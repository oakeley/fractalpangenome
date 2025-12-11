#!/usr/bin/env python3
"""
Fractal Neo4j Pangenome Builder
Builds a k-mer pangenome graph in Neo4j from FASTA files.
Single-Writer / Multi-Worker architecture to prevent deadlocks.
"""

import os
import sys
import gzip
import logging
import argparse
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from neo4j import GraphDatabase
from Bio import SeqIO

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# --- Worker Function (Pure CPU) ---
def worker_process_chunk(kmers, genome_id, contig_id, start_global_idx, is_last_chunk):
    """
    Process a chunk of k-mers:
    1. Identify unique nodes.
    2. Aggregate transitions (edges).
    Returns: (unique_nodes, edge_list, start_kmer_seq)
    """
    # 1. Unique Nodes
    unique_kmers = list(set(kmers))
    
    # 2. Unique Transitions
    # Key: (from_seq, to_seq) -> count
    transitions = {}
    
    process_limit = len(kmers) if is_last_chunk else len(kmers) - 1
    
    for i in range(process_limit):
        k1 = kmers[i]
        if i + 1 < len(kmers):
            k2 = kmers[i+1]
            key = (k1, k2)
            if key in transitions:
                transitions[key] += 1
            else:
                transitions[key] = 1
    
    # Convert transitions to list
    edge_list = []
    for (k1, k2), count in transitions.items():
        edge_list.append({
            'from': k1,
            'to': k2,
            'count': count,
            'genome': genome_id
        })
        
    start_kmer = kmers[0] if (start_global_idx == 0 and len(kmers) > 0) else None
    
    return unique_kmers, edge_list, start_kmer

# --- Writer Class (DB IO) ---
class Neo4jWriter:
    def __init__(self, uri, auth):
        self.driver = GraphDatabase.driver(uri, auth=auth)
        
    def close(self):
        self.driver.close()
        
    def init_db(self):
        with self.driver.session() as session:
            session.run("CREATE CONSTRAINT kmer_seq_unique IF NOT EXISTS FOR (k:Kmer) REQUIRE k.seq IS UNIQUE")
            session.run("CREATE INDEX genome_id_idx IF NOT EXISTS FOR (g:Genome) ON (g.id)")
            session.run("CREATE INDEX contig_id_idx IF NOT EXISTS FOR (c:Contig) ON (c.id)")
            logger.info("Database constraints and indexes initialized.")

    def create_genome_node(self, genome_id):
        with self.driver.session() as session:
            session.run("MERGE (g:Genome {id: $id})", id=genome_id)

    def create_contig_node(self, genome_id, contig_id):
        with self.driver.session() as session:
            session.run("""
                MATCH (g:Genome {id: $gid})
                MERGE (c:Contig {id: $cid, genome_id: $gid})
                MERGE (g)-[:HAS_CONTIG]->(c)
            """, gid=genome_id, cid=contig_id)

    def check_contig_status(self, genome_id, contig_id):
        with self.driver.session() as session:
            result = session.run("""
                MATCH (c:Contig {id: $cid, genome_id: $gid})
                RETURN c.status AS status
            """, gid=genome_id, cid=contig_id)
            record = result.single()
            return record["status"] if record else None

    def update_contig_status(self, genome_id, contig_id, status):
        with self.driver.session() as session:
            session.run("""
                MATCH (c:Contig {id: $cid, genome_id: $gid})
                SET c.status = $status
            """, gid=genome_id, cid=contig_id, status=status)

    def write_batch(self, nodes, edges, start_kmer_info):
        """
        Write a batch of nodes and edges.
        start_kmer_info: (contig_id, genome_id, start_seq) or None
        """
        if not nodes and not edges and not start_kmer_info:
            return

        with self.driver.session() as session:
            # 1. Nodes
            if nodes:
                session.execute_write(self._create_nodes_tx, nodes)
            
            # 2. Edges
            if edges:
                session.execute_write(self._create_edges_tx, edges)
                
            # 3. Start Link
            if start_kmer_info:
                cid, gid, seq = start_kmer_info
                session.execute_write(self._link_contig_start_tx, seq, cid, gid)

    @staticmethod
    def _create_nodes_tx(tx, batch_nodes):
        tx.run("""
            UNWIND $batch AS seq
            MERGE (k:Kmer {seq: seq})
        """, batch=batch_nodes)

    @staticmethod
    def _create_edges_tx(tx, batch_edges):
        # Note: We assume all edges in this batch belong to the same genome if we wanted to optimize
        # But the list contains 'genome' key, so we use it.
        tx.run("""
            UNWIND $batch AS item
            MATCH (k1:Kmer {seq: item.from})
            MATCH (k2:Kmer {seq: item.to})
            MERGE (k1)-[r:NEXT]->(k2)
            ON CREATE SET 
                r.count = item.count,
                r.genomes = [item.genome]
            ON MATCH SET 
                r.count = r.count + item.count,
                r.genomes = CASE 
                    WHEN r.genomes IS NULL THEN [item.genome]
                    WHEN item.genome IN r.genomes THEN r.genomes 
                    ELSE r.genomes + item.genome 
                END
        """, batch=batch_edges)

    @staticmethod
    def _link_contig_start_tx(tx, start_seq, contig_id, genome_id):
        tx.run("""
            MATCH (k:Kmer {seq: $seq})
            MERGE (c:Contig {id: $cid, genome_id: $gid})
            MERGE (c)-[:STARTS_AT]->(k)
        """, seq=start_seq, cid=contig_id, gid=genome_id)


# --- Main Builder ---
class Neo4jPangenomeBuilder:
    def __init__(self, uri, user, password, k=31, workers=16):
        self.k = k
        self.workers = workers
        self.writer = Neo4jWriter(uri, (user, password))
        self.writer.init_db()

    def close(self):
        self.writer.close()

    def add_genome(self, fasta_path, genome_id):
        logger.info(f"Adding genome {genome_id} from {fasta_path}...")
        self.writer.create_genome_node(genome_id)
        
        try:
            if fasta_path.endswith('.gz'):
                handle = gzip.open(fasta_path, 'rt')
            else:
                handle = open(fasta_path, 'r')
            
            for record in SeqIO.parse(handle, "fasta"):
                # Check if already done
                status = self.writer.check_contig_status(genome_id, record.id)
                if status == 'completed':
                    logger.info(f"Contig {record.id} already completed. Skipping.")
                    continue
                    
                self._process_contig(genome_id, record.id, str(record.seq))
                
            handle.close()
            logger.info(f"Genome {genome_id} added successfully.")
            
        except Exception as e:
            logger.error(f"Error processing genome: {e}")

    def _process_contig(self, genome_id, contig_id, sequence):
        logger.info(f"Processing contig {contig_id} (length: {len(sequence)})...")
        self.writer.create_contig_node(genome_id, contig_id)
        
        # Generate k-mers
        logger.info("Generating k-mers...")
        kmers = []
        for i in range(len(sequence) - self.k + 1):
            kmers.append(sequence[i:i+self.k])
            
        total_kmers = len(kmers)
        logger.info(f"Generated {total_kmers} k-mers. Starting processing ({self.workers} workers)...")
        
        # Chunking
        chunk_size = 50000 
        tasks = []
        for i in range(0, total_kmers, chunk_size):
            end_idx = min(i + chunk_size + 1, total_kmers)
            chunk = kmers[i : end_idx]
            is_last_chunk = (end_idx == total_kmers)
            tasks.append((chunk, genome_id, contig_id, i, is_last_chunk))
        
        logger.info(f"Created {len(tasks)} tasks.")
        
        # Execution & Writing Loop
        start_time = time.time()
        completed_kmers = 0
        
        # Buffers for batch writing
        node_buffer = []
        edge_buffer = []
        start_link_buffer = None
        
        WRITE_BATCH_SIZE = 20000 # Items to accumulate before writing
        
        with ProcessPoolExecutor(max_workers=self.workers) as executor:
            futures = [executor.submit(worker_process_chunk, *task) for task in tasks]
            
            for future in as_completed(futures):
                try:
                    nodes, edges, start_kmer = future.result()
                    
                    # Accumulate
                    node_buffer.extend(nodes)
                    edge_buffer.extend(edges)
                    if start_kmer:
                        start_link_buffer = (contig_id, genome_id, start_kmer)
                    
                    completed_kmers += len(nodes) # Approx
                    
                    # Write if buffer full
                    if len(node_buffer) >= WRITE_BATCH_SIZE or len(edge_buffer) >= WRITE_BATCH_SIZE:
                        self.writer.write_batch(node_buffer, edge_buffer, start_link_buffer)
                        # Clear buffers
                        node_buffer = []
                        edge_buffer = []
                        start_link_buffer = None 
                        
                    # Progress log
                    if completed_kmers % 500000 < 50000: # Log roughly every 500k
                         elapsed = time.time() - start_time
                         rate = completed_kmers / elapsed if elapsed > 0 else 0
                         logger.info(f"Progress: {completed_kmers}/{total_kmers} ({completed_kmers/total_kmers*100:.1f}%) - Rate: {rate:.0f} kmers/s")
                         
                except Exception as e:
                    logger.error(f"Task failed: {e}")
            
            # Flush remaining
            if node_buffer or edge_buffer or start_link_buffer:
                self.writer.write_batch(node_buffer, edge_buffer, start_link_buffer)

        total_time = time.time() - start_time
        logger.info(f"Contig {contig_id} complete. Total time: {total_time:.1f}s.")
        
        # Mark as completed
        self.writer.update_contig_status(genome_id, contig_id, 'completed')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build Neo4j Pangenome (Single-Writer)")
    parser.add_argument("--fasta", required=True, help="Path to FASTA file")
    parser.add_argument("--genome", required=True, help="Genome ID")
    parser.add_argument("--user", default="neo4j", help="Neo4j user")
    parser.add_argument("--password", default="password", help="Neo4j password")
    parser.add_argument("--workers", type=int, default=32, help="Number of worker processes")
    
    args = parser.parse_args()
    
    try:
        import neo4j
        import Bio
    except ImportError:
        print("Missing dependencies. Run: pip install neo4j biopython")
        sys.exit(1)

    builder = Neo4jPangenomeBuilder("bolt://localhost:7687", args.user, args.password, workers=args.workers)
    builder.add_genome(args.fasta, args.genome)
    builder.close()
