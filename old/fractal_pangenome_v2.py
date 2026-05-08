#!/usr/bin/env python3
"""
Fractal Pangenome Builder v2
- Pioneer Mode: Linear stream to DB.
- Explorer Mode: Alignment & Branching.
"""
import sys
import time

# IMMEDIATE FEEDBACK
print(f"[{time.strftime('%H:%M:%S')}] V2 SCRIPT LAUNCHING...")

import os
import gzip
import logging
import argparse
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from neo4j import GraphDatabase
from Bio import SeqIO

# --- Logging Setup ---
LOG_FILE = "pangenome_v2.log"
# Clear old log
if os.path.exists(LOG_FILE):
    os.remove(LOG_FILE)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Force unbuffered output
#sys.stdout.reconfigure(line_buffering=True)

print(f"[{time.strftime('%H:%M:%S')}] Modules imported. Logging configured to {LOG_FILE}", flush=True)

# --- Worker for Pioneer Mode ---
def worker_pioneer_chunk(sequence, k, start_id, genome_id):
    """
    Generates Nodes and Edges for a linear sequence chunk.
    """
    nodes = []
    edges = []
    seq_len = len(sequence)
    
    if seq_len < k:
        return [], []

    for i in range(seq_len - k + 1):
        kmer_seq = sequence[i:i+k]
        curr_id = start_id + i
        
        # Node
        # We use a compact JSON for kmers to save space/time, 
        # though for Pioneer we could just set properties directly if we wanted.
        # But sticking to the schema:
        kmer_data = [{"seq": kmer_seq, "count": 1, "genomes": [genome_id]}]
        
        nodes.append({
            "id": curr_id,
            "primary_seq": kmer_seq,
            "kmers": json.dumps(kmer_data)
        })
        
        # Edge to next
        if i < (seq_len - k):
            edges.append({
                "from": curr_id,
                "to": curr_id + 1,
                "genome": genome_id
            })
            
    return nodes, edges

# --- Main Builder Class ---
class PangenomeBuilderV2:
    def __init__(self, uri, user, password, k=31, workers=8):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))
        self.k = k
        self.workers = workers
        self.batch_size = 50000

    def close(self):
        self.driver.close()

    def check_connection(self):
        logger.info("Checking Neo4j connection...")
        try:
            with self.driver.session() as session:
                session.run("RETURN 1")
            logger.info("Neo4j is UP and responding.")
        except Exception as e:
            logger.error(f"Neo4j Connection Failed: {e}")
            sys.exit(1)

    def init_db(self):
        logger.info("Initializing Database Constraints...")
        with self.driver.session() as session:
            session.run("CREATE CONSTRAINT tile_id_unique IF NOT EXISTS FOR (t:Tile) REQUIRE t.id IS UNIQUE")
            session.run("CREATE INDEX tile_seq_index IF NOT EXISTS FOR (t:Tile) ON (t.primary_seq)")
        logger.info("Constraints verified.")

    def clear_database(self):
        logger.warning("!!! CLEARING DATABASE !!!")
        with self.driver.session() as session:
            # Check size first
            res = session.run("MATCH (n) RETURN count(n) as c")
            count = res.single()["c"]
            logger.info(f"Deleting {count} nodes...")
            
            deleted = 0
            while True:
                r = session.run("MATCH (n) WITH n LIMIT 50000 DETACH DELETE n RETURN count(n) as c")
                c = r.single()["c"]
                deleted += c
                logger.info(f"  Deleted {c} (Total: {deleted})")
                if c == 0: break
        logger.info("Database is empty.")

    def get_universe_stats(self):
        with self.driver.session() as session:
            res = session.run("MATCH (t:Tile) RETURN count(t) as c, coalesce(max(t.id), -1) as m")
            rec = res.single()
            return rec["c"], rec["m"]

    def run_pioneer(self, fasta_path, genome_id):
        logger.info(f"--- STARTING PIONEER MODE for {genome_id} ---")
        
        # 1. Stats
        u_size, max_id = self.get_universe_stats()
        logger.info(f"Current Universe Size: {u_size} Tiles")
        start_tile_id = max_id + 1
        logger.info(f"Will start creating new tiles at ID: {start_tile_id}")
        
        # 2. Read Fasta
        logger.info(f"Reading FASTA: {fasta_path}")
        if fasta_path.endswith('.gz'):
            handle = gzip.open(fasta_path, 'rt')
        else:
            handle = open(fasta_path, 'r')
            
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        logger.info(f"Loaded {len(records)} contigs.")
        
        # 3. Process Contigs
        current_id_counter = start_tile_id
        
        for rec in records:
            seq = str(rec.seq).upper()
            logger.info(f"Processing Contig: {rec.id} (Length: {len(seq)})")
            
            # Prepare chunks
            chunk_size = 50000
            chunks = []
            for i in range(0, len(seq), chunk_size):
                # We need overlap for k-mers? 
                # Actually, if we split the sequence, the boundary k-mers need care.
                # Simple approach: Overlap chunks by K-1.
                start = i
                end = min(i + chunk_size + self.k - 1, len(seq))
                chunk_seq = seq[start:end]
                
                # The ID for the first k-mer of this chunk
                chunk_start_id = current_id_counter + i
                
                chunks.append((chunk_seq, self.k, chunk_start_id, genome_id))
            
            logger.info(f"  Split into {len(chunks)} tasks. Launching workers...")
            
            # Execute
            total_nodes = 0
            with ProcessPoolExecutor(max_workers=self.workers) as executor:
                futures = [executor.submit(worker_pioneer_chunk, *c) for c in chunks]
                
                # We write as results come in
                with self.driver.session() as session:
                    for future in as_completed(futures):
                        nodes, edges = future.result()
                        if nodes:
                            self._write_batch(session, nodes, edges)
                            total_nodes += len(nodes)
                            print(f"  Written {total_nodes} tiles...", end="\r", flush=True)
            
            print("", flush=True) # Newline
            logger.info(f"  Contig {rec.id} Done. Created {total_nodes} tiles.")
            
            # Update counter for next contig (linear placement, no gap? or gap?)
            # User said "unconnected linear paths". So we just increment.
            # But we need to account for the fact that we consumed IDs.
            # The IDs used were start_id to start_id + len(seq) - k.
            current_id_counter += (len(seq) - self.k + 1)
            
            # Add a gap in IDs to ensure no accidental edge? 
            # IDs are just spatial. Edges define connectivity. 
            # So just incrementing is fine.
            current_id_counter += 1000 # Safety gap in "GPS space"
            
        logger.info("Pioneer Run Complete.")

    def _write_batch(self, session, nodes, edges):
        # Write Nodes
        session.run("""
            UNWIND $batch as row
            CREATE (t:Tile {id: row.id, primary_seq: row.primary_seq, kmers: row.kmers})
        """, batch=nodes)
        
        # Write Edges
        if edges:
            session.run("""
                UNWIND $batch as row
                MATCH (a:Tile {id: row.from})
                MATCH (b:Tile {id: row.to})
                CREATE (a)-[:NEXT {count: 1, genomes: [row.genome]}]->(b)
            """, batch=edges)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--id", required=True)
    parser.add_argument("--clear", action="store_true")
    parser.add_argument("--pioneer", action="store_true")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    
    print(f"[{time.strftime('%H:%M:%S')}] Arguments parsed. Initializing Builder...", flush=True)
    
    builder = PangenomeBuilderV2("bolt://localhost:7687", "neo4j", "password", workers=args.workers)
    builder.check_connection()
    
    if args.clear:
        builder.clear_database()
        builder.init_db()
        
    if args.pioneer:
        builder.run_pioneer(args.fasta, args.id)
    else:
        logger.info("Explorer mode not yet implemented in V2.")
        
    builder.close()
    print(f"[{time.strftime('%H:%M:%S')}] Script Finished.", flush=True)
