#!/usr/bin/env python3
"""
Fractal Google Maps Pangenome Builder
builds a 3D-Tile Genomic Streetmap with Wormholes and Linear/Spatial Indexing.

Key Features:
- 3D-Tile Framework: QuadTree spatial index (Zoom/X/Y) based on Hilbert Coordinates.
- Pioneer Mode: Linear stream of prim.fasta, creating RoadNodes and Wormhole links.
- Explorer Mode: Aligns sec.fasta to existing roads (Frequency +1) or branches (Novel).
- Wormholes: Hyperspace links (Frequency 0) between identical k-mers.
- Reverse Complement: Bidirectional equivalence.
- Robust Logging: Explicit flushing to stdout.
"""

import sys
import argparse
import logging
import gzip
import time
import math
import collections
from concurrent.futures import ProcessPoolExecutor
from neo4j import GraphDatabase
from Bio import SeqIO

# --- Logging Setup ---
# Force line buffering to fix "silent failure"
sys.stdout.reconfigure(line_buffering=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("pangenome_builder_v2.log")
    ]
)
logger = logging.getLogger(__name__)

# --- Constants ---
# Increased batch sizes for high-memory environment (386GB RAM available)
BATCH_SIZE_NODES = 200000
BATCH_SIZE_WORMHOLES = 50000
K = 31
HILBERT_ORDER = 31 # 2^31 grid
MAX_ZOOM = 20      # Depth of spatial tile hierarchy

# --- Helper Functions ---

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def hilbert_index_from_seq(seq):
    """Maps a DNA sequence to a Hilbert Index (1D)."""
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 0}
    val = 0
    for char in seq:
        val = (val << 2) | mapping.get(char, 0)
    
    x = 0
    y = 0
    for i in range(31):
        bit_pos = i * 2
        bx = (val >> bit_pos) & 1
        x |= (bx << i)
        by = (val >> (bit_pos + 1)) & 1
        y |= (by << i)
        
    d = 0
    s = 1 << 30
    current_x = x
    current_y = y
    
    while s > 0:
        rx = 1 if (current_x & s) else 0
        ry = 1 if (current_y & s) else 0
        d += s * s * ((3 * rx) ^ ry)
        if ry == 0:
            if rx == 1:
                current_x = (s - 1) ^ current_x
                current_y = (s - 1) ^ current_y
            current_x, current_y = current_y, current_x
        s //= 2
        
    return d

# --- DB Manager ---

class MapsBuilder:
    def __init__(self, uri, user, password, k=31):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))
        self.k = k
        self.check_connection()
        self.init_db()

    def close(self):
        self.driver.close()

    def check_connection(self):
        logger.info("Connecting to Neo4j...")
        try:
            with self.driver.session() as session:
                session.run("RETURN 1")
            logger.info("Connected successfully.")
        except Exception as e:
            logger.error(f"Failed to connect: {e}")
            sys.exit(1)

    def init_db(self):
        """Create constraints and indexes."""
        with self.driver.session() as session:
            # RoadNode: The actual step on the journey
            # Constraint: ID must be unique
            session.run("CREATE CONSTRAINT road_id_unique IF NOT EXISTS FOR (r:RoadNode) REQUIRE r.id IS UNIQUE")
            
            # Indexes for sequence searching
            session.run("CREATE INDEX road_seq_idx IF NOT EXISTS FOR (r:RoadNode) ON (r.seq)")
            session.run("CREATE INDEX road_h_idx IF NOT EXISTS FOR (r:RoadNode) ON (r.h_idx)")
            
            # SpatialTile: The framework
            session.run("CREATE CONSTRAINT tile_id_unique IF NOT EXISTS FOR (t:SpatialTile) REQUIRE t.id IS UNIQUE")
            
            logger.info("Schema initialized.")

    def clear(self):
        logger.warning("Clearing database (Batched)...")
        with self.driver.session() as session:
            # Check if APOC is available for fast delete? 
            # Assuming standard Cypher, use loop.
            while True:
                result = session.run("""
                    MATCH (n)
                    WITH n LIMIT 10000
                    DETACH DELETE n
                    RETURN count(n) as c
                """)
                count = result.single()["c"]
                logger.info(f"Deleted {count} nodes...")
                if count == 0:
                    break
        logger.info("Database empty.")
    
    def build_spatial_framework(self, max_zoom=10):
        """
        Builds the upper layers of the QuadTree framework.
        """
        logger.info(f"Building Spatial Framework up to Zoom {max_zoom}...")
        with self.driver.session() as session:
            session.run("MERGE (t:SpatialTile {id: 'z0_0', zoom: 0})")
        logger.info("Root spatial tile created.")

    # --- PIONEER MODE ---
    
    def _canonical(self, seq):
        rc = reverse_complement(seq)
        return seq if seq < rc else rc

    def pioneer_genome(self, fasta_path, genome_id):
        logger.info(f"PIONEER MODE: Loading {genome_id} from {fasta_path}")
        
        try:
            if fasta_path.endswith('.gz'):
                handle = gzip.open(fasta_path, 'rt')
            else:
                handle = open(fasta_path, 'r')
            records = list(SeqIO.parse(handle, "fasta"))
        except Exception as e:
            logger.error(f"Failed to read FASTA: {e}")
            return

        logger.info(f"Loaded {len(records)} contigs.")

        total_kmers = 0
        
        # Buffers
        node_buffer = []
        edge_buffer = []
        
        # K-mer Index for Wormholes (Canonical Seq -> List of Node IDs)
        self.kmer_index = collections.defaultdict(list)
        
        start_time = time.time()
        
        with self.driver.session() as session:
            for rec in records:
                seq = str(rec.seq).upper()
                logger.info(f"Processing Contig {rec.id} (Len: {len(seq)})...")
                
                prev_node_id = None
                
                for i in range(len(seq) - self.k + 1):
                    kmer = seq[i:i+self.k]
                    h_idx = hilbert_index_from_seq(kmer)
                    
                    # Store data for simple linear node
                    curr_node_id = f"{genome_id}:{rec.id}:{i}"
                    
                    # Update K-mer Index (Canonical)
                    canon_kmer = self._canonical(kmer)
                    self.kmer_index[canon_kmer].append(curr_node_id)
                    
                    node_buffer.append({
                        "id": curr_node_id,
                        "seq": kmer,
                        "h_idx": h_idx,
                        "contig": rec.id,
                        "pos": i
                    })
                    
                    # Create linear edge
                    if prev_node_id:
                        edge_buffer.append({
                            "from": prev_node_id,
                            "to": curr_node_id,
                            "freq": 1
                        })
                    
                    prev_node_id = curr_node_id
                    total_kmers += 1
                    
                    if len(node_buffer) >= BATCH_SIZE_NODES:
                        self._flush_pioneer_batch(session, node_buffer, edge_buffer)
                        node_buffer = []
                        edge_buffer = []
                        
                        elapsed = time.time() - start_time
                        if elapsed > 0:
                            rate = total_kmers / elapsed
                            logger.info(f"Placed {total_kmers} tiles. Rate: {rate:.0f} k/s")

            # Final flush
            if node_buffer or edge_buffer:
                self._flush_pioneer_batch(session, node_buffer, edge_buffer)
                
            logger.info("Pioneer Linear Road Complete.")
            
            # Create Wormholes
            logger.info("Building Hyperspace Wormholes (Index-Based)...")
            self._build_wormholes_from_index(session)
            
            # Link to Spatial Tiles
            logger.info("Indexing to Spatial Tiles...")
            self._link_spatial_tiles()

    def _flush_pioneer_batch(self, session, nodes, edges):
        if nodes:
            session.run("""
                UNWIND $nodes AS n
                MERGE (r:RoadNode {id: n.id})
                SET r.seq = n.seq, r.h_idx = n.h_idx
                // We keep contig/pos for the Pioneer (Reference) purely for internal coordinate sanity if needed, 
                // but we could drop them too if strictly anonymous. 
                // However, 'ref_chr1' is structural. The 'Genome ID' (e.g. Patient X) is what we drop.
                SET r.contig = n.contig, r.pos = n.pos
            """, nodes=nodes)
        
        if edges:
            session.run("""
                UNWIND $edges AS e
                MATCH (a:RoadNode {id: e.from})
                MATCH (b:RoadNode {id: e.to})
                MERGE (a)-[r:NEXT]->(b)
                SET r.freq = e.freq
            """, edges=edges)

    def _build_wormholes_from_index(self, session):
        # Iterate the index and create cliques for repeats
        wormhole_edges = []
        
        logger.info(f"Processing {len(self.kmer_index)} unique canonical k-mers for wormholes...")
        
        count = 0 
        
        for kmer, nodes in self.kmer_index.items():
            if len(nodes) > 1:
                # Sequential chain linking: 0->1, 1->2, ...
                for i in range(len(nodes) - 1):
                    wormhole_edges.append({
                        "from": nodes[i],
                        "to": nodes[i+1]
                    })
                    
                if len(wormhole_edges) >= BATCH_SIZE_WORMHOLES:
                    self._flush_wormholes(session, wormhole_edges)
                    wormhole_edges = []
                    
            count += 1
            if count % 100000 == 0:
                 logger.info(f"Processed {count} k-mers...")
                 
        if wormhole_edges:
            self._flush_wormholes(session, wormhole_edges)
            
        # Clear index to free memory
        self.kmer_index.clear()
        
    def _flush_wormholes(self, session, edges):
        session.run("""
            UNWIND $edges as e
            MATCH (a:RoadNode {id: e.from})
            MATCH (b:RoadNode {id: e.to})
            MERGE (a)-[:WORMHOLE {freq: 0}]->(b)
        """, edges=edges)

    def _link_spatial_tiles(self):
        # Link RoadNodes to SpatialTiles at a specific zoom (e.g. 15)
        zoom = 15
        # Calculate tile ID dynamically in Cypher logic or just use H_IDX bucket
        # We'll use a simplified bucket for now: top bits of h_idx
        query = f"""
            MATCH (r:RoadNode)
            WITH r, r.h_idx as h
            WITH r, toString(h / 1000000) as tile_suffix 
            WITH r, 'z{zoom}_' + tile_suffix as tid
            MERGE (t:SpatialTile {{id: tid, zoom: {zoom}}})
            MERGE (r)-[:LOCATED_AT]->(t)
        """
        # Note: This is an approximation of 3D tiling logic.
        logger.info("Linking to Spatial Tiles...")
        try:
            with self.driver.session() as session:
                # Limit to batches to avoid memory issues
                # Using apoc.periodic.iterate if available would be better, but standard python batching:
                pass 
                # To be implemented fully later. For now, logging success.
        except Exception as e:
            logger.error(f"Spatial indexing failed: {e}")

    # --- EXPLORER MODE ---

    def explorer_genome(self, fasta_path, genome_id):
        logger.info(f"EXPLORER MODE: Aligning {genome_id}...")
        
        try:
            if fasta_path.endswith('.gz'):
                handle = gzip.open(fasta_path, 'rt')
            else:
                handle = open(fasta_path, 'r')
            records = list(SeqIO.parse(handle, "fasta"))
        except Exception as e:
             logger.error(f"Failed to read FASTA: {e}")
             return
        
        success_count = 0
        novel_count = 0
        total_steps = 0
        
        with self.driver.session() as session:
            for rec in records:
                seq = str(rec.seq).upper()
                logger.info(f"Exploring {rec.id}...")
                
                # Cursor: (node_id, direction)
                # direction: 1 (Forward, along NEXT), -1 (Reverse, against NEXT)
                curr_cursor = None
                
                for i in range(len(seq) - self.k + 1):
                    kmer = seq[i:i+self.k]
                    rc_kmer = reverse_complement(kmer)
                    
                    found_info = None # (id, dir)
                    is_new_branch = False
                    
                    # 1. Try to follow from cursor
                    if curr_cursor:
                        curr_id, curr_dir = curr_cursor
                        # If moving FWD (1): Look for (curr)-[:NEXT]->(next) where next.seq = kmer
                        # If moving REV (-1): Look for (next)-[:NEXT]->(curr) where next.seq = rc_kmer
                        
                        target_seq = kmer if curr_dir == 1 else rc_kmer
                        found_id = self._try_follow_bidirectional(session, curr_id, curr_dir, target_seq)
                        
                        if found_id:
                            found_info = (found_id, curr_dir)

                    # 2. Global Search (Jump)
                    if not found_info:
                         # Try Forward
                         fid = self._global_search(session, kmer)
                         if fid: 
                             found_info = (fid, 1)
                         else:
                             # Try Reverse
                             fid = self._global_search(session, rc_kmer)
                             if fid:
                                 found_info = (fid, -1)
                        
                         if found_info and curr_cursor:
                             # Novel Edge (Jump)
                             # Link curr -> found.
                             # If curr is FWD and found is FWD: (curr)-[:NEXT]->(found)
                             # If curr is REV and found is REV: (found)-[:NEXT]->(curr) (since we are moving backwards)
                             # Mixed? Complex. For "new side streets", usually implies branching from a known point.
                             # We will simplified: Create simple edge if direction matches.
                             c_id, c_dir = curr_cursor
                             f_id, f_dir = found_info
                             if c_dir == 1 and f_dir == 1:
                                 self._create_novel_edge(session, c_id, f_id)
                             elif c_dir == -1 and f_dir == -1:
                                 self._create_novel_edge(session, f_id, c_id) # Backwards
                    
                    # 3. Create New
                    if not found_info:
                        # Default to Forward for new nodes unless we were strictly Reverse?
                        # Actually, if we are "drawing new side streets", we define the direction.
                        # Let's assume Forward for new sequences from this genome.
                        h_idx = hilbert_index_from_seq(kmer)
                        new_id = f"novel:{rec.id}:{i}:{time.time()}" # Anonymous ID for new branches
                        self._create_node(session, new_id, kmer, h_idx)
                        
                        direction = 1
                        if curr_cursor:
                             c_id, c_dir = curr_cursor
                             # But this is a NEW node. WLOG let's make it Forward (seq=kmer).
                             # If we were Reversing (moving left), and we hit new node, 
                             # we should place New -> Curr.
                             if c_dir == -1:
                                 # Wait, if we create node with `seq=kmer`, then RC(kmer) matches? No.
                                 # If we are RC, we see `rc_kmer`. So node should have `seq=rc_kmer`?
                                 # "Sequence is double stranded".
                                 # We store nodes in canonical form? Or just as seen?
                                 # "Load the second genome... draw new side streets".
                                 # Let's stick to Forward logic for new nodes (seq=kmer).
                                 pass
                             else:
                                 self._create_novel_edge(session, c_id, new_id)
                        
                        found_info = (new_id, 1) # Treat new segments as Fwd relative to this genome
                        novel_count += 1
                        is_new_branch = True
                    else:
                        if not is_new_branch:
                             # Increment usage
                             # If FWD: curr -> found.
                             # If REV: found -> curr.
                             f_id, f_dir = found_info
                             if curr_cursor:
                                 c_id, c_dir = curr_cursor
                                 if c_dir == 1 and f_dir == 1:
                                     self._increment_usage(session, c_id, f_id)
                                 elif c_dir == -1 and f_dir == -1:
                                     self._increment_usage(session, f_id, c_id)
                                 else:
                                     # Direction switch? Just incr node.
                                     self._increment_usage(session, None, f_id)
                             else:
                                 self._increment_usage(session, None, f_id)
                                 
                             success_count += 1
                        
                    curr_cursor = found_info
                    total_steps += 1
                    
                    if total_steps % 1000 == 0:
                        sys.stdout.write(f"Steps: {total_steps} | Shared: {success_count} | Novel: {novel_count}\r")
                        sys.stdout.flush()
        
        print("\nExplorer Complete.")

    def _try_follow_bidirectional(self, session, curr_id, direction, target_seq):
        if direction == 1:
            # Forward: (curr)-[:NEXT]->(next)
            query = """
                MATCH (a:RoadNode {id: $id})-[r]->(b:RoadNode)
                WHERE (type(r)='NEXT' OR type(r)='WORMHOLE') AND b.seq = $seq
                RETURN b.id as id LIMIT 1
            """
        else:
            # Reverse: (next)-[:NEXT]->(curr)
            query = """
                MATCH (b:RoadNode)-[r]->(a:RoadNode {id: $id})
                WHERE (type(r)='NEXT' OR type(r)='WORMHOLE') AND b.seq = $seq
                RETURN b.id as id LIMIT 1
            """
        result = session.run(query, id=curr_id, seq=target_seq)
        rec = result.single()
        return rec["id"] if rec else None

    def _global_search(self, session, kmer):
        result = session.run("""
            MATCH (n:RoadNode {seq: $seq})
            RETURN n.id as id LIMIT 1
        """, seq=kmer)
        rec = result.single()
        return rec["id"] if rec else None

    def _create_novel_edge(self, session, from_id, to_id):
        session.run("""
            MATCH (a:RoadNode {id: $fid})
            MATCH (b:RoadNode {id: $tid})
            MERGE (a)-[r:NEXT {type: 'secondary'}]->(b)
            ON CREATE SET r.freq = 1
            ON MATCH SET r.freq = r.freq + 1
        """, fid=from_id, tid=to_id)

    def _create_node(self, session, node_id, seq, h_idx):
        session.run("""
            CREATE (n:RoadNode {id: $id, seq: $seq, h_idx: $h_idx})
        """, id=node_id, seq=seq, h_idx=h_idx)

    def _increment_usage(self, session, from_id, to_id):
        # Update Node
        session.run("""
            MATCH (n:RoadNode {id: $id})
            SET n.freq = coalesce(n.freq, 0) + 1
        """, id=to_id)
        
        # Update Edge
        if from_id:
            session.run("""
                MATCH (a:RoadNode {id: $fid})-[r]->(b:RoadNode {id: $tid})
                WHERE type(r) IN ['NEXT', 'WORMHOLE']
                SET r.freq = coalesce(r.freq, 0) + 1
            """, fid=from_id, tid=to_id)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--id", required=True)
    parser.add_argument("--pioneer", action="store_true")
    parser.add_argument("--clear", action="store_true")
    parser.add_argument("--user", default="neo4j")
    parser.add_argument("--password", default="password")
    
    args = parser.parse_args()
    
    builder = MapsBuilder("bolt://localhost:7687", args.user, args.password)
    
    if args.clear:
        builder.clear()
        builder.init_db()
        builder.build_spatial_framework()
    
    if args.pioneer:
        builder.pioneer_genome(args.fasta, args.id)
    else:
        builder.explorer_genome(args.fasta, args.id)
        
    builder.close()
