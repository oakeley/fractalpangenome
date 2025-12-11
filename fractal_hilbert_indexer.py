#!/usr/bin/env python3
"""
Fractal Hilbert Indexer
Maps 31-mer sequences to a 1D Hilbert Curve Index.
Enables spatial querying and "Digital Karyotype" visualization.
"""

import sys
import logging
import argparse
import time
from neo4j import GraphDatabase

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# --- Hilbert Curve Implementation ---
def xy2d(n, x, y):
    """
    Convert (x,y) to d (Hilbert distance).
    n: order of curve (grid size is 2^n x 2^n)
    """
    d = 0
    s = n // 2
    while s > 0:
        rx = (x & s) > 0
        ry = (y & s) > 0
        d += s * s * ((3 * rx) ^ ry)
        if ry == 0:
            if rx == 1:
                x = n - 1 - x
                y = n - 1 - y
            # Swap x and y
            x, y = y, x
        s //= 2
    return d

# Optimized Hilbert calculation for 64-bit integers
# Based on standard bit-twiddling algorithms
def hilbert_index_from_seq(seq):
    """
    Maps a DNA sequence to a Hilbert Index.
    1. Convert DNA to 2-bit integer (A=00, C=01, G=10, T=11).
    2. Split into X (even bits) and Y (odd bits).
    3. Compute Hilbert distance.
    """
    # 1. DNA to Int
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 0} # Treat N as A for indexing
    val = 0
    for char in seq:
        val = (val << 2) | mapping.get(char, 0)
    
    # 2. Split into X and Y
    # 31-mer = 62 bits.
    # X = bits at 0, 2, 4...
    # Y = bits at 1, 3, 5...
    
    x = 0
    y = 0
    
    # Simple loop for splitting (can be optimized with bit masks but loop is fine for 31 iters)
    for i in range(31):
        bit_pos = i * 2
        
        # Bit for X (even pos in val)
        bx = (val >> bit_pos) & 1
        x |= (bx << i)
        
        # Bit for Y (odd pos in val)
        by = (val >> (bit_pos + 1)) & 1
        y |= (by << i)
        
    # 3. Compute Hilbert Distance
    # Order is 31 (2^31 grid)
    # Note: Python handles large integers automatically, so 64-bit overflow isn't an issue.
    # But we need a proper Hilbert function. The recursive one is slow.
    # Let's use a standard iterative approach.
    
    d = 0
    s = 1 << 30 # Start from highest bit (order-1)
    
    current_x = x
    current_y = y
    
    # We iterate from high order to low order
    # Actually, the standard algorithm iterates s = 2^(n-1) down to 1
    
    s = 1 << 30
    while s > 0:
        rx = 1 if (current_x & s) else 0
        ry = 1 if (current_y & s) else 0
        
        d += s * s * ((3 * rx) ^ ry)
        
        if ry == 0:
            if rx == 1:
                current_x = (s - 1) ^ current_x # Invert
                current_y = (s - 1) ^ current_y # Invert
            
            # Swap
            current_x, current_y = current_y, current_x
            
        s //= 2
        
    return d

class HilbertIndexer:
    def __init__(self, uri, user, password):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))
        
    def close(self):
        self.driver.close()
        
    def create_index(self):
        with self.driver.session() as session:
            session.run("CREATE INDEX kmer_hilbert_idx IF NOT EXISTS FOR (k:Kmer) ON (k.hilbert_index)")
            logger.info("Created Neo4j index on k.hilbert_index")

    def index_all_kmers(self, batch_size=10000):
        logger.info("Starting Hilbert Indexing...")
        
        # Use a dedicated session for reading to keep the stream open
        with self.driver.session() as read_session:
            # Fetch all Kmers that don't have an index yet
            query = "MATCH (k:Kmer) WHERE k.hilbert_index IS NULL RETURN k.seq as seq"
            
            result = read_session.run(query)
            
            batch = []
            count = 0
            total_processed = 0
            start_time = time.time()
            
            for record in result:
                seq = record['seq']
                h_idx = hilbert_index_from_seq(seq)
                
                batch.append({'seq': seq, 'h_idx': h_idx})
                
                if len(batch) >= batch_size:
                    # Use a SEPARATE session for writing to avoid conflict with the read stream
                    self._update_batch_safe(batch)
                    total_processed += len(batch)
                    batch = []
                    
                    if total_processed % 100000 == 0:
                        elapsed = time.time() - start_time
                        rate = total_processed / elapsed
                        logger.info(f"Indexed {total_processed} kmers. Rate: {rate:.0f} kmers/s")
            
            if batch:
                self._update_batch_safe(batch)
                total_processed += len(batch)
                
            logger.info(f"Finished. Total indexed: {total_processed}")

    def _update_batch_safe(self, batch):
        """Opens a new session to perform the write transaction safely."""
        with self.driver.session() as write_session:
            write_session.execute_write(self._update_batch_tx, batch)

    @staticmethod
    def _update_batch_tx(tx, batch):
        tx.run("""
            UNWIND $batch AS item
            MATCH (k:Kmer {seq: item.seq})
            SET k.hilbert_index = item.h_idx
        """, batch=batch)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fractal Hilbert Indexer")
    parser.add_argument("--user", default="neo4j", help="Neo4j user")
    parser.add_argument("--password", default="password", help="Neo4j password")
    
    args = parser.parse_args()
    
    indexer = HilbertIndexer("bolt://localhost:7687", args.user, args.password)
    indexer.create_index()
    indexer.index_all_kmers()
    indexer.close()
