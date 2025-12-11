#!/usr/bin/env python3
"""
Map Genes to Hilbert Intervals (Graph Source)
Generates a mapping file (Pickle) linking Gene IDs to their Hilbert Coordinate Ranges.
Queries Neo4j for Gene coordinates and uses FASTA for sequence extraction.
Used for RPKM quantification in the Graph Aligner.
"""

import sys
import gzip
import argparse
import logging
import pickle
from Bio import SeqIO
from neo4j import GraphDatabase

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# --- Hilbert Curve Logic ---
def hilbert_index_from_seq(seq):
    """Maps a DNA sequence to a Hilbert Index."""
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

class GeneMapper:
    def __init__(self, uri, user, password):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        self.driver.close()

    def generate_map(self, fasta_path, output_path, genome_id, k=31):
        logger.info(f"Loading FASTA from {fasta_path}...")
        genome_seqs = {}
        try:
            handle = gzip.open(fasta_path, 'rt') if fasta_path.endswith('.gz') else open(fasta_path, 'r')
            for record in SeqIO.parse(handle, "fasta"):
                genome_seqs[record.id] = str(record.seq)
            handle.close()
        except Exception as e:
            logger.error(f"Failed to load FASTA: {e}")
            return

        logger.info(f"Querying Genes for {genome_id} from Neo4j...")
        gene_intervals = {} # gene_id -> {'min': inf, 'max': -inf, 'length': 0}
        
        with self.driver.session() as session:
            # Fetch all genes for the given genome
            result = session.run("""
                MATCH (f:Feature {type: 'gene', genome: $genome})
                RETURN f.id AS id, f.contig AS contig, f.start AS start, f.end AS end
            """, genome=genome_id)
            
            count = 0
            for record in result:
                feat_id = record['id']
                contig = record['contig']
                start = record['start']
                end = record['end']
                
                if contig in genome_seqs:
                    # Extract sequence
                    # Ensure coordinates are within bounds
                    seq_len = len(genome_seqs[contig])
                    if start < 0 or end > seq_len:
                        logger.warning(f"Gene {feat_id} out of bounds: {start}-{end} vs {seq_len}")
                        continue
                        
                    seq = genome_seqs[contig][start:end]
                    if len(seq) < k: continue
                    
                    # Compute Hilbert Range
                    min_h = float('inf')
                    max_h = float('-inf')
                    
                    # Sample reasonably densely
                    step = max(1, len(seq) // 100) 
                    
                    for i in range(0, len(seq) - k + 1, step):
                        kmer = seq[i:i+k]
                        h_idx = hilbert_index_from_seq(kmer)
                        if h_idx < min_h: min_h = h_idx
                        if h_idx > max_h: max_h = h_idx
                    
                    # Also check exact start and end k-mers
                    start_kmer = seq[:k]
                    end_kmer = seq[-k:]
                    h_start = hilbert_index_from_seq(start_kmer)
                    h_end = hilbert_index_from_seq(end_kmer)
                    
                    min_h = min(min_h, h_start, h_end)
                    max_h = max(max_h, h_start, h_end)
                    
                    gene_intervals[feat_id] = {
                        'min': min_h,
                        'max': max_h,
                        'length': len(seq)
                    }
                    
                    count += 1
                    if count % 1000 == 0:
                        logger.info(f"Mapped {count} genes...")
        
        logger.info(f"Saving map for {len(gene_intervals)} genes to {output_path}...")
        with open(output_path, 'wb') as f:
            pickle.dump(gene_intervals, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map Genes to Hilbert Intervals (Graph Source)")
    parser.add_argument("--fasta", required=True, help="Path to FASTA file")
    parser.add_argument("--out", required=True, help="Output Pickle file")
    parser.add_argument("--genome", required=True, help="Genome ID (e.g., T2T-CHM13)")
    parser.add_argument("--user", default="neo4j", help="Neo4j user")
    parser.add_argument("--password", default="password", help="Neo4j password")
    
    args = parser.parse_args()
    
    mapper = GeneMapper("bolt://localhost:7687", args.user, args.password)
    mapper.generate_map(args.fasta, args.out, args.genome)
    mapper.close()
