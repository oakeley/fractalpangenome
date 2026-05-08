#!/usr/bin/env python3
"""
Map Genes to Hilbert Intervals
Generates a mapping file (JSON/Pickle) linking Gene IDs to their Hilbert Coordinate Ranges.
Used for RPKM quantification in the Graph Aligner.
"""

import sys
import gzip
import argparse
import logging
import json
import pickle
from Bio import SeqIO
import csv

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# --- Hilbert Curve Logic (Duplicated for standalone usage) ---
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

def generate_gene_map(gff_path, fasta_path, output_path, k=31):
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

    logger.info(f"Parsing GFF3 from {gff_path}...")
    gene_intervals = {} # gene_id -> {'min': inf, 'max': -inf, 'length': 0}
    
    try:
        handle = gzip.open(gff_path, 'rt') if gff_path.endswith('.gz') else open(gff_path, 'r')
        reader = csv.reader(handle, delimiter='\t')
        
        count = 0
        for row in reader:
            if not row or row[0].startswith('#'): continue
            if len(row) < 9: continue
            
            feature_type = row[2]
            if feature_type not in ['gene', 'mRNA', 'transcript']: continue
            
            contig = row[0]
            start = int(row[3]) - 1
            end = int(row[4])
            attributes_str = row[8]
            
            # Parse ID
            attributes = {}
            for attr in attributes_str.split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attributes[key] = val
            
            feat_id = attributes.get('ID', '')
            if not feat_id: continue
            
            # Get Sequence
            if contig in genome_seqs:
                seq = genome_seqs[contig][start:end]
                if len(seq) < k: continue
                
                # Compute Hilbert Range
                # Optimization: We don't need EVERY k-mer. 
                # Just sampling start, middle, end, plus some steps might be enough for an approximation,
                # BUT for accurate RPKM we want the true bounding box.
                # Hilbert curve preserves locality but is fractal. 
                # A linear segment might jump around in Hilbert space.
                # We should sample reasonably densely.
                
                min_h = float('inf')
                max_h = float('-inf')
                
                # Sample every K bases? Or every 10 bases?
                step = max(1, len(seq) // 100) # Sample 100 points per gene max
                
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
                    
        handle.close()
        
        logger.info(f"Saving map for {len(gene_intervals)} genes to {output_path}...")
        with open(output_path, 'wb') as f:
            pickle.dump(gene_intervals, f)
            
    except Exception as e:
        logger.error(f"Error processing GFF3: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map Genes to Hilbert Intervals")
    parser.add_argument("--gff", required=True, help="Path to GFF3 file")
    parser.add_argument("--fasta", required=True, help="Path to FASTA file")
    parser.add_argument("--out", required=True, help="Output Pickle file")
    
    args = parser.parse_args()
    
    generate_gene_map(args.gff, args.fasta, args.out)
