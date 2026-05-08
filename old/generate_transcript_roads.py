#!/usr/bin/env python3
"""
Generate Transcript Roads
Maps every k-mer of every transcript to its Hilbert Index.
Creates a "Road Map" of the transcriptome in Hilbert Space.
"""

import sys
import os
import gzip
import logging
import argparse
import csv
import pickle
from collections import defaultdict
from Bio import SeqIO

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

class TranscriptRoadGenerator:
    def __init__(self, k=31):
        self.k = k
        self.point_to_transcripts = defaultdict(set) # h_idx -> {tx_ids}
        self.transcript_paths = {} # tx_id -> [h_indices]

    def process(self, gff_path, fasta_path, output_path):
        logger.info("Loading FASTA...")
        genome_seqs = {}
        try:
            if fasta_path.endswith('.gz'):
                handle = gzip.open(fasta_path, 'rt')
            else:
                handle = open(fasta_path, 'r')
            for record in SeqIO.parse(handle, "fasta"):
                genome_seqs[record.id] = str(record.seq)
            handle.close()
        except Exception as e:
            logger.error(f"Failed to load FASTA: {e}")
            return

        logger.info("Parsing GFF3 and generating roads...")
        count = 0
        
        try:
            if gff_path.endswith('.gz'):
                handle = gzip.open(gff_path, 'rt')
            else:
                handle = open(gff_path, 'r')
            
            reader = csv.reader(handle, delimiter='\t')
            
            # We need to aggregate exons for transcripts
            transcripts = defaultdict(list) # tx_id -> [exons]
            
            for row in reader:
                if not row or row[0].startswith('#'): continue
                if len(row) < 9: continue
                
                feature_type = row[2]
                if feature_type == 'exon':
                    contig = row[0]
                    start = int(row[3]) - 1
                    end = int(row[4])
                    strand = row[6]
                    attributes = row[8]
                    
                    # Extract Parent (Transcript ID)
                    parent = None
                    for attr in attributes.split(';'):
                        if attr.startswith('Parent='):
                            parent = attr.split('=')[1]
                            break
                    
                    if parent:
                        transcripts[parent].append({
                            'contig': contig,
                            'start': start,
                            'end': end,
                            'strand': strand
                        })
            
            logger.info(f"Found {len(transcripts)} transcripts. Generating paths...")
            
            for tx_id, exons in transcripts.items():
                # Sort exons
                exons.sort(key=lambda x: x['start'])
                
                # Construct transcript sequence
                tx_seq = ""
                contig = exons[0]['contig'] # Assume all exons on same contig
                strand = exons[0]['strand']
                
                if contig not in genome_seqs: continue
                
                full_seq = ""
                for exon in exons:
                    if exon['contig'] != contig: continue # Should not happen
                    seq_part = genome_seqs[contig][exon['start']:exon['end']]
                    full_seq += seq_part
                    
                if not full_seq: continue
                
                # Handle Strand
                if strand == '-':
                    # Reverse complement
                    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
                    full_seq = "".join(complement.get(base, base) for base in reversed(full_seq))
                
                if len(full_seq) < self.k: continue
                
                # Generate Path
                path = []
                # Optimization: Stride?
                # For now, let's do every k-mer to be precise, or stride=5 to save space
                stride = 5 
                
                for i in range(0, len(full_seq) - self.k + 1, stride):
                    kmer = full_seq[i:i+self.k]
                    h_idx = hilbert_index_from_seq(kmer)
                    path.append(h_idx)
                    self.point_to_transcripts[h_idx].add(tx_id)
                
                self.transcript_paths[tx_id] = path
                
                count += 1
                if count % 1000 == 0:
                    logger.info(f"Processed {count} transcripts...")
            
            logger.info(f"Saving roads to {output_path}...")
            with open(output_path, 'wb') as f:
                pickle.dump({
                    'point_to_transcripts': self.point_to_transcripts,
                    'transcript_paths': self.transcript_paths
                }, f)
            
            logger.info("Done.")
            
        except Exception as e:
            logger.error(f"Error processing GFF3: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Transcript Roads")
    parser.add_argument("--gff", required=True, help="GFF3 file")
    parser.add_argument("--fasta", required=True, help="FASTA file")
    parser.add_argument("--out", required=True, help="Output pickle file")
    
    args = parser.parse_args()
    
    generator = TranscriptRoadGenerator()
    generator.process(args.gff, args.fasta, args.out)
