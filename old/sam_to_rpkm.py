#!/usr/bin/env python3
"""
Recover RPKM from SAM
Parses a SAM file (potentially partial) and calculates RPKM based on TX:Z tags.
Requires transcript_roads.pkl for transcript lengths.
"""

import sys
import os
import argparse
import logging
import csv
import pickle
from collections import defaultdict

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_lengths(roads_file):
    logger.info(f"Loading transcript lengths from {roads_file}...")
    lengths = {}
    try:
        with open(roads_file, 'rb') as f:
            data = pickle.load(f)
            paths = data.get('transcript_paths', {})
            for tx_id, path in paths.items():
                # Stride was 5, K=31
                lengths[tx_id] = len(path) * 5 + 31
        logger.info(f"Loaded lengths for {len(lengths)} transcripts.")
        return lengths
    except Exception as e:
        logger.error(f"Failed to load roads: {e}")
        return {}

def process_sam(sam_file, lengths, output_file):
    logger.info(f"Processing SAM file {sam_file}...")
    
    transcript_counts = defaultdict(float)
    total_reads = 0
    
    try:
        with open(sam_file, 'r') as f:
            for line in f:
                if line.startswith('@'): continue
                
                parts = line.strip().split('\t')
                if len(parts) < 12: continue
                
                # Check for TX:Z tag
                tx_tag = None
                for field in parts[11:]:
                    if field.startswith('TX:Z:'):
                        tx_tag = field[5:]
                        break
                
                if tx_tag:
                    transcripts = tx_tag.split(',')
                    if transcripts:
                        weight = 1.0 / len(transcripts)
                        for tx in transcripts:
                            transcript_counts[tx] += weight
                        total_reads += 1
                        
        logger.info(f"Counted {total_reads} mapped reads.")
        
        logger.info(f"Writing RPKM to {output_file}...")
        with open(output_file, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['TranscriptID', 'Length', 'ReadCount', 'RPKM'])
            
            # Use 1M as dummy total if 0 to avoid div/0, or use actual total
            norm_total = total_reads if total_reads > 0 else 1000000
            
            # Iterate all transcripts found in SAM + known lengths
            all_ids = set(transcript_counts.keys())
            all_ids.update(lengths.keys())
            
            for tx_id in all_ids:
                count = transcript_counts[tx_id]
                length = lengths.get(tx_id, 0)
                
                rpkm = 0.0
                if length > 0 and norm_total > 0:
                    rpkm = (count * 1e9) / (length * norm_total)
                
                writer.writerow([tx_id, length, f"{count:.2f}", f"{rpkm:.4f}"])
                
        logger.info("Done.")
        
    except Exception as e:
        logger.error(f"Error processing SAM: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Recover RPKM from SAM")
    parser.add_argument("--sam", required=True, help="Input SAM file")
    parser.add_argument("--roads", required=True, help="Transcript Roads Pickle")
    parser.add_argument("--out", required=True, help="Output RPKM TSV")
    
    args = parser.parse_args()
    
    lengths = load_lengths(args.roads)
    process_sam(args.sam, lengths, args.out)
