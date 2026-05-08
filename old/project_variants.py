#!/usr/bin/env python3
"""
Project Pangenome Variants
Maps the "World View" of variants from the Fractal Graph onto a specific
Linear Assembly (FASTA) to generate a VCF.

Usage:
    python3 project_variants.py --fasta my_assembly.fasta --out my_variants.vcf
"""

import sys
import argparse
import logging
import gzip
from collections import defaultdict
from neo4j import GraphDatabase
from Bio import SeqIO

# Logging
sys.stdout.reconfigure(line_buffering=True)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Reusing the Universal Hilbert Logic
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

class VariantProjector:
    def __init__(self, uri, user, password, k=31):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))
        self.k = k

    def close(self):
        self.driver.close()

    def process_assembly(self, fasta_path, vcf_path):
        logger.info(f"Projecting variants onto {fasta_path}...")
        
        try:
            if fasta_path.endswith('.gz'):
                handle = gzip.open(fasta_path, 'rt')
            else:
                handle = open(fasta_path, 'r')
                
            with open(vcf_path, 'w') as vcf:
                # Header
                vcf.write("##fileformat=VCFv4.2\n")
                vcf.write(f"##source=FractalPangenome\n")
                vcf.write(f"##reference={fasta_path}\n")
                vcf.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency in Pangenome Graph\">\n")
                vcf.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth (Frequency) in Graph\">\n")
                vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                
                for record in SeqIO.parse(handle, "fasta"):
                    self._process_contig(record.id, str(record.seq), vcf)
                    
            handle.close()
            logger.info("Done.")
            
        except Exception as e:
            logger.error(f"Failed: {e}")

    def _process_contig(self, chrom, seq, vcf):
        seq_len = len(seq)
        if seq_len < self.k: return
        
        logger.info(f"Processing {chrom} ({seq_len} bp)...")
        
        # We assume 1-based coordinates for VCF
        # The K-mer starts at `i` (0-indexed). The variant is roughly at `i + k` (the next base decision)
        # Actually, if we are at Node K1, the decision is "What is the next base?".
        # Node K1 corresponds to seq[i : i+k].
        # The next base is seq[i+k].
        # So at POS = i + k + 1 (1-based), the Ref allele is seq[i+k].
        # We check the graph for other outgoing edges.
        
        # Batching
        BATCH_SIZE = 1000
        
        for i in range(0, seq_len - self.k, BATCH_SIZE):
            batch_end = min(i + BATCH_SIZE, seq_len - self.k)
            
            # Prepare batch of K-mers and their expected Next Base
            batch_data = []
            for j in range(i, batch_end):
                kmer = seq[j:j+self.k]
                ref_next = seq[j+self.k]
                h_idx = hilbert_index_from_seq(kmer)
                batch_data.append({'h_idx': h_idx, 'ref': ref_next, 'pos': j+self.k+1, 'seq': kmer})
                
            # Query Graph
            variants = self._query_batch(batch_data)
            
            # Write Variants
            for item in batch_data:
                h_idx = item['h_idx']
                if h_idx in variants:
                    node_data = variants[h_idx]
                    
                    # node_data is list of (next_base, freq)
                    # Calculate Total Depth
                    total_depth = sum(x[1] for x in node_data)
                    if total_depth == 0: continue
                    
                    ref_base = item['ref']
                    
                    # Find Alt Alleles
                    alts = []
                    alt_freqs = []
                    
                    for base, freq in node_data:
                        if base != ref_base:
                            alts.append(base)
                            alt_freqs.append(freq / total_depth)
                            
                    if alts:
                        alt_str = ",".join(alts)
                        af_str = ",".join(f"{f:.4f}" for f in alt_freqs)
                        info = f"DP={total_depth};AF={af_str}"
                        
                        vcf.write(f"{chrom}\t{item['pos']}\t.\t{ref_base}\t{alt_str}\t.\tPASS\t{info}\n")
                        
            if i % 10000 == 0:
                 print(f"Scanned {i}/{seq_len}...", end='\r')

    def _query_batch(self, batch_data):
        # Extract H-indices
        h_indices = [x['h_idx'] for x in batch_data]
        
        variants = {} # h_idx -> [(base, freq)]
        
        with self.driver.session() as session:
            result = session.run("""
                UNWIND $h_indices AS h
                MATCH (n:RoadNode {h_idx: h})
                MATCH (n)-[r:NEXT]->(m:RoadNode)
                RETURN h, m.seq AS next_seq, r.freq AS freq
            """, h_indices=h_indices)
            
            for record in result:
                h = record['h']
                next_seq = record['next_seq']
                freq = record['freq']
                
                # Extract Last Base from next_seq
                # next_seq is k-mer. It overlaps n.seq by k-1. The last char is the new base.
                if not next_seq: continue
                base = next_seq[-1]
                
                if h not in variants: variants[h] = []
                variants[h].append((base, freq))
                
        return variants

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--user", default="neo4j")
    parser.add_argument("--password", default="password")
    
    args = parser.parse_args()
    
    projector = VariantProjector("bolt://localhost:7687", args.user, args.password)
    projector.process_assembly(args.fasta, args.out)
    projector.close()
