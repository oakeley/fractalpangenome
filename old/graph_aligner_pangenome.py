#!/usr/bin/env python3
"""
Custom Graph Aligner for FractalPangenome (Neo4j Version)
Aligns reads to the Neo4j pangenome graph using Universal Hilbert Coordinates.
Supports:
- Single-End and Paired-End reads
- Multiprocessing (Optimized)
- Graph-based CIGAR (Variant Calling)
- RPKM Quantification (Transcript Roads)
- Reverse Complement Alignment
- Transcript Tagging (TX:Z)
"""

import sys
import os
import logging
import argparse
import gzip
import time
import statistics
import csv
import pickle
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed, FIRST_COMPLETED, wait
from neo4j import GraphDatabase

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# --- Helper Functions ---
def reverse_complement(seq):
    """Returns the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

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

# --- Gene Quantifier ---
class GeneQuantifier:
    def __init__(self, uri, user, password, roads_file=None):
        self.point_to_transcripts = {} # h_idx -> {tx_ids}
        self.transcript_counts = defaultdict(float) # Fractional counts
        self.transcript_lengths = {}
        
        if roads_file:
            self._load_roads(roads_file)
        else:
            # Fallback to Neo4j intervals (Legacy)
            self.gene_intervals = [] 
            self._load_genes_neo4j(uri, user, password)

    def _load_roads(self, roads_file):
        logger.info(f"Loading Transcript Roads from {roads_file}...")
        try:
            with open(roads_file, 'rb') as f:
                data = pickle.load(f)
                self.point_to_transcripts = data['point_to_transcripts']
                # Estimate lengths from paths
                paths = data.get('transcript_paths', {})
                for tx_id, path in paths.items():
                    # Stride was 5, K=31
                    # Length ~ len(path) * 5 + 31
                    self.transcript_lengths[tx_id] = len(path) * 5 + 31
            logger.info(f"Loaded roads for {len(self.transcript_lengths)} transcripts.")
        except Exception as e:
            logger.error(f"Failed to load roads: {e}")

    def _load_genes_neo4j(self, uri, user, password):
        logger.info("Loading gene definitions from Neo4j (Legacy Mode)...")
        driver = GraphDatabase.driver(uri, auth=(user, password))
        with driver.session() as session:
            result = session.run("""
                MATCH (f:Feature {type: 'gene'})
                WHERE f.min_hilbert IS NOT NULL AND f.max_hilbert IS NOT NULL
                RETURN f.id AS id, f.min_hilbert AS min_h, f.max_hilbert AS max_h, (f.end - f.start) AS len
            """)
            for record in result:
                self.gene_intervals.append((
                    record['min_h'], 
                    record['max_h'], 
                    record['id'],
                    record['len']
                ))
        driver.close()
        self.gene_intervals.sort(key=lambda x: x[0])

    def count_read(self, candidates):
        """
        Increment count for transcripts matching candidates.
        candidates: list of Hilbert indices found in the read.
        """
        if self.point_to_transcripts:
            # Precise Mode
            matched_txs = set()
            for h_idx in candidates:
                if h_idx in self.point_to_transcripts:
                    matched_txs.update(self.point_to_transcripts[h_idx])
            
            if matched_txs:
                weight = 1.0 / len(matched_txs)
                for tx in matched_txs:
                    self.transcript_counts[tx] += weight
        else:
            # Legacy Interval Mode (using first candidate as pos)
            if not candidates: return
            pos = candidates[0]
            for min_h, max_h, gene_id, _ in self.gene_intervals:
                if min_h <= pos <= max_h:
                    self.transcript_counts[gene_id] += 1

    def get_overlapping_transcripts(self, candidates):
        """Returns list of transcript IDs."""
        if self.point_to_transcripts:
            matched_txs = set()
            for h_idx in candidates:
                if h_idx in self.point_to_transcripts:
                    matched_txs.update(self.point_to_transcripts[h_idx])
            return list(matched_txs)
        else:
            # Legacy
            if not candidates: return []
            pos = candidates[0]
            genes = []
            for min_h, max_h, gene_id, _ in self.gene_intervals:
                if min_h <= pos <= max_h:
                    genes.append(gene_id)
            return genes

    def write_rpkm(self, output_path, total_reads):
        logger.info(f"Writing RPKM to {output_path}...")
        with open(output_path, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['TranscriptID', 'Length', 'ReadCount', 'RPKM'])
            
            # Iterate all known transcripts
            all_ids = set(self.transcript_counts.keys())
            if self.point_to_transcripts:
                all_ids.update(self.transcript_lengths.keys())
            else:
                all_ids.update(g[2] for g in self.gene_intervals)
            
            for tx_id in all_ids:
                count = self.transcript_counts[tx_id]
                length = self.transcript_lengths.get(tx_id, 0)
                
                # Fallback length for legacy
                if length == 0 and not self.point_to_transcripts:
                     for g in self.gene_intervals:
                         if g[2] == tx_id:
                             length = g[3]
                             break
                
                rpkm = 0.0
                if length > 0 and total_reads > 0:
                    rpkm = (count * 1e9) / (length * total_reads)
                
                writer.writerow([tx_id, length, f"{count:.2f}", f"{rpkm:.4f}"])

# --- Worker Function ---
def worker_align_batch(batch_reads, uri, user, password, k):
    driver = GraphDatabase.driver(uri, auth=(user, password))
    results = []
    try:
        for r1, r2 in batch_reads:
            res1 = _align_single_read_best_strand(driver, r1, k)
            res2 = _align_single_read_best_strand(driver, r2, k) if r2 else None
            results.append((res1, res2))
    except Exception as e:
        logger.error(f"Worker failed: {e}")
    finally:
        driver.close()
    return results

def _align_single_read_best_strand(driver, read, k):
    if not read: return None
    
    res_fwd = _align_sequence(driver, read['seq'], k)
    seq_rc = reverse_complement(read['seq'])
    res_rc = _align_sequence(driver, seq_rc, k)
    
    score_fwd = res_fwd['score'] if res_fwd else -1
    score_rc = res_rc['score'] if res_rc else -1
    
    best_res = None
    is_rc = False
    
    if score_fwd >= score_rc and score_fwd > 0:
        best_res = res_fwd
        is_rc = False
    elif score_rc > score_fwd and score_rc > 0:
        best_res = res_rc
        is_rc = True
    
    if best_res:
        return {
            'id': read['id'],
            'seq': best_res['seq'], 
            'qual': read['qual'], 
            'score': best_res['score'],
            'cigar': best_res['cigar'],
            'mapq': 60,
            'rname': 'Hilbert_Space',
            'pos': best_res['pos'],
            'candidates': best_res['candidates'],
            'is_mapped': True,
            'is_rc': is_rc
        }
    else:
        return {
            'id': read['id'],
            'seq': read['seq'],
            'qual': read['qual'],
            'score': 0,
            'cigar': '*',
            'mapq': 0,
            'rname': '*',
            'pos': 0,
            'candidates': [],
            'is_mapped': False,
            'is_rc': False
        }

def _align_sequence(driver, seq, k):
    if len(seq) < k: return None
    
    read_len = len(seq)
    read_kmers = []
    for i in range(read_len - k + 1):
        read_kmers.append(seq[i:i+k])
    
    unique_read_kmers = list(set(read_kmers))
    found_indices = []
    found_kmers_set = set()
    
    BATCH_SIZE = 500
    for i in range(0, len(unique_read_kmers), BATCH_SIZE):
        batch = unique_read_kmers[i : i + BATCH_SIZE]
        indices_map = _get_kmer_indices_map(driver, batch)
        for kmer_seq, h_idx in indices_map.items():
            found_kmers_set.add(kmer_seq)
            if h_idx is None:
                h_idx = hilbert_index_from_seq(kmer_seq)
            found_indices.append(h_idx)
            
    match_count = len(found_indices)
    total_kmers = len(unique_read_kmers)
    
    if total_kmers > 0 and (match_count / total_kmers) > 0.5:
        if found_indices:
            # Use actual element for median to ensure it's a real k-mer
            found_indices.sort()
            median_idx = found_indices[len(found_indices)//2]
            
            # Select candidates for transcript matching
            # Pick a few points spread across the read
            candidates = []
            if found_indices:
                candidates.append(found_indices[0])
                candidates.append(found_indices[-1])
                candidates.append(median_idx)
                # Add quartiles
                candidates.append(found_indices[len(found_indices)//4])
                candidates.append(found_indices[3*len(found_indices)//4])
                
            cigar = _generate_graph_cigar(seq, k, found_kmers_set)
            return {
                'score': match_count,
                'cigar': cigar,
                'pos': median_idx,
                'candidates': list(set(candidates)),
                'seq': seq
            }
    return None

def _get_kmer_indices_map(driver, kmer_list):
    with driver.session() as session:
        result = session.run("""
            UNWIND $kmers AS seq
            MATCH (k:RoadNode {seq: seq})
            RETURN k.seq AS seq, k.h_idx AS h_idx LIMIT 1
        """, kmers=kmer_list)
        return {record['seq']: record['h_idx'] for record in result}

def _generate_graph_cigar(seq, k, found_kmers_set):
    read_len = len(seq)
    supported = [False] * read_len
    for i in range(read_len - k + 1):
        kmer = seq[i:i+k]
        if kmer in found_kmers_set:
            for j in range(i, i+k):
                supported[j] = True
    cigar_parts = []
    if not supported: return "*"
    current_state = supported[0]
    count = 0
    for is_supp in supported:
        if is_supp == current_state:
            count += 1
        else:
            op = 'M' if current_state else 'X'
            cigar_parts.append(f"{count}{op}")
            current_state = is_supp
            count = 1
    op = 'M' if current_state else 'X'
    cigar_parts.append(f"{count}{op}")
    return "".join(cigar_parts)

# --- Main Aligner Class ---
class Neo4jGraphAligner:
    def __init__(self, uri, user, password, k=31, cores=4, rpkm_out=None, roads_file=None):
        self.uri = uri
        self.user = user
        self.password = password
        self.k = k
        self.cores = cores
        self.rpkm_out = rpkm_out
        
        self.quantifier = None
        if rpkm_out or roads_file:
            self.quantifier = GeneQuantifier(uri, user, password, roads_file)

    def align_fastq(self, fastq1, fastq2, output_sam):
        logger.info(f"Aligning reads to {output_sam} using {self.cores} cores...")
        try:
            f1 = gzip.open(fastq1, 'rt') if fastq1.endswith('.gz') else open(fastq1, 'r')
            f2 = (gzip.open(fastq2, 'rt') if fastq2.endswith('.gz') else open(fastq2, 'r')) if fastq2 else None
            
            with open(output_sam, 'w') as f_out:
                f_out.write("@HD\tVN:1.6\tSO:unknown\n")
                f_out.write(f"@SQ\tSN:Hilbert_Space\tLN:2147483647\n") 
                
                BATCH_SIZE = 2000
                running_futures = set()
                
                with ProcessPoolExecutor(max_workers=self.cores) as executor:
                    while True:
                        while len(running_futures) < self.cores * 2:
                            chunk = self._read_chunk(f1, f2, BATCH_SIZE)
                            if not chunk: break
                            future = executor.submit(worker_align_batch, chunk, self.uri, self.user, self.password, self.k)
                            running_futures.add(future)
                        
                        if not running_futures: break
                        
                        done, _ = wait(running_futures, return_when=FIRST_COMPLETED)
                        for future in done:
                            running_futures.remove(future)
                            results = future.result()
                            for res1, res2 in results:
                                self._write_sam_record(f_out, res1, res2)
                                if self.quantifier:
                                    if res1 and res1['is_mapped']:
                                        self.quantifier.count_read(res1['candidates'])
                                    if res2 and res2['is_mapped']:
                                        self.quantifier.count_read(res2['candidates'])

            f1.close()
            if f2: f2.close()
            
            if self.quantifier and self.rpkm_out:
                self.quantifier.write_rpkm(self.rpkm_out, 1000000) 
            logger.info(f"Finished alignment.")
        except Exception as e:
            logger.error(f"Error processing FASTQ: {e}")

    def _read_chunk(self, f1, f2, size):
        chunk = []
        for _ in range(size):
            try:
                l1_1 = f1.readline(); l1_2 = f1.readline(); l1_3 = f1.readline(); l1_4 = f1.readline()
                if not l1_1: break
                r1 = {'id': l1_1.strip().split()[0][1:], 'seq': l1_2.strip(), 'qual': l1_4.strip()}
                r2 = None
                if f2:
                    l2_1 = f2.readline(); l2_2 = f2.readline(); l2_3 = f2.readline(); l2_4 = f2.readline()
                    r2 = {'id': l2_1.strip().split()[0][1:], 'seq': l2_2.strip(), 'qual': l2_4.strip()}
                chunk.append((r1, r2))
            except Exception:
                break
        return chunk

    def _write_sam_record(self, f, r1, r2):
        if r2:
            flag1 = 0x1 | 0x40
            flag2 = 0x1 | 0x80
            if not r1['is_mapped']: flag1 |= 0x4
            if not r2['is_mapped']: flag2 |= 0x4
            if not r2['is_mapped']: flag1 |= 0x8
            if not r1['is_mapped']: flag2 |= 0x8
            if r1.get('is_rc'): flag1 |= 0x10
            if r2.get('is_rc'): flag2 |= 0x10
            if r2.get('is_rc'): flag1 |= 0x20
            if r1.get('is_rc'): flag2 |= 0x20
            self._write_line(f, r1, flag1, r2)
            self._write_line(f, r2, flag2, r1)
        else:
            flag = 0
            if not r1['is_mapped']: flag |= 0x4
            if r1.get('is_rc'): flag |= 0x10
            self._write_line(f, r1, flag, None)

    def _write_line(self, f, r, flag, mate):
        rname = r['rname']
        pos = r['pos']
        mapq = r['mapq']
        cigar = r['cigar']
        rnext = '*'
        pnext = 0
        tlen = 0
        if mate:
            if mate['is_mapped']:
                rnext = '=' if mate['rname'] == rname else mate['rname']
                pnext = mate['pos']
            else:
                rnext = '*'
                pnext = 0
        
        tags = f"AS:i:{r['score']}"
        if self.quantifier and r['is_mapped']:
            txs = self.quantifier.get_overlapping_transcripts(r['candidates'])
            if txs:
                # Limit tags to avoid huge lines
                if len(txs) > 10:
                    tags += f"\tTX:Z:{','.join(txs[:10])}..."
                else:
                    tags += f"\tTX:Z:{','.join(txs)}"
        
        f.write(f"{r['id']}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{r['seq']}\t{r['qual']}\t{tags}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fractal Neo4j Graph Aligner (Hilbert)")
    parser.add_argument("--fastq", required=True, help="Input FASTQ file R1")
    parser.add_argument("--fastq2", help="Input FASTQ file R2 (Optional)")
    parser.add_argument("--out", required=True, help="Output SAM file")
    parser.add_argument("--rpkm", help="Output RPKM file (Optional)")
    parser.add_argument("--roads", help="Transcript Roads Pickle (Optional)")
    parser.add_argument("--user", default="neo4j", help="Neo4j user")
    parser.add_argument("--password", default="password", help="Neo4j password")
    parser.add_argument("--k", type=int, default=31, help="K-mer size")
    parser.add_argument("--cores", type=int, default=4, help="Number of cores")
    
    args = parser.parse_args()
    
    aligner = Neo4jGraphAligner("bolt://localhost:7687", args.user, args.password, k=args.k, cores=args.cores, rpkm_out=args.rpkm, roads_file=args.roads)
    aligner.align_fastq(args.fastq, args.fastq2, args.out)
    
    if aligner.quantifier and args.rpkm:
        aligner.quantifier.write_rpkm(args.rpkm, 1000000) 
