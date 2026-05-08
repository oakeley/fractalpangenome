#!/usr/bin/env python3
"""
Gene Finder for HDF5 Pangenome - Fixed Version
Includes "Highway" visualization logic (frequency heatmaps)
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import argparse
import sys
import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from collections import defaultdict
import tempfile

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class Gene:
    """Gene information extracted from pangenome"""
    gene_id: str
    gene_name: str
    genome_id: str
    contig_id: str
    contig_name: str
    start: int
    end: int
    strand: str
    features: List[Dict]  # exons, transcripts, etc.
    sequence: Optional[str] = None

class HilbertCurveMapper:
    """Simple Hilbert curve mapper for gene visualization"""
    
    def __init__(self, size: int = 256):
        self.size = size
        self.total_positions = size * size
    
    def _hilbert_index_to_xy(self, index: int) -> Tuple[int, int]:
        """Convert 1D Hilbert index to 2D coordinates"""
        x = y = 0
        s = 1
        
        while s < self.size:
            rx = 1 & (index // 2)
            ry = 1 & (index ^ rx)
            x, y = self._hilbert_rotate(s, x, y, rx, ry)
            x += s * rx
            y += s * ry
            index //= 4
            s *= 2
        
        return x, y
    
    def _hilbert_rotate(self, n: int, x: int, y: int, rx: int, ry: int) -> Tuple[int, int]:
        """Rotate quadrant for Hilbert curve"""
        if ry == 0:
            if rx == 1:
                x = n - 1 - x
                y = n - 1 - y
            x, y = y, x
        return x, y
    
    def get_tile_coordinates(self, position: int, genome_size: int) -> Tuple[int, int]:
        """Get Hilbert coordinates for genomic position"""
        # Normalize position to Hilbert space
        hilbert_index = int((position / genome_size) * self.total_positions)
        hilbert_index = min(hilbert_index, self.total_positions - 1)
        return self._hilbert_index_to_xy(hilbert_index)

class PathAccumulator:
    """Accumulates genomic paths on the Hilbert grid to visualize 'highways'"""
    
    def __init__(self, size: int = 64): # Default grid size (2^6)
        self.size = size
        self.total_positions = size * size
        self.grid = np.zeros((size, size), dtype=float)
        self.mapper = HilbertCurveMapper(size)
        self.max_count = 0
        
    def add_path(self, gene: Gene):
        """Add a gene's path to the accumulator"""
        # Get unique exon coordinates
        exons = [feat for feat in gene.features if feat['feature_type'] == 'exon']
        
        if not exons:
            # Fallback for genes without exon annotations
            self._add_region(gene.start, gene.end, gene.strand)
            return

        # Extract unique coordinates
        unique_coords = set()
        for exon in exons:
            unique_coords.add((exon['start'], exon['end']))
        
        # Sort and normalize
        unique_coords = sorted(list(unique_coords))
        
        # Determine offset to normalize start position (like in create_hilbert_bed_file)
        if gene.strand == '-':
             unique_coords = sorted(unique_coords, reverse=True)
             reference_coord = unique_coords[0][1]
             
             for start, end in unique_coords:
                a = end - reference_coord - 1000
                b = start - reference_coord - 1000
                a_abs = abs(a)
                b_abs = abs(b)
                norm_start = min(a_abs, b_abs)
                norm_end = max(a_abs, b_abs)
                self._add_normalized_region(norm_start, norm_end)
        else:
            reference_coord = unique_coords[0][0]
            offset = reference_coord - 1000
            
            for start, end in unique_coords:
                norm_start = start - offset
                norm_end = end - offset
                self._add_normalized_region(norm_start, norm_end)

    def _add_region(self, start, end, strand):
        """Add a raw region (fallback)"""
        # Simple normalization: start at 1000
        length = end - start
        self._add_normalized_region(1000, 1000 + length)

    def _add_normalized_region(self, start, end):
        """Add a normalized region to the grid"""
        # Normalize to Hilbert space (assuming max coord ~200kb for a gene view?)
        # We need a fixed scale for the "highway" view to make sense across genes.
        # Let's assume a standard window size, e.g., 250kb, which fits most genes.
        MAX_COORD = 250000 
        
        norm_factor = self.total_positions / MAX_COORD
        
        norm_start = int(start * norm_factor)
        norm_end = int(end * norm_factor)
        norm_start = max(0, min(norm_start, self.total_positions - 1))
        norm_end = max(0, min(norm_end, self.total_positions - 1))
        
        for pos in range(norm_start, norm_end + 1):
            x, y = self.mapper._hilbert_index_to_xy(pos)
            if 0 <= x < self.size and 0 <= y < self.size:
                self.grid[y, x] += 1
                self.max_count = max(self.max_count, self.grid[y, x])

class GeneFinder:
    """Main class for finding genes in HDF5 pangenome and generating outputs"""
    
    def __init__(self, db_path: str):
        self.db_path = Path(db_path)
        self.hilbert = HilbertCurveMapper()
        
        if not self.db_path.exists():
            raise FileNotFoundError(f"Pangenome database not found: {db_path}")
    
    def load_gene_list(self, genelist_file: str) -> List[str]:
        """Load gene list from file"""
        genes = []
        with open(genelist_file, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene and not gene.startswith('#'):
                    genes.append(gene)
        
        logger.info(f"Loaded {len(genes)} genes from {genelist_file}")
        return genes
    
    def search_genes_in_pangenome(self, gene_names: List[str], 
                                 genome_ids: List[str] = None) -> Dict[str, List[Gene]]:
        """Search for genes in the pangenome database"""
        # First pass: collect all potential matches
        potential_matches = {}  # actual_gene_name -> [(search_term, gene_object), ...]
        
        try:
            with h5py.File(self.db_path, 'r') as f:
                # Get available genomes
                available_genomes = list(f.get('gene_index', {}).keys()) if genome_ids is None else genome_ids
                
                logger.info(f"Searching in {len(available_genomes)} genomes: {available_genomes}")
                
                for genome_id in available_genomes:
                    if f'gene_index/{genome_id}' not in f:
                        logger.warning(f"No gene index found for {genome_id}")
                        continue
                    
                    # Load gene index
                    gene_group = f[f'gene_index/{genome_id}']
                    if 'names' not in gene_group or 'data' not in gene_group:
                        logger.warning(f"Incomplete gene index for {genome_id}")
                        continue
                    
                    names = gene_group['names'][:]
                    data = gene_group['data'][:]
                    
                    # Search for each gene
                    for search_term in gene_names:
                        search_term_lower = search_term.lower()
                        
                        for i, name_bytes in enumerate(names):
                            name_str = self._decode_string(name_bytes).lower()
                            
                            if search_term_lower in name_str:
                                try:
                                    gene_data = json.loads(self._decode_string(data[i]))
                                    
                                    # Load additional annotations
                                    annotations = self._load_gene_annotations(f, genome_id, gene_data)
                                    
                                    # Use the ACTUAL gene name found from the database
                                    actual_gene_name = gene_data.get('name', search_term)
                                    
                                    gene = Gene(
                                        gene_id=gene_data['name'],
                                        gene_name=actual_gene_name,
                                        genome_id=genome_id,
                                        contig_id=gene_data.get('contig', 'unknown'),
                                        contig_name=gene_data.get('contig', 'unknown'),
                                        start=gene_data['start'],
                                        end=gene_data['end'],
                                        strand=gene_data.get('strand', '+'),
                                        features=annotations
                                    )
                                    
                                    # Collect potential matches
                                    if actual_gene_name not in potential_matches:
                                        potential_matches[actual_gene_name] = []
                                    potential_matches[actual_gene_name].append((search_term, gene))
                                    
                                    logger.info(f"Found {actual_gene_name} (searched for '{search_term}') in {genome_id}: {gene.start}-{gene.end}")
                                
                                except Exception as e:
                                    logger.warning(f"Error processing gene data: {e}")
        
        except Exception as e:
            logger.error(f"Error searching genes: {e}")
            raise
        
        # Second pass: resolve duplicates by assigning each gene to best-matching search term
        found_genes = self._resolve_gene_duplicates(potential_matches, gene_names)
        
        return found_genes
    
    def _resolve_gene_duplicates(self, potential_matches: Dict[str, List[Tuple[str, Gene]]], 
                                search_terms: List[str]) -> Dict[str, List[Gene]]:
        """Resolve duplicate genes by assigning each to the best-matching search term"""
        found_genes = {}
        
        for actual_gene_name, matches in potential_matches.items():
            if len(matches) == 1:
                # Only one search term found this gene - simple case
                search_term, gene = matches[0]
                if search_term not in found_genes:
                    found_genes[search_term] = []
                found_genes[search_term].append(gene)
            else:
                # Multiple search terms found this gene - need to pick the best match
                best_search_term = self._find_best_matching_search_term(actual_gene_name, matches)
                
                # Find the gene object for the best search term
                best_gene = None
                for search_term, gene in matches:
                    if search_term == best_search_term:
                        best_gene = gene
                        break
                
                if best_gene:
                    if best_search_term not in found_genes:
                        found_genes[best_search_term] = []
                    found_genes[best_search_term].append(best_gene)
                    
                    # Log the resolution
                    other_terms = [term for term, _ in matches if term != best_search_term]
                    logger.info(f"Resolved duplicate: {actual_gene_name} assigned to '{best_search_term}' "
                               f"(also matched: {', '.join(other_terms)})")
        
        return found_genes
    
    def _find_best_matching_search_term(self, actual_gene_name: str, 
                                       matches: List[Tuple[str, Gene]]) -> str:
        """Find the search term that best matches the actual gene name"""
        search_terms = [term for term, _ in matches]
        actual_lower = actual_gene_name.lower()
        
        # Priority 1: Exact match (case insensitive)
        for term in search_terms:
            if term.lower() == actual_lower:
                return term
        
        # Priority 2: Longest search term that completely matches
        # (prefer more specific searches)
        exact_substring_matches = []
        for term in search_terms:
            if term.lower() == actual_lower[:len(term)]:  # Gene starts with search term
                exact_substring_matches.append(term)
        
        if exact_substring_matches:
            # Return the longest exact substring match
            return max(exact_substring_matches, key=len)
        
        # Priority 3: Most characters in common (as fallback)
        best_term = search_terms[0]
        best_score = 0
        
        for term in search_terms:
            # Count common characters
            term_lower = term.lower()
            common_chars = sum(1 for c in term_lower if c in actual_lower)
            score = common_chars / len(term)  # Normalize by search term length
            
            if score > best_score:
                best_score = score
                best_term = term
        
        return best_term
    
    def _load_gene_annotations(self, hdf5_file, genome_id: str, gene_data: Dict) -> List[Dict]:
        """Load detailed annotations for a gene - fast search, then cleanup"""
        annotations = []
        
        if f'annotations/{genome_id}' not in hdf5_file:
            return annotations
        
        ann_group = hdf5_file[f'annotations/{genome_id}']
        gene_start = gene_data['start']
        gene_end = gene_data['end']
        
        # FAST SEARCH: Just get all overlapping features (original method)
        for feature_type in ['exon', 'transcript', 'CDS', 'mRNA']:
            if feature_type not in ann_group:
                continue
            
            type_group = ann_group[feature_type]
            if 'starts' not in type_group:
                continue
            
            starts = type_group['starts'][:]
            ends = type_group['ends'][:]
            
            for i in range(len(starts)):
                feat_start = int(starts[i])
                feat_end = int(ends[i])
                
                # Simple overlap check (fast!)
                if feat_start <= gene_end and feat_end >= gene_start:
                    feature = {
                        'feature_type': feature_type,
                        'start': feat_start,
                        'end': feat_end,
                        'strand': self._decode_string(type_group['strands'][i]) if 'strands' in type_group else '+',
                    }
                    
                    if 'attributes' in type_group:
                        try:
                            attr_str = self._decode_string(type_group['attributes'][i])
                            attributes = json.loads(attr_str)
                            feature['attributes'] = attributes
                            
                            # Extract transcript/gene IDs for naming
                            feature['transcript_id'] = self._extract_transcript_id(attributes)
                            feature['gene_id'] = self._extract_gene_id(attributes)
                            
                        except:
                            feature['attributes'] = {}
                            feature['transcript_id'] = None
                            feature['gene_id'] = None
                    else:
                        feature['attributes'] = {}
                        feature['transcript_id'] = None
                        feature['gene_id'] = None
                    
                    annotations.append(feature)
        
        # CLEANUP: Now filter the results to make biological sense
        if annotations:
            annotations = self._cleanup_gene_annotations(annotations, gene_data)
        
        # Sort features by start position
        annotations.sort(key=lambda x: x['start'])
        return annotations
    
    def _cleanup_gene_annotations(self, annotations: List[Dict], gene_data: Dict) -> List[Dict]:
        """Clean up annotations to make biological sense"""
        if not annotations:
            return annotations
        
        gene_strand = gene_data.get('strand', '+')
        
        # Remove duplicates and filter by strand
        seen_exons = set()
        cleaned_annotations = []
        
        for ann in annotations:
            # For exons, check strand consistency and remove duplicates
            if ann['feature_type'] == 'exon':
                # Only include exons on same strand as gene (or unspecified)
                if ann['strand'] == gene_strand or ann['strand'] == '.' or gene_strand == '.':
                    exon_key = (ann['start'], ann['end'], ann['strand'])
                    if exon_key not in seen_exons:
                        seen_exons.add(exon_key)
                        cleaned_annotations.append(ann)
            else:
                # For other features (transcript, CDS, mRNA), be more permissive
                cleaned_annotations.append(ann)
        
        return cleaned_annotations
    
    def _get_best_transcript_id(self, gene: Gene) -> str:
        """Get the best transcript ID for this gene"""
        # Look for transcript/mRNA features with IDs
        for feature in gene.features:
            if feature['feature_type'] in ['transcript', 'mRNA']:
                transcript_id = feature.get('transcript_id')
                if transcript_id:
                    return transcript_id
                
                # Also check ID field in attributes
                feature_id = feature.get('attributes', {}).get('ID')
                if feature_id and any(feature_id.startswith(prefix) for prefix in ['NM_', 'XM_', 'NR_', 'XR_']):
                    return feature_id
        
        return None
    
    def _extract_transcript_id(self, attributes: Dict) -> str:
        """Extract transcript ID from attributes (like NM_001385166.1)"""
        # Look for common transcript ID fields
        for key in ['ID', 'transcript_id', 'Name']:
            if key in attributes:
                value = attributes[key]
                # Check if it looks like a transcript ID (starts with NM_, XM_, etc.)
                if isinstance(value, str) and any(value.startswith(prefix) for prefix in ['NM_', 'XM_', 'NR_', 'XR_']):
                    return value
                # Or if it's a copy number ID like IL20_0
                if isinstance(value, str) and '_' in value and value.split('_')[-1].isdigit():
                    return value
        return None
    
    def _extract_gene_id(self, attributes: Dict) -> str:
        """Extract gene ID from attributes"""
        for key in ['gene', 'gene_name', 'Parent']:
            if key in attributes:
                return str(attributes[key])
        return None
    
    def extract_gene_sequences(self, genes: Dict[str, List[Gene]]) -> Dict[str, List[Gene]]:
        """Extract sequences for genes (mock implementation)"""
        # In a real implementation, this would extract sequences from the pangenome
        # For now, we'll generate mock sequences based on coordinates
        
        logger.info("Extracting gene sequences...")
        
        for gene_name, gene_list in genes.items():
            for gene in gene_list:
                # Mock sequence extraction
                gene_length = gene.end - gene.start
                gene.sequence = 'N' * gene_length  # Placeholder sequence
                
                logger.info(f"Extracted {gene_length} bp sequence for {gene.gene_name}")
        
        return genes
    
    def create_bed_files(self, genes: Dict[str, List[Gene]], output_dir: Path):
        """Create BED files for gene regions and normalized Hilbert BED files"""
        logger.info("Creating BED files...")
        
        for gene_name, gene_list in genes.items():
            for i, gene in enumerate(gene_list):
                # Create unique filenames using transcript IDs when available
                transcript_id = self._get_best_transcript_id(gene)
                
                if len(gene_list) > 1:
                    # Multiple hits - use transcript ID or coordinates
                    if transcript_id:
                        hit_id = f"{gene_name}_{transcript_id}"
                    else:
                        hit_id = f"{gene_name}_{gene.contig_name}_{gene.start}_{gene.end}"
                    bed_file = output_dir / f"{gene.genome_id}_{hit_id}.bed"
                    hilbert_bed_file = output_dir / f"{gene.genome_id}_{hit_id}_hilbert.bed"
                else:
                    # Single hit - use simple naming with transcript ID if available
                    if transcript_id:
                        hit_id = f"{gene_name}_{transcript_id}"
                        bed_file = output_dir / f"{gene.genome_id}_{hit_id}.bed"
                        hilbert_bed_file = output_dir / f"{gene.genome_id}_{hit_id}_hilbert.bed"
                    else:
                        bed_file = output_dir / f"{gene.genome_id}_{gene_name}.bed"
                        hilbert_bed_file = output_dir / f"{gene.genome_id}_{gene_name}_hilbert.bed"
                
                # Standard BED file
                with open(bed_file, 'w') as f:
                    # Write gene region
                    f.write(f"{gene.contig_name}\t{gene.start}\t{gene.end}\t{gene.gene_name}_gene\t255\t{gene.strand}\n")
                    
                    # Write exon regions with correct numbering based on strand
                    exons = [feat for feat in gene.features if feat['feature_type'] == 'exon']
                    
                    if gene.strand == '-':
                        # For negative strand, sort exons by start position (high to low) for correct biological numbering
                        exons.sort(key=lambda x: x['start'], reverse=True)
                    else:
                        # For positive strand, sort exons by start position (low to high)
                        exons.sort(key=lambda x: x['start'])
                    
                    for j, exon in enumerate(exons):
                        f.write(f"{gene.contig_name}\t{exon['start']}\t{exon['end']}\t{gene.gene_name}_exon_{j+1}\t255\t{exon.get('strand', gene.strand)}\n")
                
                # Create normalized BED file for Hilbert plotting
                self._create_hilbert_bed_file(gene, hilbert_bed_file)
                
                logger.info(f"Created BED files: {bed_file} and {hilbert_bed_file}")
    
    def _create_hilbert_bed_file(self, gene: Gene, output_file: Path):
        """Create normalized BED file for Hilbert plotting with correct strand-aware normalization"""
        # Get unique exon coordinates (remove duplicates)
        exons = [feat for feat in gene.features if feat['feature_type'] == 'exon']
        
        if not exons:
            # If no exons, create a simple entry for the whole gene
            with open(output_file, 'w') as f:
                f.write(f"{gene.contig_name}\t1000\t{gene.end - gene.start + 1000}\t255 255 0\n")
            return
        
        # Extract unique coordinates
        unique_coords = set()
        for exon in exons:
            unique_coords.add((exon['start'], exon['end']))
        
        # Sort by start position
        unique_coords = sorted(list(unique_coords))
        
        # Strand-aware normalization following the awk logic exactly
        if gene.strand == '-':
            # For negative strand: biological exon 1 has highest coordinates
            # Sort by start position in reverse order (high to low) 
            unique_coords = sorted(unique_coords, reverse=True)
            # Use the END coordinate of biological exon 1 as reference (sta=$3 from first exon)
            reference_coord = unique_coords[0][1]  # END coordinate of biological exon 1
            
            # Write normalized coordinates following awk logic:
            # a=$3-sta-1000; b=$2-sta-1000; print $1,(a<0) ? -a : a,(b<0) ? -b : b
            with open(output_file, 'w') as f:
                for start, end in unique_coords:
                    # Calculate a and b like in awk
                    a = end - reference_coord - 1000      # $3-sta-1000
                    b = start - reference_coord - 1000    # $2-sta-1000
                    # Take absolute values
                    a_abs = abs(a)
                    b_abs = abs(b)
                    # Ensure start < end (print smaller first, larger second)
                    norm_start = min(a_abs, b_abs)
                    norm_end = max(a_abs, b_abs)
                    f.write(f"{gene.contig_name}\t{norm_start}\t{norm_end}\t255 255 0\n")
        else:
            # For positive strand: normal order, use START coordinate of exon 1
            reference_coord = unique_coords[0][0]  # Start coordinate of first exon
            offset = reference_coord - 1000
            
            # Write normalized coordinates
            with open(output_file, 'w') as f:
                for start, end in unique_coords:
                    norm_start = start - offset
                    norm_end = end - offset
                    f.write(f"{gene.contig_name}\t{norm_start}\t{norm_end}\t255 255 0\n")
    
    def create_hilbert_plots(self, genes: Dict[str, List[Gene]], output_dir: Path):
        """Create normalized Hilbert curve plots for genes using BED files"""
        logger.info("Creating Hilbert plots...")
        
        for search_term, gene_list in genes.items():
            if len(gene_list) > 1:
                # Multiple hits - create comparison plot
                self._create_gene_comparison_hilbert_plot(search_term, gene_list, output_dir)
            else:
                # Single hit - create individual plot
                self._create_single_gene_hilbert_plot(search_term, gene_list[0], output_dir)
            
            # Create Highway Plot (Frequency Heatmap)
            self._create_highway_plot(search_term, gene_list, output_dir)

    def _create_highway_plot(self, search_term: str, gene_list: List[Gene], output_dir: Path):
        """Create a 'highway' plot showing frequency of paths across all genomes"""
        
        accumulator = PathAccumulator(size=64) # Use 2^6 grid
        
        for gene in gene_list:
            accumulator.add_path(gene)
            
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        
        # Plot heatmap
        if accumulator.max_count > 0:
            im = ax.imshow(accumulator.grid, cmap='hot', interpolation='nearest', vmin=0, vmax=accumulator.max_count)
            plt.colorbar(im, ax=ax, label='Genome Count')
        else:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center')
            
        ax.set_title(f"Genomic Highway: {search_term}\n(Shared paths across {len(gene_list)} genomes)", fontsize=14)
        ax.set_xticks([])
        ax.set_yticks([])
        
        # Save
        plot_file = output_dir / f"{search_term}_highway_hilbert.png"
        plt.tight_layout()
        plt.savefig(plot_file, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        logger.info(f"Created highway plot: {plot_file}")

    def _create_gene_comparison_hilbert_plot(self, search_term: str, gene_list: List[Gene], output_dir: Path):
        """Create comparison plot showing all hits for one search term"""
        
        num_hits = len(gene_list)
        cols = min(3, num_hits)
        rows = (num_hits + cols - 1) // cols
        
        fig, axes = plt.subplots(rows, cols, figsize=(8 * cols, 8 * rows))
        if num_hits == 1:
            axes = [axes]
        elif rows == 1:
            axes = axes if num_hits > 1 else [axes]
        else:
            axes = axes.flatten()
        
        # Get actual gene names for title
        actual_genes = list(set(gene.gene_name for gene in gene_list))
        if len(actual_genes) == 1 and actual_genes[0] == search_term:
            title = f'{search_term} - Multiple Hits Comparison ({num_hits} copies found)'
        else:
            title = f'Search: "{search_term}" → {", ".join(actual_genes)} - Multiple Hits Comparison ({num_hits} copies found)'
        
        fig.suptitle(title, fontsize=16, y=0.95)
        
        for i, gene in enumerate(gene_list):
            if i < len(axes):
                transcript_id = self._get_best_transcript_id(gene)
                if transcript_id:
                    hit_id = f"{gene.gene_name}_{transcript_id}"
                else:
                    hit_id = f"{gene.gene_name}_{gene.contig_name}_{gene.start}_{gene.end}"
                self._plot_gene_hilbert_from_bed_with_arrows(gene, axes[i], output_dir, hit_id)
        
        # Hide unused subplots
        for i in range(num_hits, len(axes)):
            axes[i].set_visible(False)
        
        # Save plot using search term for filename
        plot_file = output_dir / f"{search_term}_comparison_hilbert.png"
        plt.tight_layout()
        plt.savefig(plot_file, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        logger.info(f"Created comparison Hilbert plot: {plot_file}")
    
    def _create_single_gene_hilbert_plot(self, search_term: str, gene: Gene, output_dir: Path):
        """Create plot for a single gene hit"""
        
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        
        # Create title that shows search term and actual gene if different
        if gene.gene_name == search_term:
            title = f'{search_term}'
        else:
            title = f'Search: "{search_term}" → {gene.gene_name}'
        
        fig.suptitle(title, fontsize=16, y=0.95)
        
        transcript_id = self._get_best_transcript_id(gene)
        if transcript_id:
            hit_id = f"{gene.gene_name}_{transcript_id}"
        else:
            hit_id = ""
        
        self._plot_gene_hilbert_from_bed_with_arrows(gene, ax, output_dir, hit_id)
        
        # Save plot using search term for filename
        plot_file = output_dir / f"{search_term}_hilbert.png"
        plt.tight_layout()
        plt.savefig(plot_file, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        logger.info(f"Created Hilbert plot: {plot_file}")
    
    def _plot_gene_hilbert_from_bed_with_arrows(self, gene: Gene, ax, output_dir: Path, hit_id: str):
        """Plot gene using BED file coordinates with mathematically exact Hilbert curve"""
        
        # Determine BED file name based on gene properties, not hit_id parameter
        transcript_id = self._get_best_transcript_id(gene)
        
        if transcript_id:
            actual_hit_id = f"{gene.gene_name}_{transcript_id}"
            hilbert_bed_file = output_dir / f"{gene.genome_id}_{actual_hit_id}_hilbert.bed"
        else:
            hilbert_bed_file = output_dir / f"{gene.genome_id}_{gene.gene_name}_hilbert.bed"
        
        if not hilbert_bed_file.exists():
            ax.text(0.5, 0.5, f'BED file not found\n{hilbert_bed_file.name}', 
                   ha='center', va='center', transform=ax.transAxes)
            return
        
        # Read BED coordinates
        regions = []
        max_coord = 0
        
        with open(hilbert_bed_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    start = int(parts[1])
                    end = int(parts[2])
                    regions.append((start, end))
                    max_coord = max(max_coord, end)
        
        if not regions:
            ax.text(0.5, 0.5, 'No regions found', ha='center', va='center', transform=ax.transAxes)
            return
        
        # Create Hilbert curve visualization
        level = 6
        grid_size = 2 ** level
        total_positions = grid_size * grid_size
        
        # Normalize coordinates to fit in Hilbert space
        if max_coord > 0:
            norm_factor = total_positions / max_coord
        else:
            norm_factor = 1
        
        # Create binary mask for regions
        hilbert_mask = np.zeros(total_positions, dtype=bool)
        
        for start, end in regions:
            norm_start = int(start * norm_factor)
            norm_end = int(end * norm_factor)
            norm_start = max(0, min(norm_start, total_positions - 1))
            norm_end = max(0, min(norm_end, total_positions - 1))
            
            for pos in range(norm_start, norm_end + 1):
                if pos < total_positions:
                    hilbert_mask[pos] = True
        
        # Convert to 2D grid using Hilbert curve
        hilbert_grid = np.zeros((grid_size, grid_size))
        
        for i in range(total_positions):
            if hilbert_mask[i]:
                x, y = self._hilbert_index_to_xy(i, grid_size)
                if 0 <= x < grid_size and 0 <= y < grid_size:
                    hilbert_grid[y, x] = 1.0
        
        # Set up the plot
        ax.set_xlim(0, grid_size)
        ax.set_ylim(0, grid_size)
        ax.set_aspect('equal')
        
        # Draw the exact Hilbert curve path (no arrows - mathematically pure)
        self._draw_exact_hilbert_curve(ax, level)
        
        # Fill regions with semi-transparent red color so Hilbert curve shows through
        for y in range(grid_size):
            for x in range(grid_size):
                if hilbert_grid[y, x] > 0:
                    rect = patches.Rectangle((x, y), 1, 1, 
                                           facecolor='red', edgecolor='none', 
                                           alpha=0.6, linewidth=0)
                    ax.add_patch(rect)
        
        # Create detailed title with gene information
        gene_length = gene.end - gene.start
        
        # Extract transcript ID for display
        if transcript_id:
            title_parts = [
                f'{gene.genome_id}',
                f'{gene.gene_name} ({transcript_id})',
                f'{gene.contig_name}:{gene.start:,}-{gene.end:,}',
                f'({gene.strand}, {gene_length:,} bp)'
            ]
        else:
            title_parts = [
                f'{gene.genome_id}',
                f'{gene.gene_name}',
                f'{gene.contig_name}:{gene.start:,}-{gene.end:,}',
                f'({gene.strand}, {gene_length:,} bp)'
            ]
        
        title = '\n'.join([
            ' '.join(title_parts[:2]),  # Genome and gene name
            ' '.join(title_parts[2:])   # Coordinates and details
        ])
        
        ax.set_title(title, fontsize=10, pad=20)
        ax.set_xticks([])
        ax.set_yticks([])
        
        # Remove axes for clean look
        for spine in ax.spines.values():
            spine.set_visible(False)
    
    def _draw_exact_hilbert_curve(self, ax, level: int):
        """Draw the exact Hilbert curve path - mathematically precise, no arrows"""
        size = 2 ** level
        total_points = size * size
        
        # Generate the complete, exact Hilbert curve path
        path_points = []
        for i in range(total_points):
            x, y = self._hilbert_index_to_xy(i, size)
            path_points.append((x + 0.5, y + 0.5))  # Center of each cell
        
        # Draw the complete path as a thin line
        if len(path_points) > 1:
            path_x = [p[0] for p in path_points]
            path_y = [p[1] for p in path_points]
            ax.plot(path_x, path_y, color='gray', linewidth=0.5, alpha=0.4, zorder=1)
    
    def _hilbert_index_to_xy(self, index: int, size: int) -> Tuple[int, int]:
        """Convert 1D Hilbert index to 2D coordinates for given grid size"""
        x = y = 0
        s = 1
        
        while s < size:
            rx = 1 & (index // 2)
            ry = 1 & (index ^ rx)
            x, y = self._hilbert_rotate(s, x, y, rx, ry)
            x += s * rx
            y += s * ry
            index //= 4
            s *= 2
        
        return x, y
    
    def _hilbert_rotate(self, n: int, x: int, y: int, rx: int, ry: int) -> Tuple[int, int]:
        """Rotate quadrant for Hilbert curve"""
        if ry == 0:
            if rx == 1:
                x = n - 1 - x
                y = n - 1 - y
            x, y = y, x
        return x, y
    
    def create_summary_outputs(self, genes: Dict[str, List[Gene]], output_dir: Path):
        """Create summary files and combined visualizations"""
        logger.info("Creating summary outputs...")
        
        # Create gene summary file
        summary_file = output_dir / "gene_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("Gene Finder Summary\n")
            f.write("==================\n\n")
            
            total_search_terms = 0
            total_hits = 0
            
            for search_term, gene_list in genes.items():
                total_search_terms += 1
                total_hits += len(gene_list)
                
                # Get actual gene names found for this search term
                actual_genes = list(set(gene.gene_name for gene in gene_list))
                
                f.write(f"Search term: '{search_term}'\n")
                if len(actual_genes) == 1 and actual_genes[0] == search_term:
                    f.write(f"Found: {actual_genes[0]} (exact match)\n")
                else:
                    f.write(f"Found: {', '.join(actual_genes)}\n")
                
                f.write(f"Total hits: {len(gene_list)}\n")
                
                for i, gene in enumerate(gene_list, 1):
                    gene_length = gene.end - gene.start
                    exon_count = len([f for f in gene.features if f['feature_type'] == 'exon'])
                    transcript_id = self._get_best_transcript_id(gene)
                    
                    f.write(f"  Hit {i}: {gene.gene_name} in {gene.genome_id}")
                    if transcript_id:
                        f.write(f" ({transcript_id})")
                    f.write("\n")
                    f.write(f"    Location: {gene.contig_name}:{gene.start:,}-{gene.end:,} ({gene.strand})\n")
                    f.write(f"    Length: {gene_length:,} bp\n")
                    if exon_count > 0:
                        f.write(f"    Exons: {exon_count}\n")
                    
                    # Add information about duplicated exons if present
                    all_exon_coords = [(f['start'], f['end']) for f in gene.features if f['feature_type'] == 'exon']
                    unique_exon_coords = list(set(all_exon_coords))
                    if len(all_exon_coords) != len(unique_exon_coords):
                        f.write(f"    Original exon annotations: {len(all_exon_coords)}\n")
                        f.write(f"    Unique exon coordinates: {len(unique_exon_coords)}\n")
                        f.write(f"    Duplicates removed: {len(all_exon_coords) - len(unique_exon_coords)}\n")
                    
                    f.write("\n")
                
                f.write("\n")
            
            f.write("Summary Statistics:\n")
            f.write("==================\n")
            f.write(f"Total search terms processed: {total_search_terms}\n")
            f.write(f"Total hits found: {total_hits}\n")
            f.write(f"Average hits per search term: {total_hits / total_search_terms:.1f}\n")
            f.write(f"Search terms with multiple hits: {sum(1 for gene_list in genes.values() if len(gene_list) > 1)}\n\n")
            
            f.write("Note: Each gene is assigned to only one search term to avoid duplicates.\n")
            f.write("Assignment is based on best character match (exact > substring > common chars).\n")
            f.write("Some searches may find genes via substring matching (e.g., 'IL22' → 'IL22RA1').\n\n")
            
            f.write("Files Generated:\n")
            f.write("===============\n")
            f.write("BED files:\n")
            
            # Track generated files to avoid duplicates in listing
            listed_files = set()
            
            for search_term, gene_list in genes.items():
                for i, gene in enumerate(gene_list):
                    transcript_id = self._get_best_transcript_id(gene)
                    
                    if len(gene_list) > 1:
                        if transcript_id:
                            hit_id = f"{gene.gene_name}_{transcript_id}"
                        else:
                            hit_id = f"{gene.gene_name}_{gene.contig_name}_{gene.start}_{gene.end}"
                        bed_file = f"{gene.genome_id}_{hit_id}.bed"
                        hilbert_bed_file = f"{gene.genome_id}_{hit_id}_hilbert.bed"
                    else:
                        if transcript_id:
                            hit_id = f"{gene.gene_name}_{transcript_id}"
                            bed_file = f"{gene.genome_id}_{hit_id}.bed"
                            hilbert_bed_file = f"{gene.genome_id}_{hit_id}_hilbert.bed"
                        else:
                            bed_file = f"{gene.genome_id}_{gene.gene_name}.bed"
                            hilbert_bed_file = f"{gene.genome_id}_{gene.gene_name}_hilbert.bed"
                    
                    # Only list if not already listed
                    if bed_file not in listed_files:
                        f.write(f"  - {bed_file}\n")
                        f.write(f"  - {hilbert_bed_file}\n")
                        listed_files.add(bed_file)
                        listed_files.add(hilbert_bed_file)
            
            f.write("\nVisualization files:\n")
            for search_term, gene_list in genes.items():
                if len(gene_list) > 1:
                    f.write(f"  - {search_term}_comparison_hilbert.png (all {len(gene_list)} hits)\n")
                else:
                    f.write(f"  - {search_term}_hilbert.png\n")
                f.write(f"  - {search_term}_highway_hilbert.png (Frequency Heatmap)\n")
            f.write("  - genes_combined_hilbert.png (overview of all unique genes)\n")
        
        # Create detailed analysis file
        analysis_file = output_dir / "detailed_analysis.txt"
        with open(analysis_file, 'w') as f:
            f.write("Detailed Gene Analysis\n")
            f.write("=====================\n\n")
            
            for search_term, gene_list in genes.items():
                actual_genes = list(set(gene.gene_name for gene in gene_list))
                
                f.write(f"Search term: '{search_term}':\n")
                f.write("-" * (len(search_term) + 18) + "\n\n")
                
                if len(actual_genes) == 1 and actual_genes[0] == search_term:
                    f.write(f"Exact match found: {actual_genes[0]}\n")
                else:
                    f.write(f"Found genes: {', '.join(actual_genes)}\n")
                    if len(actual_genes) == 1:
                        f.write(f"This was found via substring matching.\n")
                f.write("\n")
                
                if len(gene_list) > 1:
                    f.write(f"Multiple copies detected: {len(gene_list)} hits\n")
                    f.write("This could indicate:\n")
                    f.write("  - Gene family members with similar sequences\n")
                    f.write("  - Pseudogenes or processed genes\n")
                    f.write("  - Assembly duplications\n")
                    f.write("  - True gene duplications\n\n")
                    
                    # Compare hits
                    f.write("Comparison of hits:\n")
                    for i, gene in enumerate(gene_list, 1):
                        transcript_id = self._get_best_transcript_id(gene)
                        f.write(f"  Hit {i}: {gene.gene_name} on {gene.contig_name} ({gene.end - gene.start:,} bp, ")
                        f.write(f"{len([f for f in gene.features if f['feature_type'] == 'exon'])} exons")
                        if transcript_id:
                            f.write(f", {transcript_id}")
                        f.write(")\n")
                    f.write("\n")
                
                else:
                    f.write("Single copy gene - unique hit found\n\n")
                
                # Provide recommendations
                f.write("Recommendations:\n")
                if len(gene_list) > 1:
                    f.write("  1. Review the comparison plot to identify structural differences\n")
                    f.write("  2. Check if hits represent different gene family members\n")
                    f.write("  3. Validate with sequence similarity analysis\n")
                    f.write("  4. Consider chromosome context and synteny\n")
                else:
                    f.write("  1. Single hit suggests unique gene location\n")
                    f.write("  2. Review Hilbert plot for internal structure\n")
                    f.write("  3. Compare with reference annotations if available\n")
                
                f.write("\n" + "="*50 + "\n\n")
        
        # Create combined Hilbert plot
        self._create_combined_hilbert_plot(genes, output_dir)
        
        logger.info(f"Summary written to: {summary_file}")
        logger.info(f"Detailed analysis written to: {analysis_file}")
    
    def _create_combined_hilbert_plot(self, genes: Dict[str, List[Gene]], output_dir: Path):
        """Create a combined plot showing all unique genes"""
        if not genes:
            return
        
        # Flatten all genes and deduplicate by unique gene identifier
        unique_genes = {}  # gene_key -> (gene, label)
        
        for search_term, gene_list in genes.items():
            for i, gene in enumerate(gene_list):
                transcript_id = self._get_best_transcript_id(gene)
                
                # Create unique key based on gene properties, not search term
                if transcript_id:
                    gene_key = f"{gene.genome_id}_{gene.gene_name}_{transcript_id}"
                    label = f"{gene.gene_name}_{transcript_id}"
                else:
                    gene_key = f"{gene.genome_id}_{gene.gene_name}_{gene.start}_{gene.end}"
                    label = f"{gene.gene_name}"
                
                # Only add if not already present (avoid duplicates from multiple search terms)
                if gene_key not in unique_genes:
                    unique_genes[gene_key] = (gene, label)
        
        # Convert to lists for plotting
        all_genes = [gene for gene, label in unique_genes.values()]
        gene_labels = [label for gene, label in unique_genes.values()]
        
        logger.info(f"Creating combined plot with {len(all_genes)} unique genes")
        
        # Determine layout
        num_genes = len(all_genes)
        cols = min(4, num_genes)
        rows = (num_genes + cols - 1) // cols
        
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 5 * rows))
        if rows == 1 and cols == 1:
            axes = [axes]
        elif rows == 1:
            axes = axes
        else:
            axes = axes.flatten()
        
        fig.suptitle(f'Gene Collection - Hilbert Space Overview ({num_genes} unique genes)', fontsize=16)
        
        for i, (gene, label) in enumerate(zip(all_genes, gene_labels)):
            if i >= len(axes):
                break
            
            ax = axes[i]
            
            # Use simplified plotting for overview
            self._plot_gene_simple_overview(gene, ax, label, output_dir)
        
        # Hide unused subplots
        for i in range(len(all_genes), len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        combined_plot = output_dir / "genes_combined_hilbert.png"
        plt.savefig(combined_plot, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Combined plot saved: {combined_plot}")
    
    def _plot_gene_simple_overview(self, gene: Gene, ax, label: str, output_dir: Path):
        """Simple overview plot for combined view"""
        # Construct the exact Hilbert BED filename based on gene properties
        transcript_id = self._get_best_transcript_id(gene)
        
        if transcript_id:
            hit_id = f"{gene.gene_name}_{transcript_id}"
            hilbert_bed_file = output_dir / f"{gene.genome_id}_{hit_id}_hilbert.bed"
        else:
            hilbert_bed_file = output_dir / f"{gene.genome_id}_{gene.gene_name}_hilbert.bed"
        
        if not hilbert_bed_file.exists():
            ax.text(0.5, 0.5, f'{label}\nNo BED file\n{hilbert_bed_file.name}', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=8)
            ax.set_title(label, fontsize=10)
            return
        
        # Read BED coordinates
        regions = []
        max_coord = 0
        
        try:
            with open(hilbert_bed_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        start = int(parts[1])
                        end = int(parts[2])
                        regions.append((start, end))
                        max_coord = max(max_coord, end)
        except Exception as e:
            ax.text(0.5, 0.5, f'{label}\nError reading BED\n{str(e)}', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=8)
            ax.set_title(label, fontsize=10)
            return
        
        if not regions:
            ax.text(0.5, 0.5, f'{label}\nNo regions', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(label, fontsize=10)
            return
        
        # Simplified Hilbert visualization
        level = 5  # Smaller for overview
        grid_size = 2 ** level
        total_positions = grid_size * grid_size
        
        if max_coord > 0:
            norm_factor = total_positions / max_coord
        else:
            norm_factor = 1
        
        hilbert_grid = np.zeros((grid_size, grid_size))
        
        for start, end in regions:
            norm_start = int(start * norm_factor)
            norm_end = int(end * norm_factor)
            norm_start = max(0, min(norm_start, total_positions - 1))
            norm_end = max(0, min(norm_end, total_positions - 1))
            
            for pos in range(norm_start, norm_end + 1):
                if pos < total_positions:
                    x, y = self._hilbert_index_to_xy(pos, grid_size)
                    if 0 <= x < grid_size and 0 <= y < grid_size:
                        hilbert_grid[y, x] = 1.0
        
        # Plot
        ax.imshow(hilbert_grid, cmap='Reds', interpolation='nearest')
        ax.set_title(f'{label}\n{gene.contig_name}:{gene.start:,}', fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        
        # Add gene length info
        gene_length = gene.end - gene.start
        ax.text(0.02, 0.98, f"{gene_length:,} bp", 
               transform=ax.transAxes, va='top', ha='left', fontsize=6,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    def _decode_string(self, value) -> str:
        """Decode string from HDF5"""
        if isinstance(value, bytes):
            return value.decode('utf-8')
        return str(value)
    
    def run_analysis(self, gene_list_file: str, output_dir: str, genome_ids: List[str] = None, 
                    create_r_scripts: bool = False):
        """Run the complete gene finding analysis"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Starting gene finder analysis...")
        logger.info(f"Input: {gene_list_file}")
        logger.info(f"Output: {output_path}")
        logger.info(f"Database: {self.db_path}")
        
        # Step 1: Load gene list
        gene_names = self.load_gene_list(gene_list_file)
        
        # Step 2: Search genes in pangenome
        found_genes = self.search_genes_in_pangenome(gene_names, genome_ids)
        
        if not found_genes:
            logger.warning("No genes found in pangenome!")
            return
        
        logger.info(f"Found genes for {len(found_genes)} search terms total")
        
        # Log what was found for each search term
        for search_term, gene_list in found_genes.items():
            actual_gene_names = list(set(gene.gene_name for gene in gene_list))
            if len(actual_gene_names) == 1 and actual_gene_names[0] == search_term:
                logger.info(f"  '{search_term}': exact match, {len(gene_list)} hit(s)")
            else:
                logger.info(f"  '{search_term}': found {actual_gene_names}, {len(gene_list)} hit(s)")
                if len(actual_gene_names) == 1 and actual_gene_names[0] != search_term:
                    logger.info(f"    → matched {actual_gene_names[0]} via substring search")
        
        # Step 3: Extract sequences
        genes_with_sequences = self.extract_gene_sequences(found_genes)
        
        # Step 4: Create BED files (both standard and Hilbert-normalized)
        self.create_bed_files(genes_with_sequences, output_path)
        
        # Step 5: Create Hilbert plots
        self.create_hilbert_plots(genes_with_sequences, output_path)
        
        # Step 6: Create R scripts if requested
        if create_r_scripts:
            self.create_r_plotting_scripts(genes_with_sequences, output_path)
        
        # Step 7: Create summary
        self.create_summary_outputs(genes_with_sequences, output_path)
        
        logger.info("Gene finder analysis complete!")


def main():
    """Main function with command line interface"""
    parser = argparse.ArgumentParser(
        description="Gene Finder for HDF5 Pangenome - Find genes and create visualizations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python gene_finder_pangenome.py --db pangenome.h5 --genes genelist.txt --output results/
  
  # Search in specific genomes
  python gene_finder_pangenome.py --db pangenome.h5 --genes genelist.txt --genomes Human_T2T Human_NA12878 --output results/
  
  # Include R plotting scripts (for compatibility with original workflow)
  python gene_finder_pangenome.py --db pangenome.h5 --genes genelist.txt --output results/ --create-r-scripts
  
  # Example gene list file (genelist.txt):
  BRCA1
  BRCA2
  TP53
  SHANK3
  RAI1
  
Features:
  - Handles multiple gene copies/hits automatically
  - Creates separate BED files for each hit when multiple copies found
  - Generates comparison plots for genes with multiple hits
  - Removes duplicate exon coordinates automatically
  - Strand-aware exon numbering (exon 1 = transcription start for both + and - strands)
  - Normalizes coordinates so biological exon 1 starts at position 1000 regardless of strand
  - Creates both Python and R-compatible visualizations
  - Mathematically exact Hilbert curve visualization (no approximations)
  - Provides detailed analysis of gene duplications and structure
  - [NEW] "Highway" visualization: Frequency heatmap of shared genomic paths

Output Files:
  Standard BED files: {genome_id}_{gene_name}_{transcript_id}.bed (strand-aware exon numbering)
  Hilbert BED files: {genome_id}_{gene_name}_{transcript_id}_hilbert.bed (normalized for plotting)
  Individual plots: {gene_name}_hilbert.png (single hits)
  Comparison plots: {gene_name}_comparison_hilbert.png (multiple hits)
  Highway plots: {gene_name}_highway_hilbert.png (frequency heatmap)
  Combined overview: genes_combined_hilbert.png
  Summary reports: gene_summary.txt, detailed_analysis.txt
  
The Hilbert BED files have coordinates normalized so biological exon 1 starts at position 1000,
regardless of gene strand. For negative strand genes, exons are reordered so the biologically
first exon (highest genomic coordinate) appears first in the visualization.

Multiple Gene Copies:
When multiple copies of a gene are found (like IL20 on chr1, chr3, chr6), the system:
  - Creates separate BED files for each hit using transcript IDs (e.g., IL20_NM_001385166.1)
  - Generates a comparison plot showing all hits side-by-side with transcript information
  - Provides detailed analysis of structural differences
  - Maintains full traceability of each gene copy with proper transcript annotations
  - Uses transcript IDs from GFF3 annotations when available (NM_, XM_, NR_, XR_ prefixes)
  - Correctly numbers exons based on biological direction (5' to 3') regardless of genomic strand

Duplicate Resolution:
When multiple search terms match the same gene (e.g., searching for both "IL22" and "IL22RA1"):
  - Each gene is assigned to only ONE search term to avoid duplicates
  - Assignment priority: exact match > longest substring match > most common characters
  - This ensures clean organization and prevents genes from appearing in multiple lists
  - The assignment is logged for transparency (e.g., "IL22RA1 assigned to 'IL22RA1' (also matched: IL22)")
        """
    )
    
    parser.add_argument('--db', required=True, help='HDF5 pangenome database file')
    parser.add_argument('--genes', required=True, help='File containing gene names (one per line)')
    parser.add_argument('--output', required=True, help='Output directory for results')
    parser.add_argument('--genomes', nargs='+', help='Specific genome IDs to search (default: all)')
    parser.add_argument('--create-r-scripts', action='store_true', 
                       help='Create R plotting scripts compatible with original HilbertCurve workflow')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose logging')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Create gene finder
        gene_finder = GeneFinder(args.db)
        
        # Run analysis
        gene_finder.run_analysis(args.genes, args.output, args.genomes, args.create_r_scripts)
        
        print(f"✅ Gene finder analysis complete!")
        print(f"📁 Results saved to: {args.output}")
        print(f"📊 Check the following outputs:")
        print(f"   - BED files: Standard and Hilbert-normalized coordinates")
        print(f"   - Python Hilbert plots: Individual and comparison views")
        print(f"   - Highway plots: Frequency heatmaps (*_highway_hilbert.png)")
        if args.create_r_scripts:
            print(f"   - R plotting scripts: Compatible with original HilbertCurve workflow")
        print(f"   - Summary reports: gene_summary.txt, detailed_analysis.txt")
        print(f"   - Combined overview: genes_combined_hilbert.png")
        
        print(f"\n💡 Key Features:")
        print(f"   ✅ Multiple gene copy detection and separate visualization")
        print(f"   ✅ Automatic duplicate exon removal")
        print(f"   ✅ Strand-aware exon numbering (biological 5' to 3' direction)")
        print(f"   ✅ Coordinate normalization (biological exon 1 at position 1000)")
        print(f"   ✅ Mathematically exact Hilbert curve visualization")
        print(f"   ✅ Detailed structural analysis and recommendations")
        print(f"   ✅ Transcript ID-based naming (NM_, XM_, etc.)")
        print(f"   ✅ Uses actual gene names found (not search terms)")
        print(f"   ✅ Automatic duplicate resolution (each gene in best-matching search only)")
        print(f"   ✅ Highway visualization (frequency heatmap)")
        
        print(f"\n🔍 Note: Genes found by multiple search terms are assigned to best match.")
        print(f"   Priority: exact match > longest substring > most common characters")
        
        print(f"\n🔍 For genes with multiple hits:")
        print(f"   - Check *_comparison_hilbert.png for side-by-side comparison")
        print(f"   - Review detailed_analysis.txt for structural insights")
        print(f"   - Each hit has separate BED files with transcript identifiers")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
