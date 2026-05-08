import unittest
import numpy as np
import os
import sys
from pathlib import Path

# Add current directory to path to import the module
sys.path.append(os.getcwd())

from gene_finder_pangenome import PathAccumulator, Gene

class TestHighwayLogic(unittest.TestCase):
    def setUp(self):
        self.accumulator = PathAccumulator(size=16) # Small grid for testing

    def test_path_accumulation(self):
        # Create mock genes
        gene1 = Gene(
            gene_id="G1", gene_name="TestGene", genome_id="Genome1",
            contig_id="chr1", contig_name="chr1", start=1000, end=2000, strand="+",
            features=[{'feature_type': 'exon', 'start': 1000, 'end': 1500, 'strand': '+'}]
        )
        
        gene2 = Gene(
            gene_id="G2", gene_name="TestGene", genome_id="Genome2",
            contig_id="chr1", contig_name="chr1", start=1000, end=2000, strand="+",
            features=[{'feature_type': 'exon', 'start': 1000, 'end': 1500, 'strand': '+'}]
        )
        
        gene3 = Gene(
            gene_id="G3", gene_name="TestGene", genome_id="Genome3",
            contig_id="chr1", contig_name="chr1", start=3000, end=4000, strand="+",
            features=[{'feature_type': 'exon', 'start': 3000, 'end': 3500, 'strand': '+'}]
        )

        # Add paths
        self.accumulator.add_path(gene1)
        self.accumulator.add_path(gene2)
        self.accumulator.add_path(gene3)

        # Check max count
        # Gene 1 and 2 should overlap (count 2)
        # Gene 3 is different (count 1)
        
        print(f"Max count in grid: {self.accumulator.max_count}")
        print(f"Grid sum: {np.sum(self.accumulator.grid)}")
        
        self.assertTrue(self.accumulator.max_count >= 2, "Should have at least one cell with count 2")
        
        # Verify that some cells have 2 and some have 1 (or 0)
        unique_values = np.unique(self.accumulator.grid)
        print(f"Unique grid values: {unique_values}")
        
        self.assertIn(2.0, unique_values)
        self.assertIn(1.0, unique_values)

if __name__ == '__main__':
    unittest.main()
