#!/usr/bin/env python3
"""
Fractal Neo4j Graph Verifier
Performs comprehensive integrity checks on the Pangenome Graph.
"""

import sys
import logging
import argparse
from neo4j import GraphDatabase

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GraphVerifier:
    def __init__(self, uri, user, password):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))
        
    def close(self):
        self.driver.close()
        
    def run_checks(self):
        logger.info("Starting Graph Verification...")
        
        with self.driver.session() as session:
            # 1. General Stats
            self._log_stat(session, "Total Kmer Nodes", "MATCH (k:Kmer) RETURN count(k) as c")
            self._log_stat(session, "Total NEXT Edges", "MATCH ()-[r:NEXT]->() RETURN count(r) as c")
            self._log_stat(session, "Total Genomes", "MATCH (g:Genome) RETURN count(g) as c")
            
            # 2. Orphan Check (Kmers with NO edges)
            logger.info("Checking for Orphan Nodes (Kmers with no connections)...")
            result = session.run("""
                MATCH (k:Kmer)
                WHERE NOT (k)--()
                RETURN count(k) as c
            """)
            orphans = result.single()['c']
            if orphans > 0:
                logger.warning(f"FOUND {orphans} ORPHAN NODES! These are disconnected from everything.")
            else:
                logger.info("PASS: No orphan nodes found.")
                
            # 3. Dead End Check (Kmers with no outgoing NEXT)
            # Note: Chromosome ends are naturally dead ends.
            logger.info("Checking for Dead Ends (Kmers with no outgoing NEXT)...")
            result = session.run("""
                MATCH (k:Kmer)
                WHERE NOT (k)-[:NEXT]->()
                RETURN count(k) as c
            """)
            dead_ends = result.single()['c']
            logger.info(f"Found {dead_ends} dead ends.")
            # We expect at least one per chromosome/contig.
            
            # 4. Start Check (Kmers with no incoming NEXT)
            logger.info("Checking for Start Nodes (Kmers with no incoming NEXT)...")
            result = session.run("""
                MATCH (k:Kmer)
                WHERE NOT ()-[:NEXT]->(k)
                RETURN count(k) as c
            """)
            starts = result.single()['c']
            logger.info(f"Found {starts} start nodes.")
            
            # 5. Genome Consistency
            logger.info("Checking Genome Consistency...")
            result = session.run("""
                MATCH (g:Genome)
                WHERE NOT (g)-[:HAS_CONTIG]->()
                RETURN g.id as id
            """)
            empty_genomes = [r['id'] for r in result]
            if empty_genomes:
                logger.error(f"FAIL: Found Genomes with no Contigs: {empty_genomes}")
            else:
                logger.info("PASS: All Genomes have Contigs.")
                
            # 6. Contig Connectivity
            logger.info("Checking Contig Connectivity...")
            result = session.run("""
                MATCH (c:Contig)
                WHERE NOT (c)-[:STARTS_AT]->()
                RETURN c.id as id
            """)
            broken_contigs = [r['id'] for r in result]
            if broken_contigs:
                logger.error(f"FAIL: Found Contigs with no Start Kmer: {broken_contigs}")
            else:
                logger.info("PASS: All Contigs have a Start Kmer.")

            # 7. Journey Integrity (Sample Check)
            logger.info("Checking Journey Integrity (Sampling 100 edges)...")
            result = session.run("""
                MATCH ()-[r:NEXT]->()
                RETURN r.genomes as g
                LIMIT 100
            """)
            missing_journey = 0
            for r in result:
                if not r['g'] or not isinstance(r['g'], list):
                    missing_journey += 1
            
            if missing_journey > 0:
                logger.error(f"FAIL: Found {missing_journey}% of sampled edges missing 'genomes' property.")
            else:
                logger.info("PASS: Sampled edges have valid 'genomes' property.")

    def _log_stat(self, session, name, query):
        result = session.run(query)
        count = result.single()['c']
        logger.info(f"{name}: {count}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fractal Graph Verifier")
    parser.add_argument("--user", default="neo4j", help="Neo4j user")
    parser.add_argument("--password", default="password", help="Neo4j password")
    
    args = parser.parse_args()
    
    verifier = GraphVerifier("bolt://localhost:7687", args.user, args.password)
    verifier.run_checks()
    verifier.close()
