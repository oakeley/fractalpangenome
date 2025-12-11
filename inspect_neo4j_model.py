#!/usr/bin/env python3
"""
Inspect Neo4j Data Model
Verifies that the graph structure supports 'journeys' by checking edge properties.
"""

import sys
import logging
import argparse
from neo4j import GraphDatabase

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def inspect_model(uri, user, password):
    driver = GraphDatabase.driver(uri, auth=(user, password))
    
    try:
        with driver.session() as session:
            # 1. Check Node Count
            result = session.run("MATCH (k:Kmer) RETURN count(k) as count")
            node_count = result.single()['count']
            logger.info(f"Total Kmer Nodes: {node_count}")
            
            # 2. Check Edge Count
            result = session.run("MATCH ()-[r:NEXT]->() RETURN count(r) as count")
            edge_count = result.single()['count']
            logger.info(f"Total NEXT Edges: {edge_count}")
            
            # 3. Inspect a sample edge for 'genomes' property
            logger.info("Inspecting sample edge properties...")
            result = session.run("""
                MATCH ()-[r:NEXT]->() 
                RETURN r.count as count, r.genomes as genomes 
                LIMIT 5
            """)
            
            for record in result:
                logger.info(f"Edge: count={record['count']}, genomes={record['genomes']}")
                
                if record['genomes'] and isinstance(record['genomes'], list):
                    logger.info("SUCCESS: 'genomes' property is a list (supports journeys).")
                else:
                    logger.error("FAILURE: 'genomes' property is missing or not a list!")

    except Exception as e:
        logger.error(f"Inspection failed: {e}")
    finally:
        driver.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Inspect Neo4j Model")
    parser.add_argument("--user", default="neo4j", help="Neo4j user")
    parser.add_argument("--password", default="password", help="Neo4j password")
    
    args = parser.parse_args()
    
    inspect_model("bolt://localhost:7687", args.user, args.password)
