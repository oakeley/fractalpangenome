
import sys
from neo4j import GraphDatabase
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class HighwayInspector:
    def __init__(self, uri, user, password):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        self.driver.close()

    def inspect_kmer(self, kmer_seq):
        with self.driver.session() as session:
            # 1. Check Highway Status (Sharedness)
            # We look at outgoing edges to see which genomes traverse this node
            result = session.run("""
                MATCH (k:Kmer {seq: $seq})-[r:NEXT]->()
                RETURN r.genomes AS genomes, r.count AS count
            """, seq=kmer_seq)
            
            genomes_set = set()
            total_traversals = 0
            
            for record in result:
                g_list = record["genomes"]
                if g_list:
                    genomes_set.update(g_list)
                total_traversals += record["count"]
            
            is_highway = len(genomes_set) > 1
            
            print(f"--- K-mer Inspection: {kmer_seq[:10]}... ---")
            print(f"Found in Genomes: {list(genomes_set)}")
            print(f"Total Traversals: {total_traversals}")
            print(f"Highway Status: {'HIGHWAY' if is_highway else 'SIDE STREET'}")
            
            # 2. Check for Direct Annotations (Starts/Ends)
            result_anno = session.run("""
                MATCH (k:Kmer {seq: $seq})<-[:STARTS_AT|ENDS_AT]-(f:Feature)
                RETURN f.id AS id, f.name AS name, f.type AS type, f.genome AS genome
            """, seq=kmer_seq)
            
            print("\nDirect Annotations:")
            found_anno = False
            for record in result_anno:
                found_anno = True
                print(f" - [{record['genome']}] {record['type']}: {record['name']} ({record['id']})")
            
            if not found_anno:
                print(" - None")

    def find_highway_example(self):
        """
        Find a k-mer that is shared by multiple genomes (Highway).
        Strategy:
        1. Try highly conserved genes (GAPDH, ACTB).
        2. Fallback to direct database scan for shared edges.
        """
        with self.driver.session() as session:
            # 1. Try Conserved Genes
            for gene in ["GAPDH", "ACTB", "TP53"]:
                print(f"Searching for shared k-mer in gene: {gene}...")
                result = session.run("""
                    MATCH (f:Feature)
                    WHERE f.name = $name
                    WITH f
                    MATCH (f)-[:STARTS_AT]->(k:Kmer)
                    MATCH (k)-[r:NEXT]->()
                    WHERE size(r.genomes) > 1
                    RETURN k.seq AS seq, r.genomes AS genomes
                    LIMIT 1
                """, name=gene)
                
                record = result.single()
                if record:
                    print(f"Found shared start k-mer in {gene}! Shared by: {record['genomes']}")
                    return record["seq"]
            
            # 2. Fallback: Scan for ANY shared edge (limit scan)
            print("Specific genes not found or not shared at start. Scanning for any shared path...")
            result = session.run("""
                MATCH ()-[r:NEXT]->()
                WHERE size(r.genomes) > 1
                RETURN startNode(r).seq AS seq
                LIMIT 1
            """)
            record = result.single()
            if record:
                return record["seq"]
                
            return None

if __name__ == "__main__":
    inspector = HighwayInspector("bolt://localhost:7687", "neo4j", "password")
    
    # Try to find a real highway example
    print("Searching for a Highway k-mer...")
    highway_kmer = inspector.find_highway_example()
    
    if highway_kmer:
        inspector.inspect_kmer(highway_kmer)
    else:
        print("No highway k-mers found (yet). Try adding more genomes!")
        
    inspector.close()
