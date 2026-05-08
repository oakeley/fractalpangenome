from neo4j import GraphDatabase
import sys

def verify_privacy():
    driver = GraphDatabase.driver("bolt://localhost:7687", auth=("neo4j", "password"))
    with driver.session() as session:
        # Check Nodes for 'genome' property
        res = session.run("MATCH (n:RoadNode) WHERE n.genome IS NOT NULL RETURN count(n) as c").single()
        node_leaks = res['c']
        
        # Check Edges for 'genomes' property
        res = session.run("MATCH ()-[r]->() WHERE r.genomes IS NOT NULL RETURN count(r) as c").single()
        edge_leaks = res['c']
        
        print(f"Node Leaks: {node_leaks}")
        print(f"Edge Leaks: {edge_leaks}")
        
    driver.close()

if __name__ == "__main__":
    verify_privacy()
