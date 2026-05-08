import sys
from neo4j import GraphDatabase

def verify():
    driver = GraphDatabase.driver("bolt://localhost:7687", auth=("neo4j", "password"))
    
    with driver.session() as session:
        # 1. Count Nodes
        res = session.run("MATCH (n:RoadNode) RETURN count(n) as c").single()
        print(f"RoadNodes: {res['c']}")
        
        # 2. Check Wormholes
        res = session.run("MATCH ()-[r:WORMHOLE]->() RETURN count(r) as c").single()
        print(f"Wormholes: {res['c']}")
        
        # 3. Check Spatial Tiles
        res = session.run("MATCH (t:SpatialTile) RETURN count(t) as c").single()
        print(f"SpatialTiles: {res['c']}")
        
        # 4. Check Annotations (if any - none yet)
        
        # 5. Check Frequencies (Explorer Logic)
        # Some edges should have freq > 1 if shared
        res = session.run("MATCH ()-[r:NEXT]->() WHERE r.freq > 1 RETURN count(r) as c").single()
        print(f"Shared Edges (Freq > 1): {res['c']}")

    driver.close()

if __name__ == "__main__":
    verify()
