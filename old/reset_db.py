from neo4j import GraphDatabase
import sys

# Configuration
URI = "bolt://localhost:7687"
AUTH = ("neo4j", "password")

def reset_db():
    try:
        driver = GraphDatabase.driver(URI, auth=AUTH)
        with driver.session() as session:
            print("Deleting all nodes and relationships...")
            session.run("MATCH (n) DETACH DELETE n")
            print("Database reset complete.")
            
            # Verify
            result = session.run("MATCH (n) RETURN count(n) as count")
            count = result.single()["count"]
            print(f"Remaining nodes: {count}")
            
        driver.close()
    except Exception as e:
        print(f"Error resetting database: {e}")
        sys.exit(1)

if __name__ == "__main__":
    reset_db()
