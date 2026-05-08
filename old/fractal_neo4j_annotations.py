#!/usr/bin/env python3
"""
Fractal Neo4j Annotations V2
Maps GFF features to RoadNodes.
- Increments RoadNode frequency for Exons.
- Creates SPLICE_NEXT edges for Junctions.
- Creates Feature nodes (Gene, Transcript) and links them.
"""

import sys
import argparse
import logging
import gzip
import csv
from neo4j import GraphDatabase

# Logging
sys.stdout.reconfigure(line_buffering=True)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class Neo4jAnnotator:
    def __init__(self, uri, user, password, k=31):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))
        self.k = k
        self.init_db()

    def close(self):
        self.driver.close()

    def init_db(self):
        with self.driver.session() as session:
            session.run("CREATE INDEX feature_id_idx IF NOT EXISTS FOR (f:Feature) ON (f.id)")
            logger.info("Annotation indexes intitialized.")

    def load_annotations(self, gff_path, genome_id):
        logger.info(f"Loading annotations for {genome_id} from {gff_path}...")
        
        try:
            if gff_path.endswith('.gz'):
                handle = gzip.open(gff_path, 'rt')
            else:
                handle = open(gff_path, 'r')
            reader = csv.reader(handle, delimiter='\t')
        except Exception as e:
            logger.error(f"Failed to open GFF: {e}")
            return

        batch_features = []
        transcripts = {} # tx_id -> list of exons

        count = 0
        BATCH_SIZE = 500
        
        with self.driver.session() as session:
            for row in reader:
                if not row or row[0].startswith('#'): continue
                if len(row) < 9: continue

                contig = row[0]
                ftype = row[2]
                try:
                    start = int(row[3]) - 1 # 0-indexed
                    end = int(row[4])
                except:
                    continue
                
                # Parse attributes
                attr_map = {}
                for attr in row[8].split(';'):
                    if '=' in attr:
                        k, v = attr.split('=', 1)
                        attr_map[k] = v
                
                fid = attr_map.get('ID')
                parent = attr_map.get('Parent')
                name = attr_map.get('Name', fid)
                
                if not fid: continue

                # 1. Genes / Transcripts -> Create Feature Node
                if ftype in ['gene', 'mRNA', 'transcript']:
                    # Create Node
                    # Link to Start/End RoadNodes
                    self._create_feature(session, fid, name, ftype, contig, start, end, genome_id)
                
                # 2. Exons -> Increment Usage / Store for Splices
                if ftype == 'exon' and parent:
                    # Increment Frequency of underlying path
                    self._increment_exon_path(session, contig, start, end, genome_id)
                    
                    if parent not in transcripts: transcripts[parent] = []
                    transcripts[parent].append({'start': start, 'end': end, 'contig': contig})
                
                count += 1
                if count % 1000 == 0:
                    sys.stdout.write(f"Processed {count} features...\r")
                    sys.stdout.flush()

            # Process Splice Junctions
            logger.info("\nProcessing Splice Junctions...")
            junction_buffer = []
            
            for tx_id, exons in transcripts.items():
                exons.sort(key=lambda x: x['start'])
                if len(exons) < 2: continue
                
                for i in range(len(exons) - 1):
                    e1 = exons[i]
                    e2 = exons[i+1]
                    
                    # Splice: End of E1 -> Start of E2
                    # Nodes: E1_End_Node -> E2_Start_Node
                    # E1 End Node pos: e1['end'] - K
                    # E2 Start Node pos: e2['start']
                    
                    pos1 = e1['end'] - self.k
                    pos2 = e2['start']
                    
                    junction_buffer.append({
                        'c1': e1['contig'], 'p1': pos1,
                        'c2': e2['contig'], 'p2': pos2,
                        'gid': genome_id,
                        'tx': tx_id
                    })
                    
                    if len(junction_buffer) >= BATCH_SIZE:
                        self._flush_junctions(session, junction_buffer)
                        junction_buffer = []
            
            if junction_buffer:
                self._flush_junctions(session, junction_buffer)

        logger.info("Annotations loaded.")

    def _create_feature(self, session, fid, name, ftype, contig, start, end, genome_id):
        # Create Feature Node
        session.run("""
            MERGE (f:Feature {id: $id})
            SET f.name = $name, f.type = $type, f.genome = $gid
        """, id=fid, name=name, type=ftype, gid=genome_id)
        
        # Link to Road (Approximate Start/End)
        # Note: We link to the K-mer RoadNode at start, and one near end
        node_start_pos = start
        node_end_pos = end - self.k
        
        # We try to find the actual nodes. If they don't exist (e.g. gaps), we skip linking (or create loose link?)
        # For now, MATCH existing
        session.run("""
            MATCH (f:Feature {id: $id})
            MATCH (s:RoadNode {contig: $contig, pos: $p1})
            MERGE (f)-[:STARTS_AT]->(s)
            WITH f
            MATCH (e:RoadNode {contig: $contig, pos: $p2})
            MERGE (f)-[:ENDS_AT]->(e)
        """, id=fid, gid=genome_id, contig=contig, p1=node_start_pos, p2=node_end_pos)

    def _increment_exon_path(self, session, contig, start, end, genome_id):
        # Range of K-mers fully within exon: [start, end - K]
        # Increment Node Freq
        limit = end - self.k
        if limit < start: return # Exon shorter than K
        
        session.run("""
            MATCH (n:RoadNode {contig: $contig})
            WHERE n.pos >= $start AND n.pos <= $limit
            SET n.freq = coalesce(n.freq, 0) + 1
        """, contig=contig, start=start, limit=limit)
        
        # Also increment edges within? 
        # Ideally yes. (n)-[r:NEXT]->(m).
        session.run("""
            MATCH (n:RoadNode {contig: $contig})-[r:NEXT]->(m:RoadNode)
            WHERE n.pos >= $start AND n.pos < $limit
            SET r.freq = coalesce(r.freq, 0) + 1
        """, contig=contig, start=start, limit=limit)

    def _flush_junctions(self, session, batch):
        session.run("""
            UNWIND $batch as item
            MATCH (n1:RoadNode {contig: item.c1, pos: item.p1})
            MATCH (n2:RoadNode {contig: item.c2, pos: item.p2})
            MERGE (n1)-[r:SPLICE_NEXT]->(n2)
            ON CREATE SET r.freq = 1, r.transcripts = [item.tx]
            ON MATCH SET r.freq = r.freq + 1, r.transcripts = r.transcripts + item.tx
        """, batch=batch)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--user", default="neo4j")
    parser.add_argument("--password", default="password")
    
    args = parser.parse_args()
    
    annotator = Neo4jAnnotator("bolt://localhost:7687", args.user, args.password)
    annotator.load_annotations(args.gff, args.genome)
    annotator.close()
