#!/usr/bin/env python3
"""
Verify Graph Connectivity
Checks if the pangenome graph in HDF5 is correctly connected.
"""

import h5py
import sys
import logging
import os
from collections import defaultdict

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def verify_graph(db_path):
    logger.info(f"Verifying graph connectivity for {db_path}...")
    
    try:
        with h5py.File(db_path, 'r') as f:
            if 'path_segments' not in f:
                logger.error("No 'path_segments' group found in HDF5 file.")
                return
            
            segments_group = f['path_segments']
            total_segments = len(segments_group.keys())
            logger.info(f"Total segments found: {total_segments:,}")
            
            # 1. Collect all valid segment IDs
            valid_segment_ids = set()
            for seg_id_str in segments_group.keys():
                valid_segment_ids.add(int(seg_id_str))
            
            # 2. Check transitions
            orphan_segments = 0
            dead_end_segments = 0
            invalid_transitions = 0
            
            # Track incoming connections to find orphans
            incoming_counts = defaultdict(int)
            
            for i, seg_id_str in enumerate(segments_group.keys()):
                if i % 10000 == 0:
                    logger.info(f"Processed {i:,} segments...")
                
                segment_id = int(seg_id_str)
                segment_data = segments_group[seg_id_str]
                
                # Check outgoing transitions
                if 'transitions' in segment_data:
                    transitions = segment_data['transitions'][:]
                    if len(transitions) == 0:
                        dead_end_segments += 1
                    
                    for target_id in transitions:
                        incoming_counts[target_id] += 1
                        if target_id not in valid_segment_ids:
                            logger.warning(f"Invalid transition: Segment {segment_id} -> {target_id} (Target does not exist)")
                            invalid_transitions += 1
                else:
                    dead_end_segments += 1
            
            # Check for orphans (no incoming transitions)
            # Note: The very first segment(s) might legitimately be orphans if they are start nodes
            for seg_id in valid_segment_ids:
                if incoming_counts[seg_id] == 0:
                    orphan_segments += 1
            
            logger.info("-" * 40)
            logger.info("Graph Verification Results:")
            logger.info(f"Total Segments: {total_segments:,}")
            logger.info(f"Dead Ends (No outgoing): {dead_end_segments:,} ({dead_end_segments/total_segments*100:.2f}%)")
            logger.info(f"Orphans (No incoming): {orphan_segments:,} ({orphan_segments/total_segments*100:.2f}%)")
            logger.info(f"Invalid Transitions: {invalid_transitions:,}")
            
            if invalid_transitions == 0:
                logger.info("SUCCESS: All transitions point to valid segments.")
            else:
                logger.error("FAILURE: Found invalid transitions.")

    except Exception as e:
        logger.error(f"Error verifying graph: {e}")

if __name__ == "__main__":
    db_path = "pangenome_15ab.h5"
    if len(sys.argv) > 1:
        db_path = sys.argv[1]
    
    if not os.path.exists(db_path):
        logger.error(f"Database file not found: {db_path}")
        sys.exit(1)
        
    verify_graph(db_path)
