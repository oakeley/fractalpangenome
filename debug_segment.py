import h5py
import sys

def inspect_segment(db_path, segment_index=0):
    with h5py.File(db_path, 'r') as f:
        if 'path_segments' not in f:
            print("No path_segments group")
            return

        keys = list(f['path_segments'].keys())
        if not keys:
            print("No segments found")
            return
            
        seg_id = keys[segment_index]
        print(f"Inspecting Segment: {seg_id}")
        grp = f['path_segments'][seg_id]
        
        print("Attributes:")
        for k, v in grp.attrs.items():
            print(f"  {k}: {v}")
            
        print("Datasets:")
        for k in grp.keys():
            dset = grp[k]
            print(f"  {k}: shape={dset.shape}, dtype={dset.dtype}")
            if dset.size < 20:
                print(f"    Value: {dset[:]}")

if __name__ == "__main__":
    inspect_segment("pangenome_15ab.h5")
