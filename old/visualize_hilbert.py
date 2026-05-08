#!/usr/bin/env python3
"""
Visualize Hilbert Space Expression (Transcript Roads)
Generates PNG images of the Hilbert Curve colored by RPKM values.
Supports:
- Global View (Points)
- Zoom View (Transcript Roads/Trajectories)
"""

import sys
import os
import argparse
import logging
import csv
import math
import pickle
from collections import defaultdict
from PIL import Image, ImageDraw

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def hilbert_d2xy(n, d):
    """
    Convert Hilbert distance d to (x, y) coordinates for an n*n grid.
    n must be a power of 2.
    """
    x = 0
    y = 0
    s = 1
    t = d
    while s < n:
        rx = 1 & (t // 2)
        ry = 1 & (t ^ rx)
        
        if ry == 0:
            if rx == 1:
                x = s - 1 - x
                y = s - 1 - y
            x, y = y, x
            
        x += s * rx
        y += s * ry
        s *= 2
        t //= 4
    return x, y

class HilbertVisualizer:
    def __init__(self, rpkm_file, roads_file, resolution=1024):
        self.rpkm_data = {} # tx_id -> rpkm
        self.transcript_paths = {} # tx_id -> [h_indices]
        self.resolution = resolution
        
        self.load_rpkm(rpkm_file)
        self.load_roads(roads_file)
        
    def load_rpkm(self, rpkm_file):
        logger.info(f"Loading RPKM from {rpkm_file}...")
        try:
            with open(rpkm_file, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    tx_id = row['TranscriptID']
                    rpkm = float(row['RPKM'])
                    self.rpkm_data[tx_id] = rpkm
        except Exception as e:
            logger.error(f"Failed to load RPKM: {e}")

    def load_roads(self, roads_file):
        logger.info(f"Loading Transcript Roads from {roads_file}...")
        try:
            with open(roads_file, 'rb') as f:
                data = pickle.load(f)
                self.transcript_paths = data.get('transcript_paths', {})
            logger.info(f"Loaded paths for {len(self.transcript_paths)} transcripts.")
        except Exception as e:
            logger.error(f"Failed to load roads: {e}")

    def generate_global_plot(self, output_path):
        logger.info(f"Generating global plot ({self.resolution}x{self.resolution})...")
        img = Image.new('RGB', (self.resolution, self.resolution), "black")
        pixels = img.load()
        
        N = 2**31
        
        # Iterate all transcripts with RPKM > 0
        for tx_id, path in self.transcript_paths.items():
            rpkm = self.rpkm_data.get(tx_id, 0)
            if rpkm == 0: continue
            
            color = self._rpkm_to_color(rpkm)
            
            # Sample points from path to plot
            # Path can be long. Plot every Nth point?
            step = max(1, len(path) // 100)
            
            for i in range(0, len(path), step):
                h = path[i]
                X, Y = hilbert_d2xy(N, h)
                
                img_x = int(X / N * self.resolution)
                img_y = int(Y / N * self.resolution)
                
                if 0 <= img_x < self.resolution and 0 <= img_y < self.resolution:
                    pixels[img_x, img_y] = color
                            
        img.save(output_path)
        logger.info(f"Saved to {output_path}")

    def generate_zoom_plot(self, target_tx_id, output_path):
        if target_tx_id not in self.transcript_paths:
            logger.error(f"Transcript {target_tx_id} not found in roads.")
            return

        path = self.transcript_paths[target_tx_id]
        logger.info(f"Zooming into {target_tx_id} (Path len: {len(path)})...")
        
        # Determine bounding box of the path
        N = 2**31
        xs = []
        ys = []
        
        # Convert all points in path to XY
        # (If path is huge, maybe sample? But for zoom we want detail)
        points = []
        for h in path:
            x, y = hilbert_d2xy(N, h)
            xs.append(x)
            ys.append(y)
            points.append((x, y))
            
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        
        # Add padding
        width = max_x - min_x
        height = max_y - min_y
        padding = max(width, height) * 0.1
        
        min_x -= padding
        max_x += padding
        min_y -= padding
        max_y += padding
        
        width = max_x - min_x
        height = max_y - min_y
        scale = max(width, height)
        if scale == 0: scale = 1
        
        img = Image.new('RGB', (1024, 1024), "black")
        draw = ImageDraw.Draw(img)
        
        # Draw the main transcript path
        rpkm = self.rpkm_data.get(target_tx_id, 0)
        color = self._rpkm_to_color(rpkm)
        if rpkm == 0: color = (100, 100, 255) # Default blue if no expression
        
        prev_img_x, prev_img_y = None, None
        
        for x, y in points:
            img_x = int((x - min_x) / scale * 1000) + 12
            img_y = int((y - min_y) / scale * 1000) + 12
            
            if prev_img_x is not None:
                draw.line([(prev_img_x, prev_img_y), (img_x, img_y)], fill=color, width=3)
            
            prev_img_x, prev_img_y = img_x, img_y
            
        # Draw other transcripts in the background?
        # Maybe too expensive to check all.
        
        img.save(output_path)
        logger.info(f"Saved zoom to {output_path}")

    def _rpkm_to_color(self, rpkm):
        # Heatmap: Blue -> Green -> Red
        if rpkm <= 0: return (50, 50, 50)
        val = math.log1p(rpkm)
        norm = min(val / 10.0, 1.0)
        
        r = int(255 * norm)
        g = int(255 * (1 - abs(0.5 - norm) * 2))
        b = int(255 * (1 - norm))
        return (r, g, b)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Hilbert Visualizer (Roads)")
    parser.add_argument("--rpkm", required=True, help="RPKM TSV file")
    parser.add_argument("--roads", required=True, help="Transcript Roads Pickle")
    parser.add_argument("--out", required=True, help="Output PNG file (Global)")
    parser.add_argument("--zoom", help="Transcript ID to zoom into")
    
    args = parser.parse_args()
    
    viz = HilbertVisualizer(args.rpkm, args.roads)
    viz.generate_global_plot(args.out)
    
    if args.zoom:
        zoom_out = args.out.replace('.png', f'_{args.zoom}.png')
        viz.generate_zoom_plot(args.zoom, zoom_out)
