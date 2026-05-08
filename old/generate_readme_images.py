import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import random

def create_highway_concept_image(output_path):
    """
    Creates a schematic showing nodes/edges with varying thickness based on frequency.
    High frequency paths = Highways.
    Low frequency paths = Footpaths.
    """
    G = nx.DiGraph()
    
    # Create a main "highway" path
    highway_nodes = range(0, 6)
    for i in range(len(highway_nodes)-1):
        G.add_edge(highway_nodes[i], highway_nodes[i+1], weight=10, type='highway')
        
    # Create some "footpaths" (alternatives)
    # Detour 1
    G.add_edge(1, 6, weight=2, type='footpath')
    G.add_edge(6, 7, weight=2, type='footpath')
    G.add_edge(7, 3, weight=2, type='footpath')
    
    # Detour 2
    G.add_edge(2, 8, weight=1, type='footpath')
    G.add_edge(8, 4, weight=1, type='footpath')
    
    pos = {
        0: (0, 0), 1: (1, 0), 2: (2, 0), 3: (3, 0), 4: (4, 0), 5: (5, 0),
        6: (1.5, 1), 7: (2.5, 1),
        8: (3, -1)
    }
    
    plt.figure(figsize=(10, 6))
    
    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_color='lightblue', node_size=500)
    
    # Draw edges with varying thickness
    edges = G.edges(data=True)
    weights = [d['weight'] for u, v, d in edges]
    colors = ['darkblue' if d['type'] == 'highway' else 'gray' for u, v, d in edges]
    
    nx.draw_networkx_edges(G, pos, width=weights, edge_color=colors, arrowsize=20)
    nx.draw_networkx_labels(G, pos)
    
    plt.title("Genomic Highways vs Footpaths\nThicker edges = Higher Frequency (More Genomes)", fontsize=14)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Created {output_path}")

def create_genome_overlay_image(output_path):
    """
    Shows a new genome path overlaying existing paths.
    """
    G = nx.DiGraph()
    
    # Existing graph
    nodes = range(5)
    for i in range(4):
        G.add_edge(i, i+1, color='lightgray', style='solid', width=2)
    
    # Alternative existing path
    G.add_edge(1, 5, color='lightgray', style='solid', width=2)
    G.add_edge(5, 3, color='lightgray', style='solid', width=2)
    
    pos = {0: (0, 0), 1: (1, 0), 2: (2, 0), 3: (3, 0), 4: (4, 0), 5: (2, 1)}
    
    plt.figure(figsize=(10, 6))
    
    # Draw background graph
    edges = G.edges(data=True)
    colors = [d['color'] for u, v, d in edges]
    styles = [d['style'] for u, v, d in edges]
    widths = [d['width'] for u, v, d in edges]
    
    nx.draw(G, pos, edge_color=colors, width=widths, style=styles, with_labels=True, node_color='white', edgecolors='black')
    
    # Draw NEW genome path (overlay)
    # Path: 0 -> 1 -> 5 -> 3 -> 4
    new_path_edges = [(0, 1), (1, 5), (5, 3), (3, 4)]
    nx.draw_networkx_edges(G, pos, edgelist=new_path_edges, edge_color='red', width=4, alpha=0.7, label='New Genome Path')
    
    plt.legend(loc='upper right')
    plt.title("New Genome Overlaying Pangenome Graph\nRed path upgrades the frequency of traversed edges", fontsize=14)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Created {output_path}")

if __name__ == "__main__":
    create_highway_concept_image("highway_concept.png")
    create_genome_overlay_image("genome_overlay.png")
