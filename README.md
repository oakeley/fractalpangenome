# Genomic Streetmap (Fractal Pangenome)

**Documentation & Concepts**: [https://github.com/oakeley/fractalpangenome](https://github.com/oakeley/fractalpangenome)

A privacy-preserving, spatially-indexed pangenome graph framework that maps genomic sequences to a "GPS" coordinate system (Hilbert Curve) and visualizes them as a 3D-Tile "Google Map" of the genome.

---

## 1. The Core Metaphor: City Streets and Journeys

The central concept of this architecture is the distinction between the **infrastructure** (the anonymous map) and the **usage** (the journeys).

### The Map (Infrastructure)
*   **K-mers (RoadNodes)**: The intersections or landmarks. Each unique K-mer (sequence content) exists only once in the graph at a specific coordinate.
*   **Transitions (Edges)**: The roads connecting them.
*   **Wormholes**: "Hyperspace" links connecting identical K-mers across distant genomic regions (repeats), creating a non-linear topology.
*   **3D Tiles**: A hierarchical QuadTree index (like Google Maps tiles) allowing us to visualize the genome at different zoom levels.

### The Journey (Usage)
*   **A Genome** (e.g., T2T-CHM13, Patient_001) is a specific **route** taken through this map.
*   **Highways**: Roads traveled by many genomes become "public highways" (High Frequency), representing conserved regions.
*   **Footpaths**: Roads traveled by few are "footpaths", representing rare variants or private mutations.
*   **Privacy**: Crucially, the map stores the *amount* of traffic (Frequency) but **NOT** the license plates (Genome IDs) of the travelers. This ensures the pangenome structure is persistent and useful without storing sensitive PII in the graph core.

---

## 2. Data Model Implementation

We implement this in Neo4j using a property graph model.

### Nodes
*   **`(:RoadNode {seq: "...", h_idx: 12345, freq: 100})`**: The fundamental unit.
    *   `seq`: The K-mer sequence (content).
    *   `h_idx`: Universal Hilbert Coordinate (Location).
    *   `freq`: How many times this node has been traversed (Privacy-preserving usage count).
*   **`(:SpatialTile {id: "z10_x1_y2"})`**: A container node for spatial query optimization.
*   **`(:Feature {id: "gene1", type: "gene"})`**: An annotation overlay (metadata) linked to the roads.

### Edges
*   **`[:NEXT {freq: 50}]`**: Represents physical adjacency on a chromosome.
    *   If 50 people have the sequence `...A -> T...`, the edge frequency is 50.
*   **`[:WORMHOLE {freq: 0}]`**: A zero-cost logical link between identical K-mers, allowing the graph to "fold" repeats.
*   **`[:SPLICE_NEXT {freq: 10}]`**: A jump across an intron (splicing), linked to specific transcripts.

---

## 3. Implications for Analysis

### Graph-Aware Alignment
Instead of aligning to a linear reference, we align to the Map.
1.  **Seeding**: Find read K-mers in the graph using Hilbert Coordinates.
2.  **Extension**: Traverse `[:NEXT]` edges to find the best path.
3.  **Variant Calling**: We project the graph's "World View" frequencies onto a linear reference (using `project_variants.py`) to see where an individual diverges from the "Highway".

### RNA-seq Alignment
Aligning RNA-seq becomes a pathfinding problem.
*   A read spanning an exon junction maps to `Exon_End` and `Exon_Start`.
*   The existence of a `[:SPLICE_NEXT]` edge confirms this is a known splice site.

---

## 4. Setup & Usage

### Prerequisites
*   **Linux OS** (Ubuntu 20.04+ recommended)
*   **Anaconda** / **Miniconda**
*   **Neo4j Database** via Docker

### Installation

1.  **Clone the Repository**:
    ```bash
    git clone https://github.com/oakeley/fractalpangenome.git
    cd fractalpangenome
    ```

2.  **Create Conda Environment**:
    ```bash
    conda create -n pangenome python=3.9 numpy pandas biopython neo4j pillow matplotlib networkx
    conda activate pangenome
    ```

3.  **Start Neo4j**:
    ```bash
    docker-compose -f docker-compose-neo4j.yml up -d
    ```
    *(Tuned for high-memory systems: 64GB Heap/PageCache)*

---

## 5. Building the Map

### Step A: The Pioneer (Foundation)
Streams a reference genome to build the road network and wormholes.
```bash
python3 fractal_google_map_builder.py --fasta prim.fasta --id PRIM --pioneer --clear
```

### Step B: The Explorers (Population)
Aligns additional genomes to increment "Highways" or forge new "Footpaths".
```bash
python3 fractal_google_map_builder.py --fasta sec.fasta --id SEC
```

### Step C: Annotation
Overlays GFF3 features.
```bash
python3 fractal_neo4j_annotations.py --gff annotations.gff3.gz --fasta prim.fasta --id PRIM
```

---

## 6. Variant Projection (Graph -> Linear VCF)

Map the graph's "World View" variants back to a specific linear assembly.
```bash
python3 project_variants.py --fasta my_reference.fasta --out results.vcf
```

---

## 7. DNA Alignment

Align reads directly to the anonymous graph (SAM output).
```bash
python3 graph_aligner_pangenome.py --fastq reads_R1.fastq.gz --out alignment.sam --cores 16
```

---

## 8. RNA-seq & Expression

1.  **Generate Transcript Map**:
    ```bash
    python3 generate_transcript_roads.py --gff annotations.gff3.gz --fasta prim.fasta --out transcript_roads.pkl
    ```
2.  **Align & Quantify**:
    ```bash
    python3 graph_aligner_pangenome.py --fastq rnaseq.fastq.gz --roads transcript_roads.pkl --out rnaseq.sam --rpkm expression.tsv
    ```

---

## 9. Visualization

**Global Expression Map**:
```bash
python3 visualize_hilbert.py --rpkm expression.tsv --roads transcript_roads.pkl --out global.png
```

**Zoomed Gene View**:
```bash
python3 visualize_hilbert.py --rpkm expression.tsv --roads transcript_roads.pkl --out gene_view.png --zoom XM_00001.1
```
