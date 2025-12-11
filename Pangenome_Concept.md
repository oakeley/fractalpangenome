The Fractal Pangenome: Journeys on a K-mer Map
1. The Core Metaphor: City Streets and Journeys
The central concept of this pangenome architecture is the distinction between the infrastructure (the map) and the usage (the journeys).

The Map (Infrastructure):

K-mers (nodes) are the intersections or landmarks.
Transitions (edges) are the roads connecting them.
This map represents the total available sequence space of all observed genomes.
The Journey (Usage):

A Genome (e.g., "T2T-CHM13", "Individual_A_Paternal") is a specific route taken through this map.
A Gene or Transcript is a functional sub-journey, often with specific rules (like skipping sections via splicing).
Highways: Roads traveled by many genomes become "highways" (high edge counts), representing conserved regions.
Footpaths: Roads traveled by few genomes are "footpaths", representing rare variants or individual mutations.
2. Data Model Implementation
We implement this in Neo4j using a property graph model where edges carry the history of traversal.

Nodes
(:Kmer {seq: "..."}): The fundamental unit of the map. Unique across the entire pangenome.
(:Genome {id: "..."}): Represents an individual's genome.
(:Feature {id: "...", type: "gene"}): Represents a biological annotation.
Edges & The "Journey" Property
The key innovation is storing the list of travelers on the edges.

Genomic Path (:NEXT)
Represents physical adjacency on a chromosome.

(:Kmer)-[:NEXT {
    count: 150,                // The "Highway" score (150 genomes took this road)
    genomes: ["T2T", "HG002", ...] // The specific list of travelers
}]->(:Kmer)
New Genome Mapping: When adding a new genome, we mostly traverse existing roads, appending the new genome ID to the genomes list. We only build "new roads" (create new K-mer nodes/edges) when the genome contains a novel variant or sequence.
Transcriptomic Path (:SPLICE_NEXT)
Represents a jump across an intron (splicing).

(:Kmer)-[:SPLICE_NEXT {
    transcripts: ["ENST0001", "ENST0002"], // Isoforms taking this splice
    genomes: ["T2T"] // Genomes where this splice is valid/observed
}]->(:Kmer)
Alternative Splicing: A single K-mer node (exon end) can have multiple outgoing :SPLICE_NEXT edges pointing to different K-mer nodes (exon starts).
Path A: Exon1_End -> Exon2_Start (Isoform 1)
Path B: Exon1_End -> Exon3_Start (Isoform 2 - Exon skipping)
This naturally represents alternative splicing as branching journeys in the graph.
3. Implications for Alignment & Analysis
Graph-Aware Alignment
Instead of aligning to a linear reference, we align to the Map.

Seeding: Find read k-mers in the graph.
Extension: Traverse :NEXT edges.
Highway Scoring: Prefer edges with high 
count
 (conserved paths).
Journey Consistency: If a read maps to a specific rare variant (footpath), we can check the genomes list to see which specific population or individual has that path.
RNA-seq Alignment
Aligning RNA-seq becomes a pathfinding problem that allows traversing :SPLICE_NEXT edges.

A read spanning an exon junction will naturally map to the Exon_End k-mer and the Exon_Start k-mer.
The existence of a :SPLICE_NEXT edge confirms this is a known splice site.
Novel splice sites appear as reads bridging two k-mers that don't yet have a :SPLICE_NEXT edge (potential discovery).
Pangenome Growth
First Genome: Builds the initial city (all roads are new). Slow.
Subsequent Genomes: Mostly re-use existing roads. Fast.
Homologous Recombination: A child's genome might follow the "Mother's Journey" for a while, then switch at a junction to the "Father's Journey". The graph naturally captures this as a valid path composed of existing segments.
4. Summary
By treating biological sequences as journeys through a shared k-mer infrastructure, we create a dynamic, queryable pangenome that naturally handles:

Population variation (Highways vs. Footpaths)
Structural variants (New roads/Detours)
Alternative splicing (Branching transcript paths)
Efficient storage (Compression via shared structure)
