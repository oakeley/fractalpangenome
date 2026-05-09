# Pangenome Hilbert-Curve Spatial Indexing for Genomics

## What This Is

This tool maps DNA sequences into a three-dimensional space and builds a navigable
street map of the genomic world one that accumulates knowledge from every genome
you add, without ever erasing what was learned before.

The street map metaphor is precise. Just as OpenStreetMap records roads and their
quality without storing the names or personal data of the individuals who drove
them, this database records Hilbert-space paths and their usage without storing
the individual genome journeys that created them. Add ten reference genomes and
the map gets richer. Add a hundredth and roads shared by many genomes become
motorways while rare variants remain dirt tracks all without the database
growing proportionally, because only one copy of each street is kept regardless
of how many genomes passed through it.

---

## The Core Idea: Hilbert Cubes

### From Sequence to GPS Address

Every k-mer (a window of k consecutive DNA bases) has a unique, deterministic
GPS address in a three-dimensional cube. The mapping is:

1. The k-mer is hashed to a canonical 30-bit integer using 2-bit encoding
   (A=0, C=1, G=2, T=3). The canonical form takes the minimum of the forward
   and reverse-complement hashes, so a sequence and its complement map to the
   same address strand is irrelevant.
2. For k=15, the canonical hash is exactly 30 bits (2×15), which fills the
   30-bit Hilbert index space completely. The hash *is* the Hilbert index no further transformation is needed.
3. The 30-bit Hilbert index is decoded into normalised (x, y, z) coordinates
   in [0,1]³ using the Skilling (2004) 3D Hilbert curve algorithm. This is
   the GPS address.

**Critical property locality:** The mapping is local. Each k-mer's GPS address
depends only on its own sequence, not on any surrounding context. This means:

- One SNP at position p affects at most k trajectory points (the k sliding
  windows that overlap p). The other thousands of k-mers in a read map to
  exactly the same addresses as the reference.
- A read can start anywhere along a chromosome or transcript and still align
  correctly it does not need to start at the reference contig's beginning.
- Reads from a different individual of the same species align robustly because
  ~99% of k-mers are identical.

**Why not a cumulative walk?** Accumulated k-mer hashes into
a running position (`pos += hash.rotate_left(i)`) cumulative walk is tempting but a single
SNP at position i will shift every subsequent trajectory point, making it impossible
to align a read from an individual with per-individual variation against the current observed streetmap. 
A local hash design allows for a least worst outcome when aligning reads from an individual whos genome has not been observed previously.

### The Google Earth Analogy

Think of the Hilbert cube as a globe you can fly over in Google Earth:

| Google Earth concept | Hilbert pangenome concept |
|---|---|
| The globe | The 3D Hilbert cube (30-bit sequence space, 1024³ cells) |
| A city | A cluster of similar k-mer sequences (e.g. a gene family) |
| A motorway | A Hilbert path traversed by hundreds of genomes |
| A dirt track | A rare variant path seen in only one or two samples |
| A side street | A splice junction valid path not on the main chromosomal road |
| Zoom to fit | `zoom` command finds the tightest view around two Hilbert positions |
| GPS coordinates | Normalised (x, y, z) coordinates in the Hilbert cube |
| Navigate A>B | `route` command finds the optimal path using road quality as weight |
| Saved bookmark | `waypoint` command names a Hilbert position (gene locus, hotspot) |
| Satellite imagery | Read alignment chains k-mer seeds against known mapped roads |

---

## The Street Map: How the Database Is Built

### OpenStreetMap, Not a GPS Track Log

OpenStreetMap stores streets, not journeys. When you drive to work, your
route validates and improves the road quality of every street you used, but your
personal journey is not stored. The pangenome database works identically:

- **Tile files** (`tile_*.bin`) store Hilbert cells with usage counts and spatial
  coordinates. They contain no genome names, no sample IDs, no journey records.
- **`registry.json`** stores per-genome metadata (name, contig lengths, Hilbert
  ranges). This is a bookshelf log it is separate from the map itself.
- **`kmer_index.bin`** stores the positional k-mer index used for alignment
  (analagous to navigation by street names rather than GPS coordinates). This is a translation layer and separate from the Hilbert cell data.

When a new genome is added, its k-mer trajectory increments usage counts of
every Hilbert cell it visits. Cells visited by many genomes accumulate high usage
and are labelled high-quality roads. Rare variant cells remain low-usage.

### Road Quality Tiers

| Usage count | Road quality | Meaning |
|---|---|---|
| 0 | Desert | Cell exists but has never been traversed |
| 1-9 | Dirt Track | Seen in fewer than 10 genomes |
| 10-49 | Path | Rare but reproducible variant |
| 50-199 | Side Street | Present in a minority of samples |
| 200-999 | Road | Common variant |
| 1,000-4,999 | Main Road | Highly conserved sequence region |
| 5,000-49,999 | Motorway | Near-universal across genomes |
| 50,000+ | Pan-Continent Express | Universally conserved core genome |

These tiers emerge naturally from usage counts without any manual annotation.
Housekeeping genes and ribosomal RNA automatically become motorways. Private
mutations become dirt tracks. Cross-species conserved sequences (e.g. GAPDH
shared between human and mouse) accumulate usage from both species simultaneously.

### Chromosomes Are Processed Whole

Each chromosome is read in full and its Hilbert trajectory is computed in a
single unbroken pass. No artificial chunk boundaries are introduced. This ensures
that k-mers at positions i and i+1 are always consecutive entries in the
trajectory no gaps, no false discontinuities in the road.

### Transcripts Are First-Class Roads

When a GFF3 annotation file is supplied alongside the FASTA during `add-genome`,
each mRNA/transcript is processed as a separate contig. Its exon sequences are
extracted from the FASTA and concatenated in 5′>3′ transcript order (with minus-
strand exons reverse-complemented). The resulting spliced sequence is processed
exactly as a chromosomal contig: its k-mer trajectory is computed, voxels are
inserted into the tile grid, and minimizers are added to the kmer index.

**Splice junctions as side streets:** The sliding window at the exon-exon boundary
produces k-mers that span both exons. These k-mers have GPS addresses near both
exons' territories in the Hilbert cube. They form a "side street" connecting the
two exonic roads a valid path that RNA polymerase actually takes but that
genomic DNA does not. A spliced RNA-seq read aligns to the transcript contig's
chain; an unspliced genomic read aligns to the chromosomal contig instead.

### Structural Variation Creates New Roads Automatically

When a second individual's genome is added (`add-genome` for a second sample),
k-mers that differ from the existing road structure (due to SNPs, indels, or structural
variants) produce different canonical hashes and therefore different Hilbert
GPS addresses. Those addresses map to Hilbert cells in the same way as a city may grow over
time, new roads can be added, old ones upgraded but the GPS coordinates remain universal and
unchanging. So the tile grid inserts them as new voxels: new streets on the map.
The kmer index records these new minimizers under the new contig name.

When a read spanning the structural variant is subsequently aligned:
- Its k-mers matching the variant sequence chain strongly against the new contig
- Its k-mers matching the reference sequence chain against the reference contig
- Both are reported (up to `--max-hits`); the variant contig wins if it provides
  a longer chain for that read

No special SV-calling step is required the map accumulates variation naturally.

---

## The Tile System: Disk-Based with an LRU Cache

### Why Tiles?

A human genome at k=15 produces ~3 billion k-mers. Storing all corresponding
Hilbert cells in RAM simultaneously requires tens of gigabytes. The tile system
partitions the 3D Hilbert cube into spatial tiles (each covering a 64³ cube of
Hilbert cells, for a total of 16³ = 4096 tiles) and keeps only a working set
in RAM via an LRU cache.

### Write-Back Cache With Load-Before-Insert

The cache is write-back with a critical invariant: before inserting a k-mer
into a tile, that tile is loaded from disk if it is not already in the cache.
This ensures prior usage counts from previous genome ingestions are present in
memory before the new traversal increments them. The sequence for every k-mer is:

1. Compute the tile ID from the k-mer's Hilbert index.
2. If the tile is not in the LRU cache: read it from disk into the cache.
3. Increment the voxel's usage count in the in-memory tile.
4. After every complete chromosome: flush all dirty tiles to disk and evict them.

Flushing writes the in-memory tile directly to disk using an atomic temp-file
rename (`tile_*.tmp` > `tile_*.bin`). Because prior usage counts were loaded in
step 2, the flushed tile already contains the fully merged state no second
read-and-merge at flush time is needed (which would double-count usage).

### /dev/shm Staging for GFF3 Transcript Ingestion

Processing 181,000+ transcripts with per-transcript tile flushes against a NAS
would require ~181,000 × (read tile from disk + write tile to disk) round trips.
To avoid this, the GFF3 annotation pass stages the entire tile database into
`/dev/shm` (Linux shared memory) at the start of the pass, performs all
transcript insertions at RAM speed, then copies the updated tiles back to the
real database path. Word to the wise... don't build this stuff on your laptop!

### Atomic Tile Writes and Crash Safety

```
write tile_0_3_7_2_.tmp
rename tile_0_3_7_2_.tmp > tile_0_3_7_2_.bin < POSIX atomic
```

On POSIX filesystems `rename(2)` is atomic with respect to crashes: either the
old file survives or the new one appears, never a half-written state.

---

## The Kmer Index: Positional Index for Alignment

### Why a Separate Index?

The tile files record which Hilbert cells are occupied (the road network) but
not which contig each voxel came from or at what position within that contig.
Correct alignment requires knowing not just "this k-mer is on a road" but "this
k-mer appears at position 4721 of transcript ENST00000456328 in human chromosome
1". Without positional information it is impossible to distinguish between a
read that genuinely aligns to a contig (its k-mers appearing in the correct
order) and a read that merely shares some k-mers by coincidence.

### Structure of kmer_index.bin

The kmer index is a sorted flat binary file of `(hash:u32, contig_idx:u32,
pos:u32)` triplets, where:

- `hash` is the canonical k-mer hash (30-bit value, stored as u32)
- `contig_idx` is an index into the contig name table stored in the file header
- `pos` is the position of this k-mer in the contig's trajectory (0-based k-mer
  index, not base position)

The file is sorted by `hash` so that lookup by any k-mer hash is a binary search
(O(log n) per lookup, no HashMap overhead). For a full human genome with ~181,000
transcripts, the file is approximately 3-4 GB on disk and ~10 GB in RAM when loaded.

### Minimizers Reduce Index Size by 10×

Storing every k-mer position would require ~3 billion entries for a human genome
(~36 GB). Instead, the index stores only minimizers: for each sliding window
of `WINDOW=10` consecutive k-mers, only the position of the k-mer with the
minimum hash value is stored. Consecutive windows that share the same minimizer
position are deduplicated.

This reduces the index to ~300 million entries (~3.6 GB) while preserving enough
seeds for accurate chaining. A read from a real gene will find minimizer matches
because its k-mer sequences are identical to those used when the index was built
the same k-mer always has the same hash value.

Degenerate k-mers (canonical hash < 256) are excluded from both index construction
and alignment. This threshold filters out homopolymer runs (poly-A, poly-T, etc.)
which all hash to values near zero and would otherwise create a spurious
high-density cluster at h≈0 that attracts every read containing a poly-A tail.

### Incremental Index Updates

When a second genome is added with `add-genome`, its new minimizers are appended
to the existing `kmer_index.bin`. The existing entries are preserved; the new
contig names are assigned new indices in the header contig table. The merged file
is re-sorted and deduplicated before writing.

---

## The Alignment Algorithm: Seed-Chain (minimap2-Style)

### Overview

The alignment algorithm follows the same three-phase approach used by minimap2:
**seed > chain > score**. No reference FASTA is required at alignment time
only the `kmer_index.bin` and `registry.json`.

### Phase 1 Seed Collection

For each read, the same minimizer selection procedure used during index
construction is applied to the read's trajectory:

1. Compute the read's k-mer trajectory: for each sliding window of k bases,
   compute the canonical hash (local, not cumulative).
2. For each window of `WINDOW=10` consecutive trajectory points, find the
   position of the minimum non-degenerate hash (the read minimizer).
3. For each read minimizer, binary-search `kmer_index.bin` to retrieve all
   `(contig_idx, ref_pos)` pairs where this exact hash appears in the reference.
4. Group the resulting seeds by `contig_idx`.

A seed `(ref_pos, read_pos)` represents a k-mer that appears at position
`read_pos` in the read and at position `ref_pos` in the reference contig.

### Phase 2 Colinear Chaining

For each contig that has at least 2 seeds:

1. Sort seeds by `ref_pos` ascending.
2. Find the longest colinear chain: the largest subset of seeds where both
   `ref_pos` and `read_pos` are strictly increasing (O(n log n) patience-sort
   LIS algorithm on the `read_pos` values).
3. Record the chain length.

**Why ordering matters:** A read from gene A has minimizer k-mers that appear in
gene A's reference in the same 5′>3′ order as they appear in the read producing
a long colinear chain. A read from gene B has k-mers scattered at unrelated
`ref_pos` values in gene A's reference no long chain forms. Bacterial DNA
against a human database: almost no seeds, no chains, no alignment.

**SNP tolerance:** One SNP disrupts at most k=15 consecutive minimizer positions
in the read. For a 1000 bp read at WINDOW=10, that is at most 2 disrupted
minimizers out of ~98 total. The chain loses 2 points but still spans the full
read length in reference coordinates.

**Partial alignment:** Because each k-mer seeds independently, a read starting
mid-chromosome chains correctly against the reference starting at the corresponding
mid-chromosome position. No requirement for the read to begin at the contig start.

### Phase 3 Scoring and Contig Assignment

```
score = best_chain_length / non_degenerate_minimizer_count_in_read
```

- Score = 1.0 means every non-degenerate read minimizer is in the best chain.
- Score = 0.5 means half are in the chain (consistent with ~50% sequence identity).
- Score < `--min-score` (default 0.05): read is reported as unmapped.

The best contig is the one with the longest chain. Transcript-level contigs are
preferred over chromosomal contigs on equal chain length (shorter reference span
= more specific functional assignment).

MAPQ is computed as `-10 × log10(1 - score)`, capped at 60, mirroring standard
Phred-scaled mapping quality conventions.

### Soft-Clipping of Adapter and Homopolymer Runs

PacBio SMRTbell adapters and poly-A tails at read ends produce k-mers with
canonical hash < DEGEN (256). These are identified by finding the first and last
non-degenerate k-mer in the trajectory and marking the leading/trailing bases as
soft-clipped in the CIGAR string (e.g. `22S865M12S`). The adapter bases do not
contribute to the alignment score.

### Genome Name Resolution: O(1) Lookup

The contig name returned by chain scoring identifies which registered contig the
read aligns to. To resolve the genome name (e.g. `T2T-CHM13` vs `GRCm39`), a
`HashMap<contig_name, genome_name>` is built once at `AlignmentEngine` construction
from the registry. This makes the lookup O(1) regardless of the number of genomes
(important for xenograft samples where both human and mouse genomes are present).

### Multi-Species Alignment (e.g. Mouse-Human Xenografts)

The Hilbert k-mer address space is shared across all species simultaneously.
When a pangenome database contains both human and mouse genomes:

- A human read chains strongly against human transcripts. It also finds seeds
  in mouse orthologs but chains less well (sequence divergence reduces chain
  length proportionally to divergence).
- A mouse read chains strongly against mouse transcripts.
- A read from a conserved exon (e.g. ~90% identity between human and mouse
  GAPDH) produces two competing chains. Both are reported up to `--max-hits`.
  The `ZG:Z:` SAM tag identifies the species; the contig name identifies the
  transcript within that species.
- Cross-species conserved k-mers automatically accumulate high road quality
  (usage counts from both human and mouse traversals), becoming motorways.
  Species-specific k-mers remain lower-quality roads.

---

## Building

Requires Rust 1.75+ and Cargo.

```bash
# Build the release binary (strongly recommended for real genomes)
cargo build --release

# The binary is at:
./code/target/release/pangenome
```

---

## Command Reference

All subcommands share global options passed before the subcommand name.

```
pangenome [--db <path>] [--lru-cache <n>] <subcommand> [subcommand options]
```

### Global Options

| Option | Default | Description |
|---|---|---|
| `--db` / `-d` | `./pangenome_db` | Path to the database directory (created on first use) |
| `--lru-cache` | 4096 | Number of tiles to keep in RAM. Each tile ≈ 5-6 MB on disk; 4096 tiles ≈ 128 MB cache |

The database directory contains:
- `registry.json` genome metadata (names, contig names, lengths, Hilbert min/max)
- `tile_L_X_Y_Z_.bin` spatial tile files (the road map, 4096 files for a full human genome)
- `kmer_index.bin` sorted positional minimizer index for alignment (~3-4 GB)

---

### `add-genome` Add a Genome to the Pangenome

Loads a FASTA file, computes the local k-mer Hilbert trajectory for each
contig, inserts voxels into the tile grid (incrementing road usage counts),
builds minimizers for the positional kmer index, and saves everything to disk.
If `--gff3` is provided, transcript trajectories are also computed and indexed.

```bash
pangenome --db <db> add-genome \
    --name <genome_name> \
    --fasta <file.fa[.gz]> \
    [--kmer <k>] \
    [--gff3 <annotations.gff3>]
```

**Options:**

| Option | Default | Description |
|---|---|---|
| `--name` / `-n` | required | Logical name, e.g. `T2T-CHM13`, `GRCm39`, `sample_042` |
| `--fasta` / `-f` | | Path to FASTA file (plain or `.gz`) |
| `--seq` / `-s` | | Inline sequence string (for tests or short sequences) |
| `--kmer` / `-k` | 15 | k-mer length. Must match the k used to create the database. For k=15 the canonical hash is exactly 30 bits, filling the Hilbert space exactly |
| `--gff3` | | GFF3 annotation file. Each mRNA/transcript becomes a spliced trajectory road; splice junctions form side streets connecting exons |

**What happens internally:**

1. Each chromosome is streamed from the FASTA one at a time (only one chromosome
   in RAM at once). `seq_to_trajectory` computes the local k-mer hash trajectory.
2. For each k-mer: its tile is loaded from disk (if not already in the LRU
   cache), and the voxel's usage count is incremented.
3. After each chromosome: all dirty tiles are flushed to disk and evicted from
   the LRU cache. RAM usage drops back to approximately the LRU working set.
4. Minimizers from the chromosome are recorded in a `KmerIndexBuilder`.
5. After all chromosomes: `kmer_index.bin` is written (merging with any
   existing index from previous `add-genome` calls).
6. `registry.json` is updated with the new genome's metadata.

If `--gff3` is provided, steps 1-5 are repeated for each transcript contig
after the chromosomal pass, using `/dev/shm` staging for speed.

**Examples:**

```bash
# Add the T2T-CHM13 human reference with transcript annotations
pangenome --db ./human.k15 \
    add-genome --name T2T-CHM13 \
    --fasta T2T-CHM13.fasta \
    --gff3 T2T-CHM13.gff3 \
    --kmer 15

# Add a second human assembly shared roads increment, new variants add streets
pangenome --db ./human.k15 \
    add-genome --name HG002 \
    --fasta HG002_haplotype1.fasta

# Add a mouse genome to an existing human DB builds a mixed-species pangenome
pangenome --db ./human_mouse.k15 \
    add-genome --name GRCm39 \
    --fasta GRCm39.fasta \
    --gff3 GRCm39.gff3

# Quick test with an inline sequence
pangenome add-genome --name test --seq ACGTACGTACGTACGTACGT --kmer 8
```

**Important:** All genomes in the same database must use the same k-mer length.
The k value is fixed at database creation and stored in `registry.json`. The
`kmer_index.bin` file must be rebuilt whenever `add-genome` is run this
happens automatically. If you have an existing database built with an older
version of this tool (before the local-hash trajectory was introduced), you
must rebuild it from scratch because the Hilbert addresses stored in the tiles
and index were computed with the old cumulative-walk algorithm and are
incompatible with the current local-hash algorithm.

---

### `export-bam` Align Reads and Export SAM

Loads the positional kmer index once, then aligns all reads in parallel using
the seed-chain algorithm. Writes a SAM file incrementally (each read's result
is written as it completes the file grows in real time). Optionally writes
per-transcript and per-exon read count tables.

```bash
pangenome --db <db> export-bam \
    --reads <reads.fastq.gz> \
    --output <out.sam> \
    [--threads <n>] \
    [--gff3 <annotations.gff3>] \
    [--show-bam-cmd]
```

**Options:**

| Option | Default | Description |
|---|---|---|
| `--reads` / `-r` | | FASTA or FASTQ file of reads, plain or `.gz` |
| `--read-seq` | | Single inline read sequence |
| `--read-name` | `query_read` | Name label for an inline read |
| `--output` / `-o` | `alignments.sam` | Output SAM file path |
| `--threads` | 4 | Number of parallel alignment threads. The kmer index is loaded once and Arc-shared; each thread uses its own `AlignmentEngine` with O(1) lookups |
| `--gff3` | | If provided: reads aligning to `::transcript::` contigs are counted; `<output>.transcripts.tsv` and `<output>.exons.tsv` are written alongside the SAM |
| `--show-bam-cmd` | false | Print the samtools command to convert SAM > BAM |

**SAM custom tags:**

| Tag | Type | Description |
|---|---|---|
| `HI:i` | int | Lower 32 bits of the best-chain Hilbert index |
| `HL:i` | int | Upper 32 bits of the Hilbert index |
| `HX:f` | float | GPS x coordinate (Hilbert cube x in [0,1)) |
| `HY:f` | float | GPS y coordinate |
| `HZ:f` | float | GPS z coordinate |
| `ZG:Z` | string | Genome name of the best-matching contig |
| `AS:f` | float | Chain score (chain_length / non_degen_minimizer_count) |
| `NM:i` | int | Non-degenerate minimizers not in the best chain (approximate mismatches) |

**Transcript and exon count files:**

`<output>.transcripts.tsv` one row per transcript contig that received at
least one aligned read:
```
transcript_id    read_count    reads_per_thousand
ENST00000456328  14823         782.1
ENST00000442987   3291         173.7
```

`<output>.exons.tsv` requires `--gff3`. Read counts distributed evenly
across exons of each transcript (proportional distribution; true per-exon
sub-alignment is not performed):
```
transcript_id      exon_num  seqname  start    end      strand  read_count  reads_per_thousand
ENST00000456328    1         chr1     11869    12227    +       2471        130.4
ENST00000456328    2         chr1     12613    12721    +       2471        130.4
```

**Examples:**

```bash
# Align PacBio subreads, 24 threads, write SAM
pangenome --db ./human.k15 export-bam \
    --reads SRR36734339_subreads.fastq.gz \
    --output alignments.sam \
    --threads 24

# Align with transcript counts using the same GFF3 used at add-genome time
pangenome --db ./human.k15 export-bam \
    --reads SRR36734339_subreads.fastq.gz \
    --output alignments.sam \
    --threads 24 \
    --gff3 T2T-CHM13.gff3

# Convert to BAM
samtools view -bS alignments.sam | samtools sort -o alignments.bam
samtools index alignments.bam

# Align against a mouse-human pangenome (xenograft sample)
# Human reads chain to human contigs; mouse reads chain to mouse contigs
pangenome --db ./human_mouse.k15 export-bam \
    --reads xenograft_sample.fastq.gz \
    --output xenograft.sam \
    --threads 24 \
    --gff3 combined_human_mouse.gff3
```

---

### `align` Align Reads (Text Output)

Same alignment engine as `export-bam` but prints results to stdout in a human-
readable format rather than writing SAM.

```bash
pangenome --db <db> align \
    [--reads <file>] [--read-seq <seq>] [--read-name <name>] \
    [--max-hits <n>] [--min-score <f>]
```

**Options:**

| Option | Default | Description |
|---|---|---|
| `--reads` / `-r` | | FASTA or FASTQ file (plain or `.gz`) |
| `--read-seq` | | Single inline read |
| `--read-name` | `query_read` | Name for inline read |
| `--max-hits` / `-m` | 5 | Maximum hits to report per read |
| `--min-score` | 0.05 | Minimum chain score threshold |

```bash
# Align a single sequence and print GPS, contig, and chain score
pangenome --db ./human.k15 align \
    --read-seq ACGTACGTACGTACGTACGTACGTACGT \
    --read-name probe_1
```

---

### `annotate` Query Genomic Annotations

Loads GFF3 or BED files, maps each feature to a Hilbert position using a
deterministic coordinate hash, and queries features near a given position.

**Note:** `annotate` uses a generic hash for GPS assignment that is independent
of any specific genome. For genome-specific GPS coordinates, use the transcript
roads built by `add-genome --gff3` those are correctly positioned in Hilbert
space according to the actual k-mer sequence of each transcript.

```bash
pangenome --db <db> annotate \
    [--gff3 <file>] [--bed <file>] \
    [--seqname <chr>] [--start <n>] [--end <n>] \
    [--hilbert <idx>] [--radius <n>] \
    [--summary]
```

**Options:**

| Option | Default | Description |
|---|---|---|
| `--gff3` | | GFF3 annotation file |
| `--bed` | | BED annotation file |
| `--seqname` | | Chromosome name for linear coordinate query |
| `--start` | | Linear start position (0-based) |
| `--end` | | Linear end position (0-based, exclusive) |
| `--hilbert` | | Query near a specific Hilbert index (overrides `--seqname/--start/--end`) |
| `--radius` | 1,000,000 | Search radius in Hilbert index units for `--hilbert` queries |
| `--summary` | false | Show feature-type counts only |

```bash
# Summarise feature types in a GFF3
pangenome annotate --gff3 T2T-CHM13.gff3 --summary

# Query features at a linear coordinate range
pangenome annotate --gff3 T2T-CHM13.gff3 \
    --seqname chr1 --start 1000000 --end 2000000

# Query features near a Hilbert position from an alignment result
pangenome annotate --gff3 T2T-CHM13.gff3 \
    --hilbert 823456789 --radius 500000
```

---

### `route` Plan a Route Through the Map

Finds the optimal path between two Hilbert positions using A* search weighted
by road quality. High-usage (high-quality) roads are preferred, exactly as a
sat-nav chooses motorways over dirt tracks.

```bash
pangenome --db <db> route --from <hilbert_idx> --to <hilbert_idx> [--max-steps <n>]
```

| Option | Default | Description |
|---|---|---|
| `--from` | required | Start Hilbert index |
| `--to` | required | End Hilbert index |
| `--max-steps` | 20 | Maximum steps to display (full path is computed regardless) |

```bash
pangenome --db ./human.k15 route --from 100000000 --to 200000000
```

---

### `waypoint` Named Markers

Add and query named waypoints bookmarks for gene loci, mutation hotspots, or
regulatory elements.

```bash
pangenome --db <db> waypoint add  --name <n> --hilbert <idx> [--tag key=value ...]
pangenome --db <db> waypoint list [--tag key=value]
pangenome --db <db> waypoint near --hilbert <idx> [--radius <n>]
```

```bash
pangenome --db ./human.k15 waypoint add \
    --name BRCA1 --hilbert 823456789 \
    --tag type=gene --tag disease=breast_cancer

pangenome --db ./human.k15 waypoint list --tag type=gene
pangenome --db ./human.k15 waypoint near --hilbert 823456789 --radius 1000000
```

Note: waypoints are in-memory for the session duration and are not currently
persisted between invocations.

---

### `zoom` Zoom to Fit (Google Earth Style)

Finds the tightest Hilbert tile containing both positions, then optionally
zooms in further. Exact analogue of "zoom to fit both markers" in Google Earth.

```bash
pangenome --db <db> zoom --from <hilbert_idx> --to <hilbert_idx> [--zoom-in <n>]
```

| Option | Default | Description |
|---|---|---|
| `--from` | required | First Hilbert position |
| `--to` | required | Second Hilbert position |
| `--zoom-in` | 0 | Additional zoom-in steps after the initial fit |

```bash
pangenome --db ./human.k15 zoom --from 100000000 --to 500000000 --zoom-in 2
```

---

### `info` Database Statistics

```bash
pangenome --db <db> info
```

Prints genomes, contig lengths, Hilbert min/max ranges, cached voxel count,
and road quality distribution for all cached tiles.

---

### `demo` Built-in Demonstration

Runs a self-contained demo using five short synthetic sequences. Useful for
verifying the installation and understanding the output format.

```bash
pangenome --db /tmp/demo_db demo
```

---

## Database Corruption: Risks and Mitigations

The database has three components with different corruption profiles:

| Component | Corruption risk | Impact |
|---|---|---|
| `registry.json` | Low written once per `add-genome` | Loss of genome metadata only |
| `tile_*.bin` | Low atomic temp-file rename per tile | A corrupt tile affects all voxels in one spatial region |
| `kmer_index.bin` | Low atomic temp-file rename | Loss of alignment capability; rebuild by re-running `add-genome` |

### Risk 1: Interrupted Tile Write

Each tile is written via `tile_*.tmp` > `tile_*.bin` rename. POSIX `rename(2)` is
atomic: either the old `.bin` survives or the new one appears. A crash may leave
a `.tmp` file:

```bash
# Safe to delete never corrupts the corresponding .bin
find ./pangenome_db -name "*.tmp" -delete
```

### Risk 2: Interrupted Registry or Index Write

`registry.json` and `kmer_index.bin` are both written atomically (temp-file
rename). A crash leaves the previous version intact. Back up before large runs:

```bash
cp registry.json registry.json.bak
cp kmer_index.bin kmer_index.bin.bak
```

### Risk 3: Concurrent Writes

Two simultaneous `add-genome` processes against the same `--db` will corrupt
usage counts (the second flush overwrites the first's increments). Always run
`add-genome` sequentially:

```bash
# Correct: sequential
for f in genome_*.fa; do
    pangenome --db ./hprc add-genome --name "$(basename $f .fa)" --fasta "$f"
done

# WRONG will corrupt the database
pangenome --db ./hprc add-genome --name A --fasta A.fa &
pangenome --db ./hprc add-genome --name B --fasta B.fa &
```

### Risk 4: k-mer Length Mismatch

All genomes must use the same k. The tile voxels, kmer index, and registry are
all computed with the k stored in `registry.json`. Changing k after database
creation produces a meaningless map. Use separate directories for different k values.

### Risk 5: Database Built With Old Algorithm

Databases built before the local-hash trajectory was introduced (cumulative walk
`pos += hash.rotate_left(i)`) are incompatible with the current alignment
engine. The Hilbert addresses in the tiles and kmer index were computed with the
old algorithm and will not match the local-hash trajectory used during alignment.
Rebuild from FASTA:

```bash
rm -rf ./old.k15
pangenome --db ./new.k15 add-genome --name T2T-CHM13 \
    --fasta T2T-CHM13.fasta --gff3 T2T-CHM13.gff3 --kmer 15
```

### Risk 6: Filesystem Full During Flush

If the filesystem fills mid-write, the `.tmp` file may be partial. The original
`.bin` is not replaced until the rename succeeds, so it remains intact. Free
space estimate: ~20-80 GB of tile data per genome for the first genome; less for
subsequent genomes sharing the same sequence space.

```bash
df -h ./pangenome_db   # check before large runs
```

---

## Recommended Workflow

### Building a Human Pangenome with Transcript Annotation

```bash
# 1. Add the reference genome with transcript annotation
#    This writes tile voxels + kmer_index.bin for chromosomes + transcripts
pangenome --db ./human.k15 --lru-cache 16384 \
    add-genome --name T2T-CHM13 \
    --fasta T2T-CHM13.fasta \
    --gff3 T2T-CHM13.gff3 \
    --kmer 15

# 2. Add additional assemblies (structural variants create new side streets)
pangenome --db ./human.k15 --lru-cache 16384 \
    add-genome --name HG002_hap1 --fasta HG002.hap1.fasta

# 3. Check road quality distribution
pangenome --db ./human.k15 info

# 4. Align long reads and write SAM + transcript counts
pangenome --db ./human.k15 export-bam \
    --reads sample.fastq.gz \
    --output sample.sam \
    --threads 24 \
    --gff3 T2T-CHM13.gff3

# 5. Convert to BAM
samtools view -bS sample.sam | samtools sort -o sample.bam && samtools index sample.bam
```

### Building a Mouse-Human Pangenome for Xenograft Analysis

```bash
# Build with both species same database, same k
pangenome --db ./xenograft.k15 \
    add-genome --name T2T-CHM13 --fasta T2T-CHM13.fasta --gff3 T2T-CHM13.gff3

pangenome --db ./xenograft.k15 \
    add-genome --name GRCm39 --fasta GRCm39.fasta --gff3 GRCm39.gff3

# Align xenograft reads ZG:Z: tag identifies human vs mouse
pangenome --db ./xenograft.k15 export-bam \
    --reads xenograft.fastq.gz \
    --output xenograft.sam \
    --threads 24

# Split into human and mouse BAMs by ZG tag
samtools view -h xenograft.sam | grep -E "^@|ZG:Z:T2T" | \
    samtools view -bS - > human_reads.bam
samtools view -h xenograft.sam | grep -E "^@|ZG:Z:GRCm39" | \
    samtools view -bS - > mouse_reads.bam
```

---

## Project Structure

```
code/
  hilbert3d/        3D Hilbert curve encode/decode (Skilling 2004 algorithm)
  seqnum/           k-mer hashing, canonical encoding, local-hash trajectory
  grid/             VoxelGrid: LRU-cached disk-backed tile storage; RoadIndex
  streets/          PathStreet abstraction: Hilbert-space street segments
  query/            Hilbert range queries over the VoxelGrid
  pangenome_core/   Genome registry, FASTA/FASTQ loading, alignment engine,
                    kmer index (KmerIndexBuilder + KmerIndex), navigation
  annotation/       GFF3/BED parsing, transcript exon extraction
  bamout/           SAM export with Hilbert GPS custom tags
  pangenome/        CLI entry point
```

---

## Theoretical Background

### Sequence Space as a Metric Space

The approach treats DNA sequence space as a metric space with the following properties:

- **Locality:** The Hilbert curve preserves spatial locality k-mer hashes that are
  numerically close correspond to GPS coordinates that are spatially close. Gene
  families with similar k-mer composition cluster in the same region of the cube.

- **Determinism:** The k-mer hash is a pure function of the sequence. The same
  15-mer always maps to the same GPS address, regardless of which genome, which
  species, or which strand it came from. This makes the cube a universal coordinate
  system for all DNA sequences.

- **Anonymity:** The tile map stores only the geometry of sequence space GPS
  coordinates and usage counts. Individual genome identities are never stored in
  the tile files. This is the OpenStreetMap principle applied to genomics.

- **Accumulation:** Usage counts grow monotonically. Road quality can only increase
  as more genomes are added. The longer the project runs, the more confident the
  road quality labels become.

- **Universality:** Any DNA sequence from any species, any tissue, any sequencing
  technology maps into the same cube using the same k-mer transform. Adding a
  mouse genome populates different regions of the cube, but inter-species conserved
  sequences (ribosomal RNA, core metabolic genes) naturally emerge as shared
  high-usage roads, exactly as in OpenStreetMap where the same physical road is
  validated by drivers from many different cities.

### The Hilbert Curve (Technical)

The implementation uses the Skilling (2004) algorithm for 3D Hilbert curve
encoding and decoding. The curve has BITS=10 bits per axis, giving:

- 1024³ = 2^30 ≈ 1.07 billion addressable cells
- 30-bit Hilbert index space (exactly matching a k=15 canonical k-mer hash)
- Each cell covers a (1/1024)³ cube of normalised sequence space

The 4096 spatial tiles each cover 64³ Hilbert cells. Tile IDs are computed by
decoding the Hilbert index to (x, y, z), scaling to [0, 1024), and dividing by
64 to get integer tile coordinates (tx, ty, tz) in [0, 15]³.

### Why This Is Not a BLAST/BWA Replacement (Yet)

The current implementation does not perform base-level alignment (no Smith-
Waterman, no gap penalties, no base-quality weighting). The CIGAR string is
approximate it reports a match length equal to the read length with soft-
clipped ends for degenerate k-mer runs, but it does not report individual
mismatches or indels at base resolution. The `NM` field counts non-chained
minimizers, not base-level edit distance.

This tool is designed for:
- Rapid species classification and transcript assignment in large read sets
- Building and querying a pangenome road map for structural variant discovery
- Cross-species comparative genomics without reference-specific indexing

For base-level variant calling, use the SAM output as input to `samtools`/
`bcftools` after mapping the Hilbert coordinates back to linear reference
coordinates a `hilbert-to-linear` lookup table is a planned future extension.
