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

### From Sequence to Position

Every DNA sequence can be thought of as a journey through a three-dimensional
space. The tool works as follows:

1. A sliding window of k consecutive bases (a k-mer) is extracted from the sequence.
2. The k-mer is hashed to a canonical integer using a 2-bit encoding (A=0, C=1,
   G=2, T=3), taking the minimum of the forward and reverse-complement hashes so
   that a sequence and its complement land in the same region of space.
3. That integer is used to update a running position using a cumulative walk, much
   like a drunkard's walk through the genome's characteristic sequence space.
4. The resulting 30-bit position is decoded into a three-dimensional (x, y, z)
   coordinate using a Hilbert curve a space-filling curve that maps the
   one-dimensional integer back into 3D while preserving locality: sequences that
   are similar end up close together in space.

### Why a Hilbert Curve and Not a Simple Grid?

A naive grid would scatter similar sequences randomly. The Hilbert curve has a
critical locality property: points that are numerically close on the curve are
also spatially close in the 3D cube. This means:

- Two genome sequences that differ by a single mutation land near each other.
- Repetitive sequences (tandem repeats, transposons) cluster in specific regions
  of the cube, forming a distinctive signature.
- Evolutionary distance between species corresponds roughly to spatial distance
  in the cube.

The implementation uses a 10-bit-per-axis Hilbert curve, giving a 1024 x 1024 x 1024
cube with 2^30 ≈ 1 billion addressable cells. The Skilling (2004) algorithm encodes
and decodes these coordinates exactly and invertibly.

### The Google Earth Analogy

Think of the Hilbert cube as a globe that you can fly over in Google Earth:

| Google Earth concept | Hilbert pangenome concept |
|---|---|
| The globe | The 3D Hilbert cube (genome sequence space) |
| A city | A cluster of similar sequences (e.g. a gene family) |
| A motorway | A Hilbert path traversed by hundreds of genomes |
| A dirt track | A rare variant path seen in only one or two samples |
| Zoom to fit | `zoom` command finds the tightest view around two Hilbert positions |
| GPS coordinates | Normalised (x, y, z) coordinates in the cube |
| Navigate A→B | `route` command finds the optimal path using road quality as weight |
| Saved bookmark | `waypoint` command names a Hilbert position (gene, locus, hotspot) |
| Satellite imagery | Read alignment matches a short sequence to known mapped regions |

---

## The Street Map: How the Database Is Built

### OpenStreetMap, Not a GPS Track Log

This is the key design principle. OpenStreetMap stores streets, not journeys.
When you drive to work, your route validates and potentially improves the road
quality of every street you used, but your personal journey is not stored. Ten
thousand people taking similar commutes make those roads well-attested; a single
car that turned down an unmapped alley adds a new street to the map.

The pangenome database works identically:

- The **tile files** (`tile_*.bin`) store Hilbert cells with usage counts and
  spatial coordinates. They contain no genome names, no sample IDs, no journey
  records.
- The **registry** (`registry.json`) stores per-genome metadata (name, contig
  lengths, Hilbert ranges). This is the equivalent of a bookshelf log saying
  "genome GRCh38 was added" it is separate from the map itself.
- When a new genome is added, its path through the Hilbert cube increments the
  usage count of every cell it visits. Cells visited by many genomes accumulate
  high usage counts and are labelled as high-quality roads. Rare variant cells
  remain low-usage.

### Road Quality Tiers

| Usage count | Road quality | Meaning |
|---|---|---|
| 0 | Desert | Cell exists but has never been traversed |
| 1–9 | Dirt Track | Seen in fewer than 10 genomes |
| 10–49 | Path | Rare but reproducible variant |
| 50–199 | Side Street | Present in a minority of samples |
| 200–999 | Road | Common variant |
| 1,000–4,999 | Main Road | Highly conserved sequence region |
| 5,000–49,999 | Motorway | Near-universal across genomes |
| 50,000+ | Pan-Continent Express | Universally conserved core genome |

These tiers emerge naturally from usage counts without any manual annotation.
Housekeeping genes, ribosomal RNA, and other deeply conserved sequences
automatically become motorways. Private mutations become dirt tracks.

### Chromosomes Are Processed Whole

Each chromosome is read in full, and its Hilbert trajectory is computed in a
single unbroken pass. There are no artificial chunk boundaries that would
introduce false discontinuities into the path a k-mer at position i and a
k-mer at position i+1 always produce consecutive Hilbert cells in an unbroken
road. This matters because a break in the trajectory would appear in the map as
a missing street a gap in a road that should be continuous.

---

## The Tile System: Disk-Based with an LRU Cache

### Why Tiles?

A human genome has ~3 billion bases. With k=15, that produces ~3 billion k-mers,
each mapping to a Hilbert cell. Storing all of these in RAM simultaneously is
infeasible on most machines. The tile system solves this by partitioning the
3D Hilbert cube into spatial tiles (each covering a cube of 64x64x64 cells)
and keeping only a working set of recently used tiles in RAM.

### Write-Back Cache With Merge-on-Flush

The LRU cache is a **write-back** cache:

1. When a Hilbert cell is visited, its tile is loaded from disk into the cache
   if it is not already there. This ensures that prior usage counts from
   previous genome ingestions are present in memory before the new traversal
   is counted.
2. The cell's usage count is incremented in the in-memory tile.
3. After every complete chromosome is processed, all dirty tiles in the cache
   are flushed to disk and evicted, freeing RAM.
4. Flushing writes the in-memory tile directly to disk using an atomic
   temp-file rename (`tile_*.tmp` → `tile_*.bin`). Because prior usage counts
   were loaded during step 1, the flushed tile already contains the fully
   merged state no second read-and-merge at flush time is needed or
   performed (which would double-count usage).

This design means:
- RAM usage is bounded to roughly one chromosome's worth of touched tiles
  plus the LRU working set.
- Every genome ingestion merges cleanly with the existing map.
- Streets discovered by genome A are never erased when genome B is added later.
- Repeats and segmental duplications k-mer sequences that appear on multiple
  chromosomes or in multiple genomes correctly accumulate usage counts each
  time they are traversed.

### Atomic Writes and Crash Safety

Each tile is written using a temp-file rename pattern:

```
write tile_0_3_7_2_.tmp
rename tile_0_3_7_2_.tmp → tile_0_3_7_2_.bin
```

On POSIX filesystems (Linux, macOS) `rename(2)` is atomic with respect to
crashes: either the old file survives or the new one appears, never a
half-written intermediate state. This means a power failure during a flush
leaves at worst one tile in its previous state it does not corrupt the tile.

---

## Building

Requires Rust 1.75+ and Cargo.

```bash
# Build release binary (recommended for real genomes)
cargo build --release -p pangenome

# The binary is at:
./target/release/pangenome
```

---

## Command Reference

All subcommands share a global `--db` path argument pointing to the database
directory. If the directory does not exist it is created on first use.

```
pangenome [--db <path>] [--lru-cache <n>] <subcommand>
```

### Global Options

#### `--db <path>`
Path to the database directory. Default: `./pangenome_db`.

The directory contains:
- `registry.json` genome metadata (names, contig lengths, Hilbert ranges)
- `tile_L_X_Y_Z_.bin` spatial tile files (the actual road map)

```bash
# Use a specific database
pangenome --db /data/human_pangenome info
```

#### `--lru-cache <n>`
Number of tiles to keep in the in-memory LRU cache. Default: 4096.

Each tile is approximately 32 KB when full (2048 voxels x 32 bytes each).
4096 tiles ≈ 128 MB of tile cache. Because tiles are flushed after each
chromosome, the cache size affects speed (more cache = fewer disk reads during
a single chromosome) but not correctness.

```bash
# Use a large cache for fast loading on a 256 GB RAM server
pangenome --lru-cache 65536 add-genome --name GRCh38 --fasta GRCh38.fa

# Use a small cache on a memory-constrained machine (correct but slower)
pangenome --lru-cache 512 add-genome --name GRCh38 --fasta GRCh38.fa
```

---

### `add-genome` Add a Reference Genome

Loads a FASTA file (or inline sequence) and traces its path through the Hilbert
cube, incrementing the usage count of every cell visited. After the full genome
is loaded, all dirty tiles are flushed to disk. The genome's name and contig
metadata are written to `registry.json`.

```
pangenome add-genome --name <name> [--fasta <file>] [--seq <sequence>] [--kmer <k>]
```

**Options:**

| Option | Default | Description |
|---|---|---|
| `--name` / `-n` | required | Logical name for this genome, e.g. `GRCh38`, `CHM13`, `sample_001` |
| `--fasta` / `-f` | | Path to a FASTA file (gzip not currently supported; decompress first) |
| `--seq` / `-s` | | Inline sequence string, for testing or short sequences |
| `--kmer` / `-k` | 15 | k-mer length. Longer k = more specific, shorter k = more sensitive |

**Examples:**

```bash
# Add a human reference genome
pangenome --db ./hprc add-genome --name GRCh38 --fasta /data/GRCh38.fa

# Add a second assembly for comparison
pangenome --db ./hprc add-genome --name CHM13 --fasta /data/CHM13v2.fa

# Add a third shared roads become better attested, private variants stay rare
pangenome --db ./hprc add-genome --name HG002 --fasta /data/HG002.fa

# Quick test with an inline sequence
pangenome add-genome --name test --seq ACGTACGTACGTACGTACGT --kmer 8
```

**Notes:**
- All genomes in the same database must use the same k-mer length. The k value
  is fixed at database creation and stored in `registry.json`. If you add a
  genome with `--kmer 15` to a database created with `--kmer 8`, the new value
  is ignored and the stored value is used.
- For whole human genomes, loading takes several minutes per genome. Memory
  usage during loading is bounded to roughly the size of the largest chromosome
  plus the tile cache.
- The genome name is stored only in `registry.json`, never in the tile files.
  You can safely add genomes under pseudonyms or accession numbers.

---

### `align` Align Reads to the Pangenome

Aligns short or long reads against the map using a hierarchical Hilbert-tile
zoom strategy. This mirrors the way Google Earth narrows down to a location:
first checking coarse tiles, then descending to finer resolution inside
candidate regions.

```
pangenome align [--reads <fasta>] [--read-seq <sequence>] [--read-name <name>]
                [--max-hits <n>] [--min-score <f>]
```

**Options:**

| Option | Default | Description |
|---|---|---|
| `--reads` / `-r` | | FASTA file of reads to align |
| `--read-seq` | | Single inline read sequence |
| `--read-name` | `query_read` | Name label for an inline read |
| `--max-hits` / `-m` | 5 | Maximum alignments to report per read |
| `--min-score` | 0.05 | Minimum alignment score (0–1). Lower values report more, weaker hits |

**Alignment strategy (hierarchical zoom):**

1. The read's k-mer trajectory is computed (same transform as genome loading).
2. At coarse resolution (4 Hilbert bits), tiles that overlap the read trajectory
   are identified.
3. The search descends to medium (7 bits) and fine (10 bits = full resolution)
   resolution, retaining only candidate reference positions that continue to
   match at each level.
4. Each surviving candidate region is scored by spatial divergence (how close
   the read's Hilbert path is to the reference path) and k-mer mismatch count.
5. Forward and reverse-complement orientations are both tested.
6. Hits are ranked by combined score and the top `--max-hits` are returned.

**Output fields:**

- `GPS_start` / `GPS_end` normalised (x, y, z) coordinates in the Hilbert cube
- `Hilbert=[start..end]` the raw Hilbert index range of the aligned region
- `CIGAR` trajectory-length M operations (e.g. `150M` for a 150 bp read)
- `MAPQ` Phred-scaled mapping quality (0–60)
- `NM` approximate edit distance (k-mer mismatch count)

**Examples:**

```bash
# Align a FASTA of reads
pangenome --db ./hprc align --reads reads.fa --max-hits 3

# Align a single sequence
pangenome --db ./hprc align --read-seq ACGTACGTACGT --read-name my_probe

# Report all hits above a high-confidence threshold
pangenome --db ./hprc align --reads reads.fa --min-score 0.8 --max-hits 1
```

---

### `export-bam` Export Alignments as SAM

Runs the alignment engine and writes results to a SAM file that can be
converted to BAM with samtools.

```
pangenome export-bam [--reads <fasta>] [--read-seq <sequence>] [--read-name <name>]
                     [--output <file>] [--show-bam-cmd]
```

**Options:**

| Option | Default | Description |
|---|---|---|
| `--reads` / `-r` | | FASTA file of reads |
| `--read-seq` | | Single inline read |
| `--read-name` | `query_read` | Name for inline read |
| `--output` / `-o` | `alignments.sam` | Output SAM file path |
| `--show-bam-cmd` | false | Print the samtools command to convert SAM → BAM |

**Example:**

```bash
pangenome --db ./hprc export-bam --reads reads.fa --output hits.sam --show-bam-cmd
# prints: samtools view -bS hits.sam | samtools sort -o hits.bam && samtools index hits.bam
```

---

### `annotate` Load and Query Genomic Annotations

Loads GFF3 or BED annotation files, maps each feature to its Hilbert position,
and queries features near a given coordinate.

```
pangenome annotate [--gff3 <file>] [--bed <file>]
                   [--seqname <chr>] [--start <n>] [--end <n>]
                   [--hilbert <idx>] [--radius <n>]
                   [--summary]
```

**Options:**

| Option | Default | Description |
|---|---|---|
| `--gff3` | | GFF3 annotation file |
| `--bed` | | BED annotation file |
| `--seqname` | | Chromosome/contig name for linear coordinate query |
| `--start` | | Linear start (0-based) |
| `--end` | | Linear end (0-based, exclusive) |
| `--hilbert` | | Query near a specific Hilbert index instead of linear coordinates |
| `--radius` | 1,000,000 | Search radius in Hilbert index units for `--hilbert` queries |
| `--summary` | false | Show feature-type counts only, not individual features |

**Examples:**

```bash
# Load a GFF3 and summarise feature types
pangenome annotate --gff3 GRCh38.gff3 --summary

# Query features near a linear coordinate range
pangenome annotate --gff3 GRCh38.gff3 --seqname chr1 --start 1000000 --end 2000000

# Query features near a Hilbert position (e.g. from an alignment hit)
pangenome annotate --gff3 GRCh38.gff3 --hilbert 823456789 --radius 500000
```

---

### `route` Plan a Route Through the Map

Finds the optimal path between two Hilbert positions using A* search, weighted
by road quality. High-quality (high-usage) roads are preferred like a sat-nav
choosing motorways over dirt tracks.

```
pangenome route --from <hilbert_idx> --to <hilbert_idx> [--max-steps <n>]
```

**Options:**

| Option | Default | Description |
|---|---|---|
| `--from` | required | Start Hilbert index |
| `--to` | required | End Hilbert index |
| `--max-steps` | 20 | Maximum route steps to display (full path is still computed) |

**Example:**

```bash
# Route between two positions from an alignment result
pangenome --db ./hprc route --from 100000000 --to 200000000

# Output includes GPS coordinates and road quality at each step:
#   0: H= 100000000  GPS=(0.1234, 0.5678, 0.2345)  road=Motorway   cost=0.1000
#   1: H= 100032768  GPS=(0.1235, 0.5679, 0.2346)  road=MainRoad   cost=0.3000
#   ...
```

The route score has two components:
- **Spatial score** how directly the path moves through 3D space (1.0 = straight line)
- **Road quality score** average usage tier of cells on the path (1.0 = all motorway)
- **Total score** 60% spatial + 40% road quality

---

### `waypoint` Named Markers

Add and query named waypoints like Google Earth bookmarks for gene loci,
mutation hotspots, or regulatory elements.

```
pangenome waypoint add --name <name> --hilbert <idx> [--tag key=value ...]
pangenome waypoint list [--tag key=value]
pangenome waypoint near --hilbert <idx> [--radius <n>]
```

**Examples:**

```bash
# Mark a gene locus
pangenome --db ./hprc waypoint add --name BRCA1 --hilbert 823456789 \
    --tag type=gene --tag disease=breast_cancer

# List all waypoints tagged as genes
pangenome --db ./hprc waypoint list --tag type=gene

# Find all waypoints within 1 million Hilbert units of a position
pangenome --db ./hprc waypoint near --hilbert 823456789 --radius 1000000
```

Note: waypoints are currently held in memory for the duration of a session and
are not persisted to the database between invocations. To persist them, add
support for a `waypoints.json` file alongside `registry.json` this is a
planned future extension.

---

### `zoom` Zoom to Fit (Google Earth Style)

Finds the tightest Hilbert tile that contains both a start and end position,
then optionally zooms in further. This is the exact analogue of pressing
"zoom to fit both markers" in Google Earth.

```
pangenome zoom --from <hilbert_idx> --to <hilbert_idx> [--zoom-in <n>]
```

**Options:**

| Option | Default | Description |
|---|---|---|
| `--from` | required | First Hilbert position |
| `--to` | required | Second Hilbert position |
| `--zoom-in` | 0 | Number of additional zoom-in steps after the initial fit |

**Example:**

```bash
pangenome --db ./hprc zoom --from 100000000 --to 500000000 --zoom-in 2

# Output:
#   Zoom L0: level_bits=4  voxels=1024  tile_prefix=...
#      BBox: (0.098, 0.102, 0.201) → (0.489, 0.491, 0.600)
#   Zoom L1: level_bits=5  voxels=512   tile_prefix=...
#   Zoom L2: level_bits=6  voxels=128   tile_prefix=...
```

`level_bits` increases as you zoom in. At `level_bits=10` (maximum) you are
at single-cell resolution in the Hilbert cube.

---

### `info` Database Statistics

Shows the genomes stored, their contigs, and the road quality distribution
across all cached tiles.

```
pangenome info
```

**Example output:**

```
=== Pangenome Database: ./hprc ===
k-mer length: 15
Genomes: 3

  [GRCh38]  3099734149 bases  25 contigs  k=15
    chr1                           len= 248956422  H=[  82345678.. 923456789]
    chr2                           len= 242193529  H=[ 123456789.. 834567890]
    ...

Cached voxels: 847293

Road quality distribution:
  Dirt Track                34821
  Path                      102847
  Side Street               389201
  Road                      201843
  Main Road                  89234
  Motorway                   27891
  Pan-Continent Express       1456
```

Note: "Cached voxels" reflects only what is currently in the LRU cache. The
full database on disk may contain many more tiles. Use `find <db> -name "tile_*.bin" | wc -l`
to count all tile files on disk.

---

### `demo` Built-in Demonstration

Runs a self-contained demo using five short synthetic sequences. Useful for
verifying the installation and understanding the output format without needing
real genome files.

```
pangenome demo
```

---

## Database Corruption: Risks and Mitigations

### Understanding What Can Go Wrong

The database has two parts with different corruption profiles:

| Component | Corruption risk | Impact |
|---|---|---|
| `registry.json` | Low written once per `add-genome` call | Loss of genome metadata only; tile data is unaffected |
| `tile_*.bin` | Low for individual tiles; atomic write protects each | A corrupt tile affects all voxels in that spatial region |

### Risk 1: Interrupted Tile Write (Power Failure, Kill Signal)

**What happens:** Each tile is written by renaming a `.tmp` file over the `.bin`
file. On POSIX systems this rename is atomic. A crash during the write leaves
either the old `.bin` intact or the new one in place never a partial write.

**Residual risk:** If the process is killed after writing `tile_*.tmp` but
before the rename, a `.tmp` file may be left in the database directory.

**Mitigation:** On startup (or as a maintenance step), delete any `tile_*.tmp`
files in the database directory:

```bash
find ./pangenome_db -name "tile_*.tmp" -delete
```

These are always safe to delete they represent an uncommitted write that
has not replaced the current `.bin`.

### Risk 2: Interrupted Registry Write

**What happens:** `registry.json` is written with `std::fs::write`, which is
not atomic on all platforms. A crash during this write could leave a truncated
or malformed JSON file.

**Mitigation:** Back up `registry.json` before large ingestion runs:

```bash
cp ./pangenome_db/registry.json ./pangenome_db/registry.json.bak
pangenome add-genome --name new_genome --fasta new.fa
```

If `registry.json` is corrupted but tile files are intact, the tile data is
not lost. You can reconstruct `registry.json` by re-running `add-genome` for
each FASTA with `--name` matching the original names the tile usage counts
will simply be incremented again (doubling road quality for those genomes).
A future version will use an atomic write for the registry.

### Risk 3: Concurrent Writes to the Same Database

**What happens:** Two `add-genome` processes running simultaneously against
the same database directory will both read the same tile files, increment
usage counts independently in memory, and then both flush. The second flush
will overwrite the first, losing the first process's increments.

**Mitigation:** Never run two `add-genome` commands against the same `--db`
directory simultaneously. Process genomes sequentially:

```bash
# Correct: sequential
for f in genome_*.fa; do
    name=$(basename $f .fa)
    pangenome --db ./hprc add-genome --name "$name" --fasta "$f"
done

# WRONG: parallel runs to the same database will corrupt usage counts
pangenome --db ./hprc add-genome --name A --fasta A.fa &
pangenome --db ./hprc add-genome --name B --fasta B.fa &   # DO NOT DO THIS
wait
```

If you need to exploit parallelism, build separate databases for subsets of
genomes and merge them (a future `merge-db` command is planned).

### Risk 4: k-mer Length Mismatch

**What happens:** All genomes in a database must use the same k-mer length.
If you create a database with k=15 and then attempt to add a genome with k=8,
the stored k=15 is used silently. If you delete `registry.json` and re-create
it with a different k, the tile files contain data computed with the old k and
the new k-mer hashes will map to completely different positions in the cube —
producing a meaningless map.

**Mitigation:** Never change `--kmer` on an existing database. The k value is
fixed at creation. If you need to experiment with different k values, use
separate database directories.

```bash
pangenome --db ./hprc_k15 add-genome --name GRCh38 --fasta GRCh38.fa --kmer 15
pangenome --db ./hprc_k8  add-genome --name GRCh38 --fasta GRCh38.fa --kmer 8
# The two databases are completely independent and incompatible with each other.
```

### Risk 5: Filesystem Full During Flush

**What happens:** If the filesystem fills during a tile write, `std::fs::write`
returns an error. The `.tmp` file may be partially written. The original `.bin`
is not touched until the rename succeeds, so it remains intact.

**Mitigation:** Ensure sufficient free space before large ingestion runs.
A rough estimate: one human genome at k=15 produces roughly 3 billion k-mers.
With 32 bytes per voxel and assuming moderate deduplication across tiles,
expect 20–80 GB of tile data per genome for the first genome added; subsequent
genomes sharing the same sequence space add proportionally less.

```bash
# Check free space before a run
df -h ./pangenome_db
```

---

## Recommended Workflow for a Human Pangenome

```bash
# 1. Create the database with the first reference genome
pangenome --db ./hprc --lru-cache 16384 \
    add-genome --name GRCh38 --fasta GRCh38.fa --kmer 15

# 2. Add additional assemblies sequentially
pangenome --db ./hprc --lru-cache 16384 \
    add-genome --name CHM13v2 --fasta CHM13v2.fa

pangenome --db ./hprc --lru-cache 16384 \
    add-genome --name HG002_hap1 --fasta HG002.pat.fa

# 3. Check road quality distribution
pangenome --db ./hprc info

# 4. Align a set of reads
pangenome --db ./hprc align --reads my_reads.fa --max-hits 5 --min-score 0.1

# 5. Export alignments for downstream tools
pangenome --db ./hprc export-bam --reads my_reads.fa \
    --output alignments.sam --show-bam-cmd

# 6. Navigate between two known loci
pangenome --db ./hprc waypoint add --name BRCA1 --hilbert 823456789
pangenome --db ./hprc waypoint add --name TP53  --hilbert 412345678
pangenome --db ./hprc route --from 823456789 --to 412345678
```

---

## Project Structure

```
code/
  hilbert3d/        3D Hilbert curve encode/decode (Skilling algorithm)
  seqnum/           k-mer hashing, canonical encoding, trajectory building
  grid/             VoxelGrid: LRU-cached, disk-backed tile storage with merge
  streets/          PathStreet abstraction: Hilbert-space street segments
  query/            Hilbert range queries over the VoxelGrid
  pangenome_core/   Genome registry, FASTA loading, read alignment, navigation
  annotation/       GFF3/BED loading and feature queries
  bamout/           SAM/BAM export
  pangenome/        CLI entry point
```

---

## Theoretical Background

The approach treats DNA sequence space as a metric space with the following
properties:

- **Locality:** Similar k-mer compositions map to nearby positions. The Hilbert
  curve's locality-preservation property means that genomic neighbourhoods
  (exons, gene families, repeat classes) cluster spatially.
- **Anonymity:** The map stores only the geometry of sequence space, not the
  identity of individual genomes. This is analogous to how a physical road
  exists independent of who has driven on it.
- **Accumulation:** Usage counts grow monotonically. A voxel's road quality can
  only increase or stay the same as more genomes are added. This gives the map
  a natural directionality: the longer the project runs, the more confident the
  road quality labels become.
- **Universality:** Any DNA sequence from any species, any tissue, any
  technology can be mapped into the same cube using the same k-mer transform.
  The cube is not human-specific. Adding a mouse genome, a bacterial genome,
  or a synthetic construct simply populates different regions of the cube, and
  inter-species conserved sequences (ribosomal RNA, core metabolic genes)
  naturally emerge as shared roads.
