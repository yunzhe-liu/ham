# HAM: Hash Accelerated Matcher

**Fast, accurate sgRNA guide quantification for single-cell Perturb-seq data.**

HAM replaces k-mer pseudoalignment with a **window-restricted exact hash lookup**,
eliminating adapter-derived false-positive matches. It processes raw FASTQ pairs
and produces deduplicated MEX count matrices — directly comparable to existing
tools (simpleaf, Cell Ranger) but faster and more accurate at the per-cell level.

---

## Quick Start

```bash
# 1. Build guide hash table (one-time)
ham build-hash --fasta guides.fasta --output guide_hash.pkl

# 2. Match reads to guides
ham match \
    -1 lane01_sgRNA_R1.fastq.gz \
    -2 lane01_sgRNA_R2.fastq.gz \
    -w gex_whitelist.txt \
    -g guide_hash.pkl \
    -o hits.npz \
    -t 4

# 3. UMI deduplication + MEX matrix
ham dedup -i hits.npz -o mex_output/

# 4. Merge per-lane matrices
ham merge --lanes lane_list.tsv --out merged/ --prefix merged
```

### Python API

```python
from ham import build_guide_hash, load_whitelist, load_guide_hash
from ham.matcher import match_reads
from ham.dedup import build_count_matrix

guide_hash = build_guide_hash("guides.fasta", "guide_hash.pkl")
whitelist = load_whitelist("whitelist.txt")
match_reads("-1", "lane01_R1.fastq.gz", "-2", "lane01_R2.fastq.gz",
            whitelist=whitelist, guide_hash=guide_hash,
            output_path="hits.npz", threads=4)
build_count_matrix("hits.npz", "mex_output/")
```

### Expected Throughput

| Scenario | Speed | Memory |
|----------|:-----:|:------:|
| 4,500 guides × 3.7M whitelist, 8 FASTQ pairs, 4 threads | ~2M reads/s | 700 MB (dict) or 110 MB (low-memory) |
| 48-lane Perturb-seq (1.0 × 10⁹ read pairs) | ~8 minutes per lane | As above |

---

## Algorithm

### Overview

HAM operates in four stages, each with a single CLI command:

```
FASTQ pairs  ──→  [1. match]  ──→ hits.npz  ──→  [2. dedup]  ──→ MEX trio
guides.fasta ──→  [0. build-hash] ──→ guide_hash.pkl                │
whitelist.txt ──────────────────────────────────────────────────────┘
                                                                     │
                                          [3. merge] ← per-lane MEX ┘
```

The algorithm is built on three core ideas:

1. **2-bit encoding** — All DNA sequences (barcodes, UMIs, guides) are encoded as integers using A=00, C=01, G=10, T=11. A 16 bp barcode becomes a 32-bit integer; a 20 bp guide becomes a 40-bit integer. All comparisons are integer operations — no string parsing during the matching loop.

2. **Window-restricted extraction** — Rather than aligning the full Read2 sequence, HAM extracts a single 26 bp window (positions 28–54 of Read2) anchored against the known sgRNA construct layout. The guide sequence resides within this window. Seven sliding 20 bp sub-windows are extracted and checked against the hash table. This eliminates adapter-derived k-mer noise because the search space is physically constrained to the guide region.

3. **Pre-computed hash tables with Hamming=1 expansion** — Both the guide hash table and the cell barcode whitelist are pre-expanded to include all valid Hamming-distance-1 neighbours. During matching, a single O(1) hash lookup finds the parent guide or barcode — no iterative alignment, no scoring matrices, no seeding-and-extending.

### Stage 0: `build-hash` — Guide Hash Table Construction

**Input:** Guide FASTA file (4,532 single-guide sequences of 20 bp each).

**Output:** Python pickle dictionary with three key structures:

| Key | Type | Content |
|-----|------|---------|
| `seq_to_idx` | `dict[str, int]` | Exact guide sequence → guide index (string keys, for debugging) |
| `seq_to_idx_int` | `dict[int, int]` | Integer-encoded guide sequence → guide index (used in matching loop) |
| `idx_to_id` | `list[str]` | Guide index → guide ID string (e.g. `AAAS_-_53715438.23-P1P2`) |

**Construction steps:**

1. Parse the FASTA file, extracting guide ID (header after `>`) and 20 bp sequence.
2. For each unique guide sequence, assign an integer index and insert into both `seq_to_idx` (string key) and `seq_to_idx_int` (integer key).
3. For each guide sequence, generate all **Hamming-distance-1 variants**: for each of the 20 positions, substitute the base with the 3 alternatives (20 × 3 = 60 variants per guide). Each variant is also inserted into the hash table, pointing to the same parent guide index. Variants that collide with an existing exact guide sequence are discarded (the exact match takes priority).

**Result:** ~4,500 exact + ~268,000 Hamming-1 neighbour entries = ~275,000 total hash entries. Lookup is O(1) via Python dict. Serialised to pickle (~10 MB).

### Stage 1: `match` — Read Pair Processing

**Input:** R1 and R2 FASTQ files (gzip-compressed), cell barcode whitelist, guide hash table.

**Output:** `hits.npz` containing:
- `hits`: (N, 3) int32 array — columns `[cb_idx, umi_int, guide_idx]`
- `barcode_list`: (M,) string array — index → cell barcode string
- `idx_to_id`: (K,) string array — index → guide ID string

**Processing pipeline (per read pair):**

#### Step 1a: 2-Bit Integer Encoding

All DNA processing uses 2-bit-per-base integer encoding via pre-computed 256-entry lookup tables (`_BYTE2BITS`). This avoids string construction at every step:

| Sequence | Length | Representation | Bits |
|----------|:------:|----------------|:----:|
| Cell barcode | 16 bp | `uint32` | 32 |
| UMI | 12 bp | `uint32` | 24 |
| Guide | 20 bp | `uint64` | 40 |
| Window (R2[28:54]) | 26 bp | `uint64` | 52 |

#### Step 1b: Cell Barcode Correction (Hamming ≤ 1)

The 10x barcode whitelist (3.7M entries) is pre-expanded into a **CB hash table** containing all valid barcodes plus all their Hamming-distance-1 variants (up to 48 variants per barcode). Two modes are available:

- **Plan H (default, dict):** Python `dict[int, int]` mapping encoded 32-bit barcode integer → parent barcode index. O(1) lookup, ~700 MB. Suitable for production servers.
- **Plan G (low-memory, `--low-memory`):** Sorted `numpy` int32 arrays with `np.searchsorted` binary search. O(log N) lookup, ~110 MB. Suitable for memory-constrained environments.

During matching, the 16 bp barcode at Read1[0:16] is encoded as a `uint32` and looked up in the CB hash. If found, it is corrected to its parent whitelist barcode. Reads with barcodes not found in the hash are discarded. By default, `--cb-max-hamming 1` limits correction to Hamming ≤ 1. Use `--cb-max-hamming 2` for chemistry mismatches (e.g., sequencing Nextera-formatted barcodes against a TruSeq whitelist).

#### Step 1c: UMI Encoding

The 12 bp UMI at Read1[16:28] is encoded as a `uint32` using a loop-unrolled, 12-element LUT lookup — each byte is converted via `_BYTE2BITS` and shifted into position. The encoded UMI is stored in the hits array; string decoding is deferred to the deduplication stage.

#### Step 1d: Window-Restricted Guide Extraction

The critical design decision: **only positions 28–54 of Read2 are examined.** This window is anchored against the known 10x sgRNA construct layout:

```
Read2 layout (3' v3 chemistry, sgRNA library):
  [0:10]   — 10x TSO adapter remnant (discarded)
  [10:30]  — poly-G / linker (discarded)
  [28:48]  — guide sequence (20 bp)        ← target
  [48:54]  — poly-A tail start (6 bp)      ← provides flanking context
  [54:90]  — poly-A tail remainder (discarded)
```

The 26 bp window (positions 28–54) is encoded as a 52-bit integer using five 5 bp chunk LUT lookups for performance:

```
Window: [b0 b1 b2 b3 b4] [b5 b6 b7 b8 b9] ... [b25]
         └── 10-bit ──┘  └── 10-bit ──┘       └─ 2-bit
```

Seven 20 bp sub-windows (sliding by 1 bp, offsets 0–6) are extracted via simple bit-shift-and-mask operations on the 52-bit integer:

```python
MASK_40 = (1 << 40) - 1
for offset in range(7):
    guide_candidate = (window_int >> shift) & MASK_40
```

Each candidate 40-bit integer is checked against `guide_hash_int`:

- **Fast path (exact match):** O(1) dict lookup. Handles >99% of guide matches.
- **Slow path (Hamming=1):** For the ~1% of reads with a single sequencing error, generates all 60 Hamming-1 variants of the 20 bp guide on-the-fly and checks each against the hash. This is also O(1) per variant (60 × O(1) = constant overhead).

The first matching guide (exact or Hamming-1) is recorded. If none of the 7 windows match, the read pair is counted as `no_guide` and discarded.

#### Step 1e: Multi-File Parallelism

When multiple FASTQ pairs are provided (comma-separated lists to `-1` and `-2`), HAM uses `multiprocessing.Pool` to process each pair in parallel. The guide hash and CB hash are shared read-only via `fork()` copy-on-write semantics.

### Stage 2: `dedup` — UMI Directional Deduplication + MEX Generation

**Input:** `hits.npz` from the match stage.

**Output:** Standard MEX trio (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`).

#### Step 2a: Grouping

Hits are grouped by `(cell_barcode_idx, guide_idx)`. Each group contains a list of encoded UMI integers for that (cell, guide) pair. The encoded UMIs are decoded to 12 bp strings via `decode_umi()` (reverse of the loop-unrolled LUT encoding).

#### Step 2b: UMI-tools Directional Deduplication

For each (cell, guide) group, UMIs are deduplicated using the **UMI-tools directional algorithm**:

1. Count occurrences of each unique UMI string.
2. Sort UMIs by descending count (highest-count UMI first).
3. For each UMI in descending order:
   - If it has already been masked by a higher-count UMI, skip it.
   - Otherwise, retain it and mark all its Hamming-distance-1 neighbours as masked.

The original UMI-tools algorithm requires O(n × m) linear scanning (where n is the number of distinct UMIs and m is the number of Hamming-1 neighbours). HAM replaces this with an **O(n × 36) hash-accelerated implementation**: for each retained UMI, its 36 Hamming-1 variants (12 positions × 3 alternative bases) are generated and inserted into a local hash set. Subsequent UMIs are checked against this set in O(1) — the same algorithmic result at 10–100× the speed.

The deduplicated count for each (cell, guide) pair is the number of retained UMIs.

#### Step 2c: Sparse Matrix Construction

Deduplicated counts are assembled into a `scipy.sparse.csr_matrix` with dimensions (n_cells × n_guides) and written as a gzip-compressed Market Exchange Format (MEX) file. Only cells and guides with ≥1 retained UMI are included.

### Stage 3: `merge` — Cross-Lane Merging

**Input:** Lane list TSV (`lane_id\tmatrix_dir\tsuffix`) + per-lane MEX directories.

**Output:** Merged MEX trio with lane-suffixed barcodes (`AAACCCAAG...-L01`).

Per-lane sparse matrices are vertically stacked via `scipy.sparse.vstack`. The lane suffix (e.g. `-L01`) is appended to each barcode to ensure cross-lane uniqueness. Features (columns) must be identical across lanes — the feature set from the first lane serves as the reference.

---

## Design Rationale — Why Hash Lookup Beats k-mer Pseudoalignment

### The Problem with k-mer Approaches

k-mer pseudoalignment tools (piscem, kallisto) map reads by decomposing them into short k-mers (typically k=13 or 15) and looking up each k-mer in a pre-built coloured de Bruijn graph index. For short guide constructs on a 10x 3' v3 backbone:

- Read2 contains a TSO adapter remnant followed by a poly-G linker, then the guide itself.
- k-mers spanning the poly-G / TSO boundary can accidentally match genuine guide k-mers because poly-G and the adapter contain short subsequences shared with guide sequences (e.g., `GGGG`, `AAGCAG`).
- This produces **adapter-derived false positives**: reads where the true guide is absent but k-mers from the adapter region pseudoalign to a guide in the index.

### HAM's Solution

1. **Window-restricted extraction eliminates adapter noise entirely.** The search is confined to positions 28–54 — the adapter and poly-G linker (positions 0–27) are never examined for guide content.

2. **Integer encoding eliminates string allocation.** The entire matching loop operates on fixed-width integers. A 20 bp guide comparison is a single 64-bit dict lookup, not a string comparison.

3. **Pre-computed Hamming-1 expansion handles sequencing errors at zero runtime cost.** Since all valid single-base-error variants are already in the hash table, the fast path (exact match) catches >99% of reads with a single O(1) dict lookup. The slow path is only invoked for the rare reads where the error falls outside the pre-computed set (e.g., combined adapter-truncation + sequencing error).

---

## CLI Reference

```
ham build-hash      Build guide sequence hash table from FASTA
  --fasta PATH        Guide FASTA file
  --output PATH       Output hash table (.pkl)

ham match           Match sgRNA reads to guide reference
  -1, --read1 PATH    R1 FASTQ (CB+UMI), comma-separated for multi-file
  -2, --read2 PATH    R2 FASTQ (guide insert)
  -w, --whitelist PATH  Cell barcode whitelist (one per line)
  -g, --guide-hash PATH  Guide hash table (.pkl)
  -o, --output PATH     Output hits file (.npz)
  -t, --threads N     Parallel workers (one per FASTQ pair) [default: 1]
  --low-memory        Use numpy binary-search CB hash (~110 MB vs ~700 MB)
  --cb-max-hamming N  Max Hamming distance for CB correction [default: 1]

ham dedup           UMI deduplication + MEX matrix generation
  -i, --input PATH    Hits file (.npz) from ham match
  -o, --output-dir PATH  Output directory for MEX files
  -t, --umi-threshold N  UMI Hamming distance threshold [default: 1]

ham merge           Merge per-lane count matrices
  --lanes PATH        Lane list TSV: lane_id, matrix_dir, suffix
  --out PATH          Output directory
  --prefix PREFIX     Output file prefix [default: merged]
```

---

## Installation

```bash
pip install -e /path/to/ham
```

Requires Python ≥ 3.10, numpy, scipy. The package exposes a `ham` CLI entry point and importable Python modules (`ham.matcher`, `ham.dedup`, `ham.merge`, `ham.hash_builder`, `ham.encoding`).

