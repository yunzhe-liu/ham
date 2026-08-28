# HAM: Hash Accelerated Matcher

**Fast, accurate sgRNA guide quantification for single-cell Perturb-seq data.**

HAM replaces k-mer pseudoalignment with a **window-restricted exact hash lookup**,
eliminating adapter-derived false-positive matches. It processes raw FASTQ pairs
and produces deduplicated MEX count matrices — directly comparable to existing
tools (simpleaf, Cell Ranger) but faster and more accurate at the per-cell level.

HAM is also embedded in the [scprocess-perturb](https://github.com/yunzhe-liu/scprocess-perturb)
Snakemake workflow, which provides a unified chemistry-configuration layer and
orchestrates the full guide-extraction pipeline.

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

---

## Installation

```bash
git clone https://github.com/yunzhe-liu/ham.git
cd ham
pip install -e .
```

Requires Python ≥ 3.10, numpy, scipy. The package exposes a `ham` CLI entry
point and importable Python modules (`ham.matcher`, `ham.dedup`, `ham.merge`,
`ham.hash_builder`, `ham.encoding`).

---

## Supported 10x Chemistries

HAM supports three named chemistries plus a `custom` escape hatch. All physical
read-layout parameters are resolved from the chemistry selection.

| Chemistry | Matches these 10x kits | R1 layout | UMI | R2 guide window | Guide length | Default whitelist |
|:---|:---|:---|:---:|:---|:---:|:---|
| `10xv3` (default) | 3' v3/v3.1/v4, 3LT, Multiome | 28bp (16CB+12UMI) | 12bp | pos 28–54 (26bp, ±3bp margin around the 20bp guide) | 20bp | 3M-feb-2018 / 3M-3pgex-may-2023 |
| `10xv2-5p` | 5' v1.0, 5' v2 | 26bp (16CB+10UMI) | 10bp | pos 13–39 (26bp, ±3bp margin around the 20bp guide) | 20bp | 737K-aug-2016 |
| `10xv2-5p-12umi` | 5' v3 (GEM-X) | 28bp (16CB+12UMI) | 12bp | pos 13–39 (26bp, ±3bp margin around the 20bp guide) | 20bp | 3M-5pgex-jan-2023 |
| `custom` | Any non-standard hardware | user-defined | user-defined | user-defined | user-defined | user-provided |

The guide anchor position within R2 (position 16 for the 5' chemistries,
position 28 for `10xv3`) is validated against real data for each
chemistry — it is not derived from a generic construct diagram. Guide
length is 20bp for all three named chemistries (standard SpCas9
protospacer length; the 5' chemistries previously used an incorrect 19bp
with zero margin). All three windows now carry a ±3bp margin around their
validated anchor (`window_len - guide_len + 1` = 7 candidate offsets are
scanned automatically). If your own library prep puts the guide somewhere
else in R2 entirely, don't assume a named chemistry fits — use
`--chemistry custom` with your own `--window-start`/`--window-end`/
`--guide-len` instead.

```bash
# Standard chemistries
ham match ... --chemistry 10xv3            # 3' v3/v4, 3LT, multiome
ham match ... --chemistry 10xv2-5p         # 5' v1/v2
ham match ... --chemistry 10xv2-5p-12umi   # 5' v3 (GEM-X)

# Custom chemistry — pass each position explicitly
ham match ... --chemistry custom \
    --cb-start 0 --cb-end 14 \
    --umi-start 14 --umi-end 22 \
    --window-start 8 --window-end 28 \
    --guide-len 20
```

The `ham build-hash` and `ham dedup` stages are chemistry-independent (the
`--umi-len` flag on `ham dedup` controls UMI decoding length only).

---

## Python API

```python
from ham import build_guide_hash, load_whitelist, load_guide_hash
from ham.matcher import match_reads
from ham.dedup import build_count_matrix

guide_hash = build_guide_hash("guides.fasta", "guide_hash.pkl")
whitelist = load_whitelist("whitelist.txt")

# Standard chemistry
match_reads("lane01_R1.fastq.gz", "lane01_R2.fastq.gz",
            whitelist=whitelist, guide_hash=guide_hash,
            output_path="hits.npz", threads=4, chemistry="10xv3")
build_count_matrix("hits.npz", "mex_output/", umi_len=12)

# Custom chemistry
match_reads("lane01_R1.fastq.gz", "lane01_R2.fastq.gz",
            whitelist=whitelist, guide_hash=guide_hash,
            output_path="hits.npz", threads=4,
            chemistry="custom",
            chem_cfg={
                "cb_start": 0, "cb_end": 14,
                "umi_start": 14, "umi_end": 22,
                "window_start": 8, "window_end": 28,
                "guide_len": 20,
            })
build_count_matrix("hits.npz", "mex_output/", umi_len=8)
```

---

## Expected Throughput

HAM processes ~2 million read pairs per second per thread on typical hardware.
A genome-scale Perturb-seq dataset (~1 billion read pairs across 48 lanes)
completes in under 10 minutes per lane. Memory usage is ~700 MB in default
mode (Python dict CB hash) or ~110 MB with `--low-memory` (numpy binary-search
CB hash).

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

1. **2-bit encoding** — All DNA sequences (barcodes, UMIs, guides) are
   encoded as integers using A=00, C=01, G=10, T=11. A 16 bp barcode becomes
   a 32-bit integer; a 20 bp guide becomes a 40-bit integer. All comparisons
   are integer operations — no string parsing during the matching loop.

2. **Window-restricted extraction** — Rather than aligning the full Read2
   sequence, HAM extracts a single window anchored against the known sgRNA
   construct layout. The guide sequence resides within this window. Seven
   sliding sub-windows are extracted and checked against the hash table. This
   eliminates adapter-derived k-mer noise because the search space is
   physically constrained to the guide region.

3. **Pre-computed hash tables with Hamming=1 expansion** — Both the guide
   hash table and the cell barcode whitelist are pre-expanded to include all
   valid Hamming-distance-1 neighbours. During matching, a single O(1) hash
   lookup finds the parent guide or barcode — no iterative alignment, no
   scoring matrices, no seeding-and-extending.

### Stage 0: `build-hash` — Guide Hash Table Construction

**Input:** Guide FASTA file.

**Output:** Python pickle dictionary with three key structures:

| Key | Type | Content |
|-----|------|---------|
| `seq_to_idx` | `dict[str, int]` | Exact guide sequence → guide index (string keys, for debugging) |
| `seq_to_idx_int` | `dict[int, int]` | Integer-encoded guide sequence → guide index (used in matching loop) |
| `idx_to_id` | `list[str]` | Guide index → guide ID string |

**Construction:**

1. Parse the FASTA file, extracting guide ID (header after `>`) and sequence.
2. For each unique guide sequence, assign an integer index and insert into
   both `seq_to_idx` (string key) and `seq_to_idx_int` (integer key).
3. Generate all **Hamming-distance-1 variants** (20 positions × 3 alternative
   bases = 60 variants per guide). Each variant is also inserted into the hash
   table, pointing to the same parent guide index. Variants that collide with
   an existing exact guide sequence are discarded (exact match takes priority).

### Stage 1: `match` — Read Pair Processing

**Input:** R1 and R2 FASTQ files (gzip-compressed), cell barcode whitelist,
guide hash table.

**Output:** `hits.npz` containing:
- `hits`: (N, 3) int32 array — columns `[cb_idx, umi_int, guide_idx]`
- `barcode_list`: (M,) string array — index → cell barcode string
- `idx_to_id`: (K,) string array — index → guide ID string

**Processing per read pair:**

**Step 1a — 2-Bit Integer Encoding:** All DNA processing uses 2-bit-per-base
encoding via pre-computed 256-entry lookup tables (`_BYTE2BITS`). This avoids
string construction at every step.

**Step 1b — Cell Barcode Correction (Hamming ≤ 1):** The 10x barcode whitelist
is pre-expanded into a CB hash table containing all valid barcodes plus all
their Hamming-distance-1 variants. Two modes:
- **Plan H (default, dict):** O(1) lookup, ~700 MB. Production use.
- **Plan G (low-memory):** Sorted numpy arrays + binary search. O(log N)
  lookup, ~110 MB.

Reads with barcodes not found in the hash are discarded. Use
`--cb-max-hamming 2` for chemistry mismatches (e.g. Nextera-formatted
barcodes against a TruSeq whitelist).

**Step 1c — UMI Encoding:** The UMI region (12 bp at R1[16:28] for `10xv3`;
10 bp at R1[16:26] for `10xv2-5p`; user-specified for `custom`) is encoded
as a `uint32`. Decoding is deferred to the deduplication stage.

**Step 1d — Window-Restricted Guide Extraction:** Only a specific window of
Read2 is examined, anchored against the known construct layout:

| Chemistry | R2 window | Guide length | Margin / offsets |
|:---|:---|:---:|:---|
| `10xv3` (3' v3) | positions 28–54 | 20 bp | ±3bp, 7 offsets |
| `10xv2-5p` / `10xv2-5p-12umi` (5') | positions 13–39 | 20 bp | ±3bp, 7 offsets |
| `custom` | user-specified | user-specified | `window_len - guide_len + 1` offsets |

Example for `10xv3`:
```
Read2 layout (3' v3 chemistry, sgRNA library):
  [0:10]   — 10x TSO adapter remnant (discarded)
  [10:30]  — poly-G / linker (discarded)
  [28:48]  — guide sequence (20 bp)        ← target
  [48:54]  — poly-A tail start (6 bp)      ← flanking context
  [54:90]  — poly-A tail remainder (discarded)
```

The guide anchor within each window (position 28 for `10xv3`, position 16
for the 5' chemistries) is the validated real-data position, not derived
from a generic construct diagram; the ±3bp margin around it is a
positional-drift tolerance, matching `10xv3`'s. If your own library prep
puts the guide somewhere else in R2 entirely, use `--chemistry custom`
with your own `--window-start`/`--window-end`/`--guide-len` rather than
assuming a named chemistry fits — see "Supported 10x Chemistries" above.

The window is encoded as a big integer (2 bits/base). `window_len - guide_len
+ 1` sliding sub-windows are extracted via bit-shift-and-mask operations —
7 for all three named chemistries, given their current ±3bp margins. Each
candidate is checked against `guide_hash_int`:
- **Fast path (exact match):** O(1) dict lookup. Handles >99% of guide matches.
- **Slow path (Hamming=1):** Generates 60 Hamming-1 variants on-the-fly for
  the ~1% of reads with a sequencing error (60 × O(1) = constant overhead).

**Step 1e — Multi-File Parallelism:** When multiple FASTQ pairs are provided
(comma-separated), HAM uses `multiprocessing.Pool` to process each pair in
parallel. The guide hash and CB hash are shared read-only via `fork()`
copy-on-write.

### Stage 2: `dedup` — UMI Deduplication + MEX Generation

**Input:** `hits.npz` from the match stage.

**Output:** Standard MEX trio (`matrix.mtx.gz`, `barcodes.tsv.gz`,
`features.tsv.gz`).

Hits are grouped by `(cell_barcode_idx, guide_idx)`. Encoded UMIs are decoded
to strings with the chemistry-appropriate UMI length. For each group, UMIs are
deduplicated with a **hash-accelerated greedy count-ranked dedup**: UMIs are
visited in descending count order, and a UMI is retained unless it (or an
exact match) was already claimed by a higher-count UMI within the Hamming
threshold — claiming is done by inserting each retained UMI's 36 Hamming-1
variants into a local hash set (O(1) lookup instead of O(n × m) linear scan).
This is *not* the UMI-tools "directional" algorithm (which merges based on a
count-ratio condition over a UMI adjacency graph) and can produce different
counts from it, particularly at high UMI diversity/depth — treat it as its
own well-defined method, not a drop-in equivalent of UMI-tools directional.

Deduplicated counts are assembled into a `scipy.sparse.csr_matrix` and written
as a gzip-compressed MEX file.

### Stage 3: `merge` — Cross-Lane Merging

Per-lane sparse matrices are vertically stacked via `scipy.sparse.vstack`.
Lane suffixes (e.g. `-L01`) are appended to barcodes for cross-lane
uniqueness. Features must be identical across lanes.

---

## Design Rationale — Why Hash Lookup Beats k-mer Pseudoalignment

### The Problem with k-mer Approaches

k-mer pseudoalignment tools (piscem, kallisto) map reads by decomposing them
into short k-mers and looking up each in a pre-built coloured de Bruijn graph
index. For short guide constructs on a 10x backbone, k-mers spanning the
poly-G / TSO boundary can accidentally match genuine guide k-mers, producing
**adapter-derived false positives**.

### HAM's Solution

1. **Window-restricted extraction eliminates adapter noise entirely.** The
   search is confined to the guide region — the adapter and poly-G linker are
   never examined for guide content.
2. **Integer encoding eliminates string allocation.** A 20 bp guide comparison
   is a single 64-bit dict lookup, not a string comparison.
3. **Pre-computed Hamming-1 expansion handles sequencing errors at zero
   runtime cost.** The fast path (exact match) catches >99% of reads with a
   single O(1) dict lookup. The slow path is only invoked for rare reads where
   the error falls outside the pre-computed set.

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
  --chemistry NAME    10x chemistry: 10xv3, 10xv2-5p, or custom [default: 10xv3]

  # Custom chemistry flags (used with --chemistry custom):
  --cb-start N        CB start position in R1 [default: 0]
  --cb-end N          CB end position in R1 [default: 16]
  --umi-start N       UMI start position in R1 [default: 16]
  --umi-end N         UMI end position in R1 [default: 28]
  --window-start N    Guide window start in R2 [default: 28]
  --window-end N      Guide window end in R2 [default: 54]
  --guide-len N       Guide protospacer length in bp [default: 20]

ham dedup           UMI deduplication + MEX matrix generation
  -i, --input PATH    Hits file (.npz) from ham match
  -o, --output-dir PATH  Output directory for MEX files
  -t, --umi-threshold N  UMI Hamming distance threshold [default: 1]
  --umi-len N         UMI length in bp (12 for v3, 10 for 5' v1) [default: 12]

ham merge           Merge per-lane count matrices
  --lanes PATH        Lane list TSV: lane_id, matrix_dir, suffix
  --out PATH          Output directory
  --prefix PREFIX     Output file prefix [default: merged]
```
