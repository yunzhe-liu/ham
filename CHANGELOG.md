# Changelog

## 1.0.0 — 2026-08-28

### Fixed
- **Critical**: `extract_guides_from_bigint`'s sliding-window shift formula
  for `guide_len==20` (the default `10xv3` chemistry) was broken by a prior
  refactor (commit `00913be`) — every one of the 7 candidate offsets
  extracted a truncated, zero-padded bit fragment instead of a real 20bp
  guide substring, so exact/Hamming=1 matching against a real guide hash
  produced effectively zero hits. Restored the original, verified-correct
  formula and generalized it to arbitrary window/guide lengths so future
  chemistry additions can't reintroduce the same bug. Guarded by regression
  tests in `tests/test_encoding.py` and `tests/test_matcher.py`.
  Production `10xv3` results generated *before* this regression was
  introduced are unaffected.

### Changed
- `10xv2-5p` / `10xv2-5p-12umi` guide length corrected from 19bp to 20bp
  (standard SpCas9 protospacer length). The guide anchor position within
  R2 is unchanged (still R2 position 16 — validated against real 5' v1
  data); the window was widened from `[16:35]` (19bp, zero margin) to
  `[13:39]` (26bp, ±3bp margin on each side of the 20bp guide), mirroring
  `10xv3`'s ±3bp tolerance and giving 7 candidate offsets instead of 1. An
  earlier draft of this release moved the anchor itself to R2 position 0
  based on a generic Direct Capture construct diagram; that was reverted
  before release — the anchor position (16) was already correct and
  validated, only the guide length and margin needed correcting.
- `dedup_umis_directional`'s docstring corrected to accurately describe the
  implemented algorithm (a hash-accelerated, count-ranked greedy Hamming-1
  collapse) — it is not equivalent to UMI-tools' directional algorithm
  (which merges via a count-ratio condition over a UMI adjacency graph).
  The implementation itself is unchanged.

### Added
- `tests/` rebuilt with a real pytest suite (encoding round-trips, hash
  builder, dedup, chemistry-config validation, end-to-end synthetic FASTQ
  matching for both bugs above).
- `.github/workflows/ci.yml`: pytest across Python 3.10–3.12 on push/PR.
- `chem_cfg` validation (`_validate_chem_cfg`) shared by named and `custom`
  chemistries — catches inverted/degenerate positions and a window shorter
  than the guide length before matching starts.
- `guide_len` vs. `guide_hash.pkl`'s recorded `guide_length` consistency
  warning.
- FASTQ input validation: missing file(s) raise `FileNotFoundError`; R1/R2
  record-count mismatch or a non-`@` header raises a descriptive
  `ValueError` instead of silently producing misaligned or truncated hits.
- `hits` array in the match stage now grows by amortized doubling instead
  of a hardcoded 30,000,000-row cap.
- Named-chemistry field override: `chem_cfg` (CLI: `--cb-start`/`--cb-end`/
  `--umi-start`/`--umi-end`/`--window-start`/`--window-end`/`--guide-len`)
  can now be passed alongside a *named* chemistry (not just `custom`) to
  override individual fields on top of that chemistry's preset — e.g.
  `--chemistry 10xv2-5p --guide-len 19` for a guide library whose actual
  protospacer length differs from the 20bp default without having to
  respecify every other position via `custom`. Previously these flags were
  silently ignored unless `--chemistry custom` was set. `custom` itself is
  unchanged: it still requires all 7 keys, since there's no preset to fall
  back on.

### Internal
- `multiprocessing.set_start_method('fork', force=True)` (a process-wide
  side effect) replaced with `multiprocessing.get_context('fork').Pool(...)`
  so calling `match_reads` as a library no longer overrides the caller's
  global multiprocessing start method.

### Not changed
- `10xv3`'s window (`R2[28:54]`, guide_len=20, anchor position 28) is left
  exactly as originally shipped — validated against real 3' v3 production
  data. README now documents all three named-chemistry anchor positions as
  validated-not-generic, and points users with a different construct at
  `--chemistry custom` instead of expecting a built-in margin to cover it.
