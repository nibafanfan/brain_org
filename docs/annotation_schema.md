# Annotation Workbook Schema (metadata reference for code review)

Describes the structure of the organoid annotation/metadata workbooks **without
publishing the raw data** (the `*.xlsx` files are gitignored by policy). This is
enough for a reviewer (e.g. Codex) to reason about the metadata-analysis logic;
the actual parsing/injection code is in the tracked scripts listed at the bottom.

## Files & relationship

| File | Role |
|---|---|
| `data/brain_organoid (5).xlsx` | working iteration (has duplicate `… (1)` sheets = in-progress revisions) |
| `data/brain_organoid_GSMannotations.xlsx` | **finalized source of truth** (per-GSM annotations, cleaned encodings) |

Both share the same sheet structure. Only `*_GSMannotations.xlsx` should be treated
as authoritative for per-cell injection. The atlas was annotated from it (memory:
"finalized GSMannotations.xlsx is source of truth").

## Sheets (same in both files)

- **HNOCA Paper Datasets** (27 rows) — the papers that overlap the published HNOCA
  atlas (deposit-level, mostly used to mark HNOCA provenance / exclusion).
- **Paper Datasets** (~176 rows) — candidate-paper-level metadata (one row per paper/
  deposit), free-text, pre-GSM-expansion. Includes `Search Modality` (GEO/GPT/Self),
  `Dataset` (where data lives), `Dataset Notes` (e.g. `FASTQ Only`, `NOT PRIMARY DATA`,
  `NOT SINGLE CELL`), and the `Samples` cell that gets parsed into GSMs.
- **GSM Annotations** — **the key sheet**: one row per GSM (sample), the source for
  per-cell annotation injection. ~696–746 rows depending on file/iteration.
- **Abstract / Summary / Summary (1)** — project abstract + run-summary stats
  (e.g. "Total GSMs emitted", "Multi-Lineage? = Yes" counts); not data, just QC tallies.

## `GSM Annotations` sheet — columns (the injectable schema)

| Column | dtype | Cardinality / vocabulary | Notes |
|---|---|---|---|
| `Paper` | str | ~144 | PubMed/PMC title |
| `GEO` | str | ~146 | GSE accession (deposit id) |
| `GSM` | str | ~659–720 | GEO sample accession — **the join key to cells** |
| `GSM Label` | str | ~678 | free-text sample label (e.g. "WNTi (batch A)") |
| `Barcode ID` | (empty) | 0 | placeholder, unused |
| `Cell Type Origin` | str | `esc`, `ipsc`, `ipsc, esc`, `tissue` | → obs `cell_type_origin` |
| `Age (Days in Vitro)` | num/str | ~38–82 | days in vitro; **free-text in some sheets** ("69", "Day 100\nDay 130") → obs `age_days` |
| `Organoid Type` | str | ~39–45 | e.g. Telencephalic, Midbrain, Cortical… (trailing whitespace exists) → obs `organoid_type` |
| `Protocol` | str | ~55–58 | e.g. Conventional, Engineered Vascular Scaffold → obs `protocol` |
| `Unguided?` | bool-ish | `0/1` **or** `No/Yes` | → obs `unguided` |
| `Multi-Lineage?` | bool-ish | `0/1`, `No/Yes`, **and `"1 (sort of)"`/`"Yes (sort of)"`** | → obs `multi_lineage` |
| `Vascularized?` | bool-ish | `0/1` or `No/Yes` | → obs `vascularized` |
| `Slice?` | bool-ish | `0/1` or `No/Yes` | → obs `slice` |
| `Notes` | str | sparse | e.g. "All samples stored within the same matrix. Barcodes contain individual label…" (flags pooled/hashed deposits) |

### → atlas `obs` field map (8 finalized fields + provenance)

`GSM`→`gsm`, `Cell Type Origin`→`cell_type_origin`, `Age (Days in Vitro)`→`age_days`,
`Organoid Type`→`organoid_type`, `Protocol`→`protocol`, `Unguided?`→`unguided`,
`Multi-Lineage?`→`multi_lineage`, `Vascularized?`→`vascularized`, `Slice?`→`slice`,
plus **`annotation_level`** (provenance flag, below). Injection is keyed by `gsm`
(see `inject_finalized_annotations.py`).

## `annotation_level` semantics (gsm vs deposit)

Per-cell flag recording **how** a cell got its annotation:
- **`gsm`** — the deposit's `Samples` cell listed explicit per-GSM rows, so the cell's
  annotation came from its own GSM row ("explicit GSM path").
- **`deposit`** — the `Samples` cell was qualified "all samples" (no per-GSM split), so
  all cells in the deposit inherited the deposit-level annotation ("all samples path").

Downstream: `gsm`-level cells are higher-confidence; benchmarking can stratify by
`annotation_level` to check that conclusions aren't driven by deposit-level fallback.

## `Samples` cell grammar (what the parser handles)

Free-text, one deposit per cell. Observed forms:
- explicit list: `GSM8501109 - WNTi (batch A)\nGSM8501110 - WNTi (batch B)` (tab/space/`-`
  separated GSM + label, newline per sample) → **explicit GSM path** → `annotation_level=gsm`
- range: `GSM2665699 to GSM2666406` (expanded to the full GSM range)
- qualified bulk: `All samples`, or `All samples\nWT and Mutant are mixed…` → **all-samples
  path** → `annotation_level=deposit`
- skip conditions (logged in Summary): empty `Samples` cell, no GEO accession, descriptive-only
  (no GSMs), FASTQ-only deposits (deleted).

## Known data-quality quirks (important for any analysis)

1. **Boolean encoding is inconsistent across sheets/iterations**: `Unguided?/Vascularized?/
   Slice?/Multi-Lineage?` appear as `0/1` ints in some sheets and `No/Yes` strings in others.
   This is the root of the atlas `multi_lineage` mixed-dtype column (`'0'/'1'/'False'/'True'/'No'`)
   — normalize once: `{'1','True','Yes'}→multi`, `{'0','False','No'}→single`.
2. **`"… (sort of)"` qualifiers** (`Multi-Lineage? = "Yes (sort of)"`, `Guided (sort of)`) —
   decide a rule (we treat `"Yes (sort of)"` as not-multi unless stated).
3. **Trailing whitespace** in categoricals (`"Telencephalic "`) — strip before grouping.
4. **`Age` is free-text** at deposit level (`"Day 100\nDay 130\nDay 175"`) and sometimes at GSM
   level (`"69 Days"`) — needs unit-anchored regex (Summary notes a fix tightening this).
5. **Duplicate `… (1)` sheets** in `(5).xlsx` are older revisions — ignore; use
   `*_GSMannotations.xlsx`.
6. **Multi-value deposit cells** (newline-separated `Organoid Type`/`Protocol`) at the Paper/
   HNOCA level collapse to a single value per GSM at the GSM level.

## Tracked parsing/injection scripts (already in repo)

- `scripts/inject_finalized_annotations.py` — keys `GSM Annotations` rows to cells by `gsm`,
  writes the 8 obs fields + `annotation_level`.
- `scripts/recover_per_cell_annotations.py` — per-cell annotation recovery.
- `scripts/patch_obs_annotations.py` — in-place obs patching of existing deposits.
- `scripts/migrate_to_hnoca_schema.py` — projects obs onto the canonical HNOCA schema.

**Do not commit the `.xlsx` files** (gitignored by policy). To regenerate this schema from
the local workbooks, introspect columns/vocabularies with pandas+openpyxl (no raw rows).
