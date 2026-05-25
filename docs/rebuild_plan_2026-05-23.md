# Atlas rebuild plan v2 — post-annotation-review

**Date:** 2026-05-23
**Supersedes (for the *next* full rebuild):** `docs/rebuild_plan_2026-05-17.md`
**Reads with:** `docs/handoff_2026-05-20.md` (current state), `docs/loader_format_reference.md`

This plan governs the **next full atlas rebuild**, which is gated on the colleague
finishing the annotation review (see "Gating" below). It keeps the entire
mechanism of the 2026-05-17 plan (source-of-truth TSV, uniform QC, HNOCA gene
projection, per-deposit `data/processed/<slug>.h5ad` + manifest) and only records
what is **new or changed** since then. Do not re-read the old plan for the
unchanged parts — they still hold.

## Gating: do not rebuild until these are resolved

The rebuild should run once, cleanly, after annotation review is complete. Before
kicking it off, confirm:

1. **The 27 RAW_MISSING_GSMS deposits** (`data/raw_vs_requested_audit.tsv`,
   `category=RAW_MISSING_GSMS`) are each either corrected in
   `data/brain_organoid (2).xlsx` or confirmed for GEO re-download. See
   handoff §A.
2. **SKIP rows** in `data/rebuild_config.tsv` have been re-decided where the
   colleague has filled in the Samples cell.
3. **The pooled/HTO deposits are flagged** (see next section) so they get the
   right loader instead of the pooled-MTX loader.

## NEW — pooled + cell-hashing (HTO) deposits: genotype lives only in the Seurat object

**Lesson from gse297594 (MeCP2 / Rett):** some deposits multiplex
isogenic-control (WT) and mutant cells into shared sequencing pools via cell
hashing. The **raw GEO MTX is per-pool and cannot be split by genotype** — the
WT/MUT call exists *only* in the processed Seurat object's `meta.data`
(`condition`, `hash.ID`). The 2026-05-17 build loaded the pooled MTX
(`10x_mtx_per_gsm`) and marked all 164,189 cells `is_control=True`, silently
including thousands of Rett-mutant cells in a control-only atlas.

**The fix / reusable pattern** — `scripts/fix_gse297594_control.py`:
1. R reads each `.rds`, subsets `condition=="WT"` Singlets, writes a 10x trio +
   per-cell metadata (`orig.ident`, `genotype`, `library`, `final_annotation`).
2. Python runs the **same** `apply_qc` + `map_to_hgnc` + `attach_canonical_var`
   + `migrate_obs` + `sanitize_for_h5` as `rebuild_atlas.py`, then overrides the
   sample hierarchy (below).

Result: gse297594 = **47,844 true WT cells** (was 164,189), 8 organoids, 0 mutant.

**Action for the rebuild:** add a `seurat_hto_control` loader_type (generalize
`fix_gse297594_control.py`) and route any pooled+hashed deposit to it via
`rebuild_config.tsv`. **Do NOT** leave these on `10x_mtx_per_gsm`.

**Audit before rebuild:** sweep the other deposits for the same trap — any deposit
that (a) ships a Seurat/processed object alongside pooled raw MTX, and (b) mixes
genotypes/conditions via hashing. Candidates worth checking first: anything in
`data/raw/*/*.rds` with a `condition`/`genotype`/`hash.ID` meta column
(gse325956 NRXN1-KO, gse271118 DJ1-KO, gse219245, gse277968, …).

## NEW — `organoid_type` is PER-CELL, and one deposit can hold several

gse297594 contains two sample types in one `dataset_slug`: cerebral organoids
(**CO**, cortical, single-lineage) and telencephalic assembloids (**TA**,
cortical+MGE, multi-lineage). These are stored **per cell**, not per deposit:

| obs column | meaning | gse297594 values |
|---|---|---|
| `organoid_type` | per-cell organoid/region identity | `Cortical organoid` / `Cortical-MGE telencephalic assembloid` |
| `source` | which source object/sample-type | `CO` / `TA` |
| `multi_lineage` | per-cell single vs multi | `False` (CO) / `True` (TA) |
| `bio_sample` | **per-organoid** biological replicate | `R133C-WT-CO1` … `R255X-WT-TA3` (8) |
| `tech_sample` | the multiplexed pool/library (real batch) | `gse297594_mecp2_lib1..4` (4) |

This is correct per-cell and flows through the concat untouched:
`concatenate_atlas.py:project_deposit` does `obs=a.obs.copy()` and
`concat_on_disk(axis=0)` stacks obs row-wise. **The manifest's single
`organoid_type` string is never written onto cells** — it is path/iteration
metadata only. So CO/TA labels are assigned correctly to every cell.

**Schema convention going forward (do this in the rebuild):**
- `organoid_type` and `multi_lineage` are **per-cell** fields. A deposit MAY emit
  more than one value of each. The manifest row carries only a representative/
  display label.
- Always set `source` (the sample-type / source-object tag) and a true per-organoid
  `bio_sample` when a deposit mixes sample types.
- **Downstream tallies must group on `bio_sample`/`source`, NOT `dataset_slug`.**
  In particular the "multi-lineage vs single-lineage" benchmark count (handoff
  §6 Q1) must NOT assume one organoid_type per slug — gse297594 straddles both.
- scVI/scANVI `batch_key='tech_sample'` is unaffected (it never used
  `dataset_slug`).

## Downstream (unchanged from handoff, after per-deposit rebuild)

Concat → HVG (`cell_ranger` flavor, ~3000, `batch_key='bio_sample'`) → scVI
(`batch_key='tech_sample'`, n_latent=30, ~15 epochs) → Braun 2023 label transfer
→ benchmark questions. See `docs/handoff_2026-05-20.md` §"Next steps" for the
specifics; nothing there changes except the cell counts below.

## NEW — exclude datasets already in HNOCA (no duplication)

This atlas **extends** HNOCA; datasets HNOCA already integrated must NOT be
re-included. The finalized annotation workbook
(`data/brain_organoid_GSMannotations.xlsx`) separates these: sheet
**`HNOCA Paper Datasets`** (the ~27 HNOCA datasets, excluded) vs **`Paper Datasets`**
/ **`GSM Annotations`** (the new 2024–2025 set, included). The per-GSM source of
truth is the **`GSM Annotations`** sheet (683 control GSMs / 142 accessions).

The membership diff (`data/membership_diff_2026-05-24.tsv`) flags any built
deposit absent from the finalized set as `NOT_IN_FINALIZED`. Only two appeared:
- **`gse168323`** — confirmed (manually) to be a HNOCA dataset → **DROP** from the
  atlas/manifest in the rebuild. It is the *only* HNOCA-overlap that got built.
- `gse297594_mecp2` — expected (control defined by HTO demux, not GSM); keep.

## Bookkeeping the rebuild must reconcile

- **`data/manifest.tsv` is currently stale for gse297594_mecp2** — it still lists
  164,189 cells / 4 samples; the on-disk h5ad is now 47,844 / 8 organoids. The
  full rebuild will regenerate the manifest; until then treat the file on disk as
  truth, the manifest row as wrong.
- **`data/rebuild_config.tsv` still points gse297594_mecp2 at `10x_mtx_per_gsm`**
  on the pooled MTX — must be switched to the `seurat_hto_control` route or a full
  rebuild re-introduces the contamination.
- **Atlas total** drops by ~116k once this propagates: 4,255,476 → ~4,139,131.
  The current concatenated atlas, `data/scvi_latent_full.h5ad`, and the
  integration eval were all built on the contaminated version and need a rebuild
  to reflect the fix.

## Pointers (delta from 2026-05-17 plan)

| What | Where |
|---|---|
| Pooled/HTO control-extraction template | `scripts/fix_gse297594_control.py` |
| Current state snapshot | `docs/handoff_2026-05-20.md` |
| GSM-mismatch review queue | `data/raw_vs_requested_audit.tsv` (`category=RAW_MISSING_GSMS`) |
| Everything else (QC, gene projection, loaders, source-of-truth TSV) | `docs/rebuild_plan_2026-05-17.md` |
