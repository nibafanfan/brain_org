# Snapseed cell-type annotation

The HNOCA marker hierarchy lives at:

```
data/reference/hnoca_snapseed_markers.yaml
```

Source: `theislab/neural_organoid_atlas/supplemental_files/Data_S1_snapseed_markers.yaml`. This is the same YAML used to annotate the published HNOCA atlas — using it on our deposits gives **directly comparable cell-type labels** with HNOCA's `annot_level_1..4_rev2`.

## What's in it

- **55 total nodes**, **38 leaf cell types**, **max depth 5**.
- **10 top-level classes:** `neural_progenitor_cell`, `neuron`, `choroid_plexus_epithelium`, `astrocyte`, `oligodendrocyte_lineage`, `microglia`, `vascular_endothelial_cell`, `mesenchymal_cell`, `neural_crest`, `pns_neurons`.
- Deep neuronal subtyping under `neuron` (telencephalic dorsal/ventral, hippocampal, diencephalic, mesencephalic, hindbrain, etc.) and progenitor subtyping under `neural_progenitor_cell`.

## How Snapseed uses it

[Snapseed](https://github.com/devsystemslab/snapseed) is a marker-based annotation tool that walks the hierarchy:

1. At each level, score every cell against each candidate's marker set (`UCell`-style aggregate expression).
2. Assign the cell to the highest-scoring candidate above a confidence threshold.
3. Drop into that candidate's `subtypes` and repeat until reaching a leaf or falling below threshold.

Output: one annotation per cell at each level (`annot_level_1`, `annot_level_2`, ...). HNOCA's atlas uses level 1–4 plus separate `annot_region_rev2` (brain region) and `annot_ntt_rev2` (neurotransmitter type) annotations.

## Running it on our atlas

After the rebuild + integration:

```python
import snapseed as ss
import yaml, anndata as ad

# Load the integrated atlas
adata = ad.read_h5ad('data/processed/atlas_integrated.h5ad')  # or per-deposit

# Load hierarchy
with open('data/reference/hnoca_snapseed_markers.yaml') as f:
    markers = yaml.safe_load(f)

# Run hierarchical annotation
ss.annotate_hierarchy(
    adata,
    marker_hierarchy=markers,
    group_name='leiden',           # or any pre-computed cluster column
    layer='counts',                 # raw counts layer
)
# Results land in adata.obs as 'level_1', 'level_2', ... per cluster
```

`pip install snapseed` (or `pip install git+https://github.com/devsystemslab/snapseed`). Requires `scanpy`, `numpy`, `pandas` (already in project's `pyproject.toml`).

## Multi-lineage coverage

The hierarchy directly covers all the **added cell types** that drive Q1 (multi-lineage benchmark):

| Lineage | Node | Top marker(s) |
|---|---|---|
| Vascular endothelial | `vascular_endothelial_cell` | CLDN5 |
| Microglia | `microglia` | AIF1 |
| Astrocyte | `astrocyte` | GFAP, AQP4 |
| Oligodendrocyte (OPC + mature) | `oligodendrocyte_lineage` | OLIG1, MBP |
| Choroid plexus | `choroid_plexus_epithelium` | TTR |
| Mesenchymal | `mesenchymal_cell` | DCN |
| Neural crest | `neural_crest` | SOX10 |
| PNS neurons | `pns_neurons` | PRPH |

## Gaps to extend

These cell types appear in our 2024–2026 multi-lineage deposits but **aren't in the upstream YAML** — add them locally if needed:

- **Pericyte** — Wang neural-perivascular assembloid (`gse224346`). Suggested markers: `PDGFRB`, `RGS5`, `MCAM`, `NOTCH3`. Place as `pericyte` under root (alongside `vascular_endothelial_cell`) or as a subtype of `mesenchymal_cell`.
- **Vascular smooth muscle cell** — some vascularized protocols. Markers: `ACTA2`, `MYH11`, `TAGLN`.
- **Microglia subtypes** — if doing fine-grained microglial states (M1/M2, DAM, HAM): `TMEM119` (homeostatic), `P2RY12` (homeostatic), `CD68` / `CD86` (activated), `TREM2` (DAM).
- **Brain-region-specific neuronal types not in HNOCA** — e.g., dopaminergic midbrain (`TH`, `NR4A2`, `LMX1A`), cerebellar Purkinje (`CALB1`, `PCP2`).

Extensions should go into a separate sidecar YAML (`data/reference/extension_markers.yaml`) and be merged at annotation time, so the canonical HNOCA file stays a pure mirror of upstream.

## Alignment with HNOCA's published labels

HNOCA's atlas already carries the result of running this hierarchy: `obs['annot_level_1'..'annot_level_4_rev2']`. The `_rev2` suffix indicates they re-annotated once. When we project our deposits onto HNOCA's scPoli embedding (via reference mapping), label transfer can use these directly — Snapseed only needs to run on cells that don't get a confident HNOCA-derived label.

## When to use Snapseed vs scANVI vs label transfer

| Method | Best for | When to choose |
|---|---|---|
| **Snapseed** | Hierarchical, interpretable, marker-based | Cells with no nearby HNOCA neighbor; new lineages not in HNOCA; sanity-checking a transfer |
| **scANVI label transfer** | Soft, scvi-tools native | Cells projected into a shared latent space with HNOCA |
| **scPoli built-in label transfer** | Direct from reference mapping | Default for cells with a confident reference assignment |

Recommended workflow:
1. Project our deposits onto HNOCA's scPoli embedding → get reference-based labels for the ~80% of cells that map cleanly.
2. Run Snapseed on the remaining ~20% (low-confidence reference neighbors) using the YAML above.
3. Sanity-check by comparing Snapseed level 1 vs scPoli-transferred annot_level_1 — should agree on the easy cases.
