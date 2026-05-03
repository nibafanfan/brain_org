# brain_organoid

Benchmark of multi-lineage human brain organoid protocols, extending the
Human Neural Organoid Cell Atlas (HNOCA) framework to multi-lineage
differentiation methods published in 2024–2025.

See [`docs/proposal.md`](docs/proposal.md) for the full project proposal.

## Layout

```
configs/                   # YAML/JSON configs for integration & benchmarking runs
data/
  raw/                     # untouched downloads (gitignored)
  processed/               # harmonized AnnData (gitignored)
  external/                # third-party data (gitignored)
  reference/               # primary brain reference atlases (gitignored)
docs/                      # proposal, design notes
notebooks/                 # exploratory analyses
scripts/                   # CLI entrypoints, data fetchers
src/brain_organoid/
  integration/             # scVI / scANVI / SAE training & embeddings
  benchmarks/              # mapping accuracy, alignment, coverage metrics
  models/                  # model definitions
  utils/                   # I/O, QC, plotting helpers
tests/
results/                   # output tables, embeddings (gitignored)
figures/                   # plots (gitignored)
```

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

## Data

Primary brain references and HNOCA artifacts are downloaded by scripts under
`scripts/`. URLs and DOIs are listed in [`docs/data_sources.md`](docs/data_sources.md).
