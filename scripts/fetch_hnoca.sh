#!/usr/bin/env bash
# Download the HNOCA integrated atlas and supporting artifacts.
#
# TODO: fill in the canonical HNOCA download URLs / Zenodo DOIs before running.
# Cross-reference docs/data_sources.md and the HNOCA publication
# (He, Dony, et al., Nature 2024) for the exact records.
#
# Usage:
#   bash scripts/fetch_hnoca.sh [DEST_DIR]
# Default DEST_DIR: data/reference/hnoca

set -euo pipefail

DEST="${1:-data/reference/hnoca}"
mkdir -p "$DEST"

# --- HNOCA integrated AnnData ----------------------------------------------
# HNOCA_ATLAS_URL=""   # TODO: paste Zenodo / figshare URL for the integrated .h5ad
# HNOCA_ATLAS_OUT="$DEST/hnoca_integrated.h5ad"
# curl -L --fail -o "$HNOCA_ATLAS_OUT" "$HNOCA_ATLAS_URL"

# --- HNOCA reference embedding / model -------------------------------------
# HNOCA_MODEL_URL=""   # TODO: scVI / scANVI checkpoint, if published
# curl -L --fail -o "$DEST/hnoca_model.tar.gz" "$HNOCA_MODEL_URL"
# tar -xzf "$DEST/hnoca_model.tar.gz" -C "$DEST"

# --- Primary brain reference (developing) ----------------------------------
# PRIMARY_REF_URL=""   # TODO: pick reference (e.g. Braun 2022, Kanton 2019)
# curl -L --fail -o "$DEST/primary_dev_brain.h5ad" "$PRIMARY_REF_URL"

echo "fetch_hnoca.sh: stub — populate URLs in this script before running." >&2
exit 1
