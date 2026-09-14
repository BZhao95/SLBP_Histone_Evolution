#!/bin/bash
set -euo pipefail

# ==============================================================================
# Script: 05_map_taxonomy.sh
# Description: Maps verified SLBP accessions back to their NCBI Taxonomy IDs.
#              Prepares the dataset for clade-level evolutionary analysis.
# ==============================================================================

# --- Configuration ---
# This should be the list of hits that PASSED the RPS-BLAST check (from Script 04)
VERIFIED_ACC="results/rpsblast/slbp_hits.txt"
# The original raw BLAST output containing TaxID info (from Script 01)
RAW_BLAST="results/blastp_results.txt"

# Output Files
OUT_DIR="results/taxonomy"
TAXID_LIST="$OUT_DIR/slbp_unique_taxids.txt"
TEMP_FILTERED="$OUT_DIR/tmp_mapped_hits.tsv"

mkdir -p "$OUT_DIR"

if [[ ! -f "$VERIFIED_ACC" ]]; then
    echo "[ERROR] Verified accession file not found: $VERIFIED_ACC"
    exit 1
fi

if [[ ! -f "$RAW_BLAST" ]]; then
    echo "[ERROR] Raw BLAST results not found: $RAW_BLAST"
    exit 1
fi

echo "[INFO] Mapping verified accessions to Taxonomy IDs..."

# use fgrep (fast grep) for fixed strings since we have a large list of patterns.
# Column 12 in BLAST format is 'staxids'.
grep -F -f "$VERIFIED_ACC" "$RAW_BLAST" > "$TEMP_FILTERED"

# Extract TaxIDs, handle semi-colon separated lists (multispecies), and get unique IDs.
echo "[INFO] Extracting unique TaxIDs..."
awk '{print $13}' "$TEMP_FILTERED" | tr ';' '\n' | sort -u > "$TAXID_LIST"

TOTAL_ACC=$(wc -l < "$VERIFIED_ACC")
UNIQUE_SPECIES=$(wc -l < "$TAXID_LIST")

echo "---------------------------------------"
echo "Taxonomy Mapping Complete"
echo "---------------------------------------"
echo "Verified Accessions: $TOTAL_ACC"
echo "Unique Taxonomy IDs: $UNIQUE_SPECIES"
echo "---------------------------------------"
echo "[FINISH] TaxID list saved in $TAXID_LIST"

# Cleanup
rm "$TEMP_FILTERED"
