#!/bin/bash
set -euo pipefail

# ==============================================================================
# Script: 04_rpsblast_verification.sh
# Description: check the presence of the SLBP RNA-Binding Domain (RBD)
#              using RPS-BLAST against the NCBI CDD.
# Target Domain: Histone RNA hairpin-binding protein (PSSM: gnl|CDD|464587)
# ==============================================================================

# --- Configuration ---
# Update CDD_DB path to your local database location
CDD_DB="/Users/zhaobin/Downloads/Cdd/Cdd" #the whole CDD, downloaded from NCBI
SEQUENCES="results/retrieved_protein_sequences.fasta"

# Output Files
OUT_DIR="results/rpsblast" #all rpsblast results
RAW_OUT="$OUT_DIR/rpsblast_results_evalue1.txt" #retained the hits with evalue > 1
TOP_HITS="$OUT_DIR/rpsblast_top_hits.txt" #take the top1 hit
SLBP_HITS="$OUT_DIR/slbp_hits.txt"
SLBP_MISSED="$OUT_DIR/slbp_missed.txt"
MISSING_QUERIES="$OUT_DIR/missing_from_blast.txt"

# Temporary Tracking Files
ALL_IDS="$OUT_DIR/tmp_all_ids.txt"
FOUND_IDS="$OUT_DIR/tmp_found_ids.txt"

mkdir -p "$OUT_DIR"

# --- 1. Pre Checks ---

if [[ ! -f "${CDD_DB}.pal" ]]; then
    echo "[ERROR] CDD database not found at $CDD_DB"
    exit 1
fi

# --- 2. RPS-BLAST Execution ---
if [[ ! -f "$RAW_OUT" ]]; then
    echo "[INFO] Running RPS-BLAST (E-value threshold: 1)..."
    rpsblast -query "$SEQUENCES" \
             -db "$CDD_DB" \
             -out "$RAW_OUT" \
             -evalue 1 \
             -outfmt "6 qseqid sseqid evalue bitscore" \
             -num_threads 8
else
    echo "[INFO] Existing RPS-BLAST results found. Skipping search."
fi

# --- 3. Hit Filtering & Selection ---
echo "[INFO] Processing results to identify top hits..."

# Sort by query ID then by E-value (lowest first), then take the first instance of each query
awk '!seen[$1]++' <(sort -k1,1 -k3,3g "$RAW_OUT") > "$TOP_HITS"

# Extract all IDs from the input FASTA for comparison
grep "^>" "$SEQUENCES" | sed 's/^>//' | awk '{print $1}' | sort -u > "$ALL_IDS"

# --- 4. Reconciliation ---
# Identify which sequences returned NO results in RPS-BLAST
awk '{print $1}' "$TOP_HITS" | sort -u > "$FOUND_IDS"
comm -23 "$ALL_IDS" "$FOUND_IDS" > "$MISSING_QUERIES"

# Identify the specific SLBP RBD (PSSM 464587)
awk '$2 == "gnl|CDD|464587" {print $1}' "$TOP_HITS" | sort -u > "$SLBP_HITS"

# Identify sequences that had hits, but NOT the SLBP domain
comm -23 "$ALL_IDS" "$SLBP_HITS" > "$SLBP_MISSED"

# --- 5. Summary Report ---
EXPECTED=$(wc -l < "$ALL_IDS")
SLBP_COUNT=$(wc -l < "$SLBP_HITS")
MISS_COUNT=$(wc -l < "$SLBP_MISSED")
NO_BLAST=$(wc -l < "$MISSING_QUERIES")

echo "---------------------------------------"
echo "RPS-BLAST Verification Summary"
echo "---------------------------------------"
echo "Total Queries:    $EXPECTED"
echo "SLBP RBD Hits:    $SLBP_COUNT"
echo "Non-SLBP Hits:    $MISS_COUNT"
echo "Zero Hits Found:  $NO_BLAST"
echo "---------------------------------------"
echo "[FINISH] Results saved in $OUT_DIR"

# Cleanup temporary files
rm "$ALL_IDS" "$FOUND_IDS"
