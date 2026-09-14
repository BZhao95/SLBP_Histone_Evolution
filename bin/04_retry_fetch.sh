#!/bin/bash
set -euo pipefail

# ==============================================================================
# Script: 03_retry_fetch.sh
# Description: Retries retrieval for accessions that failed in Script 02.
# ==============================================================================

# --- Configuration ---
FAILED_ACCESSIONS="logs/failed_batches.txt"
SEQUENCES_FILE="data/not_retrieved_protein_sequences.fasta"
FINAL_FAILED="logs/final_failed_accessions.txt"
BATCH_SIZE=10 #use a smaller batch size to avoid ncbi killing the process

if [[ ! -f "$FAILED_ACCESSIONS" || ! -s "$FAILED_ACCESSIONS" ]]; then
    echo "[INFO] No failed accessions found to retry."
    exit 0
fi

echo "[INFO] Retrying failed accessions..."
rm -f "$FINAL_FAILED"

# Split using 5-character suffix for safety
split -l $BATCH_SIZE -a 5 "$FAILED_ACCESSIONS" retry_batch_

for BATCH in retry_batch_*; do
    echo "Retrying $BATCH..."

    # Logic: Direct ID fetch using comma-separated list
    if ! efetch -db protein -format fasta -id "$(paste -sd ',' "$BATCH")" >> "$SEQUENCES_FILE"; then
        cat "$BATCH" >> "$FINAL_FAILED"
        echo "[WARNING] Accessions in $BATCH failed again."
    fi

    sleep 2
done

rm -f retry_batch_*

if [[ -s "$FINAL_FAILED" ]]; then
    echo "[ALERT] Some accessions still failed. See $FINAL_FAILED."
else
    echo "[SUCCESS] All failed accessions retrieved."
fi
