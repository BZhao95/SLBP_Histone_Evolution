#!/bin/bash
set -euo pipefail

# ==============================================================================
# Script: 02_fetch_sequences.sh
# Description: Extracts unique accessions from BLAST results and retrieves 
#              FASTA sequences from NCBI using Entrez Direct (EDirect).
# ==============================================================================

# --- Configuration ---
BLAST_RESULTS="results/blastp_results.txt"
ACCESSION_FILE="results/blastp_accessions.txt"
CLEANED_ACCESSIONS="results/cleaned_accessions.txt"
SEQUENCES_FILE="results/retrieved_protein_sequences.fasta"
FAILED_BATCHES="logs/failed_batches.txt"
BATCH_SIZE=90 #NCBI limit 100

# --- pre-check ---
if [[ ! -f "$BLAST_RESULTS" ]]; then
    echo "[ERROR] BLAST results file '$BLAST_tRESULTS' not found!"
    exit 1
fi

# --- Step 1: Accession Extraction & Cleaning ---
echo "[INFO] Extracting and cleaning subject accessions..."
awk '{print $2}' "$BLAST_RESULTS" > "$ACCESSION_FILE"

# Removes prefixes like dbj| or emb| to get clean NCBI accessions
sed -E 's/^[^|]+\|([^|]+)\|?/\1/' "$ACCESSION_FILE" | sort -u > "$CLEANED_ACCESSIONS"

if [[ ! -s "$CLEANED_ACCESSIONS" ]]; then
    echo "[ERROR] No valid accessions found."
    exit 1
fi

# --- Step 2: Sequence Retrieval ---
echo "[INFO] Fetching sequences from NCBI in batches of $BATCH_SIZE..."
rm -f "$SEQUENCES_FILE" "$FAILED_BATCHES"

# Split accessions into temp batch files
split -l $BATCH_SIZE "$CLEANED_ACCESSIONS" accession_batch_

for BATCH in accession_batch_*; do
    echo "Processing $BATCH..."
    
    # Search protein DB and fetch in FASTA format
    if ! esearch -db protein -query "$(paste -sd ' ' "$BATCH")" | efetch -format fasta >> "$SEQUENCES_FILE"; then
        echo "$BATCH" >> "$FAILED_BATCHES"
        echo "[WARNING] Batch $BATCH failed. Logged for retry."
    fi
    
    sleep 2  # avoid violating NCBI rules
done

# Cleanup temporary batch files
rm -f accession_batch_*

echo "[FINISH] Retrieval complete. Sequences saved to $SEQUENCES_FILE."
