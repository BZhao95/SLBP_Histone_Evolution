#!/bin/bash

# ==============================================================================
# Script: 11_verify_histone_domains.sh
# Description: Confirms Histone-Fold Domains (HFD) using RPSTBLASTN against CDD.
#              Processes sequences individually to ensure top-hit accuracy.
# ==============================================================================

# --- Configuration ---
QUERY_FASTA="data/histones_demo.fasta" #directly download from NCBI. nucleotide CDS
BLAST_DB="/Users/zhaobin/Downloads/Cdd//Cdd" # Update to local path
OUT_DIR="results/histones"
ALL_OUT="$OUT_DIR/rpstblastn_all_results.txt"
TOP_OUT="$OUT_DIR/rpstblastn_top_hits.txt"
THREADS=8  

mkdir -p "$OUT_DIR"

# --- Execution ---
echo "[INFO] Starting RPSTBLASTN domain verification..."

# Extracting headers and looping for individual sequence processing
grep '^>' "$QUERY_FASTA" | while IFS= read -r header; do
    # Parsing Accessions
    accession=$(echo "$header" | cut -d '|' -f 2 | cut -d '_' -f 3)
    Seq_header=$(echo "$header" | cut -d ' ' -f 1 | cut -d '>' -f 2)
    
    echo "Processing: $accession"
    
    # Extract temporary sequence
    seqtk subseq "$QUERY_FASTA" <(echo "$Seq_header") > temp_query.fasta
    
    if [[ ! -s temp_query.fasta ]]; then
        echo "[WARN] No sequence data for $accession. Skipping."
        continue
    fi
    
    # RPSTBLASTN against CDD
    rpstblastn -query temp_query.fasta \
               -db "$BLAST_DB" \
               -evalue 1e-5 \
               -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
               -max_target_seqs 2 \
               -num_threads "$THREADS" \
               -out temp_result.txt

    # Recording results
    if [[ -s temp_result.txt ]]; then
        cat temp_result.txt >> "$ALL_OUT"
        head -n 1 temp_result.txt >> "$TOP_OUT"
    else
        echo -e "$accession\tNo HFD hit found" >> "$ALL_OUT"
        echo -e "$accession\tNo HFD hit found" >> "$TOP_OUT"
    fi
done

rm temp_query.fasta temp_result.txt
echo "[FINISH] Domain verification complete. Results in $OUT_DIR"
