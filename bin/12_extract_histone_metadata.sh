#!/bin/bash

# ==============================================================================
# Script: 12_extract_histone_metadata.sh
# Description: Parses NCBI CDS headers to extract protein/nucleotide IDs,
#              genomic locations, and partial sequence status.
# ==============================================================================

QUERY_FASTA="data/histones_demo.fasta"
OUTPUT_TSV="results/histones/histone_metadata.tsv"

echo "[INFO] Extracting metadata from FASTA headers..."

# Initialize file with Header
echo -e "Protein_Acc\tNucleotide_Acc\tLocation\tPartial_Info" > "$OUTPUT_TSV"

grep '^>' "$QUERY_FASTA" | while IFS= read -r header; do
    # Extract Nucleotide ID (e.g., NC_000001.1)
    nucleotide_id=$(echo "$header" | cut -d '|' -f 2 | awk -F '_cds' '{print $1}')

    # Extract Protein ID (e.g., NP_001001.1)
    protein_acc=$(echo "$header" | sed -n 's/.*protein_id=\([^]]*\).*/\1/p')

    # Extract Genomic Location string
    location=$(echo "$header" | sed -n 's/.*location=\([^]]*\).*/\1/p')

    # Identify if the sequence is partial
    partial_info="Complete"
    if [[ $header == *"partial=3'"* ]]; then
        partial_info="partial_3p"
    elif [[ $header == *"partial=5'"* ]]; then
        partial_info="partial_5p"
    fi

    echo -e "$protein_acc\t$nucleotide_id\t$location\t$partial_info" >> "$OUTPUT_TSV"
done

echo "[FINISH] Metadata extraction complete. Saved to $OUTPUT_TSV"
