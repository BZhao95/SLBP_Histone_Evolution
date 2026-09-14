#!/bin/bash

set -euo pipefail
#
#===============================================================
#script: 01_remote_blast.sh
#decription: Perform batched remote BLASTP searches and avoid
#            extensive computational warning from NCBI.
#
#===============================================================


# --- User Configuration --- 
INPUT_FASTA="./data/slbp.fa"
OUTPUT_DIR="./results/blastp_results"
FINAL_OUTPUT="./results/blastp_results.txt"
BATCH_SIZE=5 #ncbi limits
DB="nr" #change to another ncbi database like refseq if needed
#ENTREZ_QUERY="Eukaryota[Organism]" #Change it to Bacteria, Gaint Virus if needed
ENTREZ_QUERY="Cryptomycota[Organism]" 
EVALUE="1e-10"
MAX_TARGETS=100000

#BLAST Output format
FORMAT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids"

#Create output directory
mkdir -p "$OUTPUT_DIR"


# ---Functions---
# split large query files into smaller batches to prevent remote BLASTP timeouts

split_fasta() {
    local input_file=$1
    local output_prefix=$2
    local batch_size=$3
    local count=0
    local batch_number=1
    local output_file="${output_prefix}_$(printf "%03d" $batch_number).fa"

    while read -r line; do
        if [[ $line == ">"* ]]; then
            ((count++))
            if ((count>batch_size)); then
                batch_number=$((batch_number +1))
                count=1
                output_file="${output_prefix}_$(printf "%03d" $batch_number).fa"
            fi
        fi
        echo "$line" >> "$output_file"
    done < "$input_file"

}

#---Excution---

echo "[INFO] Initializating sequence batching (Size: $BATCH_SIZE)..."
split_fasta "$INPUT_FASTA" "${OUTPUT_DIR}/batch" "$BATCH_SIZE"

echo "[INFO] Connecting remote BLASTP seaches..."

for BATCH_FILE in "${OUTPUT_DIR}"/batch_*.fa; do
    OUTPUT_FILE="${BATCH_FILE%.fa}.txt"
    echo "Processing $BATCH_FILE ..."

    blastp -query "$BATCH_FILE" \
            -db "$DB" \
            -remote \
            -evalue "$EVALUE" \
            -out "$OUTPUT_FILE" \
            -outfmt "$FORMAT" \
            -entrez_query "$ENTREZ_QUERY" \
            -max_target_seqs "$MAX_TARGETS"
    

    echo "Successfully processed $(basename "$BATCH_FILE")."
    

done

#---Merge---
echo "Merging results into a single output file"
cat ${OUTPUT_DIR}/batch_*.txt > "$FINAL_OUTPUT"

echo "[FINISH] BLASTP searches completed. Results saved to $FINAL_OUTPUT"