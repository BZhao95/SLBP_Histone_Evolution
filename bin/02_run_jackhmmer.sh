#!/bin/bash
set -euo pipefail

#file path
query_file="slbp.fa"
pdb_db="pdb_seqres.fasta" #download from PDB
output_dir="jackhmmer_results"

#create output folder if it does not exit
mkdir -p "$output_dir"

#run jackhmmer
jackhmmer -N5 --cpu 4 \
    --tblout "$output_dir/jackhmmer_tblout.txt" \
    --domtblout "$output_dir/jackhmmer_dombtlout.txt" \
    -o "$output_dir/jackhmmer_out_put.txt" \
    "$query_file" "$pdb_db"

echo "jackhmmer search done!"

