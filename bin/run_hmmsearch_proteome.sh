#!/bin/bash

proteome_dir="Fungi/proteomes"
hmm_model="data/PF15247.hmm"
output_dir="Fungi/Fungi_proteome_hmm_results"

mkdir -p "$output_dir"

echo -e "Species\tHit_count" > "$output_dir/SLBP_RNA_Bind_hits.txt"

for f in "$proteome_dir"/*.faa; do
    species=$(basename "$f" .faa)
    out="$output_dir/{$species}.tbl"

    hmmsearch --tblout "$out" --cut_ga "$hmm_model" "$f"

    hit=$(grep -vc '^#' "$out")

    if [[ "$hit" -gt 0 ]]; then
        echo -e "$species\t$hit" >> "$output_dir/SLBP_RNA_bind_hits.tsv"
    else
        rm - f "$out"
    fi
done
