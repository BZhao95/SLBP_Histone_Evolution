#!/bin/bash

genome_directory="Fungi/genomes"
out_dir="Fungi/braker2"

mkdir -p $out_dir

for f in "$genome_directory"/*.fna; do
    species=$(basename "$f" .fna)
    braker.pl --genome="$f" --species"$species" --workingdir="$out_dir/$species" --etpmode

done