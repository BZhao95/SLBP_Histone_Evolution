#!/bin/bash

proteome_dir="Fungi/proteomes"
genome_dir="Fungi/genomes"
summary_file="Fungi/fungal_proteome_summary.tsv"
skipped_file="Fungi/skipped_species.txt"

mkdir -p "$proteome_dir" "$genome_dir"

echo -e "Species\tProteome_Status\tGenome_Status" > "$summary_file"
echo -e "Species\tReason" > "$skipped_file"


while IFS= read -r species_name; do

    echo "Processing $species_name ..."

    assembly_info=$(esearch -db assembly -query "\"$species_name\"" |esummary| \
        xtract -pattern DocumentSummary -element FtpPath_GenBank -element RefSeq_category | head -n 1)

    if [[ -z "$assembly_info" ]]; then
        echo -e "$species_name\tNo assembly" >> "$skipped_file"
        continue
    fi

    ftp_path=$(echo "$assembly_info" | cut -f1)

    if [[ -z "$ftp_path" || "$ftp_path" == "na" ]]; then 
        echo -e "$species_name\tInvalid FTP Path" >> "$skipped_file"
        continue
    fi

    assembly_id=$(basename "$ftp_path")

    proteome_url="${ftp_path}/${assembly_id}_protein.faa.gz"
    genome_url="${ftp_path}/${assembly_id}_genomic.fna.gz"

    proteome_file="$proteome_dir/${species_name// /_}.faa.gz"
    genome_file="$genome_dir/${species_name// /_}fna.gz"

    wget -q "$proteome_url" -O "$proteome_file"

    if [[ -s "$proteome_file" ]]; then
        gunzip -f "$proteome_file"
        echo -e "$species_name\tAvailable\t-" >> "$summary_file"
        continue
    else
        rm -f "$proteome_file"

    fi

    wget -q "$genome_url" -O "$genome_file"

    if [[ -s "$genome_file" ]]; then
    gunzip -f "$genme_file"
        echo -e "$species_name\tNo\tAvaliable" >> "$summary_file"
    else
        rm -f "$genome_file"
    fi

done < fungal_taxid_to_species.txt