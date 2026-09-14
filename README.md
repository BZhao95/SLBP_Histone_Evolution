# Evolutionary of Histone mRNA 3' End Processing Machinery in Eukaryotes

A comprehensive pipeline for identifying and analyzing Stem-Loop Binding Protein (SLBP) orthologs and replication-dependent histone processing elements across eukaryotes.

## Repository Structure
```text
Histone_SLBP_Project/
|- bin/			# scripts
|- data/		# input FASTA
|- results/		# outputs
|- logs/		# errors
|- requirements.txt	# dependencies
|- Readme.md
```
## Pipeline Overview
1. **Protein Identification**: Remote BLASTP and jackhmmer searches for SLBP, LSM10/11, and FLASH.
2. **Domain Verification**: RPS-BLAST and HMMER profiling for RBD (PF15247) and HFD domains.
3. **RNA Analysis**: Detection of histone 3' UTR Stem-Loops and PAS motifs.
4. **Evolutionary Synthesis**: Gene-species tree reconciliation and processing state classification.

## Requirements
- NCBI BLAST+
- HMMER 3.x
- RNAfold (ViennaRNA)
- Python 3.8+
- cmsearch
- MAFFT
- FastTree
- ete3
- Numpy
- networks
- matplotlib


