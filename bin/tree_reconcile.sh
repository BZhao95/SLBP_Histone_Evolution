#!/bin/bash
iqtree -s SLBP_RBD_trimal.fasta -keep-ident -T 2 -m LG -fast

#reconciliation
treerecs \
    -g slbp_gene.nwk \
    -s euk_ncbi.nwk \
    -S treerecs_smap.tsv \
    -r \
    -O nhx \
    --fevent \
    -o treerecs_output
