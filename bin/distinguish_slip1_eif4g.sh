#!/bin/bash

hmmscan --domtblout 23254.domtblout PTHR23254.hmm protein.faa
hmmscan --domtblout 23253.domtblout PTHR23253.hmm protein.faa

bash distinguish_slip1_eif4g.sh 23254.domtblout 23253.domtblout > slip1_classification.tsv

file23254=$1
file23253=$2

awk '$1 !~ /^#/ {
    id=$4; e=$7
    if (e < 1e-5 && (!(id in bestE) || e < bestE[id])) bestE[id]=e
}
END {
    for (id in bestE) print id, bestE[id]
}' "$file23254" > 23254.best

awk '$1 !~ /^#/ {
    id=$4; e=$7
    if (e < 1e-5 && (!(id in bestE) || e < bestE[id])) bestE[id]=e
}
END {
    for (id in bestE) print id, bestE[id]
}' "$file23253" > 23253.best

join -a 1 -a 2 -e "NA" -o 0,1.2,2.2 <(sort 23254.best) <(sort 23253.best) | \
awk 'BEGIN {print "ProteinID\tPTHR23254\tPTHR23253\tClassification"}
{
    e23254=$2; e23253=$3
    if (e23254=="NA" && e23253=="NA") cls="No_hit"
    else if (e23253!="NA" && e23254=="NA") cls="eIF4G"
    else if (e23254!="NA" && e23253=="NA") cls="SLIP1"
    else {
        if (e23253+0 <= e23254/1000) cls="eIF4G"
        else cls="SLIP1"
    }
    print $1"\t"e23254"\t"e23253"\t"cls
}'

rm 23254.best 23253.best