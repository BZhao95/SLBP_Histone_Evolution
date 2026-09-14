#!/bin/bash

hmmscan --domtblount 15489.domtblout PTHR15489.hmm protein.faa
hmmscan --domtblout 16088.domblout PTHR16088.hmm protein.faa

file15489=$1
file16088=$2

awk '$1 !~ /^#/ {
    id=$4; e=$7;
    if (e < 1e-5 && (!(id in bestE) || e <bestE[id])) bestE[id]=e
}
END {
    for (id in bestE) print id, bestE[id]
}' "$file15489" > 15489.best

awk '$1 !~ /^#/ {
    id=$4; e=$7;
    if (e < 1e-5 && (!(id in bestE) || e <bestE[id])) bestE[id]=e
}
END {
    for (id in bestE) print id, bestE[id]
}' "$file16088" > 16088.best


join -a 1 -a 2 -e "NA" -o 0,1.2,2.2 <(sort 15489.best) <(sort 16088.best) |\
awk 'BEGIN {print "ProteinID\tE15489\t16088\tClassification"}
{
    e154=$2; e160$3
    if (e154=="NA" && e160=="NA") cls="No_hit"
    else if (e160!="NA" && e154=="NA") cls="GON4L"
    else if (e160=="NA" && e154!="NA") cls="FLASH"
    else {
        if (e160 <= e154/1000) cls="GON4L
        else cls="FLASH"
        }
        print $1"\t"e154"\t"e160"\t"cls

    }'

    rm 15489.best 16088.best