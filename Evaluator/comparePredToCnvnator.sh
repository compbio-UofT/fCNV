#!/bin/bash

G=$1

for x in {1..22}
do
    echo ">>> Chr$x <<<"
    echo -e "PREDICTION\t\t\t\tCNVNATOR\t\t\t\tOVERLAP"
    pred_bed=/tmp/__temp.chr$x.bed
    annotation_evaluation.py *chr$x-*annotation* 30 | grep predicted | awk -v x=$x '{split($1,a,"-"); printf "chr%s\t%d\t%d\n", x, a[1], a[2]}' > $pred_bed
    intersectBed -a $pred_bed -b /dupa-filer/laci/I1/chr$x/cnvnatorCalls/$G.bed -wo
    
    echo ""
    rm $pred_bed
done
