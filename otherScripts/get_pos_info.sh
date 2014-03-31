#!/bin/bash
for pos in $@ ; do
    sed -n $pos'p' positions.txt
    chrPos=`sed -n $pos'p' positions.txt | awk '{print $1}'`
    for file in M P F; do
        echo "+++++++++++++++> info for "$file":"
        echo "  VCFftr: "`grep "[^0-9]$chrPos[^0-9]" __$file.part.snps.vcf`
        echo "  VCFraw: "`grep "[^0-9]$chrPos[^0-9]" __$file.part.genotype.vcf`
        echo "  pileup: "`samtools mpileup -r chr20:$chrPos-$chrPos __$file.part.bam 2> /dev/null`
    done
    echo "----------------------------------------------------------------------"
done
