#!/bin/bash

tmp_dir=/tmp/
if [[ -n "$TMPDIR" ]]; then
    tmp_dir=$TMPDIR
fi

plasmaFetusRate=$1
chrom=$2
phase_sites=/dupa-filer/laci/I1/$chrom/trio.phase.vcf
bamfile=/dupa-filer/laci/I1/$chrom/__plasma.part.bam
samfile=$tmp_dir/__plasma.$chrom.sam
output_bamfile=$tmp_dir/'__plasma-'$chrom'-'$plasmaFetusRate'rate.sort.bam'
output_samfile=$tmp_dir/'__plasma-'$chrom'-'$plasmaFetusRate'rate.sort.sam'

echo "sim_fetusRate.sh: Starting down-rating $chrom to $plasmaFetusRate"
date

# Copying the header
samtools view -H $bamfile -o $output_samfile
samtools view $bamfile -L /dupa-filer/laci/I1/$chrom/regions.bed -o $samfile

echo "Filtering the reads to reduce fetal DNA rate"
pypy /dupa-filer/laci/bin/fetusRateAdjuster.py $samfile $phase_sites $plasmaFetusRate >> $output_samfile
rm $samfile
samtools view -S -b $output_samfile -o $output_bamfile
rm $output_samfile

echo "Indexing"
samtools index $output_bamfile

mv $output_bamfile* /dupa-filer/laci/I1/$chrom/

echo "sim_fetusRate script - DONE."

