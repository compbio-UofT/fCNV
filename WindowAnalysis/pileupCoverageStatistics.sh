#!/bin/bash
bamfile=../Data/__plasma.part.bam
bedfile=regions.bed
samtools mpileup $bamfile -q10 -Q10 -l $bedfile | awk '{print $2,$4}'
