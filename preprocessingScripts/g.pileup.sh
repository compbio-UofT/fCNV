#!/bin/bash

samtools mpileup -Q10 -q10 -l regions.bed /dupa-filer/laci/G1/$1/__plasma.part.bam | awk '{print $2,$4}' > g.plasma.txt
