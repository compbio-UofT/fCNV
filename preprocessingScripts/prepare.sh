#!/bin/bash

# $1 - maternal reads bam file
# $2 - paternal reads bam file
# $3 - fetal reads bam file
# $4 - reference sequence .fa file
# $5 - region to extract

if [ $# -ne 5 ]; then
    echo "Wrong number of arguments: got $#, expected 5"
    exit
fi

#parse out chromosome ID
chr=`echo $5 | awk 'BEGIN {FS=":"}{print $1}'`

<<comment
# (1) remove PCR duplicates
echo "removing PCR duplicates:"
samtools view -bu $1 $5 | samtools rmdup - - > __M.part.bam &
samtools view -bu $2 $5 | samtools rmdup - - > __P.part.bam &
samtools view -bu $3 $5 | samtools rmdup - - > __F.part.bam &
wait
samtools index __M.part.bam &
samtools index __P.part.bam &
samtools index __F.part.bam &
wait
comment
echo "-------- step 1 done ----------"

# (2) genotype M, P, F, filter and phase

for file in __M.part __P.part __F.part
do
    echo "genotyping $file"
    samtools mpileup -uD -r $5 -f $4 $file.bam | bcftools view -bvcg - > $file.genotype.raw.bcf &
done
wait

for file in __M.part __P.part __F.part
do
    bcftools view $file.genotype.raw.bcf | vcfutils.pl varFilter -D100 > $file.genotype.vcf &
    # ???? what limit for depth of coverage to use?
done
wait

for file in __M.part __P.part __F.part
do
    #annotate SNPs by rs-ids from dbSNP
    echo  "annotating $file"
    java -jar ~/apps/snpEff/SnpSift.jar annotate -v dbSnp.vcf $file.genotype.vcf > $file.genotype.annot.vcf &
done
wait

for file in __M.part __P.part __F.part
do
    #extract only SNPs with reasonable quality score
    cat $file.genotype.annot.vcf | ./extract_annot_snps.awk -v qlimit=0 > $file.snps.annot.vcf &
done
wait

#----------------
# filter out SNP positions for which we don't have confident calls (or homoz. ref.) evidance in both M and P
echo "Filtering SNP positions"
time ./filter_vcfs.py __M.part.snps.annot.vcf __P.part.snps.annot.vcf

refpanels=ALL.$chr.phase1_release_v3.20101123.filt
for prefix in __M __P
do
    #convert annotated VCF to BEAGLE format (creates $file.bgl.gz, .markers, and .int)
    cat $prefix.part.snps.annot.ftr.vcf | java -jar ~/apps/jar/vcf2beagle.jar ? $prefix
    #unify alleles for the markers in $file.markers
    pypy union_markers.py $prefix.markers $refpanels.markers > $prefix.mod.markers
    #phase haplotypes
    java -Xmx10000m -jar ~/apps/jar/beagle.jar markers=$prefix.mod.markers phased=$refpanels.bgl.gz  unphased=$prefix.bgl.gz missing=? out=$prefix omitprefix=true &
done
wait
        
for prefix in __M __P
do
    #convert haplotypes from BEAGLE to VCF format
    java -jar ~/apps/jar/beagle2vcf.jar $chr $prefix.markers $prefix.bgl.gz.phased.gz ? > $prefix.phased.vcf &
done
wait

echo "-------- step 2 done ----------"
#comment

# (3) mix reads to get plasma-like reads
for gnm in M F; do
    echo "getting number of reads and coverage for $gnm"
    count_var=$gnm'_count'
    eval $count_var=$(samtools view __$gnm.part.bam | wc -l)
    
    coverage_var=$gnm'_coverage'
    eval $coverage_var=$(samtools mpileup -D -r $5 __$gnm.part.bam | awk '{sum+=$4; count+=1} END {print sum/count}')
        echo "  >> got ${!count_var} and ${!coverage_var}"
done

target_coverage=80
mixture=0.10

M_multiplier=`echo "scale=5; ($target_coverage * (1.0  - $mixture)) / $M_coverage"|bc`
F_multiplier=`echo "scale=5; ($target_coverage * $mixture) / $F_coverage"|bc`

#for debugging set multipliers by hard
#M_multiplier=0.543
#F_multiplier=0.123
echo "we need $M_multiplier x M and $F_multiplier x F reads" 

#sample reads to simulate plasma reads
temp_file='__temp.sam'
file_count=-1
samtools view -H __M.part.bam > $temp_file
samtools view -H __F.part.bam >> $temp_file

for gnm in M F; do
    frac=$gnm'_multiplier'
    frac=${!frac}
    #echo "+++> $frac"

    while [ $(bc <<< "$frac >= 1") -eq 1 ]; do
        file_count=$(($file_count+1))
        #echo $frac $gnm
        samtools view __$gnm.part.bam | awk -v f=$frac '{OFS="\t"; $1=$1"."f; print $0}' > $temp_file.$gnm.$file_count &

        frac=`echo "scale=5; $frac - 1"|bc`
        #echo $frac
    done
    if [ $(bc <<< "$frac > 0") -eq 1 ]; then
        file_count=$(($file_count+1))
        #echo $frac $gnm
        samtools view -s $frac __$gnm.part.bam | awk -v f=$frac '{OFS="\t"; $1=$1"."f; print $0}' > $temp_file.$gnm.$file_count &
    fi    
done
wait

#concatenate the temp files
#for (( c=0; c<=$file_count; c++ )); do
#    cat $temp_file$c >> $temp_file
#done
for gnm in M F; do
    samtools view -H __$gnm.part.bam > __plasma.$gnm.sam
    cat $temp_file.$gnm.* >> __plasma.$gnm.sam
done
cat $temp_file.[MF].* >> $temp_file

#sort the sam file: turn it to a bam, sort, convert back to sam
samtools view -Sbu $temp_file | samtools sort - __plasma.sort
rm $temp_file*
plasma_file='plasma.sort.sam'
samtools view -h -q 20  __plasma.sort.bam > $plasma_file
rm __plasma.sort.bam

echo "-------- step 3 done ----------"

# (4) get nucleotides counts for individual SNP positions (union of the SNP positions found in (2)
# when genotyping the parental genomes)
echo "Processing final VCFs and plasma reads to FCNV input files"
time pypy process_vcfs.py __M.phased.vcf __P.phased.vcf __F.part.snps.annot.vcf $plasma_file

echo "-------- step 4 done ----------"

# CLEAN-UP
#rm __*.*



