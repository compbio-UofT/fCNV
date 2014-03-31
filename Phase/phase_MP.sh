#!/bin/bash
# --> genotype & phase the maternal and paternal genomes

# $1 - path to maternal bam files
# $2 - path to paternal bam files
# $3 - reference sequence .fa file
# $4 - region to extract

#logfile=log_phaseMP_$4.$$.log
#exec > $logfile 2>&1

date

if [ $# -ne 4 ]; then
    echo "Wrong number of arguments: got $#, expected 4"
    exit
fi

bindir="/dupa-filer/laci/bin"

#parse out chromosome ID
region=$4
chr=`echo $region | awk 'BEGIN {FS=":"}{print $1}'`
chrno=`echo $chr | awk 'BEGIN {FS="chr"}{print $2}'`
reference=$3

<<comment
# (1) extract the region and remove PCR duplicates
echo "merging and removing PCR duplicates:"
samtools merge -R $region - $1*/*/*.bam | samtools rmdup - __M.part.bam &
samtools merge -R $region - $2*/*/*.bam | samtools rmdup - __P.part.bam &
wait
samtools index __M.part.bam &
samtools index __P.part.bam &
wait
#samtools view -h __M.part.bam > __M.part.sam &
#samtools view -h __P.part.bam > __P.part.sam &
#wait
echo "-------- step 1 done ----------"
comment

# (2) genotype M, P, filter and phase
prefix=mp
echo "genotyping"
time samtools mpileup -uDSI -C50 -r $region -f $reference __M.allreads.part.bam __P.part.bam | bcftools view -bvcg - > $prefix.genotype.raw.bcf

time bcftools view $prefix.genotype.raw.bcf | vcfutils.pl varFilter -d50 -D150 -Q20 > $prefix.genotype.vcf
# TODO: ???? what limit for depth of coverage to use?

#annotate SNPs by rs-ids from dbSNP
echo  "annotating called SNPs"
java -Xmx24000m -jar ~/apps/snpEff/SnpSift.jar annotate -v /dupa-filer/laci/dbSnp.vcf $prefix.genotype.vcf > $prefix.genotype.annot.vcf
 
#extract only SNPs with reasonable quality score TODO: change qlimit depending on nummber of samples
snpsFile=$prefix.snps.annot
cat $prefix.genotype.annot.vcf | extract_annot_snps.awk -v qlimit=100 > $snpsFile.vcf

refpanels=/dupa-filer/laci/reference_panels/$chr.1kg.ref.phase1_release_v3.20101123.vcf.gz

#exclude markers that are not in the reference
echo  "conforming markers with the reference panels"
zcat $refpanels | conform.py $snpsFile.vcf /dev/stdin
sed "s/$chr/$chrno/g" $snpsFile.ftr.vcf > $snpsFile.ftr2.vcf

#phase by Beagle
echo  "invoking Beagle for phasing w/ reference panels"
time java -Xmx24000m -jar $bindir/b4.r1128.jar gt=$snpsFile.ftr2.vcf out=$prefix.phase ref=$refpanels impute=false

zcat $prefix.phase.vcf.gz > $prefix.phase.vcf

# CLEAN-UP
#rm __*.*



