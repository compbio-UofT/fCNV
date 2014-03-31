#!/bin/bash
# --> genotype & phase the trio

# $1 - path to maternal bam files
# $2 - path to paternal bam files
# $3 - path to fetal bam files
# $4 - reference sequence .fa file
# $5 - region to extract

#logfile=log_phaseTrio_$5.$$.log
#exec > $logfile 2>&1

date

if [ $# -ne 5 ]; then
    echo "Wrong number of arguments: got $#, expected 5"
    exit
fi

bindir="/dupa-filer/laci/bin"

#parse out chromosome ID
region=$5
chr=`echo $region | awk 'BEGIN {FS=":"}{print $1}'`
chrno=`echo $chr | awk 'BEGIN {FS="chr"}{print $2}'`
reference=$4

<<comment
# (1) extract the region and remove PCR duplicates
echo "merging and removing PCR duplicates:"
samtools merge -R $region - $1*/*/*.bam | samtools rmdup - __M.part.bam &
samtools merge -R $region - $2*/*/*.bam | samtools rmdup - __P.part.bam &
samtools merge -R $region - $3*/*/*.bam | samtools rmdup - __F.part.bam &
wait
samtools index __M.part.bam &
samtools index __P.part.bam &
samtools index __F.part.bam &
wait
#samtools view -h __M.part.bam > __M.part.sam &
#samtools view -h __P.part.bam > __P.part.sam &
#wait
echo "-------- step 1 done ----------"
comment

# (2) genotype M, P, F, filter and phase
prefix=trio
echo "genotyping the trio"
time samtools mpileup -uDSI -C50 -r $region -f $reference __M.allreads.part.bam __P.part.bam __F.part.bam | bcftools view -bvcg - > $prefix.genotype.raw.bcf

time bcftools view $prefix.genotype.raw.bcf | vcfutils.pl varFilter -d60 -D200 -Q20 > $prefix.genotype.vcf
# TODO: ???? what limit for depth of coverage to use?

#annotate SNPs by rs-ids from dbSNP
echo  "annotating called SNPs"
java -Xmx24000m -jar ~/apps/snpEff/SnpSift.jar annotate -v /dupa-filer/laci/dbSnp.vcf $prefix.genotype.vcf > $prefix.genotype.annot.vcf
 
#extract only SNPs with reasonable quality score TODO: change qlimit depending on nummber of samples
snpsFile=$prefix.snps.annot
cat $prefix.genotype.annot.vcf | extract_annot_snps.awk -v qlimit=100 > $snpsFile.vcf

#compress and index
#bgzip $prefix.snps.annot.vcf
#tabix -p vcf $prefix.snps.annot.vcf.gz

refpanels=/dupa-filer/laci/reference_panels/$chr.1kg.ref.phase1_release_v3.20101123.vcf.gz

#modify our trio VCF file so that its records are consistent with the reference VCF file
#...this doesn't seems to work, rather use custom implementation - conform.py
#time java -jar ~/apps/jar/conform-gt.jar  gt=$prefix.snps.annot.vcf.gz out=conform chrom=$chr ref=$refpanels

#exclude markers that are not in the reference
echo  "conforming markers with the reference panels"
zcat $refpanels | conform.py $snpsFile.vcf /dev/stdin
sed "s/$chr/$chrno/g" $snpsFile.ftr.vcf > $snpsFile.ftr2.vcf

#phase by Beagle
echo  "invoking Beagle for trio phasing w/ reference panels"
time java -Xmx24000m -jar $bindir/b4.r1128.jar gt=$snpsFile.ftr2.vcf ped=pedigree.txt out=$prefix.phase ref=$refpanels impute=false
#time java -jar $bindir/b4.r1128.jar gt=$snpsFile.vcf.gz ped=pedigree.txt out=$prefix.phase

zcat $prefix.phase.vcf.gz > $prefix.phase.vcf

# CLEAN-UP
#rm __*.*



