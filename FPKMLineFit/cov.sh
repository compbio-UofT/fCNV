sz=1000000
begin=2300000
centroBegin=120600000
centroEnd=147000000
end=243700000

rm scaterData

for (( i=$begin; i <= $end-$sz; i+=sz ))
do
	if (( (i <= centroEnd && i >= centroBegin) || (i+sz <= centroEnd && i+sz >= centroBegin) ))
	then
		continue
	fi
	tmpEnd=$(( i+sz ))
	plasmaG1=`(samtools view /dupa-filer/laci/G1/chr1/__plasma.part.bam chr1:$i-$tmpEnd &) | wc -l`
	plasmaI1=`(samtools view /dupa-filer/laci/I1/chr1/__plasma.part.bam chr1:$i-$tmpEnd &) | wc -l`
	maternal=`(samtools view /dupa-filer/laci/I1/chr1/__M.part.bam chr1:$i-$tmpEnd &) | wc -l`
	wait
	echo $plasmaG1" "$plasmaI1" "$maternal >> scaterData
done
