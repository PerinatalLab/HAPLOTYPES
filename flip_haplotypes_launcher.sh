#!/bin/bash
set -e

## This script is just a wrapper around the flip_haplotypes.c script.
compile=true;
if [ "$compile" = true ];
then
	g++ flip_haplotypes.c -o flip_haplotypes -lz
fi

infile12=/media/local-disk2/jjuod/tsdphased/imputed_m12/
infile24=/media/local-disk2/jjuod/tsdphased/imputed_m24/
outdir=/media/local-disk2/jjuod/tsdphased/fixedhaplos/

mkdir -p ${outdir}

echo "Starting..."

# USAGE: zcat infile.gz | ./flip_haplotypes famfile logfile outfile2 | gzip > outfile1.gz
for chr in {1..22}
do
	echo "working on chromosome ${chr} in moba12"
	zcat ${infile12}${chr}.vcf.gz | ./flip_haplotypes \
		/media/local-disk/jonasbac/genetic_score/haploPrepFam_m12_n6636_19c7de9_20171201.txt \
		${outdir}log${chr}.txt \
		${outdir}haploMomUnT_m12_n6636_19c7de9_${chr}.txt.gz | gzip > ${outdir}haploPrepFam_m12_n6636_19c7de9_${chr}.vcf.gz

	echo "working on chromosome ${chr} in moba24"
	zcat ${infile24}${chr}.vcf.gz | ./flip_haplotypes \
		/media/local-disk/jonasbac/genetic_score/haploPrepFam_m24_n4271_19c7de9_20171201.txt \
		${outdir}log${chr}.txt \
		${outdir}haploMomUnT_m24_n4271_19c7de9_${chr}.txt.gz | gzip > ${outdir}haploPrepFam_m24_n4271_19c7de9_${chr}.vcf.gz

done
echo "Completed successfully."
