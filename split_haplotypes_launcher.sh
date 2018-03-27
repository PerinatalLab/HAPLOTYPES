#!/bin/bash

## This script is a wrapper around split_haplotypes.c and subsequent PLINK/GCTA usage.

###################
###   OPTIONS   ###
###################

compile=true
convert=true
makegrm=true
esth2=false

plinkfam=false
haveMothers=true
haveFathers=true

infile=~/Documents/gitrep/HAPLOTYPES/testcases/example_s3000_i6000.vcf.gz
infam=~/Documents/gitrep/HAPLOTYPES/testcases/example_s3000_i6000.fam
inpheno=~/Documents/gitrep/HAPLOTYPES/testcases/example_s3000_i6000.pheno
logfile=~/Documents/haplotypes/splitlog
outstem=~/Documents/haplotypes/split


##################
###   SCRIPT   ###
##################

set -e

## change dir to allow relative paths further.
## not entirely safe, but at least parses symlinks
scriptdir=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$scriptdir"

if [ "$compile" = true ];
then
	echo "Compiling script..."
	g++ split_haplotypes.c -o bin/split_haplotypes
	gcc divide_grm.c -o bin/divide_grm
	echo "Compilation complete."
fi

if [ "$plinkfam" = true ];
then
	echo "Converting .fam to pedigree file..."
	## drop row if this individual is ever mentioned as a parent
	## (in effect, trims multigenerational pedigrees)
	awk -F'\t' '{r[$2]=$2 FS $3 FS $4; p[$3]; m[$4]} END{for(i in r) if(!(i in p || i in m)) print r[i] }' ${infam} > ${outstem}_temp.fam
	infam=${outstem}_temp.fam
	echo "Pedigree file produced."
fi

if [ "$convert" = true ];
then
	echo "Converting VCFs to haplotype BEDs..."
	zcat ${infile} | ./bin/split_haplotypes ${infam} ${logfile} ${outstem}
	echo "Haplotype BEDs produced."
fi

if [ "$makegrm" = true ];
then
	echo "Calculating GRMs from BEDs..."

	## produce GRM matrices for each haplotype
	## (use --keep-allele-order if want to keep ref/alt assignments)
	echo "Producing maternal-transmitted GRM..."
	./bin/plink --bfile ${outstem}M1 --bim ${outstem}.bim --fam ${outstem}F.fam \
		--make-grm-bin --out ${outstem}M1
	echo "Dividing GRM M1 by 2..."
	nlines=$(< "${outstem}F.fam" wc -l)
	./bin/divide_grm ${outstem}M1.grm.bin ${outstem}M1_t.grm.bin $nlines
	mv ${outstem}M1_t.grm.bin ${outstem}M1.grm.bin

	echo "Producing paternal-transmitted GRM..."
	./bin/plink --bfile ${outstem}P1 --bim ${outstem}.bim --fam ${outstem}F.fam \
		--make-grm-bin --out ${outstem}P1
	./bin/divide_grm ${outstem}P1.grm.bin ${outstem}P1_t.grm.bin $nlines
	mv ${outstem}P1_t.grm.bin ${outstem}P1.grm.bin


	## if has full trios:
	if [ "$haveMothers" = true ];
	then
		echo "Producing maternal-untransmitted GRM..."
		./bin/plink --bfile ${outstem}M2 --bim ${outstem}.bim --fam ${outstem}M.fam \
			--make-grm-bin --out ${outstem}M2
		echo "Dividing GRM M2 by 2..."
		nlines=$(< "${outstem}M.fam" wc -l)
		./bin/divide_grm ${outstem}M2.grm.bin ${outstem}M2_t.grm.bin $nlines
		mv ${outstem}M2_t.grm.bin ${outstem}M2.grm.bin
	fi
	if [ "$haveFathers" = true ];
	then
		echo "Producing paternal-untransmitted GRM...."
		./bin/plink --bfile ${outstem}P2 --bim ${outstem}.bim --fam ${outstem}P.fam \
			--make-grm-bin --out ${outstem}P2
		echo "Dividing GRM P2 by 2..."
		nlines=$(< "${outstem}P.fam" wc -l)
		./bin/divide_grm ${outstem}P2.grm.bin ${outstem}P2_t.grm.bin $nlines
		mv ${outstem}P2_t.grm.bin ${outstem}P2.grm.bin
	fi

	echo "All GRMs produced."
fi

if [ "$esth2" = true ];
then
	echo "Estimating h2 from GRMs..."
#	./bin/gcta64 --
fi
