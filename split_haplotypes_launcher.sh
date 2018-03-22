#!/bin/bash

set -e

## This script is a wrapper around split_haplotypes.c and subsequent PLINK/GCTA usage.

###################
###   OPTIONS   ###
###################

compile=true;
convert=true;
makegrm=true;
esth2=true;

infile=~/Documents/gitrep/HAPLOTYPES/testcases/example_s3000_i18000.vcf.gz
infam=~/Documents/gitrep/HAPLOTYPES/testcases/example_s3000_i18000.fam
logfile=~/Documents/haplotypes/splitlog
outstem=~/Documents/haplotypes/split


##################
###   SCRIPT   ###
##################

set -e

if [ "$compile" = true ];
then
	echo "Compiling script..."
	g++ split_haplotypes.c -o bin/split_haplotypes
	echo "Compilation complete."
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
	echo "Producing paternal-transmitted GRM..."
	./bin/plink --bfile ${outstem}P1 --bim ${outstem}.bim --fam ${outstem}F.fam \
		--make-grm-bin --out ${outstem}P1

	## if has full trios:
	if [ "$haveMothers" = true ];
	then
		echo "Producing maternal-untransmitted GRM..."
		./bin/plink --bfile ${outstem}M2 --bim ${outstem}.bim --fam ${outstem}M.fam \
			--make-grm-bin --out ${outstem}M2
	fi
	if [ "$haveFathers" = true ];
	then
		echo "Producing paternal-untransmitted GRM...."
		./bin/plink --bfile ${outstem}P2 --bim ${outstem}.bim --fam ${outstem}P.fam \
			--make-grm-bin --out ${outstem}P2
	fi

	echo "All GRMs produced."
fi

if [ "$esth2" = true ];
then
	echo "Estimating h2 from GRMs..."
	echo "(not implemented yet)"
fi
