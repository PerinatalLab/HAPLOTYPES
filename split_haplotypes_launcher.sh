#!/bin/bash

set -e

## This script is just a wrapper around split_haplotypes.c and subsequent PLINK/GCTA usage.

compile=true;
if [ "$compile" = true ];
then
	g++ split_haplotypes.c -o bin/split_haplotypes
fi

infile = ~/Documents/haplotypes/test.vcf.gz
logfile = ~/Documents/haplotypes/splitlog
outstem = ~/Documents/haplotypes/split

echo "Starting run..."

zcat ${infile} | ./bin/split_haplotypes ${infam} ${logfile} ${outstem}

echo "Run finished."

## produce GRM matrices for each haplotype
## (use --keep-allele-order if want to keep ref/alt assignments)
./bin/plink --bfile ${outstem}M1 --bim ${outstem}.bim --fam ${outstem}F.fam \
	--make-grm-bin --out ${outstem}M1
./bin/plink --bfile ${outstem}P1 --bim ${outstem}.bim --fam ${outstem}F.fam \
	--make-grm-bin --out ${outstem}P1
./bin/plink --bfile ${outstem}M2 --bim ${outstem}.bim --fam ${outstem}M.fam \
	--make-grm-bin --out ${outstem}M2
./bin/plink --bfile ${outstem}P2 --bim ${outstem}.bim --fam ${outstem}P.fam \
	--make-grm-bin --out ${outstem}P2
