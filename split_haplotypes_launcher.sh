#!/bin/bash

## This script is a wrapper around split_haplotypes.c and subsequent PLINK/GCTA usage.

###################
###   OPTIONS   ###
###################

compile=true
convert=true
makegrm=true
esth2=true

plinkfam=false
haveMothers=true
haveFathers=true

infile=~/Documents/gitrep/HAPLOTYPES/testcases/example_s3000_i6000.vcf.gz
infam=~/Documents/gitrep/HAPLOTYPES/testcases/example_s3000_i6000.fam
inpheno=~/Documents/gitrep/HAPLOTYPES/testcases/example_s3000_i6000.pheno
logfile=~/Documents/haplotypes/splitlog
outstem=~/Documents/haplotypes/split
grmlist=~/Documents/haplotypes/grmlist
finalout=~/Documents/haplotypes/gctatest

allres=~/Documents/haplotypes/allresults.txt

##################
###   SCRIPT   ###
##################

set -e

## change dir to allow relative paths further.
## not entirely safe, but at least parses symlinks
scriptdir=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$scriptdir"
echo "" > ${grmlist}

# function for adjusting GRMs by 2
# Args:
# 1 - input GRM, 2 - n samples
function adjustGRM {
	echo "Adjusting GRM ${1}.grm.bin..."
	./bin/divide_grm ${1}.grm.bin ${1}.grm.bin_temp $2
	mv ${1}.grm.bin_temp ${1}.grm.bin
	echo "${1}" >> ${grmlist}
}


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
	adjustGRM ${outstem}M1 $nlines

	echo "Producing paternal-transmitted GRM..."
	./bin/plink --bfile ${outstem}P1 --bim ${outstem}.bim --fam ${outstem}F.fam \
		--make-grm-bin --out ${outstem}P1
	adjustGRM ${outstem}P1 $nlines


	## if has full trios:
	if [ "$haveMothers" = true ];
	then
		echo "Producing maternal-untransmitted GRM..."
		./bin/plink --bfile ${outstem}M2 --bim ${outstem}.bim --fam ${outstem}M.fam \
			--make-grm-bin --out ${outstem}M2
		nlines=$(< "${outstem}M.fam" wc -l)
		adjustGRM ${outstem}M2 $nlines
	fi
	if [ "$haveFathers" = true ];
	then
		echo "Producing paternal-untransmitted GRM...."
		./bin/plink --bfile ${outstem}P2 --bim ${outstem}.bim --fam ${outstem}P.fam \
			--make-grm-bin --out ${outstem}P2
		nlines=$(< "${outstem}P.fam" wc -l)
		adjustGRM ${outstem}P2 $nlines
	fi

	echo "All GRMs produced."
fi

if [ "$esth2" = true ];
then
	echo "Estimating h2 from GRMs..."
	./bin/gcta64 --reml --mgrm ${grmlist} --pheno ${inpheno} --out ${finalout}
	echo "analysis name: $1" >> ${allres}
	grep '/Vp' ${finalout}.hsq >> ${allres}
fi
