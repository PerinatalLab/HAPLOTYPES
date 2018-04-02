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

# constant for genotype-based GRM adjustment
divc=1

infile=~/Documents/gitrep/HAPLOTYPES/testcases/example_s3000_i6000.vcf.gz
infam=~/Documents/gitrep/HAPLOTYPES/testcases/example_s3000_i6000.fam
inpheno=~/Documents/gitrep/HAPLOTYPES/testcases/example_s3000_i6000.pheno

workdir=~/Documents/haplotypes/tmp/
allres=~/Documents/haplotypes/allresults.txt


##################
###   SCRIPT   ###
##################

set -e

## change dir to allow relative paths further.
## not entirely safe, but at least parses symlinks
scriptdir=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$scriptdir"

# initialize result file header
sed -i 's/^$/Source\tVariance\tSE\tAnalysis\tMethod/' ${allres}

mkdir -p ${workdir}
outstem=${workdir}split
logfile=${outstem}.log
grmlist=${outstem}grmlist
grmlistG=${outstem}grmlistG
finalout=${outstem}gctatest

# function for adjusting GRMs by 2
# Args:
# 1 - input GRM, 2 - n samples, 3 - division constant
function adjustGRM {
	echo "Adjusting GRM ${1}.grm.bin by a factor of ${3}..."
	./bin/divide_grm ${1}.grm.bin ${1}.grm.bin_temp $2 ${3} 
	mv ${1}.grm.bin_temp ${1}.grm.bin
}

# function for parsing GCTA log into a table
# Args:
# 1 - input file, 2 - analysis name
function parseRes {
	awk -v OFS="\t" -v NAME=$2 '/Source\tVariance/{p=1; next}
		$0~/^$/ && p==1{exit}
		p==1{print $0, NAME, "AI-REML"}' $1 > ${allres}
}


##################
###    MAIN    ###
##################

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
	nlines=$(< "${outstem}F.fam" wc -l)
	adjustGRM ${outstem}M1 $nlines 2

	echo "Producing paternal-transmitted GRM..."
	./bin/plink --bfile ${outstem}P1 --bim ${outstem}.bim --fam ${outstem}F.fam \
		--make-grm-bin --out ${outstem}P1
	adjustGRM ${outstem}P1 $nlines 2


	## produce standard GRM for fetal genotype
	echo "Producing fetal-genotype GRM..."
	./bin/plink --vcf ${infile} --keep ${outstem}F.fam \
		--make-grm-bin --out ${outstem}FG
	adjustGRM ${outstem}FG $nlines $divc 


	## if has full trios:
	if [ "$haveMothers" = true ];
	then
		echo "Producing maternal-untransmitted GRM..."
		./bin/plink --bfile ${outstem}M2 --bim ${outstem}.bim --fam ${outstem}M.fam \
			--make-grm-bin --out ${outstem}M2
		nlines=$(< "${outstem}M.fam" wc -l)
		adjustGRM ${outstem}M2 $nlines 2
	
		echo "Producing maternal-genotype GRM..."
		./bin/plink --vcf ${infile} --keep ${outstem}M.tokeep \
			--make-grm-bin --out ${outstem}MG
		adjustGRM ${outstem}MG $nlines $divc
	fi
	if [ "$haveFathers" = true ];
	then
		echo "Producing paternal-untransmitted GRM...."
		./bin/plink --bfile ${outstem}P2 --bim ${outstem}.bim --fam ${outstem}P.fam \
			--make-grm-bin --out ${outstem}P2
		nlines=$(< "${outstem}P.fam" wc -l)
		adjustGRM ${outstem}P2 $nlines 2

		echo "Producing paternal-genotype GRM..."
		./bin/plink --vcf ${infile} --keep ${outstem}P.tokeep \
			--make-grm-bin --out ${outstem}PG
		adjustGRM ${outstem}PG $nlines $divc
	fi

	echo "All GRMs produced."
fi

if [ "$esth2" = true ];
then
	echo "Estimating h2 from fetal-genotype GRM..."
	./bin/gcta64 --reml --grm ${outstem}FG --pheno ${inpheno} --out ${finalout}
	parseRes ${finalout}.log ${1}FG

	echo "Estimating h2 from transmitted haplotype GRMs..."
	echo -e "${outstem}M1\n${outstem}P1" > ${grmlist}
	./bin/gcta64 --reml --mgrm ${grmlist} --pheno ${inpheno} --out ${finalout}
	parseRes ${finalout}.log ${1}M1P1
	
	if [ "$haveMothers" = true ];
	then
		echo "Estimating h2 from maternal-genotype GRM..."
		./bin/gcta64 --reml --grm ${outstem}MG --pheno ${inpheno} --out ${finalout}
		parseRes ${finalout}.log ${1}MG

		echo "Estimating h2 from maternal haplotype GRMs..."
		echo -e "${outstem}M1\n${outstem}M2" > ${grmlist}
		./bin/gcta64 --reml --mgrm ${grmlist} --pheno ${inpheno} --out ${finalout}
		parseRes ${finalout}.log ${1}M1M2

		echo "Estimating h2 from maternal and fetal haplotype GRMs..."
		echo "${outstem}P1" >> ${grmlist}
		./bin/gcta64 --reml --mgrm ${grmlist} --pheno ${inpheno} --out ${finalout}
		parseRes ${finalout}.log ${1}M1M2P1
	fi
	if [ "$haveFathers" = true ];
	then
		echo "Estimating h2 from paternal-genotype GRM..."
		./bin/gcta64 --reml --grm ${outstem}PG --pheno ${inpheno} --out ${finalout}
		parseRes ${finalout}.log ${1}PG

		echo "Estimating h2 from saternal haplotype GRMs..."
		echo -e "${outstem}P1\n${outstem}P2" > ${grmlist}
		./bin/gcta64 --reml --mgrm ${grmlist} --pheno ${inpheno} --out ${finalout}
		parseRes ${finalout}.log ${1}P1P2

		echo "Estimating h2 from paternal and fetal haplotype GRMs..."
		echo "${outstem}M1" >> ${grmlist}
		./bin/gcta64 --reml --mgrm ${grmlist} --pheno ${inpheno} --out ${finalout}
		parseRes ${finalout}.log ${1}P1P2M1

		if [ "$haveMothers" = true ];
		then
			echo "Estimating h2 from all 4 haplotype GRMs..."
			echo "${outstem}M2" >> ${grmlist}
			./bin/gcta64 --reml --mgrm ${grmlist} --pheno ${inpheno} --out ${finalout}
			parseRes ${finalout}.log ${1}M1M2P1P2
		fi
	fi

	echo "All requested GCTA runs complete."
fi
