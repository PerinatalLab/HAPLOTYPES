#!/bin/bash
set -e
cd ~/Documents/gitrep/HAPLOTYPES/

for i in {1..10}
do
	echo "#### round $i/10 ####"
	echo "#### FETAL EFFECTS ####"
	echo "maf=0.5, phase=0.7, corrBias=0, geneff=f"
	Rscript -e 'setMaf=0.5; phaseBias=0.7; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/tmp3/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf50phase70corr0effF ~/Documents/haplotypes/tmp3/

	echo "maf=0.3, phase=0.7, corrBias=0, geneff=f"
	Rscript -e 'setMaf=0.3; phaseBias=0.7; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/tmp3/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf30phase70corr0effF ~/Documents/haplotypes/tmp3/

	echo "maf=0.1, phase=0.7, corrBias=0, geneff=f"
	Rscript -e 'setMaf=0.1; phaseBias=0.7; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/tmp3/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf10phase70corr0effF ~/Documents/haplotypes/tmp3/

	echo "maf=0.7, phase=0.7, corrBias=0, geneff=f"
	Rscript -e 'setMaf=0.7; phaseBias=0.7; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/tmp3/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf70phase70corr0effF ~/Documents/haplotypes/tmp3/

	echo "maf=0.9, phase=0.7, corrBias=0, geneff=f"
	Rscript -e 'setMaf=0.9; phaseBias=0.7; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/tmp3/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf90phase70corr0effF ~/Documents/haplotypes/tmp3/

	echo "maf=runif, phase=0.7, corrBias=0, geneff=f"
	Rscript -e 'setMaf=NULL; phaseBias=0.7; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/tmp3/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh mafUphase70corr0effF ~/Documents/haplotypes/tmp3/

	echo "#### MATERNAL EFFECTS ####"
	echo "maf=0.5, phase=0.7, corrBias=0, geneff=m"
	Rscript -e 'setMaf=0.5; phaseBias=0.7; corrBias=0; geneff="m"; workdir="~/Documents/haplotypes/tmp3/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf50phase70corr0effM ~/Documents/haplotypes/tmp3/

	echo "maf=0.3, phase=0.7, corrBias=0, geneff=m"
	Rscript -e 'setMaf=0.3; phaseBias=0.7; corrBias=0; geneff="m"; workdir="~/Documents/haplotypes/tmp3/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf30phase70corr0effM ~/Documents/haplotypes/tmp3/

	echo "maf=0.1, phase=0.7, corrBias=0, geneff=m"
	Rscript -e 'setMaf=0.1; phaseBias=0.7; corrBias=0; geneff="m"; workdir="~/Documents/haplotypes/tmp3/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf10phase70corr0effM ~/Documents/haplotypes/tmp3/

	echo "maf=0.7, phase=0.7, corrBias=0, geneff=m"
	Rscript -e 'setMaf=0.7; phaseBias=0.7; corrBias=0; geneff="m"; workdir="~/Documents/haplotypes/tmp3/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf70phase70corr0effM ~/Documents/haplotypes/tmp3/

	echo "maf=0.9, phase=0.7, corrBias=0, geneff=m"
	Rscript -e 'setMaf=0.9; phaseBias=0.7; corrBias=0; geneff="m"; workdir="~/Documents/haplotypes/tmp3/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf90phase70corr0effM ~/Documents/haplotypes/tmp3/

	echo "maf=runif, phase=0.7, corrBias=0, geneff=m"
	Rscript -e 'setMaf=NULL; phaseBias=0.7; corrBias=0; geneff="m"; workdir="~/Documents/haplotypes/tmp3/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh mafUphase70corr0effM ~/Documents/haplotypes/tmp3/
done
