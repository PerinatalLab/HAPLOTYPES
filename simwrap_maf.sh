#!/bin/bash
set -e
cd ~/Documents/gitrep/HAPLOTYPES/

for i in {1..10}
do
	echo "#### round $i/10 ####"
	echo "#### FETAL EFFECTS ####"
	echo "maf=0.5, phase=0.5, corrBias=0, geneff=f"
	Rscript -e 'setMaf=0.5; phaseBias=0.5; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/tmp/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf50phase50corr0effF ~/Documents/haplotypes/tmp/

	echo "maf=0.3, phase=0.5, corrBias=0, geneff=f"
	Rscript -e 'setMaf=0.3; phaseBias=0.5; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/tmp/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf30phase50corr0effF ~/Documents/haplotypes/tmp/

	echo "maf=0.1, phase=0.5, corrBias=0, geneff=f"
	Rscript -e 'setMaf=0.1; phaseBias=0.5; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/tmp/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf10phase50corr0effF ~/Documents/haplotypes/tmp/

	echo "maf=0.7, phase=0.5, corrBias=0, geneff=f"
	Rscript -e 'setMaf=0.7; phaseBias=0.5; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/tmp/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf70phase50corr0effF ~/Documents/haplotypes/tmp/

	echo "maf=0.9, phase=0.5, corrBias=0, geneff=f"
	Rscript -e 'setMaf=0.9; phaseBias=0.5; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/tmp/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf90phase50corr0effF ~/Documents/haplotypes/tmp/

	echo "maf=runif, phase=0.5, corrBias=0, geneff=f"
	Rscript -e 'setMaf=NULL; phaseBias=0.5; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/tmp/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh mafUphase50corr0effF ~/Documents/haplotypes/tmp/

	echo "#### MATERNAL EFFECTS ####"
	echo "maf=0.5, phase=0.5, corrBias=0, geneff=m"
	Rscript -e 'setMaf=0.5; phaseBias=0.5; corrBias=0; geneff="m"; workdir="~/Documents/haplotypes/tmp/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf50phase50corr0effM ~/Documents/haplotypes/tmp/

	echo "maf=0.3, phase=0.5, corrBias=0, geneff=m"
	Rscript -e 'setMaf=0.3; phaseBias=0.5; corrBias=0; geneff="m"; workdir="~/Documents/haplotypes/tmp/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf30phase50corr0effM ~/Documents/haplotypes/tmp/

	echo "maf=0.1, phase=0.5, corrBias=0, geneff=m"
	Rscript -e 'setMaf=0.1; phaseBias=0.5; corrBias=0; geneff="m"; workdir="~/Documents/haplotypes/tmp/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf10phase50corr0effM ~/Documents/haplotypes/tmp/

	echo "maf=0.7, phase=0.5, corrBias=0, geneff=m"
	Rscript -e 'setMaf=0.7; phaseBias=0.5; corrBias=0; geneff="m"; workdir="~/Documents/haplotypes/tmp/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf70phase50corr0effM ~/Documents/haplotypes/tmp/

	echo "maf=0.9, phase=0.5, corrBias=0, geneff=m"
	Rscript -e 'setMaf=0.9; phaseBias=0.5; corrBias=0; geneff="m"; workdir="~/Documents/haplotypes/tmp/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh maf90phase50corr0effM ~/Documents/haplotypes/tmp/

	echo "maf=runif, phase=0.5, corrBias=0, geneff=m"
	Rscript -e 'setMaf=NULL; phaseBias=0.5; corrBias=0; geneff="m"; workdir="~/Documents/haplotypes/tmp/"; source("make_pheno.R")'
	./split_haplotypes_launcher.sh mafUphase50corr0effM ~/Documents/haplotypes/tmp/
done
