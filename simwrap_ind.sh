#!/bin/bash
set -e
cd ~/Documents/gitrep/HAPLOTYPES/

for i in {7..10}
do
	echo "#### round $i/10 ####"
	echo "NO PHASING ERROR"
	echo "maf=0.5, phaseErr=0, phaseBias=0.5"
	Rscript -e 'setMaf=0.5; phaseErr=0; phaseBias=0.5; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/ind/"; source("make_pheno_ind.R")'
	./split_haplotypes_launcher.sh maf50pe0pb50effF ~/Documents/haplotypes/ind/

	echo "maf=0.3, phaseErr=0, phaseBias=0.5"
	Rscript -e 'setMaf=0.3; phaseErr=0; phaseBias=0.5; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/ind/"; source("make_pheno_ind.R")'
	./split_haplotypes_launcher.sh maf30pe0pb50effF ~/Documents/haplotypes/ind/

	echo "maf=0.1, phaseErr=0, phaseBias=0.5"
	Rscript -e 'setMaf=0.1; phaseErr=0; phaseBias=0.5; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/ind/"; source("make_pheno_ind.R")'
	./split_haplotypes_launcher.sh maf10pe0pb50effF ~/Documents/haplotypes/ind/

	echo "PHASING ERROR, BUT NO BIAS"
	echo "maf=0.5, phaseErr=1, phaseBias=0.5"
	Rscript -e 'setMaf=0.5; phaseErr=1; phaseBias=0.5; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/ind/"; source("make_pheno_ind.R")'
	./split_haplotypes_launcher.sh maf50pe100pb50effF ~/Documents/haplotypes/ind/

	echo "maf=0.3, phaseErr=1, phaseBias=0.5"
	Rscript -e 'setMaf=0.3; phaseErr=1; phaseBias=0.5; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/ind/"; source("make_pheno_ind.R")'
	./split_haplotypes_launcher.sh maf30pe100pb50effF ~/Documents/haplotypes/ind/

	echo "maf=0.1, phaseErr=1, phaseBias=0.5"
	Rscript -e 'setMaf=0.1; phaseErr=1; phaseBias=0.5; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/ind/"; source("make_pheno_ind.R")'
	./split_haplotypes_launcher.sh maf10pe100pb50effF ~/Documents/haplotypes/ind/

	echo "PHASING ERROR AND 30% BIAS"
	echo "maf=0.5, phaseErr=1, phaseBias=0.3"
	Rscript -e 'setMaf=0.5; phaseErr=1; phaseBias=0.3; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/ind/"; source("make_pheno_ind.R")'
	./split_haplotypes_launcher.sh maf50pe100pb30effF ~/Documents/haplotypes/ind/

	echo "maf=0.3, phaseErr=1, phaseBias=0.3"
	Rscript -e 'setMaf=0.3; phaseErr=1; phaseBias=0.3; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/ind/"; source("make_pheno_ind.R")'
	./split_haplotypes_launcher.sh maf30pe100pb30effF ~/Documents/haplotypes/ind/

	echo "maf=0.1, phaseErr=1, phaseBias=0.3"
	Rscript -e 'setMaf=0.1; phaseErr=1; phaseBias=0.3; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/ind/"; source("make_pheno_ind.R")'
	./split_haplotypes_launcher.sh maf10pe100pb30effF ~/Documents/haplotypes/ind/

	echo "PHASING ERROR AND 70% BIAS"
	echo "maf=0.5, phaseErr=1, phaseBias=0.7"
	Rscript -e 'setMaf=0.5; phaseErr=1; phaseBias=0.7; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/ind/"; source("make_pheno_ind.R")'
	./split_haplotypes_launcher.sh maf50pe100pb70effF ~/Documents/haplotypes/ind/

	echo "maf=0.3, phaseErr=1, phaseBias=0.7"
	Rscript -e 'setMaf=0.3; phaseErr=1; phaseBias=0.7; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/ind/"; source("make_pheno_ind.R")'
	./split_haplotypes_launcher.sh maf30pe100pb70effF ~/Documents/haplotypes/ind/

	echo "maf=0.1, phaseErr=1, phaseBias=0.7"
	Rscript -e 'setMaf=0.1; phaseErr=1; phaseBias=0.7; corrBias=0; geneff="f"; workdir="~/Documents/haplotypes/ind/"; source("make_pheno_ind.R")'
	./split_haplotypes_launcher.sh maf10pe100pb70effF ~/Documents/haplotypes/ind/
done
