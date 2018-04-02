#!/usr/bin/Rscript

## This script simulates VCF files for testing purposes.

##### PARAMETERS #####
nsnps = 3000
ninds = 6000
h2 = 0.5

# fetally or maternally determined phenotype?
# geneff="f"
# provide number to fix, or other value to generate runif()
# setMaf = 0.5
# p(heterozygote phase = 1|0)
# phaseBias = 0.5
# fraction of heterozygotes to convert to homozygotes
# corrBias = 0

# doesn't mean anything
chr = 22

# standard meta-info for VCF header
header = readLines("~/Documents/gitrep/HAPLOTYPES/vcfheader.txt")
# setwd("~/Documents/gitrep/HAPLOTYPES/testcases/")
setwd(workdir)

# allele effects
betas = rnorm(nsnps)


##### SIMULATION ######
vcf = data.frame(CHROM=chr, POS=1:nsnps,
                 ID=paste0("snp", 1:nsnps),
                 REF="A", ALT="B", QUAL=1, FILTER="PASS", INFO="INFO=1", FORMAT="GT")

if(is.numeric(setMaf)){
  mafs = rep(setMaf, nsnps)
} else {
  mafs = runif(nsnps)  
}

# simulate dosage matrix
dsm = t(sapply(mafs, function(m) rbinom(ninds, 2, m)))
# introduce confounding (reduces heterozygosity)
het = dsm==1
toHom = runif(dsm)
dsm[het & toHom < corrBias/2] = 0
dsm[het & toHom > 1-corrBias/2] = 2

# change into genotypes
gtm = matrix("0|0", nsnps, ninds)
gtm[dsm==1] = sample(c("0|1", "1|0"), sum(dsm==1), replace=T, prob=c(1-phaseBias, phaseBias))
gtm[dsm==2] = "1|1"
vcf = cbind(vcf, gtm)
colnames(vcf)[-c(1:9)] = paste0("iid", 1:ninds)


# write VCF
filestem = paste0("example_s", nsnps, "_i", ninds)
writeLines(header, paste0(filestem, ".vcf"))
cat("#", file=paste0(filestem, ".vcf"), append=T)
write.table(vcf, paste0(filestem, ".vcf"), T, col.names = T, row.names = F, quote=F, sep="\t")

# bgzip and tabix the VCF
system(paste0("bgzip -f ", filestem, ".vcf; tabix ", filestem, ".vcf.gz"))

## fam file
fam = as.data.frame(matrix(paste0("iid", 1:ninds), ncol=3))
write.table(fam, paste0(filestem, ".fam"), col.names=F, row.names=F, quote=F, sep="\t")


# pheno file
pheno = data.frame(fid = colnames(vcf)[-c(1:9)])
pheno$iid = pheno$fid
# create the genetic component of phenotype
pheno$pheno = t(dsm) %*% betas

if(geneff=="f"){
	pheno = pheno[pheno$iid %in% fam[,1],]
} else if (geneff=="m"){
	pheno = pheno[pheno$iid %in% fam[,2],]
	pheno$iid = fam[,1]
	pheno$fid = pheno$iid
}

# add environmental noise to produce the desired h2
print(var(pheno$pheno))
pheno$pheno = pheno$pheno + rnorm(nrow(pheno), 0, sqrt((1 - h2)/h2) * sd(pheno$pheno))
print(var(pheno$pheno))
 
write.table(pheno, paste0(filestem, ".pheno"), col.names = F, row.names = F, quote=F, sep="\t")

