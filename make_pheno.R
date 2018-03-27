#!/usr/bin/Rscript

## This script simulates VCF files for testing purposes.

##### PARAMETERS #####
nsnps = 3000
ninds = 6000
h2 = 0.5

# doesn't mean anything
chr = 22

# standard meta-info for VCF header
header = readLines("~/Documents/gitrep/HAPLOTYPES/vcfheader.txt")
setwd("~/Documents/gitrep/HAPLOTYPES/testcases/")

# allele effects
betas = rnorm(nsnps)


##### SIMULATION ######
vcf = data.frame(CHROM=chr, POS=1:nsnps,
                 ID=paste0("snp", 1:nsnps),
                 REF="A", ALT="B", QUAL=1, FILTER="PASS", INFO="INFO=1", FORMAT="GT")
mafs = runif(nsnps)

# simulate dosage matrix
dsm = t(sapply(mafs, function(m) rbinom(ninds, 2, m)))
# change into genotypes
gtm = matrix("0|0", nsnps, ninds)
gtm[dsm==1] = sample(c("0|1", "1|0"), sum(dsm==1), replace=T)
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

# add environmental noise to produce the desired h2
pheno$pheno = pheno$pheno + rnorm(ninds, 0, sqrt((1 - h2)/h2) * sd(pheno$pheno))

write.table(pheno, paste0(filestem, ".pheno"), col.names = F, row.names = F, quote=F, sep="\t")

