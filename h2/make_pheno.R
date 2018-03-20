#!/usr/bin/Rscript

## This script simulates VCF files for testing purposes.

##### PARAMETERS #####
nsnps = 90000
ninds = 500
h2 = 0.5

# doesn't mean anything
chr = 22

# standard meta-info for VCF header
header = readLines("vcfheader.txt")
setwd("~/Documents/haplotypes/")

# allele effects
betas = rnorm(nsnps)


##### SIMULATION ######
vcf = data.frame(CHROM=chr, POS=1:nsnps,
                 ID=paste0("snp", 1:nsnps),
                 REF="A", ALT="B", QUAL=1, FILTER="PASS", INFO="INFO=1", FORMAT="GT")
gtm = matrix(0, nsnps, ninds)

# main simulation loop
for(i in 1:ninds){
  gts = rbinom(nsnps, 2, 0.5)
  
  gtm[,i] = gts
  gts[gts==0] = "0|0"
  gts[gts==1] = ifelse(runif(sum(gts==1))>0.5, "0|1", "1|0")
  gts[gts==2] = "1|1"
  
  vcf = cbind(vcf, iid=gts)
}

# write VCF
writeLines(header, "test.vcf")
cat("#", file="test.vcf", append=T)
write.table(vcf, "test.vcf", T, col.names = T, row.names = F, quote=F, sep="\t")

# bgzip and tabix the VCF
system("bgzip -f test.vcf; tabix test.vcf.gz")


## pheno file
pheno = data.frame(fid = 1:ninds, iid=colnames(vcf)[-c(1:9)])
# create the genetic component of phenotype
pheno$pheno = t(gtm) %*% betas

# add environmental noise to produce the desired h2
pheno$pheno = pheno$pheno + rnorm(ninds, 0, sqrt((1 - h2)/h2) * sd(pheno$pheno))

write.table(pheno, "test.pheno", col.names = T, row.names = F, quote=F, sep="\t")

