#!/usr/bin/Rscript

## This script simulates VCF files for testing purposes.

##### PARAMETERS #####
nsnps = 3000
ninds = 3000
h2 = 0.5

# fetally or maternally determined phenotype?
# geneff="f"
# provide number to fix, or other value to generate runif()
# setMaf = 0.5
# p(phasing error | double het)
phaseErr = 0
# p(phasing error towards 1|0)
phaseBias = 0.5
# fraction of heterozygotes to convert to homozygotes
# confRate = 0

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

# simulate dosage matrices for each paternal haplotype
dsmM1 = t(sapply(mafs, function(m) rbinom(ninds, 1, m)))
# dsmM2 = t(sapply(mafs, function(m) rbinom(ninds/3, 1, m)))
dsmP1 = t(sapply(mafs, function(m) rbinom(ninds, 1, m)))
# dsmP2 = t(sapply(mafs, function(m) rbinom(ninds/3, 1, m)))

# change into genotypes
gtmF = matrix(paste(dsmM1, dsmP1, sep="|"), nsnps, ninds)

# introduce phasing errors and bias
toflip = which(dsmM1==0 & dsmP1==1)
het01 = toflip
toflip = toflip[runif(toflip) < phaseErr * phaseBias]
gtmF[toflip] = "1|0"

toflip = which(dsmM1==1 & dsmP1==0)
het10 = toflip
toflip = toflip[runif(toflip) < phaseErr * (1-phaseBias)]
gtmF[toflip] = "0|1"

# introduce confounding (reduces heterozygosity)
# take 01 het. children, flip equal fractions to 00 or 11
flippingIndicator = runif(het01)
toflip = het01[flippingIndicator < confRate/2]
dsmP1[toflip] = dsmM1[toflip]
gtmF[toflip] = "0|0"
toflip = het01[flippingIndicator > 1-confRate/2]
dsmM1[toflip] = dsmP1[toflip]
gtmF[toflip] = "1|1"

# repeat with 10 het. children
flippingIndicator = runif(het10)
toflip = het10[flippingIndicator < confRate/2]
dsmP1[toflip] = dsmM1[toflip]
gtmF[toflip] = "1|1"
toflip = het10[flippingIndicator > 1-confRate/2]
dsmM1[toflip] = dsmP1[toflip]
gtmF[toflip] = "0|0"

# make the other GT matrices with M1,P1 already after confounding
# gtmM = matrix(paste(dsmM1, dsmM2, sep="|"), nsnps, ninds/3)
# gtmP = matrix(paste(dsmP1, dsmP2, sep="|"), nsnps, ninds/3)


# gtms are converted to VCF and NOT USED to calculate true phenotype
vcf = cbind(vcf, gtmF)
colnames(vcf)[-c(1:9)] = paste0("iid", 1:ninds)

# write VCF
filestem = paste0("example_s", nsnps, "_i", ninds)
writeLines(header, paste0(filestem, ".vcf"))
cat("#", file=paste0(filestem, ".vcf"), append=T)
write.table(vcf, paste0(filestem, ".vcf"), T, col.names = T, row.names = F, quote=F, sep="\t")

# bgzip and tabix the VCF
system(paste0("bgzip -f ", filestem, ".vcf; tabix ", filestem, ".vcf.gz"))

## fam file
fam = data.frame(iid=paste0("iid", 1:ninds), p=0, m=0)
write.table(fam, paste0(filestem, ".fam"), col.names=F, row.names=F, quote=F, sep="\t")


# pheno file
pheno = data.frame(fid = colnames(vcf)[10:(9+ninds)])
pheno$iid = pheno$fid

# create the genetic component of phenotype
# USING DSMs (untouched by phasing errors etc.)
if(geneff=="f"){
  pheno$pheno = t(dsmM1) %*% betas + t(dsmP1) %*% betas
} else if (geneff=="m"){
  pheno$pheno = t(dsmM1) %*% betas + t(dsmM2) %*% betas
} else if (geneff=="p"){
  pheno$pheno = t(dsmP1) %*% betas + t(dsmP2) %*% betas
}

# add environmental noise to produce the desired h2
print(var(pheno$pheno))
pheno$pheno = pheno$pheno + rnorm(nrow(pheno), 0, sqrt((1 - h2)/h2) * sd(pheno$pheno))
print(var(pheno$pheno))
 
write.table(pheno, paste0(filestem, ".pheno"), col.names = F, row.names = F, quote=F, sep="\t")

