# This script checks haplotype alignment from phased files:
# 1. SHAPEITv2 with 1000G and pedigree info (--duoHMM)
# 2. same, but after Sanger imputation on HRC (w/o phasing)


### INIT
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)


### LOAD DATA
# From each dataset, extract one trio and one chromosome
# (trio selected so that all are on batch 24 and QC)
trio = c(fet="9981196059_R12C01", pat="9998714025_R08C02", mat="9998714025_R07C02")
bcfoptions = paste("-s", paste(trio, collapse=","), "-f '%POS %REF %ALT[ %GT]\n'")

in_eagle = "/mnt/HUNT/tsdphased/imputed_m24/22.vcf.gz"
in_shapeit = "/mnt/HUNT/tsdphased/m24chr22.phased.vcf.gz"
out_eagle = "/mnt/HUNT/haplotypes/test_eagle.vcf"
out_shapeit = "/mnt/HUNT/haplotypes/test_shapeit.vcf"

## extract the trio
# (sample order is same as in -s flag)
command = paste("bcftools query", in_eagle, bcfoptions, "-o", out_eagle)
print(command)
system(command)
command = paste("bcftools query", in_shapeit, bcfoptions, "-o", out_shapeit)
print(command)
system(command)

## load and split into haplotypes
eagle = read.table(out_eagle)
shapeit = read.table(out_shapeit)
colnames(eagle) = colnames(shapeit) = c("POS", "REF", "ALT", names(trio))
eagle = separate(eagle, fet, c("fetA", "fetB"), sep="\\|", convert=T) %>% 
	separate(pat, c("patA", "patB"), sep="\\|", convert=T) %>%
	separate(mat, c("matA", "matB"), sep="\\|", convert=T)
shapeit = separate(shapeit, fet, c("fetA", "fetB"), sep="\\|", convert=T) %>% 
	separate(pat, c("patA", "patB"), sep="\\|", convert=T) %>%
	separate(mat, c("matA", "matB"), sep="\\|", convert=T)

## flip (EAGLE dataset) alleles if needed
toflip = inner_join(eagle, shapeit, by=c("POS", "REF"="ALT", "ALT"="REF"))
toflip = paste(toflip$POS, toflip$REF, toflip$ALT, sep="_")
toflip = which(paste(eagle$POS, eagle$REF, eagle$ALT, sep="_") %in% toflip)
eagle[toflip, c("REF", "ALT")] = eagle[toflip, c("ALT", "REF")]
eagle[toflip, 4:9] = 1 - eagle[toflip, 4:9]

## identify transmission and merge
# (pat/mat identify the source of A hap)
phaseTable = expand.grid(fetA=0:1, fetB=0:1, pat=0:2, mat=0:2)
phaseTable$conclusion = "ME"
phaseTable$conclusion[c(1, 5, 13, 17:20, 24, 32, 36)] = "NI"
phaseTable$conclusion[c(7, 11, 14, 23, 26, 30)] = "mat"
phaseTable$conclusion[c(6, 10, 15, 22, 27, 31)] = "pat"

# for duos:
phaseTableDuoF = expand.grid(fetA=0:1, fetB=0:1, matA=0:1, matB=0:1)
phaseTableDuoF$conclusion = "NI"
phaseTableDuoF$conclusion[c(4, 13)] = "ME"
phaseTableDuoF$conclusion[c(3, 14)] = "A"
phaseTableDuoF$conclusion[c(2, 15)] = "B"

phaseTableDuoM = expand.grid(fetA=0:1, fetB=0:1, matA=0:1, matB=0:1)
phaseTableDuoM$conclusion = "NI"
phaseTableDuoM$conclusion[c(4, 13)] = "ME"
phaseTableDuoM$conclusion[c(3, 14)] = "A"
phaseTableDuoM$conclusion[c(2, 15)] = "B"


eagle = mutate(eagle, pat=patA+patB, mat=matA+matB) %>%
	inner_join(phaseTable, by=c("pat", "mat", "fetA", "fetB"))
shapeit = mutate(shapeit, pat=patA+patB, mat=matA+matB) %>%
	inner_join(phaseTable, by=c("pat", "mat", "fetA", "fetB"))
merged = full_join(eagle, shapeit, by=c("POS", "REF", "ALT"), suffix=c("e", "s"))


### PLOT
merged$conclusione = factor(merged$conclusione, c("ME", "mat", "pat", "NI"))
merged$conclusions = factor(merged$conclusions, c("ME", "mat", "pat", "NI"))
# as expected, most SNPs are non-informative (e.g. rare 00s):
filter(merged, !is.na(conclusions)) %>%
	ggplot(aes(x=POS)) +
	geom_point(aes(y="only phased", col=conclusions), pch="|") +
	geom_point(aes(y="imputed", col=conclusione), pch="|") +
	ylab(NULL) + ggtitle("origin of fetal haplotype A") +
	theme_bw()

# most of the rest seem to have paternal hap in A:
filter(merged, !is.na(conclusions), conclusions!="NI") %>%
	ggplot(aes(x=POS)) +
	geom_point(aes(y="only phased", col=conclusions), pch="|") +
	geom_point(aes(y="imputed", col=conclusione), pch="|") +
	ylab(NULL) + ggtitle("origin of fetal haplotype A") +
	theme_bw()
table(merged$conclusions)
table(merged$conclusione)


### RETEST MORE TRIOS
# get some more trios and diverse pedigrees to test this
fam = read.table("/mnt/HUNT/erc-genotypes-results/moba24/raw-data/raw-updated/raw-updated.fam")
flags = read.table("~/data/geno/harvest-aux/harvest-flag-list.txt", h=T)
flags = filter(flags, genotypesOK, phenotypesOK)

# drop duplicates and QC fails:
# (note that we ignore CORE because it includes only one sib)
fam2 = filter(fam, grepl("^9", fam$V2)) %>%
	semi_join(flags, by=c("V2"="IID"))

sort(table(fam2$V1), decreasing = T)[1:10]
# nice 3-kid family:
filter(fam2, V1==6247)
trios = filter(fam2, V1==6247, V3!=0) %>% select(V2:V4)

# 2-kid family:
filter(fam2, V1==195)
trios = filter(fam2, V1==195, V3!=0) %>% select(V2:V4) %>%
	bind_rows(trios, .)

# some random trios:
sort(table(fam2$V1), decreasing = T)[c(201, 210, 334)]
filter(fam2, V1==393)
trios = filter(fam2, V1==393, V3!=0) %>% select(V2:V4) %>%
	bind_rows(trios, .)
filter(fam2, V1==417)
trios = filter(fam2, V1==417, V3!=0) %>% select(V2:V4) %>%
	bind_rows(trios, .)
filter(fam2, V1==1078)
trios = filter(fam2, V1==1078, V3!=0) %>% select(V2:V4) %>%
	bind_rows(trios, .)

# suspicious:
filter(fam2, V1==4357)
trios = filter(fam2, V1==4357, V3!=0) %>% select(V2:V4)

colnames(trios) = c("fet", "pat", "mat")


# just a function wrapper around same script,
# using only genotyped SNPs now
in_shapeit = "/mnt/HUNT/tsdphased/m24chr22.phased.vcf.gz"
out_shapeit = "/mnt/HUNT/haplotypes/test_shapeit.vcf"

in_shapeit = "~/data/geno/imputed/fresh/22.vcf.gz"
out_shapeit = "~/Documents/haplotypes/mobatest.vcf"

processTrio = function(trio){
	print(paste("Working on trio", trio))
	bcfoptions = paste("-s", paste(trio, collapse=","), "-f '%POS %REF %ALT[ %GT]\n'")
	
	## extract the trio
	# (sample order is same as in -s flag)
	command = paste("bcftools query", in_shapeit, bcfoptions, "-o", out_shapeit)
	print(command)
	system(command)
	
	## load and split into haplotypes
	shapeit = read.table(out_shapeit)
	colnames(shapeit) = c("POS", "REF", "ALT", "fet", "pat", "mat")
	shapeit = separate(shapeit, fet, c("fetA", "fetB"), sep="\\|", convert=T) %>% 
		separate(pat, c("patA", "patB"), sep="\\|", convert=T) %>%
		separate(mat, c("matA", "matB"), sep="\\|", convert=T)
	
	## identify transmission and merge
	shapeit = mutate(shapeit, pat=patA+patB, mat=matA+matB) %>%
		inner_join(phaseTable, by=c("pat", "mat", "fetA", "fetB"))
	shapeit$conclusion = factor(shapeit$conclusion, c("ME", "mat", "pat", "NI"))
	
	return(shapeit)
}

processDuo = function(duo){
	print(paste("Working on duo", duo))
	bcfoptions = paste("-s", paste(duo, collapse=","), "-f '%POS %REF %ALT[ %GT]\n'")
	
	## extract the duo
	# (sample order is same as in -s flag)
	command = paste("bcftools query", in_shapeit, bcfoptions, "-o", out_shapeit)
	print(command)
	system(command)
	
	## load and split into haplotypes
	shapeit = read.table(out_shapeit)
	colnames(shapeit) = c("POS", "REF", "ALT", "fet", "mat")
	shapeit = separate(shapeit, fet, c("fetA", "fetB"), sep="\\|", convert=T) %>% 
		separate(mat, c("matA", "matB"), sep="\\|", convert=T)
	
	## identify transmission and merge
	shapeit = mutate(shapeit, mat=matA+matB) %>%
		inner_join(phaseTableDuo, by=c("mat", "fetA", "fetB"))
	shapeit$conclusion = factor(shapeit$conclusion, c("ME", "mat", "NI"))
	
	return(shapeit)
}


## main run
toplot = apply(trios, 1, processTrio)
names(toplot) = paste0("f", c("1a", "1b", "1c", "2a", "2b", "3", "4", "5"))
toplot = bind_rows(toplot, .id="fid")
# table(toplot$fid)


### RESULTS
# no MEs, fetA is always paternal:
table(toplot$conclusion, toplot$fid)


### ADDITIONAL CHECKS
## check which haplotype is transmitted:
toplot$patT = toplot$matT = factor("NI", levels = c("A", "B", "NI"))
toplot$patT[toplot$patA==toplot$fetA & toplot$patB!=toplot$fetA] = "A"
toplot$patT[toplot$patA!=toplot$fetA & toplot$patB==toplot$fetA] = "B"
toplot$matT[toplot$matA==toplot$fetB & toplot$matB!=toplot$fetB] = "A"
toplot$matT[toplot$matA!=toplot$fetB & toplot$matB==toplot$fetB] = "B"

## plot
pPat = filter(toplot, patT!="NI") %>%
	ggplot(aes(x=POS)) +
	geom_point(aes(y=fid, col=patT), pch="|") +
	scale_color_manual(values=c("red1", "mediumblue")) +
	ylab("father from family") + 
	theme_bw()
pMat = filter(toplot, matT!="NI") %>%
	ggplot(aes(x=POS)) +
	geom_point(aes(y=fid, col=matT), pch="|") +
	scale_color_manual(values=c("red1", "mediumblue")) +
	ylab("mother from family") + 
	theme_bw()
plot_grid(pPat, pMat, nrow = 2)


## check ALL trios
# keep only if both parents are in m24 + qc pass:
fam3 = filter(fam2, V3 %in% V2, V4 %in% V2)
famsizes = sort(table(fam3$V1), decreasing = T)
alltrios = filter(fam3, V1 %in% names(famsizes)[famsizes==1]) %>% select(V2:V4)
alltrios = apply(alltrios, 1, processTrio)
alltrios = bind_rows(alltrios, .id="fid")
# total errors: less than 1 in 6000 SNPs
table(alltrios$conclusion)
# (some trios are phased differently because children 
# that were bad duplicates were phased as singletons)
# should use /mnt/HUNT/erc-genotypes-results/postqc/to_imputation/m24-ready-for-imputation.fam
alltrios2 = alltrios
alltrios2 = filter(alltrios2, fid!=1001)
# approx. errors per genome: 5.7
table(alltrios2$conclusion)
17010/2354 * 60

## see how many trios are split across m12-m24
fam12a = read.table("/mnt/HUNT/erc-genotypes-results/moba12good-reclustered/raw-data/raw-updated/raw-updated.fam")
fam12b = read.table("/mnt/HUNT/erc-genotypes-results/moba12bad-reclustered/raw-data/raw-updated/raw-updated.fam")
# all fam12 mothers/fathers are in m12
table(fam12a$V3[fam12a$V3!=0] %in% fam$V2)
table(fam12a$V4[fam12a$V4!=0] %in% fam$V2)

table(fam12b$V3[fam12b$V3!=0] %in% fam$V2)
table(fam12b$V4[fam12b$V4!=0] %in% fam$V2)

# but 147 kids from m24 have their mothers/fathers in m12
table(fam$V3[fam$V3!=0] %in% fam$V2 & fam$V4[fam$V3!=0] %in% fam$V2)



####
#Tetyana's checks
test = read.table("/mnt/HUNT/temp_check.vcf",comment.char = "", skip=9, h=T)
test24 = read.table("/mnt/HUNT/temp_check.vcf",comment.char = "", skip=9, h=T)
testgenos = test24[,10:ncol(test24)]
testgenos = t(testgenos)
ttfet = paste(c(1,1,1,1,1,0,0,0,0,0, 1,1,1,1,1,1,1), c(1,1,1,1,1, 1,1,1,0,1,1,0,1,1,1,0,0), sep="|")
uneqs = apply(testgenos, 1, function(x) sum(x!=ttfet))
table(uneqs)
which(uneqs==1)
testgenos[which(uneqs==2),]
testgenos = testgenos[rowSums(testgenos[,c(1:5, 11, 13:15)]=="1|1")==9,]
testgenos = testgenos[rowSums(testgenos[,c(6:8, 10, 12, 16, 17)]=="0|1" |
					  	testgenos[,c(6:8, 10, 12, 16, 17)]=="1|0")==7,]
testgenos = testgenos[testgenos[,9]=="0|0",]
dim(testgenos)
testgenos
