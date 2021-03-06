# HAPLOTYPES
Collection of scripts for studies based on haplotyped trios/duos, in collaboration with CCHMC.
Basic package implements the pipeline:  
 phased VCFs -> `split_haplotypes + plink + gcta` -> heritability estimates  
Other scripts help prepare & check the phased VCFs, pedigree files, and provide simulation testing.

## Installation of basic package
1. Download **GCTA** and **PLINK 1.9**. Linux x64 binaries are provided in `bin/`, for other systems look up their provider webpages.  
2. Download `bin/split_haplotypes` and `bin/divide_grm` (for Linux x64 systems), or the source files `split_haplotypes.c` and `divide_grm.c`.  

## Usage
### Input and output formats
Main input is a *nice-phase* VCF file. It is a standard VCF with phase information (`|`-separated genotypes in `GT` field), but with alleles aligned so that the left haplotype on parents is the transmitted one, and left haplotype in child is the maternally inherited one (i.e. M|P or T|U). See `input-output_file_formats.pdf` for more details.  
Pedigree information is provided in a file with at least 3 columns, indicating child ID, father ID, and mother ID, respectively. **Duos must have 0 for the missing parent ID.** Any subsequent columns will be ignored. This format is equivalent to columns 2, 3, and 4 of PLINK's .fam files, after removing parents' rows. **All samples in the VCF file must be present in the pedigree file.**  

### Running
For basic usage, change parameters in `split_haplotype_launcher.sh` (input and output file locations) and set the `convert`, `makegrm`, and `esth2` flags to `true`. Set the `haveMothers` or `haveFathers` to `false` if you have only duos of some kind. Set `plinkfam` to `true` if you have pedigree information in a PLINK-style .fam file (correct columns will be extracted, and parents' rows removed.)


### Running only `split_haplotypes`
Input (unzipped) VCF file is streamed from STDIN, so recommend using `zcat` on gzipped files. Arguments are: pedigree input file, log output file, stem for other output files.  
So `zcat input.vcf.gz | ./bin/split_haplotypes input.fam output.log out` will produce:
1. `output.log` - log file
2. `outM1.bed` - maternal transmitted haplotypes
3. `outM2.bed` - maternal untransmitted haplotypes
4. `outP1.bed` - paternal transmitted haplotypes
5. `outP2.bed` - paternal untransmitted haplotypes
6. `out.bim` - one marker map
7. `outF.fam`, `outM.fam`, `outP.fam` - list of children corresponding to `outM1.bed` and `outP1.bed`, `outM2.bed`, `outP2.bed`, respectively
8. `outM.tokeep`, `outP.tokeep` - list of mothers and fathers (to keep from VCF)

### Downstream usage
Produced files can be directly passed to PLINK/GCTA, for example:
```
plink --bfile outM1 --bim out.bim --fam outF.fam --keep-allele-order --recodeA --out recoded
```

Note that `0/2` haplotype encoding makes PLINK overestimate GRM by a factor of 2 (i.e. expected value on the diagonal is 2), so we provide a simple `divide_grm.c` script to adjust a binary haplotype-based GRM:
```
# Arguments: inputFile outputFile nSamples divisor
./divide_grm M1.grm.bin M1adjusted.grm.bin 2000 2
```

### Advanced usage
// TODO

## Example cases
Some example files of `s` markers and `i` individuals are provided in `testcases/example_sXXX_iXXX.vcf.gz`. More can be generated by `make_pheno.R` (beware - not optimized and hogs RAM). `benchmark_testcases.sh` can help plan the runtime for large datasets.


## Current features
- speed
- PLINK- and GCTA-compatible output (.bed/.bim/.fam)
- supports mixed trio, duo, and singleton input
- supports full- and half-siblings
- proper treatment multi-allelic sites
- proper treatment of missing input genotypes

## Known limitations
- produced .bed files are aligned to ref/alt alleles, not to major/minor
- child ID must be non-missing (singletons must be coded as children)
- all samples in the VCF must be present in the input fam file at least once
- multi-generational trios will produce errors
- auto-generated .fam and .tokeep files might be wrong under weird pedigrees - recommend checking them
- half calls treated as missing
- no support for X chromosomes or other non-diploid genotypes


## Input and output formats
Described in the **input-output_file_formats.pdf**.


## Other scripts
### Phasing checker
The script **check_phasing.R** checks haplotype alignment inside trios after SHAPEIT phasing.

### Haplotype flipper
The script **flip_haplotypes.c** assigns SHAPEITv2-phased haplotypes to transmitted/untransmitted. Input must be trios, phased with duoHMM on. Then it is assumed that fetal haplotypes are paternal-left, maternal-right, and parental haplotypes are assigned transmission/non-transmission based on that. Input parental phase information is ignored.

### GRM divider
`./divide_grm` can also be used to adjust for confounding, by replacing 2 with the true variance of genotypes, i.e. `1+Fstat`.
