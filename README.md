# HAPLOTYPES
Collection of scripts for studies based on haplotyped trios/duos, in collaboration with CCHMC.
Basic package is phased VCFs -> `split_haplotypes + plink + gcta` -> heritability estimates. Other scripts help prepare & check the phased VCFs, and provide simulation testing.

## Installation of basic package
1. Download **GCTA** and **PLINK 1.9**. Linux x64 binaries are provided in `bin/`, for other systems look up their provider webpages.  
2. Download `bin/split_haplotypes` (for Linux x64 systems), or compile it yourself from `split_haplotypes.c` as `gcc split_haplotypes.c -o split_haplotypes`.  

## Usage
### Input and output formats
Main input is a *nice-phase* VCF file. It is a standard VCF with phase information (`|`-separated genotypes in `GT` field), but with alleles aligned so that the left haplotype on parents is the transmitted one, and left haplotype in child is the maternally inherited one (i.e. M|P or T|U). See `input-output_file_formats.pdf` for more details.
Pedigree information is provided in a file with at least 3 columns, indicating child ID, father ID, and mother ID, respectively. Any subsequent columns will be ignored. This format is equivalent to columns 2, 3, and 4 of PLINK's .fam files, after removing parents' rows.

### Running
`zcat input.vcf.gz | ./split_haplotypes input.fam output.log out`
will produce:
1. `output.log` - log file
2. `outM1.bed` - maternal transmitted haplotypes
3. `outM2.bed` - maternal untransmitted haplotypes
4. `outP1.bed` - paternal transmitted haplotypes
5. `outP2.bed` - paternal untransmitted haplotypes
6. `out.bim` - one marker map
7. `outF.fam` - list of children (corresponds to `outM1.bed` and `outP1.bed`)
8. `outM.fam` - list of mothers (corresponds to `outM2.bed`)
9. `outP.fam` - list of fathers (corresponds to `outP2.bed`)

For example, downstream usage of the output could be:
```
plink --bfile outM1 --bim out.bim --fam outF.fam --keep-allele-order --recodeA --out recoded
```

### Advanced usage
// TODO

## Current features
- supports mixed trio, duo, and singleton input
- supports full- and half-siblings
- proper treatment multi-allelic sites
- proper treatment of missing input genotypes

## Known limitations
- produced .bed files are aligned to ref/alt alleles, not to major/minor
- child ID must be non-missing (singletons must be coded as children)
- half calls treated as missing
- multi-generational trios will produce errors
- no support for X chromosomes or other non-diploid genotypes


## Input and output formats
Described in the **input-output_file_formats.pdf**.

## Scripts
### Phasing checker
The script **check_phasing.R** checks haplotype alignment inside trios after SHAPEIT phasing.

### Haplotype flipper
The script **flip_haplotypes.c** assigns SHAPEITv2-phased haplotypes to transmitted/untransmitted. Input must be trios, phased with duoHMM on. Then it is assumed that fetal haplotypes are paternal-left, maternal-right, and parental haplotypes are assigned transmission/non-transmission based on that. Input parental phase information is ignored.

Current features:
- supports mixed trio and duo input (resulting in half-calls)
- safe against missing input genotypes
- detects Mendelian errors
- (almost) VCF-compatible structure of main output
- runtime on the order of 1 millisecond per marker


### Haplotype splitter
The script **split_haplotypes.c** takes properly phased VCF files and produces .bed/.bim/.fam filesets for each haplotype.
'Properly phased' here means that fetal haplotypes are maternal-left, paternal-right, and both parents are aligned transmitted-left, untransmitted-right. Provided an input file with childID|fatherID|motherID columns, this script splits fetal genotypes into M1 and P1 files ("transmitted"), and takes untransmitted haplotypes (i.e. right columns) of parents into M2 and P2 files.

Current features:
- PLINK- and GCTA-compatible output (.bed/.bim/.fam)
- produces one .bim, three .fam, and four .bed files
- safe against missing and weird input genotypes
- allows mixed singleton, maternal and paternal data 
- runtime on the order of 1 millisecond per marker

Known limitations:
- all families must include a fetus, i.e. singletons must be provided as "childID 0 0"
- multigenerational families are not supported
- siblings are only allowed if the parents' columns are also repeated in the VCF (so effectively not allowed)

To include later: h1/h2/h3 haplotype splitting scripts, parent-of-origin analyses, genetic risk scores analyses, mendelian randomisation (with untransmitted haplotype).

