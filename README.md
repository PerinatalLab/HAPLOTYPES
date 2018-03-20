# HAPLOTYPES
studies based on haplotyped trios/duos, in collaboration with CCHMC

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

Known limitations:
- multigenerational trios would cause wrong results  
- no support for X chromosomes or other non-diploid genotypes

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

