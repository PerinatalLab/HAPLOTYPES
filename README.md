# HAPLOTYPES
studies based on haplotyped trios/duos, in collaboration with CCHMC

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


To include later: h1/h2/h3 haplotype splitting scripts, parent-of-origin analyses, genetic risk scores analyses, mendelian randomisation (with untransmitted haplotype).
