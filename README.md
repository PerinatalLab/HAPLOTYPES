# HAPLOTYPES
studies based on haplotyped trios/duos, in collaboration with CCHMC

### Scripts
The script check_phasing.R checks haplotype alignment inside trios after SHAPEIT phasing.

The script flip_haplotypes.c assigns SHAPEITv2-phased haplotypes to transmitted/untransmitted. Input must be trios, phased with duoHMM on. Then it is assumed that fetal haplotypes are paternal-left, maternal-right, and parental haplotypes are assigned transmission/non-transmission based on that. Input parental phase information is ignored.
Second input must be file identifying family trios, tab-separated:
`fetalID  paternalID maternalID`
Output order is
`snp_info_columns maternalT|maternalU paternalT|paternalU fetalM|fetalP`


To include later: h1/h2/h3 haplotype splitting scripts, parent-of-origin analyses, genetic risk scores analyses, mendelian randomisation (with untransmitted haplotype).
