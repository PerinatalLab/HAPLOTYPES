#!/bin/bash

# This script runs the standard testcases to benchmark performance.

echo -e "\ntime for 300 SNPs and 2000 trios:" > testcases/benchmark.txt 
{ time zcat testcases/example_s300_i6000.vcf.gz | ./bin/split_haplotypes testcases/example_s300_i6000.fam ~/Documents/haplotypes/splitlog ~/Documents/haplotypes/split; } 2>> testcases/benchmark.txt
echo -e "\ntime for 1000 SNPs and 2000 trios:" >> testcases/benchmark.txt 
{ time zcat testcases/example_s1000_i6000.vcf.gz | ./bin/split_haplotypes testcases/example_s1000_i6000.fam ~/Documents/haplotypes/splitlog ~/Documents/haplotypes/split; } 2>> testcases/benchmark.txt
echo -e "\ntime for 3000 SNPs and 6000 trios:" >> testcases/benchmark.txt 
{ time zcat testcases/example_s3000_i18000.vcf.gz | ./bin/split_haplotypes testcases/example_s3000_i18000.fam ~/Documents/haplotypes/splitlog ~/Documents/haplotypes/split; } 2>> testcases/benchmark.txt
echo -e "\ntime for 3000 SNPs and 2000 trios:" >> testcases/benchmark.txt 
{ time zcat testcases/example_s3000_i6000.vcf.gz | ./bin/split_haplotypes testcases/example_s3000_i6000.fam ~/Documents/haplotypes/splitlog ~/Documents/haplotypes/split; } 2>> testcases/benchmark.txt
echo -e "\ntime for 10000 SNPs and 2000 trios:" >> testcases/benchmark.txt 
{ time zcat testcases/example_s10000_i6000.vcf.gz | ./bin/split_haplotypes testcases/example_s10000_i6000.fam ~/Documents/haplotypes/splitlog ~/Documents/haplotypes/split; } 2>> testcases/benchmark.txt
echo -e "\ntime for 30000 SNPs and 2000 trios:" >> testcases/benchmark.txt 
{ time zcat testcases/example_s30000_i6000.vcf.gz | ./bin/split_haplotypes testcases/example_s30000_i6000.fam ~/Documents/haplotypes/splitlog ~/Documents/haplotypes/split; } 2>> testcases/benchmark.txt
