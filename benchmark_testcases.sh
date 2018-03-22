#!/bin/bash

# This script runs the standard testcases to benchmark performance.

echo "time for 300 SNPs and 2000 trios:\n" > testcases/benchmark.txt 
time -a -o testcases/benchmark.txt zcat testcases/example_s300_i6000.vcf.gz | ./bin/split_haplotypes testcases/example_s300_i6000.fam ~/Documents/haplotypes/splitlog ~/Documents/haplotypes/split
echo "time for 1000 SNPs and 2000 trios:\n" >> testcases/benchmark.txt 
time -a -o testcases/benchmark.txt zcat testcases/example_s1000_i6000.vcf.gz | ./bin/split_haplotypes testcases/example_s1000_i6000.fam ~/Documents/haplotypes/splitlog ~/Documents/haplotypes/split
echo "time for 3000 SNPs and 2000 trios:\n" >> testcases/benchmark.txt 
time -a -o testcases/benchmark.txt zcat testcases/example_s3000_i6000.vcf.gz | ./bin/split_haplotypes testcases/example_s3000_i6000.fam ~/Documents/haplotypes/splitlog ~/Documents/haplotypes/split
echo "time for 10000 SNPs and 2000 trios:\n" >> testcases/benchmark.txt 
time -a -o testcases/benchmark.txt zcat testcases/example_s10000_i6000.vcf.gz | ./bin/split_haplotypes testcases/example_s10000_i6000.fam ~/Documents/haplotypes/splitlog ~/Documents/haplotypes/split
