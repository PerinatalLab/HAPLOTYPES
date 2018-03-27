#include <stdlib.h>
#include <stdio.h>

// This script is for adjusting PLINK's GRMs
// to account for having 0/2 haplotype encoding.
// I.e. it divides by 2.

// USAGE: ./divide_grm infile.grm.bin outfile.grm.bin nsamples

int main(int argc, char *argv[]){
	FILE *infile, *outfile;
	if(argc != 4){
		fprintf(stderr, "ERROR: script needs 3 arguments, %d supplied.\n", argc-1);
		exit(1);
	}
	infile = fopen(argv[1], "rb");
	outfile = fopen(argv[2], "wb");
	if(infile == NULL || outfile == NULL){
		fprintf(stderr, "ERROR: couldn't open input/output files\n");
		exit(1);
	}
	int n = atoi(argv[3]);

	// buffer of 250 people, 4 bytes each
	int b = 250;
	long unsigned int tot = n*(n+1)/2;
	float linebuf[b];
	
	// read n full buffers
	for(int i=0; i < tot/b; i++){
		fread(linebuf, b, 4, infile);

		for(int j=0; j<b; j++){
			linebuf[j] = linebuf[j]/2;
		}
		fwrite(linebuf, b, 4, outfile);
	}
	fread(linebuf, tot % b, 4, infile);
	for(int j=0; j < tot % b; j++){
		linebuf[j] = linebuf[j]/2;
	}
	fwrite(linebuf, tot % b, 4, outfile);

	fclose(infile);
	fclose(outfile);

	return(0);
}
