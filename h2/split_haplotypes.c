#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>

using namespace std;

#define version 11


// This script shall take in a family file, format compatible with above,
// allowing trios or duos with one 0 parent,
// and a vcf.gz or vcf file, with haplotypes aligned M|P T|U T|U,
// and produce three .tped files for easy input to PLINK.
// (LATER: make it output .bedsets directly?)


// USAGE: zcat onechrom.vcf.gz | ./split_haplotypes logfile outstem
// (create .tfam yourself for now)

// replace \n with \t to fix reading of last field:
void fixEndings(char *s){
	s[strlen(s)-1] = '\t';
}

int main(int argc, char *argv[]) {
	// INIT
	// arg check
	char *logFile, *outFile, *outFile1, *outFile2;
	if(argc != 3){
		fprintf(stderr, "ERROR: script needs 2 arguments, %d supplied.\n", argc-1);
		exit(1);
	} else {
		logFile = argv[1];
		outFile = argv[2];
		outFile1 = strdupa(outFile);
		outFile2 = strdupa(outFile);
		strcat(outFile, ".tped");
		strcat(outFile1, "1.tped");
		strcat(outFile2, "2.tped");
	}

	FILE *lfp, *of, *of1, *of2;

	lfp = fopen(logFile, "w");
	of = fopen(outFile, "w");
	of1 = fopen(outFile1, "w");
	of2 = fopen(outFile2, "w");
	if(lfp == NULL){
		fprintf(stderr, "ERROR: Can't create log file!\n");
		exit(1);
	}
	if(of == NULL || of1 == NULL || of2 == NULL){
		fprintf(stderr, "ERROR: Can't create one of output files!\n");
		exit(1);
	}
	fprintf(lfp, "Starting log for run.\n");
	fprintf(lfp, "Run command: %s %s %s\n", argv[0], argv[1], argv[2]);

	// READ
	// reading variables:
	char *line, *field, *linecp;
	size_t len = 0;
	ssize_t nread;

	// sample variables:
	int nsamples;
	int nlines = 0;
	int nmiss = 0;

	// OUTPUT PRINTING STARTS HERE
	// read VCF:
	fprintf(lfp, "Haplotype splitter v1\n");
	fprintf(lfp, "Reading VCF file\n");
	fprintf(lfp, "The following meta-information was present in the VCF file:\n");
	while( (nread = getline(&line, &len, stdin)) > 0 ){
		int fieldn = 1;

		// drop meta-info lines to log
		if(strncmp(line, "##", 2) == 0){
			fprintf(lfp, "%s", line);
			continue;
		}

		// read header, parse sample ids
		if(strncmp(line, "#", 1) == 0){
			linecp = strdupa(line);
			fixEndings(linecp);

			while((field = strsep(&linecp, "\t")) != NULL){
				fieldn++;
			}

			// 9 initial fields + 1 empty after terminal \t
			nsamples = fieldn-11;
			fprintf(lfp, "Header line read.\n");
			fprintf(lfp, "In total, %d samples were seen in the VCF.\n", nsamples);

			break;
		}

		// a line doesn't start with #?
		fprintf(stderr, "ERROR: malformed VCF file?");
		return(1);
	}

	// line buffer size is not critical as it will realloc if needed
	char *linev;
	size_t linesize = nsamples * 40;
	linev = (char*)malloc(linesize * sizeof(char));
	if(linev == NULL){
		fprintf(stderr, "ERROR: unable to allocate buffer");
		exit(1);
	}

	// outbut buffer - 3 chars per person
	char outbuf[3*nsamples+1], outbuf1[3*nsamples+1], outbuf2[3*nsamples+1];
	fprintf(lfp, "Limiting output buffers to %d characters, + info fields.\n", 3*nsamples);
	for(int i=0; i<nsamples; i++){
		outbuf[3*i] = '\t';
		outbuf1[3*i] = '\t';
		outbuf2[3*i] = '\t';
	}
	outbuf[3*nsamples] = '\0';
	outbuf1[3*nsamples] = '\0';
	outbuf2[3*nsamples] = '\0';

	// continue reading VCF: genotypes
	char *snppos, *alref, *alalt;
	while( (nread = getline(&linev, &linesize, stdin)) > 0 ){
		int fieldn = 1;
		int i = 0;

		fixEndings(linev);
		char *linecp2 = linev;

		// process all GTs for this SNP
		while((field = strsep(&linecp2, "\t"))[0] != '\0'){
			if(fieldn>9){
				// store actual haplotypes:
				// (48 = ASCII '0')
				if(field[0]<48 || field[2]<48){
					// recode strange genos and missing as 0
					// NOTE: half-calls become missing!
					if(field[0]!='.' && field[2]!='.'){
						fprintf(lfp, "Warning: unusual genotype %d|%d in field %d.\n",
								field[0], field[2], fieldn);
					}
					nmiss++;

					outbuf[3*i+1] = '0';
					outbuf[3*i+2] = '0';
					outbuf1[3*i+1] = '0';
					outbuf1[3*i+2] = '0';
					outbuf2[3*i+1] = '0';
					outbuf2[3*i+2] = '0';
				} else {
					// map 0->A, 1->B
					outbuf[3*i+1] = field[0]+17;
					outbuf[3*i+2] = field[2]+17;
					outbuf1[3*i+1] = field[0]+17;
					outbuf1[3*i+2] = field[0]+17;
					outbuf2[3*i+1] = field[2]+17;
					outbuf2[3*i+2] = field[2]+17;
				}

				i++;
				if(i > nsamples){
					fprintf(stderr, "ERROR: header provided only %d samples\n", nsamples);
					exit(1);
				}

			// parse info fields ("map")
			} else if(fieldn==1){
				fprintf(of, "%s\t", field);
				fprintf(of1, "%s\t", field);
				fprintf(of2, "%s\t", field);
			} else if(fieldn==2){
				snppos = field;
			} else if(fieldn==4){
				alref = field;
			} else if(fieldn==5){
				alalt = field;
			}

			fieldn++;
		}
		if(i != nsamples){
			fprintf(stderr, "ERROR: found only %d samples\n", i);
			exit(1);
		}

		// print varID, cM position, bp position
		fprintf(of, "%s_%s_%s\t0\t%s", snppos, alref, alalt, snppos);
		fprintf(of1, "%s_%s_%s\t0\t%s", snppos, alref, alalt, snppos);
		fprintf(of2, "%s_%s_%s\t0\t%s", snppos, alref, alalt, snppos);

		// flush output buffer:
		fprintf(of, "%s\n", outbuf);
		fprintf(of1, "%s\n", outbuf1);
		fprintf(of2, "%s\n", outbuf2);

		nlines++;
	}

	fprintf(lfp, "Scan complete. Total markers read: %d\n", nlines);
	fprintf(lfp, "In total, %d genotypes were set to missing.\n", nmiss);

	free(linev);
	free(line);

	fclose(lfp);
	fclose(of);
	fclose(of1);
	fclose(of2);

	return(0);
}
