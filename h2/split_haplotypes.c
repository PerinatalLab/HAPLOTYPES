#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>


#define version 11


// This script shall take in a family file, format compatible with above,
// allowing trios or duos with one 0 parent,
// and a vcf.gz or vcf file, with haplotypes aligned M|P T|U T|U,
// and produce three .bed files for easy input to PLINK, plus a .bim and .fam file.

// USAGE: zcat onechrom.vcf.gz | ./split_haplotypes logfile outstem

// replace \n with \t to fix reading of last field:
void fixEndings(char *s){
	s[strlen(s)-1] = '\t';
}

int main(int argc, char *argv[]) {
	// INIT
	// arg check
	char *logFile, *outFile, *outFile1, *outFile2, *outFam, *outBim;
	if(argc != 3){
		fprintf(stderr, "ERROR: script needs 2 arguments, %d supplied.\n", argc-1);
		exit(1);
	} else {
		logFile = argv[1];
		outFile = strdup(argv[2]);
		outFile1 = strdup(outFile);
		outFile2 = strdup(outFile);
		outBim = strdup(argv[2]);
		outFam = strdup(argv[2]);
		strcat(outFile, ".bed");
		strcat(outFile1, "1.bed");
		strcat(outFile2, "2.bed");
		strcat(outBim, ".bim");
		strcat(outFam, ".fam");
	}

	FILE *lfp, *of, *of1, *of2, *ofB, *ofF;

	lfp = fopen(logFile, "w");
	of = fopen(outFile, "wb");
	of1 = fopen(outFile1, "wb");
	of2 = fopen(outFile2, "wb");
	ofB = fopen(outBim, "w");
	ofF = fopen(outFam, "w");

	if(lfp == NULL){
		fprintf(stderr, "ERROR: Can't create log file!\n");
		exit(1);
	}
	if(of == NULL || of1 == NULL || of2 == NULL || ofB == NULL || ofF == NULL){
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
			linecp = strdup(line);
			fixEndings(linecp);

			while((field = strsep(&linecp, "\t")) != NULL){
				if(fieldn>9 && strcmp(field, "")!=0){
					fprintf(ofF, "%s\t%s\t0\t0\t0\t0\n", field, field);
				}
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

	// outbut buffer - 1 byte per 4 people
	int nbytes = (nsamples + 3)/4;
	unsigned char outbuf[nbytes], outbuf1[nbytes], outbuf2[nbytes];
	outbuf[nbytes] = outbuf1[nbytes] = outbuf2[nbytes] = '\0';
	fprintf(lfp, "Limiting output buffers to %d bytes, + info fields.\n", nbytes);

	// continue reading VCF: genotypes
	char *snppos, *alref, *alalt;
	unsigned char byte, byte1, byte2;
       	byte = byte1 = byte2 = 0x00;
	int pos;
	char const magic[] = {0x6c, 0x1b, 0x01, 0x00};
	fprintf(of, "%s", magic);
	fprintf(of1, "%s", magic);
	fprintf(of2, "%s", magic);
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
				pos = (i % 4) * 2;
				if(field[0]<48 || field[2]<48){
					// on strange or missing genos, force missing to all outputs
					// NOTE: half-calls also become missing!
					if(field[0]!='.' && field[2]!='.'){
						fprintf(lfp, "Warning: unusual genotype %d|%d in field %d.\n",
								field[0], field[2], fieldn);
					}

					nmiss++;
					byte += 0x01 << pos;
					byte1 += 0x01 << pos;
					byte2 += 0x01 << pos;
				} else {
					// for haplotype outputs:
					if(field[0]=='0'){
						byte1 += 0x00 << pos;
					} else {
						byte1 += 0x03 << pos;
					}

					if(field[2]=='0'){
						byte2 += 0x00 << pos;
					} else {
						byte2 += 0x03 << pos;
					}

					// for main output:
					if(field[0]=='0' && field[2]=='0'){
						byte += 0x00 << pos;
					} else if (field[0]=='1' && field[2]=='1'){
						byte += 0x03 << pos;
					} else if (field[0]=='1' || field[2]=='1'){
						byte += 0x02 << pos;
					} else {
						fprintf(stderr, "ERROR: weird genotype at sample %d\n", i);
						exit(1);
					}
				}

				if(i % 4 == 3){
					// flush to bigger buffer, reset byte
					outbuf[i/4] = byte;
					outbuf1[i/4] = byte1;
					outbuf2[i/4] = byte2;
					byte = byte1 = byte2 = 0x00;
				} else if(i == nsamples-1){
					// flush to bigger buffer, reset byte
					outbuf[nbytes-1] = byte;
					outbuf1[nbytes-1] = byte1;
					outbuf2[nbytes-1] = byte2;
					byte = byte1 = byte2 = 0x00;
				}

				i++;
				if(i > nsamples){
					fprintf(stderr, "ERROR: header provided only %d samples\n", nsamples);
					exit(1);
				}

			// parse info fields ("map/bim")
			} else if(fieldn==1){
				fprintf(ofB, "%s\t%s:", field, field);
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

		// print varID, cM pos, bp pos, alleles to .bim file
		fprintf(ofB, "%s:%s:%s\t0\t%s\t%s\t%s\n",
				snppos, alref, alalt, snppos, alref, alalt);

		// flush output buffer:
		fwrite(outbuf, sizeof outbuf, 1, of);
		fwrite(outbuf1, sizeof outbuf1, 1, of1);
		fwrite(outbuf2, sizeof outbuf2, 1, of2);

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
	fclose(ofB);
	fclose(ofF);

	return(0);
}
