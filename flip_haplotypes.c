#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <zlib.h>
#include <vector>
#include <string>

using namespace std;

#define version 11

// USAGE: zcat onechrom.vcf.gz | ./flip_haplotypes famfile logfile outfile2 | gzip > out.vcf.gz
// parents.txt format MUST be:  fetid  dadid  momid  pregid  whatever, tab-separated
// missing parents encoded by 0
// will create log.txt

// replace \n with \t to fix reading of last field:
void fixEndings(char *s){
	s[strlen(s)-1] = '\t';
}
// check if field is non-empty
void checkField(char *f){
	if(f[0]=='\0'){
		fprintf(stderr, "ERROR: family file appears malformed\n");
		exit(1);
	}
}

int main(int argc, char *argv[]) {
	// INIT
	// arg check
	char *famFile, *logFile, *miniOut;
	if(argc != 4){
		fprintf(stderr, "ERROR: script needs 3 arguments, %d supplied.\n", argc-1);
		exit(1);
	} else {
		famFile = argv[1];
		logFile = argv[2];
		miniOut = argv[3];
	}

	FILE *lfp, *ffp;
	gzFile ofp;

	lfp = fopen(logFile, "w");
	ffp = fopen(famFile, "r");
	ofp = gzopen(miniOut, "w");
	int gzerr = gzbuffer(ofp, 32000);
	if(lfp == NULL){
		fprintf(stderr, "ERROR: Can't create log file!\n");
		exit(1);
	}
	if(ffp == NULL){
		fprintf(stderr, "ERROR: Can't open family info file!\n");
		exit(1);
	}
	if(ofp == NULL){
		fprintf(stderr, "ERROR: Can't create mini-result file!\n");
		exit(1);
	}
	if(gzerr != 0){
		fprintf(stderr, "ERROR: Can't resize zlib buffer!\n");
		exit(1);
	}
	fprintf(lfp, "Starting log for run.\n");
	fprintf(lfp, "Run command: %s %s %s %s\n", argv[0], argv[1], argv[2], argv[3]);

	// READ
	// reading variables:
	char *line, *field, *linecp, *linef, *fieldf;
	size_t len = 0;
	size_t lenf = 0;
	ssize_t nread, nreadf;

	// sample variables:
	vector <string> moms, dads, fets, pregs;
	int nsamples;

	// read family info file:
	nreadf = getline(&linef, &lenf, ffp);
	while( (nreadf = getline(&linef, &lenf, ffp)) > 0){
		linecp = strdupa(linef);
		fixEndings(linecp);

		fieldf = strsep(&linecp, "\t");
		checkField(fieldf);
		if(fieldf=="0"){
			fprintf(stderr, "ERROR: missing fetal IDs not allowed\n");
			exit(1);
		}
		fets.push_back(fieldf);
		
		fieldf = strsep(&linecp, "\t");
		checkField(fieldf);
		dads.push_back(fieldf);

		fieldf = strsep(&linecp, "\t");
		checkField(fieldf);
		moms.push_back(fieldf);

		fieldf = strsep(&linecp, "\t");
		checkField(fieldf);
		pregs.push_back(fieldf);
	}

	// copy sample names to fixed C-style arrays
	// and initialize arrays of trio positions in VCF
	int ntoread = moms.size();
	char *momsc[ntoread], *dadsc[ntoread], *fetsc[ntoread], *pregsc[ntoread];
	int momspos[ntoread], dadspos[ntoread], fetspos[ntoread];
	int nmatduos = 0;
	int npatduos = 0;

	fprintf(lfp, "Family info file read. The following trios were found:\n");
	for(int i=0; i<ntoread; i++){
		// sample names
		momsc[i] = new char[moms[i].size() + 1];
		dadsc[i] = new char[dads[i].size() + 1];
		fetsc[i] = new char[fets[i].size() + 1];
		pregsc[i] = new char[pregs[i].size() + 1];
		strcpy(momsc[i], moms[i].c_str());
		strcpy(dadsc[i], dads[i].c_str());
		strcpy(fetsc[i], fets[i].c_str());
		strcpy(pregsc[i], pregs[i].c_str());

		// positions of samples in VCF
		// (init to -2 as an error marker)
		dadspos[i] = -2;
		momspos[i] = -2;
		fetspos[i] = -2;

		// check and log missings
		// (-1 will tell phaser to switch to duo mode)
		if(strcmp(momsc[i], "0")==0 && strcmp(dadsc[i], "0")==0){
			fprintf(stderr, "ERROR: both parents are missing for pregnancy %s\n", pregsc[i]);
			exit(1);
		} else if (strcmp(momsc[i], "0")==0){
			fprintf(lfp, "%s is father of %s\n", dadsc[i], fetsc[i]);
			momspos[i] = -1;
			npatduos++;
		} else if (strcmp(dadsc[i], "0")==0){
			fprintf(lfp, "%s is mother of %s\n", momsc[i], fetsc[i]);
			dadspos[i] = -1;
			nmatduos++;
		} else {
			fprintf(lfp, "%s and %s are parents of %s \n", momsc[i], dadsc[i], fetsc[i]);
		}

	}
	fprintf(lfp, "Total amount of trios found: %d\n", ntoread);
	fprintf(lfp, "Total amount of maternal duos found: %d\n", nmatduos);
	fprintf(lfp, "Total amount of paternal duos found: %d\n", npatduos);


	// OUTPUT PRINTING STARTS HERE
	// read VCF:
	fprintf(lfp, "Reading VCF file\n");
	while( (nread = getline(&line, &len, stdin)) > 0 ){
		int fieldn = 1;

		// push meta-info lines straight to output
		if(strncmp(line, "##", 2) == 0){
			printf("%s", line);
			continue;
		}
		printf("##haplotype_source=flipper_v%d\n", version);

		// read header, parse sample ids
		if(strncmp(line, "#", 1) == 0){
			linecp = strdupa(line);
			fixEndings(linecp);

			while((field = strsep(&linecp, "\t")) != NULL){
				// skip info fields
				if(fieldn<10){
					fieldn++;
					continue;
				}
				
				// check each ID with 'good' ids from fam file
				// and store the corresponding column index
				// (NOTE: multigenerational structures will result in errors)
				for(int i=0; i<ntoread; i++){
					if(strcmp(momsc[i], field)==0){
						momspos[i] = fieldn-10;
					} else if(strcmp(dadsc[i], field)==0){
						dadspos[i] = fieldn-10;
					} else if(strcmp(fetsc[i], field)==0){
						fetspos[i] = fieldn-10;
					}
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

	// check if all samples were found
	// and print header line to output
	printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for(int i=0; i<ntoread; i++){
		if(momspos[i]==-2 || fetspos[i]==-2 || dadspos[i]==-2){
			fprintf(stderr, "ERROR: trio %s, %s, %s was not found in VCF.\n", momsc[i], dadsc[i], fetsc[i]);
			fprintf(stderr, "their positions are %d, %d, %d.\n", momspos[i], dadspos[i], fetspos[i]);
			exit(1);
		} else {
			printf("\t%s_%s\t%s_%s\t%s_%s", pregsc[i], fetsc[i], pregsc[i], momsc[i], pregsc[i], dadsc[i]);
		}
		delete dadsc[i];
		delete momsc[i];
		delete fetsc[i];
	}
	printf("\n");
	fprintf(lfp, "All trios were found in VCF.\n");

	// prep for reading genotypes
	int hapA[nsamples], hapB[nsamples];
	int transmM, transmP, untransmM, untransmP;
	bool missing, merr;
	int nlines = 0;
	int nmiss = 0;

	// line buffer size is not critical as it will realloc if needed
	char *linev;
	size_t linesize = nsamples * 40;
	linev = (char*)malloc(linesize * sizeof(char));
	if(linev == NULL){
		fprintf(stderr, "ERROR: unable to allocate buffer");
		exit(1);
	}

	// outbut buffer - 4 chars per person
	char outbuf[4*3*ntoread];
	// second output buffer
	char gzbuf[2*ntoread];
	fprintf(lfp, "Limiting output buffer to %d characters, + info fields.\n", 4*3*ntoread);
	fprintf(lfp, "Limiting zlib buffer to %d characters, + info fields.\n", 2*ntoread);
	for(int i=0; i<3*ntoread; i++){
		outbuf[4*i] = '\t';
		outbuf[4*i + 2] = '|';
	}
	for(int i=0; i<ntoread; i++){
		gzbuf[2*i] = '\t';
	}

	// continue reading VCF: genotypes
	while( (nread = getline(&linev, &linesize, stdin)) > 0 ){
		int fieldn = 1;

		fixEndings(linev);
		char *linecp2 = linev;

		// process all GTs for this SNP
		while((field = strsep(&linecp2, "\t"))[0] != '\0'){
			if(fieldn>9){
				// store actual haplotypes:
				// (48 = ASCII '0')
				if(field[0]<48 || field[2]<48){
					// recode missing as 9
					fprintf(lfp, "Warning: unusual genotype %d|%d in field %d.\n", field[0], field[2], fieldn);
					field[0] = 57;
					field[2] = 57;
				}
				hapA[fieldn-10] = field[0]-48;
				hapB[fieldn-10] = field[2]-48;
			} else if(fieldn==1){
				// print out info fields
				// (NOTE: could switch to sprintf buffer here)
				printf("%s", field);
			} else if(fieldn==2){
				printf("\t%s", field);
				gzprintf(ofp, "%s", field);
			} else if(fieldn<6){
				printf("\t%s", field);
				gzprintf(ofp, "\t%s", field);
			} else if(fieldn==9){
				// change FORMAT field
				printf("\tGT");
			} else {
				printf("\t%s", field);
			}
			
			fieldn++;
		}
		nlines++;

		// PHASER
		// re-phase the haplotypes, trio-by-trio
		for(int i = 0; i<ntoread; i++){
			missing = false;
			merr = false;
			transmM = hapB[fetspos[i]];
			transmP = hapA[fetspos[i]];

			// if fetal genotypes are missing, can't phase
			if(transmM + transmP > 2){
				missing = true;
			} else {
				if(momspos[i]<0){
					// switch to duo mode if a parent is missing
					untransmM = -2;
				} else {
					// phase
					untransmM = hapA[momspos[i]] + hapB[momspos[i]] - transmM;
					// check for Mendelian errors or weird input
					// NOTE: missing input will also result in this
					merr = untransmM<0 || untransmM>1;
				}
				if(dadspos[i]<0){
					untransmP = -2;
				} else {
					untransmP = hapA[dadspos[i]] + hapB[dadspos[i]] - transmP;
					merr = untransmP<0 || untransmP>1 || merr;
				}
			}

			// mask the entire trio on MEs or missing input
			// (-2 + 48 = ASCII '.')
			if(merr || missing){
				transmM = transmP = -2;
				untransmM = untransmP = -2;
				nmiss++;
			}

			// convert to ASCII char codes
			transmM += 48;
			transmP += 48;
			untransmM += 48;
			untransmP += 48;

			// print: fetM|fetp  matT|matU  patT|patU
			outbuf[i*12+1] = transmM; outbuf[i*12+3] = transmP;
			outbuf[i*12+5] = transmM; outbuf[i*12+7] = untransmM;
			outbuf[i*12+9] = transmP; outbuf[i*12+11] = untransmP;

			gzbuf[i*2+1] = untransmM; 
		}

		// flush output buffer:
		printf("%s\n", outbuf);
		gzprintf(ofp, "%s\n", gzbuf);
	}

	fprintf(lfp, "Scan complete. Total markers read: %d\n", nlines);
	fprintf(lfp, "Mendelian errors or missing input resulted in %d trios x markers set to missing.\n", nmiss);

	free(linev);

	free(line);
	free(linef);
	fclose(lfp);
	fclose(ffp);
	gzclose(ofp);

	return(0);
}
