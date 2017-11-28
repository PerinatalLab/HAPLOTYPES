#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>

using namespace std;

// USAGE: zcat onechrom.vcf.gz | ./flip_haplotypes | gzip

// replace \n with \t to fix reading of last field:
void fixEndings(char *s){
	s[strlen(s)-1] = '\t';
}
// check if field is non-empty
void checkField(char *f){
	if(f[0]=='\0'){
		printf("ERROR: family file appears malformed\n");
		exit(1);
	}
}

int main() {
	// INIT
	FILE *lfp, *ffp;
	char logFile[] = "log.txt";
	char famFile[] = "parents.txt";

	lfp = fopen(logFile, "w");
	ffp = fopen(famFile, "r");
	if(lfp == NULL){
		printf("ERROR: Can't open log file!\n");
		exit(1);
	}
	if(ffp == NULL){
		printf("ERROR: Can't open family info file!\n");
		exit(1);
	}


	// READ
	// reading variables:
	char *line, *field, *linecp, *linef, *fieldf;
	size_t len = 0;
	size_t lenf = 0;
	ssize_t nread, nreadf;

	// sample variables:
	vector <string> moms, dads, fets;
	int nsamples;

	// read family info file:
	// format MUST be:  dadid  momid  fetid
	while( (nreadf = getline(&linef, &lenf, ffp)) > 0){
		linecp = strdupa(linef);
		fixEndings(linecp);
		
		fieldf = strsep(&linecp, "\t");
		checkField(fieldf);
		dads.push_back(fieldf);

		fieldf = strsep(&linecp, "\t");
		checkField(fieldf);
		moms.push_back(fieldf);

		fieldf = strsep(&linecp, "\t");
		checkField(fieldf);
		fets.push_back(fieldf);
	}

	// C-arrays for samples:
	int ntoread = moms.size();
	char *momsc[ntoread], *dadsc[ntoread], *fetsc[ntoread];

	fprintf(lfp, "Family info file read. The following trios were found:\n");
	for(int i=0; i<ntoread; i++){
		momsc[i] = new char[moms[i].size() + 1];
		dadsc[i] = new char[dads[i].size() + 1];
		fetsc[i] = new char[fets[i].size() + 1];
		strcpy(momsc[i], moms[i].c_str());
		strcpy(dadsc[i], dads[i].c_str());
		strcpy(fetsc[i], fets[i].c_str());

		fprintf(lfp, "%s and %s are parents of %s \n", momsc[i], dadsc[i], fetsc[i]);
	}
	fprintf(lfp, "Total amount of trios found: %d\n", ntoread);

	// positions of good parents in VCF
	// (init to -1 for error catching)
	int momspos[ntoread], dadspos[ntoread], fetspos[ntoread];
	for(int i=0; i<ntoread; i++){
		dadspos[i] = -1;
		momspos[i] = -1;
		fetspos[i] = -1;
	}

	// read VCF:
	fprintf(lfp, "Reading VCF file\n");
	while( (nread = getline(&line, &len, stdin)) > 0 ){
		int fieldn = 1;

		// drop info lines
		if(strncmp(line, "##", 2) == 0){
			continue;
		}

		// read header, parse sample ids
		if(strncmp(line, "#", 1) == 0){
			linecp = strdupa(line);
			fixEndings(linecp);

			// check each ID with 'good' ids from fam file:
			while((field = strsep(&linecp, "\t")) != NULL){
				// skip info fields
				if(fieldn<10){
					fieldn++;
					continue;
				}
				
				for(int i=0; i<ntoread; i++){
					if(strcmp(momsc[i], field)==0){
						momspos[i] = fieldn-10;
						fprintf(lfp, "found sample %s in good moms\n", field);
						break;
					} else if(strcmp(dadsc[i], field)==0){
						dadspos[i] = fieldn-10;
						fprintf(lfp, "found sample %s in good dads\n", field);
						break;
					} else if(strcmp(fetsc[i], field)==0){
						fetspos[i] = fieldn-10;
						fprintf(lfp, "found sample %s in good fets\n", field);
						break;
					}
				}
				
				fieldn++;
			}
			// 9 initial fields + 1 empty after terminal \t
			nsamples = fieldn-11;
			fprintf(lfp, "Header line read.\n");
			fprintf(lfp, "In total, read %d samples from VCF\n", nsamples);

			break;
		}

		// a line doesn't start with #?
		printf("ERROR: malformed VCF file?");
		return(1);
	}

	// header of output:
	printf("#CHROM\tPOS\tID\tREF\tALT\tINFO");
	for(int i=0; i<ntoread; i++){
		if(momspos[i]<0 || fetspos[i]<0 || dadspos[i]<0){
			printf("ERROR: trio %s, %s, %s was not found in VCF.\n", momsc[i], dadsc[i], fetsc[i]);
			exit(1);
		} else {
			printf("\t%s\t%s\t%s", momsc[i], dadsc[i], fetsc[i]);
		}
		delete dadsc[i];
		delete momsc[i];
		delete fetsc[i];
	}
	printf("\n");

	// prep for reading genotypes
	int hapA[nsamples], hapB[nsamples];
	int transmM, transmP, untransmM, untransmP;
	int nlines = 0;
	int merrs = 0;

	// continue to actual genotypes
	while( (nread = getline(&line, &len, stdin)) > 0 ){
		int fieldn = 1;

		fixEndings(line);
		linecp = strdupa(line);

		// process all GTs for this SNP
		while((field = strsep(&linecp, "\t")) != NULL){
			if(fieldn>9){
				// store actual haplotypes:
				// (48 = ASCII '0')
				hapA[fieldn-10] = field[0]-48;
				hapB[fieldn-10] = field[2]-48;
			} else if(fieldn==1){
				// print out info fields:
				printf("%s", field);
			} else if(fieldn < 6 || fieldn==8){
				printf("\t%s", field);
			}
			
			fieldn++;
		}
		nlines++;

		// re-phase the haplotypes
		for(int i = 0; i<ntoread; i++){
			transmM = hapB[fetspos[i]];
			transmP = hapA[fetspos[i]];
			untransmM = hapA[momspos[i]] + hapB[momspos[i]] - transmM;
			untransmP = hapA[dadspos[i]] + hapB[dadspos[i]] - transmP;

			if(untransmM<0 || untransmM>1 || untransmP<0 || untransmP>1){
				// Mendelian error for the entire trio
				printf("\t.|.");
				printf("\t.|.");
				printf("\t.|.");
				merrs++;
			} else {
				// print: matT|matU  patT|patU  fetM|fetP
				printf("\t%d|%d", transmM, untransmM);
				printf("\t%d|%d", transmP, untransmP);
				printf("\t%d|%d", transmM, transmP);
			}
		}
		printf("\n");
	}
	fprintf(lfp, "Scan complete. Total genotypes read: %d\n", nlines);
	fprintf(lfp, "In total, %d Mendelian errors were observed.\n", merrs);

	free(line);
	free(linef);
	fclose(lfp);
	fclose(ffp);

	return(0);
}
