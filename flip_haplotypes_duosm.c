#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>

using namespace std;

// USAGE: zcat onechrom.vcf.gz | ./flip_haplotypes_duos | gzip > out.vcf.gz
// parents.txt format MUST be:  fetid  momid, tab-separated
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

int main() {
	// INIT
	FILE *lfp, *ffp;
	char logFile[] = "log.txt";
	char famFile[] = "parents.txt";

	lfp = fopen(logFile, "w");
	ffp = fopen(famFile, "r");
	if(lfp == NULL){
		fprintf(stderr, "ERROR: Can't open log file!\n");
		exit(1);
	}
	if(ffp == NULL){
		fprintf(stderr, "ERROR: Can't open family info file!\n");
		exit(1);
	}


	// READ
	// reading variables:
	char *line, *field, *linecp, *linef, *fieldf;
	size_t len = 0;
	size_t lenf = 0;
	ssize_t nread, nreadf;

	// sample variables:
	vector <string> moms, fets;
	int nsamples;

	// read family info file:
	while( (nreadf = getline(&linef, &lenf, ffp)) > 0){
		linecp = strdupa(linef);
		fixEndings(linecp);

		fieldf = strsep(&linecp, "\t");
		checkField(fieldf);
		fets.push_back(fieldf);
		
		fieldf = strsep(&linecp, "\t");
		checkField(fieldf);
		moms.push_back(fieldf);
	}

	// C-arrays for samples:
	int ntoread = moms.size();
	char *momsc[ntoread], *fetsc[ntoread];

	fprintf(lfp, "Family info file read. The following duos were found:\n");
	for(int i=0; i<ntoread; i++){
		momsc[i] = new char[moms[i].size() + 1];
		fetsc[i] = new char[fets[i].size() + 1];
		strcpy(momsc[i], moms[i].c_str());
		strcpy(fetsc[i], fets[i].c_str());

		fprintf(lfp, "%s is a parent of %s \n", momsc[i], fetsc[i]);
	}
	fprintf(lfp, "Total amount of duos found: %d\n", ntoread);

	// positions of good parents in VCF
	// (init to -1 for error catching)
	int momspos[ntoread], fetspos[ntoread];
	for(int i=0; i<ntoread; i++){
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

			// check each ID with 'good' ids from fam file
			// and store the corresponding column index:
			while((field = strsep(&linecp, "\t")) != NULL){
				// skip info fields
				if(fieldn<10){
					fieldn++;
					continue;
				}
				
				for(int i=0; i<ntoread; i++){
					// note: multigenerational structures will result in errors
					if(strcmp(momsc[i], field)==0){
						momspos[i] = fieldn-10;
					} else if(strcmp(fetsc[i], field)==0){
						fetspos[i] = fieldn-10;
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
		fprintf(stderr, "ERROR: malformed VCF file?");
		return(1);
	}

	// header of output:
	printf("#CHROM\tPOS\tID\tREF\tALT\tINFO");
	for(int i=0; i<ntoread; i++){
		if(momspos[i]<0 || fetspos[i]<0){
			fprintf(stderr, "ERROR: duo %s, %s was not found in VCF.\n", momsc[i], fetsc[i]);
			fprintf(stderr, "their positions are %d, %d.\n", momspos[i], fetspos[i]);
			exit(1);
		} else {
			fprintf(lfp, "All duos were found in VCF.\n");
			printf("\t%s\t%s", momsc[i], fetsc[i]);
		}
		delete momsc[i];
		delete fetsc[i];
	}
	printf("\n");

	// prep for reading genotypes
	int hapA[nsamples], hapB[nsamples];
	int transmM, untransmM, transmP;
	bool missing;
	int nlines = 0;
	int merrs = 0;
	int nmiss = 0;

	char *linev;
	size_t linesize = nsamples * 40;
	linev = (char*)malloc(linesize * sizeof(char));
	if(linev == NULL){
		fprintf(stderr, "ERROR: unable to allocate buffer");
		exit(1);
	}

	// continue to actual genotypes
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
				// print out info fields:
				printf("%s", field);
			} else if(fieldn < 6 || fieldn==8){
				printf("\t%s", field);
			}
			
			fieldn++;
		}
		nlines++;

		// re-phase the haplotypes, trio-by-trio
		for(int i = 0; i<ntoread; i++){
			missing = false;
			transmM = hapB[fetspos[i]];
			transmP = hapA[fetspos[i]];
			untransmM = hapA[momspos[i]] + hapB[momspos[i]];

			// if duo has more dosage than 4, one or more calls were missing
			if(transmM + transmP + untransmM > 4){
				fprintf(lfp, "Warning: duo %d had missing genotypes.\n", i);
				missing = true;
			}
			untransmM = untransmM - transmM;

			if(untransmM<0 || untransmM>1){
				// Mendelian error for the entire duo
				printf("\t.|.");
				printf("\t.|.");
				merrs++;
			} else if (missing){
				printf("\t.|.");
				printf("\t.|.");
				nmiss++;	
			} else {
				// print: matT|matU  fetM|fetP
				printf("\t%d|%d", transmM, untransmM);
				printf("\t%d|%d", transmM, transmP);
			}
		}
		printf("\n");
	}
	fprintf(lfp, "Scan complete. Total markers read: %d\n", nlines);
	fprintf(lfp, "In total, %d Mendelian errors were observed.\n", merrs);
	fprintf(lfp, "In total, %d markers x duos were missing.\n", nmiss);

	free(linev);

	free(line);
	free(linef);
	fclose(lfp);
	fclose(ffp);

	return(0);
}
