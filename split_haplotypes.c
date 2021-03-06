#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>

#define version 11

using namespace std;

// This script shall take in a family file, format compatible with above,
// allowing trios or duos with one 0 parent,
// and a vcf.gz or vcf file, with haplotypes aligned M|P T|U T|U,
// and produce three .bed files for easy input to PLINK, plus a .bim and .fam file.

// USAGE: zcat onechrom.vcf.gz | ./split_haplotypes famfile logfile outstem


// prepare I/O file streams
FILE *ffp, *lfp, *ofM1, *ofM2, *ofP1, *ofP2, *ofB, *ofFM, *ofFP, *ofFF, *ofKM, *ofKP;

// open file stream and do checks
FILE *openOut(char *file, const char *suffix, const char *mode){
	char tmpstr[256];
       	strcpy(tmpstr, file);
	strcat(tmpstr, suffix);
	FILE *stream = fopen(tmpstr, mode);
	if(stream == NULL){
		fprintf(stderr, "ERROR: Can't create output file %s!\n", tmpstr);
		exit(1);
	} else {
		fprintf(lfp, "Saving output to file %s.\n", tmpstr);
	}
	return(stream);
}

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

// check if sample was duplicated
// and return 1 if it's new
int checkDupl(int posGood, char **arrGood, char **arrBad1, char **arrBad2){
	for(int j=0; j<posGood; j++){
		if(strcmp(arrGood[posGood], arrGood[j])==0){
			fprintf(lfp, "Warning: sample %s was repeated in .fam file\n", arrGood[posGood]);
			return(0);
		}
		if(strcmp(arrGood[posGood], arrBad1[j])==0 || strcmp(arrGood[posGood], arrBad2[j])==0){
			fprintf(stderr, "ERROR: sample %s has conflicting roles in .fam file\n", arrGood[posGood]);
			exit(1);
		}
	}
	return(1);
}
// shorten the types of family relations
enum role {M=1, D, F};

int main(int argc, char *argv[]) {
	// INIT
	// arg check
	char *famFile, *logFile;
	if(argc != 4){
		fprintf(stderr, "ERROR: script needs 3 arguments, %d supplied.\n", argc-1);
		exit(1);
	} else {
		famFile = argv[1];
		logFile = argv[2];
	}

	ffp = fopen(famFile, "r");
	lfp = fopen(logFile, "w");
	if(lfp == NULL){
		fprintf(stderr, "ERROR: Can't create log file!\n");
		exit(1);
	}
	if(ffp == NULL){
	        fprintf(stderr, "ERROR: Can't open family info file!\n");
	        exit(1);
	}

	ofM1 = openOut(argv[3], "M1.bed", "wb");
	ofM2 = openOut(argv[3], "M2.bed", "wb");
	ofP1 = openOut(argv[3], "P1.bed", "wb");
	ofP2 = openOut(argv[3], "P2.bed", "wb");

	ofB = openOut(argv[3], ".bim", "w");
	ofFM = openOut(argv[3], "M.fam", "w");
	ofFP = openOut(argv[3], "P.fam", "w");
	ofFF = openOut(argv[3], "F.fam", "w");

	ofKM = openOut(argv[3], "M.tokeep", "w");
	ofKP = openOut(argv[3], "P.tokeep", "w");

	fprintf(lfp, "Haplotype splitter v1\n");
	fprintf(lfp, "Starting log for run.\n");
	fprintf(lfp, "Run command: %s %s %s %s\n", argv[0], argv[1], argv[2], argv[3]);

	// READ
	// reading variables:
	char *line = NULL;
	char *linef = NULL;
	char *field, *linecp, *fieldf;
	size_t len = 0;
	size_t lenf = 0;
	ssize_t nread, nreadf;
        vector <string> moms, dads, fets, pregs;

	// read family info file:
        while( (nreadf = getline(&linef, &lenf, ffp)) > 0){
        	linecp = strdupa(linef);
                fixEndings(linecp);

                fieldf = strsep(&linecp, "\t");
                checkField(fieldf);
                if(strcmp(fieldf, "0")==0){
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
	}

        // copy sample names to fixed C-style arrays
	// and initialize array of family roles in VCF
	int maxntrios = fets.size();
	char *momsc[maxntrios], *dadsc[maxntrios], *fetsc[maxntrios];
	int ntrios, nmatduos, npatduos, nsingletons;
	ntrios = nmatduos = npatduos = nsingletons = 0;
	int nfets, ndads, nmoms;
	nfets = ndads = nmoms = 0;

	fprintf(lfp, "Family info file read.\n");
	for(int i=0; i<maxntrios; i++){
		// sample names
		momsc[i] = new char[moms[i].size() + 1];
		dadsc[i] = new char[dads[i].size() + 1];
		fetsc[i] = new char[fets[i].size() + 1];
		strcpy(momsc[i], moms[i].c_str());
		strcpy(dadsc[i], dads[i].c_str());
		strcpy(fetsc[i], fets[i].c_str());

		// check and log missings
		// log sample repeats, check for conflicting roles
		// count total samples
		// (NOTE: multigenerational structures will result in errors)
		if(strcmp(momsc[i], "0")==0 && strcmp(dadsc[i], "0")==0){
			nsingletons++;
		} else if (strcmp(momsc[i], "0")==0){
			npatduos++;
			ndads += checkDupl(i, dadsc, momsc, fetsc);
		} else if (strcmp(dadsc[i], "0")==0){
			nmatduos++;
			nmoms += checkDupl(i, momsc, dadsc, fetsc);
		} else {
			ntrios++;
			ndads += checkDupl(i, dadsc, momsc, fetsc);
			nmoms += checkDupl(i, momsc, dadsc, fetsc);
		}
		nfets += checkDupl(i, fetsc, momsc, dadsc);
	}
	int nsamples = nmoms + ndads + nfets;
	fprintf(lfp, "Total numbers of unique samples: %d mothers, %d fathers, %d fetuses\n",
			nmoms, ndads, nfets);
	fprintf(lfp, "Total amount of pregnancies provided: %d\n", maxntrios);
	fprintf(lfp, " amount of trios found: %d\n", ntrios);
	fprintf(lfp, " amount of maternal duos found: %d\n", nmatduos);
	fprintf(lfp, " amount of paternal duos found: %d\n", npatduos);
	fprintf(lfp, " amount of singletons found: %d\n", nsingletons);


	// OUTPUT PRINTING STARTS HERE
	// NOTE: all samples in the VCF must be listed in fam file 
	// roles initialized to max expected VCF size (=all samples in .fam unique and present)
	enum role roles[nsamples];
	int nlines = 0;
	int nmiss = 0;
	nfets = ndads = nmoms = 0;

	// read VCF meta-info and header:
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
				// skip info fields
				if(fieldn < 10){
					fieldn++;
					continue;
				}
				if(fieldn-9 > nsamples && strcmp(field, "")!=0){
					fprintf(stderr, "ERROR: VCF has more samples than fam file (field %d)\n", fieldn);
					printf("%s\n", field);
					exit(1);
				}

				// check each ID with ids from fam file
				// store its role and print THE FETAL ID to corresponding .fam output
				// so all .fams will have FETAL IDs, but M.fam and P.fam can be smaller than F.fam
				// K (--keep) files will have correct parent IDs for separating them
				for(int i=0; i<maxntrios; i++){
					if(strcmp(momsc[i], field)==0){
						roles[fieldn-10] = M;
						fprintf(ofKM, "%s\t%s\t%s\t%s\n",
								momsc[i], momsc[i], fetsc[i], fetsc[i]);
						fprintf(ofFM, "%s\t%s\t0\t0\t0\t0\n", fetsc[i], fetsc[i]);
						nmoms++;
						break;
					} else if(strcmp(dadsc[i], field)==0){
						roles[fieldn-10] = D;
						fprintf(ofKP, "%s\t%s\t%s\t%s\n",
								dadsc[i], dadsc[i], fetsc[i], fetsc[i]);
						fprintf(ofFP, "%s\t%s\t0\t0\t0\t0\n", fetsc[i], fetsc[i]);
						ndads++;
						break;
					} else if(strcmp(fetsc[i], field)==0){
						roles[fieldn-10] = F;
						fprintf(ofFF, "%s\t%s\t0\t0\t0\t0\n", field, field);
						nfets++;
						break;
					}
					if(i==maxntrios-1 && strcmp(field, "")!=0){
						fprintf(stderr, "ERROR: VCF sample %s was not found in fam file\n", field);
						exit(1);
					}
				}
				fieldn++;
			}
			
			if(nmoms+nfets+ndads == nsamples){
				fprintf(lfp, "VCF header read. All %d samples from fam file were found.\n", nsamples);
			} else if(nmoms+nfets+ndads < nsamples){
				fprintf(lfp, "Warning: found only %d samples out of %d provided in fam file.\n",
						nmoms+nfets+ndads, nsamples);
			} else {
				fprintf(stderr, "ERROR: found %d samples, expected only %d?\n",
						nmoms+nfets+ndads, nsamples);
				exit(1);
			}
			break;
		}

		// a line doesn't start with #?
		fprintf(stderr, "ERROR: malformed VCF file?\n");
		exit(1);
	}

	// release memory as the .fam file is printed already
	for(int i=0; i<maxntrios; i++){
		delete[] dadsc[i];
		delete[] momsc[i];
		delete[] fetsc[i];
	}

	// double-check if all roles were assigned correctly
	nsamples = nmoms + nfets + ndads;
	for(int i=0; i<nsamples; i++){
		if(roles[i]!=M && roles[i]!=D && roles[i]!=F){
			fprintf(stderr, "ERROR: role at position %d was not assigned?\n", i+10);
			exit(1);
		}
	}
	fprintf(lfp, "All VCF samples were found in .fam file.\n");

	// input line buffer size is not critical as it will realloc if needed
	char *linev;
	size_t linesize = nsamples * 40;
	linev = (char*)malloc(linesize * sizeof(char));
	if(linev == NULL){
		fprintf(stderr, "ERROR: unable to allocate buffer\n");
		exit(1);
	}	

	// outbut buffer - 1 byte per 4 people
	nfets = (nfets + 3)/4;
	nmoms = (nmoms + 3)/4;
	ndads = (ndads + 3)/4;
	unsigned char outbufM1[nfets], outbufM2[nmoms], outbufP1[nfets], outbufP2[ndads];
	fprintf(lfp, "Limiting F/M/P output buffers to %d/%d/%d bytes, + info fields.\n", nfets, nmoms, ndads);

	// magic bytes to identify file as PLINK's .bed:
	char const magic[] = {0x6c, 0x1b, 0x01};
	fwrite(magic, 3, 1, ofM1);
	fwrite(magic, 3, 1, ofM2);
	fwrite(magic, 3, 1, ofP1);
	fwrite(magic, 3, 1, ofP2);

	// continue reading VCF: genotypes
	char *snppos, *alref, *alalt;
	unsigned char byteA, byteB, byteM1, byteM2, byteP1, byteP2;
	int pos1, posM2, posP2;

	while( (nread = getline(&linev, &linesize, stdin)) > 0 ){
		int fieldn = 1;
		int i = 0;

		fixEndings(linev);
		char *linecp2 = linev;
		pos1 = posM2 = posP2 = 0;
       		byteM1 = byteM2 = byteP1 = byteP2 = 0x00;

		// process all GTs for this SNP
		while((field = strsep(&linecp2, "\t"))[0] != '\0'){
			if(fieldn>9){
				// store actual haplotypes:
				if(field[0]<'0' || field[2]<'0' || field[0]>'1' || field[2]>'1'){
					// on strange or missing genos, force missing to all outputs
					// NOTE: half-calls also become missing!
					if(field[0]!='.' && field[2]!='.'){
						fprintf(lfp, "Warning: unusual genotype %d|%d in field %d.\n",
								field[0], field[2], fieldn);
					}

					nmiss++;
					byteA = 0x01;
					byteB = 0x01;
				} else {
					if(field[0]=='0'){
						byteA = 0x03;
					} else {
						byteA = 0x00;
					}

					if(field[2]=='0'){
						byteB = 0x03;
					} else {
						byteB = 0x00;
					}
				}

				// flush haplotype to the corresponding buffer
				if(roles[i] == F){
					// transmitted outputs:
					byteM1 += byteA << (pos1 % 4 * 2);
					byteP1 += byteB << (pos1 % 4 * 2);

					// flush to bigger buffer, reset byte
					if(pos1 % 4 == 3){
						outbufM1[pos1/4] = byteM1;
						outbufP1[pos1/4] = byteP1;
						byteM1 = byteP1 = 0x00;
					}
					pos1++;
				} else if(roles[i] == M){
					// maternal untransmitted output:
					byteM2 += byteB << (posM2 % 4 * 2);

					if(posM2 % 4 == 3){
						outbufM2[posM2/4] = byteM2;
						byteM2 = 0x00;
					}
					posM2++;
				} else if(roles[i] == D){
					// paternal untransmitted output:
					byteP2 += byteB << (posP2 % 4 * 2);

					if(posP2 % 4 == 3){
						outbufP2[posP2/4] = byteP2;
						byteP2 = 0x00;
					}
					posP2++;
				} else {
					fprintf(stderr, "ERROR: undetermined role at position %d\n", i);
					exit(1);
				}

				// on final samples - flush to bigger buffer
				if(i == nsamples-1){
					outbufM1[nfets-1] = byteM1;
					outbufP1[nfets-1] = byteM2;
					outbufM2[nmoms-1] = byteP1;
					outbufP2[ndads-1] = byteP2;
				}

				i++;
				if(i > nsamples){
					fprintf(stderr, "ERROR: expected at most %d genotypes in line %d\n", nsamples, nlines);
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
			fprintf(stderr, "ERROR: line %d had only %d genotypes\n", nlines, i);
			exit(1);
		}

		// print varID, cM pos, bp pos, alleles to .bim file
		fprintf(ofB, "%s:%s:%s\t0\t%s\t%s\t%s\n",
				snppos, alref, alalt, snppos, alalt, alref);

		// flush output buffer:
		fwrite(outbufM1, sizeof outbufM1, 1, ofM1);
		fwrite(outbufM2, sizeof outbufM2, 1, ofM2);
		fwrite(outbufP1, sizeof outbufP1, 1, ofP1);
		fwrite(outbufP2, sizeof outbufP2, 1, ofP2);

		nlines++;
	}

	fprintf(lfp, "Scan complete. Total markers read: %d\n", nlines);
	fprintf(lfp, "In total, %d genotypes were set to missing.\n", nmiss);

	free(linev);
	free(line);
	free(linecp);
	free(linef);

	fclose(ffp);
	fclose(lfp);
	fclose(ofM1);
	fclose(ofM2);
	fclose(ofP1);
	fclose(ofP2);
	fclose(ofB);
	fclose(ofFM);
	fclose(ofFP);
	fclose(ofFF);
	fclose(ofKM);
	fclose(ofKP);

	return(0);
}
