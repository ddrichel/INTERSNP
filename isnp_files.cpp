/**~ isnp_files.cpp
.---------------------------------------------------------------------------.
|  Software: INTERSNP, Genome-Wide Interaction Analysis                     |
|      Site: http://intersnp.meb.uni-bonn.de/                               |
| ------------------------------------------------------------------------- |
|      File: isnp_files.cpp                                                 |
|    Author: André Lacour, Vitalia Schüller                                 |
|   Content: source-code library on reading PLINK file-formats for INTERSNP |
| Copyright (c) 2011-2015, André Lacour                                     |
| ------------------------------------------------------------------------- |
|   License: Distributed under the General Public License (GPL)             |
|            http://www.gnu.org/licenses/                                   |
| This library is distributed in the hope that it will be useful - WITHOUT  |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or     |
| FITNESS FOR A PARTICULAR PURPOSE.                                         |
'---------------------------------------------------------------------------'
**/

using namespace std;


int update(string*, int);

void die(string str) {
    extern fstream errorfile;
    extern fstream logfile;
    errorfile.clear();
    logfile.clear();
    errorfile << "Fatal error! " << str << endl;
	logfile   << "Fatal error! " << str << endl;
	cout      << "Fatal error! " << str << "\n" << endl;
	errorfile.close();
	logfile.close();
	exit(1);
}

void error(string str) {
    extern fstream errorfile;
    extern fstream logfile;
    errorfile.clear();
    logfile.clear();
    errorfile << "Error! " << str << endl;
	logfile   << "Error! " << str << endl;
	cout      << "Error! " << str << endl;
}

void warning(string str) {
    extern fstream errorfile;
    extern fstream logfile;
    errorfile.clear();
    logfile.clear();
	errorfile << "Warning! " << str << endl;
	logfile   << "Warning! " << str << endl;
	cout      << "Warning! " << str << endl;
}

void notice(string str) {
    extern fstream errorfile;
    errorfile.clear();
	errorfile << "Notice! " << str << endl;
}

void dev(string str) {
#if DEV
    extern fstream logfile;
    logfile.clear();
	logfile << "Dev-Notice! " << str << endl;
	cout    << "Dev-Notice! " << str << endl;
#endif
}

void logg(string str) {
    extern fstream logfile;
    logfile.clear();
	logfile << str << endl;
	cout    << str << endl;
}

void die(stringstream& ss) {
    die(ss.str());
}

void error(stringstream& ss) {
    error(ss.str());
    ss.str("");
    ss.clear();
}

void warning(stringstream& ss) {
    warning(ss.str());
    ss.str("");
    ss.clear();
}

void notice(stringstream& ss) {
    notice(ss.str());
    ss.str("");
    ss.clear();
}

void dev(stringstream& ss) {
    dev(ss.str());
    ss.str("");
    ss.clear();
}

void logg(stringstream& ss) {
    logg(ss.str());
    ss.str("");
    ss.clear();
}


uint32_t getline(char*& buf, FILE* stream) {
    uint32_t size = 0;
    uint32_t len = 0;
    do {
        size += BUFSIZ;
        buf = (char*)realloc(buf, size*sizeof(char));
        if ( fgets(buf+len, BUFSIZ*sizeof(char), stream) == NULL ) { free(buf); buf=NULL; return 0; }
        len = strlen(buf);
    } while (buf[len-1]!='\n');
    buf[--len]='\0';
    if (len>0 && buf[len-1]=='\r') buf[--len]='\0';
    return len;
}


bool file_exists(string file) {
    if (FILE* fptr = fopen(file.c_str(), "r")) {
        fclose(fptr);
        return true;
    }
    return false;
}


uint32_t countLines(string file) {
    uint32_t nlines = 0;
    char buff[BUFSIZ] = "";
	FILE* fp = fopen(file.c_str(), "r");
    while ( fgets( buff, sizeof(buff), fp ) != NULL )
        if ( strchr(buff, '\n') != NULL ) nlines++;
	fclose(fp);
    return nlines;
}


uint32_t countColumns(string file) {
    stringstream sstm;
	uint32_t ncols = 0;
	char* s = NULL;
	FILE* fptr = fopen(file.c_str(), "r");
	if (fptr == NULL) { sstm << "File " << file << " does not exist."; die(sstm); }
	if (getline(s, fptr) == 0) { sstm << "First line of file " << file << " is empty! Please check infile."; die(sstm); }
	if (strtok(s, ",\t; ") != NULL) ncols++;
	while (strtok(NULL, ",\t; ") != NULL) ncols++;
	fclose(fptr);
	free(s); s=NULL;
	return ncols;
}


uint32_t countColumnsWithHeader(string file, int head, int lines) {
    stringstream sstm;
	uint32_t ncols = 0;
	uint32_t ncolst = 0;
	FILE* fptr = fopen(file.c_str(), "r");
	int skipped=0;
	int read=0;
	char* s = NULL;
	if (fptr == NULL) { sstm << "File " << file << " does not exist."; die(sstm); }
	while(skipped<head){
	  if (getline(s, fptr) == 0) { sstm << "Line "<<skipped+1<< " of file " << file << " is empty! Please check infile."; die(sstm); }
	  else skipped++;
	  free(s);s=NULL;
	}
	while(read<lines){
	  if(read>1 && ncols!=ncolst){ sstm << "Number of columns in line " << skipped+1 << " is "<<ncols<<" but in line "<<skipped+read+1<<" it is "<<ncolst; die(sstm); }
	  ncolst=0;
	  if (getline(s, fptr) == 0) { sstm << "Line "<<read+skipped+1<< " of file " << file << " is empty! Please check infile."; die(sstm); }
	  else if (read==0){
	    if (strtok(s, ",\t; ") != NULL) ncols++;
	    while (strtok(NULL, ",\t; ") != NULL) ncols++;
	    read++;
	    }
	    else{
	      if (strtok(s, ",\t; ") != NULL) ncolst++;
	      while (strtok(NULL, ",\t; ") != NULL) ncolst++;
	      read++;
	    }
	  free(s);s=NULL;
	}
	return ncols;
}


uint32_t read_fam(string famfile, bool qt, struct PERSON*& person) {
#if PARALLEL
    uint32_t thread = omp_get_thread_num();
#else
    uint32_t thread = 0;
#endif
    stringstream sstm;
	uint32_t nlinesfam = countLines(famfile);

	person = (struct PERSON *) calloc(nlinesfam, sizeof(struct PERSON));
	if (!person) die("memory allocation error in person");

    uint32_t nlines = 0;
    char* chunk = NULL;
    char* s = NULL;

    FILE* fptr = fopen(famfile.c_str(), "r");
	for (uint32_t l=1; !feof(fptr); l++) {
        if (!getline(s, fptr)) continue;

		person[nlines].qcin = true;
		person[nlines].analysis_in = true;

		chunk = strtok(s, " \t");
		if (!chunk) { sstm << "Not enough columns in file " << famfile << " line " << nlines << " column 0."; die(sstm); }
		person[nlines].fid = (char *) malloc((strlen(chunk) + 1) * sizeof(char));
		if (!person[nlines].fid ) die("memory allocation error in person[nlines].fid");
		strcpy(person[nlines].fid, chunk);

        chunk = strtok(NULL, " \t");
		if (!chunk) { sstm << "Not enough columns in file " << famfile << " line " << nlines << " column 1."; die(sstm); }
		person[nlines].pid = (char *) malloc((strlen(chunk) + 1) * sizeof(char));
		if (!person[nlines].pid ) die("memory allocation error in person[nlines].pid");
		strcpy(person[nlines].pid, chunk);

        chunk = strtok(NULL, " \t");
		if (!chunk) { sstm << "Not enough columns in file " << famfile << " line " << nlines << " column 2."; die(sstm); }
		person[nlines].vid = (char *) malloc((strlen(chunk) + 1) * sizeof(char));
		if (!person[nlines].vid ) die("memory allocation error in person[nlines].vid");
		strcpy(person[nlines].vid, chunk);

        chunk = strtok(NULL, " \t");
		if (!chunk) { sstm << "Not enough columns in file " << famfile << " line " << nlines << " column 3."; die(sstm); }
		person[nlines].mid = (char *) malloc((strlen(chunk) + 1) * sizeof(char));
		if (!person[nlines].mid ) die("memory allocation error in person[nlines].mid");
		strcpy(person[nlines].mid, chunk);

        chunk = strtok(NULL, " \t");
		if (!chunk) { sstm << "Not enough columns in file " << famfile << " line " << nlines << " column 4."; die(sstm); }
		person[nlines].sex = atoi(chunk);

        chunk = strtok(NULL, " \t");
        if (!chunk) { sstm << "Not enough columns in file " << famfile << " line " << nlines << " column 5."; die(sstm); }
		if (!qt) {
			person[nlines].aff[thread] = atoi(chunk);
        } else {
            if (!strcmp("x", chunk) || !strcmp("-", chunk) ) {
                person[nlines].aff[thread] = 0;
                person[nlines].qtaff[thread] = 0;
                person[nlines].qcin = false;
                person[nlines].analysis_in = false;
            } else {
                person[nlines].aff[thread] = 2;
                person[nlines].qtaff[thread] = atof(chunk);
            }
        }

        nlines++;
    }
	fclose(fptr);

	if (s!=NULL) { free(s); s=NULL; }
	return nlinesfam;
}


uint32_t read_map(string mapfile, struct MAP*& map) {
    stringstream sstm;
	uint32_t nlinesmap = countLines(mapfile);

	map = (struct MAP*) calloc(nlinesmap, sizeof(struct MAP));
	if (!map) die("memory allocation error in map.");

    uint32_t nlines = 0;
    char* chunk = NULL;
    char* s = NULL;

    FILE* fptr = fopen(mapfile.c_str(), "r");
	for (uint32_t l=1; !feof(fptr); l++) {
        if (!getline(s, fptr)) continue;

		map[nlines].line = l;
		map[nlines].qcin = true;

		chunk = strtok(s, " \t");
		if (!chunk) { sstm << "Not enough columns in file " << mapfile << " line " << l << " column 0."; die(sstm); }
		if (strcmp(chunk, "X") == 0) {
			map[nlines].chr[0] = '2';
			map[nlines].chr[1] = '3';
			map[nlines].chr[2] = '\0';
		} else if (strcmp(chunk,"Y") == 0) {
			map[nlines].chr[0] = '2';
			map[nlines].chr[1] = '4';
			map[nlines].chr[2] = '\0';
		} else if (strcmp(chunk,"XY") == 0) {
			map[nlines].chr[0] = '2';
			map[nlines].chr[1] = '5';
			map[nlines].chr[2] = '\0';
		} else if (strcmp(chunk,"Mt") == 0 || strcmp(chunk,"MT") == 0) {
			map[nlines].chr[0] = '2';
			map[nlines].chr[1] = '6';
			map[nlines].chr[2] = '\0';
		} else {
			map[nlines].chr[0] = chunk[0];
			if (strlen(chunk) > 1) {
				map[nlines].chr[1] = chunk[1];
				map[nlines].chr[2] = '\0';
            } else if (strlen(chunk) == 1) {
                map[nlines].chr[1] = '\0';
            } else {
                sstm << "No chromosome found in file " << mapfile << " line " << nlines << ".";
                error(sstm);
                map[nlines].qcin = false;
            }
        }

		chunk = strtok(NULL, " \t");
		if (!chunk) { sstm << "Not enough columns in file " << mapfile << " line " << l << " column 1."; die(sstm); }

		map[nlines].rs = (char *) malloc((strlen(chunk) + 1) * sizeof(char));
		if (!map[nlines].rs) die("memory allocation error in map[ĺine].rs.");
		strcpy(map[nlines].rs, chunk);

		chunk = strtok(NULL, " \t");
		if (!chunk) { sstm << "Not enough columns in file " << mapfile << " line " << l << " column 2."; die(sstm); }
		map[nlines].dist = atoi(chunk);

		chunk = strtok(NULL, " \t");
		if (!chunk) { sstm << "Not enough columns in file " << mapfile << " line " << l << " column 3."; die(sstm); }
		map[nlines].pos = atoi(chunk);

        nlines++;
	}	fclose(fptr);

	if (s!=NULL) { free(s); s=NULL; }
	return nlinesmap;
}


uint32_t read_bmap(string bmapfile, struct MAP*& map, struct ALLELECODE*& codesA) {
    stringstream sstm;
	uint32_t nlinesbmap = countLines(bmapfile);
	map = (struct MAP*) calloc(nlinesbmap, sizeof(struct MAP));
	if (!map) die("memory allocation error in map.");

    codesA = (struct ALLELECODE*) malloc(nlinesbmap * sizeof(struct ALLELECODE));
    if (!codesA) die("memory allocation error in codesA.");
    memset(codesA, '0', nlinesbmap*sizeof(struct ALLELECODE));

    uint32_t nlines = 0;

    char* chunk = NULL;
    char* s = NULL;

    FILE* fptr = fopen(bmapfile.c_str(), "r");
	for (uint32_t l=1; !feof(fptr); l++) {
        if (!getline(s, fptr)) continue;
		map[nlines].line = l;
		map[nlines].qcin = true;

		chunk = strtok(s, " \t");
		if (!chunk) { sstm << "Not enough columns in file " << bmapfile << " line " << l << " column 0."; die(sstm); }
		if (strcmp(chunk, "X") == 0) {
			map[nlines].chr[0] = '2';
			map[nlines].chr[1] = '3';
			map[nlines].chr[2] = '\0';
		} else if (strcmp(chunk,"Y") == 0) {
			map[nlines].chr[0] = '2';
			map[nlines].chr[1] = '4';
			map[nlines].chr[2] = '\0';
		} else if (strcmp(chunk,"XY") == 0) {
			map[nlines].chr[0] = '2';
			map[nlines].chr[1] = '5';
			map[nlines].chr[2] = '\0';
		} else if (strcmp(chunk,"Mt") == 0 || strcmp(chunk,"MT") == 0) {
			map[nlines].chr[0] = '2';
			map[nlines].chr[1] = '6';
			map[nlines].chr[2] = '\0';
		} else {
			map[nlines].chr[0] = chunk[0];
			if (strlen(chunk) > 1) {
				map[nlines].chr[1] = chunk[1];
				map[nlines].chr[2] = '\0';
            } else if (strlen(chunk) == 1) {
                map[nlines].chr[1] = '\0';
            } else {
                sstm << "No chromosome found in file " << bmapfile << " line " << l << ".";
                error(sstm);
                map[nlines].qcin = false;
            }
        }

		chunk = strtok(NULL, " \t");
		if (!chunk) { sstm << "Not enough columns in file " << bmapfile << " line " << l << " column 1."; die(sstm); }
		map[nlines].rs = (char *) malloc((strlen(chunk) + 1) * sizeof(char));
		if (!map[nlines].rs) die("memory allocation error in map[ĺine].rs.");
		strcpy(map[nlines].rs, chunk);

		chunk = strtok(NULL, " \t");
		if (!chunk) { sstm << "Not enough columns in file " << bmapfile << " line " << l << " column 2."; die(sstm); }
		map[nlines].dist = atoi(chunk);

		chunk = strtok(NULL, " \t");
		if (!chunk) { sstm << "Not enough columns in file " << bmapfile << " line " << l << " column 3."; die(sstm); }
		map[nlines].pos = atoi(chunk);

		chunk = strtok(NULL, " \t");
		if (!chunk) { sstm << "Not enough columns in file " << bmapfile << " line " << l << " column 4."; die(sstm); }
		codesA[nlines].a1 = (char *) malloc((strlen(chunk) + 1) * sizeof(char));
		if (!codesA[nlines].a1) die("memory allocation error in codesA[nlines].a1.");
        strcpy(codesA[nlines].a1, chunk);

		chunk = strtok(NULL, " \t");
		if (!chunk) { sstm << "Not enough columns in file " << bmapfile << " line " << l << " column 5."; die(sstm); }
		codesA[nlines].a2 = (char *) malloc((strlen(chunk) + 1) * sizeof(char));
		if (!codesA[nlines].a2) die("memory allocation error in codesA[nlines].a2.");
        strcpy(codesA[nlines].a2, chunk);

        nlines++;
	}
	fclose(fptr);

	if (s!=NULL) { free(s); s=NULL; }
	return nlinesbmap;
}


uint32_t  read_ped(string  pedfile, struct MAP* helpmap, struct PERSON* helpperson, uint32_t nsnps2read, uint32_t nppls2read, struct ALLELECODE*& codesA, struct MAP*& map, struct PERSON*& person, uint64_t***& BinSNPs) {
    stringstream sstm;
    uint64_t ncolumnsped = countColumns(pedfile);

    uint32_t nwords = (uint32_t)ceil(nppls2read/64.0);
    BinSNPs = new uint64_t**[nsnps2read];
    if (!BinSNPs) die("memory allocation error in BinSNPs.");
    for (uint32_t w=0; w<nsnps2read; w++) {
        BinSNPs[w] = new uint64_t*[nwords];
        if (!BinSNPs[w]) die("memory allocation error in BinSNPs[w]");
        BinSNPs[w][0] = new uint64_t[nwords*4];
        if (!BinSNPs[w][0]) die("memory allocation error in BinSNPs[w][0]");
        memset(BinSNPs[w][0], 0, nwords*4*sizeof(uint64_t));
        for (uint32_t n=1; n<nwords; n++) BinSNPs[w][n] = BinSNPs[w][n-1]+4;
    }

	person = (struct PERSON*) calloc(nppls2read, sizeof(struct PERSON));
	if (!person) die("memory allocation error in person.");

	map = (struct MAP*) calloc(nsnps2read, sizeof(struct MAP));
	if (!map) die("memory allocation error in map.");

    codesA = (struct ALLELECODE*) calloc(nsnps2read,sizeof(struct ALLELECODE));
    if (!codesA) die("memory allocation error in codesA.");

	char* allele=NULL;
    uint32_t pos=0,nr=0;
    uint32_t ppl = 0;
    uint32_t snp = 0;
	uint32_t pplhelp = 0;
	uint32_t snphelp = 0;
    uint32_t alls0=0,alls1=0,alls3=0;
    bool cpy = true;
    char* chunk = NULL;
    char* s = NULL;

	allele=(char*)malloc(2*sizeof(char));
	if (!allele) die("memory allocation error in allele.");
	strcpy(allele,"0");

    for (uint32_t snp=0; snp<nsnps2read; snp++) {
		codesA[snp].a1=(char*)malloc(2*sizeof(char));
		if (!codesA[snp].a1) die("memory allocation error in codesA[snp].a1.");
		strcpy(codesA[snp].a1,"0");

		codesA[snp].a2=(char*)malloc(2*sizeof(char));
		if (!codesA[snp].a2) die("memory allocation error in codesA[snp].a2.");
		strcpy(codesA[snp].a2,"0");
	}
    FILE* fptr = fopen(pedfile.c_str(), "r");
	for (uint32_t l=0; !feof(fptr); l++,pplhelp++) {
        if (!getline(s, fptr)) { pplhelp--; continue; }
		if (!helpperson[pplhelp].analysis_in) continue;

		person[ppl++] = helpperson[pplhelp];
        if (pos==64) { pos=0; nr++; }

		chunk = strtok(s, " \t");
		chunk = strtok(NULL, " \t");
		chunk = strtok(NULL, " \t");
		chunk = strtok(NULL, " \t");
		chunk = strtok(NULL, " \t");
		chunk = strtok(NULL, " \t");

        snp = 0;
		for (uint32_t j=0; j<ncolumnsped-6; j++) {
			chunk = strtok(NULL, " \t");

			if (!chunk) { sstm << "Not enough columns in file " << pedfile << " line " << l << " column " << j+6 << "."; die(sstm); }
            snphelp = j/2;
            if (!helpmap[snphelp].in) continue;
            if (cpy && ~j&1) map[snp] = helpmap[snphelp];

            if (chunk[0] != '0' && (strcmp(codesA[snp].a1, chunk)!=0) && (strcmp(codesA[snp].a2, chunk)!=0)) {
                     if (!strcmp(codesA[snp].a1,"0")){
					 	codesA[snp].a1=(char*)realloc(codesA[snp].a1,(strlen(chunk)+1 )*sizeof(char));
						if (!codesA[snp].a1) die("memory allocation error in codesA[snp].a1.");
						strcpy(codesA[snp].a1, chunk);
						}
                else if (!strcmp(codesA[snp].a2,"0")) {
						codesA[snp].a2=(char*)realloc(codesA[snp].a2,(strlen(chunk)+1 )*sizeof(char));
						if (!codesA[snp].a2) die("memory allocation error in codesA[snp].a2.");
						strcpy(codesA[snp].a2, chunk);
						}
                else if (map[snp].qcin) {
                    sstm << "Line " << map[snp].line << " SNP " << map[snp].rs << " set to qc-out: more than two alleles.";
                    notice(sstm);
                    map[snp].qcin = false;
                    alls3++;
                }
            }

			if (j&1) {
                if ((chunk[0] == '0') ^ (allele[0] == '0')) {
                    sstm << "half-missing genotypes not allowed: " << map[snp].rs << ". Set to missing.";
                    error(sstm);
                }
                if (chunk[0] == '0' || allele[0] == '0') {
                    BinSNPs[snp][nr][0] |= mask[pos];
                } else if (strcmp(allele, chunk)!=0) {
                    BinSNPs[snp][nr][2] |= mask[pos];
                } else if (strcmp(allele, codesA[snp].a1)==0) {
                    BinSNPs[snp][nr][1] |= mask[pos];
                } else if (strcmp(allele, codesA[snp].a2)==0) {
                    BinSNPs[snp][nr][3] |= mask[pos];
                } else die("Odd things happening! Watch out for ghosts!!!");
                snp++;
            } else {
				allele=(char*)realloc(allele,(strlen(chunk)+1 )*sizeof(char));
				if (!allele) die("memory allocation error in allele.");
				strcpy(allele, chunk);
			}
		}

        pos++;
        cpy = false;
	}
	fclose(fptr);

    while (64>pos) {
        for (uint32_t snp=0; snp<nsnps2read; snp++)
            BinSNPs[snp][nr][0] |= mask[pos];
        pos++;
    }

    for (uint32_t snp=0; snp<nsnps2read; snp++) {
        if (codesA[snp].a1[0] == '0' && codesA[snp].a2[0] == '0') {
            sstm << "Map line " << map[snp].line << " SNP " << map[snp].rs << " set to qc-out: missing genotypes.";
            notice(sstm);
            map[snp].qcin = false;
            alls0++;
        } else if ((codesA[snp].a1[0]=='0') ^ (codesA[snp].a2[0]=='0')) {
            alls1++;
        }
    }

    if (alls3>0) { sstm << alls3 <<   " SNPs have more than two alleles and were set to qc-out."; warning(sstm); }
    if (alls0>0) { sstm << alls0 <<       " SNPs have missing genotypes and were set to qc-out."; warning(sstm); }
    if (alls1>0) { sstm << alls1 << " SNPs are monomorphous: optionally filter by defining MAF."; warning(sstm); }

	if (s!=NULL) { free(s); s=NULL; }
	if (allele!=NULL) { free(allele); allele=NULL; }
	return nwords;
}


uint32_t read_tped(string tpedfile, struct MAP* helpmap, struct PERSON* helpperson, uint32_t nsnps2read, uint32_t nppls2read, struct ALLELECODE*& codesA, struct MAP*& map, struct PERSON*& person, uint64_t***& BinSNPs) {
    stringstream sstm;
    uint64_t ncolumnstped = countColumns(tpedfile);

    uint32_t nwords = (uint32_t)ceil(nppls2read/64.0);
    BinSNPs = new uint64_t**[nsnps2read];
    if (!BinSNPs) die("memory allocation error in BinSNPs.");
    for (uint32_t w=0; w<nsnps2read; w++) {
        BinSNPs[w] = new uint64_t*[nwords];
        if (!BinSNPs[w]) die("memory allocation error in BinSNPs[w]");
        BinSNPs[w][0] = new uint64_t[nwords*4];
        if (!BinSNPs[w][0]) die("memory allocation error in BinSNPs[w][0]");
        memset(BinSNPs[w][0], 0, nwords*4*sizeof(uint64_t));
        for (uint32_t n=1; n<nwords; n++) BinSNPs[w][n] = BinSNPs[w][n-1]+4;
    }

	person = (struct PERSON*) calloc(nppls2read, sizeof(struct PERSON));
	if (!person) die("memory allocation error in person.");

	map = (struct MAP*) calloc(nsnps2read, sizeof(struct MAP));
	if (!map) die("memory allocation error in map.");

    codesA = (struct ALLELECODE*) calloc(nsnps2read,sizeof(struct ALLELECODE));
    if (!codesA) die("memory allocation error in codesA.");

	char* allele=NULL;
    int nalleles=0;
    uint32_t pos=0,nr=0;
    uint32_t snp = 0;
    uint32_t ppl = 0;
    uint32_t snphelp = 0;
    uint32_t pplhelp = 0;
    uint32_t alls0=0,alls1=0,alls3=0;
    bool cpy = true;
    char* chunk = NULL;
    char* s = NULL;

    FILE* fptr = fopen(tpedfile.c_str(), "r");
	for (uint32_t l=1; !feof(fptr); l++,snphelp++) {
        if (!getline(s, fptr)) { snphelp--; continue; }
        if (!helpmap[snphelp].in) continue;

        map[snp] = helpmap[snphelp];

		chunk = strtok(s, " \t");
		chunk = strtok(NULL, " \t");
		chunk = strtok(NULL, " \t");
		chunk = strtok(NULL, " \t");

        nalleles = 0;
		pos=0, nr=0;

		codesA[snp].a1=(char*)malloc(2*sizeof(char));
		if (!codesA[snp].a1) die("memory allocation error in codesA[snp].a1.");
		strcpy(codesA[snp].a1,"0");

		codesA[snp].a2=(char*)malloc(2*sizeof(char));
		if (!codesA[snp].a2) die("memory allocation error in codesA[snp].a2.");
		strcpy(codesA[snp].a2,"0");

		allele=(char*)malloc(2*sizeof(char));
		if (!allele) die("memory allocation error in allele.");
		strcpy(allele,"0");

		for (uint32_t j=0; j<ncolumnstped-4; j++) {
			chunk = strtok(NULL, " \t");

			if (!chunk) { sstm << "Not enough columns in file " << tpedfile << " line " << l << " column " << j+4 << "."; die(sstm); }
			pplhelp = j/2;
			if (!helpperson[pplhelp].analysis_in) continue;
			if (cpy && ~j&1) person[ppl++] = helpperson[pplhelp];
			if (chunk[0] != '0' && (strcmp(codesA[snp].a1, chunk)!=0) && (strcmp(codesA[snp].a2, chunk)!=0)) {
                nalleles++;
                if (nalleles <= 2) {
                     if (!strcmp(codesA[snp].a1,"0"))
					 {
						codesA[snp].a1=(char*)realloc(codesA[snp].a1,(strlen(chunk)+1 )*sizeof(char));
						if (!codesA[snp].a1) die("memory allocation error in codesA[snp].a1.");
						strcpy(codesA[snp].a1, chunk);
					 }
                     else
					 {
						codesA[snp].a2=(char*)realloc(codesA[snp].a2,(strlen(chunk)+1 )*sizeof(char));
						if (!codesA[snp].a2) die("memory allocation error in codesA[snp].a2.");
						strcpy(codesA[snp].a2, chunk);
					 }
                }
            }

			if (j&1) {
                if (pos==64) { pos=0; nr++; }

                if ((chunk[0] == '0') ^ (allele[0] == '0')) {
                    sstm << "half-missing genotypes not allowed: " << map[snp].rs << ". Set to missing.";
                    error(sstm);
                }

                if (chunk[0] == '0' || allele[0] == '0') {
                    BinSNPs[snp][nr][0] |= mask[pos];
                } else if (strcmp(allele, chunk)!=0) {
                    BinSNPs[snp][nr][2] |= mask[pos];
                } else if (strcmp(allele, codesA[snp].a1)==0) {
                    BinSNPs[snp][nr][1] |= mask[pos];
                } else if (strcmp(allele, codesA[snp].a2)==0) {
                    BinSNPs[snp][nr][3] |= mask[pos];
                } else die("Odd things happening! Watch out for ghosts!!!");
                pos++;
            } else {
				allele=(char*)realloc(allele,(strlen(chunk)+1 )*sizeof(char));
				if (!allele) die("memory allocation error in allele.");
				strcpy(allele, chunk);
			}
		}
        while (64>pos) BinSNPs[snp][nr][0] |= mask[pos++];

        if (nalleles > 2) {
            sstm << "Line " << map[snp].line << " SNP " << map[snp].rs << " set to qc-out: more than two alleles.";
            notice(sstm);
            map[snp].qcin = false;
            alls3++;
        } else if (nalleles == 0) {
            sstm << "Line " << map[snp].line << " SNP " << map[snp].rs << " set to qc-out: missing genotypes.";
            notice(sstm);
            map[snp].qcin = false;
            alls0++;
        } else if (nalleles == 1) {
            alls1++;
        }

        snp++;
        cpy = false;
	}
	fclose(fptr);

    if (alls3>0) { sstm << alls3 <<   " SNPs have more than two alleles and were set to qc-out."; warning(sstm); }
    if (alls0>0) { sstm << alls0 <<       " SNPs have missing genotypes and were set to qc-out."; warning(sstm); }
    if (alls1>0) { sstm << alls1 << " SNPs are monomorphous: optionally filter by defining MAF."; warning(sstm); }

	if (s!=NULL) { free(s); s=NULL; }
	if (allele!=NULL) { free(allele); allele=NULL; }
	return nwords;
}


uint32_t read_bped(string bpedfile, struct MAP* helpmap, struct PERSON* helpperson, uint32_t nsnps2read, uint32_t nppls2read, uint32_t nlinesfam, uint32_t nlinestped, struct ALLELECODE*& codesA, struct MAP*& map, struct PERSON*& person, uint64_t***& BinSNPs) {
    stringstream sstm;

    uint32_t nwords = (uint32_t)ceil(nlinesfam/64.0);
    BinSNPs = new uint64_t**[nsnps2read];
    if (!BinSNPs) die("memory allocation error in BinSNPs.");
    for (uint32_t w=0; w<nsnps2read; w++) {
        BinSNPs[w] = new uint64_t*[nwords];
        if (!BinSNPs[w]) die("memory allocation error in BinSNPs[w]");
        BinSNPs[w][0] = new uint64_t[nwords*4];
        if (!BinSNPs[w][0]) die("memory allocation error in BinSNPs[w][0]");
        memset(BinSNPs[w][0], 0, nwords*4*sizeof(uint64_t));
        for (uint32_t n=1; n<nwords; n++) BinSNPs[w][n] = BinSNPs[w][n-1]+4;
    }

	person = (struct PERSON*) calloc(nppls2read, sizeof(struct PERSON));
	if (!person) die("memory allocation error in person.");

	map = (struct MAP*) calloc(nsnps2read, sizeof(struct MAP));
	if (!map) die("memory allocation error for map in read_bped().");

	struct ALLELECODE* codesB = (struct ALLELECODE*) calloc(nsnps2read, sizeof(struct ALLELECODE));
    if (!codesB) die("memory allocation error for codesA in read_bped().");

    const uint8_t list[4] = {1,0,2,3};
    uint16_t magicnr;
    uint8_t mode;
    uint32_t pos=0,nr=0;
    uint32_t snp=0, snphelp=0;
    uint32_t alls0=0,alls1=0;
    size_t size;
    uint8_t* s=NULL;

    FILE* fptr = fopen(bpedfile.c_str(), "rb");
    size = fread(&magicnr, 2, 1, fptr);
    if (magicnr!=0x1B6C && magicnr!=0x6C1B) { sstm << "File " << bpedfile << " is not recognized as PLINK v1.00 BPED file."; die(sstm); }
    size = fread(&mode, 1, 1, fptr);
    if (mode == 1) {
        uint32_t buffsiz = ceil(nlinesfam/4.0);
        s = (uint8_t*)malloc(buffsiz);
        while (!feof(fptr)) {
            if (!fread(s, 1, buffsiz, fptr)) break;
            if (!helpmap[snphelp].in) { snphelp++; continue; }
            map[snp] = helpmap[snphelp];
            codesB[snp] = codesA[snphelp];

            for (uint32_t b=0; b<buffsiz; b++) {
                if (pos==64) { nr++; pos=0; }
                BinSNPs[snp][nr][list[ s[b]     & 0x03]] |= mask[pos++];
                BinSNPs[snp][nr][list[(s[b]>>2) & 0x03]] |= mask[pos++];
                BinSNPs[snp][nr][list[(s[b]>>4) & 0x03]] |= mask[pos++];
                BinSNPs[snp][nr][list[(s[b]>>6) & 0x03]] |= mask[pos++];
            }
            if (uint8_t n=nlinesfam%4)
                for (uint8_t i=4-n; i>0; i--)
                    setGenotype(BinSNPs[snp][nr], pos-i, 0);
            while (64>pos) BinSNPs[snp][nr][0] |= mask[pos++];


            if ((codesB[snp].a1[0] == '0') ^ (codesB[snp].a2[0] == '0')) {
                alls1++;
            } else if (codesB[snp].a1[0] == '0' && codesB[snp].a2[0] == '0') {
                alls0++;
                sstm << "Line " << map[snp].line << " SNP " << map[snp].rs << " set to qc-out: missing genotypes.";
                notice(sstm);
                map[snp].qcin = false;
            }

            snphelp++;
            snp++;
            pos=0;
            nr=0;
        }
    } else if (mode == 0) {
        sstm << "Binary ped file " << bpedfile << " is in individual-major mode. Convert to SNP-major mode first.";
        die(sstm);
    } else die("Unknown binary file mode.");
	fclose(fptr);

    if (alls0>0) { sstm << alls0 <<       " SNPs have missing genotypes and were set to qc-out."; warning(sstm); }
    if (alls1>0) { sstm << alls1 << " SNPs are monomorphous: optionally filter by defining MAF."; warning(sstm); }

    if (nppls2read != nlinesfam) {
        uint64_t* AnalysisFlags = getAnalysisIn_bin_snp_gen(helpperson, nlinesfam, nwords, false);
        nwords = reduce_bin_gen(BinSNPs, nsnps2read, nwords, AnalysisFlags);
        delete[] AnalysisFlags; AnalysisFlags=NULL;
    }
    for (uint32_t i=0,j=0; i<nlinesfam; i++)
        if (helpperson[i].analysis_in)
            person[j++] = helpperson[i];
	for(uint32_t i = 0; i< nlinestped;i++){
		if(!helpmap[i].in){
			free(codesA[i].a1); codesA[i].a1 = NULL;
			free(codesA[i].a2); codesA[i].a2 = NULL;
		}
	}

    free(codesA);
    codesA = codesB;
    codesB = NULL;
	if (s!=NULL) { free(s); s=NULL; }
	return nwords;
}


uint32_t read_impute2(string impute2file, struct MAP* helpmap, struct PERSON* helpperson, uint32_t nsnps2read, uint32_t nppls2read, struct ALLELECODE*& codesA, struct MAP*& map, struct PERSON*& person, uint64_t***& BinSNPs, double probRequired, float*** genoWeights, int dosage) {
    stringstream sstm;
    uint64_t ncolumnsimpute2 = countColumns(impute2file);

    uint32_t nwords = (uint32_t)ceil(nppls2read/64.0);
    BinSNPs = new uint64_t**[nsnps2read];
    if (!BinSNPs) die("memory allocation error in BinSNPs.");
    for (uint32_t w=0; w<nsnps2read; w++) {
        BinSNPs[w] = new uint64_t*[nwords];
        if (!BinSNPs[w]) die("memory allocation error in BinSNPs[w]");
        BinSNPs[w][0] = new uint64_t[nwords*4];
        if (!BinSNPs[w][0]) die("memory allocation error in BinSNPs[w][0]");
        memset(BinSNPs[w][0], 0, nwords*4*sizeof(uint64_t));
        for (uint32_t n=1; n<nwords; n++) BinSNPs[w][n] = BinSNPs[w][n-1]+4;
    }

	person = (struct PERSON*) calloc(nppls2read, sizeof(struct PERSON));
	if (!person) die("memory allocation error in person.");

	map = (struct MAP*) calloc(nsnps2read, sizeof(struct MAP));
	if (!map) die("memory allocation error in map.");

    codesA = (struct ALLELECODE*) calloc(nsnps2read, sizeof(struct ALLELECODE));
    if (!codesA) die("memory allocation error in codesA.");

    uint32_t pos=0,nr=0;
    uint32_t snp = 0;
    uint32_t ppl = 0;
    uint32_t snphelp = 0;
    uint32_t pplhelp = 0;
    uint32_t alls1=0;
	uint32_t missing=0, homozyg1=0,heterozyg=0,homozyg2=0;
	uint32_t k =0;
    char* chunk = NULL;
    char* s = NULL;
	double g[3];
	double maximum;

    FILE* fptr = fopen(impute2file.c_str(), "r");
	for (uint32_t l=1; !feof(fptr); l++,snphelp++) {
        if (!getline(s, fptr)) { snphelp--; continue; }
        if (!helpmap[snphelp].in) continue;
        map[snp] = helpmap[snphelp];

		chunk = strtok(s, " \t");
		chunk = strtok(NULL, " \t");
		chunk = strtok(NULL, " \t");
		chunk = strtok(NULL, " \t");

		codesA[snp].a1=(char*)malloc((strlen(chunk)+1 )*sizeof(char));
		if (!codesA[snp].a1) die("memory allocation error in codesA[snp].a1.");
		strcpy(codesA[snp].a1, chunk);

		chunk = strtok(NULL, " \t");
		codesA[snp].a2=(char*)malloc((strlen(chunk)+1 )*sizeof(char));
		if (!codesA[snp].a2) die("memory allocation error in codesA[snp].a2.");
		strcpy(codesA[snp].a2, chunk);

		if ((codesA[snp].a1[0] == '0') ^ (codesA[snp].a2[0] == '0')) {
			sstm << "SNP " << map[snp].rs << " has one allele 0.";
            warning(sstm);
        }
		if ((codesA[snp].a1[0] == '0') && (codesA[snp].a2[0] == '0')) {
			sstm << "SNP " << map[snp].rs << " has both alleles 0.";
            warning(sstm);
        }

		pos=0, nr=0;
		missing=0, homozyg1=0,heterozyg=0,homozyg2=0;
		for (uint32_t j=0; j<ncolumnsimpute2-5; j++) {
			chunk = strtok(NULL, " \t");
			g[k]=atof(chunk);
			if (!chunk) { sstm << "Not enough columns in file " << impute2file << " line " << l << " column " << j+5 << "."; die(sstm); }
			pplhelp = j/3;
			if (!helpperson[pplhelp].analysis_in) continue;
			if(dosage==1 && ppl<nppls2read)
			{
			genoWeights[snp][ppl][k] = g[k];
			}
			k++;
			if((j+1)%3==0){
				person[ppl++] = helpperson[pplhelp];
                if (pos==64) { pos=0; nr++; }
				maximum = g[0];
				if(g[1] > maximum){maximum = g[1];}
				if(g[2] > maximum){maximum = g[2];}
			    if( maximum == 0 || maximum < probRequired)
				{
					BinSNPs[snp][nr][0] |= mask[pos];
					missing++;
				}
				else if(maximum == g[0])
				{
					BinSNPs[snp][nr][1] |= mask[pos];
					homozyg1++;
				}
				else if(maximum == g[1])
				{
					 BinSNPs[snp][nr][2] |= mask[pos];
					 heterozyg++;
				}
			    else if(maximum == g[2])
			    {
					BinSNPs[snp][nr][3] |= mask[pos];
					homozyg2++;
				}
				k=0; pos++;
			}
		}

        while (64>pos) BinSNPs[snp][nr][0] |= mask[pos++];

		if ((homozyg1!=0 && heterozyg==0 && homozyg2==0) || (homozyg1==0 && heterozyg==0 && homozyg2!=0) || (homozyg1==0 && heterozyg==0 && homozyg2==0)) {
            alls1++;
			strcpy(codesA[snp].a2,"0");
        }
        snp++; ppl=0;
	}
	fclose(fptr);

    if (alls1>0) { sstm << alls1 << " SNPs are monomorphous: optionally filter by defining MAF."; warning(sstm); }

	if (s!=NULL) { free(s); s=NULL; }
	return nwords;
}

