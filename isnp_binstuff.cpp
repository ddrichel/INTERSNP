/**~ isnp_binstuff.cpp
.---------------------------------------------------------------------------.
|  Software: INTERSNP, Genome-Wide Interaction Analysis                     |
|      Site: http://intersnp.meb.uni-bonn.de/                               |
| ------------------------------------------------------------------------- |
|      File: isnp_binstuff.cpp                                              |
|    Author: André Lacour                                                   |
|   Content: source-code library on binary manipulations for INTERSNP       |
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


double    (*getIbs)(uint64_t***, uint32_t, uint32_t, uint32_t, uint64_t*);
double  getTotalIbs(uint64_t***, uint32_t, uint32_t, uint32_t, uint64_t* =NULL);
double  getMinorIbs(uint64_t***, uint32_t, uint32_t, uint32_t, uint64_t* =NULL);
float**   (*getIbsMatrix)(uint64_t***, uint32_t, uint32_t);
float** getTotalIbsMatrix(uint64_t***, uint32_t, uint32_t);
float** getMinorIbsMatrix(uint64_t***, uint32_t, uint32_t);
uint32_t transpose_bin_gen(uint64_t***, uint32_t, uint32_t, uint64_t***&);
uint32_t    reduce_bin_gen(uint64_t***&, uint32_t, uint32_t, uint64_t*);
void chr_bin_ind_gen(uint64_t***, int*, struct MAP*, uint32_t, uint32_t, uint64_t****&, uint32_t*&);
uint64_t*       getQCIn_bin_ind_gen(struct MAP*, uint32_t, uint32_t, uint32_t);
uint64_t* getAnalysisIn_bin_ind_gen(struct MAP*, uint32_t, uint32_t, uint32_t);
uint64_t* getMatchingIn_bin_ind_gen(struct MAP*, uint32_t, uint32_t, uint32_t);
uint64_t*     getRareIn_bin_ind_gen(struct MAP*, uint32_t, uint32_t, uint32_t, double);
uint64_t*       getQCIn_bin_snp_gen(struct PERSON*, uint32_t, uint32_t);
uint64_t*     getQTaff2_bin_snp_gen(struct PERSON*, uint32_t, uint32_t, bool =true);
uint64_t* getAnalysisIn_bin_snp_gen(struct PERSON*, uint32_t, uint32_t, bool =true);
uint64_t**    getGender_bin_snp_gen(struct PERSON*, uint32_t, uint32_t, bool =true);
uint64_t**        getCC_bin_snp_gen(struct PERSON*, uint32_t, uint32_t, bool =true);
uint64_t**     updateCC_bin_snp(uint64_t**, struct PERSON*, uint32_t, uint32_t, bool =true);
uint64_t**    getCovCat_bin_snp_gen(struct PERSON*, uint32_t, uint32_t, int, bool =true);
struct PPLLOCATION* getPplLocations(struct PERSON*, uint32_t, uint32_t, bool =true);
uint8_t   getGenotype(uint64_t*, uint8_t);
uint64_t* setGenotype(uint64_t*, uint8_t, uint8_t);
uint16_t bitcount64(uint64_t);
bool getbit64(uint64_t, uint8_t);
void setbit64(uint64_t&, uint8_t);
void delbit64(uint64_t&, uint8_t);
template <typename T>string T2bin(T);
template <typename T>void delete_2dim(T**&);
template <typename T>void delete_3dim(T***&, uint32_t);
template <typename T>void delete_4dim(T****&, uint32_t);
template <typename T>void printmatrix(T**, uint32_t, uint32_t, string);
template <typename T>void printmatrix(T**, uint32_t, string);

void die(string);
void dev(string);
void dev(stringstream&);


struct PPLLOCATION {
    int32_t nr;
    int8_t pos;
};


const uint64_t mask[64] = {0x8000000000000000ull,0x4000000000000000ull,0x2000000000000000ull,0x1000000000000000ull,0x800000000000000ull,0x400000000000000ull,0x200000000000000ull,0x100000000000000ull,0x80000000000000ull,0x40000000000000ull,0x20000000000000ull,0x10000000000000ull,0x8000000000000ull,0x4000000000000ull,0x2000000000000ull,0x1000000000000ull,0x800000000000ull,0x400000000000ull,0x200000000000ull,0x100000000000ull,0x80000000000ull,0x40000000000ull,0x20000000000ull,0x10000000000ull,0x8000000000ull,0x4000000000ull,0x2000000000ull,0x1000000000ull,0x800000000ull,0x400000000ull,0x200000000ull,0x100000000ull,0x80000000ull,0x40000000ull,0x20000000ull,0x10000000ull,0x8000000ull,0x4000000ull,0x2000000ull,0x1000000ull,0x800000ull,0x400000ull,0x200000ull,0x100000ull,0x80000ull,0x40000ull,0x20000ull,0x10000ull,0x8000ull,0x4000ull,0x2000ull,0x1000ull,0x800ull,0x400ull,0x200ull,0x100ull,0x80ull,0x40ull,0x20ull,0x10ull,0x8ull,0x4ull,0x2ull,0x1ull};



double getTotalIbs(uint64_t*** People, uint32_t i, uint32_t j, uint32_t nwords, uint64_t* InFlags) {
    uint64_t equalAllels = 0;
    uint64_t   allAllels = 0;
    register uint64_t *X,*Y,Z;
    Z = 0xFFFFFFFFFFFFFFFFull;
    for (uint32_t k=0; k<nwords; k++) {
        X = People[i][k];
        Y = People[j][k];
        if (InFlags!=NULL) Z = InFlags[k];
        equalAllels += 2*bitcount64(((X[1] & Y[1]) | (X[2] & Y[2]) | (X[3] & Y[3])) & Z);
        equalAllels +=   bitcount64(((X[1] & Y[2]) | (X[2] & Y[1]) | (X[2] & Y[3]) | (X[3] & Y[2])) & Z);
          allAllels +=   bitcount64(~(X[0] | Y[0]) & Z);
    }
    return (double)equalAllels/(2.0*allAllels);
}


float** getTotalIbsMatrix(uint64_t*** People, uint32_t npeople, uint32_t nwords) {
    float** Ibs = new float*[npeople];
    if (!Ibs) die("memory allocation error for Ibs in getTotalIbs()");
    Ibs[0] = new float[npeople*npeople];
    if (!Ibs[0]) die("memory allocation error for Ibs[0] in getTotalIbs()");
    memset(Ibs[0], 0, npeople*npeople*sizeof(float));
    for (uint32_t i=1; i<npeople; i++) Ibs[i] = Ibs[i-1]+npeople;
    for (uint32_t i=0; i<npeople; i++) {
        Ibs[i][i] = 1.0;
        for (uint32_t j=i+1; j<npeople; j++) {
            Ibs[i][j] = getTotalIbs(People, i, j, nwords);
            Ibs[j][i] = Ibs[i][j];
        }
    }
    return Ibs;
}


double getMinorIbs(uint64_t*** People, uint32_t i, uint32_t j, uint32_t nwords, uint64_t* InFlags) {
    uint64_t equalAllels = 0;
    uint64_t   allAllels = 0;
    register uint64_t *X,*Y,Z;
    Z = 0xFFFFFFFFFFFFFFFFull;
    for (uint32_t k=0; k<nwords; k++) {
        X = People[i][k];
        Y = People[j][k];
        if (InFlags!=NULL) Z = InFlags[k];
        equalAllels += 2*bitcount64(((X[1] & Y[1])) & Z);
        equalAllels +=   bitcount64(((X[1] & Y[2]) | (X[2] & Y[1]) | (X[2] & Y[2])) & Z);
          allAllels +=   bitcount64(~(X[0] | Y[0]) & Z);
    }
    return (double)equalAllels/(2.0*allAllels);
}


float** getMinorIbsMatrix(uint64_t*** People, uint32_t npeople, uint32_t nwords) {
    float** Ibs = new float*[npeople];
    if (!Ibs) die("memory allocation error for Ibs in getMinorIbs()");
    Ibs[0] = new float[npeople*npeople];
    if (!Ibs[0]) die("memory allocation error for Ibs[0] in getMinorIbs()");
    memset(Ibs[0], 0, npeople*npeople*sizeof(float));
    for (uint32_t i=1; i<npeople; i++) Ibs[i] = Ibs[i-1]+npeople;
    for (uint32_t i=0; i<npeople; i++) {
        Ibs[i][i] = 1.0;
        for (uint32_t j=i+1; j<npeople; j++) {
            Ibs[i][j] = getMinorIbs(People, i, j, nwords);
            Ibs[j][i] = Ibs[i][j];
        }
    }
    return Ibs;
}


double getPartialIbs(uint64_t*** People, uint32_t i, uint32_t j, uint32_t a, uint32_t b, uint64_t* InFlags) {
    uint64_t equalAllels = 0;
    uint64_t   allAllels = 0;
    register uint64_t *X,*Y,Z;
    Z = 0xFFFFFFFFFFFFFFFFull;
    for (uint32_t k=a; k<=b; k++) {
        X = People[i][k];
        Y = People[j][k];
        if (InFlags!=NULL) Z = InFlags[k];
        equalAllels += 2*bitcount64(((X[1] & Y[1]) | (X[2] & Y[2]) | (X[3] & Y[3])) & Z);
        equalAllels +=   bitcount64(((X[1] & Y[2]) | (X[2] & Y[1]) | (X[2] & Y[3]) | (X[3] & Y[2])) & Z);
          allAllels +=   bitcount64(~(X[0] | Y[0]) & Z);
    }
    return (double)equalAllels/(2.0*allAllels);
}


uint32_t transpose_bin_gen(uint64_t*** As, uint32_t nAqc, uint32_t nBqc, uint64_t***& Bs) {
    uint64_t** buf = NULL;
    uint32_t nwordsB = (uint32_t)ceil(nAqc/64.0);
    Bs = new uint64_t**[nBqc];
    for (uint32_t nrB=0, posB=0, b=0; b<nBqc; b++, posB++) {
        if (posB==64) { posB=0; nrB++; }
        buf = new uint64_t*[nwordsB];
        if (!buf) die("memory allocation error for buf in transpose_bin_gen()");
        buf[0] = new uint64_t[nwordsB*4];
        if (!buf[0]) die("memory allocation error for buf[0] in transpose_bin_gen()");
        memset(buf[0], 0, nwordsB*4*sizeof(uint64_t));
        for (uint32_t n=1; n<nwordsB; n++) buf[n] = buf[n-1]+4;
        uint32_t nrA=0, posA=0;
        for (uint32_t a=0; a<nAqc; a++) {
            if (posA==64) { posA=0; nrA++; }
            buf[nrA][getGenotype(As[a][nrB], posB)] |= mask[posA++];
        }
        while (64>posA) buf[nrA][0] |= mask[posA++];
        Bs[b] = buf;
    }
    return nwordsB;
}


uint32_t reduce_bin_gen(uint64_t***& As, uint32_t nAqc, uint32_t nwordsA, uint64_t* setFlags) {
    uint32_t nsetBits = 0;
    for (uint32_t k=0; k<nwordsA; k++) nsetBits += bitcount64(setFlags[k]);
    uint64_t** buf = NULL;
    uint32_t nwordsB = (uint32_t)ceil(nsetBits/64.0);
    for (uint32_t row=0; row<nAqc; row++) {
        buf = new uint64_t*[nwordsB];
        if (!buf) die("memory allocation error for buf in reduce_bin_gen()");
        buf[0] = new uint64_t[nwordsB*4];
        if (!buf[0]) die("memory allocation error for buf[0] in reduce_bin_gen()");
        for (uint32_t n=1; n<nwordsB; n++) buf[n] = buf[n-1]+4;
        memset(buf[0], 0, nwordsB*4*sizeof(uint64_t));
        uint32_t nrB=0, posB=0;
        for (uint32_t nrA=0, posA=0; nrA<nwordsA-1 || posA<64; posA++) {
            if (posA==64) { posA=0; nrA++; }
            if (setFlags[nrA] & mask[posA]) {
                if (posB==64) { posB=0; nrB++; }
                buf[nrB][getGenotype(As[row][nrA], posA)] |= mask[posB++];
            }
        }
        while (64>posB) buf[nrB][0] |= mask[posB++];
        delete_2dim(As[row]);
        As[row] = buf;
        buf = NULL;
    }
    return nwordsB;
}


uint64_t* getQCIn_bin_ind_gen(struct MAP* map, uint32_t begin, uint32_t length, uint32_t nwords) {
    if (!length || !nwords) return NULL;
    uint64_t* QCInFlags = new uint64_t[nwords];
    memset(QCInFlags, 0, nwords*sizeof(uint64_t));
    for (uint32_t nr=0, pos=0, w=begin; w<begin+length; w++) {
        if (pos==64) { pos=0; nr++; }
        if (map[w].qcin) QCInFlags[nr] |= mask[pos];
        pos++;
    }
    return QCInFlags;
}


uint64_t* getAnalysisIn_bin_ind_gen(struct MAP* map, uint32_t begin, uint32_t length, uint32_t nwords) {
    if (!length || !nwords) return NULL;
    uint64_t* AnalysisInFlags = new uint64_t[nwords];
    memset(AnalysisInFlags, 0, nwords*sizeof(uint64_t));
    for (uint32_t nr=0, pos=0, w=begin; w<begin+length; w++) {
        if (!map[w].qcin) continue;
        if (pos==64) { pos=0; nr++; }
        if (map[w].analysis_in) AnalysisInFlags[nr] |= mask[pos];
        pos++;
    }
    return AnalysisInFlags;
}


uint64_t* getMatchingIn_bin_ind_gen(struct MAP* map, uint32_t begin, uint32_t length, uint32_t nwords) {
    if (!length || !nwords) return NULL;
    uint64_t* MatchingInFlags = new uint64_t[nwords];
    memset(MatchingInFlags, 0, nwords*sizeof(uint64_t));
    for (uint32_t nr=0, pos=0, w=begin; w<begin+length; w++) {
        if (!map[w].qcin) continue;
        if (pos==64) { pos=0; nr++; }
        if (map[w].matching_in) MatchingInFlags[nr] |= mask[pos];
        pos++;
    }
    return MatchingInFlags;
}


uint64_t* getRareIn_bin_ind_gen(struct MAP* map, uint32_t begin, uint32_t length, uint32_t nwords, double rarefreq) {
    if (!length || !nwords) return NULL;
    uint64_t* RareInFlags = new uint64_t[nwords];
    memset(RareInFlags, 0, nwords*sizeof(uint64_t));
    for (uint32_t nr=0, pos=0, w=begin; w<begin+length; w++) {
        if (!map[w].qcin) continue;
        if (pos==64) { pos=0; nr++; }
        if ( rarefreq+EPS>map[w].maf && map[w].maf>EPS ) RareInFlags[nr] |= mask[pos];
        pos++;
    }
    return RareInFlags;
}


uint64_t* getQCIn_bin_snp_gen(struct PERSON* person, uint32_t nlinestfam, uint32_t nwords) {
    uint64_t* QCInFlags = new uint64_t[nwords];
    memset(QCInFlags, 0, nwords*sizeof(uint64_t));
    for (uint32_t nr=0, pos=0, k=0; k<nlinestfam; k++) {
        if (pos==64) { pos=0; nr++; }
        if (person[k].qcin) QCInFlags[nr] |= mask[pos];
        pos++;
    }
    return QCInFlags;
}


uint64_t* getQTaff2_bin_snp_gen(struct PERSON* person, uint32_t nlinestfam, uint32_t nwords, bool qcflag) {
    #if PARALLEL
    uint32_t thread = omp_get_thread_num();
    #else
    uint32_t thread = 0;
    #endif
    uint64_t* QTaff2Flags = new uint64_t[nwords];
    memset(QTaff2Flags, 0, nwords*sizeof(uint64_t));
    for (uint32_t nr=0, pos=0, k=0; k<nlinestfam; k++) {
        if (qcflag && !person[k].qcin) continue;
        if (pos==64) { pos=0; nr++; }
        if (fabs(person[k].qtaff[thread]-2) < EPS) QTaff2Flags[nr] |= mask[pos];
        pos++;
    }
    return QTaff2Flags;
}


uint64_t* getAnalysisIn_bin_snp_gen(struct PERSON* person, uint32_t nlinestfam, uint32_t nwords, bool qcflag) {
    uint64_t* AnalysisInFlags = new uint64_t[nwords];
    memset(AnalysisInFlags, 0, nwords*sizeof(uint64_t));
    for (uint32_t nr=0, pos=0, k=0; k<nlinestfam; k++) {
        if (qcflag && !person[k].qcin) continue;
        if (pos==64) { pos=0; nr++; }
        if (person[k].analysis_in) AnalysisInFlags[nr] |= mask[pos];
        pos++;
    }
    return AnalysisInFlags;
}


uint64_t** getGender_bin_snp_gen(struct PERSON* person, uint32_t nlinestfam, uint32_t nwords, bool qcflag) {
    uint64_t** GenderFlags = new uint64_t*[nwords];
    if (!GenderFlags) die("memory allocation error for GenderFlags in getGender_bin_snp_gen()");
    GenderFlags[0] = new uint64_t[nwords*3];
    if (!GenderFlags[0]) die("memory allocation error for GenderFlags[0] in getGender_bin_snp_gen()");
    memset(GenderFlags[0], 0, nwords*3*sizeof(uint64_t));
    for (uint32_t n=1; n<nwords; n++) {
        GenderFlags[n] = GenderFlags[n-1]+3;
    }
    uint32_t nr=0, pos=0;
    for (uint32_t k=0; k<nlinestfam; k++) {
        if (qcflag && !person[k].qcin) continue;
        if (pos==64) { pos=0; nr++; }
        GenderFlags[nr][person[k].sex] |= mask[pos];
        pos++;
    }
    while (64>pos) GenderFlags[nr][0] |= mask[pos++];
    return GenderFlags;
}


uint64_t** getCC_bin_snp_gen(struct PERSON* person, uint32_t nlinestfam, uint32_t nwords, bool qcflag) {
    #if PARALLEL
    uint32_t thread = omp_get_thread_num();
    #else
    uint32_t thread = 0;
    #endif
    uint64_t** CCFlags = new uint64_t*[nwords];
    if (!CCFlags) die("memory allocation error for CCFlags in getCC_bin_snp_gen()");
    CCFlags[0] = new uint64_t[nwords*3];
    if (!CCFlags[0]) die("memory allocation error for CCFlags[0] in getCC_bin_snp_gen()");
    memset(CCFlags[0], 0, nwords*3*sizeof(uint64_t));
    for (uint32_t n=1; n<nwords; n++) CCFlags[n] = CCFlags[n-1]+3;
    uint32_t nr=0, pos=0;
    for (uint32_t k=0; k<nlinestfam; k++) {
        if (qcflag && !person[k].qcin) continue;
        if (pos==64) { pos=0; nr++; }
        CCFlags[nr][person[k].aff[thread]] |= mask[pos];
        pos++;
    }
    while (64>pos) CCFlags[nr][0] |= mask[pos++];
    return CCFlags;
}


uint64_t** updateCC_bin_snp(uint64_t** CCFlags, struct PERSON* person, uint32_t nlinestfam, uint32_t nwords, bool qcflag) {
    #if PARALLEL
    uint32_t thread = omp_get_thread_num();
    #else
    uint32_t thread = 0;
    #endif
    if (CCFlags==NULL) return getCC_bin_snp_gen(person, nlinestfam, nwords, qcflag);
    memset(CCFlags[0], 0, nwords*3*sizeof(uint64_t));
    uint32_t nr=0, pos=0;
    for (uint32_t k=0; k<nlinestfam; k++) {
        if (qcflag && !person[k].qcin) continue;
        if (pos==64) { pos=0; nr++; }
        CCFlags[nr][person[k].aff[thread]] |= mask[pos];
        pos++;
    }
    while (64>pos) CCFlags[nr][0] |= mask[pos++];
    return CCFlags;
}


uint64_t** getCovCat_bin_snp_gen(struct PERSON* person, uint32_t nlinestfam, uint32_t nwords, int nCovCathegories, bool qcflag) {
  uint64_t** CovCatFlags = new uint64_t*[nwords];
  if (!CovCatFlags) die("memory allocation error 1 for CovCatFlags in getCovCat_bin_snp_gen()");
  CovCatFlags[0] = new uint64_t[nwords*nCovCathegories];
  if (!CovCatFlags[0]) die("memory allocation error 2 for CovCatFlags[0] in getCovCath_bin_snp_gen()");
  memset(CovCatFlags[0], 0, nwords*nCovCathegories*sizeof(uint64_t));
  for (uint32_t n=1; n<nwords; n++) CovCatFlags[n] = CovCatFlags[n-1]+nCovCathegories;
  uint32_t nr=0, pos=0;
  for (uint32_t k=0; k<nlinestfam; k++) {
    if (qcflag && !person[k].qcin) continue;
    if (pos==64) { pos=0; nr++; }
    CovCatFlags[nr][int(person[k].cov[0])] |= mask[pos];
    pos++;
  }
  while (64>pos) CovCatFlags[nr][0] |= mask[pos++];
  return CovCatFlags;
}


struct PPLLOCATION* getPplLocations(struct PERSON* person, uint32_t nlinestfam, uint32_t npeopleqc, bool qcflag) {
    struct PPLLOCATION* PplLocation = new struct PPLLOCATION[npeopleqc];
    uint32_t nr=0, pos=0;
    for (uint32_t i=0,k=0; k<nlinestfam; k++) {
        if (!person[k].qcin) {
            if (!qcflag) {
                PplLocation[i].nr  = -1;
                PplLocation[i].pos = -1;
                i++;
            }
            continue;
        }
        if (pos==64) { pos=0; nr++; }
        PplLocation[i].nr  = nr;
        PplLocation[i].pos = pos;
        pos++;
        i++;
    }
    return PplLocation;
}


uint8_t getCC(uint64_t* SP, uint8_t pos) {
    register uint64_t m = mask[pos];
    return SP[2] & m ? 2 : SP[1] & m ? 1 : 0;
}


uint8_t getGenotype(uint64_t* SP, uint8_t pos) {
    register uint64_t m = mask[pos];
    return SP[3] & m ? 3 : SP[2] & m ? 2 : SP[1] & m ? 1 : 0;
}


uint64_t* setGenotype(uint64_t* SP, uint8_t pos, uint8_t gt) {
    register uint64_t m = ~mask[pos];
    SP[0]  &=  m;
    SP[1]  &=  m;
    SP[2]  &=  m;
    SP[3]  &=  m;
    SP[gt] |= ~m;
    return SP;
}


uint16_t bitcount64(uint64_t c) {
  static const uint8_t lookup[8192] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,
4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,8,9,9,10,9,10,10,11,9,10,10,11,10,11,11,12,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,
4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,8,9,9,10,9,10,10,11,9,10,10,11,10,11,11,12,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,
4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,8,9,9,10,9,10,10,11,9,10,10,11,10,11,11,12,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,
4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,8,9,9,10,9,10,10,11,9,10,10,11,10,11,11,12,
4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,8,9,9,10,9,10,10,11,9,10,10,11,10,11,11,12,
5,6,6,7,6,7,7,8,6,7,7,8,7,8,8,9,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,8,9,9,10,9,10,10,11,9,10,10,11,10,11,11,12,6,7,7,8,7,8,8,9,7,8,8,9,8,9,9,10,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,8,9,9,10,9,10,10,11,9,10,10,11,10,11,11,12,7,8,8,9,8,9,9,10,8,9,9,10,9,10,10,11,8,9,9,10,9,10,10,11,9,10,10,11,10,11,11,12,8,9,9,10,9,10,10,11,9,10,10,11,10,11,11,12,9,10,10,11,10,11,11,12,10,11,11,12,11,12,12,13};
  return lookup[ c      & 0x1FFF]
       + lookup[(c>>13) & 0x1FFF]
       + lookup[(c>>26) & 0x1FFF]
       + lookup[(c>>39) & 0x1FFF]
       + lookup[(c>>52) & 0x1FFF];
}


bool getbit64(uint64_t c, uint8_t pos) {
    return c & mask[pos];
}


void setbit64(uint64_t& c, uint8_t pos) {
    c |= mask[pos];
}


void delbit64(uint64_t& c, uint8_t pos) {
    c &= ~mask[pos];
}


template <typename T>
string T2bin(T number) {
    uint16_t bits = 8*sizeof(T);
    string result = number & (T)0x1<<(bits-1) ? "1" : "0";
    for (T mask=((T)0x1<<(bits-2)); mask>0; mask>>=1)
        result += (number & mask) ? "1" : "0";
    return result;
}


template <typename T>
void delete_2dim(T**& Array) {
    delete[] Array[0];
    Array[0] = NULL;
    delete[] Array;
    Array = NULL;
}


template <typename T>
void delete_3dim(T***& Array, uint32_t dim1) {
    for (uint32_t i=0; i<dim1; i++) {
        delete[] Array[i][0];
        Array[i][0] = NULL;
        delete[] Array[i];
        Array[i] = NULL;
    }
    delete[] Array;
    Array = NULL;
}


template <typename T>
void delete_4dim(T****& Array, uint32_t dim1, uint32_t dim2) {
    for (uint32_t i=0; i<dim1; i++) {
        for (uint32_t j=0; j<dim2; j++) {
            delete[] Array[i][j][0];
            Array[i][j][0] = NULL;
            delete[] Array[i][j];
            Array[i][j] = NULL;
        }
        delete[] Array[i];
        Array[i] = NULL;
    }
    delete[] Array;
    Array = NULL;
}


template <typename T>
void printmatrix(T** A, uint32_t n, uint32_t m, string name) {
  ofstream file;
  file.open(name.c_str(), ios::trunc);
  for (uint32_t i=0; i<n; i++) {
    for (uint32_t j=0; j<m; j++)
      file << "\t" << A[i][j];
    file << endl;
   }
  file << endl;
  file.close();
}


template <typename T>
void printmatrix(T** A, uint32_t n, string name) {
  printmatrix(A,n,n,name);
}

