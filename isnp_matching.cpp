/**~ isnp_matching.cpp
.---------------------------------------------------------------------------.
|  Software: INTERSNP, Genome-Wide Interaction Analysis                     |
|      Site: http://intersnp.meb.uni-bonn.de/                               |
| ------------------------------------------------------------------------- |
|      File: isnp_matching.cpp                                              |
|    Author: André Lacour                                                   |
|   Content: source-code library on structuring methods for INTERSNP        |
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


template <typename T>
void qsortrelatives(struct RELATIVES*& , int, int, struct RELATIVES&, struct RELATIVES&, T**);
uint32_t clustering_hungarian(int32_t**, uint32_t, uint32_t*&);
struct RELATIVES* maximum_weight_matching(int32_t**, uint32_t, uint32_t);
double pValueCalc(double, double);


struct LOCAL_MATCHING {
	uint8_t chr;
	uint32_t bin_start;
	uint32_t bin_stop;
	uint32_t nGroups;
    struct CLUSTERS* groups;
    uint64_t** BinSNPsCCFlags;
};


double pValNormDist(double d) {
    return erfc(abs(d)*M_SQRT1_2);
}


uint32_t** getWeightMatrix_CC(uint64_t*** BinPPLs, uint32_t* CtrlsMap, uint32_t nctrls, uint32_t* CasesMap, uint32_t ncases, int32_t* PPLMapInv, uint64_t* MatchingInFlags, uint32_t nwords) {
    uint32_t Nmax = max(ncases, nctrls);
    uint32_t** A = new uint32_t*[Nmax];
    if (!A) die("memory allocation error for A in getWeightMatrix_CC2()");
    A[0] = new uint32_t[Nmax*Nmax];
    if (!A[0]) die("memory allocation error for A[0] in getWeightMatrix_CC2()");
    memset(A[0], 0, Nmax*Nmax*sizeof(uint32_t));
    for (uint32_t i=1; i<Nmax; i++) A[i] = A[i-1]+Nmax;
    for (uint32_t i=0; i<nctrls; i++) {
        for (uint32_t j=0; j<ncases; j++) {
            A[i][j] = (uint32_t)( 128.0*nwords*getIbs(BinPPLs, PPLMapInv[CtrlsMap[i]], PPLMapInv[CasesMap[j]], nwords, MatchingInFlags) + 0.5 );
        }
    }
    return A;
}


uint32_t** getWeightMatrix_all(uint64_t*** BinPPLs, uint32_t npeople, uint64_t* MatchingInFlags, uint32_t nwords) {
    uint32_t** B = new uint32_t*[npeople];
    if (!B) die("memory allocation error for B in getWeightMatrix_all()");
    B[0] = new uint32_t[npeople*npeople];
    if (!B[0]) die("memory allocation error for B[0] in getWeightMatrix_all()");
    memset(B[0], 0, npeople*npeople*sizeof(uint32_t));
    for (uint32_t i=1; i<npeople; i++) B[i] = B[i-1]+npeople;
    for (uint32_t i=0; i<npeople; i++) {
        B[i][i] = 0.0;
        for (uint32_t j=i+1; j<npeople; j++) {
            B[i][j] = (uint32_t)( 128*nwords*getIbs(BinPPLs, i, j, nwords, MatchingInFlags) + 0.5 );
            B[j][i] = B[i][j];
        }
    }
    return B;
}


struct CLUSTERS* Affil2Cluster(uint32_t nClusters, struct PERSON* person, uint32_t nlinestfam) {
    struct CLUSTERS* Clusters = new struct CLUSTERS[nClusters];
    memset(Clusters, 0, nClusters*sizeof(struct CLUSTERS));
    #if PARALLEL
    uint32_t thread = omp_get_thread_num();
    #else
    uint32_t thread = 0;
    #endif
    for (uint32_t k=0; k<nlinestfam; k++) {
        if (!person[k].qcin || person[k].clusterAffil<0) continue;
        Clusters[person[k].clusterAffil].nppls++;
        switch (person[k].aff[thread]) {
            case 1: Clusters[person[k].clusterAffil].nctrls++; break;
            case 2: Clusters[person[k].clusterAffil].ncases++; break;
        }
    }
    for (uint32_t c=0; c<nClusters; c++) {
        Clusters[c].list = new uint32_t[Clusters[c].nppls];
        Clusters[c].nppls = 0;
    }
    for (uint32_t k=0; k<nlinestfam; k++) {
        if (!person[k].qcin || person[k].clusterAffil<0) continue;
        Clusters[person[k].clusterAffil].list[Clusters[person[k].clusterAffil].nppls] = k;
        Clusters[person[k].clusterAffil].nppls++;
    }
    return Clusters;
}


bool vicinity_condition(uint32_t ca, uint32_t cb, uint32_t T) {
    return (ca<T && cb<T) || !ca || !cb;
}


double armitageTrendGroupTestLinear(struct CLUSTERS* Groups, uint32_t nGroups, uint64_t** BinSNP, uint64_t** BinSNPsCCFlags, struct PPLLOCATION* PplLocFam, int teststat, double *inflationfactor, double *fulltests) {
    struct PPLLOCATION co;
    uint8_t gt,cc;
	int m[3][4];
	int n[4];
    int o[3];
    int p;
	double tTrend=0;
	double den;
	double var=0;
    for (uint32_t k=0; k<nGroups; k++) {
        memset(m, 0, 3*4*sizeof(int));
        memset(n, 0,   4*sizeof(int));
        memset(o, 0,   3*sizeof(int));
        p = 0;
        for (uint32_t l=0; l<Groups[k].nppls; l++) {
            co = PplLocFam[Groups[k].list[l]];
            gt = getGenotype(BinSNP[co.nr], co.pos);
            cc = getCC(BinSNPsCCFlags[co.nr], co.pos);
            m[cc][gt]++;
        }
        for (gt=1; gt<4; gt++) n[gt] = m[1][gt] + m[2][gt];
        if (!n[1] + !n[2] + !n[3] > 1) continue;
        for (cc=1; cc<3; cc++) o[cc] = m[cc][1] + m[cc][2] + m[cc][3];
        if (!o[1] || !o[2]) continue;
        p = o[1] + o[2];
        den = p*(n[2]+4*n[3]) - pow(n[2]+2*n[3], 2);
        if ( den < EPS ) continue;
        tTrend += ( o[1]*(m[2][2]+2*m[2][3])-o[2]*(m[1][2]+2*m[1][3]) )* sqrt((p-1)/(o[1]*o[2]*den));
        var++;
	}
    #pragma omp critical (INFLATIONFACTOR)
    {
      (*inflationfactor) += tTrend*tTrend/var;
	  (*fulltests)++;
    }
    if (teststat) return tTrend*tTrend/var;
    return pValNormDist(tTrend/sqrt(var));
}


double armitageTrendPairTestLinear(struct RELATIVES* Pairs, uint32_t nPairs, uint64_t** BinSNP, uint64_t** BinSNPsCCFlags, struct PPLLOCATION* PplLocFam, int teststat, double *inflationfactor, double *fulltests) {
    struct PPLLOCATION ci,cj,Ca,Co;
    uint8_t cc,g_ca,g_co;
	double tTrend=0.0;
	int var=0;
    for (uint32_t k=0; k<nPairs; k++) {
        ci = PplLocFam[Pairs[k].i];
        cj = PplLocFam[Pairs[k].j];
        cc = getCC(BinSNPsCCFlags[ci.nr], ci.pos);
        Co = cc==1 ? ci : cj;
        Ca = cc==2 ? ci : cj;
        g_co = getGenotype(BinSNP[Co.nr], Co.pos);
        g_ca = getGenotype(BinSNP[Ca.nr], Ca.pos);
        if (!g_co || !g_ca || g_co==g_ca) continue;
        tTrend += (g_ca>g_co ? -1 : g_ca==g_co ? 0 : 1);
        var++;
	}
    #pragma omp critical (INFLATIONFACTOR)
    {
      (*inflationfactor) += tTrend*tTrend/var;
	  (*fulltests)++;
    }
    if (teststat) return tTrend*tTrend/var;
    return pValNormDist(tTrend/sqrt(var));
}


double armitageTrendGroupTestSquare(struct CLUSTERS* Groups, uint32_t nGroups, uint64_t** BinSNP, uint64_t** BinSNPsCCFlags, struct PPLLOCATION* PplLocFam, int teststat, double *inflationfactor, double *fulltests) {
    struct PPLLOCATION co;
    uint8_t gt,cc;
	int m[3][4];
	int n[4];
    int o[3];
    int p;
	double tTrend=0;
	double den,buf,f_ca,f_co;
    for (uint32_t k=0; k<nGroups; k++) {
        memset(m, 0, 3*4*sizeof(int));
        memset(n, 0,   4*sizeof(int));
        memset(o, 0,   3*sizeof(int));
        p = 0;
        for (uint32_t l=0; l<Groups[k].nppls; l++) {
            co = PplLocFam[Groups[k].list[l]];
            gt = getGenotype(BinSNP[co.nr], co.pos);
            cc = getCC(BinSNPsCCFlags[co.nr], co.pos);
            m[cc][gt]++;
        }
        for (gt=1; gt<4; gt++) n[gt] = m[1][gt] + m[2][gt];
        if (!n[1] + !n[2] + !n[3] > 1) continue;
        for (cc=1; cc<3; cc++) o[cc] = m[cc][1] + m[cc][2] + m[cc][3];
        if (!o[1] || !o[2]) continue;
        p = o[1] + o[2];
        den = p*(n[2]+4*n[3]) - pow((double)(n[2]+2*n[3]), 2);
        if (den < EPS) continue;
        buf = pow(o[1]*(m[2][2]+2*m[2][3]) - o[2]*(m[1][2]+2*m[1][3]), 2) *p/(o[1]*o[2]*den);
        f_co = (2.*m[1][1]+m[1][2])/(2.*Groups[k].nctrls);
        f_ca = (2.*m[2][1]+m[2][2])/(2.*Groups[k].ncases);
        tTrend += (f_ca<f_co ? -buf : buf);
	}
    if (std::isnan(tTrend)) tTrend = 0;
	tTrend = abs(tTrend);
    #pragma omp critical (INFLATIONFACTOR)
    {
      (*inflationfactor) += tTrend;
	  (*fulltests)++;
    }
    if (teststat) return tTrend;
    return pValNormDist(sqrt(tTrend));
}


double armitageTrendPairTestSquare(struct RELATIVES* Pairs, uint32_t nPairs, uint64_t** BinSNP, uint64_t** BinSNPsCCFlags, struct PPLLOCATION* PplLocFam, int teststat, double *inflationfactor, double *fulltests) {
    struct PPLLOCATION ci,cj,Ca,Co;
    uint8_t cc,g_ca,g_co;
	double tTrend=0;
    for (uint32_t k=0; k<nPairs; k++) {
        ci = PplLocFam[Pairs[k].i];
        cj = PplLocFam[Pairs[k].j];
        cc = getCC(BinSNPsCCFlags[ci.nr], ci.pos);
        Co = cc==1 ? ci : cj;
        Ca = cc==2 ? ci : cj;
        g_co = getGenotype(BinSNP[Co.nr], Co.pos);
        g_ca = getGenotype(BinSNP[Ca.nr], Ca.pos);
        if (!g_co || !g_ca || g_co==g_ca) continue;
        tTrend += (g_ca>g_co ? -2 : g_ca==g_co ? 0 : 2);
	}
	tTrend = abs(tTrend);
    #pragma omp critical (INFLATIONFACTOR)
    {
      (*inflationfactor) += tTrend;
	  (*fulltests)++;
    }
    if (teststat) return tTrend;
    return pValNormDist(sqrt(tTrend));
}


uint32_t call_all_all_matching(uint64_t*** BinPPLs, uint32_t nwordsPPLs, struct PERSON* person, uint32_t nlinestfam, uint32_t npplqc, struct MAP* map, uint32_t begin, uint32_t length, struct CLUSTERS*& Groups) {
    uint32_t nGroups = 0;
    uint64_t* MatchingFlags = getMatchingIn_bin_ind_gen(map, begin, length, nwordsPPLs);
    uint32_t** WMAll = getWeightMatrix_all(BinPPLs, npplqc, MatchingFlags, nwordsPPLs);
    struct RELATIVES* Pairs = maximum_weight_matching((int32_t**)WMAll, npplqc, npplqc);
    for (uint32_t i=0; i<nlinestfam; i++) person[i].clusterAffil = -1;
    for (uint32_t n=0; n<npplqc; n++) {
        if (person[Pairs[n].j].clusterAffil > -1) {
            person[Pairs[n].i].clusterAffil = person[Pairs[n].j].clusterAffil;
        } else if (person[Pairs[n].i].clusterAffil > -1) {
            person[Pairs[n].j].clusterAffil = person[Pairs[n].i].clusterAffil;
        } else {
            person[Pairs[n].i].clusterAffil = nGroups;
            person[Pairs[n].j].clusterAffil = nGroups;
            nGroups++;
        }
    }
    Groups = Affil2Cluster(nGroups, person, nlinestfam);
    delete[] MatchingFlags; MatchingFlags=NULL;
    delete_2dim(WMAll);
    delete[] Pairs; Pairs=NULL;
    return nGroups;
}


uint32_t call_matching_ibs_pairs(uint64_t*** BinPPLs, uint32_t nwordsPPLs, struct PERSON* person, uint32_t nlinestfam, uint32_t nctrlsqc, uint32_t ncasesqc, struct MAP* map, uint32_t begin, uint32_t length, struct RELATIVES*& Pairs) {
    uint32_t nMatches = 0;
    #if PARALLEL
    uint32_t thread = omp_get_thread_num();
    #else
    uint32_t thread = 0;
    #endif
    stringstream sstm;
    uint32_t npplqc = nctrlsqc + ncasesqc;
    uint32_t Nmin = min(ncasesqc, nctrlsqc);
    uint32_t Nmax = max(ncasesqc, nctrlsqc);
    uint32_t* CtrlsMap = new uint32_t[nctrlsqc];
    uint32_t* CasesMap = new uint32_t[ncasesqc];
    int32_t* PPLMapInv = new int32_t[nlinestfam];
    memset(PPLMapInv, -1, nlinestfam*sizeof(int32_t));
    for (uint32_t h=0,i=0,j=0,k=0; k<nlinestfam; k++) {
        if (!person[k].qcin) continue;
        switch (person[k].aff[thread]) {
            case 1: CtrlsMap[j++] = k; break;
            case 2: CasesMap[i++] = k; break;
        }
        PPLMapInv[k] = h++;
    }
    uint64_t* MatchingFlags = getMatchingIn_bin_ind_gen(map, begin, length, nwordsPPLs);
    uint32_t** WMAll = getWeightMatrix_all(BinPPLs, npplqc, MatchingFlags, nwordsPPLs);
    delete[] MatchingFlags; MatchingFlags=NULL;
    struct RELATIVES* AllMatches = maximum_weight_matching((int32_t**)WMAll, npplqc, npplqc);
    uint64_t minIBS = ULONG_MAX;
    for (uint32_t i=0; i<npplqc; i++) if (WMAll[AllMatches[i].i][AllMatches[i].j] < minIBS) minIBS = WMAll[AllMatches[i].i][AllMatches[i].j];
    {
    float buf, mean, var;
    buf = 0;
    for (uint32_t i=0; i<npplqc; i++) buf += WMAll[AllMatches[i].i][AllMatches[i].j];
    mean = buf/(double)(128.0*nwordsPPLs*npplqc);
    buf = 0;
    for (uint32_t i=0; i<npplqc; i++) buf += (double)pow(WMAll[AllMatches[i].i][AllMatches[i].j]/(double)(128.0*nwordsPPLs)-mean, 2);
    var = sqrt(buf/(double)(npplqc-1));
    sstm << "meanPairIBS: " << mean << ", sdPairIBS: " << var << ", minPairIBS: " << minIBS/(128.0*nwordsPPLs) << " (all-all matching)";
    logg(sstm);
    }
    uint32_t** WMCC = new uint32_t*[Nmax];
    if (!WMCC) die("memory allocation error for WMCC in call_matching1()");
    WMCC[0] = new uint32_t[Nmax*Nmax];
    if (!WMCC[0]) die("memory allocation error for WMCC[0] in call_matching1()");
    for (uint32_t i=1; i<Nmax; i++) WMCC[i] = WMCC[i-1]+Nmax;
    memset(WMCC[0], 0, Nmax*Nmax*sizeof(uint32_t));
    for (uint32_t i=0; i<nctrlsqc; i++) {
        for (uint32_t j=0; j<ncasesqc; j++) {
            WMCC[i][j] = WMAll[PPLMapInv[CtrlsMap[i]]][PPLMapInv[CasesMap[j]]];
        }
    }
    struct RELATIVES* CCMatches = maximum_weight_matching((int32_t**)WMCC, nctrlsqc, ncasesqc);
    {
    float buf, mean, var;
    buf = 0;
    for (uint32_t i=0; i<Nmin; i++) buf += WMCC[CCMatches[i].i][CCMatches[i].j];
    mean = buf/(double)(128.0*nwordsPPLs*Nmin);
    buf = 0;
    for (uint32_t i=0; i<Nmin; i++) buf += (double)pow(WMCC[CCMatches[i].i][CCMatches[i].j]/(double)(128.0*nwordsPPLs)-mean, 2);
    var = sqrt(buf/(double)(Nmin-1));
    sstm << "meanPairIBS: " << mean << ", sdPairIBS: " << var << " (case-control matching)";
    logg(sstm);
    }
    for (uint32_t i=0; i<Nmin; i++) if (WMCC[CCMatches[i].i][CCMatches[i].j] >= minIBS) nMatches++;
    sstm << "Number of valid Pairs: " << nMatches;
    logg(sstm);
    sstm << "matched: " << 2*nMatches << ", unmatched: " << npplqc-2*nMatches << " (failed IBS-QC: " << 2*(Nmin-nMatches) << ")";
    notice(sstm);
    Pairs = new struct RELATIVES[nMatches];
    for (uint32_t i=0,j=0; i<Nmin; i++) {
        if (WMCC[CCMatches[i].i][CCMatches[i].j] >= minIBS) {
            Pairs[j].i = CtrlsMap[CCMatches[i].i];
            Pairs[j].j = CasesMap[CCMatches[i].j];
            j++;
        }
    }
    delete[] CtrlsMap; CtrlsMap=NULL;
    delete[] CasesMap; CasesMap=NULL;
    delete[] PPLMapInv; PPLMapInv=NULL;
    delete_2dim(WMAll);
    delete[] AllMatches; AllMatches=NULL;
    delete_2dim(WMCC);
    delete[] CCMatches; CCMatches=NULL;
    return nMatches;
}


uint32_t call_matching_ibs_groups1(uint64_t*** BinPPLs, uint32_t nwordsPPLs, struct PERSON* person, uint32_t nlinestfam, uint32_t nctrlsqc, uint32_t ncasesqc, struct MAP* map, uint32_t begin, uint32_t length, struct CLUSTERS*& Groups) {
    uint32_t nGroups = 0;
    #if PARALLEL
    uint32_t thread = omp_get_thread_num();
    #else
    uint32_t thread = 0;
    #endif
    stringstream sstm;
    uint32_t npplqc = nctrlsqc + ncasesqc;
    uint32_t Nmin = min(ncasesqc, nctrlsqc);
    uint32_t Nmax = max(ncasesqc, nctrlsqc);
    uint32_t* CtrlsMap   = new uint32_t[nctrlsqc];
    uint32_t* CasesMap   = new uint32_t[ncasesqc];
     int32_t* PPLMapInv  = new  int32_t[nlinestfam];
    memset(PPLMapInv, -1, nlinestfam*sizeof(int32_t));
    for (uint32_t h=0,i=0,j=0,k=0; k<nlinestfam; k++) {
        if (!person[k].qcin) continue;
        switch (person[k].aff[thread]) {
            case 1: CtrlsMap[j++] = k; break;
            case 2: CasesMap[i++] = k; break;
        }
        PPLMapInv[k] = h++;
    }
    uint64_t* MatchingFlags = getMatchingIn_bin_ind_gen(map, begin, length, nwordsPPLs);
    uint32_t** WMAll = getWeightMatrix_all(BinPPLs, npplqc, MatchingFlags, nwordsPPLs);
    delete[] MatchingFlags; MatchingFlags=NULL;
    struct RELATIVES* AllMatches = maximum_weight_matching((int32_t**)WMAll, npplqc, npplqc);
    uint64_t minIBS = ULONG_MAX;
    for (uint32_t i=0; i<npplqc; i++) if (WMAll[AllMatches[i].i][AllMatches[i].j] < minIBS) minIBS = WMAll[AllMatches[i].i][AllMatches[i].j];
    {
    float buf, mean, var;
    buf = 0;
    for (uint32_t i=0; i<npplqc; i++) buf += WMAll[AllMatches[i].i][AllMatches[i].j];
    mean = buf/(double)(128.0*nwordsPPLs*npplqc);
    buf = 0;
    for (uint32_t i=0; i<npplqc; i++) buf += (double)pow(WMAll[AllMatches[i].i][AllMatches[i].j]/(double)(128.0*nwordsPPLs)-mean, 2);
    var = sqrt(buf/(double)(npplqc-1));
    sstm << "meanPairIBS: " << mean << ", sdPairIBS: " << var << ", minPairIBS: " << minIBS/(128.0*nwordsPPLs) << " (all-all matching)";
    logg(sstm);
    }
    uint32_t** WMCC = new uint32_t*[Nmax];
    if (!WMCC) die("memory allocation error for WMCC in call_matching2()");
    WMCC[0] = new uint32_t[Nmax*Nmax];
    if (!WMCC[0]) die("memory allocation error for WMCC[0] in call_matching2()");
    for (uint32_t i=1; i<Nmax; i++) WMCC[i] = WMCC[i-1]+Nmax;
    memset(WMCC[0], 0, Nmax*Nmax*sizeof(uint32_t));
    for (uint32_t i=0; i<nctrlsqc; i++) {
        for (uint32_t j=0; j<ncasesqc; j++) {
            WMCC[i][j] = WMAll[PPLMapInv[CtrlsMap[i]]][PPLMapInv[CasesMap[j]]];
        }
    }
    struct RELATIVES* CCMatches = maximum_weight_matching((int32_t**)WMCC, nctrlsqc, ncasesqc);
    {
    float buf, mean, var;
    buf = 0;
    for (uint32_t i=0; i<Nmin; i++) buf += WMCC[CCMatches[i].i][CCMatches[i].j];
    mean = buf/(double)(128.0*nwordsPPLs*Nmin);
    buf = 0;
    for (uint32_t i=0; i<Nmin; i++) buf += pow(WMCC[CCMatches[i].i][CCMatches[i].j]/(double)(128.0*nwordsPPLs)-mean, 2);
    var = sqrt(buf/(double)(Nmin-1));
    sstm << "meanPairIBS: " << mean << ", sdPairIBS: " << var << " (case-control matching)";
    logg(sstm);
    }
    int32_t* GroupAffil = new int32_t[nlinestfam];
    memset(GroupAffil, 0, nlinestfam*sizeof(int32_t));
    uint32_t failed = 0;
    uint32_t unmatched = npplqc;
    for (uint32_t n=0; n<Nmin; n++) {
        if (WMCC[CCMatches[n].i][CCMatches[n].j] >= minIBS) {
            unmatched -= 2;
            nGroups++;
            GroupAffil[CtrlsMap[CCMatches[n].i]] = nGroups;
            GroupAffil[CasesMap[CCMatches[n].j]] = nGroups;
        } else if (WMCC[CCMatches[n].i][CCMatches[n].j] > 0) failed+=2;
    }
    sstm << "Number of Groups: " << nGroups;
    logg(sstm);
    sstm << "matched: " << npplqc-unmatched << ", unmatched: " << unmatched << " (failed IBS-QC: " << failed << ")";
    notice(sstm);
    delete[] CCMatches; CCMatches=NULL;
    uint32_t oldunmatched = npplqc;
    while (oldunmatched > unmatched && unmatched > 0) {
        oldunmatched = unmatched;

        for (uint32_t i=0; i<nctrlsqc; i++) {
            for (uint32_t j=0; j<ncasesqc; j++) {
                if ((bool)GroupAffil[CtrlsMap[i]] == (bool)GroupAffil[CasesMap[j]]) WMCC[i][j] = 0;
            }
        }
        CCMatches = maximum_weight_matching((int32_t**)WMCC, nctrlsqc, ncasesqc);

        failed = 0;
        for (uint32_t n=0; n<Nmin; n++) {
            if (WMCC[CCMatches[n].i][CCMatches[n].j] >= minIBS) {
                unmatched--;
                if (!GroupAffil[CtrlsMap[CCMatches[n].i]]) GroupAffil[CtrlsMap[CCMatches[n].i]] = GroupAffil[CasesMap[CCMatches[n].j]];
                if (!GroupAffil[CasesMap[CCMatches[n].j]]) GroupAffil[CasesMap[CCMatches[n].j]] = GroupAffil[CtrlsMap[CCMatches[n].i]];
            } else if (WMCC[CCMatches[n].i][CCMatches[n].j] > 0) failed++;
        }
        sstm << "matched: " << npplqc-unmatched << ", unmatched: " << unmatched << " (failed IBS-QC: " << failed << ")";
        notice(sstm);
        delete[] CCMatches; CCMatches=NULL;
    }
    for (uint32_t i=0; i<nlinestfam; i++) person[i].clusterAffil = GroupAffil[i]-1;
    delete[] GroupAffil; GroupAffil=NULL;
    Groups = Affil2Cluster(nGroups, person, nlinestfam);
    delete[] CtrlsMap; CtrlsMap=NULL;
    delete[] CasesMap; CasesMap=NULL;
    delete[] PPLMapInv; PPLMapInv=NULL;
    delete_2dim(WMAll);
    delete[] AllMatches; AllMatches=NULL;
    delete_2dim(WMCC);
    return nGroups;
}


uint32_t call_matching_ibs_groups2(uint64_t*** BinPPLs, uint32_t nwordsPPLs, struct PERSON* person, uint32_t nlinestfam, uint32_t nctrlsqc, uint32_t ncasesqc, struct MAP* map, uint32_t begin, uint32_t length, struct CLUSTERS*& Groups) {
    uint32_t nGroups = 0;
    #if PARALLEL
    uint32_t thread = omp_get_thread_num();
    #else
    uint32_t thread = 0;
    #endif
    stringstream sstm;
    uint32_t npplqc = nctrlsqc + ncasesqc;
    uint32_t Nmin = min(ncasesqc, nctrlsqc);
    uint32_t Nmax = max(ncasesqc, nctrlsqc);
    uint32_t fac = (uint32_t)ceil(((double)Nmax-EPS)/Nmin);
    uint32_t dim = fac*Nmin;
    uint32_t* CtrlsMap   = new uint32_t[nctrlsqc];
    uint32_t* CasesMap   = new uint32_t[ncasesqc];
     int32_t* PPLMapInv  = new  int32_t[nlinestfam];
    memset(PPLMapInv, -1, nlinestfam*sizeof(int32_t));
    for (uint32_t h=0,i=0,j=0,k=0; k<nlinestfam; k++) {
        if (!person[k].qcin) continue;
        switch (person[k].aff[thread]) {
            case 1: CtrlsMap[j++] = k; break;
            case 2: CasesMap[i++] = k; break;
        }
        PPLMapInv[k] = h++;
    }
    uint64_t* MatchingFlags = getMatchingIn_bin_ind_gen(map, begin, length, nwordsPPLs);
    uint32_t** WMAll = getWeightMatrix_all(BinPPLs, npplqc, MatchingFlags, nwordsPPLs);
    delete[] MatchingFlags; MatchingFlags=NULL;
    struct RELATIVES* AllMatches = maximum_weight_matching((int32_t**)WMAll, npplqc, npplqc);
    uint64_t minIBS = ULONG_MAX;
    for (uint32_t i=0; i<npplqc; i++) if (WMAll[AllMatches[i].i][AllMatches[i].j] < minIBS) minIBS = WMAll[AllMatches[i].i][AllMatches[i].j];
    {
    float buf, mean, var;
    buf = 0;
    for (uint32_t i=0; i<npplqc; i++) buf += WMAll[AllMatches[i].i][AllMatches[i].j];
    mean = buf/(double)(128.0*nwordsPPLs*npplqc);
    buf = 0;
    for (uint32_t i=0; i<npplqc; i++) buf += (double)pow(WMAll[AllMatches[i].i][AllMatches[i].j]/(double)(128.0*nwordsPPLs)-mean, 2);
    var = sqrt(buf/(double)(npplqc-1));
    sstm << "meanPairIBS: " << mean << ", sdPairIBS: " << var << ", minPairIBS: " << minIBS/(128.0*nwordsPPLs) << " (all-all matching)";
    logg(sstm);
    }
    uint32_t** WMCC = new uint32_t*[dim];
    if (!WMCC) die("memory allocation error for WMCC in call_matching3()");
    WMCC[0] = new uint32_t[dim*dim];
    if (!WMCC[0]) die("memory allocation error for WMCC[0] in call_matching3()");
    for (uint32_t i=1; i<dim; i++) WMCC[i] = WMCC[i-1]+dim;
    memset(WMCC[0], 0, dim*dim*sizeof(uint32_t));
    uint32_t ni = (nctrlsqc==Nmax ? Nmax : dim);
    uint32_t nj = (ncasesqc==Nmax ? Nmax : dim);
    for (uint32_t i=0; i<ni; i++) {
        for (uint32_t j=0; j<nj; j++) {
            WMCC[i][j] = WMAll[PPLMapInv[CtrlsMap[i%nctrlsqc]]][PPLMapInv[CasesMap[j%ncasesqc]]];
        }
    }
    for (uint32_t i=nctrlsqc; i<dim; i++) {
        for (uint32_t j=ncasesqc; j<dim; j++) {
            WMCC[i][j] = 128*nwordsPPLs + 1;
        }
    }
    struct RELATIVES* CCMatches = maximum_weight_matching((int32_t**)WMCC, ni, nj);
    for (uint32_t n=0; n<Nmax; n++) {
        CCMatches[n].i = CCMatches[n].i%nctrlsqc;
        CCMatches[n].j = CCMatches[n].j%ncasesqc;
    }
    {
    float buf, mean, var;
    buf = 0;
    for (uint32_t i=0; i<Nmax; i++) buf += WMCC[CCMatches[i].i][CCMatches[i].j];
    mean = buf/(double)(128.0*nwordsPPLs*Nmax);
    buf = 0;
    for (uint32_t i=0; i<Nmax; i++) buf += pow(WMCC[CCMatches[i].i][CCMatches[i].j]/(double)(128.0*nwordsPPLs)-mean, 2);
    var = sqrt(buf/(double)(Nmax-1));
    sstm << "meanPairIBS: " << mean << ", sdPairIBS: " << var << " (case-control matching)";
    logg(sstm);
    }
    int32_t* GroupAffil = new int32_t[nlinestfam];
    memset(GroupAffil, 0, nlinestfam*sizeof(int32_t));
    uint32_t unmatched = npplqc;
    uint32_t nPairs = 0;
    uint32_t failed = 0;
    for (uint32_t n=0; n<Nmax; n++) {
        if (WMCC[CCMatches[n].i][CCMatches[n].j] >= minIBS) {
            unmatched--;
            nPairs++;
            if (!GroupAffil[CtrlsMap[CCMatches[n].i]] && !GroupAffil[CasesMap[CCMatches[n].j]]) {
                unmatched--;
                nGroups++;
                GroupAffil[CtrlsMap[CCMatches[n].i]] = nGroups;
                GroupAffil[CasesMap[CCMatches[n].j]] = nGroups;
            }
            else if (GroupAffil[CasesMap[CCMatches[n].j]]) GroupAffil[CtrlsMap[CCMatches[n].i]] = GroupAffil[CasesMap[CCMatches[n].j]];
            else if (GroupAffil[CtrlsMap[CCMatches[n].i]]) GroupAffil[CasesMap[CCMatches[n].j]] = GroupAffil[CtrlsMap[CCMatches[n].i]];
        } else {
            if (WMCC[CCMatches[n].i][CCMatches[n].j] > 0) {
                if (!GroupAffil[CtrlsMap[CCMatches[n].i]]) failed++;
                if (!GroupAffil[CasesMap[CCMatches[n].j]]) failed++;
            }
        }
    }
    sstm << "Number of Groups: " << nGroups << " (Number of valid Pairs: " << nPairs << ")";
    logg(sstm);
    sstm << "matched: " << npplqc-unmatched << ", unmatched: " << unmatched << " (failed IBS-QC: " << failed << ")";
    notice(sstm);
    delete[] CCMatches; CCMatches=NULL;
    uint32_t oldunmatched = npplqc;
    while (oldunmatched > unmatched && unmatched > 0) {
        oldunmatched = unmatched;

        for (uint32_t i=0; i<nctrlsqc; i++) {
            for (uint32_t j=0; j<ncasesqc; j++) {
                if ((bool)GroupAffil[CtrlsMap[i]] == (bool)GroupAffil[CasesMap[j]]) WMCC[i][j] = 0;
            }
        }
        CCMatches = maximum_weight_matching((int32_t**)WMCC, nctrlsqc, ncasesqc);

        failed = 0;
        for (uint32_t n=0; n<Nmin; n++) {
            if (WMCC[CCMatches[n].i][CCMatches[n].j] >= minIBS) {
                unmatched--;
                if (!GroupAffil[CtrlsMap[CCMatches[n].i]]) GroupAffil[CtrlsMap[CCMatches[n].i]] = GroupAffil[CasesMap[CCMatches[n].j]];
                if (!GroupAffil[CasesMap[CCMatches[n].j]]) GroupAffil[CasesMap[CCMatches[n].j]] = GroupAffil[CtrlsMap[CCMatches[n].i]];
            } else if (WMCC[CCMatches[n].i][CCMatches[n].j] > 0) failed++;
        }
        sstm << "matched: " << npplqc-unmatched << ", unmatched: " << unmatched << " (failed IBS-QC: " << failed << ")";
        notice(sstm);
        delete[] CCMatches; CCMatches=NULL;
    }
    for (uint32_t i=0; i<nlinestfam; i++) person[i].clusterAffil = GroupAffil[i]-1;
    delete[] GroupAffil; GroupAffil=NULL;
    Groups = Affil2Cluster(nGroups, person, nlinestfam);
    delete[] CtrlsMap; CtrlsMap=NULL;
    delete[] CasesMap; CasesMap=NULL;
    delete[] PPLMapInv; PPLMapInv=NULL;
    delete_2dim(WMAll);
    delete[] AllMatches; AllMatches=NULL;
    delete_2dim(WMCC);
    return nGroups;
}


uint32_t call_clustering(uint64_t*** BinPPLs, uint32_t nwordsPPLs, struct PERSON* person, uint32_t nlinestfam, uint32_t nctrlsqc, uint32_t ncasesqc, struct MAP* map, uint32_t begin, uint32_t length, struct CLUSTERS*& Clusters) {
    uint32_t nClusters = 0;
    stringstream sstm;
    uint32_t npplqc = nctrlsqc + ncasesqc;

    uint64_t* MatchingFlags = getMatchingIn_bin_ind_gen(map, begin, length, nwordsPPLs);
    uint32_t** WMAll = getWeightMatrix_all(BinPPLs, npplqc, MatchingFlags, nwordsPPLs);
    delete[] MatchingFlags; MatchingFlags=NULL;
    uint32_t* Affil;
    nClusters = clustering_hungarian((int32_t**)WMAll, npplqc, Affil);
    delete_2dim(WMAll);
    sstm << "Number of Clusters: " << nClusters;
    logg(sstm);
    for (uint32_t j=0,i=0; i<nlinestfam; i++) {
        person[i].clusterAffil = person[i].qcin ? Affil[j++] : -1;
    }
    delete[] Affil; Affil=NULL;
    Clusters = Affil2Cluster(nClusters, person, nlinestfam);
    return nClusters;
}


uint32_t call_matching_cluster_pairs(uint64_t*** BinPPLs, uint32_t nwordsPPLs, struct PERSON* person, uint32_t nlinestfam, uint32_t nctrlsqc, uint32_t ncasesqc, struct MAP* map, uint32_t begin, uint32_t length, struct RELATIVES*& Pairs) {
    uint32_t nMatches = 0;
    #if PARALLEL
    uint32_t thread = omp_get_thread_num();
    #else
    uint32_t thread = 0;
    #endif
    stringstream sstm;
    uint32_t npplqc = nctrlsqc + ncasesqc;
    uint32_t Nmin = min(ncasesqc, nctrlsqc);
    uint32_t Nmax = max(ncasesqc, nctrlsqc);
    uint32_t* CtrlsMap = new uint32_t[nctrlsqc];
    uint32_t* CasesMap = new uint32_t[ncasesqc];
    int32_t* PPLMapInv = new int32_t[nlinestfam];
    memset(PPLMapInv, -1, nlinestfam*sizeof(int32_t));
    for (uint32_t h=0,i=0,j=0,k=0; k<nlinestfam; k++) {
        if (!person[k].qcin) continue;
        switch (person[k].aff[thread]) {
            case 1: CtrlsMap[j++] = k; break;
            case 2: CasesMap[i++] = k; break;
        }
        PPLMapInv[k] = h++;
    }
    uint64_t* MatchingFlags = getMatchingIn_bin_ind_gen(map, begin, length, nwordsPPLs);
    uint32_t** WMAll = getWeightMatrix_all(BinPPLs, npplqc, MatchingFlags, nwordsPPLs);
    delete[] MatchingFlags; MatchingFlags=NULL;
    uint32_t* Affil;
    uint32_t nClusters = clustering_hungarian((int32_t**)WMAll, npplqc, Affil);
    sstm << "Number of Clusters: " << nClusters;
    logg(sstm);
    uint32_t** WMCC = new uint32_t*[Nmax];
    if (!WMCC) die("memory allocation error for WMCC in call_cluster_matching1()");
    WMCC[0] = new uint32_t[Nmax*Nmax];
    if (!WMCC[0]) die("memory allocation error for WMCC[0] in call_cluster_matching1()");
    for (uint32_t i=1; i<Nmax; i++) WMCC[i] = WMCC[i-1]+Nmax;
    memset(WMCC[0], 0, Nmax*Nmax*sizeof(uint32_t));
    for (uint32_t i=0; i<nctrlsqc; i++) {
        for (uint32_t j=0; j<ncasesqc; j++) {
            if (Affil[PPLMapInv[CtrlsMap[i]]] == Affil[PPLMapInv[CasesMap[j]]]) WMCC[i][j] = WMAll[PPLMapInv[CtrlsMap[i]]][PPLMapInv[CasesMap[j]]];
        }
    }
    struct RELATIVES* CCMatches = maximum_weight_matching((int32_t**)WMCC, nctrlsqc, ncasesqc);
    for (uint32_t i=0; i<Nmin; i++) if (WMCC[CCMatches[i].i][CCMatches[i].j] > 0) nMatches++;
    sstm << "Number of valid Pairs: " << nMatches;
    logg(sstm);
    Pairs = new struct RELATIVES[nMatches];
    for (uint32_t i=0,j=0; i<Nmin; i++) {
        if (WMCC[CCMatches[i].i][CCMatches[i].j] > 0) {
            Pairs[j].i = CtrlsMap[CCMatches[i].i];
            Pairs[j].j = CasesMap[CCMatches[i].j];
            j++;
        }
    }
    delete[] Affil;
    delete[] CtrlsMap; CtrlsMap=NULL;
    delete[] CasesMap; CasesMap=NULL;
    delete[] PPLMapInv; PPLMapInv=NULL;
    delete[] MatchingFlags; MatchingFlags=NULL;
    delete[] CCMatches; CCMatches=NULL;
    delete_2dim(WMCC);
    delete_2dim(WMAll);
    return nMatches;
}


uint32_t call_matching_cluster_groups1(uint64_t*** BinPPLs, uint32_t nwordsPPLs, struct PERSON* person, uint32_t nlinestfam, uint32_t nctrlsqc, uint32_t ncasesqc, struct MAP* map, uint32_t begin, uint32_t length, struct CLUSTERS*& Groups) {
    uint32_t nGroups = 0;
    #if PARALLEL
    uint32_t thread = omp_get_thread_num();
    #else
    uint32_t thread = 0;
    #endif
    stringstream sstm;
    uint32_t npplqc = nctrlsqc + ncasesqc;
    uint32_t Nmin = min(ncasesqc, nctrlsqc);
    uint32_t Nmax = max(ncasesqc, nctrlsqc);
    uint32_t* CtrlsMap   = new uint32_t[nctrlsqc];
    uint32_t* CasesMap   = new uint32_t[ncasesqc];
     int32_t* PPLMapInv  = new  int32_t[nlinestfam];
    memset(PPLMapInv, -1, nlinestfam*sizeof(int32_t));
    for (uint32_t h=0,i=0,j=0,k=0; k<nlinestfam; k++) {
        if (!person[k].qcin) continue;
        switch (person[k].aff[thread]) {
            case 1: CtrlsMap[j++] = k; break;
            case 2: CasesMap[i++] = k; break;
        }
        PPLMapInv[k] = h++;
    }
    uint64_t* MatchingFlags = getMatchingIn_bin_ind_gen(map, begin, length, nwordsPPLs);
    uint32_t** WMAll = getWeightMatrix_all(BinPPLs, npplqc, MatchingFlags, nwordsPPLs);
    delete[] MatchingFlags; MatchingFlags=NULL;
    uint32_t* ClusterAffil;
    uint32_t nClusters = clustering_hungarian((int32_t**)WMAll, npplqc, ClusterAffil);
    sstm << "Number of Clusters: " << nClusters;
    logg(sstm);
    uint32_t** WMCC = new uint32_t*[Nmax];
    if (!WMCC) die("memory allocation error for WMCC in call_cluster_matching2()");
    WMCC[0] = new uint32_t[Nmax*Nmax];
    if (!WMCC[0]) die("memory allocation error for WMCC[0] in call_cluster_matching2()");
    for (uint32_t i=1; i<Nmax; i++) WMCC[i] = WMCC[i-1]+Nmax;
    memset(WMCC[0], 0, Nmax*Nmax*sizeof(uint32_t));
    for (uint32_t i=0; i<nctrlsqc; i++) {
        for (uint32_t j=0; j<ncasesqc; j++) {
            if (ClusterAffil[PPLMapInv[CtrlsMap[i]]] == ClusterAffil[PPLMapInv[CasesMap[j]]]) WMCC[i][j] = WMAll[PPLMapInv[CtrlsMap[i]]][PPLMapInv[CasesMap[j]]];
        }
    }
    struct RELATIVES* CCMatches = maximum_weight_matching((int32_t**)WMCC, nctrlsqc, ncasesqc);
    int32_t* GroupAffil = new int32_t[nlinestfam];
    memset(GroupAffil, 0, nlinestfam*sizeof(int32_t));
    uint32_t unmatched = npplqc;
    for (uint32_t n=0; n<Nmin; n++) {
        if (WMCC[CCMatches[n].i][CCMatches[n].j] > 0) {
            unmatched -= 2;
            nGroups++;
            GroupAffil[CtrlsMap[CCMatches[n].i]] = nGroups;
            GroupAffil[CasesMap[CCMatches[n].j]] = nGroups;
        }
    }
    sstm << "Number of Groups: " << nGroups;
    logg(sstm);
    sstm << "matched: " << npplqc-unmatched << ", unmatched: " << unmatched;
    notice(sstm);
    delete[] CCMatches; CCMatches=NULL;
    uint32_t oldunmatched = npplqc;
    while (oldunmatched > unmatched && unmatched > 0) {
        oldunmatched = unmatched;
        for (uint32_t i=0; i<nctrlsqc; i++) {
            for (uint32_t j=0; j<ncasesqc; j++) {
                if ((bool)GroupAffil[CtrlsMap[i]] == (bool)GroupAffil[CasesMap[j]]) WMCC[i][j] = 0;
            }
        }
        CCMatches = maximum_weight_matching((int32_t**)WMCC, nctrlsqc, ncasesqc);
        for (uint32_t n=0; n<Nmin; n++) {
            if (WMCC[CCMatches[n].i][CCMatches[n].j] > 0) {
                if (!GroupAffil[CtrlsMap[CCMatches[n].i]]) GroupAffil[CtrlsMap[CCMatches[n].i]] = GroupAffil[CasesMap[CCMatches[n].j]];
                if (!GroupAffil[CasesMap[CCMatches[n].j]]) GroupAffil[CasesMap[CCMatches[n].j]] = GroupAffil[CtrlsMap[CCMatches[n].i]];
                unmatched--;
            }
        }
        sstm << "matched: " << npplqc-unmatched << ", unmatched: " << unmatched;
        notice(sstm);
        delete[] CCMatches; CCMatches=NULL;
    }
    for (uint32_t i=0; i<nlinestfam; i++) person[i].clusterAffil = GroupAffil[i]-1;
    delete[] GroupAffil; GroupAffil=NULL;
    Groups = Affil2Cluster(nGroups, person, nlinestfam);
    delete[] CtrlsMap; CtrlsMap=NULL;
    delete[] CasesMap; CasesMap=NULL;
    delete[] PPLMapInv; PPLMapInv=NULL;
    delete_2dim(WMAll);
    delete_2dim(WMCC);
    delete[] ClusterAffil; ClusterAffil=NULL;
    return nGroups;
}


uint32_t call_matching_cluster_groups2(uint64_t*** BinPPLs, uint32_t nwordsPPLs, struct PERSON* person, uint32_t nlinestfam, uint32_t nctrlsqc, uint32_t ncasesqc, struct MAP* map, uint32_t begin, uint32_t length, struct CLUSTERS*& Groups) {
    uint32_t nGroups = 0;
    #if PARALLEL
    uint32_t thread = omp_get_thread_num();
    #else
    uint32_t thread = 0;
    #endif
    stringstream sstm;
    uint32_t npplqc = nctrlsqc + ncasesqc;
    int32_t* PPLMapInv = new int32_t[nlinestfam];
    memset(PPLMapInv, -1, nlinestfam*sizeof(int32_t));
    for (uint32_t h=0,k=0; k<nlinestfam; k++) {
        if (person[k].qcin) PPLMapInv[k] = h++;
    }
    uint64_t* MatchingFlags = getMatchingIn_bin_ind_gen(map, begin, length, nwordsPPLs);
    uint32_t** WMAll = getWeightMatrix_all(BinPPLs, npplqc, MatchingFlags, nwordsPPLs);
    delete[] MatchingFlags; MatchingFlags=NULL;
    uint32_t* ClusterAffil;
    uint32_t nClusters = clustering_hungarian((int32_t**)WMAll, npplqc, ClusterAffil);
    sstm << "Number of Clusters: " << nClusters;
    logg(sstm);
    for (uint32_t j=0,i=0; i<nlinestfam; i++) {
        person[i].clusterAffil = person[i].qcin ? ClusterAffil[j++] : -1;
    }
    delete[] ClusterAffil; ClusterAffil=NULL;
    struct CLUSTERS* Clusters = Affil2Cluster(nClusters, person, nlinestfam);
    int32_t* GroupAffil = new int32_t[nlinestfam];
    memset(GroupAffil, 0, nlinestfam*sizeof(int32_t));
    uint32_t nPairs = 0;
    uint32_t nmatched = 0;
    for (uint32_t c=0; c<nClusters; c++) {
        if (!Clusters[c].ncases || !Clusters[c].nctrls) continue;
        uint32_t Nmin = min(Clusters[c].ncases, Clusters[c].nctrls);
        uint32_t Nmax = max(Clusters[c].ncases, Clusters[c].nctrls);
        uint32_t fac = (uint32_t)ceil(((double)Nmax-0.5)/Nmin);
        uint32_t dim = fac*Nmin;
        uint32_t ni = (Clusters[c].nctrls==Nmax ? Nmax : dim);
        uint32_t nj = (Clusters[c].ncases==Nmax ? Nmax : dim);
        uint32_t** WMCC = new uint32_t*[dim];
        if (!WMCC) die("memory allocation error for WMCC in call_cluster_matching4()");
        WMCC[0] = new uint32_t[dim*dim];
        if (!WMCC[0]) die("memory allocation error for WMCC[0] in call_cluster_matching4()");
        for (uint32_t i=1; i<dim; i++) WMCC[i] = WMCC[i-1]+dim;
        memset(WMCC[0], 0, dim*dim*sizeof(uint32_t));
        for (uint32_t i=0,k=0; k<Clusters[c].nppls; k++) {
            if (person[Clusters[c].list[k]].aff[thread]!=1) continue;
            for (uint32_t j=0,l=0; l<Clusters[c].nppls; l++) {
                if (person[Clusters[c].list[l]].aff[thread]!=2) continue;
                WMCC[i][j] = WMAll[PPLMapInv[Clusters[c].list[k]]][PPLMapInv[Clusters[c].list[l]]];
                if (Clusters[c].nctrls!=Nmax) for (uint32_t f=1; f<fac; f++) WMCC[i+f*Clusters[c].nctrls][j] = WMCC[i][j];
                if (Clusters[c].ncases!=Nmax) for (uint32_t f=1; f<fac; f++) WMCC[i][j+f*Clusters[c].ncases] = WMCC[i][j];
                j++;
            }
            i++;
        }
        for (uint32_t i=Clusters[c].nctrls; i<dim; i++) {
            for (uint32_t j=Clusters[c].ncases; j<dim; j++) {
                WMCC[i][j] = 128*nwordsPPLs + 1;
            }
        }

        struct RELATIVES* CCMatches = maximum_weight_matching((int32_t**)WMCC, ni, nj);
        for (uint32_t n=0; n<Nmax; n++) {
            CCMatches[n].i = CCMatches[n].i%Clusters[c].nctrls;
            CCMatches[n].j = CCMatches[n].j%Clusters[c].ncases;
        }

        for (uint32_t n=0; n<Nmax; n++) {
            uint32_t i=-1,j=-1;
            if (!WMCC[CCMatches[n].i][CCMatches[n].j]) continue;
            for (int32_t a=0,b=0,k=0; k<(int32_t)Clusters[c].nppls; k++) {
                switch (person[Clusters[c].list[k]].aff[thread]) {
                    case 1: if (CCMatches[n].i==a++) i=Clusters[c].list[k]; break;
                    case 2: if (CCMatches[n].j==b++) j=Clusters[c].list[k]; break;
                }
            }
            nmatched++;
            nPairs++;
            if (!GroupAffil[i] && !GroupAffil[j]) {
                nmatched++;
                nGroups++;
                GroupAffil[i] = nGroups;
                GroupAffil[j] = nGroups;
            }
            else if (GroupAffil[j]) GroupAffil[i] = GroupAffil[j];
            else if (GroupAffil[i]) GroupAffil[j] = GroupAffil[i];
        }

        delete[] CCMatches; CCMatches=NULL;
        delete_2dim(WMCC);
    }
    sstm << "Number of Groups: " << nGroups << " (Number of valid Pairs: " << nPairs << ")";
    logg(sstm);
    sstm << "matched: " << nmatched << ", unmatched: " << npplqc-nmatched;
    notice(sstm);
    for (uint32_t i=0; i<nlinestfam; i++) person[i].clusterAffil = GroupAffil[i]-1;
    delete[] GroupAffil; GroupAffil=NULL;
    Groups = Affil2Cluster(nGroups, person, nlinestfam);
    delete[] PPLMapInv; PPLMapInv=NULL;
    delete_2dim(WMAll);
    for (uint32_t c=0; c<nClusters; c++) { delete[] Clusters[c].list; Clusters[c].list=NULL; }
    delete[] Clusters; Clusters=NULL;
    delete[] GroupAffil; GroupAffil=NULL;
    return nGroups;
}


uint32_t call_matching_vicinity_pairs(uint64_t*** BinPPLs, uint32_t nwordsPPLs, struct PERSON* person, uint32_t nlinestfam, uint32_t nctrlsqc, uint32_t ncasesqc, struct MAP* map, uint32_t begin, uint32_t length, struct RELATIVES*& Pairs) {
    extern string outputname;
    uint32_t nMatches = 0;
    #if PARALLEL
    uint32_t thread = omp_get_thread_num();
    #else
    uint32_t thread = 0;
    #endif
    stringstream sstm;
    uint32_t npplqc = ncasesqc+nctrlsqc;
    uint32_t Nmin = min(ncasesqc, nctrlsqc);
    uint32_t Nmax = max(ncasesqc, nctrlsqc);
    uint8_t T = (uint8_t)ceil(log(nctrlsqc*ncasesqc));
    uint32_t* CtrlsMap = new uint32_t[nctrlsqc];
    uint32_t* CasesMap = new uint32_t[ncasesqc];
    int32_t* PPLMapInv = new int32_t[nlinestfam];
    memset(PPLMapInv, -1, nlinestfam*sizeof(int32_t));
    for (uint32_t h=0,i=0,j=0,k=0; k<nlinestfam; k++) {
        if (!person[k].qcin) continue;
        switch (person[k].aff[thread]) {
            case 1: CtrlsMap[j++] = k; break;
            case 2: CasesMap[i++] = k; break;
        }
        PPLMapInv[k] = h++;
    }
    uint64_t* MatchingFlags = getMatchingIn_bin_ind_gen(map, begin, length, nwordsPPLs);
    uint32_t** WMCC = getWeightMatrix_CC(BinPPLs, CtrlsMap, nctrlsqc, CasesMap, ncasesqc, PPLMapInv, MatchingFlags, nwordsPPLs);
    delete[] MatchingFlags; MatchingFlags=NULL;
    struct RELATIVES* CCMatches = maximum_weight_matching((int32_t**)WMCC, nctrlsqc, ncasesqc);
    bool* pair_in = new bool[Nmin];
    memset(pair_in, 0, Nmin*sizeof(bool));
    for (uint32_t i=0; i<Nmin; i++) {
        uint32_t countera = 0;
        uint32_t counterb = 0;
        for (uint32_t k=0; k<Nmax; k++) {
            if (WMCC[CCMatches[i].i][k] > WMCC[CCMatches[i].i][CCMatches[i].j]) countera++;
            if (WMCC[k][CCMatches[i].j] > WMCC[CCMatches[i].i][CCMatches[i].j]) counterb++;
        }
        if (vicinity_condition(countera, counterb, T)) {
            pair_in[i] = true;
            nMatches++;
        }
    }
    sstm << "Number of valid Pairs: " << nMatches;
    logg(sstm);
    sstm << "matched: " << 2*nMatches << ", unmatched: " << npplqc-2*nMatches << " (failed vc-QC: " << 2*(Nmin-nMatches) << ")";
    notice(sstm);
    Pairs = new struct RELATIVES[nMatches];
    for (uint32_t i=0,j=0; i<Nmin; i++) {
        if (pair_in[i]) {
            Pairs[j].i = CtrlsMap[CCMatches[i].i];
            Pairs[j].j = CasesMap[CCMatches[i].j];
            j++;
        }
    }
    delete[] CtrlsMap; CtrlsMap=NULL;
    delete[] CasesMap; CasesMap=NULL;
    delete[] PPLMapInv; PPLMapInv=NULL;
    delete_2dim(WMCC);
    delete[] CCMatches; CCMatches=NULL;
    delete[] pair_in; pair_in=NULL;
    return nMatches;
}


uint32_t call_matching_vicinity_groups1(uint64_t*** BinPPLs, uint32_t nwordsPPLs, struct PERSON* person, uint32_t nlinestfam, uint32_t nctrlsqc, uint32_t ncasesqc, struct MAP* map, uint32_t begin, uint32_t length, struct CLUSTERS*& Groups) {
    uint32_t nGroups = 0;
    #if PARALLEL
    uint32_t thread = omp_get_thread_num();
    #else
    uint32_t thread = 0;
    #endif
    stringstream sstm;
    uint32_t npplqc = nctrlsqc + ncasesqc;
    uint32_t Nmin = min(ncasesqc, nctrlsqc);
    uint32_t Nmax = max(ncasesqc, nctrlsqc);
    uint8_t T  = (uint8_t)ceil(log(nctrlsqc*ncasesqc));
    uint32_t* CtrlsMap = new uint32_t[nctrlsqc];
    uint32_t* CasesMap = new uint32_t[ncasesqc];
    int32_t* PPLMapInv = new int32_t[nlinestfam];
    memset(PPLMapInv, -1, nlinestfam*sizeof(int32_t));
    for (uint32_t h=0,i=0,j=0,k=0; k<nlinestfam; k++) {
        if (!person[k].qcin) continue;
        switch (person[k].aff[thread]) {
            case 1: CtrlsMap[j++] = k; break;
            case 2: CasesMap[i++] = k; break;
        }
        PPLMapInv[k] = h++;
    }
    uint64_t* MatchingFlags = getMatchingIn_bin_ind_gen(map, begin, length, nwordsPPLs);
    uint32_t** WMCC = getWeightMatrix_CC(BinPPLs, CtrlsMap, nctrlsqc, CasesMap, ncasesqc, PPLMapInv, MatchingFlags, nwordsPPLs);
    delete[] MatchingFlags; MatchingFlags=NULL;
    struct RELATIVES* CCMatches = maximum_weight_matching((int32_t**)WMCC, nctrlsqc, ncasesqc);
    int32_t* GroupAffil = new int32_t[nlinestfam];
    memset(GroupAffil, 0, nlinestfam*sizeof(int32_t));
    uint32_t unmatched = npplqc;
    for (uint32_t i=0; i<Nmin; i++) {
        uint32_t countera = 0;
        uint32_t counterb = 0;
        for (uint32_t k=0; k<Nmax; k++) {
            if (WMCC[CCMatches[i].i][k] > WMCC[CCMatches[i].i][CCMatches[i].j]) countera++;
            if (WMCC[k][CCMatches[i].j] > WMCC[CCMatches[i].i][CCMatches[i].j]) counterb++;
        }
        if (vicinity_condition(countera, counterb, T)) {
            unmatched -= 2;
            nGroups++;
            GroupAffil[CtrlsMap[CCMatches[i].i]] = nGroups;
            GroupAffil[CasesMap[CCMatches[i].j]] = nGroups;
        }
    }
    sstm << "Number of Groups: " << nGroups;
    logg(sstm);
    sstm << "matched: " << npplqc-unmatched << ", unmatched: " << unmatched << " (failed vc-QC: " << 2*(Nmin-nGroups) << ")";
    notice(sstm);
    delete[] CCMatches; CCMatches=NULL;
    uint32_t** WMbak = new uint32_t*[Nmax];
    if (!WMbak) die("memory allocation error for WMbak in call_vicinity_matching2()");
    WMbak[0] = new uint32_t[Nmax*Nmax];
    if (!WMbak[0]) die("memory allocation error for WMbak[0] in call_vicinity_matching2()");
    for (uint32_t i=1; i<Nmax; i++) WMbak[i] = WMbak[i-1]+Nmax;
    memcpy(WMbak[0], WMCC[0], Nmax*Nmax*sizeof(uint32_t));
    uint32_t oldunmatched = npplqc;
    while (oldunmatched > unmatched && unmatched > 0) {

        oldunmatched = unmatched;
        for (uint32_t i=0; i<nctrlsqc; i++) {
            for (uint32_t j=0; j<ncasesqc; j++) {
                if ((bool)GroupAffil[CtrlsMap[i]] ^ (bool)GroupAffil[CasesMap[j]]) WMCC[i][j] = WMbak[i][j];
                else WMCC[i][j] = 0;
            }
        }
        CCMatches = maximum_weight_matching((int32_t**)WMCC, nctrlsqc, ncasesqc);

        uint32_t failed = 0;
        for (uint32_t n=0; n<Nmin; n++) {
            if (WMCC[CCMatches[n].i][CCMatches[n].j]==0) continue;
            uint32_t countera = 0;
            uint32_t counterb = 0;
            for (uint32_t k=0; k<Nmax; k++) {
                if (WMbak[CCMatches[n].i][k] > WMCC[CCMatches[n].i][CCMatches[n].j]) countera++;
                if (WMbak[k][CCMatches[n].j] > WMCC[CCMatches[n].i][CCMatches[n].j]) counterb++;
            }
            if (vicinity_condition(countera, counterb, T)) {
                unmatched--;
                if (!GroupAffil[CtrlsMap[CCMatches[n].i]]) GroupAffil[CtrlsMap[CCMatches[n].i]] = GroupAffil[CasesMap[CCMatches[n].j]];
                if (!GroupAffil[CasesMap[CCMatches[n].j]]) GroupAffil[CasesMap[CCMatches[n].j]] = GroupAffil[CtrlsMap[CCMatches[n].i]];
            } else failed++;
        }
        sstm << "matched: " << npplqc-unmatched << ", unmatched: " << unmatched << " (failed vc-QC: " << failed << ")";
        notice(sstm);
        delete[] CCMatches; CCMatches=NULL;
    }
    for (uint32_t i=0; i<nlinestfam; i++) person[i].clusterAffil = GroupAffil[i]-1;
    delete[] GroupAffil; GroupAffil=NULL;
    Groups = Affil2Cluster(nGroups, person, nlinestfam);
    delete[] CtrlsMap; CtrlsMap=NULL;
    delete[] CasesMap; CasesMap=NULL;
    delete[] PPLMapInv; PPLMapInv=NULL;
    delete_2dim(WMCC);
    delete_2dim(WMbak);
    return nGroups;
}


uint32_t call_matching_vicinity_groups2(uint64_t*** BinPPLs, uint32_t nwordsPPLs, struct PERSON* person, uint32_t nlinestfam, uint32_t nctrlsqc, uint32_t ncasesqc, struct MAP* map, uint32_t begin, uint32_t length, struct CLUSTERS*& Groups) {
    uint32_t nGroups = 0;
    #if PARALLEL
    uint32_t thread = omp_get_thread_num();
    #else
    uint32_t thread = 0;
    #endif
    stringstream sstm;
    uint32_t npplqc = nctrlsqc + ncasesqc;
    uint32_t Nmin = min(ncasesqc, nctrlsqc);
    uint32_t Nmax = max(ncasesqc, nctrlsqc);
    uint32_t fac = (uint32_t)ceil(((double)Nmax-EPS)/Nmin);
    uint32_t dim = fac*Nmin;
    uint8_t T  = (uint8_t)ceil(log(nctrlsqc*ncasesqc));
    uint32_t* CtrlsMap = new uint32_t[nctrlsqc];
    uint32_t* CasesMap = new uint32_t[ncasesqc];
    int32_t* PPLMapInv = new int32_t[nlinestfam];
    memset(PPLMapInv, -1, nlinestfam*sizeof(int32_t));
    for (uint32_t h=0,i=0,j=0,k=0; k<nlinestfam; k++) {
        if (!person[k].qcin) continue;
        switch (person[k].aff[thread]) {
            case 1: CtrlsMap[j++] = k; break;
            case 2: CasesMap[i++] = k; break;
        }
        PPLMapInv[k] = h++;
    }
    uint64_t* MatchingFlags = getMatchingIn_bin_ind_gen(map, begin, length, nwordsPPLs);
    uint32_t** WMCC = getWeightMatrix_CC(BinPPLs, CtrlsMap, nctrlsqc, CasesMap, ncasesqc, PPLMapInv, MatchingFlags, nwordsPPLs);
    delete[] MatchingFlags; MatchingFlags=NULL;
    uint32_t** WMmult = new uint32_t*[dim];
    if (!WMmult) die("memory allocation error for WMmult in call_vicinity_matching3()");
    WMmult[0] = new uint32_t[dim*dim];
    if (!WMmult[0]) die("memory allocation error for WMmult[0] in call_vicinity_matching3()");
    for (uint32_t i=1; i<dim; i++) WMmult[i] = WMmult[i-1]+dim;
    memset(WMmult[0], 0, dim*dim*sizeof(uint32_t));
    uint32_t ni = (nctrlsqc==Nmax ? Nmax : dim);
    uint32_t nj = (ncasesqc==Nmax ? Nmax : dim);
    for (uint32_t i=0; i<ni; i++) {
        for (uint32_t j=0; j<nj; j++) {
            WMmult[i][j] = WMCC[i%nctrlsqc][j%ncasesqc];
        }
    }
    for (uint32_t i=nctrlsqc; i<dim; i++) {
        for (uint32_t j=ncasesqc; j<dim; j++) {
            WMmult[i][j] = 128*nwordsPPLs + 1;
        }
    }
    struct RELATIVES* CCMatches = maximum_weight_matching((int32_t**)WMmult, ni, nj);
    for (uint32_t n=0; n<Nmax; n++) {
        CCMatches[n].i = CCMatches[n].i%nctrlsqc;
        CCMatches[n].j = CCMatches[n].j%ncasesqc;
    }
    int32_t* GroupAffil = new int32_t[nlinestfam];
    memset(GroupAffil, 0, nlinestfam*sizeof(int32_t));
    uint32_t nPairs = 0;
    uint32_t unmatched = npplqc;
    uint32_t failed = 0;
    for (uint32_t i=0; i<Nmax; i++) {
        uint32_t countera = 0;
        uint32_t counterb = 0;
        for (uint32_t k=0; k<Nmax; k++) {
            if (WMmult[CCMatches[i].i][k] > WMmult[CCMatches[i].i][CCMatches[i].j]) countera++;
            if (WMmult[k][CCMatches[i].j] > WMmult[CCMatches[i].i][CCMatches[i].j]) counterb++;
        }
        if (vicinity_condition(countera, counterb, T)) {
            unmatched--;
            nPairs++;
            if (!GroupAffil[CtrlsMap[CCMatches[i].i]] && !GroupAffil[CasesMap[CCMatches[i].j]]) {
                unmatched--;
                nGroups++;
                GroupAffil[CtrlsMap[CCMatches[i].i]] = nGroups;
                GroupAffil[CasesMap[CCMatches[i].j]] = nGroups;
            }
            else if (GroupAffil[CasesMap[CCMatches[i].j]]) GroupAffil[CtrlsMap[CCMatches[i].i]] = GroupAffil[CasesMap[CCMatches[i].j]];
            else if (GroupAffil[CtrlsMap[CCMatches[i].i]]) GroupAffil[CasesMap[CCMatches[i].j]] = GroupAffil[CtrlsMap[CCMatches[i].i]];
        } else {
            if (WMCC[CCMatches[i].i][CCMatches[i].j] > 0) {
                failed+=2;
            }
        }
    }
    sstm << "Number of Groups: " << nGroups << " (Number of valid Pairs: " << nPairs << ")";
    logg(sstm);
    sstm << "matched: " << npplqc-unmatched << ", unmatched: " << unmatched << " (failed vc-QC: " << failed << ")";
    notice(sstm);
    delete[] CCMatches; CCMatches=NULL;
    memset(WMmult[0], 0, dim*dim*sizeof(uint32_t));
    for (uint32_t i=0; i<nctrlsqc; i++) {
        for (uint32_t j=0; j<ncasesqc; j++) {
            WMmult[i][j] = WMCC[i][j];
        }
    }
    uint32_t oldunmatched = npplqc;
    while (oldunmatched > unmatched && unmatched > 0) {
        oldunmatched = unmatched;

        for (uint32_t i=0; i<nctrlsqc; i++) {
            for (uint32_t j=0; j<ncasesqc; j++) {
                if ((bool)GroupAffil[CtrlsMap[i]] ^ (bool)GroupAffil[CasesMap[j]]) WMCC[i][j] = WMmult[i][j];
                else WMCC[i][j] = 0;
            }
        }
        CCMatches = maximum_weight_matching((int32_t**)WMCC, nctrlsqc, ncasesqc);

        failed = 0;
        for (uint32_t n=0; n<Nmin; n++) {
            if (WMCC[CCMatches[n].i][CCMatches[n].j]==0) continue;
            uint32_t countera = 0;
            uint32_t counterb = 0;
            for (uint32_t k=0; k<Nmax; k++) {
                if (WMmult[CCMatches[n].i][k] > WMCC[CCMatches[n].i][CCMatches[n].j]) countera++;
                if (WMmult[k][CCMatches[n].j] > WMCC[CCMatches[n].i][CCMatches[n].j]) counterb++;
            }
            if (vicinity_condition(countera, counterb, T)) {
                unmatched--;
                if (!GroupAffil[CtrlsMap[CCMatches[n].i]]) GroupAffil[CtrlsMap[CCMatches[n].i]] = GroupAffil[CasesMap[CCMatches[n].j]];
                if (!GroupAffil[CasesMap[CCMatches[n].j]]) GroupAffil[CasesMap[CCMatches[n].j]] = GroupAffil[CtrlsMap[CCMatches[n].i]];
            } else failed++;
        }
        sstm << "matched: " << npplqc-unmatched << ", unmatched: " << unmatched << " (failed vc-QC: " << failed << ")";
        notice(sstm);
        delete[] CCMatches; CCMatches=NULL;
    }
    for (uint32_t i=0; i<nlinestfam; i++) person[i].clusterAffil = GroupAffil[i]-1;
    delete[] GroupAffil; GroupAffil=NULL;
    Groups = Affil2Cluster(nGroups, person, nlinestfam);
    delete[] CtrlsMap; CtrlsMap=NULL;
    delete[] CasesMap; CasesMap=NULL;
    delete[] PPLMapInv; PPLMapInv=NULL;
    delete_2dim(WMCC);
    delete_2dim(WMmult);
    return nGroups;
}


struct RELATIVES* maximum_weight_matching(int32_t** A, uint32_t nrows, uint32_t ncols) {
  uint32_t* matching_hungarian(int32_t**, uint32_t);
  uint32_t Nmax = max(ncols, nrows);
  uint32_t Nmin = min(ncols, nrows);
  uint32_t* X = matching_hungarian(A, Nmax);
  struct RELATIVES* Matches = new struct RELATIVES[Nmin];
  if (ncols > nrows) {
    uint32_t* Y = new uint32_t[Nmax];
    for (uint32_t y=0; y<Nmax; y++) Y[X[y]] = y;
    for (uint32_t x=0; x<Nmin; x++) {
      Matches[x].i = x;
      Matches[x].j = Y[x];
    }
    delete[] Y; Y=NULL;
  } else {
    for (uint32_t y=0; y<Nmin; y++) {
      Matches[y].i = X[y];
      Matches[y].j = y;
    }
  }
  delete[] X; X=NULL;
  return Matches;
}


uint32_t clustering_hungarian(int32_t** M, uint32_t ndim, uint32_t*& Affil) {
	uint32_t* matching_hungarian(int32_t**, uint32_t);
    extern uint8_t t_all;
    uint8_t T = (t_all ? t_all : (uint8_t)(floor(log(ndim))));
    uint32_t nNewClusters = ndim;
    uint32_t nOldClusters = ndim;
    uint32_t* X;
    Affil = new uint32_t[ndim];
    for (uint32_t i=0; i<ndim; i++) Affil[i] = i;
    int32_t** Dist = new int32_t*[ndim];
    if (!Dist) die("memory allocation error for Dist in clustering_hungarian()");
    Dist[0] = new int32_t[ndim*ndim];
    if (!Dist[0]) die("memory allocation error for Dist[0] in clustering_hungarian()");
    for (uint32_t i=1; i<ndim; i++) Dist[i] = Dist[i-1]+ndim;
    struct RELATIVES** c0 = new struct RELATIVES*[ndim];
    if (!c0) die("memory allocation error for c0 in clustering_hungarian()");
    c0[0] = new struct RELATIVES[ndim*ndim];
    if (!c0[0]) die("memory allocation error for c0[0] in clustering_hungarian()");
    for (uint32_t i=1; i<ndim; i++) c0[i] = c0[i-1]+ndim;
    int32_t* ClusterMap = new int32_t[ndim];
    do {
        memset(Dist[0], 0, ndim*ndim*sizeof(int32_t));
        memset(c0[0], -1, ndim*ndim*sizeof(struct RELATIVES));
        memset(ClusterMap, -1, ndim*sizeof(int32_t));
        nOldClusters = nNewClusters;
        for (uint32_t a,b,i=0; i<ndim; i++) {
            a = Affil[i];
            for (uint32_t j=0; j<ndim; j++) {
                b = Affil[j];
                if (a==b) continue;
                if (M[i][j] >= Dist[a][b]) {
                    Dist[a][b] = M[i][j];
                    c0[a][b].i = i;
                    c0[a][b].j = j;
                }
            }
        }
        for (uint32_t a=0; a<nOldClusters; a++) {
            bool found = false;
            for (uint32_t b=0; b<nOldClusters; b++) {
                if (a==b) continue;
                uint32_t countera = 0;
                uint32_t counterb = 0;
                for (uint32_t i=0; i<ndim; i++) {
                    if (Affil[i]==a && M[i][c0[a][b].i] > Dist[a][b]) countera++;
                    if (Affil[i]==b && M[c0[a][b].j][i] > Dist[a][b]) counterb++;
                }
                if (countera>T || counterb>T) Dist[a][b] = 0;
                else found = true;
            }
            if (!found) Dist[a][a] = INT_MAX;
        }
        X = matching_hungarian(Dist, nOldClusters);
        nNewClusters = 0;
        for (uint32_t b,y=0; y<nOldClusters; y++) {
            if (ClusterMap[X[y]]>-1) continue;
            ClusterMap[X[y]] = nNewClusters;
            b = y;
            do {
                for (uint32_t z=y+1; z<nOldClusters; z++) {
                    if (b==X[z]) {
                        ClusterMap[X[z]] = nNewClusters;
                        b = z;
                        break;
                    }
                }
            } while (b!=X[y]);
            nNewClusters++;
        }
        for (uint32_t i=0; i<ndim; i++) Affil[i] = ClusterMap[Affil[i]];
        delete[] X; X=NULL;
    } while (nNewClusters < nOldClusters);
    delete[] ClusterMap; ClusterMap=NULL;
    delete_2dim(c0);
    delete_2dim(Dist);
    return nNewClusters;
}


void update_slack(int32_t** W, int32_t n, int32_t* lx, int32_t* ly, int32_t* slack, int32_t* slackx, int32_t x) {
  for (int32_t y = 0; y < n; y++) {
    if (lx[x] + ly[y] - W[x][y] < slack[y]) {
      slack[y] = lx[x] + ly[y] - W[x][y];
      slackx[y] = x;
    }
  }
}


void augment(int32_t** W, int32_t n, int32_t* Y, int32_t* X, bool* S, bool* T, int32_t* lx, int32_t* ly, int32_t* slack, int32_t* slackx, int32_t* prev, int32_t* q, int32_t max_match) {
  int32_t wr = 0;
  int32_t rd = 0;
  int32_t root = 0;
  int32_t x,y;
  memset(S, false, n*sizeof(bool));
  memset(T, false, n*sizeof(bool));
  memset(prev, -1, n*sizeof(int32_t));
  for (x=0; x<n; x++) {
    if (Y[x] == -1) {
      root = x;
      q[wr++] = x;
      prev[x] = -2;
      S[x] = true;
      break;
    }
  }
  for (y=0; y<n; y++) {
    slack[y] = lx[root] + ly[y] - W[root][y];
    slackx[y] = root;
  }
  while (true) {
    while (rd < wr) {
      x = q[rd++];
      for (y=0; y<n; y++) {
        if (W[x][y] == lx[x] + ly[y] && !T[y]) {
          if (X[y] == -1) break;
          q[wr++] = X[y];
          T[y] = true;
          S[X[y]] = true;
          prev[X[y]] = x;
          update_slack(W,n,lx,ly,slack,slackx,X[y]);
        }
      }
      if (y < n) break;
    }
    if (y < n) break;
    int32_t delta = INT_MAX;
    for (int32_t i=0; i<n; i++)
      if (!T[i]) delta = min(delta, slack[i]);
    for (int32_t i=0; i<n; i++) {
      if (S[i]) lx[i] -= delta;
      if (T[i]) ly[i] += delta;
      else   slack[i] -= delta;
    }
    wr = rd = 0;
    for (y=0; y<n; y++) {
      if (!T[y] &&  slack[y] == 0) {
        if (X[y] == -1) {
          x = slackx[y];
          break;
        } else {
          T[y] = true;
          if (!S[X[y]]) {
            q[wr++] = X[y];
            S[X[y]] = true;
            prev[X[y]] = slackx[y];
            update_slack(W,n,lx,ly,slack,slackx,X[y]);
          }
        }
      }
    }
    if (y < n) break;
  }
  if (y < n) {
    max_match++;
    for (int32_t cx=x, cy=y, ty; cx!=-2; cx=prev[cx], cy=ty) {
      ty = Y[cx];
      X[cy] = cx;
      Y[cx] = cy;
    }
    if (max_match==n) return;
    augment(W,n,Y,X,S,T,lx,ly,slack,slackx,prev,q,max_match);
  }
}


uint32_t* matching_hungarian(int32_t** W, uint32_t n) {
  int32_t* X      = new int32_t[n];
  int32_t* Y      = new int32_t[n];
  bool* S         = new bool[n];
  bool* T         = new bool[n];
  int32_t* lx     = new int32_t[n];
  int32_t* ly     = new int32_t[n];
  int32_t* slack  = new int32_t[n];
  int32_t* slackx = new int32_t[n];
  int32_t* prev   = new int32_t[n];
  int32_t* q      = new int32_t[n];
  int32_t max_match = 0;
  memset(X, -1, n*sizeof(int32_t));
  memset(Y, -1, n*sizeof(int32_t));
  memset(lx, 0, n*sizeof(int32_t));
  memset(ly, 0, n*sizeof(int32_t));
  for (uint32_t x = 0; x < n; x++)
    for (uint32_t y = 0; y < n; y++)
      lx[x] = max(lx[x], W[x][y]);
  augment(W,n,Y,X,S,T,lx,ly,slack,slackx,prev,q,max_match);
  delete[] q; q=NULL;
  delete[] prev; prev=NULL;
  delete[] slackx; slackx=NULL;
  delete[] slack; slack=NULL;
  delete[] ly; ly=NULL;
  delete[] lx; lx=NULL;
  delete[] T; T=NULL;
  delete[] S; S=NULL;
  delete[] Y; Y=NULL;
  return (uint32_t*)X;
}


uint32_t call_matchingRARE(uint64_t*** BinPPLs, uint32_t nwordsPPLs, int* SNPMapInv, int* PPLMapInv, struct PERSON* person, uint32_t nlinestfam, int nctrls, int ncases, struct MAP* map, int nlinestped, struct REGION* chrPositions, double rarefreq, uint32_t minMatchingVariantsPerWindow, struct LOCAL_MATCHING*& matchingWindows) {
    if (!nwordsPPLs) return 0;
    #if PARALLEL
    uint32_t thread = omp_get_thread_num();
    #else
    uint32_t thread = 0;
    #endif
    stringstream sstm;
    sstm << "\nPerforming local matchings on common SNPs, f>" << rarefreq << ", minimum SNP per Window " << minMatchingVariantsPerWindow << ".";
    logg(sstm);
    uint64_t* RareFlags = getRareIn_bin_ind_gen(map, 0, nlinestped, nwordsPPLs, rarefreq);
    uint32_t ntotMatchingWindows = 0;
    for (uint8_t c=0; c<27; c++) {
        if (chrPositions[c].begin<0 || chrPositions[c].end<0) continue;
        uint32_t beg_nr  = SNPMapInv[chrPositions[c].begin]/64;
        uint32_t end_nr  = SNPMapInv[chrPositions[c].end  ]/64;
        uint32_t beg_pos = SNPMapInv[chrPositions[c].begin]%64;
        uint32_t end_pos = SNPMapInv[chrPositions[c].end  ]%64;
        uint32_t nMatchingVariants = 0;
        uint64_t Zbeg=0ull,Zend=0ull;
        if (beg_nr==end_nr) {
            for (uint8_t i=beg_pos; i<=end_pos; i++) Zbeg |= mask[i];
            nMatchingVariants += bitcount64(~RareFlags[beg_nr] & Zbeg);
        } else {
            for (uint8_t i=beg_pos; i<64; i++) Zbeg |= mask[i];
            for (uint8_t i=0; i<=end_pos; i++) Zend |= mask[i];
            nMatchingVariants += bitcount64(~RareFlags[beg_nr] & Zbeg);
            for (uint32_t w=beg_nr+1; w<end_nr; w++) nMatchingVariants += bitcount64(~RareFlags[w]);
            nMatchingVariants += bitcount64(~RareFlags[end_nr] & Zend);
        }
        uint32_t nMatchingWindows = nMatchingVariants/minMatchingVariantsPerWindow;
        if (!nMatchingWindows) {
          nMatchingWindows = 1;
          sstm << "Chromosome " << c+0 << " does not comprise a sufficient number of common SNPs (" << nMatchingVariants << "). Proceeding anyway with one window.";
          error(sstm);
        }
        uint32_t nMatchingVariantsPerWindow = nMatchingVariants/nMatchingWindows;
        uint32_t nRestVariantsLastWindow    = nMatchingVariants%nMatchingWindows;
        struct LOCAL_MATCHING* tmpMatchingWindows = new struct LOCAL_MATCHING[ntotMatchingWindows+nMatchingWindows];
        copy(matchingWindows, matchingWindows+ntotMatchingWindows, tmpMatchingWindows);
        memset(tmpMatchingWindows+ntotMatchingWindows, 0, nMatchingWindows*sizeof(struct LOCAL_MATCHING));
        if (matchingWindows!=NULL) { delete[] matchingWindows; matchingWindows=NULL; }
        matchingWindows = tmpMatchingWindows;
        uint32_t nr=beg_nr,pos=beg_pos;
        uint64_t Z=Zbeg;
        for (uint32_t w=ntotMatchingWindows; w<ntotMatchingWindows+nMatchingWindows; w++) {
            uint32_t i=0;
            if (w==ntotMatchingWindows+nMatchingWindows-1) nMatchingVariantsPerWindow += nRestVariantsLastWindow;
            matchingWindows[w].chr = c;
            matchingWindows[w].bin_start = 64*nr+pos;
            i += bitcount64(~RareFlags[nr] & Z);
            while (i<nMatchingVariantsPerWindow) i += bitcount64(~RareFlags[++nr]);
            pos = 63;
            Z = 0ull;
            while (i>nMatchingVariantsPerWindow) {
                if (~RareFlags[nr] & mask[pos]) i--;
                Z |= mask[pos--];
            }
            matchingWindows[w].bin_stop = 64*nr+pos;
            if (++pos==64) { pos=0; nr++; Z=-1ull; }
            if (w==ntotMatchingWindows+nMatchingWindows-1) matchingWindows[w].bin_stop = SNPMapInv[chrPositions[c].end];
        }
        ntotMatchingWindows += nMatchingWindows;
    }
    uint32_t nppl = nctrls + ncases;
    uint32_t Nmin = min(ncases, nctrls);
    uint32_t Nmax = max(ncases, nctrls);
    uint8_t T  = (uint8_t)ceil(log(nctrls*ncases));
    uint32_t* CtrlsMap = new uint32_t[nctrls];
    uint32_t* CasesMap = new uint32_t[ncases];
    for (uint32_t i=0,j=0,k=0; k<(uint32_t)nlinestfam; k++) {
        if (!person[k].qcin) continue;
        switch (person[k].aff[thread]) {
            case 1: CtrlsMap[j++] = k; break;
            case 2: CasesMap[i++] = k; break;
        }
    }
    uint32_t** WMCC = new uint32_t*[Nmax];
    if (!WMCC) die("memory allocation error for WMCC in call_matchingRARE()");
    WMCC[0] = new uint32_t[Nmax*Nmax];
    if (!WMCC[0]) die("memory allocation error for WMCC[0] in call_matchingRARE()");
    for (uint32_t i=1; i<Nmax; i++) WMCC[i] = WMCC[i-1]+Nmax;
    uint32_t** WMbak = new uint32_t*[Nmax];
    if (!WMbak) die("memory allocation error for WMbak in call_matchingsRARE()");
    WMbak[0] = new uint32_t[Nmax*Nmax];
    if (!WMbak[0]) die("memory allocation error for WMbak[0] in call_matchingsRARE()");
    for (uint32_t i=1; i<Nmax; i++) WMbak[i] = WMbak[i-1]+Nmax;
    uint64_t* MatchingFlags = new uint64_t[nwordsPPLs];
    struct RELATIVES* CCMatches = NULL;
    int32_t* GroupAffil = new int32_t[nlinestfam];
    for (uint32_t w=0; w<ntotMatchingWindows; w++) {
        sstm << "\nObtaining matching for window " << w << " (chr" << matchingWindows[w].chr+0 << ") ... ";
        logg(sstm);
        uint32_t beg_nr  = matchingWindows[w].bin_start/64;
        uint32_t end_nr  = matchingWindows[w].bin_stop /64;
        uint32_t beg_pos = matchingWindows[w].bin_start%64;
        uint32_t end_pos = matchingWindows[w].bin_stop %64;
        memset(MatchingFlags, 0, nwordsPPLs*sizeof(uint64_t));
        memset(WMCC[0], 0, Nmax*Nmax*sizeof(uint32_t));
        uint64_t Zbeg=0ull,Zend=0ull;
        for (uint8_t i=beg_pos; i<64; i++) Zbeg |= mask[i];
        for (uint8_t i=0; i<=end_pos; i++) Zend |= mask[i];
        MatchingFlags[beg_nr] = ~RareFlags[beg_nr] & Zbeg;
        for (uint32_t i=beg_nr+1; i<end_nr; i++) MatchingFlags[i] = ~RareFlags[i];
        MatchingFlags[end_nr] = ~RareFlags[end_nr] & Zend;
        for (uint32_t i=0; i<(uint32_t)nctrls; i++) {
            for (uint32_t j=0; j<(uint32_t)ncases; j++) {
                WMCC[i][j] = (uint32_t)( 128.0*nwordsPPLs*getPartialIbs(BinPPLs, PPLMapInv[CtrlsMap[i]], PPLMapInv[CasesMap[j]], beg_nr, end_nr, MatchingFlags) + 0.5 );
            }
        }
        CCMatches = maximum_weight_matching((int32_t**)WMCC, nctrls, ncases);
        uint32_t nGroups = 0;
        memset(GroupAffil, 0, nlinestfam*sizeof(int32_t));
        uint32_t unmatched = nppl;
        for (uint32_t i=0; i<Nmin; i++) {
            uint32_t countera = 0;
            uint32_t counterb = 0;
            for (uint32_t k=0; k<Nmax; k++) {
                if (WMCC[CCMatches[i].i][k] > WMCC[CCMatches[i].i][CCMatches[i].j]) countera++;
                if (WMCC[k][CCMatches[i].j] > WMCC[CCMatches[i].i][CCMatches[i].j]) counterb++;
            }
            if (vicinity_condition(countera, counterb, T)) {
                unmatched -= 2;
                nGroups++;
                GroupAffil[CtrlsMap[CCMatches[i].i]] = nGroups;
                GroupAffil[CasesMap[CCMatches[i].j]] = nGroups;
            }
        }
        sstm << "Number of Groups: " << nGroups; logg(sstm);
        delete[] CCMatches; CCMatches=NULL;
        memcpy(WMbak[0], WMCC[0], Nmax*Nmax*sizeof(uint32_t));
        uint32_t oldunmatched = nppl;
        while (oldunmatched > unmatched && unmatched > 0) {
            oldunmatched = unmatched;
            for (int32_t i=0; i<nctrls; i++) {
                for (int32_t j=0; j<ncases; j++) {
                    if ((bool)GroupAffil[CtrlsMap[i]] ^ (bool)GroupAffil[CasesMap[j]]) WMCC[i][j] = WMbak[i][j];
                    else WMCC[i][j] = 0;
                }
            }
            CCMatches = maximum_weight_matching((int32_t**)WMCC, nctrls, ncases);
            uint32_t failed = 0;
            for (uint32_t n=0; n<Nmin; n++) {
                if (WMCC[CCMatches[n].i][CCMatches[n].j]==0) continue;
                uint32_t countera = 0;
                uint32_t counterb = 0;
                for (uint32_t k=0; k<Nmax; k++) {
                    if (WMbak[CCMatches[n].i][k] > WMCC[CCMatches[n].i][CCMatches[n].j]) countera++;
                    if (WMbak[k][CCMatches[n].j] > WMCC[CCMatches[n].i][CCMatches[n].j]) counterb++;
                }
                if (vicinity_condition(countera, counterb, T)) {
                    unmatched--;
                    if (!GroupAffil[CtrlsMap[CCMatches[n].i]]) GroupAffil[CtrlsMap[CCMatches[n].i]] = GroupAffil[CasesMap[CCMatches[n].j]];
                    if (!GroupAffil[CasesMap[CCMatches[n].j]]) GroupAffil[CasesMap[CCMatches[n].j]] = GroupAffil[CtrlsMap[CCMatches[n].i]];
                } else failed++;
            }
            delete[] CCMatches; CCMatches=NULL;
        }
        for (uint32_t i=0; i<nlinestfam; i++) person[i].clusterAffil = GroupAffil[i]-1;
        matchingWindows[w].groups  = Affil2Cluster(nGroups, person, nlinestfam);
        matchingWindows[w].nGroups = nGroups;
    }
    for (uint32_t i=0; i<nlinestfam; i++) person[i].clusterAffil = 0;
    delete_2dim(WMCC);
    delete_2dim(WMbak);
    delete[] MatchingFlags; MatchingFlags=NULL;
    delete[] GroupAffil; GroupAffil=NULL;
    delete[] CasesMap; CasesMap=NULL;
    delete[] CtrlsMap; CtrlsMap=NULL;
    delete[] RareFlags; RareFlags=NULL;
    return ntotMatchingWindows;
}


struct LOCAL_MATCHING* get_matchingRARE_window(uint8_t chr, uint32_t pos, struct LOCAL_MATCHING* matchingRARE, uint32_t nMWindowsRARE) {
    int32_t dp,dm;
    for (uint32_t w=0; w<nMWindowsRARE; w++) {
        if (matchingRARE[w].chr != chr) continue;
        dp = (matchingRARE[w].bin_stop+matchingRARE[w].bin_start);
        dm = (matchingRARE[w].bin_stop-matchingRARE[w].bin_start);
        if (abs(dp-2*(int32_t)pos) <= dm) return &matchingRARE[w];
    }
    stringstream sstm;
    sstm << "No matching window found for chromosome " << chr << " bin-position " << pos;
    error(sstm);
    return NULL;
}

