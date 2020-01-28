#include <stdio.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <set>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits>

using namespace std;




struct WEIGHTS
{
  int mode;
  int betapar1;
  int betapar2;
  float logpar1;
  float logpar2;
};

double *weightvec=NULL;

struct STATplus regGeneral(int x[27], //in this function just a dummy
			   int *cov, //indicates if covar is is used;
			   int sexcov,
			   struct PERSON *person, int nlinestfam,int ncases, int ncontrols, int nrest,
			   int a, int b, int c, // here just dummies
			   double inflationfactor, //dummy
			   double casecounts5[3][3][3][3][3], //int this function just a dummy
			   double *p, double **newbeta,
			   double **X, double **Xmod, double **Xt, double **A, double **UNNT,double **VNN,
			   double *S, double **Sinv, double **A0, double **Ainv, double **AinvXt, double *Y,
			   double **Yminusp,
			   int N, //dummy
			   int alt, double *Yhelp, int xtype, int female,
			   int male,
			   struct COUNTS counts, //dummy
			   double **Yt, double **YtX, double **YtXAinv,
			   double **YtXAinvXt,
			   int haplo, //dummy
			   double **D, double **T, double **U, double **Ut, double **sumPP,
			   double **sumPJ, double **sumPK, double **MMinv, int test,int thread,
			   char *hapfile, int dohapfile, char *hapstring, FILE *file1, //dummies
			   int skip, //dummy
			   int npplqc, int* PPLMap, uint64_t*** BinSNPs, struct PPLLOCATION* PplLocations,
			   int nMc,struct STATplus result,
			   int caseOnly, int secondSystem, //dummies
			   int numberOfAllCov,int ncov,
			   int liabilityCut, int singleMarkerTest, int covNum, //dummies
			   //new parameters:
			   int nSnps, //number of snps to be included with additive term
			   int nSnpsDom, //number of snps to be included with dominance term
			   int nEnvirons, // number of non-genetic parameters (for future extensions)
			   int nInters, // number of interaction terms to be used (under construction)
			   int *snps, // list of snps to be included with additive term, indices in BinSNPs
			   int *snpsDom, // list of snps to be included with dominance term, indices in BinSNPs
			   int *environs, // list of non-genetic parameters (for future extensions)
			   int *inters, // list of interaction terms to be used (coding under construction)
			   int nx, // total number of parameters in regression model. n=nSnps+nSnpsDom+nEnvirons+nInters+sexcov
                	   int dim1, // maximum number of paramters allowed according to allocation in main
			   int collapseRare, // =1: do Morris rare-frac-method
			   int collcollapseRare, // =1: do Morris rare-coll-method
			   int qt,
			   double *weightvec,
			   int firstbinlastSNP
			   );

struct STATplus logRegGeneral(int x[27], //in this function just a dummy
			      int *cov, //indicates if covar is isused;
			      int sexcov,
			      struct PERSON *person, int nlinestfam,int ncases, int ncontrols, int nrest,
			      int a, int b, int c, // here just dummies
			      double inflationfactor, //dummy
			      double casecounts5[3][3][3][3][3], //int this function just a dummy
			      double *p, double **newbeta,
			      double **X, double **Xmod, double **Xt, double **A, double **UNNT,double **VNN,
			      double *S, double **Sinv, double **A0, double **Ainv, double **AinvXt, double *Y,
			      double **Yminusp,
			      int N, //dummy
			      int alt, double *Yhelp, int xtype, int female,
			      int male,
			      struct COUNTS counts, //dummy
			      double **Yt, double **YtX, double **YtXAinv,
			      double **YtXAinvXt,
			      int haplo, //dummy
			      double **D, double **T, double **U, double **Ut, double **sumPP,
			      double **sumPJ, double **sumPK, double **MMinv, int test,int thread,
			      char *hapfile, int dohapfile, char *hapstring, FILE *file1, //dummies
			      int skip, //dummy
			      int npplqc, int* PPLMap, uint64_t*** BinSNPs, struct PPLLOCATION* PplLocations,
			      int nMc,struct STATplus result,
			      int caseOnly, int secondSystem, //dummies
			      int numberOfAllCov,int ncov,
			      int liabilityCut, int singleMarkerTest, int covNum, //dummies
			      //new parameters:
			      int nSnps, //number of snps to be included with additive term
			      int nSnpsDom, //number of snps to be included with dominance term
			      int nEnvirons, // number of non-genetic parameters (for future extensions)
			      int nInters, // number of interaction terms to be used (under construction)
			      int *snps, // list of snps to be included with additive term, indices in BinSNPs
			      int *snpsDom, // list of snps to be included wit dominance term, indices in BinSNPs
			      int *environs, // list of non-genetic parameters (for future extensions)
			      int *inters, // list of interaction terms to be used (coding under construction)
			      int nx, // total number of parameters in regression model. n=nSnps+nSnpsDom+nEnvirons+nInters+sexcov
			      int dim1, // maximum number of paramters allowed according to allocation in main
			      int collapseRare, // =1: do Morris rare-method
			      int collcollapseRare, // =1: do Morris other rare-method
			      int qt,
			      double *weightvec,
			      int covariancematrix,  //new logRegGeneral
			      int firstbinlastSNP,
			      int collinter
			      );

void resultCopy(struct STATplus *result1, struct STATplus result2);

void init(struct STATplus *result, int haplo, int n, int numberOfAllCov, int multi, int needArrays, int dim1, int dim2);

struct STATplus logreg( int x[27], //indicates which parameters are used
			int *cov, //indicates if covar is isused;
			int sexcov, //indicates if only males or females are counted; inidcates if sex is a covariate
			struct PERSON *person, int nlinestfam, int ncases, int ncontrols, int nrest, //person info, helpful numbers
			int a, int b, int c, // data matrix; a,b,c identify snps 1, 2 and 3
			double inflationfactor,
			double casecounts5[3][3][3][3][3],double controlcounts5[3][3][3][3][3], //3-dimgenocounts
			double *p, double **newbeta, //likelihood per ind , betas
			double **X, double **Xmod, double **Xt, double **A,  double **UNNT,double **VNN, //empty matrices, vectors
			double *S, double **Sinv, double **A0, double **Ainv, double **AinvXt, double *Y, //empty matrices, vectors
			double **Yminusp,
			int N, int alt, double *Yhelp, int xtype, int female, int male, struct COUNTS counts,
			double **D, double **T, double **U, double **Ut, double **sumPP, double **sumPJ, double **sumPK, double **MMinv,int test,int thread, int npplqc, int* PPLMap, uint64_t*** BinSNPs,
			struct PPLLOCATION* PplLocations,int nMc,struct STATplus result, int numberOfAllCov, int liabilityCut, int singleMarkerTest, int covNum,int covariancematrix, int df_L1, int df_L2, int dosage, float ***genoWeights);

// Neue qtreg-Funktion
struct STATplus qtreg( int x[27], //indicates which parameters are used
		       int *cov, //indicates if covar is isused;
		       int sexcov, //indicates if only males or females are counted; inidcates if sex is a covariate
		       struct PERSON *person, int nlinestfam,int ncases, int ncontrols, int nrest, //person info, helpful numbers
		       int a, int b, int c, // data matrix; a,b,c identify snps 1, 2 and 3
		       double inflationfactor, double casecounts5[3][3][3][3][3], //3-dimgenocounts
		       double *p, double **newbeta, //likelihood per ind , betas
		       double **X, double **Xmod, double **Xt, double **A, double **UNNT,double **VNN, //empty matrices, vectors
		       double *S, double **Sinv, double **A0, double **Ainv, double **AinvXt, double *Y, //empty matrices, vectors
		       double **Yminusp, int N, int alt, double *Yhelp, int xtype, int female,
		       int male, struct COUNTS counts, double **Yt, double **YtX, double **YtXAinv, double **YtXAinvXt, int haplo,
		       double **D, double **T, double **U, double **Ut, double **sumPP, double **sumPJ, double **sumPK, double **MMinv, int test,int thread, char *hapfile, int dohapfile, char *hapstring,
		       FILE *file1,int skip, int npplqc, int* PPLMap, uint64_t*** BinSNPs, struct PPLLOCATION* PplLocations,int nMc,struct STATplus result,int caseOnly, int secondSystem,
		       int numberOfAllCov,int liabilityCut, int singleMarkerTest,int covNum,int covariancematrix, int df_L1, int df_L2);


double betai(double a, double b, double x);


void freeResult(struct STATplus result1, int n)
{
  if(n==0)
    {
      free(result1.betaNew_se);free(result1.betaNew_lcl);free(result1.betaNew_rcl);
      free(result1.oddsRatio);free(result1.lcloddsRatio);free(result1.rcloddsRatio);
      free(result1.sigma1);
    }
  free(result1.b);free(result1.bcv);free(result1.in);
};




struct WINDOW
{
  // Pretest-Permutations
  // Tests numbered as:
  // FR CMAT COLL REG FRACREG COLLREG
  int nperm[6];

  // Number of tests (still) to be conducted
  int ntests;

  // Variant data
  int **levelpos; // relative positions of SNPs in bim file, levelwise
  int n_level;
  //  int start;
  //  int end;
  //  int start1,start2,end1,end2; // For RAREINTER

  int *n_at_level; // how many SNPs per level
  double *maf_at_level; // How high is MAF at level
};

struct WINDOW *window;

struct WINDOW_INTER
{
  // Optimal Variable threshold per test
  // Number of tests (still) to be conducted

  // Variant data
  int **levelpos; // relative positions of SNPs in bim file, levelwise
  int n_level;
  int start;
  int end;
  int start1,start2,end1,end2; // For RAREINTER

  int *n_at_level; // how many SNPs per level
  double *maf_at_level; // How high is MAF at level
};

struct WINDOW_INTER *window_inter;


int *FISHERnotconverged;
int *FISHERpretestpassed;
int *REGRESSIONpretestpassed;
int *FRACREGpretestpassed;
int *COLLREGpretestpassed;
int *CMATpretestpassed;
int *COLLpretestpassed;

// Confidence Intervals
double **FISHER_CI;
double **CMAT_CI;
double **COLL_CI;
double **REGRESSION_CI;
double **FRACREG_CI;
double **COLLREG_CI;

// betas, se
double *COLLREG_beta;
double *COLLREG_se;
double *FRACREG_beta;
double *FRACREG_se;

// optimal MAF_Level
int *MAF_Level_VT_FISHER;
int *MAF_Level_VT_CMAT;
int *MAF_Level_VT_COLL;
int *MAF_Level_VT_REGRESSION;
int *MAF_Level_VT_FRACREG;
int *MAF_Level_VT_COLLREG;


double *rarefFISHER=NULL;
double *rarefREGRESSION=NULL;
double *rarefFRACREG=NULL;
double *rarefCOLLREG=NULL;
double *rarefCOLL=NULL;
double *OR_COLL_vec=NULL;
double *OR_COLL_f_vec=NULL;
double *OR_CMAT_vec=NULL;
double *rarefCMAT=NULL;
int **windowPositions=NULL;
int **doublewindowcoord=NULL; // for RAREINTER
int **intervalfileWindows=NULL;
double *pCOLLvec=NULL;
double *pCMATvec=NULL;
double *pFISHERvec=NULL;
double *pREGRESSIONvec=NULL;
double *pFRACREGvec=NULL;
double *pCOLLREGvec=NULL;
double *pFISHERvecChi=NULL;
double **limitsStatCOLL=NULL;
double **limitsStatCMAT=NULL;
double **limitsStatFISHER=NULL;
double **limitsStatREGRESSION=NULL;
double **limitsStatFRACREG=NULL;
double **limitsStatCOLLREG=NULL;
int *nCarriers=NULL;
int *nSNPsInWindow=NULL;
int nRareSNPsPerChr[26];
int *nRareSNPs=NULL;
int *nRareSNPsCOLL=NULL;
int *nRareSNPsCMAT=NULL;
int *nRareSNPsFISHER=NULL;
int *nRareSNPsREGRESSION=NULL;
int *nRareSNPsFRACREG=NULL;
int *nRareSNPsCOLLREG=NULL;

double OR_COLL=0;
double OR_COLL_f=0;
double OR_CMAT=0;

// for RARE_TESTS
bool FISHERflag=1;
bool REGRESSIONflag=0;
bool FRACREGflag=0;
bool COLLREGflag=0;
bool COLLflag=1;
bool CMATflag=1;

// for RARE_PRETEST
int rarepretest=0;
float rarepretestlimit=0;
double rareregpretest=0;
double wilsonpretest=0;
int adaptive =0;
int currentn=0;

int *FISHERcount=NULL;
int *REGRESSIONcount=NULL;
int *FRACREGcount=NULL;
int *COLLREGcount=NULL;
int *COLLcount=NULL;
int *CMATcount=NULL;

int *completedFISHER=NULL;
int *completedREGRESSION=NULL;
int *completedFRACREG=NULL;
int *completedCOLLREG=NULL;
int *completedCOLL=NULL;
int *completedCMAT=NULL;

double *FISHERstats=NULL;
double *REGRESSIONstats=NULL;
double *FRACREGstats=NULL;
double *COLLREGstats=NULL;
double *COLLstats=NULL;
double *CMATstats=NULL;

double *FISHERpermstats=NULL;
double *REGRESSIONpermstats=NULL;
double *FRACREGpermstats=NULL;
double *COLLREGpermstats=NULL;
double *COLLpermstats=NULL;
double *CMATpermstats=NULL;

int bins2test=0;

void write_setid(int setid, string SetIDfile, int intervaleditor, float raref,int nwindowssinglebin, int **windowPositions,fstream &errorfile, fstream &logfile, INTERVALS *intervals,struct MAP *map, struct WINDOW *window, string intervalfile, int *nRareLimits){
  fstream SetID;
  SetID.open(SetIDfile.c_str(), ios::out);
  for(int m=0; m<nwindowssinglebin; m++){
    for(int j=0; j<window[m].n_at_level[nRareLimits[m]-1]; j++){
      int i=window[m].levelpos[nRareLimits[m]-1][j];
      SetID<<map[i].chr<<":"<<map[windowPositions[m][0]].pos<<"-"<<map[windowPositions[m][1]].pos;
      if(intervalfile!=" " && featurecol!=-9){
     	SetID <<"("<<intervals[m].feature<<")"<<flush;
      }
      SetID << "\t"<< map[i].rs<<"\n";
    }
  }
  SetID.close();
  if(intervaleditor){
    cout << "SetID file "<<SetIDfile<<" generated."<<endl;
    logfile << "SetID file "<<SetIDfile<<" generated."<<endl;
    exit(0);
  }
};



double chiinv(double df, double p);

//struct WINDOW *window=NULL;



double pValueCalc(double a, double x);

int **calloc2Dint(int m, int n);
float **calloc2Dfloat(int m, int n);

double **calloc2Ddouble(int m, int n);
void free2Ddouble(double **matrix, int m);

void updateCounts(struct COUNTS*, int, uint64_t**, uint64_t**, uint64_t**, bool X);
void aff_permute_clusters(struct PERSON*, struct CLUSTERS*, uint32_t, int, int*, int*, int*);

// int comparator
int compareint (const void * a, const void * b){
  return ( *(int*)a - *(int*)b );
}

// float comparator
int comparefloat(const void *a, const void  *b){
  const float *da = (const float *) a;
  const float *db = (const float *) b;
  return (*da > *db) - (*da < *db);
}

// double comparator
int comparedouble(const void *a, const void  *b){
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  return (*da > *db) - (*da < *db);
}


int israre(int i, double maf, double raref, struct MAP *map){
  if(map[i].analysis_in == 1 && ((raref-maf)>EPS || fabs(maf-raref)<EPS) && maf>EPS){
    return 1;
  }
  else{
    return 0;
  }
}


int checkCovCathegories(int *cov, int nlinestfam,struct PERSON *person, int nCovCathegories){
  cout << "Attempting to use the (first) covariate for the CMAT test with covariates"<<endl;
  int dummycov[nlinestfam];
  int sdummycov[nlinestfam];
  nCovCathegories=1;

  for(int i=0; i<nlinestfam; i++){
    dummycov[i]=int(person[i].cov[0]);
    sdummycov[i]=int(person[i].cov[0]);
  }
  qsort(sdummycov,nlinestfam,sizeof(int),compareint);
  for(int i=1; i<nlinestfam; i++){
    if(sdummycov[i-1]!=sdummycov[i]){
      nCovCathegories++;
    }
  }
  cout <<nCovCathegories<<" categories found\n"<<endl;
  if(sdummycov[0]==0){
    cout<<"Among them missings!"<<endl;
    if( sdummycov[nlinestfam-1]!=nCovCathegories-1){
      cout <<"Categoric covariates have to be numbered 1 through "<<nCovCathegories<<" for CMAT test!\n"<<endl;
      exit(0);
    }
  }
  else if( sdummycov[nlinestfam-1]!=nCovCathegories ){
    cout <<"Categoric covariates have to be numbered 1 through "<<nCovCathegories<<" for CMAT test!\n"<<endl;
    exit(0);
  }
  return nCovCathegories;
}



double wilson_lower(int ncount, int pcount){
  //  float z = 3.7190165; // for alpha=10^-4 one-sided
  double z = 1.96;
  double phat = pcount / double(ncount);
  return max((double)0,((phat + z*z/(2*ncount) - z * sqrt((phat*(1-phat)+z*z/(4*ncount))/ncount))/(1+z*z/ncount)));
}
double wilson_upper(int ncount, int pcount){
  double z = 1.96;
  //  float z = 3.7190165;
  double phat = pcount / double(ncount);
  return min((double)1,((phat + z*z/(2*ncount) + z * sqrt((phat*(1-phat)+z*z/(4*ncount))/ncount))/(1+z*z/ncount)));
}

double *get_weights(int nlinestped, struct WEIGHTS weights, struct MAP *map, double *weightvec){
  int i=0;

  for(i=0; i<nlinestped; i++){
    if(weights.mode==1 && map[i].maf!=0){
      weightvec[i]=(double)1/sqrt(map[i].mafr*(1-map[i].mafr));
    }
    else if(weights.mode==2 && map[i].maf!=0){
      weightvec[i]=pow(map[i].mafr,weights.betapar1-1)*pow(1-map[i].mafr,weights.betapar2-1);
    }
    else if(weights.mode==3 && map[i].maf!=0){
      weightvec[i]=(double)1/(1+exp((-weights.logpar1+map[i].mafr)*weights.logpar2));
    }
    else if(weights.mode==0){
      weightvec[i]=1;
    }
  }

  return weightvec;
}

/// Fisher test
double calcFISHERp(double stat, double raref, int start, int end, struct COUNTS *counts, int thread, struct MAP *map, struct WINDOW *window, int l, int m){
  int i=0;
  //  int j=0;
  double p=0;

  //  for(i=start; i<=end; i++){
  //    if( ((map[i].mafr-raref)<EPS || fabs(map[i].mafr-raref)<EPS) && map[i].mafr>EPS && map[i].analysis_in == 1){
  //	  j++;
  //	}
  //    }
  p=pValueCalc((double)window[l].n_at_level[m], stat/2);
  return p;
}


double calcFISHER(int *cov, int sexcov, struct PERSON *person, int ncases, int ncontrols, int nrest, int thread, int nlinestfam, int nlinestped, int start, int end, struct COUNTS *counts,  double raref, double inflationfactor, bool fisherCorrection, struct MAP *map, double *weightvec, struct WINDOW *window, int l, int m){
  int i=0;
  double sum=0;
  double p=1;
  double stat=0;
  //  for(i=start; i<=end; i++){
  //      if(map[i].analysis_in && ((map[i].mafr-raref)<EPS || fabs(map[i].mafr-raref)<EPS) && map[i].mafr>EPS  && inflationfactor>EPS){
  for(int j=0; j<window[l].n_at_level[m]; j++){
    i=window[l].levelpos[m][j];

    p=(double)counts[i].pSingle;
    if(fisherCorrection){
      stat=chiinv(1,p);
      stat=stat/inflationfactor;
      p=pValueCalc(0.5, stat/2);
    }
    sum=sum+weightvec[i]*log(p);
    //	}
  }
  return -2*sum;
}


void find_distinct_vb_vt(struct PERSON *person,struct MAP *map, int* PPLMap, struct PPLLOCATION* PplLocations, uint64_t*** BinSNPs,  uint64_t** BinSNPsCCFlags, int nwords,  int* SNPMapInverse, struct WINDOW *window, int *nRareLimits, int **rareLimitsNCT, int **rareLimitsNCTInverse, int nwindows, int minIndInBin, int maxIndInBin,  int *nvbstartvt, int **vbstartvt, int nlinestfam, int nlinestped, int *nCarriers, int **nvblevel, int ***nchunks,  int ****chunkpos, int ****chunklen,  int optimalrare, int ****ndummyatlevel, int *****dummypos, int *****dummylevel, int ****dummyend, int ***ndummyends, int NCT, uint64_t ***BinCarriers, int ***nchunkcluster, int *****dummycluster){

  int equivstartsnps=0;
  //  nvbstart=0;
  int startcounter=0;


  int shift=0; // shift is the number of variants on previous chromosomes
  for(int l=0; l<nwindows; l++){
    nvbstartvt[l]=1;
    vbstartvt[l]=(int *)realloc(vbstartvt[l],(nvbstartvt[l])*sizeof(int));
    if(!vbstartvt[l]) die("Memory allocation error in vbstartvt[l]!");
    vbstartvt[l][0]=window[l].levelpos[nRareLimits[l]-1][0]; // index of first variant
    int i1=window[l].levelpos[nRareLimits[l]-1][0];
    uint64_t* dummy1 = new uint64_t[nwords]();
    if(!dummy1) die("Memory allocation error in dummy1!");
    int iMod = SNPMapInverse[i1];
    if (iMod==-1) continue;
    for (int p=0; p<nwords; p++){
      dummy1[p] |= (BinSNPs[iMod][p][1] | BinSNPs[iMod][p][2]);
    }
    int m1=0;
    //    cout<<"vbstartvt["<<l<<"]["<<nvbstartvt[l]-1<<"] "<<vbstartvt[l][nvbstartvt[l]-1]<<" is pos of snp "<<nvbstartvt[l]-1<<" "<<map[vbstartvt[l][nvbstartvt[l]-1]].rs<<" "<<map[vbstartvt[l][nvbstartvt[l]-1]].pos<<" "<<" nCarriers[vbstartvt[i][vbstartvt[l]["<<nvbstartvt[l]-1<<"]] "<<nCarriers[vbstartvt[l][nvbstartvt[l]-1]]<<endl;
    while(m1<window[l].n_at_level[nRareLimits[l]-1]-1){
      int m2=m1+1;
      while(m2<window[l].n_at_level[nRareLimits[l]-1]){
	int nIndiff=0; // flag for distinctivness of next variant
	int i2=window[l].levelpos[nRareLimits[l]-1][m2];
	uint64_t* dummy2 = new uint64_t[nwords]();
       if(!dummy2) die("Memory allocation error in dummy2!");
	int iMod2 = SNPMapInverse[i2];
	if (iMod2==-1) continue;
	for (int p=0; p<nwords; p++){
	  dummy2[p] |= (BinSNPs[iMod2][p][1] | BinSNPs[iMod2][p][2]);
	  if(!nIndiff && dummy1[p]!=dummy2[p]){
	    nIndiff=1;
	  }
	  dummy1[p] = dummy2[p];
	}
	m1=m2;
	if(nIndiff){
	  nvbstartvt[l]++;
	  vbstartvt[l]=(int *)realloc(vbstartvt[l],(nvbstartvt[l])*sizeof(int));
	  if(!vbstartvt[l]) die("Memory allocation error in vbstartvt[l]!");
	  vbstartvt[l][nvbstartvt[l]-1]=window[l].levelpos[nRareLimits[l]-1][m2];
	  //	  cout<<"vbstartvt["<<l<<"]["<<nvbstartvt[l]-1<<"] "<<vbstartvt[l][nvbstartvt[l]-1]<<" is pos of snp "<<nvbstartvt[l]-1<<" "<<map[vbstartvt[l][nvbstartvt[l]-1]].rs<<" "<<map[vbstartvt[l][nvbstartvt[l]-1]].pos<<" nCarriers[vbstartvt[i][vbstartvt[l]["<<nvbstartvt[l]-1<<"]] "<<nCarriers[vbstartvt[l][nvbstartvt[l]-1]]<<endl;
	}
	m2++;
	delete[] dummy2;
      }
    }
    shift+=window[l].n_at_level[nRareLimits[l]-1];
    delete[] dummy1;
  }


  // allocate memory,  compute BinCarriers
  // all arrays ending with *atlevel (except for ndummyatlevel) are temporary and replaced later with arrays without "atlevel"
  for(int l=0; l<nwindows; l++){
    nvblevel[l] = new int[nvbstartvt[l]]();
    if(!nvblevel[l]) die("Memory allocation error in nvblevel!");
    chunkpos[l] = new int**[nvbstartvt[l]];
    if(!chunkpos[l]) die("Memory allocation error in chunkpos!");
    chunklen[l] = new int**[nvbstartvt[l]];
    if(!chunklen[l]) die("Memory allocation error in chunklen!");
    nchunks[l] = new int*[nvbstartvt[l]];
    if(!nchunks[l]) die("Memory allocation error in nchunk!");
    // dummy is the actual binary bin vector. this is used for reducing the number of positions kept in memory.
    ndummyatlevel[l]=new int**[nvbstartvt[l]]; // ndummylevel is the number of dummys that are at the beginning of bins at higher levels
    if(!ndummyatlevel[l]) die("Memory allocation error in ndummyatlevel!");
    dummylevel[l]=new int***[nvbstartvt[l]];
    if(!dummylevel[l]) die("Memory allocation error in dummylevel!");
    dummypos[l]=new int***[nvbstartvt[l]];
    if(!dummypos[l]) die("Memory allocation error in dummypos!");
    dummycluster[l]=new int***[nvbstartvt[l]];
    if(!dummycluster[l]) die("Memory allocation error in dummycluster!");
    dummyend[l]=new int**[nvbstartvt[l]];
    if(!dummyend[l]) die("Memory allocation error in dummyend!");
    ndummyends[l]=new int*[nvbstartvt[l]];
    if(!ndummyends[l]) die("Memory allocation error in ndummyends!");
    nchunkcluster[l]=new int*[nvbstartvt[l]];
    if(!nchunkcluster[l]) die("Memory allocation error in nchunkcluster!");

    BinCarriers[l] = new uint64_t*[nvbstartvt[l]]; // BinCarriers[window][distinct variant][word] = binary vectors of carriers of distinct variants
    if(!BinCarriers[l]) die("Memory allocation error in BinCarriers!");
    for(int m=0; m<nvbstartvt[l]; m++){
      BinCarriers[l][m] = new uint64_t[nwords]();
      if(!BinCarriers[l][m]) die("Memory allocation error in BinCarriers!");
      int i2=vbstartvt[l][m];
      int iMod = SNPMapInverse[i2];
      if (iMod==-1) continue;
      for (int p=0; p<nwords; p++){
	BinCarriers[l][m][p] = (BinSNPs[iMod][p][1] | BinSNPs[iMod][p][2]);
      }
    }
  }

  // find chunkpos, chunklen: chunks are contigous sequences of indices, each chunk has a position and a length
  for(int l=0; l<nwindows; l++){
    for(int m1=0; m1<nvbstartvt[l]; m1++){
      int *ndummyendsatlevel=new int[nRareLimits[l]]();
      int **dummyendatlevel=new int*[nRareLimits[l]];
      int baselevel=0;
      if(optimalrare){
	baselevel=rareLimitsNCTInverse[l][nCarriers[vbstartvt[l][m1]]-1]; // baselevel is relative to chromosome. atlevelInverse[baselevel] relative to startposition is always 0
      }

      int *atlevel = new int[nRareLimits[l]](); // atlevel[level relative to startpos] = relative level, chromosome-wide
      if(!atlevel) die("Memory allocation error in atlevel!");
      int *atlevelInverse = new int[nRareLimits[l]]; // atlevelInverse[relative level, chromosome-wide] = level relative to startpos. in the beginning, set to -1, except atlevelInverse[baselevel]=0
      if(!atlevelInverse) die("Memory allocation error in atlevelInverse!");
      int *atlevelInverseStatus = new int[nRareLimits[l]](); // encoding: 0=never seen; 1=has been seen and active; 2=has been seen and inactive; 3=finished
      if(!atlevelInverseStatus) die("Memory allocation error in atlevelInverseStatus!");
      int *Ind1=new int[nRareLimits[l]]();
      if(!Ind1) die("Memory allocation error in Ind1!");
      int *Ind2=new int[nRareLimits[l]]();
      if(!Ind2) die("Memory allocation error in Ind2!");
      int *nRemaining=new int[nRareLimits[l]]();
      if(!nRemaining) die("Memory allocation error in nRemaining!");

      int *nchunksatlevel=new int[nlinestfam](); // temporary arrays
      if(!nchunksatlevel) die("Memory allocation error in nchunksatlevel!");
      int **chunkposatlevel=new int*[nRareLimits[l]];
      if(!chunkposatlevel) die("Memory allocation error in chunkposatlevel!");
      int **chunklenatlevel=new int*[nRareLimits[l]];
      if(!chunklenatlevel) die("Memory allocation error in chunklenatlevel!");
      // all variants with AND:
      uint64_t **dummy1 = new uint64_t*[nRareLimits[l]];
      if(!dummy1) die("Memory allocation error in dummy1!");
      // basevariant minus following variants, if dummy2=0 then B_{i,j}=B_{i+1,j}
      uint64_t **dummy2 = new uint64_t*[nRareLimits[l]];
      if(!dummy2) die("Memory allocation error in dummy2!");

      for(int i=0; i<nRareLimits[l]; i++){
	// dummyendatlevel[i] is the first position that is not in cluster anymore
	dummyendatlevel[i]=new int[nlinestfam]();
	if(!dummyendatlevel[i]) die("Memory allocation error in dummyendatlevel[i]!");
	chunkposatlevel[i]=new int[nlinestfam]();
	if(!chunkposatlevel[i]) die("Memory allocation error in chunkposatlevel[i]!");
	chunklenatlevel[i]=new int[nlinestfam]();
	if(!chunklenatlevel[i]) die("Memory allocation error in chunklenatlevel[i]!");
	dummy1[i]=new uint64_t[nwords]();
	if(!dummy1) die("Memory allocation error in dummy1!");
	dummy2[i]=new uint64_t[nwords]();
	if(!dummy2) die("Memory allocation error in dummy2!");
	atlevelInverse[i]=-1;
	atlevel[i]=-1;
      }

      nvblevel[l][m1]++;
      nchunksatlevel[0]++;
      chunkposatlevel[0][0]=m1;
      chunklenatlevel[0][0]=1;
      atlevel[0]=baselevel;
      if(optimalrare){
	atlevelInverse[baselevel]=0;
	atlevelInverseStatus[baselevel]=1;
      }
      else if(!optimalrare){
	atlevelInverse[0]=0;
      }

      for (int p=0; p<nwords; p++){
	dummy1[0][p] |= BinCarriers[l][m1][p];
	dummy2[0][p] = dummy1[0][p];
	Ind1[0]+=bitcount64(dummy1[0][p]);
      }
      for(int m2=m1+1; m2<nvbstartvt[l]; m2++){ // go through end positions
	int currentlevel=0;
	if(optimalrare){
	  currentlevel=rareLimitsNCTInverse[l][nCarriers[vbstartvt[l][m2]]-1]; // relative to chromsome
	}
	int newbin=0;
	// with VT, check if a new higher (>baselevel) level exists:
	if(currentlevel>baselevel){
	  int previouslevel=-2;

	  if(atlevelInverseStatus[currentlevel]==3){ // if level already satisfied row termination condition, skip variant
	    continue;
	  }
	  else if(atlevelInverseStatus[currentlevel]!=3){
	    if(atlevelInverseStatus[currentlevel]==1){ // level seen and active
	      previouslevel=atlevelInverse[currentlevel];
	    }
	    else if(atlevelInverseStatus[currentlevel]==0 || atlevelInverseStatus[currentlevel]==2){
	      for(int j=currentlevel-1; j>=baselevel; j--){ // find the largest known level < currentlevel:
		if(atlevelInverseStatus[j]==3){ // a level below satisfied row termination condition
		  for(int j1=currentlevel; j1>j; j1--){ // all levels are terminated already
		    atlevelInverseStatus[j1]=3;
		  }
		  for(int j1=currentlevel+1; j1<nRareLimits[l]; j1++){ // all levels are terminated already
		    atlevelInverseStatus[j1]=3;
		  }
		  previouslevel=-2;
		  break;
		}
		else if(atlevelInverseStatus[j]==1){ // only if seen and active
		  previouslevel=atlevelInverse[j];
		  break;
		}
	      }
	    }

	    if(previouslevel>=0){
	      uint64_t *dummy1_tmp=new uint64_t[nwords]();
	      if(!dummy1_tmp) die("Memory allocation error in dummy1_tmp!");
	      uint64_t *dummy2_tmp=new uint64_t[nwords]();
	      if(!dummy2_tmp) die("Memory allocation error in dummy2_tmp!");
	      int Ind1_tmp=Ind1[previouslevel];
	      int Ind2_tmp=0;
	      int nRemaining_tmp=0;
	      int i2=vbstartvt[l][m2];
	      int iMod2 = SNPMapInverse[i2];
	      // check if at currentlevel, B(m1,m2)==B(m1+1,m2) and  B(m1,m2)==B(m1,m2+1):
	      for (int p=0; p<nwords; p++){
		dummy1_tmp[p] = (dummy1[previouslevel][p] |  BinCarriers[l][m2][p]);
		dummy2_tmp[p] = (dummy2[previouslevel][p] & ~(BinCarriers[l][m2][p]));
		Ind2_tmp+=bitcount64(dummy1_tmp[p]);
		nRemaining_tmp+=bitcount64(dummy2_tmp[p]);
	      }
	      if(nRemaining_tmp==0 || Ind2_tmp==nlinestfam){ // this and higher levels are forbidden: row termination condition reached
		for(int i=currentlevel; i<nRareLimits[l]; i++){
		  atlevelInverseStatus[i]=3;
		}
		delete[] dummy1_tmp;
		delete[] dummy2_tmp;
		continue;
	      }
	      else if(Ind1_tmp==Ind2_tmp){
		delete[] dummy1_tmp;
		delete[] dummy2_tmp;
		continue;
	      }
	      // something new at this level
	      else if(Ind1_tmp!=Ind2_tmp){
		// check if there is a higher open level and it is equivalent to this bin

		newbin=1;
		if(atlevelInverseStatus[currentlevel]==0){ // genuinely new level
		  nvblevel[l][m1]++;
		  nchunksatlevel[nvblevel[l][m1]-1]++;
		  atlevelInverseStatus[currentlevel]=1;
		  atlevelInverse[currentlevel]=nvblevel[l][m1]-1;
		  chunkposatlevel[nvblevel[l][m1]-1][nchunksatlevel[nvblevel[l][m1]-1]-1]=m2;
		  chunklenatlevel[nvblevel[l][m1]-1][nchunksatlevel[nvblevel[l][m1]-1]-1]=1;
 		  atlevel[nvblevel[l][m1]-1]=currentlevel;
		}
		else if(atlevelInverseStatus[currentlevel]==1 || atlevelInverseStatus[currentlevel]==2){
		  if(atlevelInverseStatus[currentlevel]==2){ // if inactive, activate, because checked down already:Ind1!=Ind2, Ind1 is from previous level
		    atlevelInverseStatus[currentlevel]=1;
		  }
		  if(chunkposatlevel[atlevelInverse[currentlevel]][nchunksatlevel[atlevelInverse[currentlevel]]-1]+chunklenatlevel[atlevelInverse[currentlevel]][nchunksatlevel[atlevelInverse[currentlevel]]-1]==m2){
		    chunklenatlevel[atlevelInverse[currentlevel]][nchunksatlevel[atlevelInverse[currentlevel]]-1]++;
		  }
		  else if(chunkposatlevel[atlevelInverse[currentlevel]][nchunksatlevel[atlevelInverse[currentlevel]]-1]+chunklenatlevel[atlevelInverse[currentlevel]][nchunksatlevel[atlevelInverse[currentlevel]]-1]!=m2){
		    nchunksatlevel[atlevelInverse[currentlevel]]++;
		    chunkposatlevel[atlevelInverse[currentlevel]][nchunksatlevel[atlevelInverse[currentlevel]]-1]=m2;
		    chunklenatlevel[atlevelInverse[currentlevel]][nchunksatlevel[atlevelInverse[currentlevel]]-1]=1;
		  }
		}
		Ind1[atlevelInverse[currentlevel]]=Ind2_tmp;
		for(int p=0; p<nwords; p++){
		  dummy1[atlevelInverse[currentlevel]][p]=dummy1_tmp[p];
		  dummy2[atlevelInverse[currentlevel]][p]=dummy2_tmp[p];

		}
	      }
	      delete[] dummy1_tmp;
	      delete[] dummy2_tmp;
	    }
	  }
	}
	else if(currentlevel<=baselevel){
	  // case distinction: is this new at baselevel? atlevelInverse[basevel]=0:
	  Ind2[0]=0;
	  nRemaining[0]=0;
	  int i2=vbstartvt[l][m2];
	  int iMod2 = SNPMapInverse[i2];
	  //	  if (iMod2==-1) continue;

	  for (int p=0; p<nwords; p++){
	    dummy1[0][p] |= BinCarriers[l][m2][p];
	    dummy2[0][p] &= ~BinCarriers[l][m2][p];
	    Ind2[0]+=bitcount64(dummy1[0][p]);
	    nRemaining[0]+=bitcount64(dummy2[0][p]);
	  }
	  if(m1==0){
	    //	    cout<<"nRemaining[0]" <<nRemaining[0]<<endl;
	  }
	  if(nRemaining[0]==0 || Ind2[0]==nlinestfam){ // continue with next start position
	    break;
	  }
	  else if(Ind2[0]==Ind1[0]){ // continue with next end position
	    continue;
	  }
	  else if(Ind2[0]!=Ind1[0]){ // if something new at baselevel
	    newbin=1;
	    Ind1[0]=Ind2[0];
	    if(m2!=chunkposatlevel[0][nchunksatlevel[0]-1]+chunklenatlevel[0][nchunksatlevel[0]-1]){
	      nchunksatlevel[0]++;
	      chunkposatlevel[0][nchunksatlevel[0]-1]=m2;
	    }
	    chunklenatlevel[0][nchunksatlevel[0]-1]++;
	  }
	  currentlevel=baselevel;
	}

	// only if something new at currentlevel, check if at larger level, B(m1,m2)==B(m1+1,m2) and  B(m1,m2)==B(m1,m2+1):
	if(newbin==1){
	  int prevlevel=currentlevel; // if in level above something is found, check down at previouslevel
	  for(int m3=currentlevel+1; m3<nRareLimits[l]; m3++){
	    if(atlevelInverseStatus[m3]==3){
	      break;
	    }
	    else if(atlevelInverseStatus[m3]==1){
	      // sanity check
	      if(nCarriers[vbstartvt[l][chunkposatlevel[atlevelInverse[currentlevel]][0]]]>nCarriers[vbstartvt[l][chunkposatlevel[atlevelInverse[m3]][0]]] && atlevelInverseStatus[currentlevel]==1){
		cout<<nCarriers[vbstartvt[l][chunkposatlevel[atlevelInverse[currentlevel]][0]]]<<" "<<nCarriers[vbstartvt[l][chunkposatlevel[atlevelInverse[m3]][0]]]<<endl;
		cout<<"BUG in atlevelInverse["<<m3<<"] "<<atlevelInverse[m3]<<" "<<atlevelInverse[currentlevel]<<endl;
		exit(1);
	      }
	      //		      prevlevel=m3;
		      //		    atlevelInverseStatus[m3]==1
	      nRemaining[atlevelInverse[m3]]=0;
	      //		    Ind1[atlevelInverse[m3]]=Ind1[atlevelInverse[prevlevel]];
	      Ind2[atlevelInverse[m3]]=0;
	      int i2=vbstartvt[l][m2];
	      int iMod2 = SNPMapInverse[i2];
	      //		    if (iMod2==-1) continue;
	      for (int p=0; p<nwords; p++){
		dummy1[atlevelInverse[m3]][p] |= BinCarriers[l][m2][p];
		dummy2[atlevelInverse[m3]][p] &= ~BinCarriers[l][m2][p];
		Ind2[atlevelInverse[m3]]+=bitcount64(dummy1[atlevelInverse[m3]][p]);
		nRemaining[atlevelInverse[m3]]+=bitcount64(dummy2[atlevelInverse[m3]][p]);
	      }
	      //  B(m1,m2)==B(m1+1,m2): end of this and higher levels;
	      if(Ind2[atlevelInverse[m3]]==nlinestfam || nRemaining[atlevelInverse[m3]]==0){
		atlevelInverseStatus[m3]=3;
		for(int i=m3+1; i<nRareLimits[l]; i++){
		  //			if(atlevelInverseStatus[i]==1 || atlevelInverseStatus[i]==2){
		  atlevelInverseStatus[i]=3;
		  //			}
		}
		break;
	      }
	      else if(Ind2[atlevelInverse[m3]]!=Ind1[atlevelInverse[m3]]){ // existing open level
		// check down if differrent from prevous level
		if(Ind2[atlevelInverse[m3]]==Ind1[atlevelInverse[prevlevel]]){ // existing open level)
		  ndummyendsatlevel[atlevelInverse[m3]]++;
		  dummyendatlevel[atlevelInverse[m3]][ndummyendsatlevel[atlevelInverse[m3]]-1]=m2;
		  atlevelInverseStatus[m3]=2;
		}
		else if(Ind2[atlevelInverse[m3]]!=Ind1[atlevelInverse[prevlevel]]){ // existing open level
		  if(chunkposatlevel[atlevelInverse[m3]][nchunksatlevel[atlevelInverse[m3]]-1]+chunklenatlevel[atlevelInverse[m3]][nchunksatlevel[atlevelInverse[m3]]-1]!=m2){
		    nchunksatlevel[atlevelInverse[m3]]++;
		    chunkposatlevel[atlevelInverse[m3]][nchunksatlevel[atlevelInverse[m3]]-1]=m2;
		  }
		  chunklenatlevel[atlevelInverse[m3]][nchunksatlevel[atlevelInverse[m3]]-1]++;
		  Ind1[atlevelInverse[m3]]=Ind2[atlevelInverse[m3]];
		  prevlevel=m3;
		}
	      }
	    }
	  }
	} // new bin was found with m2 at lowest level
      } // row termination position for baselevel is reached

      // copy nchunks, chunkpos, chunklen
      nchunks[l][m1]=new int[nvblevel[l][m1]];
      if(!nchunks[l][m1]) die("Memory allocation error in nchunks[l][m1]!");
      memcpy(nchunks[l][m1], nchunksatlevel, nvblevel[l][m1]*sizeof(int));
      delete[] nchunksatlevel;
      chunkpos[l][m1] = new int*[nvblevel[l][m1]];
      if(!chunkpos[l][m1]) die("Memory allocation error in chunkpos[l][m1]!");
      chunklen[l][m1] = new int*[nvblevel[l][m1]];
      if(!chunklen[l][m1]) die("Memory allocation error in chunklen[l][m1]!");
      for(int m3=0; m3<nvblevel[l][m1]; m3++){
	chunkpos[l][m1][m3] = new int[nchunks[l][m1][m3]];
	if(!chunkpos[l][m1][m3]) die("Memory allocation error in chunkpos[l][m1][m3]!");
	chunklen[l][m1][m3] = new int[nchunks[l][m1][m3]];
	if(!chunklen[l][m1][m3]) die("Memory allocation error in chunklen[l][m1][m3]!");
	memcpy(chunkpos[l][m1][m3], chunkposatlevel[m3], nchunks[l][m1][m3]*sizeof(int));
	memcpy(chunklen[l][m1][m3], chunklenatlevel[m3], nchunks[l][m1][m3]*sizeof(int));
      }
      for(int m=0; m<nRareLimits[l]; m++){
	delete[] chunkposatlevel[m];
	delete[] chunklenatlevel[m];
      }
      delete[] chunkposatlevel;
      delete[] chunklenatlevel;

      //consistency check:
      for(int m2=0; m2<nvblevel[l][m1]; m2++){
	for(int m3=0; m3<nchunks[l][m1][m2]; m3++){
	  //  cout<<"nCarriers[vbstartvt["<<l<<"][chunkpos["<<l<<"]["<<m1<<"]["<<0<<"]["<<0<<"]]] "<<nCarriers[vbstartvt[l][chunkpos[l][m1][0][0]]]<<endl;
	  if(chunklen[l][m1][m2][m3]==0){
	    cout<<"BUG in chunklen!"<<endl;
	    exit(1);
	  }
	  for(int m4=0; m4<chunklen[l][m1][m2][m3]; m4++){
	    //    cout<<"nCarriers[vbstartvt["<<l<<"][chunkpos["<<l<<"]["<<m1<<"]["<<m2<<"]["<<m3<<"]+"<<m4<<"]] "<<nCarriers[vbstartvt[l][chunkpos[l][m1][m2][m3]+m4]]<<endl;
	    if(optimalrare && nCarriers[vbstartvt[l][chunkpos[l][m1][m2][m3]+m4]]>nCarriers[vbstartvt[l][chunkpos[l][m1][m2][0]]]){
	      cout<<nCarriers[vbstartvt[l][chunkpos[l][m1][m2][0]]]<<" "<<nCarriers[vbstartvt[l][chunkpos[l][m1][m2][m3]+m4]]<<" "<<chunklen[l][m1][m2][m3]<<endl;
	      cout<<l<<" "<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<map[vbstartvt[l][m1]].rs<<" "<<map[vbstartvt[l][chunkpos[l][m1][m2][m3]+m4]].rs<<endl;
	      cout<<"BUG in chunkpos, chunklen!"<<endl;
	      exit(1);
	    }
	  }
	}
      }

      if(nvblevel[l][m1]==0){
	cout<<"BUG in nvblevel[l][m1]"<<endl;
	exit(1);
      }

      ndummyends[l][m1]=new int[nvblevel[l][m1]]();
      if(!ndummyends[l][m1]) die("Memory allocation error in ndummyends[l][m1]!");
      nchunkcluster[l][m1]=new int[nvblevel[l][m1]]();
      if(!nchunkcluster[l][m1]) die("Memory allocation error in nchunkcluster[l][m1]!");
      dummyend[l][m1]=new int*[nvblevel[l][m1]];
      for(int m3=0; m3<nvblevel[l][m1]; m3++){
	ndummyends[l][m1][m3]=ndummyendsatlevel[m3];
	//	cout<<"ndummyendsatlevel["<<m3<<"] "<<ndummyendsatlevel[m3]<<endl;
	if(ndummyendsatlevel[m3]>0){
	  dummyend[l][m1][m3]=new int[ndummyendsatlevel[m3]]();
	  if(!dummyend[l][m1][m3]) die("Memory allocation error in ndummyends[l][m1][m3]!");
	  memcpy(dummyend[l][m1][m3],dummyendatlevel[m3], ndummyendsatlevel[m3]*sizeof(int));
	}
	else if(ndummyendsatlevel[m3]==0){
	  dummyend[l][m1][m3]=new int[1]();
	  if(!dummyend[l][m1][m3]) die("Memory allocation error in ndummyend[l][m1][m3]!");
	}
      }
      for(int i=0; i<nRareLimits[l]; i++){
	delete[] dummyendatlevel[i];
      }
      delete[] dummyendatlevel;
      delete[] ndummyendsatlevel;

      ndummyatlevel[l][m1]=new int*[nvblevel[l][m1]];
      if(!ndummyatlevel[l][m1]) die("Memory allocation error in ndummyatlevel[l][m1]!");
      dummypos[l][m1]= new int**[nvblevel[l][m1]];
      if(!dummypos[l][m1]) die("Memory allocation error in dummypos[l][m1]!");
      dummycluster[l][m1]= new int**[nvblevel[l][m1]];
      if(!dummycluster[l][m1]) die("Memory allocation error in dummycluster[l][m1]!");
      dummylevel[l][m1]= new int**[nvblevel[l][m1]];
      if(!dummylevel[l][m1]) die("Memory allocation error in dummylevel[l][m1]!");

      // count nchunkcluster: ndummyends[l][m1][m3]+1
      // dummyends: if ndummyends>0, check if chunks extend beyond last end. if no, reduce nchunkcluster
      for(int m3=0; m3<nvblevel[l][m1]; m3++){
	nchunkcluster[l][m1][m3]=ndummyends[l][m1][m3]+1;
	if(ndummyends[l][m1][m3]>0){
	  if(chunkpos[l][m1][m3][nchunks[l][m1][m3]-1]+chunklen[l][m1][m3][nchunks[l][m1][m3]-1]-1<=dummyend[l][m1][m3][ndummyends[l][m1][m3]-1]){
	    nchunkcluster[l][m1][m3]--;
	  }
	  //  else if(chunkpos[l][m1][m3][nchunks[l][m1][m3]-1]+chunklen[l][m1][m3][nchunks[l][m1][m3]-1]-1==dummyend[l][m1][m3][ndummyends[l][m1][m3]-1]){
	  //    cout<<"BUG in ndummyends!"<<endl;
	  //    exit(1);
	  //  }
	}
	if(nchunkcluster[l][m1][m3]==0){
	  cout<<"NCHUNKCLUSTER["<<l<<"]["<<m1<<"]["<<m3<<"] "<<nchunkcluster[l][m1][m3]<<endl;
	}
	//cout<<"nchunkcluster["<<l<<"]["<<m1<<"]["<<m3<<"] "<<nchunkcluster[l][m1][m3]<<endl;
	ndummyatlevel[l][m1][m3]=new int[nchunkcluster[l][m1][m3]]();
	if(!ndummyatlevel[l][m1][m3]) die("Memory allocation error in ndummyatlevel[l][m1][m3]!");
	dummypos[l][m1][m3]= new int*[nchunkcluster[l][m1][m3]];
	if(!dummypos[l][m1][m3]) die("Memory allocation error in dummypos[l][m1][m3]!");
	dummycluster[l][m1][m3]= new int*[nchunkcluster[l][m1][m3]];
	if(!dummycluster[l][m1][m3]) die("Memory allocation error in dummycluster[l][m1][m3]!");
	dummylevel[l][m1][m3]= new int*[nchunkcluster[l][m1][m3]];
	if(!dummylevel[l][m1][m3]) die("Memory allocation error in dummylevel[l][m1][m3]!");
      }

      int **startposatlevel=new int*[nvblevel[l][m1]];
      if(!startposatlevel) die("Memory allocation error in startposatlevel!");
      for(int m3=0; m3<nvblevel[l][m1]; m3++){
	//	cout<<"alloc "<<m1<<" "<<nvblevel[l][m1]<<" "<<m3<<" "<<nchunkcluster[l][m1][m3]<<endl;
	startposatlevel[m3]=new int[nchunkcluster[l][m1][m3]]();
	if(!startposatlevel[m3]) die("Memory allocation error in startposatlevel[m3]!");
      }

      // find startposatlevel[level][cluster]. startpos[level][>0] is first position after ndummyends
      /*if(m1==1966){
	cout<<"nvblevel["<<l<<"]["<<m1<<"] "<<nvblevel[l][m1]<<endl;
	} */
      for(int k=0; k<nvblevel[l][m1]; k++){
	//atlevelCarriers[k]=nCarriers[vbstartvt[l][chunkpos[l][m1][k][0]]];
	startposatlevel[k][0]=chunkpos[l][m1][k][0];
	//cout<<"chunkpos["<<l<<"]["<<m1<<"]["<<k<<"]["<<0<<"] "<<chunkpos[l][m1][k][0]<<"!!!"<<endl;
	//      cout<<"startposatlevel["<<k<<"]["<<0<<"] "<<startposatlevel[k][0]<<endl;

	if(ndummyends[l][m1][k]>0){
	  int foundstart=0;
	  for(int m4=0; m4<nchunks[l][m1][k]; m4++){
	    if(m1==2 /*m2==158 */){
	      //	     	      cout<<"chunkpos["<<l<<"]["<<m1<<"]["<<k<<"]["<<m4<<"]"<<" "<<"dummyend["<<l<<"]["<<m1<<"]["<<k<<"]["<<foundstart<<"] "<<  chunkpos[l][m1][k][m4]<<" "<<dummyend[l][m1][k][foundstart]<<endl;
		      //		      cout <<dummyend[l][m1][k][foundstart]<<" chunklen[l][m1][k][m4] "<< chunklen[l][m1][k][m4] <<endl;
	      }
	    if(chunkpos[l][m1][k][m4]>dummyend[l][m1][k][foundstart]){
	      startposatlevel[k][foundstart+1]=chunkpos[l][m1][k][m4];
	      if(m1==2){
		//		cout <<"m4 "<<m4<<" nchunks["<<l<<"]["<<m1<<"]["<<k<<"] "<< nchunks[l][m1][k]<<endl;
		//		cout<<"dummyend["<<l<<"]["<<m1<<"]["<<k<<"]["<<foundstart<<"] "<<dummyend[l][m1][k][foundstart]<<endl;
		//		cout<<"startposatlevel["<<k<<"]["<<foundstart+1<<"] "<<startposatlevel[k][foundstart+1]<<endl;
		//		cout<<"foundstart "<<foundstart<<endl;
		//		cout<<foundstart<<" "<<ndummyends[l][m1][k]<<" "<< ndummyends[l][m1][k]<<" "<<nchunkcluster[l][m1][k]<<" "<<foundstart<<" "<<ndummyends[l][m1][k]+1<<" "<< ndummyends[l][m1][k]<<" "<<nchunkcluster[l][m1][k]-1<<endl;
		}
	      foundstart++;
	      // PROBABLY DRAGONS here, does foundstart count across clusters?
	      //      cout<<foundstart<<" "<<ndummyends[l][m1][k]<<" "<< ndummyends[l][m1][k]<<" "<<nchunkcluster[l][m1][k]<<" "<<foundstart<<" "<<ndummyends[l][m1][k]+1<<" "<< ndummyends[l][m1][k]<<" "<<nchunkcluster[l][m1][k]-1<<endl;
	      if((foundstart+1==ndummyends[l][m1][k] && ndummyends[l][m1][k]==nchunkcluster[l][m1][k]) || (foundstart+1==ndummyends[l][m1][k]+1 && ndummyends[l][m1][k]==nchunkcluster[l][m1][k]-1)){
		break;
	      }
	    }
	  }
	}
      }


      for(int k=0; k<nvblevel[l][m1]; k++){
	for(int k1=0; k1<ndummyends[l][m1][k1]; k1++){
	  //	  cout<<"dummyendatlevel["<<k<<"]["<<k1<<"] "<<dummyendatlevel[k][k1]<<endl;
	}
      }

      // find ndummyatlevel: number of dummys that shall be copied from level
      // find positiions and levels for recycling of dummies for association testing
      int ***dummyposatlevel= new int**[nvblevel[l][m1]];
      if(!dummyposatlevel) die("Memory allocation error in dummyposatlevel!");
      int ***dummylevelatlevel= new int**[nvblevel[l][m1]];
      if(!dummylevelatlevel) die("Memory allocation error in dummylevelatlevel!");
      int ***dummyclusteratlevel= new int**[nvblevel[l][m1]];
      if(!dummyclusteratlevel) die("Memory allocation error in dummyclusteratlevel!");
      for(int k=0; k<nvblevel[l][m1]; k++){
	dummyposatlevel[k]= new int*[ndummyends[l][m1][k]+1];
	if(!dummyposatlevel[k]) die("Memory allocation error in dummyposatlevel[k]!");
	dummyclusteratlevel[k]= new int*[ndummyends[l][m1][k]+1];
	if(!dummyclusteratlevel[k]) die("Memory allocation error in dummyclusteratlevel[k]!");
	dummylevelatlevel[k]= new int*[ndummyends[l][m1][k]+1];
	if(!dummylevelatlevel[k]) die("Memory allocation error in dummylevelatlevel[k]!");
	for(int k1=0; k1<ndummyends[l][m1][k]+1; k1++){
	  dummyposatlevel[k][k1]= new int[nlinestped]();
	  if(!dummyposatlevel[k][k1]) die("Memory allocation error in dummyposatlevel[k][k1]!");
	  dummyclusteratlevel[k][k1]= new int[nlinestped]();
	  if(!dummyclusteratlevel[k][k1]) die("Memory allocation error in dummyclusteratlevel[k][k1]!");
	  //	  cout<<"pre "<<m1<<" "<<k<<" "<<k1<<" "<<nvblevel[l][m1]<<endl;
	  dummylevelatlevel[k][k1]= new int[nlinestped]();
	  if(!dummylevelatlevel[k][k1]) die("Memory allocation error in dummylevelatlevel[k][k1]!");
	}
	if(k==0 && ndummyends[l][m1][k]>0){ // sanity check
	  cout<<"BUG in ndummyends[l][m1][0]!"<<endl;
	}
      }

      // at lower levels, find postmp, levtmp for dummys of higher levels
      // this looks a little complicated, but it works
      //      cout<<"nvblevel["<<l<<"]["<<m1<<"] "<<nvblevel[l][m1]<<endl;
      for(int k=1; k<nvblevel[l][m1]; k++){ // begin at level above baselevel
	int *postmp=new int[nchunkcluster[l][m1][k]];
	if(!postmp) die("Memory allocation error in postmp!");
	int *chunktmp=new int[nchunkcluster[l][m1][k]]();
	if(!chunktmp) die("Memory allocation error in chunktmp!");
	int *levtmp=new int[nchunkcluster[l][m1][k]]();
	if(!levtmp) die("Memory allocation error in levtmp!");
	int clustertmp=0;
	if(nchunkcluster[l][m1][k]==0){
	  //	  cout<<"nchunkcluster["<<l<<"]["<<m1<<"]["<<k<<"] "<<nchunkcluster[l][m1][k]<<" startposatlevel["<<k<<"]["<<0<<"] "<<startposatlevel[k][0]<<endl;
	  //	  cout<<"dummyendatlevel["<<l<<"]["<<m1<<"]["<<k<<"] "<<ndummyends[l][m1][k]<<endl;
	  exit(1);
	}
	for(int k2=0; k2<nchunkcluster[l][m1][k]; k2++){

	  /*  if(m1==1966 ){

	      cout<<"nchunkcluster["<<l<<"]["<<m1<<"]["<<k<<"] "<<nchunkcluster[l][m1][k]<<" startposatlevel["<<k<<"]["<<k2<<"] "<<startposatlevel[k][k2]<<endl;
	      } */
	  postmp[k2]=startposatlevel[k][k2];
	  // note that there is a break statement at the end of the loop. position has to be found before break is reached
	  for(int lev=atlevel[k]-1; lev>=0; lev--){ // check lower levels, in decending order
	    int found=0;
	    int k3=atlevelInverse[lev]; // convert back to relative level
	    if(k3<0){ // skip non-existing level
	      continue;
	    }
	    if(nchunkcluster[l][m1][k3]==0){
	      cout<<"nchunkcluster["<<l<<"]["<<m1<<"]["<<k3<<"] "<<nchunkcluster[l][m1][k3]<<" startposatlevel["<<k<<"]["<<0<<"] "<<startposatlevel[k][0]<<endl;
	      cout<<"ndummyends["<<l<<"]["<<m1<<"]["<<k3<<"] "<<ndummyends[l][m1][k3]<<" lev "<<lev<<endl;
	      exit(1);
	    }
	    int inbetween=0;
	    levtmp[k2]=k3;
	    if(startposatlevel[k3][0]>startposatlevel[k][k2]){ // the lower level does not exist yet for this end position
	      continue;
	    }
	    // here, we check the candidate lower level that has been seen and whose first start position is smaller or equal than startposition at level k
	    /*    if(m1==1966){
		  cout<<"atlevelInverse["<<lev<<"] m1 "<<atlevelInverse[lev]<<" "<<m1<<endl;
		  } */
	    if(ndummyends[l][m1][k3] < nchunkcluster[l][m1][k3]){  // if level is open: there are variants after the last dummyend
	      // check if at previous level, startpos[k] falls between or after dummyend[k3] and startpos[k3]:
	      if(ndummyends[l][m1][k3]==0){
		clustertmp=0;
		found=1;
		//  break;
	      }
	      else if(ndummyends[l][m1][k3]>0){
		for(int k4=0; k4<nchunkcluster[l][m1][k3]-1; k4++){
		  //    cout<<"startposatlevel["<<k<<"]["<<k2<<"] dummyend["<<l<<"]["<<m1<<"]["<<k3<<"]["<<k4<<"] startposatlevel["<<k<<"]["<<k2<<"] startposatlevel["<<k3<<"]["<<k4+1<<"]  "<<startposatlevel[k][k2]<<" "<<dummyend[l][m1][k3][k4]<<" " <<startposatlevel[k][k2]<<" "<<startposatlevel[k3][k4+1]<<endl;
		  if(startposatlevel[k][k2]>dummyend[l][m1][k3][k4] && startposatlevel[k][k2]<startposatlevel[k3][k4+1]){
		    inbetween=1;
		    break;
		  }
		  else if((startposatlevel[k][k2]<=dummyend[l][m1][k3][k4] && startposatlevel[k][k2]>=startposatlevel[k3][k4]) || startposatlevel[k][k2] >= startposatlevel[k3][nchunkcluster[l][m1][k3]-1]){
		    clustertmp=k4;
		    found=1;
		    break;
		  }
		}
	      }
	      if(inbetween==1){
		continue;
	      }
	    }
	    else if(ndummyends[l][m1][k3] == nchunkcluster[l][m1][k3]){
	      if(dummyend[l][m1][k3][ndummyends[l][m1][k3]-1]<startposatlevel[k][k2]){ // if lower level already closed:
		continue;
	      }
	      for(int k4=0; k4<nchunkcluster[l][m1][k3]; k4++){
		clustertmp=k4;
		// if found in one of the chunks:
		if(startposatlevel[k][k2]<=dummyend[l][m1][k3][k4] && startposatlevel[k][k2]>=startposatlevel[k3][k4]){
		  //      for(int k5=0; k5<nchunks[l][m1][k3]; k5++){
		  // if startpos definitely not in this chunk, skip chunk:
		  //if(k5<nchunks[l][m1][k3]-1 && startposatlevel[k][k2]>chunkpos[l][m1][k3][k5+1]){
		  //  continue;
		  //}
		  //for(int k6=0; k6<chunklen[l][m1][k3][k5]; k6++){
		  //  if(chunkpos[l][m1][k3][k5]+k6 < startposatlevel[k][k2]){
		  //    postmp[k2]=chunkpos[l][m1][k3][k5]+k6;
		  //    chunktmp[k2]=k5;
		  //  }
		  //  else if(chunkpos[l][m1][k3][k5]+k6 > startposatlevel[k][k2]){
		  //    break;
		  //  }
		  //  else if(chunkpos[l][m1][k3][k5]+k6 == startposatlevel[k][k2]){
		  //    cout<<" BUG in chunkpos["<<l<<"]["<<m1<<"]["<<k3<<"]["<<k5<<"]+"<<k6<<" startposatlevel["<<k<<"]["<<k2<<"] "<<startposatlevel[k][k2]<<endl;
		  //    exit(1);
		  //  }
		  //}
		  //      }
		  found=1;
		  break;
		}
		else if(startposatlevel[k][k2]>dummyend[l][m1][k3][k4] && startposatlevel[k][k2]<startposatlevel[k3][k4+1]){ // note: k4+1 never out of bounds!
		  inbetween=1;
		  break;
		}
	      }
	      if(inbetween==1){
		continue;
	      }
	    }


	    if(found==1){
	      //      break;
	    }
	    else if(inbetween==1){
	      continue;
	    }


	    if(found==0){
	      cout<<"found=0"<<endl;
	      exit(1);
	    }
	    //    cout<<levtmp[k2]<<" levtmp["<<k2<<"] "<<endl;
	    if(m1==1966 ){
	      //      cout<<"nchunks["<<l<<"]["<<m1<<"]["<<k3<<"] "<<nchunks[l][m1][k3]<<endl;
	    }
	    for(int k4=0; k4<nchunks[l][m1][k3]; k4++){
	      int found2=0;
	      // if startpos definitely not in this chunk, skip chunk:
	      if((m1==1966 ) && levtmp[k2]==0){
		//cout<<"startposatlevel["<<k<<"]["<<k2<<"] "<<startposatlevel[k][k2]<<endl;
	      }

	      if(k4<nchunks[l][m1][k3]-1 && startposatlevel[k][k2]>=chunkpos[l][m1][k3][k4+1]){
		continue;
	      }
	      if(startposatlevel[k][k2]<chunkpos[l][m1][k3][k4]){
		break;
		//exit(1);
	      }
	      for(int k5=0; k5<chunklen[l][m1][k3][k4]; k5++){

		/*if((m1==1966 ) ){
		  cout<<"chunkpos["<<l<<"]["<<m1<<"]["<<k3<<"]["<<k4<<"]+"<<k5<<" "<< chunkpos[l][m1][k3][k4]+k5<<" startposatlevel["<<k<<"]["<<k2<<"] "<<startposatlevel[k][k2]<<endl;
		  } */
		if(chunkpos[l][m1][k3][k4]+k5 <= startposatlevel[k][k2]){
		  /*  if(m1==1966){
		      cout<<"chunkpos["<<l<<"]["<<m1<<"]["<<k3<<"]["<<k4<<"]+"<<k5<<" "<< chunkpos[l][m1][k3][k4]+k5<<" startposatlevel["<<k<<"]["<<k2<<"] "<<startposatlevel[k][k2]<<endl;
		      cout<<"!!!m1=1966 ndummyends["<<l<<"]["<<m1<<"]["<<k3<<"] == nchunkcluster["<<l<<"]["<<m1<<"]["<<k3<<"] "<<ndummyends[l][m1][k3] <<" "<< nchunkcluster[l][m1][k3]<<endl;
		      } */

		  postmp[k2]=chunkpos[l][m1][k3][k4]+k5;
		  chunktmp[k2]=k4;
		  if(m1==1966){
		    //    cout<< "postmp[k2] "<<postmp[k2]<<" chunktmp[k2] "<<chunktmp[k2]<<" nvbstartvt[l] "<<nvbstartvt[l] <<"  levtmp["<<k2<<"] "<< levtmp[k2]<<endl;
		  }
		}
		if((chunkpos[l][m1][k3][k4]+k5 >= startposatlevel[k][k2]) || (k5==chunklen[l][m1][k3][k4]-1 && k4==nchunks[l][m1][k3]-1)){
		  found2=1;

		  if((m1==1966 ) ){
		    //  cout<<"found2 "<<found2<<endl;
		    //cout<<"in break:  levtmp["<<k2<<"] "<< levtmp[k2]<<endl;
		  }

		  break;
		}
	      }
	      if((m1==1966 ) ){
		//cout<<"in break2:  levtmp["<<k2<<"] "<< levtmp[k2]<<endl;
	      }

	      if(found2==1){
		break;
	      }

	    }
	    if((m1==1966 ) ){
	      //cout<<"in break3:  levtmp["<<k2<<"] "<< levtmp[k2]<<endl;
	    }
	    break;
	  }
	  if((m1==1966 ) ){
	    //  cout<<"after break3: "<<"postmp[k2] "<<postmp[k2]<<" chunktmp[k2] "<<chunktmp[k2]<<" nvbstartvt[l] "<<nvbstartvt[l] <<"  levtmp["<<k2<<"] "<< levtmp[k2]<<endl;

	  }
	  //	  	  	  cout<<k<<" "<<startposatlevel[k][k2]<<" "<<k2<<endl;
	  //	  	  	  cout<<"clustertmp "<<clustertmp<<endl;
	  //	  	  	  cout<<"levtmp["<<k2<<"] "<<levtmp[k2]<<endl;


	  ndummyatlevel[l][m1][levtmp[k2]][clustertmp]++;

	  //			  	  	  	  cout<<"ndummyatlevel["<<l<<"]["<<m1<<"]["<<levtmp[k2]<<"]["<<clustertmp<<"] "<<ndummyatlevel[l][m1][levtmp[k2]][clustertmp]<<" "<<startposatlevel[k][k2]<<flush<< endl;//ndummyatlevel[l][m1][levtmp[k2]][chunktmp[k2]]<<endl;

	  dummylevelatlevel[levtmp[k2]][clustertmp][ndummyatlevel[l][m1][levtmp[k2]][clustertmp]-1]=k;
	  //	  cout<<" l "<<l<<" m1 "<<m1<<endl;
	  //	  cout<<"dummylevelatlevel["<<levtmp[k2]<<"]["<<clustertmp<<"]["<<ndummyatlevel[l][m1][levtmp[k2]][clustertmp-1]<<"] "<<dummylevelatlevel[levtmp[k2]][clustertmp][ndummyatlevel[l][m1][levtmp[k2]][clustertmp]-1]<<endl;
	  dummyposatlevel[levtmp[k2]][clustertmp][ndummyatlevel[l][m1][levtmp[k2]][clustertmp]-1]=postmp[k2]; // startposatlevel[k][k2];
	  dummyclusteratlevel[levtmp[k2]][clustertmp][ndummyatlevel[l][m1][levtmp[k2]][clustertmp]-1]=k2;
	  if(m1==2 && (k==20 || k==43 || k==33 || k==5)){
	  //    if(m1==1966 && levtmp[k2]==0){
	    //	    cout <<"ndummyatlevel["<<l<<"]["<<m1<<"]["<<levtmp[k2]<<"]["<<clustertmp<<"] "<<ndummyatlevel[l][m1][levtmp[k2]][clustertmp]<<endl;
	    //	    cout<<"dummylevelatlevel["<<levtmp[k2]<<"]["<<clustertmp<<"]["<<ndummyatlevel[l][m1][levtmp[k2]][clustertmp]-1<<"] "<<dummylevelatlevel[levtmp[k2]][clustertmp][ndummyatlevel[l][m1][levtmp[k2]][clustertmp]-1]<<endl;
	    //	    cout<<"dummyposatlevel["<<levtmp[k2]<<"]["<<clustertmp<<"]["<<ndummyatlevel[l][m1][levtmp[k2]][clustertmp]-1<<"] "<<dummyposatlevel[levtmp[k2]][clustertmp][ndummyatlevel[l][m1][levtmp[k2]][clustertmp]-1]<<endl;
	    //	    cout<<"dummyclusteratlevel["<<levtmp[k2]<<"]["<<clustertmp<<"]["<<ndummyatlevel[l][m1][levtmp[k2]][clustertmp]-1<<"] "<<  dummyclusteratlevel[levtmp[k2]][clustertmp][ndummyatlevel[l][m1][levtmp[k2]][clustertmp]-1]<<endl;
	  }



	}


	delete[] postmp;
	delete[] levtmp;
	delete[] chunktmp;
      }


      for(int k=0; k<nvblevel[l][m1]; k++){
	for(int k0=0; k0<nchunkcluster[l][m1][k]; k0++){
	  if(ndummyatlevel[l][m1][k][k0]>0){
	    dummypos[l][m1][k][k0]=new int[ndummyatlevel[l][m1][k][k0]]();
	    if(!dummypos[l][m1][k][k0]) die("Memory allocation error in dummypos[l][m1][k][k0]!");
	    dummylevel[l][m1][k][k0]=new int[ndummyatlevel[l][m1][k][k0]]();
	    if(!dummylevel[l][m1][k][k0]) die("Memory allocation error in dummylevel[l][m1][k][k0]!");
	    dummycluster[l][m1][k][k0]=new int[ndummyatlevel[l][m1][k][k0]]();
	    if(!dummycluster[l][m1][k][k0]) die("Memory allocation error in dummycluster[l][m1][k][k0]!");
	    int *dummypostmp=new int[ndummyatlevel[l][m1][k][k0]]();
	    if(!dummypostmp) die("Memory allocation error in dummypostmp!");
	    int *dummyclustertmp=new int[ndummyatlevel[l][m1][k][k0]]();
	    if(!dummyclustertmp) die("Memory allocation error in dummyclustertmp!");
	    int *dummyclustersrttmp=new int[ndummyatlevel[l][m1][k][k0]]();
	    if(!dummyclustersrttmp) die("Memory allocation error in dummyclustersrttmp!");
	    int *dummypossrttmp=new int[ndummyatlevel[l][m1][k][k0]]();
	    if(!dummypossrttmp) die("Memory allocation error in dummypossrttmp!");
	    int *dummyleveltmp=new int[ndummyatlevel[l][m1][k][k0]]();
	    if(!dummyleveltmp) die("Memory allocation error in dummyleveltmp!");
	    int *dummylevelsrttmp=new int[ndummyatlevel[l][m1][k][k0]]();
	    if(!dummylevelsrttmp) die("Memory allocation error in dummysrttmp!");
	    for(int k1=0; k1<ndummyatlevel[l][m1][k][k0]; k1++){
	      //cout<<"dummyposatlevel["<<k<<"]["<<k1<<"]["<<k0<<"] "<<dummyposatlevel[k][k0][k1]<<endl;
	      //	      cout<<"dummylevelatlevel["<<k<<"]["<<k0<<"]["<<k1<<"] "<<dummylevelatlevel[k][k0][k1]<<endl;
	      //	      cout<<"dummyclusteratlevel["<<k<<"]["<<k0<<"]["<<k1<<"] "<<dummyclusteratlevel[k][k0][k1]<<endl;
	      dummypossrttmp[k1]=dummyposatlevel[k][k0][k1];
	      dummypostmp[k1]=dummyposatlevel[k][k0][k1];
	      //dummypossrttmp[k1]=dummyposatlevel[k][k1][k0];
	      dummyleveltmp[k1]=dummylevelatlevel[k][k0][k1];
	      //cout<<"dummyleveltmp["<<k1<<"] "<<dummyleveltmp[k1]<<" dummylevelatlevel["<<k<<"]["<<k0<<"]["<<k1<<"] !!!"<<dummylevelatlevel[k][k0][k1]<<endl;
	      dummyclustertmp[k1]=dummyclusteratlevel[k][k0][k1];
	      if(m1==2 && k==20){
		//	              cout<<"dummyposatlevel["<<l<<"]["<<m1<<"]["<<k<<"]["<<k1<<"] "<<dummyposatlevel[k][k0][k1];
		//		      cout<<" dummylevelatlevel["<<l<<"]["<<m1<<"]["<<k<<"]["<<k1<<"] "<<dummylevelatlevel[k][k0][k1]<<endl;
		//		      cout<<"dummyclusteratlevel["<<k<<"]["<<k0<<"]["<<k1<<"] "<<dummyclusteratlevel[k][k0][k1]<<endl;
	      }
	    }
	    // dummypos[] is not necessarily in ascending order, which is required for testing. Sort dummypos:
	    //cout<<"dummypossrttmp["<<0<<"]! "<<dummypossrttmp[0]<<endl;
	    qsort(dummypossrttmp, ndummyatlevel[l][m1][k][k0], sizeof(int), compareint);
	    //cout<<"dummypossrttmp["<<0<<"]! "<<dummypossrttmp[0]<<endl;

	    for(int k1=0; k1<ndummyatlevel[l][m1][k][k0]; k1++){
	      //cout<<"dummypossrttmp["<<k1<<"]! "<<dummypossrttmp[k1]<<endl;
	      //     cout<<"ndummyatlevel["<<l<<"]["<<m1<<"]["<<k<<"] "<<  ndummyatlevel[l][m1][k]<<endl;
	      for(int k2=0; k2<ndummyatlevel[l][m1][k][k0]; k2++){
		//      cout<<"dummypossrttmp["<<k1<<"] dummypostmp["<<k2<<"] "<< dummypossrttmp[k1]<<" "<<dummypostmp[k2]<<endl;
		//  cout<<dummypossrttmp[k1]<<" "<<dummypostmp[k2]<< " pre if "<<endl;
		if(dummypossrttmp[k1]==dummypostmp[k2] && dummypostmp[k2]!=-1){
		  dummylevelsrttmp[k1]=dummyleveltmp[k2];
		  //    cout<<"dummylevelsrttmp["<<k1<<"] "<<dummylevelsrttmp[k1]<<" dummyleveltmp["<<k2<<"] "<<dummyleveltmp[k2]<<" ?!?"<<endl;
		  dummyclustersrttmp[k1]=dummyclustertmp[k2];
		  dummypostmp[k2]=-1;
		  break;
		}
	      }
	      //cout<<"dummypossrttmp["<<k1<<"] "<<dummypossrttmp[k1]<<endl;
	      dummypos[l][m1][k][k0][k1]=dummypossrttmp[k1];
	      //cout<<"dummylevelsrttmp["<<k1<<"] "<<dummylevelsrttmp[k1]<<endl;
	      dummylevel[l][m1][k][k0][k1]=dummylevelsrttmp[k1];
	      //cout<<"dummylevel["<<l<<"]["<<m1<<"]["<<k<<"]["<<k0<<"]["<<k1<<"] "<<dummylevel[l][m1][k][k0][k1]<<endl;
	      dummycluster[l][m1][k][k0][k1]=dummyclustersrttmp[k1];
	    }
	    delete[] dummypostmp;
	    delete[] dummyclustertmp;
	    delete[] dummyclustersrttmp;
	    delete[] dummypossrttmp;
	    delete[] dummyleveltmp;
	    delete[] dummylevelsrttmp;
	  }
	  else if(ndummyatlevel[l][m1][k][k0]==0){

	    dummypos[l][m1][k][k0]=new int[1]();
	    if(!dummypos[l][m1][k][k0]) die("Memory allocation error in dummypos[l][m1][k][k0]!");
	    dummylevel[l][m1][k][k0]=new int[1]();
	    if(!dummylevel[l][m1][k][k0]) die("Memory allocation error in dummylevel[l][m1][k][k0]!");
	    dummycluster[l][m1][k][k0]=new int[1]();
	    if(!dummycluster[l][m1][k][k0]) die("Memory allocation error in dummycluster[l][m1][k][k0]!");
	  }
	}
      }

      //delete[] atlevelCarriers;
      //      delete[] atlevelCarriersInverse;
      //      delete[] atlevelsorted;



      for(int k1=0; k1<nvblevel[l][m1]; k1++){
	for(int k=0; k<ndummyends[l][m1][k1]+1; k++){
	  for(int k2=0; k2<nvblevel[l][m1]; k2++){
	    //	    cout<<"dummyposatlevel["<<k1<<"]["<<k<<"]["<<k2<<"] "<<dummyposatlevel[k1][k][k2]<<endl;
	  }
	  //	  cout<<m1<<" "<<k1<<" "<<k<<endl;
	  delete[] dummyposatlevel[k1][k];
	  delete[] dummylevelatlevel[k1][k];
	  delete[] dummyclusteratlevel[k1][k];
	}
	delete[] dummyposatlevel[k1];
	delete[] dummylevelatlevel[k1];
	delete[] dummyclusteratlevel[k1];
      }
      delete[] dummyclusteratlevel;
      delete[] dummyposatlevel;
      delete[] dummylevelatlevel;

      for(int m3=0; m3<nvblevel[l][m1]; m3++){
	//	cout<<m1<<" "<<nvblevel[l][m1] <<" "<<m3<<flush<<endl;
	if(startposatlevel[m3]!=NULL)delete[] startposatlevel[m3];
      }
      delete[] startposatlevel;
      delete[] atlevelInverse;
      delete[] atlevel;

      delete[] atlevelInverseStatus;
      delete[] Ind1;
      delete[] Ind2;
      delete[] nRemaining;



      for(int m=0; m<nRareLimits[l]; m++){
	delete[] dummy1[m];
	delete[] dummy2[m];
      }
      delete[] dummy1;
      delete[] dummy2;

    } // loop over start variants
  } // loop over nwindows

    // count distint bins

  unsigned int ndistinct=0;
  for(int l=0; l<nwindows; l++){
    unsigned int ndistinctchr=0;
    int nlevels=0;
    for(int m=0; m<nvbstartvt[l]; m++){
      for(int m3=0; m3<nvblevel[l][m]; m3++){
	nlevels++;
	for(int m4=0; m4<nchunks[l][m][m3]; m4++){
	  ndistinctchr+=chunklen[l][m][m3][m4];
	  ndistinct+=chunklen[l][m][m3][m4];
	}
      }
    }
    cout<<"On chromosome "<<map[vbstartvt[l][chunkpos[l][0][0][0]]].chr<<", "<<nvbstartvt[l]<<" distinct variants result in "<<ndistinctchr<<" distinct bins at overall "<<nlevels <<" carrier threshold levels"<<endl;
  }
  cout<<"Distinct bins overall: "<<ndistinct<<endl;
}


void calcVB(int optimalrare, int ncases, int ncontrols, struct MAP *map, uint64_t*** BinSNPs, uint64_t** BinSNPsCCFlags, int nwords, int* SNPMapInverse,  int n, int *nvbstartvt, int **vbstartvt, int *nvbend, int **vbend,  float *vbmaxpermstat, int **nvblevel, int ***nchunks, int ****chunkpos, int ****chunklen, int nwindows, int *nCarriers, int ****ndummyatlevel, int *****dummypos, int *****dummylevel, uint64_t ***BinCarriers, int ****dummyend, int ***ndummyends, int ***nchunkcluster, int *****dummycluster, int vb_binwise_corr, float *****vbbinwisestat, short int *****vbbinwisecount){

  int sX, sY;
  double stat;

  for(int i=0; i<nwindows; i++){
    for(int j=0; j<nvbstartvt[i]; j++){
      // Adaptive testing
      float statmax_pos=0;
      float statmax_neg=0;
      int adaptiveVB=0; // Experimental adaptive VB testing
      int effdir=0;
      int firsteffdir=0;

      uint64_t*** dummy = new uint64_t**[nvblevel[i][j]];
      if(!dummy) die("Memory allocation error in dummy!");
      // cout<<"nchunkcluster["<<i<<"]["<<j<<"]["<<0<<"] "<<nchunkcluster[i][j][0]<<" "<<nvblevel[i][j]<<endl;
      for(int k=0; k<nvblevel[i][j]; k++){
	dummy[k]=new uint64_t*[nchunkcluster[i][j][k]];
	if(!dummy[k]) die("Memory allocation error in dummy[k]!");
	for(int k1=0; k1<nchunkcluster[i][j][k]; k1++){
	  dummy[k][k1]=new uint64_t[nwords]();
	if(!dummy[k][k1]) die("Memory allocation error in dummy[k][k1]!");
	}
      }

      int *carriersatlevel_order=NULL;
      if(optimalrare){
	int *carriersatlevel_unsort=new int[nvblevel[i][j]]();
	int *carriersatlevel=new int[nvblevel[i][j]]();
	carriersatlevel_order=new int[nvblevel[i][j]]();
	for(int k=0; k<nvblevel[i][j]; k++){
	  carriersatlevel_unsort[k]=nCarriers[vbstartvt[i][chunkpos[i][j][k][0]]];
	  carriersatlevel[k]=nCarriers[vbstartvt[i][chunkpos[i][j][k][0]]];
	}
	qsort(carriersatlevel, nvblevel[i][j], sizeof(int), compareint);
	for(int k=1; k<nvblevel[i][j]; k++){
	  for(int kk=1; kk<nvblevel[i][j]; kk++){
	    if(carriersatlevel[k]==carriersatlevel_unsort[kk]){
	      carriersatlevel_order[k]=kk;
	      //		cout<<i<<" "<<j<<" "<<kk<<endl;
	      break;
	    }
	  }
	  //	    cout<<i<<" "<<j<<" "<<k<<" "<<carriersatlevel_order[k]<<endl;
	}
	delete[] carriersatlevel_unsort;
	delete[] carriersatlevel; // to delete: carriers.levelorder
      }

      for(int kk=0; kk<nvblevel[i][j]; kk++){
	int k=kk;
	if(optimalrare){
	  k=carriersatlevel_order[kk];
	}

	int *founddummypos=new int[nchunkcluster[i][j][k]]();
	if(!founddummypos) die("Memory allocation error in founddummypos!");
	int cluster=0;
	for(int l=0; l<nchunks[i][j][k]; l++){
	  if(nchunkcluster[i][j][k]>1){
	    if(cluster<nchunkcluster[i][j][k]-1 && chunkpos[i][j][k][l]>dummyend[i][j][k][cluster]){
	      cluster++;
	    }
	  }
	  for(int m=0; m<chunklen[i][j][k][l]; m++){
	    int pos=chunkpos[i][j][k][l]+m;
	    int i2=vbstartvt[i][pos];
	    sX=0;
	    sY=0;
	    int iMod = SNPMapInverse[i2];
	    if (iMod==-1) continue;
	    for (int p=0; p<nwords; p++){
	      dummy[k][cluster][p] |= BinCarriers[i][pos][p];
	    }
	    while(ndummyatlevel[i][j][k][cluster]>0 && pos==dummypos[i][j][k][cluster][founddummypos[cluster]] && founddummypos[cluster]<ndummyatlevel[i][j][k][cluster]){
	      for (int p=0; p<nwords; p++){
		dummy[dummylevel[i][j][k][cluster][founddummypos[cluster]]][dummycluster[i][j][k][cluster][founddummypos[cluster]]][p]=dummy[k][cluster][p];
	      }
	      founddummypos[cluster]++;
	    }

	    for(int p=0; p<nwords; p++){
	      sY += bitcount64(dummy[k][cluster][p] & BinSNPsCCFlags[p][1]);
	      sX += bitcount64(dummy[k][cluster][p] & BinSNPsCCFlags[p][2]);
	    }
	    stat=((double)ncases + (double)ncontrols)*((double)sX*(double)ncontrols - (double)sY*(double)ncases)*((double)sX*(double)ncontrols-(double)sY*(double)ncases)/((double)ncases*(double)ncontrols*((double)sX+(double)sY)*((double)ncases+(double)ncontrols-(double)sX-(double)sY));

	    // Adaptive testing
	    if(adaptiveVB){
	      float propNUC=float(sY)/float(ncontrols);
	      float propNAC=float(sX)/float(ncases);
	      if((propNUC-propNAC)>EPS){ // protective
		effdir=-1;
		if((stat-statmax_neg)>EPS){
		  statmax_neg=stat;
		}
	      }
	      else if((propNAC-propNUC)>EPS){ // damaging
		effdir=1;
		if((stat-statmax_pos)>EPS){
		  statmax_pos=stat;
		}
	      }
	      if(j==0){
		firsteffdir=effdir;
	      }
	      else if(j>0){
		if(firsteffdir==1 && ((statmax_neg-statmax_pos)>EPS /* || fabs(statmax_neg-statmax_pos)<EPS */)){
		  break;
		}
		else if(firsteffdir==-1 && ((statmax_pos-statmax_neg)>EPS /* || fabs(statmax_neg-statmax_pos)<EPS */)){
		  break;
		}

	      }
	    }
	    if((ncases+ncontrols-sY-sX)==0 || (sX == 0 && sY==0)) stat=0;
	    if(n==0 && vb_binwise_corr){
	      vbbinwisestat[i][j][k][l][m]=stat;

	      //      vbstat[k][j]=stat;


	      /*if(stat>30){
		cout<<map[vbend[k][0]].chr<<"\t"<<map[vbend[k][0]].pos<<"\t"<<map[vbend[k][j]].pos<<"\t"<< pValueCalc(0.5, float(vbstat[k][j])/float(2))<<endl;
		}
		if(sX==0 || sY==0 || ncases==sX){
		OR_COLL=-9999;
		}
		else{
		OR_COLL=(double)sX*((double)ncontrols-(double)sY)/((double)sY*((double)ncases-(double)sX));
		}*/
	    }
	    else if(n>0){
	      if((stat-vbmaxpermstat[n-1])>EPS){
		vbmaxpermstat[n-1]=stat;
	      }
	      if(vb_binwise_corr){
		if(stat*(1+EPS)>vbbinwisestat[i][j][k][l][m]){
		  vbbinwisecount[i][j][k][l][m]++;
		}
	      }
	    }
	  }
	}
	if(ndummyatlevel[i][j][k][cluster]!=founddummypos[cluster]){
	  cout<<" ndummyatlevel["<<i<<"]["<<j<<"]["<<k<<"]["<<cluster<<"] "<<ndummyatlevel[i][j][k][cluster]<<" "<<founddummypos[cluster]<<endl;
	  cout<<"BUG in ndummyatlevel!"<<endl;
	  exit(1);
	}
	delete[] founddummypos;
      }
      delete[] carriersatlevel_order;

      for(int k=0; k<nvblevel[i][j]; k++){
	for(int k1=0; k1<nchunkcluster[i][j][k]; k1++){
	  delete[] dummy[k][k1];
	}
	delete[] dummy[k];
      }
      delete[] dummy;
    }
  }
}

  /// Collapsing test
  double calcCOLL_bin(int ncases, int ncontrols, int start, int end, struct MAP *map, double raref, double &OR_COLL, uint64_t*** BinSNPs, uint64_t** BinSNPsCCFlags, int nwords, int* SNPMapInverse, struct WINDOW *window, int l, int m, int n, int nwindowssinglebin, int collinter, int **doublewindowcoord){
    int sX=0,sY=0;
    double stat=0;

    uint64_t* dummy = new uint64_t[nwords];
    if(!dummy) die("Memory allocation error in dummy!");
    memset(dummy, 0, nwords*sizeof(uint64_t));

    if(collinter==2 && l>=nwindowssinglebin){
      uint64_t* dummy1 = new uint64_t[nwords];
      if(!dummy1) die("Memory allocation error in dummy1!");
      uint64_t* dummy2 = new uint64_t[nwords];
      if(!dummy2) die("Memory allocation error in dummy2!");
      memset(dummy1, 0, nwords*sizeof(uint64_t));
      memset(dummy2, 0, nwords*sizeof(uint64_t));
      int start1, end1, start2, end2, m1, m2;
      m1=doublewindowcoord[l][0];
      m2=doublewindowcoord[l][1];
      start1 = windowPositions[m1][0];
      end1 = windowPositions[m1][1];
      start2 = windowPositions[m2][0];
      end2 = windowPositions[m2][1];
      for(int j=0; j<window[l].n_at_level[m]; j++){
	int i=window[l].levelpos[m][j];
	int iMod = SNPMapInverse[i];
	if (iMod==-1) continue;
	for (int p=0; p<nwords; p++){
	  if(i>=start1 && i<=end1){
	    dummy1[p] |= (BinSNPs[iMod][p][1] | BinSNPs[iMod][p][2]);
	  }
	  if(i>=start2 && i<=end2){
	    dummy2[p] |= (BinSNPs[iMod][p][1] | BinSNPs[iMod][p][2]);
	  }
	  dummy[p] |= dummy1[p] & dummy2[p];
	}
      }
      delete[] dummy1;
      delete[] dummy2;
    }
    else{
      for(int j=0; j<window[l].n_at_level[m]; j++){
	int i=window[l].levelpos[m][j];
	int iMod = SNPMapInverse[i];
	if (iMod==-1) continue;
	for (int p=0; p<nwords; p++)
	  dummy[p] |= (BinSNPs[iMod][p][1] | BinSNPs[iMod][p][2]);
      }
    }
    for(int p=0; p<nwords; p++){
      sY += bitcount64(dummy[p] & BinSNPsCCFlags[p][1]);
      sX += bitcount64(dummy[p] & BinSNPsCCFlags[p][2]);
    }


    delete[] dummy;
    stat=((double)ncases + (double)ncontrols)*((double)sX*(double)ncontrols - (double)sY*(double)ncases)*((double)sX*(double)ncontrols-(double)sY*(double)ncases)/((double)ncases*(double)ncontrols*((double)sX+(double)sY)*((double)ncases+(double)ncontrols-(double)sX-(double)sY));
    if((ncases+ncontrols-sY-sX)==0 || (sX == 0 && sY==0)) stat=0;
    if(n==0){
      if(sX == 0 || sY==0 || ncases==sX){
	OR_COLL=-9999;
      }
      else{
	OR_COLL=(double)sX*((double)ncontrols-(double)sY)/((double)sY*((double)ncases-(double)sX));
      }
    }


    return stat;
  }


  double calcCOLL_X_bin(int ncases, int ncontrols, int start, int end, struct MAP *map, double raref, double &OR_COLL, double &OR_COLL_f, uint64_t*** BinSNPs, uint64_t** BinSNPsCCFlags, int nwords, int* SNPMapInverse, uint64_t** GenderFlags,int ncasesqcMale, int ncasesqcFemale,int ncontrolsqcMale, int ncontrolsqcFemale, struct WINDOW *window, int l, int m, int n) {
    int sX_M=0,sY_M=0;
    double stat_M=0;
    int sX_F=0,sY_F=0;
    double stat_F=0;
    uint64_t* dummy = new uint64_t[nwords];
    if(!dummy) die("Memory allocation error in dummy!");
    memset(dummy, 0, nwords*sizeof(uint64_t));
    int iMod;
    for(int j=0; j<window[l].n_at_level[m]; j++){
      int i=window[l].levelpos[m][j];

      iMod = SNPMapInverse[i];
      if (iMod==-1) continue;
      for (int p=0; p<nwords; p++)
	dummy[p] |= (BinSNPs[iMod][p][1] | BinSNPs[iMod][p][2]);
    }
    for (int p=0; p<nwords; p++) {
      sY_M += bitcount64(dummy[p] & BinSNPsCCFlags[p][1] & GenderFlags[p][1]);
      sX_M += bitcount64(dummy[p] & BinSNPsCCFlags[p][2] & GenderFlags[p][1]);
      sY_F += bitcount64(dummy[p] & BinSNPsCCFlags[p][1] & GenderFlags[p][2]);
      sX_F += bitcount64(dummy[p] & BinSNPsCCFlags[p][2] & GenderFlags[p][2]);
    }
    delete[] dummy;
    stat_M=((double)ncasesqcMale + (double)ncontrolsqcMale)*((double)sX_M*(double)ncontrolsqcMale - (double)sY_M*(double)ncasesqcMale)*((double)sX_M*(double)ncontrolsqcMale-(double)sY_M*(double)ncasesqcMale)/((double)ncasesqcMale*(double)ncontrolsqcMale*((double)sX_M+(double)sY_M)*((double)ncasesqcMale+(double)ncontrolsqcMale-(double)sX_M-(double)sY_M));
    if((ncasesqcMale+ncontrolsqcMale-sY_M-sX_M)==0 || (sX_M == 0 && sY_M==0)) stat_M=0; //??

    stat_F=((double)ncasesqcFemale + (double)ncontrolsqcFemale)*((double)sX_F*(double)ncontrolsqcFemale - (double)sY_F*(double)ncasesqcFemale)*((double)sX_F*(double)ncontrolsqcFemale-(double)sY_F*(double)ncasesqcFemale)/((double)ncasesqcFemale*(double)ncontrolsqcFemale*((double)sX_F+(double)sY_F)*((double)ncasesqcFemale+(double)ncontrolsqcFemale-(double)sX_F-(double)sY_F));
    if((ncasesqcFemale+ncontrolsqcFemale-sY_F-sX_F)==0 || (sX_F == 0 && sY_F==0)) stat_F=0;

    if(n==0){
      if(sX_M == 0 || sY_M==0 || ncasesqcMale==sX_M){
	OR_COLL=-9999;
      }
      else{
	OR_COLL=(double)sX_M*((double)ncontrolsqcMale-(double)sY_M)/((double)sY_M*((double)ncasesqcMale-(double)sX_M));
      }
      if(sX_F == 0 || sY_F==0 || ncasesqcFemale==sX_F){
	OR_COLL_f=-9999;
      }
      else{
	OR_COLL_f=(double)sX_F*((double)ncontrolsqcFemale-(double)sY_F)/((double)sY_F*((double)ncasesqcFemale-(double)sX_F));
      }
    }
    return stat_M+stat_F;
  }



  double calcCMAT(int *cov, int sexcov, struct PERSON *person, int ncases, int ncontrols, int nrest, int thread, int nlinestfam, int nlinestped, int start, int end, double &OR_CMAT, struct MAP *map, double raref,  int *nRareSNPs, uint64_t*** BinSNPs, uint64_t** BinSNPsCCFlags, int nwordsSNPs, int* SNPMapInverse, double *weightvec , struct WINDOW *window, int l, int m, int n){
    int i;
    double stat=0;

    float mA=0;   // minor/major allele counts in affected/unaffected individuals
    float mU=0;
    float MA=0;
    float MU=0;
    float sumw=0;

    //for(i = start; i <= end; i++){
    for(int j=0; j<window[l].n_at_level[m]; j++){
      i=window[l].levelpos[m][j];

      int iMod = SNPMapInverse[i];
      if (iMod==-1) continue;
      for (int p=0; p<nwordsSNPs; p++) {
	mA += weightvec[i]*(2*(bitcount64(BinSNPs[iMod][p][1] & BinSNPsCCFlags[p][2])) +   (bitcount64(BinSNPs[iMod][p][2] & BinSNPsCCFlags[p][2])));
	MA += weightvec[i]*((bitcount64(BinSNPs[iMod][p][2] & BinSNPsCCFlags[p][2])) + 2*(bitcount64(BinSNPs[iMod][p][3] & BinSNPsCCFlags[p][2])));
	mU += weightvec[i]*(2*(bitcount64(BinSNPs[iMod][p][1] & BinSNPsCCFlags[p][1])) +   (bitcount64(BinSNPs[iMod][p][2] & BinSNPsCCFlags[p][1])));
	MU += weightvec[i]*((bitcount64(BinSNPs[iMod][p][2] & BinSNPsCCFlags[p][1])) + 2*(bitcount64(BinSNPs[iMod][p][3] & BinSNPsCCFlags[p][1])));
      }
      sumw=(float)sumw+weightvec[i];
    }


    if(nRareSNPs[l]==0 || sumw==0){
      stat=-1;
    }
    else{
      stat = (((double)ncases + (double)ncontrols)/(2*(double)ncases*(double)ncontrols*(double)sumw))*((double)mA*(double)MU - (double)mU*(double)MA)* ((double)mA*(double)MU - (double)mU*(double)MA)/(((double)mA + (double)mU)*((double)MA + (double)MU));

      if(n==0 && MA!=0 && mU!=0){
	OR_CMAT=(double)MU*(double)mA/((double)MA*(double)mU);
      }
      else{
	OR_CMAT=-9999;
      }
    }


    return stat;
  };

  double calcCMATcov(int *cov,int numberOfAllCov, int sexcov, struct PERSON *person, int ncases, int ncontrols, int nrest, int thread, int nlinestfam, int nlinestped, int start, int end, double &OR_CMAT, struct MAP *map, double raref,  int *nRareSNPs, uint64_t*** BinSNPs, uint64_t** BinSNPsCCFlags, uint64_t** BinSNPsCovCatFlags, int nwordsSNPs, int* SNPMapInverse, double *weightvec , struct WINDOW *window, int l, int m, int n, int nCovCathegories){

    int i,k;
    double stat=0;
    float mA[nCovCathegories];   // minor/major allele counts in affected/unaffected individuals
    float mU[nCovCathegories];
    float MA[nCovCathegories];
    float MU[nCovCathegories];
    int NAC[nCovCathegories];
    int NUC[nCovCathegories];
    float sumw=0;

    for (int c=0; c<nCovCathegories; c++){
      mA[c]=0;
      mU[c]=0;
      MA[c]=0;
      MU[c]=0;
      NAC[c]=0;
      NUC[c]=0;
    }

    for(int j=0; j<window[l].n_at_level[m]; j++){
      i=window[l].levelpos[m][j];

      int iMod = SNPMapInverse[i];
      if (iMod==-1) continue;
      for (int c=0; c<nCovCathegories; c++){
	for (int p=0; p<nwordsSNPs; p++) {
	  mA[c] += weightvec[i]*(2*(bitcount64(BinSNPs[iMod][p][1] & BinSNPsCCFlags[p][2] & BinSNPsCovCatFlags[p][c])) + ((bitcount64(BinSNPs[iMod][p][2] & BinSNPsCCFlags[p][2])) & BinSNPsCovCatFlags[p][c]));
	  MA[c] += weightvec[i]*((bitcount64(BinSNPs[iMod][p][2] & BinSNPsCCFlags[p][2] & BinSNPsCovCatFlags[p][c])) + 2*(bitcount64(BinSNPs[iMod][p][3] & BinSNPsCCFlags[p][2] & BinSNPsCovCatFlags[p][c])));
	  mU[c] += weightvec[i]*(2*(bitcount64(BinSNPs[iMod][p][1] & BinSNPsCCFlags[p][1] & BinSNPsCovCatFlags[p][c])) +  (bitcount64(BinSNPs[iMod][p][2] & BinSNPsCCFlags[p][1] & BinSNPsCovCatFlags[p][c])));
	  MU[c] += weightvec[i]*((bitcount64(BinSNPs[iMod][p][2] & BinSNPsCCFlags[p][1] & BinSNPsCovCatFlags[p][c])) + 2*(bitcount64(BinSNPs[iMod][p][3] & BinSNPsCCFlags[p][1] & BinSNPsCovCatFlags[p][c])));
	  NAC[c]+=BinSNPsCCFlags[p][2] & BinSNPsCovCatFlags[p][c];
	  NUC[c]+=BinSNPsCCFlags[p][1] & BinSNPsCovCatFlags[p][c];
	}
      }
      sumw=(float)sumw+weightvec[i];
    }

    if(nRareSNPs[l]==0 || sumw==0){
      stat=-1;
    }
    else{
      double sum1=0;
      double sum2=0;

      for (int c=0; c<nCovCathegories; c++){
	if(NAC[c]+NUC[c]>0){
	  sum1+=(double)mA[c]-((double)NAC[c]*((double)mA[c]+(double)mU[c])/double(NAC[c]+NUC[c]));
	  sum2+=(double)NAC[c]*NUC[c]*(mA[c]+mU[c])*(MA[c]+MU[c])/(2*double(NAC[c]+NUC[c])*double(NAC[c]+NUC[c])*double(NAC[c]+NUC[c])*sumw);
	}
      }
      stat=(double)sum1*(double)sum1/(double)sum2;
      OR_CMAT=-9999;
    }


    return stat;
  };



void windowsRARE(int bin, int n, int nlinestped, int binsizeRare, struct MAP *map, struct COUNTS **counts, int *nwindows, int ***nchrwindows, int thread, fstream &errorfile, fstream &logfile, string intervalfile, double raref, string intervalfile_format, int mergelines, int minRareInBin, int maxRareInBin, int binamin, int binamax, int featurecol, string SetIDfile, int intervaleditor, int setid, int verbose,  int vb, int NCT, uint64_t ***BinSNPs, int *SNPMapInverse, int nwords)
  {
    cout    << "\n\nBegin rare variant analysis..."<<endl;
    int i, j, l, m;

    *nchrwindows=calloc2Dint(26, 2);

    for(j=0; j<26; j++){
      nRareSNPsPerChr[j]=0;
    }

    int start=0;
    int currentchr; // current chromosome number
    int oldchr=0;

    double smallestmaf=0.5;

    // Check if RARE is not too small
    for(i=0; i<nlinestped; i++){
      if(map[i].mafr>EPS && (map[i].mafr-smallestmaf)<=EPS &&  map[i].analysis_in){
	smallestmaf=map[i].mafr;
      }
    }
    if((smallestmaf-raref)>EPS){
      cout << "RARE is too small. All MAFs are larger than "<< raref<<"! The smallest MAF is "<< smallestmaf << "!" <<endl;
      logfile << "RARE is too small. All MAFs are larger than "<< raref<<"! The smallest MAF is "<< smallestmaf << "!" <<endl;
      errorfile << "RARE is too small. All MAFs are larger than "<< raref<<"! The smallest MAF is "<< smallestmaf << "!" <<endl;
      exit(1);
    }

    if(intervalfile == " "){
      if(binamax!=0){
	cout<<"For DENSE_BINNING, an interval file has to be provided."<<endl;
	logfile<<"For DENSE_BINNING, an interval file has to be provided."<<endl;
	errorfile<<"For DENSE_BINNING, an interval file has to be provided."<<endl;
	exit(1);
      }
      for(i=0; i<nlinestped; i++){
	currentchr=atoi(map[i].chr);
	if(currentchr==0 && map[i].pos==0){
	  cout << map[i].rs<< "\tignored: chr=0, position=0" <<endl ;
	  logfile << map[i].rs<< "\tignored: chr=0, position=0" <<endl ;
	  continue;
	}
	else if(currentchr==0 && map[i].pos!=0){
	  cout << map[i].rs<< "\tignored: chr=0" <<endl ;
	  logfile << map[i].rs<< "\tignored: chr=0" <<endl ;
	  continue;
	}
	else if(map[i].pos==0){
	  cout << map[i].rs<< "\tignored: position=0" <<endl ;
	  logfile << map[i].rs<< "\tignored: position=0" <<endl ;
	  continue;
	}
	if(vb==0){
	  if(((map[i].mafr-raref)<EPS || fabs(map[i].mafr-raref)<EPS) && map[i].mafr>EPS && map[i].analysis_in == 1  && map[i].pos!=0 && atoi(map[i].chr)!=0 ){
	    nRareSNPsPerChr[currentchr-1]++;
	  }
	}
	else if(vb!=0){
	  uint64_t* dummy = new uint64_t[nwords]();
	  if(!dummy) die("Memory allocation error in dummy!");
	  //memset(dummy, 0, nwords*sizeof(uint64_t));
	  int iMod = SNPMapInverse[i];
	  if (iMod==-1) continue;
	  for (int p=0; p<nwords; p++){
	    dummy[p] |= (BinSNPs[iMod][p][1] | BinSNPs[iMod][p][2]);
	    nCarriers[i]+=bitcount64(dummy[p]);
	  }
	  if(nCarriers[i]<=NCT){
	    nRareSNPsPerChr[currentchr-1]++;
	  }
	  delete[] dummy;
	}
      }
      if(verbose==2 && vb==0){
	cout <<"\nRare SNPs per chromosome:"<<endl;
	cout<<"Chr\tRareSNPs"<<endl;
	for(i=0; i<26; i++){
	  if(nRareSNPsPerChr[i]!=0){
	    cout <<i+1<<"\t"<<nRareSNPsPerChr[i]<<endl;
	    logfile <<i+1<<"\t"<< nRareSNPsPerChr[i]<<endl;
	  }
	}
      }
      for(i=0; i<nlinestped; i++){
	currentchr = atoi(map[i].chr);
	if(oldchr > currentchr){

	  cout << "Chromosomes in wrong order! Please sort your *.tped file with respect to chromosomes and positions to perform a rare variants analysis." << endl;
	  cout << "Fist two offenders: "<<map[i-1].rs<< " on chr"<<map[i-1].chr<<":"<<map[i-1].pos<<" and "<<map[i].rs<< " on chr"<<map[i].chr<<":"<<map[i].pos<<"." <<endl;
	  logfile << "Chromosomes in wrong order! Please sort your *.tped file with respect to chromosomes and positions to perform a rare variants analysis." << endl;
	  logfile << "Fist two offenders: "<<map[i-1].rs<< " on chr"<<map[i-1].chr<<":"<<map[i-1].pos<<" and "<<map[i].rs<< " on chr"<<map[i].chr<<":"<<map[i].pos<<"." <<endl;
	  errorfile << "Chromosomes in wrong order! Please sort your *.tped file with respect to chromosomes and positions to perform a rare variants analysis." << endl;
	  errorfile << "Fist two offenders: "<<map[i-1].rs<< " on chr"<<map[i-1].chr<<":"<<map[i-1].pos<<" and "<<map[i].rs<< " on chr"<<map[i].chr<<":"<<map[i].pos<<"." <<endl;
	  errorfile.close();
	  logfile.close();
	  exit(1);
	}
	if(currentchr!=0 && map[i].pos!=0){
	  (*nchrwindows)[currentchr-1][0]++; // fill nchrwindows with number of SNPs first to convert to number of windows later
	}
	oldchr = currentchr;
      }
      if(verbose==2){
	cout <<"\nSNPs per chromosome:"<<endl;
	cout<<"Chr\tSNPs"<<endl;

	for(i=0; i<26; i++){
	  if((*nchrwindows)[i][0]!=0){
	    cout <<i+1<<"\t"<< (*nchrwindows)[i][0]<<endl;
	    logfile <<i+1<<"\t"<< (*nchrwindows)[i][0]<<endl;
	  }
	}
      }

      l=0;
      for(j=0; j<26; j++){
	for(i=0; i<(*nchrwindows)[j][0]; i++){
	  l++;
	}
      }
      cout <<"\nThere are "<< l << " SNPs with well-defined positions" << endl; // with non-zero chromosome and bp position
      logfile <<"\nThere are "<< l << " SNPs with well-defined positions" << endl;

      if(vb){
	for(j=0; j<26; j++){
	  if(nRareSNPsPerChr[j]!=0){
	    cout << "Rare variants (<="<<NCT<<" carriers per variant) on chromosome " << j+1 << ": "<< nRareSNPsPerChr[j] << endl;
	    logfile << "Rare variants (<="<<NCT<<" carriers per variant) on chromosome " << j+1 << ": "<< nRareSNPsPerChr[j] << endl;
	    (*nchrwindows)[j][1]=(*nchrwindows)[j][0];
	    (*nchrwindows)[j][0]=1;
	  }
	}

      }
      else if(binsizeRare!=0){
	for(j=0; j<26; j++){
	  //	  cout<<"before "<<nRareSNPsPerChr[j]<<" "<<(*nchrwindows)[j][0]<<" "<<(*nchrwindows)[j][1]<<endl;
	  if(binsizeRare > nRareSNPsPerChr[j] &&  nRareSNPsPerChr[j]!=0 ){
	    cout << "BinsizeRare larger than number of rare SNPs on chromosome " << j+1 << ", bin adjusted to length "<< nRareSNPsPerChr[j] << " rare SNPs" << endl;
	    logfile <<  "BinsizeRare larger than number of rare SNPs on chromosome " << j+1 << ", bin adjusted to length "<<  nRareSNPsPerChr[j] << " rare SNPs"<<endl;
	    (*nchrwindows)[j][1]=nRareSNPsPerChr[j];
	    (*nchrwindows)[j][0]=1;
	  }
	  else if(nRareSNPsPerChr[j]!=0){
	    if (( nRareSNPsPerChr[j] % binsizeRare) != 0){ // if an incomplete window remains at the end
	      (*nchrwindows)[j][1] = nRareSNPsPerChr[j] % binsizeRare;
	      (*nchrwindows)[j][0]=(nRareSNPsPerChr[j]/binsizeRare)+1;
	    }
	    else{
	      (*nchrwindows)[j][0]=nRareSNPsPerChr[j]/binsizeRare;
	      (*nchrwindows)[j][1]=0;
	    }
	  }
	  else if(nRareSNPsPerChr[j]==0){
	    (*nchrwindows)[j][0]=0;
	    (*nchrwindows)[j][1]=0;
	  }
	  //	  cout<<nRareSNPsPerChr[j]<<" "<<(*nchrwindows)[j][0]<<" "<<(*nchrwindows)[j][1]<<endl;
	}
      }

      // calculate Nr bins
      *nwindows=0;
      for(j=0; j<26; j++){
	*nwindows = *nwindows + (*nchrwindows)[j][0];
      }
      if(vb==0){
      cout << "\nNr bins: "<< *nwindows << endl;
      logfile << "\nNr bins: "<< *nwindows << endl;

      }
    } // intervalfile == " "


    if(intervalfile != " "){
      if(binsizeRare!=0){
	cout << "Choose ONE binning method"<<endl;
	errorfile << "Choose ONE binning method"<<endl;
	logfile << "Choose ONE binning method"<<endl;
	exit(1);
      }

      i=0; // counts chromosomes
      j=0; // counts SNPs
      l=0; // counts bins

      *nwindows=0;
      int currentwindow=0;
      int newwindow=0;
      int snpsin=0;
      int snpsout=0;

      if((maxRareInBin-minRareInBin)<0 || maxRareInBin<1 || minRareInBin<0){
	cout <<"FILTER_SMALL_BINS"<<" "<<minRareInBin<<" and SPLIT_LARGE_BINS "<<maxRareInBin<<" is not valid."<<endl;
	logfile <<"FILTER_SMALL_BINS"<<" "<<minRareInBin<<" and SPLIT_LARGE_BINS "<<maxRareInBin<<" is not valid."<<endl;
	errorfile <<"FILTER_SMALL_BINS"<<" "<<minRareInBin<<" and SPLIT_LARGE_BINS "<<maxRareInBin<<" is not valid."<<endl;
	exit(1);
      }

      for(l=0; l<nlinestped; l++){	  // Fill mapin
	if(map[l].analysis_in == 1){
	  gffchr[atoi(map[l].chr)-1].mapin=1;
	  if(l==0){
	    gffchr[atoi(map[l].chr)-1].relmapstart=l;
	    gffchr[atoi(map[l].chr)-1].relmapend=l;
	  }
	  else if(l!=0 && atoi(map[l-1].chr) != atoi(map[l].chr)){
	    gffchr[atoi(map[l-1].chr)-1].relmapend=l-1;
	    gffchr[atoi(map[l].chr)-1].relmapstart=l;
	  }
	  else if(l==nlinestped-1){
	    gffchr[atoi(map[l].chr)-1].relmapend=l;
	  }
	}
      }
      if(verbose==2){
	cout <<"\nLines in map that were placed on intervalfile's chromosomes:"<<endl;
	for(i=0; i<26; i++){
	  if(gffchr[i].mapin==1 && gffchr[i].in==1){
	    cout<<i+1<<":\t"<<gffchr[i].relmapstart+1<<"-"<< gffchr[i].relmapend+1 <<endl;
	  }
	}
	cout<<endl;
	cout <<"\nRelative coord in merged intervals that were placed on maps's chr:"<<endl;
	for(i=0; i<26; i++){
	  if(gffchr[i].mapin==1 && gffchr[i].mergein==1 && gffchr[i].relmergestart!=-1 && gffchr[i].relmergeend!=-1){
	    cout<<i+1<<":\t"<<gffchr[i].relmergestart+1<<"-"<< gffchr[i].relmergeend+1 <<endl;
	  }
	}
	cout<<endl;
      }

      int any=0;
      for(i=0; i<26; i++){
	if(gffchr[i].mapin==1 && gffchr[i].mergein==1 && gffchr[i].relmergestart!=-1 && gffchr[i].relmergeend!=-1){
	  any=1;
	}
      }
      if(any==0){
	cout<<"Intervals cover no SNPs!"<<endl;
	logfile<<"Intervals cover no SNPs!"<<endl;
	errorfile<<"Intervals cover no SNPs!"<<endl;
	logfile.close();errorfile.close();exit(1);
      }

      i=0;
      j=0;

      int binlines=0;
      int oldbinline=0;
      int knownbin=0;

      int m=0;

      l=0;

      int mapchr=0;
      int SNPsInInterval=0;
      int tmpstart=0;
      int tmpend=0;
      int QCSNPs=0;
      int QCSNPsin=0;
      int rares=0;
      int *raremap=new int[nlinestped]();

      if(!raremap)die("Memory error in raremap!");

      // Number of rare qc-SNPs
      for(j=0; j<nlinestped; j++){
	//	cout<<j<<" "<<map[j].analysis_in<<" "<<raremap[j]<<endl;
	if(map[j].analysis_in==1){
	  raremap[j]=1;
	  if(israre(j,map[j].mafr,raref,map)){
	    QCSNPs++;
	  }
	}
	else{
	  raremap[j]=-99;
	}
      }

      //      raresnps = (struct RARES *) realloc(raresnps, (QCSNPs)* sizeof(struct RARES));
      raresnps = new struct RARES[QCSNPs];
      if(!raresnps)die("Memory error in raresnps!");

      l=0;

      for(j=0; j<nlinestped; j++){
	if(map[j].analysis_in==1){
	  if(israre(j,map[j].mafr,raref,map)==1){
	    raremap[j]=l;
	    // raresnps[l].rs = NULL;
	    raresnps[l].rs = (char *) malloc((strlen(map[j].rs)+1) * sizeof(char));
	    if(!raresnps[l].rs) die("Memory allocation error in raresnps[l].rs!");
	    //	 raresnps[l].rs = (char *) realloc(raresnps[l].rs, (strlen(map[j].rs)+ 1) * sizeof(char));
	    //	 raresnps[l].rs = (char *) realloc(raresnps[l].rs, (strlen(map[j].rs)+ 1) * sizeof(char));
	    strcpy(raresnps[l].rs,map[j].rs);
	    raresnps[l].chr = atoi(map[j].chr)-1;
	    raresnps[l].pos = map[j].pos;
	    raresnps[l].in = 0;
	    l++;
	  }
	}
      }
      l=0;
      cout <<QCSNPs <<" rare QC-SNPs found."<<endl;


      int failedMin=0;
      int failedSNPs=0;
      // Intervals are projected on bins here
      for(mapchr=0; mapchr<26; mapchr++){
	if(gffchr[mapchr].mapin==1){
	  // Go through all intervals on that chromsosome
	  for(i=gffchr[mapchr].relmergestart; i<gffchr[mapchr].relmergeend+1; i++){
	    SNPsInInterval=0;
	    // Chromosome-wise, go through map; j = line in map
	    for(j=gffchr[mapchr].relmapstart; j<gffchr[mapchr].relmapend+1; j++){
	      if(israre(j,map[j].mafr,raref,map)==1){
		if(map[j].pos>=gffmerge[i].start &&  map[j].pos<=gffmerge[i].end){
		  if(SNPsInInterval==0){
		    tmpstart=j;
		  }
		  raresnps[raremap[j]].in=1;
		  SNPsInInterval++;
		  tmpend=j;
		}
	      }
	    }
	    // INTERVAL QC-PASSED
	    if(SNPsInInterval>=minRareInBin && SNPsInInterval<=maxRareInBin){
	      intervals = (struct INTERVALS *) realloc(intervals, (binlines+1)* sizeof(struct INTERVALS));
	      if(!intervals) die("Memory allocation error in intervals!");
	      intervals[binlines].n=SNPsInInterval;
	      if(featurecol>=0){
		intervals[binlines].feature = NULL;
		intervals[binlines].feature = (char *) realloc(intervals[binlines].feature, (strlen(gffmerge[i].feature)+ 1) * sizeof(char));
		if(!intervals[binlines].feature) die("Memory allocation error in intervals[binlines].feature!");
		strcpy(intervals[binlines].feature, gffmerge[i].feature);
	      }
	      intervals[binlines].chr=gffmerge[i].chr;
	      intervals[binlines].start=tmpstart;
	      intervals[binlines].end=tmpend;
	      binlines++;
	    }
	    // MAX_RARE_IN_BIN
	    else if(SNPsInInterval>maxRareInBin){
	      int split=0;
	      int part=0;
	      int RareInInterval=0;
	      int RareInIntervalLast=0;
	      split=SNPsInInterval/maxRareInBin;
	      if(SNPsInInterval%maxRareInBin!=0){
		split++;
	      }
	      RareInInterval=SNPsInInterval/split;
	      if(SNPsInInterval%split!=0){
		RareInInterval++;
	      }
	      if(SNPsInInterval%RareInInterval!=0){
		RareInIntervalLast=SNPsInInterval%RareInInterval;
	      }
	      else{
		RareInIntervalLast=RareInInterval;
	      }
	      if(verbose==2){
		cout << SNPsInInterval<< " SNPs in bin " << binlines+1 << " from SNP " << tmpstart << " to SNP " << tmpend<<"." <<endl;
		cout <<"BIN split in "<<split << " parts, "<<RareInInterval<<" rare SNPs long, last bin has "<<RareInIntervalLast<<" rare SNPs.\n"<<endl;
	      }
	      for(l=0; l<split; l++){
		intervals = (struct INTERVALS *) realloc(intervals, (binlines+1)* sizeof(struct INTERVALS));
		if(!intervals) die("Memory allocation error in intervals!");
		if(l==0){
		  intervals[binlines].n=RareInInterval;

		  // Find net stmpstart,end
		  int candidate=0;
		  int candidatefirst=0;
		  for(m=tmpstart; candidate!=RareInInterval; m++){
		    if(israre(m,map[m].mafr,raref,map)==1){
		      if(candidate==0){
			candidatefirst=m;
		      }
		      candidate++;
		    }
		  }
		  intervals[binlines].start=candidatefirst;

		  intervals[binlines].end=m-1;
		  tmpend=m-1;
		  tmpstart=candidatefirst;
		}
		else if(l!=split-1){
		  intervals[binlines].n=RareInInterval;
		  intervals[binlines].start=intervals[binlines-1].end+1;
		  // Find net stmpstart, end
		  int candidate=0;
		  int candidatefirst=0;
		  for(m=tmpend+1; candidate!=RareInInterval; m++){
		    if(israre(m,map[m].mafr,raref,map)==1){
		      if(candidate==0){
			candidatefirst=m;
		      }
		      candidate++;
		    }
		  }
		  intervals[binlines].start=candidatefirst;
		  intervals[binlines].end=m-1;
		  tmpend=m-1;
		  tmpstart=candidatefirst;
		}
		else if(l==split-1){
		  intervals[binlines].n=RareInIntervalLast;
		  // Find net stmpstart, end
		  int candidate=0;
		  int candidatefirst=0;
		  for(m=tmpend+1; candidate!=RareInIntervalLast; m++){
		    if(israre(m,map[m].mafr,raref,map)==1){
		      if(candidate==0){
			candidatefirst=m;
		      }
		      candidate++;
		    }
		  }
		  intervals[binlines].start=candidatefirst;
		  intervals[binlines].end=m-1;
		  tmpend=m-1;
		  tmpstart=candidatefirst;
		}
		if(featurecol>=0){
		  intervals[binlines].feature = NULL;
		  intervals[binlines].feature = (char *) realloc(intervals[binlines].feature, (strlen(gffmerge[i].feature)+ 1+15) * sizeof(char));
		  if(!intervals[binlines].feature) die("Memory allocation error in intervals[binlines].feature!");
		  sprintf( intervals[binlines].feature,"%s(Part%i/%i)",gffmerge[i].feature,l+1,split);
		}
		if (!intervals){
		  cout << "Memory allocation error in intervals"<<endl;
		  exit(1);
		}
		intervals[binlines].chr=gffmerge[i].chr;
		binlines++;
	      }
	    }
	    else{
	      //	    binlines--;
	      failedMin++;
	      failedSNPs=failedSNPs+SNPsInInterval;
	    }
	    SNPsInInterval=0;
	  }
	}
      }
      if(failedMin){
	cout <<failedMin <<" bins were filtered because they were too small: "<<failedSNPs<< " eliminated SNPs."<<endl;
	logfile <<failedMin <<" bins were filtered because they were too small: "<<failedSNPs<< " eliminated SNPs."<<endl;
      }
      QCSNPsin=0;
      for(l=0;l<QCSNPs;l++){
	if(raresnps[l].in==1){
	  QCSNPsin++;
	}
      }
      delete[] raremap;

      for(l=0;l<QCSNPs;l++){
	//	free(raresnps[l].rs);
	//	delete[] raresnps[l].rs;
      }
      //      free(raresnps);
      delete[] raresnps;

      if(QCSNPsin==0){
	cout<<"No SNPs could be assigned to bins!"<<endl;
	errorfile<<"No SNPs could be assigned to bins!"<<endl;
	logfile<<"No SNPs could be assigned to bins!"<<endl;
	exit(1);
      }

      *nwindows=binlines;

      cout << binlines <<" bins assigned."<<endl;
      cout <<QCSNPsin<<" rare SNPs inside, " << QCSNPs-QCSNPsin<<" outside of bins."<<endl;

      int s1= 0;
      int s2= 0;

      if(binamax!=0){
	cout<<"DENSE_BINNING"<<endl;
	gffmerge = (struct GFFMERGE *) realloc(gffmerge, (binlines+1)* sizeof(struct GFFMERGE));
	seeds = (struct INTERVALS *) realloc(seeds, (binlines+1)* sizeof(struct INTERVALS));
	for(i=0;i<binlines;i++){

	  if(featurecol>=0){
	    gffmerge[i].feature=NULL;
	    gffmerge[i].feature = (char *) realloc(gffmerge[i].feature, (strlen(intervals[i].feature)+ 1) * sizeof(char));
	    gffmerge[i].feature=intervals[i].feature;
	  }
	  gffmerge[i].chr=intervals[i].chr;
	  gffmerge[i].start=intervals[i].start;
	  gffmerge[i].end=intervals[i].end;
	}

	int nvariations[binlines];
	int adjmin[binlines];
	int adjmax[binlines];
	int adjlines=0;

	for(i=0;i<binlines;i++){
	  if(binamin>5000){
	    adjmin[i]=binamin-10000+intervals[i].n;
	  }
	  else{
	    adjmin[i]=binamin;
	  }
	  if(binamax>5000){
	    adjmax[i]=binamax-10000+intervals[i].n;
	  }
	  else{
	    adjmax[i]=binamax;
	  }
	  if(adjmax[i]>intervals[i].n){
	    adjmax[i]=intervals[i].n;
	  }

	  nvariations[i]=(adjmax[i]-adjmin[i]+1)*(2*intervals[i].n-adjmin[i]-adjmax[i]+2)/2;
	  gffmerge[i].chr=intervals[i].chr;
	  // Find steps to rare
	  gffmerge[i].start=intervals[i].start;
	  gffmerge[i].end=intervals[i].end;

	  seeds[i].chr=intervals[i].chr;
	  seeds[i].start=intervals[i].start;
	  seeds[i].end=intervals[i].end;
	  seeds[i].n=intervals[i].n;
	  seeds[i].N=nvariations[i];

	  cout <<"Creating " << nvariations[i] <<" new bins out of interval "<<i+1 <<"."<<endl;
	  adjlines=adjlines+nvariations[i];
	}

	cout <<"Overall, " << adjlines <<" new bins have been created.\n"<<endl;

	// BINSIZE and SHIFT
	int bs=0;
	int sh=0;
	int s=0;
	int start=0;
	int end=0;
	int tstart=0;
	int tend=0;

	intervals2 = (struct INTERVALS*) realloc(intervals2, (adjlines+1)* sizeof(struct INTERVALS));
	if(!intervals2) die("Memory allocation error in intervals2!");

	int adjbinlines=0;
	for(i=0;i<binlines;i++){
	  for(bs=adjmin[i]; bs<=adjmax[i]; bs++){
	    for(sh=0; sh<intervals[i].n-bs+1; sh++){
	      start=intervals[i].start+sh;
	      end=intervals[i].start+bs+sh-1;
	      s=0;
	      tstart=intervals[i].start;
	      while(tstart<start){
		s++;
		if(israre(gffmerge[i].start+s,map[gffmerge[i].start+s].mafr, raref, map)){
		  tstart++;
		}
	      }
	      intervals2[adjbinlines].start=gffmerge[i].start+s;
	      s=0;
	      tend=0;

	      while(tend<bs-1){
		s++;
		if(israre(intervals2[adjbinlines].start+s,map[intervals2[adjbinlines].start+s].mafr, raref, map)){
		  tend++;
		}
	      }

	      intervals2[adjbinlines].end=intervals2[adjbinlines].start+s;
	      if(featurecol>=0){
		intervals2[adjbinlines].feature=NULL;
		intervals2[adjbinlines].feature = (char *) realloc(intervals[adjbinlines].feature, (strlen(gffmerge[i].feature)+ 1) * sizeof(char));
		if(!intervals2[adjbinlines].feature) die("Memory allocation error in intervals2[adjbinlines].feature!");
		intervals2[adjbinlines].feature = gffmerge[i].feature;
	      }
	      intervals2[adjbinlines].chr=intervals[i].chr;
	      intervals2[adjbinlines].s1=bs;
	      intervals2[adjbinlines].s2=sh;
	      adjbinlines++;
	    }
	  }

	  if(adjbinlines!=nvariations[0]){
	    cout<<"adjbinlines!=nvariations"<<endl;
	  }
	}

	// Copy back to interval
	intervals = (struct INTERVALS*) realloc(intervals, (adjlines+1)* sizeof(struct INTERVALS));
	if(!intervals) die("Memory allocation error in intervals!");
	for(i=0;i<adjlines;i++){
	  intervals[i].seed=adjlines+1;
	  intervals[i].start=intervals2[i].start;
	  intervals[i].end=intervals2[i].end;

	  if(featurecol>=0){
	    intervals[i].feature=NULL;
	    intervals[i].feature = (char *) realloc(intervals[i].feature, (strlen(intervals2[i].feature)+ 1) * sizeof(char));
	    if(!intervals[i].feature) die("Memory allocation error in intervals[i].feature!");
	    intervals[i].feature = intervals2[i].feature;
	  }
	  intervals[i].chr=intervals2[i].chr;
	  intervals[i].s1=intervals2[i].s1;
	  intervals[i].s2=intervals2[i].s2;
	}
	binlines=adjlines;
	*nwindows=binlines;
      }
    } // INTERVALFILE != " "

  }


void limitspositionsRARE(bool optimalrare, int n, int bin, int binsizeRare, int nlinestped, int nlinestfam, double raref, int nsim, struct MAP *map, struct COUNTS **counts, int *nRareSNPs, int **windowPositions, int *nSNPsInWindow, double **rareLimits, int **rareLimitsNCT, int **rareLimitsNCTInverse, int *nRareLimits, int **nchrwindows, int nwindows, int nwindowssinglebin, int thread, fstream &errorfile, fstream &logfile, string intervalfile, int binamin, int binamax, string SetIDfile, int intervaleditor, int setid, struct WINDOW *window, struct WINDOW_INTER *window_inter, int verbose, int **doublewindowcoord, int collinter, int nwords, int* SNPMapInverse, uint64_t*** BinSNPs, int vb, int NCT, int *nCarriers){
  if(nsim==0 && (FISHERflag!=0 ||CMATflag!=0 || (COLLflag!=0 && !vb && optimalrare==1 )|| (REGRESSIONflag!=0 && optimalrare==1 ) || (FRACREGflag!=0 && optimalrare==1 ) || (COLLREGflag!=0 && optimalrare==1 ))){
    cout<<"Rare variant analysis cannot be conducted with zero Monte-Carlo simulations! Only COLL, REG, FRACREG and COLLREG with fixed threshold MAF."<<endl;
    logfile<<"Rare variant analysis cannot be conducted with zero Monte-Carlo simulations! Only COLL, REG, FRACREG and COLLREG with fixed threshold MAF."<<endl;
    errorfile<<"Rare variant analysis cannot be conducted with zero Monte-Carlo simulations! Only COLL, REG, FRACREG and COLLREG with fixed threshold MAF."<<endl;
    exit(1);
  }

  int i, j, k, l, m, start, end;

  if(!vb){
    if(intervalfile == " "){
      if(binsizeRare!=0){
	start=0;
	end=0;
	j=0; // count chromosomes
	k=0; // count windows on chromosome
	m=0; // count global windows
	l=0; // count RareSNPs until window full
	for(i=0; i<nlinestped; i++){
	  while(map[i].pos==0 || atoi(map[i].chr)==0){
	    i++;
	  }
	  j = atoi(map[i].chr)-1;
	  if(((map[i].mafr-raref)<EPS || fabs(map[i].mafr-raref)<EPS) && map[i].mafr>EPS && map[i].analysis_in == 1){
	    l++;
	    if(l==1){
	      windowPositions[m][0]=i;
	      //k++;
	    }
	    if(l==binsizeRare){
	      windowPositions[m][1]=i;
	      nRareSNPs[m]=l;
	      l=0;
	      k++;
	      if(nchrwindows[j][0]==k){
		k=0;
	      }
	      m++;
	    }
	    else if(l==nchrwindows[j][1] && k==nchrwindows[j][0]-1){
	      windowPositions[m][1]=i;
	      nRareSNPs[m]=l;
	      l=0;
	      k=0;
	      m++;
	    }
	  }
	}
	l=0; // count SNPs
	for(m=0; m<nwindowssinglebin; m++){
	  start = windowPositions[m][0];
	  end = windowPositions[m][1];
	  if(verbose==2){
	    for(k=start; k<end+1; k++){
	      nSNPsInWindow[m]++;
	    }
	  }
	}
	cout << "\n";
      } // binsizeRare
    } // INTERVALFILE =" "
    else if(intervalfile != " "){
      for(i=0; i<nwindowssinglebin; i++){
	windowPositions[i][0]=intervals[i].start;
	windowPositions[i][1]=intervals[i].end;
      }
      for(i=0; i<nwindowssinglebin; i++){
	for(k=windowPositions[i][0]; k<windowPositions[i][1]+1; k++){
	  if(map[k].analysis_in==1){
	    if(verbose==2){
	      nSNPsInWindow[i]++;
	    }
	    if(israre(k,map[k].mafr, raref, map)==1){
	      nRareSNPs[i]++;
	    }
	  }
	}
      }
    }

    l=0; // count SNPs
    m=0; // count windows
    // Count number of DIFFERENT MAFs per window
    // Create variable-sized arrays containing DIFFERENT rare frequencies

    if(collinter!=3 && collinter!=4){
      for(m=0; m<nwindows; m++){
	if(nRareSNPs[m]==0){
	  window[m].ntests=0;
	}
      }
    }
    cout<<"Done"<<endl;


    for(m=0; m<nwindowssinglebin; m++){
      k=0; // Position of rare SNP in window
      start = windowPositions[m][0];
      end = windowPositions[m][1];


      double *mafcounter=NULL;
      mafcounter=new double[nRareSNPs[m]+1]; // nRareSNPs[m]+1 ??
      if(!mafcounter) die("Memory allocation error in mafcounter!");
      nRareLimits[m]=0;

      for(l=start; l<end+1; l++){
	if(((map[l].mafr-raref)< EPS || (fabs(map[l].mafr-raref)<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
	  int duplicate=0;
	  for(int hk=0; hk<nRareLimits[m]; hk++){
	    if(fabs(map[l].mafr-mafcounter[hk]) < EPS) duplicate=1;
	  }
	  if(duplicate==0){
	    mafcounter[nRareLimits[m]]=map[l].mafr;
	    nRareLimits[m]++;
	  }
	}
      }

      rareLimits[m]=new double[nRareLimits[m]+1];
      if(!rareLimits[m]) die("Memory allocation error in rareLimits[m]!");
      window[m].n_level=nRareLimits[m];
      if(optimalrare==1 && nRareSNPs[m]){
	window[m].n_at_level=new int[nRareLimits[m]+1];
	if(!window[m].n_at_level) die("Memory allocation error in window[m].n_at_level!");
	window[m].maf_at_level=new double[nRareLimits[m]+1];
	if(!window[m].maf_at_level) die("Memory allocation error in window[m].maf_at_level!");
	window[m].levelpos=new int*[nRareLimits[m]+1];
	if(!window[m].levelpos) die("Memory allocation error in window[m].levelpos!");
      }
      /*      else if(optimalrare==0 && nRareSNPs[m]){
	window[m].maf_at_level=new double[1];
	if(!window[m].maf_at_level) die("Memory allocation error in window[m].maf_at_level!");
	window[m].n_at_level=new int[1];
	if(!window[m].n_at_level) die("Memory allocation error in window[m].n_at_level!");
	window[m].levelpos=new int*[1];
	if(!window[m].levelpos) die("Memory allocation error in window[m].levelpos!");
      }
      */

      for(int fu=0; fu<nRareLimits[m]+1; fu++){
	rareLimits[m][fu]=mafcounter[fu];
      }
      qsort(rareLimits[m], nRareLimits[m], sizeof(double), comparedouble);


      if(optimalrare==0 && nRareSNPs[m]){
	//    rareLimits[m]=new double[1];
	window[m].n_level=0;
	window[m].maf_at_level=new double[1];
	if(!window[m].maf_at_level) die("Memory allocation error in window[m].maf_at_level!");
	window[m].n_at_level=new int[1];
	if(!window[m].n_at_level) die("Memory allocation error in window[m].n_at_level!");
	window[m].levelpos=new int*[1];
	if(!window[m].levelpos) die("Memory allocation error in window[m].levelpos!");
	window[m].levelpos[0]=new int[nRareSNPs[m]];
	//	rareLimits[m][0]=rareLimits[nRareLimits[m]-1];
	window[m].maf_at_level[0]=rareLimits[m][nRareLimits[m]-1];
	window[m].n_at_level[0]=nRareSNPs[m];
	int count=0;
	for(l=start; l<end+1; l++){
	  if(((map[l].mafr-window[m].maf_at_level[0])< EPS || (fabs(map[l].mafr-window[m].maf_at_level[0])<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
	    window[m].levelpos[0][count]=l;
	    count++;
	  }
	}
	//	rareLimits[m]=new double [1];
	//	if(!rareLimits[m]) die("Memory allocation error in rareLimits[m]!");
	rareLimits[m][0]=window[m].maf_at_level[0];
	nRareLimits[m]=1;
      }
      else if(optimalrare==1 && nRareSNPs[m]){
	for(int lev=0; lev<nRareLimits[m]; lev++){ // sorted ones
	  window[m].maf_at_level[lev]=rareLimits[m][lev];
	  window[m].n_at_level[lev]=0;
	  for(l=start; l<end+1; l++){
	    if(((map[l].mafr-window[m].maf_at_level[lev])< EPS || (fabs(map[l].mafr-window[m].maf_at_level[lev])<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
	      window[m].n_at_level[lev]+=1;
	    }
	  }
	  window[m].levelpos[lev]=new int[window[m].n_at_level[lev]];
	  if(!window[m].levelpos[lev]) die("Memory allocation error in window[m].levelpos[lev]!");
	}

	for(int lev=0; lev<nRareLimits[m]; lev++){ //Find levelpositions
	  int count=0;
	  for(l=start; l<end+1; l++){
	    if(((map[l].mafr-window[m].maf_at_level[lev])< EPS || (fabs(map[l].mafr-window[m].maf_at_level[lev])<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
	      window[m].levelpos[lev][count]=l;
	      count++;
	    }
	  }
	}
      }
      delete[] mafcounter;
    }
    cout <<"Variants at various MAF levels for rare variant analysis have been identified...\n";

    if(collinter==3 || collinter==4){
      int m1,m2,start1,end1,start2,end2;
      m=nwindowssinglebin-1;
      // for RAREINTER: compute double windows
      for(m1=0; m1<nwindowssinglebin; m1++){
	for(m2=m1+1; m2<nwindowssinglebin; m2++){
	  m++;
	  doublewindowcoord[m][0]=m1;
	  doublewindowcoord[m][1]=m2;
	}
      }
    }
    else if(collinter==1){
      int m1,m2,start1,end1,start2,end2;
      m=nwindowssinglebin-1;
      // for RAREINTER: compute double windows
      for(m1=0; m1<nwindowssinglebin; m1++){
	for(m2=m1+1; m2<nwindowssinglebin; m2++){
	  m++;
	  doublewindowcoord[m][0]=m1;
	  doublewindowcoord[m][1]=m2;
	  k=0; // Position of rare SNP in window
	  start1 = windowPositions[m1][0];
	  end1 = windowPositions[m1][1];
	  start2 = windowPositions[m2][0];
	  end2 = windowPositions[m2][1];

	  nRareSNPs[m]=nRareSNPs[m1]+nRareSNPs[m2];

	  double *mafcounter=NULL;
	  mafcounter=new double[nRareSNPs[m]+1];
	  if(!mafcounter) die("Memory allocation error in mafcounter!");
	  nRareLimits[m]=0;

	  for(l=start1; l<end1+1; l++){
	    if(((map[l].mafr-raref)< EPS || (fabs(map[l].mafr-raref)<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
	      int duplicate=0;
	      for(int hk=0; hk<nRareLimits[m]; hk++){
		if(fabs(map[l].mafr-mafcounter[hk]) < EPS) duplicate=1;
	      }
	      if(!duplicate){
		mafcounter[nRareLimits[m]]=map[l].mafr;
		nRareLimits[m]++;
	      }
	    }
	  }
	  for(l=start2; l<end2+1; l++){
	    if(((map[l].mafr-raref)< EPS || (fabs(map[l].mafr-raref)<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
	      int duplicate=0;
	      for(int hk=0; hk<nRareLimits[m]; hk++){
		if(fabs(map[l].mafr-mafcounter[hk]) < EPS) duplicate=1;
	      }
	      if(duplicate==0){
		mafcounter[nRareLimits[m]]=map[l].mafr;
		nRareLimits[m]++;
	      }
	    }
	  }
	  rareLimits[m]=new double[nRareLimits[m]+1];
	  if(!rareLimits[m]) die("Memory allocation error in rareLimits[m]!");
	  window[m].n_level=nRareLimits[m];

	  if(optimalrare==1 && nRareSNPs[m]){
	    window[m].n_at_level=new int[nRareLimits[m]+1];
	    if(!window[m].n_at_level) die("Memory allocation error in window[m].n_at_level!");
	    window[m].maf_at_level=new double[nRareLimits[m]+1];
	    if(!window[m].maf_at_level) die("Memory allocation error in window[m].maf_at_level!");
	    window[m].levelpos=new int*[nRareLimits[m]+1];
	    if(!window[m].levelpos) die("Memory allocation error in window[m].levelpos!");
	  }
	  else if(optimalrare==0 && nRareSNPs[m]){
	    window[m].maf_at_level=new double[1];
	    if(!window[m].maf_at_level) die("Memory allocation error in window[m].maf_at_level!");
	    window[m].n_at_level=new int[1];
	    if(!window[m].n_at_level) die("Memory allocation error in window[m].n_at_level!");
	    window[m].levelpos=new int*[1];
	    if(!window[m].levelpos) die("Memory allocation error in window[m].levelpos!");
	  }
	  for(int fu=0; fu<nRareLimits[m]+1; fu++){
	    rareLimits[m][fu]=mafcounter[fu];
	  }
	  qsort(rareLimits[m], nRareLimits[m], sizeof(double), comparedouble);


	  if(optimalrare==0 && nRareSNPs[m]){
	    //    rareLimits[m]=new double[1];
	    window[m].n_level=0;
	    window[m].maf_at_level=new double[1];
	    if(!window[m].maf_at_level) die("Memory allocation error in window[m].maf_at_level!");
	    window[m].n_at_level=new int[1];
	    if(!window[m].n_at_level) die("Memory allocation error in window[m].n_at_level!");
	    window[m].levelpos=new int*[1];
	    if(!window[m].levelpos) die("Memory allocation error in window[m].levelpos!");
	    window[m].levelpos[0]=new int[nRareSNPs[m]];
	    if(!window[m].levelpos[0]) die("Memory allocation error in window[m].levelpos[0]!");
	    //	rareLimits[m][0]=rareLimits[nRareLimits[m]-1];
	    window[m].maf_at_level[0]=rareLimits[m][nRareLimits[m]-1];
	    window[m].n_at_level[0]=nRareSNPs[m];
	    int count=0;
	    for(l=start1; l<end1+1; l++){
	      if(((map[l].mafr-window[m].maf_at_level[0])< EPS || (fabs(map[l].mafr-window[m].maf_at_level[0])<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
		window[m].levelpos[0][count]=l;
		count++;
	      }
	    }
	    for(l=start2; l<end2+1; l++){
	      if(((map[l].mafr-window[m].maf_at_level[0])< EPS || (fabs(map[l].mafr-window[m].maf_at_level[0])<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
		window[m].levelpos[0][count]=l;
		count++;
	      }
	    }
	    rareLimits[m]=new double [1];
	    if(!rareLimits[m]) die("Memory allocation error in rareLimits[m]!");
	    rareLimits[m][0]=window[m].maf_at_level[0];
	    nRareLimits[m]=1;
	  }
	  else if(optimalrare==1 && nRareSNPs[m]){
	    for(int lev=0; lev<nRareLimits[m]; lev++){ // sorted ones
	      window[m].maf_at_level[lev]=rareLimits[m][lev];
	      window[m].n_at_level[lev]=0;
	      for(l=start1; l<end1+1; l++){
		if(((map[l].mafr-window[m].maf_at_level[lev])< EPS || (fabs(map[l].mafr-window[m].maf_at_level[lev])<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
		  window[m].n_at_level[lev]+=1;
		}
	      }
	      for(l=start2; l<end2+1; l++){
		if(((map[l].mafr-window[m].maf_at_level[lev])< EPS || (fabs(map[l].mafr-window[m].maf_at_level[lev])<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
		  window[m].n_at_level[lev]+=1;
		}
	      }
	      window[m].levelpos[lev]=new int[window[m].n_at_level[lev]];
	      if(!window[m].levelpos[lev]) die("Memory allocation error in window[m].levelpos[lev]!");
	    }

	    for(int lev=0; lev<nRareLimits[m]; lev++){ //Find levelpositions
	      int count=0;
	      for(l=start1; l<end1+1; l++){
		if(((map[l].mafr-window[m].maf_at_level[lev])< EPS || (fabs(map[l].mafr-window[m].maf_at_level[lev])<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
		  window[m].levelpos[lev][count]=l;
		  count++;
		}
	      }
	      for(l=start2; l<end2+1; l++){
		if(((map[l].mafr-window[m].maf_at_level[lev])< EPS || (fabs(map[l].mafr-window[m].maf_at_level[lev])<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
		  window[m].levelpos[lev][count]=l;
		  count++;
		}
	      }
	    }
	  }
	  delete[] mafcounter;
	}
      }
      cout <<"Variants at various MAF levels for rare interaction analysis have been identified...\n";
    }

    if(collinter==2){
      int m1,m2,start1,end1,start2,end2,x;

      //    m=nwindowssinglebin-1;
      // for RAREINTER: compute double windows with COLL_INTER 2
      for(m=nwindowssinglebin; m<nwindows; m++){
	int nRareLimitsDoubleCarriers=0;
	m1=doublewindowcoord[m][0];
	m2=doublewindowcoord[m][1];
	k=0; // Position of rare SNP in window
	start1 = windowPositions[m1][0];
	end1 = windowPositions[m1][1];
	start2 = windowPositions[m2][0];
	end2 = windowPositions[m2][1];
	//	window_inter[m].start1=start1;
	//	window_inter[m].end1=end1;
	//	window_inter[m].start2=start2;
	//	window_inter[m].end2=end2;
	//	nRareSNPs[m]=nRareSNPs[m1]+nRareSNPs[m2];
	//	double *mafcounter=NULL;
	//      mafcounter=new double[nRareSNPs[m]+1];
	//      nRareLimits[m]=0;

	// find the larger of both smallest MAF levels in the two subbins
	double combinedlowestlevel=0;
	int ncombinedlowestlevel=0;
	if(abs(rareLimits[m1][0]-rareLimits[m2][0])<EPS){ // equal lower levels
	  combinedlowestlevel=rareLimits[m1][0];
	}
	else if((rareLimits[m1][0]-rareLimits[m2][0])>EPS){ // first lower level is larger
	  combinedlowestlevel=rareLimits[m1][0];
	}
	else if((rareLimits[m2][0]-rareLimits[m1][0])>EPS){ // second lower level is larger
	  combinedlowestlevel=rareLimits[m2][0];
	}
	// find the level number in combined bin:
	for(x=0; x<nRareLimits[m]; x++){
	  if(abs(rareLimits[m][x]-combinedlowestlevel)<EPS){
	    ncombinedlowestlevel=x;
	    //	    cout<<"ncombinedlowestlevel "<<ncombinedlowestlevel<<" of "<< nRareLimits[m] <<" m "<<m<<endl;
	    break;
	  }
	}

	// check if some carrierers are present in both bins at lowest level, otherwise find new lowest combined level
	x=ncombinedlowestlevel;
	uint64_t* dummy = new uint64_t[nwords];
	if(!dummy) die("Memory allocation error in dummy!");
	uint64_t* dummy1 = new uint64_t[nwords];
	if(!dummy1) die("Memory allocation error in dummy1!");
	uint64_t* dummy2 = new uint64_t[nwords];
	if(!dummy2) die("Memory allocation error in dummy2!");
	int doubleCarriers=0;
	int doubleCarriersPrev=0;
	double *rareLimitsDoubleCarriers=NULL;
	rareLimitsDoubleCarriers=new double[nRareLimits[m]]; // contains the number of rareLimits or less (!)
	if(!rareLimitsDoubleCarriers) die("Memory allocation error in rareLimitsDoubleCarriers!");

	nRareLimitsDoubleCarriers=0;
	while(x<nRareLimits[m]){
	  doubleCarriers=0;
	  memset(dummy, 0, nwords*sizeof(uint64_t));
	  memset(dummy1, 0, nwords*sizeof(uint64_t));
	  memset(dummy2, 0, nwords*sizeof(uint64_t));
	  for(int j=0; j<window[m].n_at_level[x]; j++){
	    int i=window[m].levelpos[x][j];
	    int iMod = SNPMapInverse[i];
	    if (iMod==-1) continue;
	    if(i>=start1 && i<=end1){
	      for (int p=0; p<nwords; p++)
		dummy1[p] |= (BinSNPs[iMod][p][1] | BinSNPs[iMod][p][2]);
	    }
	    if(i>=start2 && i<=end2){
	      for (int p=0; p<nwords; p++)
		dummy2[p] |= (BinSNPs[iMod][p][1] | BinSNPs[iMod][p][2]);
	    }
	  }
	  for(int p=0; p<nwords; p++){
	    dummy[p] |= dummy1[p] & dummy2[p];
	    doubleCarriers += bitcount64(dummy[p]);
	  }
	  if(doubleCarriers==0){ // check if level contains double carriers at all, otherwise raise lowest level
	    ncombinedlowestlevel++;
	  }
	  if(doubleCarriers>doubleCarriersPrev){
	    //	    cout <<"m "<<m<< " nRareLimitsDoubleCarriers "<< nRareLimitsDoubleCarriers<<"\n" <<flush;
	    nRareLimitsDoubleCarriers++;
	    doubleCarriersPrev=doubleCarriers;
	    rareLimitsDoubleCarriers[nRareLimitsDoubleCarriers-1]=rareLimits[m][x]; // Fill in new rareLimits
	  }
	  x++;
	}
	delete[] dummy;
	delete[] dummy1;
	delete[] dummy2;
	delete[] rareLimits[m];
	rareLimits[m]=new double[nRareLimitsDoubleCarriers];
	if(!rareLimits[m]) die("Memory allocation error in rareLimits[m]!");
	if(doubleCarriers==0){      // check if highest level contains carrieres at all, otherwise exclude bin from analysis
	  nRareSNPs[m]=0;
	  nRareLimits[m]=0;
	  nRareLimitsDoubleCarriers=0;
	}
	else{ // Fill in new rareLimits, nRareLimits
	  for(x=0; x<nRareLimitsDoubleCarriers; x++){
	    rareLimits[m][x]=rareLimitsDoubleCarriers[x];
	  }
	}
	if(rareLimitsDoubleCarriers!=NULL)delete[] rareLimitsDoubleCarriers; rareLimitsDoubleCarriers=NULL;

	nRareLimits[m]=nRareLimitsDoubleCarriers;
	if(nRareLimits[m]==0){
	  window[m].ntests=0;
	  nRareSNPs[m]=0;
	}
	// fill in: windo[m].n_at_level[level] window[m].levelpos[level][count]
	if(optimalrare==0 && nRareSNPs[m]){
	  //    rareLimits[m]=new double[1];
	  window[m].n_level=0;
	  window[m].maf_at_level=new double[1];
	  if(!window[m].maf_at_level) die("Memory allocation error in window[m].maf_at_level!");
	  window[m].n_at_level=new int[1];
	  if(!window[m].n_at_level) die("Memory allocation error in window[m].n_at_level!");
	  window[m].levelpos=new int*[1];
	  if(!window[m].levelpos) die("Memory allocation error in window[m].levelpos!");
	  window[m].levelpos[0]=new int[nRareSNPs[m]];
	  if(!window[m].levelpos[0]) die("Memory allocation error in window[m].levelpos[0]!");
	  //	rareLimits[m][0]=rareLimits[nRareLimits[m]-1];
	  window[m].maf_at_level[0]=rareLimits[m][nRareLimits[m]-1];
	  window[m].n_at_level[0]=nRareSNPs[m];
	  int count=0;
	  for(l=start1; l<end1+1; l++){
	    if(((map[l].mafr-window[m].maf_at_level[0])< EPS || (fabs(map[l].mafr-window[m].maf_at_level[0])<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
	      window[m].levelpos[0][count]=l;
	      count++;
	    }
	  }
	  for(l=start2; l<end2+1; l++){
	    if(((map[l].mafr-window[m].maf_at_level[0])< EPS || (fabs(map[l].mafr-window[m].maf_at_level[0])<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
	      window[m].levelpos[0][count]=l;
	      count++;
	    }
	  }
	  rareLimits[m]=new double [1];
	  if(!rareLimits[m]) die("Memory allocation error in rareLimits[m]!");
	  rareLimits[m][0]=window[m].maf_at_level[0];
	  nRareLimits[m]=1;
	}
	else if(optimalrare==1 && nRareSNPs[m]){
	  for(int lev=0; lev<nRareLimits[m]; lev++){ // sorted ones
	    window[m].maf_at_level[lev]=rareLimits[m][lev];
	    window[m].n_at_level[lev]=0;
	    for(l=start1; l<end1+1; l++){
	      if(((map[l].mafr-window[m].maf_at_level[lev])< EPS || (fabs(map[l].mafr-window[m].maf_at_level[lev])<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
		window[m].n_at_level[lev]+=1;
	      }
	    }
	    for(l=start2; l<end2+1; l++){
	      if(((map[l].mafr-window[m].maf_at_level[lev])< EPS || (fabs(map[l].mafr-window[m].maf_at_level[lev])<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
		window[m].n_at_level[lev]+=1;
	      }
	    }
	    window[m].levelpos[lev]=new int[window[m].n_at_level[lev]];
	    if(!window[m].levelpos[lev]) die("Memory allocation error in window[m].levelpos[lev]!");
	  }
	  for(int lev=0; lev<nRareLimits[m]; lev++){ //Find levelpositions
	    int count=0;
	    for(l=start1; l<end1+1; l++){
	      if(((map[l].mafr-window[m].maf_at_level[lev])< EPS || (fabs(map[l].mafr-window[m].maf_at_level[lev])<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
		window[m].levelpos[lev][count]=l;
		count++;
	      }
	    }
	    for(l=start2; l<end2+1; l++){
	      if(((map[l].mafr-window[m].maf_at_level[lev])< EPS || (fabs(map[l].mafr-window[m].maf_at_level[lev])<EPS)) && map[l].mafr>EPS && map[l].analysis_in==1 && map[l].pos>0 && atoi(map[l].chr)>0){
		window[m].levelpos[lev][count]=l;
		count++;
	      }
	    }
	  }
	}
      }
      cout <<"Variants at various MAF levels for rare interaction analysis have been identified...\n";
    }
  }
  else if(vb){
    int currentwindow=0;
    int currentchr=0;
    int oldchr=0;
    int **NCTcounter =NULL;
    NCTcounter = new int*[nwindows];
    if(!NCTcounter) die("Memory allocation error in NCTcounter!");
    for(i=0; i<nwindows; i++){
      NCTcounter[i] = new int[NCT]();
      if(!NCTcounter[i]) die("Memory allocation error in NCTcounter[i]!");
    }
    for(i=0; i<nlinestped; i++){ // Fill nCarriers
      //      nCarriers[i]=0;
      currentchr=atoi(map[i].chr);
      if(map[i].analysis_in){
	if(oldchr==0){
	  oldchr=atoi(map[i].chr);
	}
	else{
	  if(oldchr!=currentchr){
	    oldchr=currentchr;
	    currentwindow++;
	  }
	}
	// uint64_t* dummy = new uint64_t[nwords];
	// memset(dummy, 0, nwords*sizeof(uint64_t));
	// int iMod = SNPMapInverse[i];
	// if (iMod==-1) continue;
	// for (int p=0; p<nwords; p++){
	//   dummy[p] |= (BinSNPs[iMod][p][1] | BinSNPs[iMod][p][2]);
	//   nCarriers[i]+=bitcount64(dummy[p]);
	// }
	// delete[] dummy;
	if(map[i].analysis_in && nCarriers[i]<=NCT && nCarriers[i]!=0){
	  NCTcounter[currentwindow][nCarriers[i]-1]++;
	}
      }
    }

    if(!optimalrare){ // Find limits and positions with fixed NCT
      for(i=0; i<nwindows; i++){
	nRareLimits[i]=1;
	window[i].n_at_level=new int[1];
	if(!window[i].n_at_level) die("Memory allocation error in window[i].n_at_level!");
	window[i].levelpos=new int*[1];
	if(!window[i].levelpos) die("Memory allocation error in window[i].levelpos!");
	int carrierPerWindow = 0;
	for(j=0; j<NCT; j++){
	  carrierPerWindow+=NCTcounter[i][j];
	}
	window[i].n_at_level[0]=carrierPerWindow;
	window[i].levelpos[0]=new int[carrierPerWindow]();
	if(!window[i].levelpos[0]) die("Memory allocation error in window[i].levelpos[0]!");
      }
      int currentchr=0;
      int oldchr=0;
      int currentwindow=0;
      int carrierpos=0;
      for(i=0; i<nlinestped; i++){
	if(map[i].analysis_in){
	  currentchr=atoi(map[i].chr);
	  if(oldchr==0){
	    oldchr=atoi(map[i].chr);
	  }
	  else{
	    if(oldchr!=currentchr){
	      oldchr=currentchr;
	      currentwindow++;
	      carrierpos=0;
	    }
	  }
	  if(map[i].analysis_in && nCarriers[i]<=NCT &&  nCarriers[i]!=0){
	    window[currentwindow].levelpos[0][carrierpos]=i;
	    carrierpos++;
	  }
	}
      }
    }
    else if(optimalrare){
      for(i=0; i<nwindows; i++){
	for(j=0; j<NCT; j++){
	  if(NCTcounter[i][j]!=0){
	    nRareLimits[i]++;
	  }
	}
	window[i].n_at_level=new int[nRareLimits[i]]();
	if(!window[i].n_at_level) die("Memory allocation error in window[i].n_at_level!");
	window[i].levelpos=new int*[nRareLimits[i]];
	if(!window[i].levelpos) die("Memory allocation error in window[i].levelpos!");
	rareLimitsNCT[i]=new int[nRareLimits[i]]();
	if(!rareLimitsNCT[i]) die("Memory allocation error in rareLimitsNCT[i]!");
	rareLimitsNCTInverse[i]=new int[NCT]();
	if(!rareLimitsNCTInverse[i]) die("Memory allocation error in rareLimitsNCTInverse[i]!");
      }

      for(i=0; i<nwindows; i++){
	int currentlevel=0;
	for(j=0; j<NCT; j++){
	  int carrierPerWindowLevel=0;
	  if(NCTcounter[i][j]!=0){
	    for(int k=0; k<=j; k++){
	      carrierPerWindowLevel+=NCTcounter[i][k];
	    }
	    window[i].n_at_level[currentlevel]=carrierPerWindowLevel;
	    window[i].levelpos[currentlevel]=new int[carrierPerWindowLevel]();
	    if(!window[i].levelpos[currentlevel]) die("Memory allocation error in window[i].levelpos[currentlevel]!");
	    rareLimitsNCT[i][currentlevel]=j;
	    rareLimitsNCTInverse[i][j]=currentlevel;
	    currentlevel++;
	  }
	}
      }
      int currentwindow=0;
      int currentchr=0;
      int oldchr=0;
      int **NCTcounter2=NULL;
      NCTcounter2=new int*[nwindows];
      for(i=0; i<nwindows; i++){
	NCTcounter2[i]=new int[NCT]();
	if(!NCTcounter2[i]) die("Memory allocation error in NCTcounter2[i]!");
	//	memset(NCTcounter2[i], 0, NCT*sizeof(int));
      }
      for(i=0; i<nlinestped; i++){
	currentchr=atoi(map[i].chr);
	if(map[i].analysis_in){
	  if(oldchr==0){
	    oldchr=atoi(map[i].chr);
	  }
	  else{
	    if(oldchr!=currentchr){
	      oldchr=currentchr;
	      currentwindow++;
	    }
	  }
	  if(nCarriers[i]<=NCT && nCarriers[i]!=0){
	    window[currentwindow].levelpos[nRareLimits[currentwindow]-1][NCTcounter2[currentwindow][NCT-1]]=i;
	    //    cout<<window[currentwindow].levelpos[nRareLimits[currentwindow]-1][NCTcounter2[currentwindow][NCT-1]]<<endl;
	    //    cout<<i<<" "<<NCT<<endl;
	    NCTcounter2[currentwindow][NCT-1]++;
	    //    cout<<"NCTcounter2["<<currentwindow<<"]["<<nCarriers[i]-1<<"] "<<NCTcounter2[currentwindow][nCarriers[i]-1]<<" rareLimitsNCTInverse["<<currentwindow<<"]["<<nCarriers[i]-1<<"] "<<rareLimitsNCTInverse[currentwindow][nCarriers[i]-1]<<" "<<i<<" "<<nwindows<<" "<<NCT << endl;
	    //    cout<<"NCTcounter2[currentwindow]["<<nCarriers[i]-1<<"] "<<NCTcounter2[currentwindow][nCarriers[i]-1]<<endl;
	    //     for(j=rareLimitsNCTInverse[currentwindow][nCarriers[i]-1]; j<nRareLimits[currentwindow]; j++){
	    //       //      window[currentwindow].levelpos[rareLimitsNCTInverse[currentwindow][nCarriers[i]-1]][NCTcounter2[currentwindow][nCarriers[i]-1]]=i; // DOES NOT WORK YET??
	    //       NCTcounter2[currentwindow][rareLimitsNCT[currentwindow][nCarriers[i]-1]]++;
	    //       if(NCTcounter2[currentwindow][j]>0){
	    // window[currentwindow].levelpos[rareLimitsNCTInverse[currentwindow][j]-1][NCTcounter2[currentwindow][j]-1]=i; // DOES NOT WORK YET??
	    //       }
	    //     }

	    //    cout<<"window["<<currentwindow<<"].levelpos["<<rareLimitsNCTInverse[currentwindow][nCarriers[i]-1]<<"][NCTcounter2[currentwindow]["<<nCarriers[i]<<"-1]] "<<window[currentwindow].levelpos[rareLimitsNCTInverse[currentwindow][nCarriers[i]-1]][NCTcounter2[currentwindow][nCarriers[i]-1]]<<endl;
	    //    int k=0;
	    //    while(rareLimitsNCTInverse[currentchr][k]<=NCT){

	    //    }
	  }
	}
      }
      for(i=0; i<nwindows; i++){
	delete[] NCTcounter2[i];
	//	delete[] rareLimitsNCTInverse[i];
      }
      delete[] NCTcounter2;
      //      delete[] rareLimitsNCTInverse;
    }
    int carriersin=0;
    for(i=0; i<nwindows; i++){
      //      qsort(NCTcounter[i],NCT,sizeof(int),compareint);
      for(j=0; j<NCT; j++){
	carriersin+=NCTcounter[i][j];
      }
      delete[] NCTcounter[i];
    }
    delete[] NCTcounter;
    cout<<carriersin<<" SNPs with <= "<<NCT<<" carriers identified."<<endl;
    if(carriersin==0){
      cout<<"Not enough variants to proceed!\n";
      logfile<<"Not enough variants to proceed!\n";
      errorfile<<"Not enough variants to proceed!\n";
      exit(0);
    }

  }
}




int binRARE(int *cov, int sexcov, struct PERSON *person, int nlinestped, int nlinestfam, struct MAP *map, struct COUNTS *counts, bool optimalrare, double raref, int nsim, int n, int nwindows, int nwindowssinglebin, int *nRareSNPs, int **windowPositions, int *nRareLimits, double **rareLimits, int *MAF_Level_VT_FISHER, int *MAF_Level_VT_CMAT, int *MAF_Level_VT_COLL, int *MAF_Level_VT_REGRESSION, int *MAF_Level_VT_FRACREG, int *MAF_Level_VT_COLLREG, double *pCOLLvec, double *pCMATvec, double *pFISHERvec, double *pFISHERvecChi, double *pREGRESSIONvec, double *pFRACREG, double *pCOLLREG, double **limitsStatFISHER, double **limitsStatREGRESSION, double **limitsStatFRACREG,double **limitsStatCOLLREG, double **limitsStatCOLL, double **limitsStatCMAT, int thread, string rarefile,string intervalfile, double *rarefFISHER, double *rarefREGRESSION, double *rarefFRACREG,double *rarefCOLLREG, double *rarefCOLL, double *OR_COLL_vec, double *OR_COLL_f_vec, double *rarefCMAT, double *OR_CMAT_vec, double inflationfactor, fstream &errorfile, fstream &logfile, double wilsonpretest, double rareregpretest, int rarepretest, float rarepretestlimit, int *FISHERnotconverged, int *FISHERpretestpassed, int *REGRESSIONpretestpassed, int *FRACREGpretestpassed, int *COLLREGpretestpassed, int *CMATpretestpassed, int *COLLpretestpassed, struct WINDOW *window, double **FISHER_CI, double **CMAT_CI, double **COLL_CI, double **REGRESSION_CI, double **FRACREG_CI, double **COLLREG_CI, double *COLLREG_beta, double *COLLREG_se, double *FRACREG_beta, double *FRACREG_se, bool FISHERflag, bool REGRESSIONflag, bool FRACREGflag,bool COLLREGflag, bool COLLflag, bool CMATflag, int *FISHERcount, int *REGRESSIONcount, int *FRACREGcount,int *COLLREGcount, int *COLLcount, int *CMATcount, double *FISHERstats, double *REGRESSIONstats, double *FRACREGstats,double *COLLREGstats, double *COLLstats, double *CMATstats, double *FISHERpermstats, double *REGRESSIONpermstats,  double *FRACREGpermstats,double *COLLREGpermstats, double *COLLpermstats, double *CMATpermstats, int *completedFISHER, int *completedREGRESSION, int *completedFRACREG,int *completedCOLLREG, int *completedCOLL, int *completedCMAT, int *currentn, uint64_t*** BinSNPs, int nwordsSNPs, uint64_t** BinSNPsCCFlags,uint64_t** BinSNPsCovCatFlags, uint64_t** BinSNPsGenderFlags, int* SNPMapInverse, int* SNPMap, int* PPLMap, struct PPLLOCATION* PplLocations, int npeopleqc, int nsnpsqc, bool fisherCorrection, struct LOCAL_MATCHING* matchingRARE, uint32_t nMWindowsRARE, int rare_stratify, double fulltests, int teststat, int featcol, double casecounts5[3][3][3][3][3], double controlcounts5[3][3][3][3][3], double *p, double **newbeta, double **X, double **Xmod, double **Xt, double **A,  double **UNNT, double **VNN, double *S, double **Sinv, double **A0, double **Ainv, double **AinvXt, double *YY, double **Yt, double **YtX, double **YtXAinv, double **YtXAinvXt, double **Yminusp,  int N, int alt, double *Yhelp, int xType, int female, int male, double **D, double **T, double **U, double **Ut, double **sumPP, double **sumPJ, double **sumPK, double **MMinv,int test, int npplqc, struct STATplus result, /* struct STATplus result1_Rare, struct STATplus result2_Rare,*/ int numberOfAllCov, int ncov, int qt,int covariancematrix,  int xsinglevec[27], int zsinglevec[27], int ncasesqc, int ncontrolsqc, int nrestqc, int singleMarkerTest, int haplo, int maxthreads, int liabilityCut, int binamin, int binamax, int binadjust, int comparemode, int ncasesqcMale, int ncasesqcFemale,int ncontrolsqcMale, int ncontrolsqcFemale, int *ix, int *iy, int *iz, struct WEIGHTS weights, string SetIDfile, int intervaleditor, int setid, int verbose, int adaptive, int nCovCathegories, int *CovCathegories, int collinter, int **doublewindowcoord, int vb, int *nvbstartvt, int **vbstartvt, int *nvbend, int *vbstart, int **vbend, float *vbmaxpermstat, int dim1, int **nvblevel, int ***nchunks, int ****chunkpos, int ****chunklen, int *nCarriers, int ****ndummyatlevel, int *****dummypos, int *****dummylevel, uint64_t ***BinCarriers, float ***genoWeights, int dosage){

  int i;
  if(intervaleditor){
    write_setid(setid,SetIDfile,intervaleditor,raref,nwindowssinglebin,windowPositions,errorfile,logfile,intervals,map, window, intervalfile, nRareLimits);
  }

  if(n==0){
    bins2test=nwindows;
    if((CMATflag || FISHERflag || REGRESSIONflag || FRACREGflag || COLLREGflag)){
      weightvec=(double*)calloc(nlinestped,sizeof(double));
      if(!weightvec) die("Memory allocation error in weightvec!");
      weightvec=get_weights(nlinestped,weights,map,weightvec);
    }
  }
  //  else if(n>0 && adaptive!=0)
  struct STATplus result1,result2,resultDummy1,resultDummy2;
  int j;

  // DUMMY VARIABLES
  char *hapfile=NULL;
  char *hapstring=NULL;
  FILE *fptr10=NULL;
  int skip=0;
  int N1=1;
  int nMc=n;

  *currentn=0;

  int l, m;
  int start=0;
  int end=0;
  double stat=0;
  double dstat=0;
  double traref=0;
  double tp=0;

  int nSnps=0; //number of snps to be included with additive term
  int nSnpsDom=0; //number of snps to be included with dominance term
  int nEnvirons=0; // number of non-genetic parameters (for future extensions)
  int nInters=0; // number of interaction terms to be used (under construction)
  int *snps; // list of snps to be included with additive term, indices in BinSNPs
  int *snpsDom; // list of snps to be included with dominance term, indices in BinSNPs
  int *environs; // list of non-genetic parameters (for future extensions)
  int *inters; // list of interaction terms to be used (coding under construction)
  int nx=nSnps+nSnpsDom+nEnvirons+nInters; // total n
  if(sexcov==1){nx++;}

  snps = (int *) calloc(nSnps, sizeof(int));
  if(!snps) die("Memory allocation error in snps!");
  snpsDom = (int *) calloc(nSnpsDom, sizeof(int));
  if(!snpsDom) die("Memory allocation error in snpsDom!");
  environs = (int *) calloc(nEnvirons, sizeof(int));
  if(!environs) die("Memory allocation error in environs!");
  inters = (int *) calloc(nInters, sizeof(int));
  if(!inters) die("Memory allocation error in inters!");

  //paramters & Co for L2
  //choode all equal to zero to test full model
  int nSnps_second=0;
  int nSnpsDom_second=0;
  int nEnvirons_second=0;
  int nInters_second=0;
  int *snps_second;
  int *snpsDom_second;
  int *environs_second;
  int *inters_second;
  int nx_second=nSnps_second+nSnpsDom_second+nEnvirons_second+nInters_second;
  if(sexcov==1){nx_second++;}


  snps_second = (int *) calloc(nSnps_second, sizeof(int));
  if(!snps_second) die("Memory allocation error in snps_second!");
  //  cout <<COLLflag<<" "<<CMATflag<<" " << FISHERflag<<" "<<REGRESSIONflag<<" "<<FRACREGflag<<" "<<COLLREGflag<<endl;
  snpsDom_second = (int *) calloc(nSnpsDom_second, sizeof(int));
  if(!snpsDom_second) die("Memory allocation error in snpsDom_second!");
  environs_second = (int *) calloc(nEnvirons_second, sizeof(int));
  if(!environs_second) die("Memory allocation error in environs_second!");
  inters_second = (int *) calloc(nInters_second, sizeof(int));
  if(!inters_second) die("Memory allocation error in inters_second!");
  if(n==0){
    cout    << "\n\nRare variants: begin testing... " << nsim << " simulations" << endl;
    logfile << "\n\nRare variants: begin testing... " << nsim << " simulations" << endl;
  }


  struct COUNTS* counts_p = counts; // set to global, for rare_stratify inactive or n=0
  struct COUNTS* wCounts = NULL;
  uint64_t** BinSNPsCCFlags_local = BinSNPsCCFlags; // set to global, for rare_stratify inactive or n=0
  bool ifChrX=0;

  if (strcmp(map[start].chr, "23") == 0){ifChrX=1;}

  if(ifChrX && COLLflag){
    ncasesqcMale=0;
    ncasesqcFemale=0;
    ncontrolsqcMale=0;
    ncontrolsqcFemale=0;
    if(!qt){
      for (int i=0; i<nlinestfam; i++){
	if (person[i].qcin==1){
	  if (person[i].aff[thread]==2 &&
	      person[i].sex==1){ncasesqcMale++;}
	  else if (person[i].aff[thread]==2 &&
		   person[i].sex==2){ncasesqcFemale++;}
	  else if (person[i].aff[thread]==1 &&
		   person[i].sex==1){ncontrolsqcMale++;}
	  else if (person[i].aff[thread]==1 &&
		   person[i].sex==2){ncontrolsqcFemale++;}
	}
      }
    }
  }


  if (rare_stratify) {
    wCounts = new struct COUNTS[nlinestped];
	if(!wCounts) die("Memory allocation error in wCounts!");
	if (n>0) {
	  for (uint32_t w=0; w<nMWindowsRARE; w++) {
	    aff_permute_clusters(person, matchingRARE[w].groups, matchingRARE[w].nGroups, thread, ix, iy, iz);
	    matchingRARE[w].BinSNPsCCFlags = updateCC_bin_snp(matchingRARE[w].BinSNPsCCFlags, person, nlinestfam, nwordsSNPs);
        for (int k, kMod=matchingRARE[w].bin_start; kMod<=matchingRARE[w].bin_stop; kMod++) {
          k = SNPMap[kMod];
          updateCounts(&wCounts[k], nwordsSNPs, BinSNPs[kMod], matchingRARE[w].BinSNPsCCFlags, BinSNPsGenderFlags, (matchingRARE[w].chr==23 || matchingRARE[w].chr==24));
	    }
	    counts_p = wCounts;
      }
    }
  }


    // Bin loop
    for (l=0; l<nwindows; l++){



	if(collinter>0 && l<nwindowssinglebin) continue;

	if(collinter>0 || window[l].ntests>0){


	if (rare_stratify) {
      uint8_t chr = atoi(map[window[l].levelpos[nRareLimits[l]-1][0]].chr);
      int32_t pos = -1;
      uint32_t index = window[l].n_at_level[nRareLimits[l]-1]/2;
      while (pos<0) pos = SNPMapInverse[window[l].levelpos[nRareLimits[l]-1][index++]];
      LOCAL_MATCHING *matching = get_matchingRARE_window(chr, pos, matchingRARE, nMWindowsRARE);
      if (n>0) BinSNPsCCFlags_local = (*matching).BinSNPsCCFlags;
	}


	if(FISHERflag && ((wilsonpretest==0 && rareregpretest==0) || FISHERpretestpassed[l]==1)){
	  int countchange=0;

	  double newValue=0;
	  int helpi;
	  int i;
	  init(&(resultDummy1),0,1,numberOfAllCov,0,1,0,0);
	  init(&(resultDummy2),0,1,numberOfAllCov,0,1,0,0);

	  //	  for(i=start; i<=end; i++){
	  for(int j=0; j<window[l].n_at_level[nRareLimits[l]-1]; j++){
	    //from top level
	    i=window[l].levelpos[nRareLimits[l]-1][j];
	    int iMod = SNPMapInverse[i];

	    if (iMod==-1) continue;
	    //	    if (map[i].analysis_in == 1 && ((map[i].mafr-raref)<EPS || fabs(map[i].mafr-raref)<EPS) && map[i].mafr>EPS){

	    newValue = 1;
	    helpi = i;
	    if (singleMarkerTest == 1){
	      if (strcmp(map[i].chr, "23") == 0){
		newValue = snpTestX(counts[i], teststat, &inflationfactor, &fulltests);
	      }
	      else{
		newValue = armitageTest(counts[i], teststat, &inflationfactor, &fulltests);
	      }
	    }
	    else if (singleMarkerTest == 2){
	      if (strcmp(map[i].chr, "23") == 0){
		newValue = snpTestX(counts[i], teststat, &inflationfactor, &fulltests);
	      }
	      else{
		newValue = genotypTest(counts[i], teststat, &inflationfactor, &fulltests);
	      }
	    }
	    else if (singleMarkerTest == 3 || singleMarkerTest == 4){
	      int alt, xType;
	      double df, tLog, Fstat;
	      if (strcmp(map[i].chr, "23") == 0){
		alt=1;
		xType = 3;
	      }
	      else{
		alt=1;
		xType = 0;
	      }
	      if(!qt){
		result1 = logreg(xsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, i, 0, 0, inflationfactor, casecounts5, controlcounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[i], D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread,npplqc, PPLMap, BinSNPs, PplLocations,1,resultDummy1,numberOfAllCov,0,singleMarkerTest,-1,0,0,0,dosage,genoWeights);
		alt=0;
		if(result1.df>0){
		  result2 = logreg(zsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, i, 0, 0, inflationfactor, casecounts5, controlcounts5, p, newbeta, X, Xmod, Xt, A,  UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[i], D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread,npplqc, PPLMap, BinSNPs, PplLocations,1,resultDummy2,numberOfAllCov,0,singleMarkerTest,-1,0,0,0,dosage,genoWeights);
		}
		else{
		  result2.df=0;result2.sc=0;
		}
		tLog = 2*(result1.sc -result2.sc);
		df = result1.df -result2.df;
		if (df >= 1){
		  newValue = pValueCalc(df/2, tLog/2);
		  /*fulltests++;
		    if (singleMarkerTest == 3 || df ==2)
		    {inflationfactor+=tLog;}
		    else {inflationfactor+=2*tLog;}*/
		}
		else{
		  newValue = 1;
		}
	      }
	      else{//qt singlemarker
		newValue = 1;
		result1 = qtreg(xsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, i, 0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1,alt,Yhelp, xType, female, male, counts[i], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread,hapfile,0,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,1,resultDummy1,0,0,numberOfAllCov,0,singleMarkerTest,-1,0,0,0);

		if(result1.df>0){
		  alt=0;
		  result2 = qtreg(zsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, i, 0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[i], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread,hapfile,0,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs,  PplLocations,1,resultDummy2,0,0,numberOfAllCov,0,singleMarkerTest,-1,0,0,0);
		}
		else{
		  result2.df=0;result2.sc=0;
		}
		df = result2.df - result1.df; //!!
		Fstat = ((result2.sc -result1.sc)/(df))/(result1.sc/result1.df);

		if (df >= 1 && result2.sc > 0.000001){
		  newValue = betai(result1.df/2,df/2,result1.df/(result1.df+df*Fstat));
		  /*fulltests++;
		    if (singleMarkerTest == 3 || df ==2)
		    {inflationfactor+=Fstat;}
		    else {inflationfactor+=2*Fstat;}*/
		}
		else{
		  newValue = 1;
		}
	      }
	    }  //end singlemarker 3 + 4
	    //	      }
	    counts[i].pSingle=newValue;
	    counts_p[i].pSingle=newValue;
	  }
	  freeResult(resultDummy1,1);
	  freeResult(resultDummy2,1);

	  // FISHER calculation
	  // Test statistic with FT
	  if(optimalrare==0){
	    dstat=calcFISHER(cov, sexcov, person, ncasesqc, ncontrolsqc, nrestqc, thread, nlinestfam, nlinestped, start, end, counts_p, raref, inflationfactor, fisherCorrection, map, weightvec, window, l, 0);
	    if(n==0){
	      FISHERstats[l]=dstat;
	    }
	    else if(n!=0){
	      FISHERpermstats[l]=dstat;
	    }
	  }
	  // p as test statistic with VT
	  if(optimalrare==1){
	    for(m = 0; m < nRareLimits[l]; m++){
	      traref = (double)rareLimits[l][m];
	      limitsStatFISHER[l][m] = calcFISHER(cov, sexcov, person, ncasesqc, ncontrolsqc, nrestqc, thread, nlinestfam, nlinestped, start, end, counts_p,  traref, inflationfactor, fisherCorrection, map, weightvec, window, l, m);
	      dstat=calcFISHERp(limitsStatFISHER[l][m], traref, start, end, counts_p, thread, map, window, l, m); // l=bin, m=level

	      if(m==0){
		if(n==0){
		  FISHERstats[l]=dstat;
		  rarefFISHER[l]=traref;
		}
		else if(n!=0){
		  FISHERpermstats[l]=dstat;
		}
	      }
	      else if(m!=0){
		if(n==0){
		  if(dstat*(1-EPS)<FISHERstats[l]){
		    FISHERstats[l]=dstat;
		    rarefFISHER[l]=traref;
		    MAF_Level_VT_FISHER[l]=m;
		  }
		}
		else if(n!=0){
		  if(dstat*(1-EPS)<FISHERpermstats[l]){
		    FISHERpermstats[l]=dstat;
		  }
		}
	      }
	    }
	  }


	  ////////////////////////////////////////////////////////
	  if(!(dstat>0)){
	    int j=window[l].n_at_level[m];

	    if(n==0){ // Cannot compare to permutatoin values; finished

	      logfile << "n="<<n<<": p value of "<< limitsStatFISHER[l][m] << " with "<< 2*j <<" df is " << dstat << flush<<endl;
	      errorfile << "n="<<n<<": p value of "<< limitsStatFISHER[l][m] << " with "<< 2*j <<" df is " << dstat << flush<<endl;
	      if(wilsonpretest!=0 || rareregpretest!=0){
		FISHERpretestpassed[l]-=1;
	      }
	      window[l].ntests-=1;
	      if(window[l].ntests==0){
		bins2test--;
	      }
	    }
	    else{
#pragma omp atomic
	      FISHERnotconverged[l]++;
	      dstat=0;

	    }
	  }

	  ////////////////////////////////////////////////////////





	  if(n==0 && rareregpretest!=0){
	    if(optimalrare){
	      if(FISHERstats[l]*(1+EPS)>rareregpretest){
#pragma omp atomic
		FISHERpretestpassed[l]-=1;
#pragma omp atomic
		window[l].ntests-=1;
		if(window[l].ntests==0){
#pragma omp atomic
		  bins2test--;
		}
		pFISHERvec[l]=FISHERstats[l];
	      }
	    }
	    else if(!optimalrare){
	      double FRpre=calcFISHERp(dstat, raref, start, end, counts_p, thread, map, window, l, 0);
	      if(FRpre>rareregpretest){
#pragma omp atomic
		FISHERpretestpassed[l]-=1;
#pragma omp atomic
		window[l].ntests-=1;
		if(window[l].ntests==0){
#pragma omp atomic
		  bins2test--;
		}
		pFISHERvec[l]=FRpre;
	      }
	    }
	  } // n=0 with prescreen

	  else if(n!=0){
#pragma omp atomic
	    completedFISHER[l]++;
	    if(optimalrare==0 && (FISHERpermstats[l]*(1+EPS)>FISHERstats[l])){
#pragma omp atomic
	      FISHERcount[l]++;
#pragma omp atomic
	      countchange++;
	    }
	    else if(optimalrare == 1 && (FISHERpermstats[l]*(1-EPS)<FISHERstats[l])){
#pragma omp atomic
	      FISHERcount[l]++;
#pragma omp atomic
	      countchange++;
	    }
	    //#pragma omp critical
	    //	  {
	    if(n>0 && (wilsonpretest!=0 && FISHERpretestpassed[l]) && countchange>0){
#pragma omp atomic
	      countchange-=1;
	      double CI_low=wilson_lower(completedFISHER[l],FISHERcount[l]);
	      if(CI_low>wilsonpretest){
		pFISHERvec[l]=(double)FISHERcount[l]/completedFISHER[l];
#pragma omp atomic
		FISHERpretestpassed[l]-=1;
#pragma omp atomic
		window[l].ntests-=1;
		if(window[l].ntests==0){
#pragma omp atomic
		  bins2test--;
		}
		FISHER_CI[l][0]=CI_low;
		FISHER_CI[l][1]=wilson_upper(completedFISHER[l],FISHERcount[l]);
		window[l].nperm[0]=completedFISHER[l];
	      }
	    }
	    //	  }
	  }
	} // if FISHERflag

	if(REGRESSIONflag && ((wilsonpretest==0 && rareregpretest==0) || REGRESSIONpretestpassed[l]==1)){
	  int countchange=0;

	  double newValue=0;
	  int helpi;
	  int i;
	  int alt, xType;
	  double df, tLog, Fstat;
	  int collapseRare=0;
	  int collcollapseRare=0;
	  int firstbinlastSNP=0;
	  double newValueREGRESSION=1;
	  nSnps=0;
	  for(m = 0; m < nRareLimits[l]; m++){

	    traref=(double)rareLimits[l][m];

	    //		for(i=start; i<=end; i++){
	    for(int j=0; j<window[l].n_at_level[m]; j++){
	      i=window[l].levelpos[m][j];

	      int iMod = SNPMapInverse[i];
	      if (iMod==-1) continue;

	      //		    if (map[i].analysis_in == 1 && ((map[i].mafr-traref)<EPS || fabs(map[i].mafr-traref)<EPS) && map[i].mafr>EPS){
	      newValue = 1;
	      helpi = i;
	      nSnps++;
	      snps = (int *) realloc(snps, nSnps*sizeof(int));
	      if(!snps){
		cout<<"Memory error in snps[]!"<<endl;
		logfile<<"Memory error in snps[]!"<<endl;
		errorfile<<"Memory error in snps[]!"<<endl;
		exit(1);
	      }
	      snps[nSnps-1]=i;

	      if (strcmp(map[i].chr, "23") == 0 && nSnps==1){
		alt=1;xType = 3;
	      }
	      else if(nSnps==1){
		alt=1;xType = 0;
	      }
	    }

	    nx=nSnps+nSnpsDom+nEnvirons+nInters;
	    if(qt){
	      init(&(resultDummy1),0,1,numberOfAllCov+nx,0,1,0,0);
	      init(&(resultDummy2),0,1,numberOfAllCov+nx,0,1,0,0);
	    }
	    else{
	      init(&(resultDummy1),0,0,numberOfAllCov+nx,0,1,0,0);
	      init(&(resultDummy2),0,0,numberOfAllCov+nx,0,1,0,0);
	    }

	    if(sexcov==1){nx++;}
	    nx_second=nSnps_second+nSnpsDom_second+nEnvirons_second+nInters_second;
	    if(sexcov==1){nx_second++;}

	    alt=1;
	    collapseRare=0;
	    collcollapseRare=0;
	    newValueREGRESSION = 1;
	    firstbinlastSNP=nSnps; // information about number of snps to be used in alt=0

	    if(qt){
	      result1 = regGeneral(xsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc,
				   i, 0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv,
				   A0, Ainv, AinvXt, YY, Yminusp, N1,alt,Yhelp, xType, female, male, counts[i], Yt, YtX,
				   YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread,hapfile,0,
				   hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,1,resultDummy1,0,0,numberOfAllCov,ncov,0,
				   singleMarkerTest,-1,
				   nSnps,nSnpsDom,nEnvirons,nInters,snps,snpsDom,environs,inters,nx,dim1,collapseRare,collcollapseRare,qt,weightvec, firstbinlastSNP);
	    }
	    else{
	      result1 = logRegGeneral(xsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc,
				      i, 0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv,
				      A0, Ainv, AinvXt, YY, Yminusp, N1,alt,Yhelp, xType, female, male, counts[i], Yt, YtX,
				      YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread,hapfile,0,
				      hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,1,resultDummy1,0,0,numberOfAllCov,ncov,0,
				      singleMarkerTest,-1,
				      nSnps,nSnpsDom,nEnvirons,nInters,snps,snpsDom,environs,inters,nx,dim1,collapseRare,collcollapseRare,qt,weightvec,covariancematrix, firstbinlastSNP, collinter);

	    }
	    if(result1.df>0){
	      alt=0;
	      if(qt){
		result2 = regGeneral(zsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc,
				     i,0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv,
				     A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[i], Yt, YtX,
				     YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread,hapfile,0,
				     hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs,  PplLocations,1,resultDummy2,0,0,numberOfAllCov,ncov,
				     0,singleMarkerTest,-1,
				     nSnps_second,nSnpsDom_second,nEnvirons_second,nInters_second,snps_second,snpsDom_second,environs_second,inters_second,nx_second,dim1,collapseRare,collcollapseRare,qt,weightvec, firstbinlastSNP);
	      }
	      else{
		result2 = logRegGeneral(zsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc,
					i,0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv,
					A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[i], Yt, YtX,
					YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread,hapfile,0,
					hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs,  PplLocations,1,resultDummy2,0,0,numberOfAllCov,ncov,
					0,singleMarkerTest,-1,
					nSnps_second,nSnpsDom_second,nEnvirons_second,nInters_second,snps_second,snpsDom_second,environs_second,inters_second,nx_second,dim1,collapseRare,collcollapseRare,qt,weightvec,covariancematrix, firstbinlastSNP, collinter);
	      }
	    }
	    else{
	      result2.df=0;result2.sc=0;
	    }

	    if(qt){
	      df = result2.df - result1.df;
	      Fstat = ((result2.sc -result1.sc)/(df))/(result1.sc/result1.df);
	      if (df >= 1 && result2.sc > 0.000001){
		newValueREGRESSION = betai(result1.df/2,df/2,result1.df/(result1.df+df*Fstat));
	      }
	      else{
		newValueREGRESSION = 1;
	      }
	    }
	    else if(!qt){
	      tLog = 2*(result1.sc -result2.sc);
	      df = result1.df -result2.df;

	      if (df >= 1){
		  newValueREGRESSION = pValueCalc(df/2, tLog/2);
		}
	      else{
		newValueREGRESSION = 1;
	      }
	    }
	    nSnps=0;
	    freeResult(resultDummy1,1);
	    freeResult(resultDummy2,1);
	    if(optimalrare==0){

	      if(nsim==0){
		pREGRESSIONvec[l]=newValueREGRESSION;
	      }
	      if(n!=0){
		if((wilsonpretest==0 && rareregpretest==0) || REGRESSIONpretestpassed[l]==1){
#pragma omp atomic
		  completedREGRESSION[l]++;

		  if(!(newValueREGRESSION*(1-EPS)>0)){
		    cout <<" BUG "<<newValueREGRESSION*(1-EPS)<<" "<<newValueREGRESSION<<endl;
		    exit(1);
		  }
		  if(newValueREGRESSION*(1-EPS)<REGRESSIONstats[l]){ // p-value by permutation
#pragma omp atomic
		    REGRESSIONcount[l]++;
#pragma omp atomic
		    countchange++;
		  }
		}
	      }
	      if(n==0){
		REGRESSIONstats[l]=newValueREGRESSION;
		rarefREGRESSION[l]=traref;
		//		MAF_Level_VT_REGRESSION[l]=m;
	      }
	    }
	    else if(optimalrare==1){
	      if(n!=0){
		if((wilsonpretest==0 && rareregpretest==0) || REGRESSIONpretestpassed[l]==1){
#pragma omp atomic
		  completedREGRESSION[l]++;
		  //			     if(REGRESSIONpermstats[l]*(1-EPS)<REGRESSIONstats[l]){ // p-value by permutation
		  if(newValueREGRESSION*(1-EPS)<REGRESSIONstats[l]){ // p-value by permutation
#pragma omp atomic
		    REGRESSIONcount[l]++;
#pragma omp atomic
		    countchange++;
		    break;
		  }
		}
	      }
	      else if(n==0 && m==0){
		REGRESSIONstats[l]=newValueREGRESSION;
		rarefREGRESSION[l]=traref;
		//		MAF_Level_VT_REGRESSION[l]=m;
	      }
	      else if(n==0 && m>0){
		if(newValueREGRESSION*(1-EPS)<REGRESSIONstats[l]){ // p-value by permutation			   if()
		  REGRESSIONstats[l]=newValueREGRESSION;
		  rarefREGRESSION[l]=traref;
		  MAF_Level_VT_REGRESSION[l]=m;
		}
	      }
	      if(n!=0 && m==0){
		REGRESSIONpermstats[l]=newValueREGRESSION;
	      }
	      else if(n!=0 && m!=0){
		if(newValueREGRESSION*(1-EPS)<REGRESSIONpermstats[l]){ // p-value by permutation
		  REGRESSIONpermstats[l]=newValueREGRESSION;
		}
	      }
	    }
	  } // loop over levels
	    /*
	      if(optimalrare==1 && n!=0){
	      if(window[l].REGRESSIONpretestpassed==1){
	      if(REGRESSIONpermstats[l]*(1-EPS)<REGRESSIONstats[l]){ // p-value by permutation
	      #pragma omp atomic
	      REGRESSIONcount[l]++;
	      #pragma omp atomic
	      countchange++;
	      }
	      }
	      }
	    */
	    // PRETEST
	  if(n>0 && wilsonpretest!=0 && countchange>0 && REGRESSIONpretestpassed[l]){
#pragma omp atomic
	    countchange-=1;

	    double CI_low=wilson_lower(completedREGRESSION[l],REGRESSIONcount[l]);
	    if(CI_low>wilsonpretest){
	      pREGRESSIONvec[l]=(double)REGRESSIONcount[l]/completedREGRESSION[l];
#pragma omp atomic
	      REGRESSIONpretestpassed[l]-=1;
#pragma omp atomic
	      window[l].ntests-=1;

	      if(window[l].ntests==0){
#pragma omp atomic
		bins2test--;
	      }
	      REGRESSION_CI[l][0]=CI_low;
	      REGRESSION_CI[l][1]=wilson_upper(completedREGRESSION[l],REGRESSIONcount[l]);
	      window[l].nperm[3]=completedREGRESSION[l];
	    }
	  }
	  else if(rareregpretest!=0 && n==0){

	    if(REGRESSIONstats[l]*(1+EPS)>rareregpretest){
	      //#pragma omp atomic
	      pREGRESSIONvec[l]=REGRESSIONstats[l];
#pragma omp atomic
		REGRESSIONpretestpassed[l]-=1;
#pragma omp atomic
	      window[l].ntests-=1;
	      if(window[l].ntests==0){
#pragma omp atomic
		bins2test--;
	      }

	    } // end pretest
	  }
	} // if REGRESSION

	if(FRACREGflag && ((wilsonpretest==0 && rareregpretest==0) || FRACREGpretestpassed[l]==1)){
	  int countchange=0;
	  double newValue=0;
	  int helpi;
	  int i;
	  int alt, xType;
	  double df, tLog, Fstat;
	  int collapseRare=0;
	  int collcollapseRare=0;
	  int firstbinlastSNP=0;
	  double newValueFRACREG=1;

	  nSnps=0;
	  for(m = 0; m < nRareLimits[l]; m++)
	    {
	      if(qt){
		init(&(resultDummy1),0,0,numberOfAllCov,0,1,0,0);
		init(&(resultDummy2),0,0,numberOfAllCov,0,1,0,0);
	      }
	      else if(!qt){
		init(&(resultDummy1),0,0,numberOfAllCov,0,1,0,0);
		init(&(resultDummy2),0,0,numberOfAllCov,0,1,0,0);
	      }
	      traref = (double)rareLimits[l][m];
	      //		  for(i=start; i<=end; i++){
	      for(int j=0; j<window[l].n_at_level[m]; j++){
		i=window[l].levelpos[m][j];

		int iMod = SNPMapInverse[i];
		if (iMod==-1) continue;

		//		      if (map[i].analysis_in == 1 && ((map[i].mafr-traref)<EPS || fabs(map[i].mafr-traref)<EPS) && map[i].mafr>EPS){
		newValue = 1;
		helpi = i;
		nSnps++;
		snps = (int *) realloc(snps, nSnps*sizeof(int));
		if(!snps){
		  cout<<"Memory error in snps[]!"<<endl;
		  logfile<<"Memory error in snps[]!"<<endl;
		  errorfile<<"Memory error in snps[]!"<<endl;
		  exit(1);
		}
		snps[nSnps-1]=i;

		if (strcmp(map[i].chr, "23") == 0 && nSnps==1){
		  alt=1;xType = 3;
		}
		else if(nSnps==1){
		  alt=1;xType = 0;
		}
	      }
	      nx=nSnps+nSnpsDom+nEnvirons+nInters;
	      if(sexcov==1){nx++;}

	      nx_second=nSnps_second+nSnpsDom_second+nEnvirons_second+nInters_second;
	      if(sexcov==1){nx_second++;}

	      alt=1;
	      collapseRare=1;
	      newValueFRACREG = 1;
	      if(qt){
		result1 = regGeneral(xsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc,
					   i, 0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv,
					   A0, Ainv, AinvXt, YY, Yminusp, N1,alt,Yhelp, xType, female, male, counts[i], Yt, YtX,
					   YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread,hapfile,0,
					   hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,0,resultDummy1,0,0,numberOfAllCov,ncov,0,
					   singleMarkerTest,-1,
				     nSnps,nSnpsDom,nEnvirons,nInters,snps,snpsDom,environs,inters,nx,dim1,collapseRare,collcollapseRare,qt,weightvec, firstbinlastSNP);
	      }
	      else{
		result1 = logRegGeneral(xsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc,
					      i, 0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv,
					      A0, Ainv, AinvXt, YY, Yminusp, N1,alt,Yhelp, xType, female, male, counts[i], Yt, YtX,
					      YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread,hapfile,0,
					      hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,0,resultDummy1,0,0,numberOfAllCov,ncov,0,
					      singleMarkerTest,-1,
					nSnps,nSnpsDom,nEnvirons,nInters,snps,snpsDom,environs,inters,nx,dim1,collapseRare,collcollapseRare,qt,weightvec, covariancematrix, firstbinlastSNP, collinter);
					      }
	      if(result1.df>0){
		alt=0;
		if(qt) {
		  result2 = regGeneral(zsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc,
				       i,0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv,
				       A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[i], Yt, YtX,
				       YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread,hapfile,0,
				       hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs,  PplLocations,0,resultDummy2,0,0,numberOfAllCov,ncov,
				       0,singleMarkerTest,-1,
				       nSnps_second,nSnpsDom_second,nEnvirons_second,nInters_second,snps_second,snpsDom_second,environs_second,inters_second,nx_second,dim1,collapseRare,collcollapseRare,qt,weightvec, firstbinlastSNP);
		}
		else{
		  result2 = logRegGeneral(zsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc,
					  i,0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv,
					  A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[i], Yt, YtX,
					  YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread,hapfile,0,
					  hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs,  PplLocations,0,resultDummy2,0,0,numberOfAllCov,ncov,
					  0,singleMarkerTest,-1,
					  nSnps_second,nSnpsDom_second,nEnvirons_second,nInters_second,snps_second,snpsDom_second,environs_second,inters_second,nx_second,dim1,collapseRare,collcollapseRare,qt,weightvec, covariancematrix, firstbinlastSNP, collinter);
					  }

	      }
	      else{
		result2.df=0;result2.sc=0;
	      }
	      if(qt){
		df = result2.df - result1.df;
		Fstat = ((result2.sc -result1.sc)/(df))/(result1.sc/result1.df);

		if (df >= 1 && result2.sc > 0.000001){
		  newValueFRACREG = betai(result1.df/2,df/2,result1.df/(result1.df+df*Fstat));
		}
		else{
		  newValueFRACREG = 1;
		}
	      }
	      else if(!qt){
		tLog = 2*(result1.sc -result2.sc);
		df = result1.df -result2.df;

		if (df >= 1){
		  newValueFRACREG = pValueCalc(df/2, tLog/2);
		}
		else{
		  newValueFRACREG = 1;
		}
	      }


	    if(!optimalrare && l<nwindowssinglebin){
	      FRACREG_beta[l]=result1.b[1];
	      FRACREG_se[l]=result1.betaNew_se[1];
	    }


	      nSnps=0;
	      freeResult(resultDummy1,1);
	      freeResult(resultDummy2,1);

	      if(optimalrare==0){
		if(nsim==0){
		  pFRACREGvec[l]=newValueFRACREG;
		}
		if(n!=0){
		  if((wilsonpretest==0 && rareregpretest==0) || FRACREGpretestpassed[l]==1){
#pragma omp atomic
		    completedFRACREG[l]++;
		    if(!(newValueFRACREG*(1-EPS)>0)){
		      cout <<" BUG "<<newValueFRACREG*(1-EPS)<<" "<<newValueFRACREG<<endl;
		      exit(1);
		    }
		    if(newValueFRACREG*(1-EPS)<FRACREGstats[l]){ // p-value by permutation
#pragma omp atomic
		      FRACREGcount[l]++;
#pragma omp atomic
		      countchange++;
		    }
		  }
		}
		else if(n==0){
		  FRACREGstats[l]=newValueFRACREG;
		  rarefFRACREG[l]=traref;
		  //		  MAF_Level_VT_FRACREG[l]=m;
		}
	      }
	      else if(optimalrare==1){
		if(n!=0){
		  if((wilsonpretest==0 && rareregpretest==0) || FRACREGpretestpassed[l]==1){
#pragma omp atomic
		    completedFRACREG[l]++;
		    if(newValueFRACREG*(1-EPS)<FRACREGstats[l]){ // p-value by permutation
#pragma omp atomic
		      FRACREGcount[l]++;
#pragma omp atomic
		      countchange++;
		      break;
		    }
		  }
		}
		if(n==0 && m==0){
		  FRACREGstats[l]=newValueFRACREG;
		  rarefFRACREG[l]=traref;
		}
		else if(n==0 && m>0){
		  if(newValueFRACREG*(1-EPS)<FRACREGstats[l]){ // p-value by permutation			   if()
		    FRACREGstats[l]=newValueFRACREG;
		    rarefFRACREG[l]=traref;
		    MAF_Level_VT_FRACREG[l]=m;
		  }
		}
		if(n!=0 && m==0){
		  FRACREGpermstats[l]=newValueFRACREG;
		}
		else if(n!=0 && m!=0){
		  if(newValueFRACREG*(1-EPS)<FRACREGpermstats[l]){ // p-value by permutation			   if()
		    FRACREGpermstats[l]=newValueFRACREG;
		  }
		}
	      }
	    } // loop over levels
	  /*
	    if(optimalrare==1 && n!=0){
	    if(window[l].FRACREGpretestpassed==1){
	    if(FRACREGpermstats[l]*(1-EPS)<FRACREGstats[l]){ // p-value by permutation
	    #pragma omp atomic
	    FRACREGcount[l]++;
	    }
	    }
	    }
	  */
	  // PRETEST
	  if(n>0 && wilsonpretest!=0 && countchange>0 && FRACREGpretestpassed[l]){
#pragma omp atomic
	    countchange--;
	    double CI_low=wilson_lower(completedFRACREG[l],FRACREGcount[l]);
	    if(CI_low>wilsonpretest){
	      pFRACREGvec[l]=(double)FRACREGcount[l]/completedFRACREG[l];
	      FRACREGpretestpassed[l]=0;
#pragma omp atomic
	      window[l].ntests-=1;
	      if(window[l].ntests==0){
#pragma omp atomic
		bins2test--;
	      }
	      FRACREG_CI[l][0]=CI_low;
	      FRACREG_CI[l][1]=wilson_upper(completedFRACREG[l],FRACREGcount[l]);
	      window[l].nperm[4]=completedFRACREG[l];
	    }
	  }
	  else if(rareregpretest!=0 && n==0){
	    if(FRACREGstats[l]*(1+EPS)>rareregpretest){
	      pFRACREGvec[l]=FRACREGstats[l];
#pragma omp atomic
		FRACREGpretestpassed[l]-=1;
	      window[l].ntests-=1;
	      if(window[l].ntests==0){
#pragma omp atomic
		bins2test--;
	      }
	    }
	  } // end pretest


	} // if COLLREG
	if(COLLREGflag && ((wilsonpretest==0 && rareregpretest==0) || COLLREGpretestpassed[l]==1)){
	  int countchange=0;
	  double newValue=0;
	  int helpi;
	  int i;
	  int alt, xType;
	  double df, tLog, Fstat;
	  int collapseRare=0;
	  int collcollapseRare=0;
	  int firstbinlastSNP=0;
	  double newValueCOLLREG=1;
	  nSnps=0;

	  int nRareLimits_tmp;

	  if(collinter==3 || collinter==4){
	    nRareLimits_tmp=1;
	  }
	  else{
	    nRareLimits_tmp=nRareLimits[l];
	  }

	  for(m = 0; m < nRareLimits_tmp; m++){
	    if(qt){
	      init(&(resultDummy1),0,0,numberOfAllCov,0,1,0,0);
	      init(&(resultDummy2),0,0,numberOfAllCov,0,1,0,0);
	    }
	    else if(!qt){
	      init(&(resultDummy1),0,0,numberOfAllCov,0,1,0,0);
	      init(&(resultDummy2),0,0,numberOfAllCov,0,1,0,0);
	    }

	  if(collinter==3 || collinter==4){
	    if((rareLimits[doublewindowcoord[l][0]][m]-rareLimits[doublewindowcoord[l][1]][m])>EPS){
	      traref=(double)rareLimits[doublewindowcoord[l][0]][m];
	    }
	    else{
	      traref=(double)rareLimits[doublewindowcoord[l][1]][m];
	    }
	  }
	  else{
	    traref = (double)rareLimits[l][m];
	  }
	  if(collinter==3 || collinter==4){
	    nInters=1;
	    firstbinlastSNP = window[doublewindowcoord[l][0]].n_at_level[0]-1; // index
	    //	      firstbinlastSNP = windowPositions[doublewindowcoord[l][0]][1]; // poistion instead of index
	  }

	  int n_at_level_tmp=0;
	  if(collinter==3 || collinter==4){
	    n_at_level_tmp=window[doublewindowcoord[l][0]].n_at_level[m]+window[doublewindowcoord[l][1]].n_at_level[m];
	  }
	  else{
	    n_at_level_tmp=window[l].n_at_level[m];
	  }
	  for(int j=0; j<n_at_level_tmp; j++){
	    if(collinter==3 || collinter==4){
	      if(j<window[doublewindowcoord[l][0]].n_at_level[m]){
		i=window[doublewindowcoord[l][0]].levelpos[m][j];
	      }
	      else{
		i=window[doublewindowcoord[l][1]].levelpos[m][j-window[doublewindowcoord[l][0]].n_at_level[m]];
	      }
	    }
	    else{
	      i=window[l].levelpos[m][j];
	    }
	      int iMod = SNPMapInverse[i];
	      if (iMod==-1) continue;
	      newValue = 1;
	      helpi = i;
	      nSnps++;

	      snps = (int *) realloc(snps, nSnps*sizeof(int));
	      if(!snps)die("Memory error in snps[]!");

	      snps[nSnps-1]=i;

	      if (strcmp(map[i].chr, "23") == 0 && nSnps==1){
		alt=1;xType = 3;
	      }
	      else if(nSnps==1){
		alt=1;xType = 0;
	      }
	    }// end i-loop rare snps of bin

	    nx=nSnps+nSnpsDom+nEnvirons+nInters;
	    //if(sexcov==1){nx++;}

		if(collinter==4){
		  nSnps_second=nSnps;
		  nSnpsDom_second=nSnpsDom;
		  nEnvirons_second=nEnvirons;

            }

	    nx_second=nSnps_second+nSnpsDom_second+nEnvirons_second+nInters_second;
	    //if(sexcov==1){nx_second++;}

	    alt=1;
	    collcollapseRare=1;
	    newValueCOLLREG = 1;
	    if(qt){
	      result1= regGeneral(xsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc,
				  i, 0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv,
				  A0, Ainv, AinvXt, YY, Yminusp, N1,alt,Yhelp, xType, female, male, counts[i], Yt, YtX,
				  YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread,hapfile,0,
				  hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,0,resultDummy1,0,0,numberOfAllCov,ncov,0,
				  singleMarkerTest,-1,
				  nSnps,nSnpsDom,nEnvirons,nInters,snps,snpsDom,environs,inters,nx,dim1,collapseRare,collcollapseRare,qt,weightvec, firstbinlastSNP);
	    }
	    else{
	      result1= logRegGeneral(xsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc,
				     i, 0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv,
				     A0, Ainv, AinvXt, YY, Yminusp, N1,alt,Yhelp, xType, female, male, counts[i], Yt, YtX,
				     YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread,hapfile,0,
				     hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,0,resultDummy1,0,0,numberOfAllCov,ncov,0,
				     singleMarkerTest,-1,
				     nSnps,nSnpsDom,nEnvirons,nInters,snps,snpsDom,environs,inters,nx,dim1,collapseRare,collcollapseRare,qt,weightvec, covariancematrix, firstbinlastSNP, collinter);
				     }
	    if(result1.df>0){
	      alt=0;
	      if(qt){
		result2 = regGeneral(zsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc,
				     i,0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv,
				     A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[i], Yt, YtX,
				     YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread,hapfile,0,
				     hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs,  PplLocations,0,resultDummy2,0,0,numberOfAllCov,ncov,
				     0,singleMarkerTest,-1,
				     nSnps_second,nSnpsDom_second,nEnvirons_second,nInters_second,snps_second,snpsDom_second,environs_second,inters_second,nx_second,dim1,collapseRare,collcollapseRare,qt,weightvec, firstbinlastSNP);
	      }
	      else{
		result2 = logRegGeneral(zsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc,
					i,0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv,
					A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[i], Yt, YtX,
					YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread,hapfile,0,
					hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs,  PplLocations,0,resultDummy2,0,0,numberOfAllCov,ncov,
					0,singleMarkerTest,-1,
					nSnps_second,nSnpsDom_second,nEnvirons_second,nInters_second,snps_second,snpsDom_second,environs_second,inters_second,nx_second,dim1,collapseRare,collcollapseRare,qt,weightvec, covariancematrix, firstbinlastSNP, collinter);
	      }
	    }
	    else{
	      result2.df=0;result2.sc=0;
	    }
	    if(qt){
	      df = result2.df - result1.df;
	      Fstat = ((result2.sc -result1.sc)/(df))/(result1.sc/result1.df);

	      if (df >= 1 && result2.sc > 0.000001){
		newValueCOLLREG = betai(result1.df/2,df/2,result1.df/(result1.df+df*Fstat));
	      }
	      else{
		newValueCOLLREG = 1;
	      }
	    }
	    else if(!qt){
	      tLog = 2*(result1.sc -result2.sc);
	      df = result1.df -result2.df;
	      if (df >= 1){
		  newValueCOLLREG = pValueCalc(df/2, tLog/2);
		}
	      else{
		newValueCOLLREG = 1;
	      }
	    }

	    //	    cout<<result1.betaNew_se[1]<<endl;

	    if(!optimalrare  && l<nwindowssinglebin){
	      COLLREG_beta[l]=result1.b[1];
	      COLLREG_se[l]=result1.betaNew_se[1];
	    }

	    nSnps=0;
	    freeResult(resultDummy1,0);
	    freeResult(resultDummy2,0);

	    if(optimalrare==0){
	      if(nsim==0){

		pCOLLREGvec[l]=newValueCOLLREG;

	      }
	      if(n!=0){
		if((wilsonpretest==0 && rareregpretest==0) || COLLREGpretestpassed[l]==1){
#pragma omp atomic
		  completedCOLLREG[l]++;
		  if(!(newValueCOLLREG*(1-EPS)>0)){
		    cout <<" BUG "<<newValueCOLLREG*(1-EPS)<<" "<<newValueCOLLREG<<endl;
		    exit(1);
		  }
		  if(newValueCOLLREG*(1-EPS)<COLLREGstats[l]){ // p-value by permutation
#pragma omp atomic
		    COLLREGcount[l]++;
#pragma omp atomic
		    countchange++;
		  }
		}
	      }
	      else if(n==0){
		if(collinter!=3 && collinter!=4){
		  COLLREGstats[l]=newValueCOLLREG;
		}
		rarefCOLLREG[l]=traref;
		//		MAF_Level_VT_COLLREG[l]=m;
	      }
	    }
	    if(optimalrare==1){
	      if(n!=0){
#pragma omp atomic
		completedCOLLREG[l]++;
		if((wilsonpretest==0 && rareregpretest==0) || COLLREGpretestpassed[l]==1){
		  if(newValueCOLLREG*(1-EPS)<COLLREGstats[l]){ // p-value by permutation
#pragma omp atomic
		    COLLREGcount[l]++;
#pragma omp atomic
		    countchange++;
		    break;
		  }
		}
	      }
	      if(n==0 && m==0){
		COLLREGstats[l]=newValueCOLLREG;
		rarefCOLLREG[l]=traref;
	      }
	      else if(n==0 && m>0){
		if(newValueCOLLREG*(1-EPS)<COLLREGstats[l]){ // p-value by permutation			   if()
		  COLLREGstats[l]=newValueCOLLREG;
		  rarefCOLLREG[l]=traref;
		  MAF_Level_VT_COLLREG[l]=m;
		}
	      }
	      else if(n!=0 && m==0){
		COLLREGpermstats[l]=newValueCOLLREG;
	      }
	      else if(n!=0 && m!=0){
		if(newValueCOLLREG*(1-EPS)<COLLREGpermstats[l]){ // p-value by permutation			   if()
		  COLLREGpermstats[l]=newValueCOLLREG;
		}
	      }
	    }
	  } // loop over levels
	    /*
	      if(optimalrare==1 && n!=0){
	      if(n!=0){
	      if(window[l].COLLREGpretestpassed==1){
	      if(COLLREGpermstats[l]*(1-EPS)<COLLREGstats[l]){ // p-value by permutation
	      #pragma omp atomic
	      COLLREGcount[l]++;
	      #pragma omp atomic
	      countchange++;
	      }
	      }
	      }
	      }
	    */
	    // PRETEST
	  if(n>0 && wilsonpretest!=0 && countchange>0 && COLLREGpretestpassed[l]){
#pragma omp atomic
	    countchange--;

	    double CI_low=wilson_lower(completedCOLLREG[l],COLLREGcount[l]);
	    if(CI_low>wilsonpretest){
	      pCOLLREGvec[l]=(double)COLLREGcount[l]/completedCOLLREG[l];
	      COLLREGpretestpassed[l]=0;
#pragma omp atomic
	      window[l].ntests-=1;
	      if(window[l].ntests==0){
#pragma omp atomic
		bins2test--;
	      }
	      COLLREG_CI[l][0]=CI_low;
	      COLLREG_CI[l][1]=wilson_upper(completedCOLLREG[l],COLLREGcount[l]);
	      window[l].nperm[5]=completedCOLLREG[l];
	    }
	  }
	  else if(rareregpretest!=0 && n==0){
	    if(COLLREGstats[l]*(1+EPS)>rareregpretest){
	      pCOLLREGvec[l]=COLLREGstats[l];
#pragma omp atomic
	      COLLREGpretestpassed[l]-=1;
#pragma omp atomic
	      window[l].ntests-=1;
	      if(window[l].ntests==0){
#pragma omp atomic
		bins2test--;
	      }
	    }
	  } // pretest
	}// if COLLREG
	// CMAT statistic
	if(CMATflag && ((wilsonpretest==0 && rareregpretest==0) || CMATpretestpassed[l]==1)){
	  int countchange=0;
	  stat=0;
	  for(m = 0; m < nRareLimits[l]; m++){
	    traref = (double)rareLimits[l][m];
	    if(numberOfAllCov>0){
	      limitsStatCMAT[l][m] = (double)calcCMATcov(cov, numberOfAllCov, sexcov, person, ncasesqc, ncontrolsqc, nrestqc, thread, nlinestfam, nlinestped, start, end, OR_CMAT, map, traref, nRareSNPs, BinSNPs, BinSNPsCCFlags_local, BinSNPsCovCatFlags ,nwordsSNPs, SNPMapInverse,weightvec, window, l, m, n, nCovCathegories);
	    }
	    else{
	      limitsStatCMAT[l][m] = (double)calcCMAT(cov, sexcov, person, ncasesqc, ncontrolsqc, nrestqc, thread, nlinestfam, nlinestped, start, end, OR_CMAT, map, traref, nRareSNPs, BinSNPs, BinSNPsCCFlags_local, nwordsSNPs, SNPMapInverse,weightvec, window, l, m, n);
	    }

	    if(std::isnan(limitsStatCMAT[l][m])==1){
	      limitsStatCMAT[l][m]=0;
	    }
	    if(n==0 && m==0){
	      OR_CMAT_vec[l]=OR_CMAT;
	      rarefCMAT[l]=rareLimits[l][0];
	    }
	    if(m==0){
	      stat=limitsStatCMAT[l][0];
	    }
	    else if(limitsStatCMAT[l][m]*(1+EPS)>stat){
	      stat=limitsStatCMAT[l][m];
	      if(n==0){
		OR_CMAT_vec[l]=OR_CMAT;
		rarefCMAT[l]=rareLimits[l][m];
		MAF_Level_VT_CMAT[l]=m;
	      }
	      else if(n>0 && stat*(1+EPS)>CMATstats[l]){
		break;
	      }
	    }
	  }
	  if(n==0){
	    CMATstats[l]=stat;
	  }
	  else{
	    CMATpermstats[l]=stat;
#pragma omp atomic
	    completedCMAT[l]++;
	    if((CMATpermstats[l]*(1+EPS)>CMATstats[l])){
#pragma omp atomic
	      CMATcount[l]++;
#pragma omp atomic
	      countchange++;
	    }
	  }
	  if(n>0 && wilsonpretest!=0 && CMATpretestpassed[l] && countchange>0){
	    //#pragma omp critical
	    //	      {
#pragma omp atomic
	    countchange-=1;
	    double CI_low=wilson_lower(completedCMAT[l],CMATcount[l]);
	    if(CI_low>wilsonpretest){
	      pCMATvec[l]=(double)CMATcount[l]/completedCMAT[l];
#pragma omp atomic
	      CMATpretestpassed[l]-=1;
#pragma omp atomic
	      window[l].ntests-=1;
	      if(window[l].ntests==0){
#pragma omp atomic
		bins2test--;
	      }
	      CMAT_CI[l][0]=CI_low;
	      CMAT_CI[l][1]=wilson_upper(completedCMAT[l],CMATcount[l]);
	      window[l].nperm[1]=completedCMAT[l];
	    }
	  }
	  else if(rareregpretest!=0 && n==0){
	    double pCMATprepre=pValueCalc(0.5,(double)stat/2);
	    if(pCMATprepre*(1+EPS)>rareregpretest){
#pragma omp atomic
	      CMATpretestpassed[l]-=1;
#pragma omp atomic
	      window[l].ntests-=1;
	      if(window[l].ntests==0){
#pragma omp atomic
		bins2test--;
	      }
	      pCMATvec[l]=pCMATprepre;
	    }
	  }
	}// if CMATflag

	if(COLLflag && ((wilsonpretest==0 && rareregpretest==0) || COLLpretestpassed[l])){
	  int countchange=0;
	  if(optimalrare==0){
	    if(!ifChrX) {

	      limitsStatCOLL[l][0] = (double)calcCOLL_bin(ncasesqc, ncontrolsqc, start, end,  map, raref, OR_COLL, BinSNPs, BinSNPsCCFlags_local, nwordsSNPs, SNPMapInverse, window, l, 0, n, nwindowssinglebin, collinter, doublewindowcoord);
	    }
	    else{
	      limitsStatCOLL[l][0] = (double)calcCOLL_X_bin(ncasesqc, ncontrolsqc, start, end, map, raref, OR_COLL, OR_COLL_f, BinSNPs, BinSNPsCCFlags_local, nwordsSNPs, SNPMapInverse,BinSNPsGenderFlags,ncasesqcMale, ncasesqcFemale,ncontrolsqcMale, ncontrolsqcFemale,window, l, 0, n);
	      OR_COLL_f_vec[l]=OR_COLL_f;
	    }
	    if(n==0){
	      OR_COLL_vec[l]=OR_COLL;
	      if(nsim==0){
		pCOLLvec[l]=pValueCalc(0.5, (double)limitsStatCOLL[l][0]/2);
		//		COLLpretestpassed[l]=0;
	      }
	      else{
		COLLstats[l]=limitsStatCOLL[l][0];
	      }
	    }
	    else{
	      COLLpermstats[l]=limitsStatCOLL[l][0];
	    }
	  }
	  else if(optimalrare==1){
	    for(m = 0; m < nRareLimits[l]; m++){
	      traref = (double)rareLimits[l][m];
	      if(!ifChrX) {
		limitsStatCOLL[l][m] = (double)calcCOLL_bin(ncasesqc, ncontrolsqc, start, end, map, traref, OR_COLL, BinSNPs, BinSNPsCCFlags_local, nwordsSNPs, SNPMapInverse,window, l, m, n, nwindowssinglebin, collinter, doublewindowcoord);
	      }
	      else {
		limitsStatCOLL[l][m] = (double)calcCOLL_X_bin(ncasesqc, ncontrolsqc, start, end, map, traref, OR_COLL, OR_COLL_f, BinSNPs, BinSNPsCCFlags_local, nwordsSNPs, SNPMapInverse,BinSNPsGenderFlags,ncasesqcMale, ncasesqcFemale,ncontrolsqcMale, ncontrolsqcFemale,window, l, m, n);
	      }
	      if(m==0){
		if(n==0){
		  OR_COLL_vec[l]=OR_COLL;
		  OR_COLL_f_vec[l]=OR_COLL_f;
		  COLLstats[l]=limitsStatCOLL[l][m];
		  rarefCOLL[l]=rareLimits[l][m];
		}
		else{
		  COLLpermstats[l]=limitsStatCOLL[l][m];
		}
	      }
	      else if(m!=0 && n==0){
		if(limitsStatCOLL[l][m]*(1+EPS)>COLLstats[l]){
		  COLLstats[l]=limitsStatCOLL[l][m];
		  rarefCOLL[l]=rareLimits[l][m];
		  MAF_Level_VT_COLL[l]=m;
		  OR_COLL_vec[l]=OR_COLL;
		  OR_COLL_f_vec[l]=OR_COLL_f;
		}
	      }
	      else if(m!=0 && n!=0){
		if(limitsStatCOLL[l][m]*(1+EPS)>COLLstats[l]){
		  COLLpermstats[l]=limitsStatCOLL[l][m];
		}
	      }
	    }
	  }

	  // RARE_PRE_PRETEST
	  if(n>0){
#pragma omp atomic
	    completedCOLL[l]++;
	    if((COLLpermstats[l])*(1+EPS)>COLLstats[l]){
#pragma omp atomic
	      COLLcount[l]++;
#pragma omp atomic
	      countchange++;
	    }
	  }

	  else if(n==0 && rareregpretest!=0){
	    double pCOLLprepre=pValueCalc(0.5,(double)COLLstats[l]/2);
	    if(pCOLLprepre*(1+EPS)>wilsonpretest){
	      pCOLLvec[l]=pCOLLprepre;
#pragma omp atomic
	      COLLpretestpassed[l]-=1;
#pragma omp atomic
	      window[l].ntests-=1;
	      if(window[l].ntests==0){
#pragma omp atomic
		bins2test--;
	      }
	    }
	  }
	  if(wilsonpretest!=0 && countchange>0 && COLLpretestpassed[l]){
	    //#pragma omp critical
	    //		  {
#pragma omp atomic
	    countchange-=1;
	    double CI_low=wilson_lower(completedCOLL[l],COLLcount[l]);
	    if(CI_low>wilsonpretest){
	      pCOLLvec[l]=(double)COLLcount[l]/completedCOLL[l];
	      window[l].nperm[2]=completedCOLL[l];
#pragma omp atomic
	      COLLpretestpassed[l]-=1;
#pragma omp atomic
	      window[l].ntests-=1;
	      if(window[l].ntests==0){
#pragma omp atomic
		bins2test--;
	      }
	      COLL_CI[l][0]=CI_low;
	      COLL_CI[l][1]=wilson_upper(completedCOLL[l],COLLcount[l]);
	    }
	  }
	} // if COLLflag
      } // if window tests > 0

    }// l loop
  if (rare_stratify) {
    delete[] wCounts;
  }
  free(snps);
  free(snps_second);
  free(snpsDom);
  free(environs);
  free(inters);
  free(snpsDom_second);
  free(environs_second);
  free(inters_second);
  return bins2test;
}




  void outRare(int binsizeRare, int nwindows, int nwindowssinglebin, int **windowPositions, struct PERSON *person, struct MAP *map, struct COUNTS **counts, bool optimalrare, bool FISHERflag, bool REGRESSIONflag, bool FRACREGflag,bool COLLREGflag, bool COLLflag, bool CMATflag, int *nRareSNPsCOLL, int *nRareSNPsCMAT, int *nRareSNPsFISHER, int *nRareSNPsREGRESSION, int *nRareSNPsFRACREG,int *nRareSNPsCOLLREG, int *MAF_Level_VT_FISHER, int *MAF_Level_VT_CMAT, int *MAF_Level_VT_COLL, int *MAF_Level_VT_REGRESSION, int *MAF_Level_VT_FRACREG, int *MAF_Level_VT_COLLREG, int nsim, string rarefile, string rareTopfile, int raretop, string intervalfile, string outputname, fstream &errorfile, fstream &logfile,  int rarepretest, float rarepretestlimit, double rareregpretest,  double *rarefFISHER, double *rarefREGRESSION, double *rarefFRACREG,double *rarefCOLLREG, double *rarefCOLL, double *OR_COLL_vec, double *OR_COLL_f_vec, double *rarefCMAT, double *OR_CMAT_vec, double *pCOLLvec, double *pCMATvec, double *pFISHERvec, double *pREGRESSIONvec, double *pFRACREGvec,double *pCOLLREGvec, double *pFISHERvecChi, int *nSNPsInWindow, double raref, int thread, int *FISHERnotconverged, int *FISHERpretestpassed, int *REGRESSIONpretestpassed, int *FRACREGpretestpassed, int *COLLREGpretestpassed, int *CMATpretestpassed, int *COLLpretestpassed, struct WINDOW *window, double **FISHER_CI, double **CMAT_CI, double **COLL_CI, double **REGRESSION_CI, double **FRACREG_CI, double **COLLREG_CI, double *COLLREG_beta, double *COLLREG_se, double *FRACREG_beta, double *FRACREG_se, int *FISHERcount, int *REGRESSIONcount, int *FRACREGcount,int *COLLREGcount, int *COLLcount, int *CMATcount, int merging, int flanking, int catIntervals, int minRareInBin, int expandIntervals, int mafadjust, int nlinestped, int binamin, int binamax, struct WEIGHTS weights, int verbose, string SetIDfile, int intervaleditor, int setid)
  {

    int start, end, j, k, i, l;
    int rarepos[nwindowssinglebin][2];
    int raremapinv[nlinestped];
    int rarecount=0;

    // Find relative positions of rare SNPs in bins
    for(j=0; j<nwindowssinglebin; j++){
      for(l=0; l<2; l++){
	rarepos[j][l]=0;
      }
    }

    for(j=0; j<nlinestped; j++){
      if(israre(j, map[j].mafr, raref, map)==1){
	rarecount++;
	raremapinv[j]=rarecount;
      }
      else{
	raremapinv[j]=-99;
      }
    }

    cout << "\nOutput of rare variant analysis is being generated..." << endl;




    for(j=0; j<nwindowssinglebin; j++){
      start=windowPositions[j][0];
      end=windowPositions[j][1];

      for(l=start; l<end+1; l++){
	if(israre(l,map[l].mafr, raref, map))
	  {
	    if(rarepos[j][0]==0){
	      rarepos[j][0]=raremapinv[l];
	    }
	    rarepos[j][1]=raremapinv[l];
	  }
      }
    }


    // get permuted p

    for(j=0; j<nwindowssinglebin; j++){
      if(FISHERflag && ((wilsonpretest==0 && rareregpretest==0) || FISHERpretestpassed[j]==1)){
	pFISHERvec[j]=(double)FISHERcount[j]/nsim;
	FISHER_CI[j][0]=wilson_lower(nsim,FISHERcount[j]);
	FISHER_CI[j][1]=wilson_upper(nsim,FISHERcount[j]);
	window[j].nperm[0]=nsim;
      }
      if(CMATflag && ((wilsonpretest==0 && rareregpretest==0) || CMATpretestpassed[j]==1)){
	pCMATvec[j]=(double)CMATcount[j]/nsim;
	CMAT_CI[j][0]=wilson_lower(nsim,CMATcount[j]);
	CMAT_CI[j][1]=wilson_upper(nsim,CMATcount[j]);
	window[j].nperm[1]=nsim;
      }
      if(COLLflag && ((wilsonpretest==0 && rareregpretest==0) || COLLpretestpassed[j]==1) && nsim!=0){
	pCOLLvec[j]=(double)COLLcount[j]/nsim;
	COLL_CI[j][0]=wilson_lower(nsim,COLLcount[j]);
	COLL_CI[j][1]=wilson_upper(nsim,COLLcount[j]);
	window[j].nperm[2]=nsim;
      }
      if(REGRESSIONflag && ((wilsonpretest==0 && rareregpretest==0) || REGRESSIONpretestpassed[j]==1) && nsim!=0){
	pREGRESSIONvec[j]=(double)REGRESSIONcount[j]/nsim;
	REGRESSION_CI[j][0]=wilson_lower(nsim,REGRESSIONcount[j]);
	REGRESSION_CI[j][1]=wilson_upper(nsim,REGRESSIONcount[j]);
	window[j].nperm[3]=nsim;
      }
      if(FRACREGflag && ((wilsonpretest==0 && rareregpretest==0) || FRACREGpretestpassed[j]==1) && nsim!=0){
	pFRACREGvec[j]=(double)FRACREGcount[j]/nsim;
	FRACREG_CI[j][0]=wilson_lower(nsim,FRACREGcount[j]);
	FRACREG_CI[j][1]=wilson_upper(nsim,FRACREGcount[j]);
	window[j].nperm[4]=nsim;
      }
      if(COLLREGflag && ((wilsonpretest==0 && rareregpretest==0) || COLLREGpretestpassed[j]==1) && nsim!=0){
	pCOLLREGvec[j]=(double)COLLREGcount[j]/nsim;
	COLLREG_CI[j][0]=wilson_lower(nsim,COLLREGcount[j]);
	COLLREG_CI[j][1]=wilson_upper(nsim,COLLREGcount[j]);
	window[j].nperm[5]=nsim;
      }
    }





    // Main output file

    fstream rare;
    rare.open(rarefile.c_str(), ios::out);


    rare << "#INFO:";
    if(wilsonpretest!=0){
      rare<<"ADAPTIVE="<<wilsonpretest<<"('+'=IF_PASSED);";
    }
    if(FISHERflag){
      for(l=0; l<nwindowssinglebin;l++){
	if(FISHERnotconverged[l]>0){
	  rare<<"**=DID_NOT_CONVERGE;";
	  break;
	}
      }
    }
    rare<<"SIMULATION="<<nsim<<";";

    rare << "MAFT="<<raref<<";";
    rare << "VT="<<optimalrare<<";";
    rare << "MAF_ADJUST="<<mafadjust<<";";

    if(intervalfile == " "){
      if(binsizeRare!=0){
	rare<<"BINSIZE_RARE="<<binsizeRare<<";";
      }
    }
    if(minRareInBin>0){
      rare<<"MIN_RARE_IN_BIN="<<minRareInBin<<";";
    }
    if(maxRareInBin>0 && maxRareInBin!=1000000){
      rare<<"MAX_RARE_IN_BIN="<<maxRareInBin<<";";
    }
    if(intervalfile!=" "){
      rare<<"INTERVAL_IN="<<intervalfile<<";";
      if(merging>0){
	rare<<"MERGE="<<merging<<";";
      }
      if(flanking>0){
	rare<<"FLANKING="<<flanking<<";";
      }
      if(catIntervals>0){
	rare<<"CONCATENATE_INTERVALS="<<catIntervals<<";";
      }
      if(expandIntervals>0){
	rare<<"CLOSE_GAPS="<<expandIntervals<<";";
      }
    }
    if(weights.mode!=0){
      if(weights.mode==1){
	rare<<"WEIGHTS=1/SD;";
      }
      else if(weights.mode==2){
	rare<<"WEIGHTS=BETA;"<<weights.betapar1<<";"<<weights.betapar2<<";";
      }
      else if(weights.mode==3){
	rare<<"WEIGHTS=LOGISTIC;"<<weights.logpar1<<";"<<weights.logpar2<<";";
      }
    }

#if PARALLELN
    int maxthreads;
    maxthreads=omp_get_max_threads();
    if(maxthreads>MAXTHREAD){
      maxthreads=MAXTHREAD;
    }
    rare<<"MAXTHREADS="<<maxthreads<<";";
#endif
#if !PARALLELN
    rare<<"MAXTHREADS=1;";
#endif

    rare<<"\n";
    if(verbose==2){rare << "BinNr\tSNPs\t";}

    rare<<"chr\tstart\tend";

    if(verbose==2){rare<<"\tfrom_SNP\tto_SNP\tfrom_rare_SNP_Nr\tto_rare_SNP_Nr";}
    //  rare<<"\tMAFlevel";

    if(optimalrare==0){
      rare << "\traref";
    }

    rare<<"\tnRV";

    //FR
    if(FISHERflag==1){
      if(optimalrare==1){
	rare<<"\trarefreqFR";
	rare<<"\tnRareFR";
      }
      rare<<"\tp_FR\tCI_FR";
      if(wilsonpretest){
	rare<<"\tnPerm_FR";
      }
      if(verbose==2){rare<<"\tp_FR_Corr";}
    }
    //CMAT
    if(CMATflag==1){
      if(optimalrare==1){
	rare<<"\trarefreqCMAT";
	rare<<"\tnRareCMAT";
      }
      rare<<"\tp_CMAT\tCI_CMAT";
      if(wilsonpretest){
	rare<<"\tnPerm_CMAT";
      }
      if(verbose==2){rare<<"\tp_CMAT_Corr";}
      rare<<"\tOR_CMAT";
    }
    //COLL
    if(COLLflag==1){
      if(optimalrare==1){
	rare<<"\trarefreqCOLL";
	rare<<"\tnRareCOLL";
      }
      rare<<"\tp_COLL";
      if(nsim!=0){
	rare<<"\tCI_COLL";
	if(wilsonpretest && nsim!=0){
	  rare<<"\tnPerm_COLL";
	}
      }
      if(verbose==2){rare<<"\tp_COLL_Corr";}
      rare<<"\tOR_COLL";
    }

    // REG
    if(REGRESSIONflag==1){
      if(optimalrare==1){
	rare<<"\trarefreqREG";
	rare<<"\tnRareREG";
      }
      rare<<"\tp_REG";
      if(nsim!=0){
	rare<<"\tCI_REG";
	if(wilsonpretest && nsim!=0){
	  rare<<"\tnPerm_REG";
	}
      }
      if(verbose==2){rare<<"\tp_REG_Corr";}
    }

    // FRACREG
    if(FRACREGflag==1){
      if(optimalrare==1){
	rare<<"\trarefreqFRACREG";
	rare<<"\tnRareFRACREG";
      }
      rare<<"\tp_FRACREG";
      if(nsim!=0){
	rare<<"\tCI_FRACREG";
	if(wilsonpretest && nsim!=0){
	  rare<<"\tnPerm_FRACREG";
	}
      }
      if(verbose==2){rare<<"\tp_FRACREG_Corr";}
      if(!optimalrare){
	rare<<"\tbeta_FRACREG\tse_FRACREG";
      }
    }

    // FRACREG
    if(COLLREGflag==1){
      if(optimalrare==1){
	rare<<"\trarefreqCOLLREG";
	rare<<"\tnRareCOLLREG";
      }
      rare<<"\tp_COLLREG";
      if(nsim!=0){
	rare<<"\tCI_COLLREG";
	if(wilsonpretest && nsim!=0){
	  rare<<"\tnPerm_COLLREG";
	}
      }
      if(verbose==2){rare<<"\tp_COLLREG_Corr";}
      if(!optimalrare){
	rare<<"\tbeta_COLLREG\tse_COLLREG";
      }
    }

    if(intervalfile!=" " && featurecol>=0){
      rare<<"\tFEATURE";
    }
    if(binamax!=0){
      rare<<"\tBINSIZE\tBINSTART";
    }

    rare<<"\n";

    start=0;
    end=0;
    k=0;

    // Bonferroni-corrected p-values
    float p_CMAT_Corr;
    float p_COLL_Corr;
    float p_FISHER_Corr;
    float p_REGRESSION_Corr;
    float p_FRACREG_Corr;
    float p_COLLREG_Corr;
    for(j=0; j<nwindowssinglebin; j++){
      start=windowPositions[j][0];
      end=windowPositions[j][1];
      if(verbose==2){
	if(FISHERflag){
	  /*	if(window[j].FISHERnotconverged==1){
		pFISHERvec[j]=1;
		p_FISHER_Corr = 1;
		}*/
	  if(pFISHERvec[j]==0){
	    p_FISHER_Corr = (1/(float)nsim)*(float)nwindowssinglebin;
	  }
	  else{
	    p_FISHER_Corr = (float)pFISHERvec[j]*nwindowssinglebin;
	  }
	}
	if(p_FISHER_Corr>1){
	  p_FISHER_Corr=1;
	}
	if(CMATflag){
	  if(pCMATvec[j]==0){
	    p_CMAT_Corr = (1/(float)nsim)*(float)nwindowssinglebin;
	  }
	  else{
	    p_CMAT_Corr = (float)pCMATvec[j]*nwindowssinglebin;
	  }
	  if(p_CMAT_Corr>1){
	    p_CMAT_Corr=1;
	  }

	}
	if(COLLflag){
	  if(pCOLLvec[j]==0){
	    p_COLL_Corr = (1/(float)nsim)*(float)nwindowssinglebin;
	  }
	  else{
	    p_COLL_Corr = (float)pCOLLvec[j]*nwindowssinglebin;
	  }
	  if(p_COLL_Corr>1){
	    p_COLL_Corr=1;
	  }
	}
	if(REGRESSIONflag){
	  if(pREGRESSIONvec[j]==0){
	    p_REGRESSION_Corr = (1/(float)nsim)*(float)nwindowssinglebin;
	  }
	  else{
	    p_REGRESSION_Corr = (float)pREGRESSIONvec[j]*nwindowssinglebin;
	  }
	  if(p_REGRESSION_Corr>1){
	    p_REGRESSION_Corr=1;
	  }
	}
	if(FRACREGflag){
	  if(pFRACREGvec[j]==0){
	    p_FRACREG_Corr = (1/(float)nsim)*(float)nwindowssinglebin;
	  }
	  else{
	    p_FRACREG_Corr = (float)pFRACREGvec[j]*nwindowssinglebin;
	  }
	  if(p_FRACREG_Corr>1){
	    p_FRACREG_Corr=1;
	  }
	}
	if(COLLREGflag){
	  if(pCOLLREGvec[j]==0){
	    p_COLLREG_Corr = (1/(float)nsim)*(float)nwindowssinglebin;
	  }
	  else{
	    p_COLLREG_Corr = (float)pCOLLREGvec[j]*nwindowssinglebin;
	  }
	  if(p_COLLREG_Corr>1){
	    p_COLLREG_Corr=1;
	  }
	}
      }

      // Data output
      if(verbose==2){rare << j+1 << "\t" << nSNPsInWindow[j]<< "\t" ;}
      rare << map[start].chr<< "\t"<< map[start].pos << "\t" << map[end].pos;
      if(verbose==2){rare<< "\t" <<map[start].rs << "\t" << map[end].rs <<"\t"<<rarepos[j][0]<<"\t"<<rarepos[j][1];}

      if(optimalrare==1){
	//      rare << "\t" << window[j].n_level;
      }

      if(optimalrare==0){
	rare << "\t" << window[j].maf_at_level[0];
      }

      rare << "\t" << nRareSNPs[j];
      if(optimalrare==1 && FISHERflag==1){
	rare << "\t" << rarefFISHER[j];
	rare<< "\t"  <<window[j].n_at_level[MAF_Level_VT_FISHER[j]];
      }

      if(FISHERflag==1){
	rare <<"\t"<<pFISHERvec[j];
	if(FISHERnotconverged[j]==1){rare<<"**";}
	else if(wilsonpretest!=0 && FISHERpretestpassed[j]==1){rare<<"+";}
	rare<<"\t["<<FISHER_CI[j][0]<<","<<FISHER_CI[j][1]<<"]";
	if(wilsonpretest){
	  rare<<"\t"<<window[j].nperm[0];
	}
	if(verbose==2){rare<< "\t"<< p_FISHER_Corr;}
      }
      if(optimalrare==1 && CMATflag==1){
	rare << "\t" <<rarefCMAT[j]<<"\t";
	rare<<window[j].n_at_level[MAF_Level_VT_CMAT[j]];
      }
      if(CMATflag==1){
	rare<< "\t" << pCMATvec[j];
	if(wilsonpretest!=0 && CMATpretestpassed[j]==1){
	  rare<<"+";
	}
	rare<<"\t["<<CMAT_CI[j][0]<<","<<CMAT_CI[j][1]<<"]";
	if(wilsonpretest){
	  rare<<"\t"<<window[j].nperm[1];
	}
	if(verbose==2){ rare<< "\t" << p_CMAT_Corr;}
	rare << "\t" << OR_CMAT_vec[j];
      }
      if(optimalrare==1 && COLLflag==1){
	rare  << "\t" << rarefCOLL[j];
	rare<< "\t"<< window[j].n_at_level[MAF_Level_VT_COLL[j]];
      }
      if(COLLflag==1){
	rare <<"\t" << pCOLLvec[j];
	if(wilsonpretest!=0  && COLLpretestpassed[j]==1){
	  rare<<"+";
	}
	if(nsim!=0){
	  rare<<"\t["<<COLL_CI[j][0]<<","<<COLL_CI[j][1]<<"]";
	  if(wilsonpretest){
	    rare<<"\t"<<window[j].nperm[2];
	  }
	}
	if(verbose==2){rare << "\t" << p_COLL_Corr;}
	rare << "\t" << OR_COLL_vec[j];
	if(strcmp(map[start].chr, "23")==0){rare << "/" << OR_COLL_f_vec[j];}
      }
      // REGRESSION TESTS
      if(optimalrare==1 && REGRESSIONflag==1){
	rare << "\t" <<rarefREGRESSION[j]<<"\t";
	rare<<window[j].n_at_level[MAF_Level_VT_REGRESSION[j]];
      }
      if(REGRESSIONflag==1){
	rare<< "\t" << pREGRESSIONvec[j];
	if(wilsonpretest!=0  && REGRESSIONpretestpassed[j]==1){
	  rare<<"+";
	}
	if(nsim!=0){
	  rare<<"\t["<<REGRESSION_CI[j][0]<<","<<REGRESSION_CI[j][1]<<"]";
	  if(wilsonpretest){
	    rare<<"\t"<<window[j].nperm[3];
	  }
	}
	if(verbose==2){rare<< "\t" << p_REGRESSION_Corr;}
      }

      if(optimalrare==1 && FRACREGflag==1){
	rare << "\t" <<rarefFRACREG[j]<<"\t";
	rare<<window[j].n_at_level[MAF_Level_VT_FRACREG[j]];
      }
      if(FRACREGflag==1){
	rare<< "\t" << pFRACREGvec[j];
	if(wilsonpretest!=0 && FRACREGpretestpassed[j]==1){
	  rare<<"+";
	}
	if(nsim!=0){
	  rare<<"\t["<<FRACREG_CI[j][0]<<","<<FRACREG_CI[j][1]<<"]";
	  if(wilsonpretest){
	    rare<<"\t"<<window[j].nperm[4];
	  }
	}
	if(verbose==2){ rare<< "\t" << p_FRACREG_Corr;}
	if(!optimalrare){
	  rare<<"\t"<<FRACREG_beta[j]<<"\t"<<FRACREG_se[j];
	}
      }

      if(optimalrare==1 && COLLREGflag==1){
	rare << "\t" <<rarefCOLLREG[j]<<"\t";
	rare<<window[j].n_at_level[MAF_Level_VT_COLLREG[j]];
      }

      if(COLLREGflag==1){
	rare<< "\t" << pCOLLREGvec[j];
	if(wilsonpretest!=0 && COLLREGpretestpassed[j]==1){
	  rare<<"+";
	}
	if(nsim!=0){
	  rare<<"\t["<< COLLREG_CI[j][0]<<","<<COLLREG_CI[j][1]<<"]";
	  if(wilsonpretest){
	    rare<<"\t"<<window[j].nperm[5];
	  }
	}
	if(verbose==2){rare<< "\t" << p_COLLREG_Corr;}
	if(!optimalrare){
	  rare<<"\t"<<COLLREG_beta[j]<<"\t"<<COLLREG_se[j];
	}


      }
      if(intervalfile!=" " && featurecol>=0){
	rare<< "\t" << intervals[j].feature;
      }
      if(binamax!=0){
	rare<< "\t" << intervals[j].s1<< "\t" << intervals[j].s2;
      }

      rare << "\n";
    }

    rare.close();
    cout << "Rare variant analysis is finished!"<<endl;
    logfile << "Rare variant analysis is finished!"<<endl;
  }

void outRareInter(int binsizeRare, int nwindows, int nwindowssinglebin, int **windowPositions, struct PERSON *person, struct MAP *map, struct COUNTS **counts, bool optimalrare, bool FISHERflag, bool REGRESSIONflag, bool FRACREGflag,bool COLLREGflag, bool COLLflag, bool CMATflag, int *nRareSNPsCOLL, int *nRareSNPsCMAT, int *nRareSNPsFISHER, int *nRareSNPsREGRESSION, int *nRareSNPsFRACREG,int *nRareSNPsCOLLREG, int *MAF_Level_VT_FISHER, int *MAF_Level_VT_CMAT, int *MAF_Level_VT_COLL, int *MAF_Level_VT_REGRESSION, int *MAF_Level_VT_FRACREG, int *MAF_Level_VT_COLLREG, int nsim, string rarefileinter, string rareTopfile, int raretop, string intervalfile, string outputname, fstream &errorfile, fstream &logfile,  int rarepretest, float rarepretestlimit, double rareregpretest,  double *rarefFISHER, double *rarefREGRESSION, double *rarefFRACREG,double *rarefCOLLREG, double *rarefCOLL, double *OR_COLL_vec, double *OR_COLL_f_vec, double *rarefCMAT, double *OR_CMAT_vec, double *pCOLLvec, double *pCMATvec, double *pFISHERvec, double *pREGRESSIONvec, double *pFRACREGvec,double *pCOLLREGvec, double *pFISHERvecChi, int *nSNPsInWindow, double raref, int thread, int *FISHERnotconverged, int *FISHERpretestpassed, int *REGRESSIONpretestpassed, int *FRACREGpretestpassed, int *COLLREGpretestpassed, int *CMATpretestpassed, int *COLLpretestpassed, struct WINDOW *window, double **FISHER_CI, double **CMAT_CI, double **COLL_CI, double **REGRESSION_CI, double **FRACREG_CI, double **COLLREG_CI, double *COLLREG_beta, double *COLLREG_se, double *FRACREG_beta, double *FRACREG_se, int *FISHERcount, int *REGRESSIONcount, int *FRACREGcount,int *COLLREGcount, int *COLLcount, int *CMATcount, int merging, int flanking, int catIntervals, int minRareInBin, int expandIntervals, int mafadjust, int nlinestped, int binamin, int binamax, struct WEIGHTS weights, int verbose, string SetIDfile, int intervaleditor, int setid, int collinter)
  {

    int start, end, start1, end1, start2, end2, j, k, i, l;
    int rarepos[nwindowssinglebin][2];
    int raremapinv[nlinestped];
    int rarecount=0;

    // Find relative positions of rare SNPs in bins
    for(j=0; j<nwindowssinglebin; j++){
      for(l=0; l<2; l++){
	rarepos[j][l]=0;
      }
    }

    for(j=0; j<nlinestped; j++){
      if(israre(j, map[j].mafr, raref, map)==1){
	rarecount++;
	raremapinv[j]=rarecount;
      }
      else{
	raremapinv[j]=-99;
      }
    }

    cout << "\nOutput of rare variant interaction analysis is being generated..." << endl;

    for(j=0; j<nwindowssinglebin; j++){

      start=windowPositions[j][0];
      end=windowPositions[j][1];

      for(l=start; l<end+1; l++){
	if(israre(l,map[l].mafr, raref, map)){
	  if(rarepos[j][0]==0){
	    rarepos[j][0]=raremapinv[l];
	  }
	  rarepos[j][1]=raremapinv[l];
	}
      }
    }
    // get permuted p

    for(j=nwindowssinglebin; j<nwindows; j++){
      if(FISHERflag && ((wilsonpretest==0 && rareregpretest==0) || FISHERpretestpassed[j]==1)){
	pFISHERvec[j]=(double)FISHERcount[j]/nsim;
	FISHER_CI[j][0]=wilson_lower(nsim,FISHERcount[j]);
	FISHER_CI[j][1]=wilson_upper(nsim,FISHERcount[j]);
	window[j].nperm[0]=nsim;
      }
      if(CMATflag && ((wilsonpretest==0 && rareregpretest==0) || CMATpretestpassed[j]==1)){
	pCMATvec[j]=(double)CMATcount[j]/nsim;
	CMAT_CI[j][0]=wilson_lower(nsim,CMATcount[j]);
	CMAT_CI[j][1]=wilson_upper(nsim,CMATcount[j]);
	window[j].nperm[1]=nsim;
      }
      if(COLLflag && ((wilsonpretest==0 && rareregpretest==0) || COLLpretestpassed[j]==1) && nsim!=0){
	pCOLLvec[j]=(double)COLLcount[j]/nsim;
	COLL_CI[j][0]=wilson_lower(nsim,COLLcount[j]);
	COLL_CI[j][1]=wilson_upper(nsim,COLLcount[j]);
	window[j].nperm[2]=nsim;
      }
      if(REGRESSIONflag && ((wilsonpretest==0 && rareregpretest==0) || REGRESSIONpretestpassed[j]==1) && nsim!=0){
	pREGRESSIONvec[j]=(double)REGRESSIONcount[j]/nsim;
	REGRESSION_CI[j][0]=wilson_lower(nsim,REGRESSIONcount[j]);
	REGRESSION_CI[j][1]=wilson_upper(nsim,REGRESSIONcount[j]);
	window[j].nperm[3]=nsim;
      }
      if(FRACREGflag && ((wilsonpretest==0 && rareregpretest==0) || FRACREGpretestpassed[j]==1) && nsim!=0){
	pFRACREGvec[j]=(double)FRACREGcount[j]/nsim;
	FRACREG_CI[j][0]=wilson_lower(nsim,FRACREGcount[j]);
	FRACREG_CI[j][1]=wilson_upper(nsim,FRACREGcount[j]);
	window[j].nperm[4]=nsim;
      }
      if(COLLREGflag && ((wilsonpretest==0 && rareregpretest==0) || COLLREGpretestpassed[j]==1) && nsim!=0){
	pCOLLREGvec[j]=(double)COLLREGcount[j]/nsim;
	COLLREG_CI[j][0]=wilson_lower(nsim,COLLREGcount[j]);
	COLLREG_CI[j][1]=wilson_upper(nsim,COLLREGcount[j]);
	window[j].nperm[5]=nsim;
      }
    }

    // Main output file

    fstream rareinter;
    rareinter.open(rarefileinter.c_str(), ios::out);


    rareinter << "#INFO:";
    if(wilsonpretest!=0){
      rareinter<<"ADAPTIVE="<<wilsonpretest<<"('+'=IF_PASSED);";
    }
    if(FISHERflag){
      for(j=nwindowssinglebin; j<nwindows; j++){
	if(FISHERnotconverged[l]>0){
	  rareinter<<"**=DID_NOT_CONVERGE;";
	  break;
	}
      }
    }
    rareinter<<"SIMULATION="<<nsim<<";";

    rareinter << "MAFT="<<raref<<";";
    rareinter << "VT="<<optimalrare<<";";
    rareinter << "MAF_ADJUST="<<mafadjust<<";";

    if(intervalfile == " "){
      if(binsizeRare!=0){
	rareinter<<"BINSIZE_RARE="<<binsizeRare<<";";
      }
    }
    if(minRareInBin>0){
      rareinter<<"MIN_RARE_IN_BIN="<<minRareInBin<<";";
    }
    if(maxRareInBin>0 && maxRareInBin!=1000000){
      rareinter<<"MAX_RARE_IN_BIN="<<maxRareInBin<<";";
    }
    if(intervalfile!=" "){
      rareinter<<"INTERVAL_IN="<<intervalfile<<";";
      if(merging>0){
	rareinter<<"MERGE="<<merging<<";";
      }
      if(flanking>0){
	rareinter<<"FLANKING="<<flanking<<";";
      }
      if(catIntervals>0){
	rareinter<<"CONCATENATE_INTERVALS="<<catIntervals<<";";
      }
      if(expandIntervals>0){
	rareinter<<"CLOSE_GAPS="<<expandIntervals<<";";
      }
    }
    if(weights.mode!=0){
      if(weights.mode==1){
	rareinter<<"WEIGHTS=1/SD;";
      }
      else if(weights.mode==2){
	rareinter<<"WEIGHTS=BETA;"<<weights.betapar1<<";"<<weights.betapar2<<";";
      }
      else if(weights.mode==3){
	rareinter<<"WEIGHTS=LOGISTIC;"<<weights.logpar1<<";"<<weights.logpar2<<";";
      }
    }

#if PARALLELN
    int maxthreads;
    maxthreads=omp_get_max_threads();
    if(maxthreads>MAXTHREAD){
      maxthreads=MAXTHREAD;
    }
    rareinter<<"MAXTHREADS="<<maxthreads<<";";
#endif
#if !PARALLELN
    rareinter<<"MAXTHREADS=1;";
#endif

    rareinter<<"\n";

    rareinter << "chr1\tstart1\tend1\tchr2\tstart2\tend2";

    if(optimalrare==0){
      rareinter << "\traref";
    }

    rareinter<<"\tnRV";

    //FR
    if(FISHERflag==1){
      if(optimalrare==1){
	rareinter<<"\trarefreqFR";
	rareinter<<"\tnRareFR";
      }
      rareinter<<"\tp_FR\tCI_FR";
      if(wilsonpretest){
	rareinter<<"\tnPerm_FR";
      }
      if(verbose==2){rareinter<<"\tp_FR_Corr";}
    }
    //CMAT
    if(CMATflag==1){
      if(optimalrare==1){
	rareinter<<"\trarefreqCMAT";
	rareinter<<"\tnRareCMAT";
      }
      rareinter<<"\tp_CMAT\tCI_CMAT";
      if(wilsonpretest){
	rareinter<<"\tnPerm_CMAT";
      }
      if(verbose==2){rareinter<<"\tp_CMAT_Corr";}
      rareinter<<"\tOR_CMAT";
    }
    //COLL
    if(COLLflag==1){
      if(optimalrare==1){
	rareinter<<"\trarefreqCOLL";
	rareinter<<"\tnRareCOLL";
      }
      rareinter<<"\tp_COLL";
      if(nsim!=0){
	rareinter<<"\tCI_COLL";
	if(wilsonpretest && nsim!=0){
	  rareinter<<"\tnPerm_COLL";
	}
      }
      if(verbose==2){rareinter<<"\tp_COLL_Corr";}
      rareinter<<"\tOR_COLL";
    }

    // REG
    if(REGRESSIONflag==1){
      if(optimalrare==1){
	rareinter<<"\trarefreqREG";
	rareinter<<"\tnRareREG";
      }
      rareinter<<"\tp_REG";
      if(nsim!=0){
	rareinter<<"\tCI_REG";
	if(wilsonpretest && nsim!=0){
	  rareinter<<"\tnPerm_REG";
	}
      }
      if(verbose==2){rareinter<<"\tp_REG_Corr";}
    }

    // FRACREG
    if(FRACREGflag==1){
      if(optimalrare==1){
	rareinter<<"\trarefreqFRACREG";
	rareinter<<"\tnRareFRACREG";
      }
      rareinter<<"\tp_FRACREG";
      if(nsim!=0){
	rareinter<<"\tCI_FRACREG";
	if(wilsonpretest && nsim!=0){
	  rareinter<<"\tnPerm_FRACREG";
	}
      }
      if(verbose==2){rareinter<<"\tp_FRACREG_Corr";}
    }

    // FRACREG
    if(COLLREGflag==1){
      if(optimalrare==1){
	rareinter<<"\trarefreqCOLLREG";
	rareinter<<"\tnRareCOLLREG";
      }
      rareinter<<"\tp_COLLREG";
      if(nsim!=0){
	rareinter<<"\tCI_COLLREG";
	if(wilsonpretest && nsim!=0){
	  rareinter<<"\tnPerm_COLLREG";
	}
      }
      if(verbose==2){rareinter<<"\tp_COLLREG_Corr";}
    }

    if(intervalfile!=" " && featurecol>=0){
      rareinter<<"\tFEATURE1\tFEATURE2";
    }
    if(binamax!=0){
      rareinter<<"\tBINSIZE\tBINSTART";
    }

    rareinter<<"\n";

    start=0;
    end=0;
    k=0;

    // Bonferroni-corrected p-values
    float p_CMAT_Corr;
    float p_COLL_Corr;
    float p_FISHER_Corr;
    float p_REGRESSION_Corr;
    float p_FRACREG_Corr;
    float p_COLLREG_Corr;
    for(j=nwindowssinglebin; j<nwindows; j++){
      if(nRareSNPs[doublewindowcoord[j][0]]==0 || nRareSNPs[doublewindowcoord[j][1]]==0) continue;
      start1=windowPositions[doublewindowcoord[j][0]][0];
      end1=windowPositions[doublewindowcoord[j][0]][1];
      start2=windowPositions[doublewindowcoord[j][1]][0];
      end2=windowPositions[doublewindowcoord[j][1]][1];
      if(verbose==2){
	if(FISHERflag){
	  if(pFISHERvec[j]==0){
	    p_FISHER_Corr = (1/(float)nsim)*(float)(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	  else{
	    p_FISHER_Corr = (float)pFISHERvec[j]*(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	}
	if(CMATflag){
	  if(pCMATvec[j]==0){
	    p_CMAT_Corr = (1/(float)nsim)*(float)(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	  else{
	    p_CMAT_Corr = (float)pCMATvec[j]*(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	}
	if(COLLflag){
	  if(pCOLLvec[j]==0){
	    p_COLL_Corr = (1/(float)nsim)*(float)(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	  else{
	    p_COLL_Corr = (float)pCOLLvec[j]*(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	}
	if(REGRESSIONflag){
	  if(pREGRESSIONvec[j]==0){
	    p_REGRESSION_Corr = (1/(float)nsim)*(float)(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	  else{
	    p_REGRESSION_Corr = (float)pREGRESSIONvec[j]*(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	}
	if(FRACREGflag){
	  if(pFRACREGvec[j]==0){
	    p_FRACREG_Corr = (1/(float)nsim)*(float)(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	  else{
	    p_FRACREG_Corr = (float)pFRACREGvec[j]*(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	}
	if(COLLREGflag){
	  if(pCOLLREGvec[j]==0){
	    p_COLLREG_Corr = (1/(float)nsim)*(float)(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	  else{
	    p_COLLREG_Corr = (float)pCOLLREGvec[j]*(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	}
      }

      // Data output
      rareinter << map[start1].chr<< "\t"<< map[start1].pos << "\t" << map[end1].pos<<"\t" << map[start2].chr<< "\t"<< map[start2].pos << "\t" << map[end2].pos;

      if(optimalrare==0){
	double maf_at_level_tmp;
	if(collinter==3 || collinter==4){
	  if((window[doublewindowcoord[j][0]].maf_at_level[0]-window[doublewindowcoord[j][1]].maf_at_level[0])>EPS){
	    maf_at_level_tmp=window[doublewindowcoord[j][0]].maf_at_level[0];
	  }
	  else{
	    maf_at_level_tmp=window[doublewindowcoord[j][1]].maf_at_level[0];
	  }
	}
	else{
	  maf_at_level_tmp=window[j].maf_at_level[0];
	}
	rareinter << "\t" << maf_at_level_tmp;
      }
      if(collinter==3 || collinter==4){
	rareinter << "\t" << nRareSNPs[doublewindowcoord[j][0]]+nRareSNPs[doublewindowcoord[j][1]];
      }
      else{
	rareinter << "\t" << nRareSNPs[j];
      }
      if(optimalrare==1 && FISHERflag==1){
	rareinter << "\t" << rarefFISHER[j];
	rareinter << "\t"  <<window[j].n_at_level[MAF_Level_VT_FISHER[j]];
      }

      if(FISHERflag==1){
	rareinter <<"\t"<<pFISHERvec[j];
	if(FISHERnotconverged[j]==1){rareinter<<"**";}
	else if(wilsonpretest!=0 && FISHERpretestpassed[j]==1){rareinter<<"+";}
	rareinter<<"\t["<<FISHER_CI[j][0]<<","<<FISHER_CI[j][1]<<"]";
	if(wilsonpretest){
	  rareinter<<"\t"<<window[j].nperm[0];
	}
	if(verbose==2){rareinter<< "\t"<< p_FISHER_Corr;}
      }
      if(optimalrare==1 && CMATflag==1){
	rareinter << "\t" <<rarefCMAT[j]<<"\t";
	rareinter<<window[j].n_at_level[MAF_Level_VT_CMAT[j]];
      }
      if(CMATflag==1){
	rareinter<< "\t" << pCMATvec[j];
	if(wilsonpretest!=0 && CMATpretestpassed[j]==1){
	  rareinter<<"+";
	}
	rareinter<<"\t["<<CMAT_CI[j][0]<<","<<CMAT_CI[j][1]<<"]";
	if(wilsonpretest){
	  rareinter<<"\t"<<window[j].nperm[1];
	}
	if(verbose==2){ rareinter<< "\t" << p_CMAT_Corr;}
	rareinter << "\t" << OR_CMAT_vec[j];
      }
      if(optimalrare==1 && COLLflag==1){
	rareinter  << "\t" << rarefCOLL[j];
	rareinter<< "\t"<< window[j].n_at_level[MAF_Level_VT_COLL[j]];
      }
      if(COLLflag==1){
	rareinter <<"\t" << pCOLLvec[j];
	if(wilsonpretest!=0 && COLLpretestpassed[j]==1){
	  rareinter<<"+";
	}
	if(nsim!=0){
	  rareinter<<"\t["<<COLL_CI[j][0]<<","<<COLL_CI[j][1]<<"]";
	  if(wilsonpretest){
	    rareinter<<"\t"<<window[j].nperm[2];
	  }
	}
	if(verbose==2){rareinter << "\t" << p_COLL_Corr;}
	rareinter << "\t" << OR_COLL_vec[j];
	if(strcmp(map[start].chr, "23")==0){rareinter << "/" << OR_COLL_f_vec[j];}
      }
      // REGRESSION TESTS
      if(optimalrare==1 && REGRESSIONflag==1){
	rareinter << "\t" <<rarefREGRESSION[j]<<"\t";
	rareinter<<window[j].n_at_level[MAF_Level_VT_REGRESSION[j]];
      }
      if(REGRESSIONflag==1){
	rareinter<< "\t" << pREGRESSIONvec[j];
	if(wilsonpretest!=0 && REGRESSIONpretestpassed[j]==1){
	  rareinter<<"+";
	}
	if(nsim!=0){
	  rareinter<<"\t["<<REGRESSION_CI[j][0]<<","<<REGRESSION_CI[j][1]<<"]";
	  if(wilsonpretest){
	    rareinter<<"\t"<<window[j].nperm[3];
	  }
	}
	if(verbose==2){rareinter<< "\t" << p_REGRESSION_Corr;}
      }

      if(optimalrare==1 && FRACREGflag==1){
	rareinter << "\t" <<rarefFRACREG[j]<<"\t";
	rareinter<<window[j].n_at_level[MAF_Level_VT_FRACREG[j]];
      }
      if(FRACREGflag==1){
	rareinter<< "\t" << pFRACREGvec[j];
	if(wilsonpretest!=0 && FRACREGpretestpassed[j]==1){
	  rareinter<<"+";
	}
	if(nsim!=0){
	  rareinter<<"\t["<<FRACREG_CI[j][0]<<","<<FRACREG_CI[j][1]<<"]";
	  if(wilsonpretest){
	    rareinter<<"\t"<<window[j].nperm[4];
	  }
	}
	if(verbose==2){ rareinter<< "\t" << p_FRACREG_Corr;}
      }

      if(optimalrare==1 && COLLREGflag==1){
	rareinter << "\t" <<rarefCOLLREG[j]<<"\t";
	rareinter<<window[j].n_at_level[MAF_Level_VT_COLLREG[j]];
      }

      if(COLLREGflag==1){
	rareinter<< "\t" << pCOLLREGvec[j];
	if(wilsonpretest!=0 && COLLREGpretestpassed[j]==1){
	  rareinter<<"+";
	}
	if(nsim!=0){
	  rareinter<<"\t["<< COLLREG_CI[j][0]<<","<<COLLREG_CI[j][1]<<"]";
	  if(wilsonpretest){
	    rareinter<<"\t"<<window[j].nperm[5];
	  }
	}
	if(verbose==2){rareinter<< "\t" << p_COLLREG_Corr;}
      }
      if(intervalfile!=" " && featurecol>=0){
	rareinter<< "\t" << intervals[doublewindowcoord[j][0]].feature;
	rareinter<< "\t" << intervals[doublewindowcoord[j][1]].feature;
      }

      if(binamax!=0){
	rareinter<< "\t" << intervals[j].s1<< "\t" << intervals[j].s2;
      }

      rareinter << "\n";
    }

    rareinter.close();
    cout << "Rare variant interaction analysis is finished!"<<endl;
    logfile << "Rare variant interaction analysis is finished!"<<endl;
  }

// COLL_INTER 3,4 are conducted without permutation and output can be written on-the fly
  void outRareInterLowMem(int binsizeRare, int nwindows, int nwindowssinglebin, int **windowPositions, struct PERSON *person, struct MAP *map, struct COUNTS **counts, bool optimalrare, bool FISHERflag, bool REGRESSIONflag, bool FRACREGflag,bool COLLREGflag, bool COLLflag, bool CMATflag, int *nRareSNPsCOLL, int *nRareSNPsCMAT, int *nRareSNPsFISHER, int *nRareSNPsREGRESSION, int *nRareSNPsFRACREG,int *nRareSNPsCOLLREG, int *MAF_Level_VT_FISHER, int *MAF_Level_VT_CMAT, int *MAF_Level_VT_COLL, int *MAF_Level_VT_REGRESSION, int *MAF_Level_VT_FRACREG, int *MAF_Level_VT_COLLREG, int nsim, string rarefileinter, string rareTopfile, int raretop, string intervalfile, string outputname, fstream &errorfile, fstream &logfile,  int rarepretest, float rarepretestlimit, double rareregpretest,  double *rarefFISHER, double *rarefREGRESSION, double *rarefFRACREG,double *rarefCOLLREG, double *rarefCOLL, double *OR_COLL_vec, double *OR_COLL_f_vec, double *rarefCMAT, double *OR_CMAT_vec, double *pCOLLvec, double *pCMATvec, double *pFISHERvec, double *pREGRESSIONvec, double *pFRACREGvec,double *pCOLLREGvec, double *pFISHERvecChi, int *nSNPsInWindow, double raref, int thread, int *FISHERnotconverged, int *FISHERpretestpassed, int *REGRESSIONpretestpassed, int *FRACREGpretestpassed, int *COLLREGpretestpassed, int *CMATpretestpassed, int *COLLpretestpassed, struct WINDOW *window, double **FISHER_CI, double **CMAT_CI, double **COLL_CI, double **REGRESSION_CI, double **FRACREG_CI, double **COLLREG_CI, double *COLLREG_beta, double *COLLREG_se, double *FRACREG_beta, double *FRACREG_se, int *FISHERcount, int *REGRESSIONcount, int *FRACREGcount,int *COLLREGcount, int *COLLcount, int *CMATcount, int merging, int flanking, int catIntervals, int minRareInBin, int expandIntervals, int mafadjust, int nlinestped, int binamin, int binamax, struct WEIGHTS weights, int verbose, string SetIDfile, int intervaleditor, int setid)
  {

    int start, end, start1, end1, start2, end2, j, k, i, l;
    int rarepos[nwindowssinglebin][2];
    int raremapinv[nlinestped];
    int rarecount=0;

    // Find relative positions of rare SNPs in bins
    for(j=0; j<nwindowssinglebin; j++){
      for(l=0; l<2; l++){
	rarepos[j][l]=0;
      }
    }

    for(j=0; j<nlinestped; j++){
      if(israre(j, map[j].mafr, raref, map)==1){
	rarecount++;
	raremapinv[j]=rarecount;
      }
      else{
	raremapinv[j]=-99;
      }
    }

    cout << "\nOutput of rare variant interaction analysis is being generated..." << endl;

    for(j=0; j<nwindowssinglebin; j++){

      start=windowPositions[j][0];
      end=windowPositions[j][1];

      for(l=start; l<end+1; l++){
	if(israre(l,map[l].mafr, raref, map)){
	  if(rarepos[j][0]==0){
	    rarepos[j][0]=raremapinv[l];
	  }
	  rarepos[j][1]=raremapinv[l];
	}
      }
    }
    // get permuted p

    for(j=nwindowssinglebin; j<nwindows; j++){
      if(FISHERflag && ((wilsonpretest==0 && rareregpretest==0) || FISHERpretestpassed[j]==1)){
	pFISHERvec[j]=(double)FISHERcount[j]/nsim;
	FISHER_CI[j][0]=wilson_lower(nsim,FISHERcount[j]);
	FISHER_CI[j][1]=wilson_upper(nsim,FISHERcount[j]);
	window[j].nperm[0]=nsim;
      }
      if(CMATflag && ((wilsonpretest==0 && rareregpretest==0) || CMATpretestpassed[j]==1)){
	pCMATvec[j]=(double)CMATcount[j]/nsim;
	CMAT_CI[j][0]=wilson_lower(nsim,CMATcount[j]);
	CMAT_CI[j][1]=wilson_upper(nsim,CMATcount[j]);
	window[j].nperm[1]=nsim;
      }
      if(COLLflag && ((wilsonpretest==0 && rareregpretest==0) || COLLpretestpassed[j]==1) && nsim!=0){
	pCOLLvec[j]=(double)COLLcount[j]/nsim;
	COLL_CI[j][0]=wilson_lower(nsim,COLLcount[j]);
	COLL_CI[j][1]=wilson_upper(nsim,COLLcount[j]);
	window[j].nperm[2]=nsim;
      }
      if(REGRESSIONflag && ((wilsonpretest==0 && rareregpretest==0) || REGRESSIONpretestpassed[j]==1) && nsim!=0){
	pREGRESSIONvec[j]=(double)REGRESSIONcount[j]/nsim;
	REGRESSION_CI[j][0]=wilson_lower(nsim,REGRESSIONcount[j]);
	REGRESSION_CI[j][1]=wilson_upper(nsim,REGRESSIONcount[j]);
	window[j].nperm[3]=nsim;
      }
      if(FRACREGflag && ((wilsonpretest==0 && rareregpretest==0) || FRACREGpretestpassed[j]==1) && nsim!=0){
	pFRACREGvec[j]=(double)FRACREGcount[j]/nsim;
	FRACREG_CI[j][0]=wilson_lower(nsim,FRACREGcount[j]);
	FRACREG_CI[j][1]=wilson_upper(nsim,FRACREGcount[j]);
	window[j].nperm[4]=nsim;
      }
      if(COLLREGflag && ((wilsonpretest==0 && rareregpretest==0) || COLLREGpretestpassed[j]==1) && nsim!=0){
	pCOLLREGvec[j]=(double)COLLREGcount[j]/nsim;
	COLLREG_CI[j][0]=wilson_lower(nsim,COLLREGcount[j]);
	COLLREG_CI[j][1]=wilson_upper(nsim,COLLREGcount[j]);
	window[j].nperm[5]=nsim;
      }
    }

    // Main output file

    fstream rareinter;
    rareinter.open(rarefileinter.c_str(), ios::out);


    rareinter << "#INFO:";
    if(wilsonpretest!=0){
      rareinter<<"ADAPTIVE="<<wilsonpretest<<"('+'=IF_PASSED);";
    }
    for(j=nwindowssinglebin; j<nwindows; j++){
      if(FISHERnotconverged[j]>0){
	rareinter<<"**=DID_NOT_CONVERGE;";
	break;
      }
    }
    rareinter<<"SIMULATION="<<nsim<<";";

    rareinter << "MAFT="<<raref<<";";
    rareinter << "VT="<<optimalrare<<";";
    rareinter << "MAF_ADJUST="<<mafadjust<<";";

    if(intervalfile == " "){
      if(binsizeRare!=0){
	rareinter<<"BINSIZE_RARE="<<binsizeRare<<";";
      }
    }
    if(minRareInBin>0){
      rareinter<<"MIN_RARE_IN_BIN="<<minRareInBin<<";";
    }
    if(maxRareInBin>0 && maxRareInBin!=1000000){
      rareinter<<"MAX_RARE_IN_BIN="<<maxRareInBin<<";";
    }
    if(intervalfile!=" "){
      rareinter<<"INTERVAL_IN="<<intervalfile<<";";
      if(merging>0){
	rareinter<<"MERGE="<<merging<<";";
      }
      if(flanking>0){
	rareinter<<"FLANKING="<<flanking<<";";
      }
      if(catIntervals>0){
	rareinter<<"CONCATENATE_INTERVALS="<<catIntervals<<";";
      }
      if(expandIntervals>0){
	rareinter<<"CLOSE_GAPS="<<expandIntervals<<";";
      }
    }
    if(weights.mode!=0){
      if(weights.mode==1){
	rareinter<<"WEIGHTS=1/SD;";
      }
      else if(weights.mode==2){
	rareinter<<"WEIGHTS=BETA;"<<weights.betapar1<<";"<<weights.betapar2<<";";
      }
      else if(weights.mode==3){
	rareinter<<"WEIGHTS=LOGISTIC;"<<weights.logpar1<<";"<<weights.logpar2<<";";
      }
    }

#if PARALLELN
    int maxthreads;
    maxthreads=omp_get_max_threads();
    if(maxthreads>MAXTHREAD){
      maxthreads=MAXTHREAD;
    }
    rareinter<<"MAXTHREADS="<<maxthreads<<";";
#endif
#if !PARALLELN
    rareinter<<"MAXTHREADS=1;";
#endif

    rareinter<<"\n";

    rareinter << "chr1\tstart1\tend1\tchr2\tstart2\tend2";

    if(optimalrare==0){
      rareinter << "\traref";
    }

    rareinter<<"\tnRV";

    //FR
    if(FISHERflag==1){
      if(optimalrare==1){
	rareinter<<"\trarefreqFR";
	rareinter<<"\tnRareFR";
      }
      rareinter<<"\tp_FR\tCI_FR";
      if(wilsonpretest){
	rareinter<<"\tnPerm_FR";
      }
      if(verbose==2){rareinter<<"\tp_FR_Corr";}
    }
    //CMAT
    if(CMATflag==1){
      if(optimalrare==1){
	rareinter<<"\trarefreqCMAT";
	rareinter<<"\tnRareCMAT";
      }
      rareinter<<"\tp_CMAT\tCI_CMAT";
      if(wilsonpretest){
	rareinter<<"\tnPerm_CMAT";
      }
      if(verbose==2){rareinter<<"\tp_CMAT_Corr";}
      rareinter<<"\tOR_CMAT";
    }
    //COLL
    if(COLLflag==1){
      if(optimalrare==1){
	rareinter<<"\trarefreqCOLL";
	rareinter<<"\tnRareCOLL";
      }
      rareinter<<"\tp_COLL";
      if(nsim!=0){
	rareinter<<"\tCI_COLL";
	if(wilsonpretest && nsim!=0){
	  rareinter<<"\tnPerm_COLL";
	}
      }
      if(verbose==2){rareinter<<"\tp_COLL_Corr";}
      rareinter<<"\tOR_COLL";
    }

    // REG
    if(REGRESSIONflag==1){
      if(optimalrare==1){
	rareinter<<"\trarefreqREG";
	rareinter<<"\tnRareREG";
      }
      rareinter<<"\tp_REG";
      if(nsim!=0){
	rareinter<<"\tCI_REG";
	if(wilsonpretest && nsim!=0){
	  rareinter<<"\tnPerm_REG";
	}
      }
      if(verbose==2){rareinter<<"\tp_REG_Corr";}
    }

    // FRACREG
    if(FRACREGflag==1){
      if(optimalrare==1){
	rareinter<<"\trarefreqFRACREG";
	rareinter<<"\tnRareFRACREG";
      }
      rareinter<<"\tp_FRACREG";
      if(nsim!=0){
	rareinter<<"\tCI_FRACREG";
	if(wilsonpretest && nsim!=0){
	  rareinter<<"\tnPerm_FRACREG";
	}
      }
      if(verbose==2){rareinter<<"\tp_FRACREG_Corr";}
    }

    // FRACREG
    if(COLLREGflag==1){
      if(optimalrare==1){
	rareinter<<"\trarefreqCOLLREG";
	rareinter<<"\tnRareCOLLREG";
      }
      rareinter<<"\tp_COLLREG";
      if(nsim!=0){
	rareinter<<"\tCI_COLLREG";
	if(wilsonpretest && nsim!=0){
	  rareinter<<"\tnPerm_COLLREG";
	}
      }
      if(verbose==2){rareinter<<"\tp_COLLREG_Corr";}
    }

    if(intervalfile!=" " && featurecol>=0){
      rareinter<<"\tFEATURE1\tFEATURE2";
    }
    if(binamax!=0){
      rareinter<<"\tBINSIZE\tBINSTART";
    }

    rareinter<<"\n";

    start=0;
    end=0;
    k=0;

    // Bonferroni-corrected p-values
    float p_CMAT_Corr;
    float p_COLL_Corr;
    float p_FISHER_Corr;
    float p_REGRESSION_Corr;
    float p_FRACREG_Corr;
    float p_COLLREG_Corr;
    for(j=nwindowssinglebin; j<nwindows; j++){
      if(nRareSNPs[j]==0) continue;
      start1=windowPositions[doublewindowcoord[j][0]][0];
      end1=windowPositions[doublewindowcoord[j][0]][1];
      start2=windowPositions[doublewindowcoord[j][1]][0];
      end2=windowPositions[doublewindowcoord[j][1]][1];
      if(verbose==2){
	if(FISHERflag){
	  if(pFISHERvec[j]==0){
	    p_FISHER_Corr = (1/(float)nsim)*(float)(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	  else{
	    p_FISHER_Corr = (float)pFISHERvec[j]*(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	}
	if(CMATflag){
	  if(pCMATvec[j]==0){
	    p_CMAT_Corr = (1/(float)nsim)*(float)(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	  else{
	    p_CMAT_Corr = (float)pCMATvec[j]*(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	}
	if(COLLflag){
	  if(pCOLLvec[j]==0){
	    p_COLL_Corr = (1/(float)nsim)*(float)(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	  else{
	    p_COLL_Corr = (float)pCOLLvec[j]*(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	}
	if(REGRESSIONflag){
	  if(pREGRESSIONvec[j]==0){
	    p_REGRESSION_Corr = (1/(float)nsim)*(float)(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	  else{
	    p_REGRESSION_Corr = (float)pREGRESSIONvec[j]*(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	}
	if(FRACREGflag){
	  if(pFRACREGvec[j]==0){
	    p_FRACREG_Corr = (1/(float)nsim)*(float)(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	  else{
	    p_FRACREG_Corr = (float)pFRACREGvec[j]*(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	}
	if(COLLREGflag){
	  if(pCOLLREGvec[j]==0){
	    p_COLLREG_Corr = (1/(float)nsim)*(float)(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	  else{
	    p_COLLREG_Corr = (float)pCOLLREGvec[j]*(nwindowssinglebin*(nwindowssinglebin-1)/2);
	  }
	}
      }

      // Data output
      rareinter << map[start1].chr<< "\t"<< map[start1].pos << "\t" << map[end1].pos<<"\t" << map[start2].chr<< "\t"<< map[start2].pos << "\t" << map[end2].pos;


      if(optimalrare==0){
	rareinter << "\t" << window[j].maf_at_level[0];
      }

      rareinter << "\t" << nRareSNPs[j];
      if(optimalrare==1 && FISHERflag==1){
	rareinter << "\t" << rarefFISHER[j];
	rareinter << "\t"  <<window[j].n_at_level[MAF_Level_VT_FISHER[j]];
      }

      if(FISHERflag==1){
	rareinter <<"\t"<<pFISHERvec[j];
	if(FISHERnotconverged[j]==1){rareinter<<"**";}
	else if(wilsonpretest!=0 && FISHERpretestpassed[j]==1){rareinter<<"+";}
	rareinter<<"\t["<<FISHER_CI[j][0]<<","<<FISHER_CI[j][1]<<"]";
	if(wilsonpretest){
	  rareinter<<"\t"<<window[j].nperm[0];
	}
	if(verbose==2){rareinter<< "\t"<< p_FISHER_Corr;}
      }
      if(optimalrare==1 && CMATflag==1){
	rareinter << "\t" <<rarefCMAT[j]<<"\t";
	rareinter<<window[j].n_at_level[MAF_Level_VT_CMAT[j]];
      }
      if(CMATflag==1){
	rareinter<< "\t" << pCMATvec[j];
	if(wilsonpretest!=0 && CMATpretestpassed[j]==1){
	  rareinter<<"+";
	}
	rareinter<<"\t["<<CMAT_CI[j][0]<<","<<CMAT_CI[j][1]<<"]";
	if(wilsonpretest){
	  rareinter<<"\t"<<window[j].nperm[1];
	}
	if(verbose==2){ rareinter<< "\t" << p_CMAT_Corr;}
	rareinter << "\t" << OR_CMAT_vec[j];
      }
      if(optimalrare==1 && COLLflag==1){
	rareinter  << "\t" << rarefCOLL[j];
	rareinter<< "\t"<< window[j].n_at_level[MAF_Level_VT_COLL[j]];
      }
      if(COLLflag==1){
	rareinter <<"\t" << pCOLLvec[j];
	if(wilsonpretest!=0 && COLLpretestpassed[j]==1){
	  rareinter<<"+";
	}
	if(nsim!=0){
	  rareinter<<"\t["<<COLL_CI[j][0]<<","<<COLL_CI[j][1]<<"]";
	  if(wilsonpretest){
	    rareinter<<"\t"<<window[j].nperm[2];
	  }
	}
	if(verbose==2){rareinter << "\t" << p_COLL_Corr;}
	rareinter << "\t" << OR_COLL_vec[j];
	if(strcmp(map[start].chr, "23")==0){rareinter << "/" << OR_COLL_f_vec[j];}
      }
      // REGRESSION TESTS
      if(optimalrare==1 && REGRESSIONflag==1){
	rareinter << "\t" <<rarefREGRESSION[j]<<"\t";
	rareinter<<window[j].n_at_level[MAF_Level_VT_REGRESSION[j]];
      }
      if(REGRESSIONflag==1){
	rareinter<< "\t" << pREGRESSIONvec[j];
	if(wilsonpretest!=0 && REGRESSIONpretestpassed[j]==1){
	  rareinter<<"+";
	}
	if(nsim!=0){
	  rareinter<<"\t["<<REGRESSION_CI[j][0]<<","<<REGRESSION_CI[j][1]<<"]";
	  if(wilsonpretest){
	    rareinter<<"\t"<<window[j].nperm[3];
	  }
	}
	if(verbose==2){rareinter<< "\t" << p_REGRESSION_Corr;}
      }

      if(optimalrare==1 && FRACREGflag==1){
	rareinter << "\t" <<rarefFRACREG[j]<<"\t";
	rareinter<<window[j].n_at_level[MAF_Level_VT_FRACREG[j]];
      }
      if(FRACREGflag==1){
	rareinter<< "\t" << pFRACREGvec[j];
	if(wilsonpretest!=0 && FRACREGpretestpassed[j]==1){
	  rareinter<<"+";
	}
	if(nsim!=0){
	  rareinter<<"\t["<<FRACREG_CI[j][0]<<","<<FRACREG_CI[j][1]<<"]";
	  if(wilsonpretest){
	    rareinter<<"\t"<<window[j].nperm[4];
	  }
	}
	if(verbose==2){ rareinter<< "\t" << p_FRACREG_Corr;}
      }

      if(optimalrare==1 && COLLREGflag==1){
	rareinter << "\t" <<rarefCOLLREG[j]<<"\t";
	rareinter<<window[j].n_at_level[MAF_Level_VT_COLLREG[j]];
      }

      if(COLLREGflag==1){
	rareinter<< "\t" << pCOLLREGvec[j];
	if(wilsonpretest!=0 && COLLREGpretestpassed[j]==1){
	  rareinter<<"+";
	}
	if(nsim!=0){
	  rareinter<<"\t["<< COLLREG_CI[j][0]<<","<<COLLREG_CI[j][1]<<"]";
	  if(wilsonpretest){
	    rareinter<<"\t"<<window[j].nperm[5];
	  }
	}
	if(verbose==2){rareinter<< "\t" << p_COLLREG_Corr;}
      }
      if(intervalfile!=" " && featurecol>=0){
	rareinter<< "\t" << intervals[doublewindowcoord[j][0]].feature;
	rareinter<< "\t" << intervals[doublewindowcoord[j][1]].feature;
      }

      if(binamax!=0){
	rareinter<< "\t" << intervals[j].s1<< "\t" << intervals[j].s2;
      }

      rareinter << "\n";
    }

    rareinter.close();
    cout << "Rare variant interaction analysis is finished!"<<endl;
    logfile << "Rare variant interaction analysis is finished!"<<endl;
  }



void outRareVB(struct PERSON *person, struct MAP *map, struct COUNTS **counts, bool optimalrare, int nsim, string rarefileVB, string rarefileVBgraph, string rarefileVBpermstat, string outputname, fstream &errorfile, fstream &logfile, int rarepretest, float rarepretestlimit, double rareregpretest,  int *nSNPsInWindow, double raref, int thread, struct WINDOW *window, int minRareInBin,  int mafadjust, int nlinestped, int verbose, string SetIDfile, int intervaleditor, int setid, int *nvbstartvt, int **vbstartvt, int *nvbend, int **vbend, float *vbmaxpermstat, int **nvblevel, int ***nchunks, int ****chunkpos, int ****chunklen, int nwindows, int nwords, int* SNPMapInverse, uint64_t*** BinSNPs, uint64_t** BinSNPsCCFlagsOriginal, int ncasesqc, int ncontrolsqc, int ****ndummyatlevel, int *****dummypos, int *****dummylevel, uint64_t ***BinCarriers, int ****dummyend, int ***ndummyends, int ***nchunkcluster, int *****dummycluster, int vb_binwise_corr, int vb_print_perm_stat, float *****vbbinwisestat, short int *****vbbinwisecount){

  cout << "\nComputing and writing results of rare variant VB analysis..." << endl;

  fstream rareVB;
  rareVB.open(rarefileVB.c_str(), ios::out);


  fstream rareVBpermstat;
  rareVBpermstat.open(rarefileVBpermstat.c_str(), ios::out);
  for(int i=0; i<nsim; i++){
    rareVBpermstat<<pValueCalc(0.5, vbmaxpermstat[i]/2)<<endl;
  }
  rareVBpermstat.close();


  float minvbpermstat=0;
  if(nsim>0){
    minvbpermstat=vbmaxpermstat[0];
  }

  if(vb_print_perm_stat==1){
    for(int i=0; i<nsim; i++){
      rareVBpermstat<<vbmaxpermstat[i]<<endl;
      if((minvbpermstat-vbmaxpermstat[i])>EPS){
	minvbpermstat=vbmaxpermstat[i];
      }
    }
    rareVBpermstat.close();
  }
  if(vb_print_perm_stat==0){
    qsort(vbmaxpermstat, nsim, sizeof(float), comparefloat);
  }

  fstream rareVBgraph;
  int vbgraph=0;
  if(optimalrare==0){
    vbgraph=0;
  }
  // output graph data
  if(vbgraph==1){
    rareVBgraph.open(rarefileVBgraph.c_str(), ios::out);
        for(int i=0; i<nwindows; i++){
	  //	  for(int j=0; j<nvbstartvt[i]; j++){
	  int j=0;
	  int *carriersatlevel=new int[nvblevel[i][j]]();

	  for(int k=0; k<nvblevel[i][j]; k++){
	    carriersatlevel[k]=nCarriers[vbstartvt[i][chunkpos[i][j][k][0]]];
	    int cluster=0;
	    int *founddummypos=new int [nchunkcluster[i][j][k]]();
	    for(int l=0; l<nchunks[i][j][k]; l++){
	      if(nchunkcluster[i][j][k]>1){
		if(cluster<nchunkcluster[i][j][k]-1 && chunkpos[i][j][k][l]>dummyend[i][j][k][cluster]){
		  cluster++;
		}
	      }
	      for(int m=0; m<chunklen[i][j][k][l]; m++){
		int pos=chunkpos[i][j][k][l]+m;
		while(ndummyatlevel[i][j][k][cluster]>0  &&  pos==dummypos[i][j][k][cluster][founddummypos[cluster]] && founddummypos[cluster]<ndummyatlevel[i][j][k][cluster]){
		  int foundpos=0;
		  int pos2=0;
		  for(int l2=0; l2<nchunks[i][j][dummylevel[i][j][k][cluster][founddummypos[cluster]]]; l2++){
		    for(int m2=0; m2<chunklen[i][j][dummylevel[i][j][k][cluster][founddummypos[cluster]]][l2]; m2++){
		      pos2=chunkpos[i][j][dummylevel[i][j][k][cluster][founddummypos[cluster]]][l2]+m2;
		      if(pos2>pos){
			foundpos=1;
			break;
		      }
		    }
		    if(foundpos==1){
		      break;
		    }
		  }
		  rareVBgraph<<pos+1<<"\t"<< pos2+1<<"\t"<<nCarriers[vbstartvt[i][chunkpos[i][j][dummylevel[i][j][k][cluster][founddummypos[cluster]]][0]]]<<"\n";
		  founddummypos[cluster]++;
		}

		if(pos==chunkpos[i][j][k][nchunks[i][j][k]-1]+chunklen[i][j][k][nchunks[i][j][k]-1]-1){ // last pos
		  break;
		}
		else if(m<chunklen[i][j][k][l]-1){ // next bin in same chunk
		  rareVBgraph<<pos+1<<"\t"<<pos+2<<"\t"<<carriersatlevel[k]<<"\n";
		}
		else if(m==chunklen[i][j][k][l]-1){ // next bin in same chunk
		  rareVBgraph<<pos+1<<"\t"<<chunkpos[i][j][k][l+1]+1<<"\t"<<carriersatlevel[k]<<"\n";
		}
	      }
	    }
	  }

	}
	//	  } startpos
	rareVBgraph.close();
  }

  rareVB << "#INFO:";
  rareVB<<"SIMULATION="<<nsim<<";";

  rareVB << "MAFT="<<raref<<";";
  rareVB << "VT="<<optimalrare<<";";
  //    rareVB << "MAF_ADJUST="<<mafadjust<<";";
  if(minRareInBin>0){
    rareVB<<"MIN_RARE_IN_BIN="<<minRareInBin<<";";
  }

#if PARALLELN
  int maxthreads;
  maxthreads=omp_get_max_threads();
  if(maxthreads>MAXTHREAD){
    maxthreads=MAXTHREAD;
  }
  rareVB<<"MAXTHREADS="<<maxthreads<<";";
#endif
#if !PARALLELN
  rareVB<<"MAXTHREADS=1;";
#endif

  rareVB<<"\n";
  rareVB<<"chr\tstart\tend";
  if(verbose==2){
    rareVB<<"\tdistinct_rare_start\tdistinct_rare_end";
  }
  rareVB<<"\tnRV";
  rareVB<<"\tCarrierLevel";

  if(vb_print_perm_stat==0){
    rareVB<<"\tp_asymp";
  }
  else if(vb_print_perm_stat==1){
      rareVB<<"\tteststat_asymp";
  }

  if(nsim!=0 && vb_binwise_corr==true){
    rareVB<<"\tp_bw_corr";
  }
  if(nsim!=0 && vb_print_perm_stat==0){
    rareVB<<"\tp_gw_corr";
  }
  rareVB<<"\n";

  //  cout<<" minvbpermstat "<<minvbpermstat<<endl;
  if(vb_binwise_corr==false){
    int *carriersatlevel_order=NULL;
    double vbstat=0;
    int sX, sY;
    for(int i=0; i<nwindows; i++){
      for(int j=0; j<nvbstartvt[i]; j++){
	int carrierlevel=0;
	int *nRareVariants=new int [nvblevel[i][j]]();
	if(!nRareVariants) die("Memory allocation error in nRareVariants!");
	uint64_t*** dummy = new uint64_t**[nvblevel[i][j]];
	if(!dummy) die("Memory allocation error in dummy!");
	for(int k=0; k<nvblevel[i][j]; k++){
	  dummy[k]=new uint64_t*[nchunkcluster[i][j][k]];
	  if(!dummy[k]) die("Memory allocation error in dummy[k]!");
	  for(int k1=0; k1<nchunkcluster[i][j][k]; k1++){
	    dummy[k][k1]=new uint64_t[nwords]();
	    if(!dummy[k][k1]) die("Memory allocation error in dummy[k][k1]!");
	  }
	}
	if(!optimalrare){
	  carrierlevel=nCarriers[vbstartvt[i][chunkpos[i][j][0][0]]];
	}
	else if(optimalrare){
	  int *carriersatlevel_unsort=new int[nvblevel[i][j]]();
	  int *carriersatlevel=new int[nvblevel[i][j]]();
	  carriersatlevel_order=new int[nvblevel[i][j]]();
	  for(int k=0; k<nvblevel[i][j]; k++){
	    carriersatlevel_unsort[k]=nCarriers[vbstartvt[i][chunkpos[i][j][k][0]]];
	    carriersatlevel[k]=nCarriers[vbstartvt[i][chunkpos[i][j][k][0]]];
	  }
	  qsort(carriersatlevel, nvblevel[i][j], sizeof(int), compareint);
	  for(int k=1; k<nvblevel[i][j]; k++){
	    for(int kk=1; kk<nvblevel[i][j]; kk++){
	      if(carriersatlevel[k]==carriersatlevel_unsort[kk]){
		carriersatlevel_order[k]=kk;
		//		cout<<i<<" "<<j<<" "<<kk<<endl;
		break;
	      }
	    }
	    //	    cout<<i<<" "<<j<<" "<<k<<" "<<carriersatlevel_order[k]<<endl;
	  }
	  delete[] carriersatlevel_unsort;
	  delete[] carriersatlevel; // to delete: carriers.levelorder
	}

	for(int kk=0; kk<nvblevel[i][j]; kk++){
	  int k=kk;
	  if(optimalrare){
	    k=carriersatlevel_order[kk];
	  }

	  int *founddummypos=new int [nchunkcluster[i][j][k]]();
	  if(!founddummypos) die("Memory allocation error in founddummypos!");
	  int cluster=0;
	  for(int l=0; l<nchunks[i][j][k]; l++){
	    if(nchunkcluster[i][j][k]>1){
	      if(cluster<nchunkcluster[i][j][k]-1 && chunkpos[i][j][k][l]>dummyend[i][j][k][cluster]){
		cluster++;
	      }
	    }
	    for(int m=0; m<chunklen[i][j][k][l]; m++){
	      if(!optimalrare){
		if(carrierlevel<nCarriers[vbstartvt[i][chunkpos[i][j][0][l]+m]]){
		  carrierlevel=nCarriers[vbstartvt[i][chunkpos[i][j][0][l]+m]];
		}
	      }

	      int vbprop=nsim;
	      nRareVariants[k]++;
	      sY=0;
	      sX=0;
	      int pos=chunkpos[i][j][k][l]+m;
	      int i2=vbstartvt[i][pos];
	      int iMod = SNPMapInverse[i2];
	      if (iMod==-1) continue;

	      for (int p=0; p<nwords; p++){
		dummy[k][cluster][p] |= BinCarriers[i][pos][p];
	      }
	      while(ndummyatlevel[i][j][k][cluster]>0  &&  pos==dummypos[i][j][k][cluster][founddummypos[cluster]] && founddummypos[cluster]<ndummyatlevel[i][j][k][cluster]){
		for (int p=0; p<nwords; p++){
 		  dummy[dummylevel[i][j][k][cluster][founddummypos[cluster]]][dummycluster[i][j][k][cluster][founddummypos[cluster]]][p]=dummy[k][cluster][p];
		  nRareVariants[dummylevel[i][j][k][cluster][founddummypos[cluster]]]=nRareVariants[k];
		}
		founddummypos[cluster]++;
	      }
	      for(int p=0; p<nwords; p++){
		sY += bitcount64(dummy[k][cluster][p] & BinSNPsCCFlagsOriginal[p][1]);
		sX += bitcount64(dummy[k][cluster][p] & BinSNPsCCFlagsOriginal[p][2]);
	      }

	      vbstat=((double)ncasesqc + (double)ncontrolsqc)*((double)sX*(double)ncontrolsqc - (double)sY*(double)ncasesqc)*((double)sX*(double)ncontrolsqc-(double)sY*(double)ncasesqc)/((double)ncasesqc*(double)ncontrolsqc*((double)sX+(double)sY)*((double)ncasesqc+(double)ncontrolsqc-(double)sX-(double)sY));
	      if((ncasesqc+ncontrolsqc-sY-sX)==0 || (sX == 0 && sY==0)) vbstat=0;
	      //    cout <<"vbstat "<<vbstat<<endl;
	      if(verbose<2 && vb_print_perm_stat==0 && (vbstat-minvbpermstat)>EPS){
		for(int i2=0; i2<nsim; i2++){
		  if((vbstat-vbmaxpermstat[i2])>EPS){
		    vbprop--;
		  }
		  else{
		    break;
		  }
		}
		//      cout<<map[vbstartvt[i][0]].chr<<"\t"<<map[vbstartvt[l][chunkpos[i][j][k][0]]].pos<<"\t"<<map[pos].pos<<"\t"<<nRareVariants<<"\t"<< pValueCalc(0.5, float(vbstat)/float(2))<<"\t"<<float(vbprop)/float(nsim)<<endl;
		rareVB<<map[vbstartvt[i][0]].chr<<"\t"<<map[vbstartvt[i][chunkpos[i][j][0][0]]].pos<<"\t"<<map[pos].pos;
		if(verbose==2){
		  rareVB<<"\t"<< chunkpos[i][j][k][l]<<"\t"<<chunkpos[i][j][k][l]+m;
		}
		rareVB<<"\t"<<nRareVariants[k];
		if(optimalrare){
		  rareVB<<"\t"<<nCarriers[vbstartvt[i][chunkpos[i][j][k][0]]];
		}
		else if(!optimalrare){
		  rareVB<<"\t"<<carrierlevel;
		}
		if(vb_print_perm_stat==0){
		  rareVB<<"\t"<< pValueCalc(0.5, float(vbstat)/float(2));
		}
		else if(vb_print_perm_stat==1){
		  rareVB<<"\t"<< vbstat;
		}
		if(nsim!=0 && vb_print_perm_stat==0){
		  rareVB<<"\t"<<float(vbprop)/float(nsim);
		}
		rareVB<<endl;

	      }
	      else if(verbose==2){
		if(nsim!=0 && vb_print_perm_stat==0){
		  for(int i2=0; i2<nsim; i2++){
		    if((vbstat-vbmaxpermstat[i2])>EPS){
		      vbprop--;
		    }
		    else{
		      break;
		    }
		  }
		}
		//cout<<nvblevel[i][j]<<endl;
		//cout<<"clulster "<<cluster<<endl;
		//cout<<" pos "<<pos<<endl;

		rareVB<<map[vbstartvt[i][0]].chr<<"\t"<<map[vbstartvt[i][chunkpos[i][j][0][0]]].pos<<"\t"<<map[vbstartvt[i][pos]].pos;
		rareVB<<"\t"<< chunkpos[i][j][k][0]+1<<"\t"<<chunkpos[i][j][k][l]+m+1;
		rareVB<<"\t"<<nRareVariants[k];

		if(optimalrare){
		  rareVB<<"\t"<<nCarriers[vbstartvt[i][chunkpos[i][j][k][0]]];
		}
		else if(!optimalrare){
		  rareVB<<"\t"<<carrierlevel;
		}
		if(vb_print_perm_stat==0){
		  rareVB<<"\t"<< pValueCalc(0.5, float(vbstat)/float(2));
		}
		else if(vb_print_perm_stat==1){
		  rareVB<<"\t"<< vbstat;
		}
		if(nsim!=0 && vb_print_perm_stat==0){
		  rareVB<<"\t"<<float(vbprop)/float(nsim);
		}
		rareVB<<flush<<endl;
	      }
	    }
	  }
	  //  cout<<ndummyatlevel[i][j][k][cluster]<<" "<<founddummypos[cluster]<<endl;
	  if(ndummyatlevel[i][j][k][cluster]!=founddummypos[cluster]){
	    for(int l=0; l<ndummyatlevel[i][j][k][cluster]; l++){
	      	      cout<<"dummypos["<<i<<"]["<<j<<"]["<<k<<"]["<<cluster<<"]["<<l<<"] "<<dummypos[i][j][k][cluster][l]<<endl;
	      	      cout<<"dummylevel["<<i<<"]["<<j<<"]["<<k<<"]["<<cluster<<"]["<<l<<"] "<<dummylevel[i][j][k][cluster][l]<<endl;
	    }
	    cout<<" ndummyatlevel["<<i<<"]["<<j<<"]["<<k<<"]["<<cluster<<"] "<<ndummyatlevel[i][j][k][cluster]<<" "<<founddummypos[cluster]<<endl;
	    cout<<"BUG in ndummyatlevel!"<<endl;
	    //    exit(1);
	  }
	  delete[] founddummypos;

	}
	if(optimalrare){
	  delete[] carriersatlevel_order;
	}

	for(int k1=0; k1<nvblevel[i][j]; k1++){

	  for(int k2=0; k2<nchunkcluster[i][j][k1]; k2++){
	    delete[] dummy[k1][k2];
	  }
	  delete[] dummy[k1];
	}
	delete[] dummy;
	delete[] nRareVariants;
      }
    }
  }
  else if(vb_binwise_corr==true){ // only for FT
    double vbstat=0;
    int cluster=0;
    for(int i=0; i<nwindows; i++){
      for(int j=0; j<nvbstartvt[i]; j++){
	int carrierlevel=0;
	int *nRareVariants=new int [nvblevel[i][j]]();
	if(!nRareVariants) die("Memory allocation error in nRareVariants!");
	uint64_t*** dummy = new uint64_t**[nvblevel[i][j]];
	if(!dummy) die("Memory allocation error in dummy!");
	for(int k1=0; k1<nvblevel[i][j]; k1++){ // always 1 with vb_binwise_corr because then vt=0
	  // int distinctbins=0;
	  dummy[k1]=new uint64_t*[nchunkcluster[i][j][k1]];
	  if(!dummy[k1]) die("Memory allocation error in dummy[k1]!");
	  for(int k2=0; k2<nchunkcluster[i][j][k1]; k2++){
	    dummy[k1][k2]=new uint64_t[nwords]();
	    if(!dummy[k1][k2]) die("Memory allocation error in dummy[k1][k2]!");
	  }
	  if(!optimalrare){ // always with vb_binwise_corr
	    carrierlevel=nCarriers[vbstartvt[i][chunkpos[i][j][0][0]]];
	  }
	  for(int k=0; k<nvblevel[i][j]; k++){
	    for(int l=0; l<nchunks[i][j][k]; l++){
	      for(int m=0; m<chunklen[i][j][k][l]; m++){
		if(!optimalrare){ // always with vb_binwise_corr
		  if(carrierlevel<nCarriers[vbstartvt[i][chunkpos[i][j][0][l]+m]]){
		    carrierlevel=nCarriers[vbstartvt[i][chunkpos[i][j][0][l]+m]];
		  }
		}

		int vbprop=nsim;
		nRareVariants[k]++;
		int pos=vbstartvt[i][chunkpos[i][j][k][l]+m];

		if(verbose<2 && (vbbinwisestat[i][j][k][l][m]-minvbpermstat)>EPS){ // never because vebose=2 with VB_BINWISE_CORR
		  for(int i2=0; i2<nsim; i2++){
		    if((vbbinwisestat[i][j][k][l][m]-vbmaxpermstat[i2])>EPS){
		      vbprop--;
		    }
		    else{
		      break;
		    }
		  }
		  //      cout<<map[vbstartvt[i][0]].chr<<"\t"<<map[vbstartvt[l][chunkpos[i][j][k][0]]].pos<<"\t"<<map[pos].pos<<"\t"<<nRareVariants<<"\t"<< pValueCalc(0.5, float(vbstat)/float(2))<<"\t"<<float(vbprop)/float(nsim)<<endl;
		  rareVB<<map[vbstartvt[i][0]].chr<<"\t"<<map[vbstartvt[i][chunkpos[i][j][0][0]]].pos<<"\t"<<map[pos].pos;

		  rareVB<<"\t"<<nRareVariants[k];

		  if(!optimalrare){ // always with vb_binwise_corr
		    rareVB<<"\t"<<carrierlevel;
		  }

		  rareVB<<"\t"<< pValueCalc(0.5, float(vbbinwisestat[i][j][k][l][m])/float(2));

		  if(nsim!=0 && vb_binwise_corr==true){ // always with vb_binwise_corr
		    rareVB<<"\t"<<float(vbbinwisecount[i][j][k][l][m])/nsim;
		  }
		  if(nsim!=0){
		    rareVB<<"\t"<<float(vbprop)/float(nsim);
		  }
		  rareVB<<endl;

		}
		else if(verbose==2){
		  if(nsim!=0){
		    for(int i2=0; i2<nsim; i2++){
		      if((vbbinwisestat[i][j][k][l][m]-vbmaxpermstat[i2])>EPS){
			vbprop--;
		      }
		      else{
			break;
		      }
		    }
		  }
		  rareVB<<map[vbstartvt[i][0]].chr<<"\t"<<map[vbstartvt[i][chunkpos[i][j][0][0]]].pos<<"\t"<<map[pos].pos;
		  if(verbose==2){
		    rareVB<<"\t"<< chunkpos[i][j][k][0]+1<<"\t"<<chunkpos[i][j][k][l]+m+1;
		  }
		  rareVB<<"\t"<<nRareVariants[k];

		  if(optimalrare){ // never
		    rareVB<<"\t"<<nCarriers[vbstartvt[i][chunkpos[i][j][k][0]]];
		  }
		  else if(!optimalrare){
		    rareVB<<"\t"<<carrierlevel;
		  }

		  rareVB<<"\t"<< pValueCalc(0.5, float(vbbinwisestat[i][j][k][l][m])/float(2));

		  if(nsim!=0 && vb_binwise_corr==true){
		    rareVB<<"\t"<<float(vbbinwisecount[i][j][k][l][m])/nsim;
		  }
		  if(nsim!=0){
		    rareVB<<"\t"<<float(vbprop)/float(nsim);
		  }
		  rareVB<<endl;
		}
	      }
	    }
	  }

	  for(int k1=0; k1<nvblevel[i][j]; k1++){

	    for(int k2=0; k2<nchunkcluster[i][j][k1]; k2++){
	      delete[] dummy[k1][k2];
	    }
	    delete[] dummy[k1];
	  }
	  delete[] dummy;
	  delete[] nRareVariants;
	}
      }
    }
  }
  rareVB.close();

  cout << "Rare variant VB analysis is finished!"<<endl;
  logfile << "Rare variant VB analysis is finished!"<<endl;
}


void outRareVB_BMP(struct PERSON *person, struct MAP *map, struct COUNTS **counts, bool optimalrare, int nsim, string rarefileVB_BMP, string outputname, fstream &errorfile, fstream &logfile, int rarepretest, float rarepretestlimit, double rareregpretest,  int *nSNPsInWindow, double raref, int thread, struct WINDOW *window, int minRareInBin,  int mafadjust, int nlinestped, int verbose, string SetIDfile, int intervaleditor, int setid, int nvbstart, int *vbstart, int *nvbend, int **vbend, float *vbmaxpermstat, int ncasesqc, int ncontrolsqc, uint64_t*** BinSNPs, uint64_t** BinSNPsCCFlags, int nwords, int* SNPMapInverse)
{

  cout << "\nBMP output of rare variant VB analysis is being generated..." << endl;


  unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
  unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
  unsigned char bmppad[3] = {0,0,0};

  int w,h;
  w=nlinestped;
  h=nlinestped;


  int filesize = 54 + 3*nlinestped*nlinestped;
  bmpfileheader[ 2] = (unsigned char)(filesize    );
  bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
  bmpfileheader[ 4] = (unsigned char)(filesize>>16);
  bmpfileheader[ 5] = (unsigned char)(filesize>>24);
  bmpinfoheader[ 4] = (unsigned char)(       w    );
  bmpinfoheader[ 5] = (unsigned char)(       w>> 8);
  bmpinfoheader[ 6] = (unsigned char)(       w>>16);
  bmpinfoheader[ 7] = (unsigned char)(       w>>24);
  bmpinfoheader[ 8] = (unsigned char)(       h    );
  bmpinfoheader[ 9] = (unsigned char)(       h>> 8);
  bmpinfoheader[10] = (unsigned char)(       h>>16);
  bmpinfoheader[11] = (unsigned char)(       h>>24);

  //    fstream rareVB;


  FILE *rareVB;
  rareVB=fopen(rarefileVB_BMP.c_str(), "wb");

  fwrite(bmpfileheader,1,14,rareVB);
  fwrite(bmpinfoheader,1,40,rareVB);


  int imgstat=0;
  if(imgstat==1){


    unsigned char *img=NULL;
    img = (unsigned char *)malloc(3*w*h);
    if(!img){
      die("Memory allocation error in img!\n");
    }
    memset(img,255,3*w*h*sizeof(unsigned char));



    float **stats=NULL;
    stats = new float *[w];
    if(!stats){
      die("Memory allocation error in stats!\n");
    }
    for(int i=0; i<h; i++){
      stats[i] = new float [w]();
      if(!stats[i]){
	die("Memory allocation error in stats[i]!\n");
      }

    }

    uint64_t* dummy1 = new uint64_t[nwords];
    if(!dummy1) die("Memory allocation error in dummy1!");
    float statmax_pos=0;
    float statmax_neg=0;

    for(int j1=0; j1<nlinestped; j1++){
      dummy1=(uint64_t *)calloc(nwords, sizeof(uint64_t));
      if(!dummy1) die("Memory allocation error in dummy1!");
      for(int i1=j1; i1<nlinestped; i1++){
	int sX=0,sY=0;
	double stat=0;

	int iMod = SNPMapInverse[i1];
	if (iMod==-1) continue;
	for (int p=0; p<nwords; p++){
	  dummy1[p] |= (BinSNPs[iMod][p][1] | BinSNPs[iMod][p][2]);
	}
	for(int p=0; p<nwords; p++){
	  sY += bitcount64(dummy1[p] & BinSNPsCCFlags[p][1]);
	  sX += bitcount64(dummy1[p] & BinSNPsCCFlags[p][2]);
	}
	stat=((double)ncasesqc + (double)ncontrolsqc)*((double)sX*(double)ncontrolsqc - (double)sY*(double)ncasesqc)*((double)sX*(double)ncontrolsqc-(double)sY*(double)ncasesqc)/((double)ncasesqc*(double)ncontrolsqc*((double)sX+(double)sY)*((double)ncasesqc+(double)ncontrolsqc-(double)sX-(double)sY));
	if((ncasesqc+ncontrolsqc-sY-sX)==0 || (sX == 0 && sY==0)) stat=0;


	float propNUC=float(sY)/float(ncontrolsqc);
	float propNAC=float(sX)/float(ncasesqc);
	if((propNUC-propNAC)>EPS){ // protective
	  stats[j1][i1]=-(float)stat;
	  if((stat-statmax_neg)>EPS){
	    statmax_neg=stat;
	  }
	}
	else if((propNAC-propNUC)>EPS){ // damaging
	  stats[j1][i1]=(float)stat;
	  if((stat-statmax_pos)>EPS){
	    statmax_pos=stat;
	  }
	}
	//      cout<< stats[j1][i1]<<endl;
      }
      delete[] dummy1;
    }
    cout<<"statmax_neg "<<statmax_neg<<endl;
    cout<<"statmax_pos "<<statmax_pos<<endl;

    for(int i=0; i<w; i++)
      {
	for(int j=i; j<h; j++)
	  {
	    int y=i;
	    int x=j;
	    float r=0;
	    float g=0;
	    float b=0;
	    if(stats[i][j]>0){
	      r=255;
	      g=255-(stats[i][j]/statmax_pos)*255;
	      b=255-(stats[i][j]/statmax_pos)*255;
	    }
	    else if(stats[i][j]<0){
	      b=255;
	      r=255-(-stats[i][j]/statmax_neg)*255;
	      g=255-(-stats[i][j]/statmax_neg)*255;
	    }

	    if ((r-255)>EPS) r=255;
	    if ((g-255)>EPS) g=255;
	    if ((b-255)>EPS) b=255;

	    img[(x+y*w)*3+2] = (unsigned char)(r);
	    img[(x+y*w)*3+1] = (unsigned char)(g);
	    img[(x+y*w)*3+0] = (unsigned char)(b);

	    /*

	                  if (r<EPS) r=0;
			            if (g<EPS) g=0;
				            if (b<EPS) b=0;

					          if(fabs(stats[i][j])<EPS){
						              break;
							      } */
	  }
      }

    for(int i=0; i<h; i++)
      {
	fwrite(img+(w*(h-i-1)*3),3,w,rareVB);
	fwrite(bmppad,1,(4-(w*3)%4)%4,rareVB);
      }

    free(stats);

  }

  else if(imgstat==0){

    for(int j1=nlinestped-1; j1>=0; j1--){
      unsigned char *img=NULL;
      img = (unsigned char *)calloc(3*w, sizeof(int));
      if(!img) die("Memory allocation error in img!");
      //      memset(img,255,sizeof(img));

      uint64_t* dummy1 = new uint64_t[nwords]();
      if(!dummy1) die("Memory allocation error in dummy1!");

      for(int i1=j1; i1<nlinestped; i1++){
	int sX=0,sY=0;
	double stat=0;


	int y=j1;
	int x=i1;

	int iMod = SNPMapInverse[i1];
	if (iMod==-1) continue;
	for (int p=0; p<nwords; p++){
	  dummy1[p] |= (BinSNPs[iMod][p][1] | BinSNPs[iMod][p][2]);
	}

	int color1=0;
	int color2=0;
	int color3=0;
	uint64_t* dummycolor = new uint64_t[3]();
	if(!dummycolor) die("Memory allocation error in dummycolor!");
	for(int l=0; l<nwords/3; l++){
	  dummycolor[0]^=dummy1[l];
	}
	dummycolor[0]=dummycolor[0]%255;

	for(int m=nwords/3; m<(2*nwords/3); m++){
	  dummycolor[1]^=dummy1[m];
	}
	dummycolor[1]=dummycolor[1]%255;

	for(int n=2*nwords/3; n<nwords; n++){
	  dummycolor[2]^=dummy1[n];
	}
	dummycolor[2]=dummycolor[2]%255;

	img[(x)*3+2] = (unsigned char)(dummycolor[0]);
	img[(x)*3+1] = (unsigned char)(dummycolor[1]);
	img[(x)*3+0] = (unsigned char)(dummycolor[2]);

	delete[] dummycolor;

      }
      fwrite(img,3,w,rareVB);
      fwrite(bmppad,1,(4-(w*3)%4)%4,rareVB);

      delete[] dummy1;


      free(img);

    }


    /*    for(int i=0; i<w; i++)
	      {
	            for(int j=i; j<h; j++)
		            {
			      int y=i;
			          int x=j;
				        float r=0;
					        float g=0;
						  float b=0;
						      if(stats[i][j]>0){
						            r=255;
							            g=255-(stats[i][j]/statmax_pos)*255;
								      b=255-(stats[i][j]/statmax_pos)*255;
								          }
									        else if(stats[i][j]<0){
										        b=255;
											  r=255-(-stats[i][j]/statmax_neg)*255;
											      g=255-(-stats[i][j]/statmax_neg)*255;
											            }

												            if ((r-255)>EPS) r=255;
													      if ((g-255)>EPS) g=255;
													          if ((b-255)>EPS) b=255;

														        img[(x+y*w)*3+2] = (unsigned char)(r);
															        img[(x+y*w)*3+1] = (unsigned char)(g);
																  img[(x+y*w)*3+0] = (unsigned char)(b);


																    //      if (r<EPS) r=0;
																      //      if (g<EPS) g=0;
																        //      if (b<EPS) b=0;

																	  //      if(fabs(stats[i][j])<EPS){
																	    //      break;
																	      //      }
																	          }
																		  } */



  }

  fclose(rareVB);
  cout << "Rare variant VB BMP output is finished!"<<endl;
  logfile << "Rare variant VB BMP output is finished!"<<endl;
  exit(0);
}

/*
  void outRareVB_BMP_level(struct PERSON *person, struct MAP *map, struct COUNTS **counts, bool optimalrare, int nsim, string rarefileVB_BMP, string outputname, fstream &errorfile, fstream &logfile, int rarepretest, float rarepretestlimit, double rareregpretest,  int *nSNPsInWindow, double raref, int thread, struct WINDOW *window, int minRareInBin,  int mafadjust, int nlinestped, int verbose, string SetIDfile, int intervaleditor, int setid, int nvbstart, int *vbstart, int *nvbend, int **vbend, float **vbstat, float *vbmaxpermstat, int ncases, int ncontrols, uint64_t*** BinSNPs, uint64_t** BinSNPsCCFlags, int nwords, int* SNPMapInverse, int *nvbstartvt, int **nvblevel, int ***nchunksatlevel, int ****vbcoords){

  }

*/



