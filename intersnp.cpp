/*
 *
 *      PROGRAM: INTERSNP
 *      software for genome-wide interaction analysis (GWIA) of case-control SNP data
 *      SNPs are selected for joint analysis using a priori information
 *      Sources of information to define meaningful strategies can be statistical evidence
 *      (single marker association at a moderate level, computed from the own data)
 *      and genetic/biologic relevance (genomic location, function class or pathway information).
 *
 *      AUTHOR:         Christine Herold, Dmitriy Drichel, Andre Lacour, Tatsiana Vaitsiakhovich, Tim Becker
 *
 *	                	Institute for Medical Biometry,
 *      				Informatics, and Epidemiology (IMBIE),
 *             			University of Bonn
 *      				Sigmund-Freud-Strasse 25
 *             			53127 Bonn
 *      				Germany
 *
 *	                	German Center for Neurodegenerative Diseases, DZNE
 *						Sigmund-Freud-Strasse 25
 *             			53127 Bonn
 *      				Germany
 *
 *
 *
 *
 *      VERSION: 1.15
 *
 *      YEAR: February 2015 (11/02/2015)
 *
 *      DOCUMENTATION: http://intersnp.meb.uni-bonn.de
 *
 *
 *
 */

#define _GLIBCXX_USE_CXX11_ABI 0 

#define RARE 1
#define DEV 0
#define devREGRESSION 0

#define PARALLELN 0
#define PARALLELA 0

#if PARALLELA
 #include <omp.h>
 #define MAXTHREAD 23
 #define PARALLEL 1
#elif PARALLELN
 #include <omp.h>
 #define MAXTHREAD 23
 #define PARALLEL 1
#else
 #define MAXTHREAD 1
#endif

/// std headers
#include <cctype>
#include <climits>
#include <cmath>
#include <stdint.h>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
/// stream libraries
#include <iostream>  // wegen cout cin
#include <fstream> // wegen Dateistreamobjekt
#include <sstream>
#include <iomanip>
/// classes
#include <string> // wegen Datentyp string
#include <vector>
#include <algorithm>
#include <set> // wegen Container set



#define ITMAX 200
#define FPMIN 1.0e-30
#define EPS 1.19209290e-07F
#define MOD 0


using namespace std;


double genotypTest(struct COUNTS counts, int teststat, double *inflationfactor, double *fulltests);
double snpTestX(struct COUNTS counts, int teststat,  double *inflationfactor, double *fulltests);
double armitageTest(struct COUNTS counts, int teststat, double *inflationfactor, double *fulltests);
	int singleMarkerTest = 1;


struct TRAITAVG
{
 float snp1allele[2];
 float snp2allele[2];
 float allele[2][2];
 float snp1geno[3];
 float snp2geno[3];
 float geno[3][3];
};

struct PERSON
{
	char *fid; // Family-ID
	char *pid; // Individual-ID
	char *vid; // Paternal-ID
	char *mid; // Maternal-ID
	unsigned char sex; // Geschlecht
	unsigned char aff[MAXTHREAD]; // Affection-Status. also in thread
	int missing;
	float missingrate; // Missingrate pro Person (Spalte)
	bool qcin;
	bool analysis_in;
	double qtaff[MAXTHREAD];
	double *cov;
	unsigned char *covin;
	short unsigned int allcovin;
	float load;
	int clusterAffil;
};

struct FAMILY
{
	char *fid; // Family-ID
	short unsigned int fatherIndex;
	short unsigned int motherIndex;
	int *childIndex;
	int *personList;
	int nChildren;
	int nFamMember;
};

struct MAP
{
	char chr[3];        // Chromosome
	char *rs;           // ID
/**/double rsNum;       // ggf. Probleme mit accuracy und Länge der rs-Nummern, da [[alpha]] auch durch Ziffern ersetzt wird -> ersetze durch long
	unsigned int dist;  // genetic distance
	int pos;            // physical position on chromosome
	unsigned int line;  // line in input file
    string minorA;
    string majorA;
	bool analysis_in;   // preselection for analysis
	bool matching_in;   // preselection for matching
    bool done;
	bool qcin;          // QC
	bool in;
	int missing;        // Missings insgesamt
	float missingrate;

	double maf;          // MAF for qc (case+control)
	double controlmaf;   // MAF for showing results (controls)
	double casemaf;      // MAF for showing results (cases)
#if RARE
	double mafa;          // MAF for qc (case+control) MAF_ADJUST
	double controlmafa;   // MAF for showing results (controls) MAF_ADJUST
	double casemafa;      // MAF for showing results (cases) MAF_ADJUST

	double mafr;          // MAF for qc (case+control) RARE
	double controlmafr;   // MAF for showing results (controls) RARE
	double casemafr;      // MAF for showing results (cases) RARE
#endif
	char *gene;
	int location;
	int locationToGene;
	int codingStatus;
	double p;
	double pmod;
};

struct MAP2
{
	char *rs; // rs-Nummer
	int line; // Zeile im tped-File
};

// Family data
struct MAP3
{
	char *fid; // Family-ID
	char *pid; // Individual-ID
	char *vid; // Paternal-ID
	char *mid; // Maternal-ID
	unsigned char sex; // Geschlecht
	int line; // Zeile im tfam-File
};

struct INFOMAP
{
	char chr[3];
	char *rs;       // ID
/**/double rsNum;   // ggf. Probleme mit accuracy und Länge der rs-Nummern, da [[alpha]] auch durch Ziffern ersetzt wird -> ersetze durch long
	unsigned int pos;
	unsigned int line;
	char *gene;
	int location;
	int locationToGene;
	int codingStatus;
};

struct STATplus
{
	double sc;
	double df;
	float rsquare;
	double Fstat;
	//float b[38];
	float *b;
	float bcvsex;  //covariate sex
	//float bcv[10]; //covariate parameters
	float *bcv;
	//int in[38];
	int *in;
	float *h; // haplotype frequnecies, old: h[8]
	/*float betaNew_se[38];
	float betaNew_lcl[38];
	float betaNew_rcl[38];
	float oddsRatio[38];
	float lcloddsRatio[38];
	float rcloddsRatio[38];*/
	float *betaNew_se;
	float *betaNew_lcl;
	float *betaNew_rcl;
	float *oddsRatio;
	float *lcloddsRatio;
	float *rcloddsRatio;
	float *sigma1;
	//float *sigma2;
};

/*struct STATplus2
{

	double sc;
	double df;
	float rsquare;
	double Fstat;
	float *b;
	float bcvsex;  //covariate sex
	float *bcv; //covariate parameters
	int *in;
	float h[8]; // haplotype frequnecies, old: h[8]
	float *betaNew_se;
	float *betaNew_lcl;
	float *betaNew_rcl;
	float *oddsRatio;
	float *lcloddsRatio;
	float *rcloddsRatio;

};*/

/*void freeSTATplus2( STATplus2 result)
    {
	 free(result.b);free(result.bcv);free(result.in);
	 free(result.betaNew_se);free(result.betaNew_lcl);free(result.betaNew_rcl);
	 free(result.oddsRatio);free(result.lcloddsRatio);free(result.rcloddsRatio);
	};*/

void freeSTATplus( STATplus result)
    {
	 free(result.b);free(result.bcv);free(result.in);free(result.h);
	 free(result.betaNew_se);free(result.betaNew_lcl);free(result.betaNew_rcl);
	 free(result.oddsRatio);free(result.lcloddsRatio);free(result.rcloddsRatio);
	 free(result.sigma1);//free(result.sigma2);
	};

struct DETAILS
{
  float aCaN;
	float aCoN;
	float bCaN;
	float bCoN;
	float aCa;
	float bCa;
	float aCo;
	float bCo;
	float orA;
	float orB;
	float lclA;
	float rclA;
	float lclB;
	float rclB;
	float aaCa;
	float abCa;
	float bbCa;
	float aaCo;
	float abCo;
	float bbCo;
	double testHWE_Ca;
	double testHWE_Co;
  };

struct COUNTS // SNP (row)
{
	int AA_Ca; // number of Genotype AA cases
	int AA_Ca_male; // number of Genotype AA cases for X
	int AA_Ca_female; // number of Genotype AA cases for X
	int AB_Ca; // number of Genotype Ab cases
	int AB_Ca_male; // number of Genotype Ab cases for X
	int AB_Ca_female; // number of Genotype Ab cases for X
	int BB_Ca; // number of Genotype BB cases
	int BB_Ca_male; // number of Genotype BB cases for X
	int BB_Ca_female; // number of Genotype BB cases for X
	int AA_Co; // number of Genotype AA control
	int AA_Co_male; // number of Genotype AA control for X
	int AA_Co_female; // number of Genotype AA control for X
	int AB_Co; // number of Genotype AB control
	int AB_Co_male; // number of Genotype AB control for X
	int AB_Co_female; // number of Genotype AB control for X
	int BB_Co; // number of Genotype BB control
	int BB_Co_male; // number of Genotype BB control for X
	int BB_Co_female; // number of Genotype BB control for X
	int OO_Ca; // number of Genotype 00 cases
	int OO_Ca_male; // number of Genotype 00 cases for X
	int OO_Ca_female; // number of Genotype 00 cases for X
	int OO_Co; // number of Genotype 00 control
	int OO_Co_male; // number of Genotype 00 control for X
	int OO_Co_female; // number of Genotype 00 control for X

	struct DETAILS *det;
	struct STATplus *result1;
	double pSingle;
};

struct ALLELECODE
{
	char *a1;
	char *a2;
};

struct GENOCODE
{
	char code[8]; //z.B: 0 0, A A, A B, B A, B B
};

struct REGION
{
	char chr[3];
	signed int begin;
	signed int end;
};

struct BESTSINGLEMARKER
{
	int nr;
	double p;
	double pmod;
};

struct GENETICLIST
{
	int nr;
};

struct OVERLAPLIST
{
	int nr;
};

struct PATHWAY
{
	char name[100];
	int *list;
	int counts;
	int counts0;
	int **listsingle;
	int *listgenetic;
	int **listoverlap;
	int *singletop;
	int *singletop0;
	int genetictop;
	int *overlaptop;
	double **listp; //PWT
	double *score;
	double score0;
	double ratio0;
	int *singletopmod; //now used for number of signficant genes
	int *list_gene;
	int **listsingle_gene;
	int **listp_in; //pplus
	int **listp_in2; //pplus
	int **listp_in3; //pplus
	int ngenes;
	int *intertop;
	int *intertopmod;
};

//pplus
struct PWgenelist
{
	char *gene;
	int nr;
};

//Änderung569
struct MARKERTABLE
{
	char *rs;
	int pos;
};

struct BESTCHI5
{
	int nr1;
	int nr2;
	int nr3;
	int nr4;
	int nr5;
	double p;
	double pmod;
	double Fstat;
	double casecountsMale[3][3][3][3][3];
	double casecountsFemale[3][3][3][3][3];
	double casecounts[3][3][3][3][3];
	double controlcountsMale[3][3][3][3][3];
	double controlcountsFemale[3][3][3][3][3];
	double controlcounts[3][3][3][3][3];
	int r;
	bool X;
	struct STATplus result1;
	//struct STATplus result2;
	struct STATplus result3;
	//struct STATplus result1M;
	//struct STATplus result2M;
	//struct STATplus result1F;
	//struct STATplus result2F;
	struct TRAITAVG traitavg;
};

struct TSTAT
{
	double p;
	double pmod;
	double df;
};


struct TOPLIST
{
	int nr1;
	int nr2;
	int nr3;
	int nr4;
	int nr5;
	double p;
	int r;
};

struct INDEX
{
	int z1;
	int z2;
	int z3;
	int z4;
	int z5;
};

struct RELATIVES {
    int i;
    int j;
};

struct CLUSTERS {
  uint32_t nppls;
  uint32_t ncases;
  uint32_t nctrls;
  uint32_t* list;
};


#include "isnp_binstuff.cpp"
#include "isnp_files.cpp"
//#include "famData.cpp"
#include "isnp_matching.cpp"
#include "regression.cpp"
#include "isnp_out.cpp"  //writeOut
#if RARE
#include "intervalfile.cpp"
#include "rare.cpp"
#endif
string outputname = "IS_";
fstream errorfile; // Error message
fstream logfile; // Logfile

int mafadjust=1;
int binamin=2;
int binamax=0;
uint8_t t_all=0;
uint8_t t_cc=0;

/// Funktionen
void initTRAITavg(struct TRAITAVG *traitavg)
{
	int i,j;
	for(i=0;i<2;i++)
	{
		(*traitavg).snp1allele[i]=-1;(*traitavg).snp2allele[i]=-1;
		for(j=0;j<2;j++)
		{
			(*traitavg).allele[i][j]=-1;
		}
	}
	for(i=0;i<3;i++)
	{
		(*traitavg).snp1geno[i]=-1;(*traitavg).snp2geno[i]=-1;
		for(j=0;j<3;j++)
		{
			(*traitavg).geno[i][j]=-1;
		}
	}
}

void quickSortBESTCHI5(struct BESTCHI5 arr[], int left, int right)
{
	int i=left, j=right;
	struct BESTCHI5 tmp;
	struct BESTCHI5 pivot = arr[(left + right) / 2];
	/* partition */

	while (i <= j)

	{
		while (arr[i].p < pivot.p){i++;}
		while (arr[j].p > pivot.p){j--;}
		if (i <= j)
		{
			tmp = arr[i];
			arr[i] = arr[j];
			arr[j] = tmp;
			i++;j--;
		}
	};

	/* recursion */

	//exit(1);

	if (left < j)
	{
		quickSortBESTCHI5(arr, left, j);
	}


	if (i < right)
	{
		quickSortBESTCHI5(arr, i, right);
	}
	 //exit(1);
}


void init(struct STATplus *result, int haplo, int n, int maxIndexCov, int multi, int needArrays, int dim_L1, int dim_L2)
{
	int i;
	int size=28;

	//if(!multi){size=3;} IS_562 deactivated

	    if(needArrays)
	      {

		(*result).sc=0;
		(*result).df=0;
		(*result).bcvsex=0;

		(*result).in = (int *)calloc(size+maxIndexCov,sizeof(int));
		if (!(*result).in) die("Problem allocating (*result).in in init();");

		(*result).bcv = (float *)malloc(maxIndexCov*sizeof(float));
		if (!(*result).bcv) die("Problem allocating (*result).bcv in init();");
		for(i=0;i<maxIndexCov;i++) (*result).bcv[i]=-1;

		if(1)
		{
			(*result).b = (float *)malloc((size+maxIndexCov)*sizeof(float));
            if (!(*result).b) die("Problem allocating (*result).b in init();");
			for(i=0;i<size+maxIndexCov;i++) (*result).b[i]=-1;
		}

		if(haplo)
		{
			(*result).h = (float *)malloc(8*sizeof(float));
            if (!(*result).h) die("Problem allocating (*result).h in init();");
			for(i=0;i<8;i++) (*result).h[i]=0;
		}

	}

	if( n==0 && 1 && needArrays)
	{
		(*result).betaNew_se   = (float *)malloc((size+maxIndexCov)*sizeof(float));
		(*result).betaNew_lcl  = (float *)malloc((size+maxIndexCov)*sizeof(float));
		(*result).betaNew_rcl  = (float *)malloc((size+maxIndexCov)*sizeof(float));
		(*result).oddsRatio    = (float *)malloc((size+maxIndexCov)*sizeof(float));
		(*result).lcloddsRatio = (float *)malloc((size+maxIndexCov)*sizeof(float));
		(*result).rcloddsRatio = (float *)malloc((size+maxIndexCov)*sizeof(float));
		if (!(*result).betaNew_se  ) die("Problem allocating (*result).betaNew_se in init();");
		if (!(*result).betaNew_lcl ) die("Problem allocating (*result).betaNew_lcl in init();");
		if (!(*result).betaNew_rcl ) die("Problem allocating (*result).ibetaNew_rcl in init();");
		if (!(*result).oddsRatio   ) die("Problem allocating (*result).oddsRatio in init();");
		if (!(*result).lcloddsRatio) die("Problem allocating (*result).lcloddsRatio in init();");
		if (!(*result).rcloddsRatio) die("Problem allocating (*result).rcloddsRatio in init();");
		for(i=0;i<size+maxIndexCov;i++)
		{
			(*result).betaNew_se[i]  =-1;
			(*result).betaNew_lcl[i] =-1;
			(*result).betaNew_rcl[i] =-1;
			(*result).oddsRatio[i]   =-1;
			(*result).lcloddsRatio[i]=-1;
			(*result).rcloddsRatio[i]=-1;
		}

		(*result).sigma1   = (float *)malloc((dim_L1)*sizeof(float));
		// (*result).sigma2   = (float *)malloc((dim_L2)*sizeof(float));

	}

}

/*
void resultCopy(struct STATplus *result1,struct STATplus result2)
{
	(*result1).betaNew_se=result2.betaNew_se;
	(*result1).betaNew_lcl=result2.betaNew_lcl;
	(*result1).betaNew_rcl=result2.betaNew_rcl;
	(*result1).oddsRatio=result2.oddsRatio;
	(*result1).lcloddsRatio=result2.lcloddsRatio;
	(*result1).rcloddsRatio=result2.rcloddsRatio;
	(*result1).h=result2.h;
	(*result1).b=result2.b;
	(*result1).in=result2.in;
	(*result1).bcv=result2.bcv;
	(*result1).sigma1=result2.sigma1;
	(*result1).sigma2=result2.sigma2;
}*/


/*void Plus2ToPlus(struct STATplus *result1,struct STATplus result2, int haplo, int maxIndexCov, int printBeta) //TIM_NEW
{
	int i;

	if(haplo)
	{
	    (*result1).h = (float *)realloc((*result1).h,8*sizeof(float));
		for(i=0;i<8;i++)
		{
			(*result1).h[i]=result2.h[i];
		}
	}

	if(printBeta)
	  {

	    (*result1).betaNew_se = (float *)realloc((*result1).betaNew_se,(28+maxIndexCov)*sizeof(float));
	    (*result1).betaNew_lcl = (float *)realloc((*result1).betaNew_lcl,(28+maxIndexCov)*sizeof(float));
	    (*result1).betaNew_rcl = (float *)realloc((*result1).betaNew_rcl,(28+maxIndexCov)*sizeof(float));
	    (*result1).oddsRatio = (float *)realloc((*result1).oddsRatio,(28+maxIndexCov)*sizeof(float));
	    (*result1).lcloddsRatio = (float *)realloc((*result1).lcloddsRatio,(28+maxIndexCov)*sizeof(float));
	    (*result1).rcloddsRatio = (float *)realloc((*result1).rcloddsRatio,(28+maxIndexCov)*sizeof(float));
	    (*result1).b = (float *)realloc((*result1).b,(28+maxIndexCov)*sizeof(float));
	    (*result1).bcv = (float *)realloc((*result1).bcv,(28+maxIndexCov)*sizeof(float));
	    (*result1).in = (int *)realloc((*result1).in,(28+maxIndexCov)*sizeof(int));


		(*result1).sc=result2.sc;
		(*result1).df=result2.df;
		(*result1).rsquare=result2.rsquare;
		(*result1).Fstat=result2.Fstat;
		(*result1).bcvsex=result2.bcvsex;


		for(i=0;i<maxIndexCov;i++)
		{
			(*result1).bcv[i]=result2.bcv[i];
		}

		for(i=0;i<28+maxIndexCov;i++)
		{
			(*result1).b[i]=result2.b[i];
			(*result1).in[i]=result2.in[i];
			(*result1).betaNew_se[i]=result2.betaNew_se[i];
			(*result1).betaNew_lcl[i]=result2.betaNew_lcl[i];
			(*result1).betaNew_rcl[i]=result2.betaNew_rcl[i];
			(*result1).oddsRatio[i]=result2.oddsRatio[i];
			(*result1).lcloddsRatio[i]=result2.lcloddsRatio[i];
			(*result1).rcloddsRatio[i]=result2.rcloddsRatio[i];
		}
	}
}

*/

void PlusToPlus(struct STATplus *result1,struct STATplus result2, int haplo, int maxIndexCov, int printBeta, int multi, int covariancematrix, int covdim1,int covdim2) //TIM_NEW
{
	int i;
	int size=28;
	if(!multi){size=3;}

	if(haplo)
	{
	    (*result1).h = (float *)realloc((*result1).h,8*sizeof(float));
		for(i=0;i<8;i++)
		{
			(*result1).h[i]=result2.h[i];
		}
	}

	if(covariancematrix)
	{
	    (*result1).sigma1 = (float *)realloc((*result1).sigma1,covdim1*sizeof(float));
		// (*result1).sigma2 = (float *)realloc((*result1).sigma2,covdim2*sizeof(float));
        for(i=0;i<covdim1;i++)
		{
		 (*result1).sigma1[i]=result2.sigma1[i];
		}
        /*for(i=0;i<covdim2;i++)
		{
		 (*result1).sigma2[i]=result2.sigma2[i];
		}**/
	}

	(*result1).in = (int *)realloc((*result1).in,(28+maxIndexCov)*sizeof(int));
	for(i=0;i<28+maxIndexCov;i++)
		{
		 (*result1).in[i]=result2.in[i];
		}

	if(printBeta)
	  {
	    (*result1).betaNew_se = (float *)realloc((*result1).betaNew_se,(size+maxIndexCov)*sizeof(float));
	    (*result1).betaNew_lcl = (float *)realloc((*result1).betaNew_lcl,(size+maxIndexCov)*sizeof(float));
	    (*result1).betaNew_rcl = (float *)realloc((*result1).betaNew_rcl,(size+maxIndexCov)*sizeof(float));
	    (*result1).oddsRatio = (float *)realloc((*result1).oddsRatio,(size+maxIndexCov)*sizeof(float));
	    (*result1).lcloddsRatio = (float *)realloc((*result1).lcloddsRatio,(size+maxIndexCov)*sizeof(float));
	    (*result1).rcloddsRatio = (float *)realloc((*result1).rcloddsRatio,(size+maxIndexCov)*sizeof(float));
	    (*result1).b = (float *)realloc((*result1).b,(size+maxIndexCov)*sizeof(float));
	    (*result1).bcv = (float *)realloc((*result1).bcv,(size+maxIndexCov)*sizeof(float));

		(*result1).sc=result2.sc;
		(*result1).df=result2.df;
		(*result1).rsquare=result2.rsquare;
		(*result1).Fstat=result2.Fstat;
		(*result1).bcvsex=result2.bcvsex;


		for(i=0;i<maxIndexCov;i++)
		{
			(*result1).bcv[i]=result2.bcv[i];
		}

		for(i=0;i<size+maxIndexCov;i++)
		{
			(*result1).b[i]=result2.b[i];
			(*result1).in[i]=result2.in[i];
			(*result1).betaNew_se[i]=result2.betaNew_se[i];
			(*result1).betaNew_lcl[i]=result2.betaNew_lcl[i];
			(*result1).betaNew_rcl[i]=result2.betaNew_rcl[i];
			(*result1).oddsRatio[i]=result2.oddsRatio[i];
			(*result1).lcloddsRatio[i]=result2.lcloddsRatio[i];
			(*result1).rcloddsRatio[i]=result2.rcloddsRatio[i];
		}
	}
}


double Sidak(double p, double n)
{
    double x;
    x=1-pow(1-p,n);

    if(x > 0) {return x;}
    else {return n*p;} //if Sidak does not work: Bonferroni
}


int update(string *mystring,int a)
{
	int i=0,j;
	char help[12]="0123456789";

	if(a==200){*mystring+='2';*mystring+='0';*mystring+='0';return i;}
	if(a>=200){*mystring+='2';a-=200;}
	if(a==100){*mystring+='1';*mystring+='0';*mystring+='0';return i;}
	if(a>=100){*mystring+='1';a-=100;}

	for (i=90;i>=10;i-=10)
	{
		j=i/10;
		if(a==i){*mystring+=help[j];*mystring+='0';return i;}
		if(a>=i){*mystring+=help[j];a-=i;}
	}

	for (i=9;i>=0;i-=1)
	{
		if(a>=i){*mystring+=help[i];return i;}
	}
	return 0;
}

// Function to calculate the p-value
double lnpValueCalc(double xx)
{
	int j;
	double x, tmp, y, ser;
	double cof[14] = { 	57.1562356658629235, -59.5979603554754912, 14.1360979747417471, -0.491913816097620199, .339946499848118887e-4,
						.465236289270485756e-4, -.983744753048795646e-4, .158088703224912494e-3, -.210264441724104883e-3,
						.217439618115212643e-3, -.164318106536763890e-3, .844182239838527433e-4, -.261908384015814087e-4,
						.368991826595316234e-5 };
	if (xx <= 0)
	{
		xx = 0;
	}
	y = x = xx;
	tmp = x + 5.24218750000000000;
	tmp = (x + 0.5) * log(tmp) - tmp;
	ser = 0.999999999999997092;
	for (j = 0; j < 14; j++)
	        ser += cof[j] / ++y;
	return tmp + log(2.5066282746310005 * ser / x);
}

// Function to calculate the p-value
double pValueCalc(double a, double x)
{
	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	double gamser, gammcf, gln;

	if (x < 0.0 || a <= 0.0)
	{
	  a = 2;
	  x = 0;
	}
	if (x < (a + 1.0))
	{
		gser(&gamser, a, x, &gln);
		//		cout<<"gamser "<<1.0-gamser<<endl;
		return 1.0 - gamser;
	}
	else
	{
		gcf(&gammcf, a, x, &gln);
		//		cout<<"gammcf "<<1.0-gammcf<<endl;
		return gammcf;
	}
}

double MAX(double a, double b)
{
	if(a>=b){return a;}
	else {return b;}
}

double MIN(double a, double b)
{
    if(a<=b){return a;}
    else {return b;}
}

double SQR(double a)
{
    return a*a;
}

double rounddown(double a, int precision)
{
  return ((double)floor(a*pow(10,(double)precision))/(pow(10,(double)precision)));
}



double betacf(double a, double b, double x)
{
	int maxit=10000;
	double eps1=3.0e-10;
	double fpmin1=1.0e-30;

	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;
	qab=a+b;
	qap=a+1;
	qam=a-1;
	c=1;
	d=1-qab*x/qap;
	if (fabs(d) < fpmin1) {d=fpmin1;}
	d=1/d;
	h=d;
	for (m=1;m<maxit;m++)
	{
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1+aa*d;
		if (fabs(d) < fpmin1) {d=fpmin1;}
			c=1+aa/c;
		if (fabs(c) < fpmin1) {c=fpmin1;}
			d=1/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1+aa*d;
		if (fabs(d) < fpmin1) {d=fpmin1;}
			c=1+aa/c;
		if (fabs(c) < fpmin1) {c=fpmin1;}
			d=1/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) <= eps1) break;
	}
	return h;
}


double betaiapprox(double a, double b, double x)
{
	int j;
	double xu,t,sum,ans;
	double a1 = a-1.0, b1 = b-1.0, mu = a/(a+b);
	double lnmu=log(mu),lnmuc=log(1.-mu);
	double y[18] = {0.0021695375159141994,
					0.011413521097787704,0.027972308950302116,0.051727015600492421,
					0.082502225484340941, 0.12007019910960293,0.16415283300752470,
					0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
					0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
					0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
					0.87126389619061517, 0.95698180152629142};
	double w[18] = {0.0055657196642445571,
					0.012915947284065419,0.020181515297735382,0.027298621498568734,
					0.034213810770299537,0.040875750923643261,0.047235083490265582,
					0.053244713977759692,0.058860144245324798,0.064039797355015485,
					0.068745323835736408,0.072941885005653087,0.076598410645870640,
					0.079687828912071670,0.082187266704339706,0.084078218979661945,
					0.085346685739338721,0.085983275670394821};

	t = sqrt(a*b/(SQR(a+b)*(a+b+1.0)));
	if (x > a/(a+b))
	{
		if (x >= 1.0) return 1.0;
		xu = MIN(1.,MAX(mu + 10.*t, x + 5.0*t));
	}
	else
	{
		if (x <= 0.0) return 0.0;
		xu = MAX(0.,MIN(mu - 10.*t, x - 5.0*t));
	}
	sum = 0;
	for (j=0;j<18;j++)
	{
		t = x + (xu-x)*y[j];
		sum += w[j]*exp(a1*(log(t)-lnmu)+b1*(log(1-t)-lnmuc));
	}
	ans = sum*(xu-x)*exp(a1*lnmu-lnpValueCalc(a)+b1*lnmuc-lnpValueCalc(b)+lnpValueCalc(a+b));
	return ans>0.0? 1.0-ans : -ans;
}


double betai(double a, double b, double x)
{
	double bt;
	int SWITCH=3000;

	if (a <= 0.0 || b <= 0.0) {x=1;}
	if (x < 0.0 || x > 1.0) {x=1;}
	if (x == 0.0 || x == 1.0) return x;
	if (a > SWITCH && b > SWITCH) return betaiapprox(a,b,x);
	bt=exp(lnpValueCalc(a+b)-lnpValueCalc(a)-lnpValueCalc(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a;
	else return 1.0-bt*betacf(b,a,1.0-x)/b;
}


// Function to calculate the p-value
void gser(double *gamser, double a, double x, double *gln)
{
	double lnpValueCalc(double xx);
	int n;
	double sum, del, ap;

	*gln = lnpValueCalc(a);
	if (x <= 0.0)
	{
		//if (x < 0.0)
		      //cout << "x less than 0 in routine gser";
		*gamser = 0.0;
		return;
	}
	else
	{
		ap = a;
		del = sum = 1.0 / a;

		for (n = 1; n <= ITMAX; n++)
		{
			++ap;
			del *= x / ap;
			sum += del;
			if (fabs(del) < fabs(sum) * EPS)
			{
				*gamser = sum * exp(-x + a * log(x) - (*gln));
				return;
			}
		}
		//cout << "a too large, ITMAX too small in routine gser";
		return;
	}
}


// Function to calculate the p-value
void gcf(double *gammcf, double a, double x, double *gln)
{
	double lnpValueCalc(double xx);
	int i;
	double an, b, c, d, del, h;

	*gln = lnpValueCalc(a);
	b = x + 1.0 - a;
	c = 1.0 / FPMIN;
	d = 1.0 / b;
	h = d;
	for (i = 1; i <= ITMAX; i++)
	{
		an = -i * (i - a);
		b += 2.0;
		d = an * d + b;

		if (fabs(d) < FPMIN)
		{
			d = FPMIN;
		}

		c = b + an / c;

		if (fabs(c) < FPMIN)
		{
			c = FPMIN;
		}

		d = 1.0 / d;
		del = d * c;
		h *= del;
		if (fabs(del - 1.0) < EPS)
		{
			break;
		}
	}
	//if (i > ITMAX)
	//cout << "a too large, ITMAX too small in gcf";
	*gammcf = exp(-x + a * log(x) - (*gln)) * h;
}


// 95%-quantile
double tquantile(int n)
{
	double quantile=1.96;

	if(n<=200){quantile=1.97;}else {return quantile;}
	if (n<=150){quantile=1.98;}else {return quantile;}
	if (n<=80){quantile=1.99;}else {return quantile;}
	if (n<=60){quantile=2.00;}else {return quantile;}
	if (n<=50){quantile=2.01;}else {return quantile;}
	if (n<=44){quantile=2.02;}else {return quantile;}
	if (n<=34){quantile=2.03;}else {return quantile;}
	if (n<=30){quantile=2.04;}else {return quantile;}
	if (n<=29){quantile=2.05;}else {return quantile;}
	if (n<=27){quantile=2.06;}else {return quantile;}
	if (n<=23){quantile=2.07;}else {return quantile;}
	if (n<=21){quantile=2.08;}else {return quantile;}
	if (n<=20){quantile=2.09;}else {return quantile;}
	if (n<=18){quantile=2.10;}else {return quantile;}
	if (n<=17){quantile=2.11;}else {return quantile;}
	if (n<=16){quantile=2.12;}else {return quantile;}
	if (n<=15){quantile=2.13;}else {return quantile;}
	if (n<=14){quantile=2.14;}else {return quantile;}
	if (n<=13){quantile=2.16;}else {return quantile;}
	if (n<=12){quantile=2.18;}else {return quantile;}
	if (n<=11){quantile=2.20;}else {return quantile;}
	if (n<=10){quantile=2.23;}else {return quantile;}
	if (n<=9){quantile=2.26;}else {return quantile;}
	if (n<=8){quantile=2.31;}else {return quantile;}
	if (n<=7){quantile=2.36;}else {return quantile;}
	if (n<=6){quantile=2.45;}else {return quantile;}
	if (n<=5){quantile=2.57;}else {return quantile;}
	if (n<=4){quantile=2.78;}else {return quantile;}
	if (n<=3){quantile=3.18;}else {return quantile;}
	if (n<=2){quantile=4.30;}else {return quantile;}
	if (n<=1){quantile=12.71;}else {return quantile;}
	return quantile;
}


// HWE-test
double testHWE(int AA, int AB, int BB)
{
	double e2;
	double e1;
	double e0;
	double T;
	double p=1;
	double chisq = 0;
	double pwert;

	T = 2*(AA+AB+BB);
	if (T > 0)
	{
		p = (2*AA+AB)/T;
	}

	// expected values
	if (T > 0)
	{
		e2 = pow((p), 2) * 0.5 * T;
		e1 = 2 * p * (1 - p) * 0.5 * T;
		e0 = pow((1 - p), 2) * 0.5 * T;

		if (e2 != 0)
		{
			chisq += (pow((AA - e2), 2)) / e2;
		}
		if (e1 != 0)
		{
			chisq += (pow((AB - e1), 2)) / e1;
		}
		if (e0 != 0)
		{
			chisq += (pow((BB - e0), 2)) / e0;
		}

		//p-value chisq
		if (chisq != 0 && std::isnan(chisq) == 0)
		{
			pwert = pValueCalc(0.5, chisq / 2); // 1. value DF/2; 2. value teststat/2
		}
		else
		{
			pwert = 1;
		}
	}
	else
	{
		pwert = 1;
	}
	return pwert;
}


// Armitage-test for autosomes and pseudo-autosomal region of Y chromosome
double armitageTest(struct COUNTS counts, int teststat, double *inflationfactor, double *fulltests)
{
	double n_Ca;
	double n_Co;
	double n;
	double n1;
	double n2;
	double a;
	double b;
	double tTrend;
	double pValue;
	double n0;

	n_Ca = counts.AA_Ca + counts.AB_Ca + counts.BB_Ca;
	n_Co = counts.AA_Co + counts.AB_Co + counts.BB_Co;
	n = n_Ca + n_Co;
	n1 = counts.AB_Ca + counts.AB_Co;
	n2 = counts.BB_Ca + counts.BB_Co;
	a = counts.AB_Ca;
	b = counts.BB_Ca;
	n0 = n - n1 - n2;

	// Armitage-trendtest:
	tTrend = (n * (pow((n * (a + 2 * b) - n_Ca * (n1 + 2 * n2)), 2))) /
			(n_Ca * n_Co * (n * (n1 + 4 * n2) - pow((n1 + 2 * n2), 2)));

	//p-value chisq

	if ((n0 > 0 && n1 > 0) || ((n0 && n2) || (n1 && n2)))
	{
		if (std::isnan(tTrend) != 1)
		{
#pragma omp critical(INFLATIONFACTOR)
          {
			(*inflationfactor) += tTrend;
			(*fulltests) += 1;
          }
		}
	}

	if (teststat == 0)
	{
		if (tTrend != 0 && std::isnan(tTrend) != 1)
		{
			pValue = pValueCalc(0.5, tTrend / 2); // 1. value DF/2; 2. value teststat/2
		}
		else
		{
			pValue = 1;
		}

		return pValue;
	}
	else if (tTrend != 0 && std::isnan(tTrend) != 1)
	{
		return tTrend;
	}
	else
	{
		return 0;
	}
}


// Singlemarker analysis (chi2-test) for X chromosome
double snpTestX(struct COUNTS counts, int teststat,	double *inflationfactor, double *fulltests)
{
	// Create a 2x2 table of allele counts in cases and controls

	double xa;
	double xb;
	double n;
	double xF;
	double fCa;
	double fCo;
	double xM;
	double mCa;
	double mCo;
	double Y;
	double A;
	double P;
	double uA;
	double vM;
	double vF;
	double V;

	double tX = 0;
	double pValue;
	double nmod=0;

	xa = counts.AA_Ca_male + 2 * counts.AA_Ca_female + counts.AB_Ca_female;
	xb = counts.AA_Co_male + 2 * counts.AA_Co_female + counts.AB_Co_female;


	n = counts.AA_Ca_male + counts.BB_Ca_male + counts.AA_Co_male + counts.BB_Co_male
	    + counts.AA_Ca_female + counts.AB_Ca_female + counts.BB_Ca_female
	    + counts.AA_Co_female + counts.AB_Co_female + counts.BB_Co_female;

	nmod = counts.AA_Ca_male + counts.BB_Ca_male + counts.AA_Co_male + counts.BB_Co_male
			+ 2*counts.AA_Ca_female + 2*counts.AB_Ca_female + 2*counts.BB_Ca_female
			+ 2*counts.AA_Co_female + 2*counts.AB_Co_female + 2*counts.BB_Co_female;


	xF = counts.AA_Ca_female + counts.AB_Ca_female + counts.BB_Ca_female
	    + counts.AA_Co_female + counts.AB_Co_female + counts.BB_Co_female;

	fCa = counts.AA_Ca_female + counts.AB_Ca_female + counts.BB_Ca_female;
	fCo = xF - fCa;

	xM = n - xF;
	mCa = counts.AA_Ca_male + counts.BB_Ca_male;
	mCo = xM - mCa;

	Y = (counts.AA_Ca_male + counts.BB_Ca_male + counts.AA_Ca_female
	    + counts.AB_Ca_female + counts.BB_Ca_female)/n;

	A = (2*counts.AA_Ca_male + 2*counts.AA_Ca_female + counts.AB_Ca_female+2*counts.AA_Co_male + 2*counts.AA_Co_female + counts.AB_Co_female)/n;

	P = 1 -(xa + xb)/(nmod);

	uA = counts.AA_Ca_male * (1-Y) * 2
      	+ counts.AA_Co_male * (-Y) * 2
        + counts.AA_Ca_female * (1-Y) * 2
        + counts.AA_Co_female * (-Y) * 2
        + counts.AB_Ca_female * (1-Y)
        + counts.AB_Co_female * (-Y);

	vM = 4*P*(1-P);

	if (xF <=1)
	{
		vF = 0;
	}
	else
	{
		vF = (1/(xF-1)) * ((counts.AA_Ca_female * (2-A)* (2-A)) + (counts.AB_Ca_female * (1-A)* (1-A))
			+ (counts.BB_Ca_female * (-A)*(-A)) + (counts.AA_Co_female * (2-A)* (2-A))
			+ (counts.AB_Co_female * (1-A)* (1-A)) + (counts.BB_Co_female * (-A)* (-A)));
	}

	V = vF*(fCa*(1-Y)*(1-Y) + fCo*(-Y)*(-Y))+ vM*(mCa*(1-Y)*(1-Y)+mCo*(-Y)*(-Y));

	if(V>0){tX = uA * uA/V;}
	else
	{
		tX=0;
	}



	//calculate p-value
	if (teststat == 0)
	{
		if( xa >0 && xb >0 && V>0 && std::isnan(tX) == 0)
	    {
			(*inflationfactor) +=tX;(*fulltests)+=1;
		}


		if (tX != 0 && std::isnan(tX) == 0)
		{
			pValue = pValueCalc(0.5, tX / 2); // 1. value DF/2; 2. value teststat/2
		}
		else
		{
			pValue = 1;
		}
		return pValue;
	}
	else if (tX != 0 && std::isnan(tX) == 0)
	{
		return tX;
	}
	else
	{
		return 0;
	}
}


// Genotyp-test
double genotypTest(struct COUNTS counts, int teststat,	double *inflationfactor, double *fulltests)
{
	double expectedAA_Ca;
	double expectedAB_Ca;
	double expectedBB_Ca;
	double expectedAA_Co;
	double expectedAB_Co;
	double expectedBB_Co;
	double sumCases = 0;
	double sumControls = 0;
	double stat = 0;
	double pValue;

	double df = 0;

	sumCases = counts.AA_Ca + counts.AB_Ca + counts.BB_Ca;
	sumControls = counts.AA_Co + counts.AB_Co + counts.BB_Co;

	// expected cases and controls
	expectedAA_Ca = (sumCases * (counts.AA_Ca + counts.AA_Co)) / (sumCases + sumControls);
	expectedAB_Ca = (sumCases * (counts.AB_Ca + counts.AB_Co)) / (sumCases + sumControls);
	expectedBB_Ca = (sumCases * (counts.BB_Ca + counts.BB_Co)) / (sumCases + sumControls);

	expectedAA_Co = (sumControls * (counts.AA_Ca + counts.AA_Co)) / (sumCases + sumControls);
	expectedAB_Co = (sumControls * (counts.AB_Ca + counts.AB_Co)) / (sumCases + sumControls);
	expectedBB_Co = (sumControls * (counts.BB_Ca + counts.BB_Co)) / (sumCases + sumControls);

	if (expectedAA_Ca > 0 || expectedAA_Co > 0)
	{
		if (expectedAA_Ca > 0)
		{
			stat = stat + ((counts.AA_Ca - expectedAA_Ca) * (counts.AA_Ca - expectedAA_Ca)) / (expectedAA_Ca);
			df++;
		}
		if (expectedAA_Co > 0)
		{
			stat = stat + ((counts.AA_Co - expectedAA_Co) * (counts.AA_Co - expectedAA_Co)) / (expectedAA_Co);
		}
	}

	if (expectedAB_Ca > 0 || expectedAB_Co > 0)
	{
		if (expectedAB_Ca > 0)
		{
			stat = stat + ((counts.AB_Ca - expectedAB_Ca) * (counts.AB_Ca - expectedAB_Ca)) / (expectedAB_Ca);
			df++;
		}
		if (expectedAB_Co > 0)
		{
			stat = stat + ((counts.AB_Co - expectedAB_Co) * (counts.AB_Co - expectedAB_Co)) / (expectedAB_Co);
		}
	}

	if (expectedBB_Ca > 0 || expectedBB_Co > 0)
	{
		if (expectedBB_Ca > 0)
		{
			stat = stat + ((counts.BB_Ca - expectedBB_Ca) * (counts.BB_Ca - expectedBB_Ca)) / (expectedBB_Ca) ;
			df++;
		}
		if (expectedBB_Co > 0)
		{
			stat = stat + ((counts.BB_Co - expectedBB_Co) * (counts.BB_Co - expectedBB_Co)) / (expectedBB_Co);
		}
	}

	df -= 1;

	if (teststat == 1)
	{
		if (df == 2)
		{
			return stat;
		}
		else
		{
			return stat * 1.4;
		}
	}
	else
	{
		if(df==2){ (*inflationfactor) +=stat;(*fulltests)+=1;}
		else if(df==1){ (*inflationfactor) += (2*stat);(*fulltests)+=1;}

		if (stat > 0)
		{
			pValue = pValueCalc(df / 2, stat / 2);
		}
		else
		{
			pValue = 1;
		}
		return pValue;
	}
}


// Separate the input of  POSCHOICE, NEGCHOICE and COVARIATES
vector<string> split(string str, string delim)
{
	int pos = 0;
	vector<string> results;

	while ((pos = str.find_first_of(delim)) != str.npos)
	{
		if (pos > 0)
		{
			results.push_back(str.substr(0, pos));
		}
		str = str.substr(pos + 1);
	}
	if (str.length() > 0)
	{
		results.push_back(str);
	}
	return results;
}


// Sort best of the singlemarker values
void insertion(struct BESTSINGLEMARKER *bestsinglemarker, int lastPos)
{
	int j;
	double currentValue;
	struct BESTSINGLEMARKER currentStructValue;

	// aktuelle Position zwischenspeichern
	currentStructValue = bestsinglemarker[lastPos];
	currentValue = bestsinglemarker[lastPos].p;


	// Kleineren Wert als currentValue suchen. Schleife
	// durchläuft von aktueller Pos. von rechts nach links
	for (j = lastPos; (j > 0) && (bestsinglemarker[j - 1].p > currentValue); j--)
	{
		//Wenn Vorgänger größer als aktuelles Element in currentValue
		bestsinglemarker[j] = bestsinglemarker[j - 1];
	}
	// gespeichertes Element an neue Position -> Lücke auffüllen
	bestsinglemarker[j] = currentStructValue;
}


struct TSTAT chiTest5(double casecounts5[3][3][3][3][3], double controlcounts5[3][3][3][3][3],
                      int dim, int mWithSingletop,double helpstat)
{
	double expectedCases[3][3][3][3][3];
	double expectedControls[3][3][3][3][3];
	double sumCases = 0;
	double sumControls = 0;
	double stat = 0;
	double statNew = 0;
	int a, b, c, d, e;
	double df = 0;
	int stopC = 3;
	int stopD = 3;
	int stopE = 3;

	struct TSTAT tstat;
	tstat.p = 0;
	tstat.pmod = 0;

	if (dim == 2)
	{
		stopC = 1;
		stopD = 1;
		stopE = 1;
	}
	else if (dim == 3)
	{
		stopD = 1;
		stopE = 1;
	}
	else if (dim == 4)
	{
		stopE = 1;
	}

	// Summe Fälle und Kontrollen
	for (a = 0; a < 3; a++)
	{
		for (b = 0; b < 3; b++)
		{
			for (c = 0; c < stopC; c++)
			{
				for (d = 0; d < stopD; d++)
				{
					for (e = 0; e < stopE; e++)
					{
						// Fälle
						sumCases = sumCases + casecounts5[a][b][c][d][e];

						expectedCases[a][b][c][d][e] = 0;
						// Kontrolls
						sumControls = sumControls + controlcounts5[a][b][c][d][e];
						expectedControls[a][b][c][d][e] = 0;
					}
				}
			}
		}
	}

	// Fälle und Kontrollen erwartet
	for (a = 0; a < 3; a++)
	{
		for (b = 0; b < 3; b++)
		{
			for (c = 0; c < stopC; c++)
			{
				for (d = 0; d < stopD; d++)
				{
					for (e = 0; e < stopE; e++)
					{
						expectedCases[a][b][c][d][e] = (sumCases* (casecounts5[a][b][c][d][e] + controlcounts5[a][b][c][d][e]))/ (sumCases + sumControls);
						expectedControls[a][b][c][d][e] = (sumControls* (casecounts5[a][b][c][d][e] + controlcounts5[a][b][c][d][e]))/ (sumCases + sumControls);
					}
				}
			}
		}
	}

	for (a = 0; a < 3; a++)
	{
		for (b = 0; b < 3; b++)
		{
			for (c = 0; c < stopC; c++)
			{
				for (d = 0; d < stopD; d++)
				{
					for (e = 0; e < stopE; e++)
					{
						if (expectedCases[a][b][c][d][e] > 0 || expectedControls[a][b][c][d][e] > 0)
						{
							if (expectedCases[a][b][c][d][e] > 0)
							{
								stat= stat+ (((casecounts5[a][b][c][d][e]- expectedCases[a][b][c][d][e])* (casecounts5[a][b][c][d][e]- expectedCases[a][b][c][d][e]))/ (expectedCases[a][b][c][d][e]));
								df++;

							}
							if (expectedControls[a][b][c][d][e] > 0)
							{
								stat= stat+ (((controlcounts5[a][b][c][d][e]- expectedControls[a][b][c][d][e])* (controlcounts5[a][b][c][d][e]- expectedControls[a][b][c][d][e]))/ (expectedControls[a][b][c][d][e]));
							}
						}
						else
						{
							stat = stat + 0;
						}
					}
				}
			}
		}
	}

	if (stat != 0 && df > 0)
	{
		tstat.p = pValueCalc((df-1) / 2, stat / 2);
	}
	else
	{
		tstat.p = 1;
	}

	if (mWithSingletop == 0)
	{
		tstat.pmod = tstat.p;
	}
	else
	{
		statNew = stat - (mWithSingletop * helpstat);
		if (statNew < 0)
		{
			statNew = 0;
		}
		tstat.pmod = pValueCalc((df-1) / 2, statNew / 2);
	}
	return tstat;
}


// Chi-Test bei Interaktion für 3 Marker
struct TSTAT chiTestInter5(double casecounts5[3][3][3][3][3],double controlcounts5[3][3][3][3][3],
                          	int dim)
{
	double mAlt[3][3][3][3][3][2];
	double mNeu[3][3][3][3][3][2];
	double s1 = 0;
	double s2 = 0;
	double statAlt = 0;
	double statNeu = 0;
	int a, b, c, d, e, k, l;
	double df = 8;
	int stopC = 3;
	int stopD = 3;
	int stopE = 3;
	int merkab = 0;
	int merkac = 0;
	int merkbc = 0;
	int merka = 0;
	int merkb = 0;
	int merkc = 0;

	struct TSTAT tstat;
	tstat.p = 0;
	tstat.pmod = 0;

	if (dim == 2)
	{
		stopC = 1;
		stopD = 1;
		stopE = 1;
	}

	if (dim == 3)
	{
		stopD = 1;
		stopE = 1;
	}
	else if (dim == 4)
	{
		stopE = 1;
	}

	//update df
	for (a = 0; a < 3; a++)
	{
		for (b = 0; b < 3; b++)
		{
			if (casecounts5[a][b][0][0][0] == 0 && casecounts5[a][b][1][0][0]== 0
				&& casecounts5[a][b][2][0][0] == 0 && controlcounts5[a][b][0][0][0] == 0
				&& controlcounts5[a][b][1][0][0] == 0 && controlcounts5[a][b][2][0][0] == 0)
			{
				merkab++;
			}
		}
	}

	for (a = 0; a < 3; a++)
	{
		for (c = 0; c < stopC; c++)
		{
			if ((casecounts5[a][0][c][0][0] == 0 && casecounts5[a][1][c][0][0] == 0
				&& casecounts5[a][2][c][0][0] == 0) || (controlcounts5[a][0][c][0][0] == 0
				&& controlcounts5[a][1][c][0][0] == 0 && controlcounts5[a][2][c][0][0] == 0))
			{
				merkac++;
			}
		}
	}

	for (b = 0; b < 3; b++)
	{
		for (c = 0; c < stopC; c++)
		{
			if ((casecounts5[0][b][c][0][0] == 0 && casecounts5[1][b][c][0][0]== 0
				&& casecounts5[2][b][c][0][0] == 0) || (controlcounts5[0][b][c][0][0] == 0
				&& controlcounts5[1][b][c][0][0] == 0 && controlcounts5[2][b][c][0][0] == 0))
			{
				merkbc++;
			}
		}
	}

	merka = merkab + merkac;
	merkb = merkab + merkbc;
	merkc = merkac + merkbc;

	if (dim == 3)
	{
		if (merka >= 2 || merkb >= 2 || merkc >= 2)
		{
			df = 0;
		}
		else if ((merka + merkb + merkc) == 3)
		{
			df = 1;
		}
		else if ((merka + merkb + merkc) == 2)
		{
			df = 2;
		}
		else if ((merka + merkb + merkc) == 1)
		{
			df = 4;
		}
		else
		{
			df = 8;
		}
	}
	else if (dim == 2)
	{
		if (df - merkab == 8)
		{
			df = 4;
		}
		else if (df - merkab == 7)
		{
			df = 3;
		}
		else if  (df - merkab == 6)
		{
			df = 3;
		}
		else if  (df - merkab == 5)
		{
			df = 2;
		}
		else if  (df - merkab == 4)
		{
			df = 2;
		}
		else if  (df - merkab == 3)
		{
			df = 1;
		}
		else
		{
			df = 0;
		}
	}

	if (df > 0)
	{
		// Startvalues
		for (a = 0; a < 3; a++)
		{
			for (b = 0; b < 3; b++)
			{
				for (c = 0; c < stopC; c++)
				{
					for (d = 0; d < stopD; d++)
					{
						for (e = 0; e < stopE; e++)
						{
							for (k = 0; k < 2; k++)
							{
								mAlt[a][b][c][d][e][k] = 1;
							}//k
						}
					}
				}
			}//b
		} //end a startvalues

		for (l = 1; l < 1000; l++)
		{
			statNeu = 0;
			//Compute mNeu
			for (a = 0; a < 3; a++)
			{
				for (b = 0; b < 3; b++)
				{
					for (c = 0; c < stopC; c++)
					{
						for (d = 0; d < stopD; d++)
						{
							for (e = 0; e < stopE; e++)
							{
								for (k = 0; k < 2; k++)
								{
									if (dim == 2)
									{
										if (l % 3 == 1)
										{
											s1 = casecounts5[a][b][c][d][e] + controlcounts5[a][b][c][d][e];
											s2 = mAlt[a][b][c][d][e][0] + mAlt[a][b][c][d][e][1];
										}
										else if (l % 3 == 2)
										{
											if (k == 0)
											{
												s1 = casecounts5[a][0][c][d][e] + casecounts5[a][1][c][d][e] + casecounts5[a][2][c][d][e];
											}
											else
											{
												s1 = controlcounts5[a][0][c][d][e] + controlcounts5[a][1][c][d][e] + controlcounts5[a][2][c][d][e];
											}
											s2 = mAlt[a][0][c][d][e][k] + mAlt[a][1][c][d][e][k] + mAlt[a][2][c][d][e][k];

										}
										else if (l % 3 == 0)
										{
											if (k == 0)
											{
												s1 = casecounts5[0][b][c][d][e] + casecounts5[1][b][c][d][e] + casecounts5[2][b][c][d][e];
											}
											else
											{
												s1 = controlcounts5[0][b][c][d][e] + controlcounts5[1][b][c][d][e] + controlcounts5[2][b][c][d][e];
											}
											s2 = mAlt[0][b][c][d][e][k] + mAlt[1][b][c][d][e][k] + mAlt[2][b][c][d][e][k];

										}
									}
									else if (dim == 3)
									{
										if (l % 4 == 1)
										{
											s1 = casecounts5[a][b][c][d][e] + controlcounts5[a][b][c][d][e];
											s2 = mAlt[a][b][c][d][e][0] + mAlt[a][b][c][d][e][1];
										}
										else if (l % 4 == 2)
										{
											if (k == 0)
											{
												s1 = casecounts5[a][b][0][d][e]+ casecounts5[a][b][1][d][e]+ casecounts5[a][b][2][d][e];
											}
											else
											{
												s1 = controlcounts5[a][b][0][d][e]+ controlcounts5[a][b][1][d][e]+ controlcounts5[a][b][2][d][e];
											}
											s2 = mAlt[a][b][0][d][e][k] + mAlt[a][b][1][d][e][k] + mAlt[a][b][2][d][e][k];
										}
										else if (l % 4 == 3)
										{
											if (k == 0)
											{
												s1 = casecounts5[a][0][c][d][e] + casecounts5[a][1][c][d][e] + casecounts5[a][2][c][d][e];
											}
											else
											{
												s1 = controlcounts5[a][0][c][d][e] + controlcounts5[a][1][c][d][e] + controlcounts5[a][2][c][d][e];
											}
											s2 = mAlt[a][0][c][d][e][k] + mAlt[a][1][c][d][e][k] + mAlt[a][2][c][d][e][k];

										}
										else if (l % 4 == 0)
										{
											if (k == 0)
											{
												s1 = casecounts5[0][b][c][d][e] + casecounts5[1][b][c][d][e] + casecounts5[2][b][c][d][e];
											}
											else
											{
												s1 = controlcounts5[0][b][c][d][e] + controlcounts5[1][b][c][d][e] + controlcounts5[2][b][c][d][e];
											}
											s2 = mAlt[0][b][c][d][e][k] + mAlt[1][b][c][d][e][k] + mAlt[2][b][c][d][e][k];

										}
									}
									if (s2 > 0)
									{
										mNeu[a][b][c][d][e][k] = mAlt[a][b][c][d][e][k] * s1 / s2;
									}
									else
									{
										mNeu[a][b][c][d][e][k] = 0;
									}
								} //k
							}
						}
					}
				} //b
			} //a end update mNeu


			//update mAlt and likelihood

			for (a = 0; a < 3; a++)
			{
				for (b = 0; b < 3; b++)
				{
					for (c = 0; c < stopC; c++)
					{
						for (d = 0; d < stopD; d++)
						{
							for (e = 0; e < stopE; e++)
							{
								for (k = 0; k < 2; k++)
								{
									mAlt[a][b][c][d][e][k] = mNeu[a][b][c][d][e][k];

									if (k == 0)
									{
										if (mNeu[a][b][c][d][e][k] > 0 && casecounts5[a][b][c][d][e] > 0)
										{
											statNeu+= -2 * (casecounts5[a][b][c][d][e]* log(mNeu[a][b][c][d][e][k]) - casecounts5[a][b][c][d][e]* log(casecounts5[a][b][c][d][e]));
										}
									} //k==0
									else
									{
										if (mNeu[a][b][c][d][e][k] > 0 && controlcounts5[a][b][c][d][e] > 0)
										{
											statNeu += -2 * (controlcounts5[a][b][c][d][e]* log(mNeu[a][b][c][d][e][k])- controlcounts5[a][b][c][d][e]* log(controlcounts5[a][b][c][d][e]));

										}
									}
								}//k
							}
						}
					}
				}//b
			} //a end update mAlt

			if (fabs(statAlt - statNeu) < 0.0001)
			{
				break;
			}
			statAlt = statNeu;
		} //end l
	} //end if df > 0

	//cout << "df " << df << "stat " << statNeu << "\n";

	if (statNeu > 0 && df > 0)
	{
		tstat.p = pValueCalc(df / 2, statNeu / 2);
	}
	else
	{
		tstat.p = 1;
	}

	tstat.pmod = tstat.p;

	return tstat;

} //end chitestinter5


//NEU
struct TSTAT chiTestInterX5(double casecounts5M[3][3][3][3][3],double casecounts5F[3][3][3][3][3],
                            double controlcounts5M[3][3][3][3][3],double controlcounts5F[3][3][3][3][3],
                            int dim, int mWithSingletop,double helpstat)
{
	double mAlt[3][3][3][3][3][2];
	double mNeu[3][3][3][3][3][2];
	double s1 = 0;
	double s2 = 0;
	double statAlt = 0;
	double statNeu = 0;
	double statNeuF = 0;
	int a, b, c, d, e, k, l;
	double df = 8;
	double df_M = 8;
	double df_F = 8;
	int stopC = 3;
	int stopD = 3;
	int stopE = 3;
	int merkab = 0;
	int merkac = 0;
	int merkbc = 0;
	int merka = 0;
	int merkb = 0;
	int merkc = 0;

	struct TSTAT tstat;
	struct TSTAT tstatF;
	tstatF.p = 0;
	tstatF.pmod = 0;
	tstat.p = 0;
	tstat.pmod = 0;


	if (dim == 2)
	{
		stopC = 1;
		stopD = 1;
		stopE = 1;
	}
	else if (dim == 3)
	{
		stopD = 1;
		stopE = 1;
	}
	else if (dim == 4)
	{
		stopE = 1;
	}

	//MALE
	//update df_M
	for (a = 0; a < 3; a++)
	{
		for (b = 0; b < 3; b++)
		{
			if (casecounts5M[a][b][0][0][0] == 0
				&& casecounts5M[a][b][1][0][0] == 0
				&& casecounts5M[a][b][2][0][0] == 0
				&& controlcounts5M[a][b][0][0][0] == 0
				&& controlcounts5M[a][b][1][0][0] == 0
				&& controlcounts5M[a][b][2][0][0] == 0)
			{
				merkab++;
			}
		}
	}

	for (a = 0; a < 3; a++)
	{
		for (c = 0; c < stopC; c++)
		{
			if (casecounts5M[a][0][c][0][0] == 0
				&& casecounts5M[a][1][c][0][0] == 0
				&& casecounts5M[a][2][c][0][0] == 0
				&& controlcounts5M[a][0][c][0][0] == 0
				&& controlcounts5M[a][1][c][0][0] == 0
				&& controlcounts5M[a][2][c][0][0] == 0)
			{
				merkac++;
			}
		}
	}

	for (b = 0; b < 3; b++)
	{
		for (c = 0; c < stopC; c++)
		{
			if (casecounts5M[0][b][c][0][0] == 0
			&& casecounts5M[1][b][c][0][0] == 0
			&& casecounts5M[2][b][c][0][0] == 0
			&& controlcounts5M[0][b][c][0][0] == 0
			&& controlcounts5M[1][b][c][0][0] == 0
			&& controlcounts5M[2][b][c][0][0] == 0)
			{
				merkbc++;
			}
		}
	}

	merka = merkab + merkac;
	merkb = merkab + merkbc;
	merkc = merkac + merkbc;

	if (dim == 3)
	{
		if (merka == 2 || merkb == 2 || merkc == 2)
		{
		  df_M = 0;
		}
		else if ((merka + merkb + merkc) == 3)
		{
		  df_M = 1;
		}
		else if ((merka + merkb + merkc) == 2)
		{
		  df_M = 2;
		}
		else if ((merka + merkb + merkc) == 1)
		{
		  df_M = 4;
		}
	}
	else if (dim == 2)
	{
		if (df - merkab == 8)
		{
			df_M = 4;
		}
		else if (df - merkab == 7)
		{
			df_M = 3;
		}
		else if  (df - merkab == 6)
		{
			df_M = 3;
		}
		else if  (df - merkab == 5)
		{
			df_M = 2;
		}
		else if  (df - merkab == 4)
		{
			df_M = 2;
		}
		else if  (df - merkab == 3)
		{
			df_M = 1;
		}
		else
		{
			df_M = 0;
		}
	}


	if (df_M > 0)
	{
		// Startvalues
		for (a = 0; a < 3; a++)
		{
			for (b = 0; b < 3; b++)
			{
				for (c = 0; c < stopC; c++)
				{
					for (d = 0; d < stopD; d++)
					{
						for (e = 0; e < stopE; e++)
						{
							for (k = 0; k < 2; k++)
							{
								mAlt[a][b][c][d][e][k] = 1;
							}//k
						}
					}
				}
			}//b
		} //end a startvalues

		for (l = 1; l < 1000; l++)
		{
			statNeu = 0;

			//Compute mNeu
			for (a = 0; a < 3; a++)
			{
				for (b = 0; b < 3; b++)
				{
					for (c = 0; c < stopC; c++)
					{
						for (d = 0; d < stopD; d++)
						{
							for (e = 0; e < stopE; e++)
							{
								for (k = 0; k < 2; k++)
								{
									if (dim == 2)
									{
										if (l % 3 == 1)
										{
											s1 = casecounts5M[a][b][c][d][e] + controlcounts5M[a][b][c][d][e];
											s2 = mAlt[a][b][c][d][e][0] + mAlt[a][b][c][d][e][1];
										}
										else if (l % 3 == 2)
										{
											if (k == 0)
											{
												s1 = casecounts5M[a][0][c][d][e] + casecounts5M[a][1][c][d][e] + casecounts5M[a][2][c][d][e];
											}
											else
											{
												s1 = controlcounts5M[a][0][c][d][e] + controlcounts5M[a][1][c][d][e] + 														controlcounts5M[a][2][c][d][e];
											}
											s2 = mAlt[a][0][c][d][e][k] + mAlt[a][1][c][d][e][k] + mAlt[a][2][c][d][e][k];

										}
										else if (l % 3 == 0)
										{

											if (k == 0)
											{
												s1 = casecounts5M[0][b][c][d][e] + casecounts5M[1][b][c][d][e]
													+ casecounts5M[2][b][c][d][e];
											}
											else
											{
												s1 = controlcounts5M[0][b][c][d][e] + controlcounts5M[1][b][c][d][e] +
													controlcounts5M[2][b][c][d][e];
											}
											s2 = mAlt[0][b][c][d][e][k] + mAlt[1][b][c][d][e][k] + mAlt[2][b][c][d][e][k];

										}
									}
									if (dim == 3)
									{
										if (l % 4 == 1)
										{
											s1 = casecounts5M[a][b][c][d][e]+ controlcounts5M[a][b][c][d][e];
											s2 = mAlt[a][b][c][d][e][0]+ mAlt[a][b][c][d][e][1];
										}
										else if (l % 4 == 2)
										{
											if (k == 0)
											{
												s1 = casecounts5M[a][b][0][d][e]+ casecounts5M[a][b][1][d][e]+
													casecounts5M[a][b][2][d][e];
											}
											else
											{
												s1 = controlcounts5M[a][b][0][d][e]+ controlcounts5M[a][b][1][d][e]+
													controlcounts5M[a][b][2][d][e];
											}
											s2 = mAlt[a][b][0][d][e][k]+ mAlt[a][b][1][d][e][k]+ mAlt[a][b][2][d][e][k];
										}
										else if (l % 4 == 3)
										{
											if (k == 0)
											{
												s1 = casecounts5M[a][0][c][d][e]+ casecounts5M[a][1][c][d][e]+
													casecounts5M[a][2][c][d][e];
											}
											else
											{
												s1 = controlcounts5M[a][0][c][d][e]+ controlcounts5M[a][1][c][d][e]+
													controlcounts5M[a][2][c][d][e];
											}
											s2 = mAlt[a][0][c][d][e][k]+ mAlt[a][1][c][d][e][k]+ mAlt[a][2][c][d][e][k];
										}
										else if (l % 4 == 0)
										{
											if (k == 0)
											{
												s1 = casecounts5M[0][b][c][d][e]+ casecounts5M[1][b][c][d][e]+
													casecounts5M[2][b][c][d][e];
											}
											else
											{
												s1 = controlcounts5M[0][b][c][d][e]+ controlcounts5M[1][b][c][d][e]+
													controlcounts5M[2][b][c][d][e];
											}
											s2 = mAlt[0][b][c][d][e][k]+ mAlt[1][b][c][d][e][k]+ mAlt[2][b][c][d][e][k];
										}
									}
									if (s2 > 0)
									{
										mNeu[a][b][c][d][e][k] = mAlt[a][b][c][d][e][k] * s1 / s2;
									}
									else
									{
										mNeu[a][b][c][d][e][k] = 0;
									}
								}//k
							}
						}
					}
				}//b
			} //a end update mNeu


			//update mAlt and likelihood

			for (a = 0; a < 3; a++)
			{
				for (b = 0; b < 3; b++)
				{
					for (c = 0; c < stopC; c++)
					{
						for (d = 0; d < stopD; d++)
						{
							for (e = 0; e < stopE; e++)
							{
								for (k = 0; k < 2; k++)
								{
									mAlt[a][b][c][d][e][k] = mNeu[a][b][c][d][e][k];


									if (k == 0)
									{
										if (mNeu[a][b][c][d][e][k] > 0 && casecounts5M[a][b][c][d][e] > 0)
										{
											statNeu += -2 * (casecounts5M[a][b][c][d][e]* log(mNeu[a][b][c][d][e][k]) - casecounts5M[a][b][c][d][e]* log(casecounts5M[a][b][c][d][e]));
										}
									} //k==0
									else
									{
										if (mNeu[a][b][c][d][e][k] > 0 && controlcounts5M[a][b][c][d][e] > 0)
										{
											statNeu += -2 * (controlcounts5M[a][b][c][d][e]* log(mNeu[a][b][c][d][e][k])- controlcounts5M[a][b][c][d][e]* log(controlcounts5M[a][b][c][d][e]));
										}
									}
								}//k
							}
						}
					}
				}//b
			} //a end update mAlt

			if (fabs(statAlt - statNeu) < 0.0001)
			{
				break;
			}
			statAlt = statNeu;
		} //end l
	} //end if df_M > 0


	//FEMALE
	statNeuF = 0;
	statAlt = 0;

	//update df_F
	merkab = 0;
	merkac = 0;
	merkbc = 0;
	merka = 0;
	merkb = 0;
	merkc = 0;

	for (a = 0; a < 3; a++)
	{
		for (b = 0; b < 3; b++)
		{
			if (casecounts5F[a][b][0][0][0] == 0
				&& casecounts5F[a][b][1][0][0] == 0
				&& casecounts5F[a][b][2][0][0] == 0
				&& controlcounts5F[a][b][0][0][0] == 0
				&& controlcounts5F[a][b][1][0][0] == 0
				&& controlcounts5F[a][b][2][0][0] == 0)
			{
				merkab++;
			}
		}
	}

	for (a = 0; a < 3; a++)
	{
		for (c = 0; c < stopC; c++)
		{
			if (casecounts5F[a][0][c][0][0] == 0
			 	&& casecounts5F[a][1][c][0][0]== 0
				&& casecounts5F[a][2][c][0][0] == 0
				&& controlcounts5F[a][0][c][0][0] == 0
				&& controlcounts5F[a][1][c][0][0] == 0
				&& controlcounts5F[a][2][c][0][0] == 0)
			{
					merkac++;
			}
		}
	}

	for (b = 0; b < 3; b++)
	{
		for (c = 0; c < stopC; c++)
		{
			if (casecounts5F[0][b][c][0][0] == 0
			 	&& casecounts5F[1][b][c][0][0] == 0
				&& casecounts5F[2][b][c][0][0] == 0
				&& controlcounts5F[0][b][c][0][0] == 0
				&& controlcounts5F[1][b][c][0][0] == 0
				&& controlcounts5F[2][b][c][0][0] == 0)
			{
				merkbc++;
			}
		}
	}

	merka = merkab + merkac;
	merkb = merkab + merkbc;
	merkc = merkac + merkbc;


	if (dim == 3)
	{
		if (merka == 2 || merkb == 2 || merkc == 2)
		{
			df_F = 0;
		}
		else if ((merka + merkb + merkc) == 3)
		{
			df_F = 1;
		}
		else if ((merka + merkb + merkc) == 2)
		{
			df_F = 2;
		}
		else if ((merka + merkb + merkc) == 1)
		{
			df_F = 4;
		}
	}
	else if (dim == 2)
	{
		if (df - merkab == 8)
		{
				df_F = 4;
		}
		else if (df - merkab == 7)
		{
				df_F = 3;
		}
		else if (df - merkab == 6)
		{
				df_F = 2;
		}
		else if (df - merkab == 5)
		{
				df_F = 2;
		}
		else if (df - merkab == 4)
		{
				df_F = 1;
		}
		else
		{
				df_F = 0;
		}
	}

	if (df_F > 0)
	{
		// Startvalues
		for (a = 0; a < 3; a++)
		{
			for (b = 0; b < 3; b++)
			{
				for (c = 0; c < stopC; c++)
				{
					for (d = 0; d < stopD; d++)
					{
						for (e = 0; e < stopE; e++)
						{
							for (k = 0; k < 2; k++)
							{
								mAlt[a][b][c][d][e][k] = 1;
							}//k
						}
					}
				}
			}//	b
		} //end a startvalues

		for (l = 1; l < 1000; l++)
		{
	 		statNeuF = 0;

			//Compute mNeu
			for (a = 0; a < 3; a++)
			{
				for (b = 0; b < 3; b++)
				{
					for (c = 0; c < stopC; c++)
					{
						for (d = 0; d < stopD; d++)
						{
							for (e = 0; e < stopE; e++)
							{
								for (k = 0; k < 2; k++)
								{
									if (dim == 2)
									{
										if (l % 3 == 1)
										{
											s1 = casecounts5F[a][b][c][d][e] + controlcounts5F[a][b][c][d][e];
											s2 = mAlt[a][b][c][d][e][0] + mAlt[a][b][c][d][e][1];
										}
										else if (l % 3 == 2)
										{
											if (k == 0)
											{
												s1 = casecounts5F[a][0][c][d][e] + casecounts5F[a][1][c][d][e] + casecounts5F[a][2][c][d][e];
											}
											else
											{
												s1 = controlcounts5F[a][0][c][d][e] + controlcounts5F[a][1][c][d][e] + controlcounts5F[a][2][c][d][e];
											}
											s2 = mAlt[a][0][c][d][e][k] + mAlt[a][1][c][d][e][k] + mAlt[a][2][c][d][e][k];

										}
										else if (l % 3 == 0)
										{
											if (k == 0)
											{
												s1 = casecounts5F[0][b][c][d][e] + casecounts5F[1][b][c][d][e] + casecounts5F[2][b][c][d][e];
											}
											else
											{
												s1 = controlcounts5F[0][b][c][d][e] + controlcounts5F[1][b][c][d][e] + controlcounts5F[2][b][c][d][e];
											}
											s2 = mAlt[0][b][c][d][e][k] + mAlt[1][b][c][d][e][k] + mAlt[2][b][c][d][e][k];
										}
									}
									if (dim == 3)
									{
										if (l % 4 == 1)
										{
											s1 = casecounts5F[a][b][c][d][e] + controlcounts5F[a][b][c][d][e];
											s2 = mAlt[a][b][c][d][e][0] + mAlt[a][b][c][d][e][1];
										}

										else if (l % 4 == 2)
										{
											if (k == 0)
											{
												s1 = casecounts5F[a][b][0][d][e] + casecounts5F[a][b][1][d][e] + casecounts5F[a][b][2][d][e];
											}
											else
											{
												s1 = controlcounts5F[a][b][0][d][e]+ controlcounts5F[a][b][1][d][e] + controlcounts5F[a][b][2][d][e];
											}
											s2 = mAlt[a][b][0][d][e][k]+ mAlt[a][b][1][d][e][k]+ mAlt[a][b][2][d][e][k];
										}
										else if (l % 4 == 3)
										{
											if (k == 0)
											{
												s1 = casecounts5F[a][0][c][d][e]+ casecounts5F[a][1][c][d][e] + casecounts5F[a][2][c][d][e];
											}
											else
											{
											s1 = controlcounts5F[a][0][c][d][e] + controlcounts5F[a][1][c][d][e] + controlcounts5F[a][2][c][d][e];
											}
											s2 = mAlt[a][0][c][d][e][k] + mAlt[a][1][c][d][e][k]+ mAlt[a][2][c][d][e][k];
										}
										else if (l % 4 == 0)
										{
											if (k == 0)
											{
												s1 = casecounts5F[0][b][c][d][e] + casecounts5F[1][b][c][d][e] + casecounts5F[2][b][c][d][e];
											}
											else
											{
												s1 = controlcounts5F[0][b][c][d][e] + controlcounts5F[1][b][c][d][e] + controlcounts5F[2][b][c][d][e];
											}
										s2 = mAlt[0][b][c][d][e][k] + mAlt[1][b][c][d][e][k]+ mAlt[2][b][c][d][e][k];
										}
									}
									if (s2 > 0)
									{
										mNeu[a][b][c][d][e][k] = mAlt[a][b][c][d][e][k] * s1/ s2;
									}
									else
									{
										mNeu[a][b][c][d][e][k] = 0;
									}
								}//k
							}
						}
					}
				}//b
			} //a end update mNeu


			//update mAlt and likelihood

			for (a = 0; a < 3; a++)
			{
				for (b = 0; b < 3; b++)
				{
					for (c = 0; c < stopC; c++)
					{
						for (d = 0; d < stopD; d++)
						{
							for (e = 0; e < stopE; e++)
							{
								for (k = 0; k < 2; k++)
								{
									mAlt[a][b][c][d][e][k] = mNeu[a][b][c][d][e][k];

									if (k == 0)
									{
										if (mNeu[a][b][c][d][e][k] > 0 && casecounts5F[a][b][c][d][e]> 0)
										{
											statNeuF += -2* (casecounts5F[a][b][c][d][e]* log(mNeu[a][b][c][d][e][k]) - casecounts5F[a][b][c][d][e]* log(casecounts5F[a][b][c][d][e]));
										}
									} //k==0
									else
									{
										if (mNeu[a][b][c][d][e][k] > 0 && controlcounts5F[a][b][c][d][e]> 0)
										{
											statNeuF += -2 * (controlcounts5F[a][b][c][d][e]* log(mNeu[a][b][c][d][e][k]) - controlcounts5F[a][b][c][d][e]* log(controlcounts5F[a][b][c][d][e]));
										}
									}
								}//k
							}
						}
					}
				}//b
			} //a end update mAlt

			if (fabs(statAlt - statNeuF) < 0.0001)
			{
				break;
			}
			statAlt = statNeuF;
		} //end l
	} //end if df_F > 0


	//JOINT:

	df = df_M + df_F;
	statNeu = statNeu + statNeuF;

	if (statNeu > 0 && df > 0)
	{
		tstat.p = pValueCalc(df / 2, statNeu / 2);
	}
	else
	{
		tstat.p = 1;
	}

	tstat.pmod = tstat.p;

	return tstat;
}


struct TSTAT chiTestX5(double casecounts5M[3][3][3][3][3],double casecounts5F[3][3][3][3][3],
                  		double controlcounts5M[3][3][3][3][3],double controlcounts5F[3][3][3][3][3],

                  		int dim, int mWithSingletop,double helpstat)
{
	double expectedCasesM[3][3][3][3][3];
	double expectedCasesF[3][3][3][3][3];
	double expectedControlsM[3][3][3][3][3];
	double expectedControlsF[3][3][3][3][3];
	double sumCasesM = 0;
	double sumCasesF = 0;
	double sumControlsM = 0;
	double sumControlsF = 0;
	double stat = 0;
	double statNew = 0;
	double statM = 0;
	double statF = 0;
	double df = 0;
	double df_M = -1;
	double df_F = -1;
	int a, b, c, d, e;

	int stopC = 3;
	int stopD = 3;
	int stopE = 3;

	struct TSTAT tstat;
	tstat.p = 0;
	tstat.pmod = 0;

	if (dim == 2)
	{
		stopC = 1;
		stopD = 1;
		stopE = 1;
	}
	else if (dim == 3)
	{
		stopD = 1;
		stopE = 1;
	}
	else if (dim == 4)
	{
		stopE = 1;
	}

	// Summe Fälle und Kontrollen
	for (a = 0; a < 3; a++)
	{
		for (b = 0; b < 3; b++)
		{
			for (c = 0; c < stopC; c++)
			{
				for (d = 0; d < stopD; d++)
				{
					for (e = 0; e < stopE; e++)
					{
						// Fälle
						sumCasesM = sumCasesM + casecounts5M[a][b][c][d][e];
						sumCasesF = sumCasesF + casecounts5F[a][b][c][d][e];

						expectedCasesM[a][b][c][d][e] = 0;
						expectedCasesF[a][b][c][d][e] = 0;

						// Kontrolls
						sumControlsM = sumControlsM + controlcounts5M[a][b][c][d][e];
						sumControlsF = sumControlsF + controlcounts5F[a][b][c][d][e];

						expectedControlsM[a][b][c][d][e] = 0;
						expectedControlsF[a][b][c][d][e] = 0;
					}
				}
			}
		}
	}

	// Fälle und Kontrollen erwartet
	for (a = 0; a < 3; a++)
	{
		for (b = 0; b < 3; b++)
		{
			for (c = 0; c < stopC; c++)
			{
				for (d = 0; d < stopD; d++)
				{
					for (e = 0; e < stopE; e++)
					{
						expectedCasesM[a][b][c][d][e] = (sumCasesM * (casecounts5M[a][b][c][d][e] + controlcounts5M[a][b][c][d][e]))/ (sumCasesM + sumControlsM);
						expectedCasesF[a][b][c][d][e] = (sumCasesF * (casecounts5F[a][b][c][d][e] + controlcounts5F[a][b][c][d][e]))/(sumCasesF + sumControlsF);

						expectedControlsM[a][b][c][d][e] = (sumControlsM * (casecounts5M[a][b][c][d][e] + controlcounts5M[a][b][c][d][e]))/ (sumCasesM + sumControlsM);
						expectedControlsF[a][b][c][d][e] = (sumControlsF * (casecounts5F[a][b][c][d][e] + controlcounts5F[a][b][c][d][e]))/ (sumCasesF + sumControlsF);
					}
				}
			}
		}
	}

	for (a = 0; a < 3; a++)
	{
		for (b = 0; b < 3; b++)
		{
			for (c = 0; c < stopC; c++)
			{
				for (d = 0; d < stopD; d++)
				{
					for (e = 0; e < stopE; e++)
					{
						if (expectedCasesM[a][b][c][d][e] > 0 || expectedControlsM[a][b][c][d][e] > 0)
						{
							if (expectedCasesM[a][b][c][d][e] > 0)
							{
								statM = statM + (((casecounts5M[a][b][c][d][e]- expectedCasesM[a][b][c][d][e]) * (casecounts5M[a][b][c][d][e]- expectedCasesM[a][b][c][d][e]))/(expectedCasesM[a][b][c][d][e]));
								df_M++;
							}
							if (expectedControlsM[a][b][c][d][e] > 0)
							{
								statM = statM + (((controlcounts5M[a][b][c][d][e] - expectedControlsM[a][b][c][d][e]) * (controlcounts5M[a][b][c][d][e]- expectedControlsM[a][b][c][d][e]))/(expectedControlsM[a][b][c][d][e]));
							}
						}
						else
						{
							statM = statM + 0;
						}

						if (expectedCasesF[a][b][c][d][e] > 0 || expectedControlsF[a][b][c][d][e] > 0)
						{
							if (expectedCasesF[a][b][c][d][e] > 0)
							{
								statF = statF + (((casecounts5F[a][b][c][d][e] - expectedCasesF[a][b][c][d][e]) * (casecounts5F[a][b][c][d][e] - expectedCasesF[a][b][c][d][e]))/(expectedCasesF[a][b][c][d][e]));
								df_F++;
							}
							if (expectedControlsF[a][b][c][d][e] > 0)
							{
								statF = statF + (((controlcounts5F[a][b][c][d][e] - expectedControlsF[a][b][c][d][e]) * (controlcounts5F[a][b][c][d][e]- expectedControlsF[a][b][c][d][e]))/(expectedControlsF[a][b][c][d][e]));
							}
						}
						else
						{
							statF = statF + 0;
						}
					}
				}
			}
		}
	}

	// Add teststat and d.f.s
	stat = statM + statF;

	df = 0;
	if (df_M >= 1)
	{
		df = df + df_M;
	}

	if (df_F >= 1)
	{
		df = df + df_F;
	}

	if (stat != 0 && df > 0)
	{
		tstat.p = pValueCalc((df) / 2, stat / 2);
	}
	else
	{
		tstat.p = 1;
	}

	if (mWithSingletop == 0)
	{
		tstat.pmod = tstat.p;
	}
	else
	{
		if (df_M >= 1 && statM > mWithSingletop * helpstat)
		{
			statNew = stat - (mWithSingletop * helpstat); //conservative
		}
		else
		{
			statNew = stat;
		}

		if (df_F >= 1 && statF > mWithSingletop * helpstat)
		{
			statNew = statNew - (mWithSingletop * helpstat);
		}

		if (statNew < 0)
		{
			statNew = 0;
		}
		tstat.pmod = pValueCalc((df) / 2, statNew / 2);
	}
	return tstat;
}


void insertionChi5(struct BESTCHI5 *bestchi5, int lastPos)
{
	int j;
	double currentValue;
	struct BESTCHI5 currentStructValue;

	// store current position
	currentStructValue = bestchi5[lastPos];
	currentValue = bestchi5[lastPos].p;

	// Find smaller value than currentValue. Schleife
	// Loop goes from the current position from right to left
	for (j = lastPos; j > 0 && bestchi5[j - 1].p > currentValue; j--)
	{
		//Wenn Vorgänger größer als aktuelles Element in currentValue
		bestchi5[j] = bestchi5[j - 1];
	}
	// gespeichertes Element an neue Position -> Lücke auffüllen
	bestchi5[j] = currentStructValue;
}


void insertionTop(struct TOPLIST *toplist, int lastPos)
{
	int j;
	double currentValue;
	struct TOPLIST currentStructValue;


	// aktuelle Position zwischenspeichern
	currentStructValue = toplist[lastPos];
	currentValue = toplist[lastPos].p;

	// Kleineren Wert als currentValue suchen. Schleife
	// durchläuft von aktueller Pos. von rechts nach links
	for (j = lastPos; j > 0 && toplist[j - 1].p > currentValue; j--)
	{
		//Wenn Vorgänger größer als aktuelles Element in currentValue
		toplist[j] = toplist[j - 1];
	}
	// gespeichertes Element an neue Position -> Lücke auffüllen
	toplist[j] = currentStructValue;
}


// Funktion um Zufallszahlen zu generieren
double HillRandom(int *ix, int *iy, int *iz)
{
	double random = 0;
	double v;

	*ix = 171 * ((*ix) % 177) - 2 * (*ix / 177);
	*iy = 172 * ((*iy) % 176) - 35 * (*iy / 176);
	*iz = 170 * ((*iz) % 178) - 63 * (*iz / 178);

	if (*ix < 0)
	{
		*ix = *ix + 30269;
	}
	if (*iy < 0)
	{
		*iy = *iy + 30307;
	}
	if (*iz < 0)
	{
		*iz = *iz + 30323;
	}

	random = modf(((float) (*ix)) / 30269.0 + ((float) (*iy)) / 30307.0 + ((float) (*iz)) / 30323.0, &v);

	return random;
}


void qt_permute(int nlinestfam, double **YY, int *ix, int *iy, int *iz, struct PERSON *person, int ncasesqc,double *YY1,double *YY2)
{
	double random1;
	int i,k,l;
	int lines;

	lines=nlinestfam;

	for(i =0; i<nlinestfam; i++)
	{
		if (person[i].qcin == 1)
		{
			for(k = 0; k >= 0; k++)
			{
#pragma omp critical(HW_RANDOM)
			  random1 = HillRandom(ix, iy, iz);
				l = (int) (random1 * lines);
				if (YY2[l] == 1)
				{
					(*YY)[i] = YY1[l];
					YY1[l]=YY1[lines-1];
					YY2[l]=YY2[lines-1];

					lines--;
					break;
				}
			}//k
		} //if
	} //i
}


void qt_permute_all(struct PERSON* person, int nlinestfam, int thread, int *ix, int *iy, int *iz) {
    double random1;
    double buf_qtaff;
	double* buf_cov;
	unsigned char* buf_covin;
	for (int l,i=nlinestfam-1; i>0; i--) {
		if (person[i].qcin) {
			do {
                #pragma omp critical(HW_RANDOM)
                random1 = HillRandom(ix, iy, iz);
                l = (int)((i+1) * random1);
            } while (!person[l].qcin);
            buf_qtaff               = person[i].qtaff[thread];
            person[i].qtaff[thread] = person[l].qtaff[thread];
            person[l].qtaff[thread] = buf_qtaff;
		}
	}
}


void qt_permute_clusters(struct PERSON *person, struct CLUSTERS* Clus, uint32_t n, int thread, int *ix, int *iy, int *iz) {
    double random1;
    double buf_qtaff;
    for (uint32_t c=0; c<n; c++) {
        for (int l,m,i=Clus[c].nppls-1; i>0; i--) {
            #pragma omp critical(HW_RANDOM)
            random1 = HillRandom(ix, iy, iz);
            l = Clus[c].list[(int)((i+1) * random1)];
            m = Clus[c].list[i];
            buf_qtaff               = person[m].qtaff[thread];
            person[m].qtaff[thread] = person[l].qtaff[thread];
            person[l].qtaff[thread] = buf_qtaff;
		}
	}
}


void aff_permute_all(struct PERSON *person, int nlinestfam, int ncasesqc, int thread, int *ix, int *iy, int *iz) {
    double random1;
  	for (int k=0; k<nlinestfam; k++)
   	    if (person[k].qcin)
            person[k].aff[thread] = 1;
    int i, newcases = 0;
    while (newcases < ncasesqc) {
        #pragma omp critical(HW_RANDOM)
		random1 = HillRandom(ix, iy, iz);
        i = (int)(nlinestfam * random1);
        if (person[i].aff[thread] == 1) {
            person[i].aff[thread] = 2;
            newcases++;
        }
    }
}


void aff_permute_pairs(struct PERSON* person, struct RELATIVES* Pairs, uint32_t n, int thread) {
    uint32_t buf;
    for (uint32_t i=0; i<n; i++) {
        if (rand() & 1) {
            buf                            = person[Pairs[i].i].aff[thread];
            person[Pairs[i].i].aff[thread] = person[Pairs[i].j].aff[thread];
            person[Pairs[i].j].aff[thread] = buf;
        }
    }
}


void aff_permute_clusters(struct PERSON* person, struct CLUSTERS* Clus, uint32_t n, int thread, int *ix, int *iy, int *iz) {
    double random1;
    for (uint32_t c=0; c<n; c++) {
        for (uint32_t k=0; k<Clus[c].nppls; k++)
            person[Clus[c].list[k]].aff[thread] = 1;
        uint32_t i, newcases = 0;
        while (newcases < Clus[c].ncases) {
            #pragma omp critical(HW_RANDOM)
            random1 = HillRandom(ix, iy, iz);
            i = Clus[c].list[(int)(Clus[c].nppls * random1)];
            if (person[i].aff[thread]==1) {
                person[i].aff[thread] = 2;
                newcases++;
            }
        }
    }
}


int compareINFOMAP(struct INFOMAP info1,struct INFOMAP info2)
{

	if(strcmp(info1.rs,info2.rs)>0)
	{
		return 2;
	}
	else if(strcmp(info1.rs,info2.rs)<0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
};


void qsortbestsinglemarkerNR(struct BESTSINGLEMARKER **bestsinglemarker2,int left, int right, struct BESTSINGLEMARKER *x, struct BESTSINGLEMARKER *y)
{
	//sorts in ascending order
	int i,j;

	i=left;j=right;

	(*x)=(*bestsinglemarker2)[(left+right)/2];

	do
	{
		while( ((*bestsinglemarker2)[i].nr < (*x).nr) && (i <right)) {i++;}
		while( ((*x).nr < (*bestsinglemarker2)[j].nr) && (j>left)) {j--;}

		if (i<=j)
		{
			(*y)=(*bestsinglemarker2)[i];
			(*bestsinglemarker2)[i]=(*bestsinglemarker2)[j];
			(*bestsinglemarker2)[j]= (*y);
			i++;j--;
		}
	}while (i<=j);

	if (left <j)
	{
		qsortbestsinglemarkerNR(bestsinglemarker2,left,j, x, y);
	}
	if (i < right)
	{
		qsortbestsinglemarkerNR(bestsinglemarker2,i,right, x, y);
	}
};


void qsortbestsinglemarker(struct BESTSINGLEMARKER **bestsinglemarker,int left, int right, struct BESTSINGLEMARKER *x, struct BESTSINGLEMARKER *y)
{
	//sorts in ascending order
	int i,j;


	i=left;j=right;

	(*x)=(*bestsinglemarker)[(left+right)/2];

	do
	{
		while( ((*bestsinglemarker)[i].p < (*x).p) && (i <right)) {i++;}
		while( ((*x).p < (*bestsinglemarker)[j].p) && (j>left)) {j--;}

		if (i<=j)
		{
			(*y)=(*bestsinglemarker)[i];
			(*bestsinglemarker)[i]=(*bestsinglemarker)[j];
			(*bestsinglemarker)[j]= *y;
			i++;j--;
		}
	}while (i<=j);

	if (left <j)
	{
		qsortbestsinglemarker(bestsinglemarker,left,j, x, y);
	}
	if (i < right)
	{
		qsortbestsinglemarker(bestsinglemarker,i,right, x, y);
	}
};


void qsortbestchi(struct BESTCHI5 **bestchi,int left, int right, struct BESTCHI5 *x, struct BESTCHI5 *y)
{
	//sorts in ascending order
	int i,j;

	i=left;j=right;

	*x=(*bestchi)[(left+right)/2];

		do
		{
			while( ((*bestchi)[i].p < (*x).p) && (i <right)) {i++;}
			while( ((*x).p < (*bestchi)[j].p) && (j>left)) {j--;}

			if (i<=j)
			{
			 (*y)=(*bestchi)[i];
			 (*bestchi)[i]=(*bestchi)[j];
			 (*bestchi)[j]= (*y);
			  i++;j--;
			}
		}while (i<=j);

	if(i< left || j > right) {cout << "warning qsortbestchi!\n";}

	if (left <j)
	{
		qsortbestchi(bestchi,left,j,x,y);
	}
	if (i < right)
	{
		qsortbestchi(bestchi,i,right,x,y);
	}
	//exit(1);

};


void qsortgeneticlist(struct GENETICLIST **geneticlist2,int left, int right, struct GENETICLIST *x, struct GENETICLIST *y)
{
	//sorts in ascending order
	int i,j;

	i=left;j=right;

	(*x)=(*geneticlist2)[(left+right)/2];


	do
	{
		while( (*geneticlist2)[i].nr < (*x).nr && (i <right)) {i++;}
		while( (*x).nr < (*geneticlist2)[j].nr && (j>left)) {j--;}

		if (i<=j)
		{
			(*y)=(*geneticlist2)[i];
			(*geneticlist2)[i]=(*geneticlist2)[j];
			(*geneticlist2)[j]=(*y);
			i++;j--;
		}
	}while (i<=j);

   if (left <j)
   {
		qsortgeneticlist(geneticlist2,left,j, x, y);
   }
   if (i < right)
   {
		qsortgeneticlist(geneticlist2,i,right, x, y);
   }

};


template <typename T>
void qsortrelatives(struct RELATIVES*& relatives, int left, int right, struct RELATIVES& x, struct RELATIVES& y, T** IbsM)
{
	int i=left, j=right;
	x=relatives[(left+right)/2];
	do {
		while (( IbsM[relatives[i].i][relatives[i].j] > IbsM[x.i][x.j]) && (i<right)) i++;
		while (( IbsM[x.i][x.j] > IbsM[relatives[j].i][relatives[j].j]) && (j>left )) j--;
		if (i<=j) {
			y            = relatives[i];
            relatives[i] = relatives[j];
			relatives[j] =  y;
			i++;
			j--;
        }
    } while (i<=j);
	if (left < j) {
		qsortrelatives(relatives, left, j, x, y, IbsM);
	}
	if (i < right) {
		qsortrelatives(relatives, i, right, x, y, IbsM);
	}
}


// map2 nach rs-Nummern sortieren
void qsortmap2(struct MAP2 **map2,int left, int right, struct MAP2 *x, struct MAP2 *y)
{
	//sorts in ascending order
	int i,j;


	(*x).rs=NULL;(*x).rs=(char *)realloc((*x).rs,(strlen((*map2)[(left+right)/2].rs)+1)*sizeof(char));
	(*y).rs=NULL;

	i=left;j=right;

	strcpy((*x).rs,(*map2)[(left+right)/2].rs);
	(*x).line=(*map2)[(left+right)/2].line;

	do
	{
		while( strcmp((*map2)[i].rs,(*x).rs)< 0 && (i <right)) {i++;}
		while( strcmp((*x).rs,(*map2)[j].rs)< 0 && (j>left)) {j--;}

		if (i<=j)
		{
			(*y).rs=(char*)realloc((*y).rs,(strlen((*map2)[i].rs)+1)*sizeof(char));

			strcpy((*y).rs,(*map2)[i].rs);
			(*y).line=(*map2)[i].line;

			(*map2)[i].rs=(char*)realloc((*map2)[i].rs,(strlen((*map2)[j].rs)+1)*sizeof(char));

			strcpy((*map2)[i].rs,(*map2)[j].rs);
			(*map2)[i].line=(*map2)[j].line;

			(*map2)[j].rs=(char*)realloc((*map2)[j].rs,(strlen((*y).rs)+1)*sizeof(char));

			strcpy((*map2)[j].rs,(*y).rs);
			(*map2)[j].line=(*y).line;

			i++;j--;
		}
	}while (i<=j);

	if((*x).rs!=NULL){free((*x).rs);}

	if((*y).rs!=NULL){free((*y).rs);}
	(*x).rs=NULL;(*y).rs=NULL;

	if (left <j)
	{
		qsortmap2(map2,left,j, x, y);
	}
	if (i < right)
	{
		qsortmap2(map2,i,right, x, y);
	}

};


//PAA
void qsort(double *list,int left, int right)
{
	//sorts in ascending order
	int i,j;

	double x;
	double y;

	i=left;j=right;

	x=list[left+(right-left)/2];

	do
	{
		while( (list[i] < x) && (i<right)) {i++;}
		while( (x < list[j]) && (j>left)) {j--;}

		if (i<=j)
		{
			y=list[i];
			list[i]=list[j];
			list[j]= y;
			i++;j--;
		}
	}while (i<=j);

	if (left <j)
	{
		qsort(list,left,j);
	}
	if (i < right)
	{
		qsort(list,i,right);
	}
}


void qsort_cstring(char* *list, int left, int right) {
	int i = left;
	int j = right;
    char* P = list[left+(right-left)/2];
	while (i<j) {
		while (strcmp(list[i], P)<0) i++;
		while (strcmp(P, list[j])<0) j--;
        if (i<=j) {
            char* y = list[i];
            list[i] = list[j];
            list[j] = y;
            i++;
            j--;
        }
	}
	if (left <j) qsort_cstring(list, left , j);
	if (i<right) qsort_cstring(list, i, right);
}


int bsearch_cstring(char* *haystack, char* needle, int imin, int imax) {
  while (imax >= imin) {
      int imid = imin + ((imax - imin) / 2);
      if (strcmp(haystack[imid], needle) < 0) imin = imid + 1;
      else
      if (strcmp(haystack[imid], needle) > 0) imax = imid - 1;
      else
        return imid;
  }
  return -1;
}


int compare_cstring(const void* a, const void* b) {
    return strcmp(*(char**)a, *(char**)b);
}

int compare_int(const void* a, const void* b) {
  return *(int*)a - *(int*)b;
}



//PAA
void qsortplus(double *list,double *list2, int left, int right, int sign) //orders according to list, list2 identifies orginal ordering
{
	//sign==1: ascending
	int i,j;

	double x;
	double y;
	double z;

	i=left;j=right;

	x=list[(left+right)/2];

	do
	{
		while( (sign*list[i] < sign*x) && (i <right)) {i++;}
		while( (sign*x < sign*list[j]) && (j>left)) {j--;}

		if (i<=j)
		{
			y=list[i];
			list[i]=list[j];
			list[j]= y;


			z=list2[i];
			list2[i]=list2[j];
			list2[j]= z;

			i++;j--;
		}
	}while (i<=j);

	if (left <j)
	{
		qsortplus(list,list2,left,j,sign);
	}
	if (i < right)
	{
		qsortplus(list,list2,i,right,sign);
	}
};


int rsNumberSearch(char *rs, struct MAP2 *map2,int left, int right) // binäre Suche
{
	int l=left;
	int r=right;
	int x; // Position

	while(r >= l)
	{
		x=(l+r)/2;

		if(strcmp(rs, map2[x].rs)< 0 )/* kleiner? */
		{
	 		r=x-1;  /* Rechte Seite ist nicht mehr so interessant */
		}
		else

		{
	 		l=x+1;  /* Linke Seite ist nicht mehr so interessant */
		}

		if(strcmp(rs,map2[x].rs) == 0)
		{
	 		return map2[x].line;     /* Gefunden; x = Position*/
		}
	}
	return -1;
}


int rsNumberSearchInfo(char *rs, struct INFOMAP *infomap,int left, int right) // binäre Suche
{
	int l=left;
	int r=right;
	int x; // Position

	while(r >= l)
	{
		x=(l+r)/2;

		if(strcmp(rs, infomap[x].rs)< 0 )/* kleiner? */
		{
	 		r=x-1;  /* Rechte Seite ist nicht mehr so interessant */
		}
		else
		{
	 		l=x+1;  /* Linke Seite ist nicht mehr so interessant */
		}

		if(strcmp(rs,infomap[x].rs) == 0)
		{
	 		return infomap[x].line;     /* Gefunden; x = Position*/
		}
	}
	return -1;
}


void quicksortINFOMAP(struct INFOMAP **infomap, int left, int right, struct INFOMAP *x, struct INFOMAP *y)
{
	//sorts in ascending order
	int i,j;


	(*x).rs=NULL;(*x).rs=(char *)realloc((*x).rs,(strlen((*infomap)[(left+right)/2].rs)+1)*sizeof(char));
	if (!(*x).rs)
	{
		//errorfile << "memory allocation error in infox.rs\n";
		//cout << "memory allocation error in infox.rs\n";
		//errorfile.close();
		//exit(1);
	}
	(*x).gene=NULL;(*x).gene=(char *)realloc((*x).gene,(strlen((*infomap)[(left+right)/2].gene)+1)*sizeof(char));
	if (!(*x).gene)
	{
		//errorfile << "memory allocation error in infox.gene\n";
		//cout << "memory allocation error in infox.gene\n";
		//errorfile.close();
		//exit(1);
	}

	(*y).rs=NULL; //infoy.rs=(char *)realloc(infoy.rs,100*sizeof(char));
	(*y).gene=NULL; //infoy.gene=(char *)realloc(infoy.gene,1000*sizeof(char));

	i=left;j=right;

	strcpy((*x).rs,(*infomap)[(left+right)/2].rs);
	strcpy((*x).gene,(*infomap)[(left+right)/2].gene);
	strcpy((*x).chr,(*infomap)[(left+right)/2].chr);
	(*x).pos=(*infomap)[(left+right)/2].pos;
	(*x).location=(*infomap)[(left+right)/2].location;
	(*x).locationToGene=(*infomap)[(left+right)/2].locationToGene;
	(*x).codingStatus=(*infomap)[(left+right)/2].codingStatus;


	do
	{
		while( compareINFOMAP((*infomap)[i],(*x))==1 && (i <right)) {i++;}
		while( compareINFOMAP((*x),(*infomap)[j])==1 && (j>left)) {j--;}

		if (i<=j)
		{
			(*y).rs=(char*)realloc((*y).rs,(strlen((*infomap)[i].rs)+1)*sizeof(char));
			if (!(*y).rs)
			{
				//errorfile << "memory allocation error in infox.rs\n";
				//cout << "memory allocation error in infox.rs\n";
				//errorfile.close();
				//exit(1);
			}
			(*y).gene=(char*)realloc((*y).gene,(strlen((*infomap)[i].gene)+1)*sizeof(char));
			if (!(*y).gene)
			{
				//errorfile << "memory allocation error in infox.gene\n";
				//cout << "memory allocation error in infox.gene\n";
				//errorfile.close();
				//exit(1);
			}

			strcpy((*y).rs,(*infomap)[i].rs);
			strcpy((*y).gene,(*infomap)[i].gene);
			strcpy((*y).chr,(*infomap)[i].chr);
			(*y).pos=(*infomap)[i].pos;
			(*y).location=(*infomap)[i].location;
			(*y).locationToGene=(*infomap)[i].locationToGene;
			(*y).codingStatus=(*infomap)[i].codingStatus;

			(*infomap)[i].rs=(char*)realloc((*infomap)[i].rs,(strlen((*infomap)[j].rs)+1)*sizeof(char));
			if (!(*infomap)[i].rs)
			{
				//errorfile << "memory allocation error in infox.rs\n";
				//cout << "memory allocation error in infox.rs\n";
				//errorfile.close();
				//exit(1);
			}
			(*infomap)[i].gene=(char*)realloc((*infomap)[i].gene,(strlen((*infomap)[j].gene)+1)*sizeof(char));
			if (!(*infomap)[i].gene)
			{
				//errorfile << "memory allocation error in infox.gene\n";
				//cout << "memory allocation error in infox.gene\n";
				//errorfile.close();
				//exit(1);
			}

			strcpy((*infomap)[i].rs,(*infomap)[j].rs);
			strcpy((*infomap)[i].gene,(*infomap)[j].gene);
			strcpy((*infomap)[i].chr,(*infomap)[j].chr);
			(*infomap)[i].pos=(*infomap)[j].pos;
			(*infomap)[i].location=(*infomap)[j].location;
			(*infomap)[i].locationToGene=(*infomap)[j].locationToGene;
			(*infomap)[i].codingStatus=(*infomap)[j].codingStatus;

			(*infomap)[j].rs=(char*)realloc((*infomap)[j].rs,(strlen((*y).rs)+1)*sizeof(char));
			if (!(*infomap)[j].rs)
			{
				//errorfile << "memory allocation error in infox.rs\n";
				//cout << "memory allocation error in infox.rs\n";
				//errorfile.close();
				//exit(1);
			}
			(*infomap)[j].gene=(char*)realloc((*infomap)[j].gene,(strlen((*y).gene)+1)*sizeof(char));
			if (!(*infomap)[j].gene)
			{
				//errorfile << "memory allocation error in infox.gene\n";
				//cout << "memory allocation error in infox.gene\n";
				//errorfile.close();
				//exit(1);
			}

			strcpy((*infomap)[j].rs,(*y).rs);
			strcpy((*infomap)[j].gene,(*y).gene);
			strcpy((*infomap)[j].chr,(*y).chr);
			(*infomap)[j].pos=(*y).pos;
			(*infomap)[j].location=(*y).location;
			(*infomap)[j].locationToGene=(*y).locationToGene;
			(*infomap)[j].codingStatus=(*y).codingStatus;

			i++;j--;
		}
	}while (i<=j);

	free((*x).rs);free((*x).gene);
	if((*y).rs!=NULL){free((*y).rs);}
	if((*y).gene!=NULL){free((*y).gene);}

 	if (left <j)
 	{
		quicksortINFOMAP(infomap,left,j, x, y);

 	}

 	if (i < right)
 	{
		quicksortINFOMAP(infomap,i,right, x, y);
 	}
};


double **calloc2Ddouble(int m, int n)
{
	int i;
	double **matrix;
	matrix=(double **)calloc(m,sizeof(double *));


	if (!matrix){printf("ALLOCATION ERROR in calloc2Ddouble\n");exit(1);}
	for (i=0;i<m;i++)
	{
		matrix[i]=(double *)calloc(n, sizeof(double));
		if (!matrix[i]) {printf("ALLOCATION ERROR in calloc2Ddouble\n");exit(1);}
	}
	return (matrix);
}


float **calloc2Dfloat(int m, int n)
{
	int i;
	float **matrix;
	matrix=(float **)calloc(m,sizeof(float *));
	if (!matrix){printf("ALLOCATION ERROR in calloc2Dfloat\n");exit(1);}
	for (i=0;i<m;i++)
	{
		matrix[i]=(float *)calloc(n, sizeof(float));
		if (!matrix[i]) {printf("ALLOCATION ERROR in calloc2Dfloat\n");exit(1);}
	}
	return (matrix);
}

float ***calloc3Dfloat(int m, int n,int l)
{
	 int i,j;
     float ***matrix;
     matrix=(float ***)calloc(m,sizeof(float **));
     if (!matrix){printf("ALLOCATION ERROR in calloc3Dfloat\n");exit(1);}
     for (i=0;i<m;i++)
         {matrix[i]=(float **)calloc(n, sizeof(float *));
          if (!matrix[i]) {printf("ALLOCATION ERROR in calloc3Dfloat\n");exit(1);}
          for (j=0;j<n;j++)
              {matrix[i][j]=(float *)calloc(l, sizeof(float));
               if (!matrix[i][j]){printf("ALLOCATION ERROR in calloc3Dfloat\n");exit(1);}
              }
         }
     return (matrix);
}

int **calloc2Dint(int m, int n)
{
	int i;
	int **matrix;
	matrix=(int **)calloc(m,sizeof(int *));
	if (!matrix){printf("ALLOCATION ERROR in calloc2Dint\n");exit(1);}
	for (i=0;i<m;i++)
	{
		matrix[i]=(int *)calloc(n, sizeof(int));
		if (!matrix[i]) {printf("ALLOCATION ERROR in calloc2Dint\n");exit(1);}
	}
	return (matrix);
}


void free2Ddouble(double **matrix, int m)
{
	int i;

	for (i=0;i<m;i++)
	{
		free(matrix[i]);matrix[i]=NULL;
	}
	free(matrix);
	matrix=NULL;
}

void free2Dint(int **matrix, int m)
{
  int i;

  for (i=0;i<m;i++)
    {
      free(matrix[i]);matrix[i]=NULL;
    }
  free(matrix);
  matrix=NULL;
}


// Transformieren der Matrix
int mysvd(int n, int m, double ***U, double ***V, double **S, int la, int lv)
{
	int i,j,k,l,itr,ier, itmax=40000;
	double f, g, h, rmax, s, r1, r2, c, x, y, z, aeps, eps=1.e-20;
	double* e=NULL;

	if(n>m || n<=0 || m<=0 || n>la || n>lv) return 105;
	ier=0;

	// Reduction to Bidiagonal form using Householder transformations
	g=0.0; rmax=0.0;
	e=(double *) malloc(n*sizeof(double));
	for(i=0; i<n; ++i)
	{
		// Off-diagonal elements of bidiagonal form
		e[i]=g;
		s=0.0;

		for(j=i; j<m; ++j) s=s+pow(*((*((*U)+j))+i),2);
		if(s <= 0.0)
		{
			// transformation not required
			g=0.0;
		}
		else
		{
			f= *((*((*U)+i))+i);
			g=sqrt(s);
			if(f>=0.0) g=-g;
			h=f*g-s;

			*((*((*U)+i))+i) = f-g;

			for(j=i+1; j<n; ++j)
			{
				s=0.0;

				for(k=i; k<m; ++k) s=s+(*((*((*U)+k))+i))*(*((*((*U)+k))+j));
				f=s/h;

				for(k=i; k<m; ++k) *((*((*U)+k))+j)= *((*((*U)+k))+j)+f*(*((*((*U)+k))+i));
			}
		}

		// Diagonal elements of bidiagonal form
		(*S)[i]=g;
		s=0.0;

		for(j=i+1; j<n; ++j) s=s+pow(*((*((*U)+i))+j),2);

		if(s<= 0.0)
		{
			g=0.0;
		}
		else
		{
			f= *((*((*U)+i))+(i+1));
			g=sqrt(s);
			if(f>= 0.0) g=-g;
			h=f*g-s;

			*((*((*U)+i))+(i+1))=f-g;
			for(j=i+1; j<n; ++j) e[j]=(*((*((*U)+i))+j))/h;

			for(j=i+1; j<m; ++j)
			{
				s=0.0;
				for(k=i+1; k<n; ++k) s=s+(*((*((*U)+j))+k))*(*((*((*U)+i))+k));
				for(k=i+1; k<n; ++k) *((*((*U)+j))+k) = *((*((*U)+j))+k)+s*e[k];
			}
		}
		r1=fabs((*S)[i])+fabs(e[i]);
		if(r1 > rmax) rmax=r1;
	}

	// Accumulation of right hand transformation in array V
	for(i=n-1; i>=0; --i)
	{
		if(g != 0.0)
		{
			h=g*(*((*((*U)+i))+(i+1)));
			for(j=i+1; j<n; ++j) *((*((*V)+j))+i)=(*((*((*U)+i))+j))/h;

			for(j=i+1; j<n; ++j)
			{
				s=0.0;
				for(k=i+1; k<n; ++k) s=s+(*((*((*U)+i))+k))*(*((*((*V)+k))+j));
				for(k=i+1; k<n; ++k) *((*((*V)+k))+j)=*((*((*V)+k))+j)+s*(*((*((*V)+k))+i));
			}
		}

		for(j=i+1; j<n; ++j)
		{
			*((*((*V)+i))+j)=0.0; *((*((*V)+j))+i)=0.0;
		}
		*((*((*V)+i))+i)=1;
		g= e[i];
	}

	// Accumulation of left hand transformation overwritten on matrix A
	for(i=n-1; i>=0; --i)
	{
		g=(*S)[i];
		for(j=i+1; j<n; ++j) *((*((*U)+i))+j)=0.0;
		if(g != 0.0)
		{
			h=g*(*((*((*U)+i))+i));

			for(j=i+1; j<n; ++j)
			{
				s=0.0;
				for(k=i+1; k<m; ++k) s=s+(*((*((*U)+k))+i))*(*((*((*U)+k))+j));
				f=s/h;
				for(k=i; k<m; ++k) *((*((*U)+k))+j)=*((*((*U)+k))+j)+f*(*((*((*U)+k))+i));
			}

			for(j=i; j<m; ++j) *((*((*U)+j))+i)=(*((*((*U)+j))+i))/g;
		}
		else
		{
			for(j=i; j<m; ++j) *((*((*U)+j))+i)=0.0;
		}
		*((*((*U)+i))+i) = *((*((*U)+i))+i)+1;
	}

	// Diagonalisation of the bidiagonal form
	aeps=eps*rmax;
	// Loop over the singular values

	for(k=n-1; k>=0; --k)
	{
		// The QR transformation
		for(itr=1; itr<=itmax; ++itr)
		{
			// Test for splitting
			for(l=k; l>=0; --l)
			{
				if(fabs(e[l]) < aeps) goto split;
				if(fabs((*S)[l-1]) < aeps) break;
			}

			// cancellation of E[L] if L>1
			c=0.0; s=1.0;

			for(i=l; i<=k; ++i)
			{
				f=s*e[i];
				e[i] = c*e[i];
				if(fabs(f) < aeps) goto split;

				g=(*S)[i];
				(*S)[i]=sqrt(f*f+g*g);
				c=g/((*S)[i]);
				s=-f/((*S)[i]);

				for(j=0; j<m; ++j)
				{
					r1= *((*((*U)+j))+(l-1));
					r2= *((*((*U)+j))+i);
					*((*((*U)+j))+(l-1))=r1*c+r2*s;
					*((*((*U)+j))+i)=c*r2-s*r1;
				}
			}

			split:;  z=(*S)[k];
			if(l == k)
			{
				// QR iteration has converged
				if(z < 0.0)
				{
					(*S)[k] = -z;
					for(j=0; j<n; ++j) *((*((*V)+j))+k)=-(*((*((*V)+j))+k));
				}
				break;
			}

			if(itr==itmax) {ier=12; break;}

			// calculating shift from bottom 2x2 minor
			x=(*S)[l];
			y=(*S)[k-1];
			g=e[k-1];
			h=e[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.*h*y);
			g=sqrt(1.+f*f);
			if(f < 0.0) g=-g;
			f=((x-z)*(x+z)+h*(y/(f+g)-h))/x;

			// next QR transformation
			c=1.0; s=1.0;
			// Given's rotation
			for(i=l+1; i<=k; ++i)
			{
				g=e[i];
				y=(*S)[i];
				h=s*g;
				g=c*g;
				e[i-1]=sqrt(f*f+h*h);
				c=f/e[i-1];
				s=h/e[i-1];
				f=c*x+s*g;
				g=c*g-s*x;
				h=s*y;
				y=c*y;

				for(j=0; j<n; ++j)
				{
					x=*((*((*V)+j))+(i-1));
					z=*((*((*V)+j))+i);
					*((*((*V)+j))+(i-1))=c*x+s*z;
					*((*((*V)+j))+i)=c*z-s*x;
				}

				(*S)[i-1]=sqrt(f*f+h*h);
				if((*S)[i-1] != 0.0)
				{
					c=f/((*S)[i-1]);
					s=h/((*S)[i-1]);
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for(j=0; j<m; ++j)
				{
					y= *((*((*U)+j))+(i-1));
					z= *((*((*U)+j))+i);
					*((*((*U)+j))+(i-1)) = c*y+s*z;
					*((*((*U)+j))+i) = c*z-s*y;
				}
			}
			e[l]=0.0;
			e[k]=f;
			(*S)[k]=x;
		}
	}
       	if (e != NULL) { free(e); e=NULL; }
	return ier;
}


// Multiplizieren einer Matrix
int MatrixMult(double **A, int Adim1, int Adim2, double **B, int Bdim1, int Bdim2, double ***C)
{
	int i,j,l=0;double helfer=0;

	if (Adim2!=Bdim1){printf("Matrix Types Incompatible!\n"); exit(1);}

	for (i=0;i<Adim1;i++)
	{
		for (l=0;l<Bdim2;l++)
		{
			helfer=0;
			for (j=0;j<Adim2;j++)
			{
				helfer=helfer+(*((*(A+i))+j))*(*((*(B+j))+l));
			}
			*((*((*C)+i))+l)=helfer;
		}
	}
	return 0;
} //end MatrixMult


int MatrixMultMod(double **A, int Adim1, int Adim2, double **B, int Bdim1, int Bdim2, double **C)
{

	int i,j,l=0;double helfer=0;
	if (Adim2!=Bdim1){printf("Matrix Types Incompatible!\n"); exit(1);}

	for (i=0;i<Adim1;i++)
	{
		for (l=0;l<Bdim2;l++)
		{
			helfer=0;
			for (j=0;j<Adim2;j++)
			{
				helfer=helfer+(*((*(A+i))+j))*(*((*(B+j))+l));
			}
			// ((*((*C)+i))+l)=helfer;
			C[i][l]=helfer;
		}
	}
	return 0;
} //end MatrixMult


// Multiplizieren einer Matrix wenn bekannt ist das Ergebnis symmetrisch
int MatrixMultSym(double **A, int Adim1, int Adim2, double **B, int Bdim1, int Bdim2, double **C)
{
	int i,j,l=0;double helfer=0;
	if (Adim2!=Bdim1){printf("Matrix Types Incompatible!\n"); exit(1);}

	for (i=0;i<Adim1;i++)
	{
		for (l=0;l<Bdim2;l++)
		{
			if(l>=i)
		    {
				helfer=0;
				for (j=0;j<Adim2;j++)
				{
					helfer=helfer+(*((*(A+i))+j))*(*((*(B+j))+l));
				}
				C[i][l]=helfer;
			}
			else
			{
				C[i][l]=C[l][i];
			}
		}
	}
	return 0;
} //end MatrixMult


// Invertiern einer Matrix nach dem Dwyer Algorithmus
int DwyerInv(int n, double **MM, double **D, double **T, double **U, double **Ut, double **sumPP, double **sumPJ, double **sumPK, double **MMinv, double **Minv)
{
	int i,j,p,k;
	int dummyMatrix = 0;
	int dummy = 1;

	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
		{
			if (i == j)
			{
				D[i][j] = 1;
			}
			else
			{
				D[i][j] = 0;
			}
			T[i][j] = 0;
			U[i][j] = 0;
			sumPP[i][j] = 0;
			sumPJ[i][j] = 0;
			sumPK[i][j] = 0;
		}
	}


	for (p = 0; p< n; p++)
	{
		for (j = 0; j< n; j++)
		{
			if (p == j)
			{
				for (i = 0; i< p; i++)
				{
					sumPP[0][p] -= (T[i][p] * T[i][p]);
				}

//cout << MM[p][p] << " " << sumPP[0][p] << endl;
				if (MM[p][p] + sumPP[0][p] >= 0)
				{
					T[p][p] = sqrt(MM[p][p] + sumPP[0][p]);
				}
				else
				{
					return 0;
				}


				// Division durch Null?
				if (fabs(T[p][p]) <= 1.0e-4)
				{
					return 0;
				}
				else
				{
					dummy = 1;
				}
			}
			else if (j > p)
			{
				for (i = 0; i< p; i++)
				{
					sumPJ[p][j] -= (T[i][p] * T[i][j]);
				}
				T[p][j] = (MM[p][j] + sumPJ[p][j])/ T[p][p];
			}
			else
			{
				T[p][j] = 0;
			}
		}


		// k <= p
		for (k = 0; k<= p; k++)
		{
			for (i = 0; i< p; i++)
			{
				sumPK[p][k] -= (T[i][p]*U[i][k]);
			}
			U[p][k] = (D[p][k] + sumPK[p][k])/T[p][p];
		}
	}



	// Transponieren von U
	for (i = 0; i < n; i++)
	{
		for (j = 0; j< n; j++)
		{
			Ut[i][j] = U[j][i];
		}
	}

	// Matrixmultiplikation  von Ut*U
	dummyMatrix=MatrixMultSym(Ut,n,n,U,n,n,MMinv);

	for (i = 0; i < n; i++)
	{
		for (j = 0; j< n; j++)
		{
			//*((*((*Minv)+i))+j)=MMinv[i][j];
			Minv[i][j]=MMinv[i][j];
		}
	}

	return dummy;
}


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
						double **D, double **T, double **U, double **Ut, double **sumPP, double **sumPJ, double **sumPK, double **MMinv,int test,int thread,
						int npplqc, int* PPLMap, uint64_t*** BinSNPs, struct PPLLOCATION* PplLocations,int nMc,struct STATplus result, int maxIndexCov, int liabilityCut, int singleMarkerTest,
						int covNum, int covariancematrix,int df_L1,int df_L2, int dosage, float ***genoWeights)
{
	int indicator[28+maxIndexCov]; //which parameters (beta0-beta26) is sex a covaitate? (27) and covariates (28-37) are used?
	int numberit[28+maxIndexCov];//numbers the used indicators
	int number2indi[28+maxIndexCov]; //retrieve kth actually used parameter
	//struct STATplus result; //log-Likelihood and all betas
	double logOld,logNew;
	double betaNew[28+maxIndexCov];
	double betaOld[28+maxIndexCov];

	int i,k,j,jj,f=0;
	int withcov=0; //are there covariates?
	double eps=0.000001; //convergence criterion
	int it; //counts number of iterations
	int maxit=10; // maximum number of iterations
	int n; //number of individuals used
	double AA, AB, BB, CC=0, CD=0, DD=0, EE=0, EF=0, FF=0; //joint genocounts cases and control for SNPs 1,2,3
	double help1,help2,help3,total;
	int complete=0;
	double currentval[6]; //x-values (main effects) for current individual
	double g1,g2,g3; // indicates if genotypes for SNPs 1,2,3, respectively, are needed
	int snps=0; //1,2 or 3 (highest snp used)
	int x1; int x1D; //indicates if parameter is used
	int x2; int x2D; //indicates if parameter is used
	int x1x2; int x1x2D; int x1Dx2; int x1Dx2D; //indicates if parameter is used
	int x3; int x3D; //indicates if parameter is used
	int x1x3; int x1x3D; int x1Dx3; int x1Dx3D; //indicates if parameter is used
	int x2x3; int x2x3D; int x2Dx3; int x2Dx3D; //indicates if parameter is used
	int x1x2x3;int x1x2x3D;int x1x2Dx3;int x1x2Dx3D;int x1Dx2x3;int x1Dx2x3D;int x1Dx2Dx3;int x1Dx2Dx3D; //indicates if parameter is used
	double exponent;int dummy;
	double AA_Ca=0,AA_Co=0,AB_Ca=0,AB_Co=0,BB_Ca=0,BB_Co=0,CC_Ca=0,CC_Co=0,CD_Ca=0,CD_Co=0,DD_Ca=0,DD_Co=0,EE_Ca=0,EE_Co=0,EF_Ca=0,EF_Co=0,FF_Ca=0,FF_Co=0;
	double casecounts_ab[3][3], controlcounts_ab[3][3]; //2-dim genocounts snps 1,2
	double casecounts_ac[3][3], controlcounts_ac[3][3]; //2-dim genocounts snps 1,3
	double casecounts_bc[3][3], controlcounts_bc[3][3]; //2-dim genocounts snps 2,3
	int ncase=0; int ncontrol=0;
	double ncased=0; double ncontrold=0;
	double newf=f;
	int quick=1;
	double rss = 0;
	float liabilityOR[2][2];

	//init(&result,0,nMc);

	int dfUnchanged=1;
	currentval[0]=0;currentval[2]=0;currentval[4]=0;
	currentval[1]=0.5;currentval[3]=0.5;currentval[5]=0.5;

	x1=x[1];
	x1D=x[2];
	x2=x[3];
	x2D=x[4];
	x1x2=x[5];
	x1x2D=x[6];
	x1Dx2=x[7];
	x1Dx2D=x[8];
	x3=x[9];
	x3D=x[10];
	x1x3=x[11];
	x1x3D=x[12];
	x1Dx3=x[13];
	x1Dx3D=x[14];
	x2x3=x[15];
	x2x3D=x[16];
	x2Dx3=x[17];
	x2Dx3D=x[18];
	x1x2x3=x[19];
	x1x2x3D=x[20];
	x1x2Dx3=x[21];
	x1x2Dx3D=x[22];
	x1Dx2x3=x[23];
	x1Dx2x3D=x[24];
	x1Dx2Dx3=x[25];
	x1Dx2Dx3D=x[26];

	//printf("135 %d %d %d\n",x1,x2,x1x2);

	if(xtype>0 && !female && xtype!=3) //xtype 3 is single marker x chromosome, both sexes
	{
	  if(x[2] || x[4] || x[6] || x[7] || x[8] || x[10] || x[12] || x[13] || x[14] || x[16] || x[17] || x[18] || x[20] || x[21] || x[22] || x[23] || x[24] || x[25] || x[26])
	    {
	      dfUnchanged=0;
	    }
		x1D=0;x2D=0;x3D=0;
		x1Dx2=0;x1x2D=0;x1Dx2D=0;
		x1Dx3=0;x1x3D=0;x1Dx3D=0;
		x2Dx3=0;x2x3D=0;x2Dx3D=0;
		x1x2x3D=0;x1x2Dx3=0;x1x2Dx3D=0;x1Dx2x3=0;x1Dx2x3D=0;x1Dx2Dx3=0;x1Dx2Dx3D=0;
	}

	result.sc=0;

	if(xtype==3 && !male && !female) //NEU
	{
		sexcov=1;
	}
	else if(xtype==3) //NEU
	{
		sexcov=0;
	}
	if(xtype>0 && xtype <3)
	{
		sexcov=0;
	}


	//INITIALIZE INDICATOR
	indicator[0]=1; //always with intercept beta_0
	//exit(1);
	result.in[0]=1;

	for(i=1;i<28+maxIndexCov;i++){indicator[i]=0;result.in[i]=-1;}
	indicator[27]=sexcov; //sex as covariate yes/no
	if(indicator[27]){result.in[27]=1;}

	if(N==2) // 2-marker
	{
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				casecounts_ab[i][j]= casecounts5[i][j][0][0][0];
				controlcounts_ab[i][j]= controlcounts5[i][j][0][0][0];
			}
		}


		AA_Ca=casecounts_ab[0][0]+casecounts_ab[0][1]+casecounts_ab[0][2];
		AA_Co=controlcounts_ab[0][0]+controlcounts_ab[0][1]+controlcounts_ab[0][2];
		AB_Ca=casecounts_ab[1][0]+casecounts_ab[1][1]+casecounts_ab[1][2];
		AB_Co=controlcounts_ab[1][0]+controlcounts_ab[1][1]+controlcounts_ab[1][2];
		BB_Ca=casecounts_ab[2][0]+casecounts_ab[2][1]+casecounts_ab[2][2];
		BB_Co=controlcounts_ab[2][0]+controlcounts_ab[2][1]+controlcounts_ab[2][2];

		CC_Ca=casecounts_ab[0][0]+casecounts_ab[1][0]+casecounts_ab[2][0];
		CC_Co=controlcounts_ab[0][0]+controlcounts_ab[1][0]+controlcounts_ab[2][0];
		CD_Ca=casecounts_ab[0][1]+casecounts_ab[1][1]+casecounts_ab[2][1];
		CD_Co=controlcounts_ab[0][1]+controlcounts_ab[1][1]+controlcounts_ab[2][1];
		DD_Ca=casecounts_ab[0][2]+casecounts_ab[1][2]+casecounts_ab[2][2];
		DD_Co=controlcounts_ab[0][2]+controlcounts_ab[1][2]+controlcounts_ab[2][2];

		EE_Ca=0;EE_Co=0;EF_Ca=0;EF_Co=0;FF_Ca=0;FF_Co=0;
	}
	else if(N==3)// 3-marker
	{
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				casecounts_ab[i][j]=0;
				controlcounts_ab[i][j]=0;

				casecounts_ac[i][j]=0;
				controlcounts_ac[i][j]=0;
				casecounts_bc[i][j]=0;
				controlcounts_bc[i][j]=0;
			}
		}

		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				for(k=0;k<3;k++)
				{
					casecounts_ab[i][j]+=casecounts5[i][j][k][0][0];
					controlcounts_ab[i][j]+=controlcounts5[i][j][k][0][0];
					casecounts_ac[i][k]+=casecounts5[i][j][k][0][0];
					controlcounts_ac[i][k]+=controlcounts5[i][j][k][0][0];
					casecounts_bc[j][k]+=casecounts5[i][j][k][0][0];
					controlcounts_bc[j][k]+=controlcounts5[i][j][k][0][0];

					if(i==0)
					{
						AA_Ca+=casecounts5[i][j][k][0][0];
						AA_Co+=controlcounts5[i][j][k][0][0];
					}
					if(i==1)
					{
						AB_Ca+=casecounts5[i][j][k][0][0];
						AB_Co+=controlcounts5[i][j][k][0][0];
					}
					if(i==2)
					{
						BB_Ca+=casecounts5[i][j][k][0][0];
						BB_Co+=controlcounts5[i][j][k][0][0];
					}
					if(j==0)
					{
						CC_Ca+=casecounts5[i][j][k][0][0];
						CC_Co+=controlcounts5[i][j][k][0][0];
					}
					if(j==1)
					{
						CD_Ca+=casecounts5[i][j][k][0][0];
						CD_Co+=controlcounts5[i][j][k][0][0];
					}
					if(j==2)
					{
						DD_Ca+=casecounts5[i][j][k][0][0];
						DD_Co+=controlcounts5[i][j][k][0][0];
					}
					if(k==0)
					{
						EE_Ca+=casecounts5[i][j][k][0][0];
						EE_Co+=controlcounts5[i][j][k][0][0];
					}
					if(k==1)
					{
						EF_Ca+=casecounts5[i][j][k][0][0];
						EF_Co+=controlcounts5[i][j][k][0][0];
					}
					if(k==2)
					{
						FF_Ca+=casecounts5[i][j][k][0][0];
						FF_Co+=controlcounts5[i][j][k][0][0];
					}  // if
				} //k
			}	//j
		}   // i
	} //if N==3

	if(N>1)
	{
		AA=AA_Ca+AA_Co;AB=AB_Ca+AB_Co;BB=BB_Ca+BB_Co;
		CC=CC_Ca+CC_Co;CD=CD_Ca+CD_Co;DD=DD_Ca+DD_Co;
		EE=EE_Ca+EE_Co;EF=EF_Ca+EF_Co;FF=FF_Ca+FF_Co;
	}
	else
	{
		AA=counts.AA_Ca+counts.AA_Co;
		AB=counts.AB_Ca+counts.AB_Co;
		BB=counts.BB_Ca+counts.BB_Co;
	}

	if (liabilityCut == 0 && 0)
	{
		if( (AA==0 && AB==0) || (AA ==0 && BB==0) || (AB ==0 && BB==0) ) // SNP 1 monomorph
		{
			x1=0;x1D=0;
			x1x2=0;x1x2D=0;x1Dx2=0;x1Dx2D=0;
			x1x3=0;x1x3D=0;x1Dx3=0;x1Dx3D=0;
			x1x2x3=0;x1x2x3D=0;x1x2Dx3=0;x1x2Dx3D=0;x1Dx2x3=0;x1Dx2x3D=0;x1Dx2Dx3=0;x1Dx2Dx3D=0;
		}
		else if ( AA==0 || AB==0 || BB==0) // no dominance
		{
			x1D=0;
			x1Dx2=0;x1Dx2D=0;
			x1Dx3=0;x1Dx3D=0;
			x1Dx2x3=0;x1Dx2x3D=0;x1Dx2Dx3=0;x1Dx2Dx3D=0;
		}

		if( (CC==0 && CD==0) || (CC ==0 && DD==0) || (CD ==0 && DD==0) ) // SNP 2 monomorph
		{
			x2=0;x2D=0;
			x1x2=0;x1x2D=0;x1Dx2=0;x1Dx2D=0;
			x2x3=0;x2x3D=0;x2Dx3=0;x2Dx3D=0;
			x1x2x3=0;x1x2x3D=0;x1x2Dx3=0;x1x2Dx3D=0;x1Dx2x3=0;x1Dx2x3D=0;x1Dx2Dx3=0;x1Dx2Dx3D=0;
		}
		else if (CC==0 || CD==0 || DD==0 ) // only 2 genotypes
		{
			x2D=0;
			x1x2D=0;x1Dx2D=0;
			x2Dx3=0;x2Dx3D=0;
			x1x2Dx3=0;x1x2Dx3D=0;x1Dx2Dx3=0;x1Dx2Dx3D=0;
		}

		if( (EE==0 && EF==0) || (EE == 0 && FF==0) || (EF ==0 && FF==0) ) // SNP 3 monomorph
		{
			x3=0;x3D=0;
			x1x3=0;x1x3D=0;x1Dx3=0;x1Dx3D=0;
			x2x3=0;x2x3D=0;x2Dx3=0;x2Dx3D=0;
			x1x2x3=0;x1x2x3D=0;x1x2Dx3=0;x1x2Dx3D=0;x1Dx2x3=0;x1Dx2x3D=0;x1Dx2Dx3=0;x1Dx2Dx3D=0;
		}
		else if (EE==0 || EF==0 || FF==0) // only 2 genotypes
		{
			x3D=0;
			x2x3D=0;x2Dx3D=0;
			x1x3D=0;x1Dx3D=0;
			x1x2x3D=0;x1x2Dx3D=0;x1Dx2x3D=0;x1Dx2Dx3D=0;
		}
	}

	//now check further interaction terms:

	if(N>1 && 0)
	{
		total=0;
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				total+=casecounts_ab[i][j];
				total+=controlcounts_ab[i][j];
			}
		}

		help1=casecounts_ab[0][0]+controlcounts_ab[0][0]+casecounts_ab[2][2]+controlcounts_ab[2][2];
		help2=casecounts_ab[2][0]+controlcounts_ab[2][0]+casecounts_ab[0][2]+controlcounts_ab[0][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total)
		{
			x1x2=0;x1x2x3=0;x1x2x3D=0;
		}
		help1=casecounts_ab[0][0]+controlcounts_ab[0][0]+
		casecounts_ab[2][1]+controlcounts_ab[2][1]+casecounts_ab[0][2]+controlcounts_ab[0][2];
		help2=casecounts_ab[2][0]+controlcounts_ab[2][0]+
		casecounts_ab[0][1]+controlcounts_ab[0][1]+casecounts_ab[2][2]+controlcounts_ab[2][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total)
		{
			x1x2D=0;x1x2Dx3=0;x1x2Dx3D=0;
		}


		help1=casecounts_ab[0][0]+controlcounts_ab[0][0]+
		casecounts_ab[1][2]+controlcounts_ab[1][2]+casecounts_ab[2][0]+controlcounts_ab[2][0];
		help2=casecounts_ab[0][2]+controlcounts_ab[0][2]+
		casecounts_ab[1][0]+controlcounts_ab[1][0]+casecounts_ab[2][2]+controlcounts_ab[2][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total)
		{
			x1Dx2=0;x1Dx2x3=0;x1Dx2x3D=0;
		}


		help1=casecounts_ab[0][0]+controlcounts_ab[0][0]+
		casecounts_ab[0][2]+controlcounts_ab[0][2]+casecounts_ab[1][1]+controlcounts_ab[1][1]+
		casecounts_ab[2][0]+controlcounts_ab[2][0]+casecounts_ab[2][2]+controlcounts_ab[2][2];

		help2=total-help1;

		if(help1==total || help2==total)
		{
			x1Dx2D=0;x1Dx2Dx3=0;x1Dx2Dx3D=0;
		}
	} //end if N>1


	if(N>2 && 0)
	{
		//snps 1 and 3
		total=0;
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				total+=casecounts_ac[i][j];
				total+=controlcounts_ac[i][j];
			}
		}

		help1=casecounts_ac[0][0]+controlcounts_ac[0][0]+casecounts_ac[2][2]+controlcounts_ac[2][2];
		help2=casecounts_ac[2][0]+controlcounts_ac[2][0]+casecounts_ac[0][2]+controlcounts_ac[0][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total){x1x3=0;x1x2x3=0;x1x2Dx3=0;}


		help1=casecounts_ac[0][0]+controlcounts_ac[0][0]+ casecounts_ac[2][1]+controlcounts_ac[2][1]+casecounts_ac[0][2]+
			controlcounts_ac[0][2];
		help2=casecounts_ac[2][0]+controlcounts_ac[2][0]+ casecounts_ac[0][1]+controlcounts_ac[0][1]+casecounts_ac[2][2]+
			controlcounts_ac[2][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total){x1x3D=0;x1x2x3D=0;x1x2Dx3D=0;}


		help1=casecounts_ac[0][0]+controlcounts_ac[0][0]+ casecounts_ac[1][2]+controlcounts_ac[1][2]+casecounts_ac[2][0]+
			controlcounts_ac[2][0];
		help2=casecounts_ac[0][2]+controlcounts_ac[0][2]+ casecounts_ac[1][0]+controlcounts_ac[1][0]+casecounts_ac[2][2]+
			controlcounts_ac[2][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total){x1Dx3=0;x1Dx2x3=0;x1Dx2Dx3=0;}


		help1=casecounts_ac[0][0]+controlcounts_ac[0][0]+ casecounts_ac[0][2]+controlcounts_ac[0][2]+casecounts_ac[1][1]+
			controlcounts_ac[1][1]+
		casecounts_ac[2][0]+controlcounts_ac[2][0]+casecounts_ac[2][2]+controlcounts_ac[2][2];

		help2=total-help1;

		if(help1==total || help2==total){x1Dx3D=0;x1Dx2x3D=0;x1Dx2Dx3D=0;}


		//snps 2 and 3
		total=0;
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				total+=casecounts_bc[i][j];
				total+=controlcounts_bc[i][j];
			}
		}


		help1=casecounts_bc[0][0]+controlcounts_bc[0][0]+casecounts_bc[2][2]+controlcounts_bc[2][2];
		help2=casecounts_bc[2][0]+controlcounts_bc[2][0]+casecounts_bc[0][2]+controlcounts_bc[0][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total){x2x3=0;x1x2x3=0;x1Dx2x3=0;}


		help1=casecounts_bc[0][0]+controlcounts_bc[0][0]+ casecounts_bc[2][1]+controlcounts_bc[2][1]+casecounts_bc[0][2]+controlcounts_bc[0][2];
		help2=casecounts_bc[2][0]+controlcounts_bc[2][0]+ casecounts_bc[0][1]+controlcounts_bc[0][1]+casecounts_bc[2][2]+controlcounts_bc[2][2];
		help3=total-help2-help1;


		if(help1==total || help2==total || help3==total){x2x3D=0;x1x2x3D=0;x1Dx2x3D=0;}


		help1=casecounts_bc[0][0]+controlcounts_bc[0][0]+ casecounts_bc[1][2]+controlcounts_bc[1][2]+casecounts_bc[2][0]+controlcounts_bc[2][0];
		help2=casecounts_bc[0][2]+controlcounts_bc[0][2]+ casecounts_bc[1][0]+controlcounts_bc[1][0]+casecounts_bc[2][2]+controlcounts_bc[2][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total){x2Dx3=0;x1x2Dx3=0;x1Dx2Dx3=0;}


		help1=casecounts_bc[0][0]+controlcounts_bc[0][0]+ casecounts_bc[0][2]+controlcounts_bc[0][2]+casecounts_bc[1][1]+controlcounts_bc[1][1]+
		casecounts_bc[2][0]+controlcounts_bc[2][0]+casecounts_bc[2][2]+controlcounts_bc[2][2];

		help2=total-help1;

		if(help1==total || help2==total){x2Dx3D=0;x1x2Dx3D=0;x1Dx2Dx3D=0;}

		//now 3-fold terms
		total=0;
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				for(jj=0;jj<3;jj++)
				{
					total+=casecounts5[i][j][jj][0][0];
					total+=controlcounts5[i][j][jj][0][0];
				}
			}
		}

		if(x1x2x3) //update if always equal
		{
			help1= casecounts5[0][2][0][0][0]+casecounts5[2][0][0][0][0]+casecounts5[0][0][2][0][0]+casecounts5[2][2][2][0][0]+controlcounts5[0][2][0][0][0]+controlcounts5[2][0][0][0][0]+controlcounts5[0][0][2][0][0]+controlcounts5[2][2][2][0][0];

			help2= casecounts5[0][0][0][0][0]+casecounts5[2][2][0][0][0]+casecounts5[0][2][2][0][0]+casecounts5[2][0][2][0][0]+controlcounts5[0][0][0][0][0]+controlcounts5[2][2][0][0][0]+controlcounts5[0][2][2][0][0]+controlcounts5[2][0][2][0][0];

			help3=total-help1-help2;

			if(help1==total || help2==total || help3==total){x1x2x3=0;}
		}

		if(x1x2x3D) //update if always equal
		{
			help1= casecounts5[0][2][0][0][0]+casecounts5[2][0][0][0][0]+casecounts5[0][0][1][0][0]+casecounts5[2][2][1][0][0]+casecounts5[0][2][2][0][0]+casecounts5[2][0][2][0][0]+controlcounts5[0][2][0][0][0]+	controlcounts5[2][0][0][0][0]+controlcounts5[0][0][1][0][0]+controlcounts5[2][2][1][0][0]+controlcounts5[0][2][2][0][0]+controlcounts5[2][0][2][0][0];

			help2= casecounts5[0][0][0][0][0]+casecounts5[2][2][0][0][0]+casecounts5[0][2][1][0][0]+casecounts5[2][0][1][0][0]+
				casecounts5[0][0][2][0][0]+casecounts5[2][2][2][0][0]+controlcounts5[0][0][0][0][0]+
				controlcounts5[2][2][0][0][0]+controlcounts5[0][2][1][0][0]+controlcounts5[2][0][1][0][0]+
				controlcounts5[0][0][2][0][0]+controlcounts5[2][2][2][0][0];

			help3=total-help1-help2;

			if(help1==total || help2==total || help3==total){x1x2x3D=0;}
		}

		if(x1x2Dx3) //update if always equal
		{
			help1= casecounts5[0][1][0][0][0]+casecounts5[2][0][0][0][0]+casecounts5[2][2][0][0][0]+casecounts5[0][0][2][0][0]+
				casecounts5[0][2][2][0][0]+casecounts5[2][1][2][0][0]+controlcounts5[0][1][0][0][0]+
				controlcounts5[2][0][0][0][0]+controlcounts5[2][2][0][0][0]+controlcounts5[0][0][2][0][0]+
				controlcounts5[0][2][2][0][0]+controlcounts5[2][1][2][0][0];

			help2= casecounts5[0][0][0][0][0]+casecounts5[0][2][0][0][0]+casecounts5[2][1][0][0][0]+casecounts5[0][1][2][0][0]+
				casecounts5[2][0][2][0][0]+casecounts5[2][2][2][0][0]+controlcounts5[0][0][0][0][0]+
				controlcounts5[0][2][0][0][0]+controlcounts5[2][1][0][0][0]+controlcounts5[0][1][2][0][0]+
				controlcounts5[2][0][2][0][0]+controlcounts5[2][2][2][0][0];

			help3=total-help1-help2;

			if(help1==total || help2==total || help3==total){x1x2Dx3=0;}
		}

		if(x1x2Dx3D) //update if always equal
		{
			help1=casecounts5[0][1][0][0][0]+casecounts5[2][0][0][0][0]+casecounts5[2][2][0][0][0]+
			casecounts5[0][0][1][0][0]+casecounts5[0][2][1][0][0]+casecounts5[2][1][1][0][0]+
			casecounts5[0][1][2][0][0]+casecounts5[2][0][2][0][0]+casecounts5[2][2][2][0][0]+
			controlcounts5[0][1][0][0][0]+controlcounts5[2][0][0][0][0]+controlcounts5[2][2][0][0][0]+
			controlcounts5[0][0][1][0][0]+controlcounts5[0][2][1][0][0]+controlcounts5[2][1][1][0][0]+
			controlcounts5[0][1][2][0][0]+controlcounts5[2][0][2][0][0]+controlcounts5[2][2][2][0][0];

			help2=casecounts5[0][0][0][0][0]+casecounts5[0][2][0][0][0]+casecounts5[2][1][0][0][0]+
			casecounts5[0][1][1][0][0]+casecounts5[2][0][1][0][0]+casecounts5[2][2][1][0][0]+
			casecounts5[0][0][2][0][0]+casecounts5[0][2][2][0][0]+casecounts5[2][1][2][0][0]+
			controlcounts5[0][0][0][0][0]+controlcounts5[0][2][0][0][0]+controlcounts5[2][1][0][0][0]+
			controlcounts5[0][1][1][0][0]+controlcounts5[2][0][1][0][0]+controlcounts5[2][2][1][0][0]+
			controlcounts5[0][0][2][0][0]+controlcounts5[0][2][2][0][0]+controlcounts5[2][1][2][0][0];

			help3=total-help1-help2;

			if(help1==total || help2==total || help3==total){x1x2Dx3D=0;}
		}

		if(x1Dx2x3) //update if always equal
		{
			help1= casecounts5[0][2][0][0][0]+casecounts5[1][0][0][0][0]+casecounts5[2][2][0][0][0]+casecounts5[0][0][2][0][0]+
				casecounts5[1][2][2][0][0]+casecounts5[2][0][2][0][0]+controlcounts5[0][2][0][0][0]+
				controlcounts5[1][0][0][0][0]+controlcounts5[2][2][0][0][0]+controlcounts5[0][0][2][0][0]+
				controlcounts5[1][2][2][0][0]+controlcounts5[2][0][2][0][0];

			help2= casecounts5[0][0][0][0][0]+casecounts5[1][2][0][0][0]+casecounts5[2][0][0][0][0]+casecounts5[0][2][2][0][0]+
				casecounts5[1][0][2][0][0]+casecounts5[2][2][2][0][0]+controlcounts5[0][0][0][0][0]+
				controlcounts5[1][2][0][0][0]+controlcounts5[2][0][0][0][0]+controlcounts5[0][2][2][0][0]+
				controlcounts5[1][0][2][0][0]+controlcounts5[2][2][2][0][0];

			help3=total-help1-help2;

			if(help1==total || help2==total || help3==total){x1Dx2x3=0;}
		}

		if(x1Dx2x3D) //update if always equal
		{
			help1=casecounts5[0][2][0][0][0]+casecounts5[1][0][0][0][0]+casecounts5[2][2][0][0][0]+
				casecounts5[0][0][1][0][0]+casecounts5[1][2][1][0][0]+casecounts5[2][0][1][0][0]+
				casecounts5[0][2][2][0][0]+casecounts5[1][0][2][0][0]+casecounts5[2][2][2][0][0]+
				controlcounts5[0][2][0][0][0]+controlcounts5[1][0][0][0][0]+controlcounts5[2][2][0][0][0]+
				controlcounts5[0][0][1][0][0]+controlcounts5[1][2][1][0][0]+controlcounts5[2][0][1][0][0]+
				controlcounts5[0][2][2][0][0]+controlcounts5[1][0][2][0][0]+controlcounts5[2][2][2][0][0];

			help2=casecounts5[0][0][0][0][0]+casecounts5[1][2][0][0][0]+casecounts5[2][0][0][0][0]+
				casecounts5[0][2][1][0][0]+casecounts5[1][0][1][0][0]+casecounts5[2][2][1][0][0]+
				casecounts5[0][0][2][0][0]+casecounts5[1][2][2][0][0]+casecounts5[2][0][2][0][0]+
				controlcounts5[0][0][0][0][0]+controlcounts5[1][2][0][0][0]+controlcounts5[2][0][0][0][0]+
				controlcounts5[0][2][1][0][0]+controlcounts5[1][0][1][0][0]+controlcounts5[2][2][1][0][0]+
				controlcounts5[0][0][2][0][0]+controlcounts5[1][2][2][0][0]+controlcounts5[2][0][2][0][0];

			help3=total-help1-help2;

			if(help1==total || help2==total || help3==total){x1Dx2x3D=0;}
		}

		if(x1Dx2Dx3) //update if always equal
		{
			help1=casecounts5[0][1][0][0][0]+casecounts5[1][0][0][0][0]+casecounts5[1][2][0][0][0]+
				casecounts5[2][1][0][0][0]+casecounts5[0][0][2][0][0]+casecounts5[0][2][2][0][0]+
				casecounts5[1][1][2][0][0]+casecounts5[2][0][2][0][0]+casecounts5[2][2][2][0][0]+
				controlcounts5[0][1][0][0][0]+controlcounts5[1][0][0][0][0]+controlcounts5[1][2][0][0][0]+
				controlcounts5[2][1][0][0][0]+controlcounts5[0][0][2][0][0]+controlcounts5[0][2][2][0][0]+
				controlcounts5[1][1][2][0][0]+controlcounts5[2][0][2][0][0]+controlcounts5[2][2][2][0][0];

			help2=casecounts5[0][0][0][0][0]+casecounts5[0][2][0][0][0]+casecounts5[1][1][0][0][0]+
				casecounts5[2][0][0][0][0]+casecounts5[2][2][0][0][0]+casecounts5[0][1][2][0][0]+
				casecounts5[1][0][2][0][0]+casecounts5[1][2][2][0][0]+casecounts5[2][1][2][0][0]+
				controlcounts5[0][0][0][0][0]+controlcounts5[0][2][0][0][0]+controlcounts5[1][1][0][0][0]+
				controlcounts5[2][0][0][0][0]+controlcounts5[2][2][0][0][0]+controlcounts5[0][1][2][0][0]+
				controlcounts5[1][0][2][0][0]+controlcounts5[1][2][2][0][0]+controlcounts5[2][1][2][0][0];

			help3=total-help1-help2;

			if(help1==total || help2==total || help3==total){x1Dx2Dx3=0;}
		}

		if(x1Dx2Dx3D) //update if always equal
		{
			help1= casecounts5[0][1][0][0][0]+controlcounts5[0][1][0][0][0]+casecounts5[1][0][0][0][0]+
				controlcounts5[1][0][0][0][0]+casecounts5[1][2][0][0][0]+controlcounts5[1][2][0][0][0]+
				casecounts5[2][1][0][0][0]+controlcounts5[2][1][0][0][0]+casecounts5[0][0][1][0][0]+
				controlcounts5[0][0][1][0][0]+casecounts5[0][2][1][0][0]+controlcounts5[0][2][1][0][0]+
				casecounts5[1][1][1][0][0]+controlcounts5[1][1][1][0][0]+casecounts5[2][0][1][0][0]+

				controlcounts5[2][0][1][0][0]+casecounts5[2][2][1][0][0]+controlcounts5[2][2][1][0][0]+
				casecounts5[0][1][2][0][0]+controlcounts5[0][1][2][0][0]+casecounts5[1][0][2][0][0]+
				controlcounts5[1][0][2][0][0]+casecounts5[1][2][2][0][0]+controlcounts5[1][2][2][0][0]+
				casecounts5[2][1][2][0][0]+controlcounts5[2][1][2][0][0];

			help2=total-help1;


			if(help1==total || help2==total){x1Dx2Dx3D=0;}
		}
	} //end if snps>2

	//SET INDICATOR
	g1=x1+x1D+x1x2+x1x2D+x1Dx2+x1Dx2D+
	x1x3+x1x3D+x1Dx3+x1Dx3D+
	x1x2x3+x1x2x3D+x1x2Dx3+x1x2Dx3D+x1Dx2x3+x1Dx2x3D+x1Dx2Dx3+x1Dx2Dx3D;

	g2=x2+x2D+x1x2+x1x2D+x1Dx2+x1Dx2D+
	x2x3+x2x3D+x2Dx3+x2Dx3D+
	x1x2x3+x1x2x3D+x1x2Dx3+x1x2Dx3D+x1Dx2x3+x1Dx2x3D+x1Dx2Dx3+x1Dx2Dx3D;

	g3=x3+x3D+x1x3+x1x3D+x1Dx3+x1Dx3D+x2x3+x2x3D+x2Dx3+x2Dx3D+
	x1x2x3+x1x2x3D+x1x2Dx3+x1x2Dx3D+x1Dx2x3+x1Dx2x3D+x1Dx2Dx3+x1Dx2Dx3D;



	if(g3){snps=3;}
	else if (g2){snps=2;}
	else if (g1){snps=1;}
	//printf("snps %d\n",snps);

	help1=0;help2=0;help3=0;
	if(x1==1){indicator[1]=1;result.in[1]=1;help1++;}
	if(x1D==1){indicator[2]=1;result.in[2]=1;help1++;}

	if(snps>1)
	{
		if(x2==1){indicator[3]=1;result.in[3]=1;}
		if(x2D==1){indicator[4]=1;result.in[4]=1;}
		if(x1x2==1){indicator[5]=1;result.in[5]=1;}
		if(x1x2D==1){indicator[6]=1;result.in[6]=1;}
		if(x1Dx2==1){indicator[7]=1;result.in[7]=1;}
		if(x1Dx2D==1){indicator[8]=1;result.in[8]=1;}
	}
	if(snps>2)
	{
		if(x3==1){indicator[9]=1;result.in[9]=1;}
		if(x3D==1){indicator[10]=1;result.in[10]=1;}
		if(x1x3==1){indicator[11]=1;result.in[11]=1;}
		if(x1x3D==1){indicator[12]=1;result.in[12]=1;}
		if(x1Dx3==1){indicator[13]=1;result.in[13]=1;}
		if(x1Dx3D==1){indicator[14]=1;result.in[14]=1;}
		if(x2x3==1){indicator[15]=1;result.in[15]=1;}
		if(x2x3D==1){indicator[16]=1;result.in[16]=1;}
		if(x2Dx3==1){indicator[17]=1;result.in[17]=1;}
		if(x2Dx3D==1){indicator[18]=1;result.in[18]=1;}
		if(x1x2x3==1){indicator[19]=1;result.in[19]=1;}
		if(x1x2x3D==1){indicator[20]=1;result.in[20]=1;}
		if(x1x2Dx3==1){indicator[21]=1;result.in[21]=1;}
		if(x1x2Dx3D==1){indicator[22]=1;result.in[22]=1;}
		if(x1Dx2x3==1){indicator[23]=1;result.in[23]=1;}
		if(x1Dx2x3D==1){indicator[24]=1;result.in[24]=1;}
		if(x1Dx2Dx3==1){indicator[25]=1;result.in[25]=1;}
		if(x1Dx2Dx3D==1){indicator[26]=1;result.in[26]=1;}
	}

	for(i=28;i<28+maxIndexCov;i++)
	{
		if(cov[i-28]==1){indicator[i]=1;withcov=1;result.in[i]=1;}
	}

	//determine df
	numberit[0]=0;number2indi[0]=0;
	k=0;


	for(i=1;i<28+maxIndexCov;i++)
	{
		f+=indicator[i];//df
		if(indicator[i]==1)
		{
			k++;
			numberit[i]=k;
			number2indi[k]=i;
		}
		else
		{
			numberit[i]=-1;
		}
	}

	//moved here form below in IS 682
	for(i=0;i< f+1;i++)
	{
		result.b[number2indi[i]]=-1;  // betas
        //result.betaNew[number2indi[i]] = -1;
		if(nMc==0)
		  {
			result.betaNew_se[number2indi[i]] = -1;
			result.oddsRatio[number2indi[i]]=-1; // odds ratios
			result.lcloddsRatio[number2indi[i]] = -1;
			result.rcloddsRatio[number2indi[i]] = -1;
		  }
	}

	if(xtype==0 && help1==0 && withcov==0 && sexcov == 0 && N==1 && quick)
	{
	    result.df=f;

		ncased=counts.AA_Ca+counts.AB_Ca+counts.BB_Ca;
		ncontrold=counts.AA_Co+counts.AB_Co+counts.BB_Co;
		//printf("ncased %f ncontrold %f\n",ncased,ncontrold);

		if(ncased+ncontrold>0)
		{
			result.b[0]=log(ncased/ncontrold);
		}
		else
		{
			result.b[0]=0;
		}

		logNew=0;
		double phelp=0;
		//AABB
		exponent=result.b[0]; phelp=1/(1+exp(-exponent));
		logNew=ncased*log(phelp)+ncontrold*log(1-phelp);

		result.sc=logNew;

		//NEW
		result.b[1]=0;
		if(nMc==0)
		  {

		    result.oddsRatio[1]=1;result.lcloddsRatio[1]=-1;result.rcloddsRatio[1]=-1;
		    result.betaNew_se[1]=-1;

		    // result.oddsRatio[1]=1;result.lcloddsRatio[1]=-1;result.rcloddsRatio[1]=-1;
			//			result.betaNew_se[1]=0;
		  }

		return result;
	} //end quicker tests

	//get n, Y and X
	n=0;ncase=0;ncontrol=0;


    struct PPLLOCATION guy;
	for(int kMod=0;kMod<npplqc;kMod++)
	{
		guy = PplLocations[kMod];
		k = PPLMap[kMod];
		complete=1;

		if(alt==0 && Yhelp[k]< -0.1){continue;}

		//person in?
		if(alt==1)
		{
			if(xtype==0)
			{
				if(person[k].aff[thread]==2 && person[k].qcin==1)
				{
					Y[n]=1;
					Yhelp[k]=1;
				}
				else if(person[k].aff[thread]==1 && person[k].qcin==1)
				{
					Y[n]=0;
					Yhelp[k]=0;
				}
				else {Yhelp[k]=-1;
				continue;}
			}
			else if(xtype==1)
			{
				if(person[k].aff[thread]==2 && person[k].sex==1 && person[k].qcin==1)
				{
					Y[n]=1;
					Yhelp[k]=1;
				}
				else if(person[k].aff[thread]==1 && person[k].sex==1 && person[k].qcin==1)
				{
					Y[n]=0;
					Yhelp[k]=0;
				}
				else
				{
					Yhelp[k]=-1;
					continue;
				}
			}
			else if(xtype==2)
			{
				if(person[k].aff[thread]==2 && person[k].sex==2 && person[k].qcin==1)
				{
					Y[n]=1;
					Yhelp[k]=1;
				}
				else if(person[k].aff[thread]==1 && person[k].sex==2 && person[k].qcin==1)
				{
					Y[n]=0;
					Yhelp[k]=0;
				}
				else {Yhelp[k]=-1;
				continue;}
			}
			else if(xtype==3)
			{
				if((person[k].sex==0) || (person[k].sex==1 && getbit64(BinSNPs[a][guy.nr][2],guy.pos))){Yhelp[k]=-1;continue;}
				if(person[k].aff[thread]==2 && person[k].qcin==1)
				{
					Y[n]=1;
					Yhelp[k]=1;
				}
				else if(person[k].aff[thread]==1 && person[k].qcin==1)
				{
					Y[n]=0;
					Yhelp[k]=0;
				}
				else
				{
					Yhelp[k]=-1;
					continue;
				}
			}
		}
		else
		{
			Y[n]=Yhelp[k];
		}


//get main x values
		if(g1)
		{
			if(dosage==1){
				if(singleMarkerTest == 3){
					currentval[0]=(2*genoWeights[a][k][0]+genoWeights[a][k][1])-1;
				}
			}else{

				//if(M[a].d2[k].d3==1)
				if(getbit64(BinSNPs[a][guy.nr][1],guy.pos))
				{
					currentval[0]=1;currentval[1]=-0.5;
					if(singleMarkerTest == 7){currentval[0]=0;}

				}
				//else if(M[a].d2[k].d3==2)
				else if(getbit64(BinSNPs[a][guy.nr][2],guy.pos))
				{
					currentval[0]=0;currentval[1]=0.5;
					if(singleMarkerTest == 6 || singleMarkerTest == 7){currentval[0]=1;}

				}
				//else if(M[a].d2[k].d3==3)
				else if(getbit64(BinSNPs[a][guy.nr][3],guy.pos))
				{
					currentval[0]=-1;currentval[1]=-0.5;
					if(singleMarkerTest == 5 || singleMarkerTest == 6 || singleMarkerTest == 7){currentval[0]=0;}

				}
				else
				{
					Yhelp[k]=-1;
					continue;
				}
				if(xtype==3 /*&& person[k].sex==1*/)
				{
					//currentval[0]+=0.5;currentval[1]+=0.5;
				}

				if(covNum>=0)
				  {
				   if(covNum>0){currentval[0]=person[k].cov[covNum-1];}
				   else {currentval[0]=person[k].sex;}
				  }
			}
		}

		if(g2)
		{
			//if(M[b].d2[k].d3==1)
			if(getbit64(BinSNPs[b][guy.nr][1],guy.pos))
			{
				currentval[2]=1;currentval[3]=-0.5;
			}
			//else if(M[b].d2[k].d3==2)
			else if(getbit64(BinSNPs[b][guy.nr][2],guy.pos))
			{
				currentval[2]=0;currentval[3]=0.5;
			}
			//else if(M[b].d2[k].d3==3)
			else if(getbit64(BinSNPs[b][guy.nr][3],guy.pos))
			{
				currentval[2]=-1;currentval[3]=-0.5;
			}
			else
			{
				Yhelp[k]=-1;
				continue;
			}
		}

		if(g3)
		{
			//if(M[c].d2[k].d3==1)
			if(getbit64(BinSNPs[c][guy.nr][1],guy.pos))
			{
				currentval[4]=1;currentval[5]=-0.5;
			}

			//else if(M[c].d2[k].d3==2)
			else if(getbit64(BinSNPs[c][guy.nr][2],guy.pos))
			{
				currentval[4]=0;currentval[5]=0.5;
			}

			//else if(M[c].d2[k].d3==3)
			else if(getbit64(BinSNPs[c][guy.nr][3],guy.pos))
			{

				currentval[4]=-1;currentval[5]=-0.5;
			}
			else
			{
				Yhelp[k]=-1;
				continue;
			}
		}


		//update X
		for(j=0;j<=f;j++)
		{
			jj=number2indi[j];

			if(jj==0) //x0
			{
				X[n][j]=1;
			}
			else if(jj==1 && liabilityCut == 0) //x1, no liability model
			{

				X[n][j]=currentval[0];
			}
			else if(jj==1 && liabilityCut != 0) //x1, liability model
			{
				if(person[k].load>=liabilityCut && liabilityCut != -1)
				{
					X[n][j]=1;
				}
				else
				{
					X[n][j]=0;
				}
				if(liabilityCut == -1)
				{
					X[n][j]=pow(person[k].load/(((double)maxIndexCov)*2),3);
				}
			}
			else if(jj==2) //x1D
			{
				X[n][j]=currentval[1];
			}
			else if(jj==3) //x2
			{
				X[n][j]=currentval[2];
			}
			else if(jj==4) //x2D
			{
				X[n][j]=currentval[3];
			}
			else if(jj==5) //x1x2
			{
				X[n][j]=currentval[0]*currentval[2];
			}
			else if(jj==6) //x1x2D
			{
				X[n][j]=currentval[0]*currentval[3];
			}
			else if(jj==7) //x1Dx2
			{
				X[n][j]=currentval[1]*currentval[2];
			}
			else if(jj==8) //x1Dx2D
			{
				X[n][j]=currentval[1]*currentval[3];
			}
			else if(jj==9) //x3
			{
				X[n][j]=currentval[4];
			}
			else if(jj==10) //x3D
			{
				X[n][j]=currentval[5];
			}
			else if(jj==11) //x1x3
			{
				X[n][j]=currentval[0]*currentval[4];
			}
			else if(jj==12) //x1x3D
			{
				X[n][j]=currentval[0]*currentval[5];
			}
			else if(jj==13) //x1Dx3
			{
				X[n][j]=currentval[1]*currentval[4];
			}
			else if(jj==14) //x1Dx3D
			{
				X[n][j]=currentval[1]*currentval[5];
			}
			else if(jj==15) //x2x3
			{
				X[n][j]=currentval[2]*currentval[4];
			}
			else if(jj==16) //x2x3D
			{
				X[n][j]=currentval[2]*currentval[5];

			}
			else if(jj==17) //x2Dx3
			{
				X[n][j]=currentval[3]*currentval[4];
			}
			else if(jj==18) //x2Dx3D
			{
				X[n][j]=currentval[3]*currentval[5];
			}
			else if(jj==19) //x1x2x3
			{
				X[n][j]=currentval[0]*currentval[2]*currentval[4];
			}
			else if(jj==20) //x1x2x3D
			{
				X[n][j]=currentval[0]*currentval[2]*currentval[5];
			}
			else if(jj==21) //x1x2Dx3
			{
				X[n][j]=currentval[0]*currentval[3]*currentval[4];
			}
			else if(jj==22) //x1x2Dx3D
			{
				X[n][j]=currentval[0]*currentval[3]*currentval[5];
			}
			else if(jj==23) //x1Dx2x3
			{
				X[n][j]=currentval[1]*currentval[2]*currentval[4];
			}
			else if(jj==24) //x1Dx2x3D
			{
				X[n][j]=currentval[1]*currentval[2]*currentval[5];
			}
			else if(jj==25) //x1Dx2Dx3
			{
				X[n][j]=currentval[1]*currentval[3]*currentval[4];
			}
			else if(jj==26) //x1Dx2Dx3D
			{
				X[n][j]=currentval[1]*currentval[3]*currentval[5];
			}
			else if (jj==27)
			{
				if(person[k].sex==1){X[n][j]=1;}
				else if(person[k].sex==2){X[n][j]=0;}
				else /* if(alt==1) */ {Yhelp[k]=-1;complete=0;break;}
				//else{complete=0;break;}  never gets there
			}
			else if(jj>27)
			{
				if(person[k].covin[jj-28]==0) //NEU
				{
					Yhelp[k]=-1; complete=0;break;
				}
				else{X[n][j]=person[k].cov[jj-28];}
			}
		}  //end j-loop (paramters)

		if(complete)
		{
			n++;
			if(person[k].aff[thread]==2){ncase++;}
			else if(person[k].aff[thread]==1){ncontrol++;}
		}
	} //end k-loop (individuals)


	//MAXIMIZATION
	logNew=0;
	logOld=1;
	it=0;
	for(i=0;i<28+maxIndexCov;i++){betaNew[i]=-1;}
	for(j=0;j<=f;j++)
	{
		betaNew[number2indi[j]]=0;
	}
	betaNew[0]=0;

	int goOn=0;

	if(f>=0 && ncase>0 && ncontrol>0)
	{
	  while((fabs(logNew-logOld)>eps && it <= maxit) || (covariancematrix && it <=15) )
		{
			it++;
			logOld=logNew;
			for(j=0;j<=f;j++)
			{
				betaOld[number2indi[j]]=betaNew[number2indi[j]];
			}
			//determine p

			if(it==1)
			{
				for (i=0;i<n;i++)
				{
					exponent=0;
					for(j=0;j<=f;j++)
					{
						exponent+=betaOld[number2indi[j]]*X[i][j];
						// get Xt
						Xt[j][i]=X[i][j];
					}
					p[i]=1/(1+exp(-exponent));
					//					printf("%g ",p[i]);
				}
			} //it==1

			if(1)
			{
			//determine Xmod

				for (i=0;i<n;i++)
				{
					for(j=0;j<=f;j++)
					{
					Xmod[i][j]=X[i][j]*p[i]*(1-p[i]);
					}
				}
				dummy=MatrixMultSym(Xt,f+1,n,Xmod,n,f+1,A);

			}


			// Invertieren von A

			// Dwyer-Algorithmus
			if(maxIndexCov==0){dummy=DwyerInv(f+1, A, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, Ainv);}
			dummy=0;

			if (dummy == 0)
			{
				for(i=0;i<=f;i++)
				{
					S[i]=0;
					for(j=0;j<=f;j++)
					{
						VNN[i][j]=0;
						Sinv[i][j]=0;
					}
				}


				dummy=mysvd(f+1, f+1, &A, &VNN, &S, f+1, f+1);


				for(j=0;j<=f;j++)
				{
					for (i=0;i<=f;i++)
					{
						UNNT[j][i]=A[i][j];
					}
				}
				newf=-1;

				for (i=0;i<=f;i++)
				  {
				    if(fabs(S[i])>1.e-6){Sinv[i][i]=1/S[i];newf++;}
				    else
				      {
					   dfUnchanged=0;goOn=0;
				      }
				    if(fabs(S[i])<1.e-3 && covariancematrix){dfUnchanged=0;}
				  }

				dummy=MatrixMultMod(VNN, f+1, f+1,Sinv, f+1, f+1, A0);
				dummy=MatrixMultMod(A0, f+1, f+1,UNNT, f+1, f+1, Ainv);
			}
			else
			{
				newf=f;
			}

			//compute A^-1 * Xt
			dummy=MatrixMultMod(Ainv, f+1, f+1,Xt, f+1, n, AinvXt);

			//compute Yminusp

			for(i=0;i<n;i++){Yminusp[i][0]=Y[i]-p[i];}

			//compute AinvXt *(y-p),update betas
			dummy=MatrixMultMod(AinvXt, f+1, n, Yminusp, n, 1, newbeta);

			for(i=0;i<=f;i++)
			{
				betaNew[number2indi[i]]=betaOld[number2indi[i]]+newbeta[i][0];
			}



			//compute logNew
			logNew=0;
			for(i=0;i<n;i++)
			{
				exponent=0;
				for(j=0;j<=f;j++)
				{
					exponent+=betaNew[number2indi[j]]*X[i][j];
				}
				p[i]=1/(1+exp(-exponent));
				if(Y[i]==1)
				{
					logNew+=log(p[i]);
				}
				else if(Y[i]==0)
				{
					logNew+=log(1-p[i]);
				}
			}
		} //end while
	} //f >0
	else
	{
		f=0;logNew=0;
	}

	for(i=0;i<ncase+ncontrol;i++)
	{
		rss += (Yminusp[i][0]*Yminusp[i][0]);
	}



	result.df=newf;

	//cout << "df " << result.df << "\n";
	//cout << "df " << result.df << "\n";
	//cout << "f " << f << "\n";
	//cout << "newf " << newf << "\n";

	//moved up in IS 682
	/*for(i=0;i< f+1;i++)
	{
		result.b[number2indi[i]]=-1;  // betas

		if(nMc==0)
		  {
			result.betaNew_se[number2indi[i]] = -1;
			result.oddsRatio[number2indi[i]]=-1; // odds ratios
			result.lcloddsRatio[number2indi[i]] = -1;
			result.rcloddsRatio[number2indi[i]] = -1;
		  }
	}*/


	//NEW
	if(newf==0)
	{
		result.b[1]=0;

		if(nMc==0)
		  {
			result.oddsRatio[1]=1;result.lcloddsRatio[1]=-1;result.rcloddsRatio[1]=-1;
			result.betaNew_se[1]=-1;
		  }
	}

	else if(dfUnchanged)
	 {
		for(i=0;i< f+1;i++)
		{

			result.b[number2indi[i]]=betaNew[number2indi[i]];  // betas

			if(nMc==0 && (ncase+ncontrol-newf-1) > 0)
		    {
				if((Ainv[i][i]*rss/(ncase+ncontrol-newf-1))<0)
				{
					result.betaNew_se[number2indi[i]] = -1;
					result.lcloddsRatio[number2indi[i]] = -1;
					result.rcloddsRatio[number2indi[i]] = -1;

				}
				else
				{

				        //result.betaNew_se[number2indi[i]] = 2*sqrt(Ainv[i][i]*rss/(ncase+ncontrol-newf-1));
      				        result.betaNew_se[number2indi[i]] = sqrt(Ainv[i][i]); //changed on 09.12.2013
					result.lcloddsRatio[number2indi[i]] = exp(betaNew[number2indi[i]] -1.96*result.betaNew_se[number2indi[i]]);
					result.rcloddsRatio[number2indi[i]] = exp(betaNew[number2indi[i]] +1.96*result.betaNew_se[number2indi[i]]);
				}
				result.oddsRatio[number2indi[i]]=exp(betaNew[number2indi[i]]); // odds ratios
			}
		}
	 }
	else
	  {
	    for(i=0;i< f+1;i++)
	      {

		result.b[number2indi[i]]=-2;  // betas
		result.betaNew_se[number2indi[i]] = -2;
		result.lcloddsRatio[number2indi[i]] = -2;
		result.rcloddsRatio[number2indi[i]] = -2;
		result.oddsRatio[number2indi[i]]= -2; // odds ratios
	      }
	  }




	result.bcvsex=betaNew[27];  //covariate sex



	for(i=0;i<maxIndexCov;i++)
	{
		result.bcv[i]=betaNew[i+27]; //covariate parameters
	}
	result.sc=logNew;



	// Calculate OR and SE
	if (alt == 1 && liabilityCut != 0)
	{

		liabilityOR[0][0] = 0;  // Cases, <= load
		liabilityOR[0][1] = 0;// Cases, > load
		liabilityOR[1][0] = 0;// Controls, <= load
		liabilityOR[1][1] = 0; // Controls, > load
		result.oddsRatio[1] = 0; // odds ratio
		result.lcloddsRatio[1] = 0; // SE
		result.betaNew_lcl[1] = 0; // Cases > load
		result.betaNew_rcl[1] = 0; // Controls > load



		for (int iMod=0; iMod< npplqc; iMod++)
		{
			guy = PplLocations[iMod];
			int i = 0;
			i = PPLMap[iMod];



			if (person[i].aff[thread] == 2 &&  (person[i].load < liabilityCut))
			{
				liabilityOR[0][0]++;
			}
			else if (person[i].aff[thread] == 2 && (person[i].load >= liabilityCut))
			{
				liabilityOR[0][1]++;

			}
			else if (person[i].aff[thread] == 1 && (person[i].load < liabilityCut))
			{
				liabilityOR[1][0]++;
			}
			else if (person[i].aff[thread] == 1 && (person[i].load >= liabilityCut))
			{
				liabilityOR[1][1]++;
			}
		}

		//cout << liabilityOR[0][0] << ": " << liabilityOR[0][1] << ": " << liabilityOR[1][0] << ": " << liabilityOR[1][1] << "\n";
		if (liabilityOR[0][0] != 0 && liabilityOR[0][1] != 0 && liabilityOR[1][0] != 0 && liabilityOR[1][1] != 0)
		{
			result.oddsRatio[1] = (liabilityOR[0][1]*liabilityOR[1][0])/(liabilityOR[0][0]*liabilityOR[1][1]); // odds ratio
			result.lcloddsRatio[1] = 1.96 * sqrt(1/liabilityOR[0][0]+1/liabilityOR[0][1]+1/liabilityOR[1][0]+1/liabilityOR[1][1]); // SE
			result.betaNew_lcl[1] = liabilityOR[0][1]/ (liabilityOR[0][0] + liabilityOR[0][1]); // Cases > load
			result.betaNew_rcl[1] = liabilityOR[1][1]/ (liabilityOR[1][0] + liabilityOR[1][1]); // Controls > load

		}
		else
		{
			result.oddsRatio[1] = -99;
			result.lcloddsRatio[1] = -99;
			result.betaNew_lcl[1] = -99;
			result.betaNew_rcl[1] = -99;
		}
	}
	if(alt==1 && covariancematrix)  //SIGMA1
	  {
	    int k=0;

	    for(i=0;i<f+1;i++)
	      {
		for(j=i;j<f+1;j++)
		  {
		    if(dfUnchanged==1) {result.sigma1[k]=Ainv[i][j];} // sigma1 matrix as 1 dim array
		    else {result.sigma1[k]=-2;}
		    //cout << "f=" << f << " Sigma1[" << k << "]=" << Ainv[i][j]<< endl;
		    k=k+1;
		  }
	      }
	  }

	return result;
};



void resultReg(struct STATplus result, fstream &file, int haplo, int qt, int maxIndexCov)
{
	if(!haplo)
	{
		if (qt == 0)
		{
			if (result.in[0] == 1)
			{
				file << "OR0: " << result.oddsRatio[0] << "\t" << "OR_lcl0: " << result.lcloddsRatio[0]  << "\t" << "ORrcl0: "
				<< result.rcloddsRatio[0] << "\n";      //intercept parameter
			}
			if (result.in[1] == 1)
			{
				file << "OR1: " << result.oddsRatio[1] << "\t" << "OR_lcl1: " << result.lcloddsRatio[1]  << "\t" << "ORrcl1: "
				<< result.rcloddsRatio[1] << "\n";
			}
			if (result.in[2] == 1)
			{
				file << "OR1D: " << result.oddsRatio[2] << "\t" << "OR_lcl1D: " << result.lcloddsRatio[2]  << "\t" << "ORrcl1D: "
				<< result.rcloddsRatio[2] << "\n";
			}
			if (result.in[3] == 1)
			{
				file << "OR2: " << result.oddsRatio[3] << "\t" << "OR_lcl2: " << result.lcloddsRatio[3]  << "\t" << "ORrcl2: "
				<< result.rcloddsRatio[3] << "\n";
			}
			if (result.in[4] == 1)
			{
				file << "OR2D: " << result.oddsRatio[4] << "\t" << "OR_lcl2D: " << result.lcloddsRatio[4]  << "\t" << "ORrcl2D: "
				<< result.rcloddsRatio[4] << "\n";
			}
			if (result.in[5] == 1)
			{
				file << "OR12: " << result.oddsRatio[5] << "\t" << "OR_lcl12: " << result.lcloddsRatio[5]  << "\t" << "ORrcl12: "
				<< result.rcloddsRatio[5] << "\n"; //2-fold result.interaction parameters (snps 1 and 2)
			}
			if (result.in[6] == 1)
			{
				file << "OR12D: " << result.oddsRatio[6] << "\t" << "OR_lcl12D: " << result.lcloddsRatio[6]  << "\t" << "ORrcl12D: "
				<<result.rcloddsRatio[6] << "\n";
			}
			if (result.in[7] == 1)
			{
				file << "OR1D2: " << result.oddsRatio[7] << "\t" << "OR_lcl1D2: " << result.lcloddsRatio[7]  << "\t" << "ORrcl1D2: "
				<< result.rcloddsRatio[7] << "\n";
			}
			if (result.in[8] == 1)
			{
				file << "OR1D2D: " << result.oddsRatio[8] << "\t" << "OR_lcl1D2D: " << result.lcloddsRatio[8]  << "\t" << "ORrcl1D2D: "
				<< result.rcloddsRatio[8] << "\n";
			}
			if (result.in[9] == 1)
			{
				file << "OR3: " << result.oddsRatio[9] << "\t" << "OR_lcl3: " << result.lcloddsRatio[9]  << "\t" << "ORrcl3: "
				<< result.rcloddsRatio[9] << "\n";
			}
			if (result.in[10] == 1)
			{
				file << "OR3D: " << result.oddsRatio[10] << "\t" << "OR_lcl3D: " << result.lcloddsRatio[10]  << "\t" << "ORrcl3D: "
				<< result.rcloddsRatio[10] << "\n"; //snp3 dominant
			}
			if (result.in[11] == 1)
			{
				file << "OR13: " << result.oddsRatio[11] << "\t" << "OR_lcl13: " << result.lcloddsRatio[11]  << "\t" << "ORrcl13: "
				<< result.rcloddsRatio[11] << "\n";
			}
			if (result.in[12] == 1)
			{
				file << "OR13D: " << result.oddsRatio[12] << "\t" << "OR_lcl13D: " << result.lcloddsRatio[12]  << "\t" << "ORrcl13D: "
				<< result.rcloddsRatio[12] << "\n";
			}
			if (result.in[13] == 1)
			{
				file << "OR1D3: " << result.oddsRatio[13] << "\t" << "OR_lcl1D3: " << result.lcloddsRatio[13]  << "\t" << "ORrcl1D3: "
				<< result.rcloddsRatio[13] << "\n";
			}
			if (result.in[14] == 1)
			{
				file << "OR1D3D: " << result.oddsRatio[14] << "\t" << "OR_lcl1D3D: " << result.lcloddsRatio[14]  << "\t" << "ORrcl1D3D: "
				<< result.rcloddsRatio[14] << "\n";
			}
			if (result.in[15] == 1)
			{
				file << "OR23: " << result.oddsRatio[15] << "\t" << "OR_lcl23: " << result.lcloddsRatio[15]  << "\t" << "ORrcl23: "
				<< result.rcloddsRatio[15] << "\n";

			}

			if (result.in[16] == 1)
			{
				file << "OR23D: " << result.oddsRatio[16] << "\t" << "OR_lcl23D: " << result.lcloddsRatio[16]  << "\t" << "ORrcl23D: "
				<< result.rcloddsRatio[16] << "\n";
			}
			if (result.in[17] == 1)
			{
				file << "ORD3: " << result.oddsRatio[17] << "\t" << "OR_lclD3: " << result.lcloddsRatio[17]  << "\t" << "ORrclD3: "
				<< result.rcloddsRatio[17] << "\n";
			}
			if (result.in[18] == 1)
			{
				file << "OR2D3D: " << result.oddsRatio[18] << "\t" << "OR_lcl2D3D: " << result.lcloddsRatio[18]  << "\t" << "ORrcl2D3D: "
				<< result.rcloddsRatio[18] << "\n";
			}
			if (result.in[19] == 1)
			{
				file << "OR123: " << result.oddsRatio[19] << "\t" << "OR_lcl123: " << result.lcloddsRatio[19]  << "\t" << "ORrcl123: "
				<< result.rcloddsRatio[19] << "\n";
			}

			if (result.in[20] == 1)
			{
				file << "OR123D: " << result.oddsRatio[20] << "\t" << "OR_lcl23D: " << result.lcloddsRatio[20]  << "\t" << "ORrcl123D: "
				<< result.rcloddsRatio[20] << "\n";
			}

			if (result.in[21] == 1)
			{
				file << "OR12D3: " << result.oddsRatio[21] << "\t" << "OR_lcl12D3: " << result.lcloddsRatio[21]  << "\t" << "ORrcl12D3: "
				<< result.rcloddsRatio[21] << "\n";
			}

			if (result.in[22] == 1)
			{
				file << "OR12D3D: " << result.oddsRatio[22] << "\t" << "OR_lcl12D3D: " << result.lcloddsRatio[22]  << "\t" << "ORrcl12D3D: "
				<< result.rcloddsRatio[22] << "\n";
			}
			if (result.in[23] == 1)
			{
				file << "OR1D23: " << result.oddsRatio[23] << "\t" << "OR_lcl1D23: " << result.lcloddsRatio[23]  << "\t" << "ORrcl1D23: "
				<< result.rcloddsRatio[23] << "\n";
			}
			if (result.in[24] == 1)
			{
				file << "OR1D23D: " << result.oddsRatio[24] << "\t" << "OR_lcl1D23D: " << result.lcloddsRatio[24]  << "\t" << "ORrcl1D23D: "
				<< result.rcloddsRatio[24] << "\n";
			}
			if (result.in[25] == 1)
			{
				file << "OR1D2D3: " << result.oddsRatio[25] << "\t" << "OR_lcl1D2D3: " << result.lcloddsRatio[25]  << "\t" << "ORrcl1D2D3: "
				<< result.rcloddsRatio[25] << "\n";
			}
			if (result.in[26] == 1)
			{
				file << "OR1D2D3D: " << result.oddsRatio[25] << "\t" << "OR_lcl1D2D3D: " << result.lcloddsRatio[25]  << "\t"
				<< "ORrcl1D2D3D: " << result.rcloddsRatio[25] << "\n";
			}
			if (result.in[27] == 1)
			{
				file << "bcovariatesex: " << result.bcvsex << "\n";
			}
		}
		else
		{
			if (result.in[0] == 1)
			{
				file << "b0: " << result.b[0] << "\t" << "se0: " << result.betaNew_se[0]  << "\t" << "lcl0: " << result.betaNew_lcl[0] << "\t"
				<< "rcl0: " << result.betaNew_rcl[0] << "\n";      //intercept parameter
			}
			if (result.in[1] == 1)
			{
				file << "b1: " << result.b[1] << "\t" << "se1: " << result.betaNew_se[1]  << "\t" << "lcl1: " << result.betaNew_lcl[1] << "\t"
				<< "rcl1: " << result.betaNew_rcl[1] << "\n";      //snp1 allele
			}
			if (result.in[2] == 1)
			{
				file << "b1D: " << result.b[2] << "\t" << "se1D: " << result.betaNew_se[2]  << "\t" << "lcl1D: " << result.betaNew_lcl[2] << "\t"
				<< "rcl1D: " << result.betaNew_rcl[2] << "\n";     //snp1 dom
			}
			if (result.in[3] == 1)
			{
				file << "b2: " << result.b[3] << "\t" << "se2: " << result.betaNew_se[3]  << "\t" << "lcl2: " << result.betaNew_lcl[3] << "\t"
				<< "rcl2: " << result.betaNew_rcl[3] << "\n";     //snp2 allele
			}
			if (result.in[4] == 1)
			{
				file << "b2D: " << result.b[4] << "\t" << "se2D: " << result.betaNew_se[4]  << "\t" << "lcl2D: " << result.betaNew_lcl[4] << "\t"
				<< "rcl2D: " << result.betaNew_rcl[4] << "\n";     //snp2 dom
			}
			if (result.in[5] == 1)
			{
				file  << "b12: " << result.b[5] << "\t" << "se12: " << result.betaNew_se[5]  << "\t" << "lcl12: " << result.betaNew_lcl[5]
				<< "\t" << "rcl12: " << result.betaNew_rcl[5] << "\n";    //2-fold result.interaction parameters (snps 1 and 2)
			}
			if (result.in[6] == 1)
			{
				file << "b12D: " << result.b[6] << "\t" << "se12D: " << result.betaNew_se[6]  << "\t" << "lcl12D: " << result.betaNew_lcl[6]
				<< "\t" << "rcl12D: " << result.betaNew_rcl[6] << "\n";
			}
			if (result.in[7] == 1)
			{
				file << "b1D2: " << result.b[7] << "\t" << "se1D2: " << result.betaNew_se[7]  << "\t" << "lcl1D2: " << result.betaNew_lcl[7]
				<< "\t" << "rcl1D2: " << result.betaNew_rcl[7] << "\n";
			}
			if (result.in[8] == 1)
			{
				file << "b1D2D: " << result.b[8] << "\t" << "se1D2D: " << result.betaNew_se[8]  << "\t" << "lcl1D2D: " << result.betaNew_lcl[8]
				<< "\t" << "rcl1D2D: " << result.betaNew_rcl[8] << "\n";
			}
			if (result.in[9] == 1)
			{
				file << "b3: " << result.b[9] << "\t" << "se3: " << result.betaNew_se[9]  << "\t" << "lcl3: " << result.betaNew_lcl[9] << "\t"
				<< "rcl3: " << result.betaNew_rcl[9] << "\n";
			}
			if (result.in[10] == 1)
			{
				file << "b3D: " << result.b[10] << "\t" << "se0: " << result.betaNew_se[10]  << "\t" << "lcl0: " << result.betaNew_lcl[10]
				<< "\t" << "rcl0: " << result.betaNew_rcl[10] << "\n";     //snp3 dominant
			}
			if (result.in[11] == 1)
			{
				file << "b13: " << result.b[11] << "\t" << "se0: " << result.betaNew_se[11]  << "\t" << "lcl0: " << result.betaNew_lcl[11]
				<< "\t" << "rcl0: " << result.betaNew_rcl[11] << "\n";
			}
			if (result.in[12] == 1)
			{
				file << "b13D: " << result.b[12] << "\t" << "se13D: " << result.betaNew_se[12]  << "\t" << "lcl13D: " << result.betaNew_lcl[12]
				<< "\t" << "rcl13D: " << result.betaNew_rcl[12] << "\n";
			}
			if (result.in[13] == 1)
			{
				file << "b1D3: " << result.b[13] << "\t" << "se1D3: " << result.betaNew_se[13]  << "\t" << "lcl1D3: " << result.betaNew_lcl[13]
				<< "\t" << "rcl1D3: " << result.betaNew_rcl[13] << "\n";
			}
			if (result.in[14] == 1)
			{
				file << "b1D3D: " << result.b[14] << "\t" << "se1D3D: " << result.betaNew_se[14]  << "\t" << "lcl1D3D: "
				<< result.betaNew_lcl[14] << "\t" << "rcl1D3D: " << result.betaNew_rcl[14] << "\n";
			}
			if (result.in[15] == 1)
			{
				file << "b23: " << result.b[15] << "\t" << "se23: " << result.betaNew_se[15]  << "\t" << "lcl23: " << result.betaNew_lcl[15]
				<< "\t" << "rcl23: " << result.betaNew_rcl[15] << "\n";
			}
			if (result.in[16] == 1)
			{
				file << "b23D: " << result.b[16] << "\t" << "se23D: " << result.betaNew_se[16]  << "\t" << "lcl23D: " << result.betaNew_lcl[16]
				<< "\t" << "rcl23D: " << result.betaNew_rcl[16] << "\n";
			}
			if (result.in[17] == 1)
			{
				file << "b2D3: " << result.b[17] << "\t" << "se2D3: " << result.betaNew_se[17]  << "\t" << "lcl2D3: " << result.betaNew_lcl[17]
				<< "\t" << "rcl2D3: " << result.betaNew_rcl[17] << "\n";

				file << "ORD3: " << result.oddsRatio[17] << "\t" << "OR_lclD3: " << result.lcloddsRatio[17]  << "\t" << "ORrclD3: "
				<<result.rcloddsRatio[17] << "\n";
			}
			if (result.in[18] == 1)
			{
				file << "b2D3D: " << result.b[18] << "\t" << "se2D3D: " << result.betaNew_se[18]  << "\t" << "lcl2D3D: "
				<< result.betaNew_lcl[18] << "\t" << "rcl2D3D: " << result.betaNew_rcl[18] << "\n";
			}
			if (result.in[19] == 1)
			{
				file << "b123: " << result.b[19] << "\t" << "se123: " << result.betaNew_se[19]  << "\t" << "lcl123: " << result.betaNew_lcl[19]
				<< "\t" << "rcl123: " << result.betaNew_rcl[19] << "\n";
			}

			if (result.in[20] == 1)
			{
				file << "b123D: " << result.b[20] << "\t" << "se123D: " << result.betaNew_se[20]  << "\t" << "lcl123D: "
				<< result.betaNew_lcl[20] << "\t" << "rcl123D: " << result.betaNew_rcl[20] << "\n";
			}

			if (result.in[21] == 1)
			{
				file << "b12D3: " << result.b[21] << "\t" << "se12D3: " << result.betaNew_se[21]  << "\t" << "lcl12D3: "
				<< result.betaNew_lcl[21] << "\t" << "rcl12D3: " << result.betaNew_rcl[21] << "\n";
			}

			if (result.in[22] == 1)
			{
				file << "b12D3D: " << result.b[22] << "\t" << "se12D3D: " << result.betaNew_se[22]  << "\t" << "lcl12D3D: "
				<< result.betaNew_lcl[22] << "\t" << "rcl12D3D: " << result.betaNew_rcl[22] << "\n";
			}
			if (result.in[23] == 1)

			{

				file  << "b1D23: " << result.b[23] << "\t" << "se1D23: " << result.betaNew_se[23]  << "\t" << "lcl1D23: "
				<< result.betaNew_lcl[23] << "\t" << "rcl1D23: " << result.betaNew_rcl[23] << "\n";
			}
			if (result.in[24] == 1)
			{
				file << "b1D23D: " << result.b[24] << "\t" << "se1D23D: " << result.betaNew_se[24]  << "\t" << "lcl1D23D: "
				<< result.betaNew_lcl[24] << "\t" << "rcl1D23D: " << result.betaNew_rcl[24] << "\n";
			}
			if (result.in[25] == 1)
			{
				file << "b1D2D3: " << result.b[25] << "\t" << "se1D2D3: " << result.betaNew_se[25]  << "\t" << "lcl1D2D3: "
				<< result.betaNew_lcl[25] << "\t" << "rcl1D2D3: " << result.betaNew_rcl[25] << "\n";
			}
			if (result.in[26] == 1)
			{
				file << "b1D2D3D: " << result.b[26] << "\t" << "se1D2D3D: " << result.betaNew_se[26]  << "\t" << "lcl1D2D3D: "
				<< result.betaNew_lcl[26] << "\t" << "rcl1D2D3D: " << result.betaNew_rcl[26] << "\n";
			}
			if (result.in[27] == 1)
			{
				file << "bcovariatesex: " << result.bcvsex << "\n";
			}
		}
	} //end !haplo
	else //haplo
	{
		if (result.in[1] == 1)
		{
			file << "h1: " << result.h[0] << "\n";      //
		}
		if (result.in[2] == 1)
		{
			file << "h2: " << result.h[1] << "\n";     //
		}
		if (result.in[3] == 1)
		{
			file << "h3: " << result.h[2] << "\n";      //
		}
		if (result.in[4] == 1)
		{
			file << "h4: " << result.h[3] << "\n";     //
		}
		if (result.in[5] == 1)
		{
			file  << "h5: " << result.h[4] << "\n";
		}
		if (result.in[6] == 1)
		{
			file << "h6: " << result.h[5] << "\n";
		}
		if (result.in[7] == 1)
		{
			file << "h7: " << result.h[6] << "\n";
		}
		if (result.in[8] == 1)
		{
			file << "h8: " << result.h[7] << "\n";
		}
		if (result.in[0] == 1)
		{
				file << "b0: " << result.b[0] << "\t" << "se0: " << result.betaNew_se[0]  << "\t" << "lcl0: " << result.betaNew_lcl[0] << "\t"
				<< "rcl0: " << result.betaNew_rcl[0] << "\n";      //intercept parameter
		}
		if (result.in[1] == 1)
		{
			file << "b1: " << result.b[1] << "\t" << "se1: " << result.betaNew_se[1]  << "\t" << "lcl1: " << result.betaNew_lcl[1] << "\t"
			<< "rcl1: " << result.betaNew_rcl[1] << "\n";
		}
		if (result.in[2] == 1)
		{
			file << "b2: " << result.b[2] << "\t" << "se2: " << result.betaNew_se[2]  << "\t" << "lcl2: " << result.betaNew_lcl[2] << "\t"
			<< "rcl2: " << result.betaNew_rcl[2] << "\n";
		}
		if (result.in[3] == 1)
		{
			file << "b3: " << result.b[3] << "\t" << "se3: " << result.betaNew_se[3]  << "\t" << "lcl3: " << result.betaNew_lcl[3] << "\t"
			<< "rcl3: " << result.betaNew_rcl[3] << "\n";
		}
		if (result.in[4] == 1)
		{
			file << "b4: " << result.b[4] << "\t" << "se4: " << result.betaNew_se[4]  << "\t" << "lcl4: " << result.betaNew_lcl[4] << "\t"
			<< "rcl4: " << result.betaNew_rcl[4] << "\n";
		}
		if (result.in[5] == 1)
		{
			file  << "b5: " << result.b[5] << "\t" << "se5: " << result.betaNew_se[5]  << "\t" << "lcl5: " << result.betaNew_lcl[5]
			<< "\t" << "rcl5: " << result.betaNew_rcl[5] << "\n";
		}
		if (result.in[6] == 1)
		{
			file << "b6: " << result.b[6] << "\t" << "se6: " << result.betaNew_se[6]  << "\t" << "lcl6: " << result.betaNew_lcl[6]
			<< "\t" << "rcl6: " << result.betaNew_rcl[6] << "\n";
		}
		if (result.in[7] == 1)
		{
			file << "b7: " << result.b[7] << "\t" << "se7: " << result.betaNew_se[7]  << "\t" << "lcl7: " << result.betaNew_lcl[7]
			<< "\t" << "rcl7: " << result.betaNew_rcl[7] << "\n";
		}
		if (result.in[8] == 1)
		{
			file << "b8: " << result.b[8] << "\t" << "se8: " << result.betaNew_se[8]  << "\t" << "lcl8: " << result.betaNew_lcl[8]
			<< "\t" << "rcl8: " << result.betaNew_rcl[8] << "\n";
		}
	}

	for(int i=0;i<maxIndexCov;i++)
	{
		if (result.in[i+27] == 1)
		{
			file << "bcovariate_"<< i <<": " << result.bcv[i] << "\n";
		}

	}
}



struct HAPINFO
{
	double h[8];
	double w[10];
	double r2;
};


struct HAPINFO haplo_em(double casecounts5[3][3][3][3][3],int snps,int print, char *hapfile, int dohapfile,char *hapstring, FILE *file1)
{
	struct HAPINFO new_hapinfo;
	struct HAPINFO old_hapinfo;
	double diff=1;
	double eps=0.000001;
	int i,j,k;
	double n=0;
	int it=0;
	// FILE *file1=NULL;
	int myswitch=0;
	float r2=0;
	float f_1=0;
	float f_2=0;
	float e_11=0;
	float e_12=0;
	float e_21=0;
	float e_22=0;
	double dprime=0;
	double d[2][2];
	double dmax[2][2];

	if(snps==2)
	{
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				n+=casecounts5[i][j][0][0][0];
			}
		}


		for(i=0;i<4;i++)
		{
			old_hapinfo.h[i]=0.25;
		}
		old_hapinfo.w[0]=0.5;
		new_hapinfo.w[0]=0;

		while (diff>eps && it <100)
		{
			it++;

			new_hapinfo.h[0]=
			(2*casecounts5[0][0][0][0][0]+casecounts5[0][1][0][0][0]+casecounts5[1][0][0][0][0]
			+old_hapinfo.w[0]*casecounts5[1][1][0][0][0])/(2*n);

			new_hapinfo.h[1]=
			(2*casecounts5[0][2][0][0][0]+casecounts5[0][1][0][0][0]+casecounts5[1][2][0][0][0]
			+(1-old_hapinfo.w[0])*casecounts5[1][1][0][0][0])/(2*n);

			new_hapinfo.h[2]=
			(2*casecounts5[2][0][0][0][0]+casecounts5[2][1][0][0][0]+casecounts5[1][0][0][0][0]
			+(1-old_hapinfo.w[0])*casecounts5[1][1][0][0][0])/(2*n);

			new_hapinfo.h[3]=1-new_hapinfo.h[2]-new_hapinfo.h[1]-new_hapinfo.h[0];

			if (new_hapinfo.h[0]*new_hapinfo.h[3]+new_hapinfo.h[1]*new_hapinfo.h[2] >0)
			{
				new_hapinfo.w[0]=
				new_hapinfo.h[0]*new_hapinfo.h[3]/(new_hapinfo.h[0]*new_hapinfo.h[3]+new_hapinfo.h[1]*new_hapinfo.h[2]);
			}
			else
			{new_hapinfo.w[0]= 0;}

			diff=fabs(old_hapinfo.w[0]-new_hapinfo.w[0]);

			//new = old
			for(i=0;i<4;i++)
			{
				old_hapinfo.h[i]=new_hapinfo.h[i];
			}
			old_hapinfo.w[0]=new_hapinfo.w[0];
		}
		if(print && dohapfile)
		{
			f_1=new_hapinfo.h[0]+new_hapinfo.h[1];
			f_2=new_hapinfo.h[0]+new_hapinfo.h[2];
			e_11=f_1*f_2;e_12=f_1*(1-f_2);
			e_21=(1-f_1)*f_2;e_22=(1-f_1)*(1-f_2);


			//proxy allele
			if(new_hapinfo.h[0]>=e_11)
			{
				myswitch=0;
			}
			else
			{
				myswitch=1;
			}
            //r2
			r2=pow(new_hapinfo.h[0]*new_hapinfo.h[3]-new_hapinfo.h[1]*new_hapinfo.h[2],2);
			if(f_1*f_2*(1-f_1)*(1-f_2)>0)
			{
				r2=r2/(f_1*f_2*(1-f_1)*(1-f_2));
			}
			else if(f_1*(1-f_1)==0 && f_2*(1-f_2)==0)
			{
				r2=1;
			}
			else{r2=0;}

            //dprime
			d[0][0]=new_hapinfo.h[0]-e_11;d[0][1]=new_hapinfo.h[1]-e_12;
			d[1][0]=new_hapinfo.h[2]-e_21;d[1][1]=new_hapinfo.h[3]-e_22;

			if(d[0][0]<0)
			{
				dmax[0][0]=MIN(f_1*f_2,(1-f_1)*(1-f_2));
			}
			else
			{
				dmax[0][0]=MIN((1-f_1)*f_2,f_1*(1-f_2));
			}
			if(d[0][1]<0)
			{
				dmax[0][1]=MIN(f_1*(1-f_2),(1-f_1)*f_2);
			}
			else
			{
				dmax[0][1]=MIN((1-f_1)*(1-f_2),f_1*f_2);
			}
			if(d[1][0]<0)
			{
				dmax[1][0]=MIN((1-f_1)*f_2,f_1*(1-f_2));
			}
			else
			{
				dmax[1][0]=MIN(f_1*f_2,(1-f_1)*(1-f_2));
			}
			if(d[1][1]<0)
			{
				dmax[1][1]=MIN((1-f_1)*(1-f_2),f_1*f_2);
			}
			else
			{
				dmax[1][1]=MIN(f_1*(1-f_2),(1-f_1)*f_2);
			}

            if(dmax[0][0]>0){dprime+=e_11*abs(d[0][0])/dmax[0][0];}
            if(dmax[0][1]>0){dprime+=e_12*abs(d[0][1])/dmax[0][1];}
            if(dmax[1][0]>0){dprime+=e_21*abs(d[1][0])/dmax[1][0];}
            if(dmax[1][1]>0){dprime+=e_22*abs(d[1][1])/dmax[1][1];}

            if((f_1==0 && f_2==0) || (f_1==0 && (1-f_2)==0)	|| ((1-f_1)==0 && f_2==0)	|| ((1-f_1)==0 && (1-f_2)==0))
			{
				dprime=1;
			}
			if(dprime>1){dprime=1;}

			if(dohapfile >=2)
			{
				//file1=fopen(hapfile,"a");
				fprintf(file1,"%s %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %d\n",hapstring,f_1,1-f_1,f_2,1-f_2,new_hapinfo.h[0],new_hapinfo.h[1],new_hapinfo.h[2],new_hapinfo.h[3],dprime,r2,myswitch);
				//fclose(file1);
			}
			else if(dohapfile ==1 && r2 >= 0.5)
			{
				//file1=fopen(hapfile,"a");
				fprintf(file1,"%s %1.3f %d\n",hapstring,r2,myswitch);
				//fclose(file1);
			}
		}
	}//snp == 2
	else if(snps==3)
	{
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				for(k=0;k<3;k++)
				{
					n+=casecounts5[i][j][k][0][0];
				}
			}
		}

		//printf("n %f\n",n);

		for(i=0;i<8;i++)
		{
			old_hapinfo.h[i]=0.125;
		}
		for(i=0;i<10;i++)
		{
			old_hapinfo.w[i]=0.5;
			new_hapinfo.w[i]=0;
		}
		for(i=3;i<=6;i++)
		{
			old_hapinfo.w[i]=0.25;
			new_hapinfo.w[i]=0;
		}

		while (diff>eps && it<100)
		{
			it++;

			new_hapinfo.h[0]=
			(2*casecounts5[0][0][0][0][0]+
			casecounts5[0][0][1][0][0]+
			casecounts5[0][1][0][0][0]+
			casecounts5[1][0][0][0][0]+
			old_hapinfo.w[0]*casecounts5[1][1][0][0][0]+
			old_hapinfo.w[1]*casecounts5[1][0][1][0][0]+
			old_hapinfo.w[2]*casecounts5[0][1][1][0][0]+
			old_hapinfo.w[3]*casecounts5[1][1][1][0][0])/(2*n);

			new_hapinfo.h[1]=
			(2*casecounts5[0][0][2][0][0]+
			casecounts5[0][0][1][0][0]+
			casecounts5[1][0][2][0][0]+
			casecounts5[0][1][2][0][0]+
			old_hapinfo.w[7]*casecounts5[1][1][2][0][0]+
			(1-old_hapinfo.w[1])*casecounts5[1][0][1][0][0]+
			(1-old_hapinfo.w[2])*casecounts5[0][1][1][0][0]+
			old_hapinfo.w[4]*casecounts5[1][1][1][0][0])/(2*n);

			new_hapinfo.h[2]=
			(2*casecounts5[0][2][0][0][0]+
			casecounts5[0][2][1][0][0]+
			casecounts5[1][2][0][0][0]+
			casecounts5[0][1][0][0][0]+
			old_hapinfo.w[8]*casecounts5[1][2][1][0][0]+
			(1-old_hapinfo.w[0])*casecounts5[1][1][0][0][0]+
			(1-old_hapinfo.w[2])*casecounts5[0][1][1][0][0]+
			old_hapinfo.w[5]*casecounts5[1][1][1][0][0])/(2*n);

			new_hapinfo.h[3]=
			(2*casecounts5[0][2][2][0][0]+
			casecounts5[0][2][1][0][0]+
			casecounts5[0][1][2][0][0]+
			casecounts5[1][2][2][0][0]+
			old_hapinfo.w[2]*casecounts5[0][1][1][0][0]+
			(1-old_hapinfo.w[7])*casecounts5[1][1][2][0][0]+
			(1-old_hapinfo.w[8])*casecounts5[1][2][1][0][0]+
			old_hapinfo.w[6]*casecounts5[1][1][1][0][0])/(2*n);

			new_hapinfo.h[4]=
			(2*casecounts5[2][0][0][0][0]+
			casecounts5[2][0][1][0][0]+
			casecounts5[2][1][0][0][0]+
			casecounts5[1][0][0][0][0]+
			(1-old_hapinfo.w[0])*casecounts5[1][1][0][0][0]+
			(1-old_hapinfo.w[1])*casecounts5[1][0][1][0][0]+
			old_hapinfo.w[9]*casecounts5[2][1][1][0][0]+
			old_hapinfo.w[6]*casecounts5[1][1][1][0][0])/(2*n);

			new_hapinfo.h[5]=
			(2*casecounts5[2][0][2][0][0]+
			casecounts5[2][0][1][0][0]+
			casecounts5[2][1][2][0][0]+
			casecounts5[1][0][2][0][0]+
			old_hapinfo.w[1]*casecounts5[1][0][1][0][0]+
			(1-old_hapinfo.w[9])*casecounts5[2][1][1][0][0]+
			(1-old_hapinfo.w[7])*casecounts5[1][1][2][0][0]+
			old_hapinfo.w[5]*casecounts5[1][1][1][0][0])/(2*n);

			new_hapinfo.h[6]=
			(2*casecounts5[2][2][0][0][0]+
			casecounts5[2][2][1][0][0]+
			casecounts5[2][1][0][0][0]+
			casecounts5[1][2][0][0][0]+
			old_hapinfo.w[0]*casecounts5[1][1][0][0][0]+
			(1-old_hapinfo.w[8])*casecounts5[1][2][1][0][0]+
			(1-old_hapinfo.w[9])*casecounts5[2][1][1][0][0]+
			old_hapinfo.w[4]*casecounts5[1][1][1][0][0])/(2*n);



			new_hapinfo.h[7]=1-new_hapinfo.h[6]-new_hapinfo.h[5]-new_hapinfo.h[4]-new_hapinfo.h[3]-new_hapinfo.h[2]-new_hapinfo.h[1]-
			new_hapinfo.h[0];


			for(i=0;i<8;i++)
			{

			 if(new_hapinfo.h[i]<1/(n*100))
			   {
				new_hapinfo.h[i]=0;
			   }
			}


			if (new_hapinfo.h[0]*new_hapinfo.h[6]+new_hapinfo.h[2]*new_hapinfo.h[4] >0)
			{
				new_hapinfo.w[0]=
				new_hapinfo.h[0]*new_hapinfo.h[6]/(new_hapinfo.h[0]*new_hapinfo.h[6]+new_hapinfo.h[2]*new_hapinfo.h[4]);
			}
			else
			{new_hapinfo.w[0]= 0;}

			if (new_hapinfo.h[0]*new_hapinfo.h[5]+new_hapinfo.h[1]*new_hapinfo.h[4] >0)
			{
				new_hapinfo.w[1]=
				new_hapinfo.h[0]*new_hapinfo.h[5]/(new_hapinfo.h[0]*new_hapinfo.h[5]+new_hapinfo.h[1]*new_hapinfo.h[4]);
			}
			else
			{new_hapinfo.w[1]= 0;}

			if (new_hapinfo.h[0]*new_hapinfo.h[3]+new_hapinfo.h[2]*new_hapinfo.h[1] >0)
			{
				new_hapinfo.w[2]=
				new_hapinfo.h[0]*new_hapinfo.h[3]/(new_hapinfo.h[0]*new_hapinfo.h[3]+new_hapinfo.h[2]*new_hapinfo.h[1]);
			}
			else
			{new_hapinfo.w[2]= 0;}

			if (new_hapinfo.h[0]*new_hapinfo.h[7]+new_hapinfo.h[1]*new_hapinfo.h[6] +new_hapinfo.h[2]*new_hapinfo.h[5]+
			new_hapinfo.h[3]*new_hapinfo.h[4]>0)
			{
				new_hapinfo.w[3]=
				new_hapinfo.h[0]*new_hapinfo.h[7]/(new_hapinfo.h[0]*new_hapinfo.h[7]+new_hapinfo.h[1]*new_hapinfo.h[6]+
				new_hapinfo.h[2]*new_hapinfo.h[5]+new_hapinfo.h[3]*new_hapinfo.h[4]);
				new_hapinfo.w[4]=
				new_hapinfo.h[1]*new_hapinfo.h[6]/(new_hapinfo.h[0]*new_hapinfo.h[7]+new_hapinfo.h[1]*new_hapinfo.h[6]+
				new_hapinfo.h[2]*new_hapinfo.h[5]+new_hapinfo.h[3]*new_hapinfo.h[4]);
				new_hapinfo.w[5]=
				new_hapinfo.h[2]*new_hapinfo.h[5]/(new_hapinfo.h[0]*new_hapinfo.h[7]+new_hapinfo.h[1]*new_hapinfo.h[6]+
				new_hapinfo.h[2]*new_hapinfo.h[5]+new_hapinfo.h[3]*new_hapinfo.h[4]);
			}
			else
			{
				new_hapinfo.w[3]= 0;
				new_hapinfo.w[4]= 0;
				new_hapinfo.w[5]= 0;
			}
			new_hapinfo.w[6]=1-new_hapinfo.w[5]-new_hapinfo.w[4]-new_hapinfo.w[3];

			if (new_hapinfo.h[1]*new_hapinfo.h[7]+new_hapinfo.h[3]*new_hapinfo.h[5] >0)
			{
				new_hapinfo.w[7]=
				new_hapinfo.h[1]*new_hapinfo.h[7]/(new_hapinfo.h[1]*new_hapinfo.h[7]+new_hapinfo.h[3]*new_hapinfo.h[5]);
			}
			else
			{new_hapinfo.w[7]= 0;}

			if (new_hapinfo.h[2]*new_hapinfo.h[7]+new_hapinfo.h[3]*new_hapinfo.h[6] >0)
			{
				new_hapinfo.w[8]=
				new_hapinfo.h[2]*new_hapinfo.h[7]/(new_hapinfo.h[2]*new_hapinfo.h[7]+new_hapinfo.h[3]*new_hapinfo.h[6]);
			}
			else
			{new_hapinfo.w[8]= 0;}

			if (new_hapinfo.h[4]*new_hapinfo.h[7]+new_hapinfo.h[5]*new_hapinfo.h[6] >0)
			{
				new_hapinfo.w[9]=
				new_hapinfo.h[4]*new_hapinfo.h[7]/(new_hapinfo.h[4]*new_hapinfo.h[7]+new_hapinfo.h[5]*new_hapinfo.h[6]);
			}
			else
			{new_hapinfo.w[9]= 0;}

			diff=fabs(old_hapinfo.w[0]-new_hapinfo.w[0]);
			for(i=1;i<10;i++)
			{
				if(diff < fabs(old_hapinfo.w[i]-new_hapinfo.w[i]))
				{
					diff=fabs(old_hapinfo.w[i]-new_hapinfo.w[i]);
				}
			}

			//new = old
			for(i=0;i<8;i++)
			{
				old_hapinfo.h[i]=new_hapinfo.h[i];
			}
			for(i=0;i<10;i++)
			{
				old_hapinfo.w[i]=new_hapinfo.w[i];
			}
		}
	} //snp == 3
	new_hapinfo.r2=r2;
	return new_hapinfo;
}


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
					   double **D, double **T, double **U, double **Ut, double **sumPP, double **sumPJ, double **sumPK, double **MMinv, int test,int thread, char *hapfile, int dohapfile, char *hapstring, FILE *file1,int skip,
					   int npplqc, int* PPLMap, uint64_t*** BinSNPs, struct PPLLOCATION* PplLocations,int nMc,struct STATplus result,int caseOnly, int secondSystem, int maxIndexCov, int liabilityCut, int singleMarkerTest, int covNum,int covariancematrix,int df_L1,int df_L2)
{
	int indicator[28+maxIndexCov]; //which parameters (beta0-beta26) and covariates (28-37) are used? is sex a covaitate? (27)
	int numberit[28+maxIndexCov];//numbers the used indicators
	int number2indi[28+maxIndexCov]; //retrieve kth actually used parameter
	//struct STATplus result; //log-Likelihood and all betas
	double logNew;
	double betaNew[28+maxIndexCov];
	double betaOld[28+maxIndexCov];

	int i,k,j,jj,f=0;
	int withcov=0; //are there covariates?

	int n; //number of individuals used
	double AA=0, AB=0, BB=0, CC=0, CD=0, DD=0, EE=0, EF=0, FF=0; //joint genocounts cases and control for SNPs 1,2,3
	double help1=0,help2=0,help3=0;
	int complete=0;
	double currentval[6]; //x-values (main effects) for current individual
	double g1=0,g2=0,g3=0; // indicates if genotypes for SNPs 1,2,3, respectively, are needed
	int snps=0; //1,2 or 3 (highest snp used)
	int x1=0; int x1D=0; //indicates if parameter is used
	int x2=0; int x2D=0; //indicates if parameter is used
	int x1x2=0; int x1x2D=0; int x1Dx2=0; int x1Dx2D=0; //indicates if parameter is used
	int x3=0; int x3D; //indicates if parameter is used
	int x1x3=0; int x1x3D=0; int x1Dx3=0; int x1Dx3D=0; //indicates if parameter is used
	int x2x3=0; int x2x3D=0; int x2Dx3=0; int x2Dx3D=0; //indicates if parameter is used
	int x1x2x3=0;int x1x2x3D=0;int x1x2Dx3=0;int x1x2Dx3D=0;int x1Dx2x3=0;int x1Dx2x3D=0;int x1Dx2Dx3=0;int x1Dx2Dx3D=0; //indicates if parameter is used
	//double exponent=0;
	int dummy;
	double AA_Ca=0,AB_Ca=0,BB_Ca=0,CC_Ca=0,CD_Ca=0,DD_Ca=0,EE_Ca=0,EF_Ca=0,FF_Ca=0;
	double casecounts_ab[3][3]; //2-dim genocounts snps 1,2
	double casecounts_ac[3][3]; //2-dim genocounts snps 1,3
	double casecounts_bc[3][3]; //2-dim genocounts snps 2,3
	int ncase=0;
	//int ncontrol=0;
	//double ncased=0;
	double total=0;
	int dfUnchanged=1;

	//double ssa=0;
	double rss=0;
	double ssy=0;
	double Yavg =0;
	//double fstat=0;
	double n1=0;
	double f1=0;
	//int l=0;
	//double sigma2=0;
	double se_beta[28+maxIndexCov];
	double YtY = 0;
	double YtXAinvXtY = 0;

	struct HAPINFO hapinfo;

	int nwith=0;

	x1=x[1];
	x1D=x[2];
	x2=x[3];
	x2D=x[4];
	x1x2=x[5];
	x1x2D=x[6];
	x1Dx2=x[7];
	x1Dx2D=x[8];
	x3=x[9];
	x3D=x[10];
	x1x3=x[11];
	x1x3D=x[12];
	x1Dx3=x[13];
	x1Dx3D=x[14];
	x2x3=x[15];
	x2x3D=x[16];
	x2Dx3=x[17];
	x2Dx3D=x[18];
	x1x2x3=x[19];
	x1x2x3D=x[20];
	x1x2Dx3=x[21];
	x1x2Dx3D=x[22];
	x1Dx2x3=x[23];
	x1Dx2x3D=x[24];
	x1Dx2Dx3=x[25];
	x1Dx2Dx3D=x[26];

	if(xtype>0 && !female && xtype!=3) //xtype 3 is single marker x chromosome, both sexes
	  {
	    if(x[2] || x[4] || x[6] || x[7] || x[8] || x[10] || x[12] || x[13] || x[14] || x[16] || x[17] || x[18] || x[20] || x[21] || x[22] || x[23] || x[24] || x[25] || x[26])
	      {
		dfUnchanged=0;
	      }

		x1D=0;x2D=0;x3D=0;
		x1Dx2=0;x1x2D=0;x1Dx2D=0;
		x1Dx3=0;x1x3D=0;x1Dx3D=0;
		x2Dx3=0;x2x3D=0;x2Dx3D=0;
		x1x2x3D=0;x1x2Dx3=0;x1x2Dx3D=0;x1Dx2x3=0;x1Dx2x3D=0;x1Dx2Dx3=0;x1Dx2Dx3D=0;
	}

	result.sc=0;

	if(xtype>=3 && !male && !female)
	{
		sexcov=1;
	}
	else if(xtype>=3)
	{
		sexcov=0;
	}
	if(xtype>0 && xtype <3)
	{
		sexcov=0;
	}

	//INITIALIZE INDICATOR
	indicator[0]=1; //always with intercept beta_0
	result.in[0]=1;
	for(i=1;i<28+maxIndexCov;i++){indicator[i]=0;result.in[i]=-1;}
	indicator[27]=sexcov; //sex as covariate yes/no
	if(indicator[27]){result.in[27]=1;}

	// skip=2;


	if(!skip && 0)
	{
		if(N==2) // 2-marker
		{
			for(i=0;i<3;i++)
			{
				for(j=0;j<3;j++)
				{
					casecounts_ab[i][j]= casecounts5[i][j][0][0][0];
				}
			}


			AA_Ca=casecounts_ab[0][0]+casecounts_ab[0][1]+casecounts_ab[0][2];
			AB_Ca=casecounts_ab[1][0]+casecounts_ab[1][1]+casecounts_ab[1][2];
			BB_Ca=casecounts_ab[2][0]+casecounts_ab[2][1]+casecounts_ab[2][2];

			CC_Ca=casecounts_ab[0][0]+casecounts_ab[1][0]+casecounts_ab[2][0];
			CD_Ca=casecounts_ab[0][1]+casecounts_ab[1][1]+casecounts_ab[2][1];
			DD_Ca=casecounts_ab[0][2]+casecounts_ab[1][2]+casecounts_ab[2][2];

			EE_Ca=0;EF_Ca=0;FF_Ca=0;
		}
		else if(N==3)// 3-marker
		{
			for(i=0;i<3;i++)
			{
				for(j=0;j<3;j++)
				{
					casecounts_ab[i][j]=0;
					casecounts_ac[i][j]=0;
					casecounts_bc[i][j]=0;
				}
			}

			for(i=0;i<3;i++)
			{
				for(j=0;j<3;j++)
				{
					for(k=0;k<3;k++)
					{
						casecounts_ab[i][j]+=casecounts5[i][j][k][0][0];
						casecounts_ac[i][k]+=casecounts5[i][j][k][0][0];
						casecounts_bc[j][k]+=casecounts5[i][j][k][0][0];

						if(i==0)
						{
							AA_Ca+=casecounts5[i][j][k][0][0];
						}
						if(i==1)
						{
							AB_Ca+=casecounts5[i][j][k][0][0];
						}
						if(i==2)
						{
							BB_Ca+=casecounts5[i][j][k][0][0];
						}
						if(j==0)
						{
							CC_Ca+=casecounts5[i][j][k][0][0];
						}
						if(j==1)
						{
							CD_Ca+=casecounts5[i][j][k][0][0];
						}
						if(j==2)
						{
							DD_Ca+=casecounts5[i][j][k][0][0];
						}
						if(k==0)
						{
							EE_Ca+=casecounts5[i][j][k][0][0];
						}
						if(k==1)
						{
							EF_Ca+=casecounts5[i][j][k][0][0];
						}
						if(k==2)
						{
							FF_Ca+=casecounts5[i][j][k][0][0];
						}  // if
					} //k
				}//j
			}   // i
		} //if N==3

		if(N>1)
		{
			AA=AA_Ca;AB=AB_Ca;BB=BB_Ca;
			CC=CC_Ca;CD=CD_Ca;DD=DD_Ca;
			EE=EE_Ca;EF=EF_Ca;FF=FF_Ca;
		}
		else
		{
			AA=counts.AA_Ca;
			AB=counts.AB_Ca;
			BB=counts.BB_Ca;
		}

		if (liabilityCut == 0)
		{
			if(!haplo)
			{
				if( (AA==0 && AB==0) || (AA ==0 && BB==0) || (AB ==0 && BB==0) ) // SNP 1 monomorph
				{
					x1=0;x1D=0;
					x1x2=0;x1x2D=0;x1Dx2=0;x1Dx2D=0;
					x1x3=0;x1x3D=0;x1Dx3=0;x1Dx3D=0;
					x1x2x3=0;x1x2x3D=0;x1x2Dx3=0;x1x2Dx3D=0;x1Dx2x3=0;x1Dx2x3D=0;x1Dx2Dx3=0;x1Dx2Dx3D=0;
				}
				else if ( AA==0 || AB==0 || BB==0) // no dominance
				{
					x1D=0;
					x1Dx2=0;x1Dx2D=0;
					x1Dx3=0;x1Dx3D=0;
					x1Dx2x3=0;x1Dx2x3D=0;x1Dx2Dx3=0;x1Dx2Dx3D=0;
				}

				if( (CC==0 && CD==0) || (CC ==0 && DD==0) || (CD ==0 && DD==0) ) // SNP 2 monomorph
				{
					x2=0;x2D=0;
					x1x2=0;x1x2D=0;x1Dx2=0;x1Dx2D=0;
					x2x3=0;x2x3D=0;x2Dx3=0;x2Dx3D=0;
					x1x2x3=0;x1x2x3D=0;x1x2Dx3=0;x1x2Dx3D=0;x1Dx2x3=0;x1Dx2x3D=0;x1Dx2Dx3=0;x1Dx2Dx3D=0;
				}
				else if (CC==0 || CD==0 || DD==0 ) // only 2 genotypes
				{
					x2D=0;//printf("here x2D\n");
					x1x2D=0;x1Dx2D=0;
					x2Dx3=0;x2Dx3D=0;
					x1x2Dx3=0;x1x2Dx3D=0;x1Dx2Dx3=0;x1Dx2Dx3D=0;
				}

				if( (EE==0 && EF==0) || (EE == 0 && FF==0) || (EF ==0 && FF==0) ) // SNP 3 monomorph
				{
					x3=0;x3D=0;
					x1x3=0;x1x3D=0;x1Dx3=0;x1Dx3D=0;
					x2x3=0;x2x3D=0;x2Dx3=0;x2Dx3D=0;
					x1x2x3=0;x1x2x3D=0;x1x2Dx3=0;x1x2Dx3D=0;x1Dx2x3=0;x1Dx2x3D=0;x1Dx2Dx3=0;x1Dx2Dx3D=0;
				}
				else if (EE==0 || EF==0 || FF==0) // only 2 genotypes
				{
					x3D=0;
					x2x3D=0;x2Dx3D=0;
					x1x3D=0;x1Dx3D=0;
					x1x2x3D=0;x1x2Dx3D=0;x1Dx2x3D=0;x1Dx2Dx3D=0;
				}
			}
		}

		//now check further interaction terms:

	if(N>1 && !haplo)
	{
		total=0;
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				total+=casecounts_ab[i][j];
				if(casecounts_ab[i][j]>0){nwith++;}
			}
		}

		if(nwith<=2 || (nwith==3 && casecounts_ab[1][1]>0))
		  {
		   x1x2=0;x1x2D=0;x1Dx2=0;x1Dx2D=0;x1x3=0;
		   x1x3D=0;x1Dx3=0;x1Dx3D=0;x2x3=0;x2x3D=0;x2Dx3=0;x2Dx3D=0;
		   x1x2x3=0;x1x2x3D=0;x1x2Dx3=0;x1x2Dx3D=0;x1Dx2x3=0;x1Dx2x3D=0;x1Dx2Dx3=0;
		   x1Dx2Dx3D=0;
		  }


		help1=casecounts_ab[0][0]+casecounts_ab[2][2];
		help2=casecounts_ab[2][0]+casecounts_ab[0][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total)
		{
			x1x2=0;x1x2x3=0;x1x2x3D=0;
		}
		help1=casecounts_ab[0][0]+
		casecounts_ab[2][1]+casecounts_ab[0][2];
		help2=casecounts_ab[2][0]+
		casecounts_ab[0][1]+casecounts_ab[2][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total)
		{

			x1x2D=0;x1x2Dx3=0;x1x2Dx3D=0;

		}


		help1=casecounts_ab[0][0]+
		casecounts_ab[1][2]+casecounts_ab[2][0];
		help2=casecounts_ab[0][2]+
		casecounts_ab[1][0]+casecounts_ab[2][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total)
		{
			x1Dx2=0;x1Dx2x3=0;x1Dx2x3D=0;
		}


		help1=casecounts_ab[0][0]+
		casecounts_ab[0][2]+casecounts_ab[1][1]+
		casecounts_ab[2][0]+casecounts_ab[2][2];

		help2=total-help1;

		if(help1==total || help2==total)
		{
			x1Dx2D=0;x1Dx2Dx3=0;x1Dx2Dx3D=0;
		}
	} //end if N>1


	if(N>2 && !haplo)
	{
		//snps 1 and 3
		total=0;
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				total+=casecounts_ac[i][j];
			}
		}

		help1=casecounts_ac[0][0]+casecounts_ac[2][2];
		help2=casecounts_ac[2][0]+casecounts_ac[0][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total){x1x3=0;x1x2x3=0;x1x2Dx3=0;}


		help1=casecounts_ac[0][0]+ casecounts_ac[2][1]+casecounts_ac[0][2];
		help2=casecounts_ac[2][0]+ casecounts_ac[0][1]+casecounts_ac[2][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total){x1x3D=0;x1x2x3D=0;x1x2Dx3D=0;}


		help1=casecounts_ac[0][0]+ casecounts_ac[1][2]+casecounts_ac[2][0];
		help2=casecounts_ac[0][2]+ casecounts_ac[1][0]+casecounts_ac[2][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total){x1Dx3=0;x1Dx2x3=0;x1Dx2Dx3=0;}


		help1=casecounts_ac[0][0]+ casecounts_ac[0][2]+casecounts_ac[1][1]+
		casecounts_ac[2][0]+casecounts_ac[2][2];

		help2=total-help1;

		if(help1==total || help2==total){x1Dx3D=0;x1Dx2x3D=0;x1Dx2Dx3D=0;}


		//snps 2 and 3
		total=0;
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				total+=casecounts_bc[i][j];
			}
		}

		help1=casecounts_bc[0][0]+casecounts_bc[2][2];
		help2=casecounts_bc[2][0]+casecounts_bc[0][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total){x2x3=0;x1x2x3=0;x1Dx2x3=0;}


		help1=casecounts_bc[0][0]+ casecounts_bc[2][1] +casecounts_bc[0][2];
		help2=casecounts_bc[2][0]+ casecounts_bc[0][1] +casecounts_bc[2][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total){x2x3D=0;x1x2x3D=0;x1Dx2x3D=0;}


		help1=casecounts_bc[0][0]+ casecounts_bc[1][2]+casecounts_bc[2][0];
		help2=casecounts_bc[0][2]+ casecounts_bc[1][0]+casecounts_bc[2][2];
		help3=total-help2-help1;

		if(help1==total || help2==total || help3==total){x2Dx3=0;x1x2Dx3=0;x1Dx2Dx3=0;}


		help1=casecounts_bc[0][0]+ casecounts_bc[0][2]+casecounts_bc[1][1]+
		casecounts_bc[2][0]+casecounts_bc[2][2];

		help2=total-help1;

		if(help1==total || help2==total){x2Dx3D=0;x1x2Dx3D=0;x1Dx2Dx3D=0;}

		//now 3-fold terms
		total=0;
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				for(jj=0;jj<3;jj++)
				{
					total+=casecounts5[i][j][jj][0][0];
				}
			}
		}

		if(x1x2x3) //update if always equal
		{
			help1= casecounts5[0][2][0][0][0]+casecounts5[2][0][0][0][0]+casecounts5[0][0][2][0][0]+casecounts5[2][2][2][0][0];

			help2= casecounts5[0][0][0][0][0]+casecounts5[2][2][0][0][0]+casecounts5[0][2][2][0][0]+casecounts5[2][0][2][0][0];

			help3=total-help1-help2;

			if(help1==total || help2==total || help3==total){x1x2x3=0;}
		}

		if(x1x2x3D) //update if always equal
		{
			help1= casecounts5[0][2][0][0][0]+casecounts5[2][0][0][0][0]+casecounts5[0][0][1][0][0]+casecounts5[2][2][1][0][0]+casecounts5[0][2][2][0][0]+casecounts5[2][0][2][0][0];

			help2= casecounts5[0][0][0][0][0]+casecounts5[2][2][0][0][0]+casecounts5[0][2][1][0][0]+casecounts5[2][0][1][0][0]+
				casecounts5[0][0][2][0][0]+casecounts5[2][2][2][0][0];

			help3=total-help1-help2;

			if(help1==total || help2==total || help3==total){x1x2x3D=0;}
		}

		if(x1x2Dx3) //update if always equal
		{
			help1= casecounts5[0][1][0][0][0]+casecounts5[2][0][0][0][0]+casecounts5[2][2][0][0][0]+casecounts5[0][0][2][0][0]+
				casecounts5[0][2][2][0][0]+casecounts5[2][1][2][0][0];

			help2= casecounts5[0][0][0][0][0]+casecounts5[0][2][0][0][0]+casecounts5[2][1][0][0][0]+casecounts5[0][1][2][0][0]+
				casecounts5[2][0][2][0][0]+casecounts5[2][2][2][0][0];

			help3=total-help1-help2;


			if(help1==total || help2==total || help3==total){x1x2Dx3=0;}
		}

		if(x1x2Dx3D) //update if always equal
		{
			help1=casecounts5[0][1][0][0][0]+casecounts5[2][0][0][0][0]+casecounts5[2][2][0][0][0]+
			casecounts5[0][0][1][0][0]+casecounts5[0][2][1][0][0]+casecounts5[2][1][1][0][0]+
			casecounts5[0][1][2][0][0]+casecounts5[2][0][2][0][0]+casecounts5[2][2][2][0][0];

			help2=casecounts5[0][0][0][0][0]+casecounts5[0][2][0][0][0]+casecounts5[2][1][0][0][0]+
			casecounts5[0][1][1][0][0]+casecounts5[2][0][1][0][0]+casecounts5[2][2][1][0][0]+
			casecounts5[0][0][2][0][0]+casecounts5[0][2][2][0][0]+casecounts5[2][1][2][0][0];

			help3=total-help1-help2;

			if(help1==total || help2==total || help3==total){x1x2Dx3D=0;}
		}

		if(x1Dx2x3) //update if always equal
		{
			help1= casecounts5[0][2][0][0][0]+casecounts5[1][0][0][0][0]+casecounts5[2][2][0][0][0]+casecounts5[0][0][2][0][0]+
				casecounts5[1][2][2][0][0]+casecounts5[2][0][2][0][0];

			help2= casecounts5[0][0][0][0][0]+casecounts5[1][2][0][0][0]+casecounts5[2][0][0][0][0]+casecounts5[0][2][2][0][0]+
				casecounts5[1][0][2][0][0]+casecounts5[2][2][2][0][0];

			help3=total-help1-help2;

			if(help1==total || help2==total || help3==total){x1Dx2x3=0;}
		}

		if(x1Dx2x3D) //update if always equal
		{
			help1=casecounts5[0][2][0][0][0]+casecounts5[1][0][0][0][0]+casecounts5[2][2][0][0][0]+
				casecounts5[0][0][1][0][0]+casecounts5[1][2][1][0][0]+casecounts5[2][0][1][0][0]+
				casecounts5[0][2][2][0][0]+casecounts5[1][0][2][0][0]+casecounts5[2][2][2][0][0];

			help2=casecounts5[0][0][0][0][0]+casecounts5[1][2][0][0][0]+casecounts5[2][0][0][0][0]+
				casecounts5[0][2][1][0][0]+casecounts5[1][0][1][0][0]+casecounts5[2][2][1][0][0]+
				casecounts5[0][0][2][0][0]+casecounts5[1][2][2][0][0]+casecounts5[2][0][2][0][0];

			help3=total-help1-help2;

			if(help1==total || help2==total || help3==total){x1Dx2x3D=0;}
		}

		if(x1Dx2Dx3) //update if always equal
		{
			help1=casecounts5[0][1][0][0][0]+casecounts5[1][0][0][0][0]+casecounts5[1][2][0][0][0]+
				casecounts5[2][1][0][0][0]+casecounts5[0][0][2][0][0]+casecounts5[0][2][2][0][0]+
				casecounts5[1][1][2][0][0]+casecounts5[2][0][2][0][0]+casecounts5[2][2][2][0][0];

			help2=casecounts5[0][0][0][0][0]+casecounts5[0][2][0][0][0]+casecounts5[1][1][0][0][0]+
				casecounts5[2][0][0][0][0]+casecounts5[2][2][0][0][0]+casecounts5[0][1][2][0][0]+
				casecounts5[1][0][2][0][0]+casecounts5[1][2][2][0][0]+casecounts5[2][1][2][0][0];

			help3=total-help1-help2;

			if(help1==total || help2==total || help3==total){x1Dx2Dx3=0;}
		}

		if(x1Dx2Dx3D) //update if always equal
		{
			help1= casecounts5[0][1][0][0][0]+casecounts5[1][0][0][0][0]+
				casecounts5[1][2][0][0][0]+
				casecounts5[2][1][0][0][0]+casecounts5[0][0][1][0][0]+
				casecounts5[0][2][1][0][0]+
				casecounts5[1][1][1][0][0]+casecounts5[2][0][1][0][0]+
				casecounts5[2][2][1][0][0]+
				casecounts5[0][1][2][0][0]+casecounts5[1][0][2][0][0]+
				casecounts5[1][2][2][0][0]+
				casecounts5[2][1][2][0][0];

			help2=total-help1;

			if(help1==total || help2==total){x1Dx2Dx3D=0;}
		}
	} //end if snps>2


		//SET INDICATOR
		if(!haplo)
		{
			g1=x1+x1D+x1x2+x1x2D+x1Dx2+x1Dx2D+
			x1x3+x1x3D+x1Dx3+x1Dx3D+
			x1x2x3+x1x2x3D+x1x2Dx3+x1x2Dx3D+x1Dx2x3+x1Dx2x3D+x1Dx2Dx3+x1Dx2Dx3D;

			g2=x2+x2D+x1x2+x1x2D+x1Dx2+x1Dx2D+
			x2x3+x2x3D+x2Dx3+x2Dx3D+
			x1x2x3+x1x2x3D+x1x2Dx3+x1x2Dx3D+x1Dx2x3+x1Dx2x3D+x1Dx2Dx3+x1Dx2Dx3D;

			g3=x3+x3D+x1x3+x1x3D+x1Dx3+x1Dx3D+x2x3+x2x3D+x2Dx3+x2Dx3D+
			x1x2x3+x1x2x3D+x1x2Dx3+x1x2Dx3D+x1Dx2x3+x1Dx2x3D+x1Dx2Dx3+x1Dx2Dx3D;

			if(g3){snps=3;}
			else if (g2){snps=2;}
			else if (g1){snps=1;}
		}
		else
		{
			g1=1;g2=1;snps=2;
			if(N>2){g3=1;snps=3;}
		}
	}
	else //skip, used with pathwayTest 5

	{

		g1=x1+x1D+x1x2+x1x2D+x1Dx2+x1Dx2D+

		x1x3+x1x3D+x1Dx3+x1Dx3D+
		x1x2x3+x1x2x3D+x1x2Dx3+x1x2Dx3D+x1Dx2x3+x1Dx2x3D+x1Dx2Dx3+x1Dx2Dx3D;

		g2=x2+x2D+x1x2+x1x2D+x1Dx2+x1Dx2D+
		x2x3+x2x3D+x2Dx3+x2Dx3D+
		x1x2x3+x1x2x3D+x1x2Dx3+x1x2Dx3D+x1Dx2x3+x1Dx2x3D+x1Dx2Dx3+x1Dx2Dx3D;

		g3=x3+x3D+x1x3+x1x3D+x1Dx3+x1Dx3D+x2x3+x2x3D+x2Dx3+x2Dx3D+
		x1x2x3+x1x2x3D+x1x2Dx3+x1x2Dx3D+x1Dx2x3+x1Dx2x3D+x1Dx2Dx3+x1Dx2Dx3D;

		if(g3){snps=3;}
		else if (g2){snps=2;}
		else if (g1){snps=1;}
	}



	help1=0;help2=0;help3=0;
	if(x1==1){indicator[1]=1;result.in[1]=1;}
	if(x1D==1){indicator[2]=1;result.in[2]=1;}
	if(snps>1)
	{
		if(x2==1){indicator[3]=1;result.in[3]=1;}
		if(x2D==1){indicator[4]=1;result.in[4]=1;}
		if(x1x2==1){indicator[5]=1;result.in[5]=1;}
		if(x1x2D==1){indicator[6]=1;result.in[6]=1;}
		if(x1Dx2==1){indicator[7]=1;result.in[7]=1;}
		if(x1Dx2D==1){indicator[8]=1;result.in[8]=1;}
	}
	if(snps>2)
	{
		if(x3==1){indicator[9]=1;result.in[9]=1;}
		if(x3D==1){indicator[10]=1;result.in[10]=1;}

		if(x1x3==1){indicator[11]=1;result.in[11]=1;}
		if(x1x3D==1){indicator[12]=1;result.in[12]=1;}
		if(x1Dx3==1){indicator[13]=1;result.in[13]=1;}
		if(x1Dx3D==1){indicator[14]=1;result.in[14]=1;}
		if(x2x3==1){indicator[15]=1;result.in[15]=1;}

		if(x2x3D==1){indicator[16]=1;result.in[16]=1;}
		if(x2Dx3==1){indicator[17]=1;result.in[17]=1;}
		if(x2Dx3D==1){indicator[18]=1;result.in[18]=1;}
		if(x1x2x3==1){indicator[19]=1;result.in[19]=1;}
		if(x1x2x3D==1){indicator[20]=1;result.in[20]=1;}
		if(x1x2Dx3==1){indicator[21]=1;result.in[21]=1;}
		if(x1x2Dx3D==1){indicator[22]=1;result.in[22]=1;}
		if(x1Dx2x3==1){indicator[23]=1;result.in[23]=1;}
		if(x1Dx2x3D==1){indicator[24]=1;result.in[24]=1;}
		if(x1Dx2Dx3==1){indicator[25]=1;result.in[25]=1;}
		if(x1Dx2Dx3D==1){indicator[26]=1;result.in[26]=1;}
	}

	for(i=28;i<28+maxIndexCov;i++)
	{
		if(cov[i-28]==1){indicator[i]=1;withcov=1;result.in[i]=1;}
	}


	if(haplo && N>1) //compute haplotype frequencies
	{
		hapinfo=haplo_em(casecounts5,N,alt,hapfile,dohapfile,hapstring, file1);
	}


	//determine df
	numberit[0]=0;number2indi[0]=0;
	k=0;

	for(i=1;i<28+maxIndexCov;i++)
	{
		f+=indicator[i];//df
		if(indicator[i]==1)
		{
			k++;
			numberit[i]=k;
			number2indi[k]=i;
		}
		else
		{
			numberit[i]=-1;
		}
	}
	double newf=f;

	//moved also here from below in version 682
	for(i=0;i<f+1;i++)
	{
	    //NEW
		betaNew[number2indi[i]]=-1;
		//if(newf>0){betaNew[number2indi[i]]=YtXAinv[0][i];}
		if(nMc==0)
		{
			result.betaNew_se[number2indi[i]] = -1;
			result.betaNew_lcl[number2indi[i]] = -1;
			result.betaNew_rcl[number2indi[i]] = -1;

		}
		result.b[number2indi[i]] = -1;
		//printf("here2\n");
	}

	//get n, Y and X
	n=0;ncase=0;



	struct PPLLOCATION guy;
	for(int kMod=0;kMod<npplqc;kMod++)
	{
		guy = PplLocations[kMod];
		k = PPLMap[kMod];
		complete=1;

		if(alt==0 && Yhelp[k]< -0.1){continue;}

		//person in?
		if(alt==1)
		{

		    if(caseOnly && person[k].qtaff[thread]!=2){Yhelp[k]=-1;continue;}
			if(xtype==0)
			{
				if(person[k].qcin==1)
				{

					Y[n]=person[k].qtaff[thread];
					Yhelp[k]=1;
				}
				else {Yhelp[k]=-1;continue;}

			}
			else if(xtype==1)
			{
				if(person[k].sex==1 && person[k].qcin==1)
				{
					Y[n]=person[k].qtaff[thread];Yhelp[k]=1;
				}
				else {Yhelp[k]=-1;continue;}
			}
			else if(xtype==2)
			{
				if(person[k].sex==2 && person[k].qcin==1)
				{
					Y[n]=person[k].qtaff[thread];Yhelp[k]=1;
				}
				else {Yhelp[k]=-1;continue;}
			}
			else if(xtype>=3)
			{
				//if((person[k].sex==0) || (person[k].sex==1 && M[a].d2[k].d3==2)){Yhelp[k]=-1;continue;}
				if((person[k].sex==0) || (person[k].sex==1 && getbit64(BinSNPs[a][guy.nr][2],guy.pos))){Yhelp[k]=-1;continue;}
				if(person[k].qcin==1)
				{
					Y[n]=person[k].qtaff[thread];Yhelp[k]=1;
				}
				else {Yhelp[k]=-1;continue;}
			}

		}
		else //! alt=1
		{
			Y[n]=person[k].qtaff[thread];
		}

		//get main x values
		if(g1)
		{
			if(getbit64(BinSNPs[a][guy.nr][1],guy.pos))
			{
				currentval[0]=1;currentval[1]=-0.5;
				if(singleMarkerTest == 7){currentval[0]=0;}
			}
			else if(getbit64(BinSNPs[a][guy.nr][2],guy.pos))
			{
				currentval[0]=0;currentval[1]=0.5;
				if(singleMarkerTest == 6 || singleMarkerTest == 7){currentval[0]=1;}
			}
			else if(getbit64(BinSNPs[a][guy.nr][3],guy.pos))
			{
				currentval[0]=-1;currentval[1]=-0.5;
				if(singleMarkerTest == 5 || singleMarkerTest == 6 || singleMarkerTest == 7){currentval[0]=0;}
			}
			else
			{
				Yhelp[k]=-1;
				continue;
			}
			if(covNum>=0)
			  {
			   if(covNum>0){currentval[0]=person[k].cov[covNum-1];}
			   else {currentval[0]=person[k].sex;}
			  }
		}

		if(g2)
		{
			if(getbit64(BinSNPs[b][guy.nr][1],guy.pos))
			{
				currentval[2]=1;currentval[3]=-0.5;
			}
			else if(getbit64(BinSNPs[b][guy.nr][2],guy.pos))
			{
				currentval[2]=0;currentval[3]=0.5;
			}
			else if(getbit64(BinSNPs[b][guy.nr][3],guy.pos))
			{
				currentval[2]=-1;currentval[3]=-0.5;
			}
			else
			{
				Yhelp[k]=-1;
				continue;
			}
		}

		if(g3)
		{
			if(getbit64(BinSNPs[c][guy.nr][1],guy.pos))
			{
				currentval[4]=1;currentval[5]=-0.5;
			}
			else if(getbit64(BinSNPs[c][guy.nr][2],guy.pos))
			{
				currentval[4]=0;currentval[5]=0.5;
			}
			else if(getbit64(BinSNPs[c][guy.nr][3],guy.pos))
			{
				currentval[4]=-1;currentval[5]=-0.5;
			}
			else
			{
				Yhelp[k]=-1;
				continue;
			}
		}

		//update X
		if(haplo && N==2)
		{
			//newhaplo
			X[n][0]=1;
			if(currentval[0]==-1 && currentval[2]==-1)
			{
				X[n][1]=2;X[n][2]=0;X[n][3]=0;X[n][4]=0;
			}
			else if(currentval[0]==-1 && currentval[2]==0)
			{
				X[n][1]=1;X[n][2]=1;X[n][3]=0;X[n][4]=0;
			}
			else if(currentval[0]==-1 && currentval[2]==1)
			{
				X[n][1]=0;X[n][2]=2;X[n][3]=0;X[n][4]=0;
			}
			else if(currentval[0]==0 && currentval[2]==-1)
			{
				X[n][1]=1;X[n][2]=0;X[n][3]=1;X[n][4]=0;
			}
			else if(currentval[0]==0 && currentval[2]==0)
			{
				X[n][1]=hapinfo.w[0];X[n][2]=1-hapinfo.w[0];
				X[n][3]=1-hapinfo.w[0];X[n][4]=hapinfo.w[0];
			}
			else if(currentval[0]==0 && currentval[2]==1)
			{
				X[n][1]=0;X[n][2]=1;X[n][3]=0;X[n][4]=1;
			}
			else if(currentval[0]==1 && currentval[2]==-1)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=2;X[n][4]=0;
			}
			else if(currentval[0]==1 && currentval[2]==0)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=1;X[n][4]=1;
			}
			else if(currentval[0]==1 && currentval[2]==1)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=0;X[n][4]=2;
			}
			for(j=0;j<=f;j++)
			{
			    //newhaplo

				jj=number2indi[j];
				if (jj==27)
				{
					if(person[k].sex==1){X[n][j]=1;}
					else if(person[k].sex==2){X[n][j]=0;}
					else {Yhelp[k]=-1;complete=0;break;}
				}
				else if(jj>27)
				{
					if(person[k].covin[jj-28]==0) //NEU
					{
						Yhelp[k]=-1; complete=0;break;
					}
					else{X[n][j]=person[k].cov[jj-28];}
				}
			}
		} //haplo2
		else if(haplo && N==3)
		{
		    //double check=0;
		    //newhaplo
		    X[n][0]=1;
			if(currentval[0]==-1 && currentval[2]==-1 && currentval[4]==-1)
			{
				X[n][1]=2;X[n][2]=0;X[n][3]=0;X[n][4]=0;X[n][5]=0;X[n][6]=0;X[n][7]=0;X[n][8]=0;
			}
			else if(currentval[0]==-1 && currentval[2]==-1 && currentval[4]==0)
			{
				X[n][1]=1;X[n][2]=1;X[n][3]=0;X[n][4]=0;X[n][5]=0;X[n][6]=0;X[n][7]=0;X[n][8]=0;
			}
			else if(currentval[0]==-1 && currentval[2]==-1 && currentval[4]==1)
			{
				X[n][1]=0;X[n][2]=2;X[n][3]=0;X[n][4]=0;X[n][5]=0;X[n][6]=0;X[n][7]=0;X[n][8]=0;
			}
			else if(currentval[0]==-1 && currentval[2]==0 && currentval[4]==-1)
			{
				X[n][1]=1;X[n][2]=0;X[n][3]=1;X[n][4]=0;X[n][5]=0;X[n][6]=0;X[n][7]=0;X[n][8]=0;
			}
			else if(currentval[0]==-1 && currentval[2]==0 && currentval[4]==0)
			{
				X[n][1]=hapinfo.w[2];X[n][2]=(1-hapinfo.w[2]);X[n][3]=(1-hapinfo.w[2]);X[n][4]=hapinfo.w[2];X[n][5]=0;X[n][6]=0;X[n][7]=0;X[n][8]=0;
			}
			else if(currentval[0]==-1 && currentval[2]==0 && currentval[4]==1)
			{
				X[n][1]=0;X[n][2]=1;X[n][3]=0;X[n][4]=1;X[n][5]=0;X[n][6]=0;X[n][7]=0;X[n][8]=0;
			}
			else if(currentval[0]==-1 && currentval[2]==1 && currentval[4]==-1)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=2;X[n][4]=0;X[n][5]=0;X[n][6]=0;X[n][7]=0;X[n][8]=0;
			}
			else if(currentval[0]==-1 && currentval[2]==1 && currentval[4]==0)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=1;X[n][4]=1;X[n][5]=0;X[n][6]=0;X[n][7]=0;X[n][8]=0;
			}
			else if(currentval[0]==-1 && currentval[2]==1 && currentval[4]==1)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=0;X[n][4]=2;X[n][5]=0;X[n][6]=0;X[n][7]=0;X[n][8]=0;
			}
			else if(currentval[0]==0 && currentval[2]==-1 && currentval[4]==-1)
			{
				X[n][1]=1;X[n][2]=0;X[n][3]=0;X[n][4]=0;X[n][5]=1;X[n][6]=0;X[n][7]=0;X[n][8]=0;
			}
			else if(currentval[0]==0 && currentval[2]==-1 && currentval[4]==0)
			{
				X[n][1]=hapinfo.w[1];X[n][2]=1-hapinfo.w[1];X[n][3]=0;X[n][4]=0;X[n][5]=1-hapinfo.w[1];X[n][6]=hapinfo.w[1];X[n][7]=0;X[n][8]=0;
			}
			else if(currentval[0]==0 && currentval[2]==-1 && currentval[4]==1)
			{
				X[n][1]=0;X[n][2]=1;X[n][3]=0;X[n][4]=0;X[n][5]=0;X[n][6]=1;X[n][7]=0;X[n][8]=0;
			}
			else if(currentval[0]==0 && currentval[2]==0 && currentval[4]==-1)
			{
				X[n][1]=hapinfo.w[0];X[n][2]=0;X[n][3]=1-hapinfo.w[0];X[n][4]=0;X[n][5]=1-hapinfo.w[0];X[n][6]=0;X[n][7]=hapinfo.w[0];X[n][8]=0;
			}
			else if(currentval[0]==0 && currentval[2]==0 && currentval[4]==0)
			{
				X[n][1]=hapinfo.w[3];X[n][2]=hapinfo.w[4];X[n][3]=hapinfo.w[5];X[n][4]=hapinfo.w[6];X[n][5]=hapinfo.w[6];X[n][6]=hapinfo.w[5];X[n][7]=hapinfo.w[4];X[n][8]=hapinfo.w[3];
			}
			else if(currentval[0]==0 && currentval[2]==0 && currentval[4]==1)
			{
				X[n][1]=0;X[n][2]=hapinfo.w[7];X[n][3]=0;X[n][4]=1-hapinfo.w[7];X[n][5]=0;X[n][6]=1-hapinfo.w[7];X[n][7]=0;X[n][8]=hapinfo.w[7];
			}
			else if(currentval[0]==0 && currentval[2]==1 && currentval[4]==-1)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=1;X[n][4]=0;X[n][5]=0;X[n][6]=0;X[n][7]=1;X[n][8]=0;
			}
			else if(currentval[0]==0 && currentval[2]==1 && currentval[4]==0)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=hapinfo.w[8];X[n][4]=1-hapinfo.w[8];X[n][5]=0;X[n][6]=0;X[n][7]=1-hapinfo.w[8];X[n][8]=hapinfo.w[8];
			}
			else if(currentval[0]==0 && currentval[2]==1 && currentval[4]==1)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=0;X[n][4]=1;X[n][5]=0;X[n][6]=0;X[n][7]=0;X[n][8]=1;
			}
			else if(currentval[0]==1 && currentval[2]==-1 && currentval[4]==-1)

			{
				X[n][1]=0;X[n][2]=0;X[n][3]=0;X[n][4]=0;X[n][5]=2;X[n][6]=0;X[n][7]=0;X[n][8]=0;
			}
			else if(currentval[0]==1 && currentval[2]==-1 && currentval[4]==0)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=0;X[n][4]=0;X[n][5]=1;X[n][6]=1;X[n][7]=0;X[n][8]=0;
			}
			else if(currentval[0]==1 && currentval[2]==-1 && currentval[4]==1)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=0;X[n][4]=0;X[n][5]=0;X[n][6]=2;X[n][7]=0;X[n][8]=0;
			}
			else if(currentval[0]==1 && currentval[2]==0 && currentval[4]==-1)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=0;X[n][4]=0;X[n][5]=1;X[n][6]=0;X[n][7]=1;X[n][8]=0;
			}
			else if(currentval[0]==1 && currentval[2]==0 && currentval[4]==0)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=0;X[n][4]=0;X[n][5]=hapinfo.w[9];X[n][6]=1-hapinfo.w[9];X[n][7]=1-hapinfo.w[9];X[n][8]=hapinfo.w[9];
			}
			else if(currentval[0]==1 && currentval[2]==0 && currentval[4]==1)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=0;X[n][4]=0;X[n][5]=0;X[n][6]=1;X[n][7]=0;X[n][8]=1;
			}
			else if(currentval[0]==1 && currentval[2]==1 && currentval[4]==-1)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=0;X[n][4]=0;X[n][5]=0;X[n][6]=0;X[n][7]=2;X[n][8]=0;
			}
			else if(currentval[0]==1 && currentval[2]==1 && currentval[4]==0)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=0;X[n][4]=0;X[n][5]=0;X[n][6]=0;X[n][7]=1;X[n][8]=1;
			}
			else if(currentval[0]==1 && currentval[2]==1 && currentval[4]==1)
			{
				X[n][1]=0;X[n][2]=0;X[n][3]=0;X[n][4]=0;X[n][5]=0;X[n][6]=0;X[n][7]=0;X[n][8]=2;
			}


			for(j=0;j<=f;j++)
			{
				jj=number2indi[j];
				if (jj==27)
				{
					if(person[k].sex==1){X[n][j]=1;}
					else if(person[k].sex==2){X[n][j]=0;}
					else {Yhelp[k]=-1;complete=0;break;}
				}
				else if(jj>27)
				{
					if(person[k].covin[jj-28]==0) //NEU
					{
						Yhelp[k]=-1; complete=0;break;
					}
					else{X[n][j]=person[k].cov[jj-28];}
				}
			}
		}  //haplo3
		else //!haplo:
		{
			for(j=0;j<=f;j++)
			{
				jj=number2indi[j];
				if(jj==0) //x0
				{
					X[n][j]=1;
				}
				else if(jj==1 && liabilityCut == 0) //x1, no liability model
				{
					X[n][j]=currentval[0];
				}
				else if(jj==1 && liabilityCut != 0) //x1, liability model
				{
					if(person[k].load>=liabilityCut && liabilityCut != -1)
					{
						X[n][j]=1;
					}
					else
					{
						X[n][j]=0;
					}
					if(liabilityCut == -1)
					{
						X[n][j]=pow(person[k].load/(((double)maxIndexCov)*2),3);
					}
				}
				else if(jj==2) //x1D
				{
					X[n][j]=currentval[1];
				}
				else if(jj==3) //x2
				{
					X[n][j]=currentval[2];
				}
				else if(jj==4) //x2D
				{
					X[n][j]=currentval[3];
				}
				else if(jj==5) //x1x2
				{
					X[n][j]=currentval[0]*currentval[2];
				}
				else if(jj==6) //x1x2D
				{
					X[n][j]=currentval[0]*currentval[3];
				}
				else if(jj==7) //x1Dx2
				{
					X[n][j]=currentval[1]*currentval[2];
				}
				else if(jj==8) //x1Dx2D
				{
					X[n][j]=currentval[1]*currentval[3];
				}
				else if(jj==9) //x3
				{
					X[n][j]=currentval[4];
				}
				else if(jj==10) //x3D
				{
					X[n][j]=currentval[5];
				}
				else if(jj==11) //x1x3
				{
					X[n][j]=currentval[0]*currentval[4];
				}
				else if(jj==12) //x1x3D
				{
					X[n][j]=currentval[0]*currentval[5];
				}
				else if(jj==13) //x1Dx3
				{
					X[n][j]=currentval[1]*currentval[4];
				}
				else if(jj==14) //x1Dx3D
				{
					X[n][j]=currentval[1]*currentval[5];
				}
				else if(jj==15) //x2x3
				{
					X[n][j]=currentval[2]*currentval[4];
				}
				else if(jj==16) //x2x3D
				{
					X[n][j]=currentval[2]*currentval[5];
				}
				else if(jj==17) //x2Dx3
				{
					X[n][j]=currentval[3]*currentval[4];
				}
				else if(jj==18) //x2Dx3D
				{
					X[n][j]=currentval[3]*currentval[5];
				}
				else if(jj==19) //x1x2x3
				{
					X[n][j]=currentval[0]*currentval[2]*currentval[4];
				}
				else if(jj==20) //x1x2x3D
				{
					X[n][j]=currentval[0]*currentval[2]*currentval[5];
				}
				else if(jj==21) //x1x2Dx3
				{
					X[n][j]=currentval[0]*currentval[3]*currentval[4];
				}
				else if(jj==22) //x1x2Dx3D
				{
					X[n][j]=currentval[0]*currentval[3]*currentval[5];
				}
				else if(jj==23) //x1Dx2x3
				{
					X[n][j]=currentval[1]*currentval[2]*currentval[4];
				}
				else if(jj==24) //x1Dx2x3D
				{
					X[n][j]=currentval[1]*currentval[2]*currentval[5];
				}
				else if(jj==25) //x1Dx2Dx3
				{
					X[n][j]=currentval[1]*currentval[3]*currentval[4];
				}
				else if(jj==26) //x1Dx2Dx3D
				{
					X[n][j]=currentval[1]*currentval[3]*currentval[5];
				}
				else if (jj==27)
				{
					if(person[k].sex==1){X[n][j]=1;}
					else if(person[k].sex==2){X[n][j]=0;}
					else /* if(alt==1) */ {Yhelp[k]=-1;complete=0;break;}
					//else{complete=0;break;}  never gets there
				}
				else if(jj>27)
			    {
					if(person[k].covin[jj-28]==0) //NEU
					{
						Yhelp[k]=-1; complete=0;break;
					}
					else
					{
						X[n][j]=person[k].cov[jj-28];
					}

					if (jj < maxIndexCov)
					{
						X[n][j]=person[k].cov[jj-28];
					}
			    }
			}  //end j-loop, nothaplo (paramters)
		} //else, not haplo

		if(complete)
		{
			n++;
			ncase++;
		}

	} //end k-loop (individuals)


	for(i=0;i<28+maxIndexCov;i++){betaNew[i]=-1;se_beta[i]=-1;betaOld[i]=-1;}
	for(j=0;j<=f;j++)
	{
		betaNew[number2indi[j]]=-1;se_beta[number2indi[j]]=-1;
	}
	betaNew[0]=0;se_beta[0]=0;

	n1=n;
	f1=f;


	//CASE_ONLY: change matrix
	  if(caseOnly)
	    {
		 if(test==17)  //caseOnly1DF
		   {
		    for (i=0;i<n;i++)
    	        {
    	         Y[i]=X[i][1];
				 X[i][1]=1;
		        }
		   }
		else if(test==18)  //caseOnly1DF
		   {
		   if(secondSystem==0)
		     {
				for (i=0;i<n;i++)
					{
					 Y[i]=X[i][1];
					 X[i][1]=1;
					}
			 }
            else
		     {
				for (i=0;i<n;i++)
					{
					 Y[i]=X[i][2];
					 X[i][2]=1;
					}
			 }
		   } //18
		}
	// END CASE_ONLY

	for(i=0;i<n;i++)
	{
		Yavg+=Y[i];
	}

	if(n>0){Yavg=Yavg/n;}


	if(f>0 && ncase>0)
	{
		for(j=0;j<=f;j++)
		{
			betaOld[number2indi[j]]=betaNew[number2indi[j]];
    	}

      	for (i=0;i<n;i++)
       	{
			//exponent=0;
			for(j=0;j<=f;j++)
			{
				Xt[j][i]=X[i][j];
			}
		}

		if(1)
		{
		//new XtX

		if (sexcov==0 && withcov==0 && N==2)
		{
			if((test==5 || (test==3 && alt==1)) && x1==1 && x2==1 && x1x2==1)
		    {
				A[0][0] = n;

				A[0][1] =  casecounts5[0][0][0][0][0] + casecounts5[0][1][0][0][0] + casecounts5[0][2][0][0][0]
				- casecounts5[2][0][0][0][0] - casecounts5[2][1][0][0][0] - casecounts5[2][2][0][0][0];

				A[0][2]=  casecounts5[0][0][0][0][0] + casecounts5[1][0][0][0][0] + casecounts5[2][0][0][0][0]
				- casecounts5[0][2][0][0][0] - casecounts5[1][2][0][0][0] - casecounts5[2][2][0][0][0];


				A[1][2]=  +casecounts5[0][0][0][0][0] - casecounts5[2][0][0][0][0]
							  - casecounts5[0][2][0][0][0] + casecounts5[2][2][0][0][0];

				A[1][0]=A[0][1];

				A[1][1] = casecounts5[0][0][0][0][0] + casecounts5[0][1][0][0][0] + casecounts5[0][2][0][0][0]
						  + casecounts5[2][0][0][0][0] + casecounts5[2][1][0][0][0] + casecounts5[2][2][0][0][0];

				A[2][0]=A[0][2];

				A[2][1]=A[1][2];

				A[2][2]= casecounts5[0][0][0][0][0] + casecounts5[1][0][0][0][0] + casecounts5[2][0][0][0][0] + casecounts5[0][2][0][0][0] +

				casecounts5[1][2][0][0][0] + casecounts5[2][2][0][0][0];

				if(alt==1)
				{
					A[0][3]=A[1][2];
					A[1][3]= casecounts5[0][0][0][0][0] - casecounts5[2][2][0][0][0] - casecounts5[0][2][0][0][0] + casecounts5[2][0][0][0][0];

					A[2][3]= +casecounts5[0][0][0][0][0] - casecounts5[2][2][0][0][0] + casecounts5[0][2][0][0][0] - casecounts5[2][0][0][0][0];

					A[3][0]=A[0][3];

					A[3][1]=A[1][3];

					A[3][2]= A[2][3];

					A[3][3]= casecounts5[0][0][0][0][0] + casecounts5[2][2][0][0][0] + casecounts5[0][2][0][0][0] + casecounts5[2][0][0][0][0];
				}
			} //test 3,5
			//else if (test==1 || test==2 || test==6)
			//{
			else if((test==6 || (test==4 && alt==1)) && x1==1 && x2==1 && x1D==1 && x2D==1 && x1x2==1 && x1x2D==1 && x1Dx2==1 && x1Dx2D==1)
			{
				A[0][0]=n;
				//if(alt==1){cout << "A[0][0]" << A[0][0] << "\n";}

				A[0][1]= casecounts5[0][0][0][0][0] + casecounts5[0][1][0][0][0] + casecounts5[0][2][0][0][0] - casecounts5[2][0][0][0][0] -
				casecounts5[2][1][0][0][0] - casecounts5[2][2][0][0][0];
				//if(alt==1){cout << "A[0][1]" << A[0][1] << "\n";}

				A[0][2]= -(casecounts5[0][0][0][0][0]*0.5) - (casecounts5[0][1][0][0][0]*0.5) - (casecounts5[0][2][0][0][0]*0.5) +
				(casecounts5[1][0][0][0][0]*0.5) + (casecounts5[1][1][0][0][0]*0.5) + (casecounts5[1][2][0][0][0]*0.5) -
				(casecounts5[2][0][0][0][0]*0.5) - (casecounts5[2][1][0][0][0]*0.5) - (casecounts5[2][2][0][0][0]*0.5);
				//if(alt==1){cout << "A[0][2]=" << A[0][2] << "\n";}

				A[0][3]= casecounts5[0][0][0][0][0] - casecounts5[0][2][0][0][0] + casecounts5[1][0][0][0][0] - casecounts5[1][2][0][0][0] +
				casecounts5[2][0][0][0][0]   - casecounts5[2][2][0][0][0];
				//if(alt==1){cout << "A[0][3]=" << A[0][3] << "\n";}

				A[0][4]= -(casecounts5[0][0][0][0][0]*0.5) + (casecounts5[0][1][0][0][0]*0.5) - (casecounts5[0][2][0][0][0]*0.5) -
				(casecounts5[1][0][0][0][0]*0.5) + (casecounts5[1][1][0][0][0]*0.5) - (casecounts5[1][2][0][0][0]*0.5) -
				(casecounts5[2][0][0][0][0]*0.5) + (casecounts5[2][1][0][0][0]*0.5) - (casecounts5[2][2][0][0][0]*0.5);
				//if (alt==1){cout << "A[0][4]=" << A[0][4] << "\n";}

				A[0][5]= casecounts5[0][0][0][0][0] - casecounts5[0][2][0][0][0] - casecounts5[2][0][0][0][0] + casecounts5[2][2][0][0][0];
				//if(alt==1){cout << "A[0][5]=" << A[0][5] << "\n";}

				A[0][6]= -(casecounts5[0][0][0][0][0]*0.5) + (casecounts5[0][1][0][0][0]*0.5) - (casecounts5[0][2][0][0][0]*0.5) +
				(casecounts5[2][0][0][0][0]*0.5) - (casecounts5[2][1][0][0][0]*0.5) + (casecounts5[2][2][0][0][0]*0.5);
				//if(alt==1){cout << "A[0][6]=" << A[0][6] << "\n";}

				A[0][7]= -(casecounts5[0][0][0][0][0]*0.5) + (casecounts5[0][2][0][0][0]*0.5) + (casecounts5[1][0][0][0][0]*0.5) -
				(casecounts5[1][2][0][0][0]*0.5) - (casecounts5[2][0][0][0][0]*0.5) + (casecounts5[2][2][0][0][0]*0.5);
				//if(alt==1){cout << "A[0][7]=" << A[0][7] << "\n";}

				A[1][0]= A[0][1];
				//if(alt==1){cout << "A[1][0]=" << A[1][0] << "\n";}

				A[1][1]= casecounts5[0][0][0][0][0] + casecounts5[0][1][0][0][0] + casecounts5[0][2][0][0][0] +
				casecounts5[2][0][0][0][0] + casecounts5[2][1][0][0][0] + casecounts5[2][2][0][0][0];
				//if(alt==1){cout << "A[1][1]=" << A[1][1] << "\n";}

				A[1][2]= -(casecounts5[0][0][0][0][0]*0.5) - (casecounts5[0][1][0][0][0]*0.5) - (casecounts5[0][2][0][0][0]*0.5) +
				(casecounts5[2][0][0][0][0]*0.5) + (casecounts5[2][1][0][0][0]*0.5) + (casecounts5[2][2][0][0][0]*0.5);
				//if(alt==1){cout << "A[1][2]=" << A[1][2] << "\n";}

				A[1][3]= A[0][5];
				//if(alt==1){cout << "A[1][3]=" << A[1][3] << "\n";}
				A[1][4]= A[0][6];
				//if(alt==1){cout << "A[1][4]=" << A[1][4] << "\n";}

				A[1][5]= casecounts5[0][0][0][0][0] - casecounts5[0][2][0][0][0] + casecounts5[2][0][0][0][0] - casecounts5[2][2][0][0][0];
				//if(alt==1){cout << "A[1][5]=" << A[1][5] << "\n";}

				A[1][6]= -(casecounts5[0][0][0][0][0]*0.5) + (casecounts5[0][1][0][0][0]*0.5) - (casecounts5[0][2][0][0][0]*0.5) -
				(casecounts5[2][0][0][0][0]*0.5) + (casecounts5[2][1][0][0][0]*0.5) - (casecounts5[2][2][0][0][0]*0.5);
				//if(alt==1){cout << "A[1][6]=" << A[1][6] << "\n";}

				A[1][7]= -(casecounts5[0][0][0][0][0]*0.5) + (casecounts5[0][2][0][0][0]*0.5) + (casecounts5[2][0][0][0][0]*0.5) -
				(casecounts5[2][2][0][0][0]*0.5);
				//if(alt==1){cout << "A[1][7]=" << A[1][7] << "\n";}



				A[2][0]= A[0][2];
				//if(alt==1){cout << "A[2][0]=" << A[2][0] << "\n";}
				A[2][1]= A[1][2];
				//if(alt==1){cout << "A[2][1]=" << A[2][1] << "\n";}

				A[2][2]= (casecounts5[0][0][0][0][0]*0.25) + (casecounts5[0][1][0][0][0]*0.25) + (casecounts5[0][2][0][0][0]*0.25) +
				(casecounts5[1][0][0][0][0]*0.25) + (casecounts5[1][1][0][0][0]*0.25) + (casecounts5[1][2][0][0][0]*0.25) +
				(casecounts5[2][0][0][0][0]*0.25) + (casecounts5[2][1][0][0][0]*0.25) + (casecounts5[2][2][0][0][0]*0.25);
				//if(alt==1){cout << "A[2][2]=" << A[2][2] << "\n";}

				A[2][3]= A[0][7];
				//if(alt==1){cout << "A[2][3]=" << A[2][3] << "\n";}

				A[2][4]= (casecounts5[0][0][0][0][0]*0.25) - (casecounts5[0][1][0][0][0]*0.25) + (casecounts5[0][2][0][0][0]*0.25) -
				(casecounts5[1][0][0][0][0]*0.25) + (casecounts5[1][1][0][0][0]*0.25) - (casecounts5[1][2][0][0][0]*0.25) +
				(casecounts5[2][0][0][0][0]*0.25) - (casecounts5[2][1][0][0][0]*0.25) + (casecounts5[2][2][0][0][0]*0.25);
				//if(alt==1){cout << "A[2][4]=" << A[2][4] << "\n";}

				A[2][5]= A[1][7];
				//if(alt==1){cout << "A[2][5]=" << A[2][5] << "\n";}

				A[2][6]= (casecounts5[0][0][0][0][0]*0.25) - (casecounts5[0][1][0][0][0]*0.25) + (casecounts5[0][2][0][0][0]*0.25) -
				(casecounts5[2][0][0][0][0]*0.25) + (casecounts5[2][1][0][0][0]*0.25) - (casecounts5[2][2][0][0][0]*0.25);
				//if(alt==1){cout << "A[2][6]=" << A[2][6] << "\n";}

				A[2][7]= (casecounts5[0][0][0][0][0]*0.25) - (casecounts5[0][2][0][0][0]*0.25) + (casecounts5[1][0][0][0][0]*0.25) -
				(casecounts5[1][2][0][0][0]*0.25) + (casecounts5[2][0][0][0][0]*0.25) - (casecounts5[2][2][0][0][0]*0.25);
				//if(alt==1){cout << "A[2][7]=" << A[2][7] << "\n";}



				A[3][0]= A[0][3];
				//if(alt==1){cout << "A[3][0]=" << A[3][0] << "\n";}
				A[3][1]= A[1][3];
				//if(alt==1){cout << "A[3][1]=" << A[3][1] << "\n";}
				A[3][2]= A[2][3];
				//if(alt==1){cout << "A[3][2]=" << A[3][2] << "\n";}

				A[3][3]= casecounts5[0][0][0][0][0] + casecounts5[0][2][0][0][0] + casecounts5[1][0][0][0][0] +
				casecounts5[1][2][0][0][0] + casecounts5[2][0][0][0][0] + casecounts5[2][2][0][0][0];
				//if(alt==1){cout << "A[3][3]=" << A[3][3] << "\n";}

				A[3][4]= -(casecounts5[0][0][0][0][0]*0.5) + (casecounts5[0][2][0][0][0]*0.5) - (casecounts5[1][0][0][0][0]*0.5) +
				(casecounts5[1][2][0][0][0]*0.5) - (casecounts5[2][0][0][0][0]*0.5) + (casecounts5[2][2][0][0][0]*0.5);
				//if(alt==1){cout << "A[3][4]=" << A[3][4] << "\n";}

				A[3][5]= casecounts5[0][0][0][0][0] + casecounts5[0][2][0][0][0]- casecounts5[2][0][0][0][0] - casecounts5[2][2][0][0][0];
				//if(alt==1){cout << "A[3][5]=" << A[3][5] << "\n";}

				A[3][6]= A[1][7];
				//if(alt==1){cout << "A[3][6]=" << A[3][6] << "\n";}

				A[3][7]= -(casecounts5[0][0][0][0][0]*0.5) - (casecounts5[0][2][0][0][0]*0.5) + (casecounts5[1][0][0][0][0]*0.5) +
				(casecounts5[1][2][0][0][0]*0.5) - (casecounts5[2][0][0][0][0]*0.5) - (casecounts5[2][2][0][0][0]*0.5);
				//if(alt==1){cout << "A[3][7]=" << A[3][7] << "\n";}



				A[4][0]= A[0][4];
				//if(alt==1){cout << "A[4][0]=" << A[4][0] << "\n";}
				A[4][1]= A[1][4];
				//if(alt==1){cout << "A[4][1]=" << A[4][1] << "\n";}
				A[4][2]= A[2][4];
				//if(alt==1){cout << "A[4][2]=" << A[4][2] << "\n";}
				A[4][3]= A[3][4];
				//if(alt==1){cout << "A[4][3]=" << A[4][3] << "\n";}
				A[4][4]= A[2][2];
				//if(alt==1){cout << "A[4][4]=" << A[4][4] << "\n";}
				A[4][5]= A[1][7];
				//if(alt==1){cout << "A[4][5]=" << A[4][5] << "\n";}

				A[4][6]= (casecounts5[0][0][0][0][0]*0.25) + (casecounts5[0][1][0][0][0]*0.25) + (casecounts5[0][2][0][0][0]*0.25) -
				(casecounts5[2][0][0][0][0]*0.25) - (casecounts5[2][1][0][0][0]*0.25) - (casecounts5[2][2][0][0][0]*0.25);
				//if(alt==1){cout << "A[4][6]=" << A[4][6] << "\n";}

				A[4][7]= (casecounts5[0][0][0][0][0]*0.25) - (casecounts5[0][2][0][0][0]*0.25) - (casecounts5[1][0][0][0][0]*0.25) +
				(casecounts5[1][2][0][0][0]*0.25) + (casecounts5[2][0][0][0][0]*0.25) - (casecounts5[2][2][0][0][0]*0.25);
				//if(alt==1){cout << "A[4][7]=" << A[4][7] << "\n";}



				A[5][0]= A[0][5];
				//if(alt==1){cout << "A[5][0]=" << A[5][0] << "\n";}
				A[5][1]= A[1][5];
				//if(alt==1){cout << "A[5][1]=" << A[5][1] << "\n";}
				A[5][2]= A[2][5];
				//if(alt==1){cout << "A[5][2]=" << A[5][2] << "\n";}
				A[5][3]= A[3][5];
				//if(alt==1){cout << "A[5][3]=" << A[5][3] << "\n";}
				A[5][4]= A[4][5];
				//if(alt==1){cout << "A[5][4]=" << A[5][4] << "\n";}

				A[5][5]= casecounts5[0][0][0][0][0] + casecounts5[0][2][0][0][0] + casecounts5[2][0][0][0][0] + casecounts5[2][2][0][0][0];
				//if(alt==1){cout << "A[5][5]=" << A[5][5] << "\n";}

				A[5][6]= -(casecounts5[0][0][0][0][0]*0.5) + (casecounts5[0][2][0][0][0]*0.5) - (casecounts5[2][0][0][0][0]*0.5) +
				(casecounts5[2][2][0][0][0]*0.5);
				//if(alt==1){cout << "A[5][6]=" << A[5][6] << "\n";}

				A[5][7]= -(casecounts5[0][0][0][0][0]*0.5) - (casecounts5[0][2][0][0][0]*0.5) + (casecounts5[2][0][0][0][0]*0.5) +
				(casecounts5[2][2][0][0][0]*0.5);
				//if(alt==1){cout << "A[5][7]=" << A[5][7] << "\n";}

				A[6][0]= A[0][6];
				//if(alt==1){cout << "A[6][0]=" << A[6][0] << "\n";}
				A[6][1]= A[1][6];
				//if(alt==1){cout << "A[6][1]=" << A[6][1] << "\n";}
				A[6][2]= A[2][6];
				//if(alt==1){cout << "A[6][2]=" << A[6][2] << "\n";}
				A[6][3]= A[3][6];
				//if(alt==1){cout << "A[6][3]=" << A[6][3] << "\n";}
				A[6][4]= A[4][6];
				//if(alt==1){cout << "A[6][4]=" << A[6][4] << "\n";}
				A[6][5]= A[5][6];
				//if(alt==1){cout << "A[6][5]=" << A[6][5] << "\n";}

				A[6][6]= (casecounts5[0][0][0][0][0]*0.25) + (casecounts5[0][1][0][0][0]*0.25) + (casecounts5[0][2][0][0][0]*0.25) +
				(casecounts5[2][0][0][0][0]*0.25) + (casecounts5[2][1][0][0][0]*0.25) + (casecounts5[2][2][0][0][0]*0.25);
				//if(alt==1){cout << "A[6][6]=" << A[6][6] << "\n";}

				A[6][7]= (casecounts5[0][0][0][0][0]*0.25) - (casecounts5[0][2][0][0][0]*0.25) - (casecounts5[2][0][0][0][0]*0.25) +
				(casecounts5[2][2][0][0][0]*0.25);
				//if(alt==1){cout << "A[6][7]=" << A[6][7] << "\n";}



				A[7][0]= A[0][7];
				//if(alt==1){cout << "A[7][0]=" << A[7][0] << "\n";}
				A[7][1]= A[1][7];
				//if(alt==1){cout << "A[7][1]=" << A[7][1] << "\n";}
				A[7][2]= A[2][7];
				//if(alt==1){cout << "A[7][2]=" << A[7][2] << "\n";}
				A[7][3]= A[3][7];
				//if(alt==1){cout << "A[7][3]=" << A[7][3] << "\n";}
				A[7][4]= A[4][7];
				//if(alt==1){cout << "A[7][4]=" << A[7][4] << "\n";}
				A[7][5]= A[5][7];
				//if(alt==1){cout << "A[7][5]=" << A[7][5] << "\n";}
				A[7][6]= A[6][7];
				//if(alt==1){cout << "A[7][6]=" << A[7][6] << "\n";}

				A[7][7]= (casecounts5[0][0][0][0][0]*0.25) + (casecounts5[0][2][0][0][0]*0.25) + (casecounts5[1][0][0][0][0]*0.25) +
				(casecounts5[1][2][0][0][0]*0.25) + (casecounts5[2][0][0][0][0]*0.25) + (casecounts5[2][2][0][0][0]*0.25);
				//if(alt==1){cout << "A[7][7]=" << A[7][7] << "\n";}

				if(alt==1)
				{


					A[0][8]= (casecounts5[0][0][0][0][0]*0.25) - (casecounts5[0][1][0][0][0]*0.25) + (casecounts5[0][2][0][0][0]*0.25) -
					(casecounts5[1][0][0][0][0]*0.25) + (casecounts5[1][1][0][0][0]*0.25) - (casecounts5[1][2][0][0][0]*0.25) +
					(casecounts5[2][0][0][0][0]*0.25) - (casecounts5[2][1][0][0][0]*0.25) + (casecounts5[2][2][0][0][0]*0.25);
					//if(alt==1){cout << "A[0][8]=" << A[0][8] << "\n";}

					A[1][8]= (casecounts5[0][0][0][0][0]*0.25) - (casecounts5[0][1][0][0][0]*0.25) + (casecounts5[0][2][0][0][0]*0.25) -
					(casecounts5[2][0][0][0][0]*0.25) + (casecounts5[2][1][0][0][0]*0.25) - (casecounts5[2][2][0][0][0]*0.25);
					//if(alt==1){cout << "A[1][8]=" << A[1][8] << "\n";}

					A[2][8]= -(casecounts5[0][0][0][0][0]*0.125) + (casecounts5[0][1][0][0][0]*0.125) - (casecounts5[0][2][0][0][0]*0.125) -
					(casecounts5[1][0][0][0][0]*0.125) + (casecounts5[1][1][0][0][0]*0.125) - (casecounts5[1][2][0][0][0]*0.125) -
					(casecounts5[2][0][0][0][0]*0.125) + (casecounts5[2][1][0][0][0]*0.125) - (casecounts5[2][2][0][0][0]*0.125);
					//if(alt==1){cout << "A[2][8]=" << A[2][8] << "\n";}

					A[3][8]= (casecounts5[0][0][0][0][0]*0.25) - (casecounts5[0][2][0][0][0]*0.25) - (casecounts5[1][0][0][0][0]*0.25) +
					(casecounts5[1][2][0][0][0]*0.25) + (casecounts5[2][0][0][0][0]*0.25) - (casecounts5[2][2][0][0][0]*0.25);
					//if(alt==1){cout << "A[3][8]=" << A[3][8] << "\n";}

					A[4][8]= -(casecounts5[0][0][0][0][0]*0.125) - (casecounts5[0][1][0][0][0]*0.125) - (casecounts5[0][2][0][0][0]*0.125) +
					(casecounts5[1][0][0][0][0]*0.125) + (casecounts5[1][1][0][0][0]*0.125) + (casecounts5[1][2][0][0][0]*0.125) -
					(casecounts5[2][0][0][0][0]*0.125) - (casecounts5[2][1][0][0][0]*0.125) - (casecounts5[2][2][0][0][0]*0.125);
					//if(alt==1){cout << "A[4][8]=" << A[4][8] << "\n";}

					A[5][8]= (casecounts5[0][0][0][0][0]*0.25) - (casecounts5[0][2][0][0][0]*0.25) - (casecounts5[2][0][0][0][0]*0.25) +
					(casecounts5[2][2][0][0][0]*0.25);
					//if(alt==1){cout << "A[5][8]=" << A[5][8] << "\n";}

					A[6][8]= -(casecounts5[0][0][0][0][0]*0.125) - (casecounts5[0][1][0][0][0]*0.125) - (casecounts5[0][2][0][0][0]*0.125) +
					(casecounts5[2][0][0][0][0]*0.125) + (casecounts5[2][1][0][0][0]*0.125) + (casecounts5[2][2][0][0][0]*0.125);
					//if(alt==1){cout << "A[6][8]=" << A[6][8] << "\n";}

					A[7][8]= -(casecounts5[0][0][0][0][0]*0.125) + (casecounts5[0][2][0][0][0]*0.125) - (casecounts5[1][0][0][0][0]*0.125) +
					(casecounts5[1][2][0][0][0]*0.125) - (casecounts5[2][0][0][0][0]*0.125) + (casecounts5[2][2][0][0][0]*0.125);
					//if(alt==1){cout << "A[7][8]=" << A[7][8] << "\n";}

					A[8][0]= A[0][8];
					//if(alt==1){cout << "A[8][0]=" << A[8][0] << "\n";}
					A[8][1]= A[1][8];
					//if(alt==1){cout << "A[8][1]=" << A[8][1] << "\n";}
					A[8][2]= A[2][8];
					//if(alt==1){cout << "A[8][2]=" << A[8][2] << "\n";}
					A[8][3]= A[3][8];
					//if(alt==1){cout << "A[8][3]=" << A[8][3] << "\n";}
					A[8][4]= A[4][8];
					//if(alt==1){cout << "A[8][4]=" << A[8][4] << "\n";}
					A[8][5]= A[5][8];
					//if(alt==1){cout << "A[8][5]=" << A[8][5] << "\n";}
					A[8][6]= A[6][8];
					//if(alt==1){cout << "A[8][6]=" << A[8][6] << "\n";}
					A[8][7]= A[7][8];
					//if(alt==1){cout << "A[8][7]=" << A[8][7] << "\n";}

					A[8][8]= (casecounts5[0][0][0][0][0]*0.0625) + (casecounts5[0][1][0][0][0]*0.0625) + (casecounts5[0][2][0][0][0]*0.0625) +
					(casecounts5[1][0][0][0][0]*0.0625) + (casecounts5[1][1][0][0][0]*0.0625) + (casecounts5[1][2][0][0][0]*0.0625) +
					(casecounts5[2][0][0][0][0]*0.0625) + (casecounts5[2][1][0][0][0]*0.0625) + (casecounts5[2][2][0][0][0]*0.0625);
					//if(alt==1){cout << "A[8][8]=" << A[8][8] << "\n";}
				}
			}
			//} //test 1,2
            else
            {
				dummy=MatrixMultSym(Xt,f+1,n,X,n,f+1,A);
            }
		}
		else //sexcov & Co
		{
			dummy=MatrixMultSym(Xt,f+1,n,X,n,f+1,A);
		}

      } // alt !=0
		// Invertieren von A


		// Dwyer-Algorithmus
		//dummy=DwyerInv(f+1, A, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, Ainv);
		dummy=0;

		if (dummy == 0)
		{
			for(i=0;i<=f;i++)
			{
				S[i]=0;
				for(j=0;j<=f;j++)
				{
					VNN[i][j]=0;
					Sinv[i][j]=0;
				}
			}
			dummy=mysvd(f+1, f+1, &A, &VNN, &S, f+1, f+1);

			for(j=0;j<=f;j++)
			{
				for (i=0;i<=f;i++)
				{
					UNNT[j][i]=A[i][j];
				}
			}

			newf=-1;
			for (i=0;i<=f;i++)
			  {
			    if(fabs(S[i])>1.e-6){Sinv[i][i]=1/S[i];newf++;}
			    else {dfUnchanged=0;}
			  }

			/*dummy=MatrixMult(VNN, f+1, f+1,Sinv, f+1, f+1, &A0);
			  dummy=MatrixMult(A0, f+1, f+1,UNNT, f+1, f+1, &Ainv);*/
			dummy=MatrixMultMod(VNN, f+1, f+1,Sinv, f+1, f+1, A0);
			dummy=MatrixMultMod(A0, f+1, f+1,UNNT, f+1, f+1, Ainv);
		}
		else
		{
				newf=f;
		}

		//from here different from logreg
		//rss = Y'Y - Y'XAinvX'Y berechnen
		for (i=0;i<n;i++)
    	{
    	     Yt[0][i]=Y[i];
		}

		for (i=0;i<n;i++)
    	{
			YtY+=Yt[0][i]*Y[i];
		}

		/*dummy=MatrixMult(Yt,1,n,X,n,f+1,&YtX);
		dummy=MatrixMult(YtX, 1, f+1,Ainv, f+1, f+1, &YtXAinv);
		dummy=MatrixMult(YtXAinv, 1, f+1,Xt,f+1,n, &YtXAinvXt);*/

		dummy=MatrixMultMod(Yt,1,n,X,n,f+1,YtX);
		dummy=MatrixMultMod(YtX, 1, f+1,Ainv, f+1, f+1, YtXAinv);
		dummy=MatrixMultMod(YtXAinv, 1, f+1,Xt,f+1,n, YtXAinvXt);

		for (i=0;i<n;i++)
    	{
			YtXAinvXtY+=YtXAinvXt[0][i]*Y[i];
			ssy+=(Y[i]-Yavg)*(Y[i]-Yavg);
		}

		rss =YtY-YtXAinvXtY;
	} // f>0
	else
	{
		rss=0;
		f=0;
		logNew=0;
		newf=0;

		Yavg=0;
		for(i=0;i<n;i++)
		{
			Yavg+=Y[i];
		}

		if(n>0){Yavg=Yavg/n;}
		for(i=0;i<n;i++)
		{
			rss+=(Y[i]-Yavg)*(Y[i]-Yavg);
			ssy+=(Y[i]-Yavg)*(Y[i]-Yavg);
		}
		//NEW
		if(nMc==0)
		{
			result.betaNew_se[1] = -1;
			result.betaNew_lcl[1] = -1;
			result.betaNew_rcl[1] = -1;
		}

		result.b[1] = 0;
	}


	//moved up in version 682
	for(i=0;i<f+1;i++)
	{
	    //NEW
		betaNew[number2indi[i]]=0;
		if(newf>0){betaNew[number2indi[i]]=YtXAinv[0][i];}
		if(nMc==0)
		{
			result.betaNew_se[number2indi[i]] = -1;
			result.betaNew_lcl[number2indi[i]] = -1;
			result.betaNew_rcl[number2indi[i]] = -1;

		}
		result.b[number2indi[i]] = -1;
		//printf("here2\n");
	}

	if(newf>0 && dfUnchanged)
	{
		for(i=0;i<f+1;i++)
		{
			//NEW
			if (n1-newf-1 > 0 && f > 0 && nMc==0)
			{
				if((Ainv[i][i]*rss/(n1-newf-1))<0)
				{
					result.betaNew_se[number2indi[i]] = -1;
					result.betaNew_lcl[number2indi[i]] = -1;
					result.betaNew_rcl[number2indi[i]] = -1;
				}
				else
				{



					result.betaNew_se[number2indi[i]]=sqrt(Ainv[i][i]*rss/(n1-newf-1));
					result.betaNew_lcl[number2indi[i]]= betaNew[number2indi[i]]-result.betaNew_se[number2indi[i]]*tquantile((int)(n1-newf-1));
					result.betaNew_rcl[number2indi[i]]= betaNew[number2indi[i]]+result.betaNew_se[number2indi[i]]*tquantile((int)(n1-newf-1));
				}
			}
		}

	}
	else if (!dfUnchanged)
	  {
	    for(i=0;i<f+1;i++)
	      {
		result.betaNew_se[number2indi[i]] = -2;
		result.betaNew_lcl[number2indi[i]] = -2;
		result.betaNew_rcl[number2indi[i]] = -2;
	      }
	  }
	else
	{
		betaNew[1]=0;

		if(nMc==0)
		{
			result.betaNew_se[1] = -1;
			result.betaNew_lcl[1] = -1;
			result.betaNew_rcl[1] = -1;
		}
		result.b[1] = -1;
	}


	result.df=n1-newf-1;

	for(i=0;i<f+1;i++)
	  {
	    if(dfUnchanged) {result.b[number2indi[i]]=betaNew[number2indi[i]];}
	    else {result.b[number2indi[i]]=-2;}
	  }

	result.bcvsex=betaNew[27];  //covariate sex


	if(haplo && alt==1)
	{
		for(i=0;i<8;i++)
		{
			result.h[i]=hapinfo.h[i];
		}
	}


	for(i=0;i<maxIndexCov;i++)
	{
		result.bcv[i]=betaNew[i+27]; //covariate parameters
	}
	result.sc=rss;
	result.rsquare=1-(rss/ssy); // calculate rsquare

    //exit(1);
	//exit(1);

	if(alt==1 && covariancematrix)  //SIGMA1
	  {
	    int k=0;

	    for(i=0;i<f+1;i++)
	      {
		for(j=i;j<f+1;j++)
		  {
		    if((n-f-1)!=0)
		      {
			if(dfUnchanged==1) {result.sigma1[k]=Ainv[i][j]*rss/(n-f-1);} // sigma1 matrix as 1 dim array
			else {result.sigma1[k]=-2;}
			//cout << "f=" << f << " Sigma1[" << k << "]=" << Ainv[i][j]*rss/(n-f-1) << endl;
			k=k+1;
		      }
		    else
		      {
			result.sigma1[k]=-99999;
		      }
		  }
	      }
	  }

	return result;
};


// Neue anova1df-Funktion
double anova1df(struct PERSON *person, int nlinestfam, int a, int b, int c, int thread, struct TRAITAVG *traitavg,
				int npplqc, int* PPLMap, uint64_t*** BinSNPs, struct PPLLOCATION* PplLocations)
{
	int i = 0;
	//int j = 0;
	double n=0;
	double n1_=0;
	double n2_=0;
	double n_1=0;
	double n_2=0;
	double n11=0;
	double n12=0;
	double n21=0;
	double n22=0;
	double x=0; //overall sum
	double xsq=0; //overall sum of squares
	double x1_=0; //sum allele 1 SNP 1
	double x2_=0; //sum allele 2 SNP 1
	double x_1=0; //sum allele 1 SNP 2
	double x_2=0; //sum allele 2 SNP 2
	double x11=0; //sum allele combi 1-1
	double x12=0; //sum allele combi 1-2
	double x21=0; //sum allele combi 2-1
	double x22=0; //sum allele combi 2-2
    double x11sq=0; //sum sq allele combi 1-1
	double x12sq=0; //sum sq allele combi 1-2
	double x21sq=0; //sum sq allele combi 2-1
	double x22sq=0; //sum sq allele combi 2-2

	double SAQ=0; //total
	double SAQa=0; //first factor
	double SAQb=0; //second factor
	double SAQin=0; //inside classes
	double SAQab=0; //interaction
	double F=0;
	//double check=0;


	initTRAITavg(traitavg);


	 struct PPLLOCATION guy;

	for (int iMod=0; iMod< npplqc; iMod++)
	{
		guy = PplLocations[iMod];
		i = PPLMap[iMod];

		//if (person[i].qcin == 1 && M[a].d2[i].d3 !=0 && M[b].d2[i].d3 != 0)
		if (person[i].qcin == 1 && !getbit64(BinSNPs[a][guy.nr][0],guy.pos) && !getbit64(BinSNPs[b][guy.nr][0],guy.pos))
		{
			n+=2;
			x+=2*person[i].qtaff[thread];
			xsq+=2*person[i].qtaff[thread]*person[i].qtaff[thread];

			//if (M[a].d2[i].d3 == 1 && M[b].d2[i].d3 == 1) // 11 11
			if (getbit64(BinSNPs[a][guy.nr][1],guy.pos) && getbit64(BinSNPs[b][guy.nr][1],guy.pos))
			{
			 x1_+=2*person[i].qtaff[thread];
			 x_1+=2*person[i].qtaff[thread];
			 x11+=2*person[i].qtaff[thread];n11+=2;
             x11sq+=2*(person[i].qtaff[thread]*person[i].qtaff[thread]);
			}
			//else if (M[a].d2[i].d3 == 1 && M[b].d2[i].d3 == 2) // 11 12
			else if (getbit64(BinSNPs[a][guy.nr][1],guy.pos)  && getbit64(BinSNPs[b][guy.nr][2],guy.pos)) // 11 12
			{
			 x1_+=2*person[i].qtaff[thread];
			 x_1+=person[i].qtaff[thread];
			 x_2+=person[i].qtaff[thread];
			 x11+=person[i].qtaff[thread];n11++;
             x12+=person[i].qtaff[thread];n12++;
             x11sq+=person[i].qtaff[thread]*person[i].qtaff[thread];
             x12sq+=person[i].qtaff[thread]*person[i].qtaff[thread];
			}
			//else if (M[a].d2[i].d3 == 1 && M[b].d2[i].d3 == 3) // 11 22
			else if (getbit64(BinSNPs[a][guy.nr][1],guy.pos) && getbit64(BinSNPs[b][guy.nr][3],guy.pos)) // 11 22
			{
			 x1_+=2*person[i].qtaff[thread];
			 x_2+=2*person[i].qtaff[thread];
             x12+=2*person[i].qtaff[thread];n12+=2;
             x12sq+=2*person[i].qtaff[thread]*person[i].qtaff[thread];
			}
			//else if (M[a].d2[i].d3 == 2 && M[b].d2[i].d3 == 1) // 12 11
			else if (getbit64(BinSNPs[a][guy.nr][2],guy.pos) && getbit64(BinSNPs[b][guy.nr][1],guy.pos)) // 12 11
			{
			 x1_+=person[i].qtaff[thread];
			 x2_+=person[i].qtaff[thread];
			 x_1+=2*person[i].qtaff[thread];
             x11+=person[i].qtaff[thread];n11++;
			 x21+=person[i].qtaff[thread];n21++;
			 x11sq+=person[i].qtaff[thread]*person[i].qtaff[thread];
			 x21sq+=person[i].qtaff[thread]*person[i].qtaff[thread];
			}
			//else if (M[a].d2[i].d3 == 2 && M[b].d2[i].d3 == 2) // 12 12
			else if (getbit64(BinSNPs[a][guy.nr][2],guy.pos) && getbit64(BinSNPs[b][guy.nr][2],guy.pos)) // 12 12
			{
			 x1_+=person[i].qtaff[thread];
			 x2_+=person[i].qtaff[thread];
			 x_1+=person[i].qtaff[thread];
             x_2+=person[i].qtaff[thread];
             x11+=0.5*person[i].qtaff[thread];n11+=0.5;
			 x12+=0.5*person[i].qtaff[thread];n12+=0.5;
             x21+=0.5*person[i].qtaff[thread];n21+=0.5;
             x22+=0.5*person[i].qtaff[thread];n22+=0.5;
             x11sq+=0.5*person[i].qtaff[thread]*person[i].qtaff[thread];
			 x12sq+=0.5*person[i].qtaff[thread]*person[i].qtaff[thread];
             x21sq+=0.5*person[i].qtaff[thread]*person[i].qtaff[thread];
             x22sq+=0.5*person[i].qtaff[thread]*person[i].qtaff[thread];
			}
			//else if (M[a].d2[i].d3 == 2 && M[b].d2[i].d3 == 3) // 12 22
			else if (getbit64(BinSNPs[a][guy.nr][2],guy.pos) && getbit64(BinSNPs[b][guy.nr][3],guy.pos)) // 12 22
			{
			 x1_+=person[i].qtaff[thread];
			 x2_+=person[i].qtaff[thread];
             x_2+=2*person[i].qtaff[thread];
			 x12+=person[i].qtaff[thread];n12++;
             x22+=person[i].qtaff[thread];n22++;
             x12sq+=person[i].qtaff[thread]*person[i].qtaff[thread];
             x22sq+=person[i].qtaff[thread]*person[i].qtaff[thread];
			}
			//else if (M[a].d2[i].d3 == 3 && M[b].d2[i].d3 == 1) // 22 11
			else if (getbit64(BinSNPs[a][guy.nr][3],guy.pos) && getbit64(BinSNPs[b][guy.nr][1],guy.pos)) // 22 11
			{
			 x2_+=2*person[i].qtaff[thread];
             x_1+=2*person[i].qtaff[thread];
			 x21+=2*person[i].qtaff[thread];n21+=2;
             x21sq+=2*person[i].qtaff[thread]*person[i].qtaff[thread];
			}
			//else if (M[a].d2[i].d3 == 3 && M[b].d2[i].d3 == 2) // 22 12
			else if (getbit64(BinSNPs[a][guy.nr][3],guy.pos) && getbit64(BinSNPs[b][guy.nr][2],guy.pos)) // 22 12
			{
			 x2_+=2*person[i].qtaff[thread];
             x_1+=person[i].qtaff[thread];
             x_2+=person[i].qtaff[thread];
			 x21+=person[i].qtaff[thread];n21+=1;

             x22+=person[i].qtaff[thread];n22+=1;
             x21sq+=person[i].qtaff[thread]*person[i].qtaff[thread];
             x22sq+=person[i].qtaff[thread]*person[i].qtaff[thread];
			}
			//else if (M[a].d2[i].d3 == 3 && M[b].d2[i].d3 == 3) // 22 22
			else if (getbit64(BinSNPs[a][guy.nr][3],guy.pos) && getbit64(BinSNPs[b][guy.nr][3],guy.pos)) // 22 22
			{
			 x2_+=2*person[i].qtaff[thread];
             x_2+=2*person[i].qtaff[thread];
			 x22+=2*person[i].qtaff[thread];n22+=2;
			 x22sq+=2*person[i].qtaff[thread]*person[i].qtaff[thread];
			}
			else
			{
			 printf("error anova1df\n");exit(1);
			}
		} //person in
	} // loop person

	if(n<=4)
	  {
	   return 1;
	  }
	else
      {
	   SAQ=xsq-(x*x/n);
	   //SAQin=x11sq-(x11*x11/n11)+x12sq-(x12*x12/n12)+x21sq-(x21*x21/n21)+x22sq-(x22*x22/n22);

	   SAQin=0;
	   if(n11>0){SAQin=SAQin+x11sq-(x11*x11/n11);}
	   if(n12>0){SAQin=SAQin+x12sq-(x12*x12/n12);}
	   if(n21>0){SAQin=SAQin+x21sq-(x21*x21/n21);}
       if(n22>0){SAQin=SAQin+x22sq-(x22*x22/n22);}

	   //cout << "SAQin " << SAQin << "\n";

	   n1_=n11+n12;n2_=n21+n22;
	   n_1=n11+n21;n_2=n12+n22;
	   x=x/n;
	   if(n1_>0){x1_=x1_/n1_;(*traitavg).snp1allele[0]=x1_;}
	   if(n2_>0){x2_=x2_/n2_;(*traitavg).snp1allele[1]=x2_;}
	   if(n_1>0){x_1=x_1/n_1;(*traitavg).snp2allele[0]=x_1;}
	   if(n_2>0){x_2=x_2/n_2;(*traitavg).snp2allele[1]=x_2;}
	   if(n11>0){x11=x11/n11;(*traitavg).allele[0][0]=x11;}
	   if(n12>0){x12=x12/n12;(*traitavg).allele[0][1]=x12;}
	   if(n21>0){x21=x21/n21;(*traitavg).allele[1][0]=x21;}
	   if(n22>0){x22=x22/n22;(*traitavg).allele[1][1]=x22;}

	   SAQa=n1_*(x1_-x)*(x1_-x)+n2_*(x2_-x)*(x2_-x);
	   SAQb=n_1*(x_1-x)*(x_1-x)+n_2*(x_2-x)*(x_2-x);



	   SAQab=SAQ-SAQa-SAQb-SAQin;

	   F=(SAQab/SAQin)*(n-4);
	   if(F<=0){return 1;}

	   //cout << "Anova " << betai((n-4)/2,0.5,(n-4)/((n-4)+1*F)) << "\n";
	   return betai((n-4)/2,0.5,(n-4)/((n-4)+1*F));
      }


} //end anova1df


// Neue anova4df-Funktion
double anova4df(struct PERSON *person, int nlinestfam, int a, int b, int c, int thread, struct TRAITAVG *traitavg,
				int npplqc, int* PPLMap, uint64_t*** BinSNPs, struct PPLLOCATION* PplLocations)
{
	int i = 0;
	int j = 0;
	double n=0;
	double n1[3]; //genotypes first SNP
	double n2[3]; //genotypes second SNP
	double N[3][3];// genotypes
	double x1[3]; //values first SNP
	double x2[3]; //values second SNP
	double X[3][3];// values
	double Xsq[3][3];// values
	double x=0; //overall sum
	double xsq=0; //overall sum of squares
	double SAQ=0; //total
	double SAQa=0; //first factor
	double SAQb=0; //second factor
	double SAQin=0; //inside classes
	double SAQab=0; //interaction
	double F=0;
	//double check=0;

	initTRAITavg(traitavg);

	for(i=0;i<3;i++)
	   {
	    n1[i]=0;n2[i]=0;x1[i]=0;x2[i]=0;
		for(j=0;j<3;j++)
		   {
			N[i][j]=0;X[i][j]=0;Xsq[i][j]=0;
		   }
	   }

	struct PPLLOCATION guy;

	for (int iMod=0; iMod< npplqc; iMod++)
	{
		guy = PplLocations[iMod];
		i = PPLMap[iMod];


		//if (person[i].qcin == 1 && M[a].d2[i].d3 !=0 && M[b].d2[i].d3 != 0)
		if (person[i].qcin == 1 && !getbit64(BinSNPs[a][guy.nr][0],guy.pos) && !getbit64(BinSNPs[b][guy.nr][0],guy.pos))
		{
			n+=1;
			x+=person[i].qtaff[thread];
			xsq+=person[i].qtaff[thread]*person[i].qtaff[thread];

			//if (M[a].d2[i].d3 == 1 && M[b].d2[i].d3 == 1) // 11 11
			if (getbit64(BinSNPs[a][guy.nr][1],guy.pos) && getbit64(BinSNPs[b][guy.nr][1],guy.pos)) // 11 11
			{
			 x1[0]+=person[i].qtaff[thread];
			 x2[0]+=person[i].qtaff[thread];
			 X[0][0]+=person[i].qtaff[thread];N[0][0]++;
             Xsq[0][0]+=person[i].qtaff[thread]*person[i].qtaff[thread];
			}
			//else if (M[a].d2[i].d3 == 1 && M[b].d2[i].d3 == 2) // 11 12
			else if (getbit64(BinSNPs[a][guy.nr][1],guy.pos) && getbit64(BinSNPs[b][guy.nr][2],guy.pos)) // 11 12
			{
			 x1[0]+=person[i].qtaff[thread];
			 x2[1]+=person[i].qtaff[thread];
			 X[0][1]+=person[i].qtaff[thread];N[0][1]++;
			 Xsq[0][1]+=person[i].qtaff[thread]*person[i].qtaff[thread];
			}
			//else if (M[a].d2[i].d3 == 1 && M[b].d2[i].d3 == 3) // 11 22
			else if (getbit64(BinSNPs[a][guy.nr][1],guy.pos) && getbit64(BinSNPs[b][guy.nr][3],guy.pos)) // 11 22
			{
			 x1[0]+=person[i].qtaff[thread];
			 x2[2]+=person[i].qtaff[thread];
			 X[0][2]+=person[i].qtaff[thread];N[0][2]++;
			 Xsq[0][2]+=person[i].qtaff[thread]*person[i].qtaff[thread];
			}
			//else if (M[a].d2[i].d3 == 2 && M[b].d2[i].d3 == 1) // 12 11
			else if (getbit64(BinSNPs[a][guy.nr][2],guy.pos) && getbit64(BinSNPs[b][guy.nr][1],guy.pos)) // 12 11
			{
			 x1[1]+=person[i].qtaff[thread];
			 x2[0]+=person[i].qtaff[thread];
			 X[1][0]+=person[i].qtaff[thread];N[1][0]++;
			 Xsq[1][0]+=person[i].qtaff[thread]*person[i].qtaff[thread];
			}
			//else if (M[a].d2[i].d3 == 2 && M[b].d2[i].d3 == 2) // 12 12
			else if (getbit64(BinSNPs[a][guy.nr][2],guy.pos) && getbit64(BinSNPs[b][guy.nr][2],guy.pos)) // 12 12
			{
			 x1[1]+=person[i].qtaff[thread];
			 x2[1]+=person[i].qtaff[thread];
			 X[1][1]+=person[i].qtaff[thread];N[1][1]++;
			 Xsq[1][1]+=person[i].qtaff[thread]*person[i].qtaff[thread];
			}
			//else if (M[a].d2[i].d3 == 2 && M[b].d2[i].d3 == 3) // 12 22
			else if (getbit64(BinSNPs[a][guy.nr][2],guy.pos) && getbit64(BinSNPs[b][guy.nr][3],guy.pos)) // 12 22
			{
			 x1[1]+=person[i].qtaff[thread];
			 x2[2]+=person[i].qtaff[thread];
			 X[1][2]+=person[i].qtaff[thread];N[1][2]++;
			 Xsq[1][2]+=person[i].qtaff[thread]*person[i].qtaff[thread];
			}
			//else if (M[a].d2[i].d3 == 3 && M[b].d2[i].d3 == 1) // 22 11
			else if (getbit64(BinSNPs[a][guy.nr][3],guy.pos) && getbit64(BinSNPs[b][guy.nr][1],guy.pos)) // 22 11
			{
			 x1[2]+=person[i].qtaff[thread];
			 x2[0]+=person[i].qtaff[thread];
			 X[2][0]+=person[i].qtaff[thread];N[2][0]++;
			 Xsq[2][0]+=person[i].qtaff[thread]*person[i].qtaff[thread];
			}
			//else if (M[a].d2[i].d3 == 3 && M[b].d2[i].d3 == 2) // 22 12
			else if (getbit64(BinSNPs[a][guy.nr][3],guy.pos) && getbit64(BinSNPs[b][guy.nr][2],guy.pos)) // 22 12
			{
			 x1[2]+=person[i].qtaff[thread];
			 x2[1]+=person[i].qtaff[thread];
			 X[2][1]+=person[i].qtaff[thread];N[2][1]++;
			 Xsq[2][1]+=person[i].qtaff[thread]*person[i].qtaff[thread];
			}

			else if (getbit64(BinSNPs[a][guy.nr][3],guy.pos) && getbit64(BinSNPs[b][guy.nr][3],guy.pos)) // 22 22
			{
			 x1[2]+=person[i].qtaff[thread];
			 x2[2]+=person[i].qtaff[thread];
			 X[2][2]+=person[i].qtaff[thread];N[2][2]++;
			 Xsq[2][2]+=person[i].qtaff[thread]*person[i].qtaff[thread];
			}
			else
			{
			 printf("error anova1df\n");exit(1);
			}
		} //person in
	} // loop person

	if(n<=8)
	  {
	   return 1;
	  }
	else
      {
	   SAQ=xsq-(x*x/n);

	   SAQin=0;

	   for(i=0;i<3;i++)
	   {
		for(j=0;j<3;j++)
		   {
		    n1[i]+=N[i][j];
		    n2[j]+=N[i][j];
			if(N[i][j]>0){SAQin+=Xsq[i][j]-(X[i][j]*X[i][j]/N[i][j]);}
		   }
	   }

	   x=x/n;


	   for(i=0;i<3;i++)
	   {
        if(n1[i]>0){x1[i]=x1[i]/n1[i];(*traitavg).snp1geno[i]=x1[i];}

         if(n2[i]>0){x2[i]=x2[i]/n2[i];(*traitavg).snp2geno[i]=x2[i];}

		for(j=0;j<3;j++)
		   {
		    if(N[i][j]>0){X[i][j]=X[i][j]/N[i][j];(*traitavg).geno[i][j]=X[i][j];}

		   }
	   }

       for(i=0;i<3;i++)
	   {
        SAQa+=n1[i]*(x1[i]-x)*(x1[i]-x);
		SAQb+=n2[i]*(x2[i]-x)*(x2[i]-x);
	   }

	   SAQab=SAQ-SAQa-SAQb-SAQin;

	   F=(SAQab/SAQin)*(n-8)/4;

	   if(F<=0){return 1;}

	   return betai((n-8)/2,2,(n-8)/((n-8)+4*F));
      }


} //end anova4df


//PWT: Bewertungsfunktion aufrufen ratio pro pathway berechnen (pw.singletop/pw.counts)
double snpratio(struct PATHWAY pathway, int thread)
{
	if (pathway.counts == 0)
	{
		return 0;
	}
	else
	{
	    //printf("inside %f\n",double(pathway.singletop[thread])/double(pathway.counts));
		return (double(pathway.singletop[thread])/double(pathway.counts));
	}
}


// interactionRatio // Pathwaytest 6
double interactionratio(struct PATHWAY pathway, int n, int thread)
{
	if (pathway.counts == 0)
	{
		return 0;
	}
	else
	{
	    if (n == 0)
		{
			return (double(pathway.intertop[thread])/(double(pathway.counts)*(double(pathway.counts)-1)/2));
		}
		else
		{
			return (double(pathway.intertopmod[thread])/(double(pathway.counts)*(double(pathway.counts)-1)/2));
		}
	}
}


//PWT_TIM
double generatio(struct PATHWAY pathway, int thread)
{
	if (pathway.counts == 0 || pathway.ngenes == 0)
	{
		return 0;
	}
	else
	{
		return (double(pathway.singletopmod[thread])/double(pathway.ngenes));
	}
}


//PWT_TIM
double gammp(double a, double x)
{
	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0)
	{
		printf("Invalid arguments in routine gammp");exit(1);
	}
	if (x < (a+1))
	{
		gser(&gamser,a,x,&gln);
		return gamser;
	}
	else
	{
		gcf(&gammcf,a,x,&gln);
		return 1-gammcf;
	}
}


double invgammap(double p, double df)
{
	int j;
	double gln, x,err,t,u,pp,lna1,afac,a1,min;
	double eps=1.e-8;

	a1=df-1;
	gln=lnpValueCalc(df);

	if (df <= 0.0)
	{
		printf("df must be pos in invgammap.\n"); exit(1);
	}
	if (p >= 1)
	{
		return df + 100*sqrt(df);
	}
	if (p <= 0) {return 0;}
	if (df > 1)
	{
		lna1=log(a1);
		afac = exp(a1*(lna1-1)-gln);
		pp = (p < 0.5)? p : 1 - p;
		t = sqrt(-2*log(pp));
		x = (2.30753+t*0.27061)/(1+t*(0.99229+t*0.04481)) - t;
		if (p < 0.5) {x = -x;}

		x = df*pow(1-1/(9*df)-x/(3*sqrt(df)),3);
		if(df <= 1.e-3){x=1.e-3;}
  	}
    else
    {
		t = 1.0 - df*(0.253+df*0.12);
		if (p < t) {x = pow(p/t,1/df);}
		else {x = 1-log(1-(p-t)/(1-t));}
	}
 	for (j=0;j<12;j++)
	{
		if (x <= 0.0) {return 0.0;}
		err = gammp(df,x) - p;
		if (df > 1) {t = afac*exp(-(x-a1)+a1*(log(x)-lna1));}
		else {t = exp(-x+a1*log(x)-gln);}
		u = err/t;

		min=1;
		if (( u*( (df-1)/x - 1)) < 1)
		{
			min=u*( (df-1)/x - 1);
		}

		x -= (t = u/(1-0.5*min));

		if (x <= 0) {x = 0.5*(x + t);}
		if (fabs(t) < eps*x ) {break;}
	} //end j
	return x;
}


double chiinv(double df, double p)
{
 if(df <=0){return 0;}
 return 2*invgammap(1-p, 0.5*df);
}


double fisher2(struct PATHWAY pathway, double inflationfactor,int thread, int adjust)
{
	int i=0;
	double sum=0;
	double p=1;
	double stat=0;

	if (pathway.counts == 0)
	{
		return 0;
	}
	else
	{
		for(i=0;i<pathway.singletop[thread];i++)
		{
			p= pathway.listp[thread][i];

			if(adjust)
			  {
			   stat=chiinv(1,p);
			   stat=stat/inflationfactor;
			   p=pValueCalc(0.5, stat / 2);
              }
			sum+=log(p);
		}

		return -2*sum;
	}
}


double maxT(struct PATHWAY pathway, double inflationfactor, int thread, int adjust)
{
	int i=0;
	double sum=0;
	double p=1;
	double stat=0;
	if (pathway.counts == 0)
	{
		return 0;
	}
	else
	{
		for(i=0;i<pathway.singletop[thread];i++)
		{
			if(pathway.listp_in[thread][i])
			{
				p= pathway.listp[thread][i];
				if(adjust)
				  {
					stat=chiinv(1,p);
					stat=stat/inflationfactor;
					p=pValueCalc(0.5, stat / 2);
				  }
				sum+=log(p);
			}
		}
		return -2*sum;
	}
}


double maxTplus(struct PATHWAY pathway, double inflationfactor, int thread, double cutoff, int adjust)
{
	int i=0;
	double sum=0;
	double p=1;
	double stat=0;
	if (pathway.counts == 0)
	{
		return 0;
	}
	else
	{
		for(i=0;i<pathway.singletop[thread];i++)
		{

			if(pathway.listp_in[thread][i])
			{

				p= pathway.listp[thread][i];
				if(adjust)
				  {
					stat=chiinv(1,p);
					stat=stat/inflationfactor;
					p=pValueCalc(0.5, stat / 2);
				  }
				sum+=log(p);
			}
			else if(pathway.listp_in2[thread][i])
			{

				p= pathway.listp[thread][i];
				if(adjust)
				  {
					stat=chiinv(1,p);
					stat=stat/inflationfactor;
					p=pValueCalc(0.5, stat / 2);
				  }
				if(p<=cutoff) {sum+=log(p);}
			}
			else if(pathway.listp_in3[thread][i])
			{

				p= pathway.listp[thread][i];
				if(adjust)
				  {
					stat=chiinv(1,p);
					stat=stat/inflationfactor;
					p=pValueCalc(0.5, stat / 2);
				  }
				if(p<=cutoff) {sum+=log(p);}
			}
		}
		return -2*sum;
	}
}


// pre-test allelic interaction
struct TSTAT pretest1df(double casecounts5[3][3][3][3][3],double controlcounts5[3][3][3][3][3],
                          	int dim, int mWithSingletop,double helpstat)
{
	double mAlt[3][3][3][3][3][2];
	double mNeu[3][3][3][3][3][2];
	double s1 = 0;
	double s2 = 0;
	double statAlt = 0;
	double statNeu = 0;
	int a, b, c, d, e, k, l;
	double df = 1;
	int stopC = 2;
	int stopD = 2;
	int stopE = 2;

	double casecountsAllele[2][2][2];
	double controlcountsAllele[2][2][2];
	double countsAllele[2][2][2];

	struct TSTAT tstat;
	double weight=0.5;

	tstat.p = 0;
	tstat.pmod = 0;

	if (dim == 2)
	{
	 stopC = 1;
	 stopD = 1;
	 stopE = 1;
	}

	if (dim == 3)
	{
	  stopD = 1;
	  stopE = 1;
	}
	else if (dim == 4)
	{
		stopE = 1;
	}

	if(dim > 2){printf("not yet implemented\n");exit(1);}


	//update df

	if(dim==2)
	{
	   casecountsAllele[0][0][0]=
	   2*casecounts5[0][0][0][0][0]+casecounts5[0][1][0][0][0]+casecounts5[1][0][0][0][0]+weight*casecounts5[1][1][0][0][0];
	   casecountsAllele[0][1][0]=
	   2*casecounts5[0][2][0][0][0]+casecounts5[0][1][0][0][0]+casecounts5[1][2][0][0][0]+weight*casecounts5[1][1][0][0][0];
	   casecountsAllele[1][0][0]=
	   2*casecounts5[2][0][0][0][0]+casecounts5[1][0][0][0][0]+casecounts5[2][1][0][0][0]+weight*casecounts5[1][1][0][0][0];
	   casecountsAllele[1][1][0]=
	   2*casecounts5[2][2][0][0][0]+casecounts5[1][2][0][0][0]+casecounts5[2][1][0][0][0]+weight*casecounts5[1][1][0][0][0];

	   controlcountsAllele[0][0][0]=
	   2*controlcounts5[0][0][0][0][0]+controlcounts5[0][1][0][0][0]+controlcounts5[1][0][0][0][0]+weight*controlcounts5[1][1][0][0][0];
	   controlcountsAllele[0][1][0]=
	   2*controlcounts5[0][2][0][0][0]+controlcounts5[0][1][0][0][0]+controlcounts5[1][2][0][0][0]+weight*controlcounts5[1][1][0][0][0];
	   controlcountsAllele[1][0][0]=
	   2*controlcounts5[2][0][0][0][0]+controlcounts5[1][0][0][0][0]+controlcounts5[2][1][0][0][0]+weight*controlcounts5[1][1][0][0][0];
	   controlcountsAllele[1][1][0]=
	   2*controlcounts5[2][2][0][0][0]+controlcounts5[1][2][0][0][0]+controlcounts5[2][1][0][0][0]+weight*controlcounts5[1][1][0][0][0];


	   casecountsAllele[0][0][0]*=2;casecountsAllele[0][1][0]*=2; // motivated by reasonning on plink page
	   casecountsAllele[1][0][0]*=2;casecountsAllele[1][1][0]*=2;
	   controlcountsAllele[0][0][0]*=2;controlcountsAllele[0][1][0]*=2;
	   controlcountsAllele[1][0][0]*=2;controlcountsAllele[1][1][0]*=2;


	   countsAllele[0][0][0]=casecountsAllele[0][0][0]+controlcountsAllele[0][0][0];
	   countsAllele[0][1][0]=casecountsAllele[0][1][0]+controlcountsAllele[0][1][0];
	   countsAllele[1][0][0]=casecountsAllele[1][0][0]+controlcountsAllele[1][0][0];
	   countsAllele[1][1][0]=casecountsAllele[1][1][0]+controlcountsAllele[1][1][0];

	   df=4;
	   if (countsAllele[0][0][0]==0){df--;}
	   if (countsAllele[0][1][0]==0){df--;}
	   if (countsAllele[1][0][0]==0){df--;}
	   if (countsAllele[1][1][0]==0){df--;}
	   if(df>=3){df=1;}
	   else{df=0;}
	}

	if (df > 0)
	{
		// Startvalues
		for (a = 0; a < 2; a++)
		{

			for (b = 0; b < 2; b++)
			{

				for (c = 0; c < stopC; c++)
				{

					for (d = 0; d < stopD; d++)
					{
						for (e = 0; e < stopE; e++)
						{
							for (k = 0; k < 2; k++)
							{
								mAlt[a][b][c][d][e][k] = 1;
							}//k
						}
					}
				}
			}//b
		} //end a startvalues

		for (l = 1; l < 1000; l++)
		{
			statNeu = 0;
			//Compute mNeu
			for (a = 0; a < 2; a++)
			{
				for (b = 0; b < 2; b++)
				{
					for (c = 0; c < stopC; c++)
					{
						for (d = 0; d < stopD; d++)
						{
							for (e = 0; e < stopE; e++)
							{
								for (k = 0; k < 2; k++)
								{
									if (dim == 2)
									{
										if (l % 3 == 1)
										{
											s1 = casecountsAllele[a][b][c] + controlcountsAllele[a][b][c];
											s2 = mAlt[a][b][c][d][e][0] + mAlt[a][b][c][d][e][1];
										}
										else if (l % 3 == 2)
										{
											if (k == 0)
											{
												s1 = casecountsAllele[a][0][c] + casecountsAllele[a][1][c];
											}
											else
											{
												s1 = controlcountsAllele[a][0][c] + controlcountsAllele[a][1][c];
											}
											s2 = mAlt[a][0][c][d][e][k] + mAlt[a][1][c][d][e][k];

										}
										else if (l % 3 == 0)
										{
											if (k == 0)
											{
												s1 = casecountsAllele[0][b][c] + casecountsAllele[1][b][c];
											}
											else
											{
												s1 = controlcountsAllele[0][b][c] + controlcountsAllele[1][b][c];
											}
											s2 = mAlt[0][b][c][d][e][k] + mAlt[1][b][c][d][e][k];

										}
									}

									if (s2 > 0)
									{
										mNeu[a][b][c][d][e][k] = mAlt[a][b][c][d][e][k] * s1 / s2;
									}
									else
									{
										mNeu[a][b][c][d][e][k] = 0;
									}
								} //k
							}
						}
					}
				} //b
			} //a end update mNeu


			//update mAlt and likelihood

			for (a = 0; a < 2; a++)
			{
				for (b = 0; b < 2; b++)
				{
					for (c = 0; c < stopC; c++)
					{
						for (d = 0; d < stopD; d++)
						{
							for (e = 0; e < stopE; e++)
							{
								for (k = 0; k < 2; k++)
								{
									mAlt[a][b][c][d][e][k] = mNeu[a][b][c][d][e][k];

									if (k == 0)
									{
										if (mNeu[a][b][c][d][e][k] > 0 && casecountsAllele[a][b][c]> 0)
										{
											statNeu+= -2 * (casecountsAllele[a][b][c]* log(mNeu[a][b][c][d][e][k]) - casecountsAllele[a][b][c]* log(casecountsAllele[a][b][c]));
										}
									} //k==0
									else
									{
										if (mNeu[a][b][c][d][e][k] > 0 && controlcountsAllele[a][b][c] > 0)
										{
											statNeu += -2 * (controlcountsAllele[a][b][c]* log(mNeu[a][b][c][d][e][k])- controlcountsAllele[a][b][c]* log(controlcountsAllele[a][b][c]));

										}
									}
								}//k
							}
						}
					}
				}//b
			} //a end update mAlt

			if (fabs(statAlt - statNeu) < 0.0001)
			{
				break;
			}
			statAlt = statNeu;
		} //end l
	} //end if df > 0


	//printf("df %f statNeu %f\n",df,statNeu);

	if (statNeu > 0 && df > 0)
	{
		tstat.p = pValueCalc(df / 2, statNeu / 2);
	}
	else
	{
		tstat.p = 1;
	}

	tstat.pmod = tstat.p;

	tstat.df=df;

	return tstat;

} //end pretest


int numberCharacter(string file)
{
	int i = 0;
	char c;
	FILE *fptr = NULL;
	fptr = fopen(file.c_str(), "r");

	while(1)
	{
		c = fgetc(fptr);
		if (c == '\n')
		{
			fclose(fptr);break;
		}
		i++;

    }

	return i;
}


double adjustP(double p, double inflationfactor, int singleMarkerTest)
{
	double stat=0,newp=1;
	double smadjust=1;

	if(singleMarkerTest==2 || singleMarkerTest==4){smadjust=0.5;}

	if(p>0)
	{
		stat=chiinv(1,p);
		stat=stat/(inflationfactor*smadjust);
		newp=pValueCalc(0.5, stat / 2);
		return newp;
	}
	else {return 0;}
}

void initCounts(struct COUNTS* counts, int n, int regression)
{
	(*counts).AA_Co = 0;
	(*counts).AA_Co_male = 0;
	(*counts).AA_Co_female = 0;
	(*counts).AA_Ca = 0;
	(*counts).AA_Ca_male = 0;
	(*counts).AA_Ca_female = 0;
	(*counts).BB_Co = 0;
	(*counts).BB_Co_male = 0;
	(*counts).BB_Co_female = 0;
	(*counts).BB_Ca = 0;
	(*counts).BB_Ca_male = 0;
	(*counts).BB_Ca_female = 0;
	(*counts).AB_Co = 0;
	(*counts).AB_Co_male = 0;
	(*counts).AB_Co_female = 0;
	(*counts).AB_Ca = 0;
	(*counts).AB_Ca_male = 0;
	(*counts).AB_Ca_female = 0;
	(*counts).pSingle = 1;

	(*counts).OO_Co = 0;
	(*counts).OO_Co_male = 0;
	(*counts).OO_Co_female = 0;
	(*counts).OO_Ca = 0;
	(*counts).OO_Ca_male = 0;
	(*counts).OO_Ca_female = 0;


	if(n==0)
	{
//		(*counts).det = (struct DETAILS *) calloc(1, sizeof(struct DETAILS));
		(*counts).det = (struct DETAILS *) realloc((*counts).det, sizeof(struct DETAILS));
		(*counts).det->aCaN = 0;
		(*counts).det->aCoN = 0;
		(*counts).det->bCaN = 0;
		(*counts).det->bCoN = 0;
		(*counts).det->aCa = 0;
		(*counts).det->bCa = 0;
		(*counts).det->aCo = 0;
		(*counts).det->bCo = 0;
		(*counts).det->orA = 0;
		(*counts).det->orB = 0;
		(*counts).det->lclA = 0;
		(*counts).det->rclA = 0;
		(*counts).det->lclB = 0;
		(*counts).det->rclB = 0;
		(*counts).det->aaCa = 0;
		(*counts).det->abCa = 0;
		(*counts).det->bbCa = 0;
		(*counts).det->aaCo = 0;
		(*counts).det->abCo = 0;
		(*counts).det->bbCo = 0;
		(*counts).det->testHWE_Ca = 0;
		(*counts).det->testHWE_Co = 0;
		if(regression)
		  {
		    ((*counts).result1)=(struct STATplus *) malloc(1*sizeof(struct STATplus));
		  }
	}
};

//binäre Kodierung
void updateCounts(struct COUNTS* counts, int nwordsSNPs, uint64_t** BinSNPs, uint64_t** BinSNPsCCFlags, uint64_t** BinSNPsGenderFlags, bool gender)
{
	for (int p=0; p<nwordsSNPs; p++)
	{
		// Kontrollen
		(*counts).AA_Co += bitcount64(BinSNPs[p][1] & BinSNPsCCFlags[p][1]);
		(*counts).AB_Co += bitcount64(BinSNPs[p][2] & BinSNPsCCFlags[p][1]);
		(*counts).BB_Co += bitcount64(BinSNPs[p][3] & BinSNPsCCFlags[p][1]);
		(*counts).OO_Co += bitcount64(BinSNPs[p][0] & BinSNPsCCFlags[p][1]);

		// Fälle
		(*counts).AA_Ca += bitcount64(BinSNPs[p][1] & BinSNPsCCFlags[p][2]);
		(*counts).AB_Ca += bitcount64(BinSNPs[p][2] & BinSNPsCCFlags[p][2]);
		(*counts).BB_Ca += bitcount64(BinSNPs[p][3] & BinSNPsCCFlags[p][2]);
		(*counts).OO_Ca += bitcount64(BinSNPs[p][0] & BinSNPsCCFlags[p][2]);
	}
	if (gender)
	{
        for (int p=0; p<nwordsSNPs; p++)
        {
            // Kontrollen
			(*counts).AA_Co_male += bitcount64(BinSNPs[p][1] & BinSNPsCCFlags[p][1] & BinSNPsGenderFlags[p][1]);
			(*counts).AB_Co_male += bitcount64(BinSNPs[p][2] & BinSNPsCCFlags[p][1] & BinSNPsGenderFlags[p][1]);
			(*counts).BB_Co_male += bitcount64(BinSNPs[p][3] & BinSNPsCCFlags[p][1] & BinSNPsGenderFlags[p][1]);
			(*counts).OO_Co_male += bitcount64(BinSNPs[p][0] & BinSNPsCCFlags[p][1] & BinSNPsGenderFlags[p][1]);

			(*counts).AA_Co_female += bitcount64(BinSNPs[p][1] & BinSNPsCCFlags[p][1] & BinSNPsGenderFlags[p][2]);
			(*counts).AB_Co_female += bitcount64(BinSNPs[p][2] & BinSNPsCCFlags[p][1] & BinSNPsGenderFlags[p][2]);
			(*counts).BB_Co_female += bitcount64(BinSNPs[p][3] & BinSNPsCCFlags[p][1] & BinSNPsGenderFlags[p][2]);
			(*counts).OO_Co_female += bitcount64(BinSNPs[p][0] & BinSNPsCCFlags[p][1] & BinSNPsGenderFlags[p][2]);

            // Fälle
			(*counts).AA_Ca_male += bitcount64(BinSNPs[p][1] & BinSNPsCCFlags[p][2] & BinSNPsGenderFlags[p][1]);
			(*counts).AB_Ca_male += bitcount64(BinSNPs[p][2] & BinSNPsCCFlags[p][2] & BinSNPsGenderFlags[p][1]);
			(*counts).BB_Ca_male += bitcount64(BinSNPs[p][3] & BinSNPsCCFlags[p][2] & BinSNPsGenderFlags[p][1]);
			(*counts).OO_Ca_male += bitcount64(BinSNPs[p][0] & BinSNPsCCFlags[p][2] & BinSNPsGenderFlags[p][1]);

			(*counts).AA_Ca_female += bitcount64(BinSNPs[p][1] & BinSNPsCCFlags[p][2] & BinSNPsGenderFlags[p][2]);
			(*counts).AB_Ca_female += bitcount64(BinSNPs[p][2] & BinSNPsCCFlags[p][2] & BinSNPsGenderFlags[p][2]);
			(*counts).BB_Ca_female += bitcount64(BinSNPs[p][3] & BinSNPsCCFlags[p][2] & BinSNPsGenderFlags[p][2]);
			(*counts).OO_Ca_female += bitcount64(BinSNPs[p][0] & BinSNPsCCFlags[p][2] & BinSNPsGenderFlags[p][2]);
		}
	}
}


int filledCells(double counts5[3][3][3][3][3], int cut)
{
	 int i,j,filled=0;

	 for(i=0;i<3;i++)
	    {
		 for(j=0;j<3;j++)
			{
			 if(counts5[i][j][0][0][0]>=cut){filled++;}
			}
		}

	 return filled;
	}


void updateCountsMulti(double casecounts5[3][3][3][3][3], double casecounts5M[3][3][3][3][3],double casecounts5F[3][3][3][3][3],double controlcounts5[3][3][3][3][3], double controlcounts5M[3][3][3][3][3],double controlcounts5F[3][3][3][3][3],int aMod, int bMod, int nwordsSNPs, uint64_t*** BinSNPs, uint64_t** BinSNPsCCFlags, uint64_t** BinSNPsGenderFlags, bool X, int caseOnly, uint64_t* BinSNPsQTaffFlags,int test)
{

	int p;
	int aa=0;
	int bb=0;

	if(!caseOnly && !X)
	  {
	   for (p=0; p<nwordsSNPs; p++)
			{
				for(aa=0; aa<3; aa++)
				{
					for( bb=0; bb<3; bb++)
					{
						//Fälle
						 casecounts5[aa][bb][0][0][0] += bitcount64(BinSNPs[aMod][p][aa+1] & BinSNPs[bMod][p][bb+1] & BinSNPsCCFlags[p][2]);

						//Kontrollen
						 controlcounts5[aa][bb][0][0][0] += bitcount64(BinSNPs[aMod][p][aa+1] & BinSNPs[bMod][p][bb+1] & BinSNPsCCFlags[p][1]);
					}
				}
			}
	  }
	else if (!caseOnly && X)
	  {
	   for (p=0; p<nwordsSNPs; p++)
		{
			for(aa=0; aa<3; aa++)
			{
				for( bb=0; bb<3; bb++)
				{
					//Fälle

						casecounts5[aa][bb][0][0][0] += bitcount64(((BinSNPs[aMod][p][aa+1]) & (BinSNPsCCFlags[p][2])) & ((BinSNPs[bMod][p][bb+1]) ));
						casecounts5M[aa][bb][0][0][0] += bitcount64((((BinSNPs[aMod][p][aa+1]) & (BinSNPsCCFlags[p][2])) & ((BinSNPs[bMod][p][bb+1]) )) & BinSNPsGenderFlags[p][1]);
						casecounts5F[aa][bb][0][0][0] += bitcount64((((BinSNPs[aMod][p][aa+1]) & (BinSNPsCCFlags[p][2])) & ((BinSNPs[bMod][p][bb+1]) )) & BinSNPsGenderFlags[p][2]);


					//Kontrollen

					   controlcounts5[aa][bb][0][0][0] += bitcount64(((BinSNPs[aMod][p][aa+1]) & (BinSNPsCCFlags[p][1])) & ((BinSNPs[bMod][p][bb+1])  ));
					   controlcounts5M[aa][bb][0][0][0] += bitcount64((((BinSNPs[aMod][p][aa+1]) & (BinSNPsCCFlags[p][1])) & ((BinSNPs[bMod][p][bb+1])  )) & BinSNPsGenderFlags[p][1]);
					   controlcounts5F[aa][bb][0][0][0] += bitcount64((((BinSNPs[aMod][p][aa+1]) & (BinSNPsCCFlags[p][1])) & ((BinSNPs[bMod][p][bb+1]) )) & BinSNPsGenderFlags[p][2]);

				}
			}
		}
	  }
	else if (caseOnly && test!=17 && test !=18 && !X)
	  {
	   for (p=0; p<nwordsSNPs; p++)
		{
			for(aa=0; aa<3; aa++)
			{
				for( bb=0; bb<3; bb++)
				{
					//Fälle
					  casecounts5[aa][bb][0][0][0] += bitcount64(((BinSNPs[aMod][p][aa+1]) & (BinSNPsCCFlags[p][2])) & ((BinSNPs[bMod][p][bb+1]) ));
				}
			}
		}
	  }
    else if (caseOnly && test!=17 && test !=18 && X)
	  {
	   for (p=0; p<nwordsSNPs; p++)
		{
			for(aa=0; aa<3; aa++)
			{
				for( bb=0; bb<3; bb++)
				{
					//Fälle
					  casecounts5[aa][bb][0][0][0] += bitcount64(((BinSNPs[aMod][p][aa+1]) & (BinSNPsCCFlags[p][2])) & ((BinSNPs[bMod][p][bb+1]) ));
					  casecounts5M[aa][bb][0][0][0] += bitcount64((((BinSNPs[aMod][p][aa+1]) & (BinSNPsCCFlags[p][2])) & ((BinSNPs[bMod][p][bb+1]) )) & BinSNPsGenderFlags[p][1]);
					  casecounts5F[aa][bb][0][0][0] += bitcount64((((BinSNPs[aMod][p][aa+1]) & (BinSNPsCCFlags[p][2])) & ((BinSNPs[bMod][p][bb+1]) )) & BinSNPsGenderFlags[p][2]);

				}
			}
		}
	  }
    else if (caseOnly && (test==17 || test ==18) && !X)// caseOnly 17,18
      {
	   for (p=0; p<nwordsSNPs; p++)
		{
			for(aa=0; aa<3; aa++)
			{
				for( bb=0; bb<3; bb++)
				{
					//Fälle
					   casecounts5[aa][bb][0][0][0] += bitcount64(((BinSNPs[aMod][p][aa+1]) ) & ((BinSNPs[bMod][p][bb+1]) & (BinSNPsQTaffFlags[p])));
				}
			}
		}
      }
     else if (caseOnly && (test==17 || test ==18) && X)// caseOnly 17,18
      {
	   for (p=0; p<nwordsSNPs; p++)
		{
			for(aa=0; aa<3; aa++)
			{
				for( bb=0; bb<3; bb++)
				{
					//Fälle
					   casecounts5[aa][bb][0][0][0] += bitcount64(((BinSNPs[aMod][p][aa+1]) ) & ((BinSNPs[bMod][p][bb+1]) & (BinSNPsQTaffFlags[p])));
					   casecounts5M[aa][bb][0][0][0] += bitcount64((((BinSNPs[aMod][p][aa+1]) ) & ((BinSNPs[bMod][p][bb+1]) & (BinSNPsQTaffFlags[p]))) & BinSNPsGenderFlags[p][1]);
					   casecounts5F[aa][bb][0][0][0] += bitcount64((((BinSNPs[aMod][p][aa+1]) ) & ((BinSNPs[bMod][p][bb+1]) & (BinSNPsQTaffFlags[p]))) & BinSNPsGenderFlags[p][2]);

				}
			}
		}
      }

}

void updateCountsMulti3(double casecounts5[3][3][3][3][3], double casecounts5M[3][3][3][3][3], double casecounts5F[3][3][3][3][3], double controlcounts5[3][3][3][3][3], double controlcounts5M[3][3][3][3][3], double controlcounts5F[3][3][3][3][3], int aMod, int bMod, int cMod, int nwordsSNPs, uint64_t*** BinSNPs, uint64_t** BinSNPsCCFlags, uint64_t** BinSNPsGenderFlags, bool X, int caseOnly)
{
	int p;
	int aa = 0;
	int bb = 0;
	int cc = 0;

	for(p=0; p<nwordsSNPs;p++)
	{
		for(aa=0; aa<3; aa++)
		{
			for(bb=0; bb<3; bb++)
			{
				for(cc=0; cc<3; cc++)
				{
				//Fälle
					casecounts5[aa][bb][cc][0][0] += bitcount64(((BinSNPs[aMod][p][aa+1]) & (BinSNPsCCFlags[p][2])) & ((BinSNPs[bMod][p][bb+1]) & (BinSNPsCCFlags[p][2])) & ((BinSNPs[cMod][p][cc+1]) & (BinSNPsCCFlags[p][2])));

					if(X)
					{
						casecounts5M[aa][bb][cc][0][0] += bitcount64((((BinSNPs[aMod][p][aa+1]) & (BinSNPsCCFlags[p][2])) & ((BinSNPs[bMod][p][bb+1]) & (BinSNPsCCFlags[p][2])) & ((BinSNPs[cMod][p][cc+1]) & (BinSNPsCCFlags[p][2]))) & BinSNPsGenderFlags[p][1]);
						casecounts5F[aa][bb][cc][0][0] += bitcount64((((BinSNPs[aMod][p][aa+1]) & (BinSNPsCCFlags[p][2])) & ((BinSNPs[bMod][p][bb+1]) & (BinSNPsCCFlags[p][2])) & ((BinSNPs[cMod][p][cc+1]) & (BinSNPsCCFlags[p][2]))) & BinSNPsGenderFlags[p][2]);
					}

					if(!caseOnly)
					{
						//Kontrollen
						controlcounts5[aa][bb][cc][0][0] += bitcount64(((BinSNPs[aMod][p][aa+1]) & (BinSNPsCCFlags[p][1])) & ((BinSNPs[bMod][p][bb+1]) & (BinSNPsCCFlags[p][1])) & ((BinSNPs[cMod][p][cc+1]) & (BinSNPsCCFlags[p][1])));

						if(X)
						{
							controlcounts5M[aa][bb][cc][0][0] += bitcount64((((BinSNPs[aMod][p][aa+1]) & (BinSNPsCCFlags[p][1])) & ((BinSNPs[bMod][p][bb+1]) & (BinSNPsCCFlags[p][1])) & ((BinSNPs[cMod][p][cc+1]) & (BinSNPsCCFlags[p][1]))) & BinSNPsGenderFlags[p][1]);
							controlcounts5F[aa][bb][cc][0][0] += bitcount64((((BinSNPs[aMod][p][aa+1]) & (BinSNPsCCFlags[p][1])) & ((BinSNPs[bMod][p][bb+1]) & (BinSNPsCCFlags[p][1])) & ((BinSNPs[cMod][p][cc+1]) & (BinSNPsCCFlags[p][1]))) & BinSNPsGenderFlags[p][2]);
						}
					}
				}
			}
		}
	}
}


void wait(int ms) {
#if defined __linux__ || defined __unix__
    struct timespec delay;
    delay.tv_sec = 0;
    delay.tv_nsec = ms*1e6;
    nanosleep(&delay, NULL);
#else
    clock_t endwait = clock() + ms * CLOCKS_PER_SEC/1000;
    while (clock() < endwait) {}
#endif
}


struct TSTAT caseOnly4DF(double casecounts5[3][3][3][3][3])
{
	struct TSTAT tstat;
	tstat.p = 1;
	tstat.pmod = 1;
	tstat.df =0;


	double CasesExpected[4][4];
	double CasesObserved[4][4];

	double summaLine[]={0,0,0};
	double summaColumn[]={0,0,0};
	double summaTable=0;

	double summaLikelihood=0;

	double df=0;
	double dfSNP1=0;
	double dfSNP2=0;


//	int merkab = 0;


	for(int a=0; a<3; a++)
	{
		for(int b=0; b<3; b++)
		{
			CasesObserved[a][b]=casecounts5[a][b][0][0][0];
			//cout << "casecounts" << "\t" <<  casecounts5[a][b][0][0][0] << "CasesObserved" << "\t" << CasesObserved[a][b] << endl;
		}

	}


	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			summaLine[i]=summaLine[i]+CasesObserved[i][j];
			CasesObserved[i][3]=summaLine[i];
			summaColumn[i]=summaColumn[i]+CasesObserved[j][i];
			CasesObserved[3][i]=summaColumn[i];
			summaTable=summaTable+CasesObserved[i][j];
			CasesObserved[3][3]=summaTable;
		}
	}

// Cases Expected

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			if(CasesObserved[3][3]!=0)
				CasesExpected[i][j]=CasesObserved[i][3]*CasesObserved[3][j]/CasesObserved[3][3];
			else
				CasesExpected[i][j]=0;
		}
	}

	for(int i=0; i<=2; i++)
	{
		summaLine[i]=0;
		summaColumn[i]=0;
	}

	summaTable=0;

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			summaLine[i]=summaLine[i]+CasesExpected[i][j];
			CasesExpected[i][3]=summaLine[i];
			summaColumn[i]=summaColumn[i]+CasesExpected[j][i];
			CasesExpected[3][i]=summaColumn[i];
			summaTable=summaTable+CasesExpected[i][j];
			CasesExpected[3][3]=summaTable;
		}
	}

	for(int i=0; i<3; i++)
	{
		if(CasesObserved[i][3] != 0)
			dfSNP1=dfSNP1+1;
	}

	for(int i=0; i<3; i++)
	{
		if(CasesObserved[3][i] != 0)
			dfSNP2=dfSNP2+1;
	}

	df=(dfSNP1-1)*(dfSNP2-1);


// Likelihood-Ratio tets

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			if(CasesExpected[i][j]>0 && CasesObserved[i][j]>0)
				summaLikelihood=summaLikelihood+2*CasesObserved[i][j]*log(CasesObserved[i][j]/CasesExpected[i][j]);
		}
	}

	if(df>0)
	{
		tstat.p=pValueCalc(df/2, summaLikelihood/2);
	}
	else
		tstat.p=1;

	tstat.df=df;

return (tstat);
}

struct TSTAT caseOnly4DF_X(double casecounts5M[3][3][3][3][3], double casecounts5F[3][3][3][3][3])
{
	struct TSTAT tstat;
	tstat.p = 1;
	tstat.pmod = 1;
	tstat.df =0;


	double CasesObservedM[4][4];
	double CasesObservedF[4][4];

	double CasesExpectedM[4][4];
	double CasesExpectedF[4][4];


	double summaLine[]={0,0,0};
	double summaColumn[]={0,0,0};
	double summaTable=0;

	double summaLikelihoodM=0;
	double summaLikelihoodF=0;

	double dfM=0;
	double dfSNP1M=0;
	double dfSNP2M=0;

	double dfF=0;
	double dfSNP1F=0;
	double dfSNP2F=0;

	// int merkab = 0;



	//Cases Observed Males

	CasesObservedM[0][0]=casecounts5M[0][0][0][0][0];
	CasesObservedM[0][1]=0;
	CasesObservedM[0][2]=casecounts5M[0][2][0][0][0];
	CasesObservedM[1][0]=0;
	CasesObservedM[1][1]=0;
	CasesObservedM[1][2]=0;
	CasesObservedM[2][0]=casecounts5M[2][0][0][0][0];
	CasesObservedM[2][1]=0;
	CasesObservedM[2][2]=casecounts5M[2][2][0][0][0];

	//Cases Observed Females

	for(int a=0; a<3; a++)
	{
		for(int b=0; b<3; b++)
		{
			//CasesObservedM[a][b]=casecounts5M[a][b][0][0][0];
			CasesObservedF[a][b]=casecounts5F[a][b][0][0][0];
		}

	}


	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			summaLine[i]=summaLine[i]+CasesObservedM[i][j];
			CasesObservedM[i][3]=summaLine[i];
			summaColumn[i]=summaColumn[i]+CasesObservedM[j][i];
			CasesObservedM[3][i]=summaColumn[i];
			summaTable=summaTable+CasesObservedM[i][j];
			CasesObservedM[3][3]=summaTable;
		}
	}

	for(int i=0; i<=2; i++)
	{
		summaLine[i]=0;
		summaColumn[i]=0;
	}
	summaTable=0;



	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			summaLine[i]=summaLine[i]+CasesObservedF[i][j];
			CasesObservedF[i][3]=summaLine[i];
			summaColumn[i]=summaColumn[i]+CasesObservedF[j][i];
			CasesObservedF[3][i]=summaColumn[i];
			summaTable=summaTable+CasesObservedF[i][j];
			CasesObservedF[3][3]=summaTable;
		}
	}


// Cases Expected Males

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			if(CasesObservedM[3][3]!=0)
				CasesExpectedM[i][j]=CasesObservedM[i][3]*CasesObservedM[3][j]/CasesObservedM[3][3];
			else
				CasesExpectedM[i][j]=0;
		}
	}


	for(int i=0; i<=2; i++)
	{
		summaLine[i]=0;
		summaColumn[i]=0;
	}
	summaTable=0;

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			summaLine[i]=summaLine[i]+CasesExpectedM[i][j];
			CasesExpectedM[i][3]=summaLine[i];
			summaColumn[i]=summaColumn[i]+CasesExpectedM[j][i];
			CasesExpectedM[3][i]=summaColumn[i];
			summaTable=summaTable+CasesExpectedM[i][j];
			CasesExpectedM[3][3]=summaTable;
		}
	}


// Cases Expected Females

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			if(CasesObservedF[3][3]!=0)
				CasesExpectedF[i][j]=CasesObservedF[i][3]*CasesObservedF[3][j]/CasesObservedF[3][3];
			else
				CasesExpectedF[i][j]=0;
		}
	}

	for(int i=0; i<=2; i++)
	{
		summaLine[i]=0;
		summaColumn[i]=0;
	}

	summaTable=0;

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			summaLine[i]=summaLine[i]+CasesExpectedF[i][j];
			CasesExpectedF[i][3]=summaLine[i];
			summaColumn[i]=summaColumn[i]+CasesExpectedF[j][i];
			CasesExpectedF[3][i]=summaColumn[i];
			summaTable=summaTable+CasesExpectedF[i][j];
			CasesExpectedF[3][3]=summaTable;
		}
	}


	for(int i=0; i<3; i++)
	{
		if(CasesObservedM[i][3] != 0)
			dfSNP1M=dfSNP1M+1;
	}

	for(int i=0; i<3; i++)
	{
		if(CasesObservedM[3][i] != 0)
			dfSNP2M=dfSNP2M+1;
	}

	if(CasesObservedM[3][3]!=0)
		dfM=(dfSNP1M-1)*(dfSNP2M-1);

	for(int i=0; i<3; i++)
	{
		if(CasesObservedF[i][3] != 0)
			dfSNP1F=dfSNP1F+1;
	}

	for(int i=0; i<3; i++)
	{
		if(CasesExpectedF[3][i] != 0)
			dfSNP2F=dfSNP2F+1;
	}

	if(CasesObservedF[3][3]!=0)
		dfF=(dfSNP1F-1)*(dfSNP2F-1);


	// Likelihood-Ratio test

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			if(CasesExpectedM[i][j]>0 && CasesObservedM[i][j]>0)
				summaLikelihoodM=summaLikelihoodM+2*CasesObservedM[i][j]*log(CasesObservedM[i][j]/CasesExpectedM[i][j]);
		}
	}

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			if(CasesExpectedF[i][j]>0 && CasesObservedF[i][j]>0)
				summaLikelihoodF=summaLikelihoodF+2*CasesObservedF[i][j]*log(CasesObservedF[i][j]/CasesExpectedF[i][j]);
		}
	}



	if(dfM>0)
	{
		if(dfF>0)
		{
			tstat.p=pValueCalc((dfM+dfF)/2, (summaLikelihoodM+summaLikelihoodF)/2);
			tstat.df=dfM+dfF;
		}
		else
		{
			tstat.p=pValueCalc(dfM/2, summaLikelihoodM/2);
			tstat.df=dfM;
		}
	}
	else
	{
		if(dfF>0)
		{
			tstat.p=pValueCalc(dfF/2, summaLikelihoodF/2);
			tstat.df=dfF;
		}
		else
			tstat.p=1;
	}


return (tstat);
}

struct TSTAT caseOnly1DF(double casecounts5[3][3][3][3][3])
{
	struct TSTAT tstat;
	tstat.p = 1;
	tstat.pmod = 1;
	tstat.df =0;

	double CasesExpected[4][4];
	double CasesObserved[4][4];

	double CasesObservedAlleles[3][3];
	double CasesExpectedAlleles[3][3];

	double summaLine[]={0,0,0};
	double summaColumn[]={0,0,0};

	double summaTable=0;

	double summaLikelihood=0;

	double df;

	int merkab = 0;


	for(int a=0; a<3; a++)
	{
		for(int b=0; b<3; b++)
		{
			CasesObserved[a][b]=casecounts5[a][b][0][0][0];
			//cout << "casecounts" << "\t" <<  casecounts5[a][b][0][0][0] << "CasesObserved" << "\t" << CasesObserved[a][b] << endl;
		}

	}

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			summaLine[i]=summaLine[i]+CasesObserved[i][j];
			CasesObserved[i][3]=summaLine[i];
			summaColumn[i]=summaColumn[i]+CasesObserved[j][i];
			CasesObserved[3][i]=summaColumn[i];
			summaTable=summaTable+CasesObserved[i][j];
			CasesObserved[3][3]=summaTable;
		}
	}


	// Observed Alleles Table

	CasesObservedAlleles[0][0]=4*CasesObserved[0][0]+2*CasesObserved[0][1]+2*CasesObserved[1][0]+CasesObserved[1][1];
	CasesObservedAlleles[0][1]=4*CasesObserved[0][2]+2*CasesObserved[0][1]+2*CasesObserved[1][2]+CasesObserved[1][1];
	CasesObservedAlleles[1][0]=4*CasesObserved[2][0]+2*CasesObserved[1][0]+2*CasesObserved[2][1]+CasesObserved[1][1];
	CasesObservedAlleles[1][1]=4*CasesObserved[2][2]+2*CasesObserved[1][2]+2*CasesObserved[2][1]+CasesObserved[1][1];

	for(int i=0; i<=2; i++)
	{
		summaLine[i]=0;
		summaColumn[i]=0;
	}

	summaTable=0;

	for(int i=0; i<2; i++)
	{
		for(int j=0; j<2; j++)
		{
			summaLine[i]=summaLine[i]+CasesObservedAlleles[i][j];
			CasesObservedAlleles[i][2]=summaLine[i];
			summaColumn[i]=summaColumn[i]+CasesObservedAlleles[j][i];
			CasesObservedAlleles[2][i]=summaColumn[i];
			summaTable=summaTable+CasesObservedAlleles[i][j];
			CasesObservedAlleles[2][2]=summaTable;
		}
	}



	// Alleles Expected

	for(int i=0; i<2; i++)
	{
		for(int j=0; j<2; j++)
		{
			if(CasesObservedAlleles[2][2]!=0)
				CasesExpectedAlleles[i][j]=CasesObservedAlleles[i][2]*CasesObservedAlleles[2][j]/CasesObservedAlleles[2][2];
			else
				CasesExpectedAlleles[i][j]=0;
		}
	}

	for(int i=0; i<=2; i++)
	{
		summaLine[i]=0;
		summaColumn[i]=0;
	}
	summaTable=0;


	for(int i=0; i<2; i++)
	{
		for(int j=0; j<2; j++)
		{
			summaLine[i]=summaLine[i]+CasesExpectedAlleles[i][j];
			CasesExpectedAlleles[i][2]=summaLine[i];
			summaColumn[i]=summaColumn[i]+CasesExpectedAlleles[j][i];
			CasesExpectedAlleles[2][i]=summaColumn[i];
			summaTable=summaTable+CasesExpectedAlleles[i][j];
			CasesExpectedAlleles[2][2]=summaTable;
		}
	}


	df=1;


	if(CasesObservedAlleles[0][2]==0 || CasesObservedAlleles[1][2]==0 || CasesObservedAlleles[2][0]==0 || CasesObservedAlleles[2][1]==0)
		df=0;

// Likelihood-Ratio test

	for(int i=0; i<2; i++)
	{
		for(int j=0; j<2; j++)
		{
			if(CasesExpectedAlleles[i][j]>0 && CasesObservedAlleles[i][j]>0)
				summaLikelihood=summaLikelihood+2*CasesObservedAlleles[i][j]*log(CasesObservedAlleles[i][j]/CasesExpectedAlleles[i][j]);
		}
	}

	if(df==1)
		{
			tstat.p=pValueCalc(df/2, summaLikelihood/2);
		}
	else
		tstat.p=1;

	tstat.df=df;

return (tstat);
}

struct TSTAT caseOnly1DF_X(double casecounts5M[3][3][3][3][3], double casecounts5F[3][3][3][3][3])
{
	struct TSTAT tstat;
	tstat.p = 1;
	tstat.pmod = 1;
	tstat.df =0;


	double CasesObservedM[4][4];
	double CasesObservedF[4][4];

	double CasesObservedMAlleles[3][3];
	double CasesObservedFAlleles[3][3];

	double CasesExpectedM[4][4];
	double CasesExpectedF[4][4];

	double CasesExpectedMAlleles[3][3];
	double CasesExpectedFAlleles[3][3];


	double summaLine[]={0,0,0};
	double summaColumn[]={0,0,0};

	double summaTable=0;

	double summaLikelihoodM=0;
	double summaLikelihoodF=0;

	double dfM;
	double dfF;

	//int merkab = 0;


	//Cases Observed Males

	CasesObservedM[0][0]=casecounts5M[0][0][0][0][0];
	CasesObservedM[0][1]=0;
	CasesObservedM[0][2]=casecounts5M[0][2][0][0][0];
	CasesObservedM[1][0]=0;
	CasesObservedM[1][1]=0;
	CasesObservedM[1][2]=0;
	CasesObservedM[2][0]=casecounts5M[2][0][0][0][0];
	CasesObservedM[2][1]=0;
	CasesObservedM[2][2]=casecounts5M[2][2][0][0][0];

	//Cases Observed Females

	for(int a=0; a<3; a++)
	{
		for(int b=0; b<3; b++)
		{
			CasesObservedF[a][b]=casecounts5F[a][b][0][0][0];
		}
	}

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			summaLine[i]=summaLine[i]+CasesObservedM[i][j];
			CasesObservedM[i][3]=summaLine[i];
			summaColumn[i]=summaColumn[i]+CasesObservedM[j][i];
			CasesObservedM[3][i]=summaColumn[i];
			summaTable=summaTable+CasesObservedM[i][j];
			CasesObservedM[3][3]=summaTable;
		}
	}

	for(int i=0; i<=2; i++)
	{
		summaLine[i]=0;
		summaColumn[i]=0;
	}
	summaTable=0;

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			summaLine[i]=summaLine[i]+CasesObservedF[i][j];
			CasesObservedF[i][3]=summaLine[i];
			summaColumn[i]=summaColumn[i]+CasesObservedF[j][i];
			CasesObservedF[3][i]=summaColumn[i];
			summaTable=summaTable+CasesObservedF[i][j];
			CasesObservedF[3][3]=summaTable;
		}
	}



// Observed Alleles Table (males)

	CasesObservedMAlleles[0][0]=CasesObservedM[0][0];
	CasesObservedMAlleles[0][1]=CasesObservedM[0][2];
	CasesObservedMAlleles[1][0]=CasesObservedM[2][0];
	CasesObservedMAlleles[1][1]=CasesObservedM[2][2];

	for(int i=0; i<=2; i++)
	{
		summaLine[i]=0;
		summaColumn[i]=0;
	}
	summaTable=0;

	for(int i=0; i<2; i++)
	{
		for(int j=0; j<2; j++)
		{
			summaLine[i]=summaLine[i]+CasesObservedMAlleles[i][j];
			CasesObservedMAlleles[i][2]=summaLine[i];
			summaColumn[i]=summaColumn[i]+CasesObservedMAlleles[j][i];
			CasesObservedMAlleles[2][i]=summaColumn[i];
			summaTable=summaTable+CasesObservedMAlleles[i][j];
			CasesObservedMAlleles[2][2]=summaTable;
		}
	}


// Observed Alleles Table (females)

	CasesObservedFAlleles[0][0]=4*CasesObservedF[0][0]+2*CasesObservedF[0][1]+2*CasesObservedF[1][0]+CasesObservedF[1][1];
	CasesObservedFAlleles[0][1]=4*CasesObservedF[0][2]+2*CasesObservedF[0][1]+2*CasesObservedF[1][2]+CasesObservedF[1][1];
	CasesObservedFAlleles[1][0]=4*CasesObservedF[2][0]+2*CasesObservedF[1][0]+2*CasesObservedF[2][1]+CasesObservedF[1][1];
	CasesObservedFAlleles[1][1]=4*CasesObservedF[2][2]+2*CasesObservedF[1][2]+2*CasesObservedF[2][1]+CasesObservedF[1][1];


	for(int i=0; i<=2; i++)
	{
		summaLine[i]=0;
		summaColumn[i]=0;
	}
	summaTable=0;

	for(int i=0; i<2; i++)
	{
		for(int j=0; j<2; j++)
		{
			summaLine[i]=summaLine[i]+CasesObservedFAlleles[i][j];
			CasesObservedFAlleles[i][2]=summaLine[i];
			summaColumn[i]=summaColumn[i]+CasesObservedFAlleles[j][i];
			CasesObservedFAlleles[2][i]=summaColumn[i];
			summaTable=summaTable+CasesObservedFAlleles[i][j];
			CasesObservedFAlleles[2][2]=summaTable;
		}
	}



// Expected Alleles (males)

	for(int i=0; i<2; i++)
	{
		for(int j=0; j<2; j++)
		{
			if(CasesObservedMAlleles[2][2]!=0)
				CasesExpectedMAlleles[i][j]=CasesObservedMAlleles[i][2]*CasesObservedMAlleles[2][j]/CasesObservedMAlleles[2][2];
			else
				CasesExpectedMAlleles[i][j]=0;
		}
	}

// Expected Alleles (females)

	for(int i=0; i<2; i++)
	{
		for(int j=0; j<2; j++)
		{
			if(CasesObservedFAlleles[2][2]!=0)
				CasesExpectedFAlleles[i][j]=CasesObservedFAlleles[i][2]*CasesObservedFAlleles[2][j]/CasesObservedFAlleles[2][2];
			else
			CasesExpectedFAlleles[i][j]=0;
		}
	}


	for(int i=0; i<=2; i++)
	{
		summaLine[i]=0;
		summaColumn[i]=0;
	}
	summaTable=0;

	for(int i=0; i<2; i++)
	{
		for(int j=0; j<2; j++)
		{
			summaLine[i]=summaLine[i]+CasesExpectedMAlleles[i][j];
			CasesExpectedMAlleles[i][2]=summaLine[i];
			summaColumn[i]=summaColumn[i]+CasesExpectedMAlleles[j][i];
			CasesExpectedMAlleles[2][i]=summaColumn[i];
			summaTable=summaTable+CasesExpectedMAlleles[i][j];
			CasesExpectedMAlleles[2][2]=summaTable;
		}
	}



	for(int i=0; i<=2; i++)
	{
		summaLine[i]=0;
		summaColumn[i]=0;
	}
	summaTable=0;

	for(int i=0; i<2; i++)
	{
		for(int j=0; j<2; j++)
		{
			summaLine[i]=summaLine[i]+CasesExpectedFAlleles[i][j];
			CasesExpectedFAlleles[i][2]=summaLine[i];
			summaColumn[i]=summaColumn[i]+CasesExpectedFAlleles[j][i];
			CasesExpectedFAlleles[2][i]=summaColumn[i];
			summaTable=summaTable+CasesExpectedFAlleles[i][j];
			CasesExpectedFAlleles[2][2]=summaTable;
		}
	}

// dfM

	dfM=1;

	if(CasesObservedMAlleles[0][2]==0 || CasesObservedMAlleles[1][2]==0 || CasesObservedMAlleles[2][0]==0 || CasesObservedMAlleles[2][1]==0)
		dfM=0;

// dfF

	dfF=1;

	if(CasesObservedFAlleles[0][2]==0 || CasesObservedFAlleles[1][2]==0 || CasesObservedFAlleles[2][0]==0 || CasesObservedFAlleles[2][1]==0)
		dfF=0;


// Likelihood-Ratio test

	for(int i=0; i<2; i++)
	{
		for(int j=0; j<2; j++)
		{
			if(CasesExpectedMAlleles[i][j]>0 && CasesObservedMAlleles[i][j]>0)
				summaLikelihoodM=summaLikelihoodM+2*CasesObservedMAlleles[i][j]*log(CasesObservedMAlleles[i][j]/CasesExpectedMAlleles[i][j]);
		}
	}



	for(int i=0; i<2; i++)
	{
		for(int j=0; j<2; j++)
		{
			if(CasesExpectedFAlleles[i][j]>0 && CasesObservedFAlleles[i][j]>0)
				summaLikelihoodF=summaLikelihoodF+2*CasesObservedFAlleles[i][j]*log(CasesObservedFAlleles[i][j]/CasesExpectedFAlleles[i][j]);
		}
	}

	tstat.p=1;

	if(dfM==1)
	{
		if(dfF==1)
		{
			tstat.p=pValueCalc((dfM+dfF)/2, (summaLikelihoodM+summaLikelihoodF)/2);
			tstat.df=dfM+dfF;
		}
		else
		{
			tstat.p=pValueCalc(dfM/2, summaLikelihoodM/2);
			tstat.df=dfM;
		}
	}
	else
	{
		if(dfF==1)
		{
			tstat.p=pValueCalc(dfF/2, summaLikelihoodF/2);
			tstat.df=dfF;
		}
	}

	return (tstat);
}

struct TSTAT CaseOnly3SNPs_X(double casecounts5M[3][3][3][3][3], double casecounts5F[3][3][3][3][3])
{
	struct TSTAT tstat;
	tstat.p = 1;
	tstat.pmod = 1;
	tstat.df =0;

	double CasesExpectedM[3][3][3];
	double CasesExpectedF[3][3][3];

	double CasesObservedM[3][3][3];
	double CasesObservedF[3][3][3];

	double summaAllLine[3]={0,0,0};
	double summaAllColumn[3]={0,0,0};
	double summaAllApplic[3]={0,0,0};
	double summaTable=0;

	double dfM=0;
	double dfF=0;
	double dfGes=0;
	double dfSNP1=0;
	double dfSNP2=0;
	double dfSNP3=0;

	double summaLikelihoodM=0;
	double summaLikelihoodF=0;


	//Cases Observed Males

	CasesObservedM[0][0][0]=casecounts5M[0][0][0][0][0];
	CasesObservedM[0][1][0]=0;
	CasesObservedM[0][2][0]=casecounts5M[0][2][0][0][0];
	CasesObservedM[1][0][0]=0;
	CasesObservedM[1][1][0]=0;
	CasesObservedM[1][2][0]=0;
	CasesObservedM[2][0][0]=casecounts5M[2][0][0][0][0];
	CasesObservedM[2][1][0]=0;
	CasesObservedM[2][2][0]=casecounts5M[2][2][0][0][0];

	CasesObservedM[0][0][1]=0;
	CasesObservedM[0][1][1]=0;
	CasesObservedM[0][2][1]=0;
	CasesObservedM[1][0][1]=0;
	CasesObservedM[1][1][1]=0;
	CasesObservedM[1][2][1]=0;
	CasesObservedM[2][0][1]=0;
	CasesObservedM[2][1][1]=0;
	CasesObservedM[2][2][1]=0;

	CasesObservedM[0][0][2]=casecounts5M[0][0][2][0][0];
	CasesObservedM[0][1][2]=0;
	CasesObservedM[0][2][2]=casecounts5M[0][2][2][0][0];
	CasesObservedM[1][0][2]=0;
	CasesObservedM[1][1][2]=0;
	CasesObservedM[1][2][2]=0;
	CasesObservedM[2][0][2]=casecounts5M[2][0][2][0][0];
	CasesObservedM[2][1][2]=0;
	CasesObservedM[2][2][2]=casecounts5M[2][2][2][0][0];

	//Cases Observed Females

	for(int a=0; a<3; a++)
	{
		for(int b=0; b<3; b++)
		{
			for(int c=0; c<3; c++)
			{
				CasesObservedF[a][b][c]=casecounts5F[a][b][c][0][0];
			}
		}
	}



	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				summaAllLine[i]=summaAllLine[i]+CasesObservedM[j][i][k];
				summaAllColumn[i]=summaAllColumn[i]+CasesObservedM[j][k][i];
				summaAllApplic[i]=summaAllApplic[i]+CasesObservedM[i][j][k];
				summaTable=summaTable+CasesObservedM[i][j][k];
			}
		}
	}

	/*cout << "Cases observed table (males)" << "\n";

	for(int i=0; i<3; i++)
		{
			for(int j=0; j<3; j++)
				{
					for(int k=0; k<3; k++)
						{
							cout << "CasesObservedM" << "\t" << i << "\t" << j << "\t" << k << "\t" << CasesObservedM[i][j][k] << "\t";
							cout << endl;
						}
				}
			cout << endl;
		}
	cout << endl;*/


	// Cases Expected Males

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				if(summaTable!=0)
					CasesExpectedM[i][j][k]= summaAllApplic[i]*summaAllLine[j]*summaAllColumn[k]/(summaTable*summaTable);
				else
					CasesExpectedM[i][j][k]=0;
			}
		}
	}

	for(int i=0; i<=2; i++)
	{
		summaAllLine[i]=0;
		summaAllColumn[i]=0;
		summaAllApplic[i]=0;
	}

	summaTable=0;


	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				summaAllLine[i]=summaAllLine[i]+CasesExpectedM[j][i][k];
				summaAllColumn[i]=summaAllColumn[i]+CasesExpectedM[j][k][i];
				summaAllApplic[i]=summaAllApplic[i]+CasesExpectedM[i][j][k];
				summaTable=summaTable+CasesExpectedM[i][j][k];
			}
		}
	}

	/*cout << "Cases expected table males" << "\n";

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
			{
				for(int k=0; k<3; k++)
					{
						cout << "CasesExpectedM" << "\t" << i << "\t" << j << "\t" << k << "\t" << CasesExpectedM[i][j][k] << "\t";
						cout << endl;
					}
			}
		cout << endl;
	}
	cout << endl;


	cout << "Summa Table Males = " << summaTable << "\n";*/


	// determine degrees of freedom males
	for(int i=0; i<3;i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				if(CasesObservedM[i][j][k] != 0)
					dfGes=dfGes+1;
			}
		}
	}

	dfGes=dfGes-1;
	//cout << " dfGesM " << dfGes << "\n";

	for(int i=0; i<3; i++)
	{
		if(summaAllLine[i] != 0)
			dfSNP1=dfSNP1+1;

	}
	dfSNP1=dfSNP1-1;
	//cout << " dfSNP1M " << dfSNP1 << "\n";


	for(int i=0; i<3; i++)
	{
		if(summaAllColumn[i] != 0)
			dfSNP2=dfSNP2+1;
	}
	dfSNP2=dfSNP2-1;
	//cout << " dfSNP2M " << dfSNP2 << "\n";

	for(int i=0; i<3;i++)
	{
		if(summaAllApplic[i] != 0)
			dfSNP3=dfSNP3+1;
	}

	dfSNP3=dfSNP3-1;
	//cout << " dfSNP3M " << dfSNP3 << "\n";

	dfM=dfGes-dfSNP1-dfSNP2-dfSNP3;
	//cout << " dfM = " << dfM << "\n";


	// Likelihood-Ratio tets males

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				if(CasesExpectedM[i][j][k]>0 && CasesObservedM[i][j][k]>0)
					summaLikelihoodM=summaLikelihoodM+2*CasesObservedM[i][j][k]*log(CasesObservedM[i][j][k]/CasesExpectedM[i][j][k]);
			}
		}
	}

	//cout << "Stat Males = " << summaLikelihoodM << "\n";

	///Females

	// Cases Observed Females

	for(int i=0; i<=2; i++)
	{
		summaAllLine[i]=0;
		summaAllColumn[i]=0;
		summaAllApplic[i]=0;
	}

	summaTable=0;


	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				summaAllLine[i]=summaAllLine[i]+CasesObservedF[j][i][k];
				summaAllColumn[i]=summaAllColumn[i]+CasesObservedF[j][k][i];
				summaAllApplic[i]=summaAllApplic[i]+CasesObservedF[i][j][k];
				summaTable=summaTable+CasesObservedF[i][j][k];
			}
		}
	}

	/*cout << "Cases observed table (females)" << "\n";

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				cout << "CasesObservedF" << "\t" << i << "\t" << j << "\t" << k << "\t" << CasesObservedF[i][j][k] << "\t";
				cout << endl;
			}
		}
		cout << endl;
	}
	cout << endl;*/


	// Cases Expected Females

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				if(summaTable!=0)
					CasesExpectedF[i][j][k]= summaAllApplic[i]*summaAllLine[j]*summaAllColumn[k]/(summaTable*summaTable);
				else
					CasesExpectedF[i][j][k]=0;
			}
		}
	}

	for(int i=0; i<=2; i++)
	{
		summaAllLine[i]=0;
		summaAllColumn[i]=0;
		summaAllApplic[i]=0;
	}

	summaTable=0;


	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				summaAllLine[i]=summaAllLine[i]+CasesExpectedF[j][i][k];
				summaAllColumn[i]=summaAllColumn[i]+CasesExpectedF[j][k][i];
				summaAllApplic[i]=summaAllApplic[i]+CasesExpectedF[i][j][k];
				summaTable=summaTable+CasesExpectedF[i][j][k];
			}
		}
	}

	/*cout << "Cases expected table females" << "\n";

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				cout << "CasesExpectedF" << "\t" << i << "\t" << j << "\t" << k << "\t" << CasesExpectedF[i][j][k] << "\t";
				cout << endl;
			}
		}
		cout << endl;
	}
	cout << endl;


	cout << "Summa Table Females = " << summaTable << "\n";		*/

	// determine degrees of freedom females

	dfGes=0;
	dfSNP1=0;
	dfSNP2=0;
	dfSNP3=0;


	for(int i=0; i<3;i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				if(CasesObservedF[i][j][k] != 0)
					dfGes=dfGes+1;
			}
		}
	}

	dfGes=dfGes-1;
	//cout << " dfGesF " << dfGes << "\n";

	for(int i=0; i<3; i++)
	{
		if(summaAllLine[i] != 0)
			dfSNP1=dfSNP1+1;

	}
	dfSNP1=dfSNP1-1;
	//cout << " dfSNP1F " << dfSNP1 << "\n";


	for(int i=0; i<3; i++)

	{
		if(summaAllColumn[i] != 0)
			dfSNP2=dfSNP2+1;
	}
	dfSNP2=dfSNP2-1;
	//cout << " dfSNP2F " << dfSNP2 << "\n";

	for(int i=0; i<3;i++)
	{
		if(summaAllApplic[i] != 0)
			dfSNP3=dfSNP3+1;
	}

	dfSNP3=dfSNP3-1;
	//cout << " dfSNP3F " << dfSNP3 << "\n";

	dfF=dfGes-dfSNP1-dfSNP2-dfSNP3;
	//cout << " dfF = " << dfF << "\n";


	// Likelihood-Ratio tets females

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				if(CasesExpectedF[i][j][k]>0 && CasesObservedF[i][j][k]>0)
					summaLikelihoodF=summaLikelihoodF+2*CasesObservedF[i][j][k]*log(CasesObservedF[i][j][k]/CasesExpectedF[i][j][k]);
			}
		}
	}

	//cout << "Stat Females = " << summaLikelihoodF << "\n";

	/// Join

	if(dfM>0)
	{
		if(dfF>0)
		{
			tstat.p=pValueCalc((dfM+dfF)/2, (summaLikelihoodM+summaLikelihoodF)/2);
			tstat.df = dfM+dfF;
			//cout << "df=" << dfM+dfF << endl;
			//cout << "Stat = " << summaLikelihoodM+summaLikelihoodF << endl;
		}
		else
		{
			tstat.p=pValueCalc(dfM/2, summaLikelihoodM/2);
			tstat.df = dfM;
			//cout << "df=" << dfM<< endl;
			//cout << "Stat = " << summaLikelihoodM << endl;
		}
	}
	else if(dfF>0)
	{
		tstat.p=pValueCalc(dfF/2, summaLikelihoodF/2);
		tstat.df = dfF;
		//cout << "df=" << dfF << endl;
		//cout << "Stat = " << summaLikelihoodF << endl;
	}

	return (tstat);
}


struct TSTAT caseOnly3_20DF(double casecounts5[3][3][3][3][3])
{
	struct TSTAT tstat;
	tstat.p = 1;
	tstat.pmod = 1;
	tstat.df =20;

	double CasesExpected[4][4][4];
	double CasesObserved[4][4][4];

	double summaLine[3][3]={{0,0,0},{0,0,0}};
	double summaColumn[3][3]={{0,0,0},{0,0,0}};
	double summa3rdDim[3]={0,0,0};
	double summaTable=0;
	double summaAllLine[3]={0,0,0};
	double summaAllColumn[3]={0,0,0};

	double summaLikelihood=0;

	double dfGes=0;
	double dfSNP1=0;
	double dfSNP2=0;
	double dfSNP3=0;
	double dfPair12=0;
	double dfPair13=0;
	double dfPair23=0;
	double dfTriple=0;
	double dfInteraction=0;


	for(int a=0; a<3; a++)
	{
		for(int b=0; b<3; b++)
		{
			for(int c=0; c<3; c++)
			{
				CasesObserved[a][b][c]=casecounts5[a][b][c][0][0];
			}
		}

	}

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				summaLine[i][j]=summaLine[i][j]+CasesObserved[i][j][k];
				summaColumn[i][j]=summaColumn[i][j]+CasesObserved[i][k][j];
				summa3rdDim[j]=summa3rdDim[j]+CasesObserved[j][k][i];
				summaTable=summaTable+CasesObserved[i][j][k];
			}

	//	CasesObserved[i][j][3]=summaLine[i][j];
	//	CasesObserved[i][3][j]=summaColumn[i][j];
	//	CasesObserved[j][3][3]=summa3rdDim[j];
		}
	}

	// CasesObserved[3][3][3]=summaTable;

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			summaAllLine[i]= summaAllLine[i]+summaLine[j][i];
			summaAllColumn[i]=summaAllColumn[i]+summaColumn[j][i];

		}
	}


	/*cout << "Cases observed table" << "\n";

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				cout << "CasesObserved" << "\t" << i << "\t" << j << "\t" << k << "\t" << CasesObserved[i][j][k] << "\t";
				cout << endl;
			}
			cout << endl;
		}
		cout << endl;
	*/


	// determine degrees of freedom
		for(int i=0; i<3;i++)
		{
			for(int j=0; j<3; j++)
			{
				for(int k=0; k<3; k++)
				{
					if(CasesObserved[i][j][k] != 0)
						dfGes=dfGes+1;
				}
			}
		}

	dfGes=dfGes-1;
	//cout << " dfGes " << dfGes << "\n";

	for(int i=0; i<3; i++)
	{
		if(summaAllLine[i] != 0)
			dfSNP1=dfSNP1+1;

	}
	dfSNP1=dfSNP1-1;
	//cout << " dfSNP1 " << dfSNP1 << "\n";


	for(int i=0; i<3; i++)
	{
		if(summaAllColumn[i] != 0)
			dfSNP2=dfSNP2+1;
	}
	dfSNP2=dfSNP2-1;
	//cout << " dfSNP2 " << dfSNP2 << "\n";

	for(int i=0; i<3;i++)
	{
		if(summa3rdDim[i] != 0)
			dfSNP3=dfSNP3+1;
	}

	dfSNP3=dfSNP3-1;
	//cout << " dfSNP3 " << dfSNP3 << "\n";

	/*
	dfPair12=dfSNP1*dfSNP2;
	cout << " dfPair12 " << dfPair12 << "\n";

	dfPair13=dfSNP1*dfSNP3;
	cout << " dfPair13 " << dfPair13 << "\n";

	dfPair23=dfSNP2*dfSNP3;
	cout << " dfPair23 " << dfPair23 << "\n";

	dfTriple=dfSNP1*dfSNP2*dfSNP3;
	cout << " dfTriple " << dfTriple << "\n";

	dfInteraction= dfTriple+dfPair12+dfPair13+dfPair23;
	cout << " dfInteraction " << dfInteraction << "\n";
	*/

	tstat.df=dfGes-dfSNP1-dfSNP2-dfSNP3;
	//cout << " df " << tstat.df << "\n";

// Cases Expected

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				//if(CasesObserved[i][3][3]!=0)
				if(summaTable!=0)
					//CasesExpected[i][j][k]=CasesObserved[i][j][3]*CasesObserved[i][3][j]*CasesObserved[j][3][3]/(CasesObserved[3][3][3]*CasesObserved[3][3][3]);
					CasesExpected[i][j][k]= summaAllLine[j]*summaAllColumn[k]*summa3rdDim[i]/(summaTable*summaTable);
				else
					CasesExpected[i][j][k]=0;
			}
		}
	}

	for(int i=0; i<=2; i++)
	{
		for(int j=0; j<=2; j++)
		{
			summaLine[i][j]=0;
			summaColumn[i][j]=0;

		}

		summaAllLine[i]=0;
		summaAllColumn[i]=0;
		summa3rdDim[i]=0;
	}

	summaTable=0;


	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				summaLine[i][j]=summaLine[i][j]+CasesExpected[i][j][k];
				summaColumn[i][j]=summaColumn[i][j]+CasesExpected[i][k][j];
				summa3rdDim[j]=summa3rdDim[j]+CasesExpected[j][k][i];
				summaTable=summaTable+CasesExpected[i][j][k];
			}

			//CasesExpected[i][j][3]=summaLine[i][j];
			//CasesExpected[i][3][j]=summaColumn[i][j];
			//CasesExpected[j][3][3]=summa3rdDim[j];
		}

	}

	//CasesExpected[3][3][3]=summaTable;

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			summaAllLine[i]= summaAllLine[i]+summaLine[j][i];
			summaAllColumn[i]=summaAllColumn[i]+summaColumn[j][i];
		}
	}


	//tstat.df=20;

	/*cout << "Cases expected table" << "\n";

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				cout << "CasesExpected" << "\t" << i << "\t" << j << "\t" << k << "\t" << CasesExpected[i][j][k] << "\t";
				cout << endl;
			}
		}
		cout << endl;
	}
	cout << endl;*/

// Likelihood-Ratio tets

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				if(CasesExpected[i][j][k]>0 && CasesObserved[i][j][k]>0)
					summaLikelihood=summaLikelihood+2*CasesObserved[i][j][k]*log(CasesObserved[i][j][k]/CasesExpected[i][j][k]);
			}
		}
	}

	if(tstat.df>0)
	{
		tstat.p=pValueCalc(tstat.df/2, summaLikelihood/2);
	}
	return (tstat);
}


double ComputePlusSingle(bool i, bool j, bool k, double dfInter, int test, bool ifChrX, double pSingleFirst, double pSingleSecond, double pSingleThird, double pInter, int dfSingle)
{
	double pValue=1;

	if (test!=19)
	{
		pValue=pValueCalc(0.5*(dfSingle*(i+j)+dfInter), 0.5*(chiinv(dfInter, pInter)+i*chiinv(dfSingle, pSingleFirst)+j*chiinv(dfSingle, pSingleSecond)));
	}
	else if (test==19)
	{
		pValue=pValueCalc(0.5*(dfSingle*(i+j+k)+dfInter), 0.5*(chiinv(dfInter, pInter)+i*chiinv(dfSingle, pSingleFirst)+j*chiinv(dfSingle, pSingleSecond)+k*chiinv(dfSingle, pSingleThird)));
	}

return(pValue);
}


void calcMaf(struct MAP* map, struct COUNTS oneCount) {
    if (!strcmp((*map).chr, "23")) {
		(*map).maf        = (double)(2*oneCount.AA_Ca_female+oneCount.AB_Ca_female+2*oneCount.AA_Co_female+oneCount.AB_Co_female+oneCount.AA_Ca_male+oneCount.AA_Co_male)
		                  / (double)(2*oneCount.AA_Ca_female+oneCount.AB_Ca_female+2*oneCount.AA_Co_female+2*oneCount.AB_Co_female+2*oneCount.BB_Ca_female+oneCount.AB_Ca_female+2*oneCount.BB_Co_female+oneCount.AA_Ca_male+oneCount.AA_Co_male+oneCount.BB_Ca_male+oneCount.BB_Co_male);
		(*map).controlmaf = (double)(2*oneCount.AA_Co_female+oneCount.AB_Co_female+oneCount.AA_Co_male)/(2*oneCount.AA_Co_female+2*oneCount.AB_Co_female+2*oneCount.BB_Co_female+oneCount.AA_Co_male+oneCount.BB_Co_male);
        (*map).casemaf    = (double)(2*oneCount.AA_Ca_female+oneCount.AB_Ca_female+oneCount.AA_Ca_male)/(2*oneCount.AA_Ca_female+2*oneCount.AB_Ca_female+2*oneCount.BB_Ca_female+oneCount.AA_Ca_male+oneCount.BB_Ca_male);
    } else
    if (!strcmp((*map).chr, "24")) {
        (*map).maf        = (double)(oneCount.AA_Ca_male+oneCount.AA_Co_male)/(oneCount.AA_Ca_male+oneCount.AA_Co_male+oneCount.BB_Ca_male+oneCount.BB_Co_male);
		(*map).controlmaf = (double)(oneCount.AA_Co_male)/(oneCount.AA_Co_male+oneCount.BB_Co_male);
		(*map).casemaf    = (double)(oneCount.AA_Ca_male)/(oneCount.AA_Ca_male+oneCount.BB_Ca_male);
    } else
    if (!strcmp((*map).chr, "26")) {
        (*map).maf        = (double)(oneCount.AA_Ca+oneCount.AA_Co)/(oneCount.AA_Ca+oneCount.BB_Ca+oneCount.AA_Co+oneCount.BB_Co);
		(*map).controlmaf = (double)(oneCount.AA_Co)/(oneCount.AA_Co+oneCount.BB_Co);
		(*map).casemaf    = (double)(oneCount.AA_Ca)/(oneCount.AA_Ca+oneCount.BB_Ca);
    } else {
        (*map).maf        = (double)(2*oneCount.AA_Ca+oneCount.AB_Ca+2*oneCount.AA_Co+oneCount.AB_Co)/(2*(oneCount.AA_Ca+oneCount.AB_Ca+oneCount.BB_Ca+oneCount.AA_Co+oneCount.AB_Co+oneCount.BB_Co));
        (*map).controlmaf = (double)(2*oneCount.AA_Co+oneCount.AB_Co)/(2*(oneCount.AA_Co+oneCount.AB_Co+oneCount.BB_Co));
		(*map).casemaf    = (double)(2*oneCount.AA_Ca+oneCount.AB_Ca)/(2*(oneCount.AA_Ca+oneCount.AB_Ca+oneCount.BB_Ca));
    }
}


void calcMaf(struct MAP* map, uint64_t** BinSNPs, uint64_t** BinSNPsCCFlags, uint64_t** BinSNPsGenderFlags, uint32_t nwordsSNPs) {
    uint32_t maf1num=0, maf1den=0, maf2num=0, maf2den=0, maf3num=0, maf3den=0;
    uint8_t a,b,c,d,e;
#if RARE
    uint32_t maf1numa=0, maf1dena=0, maf2numa=0, maf2dena=0, maf3numa=0, maf3dena=0;
    uint8_t aa,ba,ca,da,ea,fa,ga;
#endif
    if (!strcmp((*map).chr, "23")) {
        for (uint32_t p=0; p<nwordsSNPs; p++) {
            a = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][2] & ~BinSNPsCCFlags[p][0]);
            b = bitcount64(BinSNPs[p][2] & BinSNPsGenderFlags[p][2] & ~BinSNPsCCFlags[p][0]);
            c = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][2] & ~BinSNPsCCFlags[p][0]);
            d = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][1] & ~BinSNPsCCFlags[p][0]);
            e = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][1] & ~BinSNPsCCFlags[p][0]);
            maf1num += 2*a+b+d;
            maf1den += 2*(a+b+c)+d+e;
            a = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][2] &  BinSNPsCCFlags[p][1]);
            b = bitcount64(BinSNPs[p][2] & BinSNPsGenderFlags[p][2] &  BinSNPsCCFlags[p][1]);
            c = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][2] &  BinSNPsCCFlags[p][1]);
            d = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][1]);
            e = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][1]);
            maf2num += 2*a+b+d;
            maf2den += 2*(a+b+c)+d+e;
            a = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][2] &  BinSNPsCCFlags[p][2]);
            b = bitcount64(BinSNPs[p][2] & BinSNPsGenderFlags[p][2] &  BinSNPsCCFlags[p][2]);
            c = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][2] &  BinSNPsCCFlags[p][2]);
            d = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][2]);
            e = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][2]);
            maf3num += 2*a+b+d;
            maf3den += 2*(a+b+c)+d+e;
#if RARE
             aa = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][2] & ~BinSNPsCCFlags[p][0]);
             ba = bitcount64(BinSNPs[p][2] & BinSNPsGenderFlags[p][2] & ~BinSNPsCCFlags[p][0]);
             ca = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][2] & ~BinSNPsCCFlags[p][0]);
             da = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][1] & ~BinSNPsCCFlags[p][0]);
             ea = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][1] & ~BinSNPsCCFlags[p][0]);
             fa = bitcount64(BinSNPs[p][0] & BinSNPsGenderFlags[p][1] & ~BinSNPsCCFlags[p][0]);
             ga = bitcount64(BinSNPs[p][0] & BinSNPsGenderFlags[p][2] & ~BinSNPsCCFlags[p][0]);
             maf1numa += 2*aa+ba+da;
             maf1dena += 2*(aa+ba+ca+ga)+da+ea+fa;
             aa = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][2] &  BinSNPsCCFlags[p][1]);
             ba = bitcount64(BinSNPs[p][2] & BinSNPsGenderFlags[p][2] &  BinSNPsCCFlags[p][1]);
             ca = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][2] &  BinSNPsCCFlags[p][1]);
             da = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][1]);
             ea = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][1]);
             fa = bitcount64(BinSNPs[p][0] & BinSNPsGenderFlags[p][1] & ~BinSNPsCCFlags[p][1]);
             ga = bitcount64(BinSNPs[p][0] & BinSNPsGenderFlags[p][2] & ~BinSNPsCCFlags[p][1]);
             maf2numa += 2*aa+ba+da;
             maf2dena += 2*(aa+ba+ca+ga)+da+ea+fa;
             aa = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][2] &  BinSNPsCCFlags[p][2]);
             ba = bitcount64(BinSNPs[p][2] & BinSNPsGenderFlags[p][2] &  BinSNPsCCFlags[p][2]);
             ca = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][2] &  BinSNPsCCFlags[p][2]);
             da = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][2]);
             ea = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][2]);
             fa = bitcount64(BinSNPs[p][0] & BinSNPsGenderFlags[p][1] & ~BinSNPsCCFlags[p][2]);
             ga = bitcount64(BinSNPs[p][0] & BinSNPsGenderFlags[p][2] & ~BinSNPsCCFlags[p][2]);
             maf3numa += 2*aa+ba+da;
             maf3dena += 2*(aa+ba+ca+ga)+da+ea+fa;
#endif
        }
    }
      else
    if (!strcmp((*map).chr, "24")) {
        for (uint32_t p=0; p<nwordsSNPs; p++) {
            a = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][1] & ~BinSNPsCCFlags[p][0]);
            c = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][1] & ~BinSNPsCCFlags[p][0]);
            maf1num += a;
            maf1den += a+c;
            a = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][1]);
            c = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][1]);
            maf2num += a;
            maf2den += a+c;
            a = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][2]);
            c = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][2]);
            maf3num += a;
            maf3den += a+c;
#if RARE
            aa = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][1] & ~BinSNPsCCFlags[p][0]);
            ca = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][1] & ~BinSNPsCCFlags[p][0]);
            fa = bitcount64(BinSNPs[p][0] & BinSNPsGenderFlags[p][1] & ~BinSNPsCCFlags[p][0]);
            maf1numa += aa;
            maf1dena += aa+ca;
            aa = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][1]);
            ca = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][1]);
            fa = bitcount64(BinSNPs[p][0] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][1]);
            maf2numa += aa;
            maf2dena += aa+ca+fa;
            aa = bitcount64(BinSNPs[p][1] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][2]);
            ca = bitcount64(BinSNPs[p][3] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][2]);
            fa = bitcount64(BinSNPs[p][0] & BinSNPsGenderFlags[p][1] &  BinSNPsCCFlags[p][2]);
            maf3numa += aa;
            maf3dena += aa+ca+fa;
#endif
        }
    } else
    if (!strcmp((*map).chr, "26")) {
        for (uint32_t p=0; p<nwordsSNPs; p++) {
            a = bitcount64(BinSNPs[p][1] & ~BinSNPsCCFlags[p][0]);
            c = bitcount64(BinSNPs[p][3] & ~BinSNPsCCFlags[p][0]);
            maf1num += a;
            maf1den += a+c;
            a = bitcount64(BinSNPs[p][1] &  BinSNPsCCFlags[p][1]);
            c = bitcount64(BinSNPs[p][3] &  BinSNPsCCFlags[p][1]);
            maf2num += a;
            maf2den += a+c;
            a = bitcount64(BinSNPs[p][1] &  BinSNPsCCFlags[p][2]);
            c = bitcount64(BinSNPs[p][3] &  BinSNPsCCFlags[p][2]);
            maf3num += a;
            maf3den += a+c;
#if RARE
            aa = bitcount64(BinSNPs[p][1] & ~BinSNPsCCFlags[p][0]);
            ca = bitcount64(BinSNPs[p][3] & ~BinSNPsCCFlags[p][0]);
            fa = bitcount64(BinSNPs[p][0] & ~BinSNPsCCFlags[p][0]);
            maf1numa += aa;
            maf1dena += aa+ca+fa;
            aa = bitcount64(BinSNPs[p][1] &  BinSNPsCCFlags[p][1]);
            ca = bitcount64(BinSNPs[p][3] &  BinSNPsCCFlags[p][1]);
            fa = bitcount64(BinSNPs[p][0] &  BinSNPsCCFlags[p][1]);
            maf2numa += aa;
            maf2dena += aa+ca+fa;
            aa = bitcount64(BinSNPs[p][1] &  BinSNPsCCFlags[p][2]);
            ca = bitcount64(BinSNPs[p][3] &  BinSNPsCCFlags[p][2]);
            fa = bitcount64(BinSNPs[p][0] &  BinSNPsCCFlags[p][2]);
            maf3numa += aa;
            maf3dena += aa+ca+fa;
#endif
        }
    } else {
        for (uint32_t p=0; p<nwordsSNPs; p++) {
            a = bitcount64(BinSNPs[p][1] & ~BinSNPsCCFlags[p][0]);
            b = bitcount64(BinSNPs[p][2] & ~BinSNPsCCFlags[p][0]);
            c = bitcount64(BinSNPs[p][3] & ~BinSNPsCCFlags[p][0]);
            maf1num += 2*a+b;
            maf1den += 2*(a+b+c);
            a = bitcount64(BinSNPs[p][1] &  BinSNPsCCFlags[p][1]);
            b = bitcount64(BinSNPs[p][2] &  BinSNPsCCFlags[p][1]);
            c = bitcount64(BinSNPs[p][3] &  BinSNPsCCFlags[p][1]);
            maf2num += 2*a+b;
            maf2den += 2*(a+b+c);
            a = bitcount64(BinSNPs[p][1] &  BinSNPsCCFlags[p][2]);
            b = bitcount64(BinSNPs[p][2] &  BinSNPsCCFlags[p][2]);
            c = bitcount64(BinSNPs[p][3] &  BinSNPsCCFlags[p][2]);
            maf3num += 2*a+b;
            maf3den += 2*(a+b+c);
#if RARE
            aa = bitcount64(BinSNPs[p][1] & ~BinSNPsCCFlags[p][0]);
            ba = bitcount64(BinSNPs[p][2] & ~BinSNPsCCFlags[p][0]);
            ca = bitcount64(BinSNPs[p][3] & ~BinSNPsCCFlags[p][0]);
            fa = bitcount64(BinSNPs[p][0] & ~BinSNPsCCFlags[p][0]);
            maf1numa += 2*aa+ba;
            maf1dena += 2*(aa+ba+ca+fa);
            aa = bitcount64(BinSNPs[p][1] &  BinSNPsCCFlags[p][1]);
            ba = bitcount64(BinSNPs[p][2] &  BinSNPsCCFlags[p][1]);
            ca = bitcount64(BinSNPs[p][3] &  BinSNPsCCFlags[p][1]);
            fa = bitcount64(BinSNPs[p][0] &  BinSNPsCCFlags[p][1]);
            maf2numa += 2*aa+ba;
            maf2dena += 2*(aa+ba+ca+fa);
            aa = bitcount64(BinSNPs[p][1] &  BinSNPsCCFlags[p][2]);
            ba = bitcount64(BinSNPs[p][2] &  BinSNPsCCFlags[p][2]);
            ca = bitcount64(BinSNPs[p][3] &  BinSNPsCCFlags[p][2]);
            fa = bitcount64(BinSNPs[p][0] &  BinSNPsCCFlags[p][2]);
            maf3numa += 2*aa+ba;
            maf3dena += 2*(aa+ba+ca+fa);
#endif
        }
    }
     (*map).maf        =  (double)maf1num/maf1den;
     (*map).controlmaf =  (double)maf2num/maf2den;
     (*map).casemaf    =  (double)maf3num/maf3den;
#if RARE
     (*map).mafa        =  (double)maf1numa/maf1dena;
     (*map).controlmafa =  (double)maf2numa/maf2dena;
     (*map).casemafa    =  (double)maf3numa/maf3dena;

     if(mafadjust==1)
       {
	 (*map).mafr        =  rounddown((double)maf1numa/maf1dena,7);
	 (*map).controlmafr =  rounddown((double)maf2numa/maf2dena,7);
	 (*map).casemafr    =  rounddown((double)maf3numa/maf3dena,7);
       }
     else if(mafadjust==0)
       {
	 (*map).mafr        =  rounddown((double)maf1num/maf1den,7);
	 (*map).controlmafr =  rounddown((double)maf2num/maf2den,7);
	 (*map).casemafr    =  rounddown((double)maf3num/maf3den,7);
       }
#endif
}


void calcMR(struct MAP* map, uint64_t** BinSNPs, uint64_t** BinSNPsCCFlags, uint64_t** BinSNPsGenderFlags, uint32_t nwordsSNPs) {
  int missing=0,npplAff=0;

    if (!strcmp((*map).chr, "24")) {
        for (uint32_t p=0; p<nwordsSNPs; p++) {
            missing += bitcount64(BinSNPs[p][0] & ~BinSNPsCCFlags[p][0] & BinSNPsGenderFlags[p][1]);
            npplAff += bitcount64(~BinSNPsCCFlags[p][0] & BinSNPsGenderFlags[p][1]);
        }
    } else {
        for (uint32_t p=0; p<nwordsSNPs; p++) {
            missing += bitcount64(BinSNPs[p][0] & ~BinSNPsCCFlags[p][0]);
            npplAff += bitcount64(~BinSNPsCCFlags[p][0]);
        }
	}
    (*map).missing     = missing;
    (*map).missingrate = (double)(missing)/(npplAff);
}

// SNP-Covariates
void ReadSnpCovariates(int thread,int npplqc, uint64_t*** BinSNPs, int &maxIndexCov, int &numberOfSNPCov, int nlinestfam, int nsnpqc, struct PERSON *person, struct MAP *map, vector<string> SNPcovar, double missingpheno, int qt,int nlinestped, int liability,int dosage, float ***genoWeights)
{
	int l=0;
	int k = 0;
	int in = 0;
	struct PPLLOCATION guy;

	for (int i=0; i< nlinestfam; i++)
	{

		l = 0;


		person[i].allcovin = maxIndexCov;

		person[i].covin = (unsigned char *) realloc(person[i].covin,(maxIndexCov+numberOfSNPCov)*sizeof(unsigned char));
		if (!person[i].covin) die("memory allocation error in person[i].covin");

		person[i].cov= (double *) realloc(person[i].cov,(maxIndexCov+numberOfSNPCov)*sizeof(double));
		if (!person[i].cov) die("memory allocation error in person[i].cov");

		for (int j=0; j < numberOfSNPCov; j++)
		{
			for (int k = 0; k < nlinestped; k++)
			{
				if ((liability == 1 || liability ==2) && strcmp(map[k].chr,"23") && map[k].rs == SNPcovar[j])
				{
					if (getbit64(BinSNPs[k][i/64][1],i%64) == 1) // 11
					{
						person[i].cov[maxIndexCov+l] = 1.00;
						person[i].allcovin += 1;
						person[i].covin[maxIndexCov+l] = 1;
					}
					else if (getbit64(BinSNPs[k][i/64][2],i%64) == 1) // 12
					{
						person[i].cov[maxIndexCov+l] = 0.00;
						person[i].allcovin += 1;
						person[i].covin[maxIndexCov+l] = 1;
					}
					else if (getbit64(BinSNPs[k][i/64][3],i%64) == 1) // 22
					{
						person[i].cov[maxIndexCov+l] = -1.00;
						person[i].allcovin += 1;
						person[i].covin[maxIndexCov+l] = 1;
					}
					else
					{
						person[i].allcovin=0;
						person[i].qcin = 0;
						person[i].aff[thread] = 0;
						person[i].qtaff[thread] = 0;
						person[i].covin[maxIndexCov+l] = 0;
					}
					l++;

					break;
				}
				else if ((liability == 1 || liability ==2) && !strcmp(map[k].chr,"23") && map[k].rs == SNPcovar[j])
				{
					cout << "\n\nLiability model:\n" << map[k].rs << ": Risk SNPs on X chromosome can not be analyzed. Please remove these SNPs from the SNP_COVARIATES list\n";
					exit(1);
				}
				else if (!dosage && liability == 0 && map[k].rs == SNPcovar[j])
				{
					if (getbit64(BinSNPs[k][i/64][1],i%64) == 1) // 11
					{
						person[i].cov[maxIndexCov+l] = 1.00;
						person[i].allcovin += 1;
						person[i].covin[maxIndexCov+l] = 1;
					}
					else if (getbit64(BinSNPs[k][i/64][2],i%64) == 1) // 12
					{
						person[i].cov[maxIndexCov+l] = 0.00;
						person[i].allcovin += 1;
						person[i].covin[maxIndexCov+l] = 1;
					}
					else if (getbit64(BinSNPs[k][i/64][3],i%64) == 1) // 22
					{
						person[i].cov[maxIndexCov+l] = -1.00;
						person[i].allcovin += 1;
						person[i].covin[maxIndexCov+l] = 1;
					}
					else
					{
						person[i].allcovin=0;
						person[i].qcin = 0;
						person[i].aff[thread] = 0;
						person[i].qtaff[thread] = 0;
						person[i].covin[maxIndexCov+l] = 0;
					}
					l++;

					break;
				}
				else if (dosage==1 && liability == 0 && map[k].rs == SNPcovar[j])
				{
					if (1) // assumes dosage data has no missings
					{
						person[i].cov[maxIndexCov+l] = (2*genoWeights[k][i][0]+genoWeights[k][i][1])-1; //index? k=person, i=index
						person[i].allcovin += 1;
						person[i].covin[maxIndexCov+l] = 1;
					}
					else
					{
						person[i].allcovin=0;
						person[i].qcin = 0;
						person[i].aff[thread] = 0;
						person[i].qtaff[thread] = 0;
						person[i].covin[maxIndexCov+l] = 0;
					}
					l++;

					break;
				}
			}
		}
	}
	numberOfSNPCov = l;
	maxIndexCov = maxIndexCov + numberOfSNPCov;

	for (int i=0; i< nlinestfam; i++)
	{
		if (person[i].qcin == 1)
		{
			in++;
		}
	}

	cout << "Persons with all SNP covariates: " << in << "\n";
	//sstm << "After reading snp-covariate-file: ";
	//if (qt) sstm << "individuals: " << nlinestfam-nrest << ", notUsed: " << nrest << "\n";
	//else    sstm << "cases: " << ncases << ", controls: " << ncontrols << ", notUsed: " << nrest << "\n";
	//logg(sstm);

}


// Liability
void setLoad(int thread_nloop,int npplqc, int* PPLMap, int* SNPMap, uint64_t*** BinSNPs, struct PPLLOCATION* PplLocations, int maxIndexCov, int numberOfSNPCov, int nsnpqc, struct PERSON *person,struct MAP *map, struct COUNTS **counts,vector<string> SNPcovar, vector<string> SNPrisk, struct ALLELECODE *codesA)
{
	struct PPLLOCATION guy;
	int k = 0;
	int i = 0;
	string risk = "";
	string alleleA = "";
	string alleleB = "";

	for (int iMod=0; iMod< npplqc; iMod++)
	{
		i = PPLMap[iMod];
		guy = PplLocations[iMod];
		person[i].load = 0;
	}

	for (int j=0; j < numberOfSNPCov; j++)
	{
		int found=0;

		for (int kMod = 0; kMod < nsnpqc; kMod++)
		{
			k = SNPMap[kMod];

			alleleA = codesA[k].a1;
			alleleB = codesA[k].a2;

			if (SNPrisk[j] == "R")
			{
				if (counts[thread_nloop][k].det->orA > 1) // OR_A > 1?
				{
					risk = alleleA;
				}
				else
				{
					risk = alleleB;
				}
			}
			else
			{
				risk = SNPrisk[j]; // risk allele is read from selection file
			}

			if (map[k].rs == SNPcovar[j])
			{
				for (int iMod=0; iMod< npplqc; iMod++)
				{
					i = PPLMap[iMod];
					guy = PplLocations[iMod];

					if (getbit64(BinSNPs[kMod][guy.nr][1],guy.pos) == 1) // 11
					{
						if (risk == alleleA)
						{
							person[i].load += 2;
						}
					}
					else if (getbit64(BinSNPs[kMod][guy.nr][2],guy.pos)) // 12
					{
						person[i].load += 1;
					}
					else if (getbit64(BinSNPs[kMod][guy.nr][3],guy.pos) == 1) // 22
					{
						if (risk == alleleB)
						{
							person[i].load += 2;
						}
					}
				}
				if (j == 0)
				{
					cout << "\nrs_No" << "\t" << "minor allele (data)" << "\t" << "major allele (data)" << "\t" << "risk allele" << "\n";
				}

				if (SNPrisk[j] == "R")
				{
					cout << map[k].rs << "\t" << codesA[k].a1 << "\t" << codesA[k].a2 << "\t" << risk << " (data)\n";
				}
				else
				{
					cout << map[k].rs << "\t" << codesA[k].a1 << "\t" << codesA[k].a2 << "\t" << risk << " (user specified)\n";
				}

				// error message
				if ((SNPrisk[j] != alleleA) && (SNPrisk[j] != alleleB) && (SNPrisk[j] != "R"))
				{
					cout << "\nRisk allele " << SNPrisk[j] << " from " << SNPcovar[j] << " doesn't exit!\n";
					exit(1);
				}
				found=1;break;
			}
		}
		if(!found)
		{
			cout << "\nSNP-Covariate " << j+1 << " was not found or is qc-out. Not used for liability test.\nNote: Individuals which have missing genotype for that SNP are nevertheless excluded from the analysis.\n";
		}
	}
}


int main(int argc, char *argv[])
{
  srand((int)time(NULL));

  int writeOut=0; //writeOut

  int dim1=300; // number of parameters allowed in regGeneral. used to be 27+sexcov=28 until IS 582
  int poschoice_Selected=0;
  int negchoice_Selected=0;
  int snp1_Selected=0;
  int snp2_Selected=0;
  int snp3_Selected=0;
  bool fisherCorrection=0;
  int intervaleditor=0;
  int verbose=1;

  int dim_Cov1=0;int dim_Cov2=0; //dimensions of covariance matrix
  int dim_Cov1_Single=0;int dim_Cov2_Single=0;
  int df_L1=0;
  int df_L2=0;
  int df_L1_Single=0;
  int df_L2_Single=0;

  int *SigmaAuxVector=NULL;
  int *SigmaAuxVector_Single=NULL;

  // switches
  bool bin = false;
  bool vb = false;
  bool vb_binwise_corr = false;
  int stratify = 0;
  string stratify_struct = "";
  string stratify_valid = "VICINITY";
  bool familyCluster = false;
  bool cluster_covars = false;
  int group_test = 0;
  int sim_limit = false;
  bool allAllMatching = false;
  bool annotate = false;
  bool pathwayAnalysis = false;
  bool switchCaseControl = false;
  bool SM_showresults = true;
  bool skip_it_qc = false;

  int lambdaAdjust=0;
  int printBeta=0;
  int choosePrintBeta= 1;
  int covariancematrix=0;

  int filled;
  int passedPreTest=1;
  int dfSingle=1;

  int doIbs = 0;  // IBS handling
  bool ibs_minor_only = false;
  float ibs_SD_relatives = 4.;
  float ibs_SD_outlier   = 4.;
  int ibs_SD_Selected=0;
  int caseOnly=0;
  bool plusSingleFirst=0;
  bool plusSingleSecond=0;
  bool plusSingleThird=0;
  int beginPlusSingle = 0;
  int endPlusSingle = 0;

  double pHelpCut=0.01;
  struct BESTCHI5 x1;
  struct BESTCHI5 y1;
  struct BESTSINGLEMARKER x1sm;
  struct BESTSINGLEMARKER y1sm;
  struct INFOMAP x1infomap;
  struct INFOMAP y1infomap;
  struct MAP2 x1map2;
  struct MAP2 y1map2;
  struct GENETICLIST x1genetic;
  struct GENETICLIST y1genetic;
  struct BESTSINGLEMARKER x1smNR;
  struct BESTSINGLEMARKER y1smNR;
  double pHelp=1;
  struct TRAITAVG traitavg;
  int genecol=5;
  int with25=0;

  int maxposSingle=0;
  int maxposMulti=0;
  double lastpSingle=1.1;
  double lastpMulti=1.1;

  double nbetter=0;
  int pretest=0;
  int kkk=0;
  unsigned int lll=0;
  int ll2=0;
  int jjj=0;
  int mmm=0;
  char hapstring[1000];
  int xsinglevec2[27];
  int zsinglevec2[27];
  int xsinglevec3[27];
  int zsinglevec3[27];
  double cutoff=0.05;
  double p_2nd=1;
  double p_3rd=1;
  struct STATplus result0;
  int skip=0;
  int helpj=0;

  struct MAP *map = NULL;
  struct MAP2 *map2 = NULL;
  struct MAP3 *map3 = NULL;
  struct INFOMAP *infomap = NULL;
  struct COUNTS **counts = NULL;
  struct PERSON *person = NULL;
  struct FAMILY *family = NULL;
  struct REGION *region = NULL;
  struct REGION *regionM = NULL;
  struct BESTSINGLEMARKER *bestsinglemarker = NULL;
  struct BESTSINGLEMARKER *bestsinglemarker2 = NULL; //bestsinglemarker ordered by nr
  struct BESTCHI5 *bestchi3 = NULL;
  struct TSTAT tstat;
  struct TSTAT tstat2; //CASEONLY18
  struct STATplus result1[MAXTHREAD];
  struct STATplus result2[MAXTHREAD];
  struct STATplus result1Single[MAXTHREAD];
  struct STATplus result2Single[MAXTHREAD];
  //struct STATplus result1M[MAXTHREAD];
  //struct STATplus result2M[MAXTHREAD];
  //struct STATplus result1F[MAXTHREAD];
  //struct STATplus result2F[MAXTHREAD];
  struct STATplus result3[MAXTHREAD];
  //struct STATplus result3M[MAXTHREAD];
  //struct STATplus result3F[MAXTHREAD];
  struct STATplus result4[MAXTHREAD];
  //struct STATplus result4M[MAXTHREAD];
  //struct STATplus result4F[MAXTHREAD];
  struct TOPLIST *toplist = NULL;
  struct REGION* chrPositions = NULL;

  struct STATplus resultSingle[MAXTHREAD];
  struct STATplus resultMulti[MAXTHREAD];
  //struct STATplus resultDummy[MAXTHREAD];

  double *bestofSim = NULL;
  struct GENETICLIST *geneticlist = NULL;
  struct GENETICLIST *geneticlist2 = NULL;//ordered by number
  struct OVERLAPLIST *overlaplist = NULL;
  struct PATHWAY *pathway = NULL;
  struct PWgenelist *pwgenelist=NULL; //pplus
  char* *selectedSNP=NULL; // Änderung
  char* *selectedSNPinfo=NULL;
  char* *INFOScore=NULL;
  struct MARKERTABLE **markerTable=NULL;// Änderung569
  int ngenes=0; //pplus
  int genefound=0; //pplus
  int current=0; //pplus
  struct ALLELECODE *codesA = NULL; // Speichern der Allele

  int firstCallThread[MAXTHREAD];

  int maxSize=30; // max strlen of snps in combilist
  int maxSizeNew=30;
  char dummyString[10]="dummy";

  FILE *fptr3 = NULL; // Selectionfile
  FILE *fptr4 = NULL; // Annotationfile
  FILE *fptr5 = NULL; // Pathwayfile
  FILE *fptr6 = NULL; // Covariatefile
  FILE *fptr7 = NULL; // Modelfile
  FILE *fptr8 = NULL; // Combifile
  FILE *fptr9 = NULL; // SNPfile
  FILE *fptr10 = NULL; //Hapfile
  FILE *fptr11 = NULL; // Personfile
  FILE *fptr12 = NULL; // Infofile
  fstream singlemarker; // Singlemarkerfile
  fstream singlemarkerTop; // SinglemarkerTopfile
  fstream rareTop;
  fstream deletedPerson;
  fstream deletedSnps;
  fstream bestMarkerCombi2;
  fstream bestMarkerCombi3;
  fstream bestMarkerCombi2Details;
  fstream bestMarkerCombi3Details;
  fstream toplistMC;
  fstream liabilityResult;
  fstream rare;
  fstream rarepassed;
  fstream missingCombis;
  fstream bestMarkerCombiParallel[MAXTHREAD]; //writeOut
  int raretop=0;
  int bins2testt=0;

  int nlinestped = 0;
  int nlinesy = 0; // nlinestped ohne die Marker auf Y Chromosom
  int ylines = 0;
  int nlinestfam = 0;
  int nlines2Aff = 0;
  int nlines2AffMale = 0;
  int nlinesSelectionfile = 0;
  int ncolumnsSelectionfile = 0;
  int nlinesinfo = 0;
  int ncolumnsinfo = 0;
  int nlinespathway = 0;
  int ncolumnspathway = 0;
  int nlinescovariate = 0;
  int ncolumnscovariate = 0;
  int nlinesmodel = 0;
  int ncolumnsmodel = 0;
  int nlinescombi = 0;
  int ncolumnscombi = 0;
  int nlinesSNP = 0;
  int nlinesSNPinfo = 0;
  int ncolumnsSNP = 0;
  int ncolumnsSNPinfo = 0;
  int nlinesperson = 0;
  int ncolumnsperson = 0;
  int i=0;
  int j=0, k, l;
  int ii = 0;
  int reg = 0;
  int regM = 0;
  int change = 1;
  int choice = 0; // 1.Fall: weder poschoice noch negchoice, 2.Fall: poschoice, 3.Fall: negchoice
  int choiceM = 0;
  int nlinestpedMod = 0;
  int nlinestfamMod = 0;
  int mod = 0;
  int mod2 = 0;
  unsigned int y = 0;
  unsigned int z = 0;
  int a=0, b=0, c=0, d, e;
  int aa, bb, cc, dd, ee;
  int ll, kk;
  int singletop = 10; // NEU singletop wurde auf 10 gesetzt
  int mWithSingletop = 0;
  int markercombi2 = 0;
  int markercombi3 = 0;
  int x = 0;
  int oldCombi = 0;
  int aMod,bMod,cMod;

  int ncases = 0;
  int ncontrols = 0;
  int nrest = 0;
  int ncasesqc = 0;
  int ncontrolsqc = 0;
  int nrestqc = 0;
  uint32_t npplqc = 0;
  int nsnpqc = 0;
  int nsnpAnalysisqc = 0;
  int nsnpMatchingqc = 0;
  int nsnpAnalysisMatchingqc = 0;
  long long int ntests = 0;
  int stopA = 0;
  int stopB = 0;
  int stopC = 0;
  int stopD = 0;
  int stopE = 0;
  int stopCC = 0;
  int stopDD = 0;
  int stopEE = 0;
  int marker1 = 0;
  int marker2 = 0;
  int marker3 = 0;
  int startA = 0;
  int startB = 0;
  int startC = 0;
  int a2 = 0;
  int b2 = 0;
  int c2 = 0;
  int selectedSnp = 0;
  int snp1Pos = -1;
  int snp2Pos = -1;
  int snp3Pos = -1;
  int teststat = 0;
  int test = -1;
  int storeall = 1;
  int nsnps = 0;
  int nsim = 0;
  int ix = 2505;
  int iy = 11831;
  int iz = 23492;
  int helpk = 0;
  int geneticImpact = 0;
  int pathwayImpact = 0;
  int mWithGenetictop = 0;
  int genetictop = 0;
  int jstart = 0;
  int mapsorted = 1;
  int overlaptop = 0;
  int a_list = 0;
  int b_list = 0;
  int c_list = 0;
  int d_list = 0;
  int e_list = 0;
  int r;
  int firstc = 0;
  int storeda = -1;
  int storedb = -1;
  int nlist = 0;
  int snpInPathway1 = 0;
  int snpInPathway2 = 0;
  int snpInPathway3 = 0;
  int snpInPathway = 1;
  int rstart = 0;
  int covariate = 0;
  int xvec[27];
  int zvec[27];
  int xvec2[27]; //CASEONLY18
  int zvec2[27]; //CASEONLY18
  int xsinglevec[27];
  int zsinglevec[27];
  int beginCovar = 0;
  int endCovar = 0;
  int male = 0;
  int female = 0;
  int MCWithSM = 0;
  int readmodel = 0;
  int regression = 0;
  int collapseRare=0;
  int collcollapseRare=0;
  int firstbinlastSNP=0;
  int sexcov = 0;
  int N = 2;
  int N1 = 1;
  int withmarginal = 0;
  int alt;
  int xType = 0;
  int found = 0;
  int rsline = 0;
  int jj = 0;
  int mm = 0;
  int overlapneeded = 0;
  int sum = 0;
  int icov= 0;
  int singleMarkerTest = 1;
  int allequal = 1;
  int combilist = 0;
  int q = 0;
  int qstart = 1;
  int pathwayTest = 2;
  double p_PAA=-1;
  int printtop=100;
  int countA = 0;
  int countB = 0;
  int countC = 0;
  int helpi=0;
  int mlist = 0;
  int pos1=-1;
  int pos2=-1;
  int pos3=-1;
  int snplist = 0;
  int infolist = 0;
  int dosage = 0;
  int setid = 0;
  int personlist= 0;
  int list = 0;
  int list1 = 0;
  int dohapfile = 0;
  int binsizeRare = 0;
  int familyData = 0;
  int snpCov = 0;
  int liabilityCut = 0;
  int liability = 0;
  int riskSNP = 0;
  int corrL = 0;

  double bestL = 1;
  //double snp1 = -1;
  char *snp1=NULL; // Änderung569
  //double snp2 = -1;
  char *snp2=NULL; // Änderung569
  //double snp3 = -1;
  char *snp3=NULL; // Änderung569
  //double mastersnp1 = -1;
  char *mastersnp1=NULL; // Änderung569
  //double mastersnp2 = -1;
  char *mastersnp2=NULL; // Änderung569
  //double mastersnp3 = -1;
  char *mastersnp3=NULL; // Änderung569
  //double refSnp = -1;
  char *refSnp=NULL; // Änderung
  //double **markerTable= NULL; // Änderung569
  //int **markerPosTable= NULL; // Änderung
  int *cov =NULL;

  double mrdiff = 1;
  double hweCa = 0.00;
  double hweCo = 0.00;
  double testHweCa;
  double testHweCo;
  double p_interratio = 0.05;
  double meanMissingPerson = 0;
  double meanMissingSNP = 0;
  double meanMissingNew = 0;
  double newValue = 0;
  double newValueCorr = 0;
  double newValueTop = 0;
  double helpstat = 0;
  double casecounts5[3][3][3][3][3];
  double casecounts5M[3][3][3][3][3];
  double casecounts5F[3][3][3][3][3];
  double controlcounts5[3][3][3][3][3];
  double controlcounts5M[3][3][3][3][3];
  double controlcounts5F[3][3][3][3][3];
  double inflationfactor = 0; // Durchschnitt der Teststatistik vom Singlemarker (bei df = 1 ~ 1)
  double fulltests = 0;
  double val1 = -999.0;
  double val2 = -999.0;
  double ntestsPlus= 0;
  double *countsP=NULL;
  double ntestsUser = 0;
  double pfilter = 1;
  double probRequired = 0;
  double infoScore = 0.4;

  // Logistische Regression
  double *p = NULL;
  double **X = NULL;
  double **Xt = NULL;
  double **A = NULL;
  double **VNN = NULL;
  double *S = NULL;
  double **Sinv = NULL;
  double **A0 = NULL;
  double **UNNT = NULL;
  double **Ainv = NULL;
  double **AinvXt = NULL;
  double **Yminusp = NULL;
  double **newbeta = NULL;
  double **Xmod = NULL;
  double *YY = NULL;
  double *Yhelp=NULL;
  double **Yt=NULL;
  double **YtX=NULL;
  double **YtXAinv=NULL;
  double **YtXAinvXt=NULL;
  double **D=NULL;
  double **T=NULL;
  double **U=NULL;
  double **Ut=NULL;
  double **sumPP=NULL;
  double **sumPJ=NULL;
  double **sumPK=NULL;
  double **MMinv=NULL;

  //MC-procdure
  int **Y = NULL;
  double **Yd = NULL;
  double **YY1 = NULL;
  double **YY2 = NULL;

  double tLog = 1;
  double tLogMale = 1;
  double tLogFemale = 1;
  double df = 0;
  double dfMale = 0;
  double dfFemale = 0;

  char s3[100000];
  char s4[100000];
  char s5[100000];
  char s6[100000];
  char s7[100000];
  char s8[100000];
  char s9[100000];
  char s10[100000];
  char s11[100000];
  char character;
  char *str3 = NULL;
  char *str44 = NULL;
  char *str55 = NULL;
  char *str66 = NULL;
  char *str77 = NULL;
  char *str88 = NULL;
  char *str99 = NULL;
  char *str33 = NULL;
  char *str22 = NULL;
  char *str11 = NULL;
  char *genestring = NULL; //pplus
  char *snpstring=NULL;
  int newgene=0; //pplus
  char *str_fid=NULL; //NEU
  char *str_pid=NULL; //NEU
  char xstring[2]="x"; //NEU
  char minstring[2]="-"; //NEU

  string singlemarkerfile = " ";
  string singlemarkerTopfile = " ";
  string rareTopfile = " ";
  string liabilityfile = " ";
  string markerCombi2file = " ";
  string markerCombifileParallel[MAXTHREAD]; //writeOut
  string markerCombi2Detailsfile = " ";
  string markerCombi3file = " ";
  string markerCombi3Detailsfile = " ";
  string errorname = "errorMessage.txt";
  string logname = "logfile.txt";

  string delPers = " ";
  string delSnps = " ";
  string toplistMCfile = " ";
  string rarefile = " ";
  string rarefileinter = " ";   // for RAREINTER
  string rarefileVB = " ";
  string rarefileVBpermstat = " ";
  string rarefileVBgraph = " ";
  string rarefileVB_BMP = " ";
  string SetIDfile = " ";
  string keyword;
  string info;
  string missCombis = " ";

  string famfile  = "";
  string mapfile  = "";
  string bmapfile = "";
  string pedfile  = "";
  string tpedfile = "";
  string bpedfile = "";
  string impute2file = "";
  uint8_t filemode = 0;

  string annotationfile = " ";
  string pathwayfile = " ";
  string covariatefile = " ";
  string modelfile = " ";
  string combifile = " ";
  string SNPfile = " ";
  string personfile = " ";
  string infofile = " ";
  string str4;
  string snp1help;
  string snp2help;
  string snp3help;
  string maprs;
  string infomaprs;
  string emptystring = " ";
  string rsHelpmap;
  string rsStr = " ";
  char *hapfile=NULL;
  string intervalfile= " ";

  // GFF
  string intervalfile_format= " ";
  string gffname= " ";
  vector<string> selectedChr;
  vector<string> selectedregion;
  vector<string> selectedregion2;
  vector<string> selectedCovar;
  vector<string> covar;
  vector<string> SNPcovar;
  vector<string> SNPrisk;
  vector<string> selectedPlussingle;
  vector<string> plussingleSNPs;
  vector<string> selectedRareTests;
  vector<string> weightParams;

  bool ifChrX = 0;

  int qt=0;
  double missingpheno=-9;
  double Fstat=0;


  int ncov=0;
  int maxIndexCov = 0;
  int numberOfSNPCov = 0;
  int haplo=0;
  double distance=30000;
  int LDdistance=50000000;

  int plusSingle = 0;

  //PAA
  double **Ttable=NULL;
  double **hTable=NULL; //one row index, one row order

  //MAF
  double maf=0;
  double thismaf=0;
  double thiscontrolmaf=0;
  double thiscasemaf=0;

  double thismafa=0;
  double thiscontrolmafa=0;
  double thiscasemafa=0;

  double thismafr=0;
  double thiscontrolmafr=0;
  double thiscasemafr=0;

  int c_in=0;
  int d_in=0;
  int e_in=0;

  //PLL
  int plla2 = 0;
  int plln = 0;
  long int n=0;
  int maxthreads=1;
  int thread=0;
  int thread_nloop=0;
  int thread_aloop=0;

  stringstream sstm;
  uint32_t nwordsSNPs = 0;
  uint32_t nwordsPPLs = 0;
  uint64_t*** BinSNPs = NULL;
  uint64_t*** BinPPLs = NULL;
  uint64_t**  BinSNPsGenderFlags = NULL;
  uint64_t**  BinSNPsCCFlags     = NULL;
  uint64_t*   BinSNPsQTaffFlags  = NULL;
  uint64_t*   BinSNPsQCinFlags   = NULL;
  int* SNPMap = NULL;
  int* PPLMap = NULL;
  int* SNPMapInverse = NULL;
  int* PPLMapInverse = NULL;
  struct PPLLOCATION* PplLocations = NULL;
  struct PPLLOCATION* PplLocFam = NULL;
  struct RELATIVES* MatchedPairs = NULL;
  struct CLUSTERS* Clusters = NULL;
  uint32_t nMatches  = 0;
  uint32_t nClusters = 0;
  uint32_t nFamilies = 0;
  getIbsMatrix = getTotalIbsMatrix;
  getIbs = getTotalIbs;

  // impute2
  float ***genoWeights = NULL;

#if PARALLEL
  struct INDEX index;
  double stat=0;
#endif

#if RARE
  struct WEIGHTS weights;
  weights.mode=0;
  bool FISHERtest=0;
  bool REGRESSIONtest=0;
  bool FRACREGtest=0;
  bool COLLREGtest=0;
  bool COLLtest=0;
  bool CMATtest=0;

  int nCovCathegories=0; // for CMAT with covariates
  int *CovCathegories=NULL;// for CMAT with covariates

  bool optimalrare = false;
  bool vb_bmp = false;
  bool vb_print_perm_stat = false;
  bool vb_bmp_level = false;
  int nvbstart=0;
  int *nvbend=NULL;
  int *vbstart=NULL;
  int **vbend=NULL;
  float *vbmaxpermstat=NULL;
  int **nvblevel=NULL;
  int ***nchunks=NULL;
  int ****chunkpos=NULL;
  int ****chunklen=NULL;

  int *vbnchunks=NULL;
  int **vbchunkpos=NULL;
  int **vbchunklen=NULL;

  float *****vbbinwisestat=NULL;
  short int *****vbbinwisecount=NULL;

  uint64_t**  BinSNPsCCFlagsOriginal = NULL;

  int *nvbstartvt=NULL;
  int **vbstartvt=NULL;
  int ****ndummyatlevel=NULL;
  int *****dummylevel=NULL;
  int *****dummypos=NULL;
  int *****dummycluster=NULL;

  int ****dummyend=NULL;
  int ***ndummyends=NULL;
  int ***nchunkcluster=NULL;

  uint64_t ***BinCarriers=NULL;

  int ****vbendvt=NULL;
  int **rareLimitsNCT=NULL;
  int **rareLimitsNCTInverse=NULL;
  //  float *vbstatvt=NULL;
  //  float *vbmaxpermstatvt=NULL;


  int minIndInBin=20;
  int maxIndInBin=ncases+ncontrols-20;


  int binadjust=0;
  int collinter=0;
  int comparemode=0;
  bool rare_stratify = false;
  uint32_t rare_strat_min = 10000;
  float rare_strat_maft = -1.0;
  double **limitsStatCOLL=NULL;
  double **limitsStatCMAT=NULL;
  double **limitsStatFISHER=NULL;
  double **limitsStatREGRESSION=NULL;
  double **limitsStatFRACREG=NULL;
  double **limitsStatCOLLREG=NULL;

  int *nRareLimits=NULL;
  double **rareLimits=NULL;

  double raref = 0;
  int NCT = 0;
  int nwindows=0;
  int nwindowssinglebin=0;
  int **nchrwindows=NULL;
  LOCAL_MATCHING* MatchingRARE = NULL;
  uint32_t nMWindowsRARE = 0;
#endif


#if PARALLELN
  plln=1;
  maxthreads=omp_get_max_threads();
  if(maxthreads>MAXTHREAD)
    {
      maxthreads=MAXTHREAD;
    }
  omp_set_num_threads(maxthreads);
  cout << "maxthreads " << maxthreads << "\n";
  thread=omp_get_thread_num();
  cout << "thread " << thread << "\n";
#endif

#if PARALLELA
  plla2=1;
  maxthreads=omp_get_max_threads();
  if(maxthreads>MAXTHREAD)
    {
      maxthreads=MAXTHREAD;
    }
  omp_set_num_threads(maxthreads);
  cout << "maxthreads " << maxthreads << "\n";
  thread=omp_get_thread_num();
  cout << "thread " << thread << "\n";
#endif


  if (argc < 2)
    {
      errorfile.open("errorMessage.txt", ios::out);
      logfile.open("logfile.txt", ios::out);
      errorfile << "Please specify selectionfile!\n";
      logfile << "Please specify selectionfile!\n";
      cout << "Please specify selectionfile!\n";
      errorfile.close();logfile.close();exit(1);
    }

  fptr3 = fopen(argv[1], "r"); // Öffnen der Datei
  if (fptr3 == NULL)
    {
      errorfile.open("errorMessage.txt", ios::out);
      logfile.open("logfile.txt", ios::out);
      errorfile << "The selection file " << argv[1] << " does not exist\n";
      logfile << "The selection file " << argv[1] << " does not exist\n";
      cout << "The selection file " << argv[1] << " does not exist\n";
      errorfile.close();logfile.close();exit(1);
    }
  fclose(fptr3);

  ncolumnsSelectionfile = countColumns(argv[1]);

  if (ncolumnsSelectionfile == -1)
    {

      errorfile << "File " << argv[1] << " does not exist\n";
      logfile << "File " << argv[1] << " does not exist\n";
      cout << "File " << argv[1] << " does not exist\n";
      errorfile.close();
      logfile.close();
      exit(1);
    }
  else if (ncolumnsSelectionfile == -2)
    {
      errorfile << "First line of " <<  argv[1] << " is empty! Please check infile.\n";
      logfile << "First line of " <<  argv[1] << " is empty! Please check infile.\n";
      cout << "First line of " <<  argv[1] << " is empty! Please check infile.\n";
      errorfile.close();
      logfile.close();
      exit(1);
    }

  hapfile = (char *) calloc(15,sizeof(char));
  sprintf(hapfile,"IS_hapfile.txt");

  /// read selection-file for outputname
  character = '0';
  i = 0;
  choice = 0;
  fptr3 = fopen(argv[1], "r");
  while ((feof(fptr3)) == 0)
    {
      character = fgetc(fptr3);

      if (character == '\n' && i > 1)
	{
	  s3[i] = '\0';
	  i = 0;
	  nlinesSelectionfile++;

	  for (j = 0; j < ncolumnsSelectionfile; j++) // Speichern der einzelnen Spalten
	    {
	      if (j == 0) // Schlüsselwort
		{
		  str3 = strtok(s3, " \t");
		  keyword = str3;
		  //printf("str3 %s\n", str3);
		  if (!str3)
		    {
		      errorfile.open("errorMessage.txt", ios::out);
		      logfile.open("logfile.txt", ios::out);
		      errorfile << "Selection file: Not enough columns in line " << nlinesSelectionfile << "\n";
		      logfile << "Selection file: Not enough columns in line " << nlinesSelectionfile << "\n";
		      cout << "Selection file: Not enough columns in line " << nlinesSelectionfile << "\n";
		      errorfile.close();logfile.close();exit(1);
		    }
		}
	      else
		{

		  str3 = strtok(NULL, " \t");

		  if (str3 != NULL && str3[0]!=' ' && str3[0]!='\t' && str3[0]!='\r' )
		    {
		      if (j == 1) // Parameter
			{
			  if (keyword == "OUTPUTNAME"||keyword == "OUTPUT")
			    {
			      outputname = str3;
			      hapfile = (char *) realloc(hapfile, (strlen(str3)+15)*sizeof(char));
			      sprintf(hapfile,"%shapfile.txt",str3);

			      if (outputname == "//")
				{
				  outputname = "";
				  sprintf(hapfile,"hapfile.txt");
				}
			    }
			}
		    }//j>0
		}
	    }// end j
	} //neue zeile
      else
	{
	  if (character != '\r')
	    {
	      s3[i++] = character;
	    }
	}
    } //end while
  fclose(fptr3);


  errorname = outputname + "ERRORS.txt";
  logname = outputname + "log.txt";
  delPers = outputname + "deletedIndividuals.txt";
  delSnps = outputname + "deletedSnps.txt";
  missCombis = outputname + "missingCombis.txt";

  errorfile.open(errorname.c_str(), ios::out);
  logfile.open(logname.c_str(), ios::out);
  errorfile.clear();
  logfile.clear();

  sstm << "\n"
       << "*---------------------------------------------------*\n"
#if !RARE
       << "|    INTERSNP     |   v1.15    |    22/01/2015      |\n"
#endif
#if RARE
       << "|  INTERSNP-RARE  |   v1.0     |    22/01/2015      |\n"
#endif
       << "*---------------------------------------------------*\n"
    //         << "| (c) Tim Becker                                    |\n"
    //         << "|     Dmitriy Drichel                               |\n"
    //         << "|     Christine Herold                              |\n"
    //         << "|     Andre Lacour                                  |\n"
    //         << "|     Tatsiana Vaitsiakhovich                       |\n"
    //         << "|     Vitalia Schueller                             |\n"
    //         << "*---------------------------------------------------*\n"
    ;
  logg(sstm);

  if (maxthreads > MAXTHREAD)
    {
      errorfile << "Too many threads detected with omp_get_num_thread. Alter define MAXTHREAD statement and recompile!\n";
      logfile << "Too many threads with omp_get_num_thread. Alter define MAXTHREAD statement and recompile!\n";
      cout << "TToo many threads with omp_get_num_thread. Alter define MAXTHREAD statement and recompile!\n";
      errorfile.close();logfile.close();exit(1);
    }

#if defined __WINDOWS__ || defined __DOS__ || defined _WIN32 || defined _WIN64 || defined _WINDOWS || defined _MSDOS || defined _DOS
  cout    << "WARNING: If you use INTERSNP on DOS/WINDOWS systems together with PLINK files, proceed the files with the PLINK DOS version first!!\n";
  logfile << "WARNING: If you use INTERSNP on DOS/WINDOWS systems together with PLINK files, proceed the files with the PLINK DOS version first!!\n";
#endif

  fptr3 = fopen(argv[1], "r"); // Öffnen der Datei
  if (fptr3 == NULL)
    {
      errorfile << "The selection file " << argv[1] << " does not exist\n";
      logfile << "The selection file " << argv[1] << " does not exist\n";
      cout << "The selection file " << argv[1] << " does not exist\n";
      errorfile.close();
      logfile.close();
      exit(1);
    }

  sstm << "Reading selection-file " << argv[1] << " ...";
  logg(sstm);
  logg("Logfile is " + logname + " .");
  fclose(fptr3);


  /// read selectionfile again
  character = '0';
  i = 0;
  choice = 0;
  nlinesSelectionfile=0;
  fptr3 = fopen(argv[1], "r");
  while ((feof(fptr3)) == 0)
    {
      character = fgetc(fptr3);

      if (character == '\n' && i > 1)
	{
	  s3[i] = '\0';
	  i = 0;
	  nlinesSelectionfile++;

	  for (j = 0; j < ncolumnsSelectionfile; j++) // Speichern der einzelnen Spalten
	    {
	      if (j == 0) // Schlüsselwort
		{
		  str3 = strtok(s3, " \t");
		  keyword = str3;
		  if (!str3)
		    {
		      errorfile << "Selection file: Not enough columns in line " << nlinesSelectionfile << "\n";
		      logfile << "Selection file: Not enough columns in line " << nlinesSelectionfile << "\n";
		      cout << "Selection file: Not enough columns in line " << nlinesSelectionfile << "\n";
		      errorfile.close();
		      logfile.close();
		      exit(1);
		    }
		}
	      else
		{
		  str3 = strtok(NULL, " \t");

		  if (str3 != NULL && str3[0]!=' ' && str3[0]!='\t' && str3[0]!='\r' )
		    {
		      if (j == 1) // Parameter
			{
			  if (keyword == "FILE")
			    {
			      string str = str3;
			      if (mapfile=="") mapfile = str + ".map";
			      if (pedfile=="") pedfile = str + ".ped";
			      logfile << keyword << "\t" << str << "\n";
			      filemode = 1;
			    }
			  if (keyword == "TFILE")
			    {
			      string str = str3;
			      if ( famfile=="")  famfile = str + ".tfam";
			      if (tpedfile=="") tpedfile = str + ".tped";
			      logfile << keyword << "\t" << str << "\n";
			      filemode = 2;
			    }
			  if (keyword == "BFILE")
			    {
			      string str = str3;
			      if ( famfile=="")  famfile = str + ".fam";
			      if (bmapfile=="") bmapfile = str + ".bim";
			      if (bpedfile=="") bpedfile = str + ".bed";
			      logfile << keyword << "\t" << str << "\n";
			      filemode = 3;
			    }
			  if (keyword == "IFILE")
			    {
			      string str = str3;
			      if (mapfile=="") mapfile = str + ".map";
			      if ( famfile=="")  famfile = str + ".fam";
			      if (impute2file=="") impute2file = str + ".impute2";
			      logfile << keyword << "\t" << str << "\n";
			      filemode = 4;
			    }
			  else if (keyword == "FAM" || keyword == "TFAM" || keyword == "BFAM")
			    {
			      famfile = str3;
			      logfile << keyword << "\t" << famfile.c_str() << "\n";
			    }
			  else if (keyword == "MAP")
			    {
			      mapfile = str3;
			      logfile << keyword << "\t" << mapfile.c_str() << "\n";
			    }
			  else if (keyword == "BMAP" || keyword == "BIM")
			    {
			      bmapfile = str3;
			      logfile << keyword << "\t" << bmapfile.c_str() << "\n";
			    }
			  else if (keyword == "PED")
			    {
			      pedfile = str3;
			      logfile << keyword << "\t" << pedfile.c_str() << "\n";
			      filemode = 1;
			    }
			  else if (keyword == "TPED")
			    {
			      tpedfile = str3;
			      logfile << keyword << "\t" << tpedfile.c_str() << "\n";
			      filemode = 2;
			    }
			  else if (keyword == "BPED" || keyword == "BED")
			    {
			      bpedfile = str3;
			      logfile << keyword << "\t" << bpedfile.c_str() << "\n";
			      filemode = 3;
			    }
			  else if (keyword == "IMPUTE2")
			    {
			      impute2file = str3;
			      logfile << keyword << "\t" << impute2file.c_str() << "\n";
			      filemode = 4;
			    }
			  else if (keyword == "ANNOTATIONFILE")
			    {
			      annotationfile = str3;
			      if (annotationfile == "//")
				{
				  annotationfile = "";
				}
			      logfile << "ANNOTATIONFILE\t" << annotationfile.c_str() << "\n";

			    }
			  else if (keyword == "PATHWAYFILE")
			    {
			      pathwayfile = str3;
			      if (pathwayfile == "//")
				{
				  pathwayfile = "";
				}
			      logfile << "PATHWAYFILE\t" << pathwayfile.c_str() << "\n";
			    }
			  else if (keyword == "COVARIATEFILE")
			    {
			      covariatefile = str3;
			      if (covariatefile == "//")
				{
				  covariatefile = "";
				}
			      logfile << "COVARIATEFILE\t" << covariatefile.c_str() << "\n";
			    }
			  else if (keyword == "MODELFILE")
			    {
			      modelfile = str3;
			      if (modelfile == "//")
				{
				  modelfile = "";
				}
			      logfile << "MODELFILE\t" << modelfile.c_str() << "\n";
			    }
			  else if (keyword == "COMBIFILE")
			    {
			      combifile = str3;
			      if (combifile == "//")
				{
				  combifile = "";
				}
			      logfile << "COMBIFILE\t" << combifile.c_str() << "\n";

			    }
			  else if (keyword == "SNPFILE")
			    {
			      SNPfile = str3;
			      if (SNPfile == "//")
				{
				  SNPfile = "";
				}
			      logfile << "SNPFILE\t" << SNPfile.c_str() << "\n";
			    }
			  else if (keyword == "PERSONFILE")
			    {
			      personfile = str3;
			      if (personfile == "//")
				{
				  personfile = "";
				}
			      logfile << "PERSONFILE\t" << personfile.c_str() << "\n";
			    }
			  else if (keyword == "INFOFILE")
			    {
			      infofile = str3;
			      if (infofile == "//")
				{
				  infofile = "";
				}
			      logfile << "INFOFILE\t" << infofile.c_str() << "\n";
			    }
			  else if (keyword == "ONLY_MALE")
			    {
			      male = atoi(str3);
			      logfile << "ONLY_MALE\t" << male << "\n";
			      if (male > 1)
				{
				  errorfile << "For parameter ONLY_MALE only 0 or 1 is allowed.\n";
				  logfile << "For parameter ONLY_MALE only 0 or 1 is allowed.\n";
				  cout << "For parameter ONLY_MALE only 0 or 1 is allowed.\n";
				  logfile.close();
				  errorfile.close();
				  exit(1);
				}
			    }
			  else if (keyword == "WRITE_OUT") //writeOut
			    {
			      writeOut = atoi(str3);
			      logfile << "WRITE_OUT\t" << writeOut << "\n";
			      if (writeOut > 1)
				{
				  errorfile << "For parameter WRITE_OUT only 0 or 1 is allowed.\n";
				  logfile << "For parameter WRITE_OUT only 0 or 1 is allowed.\n";
				  cout << "For parameter WRITE_OUT only 0 or 1 is allowed.\n";
				  logfile.close();errorfile.close();exit(1);
				}
			    }
			  else if (keyword == "ONLY_FEMALE")
			    {
			      female = atoi(str3);
			      logfile << "ONLY_FEMALE\t" << female << "\n";
			      if (female > 1)
				{
				  errorfile << "For parameter ONLY_FEMALE only 0 or 1 is allowed.\n";
				  logfile << "For parameter ONLY_FEMALE only 0 or 1 is allowed.\n";
				  cout << "For parameter ONLY_FEMALE only 0 or 1 is allowed\n";
				  errorfile.close();
				  logfile.close();
				  exit(1);
				}
			    }
			  else if (keyword == "COMBILIST")
			    {
			      combilist = atoi(str3);
			      if(combilist>0){combilist=1;}
			      logfile << "COMBILIST\t" << combilist << "\n";
			    }
			  else if (keyword == "SNPLIST")
			    {
			      snplist = atoi(str3);
			      if (snplist > 0)
				{
				  choice = 3;
				}
			      if(setid == 1)
				{
				  cout<<"SETID is not compatible with SNPLIST. SNPLIST will be ignored.\n";
				  logfile<<"SETID is not compatible with SNPLIST. SNPLIST will be ignored.\n";
				}
			      logfile << "SNPLIST\t" << snplist << "\n";
			    }
			  else if (keyword == "INFOLIST")
			    {
			      infolist = atoi(str3);
			      if (infolist > 0) infolist=1;
			      logfile << "INFOLIST\t" << infolist << "\n";
			    }
			  else if (keyword == "INFO_SCORE")
			    {
			      infoScore = atof(str3);
			      logfile << "INFO_SCORE\t" << infoScore << "\n";
			      if (infoScore < 0 || infoScore >1)
				{
				  errorfile << "INFO_SCORE takes only values between 0 and 1!\n";
				  logfile << "INFO_SCORE takes only values between 0 and 1!\n";
				  cout << "INFO_SCORE takes only values between 0 and 1!\n";
				  errorfile.close();
				  logfile.close();
				  exit(1);

				}
			    }
			  else if (keyword == "SETID")
			    {
			      if(str3 != "//" && str3!="")
				{
				  SNPfile = str3;
				  SetIDfile=str3;
				  if(!file_exists(SetIDfile)){
				    cout<<SetIDfile<<" Not found!"<<endl;
				    logfile<<SetIDfile<<" Not found!"<<endl;
				    errorfile<<SetIDfile<<" Not found!"<<endl;
				    exit(1);
				  }
				  else{
				    snplist=1;
				    setid=1;
				  }
				}
			      choice = 3;
			      logfile << "SETID\t" << str3 << "\n";
			    }
			  else if (keyword == "PERSONLIST")
			    {
			      personlist = atoi(str3);

			      logfile << "PERSONLIST\t" << personlist << "\n";
			    }
			  else if (keyword == "POSCHOICE" && str3[0]!='/')
			    {
			      poschoice_Selected=1;
			      if (snplist == 1)
				{
				  errorfile << "Cannot select SNPLIST and POSCHOICE.\n";
				  logfile << "Cannot select SNPLIST and POSCHOICE.\n";
				  cout << "Cannot select SNPLIST and POSCHOICE.\n";
				  errorfile.close();logfile.close();exit(1);
				}

			      logfile << "POSCHOICE\t" << str3 << "\n";
			      if (str3 != NULL)
				{
				  choice = 1;

				  // Eingabe nach Simikolon trennen
				  str4 = str3;
				  selectedChr = split(str4, ";");


				  for (y = 0; y < selectedChr.size(); y++)
				    {
				      if (selectedChr[y] == "X")
					{
					  selectedChr[y] = "23";
					}
				      else if (selectedChr[y] == "Y")
					{
					  selectedChr[y] = "24";
					}
				      else if (selectedChr[y] == "YX")
					{
					  selectedChr[y] = "25";
					}
				      else if (selectedChr[y] == "Mt" || selectedChr[y] == "MT" )
					{
					  selectedChr[y] = "26";
					}
				      // Speicher anlegen
				      region = (struct REGION *) realloc(region, (reg + 1) * sizeof(struct REGION));
				      if (!region)
					{
					  errorfile << "memory allocation error in region\n";
					  logfile << "memory allocation error in region\n";
					  cout << "memory allocation error in region\n";
					  errorfile.close();
					  logfile.close();
					  exit(1);
					}

				      // Teilstring nach "," trennen
				      selectedregion = split(selectedChr[y],",");

				      if (selectedregion.size() == 1)
					{
					  if (selectedregion[0].size() == 4)
					    {
					      region[reg].chr[0] = selectedregion[0][3];
					      region[reg].chr[1] = '\0';
					    }
					  else if (selectedregion[0].size()== 5)
					    {
					      region[reg].chr[0] = selectedregion[0][3];
					      region[reg].chr[1] = selectedregion[0][4];
					      region[reg].chr[2] = '\0';
					    }
					  region[reg].begin = 0;
					  region[reg].end = 1000000000;

					  reg++;


					}
				      else if (selectedregion.size() > 1)
					{
					  for (z = 1; z < selectedregion.size(); z++)
					    {
					      // Teilstring nach "-" trennen
					      selectedregion2 = split(selectedregion[z], "-");

					      if (selectedregion2.size() == 1)
						{
						  cout << "Usage POSCHOICE:\nchr1;\nchr1;chr7;\nchr1,100-20000;\nchr1,100-20000;chr7,8000-200000;\nchr1,100-20000;chr7,8000-200000;chr7\n";
						  errorfile << "Usage POSCHOICE:\nchr1;\nchr1;chr7;\nchr1,100-20000;\nchr1,100-20000;chr7,8000-200000;\nchr1,100-20000;chr7,8000-200000;chr7\n";
						  logfile << "Usage POSCHOICE:\nchr1;\nchr1;chr7;\nchr1,100-20000;\nchr1,100-20000;chr7,8000-200000;\nchr1,100-20000;chr7,8000-200000;chr7\n";
						  errorfile.close();
						  logfile.close();
						  exit(1);
						}

					      if (selectedregion[0].size()== 4)
						{
						  region[reg].chr[0] = selectedregion[0][3];
						  region[reg].chr[1] = '\0';
						}
					      else if (selectedregion[0].size()== 5)
						{
						  region[reg].chr[0] = selectedregion[0][3];
						  region[reg].chr[1] = selectedregion[0][4];
						  region[reg].chr[2] = '\0';
						}
					      region[reg].begin = (int) atof(selectedregion2[0].c_str());
					      region[reg].end = (int) atof(selectedregion2[1].c_str());
					      reg++;
					    }
					}
				    }
				}

			    }
			  else if (keyword == "NEGCHOICE" && str3[0]!='/')
			    {
			      negchoice_Selected=1;
			      if (snplist == 1)
				{
				  errorfile << "Cannot select SNPLIST and NEGCHOICE.\n";
				  logfile << "Cannot select SNPLIST and NEGCHOICE.\n";
				  cout << "Cannot select SNPLIST and NEGCHOICE.\n";
				  errorfile.close();logfile.close();exit(1);
				}


			      logfile << "NEGCHOICE\t" << str3 << "\n";
			      choice = 2;

			      // Eingabe nach Simikolon trennen
			      str4 = str3;
			      selectedChr = split(str4, ";");

			      for (y = 0; y < selectedChr.size(); y++)
				{
				  if (selectedChr[y] == "X")
				    {
				      selectedChr[y] = "23";
				    }
				  else if (selectedChr[y] == "Y")
				    {
				      selectedChr[y] = "24";
				    }
				  else if (selectedChr[y] == "YX")
				    {
				      selectedChr[y] = "25";
				    }
				  else if (selectedChr[y] == "Mt" || selectedChr[y] == "MT" )
				    {
				      selectedChr[y] = "26";
				    }

				  // Speicher anlegen
				  region = (struct REGION *) realloc(region,
								     (reg + 1) * sizeof(struct REGION));
				  if (!region)
				    {
				      errorfile << "memory allocation error in region\n";
				      logfile << "memory allocation error in region\n";
				      cout << "memory allocation error in region\n";
				      errorfile.close();
				      logfile.close();
				      exit(1);
				    }

				  // Teilstring nach "," trennen
				  selectedregion = split(selectedChr[y], ",");

				  if (selectedregion.size() == 1)
				    {
				      if (selectedregion[0].size() == 4)
					{
					  region[reg].chr[0] = selectedregion[0][3];
					  region[reg].chr[1] = '\0';
					}
				      else if (selectedregion[0].size()== 5)
					{
					  region[reg].chr[0] = selectedregion[0][3];
					  region[reg].chr[1] = selectedregion[0][4];
					  region[reg].chr[2] = '\0';
					}
				      region[reg].begin = 0;
				      region[reg].end = 1000000000;
				      reg++;
				    }
				  else if (selectedregion.size() > 1)
				    {
				      for (z = 1; z < selectedregion.size(); z++)
					{
					  // Teilstring nach "-" trennen
					  selectedregion2 = split(selectedregion[z], "-");

					  if (selectedregion2.size() == 1)
					    {
					      cout << "Usage NEGCHOICE:\nchr1;\nchr1;chr7;\nchr1,100-20000;\nchr1,100-20000;chr7,8000-200000;\nchr1,100-20000;chr7,8000-200000;chr7\n";
					      logfile << "Usage NEGCHOICE:\nchr1;\nchr1;chr7;\nchr1,100-20000;\nchr1,100-20000;chr7,8000-200000;\nchr1,100-20000;chr7,8000-200000;chr7\n";
					      errorfile << "Usage NEGCHOICE:\nchr1;\nchr1;chr7;\nchr1,100-20000;\nchr1,100-20000;chr7,8000-200000;\nchr1,100-20000;chr7,8000-200000;chr7\n";
					      errorfile.close();
					      logfile.close();
					      exit(1);
					    }
					  if (selectedregion[0].size() == 4)
					    {
					      region[reg].chr[0] = selectedregion[0][3];
					      region[reg].chr[1] = '\0';
					    }
					  else if (selectedregion[0].size() == 5)
					    {
					      region[reg].chr[0] = selectedregion[0][3];
					      region[reg].chr[1] = selectedregion[0][4];
					      region[reg].chr[2] = '\0';
					    }
					  region[reg].begin = (int) atof(selectedregion2[0].c_str());
					  region[reg].end = (int) atof(selectedregion2[1].c_str());
					  reg++;
					}
				    }
				}
			    }
			  else if (keyword == "MATCHINGCHOICE" && str3[0]!='/')
			    {
			      logfile << "MATCHINGCHOICE\t" << str3 << "\n";
			      if (str3 != NULL)
				{
				  choiceM = 1;
				  // Eingabe nach Simikolon trennen
				  str4 = str3;
				  selectedChr = split(str4, ";");

				  for (y = 0; y < selectedChr.size(); y++)
				    {
				      if (selectedChr[y] == "X")
					{
					  selectedChr[y] = "23";
					}
				      else if (selectedChr[y] == "Y")
					{
					  selectedChr[y] = "24";
					}
				      else if (selectedChr[y] == "YX")
					{
					  selectedChr[y] = "25";
					}
				      else if (selectedChr[y] == "Mt" || selectedChr[y] == "MT" )
					{
					  selectedChr[y] = "26";
					}
				      // Speicher anlegen
				      regionM = (struct REGION*) realloc(regionM, (regM + 1) * sizeof(struct REGION));
				      if (!regionM)
					{
					  errorfile << "memory allocation error in regionM\n";
					  logfile << "memory allocation error in regionM\n";
					  cout << "memory allocation error in regionM\n";
					  errorfile.close();
					  logfile.close();
					  exit(1);
					}
				      // Teilstring nach "," trennen
				      selectedregion = split(selectedChr[y],",");

				      if (selectedregion.size() == 1)
					{
					  if (selectedregion[0].size() == 4)
					    {
					      regionM[regM].chr[0] = selectedregion[0][3];
					      regionM[regM].chr[1] = '\0';
					    }
					  else if (selectedregion[0].size()== 5)
					    {
					      regionM[regM].chr[0] = selectedregion[0][3];
					      regionM[regM].chr[1] = selectedregion[0][4];
					      regionM[regM].chr[2] = '\0';
					    }
					  regionM[regM].begin = 0;
					  regionM[regM].end = 1000000000;
					  regM++;
					}
				      else if (selectedregion.size() > 1)
					{
					  for (z = 1; z < selectedregion.size(); z++)
					    {
					      // Teilstring nach "-" trennen
					      selectedregion2 = split(selectedregion[z], "-");

					      if (selectedregion2.size() == 1)
						{
						  cout << "Usage POSCHOICE:\nchr1;\nchr1;chr7;\nchr1,100-20000;\nchr1,100-20000;chr7,8000-200000;\nchr1,100-20000;chr7,8000-200000;chr7\n";
						  errorfile << "Usage POSCHOICE:\nchr1;\nchr1;chr7;\nchr1,100-20000;\nchr1,100-20000;chr7,8000-200000;\nchr1,100-20000;chr7,8000-200000;chr7\n";
						  logfile << "Usage POSCHOICE:\nchr1;\nchr1;chr7;\nchr1,100-20000;\nchr1,100-20000;chr7,8000-200000;\nchr1,100-20000;chr7,8000-200000;chr7\n";
						  errorfile.close();
						  logfile.close();
						  exit(1);
						}

					      if (selectedregion[0].size()== 4)
						{
						  regionM[regM].chr[0] = selectedregion[0][3];
						  regionM[regM].chr[1] = '\0';
						}
					      else if (selectedregion[0].size()== 5)
						{
						  regionM[regM].chr[0] = selectedregion[0][3];
						  regionM[regM].chr[1] = selectedregion[0][4];
						  regionM[regM].chr[2] = '\0';
						}
					      regionM[regM].begin = (int) atof(selectedregion2[0].c_str());
					      regionM[regM].end = (int) atof(selectedregion2[1].c_str());
					      regM++;
					    }
					}
				    }
				}

			    }
			  else if (keyword == "HWE_P_CASE")
			    {
			      hweCa = atof(str3);
			      if (strcmp (str3, "//") == 0)
				{
				  hweCa = 0.00;
				}
			      logfile << "HWE_P_CASE\t" << hweCa << "\n";
			    }
			  else if (keyword == "HWE_P_CONTROL")
			    {
			      hweCo = atof(str3);
			      if (strcmp (str3, "//") == 0)
				{
				  hweCo = 0.00;
				}
			      logfile << "HWE_P_CONTROL\t" << hweCo << "\n";
			    }
			  else if (keyword == "MAF")
			    {
			      maf = atof(str3);
			      if (strcmp (str3, "//") == 0)
				{
				  maf = 0;
				}
			      logfile << "MAF\t" << maf << "\n";
			    }
			  else if (keyword == "ANNOTATE")
			    {
			      annotate = atoi(str3);
			      if (strcmp (str3, "//") == 0)
				{
				  annotate = 0;
				}

			      logfile << "ANNOTATE\t" << annotate << "\n";
			    }
			  else if (keyword == "GENECOL")
			    {
			      genecol = atoi(str3);

			      if (strcmp (str3, "//") == 0)
				{
				  genecol = 5;
				}

			      if(genecol<3)
				{
				  errorfile << "GENECOL: selection not valid. Must be <=3.\n";
				  logfile << "GENECOL: selection not valid. Must be <=3.\n";
				  cout << "GENECOL: selection not valid. Must be <=3.\n";
				  errorfile.close();logfile.close();exit(1);
				}

			      logfile << "GENECOL\t" << genecol << "\n";
			    }
			  else if (keyword == "MRDIFF")
			    {
			      mrdiff = atof(str3);
			      if (strcmp (str3, "//") == 0)
				{
				  mrdiff = 1;
				}

			      logfile << "MRDIFF\t" << mrdiff << "\n";
			    }
			  else if (keyword == "SINGLE_MARKER" || keyword == "SINGLEMARKER")
			    {
			      singleMarkerTest = atoi(str3);
			      logfile << "SINGLE_MARKER\t" << singleMarkerTest << "\n";
			      if(singleMarkerTest>2){regression=1;}
			      if(singleMarkerTest==2 || singleMarkerTest==4){dfSingle=2;}
			      if (singleMarkerTest>7 || singleMarkerTest <0)
				{
				  errorfile << "SINGLE_MARKER: selection not valid.\n";
				  logfile << "SINGLE_MARKER: selection not valid.\n";
				  cout << "SINGLE_MARKER: selection not valid.\n";
				  errorfile.close();
				  logfile.close();
				  exit(1);
				}
			    }
			  else if (keyword == "PRETEST")
			    {
			      pretest = atoi(str3);
			      logfile << "PRETEST\t" << pretest << "\n";

			      if (pretest !=1 && pretest !=0)
				{
				  errorfile << "PRETEST allowed values: {0,1}\n";
				  logfile << "PRETEST allowed values: {0,1}\n";
				  cout << "PRETEST allowed values: {0,1}\n";
				  errorfile.close();logfile.close();exit(1);
				}
			    }
			  else if (keyword == "PRETEST_CUTOFF")
			    {
			      pHelpCut = atof(str3);
			      logfile << "PRETEST_CUTOFF\t" << pHelpCut << "\n";

			      if (pHelpCut >1 || pHelpCut <0)
				{

				  logfile << "PRETEST_CUTOFF changed to 0.01.\n";
				  cout << "PRETEST_CUTOFF changed to 0.01.\n";
				  logfile.close();
				}
			    }
			  else if (keyword == "TWO_MARKER")
			    {
			      markercombi2 = atoi(str3);
			      logfile << "TWO_MARKER\t" << markercombi2 << "\n";

			      if (markercombi2 > 1)
				{
				  errorfile << "For parameter TWO_MARKER only 0 or 1 is allowed.\n";
				  logfile << "For parameter TWO_MARKER only 0 or 1 is allowed.\n";
				  cout << "For parameter TWO_MARKER only 0 or 1 is allowed.\n";
				  errorfile.close();
				  logfile.close();
				  exit(1);
				}
			    }
			  else if (keyword == "THREE_MARKER")
			    {
			      markercombi3 = atoi(str3);
			      logfile << "THREE_MARKER\t" << markercombi3 << "\n";

			      if (markercombi3 > 1)
				{
				  errorfile << "For parameter THREE_MARKER only 0 or 1 is allowed.\n";
				  logfile << "For parameter THREE_MARKER only 0 or 1 is allowed.\n";
				  cout << "For parameter THREE_MARKER only 0 or 1 is allowed.\n";
				  errorfile.close();
				  logfile.close();
				  exit(1);
				}
			    }
			  else if (keyword == "TEST")
			    {
			      logfile << "TEST\t" << atoi(str3) << "\n";

			      if (strcmp(str3,"M") == 0)
				{
				  readmodel = 1;test=77;
				}
			      else
				{
				  test = atoi(str3);
				}
			      if (test > 2 && test!=15 && test!=16 && test!=19)  // regression
				{
				  regression = 1;
				  if(test==17 || test==18){qt=1;}
				  if(test ==18)
				    {
				      errorfile << "Test 18 is under construction.\n";
				      logfile << "Test 18 is under construction.\n";
				      cout << "Test 18 is under construction.\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				}
			    }
			  else if (keyword == "COVARIATES" && str3!=NULL && str3[0]!='/' )
			    {
			      logfile << "COVARIATES\t" << str3 << "\n";
			      str4 = strcat(str3,";");

			      covar = split(str4, ";");

			      if (str4.find(';') != string::npos)
				{
				  covariate = 1;
				}
			      else
				{
				  covariate = 0;
				}


			      for (int i=0; i<covar.size(); i++)
				{
				  selectedCovar = split(covar[i], "-");

				  if (selectedCovar.size() == 1)
				    {
				      if (atoi(covar[i].c_str()) > maxIndexCov)
					{
					  maxIndexCov = atoi(covar[i].c_str());
					}
				    }
				  else if (selectedCovar.size() == 2)
				    {
				      if (atoi(selectedCovar[1].c_str()) > maxIndexCov)
					{
					  maxIndexCov = atoi(selectedCovar[1].c_str());
					}
				    }
				}

			      if (maxIndexCov <= 0) die("Invalid number of covariates. Non-positive count is not allowed.");
			      cov = (int *) calloc(maxIndexCov, sizeof(int));
			      for (int y = 0; y < covar.size(); y++) {
				selectedCovar = split(covar[y], "-");
				if (selectedCovar.size() == 1) {
				  if (atoi(covar[y].c_str()) <= 0) die("Invalid covariate selected. Non-positive values is not allowed.");
				  cov[atoi(covar[y].c_str())-1] = 1;
				} else if (selectedCovar.size() == 2) {
				  beginCovar = atoi(selectedCovar[0].c_str());
				  endCovar   = atoi(selectedCovar[1].c_str());
				  if (beginCovar <= 0) die("Invalid covariate selected. Non-positive values is not allowed.");
				  for (int kkk = beginCovar; kkk < endCovar+1; kkk++)cov[kkk-1] = 1;
				} else die("Odd things happenings in definition of covariate list.");
			      }
			    }
			  else if (keyword == "SEXCOV")
			    {
			      //logfile.close();exit(1);
			      sexcov= atoi(str3);
			      logfile << "SEXCOV\t" << sexcov << "\n";
			      //printf("sexcov %d\n",sexcov);
			      if (sexcov > 1)
				{
				  errorfile << "For parameter SEXCOV only 0 or 1 is allowed.\n";
				  logfile << "For parameter SEXCOV only 0 or 1 is allowed.\n";
				  cout << "For parameter SEXCOV only 0 or 1 is allowed.\n";
				  errorfile.close();logfile.close();exit(1);
				}

			    }
			  else if (keyword == "PRINT_BETA") //TIM_NEW
			    {
			      choosePrintBeta = atoi(str3);
			      logfile << "PRINT_BETA\t" << choosePrintBeta << "\n";
			    }
			  else if (keyword == "COVARIANCE_MATRIX" || keyword == "COVARIANCEMATRIX")
			    {
			      covariancematrix = atoi(str3);
			      logfile << "COVARIANCE_MATRIX\t" << covariancematrix << "\n";
			    }
			  else if (keyword == "SINGLETOP")
			    {
			      singletop = atoi(str3);
			      logfile << "SINGLETOP\t" << singletop << "\n";
			      //printf("singletop %d\n",singletop);exit(1);
			    }
			  else if (keyword == "M_WITH_SINGLETOP")
			    {
			      mWithSingletop = atoi(str3);
			      logfile << "M_WITH_SINGLETOP\t" << mWithSingletop << "\n";
			      if (markercombi2 == 1 && mWithSingletop > 2)
				{
				  errorfile << "TWO_MARKER: Only M_WITH_SINGLETOP 0,1 or 2 is allowed.\n";
				  logfile << "TWO_MARKER: Only M_WITH_SINGLETOP 0,1 or 2 is allowed.\n";
				  cout << "TWO_MARKER: Only M_WITH_SINGLETOP 0,1 or 2 is allowed.\n";
				  errorfile.close();
				  logfile.close();
				  exit(1);
				}
			      else if (markercombi3 == 1 && mWithSingletop > 3)
				{
				  errorfile << "THREE_MARKER: Only M_WITH_SINGLETOP 0,1,2 or 3 is allowed.\n";
				  logfile << "THREE_MARKER: Only M_WITH_SINGLETOP 0,1,2 or 3 is allowed.\n";
				  cout << "THREE_MARKER: Only M_WITH_SINGLETOP 0,1,2 or 3 is allowed.\n";
				  errorfile.close();
				  logfile.close();
				  exit(1);
				}
			    }
			  else if (keyword == "GENETIC_IMPACT")
			    {
			      geneticImpact = atoi(str3);
			      logfile << "GENETIC_IMPACT\t" << geneticImpact << "\n";
			    }
			  else if (keyword == "M_WITH_GENETIC_IMPACT")
			    {
			      mWithGenetictop = atoi(str3);
			      logfile << "M_WITH_GENETIC_IMPACT\t" << mWithGenetictop << "\n";

			      if (mWithGenetictop > 3)
				{
				  errorfile << "For parameter M_WITH_GENETIC_IMPACT only 0,1,2 or 3 is allowed.\n";
				  logfile << "For parameter M_WITH_GENETIC_IMPACT only 0,1,2 or 3 is allowed.\n";
				  cout << "For parameter M_WITH_GENETIC_IMPACT only 0,1,2 or 3 is allowed.\n";
				  errorfile.close();logfile.close();exit(1);
				}
			    }
			  else if (keyword == "SNP1" && str3[0]!='/')
			    {
			      snp1=(char *) realloc(snp1,(strlen(str3)+1)*sizeof(char));
			      strcpy(snp1,str3);
			      logfile << "SNP1\t" << str3 << "\n";
			      nsnps++;snp1_Selected=1;
			    }
			  else if (keyword == "SNP2" && str3[0]!='/')
			    {
			      snp2=(char *) realloc(snp2,(strlen(str3)+1)*sizeof(char));
			      strcpy(snp2,str3);
			      logfile << "SNP2\t" << str3 << "\n";
			      nsnps++;snp2_Selected=1;
			    }
			  else if (keyword == "SNP3" && str3[0]!='/')
			    {
			      snp3=(char *) realloc(snp3,(strlen(str3)+1)*sizeof(char));
			      strcpy(snp3,str3);
			      logfile << "SNP3\t" << str3 << "\n";
			      nsnps++;snp3_Selected=1;
			    }
			  else if (keyword == "PATHWAY")
			    {
			      pathwayImpact = atoi(str3);
			      logfile << "PATHWAY\t" << pathwayImpact << "\n";

			      if (pathwayImpact > 1)
				{
				  errorfile << "For parameter PATHWAY only 0 or 1 is allowed.\n";
				  logfile << "For parameter PATHWAY only 0 or 1 is allowed.\n";
				  cout << "For parameter PATHWAY only 0 or 1 is allowed.\n";
				  errorfile.close();logfile.close();exit(1);
				}
			    }
			  else if (keyword == "SIMULATION" || keyword == "SIMULATIONS" || keyword == "PERMUTATIONS" || keyword == "PERMUTATION")
			    {
			      nsim = atoi(str3);
			      logfile << "SIMULATION\t" << nsim << "\n";
			      if(nsim>0 && plla2)
				{
				  die("Simulations are not compatible with PARALLELA. Use PARALLELN for this analysis instead.");
				}
			    }
			  else if (keyword == "MC_WITH_SM" || keyword == "SM_WITH_MC")
			    {
			      MCWithSM = atoi(str3);
			      logfile << "MC_With_SM\t" << MCWithSM << "\n";
			      if (MCWithSM > 1)
				{
				  errorfile << "For parameter MC_WITH_SM only 0 or 1 is allowed.\n";
				  logfile << "For parameter MC_WITH_SM only 0 or 1 is allowed.\n";
				  cout << "For parameter MC_WITH_SM only 0 or 1 is allowed.\n";
				  errorfile.close();
				  logfile.close();
				  exit(1);
				}
			    }
			  else if (keyword == "LD_DISTANCE" || keyword == "LD_DIST")
			    {
			      LDdistance = atof(str3);
			      if(LDdistance<=0)
				{
				  errorfile << "LD_DISTANCE must be greater than zero or remove keyword if not needed.\n";
				  logfile << "LD_DISTANCE must be greater than zero or remove keyword if not needed.\n";
				  cout << "LD_DISTANCE must be greater than zero or remove keyword if not needed.\n";
				  errorfile.close();logfile.close();exit(1);
				}
			      logfile << "LD_DISTANCE\t" << LDdistance << "\n";
			    }
			  else if (keyword == "PLUSSINGLE" || keyword == "PLUS_SINGLE")
			    {


			      str4 = str3;

			      plussingleSNPs = split(str4, ";");

			      if (str4.find(';') != string::npos)
				{
				  plusSingle = 1;
				}
			      else
				{
				  plusSingle = 0;
				}
			      logfile << "PLUSSINGLE\t" << plusSingle << "\n";

			      for (int a = 0; a < plussingleSNPs.size(); a++)
				{
				  // Teilstring nach "-" trennen
				  selectedPlussingle = split(plussingleSNPs[a], "-");

				  if (selectedPlussingle.size() == 1)
				    {
				      if (plussingleSNPs[a] == "1")
					{
					  plusSingleFirst = 1;
					}
				      else if (plussingleSNPs[a] == "2")
					{

					  plusSingleSecond = 1;
					}
				      else if (plussingleSNPs[a] == "3")
					{
					  plusSingleThird = 1;
					}

				      else if( atoi(plussingleSNPs[a].c_str())>3)
					{
					  errorfile << "Invalid number of SNPs selected. Not more than 3 SNPs allowed.\n";
					  logfile << "Invalid number of SNPs selected. Not more than 3 SNPs allowed.\n";
					  cout << "Invalid number of SNPs selected. Not more than 3 SNPs allowed.\n";
					  errorfile.close();
					  logfile.close();
					  exit(1);
					}
				    }
				  else if (selectedPlussingle.size() > 1)
				    {
				      beginPlusSingle = atoi(selectedPlussingle[0].c_str());
				      endPlusSingle = atoi(selectedPlussingle[1].c_str());

				      if( endPlusSingle > 3)
					{
					  errorfile << "Invalid number of SNPs selected. Not more than 3 SNPs allowed.\n";
					  logfile << "Invalid number of SNPs selected. Not more than 3 SNPs allowed.\n";
					  cout << "Invalid number of SNPs selected. Not more than 3 SNPs allowed.\n";
					  errorfile.close();
					  logfile.close();
					  exit(1);
					}

				      for (int b = beginPlusSingle; b < endPlusSingle+1; b++)
					{
					  if (b == 1)
					    {
					      plusSingleFirst = 1;
					    }
					  else if (b == 2)
					    {

					      plusSingleSecond = 1;
					    }
					  else if (b == 3)
					    {
					      plusSingleThird = 1;
					    }
					}
				    }
				}
			    }
			  else if (keyword == "ALLALLMATCHING") {
			      allAllMatching = (bool)atoi(str3);
			      logfile << keyword << "\t" << allAllMatching << "\n";
			    }
			  else if (keyword == "STRATIFY" || keyword == "DESTRATIFY") {
			      stratify = atoi(str3);
			      logfile << keyword << "\t" << stratify << "\n";
			    }
			  else if (keyword == "STRATIFY_STRUCT") {
			      if (strcmp(str3, "PAIR") && strcmp(str3, "GROUP") && strcmp(str3, "GROUP2") && strcmp(str3, "CLUSTER")) { sstm << "Do not unterstand STRATIFY_STRUCT " << str3 << "."; logg(sstm); }
			      else {
				stratify_struct = str3;
				logfile << keyword << "\t" << stratify_struct << "\n";
			      }
			    }
			  else if (keyword == "STRATIFY_VALID") {
			      if (strcmp(str3, "IBS") && strcmp(str3, "VICINITY") && strcmp(str3, "CLUSTER")) { sstm << "Do not unterstand STRATIFY_VALID " << str3 << "."; logg(sstm); }
			      else {
				    stratify_valid = str3;
				    logfile << keyword << "\t" << stratify_valid << "\n";
			      }
			  }
			  else if (keyword == "T_ALL") {
			      t_all = atoi(str3);
			      logfile << keyword << "\t" << t_all << "\n";
            }
			  else if (keyword == "T_CC") {
			      t_cc = atoi(str3);
			      logfile << keyword << "\t" << t_cc << "\n";
			  }
			  else if (keyword == "CLUSTER_COVARS") {
			      cluster_covars = (bool)atoi(str3);
			      logfile << keyword << "\t" << cluster_covars << "\n";
			  }
			  else if (keyword == "GROUP_TEST") {
			      group_test = 1;
			      if (!strcmp(str3, "LINEAR")) group_test = 1;
			      if (!strcmp(str3, "SQUARE")) group_test = 2;
			      logfile << keyword << "\t" << group_test << "\n";
			  }
			  else if (keyword == "FAMCLUSTER") {
			      familyCluster = atoi(str3);
			      logfile << keyword << "\t" << familyCluster << "\n";
			  }
			  else if (keyword == "SIM_DELTA") {
			      sim_limit = pow(10,atoi(str3));
			      logfile << keyword << "\t" << atoi(str3) << "\n";
              }
              else if (keyword == "SKIP_IT_QC") {
			      skip_it_qc = (bool)atoi(str3);
			      logfile << keyword << "\t" << skip_it_qc << "\n";
              }
			  else if (keyword == "PRINTTOP") {
			      printtop = atoi(str3);
			      logfile << "PRINTTOP\t" << printtop << "\n";
			      if(printtop<=0){printtop=10;}
			      if(printtop>=100000)
				{
				  cout << "Warning! Large PRINTTOP may require too much running time and computer memory.\n";
				  logfile << "Warning! Large PRINTTOP may require too much running time and computer memory.\n";

				}
			    }
			  else if (keyword == "PATHWAYANALYSIS")
			    {
			      pathwayAnalysis = atoi(str3);
			      if(pathwayAnalysis>1)
				{
				  pathwayAnalysis=1;
				}

			      logfile << "PATHWAYANALYSIS\t" << pathwayAnalysis << "\n";
			    }
			  else if (keyword == "PATHWAYTEST")
			    {
			      pathwayTest = atoi(str3);
			      if (pathwayTest!=1 && pathwayTest!=2 && pathwayTest!=3 && pathwayTest!=4 && pathwayTest!=5 && pathwayTest!=6)
				{
				  pathwayTest=2;
				}
			      logfile << "PATHWAYTEST\t" << pathwayTest << "\n";
			    }
			  else if (keyword == "P_INTERRATIO")
			    {
			      p_interratio = atof(str3);
			      if((strcmp (str3, "//") == 0))
				{
				  p_interratio=0.05;
				}
			      logfile << "P_INTERRATIO\t" << p_interratio << "\n";
			    }
			  else if (keyword == "P_PAA")
			    {
			      p_PAA = atof(str3);
			      if((strcmp (str3, "//") == 0))
				{
				  p_PAA=0.05;
				}
			      logfile << "P_PAA\t" << p_PAA << "\n";
			    }
			  else if (keyword == "QT")
			    {
			      qt = atoi(str3);

			      if(qt!=0){qt=1;}

			      logfile << "QT\t" << qt << "\n";
			    }
			  else if (keyword == "MISSING_PHENO")
			    {
			      missingpheno = atof(str3);
			      logfile << "MISSING_PHENO\t" << missingpheno << "\n";
			    }
			  else if (keyword == "HAPLO")
			    {
			      haplo = atoi(str3);
			      if(haplo>=1)
				{
				  haplo=1;
				  regression=1;
				  qt=1;
				}
			      else
				{
				  haplo=0;
				}
			      logfile << "HAPLO\t" << haplo << "\n";
			      cout << "Haplotype test follows haplotype trend regression (HTR) suggested by Zaykin et. al, Hum Hered 53:79-91, will be performed. Internally handled with QT set to 1." << "\n";
			      logfile << "Haplotype test follows haplotype trend regression (HTR) suggested by Zaykin et. al, Hum Hered 53:79-91, will be performed. Internally handled with QT set to 1." << "\n";

			    }
			  else if (keyword == "HAPLO_DIST" || keyword == "HAPLO_DISTANCE")
			    {
			      distance = atof(str3);
			      if(distance<=0)
				{
				  errorfile << "HAPLO_DIST must be greater than zero or remove keyword if not needed.\n";
				  logfile << "HAPLO_DIST must be greater than zero or remove keyword if not needed.\n";
				  cout << "HAPLO_DIST must be greater than zero or remove keyword if not needed.\n";
				  errorfile.close();logfile.close();exit(1);
				}
			      logfile << "HAPLO_DIST\t" << distance << "\n";
			    }
			  else if (keyword == "DOHAPFILE")
			    {
			      dohapfile = atoi(str3);
			      logfile << "DOHAPFILE\t" << dohapfile << "\n";
			    }
			  else if (keyword == "DOIBS")
			    {
			      doIbs = atoi(str3);
			      logfile << "DOIBS\t" << doIbs << "\n";
			    }
			  else if (keyword == "IBS_MINOR_ONLY")
			    {
			      ibs_minor_only = (bool)atoi(str3);
			      ibs_SD_Selected=1;
			      if (ibs_minor_only) {
				getIbs       = getMinorIbs;
				getIbsMatrix = getMinorIbsMatrix;
			      }
			      logfile << "IBS_MINOR_ONLY\t" << ibs_minor_only << "\n";
			    }
			  else if (keyword == "IBS_SD_RELATIVES")
			    {
			      ibs_SD_relatives = atof(str3);
			      ibs_SD_Selected=1;
			      logfile << "IBS_SD_RELATIVES\t" << ibs_SD_relatives << "\n";
			    }
			  else if (keyword == "IBS_SD_OUTLIER")
			    {
			      ibs_SD_outlier = atof(str3);
			      ibs_SD_Selected=1;
			      logfile << "IBS_SD_OUTLIER\t" << ibs_SD_outlier << "\n";
			    }
			  else if (keyword == "PFILTER")
			    {
			      pfilter = atof(str3);
			      logfile << "PFILTER\t" << pfilter << "\n";
			    }
			  else if (keyword == "FAMDATA")
			    {
			      familyData = atoi(str3);
			      logfile << "FAMDATA\t" << familyData << "\n";
			    }
			  else if (keyword == "SWITCH_CASE_CONTROL")
			    {
			      switchCaseControl = atoi(str3);
			      logfile << "SWITCH_CASE_CONTROL\t" << switchCaseControl << "\n";
			    }
			  else if (keyword == "SNP_COVARIATES")
			    {
			      logfile << "SNP_COVARIATES\t" << str3 << "\n";

			      str4 = str3;
			      snpCov = 1;

			      SNPcovar = split(str4, ";");
			      numberOfSNPCov = SNPcovar.size();
			      cout << "SNPcovariates: \n";
			      logfile << "SNPcovariates: \n";
			      for (int y = 0; y < numberOfSNPCov; y++)
				{
				  cout << "SNPcovar: " << SNPcovar[y] << "\n";
				  logfile << "SNPcovar: " << SNPcovar[y] << "\n";
				}
			      cout << "Note: snp-covariates are used for regression tests even when they are qc-out. Only exception: liability analysis.\n";
			      logfile << "Note: snp-covariates are used for regression tests even when they are qc-out. Only exception: liability analysis.\n";
			    }
			  else if (keyword == "SNP_RISK")
			    {
			      logfile << "SNP_RISK\t" << str3 << "\n";
			      riskSNP = 1;

			      str4 = str3;

			      SNPrisk = split(str4, ";");
			      cout << "Risk alleles: \n";
			      logfile << "Risk alleles: \n";
			      for (int y = 0; y < SNPrisk.size(); y++)
				{
				  cout << "SNPrisk: " << SNPrisk[y] << "\n";
				  logfile << "SNPrisk: " << SNPrisk[y] << "\n";
				}
			    }
			  else if (keyword == "LIABILITY")
			    {
			      logfile << "LIABILITY\t" << atoi(str3) << "\n";

			      liability = atoi(str3);
			    }
			  else if(keyword == "LAMBDA_ADJUST_PW")
			    {
			      logfile << "LAMBDA_ADJUST_PW\t" << atoi(str3) << "\n";
			      lambdaAdjust = atoi(str3);
			    }
			  else if(keyword == "VERBOSE" || keyword == "VERBOSE_OUTPUT")
			    {
			      logfile << "VERBOSE_OUTPUT\t" << atoi(str3) << "\n";
			      verbose = atoi(str3);
			      if(verbose!=0 && verbose!=1 && verbose!=2){
				cout <<"VERBOSE_OUTPUT takes only values 0, 1 and 2!"<<endl;
				logfile <<"VERBOSE_OUTPUT takes only values 0, 1 and 2!"<<endl;
				errorfile <<"VERBOSE_OUTPUT takes only values 0, 1 and 2!"<<endl;
				verbose=1;
			      }
			    }
			  else if(keyword == "PROB_REQUIRED")
			    {
			      logfile << "PROB_REQUIRED\t" << atof(str3) << "\n";
			      probRequired = atof(str3);
			      if(probRequired>1 || probRequired<0)
				{
				  cout <<"PROB_REQUIRED takes only values between 0 and 1!"<<endl;
				  logfile <<"PROB_REQUIRED takes only values between 0 and 1!"<<endl;
				  errorfile <<"PROB_REQUIRED takes only values between 0 and 1!"<<endl;
				  probRequired=0;
				}
			    }
			  else if (keyword == "DOSAGE")
			    {
			      dosage = atoi(str3);
			      if (dosage > 0) dosage=1;
			      if(dosage==1 && filemode!=4)
				{
				  errorfile << "It isn't possible to choose DOSAGE without IMPUTE2 or IFILE.\n";
				  logfile << "It isn't possible to choose DOSAGE without IMPUTE2 or IFILE.\n";
				  cout << "It isn't possible to choose DOSAGE without IMPUTE2 or IFILE.\n";
				  errorfile.close();logfile.close();exit(1);
				}
			      logfile << "DOSAGE\t" << dosage << "\n";
			    }
#if RARE  // intersnpRARE
			  else if(keyword == "BINSIZE_SNP" || keyword == "BINSIZE_DIST") // number of SNPs
			    {
			      cout<<keyword<<" is deprecated. Use BINSIZE_RARE or INTERVALFILE."<<endl;
			      logfile<<keyword<<" is deprecated. Use BINSIZE_RARE or INTERVALFILE."<<endl;
			      errorfile<<keyword<<" is deprecated. Use BINSIZE_RARE or INTERVALFILE."<<endl;
			      exit(1);
			    }
			  else if (keyword == "BINSIZE_RARE") // number of SNPs
			    {
			      binsizeRare = atoi(str3);
			      logfile << "BINSIZE_RARE\t" << binsizeRare << "\n";
			    }
			  else if (keyword == "RARE" || keyword=="MAFT")
			    {
			      raref = (double)atof(str3);
			      if(raref<=EPS)
				{
				  cout << "MAFT too small!"<<endl;
				  logfile << "MAFT too small!"<<endl;
				  errorfile << "MAFT too small!"<<endl;
				  exit(1);
				}
			      if((raref-0.5)>=EPS)
				{
				  cout << "MAFT! MAF is always between 0 and 0.5"<<endl;
				  cout << "Setting MAFT to 0.5"<<endl;
				  logfile << "MAFT too large! MAF is always between 0 and 0.5"<<endl;
				  logfile << "Setting MAFT to 0.5"<<endl;
				  raref=0.5;
				}

			      logfile << "MAFT\t" << raref <<"\n";
			      bin=1;
			    }
			  else if (keyword == "NCT")
			    {
			      NCT = (double)atof(str3);
			      logfile << "NCT\t" << NCT <<"\n";
			      bin=1;
			    }
			  else if(keyword == "OPTIMAL" || keyword == "VT")
			    {
			      if(isdigit(str3[0]))
				{
				  optimalrare = atoi(str3);
				  logfile <<"VT\t"<<optimalrare<<endl;
				  if(optimalrare!=0 && optimalrare!=1)
				    {
				      cout << "VT can be either 0 or 1\n";
				      logfile << "VT can be either 0 or 1\n";
				      errorfile << "VT can be either 0 or 1\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				}
			    }
			  else if(keyword == "VB")
			    {
			      if(isdigit(str3[0]))
				{
				  vb = atoi(str3);
				  logfile <<"VB\t"<<vb<<endl;
				  if(vb!=0 && vb!=1)
				    {
				      cout << "VB can be either 0 or 1\n";
				      logfile << "VB can be either 0 or 1\n";
				      errorfile << "VB can be either 0 or 1\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				  if(vb==1){
				    bin=1;
				  }
				}
			    }
			  else if(keyword == "VB_PRINT_PERM_STAT"  || keyword == "VB_PRINT_PERM_P")
			    {
			      if(isdigit(str3[0]))
				{
				  vb_print_perm_stat = atoi(str3);
				  logfile <<"VB_PRINT_PERM_STAT\t"<<vb_print_perm_stat<<endl;
				  if(vb_print_perm_stat!=0 && vb_print_perm_stat!=1)
				    {
				      cout << "VB_PRINT_PERM_STAT can be either 0 or 1\n";
				      logfile << "VB_PRINT_PERM_STAT can be either 0 or 1\n";
				      errorfile << "VB_PRINT_PERM_STAT can be either 0 or 1\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				}
			    }

			  else if(keyword == "VB_BINWISE_CORR" || keyword == "VB_BINSIZE_CORR")
			    {
			      if(isdigit(str3[0]))
				{
				  vb_binwise_corr = atoi(str3);
				  logfile <<"VB_BINWISE_CORR\t"<<vb_binwise_corr<<endl;
				  if(vb_binwise_corr!=0 && vb_binwise_corr!=1)
				    {
				      cout << "VB_BINWISE_CORR can be either 0 or 1\n";
				      logfile << "VB_BINWISE_CORR can be either 0 or 1\n";
				      errorfile << "VB_BINWISE_CORR can be either 0 or 1\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				}
			    }
			  else if(keyword == "VB_BMP")
			    {
			      if(isdigit(str3[0]))
				{
				  vb_bmp = atoi(str3);
				  logfile <<"VB_BMP\t"<<vb_bmp<<endl;
				  if(vb_bmp!=0 && vb_bmp!=1)
				    {
				      cout << "VB_BMP can be either 0 or 1\n";
				      logfile << "VB_BMP can be either 0 or 1\n";
				      errorfile << "VB_BMP can be either 0 or 1\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				}
			    }
			  else if (keyword == "RARE_STRATIFY") {
			      rare_stratify = (bool)atoi(str3);
			      logfile << keyword << "\t" << rare_stratify << "\n";
			    }
			  else if (keyword == "RARE_STRAT_MIN_SNPS") {
			      rare_strat_min = atoi(str3);
			      logfile << keyword << "\t" << rare_strat_min << "\n";
			    }
			  else if (keyword == "RARE_STRAT_MAFT") {
			      rare_strat_maft = atof(str3);
			      if(rare_strat_maft < 0.0 || 0.5 < rare_strat_maft) die("RARE_STRAT_MAFT has to be in [0;0.5]!");
			      logfile << keyword << "\t" << rare_strat_maft << "\n";
			    }
			  else if(keyword == "WEIGHTS")
			    {
			      logfile << "WEIGHTS\t" << str3 << "\n";
			      if (str3 != NULL && str3[0]!='/')
				{
				  str4=str3;
				  weightParams=split(str4, ";");
				  if(weightParams.size()==1 && (weightParams[0]=="MAF" || weightParams[0]=="SD" ||  weightParams[0]=="1/SD"))
				    {
				      weights.mode=1;
				      cout << "MAF-weighted variants (Madsen and Browning, 2009)"<<endl;
				      logfile << "MAF-weighted variants (Madsen and Browning, 2009)"<<endl;
				    }
				  if(weightParams.size()==3 && weightParams[0]=="BETA")
				    {
				      weights.mode=2;
				      weights.betapar1=atoi(weightParams[1].c_str());
				      weights.betapar2=atoi(weightParams[2].c_str());
				      cout << "Beta distribution-weighted variants with parameters "<<weights.betapar1<<", "<<weights.betapar2<<""<<endl;
				      logfile << "Beta distribution-weighted variants with parameters "<<weights.betapar1<<", "<<weights.betapar2<<" "<<endl;
				    }

				  if(weightParams.size()==3 && weightParams[0]=="LOGISTIC")
				    {
				      weights.mode=3;
				      weights.logpar1=atof(weightParams[1].c_str());
				      weights.logpar2=atof(weightParams[2].c_str());
				      cout << "Logfistic weights with parameters "<<weights.logpar1<<", "<<weights.logpar2<<""<<endl;
				      logfile << "Logistic weights with parameters "<<weights.logpar1<<", "<<weights.logpar2<<""<<endl;
				    }// Eingabe nach Simikolon trennen
				}
			    }
			  else if(keyword == "RARE_TESTS" || keyword == "RARE_TEST")
			    {
			      logfile << "RARE_TESTS\t" << str3 << "\n";
			      if (str3 != NULL && str3[0]!='/')
				{
				  // Eingabe nach Simikolon trennen
				  str4 = str3;
				  selectedRareTests = split(str4, ";");


				  for (y = 0; y < selectedRareTests.size(); y++)
				    {

				      if(selectedRareTests[y]=="FR" || selectedRareTests[y]=="FISHER-RARE")
					{
					  FISHERtest=1;
					}
				      if(selectedRareTests[y]=="COLL")
					{
					  COLLtest=1;
					}
				      if(selectedRareTests[y]=="CMAT")
					{
					  CMATtest=1;
					}
				      if(selectedRareTests[y]=="REG" || selectedRareTests[y]=="REGRESSION")
					{
					  REGRESSIONtest=1;
					  regression=1;
					}
				      if(selectedRareTests[y]=="FRACREG")
					{
					  FRACREGtest=1;
					  regression=1;
					  collapseRare=1;
					}
				      if(selectedRareTests[y]=="COLLREG")
					{
					  COLLREGtest=1;
					  regression=1;
					  collcollapseRare=1;
					}
				    }
				}
			    }
			  else if (keyword =="PRESCREEN" || keyword == "RARE_REG_PRETEST"||keyword == "RARE_PRE_PRETEST")
			    {
			      if(str3!=NULL && str3[0]!='/')
				{

				  if(isdigit(str3[0]))
				    {
				      logfile << "PRESCREEN\t" << str3<<endl;
				      rareregpretest=atof(str3);

				      if(rareregpretest==0)
					{
					  cout   <<"Argument of PRESCREEN is ZERO! No prescreening will be conducted."<<endl;
					  logfile<<"Argument of PRESCREEN is ZERO! No prescreening will be conducted."<<endl;
					}
				      //string str4;

				    }

				}
			    }
			  else if((keyword == "GENERATE_INTERVALS" || keyword == "GENERATE_SETID" ||  keyword == "GENERATE_SETIDS" || keyword == "INTERVAL_EDITOR" || keyword == "INTERVALEDITOR") && str3 != NULL && str3[0]!='/')
			    {
			      logfile << "GENERATE_SETID\t" << str3 << "\n";
			      intervaleditor=atoi(str3);
			      //			      setid=1;
			    }
			  else if (keyword == "ADAPTIVE")
			    {
			      if(str3!=NULL && str3[0]!='/')
				{

				  if(isdigit(str3[0]))
				    {
				      logfile << "ADAPTIVE\t" << str3<<endl;
				      adaptive=atoi(str3);

				      if(adaptive!=0 && adaptive!=1)
					{
					  cout   <<"ADAPTIVE is 0 or 1! Default: 0"<<endl;
					  logfile<<"ADAPTIVE is 0 or 1! Default: 0."<<endl;
					}
				      //string str4;

				    }

				}
			    }
			  else if (keyword == "WILSON_PRETEST")
			    {
			      //				cout << "WILSON_PRETEST is obsolete! Use ADAPTIVE instead"<<endl;
			      //				logfile << "WILSON_PRETEST is obsolete! Use ADAPTIVE instead"<<endl;
			      //				errorfile << "WILSON_PRETEST is obsolete! Use ADAPTIVE instead"<<endl;
			      //				exit(1);
			      if(str3!=NULL && str3[0]!='/')
				{

				  if(isdigit(str3[0]))
				    {
				      logfile << "WILSON_PRETEST\t" << str3<<endl;
				      wilsonpretest=atof(str3);

				      if(wilsonpretest==0)
					{
					  cout   <<"Argument of WILSON_PRETEST is ZERO! No Wilson pretest will be conducted."<<endl;
					  logfile<<"Argument of WILSON_PRETEST is ZERO! No Wilson pretest will be conducted."<<endl;
					}
				      //string str4;

				    }

				}
			    }
			  else if (keyword == "RARE_PRETEST")
			    {
			      cout << "RARE_PRETEST is obsolete! Use ADAPTIVE instead"<<endl;
			      logfile << "RARE_PRETEST is obsolete! Use ADAPTIVE instead"<<endl;
			      errorfile << "RARE_PRETEST is obsolete! Use ADAPTIVE instead"<<endl;
			      exit(1);
			      if(str3!=NULL && str3[0]!='/')
				{

				  if(isdigit(str3[0]))
				    {
				      logfile << "RARE_PRETEST\t" << str3<<endl;
				      str3 = strtok(str3, "; \t");
				      rarepretest=atoi(str3);

				      if(rarepretest==0)
					{
					  cout   <<"RARE_PRETEST is ZERO! No pretest will be conducted."<<endl;
					  logfile<<"RARE_PRETEST is ZERO! No pretest will be conducted."<<endl;
					}
				      cout   <<"RARE_PRETEST is obsolete and not supoted anymore. Use WILSON_PRETEST instead."<<endl;
				      logfile   <<"RARE_PRETEST is obsolete and not supoted anymore. Use WILSON_PRETEST instead."<<endl;

				      //string str4;

				      //string str4;
				      /*
					cout <<str3<<flush<<endl;
					str3 = strtok(NULL, "; \t");
					cout <<str3<<flush<<endl;
					if(!isdigit(str3[0]))
					{
					cout <<"Can not interpret p-limit for pretest! Modify your selectionfile accordingly."<<endl;
					errorfile <<"Can not interpret p-limit for pretest! Modify your selectionfile accordingly."<<endl;    logfile <<"Can not interpret p-limit for pretest! Modify your selectionfile accordingly."<<endl;
					logfile.close();errorfile.close();exit(1);
					}
					else if((1/(float)rarepretest-atof(str3))>EPS)
					{
					cout <<"It is impossible not to pass RARE_PRETEST! With " << rarepretest << " simulations, we need AT LEAST p-value limit of " << 1/(float)rarepretest << " for RARE_PRETEST to be effective! Modify selectionfile accordingly."<<endl;
					errorfile <<"It is impossible not to pass RARE_PRETEST! With " << rarepretest << " simulations, we need AT LEAST p-value limit of " << 1/(float)rarepretest << " for RARE_PRETEST to be effective! Modify selectionfile accordingly."<<endl;
					logfile <<"It is impossible not to pass RARE_PRETEST! With " << rarepretest << " simulations, we need AT LEAST p-value limit of " << 1/(float)rarepretest << " for RARE_PRETEST to be effective! Modify selectionfile accordingly."<<endl;
					logfile.close();errorfile.close();exit(1);exit(1);
					}
					else
					{
					rarepretestlimit=atof(str3);
					if(rarepretestlimit<=0 || rarepretestlimit>=1)
					{
					cout << "p-Limit for pretest is " << rarepretestlimit <<", but it has to be a number between 0 and 1!"<<endl;
					logfile << "p-Limit for pretest is " << rarepretestlimit <<", but it has to be a number between 0 and 1!"<<endl;
					errorfile << "p-Limit for pretest is " << rarepretestlimit <<", but it has to be a number between 0 and 1!"<<endl;
					logfile.close();errorfile.close();exit(1);

					}

					}
					//continue;
					*/
				    }

				}
			    }
			  else if(keyword == "LAMBDA_ADJUST" && str3 != NULL && str3[0]!='/')
			    {
			      logfile << "LAMBDA_ADJUST\t" << str3 << "\n";
			      if (str3[0]=='1')
				{
				  fisherCorrection=1;
				  cout <<"Fisher lambda-correction was selected. Be careful not to apply to small datasets." <<endl;
				  logfile <<"Fisher lambda-correction was selected. Be careful not to apply to small datasets." <<endl;
				}
			    }
			  else if(keyword == "MAF_ADJUST" && str3 != NULL && str3[0]!='/')
			    {
			      logfile << "MAF_ADJUST\t" << str3 << "\n";
			      if (str3[0]=='0')
				{
				  mafadjust=0;
				}
			    }
			  else if((keyword == "BIN_EXPLORE" || keyword == "DENSE_BINNING") && str3 != NULL && str3[0]!='/')
			    {
			      logfile << "DENSE_BINNING\t" << str3 << "\n";
			      cout << "DENSE_BINNING\t" << str3 << "\n";
			      char *str4=NULL;

			      str4=strtok(str3,";");

			      cout <<str4<<endl;

			      string tempchr=" ";
			      cout<<str4<<" "<<str4[1] <<endl;
			      tempchr.assign(str4);
			      if(str4[0]=='N')
				{
				  if(tempchr.length()<3&&tempchr.length()!=1)
				    {
				      cout <<"Usage: DENSE_BINNING MIN_RARE MAX_RARE."<<endl;
				      logfile <<"Usage: DENSE_BINNING MIN_RARE MAX_RARE."<<endl;
				      errorfile <<"Usage: DENSE_BINNING MIN_RARE MAX_RARE."<<endl;
				      exit(1);
				    }
				  if(tempchr.length()==1)
				    {
				      binamin=0;
				    }
				  else if(str4[1]=='+'||str4[1]=='-')
				    {
				      tempchr.assign(tempchr.substr(2,tempchr.length()));

				      // Add 10000 to distinguish between fixed constant values
				      binamin=atoi(tempchr.c_str());
				      cout<<str4<<" "<<str4[1] <<endl;
				      if(str4[1]=='+')
					{
					  cout <<"MIN_RARE is too large."<<endl;
					  logfile <<"MIN_RARE is too large."<<endl;
					  errorfile <<"MIN_RARE is too large."<<endl;
					  exit(1);
					}
				      else if(str4[1]=='-')
					{
					  binamin=-binamin+10000;
					}

				    }
				}
			      else
				{
				  binamin=atoi(str4);
				}
			      cout<<binamin<<endl;
			      str4=strtok(NULL,";");
			      cout<<binamin<<endl;
			      if(str4==NULL)
				{
				  cout <<"Usage: DENSE_BINNING MIN_RARE MAX_RARE."<<endl;
				  logfile <<"Usage: DENSE_BINNING MIN_RARE MAX_RARE."<<endl;
				  errorfile <<"Usage: DENSE_BINNING MIN_RARE MAX_RARE."<<endl;
				  exit(1);
				}
			      cout<<str4<<endl;
			      tempchr.assign(str4);
			      if(str4[0]=='N')
				{
				  if(tempchr.length()<3)
				    {
				      binamax=10000;
				    }
				  else
				    {
				      tempchr.assign(tempchr.substr(1,tempchr.length()));
				      // Add 10000 to distinguish between fixed constant values
				      binamax=atoi(tempchr.c_str());
				      binamax=binamax+10000;
				    }
				}
			      else
				{
				  binamax=atoi(str4);
				}

			      cout <<"DENSE_BINNING "<<binamin <<"\t"<<binamax<<"." <<endl;
			      logfile <<"DENSE_BINNING "<<binamin <<"\t"<<binamax<<"." <<endl;

			      if (binamin<2)
				{
				  cout << "MIN_RARE too small. Set to 2."<<endl;
				  logfile << "MIN_RARE too small. Set to 2."<<endl;
				}
			      if (binamax-binamin<2)
				{

				  cout << "MAX_RARE-MIN_RARE="<<binamax-binamin<<"."<<endl;
				  logfile << "MAX_RARE-MIN_RARE="<<binamax-binamin<<"."<<endl;
				  errorfile << "MAX_RARE-MIN_RARE="<<binamax-binamin<<"."<<endl;
				  exit(1);
				}


			    }
			  else if ((keyword == "INTERVALFILE" || keyword == "INTERVALFILE_IN" ||keyword == "INTERVAL_IN") && str3!=NULL && (str3[0]!='/' || str3[1]!='/'))
			    {
			      intervalfile = str3;
			      bin=1;
			      if(!file_exists(intervalfile)){
				cout<<intervalfile<<" Not found!"<<endl;
				logfile<<intervalfile<<" Not found!"<<endl;
				errorfile<<intervalfile<<" Not found!"<<endl;
				exit(1);
			      }
			      cout << "Input intervalfile has been selected!" << endl;
			      binsizeRare=0;
			      logfile << "INTERVAL_IN\t" << intervalfile << "\n";
			    }
			  else if ((keyword == "MERGE_INTERVALS" || keyword == "MERGING" || keyword == "MERGE") && str3!=NULL &&  str3[0]!='/')
			    {
			      if(atoi(str3)==1)
				{
				  merging = 1;
				  cout << "Bins from intervalfile will be merged, if they overlap."<<endl;
				}
			      else if (bin)
				{
				  cout<< "No merging for bins from intervalfile"<<endl;
				}

			      logfile << "MERGE_INTERVALS\t" << str3 << "\n";

			    }
			  else if (keyword == "FLANKING" && str3!=NULL && str3[0]!='/')
                            {

			      if(atoi(str3)>=0)
				{
				  flanking = atoi(str3);
				}
			      else
				{
				  cout << "Can not recognize FLANKING! Must be a positive integer"<<endl;
				  logfile << "Can not recognize FLANKING! Must be a positive integer"<<endl;
				  errorfile << "Can not recognize FLANKING! Must be a positive integer"<<endl;
				  exit(1);
				}
			      logfile << "FLANKING\t"<<str3<<endl;

			    }
			  else if ((keyword == "INTERVALFILE_FORMAT"||keyword == "INTERVAL_FORMAT") && str3!=NULL && str3[0]!='/')
			    {
			      intervalfile_format = str3;
			      cout << "Intervalfile format is encoded in string "<< intervalfile_format<<"."<<endl;

			      logfile << "INTERVAL_FORMAT\t" << intervalfile_format << "\n";

			    }
			  else if ((keyword =="CLOSE_GAPS"|| keyword =="FILL_GAPS"||keyword =="FILLGAPS") && str3!=NULL && str3[0]!='/')
			    {
			      if(atoi(str3)==1)
				{
				  expandIntervals=1;
				  cout << "FILL_GAPS: Intervals will be auto-expanded."<< endl;
				}
			      else if(atoi(str3)!=0)
				{
				  cout << "Can not interpret FILL_GAPS "<< expandIntervals << endl;
				}
			      logfile << "FILL_GAPS\t" << str3 << "\n";

			    }
			  else if ((keyword =="CONCATENATE_INTERVALS" || keyword =="CONCATENATE") && str3!=NULL && str3[0]!='/')
			    {

			      if(atoi(str3)>=0)
				{
				  catIntervals=atoi(str3);
				  cout << "This many consecutive intervals will be concatenated: "<< catIntervals<< endl;
				}
			      else
				{
				  cout << "Can not interpret CONCATENATE: "<< catIntervals << endl;
				}
			      if(catIntervals==1)
				{
				  cout <<"CONCATENENATE 1 has no effect! Set to 0."<<endl;
				  logfile <<"CONCATENENATE 1 has no effect! Set to 0."<<endl;
				  catIntervals=0;
				}
			      logfile << "CONCATENATE\t" << str3 << "\n";

			    }
			  else if ((keyword =="MIN_RARE_IN_BIN"||keyword =="FILTER_SMALL_BINS") && str3!=NULL && str3[0]!='/')
			    {

			      if(atoi(str3))
				{
				  minRareInBin=atoi(str3);
				  cout << "MIN_RARE_IN_BIN: a bin requires at least "<< minRareInBin <<" rare variants."<< endl;
				}
			      else
				{
				  cout << "Can not interpret MIN_RARE_IN_BIN: "<< minRareInBin << endl;
				}
			      logfile << "MIN_RARE_IN_BIN\t" << str3 << "\n";

			    }
			  else if ( (keyword =="SPLIT_LARGE_BINS"||keyword =="MAX_RARE_IN_BIN") && str3!=NULL && str3[0]!='/')
			    {

			      if(atoi(str3))
				{
				  maxRareInBin=atoi(str3);
				  cout << "MAX_RARE_IN_BIN: a bin can have "<< maxRareInBin <<" rare variants at most."<< endl;
				}
			      else
				{
				  cout << "Can not interpret MAX_RARE_IN_BIN: "<< maxRareInBin << endl;
				}
			      logfile << "MAX_RARE_IN_BIN\t" << str3 << "\n";

			    }
			  else if(keyword == "BIN_ADJUST")
			    {
			      if(isdigit(str3[0]))
				{
				  binadjust = atoi(str3);
				  logfile <<"BIN_ADJUST\t"<<binadjust<<endl;
				  if(binadjust!=0 && binadjust!=1)
				    {
				      cout << "BIN_ADJUST can be either 0 or 1\n";
				      logfile << "BIN_ADJUST can be either 0 or 1\n";
				      errorfile << "BIN_ADJUST can be either 0 or 1\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				}
			    }
			  else if(keyword == "COLL_INTER")
			    {
			      if(isdigit(str3[0]))
				{
				  collinter = atoi(str3);
				  logfile <<"COLL_INTER\t"<<collinter<<endl;
				  if(collinter!=0 && collinter!=1 && collinter!=2 && collinter!=3 && collinter!=4)
				    {
				      cout << "COLL_INTER only takes values between 0 and 4!\n";
				      logfile << "COLL_INTER only takes values between 0 and 4!\n";
				      errorfile << "COLL_INTER only takes values between 0 and 4!\n";
				      errorfile.close();logfile.close();exit(1);
				    }

				  else if(collinter>0 && collinter <=2){
				    COLLflag=1;
				    CMATflag=0;
				    FISHERflag=0;
				    REGRESSIONflag=0;
				    FRACREGflag=0;
				    COLLREGflag=0;
				  }
				  else if(collinter==3 || collinter==4){
				    regression=1;
				    COLLflag=0;
				    CMATflag=0;
				    FISHERflag=0;
				    REGRESSIONflag=0;
				    FRACREGflag=0;
				    COLLREGflag=1;
				  }
				}
			    }
#endif // RARE
			  else if (keyword == "OUTPUTNAME" && str3[0]!='/' && str3!=NULL)
			    {
			      outputname = str3;

			      if (outputname == "//")
				{
				  outputname = "";
				}
			      logfile << "OUTPUTNAME\t" << outputname << "\n";
			    }
			  else if (keyword[0] != '#' && strcmp(keyword.c_str(),"FILE")!=0 && strcmp(keyword.c_str(),"TFILE")!=0 && strcmp(keyword.c_str(),"BFILE")!=0 && strcmp(keyword.c_str(),"IFILE")!=0 && strcmp(keyword.c_str(),"FAM")!=0 && strcmp(keyword.c_str(),"TFAM")!=0 && strcmp(keyword.c_str(),"BFAM")!=0 && strcmp(keyword.c_str(),"MAP")!=0 && strcmp(keyword.c_str(),"BMAP")!=0 && strcmp(keyword.c_str(),"BIM")!=0 && strcmp(keyword.c_str(),"PED")!=0 && strcmp(keyword.c_str(),"TPED")!=0 && strcmp(keyword.c_str(),"BPED")!=0 && strcmp(keyword.c_str(),"BED")!=0 && strcmp(keyword.c_str(),"IMPUTE2")!=0 && strcmp(keyword.c_str(),"OUTPUTNAME")!=0) warning("Keyword " + keyword + " not recognized.");
			}
		    }//j>0
		}
	    }// end j
	} //neue zeile
      else
	{
	  if (character != '\r' && character != '\n')
	    {
	      s3[i++] = character;
	    }
	}
    } //end while
  fclose(fptr3);
  //exit(1);

  /***** Fehlermeldungen *****/
  {
    switch (filemode) {
    case 1:
      logg("File format is ped/map.");
      if (!file_exists(pedfile)) die("ped-file " + pedfile + " does not exist! Please check input files and/or selection-file!\n");
      if (!file_exists(mapfile)) die("map-file " + mapfile + " does not exist! Please check input files and/or selection-file!\n");
      if (countColumns(pedfile) != 2*countLines(mapfile)+6) die("Number of SNPs in ped-file differs from the number of SNPs in map-file!\n");
      if (pedfile.substr( pedfile.size()-3,3) != "ped") warning("ped-file " + pedfile + " has not the proper filetype. Check if file is correct!");
      if (mapfile.substr( mapfile.size()-3,3) != "map") warning("map-file " + mapfile + " has not the proper filetype. Check if file is correct!");
      break;
    case 2:
      logg("File format is tped/tfam.");
      if (!file_exists(tpedfile)) die("tped-file " + tpedfile + " does not exist! Please check input files and/or selection-file!\n");
      if (!file_exists( famfile)) die("tfam-file " +  famfile + " does not exist! Please check input files and/or selection-file!\n");
      if (countColumns(tpedfile) != 2*countLines(famfile)+4) die("Number of persons in tped-file differs from the number of persons in tfam-file!\n");
      if (tpedfile.substr(tpedfile.size()-4,4) != "tped") warning("tped-file " + tpedfile + " has not the proper filetype. Check if file is correct!");
      if ( famfile.substr( famfile.size()-4,4) != "tfam") warning("tfam-file " +  famfile + " has not the proper filetype. Check if file is correct!");
      break;
    case 3:
      logg("File format is bed/bim/fam.");
      if (!file_exists(bpedfile)) die("bped-file " + bpedfile + " does not exist! Please check input files and/or selection-file!\n");
      if (!file_exists(bmapfile)) die("bmap-file " + bmapfile + " does not exist! Please check input files and/or selection-file!\n");
      if (!file_exists( famfile)) die("fam-file "  +  famfile + " does not exist! Please check input files and/or selection-file!\n");
      if (bpedfile.substr(bpedfile.size()-3,3) != "bed") warning("bed-file " + bpedfile + " has not the proper filetype. Check if file is correct!");
      if (bmapfile.substr(bmapfile.size()-3,3) != "bim") warning("bim-file " + bmapfile + " has not the proper filetype. Check if file is correct!");
      if ( famfile.substr( famfile.size()-3,3) != "fam") warning("fam-file " +  famfile + " has not the proper filetype. Check if file is correct!");
      break;
    case 4:
      logg("File format is impute2/fam/map.");
      if (!file_exists(impute2file)) die("impute2-file " + impute2file + " does not exist! Please check input files and/or selection-file!\n");
      if (!file_exists( famfile)) die("fam-file "  +  famfile + " does not exist! Please check input files and/or selection-file!\n");
      if (!file_exists(mapfile)) die("map-file " + mapfile + " does not exist! Please check input files and/or selection-file!\n");
      if (countColumns(impute2file) != 3*countLines(famfile)+5) die("Number of persons in impute2-file differs from the number of persons in fam-file!\n");
      if (countLines(impute2file) != countLines(mapfile)) die("Number of SNPs in impute2-file differs from the number of SNPs in map-file!\n");
      if (impute2file.substr(impute2file.size()-7,7) != "impute2") warning("impute2-file " + impute2file + " has not the proper filetype. Check if file is correct!");
      if ( famfile.substr( famfile.size()-3,3) != "fam") warning("fam-file " +  famfile + " has not the proper filetype. Check if file is correct!");
      if (mapfile.substr(mapfile.size()-3,3) != "map") warning("map-file " + mapfile + " has not the proper filetype. Check if file is correct!");
      break;
    default:
      die("No input file format specified!\n");
    }

    if(liability && singleMarkerTest!=3)
      {
	singleMarkerTest=3;regression=1;
	logfile << "Note: SINGLE_MARKER was changed to 3, because of liability test.\n";
	cout << "Note: SINGLE_MARKER was changed to 3, because of liability test.\n";
      }

    if(negchoice_Selected && poschoice_Selected ){die("NEGCHOICE and POSCHOICE are not compatible!\n");}
    if(negchoice_Selected && setid ){die("NEGCHOICE and SETID are not compatible!\n");}
    if(negchoice_Selected && snplist ){die("NEGCHOICE and SNPLIST are not compatible!\n");}
    if(poschoice_Selected && setid ){die("POSCHOICE and SETID are not compatible!\n");}
    if(poschoice_Selected && snplist ){die("POSCHOICE and SNPLIST are not compatible!\n");}

    if(snp2_Selected && !snp1_Selected){die("You have to specify SNP1 as well!\n");}
    if(snp3_Selected && !snp1_Selected && !snp2_Selected){die("You have to specify SNP1 and SNP2 as well!\n");}
    if(snp3_Selected && !snp2_Selected){die("You have to specify SNP2 as well!\n");}
    if(snp3_Selected && !snp1_Selected){die("You have to specify SNP1 as well!\n");}
    if(nsnps && combilist){die("COMBILIST and SNP1,SNP2,SNP3 are mutually exclusive!\n");}

    if (markercombi2 == 1 && markercombi3 == 1)
      {
	errorfile << "Cannot select TWO_MARKER and THREE_MARKER.\n";
	logfile << "Cannot select TWO_MARKER and THREE_MARKER.\n";
	cout << "Cannot select TWO_MARKER and THREE_MARKER.\n";
	errorfile.close();logfile.close();exit(1);
      }

    if(haplo)
      {
	qt=1;
	regression=1;
	if(markercombi2==1)
	  {
	    test=13;
	  }
	else if(markercombi3==1)
	  {
	    test=14;
	  }
	else
	  {
	    logfile << "Haplotype analysis will not be conducted. Choose two or three-marker analysis.\n";
	    cout << "Haplotype analysis will not be conducted. Choose two or three-marker analysis.\n";
	    errorfile << "Haplotype analysis will not be conducted. Choose two or three-marker analysis.\n";
	  }
      }

    if(pathwayAnalysis && pathwayTest==6 && !markercombi2)
      {
	markercombi2=1;markercombi3=0;
	logfile << "TWO_MARKER internally set to 1 for pathwaytest 6.\n";
      }

    if ((markercombi2 == 1 || markercombi3 == 1) && test == -1)
      {
	errorfile << "You have to specify number of multi-marker test with keyword TEST.\n";
	logfile << "You have to specify number of multi-marker test with keyword TEST.\n";
	cout << "You have to specify number of multi-marker test with keyword TEST.\n";
	errorfile.close();logfile.close();exit(1);
      }


    if ((markercombi2 == 1 || markercombi3 == 1) && (test==15 || test==16 || test==19))
      {
	if(qt==1)
	  {
	    errorfile << "\nTests 15,16,19 not compatible with QT 1.\n";
	    logfile << "\nTests 15,16,19 not compatible with QT 1.\n";
	    cout << "\nTests 15,16,19 not compatible with QT 1.\n";
	    errorfile.close();logfile.close();exit(1);
	  }
	logfile << "\n\nCase-only analysis: individuals with affection status '2' are used for interaction test.\n\n";
	cout << "\n\nCase-only analysis: individuals with affection status '2' are used for interaction test.\n\n";
      }
    if ((markercombi2 == 1 || markercombi3 == 1) && (test==17 || test==18))
      {
	qt=1;
	logfile << "\n\nNote: Tests 17, 18 are handled with linear regression, internally QT 1 is set.\n";
	logfile << "Case-only analysis: individuals with affection status '2' are used for interaction test.\n\n";
	cout << "\n\nNote: Tests 17, 18 are handled with linear regression, internally QT 1 is set.\n";
	cout << "Case-only analysis: individuals with affection status '2' are used for interaction test.\n\n";
      }

    if (plusSingleThird == 1 && markercombi2 ==1)
      {
	errorfile << "The third SNP can only be used if three-marker-analysis is selected.\n";
	logfile << "The third SNP can only be used if three-marker-analysis is selected.\n";
	cout << "The third SNP can only be used if three-marker-analysis is selected.\n";
	errorfile.close();logfile.close();exit(1);
      }

    logfile << "Missing phenotype value is " << missingpheno << ".\n";
    if(!qt) {logfile << "Missing phenotype value is also 0.\n";}


    if (plla2 && (combilist || nsnps> 0)) // Änderung569
      {
	errorfile << "PARALLELA is not compatible with any of the parameters [COMBILIST], [SNP1], [SNP2], [SNP3].\n";
	logfile << "PARALLELA is not compatible with any of the parameters [COMBILIST], [SNP1], [SNP2], [SNP3].\n";
	cout << "PARALLELA is not compatible with any of the parameters [COMBILIST], [SNP1], [SNP2], [SNP3].\n";
	errorfile.close();logfile.close();exit(1);
      }

    if (dohapfile && (PARALLELN || PARALLELA))
      {
	errorfile << "DOHAPIFLE not possible in parallel mode! Change DEFINE for PARALLELN and(or) PARALLELA Recompile!\n";
	logfile << "DOHAPIFLE not possible in parallel mode! Change DEFINE for PARALLELN and(or) PARALLELA Recompile!\n";
	cout << "DOHAPIFLE not possible in parallel mode! Change DEFINE for PARALLELN and(or) PARALLELA Recompile!\n";
        errorfile.close();logfile.close();exit(1);
      }

    if(dohapfile)
      {
	fptr10=fopen(hapfile,"w");
      }
    if(bin==1 && binsizeRare == 0 && intervalfile == " " && !vb)
      {
	cout << "For rare variant analysis, a binning method has to be chosen via either INTERVALFILE or BINSIZE_RARE!" << endl;
	logfile << "For rare variant analysis, a binning method has to be chosen via either INTERVALFILE or BINSIZE_RARE!" << endl;
	errorfile << "For rare variant analysis, a binning method has to be chosen via either INTERVALFILE or BINSIZE_RARE!" << endl;
	errorfile.close();logfile.close();exit(1);
      }
    else if(bin==1 && (binsizeRare != 0 || intervalfile != " ") && vb)
      {
	cout << "For VB analysis, all fixed-binning methods will be ignored!" << endl;
	logfile << "For VB analysis, all fixed-binning methods will be ignored!" << endl;
	binsizeRare=0;
	intervalfile=" ";
      }

    if(pathwayAnalysis)
      {
	// Es wird nur der Singlemarker für Pathwaytest 1-5 berechnet
	markercombi2 = 0;
	markercombi3 = 0;

	if(pathwayTest==5)
	  {
	    if(!qt || !regression)
	      {
		regression=1;qt=1;
		logfile << "pathwaytest 5 requires linear regression. parameters were changed accordingly.\n";
	      }
	  }
	else if(pathwayTest==6) // Pathwaytest 6
	  {
	    pathwayImpact = 1;
	    markercombi2 = 1;
	  }
      }

    if(haplo && combilist)
      {
	errorfile << "Warning. Haplotype analysis and combilist may be incompatibel. Under construction\n";
	logfile << "Warning. Haplotype analysis and combilist may be incompatibel. Under construction\n";
	cout << "Warning. Haplotype analysis and combilist may be incompatibel. Under construction\n";
      }

    if (markercombi2 == 1 && test > 8 && test != 77 && !haplo && test !=15 && test !=16 && test!=17 && test!=18)
      {
	errorfile << "This test can only be used for three-marker-analysis.\n";
	logfile << "This test can only be used for three-marker-analysis.\n";
	cout << "This test can only be used for three-marker-analysis.\n";
	errorfile.close();logfile.close();exit(1);
      }
    else if (markercombi3 == 1 && test < 9 && test > 3 && test != 77)
      {
	errorfile << "This test can only be used for two-marker-analysis.\n";
	logfile << "This test can only be used for two-marker-analysis.\n";
	cout << "This test can only be used for two-marker-analysis.\n";
	errorfile.close();
	logfile.close();
	exit(1);
      }

    if(pathwayAnalysis && haplo==1)
      {
	errorfile << "Pathwayanalysis and haplotype analysis are mutually exclusive!\n";
	logfile << "Pathwayanalysis and haplotype analysis are mutually exclusive!\n";
	cout << "Pathwayanalysis and haplotype analysis are mutually exclusive!\n";
	errorfile.close();logfile.close();exit(1);
      }

    if(haplo && (mWithGenetictop>0 || mWithSingletop>0 ))
      {
	errorfile << "Haplotype analysis is not compatible with M_WITH_SINGLETOP > 0 or M_WITH_GENETIC_IMPACT > 0 !!!\n";
	logfile << "Haplotype analysis is not compatible with M_WITH_SINGLETOP > 0 or M_WITH_GENETIC_IMPACT > 0 !!!\n";
	cout << "Haplotype analysis is not compatible with M_WITH_SINGLETOP > 0 or M_WITH_GENETIC_IMPACT > 0 !!!\n";
	errorfile.close();logfile.close();exit(1);
      }

    if(qt==1)
      {
	regression=1;
	if(singleMarkerTest==1)
	  {
	    singleMarkerTest=3;
	    logfile << "Quantitative trait analysis: singleMarkerTest was changed to " << singleMarkerTest << ".\n";
	  }
	if(singleMarkerTest==2)
	  {
	    singleMarkerTest=4;
	    logfile << "Quantitative trait analysis: singleMarkerTest was changed to " << singleMarkerTest << ".\n";
	  }
	if(markercombi2==1 && test==1)
	  {
	    test=4;
	    regression=1;
	    logfile << "Quantitative trait analysis: 2-Marker-Test was changed to " << test << ".\n";
	  }
	if(markercombi2==1 && test==2)
	  {
	    test=6;
	    regression=1;
	    logfile << "Quantitative trait analysis: 2-Marker-Test was changed to " << test << ".\n";
	  }

	if(markercombi3==1 && test<=2)
	  {
	    errorfile << "Tests 1 and 2 cannot be combined with quantitative traits and 3-marker-analysis.\n";
	    logfile << "Tests 1 and 2 cannot be combined with quantitative traits and 3-marker-analysis.\n";
	    cout << "Tests 1 and 2 cannot be combined with quantitative traits and 3-marker-analysis.\n";
	    errorfile.close();
	    logfile.close();
	    exit(1);
	  }
      }

    if (male == 1 && female == 1) {
      male = 0;
      female = 0;
      sstm << "Cannot select MALE only and FEMALE only. Selection ignored!";
      error(sstm);
    }

    if (((male == 0 && female == 1) || (male == 1 && female == 0)) && sexcov==1)
      {
	sexcov = 0;
	cout << "Note: SEXCOV is ignored because ONLY_MALE or ONLY_FEMALE is selected\n";
	logfile << "Note: SEXCOV is ignored because ONLY_MALE or ONLY_FEMALE is selected\n";
	errorfile << "Note: SEXCOV is ignored because ONLY_MALE or ONLY_FEMALE is selected\n";
      }

    if (markercombi2 == 0 && markercombi3 == 0 && pathwayImpact == 1)
      {
	errorfile << "To use the PATHWAY option TWO_MARKER or THREE_MARKER shall be selected.\n";
	logfile << "To use the PATHWAY option TWO_MARKER or THREE_MARKER shall be selected.\n";
	cout << "To use the PATHWAY option TWO_MARKER or THREE_MARKER shall be selected.\n";
	errorfile.close();logfile.close();exit(1);
      }

    if (combilist==1 && (mWithSingletop >0 || mWithGenetictop >0))
      {
	errorfile << "Warning!! Priorities M_WITH_SINGLETOP M_WITH_GENETICTOP are not in effect when combined with COMBILIST!!!\n\n";
	logfile << "Warning!! Priorities M_WITH_SINGLETOP M_WITH_GENETICTOP are not in effect when combined with COMBILIST!!!\n\n";
	cout << "Warning!! Priorities M_WITH_SINGLETOP M_WITH_GENETICTOP are not in effect when combined with COMBILIST!!!\n\n";
	mWithSingletop =0; mWithGenetictop =0;
	//errorfile.close();
      }

    if (combilist==1 && pathwayImpact == 1 )
      {
	errorfile << "PATHWAY 1 cannot be combined with COMBILIST!!!\n";
	logfile << "PATHWAY 1 cannot be combined with COMBILIST!!!\n";
	cout << "PATHWAY 1 cannot be combined with COMBILIST!!!\n";
	errorfile.close();logfile.close();exit(1);
      }

    if (combilist==1 && pathwayAnalysis == 1 )
      {
	errorfile << "COMBILIST and PATHWAYANALYSIS are mutually exclusive!!!\n";
	logfile << "COMBILIST and PATHWAYANALYSIS are mutually exclusive!!!\n";
	cout << "COMBILIST and PATHWAYANALYSIS are mutually exclusive!!!\n";
	errorfile.close();logfile.close();exit(1);
      }
    if (markercombi2==1 && nsnps+mWithSingletop>2 || markercombi3==1 && nsnps+mWithSingletop>3)
      {
	errorfile << "Error! Too many snps selected via SNP1, SNP2, SNP3 and M_WITH_SINGLETOP.\n";
	logfile << "Error! Too many snps selected via SNP1, SNP2, SNP3 and M_WITH_SINGLETOP.\n";
	cout << "Error! Too many snps selected via SNP1, SNP2, SNP3 and M_WITH_SINGLETOP.\n";
	errorfile.close();logfile.close();exit(1);
      }
    if (sexcov && singleMarkerTest<3) warning("SEXCOV is ignored. Choose a singleMarkerTest > 2!");
    if (ibs_SD_Selected && !doIbs){error("You have to select DOIBS to activate ibs_SD_outlier and ibs_SD_relatives.");}

#if RARE
    if (bin && rare_stratify && stratify) {
      die("You cannot choose genome-wide stratification- and local stratification-correcting methods at once. Please deactivate one of those options.");
    }
    if(rarepretest>nsim)
      {
	cout  <<"Argument to RARE_PRETEST can not be greater than argument to SIMULATION."<<endl;
	logfile  <<"Argument to RARE_PRETEST can not be greater than argument to SIMULATION."<<endl;
	errorfile<<"Argument to RARE_PRETEST can not be greater than argument to SIMULATION."<<endl;
	logfile.close();errorfile.close();exit(1);
      }
    if(rarepretest==nsim && bin && nsim!=0)
      {
	cout   <<"Argument of RARE_PRETEST is equal to argument of SIMULATIONS, pretest has no effect."<<endl;
	logfile<<"Argument of RARE_PRETEST is equal to argument of SIMULATIONS, pretest has no effect."<<endl;
      }
#endif

    if(choosePrintBeta > 0 && regression)
      {
	printBeta=1;
      }

    for(i=0;i<maxIndexCov;i++)
      {
	ncov+=cov[i];
      }
    if(ncov>0)
      {

	printf("%d covariates selected.\n",ncov);
	logfile << ncov << " covariates selected.\n";

	if((markercombi2==1 || markercombi3 ==1) && test <=2)
	  {
	    errorfile << "Warning! Selected multi-marker test does not use covariates!\n";
	    logfile << "Warning! Selected multi-marker test does not use covariates!\n";
	    cout << "Warning! Selected multi-marker test does not use covariates!\n";
	  }
	if(singleMarkerTest <=2)
	  {
	    errorfile << "Warning! Selected single-marker test does not use covariates!\n";
	    logfile << "Warning! Selected single-marker test does not use covariates!\n";
	    cout << "Warning! Selected single-marker test does not use covariates!\n";
	  }
      }
    if (stratify_struct != "") {
      if (stratify_struct == "CLUSTER") stratify = 4;
      else {
	    if (stratify_valid == "VICINITY")     stratify = 8;
        else if (stratify_valid == "CLUSTER") stratify = 5;
        else if (stratify_valid == "IBS")     stratify = 1;
	    if (stratify_struct == "GROUP")       stratify++;
        else if (stratify_struct == "GROUP2") stratify+=2;
      }
    }
    if (cluster_covars) {
      //        if (covariate || snpCov) die("Inputfile- or SNP-covariates cannot be used simultaneously with cluster covariates!\n");
      if (singleMarkerTest < 3) error("Selected single-marker test does not use covariates!\n");
      if (!stratify && !familyCluster) warning("CLUSTER_COVARS has no effect without usage of STRATIFYing method!");
    }
    if (familyCluster && stratify) {
      stratify = 0;
      warning("FAMCLUSTER overrides STRATIFY!");
    }
    if (group_test) {
      if (singleMarkerTest!=1) notice("Changed SINGLEMARKER to 1.");
      singleMarkerTest = 1;
    }
    if (group_test==2) {
      if (!nsim) error("No SIMULATION given for GROUP_TEST 2. Singlemarker analysis cannot be performed.");
      notice("Switching off output of asymptotic singlemarker results for MCAT(2)-test.");
      SM_showresults = false;
    }
#if RARE
    if (rare_stratify) {
      if (rare_strat_maft < 0.0) rare_strat_maft = raref;
    }
#endif
  }

  /// Ausgabefiles erstellen: OUTPUTNAME_*
  singlemarkerfile = outputname + "Singlemarker.txt";
  singlemarkerTopfile = outputname + "SinglemarkerTop.txt";
  markerCombi2file = outputname + "BestMarkerCombi2.txt";
  markerCombi2Detailsfile = outputname + "BestMarkerCombi2Details.txt";
  markerCombi3file = outputname + "BestMarkerCombi3.txt";
  markerCombi3Detailsfile = outputname + "BestMarkerCombi3Details.txt";
  toplistMCfile = outputname + "ToplistMC.txt";
  rarefile = outputname + "Rare.txt";
  rarefileinter = outputname + "RareInter.txt";
  if(optimalrare==0){
    rarefileVB = outputname + "RareVB.txt";
    rarefileVBpermstat = outputname + "RareVB.perm_p";
    rarefileVBgraph = outputname + "RareVB.graph";
    rarefileVB_BMP = outputname + "RareVB.bmp";
  }
  else if(optimalrare==1){
    rarefileVB = outputname + "RareVBVT.txt";
    rarefileVBpermstat = outputname + "RareVBVT.perm_p";
    rarefileVBgraph = outputname + "RareVBVT.graph";
    rarefileVB_BMP = outputname + "RareVBVT.bmp";
  }
  SetIDfile = outputname + ".SetID";
  if(intervaleditor==1){
    if(strcmp(intervalfile.c_str()," ")!=0){
      SetIDfile = outputname +"_"+intervalfile.substr(intervalfile.find_last_of('/')+1)+".SetID";
    }
    else{
      SetIDfile = outputname + ".SetID";
    }
    }
  rareTopfile = outputname + "RareTopfile.txt";
  liabilityfile = outputname + "Liability.txt";


  /***** Personen einlesen *****/
  struct PERSON *helpperson = NULL;
  uint32_t nperson2read = 0;
  {
    switch (filemode) {
    case 1:
      logg("\nReading pedigree data from " + pedfile + " ...");
      nlinestfam = read_fam(pedfile, qt, helpperson);
      sstm <<  pedfile << " contains " << nlinestfam << " individuals."; logg(sstm);
      break;
    case 2:
      logg("\nReading pedigree data from " + famfile + " ...");
      nlinestfam = read_fam(famfile, qt, helpperson);
      sstm <<  famfile << " contains " << nlinestfam << " individuals."; logg(sstm);
      break;
    case 3:
      logg("\nReading pedigree data from " + famfile + " ...");
      nlinestfam = read_fam(famfile, qt, helpperson);
      sstm <<  famfile << " contains " << nlinestfam << " individuals."; logg(sstm);
      break;
    case 4:
      logg("\nReading pedigree data from " + famfile + " ...");
      nlinestfam = read_fam(famfile, qt, helpperson);
      sstm <<  famfile << " contains " << nlinestfam << " individuals."; logg(sstm);
      break;
    }

    /** cross-checks on sex and aff **/
    for (int i = 0; i < nlinestfam; i++) {
      if (helpperson[i].sex != 1 && helpperson[i].sex != 2) {
	helpperson[i].sex = 0;
      }
      if ( (helpperson[i].aff[thread] != 2 && helpperson[i].aff[thread] != 1)
	   || (qt && helpperson[i].qtaff[thread]==missingpheno)
	   || (male   == 1 && helpperson[i].sex != 1)
	   || (female == 1 && helpperson[i].sex != 2)
	   || (test > 2 && sexcov == 1 && helpperson[i].sex == 0) ) {
	helpperson[i].aff[thread] = 0;
	helpperson[i].qtaff[thread] = 0;
	helpperson[i].qcin = false;
      }
    }
  }

  /***** Personen selektieren *****/
  {
    if (personlist == 1) {  /// personlist einlesen und selektieren
      ncolumnsperson = countColumns(personfile);
      if (ncolumnsperson == -1) die("File " + personfile + " does not exist!");
      if (ncolumnsperson == -2) die("First line of " + personfile + " is empty! Please check infile!");
      if (ncolumnsperson <   2) die("File " + personfile + " needs at least 2 columns!");
      character = '0';
      i = 0;
      fptr11 = fopen(personfile.c_str(), "r");
      for (int k=0; k<nlinestfam; k++) helpperson[k].analysis_in = 0;
      while ((feof(fptr11)) == 0) {
	character = fgetc(fptr11);
	if (character == '\n' && i > 1){
	  bool found = false;
	  s10[i] = '\0';
	  i = 0;
	  nlinesperson++;
	  for (j = 0; j < ncolumnsperson; j++) {
	    if (j == 0){
	      str33 = strtok(s10, ",\t; ");
	      str_fid = (char *) realloc(str_fid, (strlen(str33) + 1) * sizeof(char));
	      strcpy(str_fid,str33);
	      if (!str33) { sstm << "\nPerson file: Not enough columns in line " << nlinesperson;  die(sstm); }
	    } else {
	      str33 = strtok(NULL, ",\t; ");
	      if (!str33) { sstm << "\nPerson file: Not enough columns in line " << nlinesperson; die(sstm); }
	      if(j==1) {
		str_pid = (char *) realloc(str_pid, (strlen(str33) + 1) * sizeof(char));
		strcpy(str_pid,str33);
		for (int k=0; k<nlinestfam; k++) {
		  if (strcmp(helpperson[k].fid, str_fid) == 0 && strcmp(helpperson[k].pid, str_pid) == 0){
		    helpperson[k].analysis_in = 1;
		    nperson2read++;
		    found = true;
		    break;
		  }
		}
	      }
	    }
	  }
	  if (!found) {
	    sstm << "Person fid:" << str_fid << ", pid:" << str_pid << " not found in sample!";
	    warning(sstm);
	  }
	} else if (character != '\r') s10[i++] = character;
      } /// while
      fclose(fptr11);
    } else nperson2read = nlinestfam;
    if (nperson2read == 0) die("After specification (PERSONLIST) no persons are left.");
    if (nperson2read > (uint32_t)nlinestfam) die("Odd things happening! Watch out for ghosts!!!");
    if (nperson2read < (uint32_t)nlinestfam) { sstm << nperson2read << " persons were pre-selected by PERSONLIST " << personfile << "."; logg(sstm); }

    /** numbers **/
    ncases=0; ncontrols=0; nrest=0;
    int nmales=0;
    int nfemales=0;
    int nnosex=0;
    for (int i = 0; i < nlinestfam; i++) {
      if (!helpperson[i].analysis_in) continue;
      switch (helpperson[i].aff[thread]) {
      case 1:  ncontrols++; break;
      case 2:  ncases++;    break;
      default: nrest++;
      }
      switch (helpperson[i].sex) {
      case 1:  nmales++;   break;
      case 2:  nfemales++; break;
      default: nnosex++;
      }
    }
    if (qt) sstm << "individuals: " << nperson2read-nrest << ", notUsed: " << nrest;
    else    sstm << "cases: " << ncases << ", controls: " << ncontrols << ", notUsed: " << nrest;
    sstm << "\nmales: " << nmales << ", females: " << nfemales << ", nosex: " << nnosex;
    logg(sstm);
  }


  /***** SNPs einlesen *****/
  struct MAP* helpmap = NULL;
  uint32_t nsnp2read = 0;
  {
    switch (filemode) {
    case 1:
      logg("\nReading SNP data from " + mapfile + " ...");
      nlinestped = read_map(mapfile, helpmap);
      sstm <<  mapfile << " contains " << nlinestped << " SNPs."; logg(sstm);
      break;
    case 2:
      logg("\nReading SNP data from " + tpedfile + " ...");
      nlinestped = read_map(tpedfile, helpmap);
      sstm << tpedfile << " contains " << nlinestped << " SNPs."; logg(sstm);
      break;
    case 3:

      logg("\nReading SNP data from " + bmapfile + " ...");
      nlinestped = read_bmap(bmapfile, helpmap, codesA);
      sstm << bmapfile << " contains " << nlinestped << " SNPs."; logg(sstm);
      break;
    case 4:
      logg("\nReading SNP data from " + mapfile + " ...");
      nlinestped = read_map(mapfile, helpmap);
      sstm << mapfile << " contains " << nlinestped << " SNPs."; logg(sstm);
      break;
    }
  }


  /***** SNPs selektieren *****/
  {
    /** SNPlist einlesen **/

    if (snplist == 1 && setid==0){
      ncolumnsSNP = countColumns(SNPfile);
      if (ncolumnsSNP == -1) { sstm << "File " << SNPfile << " does not exist\n"; die(sstm); }
      if (ncolumnsSNP == -2) { sstm << "First line of " << SNPfile << " is empty! Please check infile.\n"; die(sstm); }
      character = '0';
      i = 0;

      selectedSNP = (char**) calloc(countLines(SNPfile), sizeof(char*)); //Änderung
      fptr9 = fopen(SNPfile.c_str(), "r");
      while ((feof(fptr9)) == 0)
	{
	  character = fgetc(fptr9);
	  if (character == '\n' && i > 1)
	    {
	      s9[i] = '\0';
	      i = 0;
	      nlinesSNP++;

	      for (j = 0; j < ncolumnsSNP; j++)
		{
		  if (j == 0)
		    {
		      str99 = strtok(s9, ",\t; ");

		      if (!str99)
			{
			  errorfile << "SNP file: Not enough columns in line " << nlinesSNP << " column " << j << "\n";
			  logfile << "SNP file: Not enough columns in line " << nlinesSNP  << " column " << j << "\n";
			  cout << "SNP file: Not enough columns in line " << nlinesSNP  << " column " << j << "\n";
			  errorfile.close();
			  logfile.close();
			  exit(1);
			}

		      //Änderung
		      //rsHelpmap = "";
		      selectedSNP[list] = NULL;
		      //selectedSNP[list] = (char *) realloc(selectedSNP[list], (strlen(rsHelpmap.c_str())+1)*sizeof(char));

		      selectedSNP[list] = (char *) realloc(selectedSNP[list], (strlen(str99)+1)*sizeof(char));

		      //for (lll = 0; lll < strlen(str99); lll++)
		      //{
		      strcpy(selectedSNP[list], str99);
		      //}
		    }
		}
	      list++;
	    }
	  else
	    {
	      if (character != '\r' || character != '\n')
		{
		  s9[i++] = character;
		}
	    }
	} //end while
      fclose(fptr9);
      //        qsort_cstring(selectedSNP, 0, nlinesSNP-1);
      qsort(selectedSNP, nlinesSNP, sizeof(char*), compare_cstring);
    }
    else if (setid == 1 && intervaleditor==0) {
      ncolumnsSNP = countColumns(SNPfile);
      if (ncolumnsSNP == -1) { sstm << "File " << SNPfile << " does not exist\n"; die(sstm); }
      if (ncolumnsSNP == -2) { sstm << "First line of " << SNPfile << " is empty! Please check infile.\n"; die(sstm); }
      character = '0';
      i = 0;

      selectedSNP = (char**) calloc(countLines(SNPfile), sizeof(char*));
      fptr9 = fopen(SNPfile.c_str(), "r");
      while ((feof(fptr9)) == 0)
	{
	  character = fgetc(fptr9);

	  if (character == '\n' && ((setid==1 && i>2) || (setid==0 && i> 1)))
	    {
	      s9[i] = '\0';
	      i = 0;
	      nlinesSNP++;

	      for (j = 0; j < ncolumnsSNP; j++)
		{
		  if ((setid==0 && j == 0) || (setid == 1 && j==1))
		    {
		      str99 = strtok(s9, ",\t; ");
		      if(setid==1){
			str99 = strtok(NULL, ",\t; ");
		      }
		      if (!str99)
			{
			  errorfile << "SNP file: Not enough columns in line " << nlinesSNP << " column " << j << "\n";
			  logfile << "SNP file: Not enough columns in line " << nlinesSNP  << " column " << j << "\n";
			  cout << "SNP file: Not enough columns in line " << nlinesSNP  << " column " << j << "\n";
			  errorfile.close();
			  logfile.close();
			  exit(1);
			}

		      selectedSNP[list] = NULL;
		      //selectedSNP[list] = (char *) realloc(selectedSNP[list], (strlen(rsHelpmap.c_str())+1)*sizeof(char));
		      selectedSNP[list] = (char *) realloc(selectedSNP[list], (strlen(str99)+1)*sizeof(char));

		      // for (lll = 0; lll < strlen(str99); lll++)
		      //{
		      strcpy(selectedSNP[list], str99);
		      //}
		    }
		}
	      list++;
	    }
	  else
	    {
	      if (character != '\r' || character != '\n')
		{
		  s9[i++] = character;
		}
	    }
	} //end while
      fclose(fptr9);
      //        qsort_cstring(selectedSNP, 0, nlinesSNP-1);
      qsort(selectedSNP, nlinesSNP, sizeof(char*), compare_cstring);
    }

    /** INFOlist einlesen **/
    if (infolist == 1){
      ncolumnsSNPinfo = countColumns(infofile);
      if (ncolumnsSNPinfo == -1) { sstm << "File " << infofile << " does not exist\n"; die(sstm); }
      if (ncolumnsSNPinfo == -2) { sstm << "First line of " << infofile << " is empty! Please check infile.\n"; die(sstm); }
      character = '0';
      i = 0;
      selectedSNPinfo = (char**) calloc(countLines(infofile), sizeof(char*)); //Änderung
      INFOScore = (char**) calloc(countLines(infofile), sizeof(char*));
      fptr12 = fopen(infofile.c_str(), "r");
      while ((feof(fptr12)) == 0)
	{
	  character = fgetc(fptr12);
	  if (character == '\n' && i > 1)
	    {
	      s11[i] = '\0';
	      i = 0;
	      nlinesSNPinfo++;
	      for (j = 0; j < ncolumnsSNPinfo; j++)
		{
		  if (j == 0)
		    {
		      str11 = strtok(s11, ",\t; ");
		      if (!str11){ sstm << "\nINFO file: Not enough columns in line " << nlinesSNPinfo;  die(sstm); }
		    }
		  else
		    {
		      str22 = strtok(NULL, ",\t; ");
		      if (!str22){ sstm << "\nINFO file: Not enough columns in line " << nlinesSNPinfo;  die(sstm); } // && j < ncolumnsSNPinfo -1
		    }
		}

	      //Änderung
	      //rsHelpmap = "";
	      selectedSNPinfo[list1] = NULL;
	      //selectedSNPinfo[list1] = (char *) realloc(selectedSNPinfo[list1], (strlen(rsHelpmap.c_str())+1)*sizeof(char));
	      selectedSNPinfo[list1] = (char *) realloc(selectedSNPinfo[list1], (strlen(str11)+1)*sizeof(char));

	      INFOScore[list1] = NULL;
	      //INFOScore[list1] = (char *) realloc(INFOScore[list1], (strlen(rsHelpmap.c_str())+1)*sizeof(char));
	      INFOScore[list1] = (char *) realloc(INFOScore[list1], (strlen(str22)+1)*sizeof(char));

	      //for (lll = 0; lll < strlen(str11); lll++)
	      //{
	      strcpy(INFOScore[list1],str22);
	      if(atof(INFOScore[list1])>=infoScore) strcpy(selectedSNPinfo[list1], str11);
	      //}
	      list1++;
	    }
	  else
	    {
	      if (character != '\r' || character != '\n')
		{
		  s11[i++] = character;
		}
	    }

	} //end while
      fclose(fptr12);
      //        qsort_cstring(selectedSNP, 0, nlinesSNP-1);
      qsort(selectedSNPinfo, nlinesSNPinfo, sizeof(char*), compare_cstring);
    }



    /** selektieren **/
    /// analysis-in selection (matching follows analysis, if not separately defined)
    if (choice == 1) {  /// poschoice
      for (int i=0; i<nlinestped; i++) {
	helpmap[i].analysis_in = false;
	helpmap[i].matching_in = false;
	for (int j = 0; j < reg; j++) {
	  if (!strcmp(helpmap[i].chr, region[j].chr) && helpmap[i].pos >= region[j].begin && helpmap[i].pos <= region[j].end) {
	    helpmap[i].analysis_in = true;
	    helpmap[i].matching_in = true;
	    break;
	  }
	}
      }
    } else if (choice == 2) { /// negchoice
      for (int i=0; i<nlinestped; i++) {
	helpmap[i].analysis_in = true;
	helpmap[i].matching_in = true;
	for (int j = 0; j < reg; j++) {
	  if (!strcmp(helpmap[i].chr, region[j].chr) && helpmap[i].pos >= region[j].begin && helpmap[i].pos <= region[j].end) {
	    helpmap[i].analysis_in = false;
	    helpmap[i].matching_in = false;
	    break;
	  }
	}
      }
    } else if (choice == 3) {  /// SNP list
      for (int i=0; i<nlinestped; i++) {
	helpmap[i].analysis_in = false;
	helpmap[i].matching_in = false;
	//			if (bsearch_cstring(selectedSNP, helpmap[i].rs, 0, nlinesSNP-1) > -1) {
	if (bsearch(&helpmap[i].rs, selectedSNP, nlinesSNP, sizeof(char*), compare_cstring)) {
	  helpmap[i].analysis_in = true;
	  helpmap[i].matching_in = true;
	}
      }
    } else { // no choice
      for (int i=0; i<nlinestped; i++) {
	helpmap[i].analysis_in = true;
	helpmap[i].matching_in = true;
      }
    }
    if (infolist == 1) {  /// Info list
      for (int i=0; i<nlinestped; i++) {
	if (!bsearch(&helpmap[i].rs, selectedSNPinfo, nlinesSNPinfo, sizeof(char*), compare_cstring)) {
	  helpmap[i].analysis_in = false;
	  helpmap[i].matching_in = false;
	}
      }
    }
    /// matching-in selection
    if (choiceM == 1) {  /// poschoice
      for (int i=0; i<nlinestped; i++) {
	helpmap[i].matching_in = false;
	for (int j = 0; j < regM; j++) {
	  if (!strcmp(helpmap[i].chr, regionM[j].chr) && helpmap[i].pos >= regionM[j].begin && helpmap[i].pos <= regionM[j].end) {
	    helpmap[i].matching_in = true;
	    break;
	  }
	}
      }
    }
    /// total selection
    for (int i=0; i<nlinestped; i++) {
      if (helpmap[i].analysis_in || helpmap[i].matching_in)	{
	helpmap[i].in = 1;
	nsnp2read++;
      } else {
	helpmap[i].in = 0;
      }
    }
    if (nsnp2read == 0) die("After specification (POSCHOICE, NEGCHOICE, SNPLIST, SETID, INFOLIST) no SNPs are left.");
    if (nsnp2read > (uint32_t)nlinestped) die("Odd things happening! Watch out for ghosts!!!");
    if (nsnp2read < (uint32_t)nlinestped) { sstm << nsnp2read << " SNPs were pre-selected by " << ( choice==1&&infolist==0  ? "POSCHOICE" : choice==2&&infolist==0  ? "NEGCHOICE" : choice==3&&setid==0&&infolist==0  ? "SNPfile "+SNPfile : choice==3&&setid==1&&infolist==0  ? "SETID "+SetIDfile : choice==0&&infolist==1 ? "INFOfile "+infofile : choice==1&&infolist==1 ? "POSCHOICE and INFOfile "+infofile : choice==2&&infolist==1 ? "NEGCHOICE and INFOfile "+infofile : choice==3&&setid==0&&infolist==1 ? "SNPfile "+SNPfile +" and INFOfile "+infofile : choice==3&&setid==1&&infolist==1 ? "SETID "+SetIDfile +" and INFOfile "+infofile : 0) << "."; logg(sstm); }

    free(region); region = NULL;
    free(regionM); regionM = NULL;
    free(selectedSNP); selectedSNP = NULL;
    free(selectedSNPinfo); selectedSNPinfo = NULL;
    free(INFOScore); INFOScore = NULL;
  }

  if(dosage){
    genoWeights = calloc3Dfloat(nsnp2read,nperson2read,3);
  }
  /***** Genotypen einlesen *****/
  {
    /** read genotypes from TPED-file to binary coding **/
    switch (filemode) {
    case 1:
      sstm << "\nReading genotypes from " <<  pedfile << " ..."; logg(sstm);
      nwordsSNPs =  read_ped( pedfile, helpmap, helpperson, nsnp2read, nperson2read, codesA, map, person, BinSNPs);
      break;
    case 2:
      sstm << "\nReading genotypes from " << tpedfile << " ..."; logg(sstm);
      nwordsSNPs = read_tped(tpedfile, helpmap, helpperson, nsnp2read, nperson2read, codesA, map, person, BinSNPs);
      break;
    case 3:
      sstm << "\nReading genotypes from " << bpedfile << " ..."; logg(sstm);
      nwordsSNPs = read_bped(bpedfile, helpmap, helpperson, nsnp2read, nperson2read, nlinestfam, nlinestped, codesA, map, person, BinSNPs);
      break;
    case 4:
      sstm << "\nReading genotypes from " << impute2file << " ..."; logg(sstm);
      nwordsSNPs = read_impute2(impute2file, helpmap, helpperson, nsnp2read, nperson2read, codesA, map, person, BinSNPs, probRequired,genoWeights,dosage);
    }

    /** delete helpmap **/
    for (int i=0; i<nlinestped; i++) {
      if (!helpmap[i].in) {
	free(helpmap[i].rs); helpmap[i].rs=NULL;
      }
    }
    free(helpmap); helpmap = NULL;
    nlinestped = nsnp2read;

    /** delete helpperson **/
    for (int i=0; i<nlinestfam; i++) {
      if (!helpperson[i].analysis_in) {
	free(helpperson[i].fid); helpperson[i].fid=NULL;
	free(helpperson[i].pid); helpperson[i].pid=NULL;
	free(helpperson[i].vid); helpperson[i].vid=NULL;
	free(helpperson[i].mid); helpperson[i].mid=NULL;
      }
    }
    free(helpperson); helpperson = NULL;
    nlinestfam = nperson2read;

    /** put minor allele first **/
    BinSNPsGenderFlags = getGender_bin_snp_gen(person, nlinestfam, nwordsSNPs, false);
    BinSNPsCCFlags     =     getCC_bin_snp_gen(person, nlinestfam, nwordsSNPs, false);
    struct COUNTS oneCount;
    for (int k=0; k<nlinestped; k++) {
      initCounts(&oneCount, -1, regression);
      updateCounts(&oneCount, nwordsSNPs, BinSNPs[k], BinSNPsCCFlags, BinSNPsGenderFlags, false);
      if (oneCount.AA_Ca+oneCount.AA_Co > oneCount.BB_Ca+oneCount.BB_Co) {
	char *char1;
	uint64_t wordhelp;
	char1 = codesA[k].a1;
	codesA[k].a1 = codesA[k].a2;
	codesA[k].a2 = char1;
	for (uint32_t j=0; j<nwordsSNPs; j++) {
	  wordhelp = BinSNPs[k][j][1];
	  BinSNPs[k][j][1] = BinSNPs[k][j][3];
	  BinSNPs[k][j][3] = wordhelp;
	}
      }  // if
    }  // for k
    delete_2dim(BinSNPsGenderFlags);
    delete_2dim(BinSNPsCCFlags);

    /** check if map is sorted **/
    for (int i=1; i<nlinestped; i++) {
      if ( atoi(map[i].chr) < atoi(map[i-1].chr) || ( atoi(map[i].chr)==atoi(map[i-1].chr) && map[i].pos<map[i-1].pos ) ) {
	mapsorted = false;
	break;
      }
    }
    if (!mapsorted) {
      if (haplo) die("\nInput files have to be ordered by chromosome and position for haplotype analysis!");
      if (bin)   die("\nInput files have to be ordered by chromosome and position for rare variant analysis!");
    } else {
      chrPositions = new struct REGION[27];
      for (uint8_t i=0; i<27; i++) {
	    sprintf(chrPositions[i].chr, "%d", i);
	    chrPositions[i].begin = -1;
	    chrPositions[i].end   = -1;
      }
      for (int32_t chr,p=0; p<nlinestped; p++) {
	    chr = atoi(map[p].chr);
	    if (chrPositions[chr].begin==-1) chrPositions[chr].begin = p;
	    chrPositions[chr].end = p;
      }
    }
  }


  /// Covariatefile einlesen ///
  if (covariate == 1)
    {
      sstm << "\nReading covariate file " << covariatefile << " ...";
      logg(sstm);
      cout << "Note: It is recommended that covariate file has column headings.\n";

      ncolumnscovariate = countColumns(covariatefile);
      if (ncolumnscovariate == -1)
	{

	  errorfile << "Covariatefile " << covariatefile << " does not exist!\n";
	  logfile << "Covariatefile " << covariatefile << " does not exist!\n";
	  cout << "Covariatefile " << covariatefile << " does not exist!\n";
	  errorfile.close();
	  logfile.close();
	  exit(1);
	}
      else if (ncolumnscovariate == -2)
	{
	  errorfile << "First line of " <<  covariatefile << " is empty! Please check infile.\n";
	  logfile << "First line of " <<  covariatefile << " is empty! Please check infile.\n";
	  cout << "First line of " <<  covariatefile << " is empty! Please check infile.\n";
	  errorfile.close();
	  logfile.close();
	  exit(1);
	}


      icov=0;
      while(icov<maxIndexCov)
	{
	  sum +=cov[icov];

	  if(cov[icov]==1 && (icov+1 > ncolumnscovariate-2))
	    {
	      errorfile << "Invalid covariate " << icov+1 << " chosen. Covariate file has only " << ncolumnscovariate-2 << " data columns\n";
	      logfile << "Invalid covariate " << icov+1 << " chosen. Covariate file has only " << ncolumnscovariate-2 << " data columns\n";
	      cout << "Invalid covariate " << icov+1 << " chosen. Covariate file has only " << ncolumnscovariate-2 << " data columns\n";
	      errorfile.close();logfile.close();exit(1);
	    }
	  icov++;
	}

      if (sum > (ncolumnscovariate-2))
	{
	  errorfile << "Covariate file: Too many covariates are chosen\n";
	  logfile << "Covariate file: Too many covariates are chosen\n";
	  cout << "Covariate file: Too many covariates are chosen\n";
	  errorfile.close();logfile.close();exit(1);

	}


      //read Covariate-File
      character = '0';
      i = 0;
      int nextk=0; int order=1;//NEU


      fptr6 = fopen(covariatefile.c_str(), "r");

      while ((feof(fptr6)) == 0)
	{
	  character = fgetc(fptr6);
	  if (character == '\n' && i > 1)
	    {
	      s6[i] = '\0';
	      i = 0;
	      nlinescovariate++;
	      for (j = 0; j < maxIndexCov+2; j++) // Speichern der einzelnen Spalten
		{
		  //cout << j << "\n";
		  if (j == 0)
		    {
		      str66 = strtok(s6, ",\t; ");
		      str_fid = (char *) realloc(str_fid, (strlen(str66) + 1) * sizeof(char)); //NEU
		      strcpy(str_fid,str66); //NEU

		      if (!str66)
			{
			  errorfile << "Covariate file: Not enough columns in line " << nlinescovariate << " column " << j << "\n";
			  logfile << "Covariate file: Not enough columns in line " << nlinescovariate << " column " << j << "\n";
			  cout << "Covariatefile: Not enough columns in line " << nlinescovariate << " column " << j << "\n";
			  errorfile.close();logfile.close();exit(1);
			}
		    }
		  else
		    {
		      str66 = strtok(NULL, ",\t; ");
		      if(j==1) //NEU
			{
			  str_pid = (char *) realloc(str_pid, (strlen(str66) + 1) * sizeof(char)); //NEU
			  strcpy(str_pid,str66);
			} //NEU


		      if (!str66 && j < maxIndexCov+2 -1)
			{
			  errorfile << "Covariate file: Not enough columns in line " << nlinescovariate << " column " << j << "\n";
			  logfile << "Covariate file: Not enough columns in line " << nlinescovariate << " column " << j << "\n";
			  cout << "Covariate file: Not enough columns in line " << nlinescovariate << " column " << j << "\n";
			  errorfile.close();logfile.close();exit(1);
			}

		      if (j > 1)
			{
			  found=0; //NEU
			  for(k=nextk;k<nlinestfam;k++) //NEU
			    {
			      person[k].covin= (unsigned char *) realloc(person[k].covin, maxIndexCov*sizeof(unsigned char));
			      if (!person[k].covin) die("memory allocation error in person[k].covin");
			      person[k].cov= (double *) realloc(person[k].cov, maxIndexCov*sizeof(double));
			      if (!person[k].cov) die("memory allocation error in person[k].cov");

			      if( strcmp(person[k].fid,str_fid)== 0 && strcmp(person[k].pid,str_pid)== 0 && person[k].qcin == 1)
				{
				  if( strcmp(minstring,str66)!= 0 && strcmp(xstring,str66)!= 0 && atof(str66)!= missingpheno)
				    {
				      person[k].cov[j-2] = atof(str66);
				      person[k].covin[j-2] = 1;
				      if(cov[j-2]==1){person[k].allcovin += 1;}
				      found=1;
				      break;
				    }
				} //NEU
			    }//end k
			  if(1) //always with heading
			    {
			      if(k==nlinescovariate-1 && order==1 && found==1)
				{
				  nextk++;
				}
			      else{order=0;nextk=0;}
			    }
			} // (j > 1 && nlinescovariate != 1)
		    } //if j==0
		} //j-loop
	    } //if end of line
	  else
	    {
	      if (character != '\r')
		{
		  s6[i++] = character;
		}
	    }
	} //end while

      free(str_fid);

      fclose(fptr6);
      //		printf("%s has %d columns.\n", covariatefile.c_str(), ncolumnscovariate);
      //		printf("%s has %d lines (including heading).\n", covariatefile.c_str(), nlinescovariate);
      //		logfile << covariatefile.c_str() << " has " << countColumns(covariatefile) << " columns.\n";
      //		logfile << covariatefile.c_str() << " has " << nlinescovariate << " lines (including heading).\n";


      cout    << "#covs:" << ncov << ",  ";
      logfile << "#covs:" << ncov << ",  ";
      for (i=0;i<maxIndexCov;i++)
	{
	  cout    << "cov" << i+1 << ":" << cov[i] << "  ";
	  logfile << "cov" << i+1 << ":" << cov[i] << "  ";
	}
      cout    << endl;
      logfile << endl;


      for(k=0;k<nlinestfam;k++) //NEU
	{
	  if(person[k].allcovin==ncov)
	    {
	      person[k].allcovin=1;
	    }
	  else if(person[k].allcovin<ncov)
	    {
	      person[k].allcovin=0;
	    }
	  else
	    {
	      cout << "\nError in " << covariatefile.c_str() << ": At least one of the person IDs appears multiple times\n\n";
	      logfile << "\nError in " << covariatefile.c_str() << ": At least one of the person IDs appears multiple times\n\n";
	      exit(1);
	    }

	}


      for (i = 0; i < nlinestfam; i++) //NEUNEU
	{
	  if(person[i].allcovin==0 && person[i].qcin==1)
	    {
	      if(!qt)
		{
		  if(person[i].aff[thread]==1)
		    {
		      ncontrols--;
		    }
		  else
		    {
		      ncases--;
		    }
		}
	      else //qt
		{
		  ncases--;
		}
	      nrest++;
	      person[i].qcin = 0;
	      person[i].aff[thread] = 0;
	      person[i].qtaff[thread] = 0;
	    }
	}

      if (qt) sstm << "individuals: " << nlinestfam-nrest << ", notUsed: " << nrest;
      else    sstm << "cases: " << ncases << ", controls: " << ncontrols << ", notUsed: " << nrest;
      logg(sstm);
    }


  /// Annotationfile einlesen ///
  for (i = 0; i < nlinestped; i++)  // neccessary for output
    {
      map[i].gene=NULL;
      map[i].gene= (char *) realloc(map[i].gene,3*sizeof(char));
      map[i].gene[0]='-';
      map[i].gene[1]='\0';
    }


  if (annotationfile != " " && ((mWithGenetictop != 0 && (markercombi2 >=1 || markercombi3 >=1 )) || annotate))
    {


      ncolumnsinfo = countColumns(annotationfile);
      if (ncolumnsinfo == -1)
	{

	  errorfile << "File " << annotationfile << " does not exist\n";
	  logfile << "File " << annotationfile << " does not exist\n";
	  cout << "File " << annotationfile << " does not exist\n";
	  errorfile.close();
	  logfile.close();
	  exit(1);
	}
      else if (ncolumnsinfo == -2)
	{
	  errorfile << "First line of " <<  annotationfile << " is empty! Please check infile.\n";
	  logfile << "First line of " <<  annotationfile << " is empty! Please check infile.\n";
	  cout << "First line of " <<  annotationfile << " is empty! Please check infile.\n";
	  errorfile.close();
	  logfile.close();
	  exit(1);
	}

      //read Annotation-File
      character = '0';
      i = 0;
      fptr4 = fopen(annotationfile.c_str(), "r");
      while ((feof(fptr4)) == 0)

	{
	  character = fgetc(fptr4);
	  if (character == '\n' && i > 1)
	    {
	      s4[i] = '\0';
	      i = 0;
	      nlinesinfo++;

	      // Speicher anlegen
	      infomap = (struct INFOMAP *) realloc(infomap, nlinesinfo* sizeof(struct INFOMAP));
	      if (!infomap)
		{
		  errorfile << "memory allocation error in infomap\n";
		  logfile << "memory allocation error in infomap\n";
		  cout << "memory allocation error in infomap\n";
		  errorfile.close();logfile.close();exit(1);
		}

	      infomap[nlinesinfo - 1].pos = 0;

	      infomap[nlinesinfo - 1].rs = 0;
	      infomap[nlinesinfo - 1].gene = 0;
	      infomap[nlinesinfo - 1].location = 0;
	      infomap[nlinesinfo - 1].locationToGene = 0;
	      infomap[nlinesinfo - 1].codingStatus = 0;
	      infomap[nlinesinfo - 1].line = nlinesinfo-1;

	      for (j = 0; j < ncolumnsinfo; j++) // Speichern der einzelnen Spalten
		{
		  if (j == 0)
		    {
		      str44 = strtok(s4, ",\t; ");
		      info = str44;

		      infomap[nlinesinfo - 1].rs = NULL;
		      infomap[nlinesinfo - 1].rs = (char *) realloc(infomap[nlinesinfo - 1].rs, (strlen(str44) + 1)* sizeof(char));
		      if (!infomap[nlinesinfo - 1].rs)
			{
			  errorfile << "memory allocation error in infomap[nlinesinfo - 1].rs\n";
			  logfile << "memory allocation error in infomap[nlinesinfo - 1].rs\n";
			  cout << "memory allocation error in infomap[nlinesinfo - 1].rs\n";
			  errorfile.close();logfile.close();exit(1);
			}
		      strcpy(infomap[nlinesinfo - 1].rs, str44);

		      if (!str44)
			{
			  errorfile << "Annotation file: Not enough columns in line " << nlinesinfo << "\n";
			  logfile << "Annotation file: Not enough columns in line " << nlinesinfo << "\n";
			  cout << "Annotation file: Not enough columns in line " << nlinesinfo << "\n";
			  errorfile.close();logfile.close();exit(1);
			}
		    }
		  else
		    {
		      str44 = strtok(NULL, ",\t; ");
		      info = str44;
		      if (!str44 && j < ncolumnsinfo - 1)
			{
			  errorfile << "Annotation file: Not enough columns in line " << nlinesinfo << " column " << j << "\n";
			  logfile << "Annotation file: Not enough columns in line " << nlinesinfo << " column " << j << "\n";
			  cout << "Annotation file: Not enough columns in line " << nlinesinfo << " column " << j << "\n";
			  errorfile.close();logfile.close();exit(1);
			}

		      if (j == 1)
			{
			  if (info == "X")
			    {
			      infomap[nlinesinfo - 1].chr[0] = '2'; // Chromosom
			      infomap[nlinesinfo - 1].chr[1] = '3';
			      infomap[nlinesinfo - 1].chr[2] = '\0';
			    }
			  else if (info == "Y")
			    {
			      infomap[nlinesinfo - 1].chr[0] = '2'; // Chromosom
			      infomap[nlinesinfo - 1].chr[1] = '4';
			      infomap[nlinesinfo - 1].chr[2] = '\0';
			    }
			  else if (info == "XY")
			    {
			      infomap[nlinesinfo - 1].chr[0] = '2'; // Chromosom
			      infomap[nlinesinfo - 1].chr[1] = '5';
			      infomap[nlinesinfo - 1].chr[2] = '\0';
			    }
			  else if (info == "Mt" || info == "MT")
			    {
			      infomap[nlinesinfo - 1].chr[0] = '2'; // Chromosom
			      infomap[nlinesinfo - 1].chr[1] = '6';
			      infomap[nlinesinfo - 1].chr[2] = '\0';
			    }
			  else
			    {
			      infomap[nlinesinfo - 1].chr[0] = str44[0]; // Chromosom
			      if (strlen(str44) > 1) // zweistellig
				{
				  infomap[nlinesinfo - 1].chr[1] = str44[1];
				  infomap[nlinesinfo - 1].chr[2] = '\0';
				}
			      else if (strlen(str44) == 1) // einstellig
				{
				  infomap[nlinesinfo - 1].chr[1] = '\0';
				}
			    }
			}
		      else if (j == 2)
			{
			  infomap[nlinesinfo - 1].pos = (int) atof(str44);
			}
		      else if (j == (genecol-1)) // gene_symbol (char)
			{
			  infomap[nlinesinfo - 1].gene = NULL;
			  infomap[nlinesinfo - 1].gene = (char *) realloc(infomap[nlinesinfo - 1].gene, (strlen(str44)+ 1) * sizeof(char));
			  if (!infomap[nlinesinfo - 1].gene)
			    {
			      errorfile << "memory allocation error in infomap[nlinesinfo - 1].gene\n";
			      logfile << "memory allocation error in infomap[nlinesinfo - 1].gene\n";
			      cout << "memory allocation error in infomap[nlinesinfo - 1].gene\n";
			      errorfile.close();logfile.close();exit(1);
			    }
			  strcpy(infomap[nlinesinfo - 1].gene, str44);
			}
		      else if (j == 7) // location (int)
			{

			  if (info == "coding" || info == "CDS")
			    {
			      infomap[nlinesinfo - 1].location = 1;
			      //printf("here"); exit(1);
			    }
			  else if (info == "intron")

			    {
			      infomap[nlinesinfo - 1].location = 2;
			    }
			  else if (info == "3UTR")
			    {
			      infomap[nlinesinfo - 1].location = 3;
			    }
			  else if (info == "5UTR")
			    {
			      infomap[nlinesinfo - 1].location = 4;
			    }
			  else if (info == "UTR")
			    {
			      infomap[nlinesinfo - 1].location = 5;
			    }
			  else if (info == "flanking_3UTR")
			    {
			      infomap[nlinesinfo - 1].location = 6;
			    }
			  else if (info == "flanking_5UTR")
			    {
			      infomap[nlinesinfo - 1].location = 7;
			    }
			  else
			    {
			      infomap[nlinesinfo - 1].location = 0;
			    }
			}
		      else if (j == 8) // location_relative_to_gene (int)

			{
			  if (str44[0] == '[') // Zeichen ist nicht alphanumerisch
			    {
			      infomap[nlinesinfo - 1].locationToGene = 0;
			    }
			  else
			    {
			      infomap[nlinesinfo - 1].locationToGene = (int) atof(str44);
			    }
			}
		      else if (j == 9) // coding_status (int)
			{
			  if (info == "-1")
			    {
			      infomap[nlinesinfo - 1].codingStatus = -1;
			    }
			  if (info == "NULL")
			    {
			      infomap[nlinesinfo - 1].codingStatus = 0;

			    }
			  if (info == "SYNON")
			    {
			      infomap[nlinesinfo - 1].codingStatus = 1;
			    }
			  if (info == "[complex]")
			    {
			      infomap[nlinesinfo - 1].codingStatus = 2;
			    }
			  if (info == "NONSYN" || info == "NONSYNON")
			    {
			      infomap[nlinesinfo - 1].codingStatus = 3;
			    }
			}
		    } //j>0
		}// end j
	    } //neue zeile
	  else
	    {
	      if (character != '\r')
		{
		  s4[i++] = character;
		}
	    }
	} //end while
      fclose(fptr4);
      printf("%s has %d columns.\n", annotationfile.c_str(), ncolumnsinfo);
      printf("%s has %d lines.\n\n", annotationfile.c_str(), nlinesinfo);
      logfile << annotationfile.c_str() << " has " << ncolumnsinfo << " columns.\n";
      logfile << annotationfile.c_str() << " has " << nlinesinfo << " lines.\n\n";



      // Sortieren:
      if(mapsorted==1 || 1)
	{
	  quicksortINFOMAP(&infomap,0,nlinesinfo-1, &x1infomap, &y1infomap);

	  /*for (j = 0; j < nlinesinfo; j++)
	    {
	    cout << infomap[j].rs << "\n";
	    }*/
	}

      jstart=0;
      for (i = 0; i < nlinestped; i++)
	{
	  /*for (j = jstart; j < nlinesinfo; j++)
	    {
	    //if (strcmp(map[i].rs, infomap[j].rs) == 0)
	    if (strcmp(map[i].rs, infomap[j].rs) == 0)
	    {
	    strcpy(map[i].chr, infomap[j].chr);
	    map[i].pos = infomap[j].pos;
	    map[i].gene=(char *) calloc (strlen(infomap[j].gene)+1,sizeof(char));
	    strcpy(map[i].gene,infomap[j].gene);
	    map[i].location = infomap[j].location;
	    map[i].locationToGene = infomap[j].locationToGene;
	    map[i].codingStatus = infomap[j].codingStatus;

	    if(mapsorted==1)
	    {
	    jstart=j+1;
	    }
	    break;
	    }
	    //verbesserung auch wenn mapsorted=0
	    if (atoi(map[i].chr) < atoi(infomap[j].chr))
	    {
	    break;
	    }
	    else if (atoi(map[i].chr) == atoi(infomap[j].chr) && map[i].pos < infomap[j].pos)
	    {
	    break;
	    }
	    }*/

	  //TODO

	  j=rsNumberSearchInfo(map[i].rs, infomap,0,nlinesinfo-1); //ACHTUNG infomap hat nicht den richtigen Typ!! Neue Funktion rsNumberSearchINFOMAP nötig!

	  if(j>=0 && j <nlinesinfo)
	    {

	      //map[i].pos = infomap[j].pos;
	      map[i].gene=(char *) realloc (map[i].gene,(strlen(infomap[j].gene)+1)*sizeof(char));
	      strcpy(map[i].gene,infomap[j].gene);
	      map[i].location = infomap[j].location;
	      map[i].locationToGene = infomap[j].locationToGene;
	      map[i].codingStatus = infomap[j].codingStatus;
	    }
	}
      for (i=0; i<nlinesinfo; i++) {
	free(infomap[i].rs); infomap[i].rs=NULL;
	free(infomap[i].gene); infomap[i].gene=NULL;
      }
      free(infomap); infomap=NULL;
    }



  /// Pathwayfile einlesen ///
  if (pathwayImpact == 1 || pathwayAnalysis == 1) //or PWT ==1
    {
      /** create map2 **/
      map2 = (struct MAP2 *) calloc(nlinestped, sizeof(struct MAP2));
      if (!map2) die("memory allocation error in map2");
      for (i=0; i<nlinestped; i++) {
	map2[i].rs = (char *) calloc((strlen(map[i].rs)+1), sizeof(char));
	if (!map2[i].rs) die("memory allocation error in map2[i].rs");
	strcpy(map2[i].rs, map[i].rs);
	//map2[i].line = map[i].line;
	map2[i].line = i;
      }
      qsortmap2(&map2,0,nlinestped-1, &x1map2, &y1map2);

      ncolumnspathway = countColumns(pathwayfile);
      if (ncolumnspathway == -1)
	{

	  errorfile << "File " << pathwayfile << " does not exist\n";
	  logfile << "File " << pathwayfile << " does not exist\n";
	  cout << "File " << pathwayfile << " does not exist\n";
	  errorfile.close();
	  logfile.close();
	  exit(1);
	}
      else if (ncolumnspathway == -2)
	{
	  errorfile << "First line of " <<  pathwayfile << " is empty! Please check infile.\n";
	  logfile << "First line of " <<  pathwayfile << " is empty! Please check infile.\n";
	  cout << "First line of " <<  pathwayfile << " is empty! Please check infile.\n";
	  errorfile.close();
	  logfile.close();
	  exit(1);
	}



      //read Pathway-File
      character = '0';
      i = 0;

      fptr5 = fopen(pathwayfile.c_str(), "r");
      while ((feof(fptr5)) == 0)
	{
	  character = fgetc(fptr5);
	  if (character == '\n' && i > 1)
	    {
	      s5[i] = '\0';
	      i = 0;
	      nlinespathway++;

	      str55 = strtok(s5, ",\t; ");

	      snpstring = strtok(NULL, ",\t; ");

	      if (!snpstring)
		{
		  errorfile << "Pathway file: Not enough columns in line " << nlinespathway << " column " << j << "\n";
		  logfile << "Pathway file: Not enough columns in line " << nlinespathway << " column " << j << "\n";
		  cout << "Pathway file: Not enough columns in line " << nlinespathway<< " column " << j << "\n";
		  errorfile.close();logfile.close();exit(1);
		}

	      current=0;
	      newgene=1;
	      if(ncolumnspathway>2)
		{
		  genestring=strtok(NULL, ",\t; ");
		  genefound=0;
		  for(k=0;k<ngenes;k++)
		    {
		      if(strcmp(genestring,pwgenelist[k].gene)==0)
			{
			  genefound=1;
			  current=pwgenelist[k].nr;
			  break;
			}
		    }

		  if(!genefound)
		    {
		      ngenes++;
		      pwgenelist= (struct PWgenelist *) realloc(pwgenelist,ngenes*sizeof(struct PWgenelist));
		      if (!pwgenelist)
			{
			  errorfile << "memory allocation error in pwgenelist.\n";
			  logfile << "memory allocation error in pwgenelist.\n";
			  cout << "memory allocation error in pwgenelist.\n";
			  errorfile.close();logfile.close();exit(1);
			}

		      pwgenelist[ngenes-1].nr=ngenes;
		      pwgenelist[ngenes-1].gene=NULL;
		      pwgenelist[ngenes-1].gene=(char *) realloc(pwgenelist[ngenes-1].gene,(strlen(genestring)+1)*sizeof(char));
		      strcpy(pwgenelist[ngenes-1].gene,genestring);
		      current=ngenes;
		    }
		}		//ncolumnspathway > 2



	      if (nlist > 0 && strcmp(pathway[nlist-1].name,str55) == 0)
		{
		  rsline = rsNumberSearch(snpstring, map2, 0, nlinestped-1); // rs-Nummern werden durch die Zeilennummern ersetzt


		  if (rsline != -1)
		    {
		      pathway[nlist-1].counts++;

		      pathway[nlist-1].list= (int *) realloc(pathway[nlist-1].list,pathway[nlist-1].counts*sizeof(int));
		      if (!pathway[nlist-1].list)
			{
			  errorfile << "memory allocation error in pathway[nlist-1].list\n";
			  logfile << "memory allocation error in pathway[nlist-1].list\n";
			  cout << "memory allocation error in pathway[nlist-1].list\n";
			  errorfile.close();logfile.close();exit(1);
			}
		      pathway[nlist-1].list[pathway[nlist-1].counts-1] = rsline-1;


		      pathway[nlist-1].list_gene= (int *) realloc(pathway[nlist-1].list_gene,pathway[nlist-1].counts*sizeof(int));
		      if (!pathway[nlist-1].list_gene)
			{
			  errorfile << "memory allocation error in pathway[nlist-1].list_gene\n";
			  logfile << "memory allocation error in pathway[nlist-1].list_gene\n";
			  cout << "memory allocation error in pathway[nlist-1].list_gene\n";
			  errorfile.close();logfile.close();exit(1);
			}
		      pathway[nlist-1].list_gene[pathway[nlist-1].counts-1] = current;

		      for(k=0;k<pathway[nlist-1].counts-1;k++)
			{
			  if(current==pathway[nlist-1].list_gene[k])
			    {
			      newgene=0;break;
			    }
			}
		      if(newgene)
			{
			  pathway[nlist-1].ngenes++;
			}
		    }
		}
	      else
		{
		  //neuer Pathway: Test, ob dieser bereits vorhanden ist
		  found = 0;
		  for (j = 0; j<nlist; j++)
		    {
		      if (strcmp(str55, pathway[j].name) == 0 )
			{

			  rsline = rsNumberSearch(snpstring, map2, 0, nlinestped-1); // rs-Nummern werden durch die Zeilennummern ersetzt


			  if (rsline != -1)
			    {
			      pathway[j].counts++;

			      pathway[j].list= (int *) realloc(pathway[j].list,pathway[j].counts*sizeof(int));
			      if (!pathway[j].list)
				{
				  errorfile << "memory allocation error in pathway[j].list\n";
				  logfile << "memory allocation error in pathway[j].list\n";
				  cout << "memory allocation error in pathway[j].list\n";
				  errorfile.close();logfile.close();exit(1);
				}
			      pathway[j].list[pathway[j].counts-1] = rsline-1;

			      pathway[j].list_gene= (int *) realloc(pathway[j].list_gene,pathway[j].counts*sizeof(int));
			      if (!pathway[j].list_gene)
				{
				  errorfile << "memory allocation error in pathway[j].list_gene\n";
				  logfile << "memory allocation error in pathway[j].list_gene\n";
				  cout << "memory allocation error in pathway[j].list_gene\n";
				  errorfile.close();logfile.close();exit(1);
				}
			      pathway[j].list_gene[pathway[j].counts-1] = current;

			      for(k=0;k<pathway[nlist-1].counts-1;k++)
				{
				  if(current==pathway[nlist-1].list_gene[k])
				    {
				      newgene=0;break;
				    }
				}
			      if(newgene)
				{
				  pathway[nlist-1].ngenes++;
				}
			    }
			  found = 1;
			  break;
			}
		    }

		  if (!found) // neuer Pathway, der vorher noch nicht aufgetaucht ist
		    {
		      nlist++;
		      pathway = (struct PATHWAY *) realloc(pathway, nlist* sizeof(struct PATHWAY));
		      if (!pathway)
			{
			  errorfile << "memory allocation error in pathway\n";
			  logfile << "memory allocation error in pathway\n";
			  cout << "memory allocation error in pathway\n";
			  errorfile.close();logfile.close();exit(1);
			}

		      pathway[nlist-1].counts = 0;
		      pathway[nlist-1].counts0 = 0;
		      pathway[nlist-1].ngenes = 0;
		      pathway[nlist-1].list=NULL;
		      pathway[nlist-1].list_gene=NULL;

		      strcpy(pathway[nlist-1].name, str55);

		      rsline = rsNumberSearch(snpstring, map2, 0, nlinestped-1); // rs-Nummern werden durch die Zeilennummern ersetzt


		      if (rsline != -1)
			{
			  pathway[nlist-1].counts++;

			  pathway[nlist-1].list= (int *) realloc(pathway[nlist-1].list,pathway[nlist-1].counts*sizeof(int));
			  if (!pathway[nlist-1].list)
			    {
			      errorfile << "memory allocation error in pathway[nlist-1].list\n";
			      logfile << "memory allocation error in pathway[nlist-1].list\n";
			      cout << "memory allocation error in pathway[nlist-1].list\n";
			      errorfile.close();logfile.close();exit(1);
			    }
			  pathway[nlist-1].list[pathway[nlist-1].counts-1] = rsline-1;

			  pathway[nlist-1].list_gene= (int *) realloc(pathway[nlist-1].list_gene,pathway[nlist-1].counts*sizeof(int));

			  if (!pathway[nlist-1].list_gene)
			    {
			      errorfile << "memory allocation error in pathway[nlist-1].list_gene\n";
			      logfile << "memory allocation error in pathway[nlist-1].list_gene\n";
			      cout << "memory allocation error in pathway[nlist-1].list_gene\n";
			      errorfile.close();logfile.close();exit(1);
			    }
			  pathway[nlist-1].list_gene[pathway[nlist-1].counts-1] = current;
			  if(current!=0)
			    {
			      pathway[nlist-1].ngenes = 1;
			    }
			}
		    }
		}
	    } //neue zeile
	  else
	    {
	      if (character != '\r')
		{
		  s5[i++] = character;
		}
	    }
	} //end while

      fclose(fptr5);
      printf("%s has %d columns.\n", pathwayfile.c_str(), ncolumnspathway);
      printf("%s has %d lines.\n\n", pathwayfile.c_str(), nlinespathway);
      logfile << pathwayfile.c_str() << " has " << ncolumnspathway << " columns\n";
      logfile << pathwayfile.c_str() << " has " << nlinespathway << " lines\n\n";
      cout << "Number of different pathways: " << nlist << "\n\n";
      logfile << "Number of different pathways: " << nlist<< "\n\n";

      for (i=0; i<nlinestped; i++) {
	free(map2[i].rs); map2[i].rs = NULL;
      }
      free(map2); map2 = NULL;
      for(k=0;k<ngenes;k++) {
	free(pwgenelist[k].gene); pwgenelist[k].gene = NULL;
      }
      free(pwgenelist); pwgenelist = NULL;
    }

  /// PWT
  if (pathwayAnalysis == 1)
    {
      printtop = nlist;
      mWithSingletop = 0;
      mWithGenetictop = 0;
      genetictop = 0;

    }


  /// Modelfile einlesen ///
  //x[1] = x1 ,x[2] = x1D ,x[3] = x2 ,x[4] = x2D ,x[5] = x1x2 ,
  //x[6] = x1x2D,x[7] = x1Dx2,x[8] = x1Dx2D,x[9] = x3,x[10] = x3D,
  //x[11] = x1x3,x[12] = x1x3D ,x[13] = x1Dx3,x[14] = x1Dx3D ,x[15] = x2x3 ,
  //x[16] = x2x3D,x[17] = x2Dx3,x[18] = x2Dx3D,x[19] = x1x2x3, x[20] = x1x2x3D,
  //x[21] = x1x2Dx3 ,x[22] = x1x2Dx3D ,x[23] = x1Dx2x3 ,x[24] = x1Dx2x3D ,x[25] = x1Dx2Dx3 ,x[26] = x1Dx2Dx3D
  if (readmodel == 1)
    {
      xvec[0] = 1;
      zvec[0] = 1;
      for (i=1; i < 27; i++)
	{
	  xvec[i] = 0;
	  zvec[i] = 0;
	}

      xvec2[0] = 1; //CASEONLY18
      zvec2[0] = 1;
      for (i=1; i < 27; i++)
	{
	  xvec2[i] = 0;
	  zvec2[i] = 0;
	}

      fptr7 = fopen(modelfile.c_str(), "r"); // Öffnen der Datei

      i = 0;
      if (fptr7 == NULL)
	{
	  errorfile << "Model file " << modelfile.c_str()<< " does not exist\n";
	  logfile << "Model file " << modelfile.c_str()<< " does not exist\n";
	  cout << "Model file " << modelfile.c_str()<< " does not exist\n";
	  errorfile.close();logfile.close();exit(1);
	}

      while ((character = fgetc(fptr7)) != '\n') // Einlesen der Zeilen
	{
	  s7[i++] = character;
	  if (i > 9998)
	    {
	      errorfile << "First line of model file is too long!\n";
	      logfile << "First line of model file is too long!\n";
	      cout << "First line of model file is too long!\n";
	      errorfile.close();logfile.close();exit(1);
	    }
	}
      if (character == '\n')
	s7[i] = '\0';

      if (i == 1)
	{
	  errorfile << "First line of model file is empty! Please check infile.\n";
	  logfile << "First line of model file is empty! Please check infile.\n";
	  cout << "First line of model file is empty! Please check infile.\n";
	  errorfile.close();logfile.close();exit(1);
	}

      str77 = strtok(s7, ",\t; "); // Zeile an Trennzeichen spalten

      ncolumnsmodel++;

      do {
	str77 = strtok(NULL, ",\t; ");

	if (str77 != NULL)
	  {
	    ncolumnsmodel++;
	  } // Anzahl der Spalten ermitteln
      } while (str77);

      fclose(fptr7);


      //read Modelfile
      character = '0';
      i = 0;
      fptr7 = fopen(modelfile.c_str(), "r");
      while ((feof(fptr7)) == 0)
	{
	  character = fgetc(fptr7);
	  if (character == '\n' && i > 1)
	    {
	      s7[i] = '\0';
	      i = 0;
	      nlinesmodel++;

	      for (j = 0; j < ncolumnsmodel; j++) // Speichern der einzelnen Spalten
		{
		  if (j == 0)
		    {
		      str77 = strtok(s7, ",\t; ");

		      if (!str77)
			{
			  errorfile << "Model file: Not enough columns in line " << nlinesmodel << " column " << j << "\n";
			  logfile << "Model file: Not enough columns in line " << nlinesmodel << " column " << j << "\n";
			  cout << "Model file: Not enough columns in line " << nlinesmodel << " column " << j << "\n";
			  errorfile.close();logfile.close();exit(1);
			}
		    }
		  else
		    {
		      str77 = strtok(NULL, ",\t; ");

		      if (!str77 && j < ncolumnsmodel -1)
			{
			  errorfile << "Model file: Not enough columns in line " << nlinesmodel << " column " << j << "\n";
			  logfile << "Model file: Not enough columns in line " << nlinesmodel << " column " << j << "\n";
			  cout << "Model file: Not enough columns in line " << nlinesmodel << " column " << j << "\n";
			  errorfile.close();logfile.close();exit(1);
			}

		      if (j == 1)
			{
			  xvec[nlinesmodel-1] = atoi(str77);
			}
		      else if (j == 2)
			{
			  zvec[nlinesmodel-1] = atoi(str77);
			}
		    }
		}
	    }
	  else
	    {
	      if (character != '\r')
		{
		  s7[i++] = character;
		}
	    }
	} //end while
      fclose(fptr7);
      printf("%s has %d columns.\n", modelfile.c_str(), ncolumnsmodel);
      printf("%s has %d lines.\n\n", modelfile.c_str(), nlinesmodel);
      logfile << modelfile.c_str() << " has " << ncolumnsmodel << " columns.\n";
      logfile << modelfile.c_str() << " has " << nlinesmodel << " lines.\n\n";
    }
  else
    {
      xvec[0] = 1;
      zvec[0] = 1;
      for (i=1; i < 27; i++)
	{
	  xvec[i] = 0;
	  zvec[i] = 0;
	}
    }

  /// Combifile einlesen ///
  if (combilist == 1)
    {
      missingCombis.open(missCombis.c_str(), ios::out);

      string snpA = " ";string snpB = " ";
      int firstNotFound=1;
      /** create map2 **/
      map2 = (struct MAP2 *) calloc(nlinestped, sizeof(struct MAP2));
      if (!map2) die("memory allocation error in map2");
      for (i=0; i<nlinestped; i++) {
	map2[i].rs = (char *) calloc((strlen(map[i].rs)+1), sizeof(char));
	if (!map2[i].rs) die("memory allocation error in map2[i].rs");
	strcpy(map2[i].rs, map[i].rs);
	//map2[i].line = map[i].line;
	map2[i].line = i;
      }

      qsortmap2(&map2,0,nlinestped-1, &x1map2, &y1map2);

      fptr8 = fopen(combifile.c_str(), "r"); // Öffnen der Datei


      i = 0;
      if (fptr8 == NULL)
	{
	  errorfile << "Combi file " << combifile.c_str()<< " does not exist\n";
	  logfile << "Combi file " << combifile.c_str()<< " does not exist\n";
	  cout << "Combi file " << combifile.c_str()<< " does not exist\n";
	  errorfile.close();logfile.close();exit(1);
	}

      while ((character = fgetc(fptr8)) != '\n') // Einlesen der Zeilen
	{
	  s8[i++] = character;
	  if (i > 9998)
	    {
	      errorfile << "First line of marker file is too long!\n";
	      logfile << "First line of combi file is too long!\n";
	      cout << "First line of combi file is too long!\n";
	      errorfile.close();logfile.close();exit(1);
	    }
	}
      if (character == '\n')
	s8[i] = '\0';
      if (i == 1)
	{
	  errorfile << "First line of combi file is empty! Please check infile.\n";
	  logfile << "First line of combi file is empty! Please check infile.\n";
	  cout << "First line of combi file is empty! Please check infile.\n";
	  errorfile.close();logfile.close();exit(1);
	}

      str88 = strtok(s8, ",\t; "); // Zeile an Trennzeichen spalten
      ncolumnscombi++;

      do {
	str88 = strtok(NULL, ",\t; ");
	if (str88 != NULL)
	  {
	    ncolumnscombi++;
	  } // Anzahl der Spalten ermitteln
      } while (str88);

      fclose(fptr8);

      if (ncolumnscombi == 2 && markercombi3 == 1)
	{
	  cout << "To use the three-marker-analysis COMBIFILE must have three columns\n";
	  logfile << "To use the three-marker-analysis COMBIFILE must have three columns\n";
	  errorfile << "To use the three-marker-analysis COMBIFILE must have three columns\n";
	}

      if (ncolumnscombi == 3 && markercombi2 == 1)
	{
	  cout << "To use the two-marker-analysis COMBIFILE must have two columns\n";
	  logfile << "To use the two-marker-analysis COMBIFILE must have two columns\n";
	  errorfile << "To use the two-marker-analysis COMBIFILE must have two columns\n";
	}


      //read Combifile
      character = '0';
      i = 0;

      int combiComplete=1;
      double combifileLines = 0;
      combifileLines = countLines(combifile); // Aenderung 712


      // Änderung569
      markerTable = (struct MARKERTABLE ** )calloc(countLines(combifile),sizeof(struct MARKERTABLE *));
      for(int x=0; x<combifileLines;x++)
	{
	  markerTable[x]=(struct MARKERTABLE *)calloc(3,sizeof(struct MARKERTABLE));
	  markerTable[x][0].rs=(char *) calloc(maxSize,sizeof(char));
	  markerTable[x][1].rs=(char *) calloc(maxSize,sizeof(char));
	  markerTable[x][2].rs=(char *) calloc(maxSize,sizeof(char));
	}

      //markerPosTable = calloc2Dint(countLines(combifile),3);

      fptr8 = fopen(combifile.c_str(), "r");
      while ((feof(fptr8)) == 0)
	{
	  character = fgetc(fptr8);
	  if (character == '\n' && i > 1)
	    {
	      s8[i] = '\0';
	      i = 0;
	      nlinescombi++;

	      combiComplete=1;
	      for (j = 0; j < ncolumnscombi; j++) // Speichern der einzelnen Spalten
		{
		  if (j == 0)
		    {
		      str88 = strtok(s8, ",\t; ");

		      if (!str88)
			{
			  errorfile << "Combifile: Not enough columns in line " << nlinescombi << " column " << j << "\n";
			  logfile << "Combifile: Not enough columns in line " << nlinescombi << " column " << j << "\n";
			  cout << "Combifile: Not enough columns in line " << nlinescombi << " column " << j << "\n";
			  errorfile.close();logfile.close();exit(1);
			}

		      // Abgleich, ob in tped-File vorhanden
		      strcpy(markerTable[mlist][j].rs,dummyString);
		      markerTable[mlist][j].pos =-1; // Änderung569
		      rsline = -1;

		      rsline=rsNumberSearch(str88, map2, 0, nlinestped-1);

		      snpA=str88;
		      if (rsline == -1)
			{
			  combiComplete=0; //missingCombis << str88 << "\t";

			  if(firstNotFound)
			    {
			      errorfile << "There are SNPs in the combifile that are not available. See file CombisNotFound.\n\n";
			      cout << "There are SNPs in the combifile that are not available. See file CombisNotFound.\n\n";
			      logfile << "There are SNPs in the combifile that are not available. See file CombisNotFound.\n\n";
			    }
			  firstNotFound=0;
			}
		      else
			{
			  if(strlen(str88)+1>maxSize)
			    {
			      if( (strlen(str88)+1) > maxSizeNew) {maxSizeNew=strlen(str88)+1;}
			      markerTable[mlist][j].rs=(char *) realloc(markerTable[mlist][j].rs,(strlen(str88)+1)*sizeof(char));
			    }
			  strcpy(markerTable[mlist][j].rs,str88);
			  markerTable[mlist][j].pos=rsline; // Änderung569
			}
		    }
		  else
		    {
		      str88 = strtok(NULL, ",\t; ");

		      if (!str88 && j < ncolumnscombi -1)
			{
			  errorfile << "Combifile: Not enough columns in line " << nlinescombi << " column " << j << "\n";
			  logfile << "Combifile: Not enough columns in line " << nlinescombi << " column " << j << "\n";
			  cout << "Combifile: Not enough columns in line " << nlinescombi << " column " << j << "\n";
			  errorfile.close();logfile.close();exit(1);
			}

		      if (j == 1)
			{
			  // Abgleich, ob in tped-File vorhanden
			  strcpy(markerTable[mlist][j].rs,dummyString);
			  markerTable[mlist][j].pos = -1; // Änderung569

			  rsline=rsNumberSearch(str88, map2, 0, nlinestped-1);

			  snpB=str88;
			  if (rsline == -1)
			    {
			      combiComplete=0;

			      //if(markercombi2) {missingCombis << str88 << "\n";}
			      //else {missingCombis << str88 << "\t";}
			      if(firstNotFound)
				{
				  errorfile << "There are SNPs in the combifile that are not available.\n\n";
				  cout << "There are SNPs in the combifile that are not available.\n\n";
				  logfile << "There are SNPs in the combifile that are not available.\n\n";
				}
			      firstNotFound=0;
			    }
			  else
			    {
			      if(strlen(str88)+1>maxSize)
				{
				  if( (strlen(str88)+1) > maxSizeNew) {maxSizeNew=strlen(str88)+1;}
				  markerTable[mlist][j].rs=(char *) realloc(markerTable[mlist][j].rs,(strlen(str88)+1)*sizeof(char));
				}
			      strcpy(markerTable[mlist][j].rs,str88); // Änderung569
			      markerTable[mlist][j].pos=rsline; // Änderung569
			    }

			}
		      if (j == 2 && markercombi3 == 1)
			{
			  // Abgleich, ob in tped-File vorhanden
			  strcpy(markerTable[mlist][j].rs,dummyString);
			  markerTable[mlist][j].pos = -1; // Änderung569

			  rsline=rsNumberSearch(str88, map2, 0, nlinestped-1);

			  if (rsline == -1)
			    {
			      combiComplete=0;
			      //missingCombis << str88 << "\n";
			      if(firstNotFound)
				{
				  errorfile << "There are SNPs in the combifile that are not available.\n\n";
				  cout << "There are SNPs in the combifile that are not available.\n\n";
				  logfile << "There are SNPs in the combifile that are not available.\n\n";
				}
			      firstNotFound=0;
			    }
			  else
			    {
			      if(strlen(str88)+1>maxSize)
				{
				  if( (strlen(str88)+1) > maxSizeNew) {maxSizeNew=strlen(str88)+1;}
				  markerTable[mlist][j].rs=(char *) realloc(markerTable[mlist][j].rs,(strlen(str88)+1)*sizeof(char));
				}
			      strcpy(markerTable[mlist][j].rs,str88);
			      markerTable[mlist][j].pos=rsline; // Änderung569
			    }
			}
		    }
		}
	      mlist++;
	      if(!combiComplete)
		{
		  if(markercombi2){missingCombis << snpA << "\t" << snpB << "\n";}
		  else if(markercombi3){missingCombis << snpA << "\t" << snpB << "\t" << str88 << "\n";}
		}
	    }
	  else
	    {
	      if (character != '\r')
		{
		  s8[i++] = character;
		}
	    }

	} //end while
      fclose(fptr8);
      printf("%s has %d columns.\n", combifile.c_str(), ncolumnscombi);
      printf("%s has %d lines.\n\n", combifile.c_str(), nlinescombi);
      logfile << combifile.c_str() << " has " << ncolumnscombi << "columns.\n";
      logfile << combifile.c_str() << " has " << nlinescombi << " lines.\n\n";

      for (i=0; i<nlinestped; i++)
	{
	  free(map2[i].rs); map2[i].rs = NULL;
	}
      free(map2); map2 = NULL;
      missingCombis.close();
    }


  /// SNP-Covariates
  if (snpCov == 1)
    {
      ReadSnpCovariates(thread,npplqc, BinSNPs, maxIndexCov, numberOfSNPCov, nlinestfam, nsnpqc, person, map, SNPcovar, missingpheno, qt, nlinestped, liability,dosage,genoWeights);

      cov = (int *) realloc(cov, maxIndexCov*sizeof(int));

      for (int icov=maxIndexCov-numberOfSNPCov;icov<maxIndexCov;icov++)
	{
	  cov[icov] = 1;
	  ncov += 1;
	}
      cout << "After reading snp-covariate-file: ";
      cout << "Number of covariates and SNP covariates: " <<  ncov << "\n";
      //cout << "maxIndexCov: " <<  maxIndexCov << "\n";
    }

#if RARE
  if(bin){
    if(collinter==0){
      FISHERflag=0;
      CMATflag=0;
      COLLflag=0;
      REGRESSIONflag=0;
      FRACREGflag=0;
      COLLREGflag=0;
      if(vb==0 && FISHERtest==0 && COLLtest==0 && CMATtest==0 && REGRESSIONtest==0 && FRACREGtest==0 &&  COLLREGtest==0)
	{
	  cout << "RARE_TESTS has no valid argumens: FR, CMAT, COLL will be conducted!"<<endl;
	  logfile << "RARE_TESTS has no valid argumens: FR, CMAT, COLL will be conducted!"<<endl;
	  FISHERtest=1;
	  CMATtest=1;
	  COLLtest=1;
	  REGRESSIONtest=0;
	  FRACREGtest=0;
	  COLLREGtest=0;
	}
      cout << "\nThe following rare variant tests will be conducted:";
      logfile << "\nThe following rare variant tests will be conducted:";
    }
    if(vb){
      FISHERtest=0;
      CMATtest=0;
      COLLtest=0;
      REGRESSIONtest=0;
      FRACREGtest=0;
      COLLREGtest=0;
      if(!optimalrare){
	cout << " COLL with the VB and FT methods";
	logfile << " COLL with the VB and FT methods";
      }
      else if(optimalrare){
	cout << " COLL with the VB and VT methods";
	logfile << " COLL with the VB and VT methods";
      }
    }
    else{
      if(REGRESSIONtest==1){
	cout << " REG";
	logfile << "REG";
	REGRESSIONflag=1;
      }
      if(FRACREGtest==1){
	cout << " FRACREG";
	logfile << " FRACREG";
	FRACREGflag=1;
      }
      if(COLLREGtest==1){
	cout << " COLLREG";
	logfile << " COLLREG";
	COLLREGflag=1;
      }
      if(FISHERtest==1){
	cout << " FR";
	logfile << " FR";
	FISHERflag=1;
      }
      if(CMATtest==1){
	cout << " CMAT";
	logfile << " CMAT";
	CMATflag=1;
      }
      if(COLLtest==1){
	cout << " COLL";
	logfile << " COLL";
	COLLflag=1;
      }
    }
    cout << "."<<endl;
    logfile <<"."<<endl;

    if(nsim==0 && wilsonpretest!=0){
      cout<<"With SIMULATION 0, ADAPTIVE set to 0."<<endl;
      logfile<<"With SIMULATION 0, ADAPTIVE set to 0."<<endl;
      wilsonpretest=0;
    }
    if(nsim==0 && rareregpretest!=0){
      cout<<"With SIMULATION 0, PRESCREEN set to 0."<<endl;
      logfile<<"With SIMULATION 0, PRESCREEN set to 0."<<endl;
      rareregpretest=0;
    }
  }
#endif


  /***** Quality Control *****/
  {
    logg("\nPerforming Quality Control ...");

    nlines2Aff = nlinestfam - nrest; // nur Personen mit Affectionstatus
    nlines2AffMale=0;
    for (i = 0; i < nlinestfam; i++) {
      if(person[i].qcin && person[i].sex==1) nlines2AffMale++;
    }
    for (i = 0; i < nlinestped; i++) {
      if (!strcmp(map[i].chr, "24")) ylines++;
    }

    /// check haploids
    uint64_t nhhg=0;
    uint64_t nfyg=0;
    for (i=0; i < nlinestfam; i++) {
      uint64_t nr  = i/64;
      uint64_t pos = i%64;
      for (int j=0; j < nlinestped; j++) {
	if (!strcmp(map[j].chr, "23")) {
	  if (person[i].sex==1 && getbit64(BinSNPs[j][nr][2], pos)) {
	    setGenotype(BinSNPs[j][nr], pos, 0);
	    sstm << "Person " << person[i].pid << ":" << person[i].fid << " SNP " << map[j].chr  << ":" << map[j].rs << " set to missing (ChrX, Male, heterozygous).";
	    notice(sstm);
	    nhhg++;
	  }
	} else
	  if (!strcmp(map[j].chr, "24")) {
	    if (person[i].sex==1 && getbit64(BinSNPs[j][nr][2], pos)) {
	      setGenotype(BinSNPs[j][nr], pos, 0);
	      sstm << "Person " << person[i].pid << ":" << person[i].fid << " SNP " << map[j].chr  << ":" << map[j].rs << " set to missing (ChrY, Male, heterozygous).";
	      notice(sstm);
	      nhhg++;
	    } else
	      if (person[i].sex==2 && getGenotype(BinSNPs[j][nr], pos)) {
		setGenotype(BinSNPs[j][nr], pos, 0);
		sstm << "Person " << person[i].pid << ":" << person[i].fid << " SNP " << map[j].chr  << ":" << map[j].rs << " set to missing (ChrY, Female, non-missing).";
		notice(sstm);
		nfyg++;
	      }
	  } else
            if (!strcmp(map[j].chr, "26")) {
	      if (getbit64(BinSNPs[j][nr][2], pos)) {
		setGenotype(BinSNPs[j][nr], pos, 0);
		sstm << "Person " << person[i].pid << ":" << person[i].fid << " SNP " << map[j].chr  << ":" << map[j].rs << " set to missing (MT, heterozygous).";
		notice(sstm);
		nhhg++;
	      }
            }
      }
    }
    if (nhhg > 0) { sstm << nhhg << " heterozygous haploid genotypes; set to missing."; logg(sstm); }
    if (nfyg > 0) { sstm << nfyg << " non-missing female genotypes on Y-chromosome; set to missing."; logg(sstm); }

    /// count missings
    for (i = 0; i < nlinestfam; i++) {
      for (int j = 0; j < nlinestped; j++) {
	if ( getbit64(BinSNPs[j][i/64][0], i%64) ) {
	  if ( !strcmp(map[j].chr, "24") )  { // ChrY
	    if (person[i].sex==1) {
	      //                        person[i].missing += 1;
	      if ( person[i].aff[thread]==1 || person[i].aff[thread]==2 ) {
		map[j].missing += 1;
	      }
	    }
	  } else {
	    person[i].missing += 1;
	    if ( person[i].aff[thread]==1 || person[i].aff[thread]==2 ) {
	      map[j].missing += 1;
	    }
	  }
	}
      }
    }

    if (switchCaseControl == 1)
      {
	for (i = 0; i < nlinestfam; i++)
	  {
	    if (person[i].aff[thread] == 2)
	      {
		person[i].aff[thread] = 1;
	      }
	    else if (person[i].aff[thread] == 1)
	      {
		person[i].aff[thread] = 2;
	      }
	  }
	cout << "!!!Cases and Controls are switched!!!\n";
      }





    /// Iterative QC-algorithm
    /// starting point: average genotype missing rate taken over all SNPs and individuals
    /// In every iteration, alternately either SNPs or individuals with a missing rate worse than
    /// the average missing rate plus a user-defined missing rate difference (mrdiff) are discarded.
    /// Then, the new average missing rate is calculated and further SNPs or individuals are discarded,
    /// when their missing rate is higher than the new missing rate plus mrdiff.
    /// The algorithm terminates when there are no SNPs or individuals left that have to be discarded.
    if (skip_it_qc) logg("Skipping iterative QC algorithm.");
    else {
    if (storeall && mrdiff >= 0 && mrdiff <= 1 )
      {

	/// Missingrate pro Person
	for (i = 0; i < nlinestfam; i++)
	  {
	    nlinesy = nlinestped - ylines;
	    person[i].missingrate = (double)person[i].missing / nlinesy;
	    meanMissingPerson += person[i].missingrate;
	  }
	meanMissingPerson = meanMissingPerson / nlines2Aff;

	/// Missing rate per SNP
	for (j = 0; j < nlinestped; j++)
	  {
	    if ( !strcmp(map[j].chr, "24") )
	      {
		map[j].missingrate = (double)(map[j].missing) / (nlines2AffMale);
	      }
	    else
	      {
		//map[j].missings_co = (double)map[j].missings_co / nlines2Aff; // Control
		//map[j].missings_ca = (double)map[j].missings_ca / nlines2Aff; // Case
		map[j].missingrate = (double)(map[j].missing) / (nlines2Aff); // total Missing rate
		meanMissingSNP += map[j].missingrate;
	      }
	  }

	meanMissingSNP = meanMissingSNP / nlinesy;
	meanMissingNew = meanMissingSNP; // starting point

	nlinestpedMod = 0;
	nlinestfamMod = 0;
	while (change == 1 && nlinesy > 0)
	  {
	    change = 0;

	    /// Markierung der SNPs
	    for (j = 0; j < nlinestped; j++) {
	      if (!strcmp(map[j].chr, "24")) continue;  // ChrY
	      /// Markieren, wenn Missingrate ist > meanMissingSNP + MRDIFF
	      if (map[j].missingrate > (meanMissingNew + mrdiff) && map[j].qcin)
		{
		  change = 1;
		  nlinestpedMod++;
		  map[j].qcin = 0;
		}
	    }

	    /// Durchschnitt neu berechnen für Personen, aufgrund aktualisierter SNP-Liste
	    meanMissingPerson = 0;
	    for (i = 0; i < nlinestfam; i++) {
	      person[i].missingrate = 0;

	      for (j = 0; j < nlinestped; j++) // SNPs
		{
		  if (!strcmp(map[j].chr, "24")) continue;  // ChrY
		  if ( !map[j].qcin && getbit64(BinSNPs[j][i/64][0], i%64) )
		    {
		      mod++;
		    }
		}
	      if (person[i].qcin && person[i].missing > mod)
		{
		  person[i].missingrate = ((double)person[i].missing - mod)/ (nlinesy - nlinestpedMod);
		}
	      else
		{
		  person[i].missingrate = 0;
		}
	      meanMissingPerson += person[i].missingrate;
	      mod = 0;
	    }
	    meanMissingNew = meanMissingPerson / (nlines2Aff - nlinestfamMod);
	    sstm << "meanMissingNew Person: " << meanMissingNew;
	    logg(sstm);

	    /// Markierung der Personen
	    for (i = 0; i < nlinestfam; i++) {
	      /// Markieren, wenn Missingrate ist > meanMissingSNP + MRDIFF
	      if ((person[i].missingrate > meanMissingNew + mrdiff) && person[i].qcin)
		{
		  if (person[i].aff[thread] == 2)
		    {
		      ncases--;
		    }
		  else if (person[i].aff[thread] == 1)
		    {
		      ncontrols--;
		    }
		  nrest++;
		  change = 1;
		  nlinestfamMod++;
		  person[i].qcin = 0;
		  person[i].aff[thread] = 0;
		  person[i].qtaff[thread] = 0;
		}
	    }

	    /// Durchschnitt neu berechnen für SNPs, aufgrund aktualisierter Personen-Liste
	    meanMissingSNP = 0;
	    for (j = 0; j < nlinestped; j++) {
	      if (!strcmp(map[j].chr, "24")) continue;  // ChrY
	      for (i = 0; i < nlinestfam; i++) // Person
		{
		  if ( !person[i].qcin && getbit64(BinSNPs[j][i/64][0], i%64))
		    {
		      mod2++;
		    }
		}
	      if (map[j].qcin && (map[j].missing > mod2))
		{
		  map[j].missingrate = ((double)map[j].missing - mod2)/ (nlines2Aff - nlinestfamMod);
		}
	      else
		{
		  map[j].missingrate = 0;
		}

	      meanMissingSNP += map[j].missingrate;
	      mod2 = 0;
	    }
	    meanMissingNew = meanMissingSNP / (nlinesy - nlinestpedMod);
	    sstm << "meanMissingNew SNP: " << meanMissingNew;
	    logg(sstm);

	  }
      } //end if storeall
    }

	/// SNPs which are not in HWE in either cases or controls are removed
    BinSNPsCCFlags     =     getCC_bin_snp_gen(person, nlinestfam, nwordsSNPs, false);
    BinSNPsGenderFlags = getGender_bin_snp_gen(person, nlinestfam, nwordsSNPs, false);
    struct COUNTS oneCount;
    for (int j=0; j<nlinestped; j++) {
      initCounts(&oneCount, -1, regression);
      updateCounts(&oneCount, nwordsSNPs, BinSNPs[j], BinSNPsCCFlags, BinSNPsGenderFlags, !(strcmp(map[j].chr, "23") && strcmp(map[j].chr, "24")));
      if (!(strcmp(map[j].chr, "24") == 0) && !(strcmp(map[j].chr, "26") == 0))
	{
	  if (!(strcmp(map[j].chr, "23") == 0))
	    {
	      testHweCo = testHWE(oneCount.AA_Co, oneCount.AB_Co, oneCount.BB_Co);
	      testHweCa = testHWE(oneCount.AA_Ca, oneCount.AB_Ca, oneCount.BB_Ca);
	      //MAF
	      thismaf= (double)(2*oneCount.AA_Ca+oneCount.AB_Ca+2*oneCount.AA_Co+oneCount.AB_Co)/(2*oneCount.AA_Ca+2*oneCount.AB_Ca+2*oneCount.BB_Ca+2*oneCount.AA_Co+2*oneCount.AB_Co+2*oneCount.BB_Co);
	      // With missings assumed to be common:
	      //			    thismafa= (double)(2*oneCount.AA_Ca+oneCount.AB_Ca+2*oneCount.AA_Co+oneCount.AB_Co)/(2*oneCount.AA_Ca+2*oneCount.AB_Ca+2*oneCount.BB_Ca+2*oneCount.AA_Co+2*oneCount.AB_Co+2*oneCount.BB_Co+2*oneCount.OO_Co+2*oneCount.OO_Ca);
	      int Nalleles=(2*oneCount.AA_Ca+2*oneCount.AB_Ca+2*oneCount.BB_Ca+2*oneCount.AA_Co+2*oneCount.AB_Co+2*oneCount.BB_Co+2*oneCount.OO_Co+2*oneCount.OO_Ca);
	      int Nalleles_Co=(2*oneCount.AA_Co+2*oneCount.AB_Co+2*oneCount.BB_Co+2*oneCount.OO_Co);
	      int Nalleles_Ca=(2*oneCount.AA_Ca+2*oneCount.AB_Ca+2*oneCount.BB_Ca+2*oneCount.OO_Ca);
	      thismafa= 1/(double)Nalleles*floor((double)((2*oneCount.AA_Ca+oneCount.AB_Ca+2*oneCount.AA_Co+oneCount.AB_Co)*Nalleles)/(Nalleles-2*oneCount.OO_Co-2*oneCount.OO_Ca));
	      if(!qt)
		{
		  thiscontrolmaf=(double)(2*oneCount.AA_Co+oneCount.AB_Co)/(2*oneCount.AA_Co+2*oneCount.AB_Co+2*oneCount.BB_Co);
		  thiscasemaf=(double)(2*oneCount.AA_Ca+oneCount.AB_Ca)/(2*oneCount.AA_Ca+2*oneCount.AB_Ca+2*oneCount.BB_Ca);

		  //thiscontrolmafa=(double)(2*oneCount.AA_Co+oneCount.AB_Co)/(2*oneCount.AA_Co+2*oneCount.AB_Co+2*oneCount.BB_Co+2*oneCount.OO_Co);
		  thiscontrolmafa=1/(double)Nalleles_Co*floor(((double)(2*oneCount.AA_Co+oneCount.AB_Co)*Nalleles_Co)/(2*oneCount.AA_Co+2*oneCount.AB_Co+2*oneCount.BB_Co-2*oneCount.OO_Co));

		  // thiscasemafa=(double)(2*oneCount.AA_Ca+oneCount.AB_Ca)/(2*oneCount.AA_Ca+2*oneCount.AB_Ca+2*oneCount.BB_Ca+2*oneCount.OO_Ca);
		  thiscasemafa=1/(double)Nalleles_Ca*floor(((double)(2*oneCount.AA_Ca+oneCount.AB_Ca)*Nalleles_Ca)/(2*oneCount.AA_Ca+2*oneCount.AB_Ca+2*oneCount.BB_Ca-2*oneCount.OO_Ca));
		}
	      else
		{
		  thiscontrolmaf=thismaf;
		  thiscasemaf=thismaf;
		  thiscontrolmafa=thismafa;
		  thiscasemafa=thismafa;
		}
	    }
	  else
	    {

	      testHweCo = testHWE(oneCount.AA_Co_female,oneCount.AB_Co_female, oneCount.BB_Co_female);
	      testHweCa = testHWE(oneCount.AA_Ca_female,oneCount.AB_Ca_female, oneCount.BB_Ca_female);
	      //MAF
	      int Nalleles=(2*oneCount.AA_Ca_female+oneCount.AB_Ca_female+2*oneCount.AA_Co_female+oneCount.AB_Co_female+2*oneCount.BB_Ca_female+oneCount.AB_Ca_female+2*oneCount.BB_Co_female+oneCount.AB_Co_female+oneCount.AA_Ca_male+oneCount.AA_Co_male+oneCount.BB_Ca_male+oneCount.BB_Co_male+2*oneCount.OO_Ca_female+2*oneCount.OO_Co_female+oneCount.OO_Ca_male+oneCount.OO_Co_male);
	      int Nalleles_Co=(2*oneCount.AA_Co_female+oneCount.AB_Co_female+2*oneCount.BB_Co_female+oneCount.AB_Co_female+oneCount.AA_Co_male+oneCount.BB_Co_male+2*oneCount.OO_Co_female+oneCount.OO_Co_male);
	      int Nalleles_Ca=(2*oneCount.AA_Ca_female+oneCount.AB_Ca_female+2*oneCount.BB_Ca_female+oneCount.AB_Ca_female+oneCount.AA_Ca_male+oneCount.BB_Ca_male+2*oneCount.OO_Ca_female+oneCount.OO_Ca_male);

	      thismaf= (double)(2*oneCount.AA_Ca_female+oneCount.AB_Ca_female+2*oneCount.AA_Co_female+oneCount.AB_Co_female+oneCount.AA_Ca_male+oneCount.AA_Co_male)/(2*oneCount.AA_Ca_female+oneCount.AB_Ca_female+2*oneCount.AA_Co_female+oneCount.AB_Co_female+2*oneCount.BB_Ca_female+oneCount.AB_Ca_female+2*oneCount.BB_Co_female+oneCount.AB_Co_female+oneCount.AA_Ca_male+oneCount.AA_Co_male+oneCount.BB_Ca_male+oneCount.BB_Co_male);

	      //				thismafa= (double)(2*oneCount.AA_Ca_female+oneCount.AB_Ca_female+2*oneCount.AA_Co_female+oneCount.AB_Co_female+oneCount.AA_Ca_male+oneCount.AA_Co_male)/(2*oneCount.AA_Ca_female+oneCount.AB_Ca_female+2*oneCount.AA_Co_female+oneCount.AB_Co_female+2*oneCount.BB_Ca_female+oneCount.AB_Ca_female+2*oneCount.BB_Co_female+oneCount.AB_Co_female+oneCount.AA_Ca_male+oneCount.AA_Co_male+oneCount.BB_Ca_male+oneCount.BB_Co_male+2*oneCount.OO_Ca_female+2*oneCount.OO_Co_female+oneCount.OO_Ca_male+oneCount.OO_Co_male);
	      thismafa= 1/(double)Nalleles*floor((double)(2*oneCount.AA_Ca_female+oneCount.AB_Ca_female+2*oneCount.AA_Co_female+oneCount.AB_Co_female+oneCount.AA_Ca_male+oneCount.AA_Co_male)*Nalleles/(2*oneCount.AA_Ca_female+oneCount.AB_Ca_female+2*oneCount.AA_Co_female+oneCount.AB_Co_female+2*oneCount.BB_Ca_female+oneCount.AB_Ca_female+2*oneCount.BB_Co_female+oneCount.AB_Co_female+oneCount.AA_Ca_male+oneCount.AA_Co_male+oneCount.BB_Ca_male+oneCount.BB_Co_male-2*oneCount.OO_Ca_female-2*oneCount.OO_Co_female-oneCount.OO_Ca_male-oneCount.OO_Co_male));
	      if(!qt)
		{
		  thiscontrolmaf=(double)(2*oneCount.AA_Co_female+oneCount.AB_Co_female+oneCount.AA_Co_male)/(2*oneCount.AA_Co_female+oneCount.AB_Co_female+2*oneCount.BB_Co_female+oneCount.AB_Co_female+oneCount.AA_Co_male+oneCount.BB_Co_male);
		  thiscasemaf= (double)(2*oneCount.AA_Ca_female+oneCount.AB_Ca_female+oneCount.AA_Ca_male)/(2*oneCount.AA_Ca_female+oneCount.AB_Ca_female+2*oneCount.BB_Ca_female+oneCount.AB_Ca_female+oneCount.AA_Ca_male+oneCount.BB_Ca_male);

		  thiscontrolmafa= 1/(double)Nalleles_Co*floor((double)(2*oneCount.AA_Co_female+oneCount.AB_Co_female+oneCount.AA_Co_male)*Nalleles_Co/(2*oneCount.AA_Co_female+oneCount.AB_Co_female+2*oneCount.BB_Co_female+oneCount.AB_Co_female+oneCount.AA_Co_male+oneCount.BB_Co_male-2*oneCount.OO_Co_female-oneCount.OO_Co_male));
		  thiscasemafa=  1/(double)Nalleles_Ca*floor((double)(2*oneCount.AA_Ca_female+oneCount.AB_Ca_female+oneCount.AA_Ca_male)*Nalleles_Ca/(2*oneCount.AA_Ca_female+oneCount.AB_Ca_female+2*oneCount.BB_Ca_female+oneCount.AB_Ca_female+oneCount.AA_Ca_male+oneCount.BB_Ca_male-2*oneCount.OO_Ca_female-oneCount.OO_Ca_male));
		}
	      else{thiscontrolmaf=thismaf;thiscasemaf=thismaf;thiscontrolmafa=thismafa;thiscasemafa=thismafa;}

	    }

	  if (testHweCo < hweCo || testHweCa < hweCa)
	    {
	      map[j].qcin = 0;
	    }
	}
      else if (strcmp(map[j].chr, "24") == 0)
	{
	  if(map[j].missingrate > mrdiff){map[j].qcin = 0;}

	  int Nalleles=oneCount.AA_Ca_male+oneCount.AA_Co_male+oneCount.BB_Ca_male+oneCount.BB_Co_male+oneCount.OO_Co_male+oneCount.OO_Ca_male;
	  int Nalleles_Co=oneCount.AA_Co_male+oneCount.BB_Co_male+oneCount.OO_Co_male;
	  int Nalleles_Ca=oneCount.AA_Ca_male+oneCount.BB_Ca_male+oneCount.OO_Ca_male;

	  //MAF
	  thismaf=(double)(oneCount.AA_Ca_male+oneCount.AA_Co_male)/(oneCount.AA_Ca_male+oneCount.AA_Co_male+oneCount.BB_Ca_male+oneCount.BB_Co_male);
	  //thismafa=(double)(oneCount.AA_Ca_male+oneCount.AA_Co_male)/(oneCount.AA_Ca_male+oneCount.AA_Co_male+oneCount.BB_Ca_male+oneCount.BB_Co_male+oneCount.OO_Co_male+oneCount.OO_Ca_male);
	  thismafa=1/(double)Nalleles*floor((double)(oneCount.AA_Ca_male+oneCount.AA_Co_male)*Nalleles/(oneCount.AA_Ca_male+oneCount.AA_Co_male+oneCount.BB_Ca_male+oneCount.BB_Co_male-oneCount.OO_Co_male-oneCount.OO_Ca_male));

	  if(!qt)
	    {
	      thiscontrolmaf= (double)(oneCount.AA_Co_male)/(oneCount.AA_Co_male+oneCount.BB_Co_male);
	      thiscasemaf= (double)(oneCount.AA_Ca_male)/(oneCount.AA_Ca_male+oneCount.BB_Ca_male);

	      // thiscontrolmafa= (double)(oneCount.AA_Co_male)/(oneCount.AA_Co_male+oneCount.BB_Co_male+oneCount.OO_Co_male);
	      // thiscasemafa=(double)(oneCount.AA_Ca_male)/(oneCount.AA_Ca_male+oneCount.BB_Ca_male+oneCount.OO_Ca_male);

	      thiscontrolmafa=1/(double)Nalleles_Co*floor((double)oneCount.AA_Co_male*Nalleles_Co/(oneCount.AA_Co_male+oneCount.BB_Co_male-oneCount.OO_Co_male));
	      thiscasemafa=1/(double)Nalleles_Ca*floor((double)oneCount.AA_Ca_male*Nalleles_Ca/(oneCount.AA_Ca_male+oneCount.BB_Ca_male-oneCount.OO_Ca_male));

	    }
	  else
	    {
	      thiscontrolmaf=thismaf;
	      thiscasemaf=thismaf;

	      thiscontrolmafa=thismafa;
	      thiscasemafa=thismafa;

	    }
	} //MAF
      else if (strcmp(map[j].chr, "26") == 0)
	{
	  //MAF
	  int Nalleles=(2*oneCount.AA_Ca+2*oneCount.AB_Ca+2*oneCount.BB_Ca+2*oneCount.AA_Co+2*oneCount.AB_Co+2*oneCount.BB_Co+2*oneCount.OO_Co+2*oneCount.OO_Ca);
	  int Nalleles_Co=(2*oneCount.AA_Co+2*oneCount.AB_Co+2*oneCount.BB_Co+2*oneCount.OO_Co);
	  int Nalleles_Ca=(2*oneCount.AA_Ca+2*oneCount.AB_Ca+2*oneCount.BB_Ca+2*oneCount.OO_Ca);

	  thismaf= (double)(2*oneCount.AA_Ca+oneCount.AB_Ca+2*oneCount.AA_Co+oneCount.AB_Co)/(2*oneCount.AA_Ca+2*oneCount.AB_Ca+2*oneCount.BB_Ca+2*oneCount.AA_Co+2*oneCount.AB_Co+2*oneCount.BB_Co);

	  // thismafa= (double)(2*oneCount.AA_Ca+oneCount.AB_Ca+2*oneCount.AA_Co+oneCount.AB_Co)/(2*oneCount.AA_Ca+2*oneCount.AB_Ca+2*oneCount.BB_Ca+2*oneCount.AA_Co+2*oneCount.AB_Co+2*oneCount.BB_Co+2*oneCount.OO_Co+2*oneCount.OO_Ca);
	  thismafa= 1/(double)Nalleles*floor((double)(2*oneCount.AA_Ca+oneCount.AB_Ca+2*oneCount.AA_Co+oneCount.AB_Co)*Nalleles/(2*oneCount.AA_Ca+2*oneCount.AB_Ca+2*oneCount.BB_Ca+2*oneCount.AA_Co+2*oneCount.AB_Co+2*oneCount.BB_Co-2*oneCount.OO_Co-2*oneCount.OO_Ca));
	  if(!qt)
	    {
	      thiscontrolmaf=
		(double)(2*oneCount.AA_Co+oneCount.AB_Co)/(2*oneCount.AA_Co+2*oneCount.AB_Co+2*oneCount.BB_Co);
	      thiscasemaf=
		(double)(2*oneCount.AA_Ca+oneCount.AB_Ca)/(2*oneCount.AA_Ca+2*oneCount.AB_Ca+2*oneCount.BB_Ca);

	      //			   thiscontrolmafa=
	      //			   (double)(2*oneCount.AA_Co+oneCount.AB_Co)/(2*oneCount.AA_Co+2*oneCount.AB_Co+2*oneCount.BB_Co+2*oneCount.OO_Co);
	      //			   thiscasemafa=
	      //			   (double)(2*oneCount.AA_Ca+oneCount.AB_Ca)/(2*oneCount.AA_Ca+2*oneCount.AB_Ca+2*oneCount.BB_Ca+2*oneCount.OO_Ca);

	      thiscontrolmafa=
		1/(double)Nalleles_Co*floor((double)(2*oneCount.AA_Co+oneCount.AB_Co)*Nalleles_Co/(2*oneCount.AA_Co+2*oneCount.AB_Co+2*oneCount.BB_Co-2*oneCount.OO_Co));
	      thiscasemafa=
		1/(double)Nalleles_Ca*floor((double)(2*oneCount.AA_Ca+oneCount.AB_Ca)*Nalleles_Ca/(2*oneCount.AA_Ca+2*oneCount.AB_Ca+2*oneCount.BB_Ca-2*oneCount.OO_Ca));


	    }
	  else
	    {
	      thiscontrolmaf=thismaf;
	      thiscasemaf=thismaf;

	      thiscontrolmafa=thismafa;
	      thiscasemafa=thismafa;
	    }
	}

      /// MAF
      map[j].maf=thismaf;
      map[j].controlmaf=thiscontrolmaf;
      map[j].casemaf=thiscasemaf;
#if RARE
      map[j].mafa=thismafa;
      map[j].controlmafa=thiscontrolmafa;
      map[j].casemafa=thiscasemafa;
      if(mafadjust==1)
	{
	  map[j].mafr=rounddown(thismafa,7);
	  map[j].controlmafr=rounddown(thiscontrolmafa,7);
	  map[j].casemafr=rounddown(thiscasemafa,7);
	}
      if(mafadjust==0)
	{

	  map[j].mafr= rounddown(thismaf,7);
	  map[j].controlmafr= rounddown(thiscontrolmaf,7);
	  map[j].casemafr= rounddown(thiscasemaf,7);
	}
#endif
      if(map[j].maf<maf)
	{
	  map[j].qcin = 0;
	}

      /// qcin
      if (!map[j].qcin) {
	map[j].analysis_in = 0;
	map[j].matching_in = 0;
      }

      /// set counters for snps
      if (map[j].qcin == 1)
	{
	  nsnpqc++;
	}
      if (map[j].analysis_in == 1)
	{
	  nsnpAnalysisqc++;
	}
      if (map[j].matching_in == 1)
	{
	  nsnpMatchingqc++;
	}
      if(map[j].analysis_in == 1 && map[j].matching_in == 1)
	{
	  nsnpAnalysisMatchingqc++;
	}
      if (!strcmp(map[j].chr, "23") || !strcmp(map[j].chr, "24"))
	{
	  with25=1;
	}
    }

    if (!nsnpAnalysisqc)               die("No SNPs are left for Analysis after QC!");
    if (!nsnpMatchingqc && stratify) error("No SNPs are left for Stratifying after QC!");
    if (with25 && (!sexcov || singleMarkerTest <3))	logg("WARNING! Chromosomes 23 & 24 should be analyzed using regression with sex as a covariate!\n");

    /// set counters for ppls
    for (i = 0; i < nlinestfam; i++) {
      switch (person[i].aff[thread]) {
      case 1:  ncontrolsqc++; break;
      case 2:  ncasesqc++;    break;
      default: nrestqc++;
      }
    }
    npplqc = ncasesqc + ncontrolsqc;
    ncases=0; ncontrols=0; nrest=0;

    if (qt) sstm << "individuals: " << npplqc << ", notUsed: " << nrestqc;
    else    sstm << "cases: " << ncasesqc << ", controls: " << ncontrolsqc << ", notUsed: " << nrestqc;
    sstm << "\nsnps: " << nsnpqc << ", " << "analysis-snps: " << nsnpAnalysisqc;
    if (stratify) sstm << ", " << "matching-snps: " << nsnpMatchingqc << ", " << "analysis-matching-snps: " << nsnpAnalysisMatchingqc;
    logg(sstm);

    cutoff=((double)singletop)/((double)nsnpqc);
    if(pathwayAnalysis && p_PAA != -1)
      {
	singletop=p_PAA*((double)nsnpqc);if(singletop<1){singletop=1;}
	logfile << "SINGLETOP was adjusted to " << singletop << " because of option P_PAA\n";
	cutoff=p_PAA;
      }


    /// remove SNPs from binary coding and write removed to file
    deletedSnps.open(delSnps.c_str(), ios::out);
    deletedSnps << "rs_No\tRow_No\n";
    if (nsnpqc < nlinestped) {
      uint64_t*** oldBinSNPs = BinSNPs;
      BinSNPs = new uint64_t**[nsnpqc];
      for (int i=0, j=0; j<nlinestped; j++) {
	if (map[j].qcin) {
	  BinSNPs[i++] = oldBinSNPs[j];
	} else {
	  delete_2dim(oldBinSNPs[j]);
	  deletedSnps << map[j].rs << "\t" << map[j].line << "\n";
	}
      }
      delete[] oldBinSNPs; oldBinSNPs = NULL;
    }
    deletedSnps.close();

    /// remove persons from binary coding and write removed to file

    deletedPerson.open(delPers.c_str(), ios::out);
    deletedPerson << "FID\tPID\n";
    if (npplqc < (uint32_t)nlinestfam) {
      BinSNPsQCinFlags = getQCIn_bin_snp_gen(person, nlinestfam, nwordsSNPs);
      nwordsSNPs = reduce_bin_gen(BinSNPs, nsnpqc, nwordsSNPs, BinSNPsQCinFlags);
      delete[] BinSNPsQCinFlags; BinSNPsQCinFlags=NULL;
      for (int i = 0; i < nlinestfam; i++) {
	if (!person[i].qcin) {
	  deletedPerson << person[i].fid << "\t" << person[i].pid << "\n";
	}
      }
    }
    deletedPerson.close();

    delete_2dim(BinSNPsCCFlags);
    delete_2dim(BinSNPsGenderFlags);
    BinSNPsCCFlags     =     getCC_bin_snp_gen(person, nlinestfam, nwordsSNPs);
    BinSNPsGenderFlags = getGender_bin_snp_gen(person, nlinestfam, nwordsSNPs);
  }

  /***** IBS determination *****/
  if (doIbs) {
    if (BinPPLs == NULL) {
      logg("transposing binary coding to people major-mode ...");
      nwordsPPLs = transpose_bin_gen(BinSNPs, nsnpqc, npplqc, BinPPLs);
    }
    logg("\nPerforming IBS-analysis ...");

    struct IBS {
      double mean;
      double n;
      double squaresum;
      double sd;
    };
    fstream ibsRelativesFile;
    fstream ibsOutlierFile;
    string ibsRelatives  = outputname + "ibsRelatives.txt";
    string ibsOutlier    = outputname + "ibsOutlier.txt";
    struct IBS ibs = {0,0,0,0};
    struct IBS lowIbs = {0,0,0,0};
    uint32_t nRelatives=0;
    uint32_t nOutlier=0;
    float ibsHelp=0;

    struct RELATIVES *relatives = NULL;     // list of relative pairs
    int *outlier = NULL;                    // list of outliers
    int* lowIbsCounts = new int[npplqc]; // nr of low IBS per person
    memset(lowIbsCounts, 0, npplqc*sizeof(int));

    /** get IBS-matrix **/
    float** IbsM = getIbsMatrix(BinPPLs, npplqc, nwordsPPLs);
#if DEV
    printmatrix(IbsM, npplqc, outputname+"mibs.txt");
#endif

    int* People = new int[npplqc];                 // saves line-nr of tfam-file
    for (int i=0,k=0; k < nlinestfam; k++) {
      if (person[k].qcin) People[i++] = k;
    }

    /** calculate mean IBS **/
    for (uint32_t i = 0; i < npplqc; i++) {
      for (uint32_t j = i+1; j < npplqc; j++) {
	ibsHelp = IbsM[i][j];
	if (0 <= ibsHelp && ibsHelp <= 1) {
	  ibs.n++;
	  ibs.mean+=ibsHelp;
	  ibs.squaresum+=ibsHelp*ibsHelp;
	}
      }	 // j
    } // i
    ibs.mean=ibs.mean/ibs.n;
    ibs.sd=sqrt((ibs.squaresum-ibs.n*ibs.mean*ibs.mean)/(ibs.n-1));

    /** determine high/low ibs **/
    for (uint32_t i = 0; i < npplqc; i++) {
      for (uint32_t j = i+1; j < npplqc; j++) {
	ibsHelp = IbsM[i][j];
	if (fabs(ibsHelp-ibs.mean) > ibs_SD_relatives*ibs.sd) {
	  if (ibsHelp>ibs.mean) {
	    nRelatives++;
	    relatives = (struct RELATIVES *) realloc(relatives, nRelatives*sizeof(struct RELATIVES));
	    relatives[nRelatives-1].i = i;
	    relatives[nRelatives-1].j = j;
	  } else {
	    lowIbsCounts[i]++;
	    lowIbsCounts[j]++;
	  }
	}
      } // j
    } // i

    struct RELATIVES xr, yr;
    if (nRelatives>1) qsortrelatives(relatives, 0, nRelatives-1, xr, yr, IbsM);

    ibsRelativesFile.open(ibsRelatives.c_str(), ios::out);
    for (uint32_t i,j,k=0; k < nRelatives; k++) {
      i = People[relatives[k].i];
      j = People[relatives[k].j];
      ibsHelp = IbsM[relatives[k].i][relatives[k].j];
      ibsRelativesFile << setprecision(4);
      if (person[i].missing > person[j].missing)
	ibsRelativesFile << person[i].fid << "\t" << person[i].pid << "\t" << person[j].fid << "\t" << person[j].pid;
      else
	ibsRelativesFile << person[j].fid << "\t" << person[j].pid << "\t" << person[i].fid << "\t" << person[i].pid;
      ibsRelativesFile << "\t" << ibsHelp << "\t" << (ibsHelp-ibs.mean)/ibs.sd << "\n";
    }
    ibsRelativesFile.close();
    sstm << "number of found relative-pairs " << nRelatives << " (meanIBS " << ibs.mean << ", sdIBS " << ibs.sd << ", threshold " << ibs_SD_relatives << "sd)";
    logg(sstm);

    /** determine outlier **/
    for (uint32_t i = 0; i < npplqc; i++) {
      lowIbs.n++;
      lowIbs.mean+=lowIbsCounts[i];
      lowIbs.squaresum+=lowIbsCounts[i]*lowIbsCounts[i];
    }
    lowIbs.mean=lowIbs.mean/lowIbs.n;
    lowIbs.sd=sqrt((lowIbs.squaresum-lowIbs.n*lowIbs.mean*lowIbs.mean)/(lowIbs.n-1));
    ibsOutlierFile.open(ibsOutlier.c_str(), ios::out);
    for (uint32_t i = 0; i < npplqc; i++) {
      if( (lowIbsCounts[i]-lowIbs.mean) > ibs_SD_outlier*lowIbs.sd) {
	nOutlier++;
	ibsOutlierFile << person[People[i]].fid << "\t" << person[People[i]].pid << "\t" << lowIbsCounts[i] << "\t" << (lowIbsCounts[i]-lowIbs.mean)/lowIbs.sd << "\n";
      }
    }
    ibsOutlierFile.close();
    sstm << "number of found outliers " << nOutlier << " (meanLowIBS " << lowIbs.mean << ", sdLowIBS " << lowIbs.sd << ", threshold " << ibs_SD_outlier << "sd)";
    logg(sstm);

    free(outlier); outlier=NULL;
    free(relatives); relatives=NULL;
    delete[] lowIbsCounts; lowIbsCounts=NULL;
    delete[] People; People=NULL;
    delete_2dim(IbsM);

    if (doIbs < 2) {
      errorfile.close();
      logfile.close();
      exit(0);
    }
  }


  if ((mWithSingletop == 2) && (singletop < 2))
    {
      singletop = 2;
    }
  else if ((mWithSingletop == 3) && (singletop < 3))
    {
      singletop = 3;
    }
  else if (singletop < 1)
    {
      singletop = 10;
    }

  // Anzahl der SNPs ist kleiner singletop
  if (nsnpAnalysisqc < singletop)
    {
      cout << "Singletop was reduced: old singletop > number of SNPs after QC\n";
      logfile << "Singletop was reduced: old singletop > number of SNPs after QC\n";
      singletop = nsnpAnalysisqc;
    }

  j=0;
  // Erstellen der Geneticlist
  if(nlinesinfo>0 && (mWithGenetictop != 0 && (markercombi2 >=1 || markercombi3 >=1 )) )
    {
      geneticlist = (struct GENETICLIST *) realloc(geneticlist, (nlinestped) * sizeof(struct GENETICLIST));
      if (!geneticlist)
	{
	  errorfile << "memory allocation error in geneticlist\n";
	  logfile << "memory allocation error in geneticlist\n";
	  cout << "memory allocation error in geneticlist\n";
	  //cout << nlinestped << "\n";
	  errorfile.close();logfile.close();exit(1);
	}
      for (i=0; i< nlinestped; i++)
	{
	  geneticlist[j].nr = 99;
	  if (map[i].qcin==1 && ((geneticImpact == 4 && map[i].codingStatus == 3) || (geneticImpact == 3 && map[i].location == 1)
				 || (geneticImpact == 2 && map[i].locationToGene == 0) || (geneticImpact == 1 && (abs(map[i].locationToGene) < 100000 || map[i].location == 2))))
	    {
	      geneticlist[j].nr = i;
	      j++;
	    }
	}

      genetictop = j;
      if (genetictop==0)
	{
	  errorfile << "No SNPs found that meet GENETIC_IMPACT criterion\n";
	  logfile << "No SNPs found that meet GENETIC_IMPACT criterion\n";
	  cout << "No SNPs found that meet GENETIC_IMPACT criterion\n";
	  errorfile.close();logfile.close();exit(1);
	}
      geneticlist = (struct GENETICLIST *) realloc(geneticlist, (genetictop) * sizeof(struct GENETICLIST));
      if (!geneticlist)
	{
	  errorfile << "memory allocation error in geneticlist\n";
	  logfile << "memory allocation error in geneticlist\n";
	  cout << "memory allocation error in geneticlist\n";
	  errorfile.close();logfile.close();exit(1);
	}
    }
  genetictop = j;


  // combilist
  if (combilist == 0)
    {
      qstart = 0;
      nlinescombi = 0;
    }
  else
    {
      qstart = 1;
      printtop = nlinescombi;
      mWithSingletop = 0;
      mWithGenetictop = 0;
    }
  countsP=(double *)realloc(countsP,printtop*sizeof(double));

  if (!countsP)
    {
      errorfile << "memory allocation error in countsP\n";
      logfile << "memory allocation error in countsP\n";
      cout << "memory allocation error in countsP\n";
      errorfile.close();logfile.close();exit(1);
    }

  for (i = 0; i < printtop; i++)
    {
      countsP[i] = 0;
    }

  uint32_t* countsMC = (uint32_t*)calloc(singletop, sizeof(uint32_t));
  if (!countsMC) die("memory allocation error in countsMC");

  if (pathwayImpact == 0 && pathwayAnalysis == 0)
    {
      rstart = 0;
      nlist = 0;
    }
  else
    {
      rstart = 1;
    }

  //PAA
  if (pathwayAnalysis==1)
    {
      Ttable = calloc2Ddouble(nlist+1,nsim+2);
      hTable = calloc2Ddouble(2,nsim+2);
    }

  /// betrifft singlemarker, Anzahl der SNPs ist kleiner printtop
  if ( !pathwayImpact && !markercombi2 && !markercombi3 && nsnpAnalysisqc<printtop )	{
    logg("Printtop was reduced: old printtop > number of SNPs after QC.");
    printtop = nsnpAnalysisqc;
  }


  /***** Case-Control IBS-stratifying *****/
  {
    if (stratify && stratify!=4 && (ncasesqc < 1 || ncontrolsqc < 1)) {
      stratify = 0;
      logg("\nSwitched off case-control IBS-stratifying due to lack of cases/controls.");
    }
    if (stratify) {
      logg("");
      if (BinPPLs == NULL) {
        logg("transposing binary coding to people major-mode ...");
        nwordsPPLs = transpose_bin_gen(BinSNPs, nsnpqc, npplqc, BinPPLs);
      }
      sstm << "Performing IBS-stratifying mode " << stratify << " ...";
      logg(sstm);
      switch (stratify) {
      case 1:
	nMatches  = call_matching_ibs_pairs(BinPPLs, nwordsPPLs, person, nlinestfam, ncontrolsqc, ncasesqc, map, 0, nlinestped, MatchedPairs);
	break;
      case 2:
	nClusters = call_matching_ibs_groups1(BinPPLs, nwordsPPLs, person, nlinestfam, ncontrolsqc, ncasesqc, map, 0, nlinestped, Clusters);
	break;
      case 3:
	nClusters = call_matching_ibs_groups2(BinPPLs, nwordsPPLs, person, nlinestfam, ncontrolsqc, ncasesqc, map, 0, nlinestped, Clusters);
	break;
      case 4:
	nClusters = call_clustering(BinPPLs, nwordsPPLs, person, nlinestfam, ncontrolsqc, ncasesqc, map, 0, nlinestped, Clusters);
	break;
      case 5:
	nMatches  = call_matching_cluster_pairs(BinPPLs, nwordsPPLs, person, nlinestfam, ncontrolsqc, ncasesqc, map, 0, nlinestped, MatchedPairs);
	break;
      case 6:
	nClusters = call_matching_cluster_groups1(BinPPLs, nwordsPPLs, person, nlinestfam, ncontrolsqc, ncasesqc, map, 0, nlinestped, Clusters);
	break;
      case 7:
	nClusters = call_matching_cluster_groups2(BinPPLs, nwordsPPLs, person, nlinestfam, ncontrolsqc, ncasesqc, map, 0, nlinestped, Clusters);
	break;
      case 8:
	nMatches  = call_matching_vicinity_pairs(BinPPLs, nwordsPPLs, person, nlinestfam, ncontrolsqc, ncasesqc, map, 0, nlinestped, MatchedPairs);
	break;
      case 9:
	nClusters = call_matching_vicinity_groups1(BinPPLs, nwordsPPLs, person, nlinestfam, ncontrolsqc, ncasesqc, map, 0, nlinestped, Clusters);
	break;
      case 10:
	nClusters = call_matching_vicinity_groups2(BinPPLs, nwordsPPLs, person, nlinestfam, ncontrolsqc, ncasesqc, map, 0, nlinestped, Clusters);
	break;
      default:
	die("Invalid STRATIFY value!");
      }
      uint32_t oldpplqc = npplqc;
      switch (stratify) {
      case 1: case 5: case 8:
	for (int i=0; i<nlinestfam; i++) {
	  person[i].aff[thread] = 0;
	  person[i].analysis_in = false;
	}
	for (uint32_t i=0; i<nMatches; i++) {
	  person[MatchedPairs[i].i].aff[thread] = 1;
	  person[MatchedPairs[i].i].analysis_in = true;
	  person[MatchedPairs[i].j].aff[thread] = 2;
	  person[MatchedPairs[i].j].analysis_in = true;
	}
	ncontrolsqc = nMatches;
	ncasesqc    = nMatches;
	npplqc = 2*nMatches;
	nrestqc = nlinestfam - npplqc;
	//  string clusterFilename = outputname + "cluster.txt";
	//  fstream clusterFile;
	//	clusterFile.open(clusterFilename.c_str(), ios::out);
	//	for (uint32_t i=0; i<nMatches; i++) {
	//	  clusterFile << i << ":\t" << person[MatchedPairs[i].i].fid << " " << person[MatchedPairs[i].i].pid << "\t" << person[MatchedPairs[i].j].fid << " " << person[MatchedPairs[i].j].pid << endl;
	//	}
	//  clusterFile.close();
	break;
      case 2: case 3: case 4: case 6: case 7: case 9: case 10:
	for (int i=0; i<nlinestfam; i++) {
	  if (person[i].clusterAffil == -1) {
	    person[i].analysis_in = false;
	    switch (person[i].aff[thread]) {
	    case 1: ncontrolsqc--; break;
	    case 2:    ncasesqc--; break;
	    }
	    person[i].aff[thread] = 0;
	  }
	}
	nrestqc = nlinestfam-ncontrolsqc-ncasesqc;
	npplqc = ncontrolsqc + ncasesqc;
	//  string clusterFilename = outputname + "cluster.txt";
	//  fstream clusterFile;
	//	clusterFile.open(clusterFilename.c_str(), ios::out);
	//	for (uint32_t i=0; i<nClusters; i++) {
	//	  clusterFile << i << ":\t";
	//	  for (uint32_t j=0; j<Clusters[i].nppls; j++) {
	//	    clusterFile << person[Clusters[i].list[j]].fid << " " << person[Clusters[i].list[j]].pid << "\t";
	//	  }
	//	  clusterFile << endl;
	//	}
	//  clusterFile.close();
	break;
      }

      /// reduce sample
      if (oldpplqc > npplqc) {
	delete_3dim(BinPPLs, oldpplqc);
	uint64_t* AnalysisFlags = getAnalysisIn_bin_snp_gen(person, nlinestfam, nwordsSNPs);
	nwordsSNPs = reduce_bin_gen(BinSNPs, nsnpqc, nwordsSNPs, AnalysisFlags);
	delete[] AnalysisFlags; AnalysisFlags=NULL;
	for (int i=0; i<nlinestfam; i++) person[i].qcin = person[i].analysis_in;
	if (BinSNPsCCFlags    !=NULL) delete_2dim(BinSNPsCCFlags);
	if (BinSNPsGenderFlags!=NULL) delete_2dim(BinSNPsGenderFlags);
	BinSNPsCCFlags     =     getCC_bin_snp_gen(person, nlinestfam, nwordsSNPs);
	BinSNPsGenderFlags = getGender_bin_snp_gen(person, nlinestfam, nwordsSNPs);
	for (int j=0,i=0; i<nlinestped; i++) {
	  if (map[i].qcin) {
	    calcMaf(&map[i], BinSNPs[j], BinSNPsCCFlags, BinSNPsGenderFlags, nwordsSNPs);
	    calcMR(&map[i], BinSNPs[j], BinSNPsCCFlags, BinSNPsGenderFlags, nwordsSNPs);
	    j++;
	  }
	}
      }

    }  // if (stratify)

    if (familyCluster) {
      if (!familyData) {
	char* famID;
	nClusters = 0;
	for (int n=0; n<nlinestfam; n++) {
	  if (!person[n].qcin) { person[n].clusterAffil=-1; continue; }
	  famID = NULL;
	  for (int o=0; o<n; o++) {
	    if (!person[o].qcin) continue;
	    if (!strcmp(person[n].fid, person[o].fid)) {
	      person[n].clusterAffil = person[o].clusterAffil;
	      famID = person[o].fid;
	      break;
	    }
	  }
	  if (famID==NULL) {
	    person[n].clusterAffil = nClusters;
	    nClusters++;
	  }
	}
      }
      Clusters = Affil2Cluster(nClusters, person, nlinestfam);
      stratify = 4;
      sstm << "\nSetting number of family clusters: " << nClusters;
      logg(sstm);
    }  // if (familyCluster)

    if (cluster_covars) {
      logg("Setting cluster affiliation as indicator covariates.");
      for (int i=0; i<nlinestfam; i++) {
	if (person[i].qcin) {
	  person[i].covin = (unsigned char*)realloc(person[i].covin, (maxIndexCov+nClusters)*sizeof(unsigned char));
	  if (!person[i].covin) die("Memory re-allocation error in person[i].covin!");
	  memset(person[i].covin+maxIndexCov, 1, nClusters*sizeof(unsigned char));

	  person[i].cov = (double*)realloc(person[i].cov, (maxIndexCov+nClusters)*sizeof(double));
	  if (!person[i].cov) die("Memory re-allocation error in person[i].cov!");
	  memset(person[i].cov+maxIndexCov, 0, nClusters*sizeof(double));
	  person[i].cov[maxIndexCov+person[i].clusterAffil] = 1;
	}
      }
      cov = (int*)realloc(cov, (maxIndexCov+nClusters)*sizeof(int));
      if (!cov) die("Memory allocation error in cov!");
      for (int i=0; i<nClusters; i++) cov[maxIndexCov+i] = 1;

      maxIndexCov += nClusters;
      ncov += nClusters;
    }  // if (cluster_covar)

    if (stratify) {
      if (qt) sstm << "individuals: " << npplqc << ", notUsed: " << nrestqc;
      else    sstm << "cases: " << ncasesqc << ", controls: " << ncontrolsqc << ", notUsed: " << nrestqc;
      logg(sstm);
    }

    if (allAllMatching) {
      logg("");
      if (BinPPLs == NULL) {
        logg("transposing binary coding to people major-mode ...");
        nwordsPPLs = transpose_bin_gen(BinSNPs, nsnpqc, npplqc, BinPPLs);
      }
	  string clusterFilename = outputname + "cluster.txt";
      logg("Performing all-all-matching, filename [ " + clusterFilename + " ] ...");
      uint32_t nClusters = call_all_all_matching(BinPPLs, nwordsPPLs, person, nlinestfam, npplqc, map, 0, nlinestped, Clusters);
	  fstream clusterFile;
      clusterFile.open(clusterFilename.c_str(), ios::out);
      for (uint32_t i=0; i<nClusters; i++) {
		clusterFile << i << "\t" << Clusters[i].nppls;
		for (uint32_t j=0; j<Clusters[i].nppls; j++) {
		  clusterFile << "\t" << person[Clusters[i].list[j]].fid << " " << person[Clusters[i].list[j]].pid << " " << ( qt ? person[Clusters[i].list[j]].qtaff[0] : person[Clusters[i].list[j]].aff[0] );
		}
		clusterFile << endl;
      }
	  clusterFile.close();
    }

  }


  /** classic <-> binary mappings **/
  {
    SNPMap        = new int[nsnpqc];      /// maps QCin position to line in MAP-file
    SNPMapInverse = new int[nlinestped];  /// maps line in MAP-file to QCin position
    memset(SNPMapInverse, -1, nlinestped*sizeof(int));
    for (int j=0,i=0; i<nlinestped; i++) {
      if (map[i].qcin) {
        SNPMap[j] = i;
        SNPMapInverse[i] = j;
        j++;
      }
    }

    PPLMap        = new int[npplqc];      /// maps QCin position to line in FAM-file
    PPLMapInverse = new int[nlinestfam];  /// maps line in FAM-file to QCin position
    memset(PPLMapInverse, -1, nlinestfam*sizeof(int));
    for (int j=0,i=0; i<nlinestfam; i++) {
      if (person[i].qcin) {
        PPLMap[j] = i;
        PPLMapInverse[i] = j;
	    j++;
      }
    }

    PplLocations = getPplLocations(person, nlinestfam, npplqc);             /// peoples' coordinates in the binary coding
    PplLocFam    = getPplLocations(person, nlinestfam, nlinestfam, false);  /// peoples' coordinates in the binary coding from fam-list
  }


  /** create counts **/
  counts = (struct COUNTS **)malloc(maxthreads * sizeof(struct COUNTS *));
  if (!counts) die("memory allocation error in counts.");
  for (uint16_t k=0; k<maxthreads; k++) {
    counts[k] = (struct COUNTS *)calloc(nlinestped, sizeof(struct COUNTS));
    if (!counts[k]) die("memory allocation error in counts[k].");
  }
  /** initialize counts **/
  for (int k=0,j=0; j<nlinestped; j++) {
    initCounts(&counts[thread][j], -1, regression);
    if (map[j].qcin) {
      updateCounts(&counts[thread][j], nwordsSNPs, BinSNPs[k], BinSNPsCCFlags, BinSNPsGenderFlags, !(strcmp(map[j].chr, "23") && strcmp(map[j].chr, "24")));
      k++;
    }
  }

  int ncasesqcMale=0;
  int ncasesqcFemale=0;
  int ncontrolsqcMale=0;
  int ncontrolsqcFemale=0;

  if(!qt)
    {
      for (int i=0; i<nlinestfam; i++)
	{
	  if (person[i].qcin==1)
	    {
	      if (person[i].aff[thread]==2 && person[i].sex==1){ncasesqcMale++;}
	      else if (person[i].aff[thread]==2 && person[i].sex==2){ncasesqcFemale++;}
	      else if (person[i].aff[thread]==1 && person[i].sex==1){ncontrolsqcMale++;}
	      else if (person[i].aff[thread]==1 && person[i].sex==2){ncontrolsqcFemale++;}
	    }
	}
      cout << "\ncasesMale: " << ncasesqcMale << ", casesFemale: " << ncasesqcFemale << endl;
      cout << "controlsMale: " << ncontrolsqcMale << ", controlsFemale: " << ncontrolsqcFemale << endl;
    }


#if RARE
  if(bin){
    if((collinter==3 || collinter==4) && optimalrare){
      cout<<"COLL_INTER 3 and 4 is not compatible with the Variable Threshold (VT) method!"<<endl;
      logfile<<"COLL_INTER 3 and 4 is not compatible with the Variable Threshold (VT) method!"<<endl;
      errorfile<<"COLL_INTER 3 and 4 is not compatible with the Variable Threshold (VT) method!"<<endl;
      logfile.close(); errorfile.close(); exit(1);
    }
    if((collinter==1 || collinter==2) && (FISHERflag || CMATflag || COLLREGflag || FRACREGflag || REGRESSIONflag)){
      cout<<"COLL_INTER is not compatible with RARE_TESTS!"<<endl;
      logfile<<"COLL_INTER is not compatible with RARE_TESTS!"<<endl;
      errorfile<<"COLL_INTER is not compatible with RARE_TESTS!"<<endl;
      logfile.close(); errorfile.close(); exit(1);
    }
    else if((collinter==3 || collinter==4) && (FISHERflag || CMATflag || COLLflag || FRACREGflag || REGRESSIONflag)){
      cout<<"COLL_INTER 3 and 4 is not compatible with RARE_TESTS!"<<endl;
      logfile<<"COLL_INTER 3 and 4 is not compatible with RARE_TESTS!"<<endl;
      errorfile<<"COLL_INTER 3 and 4 is not compatible with RARE_TESTS!"<<endl;
      logfile.close(); errorfile.close(); exit(1);
    }
    if(nsim>0 && (collinter==3 || collinter==4)){
      cout   <<"COLL_INTER 3 and 4 does not work with SIMULATIONS>0."<<endl;
      logfile   <<"COLL_INTER 3 and 4 does not work with SIMULATIONS>0."<<endl;
      errorfile   <<"COLL_INTER 3 and 4 does not work with SIMULATIONS>0."<<endl;
      logfile.close();errorfile.close();exit(1);
    }
    // INTERSNPRare
    // 1) count rare SNPs nRareSNPs[nwindows]
    // 2) write mafs of rare SNPs in arrays mafRareSNPs[nwindows][nRareSNPs[j]]
    // 3) calculate nRareLimits[nwindows] and rareLimits[nwindows][nRarelimits[j]]

    if (intervalfile !=" "){
      mergelines=get_gff_data(intervalfile,flanking, merging, expandIntervals, catIntervals, intervalfile_format, featurecol, gffchr, errorfile, logfile, outputname, intervaleditor, verbose);
    }
    if(intervalfile!=" "){
      featurecol=*pfeaturecol;
    }

    if(vb){
      if(NCT==0){
	cout<<"Conversion from MAFT to Carrier Frequency Threshold...\n";
	logfile<<"Conversion from MAFT to Carrier Frequency Threshold...\n";
	float NCFT=2*raref-raref*raref;
	NCT=int(NCFT*npplqc);
	cout<<"MAF of "<<raref<<" corresponds to a Carrier Frequency of "<<NCFT<<". SNPs with up to "<<NCT<<" carriers will be included.\n";
	logfile<<"MAF of "<<raref<<" corresponds to a Carrier Frequency of "<<NCFT<<". SNPs with up to "<<NCT<<" carriers will be included.\n";
      }
      nCarriers = new int[nlinestped]();
      if(!nCarriers)die("Memory allocation error in nCarriers!");
    }
    if(raref==0 && NCT!=0){
      cout<<"Conversion from Carrier Frequency Threshold to MAFT...\n";
      logfile<<"Conversion from Carrier Frequency Threshold to MAFT...\n";
      raref=1-sqrt(1-float(float(NCT)/float(npplqc)));
      cout<<"MAF of "<<raref<<" corresponds to a Carrier Number of "<<NCT<<". SNPs with up to "<<NCT<<" carriers will be included.\n";
      logfile<<"MAF of "<<raref<<" corresponds to a Carrier Number of "<<NCT<<". SNPs with up to "<<NCT<<" carriers will be included.\n";

    }
    windowsRARE(bin, n, nlinestped, binsizeRare, map, counts, &nwindows, &nchrwindows, thread_nloop, errorfile, logfile, intervalfile, raref, intervalfile_format, mergelines, minRareInBin, maxRareInBin, binamin, binamax, featurecol, SetIDfile,intervaleditor,setid, verbose, vb, NCT, BinSNPs, SNPMapInverse, nwordsSNPs);
    // nwindows exists here

    if(vb && optimalrare){
      rareLimitsNCT=new int*[nwindows];
      if (!rareLimitsNCT) die("Memory allocation error in rareLimitsNCT");
      rareLimitsNCTInverse=new int*[nwindows];
      if (!rareLimitsNCTInverse) die("Memory allocation error in rareLimitsNCTInverse");
    }

    // nwindows exists here
    nwindowssinglebin=nwindows; // in case of no rare interaction analysis

    if(collinter>0){
      if(nwindows==1){
	cout << "Only 1 bin, not enough for rare interaction analysis!"<<endl;
	logfile << "Only 1 bin, not enough for rare interaction analysis!"<<endl;
	errorfile<< "Only 1 bin, not enough for rare interaction analysis!"<<endl;
	exit(1);
      }
      nwindows=nwindows*(nwindows+1)/2;
      doublewindowcoord=calloc2Dint(nwindows,2);
    }

    bins2testt=nwindows;

    if(vb==0){
      rareLimits=new double*[nwindows];
      if (!rareLimits) die("Memory allocation error in rareLimits!");
    }

    windowPositions=calloc2Dint(nwindowssinglebin, 2);

    if(verbose==2){
      nSNPsInWindow=new int[nwindowssinglebin]();
      if (!nSNPsInWindow) die("Memory allocation error in nSNPsInWindow!");
    }
    if(collinter==3 || collinter==4){
      nRareSNPs = new int[nwindowssinglebin]();
      if (!nRareSNPs) die("Memory allocation error in nRareSNPs!");
    }
    else{
      nRareSNPs = new int[nwindows]();
      if (!nRareSNPs) die("Memory allocation error in nRareSNPs!");
    }
    if(collinter==3 || collinter==4){
      nRareLimits = new int[nwindowssinglebin]();
      if (!nRareLimits) die("Memory allocation error in nRareLimits!");
     }
    else{
      nRareLimits = new int[nwindows]();
      if (!nRareLimits) die("Memory allocation error in nRareLimits!");
     }


    if(collinter==3 || collinter==4){
      window=new struct WINDOW[nwindowssinglebin];
      if(!window)die("Memory allocation error in window struct!");
    }
    else{
      window=new struct WINDOW[nwindows];
      if(!window)die("Memory allocation error in window struct!");
    }


    if(!optimalrare && vb==1){
      for(l=0; l<nwindows; l++){
	nRareLimits[l]=1;
      }
    }

    if(vb_bmp){
      outRareVB_BMP(person, map, counts, optimalrare, nsim, rarefileVB_BMP, outputname, errorfile, logfile, rarepretest, rarepretestlimit, rareregpretest,  nSNPsInWindow, raref, thread, window, minRareInBin,  mafadjust, nlinestped, verbose, SetIDfile, intervaleditor, setid, nvbstart, vbstart, nvbend, vbend, vbmaxpermstat, ncasesqc, ncontrolsqc, BinSNPs, BinSNPsCCFlags, nwordsSNPs, SNPMapInverse);
    }

    limitspositionsRARE(optimalrare, n, bin, binsizeRare, nlinestped, nlinestfam, raref,  nsim,  map, counts, nRareSNPs, windowPositions, nSNPsInWindow, rareLimits, rareLimitsNCT,rareLimitsNCTInverse, nRareLimits, nchrwindows, nwindows, nwindowssinglebin, thread_nloop, errorfile, logfile, intervalfile, binamin, binamax, SetIDfile, intervaleditor,setid, window, window_inter, verbose, doublewindowcoord, collinter, nwordsSNPs, SNPMapInverse, BinSNPs, vb, NCT, nCarriers);

    if(vb==1){
      vbmaxpermstat=new float[nsim]();
      if(!vbmaxpermstat)die("Memory allocation error in vbmaxpermstat!");
      nvbstartvt=new int[nwindows]();
      if(!nvbstartvt)die("Memory allocation error in nvbstartvt!");
      vbstartvt=new int*[nwindows];
      if(!vbstartvt)die("Memory allocation error in vbstartvt!");
      for(int l=0; l<nwindows; l++){
	vbstartvt[l] = (int *)calloc(1, sizeof(int));
	if(!vbstartvt[l])die("Memory allocation error in vbstartvt[l]!");
      }
      nvblevel=new int*[nwindows];
      if(!nvblevel)die("Memory allocation error in nvblevel!");
      nchunks=new int**[nwindows];
      if(!nchunks)die("Memory allocation error in nchunks!");
      chunkpos=new int***[nwindows];
      if(!chunkpos)die("Memory allocation error in chunkpos!");
      chunklen=new int***[nwindows];
      if(!chunklen)die("Memory allocation error in chunkpos!");


      BinSNPsCCFlagsOriginal=new uint64_t*[nwordsSNPs];
      if(!BinSNPsCCFlagsOriginal)die("Memory allocation error in BinSNPsCCFlagsOriginal!");
      for (int p=0; p<nwordsSNPs; p++){
	BinSNPsCCFlagsOriginal[p]=new uint64_t[3];
	if(!BinSNPsCCFlagsOriginal[p])die("Memory allocation error in BinSNPsCCFlagsOriginal[p]!");
	for(int i=0; i<3; i++){
	  BinSNPsCCFlagsOriginal[p][i]=BinSNPsCCFlags[p][i];
	  //	  if(!BinSNPsCCFlagsOriginal[p][i])die("Memory allocation error in BinSNPsCCFlagsOriginal[p][i]!");
	}
      }

      ndummyatlevel=new int***[nwindows];
      if(!ndummyatlevel)die("Memory allocation error in ndummyatlevel!");
      dummylevel=new int****[nwindows];
      if(!dummylevel)die("Memory allocation error in dummylevel!");
      dummypos=new int****[nwindows];
      if(!dummypos)die("Memory allocation error in dummypos!");
      dummycluster=new int****[nwindows];
      if(!dummycluster)die("Memory allocation error in dummycluster!");
      dummyend=new  int***[nwindows];
      if(!dummyend)die("Memory allocation error in dummyend!");
      ndummyends=new  int**[nwindows];
      if(!ndummyends)die("Memory allocation error in ndummyends!");
      nchunkcluster=new  int**[nwindows];
      if(!nchunkcluster)die("Memory allocation error in nchunkcluster!");
      BinCarriers = new uint64_t**[nwindows];
      if(!BinCarriers)die("Memory allocation error in BinCarriers!");

      if(vb_binwise_corr){
	vbbinwisestat=new float****[nwindows];
	if(!vbbinwisestat)die("Memory allocation error in vbbinwisestat!");
	vbbinwisecount=new short int****[nwindows];
	if(!vbbinwisecount)die("Memory allocation error in vbbinwisecount!");
      }

      find_distinct_vb_vt(person, map, PPLMap, PplLocations, BinSNPs, BinSNPsCCFlags, nwordsSNPs, SNPMapInverse, window, nRareLimits, rareLimitsNCT, rareLimitsNCTInverse, nwindows, minIndInBin, maxIndInBin, nvbstartvt, vbstartvt, nlinestfam, nlinestped, nCarriers, nvblevel, nchunks, chunkpos, chunklen, optimalrare, ndummyatlevel, dummypos, dummylevel, dummyend, ndummyends, NCT, BinCarriers, nchunkcluster, dummycluster);

      if(vb_binwise_corr){
	for(int l=0; l<nwindows; l++){
	  vbbinwisestat[l]=new float***[nvbstartvt[l]];
	  if(!vbbinwisestat[l])die("Memory allocation error in vbbinwisestat[l]!");
	  vbbinwisecount[l]=new short int***[nvbstartvt[l]];
	  if(!vbbinwisecount[l])die("Memory allocation error in vbbinwisecount[l]!");
	  for(int m=0; m<nvbstartvt[l]; m++){
	    vbbinwisestat[l][m]=new float**[nvblevel[l][m]];
	    if(!vbbinwisestat[l][m])die("Memory allocation error in vbbinwisestat[l][m]!");
	    vbbinwisecount[l][m]=new short int**[nvblevel[l][m]];
	    if(!vbbinwisecount[l][m])die("Memory allocation error in vbbinwisecount[l][m]!");
	    for(int m2=0; m2<nvblevel[l][m]; m2++){
	      vbbinwisestat[l][m][m2]=new float*[nchunks[l][m][m2]];
	      if(!vbbinwisestat[l][m][m2])die("Memory allocation error in vbbinwisestat[l][m][m2]!");
	      vbbinwisecount[l][m][m2]=new short int*[nchunks[l][m][m2]];
	      if(!vbbinwisecount[l][m][m2])die("Memory allocation error in vbbinwisecount[l][m][m2]!");
	      for(int m3=0; m3<nchunks[l][m][m2]; m3++){
		vbbinwisestat[l][m][m2][m3]=new float[chunklen[l][m][m2][m3]]();
		if(!vbbinwisestat[l][m][m2])die("Memory allocation error in vbbinwisestat[l][m][m2][m3]!");
		vbbinwisecount[l][m][m2][m3]=new short int[chunklen[l][m][m2][m3]]();
		if(!vbbinwisecount[l][m][m2][m3])die("Memory allocation error in vbbinwisecount[l][m][m2][m3]!");
	      }
	    }
	  }
	}
      }


      if(vb_bmp_level){
	//  outRareVB_BMP_level(person, map, counts, optimalrare, nsim, rarefileVB_BMP, outputname, errorfile, logfile, rarepretest, rarepretestlimit, rareregpretest,  nSNPsInWindow, raref, thread, window, minRareInBin,  mafadjust, nlinestped, verbose, SetIDfile, intervaleditor, setid, nvbstart, vbstart, nvbend, vbend, vbstat, vbmaxpermstat, ncasesqc, ncontrolsqc, BinSNPs, BinSNPsCCFlags, nwordsSNPs, SNPMapInverse, nvbstartvt, nvblevel, nchunksatlevel, vbcoords);
      }

      for(int l=0; l<nwindows; l++){
	for(int m=0; m<nvbstartvt[l]; m++){
	  //      cout <<"nvblevel["<<l<<"]["<<m<<"] "<<nvblevel[l][m]<<endl;;
	  for(int m3=0; m3<nvblevel[l][m]; m3++){
	    for(int m4=0; m4<nchunks[l][m][m3]; m4++){
	    }
	  }
	  //    delete[] nchunks[l][m];
	}
	//  delete[] nvblevel[l];
      }
      //delete[] nvblevel;
      //delete[] nchunks;

      if(vb_bmp_level){
	cout<<"*bmp of SNP with most vb Carrier Threshold levels generated!"<<endl;
	logfile<<"*bmp of SNP with most vb Carrier Threshold levels generated!"<<endl;
	exit(0);
      }
      //      }
    }


    if(vb==0){
      if(FISHERflag){
	rarefFISHER= new double[nwindows]();
	pFISHERvec = new double[nwindows]();
	pFISHERvecChi = new double[nwindows]();
	FISHERcount=new int[nwindows]();
	completedFISHER=new int[nwindows]();
	FISHERstats=new double[nwindows]();
	FISHERpermstats=new double[nwindows]();
	FISHERnotconverged=new int[nwindows]();
	if(nsim!=0){
	  FISHER_CI=calloc2Ddouble(nwindows,2);
	}
	if(optimalrare){
	  MAF_Level_VT_FISHER=new int[nwindows]();
	}
	if(wilsonpretest!=0 || rareregpretest!=0){
	  FISHERpretestpassed=new int[nwindows]();
	}
      }
      if(REGRESSIONflag){
	rarefREGRESSION= new double[nwindows]();
	pREGRESSIONvec = new double[nwindows]();
	REGRESSIONcount=new int[nwindows]();
	completedREGRESSION=new int[nwindows]();
	REGRESSIONstats=new double[nwindows]();
	REGRESSIONpermstats=new double[nwindows]();
	if(optimalrare){
	  MAF_Level_VT_REGRESSION=new int[nwindows]();
	}
	if(nsim!=0){
	  REGRESSION_CI=calloc2Ddouble(nwindows,2);
	}
	if(wilsonpretest!=0 || rareregpretest!=0){
	  REGRESSIONpretestpassed=new int[nwindows]();
	}
      }
      if(FRACREGflag){
	rarefFRACREG= new double[nwindows]();
	pFRACREGvec= new double[nwindows]();
	FRACREGcount=new int[nwindows]();
	completedFRACREG=new int[nwindows]();
	FRACREGstats=new double[nwindows]();
	FRACREGpermstats=new double[nwindows]();
	FRACREG_beta=new double[nwindows]();
	FRACREG_se=new double[nwindows]();
	if(optimalrare){
	  MAF_Level_VT_FRACREG=new int[nwindows]();
	}
	if(nsim!=0){
	  FRACREG_CI=calloc2Ddouble(nwindows,2);
	}
	if(wilsonpretest!=0 || rareregpretest!=0){
	  FRACREGpretestpassed=new int[nwindows]();
	}
      }
      if(COLLREGflag){
	pCOLLREGvec = new double[nwindows]();
	rarefCOLLREG= new double[nwindows]();
	if(collinter!=3 && collinter!=4){
	  COLLREGcount=new int[nwindows]();
	  completedCOLLREG=new int[nwindows]();
	  COLLREGstats=new double[nwindows]();
	  COLLREGpermstats=new double[nwindows]();
	  COLLREG_beta=new double[nwindows]();
	  COLLREG_se=new double[nwindows]();
	}
	if(optimalrare){
	  MAF_Level_VT_COLLREG=new int[nwindows]();
	}
	if(nsim!=0){
	  COLLREG_CI=calloc2Ddouble(nwindows,2);
	}
	if(wilsonpretest!=0 || rareregpretest!=0){
	  COLLREGpretestpassed=new int[nwindows]();
	}
      }
      if(CMATflag){
	rarefCMAT= new double[nwindows]();
	OR_CMAT_vec=new double[nwindows]();
	pCMATvec= new double[nwindows]();
	CMATcount=new int[nwindows]();
	completedCMAT=new int[nwindows]();
	CMATstats=new double[nwindows]();
	CMATpermstats=new double[nwindows]();
	CMAT_CI=calloc2Ddouble(nwindows,2);
	if(optimalrare){
	  MAF_Level_VT_CMAT=new int[nwindows]();
	}
	if(wilsonpretest!=0 || rareregpretest!=0){
	  CMATpretestpassed=new int[nwindows]();
	}
      }
      if(COLLflag){
	rarefCOLL= new double[nwindows]();
	OR_COLL_vec=new double[nwindows]();
	OR_COLL_f_vec=new double[nwindows]();
	pCOLLvec= new double[nwindows]();
	COLLcount=new int[nwindows]();
	completedCOLL=new int[nwindows]();
	COLLstats=new double[nwindows]();
	COLLpermstats=new double[nwindows]();
	if(optimalrare){
	  MAF_Level_VT_COLL=new int[nwindows]();
	}
	if(nsim!=0){
	  COLL_CI=calloc2Ddouble(nwindows,2);
	}
	if(wilsonpretest!=0 || rareregpretest!=0){
	  COLLpretestpassed=new int[nwindows]();
	}
      }


      // struct WINDOW *window;

      for(l=0; l<nwindows; l++)
	{
	  if(FISHERflag){
	    rarefFISHER[l]=raref;
	    pFISHERvec[l]=0;
	    pFISHERvecChi[l]=0;
	    FISHERnotconverged[l]=0;
	    if(wilsonpretest!=0 || rareregpretest!=0){
	      FISHERpretestpassed[l]=1;
	    }
	    FISHER_CI[l][1]=1;
	  }
	  if(REGRESSIONflag){
	    rarefREGRESSION[l]=raref;
	    pREGRESSIONvec[l]=0;
	    if(wilsonpretest!=0 || rareregpretest!=0){
	      REGRESSIONpretestpassed[l]=1;
	    }
	    if(nsim!=0){
	      REGRESSION_CI[l][1]=1;
	    }
	  }
	  if(FRACREGflag){
	    rarefFRACREG[l]=raref;
	    pFRACREGvec[l]=0;
	    if(wilsonpretest!=0 || rareregpretest!=0){
	      FRACREGpretestpassed[l]=1;
	    }
	    if(nsim!=0){
	      FRACREG_CI[l][1]=1;
	    }
	  }
	  if(COLLREGflag){
	    rarefCOLLREG[l]=raref;
	    pCOLLREGvec[l]=0;
	    if(wilsonpretest!=0 || rareregpretest!=0){
	      COLLREGpretestpassed[l]=1;
	    }
	    if(nsim!=0){
	      COLLREG_CI[l][1]=1;
	    }
	  }
	  if(CMATflag){
	    rarefCMAT[l]=raref;
	    OR_CMAT_vec[l]=0;
	    pCMATvec[l]=0;
	    if(wilsonpretest!=0 || rareregpretest!=0){
	      CMATpretestpassed[l]=1;
	    }
	    CMAT_CI[l][1]=1;
	  }
	  if(COLLflag){
	    rarefCOLL[l]=raref;
	    OR_COLL_vec[l]=0;
	    OR_COLL_f_vec[l]=0;
	    pCOLLvec[l]=0;
	    if(wilsonpretest!=0 || rareregpretest!=0){
	      COLLpretestpassed[l]=1;
	    }
	    if(nsim!=0){
	      COLL_CI[l][1]=1;
	    }
	  }
	  if(collinter!=3 && collinter!=4){
	    window[l].ntests=0;
	    for(int k=0; k<6; k++){
	      window[l].nperm[k]=0;
	    }
	    if(nRareSNPs!=0){
	      if(FISHERflag){
		window[l].ntests+=1;
	      }
	      if(CMATflag){
		window[l].ntests+=1;
	      }
	      if(COLLflag){
		window[l].ntests+=1;
	      }
	      if(REGRESSIONflag){
		window[l].ntests+=1;
	      }
	      if(FRACREGflag){
		window[l].ntests+=1;
	      }
	      if(COLLREGflag){
		window[l].ntests+=1;
	      }
	    }
	  }
	}
      if(!vb && COLLflag){
	limitsStatCOLL=new double*[nwindows];
	for(l=0; l<nwindows; l++){
	  limitsStatCOLL[l]=new double [nRareLimits[l]]();
	}
      }
      if(CMATflag){
	limitsStatCMAT=new double*[nwindows];
	for(l=0; l<nwindows; l++){
	  limitsStatCMAT[l]=new double [nRareLimits[l]]();
	}
      }
      if(FISHERflag){
	limitsStatFISHER=new double*[nwindows];
	for(l=0; l<nwindows; l++){
	  limitsStatFISHER[l]=new double [nRareLimits[l]]();
	}
      }
      if(REGRESSIONflag){
	//	limitsStatREGRESSION=new double*[nwindows];
	for(l=0; l<nwindows; l++){
	  //	  limitsStatREGRESSION[l]=new double [nRareLimits[l]]();
	}
      }
      if(FRACREGflag){
	//	limitsStatFRACREG=new double*[nwindows];
	for(l=0; l<nwindows; l++){
	  //	  limitsStatFRACREG[l]=new double [nRareLimits[l]]();
	}
      }
      if(COLLREGflag){
	//	limitsStatCOLLREG=new double*[nwindows];
	for(l=0; l<nwindows; l++){
	  //	  limitsStatCOLLREG[l]=new double [nRareLimits[l]]();
	}
      }
    }
    /// RARE matching
    if (rare_stratify) {
      if (BinPPLs == NULL) {
	    notice("transposing binary coding to people major-mode.");
	    nwordsPPLs = transpose_bin_gen(BinSNPs, nsnpqc, npplqc, BinPPLs);
      }
      nMWindowsRARE = call_matchingRARE(BinPPLs, nwordsPPLs, SNPMapInverse, PPLMapInverse, person, nlinestfam, ncontrolsqc, ncasesqc, map, nlinestped, chrPositions, rare_strat_maft, rare_strat_min, MatchingRARE);
      string matchingFilename = outputname + "matchingRARE.txt";
      fstream matchingFile;
      matchingFile.open(matchingFilename.c_str(), ios::out);
      for (int w=0; w<nMWindowsRARE; w++) {
	    matchingFile << "========== Window " << w << " ==========" << endl;
	    for (int i=0; i<MatchingRARE[w].nGroups; i++) {
	      matchingFile << i;
	      for (int j=0; j<MatchingRARE[w].groups[i].nppls; j++)
	        matchingFile << "\t" << person[MatchingRARE[w].groups[i].list[j]].pid << " " << person[MatchingRARE[w].groups[i].list[j]].fid;
	      matchingFile << endl;
	    }
      }
      matchingFile.close();
    }

  }	// end bin
#endif


  /// Logistic regression
  {
    // Logistic regression framework introcuded by Cordell and Clayton (2002)
    // Within this framework, it is possible to include or exclude marginal effects,
    // distinguish allelic and genotypic tests and adjust for covariates


    // Singlemarker mit logistischer Regression

    xsinglevec[0] = 1;
    zsinglevec[0] = 1;

    for (i = 1; i < 27; i++)
      {
	xsinglevec[i] = 0;
	zsinglevec[i] = 0;
	xsinglevec2[i] = 0;
	zsinglevec2[i] = 0;
	xsinglevec3[i] = 0;
	zsinglevec3[i] = 0;
      }


    if (singleMarkerTest == 3 || singleMarkerTest >= 5)
      {
	xsinglevec[1] = 1;
      }
    else if (singleMarkerTest == 4)
      {
	xsinglevec[1] = 1;
	xsinglevec[2] = 1;
      }

    if(pathwayAnalysis && pathwayTest==5)
      {
	if ((singleMarkerTest==3 || singleMarkerTest >=5) || singleMarkerTest == 1)
	  {
	    xsinglevec2[1] = 1;
	    xsinglevec2[3] = 1;
	    zsinglevec2[1] = 1;

	    xsinglevec3[1] = 1;
	    xsinglevec3[3] = 1;
	    xsinglevec3[9] = 1;
	    zsinglevec3[1] = 1;
	    zsinglevec3[3] = 1;
	  }
	else if (singleMarkerTest == 4 || singleMarkerTest == 2)
	  {
	    xsinglevec2[1] = 1;
	    xsinglevec2[2] = 1;
	    xsinglevec2[3] = 1;
	    xsinglevec2[4] = 1;
	    zsinglevec2[1] = 1;
	    zsinglevec2[2] = 1;

	    xsinglevec3[1] = 1;
	    xsinglevec3[2] = 1;
	    xsinglevec3[3] = 1;
	    xsinglevec3[4] = 1;
	    xsinglevec3[9] = 1;
	    xsinglevec3[10] = 1;

	    zsinglevec3[1] = 1;
	    zsinglevec3[2] = 1;
	    zsinglevec3[3] = 1;
	    zsinglevec3[4] = 1;
	  }

      }

    //Überprüfung, ob Covariatenauswahl möglich
    for (l=0; l<0;l++)
      {
	if (cov[l] == 1)
	  {
	    allequal = 1;
	    val1 = -999.0;
	    val2 = -999.0;

	    for (k=0; k<nlinestfam;k++)
	      {
		if (person[k].qcin == 1 && person[k].aff[thread] != 0
		    && person[k].covin[l] == 1 && val1 == -999.0)
		  {
		    val1 = person[k].cov[l];
		  }
		else if (person[k].qcin == 1 && person[k].aff[thread] != 0
			 && person[k].covin[l] == 1 && val2 != -999.0)
		  {
		    val2 = person[k].cov[l];
		  }

		if (val1 != -999.0 && val1 != val2)
		  {
		    allequal = 0;
		    break;


		  }
	      } // end k

	    if (allequal == 1)
	      {
		cov[l] = 0;
	      }
	  }
      } // end l


    if (regression || singleMarkerTest>2 || qt==1)
      {
	if (markercombi2 == 1 && markercombi3 == 1)
	  {
	    cout << "Simultaneous 2-marker and 3-marker analysis not possible with regression\n";
	    logfile << "Simultaneous 2-marker and 3-marker analysis not possible with regression\n";
	    errorfile << "Simultaneous 2-marker and 3-marker analysis not possible with regression\n";
	    errorfile.close();logfile.close();exit(1);
	  }
	//	if(qt || regression)
	if(qt)
	  {
	    YY1 = calloc2Ddouble(maxthreads,nlinestfam);
	    YY2 = calloc2Ddouble(maxthreads,nlinestfam);
	    Yd = calloc2Ddouble(maxthreads,nlinestfam);

	  }
	if(!qt)
	  {
	    Y = calloc2Dint(maxthreads,nlinestfam);

	  }

	//if(!plla2 && !plln) IS_561 deactivated
	{

	  p = (double *) calloc(nlinestfam,sizeof(double));
	  //	  cout<<dim1+maxIndexCov<< " dim1+maxIndexCov"<<endl;
	  S = (double *) calloc(dim1+maxIndexCov,sizeof(double));
	  YY = (double *) calloc(nlinestfam,sizeof(double));
	  Yhelp = (double *) calloc(nlinestfam,sizeof(double));


	  if (!p)
	    {
	      errorfile << "memory allocation error in p\n";
	      logfile << "memory allocation error in p\n";
	      cout << "memory allocation error in p\n";
	      errorfile.close();logfile.close();exit(1);
	    }
	  if (!S)
	    {
	      errorfile << "memory allocation error in S\n";
	      logfile << "memory allocation error in S\n";
	      cout << "memory allocation error in S\n";
	      errorfile.close();logfile.close();exit(1);
	    }
	  if (!YY)
	    {
	      errorfile << "memory allocation error in YY\n";
	      logfile << "memory allocation error in YY\n";
	      cout << "memory allocation error in YY\n";
	      errorfile.close();logfile.close();exit(1);
	    }
	  if (!Yhelp)
	    {
	      errorfile << "memory allocation error in Yhelp\n";
	      logfile << "memory allocation error in Yhelp\n";
	      cout << "memory allocation error in Yhelp\n";
	      errorfile.close();logfile.close();exit(1);
	    }

	  X = calloc2Ddouble(nlinestfam, dim1+maxIndexCov);
	  Xmod = calloc2Ddouble(nlinestfam, dim1+maxIndexCov);
	  Xt = calloc2Ddouble(dim1+maxIndexCov, nlinestfam);
	  A = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
	  VNN = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
	  Sinv = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
	  A0 = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
	  UNNT = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
	  Ainv = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
	  AinvXt = calloc2Ddouble(dim1+maxIndexCov, nlinestfam);
	  Yminusp = calloc2Ddouble(nlinestfam, 1);
	  newbeta = calloc2Ddouble(dim1+maxIndexCov, 1);
	  D = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
	  T = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
	  U = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
	  Ut = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
	  sumPP = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
	  sumPJ = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
	  sumPK = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
	  MMinv = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
	  if(qt)
	    {
	      Yt = calloc2Ddouble(1, nlinestfam);
	      YtX = calloc2Ddouble(1, dim1+maxIndexCov);
	      YtXAinv = calloc2Ddouble(1, dim1+maxIndexCov);
	      YtXAinvXt = calloc2Ddouble(1, nlinestfam);
	    }
	} // !plln !plla


	if (readmodel == 0)
	  {
	    // Parameterauswahl für die verschiedenen Modelle
	    //x[1] = x1 ,x[2] = x1D ,x[3] = x2 ,x[4] = x2D ,x[5] = x1x2 ,
	    //x[6] = x1x2D ,x[7] = x1Dx2  ,x[8] = x1Dx2D ,x[9] = x3  ,x[10] = x3D  ,
	    //x[11] = x1x3 ,x[12] = x1x3D ,x[13] = x1Dx3 ,x[14] = x1Dx3D ,x[15] = x2x3 ,
	    //x[16] = x2x3D ,x[17] = x2Dx3 ,x[18] = x2Dx3D ,x[19] = x1x2x3 ,x[20] = x1x2x3D ,
	    //x[21] = x1x2Dx3 ,x[22] = x1x2Dx3D ,x[23] = x1Dx2x3 ,x[24] = x1Dx2x3D ,x[25] = x1Dx2Dx3 ,x[26] = x1Dx2Dx3D

	    //xvec[0] = 1;
	    //zvec[0] = 1;
	    //for (i=1; i < 27; i++)
	    //{
	    //	xvec[i] = 0;
	    //	zvec[i] = 0;
	    //}

	    if (liability == 1 || liability == 2)
	      {
		xvec[1] = 1;
	      }
	    else if (test == 3) // full allelic test with 3 d.f.
	      {
		xvec[1] = 1;
		xvec[3] = 1;
		xvec[5] = 1;
	      }
	    else if (test == 4) // full genotpye test with 8 d.f.
	      {
		xvec[1] = 1;
		xvec[2] = 1;

		xvec[3] = 1;
		xvec[4] = 1;
		xvec[5] = 1;
		xvec[6] = 1;

		xvec[7] = 1;
		xvec[8] = 1;
	      }
	    else if (test == 5) //  test for additive interaction with 1 d.f.
	      {
		xvec[1] = 1;
		xvec[3] = 1;
		xvec[5] = 1;
		zvec[1] = 1;
		zvec[3] = 1;
	      }
	    else if (test == 6) // test for genotypic interaction with 4 d.f.
	      {
		xvec[1] = 1;
		xvec[2] = 1;
		xvec[3] = 1;
		xvec[4] = 1;
		xvec[5] = 1;
		xvec[6] = 1;
		xvec[7] = 1;
		xvec[8] = 1;
		zvec[1] = 1;
		zvec[2] = 1;
		zvec[3] = 1;
		zvec[4] = 1;
	      }
	    else if (test == 7) // Additional allelic effect of SNP2 with 2 d.f.
	      {
		xvec[1] = 1;
		xvec[3] = 1;
		xvec[5] = 1;
		zvec[1] = 1;
	      }
	    else if (test == 8) // Additional genotypic effect of SNP2 with 6 d.f.
	      {
		xvec[1] = 1;
		xvec[2] = 1;
		xvec[3] = 1;
		xvec[4] = 1;
		xvec[5] = 1;
		xvec[6] = 1;
		xvec[7] = 1;
		xvec[8] = 1;
		zvec[1] = 1;
		zvec[2] = 1;
	      }
	    else if (test == 9) // Full additive test with 7 d.f.
	      {
		xvec[1] = 1;
		xvec[3] = 1;
		xvec[9] = 1;
		xvec[5] = 1;
		xvec[11] = 1;
		xvec[15] = 1;
		xvec[19] = 1;
	      }
	    else if (test == 10) // test for allelic interaction with 4 d.f.
	      {
		xvec[1] = 1;
		xvec[3] = 1;
		xvec[9] = 1;
		xvec[5] = 1;
		xvec[11] = 1;
		xvec[15] = 1;
		xvec[19] = 1;
		zvec[1] = 1;
		zvec[3] = 1;
		zvec[9] = 1;
	      }
	    else if (test == 11) // test for 3-fold allelic interaction with 1 d.f.
	      {
		xvec[1] = 1;
		xvec[3] = 1;
		xvec[9] = 1;
		xvec[5] = 1;
		xvec[11] = 1;
		xvec[15] = 1;
		xvec[19] = 1;
		zvec[1] = 1;
		zvec[3] = 1;
		zvec[9] = 1;
		zvec[5] = 1;
		zvec[11] = 1;
		zvec[15] = 1;
	      }
	    else if (test == 12) // Additional allelic effect of third locus with 4 d.f.
	      {
		xvec[1] = 1;
		xvec[3] = 1;
		xvec[9] = 1;
		xvec[5] = 1;
		xvec[11] = 1;
		xvec[15] = 1;
		xvec[19] = 1;
		zvec[1] = 1;
		zvec[3] = 1;
		zvec[5] = 1;
	      }
	    else if (test == 13) // two-marker haplotype test
	      {
		xvec[1] = 1;
		xvec[2] = 1;
		xvec[3] = 1;
		xvec[4] = 1;
	      }
	    else if (test == 14) // three-marker haplotype test
	      {
		xvec[1] = 1;
		xvec[2] = 1;
		xvec[3] = 1;
		xvec[4] = 1;
		xvec[5] = 1;
		xvec[6] = 1;
		xvec[7] = 1;
		xvec[8] = 1;
	      }
	    else if (test == 17) // caseOnly1DF regression:
	      {
		xvec[1] = 1;
		xvec[3] = 1;
		zvec[1] = 1;
	      }
	    //CASE_ONLY
	    else if (test == 18) // CASEONLY18 4DF regression
	      {
		xvec[1] = 1;
		xvec[3] = 1;
		xvec[4] = 1;
		zvec[1] = 1;

		xvec2[1] = 1;
		xvec2[2] = 1;
		xvec2[3] = 1;
		xvec2[4] = 1;
		zvec2[1] = 1;
		zvec2[2] = 1;

	      }
	  }
      }
    else
      {
	Y = calloc2Dint(maxthreads,nlinestfam);
      }

    if (markercombi2 == 1 && test > 2)
      {
	for (i=9; i < 27; i++)
	  {
	    if (xvec[i] == 1)
	      {
		cout << "Do not choose a third SNP in logistic model with 2-marker analysis. "
		     << "Select TWO-MARKER = 0 and THREE-MARKER = 1 or modify logistic model\n";

		errorfile << "Do not choose a third SNP in logistic model with 2-marker analysis. "
			  << "Select TWO-MARKER = 0 and THREE-MARKER = 1 or modify logistic model\n";

		errorfile.close();

		logfile << "Do not choose a third SNP in logistic model with 2-marker analysis. "
			<< "Select TWO-MARKER = 0 and THREE-MARKER = 1 or modify logistic model\n";
		logfile.close();
		exit(1);
	      }
	  }
      }

    // Topliste
    int topListSize=singletop;
    if(printtop> singletop){topListSize=printtop;}
    toplist = (struct TOPLIST *) calloc(topListSize, sizeof(struct TOPLIST));

    if (!toplist)
      {
	errorfile << "memory allocation error in toplist\n";
	logfile << "memory allocation error in toplist\n";
	cout << "memory allocation error in toplist\n";
	errorfile.close();logfile.close();exit(1);
      }

    bestofSim = (double *) calloc(nsim+1,sizeof(double));

    if(!bestofSim)
      {
	errorfile << "memory allocation error in bestofSim\n";

	logfile << "memory allocation error in bestofSim\n";
	cout << "memory allocation error in bestofSim\n";
	errorfile.close();logfile.close();exit(1);
      }


    for (k=0; k<maxthreads;k++)
      {
	if(!qt)
	  {
	    for (i = 0; i < nlinestfam; i++)
	      {
		person[i].aff[k]=person[i].aff[thread];
		Y[k][i]=person[i].aff[thread];
	      }
	  }
	else
	  {
	    for (i = 0; i < nlinestfam; i++)
	      {
		person[i].aff[k]=person[i].aff[thread];
		person[i].qtaff[k]=person[i].qtaff[thread];

		Yd[k][i]  = person[i].qtaff[thread];
		YY1[k][i] = person[i].qtaff[thread];
		YY2[k][i] = person[i].qcin;
	      }
	  }
      } //k-llop maxthreads


    if(pathwayAnalysis==1 || pathwayImpact==1)
      {
	// Erstellen der Geneticlist in Kombination mit den Pathways
	for (r = rstart; r <= nlist; r++)
	  {
	    pathway[r-1].singletop=NULL;
	    pathway[r-1].singletop0=NULL;
	    pathway[r-1].overlaptop=NULL;
	    pathway[r-1].singletopmod=NULL;
	    pathway[r-1].intertop=NULL;
	    pathway[r-1].intertopmod=NULL;
	    pathway[r-1].score=NULL;
	    pathway[r-1].singletop = (int *) realloc(pathway[r-1].singletop, maxthreads*sizeof(int));
	    pathway[r-1].singletopmod = (int *) realloc(pathway[r-1].singletopmod, maxthreads*sizeof(int));
	    pathway[r-1].singletop0 = (int *) realloc(pathway[r-1].singletop0, maxthreads*sizeof(int));
	    pathway[r-1].overlaptop = (int *) realloc(pathway[r-1].overlaptop, maxthreads*sizeof(int));
	    pathway[r-1].score = (double *) realloc(pathway[r-1].score, maxthreads*sizeof(double)); //TIM 1.0.9
	    pathway[r-1].intertop = (int *) realloc(pathway[r-1].intertop, maxthreads*sizeof(int));
	    pathway[r-1].intertopmod = (int *) realloc(pathway[r-1].intertopmod, maxthreads*sizeof(int));


	    if (!pathway[r-1].singletop || !pathway[r-1].singletopmod || !pathway[r-1].singletop0 || !pathway[r-1].overlaptop || !pathway[r-1].score)
	      {
		errorfile << "memory allocation error in pathway[r-1]n";
		logfile << "memory allocation error in pathway[r-1]\n";
		cout << "memory allocation error in pathway[r-1]\n";
		errorfile.close();logfile.close();exit(1);
	      }

	    pathway[r-1].listsingle = NULL;
	    pathway[r-1].listgenetic = NULL;
	    pathway[r-1].listoverlap = NULL;
	    pathway[r-1].listp = NULL;
	    pathway[r-1].listsingle_gene = NULL;
	    pathway[r-1].listp_in = NULL;
	    pathway[r-1].listp_in2 = NULL;
	    pathway[r-1].listp_in3 = NULL;

	    pathway[r-1].listsingle = (int **) realloc(pathway[r-1].listsingle, maxthreads*sizeof(int *));
	    pathway[r-1].listoverlap = (int **) realloc(pathway[r-1].listoverlap, maxthreads*sizeof(int *));
	    pathway[r-1].listp = (double **) realloc(pathway[r-1].listp, maxthreads*sizeof(double *));
	    pathway[r-1].listp_in = (int **) realloc(pathway[r-1].listp_in, maxthreads*sizeof(int *));
	    pathway[r-1].listp_in2 = (int **) realloc(pathway[r-1].listp_in2, maxthreads*sizeof(int *));
	    pathway[r-1].listp_in3 = (int **) realloc(pathway[r-1].listp_in3, maxthreads*sizeof(int *));
	    pathway[r-1].listsingle_gene = (int **) realloc(pathway[r-1].listsingle_gene, maxthreads*sizeof(int *));

	    if (!pathway[r-1].listsingle || !pathway[r-1].listoverlap || !pathway[r-1].listp || !pathway[r-1].listp_in || !pathway[r-1].listp_in2 || !pathway[r-1].listp_in3 || !pathway[r-1].listsingle_gene)
	      {
		errorfile << "memory allocation error in pathway[r-1]n";
		logfile << "memory allocation error in pathway[r-1]\n";
		cout << "memory allocation error in pathway[r-1]\n";
		errorfile.close();logfile.close();exit(1);
	      }

	    for(k=0;k<maxthreads;k++)
	      {
		pathway[r-1].listsingle[k]=NULL;
		pathway[r-1].listoverlap[k]=NULL;
		pathway[r-1].listp[k]=NULL;
		pathway[r-1].listp_in[k]=NULL;
		pathway[r-1].listp_in2[k]=NULL;
		pathway[r-1].listp_in3[k]=NULL;
		pathway[r-1].listsingle_gene[k]=NULL;
	      }


	    pathway[r-1].genetictop = 0;

	    if (mWithGenetictop != 0)
	      {
		for (k = 0; k < genetictop; k++)
		  {
		    for(l = 0; l < pathway[r-1].counts; l++)
		      {
			if (geneticlist[k].nr == pathway[r-1].list[l])
			  {
			    pathway[r-1].genetictop++;
			    pathway[r-1].listgenetic = (int*) realloc(pathway[r-1].listgenetic, pathway[r-1].genetictop*sizeof(int));
			    if (!pathway[r-1].listgenetic)
			      {
				errorfile << "memory allocation error in pathway[r-1].listgenetic\n";
				logfile << "memory allocation error in pathway[r-1].listgenetic\n";
				cout << "memory allocation error in pathway[r-1].listgenetic\n";
				errorfile.close();logfile.close();exit(1);
			      }
			    pathway[r-1].listgenetic[pathway[r-1].genetictop-1] = geneticlist[k].nr;
			  }
		      }
		  }
	      }
	  } //end r
      } //


    if(!markercombi2 && !markercombi3 && !combilist && !pathwayAnalysis)
      {
	ntestsPlus=nsnpqc;
	ntestsUser=nsnpqc;
      }

    if(nsnps>=1)
      {
	mastersnp1=(char *) realloc(mastersnp1,(strlen(snp1)+1)*sizeof(char));
	strcpy(mastersnp1,snp1);
      }
    if(nsnps>=2)
      {
	mastersnp2=(char *) realloc(mastersnp2,(strlen(snp2)+1)*sizeof(char));
	strcpy(mastersnp2,snp2);
      }
    if(nsnps>=3)
      {
	mastersnp3=(char *) realloc(mastersnp3,(strlen(snp3)+1)*sizeof(char));
	strcpy(mastersnp3,snp3);
      }



    if(covariancematrix)
      {
	for (int i=1; i< 27;i++)
	  {
	    if(xvec[i]==1){df_L1++;}
	    if(zvec[i]==1){df_L2++;}
	    if(xsinglevec[i]==1){df_L1_Single++;}
	    if(zsinglevec[i]==1){df_L2_Single++;}
	  }
	/* This part can be used, if the need to get the elements of the covariance matrix corresponding to all the covariates will appear (Tanya, 24.01.2014)

	   int df_L1_param=df_L1+1; //+1 wegen beta_0
	   int df_L2_param=df_L2+1; //+1 wegen beta_0
	   int df_L1_Single_param=df_L1_Single+1; //+1 wegen beta_0
	   int df_L2_Single_param=df_L2_Single+1; //+1 wegen beta_0


	   int SigmaAuxArray[df_L1_param][df_L1_param];
	   int SigmaAuxVector[df_L1_param*(df_L1_param+1)/2];

	   for(int J=0; J<df_L1_param; J++)
	   {
	   int N=J;
	   int i=0;
	   //int SigmaAuxArray[df_L1_param][J];
	   SigmaAuxArray[0][J]=N;

	   for(int j=0; j<df_L1_param-1;j++)
	   {
	   if(j<J)
	   {
	   N=N+df_L1_param-1-i+ncov;
	   if(sexcov==1)
	   {
	   N=N+1;
	   }
	   }
	   else
	   {
	   N=N+1;
	   }
	   i+=1;
	   SigmaAuxArray[i][J]=N;
	   }
	   }


	   int VV=0;

	   for(int i=0; i<df_L1_param;i++)
	   {
	   for(int j=i;j<df_L1_param;j++)
	   {
	   SigmaAuxVector[VV]=SigmaAuxArray[i][j];
	   cout << "["<< i << "][" << j << "] " << SigmaAuxArray[i][j] << " [" << VV << "]" << SigmaAuxVector[VV] << endl;
	   VV+=1;
	   }
	   }


	   int SigmaAuxArray_Single[df_L1_Single_param][df_L1_Single_param];
	   int SigmaAuxVector_Single[df_L1_Single_param*(df_L1_Single_param+1)/2];

	   for(int J=0; J<df_L1_Single_param; J++)
	   {
	   int N=J;
	   int i=0;
	   //int SigmaAuxArray[df_L1_param][J];
	   SigmaAuxArray_Single[0][J]=N;

	   for(int j=0; j<df_L1_Single_param-1;j++)
	   {
	   if(j<J)
	   {
	   N=N+df_L1_Single_param-1-i+ncov;
	   if(sexcov==1)
	   {
	   N=N+1;
	   }
	   }
	   else
	   {
	   N=N+1;
	   }
	   i+=1;
	   SigmaAuxArray_Single[i][J]=N;
	   }
	   }


	   int VV_Single=0;

	   for(int i=0; i<df_L1_Single_param;i++)
	   {
	   for(int j=i;j<df_L1_Single_param;j++)
	   {
	   SigmaAuxVector_Single[VV_Single]=SigmaAuxArray_Single[i][j];
	   cout << "["<< i << "][" << j << "] " << SigmaAuxArray_Single[i][j] << " [" << VV_Single << "]" << SigmaAuxVector_Single[VV_Single] << endl;
	   VV_Single+=1;
	   }
	   }



	   //df_L1+=numberOfAllCov+1; //+1 wegen beta_0
	   //df_L2+=numberOfAllCov+1;
	   //df_L1_Single+=numberOfAllCov+1; //+1 wegen beta_0
	   //df_L2_Single+=numberOfAllCov+1;

	   // ncov includes the number of the covariates from the COVARIATEFILE, numberOfSNPCov and clusters covariates
	   // ncov does not include sexcov


	   df_L1+=ncov+1; //+1 wegen beta_0
	   df_L2+=ncov+1;
	   df_L1_Single+=ncov+1; //+1 wegen beta_0
	   df_L2_Single+=ncov+1;

	   if(sexcov==1)
	   {
	   df_L1+=1;
	   df_L2+=1;
	   df_L1_Single+=1;
	   df_L2_Single+=1;
	   }

	*/


	df_L1=df_L1+1; //+1 wegen beta_0
	df_L2=df_L2+1; //+1 wegen beta_0
	df_L1_Single=df_L1_Single+1; //+1 wegen beta_0
	df_L2_Single=df_L2_Single+1; //+1 wegen beta_0


	int SigmaAuxArray[df_L1][df_L1];

	SigmaAuxVector=(int *)calloc(df_L1*(df_L1+1)/2, sizeof(int));

	for(int J=0; J<df_L1; J++)
	  {
	    int N=J;
	    int i=0;
	    SigmaAuxArray[0][J]=N;

	    for(int j=0; j<df_L1-1;j++)
	      {
		if(j<J)
		  {
		    N=N+df_L1-1-i+ncov;
		    if(sexcov==1)
		      {
			N=N+1;
		      }
		  }
		else
		  {
		    N=N+1;
		  }
		i+=1;
		SigmaAuxArray[i][J]=N;
	      }
	  }


	int VV=0;

	for(int i=0; i<df_L1;i++)
	  {
	    for(int j=i;j<df_L1;j++)
	      {
		SigmaAuxVector[VV]=SigmaAuxArray[i][j];
		//cout << "["<< i << "][" << j << "] " << SigmaAuxArray[i][j] << " SigmaAuxVector[" << VV << "]=" << SigmaAuxVector[VV] << endl;
		VV+=1;
	      }
	  }


	int SigmaAuxArray_Single[df_L1_Single][df_L1_Single];

	SigmaAuxVector_Single=(int *)calloc(df_L1_Single*(df_L1_Single+1)/2, sizeof(int));

	for(int J=0; J<df_L1_Single; J++)
	  {
	    int N=J;
	    int i=0;
	    SigmaAuxArray_Single[0][J]=N;

	    for(int j=0; j<df_L1_Single-1;j++)
	      {
		if(j<J)
		  {
		    N=N+df_L1_Single-1-i+ncov;
		    if(sexcov==1)
		      {
			N=N+1;
		      }
		  }
		else
		  {
		    N=N+1;
		  }
		i+=1;
		SigmaAuxArray_Single[i][J]=N;
	      }
	  }


	int VV_Single=0;

	for(int i=0; i<df_L1_Single;i++)
	  {
	    for(int j=i;j<df_L1_Single;j++)
	      {
		SigmaAuxVector_Single[VV_Single]=SigmaAuxArray_Single[i][j];
		//cout << "["<< i << "][" << j << "] " << SigmaAuxArray_Single[i][j] << " SigmaAuxVector_Single[" << VV_Single << "]=" << SigmaAuxVector_Single[VV_Single] << endl;
		VV_Single+=1;
	      }
	  }


	dim_Cov1=(df_L1+ncov)*(df_L1+ncov+1)/2;
	dim_Cov2=(df_L2+ncov)*(df_L2+ncov+1)/2;
	dim_Cov1_Single=(df_L1_Single+ncov)*(df_L1_Single+ncov+1)/2;
	dim_Cov2_Single=(df_L2_Single+ncov)*(df_L2_Single+ncov+1)/2;

	if(sexcov==1)
	  {
	    dim_Cov1=(df_L1+ncov+1)*(df_L1+ncov+2)/2;
	    dim_Cov2=(df_L2+ncov+1)*(df_L2+ncov+2)/2;
	    dim_Cov1_Single=(df_L1_Single+ncov+1)*(df_L1_Single+ncov+2)/2;
	    dim_Cov2_Single=(df_L2_Single+ncov+1)*(df_L2_Single+ncov+2)/2;
	  }

      }


    for(j=0;j<maxthreads;j++)
      {
	init(&(resultSingle[j]),0,0,maxIndexCov,1,regression,dim_Cov1_Single,dim_Cov2_Single); // next ten lines: do not use var printBeta instead of regression!
	init(&(resultMulti[j]),haplo,0,maxIndexCov,1,regression,dim_Cov1,dim_Cov2);
      }



    for(j=0;j<maxthreads;j++)
      {
	init(&(result1[j]),haplo,0,maxIndexCov,1,regression,dim_Cov1,dim_Cov2);
	init(&(result2[j]),haplo,0,maxIndexCov,1,regression,dim_Cov1,dim_Cov2);
	init(&(result1Single[j]),haplo,0,maxIndexCov,1,regression,dim_Cov1_Single,dim_Cov2_Single);
	init(&(result2Single[j]),haplo,0,maxIndexCov,1,regression,dim_Cov1_Single,dim_Cov2_Single);
	//init(&(result1M[j]),haplo,0,maxIndexCov,1,regression);
	//init(&(result2M[j]),haplo,0,maxIndexCov,1,regression);
	//init(&(result1F[j]),haplo,0,maxIndexCov,1,regression);
	//init(&(result2F[j]),haplo,0,maxIndexCov,1,regression);
	init(&(result3[j]),haplo,0,maxIndexCov,1,regression,dim_Cov1,dim_Cov2);
	init(&(result4[j]),haplo,0,maxIndexCov,1,regression,dim_Cov1,dim_Cov2);
	//init(&(result3M[j]),haplo,0,maxIndexCov,1,regression);
	//init(&(result4M[j]),haplo,0,maxIndexCov,1,regression);
	//init(&(result3F[j]),haplo,0,maxIndexCov,1,regression);
	//init(&(result4F[j]),haplo,0,maxIndexCov,1,regression);
	firstCallThread[j]=0;
      }


    /// test covariates for signficance
    if (0 && regression && (maxIndexCov || sexcov) && !plla2 && !plln)
      {
	logg("");
	struct STATplus result1Dummy,result2Dummy;
	init(&(result1Dummy),0,0,maxIndexCov,1,regression,dim_Cov1,dim_Cov2);
	init(&(result2Dummy),0,0,maxIndexCov,1,regression,dim_Cov1,dim_Cov2);
	int dummyIndex=0;

	//get first valid SNP aus dummy index
	for (int iMod=0; iMod < nsnpqc ; iMod++)
	  {
	    dummyIndex=SNPMap[iMod];
	    if (map[dummyIndex].analysis_in == 1) {break;}
	  }

	for (int i=0;i<maxIndexCov+1;i++)
	  {
	    if (i==0 && !sexcov) continue;
	    newValue=1;
	    if(qt)
	      {
		alt=1;
		result1Dummy = qtreg(xsinglevec, cov, 0, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, 0, 0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1,alt,Yhelp, xType, female, male, counts[thread_nloop][dummyIndex], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,n,result1Dummy,0,0,0,0,3,i,0,0,0);

		if(result1Dummy.df>0)
		  {
		    alt=0;
		    result2Dummy = qtreg(zsinglevec, cov, 0, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, 0, 0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[thread_nloop][dummyIndex], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,n,result2Dummy,0,0,0,0,3,i,0,0,0);
		  }
		else
		  {
		    result2Dummy.df=0;result2Dummy.sc=0;
		  }

		df = result2Dummy.df -result1Dummy.df; //!!
		Fstat = ((result2Dummy.sc -result1Dummy.sc)/(df))/(result1Dummy.sc/result1Dummy.df);

		if (df >= 1 && result2Dummy.sc > 0.000001)
		  {
		    newValue = betai(result1Dummy.df/2,df/2,result1Dummy.df/(result1Dummy.df+df*Fstat));
		  }
	      }
	    else
	      {
		alt=1;
		result1Dummy = logreg(xsinglevec, cov, 0, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, 0, 0, 0, inflationfactor, casecounts5, controlcounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[0][dummyIndex], D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,npplqc, PPLMap, BinSNPs, PplLocations,n,result1Dummy,0,0,3,i,0,0,0,dosage,genoWeights);

		alt=0;
		if(result1Dummy.df>0)
		  {
		    result2Dummy = logreg(zsinglevec, cov, 0, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, 0, 0, 0, inflationfactor, casecounts5, controlcounts5, p, newbeta, X, Xmod, Xt, A,  UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[0][dummyIndex], D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,npplqc, PPLMap, BinSNPs, PplLocations,n,result2Dummy,0,0,3,i,0,0,0,dosage,genoWeights);
		  }
		else
		  {
		    result2Dummy.df=0;result2Dummy.sc=0;
		  }

		tLog = 2*(result1Dummy.sc -result2Dummy.sc);
		df = result1Dummy.df -result2Dummy.df;

		if (df >= 1)
		  {
		    newValue = pValueCalc(df/2, tLog/2);
		  }
	      }	// !qt

	    if (i==0) {
	      cout    << "P-value for sex covariate is " << newValue << endl;
	      logfile << "P-value for sex covariate is " << newValue << endl;
	    } else if (cov[i-1] == 1) {
	      cout    << "P-value for covariate " << i << " is " << newValue << endl;
	      logfile << "P-value for covariate " << i << " is " << newValue << endl;
	    }
	  } //i-loop maxIndexCov

	if(plla2 || plln)
	  {
	    free(p);p=NULL;free(S);S=NULL;
	    free(YY);YY=NULL;free(Yhelp);Yhelp = NULL;

	    free2Ddouble(X,nlinestfam);free2Ddouble(Xmod,nlinestfam);
	    free2Ddouble(Xt,dim1+maxIndexCov);free2Ddouble(A,dim1+maxIndexCov);
	    free2Ddouble(VNN,dim1+maxIndexCov);free2Ddouble(Sinv,dim1+maxIndexCov);
	    free2Ddouble(A0,dim1+maxIndexCov);free2Ddouble(UNNT,dim1+maxIndexCov);
	    free2Ddouble(Ainv,dim1+maxIndexCov);free2Ddouble(AinvXt,dim1+maxIndexCov);
	    free2Ddouble(Yminusp,nlinestfam);free2Ddouble(newbeta,dim1+maxIndexCov);
	    free2Ddouble(D,dim1+maxIndexCov);free2Ddouble(T,dim1+maxIndexCov);
	    free2Ddouble(U,dim1+maxIndexCov);free2Ddouble(Ut,dim1+maxIndexCov);
	    free2Ddouble(sumPP,dim1+maxIndexCov);free2Ddouble(sumPJ,dim1+maxIndexCov);
	    free2Ddouble(sumPK,dim1+maxIndexCov);free2Ddouble(MMinv,dim1+maxIndexCov);

	    if(qt)
	      {
		free2Ddouble(Yt,1);free2Ddouble(YtX,1);
		free2Ddouble(YtXAinv,1);free2Ddouble(YtXAinvXt,1);
	      }
	  }


      } // if regression compute sig of cov

    if(test>=15 && test<=19){caseOnly=1;}

    if(plusSingle && (!caseOnly || !markercombi2 && !markercombi3))
      {
	logfile << "PLUS_SINGLE was ignored. Active only with case-only interaction analysis.\n";
	cout << "PLUS_SINGLE was ignored. Active only with case-only interaction analysis.\n";
      }

    if (!pretest){pHelpCut=1;passedPreTest=1;}
  }


  uint64_t** BinSNPsCovCatFlags = NULL;
#if RARE
  if(bin && CMATflag && maxIndexCov>0){ // check once if cathegorical covariates are ok
    nCovCathegories=checkCovCathegories(cov, nlinestfam, person, nCovCathegories);
    BinSNPsCovCatFlags=getCovCat_bin_snp_gen(person, nlinestfam, nwordsSNPs, nCovCathegories);
  }

  /*
    #if devREGRESSION

    {
    if(REGRESSIONflag==1 || FRACREGflag==1 || COLLREGflag==1){regression=1;}
    if(FRACREGflag==1){collapseRare=1;}
    if(COLLREGflag==1){collcollapseRare=1;}
    if(regression || CMATflag || FISHERflag)
    {
    double weightvec[nlinestped];
    get_weights(nlinestped,weights,map,weightvec);
    YtXAinv = calloc2Ddouble(1, dim1+maxIndexCov);
    YtXAinvXt = calloc2Ddouble(1, nlinestfam);
    YtX = calloc2Ddouble(1, dim1+maxIndexCov);
    Yt = calloc2Ddouble(1, nlinestfam);
    Ainv = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
    A0 = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
    UNNT = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
    Sinv = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
    VNN = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
    S = (double *) calloc(dim1+maxIndexCov,sizeof(double));
    YY = (double *) calloc(nlinestfam,sizeof(double));
    Yhelp = (double *) calloc(nlinestfam,sizeof(double));
    X = calloc2Ddouble(nlinestfam, dim1+maxIndexCov);
    Xt = calloc2Ddouble(dim1+maxIndexCov, nlinestfam);
    A = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);

    xType=0;
    struct STATplus result1Dummy,result2Dummy;
    int collapseRare =1; // 1=do Morris rare-method
    int collcollapseRare =1; // 1=do Morris other rare-method
    init(&(result1Dummy),0,0,maxIndexCov,1,regression,dim_Cov1,dim_Cov2);
    init(&(result2Dummy),0,0,maxIndexCov,1,regression,dim_Cov1,dim_Cov2);
    result1Dummy.betaNew_se[0]=0;

    //paramters & Co for L1
    int nSnps=100; //number of snps to be included with additive term
    if(nSnps > nsnpqc)
    {
    //    cout << "\nError. nSnps is larger than nsnpsqc.\n";exit(1);
    nSnps=nsnpqc;
    }
    int nSnpsDom=0; //number of snps to be included with dominance term
    int nEnvirons=0; // number of non-genetic parameters (for future extensions)
    int nInters=0; // number of interaction terms to be used (under construction)
    int *snps; // list of snps to be included with additive term, indices in BinSNPs
    int *snpsDom; // list of snps to be included wit dominance term, indices in BinSNPs
    int *environs; // list of non-genetic parameters (for future extensions)
    int *inters; // list of interaction terms to be used (coding under construction)
    int nx=nSnps+nSnpsDom+nEnvirons+nInters; // total n
    if(sexcov==1){nx++;}

    snps = (int *) calloc(nSnps, sizeof(int));
    snpsDom = (int *) calloc(nSnpsDom, sizeof(int));
    environs = (int *) calloc(nEnvirons, sizeof(int));
    inters = (int *) calloc(nInters, sizeof(int));

    if(nSnps > nsnpqc)
    {
    cout << "\nError. nSnps "<<nSnps<< " is larger than nsnpqc "<<nsnpqc<< ".\n";exit(1);
    }

    for(int h=0; h<nSnps;h++)
    {
    snps[h]=h; //index of snps
    //snpsDom[h]=snps[h];
    }

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
    snpsDom_second = (int *) calloc(nSnpsDom_second, sizeof(int));
    environs_second = (int *) calloc(nEnvirons_second, sizeof(int));
    inters_second = (int *) calloc(nInters_second, sizeof(int));

    //snps_second[0]=0; //index of snps

    alt=1;
    if(qt){
    result1Dummy = regGeneral(xsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, 0, 0, 0,
    inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY,
    Yminusp, N1,alt,Yhelp, xType, female, male, counts[0][0], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut,
    sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap,
    BinSNPs, PplLocations,n,result1Dummy,0,0,maxIndexCov,0,singleMarkerTest,-1,
    nSnps,nSnpsDom,nEnvirons,nInters,snps,snpsDom,environs,inters,nx,dim1,collapseRare,collcollapseRare,qt,weightvec);
    }
    else{
    result1Dummy = logRegGeneral(xsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, 0, 0, 0,
    inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY,
    Yminusp, N1,alt,Yhelp, xType, female, male, counts[0][0], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut,
    sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap,
    BinSNPs, PplLocations,n,result1Dummy,0,0,maxIndexCov,0,singleMarkerTest,-1,
    nSnps,nSnpsDom,nEnvirons,nInters,snps,snpsDom,environs,inters,nx,dim1,collapseRare,collcollapseRare,qt,weightvec,covariancematrix);
    }

    cout << result1Dummy.df << " df1\n";
    cout << result1Dummy.sc << " sc1\n";
    if(result1Dummy.df>0)
    {
    alt=0;
    cout << nx_second << " nx_second\n";
    if(qt){
    result2Dummy = regGeneral(zsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, 0, 0, 0,
    inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY,
    Yminusp, N1, alt, Yhelp, xType, female, male, counts[0][0], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut,
    sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap,
    BinSNPs,  PplLocations,n,result2Dummy,0,0,maxIndexCov,0,singleMarkerTest,-1,
    nSnps_second,nSnpsDom_second,nEnvirons_second,nInters_second,snps_second,snpsDom_second,
    environs_second,inters_second,nx_second,dim1,collapseRare,collcollapseRare,qt,weightvec);
    }
    else{
    result2Dummy = logRegGeneral(zsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, 0, 0, 0,
    inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY,
    Yminusp, N1, alt, Yhelp, xType, female, male, counts[0][0], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut,
    sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap,
    BinSNPs,  PplLocations,n,result2Dummy,0,0,maxIndexCov,0,singleMarkerTest,-1,
    nSnps_second,nSnpsDom_second,nEnvirons_second,nInters_second,snps_second,snpsDom_second,
    environs_second,inters_second,nx_second,dim1,collapseRare,collcollapseRare,qt,weightvec,covariancematrix);
    }

    cout << result2Dummy.df << " df2\n";
    cout << result2Dummy.sc << " sc2\n";
    }
    else
    {
    result2Dummy.df=0;result2Dummy.sc=0;
    }

    df = result2Dummy.df -result1Dummy.df;
    Fstat = ((result2Dummy.sc -result1Dummy.sc)/(df))/(result1Dummy.sc/result1Dummy.df);

    if (df >= 1 && result2Dummy.sc > 0.000001)
    {
    newValue = betai(result1Dummy.df/2,df/2,result1Dummy.df/(result1Dummy.df+df*Fstat));
    }
    else {newValue=1;}

    cout << "\nTest output p-value from function regGeneral. A model with " << nx << " parameters is tested. ";

    cout << "p-value: " << newValue << ".\n\n";

    }
    }
    #endif
  */
#endif


  cout << "\nPerforming analysis. This might take a while ..." << endl;
  bool barrier = true;
  long int m = 1;

  int multimarker=0;
  if(markercombi2 || markercombi3){multimarker=1;}

  char number[5]="0000";
  if(writeOut && multimarker) //writeOut
    {
      for(int h=0;h<MAXTHREAD;h++) //writeout
	{
	  sprintf(number,"%d_",h+1);
	  markerCombifileParallel[h] = outputname + number + "BestMarkerCombi.txt";
	  out2header(markerCombifileParallel[h],pathwayImpact, haplo, test, covariancematrix, df_L1);
	}
      if(pfilter>0.00001)
	{
	  errorfile << "\nWarning!! Writeout is active. PFILTER has to be chosen reasonnably low to avoid huge output file size in pair-wise analysis!\n\n";
	  logfile << "\nWarning!! Writeout is active. PFILTER has to be chosen reasonnably low to avoid huge output file size in pair-wise analysis!\n\n";
	  cout << "\nWarning!! Writeout is active. PFILTER has to be chosen reasonnably low to avoid huge output file size in pair-wise analysis!\n\n";
	}

    }



  bool abortnsim=false;
#if PARALLELN
#pragma omp parallel
  {
    srand(int(time(NULL)) ^ omp_get_thread_num());

#pragma omp for schedule(dynamic,1) private(thread_aloop,aMod,bMod,cMod,traitavg,pHelp,snp1,snp2,snp3,pos1,pos2,pos3,maxposSingle,maxposMulti,lastpSingle,lastpMulti,l,skip,result0,kkk,jjj,lll,mmm,p_2nd,p_3rd,ntests,ntestsUser,n,a2,helpi,overlapneeded,selectedSnp,j,teststat,helpstat,newValue,newValueCorr,i,stopD,stopE,N,bestchi3,fulltests,inflationfactor,r,refSnp,snp1Pos,snp2Pos,snp3Pos,helpk,helpj,k,overlaptop,startA,stopA,stopB,q,bestsinglemarker2,geneticlist2,overlaplist,jstart,bestsinglemarker,thread_nloop,b2,c2,d,e,countA,startB,a,a_list,countB,b,startC,b_list,firstc,stopC,countC,c,c_list,d_list,e_list,x,oldCombi,ifChrX,c_in,d_in,e_in,tstat,stopCC,stopDD,stopEE,aa,bb,cc,dd,ee,controlcounts5,controlcounts5M,controlcounts5F,casecounts5,casecounts5M,casecounts5F,index,alt,xType,tLog,df,Fstat,tLogMale,dfMale,tLogFemale,dfFemale,marker1,marker2,marker3,ii,jj,kk,ll,mm,p,S,YY,Yhelp,X,Xmod,Xt,A,VNN,Sinv,A0,UNNT,Ainv,AinvXt,Yminusp,newbeta,D,T,U,Ut,sumPP,sumPJ,sumPK,MMinv,Yt,YtX,YtXAinv,YtXAinvXt,ll2,stat,tstat2,filled,passedPreTest)
#endif
    for (n = 0; n <= nsim; n++) // MC-simulations
      {
#if PARALLELN
	if (n) while (barrier) wait(10);
#endif

#pragma omp flush (abortnsim)
	if(bin && bins2testt==0){
	  abortnsim=true;
	}
	if(!abortnsim){
	  pHelp=1;
	  maxposSingle=0;
	  maxposMulti=0;
	  lastpSingle=1.1;
	  lastpMulti=1.1;
	  selectedSnp=nsnps;
	  if(combilist && markercombi2)
	    {
	      selectedSnp=2;
	    }
	  else if(combilist && markercombi3)
	    {
	      selectedSnp=3;
	    }


	  N=1; //to be on the safe side

	  if(plln)
	    {
	      if(nsnps>=1) {snp1= (char *) calloc(strlen(mastersnp1)+1,sizeof(char));}
	      if(nsnps>=2) {snp2= (char *) calloc(strlen(mastersnp2)+1,sizeof(char));}
	      if(nsnps>=3) {snp3= (char *) calloc(strlen(mastersnp3)+1,sizeof(char));}
	    }

	  if(combilist)
	    {
	      if(!plln && n==0 || plln)
		{
		  snp1= (char *) calloc(maxSizeNew,sizeof(char));
		  snp2= (char *) calloc(maxSizeNew,sizeof(char));
		  snp3= (char *) calloc(maxSizeNew,sizeof(char));
		}

	    }
	  else
	    {
	      if(nsnps>=1) {strcpy(snp1,mastersnp1);}
	      if(nsnps>=2) {strcpy(snp2,mastersnp2);}
	      if(nsnps>=3) {strcpy(snp3,mastersnp3);}
	    }


	  pos1=-1;pos2=-1;pos3=-1;

	  helpi=0;ntestsUser=0;skip=0;ntests=0;

	  marker1=0;marker2=0;marker3=0;

	  startA=0;stopA=0;startB=0;stopB=0;
	  thread_nloop=0;overlapneeded=0;
	  overlaptop=0;
#if PARALLELN
	  thread_nloop=omp_get_thread_num();
	  if(thread_nloop > maxthreads){printf("error\n");exit(1);}
	  //YY=NULL;
#endif
	  if(firstCallThread[thread_nloop]==0){firstCallThread[thread_nloop]++;}
	  //bench=1;
	  fulltests = 0;
	  inflationfactor = 0;

	  //cout << thread_nloop << " thread\n";

	  if(n==0 || plln)
	    {
	      bestsinglemarker = (struct BESTSINGLEMARKER *) calloc(singletop, sizeof(struct BESTSINGLEMARKER));
	      if(!bestsinglemarker)
		{
		  errorfile << "memory allocation error in bestsinglemarker\n";
		  logfile << "memory allocation error in bestsinglemarker\n";
		  cout << "memory allocation error in bestsinglemarker\n";
		  errorfile.close();logfile.close(); exit(1);
		}
	      if(multimarker)
		{
		  bestchi3 = (struct BESTCHI5 *) calloc(printtop, sizeof(struct BESTCHI5));

		  if (!bestchi3)
		    {

		      errorfile << "memory allocation error in bestchi3\n";
		      logfile << "memory allocation error in bestchi3\n";
		      cout << "memory allocation error in bestchi3\n";
		      errorfile.close();logfile.close();exit(1);
		    }
		}
	    }

	  if (n > 0)
	    {
	      // Affectionstatus permutieren
	      if(!qt)
		{
		  switch (stratify) {
		  case 1:
		  case 5:
		  case 8:
		    aff_permute_pairs(person, MatchedPairs, nMatches, thread_nloop);
		    break;
		  case 2:
		  case 3:
		  case 4:
		  case 6:
		  case 7:
		  case 9:
		  case 10:
		    aff_permute_clusters(person, Clusters, nClusters, thread_nloop, &ix, &iy, &iz);
		    break;
		  default:
		    aff_permute_all(person, nlinestfam, ncasesqc, thread_nloop, &ix, &iy, &iz);
		  }
		}
	      else // qt
		{
		  switch (stratify) {
		  case 4:
		    qt_permute_clusters(person, Clusters, nClusters, thread_nloop, &ix, &iy, &iz);
		    break;
		  default:
		    //                            qt_permute_all(person, nlinestfam, thread, &ix, &iy, &iz);
		    qt_permute(nlinestfam, &(Yd[thread_nloop]), &ix, &iy, &iz, person, ncasesqc, YY1[thread_nloop], YY2[thread_nloop]);
		    for (i = 0; i < nlinestfam; i++)  // write permuted dummy into person
		      {
			person[i].qtaff[thread_nloop] = Yd[thread_nloop][i];
		      }
		    for (i = 0; i < nlinestfam; i++)  // reset of dummy-arrays, because of side-effects
		      {
			Yd[thread_nloop][i]=person[i].qtaff[thread_nloop];
			YY1[thread_nloop][i] = person[i].qtaff[thread_nloop];
			if(person[i].qcin ==1) { YY2[thread_nloop][i] = 1; }
			else {YY2[thread_nloop][i] = 0;}
		      }
		  }

		}   //qt
	    } // if n >0
	  uint64_t** BinSNPsCCFlags = getCC_bin_snp_gen(person, nlinestfam, nwordsSNPs);
	  if(caseOnly && markercombi2 && (test==17 || test==18)) {BinSNPsQTaffFlags = getQTaff2_bin_snp_gen(person, nlinestfam, nwordsSNPs);}


	  for (r = rstart; r <= nlist; r++) // Pathways
	    {

	      if (regression == 1 && r==rstart && (plln || (plla2 && r<=1))) //plla2 wegen singlemarker // TIM 1.0.9
		{
		  //		cout << "r n " << r << " " << n << "\n";
		  p = (double *) calloc(nlinestfam,sizeof(double));
		  if (!p)
		    {
		      errorfile << "memory allocation error in p\n";
		      logfile << "memory allocation error in p\n";
		      cout << "memory allocation error in p\n";
		      errorfile.close();logfile.close();exit(1);
		    }

		  S = (double *) calloc(dim1+maxIndexCov,sizeof(double));
		  if (!S)
		    {
		      errorfile << "memory allocation error in S\n";
		      logfile << "memory allocation error in S\n";
		      cout << "memory allocation error in S\n";
		      errorfile.close();logfile.close();exit(1);
		    }
		  //YY = (double *) realloc(YY,0*sizeof(double));
		  YY = (double *) calloc(nlinestfam,sizeof(double));
		  if (!YY)
		    {
		      errorfile << "memory allocation error in YY\n";
		      logfile << "memory allocation error in YY\n";
		      cout << "memory allocation error in YY\n";
		      errorfile.close();logfile.close();exit(1);
		    }
		  //Yhelp = (double *) realloc(Yhelp,0*sizeof(double));
		  Yhelp = (double *) calloc(nlinestfam,sizeof(double));
		  if (!Yhelp)
		    {
		      errorfile << "memory allocation error in Yhelp\n";
		      logfile << "memory allocation error in Yhelp\n";
		      cout << "memory allocation error in Yhelp\n";
		      errorfile.close();logfile.close();exit(1);
		    }

		  X = calloc2Ddouble(nlinestfam, dim1+maxIndexCov);
		  Xmod = calloc2Ddouble(nlinestfam, dim1+maxIndexCov);
		  Xt = calloc2Ddouble(dim1+maxIndexCov, nlinestfam);
		  A = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
		  VNN = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
		  Sinv = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
		  A0 = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
		  UNNT = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
		  Ainv = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
		  AinvXt = calloc2Ddouble(dim1+maxIndexCov, nlinestfam);
		  Yminusp = calloc2Ddouble(nlinestfam, 1);
		  newbeta = calloc2Ddouble(dim1+maxIndexCov, 1);
		  D = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
		  T = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
		  U = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
		  Ut = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
		  sumPP = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
		  sumPJ = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
		  sumPK = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
		  MMinv = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);

		  if(qt)
		    {
		      Yt = calloc2Ddouble(1, nlinestfam);
		      YtX = calloc2Ddouble(1, dim1+maxIndexCov);
		      YtXAinv = calloc2Ddouble(1, dim1+maxIndexCov);
		      YtXAinvXt = calloc2Ddouble(1, nlinestfam);
		    }

		} //plla2 plln

	      for (q = qstart; q <= nlinescombi; q++) // combilist
		{

		  if (combilist == 1)
		    {
		      if (markercombi2 == 1)
			{

			  strcpy(snp1,markerTable[q-1][0].rs);
			  //exit(1);
			  strcpy(snp2,markerTable[q-1][1].rs);

			  pos1 = markerTable[q-1][0].pos; // Änderung569
			  pos2 = markerTable[q-1][1].pos; // Änderung569
			  pos3 = -1;
			  storeda = -1;
			  storedb = -1;
			}
		      else if (markercombi3 == 1)
			{
			  strcpy(snp1,markerTable[q-1][0].rs);
			  strcpy(snp2,markerTable[q-1][1].rs);
			  strcpy(snp3,markerTable[q-1][2].rs);
			  pos1 = markerTable[q-1][0].pos; // Änderung569
			  pos2 = markerTable[q-1][1].pos; // Änderung569
			  pos3 = markerTable[q-1][2].pos; // Änderung569
			}
		    }

		  if (r == rstart && q<=1)
		    {

		      if (markercombi2 == 1)
			{
			  stopC = 1;
			  stopD = 1;
			  stopE = 1;
			  N=2;
			}
		      else if (markercombi3 == 1)
			{
			  stopC = 3;
			  stopD = 1;
			  stopE = 1;
			  N = 3;
			}

		      if(multimarker)
			{
			  for (i = 0; i < printtop; i++)
			    {
			      bestchi3[i].p = 1.01;
			      bestchi3[i].pmod = 1.01;
			      bestchi3[i].nr1 = -1;
			      bestchi3[i].nr2 = -1;
			      bestchi3[i].nr3 = -1;
			      bestchi3[i].nr4 = -1;
			      bestchi3[i].nr5 = -1;
			      bestchi3[i].r = -1;
			      bestchi3[i].X = 0;
			      initTRAITavg(&(bestchi3[i].traitavg));

			      for (ii= 0; ii < 3; ii++)
				{
				  for (jj= 0; jj < 3; jj++)
				    {
				      for (kk = 0; kk < 3; kk++)
					{
					  for (ll= 0; ll < 1; ll++)
					    {
					      for (mm = 0; mm < 1; mm++)
						{
						  bestchi3[i].casecounts[ii][jj][kk][ll][mm] = -1;
						  bestchi3[i].casecountsMale[ii][jj][kk][ll][mm] = -1;
						  bestchi3[i].casecountsFemale[ii][jj][kk][ll][mm] = -1;
						  bestchi3[i].controlcounts[ii][jj][kk][ll][mm] = -1;
						  bestchi3[i].controlcountsMale[ii][jj][kk][ll][mm] = -1;
						  bestchi3[i].controlcountsFemale[ii][jj][kk][ll][mm] = -1;
						}
					    }
					}
				    }
				}
			    }
			}
		    }

		  if (r <=1) // only once!!!
		    {
		      // Speichern der besten 1000 singlemarker-Werte und Singlemarker-Ausgabe-File
		      if (q == qstart)
			{
			  for (i = 0; i < singletop; i++)
			    {
			      bestsinglemarker[i].p = 1.01;
			      bestsinglemarker[i].pmod = 1.01;
			      bestsinglemarker[i].nr = -1;
			    }

			}


		      if(q==qstart)
			{   //TIM

			  for (int kMod = 0; kMod < nsnpqc; kMod++)
			    {

			      k = SNPMap[kMod];
			      if (!map[k].analysis_in) continue;

			      initCounts(&counts[thread_nloop][k],n,regression);
			      //counts[thread_nloop][k].result1.b[0] = 0;
			      //counts[thread_nloop][k].result1.b[1] = 0;
			      //counts[thread_nloop][k].result1.b[2] = 0;
			      if(n==0 || firstCallThread[thread_nloop]==0)
				{
				  init(&(*(counts[thread_nloop][k].result1)),0,0,maxIndexCov,0,regression,dim_Cov1,dim_Cov2); // do not replace regression with printBeta here
				  if(regression)
				    {
				      (*counts[thread_nloop][k].result1).betaNew_se[1] = 0;
				      (*counts[thread_nloop][k].result1).oddsRatio[1] = 0;
				      (*counts[thread_nloop][k].result1).lcloddsRatio[1] = 0;
				      (*counts[thread_nloop][k].result1).rcloddsRatio[1] = 0;
				      (*counts[thread_nloop][k].result1).betaNew_lcl[1] = 0;
				      (*counts[thread_nloop][k].result1).betaNew_rcl[1] = 0;

				      (*counts[thread_nloop][k].result1).betaNew_se[2] = 0;
				      (*counts[thread_nloop][k].result1).oddsRatio[2] = 0;
				      (*counts[thread_nloop][k].result1).lcloddsRatio[2] = 0;
				      (*counts[thread_nloop][k].result1).rcloddsRatio[2] = 0;
				      (*counts[thread_nloop][k].result1).betaNew_lcl[2] = 0;
				      (*counts[thread_nloop][k].result1).betaNew_rcl[2] = 0;

				      (*counts[thread_nloop][k].result1).b[0] = 0;
				      (*counts[thread_nloop][k].result1).b[1] = 0;
				      (*counts[thread_nloop][k].result1).b[2] = 0;
				    }
				}

			      helpk = k;

			      updateCounts(&counts[thread_nloop][k], nwordsSNPs, BinSNPs[kMod], BinSNPsCCFlags, BinSNPsGenderFlags, !(strcmp(map[k].chr, "23") && strcmp(map[k].chr, "24")));

			    } // end kMod
			} //endif qstart



		      for (int iMod=0; iMod < nsnpqc ; iMod++)
			{
			  i=SNPMap[iMod];
			  if (!map[i].analysis_in || map[i].done) { continue; }
			  teststat = 0;

			  if(q==qstart) //TIM
			    {
			      newValue = 1;
			      helpi = i;

			      // Singlemarker NEU

			      // HWE und Singlemarker-Ausgabe-File
			      if (singleMarkerTest == 1)
				{
				  if (strcmp(map[i].chr, "23") == 0)
				    {
				      newValue = snpTestX(counts[thread_nloop][i], teststat,&inflationfactor, &fulltests);
				    }
				  else {
				    if (group_test) {
				      switch (stratify) {
				      case 1:
				      case 5:
				      case 8:
					if (group_test==1) {
					  newValue = armitageTrendPairTestLinear(MatchedPairs, nMatches, BinSNPs[iMod], BinSNPsCCFlags, PplLocFam, teststat, &inflationfactor, &fulltests);
					} else if (group_test==2) {
					  newValue = armitageTrendPairTestSquare(MatchedPairs, nMatches, BinSNPs[iMod], BinSNPsCCFlags, PplLocFam, teststat, &inflationfactor, &fulltests);
					}
					break;
				      case 2:
				      case 3:
				      case 4:
				      case 6:
				      case 7:
				      case 9:
				      case 10:
					if (group_test==1) {
					  newValue = armitageTrendGroupTestLinear(Clusters, nClusters, BinSNPs[iMod], BinSNPsCCFlags, PplLocFam, teststat, &inflationfactor, &fulltests);
					} else if (group_test==2) {
					  newValue = armitageTrendGroupTestSquare(Clusters, nClusters, BinSNPs[iMod], BinSNPsCCFlags, PplLocFam, teststat, &inflationfactor, &fulltests);
					}
					break;
				      default:
					newValue = armitageTest(counts[thread_nloop][i], teststat, &inflationfactor, &fulltests);
				      }
				    } else {
				      newValue = armitageTest(counts[thread_nloop][i], teststat, &inflationfactor, &fulltests);
				    }
				  }
				}
			      else if (singleMarkerTest == 2)
				{
				  if (strcmp(map[i].chr, "23") == 0)
				    {
				      newValue = snpTestX(counts[thread_nloop][i], teststat,&inflationfactor, &fulltests);
				    }
				  else
				    {
				      newValue = genotypTest(counts[thread_nloop][i], teststat,&inflationfactor, &fulltests);
				    }
				}
			      else if ((singleMarkerTest==3 || singleMarkerTest >=5) || singleMarkerTest == 4)
				{
				  if (strcmp(map[i].chr, "23") == 0)
				    {
				      alt=1;
				      xType = 3;
				    }
				  else
				    {
				      alt=1;
				      xType = 0;
				    }
				  if(!qt)
				    {
				      //counts[thread_nloop][i].result1 = logreg(xsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, iMod, 0, 0, inflationfactor, casecounts5, controlcounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[thread_nloop][i], D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,npplqc, PPLMap, BinSNPs, PplLocations,n,resultSingle[thread_nloop],maxIndexCov,liabilityCut);
				      result1Single[thread_nloop] = logreg(xsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, iMod, 0, 0, inflationfactor, casecounts5, controlcounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[thread_nloop][i], D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,npplqc, PPLMap, BinSNPs, PplLocations,n,resultSingle[thread_nloop],maxIndexCov,0, singleMarkerTest,-1,covariancematrix,df_L1_Single,df_L2_Single,dosage,genoWeights);
				      PlusToPlus(&(*(counts[thread_nloop][i].result1)),result1Single[thread_nloop],0,maxIndexCov,regression,0,covariancematrix,dim_Cov1_Single,dim_Cov2_Single);
				      alt=0;
				      if(result1Single[thread_nloop].df>0)
					{

					  // resultCopy(&(resultSingle[thread_nloop]),result2Single[thread_nloop]);
					  //resultCopy(&(resultSingle_Rare[thread_nloop]),result2Single[thread_nloop]);
					  result2Single[thread_nloop] = logreg(zsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, iMod, 0, 0, inflationfactor, casecounts5, controlcounts5, p, newbeta, X, Xmod, Xt, A,  UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[thread_nloop][i], D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,npplqc, PPLMap, BinSNPs, PplLocations,n,resultSingle[thread_nloop],maxIndexCov,0, singleMarkerTest,-1,covariancematrix,df_L1_Single,df_L2_Single,dosage,genoWeights);
					}
				      else
					{
					  result2Single[thread_nloop].df=0;result2Single[thread_nloop].sc=0;
					}

				      tLog = 2*(result1Single[thread_nloop].sc -result2Single[thread_nloop].sc);
				      df = result1Single[thread_nloop].df -result2Single[thread_nloop].df;

				      if (df >= 1)
					{
					  newValue = pValueCalc(df/2, tLog/2);
					  fulltests++;
					  if ((singleMarkerTest==3 || singleMarkerTest >=5) || df ==2)
					    {inflationfactor+=tLog;}
					  else {inflationfactor+=2*tLog;}
					}
				      else
					{
					  newValue = 1;
					}

				    }
				  else //qt singlemarker
				    {
				      newValue = 1;

				      result1Single[thread_nloop] = qtreg(xsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, iMod, 0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1,alt,Yhelp, xType, female, male, counts[thread_nloop][i], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,n,resultSingle[thread_nloop],0,0,maxIndexCov,0,singleMarkerTest,-1,covariancematrix,df_L1_Single,df_L2_Single);
				      PlusToPlus(&(*(counts[thread_nloop][i].result1)),result1Single[thread_nloop],0,maxIndexCov,regression,0,covariancematrix,dim_Cov1_Single,dim_Cov2_Single);

				      if(result1Single[thread_nloop].df>0)
					{
					  alt=0;
					  result2Single[thread_nloop] = qtreg(zsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, iMod, 0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[thread_nloop][i], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs,  PplLocations,n,resultSingle[thread_nloop],0,0,maxIndexCov,0,singleMarkerTest,-1,covariancematrix,df_L1_Single,df_L2_Single);
					}
				      else
					{
					  result2Single[thread_nloop].df=0;result2Single[thread_nloop].sc=0;
					}

				      df = result2Single[thread_nloop].df -result1Single[thread_nloop].df; //!!
				      Fstat = ((result2Single[thread_nloop].sc -result1Single[thread_nloop].sc)/(df))/(result1Single[thread_nloop].sc/result1Single[thread_nloop].df);

				      if (df >= 1 && result2Single[thread_nloop].sc > 0.000001)
					{
					  newValue = betai(result1Single[thread_nloop].df/2,df/2,result1Single[thread_nloop].df/(result1Single[thread_nloop].df+df*Fstat));
					  fulltests++;
					  if ((singleMarkerTest==3 || singleMarkerTest >=5) || df ==2)
					    {inflationfactor+=Fstat;}
					  else {inflationfactor+=2*Fstat;}
					}
				      else
					{
					  newValue = 1;
					}
				    }
				}

			      if(n==0){newValueCorr = Sidak(newValue, nsnpAnalysisqc);}

			      // HWE und Counts
			      if (n == 0)
				{

				  if (strcmp(map[i].chr, "23") == 0)
				    {
				      // HWE
				      counts[thread_nloop][i].det->testHWE_Ca = testHWE(counts[thread_nloop][i].AA_Ca_female, counts[thread_nloop][i].AB_Ca_female,counts[thread_nloop][i].BB_Ca_female);
				      counts[thread_nloop][i].det->testHWE_Co = testHWE(counts[thread_nloop][i].AA_Co_female, counts[thread_nloop][i].AB_Co_female,counts[thread_nloop][i].BB_Co_female);

				      counts[thread_nloop][i].det->aCaN = counts[thread_nloop][i].AA_Ca_male + 2* counts[thread_nloop][i].AA_Ca_female+ counts[thread_nloop][i].AB_Ca_female;
				      counts[thread_nloop][i].det->aCoN = counts[thread_nloop][i].AA_Co_male + 2* counts[thread_nloop][i].AA_Co_female+ counts[thread_nloop][i].AB_Co_female;
				      counts[thread_nloop][i].det->bCaN = counts[thread_nloop][i].BB_Ca_male + 2* counts[thread_nloop][i].BB_Ca_female+ counts[thread_nloop][i].AB_Ca_female;
				      counts[thread_nloop][i].det->bCoN = counts[thread_nloop][i].BB_Co_male + 2* counts[thread_nloop][i].BB_Co_female+ counts[thread_nloop][i].AB_Co_female;
				    }
				  else if ((strcmp(map[i].chr, "24") == 0) || (strcmp(map[i].chr, "26") == 0))
				    {

				      // HWE
				      counts[thread_nloop][i].det->testHWE_Ca = 1;
				      counts[thread_nloop][i].det->testHWE_Co = 1;

				      counts[thread_nloop][i].det->aCaN = counts[thread_nloop][i].AA_Ca + counts[thread_nloop][i].AB_Ca;
				      counts[thread_nloop][i].det->aCoN = counts[thread_nloop][i].AA_Co + counts[thread_nloop][i].AB_Co;
				      counts[thread_nloop][i].det->bCaN = counts[thread_nloop][i].BB_Ca + counts[thread_nloop][i].AB_Ca;
				      counts[thread_nloop][i].det->bCoN = counts[thread_nloop][i].BB_Co + counts[thread_nloop][i].AB_Co;
				    }
				  else
				    {
				      // HWE
				      counts[thread_nloop][i].det->testHWE_Ca = testHWE(counts[thread_nloop][i].AA_Ca,counts[thread_nloop][i].AB_Ca, counts[thread_nloop][i].BB_Ca);
				      counts[thread_nloop][i].det->testHWE_Co = testHWE(counts[thread_nloop][i].AA_Co,counts[thread_nloop][i].AB_Co, counts[thread_nloop][i].BB_Co);

				      counts[thread_nloop][i].det->aCaN = 2 * counts[thread_nloop][i].AA_Ca + counts[thread_nloop][i].AB_Ca;
				      counts[thread_nloop][i].det->aCoN = 2 * counts[thread_nloop][i].AA_Co + counts[thread_nloop][i].AB_Co;
				      counts[thread_nloop][i].det->bCaN = 2 * counts[thread_nloop][i].BB_Ca + counts[thread_nloop][i].AB_Ca;
				      counts[thread_nloop][i].det->bCoN = 2 * counts[thread_nloop][i].BB_Co + counts[thread_nloop][i].AB_Co;
				    }

				  // Berechnen der Oddsratios, Konfidenzintervalle
				  if (singleMarkerTest == 1 || (singleMarkerTest==3 || singleMarkerTest >=5))

				    {
				      if(singleMarkerTest==5)
					{
					  counts[thread_nloop][i].det->aCaN = counts[thread_nloop][i].AA_Ca;
					  counts[thread_nloop][i].det->aCoN = counts[thread_nloop][i].AA_Co;
					  counts[thread_nloop][i].det->bCaN = counts[thread_nloop][i].BB_Ca + counts[thread_nloop][i].AB_Ca;
					  counts[thread_nloop][i].det->bCoN = counts[thread_nloop][i].BB_Co + counts[thread_nloop][i].AB_Co;
					}

				      if(singleMarkerTest==6)
					{
					  counts[thread_nloop][i].det->aCaN = counts[thread_nloop][i].AA_Ca + counts[thread_nloop][i].AB_Ca;
					  counts[thread_nloop][i].det->aCoN = counts[thread_nloop][i].AA_Co + counts[thread_nloop][i].AB_Co;
					  counts[thread_nloop][i].det->bCaN = counts[thread_nloop][i].BB_Ca;
					  counts[thread_nloop][i].det->bCoN = counts[thread_nloop][i].BB_Co;
					}

				      if(singleMarkerTest==7)
					{
					  counts[thread_nloop][i].det->aCaN = counts[thread_nloop][i].AB_Ca;
					  counts[thread_nloop][i].det->aCoN = counts[thread_nloop][i].AB_Co;
					  counts[thread_nloop][i].det->bCaN = counts[thread_nloop][i].AA_Ca + counts[thread_nloop][i].BB_Ca;
					  counts[thread_nloop][i].det->bCoN = counts[thread_nloop][i].AA_Co + counts[thread_nloop][i].BB_Co;
					}
				      if ((counts[thread_nloop][i].det->aCaN + counts[thread_nloop][i].det->bCaN) > 0)
					{
					  counts[thread_nloop][i].det->aCa = counts[thread_nloop][i].det->aCaN / (counts[thread_nloop][i].det->aCaN + counts[thread_nloop][i].det->bCaN);
					  counts[thread_nloop][i].det->bCa = 1 - counts[thread_nloop][i].det->aCa;
					}
				      else
					{
					  counts[thread_nloop][i].det->aCa = 0;
					  counts[thread_nloop][i].det->bCa = 0;
					}

				      if ((counts[thread_nloop][i].det->aCoN + counts[thread_nloop][i].det->bCoN) > 0)
					{
					  counts[thread_nloop][i].det->aCo = counts[thread_nloop][i].det->aCoN / (counts[thread_nloop][i].det->aCoN + counts[thread_nloop][i].det->bCoN);
					  counts[thread_nloop][i].det->bCo = 1 - counts[thread_nloop][i].det->aCo;
					}
				      else
					{
					  counts[thread_nloop][i].det->aCo = 0;
					  counts[thread_nloop][i].det->bCo = 0;
					}


				      if (counts[thread_nloop][i].det->aCaN > 0 && counts[thread_nloop][i].det->bCaN > 0 && counts[thread_nloop][i].det->aCoN > 0 && counts[thread_nloop][i].det->bCoN > 0)
					{
					  counts[thread_nloop][i].det->orA = (counts[thread_nloop][i].det->aCaN * counts[thread_nloop][i].det->bCoN) / (counts[thread_nloop][i].det->bCaN * counts[thread_nloop][i].det->aCoN);
					  counts[thread_nloop][i].det->orB = 1/counts[thread_nloop][i].det->orA;
					  counts[thread_nloop][i].det->lclA = exp(sqrt((1/counts[thread_nloop][i].det->aCaN) + (1/counts[thread_nloop][i].det->bCaN)+ (1/counts[thread_nloop][i].det->aCoN) + (1/counts[thread_nloop][i].det->bCoN)) * 1.96 * (-1))* counts[thread_nloop][i].det->orA;
					  counts[thread_nloop][i].det->rclA = exp(sqrt((1/counts[thread_nloop][i].det->aCaN) + (1/counts[thread_nloop][i].det->bCaN)+ (1/counts[thread_nloop][i].det->aCoN) + (1/counts[thread_nloop][i].det->bCoN)) * 1.96 * (1))* counts[thread_nloop][i].det->orA;
					  counts[thread_nloop][i].det->lclB = exp(sqrt((1/counts[thread_nloop][i].det->aCaN) + (1/counts[thread_nloop][i].det->bCaN)+ (1/counts[thread_nloop][i].det->aCoN) + (1/counts[thread_nloop][i].det->bCoN)) * 1.96 * (-1))* counts[thread_nloop][i].det->orB;
					  counts[thread_nloop][i].det->rclB = exp(sqrt((1/counts[thread_nloop][i].det->aCaN) + (1/counts[thread_nloop][i].det->bCaN)+ (1/counts[thread_nloop][i].det->aCoN) + (1/counts[thread_nloop][i].det->bCoN)) * 1.96 * (1))* counts[thread_nloop][i].det->orB;
					}
				      else
					{
					  counts[thread_nloop][i].det->lclA = -1;
					  counts[thread_nloop][i].det->rclA = -1;
					  counts[thread_nloop][i].det->lclB = -1;
					  counts[thread_nloop][i].det->rclB = -1;
					  counts[thread_nloop][i].det->orA = -1;
					  counts[thread_nloop][i].det->orB = -1;
					}
				    }
				  else if (singleMarkerTest == 2 || singleMarkerTest == 4)
				    {
				      if ((counts[thread_nloop][i].AA_Ca + counts[thread_nloop][i].AB_Ca + counts[thread_nloop][i].BB_Ca) > 0)
					{
					  counts[thread_nloop][i].det->aaCa = double(counts[thread_nloop][i].AA_Ca) / double(counts[thread_nloop][i].AA_Ca + counts[thread_nloop][i].AB_Ca + counts[thread_nloop][i].BB_Ca);
					  counts[thread_nloop][i].det->abCa = double(counts[thread_nloop][i].AB_Ca) / double(counts[thread_nloop][i].AA_Ca + counts[thread_nloop][i].AB_Ca + counts[thread_nloop][i].BB_Ca);
					  counts[thread_nloop][i].det->bbCa = double(counts[thread_nloop][i].BB_Ca) / double(counts[thread_nloop][i].AA_Ca + counts[thread_nloop][i].AB_Ca + counts[thread_nloop][i].BB_Ca);
					}
				      else
					{
					  counts[thread_nloop][i].det->aaCa = -1;
					  counts[thread_nloop][i].det->abCa = -1;
					  counts[thread_nloop][i].det->bbCa = -1;
					}

				      if 	((counts[thread_nloop][i].AA_Co + counts[thread_nloop][i].AB_Co + counts[thread_nloop][i].BB_Co) > 0)
					{
					  counts[thread_nloop][i].det->aaCo = double(counts[thread_nloop][i].AA_Co) / double(counts[thread_nloop][i].AA_Co + counts[thread_nloop][i].AB_Co + counts[thread_nloop][i].BB_Co);
					  counts[thread_nloop][i].det->abCo = double(counts[thread_nloop][i].AB_Co) / double(counts[thread_nloop][i].AA_Co + counts[thread_nloop][i].AB_Co + counts[thread_nloop][i].BB_Co);
					  counts[thread_nloop][i].det->bbCo = double(counts[thread_nloop][i].BB_Co) / double(counts[thread_nloop][i].AA_Co + counts[thread_nloop][i].AB_Co + counts[thread_nloop][i].BB_Co);
					}
				      else
					{
					  counts[thread_nloop][i].det->aaCo = -1;
					  counts[thread_nloop][i].det->abCo = -1;
					  counts[thread_nloop][i].det->bbCo = -1;
					}
				    }
				  map[i].p = newValue;
				  map[i].pmod = newValueCorr;
				} // n==0


			      counts[thread_nloop][i].pSingle=newValue;

			      if (newValue <= lastpSingle)
				{
				  //lastpSingle=newValue;
				  if(lastpSingle<1.1) // list is full
				    {
				      bestsinglemarker[maxposSingle].p = newValue;
				      bestsinglemarker[maxposSingle].pmod = newValueCorr;
				      bestsinglemarker[maxposSingle].nr = i;

				      //get new maxpos
				      lastpSingle=-0.1;
				      for(ii=0;ii<singletop;ii++)
					{
					  if(bestsinglemarker[ii].p>lastpSingle)
					    {
					      lastpSingle=bestsinglemarker[ii].p;
					      maxposSingle=ii;
					    }
					}
				    }
				  else if (maxposSingle < singletop -1 )
				    {
				      bestsinglemarker[maxposSingle].p = newValue;
				      bestsinglemarker[maxposSingle].pmod = newValueCorr;
				      bestsinglemarker[maxposSingle].nr = i;
				      maxposSingle++;
				    }
				  else if (maxposSingle == singletop -1 )
				    {
				      bestsinglemarker[maxposSingle].p = newValue;
				      bestsinglemarker[maxposSingle].pmod = newValueCorr;
				      bestsinglemarker[maxposSingle].nr = i;
				      //get new maxpos
				      lastpSingle=-0.1;
				      for(ii=0;ii<singletop;ii++)
					{
					  if(bestsinglemarker[ii].p>lastpSingle)
					    {
					      lastpSingle=bestsinglemarker[ii].p;
					      maxposSingle=ii;
					    }
					}
				    }
				} //newValue <=
			    } //end TIM
			  newValue = 1.1;
			} // end i

		      // update countsMC
		      if (nsim && n && !pathwayAnalysis && !markercombi2 && !markercombi3) {
			for (int z=0; z<singletop; z++) {
			  int *i = (int*)bsearch(&toplist[z].nr1, SNPMap, nsnpqc, sizeof(int), compare_int);
#pragma omp critical (COUNTS_MC)
			  if (!map[*i].done && counts[thread_nloop][*i].pSingle < toplist[z].p*(1+EPS)) {
			    countsMC[z]++;
			    if (0<sim_limit && sim_limit<countsMC[z]) {
			      countsMC[z] = (int)(0.5 + countsMC[z]*(double)nsim/(double)m);
			      map[*i].done = 1;
			    }
			  }
			}
#pragma omp critical (COUNTS_MC)
			m++;
		      }

		      qsortbestsinglemarker(&bestsinglemarker,0,singletop-1, &x1sm, &y1sm);
		      // exit(1);
		      //exit(1);
		      if (q == qstart )
			{
			  inflationfactor/=fulltests;
			}

		      if(n ==0 || plln)
			{

			  // Overlaplist
			  if(mWithGenetictop != 0 && mWithSingletop != 0)
			    {
			      if (singletop > genetictop && genetictop>0)
				{
				  overlaplist = (struct OVERLAPLIST *) calloc(singletop,sizeof(struct OVERLAPLIST));
				  if (!overlaplist)
				    {
				      errorfile << "memory allocation error in overlaplist\n";
				      logfile << "memory allocation error in overlaplist\n";
				      cout << "memory allocation error in overlaplist\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				}
			      else if(genetictop>0)
				{
				  overlaplist = (struct OVERLAPLIST *) calloc(genetictop,sizeof(struct OVERLAPLIST));
				  if (!overlaplist)
				    {
				      errorfile << "memory allocation error in overlaplist\n";
				      logfile << "memory allocation error in overlaplist\n";
				      cout << "memory allocation error in overlaplist\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				}
			      k =0;

			      if(genetictop>0)
				{
				  bestsinglemarker2 = (struct BESTSINGLEMARKER *) calloc(singletop,sizeof(struct BESTSINGLEMARKER));
				  if (!bestsinglemarker2)
				    {
				      errorfile << "memory allocation error in bestsinglemarker2\n";
				      logfile << "memory allocation error in bestsinglemarker2\n";
				      cout << "memory allocation error in bestsinglemarker2\n";
				      errorfile.close();logfile.close();exit(1);
				    }

				  geneticlist2 = (struct GENETICLIST *) calloc(genetictop,sizeof(struct GENETICLIST));
				  if (!geneticlist2)
				    {

				      errorfile << "memory allocation error in geneticlist2\n";
				      logfile << "memory allocation error in geneticlist2\n";
				      cout << "memory allocation error in geneticlist2\n";
				      errorfile.close();logfile.close();exit(1);
				    }

				  for(i=0;i<singletop;i++)
				    {
				      bestsinglemarker2[i]=bestsinglemarker[i];
				    }

				  for(i=0;i<genetictop;i++)
				    {
				      geneticlist2[i]=geneticlist[i];
				    }

				  qsortbestsinglemarkerNR(&bestsinglemarker2,0,singletop-1, &x1smNR, &y1smNR);
				  qsortgeneticlist(&geneticlist2,0,genetictop-1, &x1genetic, &y1genetic);


				  for (i=0; i < singletop; i++)
				    {
				      for (j=jstart; j < genetictop; j++)
					{
					  overlaplist[k].nr = 99;
					  if (bestsinglemarker2[i].nr == geneticlist2[j].nr)
					    {
					      overlaplist[k].nr = bestsinglemarker2[i].nr;
					      k++;
					      jstart=j+1;
					      break;
					    }
					  if(geneticlist2[j].nr>bestsinglemarker2[i].nr){break;}
					}//j
				    }//i
				}//if genetictop >0
			      overlaptop = k;
			    }
			} //n==0 || plln

		      if (n == 0 && q == qstart)
			{

			  if (ncasesqc == 0 && !qt)
			    {
			      errorfile << "\nAfter QC no cases are left.\n";
			      logfile << "\nAfter QC no cases are left.\n";
			      cout << "\nAfter QC no cases are left.\n";
			      errorfile.close();logfile.close();exit(1);
			    }
			  else if (ncontrolsqc == 0 && !qt && !caseOnly)
			    {
			      errorfile << "\nAfter QC no controls are left.\n";
			      logfile << "\nAfter QC no controls are left.\n";
			      cout << "\nAfter QC no controls are left.\n";
			      errorfile.close();logfile.close();exit(1);
			    }
			}

		      teststat = 1;

		      if(q==qstart) //TIM
			{
			  if (bestsinglemarker[singletop - 1].nr !=-1 && strcmp(map[bestsinglemarker[singletop - 1].nr].chr, "23") == 0)
			    {
			      helpstat = snpTestX(counts[thread_nloop][bestsinglemarker[singletop - 1].nr],teststat,&inflationfactor, &fulltests);
			      helpstat *= 1.4;
			    }
			  else if (bestsinglemarker[singletop - 1].nr !=-1)
			    {
			      helpstat = genotypTest(counts[thread_nloop][bestsinglemarker[singletop - 1].nr],teststat,&inflationfactor, &fulltests);
			    }
			  else if (strcmp(map[bestsinglemarker[0].nr].chr, "23") == 0)
			    {
			      helpstat = snpTestX(counts[thread_nloop][bestsinglemarker[singletop - 1].nr],teststat,&inflationfactor, &fulltests);
			      helpstat *= 1.4;
			    }
			  else if (bestsinglemarker[0].nr !=-1)
			    {
			      helpstat = genotypTest(counts[thread_nloop][bestsinglemarker[0].nr],teststat,&inflationfactor, &fulltests);
			    }
			  else {
			    helpstat = 0;
			  }
			  teststat=0;


			} //endif qstart TIM

		      if (SM_showresults && n == 0 && q <= 1 && r <=1)
			{
			  //							cout << "\n";
			  //							cout << "singletop: " << singletop << "\n";
			  //							cout << "genetictop: " << genetictop << "\n";
			  //							cout << "overlaptop: " << overlaptop << "\n";
			  logfile << "\n";
			  logfile << "singletop: " << singletop << "\n";
			  logfile << "genetictop: " << genetictop << "\n";
			  logfile << "overlaptop: " << overlaptop << "\n";

			  if (fulltests > 0)
			    {
			      cout    << "inflation factor: " << inflationfactor << endl;
			      logfile << "inflation factor: " << inflationfactor << endl;
			    }

			  // 10 besten P-Werte werden ausgeben
			  cout << "\ntop single-marker p-values:\n";
			  cout << "Chr rs_No Position P_Single-marker P_Corr\n";

			  singlemarkerTop.open(singlemarkerTopfile.c_str(), ios::out);
			  singlemarker.open(singlemarkerfile.c_str(), ios::out);

			  if (singleMarkerTest == 1)
			    {
			      singlemarkerTop <<"No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tA_Ca_N\tB_Ca_N\tA_Co_N\tB_Co_N\tA_Ca\tB_Ca\tA_Co\tB_Co\tOR_A\tLCL_A\tRCL_A\tOR_B\tLCL_B\tRCL_B\n";

			      singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tA_Ca_N\tB_Ca_N\tA_Co_N\tB_Co_N\tA_Ca\tB_Ca\tA_Co\tB_Co\tOR_A\tLCL_A\tRCL_A\tOR_B\tLCL_B\tRCL_B\n";
			    }
			  else if (singleMarkerTest == 2)
			    {
			      singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAA_Ca_N\tAB_Ca_N\tBB_Ca_N\tAA_Co_N\tAB_Co_N\tBB_Co_N\tAA_Ca\tAB_Ca\tBB_Ca\tAA_Co\tAB_Co\tBB_Co\n";

			      singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAA_Ca_N\tAB_Ca_N\tBB_Ca_N\tAA_Co_N\tAB_Co_N\tBB_Co_N\tAA_Ca\tAB_Ca\tBB_Ca\tAA_Co\tAB_Co\tBB_Co\n";
			    }
			  else if (singleMarkerTest==3)
			    {
			      if (qt == 1)
				{
				  //singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tN_minor\tN_major\tFreq_minor\tFreq_major\tbeta1\tse1\tLCL1\tRLC1\n";

				  singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tN_minor\tN_major\tFreq_minor\tFreq_major\tbeta1\tse1\tLCL1\tRLC1";
				  if(covariancematrix)
				    {
				      for(int m=0; m<df_L1_Single; m++)
					{
					  for(int s=m; s<df_L1_Single; s++)
					    {
					      singlemarkerTop << "\tsigma1[" << m <<"]" << "[" << s << "]";
					    }
					}
				    }

				  singlemarkerTop << "\n";


				  //singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tN_minor\tN_major\tFreq_minor\tFreq_major\tbeta1\tse1\tLCL1\tRLC1\n";

				  singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tN_minor\tN_major\tFreq_minor\tFreq_major\tbeta1\tse1\tLCL1\tRLC1";



				  if(covariancematrix)
				    {
				      for(int m=0; m<df_L1_Single; m++)
					{
					  for(int s=m; s<df_L1_Single; s++)
					    {
					      singlemarker << "\tsigma1[" << m <<"]" << "[" << s << "]";
					    }
					}
				    }
				  singlemarker << "\n";
				}
			      else
				{
				  singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tA_Ca_N\tB_Ca_N\tA_Co_N\tB_Co_N\tA_Ca\tB_Ca\tA_Co\tB_Co\tOR_A\tLCL_A\tRCL_A\tOR_B\tLCL_B\tRCL_B\tbeta1\tse1\tOR1\tLCL1\tRLC1";
				  if(covariancematrix)
				    {
				      for(int m=0; m<df_L1_Single; m++)
					{
					  for(int s=m; s<df_L1_Single; s++)
					    {
					      singlemarkerTop << "\tsigma1[" << m <<"]" << "[" << s << "]";
					    }
					}
				    }

				  singlemarkerTop << "\n";


				  //singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tA_Ca_N\tB_Ca_N\tA_Co_N\tB_Co_N\tA_Ca\tB_Ca\tA_Co\tB_Co\tOR_A\tLCL_A\tRCL_A\tOR_B\tLCL_B\tRCL_B\tbeta1\tse1\tOR1\tLCL1\tRLC1\n";

				  singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tA_Ca_N\tB_Ca_N\tA_Co_N\tB_Co_N\tA_Ca\tB_Ca\tA_Co\tB_Co\tOR_A\tLCL_A\tRCL_A\tOR_B\tLCL_B\tRCL_B\tbeta1\tse1\tOR1\tLCL1\tRLC1";


				  if(covariancematrix)
				    {
				      for(int m=0; m<df_L1_Single; m++)
					{
					  for(int s=m; s<df_L1_Single; s++)
					    {
					      singlemarker << "\tsigma1[" << m <<"]" << "[" << s << "]";
					    }
					}
				    }
				  singlemarker << "\n";

				}
			    }
			  else if (singleMarkerTest == 4)
			    {
			      if (qt == 1)
				{
				  singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tAA_N\tAB_N\tBB_N\tAA_Freq\tAB_Freq\tBB_Freq\tbeta1\tse1\tLCL1\tRLC1\tbeta1D\tse1D\tLCL1D\tRLC1D";



				  if(covariancematrix)
				    {
				      for(int m=0; m<df_L1_Single; m++)
					{
					  for(int s=m; s<df_L1_Single; s++)
					    {
					      singlemarkerTop << "\tsigma1[" << m <<"]" << "[" << s << "]";
					    }
					}
				    }
				  singlemarkerTop << "\n";


				  singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tAA_N\tAB_N\tBB_N\tAA_Freq\tAB_Freq\tBB_Freq\tbeta1\tse1\tLCL1\tRLC1\tbeta1D\tse1D\tLCL1D\tRLC1D";


				  if(covariancematrix)
				    {
				      for(int m=0; m<df_L1_Single; m++)
					{
					  for(int s=m; s<df_L1_Single; s++)
					    {
					      singlemarker << "\tsigma1[" << m <<"]" << "[" << s << "]";
					    }
					}
				    }
				  singlemarker << "\n";

				}
			      else
				{
				  //singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAA_Ca_N\tAB_Ca_N\tBB_Ca_N\tAA_Co_N\tAB_Co_N\tBB_Co_N\tAA_Ca\tAB_Ca\tBB_Ca\tAA_Co\tAB_Co\tBB_Co\tbeta1\tse1\tOR1\tLCL1\tRLC1\tbeta1D\tse1D\tOR1D\tLCL1D\tRLC1D\n";

				  singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAA_Ca_N\tAB_Ca_N\tBB_Ca_N\tAA_Co_N\tAB_Co_N\tBB_Co_N\tAA_Ca\tAB_Ca\tBB_Ca\tAA_Co\tAB_Co\tBB_Co\tbeta1\tse1\tOR1\tLCL1\tRLC1\tbeta1D\tse1D\tOR1D\tLCL1D\tRLC1D";


				  if(covariancematrix)
				    {
				      for(int m=0; m<df_L1_Single; m++)
					{
					  for(int s=m; s<df_L1_Single; s++)
					    {
					      singlemarkerTop << "\tsigma1[" << m <<"]" << "[" << s << "]";
					    }
					}
				    }

				  singlemarkerTop << "\n";


				  //singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAA_Ca_N\tAB_Ca_N\tBB_Ca_N\tAA_Co_N\tAB_Co_N\tBB_Co_N\tAA_Ca\tAB_Ca\tBB_Ca\tAA_Co\tAB_Co\tBB_Co\tbeta1\tse1\tOR1\tLCL1\tRLC1\tbeta1D\tse1D\tOR1D\tLCL1D\tRLC1D\n";

				  singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAA_Ca_N\tAB_Ca_N\tBB_Ca_N\tAA_Co_N\tAB_Co_N\tBB_Co_N\tAA_Ca\tAB_Ca\tBB_Ca\tAA_Co\tAB_Co\tBB_Co\tbeta1\tse1\tOR1\tLCL1\tRLC1\tbeta1D\tse1D\tOR1D\tLCL1D\tRLC1D";


				  if(covariancematrix)
				    {
				      for(int m=0; m<df_L1_Single; m++)
					{
					  for(int s=m; s<df_L1_Single; s++)
					    {
					      singlemarker << "\tsigma1[" << m <<"]" << "[" << s << "]";
					    }
					}
				    }
				  singlemarker << "\n";
				}

			    }
			  else if(singleMarkerTest >=5)
			    {
			      if (qt == 1)
				{
				  if(singleMarkerTest==5)
				    {
				      //singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tAA_N\tnotAA_N\tFreq_AA\tFreq_notAA\tbeta1\tse1\tLCL1\tRLC1\n";
				      singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tAA_N\tnotAA_N\tFreq_AA\tFreq_notAA\tbeta1\tse1\tLCL1\tRLC1";

				      if(covariancematrix)
					{
					  for(int m=0; m<df_L1_Single;m++)
					    {
					      for(int s=m; s<df_L1_Single;s++)
						{
						  singlemarkerTop << "\tsigma1[" << m <<"]" << "[" << s << "]";
						}
					    }
					}
				      singlemarkerTop << "\n";

				      //singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tAA_N\tnotAA_N\tFreq_AA\tFreq_notAA\tbeta1\tse1\tLCL1\tRLC1\n";
				      singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tAA_N\tnotAA_N\tFreq_AA\tFreq_notAA\tbeta1\tse1\tLCL1\tRLC1";
				      if(covariancematrix)
					{
					  for(int m=0; m<df_L1_Single;m++)
					    {
					      for(int s=m; s<df_L1_Single;s++)
						{
						  singlemarker << "\tsigma1[" << m <<"]" << "[" << s << "]";
						}
					    }
					}
				      singlemarker << "\n";
				    }
				  if(singleMarkerTest==6)
				    {
				      //singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tAA/AB_N\tBB_N\tFreq_AA/AB\tFreq_BB\tbeta1\tse1\tLCL1\tRLC1\n";
				      singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tAA/AB_N\tBB_N\tFreq_AA/AB\tFreq_BB\tbeta1\tse1\tLCL1\tRLC1";
				      if(covariancematrix)
					{
					  for(int m=0; m<df_L1_Single;m++)
					    {
					      for(int s=m; s<df_L1_Single;s++)
						{
						  singlemarkerTop << "\tsigma1[" << m <<"]" << "[" << s << "]";
						}
					    }
					}
				      singlemarkerTop << "\n";

				      //singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tAA/AB_N\tBB_N\tFreq_AA/AB\tFreq_BB\tbeta1\tse1\tLCL1\tRLC1\n";
				      singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tAA/AB_N\tBB_N\tFreq_AA/AB\tFreq_BB\tbeta1\tse1\tLCL1\tRLC1";
				      if(covariancematrix)
					{
					  for(int m=0; m<df_L1_Single;m++)
					    {
					      for(int s=m; s<df_L1_Single;s++)
						{
						  singlemarker << "\tsigma1[" << m <<"]" << "[" << s << "]";
						}
					    }
					}
				      singlemarker << "\n";
				    }
				  if(singleMarkerTest==7)
				    {
				      //singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tAB_N\tnotAB_N\tFreq_AB\tFreq_notAB\tbeta1\tse1\tLCL1\tRLC1\tGeno1\tGeno2\n";
				      singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tAB_N\tnotAB_N\tFreq_AB\tFreq_notAB\tbeta1\tse1\tLCL1\tRLC1\tGeno1\tGeno2";
				      if(covariancematrix)
					{
					  for(int m=0; m<df_L1_Single;m++)
					    {
					      for(int s=m; s<df_L1_Single;s++)
						{
						  singlemarkerTop << "\tsigma1[" << m <<"]" << "[" << s << "]";
						}
					    }
					}

				      singlemarkerTop << "\n";

				      //singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tAB_N\tnotAB_N\tFreq_AB\tFreq_notAB\tbeta1\tse1\tLCL1\tRLC1\tGeno1\tGeno2\n";
				      singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF\tSNP_MR\tHWE\tP_Single-marker\tP_Corr\tAB_N\tnotAB_N\tFreq_AB\tFreq_notAB\tbeta1\tse1\tLCL1\tRLC1\tGeno1\tGeno2";
				      if(covariancematrix)
					{
					  for(int m=0; m<df_L1_Single;m++)
					    {
					      for(int s=m; s<df_L1_Single;s++)
						{
						  singlemarker << "\tsigma1[" << m <<"]" << "[" << s << "]";
						}
					    }
					}
				      singlemarker << "\n";
				    }
				}
			      else
				{
				  if(singleMarkerTest==5)
				    {
				      //singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAA_Ca_N\tnotAA_Ca_N\tAA_Co_N\tAA_Co_N\tAA_Ca\tnotAA_Ca\tAA_Co\tnotAA_Co\tOR_AA\tLCL_AA\tRCL_AA\tOR_notAA\tLCL_notAA\tRCL_notAA\tbeta1\tse1\tOR1\tLCL1\tRLC1\n";
				      singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAA_Ca_N\tnotAA_Ca_N\tAA_Co_N\tAA_Co_N\tAA_Ca\tnotAA_Ca\tAA_Co\tnotAA_Co\tOR_AA\tLCL_AA\tRCL_AA\tOR_notAA\tLCL_notAA\tRCL_notAA\tbeta1\tse1\tOR1\tLCL1\tRLC1";
				      if(covariancematrix)
					{
					  for(int m=0; m<df_L1_Single;m++)
					    {
					      for(int s=m; s<df_L1_Single;s++)
						{
						  singlemarkerTop << "\tsigma1[" << m <<"]" << "[" << s << "]";
						}
					    }
					}
				      singlemarkerTop << "\n";


				      //singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAA_Ca_N\tnotAA_Ca_N\tAA_Co_N\tnotAA_Co_N\tAA_Ca\tnotAA_Ca\tAA_Co\tnotAA_Co\tOR_AA\tLCL_AA\tRCL_AA\tOR_notAA\tLCL_notAA\tRCL_notAA\tbeta1\tse1\tOR1\tLCL1\tRLC1\n";
				      singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAA_Ca_N\tnotAA_Ca_N\tAA_Co_N\tnotAA_Co_N\tAA_Ca\tnotAA_Ca\tAA_Co\tnotAA_Co\tOR_AA\tLCL_AA\tRCL_AA\tOR_notAA\tLCL_notAA\tRCL_notAA\tbeta1\tse1\tOR1\tLCL1\tRLC1";
				      if(covariancematrix)
					{
					  for(int m=0; m<df_L1_Single;m++)
					    {
					      for(int s=m; s<df_L1_Single;s++)
						{
						  singlemarker << "\tsigma1[" << m <<"]" << "[" << s << "]";
						}
					    }
					}
				      singlemarker << "\n";
				    }

				  if(singleMarkerTest==6)
				    {
				      //singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAA/AB_Ca_N\tBB_Ca_N\tAA/AB_Co_N\tAA/AB_Co_N\tAA/AB_Ca\tBB_Ca\tAA/AB_Co\tBB_Co\tOR_AA/AB\tLCL_AA/AB\tRCL_AA/AB\tOR_BB\tLCL_BB\tRCL_BB\tbeta1\tse1\tOR1\tLCL1\tRLC1\n";
				      singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAA/AB_Ca_N\tBB_Ca_N\tAA/AB_Co_N\tAA/AB_Co_N\tAA/AB_Ca\tBB_Ca\tAA/AB_Co\tBB_Co\tOR_AA/AB\tLCL_AA/AB\tRCL_AA/AB\tOR_BB\tLCL_BB\tRCL_BB\tbeta1\tse1\tOR1\tLCL1\tRLC1";
				      if(covariancematrix)
					{
					  for(int m=0; m<df_L1_Single;m++)
					    {
					      for(int s=m; s<df_L1_Single;s++)
						{
						  singlemarkerTop << "\tsigma1[" << m <<"]" << "[" << s << "]";
						}
					    }
					}
				      singlemarkerTop << "\n";

				      //singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAA/AB_Ca_N\tBB_Ca_N\tAA/AB_Co_N\tBB_Co_N\tAA/AB_Ca\tBB_Ca\tAA/AB_Co\tBB_Co\tOR_AA/AB\tLCL_AA/AB\tRCL_AA/AB\tOR_BB\tLCL_BB\tRCL_BB\tbeta1\tse1\tOR1\tLCL1\tRLC1\n";
				      singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAA/AB_Ca_N\tBB_Ca_N\tAA/AB_Co_N\tBB_Co_N\tAA/AB_Ca\tBB_Ca\tAA/AB_Co\tBB_Co\tOR_AA/AB\tLCL_AA/AB\tRCL_AA/AB\tOR_BB\tLCL_BB\tRCL_BB\tbeta1\tse1\tOR1\tLCL1\tRLC1";
				      if(covariancematrix)
					{
					  for(int m=0; m<df_L1_Single;m++)
					    {
					      for(int s=m; s<df_L1_Single;s++)
						{
						  singlemarker << "\tsigma1[" << m <<"]" << "[" << s << "]";
						}
					    }
					}
				      singlemarker << "\n";
				    }

				  if(singleMarkerTest==7)
				    {
				      //singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAB_Ca_N\tnotAB_Ca_N\tAB_Co_N\tAB_Co_N\tAB_Ca\tnotAB_Ca\tAB_Co\tnotAB_Co\tOR_AB\tLCL_AB\tRCL_AB\tOR_notAB\tLCL_notAB\tRCL_notAB\tbeta1\tse1\tOR1\tLCL1\tRLC1\tGeno1\tGeno2\n";
				      singlemarkerTop << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAB_Ca_N\tnotAB_Ca_N\tAB_Co_N\tAB_Co_N\tAB_Ca\tnotAB_Ca\tAB_Co\tnotAB_Co\tOR_AB\tLCL_AB\tRCL_AB\tOR_notAB\tLCL_notAB\tRCL_notAB\tbeta1\tse1\tOR1\tLCL1\tRLC1\tGeno1\tGeno2";
				      if(covariancematrix)
					{
					  for(int m=0; m<df_L1_Single;m++)
					    {
					      for(int s=m; s<df_L1_Single;s++)
						{
						  singlemarkerTop << "\tsigma1[" << m <<"]" << "[" << s << "]";
						}
					    }
					}
				      singlemarkerTop << "\n";

				      //singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAB_Ca_N\tnotAB_Ca_N\tAB_Co_N\tnotAB_Co_N\tAB_Ca\tnotAB_Ca\tAB_Co\tnotAB_Co\tOR_AB\tLCL_AB\tRCL_AB\tOR_notAB\tLCL_notAB\tRCL_notAB\tbeta1\tse1\tOR1\tLCL1\tRLC1\tGeno1\tGeno2\n";
				      singlemarker << "No\tChr\trs_No\tPosition\tGene\tminor\tmajor\tMAF(All/Co/Ca)\tSNP_MR\tHWE_Ca\tHWE_Co\tP_Single-marker\tP_Corr\tAB_Ca_N\tnotAB_Ca_N\tAB_Co_N\tnotAB_Co_N\tAB_Ca\tnotAB_Ca\tAB_Co\tnotAB_Co\tOR_AB\tLCL_AB\tRCL_AB\tOR_notAB\tLCL_notAB\tRCL_notAB\tbeta1\tse1\tOR1\tLCL1\tRLC1\tGeno1\tGeno2";
				      if(covariancematrix)
					{
					  for(int m=0; m<df_L1_Single;m++)
					    {
					      for(int s=m; s<df_L1_Single;s++)
						{
						  singlemarker << "\tsigma1[" << m <<"]" << "[" << s << "]";
						}
					    }
					}
				      singlemarker << "\n";
				    }
				}
			    }

			  logfile << "\ntop single-marker p-values:\n";
			  logfile << "Chr rs_No Position P_Single-marker P_Corr\n";

			  for (k = 0; k < singletop; k++)
			    {
			      if (bestsinglemarker[k].nr < 0)
				{
				  cout    << k << ": Warning: bestsinglemarker[k].nr < 0!\nexiting!";
				  logfile << k << ": Warning: bestsinglemarker[k].nr < 0!\nexiting!";
				  exit(1);
				}
			      if (k < 10)
				{
				  cout << map[bestsinglemarker[k].nr].chr << " " << map[bestsinglemarker[k].nr].rs << " " << map[bestsinglemarker[k].nr].pos << " "
				       << bestsinglemarker[k].p << " " << bestsinglemarker[k].pmod << "\n";

				  logfile << map[bestsinglemarker[k].nr].chr << " " << map[bestsinglemarker[k].nr].rs << " " << map[bestsinglemarker[k].nr].pos << " "
					  << bestsinglemarker[k].p << " " << bestsinglemarker[k].pmod << "\n";
				}

			      if (singleMarkerTest == 1)
				{
				  singlemarkerTop << k + 1 << "\t" << map[bestsinglemarker[k].nr].chr << "\t" << map[bestsinglemarker[k].nr].rs << "\t"
						  << map[bestsinglemarker[k].nr].pos << "\t" << map[bestsinglemarker[k].nr].gene << "\t" << codesA[bestsinglemarker[k].nr].a1 << "\t"
						  << codesA[bestsinglemarker[k].nr].a2 << "\t" << map[bestsinglemarker[k].nr].maf <<"/"<< map[bestsinglemarker[k].nr].controlmaf <<"/"<< map[bestsinglemarker[k].nr].casemaf << "\t"
						  << map[bestsinglemarker[k].nr].missingrate << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->testHWE_Ca << "\t"
						  << counts[thread_nloop][bestsinglemarker[k].nr].det->testHWE_Co << "\t" << bestsinglemarker[k].p << "\t" << bestsinglemarker[k].pmod << "\t"
						  << counts[thread_nloop][bestsinglemarker[k].nr].det->aCaN << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->bCaN << "\t"
						  << counts[thread_nloop][bestsinglemarker[k].nr].det->aCoN << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->bCoN << "\t"
						  << counts[thread_nloop][bestsinglemarker[k].nr].det->aCa << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->bCa << "\t"
						  << counts[thread_nloop][bestsinglemarker[k].nr].det->aCo << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->bCo << "\t"
						  << counts[thread_nloop][bestsinglemarker[k].nr].det->orA << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->lclA << "\t"
						  << counts[thread_nloop][bestsinglemarker[k].nr].det->rclA << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->orB << "\t"
						  << counts[thread_nloop][bestsinglemarker[k].nr].det->lclB << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->rclB << "\n";
				}
			      else if (singleMarkerTest == 2)
				{
				  singlemarkerTop << k + 1 << "\t" << map[bestsinglemarker[k].nr].chr << "\t" << map[bestsinglemarker[k].nr].rs << "\t"
						  << map[bestsinglemarker[k].nr].pos << "\t" << map[bestsinglemarker[k].nr].gene << "\t" << codesA[bestsinglemarker[k].nr].a1 << "\t"
						  << codesA[bestsinglemarker[k].nr].a2 << "\t" << map[bestsinglemarker[k].nr].maf <<"/"<< map[bestsinglemarker[k].nr].controlmaf <<"/"<< map[bestsinglemarker[k].nr].casemaf << "\t"
						  << map[bestsinglemarker[k].nr].missingrate << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->testHWE_Ca << "\t"
						  << counts[thread_nloop][bestsinglemarker[k].nr].det->testHWE_Co << "\t" << bestsinglemarker[k].p << "\t" << bestsinglemarker[k].pmod << "\t"
						  << counts[thread_nloop][bestsinglemarker[k].nr].AA_Ca << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].AB_Ca << "\t"
						  << counts[thread_nloop][bestsinglemarker[k].nr].BB_Ca << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].AA_Co << "\t"
						  << counts[thread_nloop][bestsinglemarker[k].nr].AB_Co << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].BB_Co << "\t"
						  << counts[thread_nloop][bestsinglemarker[k].nr].det->aaCa << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->abCa << "\t"
						  << counts[thread_nloop][bestsinglemarker[k].nr].det->bbCa << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->aaCo << "\t"
						  << counts[thread_nloop][bestsinglemarker[k].nr].det->abCo << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->bbCo << "\n";
				}
			      else if (singleMarkerTest==3 || singleMarkerTest >=5)
				{
				  if (qt == 1)
				    {
				      singlemarkerTop << k + 1 << "\t" << map[bestsinglemarker[k].nr].chr << "\t"
						      << map[bestsinglemarker[k].nr].rs << "\t" << map[bestsinglemarker[k].nr].pos << "\t" << map[bestsinglemarker[k].nr].gene << "\t"
						      << codesA[bestsinglemarker[k].nr].a1 << "\t" << codesA[bestsinglemarker[k].nr].a2 << "\t" << map[bestsinglemarker[k].nr].maf <<"/"<< map[bestsinglemarker[k].nr].controlmaf <<"/"<< map[bestsinglemarker[k].nr].casemaf<< "\t" << map[bestsinglemarker[k].nr].missingrate << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].det->testHWE_Ca << "\t" << bestsinglemarker[k].p << "\t"
						      << bestsinglemarker[k].pmod << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->aCaN << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].det->bCaN << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->aCa << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].det->bCa << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).b[1] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).betaNew_se[1] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).betaNew_lcl[1] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).betaNew_rcl[1];
				      if(singleMarkerTest ==7) {singlemarkerTop << "\tH\th";}
				      //if(singleMarkerTest ==3)
				      //{
				      if(covariancematrix)
					{
					  for(int m=0; m <(df_L1_Single*(df_L1_Single+1))/2; m++)
					    {
					      singlemarkerTop << "\t" << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).sigma1[SigmaAuxVector_Single[m]];
					      //singlemarkerTop << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].result1.sigma1[m]; // output of sigma1 elements in a singlemarkerTop file
					    }
					}
				      //}
				      singlemarkerTop << "\n";
				    }
				  else
				    {
				      singlemarkerTop << k + 1 << "\t" << map[bestsinglemarker[k].nr].chr << "\t"
						      << map[bestsinglemarker[k].nr].rs << "\t" << map[bestsinglemarker[k].nr].pos << "\t" << map[bestsinglemarker[k].nr].gene << "\t"
						      << codesA[bestsinglemarker[k].nr].a1 << "\t" << codesA[bestsinglemarker[k].nr].a2 << "\t" << map[bestsinglemarker[k].nr].maf <<"/"<< map[bestsinglemarker[k].nr].controlmaf <<"/"<< map[bestsinglemarker[k].nr].casemaf<< "\t" << map[bestsinglemarker[k].nr].missingrate << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].det->testHWE_Ca << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->testHWE_Co << "\t"
						      << bestsinglemarker[k].p << "\t"<< bestsinglemarker[k].pmod << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->aCaN << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].det->bCaN << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->aCoN << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].det->bCoN << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->aCa << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].det->bCa << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->aCo << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].det->bCo << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->orA << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].det->lclA << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->rclA << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].det->orB << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->lclB << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].det->rclB << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).b[1] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).betaNew_se[1] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).oddsRatio[1]  << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).lcloddsRatio[1] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).rcloddsRatio[1];
				      if(singleMarkerTest ==7) {singlemarkerTop << "\tH\th";}
				      //if(singleMarkerTest ==3)
				      //{
				      if(covariancematrix)
					{
					  for(int m=0; m < (df_L1_Single*(df_L1_Single+1))/2; m++)
					    {
					      singlemarkerTop << "\t" << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).sigma1[SigmaAuxVector_Single[m]];
					    }
					}
				      //}
				      singlemarkerTop << "\n";
				    }
				}
			      else if (singleMarkerTest == 4 )
				{
				  if (qt == 1)
				    {
				      singlemarkerTop << k + 1 << "\t" << map[bestsinglemarker[k].nr].chr << "\t" << map[bestsinglemarker[k].nr].rs << "\t" << map[bestsinglemarker[k].nr].pos << "\t" << map[bestsinglemarker[k].nr].gene << "\t" << codesA[bestsinglemarker[k].nr].a1 << "\t"
						      << codesA[bestsinglemarker[k].nr].a2 << "\t" << map[bestsinglemarker[k].nr].maf <<"/"<< map[bestsinglemarker[k].nr].controlmaf <<"/"<< map[bestsinglemarker[k].nr].casemaf << "\t" << map[bestsinglemarker[k].nr].missingrate << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->testHWE_Ca << "\t"
						      << map[bestsinglemarker[k].nr].p << "\t" << map[bestsinglemarker[k].nr].pmod << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].AA_Ca << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].AB_Ca << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].BB_Ca << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].det->aaCa << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->abCa << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->bbCa << "\t" << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).b[1] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).betaNew_se[1] << "\t"  << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).betaNew_lcl[1] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).betaNew_rcl[1] << "\t" << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).b[2] << "\t" << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).betaNew_se[2] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).betaNew_lcl[2] << "\t" << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).betaNew_rcl[2];

				      if(covariancematrix)
					{
					  for(int m=0; m < (df_L1_Single*(df_L1_Single+1))/2; m++)
					    {
					      singlemarkerTop << "\t" << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).sigma1[SigmaAuxVector_Single[m]];
					    }

					}
				      singlemarkerTop << "\n";
				    }
				  else
				    {
				      singlemarkerTop << k + 1 << "\t" << map[bestsinglemarker[k].nr].chr << "\t"
						      << map[bestsinglemarker[k].nr].rs << "\t" << map[bestsinglemarker[k].nr].pos << "\t" << map[bestsinglemarker[k].nr].gene << "\t"
						      << codesA[bestsinglemarker[k].nr].a1 << "\t" << codesA[bestsinglemarker[k].nr].a2 << "\t" << map[bestsinglemarker[k].nr].maf <<"/"<< map[bestsinglemarker[k].nr].controlmaf <<"/"<< map[bestsinglemarker[k].nr].casemaf << "\t" << map[bestsinglemarker[k].nr].missingrate << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].det->testHWE_Ca << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->testHWE_Co << "\t"
						      << bestsinglemarker[k].p << "\t" << bestsinglemarker[k].pmod << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].AA_Ca << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].AB_Ca << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].BB_Ca << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].AA_Co << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].AB_Co << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].BB_Co << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].det->aaCa << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->abCa << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].det->bbCa << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->aaCo << "\t"
						      << counts[thread_nloop][bestsinglemarker[k].nr].det->abCo << "\t" << counts[thread_nloop][bestsinglemarker[k].nr].det->bbCo << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).b[1] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).betaNew_se[1] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).oddsRatio[1]  << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).lcloddsRatio[1] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).rcloddsRatio[1] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).b[2] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).betaNew_se[2] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).oddsRatio[2] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).lcloddsRatio[2] << "\t"
						      << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).rcloddsRatio[2];

			  if(covariancematrix)
			    {
			      for(int m=0; m < (df_L1_Single*(df_L1_Single+1))/2; m++)
				{
				  singlemarkerTop << "\t" << (*counts[thread_nloop][bestsinglemarker[k].nr].result1).sigma1[SigmaAuxVector_Single[m]];
				}
			    }

			  singlemarkerTop << "\n";
			}
		    }
		}

			  // Ausgabe Singlemarker
			  for (i = 0; i < nlinestped; i++)
			    {
			      if(map[i].analysis_in==1)
				{
				  if (singleMarkerTest == 1)
				    {
				      singlemarker << i + 1 << "\t" << map[i].chr << "\t" << map[i].rs << "\t" << map[i].pos << "\t" << map[i].gene << "\t" << codesA[i].a1 << "\t" << codesA[i].a2 << "\t" << map[i].maf <<"/"<< map[i].controlmaf <<"/"<< map[i].casemaf << "\t" << map[i].missingrate << "\t" << counts[thread_nloop][i].det->testHWE_Ca << "\t" << counts[thread_nloop][i].det->testHWE_Co << "\t" << map[i].p << "\t" << map[i].pmod << "\t" << counts[thread_nloop][i].det->aCaN << "\t" << counts[thread_nloop][i].det->bCaN << "\t" << counts[thread_nloop][i].det->aCoN << "\t" << counts[thread_nloop][i].det->bCoN << "\t" << counts[thread_nloop][i].det->aCa << "\t" << counts[thread_nloop][i].det->bCa << "\t" << counts[thread_nloop][i].det->aCo << "\t" << counts[thread_nloop][i].det->bCo << "\t" << counts[thread_nloop][i].det->orA << "\t" << counts[thread_nloop][i].det->lclA << "\t" << counts[thread_nloop][i].det->rclA << "\t" << counts[thread_nloop][i].det->orB << "\t" << counts[thread_nloop][i].det->lclB << "\t" << counts[thread_nloop][i].det->rclB << "\n";
				    }
				  else if (singleMarkerTest == 2)
				    {
				      singlemarker << i + 1 << "\t" << map[i].chr << "\t" << map[i].rs << "\t" << map[i].pos << "\t" << map[i].gene << "\t" << codesA[i].a1 << "\t" << codesA[i].a2 << "\t" << map[i].maf <<"/"<< map[i].controlmaf <<"/"<< map[i].casemaf << "\t" << map[i].missingrate << "\t" << counts[thread_nloop][i].det->testHWE_Ca << "\t" << counts[thread_nloop][i].det->testHWE_Co << "\t" << map[i].p << "\t" << map[i].pmod << "\t" << counts[thread_nloop][i].AA_Ca << "\t" << counts[thread_nloop][i].AB_Ca << "\t" << counts[thread_nloop][i].BB_Ca << "\t" << counts[thread_nloop][i].AA_Co << "\t" << counts[thread_nloop][i].AB_Co << "\t" << counts[thread_nloop][i].BB_Co << "\t" << counts[thread_nloop][i].det->aaCa << "\t" << counts[thread_nloop][i].det->abCa << "\t" << counts[thread_nloop][i].det->bbCa << "\t" << counts[thread_nloop][i].det->aaCo << "\t" << counts[thread_nloop][i].det->abCo << "\t" << counts[thread_nloop][i].det->bbCo << "\n";
				    }
				  else if (singleMarkerTest==3 || singleMarkerTest >=5)
				    {
				      if (qt == 1)
					{
					  singlemarker << i + 1 << "\t" << map[i].chr << "\t" << map[i].rs << "\t" << map[i].pos << "\t" << map[i].gene << "\t" << codesA[i].a1 << "\t" << codesA[i].a2 << "\t" << map[i].maf <<"/"<< map[i].controlmaf <<"/"<< map[i].casemaf << "\t" << map[i].missingrate << "\t" << counts[thread_nloop][i].det->testHWE_Ca << "\t"<< map[i].p << "\t" << map[i].pmod << "\t" << counts[thread_nloop][i].det->aCaN << "\t" << counts[thread_nloop][i].det->bCaN << "\t" << counts[thread_nloop][i].det->aCa << "\t" << counts[thread_nloop][i].det->bCa << "\t" << (*counts[thread_nloop][i].result1).b[1] << "\t" << (*counts[thread_nloop][i].result1).betaNew_se[1] << "\t" << (*counts[thread_nloop][i].result1).betaNew_lcl[1] << "\t" << (*counts[thread_nloop][i].result1).betaNew_rcl[1];
					  if(singleMarkerTest ==7) {singlemarker << "\tH\th";}
					  //if(singleMarkerTest ==3)
					  //{
					  if(covariancematrix)
					    {
					      for(int m=0; m < (df_L1_Single*(df_L1_Single+1))/2; m++)
						{
						  singlemarker << "\t" << (*counts[thread_nloop][i].result1).sigma1[SigmaAuxVector_Single[m]];
						}
					    }
					  //}
					  singlemarker<< "\n";
					}
				      else
					{
					  singlemarker << i + 1 << "\t" << map[i].chr << "\t" << map[i].rs << "\t" << map[i].pos << "\t" << map[i].gene << "\t" << codesA[i].a1 << "\t" << codesA[i].a2 << "\t" << map[i].maf <<"/"<< map[i].controlmaf <<"/"<< map[i].casemaf << "\t" << map[i].missingrate << "\t" << counts[thread_nloop][i].det->testHWE_Ca << "\t" << counts[thread_nloop][i].det->testHWE_Co << "\t" << map[i].p << "\t" << map[i].pmod << "\t"<< counts[thread_nloop][i].det->aCaN << "\t" << counts[thread_nloop][i].det->bCaN << "\t" << counts[thread_nloop][i].det->aCoN << "\t" << counts[thread_nloop][i].det->bCoN << "\t" << counts[thread_nloop][i].det->aCa << "\t" << counts[thread_nloop][i].det->bCa << "\t" << counts[thread_nloop][i].det->aCo << "\t" << counts[thread_nloop][i].det->bCo << "\t" << counts[thread_nloop][i].det->orA << "\t" << counts[thread_nloop][i].det->lclA << "\t" << counts[thread_nloop][i].det->rclA << "\t" << counts[thread_nloop][i].det->orB << "\t" << counts[thread_nloop][i].det->lclB << "\t" << counts[thread_nloop][i].det->rclB << "\t" << (*counts[thread_nloop][i].result1).b[1] << "\t" << (*counts[thread_nloop][i].result1).betaNew_se[1]<< "\t" << (*counts[thread_nloop][i].result1).oddsRatio[1] << "\t" << (*counts[thread_nloop][i].result1).lcloddsRatio[1] << "\t" << (*counts[thread_nloop][i].result1).rcloddsRatio[1];
					  if(singleMarkerTest ==7) {singlemarker << "\tH\th";}
					  //if(singleMarkerTest ==3)
					  //{
					  if(covariancematrix)
					    {
					      for(int m=0; m < (df_L1_Single*(df_L1_Single+1))/2; m++)
						{
						  singlemarker << "\t" << (*counts[thread_nloop][i].result1).sigma1[SigmaAuxVector_Single[m]];
						}
					    }
					  //}
					  singlemarker << "\n";
					}

				    }
				  else if (singleMarkerTest == 4 )
				    {
				      if (qt == 1)
					{
					  singlemarker << i + 1 << "\t" << map[i].chr << "\t" << map[i].rs << "\t" << map[i].pos << "\t"
						       << map[i].gene << "\t" << codesA[i].a1 << "\t" << codesA[i].a2 << "\t"
						       << map[i].maf <<"/"<< map[i].controlmaf <<"/"<< map[i].casemaf << "\t" << map[i].missingrate
						       << "\t" << counts[thread_nloop][i].det->testHWE_Ca << "\t" << map[i].p << "\t" << map[i].pmod
						       << "\t" << counts[thread_nloop][i].AA_Ca << "\t" << counts[thread_nloop][i].AB_Ca << "\t"
						       << counts[thread_nloop][i].BB_Ca << "\t" << counts[thread_nloop][i].det->aaCa << "\t"
						       << counts[thread_nloop][i].det->abCa << "\t" << counts[thread_nloop][i].det->bbCa
						       << "\t" << (*counts[thread_nloop][i].result1).b[1] << "\t" << (*counts[thread_nloop][i].result1).betaNew_se[1]
						       << "\t"  << (*counts[thread_nloop][i].result1).betaNew_lcl[1] << "\t" << (*counts[thread_nloop][i].result1).betaNew_rcl[1]
						       << "\t" << (*counts[thread_nloop][i].result1).b[2] << "\t" << (*counts[thread_nloop][i].result1).betaNew_se[2]
						       << "\t" << (*counts[thread_nloop][i].result1).betaNew_lcl[2] << "\t" << (*counts[thread_nloop][i].result1).betaNew_rcl[2];

					  if(covariancematrix)
					    {
					      for(int m=0; m < (df_L1_Single*(df_L1_Single+1))/2; m++)
						{
						  singlemarker << "\t" << (*counts[thread_nloop][i].result1).sigma1[SigmaAuxVector_Single[m]];
						}
					    }
					  singlemarker << "\n";

					}
				      else
					{
					  //singlemarker << i + 1 << "\t" << map[i].chr << "\t" << map[i].rs << "\t" << map[i].pos << "\t" << map[i].gene << "\t" << codesA[i].a1 << "\t" << codesA[i].a2 << "\t" << map[i].maf <<"/"<< map[i].controlmaf <<"/"<< map[i].casemaf << "\t" << map[i].missingrate << "\t" << counts[thread_nloop][i].det->testHWE_Ca << "\t" << counts[thread_nloop][i].det->testHWE_Co << "\t" << map[i].p << "\t" << map[i].pmod << "\t" << counts[thread_nloop][i].AA_Ca << "\t" << counts[thread_nloop][i].AB_Ca << "\t" << counts[thread_nloop][i].BB_Ca << "\t" << counts[thread_nloop][i].AA_Co << "\t" << counts[thread_nloop][i].AB_Co << "\t" << counts[thread_nloop][i].BB_Co << "\t" << counts[thread_nloop][i].det->aaCa << "\t" << counts[thread_nloop][i].det->abCa << "\t" << counts[thread_nloop][i].det->bbCa << "\t" << counts[thread_nloop][i].det->aaCo << "\t" << counts[thread_nloop][i].det->abCo << "\t" << counts[thread_nloop][i].det->bbCo << "\t" << counts[thread_nloop][i].result1.b[1] << "\t" << counts[thread_nloop][i].result1.betaNew_se[1] << "\t" << counts[thread_nloop][i].result1.oddsRatio[1] << "\t" << counts[thread_nloop][i].result1.lcloddsRatio[1] << "\t" << counts[thread_nloop][i].result1.rcloddsRatio[1] << "\t" << counts[thread_nloop][i].result1.b[2] << "\t" << counts[thread_nloop][i].result1.betaNew_se[2] << "\t" << counts[thread_nloop][i].result1.oddsRatio[2] << "\t" << counts[thread_nloop][i].result1.lcloddsRatio[2] << "\t" << counts[thread_nloop][i].result1.rcloddsRatio[2] << "\n";
					  singlemarker << i + 1 << "\t" << map[i].chr << "\t" << map[i].rs << "\t" << map[i].pos << "\t" << map[i].gene << "\t" << codesA[i].a1 << "\t" << codesA[i].a2 << "\t" << map[i].maf <<"/"<< map[i].controlmaf <<"/"<< map[i].casemaf << "\t" << map[i].missingrate << "\t" << counts[thread_nloop][i].det->testHWE_Ca << "\t" << counts[thread_nloop][i].det->testHWE_Co << "\t" << map[i].p << "\t" << map[i].pmod << "\t" << counts[thread_nloop][i].AA_Ca << "\t" << counts[thread_nloop][i].AB_Ca << "\t" << counts[thread_nloop][i].BB_Ca << "\t" << counts[thread_nloop][i].AA_Co << "\t" << counts[thread_nloop][i].AB_Co << "\t" << counts[thread_nloop][i].BB_Co << "\t" << counts[thread_nloop][i].det->aaCa << "\t" << counts[thread_nloop][i].det->abCa << "\t" << counts[thread_nloop][i].det->bbCa << "\t" << counts[thread_nloop][i].det->aaCo << "\t" << counts[thread_nloop][i].det->abCo << "\t" << counts[thread_nloop][i].det->bbCo << "\t" << (*counts[thread_nloop][i].result1).b[1] << "\t" << (*counts[thread_nloop][i].result1).betaNew_se[1] << "\t" << (*counts[thread_nloop][i].result1).oddsRatio[1] << "\t" << (*counts[thread_nloop][i].result1).lcloddsRatio[1] << "\t" << (*counts[thread_nloop][i].result1).rcloddsRatio[1] << "\t" << (*counts[thread_nloop][i].result1).b[2] << "\t" << (*counts[thread_nloop][i].result1).betaNew_se[2] << "\t" << (*counts[thread_nloop][i].result1).oddsRatio[2] << "\t" << (*counts[thread_nloop][i].result1).lcloddsRatio[2] << "\t" << (*counts[thread_nloop][i].result1).rcloddsRatio[2];

					  if(covariancematrix)
					    {
					      for(int m=0; m < (df_L1_Single*(df_L1_Single+1))/2; m++)
						{
						  singlemarker << "\t" << (*counts[thread_nloop][i].result1).sigma1[SigmaAuxVector_Single[m]];
						}

					    }
					  singlemarker << "\n";


					}
				    }
				}
			    }
			  if(q==qstart)
			    {
			      singlemarkerTop.close();
			      singlemarker.close();
			    }

			} // if n == 0 qstart etc




		      // Liability Model
		      if (liability == 1 || liability == 2)
			{
			  if (riskSNP == 0)
			    {
			      for (int a = 0; a < SNPcovar.size(); a++)
				{
				  SNPrisk.push_back("R");
				}
			    }
			  else
			    {
			      if (SNPrisk.size() != SNPcovar.size())
				{
				  cout << "\n\nLiability: Number of risk SNPs differs from the number of risk alleles.\n";
				  logfile << "\n\nLiability: Number of risk SNPs differs from the number of risk alleles.\n";
				  errorfile << "\n\nLiability: Number of risk SNPs differs from the number of risk alleles.\n";
				  exit(1);
				}
			    }

			  cout << "\n\nnumberOfAllCov: " << maxIndexCov << "\n";
			  cout << "numberOfSNPCov: " << numberOfSNPCov << "\n";

			  setLoad(thread_nloop,npplqc, PPLMap, SNPMap, BinSNPs, PplLocations,maxIndexCov, numberOfSNPCov, nsnpqc, person, map, counts, SNPcovar, SNPrisk, codesA);

			  cout << "\n\n";
			  struct PPLLOCATION guy;
			  liabilityResult.open(liabilityfile.c_str(), ios::out);

			  if (liability == 1)
			    {
			      for (int l = 1; l < 2*numberOfSNPCov+1; l++)
				{

				  liabilityCut = l;

				  if (qt != 1)
				    {
				      alt = 1;
				      xType = 0;

				      result1Single[thread_nloop] = logreg(xvec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, 0, 0, 0, inflationfactor, casecounts5, controlcounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[thread_nloop][i], D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, 3,thread_nloop,npplqc, PPLMap, BinSNPs, PplLocations,n,resultSingle[thread_nloop],maxIndexCov,liabilityCut,singleMarkerTest,-1,0,0,0,dosage,genoWeights);

				      alt=0;
				      if(result1Single[thread_nloop].df>0)
					{
					  result2Single[thread_nloop] = logreg(zvec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, 0, 0, 0, inflationfactor, casecounts5, controlcounts5, p, newbeta, X, Xmod, Xt, A,  UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[thread_nloop][i], D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, 3,thread_nloop,npplqc, PPLMap, BinSNPs, PplLocations,n,resultSingle[thread_nloop],maxIndexCov,liabilityCut,singleMarkerTest,-1,0,0,0,dosage,genoWeights);
					}
				      else
					{
					  result2Single[thread_nloop].df=0;result2Single[thread_nloop].sc=0;
					}

				      tLog = 2*(result1Single[thread_nloop].sc -result2Single[thread_nloop].sc);
				      df = result1Single[thread_nloop].df -result2Single[thread_nloop].df;


				      if (df >= 1)
					{
					  newValue = pValueCalc(df/2, tLog/2);
					}
				      else
					{
					  newValue = 1;
					}
				    }
				  else //qt
				    {
				      newValue = 1;

				      result1Single[thread_nloop] = qtreg(xsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, 0, 0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1,alt,Yhelp, xType, female, male, counts[thread_nloop][i], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, 3,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,n,resultSingle[thread_nloop],0,0,maxIndexCov,liabilityCut,singleMarkerTest,-1,covariancematrix,df_L1_Single,df_L2_Single);

				      if(result1Single[thread_nloop].df>0)
					{
					  alt=0;
					  result2Single[thread_nloop] = qtreg(zsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, 0, 0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[thread_nloop][i], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, 3,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs,  PplLocations,n,resultSingle[thread_nloop],0,0,maxIndexCov,liabilityCut,singleMarkerTest,-1,covariancematrix,df_L1_Single,df_L2_Single);
					}
				      else
					{
					  result2Single[thread_nloop].df=0;result2Single[thread_nloop].sc=0;
					}

				      df = result2Single[thread_nloop].df -result1Single[thread_nloop].df; //!!
				      Fstat = ((result2Single[thread_nloop].sc -result1Single[thread_nloop].sc)/(df))/(result1Single[thread_nloop].sc/result1Single[thread_nloop].df);

				      if (df >= 1 && result2Single[thread_nloop].sc > 0.000001)
					{
					  newValue = betai(result1Single[thread_nloop].df/2,df/2,result1Single[thread_nloop].df/(result1Single[thread_nloop].df+df*Fstat));
					}
				      else
					{
					  newValue = 1;
					}
				    }



				  // Results

				  // Correction Sidak
				  if (newValue != 1)
				    {
				      corrL++;
				    }
				  if (l == 2*numberOfSNPCov && corrL==0)
				    {
				      corrL=1;
				    }


				  if (newValue < bestL)
				    {
				      bestL = newValue;
				    }


				  if (n == 0)
				    {
				      if (!qt)
					{
					  if (l == 1)
					    {
					      liabilityResult << "Allele_Load\tp_Value\tbeta1\tse1\tOR_unadjusted\tseOR\tFreq_Cases_Loads\tFreq_Control_Loads\n";
					      liabilityResult << l << "\t" << newValue << "\t" << result1Single[thread_nloop].b[1] << "\t" <<result1Single[thread_nloop].betaNew_se[1] << "\t" << result1Single[thread_nloop].oddsRatio[1] << "\t" << result1Single[thread_nloop].lcloddsRatio[1] << "\t" << result1Single[thread_nloop].betaNew_lcl[1] << "\t" << result1Single[thread_nloop].betaNew_rcl[1] << "\n";
					      cout << "Allele_Load\tp_Value\tbeta1\tse1\tOR_unadjusted\tseOR\tFreq_Cases_Loads\tFreq_Control_Loads\n";
					      cout << l << "\t" << newValue << "\t" << result1Single[thread_nloop].b[1] << "\t" <<result1Single[thread_nloop].betaNew_se[1] << "\t" << result1Single[thread_nloop].oddsRatio[1] << "\t" << result1Single[thread_nloop].lcloddsRatio[1] << "\t" << result1Single[thread_nloop].betaNew_lcl[1] << "\t" << result1Single[thread_nloop].betaNew_rcl[1] << "\n";
					    }
					  else if (l == 2*numberOfSNPCov)
					    {
					      newValueCorr = Sidak(bestL,corrL);
					      liabilityResult << l << "\t" << newValue << "\t" << result1Single[thread_nloop].b[1] << "\t" << result1Single[thread_nloop].betaNew_se[1] << "\t" << result1Single[thread_nloop].oddsRatio[1] << "\t" << result1Single[thread_nloop].lcloddsRatio[1] << "\t" << result1Single[thread_nloop].betaNew_lcl[1] << "\t" << result1Single[thread_nloop].betaNew_rcl[1] << "\n";
					      liabilityResult << "Total p-Value: " << newValueCorr << "\n";
					      cout << l << "\t" << newValue << "\t" << result1Single[thread_nloop].b[1] << "\t" << result1Single[thread_nloop].betaNew_se[1] << "\t" << result1Single[thread_nloop].oddsRatio[1] << "\t" << result1Single[thread_nloop].lcloddsRatio[1] << "\t" << result1Single[thread_nloop].betaNew_lcl[1] << "\t" << result1Single[thread_nloop].betaNew_rcl[1] << "\n";
					      cout << "Total p-Value: " << newValueCorr << "\n";
					    }
					  else
					    {
					      liabilityResult << l << "\t" << newValue << "\t" << result1Single[thread_nloop].b[1] << "\t" << result1Single[thread_nloop].betaNew_se[1]  << "\t" << result1Single[thread_nloop].oddsRatio[1] << "\t" << result1Single[thread_nloop].lcloddsRatio[1] << "\t" << result1Single[thread_nloop].betaNew_lcl[1] << "\t" << result1Single[thread_nloop].betaNew_rcl[1] << "\n";
					      cout << l << "\t" << newValue << "\t" << result1Single[thread_nloop].b[1] << "\t" << result1Single[thread_nloop].betaNew_se[1]  << "\t" << result1Single[thread_nloop].oddsRatio[1] << "\t" << result1Single[thread_nloop].lcloddsRatio[1] << "\t" << result1Single[thread_nloop].betaNew_lcl[1] << "\t" << result1Single[thread_nloop].betaNew_rcl[1] << "\n";
					    }
					}
				      else
					{
					  if (l == 1)
					    {
					      liabilityResult << "Allele_Load\tp_Value\tbeta1\tse1\n";
					      liabilityResult << l << "\t" << newValue << "\t" << result1Single[thread_nloop].b[1] << "\t" <<result1Single[thread_nloop].betaNew_se[1] << "\n";
					      cout << "Allele_Load\tp_Value\tbeta1\tse1\n";
					      cout << l << "\t" << newValue << "\t" << result1Single[thread_nloop].b[1] << "\t" <<result1Single[thread_nloop].betaNew_se[1] << "\n";
					    }
					  else if (l == 2*numberOfSNPCov)
					    {
					      newValueCorr = Sidak(bestL,corrL);
					      liabilityResult << l << "\t" << newValue << "\t" << result1Single[thread_nloop].b[1] << "\t" << result1Single[thread_nloop].betaNew_se[1] << "\n";
					      liabilityResult << "Total p-Value: " << newValueCorr << "\n";
					      cout << l << "\t" << newValue << "\t" << result1Single[thread_nloop].b[1] << "\t" << result1Single[thread_nloop].betaNew_se[1] << "\n";
					      cout << "Total p-Value: " << newValueCorr << "\n";
					    }
					  else
					    {
					      liabilityResult << l << "\t" << newValue << "\t" << result1Single[thread_nloop].b[1] << "\t" << result1Single[thread_nloop].betaNew_se[1]  << "\n";
					      cout << l << "\t" << newValue << "\t" << result1Single[thread_nloop].b[1] << "\t" << result1Single[thread_nloop].betaNew_se[1]  << "\n";
					    }

					}
				    }
				}
			    }
			  if(liability == 2)
			    {
			      liabilityCut = -1;

			      if (qt != 1)
				{
				  alt = 1;
				  xType = 0;
				  result1Single[thread_nloop] = logreg(xvec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, 0, 0, 0, inflationfactor, casecounts5, controlcounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[thread_nloop][i], D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, 3,thread_nloop,npplqc, PPLMap, BinSNPs, PplLocations,n,resultSingle[thread_nloop],maxIndexCov,liabilityCut,singleMarkerTest,-1,0,0,0,dosage,genoWeights);
				  alt=0;
				  if(result1Single[thread_nloop].df>0)
				    {
				      result2Single[thread_nloop] = logreg(zvec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, 0, 0, 0, inflationfactor, casecounts5, controlcounts5, p, newbeta, X, Xmod, Xt, A,  UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[thread_nloop][i], D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, 3,thread_nloop,npplqc, PPLMap, BinSNPs, PplLocations,n,resultSingle[thread_nloop],maxIndexCov,liabilityCut,singleMarkerTest,-1,0,0,0,dosage,genoWeights);
				    }
				  else
				    {
				      result2Single[thread_nloop].df=0;result2Single[thread_nloop].sc=0;
				    }

				  tLog = 2*(result1Single[thread_nloop].sc -result2Single[thread_nloop].sc);
				  df = result1Single[thread_nloop].df -result2Single[thread_nloop].df;

				  if (df >= 1)
				    {
				      newValue = pValueCalc(df/2, tLog/2);
				    }
				  else
				    {
				      newValue = 1;
				    }
				  cout << "p-wert " << newValue << "\n";
				}
			      else //qt
				{
				  newValue = 1;

				  result1Single[thread_nloop] = qtreg(xsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, 0, 0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1,alt,Yhelp, xType, female, male, counts[thread_nloop][i], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, 3,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,n,resultSingle[thread_nloop],0,0,maxIndexCov,liabilityCut,singleMarkerTest,-1,covariancematrix,df_L1_Single,df_L2_Single);

				  if(result1Single[thread_nloop].df>0)
				    {
				      alt=0;
				      result2Single[thread_nloop] = qtreg(zsinglevec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, 0, 0, 0, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N1, alt, Yhelp, xType, female, male, counts[thread_nloop][i], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, 3,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs,  PplLocations,n,resultSingle[thread_nloop],0,0,maxIndexCov,liabilityCut,singleMarkerTest,-1,covariancematrix,df_L1_Single,df_L2_Single);
				    }
				  else
				    {
				      result2Single[thread_nloop].df=0;result2Single[thread_nloop].sc=0;
				    }

				  df = result2Single[thread_nloop].df -result1Single[thread_nloop].df; //!!
				  Fstat = ((result2Single[thread_nloop].sc -result1Single[thread_nloop].sc)/(df))/(result1Single[thread_nloop].sc/result1Single[thread_nloop].df);

				  if (df >= 1 && result2Single[thread_nloop].sc > 0.000001)
				    {
				      newValue = betai(result1Single[thread_nloop].df/2,df/2,result1Single[thread_nloop].df/(result1Single[thread_nloop].df+df*Fstat));
				    }
				  else
				    {
				      newValue = 1;
				    }
				}

			      if (n == 0)
				{
				  if (!qt)
				    {

				      liabilityResult << "p_Value\tbeta1\tse1\n";
				      liabilityResult  << newValue << "\t" << result1Single[thread_nloop].b[1] << "\t" <<result1Single[thread_nloop].betaNew_se[1] << "\n";
				      cout << "p_Value\tbeta1\tse1\n";
				      cout <<  newValue << "\t" << result1Single[thread_nloop].b[1] << "\t" <<result1Single[thread_nloop].betaNew_se[1] << "\n";


				    }
				  else
				    {

				      liabilityResult << "p_Value\tbeta1\tse1\n";
				      liabilityResult << newValue << "\t" << result1Single[thread_nloop].b[1] << "\t" <<result1Single[thread_nloop].betaNew_se[1] << "\n";
				      cout << "p_Value\tbeta1\tse1\n";
				      cout << newValue << "\t" << result1Single[thread_nloop].b[1] << "\t" <<result1Single[thread_nloop].betaNew_se[1] << "\n";


				    }
				}
			    }
			  liabilityResult.close();
			}


		      // Einzelne SNPs auswählen
		      if (n == 0 || combilist == 1)
			{

			  snp1Pos = -1;
			  snp2Pos = -1;
			  snp3Pos = -1;

			  if(combilist!=1 && nsnps>0)
			    {
			      for (i = 0; i < nlinestped; i++)
				{
				  if (map[i].analysis_in == 1)
				    {

				      if (nsnps>=1 && strcmp(snp1,map[i].rs)==0)
					{
					  snp1Pos = i;
					}

				      if (nsnps>=2 && strcmp(snp2,map[i].rs)==0)
					{

					  snp2Pos = i;
					}

				      if (nsnps>=3 && strcmp(snp3,map[i].rs)==0)
					{
					  snp3Pos = i;
					}
				    }
				}
			    }


			  if (combilist == 1)
			    {
			      if(pos1 != -1 && map[pos1].analysis_in==1) {snp1Pos = pos1;}
			      if(pos2 != -1 && map[pos2].analysis_in==1) {snp2Pos = pos2;}
			      if(markercombi3 && pos3 != -1 && map[pos3].analysis_in==1) {snp3Pos = pos3;}

			      if (snp1Pos == -1)
				{
				  //errorfile << "SNP1 is not left after QC\n";
				  //logfile << "SNP1 is not left after QC\n";
				  //cout << "SNP1 is not left after QC\n";
				}
			      if (snp2Pos == -1)
				{
				  //errorfile << "SNP2 is not left after QC\n";
				  //logfile << "SNP2 is not left after QC\n";
				  //cout << "SNP2 is not left after QC\n";
				}
			      if (markercombi3 == 1 && snp3Pos == -1)
				{
				  //errorfile << "SNP3 is not left after QC\n";
				  //logfile << "SNP3 is not left after QC\n";
				  //cout << "SNP3 is not left after QC\n";
				}
			    }
			  else // not combilist
			    {
			      // Ausgewählter SNP ist nicht in der Liste verhanden
			      if ((snp1_Selected && snp1Pos == -1) || (snp2_Selected && snp2Pos == -1) || (snp3_Selected && snp3Pos == -1))
				{
				  errorfile << "Selected SNP is not found in the tpedfiles\n";
				  logfile << "Selected SNP is not found in the tpedfiles\n";
				  cout << "Selected SNP is not found in the tpedfiles\n";
				  errorfile.close();logfile.close();exit(1);
				}


			      if (nsnps>0)
				{
				  cout << "\nselected snp1: " << map[snp1Pos].rs << "\n";
				  logfile << "\nselected snp1: " << map[snp1Pos].rs << "\n";

				  cout << "snp1Pos: " << snp1Pos << "\n";
				  logfile << "snp1Pos: " << snp1Pos << "\n";

				}
			      if (nsnps>1)// Änderung569
				{
				  cout << "selected snp2: " << map[snp2Pos].rs << "\n";
				  logfile << "selected snp2: " << map[snp2Pos].rs << "\n";

				  cout << "snp2Pos: " << snp2Pos << "\n";
				  logfile << "snp2Pos: " << snp2Pos << "\n";
				}
			      if (nsnps>2) // Änderung569
				{
				  cout << "selected snp3: " << map[snp3Pos].rs << "\n";
				  logfile << "selected snp3: " << map[snp3Pos].rs << "\n";

				  cout << "snp3Pos: " << snp3Pos << "\n";
				  logfile << "snp3Pos: " << snp3Pos << "\n";
				}
			    }
			} // end if n == 0 || combilist == 1
		    } // nur einmal r < 1

		  if(markercombi2 ==1)
		    {
		      if (n == 0 && r <= 1 && combilist!=1) // keine Pathwayinformation oder erster Pathway TIM
			{
			  // Anzahl der ausgewählten Tests ermitteln: bei Pathway ist es Obergrenze"
			  ntestsUser = 0;
			  ntestsPlus = 0;
			  if (selectedSnp == 0)
			    {
			      if (mWithSingletop == 0)
				{
				  if (mWithGenetictop == 0)
				    {
				      ntestsUser = (((double)nsnpqc * ((double)nsnpqc - 1)) / 2);
				    }
				  else if (mWithGenetictop == 1)
				    {
				      ntestsUser = ((double)nsnpqc - (double)genetictop) * (double)genetictop + (double)genetictop *((double)genetictop-1)/2;
				    }
				  else if (mWithGenetictop == 2)
				    {
				      ntestsUser = (double)genetictop * ((double)genetictop-1)/2;
				    }
				  ntestsPlus = ntestsUser;
				}
			      else if (mWithSingletop == 1)
				{
				  if (mWithGenetictop == 0)
				    {
				      ntestsUser = ((double)nsnpqc - (double)singletop) * (double)singletop + (double)singletop *((double)singletop-1)/2;
				      ntestsPlus = (((double)nsnpqc * ((double)nsnpqc - 1)) / 2);
				    }
				  else if (mWithGenetictop == 1)
				    {

				      ntestsUser = (double)singletop*(double)genetictop -((double)overlaptop*((double)overlaptop-1)/2);
				      ntestsPlus = ((double)nsnpqc - (double)genetictop) * (double)genetictop + (double)genetictop *((double)genetictop-1)/2;
				    }
				  else if (mWithGenetictop == 2)
				    {
				      ntestsUser = (double)overlaptop*(double)genetictop -((double)overlaptop*((double)overlaptop-1)/2);
				      ntestsPlus = (double)genetictop * ((double)genetictop-1)/2;
				    }
				}
			      else if (mWithSingletop == 2)
				{
				  if (mWithGenetictop == 0)
				    {
				      ntestsUser = (double)singletop * ((double)singletop-1)/2;
				      ntestsPlus = (((double)nsnpqc * ((double)nsnpqc - 1)) / 2);
				    }
				  else if (mWithGenetictop == 1)
				    {
				      ntestsUser = (double)overlaptop*(double)singletop -((double)overlaptop*((double)overlaptop-1)/2);
				      ntestsPlus = ((double)nsnpqc - (double)genetictop) * (double)genetictop + (double)genetictop *((double)genetictop-1)/2;
				    }
				  else if (mWithGenetictop == 2)
				    {
				      ntestsUser = (double)overlaptop *((double)overlaptop-1)/2;
				      ntestsPlus = (double)genetictop * ((double)genetictop-1)/2;
				    }
				}
			    }

			  else if (selectedSnp == 1)
			    {
			      if (mWithSingletop == 0)
				{
				  if (mWithGenetictop == 0)

				    {
				      ntestsUser = (double)nsnpqc-1;
				    }
				  else if (mWithGenetictop == 1)

				    {
				      ntestsUser = (double)genetictop-1;
				    }
				  ntestsPlus = ntestsUser;
				}
			      else if (mWithSingletop == 1)
				{
				  if (mWithGenetictop == 0)
				    {
				      ntestsUser = (double)singletop-1;
				      ntestsPlus = (double)nsnpqc-1;
				    }
				  else if (mWithGenetictop == 1)
				    {
				      ntestsUser = (double)overlaptop-1;
				      ntestsPlus = (double)genetictop-1;
				    }
				}
			    }
			  else if (selectedSnp == 2)
			    {
			      ntestsUser = 1;
			      ntestsPlus = 1;
			    }

			}
		    } // end if markercombi2
		  else if (markercombi3 ==1 && n==0 && r<=1 && combilist!=1) //TIM
		    {
		      ntestsUser = 0;
		      ntestsPlus = 0;
		      marker1 = 0;
		      marker2 = 0;
		      marker3 = 0;

		      if (selectedSnp == 0)
			{
			  if (mWithSingletop == 0)
			    {
			      if (mWithGenetictop == 0)
				{
				  ntestsUser = ((((double)nsnpqc-2)*((double)nsnpqc-1)*(double)nsnpqc)/6);
				}
			      else if (mWithGenetictop == 1)
				{
				  ntestsUser = (double)genetictop*((double)nsnpqc-(double)genetictop)*((double)nsnpqc- (double)genetictop - 1)/2
				    + (double)genetictop *((double)genetictop-1)*((double)genetictop-2)/6
				    +(double)genetictop*((double)genetictop-1)/2*((double)nsnpqc-(double)genetictop);
				}
			      else if (mWithGenetictop == 2)
				{
				  ntestsUser = (((double)nsnpqc - (double)genetictop) * ((double)genetictop - 1)* (double)genetictop)/2
				    + (double)genetictop * ((double)genetictop - 1)* ((double)genetictop - 2) / 6;
				}
			      else if (mWithGenetictop == 3)
				{
				  ntestsUser = ((((double)genetictop - 2)*((double)genetictop - 1)* (double)genetictop)/6);
				}
			      ntestsPlus = ntestsUser;
			    }
			  else if (mWithSingletop == 1)
			    {
			      if (mWithGenetictop == 0)
				{
				  ntestsUser = (double)singletop*((double)nsnpqc-(double)singletop)*((double)nsnpqc- (double)singletop - 1)/2
				    + (double)singletop*((double)singletop-1)*((double)singletop-2)/6
				    +(double)singletop*((double)singletop-1)/2*((double)nsnpqc-(double)singletop);
				  ntestsPlus = ((((double)nsnpqc-2)*((double)nsnpqc-1)*(double)nsnpqc)/6);
				}
			      else if (mWithGenetictop == 1)
				{
				  ntestsUser = overlaptop*((double)singletop-(double)overlaptop)*((double)singletop-(double)overlaptop-1)/2
				    + (double)overlaptop*((double)overlaptop-1)/2*((double)singletop-(double)overlaptop)
				    + (double)overlaptop*((double)overlaptop-1)*((double)overlaptop-2)/6
				    + (double)overlaptop*((double)overlaptop-1)/2*((double)genetictop-(double)overlaptop)
				    + ((double)singletop-(double)overlaptop)*((double)singletop-(double)overlaptop-1)/2
				    *((double)genetictop-(double)overlaptop)
				    + (double)overlaptop*((double)overlaptop-1)/2*((double)nsnpqc-(double)singletop-(double)genetictop
										   +(double)overlaptop)
				    + (double)overlaptop*((double)singletop-(double)overlaptop)*((double)genetictop-(double)overlaptop)
				    + (double)overlaptop*((double)singletop-(double)overlaptop)* ((double)nsnpqc-(double)singletop
												  -(double)genetictop+(double)overlaptop)
				    + (double)overlaptop*((double)genetictop-(double)overlaptop)*((double)genetictop
												  *(double)overlaptop-1)/2
				    + (double)overlaptop*((double)genetictop-(double)overlaptop)
				    + ((double)singletop- (double)overlaptop)* ((double)genetictop-(double)overlaptop)
				    *((double)genetictop-(double)overlaptop-1)/2
				    + ((double)singletop-(double)overlaptop)*((double)genetictop-(double)overlaptop)*((double)nsnpqc
														      -(double)singletop-(double)genetictop+(double)overlaptop);

				  ntestsPlus = (double) genetictop*((double) nsnpqc-(double) genetictop)*((double) nsnpqc- (double) genetictop
													  - 1)/2
				    + (double) genetictop *((double) genetictop-1)*((double) genetictop-2)/6+(double) genetictop*((double) genetictop-1)/2*((double) nsnpqc-(double) genetictop);

				}
			      else if (mWithGenetictop == 2)
				{
				  ntestsUser = (double) overlaptop*((double) overlaptop-1)*((double) overlaptop-2)/6
				    + (double) overlaptop*((double) overlaptop-1)/2*((double) genetictop-(double) overlaptop)
				    + (double) overlaptop*((double) overlaptop-1)/2*(singletop-(double) overlaptop)
				    + (double) overlaptop*((double) genetictop-(double) overlaptop)*((double) genetictop-(double) overlaptop-1)/2
				    + (double) overlaptop*((double) genetictop-(double) overlaptop)*(singletop-(double) overlaptop)
				    + ((double) genetictop-(double) overlaptop)*((double) genetictop-(double) overlaptop-1)/2 * (singletop-(double) overlaptop);

				  ntestsPlus = (((double) nsnpqc - (double) genetictop) * ((double) genetictop - 1)* (double) genetictop)/2
				    + (double) genetictop * ((double) genetictop - 1)* ((double) genetictop - 2) / 6;
				}
			      else if (mWithGenetictop == 3)
				{
				  ntestsUser = (double) overlaptop*((double) genetictop-(double) overlaptop)*((double) genetictop-(double) overlaptop-1)/2
				    + (double) overlaptop*((double) overlaptop-1)/2*((double) genetictop-(double) overlaptop)
				    + (double) overlaptop*((double) overlaptop-1)*((double) overlaptop-2)/6;

				  ntestsPlus = ((((double) genetictop - 2)*((double) genetictop - 1)* (double) genetictop)/6);
				}
			    }
			  else if (mWithSingletop == 2)
			    {
			      if (mWithGenetictop == 0)
				{
				  ntestsUser = (((double) nsnpqc - singletop) * (singletop - 1)* singletop) / 2
				    + singletop * (singletop - 1)* (singletop - 2) / 6;

				  ntestsPlus = ((((double) nsnpqc-2)*((double) nsnpqc-1)*(double) nsnpqc)/6);
				}
			      else if (mWithGenetictop == 1)
				{
				  ntestsUser = (double) overlaptop*((double) overlaptop-1)*((double) overlaptop-2)/6
				    + (double) overlaptop*((double) overlaptop-1)/2*(singletop-(double) overlaptop)
				    + (double) overlaptop*((double) overlaptop-1)/2*((double) genetictop-(double) overlaptop)
				    + (double) overlaptop*(singletop-(double) overlaptop)*(singletop-(double) overlaptop-1)/2
				    + (double) overlaptop*(singletop-(double) overlaptop)*((double) genetictop-(double) overlaptop)
				    + (singletop-(double) overlaptop)*(singletop-(double) overlaptop-1)/2 * ((double) genetictop-(double) overlaptop);

				  ntestsPlus = (double) genetictop*((double) nsnpqc-(double) genetictop)*((double) nsnpqc- (double) genetictop - 1)/2
				    + (double) genetictop *((double) genetictop-1)*((double) genetictop-2)/6+(double) genetictop*((double) genetictop-1)/2*((double) nsnpqc-(double) genetictop);
				}
			      else if (mWithGenetictop == 2)
				{
				  ntestsUser = (double) overlaptop*((double) overlaptop-1)*((double) overlaptop-2)/6
				    + (double) overlaptop*((double) overlaptop-1)/2*(singletop-(double) overlaptop)
				    + (double) overlaptop*((double) overlaptop-1)/2*((double) genetictop-(double) overlaptop)
				    + (double) overlaptop*(singletop-(double) overlaptop)*((double) genetictop-(double) overlaptop);

				  ntestsPlus = (((double) nsnpqc - (double) genetictop) * ((double) genetictop - 1)* (double) genetictop)/2
				    + (double) genetictop * ((double) genetictop - 1)* ((double) genetictop - 2) / 6;
				}
			      else if (mWithGenetictop == 3)
				{
				  ntestsUser = (double) overlaptop*((double) overlaptop-1)*((double) overlaptop-2)/6
				    + (double) overlaptop*((double) overlaptop-1)/2*((double) genetictop-(double) overlaptop);

				  ntestsPlus = ((((double) genetictop - 2)*((double) genetictop - 1)* (double) genetictop)/6);
				}
			    }
			  else if (mWithSingletop == 3)
			    {
			      if (mWithGenetictop == 0)
				{
				  ntestsUser = (((singletop - 2)*(singletop - 1)* singletop)/6);

				  ntestsPlus = ((((double) nsnpqc-2)*((double) nsnpqc-1)*(double) nsnpqc)/6);
				}
			      else if (mWithGenetictop == 1)
				{
				  ntestsUser = (double) overlaptop*(singletop-(double) overlaptop)*(singletop-(double) overlaptop-1)/2
				    + (double) overlaptop*((double) overlaptop-1)/2*(singletop-(double) overlaptop)
				    + (double) overlaptop*((double) overlaptop-1)*((double) overlaptop-2)/6;

				  ntestsPlus = (double) genetictop*((double) nsnpqc-(double) genetictop)*((double) nsnpqc- (double) genetictop - 1)/2

				    + (double) genetictop *((double) genetictop-1)*((double) genetictop-2)/6+(double) genetictop*((double) genetictop-1)/2*((double) nsnpqc-(double) genetictop);
				}
			      else if (mWithGenetictop == 2)
				{
				  ntestsUser = (double) overlaptop*((double) overlaptop-1)/2*(singletop-(double) overlaptop)
				    + (double) overlaptop*((double) overlaptop-1)*((double) overlaptop-2)/6;

				  ntestsPlus = (((double) nsnpqc - (double) genetictop) * ((double) genetictop - 1)* (double) genetictop)/2
				    + (double) genetictop * ((double) genetictop - 1)* ((double) genetictop - 2) / 6;
				}
			      else if (mWithGenetictop == 3)
				{
				  ntestsUser = ((double) overlaptop-2)*((double) overlaptop-1)*(double) overlaptop/6;
				  ntestsPlus = ((((double) genetictop - 2)*((double) genetictop - 1)* (double) genetictop)/6);
				}
			    }
			}
		      else if (selectedSnp == 1)
			{
			  if (mWithSingletop == 0)
			    {
			      if (mWithGenetictop == 0)
				{
				  ntestsUser = ((((double) nsnpqc-1)*(double) nsnpqc)/2);
				}
			      else if (mWithGenetictop == 1)
				{

				  ntestsUser = ((double) nsnpqc - (double) genetictop) * (double) genetictop + (double) genetictop* ((double) genetictop - 1)/2;
				}
			      else if (mWithGenetictop == 2)
				{
				  ntestsUser = ((((double) genetictop - 1) * (double) genetictop) / 2);
				}
			      ntestsPlus = ntestsUser;
			    }
			  else if (mWithSingletop == 1)
			    {
			      if (mWithGenetictop == 0)
				{
				  ntestsUser = ((double) nsnpqc - singletop) * singletop + singletop*(singletop - 1)/2;
				  ntestsPlus = ((((double) nsnpqc-1)*(double) nsnpqc)/2);
				}
			      else if (mWithGenetictop == 1)
				{
				  ntestsUser = (double) genetictop*singletop -((double) overlaptop*((double) overlaptop-1)/2);
				  ntestsPlus = ((double) nsnpqc - (double) genetictop) * (double) genetictop + (double) genetictop* ((double) genetictop - 1)/2;
				}
			      else if (mWithGenetictop == 2)
				{
				  ntestsUser = (double) overlaptop*(double) genetictop -((double) overlaptop*((double) overlaptop-1)/2);
				  ntestsPlus = ((((double) genetictop - 1) * (double) genetictop) / 2);
				}
			    }
			  else if (mWithSingletop == 2)
			    {
			      if (mWithGenetictop == 0)
				{
				  ntestsUser = (((singletop - 1) * singletop) / 2);
				  ntestsPlus = ((((double) nsnpqc-1)*(double) nsnpqc)/2);
				}
			      else if (mWithGenetictop == 1)
				{
				  ntestsUser = (double) overlaptop*singletop -((double) overlaptop*((double) overlaptop-1)/2);
				  ntestsPlus = ((double) nsnpqc - (double) genetictop) * (double) genetictop + (double) genetictop* ((double) genetictop - 1)/2;
				}
			      else if (mWithGenetictop == 2)
				{
				  ntestsUser = (double) overlaptop *((double) overlaptop-1)/2;
				  ntestsPlus = ((((double) genetictop - 1) * (double) genetictop) / 2);
				}
			    }
			}
		      else if (selectedSnp == 2)
			{
			  if (mWithSingletop == 0)
			    {
			      if (mWithGenetictop == 0)
				{
				  ntestsUser = (double) nsnpqc-2;
				}
			      else if (mWithGenetictop == 1)
				{
				  ntestsUser = (double) genetictop;
				}
			      ntestsPlus = ntestsUser;
			    }

			  else if (mWithSingletop == 1)
			    {
			      if (mWithGenetictop == 0)
				{
				  ntestsUser = singletop;
				  ntestsPlus = (double) nsnpqc-2;
				}
			      else if (mWithGenetictop == 1)
				{
				  ntestsUser = (double) overlaptop;
				  ntestsPlus = (double) genetictop;
				}
			    }
			}
		      else if (selectedSnp == 3)
			{
			  ntestsUser = 1;
			  ntestsPlus = 1;
			}
		    }//markercomb3 anzahl tests

		  if (combilist == 1 && n==0)
		    {
		      ntestsUser = nlinescombi;
		      ntestsPlus = nlinescombi;
		    }


		  if (r > 0)
		    {
		      // Overlapneeded
		      if (markercombi2 == 1)
			{
			  if (mWithSingletop == 2 && mWithGenetictop == 2)
			    {
			      overlapneeded = 2;
			    }
			  else if ((mWithSingletop == 1 && mWithGenetictop == 2) || (mWithSingletop == 2 && mWithGenetictop == 1))
			    {
			      overlapneeded = 1;
			    }
			  else
			    {
			      overlapneeded = 0;
			    }
			}
		      else if (markercombi3 == 1)
			{
			  if (mWithSingletop == 3 && mWithGenetictop == 3)
			    {
			      overlapneeded = 3;
			    }
			  else if ((mWithSingletop == 3 && mWithGenetictop == 2) || (mWithSingletop == 2 && mWithGenetictop == 3))
			    {
			      overlapneeded = 2;
			    }
			  else if ((mWithSingletop == 3 && mWithGenetictop == 1) || (mWithSingletop == 1 && mWithGenetictop == 3)
				   || (mWithSingletop == 2 && mWithGenetictop == 2))
			    {
			      overlapneeded = 1;
			    }
			  else
			    {
			      overlapneeded = 0;
			    }
			}

		      pathway[r-1].singletop[thread_nloop] = 0;
		      pathway[r-1].singletop0[thread_nloop] = 0;
		      pathway[r-1].overlaptop[thread_nloop] = 0;
		      pathway[r-1].singletopmod[thread_nloop] = 0;
		      pathway[r-1].intertop[thread_nloop] = 0;
		      pathway[r-1].intertopmod[thread_nloop] = 0;

		      if (pathway[r-1].genetictop < mWithGenetictop)
			{
			  continue;
			}

		      // Singletopliste in Kombination mit dem Pathway
		      if (mWithSingletop != 0 || pathwayAnalysis == 1) //oder PWT
			{
			  for (k = 0; k < singletop; k++)
			    {
			      for(l = 0; l < pathway[r-1].counts; l++)
				{
				  if (bestsinglemarker[k].nr == pathway[r-1].list[l])
				    {
				      pathway[r-1].singletop[thread_nloop]++;
				      pathway[r-1].singletopmod[thread_nloop]++;
				      pathway[r-1].listsingle[thread_nloop] = (int*) realloc(pathway[r-1].listsingle[thread_nloop], pathway[r-1].singletop[thread_nloop]*sizeof(int));
				      pathway[r-1].listsingle_gene[thread_nloop] = (int*) realloc(pathway[r-1].listsingle_gene[thread_nloop], pathway[r-1].singletop[thread_nloop]*sizeof(int));

				      if (!pathway[r-1].listsingle || ! pathway[r-1].list_gene)
					{
					  errorfile << "memory allocation error in pathway[r-1].listsingle\n";
					  logfile << "memory allocation error in pathway[r-1].listsingle\n";
					  cout << "memory allocation error in pathway[r-1].listsingle\n";
					  errorfile.close();logfile.close();exit(1);
					}



				      pathway[r-1].listsingle[thread_nloop][pathway[r-1].singletop[thread_nloop]-1] = bestsinglemarker[k].nr;
				      pathway[r-1].listsingle_gene[thread_nloop][pathway[r-1].singletop[thread_nloop]-1] = pathway[r-1].list_gene[l];

				      if (pathwayAnalysis == 1)
					{
					  pathway[r-1].listp[thread_nloop] = (double*) realloc(pathway[r-1].listp[thread_nloop], pathway[r-1].singletop[thread_nloop]*sizeof(double));


					  if (!pathway[r-1].listp[thread_nloop])
					    {
					      errorfile << "memory allocation error in pathway[r-1].listp[thread_nloop]\n";
					      logfile << "memory allocation error in pathway[r-1].listp[thread_nloop]\n";
					      cout << "memory allocation error in pathway[r-1].listp[thread_nloop]\n";
					      errorfile.close();logfile.close();exit(1);
					    }

					  //PWT: listp anlegen, in listp bestsingelemarker[k].p schreiben
					  pathway[r-1].listp[thread_nloop][pathway[r-1].singletop[thread_nloop]-1] = bestsinglemarker[k].p;

					  //pplus: falls gen schon drin: wieder raus
					  if (pathwayTest >=3)
					    {
					      pathway[r-1].listp_in[thread_nloop] = (int*) realloc(pathway[r-1].listp_in[thread_nloop], pathway[r-1].singletop[thread_nloop]*sizeof(int));
					      pathway[r-1].listp_in[thread_nloop][pathway[r-1].singletop[thread_nloop]-1] = 1;
					      if (!pathway[r-1].listp_in[thread_nloop])
						{
						  errorfile << "memory allocation error in pathway[r-1].listp[thread_nloop]\n";
						  logfile << "memory allocation error in pathway[r-1].listp[thread_nloop]\n";
						  cout << "memory allocation error in pathway[r-1].listp[thread_nloop]\n";
						  errorfile.close();logfile.close();exit(1);
						}

					      if (pathwayTest ==5) // NEU Pathwaytest
						{
						  pathway[r-1].listp_in2[thread_nloop] = (int*) realloc(pathway[r-1].listp_in2[thread_nloop], pathway[r-1].singletop[thread_nloop]*sizeof(int));
						  pathway[r-1].listp_in2[thread_nloop][pathway[r-1].singletop[thread_nloop]-1] = 0;

						  pathway[r-1].listp_in3[thread_nloop] = (int*) realloc(pathway[r-1].listp_in3[thread_nloop], pathway[r-1].singletop[thread_nloop]*sizeof(int));
						  pathway[r-1].listp_in3[thread_nloop][pathway[r-1].singletop[thread_nloop]-1] = 0;

						  if (!pathway[r-1].listp_in2[thread_nloop] || !pathway[r-1].listp_in3[thread_nloop])
						    {
						      errorfile << "memory allocation error in pathway[r-1].listp[thread_nloop]\n";
						      logfile << "memory allocation error in pathway[r-1].listp[thread_nloop]\n";
						      cout << "memory allocation error in pathway[r-1].listp[thread_nloop]\n";
						      errorfile.close();logfile.close();exit(1);
						    }
						}

					      for(kkk=0;kkk<pathway[r-1].singletop[thread_nloop]-1;kkk++)
						{
						  if(pathway[r-1].listp_in[thread_nloop][kkk] == 1 && pathway[r-1].listsingle_gene[thread_nloop][kkk]==pathway[r-1].list_gene[l])
						    {

						      pathway[r-1].singletopmod[thread_nloop]--;

						      if(pathway[r-1].listp[thread_nloop][pathway[r-1].singletop[thread_nloop]-1] < pathway[r-1].listp[thread_nloop][kkk])
							{
							  pathway[r-1].listp_in[thread_nloop][kkk] = 0;
							}
						      else
							{
							  pathway[r-1].listp_in[thread_nloop][pathway[r-1].singletop[thread_nloop]-1] = 0;
							}
						      break;
						    }
						}

					    } //pathwayTest >=3
					} // if PAA
				    } //if (update singletop in pw)
				} //l pathway counts
			    } // k singletop


			  //PWT: Bewertungsfunktion aufrufen ratio pro pathway berechnen (pw.singletop/pw.counts)
			  //PWT pw-score abspeichern


			  if (pathwayAnalysis == 1)
			    {
			      pathway[r-1].score[thread_nloop] = 0;
			      //PWT_TIM

			      if (pathwayTest==1)
				{
				  pathway[r-1].score[thread_nloop] = snpratio(pathway[r-1],thread_nloop);

				  if (n == 0)
				    {
				      pathway[r-1].ratio0 = snpratio(pathway[r-1],thread_nloop);
				      pathway[r-1].counts0 = pathway[r-1].counts;
				    }
				  //if(r<5) {printf("n%d\tr\t%d\t%f\n",n,r,pathway[r-1].score[thread_nloop]);}
				}

			      else if (pathwayTest==2)
				{
				  pathway[r-1].score[thread_nloop] = fisher2(pathway[r-1],inflationfactor,thread_nloop,lambdaAdjust);

				  if (n == 0)
				    {
				      pathway[r-1].ratio0 = snpratio(pathway[r-1],thread_nloop);
				      pathway[r-1].score0 = fisher2(pathway[r-1],inflationfactor,thread_nloop,lambdaAdjust);
				      pathway[r-1].counts0 = pathway[r-1].counts;
				    }
				}
			      else if (pathwayTest==3)
				{
				  pathway[r-1].score[thread_nloop] = generatio(pathway[r-1],thread_nloop);

				  if (n == 0)
				    {
				      pathway[r-1].ratio0 = generatio(pathway[r-1],thread_nloop);
				      pathway[r-1].counts0 = pathway[r-1].counts;
				    }
				}

			      else if (pathwayTest==4)
				{
				  pathway[r-1].score[thread_nloop] = maxT(pathway[r-1],inflationfactor,thread_nloop,lambdaAdjust);

				  if (n == 0)
				    {
				      pathway[r-1].ratio0 = generatio(pathway[r-1],thread_nloop);
				      pathway[r-1].score0 = maxT(pathway[r-1],inflationfactor,thread_nloop,lambdaAdjust);
				      pathway[r-1].counts0 = pathway[r-1].counts;
				    }

				}
			      else if (pathwayTest==5)
				{
				  //2ndbest
				  skip=1;
				  for(kkk=0;kkk<pathway[r-1].singletop[thread_nloop];kkk++)
				    {
				      if(pathway[r-1].listp_in[thread_nloop][kkk] == 1)
					{
					  for(ll2=0;ll2<pathway[r-1].singletop[thread_nloop];ll2++)
					    {
					      //from same gene:
					      if(pathway[r-1].listsingle_gene[thread_nloop][kkk]==pathway[r-1].listsingle_gene[thread_nloop][ll2] && kkk!=ll2)
						{
						  //do regression
						  p_2nd=1;
						  alt=1;
						  xType=0;
						  helpi=pathway[r-1].listsingle[thread_nloop][kkk];
						  helpj=pathway[r-1].listsingle[thread_nloop][ll2];
						  helpk=0;
						  if(helpi<0 || helpj <0 || helpi >=nlinestped || helpj >=nlinestped)
						    {
						      printf("out\n");exit(1);
						    }
						  if (strcmp(map[helpi].chr, "23") == 0 || strcmp(map[helpj].chr, "23") == 0 )
						    {
						      xType=3;
						    }

						  //resultCopy(&(resultSingle[thread_nloop]),result0);
						  result0 = qtreg(xsinglevec2, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, helpi, helpj, helpk, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, 2,alt,Yhelp, xType, female, male, counts[thread_nloop][helpi], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,n,resultSingle[thread_nloop],0,0,maxIndexCov,liabilityCut,singleMarkerTest,-1,covariancematrix,df_L1_Single,df_L2_Single);

						  if(result0.df>0)
						    {
						      alt=0;
						      result2Single[thread_nloop] = qtreg(zsinglevec2, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, helpi, helpj, helpk, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, 1, alt, Yhelp, xType, female, male, counts[thread_nloop][helpi], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs,PplLocations,n,resultSingle[thread_nloop],0,0,maxIndexCov,liabilityCut,singleMarkerTest,-1,covariancematrix,df_L1_Single,df_L2_Single);
						    }
						  else
						    {
						      result2Single[thread_nloop].df=0;result2Single[thread_nloop].sc=0;
						    }
						  df = result2Single[thread_nloop].df - result0.df;


						  Fstat = ((result2Single[thread_nloop].sc - result0.sc)/(df))/(result0.sc/result0.df);

						  if (df >= 1 && result2Single[thread_nloop].sc > 0.000001)
						    {
						      p_2nd = betai(result0.df/2,df/2,result0.df/(result0.df+df*Fstat));
						    }
						  else
						    {
						      p_2nd = 1;
						    }

						  //if result < cutoff
						  //p_2nd=1;
						  if(p_2nd <= cutoff)
						    {
						      //loop over listp_in2
						      for(mmm=0;mmm<pathway[r-1].singletop[thread_nloop];mmm++)
							{
							  if(pathway[r-1].listp_in2[thread_nloop][mmm]==1)
							    {
							      //vergleich mit besten bisher
							      if(p_2nd < pathway[r-1].listp[thread_nloop][mmm])
								{
								  //printf("here\n");
								  pathway[r-1].listp_in2[thread_nloop][mmm]=0;
								  pathway[r-1].listp_in2[thread_nloop][ll2]=1;
								  pathway[r-1].listp[thread_nloop][ll2]=p_2nd;
								}
							      //update
							      break;
							    }
							}
						      //not found:
						      if(mmm==pathway[r-1].singletop[thread_nloop])
							{
							  pathway[r-1].listp_in2[thread_nloop][ll2]=1;
							  //overwrite
							  pathway[r-1].listp[thread_nloop][ll2]=p_2nd;
							}
						    }
						}
					    }
					}
				    } //2ndbest kkk

				  //3rdbest
				  for(kkk=0;kkk<pathway[r-1].singletop[thread_nloop];kkk++)
				    {
				      if(pathway[r-1].listp_in[thread_nloop][kkk] == 1)
					{
					  //printf("here3 kkk %d\n",kkk);
					  for(jjj=0;jjj<pathway[r-1].singletop[thread_nloop];jjj++)
					    {
					      if(pathway[r-1].listp_in2[thread_nloop][jjj] == 1 && pathway[r-1].listsingle_gene[thread_nloop][kkk]==pathway[r-1].listsingle_gene[thread_nloop][jjj] && jjj!=kkk)
						{
						  for(ll2=0;ll2<pathway[r-1].singletop[thread_nloop];ll2++)
						    {
						      //from same gene:
						      if(pathway[r-1].listsingle_gene[thread_nloop][kkk]==pathway[r-1].listsingle_gene[thread_nloop][ll2] && kkk!=ll2 && jjj!=ll2)
							{
							  //do regression
							  //p_3rd=0.01;

							  alt=1;xType=0;
							  //get index
							  helpi=pathway[r-1].listsingle[thread_nloop][kkk];
							  helpj=pathway[r-1].listsingle[thread_nloop][jjj];
							  helpk=pathway[r-1].listsingle[thread_nloop][ll2];

							  if(helpi<0 || helpj <0 || helpi >=nlinestped || helpj >=nlinestped)
							    {
							      printf("out\n");exit(1);
							    }

							  if(helpi<0 || helpj <0 ||  helpk <0 || helpi >=nlinestped || helpj >=nlinestped || helpk >=nlinestped)
							    {
							      printf("out\n");exit(1);
							    }
							  if (strcmp(map[helpi].chr, "23") == 0 || strcmp(map[helpj].chr, "23") == 0 || strcmp(map[helpk].chr, "23") == 0)
							    {
							      xType=3;
							    }

							  result0 = qtreg(xsinglevec3, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, helpi, helpj, helpk, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, 3,alt,Yhelp, xType, female, male, counts[thread_nloop][helpi], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,n,resultSingle[thread_nloop],0,0,maxIndexCov,liabilityCut,singleMarkerTest,-1,covariancematrix,df_L1_Single,df_L2_Single);

							  if(result0.df>0)
							    {
							      alt=0;
							      result2Single[thread_nloop] = qtreg(zsinglevec3, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, helpi, helpj, helpk, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, 2, alt, Yhelp, xType, female, male, counts[thread_nloop][helpi], Yt, YtX, YtXAinv, YtXAinvXt, 0, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,n,resultSingle[thread_nloop],0,0,maxIndexCov,liabilityCut,singleMarkerTest,-1,covariancematrix,df_L1_Single,df_L2_Single);
							    }
							  else
							    {
							      result2Single[thread_nloop].df=0;result2Single[thread_nloop].sc=0;
							    }

							  df = result2Single[thread_nloop].df - result0.df;
							  Fstat = ((result2Single[thread_nloop].sc - result0.sc)/(df))/(result0.sc/result0.df);

							  if (df >= 1 && result2Single[thread_nloop].sc > 0.000001)
							    {
							      p_3rd = betai(result0.df/2,df/2,result0.df/(result0.df+df*Fstat));
							    }
							  else
							    {
							      p_3rd = 1;
							    }

							  //if result < cutoff
							  if(p_3rd <= cutoff)
							    {
							      //loop over listp_in2
							      for(mmm=0;mmm<pathway[r-1].singletop[thread_nloop];mmm++)
								{
								  if(pathway[r-1].listp_in3[thread_nloop][mmm]==1)
								    {
								      //vergleich mit besten bisher
								      if(p_3rd < pathway[r-1].listp[thread_nloop][mmm])
									{
									  //printf("here\n");
									  pathway[r-1].listp_in3[thread_nloop][mmm]=0;
									  pathway[r-1].listp_in3[thread_nloop][ll2]=1;
									  pathway[r-1].listp[thread_nloop][ll2]=p_3rd;
									}
								      //update
								      break;
								    }
								}
							      //not found:
							      if(mmm==pathway[r-1].singletop[thread_nloop])
								{
								  pathway[r-1].listp_in3[thread_nloop][ll2]=1;
								  //overwrite
								  pathway[r-1].listp[thread_nloop][ll2]=p_3rd;
								}
							    }
							}
						    }  ///ll2
						} //if
					    } //jjj
					}
				    } //3rdbest kkk

				  skip=0;xType=0;alt=1;
				  pathway[r-1].score[thread_nloop] = maxTplus(pathway[r-1],inflationfactor,thread_nloop,cutoff,lambdaAdjust);

				  if (n == 0)
				    {
				      pathway[r-1].ratio0 = generatio(pathway[r-1],thread_nloop);
				      pathway[r-1].score0 = maxTplus(pathway[r-1],inflationfactor,thread_nloop,cutoff,lambdaAdjust);
				      pathway[r-1].counts0 = pathway[r-1].counts;
				    }
				} //pathwaytest 5
			      else if (pathwayTest == 6)
				{
				  if (n == 0)
				    {
				      pathway[r-1].counts0 = pathway[r-1].counts;
				    }

				}

			      if (pathwayTest < 6)
				{
				  Ttable[r-1][n]=pathway[r-1].score[thread_nloop];
				}
			    } //PAA


			  if (pathway[r-1].singletop[thread_nloop] < mWithSingletop)
			    {
			      continue;
			    }
			} //pathwaylists and PAA

		      // Overlapliste in Kombination mit dem Pathway
		      if (mWithSingletop != 0 && mWithGenetictop != 0)
			{
			  for (k = 0; k < overlaptop; k++)
			    {
			      for(l = 0; l < pathway[r-1].counts; l++)
				{
				  if (overlaplist[k].nr == pathway[r-1].list[l])
				    {
				      pathway[r-1].overlaptop[thread_nloop]++;
				      pathway[r-1].listoverlap[thread_nloop] = (int*) realloc(pathway[r-1].listoverlap[thread_nloop], pathway[r-1].overlaptop[thread_nloop]*sizeof(int));
				      if (!pathway[r-1].listoverlap[thread_nloop])
					{
					  errorfile << "memory allocation error in pathway[r-1].listoverlap[thread_nloop]\n";
					  logfile << "memory allocation error in pathway[r-1].listoverlap[thread_nloop]\n";
					  cout << "memory allocation error in pathway[r-1].listlistoverlap[thread_nloop]\n";
					  errorfile.close();logfile.close();exit(1);
					}
				      pathway[r-1].listoverlap[thread_nloop][pathway[r-1].overlaptop[thread_nloop]-1] = overlaplist[k].nr;
				    }
				}
			    }
			  if (pathway[r-1].overlaptop[thread_nloop] < overlapneeded)
			    {
			      continue;
			    }
			}

		      if(markercombi2>=1 || markercombi3>=1)
			{
			  // Prüfen, ob ausgewählte SNPs im Pathway sind
			  snpInPathway1 = 0;
			  snpInPathway2 = 0;
			  snpInPathway3 = 0;

			  for (i = 0; i < pathway[r-1].counts; i++)
			    {
			      if (snp1Pos == pathway[r-1].list[i])
				{
				  snpInPathway1 = 1;
				}
			      else if (snp2Pos == pathway[r-1].list[i])
				{
				  snpInPathway2 = 0;
				}
			      else if (snp3Pos == pathway[r-1].list[i])
				{
				  snpInPathway3 = 0;

				}
			    }
			}
		    } // end r > 0

		  startA = 0;
		  startB = 0;
		  startC = 0;

		  //cout << "mWithSingletop: " << mWithSingletop << "\n";
		  //cout << "mWithGenetictop: " << mWithGenetictop << "\n";



		  if (markercombi3 == 1 || markercombi2 == 1) // 3-Marker
		    {
		      if (markercombi2 == 1)
			{
			  if (r < 1 && n<6)
			    {
			      if (combilist != 1 && haplo==0)
				{
				  if(n==0) {cout << "\nNumber of selected 2-marker-tests: " << ntestsUser << "\n\n";}
				  if(n==0) {logfile << "\nNumber of selected 2-marker-tests: " << ntestsUser << "\n\n";}
				}
			      if (ntestsUser == 0 && n==0)
				{
				  errorfile << "Number of selected 2-marker-tests: " << ntestsUser << "\n\n";
				  logfile << "Number of selected 2-marker-tests: " << ntestsUser << "\n\n";
				  errorfile.close();logfile.close();exit(1);
				}
			    }
			}
		      else if (markercombi3 == 1)
			{
			  if (r < 1 && n< 6)
			    {
			      if (combilist != 1 && haplo==0)
				{

				  if(n==0)
				    {
				      cout << "\nNumber of selected 3-marker-tests: " << ntestsUser << "\n\n";
				      logfile << "\nNumber of selected 3-marker-tests: " << ntestsUser << "\n\n";
				    }
				}
			      if (ntestsUser == 0 && n==0)
				{
				  errorfile << "Number of selected 3-marker-tests: " << ntestsUser << "\n\n";
				  logfile << "Number of selected 3-marker-tests: " << ntestsUser << "\n\n";
				  errorfile.close();logfile.close();exit(1);
				}
			    }
			  oldCombi = 0;
			}

		      if (r < 1)
			{
			  if(combilist == 1)
			    {
			      startA = snp1Pos;
			      stopA = snp1Pos + 1;
			      startB = snp2Pos;
			      stopB = snp2Pos + 1;
			      startC = snp3Pos;
			      stopC = snp3Pos + 1;

			    }
			  else if (selectedSnp == 0)
			    {
			      if (mWithSingletop == 0 && mWithGenetictop == 0)
				{
				  stopA = nlinestped;
				  stopB = nlinestped;
				  stopC = nlinestped;
				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 0)
				{
				  stopA = singletop;
				  stopB = nlinestped;
				  stopC = nlinestped;
				}
			      else if (mWithSingletop == 2 && mWithGenetictop == 0)
				{
				  stopA = singletop;
				  stopB = singletop;
				  stopC = nlinestped;
				}
			      else if (mWithSingletop == 3 && mWithGenetictop == 0)
				{
				  stopA = singletop;
				  stopB = singletop;
				  stopC = singletop;
				}
			      else if (mWithSingletop == 0 && mWithGenetictop == 1)
				{
				  stopA = genetictop;
				  stopB = nlinestped;
				  stopC = nlinestped;
				}
			      else if (mWithSingletop == 0 && mWithGenetictop == 2)
				{
				  stopA = genetictop;
				  stopB = genetictop;
				  stopC = nlinestped;
				}
			      else if (mWithSingletop == 0 && mWithGenetictop == 3)
				{
				  stopA = genetictop;
				  stopB = genetictop;
				  stopC = genetictop;
				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 1)
				{
				  stopA = singletop;
				  stopB = genetictop;
				  stopC = nlinestped;
				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 2)
				{
				  if (markercombi2 == 1)
				    {
				      stopA = genetictop;
				      stopB = overlaptop;
				    }
				  else
				    {
				      stopA = singletop;
				      stopB = genetictop;
				      stopC = genetictop;
				    }
				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 3)

				{
				  stopA = genetictop;
				  stopB = genetictop;
				  stopC = overlaptop;
				}
			      else if (mWithSingletop == 3 && mWithGenetictop == 1)
				{
				  stopA = singletop;
				  stopB = singletop;
				  stopC = overlaptop;
				}
			      else if (mWithSingletop == 2 && mWithGenetictop == 1)
				{
				  if (markercombi2 == 1)
				    {
				      stopA = singletop;
				      stopB = overlaptop;
				    }
				  else
				    {
				      stopA = singletop;
				      stopB = singletop;
				      stopC = genetictop;
				    }
				}
			      else if (mWithSingletop == 2 && mWithGenetictop == 2) // if-Abfrage, ob ein singletop auch genetictop
				{
				  if (markercombi2 == 1)
				    {
				      stopA = overlaptop;
				      stopB = overlaptop;
				    }
				  else
				    {
				      stopA = singletop;
				      stopB = genetictop;
				      stopC = overlaptop;
				    }
				}
			      else if (mWithSingletop == 3 && mWithGenetictop == 2) // if-Abfrage, ob 2 singletop auch genetictop
				{
				  stopA = singletop;
				  stopB = overlaptop;
				  stopC = overlaptop;
				}
			      else if (mWithSingletop == 2 && mWithGenetictop == 3) // if-Abfrage, ob 2 genetictop auch singletop
				{
				  stopA = genetictop;
				  stopB = overlaptop;
				  stopC = overlaptop;
				}
			      else if (mWithSingletop == 3 && mWithGenetictop == 3) // if-Abfrage, ob auch genetictop
				{
				  stopA = overlaptop;
				  stopB = overlaptop;
				  stopC = overlaptop;
				}
			    }
			  else if (selectedSnp == 1)
			    {
			      if (mWithSingletop == 0 && mWithGenetictop == 0)
				{
				  startA = snp1Pos;
				  stopA = snp1Pos + 1;
				  stopB = nlinestped;
				  stopC = nlinestped;
				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 0)
				{
				  startA = snp1Pos;
				  stopA = snp1Pos + 1;
				  stopB = singletop;
				  stopC = nlinestped;
				}
			      else if (mWithSingletop == 2 && mWithGenetictop == 0)
				{
				  startA = snp1Pos;
				  stopA = snp1Pos + 1;
				  stopB = singletop;
				  stopC = singletop;
				}
			      else if (mWithSingletop == 0 && mWithGenetictop == 1)
				{
				  startA = snp1Pos;
				  stopA = snp1Pos + 1;
				  stopB = genetictop;
				  stopC = nlinestped;
				}
			      else if (mWithSingletop == 0 && mWithGenetictop == 2)
				{
				  startA = snp1Pos;
				  stopA = snp1Pos + 1;
				  stopB = genetictop;
				  stopC = genetictop;
				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 1)
				{
				  if (markercombi2 == 1)
				    {
				      startA = snp1Pos;
				      stopA = snp1Pos + 1;
				      stopB = overlaptop;
				    }
				  else
				    {
				      startA = snp1Pos;
				      stopA = snp1Pos + 1;
				      stopB = singletop;
				      stopC = genetictop;
				    }

				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 2) // if-Abfrage, ob einer auch in singletop
				{
				  startA = snp1Pos;
				  stopA = snp1Pos + 1;
				  stopB = genetictop;
				  stopC = overlaptop;
				}
			      else if (mWithSingletop == 2 && mWithGenetictop == 1) // if-Abfrage, ob einer auch in genetictop
				{
				  startA = snp1Pos;
				  stopA = snp1Pos + 1;
				  stopB = singletop;
				  stopC = overlaptop;
				}
			      else if (mWithSingletop == 2 && mWithGenetictop == 2) // if-Abfrage, ob einer auch in genetictop
				{
				  startA = snp1Pos;
				  stopA = snp1Pos + 1;
				  stopB = overlaptop;
				  stopC = overlaptop;
				}
			      else
				{
				  if (markercombi2 == 1)
				    {
				      cout << "Combination is not possible if one SNP is fixed. Cannot selcet more than 1 marker with GENETIC_IMPACT or/and in SINGLETOP\n";
				      logfile << "Combination is not possible if one SNP is fixed. Cannot selcet more than 1 marker with GENETIC_IMPACT or/and in SINGLETOP\n";
				      errorfile << "Combination is not possible if one SNP is fixed. Cannot selcet more than 1 marker with GENETIC_IMPACT or/and in SINGLETOP\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				  else
				    {
				      cout << "Combination is not possible if one SNP is fixed. Cannot selcet more than 2 markers with GENETIC_IMPACT or/and in SINGLETOP\n";
				      logfile << "Combination is not possible if one SNP is fixed. Cannot selcet more than 2 markers with GENETIC_IMPACT or/and in SINGLETOP\n";
				      errorfile << "Combination is not possible if one SNP is fixed. Cannot selcet more than 2 markers with GENETIC_IMPACT or/and in SINGLETOP\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				}
			    }
			  else if (selectedSnp == 2)
			    {
			      if (mWithSingletop == 0 && mWithGenetictop == 0)
				{
				  startA = snp1Pos;
				  stopA = snp1Pos + 1;
				  startB = snp2Pos;
				  stopB = snp2Pos + 1;
				  stopC = nlinestped;
				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 0)
				{
				  startA = snp1Pos;
				  stopA = snp1Pos + 1;
				  startB = snp2Pos;
				  stopB = snp2Pos + 1;
				  stopC = singletop;
				}
			      else if (mWithSingletop == 0 && mWithGenetictop == 1)
				{
				  startA = snp1Pos;
				  stopA = snp1Pos + 1;
				  startB = snp2Pos;
				  stopB = snp2Pos + 1;
				  stopC = genetictop;
				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 1) // Abfrage, ob auch in genetictop
				{
				  startA = snp1Pos;
				  stopA = snp1Pos + 1;
				  startB = snp2Pos;
				  stopB = snp2Pos + 1;
				  stopC = overlaptop;
				}
			      else
				{
				  if (markercombi2 == 1)
				    {
				      cout << "Combination is not possible if 2 SNPs are fixed.\n";
				      logfile << "Combination is not possible if 2 SNPs are fixed.\n";
				      errorfile << "Combination is not possible if 2 SNPs are fixed.\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				  else
				    {
				      cout << "Combination is not possible if 2 SNPs are fixed. Cannot selcet more than 1 marker with GENETIC_IMPACT or/and in SINGLETOP\n";
				      logfile << "Combination is not possible if 2 SNPs are fixed. Cannot selcet more than 1 marker with GENETIC_IMPACT or/and in SINGLETOP\n";
				      errorfile << "Combination is not possible if 2 SNPs are fixed. Cannot selcet more than 1 marker with GENETIC_IMPACT or/and in SINGLETOP\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				}
			    }
			  else if (selectedSnp == 3)
			    {
			      if (mWithSingletop == 0 && mWithGenetictop == 0)
				{

				  startA = snp1Pos;
				  stopA = snp1Pos + 1;
				  startB = snp2Pos;
				  stopB = snp2Pos + 1;
				  startC = snp3Pos;
				  stopC = snp3Pos + 1;
				}
			    }
			  else
			    {
			      cout << "Combination is not possible if 3 SNPs are fixed\n";
			      logfile << "Combination is not possible if 3 SNPs are fixed\n";
			      errorfile << "Combination is not possible if 3 SNPs are fixed\n";
			      errorfile.close();logfile.close();exit(1);
			    }
			}
		      else // mit Pathwayinformation
			{
			  if (selectedSnp == 0)
			    {
			      if (mWithSingletop == 0 && mWithGenetictop == 0)
				{
				  stopA = pathway[r-1].counts;
				  stopB = pathway[r-1].counts;
				  stopC = pathway[r-1].counts;
				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 0)
				{
				  stopA = pathway[r-1].singletop[thread_nloop];
				  stopB = pathway[r-1].counts;
				  stopC = pathway[r-1].counts;
				}
			      else if (mWithSingletop == 2 && mWithGenetictop == 0)
				{
				  stopA = pathway[r-1].singletop[thread_nloop];
				  stopB = pathway[r-1].singletop[thread_nloop];
				  stopC = pathway[r-1].counts;
				}
			      else if (mWithSingletop == 3 && mWithGenetictop == 0)
				{
				  stopA = pathway[r-1].singletop[thread_nloop];
				  stopB = pathway[r-1].singletop[thread_nloop];
				  stopC = pathway[r-1].singletop[thread_nloop];
				}
			      else if (mWithSingletop == 0 && mWithGenetictop == 1)
				{
				  stopA = pathway[r-1].genetictop;
				  stopB = pathway[r-1].counts;
				  stopC = pathway[r-1].counts;
				}
			      else if (mWithSingletop == 0 && mWithGenetictop == 2)
				{
				  stopA = pathway[r-1].genetictop;
				  stopB = pathway[r-1].genetictop;
				  stopC = pathway[r-1].counts;
				}
			      else if (mWithSingletop == 0 && mWithGenetictop == 3)
				{
				  stopA = pathway[r-1].genetictop;
				  stopB = pathway[r-1].genetictop;
				  stopC = pathway[r-1].genetictop;
				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 1)
				{
				  stopA = pathway[r-1].singletop[thread_nloop];
				  stopB = pathway[r-1].genetictop;
				  stopC = pathway[r-1].counts;
				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 2)
				{
				  if (markercombi2 == 1)
				    {
				      stopA = pathway[r-1].genetictop;
				      stopB = pathway[r-1].overlaptop[thread_nloop];
				    }
				  else
				    {
				      stopA = pathway[r-1].singletop[thread_nloop];
				      stopB = pathway[r-1].genetictop;
				      stopC = pathway[r-1].genetictop;
				    }

				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 3)
				{
				  stopA = pathway[r-1].genetictop;
				  stopB = pathway[r-1].genetictop;
				  stopC = pathway[r-1].overlaptop[thread_nloop];
				}
			      else if (mWithSingletop == 3 && mWithGenetictop == 1)
				{
				  stopA = pathway[r-1].singletop[thread_nloop];
				  stopB = pathway[r-1].singletop[thread_nloop];
				  stopC = pathway[r-1].overlaptop[thread_nloop];
				}
			      else if (mWithSingletop == 2 && mWithGenetictop == 1)
				{
				  if (markercombi2 == 1)
				    {
				      stopA = pathway[r-1].singletop[thread_nloop];
				      stopB = pathway[r-1].overlaptop[thread_nloop];
				    }
				  else
				    {
				      stopA = pathway[r-1].singletop[thread_nloop];
				      stopB = pathway[r-1].singletop[thread_nloop];
				      stopC = pathway[r-1].genetictop;
				    }
				}
			      else if (mWithSingletop == 2 && mWithGenetictop == 2) // if-Abfrage, ob ein singletop auch genetictop
				{
				  if (markercombi2 == 1)
				    {
				      stopA = pathway[r-1].overlaptop[thread_nloop];
				      stopB = pathway[r-1].overlaptop[thread_nloop];
				    }
				  else
				    {
				      stopA = pathway[r-1].singletop[thread_nloop];
				      stopB = pathway[r-1].genetictop;
				      stopC = pathway[r-1].overlaptop[thread_nloop];
				    }
				}
			      else if (mWithSingletop == 3 && mWithGenetictop == 2) // if-Abfrage, ob 2 singletop auch genetictop
				{
				  stopA = pathway[r-1].singletop[thread_nloop];
				  stopB = pathway[r-1].overlaptop[thread_nloop];
				  stopC = pathway[r-1].overlaptop[thread_nloop];
				}
			      else if (mWithSingletop == 2 && mWithGenetictop == 3) // if-Abfrage, ob 2 genetictop auch singletop
				{
				  stopA = pathway[r-1].genetictop;
				  stopB = pathway[r-1].overlaptop[thread_nloop];
				  stopC = pathway[r-1].overlaptop[thread_nloop];
				}
			      else if (mWithSingletop == 3 && mWithGenetictop == 3) // if-Abfrage, ob auch genetictop
				{
				  stopA = pathway[r-1].overlaptop[thread_nloop];
				  stopB = pathway[r-1].overlaptop[thread_nloop];
				  stopC = pathway[r-1].overlaptop[thread_nloop];
				}
			    }
			  else if (selectedSnp == 1)
			    {
			      if (mWithSingletop == 0 && mWithGenetictop == 0)
				{
				  if (snpInPathway1 == 1)
				    {
				      startA = snp1Pos;
				      stopA = snp1Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }
				  stopB = pathway[r-1].counts;
				  stopC = pathway[r-1].counts;
				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 0)
				{
				  if (snpInPathway1 == 1)
				    {
				      startA = snp1Pos;
				      stopA = snp1Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;

				    }
				  stopB = pathway[r-1].singletop[thread_nloop];
				  stopC = pathway[r-1].counts;
				}
			      else if (mWithSingletop == 2 && mWithGenetictop == 0)
				{
				  if (snpInPathway1 == 1)
				    {
				      startA = snp1Pos;
				      stopA = snp1Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }
				  stopB = pathway[r-1].singletop[thread_nloop];
				  stopC = pathway[r-1].singletop[thread_nloop];
				}
			      else if (mWithSingletop == 0 && mWithGenetictop == 1)
				{
				  if (snpInPathway1 == 1)
				    {
				      startA = snp1Pos;
				      stopA = snp1Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }
				  stopB = pathway[r-1].genetictop;
				  stopC = pathway[r-1].counts;
				}
			      else if (mWithSingletop == 0 && mWithGenetictop == 2)
				{
				  if (snpInPathway1 == 1)
				    {
				      startA = snp1Pos;
				      stopA = snp1Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }
				  stopB = pathway[r-1].genetictop;
				  stopC = pathway[r-1].genetictop;
				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 1)
				{
				  if (snpInPathway1 == 1)
				    {
				      startA = snp1Pos;
				      stopA = snp1Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }

				  if (markercombi2 == 1)
				    {
				      stopB = pathway[r-1].overlaptop[thread_nloop];
				    }
				  else
				    {
				      stopB = pathway[r-1].singletop[thread_nloop];
				      stopC = pathway[r-1].genetictop;
				    }
				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 2) // if-Abfrage, ob einer auch in singletop
				{
				  if (snpInPathway1 == 1)
				    {
				      startA = snp1Pos;
				      stopA = snp1Pos + 1;
				    }
				  else
				    {

				      snpInPathway = 0;
				      continue;
				    }
				  stopB = pathway[r-1].genetictop;
				  stopC = pathway[r-1].overlaptop[thread_nloop];
				}
			      else if (mWithSingletop == 2 && mWithGenetictop == 1) // if-Abfrage, ob einer auch in genetictop
				{
				  if (snpInPathway1 == 1)
				    {
				      startA = snp1Pos;
				      stopA = snp1Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }
				  stopB = pathway[r-1].singletop[thread_nloop];
				  stopC = pathway[r-1].overlaptop[thread_nloop];

				}
			      else if (mWithSingletop == 2 && mWithGenetictop == 2) // if-Abfrage, ob einer auch in genetictop
				{
				  if (snpInPathway1 == 1)
				    {
				      startA = snp1Pos;
				      stopA = snp1Pos + 1;
				    }
				  else

				    {
				      snpInPathway = 0;
				      continue;
				    }
				  stopB = pathway[r-1].overlaptop[thread_nloop];
				  stopC = pathway[r-1].overlaptop[thread_nloop];
				}
			      else
				{
				  if (markercombi2 == 1)
				    {
				      cout << "Combination is not possible if 1 SNP is fixed. Cannot selcet more than 1 marker with GENETIC_IMPACT or/and in SINGLETOP\n";
				      logfile << "Combination is not possible if 1 SNP is fixed. Cannot selcet more than 1 marker with GENETIC_IMPACT or/and in SINGLETOP\n";
				      errorfile << "Combination is not possible if 1 SNP is fixed. Cannot selcet more than 1 marker with GENETIC_IMPACT or/and in SINGLETOP\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				  else
				    {
				      cout << "Combination is not possible if 1 SNP is fixed. Cannot selcet more than 2 markers with GENETIC_IMPACT or/and in SINGLETOP\n";
				      logfile << "Combination is not possible if 1 SNP is fixed. Cannot selcet more than 2 markers with GENETIC_IMPACT or/and in SINGLETOP\n";
				      errorfile << "Combination is not possible if 1 SNP is fixed. Cannot selcet more than 2 markers with GENETIC_IMPACT or/and in SINGLETOP\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				}
			    }
			  else if (selectedSnp == 2)
			    {
			      if (mWithSingletop == 0 && mWithGenetictop == 0)
				{
				  if (snpInPathway1 == 1)
				    {
				      startA = snp1Pos;
				      stopA = snp1Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }
				  if (snpInPathway2 == 1)
				    {
				      startB = snp2Pos;
				      stopB = snp2Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }
				  stopC = pathway[r-1].counts;
				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 0)
				{
				  if (snpInPathway1 == 1)
				    {
				      startA = snp1Pos;
				      stopA = snp1Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }
				  if (snpInPathway2 == 1)
				    {
				      startB = snp2Pos;
				      stopB = snp2Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }
				  stopC = pathway[r-1].singletop[thread_nloop];
				}
			      else if (mWithSingletop == 0 && mWithGenetictop == 1)
				{
				  if (snpInPathway1 == 1)
				    {
				      startA = snp1Pos;
				      stopA = snp1Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }
				  if (snpInPathway2 == 1)
				    {
				      startB = snp2Pos;
				      stopB = snp2Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }
				  stopC = pathway[r-1].genetictop;
				}
			      else if (mWithSingletop == 1 && mWithGenetictop == 1) // Abfrage, ob auch in genetictop
				{
				  if (snpInPathway1 == 1)
				    {
				      startA = snp1Pos;
				      stopA = snp1Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }
				  if (snpInPathway2 == 1)
				    {
				      startB = snp2Pos;
				      stopB = snp2Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }
				  stopC = pathway[r-1].overlaptop[thread_nloop];
				}
			      else
				{
				  if (markercombi2 == 1)
				    {
				      cout << "Combination is not possible if 2 SNPs are fixed.\n";
				      logfile << "Combination is not possible if 2 SNPs are fixed.\n";
				      errorfile << "Combination is not possible if 2 SNPs are fixed.\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				  else
				    {
				      cout << "Combination is not possible if 2 SNPs are fixed. Cannot selcet more than 1 marker with GENETIC_IMPACT or/and in SINGLETOP\n";
				      logfile << "Combination is not possible if 2 SNPs are fixed. Cannot selcet more than 1 marker with GENETIC_IMPACT or/and in SINGLETOP\n";
				      errorfile << "Combination is not possible if 2 SNPs are fixed. Cannot selcet more than 1 marker with GENETIC_IMPACT or/and in SINGLETOP\n";
				      errorfile.close();logfile.close();exit(1);
				    }
				}
			    }
			  else if (selectedSnp == 3)
			    {
			      if (mWithSingletop == 0 && mWithGenetictop == 0)
				{
				  if (snpInPathway1 == 1)
				    {
				      startA = snp1Pos;
				      stopA = snp1Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }
				  if (snpInPathway2 == 1)
				    {
				      startB = snp2Pos;
				      stopB = snp2Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }
				  if (snpInPathway3 == 1)
				    {
				      startC = snp3Pos;
				      stopC = snp3Pos + 1;
				    }
				  else
				    {
				      snpInPathway = 0;
				      continue;
				    }
				}
			      else
				{
				  cout << "Combination is not possible if 3 SNPs are fixed\n";
				  logfile << "Combination is not possible if 3 SNPs are fixed\n";
				  errorfile << "Combination is not possible if 3 SNPs are fixed\n";
				  errorfile.close();logfile.close();exit(1);
				}
			    }
			}



		      if (markercombi2 == 1)
			{
			  startC = 0;
			  stopC = 1;
			}

		      //cout << "stopA: " << stopA << "\n";
		      //cout << "stopB: " << stopB << "\n";
		      //cout << "stopC: " << stopC << "\n";
		      //Überprüfen, ob stopA, stopB und stopC nicht 0 sind d.h. die augewählten Listen leer
		      if (stopA == 0 || stopB == 0 || stopC == 0)
			{
			  continue;
			}
		    }

		  a=0;b=0;c=0;

		  if(regression && plla2 && r <=1) //TIM 1.0.9
		    {
		      //delete stuff
		      free(p);p=NULL;free(S);S=NULL;
		      free(YY);YY=NULL;free(Yhelp);Yhelp = NULL;

		      free2Ddouble(X,nlinestfam);free2Ddouble(Xmod,nlinestfam);
		      free2Ddouble(Xt,dim1+maxIndexCov);free2Ddouble(A,dim1+maxIndexCov);
		      free2Ddouble(VNN,dim1+maxIndexCov);free2Ddouble(Sinv,dim1+maxIndexCov);
		      free2Ddouble(A0,dim1+maxIndexCov);free2Ddouble(UNNT,dim1+maxIndexCov);
		      free2Ddouble(Ainv,dim1+maxIndexCov);free2Ddouble(AinvXt,dim1+maxIndexCov);
		      free2Ddouble(Yminusp,nlinestfam);free2Ddouble(newbeta,dim1+maxIndexCov);
		      free2Ddouble(D,dim1+maxIndexCov);free2Ddouble(T,dim1+maxIndexCov);
		      free2Ddouble(U,dim1+maxIndexCov);free2Ddouble(Ut,dim1+maxIndexCov);
		      free2Ddouble(sumPP,dim1+maxIndexCov);free2Ddouble(sumPJ,dim1+maxIndexCov);
		      free2Ddouble(sumPK,dim1+maxIndexCov);free2Ddouble(MMinv,dim1+maxIndexCov);

		      if(qt)
			{
			  free2Ddouble(Yt,1);free2Ddouble(YtX,1);
			  free2Ddouble(YtXAinv,1);free2Ddouble(YtXAinvXt,1);
			}
		    }

#if PARALLELA
# pragma omp parallel // if (nsim == 0)
		  {
		    srand(int(time(NULL)) ^ omp_get_thread_num());

# pragma omp for schedule(dynamic,100) private(aMod,bMod,cMod,thread_aloop,traitavg,i,b2,c2,d,e,countA,startB,a,a_list,countB,b,startC,b_list,firstc,countC,c,c_list,d_list,e_list,x,oldCombi,ifChrX,c_in,d_in,e_in,tstat,stopCC,stopDD,stopEE,aa,bb,cc,dd,ee,controlcounts5,controlcounts5M,controlcounts5F,casecounts5,casecounts5M,casecounts5F,index,alt,xType,tLog,df,Fstat,tLogMale,dfMale,tLogFemale,dfFemale,marker1,marker2,marker3,ii,jj,kk,ll,mm,p,S,YY,Yhelp,X,Xmod,Xt,A,VNN,Sinv,A0,UNNT,Ainv,AinvXt,Yminusp,newbeta,D,T,U,Ut,sumPP,sumPJ,sumPK,MMinv,Yt,YtX,YtXAinv,YtXAinvXt,ll2,stat,tstat2,filled,passedPreTest)

#endif
		    for (a2 = startA; a2 < stopA; a2++)
		      {
			thread_aloop=0;
			if(plln){thread_aloop=thread_nloop;}
#if PARALLELA
			thread_aloop=omp_get_thread_num();
			if(thread_aloop > maxthreads){printf("error\n");exit(1);}
			//YY=NULL;
#endif
			filled=1;passedPreTest=1;

			if(pathwayAnalysis && pathwayTest < 6){cout << "error_code45\n";exit(1);}
			oldCombi=0;
			if (regression == 1 && (plla2) ) //TIM 1.0.9
			  {
			    p = (double *) calloc(nlinestfam,sizeof(double));
			    if (!p)
			      {
				errorfile << "memory allocation error in p\n";
				logfile << "memory allocation error in p\n";

				cout << "memory allocation error in p\n";
				errorfile.close();logfile.close();exit(1);
			      }

			    S = (double *) calloc(dim1+maxIndexCov,sizeof(double));
			    if (!S)
			      {
				errorfile << "memory allocation error in S\n";
				logfile << "memory allocation error in S\n";
				cout << "memory allocation error in S\n";
				errorfile.close();logfile.close();exit(1);
			      }
			    YY = (double *) calloc(nlinestfam,sizeof(double));
			    if (!YY)
			      {
				errorfile << "memory allocation error in YY\n";
				logfile << "memory allocation error in YY\n";
				cout << "memory allocation error in YY\n";
				errorfile.close();logfile.close();exit(1);
			      }
			    Yhelp = (double *) calloc(nlinestfam,sizeof(double));
			    if (!Yhelp)
			      {
				errorfile << "memory allocation error in Yhelp\n";
				logfile << "memory allocation error in Yhelp\n";
				cout << "memory allocation error in Yhelp\n";
				errorfile.close();logfile.close();exit(1);
			      }

			    X = calloc2Ddouble(nlinestfam, dim1+maxIndexCov);
			    Xmod = calloc2Ddouble(nlinestfam, dim1+maxIndexCov);
			    Xt = calloc2Ddouble(dim1+maxIndexCov, nlinestfam);
			    A = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
			    VNN = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
			    Sinv = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
			    A0 = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
			    UNNT = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
			    Ainv = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
			    AinvXt = calloc2Ddouble(dim1+maxIndexCov, nlinestfam);
			    Yminusp = calloc2Ddouble(nlinestfam, 1);
			    newbeta = calloc2Ddouble(dim1+maxIndexCov, 1);
			    D = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
			    T = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
			    U = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
			    Ut = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
			    sumPP = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
			    sumPJ = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
			    sumPK = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);
			    MMinv = calloc2Ddouble(dim1+maxIndexCov, dim1+maxIndexCov);

			    if(qt)
			      {
				Yt = calloc2Ddouble(1,nlinestfam);
				YtX = calloc2Ddouble(1, dim1+maxIndexCov);
				YtXAinv = calloc2Ddouble(1, dim1+maxIndexCov);
				YtXAinvXt = calloc2Ddouble(1, nlinestfam);
			      }
			  } //plla2


			countA++;
			if (selectedSnp == 0)
			  {
			    if (mWithSingletop == 0)
			      {
				if (mWithGenetictop == 0)
				  {
				    startB = a2+1;
				    if (r == 0)
				      {
					a = a2;
				      }
				    else
				      {
					a = pathway[r-1].list[a2];
				      }
				  }
				else if (mWithGenetictop == 1)
				  {
				    startB = 0;
				    if (r == 0)
				      {
					a = geneticlist[a2].nr;
				      }
				    else
				      {
					a = pathway[r-1].listgenetic[a2];
				      }
				  }
				else
				  {
				    startB = a2+1;
				    if (r == 0)
				      {
					a = geneticlist[a2].nr;
				      }
				    else
				      {
					a = pathway[r-1].listgenetic[a2];
				      }
				  }
			      }
			    else if (mWithSingletop == 1)
			      {
				if (mWithGenetictop == 2)
				  {
				    startB = 0;
				    if (r == 0)
				      {
					a = geneticlist[a2].nr;
				      }
				    else
				      {
					a = pathway[r-1].listgenetic[a2];
				      }
				  }
				else if (mWithGenetictop == 3)
				  {
				    startB = a2+1;
				    if (r == 0)
				      {
					a = geneticlist[a2].nr;
				      }
				    else
				      {
					a = pathway[r-1].listgenetic[a2];
				      }
				  }
				else
				  {
				    startB = 0;
				    if (r == 0)
				      {
					a = bestsinglemarker[a2].nr;
				      }
				    else
				      {
					a = pathway[r-1].listsingle[thread_nloop][a2];
				      }

				  }
			      }
			    else if (mWithSingletop == 2)
			      {
				if (mWithGenetictop == 0 || (mWithGenetictop == 1 && markercombi3 == 1))
				  {
				    startB = a2+1;
				    if (r == 0)
				      {
					a = bestsinglemarker[a2].nr;
				      }
				    else
				      {
					a = pathway[r-1].listsingle[thread_nloop][a2];
				      }
				  }
				else if ((mWithGenetictop == 1 && markercombi2 == 1) || (mWithGenetictop == 2 && markercombi3 == 1))
				  {
				    startB = 0;
				    if (r == 0)
				      {
					a = bestsinglemarker[a2].nr;
				      }
				    else
				      {
					a = pathway[r-1].listsingle[thread_nloop][a2];
				      }
				  }
				else if (mWithGenetictop == 2 && markercombi2 == 1)
				  {
				    startB = a2+1;
				    if (r == 0)
				      {
					a = overlaplist[a2].nr;
				      }
				    else
				      {
					a = pathway[r-1].listoverlap[thread_nloop][a2];
				      }

				  }

				else if (mWithGenetictop == 3)
				  {
				    startB = 0;
				    if (r == 0)
				      {
					a = geneticlist[a2].nr;
				      }
				    else
				      {
					a = pathway[r-1].listgenetic[a2];
				      }
				  }
			      }
			    else if (mWithSingletop == 3)
			      {
				if (mWithGenetictop < 2 )
				  {
				    startB = a2+1;
				    if (r == 0)
				      {
					a = bestsinglemarker[a2].nr;
				      }
				    else
				      {
					a = pathway[r-1].listsingle[thread_nloop][a2];
				      }
				  }
				else if (mWithGenetictop == 2)
				  {
				    startB = 0;
				    if (r == 0)
				      {
					a = bestsinglemarker[a2].nr;
				      }
				    else
				      {
					a = pathway[r-1].listsingle[thread_nloop][a2];
				      }
				  }
				else if (mWithGenetictop == 3)
				  {
				    startB = a2+1;
				    if (r == 0)
				      {
					a = overlaplist[a2].nr;
				      }
				    else
				      {
					a = pathway[r-1].listoverlap[thread_nloop][a2];
				      }
				  }
			      }
			  }
			else if (selectedSnp == 2 || selectedSnp == 3)
			  {
			    if (r == 0)
			      {
				a = snp1Pos;
			      }
			  }
			else
			  {
			    startB = 0;
			    if (r == 0)
			      {
				a = a2;
			      }
			    else
			      {
				a = pathway[r-1].list[a2];
			      }
			  }
			a_list = a;

			for (b2 = startB; b2 < stopB; b2++)
			  {
			    countB++;
			    if (selectedSnp == 0)
			      {
				if (mWithSingletop == 0)
				  {
				    if (mWithGenetictop < 2)
				      {
					startC = b2+1;
					if (r == 0)
					  {
					    b = b2;
					  }
					else
					  {
					    b = pathway[r-1].list[b2];
					  }
				      }
				    else if (mWithGenetictop == 2)
				      {
					startC = 0;
					if (r == 0)
					  {
					    b = geneticlist[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listgenetic[b2];
					  }
				      }
				    else if (mWithGenetictop == 3)
				      {
					startC = b2+1;
					if (r == 0)
					  {
					    b = geneticlist[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listgenetic[b2];
					  }
				      }
				  }
				else if (mWithSingletop == 1)
				  {
				    if (mWithGenetictop == 0 && markercombi2 == 1)
				      {
					startC = 0;
					if (r == 0)
					  {
					    b = b2;
					  }
					else
					  {
					    b = pathway[r-1].list[b2];
					  }
				      }
				    if (mWithGenetictop == 0 && markercombi3 == 1)
				      {
					startC = b2+1;
					if (r == 0)
					  {
					    b = b2;
					  }

					else
					  {
					    b = pathway[r-1].list[b2];
					  }
				      }
				    else if (mWithGenetictop == 1 || mWithGenetictop == 3)
				      {
					startC = 0;
					if (r == 0)
					  {
					    b = geneticlist[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listgenetic[b2];
					  }
				      }
				    else if (mWithGenetictop == 2 && markercombi2 == 1)
				      {
					startC = 0;
					if (r == 0)

					  {
					    b = overlaplist[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listoverlap[thread_nloop][b2];
					  }
				      }
				    else if (mWithGenetictop == 2 && markercombi3 == 1 )
				      {
					startC = b2+1;
					if (r == 0)
					  {
					    b = geneticlist[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listgenetic[b2];
					  }
				      }
				  }
				else if (mWithSingletop == 2)
				  {
				    if (mWithGenetictop == 0)
				      {
					startC = 0;
					if (r == 0)
					  {
					    b = bestsinglemarker[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listsingle[thread_nloop][b2];
					  }
				      }
				    else if ((mWithGenetictop == 1 && markercombi2 == 1)  || (mWithGenetictop == 2 && markercombi2 == 1))
				      {
					startC = 0;
					if (r == 0)
					  {
					    b = overlaplist[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listoverlap[thread_nloop][b2];
					  }
				      }
				    else if (mWithGenetictop == 1 && markercombi3 == 1)
				      {
					startC = 0;
					if (r == 0)
					  {
					    b = bestsinglemarker[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listsingle[thread_nloop][b2];
					  }
				      }
				    else if (mWithGenetictop == 2 && markercombi3 == 1)
				      {
					startC = 0;
					if (r == 0)
					  {
					    b = geneticlist[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listgenetic[b2];
					  }
				      }
				    else if (mWithGenetictop == 3)
				      {
					startC = b2+1;
					if (r == 0)
					  {
					    b = overlaplist[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listoverlap[thread_nloop][b2];
					  }
				      }
				  }
				else if (mWithSingletop == 3)
				  {
				    if (mWithGenetictop == 0)
				      {
					startC = b2+1;
					if (r == 0)
					  {
					    b = bestsinglemarker[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listsingle[thread_nloop][b2];
					  }
				      }
				    else if (mWithGenetictop == 1)
				      {
					startC = 0;
					if (r == 0)
					  {
					    b = bestsinglemarker[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listsingle[thread_nloop][b2];
					  }
				      }
				    else if (mWithGenetictop == 2 || mWithGenetictop == 3)
				      {
					startC = b2+1;
					if (r == 0)
					  {
					    b = overlaplist[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listoverlap[thread_nloop][b2];
					  }
				      }
				  }
			      }
			    else if (selectedSnp == 1)
			      {
				if (mWithSingletop == 0)
				  {
				    if (mWithGenetictop == 0)
				      {
					startC = b2+1;
					if (r == 0)
					  {
					    b = b2;
					  }
					else
					  {
					    b = pathway[r-1].list[b2];
					  }

				      }
				    else if (mWithGenetictop == 1)
				      {
					startC = 0;
					if (r == 0)
					  {
					    b = geneticlist[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listgenetic[b2];
					  }
				      }
				    else if (mWithGenetictop == 2)
				      {
					startC = b2+1;
					if (r == 0)
					  {
					    b = geneticlist[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listgenetic[b2];
					  }
				      }
				  }
				else if (mWithSingletop == 1)
				  {
				    if ((mWithGenetictop < 2 && markercombi3 == 1) || (mWithGenetictop == 0 && markercombi2 == 1))
				      {
					startC = 0;

					if (r == 0)
					  {
					    b = bestsinglemarker[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listsingle[thread_nloop][b2];
					  }
				      }
				    else if (mWithGenetictop == 1 && markercombi2 == 1)
				      {
					startC = 0;
					if (r == 0)
					  {
					    b = overlaplist[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listoverlap[thread_nloop][b2];
					  }
				      }
				    else if (mWithGenetictop == 2)
				      {
					startC = 0;
					if (r == 0)
					  {
					    b = geneticlist[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listgenetic[b2];
					  }
				      }

				  }
				else if (mWithSingletop == 2)
				  {
				    if (mWithGenetictop == 0)
				      {
					startC = b2+1;
					if (r == 0)
					  {
					    b = bestsinglemarker[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listsingle[thread_nloop][b2];
					  }
				      }
				    else if (mWithGenetictop == 1)
				      {
					startC = 0;
					if (r == 0)
					  {
					    b = bestsinglemarker[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listsingle[thread_nloop][b2];
					  }
				      }
				    else if (mWithGenetictop == 2)
				      {
					startC = b2+1;
					if (r == 0)
					  {
					    b = overlaplist[b2].nr;
					  }
					else
					  {
					    b = pathway[r-1].listoverlap[thread_nloop][b2];
					  }
				      }
				  }
			      }
			    else if (selectedSnp == 2)
			      {
				startC = 0;
				if (r == 0)
				  {
				    b = snp2Pos;
				  }
			      }
			    else if (selectedSnp == 3)
			      {
				if (r == 0)
				  {
				    b = snp2Pos;
				  }
			      }

			    b_list = b;

			    firstc = 0;

			    if (markercombi2 == 1)
			      {
				startC = 0;
			      }


			    if(haplo && (markercombi2>=1 || markercombi3>=1) && !combilist)
			      {
				if(strcmp(map[a_list].chr,map[b_list].chr)!=0
				   || (strcmp(map[a_list].chr,map[b_list].chr)==0 && ((map[b_list].pos - map[a_list].pos) > distance)))
				  {
				    break;
				  }
			      }

			    for (c2 = startC; c2 < stopC; c2++)
			      {
				if (markercombi3 == 1)
				  {
				    countC++;
				    if (selectedSnp == 0)
				      {
					if (mWithSingletop == 0)
					  {
					    if (mWithGenetictop < 3)
					      {
						if (r == 0)
						  {
						    c = c2;
						  }
						else
						  {
						    c = pathway[r-1].list[c2];
						  }
					      }
					    else if (mWithGenetictop == 3)
					      {
						if (r == 0)
						  {
						    c = geneticlist[c2].nr;
						  }
						else
						  {
						    c = pathway[r-1].listgenetic[c2];
						  }
					      }
					  }
					else if (mWithSingletop == 1)
					  {
					    if (mWithGenetictop < 2)
					      {
						if (r == 0)
						  {
						    c = c2;
						  }
						else
						  {
						    c = pathway[r-1].list[c2];
						  }
					      }
					    else if (mWithGenetictop == 2)
					      {
						if (r == 0)
						  {
						    c = geneticlist[2].nr;
						  }
						else
						  {
						    c = pathway[r-1].listgenetic[c2];
						  }
					      }
					    else if (mWithGenetictop == 3)
					      {
						if (r == 0)
						  {
						    c = overlaplist[c2].nr;
						  }
						else
						  {
						    c = pathway[r-1].listoverlap[thread_nloop][c2];
						  }
					      }
					  }
					else if (mWithSingletop == 2)
					  {
					    if (mWithGenetictop == 0)
					      {
						if (r == 0)
						  {
						    c = c2;
						  }
						else
						  {
						    c = pathway[r-1].list[c2];
						  }
					      }
					    else if (mWithGenetictop == 1)
					      {
						if (r == 0)
						  {
						    c = geneticlist[c2].nr;
						  }
						else
						  {
						    c = pathway[r-1].listgenetic[c2];
						  }
					      }
					    else if (mWithGenetictop > 1)
					      {
						if (r == 0)
						  {
						    c = overlaplist[c2].nr;
						  }
						else
						  {
						    c = pathway[r-1].listoverlap[thread_nloop][c2];
						  }
					      }
					  }
					else if (mWithSingletop == 3)
					  {
					    if (mWithGenetictop == 0)
					      {
						if (r == 0)
						  {
						    c = bestsinglemarker[c2].nr;
						  }
						else
						  {
						    c = pathway[r-1].listsingle[thread_nloop][c2];
						  }
					      }
					    else if (mWithGenetictop > 0)
					      {
						if (r == 0)
						  {
						    c = overlaplist[c2].nr;
						  }
						else
						  {
						    c = pathway[r-1].listoverlap[thread_nloop][c2];
						  }
					      }
					  }
				      }
				    else if (selectedSnp == 1)
				      {
					if (mWithSingletop == 0)
					  {
					    if (mWithGenetictop < 2)
					      {
						if (r == 0)
						  {
						    c = c2;
						  }
						else
						  {
						    c = pathway[r-1].list[c2];
						  }
					      }
					    else if (mWithGenetictop == 2)
					      {
						if (r == 0)
						  {
						    c = geneticlist[c2].nr;
						  }
						else
						  {
						    c = pathway[r-1].listgenetic[c2];
						  }
					      }
					  }
					else if (mWithSingletop == 1)
					  {
					    if (mWithGenetictop == 0)
					      {
						if (r == 0)
						  {
						    c = c2;
						  }
						else
						  {
						    c = pathway[r-1].list[c2];
						  }
					      }
					    else if (mWithGenetictop == 1)
					      {
						if (r == 0)
						  {
						    c = geneticlist[c2].nr;
						  }
						else
						  {
						    c = pathway[r-1].listgenetic[c2];
						  }
					      }
					    else if (mWithGenetictop == 2)
					      {
						if (r == 0)
						  {
						    c = overlaplist[c2].nr;
						  }
						else
						  {
						    c = pathway[r-1].listoverlap[thread_nloop][c2];
						  }
					      }
					  }
					else if (mWithSingletop == 2)
					  {
					    if (mWithGenetictop == 0)
					      {
						if (r == 0)
						  {
						    c = bestsinglemarker[c2].nr;
						  }
						else
						  {
						    c = pathway[r-1].listsingle[thread_nloop][c2];
						  }
					      }
					    else if (mWithGenetictop > 0)
					      {
						if (r == 0)
						  {
						    c = overlaplist[c2].nr;
						  }
						else
						  {
						    c = pathway[r-1].listoverlap[thread_nloop][c2];
						  }
					      }
					  }
				      }
				    else if (selectedSnp == 2)
				      {
					if (mWithSingletop == 0)
					  {
					    if (mWithGenetictop == 0)
					      {
						if (r == 0)
						  {
						    c = c2;
						  }
						else
						  {
						    c = pathway[r-1].list[c2];
						  }
					      }
					    else if (mWithGenetictop == 1)
					      {
						if (r == 0)
						  {
						    c = geneticlist[c2].nr;
						  }
						else
						  {
						    c = pathway[r-1].listgenetic[c2];
						  }
					      }
					  }
					else if (mWithSingletop == 1)
					  {
					    if (mWithGenetictop == 0)
					      {
						if (r == 0)
						  {
						    c = bestsinglemarker[c2].nr;
						  }
						else
						  {
						    c = pathway[r-1].listsingle[thread_nloop][c2];
						  }
					      }
					    else if (mWithGenetictop == 1)
					      {
						if (r == 0)
						  {
						    c = overlaplist[c2].nr;
						  }
						else
						  {
						    c = pathway[r-1].listoverlap[thread_nloop][c2];
						  }
					      }
					  }
				      }
				    else if (selectedSnp == 3)
				      {
					if (r == 0)
					  {
					    c = snp3Pos;
					  }
				      }
				    c_list = c;
				  }
				else
				  {
				    c_list = 0;c=0;
				  }

				if(haplo && (markercombi3>=1) && !combilist)
				  {
				    if(strcmp(map[a_list].chr,map[b_list].chr)!=0
				       || (strcmp(map[a_list].chr,map[b_list].chr)==0 && ((map[b_list].pos - map[a_list].pos) > distance)))
				      {
					break;
				      }


				    if(strcmp(map[a_list].chr,map[c_list].chr)!=0 || (strcmp(map[a_list].chr,map[c_list].chr)==0 && ((map[c_list].pos - map[a_list].pos) > distance)))
				      {
					break;
				      }
				  }

				for (d = 0; d < 1; d++)
				  {
				    d_list = d;
				    for (e = 0; e < 1; e++)
				      {
					e_list = e;

					if (r==0)
					  {
					    oldCombi=0;
					    if (markercombi2 == 1 && !combilist)
					      {
#pragma omp critical (OLDCOMBI)
						{
						  oldCombi=0;
						  if ( !((mWithSingletop==0 && mWithGenetictop==0) || ((mWithSingletop==2 && mWithGenetictop==0))))
						    {
						      for (x = 0; x < printtop; x++)
							{
							  if (bestchi3[x].nr1 == b_list && bestchi3[x].nr2 == a_list)
							    {
							      oldCombi = 1;
							      break;
							    }
							}
						    }

						  if(mWithSingletop==1 && mWithGenetictop==0)
						    {
						      for (x = 0; x < singletop; x++)
							{
							  if (bestsinglemarker[x].nr == b_list)
							    {
							      if(selectedSnp ==0 && (bestsinglemarker[a2].nr > b_list))
								{
								  oldCombi = 1;
								  break;
								}
							    }
							}
						    }
						}
					      }
					    else if (!combilist)
					      {
#pragma omp critical (OLDCOMBI)
						{
						  oldCombi=0;
						  for (x = 0; x < printtop; x++)
						    {
						      if ((bestchi3[x].nr1 == a_list && bestchi3[x].nr2 == b_list && bestchi3[x].nr3 == c_list)
							  || (bestchi3[x].nr1 == a_list && bestchi3[x].nr2 == c_list && bestchi3[x].nr3 == b_list)
							  || (bestchi3[x].nr1 == b_list && bestchi3[x].nr2 == a_list && bestchi3[x].nr3 == c_list)
							  || (bestchi3[x].nr1 == b_list && bestchi3[x].nr2 == c_list && bestchi3[x].nr3 == a_list)
							  || (bestchi3[x].nr1 == c_list && bestchi3[x].nr2 == a_list && bestchi3[x].nr3 == b_list)
							  || (bestchi3[x].nr1 == c_list && bestchi3[x].nr2 == b_list && bestchi3[x].nr3 == a_list))
							{
							  oldCombi = 1;
							  break;
							}
						    }
						}
					      }
					  }


					if (markercombi2 == 1)
					  {
					    //NEW3
					    c_in=1;d_in=1;e_in=1;

					  }
					else
					  {
					    //NEW3
					    c_in=map[c_list].analysis_in;
					    d_in=1;e_in=1;
					  }

					if ((strcmp(map[a_list].chr, "23") == 0)
					    || (strcmp(map[b_list].chr, "23") == 0)
					    || (strcmp(map[c_list].chr, "23") == 0))
					  {
					    ifChrX = 1;
					  }
					else
					  {
					    ifChrX = 0;
					  }

					//NEW3
					if (map[a_list].analysis_in == 1 && map[b_list].analysis_in == 1 && c_in == 1 && d_in == 1 && e_in == 1 && oldCombi != 1)
					  {
					    if ((markercombi2 == 1 && a_list != b_list) || (markercombi3 == 1 && a_list != b_list && a_list != c_list && b_list != c_list))
					      {
						tstat.p = 1;
						tstat.pmod = 1;

						//CASEONLY_18
						tstat2.p = 1;
						tstat2.pmod = 1;

						stopCC = 3;
						stopDD = 1;
						stopEE = 1;


						for (aa = 0; aa < 3; aa++)
						  {
						    for (bb = 0; bb < 3; bb++)
						      {
							for (cc = 0; cc < stopCC ; cc++) // Achtung Änderung!!!!!! für Marker 2
							  {
							    for (dd = 0; dd < stopDD; dd++)
							      {
								for (ee = 0; ee < stopEE; ee++)
								  {
								    controlcounts5[aa][bb][cc][dd][ee] = 0;
								    controlcounts5M[aa][bb][cc][dd][ee] = 0;
								    controlcounts5F[aa][bb][cc][dd][ee] = 0;
								    casecounts5[aa][bb][cc][dd][ee] = 0;
								    casecounts5M[aa][bb][cc][dd][ee] = 0;
								    casecounts5F[aa][bb][cc][dd][ee] = 0;

								  }
							      }
							  }
						      }
						  }

						//int aMod;
						//int bMod;
						//int cMod;

						//CASE_ONLY
						if(caseOnly)
						  {
						    if((test == 15 || test == 16 || test == 17 || test == 18) && (strcmp(map[a_list].chr,map[b_list].chr)!=0 || (strcmp(map[a_list].chr,map[b_list].chr)==0 && (abs(map[b_list].pos - map[a_list].pos) > LDdistance))))
						      {

						      }
						    else if((test == 19) && (strcmp(map[a_list].chr,map[b_list].chr)!=0 || (strcmp(map[a_list].chr,map[b_list].chr)==0 && (abs(map[b_list].pos - map[a_list].pos) > LDdistance)))
							    && (strcmp(map[a_list].chr,map[c_list].chr)!=0 || (strcmp(map[a_list].chr,map[c_list].chr)==0&& (abs(map[c_list].pos - map[a_list].pos) > LDdistance)))
							    && (strcmp(map[b_list].chr,map[c_list].chr)!=0 || (strcmp(map[b_list].chr,map[c_list].chr)==0&& (abs(map[c_list].pos - map[b_list].pos) > LDdistance))))
						      {

						      }
						    else
						      {
							continue;
						      }
						  }
						//CASE_ONLY

						aMod= SNPMapInverse[a];
						bMod= SNPMapInverse[b];
						cMod= SNPMapInverse[c];

#pragma omp critical (NTESTS)
						{
						  if (!caseOnly)
						    {
						      ntests++;
						    }
						  //CASE_ONLY
						  else if((test == 15 || test == 16 || test == 17 || test == 18) && (strcmp(map[a_list].chr,map[b_list].chr)!=0 || (strcmp(map[a_list].chr,map[b_list].chr)==0 && (abs(map[b_list].pos - map[a_list].pos) > LDdistance))))
						    {
						      ntests++;
						    }
						  //CASE_ONLY
						  else if((test == 19) && (strcmp(map[a_list].chr,map[b_list].chr)!=0 || (strcmp(map[a_list].chr,map[b_list].chr)==0 && (abs(map[b_list].pos - map[a_list].pos) > LDdistance)))
							  && (strcmp(map[a_list].chr,map[c_list].chr)!=0 || (strcmp(map[a_list].chr,map[c_list].chr)==0 && (abs(map[c_list].pos - map[a_list].pos) > LDdistance)))
							  && (strcmp(map[b_list].chr,map[c_list].chr)!=0 || (strcmp(map[b_list].chr,map[c_list].chr)==0 && (abs(map[c_list].pos - map[b_list].pos) > LDdistance))))
						    {
						      ntests++;
						    }

						  if (n == 0)
						    {
						      if ((ntests % 1000000 == 0 && ntests > 0) || (ntests % 100000 == 0 && ntests > 0 && ntests < 400000))
							{
							  cout << ntests << " tests done\n";
							}
						    }
						}


						if(markercombi2 == 1)
						  {
						    updateCountsMulti(casecounts5, casecounts5M, casecounts5F,controlcounts5, controlcounts5M, controlcounts5F, aMod, bMod, nwordsSNPs, BinSNPs ,BinSNPsCCFlags, BinSNPsGenderFlags,ifChrX,caseOnly,BinSNPsQTaffFlags,test);
						  }
						else if(markercombi3 == 1)
						  {
						    updateCountsMulti3(casecounts5, casecounts5M, casecounts5F, controlcounts5, controlcounts5M, controlcounts5F, aMod, bMod, cMod, nwordsSNPs, BinSNPs, BinSNPsCCFlags, BinSNPsGenderFlags,ifChrX,caseOnly);
						  }


						//multimarker tests
						if (ifChrX == 0  || regression) // not X chromosome, or regression with sexcov
						  {
						    if (test == 1)
						      {
							if (markercombi2 == 1)
							  {
							    tstat = chiTest5(casecounts5, controlcounts5, 2, mWithSingletop, helpstat);
							  }
							else if (markercombi3 == 1)
							  {
							    tstat = chiTest5(casecounts5, controlcounts5, 3, mWithSingletop, helpstat);
							  }
						      }
						    else if (test == 2)
						      {
							if (markercombi2 == 1)
							  {
							    if(!qt)
							      {
								tstat = chiTestInter5(casecounts5, controlcounts5, 2);
							      }

							  }
							else if (markercombi3 == 1)
							  {
							    tstat = chiTestInter5(casecounts5, controlcounts5, 3);
							  }
						      }
						    else if (test == 15)
						      {
							if (markercombi2 == 1)
							  {
							    if(!qt)
							      {
								if(strcmp(map[a_list].chr,map[b_list].chr)!=0 || (strcmp(map[a_list].chr,map[b_list].chr)==0 && (abs(map[b_list].pos - map[a_list].pos) > LDdistance)))
								  {
								    tstat = caseOnly1DF(casecounts5);
								    if (plusSingle == 1)
								      {
									tstat.p=ComputePlusSingle(plusSingleFirst, plusSingleSecond, plusSingleThird, tstat.df, test, ifChrX, counts[thread_nloop][a_list].pSingle, counts[thread_nloop][b_list].pSingle, 0, tstat.p,dfSingle);
								      }
								  }
							      }
							  }
						      }
						    else if (test == 16)
						      {
							if (markercombi2 == 1)
							  {
							    if(!qt)
							      {
								if(strcmp(map[a_list].chr,map[b_list].chr)!=0 || (strcmp(map[a_list].chr,map[b_list].chr)==0 && (abs(map[b_list].pos - map[a_list].pos) > LDdistance)))
								  {
								    tstat = caseOnly4DF(casecounts5);
								    if (plusSingle == 1)
								      {
									tstat.p=ComputePlusSingle(plusSingleFirst, plusSingleSecond, plusSingleThird, tstat.df, test, ifChrX, counts[thread_nloop][a_list].pSingle, counts[thread_nloop][b_list].pSingle, 0, tstat.p,dfSingle);
								      }
								  }
							      }
							  }
						      }
						    else if (test == 19)
						      {
							if (markercombi3 == 1)
							  {
							    if(!qt)
							      {
								if((strcmp(map[a_list].chr,map[b_list].chr)!=0 || (strcmp(map[a_list].chr,map[b_list].chr)==0 && (abs(map[b_list].pos - map[a_list].pos) > LDdistance)))
								   && (strcmp(map[a_list].chr,map[c_list].chr)!=0 || (strcmp(map[a_list].chr,map[c_list].chr)==0 && (abs(map[c_list].pos - map[a_list].pos) > LDdistance)))
								   && (strcmp(map[b_list].chr,map[c_list].chr)!=0 || (strcmp(map[b_list].chr,map[c_list].chr)==0 && (abs(map[c_list].pos - map[b_list].pos) > LDdistance))))
								  {
								    tstat = caseOnly3_20DF(casecounts5);
								    if (plusSingle == 1)
								      {
									tstat.p=ComputePlusSingle(plusSingleFirst, plusSingleSecond, plusSingleThird, tstat.df, test, ifChrX, counts[thread_nloop][a_list].pSingle, counts[thread_nloop][b_list].pSingle, counts[thread_nloop][c_list].pSingle, tstat.p,dfSingle);
								      }
								  }
							      }
							  }
						      }
						    else // regression multimarker
						      {
							alt = 1;
							xType=0; if(ifChrX){xType=3;}

							//TODO: pre-test X

							if(!qt)
							  {
							    tstat.p = 0;
							    if(pretest)
							      {
								if(!ifChrX)
								  {
								    if(test==3 || test==5)
								      {
									if (markercombi2 == 1)
									  {
									    tstat = pretest1df(casecounts5, controlcounts5, 2, mWithSingletop, helpstat);
									  }
									else if (markercombi3 == 1)
									  {
									    tstat = pretest1df(casecounts5, controlcounts5, 3, mWithSingletop, helpstat);
									  }

									if(test==3)
									  {
									    tstat.p=ComputePlusSingle(1, 1, 0, tstat.df, test, ifChrX, counts[thread_nloop][a_list].pSingle, counts[thread_nloop][b_list].pSingle, 0, tstat.p,dfSingle);
									  }

								      }
								    else if(test==4 || test==6)
								      {
									if (markercombi2 == 1)
									  {
									    tstat = chiTestInter5(casecounts5, controlcounts5, 2);
									  }
									else if (markercombi3 == 1)
									  {
									    tstat = chiTestInter5(casecounts5, controlcounts5, 3);
									  }
								      }
								    else if(test==17)
								      {
									tstat = caseOnly1DF(casecounts5);
								      }
								    else if(test==18)
								      {
									tstat = caseOnly4DF(casecounts5);
								      }
								  }
								else //pre-test X
								  {
								    if(test==17)
								      {
									tstat = caseOnly1DF_X(casecounts5M, casecounts5F);
								      }
								    else if(test==18)
								      {
									tstat = caseOnly4DF_X(casecounts5M, casecounts5F);
								      }
								    else if (markercombi2 == 1)
								      {
									tstat = chiTestInterX5(casecounts5M, casecounts5F, controlcounts5M,controlcounts5F, 2, mWithSingletop, helpstat);
								      }
								    else if (markercombi3 == 1)
								      {
									tstat = chiTestInterX5(casecounts5M, casecounts5F, controlcounts5M,controlcounts5F, 3, mWithSingletop, helpstat);
								      }
								  }
							      }

							    if((tstat.p <=pHelpCut))
							      {
								passedPreTest=1;
								result1[thread_aloop] = logreg(xvec, cov, sexcov || ifChrX, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, aMod, bMod, cMod,inflationfactor, casecounts5, controlcounts5, p, newbeta, X, Xmod, Xt, A,  UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N,alt,Yhelp, xType, female, male,counts[thread_nloop][0],D, T, U, Ut, sumPP, sumPJ, sumPK,MMinv,test,thread_nloop,npplqc, PPLMap, BinSNPs, PplLocations,n,result1[thread_aloop],maxIndexCov,liabilityCut,singleMarkerTest,-1,covariancematrix,df_L1,df_L2,dosage,genoWeights);

								if(result1[thread_aloop].df>0)
								  {
								    alt = 0;

								    result2[thread_aloop] = logreg(zvec, cov, sexcov || ifChrX, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, aMod, bMod, cMod,inflationfactor, casecounts5, controlcounts5, p, newbeta, X, Xmod, Xt, A,  UNNT, VNN, S,Sinv, A0, Ainv, AinvXt, YY, Yminusp, N,alt,Yhelp, xType, female, male, counts[thread_nloop][0],D, T, U, Ut, sumPP, sumPJ, sumPK,MMinv,test,thread_nloop,npplqc, PPLMap, BinSNPs, PplLocations,n,result2[thread_aloop],maxIndexCov,liabilityCut,singleMarkerTest,-1,covariancematrix,df_L1,df_L2,dosage,genoWeights);
								  }
								else
								  {
								    result2[thread_aloop].df=0;result2[thread_aloop].sc=0;
								  }
								tLog = 2*(result1[thread_aloop].sc -result2[thread_aloop].sc);
								df = result1[thread_aloop].df -result2[thread_aloop].df;


								if(MOD)
								  {
								    cout << "df " << df << "\t" << tLog << "\t" << pValueCalc(df/2, tLog/2) << "\t" << thread_aloop << " thread \n";

								    //if(n==0){exit(1);}
								  }

								if (df >= 1)
								  {
								    tstat.p = pValueCalc(df/2, tLog/2);
								    tstat.pmod = tstat.p;
								  }
								else
								  {
								    tstat.p = 1;
								    tstat.pmod = 1;
								  }
							      }
							    else {
							      tstat.p = 1;
							      tstat.pmod = 1;
							      passedPreTest=0;
							    }
							  }
							else //qt
							  {
							    if(haplo && markercombi2==1 && dohapfile && n==0) //newhaplo
							      {
								if(dohapfile >=2)
								  {
								    sprintf(hapstring,"%s %s %s %s %s %s %s %d",map[a_list].chr,map[a_list].rs,codesA[a_list].a1,codesA[a_list].a2,map[b_list].rs,codesA[b_list].a1,codesA[b_list].a2,map[b_list].pos-map[a_list].pos);
								  }
								else if(dohapfile ==1)
								  {
								    sprintf(hapstring,"%s %s %s %s %s %s %s %d",map[a_list].chr,map[a_list].rs,codesA[a_list].a1,codesA[a_list].a2,map[b_list].rs,codesA[b_list].a1,codesA[b_list].a2,map[b_list].pos-map[a_list].pos);
								  }
							      }

							    pHelp=0;
							    if(markercombi2 && pretest)
							      {
								if(test==3 || test==5 || ifChrX)
								  {
								    pHelp=anova1df(person,nlinestfam,aMod,bMod,cMod,thread_nloop,&traitavg,npplqc,PPLMap, BinSNPs, PplLocations);

								  }
								else if(test==4 || test==6)
								  {
								    pHelp=anova4df(person,nlinestfam,aMod,bMod,cMod,thread_nloop,&traitavg,npplqc,PPLMap, BinSNPs, PplLocations);
								  }
							      }

							    if(pHelp<=pHelpCut) //survived pretest
							      {
								passedPreTest=1;
								result1[thread_aloop] = qtreg(xvec, cov, sexcov || ifChrX, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, aMod, bMod, cMod, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N,alt,Yhelp, xType, female, male,counts[thread_nloop][0], Yt,YtX, YtXAinv, YtXAinvXt,haplo,D, T, U, Ut, sumPP, sumPJ, sumPK,MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,n,result1[thread_aloop],caseOnly,0,maxIndexCov,liabilityCut,singleMarkerTest,-1,covariancematrix,df_L1,df_L2);

								if(result1[thread_aloop].df>0)
								  {
								    alt = 0;
								    result2[thread_aloop] = qtreg(zvec, cov, sexcov || ifChrX, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, aMod, bMod, cMod, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N,alt,Yhelp, xType, female, male,counts[thread_nloop][0],Yt,YtX, YtXAinv, YtXAinvXt,haplo,D, T, U, Ut, sumPP, sumPJ, sumPK,MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,n,result2[thread_aloop],caseOnly,0,maxIndexCov,liabilityCut,singleMarkerTest,-1,covariancematrix,df_L1,df_L2);
								  }
								else
								  {
								    result2[thread_aloop].df=0;result2[thread_aloop].sc=0;
								  }



								df = result2[thread_aloop].df -result1[thread_aloop].df; //!!
								Fstat = ((result2[thread_aloop].sc -result1[thread_aloop].sc)/(df))/(result1[thread_aloop].sc/result1[thread_aloop].df);

								if (df >= 1 && result2[thread_aloop].sc > 0.000001)
								  {
								    tstat.p = betai(result1[thread_aloop].df/2,df/2,result1[thread_aloop].df/(result1[thread_aloop].df+df*Fstat));
								    tstat.pmod = tstat.p;

								  }
								else
								  {
								    tstat.p = 1;
								    tstat.pmod = 1;
								  }

								if(test == 18 ) //CASEONLY 18, do again for second system
								  {
								    filled=filledCells(casecounts5,0);

								    if(filled==9)
								      {
									alt=1;
									result3[thread_aloop] = qtreg(xvec2, cov, sexcov || ifChrX, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, aMod, bMod, cMod, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N,alt,Yhelp, xType, female, male,counts[thread_nloop][0], Yt,YtX, YtXAinv, YtXAinvXt,haplo,D, T, U, Ut, sumPP, sumPJ, sumPK,MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,n,result3[thread_aloop],caseOnly,1,maxIndexCov,liabilityCut,singleMarkerTest,-1,covariancematrix,df_L1,df_L2);

									if(result3[thread_aloop].df>0)
									  {
									    alt = 0;
									    result4[thread_aloop] = qtreg(zvec2, cov, sexcov || ifChrX, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, aMod, bMod, cMod, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N,alt,Yhelp, xType, female, male,counts[thread_nloop][0],Yt,YtX, YtXAinv, YtXAinvXt,haplo,D, T, U, Ut, sumPP, sumPJ, sumPK,MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,n,result4[thread_aloop],caseOnly,1,maxIndexCov,liabilityCut,singleMarkerTest,-1,covariancematrix,df_L1,df_L2);
									  }
									else
									  {
									    result4[thread_aloop].df=0;result4[thread_aloop].sc=0;
									  }

									df = result4[thread_aloop].df -result3[thread_aloop].df; //!!
									Fstat = ((result4[thread_aloop].sc -result3[thread_aloop].sc)/(df))/(result3[thread_aloop].sc/result3[thread_aloop].df);


									if (df >= 1 && result4[thread_aloop].sc > 0.000001)
									  {
									    tstat2.p = betai(result3[thread_aloop].df/2,df/2,result3[thread_aloop].df/(result3[thread_aloop].df+df*Fstat));
									    tstat2.pmod = tstat2.p;
									  }
									else
									  {
									    tstat2.p = 1;
									    tstat2.pmod = 1;
									  }

									tstat.p = pValueCalc(2, -2*(log(tstat.p)+log(tstat2.p)) / 2);
								      }

								    tstat.pmod=tstat.p;

								  } //end test ==18
							      }	//pretest
							    else
							      {
								tstat.p = 1;
								tstat.pmod = 1;
								passedPreTest=0;
							      }
							  } 	//else qt
						      } //else regression
						  } //not X
						else // multimarker X chromosom not regression
						  {
						    if (test == 1)
						      {
							if (markercombi2 == 1)
							  {
							    tstat = chiTestX5(casecounts5M, casecounts5F, controlcounts5M, controlcounts5F, 2, mWithSingletop, helpstat);
							  }
							else if (markercombi3 == 1)
							  {
							    tstat = chiTestX5(casecounts5M, casecounts5F, controlcounts5M, controlcounts5F, 3, mWithSingletop, helpstat);
							  }
						      }
						    else if (test == 2)
						      {
							if (markercombi2 == 1)
							  {
							    tstat = chiTestInterX5(casecounts5M, casecounts5F, controlcounts5M,controlcounts5F, 2, mWithSingletop, helpstat);
							  }
							else if (markercombi3 == 1)
							  {
							    tstat = chiTestInterX5(casecounts5M, casecounts5F, controlcounts5M,controlcounts5F, 3, mWithSingletop, helpstat);
							  }
						      }
						    else if (test == 15)
						      {
							if (markercombi2 == 1)
							  {
							    if(!qt)
							      {
								if(strcmp(map[a_list].chr,map[b_list].chr)!=0 || (strcmp(map[a_list].chr,map[b_list].chr)==0 && (abs(map[b_list].pos - map[a_list].pos) > LDdistance)))
								  {
								    tstat = caseOnly1DF_X(casecounts5M, casecounts5F);
								    if (plusSingle == 1)
								      {
									tstat.p=ComputePlusSingle(plusSingleFirst, plusSingleSecond, plusSingleThird, tstat.df, test, ifChrX, counts[thread_nloop][a_list].pSingle, counts[thread_nloop][b_list].pSingle, 0, tstat.p,dfSingle);
								      }
								  }

							      }
							  }
						      }
						    else if (test == 16)
						      {
							if (markercombi2 == 1)
							  {
							    if(!qt)
							      {
								if(strcmp(map[a_list].chr,map[b_list].chr)!=0 || (strcmp(map[a_list].chr,map[b_list].chr)==0 && (abs(map[b_list].pos - map[a_list].pos) > LDdistance)))
								  {
								    tstat = caseOnly4DF_X(casecounts5M, casecounts5F);
								    if (plusSingle == 1)
								      {
									tstat.p=ComputePlusSingle(plusSingleFirst, plusSingleSecond, plusSingleThird, tstat.df, test, ifChrX, counts[thread_nloop][a_list].pSingle, counts[thread_nloop][b_list].pSingle, 0, tstat.p,dfSingle);
								      }
								  }

							      }
							  }
						      }
						    else if (test == 19)
						      {
							if (markercombi3 == 1)
							  {
							    if(!qt)
							      {
								if((strcmp(map[a_list].chr,map[b_list].chr)!=0 || (strcmp(map[a_list].chr,map[b_list].chr)==0 && (abs(map[b_list].pos - map[a_list].pos) > LDdistance)))
								   && (strcmp(map[a_list].chr,map[c_list].chr)!=0 || (strcmp(map[a_list].chr,map[c_list].chr)==0 && (abs(map[c_list].pos - map[a_list].pos) > LDdistance)))
								   && (strcmp(map[b_list].chr,map[c_list].chr)!=0 || (strcmp(map[b_list].chr,map[c_list].chr)==0 && (abs(map[c_list].pos - map[b_list].pos) > LDdistance))))
								  {
								    tstat = CaseOnly3SNPs_X(casecounts5M, casecounts5F);
								    if (plusSingle == 1)
								      {
									tstat.p=ComputePlusSingle(plusSingleFirst, plusSingleSecond, plusSingleThird, tstat.df, test, ifChrX, counts[thread_nloop][a_list].pSingle, counts[thread_nloop][b_list].pSingle, counts[thread_nloop][c_list].pSingle, tstat.p,dfSingle);
								      }
								  }
							      }
							  }
						      }
						    /*else separate call for regression X chromosome, no longer active
						      {
						      if(!qt)
						      {
						      alt = 1;
						      xType = 1;


						      result1M[thread_aloop] = logreg(xvec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, aMod, bMod, cMod,inflationfactor, casecounts5M, controlcounts5M, p, newbeta, X, Xmod, Xt, A,  UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N,alt,Yhelp,xType, female, male,counts[thread_nloop][0],D, T, U, Ut, sumPP, sumPJ, sumPK,MMinv,test,thread_nloop,npplqc, PPLMap, BinSNPs, PplLocations,n,result1M[thread_aloop],maxIndexCov,liabilityCut);

						      if(result1M[thread_aloop].df>0)
						      {
						      alt = 0;

						      //resultCopy(&(resultMulti[thread_aloop]),result2M[thread_aloop]);
						      result2M[thread_aloop] = logreg(zvec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, aMod, bMod, cMod, inflationfactor, casecounts5M, controlcounts5M, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N,alt,Yhelp,xType, female, male, counts[thread_nloop][0],D, T, U, Ut, sumPP, sumPJ, sumPK,MMinv,test,thread_nloop,npplqc, PPLMap, BinSNPs, PplLocations,n,result2M[thread_aloop],maxIndexCov,liabilityCut);
						      }
						      else
						      {
						      result2M[thread_aloop].df=0;result2M[thread_aloop].sc=0;
						      }

						      tLogMale = 2*(result1M[thread_aloop].sc -result2M[thread_aloop].sc);
						      dfMale = result1M[thread_aloop].df -result2M[thread_aloop].df;

						      alt = 1;
						      xType = 2;

						      result1F[thread_aloop] = logreg(xvec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, aMod, bMod,
						      cMod, inflationfactor, casecounts5F, controlcounts5F, p, newbeta, X, Xmod, Xt, A,  UNNT, VNN,
						      S,Sinv, A0, Ainv, AinvXt, YY, Yminusp, N,alt,Yhelp,xType, female, male,counts[thread_nloop][0],D, T, U,
						      Ut, sumPP, sumPJ, sumPK,MMinv,test,thread_nloop,npplqc, PPLMap, BinSNPs, PplLocations,n,result1F[thread_aloop],maxIndexCov,liabilityCut);

						      if(result1F[thread_aloop].df>0)
						      {
						      alt = 0;

						      //resultCopy(&(resultMulti[thread_aloop]),result2F[thread_aloop]);
						      result2F[thread_aloop] = logreg(zvec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, aMod, bMod, cMod,inflationfactor, casecounts5F, controlcounts5F, p, newbeta, X, Xmod, Xt, A,UNNT, VNN, S,Sinv, A0, Ainv, AinvXt, YY, Yminusp, N,alt,Yhelp,xType, female, male, counts[thread_nloop][0],D, T, U, Ut, sumPP, sumPJ, sumPK,MMinv,test,thread_nloop,npplqc, PPLMap, BinSNPs, PplLocations,n,result2F[thread_aloop],maxIndexCov,liabilityCut);
						      }
						      else
						      {
						      result2F[thread_aloop].df=0;result2F[thread_aloop].sc=0;
						      }

						      tLogFemale = 2*(result1F[thread_aloop].sc -result2F[thread_aloop].sc);
						      dfFemale = result1F[thread_aloop].df -result2F[thread_aloop].df;

						      df = dfMale + dfFemale;

						      if (df >= 1)
						      {
						      tstat.p = pValueCalc((dfMale+dfFemale)/2, (tLogMale+tLogFemale)/2);
						      tstat.pmod = tstat.p;
						      }
						      else
						      {
						      tstat.p = 1;
						      tstat.pmod = 1;
						      }
						      }
						      else //qt chr x
						      {
						      alt = 1;
						      xType = 4;

						      if(haplo && markercombi2==1 && dohapfile && n==0) //newhaplo
						      {
						      if(dohapfile >=2)
						      {

						      sprintf(hapstring,"%s %s %c %c %s %c %c %d",map[a_list].chr,map[a_list].rs,codesA[a_list].a1,codesA[a_list].a2,map[b_list].rs,codesA[b_list].a1,codesA[b_list].a2,map[b_list].pos-map[a_list].pos);
						      }
						      else if(dohapfile ==1)
						      {
						      sprintf(hapstring,"%s %s %c %c %s %c %c %d",map[a_list].chr,map[a_list].rs,codesA[a_list].a1,codesA[a_list].a2,map[b_list].rs,codesA[b_list].a1,codesA[b_list].a2,map[b_list].pos-map[a_list].pos);
						      }
						      }


						      result1[thread_aloop] = qtreg(xvec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, aMod, bMod, cMod, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N,alt,Yhelp, xType, female, male,counts[thread_nloop][0], Yt,YtX, YtXAinv, YtXAinvXt,haplo,D, T, U, Ut, sumPP, sumPJ, sumPK,MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs, PplLocations,n,result1[thread_aloop],caseOnly,0,maxIndexCov,liabilityCut);

						      if(result1[thread_aloop].df>0)
						      {
						      alt = 0;
						      result2[thread_aloop] = qtreg(zvec, cov, sexcov, person, nlinestfam, ncasesqc, ncontrolsqc, nrestqc, aMod, bMod, cMod, inflationfactor, casecounts5, p, newbeta, X, Xmod, Xt, A, UNNT, VNN, S, Sinv, A0, Ainv, AinvXt, YY, Yminusp, N,alt,Yhelp, xType, female, male, counts[thread_nloop][0],Yt, YtX, YtXAinv, YtXAinvXt,haplo,D, T, U, Ut, sumPP, sumPJ, sumPK,MMinv, test,thread_nloop,hapfile,dohapfile,hapstring,fptr10,skip,npplqc, PPLMap, BinSNPs,  PplLocations,n,result2[thread_aloop],caseOnly,0,maxIndexCov,liabilityCut);
						      }
						      else
						      {
						      result2[thread_aloop].df=0;result2[thread_aloop].sc=0;
						      }

						      df = result2[thread_aloop].df -result1[thread_aloop].df; //!!
						      Fstat = ((result2[thread_aloop].sc -result1[thread_aloop].sc)/(df))/(result1[thread_aloop].sc/result1[thread_aloop].df);

						      if (df >= 1 && result2[thread_aloop].sc > 0.000001)
						      {
						      tstat.p = betai(result1[thread_aloop].df/2,df/2,result1[thread_aloop].df/(result1[thread_aloop].df+df*Fstat));

						      tstat.pmod = tstat.p;
						      }
						      else
						      {
						      tstat.p = 1;
						      tstat.pmod = 1;
						      }
						      }
						      }
						    */
						  }

						// Speichern der besten printtop Chi-Werte
						if (pathwayAnalysis == 1 && pathwayTest == 6)
						  {
#pragma omp critical (PATHWAY)
						    if (adjustP(tstat.p, inflationfactor, singleMarkerTest) <= p_interratio)
						      {
							if (n==0)
							  {
							    pathway[r-1].intertop[thread_nloop]++;
							  }
							else
							  {
							    pathway[r-1].intertopmod[thread_nloop]++;
							  }
						      }
						  }


						if(writeOut && tstat.p<=pfilter && passedPreTest) //writeOut
						  {
						    out2(markerCombifileParallel[thread_aloop], pathwayImpact, haplo, test, covariancematrix, df_L1, 1,
							 map, a_list, b_list, tstat.p, result1[thread_aloop],codesA,xvec,zvec,SigmaAuxVector); //writeOut
						  }


#pragma omp critical (BESTCHI)
						if (tstat.p <= pfilter && tstat.p <= lastpMulti && passedPreTest)
						  {
						    if(combilist == 1)
						      {
							maxposMulti=ntests-1;

							if(maxposMulti>(printtop-1) || maxposMulti < 0){cout << maxposMulti << " " << printtop << " warning!!!\n";exit(1);}

							bestchi3[maxposMulti].p = tstat.p;
							bestchi3[maxposMulti].pmod = tstat.pmod;
							bestchi3[maxposMulti].Fstat = Fstat;
							bestchi3[maxposMulti].nr1 = a_list;
							bestchi3[maxposMulti].nr2 = b_list;
							bestchi3[maxposMulti].nr3 = c_list;
							bestchi3[maxposMulti].nr4 = d_list;
							bestchi3[maxposMulti].nr5 = e_list;
							bestchi3[maxposMulti].r = r;
							bestchi3[maxposMulti].X = ifChrX;

							for (ii = 0; ii < 3; ii++)
							  {
							    for (jj = 0; jj < 3; jj++)
							      {
								for (kk = 0; kk < 3; kk++)
								  {
								    for (ll = 0; ll < 1; ll++)
								      {
									for (mm = 0; mm < 1; mm++)
									  {
									    bestchi3[maxposMulti].casecounts[ii][jj][kk][ll][mm] = casecounts5[ii][jj][kk][ll][mm];
									    bestchi3[maxposMulti].casecountsMale[ii][jj][kk][ll][mm] = casecounts5M[ii][jj][kk][ll][mm];
									    bestchi3[maxposMulti].casecountsFemale[ii][jj][kk][ll][mm] = casecounts5F[ii][jj][kk][ll][mm];
									    bestchi3[maxposMulti].controlcounts[ii][jj][kk][ll][mm] = controlcounts5[ii][jj][kk][ll][mm];
									    bestchi3[maxposMulti].controlcountsMale[ii][jj][kk][ll][mm] = controlcounts5M[ii][jj][kk][ll][mm];
									    bestchi3[maxposMulti].controlcountsFemale[ii][jj][kk][ll][mm] = controlcounts5F[ii][jj][kk][ll][mm];
									  }
								      }
								  }
							      }
							  }

							if(regression)
							  {
							    PlusToPlus(&(bestchi3[maxposMulti].result1),result1[thread_aloop],haplo,maxIndexCov,printBeta,1,covariancematrix,dim_Cov1,dim_Cov2);
							    if(test==18) {PlusToPlus(&(bestchi3[maxposMulti].result3),result3[thread_aloop],haplo,maxIndexCov,printBeta,1,covariancematrix,dim_Cov1,dim_Cov2);}
							  }
							bestchi3[maxposMulti].traitavg=traitavg;

						      }
						    else
						      {
							if(maxposMulti>(printtop-1) || maxposMulti < 0){cout << maxposMulti << " " << printtop << " warning!!!\n";exit(1);}


							if (tstat.p <= lastpMulti)
							  {
							    if(lastpMulti<1.1) // list is full
							      {
								bestchi3[maxposMulti].p = tstat.p;
								bestchi3[maxposMulti].pmod = tstat.pmod;
								bestchi3[maxposMulti].Fstat = Fstat;
								bestchi3[maxposMulti].nr1 = a_list;
								bestchi3[maxposMulti].nr2 = b_list;
								bestchi3[maxposMulti].nr3 = c_list;
								bestchi3[maxposMulti].nr4 = d_list;
								bestchi3[maxposMulti].nr5 = e_list;
								bestchi3[maxposMulti].r = r;
								bestchi3[maxposMulti].X = ifChrX;



								for (ii = 0; ii < 3; ii++)
								  {
								    for (jj = 0; jj < 3; jj++)
								      {
									for (kk = 0; kk < 3; kk++)
									  {
									    for (ll = 0; ll < 1; ll++)
									      {
										for (mm = 0; mm < 1; mm++)
										  {
										    bestchi3[maxposMulti].casecounts[ii][jj][kk][ll][mm] = casecounts5[ii][jj][kk][ll][mm];
										    bestchi3[maxposMulti].casecountsMale[ii][jj][kk][ll][mm] = casecounts5M[ii][jj][kk][ll][mm];
										    bestchi3[maxposMulti].casecountsFemale[ii][jj][kk][ll][mm] = casecounts5F[ii][jj][kk][ll][mm];
										    bestchi3[maxposMulti].controlcounts[ii][jj][kk][ll][mm] = controlcounts5[ii][jj][kk][ll][mm];
										    bestchi3[maxposMulti].controlcountsMale[ii][jj][kk][ll][mm] = controlcounts5M[ii][jj][kk][ll][mm];
										    bestchi3[maxposMulti].controlcountsFemale[ii][jj][kk][ll][mm] = controlcounts5F[ii][jj][kk][ll][mm];
										  }
									      }
									  }
								      }
								  }

								if(regression)
								  {
								    PlusToPlus(&(bestchi3[maxposMulti].result1),result1[thread_aloop],haplo,maxIndexCov,printBeta,1,covariancematrix,dim_Cov1,dim_Cov2);
								    if(test==18) PlusToPlus(&(bestchi3[maxposMulti].result3),result3[thread_aloop],haplo,maxIndexCov,printBeta,1,covariancematrix,dim_Cov1,dim_Cov2);
								  }
								bestchi3[maxposMulti].traitavg=traitavg;

								//get new maxpos
								lastpMulti=-0.1;
								for(ii=0;ii<printtop;ii++)
								  {
								    if(bestchi3[ii].p>lastpMulti)
								      {
									lastpMulti=bestchi3[ii].p;
									maxposMulti=ii;
								      }
								  }
							      }
							    else if (maxposMulti < printtop -1 )
							      {
								bestchi3[maxposMulti].p = tstat.p;
								bestchi3[maxposMulti].pmod = tstat.pmod;
								bestchi3[maxposMulti].Fstat = Fstat;
								bestchi3[maxposMulti].nr1 = a_list;
								bestchi3[maxposMulti].nr2 = b_list;
								bestchi3[maxposMulti].nr3 = c_list;
								bestchi3[maxposMulti].nr4 = d_list;
								bestchi3[maxposMulti].nr5 = e_list;
								bestchi3[maxposMulti].r = r;
								bestchi3[maxposMulti].X = ifChrX;



								for (ii = 0; ii < 3; ii++)
								  {
								    for (jj = 0; jj < 3; jj++)
								      {
									for (kk = 0; kk < 3; kk++)
									  {
									    for (ll = 0; ll < 1; ll++)
									      {
										for (mm = 0; mm < 1; mm++)
										  {
										    bestchi3[maxposMulti].casecounts[ii][jj][kk][ll][mm] = casecounts5[ii][jj][kk][ll][mm];
										    bestchi3[maxposMulti].casecountsMale[ii][jj][kk][ll][mm] = casecounts5M[ii][jj][kk][ll][mm];
										    bestchi3[maxposMulti].casecountsFemale[ii][jj][kk][ll][mm] = casecounts5F[ii][jj][kk][ll][mm];
										    bestchi3[maxposMulti].controlcounts[ii][jj][kk][ll][mm] = controlcounts5[ii][jj][kk][ll][mm];
										    bestchi3[maxposMulti].controlcountsMale[ii][jj][kk][ll][mm] = controlcounts5M[ii][jj][kk][ll][mm];
										    bestchi3[maxposMulti].controlcountsFemale[ii][jj][kk][ll][mm] = controlcounts5F[ii][jj][kk][ll][mm];
										  }
									      }
									  }
								      }
								  }

								if(regression)
								  {
								    PlusToPlus(&(bestchi3[maxposMulti].result1),result1[thread_aloop],haplo,maxIndexCov,printBeta,1,covariancematrix,dim_Cov1,dim_Cov2);
								    PlusToPlus(&(bestchi3[maxposMulti].result3),result3[thread_aloop],haplo,maxIndexCov,printBeta,1,covariancematrix,dim_Cov1,dim_Cov2);
								  }
								bestchi3[maxposMulti].traitavg=traitavg;

								maxposMulti++;
							      }
							    else if (maxposMulti == printtop -1 )
							      {
								bestchi3[maxposMulti].p = tstat.p;
								bestchi3[maxposMulti].pmod = tstat.pmod;
								bestchi3[maxposMulti].Fstat = Fstat;
								bestchi3[maxposMulti].nr1 = a_list;
								bestchi3[maxposMulti].nr2 = b_list;
								bestchi3[maxposMulti].nr3 = c_list;
								bestchi3[maxposMulti].nr4 = d_list;
								bestchi3[maxposMulti].nr5 = e_list;
								bestchi3[maxposMulti].r = r;
								bestchi3[maxposMulti].X = ifChrX;

								for (ii = 0; ii < 3; ii++)
								  {
								    for (jj = 0; jj < 3; jj++)
								      {
									for (kk = 0; kk < 3; kk++)
									  {
									    for (ll = 0; ll < 1; ll++)
									      {
										for (mm = 0; mm < 1; mm++)
										  {
										    bestchi3[maxposMulti].casecounts[ii][jj][kk][ll][mm] = casecounts5[ii][jj][kk][ll][mm];
										    bestchi3[maxposMulti].casecountsMale[ii][jj][kk][ll][mm] = casecounts5M[ii][jj][kk][ll][mm];
										    bestchi3[maxposMulti].casecountsFemale[ii][jj][kk][ll][mm] = casecounts5F[ii][jj][kk][ll][mm];
										    bestchi3[maxposMulti].controlcounts[ii][jj][kk][ll][mm] = controlcounts5[ii][jj][kk][ll][mm];
										    bestchi3[maxposMulti].controlcountsMale[ii][jj][kk][ll][mm] = controlcounts5M[ii][jj][kk][ll][mm];
										    bestchi3[maxposMulti].controlcountsFemale[ii][jj][kk][ll][mm] = controlcounts5F[ii][jj][kk][ll][mm];
										  }
									      }
									  }
								      }
								  }

								if(regression)
								  {
								    PlusToPlus(&(bestchi3[maxposMulti].result1),result1[thread_aloop],haplo,maxIndexCov,printBeta,1,covariancematrix,dim_Cov1,dim_Cov2);
								    PlusToPlus(&(bestchi3[maxposMulti].result3),result3[thread_aloop],haplo,maxIndexCov,printBeta,1,covariancematrix,dim_Cov1,dim_Cov2);
								  }
								bestchi3[maxposMulti].traitavg=traitavg;

								//get new maxpos
								lastpMulti=-0.1;
								for(ii=0;ii<printtop;ii++)
								  {
								    if(bestchi3[ii].p>lastpMulti)
								      {
									lastpMulti=bestchi3[ii].p;
									maxposMulti=ii;
								      }
								  }
							      }
							  } //newValue <=
						      }
						    marker1 = 0;
						    marker2 = 0;
						    marker3 = 0;
						  }
					      }
					  }
					oldCombi = 0;
				      } // end e
				  } // end d
			      } // end c
			  } // end b

			//nur nötig bei omp:

			if (regression == 1 && (plla2==1 /*|| plln*/)) //TIM 1.0.9
			  {
			    free(p);p=NULL;free(S);S=NULL;
			    free(YY);YY=NULL;free(Yhelp);Yhelp = NULL;

			    free2Ddouble(X,nlinestfam);free2Ddouble(Xmod,nlinestfam);
			    free2Ddouble(Xt,dim1+maxIndexCov);free2Ddouble(A,dim1+maxIndexCov);
			    free2Ddouble(VNN,dim1+maxIndexCov);free2Ddouble(Sinv,dim1+maxIndexCov);
			    free2Ddouble(A0,dim1+maxIndexCov);free2Ddouble(UNNT,dim1+maxIndexCov);
			    free2Ddouble(Ainv,dim1+maxIndexCov);free2Ddouble(AinvXt,dim1+maxIndexCov);
			    free2Ddouble(Yminusp,nlinestfam);free2Ddouble(newbeta,dim1+maxIndexCov);
			    free2Ddouble(D,dim1+maxIndexCov);free2Ddouble(T,dim1+maxIndexCov);
			    free2Ddouble(U,dim1+maxIndexCov);free2Ddouble(Ut,dim1+maxIndexCov);
			    free2Ddouble(sumPP,dim1+maxIndexCov);free2Ddouble(sumPJ,dim1+maxIndexCov);
			    free2Ddouble(sumPK,dim1+maxIndexCov);free2Ddouble(MMinv,dim1+maxIndexCov);

			    if(qt)
			      {
				free2Ddouble(Yt,1);free2Ddouble(YtX,1);
				free2Ddouble(YtXAinv,1);free2Ddouble(YtXAinvXt,1);
			      }
			  } //free plla2
		      } // end a
#if PARALLELA
		  } //end pragma a-loop
#endif
		} // end q

	      if (regression == 1 && r==nlist && (plln) && (!bin)) //TIM 1.0.9
		{

		  free(p);p=NULL;free(S);S=NULL;
		  free(YY);YY=NULL;free(Yhelp);Yhelp = NULL;

		  free2Ddouble(X,nlinestfam);free2Ddouble(Xmod,nlinestfam);
		  free2Ddouble(Xt,dim1+maxIndexCov);free2Ddouble(A,dim1+maxIndexCov);
		  free2Ddouble(VNN,dim1+maxIndexCov);free2Ddouble(Sinv,dim1+maxIndexCov);
		  free2Ddouble(A0,dim1+maxIndexCov);free2Ddouble(UNNT,dim1+maxIndexCov);
		  free2Ddouble(Ainv,dim1+maxIndexCov);free2Ddouble(AinvXt,dim1+maxIndexCov);
		  free2Ddouble(Yminusp,nlinestfam);free2Ddouble(newbeta,dim1+maxIndexCov);
		  free2Ddouble(D,dim1+maxIndexCov);free2Ddouble(T,dim1+maxIndexCov);
		  free2Ddouble(U,dim1+maxIndexCov);free2Ddouble(Ut,dim1+maxIndexCov);
		  free2Ddouble(sumPP,dim1+maxIndexCov);free2Ddouble(sumPJ,dim1+maxIndexCov);
		  free2Ddouble(sumPK,dim1+maxIndexCov);free2Ddouble(MMinv,dim1+maxIndexCov);

		  if(qt)
		    {
		      free2Ddouble(Yt,1);free2Ddouble(YtX,1);
		      free2Ddouble(YtXAinv,1);free2Ddouble(YtXAinvXt,1);
		    }
		} //free plla2

	      if (pathwayAnalysis == 1 && pathwayTest == 6)
		{
		  if (n==0)
		    {
		      pathway[r-1].ratio0 = interactionratio(pathway[r-1], n, thread_nloop);
		      Ttable[r-1][n]=pathway[r-1].ratio0;
		    }
		  else
		    {
		      pathway[r-1].score[thread_nloop] = interactionratio(pathway[r-1], n, thread_nloop);
		      Ttable[r-1][n]=pathway[r-1].score[thread_nloop];
		    }
		}

	    } // end r

#pragma omp critical(QSORTBESTCHI)
	  if(multimarker){qsortbestchi(&bestchi3,0,printtop - 1,&x1,&y1);}
	  //			#pragma omp end critical

	  if (n == 0) // keine Simulation
	    {
	      if (markercombi2 == 1)
		{
		  cout << "\n\nNumber of effective 2-marker-tests: " << ntests << "\n";
		  logfile << "\n\nNumber of effective 2-marker-tests: " << ntests << "\n";
		  if (ntests == 0)

		    {
		      bestMarkerCombi2.open(markerCombi2file.c_str(), ios::out);
		      errorfile << "Number of effective 2-marker-tests: " << ntests << "\n";
		      logfile << "Number of effective 2-marker-tests: " << ntests << "\n";
		      errorfile.close();logfile.close();bestMarkerCombi2.close();exit(0);
		    }
		}
	      else if (markercombi3 == 1)
		{
		  cout << "\n\nNumber of effective 3-marker-tests: " << ntests << "\n";
		  logfile << "\n\nNumber of effective 3-marker-tests: " << ntests << "\n";

		  if (ntests == 0)
		    {
		      bestMarkerCombi3.open(markerCombi2file.c_str(), ios::out);
		      errorfile << "Number of effective 3-marker-tests: " << ntests << "\n";
		      logfile << "Number of effective 3-marker-tests: " << ntests << "\n";
		      errorfile.close();logfile.close();bestMarkerCombi3.close();exit(0);
		    }
		}


	      // Singlemarkertest
	      //PWT als sonderfall
	      if (pathwayAnalysis == 1)
		{
		  for (i = 0; i < printtop; i++)
		    {
		      toplist[i].p = pathway[i].score[thread_nloop];
		      toplist[i].nr1 = i;
		      toplist[i].nr2 = -1;
		      toplist[i].nr3 = -1;
		      toplist[i].r = -1;
		    }

		}
	      else if ((markercombi2 == 0 && markercombi3 == 0) || MCWithSM == 1)
		{
		  for (i = 0; i < singletop; i++)
		    {
		      toplist[i].p = bestsinglemarker[i].p;
		      toplist[i].nr1 = bestsinglemarker[i].nr;
		      toplist[i].nr2 = -1;
		      toplist[i].nr3 = -1;
		      toplist[i].r = -1;
		      if(i==printtop-1){break;}
		    }


		  for (i = singletop; i < printtop; i++)
		    {
		      toplist[i].p = 1.01;
		      toplist[i].nr1 = -1;
		      toplist[i].nr2 = -1;
		      toplist[i].nr3 = -1;
		      toplist[i].r = -1;
		    }
		}
	      else
		{
		  for (i = 0; i < printtop; i++)
		    {
		      toplist[i].p = 1.01;
		      toplist[i].nr1 = -1;
		      toplist[i].nr2 = -1;
		      toplist[i].nr3 = -1;
		      toplist[i].r = -1;
		    }
		}


	      // Chitest für 2 Marker
	      if (markercombi2 == 1)
		{
		  for (i = 0; i < printtop; i++)
		    {
		      newValueTop = bestchi3[i].p;
		      if (newValueTop <= toplist[printtop - 1].p)
			{
			  toplist[printtop - 1].p = newValueTop;
			  toplist[printtop - 1].nr1 = bestchi3[i].nr1;
			  toplist[printtop - 1].nr2 = bestchi3[i].nr2;
			  toplist[printtop - 1].nr3 = -1;
			  toplist[printtop - 1].r = bestchi3[i].r;

			  insertionTop(toplist, printtop - 1);
			}
		    }
		}

	      // Chitest für 3 Marker
	      if (markercombi3 == 1)
		{
		  for (i = 0; i < printtop; i++)
		    {
		      newValueTop = bestchi3[i].p;
		      if (newValueTop <= toplist[printtop - 1].p)
			{
			  toplist[printtop - 1].p = newValueTop;

			  toplist[printtop - 1].nr1 = bestchi3[i].nr1;
			  toplist[printtop - 1].nr2 = bestchi3[i].nr2;
			  toplist[printtop - 1].nr3 = bestchi3[i].nr3;
			  toplist[printtop - 1].r = bestchi3[i].r;

			  insertionTop(toplist, printtop - 1);
			}
		    }
		}

#if RARE
	      if(bin)
		{
		  if((qt==1 && COLLflag==0 && CMATflag==0)||qt==0)
		    {
		      if(vb==0){
			bins2testt=binRARE(cov, sexcov, person, nlinestped, nlinestfam, map, counts[thread_nloop], optimalrare, raref, nsim, n, nwindows, nwindowssinglebin, nRareSNPs, windowPositions, nRareLimits, rareLimits, MAF_Level_VT_FISHER, MAF_Level_VT_CMAT, MAF_Level_VT_COLL, MAF_Level_VT_REGRESSION, MAF_Level_VT_FRACREG, MAF_Level_VT_COLLREG, pCOLLvec, pCMATvec, pFISHERvec, pFISHERvecChi, pREGRESSIONvec, pFRACREGvec, pCOLLREGvec, limitsStatFISHER,limitsStatREGRESSION, limitsStatFRACREG,limitsStatCOLLREG, limitsStatCOLL, limitsStatCMAT, thread_nloop, rarefile, intervalfile, rarefFISHER, rarefREGRESSION, rarefFRACREG,rarefCOLLREG, rarefCOLL, OR_COLL_vec, OR_COLL_f_vec, rarefCMAT, OR_CMAT_vec, inflationfactor, errorfile, logfile, wilsonpretest, rareregpretest, rarepretest, rarepretestlimit, FISHERnotconverged,FISHERpretestpassed, REGRESSIONpretestpassed, FRACREGpretestpassed, COLLREGpretestpassed, CMATpretestpassed, COLLpretestpassed, window, FISHER_CI, CMAT_CI, COLL_CI, REGRESSION_CI, FRACREG_CI, COLLREG_CI, COLLREG_beta, COLLREG_se, FRACREG_beta, FRACREG_se, FISHERflag, REGRESSIONflag, FRACREGflag,COLLREGflag, COLLflag, CMATflag, FISHERcount, REGRESSIONcount, FRACREGcount,COLLREGcount, COLLcount, CMATcount, FISHERstats, REGRESSIONstats, FRACREGstats,COLLREGstats, COLLstats, CMATstats, FISHERpermstats, REGRESSIONpermstats, FRACREGpermstats,COLLREGpermstats, COLLpermstats, CMATpermstats, completedFISHER, completedREGRESSION, completedFRACREG,completedCOLLREG, completedCOLL, completedCMAT, &currentn, BinSNPs, nwordsSNPs, BinSNPsCCFlags,BinSNPsCovCatFlags, BinSNPsGenderFlags, SNPMapInverse, SNPMap, PPLMap, PplLocations, npplqc, nsnpqc, fisherCorrection, MatchingRARE, nMWindowsRARE, rare_stratify, fulltests, teststat, featurecol, casecounts5, controlcounts5, p, newbeta, X, Xmod, Xt, A,  UNNT, VNN,  S, Sinv, A0, Ainv, AinvXt, YY, Yt, YtX, YtXAinv, YtXAinvXt, Yminusp, N, alt, Yhelp, xType, female, male, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test, npplqc, resultSingle[thread_nloop], maxIndexCov, ncov, qt, covariancematrix, xsinglevec, zsinglevec, ncasesqc, ncontrolsqc, nrestqc, singleMarkerTest,haplo,maxthreads, liabilityCut, binamin, binamax, binadjust, comparemode,ncasesqcMale, ncasesqcFemale,ncontrolsqcMale, ncontrolsqcFemale, &ix, &iy, &iz, weights,SetIDfile,intervaleditor,setid,verbose,adaptive, nCovCathegories, CovCathegories, collinter, doublewindowcoord, vb, nvbstartvt, vbstartvt, nvbend, vbstart, vbend, vbmaxpermstat, dim1, nvblevel, nchunks, chunkpos, chunklen, nCarriers, ndummyatlevel, dummypos, dummylevel, BinCarriers, genoWeights, dosage);
		      }
		      else if(vb==1 && vb_binwise_corr){
			if(optimalrare==1){
			  die("VB+VT is not compatible with VB_BINWISE_CORR!");
			}
			else if(vb_binwise_corr==1 && vb_print_perm_stat==1){
			  die("PRINT_PERM_STAT is not compatible with VB_BINWISE_CORR!");
			}
			else if(vb_binwise_corr==1 && verbose<2){
			  die("VB_BINWISE_CORR works obly with VERBOSE 2!");
			}

			calcVB(optimalrare, ncasesqc, ncontrolsqc, map, BinSNPs, BinSNPsCCFlags, nwordsSNPs, SNPMapInverse, n, nvbstartvt, vbstartvt, nvbend, vbend, vbmaxpermstat, nvblevel, nchunks, chunkpos, chunklen, nwindows, nCarriers, ndummyatlevel, dummypos, dummylevel, BinCarriers, dummyend, ndummyends, nchunkcluster, dummycluster, vb_binwise_corr, vbbinwisestat, vbbinwisecount);
		      }
		      if (regression == 1 /*&& r==nlist*/ && (plln)) //TIM 1.0.9 + 507
			{

			  free(p);p=NULL;free(S);S=NULL;
			  free(YY);YY=NULL;free(Yhelp);Yhelp = NULL;

			  free2Ddouble(X,nlinestfam);free2Ddouble(Xmod,nlinestfam);
			  free2Ddouble(Xt,dim1+maxIndexCov);free2Ddouble(A,dim1+maxIndexCov);
			  free2Ddouble(VNN,dim1+maxIndexCov);free2Ddouble(Sinv,dim1+maxIndexCov);
			  free2Ddouble(A0,dim1+maxIndexCov);free2Ddouble(UNNT,dim1+maxIndexCov);
			  free2Ddouble(Ainv,dim1+maxIndexCov);free2Ddouble(AinvXt,dim1+maxIndexCov);
			  free2Ddouble(Yminusp,nlinestfam);free2Ddouble(newbeta,dim1+maxIndexCov);
			  free2Ddouble(D,dim1+maxIndexCov);free2Ddouble(T,dim1+maxIndexCov);
			  free2Ddouble(U,dim1+maxIndexCov);free2Ddouble(Ut,dim1+maxIndexCov);
			  free2Ddouble(sumPP,dim1+maxIndexCov);free2Ddouble(sumPJ,dim1+maxIndexCov);
			  free2Ddouble(sumPK,dim1+maxIndexCov);free2Ddouble(MMinv,dim1+maxIndexCov);

			  if(qt)
			    {
			      free2Ddouble(Yt,1);free2Ddouble(YtX,1);
			      free2Ddouble(YtXAinv,1);free2Ddouble(YtXAinvXt,1);
			    }
			} //free plla2

		    }
		  else
		    {
		      cout << "COLL/CMAT can not be used with quantitative traits!"<<endl;
		      logfile << "COLL/CMAT can not be used with quantitative traits!"<<endl;
		      errorfile << "COLL/CMAT can not be used with quantitative traits!"<<endl;
		      exit(1);
		    }
		}
#endif
	      if (nsim && !bin) {
		sstm << "\nMonte-Carlo on P-values: " << nsim << " simulations, IBS-stratifying mode " << stratify;
		if (singleMarkerTest==1) sstm << " (" << (group_test==2 ? "MCAT(2)" : group_test==1 ? "MCAT(1)" : "CAT") << "-test).";
		logg(sstm);
	      }

	      barrier = false;
	    }  // if (n==0)

	  if (nsim > 0) // Simulation
	    {
	      if (!thread_nloop) {
		cout << "\rn = " << setw(ceil(log10(nsim))) << n << " of " << nsim << "\t\t\t\t\0" << flush;
	      }
	      if (n > 0)
		{

#if RARE
		  if(bin){
		    if(vb==0){
		      bins2testt=binRARE(cov, sexcov, person, nlinestped, nlinestfam, map, counts[thread_nloop], optimalrare, raref, nsim, n, nwindows, nwindowssinglebin, nRareSNPs, windowPositions, nRareLimits, rareLimits, MAF_Level_VT_FISHER, MAF_Level_VT_CMAT, MAF_Level_VT_COLL, MAF_Level_VT_REGRESSION, MAF_Level_VT_FRACREG, MAF_Level_VT_COLLREG, pCOLLvec, pCMATvec, pFISHERvec, pFISHERvecChi, pREGRESSIONvec, pFRACREGvec, pCOLLREGvec, limitsStatFISHER, limitsStatREGRESSION, limitsStatFRACREG, limitsStatCOLLREG, limitsStatCOLL, limitsStatCMAT, thread_nloop, rarefile, intervalfile, rarefFISHER, rarefREGRESSION, rarefFRACREG,rarefCOLLREG, rarefCOLL, OR_COLL_vec, OR_COLL_f_vec, rarefCMAT, OR_CMAT_vec, inflationfactor, errorfile, logfile, wilsonpretest,rareregpretest,rarepretest, rarepretestlimit, FISHERnotconverged,FISHERpretestpassed, REGRESSIONpretestpassed, FRACREGpretestpassed, COLLREGpretestpassed, CMATpretestpassed, COLLpretestpassed, window, FISHER_CI, CMAT_CI, COLL_CI, REGRESSION_CI, FRACREG_CI, COLLREG_CI, COLLREG_beta, COLLREG_se, FRACREG_beta, FRACREG_se, FISHERflag, REGRESSIONflag, FRACREGflag,COLLREGflag, COLLflag, CMATflag, FISHERcount, REGRESSIONcount, FRACREGcount,COLLREGcount, COLLcount, CMATcount, FISHERstats, REGRESSIONstats, FRACREGstats,COLLREGstats, COLLstats, CMATstats, FISHERpermstats, REGRESSIONpermstats, FRACREGpermstats,COLLREGpermstats, COLLpermstats, CMATpermstats, completedFISHER, completedREGRESSION, completedFRACREG,completedCOLLREG, completedCOLL, completedCMAT, &currentn, BinSNPs, nwordsSNPs, BinSNPsCCFlags, BinSNPsCovCatFlags, BinSNPsGenderFlags, SNPMapInverse, SNPMap, PPLMap, PplLocations, npplqc, nsnpqc, fisherCorrection, MatchingRARE, nMWindowsRARE, rare_stratify, fulltests, teststat, featurecol, casecounts5, controlcounts5, p, newbeta, X, Xmod, Xt, A,  UNNT, VNN,  S, Sinv, A0, Ainv, AinvXt, YY, Yt, YtX, YtXAinv, YtXAinvXt, Yminusp, N, alt, Yhelp, xType, female, male, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, test, npplqc, resultSingle[thread_nloop], maxIndexCov, ncov, qt, covariancematrix, xsinglevec, zsinglevec, ncasesqc, ncontrolsqc, nrestqc, singleMarkerTest,haplo,maxthreads, liabilityCut, binamin, binamax, binadjust, comparemode,ncasesqcMale, ncasesqcFemale,ncontrolsqcMale, ncontrolsqcFemale, &ix, &iy, &iz, weights,SetIDfile,intervaleditor,setid,verbose,adaptive, nCovCathegories, CovCathegories, collinter, doublewindowcoord, vb, nvbstartvt, vbstartvt, nvbend, vbstart, vbend, vbmaxpermstat, dim1,  nvblevel, nchunks, chunkpos, chunklen, nCarriers, ndummyatlevel, dummypos, dummylevel, BinCarriers, genoWeights, dosage);
		    }
		    else if(vb==1){
		      calcVB(optimalrare, ncasesqc, ncontrolsqc, map, BinSNPs, BinSNPsCCFlags, nwordsSNPs, SNPMapInverse, n, nvbstartvt, vbstartvt, nvbend, vbend, vbmaxpermstat, nvblevel, nchunks, chunkpos, chunklen, nwindows, nCarriers, ndummyatlevel, dummypos, dummylevel, BinCarriers, dummyend, ndummyends, nchunkcluster, dummycluster, vb_binwise_corr, vbbinwisestat, vbbinwisecount);
		      //		      if(n==nsim){
		      // for(int l=0; l<nwindows; l++){
		      //   for(int k=0; k<nvbstartvt[l]; k++){
		      //   }
		      // }
		      //			for(int k=0; k<nsim; k++){
		      //		  cout <<vbmaxpermstat[k]<<endl;
		      //			}
		      //		      }
		    }

		    if (regression == 1 /*&& r==nlist*/ && plln) //TIM 1.0.9
		      {
			free(p);p=NULL;free(S);S=NULL;
			free(YY);YY=NULL;free(Yhelp);Yhelp = NULL;

			free2Ddouble(X,nlinestfam);free2Ddouble(Xmod,nlinestfam);
			free2Ddouble(Xt,dim1+maxIndexCov);free2Ddouble(A,dim1+maxIndexCov);
			free2Ddouble(VNN,dim1+maxIndexCov);free2Ddouble(Sinv,dim1+maxIndexCov);
			free2Ddouble(A0,dim1+maxIndexCov);free2Ddouble(UNNT,dim1+maxIndexCov);
			free2Ddouble(Ainv,dim1+maxIndexCov);free2Ddouble(AinvXt,dim1+maxIndexCov);
			free2Ddouble(Yminusp,nlinestfam);free2Ddouble(newbeta,dim1+maxIndexCov);
			free2Ddouble(D,dim1+maxIndexCov);free2Ddouble(T,dim1+maxIndexCov);
			free2Ddouble(U,dim1+maxIndexCov);free2Ddouble(Ut,dim1+maxIndexCov);
			free2Ddouble(sumPP,dim1+maxIndexCov);free2Ddouble(sumPJ,dim1+maxIndexCov);
			free2Ddouble(sumPK,dim1+maxIndexCov);free2Ddouble(MMinv,dim1+maxIndexCov);

			if(qt)
			  {
			    free2Ddouble(Yt,1);free2Ddouble(YtX,1);
			    free2Ddouble(YtXAinv,1);free2Ddouble(YtXAinvXt,1);
			  }
		      } //free plla2

		  }
#endif

		  bestofSim[n]= 1.01;
		  //bestofSim2[n]= 1.01;
		  //bestofSim3[n]= 1.01;
		  // Topliste
		  // Singlemarkertest

		  if (pathwayAnalysis == 1)
		    {

		    }
		  else if ((markercombi2 == 0 && markercombi3 == 0) || MCWithSM == 1)
		    {
		      bestofSim[n]= bestsinglemarker[0].p;
		    }


		  // Chitest für 2 oder 3 Marker
		  if (markercombi2 == 1 || markercombi3 == 1)
		    {
		      if(!MCWithSM || bestchi3[0].p <= bestofSim[n])
			{
			  bestofSim[n] = bestchi3[0].p;
			}
		    }
		} // end n > 0

	    } // end nsim>0

	  if (n == 0 )
	    {
	      // Ausgabe auf dem Bildschirm und Schreiben in bestMarkerCombi2 bzw. bestMarkerCombi3
	      if (markercombi2 == 1)
		{

		  bestMarkerCombi2.open(markerCombi2file.c_str(), ios::out);

		  bestMarkerCombi2Details.open(markerCombi2Detailsfile.c_str(), ios::out);
		  // 10 besten P-Werte werden ausgeben und  besten in eine Datei
		  cout << "\n2-marker: top p-values multi-marker analysis:\n";
		  logfile << "\n2-marker: top p-values multi-marker analysis:\n";

		  if (pathwayImpact == 1)
		    {
		      if (test < 3)
			{
			  cout << "No Chr_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_2 rs_No_2 Pos_No_2 p-Single-marker_2 p-value p-value_corr Pathway\n";
			  logfile << "No Chr_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_2 rs_No_2 Pos_No_2 p-Single-marker_2 p-value p-value_corr Pathway\n";
			  bestMarkerCombi2 <<"No\tChr_1\trs_No_1\tPos_No_1\tp-Single-marker_1\tChr_2\trs_No_2\tPos_Pos_2\tp-Single-marker_2\tp-value\tp-value_corr\tPathway\n";
			}
		      else
			{
			  cout << "No Chr_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_2 rs_No_2 Pos_No_2 p-Single-marker_2 p-value p-value_corr Pathway\n";
			  logfile << "No Chr_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_2 rs_No_2 Pos_No_2 p-Single-marker_2 p-value p-value_corr Pathway\n";
			  bestMarkerCombi2 <<"No\tChr_1\trs_No_1\tPos_No_1\tp-Single-marker_1\tChr_2\trs_No_2\tPos_Pos_2\tp-Single-marker_2\tp-value\tp-value_corr\tPathway\tSNP1_Allele_A\tSNP1_Allele_B\tSNP2_Allele_A\tSNP2_Allele_B\tbeta_x1\tse_x1\tbeta_x1D\tse_x1D\tbeta_x2\tse_x2\tbeta_x2D\tse_x2D\tbeta_x1x2\tse_x1x2\tbeta_x1x2D\tse_x1x2D\tbeta_x1Dx2\tse_x1Dx2\tbeta_x1Dx2D\tse_x1Dx2D\n";
			}
		    }
		  else
		    {

		      if (test < 3 || test==15 || test==16)
			{
			  cout << "No Chr_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_2 rs_No_2 Pos_No_2 p-Single-marker_2 p-value p-value_corr\n";
			  logfile << "No Chr_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_2 rs_No_2 Pos_No_2 p-Single-marker_2 p-value p-value_corr\n";

			  bestMarkerCombi2 <<
			    "No\tChr_1\trs_No_1\tPos_No_1\tp-Single-marker_1\tChr_2\trs_No_2\tPos_No_2\tp-Single-marker_2\tp-value\tp-value_corr\n";


			}
		      else if(!haplo)

			{
			  cout << "No Chr_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_2 rs_No_2 Pos_No_2 p-Single-marker_2 p-value p-value_corr\n";
			  logfile << "No Chr_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_2 rs_No_2 Pos_No_2 p-Single-marker_2 p-value p-value_corr\n";

			  bestMarkerCombi2 <<
			    "No\tChr_1\trs_No_1\tPos_No_1\tp-Single-marker_1\tChr_2\trs_No_2\tPos_No_2\tp-Single-marker_2\tp-value\tp-value_corr\tSNP1_Allele_A\tSNP1_Allele_B\tSNP2_Allele_A\tSNP2_Allele_B\tbeta_x1\tse_x1\tbeta_x1D\tse_x1D\tbeta_x2\tse_x2\tbeta_x2D\tse_x2D\tbeta_x1x2\tse_x1x2\tbeta_x1x2D\tse_x1x2D\tbeta_x1Dx2\tse_x1Dx2\tbeta_x1Dx2D\tse_x1Dx2D";

			  if(covariancematrix)
			    {
			      for(int m=0; m<df_L1; m++)
				//for(int m=0; m<df_L1-maxIndexCov; m++)
				{
				  for(int s=m; s<df_L1; s++)
				    //for(int s=m; s<df_L1-maxIndexCov; s++)
				    {
				      bestMarkerCombi2 << "\tsigma1[" << m <<"]" << "[" << s << "]";
				    }
				}

			    }
			  bestMarkerCombi2 << "\n";

			}
		      else if (haplo)

			{
			  cout << "No Chr_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_2 rs_No_2 Pos_No_2 p-Single-marker_2 p-value p-value_corr\n";
			  logfile << "No Chr_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_2 rs_No_2 Pos_No_2 p-Single-marker_2 p-value p-value_corr\n";

			  bestMarkerCombi2 <<
			    "No\tChr_1\trs_No_1\tPos_No_1\tp-Single-marker_1\tChr_2\trs_No_2\tPos_No_2\tp-Single-marker_2\tp-value\tp-value_corr\tSNP1_Allele_A\tSNP1_Allele_B\tSNP2_Allele_A\tSNP2_Allele_B\tbeta_h1\tse_h1\tbeta_h2\tse_h2\tbeta_h3\tse_h3\tbeta_h4\tse_h4\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\n";

			}
		    }

		  withmarginal = xvec[1]+xvec[2]+xvec[3]+xvec[4]+xvec[9]+xvec[10]-zvec[1]-zvec[2]-zvec[3]-zvec[4]-zvec[9]-zvec[10];

		  if(mWithSingletop==1)
		    {
		      withmarginal=withmarginal-(xvec[3]-zvec[3])-(xvec[4]-zvec[4]);
		    }

		  if(haplo){withmarginal=0;}

		  for (k = 0; k < printtop; k++)
		    {

		      //Korrektur
		      if (mWithSingletop == 0 || combilist == 1 || withmarginal==0 || haplo || caseOnly && !(plusSingle && mWithSingletop>0) )
			{
			  bestchi3[k].pmod = Sidak(bestchi3[k].p,ntests);

			}
		      /*else if (mWithSingletop == 0 || test==2 || (test>2 && withmarginal==0) || haplo)
			{
			if (pathwayImpact == 0 && !haplo)
			{
			//bestchi3[k].pmod = bestchi3[k].p * ntestsUser;
			bestchi3[k].pmod = Sidak(bestchi3[k].p,ntestsUser);

			}
			else
			{
			//bestchi3[k].pmod = bestchi3[k].p * ntests;
			bestchi3[k].pmod = Sidak(bestchi3[k].p,ntests);
			}
			}
			else if (test == 1 && pathwayImpact == 0)
			{
			//bestchi3[k].pmod = bestchi3[k].pmod * ntestsUser;
			bestchi3[k].pmod = Sidak(bestchi3[k].pmod,ntestsUser);
			}*/
		      else
			{
			  //bestchi3[k].pmod = bestchi3[k].p * ntestsPlus;
			  bestchi3[k].pmod = Sidak(bestchi3[k].p,ntestsPlus);
			}

		      if (bestchi3[k].nr1 != -1 && bestchi3[k].nr2 != -1)
			{
			  //cout << bestchi3[k].nr1 << " " << bestchi3[k].nr2 << "\n";

			  if (pathwayImpact == 1)
			    {
			      if (k < 10)
				{
				  cout << k + 1 << " " << map[bestchi3[k].nr1].chr << " " << map[bestchi3[k].nr1].rs << " "
				       << map[bestchi3[k].nr1].pos  << " " << map[bestchi3[k].nr1].p << " " << map[bestchi3[k].nr2].chr << " "
				       << map[bestchi3[k].nr2].rs << " " << map[bestchi3[k].nr2].pos << " "
				       << map[bestchi3[k].nr2].p << " " << bestchi3[k].p << " " << bestchi3[k].pmod << " "
				       << pathway[bestchi3[k].r-1].name << "\n";

				  logfile << k + 1 << " " << map[bestchi3[k].nr1].chr << " " << map[bestchi3[k].nr1].rs << " "
					  << map[bestchi3[k].nr1].pos  << " " << map[bestchi3[k].nr1].p << " "
					  << map[bestchi3[k].nr2].chr << " " << map[bestchi3[k].nr2].rs << " " << map[bestchi3[k].nr2].pos << " "
					  << map[bestchi3[k].nr2].p << " " << bestchi3[k].p << " " << bestchi3[k].pmod << " "
					  << pathway[bestchi3[k].r-1].name << "\n";
				}
			      bestMarkerCombi2 << k + 1 << "\t" << map[bestchi3[k].nr1].chr << "\t" << map[bestchi3[k].nr1].rs << "\t"
					       << map[bestchi3[k].nr1].pos  << "\t" << map[bestchi3[k].nr1].p
					       << "\t" << map[bestchi3[k].nr2].chr << "\t" << map[bestchi3[k].nr2].rs << "\t" << map[bestchi3[k].nr2].pos
					       << "\t" << map[bestchi3[k].nr2].p << "\t" << bestchi3[k].p << "\t" << bestchi3[k].pmod << "\t"
					       << pathway[bestchi3[k].r-1].name << "\n";
			    }
			  else
			    {
			      if (k < 10)
				{
				  cout << k + 1 << " " << map[bestchi3[k].nr1].chr << " " << map[bestchi3[k].nr1].rs << " "
				       << map[bestchi3[k].nr1].pos  << " " << map[bestchi3[k].nr1].p
				       << " " << map[bestchi3[k].nr2].chr << " " << map[bestchi3[k].nr2].rs << " " << map[bestchi3[k].nr2].pos << " "
				       << map[bestchi3[k].nr2].p << " " << bestchi3[k].p << " " << bestchi3[k].pmod << "\n";

				  logfile << k + 1 << " " << map[bestchi3[k].nr1].chr << " " << map[bestchi3[k].nr1].rs << " "
					  << map[bestchi3[k].nr1].pos  << " " << map[bestchi3[k].nr1].p << " " << map[bestchi3[k].nr2].chr << " "
					  << map[bestchi3[k].nr2].rs << " " << map[bestchi3[k].nr2].pos << " " << map[bestchi3[k].nr2].p
					  << " " << bestchi3[k].p << " " << bestchi3[k].pmod  << "\n";
				}


			      bestMarkerCombi2 << k + 1 << "\t" << map[bestchi3[k].nr1].chr << "\t" << map[bestchi3[k].nr1].rs << "\t"
					       << map[bestchi3[k].nr1].pos  << "\t" << map[bestchi3[k].nr1].p << "\t" << map[bestchi3[k].nr2].chr << "\t"
					       << map[bestchi3[k].nr2].rs << "\t" << map[bestchi3[k].nr2].pos << "\t" << map[bestchi3[k].nr2].p
					       << "\t" << bestchi3[k].p << "\t" << bestchi3[k].pmod;
			      if (test < 3 || test==15 || test==16)
				{
				  bestMarkerCombi2 << "\n";
				}
			      else
				{
				  bestMarkerCombi2 << "\t" << codesA[bestchi3[k].nr1].a1 << "\t" << codesA[bestchi3[k].nr1].a2 << "\t" << codesA[bestchi3[k].nr2].a1 << "\t" << codesA[bestchi3[k].nr2].a2;

				  if (test != 17 && test != 18)
				    {
				      for (i =1; i<9; i++)
					{
					  if (xvec[i] != zvec[i])
					    {
					      bestMarkerCombi2 << "\t" << bestchi3[k].result1.b[i] << "\t" << bestchi3[k].result1.betaNew_se[i];
					    }
					  else
					    {
					      bestMarkerCombi2 << "\t-\t-";
					    }
					}

				      if(covariancematrix)
					{
					  for(int m=0; m < ((df_L1)*(df_L1+1))/2; m++)
					    {
					      bestMarkerCombi2 << "\t" << bestchi3[k].result1.sigma1[SigmaAuxVector[m]];
					    }
					}

				      bestMarkerCombi2 << "\n";
				    }
				  else
				    {
				      if (test == 17)
					{
					  bestMarkerCombi2 << "\t-\t-\t-\t-\t-\t-\t-\t-\t" << bestchi3[k].result1.b[3] << "\t" << bestchi3[k].result1.betaNew_se[3] << "\t-\t-\t-\t-\t-\t-\n";
					}
				      if (test == 18)
					{
					  bestMarkerCombi2 << "\t-\t-\t-\t-\t-\t-\t-\t-\t" << bestchi3[k].result1.b[3] << "\t" << bestchi3[k].result1.betaNew_se[3] << "\t" << bestchi3[k].result1.b[4] << "\t" << bestchi3[k].result1.betaNew_se[4] << "\t" << bestchi3[k].result3.b[3] << "\t" << bestchi3[k].result3.betaNew_se[3] << "\t" << bestchi3[k].result3.b[4] << "\t" << bestchi3[k].result3.betaNew_se[4] << "\n";
					}
				    }

				}
			    }

			  // Details
			  if (bestchi3[k].X == 0)
			    {
			      if (!qt)
				{
				  // Ausgabe logistische Regression
				  bestMarkerCombi2Details << k+1 << ": SNP1: " << map[bestchi3[k].nr1].rs << "\tSNP2: " << map[bestchi3[k].nr2].rs
							  << "\n\n";
				  bestMarkerCombi2Details << "cases:\t\t\t\tcontrols:\t\t\n";
				  bestMarkerCombi2Details << bestchi3[k].casecounts[0][0][0][0][0] << "\t" << bestchi3[k].casecounts[0][1][0][0][0]
							  << "\t" << bestchi3[k].casecounts[0][2][0][0][0] << "\t\t\t" << bestchi3[k].controlcounts[0][0][0][0][0] << "\t"

							  << bestchi3[k].controlcounts[0][1][0][0][0] << "\t" << bestchi3[k].controlcounts[0][2][0][0][0] << "\n"
							  << bestchi3[k].casecounts[1][0][0][0][0] << "\t" << bestchi3[k].casecounts[1][1][0][0][0] << "\t"
							  << bestchi3[k].casecounts[1][2][0][0][0]  << "\t\t\t" << bestchi3[k].controlcounts[1][0][0][0][0] << "\t"
							  << bestchi3[k].controlcounts[1][1][0][0][0] << "\t" << bestchi3[k].controlcounts[1][2][0][0][0] << "\n"
							  << bestchi3[k].casecounts[2][0][0][0][0] << "\t" << bestchi3[k].casecounts[2][1][0][0][0] << "\t"
							  << bestchi3[k].casecounts[2][2][0][0][0]  << "\t\t\t" << bestchi3[k].controlcounts[2][0][0][0][0] << "\t"
							  << bestchi3[k].controlcounts[2][1][0][0][0] << "\t" << bestchi3[k].controlcounts[2][2][0][0][0] << "\n\n";

				  if (test > 2)
				    {
				      bestMarkerCombi2Details << "\nL1:\n";
				      if(printBeta) resultReg(bestchi3[k].result1, bestMarkerCombi2Details,haplo, qt,maxIndexCov);
				      //if(printBeta && test==18) {bestMarkerCombi2Details << "\nL2:\n";}
				      if(printBeta && test==18) resultReg(bestchi3[k].result3, bestMarkerCombi2Details,haplo, qt,maxIndexCov);
				      bestMarkerCombi2Details << "\n\n";
				    }
				}

			      else // qt
				{
				  // Ausgabe lineare Regression
				  bestMarkerCombi2Details << k+1 << ": SNP1: " << map[bestchi3[k].nr1].rs << "\tSNP2: " << map[bestchi3[k].nr2].rs
							  << "\n\n";
				  bestMarkerCombi2Details << "counts:\n";
				  bestMarkerCombi2Details << bestchi3[k].casecounts[0][0][0][0][0] << "\t" << bestchi3[k].casecounts[0][1][0][0][0]
							  << "\t" << bestchi3[k].casecounts[0][2][0][0][0] << "\n"
							  << bestchi3[k].casecounts[1][0][0][0][0] << "\t" << bestchi3[k].casecounts[1][1][0][0][0] << "\t"
							  << bestchi3[k].casecounts[1][2][0][0][0] << "\n"
							  << bestchi3[k].casecounts[2][0][0][0][0] << "\t" << bestchi3[k].casecounts[2][1][0][0][0] << "\t"
							  << bestchi3[k].casecounts[2][2][0][0][0] << "\n\n";

				  if (test > 2)
				    {

				      bestMarkerCombi2Details << "\nL1:\n";
				      if(printBeta) resultReg(bestchi3[k].result1, bestMarkerCombi2Details,haplo, qt,maxIndexCov);
				      //bestMarkerCombi2Details << "\nL2:\n";
				      //if(printBeta) resultReg(bestchi3[k].result2, bestMarkerCombi2Details,haplo, qt,maxIndexCov);
				      bestMarkerCombi2Details << "\n";
				      bestMarkerCombi2Details << "L1 rsquare: " << bestchi3[k].result1.rsquare << "\t" /*<< "L2 rsquare: " << bestchi3[k].result2.rsquare << "\t"*/ << "Fstat: " << bestchi3[k].Fstat << "\n\n";

				    }

				  if(pretest && (test ==3 || test==5))
				    {
				      bestMarkerCombi2Details << "avgSnp1\t" << bestchi3[k].traitavg.snp1allele[0] << "\t" << bestchi3[k].traitavg.snp1allele[1] << "\n";
				      bestMarkerCombi2Details << "avgSnp2\t" << bestchi3[k].traitavg.snp2allele[0] << "\t" << bestchi3[k].traitavg.snp2allele[1] << "\n";
				      bestMarkerCombi2Details << "avgAlleleCombinations\t" << bestchi3[k].traitavg.allele[0][0] << "\t" << bestchi3[k].traitavg.allele[0][1] << "\t" << bestchi3[k].traitavg.allele[1][0] << "\t" << bestchi3[k].traitavg.allele[1][1] << "\n";
				    }
				  else if(pretest && (test ==4 || test==6))
				    {
				      bestMarkerCombi2Details << "avgSnp1\t" << bestchi3[k].traitavg.snp1geno[0] << "\t" << bestchi3[k].traitavg.snp1geno[1] << "\t" << bestchi3[k].traitavg.snp1geno[2] << "\n";
				      bestMarkerCombi2Details << "avgSnp2\t" << bestchi3[k].traitavg.snp2geno[0] << "\t" << bestchi3[k].traitavg.snp2geno[1] << "\t" << bestchi3[k].traitavg.snp2geno[2] << "\n";
				      bestMarkerCombi2Details << "avgGenoCombinations\t" << bestchi3[k].traitavg.geno[0][0];
				      bestMarkerCombi2Details << "\t" << bestchi3[k].traitavg.geno[0][1];
				      bestMarkerCombi2Details << "\t" << bestchi3[k].traitavg.geno[0][2];
				      bestMarkerCombi2Details << "\t" << bestchi3[k].traitavg.geno[1][0];
				      bestMarkerCombi2Details << "\t" << bestchi3[k].traitavg.geno[1][1];
				      bestMarkerCombi2Details << "\t" << bestchi3[k].traitavg.geno[1][2];
				      bestMarkerCombi2Details << "\t" << bestchi3[k].traitavg.geno[2][0];
				      bestMarkerCombi2Details << "\t" << bestchi3[k].traitavg.geno[2][1];
				      bestMarkerCombi2Details << "\t" << bestchi3[k].traitavg.geno[2][2];
				      bestMarkerCombi2Details << "\n\n";
				    }
				}
			    }
			  else
			    {
			      if (!qt)
				{
				  // male
				  bestMarkerCombi2Details << k+1 << ": SNP1: " << map[bestchi3[k].nr1].rs << "\tSNP2: " << map[bestchi3[k].nr2].rs
							  << "\n\n";
				  bestMarkerCombi2Details << "male: \n";
				  bestMarkerCombi2Details << "cases:\t\t\t\tcontrols:\t\t\n";
				  bestMarkerCombi2Details << bestchi3[k].casecountsMale[0][0][0][0][0] << "\t"
							  << bestchi3[k].casecountsMale[0][1][0][0][0] << "\t" << bestchi3[k].casecountsMale[0][2][0][0][0] << "\t\t\t"
							  << bestchi3[k].controlcountsMale[0][0][0][0][0] << "\t" << bestchi3[k].controlcountsMale[0][1][0][0][0] << "\t"
							  << bestchi3[k].controlcountsMale[0][2][0][0][0] << "\n" << bestchi3[k].casecountsMale[1][0][0][0][0] << "\t"
							  << bestchi3[k].casecountsMale[1][1][0][0][0] << "\t" << bestchi3[k].casecountsMale[1][2][0][0][0] << "\t\t\t"
							  << bestchi3[k].controlcountsMale[1][0][0][0][0] << "\t" << bestchi3[k].controlcountsMale[1][1][0][0][0] << "\t"
							  << bestchi3[k].controlcountsMale[1][2][0][0][0] << "\n" << bestchi3[k].casecountsMale[2][0][0][0][0] << "\t"
							  << bestchi3[k].casecountsMale[2][1][0][0][0] << "\t" << bestchi3[k].casecountsMale[2][2][0][0][0] << "\t\t\t"
							  << bestchi3[k].controlcountsMale[2][0][0][0][0] << "\t" << bestchi3[k].controlcountsMale[2][1][0][0][0] << "\t"
							  << bestchi3[k].controlcountsMale[2][2][0][0][0] << "\n";

				  // female
				  bestMarkerCombi2Details << "female: \n";
				  bestMarkerCombi2Details << "cases:\t\t\t\tcontrols:\t\t\n";
				  bestMarkerCombi2Details << bestchi3[k].casecountsFemale[0][0][0][0][0] << "\t"
							  << bestchi3[k].casecountsFemale[0][1][0][0][0] << "\t" << bestchi3[k].casecountsFemale[0][2][0][0][0] << "\t\t\t"
							  << bestchi3[k].controlcountsFemale[0][0][0][0][0] << "\t" << bestchi3[k].controlcountsFemale[0][1][0][0][0] << "\t" << bestchi3[k].controlcountsFemale[0][2][0][0][0] << "\n" << bestchi3[k].casecountsFemale[1][0][0][0][0] << "\t"
							  << bestchi3[k].casecountsFemale[1][1][0][0][0] << "\t" << bestchi3[k].casecountsFemale[1][2][0][0][0]  << "\t\t\t"
							  << bestchi3[k].controlcountsFemale[1][0][0][0][0] << "\t" << bestchi3[k].controlcountsFemale[1][1][0][0][0] << "\t" << bestchi3[k].controlcountsFemale[1][2][0][0][0] << "\n" << bestchi3[k].casecountsFemale[2][0][0][0][0] << "\t"
							  << bestchi3[k].casecountsFemale[2][1][0][0][0] << "\t" << bestchi3[k].casecountsFemale[2][2][0][0][0]  << "\t\t\t"
							  << bestchi3[k].controlcountsFemale[2][0][0][0][0] << "\t" << bestchi3[k].controlcountsFemale[2][1][0][0][0] << "\t" << bestchi3[k].controlcountsFemale[2][2][0][0][0] << "\n";

				  if (test > 2)
				    {
				      bestMarkerCombi2Details << "\nL1:\n";
				      if(printBeta) resultReg(bestchi3[k].result1, bestMarkerCombi2Details,haplo, qt,maxIndexCov);
				      //if(printBeta && test==18) {bestMarkerCombi2Details << "\nL2:\n";}
				      if(printBeta && test==18) resultReg(bestchi3[k].result3, bestMarkerCombi2Details,haplo, qt,maxIndexCov);
				      bestMarkerCombi2Details << "\n\n";
				    }
				}
			      else
				{
				  // Ausgabe lineare Regression
				  bestMarkerCombi2Details << k+1 << ": SNP1: " << map[bestchi3[k].nr1].rs << "\tSNP2: " << map[bestchi3[k].nr2].rs
							  << "\n\n";
				  bestMarkerCombi2Details << "counts:\n";
				  bestMarkerCombi2Details << bestchi3[k].casecounts[0][0][0][0][0] << "\t" << bestchi3[k].casecounts[0][1][0][0][0]
							  << "\t" << bestchi3[k].casecounts[0][2][0][0][0] << "\n"
							  << bestchi3[k].casecounts[1][0][0][0][0] << "\t" << bestchi3[k].casecounts[1][1][0][0][0] << "\t"
							  << bestchi3[k].casecounts[1][2][0][0][0] << "\n"
							  << bestchi3[k].casecounts[2][0][0][0][0] << "\t" << bestchi3[k].casecounts[2][1][0][0][0] << "\t"
							  << bestchi3[k].casecounts[2][2][0][0][0] << "\n\n";

				  if (test > 2)
				    {

				      bestMarkerCombi2Details << "\nL1:\n";
				      if(printBeta) resultReg(bestchi3[k].result1, bestMarkerCombi2Details,haplo, qt,maxIndexCov);
				      //bestMarkerCombi2Details << "\nL2:\n";
				      //if(printBeta) resultReg(bestchi3[k].result2, bestMarkerCombi2Details,haplo, qt,maxIndexCov);
				      bestMarkerCombi2Details << "\n";
				      bestMarkerCombi2Details << "L1 rsquare: " << bestchi3[k].result1.rsquare << "\t" /* << "L2 rsquare: " << bestchi3[k].result2.rsquare << "\t" */ << "Fstat: " << bestchi3[k].Fstat << "\n\n";

				    }

				  if(pretest && (test ==3 || test==5))
				    {
				      bestMarkerCombi2Details << "avgSnp1\t" << bestchi3[k].traitavg.snp1allele[0] << "\t" << bestchi3[k].traitavg.snp1allele[1] << "\n";
				      bestMarkerCombi2Details << "avgSnp2\t" << bestchi3[k].traitavg.snp2allele[0] << "\t" << bestchi3[k].traitavg.snp2allele[1] << "\n";
				      bestMarkerCombi2Details << "avgAlleleCombinations\t" << bestchi3[k].traitavg.allele[0][0] << "\t" << bestchi3[k].traitavg.allele[0][1] << "\t" << bestchi3[k].traitavg.allele[1][0] << "\t" << bestchi3[k].traitavg.allele[1][1] << "\n";
				    }
				  else if(pretest && (test ==4 || test==6))
				    {
				      bestMarkerCombi2Details << "avgSnp1\t" << bestchi3[k].traitavg.snp1geno[0] << "\t" << bestchi3[k].traitavg.snp1geno[1] << "\t" << bestchi3[k].traitavg.snp1geno[2] << "\n";
				      bestMarkerCombi2Details << "avgSnp2\t" << bestchi3[k].traitavg.snp2geno[0] << "\t" << bestchi3[k].traitavg.snp2geno[1] << "\t" << bestchi3[k].traitavg.snp2geno[2] << "\n";
				      bestMarkerCombi2Details << "avgGenoCombinations\t" << bestchi3[k].traitavg.geno[0][0];
				      bestMarkerCombi2Details << "\t" << bestchi3[k].traitavg.geno[0][1];
				      bestMarkerCombi2Details << "\t" << bestchi3[k].traitavg.geno[0][2];
				      bestMarkerCombi2Details << "\t" << bestchi3[k].traitavg.geno[1][0];
				      bestMarkerCombi2Details << "\t" << bestchi3[k].traitavg.geno[1][1];
				      bestMarkerCombi2Details << "\t" << bestchi3[k].traitavg.geno[1][2];
				      bestMarkerCombi2Details << "\t" << bestchi3[k].traitavg.geno[2][0];
				      bestMarkerCombi2Details << "\t" << bestchi3[k].traitavg.geno[2][1];
				      bestMarkerCombi2Details << "\t" << bestchi3[k].traitavg.geno[2][2];
				      bestMarkerCombi2Details << "\n\n";
				    }
				}
			    }
			}
		    }
		  //if(combilist){bestMarkerCombi2 << "\n\n Note: Each combination is listed only once, even if it occurs multiple times in the combilist\n e.g.: (SNP1, SNP2) = (SNP2, SNP1)\n";}
		  bestMarkerCombi2.close();
		  bestMarkerCombi2Details.close();

		}
	      else if (markercombi3 == 1)
		{
		  bestMarkerCombi3.open(markerCombi3file.c_str(),ios::out);
		  bestMarkerCombi3Details.open(markerCombi3Detailsfile.c_str(), ios::out);

		  // 10 besten P-Werte werden ausgeben
		  cout << "\n3-marker: top p-values multi-marker analysis:\n";
		  logfile << "\n3-marker: top p-values multi-marker analysis:\n";

		  if (pathwayImpact == 1)
		    {
		      if (test < 3)
			{
			  cout << "No Chr_No_1 rs_No_1 SNP_Pos_1 p-Single-marker_1 Chr_No_2 rs_No_2 Pos_No_2 p-Single-marker_2 Chr_No_3 rs_No_3 Pos_No_3 p-Single-marker_3 p-value p-value_corr Pathway\n";
			  logfile << "No Chr_No_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_No_2 rs_No_2 Pos_No_2 p-Single-marker_2 Chr_No_3 rs_No_3 Pos_No_3 p-Single-marker_3 p-value p-value_corr Pathway\n";
			  bestMarkerCombi3 << "No\tChr_No_1\trs_No_1\tPos_No_1\tp-Single-marker_1\tChr_No_2\trs_No_2\tPos_No_2\tp-Single-marker_2\tChr_No_3\trs_No_3\tPos_No_3\tp-Single-marker_3\tp-value\tp-value_corr\tPathway\n";
			}
		      else
			{
			  cout << "No Chr_No_1 rs_No_1 SNP_Pos_1 p-Single-marker_1 Chr_No_2 rs_No_2 Pos_No_2 p-Single-marker_2 Chr_No_3 rs_No_3 Pos_No_3 p-Single-marker_3 p-value p-value_corr Pathway\n";
			  logfile << "No Chr_No_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_No_2 rs_No_2 Pos_No_2 p-Single-marker_2 Chr_No_3 rs_No_3 Pos_No_3 p-Single-marker_3 p-value p-value_corr Pathway\n";
			  bestMarkerCombi3 << "No\tChr_No_1\trs_No_1\tPos_No_1\tp-Single-marker_1\tChr_No_2\trs_No_2\tPos_No_2\tp-Single-marker_2\tChr_No_3\trs_No_3\tPos_No_3\tp-Single-marker_3\tp-value\tp-value_corr\tPathway\tSNP1_Allele_A\tSNP1_Allele_B\tSNP2_Allele_A\tSNP2_Allele_B\tSNP3_Allele_A\tSNP3_Allele_B\tbeta_x1\tse_x1\tbeta_x1D\tse_x1D\tbeta_x2\tse_x2\tbeta_x2D\tse_x2D\tbeta_x1x2\tse_x1x2\tbeta_x1x2D\tse_x1x2D\tbeta_x1Dx2\tse_x1Dx2\tbeta_x1Dx2D\tse_x1Dx2D\tbeta_x3\tse_x3\tbeta_x3D\tse_x3D\tbeta_x1x3\tse_x1x3\tbeta_x1x3D\tse_x1x3D\tbeta_x1Dx3\tse_x1Dx3\tbeta_x1Dx3D\tse_x1Dx3D\tbeta_x2x3\tse_x2x3\tbeta_x2x3D\tse_x2x3D\tbeta_x2Dx3\tse_x2Dx3\tbeta_x2Dx3D\tse_x2Dx3D\tbeta_x1x2x3\tse_x1x2x3\tbeta_x1x2x3D\tse_x1x2x3D\tbeta_x1x2Dx3\tse_x1x2Dx3\tbeta_x1x2Dx3D\tse_x1x2Dx3D\tbeta_x1Dx2x3\tse_x1Dx2x3\tbeta_x1Dx2x3D\tse_x1Dx2x3D\tbeta_x1Dx2Dx3\tse_x1Dx2Dx3\tbeta_x1Dx2Dx3D\tse_x1Dx2Dx3D\n";
			}
		    }
		  else
		    {
		      if (test < 3)
			{
			  cout << "No Chr_No_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_No_2 rs_No_2 Pos_No_2 p-Single-marker_2 Chr_No_3 rs_No_3 Pos_No_3 p-Single-marker_3 p-value p-value_corr\n";
			  logfile << "No Chr_No_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_No_2 rs_No_2 Pos_No_2 p-Single-marker_2 Chr_No_3 rs_No_3 Pos_No_3 p-Single-marker_3 p-value p-value_corr\n";
			  bestMarkerCombi3 << "No\tChr_No_1\trs_No_1\tPos_No_1\tp-Single-marker_1\tChr_No_2\trs_No_2\tPos_No_2\tp-Single-marker_2\tChr_No_3\trs_No_3\tPos_No_3\tp-Single-marker_3\tp-value\tp-value_corr\n";
			}
		      else if(!haplo)
			{
			  cout << "No Chr_No_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_No_2 rs_No_2 Pos_No_2 p-Single-marker_2 Chr_No_3 rs_No_3 Pos_No_3 p-Single-marker_3 p-value p-value_corr\n";
			  logfile << "No Chr_No_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_No_2 rs_No_2 Pos_No_2 p-Single-marker_2 Chr_No_3 rs_No_3 Pos_No_3 p-Single-marker_3 p-value p-value_corr\n";
			  bestMarkerCombi3 << "No\tChr_No_1\trs_No_1\tPos_No_1\tp-Single-marker_1\tChr_No_2\trs_No_2\tPos_No_2\tp-Single-marker_2\tChr_No_3\trs_No_3\tPos_No_3\tp-Single-marker_3\tp-value\tp-value_corr\tSNP1_Allele_A\tSNP1_Allele_B\tSNP2_Allele_A\tSNP2_Allele_B\tSNP3_Allele_A\tSNP3_Allele_B\tbeta_x1\tse_x1\tbeta_x1D\tse_x1D\tbeta_x2\tse_x2\tbeta_x2D\tse_x2D\tbeta_x1x2\tse_x1x2\tbeta_x1x2D\tse_x1x2D\tbeta_x1Dx2\tse_x1Dx2\tbeta_x1Dx2D\tse_x1Dx2D\tbeta_x3\tse_x3\tbeta_x3D\tse_x3D\tbeta_x1x3\tse_x1x3\tbeta_x1x3D\tse_x1x3D\tbeta_x1Dx3\tse_x1Dx3\tbeta_x1Dx3D\tse_x1Dx3D\tbeta_x2x3\tse_x2x3\tbeta_x2x3D\tse_x2x3D\tbeta_x2Dx3\tse_x2Dx3\tbeta_x2Dx3D\tse_x2Dx3D\tbeta_x1x2x3\tse_x1x2x3\tbeta_x1x2x3D\tse_x1x2x3D\tbeta_x1x2Dx3\tse_x1x2Dx3\tbeta_x1x2Dx3D\tse_x1x2Dx3D\tbeta_x1Dx2x3\tse_x1Dx2x3\tbeta_x1Dx2x3D\tse_x1Dx2x3D\tbeta_x1Dx2Dx3\tse_x1Dx2Dx3\tbeta_x1Dx2Dx3D\tse_x1Dx2Dx3D\n";
			}
		      else if(haplo)
			{
			  cout << "No Chr_No_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_No_2 rs_No_2 Pos_No_2 p-Single-marker_2 Chr_No_3 rs_No_3 Pos_No_3 p-Single-marker_3 p-value p-value_corr\n";
			  logfile << "No Chr_No_1 rs_No_1 Pos_No_1 p-Single-marker_1 Chr_No_2 rs_No_2 Pos_No_2 p-Single-marker_2 Chr_No_3 rs_No_3 Pos_No_3 p-Single-marker_3 p-value p-value_corr\n";
			  bestMarkerCombi3 << "No\tChr_No_1\trs_No_1\tPos_No_1\tp-Single-marker_1\tChr_No_2\trs_No_2\tPos_No_2\tp-Single-marker_2\tChr_No_3\trs_No_3\tPos_No_3\tp-Single-marker_3\tp-value\tp-value_corr\tSNP1_Allele_A\tSNP1_Allele_B\tSNP2_Allele_A\tSNP2_Allele_B\tSNP3_Allele_A\tSNP3_Allele_B\tbeta_h1\tse_h1\tbeta_h2\tse_h2\tbeta_h3\tse_h3\tbeta_h4\tse_h4\tbeta_h5\tse_h5\tbeta_h6\tse_h6\tbeta_h7\tse_h7\tbeta_h8\tse_h8\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse\n";
			}
		    }

		  withmarginal = xvec[1]+xvec[2]+xvec[3]+xvec[4]+xvec[9]+xvec[10]-zvec[1]-zvec[2]-zvec[3]-zvec[4]-zvec[9]-zvec[10];

		  if(mWithSingletop==1)
		    {
		      withmarginal=withmarginal-(xvec[3]-zvec[3])-(xvec[4]-zvec[4])-(xvec[9]-zvec[9])-(xvec[10]-zvec[10]);
		    }
		  else if(mWithSingletop==2)
		    {
		      withmarginal=withmarginal-(xvec[9]-zvec[9])-(xvec[10]-zvec[10]);
		    }

		  if(haplo){withmarginal=0;}

		  for (k = 0; k < printtop; k++)
		    {
		      // Korrektur

		      if (mWithSingletop == 0 || combilist == 1 || withmarginal==0 || haplo || caseOnly && !(plusSingle && mWithSingletop>0) )
			{
			  bestchi3[k].pmod = Sidak(bestchi3[k].p,ntests);

			}
		      else
			{
			  bestchi3[k].pmod = Sidak(bestchi3[k].p,ntestsPlus);
			}

		      if (bestchi3[k].nr1 != -1 && bestchi3[k].nr2 != -1
			  && bestchi3[k].nr3 != -1 && bestchi3[k].nr4 != -1
			  && bestchi3[k].nr5 != -1)
			{
			  if (pathwayImpact == 1)
			    {
			      if (k < 10)
				{
				  cout << k + 1 << " " <<  map[bestchi3[k].nr1].chr << " " << map[bestchi3[k].nr1].rs << " " << map[bestchi3[k].nr1].pos << " "
				       << map[bestchi3[k].nr1].p << " " << map[bestchi3[k].nr2].chr << " " << map[bestchi3[k].nr2].rs << " "
				       << map[bestchi3[k].nr2].pos << " " << map[bestchi3[k].nr2].p << " " << map[bestchi3[k].nr3].chr << " "
				       << map[bestchi3[k].nr3].rs << " " << map[bestchi3[k].nr3].pos << " " << map[bestchi3[k].nr3].p << " "
				       << bestchi3[k].p << " " << bestchi3[k].pmod << " " << pathway[bestchi3[k].r-1].name << "\n";

				  logfile << k + 1 << " " <<  map[bestchi3[k].nr1].chr << " " << map[bestchi3[k].nr1].rs << " "
					  << map[bestchi3[k].nr1].pos << " " << map[bestchi3[k].nr1].p << " " << map[bestchi3[k].nr2].chr << " "
					  << map[bestchi3[k].nr2].rs << " " << map[bestchi3[k].nr2].pos << " " << map[bestchi3[k].nr2].p << " "
					  << map[bestchi3[k].nr3].chr << " " << map[bestchi3[k].nr3].rs << " " << map[bestchi3[k].nr3].pos << " "
					  << map[bestchi3[k].nr3].p << " " << bestchi3[k].p << " " << bestchi3[k].pmod << " " << pathway[bestchi3[k].r-1].name << "\n";
				}

			      bestMarkerCombi3 << k + 1 << "\t" <<  map[bestchi3[k].nr1].chr << "\t" << map[bestchi3[k].nr1].rs << "\t"
					       << map[bestchi3[k].nr1].pos << "\t" << map[bestchi3[k].nr1].p << "\t" << map[bestchi3[k].nr2].chr << "\t"
					       << map[bestchi3[k].nr2].rs << "\t" << map[bestchi3[k].nr2].pos << "\t" << map[bestchi3[k].nr2].p << "\t"
					       << map[bestchi3[k].nr3].chr << "\t" << map[bestchi3[k].nr3].rs << "\t" << map[bestchi3[k].nr3].pos << "\t"
					       << map[bestchi3[k].nr3].p << "\t" << bestchi3[k].p << "\t" << bestchi3[k].pmod << "\t" << pathway[bestchi3[k].r-1].name << "\n";
			    }
			  else
			    {
			      if (k < 10)
				{
				  cout << k + 1 << " "
				       <<  map[bestchi3[k].nr1].chr << " "
				       << map[bestchi3[k].nr1].rs << " "
				       << map[bestchi3[k].nr1].pos << " "
				       << map[bestchi3[k].nr1].p << " "
				       << map[bestchi3[k].nr2].chr << " "
				       << map[bestchi3[k].nr2].rs << " "
				       << map[bestchi3[k].nr2].pos << " "
				       << map[bestchi3[k].nr2].p << " "
				       << map[bestchi3[k].nr3].chr << " "
				       << map[bestchi3[k].nr3].rs << " "
				       << map[bestchi3[k].nr3].pos << " "
				       << map[bestchi3[k].nr3].p << " "
				       << bestchi3[k].p << " "
				       << bestchi3[k].pmod << " "
				       << "\n";

				  logfile << k + 1 << " "
					  <<  map[bestchi3[k].nr1].chr << " "
					  << map[bestchi3[k].nr1].rs << " "
					  << map[bestchi3[k].nr1].pos << " "
					  << map[bestchi3[k].nr1].p << " "
					  << map[bestchi3[k].nr2].chr << " "
					  << map[bestchi3[k].nr2].rs << " "
					  << map[bestchi3[k].nr2].pos << " "
					  << map[bestchi3[k].nr2].p << " "
					  << map[bestchi3[k].nr3].chr << " "
					  << map[bestchi3[k].nr3].rs << " "
					  << map[bestchi3[k].nr3].pos << " "
					  << map[bestchi3[k].nr3].p << " "
					  << bestchi3[k].p << " "
					  << bestchi3[k].pmod << " "
					  << "\n";
				}


			      bestMarkerCombi3 << k + 1 << "\t"
					       <<  map[bestchi3[k].nr1].chr << "\t" << map[bestchi3[k].nr1].rs << "\t" << map[bestchi3[k].nr1].pos << "\t"
					       << map[bestchi3[k].nr1].p << "\t" << map[bestchi3[k].nr2].chr << "\t" << map[bestchi3[k].nr2].rs << "\t"
					       << map[bestchi3[k].nr2].pos << "\t" << map[bestchi3[k].nr2].p << "\t" << map[bestchi3[k].nr3].chr << "\t"
					       << map[bestchi3[k].nr3].rs << "\t" << map[bestchi3[k].nr3].pos << "\t" << map[bestchi3[k].nr3].p << "\t"
					       << bestchi3[k].p << "\t" << bestchi3[k].pmod;

			      if (test < 3)
				{
				  bestMarkerCombi3 << "\n";
				}
			      else
				{
				  bestMarkerCombi3 << "\t" << codesA[bestchi3[k].nr1].a1 << "\t" << codesA[bestchi3[k].nr1].a2 << "\t" << codesA[bestchi3[k].nr2].a1 << "\t" << codesA[bestchi3[k].nr2].a2 << "\t" << codesA[bestchi3[k].nr3].a1 << "\t" << codesA[bestchi3[k].nr3].a2;

				  for (i =1; i<=26; i++)
				    {
				      if (xvec[i] != zvec[i] /*&& i != 10 && i != 12 && i != 13 && i != 14 && i != 16 && i != 17 && i != 18*/)
					{
					  bestMarkerCombi3 << "\t" << bestchi3[k].result1.b[i] << "\t" << bestchi3[k].result1.betaNew_se[i];
					}
				      else //if (i != 10 && i != 12 && i != 13 && i != 14 && i != 16 && i != 17 && i != 18)
					{
					  bestMarkerCombi3 << "\t-\t-";
					}

				    }
				  bestMarkerCombi3 << "\n";
				}
			    }


			  // Details
			  if (bestchi3[k].X == 0)
			    {
			      if (!qt)
				{
				  bestMarkerCombi3Details << k+1 << ": SNP1: " << map[bestchi3[k].nr1].rs << "\tSNP2: " << map[bestchi3[k].nr2].rs
							  << "\tSNP3: " << map[bestchi3[k].nr3].rs << "\n\n";
				  bestMarkerCombi3Details << "cases:\t\t\t\tcontrols:\t\t\n";
				  bestMarkerCombi3Details << bestchi3[k].casecounts[0][0][0][0][0] << "\t" << bestchi3[k].casecounts[0][0][1][0][0]
							  << "\t" << bestchi3[k].casecounts[0][0][2][0][0] << "\t\t\t" << bestchi3[k].controlcounts[0][0][0][0][0] << "\t"
							  << bestchi3[k].controlcounts[0][0][1][0][0] << "\t" << bestchi3[k].controlcounts[0][0][2][0][0] << "\n"
							  << bestchi3[k].casecounts[0][1][0][0][0] << "\t" << bestchi3[k].casecounts[0][1][1][0][0] << "\t"
							  << bestchi3[k].casecounts[0][1][2][0][0] << "\t\t\t" << bestchi3[k].controlcounts[0][1][0][0][0] << "\t"
							  << bestchi3[k].controlcounts[0][1][1][0][0] << "\t" << bestchi3[k].controlcounts[0][1][2][0][0] << "\n"
							  << bestchi3[k].casecounts[0][2][0][0][0] << "\t" << bestchi3[k].casecounts[0][2][1][0][0] << "\t"
							  << bestchi3[k].casecounts[0][2][2][0][0] << "\t\t\t" << bestchi3[k].controlcounts[0][2][0][0][0] << "\t"

							  << bestchi3[k].controlcounts[0][2][1][0][0] << "\t" << bestchi3[k].controlcounts[0][2][2][0][0] << "\n"
							  << bestchi3[k].casecounts[1][0][0][0][0] << "\t" << bestchi3[k].casecounts[1][0][1][0][0] << "\t"
							  << bestchi3[k].casecounts[1][0][2][0][0] << "\t\t\t" << bestchi3[k].controlcounts[1][0][0][0][0] << "\t"
							  << bestchi3[k].controlcounts[1][0][1][0][0] << "\t" << bestchi3[k].controlcounts[1][0][2][0][0] << "\n"
							  << bestchi3[k].casecounts[1][1][0][0][0] << "\t" << bestchi3[k].casecounts[1][1][1][0][0] << "\t"
							  << bestchi3[k].casecounts[1][1][2][0][0] << "\t\t\t" << bestchi3[k].controlcounts[1][1][0][0][0] << "\t"
							  << bestchi3[k].controlcounts[1][1][1][0][0] << "\t" << bestchi3[k].controlcounts[1][1][2][0][0] << "\n"
							  << bestchi3[k].casecounts[1][2][0][0][0] << "\t" << bestchi3[k].casecounts[1][2][1][0][0] << "\t"
							  << bestchi3[k].casecounts[1][2][2][0][0] << "\t\t\t" << bestchi3[k].controlcounts[1][2][0][0][0] << "\t"
							  << bestchi3[k].controlcounts[1][2][1][0][0] << "\t" << bestchi3[k].controlcounts[1][2][2][0][0] << "\n"
							  << bestchi3[k].casecounts[2][0][0][0][0] << "\t" << bestchi3[k].casecounts[2][0][1][0][0] << "\t"
							  << bestchi3[k].casecounts[2][0][2][0][0] << "\t\t\t" << bestchi3[k].controlcounts[2][0][0][0][0] << "\t"
							  << bestchi3[k].controlcounts[2][0][1][0][0] << "\t"  << bestchi3[k].controlcounts[2][0][2][0][0] << "\n"
							  << bestchi3[k].casecounts[2][1][0][0][0] << "\t" << bestchi3[k].casecounts[2][1][1][0][0] << "\t"
							  << bestchi3[k].casecounts[2][1][2][0][0] << "\t\t\t" << bestchi3[k].controlcounts[2][1][0][0][0] << "\t"
							  << bestchi3[k].controlcounts[2][1][1][0][0] << "\t" << bestchi3[k].controlcounts[2][1][2][0][0] << "\n"
							  << bestchi3[k].casecounts[2][2][0][0][0] << "\t" << bestchi3[k].casecounts[2][2][1][0][0] << "\t"
							  << bestchi3[k].casecounts[2][2][2][0][0] << "\t\t\t" << bestchi3[k].controlcounts[2][2][0][0][0] << "\t"
							  << bestchi3[k].controlcounts[2][2][1][0][0] << "\t" << bestchi3[k].controlcounts[2][2][2][0][0] << "\n\n";

				  if (test > 2)
				    {
				      bestMarkerCombi3Details << "\nL1:\n";
				      if(printBeta) resultReg(bestchi3[k].result1, bestMarkerCombi3Details,haplo, qt,maxIndexCov);
				      //bestMarkerCombi3Details << "\nL2:\n";
				      //resultReg(bestchi3[k].result2, bestMarkerCombi3Details,haplo, qt);
				      bestMarkerCombi3Details << "\n\n";
				    }
				}
			      else
				{
				  bestMarkerCombi3Details << k+1 << ": SNP1: " << map[bestchi3[k].nr1].rs << "\tSNP2: " << map[bestchi3[k].nr2].rs
							  << "\tSNP3: rs" << map[bestchi3[k].nr3].rs << "\n\n";
				  bestMarkerCombi3Details << "counts:\n";
				  bestMarkerCombi3Details << bestchi3[k].casecounts[0][0][0][0][0] << "\t" << bestchi3[k].casecounts[0][0][1][0][0]
							  << "\t" << bestchi3[k].casecounts[0][0][2][0][0] << "\n"
							  << bestchi3[k].casecounts[0][1][0][0][0] << "\t" << bestchi3[k].casecounts[0][1][1][0][0] << "\t"
							  << bestchi3[k].casecounts[0][1][2][0][0] << "\n"
							  << bestchi3[k].casecounts[0][2][0][0][0] << "\t" << bestchi3[k].casecounts[0][2][1][0][0] << "\t"
							  << bestchi3[k].casecounts[0][2][2][0][0] << "\n"

							  << bestchi3[k].casecounts[1][0][0][0][0] << "\t" << bestchi3[k].casecounts[1][0][1][0][0] << "\t"
							  << bestchi3[k].casecounts[1][0][2][0][0] << "\n"
							  << bestchi3[k].casecounts[1][1][0][0][0] << "\t" << bestchi3[k].casecounts[1][1][1][0][0] << "\t"
							  << bestchi3[k].casecounts[1][1][2][0][0] << "\n"
							  << bestchi3[k].casecounts[1][2][0][0][0] << "\t" << bestchi3[k].casecounts[1][2][1][0][0] << "\t"

							  << bestchi3[k].casecounts[1][2][2][0][0] << "\n"
							  << bestchi3[k].casecounts[2][0][0][0][0] << "\t" << bestchi3[k].casecounts[2][0][1][0][0] << "\t"
							  << bestchi3[k].casecounts[2][0][2][0][0] << "\n"
							  << bestchi3[k].casecounts[2][1][0][0][0] << "\t" << bestchi3[k].casecounts[2][1][1][0][0] << "\t"
							  << bestchi3[k].casecounts[2][1][2][0][0] << "\n"
							  << bestchi3[k].casecounts[2][2][0][0][0] << "\t" << bestchi3[k].casecounts[2][2][1][0][0] << "\t"
							  << bestchi3[k].casecounts[2][2][2][0][0] << "\n\n";

				  if (test > 2)
				    {
				      bestMarkerCombi3Details << "\nL1:\n";
				      if(printBeta) resultReg(bestchi3[k].result1, bestMarkerCombi3Details,haplo, qt,maxIndexCov);
				      //bestMarkerCombi3Details << "\nL2:\n";
				      //resultReg(bestchi3[k].result2, bestMarkerCombi3Details,haplo, qt);
				      bestMarkerCombi3Details << "\n\n";
				    }
				}
			    }
			  else
			    {
			      // male
			      bestMarkerCombi3Details << k+1 << ": SNP1: " << map[bestchi3[k].nr1].rs << "\tSNP2: " << map[bestchi3[k].nr2].rs
						      << "\tSNP3: " << map[bestchi3[k].nr3].rs << "\n\n";
			      bestMarkerCombi3Details << "male: \n";
			      bestMarkerCombi3Details << "cases:\t\t\t\tcontrols:\t\t\n";
			      bestMarkerCombi3Details << bestchi3[k].casecountsMale[0][0][0][0][0] << "\t"
						      << bestchi3[k].casecountsMale[0][0][1][0][0] << "\t" << bestchi3[k].casecountsMale[0][0][2][0][0] << "\t\t\t"
						      << bestchi3[k].controlcountsMale[0][0][0][0][0] << "\t" << bestchi3[k].controlcountsMale[0][0][1][0][0] << "\t"
						      << bestchi3[k].controlcountsMale[0][0][2][0][0] << "\n" << bestchi3[k].casecountsMale[0][1][0][0][0] << "\t"
						      << bestchi3[k].casecountsMale[0][1][1][0][0] << "\t" << bestchi3[k].casecountsMale[0][1][2][0][0] << "\t\t\t"
						      << bestchi3[k].controlcountsMale[0][1][0][0][0] << "\t" << bestchi3[k].controlcountsMale[0][1][1][0][0] << "\t"
						      << bestchi3[k].controlcountsMale[0][1][2][0][0] << "\n" << bestchi3[k].casecountsMale[0][2][0][0][0] << "\t"
						      << bestchi3[k].casecountsMale[0][2][1][0][0] << "\t" << bestchi3[k].casecountsMale[0][2][2][0][0] << "\t\t\t"
						      << bestchi3[k].controlcountsMale[0][2][0][0][0] << "\t" << bestchi3[k].controlcountsMale[0][2][1][0][0] << "\t"
						      << bestchi3[k].controlcountsMale[0][2][2][0][0] << "\n" << bestchi3[k].casecountsMale[1][0][0][0][0] << "\t"
						      << bestchi3[k].casecountsMale[1][0][1][0][0] << "\t" << bestchi3[k].casecountsMale[1][0][2][0][0] << "\t\t\t"
						      << bestchi3[k].controlcountsMale[1][0][0][0][0] << "\t" << bestchi3[k].controlcountsMale[1][0][1][0][0] << "\t"
						      << bestchi3[k].controlcountsMale[1][0][2][0][0] << "\n" << bestchi3[k].casecountsMale[1][1][0][0][0] << "\t"
						      << bestchi3[k].casecountsMale[1][1][1][0][0] << "\t" << bestchi3[k].casecountsMale[1][1][2][0][0] << "\t\t\t"
						      << bestchi3[k].controlcountsMale[1][1][0][0][0] << "\t" << bestchi3[k].controlcountsMale[1][1][1][0][0] << "\t"
						      << bestchi3[k].controlcountsMale[1][1][2][0][0] << "\n" << bestchi3[k].casecountsMale[1][2][0][0][0] << "\t"
						      << bestchi3[k].casecountsMale[1][2][1][0][0] << "\t" << bestchi3[k].casecountsMale[1][2][2][0][0] << "\t\t\t"
						      << bestchi3[k].controlcountsMale[1][2][0][0][0] << "\t" << bestchi3[k].controlcountsMale[1][2][1][0][0] << "\t"
						      << bestchi3[k].controlcountsMale[1][2][2][0][0] << "\n" << bestchi3[k].casecountsMale[2][0][0][0][0] << "\t"
						      << bestchi3[k].casecountsMale[2][0][1][0][0] << "\t" << bestchi3[k].casecountsMale[2][0][2][0][0] << "\t\t\t"
						      << bestchi3[k].controlcountsMale[2][0][0][0][0] << "\t" << bestchi3[k].controlcountsMale[2][0][1][0][0] << "\t"
						      << bestchi3[k].controlcountsMale[2][0][2][0][0] << "\n" << bestchi3[k].casecountsMale[2][1][0][0][0] << "\t"
						      << bestchi3[k].casecountsMale[2][1][1][0][0] << "\t" << bestchi3[k].casecountsMale[2][1][2][0][0] << "\t\t\t"
						      << bestchi3[k].controlcountsMale[2][1][0][0][0] << "\t" << bestchi3[k].controlcountsMale[2][1][1][0][0] << "\t"
						      << bestchi3[k].controlcountsMale[2][1][2][0][0] << "\n" << bestchi3[k].casecountsMale[2][2][0][0][0] << "\t"
						      << bestchi3[k].casecountsMale[2][2][1][0][0] << "\t" << bestchi3[k].casecountsMale[2][2][2][0][0] << "\t\t\t"
						      << bestchi3[k].controlcountsMale[2][2][0][0][0] << "\t" << bestchi3[k].controlcountsMale[2][2][1][0][0] << "\t"
						      << bestchi3[k].controlcountsMale[2][2][2][0][0] << "\n\n";


			      // female
			      bestMarkerCombi3Details << "female: \n";
			      bestMarkerCombi3Details << "cases:\t\t\t\tcontrols:\t\t\n";
			      bestMarkerCombi3Details << bestchi3[k].casecountsFemale[0][0][0][0][0] << "\t" << bestchi3[k].casecountsFemale[0][0][1][0][0] << "\t" << bestchi3[k].casecountsFemale[0][0][2][0][0] << "\t\t\t" << bestchi3[k].controlcountsFemale[0][0][0][0][0] << "\t" << bestchi3[k].controlcountsFemale[0][0][1][0][0] << "\t" << bestchi3[k].controlcountsFemale[0][0][2][0][0] << "\n" << bestchi3[k].casecountsFemale[0][1][0][0][0] << "\t" << bestchi3[k].casecountsFemale[0][1][1][0][0] << "\t" << bestchi3[k].casecountsFemale[0][1][2][0][0] << "\t\t\t" << bestchi3[k].controlcountsFemale[0][1][0][0][0] << "\t" << bestchi3[k].controlcountsFemale[0][1][1][0][0] << "\t" << bestchi3[k].controlcountsFemale[0][1][2][0][0] << "\n" << bestchi3[k].casecountsFemale[0][2][0][0][0] << "\t" << bestchi3[k].casecountsFemale[0][2][1][0][0] << "\t" << bestchi3[k].casecountsFemale[0][2][2][0][0] << "\t\t\t" << bestchi3[k].controlcountsFemale[0][2][0][0][0] << "\t" << bestchi3[k].controlcountsFemale[0][2][1][0][0] << "\t" << bestchi3[k].controlcountsFemale[0][2][2][0][0] << "\n" << bestchi3[k].casecountsFemale[1][0][0][0][0] << "\t" << bestchi3[k].casecountsFemale[1][0][1][0][0] << "\t" << bestchi3[k].casecountsFemale[1][0][2][0][0] << "\t\t\t" << bestchi3[k].controlcountsFemale[1][0][0][0][0] << "\t" << bestchi3[k].controlcountsMale[1][0][1][0][0] << "\t" << bestchi3[k].controlcountsMale[1][0][2][0][0] << "\n" << bestchi3[k].casecountsFemale[1][1][0][0][0] << "\t" << bestchi3[k].casecountsFemale[1][1][1][0][0] << "\t" << bestchi3[k].casecountsFemale[1][1][2][0][0] << "\t\t\t" << bestchi3[k].controlcountsFemale[1][1][0][0][0] << "\t" << bestchi3[k].controlcountsFemale[1][1][1][0][0] << "\t" << bestchi3[k].controlcountsFemale[1][1][2][0][0] << "\n" << bestchi3[k].casecountsFemale[1][2][0][0][0] << "\t" << bestchi3[k].casecountsFemale[1][2][1][0][0] << "\t" << bestchi3[k].casecountsFemale[1][2][2][0][0] << "\t\t\t" << bestchi3[k].controlcountsFemale[1][2][0][0][0] << "\t" << bestchi3[k].controlcountsFemale[1][2][1][0][0] << "\t" << bestchi3[k].controlcountsMale[1][2][2][0][0] << "\n" << bestchi3[k].casecountsFemale[2][0][0][0][0] << "\t" << bestchi3[k].casecountsFemale[2][0][1][0][0] << "\t" << bestchi3[k].casecountsFemale[2][0][2][0][0] << "\t\t\t" << bestchi3[k].controlcountsFemale[2][0][0][0][0] << "\t" << bestchi3[k].controlcountsFemale[2][0][1][0][0] << "\t" << bestchi3[k].controlcountsFemale[2][0][2][0][0] << "\n" << bestchi3[k].casecountsFemale[2][1][0][0][0] << "\t" << bestchi3[k].casecountsFemale[2][1][1][0][0] << "\t" << bestchi3[k].casecountsFemale[2][1][2][0][0] << "\t\t\t" << bestchi3[k].controlcountsFemale[2][1][0][0][0] << "\t" << bestchi3[k].controlcountsFemale[2][1][1][0][0] << "\t" << bestchi3[k].controlcountsFemale[2][1][2][0][0] << "\n" << bestchi3[k].casecountsFemale[2][2][0][0][0] << "\t" << bestchi3[k].casecountsFemale[2][2][1][0][0] << "\t" << bestchi3[k].casecountsFemale[2][2][2][0][0] << "\t\t\t" << bestchi3[k].controlcountsFemale[2][2][0][0][0] << "\t" << bestchi3[k].controlcountsFemale[2][2][1][0][0] << "\t" << bestchi3[k].controlcountsFemale[2][2][2][0][0] << "\n\n";

			      if (test > 2)
				{
				  bestMarkerCombi3Details << "\nmale: L1:\n";
				  if(printBeta) resultReg(bestchi3[k].result1, bestMarkerCombi3Details,haplo, qt,maxIndexCov);
				  //bestMarkerCombi3Details << "\nmale: L2:\n";
				  //resultReg(bestchi3[k].result2M, bestMarkerCombi3Details,haplo, qt);
				  bestMarkerCombi3Details << "\n";

				  //bestMarkerCombi3Details << "\nfemale: L1:\n";
				  //if(printBeta) resultReg(bestchi3[k].result1F, bestMarkerCombi3Details,haplo, qt,maxIndexCov);
				  //bestMarkerCombi3Details << "\nfemale: L2:\n";
				  //resultReg(bestchi3[k].result2F, bestMarkerCombi3Details,haplo, qt);
				  bestMarkerCombi3Details << "\n\n";
				}
			    }
			}
		    }
		  bestMarkerCombi3 << "\n\n Note: Each combination is listed only once, even if it occurs multiple times in the combilist\n e.g.: (SNP1, SNP2, SNP3) = (SNP2, SNP1, SNP3)\n";
		  bestMarkerCombi3.close();
		  if (qt == 1)
		    {
		      bestMarkerCombi3Details << "\n\n Note: The output for linear regression is still under construction. Parameter estimates will be provided in the next intersnp version.\n";
		    }

		  bestMarkerCombi3Details.close();
		}
	    }
	  if(plln)
	    {
	      free(bestsinglemarker);bestsinglemarker=NULL;

	      if (genetictop > 0)

		{
		  free(bestsinglemarker2);bestsinglemarker2=NULL;
		  free(geneticlist2);geneticlist2=NULL;
		}

	      if(multimarker && bestchi3)
		{
		  for (int z=0;z<printtop;z++)
		    {
		      freeSTATplus(bestchi3[z].result1);freeSTATplus(bestchi3[z].result3);
		      //freeSTATplus(bestchi3[z].result1M);freeSTATplus(bestchi3[z].result2M);
		      //freeSTATplus(bestchi3[z].result1F);freeSTATplus(bestchi3[z].result2F);
		    }
		  free(bestchi3);bestchi3=NULL;
		}
	    }
	  if(dohapfile && n==0)
	    {
	      fclose(fptr10);
	    }

	  if (BinSNPsCCFlags != NULL) delete_2dim(BinSNPsCCFlags);

	  if(plln)
	    {
	      if(combilist)
		{
		  free(snp1);free(snp2);free(snp3);
		}
	      else
		{
		  if(nsnps>=1){free(snp1);}
		  if(nsnps>=2){free(snp2);}
		  if(nsnps>=3){free(snp3);}
		}
	    }
	} // abortnsim

      } // end n

#if PARALLELN


  } //end pragma n
#endif


  if (nsim>0) {
    cout << "\rn = " << nsim << " of " << nsim << endl;

  }




  //PAA
  if (pathwayAnalysis == 1)
    {

      for(i=0;i<nlist+1;i++)
	{
	  Ttable[i][nsim+1]=1;
	}
      for(i=0;i<nsim+2;i++)
	{
	  Ttable[nlist][i]=0;
	}

      //jede zeile sortieren und in p-wert umwandeln
      for(i=0;i<nlist+1;i++)
	{
	  nbetter=0;
	  for(j=0;j<nsim+2;j++)
	    {
	      hTable[0][j]=j; //index
	      hTable[1][j]=Ttable[i][j]; //value of current row
	    }
	  qsortplus(hTable[1],hTable[0],0,nsim,-1); //want to order hTable[1]!

	  for(j=0;j<=nsim;j++)
	    {
	      for(k=j+1;k<=nsim;k++)
		{
		  if(hTable[1][k]<hTable[1][j])
		    {
		      break;
		    }
		}
	      nbetter=k-1;
	      if(nbetter>nsim){nbetter=nsim;}

	      for(kk=j;kk<k;kk++)
		{
		  hTable[1][kk]=(double)(nbetter/(double)nsim);
		}
	      j=k-1;
	    }

	  hTable[1][nsim]=1;

	  //sort back
	  qsortplus(hTable[0],hTable[1],0,nsim,1); //want to order according to hTable[0]!

	  //einfügen in ursprüngliche tabelle
	  for(j=0;j<nsim+1;j++)
	    {
	      Ttable[i][j]=hTable[1][j];
	    }
	  // cout << nbetter << " " << Ttable[i][0] << "\n";
	} // i

      //berechnen von bestp per simulation
      for(k=0;k<=nsim;k++)
	{
	  Ttable[nlist][k]=1;
	  for(j=0;j<nlist;j++)
	    {
	      if (Ttable[j][k]<Ttable[nlist][k])
		{
		  Ttable[nlist][k]=Ttable[j][k];
		}
	    }
	}
      //berechnen von korrp per pathway
      for(i=0;i<nlist;i++)
	{
	  Ttable[i][nsim+1]=0;
	  for(j=1;j<=nsim;j++)
	    {
	      if(Ttable[nlist][j]<=Ttable[i][0])
		{
		  Ttable[i][nsim+1]++;
		}
	    }
	  Ttable[i][nsim+1]=Ttable[i][nsim+1]/((double)(nsim));
	}
    }	//end PAA


  if (nsim > 0 && !bin)  //ausgabe PWT anpassen
    {
      cout << "\n\n";


      /// calculate corrected MC P-values
      for (i = 0; i < printtop; i++) {
	for(k=1; k<=nsim; k++) {
	  if (bestofSim[k] < toplist[i].p*(1+EPS)) countsP[i] += 1;
	}
	countsP[i] = countsP[i] / ((double)nsim);
      }

      /// estimation of MC inflation factor
      if (!pathwayAnalysis && !pathwayImpact && !markercombi2 && !markercombi3) {
	double pVal, inflationfactorMC = 0;
	for (int c=0,i=0; i < singletop; i++) {
	  if (countsMC[i]) pVal = countsMC[i]/(double)nsim;
	  else {
	    pVal = (++c-0.5)/(double)nsnpAnalysisqc;
	    if (pVal > 1./nsim) pVal = exp(-1.0)/nsim;
	  }
	  inflationfactorMC += chiinv(1, pVal)/nsnpAnalysisqc;
	}
	sstm << "Monte-Carlo inflation factor: " << inflationfactorMC << "\n";
	logg(sstm);
      }

      toplistMC.open(toplistMCfile.c_str(), ios::out);

      for (i = 0; i < printtop; i++) // i läuft durch die Pathways
	{
	  if (toplist[i].nr1 != -1 && nsim > 0)
	    {

	      if (i == 0)
		{
		  //PWT
		  if (pathwayAnalysis == 1)
		    {
		      if (pathwayTest == 1)
			{
			  toplistMC << "Pathway\tMC-p-value\tp-corr\tGENEcounts\tSNPcounts\tSNPratio\n";
			  cout << "Pathway MC-p-value p-corr GENEcounts SNPcounts SNPratio\n";
			  logfile << "Pathway\tMC-p-value\tp-corr\tGENEcounts\tSNPcounts\tSNPratio\n";
			}
		      if (pathwayTest == 2)
			{
			  toplistMC << "Pathway\tMC-p-value\tp-corr\tGENEcounts\tSNPcounts\tSNPratio\tFisherStat\n";
			  cout << "Pathway MC-p-value p-corr GENEcounts SNPcounts SNPratio FisherStat\n";
			  logfile << "Pathway\tMC-p-value\tGENEcounts\tSNPcounts\tSNPratio\tFisherStat\n";
			}
		      if (pathwayTest == 3)
			{
			  toplistMC << "Pathway\tMC-p-value\tp-corr\tGENEcounts\tSNPcounts\tGENEratio\n";
			  cout << "Pathway MC-p-value p-corr GENEcounts SNPcounts GENEratio\n";
			  logfile << "Pathway\tMC-p-value\tp-corr\tGENEcounts\tSNPcounts\tGENEratio\n";
			}
		      if (pathwayTest == 4)
			{
			  toplistMC << "Pathway\tMC-p-value\tp-corr\tGENEcounts\tSNPcounts\tGENEratio\tFisherMaxStat\n";
			  cout << "Pathway MC-p-value p-corr GENEcounts SNPcounts GENEratio FisherMaxStat\n";
			  logfile << "Pathway\tMC-p-value\tp-corr\tGENEcounts\tSNPcounts\tGENEratio\tFisherMaxStat\n";
			}
		      if (pathwayTest == 5)
			{
			  toplistMC << "Pathway\tMC-p-value\tp-corr\tGENEcounts\tSNPcounts\tGENEratio\tFisherMaxStatPlusn";
			  cout << "Pathway MC-p-value p-corr GENEcounts SNPcounts GENEratio FisherMaxStatPlus\n";


			  logfile << "Pathway\tMC-p-value\tp-corr\tGENEcounts\tSNPcounts\tGENEratio\tFisherMaxStatPlus\n";
			}
		      if (pathwayTest == 6)
			{
			  toplistMC << "Pathway\tMC-p-value\tp-corr\tGENEcounts\tSNPcounts\tInteractionRatio\n";
			  cout << "Pathway MC-p-value p-corr GENEcounts SNPcounts InteractionRatio\n";
			  logfile << "Pathway\tMC-p-value\tp-corr\tGENEcounts\tSNPcounts\tInteractionRatio\n";
			}

		    }
		  else if (pathwayImpact == 1)
		    {
		      toplistMC << "No\tChr_No_1\trs_No_1\tChr_No_2\trs_No_2\tChr_No_3\trs_No_3\tPos_No_1\tPos_No_2\tPos_No_3\tp-value\tMC-p-value_Corr\tPathwy\n";
		      cout << "No Chr_No_1 rs_No_1 Chr_No_2 rs_No_2 Chr_No_3 rs_No_3 Pos_No_1 Pos_No_2 Pos_No_3 p-value MC-p-value_Corr Pathway\n";
		      logfile << "No Chr_No_1 rs_No_1 Chr_No_2 rs_No_2 Chr_No_3 rs_No_3 Pos_No_1 Pos_No_2 Pos_No_3 p-value MC-p-value_Corr Pathway\n";
		    }
		  else if(!markercombi2 && !markercombi3)
		    {
		      toplistMC << "No\tChr_No_1\trs_No_1\tChr_No_2\trs_No_2\tChr_No_3\trs_No_3\tPos_No_1\tPos_No_2\tPos_No_3\tP-value\tMC-P-value\tMC-P-value_Corr\n";
		      logfile   << "No Chr_No_1 rs_No_1 Chr_No_2 rs_No_2 Chr_No_3 rs_No_3 Pos_No_1 Pos_No_2 Pos_No_3 P-value MC-P-value MC-P-value_Corr\n";
		      cout      << "No Chr_No_1 rs_No_1 Chr_No_2 rs_No_2 Chr_No_3 rs_No_3 Pos_No_1 Pos_No_2 Pos_No_3 P-value MC-P-value MC-P-value_Corr\n";
		    }
		  else
		    {
		      toplistMC << "No\tChr_No_1\trs_No_1\tChr_No_2\trs_No_2\tChr_No_3\trs_No_3\tPos_No_1\tPos_No_2\tPos_No_3\tP-value\tMC-P-value_Corr\n";
		      logfile   << "No Chr_No_1 rs_No_1 Chr_No_2 rs_No_2 Chr_No_3 rs_No_3 Pos_No_1 Pos_No_2 Pos_No_3 P-value MC-P-value_Corr\n";
		      cout      << "No Chr_No_1 rs_No_1 Chr_No_2 rs_No_2 Chr_No_3 rs_No_3 Pos_No_1 Pos_No_2 Pos_No_3 P-value MC-P-value_Corr\n";
		    }
		}

	      //PWT
	      if (pathwayAnalysis == 1)
		{
		  if (pathwayTest == 1 || pathwayTest == 3 || pathwayTest == 6)
		    {
		      toplistMC << pathway[i].name << "\t" <<  Ttable[i][0] << "\t" << Ttable[i][nsim+1] << "\t" << pathway[i].ngenes << "\t" << pathway[i].counts0 <<  "\t" <<  pathway[i].ratio0 << "\n";
		      if (i < 10)
			{
			  cout << pathway[i].name << " " <<  Ttable[i][0] << " " << Ttable[i][nsim+1] << " " << pathway[i].ngenes << " " << pathway[i].counts0 <<  " " <<  pathway[i].ratio0 << "\n";
			  logfile << pathway[i].name << "\t" <<  Ttable[i][0] << "\t" << Ttable[i][nsim+1] << "\t" << pathway[i].ngenes << "\t" << pathway[i].counts0 <<  "\t"  <<  pathway[i].ratio0 << "\n";
			}
		    }
		  else if (pathwayTest == 2 || pathwayTest == 4 || pathwayTest == 5 )
		    {
		      toplistMC << pathway[i].name << "\t" <<  Ttable[i][0] << "\t" << Ttable[i][nsim+1] << "\t" << pathway[i].ngenes << "\t" << pathway[i].counts0 <<  "\t" <<  pathway[i].ratio0 << "\t" <<  pathway[i].score0 << "\n";
		      if (i < 10)
			{
			  cout << pathway[i].name << " " <<  Ttable[i][0] << " " << Ttable[i][nsim+1] << " " << pathway[i].ngenes << " " << pathway[i].counts0 <<  " " <<  pathway[i].ratio0 << " " <<  pathway[i].score0 << "\n";
			  logfile << pathway[i].name << "\t" <<  Ttable[i][0] << "\t" << Ttable[i][nsim+1] << "\t" << pathway[i].ngenes << "\t" << pathway[i].counts0 <<  "\t" <<  pathway[i].ratio0 << "\t" <<  pathway[i].score0<< "\n";
			}
		    }
		}
	      else if (pathwayImpact == 1)
		{
		  if (toplist[i].nr2 != -1 && toplist[i].nr3 != -1)
		    {
		      toplistMC << i + 1 << "\t" << map[toplist[i].nr1].chr << "\t" << map[toplist[i].nr1].rs << "\t"
				<< map[toplist[i].nr2].chr << "\t" << map[toplist[i].nr2].rs << "\t" << map[toplist[i].nr3].chr
				<< "\t" << map[toplist[i].nr3].rs << "\t" << map[toplist[i].nr1].pos << "\t" << map[toplist[i].nr2].pos
				<< "\t" << map[toplist[i].nr3].pos << "\t" << toplist[i].p << "\t" << countsP[i]  << "\t"
				<< pathway[toplist[i].r-1].name << "\n";

		      if (i < 10)
			{
			  cout << i + 1 << " " << map[toplist[i].nr1].chr << " " << map[toplist[i].nr1].rs << " "
			       << map[toplist[i].nr2].chr << " " << map[toplist[i].nr2].rs << " " << map[toplist[i].nr3].chr
			       << " " << map[toplist[i].nr3].rs << " " << map[toplist[i].nr1].pos << " " << map[toplist[i].nr2].pos
			       << " " << map[toplist[i].nr3].pos + 1 << " " << toplist[i].p << " " << countsP[i]  << " "
			       << pathway[toplist[i].r-1].name << "\n";

			  logfile << i + 1 << " " << map[toplist[i].nr1].chr << " " << map[toplist[i].nr1].rs << " "
				  << map[toplist[i].nr2].chr << " " << map[toplist[i].nr2].rs << " " << map[toplist[i].nr3].chr
				  << " " << map[toplist[i].nr3].rs << " " << map[toplist[i].nr1].pos << " " << map[toplist[i].nr2].pos
				  << " " << map[toplist[i].nr3].pos << " " << toplist[i].p << " " << countsP[i]  << " "
				  << pathway[toplist[i].r-1].name << "\n";

			}
		    }
		  else if (toplist[i].nr2 != -1)
		    {
		      toplistMC << i + 1 << "\t" << map[toplist[i].nr1].chr << "\t" << map[toplist[i].nr1].rs << "\t"
				<< map[toplist[i].nr2].chr << "\t" << map[toplist[i].nr2].rs << "\t" << "-"
				<< "\t" << "-" << "\t" << map[toplist[i].nr1].pos << "\t" << map[toplist[i].nr2].pos
				<< "\t" << "-" << "\t" << toplist[i].p << "\t" << countsP[i]  << "\t" << pathway[toplist[i].r-1].name << "\n";

		      if (i < 10)
			{
			  cout << i + 1 << " " << map[toplist[i].nr1].chr << " " << map[toplist[i].nr1].rs << " "
			       << map[toplist[i].nr2].chr << " " << map[toplist[i].nr2].rs << " " << "-"
			       << " " << "-" << " " << map[toplist[i].nr1].pos << " " << map[toplist[i].nr2].pos
			       << " " << "-" << " " << toplist[i].p << " " << countsP[i]  << " " << pathway[toplist[i].r-1].name << "\n";

			  logfile << i + 1 << " " << map[toplist[i].nr1].chr << " " << map[toplist[i].nr1].rs << " "
				  << map[toplist[i].nr2].chr << " " << map[toplist[i].nr2].rs << " " << "-"
				  << " " << "-" << " " << map[toplist[i].nr1].pos << " " << map[toplist[i].nr2].pos
				  << " " << "-" << " " << toplist[i].p << " " << countsP[i]  << " " << pathway[toplist[i].r-1].name << "\n";
			}
		    }
		  else
		    {
		      toplistMC << i + 1 << "\t" << map[toplist[i].nr1].chr << "\t" << map[toplist[i].nr1].rs << "\t"
				<< "-" << "\t" << "-" << "\t" << "-" << "\t" << "-" << "\t" << map[toplist[i].nr1].pos << "\t" << "-"
				<< "\t" << "-" << "\t" << toplist[i].p << "\t" << countsP[i]  << "\t" << pathway[toplist[i].r-1].name << "\n";

		      if (i < 10)
			{
			  cout << i + 1 << " " << map[toplist[i].nr1].chr << " " << map[toplist[i].nr1].rs << " "
			       << "-" << " " << "-" << " " << "-" << " " << "-" << " " << map[toplist[i].nr1].pos << " " << "-"
			       << " " << "-" << " " << toplist[i].p << " " << countsP[i] << " " << pathway[toplist[i].r-1].name << "\n";

			  logfile << i + 1 << " " << map[toplist[i].nr1].chr << " " << map[toplist[i].nr1].rs << " "
				  << "-" << " " << "-" << " " << "-" << " " << "-" << " " << map[toplist[i].nr1].pos << " " << "-"
				  << " " << "-" << " " << toplist[i].p << " " << countsP[i] << " " << pathway[toplist[i].r-1].name << "\n";
			}
		    }
		}
	      else
		{
		  if (toplist[i].nr2 != -1 && toplist[i].nr3 != -1)
		    {
		      toplistMC << i + 1 << "\t" << map[toplist[i].nr1].chr << "\t" << map[toplist[i].nr1].rs << "\t"
				<< map[toplist[i].nr2].chr << "\t" << map[toplist[i].nr2].rs << "\t" << map[toplist[i].nr3].chr
				<< "\t" << map[toplist[i].nr3].rs << "\t" << map[toplist[i].nr1].pos << "\t" << map[toplist[i].nr2].pos
				<< "\t" << map[toplist[i].nr3].pos << "\t" << toplist[i].p << "\t" << countsP[i]  << "\n";

		      if (i < 10)
			{
			  cout << i + 1 << " " << map[toplist[i].nr1].chr << " " << map[toplist[i].nr1].rs << " "

			       << map[toplist[i].nr2].chr << " " << map[toplist[i].nr2].rs << " " << map[toplist[i].nr3].chr
			       << " " << map[toplist[i].nr3].rs << " " << map[toplist[i].nr1].pos << " " << map[toplist[i].nr2].pos
			       << " " << map[toplist[i].nr3].pos << " " << toplist[i].p << " " << countsP[i]  << "\n";

			  logfile << i + 1 << " " << map[toplist[i].nr1].chr << " " << map[toplist[i].nr1].rs << " "
				  << map[toplist[i].nr2].chr << " " << map[toplist[i].nr2].rs << " " << map[toplist[i].nr3].chr
				  << " " << map[toplist[i].nr3].rs << " " << map[toplist[i].nr1].pos << " " << map[toplist[i].nr2].pos
				  << " " << map[toplist[i].nr3].pos << " " << toplist[i].p << " " << countsP[i] << "\n";
			}
		    }
		  else if (toplist[i].nr2 != -1)
		    {
		      toplistMC << i + 1 << "\t" << map[toplist[i].nr1].chr << "\t" << map[toplist[i].nr1].rs << "\t"
				<< map[toplist[i].nr2].chr << "\t" << map[toplist[i].nr2].rs << "\t" << "-"
				<< "\t" << "-" << "\t" << map[toplist[i].nr1].pos << "\t" << map[toplist[i].nr2].pos
				<< "\t" << "-" << "\t" << toplist[i].p << "\t" << countsP[i]  << "\n";

		      if (i < 10)
			{
			  cout << i + 1 << " " << map[toplist[i].nr1].chr << " " << map[toplist[i].nr1].rs << " "
			       << map[toplist[i].nr2].chr << " " << map[toplist[i].nr2].rs << " " << "-"
			       << " " << "-" << " " << map[toplist[i].nr1].pos << " " << map[toplist[i].nr2].pos
			       << " " << "-" << " " << toplist[i].p << " " << countsP[i]  << "\n";

			  logfile << i + 1 << " " << map[toplist[i].nr1].chr << " " << map[toplist[i].nr1].rs << " "
				  << map[toplist[i].nr2].chr << " " << map[toplist[i].nr2].rs << " " << "-"
				  << " " << "-" << " " << map[toplist[i].nr1].pos << " " << map[toplist[i].nr2].pos
				  << " " << "-" << " " << toplist[i].p << " " << countsP[i]  << "\n";
			}
		    }
		  else  /// single-marker ToplistMC output
		    {
		      if (group_test) {
			toplistMC << i + 1 << "\t" << map[toplist[i].nr1].chr << "\t" << map[toplist[i].nr1].rs << "\t"
				  << "-" << "\t" << "-" << "\t" << "-" << "\t" << "-" << "\t" << map[toplist[i].nr1].pos << "\t" << "-"
				  << "\t" << "-" << "\t" << "-" << "\t" << countsMC[i]/(double)nsim << "\t" << countsP[i]  << "\n";
			if (i < 10) {
			  cout << i + 1 << " " << map[toplist[i].nr1].chr << " " << map[toplist[i].nr1].rs << " "
			       << "-" << " " << "-" << " " << "-" << " " << "-" << " " << map[toplist[i].nr1].pos << " " << "-"
			       << " " << "-" << " " << "-" << " " << countsMC[i]/(double)nsim << " " << countsP[i] << "\n";
			  logfile << i + 1 << " " << map[toplist[i].nr1].chr << " " << map[toplist[i].nr1].rs << " "
				  << "-" << " " << "-" << " " << "-" << " " << "-" << " " << map[toplist[i].nr1].pos << " " << "-"
				  << " " << "-" << " " << "-" << " " << countsMC[i]/(double)nsim << " " << countsP[i] << "\n";
			}
		      } else {
			toplistMC << i + 1 << "\t" << map[toplist[i].nr1].chr << "\t" << map[toplist[i].nr1].rs << "\t"
				  << "-" << "\t" << "-" << "\t" << "-" << "\t" << "-" << "\t" << map[toplist[i].nr1].pos << "\t" << "-"
				  << "\t" << "-" << "\t" << toplist[i].p << "\t" << countsMC[i]/(double)nsim << "\t" << countsP[i]  << "\n";
			if (i < 10) {
			  cout << i + 1 << " " << map[toplist[i].nr1].chr << " " << map[toplist[i].nr1].rs << " "
			       << "-" << " " << "-" << " " << "-" << " " << "-" << " " << map[toplist[i].nr1].pos << " " << "-"
			       << " " << "-" << " " << toplist[i].p << " " << countsMC[i]/(double)nsim << " " << countsP[i] << "\n";
			  logfile << i + 1 << " " << map[toplist[i].nr1].chr << " " << map[toplist[i].nr1].rs << " "
				  << "-" << " " << "-" << " " << "-" << " " << "-" << " " << map[toplist[i].nr1].pos << " " << "-"
				  << " " << "-" << " " << toplist[i].p << " " << countsMC[i]/(double)nsim << " " << countsP[i] << "\n";
			}
		      }
		    }
		}
	    }
	}
      toplistMC.close();



    } // nsim > 0


#if RARE
  if(bin){
    if(collinter==0 && vb==0){
      outRare(binsizeRare, nwindows, nwindowssinglebin, windowPositions, person, map, counts, optimalrare, FISHERflag, REGRESSIONflag, FRACREGflag,COLLREGflag, COLLflag, CMATflag, nRareSNPsCOLL, nRareSNPsCMAT, nRareSNPsFISHER, nRareSNPsREGRESSION, nRareSNPsFRACREG,nRareSNPsCOLLREG,  MAF_Level_VT_FISHER, MAF_Level_VT_CMAT, MAF_Level_VT_COLL, MAF_Level_VT_REGRESSION, MAF_Level_VT_FRACREG, MAF_Level_VT_COLLREG, nsim, rarefile,rareTopfile, raretop, intervalfile, outputname, errorfile, logfile,  rarepretest, rarepretestlimit, rareregpretest, rarefFISHER, rarefREGRESSION, rarefFRACREG,rarefCOLLREG, rarefCOLL, OR_COLL_vec, OR_COLL_f_vec, rarefCMAT, OR_CMAT_vec, pCOLLvec, pCMATvec, pFISHERvec, pREGRESSIONvec, pFRACREGvec, pCOLLREGvec, pFISHERvecChi, nSNPsInWindow, raref, thread_nloop, FISHERnotconverged,FISHERpretestpassed, REGRESSIONpretestpassed, FRACREGpretestpassed, COLLREGpretestpassed, CMATpretestpassed, COLLpretestpassed, window, FISHER_CI, CMAT_CI, COLL_CI, REGRESSION_CI, FRACREG_CI, COLLREG_CI, COLLREG_beta, COLLREG_se, FRACREG_beta, FRACREG_se, FISHERcount, REGRESSIONcount, FRACREGcount,COLLREGcount, COLLcount, CMATcount, merging, flanking, catIntervals, minRareInBin, expandIntervals, mafadjust, nlinestped, binamin, binamax, weights,verbose, SetIDfile, intervaleditor,setid);
    }
    else if((collinter==1 || collinter==2 || collinter==3 || collinter== 4) && vb==0){
      outRareInter(binsizeRare, nwindows, nwindowssinglebin, windowPositions, person, map, counts, optimalrare, FISHERflag, REGRESSIONflag, FRACREGflag,COLLREGflag, COLLflag, CMATflag, nRareSNPsCOLL, nRareSNPsCMAT, nRareSNPsFISHER, nRareSNPsREGRESSION, nRareSNPsFRACREG,nRareSNPsCOLLREG,  MAF_Level_VT_FISHER, MAF_Level_VT_CMAT, MAF_Level_VT_COLL, MAF_Level_VT_REGRESSION, MAF_Level_VT_FRACREG, MAF_Level_VT_COLLREG, nsim, rarefileinter,rareTopfile, raretop, intervalfile, outputname, errorfile, logfile,  rarepretest, rarepretestlimit, rareregpretest, rarefFISHER, rarefREGRESSION, rarefFRACREG,rarefCOLLREG, rarefCOLL, OR_COLL_vec, OR_COLL_f_vec, rarefCMAT, OR_CMAT_vec, pCOLLvec, pCMATvec, pFISHERvec, pREGRESSIONvec, pFRACREGvec, pCOLLREGvec, pFISHERvecChi, nSNPsInWindow, raref, thread_nloop, FISHERnotconverged,FISHERpretestpassed, REGRESSIONpretestpassed, FRACREGpretestpassed, COLLREGpretestpassed, CMATpretestpassed, COLLpretestpassed, window, FISHER_CI, CMAT_CI, COLL_CI, REGRESSION_CI, FRACREG_CI, COLLREG_CI, COLLREG_beta, COLLREG_se, FRACREG_beta, FRACREG_se, FISHERcount, REGRESSIONcount, FRACREGcount,COLLREGcount, COLLcount, CMATcount, merging, flanking, catIntervals, minRareInBin, expandIntervals, mafadjust, nlinestped, binamin, binamax, weights,verbose, SetIDfile, intervaleditor,setid, collinter);
    }
    else if(vb==1){
      outRareVB(person, map, counts, optimalrare, nsim, rarefileVB, rarefileVBgraph, rarefileVBpermstat, outputname, errorfile, logfile,  rarepretest, rarepretestlimit, rareregpretest, nSNPsInWindow, raref, thread_nloop, window, minRareInBin,  mafadjust, nlinestped, verbose, SetIDfile, intervaleditor,setid, nvbstartvt, vbstartvt, nvbend, vbend, vbmaxpermstat, nvblevel, nchunks, chunkpos, chunklen, nwindows, nwordsSNPs, SNPMapInverse, BinSNPs, BinSNPsCCFlagsOriginal, ncasesqc, ncontrolsqc, ndummyatlevel, dummypos, dummylevel, BinCarriers, dummyend, ndummyends,  nchunkcluster, dummycluster, vb_binwise_corr, vb_print_perm_stat, vbbinwisestat, vbbinwisecount);
      // delete
      if(vb_binwise_corr){
	for(int l=0; l<nwindows; l++){
	  for(int m=0; m<nvbstartvt[l]; m++){
	    for(int m2=0; m2<nvblevel[l][m]; m2++){
	      for(int m3=0; m3<nchunks[l][m][m2]; m3++){
		delete[] vbbinwisestat[l][m][m2][m3];
		delete[] vbbinwisecount[l][m][m2][m3];
	      }
	      delete[] vbbinwisestat[l][m][m2];
	      delete[] vbbinwisecount[l][m][m2];
	    }
	    delete[] vbbinwisestat[l][m];
	    delete[] vbbinwisecount[l][m];
	  }
	  delete[] vbbinwisestat[l];
	  delete[] vbbinwisecount[l];
	}
	delete[] vbbinwisestat;
	delete[] vbbinwisecount;
      }

      delete[] nCarriers;
      for(int i=0; i<nwordsSNPs; i++){
	delete[] BinSNPsCCFlagsOriginal[i];
      }
      delete[] BinSNPsCCFlagsOriginal;
      delete[] vbmaxpermstat;
      for(int i=0; i<nwindows; i++){
	for(int j=0; j<nvbstartvt[i]; j++){
	  for(int k=0; k<nvblevel[i][j]; k++){
	    delete[] chunkpos[i][j][k];
	    delete[] chunklen[i][j][k];
	  }
	  delete[] BinCarriers[i][j];
	  delete[] chunkpos[i][j];
	  delete[] chunklen[i][j];
	  delete[] nchunks[i][j];
	}
	delete[] chunkpos[i];
	delete[] chunklen[i];
	delete[] nchunks[i];
	if(optimalrare){
	  delete[] rareLimitsNCT[i];
	  delete[] rareLimitsNCTInverse[i];
	}
	free(vbstartvt[i]);
	delete[] BinCarriers[i];
      }
      delete[] BinCarriers;
      if(optimalrare){
	delete[] rareLimitsNCT;
	delete[] rareLimitsNCTInverse;
      }
      delete[] chunkpos;
      delete[] chunklen;
      delete[] nchunks;
      delete[] vbstartvt;

      for(int i=0; i<nwindows; i++){
	for(int j=0; j<nvbstartvt[i]; j++){
	  for(int k=0; k<nvblevel[i][j]; k++){
	    for(int l=0; l<nchunkcluster[i][j][k]; l++){
	      //      for(int m=0; m<ndummyatlevel[i][j][k][l]; m++){
	      //    for(int l=0; l<ndummyends[i][j][k]+1; l++){
	      delete[] dummylevel[i][j][k][l];
	      delete[] dummypos[i][j][k][l];
	    }
	    delete[] ndummyatlevel[i][j][k];
	    delete[] dummyend[i][j][k];
	    delete[] dummylevel[i][j][k];
	    delete[] dummypos[i][j][k];
	    delete[] dummycluster[i][j][k];
	    //  }
	  }
	  delete[] ndummyatlevel[i][j];
	  delete[] dummyend[i][j];
	  delete[] ndummyends[i][j];
	  delete[] nchunkcluster[i][j];
	  delete[] dummylevel[i][j];
	  delete[] dummypos[i][j];
	  delete[] dummycluster[i][j];
	}
	delete[] ndummyatlevel[i];
	delete[] dummyend[i];
	delete[] ndummyends[i];
	delete[] nchunkcluster[i];
	delete[] dummylevel[i];
	delete[] dummypos[i];
	delete[] dummycluster[i];
	delete[] nvblevel[i];
      }
      delete[] ndummyends;
      delete[] dummyend;
      delete[] nchunkcluster;
      delete[] nvbstartvt;
      delete[] dummylevel;
      delete[] dummypos;
      delete[] dummycluster;
      delete[] ndummyatlevel;
      delete[] nvblevel;
      //      delete[] window;
    }
    else if(vb==0){
      if(FISHERflag){
	delete[] rarefFISHER;
	delete[] pFISHERvec;
	delete[] pFISHERvecChi;
	delete[] FISHERcount;
	delete[] completedFISHER;
	delete[] FISHERstats;
	delete[] FISHERpermstats;
	delete[] FISHERnotconverged;
	if(optimalrare){
	  delete[] MAF_Level_VT_FISHER;
	}
	if(wilsonpretest!=0 || rareregpretest!=0){
	  delete[] FISHERpretestpassed;
	}
	if(nsim!=0){
	  free2Ddouble(FISHER_CI,nwindows);
	}
      }
      if(REGRESSIONflag){
	delete[] rarefREGRESSION;
	delete[] pREGRESSIONvec;
	delete[] REGRESSIONcount;
	delete[] completedREGRESSION;
	delete[] REGRESSIONstats;
	delete[] REGRESSIONpermstats;
	if(optimalrare){
	  delete[] MAF_Level_VT_REGRESSION;
	}
	if(wilsonpretest!=0 || rareregpretest!=0){
	  delete[] REGRESSIONpretestpassed;
	}
	if(nsim!=0){
	  free2Ddouble(REGRESSION_CI,nwindows);
	}
      }
      if(FRACREGflag){
	delete[] rarefFRACREG;
	delete[] pFRACREGvec;
	delete[] FRACREGcount;
	delete[] completedFRACREG;
	delete[] FRACREGstats;
	delete[] FRACREGpermstats;
	delete[] FRACREG_beta;
	delete[] FRACREG_se;
	if(optimalrare){
	  delete[] MAF_Level_VT_FRACREG;
	}
	if(wilsonpretest!=0 || rareregpretest!=0){
	  delete[] FRACREGpretestpassed;
	}
	if(nsim!=0){
	  free2Ddouble(FRACREG_CI,nwindows);
	}
      }
      if(COLLREGflag){
	delete[] rarefCOLLREG;
	delete[] pCOLLREGvec;
	if(collinter!=3 && collinter!=4){
	  delete[] COLLREGcount;
	  delete[] completedCOLLREG;
	  delete[] COLLREGstats;
	  delete[] COLLREGpermstats;
	  delete[] COLLREG_beta;
	  delete[] COLLREG_se;
	}
	if(optimalrare){
	  delete[] MAF_Level_VT_COLLREG;
	}
	if(wilsonpretest!=0 || rareregpretest!=0){
	  delete[] COLLREGpretestpassed;
	}
	if(nsim!=0){
	  free2Ddouble(COLLREG_CI,nwindows);
	}
      }
      if(CMATflag){
	delete[] rarefCMAT;
	delete[] OR_CMAT_vec;
	delete[] pCMATvec;
	delete[] CMATcount;
	delete[] completedCMAT;
	delete[] CMATstats;
	delete[] CMATpermstats;
	if(optimalrare){
	  delete[] MAF_Level_VT_CMAT;
	}
	if(wilsonpretest!=0 || rareregpretest!=0){
	  delete[] CMATpretestpassed;
	}
	if(nsim!=0){
	  free2Ddouble(CMAT_CI,nwindows);
	}
      }
      if(COLLflag){
	delete[] rarefCOLL;
	delete[] OR_COLL_vec;
	delete[] OR_COLL_f_vec;
	delete[] pCOLLvec;
	delete[] COLLcount;
	delete[] completedCOLL;
	delete[] COLLstats;
	delete[] COLLpermstats;
	if(optimalrare){
	  delete[] MAF_Level_VT_COLL;
	}
	if(wilsonpretest!=0 || rareregpretest!=0){
	  delete[] COLLpretestpassed;
	}
	if(nsim!=0){
	  free2Ddouble(COLL_CI,nwindows);
	}
      }

      if(!vb && COLLflag){
	for(l=0; l<nwindows; l++){
	  delete[] limitsStatCOLL[l];
	}
	delete[] limitsStatCOLL;

      }
      if(CMATflag){
	for(l=0; l<nwindows; l++){
	  delete[] limitsStatCMAT[l];
	}
	delete[] limitsStatCMAT;

      }
      if(FISHERflag){
	for(l=0; l<nwindows; l++){
	  delete[] limitsStatFISHER[l];
	}
	delete[] limitsStatFISHER;;
      }
      if(REGRESSIONflag){
	for(l=0; l<nwindows; l++){
	  //	  delete[] limitsStatREGRESSION[l];
	}
	//	delete[] limitsStatREGRESSION;
      }
      if(FRACREGflag){
	for(l=0; l<nwindows; l++){
	  //	  delete[] limitsStatFRACREG[l];
	}
	//	delete[] limitsStatFRACREG;
      }
      if(COLLREGflag){
	for(l=0; l<nwindows; l++){
	  //	  delete[] limitsStatCOLLREG[l];
	}
	//	delete[] limitsStatCOLLREG;
      }
    }
    if(optimalrare && vb==0){
      for(j=0; j<nwindows; j++){
	delete[] rareLimits[j];
      }
      delete[] rareLimits;
    }
    if(collinter!=3 && collinter!=4){
      for(j=0; j<nwindows; j++){
	for(int k=0; k<nRareLimits[j]; k++){
	  delete[] window[j].levelpos[k];
	}
	delete[] window[j].levelpos;
	delete[] window[j].n_at_level;
      }
      delete[] window;
    }
    if(verbose==2){
      delete[] nSNPsInWindow;
    }
    delete[] nRareSNPs;
    delete[] nRareLimits;
    free2Dint(windowPositions,nwindowssinglebin);
    free2Dint(nchrwindows,26);
    if (rare_stratify) {
      for (uint32_t w=0; w<nMWindowsRARE; w++) {
        for (int i=0; i<MatchingRARE[w].nGroups; i++) delete[] MatchingRARE[w].groups[i].list;
        delete[] MatchingRARE[w].groups;
        if (MatchingRARE[w].BinSNPsCCFlags!=NULL) delete_2dim(MatchingRARE[w].BinSNPsCCFlags);
      }
      delete[] MatchingRARE; MatchingRARE=NULL;
    }
    free(intervals);
  }
#endif

  /***** free heap *****/
  if (BinSNPs            != NULL) delete_3dim(BinSNPs, nsnpqc);
  if (BinPPLs            != NULL) delete_3dim(BinPPLs, npplqc);
  if (BinSNPsGenderFlags != NULL) delete_2dim(BinSNPsGenderFlags);
  if (BinSNPsCCFlags     != NULL) delete_2dim(BinSNPsCCFlags);
  if (SNPMap        !=NULL) { delete[] SNPMap;        SNPMap       =NULL; }
  if (PPLMap        !=NULL) { delete[] PPLMap;        PPLMap       =NULL; }
  if (SNPMapInverse !=NULL) { delete[] SNPMapInverse; SNPMapInverse=NULL; }
  if (PPLMapInverse !=NULL) { delete[] PPLMapInverse; PPLMapInverse=NULL; }
  if (PplLocations  !=NULL) { delete[] PplLocations;  PplLocations =NULL; }
  if (PplLocFam     !=NULL) { delete[] PplLocFam;     PplLocFam    =NULL; }
  if (chrPositions  !=NULL) { delete[] chrPositions;  chrPositions =NULL; }
  if (MatchedPairs  !=NULL) { delete[] MatchedPairs;  MatchedPairs =NULL; }
  if (Clusters      !=NULL) {
    for (int i=0; i<nClusters; i++) { delete[] Clusters[i].list; Clusters[i].list=NULL; }
    delete[] Clusters; Clusters=NULL;
  }
  if(dosage){
    for(int i=0;i<nlinestped;i++){
      for(int j=0;j<nlinestfam;j++){
	free(genoWeights[i][j]);
      }
      free(genoWeights[i]);
    }
    free(genoWeights);
  }
  for (int i=0; i<nlinestfam; i++) {
    free(person[i].fid);
    free(person[i].pid);
    free(person[i].vid);
    free(person[i].mid);
  }
  free(person);
  for (int i=0; i<nlinestped; i++) {
    free(map[i].rs);
    free(map[i].gene);
  }
  free(map);


  for (int i=0; i<maxthreads; i++) {
    if(!qt) {free(Y[i]);} //provisional fix IS691
    for (int j=0; j<nlinestped; j++) free(counts[i][j].det);
    free(counts[i]);
  }
  //exit(1);
  free(Y);

  free(counts);
  for(int i = 0; i<nlinestped;i++){
    free(codesA[i].a1); codesA[i].a1 = NULL;
    free(codesA[i].a2); codesA[i].a2 = NULL;
  }
  free(codesA);
  free(toplist);
  free(hapfile);
  free(bestofSim);
  free(countsMC);
  free(countsP);
  free(bestsinglemarker);
  free(SigmaAuxVector);
  free(SigmaAuxVector_Single);


  errorfile.close();
  logfile.close();
  return 0;
} // end main

