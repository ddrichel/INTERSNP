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

int MatrixMultSym(double **A, int Adim1, int Adim2, double **B, int Bdim1, int Bdim2, double **C);
int DwyerInv(int n, double **MM, double **D, double **T, double **U, double **Ut, double **sumPP, double **sumPJ, double **sumPK, double **MMinv, double **Minv);
int mysvd(int n, int m, double ***U, double ***V, double **S, int la, int lv);
int MatrixMultMod(double **A, int Adim1, int Adim2, double **B, int Bdim1, int Bdim2, double **C);
double tquantile(int n);

struct STATplus regGeneral(int x[27], //in this function just a dummy
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
			   int numberOfAllCov, int ncov,
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
			   int firstbinlastSNP
			   )
{
  int covIndex=0;
  double logNew;
  double betaNew[nx+ncov];
  double betaOld[nx+ncov];
  int i,k,j,jj,f=0;
  int complete=0;
  int dummy;
  int ncase=0;
  double rss=0;
  double ssy=0;
  double Yavg =0;
  double n1=0;
  double f1=0;
  double se_beta[nx+ncov];
  double YtY = 0;
  double YtXAinvXtY = 0;
  double newf=-1;
  int shift=0;
  int n=0;


  result.sc=0;

  if(xtype>=3 && !male && !female)
    {
      if(!sexcov){sexcov=1;nx++;}
    }
  else if(xtype>=3)
    {
      if(sexcov){sexcov=0;nx--;}
    }

  if(nx>dim1)
    {
      cout << "\nWarning! Model cannot be computed because matrix dimension ("<<dim1<<") is too small.\n";
      exit(1);
      //realloc matrices:
      result.df=0;return result;
    }

  if(collapseRare==1 or (collcollapseRare==1 && nInters==0))
    {
      if(nSnpsDom>0 || nEnvirons>0)
	{
	  cout << "Error from regGeneral. FRACREG and COLLREG are implemented only for additive snp-paramters.\n";
	  exit(1);
	}
    }

  //determine df
  f=nx+ncov;
  k=0;

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
	  if(xtype==0)
	    {
	      if(person[k].qcin==1)
		{
		  if(qt)
		    {
		      Y[n]=person[k].qtaff[thread];
		    }
		  if(!qt)
		    {
		      Y[n]=person[k].aff[thread];
		    }
		  Yhelp[k]=1;
		}
	      else {Yhelp[k]=-1;continue;}
	    }
	  else if(xtype>=3)
	    {
	      if(person[k].sex==0 ){Yhelp[k]=-1;continue;}
	      if(person[k].qcin==1)
		{
		  if(qt)
		    {
		      Y[n]=person[k].qtaff[thread];
		    }
		  if(!qt)
		    {
		      Y[n]=person[k].aff[thread];
		    }
		  Yhelp[k]=1;
		}
	      else {Yhelp[k]=-1;continue;}
	    }
	}
      else //! alt=1
	{
	  if(qt)
	    {
	      Y[n]=person[k].qtaff[thread];
	    }
	  if(!qt)
	    {
	      Y[n]=person[k].aff[thread];
	    }
	}

      //intercept
      X[n][0]=1;

      //additve snp parameters
      for (int h=0;h<nSnps;h++)
	{
	  if(getbit64(BinSNPs[snps[h]][guy.nr][1],guy.pos))
	    {
	      X[n][h+1]=1*weightvec[snps[h]];
	    }
	  else if(getbit64(BinSNPs[snps[h]][guy.nr][2],guy.pos))
	    {
	      X[n][h+1]=0;
	    }
	  else if(getbit64(BinSNPs[snps[h]][guy.nr][3],guy.pos))
	    {
	      X[n][h+1]=-1*weightvec[snps[h]];
	    }
	  else
	    {
	      Yhelp[k]=-1;
	      goto nextInd;
	    }
	}

      //dominant snp parameters
      shift=nSnps;
      for (int h=0;h<nSnpsDom;h++)
	{
	  if(getbit64(BinSNPs[snpsDom[h]][guy.nr][1],guy.pos))
	    {
	      X[n][h+shift+1]=-0.5*weightvec[snps[h]];
	    }
	  else if(getbit64(BinSNPs[snpsDom[h]][guy.nr][2],guy.pos))
	    {
	      X[n][h+shift+1]=0.5*weightvec[snps[h]];
	    }
	  else if(getbit64(BinSNPs[snpsDom[h]][guy.nr][3],guy.pos))
	    {
	      X[n][h+shift+1]=-0.5*weightvec[snps[h]];
	    }
	  else
	    {
	      Yhelp[k]=-1;
	      goto nextInd;
	    }
	}

      //non-genetic parameters, under concstruction
      shift+=nSnpsDom;
      for (int h=0;h<nEnvirons;h++)
	{
	  if(1)
	    {
	      X[n][h+shift+1]=1;
	    }
	  else
	    {
	      Yhelp[k]=-1;
	      goto nextInd;
	    }
	}
      //interaction parameters, under construction
      //shift+=nEnvirons;
      for (int h=0;h<nInters;h++)
	{
	  if(1)
	    {
	      X[n][h+shift+1]=1;
	    }
	  else
	    {
	      Yhelp[k]=-1;
	      goto nextInd;
	    }
	}

      //covariates
      //shift+=nInters;

      covIndex=0;
      for (int h=0;h<numberOfAllCov;h++){
	  if(cov[h]){
	    if(person[k].covin[h]==0){
		Yhelp[k]=-1; complete=0;goto nextInd;
	      }
	    else{
		X[n][covIndex+shift+1]=person[k].cov[h];
	      }
	    covIndex++;
	  }
	}
      if(sexcov)
	    {
	      if(person[k].sex==1){X[n][f]=1;}
	      else if(person[k].sex==2){X[n][f]=0;}
	      else {Yhelp[k]=-1;complete=0;goto nextInd;}
	    }

	  if(complete)
	    {
	      n++;
	      ncase++;
	    }
	nextInd:;

	} //end k-loop (individuals)

      if((collapseRare==1 || collcollapseRare==1) && alt==1 && nInters==0)
	{
	  if(collapseRare==1)
	    {
	      double thisX=0;
	      double weightsum=0;
	      for (int h=0;h<nSnps;h++)
		{
		  weightsum=weightsum+weightvec[snps[h]];
		}
	      for (int i=0;i<n;i++)
		{
		  thisX=0;
		  for (int h=0;h<nSnps;h++)
		    {
		      if(X[i][h+1]>=0){
			thisX+=weightvec[snps[h]];
		      }
		    }
		  thisX=thisX/(double(weightsum));

		  X[i][1]=thisX;
		}
	    }
	  if(collcollapseRare==1)
	    {
	      double thisX=0;
	      for (int i=0;i<n;i++)
		{
		  thisX=0;
		  for (int h=0;h<nSnps;h++)
		    {
		      if(X[i][h+1]>=0){
			thisX=1; // at least one rare allele
			break;
		      }
		    }
		  X[i][1]=thisX;
		}
	    }
	  //shift back covariates
	  for (int i=0;i<n;i++)
	    {
	      for (int h=0;h<ncov;h++)
		{
		  X[i][h+2]=X[i][h+shift+1];
		}
	    }
	  if(sexcov)
	    {
	      for (int i=0;i<n;i++)
		{
		  X[i][2+ncov]=X[i][f];
		}
	    }
	  f=f-nSnps+1;
	} //end matrix modification collapseRare

      //      for(int i=0; i<=f; i++)
      //      	{
      //	       		cout << "postmod: X_Init[10][" << i << "]=" << X[10][i] << endl;
      //      	}


      //	  for(int i=0; i<=f; i++){
	    //   	    	    cout<<"postmod X[n-1]["<<i<<"] "<<X[nlinestfam-1][i]<<endl;
	    //	  }
	  //	  cout<<n <<" " <<nlinestfam<<endl;

      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      /*
	cout << "nInters=" << nInters << endl;

	if(collcollapseRare==1 && alt==1 && nInters==1 ) // COLL_INTER 3
	{
	double thisX_Bin1=0;
	double thisX_Bin2=0;

	for (int i=0;i<n;i++)
	{
	thisX_Bin1=0;
	thisX_Bin2=0;
	for (int h=0;h<firstbinlastSNP;h++)
	{
	if(X[i][h+1]>=0)
	{
	thisX_Bin1=1; // at least one rare allele
	break;
	}
	}
	X[i][1]=thisX_Bin1;
	//cout << "X[" << i << "][1]" << X[i][1] << endl;

	for (int h=firstbinlastSNP;h<nSnps;h++)
	{
	if(X[i][h+1]>=0)
	{
	thisX_Bin2=1; // at least one rare allele
	break;
	}
	}
	X[i][2]=thisX_Bin2;

	}

	//non-genetic parameters, under concstruction
	//shift+=nSnpsDom;
	for (int i=0;i<n;i++)
	{
	for (int h=0;h<nEnvirons;h++)
	{
	X[i][3+h]=1;
	}
	}



	//interaction parameter, under construction

	for (int i=0;i<n;i++)
	{
	X[i][3+nEnvirons]=(X[i][1])*(X[i][2]);
	}


	//covariates

	for (int i=0;i<n;i++)
	{
	for (int h=0;h<numberOfAllCov;h++)
	{
	X[i][3+nEnvirons+1+h]=X[i][h+1+nSnps+nSnpsDom+nEnvirons+nInters];
	}
	}

	//sexcov

	for (int i=0;i<n;i++)
	{
	X[i][3+nEnvirons+numberOfAllCov+1]=X[i][f];
	}

	cout << "nEnvirons=" << nEnvirons <<endl;
	cout << "numberOfAllCov=" << numberOfAllCov << endl;
	cout << "sexcov=" << sexcov << endl;

	f=f-nSnps-nSnpsDom+2;


	for(int i=0; i<4+nEnvirons+numberOfAllCov+sexcov;i++)
	{
	cout << "X_new[10][" << i << "]=" << X[10][i] << endl;

	}

	for(int i=0; i<4+nEnvirons+numberOfAllCov+sexcov;i++)
	{
	cout << "X_new[1000][" << i << "]=" << X[1000][i] << endl;
	}


	} //end X matrix modification for COLL_INTER 3



	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      */
      for(i=0;i<nx+ncov;i++){betaNew[i]=-1;se_beta[i]=-1;}
      for(j=0;j<=f;j++)
	{
	  betaNew[j]=0;se_beta[j]=0;
	}

      n1=n;
      f1=f;

      for(i=0;i<n;i++)
	{
	  Yavg+=Y[i];
	}

      if(n>0){Yavg=Yavg/n;}

      if(f>0 && ncase>0)
	{
	  for(j=0;j<=f;j++)
	    {
	      betaOld[j]=betaNew[j];
	    }

	  for (i=0;i<n;i++)
	    {
	      //exponent=0;
	      for(j=0;j<=f;j++)
		{
		  Xt[j][i]=X[i][j];
		}
	    }

	  dummy=MatrixMultSym(Xt,f+1,n,X,n,f+1,A);

	  // Invertieren von A
	  // first try Dwyer-Algorithmus
	  //	dummy=DwyerInv(f+1, A, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, Ainv);
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
		}
	      dummy=MatrixMultMod(VNN, f+1, f+1,Sinv, f+1, f+1, A0);
	      dummy=MatrixMultMod(A0, f+1, f+1,UNNT, f+1, f+1, Ainv);
	    }
	  else
	    {
	      newf=f;
	    }

	  //rss = Y'Y - Y'XAinvX'Y berechnen
	  for (i=0;i<n;i++)
	    {
	      Yt[0][i]=Y[i];
	    }

	  for (i=0;i<n;i++)
	    {
	      YtY+=Yt[0][i]*Y[i];
	    }
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
	  if(nMc==0)
	    {
	      result.betaNew_se[1] = -1;result.betaNew_lcl[1] = -1;result.betaNew_rcl[1] = -1;
	    }
	  result.b[1] = 0;
	}
      for(i=0;i<f+1;i++)
	{
	  betaNew[i]=0;
	  if(newf>0){betaNew[i]=YtXAinv[0][i];}
	  if(nMc==0)
	    {
	      result.betaNew_se[i] = -1;result.betaNew_lcl[i] = -1;result.betaNew_rcl[i] = -1;
	    }
	  result.b[i] = -1;
	}

      if(newf>0)
	{
	  for(i=0;i<f+1;i++)
	    {
	      if (n1-newf-1 > 0 && f > 0 && nMc==0)
		{
		  if((Ainv[i][i]*rss/(n1-newf-1))<0)
		    {
		      result.betaNew_se[i] = -1;result.betaNew_lcl[i] = -1;result.betaNew_rcl[i] = -1;
		    }
		  else
		    {
		      result.betaNew_se[i]=sqrt(Ainv[i][i]*rss/(n1-newf-1));
		      result.betaNew_lcl[i]= betaNew[i]-result.betaNew_se[i]*tquantile((int)(n1-newf-1));
		      result.betaNew_rcl[i]= betaNew[i]+result.betaNew_se[i]*tquantile((int)(n1-newf-1));
		    }
		}
	    }

	}
      else
	{
	  betaNew[1]=0;

	  if(nMc==0)
	    {
	      result.betaNew_se[1] = -1;result.betaNew_lcl[1] = -1;result.betaNew_rcl[1] = -1;
	    }
	  result.b[1] = -1;
	}

      result.df=n1-newf-1;
      for(i=0;i<f+1;i++)
	{
	  result.b[i]=betaNew[i];

	}
      result.bcvsex=betaNew[f];  //covariate sex

      for(i=0;i<ncov;i++)
	{
	  result.bcv[i]=betaNew[i+nx]; //covariate parameters
	}
      result.sc=rss;
      result.rsquare=1-(rss/ssy); // calculate rsquare

      return result;
    };


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
			      int numberOfAllCov, int ncov,
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
			      )
{
  int covIndex=0;
  double logNew;
  double betaNew[nx+ncov];
  double betaOld[nx+ncov];
  int i,k,j,jj,f=0;
  int complete=0;
  int dummy;
  int ncase=0;
  int ncontrol=0;
  double rss=0;
  double ssy=0;
  double Yavg =0;
  double n1=0;
  double f1=0;
  double se_beta[nx+ncov];
  double YtY = 0;
  double YtXAinvXtY = 0;
  double newf=-1;
  int shift=0;
  int n=0;


  result.sc=0;

  if(xtype>=3 && !male && !female)
    {
      if(!sexcov){sexcov=1;nx++;}
    }
  else if(xtype>=3)
    {
      if(sexcov){sexcov=0;nx--;}
    }

  if(nx>dim1)
    {
      cout << "\nWarning! Model cannot be computed because matrix dimension ("<<dim1<<") is too small.\n";
      exit(1);
      //realloc matrices:
      result.df=0;return result;
    }

  if(collapseRare==1 or (collcollapseRare==1 && nInters==0))
    {
      if(nSnpsDom>0 || nEnvirons>0)
	{
	  cout << "Error from logRegGeneral. FRACREG and COLLREG are implemented only for additive snp-paramters.\n";
	  exit(1);
	}
    }

  //determine df
  f=nx+ncov;
  k=0;

  //get n, Y and X
  n=0;ncase=0;ncontrol=0;

  struct PPLLOCATION guy;

  for(int kMod=0;kMod<npplqc;kMod++)
    {
      guy = PplLocations[kMod];
      k = PPLMap[kMod];
      complete=1;

      if(alt==0 && Yhelp[k]< -0.1){continue;}

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


	if(alt==1)
	{
	  //intercept
      X[n][0]=1;

      //additve snp parameters
      for (int h=0;h<nSnps;h++)
	{
	  if(getbit64(BinSNPs[snps[h]][guy.nr][1],guy.pos))
	    {
	      X[n][h+1]=1*weightvec[snps[h]];
	    }
	  else if(getbit64(BinSNPs[snps[h]][guy.nr][2],guy.pos))
	    {
	      X[n][h+1]=0;
	    }
	  else if(getbit64(BinSNPs[snps[h]][guy.nr][3],guy.pos))
	    {
	      X[n][h+1]=-1*weightvec[snps[h]];
	    }
	  else
	    {
	      Yhelp[k]=-1;
	      goto nextInd;
	    }
	}

      //dominant snp parameters
      shift=nSnps;
      for (int h=0;h<nSnpsDom;h++)
	{
	  if(getbit64(BinSNPs[snpsDom[h]][guy.nr][1],guy.pos))
	    {
	      X[n][h+shift+1]=-0.5*weightvec[snps[h]];
	    }
	  else if(getbit64(BinSNPs[snpsDom[h]][guy.nr][2],guy.pos))
	    {
	      X[n][h+shift+1]=0.5*weightvec[snps[h]];
	    }
	  else if(getbit64(BinSNPs[snpsDom[h]][guy.nr][3],guy.pos))
	    {
	      X[n][h+shift+1]=-0.5*weightvec[snps[h]];
	    }
	  else
	    {
	      Yhelp[k]=-1;
	      goto nextInd;
	    }
	}

      //non-genetic parameters, under construction
      shift+=nSnpsDom;
      for (int h=0;h<nEnvirons;h++)
	{
	  if(1)
	    {
	      X[n][h+shift+1]=1;
	    }
	  else
	    {
	      Yhelp[k]=-1;
	      goto nextInd;
	    }
	}
      //interaction parameters, under construction
      shift+=nEnvirons;
      for (int h=0;h<nInters;h++)
	{
	  if(nInters==1)
	    {
	      X[n][h+shift+1]=1;
	    }
	  else
	    {
	      Yhelp[k]=-1;
	      goto nextInd;
	    }
	}

      //covariates
      shift+=nInters;

		covIndex=0;
		for(int h=0;h<numberOfAllCov;h++){
		if(cov[h])
			{
				if(person[k].covin[h]==0)
					{
						Yhelp[k]=-1; complete=0;goto nextInd;
					}
				else
					{
						X[n][covIndex+shift+1]=person[k].cov[h];
					}
		covIndex++;
			}
    }

      //sexcov
      if(sexcov)
	{
	  if(person[k].sex==1){X[n][f+nInters]=1;}
	  else if(person[k].sex==2){X[n][f+nInters]=0;}
	  else {Yhelp[k]=-1;complete=0;goto nextInd;}
	}
    }
      if(complete)
	{
	  n++;
	  //ncase++;
	  if(person[k].aff[thread]==2){ncase++;}
	  else if(person[k].aff[thread]==1){ncontrol++;}
	}
    nextInd:;
    } //end k-loop (individuals)


  if((collapseRare==1 || collcollapseRare==1) && alt==1 && nInters==0) // matrix modification collapseRare
    {
      if(collapseRare==1)
	{
	  double thisX=0;
	  double weightsum=0;
	  for (int h=0;h<nSnps;h++)
	    {
	      weightsum=weightsum+weightvec[snps[h]];
	    }
	  for (int i=0;i<n;i++)
	    {
	      thisX=0;
	      for (int h=0;h<nSnps;h++)
		{
		  if(X[i][h+1]>=0)
		    {
		      thisX+=weightvec[snps[h]];
		    }
		}
	      thisX=thisX/(double(weightsum));

	      X[i][1]=thisX;
	    }
	}
      if(collcollapseRare==1)
	{
	  double thisX=0;
	  for (int i=0;i<n;i++)
	    {
	      thisX=0;
	      for (int h=0;h<nSnps;h++)
		{
		  if(X[i][h+1]>=0){
		    thisX=1; // at least one rare allele
		    break;
		  }
		}
	      X[i][1]=thisX;
	    }
	}
      //shift back covariates
      for (int i=0;i<n;i++)
	{
	  for (int h=0;h<ncov;h++)
	    {
	      X[i][h+2]=X[i][h+shift+1];
	    }
	}
      if(sexcov)
	{
	  for (int i=0;i<n;i++)
	    {
	      X[i][2+ncov]=X[i][f];
	    }
	}
      f=f-nSnps+1;
    } //end matrix modification collapseRare

	/*
	cout << "collinter=" << collinter << endl;
	cout << "alt=" << alt << endl;
	cout << "nSnps=" << nSnps << endl;
	cout << "nSnpsDom=" << nSnpsDom << endl;
	cout << "nEnvirons=" << nEnvirons << endl;
	cout << "nInters=" << nInters << endl;
	cout << "nx=" << nx << endl;
	cout << "sexcov=" << sexcov << endl;
	cout << "ncov=" << ncov << endl;

	cout << "X_initial" << endl;
	for(int t=0;t<10;t++)
			{
				for(int j=0; j<2+nSnps+nSnpsDom+nEnvirons+ncov+sexcov+10;j++)
					{
						cout << X[t][j] << "\t";
					}
				cout << endl;
			}


	ofstream file3;
	file3.open("Matrix_initial.txt");

	for(int t=0;t<n;t++)
			{
				for(int j=0; j<10;j++)
					{

						file3 << X[t][j] << "\t";

					}
				file3 << endl;
			}
	file3.clear();
	file3.close();
	*/

	if((collinter==3 || collinter==4) && alt==1) // COLL_INTER 3 and 4, alt==1
	{

		double thisX_Bin1=0;
		double thisX_Bin2=0;

		for (int i=0;i<n;i++)
			{
				thisX_Bin1=0;
				thisX_Bin2=0;



				for (int h=0;h<=firstbinlastSNP;h++)
					{
						if(X[i][h+1]>=0)
							{
								thisX_Bin1=1; // at least one rare allele
								break;
							}
					}

				for (int h=firstbinlastSNP+1;h<nSnps;h++)
					{
						 if(X[i][h+1]>=0)
							{
								thisX_Bin2=1; // at least one rare allele
								break;
							}
					}

				X[i][1]=thisX_Bin1;
				X[i][2]=thisX_Bin2;
			}

		//non-genetic parameters, under concstruction

		for (int i=0;i<n;i++)
			{
				for (int h=0;h<nEnvirons;h++)
					{
						X[i][3+h]=1;
					}

			}

		//interaction parameter

		for (int i=0;i<n;i++)
			{
				X[i][3+nEnvirons]=(X[i][1])*(X[i][2]);
			}


		//covariates

		for (int i=0;i<n;i++)
			{
				for (int h=0;h<ncov;h++)
					{
						X[i][3+nEnvirons+1+h]=X[i][h+1+nSnps+nSnpsDom+nEnvirons+nInters];
					}
			}

		//sexcov

		if(sexcov)
			{
				for (int i=0;i<n;i++)
					{
						//X[i][3+nEnvirons+ncov+1]=X[i][f+sexcov];
						X[i][3+nEnvirons+ncov+1]=X[i][nSnps+nSnpsDom+nEnvirons+nInters+ncov+sexcov];
					}
			}

		f=f-nSnps-nSnpsDom+2;
		//cout << "f=" << f << endl;

		/*
		cout << "X_new" << endl;
		for(int tt=0;tt<10;tt++) // 10 first individuals
			{

				for(int j=0; j<2+nEnvirons+nInters+ncov+sexcov+1;j++)
				{
					cout << X[tt][j] << "\t";
				}
				cout << endl;
			}


		ofstream file4;
		file4.open("Matrix_new_alt1.txt");

		for(int tt=0;tt<n;tt++)
			{

				for(int j=0; j<2+nEnvirons+nInters+ncov+sexcov+1;j++)
				{
					file4 << X[tt][j] << "\t";
				}
				file4 << endl;
			}
		file4.clear();
		file4.close();
		*/

    } //end X matrix modification for COLL_INTER 3 and 4, alt==1

	if(collinter==3 && alt==0) // COLL_INTER 3, alt==0
	{

		for(int i=0;i<n;i++)
			{

				for(int h=0;h<nEnvirons;h++)
					{
						X[i][h+1]=X[i][h+3];
					}


				for (int h=0;h<ncov;h++)
					{
						X[i][1+nEnvirons+h]=X[i][4+nEnvirons+h];
					}

				if(sexcov==1)
					{
						X[i][1+nEnvirons+ncov]=X[i][3+nEnvirons+ncov+sexcov];
					}

			}

		f=f-nSnps-nSnpsDom-nInters;
		//cout << "f=" << f << endl;

		/*
		cout << "X_new" << endl;
		for(int ttt=0;ttt<10;ttt++)
			{

				for(int j=0; j<1+nEnvirons+ncov+sexcov;j++)
					{
						cout << X[ttt][j] << "\t";

					}
				cout << endl;
			}
		*/
	}


	else if(collinter==4 && alt==0) // COLL_INTER 4, alt==0
	{
		for(int i=0;i<n;i++)
			{

				for (int h=0;h<ncov;h++)
					{
						X[i][3+nEnvirons+h]=X[i][4+nEnvirons+h];
					}

				if(sexcov==1)
					{
						X[i][3+nEnvirons+ncov]=X[i][3+nEnvirons+ncov+sexcov];
					}

			}

		f=f-nSnps-nSnpsDom+2;
		//cout << "f=" << f << endl;

		/*
		cout << "X_new" << endl;
		for(int ttt=0;ttt<10;ttt++)
			{

				for(int j=0; j<3+nEnvirons+ncov+sexcov;j++)
					{
						cout << X[ttt][j] << "\t";
					}
				cout << endl;
			}


		ofstream file5;
		file5.open("Matrix_new_alt0.txt");

		for(int ttt=0;ttt<n;ttt++)
			{

				for(int j=0; j<3+nEnvirons+ncov+sexcov;j++)
					{
						file5 << X[ttt][j] << "\t";
					}
				file5 << endl;
			}
		file5.clear();
		file5.close();
		*/
	}

	else if(collinter==0 && collapseRare==0 && collcollapseRare==0 && alt==0) // REG test
	{
		for(int i=0;i<n;i++)
			{

				for (int h=0;h<ncov;h++)
					{
						X[i][1+nEnvirons+h]=X[i][1+firstbinlastSNP+nEnvirons+h];
					}
				if(sexcov==1)
					{
						X[i][1+nEnvirons+ncov]=X[i][1+firstbinlastSNP+nEnvirons+ncov+sexcov];
					}

			}

		f=f-nSnps-nSnpsDom;
	}

	else if(collinter==0 && (collapseRare==1 || collcollapseRare==1) && alt==0) // COLLREG, FRACREG test
	{
		for(int i=0;i<n;i++)
			{

				for (int h=0;h<ncov;h++)
					{
						X[i][1+nEnvirons+h]=X[i][2+nEnvirons+h];
					}
				if(sexcov==1)
					{
						X[i][1+nEnvirons+ncov]=X[i][2+nEnvirons+ncov+sexcov];
					}

			}

		f=f-nSnps-nSnpsDom;
	}




	for(i=0;i<nx+ncov;i++){betaNew[i]=-1;se_beta[i]=-1;}
	for(j=0;j<=f;j++)
		{
			betaNew[j]=0;se_beta[j]=0;
		}

	n1=n;
	f1=f;

	//cout << "f1=" << f1 << endl;
	/*

	cout<<"alt "<<alt<<endl;
	for (i=0;i<1;i++)
	  {
	    cout<<person[i].pid<< " ";


	    for(j=0;j<=f;j++)
	      {
		cout<<X[i][j]<<" ";
	      }
	    cout<<endl;
	  }

	*/




  //MAXIMIZATION aus funktion logReg in intersnp.cpp
  logNew=0;
  double logOld=1;
  double exponent=0;
  int dfUnchanged=1;
  int maxit=10;
  int it=0;
  double eps=0.000001; //convergence criterion

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
	      betaOld[j]=betaNew[j];
	    }
	  //determine p

	  if(it==1)
	    {
	      for (i=0;i<n;i++)
		{
		  exponent=0;
		  for(j=0;j<=f;j++)
		    {
		      exponent+=betaOld[j]*X[i][j];
		      Xt[j][i]=X[i][j];
		    }

		  p[i]=1/(1+exp(-exponent));
		}
	    }

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
	  //if(maxIndexCov==0){dummy=DwyerInv(f+1, A, D, T, U, Ut, sumPP, sumPJ, sumPK, MMinv, Ainv);}
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
		  //		  cout<<" S["<<i<<"] "<<S[i]<<endl;
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
	      betaNew[i]=betaOld[i]+newbeta[i][0];
	    }

	  //compute logNew
	  logNew=0;
	  for(i=0;i<n;i++)
	    {
	      exponent=0;

	      for(j=0;j<=f;j++)
		{
		  exponent+=betaNew[j]*X[i][j];
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
  //cout << "df=" << newf << endl;

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
	  result.b[i]=betaNew[i];  // betas

	  if(nMc==0 && (ncase+ncontrol-newf-1) > 0)
	    {
	      if((Ainv[i][i]*rss/(ncase+ncontrol-newf-1))<0)
		{
		  result.betaNew_se[i] = -1;
		  result.lcloddsRatio[i] = -1;
		  result.rcloddsRatio[i] = -1;

		}
	      else
		{
		  //result.betaNew_se[number2indi[i]] = 2*sqrt(Ainv[i][i]*rss/(ncase+ncontrol-newf-1));
		  result.betaNew_se[i] = sqrt(Ainv[i][i]); //changed on 09.12.2013
		  result.lcloddsRatio[i] = exp(betaNew[i] -1.96*result.betaNew_se[i]);
		  result.rcloddsRatio[i] = exp(betaNew[i] +1.96*result.betaNew_se[i]);
		}
	      result.oddsRatio[i]=exp(betaNew[i]); // odds ratios
	    }
	}
    }
  else
    {
      for(i=0;i< f+1;i++)
	{
	  result.b[i]=-2;  // betas
	  result.betaNew_se[i] = -2;
	  result.lcloddsRatio[i] = -2;
	  result.rcloddsRatio[i] = -2;
	  result.oddsRatio[i]= -2; // odds ratios
	}
    }
  /*result.bcvsex=betaNew[27];  //covariate sex
    for(i=0;i<maxIndexCov;i++)
    {
    result.bcv[i]=betaNew[i+27]; //covariate parameters
    }*/
  result.sc=logNew;

  //cout << "alt=" << alt << "result.sc=" << result.sc << endl;


  if(alt==1 && covariancematrix)  //SIGMA1
    {
      int k=0;

      for(i=0;i<f+1;i++)
	{
	  for(j=i;j<f+1;j++)
	    {
	      if(dfUnchanged==1) {result.sigma1[k]=Ainv[i][j];} // sigma1 matrix as 1 dim array
	      else {result.sigma1[k]=-2;}
	      k=k+1;
	    }
	}
    }
  return result;
};
