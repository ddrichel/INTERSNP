/**~ isnp_files.cpp
.---------------------------------------------------------------------------.
|  Software: INTERSNP, Genome-Wide Interaction Analysis                     |
|      Site: http://intersnp.meb.uni-bonn.de/                               |
| ------------------------------------------------------------------------- |
|      File: isnp_out.cpp                                                   |
|    Author:                                                                |
|   Content: source-code library on reading PLINK file-formats for INTERSNP |
| Copyright (c) 2014-,                                                      |
| ------------------------------------------------------------------------- |
|   License: Distributed under the General Public License (GPL)             |
|            http://www.gnu.org/licenses/                                   |
| This library is distributed in the hope that it will be useful - WITHOUT  |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or     |
| FITNESS FOR A PARTICULAR PURPOSE.                                         |
'---------------------------------------------------------------------------'
**/
using namespace std;

void out2(string name, int pathwayImpact, int haplo, int test, int covariancematrix, int df_L1, int lfd, struct MAP *map, int nr1, int nr2, double pVal, 
          struct STATplus result, struct ALLELECODE *codesA, int *xvec, int *zvec, int *SigmaAuxVector) {
    
	fstream Combi2;
	//Combi2.open(name.c_str(), ios::out);
	//Combi2 << "try\n";
	FILE *fptr=NULL;
	fptr=fopen(name.c_str(),"a");
		
	//Combi2 << "No\tChr_1\trs_No_1\tPos_No_1\tp-Single-marker_1\tChr_2\trs_No_2\tPos_No_2\tp-Single-marker_2\tp-value\tp-value_corr";
	fprintf(fptr,"%d\t %s\t%s\t%d\t%g \t%s\t%s\t%d\t%g\t%g\txxx",lfd,map[nr1].chr,map[nr1].rs,map[nr1].pos,map[nr1].p,map[nr2].chr,map[nr2].rs,map[nr2].pos,map[nr2].p,pVal);
		 
		 if (test < 3 || test==15 || test==16 || test==17)
			{
			  //do nth
			}
		 else
			{
			 //Combi2 << "\tSNP1_Allele_A\tSNP1_Allele_B\tSNP2_Allele_A\tSNP2_Allele_B\tbeta_x1\tse_x1\tbeta_x1D\tse_x1D\tbeta_x2\tse_x2\tbeta_x2D\tse_x2D\tbeta_x1x2\tse_x1x2\tbeta_x1x2D\tse_x1x2D\tbeta_x1Dx2\tse_x1Dx2\tbeta_x1Dx2D\tse_x1Dx2D";
             fprintf(fptr,"\t%s\t%s\t%s\t%s",codesA[nr1].a1,codesA[nr1].a2,codesA[nr2].a1,codesA[nr2].a2);
			 
			 for (int i =1; i<9; i++)
				{
				 if (xvec[i] != zvec[i])
					{
					  fprintf(fptr,"\t%g\t%g", result.b[i],result.betaNew_se[i]);
					}
				 else
					{
					 fprintf(fptr,"\t-\t-");
					}
				}
			 
			 if(covariancematrix)
			    {
			      for(int m=0; m < ((df_L1)*(df_L1+1))/2; m++)
					 {
					  fprintf(fptr,"\t%g",result.sigma1[SigmaAuxVector[m]]);
					 }
			    }			  
			}		 
		
	if (pathwayImpact == 1)
		{		  
		 fprintf(fptr,"pathway\n");
		}
	else
		{		  
		 fprintf(fptr,"\n");
		}
		
	//Combi2.close();
	fclose(fptr);

}


void out2header(string name, int pathwayImpact, int haplo, int test, int covariancematrix, int df_L1) {
    
	fstream Combi2;
	Combi2.open(name.c_str(), ios::out);
    //Combi2 << "Try\n";
	
		 if (test < 3 || test==15 || test==16)
			{
			  Combi2 << "No\tChr_1\trs_No_1\tPos_No_1\tp-Single-marker_1\tChr_2\trs_No_2\tPos_No_2\tp-Single-marker_2\tp-value\tp-value_corr";
			}
		 else if(!haplo)
			{
			 Combi2 << "No\tChr_1\trs_No_1\tPos_No_1\tp-Single-marker_1\tChr_2\trs_No_2\tPos_No_2\tp-Single-marker_2\tp-value\tp-value_corr\tSNP1_Allele_A\tSNP1_Allele_B\tSNP2_Allele_A\tSNP2_Allele_B\tbeta_x1\tse_x1\tbeta_x1D\tse_x1D\tbeta_x2\tse_x2\tbeta_x2D\tse_x2D\tbeta_x1x2\tse_x1x2\tbeta_x1x2D\tse_x1x2D\tbeta_x1Dx2\tse_x1Dx2\tbeta_x1Dx2D\tse_x1Dx2D";

			 if(covariancematrix)
			    {
			      for(int m=0; m<df_L1; m++)
				     {
				      for(int s=m; s<df_L1; s++)				    
				        {
				         Combi2 << "\tsigma1[" << m <<"]" << "[" << s << "]";
				        }
				     }

			    }			  
			}
		  else if (haplo)
			      {
			       Combi2 << "No\tChr_1\trs_No_1\tPos_No_1\tp-Single-marker_1\tChr_2\trs_No_2\tPos_No_2\tp-Single-marker_2\tp-value\tp-value_corr\tSNP1_Allele_A\tSNP1_Allele_B\tSNP2_Allele_A\tSNP2_Allele_B\tbeta_h1\tse_h1\tbeta_h2\tse_h2\tbeta_h3\tse_h3\tbeta_h4\tse_h4\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_\tbeta_\tse_";
			      }
		
	if (pathwayImpact == 1)
		{		  
		 Combi2 << "\tPathway\n";
		}
	else
		{		  
		 Combi2 << "\n";
		}
		
	Combi2.close();

}
