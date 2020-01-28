#include <stdio.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <set>
#include <map>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits>
#include "intervalfile.h"

using namespace std;

// qsort chromosomes

set<string> recognized_chr;
set<string> unrecognized_chr;

//vector<string> features;


int qs_struct_chr(struct GFF items[], int left, int right);
int qs_struct_chr(struct GFF items[], int left, int right)
  {
    register int a, b;
    int x;
    string tmp;
    struct GFF temp;

    a = left; b = right;
    int a2=0;
    int b2=0;
    if(a2!=a || b2!=b)
    {
      a2=a;
      b2=b;
    }
    x = items[(left+right)/2].chr;

    do
    {
      while(items[a].chr<x && (a < right)) a++;
      while(items[b].chr>x && (b > left)) b--;
      if(a <= b)
      {
        temp = items[a];
        items[a] = items[b];
        items[b] = temp;
        a++; b--;
      }
    } while(a <= b);

    if(left < b) qs_struct_chr(items, left, b);
    if(a < right) qs_struct_chr(items, a, right);
    return 1;
  }

void quick_struct_chr(struct GFF items[], int count);

void quick_struct_chr(struct GFF items[], int count)
  {
   qs_struct_chr(items,0,count-1);
  }

// qsort GFFMERGE

int qs_struct_chr_merge(struct GFFMERGE items[], int left, int right);
int qs_struct_chr_merge(struct GFFMERGE items[], int left, int right)
  {

    register int a, b;
    int x;
    string tmp;
    struct GFFMERGE temp;

    a = left; b = right;
    int a2=0;
    int b2=0;
    if(a2!=a || b2!=b)
    {
      a2=a;
      b2=b;
    }
    x = items[(left+right)/2].chr;

    do
    {
      while(items[a].chr<x && (a < right)) a++;
      while(items[b].chr>x && (b > left)) b--;
      if(a <= b)
      {
        temp = items[a];
        items[a] = items[b];
        items[b] = temp;
        a++; b--;
      }
    } while(a <= b);

    if(left < b) qs_struct_chr_merge(items, left, b);
    if(a < right) qs_struct_chr_merge(items, a, right);
    return 1;
 }

void quick_struct_chr_merge(struct GFFMERGE items[], int count);

void quick_struct_chr_merge(struct GFFMERGE items[], int count)
  {
   qs_struct_chr_merge(items,0,count-1);
  }



//qsort positions
int qs_struct_pos(struct GFF items[], int left, int right);

int qs_struct_pos(struct GFF items[], int left, int right)
  {

    register int a, b;
    int x;
    struct GFF temp;

    a = left; b = right;
    x = items[(left+right)/2].start;

    do
    {
      while( (items[a].start < x) && (a < right)) a++;
      while( (items[b].start > x) && (b > left)) b--;
      if(a <= b)
      {
        temp = items[a];
        items[a] = items[b];
        items[b] = temp;
        a++; b--;
      }
    } while(a <= b);

    if(left < b) qs_struct_pos(items, left, b);
    if(a < right) qs_struct_pos(items, a, right);
    return 1;
  }

void quick_struct_pos(struct GFF items[], int begin, int count);

void quick_struct_pos(struct GFF items[], int begin, int count)
  {
    qs_struct_pos(items,begin,count-1);
  }

//qsort MERGE positions
int qs_struct_pos_merge(struct GFFMERGE items[], int left, int right);

int qs_struct_pos_merge(struct GFFMERGE items[], int left, int right)
  {

    register int a, b;
    int x;
    struct GFFMERGE temp;

    a = left; b = right;
    x = items[(left+right)/2].start;

    do
    {
      while( (items[a].start < x) && (a < right)) a++;
      while( (items[b].start > x) && (b > left)) b--;
      if(a <= b)
      {
        temp = items[a];
        items[a] = items[b];
        items[b] = temp;
        a++; b--;
      }
    } while(a <= b);

    if(left < b) qs_struct_pos_merge(items, left, b);
    if(a < right) qs_struct_pos_merge(items, a, right);
    return 1;
  }

void quick_struct_pos_merge(struct GFFMERGE items[], int begin, int count);

void quick_struct_pos_merge(struct GFFMERGE items[], int begin, int count)
  {
    qs_struct_pos_merge(items,begin,count-1);
  }




// Check chromosome; chromosome name recognition

int guess_chr(char *str1);

int guess_chr(char *str1)
{
  string tempchr=" ";
  int guess=0;

  tempchr.assign(str1);
    if(tempchr.length()>3 && (tempchr.substr(0,3)=="CHR" || tempchr.substr(0,3)=="Chr"))
      {
	tempchr.assign(tempchr.substr(3,tempchr.length()));
      }
    if(((tempchr.length()==1 || tempchr.length()==2) && (atoi(tempchr.c_str()) > 0 && atoi(tempchr.c_str()) < 27 && isdigit(tempchr[0]))))
      {
	guess=atoi(tempchr.c_str());
      }
    // Non-autosomal chromosomes
    else if(tempchr.length()==1 && (tempchr.substr(0,1)=="X" || tempchr.substr(0,1)=="x"))
      {
	guess=23;
      }
    else if(tempchr.length()==1 &&  (tempchr.substr(0,1)=="Y" || tempchr.substr(0,1)=="y"))
      {
	guess=24;
      }
    else if(tempchr.length()==2 && (tempchr.substr(0,2)=="XY" || tempchr.substr(0,2)=="xy"))
      {
	guess=25;
      }
    else if((tempchr.length()==2 && (tempchr.substr(0,2)=="MT" || tempchr.substr(0,2)=="mt" || tempchr.substr(0,2)=="Mt")) || (tempchr.length()==1 && (tempchr.substr(0,1)=="M" || tempchr.substr(0,1)=="m") ))
      {
	guess=26;
      }

    if(guess!=0)
      {
		recognized_chr.insert(tempchr);
      }
     else if(guess==0)
      {
	unrecognized_chr.insert(tempchr);
      }
    return guess;
}


int get_gff_data(string intervalfile, int flanking, int merging, int expandIntervals, int catIntervals, string intervalfile_format, int featurecol, struct GFFCHR * gffchr, fstream &errorfile, fstream &logfile, string outputname,int intervaleditor, int verbose);

int get_gff_data(string intervalfile, int flanking, int merging, int expandIntervals, int catIntervals, string intervalfile_format, int featurecol,  struct GFFCHR * gffchr, fstream &errorfile, fstream &logfile, string outputname,int intervaleditor, int verbose)
{
int i,j,k;
int lines=0;
int chrlines=0;
int tmpchrom=0;

int cols=3;
int chrcol=0;
int startcol=1;
int endcol=2;

 char *str2=NULL;
 char *str3=NULL;



 for(i=0; i<26; i++)
 {
   gffchr[i].relmapstart=0; // relative start position in *.map
   gffchr[i].relmapend=0; // relative end position in *.map
   gffchr[i].relmergestart=-1; // relative start position in gffmerge
   gffchr[i].relmergeend=-1;  // relative end position in gffmerge
   gffchr[i].mergein=0;
   gffchr[i].mapin=0;
 }


// Fill gffchr[i],absstart/end;

char s[500000]; //First line of input file should not have more characters than this!!!!!
char s2[500000];
char c;
char *str1=NULL;
// char *str2=NULL;

//bool commentline;
bool chrfound;
string uniqtemp=" ";
string tempchr=" "  ;
//char commentchar='#';

int l=0;


FILE *gff=NULL;

int chr=0;
int start=0;
int end=0;
char *feature=NULL;
char *attribute=NULL;
char *gene_id=NULL;
char *gene_name=NULL;
char *sep=NULL;

 char *param=NULL;
 char **mparam=NULL;
 int nrparams=0;
 int headerlines=0;
 int tmpheaderlines=0;
 vector<char*> params;

    char *bf=NULL;

if(intervalfile!=" " && intervalfile_format==" ")
   {
     cout<<"Proceed with default parameters 'CHR'=1,'START'=2,END='3';";
     cout<<"\nIf necessary, modify with INTERVAL_FORMAT CHR=\"[1-N]\";START=\"[1-N]\";END=\"[1-N]\";FEATURE=\"[1-N]\";"<<endl;
   }

   //   commentchar='#';
   chrcol=0;
   startcol=1;
   endcol=2;
   featurecol=-9;
   if(intervalfile_format!=" "){
   cout << "\nInterpreting INTERVALFILE_FORMAT string "<<intervalfile_format<<"\n"<<endl;
   logfile<< "\nInterpreting INTERVALFILE_FORMAT string "<<intervalfile_format<<"\n"<<endl;
   bf=(char*)intervalfile_format.c_str();

	str1=strtok(bf,";");

	params.push_back(str1);
	nrparams++;

	while(str1!=NULL)
	  {
     	    str1=strtok(NULL,";");

	    if(str1!=NULL)
	      {

	    params.push_back(str1);
		nrparams++;
	      }
	  }

	for(i=0; i<nrparams; i++)
	  {
	    str2=params[i];
	    // str2=mparam[i];
	    str2=strtok(str2,"=");
 	    param=str2;
	    str2=strtok(NULL,"=");
 	    str2=strtok(str2,"\"\'");
 	    if(strcmp(param,"HEADER")==0||strcmp(param,"HEAD")==0)
	      {
		headerlines=atoi(str2);
		cout <<param<<":\t\t "<<headerlines << " headerlines expected."<<endl;
		logfile <<param<<":\t\t "<<headerlines << " headerlines expected."<<endl;
	      }
	    else if(strcmp(param,"CHR")==0||strcmp(param,"CHROMOSOME")==0||strcmp(param,"CHRCOL")==0)
	      {
		chrcol=atoi(str2)-1;
		cout <<param<<":\t\t chromosome will be read from column nr. "<<chrcol+1 << "."<<endl;
		logfile <<param<<":\t\t chromosome will be read from column nr. "<<chrcol+1 << "."<<endl;
	      }
 else if(strcmp(param,"START")==0||strcmp(param,"STARTCOL")==0)
	      {
		startcol=atoi(str2)-1;
		cout <<param<<":\t\t interval start will be read from column nr. "<<startcol+1 << "."<<endl;
		logfile <<param<<":\t\t interval start will be read from column nr. "<<startcol+1 << "."<<endl;
	      }
 else if(strcmp(param,"END")==0||strcmp(param,"ENDCOL")==0)
	      {
		endcol=atoi(str2)-1;
		cout <<param<<":\t\t interval end will be read from column nr. "<<endcol+1 << "."<<endl;
		logfile <<param<<":\t\t interval end will be read from column nr. "<<endcol+1 << "."<<endl;
	      }
 else if(strcmp(param,"FEAT")==0||strcmp(param,"FEATURE")==0||strcmp(param,"FEATURECOL")==0)
	      {
		featurecol=atoi(str2)-1;
		if(featurecol>=0)
		  {
		    cout <<param<<":\t\t feature will be read from column nr. "<<featurecol+1 << "."<<endl;
		    logfile <<param<<":\t\t feature will be read from column nr. "<<featurecol+1 << "."<<endl;
		  }
	        else
		  {
		    cout <<"FEATURE:\t\t feature was not provided."<<endl;
		    logfile <<"FEATURE:\t\t feature was not provided."<<endl;
		    featurecol=-9;
		  }
	      }
	    str2=strtok(NULL,"\"");
	  }
      }
    c='0';





    //check  lines and columns
     cout << "Opening intervalfile "<<intervalfile<<"...\n";
     logfile << "Opening intervalfile "<<intervalfile<<"...\n";

     int readlines=countLines(intervalfile);

     lines=readlines-headerlines;
     cols=countColumnsWithHeader(intervalfile,headerlines, lines);

     cout << cols << " columns in " << intervalfile<<endl;

     cout <<"\n"<< headerlines<<" header lines + "<<lines << " interval lines in " << intervalfile;





     gff=fopen(intervalfile.c_str(),"r");

    // Parameter-QC
    if(intervalfile_format!=" " || cols<3){
    if(cols<3 || cols > 100)
      {
	cout << "Invalid number of columns ("<<cols<<") in INTERVALFILE"<<endl;
	logfile << "Invalid number of columns ("<<cols<<") in INTERVALFILE"<<endl;
	errorfile << "Invalid number of columns ("<<cols<<") in INTERVALFILE"<<endl;
	exit(1);
      }
    if(chrcol>=cols)
      {
	cout << "Chromosome column larger than number of cols!"<<endl;
	logfile << "Chromosome column larger than number of cols!"<<endl;
	errorfile << "Chromosome column larger than number of cols!"<<endl;
	exit(1);
      }
    if(startcol>=cols)
      {
	cout << "Start column larger than number of cols!"<<endl;
	logfile << "Start column larger than number of cols!"<<endl;
	errorfile << "Start column larger than number of cols!"<<endl;
	exit(1);
      }
    if(endcol>=cols)
      {
	cout << "End column larger than number of cols!"<<endl;
	logfile << "End column larger than number of cols!"<<endl;
	errorfile << "End column larger than number of cols!"<<endl;
	exit(1);
      }
    if(featurecol!=-9 && featurecol>=cols)
      {
	cout << "Feature column larger than number of cols!"<<endl;
	logfile << "Feature column larger than number of cols!"<<endl;
	errorfile << "Feature column larger than number of cols!"<<endl;
	exit(1);
      }
    cout<<"\n";
}
    i=0;


    // read intervalfile

    while(!feof(gff))
    {
        c=fgetc(gff);
	if(headerlines!=0)
         {
	   for(tmpheaderlines=0; tmpheaderlines<headerlines; tmpheaderlines++)
	     {
	       c=fgetc(gff);
 	       if(feof(gff)){
		   cout<<"File ended while reading header."<<endl;
		   logfile<<"File ended while reading header."<<endl;
		   errorfile<<"File ended while reading header."<<endl;
		   exit(1);
		 }
	       while(c!='\n'){
		   c=fgetc(gff);
		   if(feof(gff)){
		       cout<<"File ended while reading header."<<endl;
		       logfile<<"File ended while reading header."<<endl;
		       errorfile<<"File ended while reading header."<<endl;
		       exit(1);
		     }

		 }
	             lines++;
	     }
	   headerlines=0;
	   if(feof(gff))
	     {
	       cout<<"File unexpectedly ended."<<endl;
	       logfile<<"File unexpectedly ended."<<endl;
	       errorfile<<"File unexpectedly ended."<<endl;
	       exit(1);
	     }
         }

	else if(c=='\n')
        {

	  lines++;
	  gffdata = (struct GFF *) realloc(gffdata, (chrlines+1) * sizeof(struct GFF));
 	  if(featurecol<0)
	    {
	      gffdata[chrlines].feature=NULL;
	    }

          chrfound=0;
          s2[i]='\0';
          i=0;

	  str1 = strtok(s2, "\t ");
	  j=0;

	  while(j<cols)
            {

	      if (j==chrcol)
                {

		  chr=guess_chr(str1);
		   if(chr==0)
		    {
		      break;
		    }
		  gffdata[chrlines].chr = 0;
		  gffdata[chrlines].chr=chr-1;
		}
	      else if(j==startcol)
                {
		  start=atoi(str1);
		   if(start==0)
		    {
		      break;
		    }
		  gffdata[chrlines].start=start;

                }
	      else if(j==endcol)
                {
		  end=atoi(str1);
		   if(end==0)
		    {
		      break;
		    }
		  gffdata[chrlines].end=end;
                }

 	      else if (featurecol>=0 && j==featurecol)
 		{
		  gffdata[chrlines].feature = NULL;
		  gffdata[chrlines].feature=(char *) realloc(gffdata[chrlines].feature, (strlen(str1)+ 1) * sizeof(char));
		  strcpy(gffdata[chrlines].feature,str1);
 		}
	      if(j<cols)
	      	{
		  str1 = strtok(NULL, "\t ");
		  if(str1==NULL)
		    {
		      str1=(char*)"-";
		    }
		}

	if(chr!=0)
	   {
	     gffchr[chr-1].in=1;
	   }
       	j++;
	    }

	  if(chr!=0 && start!=0 && end!=0)
	    {
	      // Check Identicals
	      if (chrlines!=0)
		{
	      if(gffdata[chrlines].chr==gffdata[chrlines-1].chr && gffdata[chrlines].start==gffdata[chrlines-1].start && gffdata[chrlines].end==gffdata[chrlines-1].end)
		    {

		      cout << "WARNING! known interval encountered and dismissed in line "<<lines<<" ";
		      cout << " chr:"<<gffdata[chrlines].chr+1<<", start:"<<gffdata[chrlines-1].start <<", end:"<<gffdata[chrlines-1].end <<"."<<endl;
		      logfile << "WARNING! Identical/invalid interval detected and dismissed in line "<<lines<<" ";
		      logfile << " chr:"<<gffdata[chrlines].chr+1<<", start:"<<gffdata[chrlines-1].start <<", end:"<<gffdata[chrlines-1].end <<"."<<endl;
		      errorfile << "WARNING! Identical/invalid interval detected and dismissed in line "<<lines<<" ";
		      errorfile << " chr:"<<gffdata[chrlines].chr+1<<", start:"<<gffdata[chrlines-1].start <<", end:"<<gffdata[chrlines-1].end <<"."<<endl;

		      chrlines--;
		    }
		}
	      chrlines++;

	    }
	}
        else
        {

            s2[i++]=c;
        }

    } //end while

    fclose(gff);
    int chrqcin=0;
    cout << chrlines << " accepted lines.\n" << endl;
    logfile <<"\n"<< lines << " lines found in " << intervalfile << "."<<endl;
    logfile << chrlines << " accepted lines.\n" << endl;

if(recognized_chr.size()!=0)
  {
cout<< "Recognized chromosomes in intervalfile: ";
logfile<< "Recognized chromosomes in intervalfile: ";
 for (set<string>::iterator it=recognized_chr.begin(); it!=recognized_chr.end(); it++)
   {
    cout << " " << *it<< " ";
    logfile << " " << *it<< " ";
   }
  }

cout <<endl;

if(unrecognized_chr.size()!=0)
  {
cout<< "\nCould not interpret element: ";
logfile<< "\nCould not interpret element: ";
logfile<< "\nWARNING: Could not interpret element: ";
 for (set<string>::iterator it=unrecognized_chr.begin(); it!=unrecognized_chr.end(); it++)
   {
    cout << " " << *it<< " ";
    logfile << " " << *it<< " ";
    errorfile << " " << *it<< " ";
   }

cout <<endl;
logfile <<endl;
errorfile <<endl;
  }

 recognized_chr.clear();
 unrecognized_chr.clear();

    int temppos;
    int nrswitches=0;

    // Positions Check
if(chrlines>2)
  {

    cout << "\nChecking if start and end positions of intervals are in the correct order...";
    logfile << "\nChecking if start and end positions of intervals are in the correct order...";
    for(i=0; i<chrlines; i++)
    {
    if ((gffdata[i].end-gffdata[i].start)<0)
    {
      temppos=gffdata[i].start;
      gffdata[i].start=gffdata[i].end;
      gffdata[i].end=temppos;
      nrswitches++;
    }
    }
    if(nrswitches!=0)
    {
    errorfile << "\nChecking if start and end positions of intervals are in the correct order...";
    cout<<"\n"<<nrswitches << " start and end positions have been switched on recognized chromosomes.\n"<<endl;
    errorfile<<"\n"<<nrswitches << " start and end positions have been switched on recognized chromosomes.\n"<<endl;
    }
    else
    {
    	cout << " Ok!\n"<<endl;
    }
  }


// Check Positions, fix PAR


if(chrlines>2)
  {
cout << "\nSorting chromosomes...\t";
quick_struct_chr(gffdata, chrlines);

int tmpchr=-1;

for(i = 0; i<chrlines; i++)
{
if(tmpchr!=gffdata[i].chr)
{
cout << gffdata[i].chr+1 <<" ";
tmpchr=gffdata[i].chr;
}
}

cout << "\n"<<endl;
  }
i=0;
j=0;

// Find relative positions of chromosomes

for(i=1; i<chrlines; i++)
{
  if(i==1)
    {
      gffchr[gffdata[i-1].chr].start=0;
      if(gffdata[i-1].chr!=gffdata[i].chr)
	{
	  gffchr[gffdata[i-1].chr].end=i;
	}
      else if(gffdata[i-1].chr==gffdata[i].chr)
	{
	  gffchr[gffdata[i-1].chr].end=i;
	}
    }
  else if(gffdata[i].chr!=gffdata[i-1].chr)
    {
      gffchr[gffdata[i].chr].start=i;
      gffchr[gffdata[i-1].chr].end=i-1;
    }
  else if(i==chrlines-1)
    {
      gffchr[gffdata[i].chr].end=i;
      if(gffdata[i-1].chr!=gffdata[i].chr)
	{
	  gffchr[gffdata[i].chr].start=i;
	}
    }
  // cout<<gffdata[i-1].start<<" " <<gffdata[i-1].end<<endl;
  // cout<<gffdata[i].start<<" " <<gffdata[i].end<<endl;
}


// Print chromosome start- and endpoints
if(verbose==2 && chrlines>2)
  {
 cout << "Relative positions of chromosomes in interval list:"<<endl;
 cout << "Chr\tfrom\tto"<<endl;
for(i=0; i<26; i++)
{
if(gffchr[i].in==1)
{
cout << i+1<<"\t"<< gffchr[i].start+1 <<"\t" << gffchr[i].end+1 <<endl;
}
}
cout <<endl;

  }

// Sort according to positions
if(chrlines>2)
  {
cout << "Sorting positions on chr ";
for(i=0; i< 26; i++)
{
if(gffchr[i].in==1)
{
cout << i+1<<" ";
quick_struct_pos(gffdata, gffchr[i].start, gffchr[i].end+1);
}
}
cout <<" \n" <<endl;
  }

// Print chromosome spans

// cout << "Interval coverage on individual chromosomes (not counting uncovered start and end segments of chromosomes):\n"<<endl;
// cout << "Chromosome\tSpan(bp)\t\tSumIntervalLengths\t\tCoverage"<<endl;
//
// for(i=0; i<26; i++)
// {
//   int intervalLengths=0;
//   int intervalSpan=0;
// if(gffchr[i].in==1)
// {
//   for(j=gffchr[i].start; j<gffchr[i].end+1; j++)
//     {
//       intervalLengths+=gffdata[j].end-gffdata[j].start+1;
//     }
//
//   intervalSpan=gffdata[gffchr[i].end].end - gffdata[gffchr[i].start].start + 1;
//
//   if(i==24)
//     {
//       intervalSpan=gffdata[gffchr[i].end].end-gffdata[gffchr[i].start].start+1-(par[0].endgap-par[0].startgap+1);
//     }
//   cout << i+1<<"\t\t"<< intervalSpan <<"\t\t"<< intervalLengths <<"\t\t"<<(float)intervalLengths/(float)intervalSpan<<endl;
// }
// }
//
// cout <<"\n";
//
mergelines=chrlines; // Holds if there is no merging


// Check if all starts are sorted
if(chrlines>2)
{
cout << "Checking order of start positions...";
logfile << "Checking order of start positions...";

 for(i=1; i<chrlines; i++)
   {

     gffchr[gffdata[i-1].chr].end=i-1;
     if(gffdata[i-1].chr==gffdata[i].chr)
       {
	 if(gffdata[i-1].start>gffdata[i].start)
	   {
	     cout <<"\nProblem with sorting of positions on chr "<< gffdata[i].chr+1 <<":" <<gffdata[i].start <<"-" << gffdata[i].end<<endl;
	     logfile <<"\nProblem with sorting of positions on chr "<< gffdata[i].chr+1 <<":" <<gffdata[i].start <<"-" << gffdata[i].end<<endl;
	     errorfile <<"\nProblem with sorting of positions on chr "<< gffdata[i].chr+1 <<":" <<gffdata[i].start <<"-" << gffdata[i].end<<endl;
	     exit(1);
	   }
       }
   }
cout<<"\n";
 }

// MERGING

// Check if merging necessary

// Check if merging necessary
if(chrlines>2){
    cout << "Checking if there is an overlap between intervals...";
    logfile << "Checking if there is an overlap between intervals...";
  }

int overlappingIntervals=0;

for(i=1; i<chrlines; i++)
  {
    if(gffdata[i-1].chr==gffdata[i].chr)
      {
	if(gffdata[i-1].end>=gffdata[i].start)
	  {
	    cout <<" Yes!\n"<<endl;
	    logfile <<" Yes!\n"<<endl;
	    overlappingIntervals=1;
	    break;
	  }
      }
  }

if(overlappingIntervals==0 && chrlines>2)
  {
    cout <<" No!\n"<<endl;
    logfile <<" No!\n"<<endl;
  }
if(merging==1 && overlappingIntervals==0)
  {
    cout << "MERGE_INTERVALS is not necessary! There is no overlap between intervals."<<endl;
    logfile << "MERGE_INTERVALS is not necessary! There is no overlap between intervals."<<endl;
  }
 else if(merging==1 && overlappingIntervals==1)
   {
     cout << "Merging intervalfile intervals...\n"<<endl;
     logfile << "Merging intervalfile intervals...\n"<<endl;

     // Intervals are extended to the limits of their overlap; their number is NOT reduced

     int overlapBins=0;
     char *tmpfeature=NULL;

     int oldbin=0;
     int newchr=1;
     int stepback=1;
     int i2=0;
     string tmpstr=" ";
     // Number of lines in future merged struct

     mergelines=0;
     i2=0; // future mergelines;

     for(i=0; i<chrlines; i++)
       {

	 gfftmp = (struct GFFMERGE *) realloc(gfftmp, (i2+1) * sizeof(struct GFFMERGE));

	 gfftmp[i2].chr=gffdata[i].chr;
	 gfftmp[i2].end=gffdata[i].end;
	 gfftmp[i2].start=gffdata[i].start;
	 if(featurecol>=0)
	   {
	     gfftmp[i2].feature = NULL;
	     gfftmp[i2].feature = (char*) realloc(gfftmp[i2].feature, (strlen(gffdata[i].feature)+ 1) * sizeof(char));
	     gfftmp[i2].feature = gffdata[i].feature;
	   }

	 merged=0;
	 mergeend=gffdata[i].end;
	 mergestart=gffdata[i].start;
	 if(featurecol>=0)
	   {
	     tmpfeature=NULL;
	     tmpfeature=gffdata[i].feature;
	   }
	 while(gffdata[i].chr==gffdata[i+1].chr && i<chrlines-1 && gffdata[i+1].start <= mergeend)
	   {

	     merged=1;
	     if(featurecol>=0)
	       {

		 //  tmpfeature=gffdata[i+1].feature;

		 sprintf(tmpfeature,"%s,%s",tmpfeature,gffdata[i+1].feature);



	       }
	     if(gffdata[i+1].end >= mergeend)
	       {
		 mergeend=gffdata[i+1].end;
	       }

	     gfftmp[i2].end=mergeend;
	     i++;
	   }


	 if(merged==1)
	   {
	     if(featurecol>=0)
	       {
		 gfftmp[i2].feature = NULL;
		 gfftmp[i2].feature = (char*) realloc(gfftmp[i2].feature, (strlen(tmpfeature)+ 1) * sizeof(char));
		 gfftmp[i2].feature = tmpfeature;
	       }
	   }
	 merged=0;
	 i2++;
       }
     mergelines=i2;
     cout <<"After MERGE_INTERVALFILE_INTERVALS: "<<mergelines<<" intervals obtained from "<<chrlines<<" intervals."<<endl;
     logfile <<"After MERGE_INTERVALFILE_INTERVALS: "<<mergelines<<" intervals obtained from "<<chrlines<<" intervals."<<endl;

   }
cout <<"\nResulting intervals: "<<mergelines<<endl;
logfile <<"\nResulting intervals: "<<mergelines<<endl;

gffmerge = (struct GFFMERGE *) realloc(gffmerge, (mergelines) * sizeof(struct GFFMERGE));

if(merging==1 && overlappingIntervals ==1)
  {
    for(i=0; i<mergelines; i++)
      {
	gffmerge[i].chr=gfftmp[i].chr;
	gffmerge[i].start=gfftmp[i].start;
	gffmerge[i].end=gfftmp[i].end;
	if(featurecol>=0)
	  {
	    gffmerge[i].feature = NULL;
	    gffmerge[i].feature = (char*) realloc(gffmerge[i].feature, (strlen(gfftmp[i].feature)+ 1) * sizeof(char));
	    //	    strcpy( gffmerge[i].feature, gfftmp[i].feature);
	    gffmerge[i].feature = gfftmp[i].feature;
	  }
      }
    overlappingIntervals=0;
  }
// else if( overlappingIntervals == 0)
 else
   {
     for(i=0; i<mergelines; i++)
       {
	 gffmerge[i].chr=gffdata[i].chr;
	 gffmerge[i].start=gffdata[i].start;
	 gffmerge[i].end=gffdata[i].end;

	 if(featurecol>=0)
	   {
	     gffmerge[i].feature = NULL;
	     gffmerge[i].feature = (char*) realloc(gffmerge[i].feature, (strlen(gffdata[i].feature)+ 1) * sizeof(char));
	     strcpy(gffmerge[i].feature,gffdata[i].feature);
	   }
       }
   }
 free(gffdata);

i=0;

// Merging finished

// Find NEW relative positions of chromosomes

for(i=1; i<mergelines; i++)
{
  if(i==1)
    {
      gffchr[gffmerge[i-1].chr].start=0;
      if(gffmerge[i-1].chr!=gffmerge[i].chr)
	{
	  gffchr[gffmerge[i-1].chr].end=i;
	}
      else if(gffmerge[i-1].chr==gffmerge[i].chr)
	{
	  gffchr[gffmerge[i-1].chr].end=i-1;
	}
   }
  else if(gffmerge[i].chr!=gffmerge[i-1].chr)
    {
      gffchr[gffmerge[i].chr].start=i;
      gffchr[gffmerge[i-1].chr].end=i-1;
    }
  else if(i==mergelines-1)
    {
      gffchr[gffmerge[i].chr].end=i;
      if(gffmerge[i-1].chr!=gffmerge[i].chr)
	{
	  gffchr[gffmerge[i].chr].start=i;
	}
    }
}


// FLANKING

  int y1=0;
  int y2=0;

if(flanking>=1)
  {
    cout << "FLANKING each interval with "<<flanking<<" bp"<<endl;
    logfile << "FLANKING each interval with "<<flanking<<" bp"<<endl;

 for(i=0; i<mergelines; i++)
   {

	 if(gffmerge[i].start-flanking <1 )
	   {
	     gffmerge[i].start=1;
	   }
	 else
	   {
	     gffmerge[i].start=gffmerge[i].start-flanking;
	   }
	 gffmerge[i].end=gffmerge[i].end+flanking;
   }
 cout <<"\n"<<endl;
 logfile <<"\n"<<endl;


// Find NEW gffchr start, end again

for(i=1; i<mergelines; i++)
  {
    if(i==1)
      {
	gffchr[gffmerge[i-1].chr].start=0;
	if(gffmerge[i-1].chr!=gffmerge[i].chr)
	{
	  gffchr[gffmerge[i-1].chr].end=0;
	}
      }
    else if(gffmerge[i].chr!=gffmerge[i-1].chr)
      {
	gffchr[gffmerge[i].chr].start=i;
	gffchr[gffmerge[i-1].chr].end=i-1;
      }
    else if(i==mergelines-1)
      {
	gffchr[gffmerge[i].chr].end=i;
	if(gffmerge[i-1].chr!=gffmerge[i].chr)
	{
	  gffchr[gffmerge[i].chr].start=i;
	}
      }
  }
}


// FILL_GAPS

if(expandIntervals==1)
{

cout << "\nExecuting FILL_GAPS..."<<endl;
cout << "Checking if there is an overlap between intervals... ";
logfile << "\nExecuting FILL_GAPS..."<<endl;
logfile << "Checking if there is an overlap between intervals... ";

overlappingIntervals=0;
for(i=1; i<mergelines; i++)
   {
     if(gffmerge[i-1].chr==gffmerge[i].chr)
       {
	 if(gffmerge[i-1].end>=gffmerge[i].start)
	   {
	     cout <<" Yes!\n"<<endl;
	     logfile <<" Yes!\n"<<endl;
	     if(flanking==0 && merging==0)
	       {
		 cout << "MERGE_INTERVALS is required for FILL_GAPS." <<endl;
		 logfile << "MERGE_INTERVALS is required for FILL_GAPS." <<endl;
		 errorfile << "MERGE_INTERVALS is required for FILL_GAPS." <<endl;
		 overlappingIntervals=1;
		 exit(1);
	       }
	     if(flanking!=0 && merging==0)
	       {
		 cout << "Non-overlapping intervals are required for FILL_GAPS. Set FLANKING to 0." <<endl;
		 logfile << "Non-overlapping intervals are required for FILL_GAPS. Set FLANKING to 0." <<endl;
		 errorfile << "Non-overlapping intervals are required for FILL_GAPS. Set FLANKING to 0." <<endl;
		 overlappingIntervals=1;
		 exit(1);
	       }
	     if(flanking==0 && merging!=0)
	       {
		 cout << "MERGE_INTERVALS is required for FILL_GAPS." <<endl;
		 logfile << "MERGE_INTERVALS is required for FILL_GAPS." <<endl;
		 errorfile << "MERGE_INTERVALS is required for FILL_GAPS." <<endl;
		 overlappingIntervals=1;
		 exit(1);
	       }
	     if(flanking!=0 && merging!=0)
	       {
		 cout << "Non-overlapping intervals are required for FILL_GAPS. Set MERGE_INTERVALS to 1 and FLANKING to 0." <<endl;
		 logfile << "Non-overlapping intervals are required for FILL_GAPS. Set MERGE_INTERVALS to 1 and FLANKING to 0." <<endl;
		 errorfile << "Non-overlapping intervals are required for FILL_GAPS. Set MERGE_INTERVALS to 1 and FLANKING to 0." <<endl;
		 overlappingIntervals=1;
		 exit(1);
	       }


	   }
       }
  }

 if(overlappingIntervals==0)
   {
	     cout <<" No. Proceed to FILL_GAPS."<<endl;
	     logfile <<" No. Proceed to FILL_GAPS."<<endl;
   }


for(i=0; i<mergelines-1; i++)
{
  // Chromososome start/ends
      if(i==0)
	{
	  gffmerge[i].start=gffchr[gffmerge[i].chr].absstart;
	}
      else if(gffmerge[i].chr!=gffmerge[i+1].chr)
	{
	  if(gffmerge[i].chr!=25-1 && gffmerge[i+1].chr!=25-1)
	    {
	      gffmerge[i].end=gffchr[gffmerge[i].chr].absend;
	      gffmerge[i+1].start=gffchr[gffmerge[i+1].chr].absstart;
	    }
	}
	  // If PAR
      else if(gffmerge[i].chr==25-1)
	{
	  if(gffmerge[i].end<=par[0].xabsend && gffmerge[i+1].start>=par[1].xabsstart)
	    {
	      gffmerge[i].end=par[0].xabsend;
	      gffmerge[i+1].start=par[1].xabsstart;
	    }
	}
      else if(i==mergelines-2)
	{
	  gffmerge[i+1].end=gffchr[gffmerge[i+1].chr].absend;
	}
      // Maximal extension for chromosome ends
      int ext=30000000;
      if(i==0)
	{
	  gffmerge[i].start=1;
	}
      else if(gffmerge[i].chr!=gffmerge[i+1].chr)
	{
	  if(gffmerge[i].chr!=25-1 && gffmerge[i+1].chr!=25-1)
	    {
	      gffmerge[i].end=gffmerge[i].end+ext;
	      gffmerge[i+1].start=1;
	    }
	}
      else if(i==mergelines-2)
	{
	  gffmerge[i+1].end=gffmerge[i+1].end+ext;
	}
 }

int randn=0;
int diff=0;
for(i=0; i<mergelines-1; i++)
  {
  if(gffmerge[i].chr==gffmerge[i+1].chr && gffmerge[i].end!=gffmerge[i+1].start-1)
    {
      diff=gffmerge[i+1].start-gffmerge[i].end-1;
	  gffmerge[i+1].start=gffmerge[i+1].start-diff/2;
	  gffmerge[i].end=gffmerge[i].end+diff/2;
	  if(diff%2==1)
	    {
	      gffmerge[i].end=gffmerge[i].end+1;
	    }
    }
  }
cout << "\n"<<endl;
logfile<< "\n"<<endl;
} // FILL_GAPS

// Concatenate intervals
if(verbose==2 && chrlines>2)
  {
 cout << "Relative positions of chromosomes in interval list:"<<endl;
 cout << "Chr\tfrom\tto"<<endl;
for(i=0; i<26; i++)
{
if(gffchr[i].in==1)
{
cout << i+1<<"\t"<< gffchr[i].start+1 <<"\t" << gffchr[i].end+1 <<endl;
}
}
cout <<endl;
  }
if(catIntervals!=0)
  {
    cout<<"Concatenating intervals..."<<endl;
    logfile<<"Concatenating intervals..."<<endl;
    int i2=0;
    int i3=0;
    char *tmpfeature=NULL;
    char *tmpfeature2=NULL;
    int catIntervals2=0;
    for(i=0; i<26; i++)
      {
	gffchr[i].catint=0;
	gffchr[i].catintlast=0;
	gffchr[i].interv=0;
      }

    if(expandIntervals==0)
      {
	cout<<"It is sometimes recommended to use FILL_GAPS with CONCATENATE_INTERVALS!"<<endl;
      }

    tmpchrom=0;
for(i=0; i<mergelines; i++)
  {
    gffchr[gffmerge[i].chr].mergein=1;
    gffchr[gffmerge[i].chr].relmergestart=-1;
    gffchr[gffmerge[i].chr].relmergeend=-1;
  }
for(i=0; i<mergelines; i++)
  {
    if(tmpchrom!=gffmerge[i].chr)
      {
	tmpchrom=gffmerge[i].chr;
      }
    if(gffchr[gffmerge[i].chr].mergein==1)
      {
	if(gffchr[gffmerge[i].chr].relmergestart==-1)
	  {
	    gffchr[gffmerge[i].chr].relmergestart=i;
	  }
	gffchr[gffmerge[i].chr].relmergeend=i;
      }
  }


    int catlines=0;

    if(verbose==2){
    cout << "\nChr\tOldIntervals\tNewIntervals\n"<<endl;
    }
     for(i=0; i<26; i++)
      {
	if(gffchr[i].mergein==1)
	  {
	    gffchr[i].interv=gffchr[i].end-gffchr[i].start+1;
	    gffchr[i].catint=gffchr[i].interv/catIntervals;

	    if(gffchr[i].catint==0 && gffchr[i].interv!=0)
	      {
		gffchr[i].catint=1;
		gffchr[i].catintlast=gffchr[i].interv;
	      }
	    else if(gffchr[i].catint!=0)
	      {
		if(gffchr[i].interv%catIntervals==0)
		  {
		    gffchr[i].catintlast=catIntervals;
		  }
		else if(gffchr[i].interv%catIntervals!=0)
		  {
		    gffchr[i].catintlast=gffchr[i].interv%catIntervals;
		    gffchr[i].catint++;
		  }
	      }
	    if(verbose==2){
	    cout <<i+1<<"\t"<<gffchr[i].interv<<"\t"<<gffchr[i].catint<<endl;
	    }
	    catlines=catlines+gffchr[i].catint;

	  }
      }
    cout <<"Concatenating resulted in "<<catlines<<" new intervals.\n"<<endl;
    logfile <<"Concatenating resulted in "<<catlines<<" new intervals.\n"<<endl;





    int catpar1=0;
    int catpar2=0;
    int intervpar1=0;
    int intervpar2=0;
    int catparlast1=0;
    int catparlast2=0;


    gfftmp = (struct GFFMERGE *) realloc(gfftmp, catlines * sizeof(struct GFFMERGE));
    //   i2 tracks chromosome catlines
    //   i3 tracks mergelines;
    //   i4 tracks catines;
    //   i5 tracks lines in mergeinterval;
    int i4=0;
    int i5=0;
    i3=0;
    mergelines=catlines;
    int catrepair=0;
    int parcatint=gffchr[25-1].catint;
    gffchr[i].catpar1=catpar1;
    gffchr[i].catpar2=catpar2;
    gffchr[i].catparlast1=catparlast1;
    gffchr[i].catparlast2=catparlast2;
    for(i=0; i<26; i++)
      {
	if(i!=25-1 || ((catpar1 ==0) && (catpar2==0)))
	  {
	    for(i2=0; i2<gffchr[i].catint; i2++)
	      {
		if(featurecol>0)
		  {
		    gfftmp[i4].feature=NULL;
		    gfftmp[i4].feature = (char*) realloc(gfftmp[i4].feature, 50000 * sizeof(char));
		    strcpy(gfftmp[i4].feature,gffmerge[i3].feature);
		  }

	    if(i2!=gffchr[i].catint-1)
	      {
		gfftmp[i4].chr=i;

		    gfftmp[i4].start=gffmerge[i3].start;
		    gfftmp[i4].end=gffmerge[i3+catIntervals-1].end;

		if(featurecol>0)
		  {
		    for(i5=i3+1; i5<i3+catIntervals;i5++)
		      {
			strcat(gfftmp[i4].feature,",");
			strcat(gfftmp[i4].feature,gffmerge[i5].feature);
		      }
		  }
		i3=i3+catIntervals;
		i4++;
	      }
	    else if(i2==gffchr[i].catint-1)
	      {
		gfftmp[i4].chr=i;
		gfftmp[i4].start=gffmerge[i3].start;
		gfftmp[i4].end=gffmerge[i3+gffchr[i].catintlast-1].end;
		if(featurecol>0)
		  {
		    strcpy(gfftmp[i4].feature,gffmerge[i3].feature);
		    for(i5=i3+1; i5<i3+gffchr[i].catintlast;i5++)
		      {
			strcat(gfftmp[i4].feature,",");
			strcat(gfftmp[i4].feature,gffmerge[i5].feature);
		      }
 		  }
		i4++;
		i3=i3+gffchr[i].catintlast;

	      }
	  }
   }
      }

    //    gffchr[25-1].catint=parcatint;

    for(i=1; i<mergelines; i++)
      {

	if(i==1)
	  {
	    gffchr[gfftmp[i-1].chr].start=0;
	    if(gfftmp[i-1].chr!=gfftmp[i].chr)
	      {
		gffchr[gfftmp[i-1].chr].end=0;
	      }
	  }
	else if(gfftmp[i].chr!=gfftmp[i-1].chr)
	  {
	    gffchr[gfftmp[i].chr].start=i;
	    gffchr[gfftmp[i-1].chr].end=i-1;
	  }
	else if(i==mergelines-1)
	  {
	    gffchr[gfftmp[i].chr].end=i;
	    if(gfftmp[i-1].chr!=gfftmp[i].chr)
	      {
		gffchr[gfftmp[i].chr].start=i;
	      }

	  }
      }

    gffmerge=(struct GFFMERGE *) realloc(gffmerge,(catlines) * sizeof(struct GFFMERGE));
    i=0;

    for(i2=0; i2<catlines; i2++)
      {

	gffmerge[i2].chr=gfftmp[i2].chr;
	gffmerge[i2].start=gfftmp[i2].start;
	gffmerge[i2].end=gfftmp[i2].end;
	if(featurecol>0)
	  {
	    gffmerge[i2].feature=NULL;
	    gffmerge[i2].feature = (char *) realloc( gffmerge[i2].feature,strlen(gfftmp[i2].feature+1) * sizeof(char));
	    gffmerge[i2].feature=gfftmp[i2].feature;
	  }
      }
    cout <<"... Done!"<<endl;
    mergelines=catlines;
  }



for(i=0; i<26; i++)
  {
    if(gffchr[i].mergein==1)
      {
	gffchr[i].relmergestart=-1;
	gffchr[i].relmergeend=-1;
      }
  }


tmpchrom=0;
for(i=0; i<mergelines; i++)
  {
    gffchr[gffmerge[i].chr].mergein=1;
  }
for(i=0; i<mergelines; i++)
  {
    if(tmpchrom!=gffmerge[i].chr)
      {
	tmpchrom=gffmerge[i].chr;
      }
    if(gffchr[gffmerge[i].chr].mergein==1)
      {
	if(gffchr[gffmerge[i].chr].relmergestart==-1)
	  {
	    gffchr[gffmerge[i].chr].relmergestart=i;
	  }
	gffchr[gffmerge[i].chr].relmergeend=i;
      }
  }

if(verbose==2 && chrlines>2)
  {
cout<<"\nIntervals that could be assigned to chromosomes:\n" <<endl;
cout << "Chr\tfrom\tto\n";
for(i=0; i<26; i++)
  {
    if(gffchr[i].mergein==1)
      {
	cout <<i+1 <<"\t"<<gffchr[i].relmergestart+1<<"\t"<<gffchr[i].relmergeend+1<<endl;
      }
  }

  }


int deep=0;
const char *intervalfull=NULL;


if(intervaleditor){

string intervalbase=intervalfile.substr(intervalfile.find_last_of('/')+1);
string intervalfileout=outputname+"_"+intervalbase+"_modified.txt";
const char *modified="_modified.txt";

cout <<"\nModified intervalfile will be written to "<< intervalfileout<<endl;
logfile <<"\nModified intervalfile will be written to "<< intervalfileout<<endl;

fstream bf_out;

bf_out.open(intervalfileout.c_str(), ios::out);

bf_out <<"#Intervalfile_modified_from_"<<intervalfile<<":MERGE_INTERVALS="<<merging <<";FLANKING="<<flanking<<";FILL_GAPS="<<expandIntervals<<";CONCATENATE_INTERVALS="<<catIntervals<<endl;

for(i=0; i<mergelines; i++)
  {
    bf_out<<gffmerge[i].chr+1<<"\t"<<gffmerge[i].start<<"\t"<<gffmerge[i].end;

    if(featurecol>=0)
      {
	bf_out<<"\t"<<gffmerge[i].feature;
      }
    bf_out<<"\n";
  }

bf_out.close();
 }
pfeaturecol=&featurecol;

return mergelines;

}
