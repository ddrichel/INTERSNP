#include <stdio.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream> // wegen Dateistreamobjekt
#include <set> // wegen Container setn
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits>




using namespace std;


int rsNumberSearchFam(char *pid, struct MAP3 *map3,int left, int right) // binäre Suche
{
	int l=left;
	int r=right;
	int x; //

	while(r >= l)
	{
		x=(l+r)/2;

		if(strcmp(pid, map3[x].pid) < 0)/* kleiner? */
		{
	 		r=x-1;  /* Rechte Seite ist nicht mehr so interessant */
		}
		else

		{
	 		l=x+1;  /* Linke Seite ist nicht mehr so interessant */
		}

		if(strcmp(pid, map3[x].pid) == 0)
		{
			return map3[x].line;     /* Gefunden; x = Position*/
		}
	}
	return -1;
}

// map3 nach family-IDs sortieren
void qsortmap3(struct MAP3 **map3,int left, int right, struct MAP3 *x, struct MAP3 *y)
{
	//sorts in ascending order
	int i,j;


	(*x).fid=NULL;(*x).fid=(char *)realloc((*x).fid,(strlen((*map3)[(left+right)/2].fid)+1)*sizeof(char));
	(*x).pid=NULL;(*x).pid=(char *)realloc((*x).pid,(strlen((*map3)[(left+right)/2].pid)+1)*sizeof(char));
	(*x).vid=NULL;(*x).vid=(char *)realloc((*x).vid,(strlen((*map3)[(left+right)/2].vid)+1)*sizeof(char));
	(*x).mid=NULL;(*x).mid=(char *)realloc((*x).mid,(strlen((*map3)[(left+right)/2].mid)+1)*sizeof(char));
	(*y).fid=NULL;
	(*y).pid=NULL;
	(*y).vid=NULL;
	(*y).mid=NULL;

	i=left;j=right;

	strcpy((*x).fid,(*map3)[(left+right)/2].fid);
	strcpy((*x).pid,(*map3)[(left+right)/2].pid);
	strcpy((*x).vid,(*map3)[(left+right)/2].vid);
	strcpy((*x).mid,(*map3)[(left+right)/2].mid);
	(*x).line=(*map3)[(left+right)/2].line;

	do
	{
		while( strcmp((*map3)[i].pid,(*x).pid)< 0 && (i <right)) {i++;}
		while( strcmp((*x).pid,(*map3)[j].pid)< 0 && (j>left)) {j--;}

		if (i<=j)
		{
			(*y).fid=(char*)realloc((*y).fid,(strlen((*map3)[i].fid)+1)*sizeof(char));
			(*y).pid=(char*)realloc((*y).pid,(strlen((*map3)[i].pid)+1)*sizeof(char));
			(*y).vid=(char*)realloc((*y).vid,(strlen((*map3)[i].vid)+1)*sizeof(char));
			(*y).mid=(char*)realloc((*y).mid,(strlen((*map3)[i].mid)+1)*sizeof(char));

			strcpy((*y).fid,(*map3)[i].fid);
			strcpy((*y).pid,(*map3)[i].pid);
			strcpy((*y).vid,(*map3)[i].vid);
			strcpy((*y).mid,(*map3)[i].mid);
			(*y).line=(*map3)[i].line;

			(*map3)[i].fid=(char*)realloc((*map3)[i].fid,(strlen((*map3)[j].fid)+1)*sizeof(char));
			(*map3)[i].pid=(char*)realloc((*map3)[i].pid,(strlen((*map3)[j].pid)+1)*sizeof(char));
			(*map3)[i].vid=(char*)realloc((*map3)[i].vid,(strlen((*map3)[j].vid)+1)*sizeof(char));
			(*map3)[i].mid=(char*)realloc((*map3)[i].mid,(strlen((*map3)[j].mid)+1)*sizeof(char));

			strcpy((*map3)[i].fid,(*map3)[j].fid);
			strcpy((*map3)[i].pid,(*map3)[j].pid);
			strcpy((*map3)[i].vid,(*map3)[j].vid);
			strcpy((*map3)[i].mid,(*map3)[j].mid);
			(*map3)[i].line=(*map3)[j].line;

			(*map3)[j].fid=(char*)realloc((*map3)[j].fid,(strlen((*y).fid)+1)*sizeof(char));
			(*map3)[j].pid=(char*)realloc((*map3)[j].pid,(strlen((*y).pid)+1)*sizeof(char));
			(*map3)[j].vid=(char*)realloc((*map3)[j].vid,(strlen((*y).vid)+1)*sizeof(char));
			(*map3)[j].mid=(char*)realloc((*map3)[j].mid,(strlen((*y).mid)+1)*sizeof(char));

			strcpy((*map3)[j].fid,(*y).fid);
			strcpy((*map3)[j].pid,(*y).pid);
			strcpy((*map3)[j].vid,(*y).vid);
			strcpy((*map3)[j].mid,(*y).mid);
			(*map3)[j].line=(*y).line;

			i++;j--;
		}
	}while (i<=j);

	if((*x).fid!=NULL){free((*x).fid);}
	if((*x).pid!=NULL){free((*x).pid);}
	if((*x).vid!=NULL){free((*x).vid);}
	if((*x).mid!=NULL){free((*x).mid);}

	if((*y).fid!=NULL){free((*y).fid);}
	if((*y).pid!=NULL){free((*y).pid);}
	if((*y).vid!=NULL){free((*y).vid);}
	if((*y).mid!=NULL){free((*y).mid);}
	(*x).fid=NULL;(*y).fid=NULL;
	(*x).pid=NULL;(*y).pid=NULL;
	(*x).vid=NULL;(*y).vid=NULL;
	(*x).mid=NULL;(*y).mid=NULL;

	if (left <j)
	{
		qsortmap3(map3,left,j, x, y);
	}
	if (i < right)
	{
		qsortmap3(map3,i,right, x, y);
	}

};
int ReadFamData(struct FAMILY*& family, struct PERSON *person, struct MAP *map, int nlinestfam, fstream &errorfile, fstream &logfile)
{
	int nlist = 0;
	int found = 0;
	int pline = 0;
	struct MAP3 *map3 = NULL;
	struct MAP3 x1map3;
	struct MAP3 y1map3;
	unsigned short int singleperson = 0;
	unsigned short int trios = 0;
	unsigned short int nuclearFam2 = 0;
	unsigned short int nuclearFam3 = 0;
	unsigned short int nuclearFam4 = 0;
	unsigned short int nuclearFamMoreThan4 = 0;
	unsigned short int childFather = 0;
	unsigned short int childMother = 0;
	unsigned short int nChildrenTotal = 0;
	unsigned short int nfamilies = 0;


	/** create map3 **/
	map3 = (struct MAP3 *) calloc(nlinestfam, sizeof(struct MAP3));
	if (!map3) die("memory allocation error in map3");
	for (int i=0; i<nlinestfam; i++) {
		map3[i].fid = (char *) calloc((strlen(person[i].fid)+1), sizeof(char));
		if (!map3[i].fid) die("memory allocation error in map3[i].fid");
		map3[i].pid = (char *) calloc((strlen(person[i].pid)+1), sizeof(char));
		if (!map3[i].fid) die("memory allocation error in map3[i].pid");
		map3[i].vid = (char *) calloc((strlen(person[i].vid)+1), sizeof(char));
		if (!map3[i].vid) die("memory allocation error in map3[i].vid");
		map3[i].mid = (char *) calloc((strlen(person[i].mid)+1), sizeof(char));
		if (!map3[i].mid) die("memory allocation error in map3[i].mid");

		strcpy(map3[i].fid, person[i].fid);
		strcpy(map3[i].pid, person[i].pid);
		strcpy(map3[i].vid, person[i].vid);
		strcpy(map3[i].mid, person[i].mid);
		map3[i].line = i;
	}

	qsortmap3(&map3,0,nlinestfam-1, &x1map3, &y1map3);

	for(int i=0;i<nlinestfam;i++)
	{
		if (nlist > 0 && (strcmp(family[nlist-1].fid,map3[i].fid)==0))
		{
			pline = rsNumberSearchFam(map3[i].pid, map3, 0, nlinestfam-1); // family-ID wird durch die Zeilennummern ersetzt

			if (pline != -1)
			{
				family[nlist-1].nFamMember++;

				family[nlist-1].personList= (int *) realloc(family[nlist-1].personList,family[nlist-1].nFamMember*sizeof(int));
				if (!family[nlist-1].personList)
				{
					errorfile << "memory allocation error in family[nlist-1].personList\n";
					logfile << "memory allocation error in family[nlist-1].personList\n";
					cout << "memory allocation error in family[nlist-1].personList\n";
					errorfile.close();logfile.close();exit(1);
				}
				family[nlist-1].personList[family[nlist-1].nFamMember-1] = pline;

				if (strcmp(map3[i].vid, "0") != 0 && strcmp(map3[i].mid, "0") != 0)
				{
					family[nlist-1].nChildren++;

					family[nlist-1].childIndex= (int *) realloc(family[nlist-1].childIndex,family[nlist-1].nChildren*sizeof(int));
					if (!family[nlist-1].childIndex)
					{
						errorfile << "memory allocation error in family[nlist-1].childIndex\n";
						logfile << "memory allocation error in family[nlist-1].childIndex\n";
						cout << "memory allocation error in family[nlist-1].childIndex\n";
						errorfile.close();logfile.close();exit(1);
					}

					family[nlist-1].childIndex[family[nlist-1].nChildren-1] = pline;

					pline = rsNumberSearchFam(map3[i].vid, map3, 0, nlinestfam-1);

					if (pline != -1)
					{
						if (person[pline].sex == 1)
						{
							family[nlist-1].fatherIndex = pline;
						}
						else
						{
							family[nlist-1].fatherIndex = pline;
//							cout << "Father " << map3[i].vid << " of child " << map3[i].pid << ": Sex option must be wrong\n";
						}
					}
					else
					{
//						cout << "Father of child " << map3[i].pid << " isn't found in the file\n";
					}
					pline = rsNumberSearchFam(map3[i].mid, map3, 0, nlinestfam-1);

					if (pline != -1 )
					{
						if (person[pline].sex == 2)
						{
							family[nlist-1].motherIndex = pline;
						}
						else
						{
							family[nlist-1].motherIndex = pline;
//							cout << "Mother " << map3[i].mid << " of child " << map3[i].pid << ": Sex option must be wrong\n";
						}
					}
					else
					{
//						cout << "Mother of child " << map3[i].pid << " isn't found in the file\n";
					}
				}
			}
		}
		else
		{
			//neue Familie: Test, ob diese bereits vorhanden ist
			found = 0;
			for (int j = 0; j<nlist; j++)
			{
				if (strcmp(family[j].fid,map3[i].fid)==0)
				{
					pline = rsNumberSearchFam(map3[i].pid, map3, 0, nlinestfam-1); // Person-ID wird durch die Zeilennummern ersetzt

					if (pline != -1)
					{
						family[j].nFamMember++;

						family[j].personList= (int *) realloc(family[j].personList,family[j].nFamMember*sizeof(int));
						if (!family[j].personList)
						{
							errorfile << "memory allocation error in family[j].personList\n";
							logfile << "memory allocation error in family[j].personList\n";
							cout << "memory allocation error in family[j].personList\n";
							errorfile.close();logfile.close();exit(1);
						}
						family[j].personList[family[j].nFamMember-1] = pline;

						if (strcmp(map3[i].vid, "0") != 0 && strcmp(map3[i].mid, "0") != 0)
						{
							family[j].nChildren++;

							family[j].childIndex= (int *) realloc(family[j].childIndex,family[j].nChildren*sizeof(int));
							if (!family[j].childIndex)
							{
								errorfile << "memory allocation error in family[nlist-1].childIndex\n";
								logfile << "memory allocation error in family[nlist-1].childIndex\n";
								cout << "memory allocation error in family[nlist-1].childIndex\n";
								errorfile.close();logfile.close();exit(1);
							}
							family[j].childIndex[family[j].nChildren-1] = pline;

							pline = rsNumberSearchFam(map3[i].vid, map3, 0, nlinestfam-1);
							if (pline != -1)
							{
								if (person[pline].sex == 1)
								{
									family[j].fatherIndex = pline;
								}
								else
								{
									family[j].fatherIndex = pline;
//									cout << "Father " << map3[i].vid << " of child " << map3[i].pid << ": Sex option must be wrong\n";
								}
							}
							else
							{
								cout << "Father of child " << map3[i].pid << " isn't found in the file\n";
							}
							pline = rsNumberSearchFam(map3[i].mid, map3, 0, nlinestfam-1);
							if (pline != -1 )
							{
								if (person[pline].sex == 2)
								{
									family[j].motherIndex = pline;
								}
								else
								{
									family[j].motherIndex = pline;
//									cout << "Mother " << map3[i].mid << " of child " << map3[i].pid << ": Sex option must be wrong\n";
								}
							}
							else
							{
//								cout << "Mother of child " << map3[i].pid << " isn't found in the file\n";
							}
						}
					}
					found = 1;
					break;
				}
			}

			if (!found) // neue Familie, die vorher noch nicht aufgetaucht ist
			{
				nlist++;
				family = (struct FAMILY *) realloc(family, nlist* sizeof(struct FAMILY));
				if (!family)
				{
					errorfile << "memory allocation error in family\n";
					logfile << "memory allocation error in family\n";
					cout << "memory allocation error in family\n";
					errorfile.close();logfile.close();exit(1);
				}

				family[nlist-1].fid=NULL;
				family[nlist-1].personList=NULL;
				family[nlist-1].nFamMember = 0;
				family[nlist-1].nChildren = 0;
				family[nlist-1].fatherIndex = -1;
				family[nlist-1].motherIndex = -1;
				family[nlist-1].childIndex = NULL;

				family[nlist-1].fid = (char *) calloc((strlen(map3[i].fid)+1), sizeof(char));
				if (!family[nlist-1].fid) die("memory allocation error in map3[i].fid");

				strcpy(family[nlist-1].fid, map3[i].fid);

				pline = rsNumberSearchFam(map3[i].pid, map3, 0, nlinestfam-1); // person-ID wird durch die Zeilennummern ersetzt

				if (pline != -1)
				{
					family[nlist-1].nFamMember++;

					family[nlist-1].personList= (int *) realloc(family[nlist-1].personList,family[nlist-1].nFamMember*sizeof(int));
					if (!family[nlist-1].personList)
					{
						errorfile << "memory allocation error in family[nlist-1].personList\n";
						logfile << "memory allocation error in family[nlist-1].personList\n";
						cout << "memory allocation error in family[nlist-1].personList\n";
						errorfile.close();logfile.close();exit(1);
					}
					family[nlist-1].personList[family[nlist-1].nFamMember-1] =pline;

					if (strcmp(map3[i].vid, "0") != 0 && strcmp(map3[i].mid, "0") != 0)
					{
						family[nlist-1].nChildren++;

						family[nlist-1].childIndex= (int *) realloc(family[nlist-1].childIndex,family[nlist-1].nChildren*sizeof(int));
						if (!family[nlist-1].childIndex)
						{
							errorfile << "memory allocation error in family[nlist-1].childIndex\n";
							logfile << "memory allocation error in family[nlist-1].childIndex\n";
							cout << "memory allocation error in family[nlist-1].childIndex\n";
							errorfile.close();logfile.close();exit(1);
						}
						family[nlist-1].childIndex[family[nlist-1].nChildren-1] = pline;
						pline = rsNumberSearchFam(map3[i].vid, map3, 0, nlinestfam-1);
						if (pline != -1)
						{
							if (person[pline].sex == 1)
							{
								family[nlist-1].fatherIndex = pline;
							}
							else
							{
								family[nlist-1].fatherIndex = pline;
//								cout << "Father " << map3[i].vid << " of child " << map3[i].pid << ": Sex option must be wrong\n";
							}
						}
						else
						{
//							cout << "Father of child " << map3[i].pid << " isn't found in the file\n";
						}
						pline = rsNumberSearchFam(map3[i].mid, map3, 0, nlinestfam-1);
						if (pline != -1 )
						{
							if (person[pline].sex == 2)
							{
								family[nlist-1].motherIndex = pline;
							}
							else
							{
								family[nlist-1].motherIndex = pline;
//								cout << "Mother " << map3[i].mid << " of child " << map3[i].pid << ": Sex option must be wrong\n";
							}
						}
						else
						{
//							cout << "Mother of child " << map3[i].pid << " isn't found in the file\n";
						}
					}
				}
			}
		}
	}


	for (int i=0;i< nlist;i++)
	{
		if (family[i].nFamMember >= 3)
		{
			if (family[i].nChildren >= 1)
			{
				switch(family[i].nChildren)
				{
					case 1:
						trios++;
						break;
					case 2:
						nuclearFam2++;
						break;
					case 3:
						nuclearFam3++;
						break;
					case 4:
						nuclearFam4++;
						break;
					default:
						nuclearFamMoreThan4++;
						break;
				}

				nChildrenTotal += family[i].nChildren;
			}
			else
			{
//				cout << "Family without children\n";
				singleperson += 2;
			}
		}
		else if (family[i].nFamMember == 2)
		{

			if (family[i].nChildren > 1)
			{
//				cout << "Children without partens\n";
			}
			else if (family[i].nChildren == 1 && family[i].fatherIndex !=-1 && family[i].motherIndex == -1)
			{
				childFather++;
				nChildrenTotal += family[i].nChildren;
			}
			else if (family[i].nChildren == 1 && family[i].fatherIndex ==-1 && family[i].motherIndex != -1)
			{
				childMother++;
				nChildrenTotal += family[i].nChildren;
			}
			else if (family[i].nChildren == 1 && family[i].fatherIndex !=-1 && family[i].motherIndex != -1)
			{
				for (int j=0; j<family[i].nFamMember;j++)
				{
					if (family[i].fatherIndex == family[i].personList[j])
					{
						childFather++;
						nChildrenTotal += family[i].nChildren;
//						cout << "Child without a mother\n";
						break;
					}
					else if (family[i].motherIndex == family[i].personList[j])
					{
						childMother++;
						nChildrenTotal += family[i].nChildren;
//						cout << "Child without a father\n";
						break;
					}
				}
			}
			else
			{
//				cout << "Family without children\n";
			}
		}
		else if (family[i].nFamMember == 1 && family[i].fatherIndex == -1 && family[i].motherIndex == -1)
		{
			singleperson++;
		}
		nfamilies = nlist;
	}


cout << "\nFamily data:\n";
cout << nlinestfam-nChildrenTotal << " founders\n";
cout << nChildrenTotal << " non-founders\n";
cout << nfamilies << " nuclear families\n";
cout << singleperson << " founder singletons\n";
cout << trios << " trios\n";


	if (nuclearFam2 != 0)
	{
//		cout << nuclearFam2 << " nuclear families with 2 children\n";
	}
	if (nuclearFam3 != 0)
	{
//		cout << nuclearFam3 << " nuclear families with 3 children\n";
	}
	if (nuclearFam4 != 0)
	{
//		cout << nuclearFam4 << " nuclear families with 4 children\n";
	}
	if (nuclearFamMoreThan4 != 0)
	{
//		cout << nuclearFamMoreThan4 << " nuclear families with more than 4 children\n";
	}
	if (childFather != 0)
	{
//		cout << childFather << " child + father\n";
	}
	if (childMother != 0)
	{
//		cout << childMother << " child + mother\n";
	}

	for (int i=0; i<nlinestfam; i++) {
		free(map3[i].fid); map3[i].fid = NULL;
		free(map3[i].pid); map3[i].pid = NULL;
		free(map3[i].vid); map3[i].vid = NULL;
		free(map3[i].mid); map3[i].mid = NULL;
	}
	free(map3); map3 = NULL;

	return nfamilies;
}

double tdt (struct PERSON *person, struct FAMILY *family, int npplqc, int* PPLMap, uint64_t*** BinSNPs, struct PPLLOCATION* PplLocations, int thread, int nfamilies, int snpPos)
{
	struct PPLLOCATION guy;
	int a=0, b=0, c=0, d=0, e=0, f=0;
	int ft=0, fnt=0, mt=0, mnt=0;
	unsigned short int T[2];
	unsigned short int NT[2];
	double stat = 0.00;
	double pValue = 1;
	int df = 1;
	int lMod = 0;

	T[0] = 0;
	T[1] = 0;
	NT[0] = 0;
	NT[1] = 0;
	int aff = 0;


	for (int i=0;i< nfamilies;i++) // Schleife über alle Familien
	{
		for (int j=2; j< family[i].nFamMember;j++)
		{
			for(int l=0;l<npplqc;l++) // persons
			{
				guy = PplLocations[l];
				lMod = PPLMap[l];

				if (person[lMod].qcin == 1)
				{
					if (family[i].childIndex[j-2] == lMod) // Genotypes Child
					{

						if (getbit64(BinSNPs[snpPos][guy.nr][1],guy.pos) == 1)
						{
							e = 1;
							f = 1;
						}
						else if (getbit64(BinSNPs[snpPos][guy.nr][2],guy.pos) == 1)
						{
							e = 1;
							f = 2;
						}
						else if (getbit64(BinSNPs[snpPos][guy.nr][3],guy.pos) == 1)
						{
							e = 2;
							f = 2;
						}
						else
						{
							e = 0;
							f = 0;
						}
					}
					else if(family[i].fatherIndex == lMod) // Genotypes Father
					{
						if (getbit64(BinSNPs[snpPos][guy.nr][1],guy.pos) == 1)
						{
							a = 1;
							b = 1;
						}
						else if (getbit64(BinSNPs[snpPos][guy.nr][2],guy.pos) == 1)
						{
							a = 1;
							b = 2;
						}
						else if (getbit64(BinSNPs[snpPos][guy.nr][3],guy.pos) == 1)
						{
							a = 2;
							b = 2;
						}
						else
						{
							a = 0;
							b = 0;
						}
					}
					else if(family[i].motherIndex == lMod) // Genotypes Mother
					{
						if (getbit64(BinSNPs[snpPos][guy.nr][1],guy.pos) == 1)
						{
							c = 1;
							d = 1;
						}
						else if (getbit64(BinSNPs[snpPos][guy.nr][2],guy.pos) == 1)
						{
							c = 1;
							d = 2;
						}
						else if (getbit64(BinSNPs[snpPos][guy.nr][3],guy.pos) == 1)
						{
							c = 2;
							d = 2;
						}
						else
						{
							c = 0;
							d = 0;
						}
					}
				}
			}

			if (e*f != 0 && (((!(a*e*(a!=e))) && (!(c*f*(c!=f)))) || ((!(a*f*(a!=f))) && (!(c*e*(c!=e))))
			|| ((!(a*e*(a!=e))) && (!(d*f*(d!=f)))) || ((!(a*f*(a!=f))) && (!(d*e*(d!=e))))
			|| ((!(b*e*(b!=e))) && (!(c*f*(c!=f)))) || ((!(b*f*(b!=f))) && (!(c*e*(c!=e))))
			|| ((!(b*e*(b!=e))) && (!(d*f*(d!=f)))) || ((!(b*f*(b!=f))) && (!(d*e*(d!=e)))))) // Überprüfung, ob Mendel erfüllt
			{
				if (a*b*c*d*e*f!=0  && person[lMod].aff[thread] == 2) //if child affected and not missing
				{
					//tdt
					if (a==e && c==f)
					{
						ft = a;
						fnt = b;
						mt = c;
						mnt = d;
					}
					else if (a==e && d==f)
					{
						ft = a;
						fnt = b;
						mt = d;
						mnt = c;
					}
					else if (b==e && c==f)
					{
						ft = b;
						fnt = a;
						mt = c;
						mnt = d;
					}
					else if (b==e && d==f)
					{
						ft = b;
						fnt = a;
						mt = d;
						mnt = c;
					}
					else if (a==f && c==e)
					{
						ft = a;
						fnt = b;
						mt = c;
						mnt = d;
					}
					else if (a==f && d==e)
					{
						ft = a;
						fnt = b;
						mt = d;
						mnt = c;
					}
					else if (b==f && c==e)
					{
						ft = b;
						fnt = a;
						mt = c;
						mnt = d;
					}
					else if (b==f && d==e)
					{
						ft = b;
						fnt = a;
						mt = d;
						mnt = c;
					}

					if (ft != fnt)
					{
						if (ft == 1)
						{
							T[0] += 1;
						}
						if (fnt == 1)
						{
							NT[0] += 1;
						}
						if (ft == 2)
						{
							T[1] += 1;
						}
						if (fnt == 2)
						{
							NT[1] += 1;
						}
					}

					if (mt != mnt)
					{
						if (mt == 1)
						{
							T[0] += 1;
						}
						if (mnt == 1)
						{
							NT[0] += 1;
						}
						if (mt == 2)
						{
							T[1] += 1;
						}
						if (mnt == 2)
						{
							NT[1] += 1;
						}
					}
				}
			}
			else
			{
				if (a*b*c*d == 0)
				{
					// Mendelfehler
				}
			}
			ft=0;
			fnt=0;
			mt=0;
			mnt=0;
			a = 0;
			b = 0;
			c = 0;
			d = 0;
		}
	}

	// Teststatistik
	for (int k = 0; k<2; k++)
	{
		if((T[k]+NT[k]) > 0)
		{
			stat += pow((double)T[k]-NT[k],2)/(T[k]+NT[k]);
		}
	}
	cout << snpPos+1 << "\t" << T[0] << "\t" << T[1] << "\t" << "\n";
	cout << snpPos+1 << "\t" << NT[0] << "\t" << NT[1] << "\t" << "\n";
	return stat;
}
