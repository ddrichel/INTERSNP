using namespace std;

// Gene

struct RARES
{
  int in;
  char *rs;
  int chr;
  int pos;
};
struct RARES *raresnps = NULL;

struct GENE
{
  char *gene_name;
  char *gene_id;
  int start;
  int end;
  int relstart;
  int relend;
};

struct GENE *gene = NULL;


// GFF

struct GFF
{
  char *gene_id;
  char *gene_name;

  char *source;
  char *feature;
  char *attr;                                                    
  
  int chr; // Chromosome                                                       
  int start;
  int end;
    

  char qcin;                                                   
  //  char *rs; // rs-Nummer                                                      
   
};

struct GFF *gffdata = NULL; //reference legend file


struct GFFMERGE
{
  char *feature;
  char *attr;                                                    
  
  int start;
  int end;
                                                       
  int relmapstart;
  int relmapend;
  int relstart; // number of line where interval begins
  int relend;
  int merged;
    
  int chr; // Chromosome
  int mapin;                                                   
  //  char *rs; // rs-Nummer                                                      

  int hit;
};
struct GFFMERGE *gffmerge = NULL; //reference legend file
struct GFFMERGE *gfftmp = NULL; //reference legend file
struct GFFMERGE *gfftmp2 = NULL; //reference legend file
 
struct GFFCHR
{
    
  bool in;
  bool mapin;
  bool mergein;
  bool hit;
    int start;
    int end;
    int relmapstart;
    int relmapend;
    int relmergestart;
    int relmergeend;
 
  int absstart;
  int absend;
  int relstart; // number of line in gff where chromosome begins
  int relend;

  int interv;
  int catint;
  int catintlast;
  int catparlast1;
  int catparlast2;
  int catpar1;
  int catpar2;


};

struct GFFCHR gffchr[26]; //reference legend file

struct PAR
{
//     bool in;
//   bool mapin;
//   bool mergein;
//   bool hit;
//     int start;
//     int end;
//     int relmapstart;
//     int relmapend;
//     int relmergestart;
//     int relmergeend;
//  
  int startgap;
  int endgap;



  int yabsstart;
  int yabsend;
  int xabsstart;
  int xabsend;


  //  int relstart; // number of line in gff where chromosome begins
  //  int relend;
};

struct PAR par[2]; 


int mergelines=0;

string chrname=" ";    

struct INTERVALS
{
  int relmapstart; 
  int relmapend;
  int start;
  int end;
  int chr;
  char *feature;

  // Number of Rare variants in Interval; max BINSIZE
  int n;
  // For seed intervals: number of dense bined intervals
  int N;
  // BINSIZE
  int s1;
  // BINSHIFT
  int s2;
  // SEED interval number for BIN_ADJUS,
  int seed;

};

struct INTERVALS *intervals, *intervalstmp, *intervals2, *seeds;

int mergestart;
int mergeend;
int merged=1;

int merging=0;
int fixChromosomes=1;
int flanking=0;

int expandIntervals=0;
int catIntervals=0;
int minRareInBin=1;
int maxRareInBin=1000000;
string intervalfile=" ";
string intervalfile_out=" ";

string gene_ids=" ";
string gene_names=" ";

int nrgenes=0;

int attr=0;

int *nvariations;
int *seedintervals;

int featurecol=-9;

int *pfeaturecol;

