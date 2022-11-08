/*
 Goal: calculate covariance matrix between all possible gene combinations of two chromosome files
Consider only samples with genes with a frequency below a bound of Bound (=0.9)

# Input: annovar freqs
# Output: User definition --> Example: All.chr1.exonic_variant_function_paired_frequencies
*/


#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <math.h> 
#include <sys/types.h>
#include <unistd.h>
#include <time.h>

#define max(a,b) (a<b?b:a)
#define min(a,b) (a<b?a:b)

struct GENE {
  char name[100];
  int numsamples;
  int *freq;               // list of frequencies 
  float average,std;
  int rares;
};

struct CHROM
{
 char fname[128];
 int Ngenes;
 struct GENE *genes[200000];  // Max 100000 genes per chrom
};

struct GENE *Malloc_gene( int Nsamples){

   struct GENE *g;
   g=malloc(sizeof(struct GENE));
   g->freq=malloc(sizeof (int)*(Nsamples+1));
   return(g);
}
int Maxsamples = 1000;
int minsamples = 10; 
int below10 = 0;
int numcells =0;
FILE *outpairs;
/*
struct CHROM *Malloc_chrom(){
    return(malloc(sizeof(struct CHROM)));
}
*/

int *muttmp;

float Bound1,Bound2;



read_chrom_row( char *name, int *Nsamples, struct CHROM *chr){
/*
# Example without any ","
GRHPR 0 0.07396449320577009 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.40437687695736646 0 0 0  0 0 0 0 0 0.10699008866159604TUCS ...

*/

   int Ngenes=0, rares,first;
   FILE *fp;
   int freq;
   float sumfreq,sumfreqsq,average;
   char  gene[100],c;
   int numsamples;

   fp = fopen(name,"r");
   first=1;

   while (!feof(fp)) {
    fscanf(fp,"%s",gene);
    if (feof(fp)) break;
    fscanf(fp,"%c",&c);   // " "

    chr->genes[Ngenes]=Malloc_gene(*Nsamples);
    strcpy(chr->genes[Ngenes]->name,gene);
    numsamples=0;
    rares=0;
    sumfreq=0;sumfreqsq=0;
    while (c==' ' && !feof(fp)) {
      int numinp = fscanf(fp,"%d",&freq);
      if (numinp == 0)
          break;
      chr->genes[Ngenes]->freq[numsamples]=freq;
      sumfreq+=freq;
      sumfreqsq+=freq*freq;
      if (freq < 0.5) rares++;
      numsamples++;

      fscanf(fp,"%c",&c);   // " "  or \n
    }
    chr->genes[Ngenes]->freq[numsamples]=-1000;   // Mark the last element
    chr->genes[Ngenes]->numsamples=numsamples;
    if (first) {
       if (numsamples>*Nsamples) {
            printf("Too many samples. Use -N ##\n");
            exit(1);
       }
       *Nsamples = numsamples;
       first=0;
    }
    if (numsamples != *Nsamples)
         printf("Error: the number of samples in %s %d is not the expected %d.\n", gene, numsamples, *Nsamples);
    if (numsamples==0)
         printf("Error: the number of samples in %s %d is not the expected %d.\n", gene, numsamples, *Nsamples);
    average=chr->genes[Ngenes]->average=sumfreq/numsamples;
    chr->genes[Ngenes]->std=sqrt((double) ( ( sumfreqsq - sumfreq*sumfreq/numsamples)/numsamples));
    chr->genes[Ngenes]->rares=rares;
//    if (chr->genes[Ngenes]->std==(float)0.0)
//         printf("Error: %s:  nul standard deviation.\n", gene);

//    if (chr->genes[Ngenes]->rares > Nsamples/5 )
//         printf("GeneRare: %s %d\n", gene, rares);
    Ngenes++;


   }
   chr->Ngenes=Ngenes;
}


void sub(struct GENE *g1, struct GENE *g2){

int i,j,max;
float statistic=0; //Gtest;

float  Oold[3][3], O[3][3], marg[2][3]; // first index, snp1, more mutated in tumor;  second index, snp2
float  E[3][3]; // expected numbers

   for (i=0;i<3;i++)
    for (j=0;j<3;j++){
     O[i][j]=0.0;
   }
// input code 0=none, 1=somatic, 2=germline


int *fp1,*fp2; // normals, first;  tumor, afterwards
max = g1->numsamples;

for (fp1=g1->freq, fp2=g2->freq;*fp1 != -1000; fp1++, fp2++) {
     O[*fp1][*fp2]++;
 }

/* The O[2][2] should not be too big because subjects with bb have germline mutation in both genes in all tissues so it is not a driver pair. */
/* More evidence in favor than against */

int a=O[1][1];
int b=O[1][2];
int c=O[2][1];

//if (O[2][2] > (a+b+c)) return;


 for (i=0;i<3;i++){
  marg[0][i]=0;
  for (j=0; j<3; j++)
    marg[0][i]+=O[i][j];
 }

 for (j=0;j<3;j++){
  marg[1][j]=0;
  for (i=0; i<3; i++)
    marg[1][j]+=O[i][j];
 }

 for (i=0;i<3;i++)
    for (j=0;j<3;j++){
     E[i][j]=(float)marg[0][i]*marg[1][j]/max;
    }

/* If any of the three cells has more observations than expected then we consider it. Otherwise, not considered.*/

if (E[1][1] >= O[1][1] && E[1][2] >= O[1][2] && E[2][1] >= O[2][1]) return;

/* If the mutation in normal tissue is too big it could not be driver  */

//if (E[2][2] < O[2][2] && O[2][2] > 2*(a+b+c) ) return;

//if (O[2][2] > 0.1*max) return;

//if (marg[0][2] > 0.2*max && marg[1][2] > 0.2*max) return;

/* Calculate the x=O[2][2] for which the G test statistic is minimum.
   Keeping the 3 cells, it gives the degree of dependency created only by those cells and freeing the others as much as possible. */

//Make a copy of Observations 
for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      Oold[i][j]=O[i][j];

float x=-1;

//  Constants


int f=marg[1][0];
int e=marg[1][1];
int d=marg[1][2];
int g=marg[0][0];
int h=marg[0][1];
int ii=marg[0][2];

float c1=f-h+a+b-ii +c;      // n_w - m_s + o_{ss} + o_{sb} - m_b + o_{bs} = n_w -o_ws + o_sb -m_b  = o_ww + o_sw -o_ws - o_sb - o_bs - o_bb

float c3=d-b;
float c4=ii-c;

if (c1+c3+c4 != 0) 
  x=(c3*c4)/(c1+c3+c4);  // denom = f - h + a  + d  = n_w - m_s +o_ss + n_b
else 
  x=0;

if ( (x >= 0) && ((c1 + x) >= 0) && (x <= c3) && (x <= c4))  {  // always TRUE, see the paper
/*G test for independence*/

O[0][0] = c1 + x;
O[0][2] = c3 - x;
O[2][0] = c4 - x;
O[2][2] = x;

statistic=0;
 for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      if (0 < O[i][j])
        statistic += O[i][j]*logf(O[i][j]/E[i][j]);

statistic *=2; // Only for Gtest

if (statistic<= 20.0) 
    below10++;
else 
{
printf("%8s %8s ",g1->name, g2->name);
   for (i=0;i<3;i++)
     for (j=0;j<3;j++) 
       printf("%3.0f ",Oold[i][j]);
printf(" %3.3f\n",statistic );
}


}  // always valid 


}

int samegene(char *g1, char *g2){
char *p1, *p2;

  p1 = g1 + (strlen(g1) -1);
  p2 = g2 + (strlen(g2) -1);

  while (*p1 == *p2 && *p1 != ':'){
    p1--; p2--;

  }

return(*p1 == *p2);
}


main(int argc, char**argv){

  struct CHROM chr1, chr2;
  int i,j,k;
  int Nsamples;
  int nosame,rowformat,numfile,L;
  char listname[100],nomfiles[2][100];
  rowformat=0;
  Bound2=1.1;
  Bound1=0.0;
  Nsamples=Maxsamples;
  numfile=0;
  k = 1;
  while (k<argc)
  {
    if (argv[k][0]=='-')
    { 
     L = strlen(argv[k]); 

     if ((argv[k][1]=='B') && (argv[k][2]=='1') &&(L==3)){ ++k; Bound1=atof(argv[k]); 
                                         }
     else if ((argv[k][1]=='B') && (argv[k][2]=='2') &&(L==3)){ ++k; Bound2=atof(argv[k]); 
                                         }

     else if ((argv[k][1]=='R')&&(L==2)) {rowformat=1;}
     else if ((argv[k][1]=='N')&&(L==2)) { ++k; Nsamples=atoi(argv[k]);}
     else if ((argv[k][1]=='l')&&(L==2)) { ++k; sprintf(listname,"%s",argv[k]); 
                                        } 
     else if ((argv[k][1]=='k')&&(L==2)) { ++k; minsamples=atoi(argv[k]);
                                         }
     else {
       printf("Error in -options.\n");
       exit(1);
     }                          
    }
    else {
       sprintf(nomfiles[numfile],"%s",argv[k]);
       numfile++;
    }
    ++k;
   }
   if (numfile !=2 || Nsamples==0 ) {
       printf("Error in parameters.\n");
       exit(1);
   }
 
   read_chrom_row(nomfiles[0],&Nsamples, &chr1);
   read_chrom_row(nomfiles[1],&Nsamples, &chr2);

   nosame=strcmp(nomfiles[1],nomfiles[0]);
 

   for (i=0; i<chr1.Ngenes; i++) {
     for (j=(nosame?0:i+1); j<chr2.Ngenes; j++) {
	 if (samegene(chr1.genes[i]->name,chr2.genes[j]->name)==0)
          sub(chr1.genes[i],chr2.genes[j]);
//	 else
//           printf("%s %s\n",chr1.genes[i]->name,chr2.genes[j]->name); // test
        }
   }
char sizefilename[200];
FILE *sizefile;
    strcpy(sizefilename,nomfiles[0]);
    int namesize= strlen(nomfiles[1]);
    strcat(strcat(sizefilename,nomfiles[1]+(namesize-4)),"Pairs.size");
    sizefile=fopen(sizefilename,"w");
    fprintf(sizefile,"%d",below10);
    fclose(sizefile);
}
