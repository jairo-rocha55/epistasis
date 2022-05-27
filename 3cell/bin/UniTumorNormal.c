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

#define Maxgenes  400000
int Maxsamples = 2000;
int minsamples = 10;
int *muttmp;

struct GENE {
  char name[100];
  int numsamples;
  float *freq;               // list of frequencies 
  float average,std;
  int rares;
};

struct CHROM
{
 char fname[128];
 int Ngenes;
 struct GENE *genes[Maxgenes];  // Max 100000 genes per chrom
};

struct GENE *Malloc_gene( int Nsamples){

   struct GENE *g;
   g=malloc(sizeof(struct GENE));
   g->freq=malloc(sizeof (float)*(Nsamples+1));
   return(g);
}

/*
struct CHROM *Malloc_chrom(){
    return(malloc(sizeof(struct CHROM)));
}
*/

read_chrom_row( char *name, int *Nsamples, struct CHROM *chr){
/*
# Example without any ","
GRHPR 0 0.07396449320577009 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.40437687695736646 0 0 0  0 0 0 0 0 0.10699008866159604
TUCS ...

*/
   
   int Ngenes=0, rares,first;
   FILE *fp;
   float freq,sumfreq,sumfreqsq,average;
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
      fscanf(fp,"%f",&freq);
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

void sub(struct GENE *g1){
int Nsamples=g1->numsamples;
int i,max;
int conttumor=0, contnormal=0;
int vfp1, vn1, mut,out=0;

float *fp1, *n1  ; // normals, first;  tumor, afterwards
int *tmp, *tmpn;

max = g1->numsamples/2;
   for (fp1=g1->freq, n1=fp1+max, tmp= muttmp, tmpn=tmp+max;*n1 != -1000; fp1++, n1++, tmp++, tmpn++)
{
    if (*fp1 != 0 && *n1 != 0)  mut = 2;
    else if  (*fp1 == 0 && *n1 == 0) mut = 0;
    else if  (*fp1<*n1) mut = 1;
    else {mut = 3; out = 1;}
    
    *tmp=mut ;
}
 /*
{
     if (*fp1<*n1) {vfp1=0; vn1=1; mut = 1;}
     else if (*fp1 > *n1) 
          if (*fp1 == 2 ) { vfp1=1; vn1=1; mut=2;}
          else { vfp1=1; vn1=0; mut=3;out=1;}
     else if ( *fp1 >= 1 )  {  vfp1=1; vn1=1; mut=2;}
     else { vfp1=0; vn1=0; mut=0;}
     *tmp=mut;
  //  *tmpn=vn1;
     if (vfp1==1) contnormal++;
     if (vn1==1) conttumor++;
 }
*/
// 0 = n ,  1 = som,  2 = ger, 3 = mut nor, no tumor


   if (out==0 ) { // && conttumor >= minsamples   && contnormal <= 0.8*max ) { 
    *tmp=-1000;
    *tmpn=-1000;
    printf("%s",g1->name);

    for (tmp=muttmp; *tmp != -1000; tmp++) {
     printf(" %d",*tmp);
    }
    printf("\n");


   }   

}


main(int argc, char**argv){

  struct CHROM chr1, chr2;
  int i,j,k;
  int Nsamples;
  int nosame,rowformat,numfile,L;
  Nsamples=atoi(argv[3]); 
  char nomfiles[2][100];
  rowformat=0;
  Nsamples=Maxsamples;
  numfile=0;
  k = 1;
  while (k<argc)
  {
    if (argv[k][0]=='-')
    { 
     L = strlen(argv[k]); 


     if ((argv[k][1]=='N')&&(L==2)) { ++k; Nsamples=atoi(argv[k]);
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
   if (numfile !=1  ) {
       printf("Error in parameters.\n");
       exit(1);
   }
 
   read_chrom_row(nomfiles[0],&Nsamples, &chr1);
   muttmp = malloc(sizeof (float)*(Nsamples+1));

   
   for (i=0; i<chr1.Ngenes; i++)
          sub(chr1.genes[i]);  




}
