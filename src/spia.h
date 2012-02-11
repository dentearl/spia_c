#ifndef SPIA_H_
#define SPIA_H_

#include <uthash.h>    /* hashes! */
#include <stdlib.h>

#define MAX_PATH_LENGTH 512 /* file system path names  */
#define MAX_ID_LENGTH 100   /* maximum gene name (id) length */
#define MAX_PATHWAY 1000    /* largest pathway in terms of unique genes*/
#define NUM_REL 28 /* the number of relationships possible. length of 
                      enum relationType.
                   */
#define TRUE 1
/* 
   techincal note here, the actual binding/association relationships use `/'
   in their names, but c complains when we do this, so we have substituted
   `_' in their place.
*/ 
enum relationType {activation, compound, binding_association, expression,
                   inhibition, activation_phosphorylation, phosphorylation,
                   indirect, inhibition_phosphorylation, dephosphorylation_inhibition,
                   dissociation, dephosphorylation, activation_dephosphorylation,
                   state, activation_indirect, inhibition_ubiquination, ubiquination,
                   expression_indirect, indirect_inhibition, repression,
                   binding_association_phosphorylation, dissociation_phosphorylation,
                   indirect_phosphorylation, family_membership,
                   transcriptional_activation, transcriptional_inhibition,
                   process_activation, process_inhibition};
typedef enum relationType relationType;


/* STEVE BENZ SAYS:
  Yea I can give you some coefficients:
  member> family_membership       0
  -t> transcriptional_activation  1
  -t| transcriptional_inhibition -1
  -ap> process_activation         1
  -ap| process_inhibition        -1
*/


typedef struct downList downList;
struct downList {
  char id[MAX_ID_LENGTH];
  relationType relation;
  downList *next;
};

struct upstreamGene {
  char id[MAX_ID_LENGTH];                /* key */
  downList *down;
  UT_hash_handle hh;     /* this struct is hashable */
};
typedef struct upstreamGene upstreamGene;

struct geneItem {
  char id[MAX_ID_LENGTH];    /* key */
  int order; 
  UT_hash_handle hh;
};
typedef struct geneItem geneItem;

typedef struct diffE diffE;
struct diffE{
  char id[MAX_ID_LENGTH];
  double expr;
  int count;
  UT_hash_handle hh;
};

struct allGene {
  char id[MAX_ID_LENGTH];    /* key */
  UT_hash_handle hh;
};
typedef struct allGene allGene;

typedef struct isectList isectList;
struct isectList{
  char id[MAX_ID_LENGTH];
  int origOrder;
  int newOrder;
  isectList *next;
};

typedef struct pGlobal pGlobal;
struct pGlobal{
  double p;
  double bonf;
  double fdr;
  char path[MAX_ID_LENGTH];
  int pathSize;
  int NDE;
  double tAc;
  double pPERT;
  double pNDE;
  UT_hash_handle hh;
};

/* General functions */
void gatherOptions(int argc, char **argv, char **dir, char **de, char **ar, char **spf, char **betaCoFile);
void usage(void);
int endsIn_tab(char *filename);
int readOldPathway(char *filename); // charlie parse format
int readNewPathway(char *filename); // cancer browser format
void readBetaCoeffFile(char *filename);
double processPathway(int *status);
void readDETab(char *filename);
void readArrayTab(char *filename);
int isRelationship(char *rel, relationType *relType_ptr);
void addItemToIntersect(char *id, int orgOrd, int newOrd);
void printIsectList(void);
void deleteIntersect(void);
void printBetaCoeffs(void);
void cleanup(void);
void* de_malloc(size_t n);
void verbose(char const *s, ...);
void debug(char const *s, ...);
void message(char const *t, char const *s, ...);

/* pathway hash and list functions */
int  addGenePath(char *geneID);
int  addInteraction(char *upID, char *downID, relationType rel);
upstreamGene* findGenePath(char *geneID);
void deleteAllPath(void);
void printPathway(void);
int  id_sortPath(upstreamGene *a, upstreamGene *b);
void sort_by_idPath(void);
int  countDowns(char *geneID);
void addGenePathAll(char *geneID);
allGene* findGenePathAll(char *geneID);
void deleteAllPathAll(void);
void printPathwayAll(void);
void printPathwayAllArray(void); // alternative method, using the global array of ptrs
int countIntersect_array_path(void);

/* gene order hash and list functions */
int  addGeneOrder(char *geneID);
void deleteAllOrder(void);
geneItem* findGeneOrder(char *geneID);
int  id_sortOrder(geneItem *a, geneItem *b);
void sort_by_idOrder(void);
void printOrder(void);

/* Differential expression hash and list functions*/
int addDiffExprsGene(char *id, double de); //
int addDiffExprsGeneEntry(char *id, double de); // 
diffE* findDiffExpr(char *geneID);
void printDiffExpr(void);
int countIntersect_de_path(void);

/* Bootstrap related functions*/
void deleteAllBootGenes(void); // kill the boot hash
void populateBootGenes(int n); // populate the boothash once
void printBootGenes(void);
int addBootGene(char *id); //
diffE* findBootGene(char *geneID);
double randomDEvalue(void);
char** randomPathwayGene(void);
double calcBootPerterb(void);
void fillBootDEVec(double *vecNew, int n);
double pPERTcalc(double a, double *b, int n, int opt);

/* All (tested) Genes hash and list functions*/
int addAllGeneEntry(char *id);
allGene* findAllGene(char *geneID);
void printAllGene(void);

/* Linear Algebra functions */
double** buildBeta(relationType rel);
double** zeros(int n);
double* zerosVec(int n);
void printMatrix(double **a, int n);
/*void prettyPrintMatrix(double **a, int n);*/ /* unused, I'm going to remove this soon  */
void printVector(double *a, int n);
void printNamedVector(double *a, int n);       /* uses the ordered hash to print gene names */
void transposeMatrix(double **a, int n);       /* for use with lapack / fortran */
void colNorm(double **a, int n);               /* normalize a matrix by columns*/
double* colSum(double **a, int n);             /* sum of matrix columns. returns vector*/
void matAdd(double **a, double **b, int n);    /* matrix - matrix addition*/
void matScalMult(double **a, int n, double b); /* matrix - scalar multiplication*/
void fillDEVec(double *vecNew, int n);
int  isMatrixEmpty(double **a, int n);
double determinant(double **a, int n);         /* Gaussian elimination determinant of the matrix*/
void subtractIdent(double **a, int n);         /* a = I - a*/
void invert(double **a, int n);                /* a = a^-1 MAKES CALL TO LAPACK*/
void matVecMultiply(double **a, int n, double *b, double *c);  /* result stored in c */
void matMatMultiply(double **a, int n, double **b, double **c); /* result stored in c */
double sumVec(double *a, int n);
int  compare_dbls(const void *a, const void *b);     /* qsort(), ahoy!*/
int  compare_dbls_rev(const void *a, const void *b); /* qsort(), ahoy!*/
double correctTA(double *a, int n);                /* takes sorted *a, finds median, subtracts median*/
void copyMatrix(double **a, double **b, int n); /* copy matrix a into matrix b. */
void solveForPF(double *a, double *b, double *c, int n);

/* Probability functions ... Part of this project was cancelled in mid development.*/
double  probCDFHyper(int x, int m, int n, int k); 
double  probPDFHyper(int x, int m, int n, int k); // unused
double  logOfFactorial(int n);         // unused
double  sigmaLog(int i, int n);
double  combPValue(double a, double b);
double* bonferroni(double *a, int b); // unused
double* falseDR(double *a, int b); // unimplemented

/* Probability values hashes */
int  addPGlobal(double p, char *c, int pathSize, int Nde, double t, double pPERT, double pNDE);
void sort_by_pValue(void);
void rev_sort_by_pValue(void);
void printPValues(void);
void bonferrPGlobal(void);
void fdrPGlobal(void);

#endif // SPIA_H_
