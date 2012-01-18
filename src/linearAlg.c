/*
 * Copyright (C) 2009-2012 by
 * Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
 *
 * ... and other members of Josh Stuart's lab (BME Dept. UCSC)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <assert.h>
/*#include <gsl_linalg.h> // gnu sci lib, linear alg. for det()... maybe*/
#include "spia.h"

#define f77_dgetrf dgetrf_
extern void f77_dgetrf (
                        /*
                          dgetrf will fill in the ipiv, the pivot table. 
                          We need this in order to call dgetri.
                        */
                        int *m, int *n, double *a, int *lda, 
                        int *ipiv, int *info
                        );
#define f77_dgetri dgetri_
extern void f77_dgetri (
                        /*
                          dgetri will calculate the inverse of our matrix
                          a, using an LU decomposition approach. Yay.
                        */
                        int *n, double *a, int *lda, int *ipiv,
                        double *workspace, int *lworkspace, int *info
                        );

double** buildBeta(relationType rel){
    /* 
       Given an enum relationType, this function will build the Beta_rel
       matrix by stepping through the pathway struct and hand back a
       pointer to the Beta_rel matrix. You must free the Beta_rel matrix.
       betaCoefs is a vector of beta coefficients.
    */
    extern upstreamGene *pathway_up;
    extern geneItem     *geneOrder;
    //  extern double *betaCoefs;
    upstreamGene *cur;
    geneItem *ord;
    cur = pathway_up;
    ord = geneOrder;
    unsigned int n;
    n = HASH_COUNT(geneOrder);
    double **beta;
    beta = zeros(n);
    int i;
    for (i = 0; i < n; ++i){
        /* look down every column*/
        upstreamGene *g = NULL;
        geneItem     *o = NULL;
        g = findGenePath(ord->id);
        if(g == NULL){
            ord = ord->hh.next;
            continue;
            /* some genes are only downstream and have no
               upstream component */
        }
        downList *d = g->down;
        while(d != NULL){
            if((int) rel == (int) d->relation){
                o = findGeneOrder(d->id);
                assert(o != NULL);
                beta[o->order][i] = 1;
            }
            d = d->next;
        }
        ord = ord->hh.next;
    }
    return beta;
}

void printMatrix(double **a, int n){
    /* for debug purposes only */
    int i, j;
    for(i = 0; i < n; ++i)
        for(j = 0; j< n; ++j)
            fprintf(stderr, "%7.4f%s", a[i][j], (j == (n-1))?"\n":" ");
}

/*
  void prettyPrintMatrix(double **a, int n){
  // for debug purposes only 
  int i, j;
  for(i = 0; i < n; ++i)
  for(j = 0; j< n; ++j)
  printf("%7.4f%s", a[i][j], (j == (n-1))?"\n":" ");
  }
*/

void printVector(double *a, int n){
    int i;
    for(i = 0; i < n; ++i){
        printf(" %d: ",i);
        printf("%7.4f%s", a[i], (i == (n-1))?"\n":" ");
    }
}
void printNamedVector(double *a, int n){
    extern geneItem     *geneOrder;
    geneItem *ord;
    ord = geneOrder;  
    int i;
    for(i = 0; i < n; ++i){
        printf(" %s: ",ord->id);
        printf("%7.4f%s", a[i], (i == (n-1))?"\n":" ");
        ord = ord->hh.next;
    }
}

void printBetaCoeffs(void){
    extern char *relationTypeStr[];
    extern double betaCoefs[];
    int i;
    //  char *tmp;
    for(i=0; i < NUM_REL; ++i)
        printf("%s\t%lf\n",relationTypeStr[i], betaCoefs[i]);
}


void transposeMatrix(double **a, int n){
    /* leaves the diagonals alone.  */
    int i, j;
    double tmp;
    for (i = 0; i < n; ++i){
        for (j = i + 1; j < n; ++j){
            tmp = a[i][j];
            a[i][j] = a[j][i];
            a[j][i] = tmp;
        }
    }
}

void colNorm(double **a, int n){
    double *sumVec = NULL;
    int i,j;
    sumVec = colSum(a, n);
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            if (sumVec[j] !=0)
                a[i][j] /= sumVec[j];
    free(sumVec);
}

double* colSum(double **a, int n){
    /* a is a ptr to a matrix of size n */
    double *b = NULL;
    int i,j;
    b = malloc(n * sizeof(double));
    assert(b != NULL);
    for (i = 0; i < n; ++i)
        b[i] = 0;
    for (i = 0; i < n; ++i){
        for (j = 0; j < n; ++j){
            b[j] += a[i][j];
        }
    }
    return b;
}

void matAdd(double **a, double **b, int n){
    int i,j;
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            a[i][j] += b[i][j];

}

void matScalMult(double **a, int n, double b){
    int i,j;
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            a[i][j] *= b;
}

double** zeros(int n){
    double *p = NULL;
    double **a = NULL;
    int i, j;
    p = malloc(n * n * sizeof(double));
    a = malloc(n * sizeof(double *));
    //a = calloc(n, sizeof(double *));
    assert(p != NULL);
    assert(a != NULL);
    for (i = 0; i < n; ++i)
        a[i] = p + i * n;
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            a[i][j] = 0;
    /*   for (i = 0; i < n; ++i) */
    /*      a[i] = calloc(n, sizeof(double)); */
    return a;
}

double* zerosVec(int n){
    double *a = NULL;
    //  int i;
    a = calloc(n, sizeof(double));
    assert(a != NULL);
    //  for (i = 0; i < n; ++i)
    //    a[i] = 0;
    return a;
}

void fillDEVec(double *vecNew, int n){
    extern geneItem *geneOrder;
    extern int debug_flag;
    int index[n], i; /* index contains intersected genes' positions in big beta matrix*/
    diffE *de = NULL;
    geneItem *ord   = NULL;
    ord = geneOrder;
    i = 0;
    /* this loop should build the index array
       and build the ordered intersected diff expr vector
    */
    while(ord != NULL){
        de = findDiffExpr(ord->id);
        if(de != NULL){
            vecNew[i]  = de->expr;
            //      addItemToIntersect(de->id, ord->order, i);
            index[i++] = ord->order;
            if(debug_flag)
                fprintf(stderr,"\t[%d]-%s",ord->order,de->id);
        }else{
            vecNew[i] = 0.0;
            //      addItemToIntersect(ord->id, ord->order, i);
            index[i++] = ord->order;
            if(debug_flag)
                fprintf(stderr,"\t[%d]-NA(%s)",ord->order,ord->id);
        }
        ord = ord->hh.next;
    }
    if(debug_flag)
        fprintf(stderr,"\n");
}

int  isMatrixEmpty(double **a, int n){
    int i,j;
    for (i = 0; i < n; ++i){
        for (j = 0; j < n; ++j){
            if(a[i][j] != 0)
                return 0;
        }
    }
    return 1;
}

double determinant(double **a, int n){
    /* Gaussian elimination method to calculate determinant of matrix*/
    /* This is not an effifient way of doing this, but it's simple 
     * enough for our needs. */
    int i, j, k, m, orig, cnt = 0;
    double **b, ans, pivot, multiplier, absmax, *tmp = NULL, *row = NULL;
    // make a copy of a in b, since we're going to be adjusting the contents
    b = zeros(n);
    copyMatrix(a, b, n);

    for(k = 0; k < n-1; ++k){
        // find the pivot element
        absmax = b[k][k] * b[k][k];
        for(m = k + 1; m < n; ++m){
            if((b[m][k] * b[m][k]) > absmax){
                row = b[m];
                orig = m;
                absmax = b[m][k];
            }
        }
        if(row != NULL){
            ++cnt;
            tmp = b[k];
            b[k] = row;
            b[orig] = tmp;
        }
        if(b[k][k] == 0.0){
            return 0.0;
        }
        //
        pivot = b[k][k];
        for(i = k + 1; i < n; ++i){
            multiplier = -b[i][k] * (1.0 / pivot);
            for (j = k; j < n; ++j){
                b[i][j] += multiplier * b[k][j];
            }
        }
    }
    if(cnt % 2 == 0)
        ans = 1.0;
    else
        ans = -1.0;
    for (i = 0; i < n; ++i){
        ans *= b[i][i];
    }
    return ans;
}


void subtractIdent(double **a, int n){
    /* 
       we're subtracting the beta pathway matrix
       from the identity matrix of the same degree
    */
    int i,j;
    for (i=0; i < n; ++i){
        for (j = 0; j < n; ++j){
            if(i == j)
                a[i][j] = 1.0 - a[i][j];
            else
                a[i][j] = 0.0 - a[i][j];
        }
    }
}

void invert(double **a, int n){
    transposeMatrix(a, n);
    /* now we come to the heart of the matter... First get the pivot table, ipiv,
       then we'll calculate the inverse, then it's high-fives and root-beers.
    */
    double **workspace;
    int lda, lworkspace, info, ipiv[n];
    lda = lworkspace = n;
    info = 0;
    workspace = zeros(n);
    f77_dgetrf(&n, &n, *a, &lda, ipiv, &info);
    if(info != 0){
        printf("f77_dgetrf status:[%d]\n", info);
        exit(1);
    }
    f77_dgetri(&n, *a, &lda, ipiv, *workspace, &lworkspace, &info);
    if(info !=0){
        printf("f77_dgetri status:[%d]\n", info);
        exit(1);
    }
    transposeMatrix(a, n);
}

void matVecMultiply(double **a, int n, double *b, double *c)
{
    int i, j;
    double sum = 0;
    for (i = 0; i < n; ++i){
        sum = 0;
        for (j = 0 ; j < n; ++j){
            sum += a[i][j] * b[j];
        }
        c[i] = sum;
    }
}

void matMatMultiply(double **a, int n, double **b, double **c)
{
    int i, j, m;
    double sum = 0;
    for (i = 0; i < n; ++i){
        for (j = 0 ; j < n; ++j){
            sum = 0;
            for (m = 0; m < n; ++m){
                sum += a[i][m] * b[m][j];
            }
            c[i][j] = sum;
        }
    }
}

double sumVec(double *a, int n)
{
    int i;
    double sum = 0.0;
    for (i = 0; i < n; ++i)
        sum += a[i];
    return sum;
}

int compare_dbls(const void *a, const void *b){
    // qsort(array, size, sizeof(double), compare_dbls);
    const double *da = (const double *) a;
    const double *db = (const double *) b;
    return (*da > *db) - (*da < *db);
}
int compare_dbls_rev(const void *a, const void *b){
    // qsort(array, size, sizeof(double), compare_dbls_rev);
    const double *da = (const double *) a;
    const double *db = (const double *) b;
    return (*da < *db) - (*da > *db);
}

double correctTA(double *a, int n){
    // finds the median of the pre-sorted TA array (*a),
    // and subtracts the median from all values in the array.
    int i;
    double med;
    if((n - 1) % 2){
        i = n / 2;
        med = (a[i] + a[i - 1]) / 2.0;
        //      printf("it was even, median = %f\n", med);
    }else{
        i = (n - 1)/2;
        med = a[i];
        //      printf("it was odd, median = %f\n", med);
    }
    for(i = 0; i < n; ++i)
        a[i] -= med;
    //  printf("med = %f\n",med);
    return med;
}

void copyMatrix(double **a, double **b, int n){
    // The matrices may not necessarily be in contiguous
    // chunks, which makes the use of memcpy not ideal.
    int i,j;
    for(i = 0; i < n; ++i){
        for(j = 0; j < n; ++j){
            b[i][j] = a[i][j];
        }
    }
}

void solveForPF(double *a, double *b, double *c, int n){
    int i;
    for(i = 0; i < n; ++i)
        c[i] = b[i] + a[i];
}
