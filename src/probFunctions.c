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
#include <assert.h>
#include <math.h>
#include "spia.h"

double probCDFHyper(int x, int m, int n, int k){
    // this isn't so bad if we simplify first, lets give this a shot ourselves...
    // x = black chosen, m = black total, n = white total, k = chosen total
    if(x==0)
        return 0;
    double ans=0;
    double sum=0;
    int i;
    for(i=0; i <= x; ++i){
        ans = sigmaLog(1+m-i, m) + sigmaLog(1+n-(k-i), n) +
            sigmaLog(1+k-i, k) - sigmaLog(1, i) - sigmaLog((1+n+m)-k, n+m);
        sum += exp(ans);
    }
    return (sum);
}

double probPDFHyper(int x, int m, int n, int k){
    // x = black chosen, m = black total, n = white total, k = chosen total
    double ans=0;
	ans = sigmaLog(1+m-x, m) + sigmaLog(1+n-(k-x), n) +
        sigmaLog(1+k-x, k) - sigmaLog(1, x) - sigmaLog((1+n+m)-k, n+m);
	ans= exp(ans);
    return (ans);
}

/*** MARKED FOR REMOVAL  ***/
/* double logOfFactorial(int n){ */
/*   double ans = 0; */
/*   int i; */
/*   for(i=1; i<=n; ++i){ */
/* 	ans += log(i); */
/*   } */
/* } */

double sigmaLog(int i, int n){
    /* calcuates \sum_{i}^n \log(i) */
    double ans = 0;
    for( ; i <= n; ++i)
        ans += log(i);
    return ans;
}

double combPValue(double a, double b){
    double k;
    if(a < 0)
        return -1;
    if(b < 0)
        return -1;
    k = a * b;
    k -= k*log(k);
    return k;
}

double* bonferroni(double *a, int b){
    // we're going to multiply b, the number of trials,
    // by each element of a and return the least of these.
    int i;
    double *ans;
    ans = zerosVec(b);
    for(i = 0; i < b; ++i){
        ans[i] = a[i] * b;
        if(ans[i] > 1.0)
            ans[i] = 1.0;
    }
    return(ans);
}

double* falseDR(double *a, int b){
    // THIS IS NOT COMPLETE AND WAS CANCELED.
    // IT IS LEFT HERE IN CASE THE FEATURE IS RE-REQUESTED
    /***
        This is the subsection of the R function p.adjust:
        i <- n:1
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummin( (n / i) * p[o] ))[ro]
    ***/
    int i;
    double *ord, *revOrd;
    ord = zerosVec(b);
    revOrd = zerosVec(b);
    for(i = 0; i < b; ++i){
        ord[i] = a[i];
        revOrd[i] = a[i];
    }
    double *index;
    index = zerosVec(b);
    for(i = 0; i < b; ++i){
        index[i] = b - i;
    }
    return ord;
    //  qsort(ord, size, sizeof(double), compare_dbls);
    //  qsort(revOrd, size, sizeof(double), compare_dbls_rev);

    /*** This is a secord R function, maybe it will help?:
         xofunction [pID,pN] = FDR(p,q)
         % FORMAT pt = FDR(p,q)
         % 
         % p   - vector of p-values
         % q   - False Discovery Rate level
         %
         % pID - p-value threshold based on independence or positive dependence
         % pN  - Nonparametric p-value threshold
         %______________________________________________________________________________
         % @(#)FDR.m1.3 Tom Nichols 02/01/18


         p = sort(p(:));
         V = length(p);
         I = (1:V)';

         cVID = 1;
         cVN = sum(1./(1:V));

         pID = p(max(find(p<=I/V*q/cVID)));
         pN = p(max(find(p<=I/V*q/cVN)));

    ***/
}


double pPERTcalc(double a, double *b, int n, int opt){
    // b comes in as a sorted vector.
    int i;
    i = (n/2) - 1;
    double ans;
    //printVector(b, n);
    if(opt == 1){ // for positive perturbation...
        //	fprintf(stderr,"%d: %f, %f\n",i, a, b[i]);
        if(b[i] < a){
            while((b[i] < a) && (i < n)){
                //fprintf(stderr,"%d: %f %f\n",i, a, b[i]);
                ++i;
            }
            if(i == n)
                return (double) 1 / n / 100; // from original code, gives a bump to best possible score.
            else
                ans =  ( n - i) / (double) n ;
        }else if(b[i] > a){
            while((b[i] > a) && (i > 0)){
                --i;
            }
            ans =  ( n - i) / (double) n ;
        }else{
            while((b[i] > a) && (i > 0)){
                --i;
            }
            ++i;
            ans =  (n - i) / (double) n;
        }
    }else{ // for negative perturbation
	
        if(b[i] < a){
            while((b[i] < a) && (i < n))
                ++i;
            if(i == n)
                return (double) 1 / n / 100;  // from original code, gives a bump to best possible score.
            else
                ans = (i + 1)/ (double) n ;
        }else if(b[i] > a){
            while((b[i] > a) &&(i >= 0))
                --i;
            ans = (i + 1)/ (double) n ;
        }else{
            while((b[i] > a) && (i >= 0))
                --i;
            ++i;
            ans =  (i + 1)/ (double) n ;
        }
    }
    //  fprintf(stderr,"i = %d, n = %d, ans = %f\n", i, n, ans);
    return 2.0 * ans;
}
