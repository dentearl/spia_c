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
#include <stdarg.h> // va_list
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include "spia.h"

int debug_flag; /* an extern that gets set in gatherOptions() */
int verbose_flag; /* an extern that gets set in gatherOptions() */
int showNetAcc_flag = 1; /* an extern that gets set in gatherOptions() */
int quietNetAcc_flag = 0; /* an extern that gets set in gatherOptions() */
int nBoots = 0; /* The number of bootstraps to be run, per pathway */

int main(int argc, char *argv[]){
    extern double probNDE;
    char *oldFormatDir = NULL; 
    char *newPathFormatDir = NULL; /* alternate pathway input mode*/
    char dirName[MAX_PATH_LENGTH];
    char *deName  = NULL; /* differential expression, file  */
    char *arrName = NULL; /* list of all genes tested, file */
    char *betaCoeffFile = NULL; /* user specified beta coefficients */

    srand(time(NULL));
    gatherOptions(argc, argv, &oldFormatDir, &deName, &arrName, &newPathFormatDir, &betaCoeffFile);
    if(betaCoeffFile){
        readBetaCoeffFile(betaCoeffFile);
        //    printBetaCoeffs();
    }
    verbose("Opening differetially expressed (DE) genes file: `%s'\n", deName);

    readDETab(deName);
    verbose("Opening array file: `%s'\n", arrName);
    readArrayTab(arrName);
    /*** Get the pathway dir name, regardless of the format  ***/
    DIR *dip;
    struct dirent *dit;
    if(oldFormatDir){
        verbose("Opening pathway dir (oldFormat): `%s'\n", oldFormatDir);
        strcpy(dirName, oldFormatDir);
    }
    if(newPathFormatDir){
        verbose("Opening pathway dir (newFormat): `%s'\n", newPathFormatDir);
        strcpy(dirName, newPathFormatDir);
    }
    int i = 0;
    int leng = 0;
    char fullPath[MAX_PATH_LENGTH];
    char tmpPath[MAX_PATH_LENGTH];
    struct stat buff;
    strcpy(fullPath, dirName);
    if((dip = opendir(dirName)) == NULL){
        fprintf(stderr, "Error while openning directory `%s'!\n", dirName);
        exit(1);
    }
    double tA; // this is the sum of pert factors for our observed data.
    double tAc; // median corrected pert factor value.
    double pPERT; // probability of perturbation.
    int szMat; // this is the size (# uniq genes) of the pathway.
    int Nde;   // number of diff expr genes that are in the pathway
    int status = 0;
    while((dit = readdir(dip)) != NULL){
        i++;
        strcpy(tmpPath, fullPath);
        strcat(tmpPath, dit->d_name);
        leng = strlen(dit->d_name);
        stat(tmpPath, &buff);
        cleanup();
        if(!S_ISREG(buff.st_mode)){       // Is it a regular file?
            continue; // if not, bail on this file
        }
        if(!endsIn_tab(tmpPath)){
            continue; // if it does not end in .tab, bail on this file
        }
        // printf("%s\n", dit->d_name); // print name of pathway
        if(oldFormatDir)
            szMat = readOldPathway(tmpPath);
        if(newPathFormatDir)
            szMat = readNewPathway(tmpPath);

        if(szMat == 0){
            fprintf(stderr, "Unable to process pathway `%s', nothing read from pathway file.\n", tmpPath);
            continue;
        }
        status = 0;
        tA = processPathway(&status);
        if(status != 0){
            if(status == -1)
                fprintf(stderr, "Unable to process pathway `%s', beta matrix is "
                        "empty (this is a pathway issue).\n", tmpPath);
            if(status == -2)
                fprintf(stderr, "Unable to process pathway `%s', beta matrix is "
                        "singular (this as a pathway issue).\n", tmpPath);
            continue;
        }
        if(nBoots){ // if bootstrapping is turned 'on'
            Nde = countIntersect_de_path();
            double *totalAcc; // T_A from Tarca et al.
            totalAcc = zerosVec(nBoots);
            int k;
            if(Nde == 0){
                fprintf(stderr, "Unable to process pathway `%s', number of "
                        "differentially expressed genes, Nde = 0 (None of the genes "
                        "in your DE file are present in this pathway).\n", tmpPath);
            }else{
                for(k = 0; k < nBoots; ++k){
                    populateBootGenes(Nde);
                    totalAcc[k] = calcBootPerterb();
                    deleteAllBootGenes();
                }
                qsort(totalAcc, nBoots, sizeof(double), compare_dbls);
                double median = correctTA(totalAcc, nBoots);
                tAc = tA - median;
                printf("%6s = %f\n", "t_Ac", tAc);
                if(tAc >= 0.)
                    pPERT = pPERTcalc(tAc, totalAcc, nBoots, 1);
                else
                    pPERT = pPERTcalc(tAc, totalAcc, nBoots, 0);
                printf("%6s = %e\n", "pPERT", pPERT);
                double pG = combPValue(pPERT, probNDE);
                if((probNDE > -1) && pG != -1 ){
                    printf("%6s = %e\n", "pG", pG);
                    int pSize = countIntersect_array_path();
                    addPGlobal(pG, dit->d_name, pSize , Nde, tAc, pPERT, probNDE);
                }else
                    printf("%6s = NaN\n", "pG");
            }
        } // end bootstrapping loop
    } // end file reading while loop
    if(closedir(dip) != 0){
        fprintf(stderr, "Error while closing directory `%s'!\n", dirName);
        exit(1);
    }

    /*** Now, having processed all pathways, we take the results, sort them, 
         perform the appropriate corrections, and report them. ta-da!
    ***/
    extern pGlobal *pGlist;
    if(HASH_COUNT(pGlist) > 0){
        sort_by_pValue();
        bonferrPGlobal();
        fdrPGlobal();
        printPValues();
    }
    exit(0);
}

void* daemalloc(size_t n)
{
    void *i;
    i = malloc(n);
    if (i == NULL){
        fprintf(stderr, "malloc failed on a request for %zu bytes\n", n);
        exit(1);
    }
    return i;
}

void verbose(char const *fmt, ...)
{
    extern int verbose_flag;
    char str[512];
    va_list args;
    va_start(args, fmt);
    if(verbose_flag){
        int n = vsprintf(str, fmt, args);
        if(n >= 512){
            fprintf(stderr, "Error, failure in verbose(), n = %d\n", n);
            exit(1);
        }
        message("Verbose", str, args);
    }
    va_end(args);
}

void debug(char const *fmt, ...)
{
    extern int debug_flag;
    char str[512];
    va_list args;
    va_start(args, fmt);
    if(debug_flag){
        int n = vsprintf(str, fmt, args);
        if(n >= 512){
            fprintf(stderr, "Error, failure in debug(), n = %d\n", n);
            exit(1);
        }
        message("Debug", str, args);
    }
    va_end(args);
}

void message(char const *type, char const *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    if(verbose_flag){
        fprintf(stderr, "%s: ", type);
        vfprintf(stderr, fmt, args);
    }
    va_end(args);
}
