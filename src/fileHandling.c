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
// #define _GNU_SOURCE // getline()
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <getopt.h>    /* get command line options */
#include <sys/types.h>
#include <sys/stat.h>  /* gcc test for regular file types when opening things */
#include <assert.h>
#include "spia.h"

/*
  Functions to read directories, parse tab delimited files, populate hashes and lists.
  as this is written, this is like the main() function, calling the others.
*/

upstreamGene *pathway_up = NULL; /* our pathway hash, upstream genes as keys */
allGene *pathway_all = NULL; /* our pathway_all hash, all pathway genes as keys */
geneItem *geneOrder = NULL; /* our hash to report proper gene order (for matrix access)  */
diffE *diffGeneExp = NULL; /* our hash to store differential gene expression  */
allGene *allGenesTested = NULL; /* our hash to store the list of all genes tested*/
diffE *bootGenes = NULL; /* our hash to store a bootstrap set of genes */
double betaCoefs[] = {1., 0., 0., 1., -1.,
                      1., 0., 0., -1., -1.,
                      0., 0., 1., 0.,  1.,
                      -1., 0. , 1., -1., -1.,
                      0., 0., 0., 0.,  1.,
                      -1., 1., -1.};
/* check spia.h and the enum relationType for context of these values. */
char *relationTypeStr[48] = {"activation", "compound", "binding_association", "expression",
                             "inhibition", "activation_phosphorylation", "phosphorylation",
                             "indirect", "inhibition_phosphorylation", "dephosphorylation_inhibition",
                             "dissociation", "dephosphorylation", "activation_dephosphorylation",
                             "state", "activation_indirect", "inhibition_ubiquination", "ubiquination",
                             "expression_indirect", "indirect_inhibition", "repression",
                             "binding_association_phosphorylation", "dissociation_phosphorylation",
                             "indirect_phosphorylation", "family_membership",
                             "transcriptional_activation", "transcriptional_inhibition",
                             "process_activation", "process_inhibition"};

double *all_de_values = NULL; /* to make boot straps faster when grabbing random de values*/
char *all_pathway_ids[MAX_PATHWAY]; /* to make boot straps faster when grabbing random pathway genes*/
double **beta2 = NULL;
double probNDE = -1.0;
pGlobal *pGlist  = NULL; /* our hash to store all the global p values for all pathways tested*/ 

void gatherOptions(int argc, char **argv, char **oldPathFormatDir, char **de, 
                   char **ar, char **newPathFormatDir, char **betaFile){
    /* 
       gather command line options and send them back to main
    */
    extern int debug_flag;
    extern int verbose_flag;
    extern int showNetAcc_flag;
    extern int quietNetAcc_flag;
    extern int nBoots;
    int c;
    while (1)
        {
            static struct option long_options[] =
                {
                    {"debug", no_argument, &debug_flag, 1},
                    {"verbose", no_argument, NULL, 'v'},
                    {"printNetAcc", no_argument, &showNetAcc_flag, 1},
                    {"quietNetAcc", no_argument, &quietNetAcc_flag, 1},
                    {"help",  no_argument, 0, 'h'},
                    /* These options don't set a flag.
                       We distinguish them by their indices. */
                    {"dir",  required_argument, 0, 'd'},
                    {"de", required_argument, 0, 'e'},
                    {"nBoots", required_argument, 0,'b'},
                    {"array", required_argument, 0,'a'},
                    {"pathFiles", required_argument, 0, 'p'},
                    {"betaCoFile", required_argument, 0, 'c'},
                    {0, 0, 0, 0}
                };
            /* getopt_long stores the option index here. */
            int option_index = 0;
            c = getopt_long(argc, argv, "d:e:b:a:p:c:vh",
                            long_options, &option_index);
            /* Detect the end of the options. */
            if (c == -1){
                break;
            }
            switch (c)
                {
                case 0:
                    break;
                case 'd':
                    *oldPathFormatDir = optarg;
                    break;
                case 'e':
                    *de = optarg;
                    break;
                case 'c':
                    *betaFile = optarg;
                    break;
                case 'b':
                    sscanf(optarg, "%d", &nBoots);
                    break;
                case 'a':
                    *ar = optarg;
                    break;
                case 'p':
                    *newPathFormatDir = optarg;
                    break;
                case 'v':
                    verbose_flag++;
                    break;
                case 'h':
                case '?':
                    usage();
                    /* getopt_long already printed an error message. */
                    break;
                default:
                    abort ();
                }
        }
    if((*oldPathFormatDir == NULL && *newPathFormatDir == NULL) || *de == NULL || *ar == NULL)
        usage();
}


size_t de_getline(char **s, size_t *n, FILE *f) {
    register size_t nMinus1 = ((*n) - 1), i = 0;
    char *s2 = *s;

    while(TRUE) {
        register size_t ch = (char) getc(f);

        if(ch == '\r') {
            ch = getc(f);
        }

        if(i == nMinus1) {
            *n = 2 * (*n) + 1;
            *s = realloc(*s, (*n + 1) * sizeof(char));
            assert(*s != NULL);
            s2 = *s + i;
            nMinus1 = ((*n) - 1);
        }

        if((ch == '\n') || (ch == EOF)) {
            *s2 = '\0';
            return(feof(f) ? 0 : i);
        }else{
            *s2 = ch;
            s2++;
        }
        ++i;
    }
}



int endsIn_tab(char *filename){
    /* endsIn should return a 0 if the name does not end
       in `.tab', and a 1 if it does */
    char nameCopy[MAX_PATH_LENGTH];
    strcpy(nameCopy, filename);
    char *tmp;
    int status = 0;
    /* We'll use strtok to cut up the filename and see if .tab is on the end  */
    for(tmp=strtok(nameCopy,"."); tmp != NULL; tmp = strtok(NULL, ".")){
        if(strcmp(tmp, "tab") == 0)
            status = 1;
    }
    return status;
}

void usage(void){
    fprintf(stderr, "Usage: --dir <pathway directory> --de <Diff Exp File> "
            "--array <Entire Test Set File>  --nBoots <int> --quietNetAcc [options]\n\n");
    fprintf(stderr, "TERMS\n");
    fprintf(stderr, "NDE: Number of differentially expressend genes in pathway.\n");
    fprintf(stderr, "Acc: net perturbation acumulation at a gene.\n");
    fprintf(stderr, "PF: perturbation factor at a gene.\n");
    fprintf(stderr, "t_A: total net accumulated perturbation in the paythway.\n");
    fprintf(stderr, "pPERT: The probability that the total accumulated perturbation of "
            "the pathway, as a random variable of the observed tA is greater than the "
            "observed. i.e. Pr(TA >= t_A | H0). Calculated via bootstrapping.\n");
    fprintf(stderr, "pNDE: Probability based on a standard overrepresentation analysis "
            "using the hypergeometric distribution.\n");
    fprintf(stderr, "pGlobal: Combined probabilities of P_NDE and P_PERT. See Tarca et "
            "al. supplemental for details.\n");
    exit(2);
}

void readDETab(char *filename){
    /* readDETab should be able to take a pathname to a tab file,
       read out the good bits, and then stuff them in a hash.*/
    FILE *ifp = NULL;
    
    size_t nbytes = 512;
    size_t bytes_read = 1;
    char *line = (char *) daemalloc((nbytes + 1));

    extern double *all_de_values;
    extern int debug_flag;
    extern int verbose_flag;
    int i = 0;
    char *id = (char *) daemalloc(MAX_ID_LENGTH + 1);
    double de;
    diffE *g;
    int nArgs;
    ifp = fopen(filename, "r");
    verbose("Reading `%s'\n", filename);
    if(ifp == NULL){
        fprintf(stderr, "ERROR, unable to open `%s', is path correct?\n", filename);
        exit(1);
    }
    while(bytes_read > 0){
        bytes_read = de_getline(&line, &nbytes, ifp);
        nArgs = sscanf(line, "%s\t%lf", id, &de);
        if(nArgs < 2)
            nArgs = sscanf(line, "%s %lf", id, &de);
        g = findDiffExpr(id);
        if(nArgs == 2){
            if(g != NULL){
                addDiffExprsGeneEntry(id, de);
            }else{
                addDiffExprsGene(id, de);
            }
        }
    }
    fclose(ifp);
    verbose("Read complete. Populated hash, now creating array with DE values.\n");
    // now that all of the de_genes and values are read in, populate our
    // extern double array, 
    int n;
    diffE *d;
    d = diffGeneExp;
    n = HASH_COUNT(diffGeneExp);
    verbose("Verbose: DE hash count: %d\n", n);
    all_de_values = zerosVec(n);
    for(i = 0; i < n; ++i){
        all_de_values[i] = d->expr;
        d = d->hh.next;
    }
    free(line);
}

void readArrayTab(char *filename){
    /* readAllTab should be able to take a pathname to a tab file,
       read out the good bits, and then stuff them in a hash.
       reads in the array file.
    */
    FILE *ifp = NULL;
    size_t nbytes = 255;
    size_t bytes_read = 1;
    char *line = (char *) daemalloc((nbytes + 1));
    extern int debug_flag;
    extern int verbose_flag;
    char *id = (char *) daemalloc(MAX_ID_LENGTH + 1);
    allGene *g;
    int nArgs,i;
    ifp = fopen(filename, "r");
    verbose("Reading Array file:`%s'\n", filename);
    if(ifp == NULL){
        fprintf(stderr, "ERROR, unable to open `%s', is path correct?\n", filename);
        exit(1);
    }
    while(bytes_read > 0){
        bytes_read = de_getline(&line, &nbytes, ifp);
        nArgs = sscanf(line, "%d\t%s", &i, id);
        // fprintf(stdout, "%s\n",id);
        g = findAllGene(id);
        if(nArgs == 2){
            if(g == NULL){
                addAllGeneEntry(id);
            }
        }
    }
    int n;
    allGene *d;
    n = HASH_COUNT(d);
    verbose("array hash count: %d\n", n);
    fclose(ifp);
    free(line);
}

int readOldPathway(char *filename){
    /* readOldPathwayTab should be able to take a pathname to a tab file,
       read out the good bits, and then stuff them in a hash.*/
    extern int debug_flag;
    extern int verbose_flag;
    extern int showNetAcc_flag;
    extern int quietNetAcc_flag;
    extern int nBoots;
    extern geneItem *geneOrder;
    FILE *ifp = NULL;
    size_t nbytes = 512;
    size_t bytes_read = 1;
    char *ups, *downs, *pathname, *relType, *relName, *relSymb, *descrip;
    char *line = (char *) daemalloc((nbytes + 1));

    upstreamGene *g;
    relationType enumRelType;
    int nArgs;
    ups = (char *) daemalloc(nbytes + 1);
    downs = (char *) daemalloc(nbytes + 1);
    pathname = (char *) daemalloc(nbytes + 1);
    relType = (char *) daemalloc(nbytes + 1);
    relName = (char *) daemalloc(nbytes + 1);
    relSymb = (char *) daemalloc(nbytes + 1);
    descrip = (char *) daemalloc(nbytes + 1);
    ifp = fopen(filename, "r");
    if(ifp == NULL){
        fprintf(stderr, "ERROR, unable to open `%s', is path correct?\n", filename);
        exit(1);
    }
    verbose("READING `%s'\n", filename);
    while(bytes_read > 0){
        bytes_read = de_getline(&line, &nbytes, ifp);
        if(bytes_read <= 0)
            continue;
        nArgs = sscanf(line, "hsa:%s hsa:%s %s %s %s path:%s %s", ups, downs, relType, 
                       relName, relSymb, pathname, descrip);
        g = findGenePath(ups);
        if(isRelationship(relName, &enumRelType) && (nArgs == 7)){
            addGenePathAll(ups);
            addGenePathAll(downs);
            if(g != NULL){
                verbose("adding %s to existing gene %s\n", downs, ups);
                addInteraction(ups, downs, enumRelType);
            }else{
                verbose("adding new gene %s\n", ups);
                addGenePath(ups);
                addInteraction(ups, downs, enumRelType);
            }
        }
    }
    verbose("%d items stored in the geneOrder hash after reading %s.\n",
            HASH_COUNT(geneOrder), filename);
    fclose(ifp);
    free(line);
    free(ups);
    free(downs);
    free(pathname);
    free(relType);
    free(relName);
    free(relSymb);
    free(descrip);
    /*** PATHWAY STORED IN HASHES, BEGIN POST PROCESSING ***/
    return HASH_COUNT(geneOrder);
}

int readNewPathway(char *filename){
    /* readNewPathwayTab should be able to take a pathname to a tab file,
       read out the good bits, and then stuff them in a hash.*/
    extern int debug_flag;
    extern int verbose_flag;
    extern int showNetAcc_flag;
    extern int quietNetAcc_flag;
    extern int nBoots;
    extern geneItem *geneOrder;
    FILE *ifp = NULL;
    size_t nbytes = 200;
    size_t bytes_read = 1;
    char *itemA, *itemB, *interact;
    char *line = (char *) daemalloc((nbytes + 1));
    upstreamGene *g;
    relationType enumRelType;
    int nArgs;
    line = (char *) daemalloc(nbytes + 1);
    itemA = (char *) daemalloc(nbytes + 1);
    itemB = (char *) daemalloc(nbytes + 1);
    interact = (char *) daemalloc(nbytes + 1);
    ifp = fopen(filename, "r");
    verbose("Reading new pathway format file `%s'\n", filename);
    if(ifp == NULL){
        fprintf(stderr, "ERROR, unable to open `%s', is path correct?\n", filename);
        exit(1);
    }
    while(bytes_read > 0){
        bytes_read = de_getline(&line, &nbytes, ifp);
        if(bytes_read <= 0)
            continue;
        nArgs = sscanf(line, "%s\t%s\t%s", itemA, itemB, interact);
        g = findGenePath(itemA);
        if(isRelationship(interact, &enumRelType) && (nArgs == 3)){
            addGenePathAll(itemA);
            addGenePathAll(itemB);
            if(g != NULL){
                debug("adding interaction %s (new) on %s,\n", itemA, itemB);
                addInteraction(itemA, itemB, enumRelType);
            }else{
                debug("adding interaction %s on %s,\n", itemA, itemB);
                addGenePath(itemA);
                addInteraction(itemA, itemB, enumRelType);
            }
        }
    }
    fclose(ifp);
    free(line);
    free(itemA);
    free(itemB);
    free(interact);
    verbose("%d items stored in the geneOrder hash after reading %s.\n",
            HASH_COUNT(geneOrder), filename);
    /*** PATHWAY STORED IN HASHES, BEGIN POST PROCESSING ***/
    return HASH_COUNT(geneOrder);
}

void readBetaCoeffFile(char *filename){
    extern int debug_flag;
    extern int verbose_flag;
    // extern char *relationTypeStr[];
    extern double betaCoefs[];
    FILE *ifp = NULL;
    size_t nbytes = 200;
    size_t bytes_read = 1;
    char *line = (char *) daemalloc((nbytes + 1) * sizeof(char));
    char *relName;
    relationType enumRelType;
    int nArgs;
    relName  = (char *) daemalloc(nbytes + 1);
    double beta;
  
    ifp = fopen(filename, "r");
    verbose("Reading beta coef file `%s'\n", filename);
    if(ifp == NULL){
        fprintf(stderr, "ERROR, unable to open `%s', is path correct?\n", filename);
        exit(1);
    }
    // printBetaCoeffs();
    while(bytes_read > 0){
        bytes_read = de_getline(&line, &nbytes, ifp);
        if(bytes_read <= 0)
            continue;
        nArgs = sscanf(line, "%s\t%lf", relName, &beta);
        // printf("I see %s, %lf")
        if(isRelationship(relName, &enumRelType) && (nArgs == 2)){
            // char *tmp = relationTypeStr[enumRelType];
            // printf("I see you want to use %s, [%d] with value %lf\n", tmp, enumRelType, beta);
            betaCoefs[enumRelType] = beta;
            if((beta > 1.0) || (beta < -1.0))
                fprintf(stderr, "Warning, your beta coefficient for %s, %e, is not in "
                        "[-1.0 < x < 1.0] which may yield beta matrix singularity "
                        "(i.e. determinant = 0) which would prevent matrix inversion.\n", 
                        relName, beta);
        }
    }
    //  printBetaCoeffs();
    fclose(ifp);
    free(line);
    free(relName);
}

double processPathway(int *status){
    extern int debug_flag;
    extern int verbose_flag;
    extern int showNetAcc_flag;
    extern int quietNetAcc_flag;
    extern int nBoots;
    extern allGene *pathway_all;
    extern geneItem *geneOrder;
    extern diffE *diffGeneExp;
    extern allGene *allGenesTested;
    extern double *all_de_values;
    extern double **beta2;
    extern double probNDE;
    if(beta2 != NULL){
        free(beta2);
        beta2 = NULL;
    }
    unsigned int all_genes_tested;
    unsigned int pathSize;
    unsigned int all_de;
    all_genes_tested = HASH_COUNT(allGenesTested);
    pathSize = HASH_COUNT(pathway_all); //changed this from pathway_de. 17 June 2009 dae
    all_de = HASH_COUNT(diffGeneExp);
    sort_by_idPath();
    sort_by_idOrder();

    /*** HERE IS WHERE THE BETA MATRIX WAS BEING SHRUNK ***/
    unsigned int szMat;
    szMat = HASH_COUNT(geneOrder);
    double **beta = NULL;
    beta = zeros(szMat);
    double **tmp = NULL;
    int i;
    for (i = 0; i < NUM_REL; ++i){
        tmp = buildBeta((relationType) i);
        colNorm(tmp, szMat);
        matScalMult(tmp, szMat, betaCoefs[i]);
        matAdd(beta, tmp, szMat);
    }
    if(debug_flag){
        debug("beta matrix = \n");
        printMatrix(beta, szMat);
    }
    int Nde = countIntersect_de_path(); // number of diff exp genes on pathway
    verbose("There are %d intersections between genes in the pathway "
            "and your DE genes\n", Nde);

    double *deVec = zerosVec(szMat);
    fillDEVec(deVec, szMat);
    if(isMatrixEmpty(beta, szMat)){
        // fprintf(stderr, "Oh noes! Beta matrix is empty!\n");
        free(beta);
        free(tmp);
        free(deVec);
        *status = -1;
        return -1.0;
    }
    double **betaOrig;
    betaOrig = zeros(szMat);
    copyMatrix(beta, betaOrig, szMat);
    subtractIdent(beta, szMat);
    if(debug_flag){
        debug("I - B = \n");
        printMatrix(beta,szMat);
    }
    double det = determinant(beta, szMat); // check for matrix singularity
    if(abs(det) < 1e-7){
        verbose("abs(Determinant), %e, is less than 1e-7.\n", det);
        if(verbose_flag){
            printMatrix(beta, szMat);
        }
        free(beta);
        free(betaOrig);
        free(tmp);
        free(deVec);
        *status = -2;
        return -1.0;
    }
    invert(beta, szMat);
    if(debug_flag){
        debug("inv(I - B) = \n");
        printMatrix(beta,szMat);
    }
    //double **beta2;
    beta2 = zeros(szMat);
    matMatMultiply(betaOrig, szMat, beta, beta2);
    if(debug_flag){
        debug("B * (inv(I - B)) = \n");
        printMatrix(beta2,szMat);
    }

    double *netAcc = zerosVec(szMat);
    double sumNetAcc;
    assert(netAcc != NULL);
    matVecMultiply(beta2, szMat, deVec, netAcc);
    double *pertFact = zerosVec(szMat);
    solveForPF(deVec, netAcc, pertFact, szMat);
    if(!quietNetAcc_flag){
        printf("%6s = ", "Acc");
        // printVector(netAcc, szMat);
        printNamedVector(netAcc, szMat);
        printf("%6s = ", "PF");
        // printVector(pertFact, szMat);
        printNamedVector(pertFact, szMat);
    }
    sumNetAcc = sumVec(netAcc, szMat);
    printf("%6s = %d\n", "pSize", pathSize);
    printf("%6s = %d\n", "NDE", Nde);
    printf("%6s = %f\n", "t_A", sumNetAcc);

    /* PROBABILITY OF SEEING THIS MANY DIFF EXPR GENES TESTING VIA 
     * THE HYPER GEOMETRIC
    */  

    int pathArrayIntersect = countIntersect_array_path();
    verbose("1 - hygecdf(x = %d, m = %d, n = %d, k = %d) = pNDE\n", Nde - 1, 
            pathArrayIntersect, all_genes_tested-pathArrayIntersect, all_de);
    if((Nde - 1 < pathArrayIntersect) && (pathArrayIntersect < all_genes_tested) && (Nde - 1 < all_de)){
        double ans;
        ans = 1 - probCDFHyper(Nde - 1, pathArrayIntersect, all_genes_tested-pathArrayIntersect, all_de);  
        printf("%6s = %e\n", "pNDE", ans);
        probNDE = ans;
    }else{
        printf("%6s = NA\n", "pNDE");
        probNDE = -1;
    }
    free(pertFact);
    free(netAcc);
    free(beta);
    free(betaOrig);
    free(tmp);
    free(deVec);
    return sumNetAcc;
}

int  isRelationship(char *rel, relationType *relType_ptr){
    /* 
       Truly, this is an ugly way of accomplishing this.
       This Function tests: does the string represent a known
       relationship type? returns 0 or 1
    */
    if(strcmp(rel, "activation")==0){
        *relType_ptr = activation;
        return 1;
    }else if(strcmp(rel, "compound")==0){
        *relType_ptr = compound;
        return 1;
    }else if(strcmp(rel, "binding/association")==0){
        *relType_ptr = binding_association;
        return 1;
    }else if(strcmp(rel, "expression")==0){
        *relType_ptr = expression;
        return 1;
    }else if(strcmp(rel, "inhibition")==0){
        *relType_ptr = inhibition;
        return 1;
    }else if(strcmp(rel, "activation_phosphorylation")==0){
        *relType_ptr = activation_phosphorylation;
        return 1;
    }else if(strcmp(rel, "phosphorylation")==0){
        *relType_ptr = phosphorylation;
        return 1;
    }else if(strcmp(rel, "indirect")==0){
        *relType_ptr = indirect;
        return 1;
    }else if(strcmp(rel, "inhibition_phosphorylation")==0){
        *relType_ptr = inhibition_phosphorylation;
        return 1;
    }else if(strcmp(rel, "dephosphorylation_inhibition")==0){
        *relType_ptr = dephosphorylation_inhibition;
        return 1;
    }else if(strcmp(rel, "dissociation")==0){
        *relType_ptr = dissociation;
        return 1;
    }else if(strcmp(rel, "dephosphorylation")==0){
        *relType_ptr = dephosphorylation;
        return 1;
    }else if(strcmp(rel, "activation_dephosphorylation")==0){
        *relType_ptr = activation_dephosphorylation;
        return 1;
    }else if(strcmp(rel, "state")==0){
        *relType_ptr = state;
        return 1;
    }else if(strcmp(rel, "activation_indirect")==0){
        *relType_ptr = activation_indirect;
        return 1;
    }else if(strcmp(rel, "inhibition_ubiquination")==0){
        *relType_ptr = inhibition_ubiquination;
        return 1;
    }else if(strcmp(rel, "ubiquination")==0){
        *relType_ptr = ubiquination;
        return 1;
    }else if(strcmp(rel, "expression_indirect")==0){
        *relType_ptr = expression_indirect;
        return 1;
    }else if(strcmp(rel, "indirect_inhibition")==0){
        *relType_ptr = indirect_inhibition;
        return 1;
    }else if(strcmp(rel, "repression")==0){
        *relType_ptr = repression;
        return 1;
    }else if(strcmp(rel, "binding/association_phosphorylation")==0){
        *relType_ptr = binding_association_phosphorylation;
        return 1;
    }else if(strcmp(rel, "dissociation_phosphorylation")==0){
        *relType_ptr = dissociation_phosphorylation;
        return 1;
    }else if(strcmp(rel, "indirect_phosphorylation")==0){
        *relType_ptr = indirect_phosphorylation;
        return 1;
        /* What follows are steve's relation types!  */
    }else if(strcmp(rel, "family_membership")==0){
        *relType_ptr = family_membership;
        return 1;
    }else if(strcmp(rel, "transcriptional_activation")==0){
        *relType_ptr = transcriptional_activation;
        return 1;
    }else if(strcmp(rel, "transcriptional_inhibition")==0){
        *relType_ptr = transcriptional_inhibition;
        return 1;
    }else if(strcmp(rel, "process_activation")==0){
        *relType_ptr = process_activation;
        return 1;
    }else if(strcmp(rel, "process_inhibition")==0){
        *relType_ptr = process_inhibition;
        return 1;
    }else if(strcmp(rel, "-p>")==0){
        *relType_ptr = activation;
        return 1;
    }else if(strcmp(rel, "-p|")==0){
        *relType_ptr = inhibition;
        return 1;
    }else if(strcmp(rel, "component>")==0){
        *relType_ptr = compound;
        return 1;
    }else if(strcmp(rel, "member>")==0){
        *relType_ptr = family_membership;
        return 1;
    }else if(strcmp(rel, "-t>")==0){
        *relType_ptr = transcriptional_activation;
        return 1;
    }else if(strcmp(rel, "-t|")==0){
        *relType_ptr = transcriptional_inhibition;
        return 1;
    }else if(strcmp(rel, "-ap>")==0){
        *relType_ptr = process_activation;
        return 1;
    }else if(strcmp(rel, "-ap|")==0){
        *relType_ptr = process_inhibition;
        return 1;
    }
    return 0;
}
