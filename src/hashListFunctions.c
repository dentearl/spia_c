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
#include <string.h>
#include "spia.h"

/*
  Functions that deal with either hashes or with linked lists
*/

int  addGenePath(char *geneID){
    /* front add  */
    extern upstreamGene *pathway_up;
    upstreamGene *g = NULL;
    g = daemalloc(sizeof(upstreamGene));
    strcpy(g->id, geneID);
    g->down = NULL;
    HASH_ADD_STR(pathway_up, id, g);
    return 0;
}

upstreamGene* findGenePath(char *geneID){
    extern upstreamGene *pathway_up;
    upstreamGene *g = NULL;
    HASH_FIND_STR(pathway_up, geneID, g);
    return g;
}

void addGenePathAll(char *geneID){
    extern allGene *pathway_all;
    extern char *all_pathway_ids[];
    int i=0;
    allGene *g = NULL;
    allGene *f = NULL;
    f = findGenePathAll(geneID);
    if(f != NULL)
        return;
    g = daemalloc(sizeof(allGene));
    strcpy(g->id, geneID);
    HASH_ADD_STR(pathway_all, id, g);
    while(all_pathway_ids[i] != NULL)
        ++i;
    all_pathway_ids[i] = calloc(MAX_ID_LENGTH+1, sizeof(char));
    strcpy(all_pathway_ids[i], geneID);
    //  fprintf(stderr,"adding %s @%d.\n",all_pathway_ids[i],i);
  
}
allGene* findGenePathAll(char *geneID){
    extern allGene *pathway_all;
    allGene *g = NULL;
    HASH_FIND_STR(pathway_all, geneID, g);
    return g;
}
void printPathwayAll(void){
    extern allGene *pathway_all;
    allGene *current;
    current = pathway_all;
    while(current != NULL){
        printf("%s\n", current->id);
        current = current->hh.next;
    }
}
void printPathwayAllArray(void){
    extern allGene *pathway_all;
    extern char *all_pathway_ids[];
    int i,n;
    n = HASH_COUNT(pathway_all);
    for(i = 0; i < n; ++i){
        printf("%s\n", all_pathway_ids[i]);
    }
}


int countDowns(char *geneID){
    upstreamGene *g = NULL;
    downList *d = NULL;
    int count=0;
    g = findGenePath(geneID);
    if(g != NULL){
        if(g->down != NULL){
            d = g->down;
            while(d != NULL){
                ++count;
                d=d->next;
            }
        }
    }
    return count;
}

int addInteraction(char *upID, char *downID, relationType rel){
    // extern geneItem *geneOrder;
    upstreamGene *g = NULL;
    downList     *d = NULL;
    geneItem     *i = NULL;
    g = findGenePath(upID);
    if(g == NULL){
        return -1;
    }
    d = daemalloc(sizeof(downList));
    strcpy(d->id, downID);
    d->relation = rel;
    d->next = g->down;
    g->down = d;
    /*
      check to see if up or down are unqiue to our geneOrder hash.
      if so, add them, if not, carry on.
    */
    i = findGeneOrder(upID);
    if(i == NULL){
        addGeneOrder(upID);
    }
    i = findGeneOrder(downID);
    if(i == NULL){
        addGeneOrder(downID);
    }
    return 0;
}

void deleteAllPath(void){
    upstreamGene *current = NULL;
    extern upstreamGene *pathway_up;
    while(pathway_up){
        current = pathway_up;
        HASH_DEL(pathway_up, current);
        free(current);
    }
}

void printPathway(void){
    upstreamGene *current;
    extern upstreamGene *pathway_up;
    extern char *relationTypeStr[];
    downList *d = NULL;
    char *tmp;
    current = pathway_up;
    while(current != NULL){
        printf("%s", current->id);
        d = current->down;
        while(d != NULL){
            tmp = relationTypeStr[d->relation];
            printf("\n\t%s [%s]", d->id, tmp);
            d = d->next;
        }
        printf("\n");
        current = current->hh.next;
    }
}

int countIntersect_array_path(void){
    extern geneItem *geneOrder;
    int n = 0;
    geneItem *order = NULL;
    allGene *ag = NULL;
    order = geneOrder;
    while(order != NULL){
        ag = findAllGene(order->id);
        if(ag != NULL){
            ++n;
        }
        order = order->hh.next;
    }
    return n;
}

int id_sortPath(upstreamGene *a, upstreamGene *b){
    return strcmp(a->id, b->id);
    /* return(a->id - b->id);  back when we used int ids*/
}

void sort_by_idPath(void){
    extern upstreamGene *pathway_up;
    HASH_SORT(pathway_up, id_sortPath);
}

geneItem* findGeneOrder(char *geneID){
    extern geneItem *geneOrder;
    geneItem *g = NULL;
    HASH_FIND_STR(geneOrder, geneID, g);
    return g;
}

int addGeneOrder(char *geneID){
    /* front add  */
    extern geneItem *geneOrder;
    geneItem *g = NULL;
    g = daemalloc(sizeof(geneItem)); // changed from upstreamGene
    strcpy(g->id, geneID);
    g->order =-1;
    HASH_ADD_STR(geneOrder, id, g);
    return 0;
}

void deleteAllOrder(void){
    geneItem *current = NULL;
    extern geneItem *geneOrder;
    while(geneOrder){
        current = geneOrder;
        HASH_DEL(geneOrder, current);
        free(current);
    }
}

int id_sortOrder(geneItem *a, geneItem *b){
    return strcmp(a->id, b->id);
}

void sort_by_idOrder(void){
    extern geneItem *geneOrder;
    HASH_SORT(geneOrder, id_sortOrder);
    geneItem *current = NULL;
    int i=0;
    current = geneOrder;
    while(current != NULL){
        current->order = i++;
        current = current->hh.next;
    }
}

void printOrder(void){
    geneItem *current = NULL;
    extern geneItem *geneOrder;
    current = geneOrder;
    while(current != NULL){
        /*printf("%2d %s\n", current->order, current->id);*/
        printf("%d:%s%s",current->order,current->id, (current->hh.next==NULL)?"\n":" ");
        current = current->hh.next;
    }
}


int addDiffExprsGene(char *geneID, double expr){
    /* front add  */
    extern diffE *diffGeneExp;
    diffE *d = NULL;
    d = daemalloc(sizeof(diffE));
    strcpy(d->id, geneID);
    d->expr = expr;
    d->count = 1;
    HASH_ADD_STR(diffGeneExp, id, d);
    return 0;
}

int addDiffExprsGeneEntry(char *geneID, double expr){
    /* 
       WE TAKE THE MEAN DIFF EXPRS FOR MULTIPLE
       ENTRIES WITH THE SAME ENTREZ GENE ID.
    */
    extern diffE *diffGeneExp;
    diffE *d = NULL;
    int c;
    double e;
    d = findDiffExpr(geneID);
    assert(d != NULL);
    e = d->expr;
    c = d->count;
    d->expr = (expr - e) / (double) c;
    d->count= ++c;
    return 0;
}

diffE* findDiffExpr(char *geneID){
    extern diffE *diffGeneExp;
    diffE *d;
    HASH_FIND_STR(diffGeneExp, geneID, d);
    return d;
}

void printDiffExpr(void){
    extern diffE *diffGeneExp;
    diffE *current = NULL;
    current = diffGeneExp;
    while(current != NULL){
        printf("[%s] %6.4f\n", current->id, current->expr);
        current = current->hh.next;
    }
}

int countIntersect_de_path(void){
    extern diffE *diffGeneExp;
    extern geneItem *geneOrder;
    int n = 0;
    geneItem *order = NULL;
    diffE *de = NULL;
    order = geneOrder;
    while(order != NULL){
        //	fprintf(stderr,"checking path gene [%s]...", order->id);
        de = findDiffExpr(order->id);
        if(de != NULL){
            ++n;
            //fprintf(stderr,"\t[%s]", de->id);
        }
        order = order->hh.next;
        //fprintf(stderr,"\n");
    }
    return n;
}

/* void addItemToIntersect(char *id, int orgOrd, int newOrd){ */
/*   extern isectList *isect; */
/*   isectList *head = isect; */
/*   isectList *ilist = isect; */
/*   isectList *new = isect; */
/*   int cnt=0; */
/*   new = (isectList *) malloc(sizeof(isectList)); */
/*   assert(new != NULL); */
/*   strcpy(new->id , id); */
/*   new->origOrder = orgOrd; */
/*   new->newOrder  = newOrd; */
/*   new->next = NULL; */
/*   if(head == NULL){ */
/* 	head = new; */
/*   }else{ */
/* 	for( ; ilist->next != NULL; ilist = ilist->next){ */
/* 	  ++cnt; */
/* 	} */
/* 	ilist->next = new; */
/*   } */
/* } */
/* void deleteIntersect(void){ */
/*   printIsectList(); */
/*   fprintf(stderr,"free intersect -- "); */
/*   extern isectList *isect; */
/*   if(isect == NULL) */
/* 	return; */
/*   isectList *cur = isect; */
/*   isectList *nxt = cur->next; */
/*   while(nxt!= NULL){ */
/* 	fprintf(stderr,"."); */
/* 	free(cur); */
/* 	cur = nxt; */
/* 	nxt = nxt->next; */
/*   } */
/*   free(cur); */
/*   fprintf(stderr," success\n"); */
/* } */

/* void printIsectList(void){ */
/*   extern isectList *isect; */
/*   isectList *cur; */
/*   cur = isect; */
/*   fprintf(stderr,"\npreparing to print interesect list:\n"); */
/*   while(cur != NULL){ */
/* 	fprintf(stderr,"[%d %d] %s\n", cur-> newOrder, cur->origOrder, cur->id); */
/* 	cur = cur->next; */
/*   } */
/* } */


int addAllGeneEntry(char *geneID){
    /* front add  */
    extern allGene *allGenesTested;
    allGene *a = NULL;
    a = daemalloc(sizeof(allGene));
    strcpy(a->id, geneID);
    HASH_ADD_STR(allGenesTested, id, a);
    return 0;
}

allGene* findAllGene(char *geneID){
    extern allGene *allGenesTested;
    allGene *a = NULL;
    HASH_FIND_STR(allGenesTested, geneID, a);
    return a;
}

void printAllGene(void){
    extern allGene *allGenesTested;
    allGene *current;
    current = allGenesTested;
    while(current != NULL){
        printf("[%s]\n", current->id);
        current = current->hh.next;
    }
}

void deleteAllBootGenes(void){
    diffE *current = NULL;
    extern diffE *bootGenes;
    //fprintf(stderr,"deleteAllBootGenes -- ");
    while(bootGenes){
        current = bootGenes;
        HASH_DEL(bootGenes, current);
        free(current);
    }
    //fprintf(stderr,"success\n");
}

void populateBootGenes(int n){
    extern diffE *bootGenes;
    diffE *d = NULL;
    int i=0;
    char **gene;
    while(i <n){
        d = NULL;
        gene = randomPathwayGene(); // pick a random pathway gene
        d = findBootGene(*gene);
        if (d == NULL){
            addBootGene(*gene);
            ++i;
        }
    }
}

int addBootGene(char *geneID){
    /* front add  */
    extern diffE *bootGenes;
    diffE *d;
    d = daemalloc(sizeof(diffE));
    strcpy(d->id, geneID);
    d->expr = randomDEvalue(); // random expression value, with resampling
    HASH_ADD_STR(bootGenes, id, d);
    return 0;
}

diffE* findBootGene(char *geneID){
    extern diffE *bootGenes;
    diffE *d = NULL;
    HASH_FIND_STR(bootGenes, geneID, d);
    return d;
}

double randomDEvalue(void){
    extern diffE *diffGeneExp;
    extern double *all_de_values;
    int n = HASH_COUNT(diffGeneExp);
    assert(n != 0);
    int i = rand()%n;
    // printf("DE I picked: %f\n",all_de_values[i]);
    return all_de_values[i];
}

char** randomPathwayGene(void){
    // returns a random char* with the gene ID for one
    // gene in the pathway.
    extern char *all_pathway_ids[];
    extern allGene *pathway_all;
    int n = HASH_COUNT(pathway_all);
    assert(n != 0);
    int i = rand()%n;
    //  fprintf(stderr,"%s\n", all_pathway_ids[i]);
    return &all_pathway_ids[i];
}
void deleteAllPathAll(void){
    extern allGene *pathway_all;
    allGene *current = NULL;
    //  fprintf(stderr,"deleteAllPathAll -- ");
    while(pathway_all){
        current = pathway_all;
        HASH_DEL(pathway_all, current);
        free(current);
    }
    //fprintf(stderr,"success\n");
}

void printBootGenes(void){
    extern diffE *bootGenes;
    diffE *b = NULL;
    b = bootGenes;
    while(b != NULL){
        printf("[%s] %f\n", b->id, b->expr);
        b = b->hh.next;
    }
}

double calcBootPerterb(void){
    extern diffE *bootGenes;
    extern double **beta2;
    extern geneItem *geneOrder;    
    diffE *b = NULL;
    b = bootGenes;
    assert(b != NULL);
    assert(beta2 != NULL);
    int szMat = HASH_COUNT(geneOrder);
    double *netAcc = zerosVec(szMat);
    double sumNetAcc;
    assert(netAcc != NULL);
    double *deVec = zerosVec(szMat);
    fillBootDEVec(deVec, szMat);
    //  printVector(deVec, szMat);
    //  printf("\n");
    matVecMultiply(beta2, szMat, deVec, netAcc);
    double *pertFact = zerosVec(szMat);
    solveForPF(deVec, netAcc, pertFact, szMat);
    sumNetAcc = sumVec(netAcc, szMat); 
    //  fprintf(stdout, "T_A  = %f\n",sumNetAcc);
    return sumNetAcc;
}

void fillBootDEVec(double *vecNew, int n){
    extern diffE *bootGenes;
    extern geneItem *geneOrder;
    int index[n], i; /* index contains intersected genes' positions in big beta matrix*/
    diffE *de = NULL;
    geneItem *ord   = NULL;
    ord = geneOrder;
    i = 0;
    /* this loop should build the index array
       and build the ordered intersected diff expr vector
    */
    while(ord != NULL){
        de = findBootGene(ord->id);
        if(de != NULL){
            vecNew[i]  = de->expr;
            index[i++] = ord->order;
        }else{
            vecNew[i] = 0.0;
            index[i++] = ord->order;
        }
        ord = ord->hh.next;
    }
}

int addPGlobal(double p, char *c, int pathSize, int Nde, double t, double pPERT, double pNDE){
    /* front add  */
    extern pGlobal *pGlist;
    pGlobal *pGlob = NULL;
    pGlob = daemalloc(sizeof(pGlobal));
    strcpy(pGlob->path, c);
    pGlob->p = p;
    pGlob->bonf = -1;
    pGlob->fdr  = -1;
    pGlob->pathSize = pathSize;
    pGlob->tAc  = t;
    pGlob->pPERT = pPERT;
    pGlob->pNDE  = pNDE;
    pGlob->NDE   = Nde;
    HASH_ADD_STR(pGlist, path, pGlob);
    return 0;
}

int value_sort(pGlobal *a, pGlobal *b){
    return (1000000* (a->p - b->p));
}
int value_rev_sort(pGlobal *a, pGlobal *b){
    return (1000000* (b->p - a->p));
}
void sort_by_pValue(void){
    extern pGlobal *pGlist;
    HASH_SORT(pGlist, value_sort);
}
void rev_sort_by_pValue(void){
    extern pGlobal *pGlist;
    HASH_SORT(pGlist, value_rev_sort);
}
void bonferrPGlobal(void){
    extern pGlobal *pGlist;
    pGlobal *cur = NULL;
    cur = pGlist;
    int n = HASH_COUNT(pGlist);
    while(cur != NULL){
        cur->bonf = cur->p * n;
        if(cur->bonf > 1)
            cur->bonf = 1;
        cur = cur->hh.next;
    }
}
void fdrPGlobal(void){
    extern pGlobal *pGlist;
    pGlobal *cur = NULL;
    int i,n = HASH_COUNT(pGlist);
    rev_sort_by_pValue();
    cur = pGlist;
    double cummin = 1.0;
    printf("\n");
    for(i = n;i > 0; --i){
        //	printf("(%d / %d) * %e = %e < %e ...\n ", n, i, cur->p, (n/(double)i)*cur->p, cummin);
        if(((n/(double)i) * cur->p) < cummin){
            //printf("yes, %e < %e\n",(n/(double)i) * cur->p, cummin );
            cummin = (n /(double)i) * cur->p;
        }
        cur->fdr = cummin;
        //	printf("... %e\n",cur->fdr);
        cur = cur->hh.next;
    }
}
void printPValues(void){
    extern pGlobal *pGlist;
    pGlobal *cur = NULL;
    sort_by_pValue();
    cur = pGlist;
    printf("Path\tpSize\tNDE\ttAc\tpNDE\tpPERT\tpGlobal\tBonferr\tFDR\n");
    while(cur != NULL){
        printf("%s\t%d\t%d\t%e\t%e\t%e\t%e",cur->path, cur->pathSize, cur->NDE, cur->tAc, cur->pNDE, cur->pPERT ,cur->p);
        if(cur->bonf > -1)
            printf("\t%e",cur->bonf);
        if(cur->fdr > -1)
            printf("\t%e",cur->fdr);
        printf("\n");
        cur = cur->hh.next;
    }
}
void cleanup(void){
    extern char *all_pathway_ids[];
    int i=0;
    deleteAllPath();
    deleteAllOrder();
    deleteAllPathAll();
    while(all_pathway_ids[i] != NULL){
        free(all_pathway_ids[i]);
        all_pathway_ids[i] = NULL;
        ++i;
    }
}
