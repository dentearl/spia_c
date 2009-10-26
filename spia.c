#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "spia.h"

int debug_flag;   /* a nice extern that gets set in gatherOptions()*/
int verbose_flag; /* a nice extern that gets set in gatherOptions()*/
int showNetAcc_flag=1; /* a nice extern that gets set in gatherOptions()*/
int quietNetAcc_flag=0; /* a nice extern that gets set in gatherOptions()*/
int nBoots = 0; /* The number of bootstraps to be run, per pathway */

int main(int argc, char *argv[]){
  extern double probNDE;
  char *dirName = NULL; 
  char *deName  = NULL; /* differential expression, file  */
  char *arrName = NULL; /* list of all genes tested, file*/
  char *pathFormatName = NULL; /* alternate pathway input mode, one pathway file*/
  srand(time(NULL));
  gatherOptions(argc, argv, &dirName, &deName, &arrName, &pathFormatName);
  if(verbose_flag)
    fprintf(stdout,"Opening differetially expressed (DE) genes file: `%s'\n",deName);
  readDETab(deName);
  readArrayTab(arrName);
  if(dirName){
    if(verbose_flag)
      fprintf(stdout,"Opening pathway dir: `%s'\n",dirName);
    DIR *dip;
    struct dirent *dit;
    int i = 0;
    int leng = 0;
    char fullPath[MAX_PATH_LENGTH];
    char tmpPath[MAX_PATH_LENGTH];
    struct stat buff;
    strcpy(fullPath, dirName);
    if((dip = opendir(dirName))==NULL){
      fprintf(stderr,"Error while openning directory `%s'!\n",dirName);
      exit(-1);
    }
    double tA; // this is the sum of pert factors for our observed data.
    double tAc; // median corrected pert factor value.
    double pPERT; // probability of perturbation.
    int szMat; // this is the size (# uniq genes) of the pathway.
    int Nde;   // number of diff expr genes that are in the pathway
	int status = 0;
    while((dit = readdir(dip)) != NULL){
      i++;
      strcpy(tmpPath,fullPath);
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
	  printf("%s\n", dit->d_name); // print name of pathway
	  szMat = readPathway(tmpPath);
	  if(szMat == 0){
		printf("Unable to process pathway, nothing read from pathway file.\n");
		continue;
	  }
	  status = 0;
	  tA = processPathway(&status);
	  if(status != 0){
		if(status == -1)
		  printf("Unable to process pathway, beta empty.\n");
		if(status == -2)
		  printf("Unable to process pathway, beta singularity.\n");
		continue;
	  }
	  if(nBoots){ // if bootstrapping is turned 'on'
		Nde = countIntersect_de_path();
		double *totalAcc; // T_A from Tarca et al.
		totalAcc = zerosVec(nBoots);
		int k;
		if(Nde == 0){
		  printf("Unable to process pathway, Nde == 0.\n");
		}else{
		  for(k = 0; k < nBoots; ++k){
			populateBootGenes(Nde);
			totalAcc[k] = calcBootPerterb();
			deleteAllBootGenes();
		  }
		  qsort(totalAcc, nBoots, sizeof(double), compare_dbls);
		  double median = correctTA(totalAcc, nBoots);
		  tAc = tA - median;
		  printf("t_Ac  = %f\n",tAc);
		  if(tAc >= 0)
			pPERT = pPERTcalc(tAc, totalAcc, nBoots, 1);
		  else
			pPERT = pPERTcalc(tAc, totalAcc, nBoots, 0);
		  printf("pPERT = %e\n",pPERT);
		  double pG = combPValue(pPERT,probNDE);
		  if((probNDE > -1) && pG != -1 ){
			printf("pG    = %e\n", pG);
			int pSize = countIntersect_array_path();
			addPGlobal(pG, dit->d_name, pSize ,Nde, tAc, pPERT, probNDE);
		  }else
			printf("pG    = NaN\n");
		}
	  } // end bootstrapping loop
    } // end file reading while loop
    if(closedir(dip) != 0){
      fprintf(stderr, "Error while closing directory `%s'!\n",dirName);
      exit(1);
    }
  }else if(singlePathFormatName){
    printf("So, you want a single path file, huh?\nLet me get back to you on that.\n");
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
