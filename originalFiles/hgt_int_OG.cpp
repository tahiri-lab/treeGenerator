//================================================================================================
//=  HGT-DETECTION v3.3.1
//=  Modified by : Nadia Tahiri
//=  Date of modification : March 2014
//=  Authors : Nadia Tahiri, Alix Boc and Vladimir Makarenkov
//=  Date : November 2009
//=
//=  Description : This program detect horizontal gene transfer (HGT). As input it takes 2
//=  trees: a species tree and a gene tree. the goal is to transform the species tree
//=  into the gene tree following a transfer scenario. There are 3 criteria : the robinson and
//=  Foulds distance, the least-square criterion and the bipartition distance. We also use the
//=  subtree constraint. With this version we can now perform simulation.
//=
//=	 input   : file with species tree and gene tree in the newick format.
//=			   In case of simulation, the species tree and all the gene trees in the same file in
//=            the phylip format or newick string
//=  output  : a list of HGT and the criteria values for each one.
//=	 options :
//=
//================================================================================================


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <signal.h>

#include <iostream>
using namespace std;
#pragma warning(disable:4996)

#include "structures.h"
#include "utils_tree.cpp"
#include "fonctions.cpp"

#define binaireSpecies 0
#define binaireGene    0 //VM bilo 1
int compteur=0;

void traiterSignal(int sig){
	printf("\nMESSAGE : SEGMENTATION FAULT #%d DETECTED",sig);
	printf("\nUse valgrind or gdb to fix the problem");
	printf("\n");
	exit(-1);
}

//========================================================================================================
//============================================ MAIN ======================================================
//========================================================================================================
void main_hgt(string tree1, string tree2, double *distances){
	double alpha = distances[5];
	struct InputTree SpeciesTree;				    //== initial species tree
	struct InputTree GeneTree;					    //== initial gene tree
	struct InputTree SpeciesTreeReduce;				    //== initial species tree reduit
	struct InputTree GeneTreeReduce;					    //== initial gene tree reduit
	struct InputTree SpeciesTreeRed;			  //== reduced species tree
	struct InputTree GeneTreeRed;				    //== reduced gene tree
	int cpt_hgt,i,j,nb_same_espece, nb_leaves, nbTree=0;
	int min_diff = 0; // difference minimum of species between T1 and nb_same_espece or between T2 and nb_same_espece
	int bootstrap = 0;
	int multigene = 0;
	int nbHgtFound = 0;
	struct CRITERIA aCrit;						       //== struture of all the criteria
	struct Parameters param;

	//== read parameters
	if(readParameters(&param)==-1){
		printf("\nhgt : no options specified, see the README file for more details\n");
		exit(-1);
	}

	rand_bootstrap = param.rand_bootstrap;

	if(strcmp(param.speciesroot,"file") == 0){
		if(!file_exists(param.speciesRootfileLeaves) && !file_exists(param.speciesRootfile)){
			printf("\nhgt : The file %s does not exist",param.speciesRootfileLeaves);
			exit(-1);
		}
	}
	if(strcmp(param.generoot,"file") == 0){
		if(!file_exists(param.geneRootfileLeaves) && !file_exists(param.geneRootfile)){
			printf("\nhgt : The file %s does not exist",param.geneRootfileLeaves);
			exit(-1);
		}
	}

//==============================================================================
//============================= LECTURE DES ARBRES =============================
//==============================================================================
	initInputTree(&SpeciesTree);
	initInputTree(&GeneTree);
	//printf("TEST");
	nb_same_espece = readInputFile(tree1,tree2, param.input,&SpeciesTree,&GeneTree,param.errorFile);

	if(nb_same_espece<0) {
		nb_same_espece=0;
	}

	initInputTree(&SpeciesTreeReduce);
	initInputTree(&GeneTreeReduce);
	//== lecture des matrices ou chaines newick en entree
	if(readInput(SPECIE,param.input,&SpeciesTreeReduce) == -1)
	{
		printf("\nError in species tree\n");
		exit(-1);
	}
	if(readInput(SPECIE,param.input,&GeneTreeReduce) == -1)
	{
		printf("\nError in gene tree\n");
		getchar();
		exit(-1);
	}


//== lecture des matrices ou chaines newick en entree VM
	TrierMatrices(GeneTreeReduce.Input,GeneTreeReduce.SpeciesName,SpeciesTreeReduce.SpeciesName,SpeciesTreeReduce.size);

     //TrierMatrices(SpeciesTreeReduce.ADD,GeneTreeReduce.SpeciesName,SpeciesTreeReduce.SpeciesName,SpeciesTreeReduce.size);
     //TrierMatrices(GeneTreeReduce.ADD,GeneTreeReduce.SpeciesName,SpeciesTreeReduce.SpeciesName,SpeciesTreeReduce.size);
    //printf ("\nSpeciesTreeSize =%d, GeneTreeSize=%d\n",SpeciesTreeReduce.size, GeneTreeReduce.size);
    for(i=1;i<=SpeciesTreeReduce.size;i++){
		for(j=1;j<=SpeciesTreeReduce.size;j++){		
			SpeciesTreeReduce.ADD[i][j] = SpeciesTreeReduce.Input[i][j];
			GeneTreeReduce.ADD[i][j] = GeneTreeReduce.Input[i][j];
		}
	}

	//NJ(SpeciesTreeReduce.Input,SpeciesTreeReduce.ADD,SpeciesTreeReduce.size);
	//NJ(GeneTreeReduce.Input,GeneTreeReduce.ADD,GeneTreeReduce.size);

	//== construction des differentes repr�sentation des arbres (adjacence,aretes,longueur,degre)
	CreateSubStructures(&SpeciesTreeReduce,1,binaireSpecies);
	CreateSubStructures(&GeneTreeReduce,1,binaireGene);

	InitCriteria(&aCrit,SpeciesTreeReduce.size);
	computeCriteria(SpeciesTreeReduce.ADD,GeneTreeReduce.ADD,SpeciesTreeReduce.size,&aCrit,SpeciesTreeReduce.LONGUEUR,SpeciesTreeReduce.ARETE,GeneTreeReduce.LONGUEUR,GeneTreeReduce.ARETE);
	distances[0]=aCrit.RF;
	//printf("RF=%lf\n",distances[0]);
	if(SpeciesTree.size-nb_same_espece>GeneTree.size-nb_same_espece){
		min_diff = GeneTree.size-nb_same_espece;
	}else{
		min_diff = SpeciesTree.size-nb_same_espece;
	}
	min_diff = fabs(min_diff);
	min_diff = min_diff * min_diff;

	nb_leaves = SpeciesTree.size;
	distances[0] = floor(distances[0]*pow(10,3)+0.5)/(1.0*pow(10,3));

	if(nb_same_espece<=3){
		distances[0]=-((alpha+0.000001)*((min(SpeciesTree.size,GeneTree.size)-(1.0*nb_same_espece))/(1.0*min(SpeciesTree.size,GeneTree.size))));
		//distances[0]=0;
		// distances[0]=1.0+((alpha)*((min(SpeciesTree.size,GeneTree.size)-(1.0*nb_same_espece))/(1.0*min(SpeciesTree.size,GeneTree.size))));
	}else{
		distances[0]=distances[0]/((2.0*nb_same_espece)-6.0)+alpha*((min(SpeciesTree.size,GeneTree.size)-(1.0*nb_same_espece))/(1.0*min(SpeciesTree.size,GeneTree.size))); //normalisation par le nombre d'espèces communes entre les deux arbres
	}

	distances[1]=aCrit.LS;
	distances[2]=aCrit.BD;
	distances[3]=nb_same_espece;

	distances[4]=nb_leaves;

	printf("\n RF=%lf \t nb species in common %d",distances[0], nb_same_espece);

	FreeCriteria(&aCrit,SpeciesTreeReduce.size);
	freeInputTree(&SpeciesTree,SpeciesTree.size);
	freeInputTree(&GeneTree,GeneTree.size);
	freeReducedTree(&SpeciesTreeReduce,SpeciesTreeReduce.size);
	freeReducedTree(&GeneTreeReduce,GeneTreeReduce.size);
}
