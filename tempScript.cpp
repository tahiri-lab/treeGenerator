#include <iostream>
#include <fstream>
#include <vector>
#include <string>
//#include <boost/regex.hpp>

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sstream>
#include <algorithm>
//#include "structures.h"
#include "hgt_int.cpp"


using namespace std;

#define INFINI 999999.99
#define MaxRF 0
double MaxLong = 0;
double MinLong = INFINI;
double seuil;
double epsilona = 0.0005;
int n;
string myreplace(string &s, string toReplace,string replaceWith);
string intToString(int i);

//
// Temporary script to improve / modify the tree generator
//
string modifTree (int nbSpecies, string newTree, int i, int j){
	int nbSpecie1 = 0;
	int nbSpecie2 = 0;
	string specie1 = "";
	string specie2 = "";
	
	do{
		//Generate randomly the number of specie I remove in the species tree (This number must be between 1 and nbSpecies
		nbSpecie1 = rand() % nbSpecies + 1;
		specie1 = intToString(nbSpecie1);
		
		nbSpecie2 = nbSpecie1;
		if(nbSpecie2!=i && nbSpecie2!=j){
			specie1 = specie1 + ":";
		
			nbSpecie2 = nbSpecie2 + nbSpecies;
			specie2 = intToString(nbSpecie2);
			specie2 = specie2 + ":";
		}
		
	}while((newTree.find(specie1)>newTree.length() || newTree.find(specie1)==-1) || nbSpecie2==i || nbSpecie2==j);
	//printf("\n mais on répète un autre truc là non ? %d %d", nbSpecie1, nbSpecie2);
	newTree = myreplace(newTree, specie1 ,specie2);

	return newTree;
}


string intToString(int i){
     ostringstream oss;
     oss << i;
     return oss.str();
}

//=============================================
//
//=============================================
int floor1(double x)
{  
	int i;
  
	if (ceil(x)-floor(x)==2) { i = (int)x; } 
	else if (fabs(x-floor(x)) > fabs(x-ceil(x))) { i = (int)ceil(x); }
	else { i = (int)floor(x); }

	return i;
} 

//=============================================
//
//=============================================
void odp1ct(double **D, int *X, int *i1, int *j1, int n)
{
	double S1, S;
	int i, j, k, a, *Y1;
 
	Y1 = (int *) malloc((n+1)*sizeof(int));
 
	for(i = 1; i <= n; i++)
	{ Y1[i] = 1; }
  
	X[1] = *i1;
	X[n] = *j1;
	if (n==2) { return; }
	Y1[*i1] = 0;
	Y1[*j1] = 0;

	for(i = 0; i <= n-3; i++)
	{
		a = 2;
		S = 0;
		for(j = 1; j <= n; j++)
		{
			if (Y1[j]>0)
			{
				S1 = D[X[n-i]][X[1]] - D[j][X[1]] + D[X[n-i]][j];
				if ((a==2)||(S1<=S))
				{
					S = S1;
					a = 1;
					X[n-i-1] = j;
					k = j;
				}
			}
		}
		Y1[k] = 0;
	}     
   free(Y1);
}


//=============================================
//
//=============================================
void SAVEASNewick2(double *LONGUEUR, long int *ARETE,int nn,const char* t, string& newick) 
{
	int n,root,a;
	int Ns;
	int i, j, sd, sf, *Suc, *Fre, *Tree, *degre, *Mark;
	double *Long;
	int *boot; 
	char *string = (char*)malloc(100000);	
	n = nn;
	Ns = 2*n-3;
	
	double * bootStrap = NULL;
	
	Suc = (int*) malloc((2*n) * sizeof(int));
	Fre = (int*) malloc((2*n) * sizeof(int));
	degre = (int*) malloc((2*n) * sizeof(int));
	Long = (double*) malloc((2*n) * sizeof(double));	
	boot = (int*) malloc((2*n) * sizeof(int));
	Tree = (int*) malloc((2*n) * sizeof(int));
	Mark = (int*) malloc((2*n) * sizeof(int));
	
	if ((degre==NULL)||(Mark==NULL)||(string==NULL)||(Suc==NULL)||(Fre==NULL)||(Long==NULL)||(Tree==NULL)||(ARETE==NULL)||(LONGUEUR==NULL))	
		{ printf("Tree is too large to be saved"); return;} 
	
	for (i = 1; i <= 2*n-3; i++)
	{ 	
		if (i <= n) { degre[i] = 1; }
		else { degre[i] = 3; }
	}
  degre[2*n-2] = 3;
	
	root = Ns+1;
	
	int cpt = 0;
	
	for (;;)
	{
		cpt++;
		if(cpt > 1000000) { exit(1); }
		a = 0;
		a++;
		for (j = 1; j <= 2*n-2; j++)
			{ Mark[j] = 0; }
		
		for (i = 1; i <= 2*n-3; i++)
		{ 	  									
			if ((degre[ARETE[2*i-2]]==1)&&(degre[ARETE[2*i-1]]>1)&&(Mark[ARETE[2*i-1]]==0)&&(Mark[ARETE[2*i-2]]==0))
			{
				Tree[ARETE[2*i-2]] = ARETE[2*i-1];
				degre[ARETE[2*i-1]]--;
				degre[ARETE[2*i-2]]--;
				Mark[ARETE[2*i-1]] = 1;
				Mark[ARETE[2*i-2]] = 1;
				Long[ARETE[2*i-2]] = LONGUEUR[i-1];
				if(bootStrap != NULL) { boot[ARETE[2*i-2]] = (int) bootStrap[i-1]; }
				
			}

			else if ((degre[ARETE[2*i-1]]==1)&&(degre[ARETE[2*i-2]]>1)&&(Mark[ARETE[2*i-1]]==0)&&(Mark[ARETE[2*i-2]]==0))
			{
				Tree[ARETE[2*i-1]] = ARETE[2*i-2];
				degre[ARETE[2*i-1]]--;
				degre[ARETE[2*i-2]]--;
				Mark[ARETE[2*i-1]] = 1;
				Mark[ARETE[2*i-2]] = 1;
				Long[ARETE[2*i-1]] = LONGUEUR[i-1];
				if(bootStrap != NULL) { boot[ARETE[2*i-1]] = (int) bootStrap[i-1]; }
				
			}

			else if ((degre[ARETE[2*i-2]]==1)&&(degre[ARETE[2*i-1]]==1)&&(Mark[ARETE[2*i-2]]==0)&&(Mark[ARETE[2*i-1]]==0))
			{
				Tree[ARETE[2*i-2]] = ARETE[2*i-1];
				root = ARETE[2*i-1];
				degre[ARETE[2*i-1]]--;
				degre[ARETE[2*i-2]]--;
				a =-1;
				Long[ARETE[2*i-2]] = LONGUEUR[i-1];
				if(bootStrap != NULL) { boot[ARETE[2*i-2]] = (int) bootStrap[i-1]; }
			}

			if (a==-1) { break; }
		}
		if (a==-1) { break; }
	}
	
	
	/*  On decale et on complete la structure d'arbre avec Successeurs et Freres  */
	for (i = Ns+1; i > 0; i--)
	{
		Fre[i] = 0;
		Suc[i] = 0;
		//Tree[i] = Tree[i-1]+1;
		//Long[i] = Long[i-1];
	}
	Tree[root] = 0;/*Tree[Ns+1]=0;*/
	
	for (i = 1; i <= Ns+1/*Ns*/; i++)
	{	
		if (i!=root) 
		{
			sd = i;
			sf = Tree[i];

			if (Suc[sf]==0) { Suc[sf] = sd; }
			else
			{	
				sf = Suc[sf];
				while (Fre[sf]>0) { sf = Fre[sf]; }
				Fre[sf] = sd;
			}		 
		}
	}
	
	
	/* On compose la chaine parenthesee */
	string[0] = 0;
	i = root;/*i=Ns+1;*/
	cpt = 0;
	for (;;)
	{	
		if(cpt > 1000000) exit(1);
		
		if (Suc[i]>0)
		{
			sprintf(string,"%s(",string);
			Suc[i] =- Suc[i];
			i =- Suc[i];
		}

		else if (Fre[i]!=0)
		{
			if (Suc[i]==0) { sprintf(string,"%s%d:%.4f,",string,i,Long[i]); }
			else
			{
				if(bootStrap != NULL) { sprintf(string,"%s%d:%.4f,",string,boot[i],Long[i]); }
				else { sprintf(string,"%s:%.4f,",string,Long[i]); }
			}	
			i = Fre[i];
    	}

		else if (Tree[i]!=0)
		{
			if (Suc[i]==0) { sprintf(string,"%s%d:%.4f)",string,i,Long[i]); }
			else
			{
				if(bootStrap != NULL) { sprintf(string,"%s%d:%.4f)",string,boot[i],Long[i]); }
				else { sprintf(string,"%s:%.4f)",string,Long[i]); }
			}
			i = Tree[i];
		}
		else { break; }
	}	

	strcat(string,";");
	//printf("\n %s \n", string);
	newick = string;
	//FILE *pt_t = fopen(t,"a+");
	//fprintf(pt_t,"%s\n",string);
	//fclose(pt_t);
	
	free(Suc);
	free(Fre);
	free(Tree);
	free(Long);
	free(degre);
	free(Mark);
	free(string);
}

//=============================================
//
//=============================================
void tree_generation(double **DA, double **DI, int n, double Sigma)
{
	struct TABLEAU { int V; } **NUM, **A;
	int i, j, k, p, a, a1, a2, *L, *L1, n1;
	double *LON,X0,X,U;
   //time_t t;
   
	n1 = n*(n-1)/2;
    L = (int *)malloc((2*n-2)*sizeof(int));
	L1 = (int *)malloc((2*n-2)*sizeof(int));
	LON = (double *)malloc((2*n-2)*sizeof(double));
    
	NUM = (TABLEAU **)malloc((2*n-2)*sizeof(TABLEAU*));
	A = (TABLEAU **)malloc((n1+1)*sizeof(TABLEAU*));
    
	for (i = 0; i <= n1; i++)
	{
		A[i] = (TABLEAU*)malloc((2*n-2)*sizeof(TABLEAU));
		if (i <= 2*n-3) { NUM[i] = (TABLEAU*)malloc((n+1)*sizeof(TABLEAU)); }
      
		if ((A[i]==NULL)||((i<=2*n-3)&&(NUM[i]==NULL))) 
		{
			printf("\nData matrix is too large\n");
			exit(1);
		}
	}  

	/* Generation of a random additive tree topology T*/

	for (j = 1; j <= 2*n-3; j++)
	{
		for (i = 1; i <= n; i++)
		{
			A[i][j].V = 0;
			NUM[j][i].V = 0;
    	}
    	for (i = n+1; i <= n1; i++)
		{ A[i][j].V = 0; }
	}
	A[1][1].V = 1;
	L[1] = 1;
	L1[1] = 2;
	NUM[1][1].V = 1;
	NUM[1][2].V = 0;
     
   //srand((unsigned) time(&t));

	for (k = 2; k <= n-1; k++)
	{
		p = (rand() % (2*k-3))+1;

		for (i = 1; i <= (n*(k-2)-(k-1)*(k-2)/2+1); i++)
		{ A[i][2*k-2].V = A[i][p].V; }

	    for (i = 1; i <= k; i++)
		{
			a = n*(i-1)-i*(i-1)/2+k+1-i;
			
			if (NUM[p][i].V==0) { A[a][2*k-2].V = 1; }
			else { A[a][p].V = 1; }
		}
		
		for (i = 1; i <= k; i++)
		{
			a = n*(i-1) - i*(i-1)/2 +k+1-i;
			A[a][2*k-1].V = 1;
    	}

    	for (j = 1; j <= k; j++)
    	{
			if (j==L[p])
			{
				for (i = 1; i <= 2*k-3; i++)
				{
					if (i!=p)
					{
						if (L1[p]>L[p]) { a = floor1((n-0.5*L[p]) * (L[p]-1)+L1[p]-L[p]); }
						else { a = floor1((n-0.5*L1[p]) * (L1[p]-1)+L[p]-L1[p]); }

						if (A[a][i].V==1)
						{
							if (NUM[i][L[p]].V==0) { a = floor1((n-0.5*L[p]) * (L[p]-1)+k+1-L[p]); }
							else { a = floor1((n-0.5*L1[p]) * (L1[p]-1)+k+1-L1[p]); }
							A[a][i].V = 1;
						}
					}
				}
			}

			else if (j!=L1[p])
			{
				a = floor1((n-0.5*j)*(j-1)+k+1-j);
        		if (j<L[p]) { a1 = floor1((n-0.5*j)*(j-1)+L[p]-j); }
        		else { a1 = floor1((n-0.5*L[p])*(L[p]-1)+j-L[p]); }

        		if (j<L1[p]) { a2 = floor1((n-0.5*j)*(j-1)+L1[p]-j); }
        		else { a2 = floor1((n-0.5*L1[p])*(L1[p]-1)+j-L1[p]); }
        
        		for(i=1;i<=2*k-3;i++)
        		{
        			if ((i!=p)&&((A[a1][i].V+A[a2][i].V==2)||((NUM[i][j].V+NUM[i][L[p]].V==0)&&(A[a2][i].V==1))||((NUM[i][j].V+NUM[i][L1[p]].V==0)&&(A[a1][i].V==1))))
            		{ A[a][i].V = 1; }
       			}
    		}
    	}

    	for(i = 1; i <= k; i++)
    	{ NUM[2*k-2][i].V = NUM[p][i].V; }

    	NUM[2*k-2][k+1].V = 1;

    	for(i = 1; i <= k; i++)
    	{ NUM[2*k-1][i].V = 1; }

    	for(i = 1; i <= 2*k-3; i++)
    	{
    		if (((NUM[i][L[p]].V+NUM[i][L1[p]].V)!=0)&&(i!=p)) { NUM[i][k+1].V = 1; }
    	}

    	L[2*k-2] = k+1;
    	L1[2*k-2] = L1[p];
    	L[2*k-1] = L1[p];
    	L1[2*k-1] = k+1;
    	L1[p] = k+1;
	}

   //srand((unsigned) time(&t));
   
 // Exponential distribution generator with expectation 1/(2n-3)  
 // 1. Generate U~U(0,1)
 // 2. Compute X = -ln(U)/lambda
 // 3. To obtain an exponential distribution with theoretical mean (= "expected
 // value") of 1/(2n-3), use lambda = 1 and multiply X by 1/(2n-3):     
 // Thus, X = -(1/(2n-3))*ln(U)
 // This formula is given in Numerical Recipes, p.200.*/ 
 // Every branch length of the tree T is then multiplied by 1+aX,
 // where X followed the standard exponential distribution (P(X>n)=exp(-n)) and 
 //  "a" is a tuning factor to adjust the deviation from the molecular clock; "a" is set at 0.8. 
 // 4. LON[i] = LON[i]*(1+0.8*-ln(U))
 
	U = 0.1;
	while (U<0.1) { U = 1.0*rand()/RAND_MAX; }
   
	for(i = 1; i <= 2*n-3; i++)
	{ LON[i] = -1.0/(2*n-3)*log(U); } 
   
   //for (i=1;i<=2*n-3;i++)
	i = 1;
	while (i<=2*n-3)
	{
    	U = 1.0*rand()/RAND_MAX;
    	LON[i] = 1.0*LON[i]*(1.0+0.8*(-log(U)));
    	LON[i] = LON[i]*0.8; ///5.0 ;//0.8
    	if (LON[i]>2*epsilona) { i++; }/*printf("U=%f LON[%d]=%f \n",U,i,LON[i]);*/
	}
      
 // Computation of a tree distance matrix (tree metric matrix)
	for(i = 1; i <= n; i++)
	{
    	DA[i][i] = 0;
    	for(j = i+1; j <= n; j++)
    	{
    		DA[i][j] = 0;
    		a = floor1((n-0.5*i)*(i-1)+j-i);

    		for(k = 1; k <= 2*n-3; k++)
			{
        		if (A[a][k].V==1) { DA[i][j] = DA[i][j]+LON[k]; }
    		}
    		DA[j][i] = DA[i][j];
    	}
	}

 // Normalisation of the tree distance matrix
 /*double m=0.0, var=0.0;
 for (i=1;i<=n;i++)
 {
  for (j=i+1;j<=n;j++)
  { 
   m=m+2.0*DA[i][j]/n/(n-1);
   var=var+2.0*DA[i][j]*DA[i][j]/n/(n-1);
  }
 }
 var=var-m*m;
 for (i=1;i<=n;i++)
 {
  for (j=1;j<=n;j++)
      DA[i][j]=DA[i][j]/sqrt(var); 
 }*/

 	/* Addition of noise to the tree metric */
  //srand((unsigned) time(&t));
  
	for(i = 1; i <= n; i++)
	{
    	DI[i][i] = 0.0;

    	for(j = i+1; j <= n; j++)
    	{
    		X = 0.0;

    		for(k = 1; k <= 5; k++)
    		{
	    		X0 = 1.0*rand()/RAND_MAX;
        		X = X+0.0001*X0;
    		}

    		X = 2*sqrt(0.6)*(X-2.5);
    		U = X-0.01*(3*X-X*X*X);
    		DI[i][j] = DA[i][j]+Sigma*U;

    		if (DI[i][j]<0) { DI[i][j] = 0.01; }

    		DI[j][i] = DI[i][j];
    	}
	}
   
	free (L);
	free(L1);
	free(LON);    
	for(i = 0; i <= n1; i++)
	{
    	free(A[i]);
    	if (i<=2*n-3) { free(NUM[i]); }
	}
	free(NUM);
	free(A);
}

//=============================================
//
//=============================================
void Tree_edges (double **DI, long int *ARETE, double *LONGUEUR, int n)
{ 
	struct EDGE { unsigned int U; unsigned int V; double LN;};
	struct EDGE *Path,*Tree;
	int i,j,k,p,P,*X;
	double S,DIS,DIS1,*L,**D;
 
	X = (int *)malloc((n+1)*sizeof(int));  
	L = (double *)malloc((n+1)*sizeof(double));
	Tree = (EDGE *)malloc((2*n-2)*sizeof(EDGE));
	Path = (EDGE *)malloc((n+2)*sizeof(EDGE));

	D = (double **) malloc((n+1)*sizeof(double*));

	for (i=0;i<=n;i++)
	{
		D[i] = (double*)malloc((n+1)*sizeof(double)); 
  
		if (D[i]==NULL)
		{
			printf("Data matrix is too large"); exit(1); 
		}
	}

	i = 1; j = n;
	odp1ct(DI,X,&i,&j,n);

	for (i = 1; i <= n; i++)
	{ 
		for (j = 1; j <= n; j++) 
			{D[i][j] = DI[i][j];}
	}  
   
 /* Verification de l'equivalence des topologies */
	L[1] = D[X[1]][X[2]];
	Path[1].U = X[1];
	Path[1].V = X[2];
	p = 0;
	P = 1;

	for(k = 2; k <= n-1; k++)
	{
		DIS = (D[X[1]][X[k]]+D[X[1]][X[k+1]]-D[X[k]][X[k+1]])/2;
		DIS1 = (D[X[1]][X[k+1]]+D[X[k]][X[k+1]]-D[X[1]][X[k]])/2;

		S = 0.0;
		i = 0;

		if (DIS>2*epsilona)
		{
    		while (S<DIS-epsilona)
    		{
    			i = i+1;
    			S = S+L[i];
			}
		}
		else { DIS = 0; i = 1; }
  
		Tree[p+1].U = n+k-1;
		Tree[p+1].V = Path[i].V;
		Tree[p+1].LN = S-DIS;
		if (Tree[p+1].LN<2*epsilona) { Tree[p+1].LN = epsilona;  }
  
		for (j = i+1; j <= P; j++)
		{
			Tree[p+j-i+1].U = Path[j].U;
			Tree[p+j-i+1].V = Path[j].V;
			Tree[p+j-i+1].LN = L[j];
			if (L[j]<2*epsilona) { L[j] = 2*epsilona; }
		}
		
		p = p+P-i+1;
		Path[i].V = n+k-1;
		Path[i+1].U = n+k-1;
		Path[i+1].V = X[k+1];
		L[i] = L[i]+DIS-S;
		L[i+1] = DIS1;
		P = i+1;
	}
 
	for (i = 1; i <= P; i++) 
	{
		Tree[p+i].U = Path[i].U;
		Tree[p+i].V = Path[i].V;
		Tree[p+i].LN = L[i];
	}
 
	for (i = 1; i <= 2*n-3; i++)
	{
		if (fabs(Tree[i].LN-epsilona) <= 2*epsilona) { Tree[i].LN = 0.0; }

		ARETE[2*i-2] = Tree[i].U;
		ARETE[2*i-1] = Tree[i].V;
		LONGUEUR[i-1] = Tree[i].LN;
		
		if (LONGUEUR[i-1]<2*epsilona) { LONGUEUR[i-1] = 2*epsilona; }
	} 
 
	free(X);
	free(Tree);
	free(L);
	free(Path);

	for (i=0;i<=n;i++)    
	{
		free(D[i]);
	}
    
	free(D);
}


//=============================================
//
//=============================================
// CONSIGNE_1 Loop plusieurs fois et supprime création outputfile, remplace par structure -> str
void createTree2 (int nbSpecies,double l, const char * outputfilename, string& newickRef){	 
    int i,j;
    double *LONGUEUR,Sigma;        //= Longueurs des aretes
    long int *ARETE;               //= Liste des aretes de l'arbre
    double **DA,**DI;		       //= Distance d'arbre
    
    int n = nbSpecies;
    
	// Allocation dynamique pour les matrices
	DI = (double **) malloc((n+1)*sizeof(double*));
	DA = (double **) malloc((2*n-1)*sizeof(double*));

	for (i = 0; i <= n; i++){
		DI[i] = (double*)malloc((n+1)*sizeof(double));
		 
		if (DI[i]==NULL)
		{
			printf(" Data matrix is too large(2)\n"); 
			exit(1);
		}
	}

    for(i = 0; i <= 2*n-2; i++){
        DA[i] = (double*)malloc((2*n-1)*sizeof(double));
    }	 

    ARETE = (long int *) malloc((4*n)*sizeof(long int));
    LONGUEUR = (double *) malloc((2*n)*sizeof(double));

	//===== CA COMMENCE ICI =====
	
	
	//==================================
	//= CREATION D'UN ARBRE ALEATOIRE
	//==================================
	tree_generation(DA, DI, n, Sigma = 0.0);
	Tree_edges(DA, ARETE, LONGUEUR, n);
	
	float MaxLong = 0;
	int MinLong = INFINI;
    
	for(int j = 1; j <= 2*n-3; j++){
		MaxLong += LONGUEUR[j-1];
	}
    
    double facteur = 1.0/(MaxLong/(2*n-3));
	//printf("\n check parameters %f %f", facteur, MaxLong);
    MaxLong = 0;
    
    for(int j = 1; j <= 2*n-3; j++){
        LONGUEUR[j-1] = (LONGUEUR[j-1]*facteur)*l;
		MaxLong += LONGUEUR[j-1];
	} 

	/*printf("\n check parameters 2 %f %f", facteur, MaxLong);
	for(i=0;i<=2*n;i++){
         printf("\n Longueur %f et Aretes %ld & %ld", LONGUEUR[i], ARETE[2*i], ARETE[2*i+1]);

		 //for (j=1; j<= n; j++) { printf("test %f \t", DI[i][j]); }
    }
	printf("\nl= %lf,facteur = %lf,longueur moyenne des branches = %f\n",l,facteur,MaxLong/(2*n-3));*/
    //=============================================================================
	//= SAUVEGARDE DE L'ARBRE
	//=============================================================================
    SAVEASNewick2(LONGUEUR,ARETE,n, outputfilename, newickRef);
	
	for (i = 0; i <= n; i++){
		free(DI[i]);
	}
	free(DI);
    for(i = 0; i <= 2*n-2; i++){
        free(DA[i]);
    }	
	free(DA);
	free(ARETE);
	free(LONGUEUR);

	
}

//=============================================
//	Function creating a variant tree and computing the RF distance
//=============================================
double swapLeafComputeRF(string treeRef, string& secTree, int sp1, int sp2)
{
	string trueSpecie1, fakeSpecie1, trueSpecie2, fakeSpecie2;
	string espece11, espece12, espece21, espece22;
	string filling = "XXXXX", espece1 = intToString(sp1), espece2 = intToString(sp2);
	int start_pos1, start_pos2;
	double *mat_distances = new double[4];

	
	espece11 = "(" + espece1 + ":";
	espece12 = "," + espece1 + ":";
	start_pos1 = secTree.find(espece11);

	if(start_pos1 <= secTree.length() && start_pos1 >= 0)
	{
		trueSpecie1 = espece11;
		trueSpecie2 = espece11[0] + espece2 + ":";
		fakeSpecie1 = espece12;
		fakeSpecie2 = espece12[0] + espece2 + ":";
	}
	else
	{
		trueSpecie1 = espece12;
		trueSpecie2 = espece12[0] + espece2 + ":";
		fakeSpecie1 = espece11;
		fakeSpecie2 = espece11[0] + espece2 + ":";

	}
	secTree = myreplace(secTree, trueSpecie1, filling);
	start_pos2 = secTree.find(trueSpecie2);

	if(start_pos2 <= secTree.length() && start_pos2 >= 0)
	{ secTree = myreplace(secTree, trueSpecie2, trueSpecie1); }
	else
	{ secTree = myreplace(secTree, fakeSpecie2, fakeSpecie1); }

	secTree = myreplace(secTree, filling, trueSpecie2);
	main_hgt(treeRef, secTree, mat_distances);

	return mat_distances[0];

}


string myreplace(string &s, string toReplace,string replaceWith){
    return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
}

void createSubstituteTrees2(string treeRef, int nSpecies, int nbArbres)
{
	int nbSpecies = nSpecies;
	//double distMax = 2*nbSpecies-6;
	double dix = 0.10, vingtCinq = 0.25, cinquante = 0.50, soixanteQuinze = 0.75;
	int quota0 = 0, quota10 = 0, quota25 = 0, quota50 = 0, quota75 = 0, quotaMax = 0;
	int nbEspPart = 9;
	double pLeavesAbsent = 0.25;
	double *mat_distances = new double[4];
	double distancesRF;
	double nbEspSupp = (double(nbSpecies)*double(pLeavesAbsent))/100.0;
	nbEspSupp = ceil(nbEspSupp);

	vector <string> mesTrees0;
	vector <string> mesTrees10;
	vector <string> mesTrees25;
	vector <string> mesTrees50;
	vector <string> mesTrees75;

	vector <string> mesTrees0_tmp;			
	vector <string> mesTrees10_tmp;
	vector <string> mesTrees25_tmp;
	vector <string> mesTrees50_tmp;
	vector <string> mesTrees75_tmp;

	string espece1, espece11, espece12;
	string espece2;
	
	string specie1 = "";
	string specie2 = "";
			
	int nbSpecie1 = 0;
	int nbSpecie2 = 0;
	
	string trueSpecie1, fakeSpecie1, trueSpecie2, fakeSpecie2, tree2; 
			
	string newTree, nvlArbre;
	int start_pos1;
	int start_pos2;

	for (int i = 1; i < nbSpecies; i++)
	{
		//distancesRF = 0;
		printf("\n");
								
									//Pour savoir si les quotas ne sont pas complétés
		if(quota0 <= nbEspPart || quota10 <= nbEspPart || quota25 <= nbEspPart || quota50 <= nbEspPart || quota75 <= nbEspPart)
		{
			//Les deux possibilitées concatener l'espece1
			espece1 = intToString(i);
			espece11 = "(" + espece1 + ":";
			espece12 = "," + espece1 + ":";
			start_pos1 = treeRef.find(espece11);
			if(start_pos1 <= treeRef.length() && start_pos1 >=0){
				trueSpecie1 = espece11;
				fakeSpecie1 = espece12;
			}
			else{
				trueSpecie1 = espece12;
				fakeSpecie1 = espece11;
			}
															
			for (int j = (i+1); j <= nbSpecies; j++)
			{
				distancesRF = 0;
				//printf("\n normalement on rentre ici");
				if(quota0 <= nbEspPart || quota10 <= nbEspPart || quota25 <= nbEspPart || quota50 <= nbEspPart || quota75 <= nbEspPart)
				{
					newTree = treeRef;
					string filling = "XXXXXX";
					newTree = myreplace(newTree, trueSpecie1 ,filling);
					espece2 = intToString(j);				
					trueSpecie2 = trueSpecie1[0] + espece2 + ":";
					fakeSpecie2 = fakeSpecie1[0] + espece2 + ":";
													
					start_pos2 = newTree.find(trueSpecie2);
					int treeSize = newTree.length();

					if(start_pos2 <= treeSize && start_pos2 >=0){
						newTree = myreplace(newTree, trueSpecie2 ,trueSpecie1);
					}
					else if (start_pos2 < 0)
					{
						newTree = myreplace(newTree, fakeSpecie2 ,fakeSpecie1);
					}
					//nvlArbre = treeRef;
					distancesRF = swapLeafComputeRF(treeRef, newTree, i, j);
					//newTree = nvlArbre;
					newTree = myreplace(newTree, filling ,trueSpecie2);
					main_hgt(treeRef, newTree, mat_distances);
					distancesRF = mat_distances[0];
					if(distancesRF==0){
						distancesRF = 0.0;
					}

					if(quota0<=nbEspPart){
						mesTrees0.push_back(treeRef);
						mesTrees0_tmp.push_back(treeRef);
						quota0++;
					}
													
					if(distancesRF<soixanteQuinze)
					{
						tree2 = newTree;
						for(int nbModif = 0; nbModif<nbEspSupp; nbModif++){
							tree2 = modifTree (nbSpecies,tree2,i,j);
						}
						if(distancesRF>0 && distancesRF<=dix && quota10<=nbEspPart)
						{
							mesTrees10.push_back(tree2);
							mesTrees10_tmp.push_back(tree2);
							quota10++;
						}								
					
						else if(distancesRF>dix && distancesRF<=vingtCinq && quota25<=nbEspPart)
						{		
							mesTrees25.push_back(tree2);
							mesTrees25_tmp.push_back(tree2);
							quota25++;
						}
						else if(distancesRF>vingtCinq && distancesRF<=cinquante && quota50<=nbEspPart)
						{										
							mesTrees50.push_back(tree2);
							mesTrees50_tmp.push_back(tree2);
							quota50++;
						}
						else if(distancesRF>cinquante && distancesRF<=soixanteQuinze && quota75<=nbEspPart)
						{										
							mesTrees75.push_back(tree2);
							mesTrees75_tmp.push_back(tree2);
							quota75++;
						}
					}
					else if(distancesRF>soixanteQuinze){
						//tree75Noise = newTree;
						quotaMax++;
					}	
				}
			}
		}					
	}
	if(std::find(mesTrees10.begin(), mesTrees10.end(), treeRef) != mesTrees10.end()) {printf("\nhaha");}
	
	//printf("\n okay donc les résultats sont : Quota0=%d; Quota10=%d; Quota25=%d; Quota50=%d; Quota75=%d; QuotaMax=%d", quota0, quota10, quota25, quota50, quota75, quotaMax);
}

void createClusters(string treeRef, int nbSpecies, int nbArbres, double lowLimit, double highLimit, vector <string>& allTheTrees)
{
	int nbEspPart = 9;
	double pLeavesAbsent = 0.25;
	double *matrices = new double[4];
	double nbEspSupp = (double(nbSpecies)*double(pLeavesAbsent))/100.0;
	nbEspSupp = ceil(nbEspSupp);

	string newTree, nvlArbre;
	vector <string> mesArbres;
	vector <string> arbresAux;
	vector <string> arbresTemp;
	int quota = 0, quotaT = 1, x, y, i, j, nb = 0, first = 0, compt;
	double distancesRF;

	mesArbres.push_back(treeRef);
	arbresAux.push_back(treeRef);
	allTheTrees.push_back(treeRef);
	arbresTemp.push_back(treeRef);
	x = 0; y = 1;
	int stop = 5, count = 0;
	printf("\n\n NEW CLUSTER ! \n");

	while(quota < nbArbres)
	{
		quotaT = 1;
		if(quota == x)
		{
			printf("\n check x %d\n", x);
			newTree = arbresAux[x];
			printf("\n arbre de X %s", newTree.c_str());
			x++;
		}
		else if(y <= quota)
		{
			printf("\n check y \n");
			newTree = mesArbres[y];
			y++;
		}
		else {
			printf("\nProblème quelque part :( \n");
			exit(1);
		}
		printf("\n yellow");
		for(i = 1; i < nbSpecies; i++)
		{
			for(j = i + 1; j <= nbSpecies; j++)
			{
				nvlArbre = newTree;
				distancesRF = swapLeafComputeRF(treeRef, nvlArbre, i, j);
				
				if(distancesRF > lowLimit && distancesRF <= highLimit && std::find(arbresTemp.begin(), arbresTemp.end(), nvlArbre) == arbresTemp.end())
				{
					//mesArbres.push_back(nvlArbre);
					//allTheTrees.push_back(nvlArbre);
					arbresTemp.push_back(nvlArbre);
					quotaT++;
					if(quotaT == nbArbres) { i = j = nbSpecies; }
				}
				else
				{
					arbresAux.push_back(nvlArbre);
				}
			}
		}

		compt = 1;
		/*while(quotaT > 0)
		{
			printf("\n%s\n%s\n\n", mesArbres[compt-1].c_str(), arbresTemp[quotaT-1].c_str());
			//printf("\n On rentre dans la boucle ?");
			main_hgt(mesArbres[compt-1].c_str(), arbresTemp[quotaT-1].c_str(), matrices);
			if(matrices[0] == 0 || matrices[0] > highLimit || matrices[0] < lowLimit){
				arbresAux.push_back(arbresTemp[quotaT-1]);
				compt = 1;
				quotaT = quotaT-1;
				printf("BAD ? %f  %d\t\n", matrices[0], quotaT);
			}
			else if(matrices[0]<=highLimit){
				compt++;
				printf("GOOD ? %f %d \t %s \n", matrices[0], compt, "hello");
				if(compt >= quota){
					mesArbres.push_back(arbresTemp[quotaT-1].c_str());
					allTheTrees.push_back(arbresTemp[quotaT-1].c_str());

					//printf("\n%s\n", arbresTemp[quotaT-1].c_str());
					quota++;
					compt = 1;
					printf("QUOTA=%d \n\n", quotaT);
					quotaT = quotaT-1;
					if(quota == nbArbres){quotaT = 0;}
				}
			}
		}*/
		printf("\n on fonctionne ici ?");
		for(i = 0; i < quotaT; i++)
		{
			//printf("\n Start first for loop\n");
			int correct = 1;
			for(j = 0; j <= quota; j++)
			{
				printf("\t Start second for loop \t");
				printf("\n cello");
				main_hgt(mesArbres[j].c_str(), arbresTemp[i].c_str(), matrices);
				if(matrices[0] > lowLimit && matrices[0] <= highLimit)
				{ printf("\n At least one good tree"); }
				else{
					correct = 0;
					j = quota +1;
				}
			}

			if(correct == 1)
			{
				mesArbres.push_back(arbresTemp[i].c_str());
				allTheTrees.push_back(arbresTemp[i].c_str());
				quota++;
				if(quota == nbArbres) { i = quotaT;}
			}
			else { arbresAux.push_back(arbresTemp[i].c_str()); 
				printf("\n on ajoute un mauvais arbre à arbresAux");}

			printf("\n Ok le nouveau quota de bons arbres est : %d\n", quota);
		}
		printf("\n QuotaT vaut %d donc on est sorti :) \n", quotaT);
		arbresTemp.clear();
	}
	printf("\n\t\t quota final = %f   %d\n", matrices[0], quota);
	printf("\n====================================================\n====================================================");

}

int main(int nargs,char ** argv)
{
	//if(nargs == 1){ printf("nope !"); exit(-1);}
	//int nb_Clusters = atoi(argv[1]), nb_Feuilles = atoi(argv[2]), noiseLvl = atoi(argv[3]), nbTrees = atoi(argv[4]);
	//int nbTreesTot = (nbTrees+1) * nb_Clusters;
	int nbarbres[10];
	int nbFeuilles[10];
	int nbClusters[10];
	int nb_trees[10];
	nbarbres[0] = 19, nbarbres[1] = 9, nbarbres[2] = 6, nbarbres[3] = 4, nbarbres[4] = 3, nbarbres[5] = 1;
	nbFeuilles[0] = 8, nbFeuilles[1] = 16, nbFeuilles[2] = 32, nbFeuilles[3] = 64;
	nbClusters[0] = 1, nbClusters[1] = 2, nbClusters[2] = 3, nbClusters[3] = 4, nbClusters[4] = 5, nbClusters[5] = 10;
	nb_trees[0] = 19, nb_trees[1] = 9, nb_trees[2] = 6, nb_trees[3] = 4, nb_trees[4] = 3, nb_trees[5] = 1; 
	float volNoise[10];
	volNoise[0] = 0.0, volNoise[1] = 0.1, volNoise[2] = 0.25, volNoise[3] = 0.5, volNoise[4] = 0.75;
	int startL = 0, endL, startN = 1, endN, startC, endC, limSup;
	string refTree;
	double *mat_dist = new double[4];
	double lowNoise, highNoise, RF;
	vector <string> allTrees;
	FILE * outfile = fopen("test_auto2.txt","w");
	//printf("On arrive ici au moins ?");
	char okay = argv[1][1];
	
	//Pour la génération de plusieurs matrices
	int nbTreesTot = 15;
	if(okay == 'm')
	{
		// enchainement de cout/cin pour récupérer les infos
		int nb_Clusters, leaves, noiseLvl;
		int nb_repet = atoi(argv[2]);
		cout<<"How many clusters ?";
		cin>>nb_Clusters;
		cout<<"How many leaves ?";
		cin>>leaves;
		cout<<"Level of noise ?  (1 = [0-10], 2 = [10-25], 3 = [25-50], 4 = [50-75]) \nYou cannot have [0-10] noise with 8 leaves.\n";
		cin>>noiseLvl;
		lowNoise = volNoise[noiseLvl-1];
		highNoise = volNoise[noiseLvl];
		for(int m = 1; m <= nb_repet; m++)
		{
			//cout<<"\nHello :)";
			//printf("\t what ? %d", nb_Clusters);
			for(int i = 1; i <= nb_Clusters; i++)
			{
				createTree2(leaves, 1, "something", refTree);
				printf("\n %s \n", refTree.c_str());
				createClusters(refTree, leaves, 2, lowNoise, highNoise, allTrees);
			}
			printf("\n\n");
			fprintf(outfile, "%d   %d   %d   0   %d", nbTreesTot, leaves, nb_Clusters, noiseLvl);
			for( int i = 1; i <= nbTreesTot; i++)
			{
				fprintf(outfile, "\n");
				printf("\n%s", allTrees[i-1].c_str());
				for(int j = 0; j < nbTreesTot; j++)
				{
					main_hgt(allTrees[i-1].c_str(), allTrees[j].c_str(), mat_dist);
					RF = mat_dist[0];
					fprintf(outfile, "%f   ", RF);
				}
			}
			allTrees.clear();
			fprintf(outfile, "\n");
		}
	}

	else{
		if(okay == 'L')
		{
			nbFeuilles[0] = atoi(argv[2]);
			endL = 1;
			endN = 3;
			endC = 6;
		}
		else if(argv[1] == "C")
		{
			nbClusters[0] = atoi(argv[2]);
			endC = 1;
			endL = 4;
			endN = 5;
		}
		else if(argv[1] == "-N")
		{
			limSup = atoi(argv[2]);
			volNoise[0] = volNoise[limSup-1];
			volNoise[1] = volNoise[limSup];

			if(limSup == 1){
				startL = 1; 
			}
			endN = 1;
			endL = 4;
			endC = 6;
		}
		//printf("\n on va checker les arbres de %d à  %d \n", startL, endL);

		for(int a = startL; a < endL; a++)
		{
			int nb_Feuilles = nbFeuilles[a];
			if (nb_Feuilles == 8){
				startN = 2;
				startC = 3;
			}
			else if (nb_Feuilles == 16){
				startC = 2;
				startN = 1;
			}
			else{
				startC = 0;
				startN = 1;
			}
			for(int b = startN; b < endN; b++)
			{
				float lowNoise = volNoise[b-1];
				float highNoise = volNoise[b];
				for(int c = 2; c < 5; c++)
				{
					int nb_Clusters = nbClusters[c];
					int nbTrees = nb_trees[c];
					for(int d = 0; d < nb_Clusters; d++)
					{
						createTree2(nb_Feuilles, 0.5, "something", refTree);
						createClusters(refTree, nb_Feuilles, nbTrees, lowNoise, highNoise, allTrees);
					}
					fprintf(outfile, "%d   %d   %d   0   %f", nbTreesTot, nb_Feuilles, nb_Clusters, highNoise*100);
					for( int i = 1; i <= nbTreesTot; i++)
					{
						fprintf(outfile, "\n");
						//printf("\n\n%s", allTrees[i-1].c_str());
						for(int j = 1; j <= nbTreesTot; j++)
						{
							main_hgt(allTrees[i-1].c_str(), allTrees[j-1].c_str(), mat_dist);
							RF = mat_dist[0];
							fprintf(outfile, "%f   ", RF);
						}
					}
					allTrees.clear();
					fprintf(outfile, "\n");
					//fprintf(outfile, "\n %s", refTree.c_str());
				}
			}
		}
	}
	


	/*if(noiseLvl == 10)
	{
		lowNoise = 0;
		highNoise = 0.10;
	}
	else if(noiseLvl == 25)
	{
		lowNoise = 0.10;
		highNoise = 0.25;
	}
	else if(noiseLvl == 50)
	{
		lowNoise = 0.25;
		highNoise = 0.50;
	}
	else if(noiseLvl == 75)
	{
		lowNoise = 0.50;
		highNoise = 0.75;
	}

	fprintf(outfile, "%d   %d   %d   0   %d", nbTreesTot, nbSpecies, nb_Clusters, noiseLvl);

	for(int i = 1; i <= nb_Clusters; i++) {
		createTree2(nbSpecies, 0.5, "something", refTree);
		createClusters(refTree, nbSpecies, nbTrees, lowNoise, highNoise, allTrees);
	}

	for( int i = 1; i <= nbTreesTot; i++)
	{
		fprintf(outfile, "\n");
		//printf("\n\n%s", allTrees[i-1].c_str());
		for(int j = 1; j <= nbTreesTot; j++)
		{
			main_hgt(allTrees[i-1].c_str(), allTrees[j-1].c_str(), mat_dist);
			RF = mat_dist[0];
			fprintf(outfile, "%f   ", RF);
		}
	}*/
	fprintf(outfile, "\n");
	fclose(outfile);
	printf("\n");
}