#include <R.h>
#include "Rmath.h"
#include <iostream>
using namespace std;

//#include "others/combn.cpp"
//#include "others/Cnp.cpp"

extern"C" {
  
  void Dcov1C(double *X, int *vecd, double *a, int *n, int *p, int *weightchoice, double *res) {

    double gammajl(int j, int l, double *X, double *a, int *n, int *vecd, int *weightchoice);
    double gammajjprimel(int j, int jprime, int l, double *X, double *a, int *n, int *vecd, int *weightchoice);
    double betajjprimel(int j, int jprime, int l, double *X, double *a, int *n, int *vecd, int *weightchoice);
    void combn(int *combmat, int *n, int *m);
    double Cnp(int n, int k);

    double xijjprimel(int j, int jprime, int l, double *X, double *a, int *n, int *vecd, int *weightchoice);
    double term1=0, term2=0, term3=1, somme;
    int j, jprime, l, lprime, n2=n[0]*n[0], *indB1, *indBprime1, indB2, indBprime2, *combmatB, *combmatBprime, *B, *Bprime, tB, tBprime, cmpt;
    double CpindB, CpindBprime, prod1, prod2, prod3;

    indB1 = new int[1];
    indBprime1 = new int[1];

    for (indB1[0] = 1; indB1[0] <= p[0]; indB1[0]++) {
      CpindB = Cnp(p[0],indB1[0]);
      combmatB = new int[(int)CpindB*indB1[0]];
      combn(combmatB,p,indB1);
      for (indB2 = 1; indB2 <= (int)CpindB; indB2++) {
	B = new int[indB1[0]];
	for (tB = 1; tB <= indB1[0]; tB++) B[tB-1] = combmatB[(indB2-1)*indB1[0]+(tB-1)];
	for (indBprime1[0] = 1; indBprime1[0] <= p[0]; indBprime1[0]++) {
	  CpindBprime = Cnp(p[0],indBprime1[0]);
	  combmatBprime = new int[(int)CpindBprime*indBprime1[0]];
	  combn(combmatBprime,p,indBprime1);
	  for (indBprime2 = 1; indBprime2 <= (int)CpindBprime; indBprime2++) {
	    Bprime = new int[indBprime1[0]];
	    for (tBprime = 1; tBprime <= indBprime1[0]; tBprime++) Bprime[tBprime-1] = combmatBprime[(indBprime2-1)*indBprime1[0]+(tBprime-1)];
	    term1 = 0.0;
	    for (j = 1; j <= n[0]; j++) {
	      for (jprime = 1; jprime <= n[0]; jprime++) {
		prod1 = 1.0;
		for (l = 1; l <= indB1[0]; l++) {
		  for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
		    if (B[l-1] == Bprime[lprime-1]) prod1 = prod1 * betajjprimel(j,jprime,B[l-1],X,a,n,vecd,weightchoice);
		  }
		}
		prod2 = 1.0;
		for (l = 1; l <= indB1[0]; l++) {
		  cmpt = 0;
		  for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
		    if (B[l-1] == Bprime[lprime-1]) cmpt = cmpt + 1;
		  }
		  if (cmpt == 0) prod2 = prod2 * gammajl(j,B[l-1],X,a,n,vecd,weightchoice);
		}
		prod3 = 1.0;
		for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
		  cmpt = 0;
		  for (l = 1; l <= indB1[0]; l++) {
		    if (Bprime[lprime-1] == B[l-1]) cmpt = cmpt + 1;
		  }
		  if (cmpt == 0) prod3 = prod3 * gammajl(jprime,Bprime[lprime-1],X,a,n,vecd,weightchoice);
		}
		term1 = term1 + prod1 * prod2 * prod3;
	      }
	    }
	    term1 = term1 / (double)n2;
	    term2 = 0.0;
	    for (j = 1; j <= n[0]; j++) {
	      prod1 = 1.0;
		for (l = 1; l <= indB1[0]; l++) {
		  for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
		    if (B[l-1] == Bprime[lprime-1]) {
			somme = 0.0;
			for (jprime = 1; jprime <= n[0]; jprime++) somme = somme + betajjprimel(j,jprime,B[l-1],X,a,n,vecd,weightchoice);
			somme = somme / (double)n[0];
			prod1 = prod1 * somme;
		      }
		  }
		}
		prod2 = 1.0;
		for (l = 1; l <= indB1[0]; l++) {
		  cmpt = 0;
		  for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
		    if (B[l-1] == Bprime[lprime-1]) cmpt = cmpt + 1;
		  }
		  if (cmpt == 0) prod2 = prod2 * gammajl(j,B[l-1],X,a,n,vecd,weightchoice);
		}
		prod3 = 1.0;
		for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
		  cmpt = 0;
		  for (l = 1; l <= indB1[0]; l++) {
		    if (Bprime[lprime-1] == B[l-1]) cmpt = cmpt + 1;
		  }
		  if (cmpt == 0) {
		    somme = 0.0;
		    for (jprime = 1; jprime <= n[0]; jprime++) somme = somme + gammajl(jprime,Bprime[lprime-1],X,a,n,vecd,weightchoice);
		    somme = somme / (double)n[0];
		    prod3 = prod3 * somme;
		  }
		}
		term2 = term2 + prod1 * prod2 * prod3;
	      }
	    term2 = -2.0 * term2 / (double)n[0];
	    prod1 = 1.0;
	    for (l = 1; l <= indB1[0]; l++) {
	      for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
		if (B[l-1] == Bprime[lprime-1]) {
		  somme = 0.0;
		  for (j = 1; j <= n[0]; j++) {
		    for (jprime = 1; jprime <= n[0]; jprime++) {
		      somme = somme + betajjprimel(j,jprime,B[l-1],X,a,n,vecd,weightchoice);
		    }
		  }
		  somme = somme / (double)n2;
		  prod1 = prod1 * somme;
		}
	      }
	    }
	    prod2 = 1.0;
	    for (l = 1; l <= indB1[0]; l++) {
	      cmpt = 0;
	      for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
		if (B[l-1] == Bprime[lprime-1]) cmpt = cmpt + 1;
	      }
	      if (cmpt == 0) {
		somme = 0.0;
		for (jprime = 1; jprime <= n[0]; jprime++) {
		  somme = somme + gammajl(jprime,B[l-1],X,a,n,vecd,weightchoice);
		}
		somme = somme / (double)n[0];
		prod2 = prod2 * somme;
	      }
	    }
	    prod3 = 1.0;
	    for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
	      cmpt = 0;
	      for (l = 1; l <= indB1[0]; l++) {
		if (Bprime[lprime-1] == B[l-1]) cmpt = cmpt + 1;
	      }
	      if (cmpt == 0) {
		somme = 0.0;
		for (jprime = 1; jprime <= n[0]; jprime++) {
		  somme = somme + gammajl(jprime,Bprime[lprime-1],X,a,n,vecd,weightchoice);
		}
		somme = somme / (double)n[0];
		prod3 = prod3 * somme;
	      }
	    }
	    term3 = prod1 * prod2 * prod3;
	    
	    res[0] = res[0] + R_pow(-1.0,(double)(indB1[0]+indBprime1[0])) * (term1 + term2 + term3);
	      
	    delete[] Bprime;
	  }
	  delete[] combmatBprime;

	}
	delete[] B;
      }
      delete[] combmatB;
    }
    
    delete[] indB1;
    delete[] indBprime1;
    
  }
  void Dcov1Cnormed(double *X, int *vecd, double *a, int *n, int *p, int *weightchoice, double *res, double *denom1) {

    double gammajl(int j, int l, double *X, double *a, int *n, int *vecd, int *weightchoice);
    double gammajjprimel(int j, int jprime, int l, double *X, double *a, int *n, int *vecd, int *weightchoice);
    double betajjprimel(int j, int jprime, int l, double *X, double *a, int *n, int *vecd, int *weightchoice);
    void combn(int *combmat, int *n, int *m);
    double Cnp(int n, int k);
    int Union(int *arr1, int *arr2, int m, int n, int *arr3);
    int Intersection(int *arr1, int *arr2, int m, int n, int *arr3);
    int Setdiff(int *arr1, int *arr2, int m, int n, int p, int *arr3);

    double xijjprimel(int j, int jprime, int l, double *X, double *a, int *n, int *vecd, int *weightchoice);
    double term1=0, term2=0, term3=1, somme;
    int j, jprime, l, lprime, n2=n[0]*n[0], *indB1, *indBprime1, indB2, indBprime2, *combmatB, *combmatBprime, *B, *Bprime, tB, tBprime, cmpt;
    int cardBuBprime, *BuBprime, cardBnBprime, *BnBprime, cardSymdif, *Symdif;
    double CpindB, CpindBprime, prod1, prod2, prod3, *blvec, *glvec, prodtmp;

    denom1[0] = 0.0;

    blvec = new double[p[0]];
    glvec = new double[p[0]];

    indB1 = new int[1];
    indBprime1 = new int[1];


    for (l = 1; l <= p[0]; l++) {
      blvec[l-1] = 0.0;
      glvec[l-1] = 0.0;
      for (j = 1; j <= n[0]; j++) {
	for (jprime = 1; jprime <= n[0]; jprime++) {
	  blvec[l-1] =  blvec[l-1] + betajjprimel(j,jprime,l,X,a,n,vecd,weightchoice);
	}
	glvec[l-1] = glvec[l-1] + gammajl(j,l,X,a,n,vecd,weightchoice);
      }
      blvec [l-1] = blvec[l-1] / n2;
      glvec[l-1] = glvec[l-1] / n[0];
    }
  

    for (indB1[0] = 1; indB1[0] <= p[0]; indB1[0]++) {
      CpindB = Cnp(p[0],indB1[0]);
      combmatB = new int[(int)CpindB*indB1[0]];
      combn(combmatB,p,indB1);
      for (indB2 = 1; indB2 <= (int)CpindB; indB2++) {
	B = new int[indB1[0]];
	for (tB = 1; tB <= indB1[0]; tB++) B[tB-1] = combmatB[(indB2-1)*indB1[0]+(tB-1)];
	for (indBprime1[0] = 1; indBprime1[0] <= p[0]; indBprime1[0]++) {
	  CpindBprime = Cnp(p[0],indBprime1[0]);
	  combmatBprime = new int[(int)CpindBprime*indBprime1[0]];
	  combn(combmatBprime,p,indBprime1);
	  for (indBprime2 = 1; indBprime2 <= (int)CpindBprime; indBprime2++) {
	    Bprime = new int[indBprime1[0]];
	    for (tBprime = 1; tBprime <= indBprime1[0]; tBprime++) Bprime[tBprime-1] = combmatBprime[(indBprime2-1)*indBprime1[0]+(tBprime-1)];
	    BuBprime = new int[indB1[0] + indBprime1[0]];
	    cardBuBprime = Union(B,Bprime,indB1[0],indBprime1[0], BuBprime);
	    if (cardBuBprime != 1) {
	      BnBprime = new int[indB1[0]];
	      cardBnBprime = Intersection(B,Bprime,indB1[0],indBprime1[0],BnBprime);
	      Symdif = new int[cardBuBprime];
	      cardSymdif = Setdiff(BuBprime,BnBprime,cardBuBprime,cardBnBprime,p[0],Symdif);
	      prodtmp = 1.0;
	      for (l = 0; l <cardBnBprime; l++) prodtmp = prodtmp * blvec[BnBprime[l]-1];
	      for (l = 0; l <cardSymdif; l++) prodtmp = prodtmp * (-1.0) * glvec[Symdif[l]-1]; 
	      denom1[0] = denom1[0] + (cardBuBprime - 1) * prodtmp;
	      delete[] BnBprime;
	      delete[] Symdif;
	    }
	    delete[] BuBprime;
	    term1 = 0.0;
	    for (j = 1; j <= n[0]; j++) {
	      for (jprime = 1; jprime <= n[0]; jprime++) {
		prod1 = 1.0;
		for (l = 1; l <= indB1[0]; l++) {
		  for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
		    if (B[l-1] == Bprime[lprime-1]) prod1 = prod1 * betajjprimel(j,jprime,B[l-1],X,a,n,vecd,weightchoice);
		  }
		}
		prod2 = 1.0;
		for (l = 1; l <= indB1[0]; l++) {
		  cmpt = 0;
		  for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
		    if (B[l-1] == Bprime[lprime-1]) cmpt = cmpt + 1;
		  }
		  if (cmpt == 0) prod2 = prod2 * gammajl(j,B[l-1],X,a,n,vecd,weightchoice);
		}
		prod3 = 1.0;
		for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
		  cmpt = 0;
		  for (l = 1; l <= indB1[0]; l++) {
		    if (Bprime[lprime-1] == B[l-1]) cmpt = cmpt + 1;
		  }
		  if (cmpt == 0) prod3 = prod3 * gammajl(jprime,Bprime[lprime-1],X,a,n,vecd,weightchoice);
		}
		term1 = term1 + prod1 * prod2 * prod3;
	      }
	    }
	    term1 = term1 / (double)n2;
	    term2 = 0.0;
	    for (j = 1; j <= n[0]; j++) {
	      prod1 = 1.0;
		for (l = 1; l <= indB1[0]; l++) {
		  for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
		    if (B[l-1] == Bprime[lprime-1]) {
			somme = 0.0;
			for (jprime = 1; jprime <= n[0]; jprime++) somme = somme + betajjprimel(j,jprime,B[l-1],X,a,n,vecd,weightchoice);
			somme = somme / (double)n[0];
			prod1 = prod1 * somme;
		      }
		  }
		}
		prod2 = 1.0;
		for (l = 1; l <= indB1[0]; l++) {
		  cmpt = 0;
		  for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
		    if (B[l-1] == Bprime[lprime-1]) cmpt = cmpt + 1;
		  }
		  if (cmpt == 0) prod2 = prod2 * gammajl(j,B[l-1],X,a,n,vecd,weightchoice);
		}
		prod3 = 1.0;
		for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
		  cmpt = 0;
		  for (l = 1; l <= indB1[0]; l++) {
		    if (Bprime[lprime-1] == B[l-1]) cmpt = cmpt + 1;
		  }
		  if (cmpt == 0) {
		    somme = 0.0;
		    for (jprime = 1; jprime <= n[0]; jprime++) somme = somme + gammajl(jprime,Bprime[lprime-1],X,a,n,vecd,weightchoice);
		    somme = somme / (double)n[0];
		    prod3 = prod3 * somme;
		  }
		}
		term2 = term2 + prod1 * prod2 * prod3;
	      }
	    term2 = -2.0 * term2 / (double)n[0];
	    prod1 = 1.0;
	    for (l = 1; l <= indB1[0]; l++) {
	      for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
		if (B[l-1] == Bprime[lprime-1]) {
		  somme = 0.0;
		  for (j = 1; j <= n[0]; j++) {
		    for (jprime = 1; jprime <= n[0]; jprime++) {
		      somme = somme + betajjprimel(j,jprime,B[l-1],X,a,n,vecd,weightchoice);
		    }
		  }
		  somme = somme / (double)n2;
		  prod1 = prod1 * somme;
		}
	      }
	    }
	    prod2 = 1.0;
	    for (l = 1; l <= indB1[0]; l++) {
	      cmpt = 0;
	      for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
		if (B[l-1] == Bprime[lprime-1]) cmpt = cmpt + 1;
	      }
	      if (cmpt == 0) {
		somme = 0.0;
		for (jprime = 1; jprime <= n[0]; jprime++) {
		  somme = somme + gammajl(jprime,B[l-1],X,a,n,vecd,weightchoice);
		}
		somme = somme / (double)n[0];
		prod2 = prod2 * somme;
	      }
	    }
	    prod3 = 1.0;
	    for (lprime = 1; lprime <= indBprime1[0]; lprime++) {
	      cmpt = 0;
	      for (l = 1; l <= indB1[0]; l++) {
		if (Bprime[lprime-1] == B[l-1]) cmpt = cmpt + 1;
	      }
	      if (cmpt == 0) {
		somme = 0.0;
		for (jprime = 1; jprime <= n[0]; jprime++) {
		  somme = somme + gammajl(jprime,Bprime[lprime-1],X,a,n,vecd,weightchoice);
		}
		somme = somme / (double)n[0];
		prod3 = prod3 * somme;
	      }
	    }
	    term3 = prod1 * prod2 * prod3;
	    
	    res[0] = res[0] + R_pow(-1.0,(double)(indB1[0]+indBprime1[0])) * (term1 + term2 + term3);
	      
	    delete[] Bprime;
	  }
	  delete[] combmatBprime;

	}
	if (indB1[0] != 1) {
	  prodtmp = 1.0;
	  for (l = 0; l <indB1[0]; l++) prodtmp = prodtmp * (-1.0) * glvec[B[l]-1];
	  denom1[0] = denom1[0] + 2 * (indB1[0] - 1) * prodtmp;
	}
	delete[] B;
      }
      delete[] combmatB;
    }

    delete[] indB1;
    delete[] indBprime1;
    delete[] blvec;
    delete[] glvec;
    
  }
  double gammajl(int j, int l, double *X, double *a, int *n, int *vecd, int *weightchoice) {
    int indbeginblocl = 1;
    int i, k;
    double prod = 1.0, somme = 0.0, tmp = 0.0;
    
    if (l != 1) {
      for (i = 1; i <= (l-1); i++) {
	indbeginblocl = indbeginblocl + vecd[i-1];
      }
    }
    
    if (weightchoice[0] == 1) {
      for (k = 1; k <= vecd[l-1]; k++) {
	prod = prod * exp(-R_pow(a[0] * X[(indbeginblocl+k-2)*n[0] + (j-1)],2.0) / 2.0);
	  }
    }
    if (weightchoice[0] == 2) {
      for (k = 1; k <= vecd[l-1]; k++) {
	prod = prod / (1.0 + R_pow(a[0] * (X[(indbeginblocl+k-2)*n[0] + (j-1)]),2.0));
	  }
    }
    if (weightchoice[0] == 3) {
      for (k = 1; k <= vecd[l-1]; k++) {
	somme = somme + R_pow(X[(indbeginblocl+k-2)*n[0] + (j-1)],2.0);
      }
      prod = R_pow(somme,a[0]/2.0);
    }
    if (weightchoice[0] == 4) {
      for (k = 1; k <= vecd[l-1]; k++) {
	prod = prod * (-2.0*fabs(X[(indbeginblocl+k-2)*n[0] + (j-1)]) + fabs(X[(indbeginblocl+k-2)*n[0] + (j-1)] - 2.0*a[0]) + fabs(X[(indbeginblocl+k-2)*n[0] + (j-1)] + 2.0*a[0])) / (4.0*fabs(a[0]));
	  }
    }
    if (weightchoice[0] == 5) {
      for (k = 1; k <= vecd[l-1]; k++) {
	tmp = X[(indbeginblocl+k-2)*n[0] + (j-1)];
	prod = prod * (1.0 - R_pow(a[0]*tmp,2.0)) * exp(- R_pow(a[0]*tmp,2.0) / 2.0);
	  }
    }
    return 1.0 - prod;
  }

  double gammajjprimel(int j, int jprime, int l, double *X, double *a, int *n, int *vecd, int *weightchoice) {
    int indbeginblocl = 1;
    int i, k;
    double prod = 1.0, somme = 0.0, tmp = 0.0;
    
    if (l != 1) {
      for (i = 1; i <= (l-1); i++) {
	indbeginblocl = indbeginblocl + vecd[i-1];
      }
    }
    
    if (weightchoice[0] == 1) {
      for (k = 1; k <= vecd[l-1]; k++) {
	prod = prod * exp(-R_pow(a[0] * (X[(indbeginblocl+k-2)*n[0] + (j-1)] - X[(indbeginblocl+k-2)*n[0] + (jprime-1)]),2.0) / 2.0);
	  }
    }
    if (weightchoice[0] == 2) {
      for (k = 1; k <= vecd[l-1]; k++) {
	prod = prod / (1.0 + R_pow(a[0] * (X[(indbeginblocl+k-2)*n[0] + (j-1)] - X[(indbeginblocl+k-2)*n[0] + (jprime-1)]),2.0));
	  }
    }
    if (weightchoice[0] == 3) {
      for (k = 1; k <= vecd[l-1]; k++) {
	somme = somme + R_pow(X[(indbeginblocl+k-2)*n[0] + (j-1)] - X[(indbeginblocl+k-2)*n[0] + (jprime-1)],2.0);
      }
      prod = R_pow(somme,a[0]/2.0);
    }
    if (weightchoice[0] == 4) {
      for (k = 1; k <= vecd[l-1]; k++) {
	prod = prod * (-2.0*fabs(X[(indbeginblocl+k-2)*n[0] + (j-1)] - X[(indbeginblocl+k-2)*n[0] + (jprime-1)]) + fabs(X[(indbeginblocl+k-2)*n[0] + (j-1)] - X[(indbeginblocl+k-2)*n[0] + (jprime-1)] - 2.0*a[0]) + fabs(X[(indbeginblocl+k-2)*n[0] + (j-1)] - X[(indbeginblocl+k-2)*n[0] + (jprime-1)] + 2.0*a[0])) / (4.0*fabs(a[0]));
	  }
    }
    if (weightchoice[0] == 5) {
      for (k = 1; k <= vecd[l-1]; k++) {
	tmp = X[(indbeginblocl+k-2)*n[0] + (j-1)] - X[(indbeginblocl+k-2)*n[0] + (jprime-1)];
	prod = prod * (1.0 - R_pow(a[0]*tmp,2.0)) * exp(- R_pow(a[0]*tmp,2.0) / 2.0);
	  }
    }
    return 1.0 - prod;

  }

  double betajjprimel(int j, int jprime, int l, double *X, double *a, int *n, int *vecd, int *weightchoice) {
    return gammajl(j,l,X,a,n,vecd,weightchoice) + gammajl(jprime,l,X,a,n,vecd,weightchoice) - gammajjprimel(j,jprime,l,X,a,n,vecd,weightchoice);
  }


  void Dcov2C(double *X, int *vecd, double *a, int *n, int *p, int *weightchoice, double *res) {
    
    double xijjprimel(int j, int jprime, int l, double *X, double *a, int *n, int *vecd, int *weightchoice);
    double term1=0, term2=0, term3=1, somme, prod;
    int j, jprime, l, n2=n[0]*n[0];


      for (j = 1; j <= n[0]; j++) {
	for (jprime = 1; jprime <= n[0]; jprime++) {
	  prod = 1.0;
	    for (l = 1; l <= p[0]; l++) {
	      prod = prod * xijjprimel(j,jprime,l,X,a,n,vecd,weightchoice);
	    }
	  term1 = term1 + prod;
	}
      }
      term1 = term1 / (double)n2;

      for (j = 1; j <= n[0]; j++) {
	prod = 1.0;
	for (l = 1; l <= p[0]; l++) {
	  somme = 0.0;
	  for (jprime = 1; jprime <= n[0]; jprime++) {
	    somme = somme + xijjprimel(j,jprime,l,X,a,n,vecd,weightchoice);
	  }
	  somme = somme / (double)(n[0]);
	  prod = prod * somme;
	}
	term2 = term2 + prod;
      }
      term2 = -2.0 * term2 / (double)(n[0]);

      for (l = 1; l <= p[0]; l++) {
	somme = 0.0;
	for (j = 1; j <= n[0]; j++) {
	  for (jprime = 1; jprime <= n[0]; jprime++) {
	    somme = somme + xijjprimel(j,jprime,l,X,a,n,vecd,weightchoice);
	  }
	}
	somme = somme / (double)n2;
	term3 = term3 * somme;
      }

      res[0] = term1 +  term2 + term3;
  }

  void Dcov2Cnormed(double *X, int *vecd, double *a, int *n, int *p, int *weightchoice, double *res, double *denom2) {
    
    double xijjprimel(int j, int jprime, int l, double *X, double *a, int *n, int *vecd, int *weightchoice);
    double term1=0.0, term2=0.0, somme, prod, xl;
    int j, jprime, l, n2=n[0]*n[0];


      for (j = 1; j <= n[0]; j++) {
	for (jprime = 1; jprime <= n[0]; jprime++) {
	  prod = 1.0;
	    for (l = 1; l <= p[0]; l++) {
	      prod = prod * xijjprimel(j,jprime,l,X,a,n,vecd,weightchoice);
	    }
	  term1 = term1 + prod;
	}
      }
      term1 = term1 / (double)n2;

      for (j = 1; j <= n[0]; j++) {
	prod = 1.0;
	for (l = 1; l <= p[0]; l++) {
	  somme = 0.0;
	  for (jprime = 1; jprime <= n[0]; jprime++) {
	    somme = somme + xijjprimel(j,jprime,l,X,a,n,vecd,weightchoice);
	  }
	  somme = somme / (double)(n[0]);
	  prod = prod * somme;
	}
	term2 = term2 + prod;
      }
      term2 = -2.0 * term2 / (double)(n[0]);

      prod = 1.0;
      somme = 0.0;
      for (l = 1; l <= p[0]; l++) {
	xl = 0.0;
	for (j = 1; j <= n[0]; j++) {
	  for (jprime = 1; jprime <= n[0]; jprime++) {
	    xl = xl + xijjprimel(j,jprime,l,X,a,n,vecd,weightchoice);
	  }
	}
	xl = xl / (double)n2;
	somme = somme + 1.0 / xl;
	prod = prod * xl;
      }

      denom2[0] = 1.0 - (1.0 + somme - p[0]) * prod;
      res[0] = term1 +  term2 + prod;
  }

  double xijjprimel(int j, int jprime, int l, double *X, double *a, int *n, int *vecd, int *weightchoice) {

    int indbeginblocl = 1;
    int i, k;
    double prod = 1.0, tmp = 0.0;

    if (l != 1) {
      for (i = 1; i <= (l-1); i++) {
	indbeginblocl = indbeginblocl + vecd[i-1];
      }
    }

    if (weightchoice[0] == 1) {
      for (k = 1; k <= vecd[l-1]; k++) {
	prod = prod * exp(-R_pow(a[0] * (X[(indbeginblocl+k-2)*n[0] + (j-1)] - X[(indbeginblocl+k-2)*n[0] + (jprime-1)]),2.0) / 2.0 );
	  }
    }
    if (weightchoice[0] == 2) {
      for (k = 1; k <= vecd[l-1]; k++) {
	prod = prod / (1.0 + R_pow((X[(indbeginblocl+k-2)*n[0] + (j-1)] - X[(indbeginblocl+k-2)*n[0] + (jprime-1)]) / a[0],2.0));
	  }
    }
    if (weightchoice[0] == 3) {
      prod = R_PosInf;
    }
    if (weightchoice[0] == 4) {
      for (k = 1; k <= vecd[l-1]; k++) {
	prod = prod * (-2.0*fabs(X[(indbeginblocl+k-2)*n[0] + (j-1)] - X[(indbeginblocl+k-2)*n[0] + (jprime-1)]) + fabs(X[(indbeginblocl+k-2)*n[0] + (j-1)] - X[(indbeginblocl+k-2)*n[0] + (jprime-1)] - 2.0*a[0]) + fabs(X[(indbeginblocl+k-2)*n[0] + (j-1)] - X[(indbeginblocl+k-2)*n[0] + (jprime-1)] + 2.0*a[0])) / (4.0*fabs(a[0]));
	  }
    }
    if (weightchoice[0] == 5) {
      for (k = 1; k <= vecd[l-1]; k++) {
	tmp = X[(indbeginblocl+k-2)*n[0] + (j-1)] - X[(indbeginblocl+k-2)*n[0] + (jprime-1)];
	prod = prod * (1.0 - R_pow(a[0]*tmp,2.0)) * exp(- R_pow(a[0]*tmp,2.0) / 2.0);
	  }
    }
    return prod;
  }



/* Function outputs union of arr1[] and arr2[]
   m is the number of elements in arr1[]
   n is the number of elements in arr2[] 
   arr3 should be of length n+m in input, arr3[0..k-1] will contain the union in output
   k is the card of the union
*/

int Union(int *arr1, int *arr2, int m, int n, int *arr3) {
  int i = 0, j = 0, k = 0;
  while(i < m && j < n) {
    if(arr1[i] < arr2[j]) {
      arr3[k] = arr1[i++];
      k = k + 1;
    } else if(arr2[j] < arr1[i]) {
      arr3[k] = arr2[j++];
      k = k + 1;
    } else {
      arr3[k] =arr2[j++];
      k = k + 1;
      i++;
    }
  }
 
  /* remaining elements of the larger array */
  while(i < m) {
    arr3[k] = arr1[i++];
    k = k + 1;
  }
  while(j < n) {
    arr3[k] = arr2[j++];
    k = k + 1;
  }
  return(k);
}

/* Function outputs Intersection of arr1[] and arr2[]
   m is the number of elements in arr1[]
   n is the number of elements in arr2[] 
   arr3 should be of length n in input, arr3[0..k-1] will contain the intersection in output
   k is the card of the intersection
*/
int Intersection(int *arr1, int *arr2, int m, int n, int *arr3) {
  int i = 0, j = 0, k = 0;
  while(i < m && j < n) {
    if(arr1[i] < arr2[j])
      i++;
    else if(arr2[j] < arr1[i])
      j++;
    else {/* if arr1[i] == arr2[j] */
      arr3[k] = arr2[j++];
      k = k + 1;
      i++;
    }
  }
  return(k);
}

/* Function outputs Setc difference of arr1[] and arr2[]
   m is the number of elements in arr1[]
   n is the number of elements in arr2[] , n<m
   arr3 should be of length m in input, arr3[0..k-1] will contain the symmetric difference in output
   j is the card of the set difference
*/
int Setdiff(int *arr1, int *arr2, int m, int n, int p, int *arr3) {
  int i = 0, j = 0; // Because elements in arr1 and arr2 are supposed sorted, we should be able to improve the code below!!
  int *tmp;
  tmp = new int[p];
  for (i = 0; i <p; i++) tmp[i] = 0;
  for (i = 0; i <m; i++) tmp[arr1[i]-1] = 1;
  for (i = 0; i <n; i++) tmp[arr2[i]-1] = 0;
  for (i = 0; i <p; i++) {
    if (tmp[i] == 1) arr3[j++] = i + 1;
  }
  delete[] tmp;
  return(j);
}

}
