#include <R.h>
#include "Rmath.h"
#include "mycomplex.h"
#include <iostream>
using namespace std;

extern"C" {

  void CnhatC(double *vecs, double *vect, double *X, int *n, int *q, int *p, int *vecd, std::complex<double> *res) {

    void phinhatReturn(double *vect1, double *vect2, double *vect3, double *X, int *q, int *n, cmplx *res1, cmplx *res2, cmplx *res3);
    cmplx prod1 = C(1.0, 0.0), prod2 = C(1.0, 0.0), somme = C(0.0, 0.0);
    cmplx tmp2;
    int indbeginblocl, l, i, *vecdl, k, j;
    double *vecsl, *mvectl, *diffvecsvectl, *Xl;
    cmplx *res1, *res2, *res3;
    vecdl = new int[1];
    res1 =  new cmplx[1];
    res2 =  new cmplx[1];
    res3 =  new cmplx[1];

    for (l = 0; l <= (p[0]-1); l++) {
      indbeginblocl = 1;
      if (l != 0) {
	for (i = 0; i <= (l-1); i++) {
	  indbeginblocl = indbeginblocl + vecd[i];
	}
      }
      vecdl[0] = vecd[l];
      vecsl = new double[vecdl[0]];
      for (k = 0; k <= (vecdl[0]-1); k++) vecsl[k] = vecs[indbeginblocl+k-1];
      mvectl = new double[vecdl[0]];
      diffvecsvectl = new double[vecdl[0]];
      for (k = 0; k <= (vecdl[0]-1); k++) {
	mvectl[k] = - vect[indbeginblocl+k-1];
	diffvecsvectl[k] = vecsl[k] + mvectl[k];
      }
      Xl = new double[n[0]*vecdl[0]];
      for (j = 0; j <= (n[0]-1); j++) {
	for (k = 0; k <= (vecdl[0]-1); k++) {
	  Xl[k*n[0] + j] = X[(indbeginblocl+k-1)*n[0] + j];
	}
      }
      phinhatReturn(diffvecsvectl,vecsl,mvectl,Xl,vecdl,n,res1,res2,res3);

      prod1 = prod1 * res1[0];
      tmp2 = res2[0] * res3[0];
      prod2 = prod2 * tmp2;
      somme = somme + res1[0] / tmp2;

      delete[] vecsl;
      delete[] mvectl;
      delete[] diffvecsvectl;
      delete[] Xl;
    }
    delete[] vecdl;
    delete[] res1;
    delete[] res2;
    delete[] res3;
    res[0] = prod1 - prod2 * (cmplx)(1 - (double)p[0] + somme);

  }

  void CnhatmatC(double *yMat, int *N, double *X, int *n, int *q, int *p, int *vecd, cmplx *res) {

    int i, j, k, qm = q[0] - 1, Nm = N[0] - 1;
    double *vect, *vecs;
    cmplx *restmp;
    restmp = new cmplx[1];
    vecs = new double[q[0]];
    vect = new double[q[0]];

    for (i = 0; i <= Nm; i++) {
      for (j = 0; j <= i; j++) {
	for (k = 0; k <= qm; k++) {
	  vecs[k] = yMat[k*N[0] + i];
	  vect[k] = yMat[k*N[0] + j];
	}
	CnhatC(vecs,vect,X,n,q,p,vecd,restmp);
	res[j*N[0] + i] = restmp[0];
      }

    }

    for (i = 0; i <= (N[0]-2); i++) {
      for (j = i+1; j <= Nm; j++) {
	res[j*N[0] + i] = conj(res[i*N[0] + j]);
      }
    }

    delete[] restmp;
    delete[] vecs;
    delete[] vect;
  }

  void CnhatmatClower(double *yMat, int *N, double *X, int *n, int *q, int *p, int *vecd, cmplx *res, double *sqrtweights) {

    int i, j, k, l, qm = q[0] - 1, Nm = N[0] - 1;
    double *vect, *vecs;
    cmplx *restmp;
    restmp = new cmplx[1];
    vecs = new double[q[0]];
    vect = new double[q[0]];

    l = 0;
    for (j = 0; j <= Nm; j++) {
      for (k = 0; k <= qm; k++) {
	vect[k] = yMat[k*N[0] + j];
      }
      for (i = j; i <= Nm; i++) {
	for (k = 0; k <= qm; k++) {
	  vecs[k] = yMat[k*N[0] + i];
	}
	CnhatC(vecs,vect,X,n,q,p,vecd,restmp);
	res[l] = restmp[0] * sqrtweights[i] * sqrtweights[j];
	l = l + 1;
      }
    }

    delete[] restmp;
    delete[] vecs;
    delete[] vect;
  }
  
  void phinhatReturn(double *vect1, double *vect2, double *vect3, double *X, int *q, int *n, cmplx *res1, cmplx *res2, cmplx *res3) {
    int j, col, qm = q[0]-1;
    double tmp1, tmp2, tmp3, tmp;
    cmplx I = C(0.0, 1.0);
    res1[0] = C(0.0, 0.0);
    res2[0] = C(0.0, 0.0);
    res3[0] = C(0.0, 0.0);

    for (j = 0; j<= (n[0]-1); j++) {
      tmp1 = 0.0;
      tmp2 = 0.0;
      tmp3 = 0.0;
      for (col = 0; col <= qm; col++) {
	tmp = X[col * n[0] + j];
	tmp1 = tmp1 + vect1[col] * tmp;
	tmp2 = tmp2 + vect2[col] * tmp;
	tmp3 = tmp3 + vect3[col] * tmp;
      }      
      res1[0] = res1[0] + cexp(I * tmp1);
      res2[0] = res2[0] + cexp(I * tmp2);
      res3[0] = res3[0] + cexp(I * tmp3);
    }
    res1[0] = res1[0] / (double)n[0];
    res2[0] = res2[0] / (double)n[0];
    res3[0] = res3[0] / (double)n[0];
  }
  

}
