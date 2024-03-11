#include <R.h>
#include "Rmath.h"
#include "mycomplex.h"

extern"C" {

  void phinhatC(double *vect, double *X, int *q, int *n, std::complex<double> *res) {

    int j, col;
    double tmp = 0.0;
    cmplx somme = C(0.0, 0.0);
    cmplx I = C(0.0, 1.0);

    for (j = 1; j<= n[0]; j++) {
      tmp = 0.0;
      for (col = 1; col <= q[0]; col++) {
	tmp = tmp + vect[col-1] * X[(col-1) * n[0] + (j-1)];
      }
      somme = somme + cexp(I * tmp);
    }
    res[0] = somme / (double)n[0];
  }

}



