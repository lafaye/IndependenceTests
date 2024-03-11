#define pi 3.141592653589793115998
         
void gauher(double b, int N, double prec, int MAXIT, double *weightsandpoints) {

/*

Returns a N x 2 matrix (via weightsandpoints pointer).
First column contains points.
Second colonne contains  weights.

*/

  double EPS, pim4, z, p1, p2, p3, pp, z1;
  int its, i, j;
  
  //  double var = 0.0;
  z = 0.0;
  pp = 0.0;
  EPS = prec;
  pim4 = 1.0/(pow(pi,0.25));
        
  double *x, *w;
  x = new double[N];
  w = new double[N];
  for (i=0; i<N; i++) {x[i] = 0.0; w[i] = 0;}
        
  int m = (int)((N+1) / 2.0);

  for (i=1; i<=m; i=i+1) {

      if (i == 1) { z = (sqrt((double)(2 * N + 1)) - 1.85575 * (pow((double)(2 * N + 1),(double)(-1.0 / 6.0))));}
      else { if (i == 2) {z = (z - 1.14 * (pow(N,0.426)) / z);}
      else { if (i == 3) {z = (1.86 * z - 0.86 * x[0]);}
      else { if (i == 4) {z = (1.91 * z - 0.91 * x[1]);}
      else { z = (2 * z - x[i-3]);} 
      }}}

      //var = 0.0;
        
      for (its=1; its<=MAXIT; its=its+1) {
	  p1 = pim4;
	  p2 = 0.0;
	  for (j=1; j<=N; j=j+1) {
	    p3 = p2;
	    p2 = p1;
	    p1 = z * (sqrt(2.0 / j)) * p2 - sqrt((j - 1.0) / (double)j) * p3;
	  }

	  pp = sqrt((double)(2 * N)) * p2;
	  z1 = z;
	  z = z1 - p1 / pp;
	  if (fabs(z-z1) < EPS) {//var = 1.0;
                          break;}
	}
                
      //      if (var == 0.0) cout << "Too many iterations\n";

      x[i-1] = z;
      x[N-i] = -z;
      w[i-1] = 2.0 / (pp*pp);
      w[N-i] = w[i-1];
    }

  for (i=0; i<N; i++) {
    weightsandpoints[i] = sqrt(2.0) * b * x[i];
    weightsandpoints[N + i] = w[i] / sqrt(pi);
  }
        

  delete[] x;
  delete[] w;

  return;
}
