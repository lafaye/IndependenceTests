/*

Nom de la fonction: Cnp
-------------------

Entrées:
--------

 int n et int k

Sorties:
--------

 Returns the binomial coefficient Cnk as a double number.

Fonctions extérieures appelées:
-------------------------------

gammln, factrl, factln, 

Références:
-----------
Numerical recipies in C chapter 6.1

*/

extern"C" {

#include <math.h>
double gammln(double xx)
     // Returns the value ln[Gamma(xx)] for xx > 0.
{
  // Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure accuracy is good enough.
double x,y,tmp,ser;
static double cof[6]={76.18009172947146,-86.50532032941677,
24.01409824083091,-1.231739572450155,
0.1208650973866179e-2,-0.5395239384953e-5};
int j;
y=x=xx;
tmp=x+5.5;
tmp -= (x+0.5)*log(tmp);
ser=1.000000000190015;
for (j=0;j<=5;j++) ser += cof[j]/++y;
return -tmp+log(2.5066282746310005*ser/x);
}

double factrl(int n)
     // Returns the value n! as a double number.
{
double gammln(double xx);
static int ntop=4;
 static double a[33]={1.0,1.0,2.0,6.0,24.0}; //Fill in table only as required.
int j;
if (n > 32) return exp(gammln(n+1.0));
// Larger value than size of table is required. Actually, this big a value is going to overflow
// on many computers, but no harm in trying.
 while (ntop<n) { // Fill in table up to desired value.
j=ntop++;
a[ntop]=a[j]*ntop;
}
return a[n];
}

double factln(int n)
     // Returns ln(n!).
{
double gammln(double xx);
 static double a[101]; // A static array is automatically initialized to zero.
if (n <= 1) return 0.0;
 if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0)); // In range of table.
							   else return gammln(n+1.0); // Out of range of table.
}

double Cnp(int n, int k)
     // Returns the binomial coefficient Cnk as a double number.
{
double factln(int n);
return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
// The floor function cleans up roundoff error for smaller values of n and k.
}

}
