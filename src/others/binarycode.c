/* Beginning comments

Function name: binarycode

Entries: i integer, L integer, *y pointer

Output: void

Description:
Translate into binary over L bits the integer i.
Results are retrieved via pointer *y

Author: Pierre Lafaye de Micheaux

Date: 30/05/2014

End comments */


void binarycode(int i,int L, int *y) {

   int n, j;
   int R, q, *x;
   
   x = new int [L];
   q = i;
   n = 0;

   do {
     R = q - 2 * (int)floor(q / 2.0);
     q = (int)floor(q / 2.0);
     for (j=L-1; j>0; j=j-1) x[j] = x[j - 1];
     x[0] = R;
     n = n+1;
   }

   while (q>0);

   if (n == L) 
     for (j=0; j<L; j=j+1) y[j] = x[j];
   else if (L-n > 0) {
    for (j=0; j<L-n; j=j+1) y[j] = 0;
    for (j=L-n; j<L; j=j+1) y[j] = x[n - L + j];
   }
   else {
     // cout << "Choose a larger value for L.\n"; 
     y[0] = L - n;
   }
 
}
