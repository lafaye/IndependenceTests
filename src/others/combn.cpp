/* Debut des commentaires

Nom de la fonction: combn
-------------------

Auteur: Pierre Lafaye de Micheaux
-------

Date: 09/11/2005
-----

Entrées:
--------

    int *combmat: 
    int *n: 
    int *m: 

Sorties:
--------

    Le pointeur modifié est *combmat. La fonction ne renvoie rien.
    *combmat contiendra le vecteur constitué de la concaténation des colonnes de la matrice m x Cnm des combinaisons des elements de seq(n) pris m à la fois.

Fonctions extérieures appelées:
-------------------------------

    Cnp

Description:
------------

    Generate all combinations of the elements of seq(n) taken m at a time. C'est-à-dire une matrice de taille m x Cnm.


Références:
-----------

    Nijenhuis, A. and Wilf, H.S. (1978) Combinatorial Algorithms for 
    Computers and Calculators.  NY:  Academic Press.
 
Exemples:
---------

    combn(5,3)

Equivalent en R:
----------------

    require(combinat)
    combn(5,3)

Instructions de compilation pour utilisation dans un terminal:
--------------------------------------------------------------

g++ -fPIC  -O2 -march=i686 -fomit-frame-pointer combn.cpp 
./a.out

Instructions de compilation pour utilisation depuis R:
------------------------------------------------------

g++ -I/usr/lib/R/include  -I/usr/local/include  -fPIC  -O2 -march=i686 -fomit-frame-pointer -c combn.cpp -o combn.o
g++ -shared -L/usr/local/lib -o combn.so combn.o

Pour utiliser dans R, taper source("combn.R") où le fichier normeRnAwR.R contient le code R suivant:

combn <- function(n,m) {

combmat<-matrix(0,nrow=m,ncol=choose(n,m))

dyn.load(paste("combn", .Platform$dynlib.ext, sep="")) 
out <- .C("combn",res=as.integer(combmat),as.integer(n),as.integer(m))
result<-matrix(out$res,ncol=m,byrow=F)
return(result)
dyn.unload(paste("combn", .Platform$dynlib.ext, sep="")) 
}


Utilisation du debugger gdb avec electric fence:
------------------------------------------------

g++ -g combn.cpp -lm -u malloc -lefence
gdb ./a.out
run

Utilisation du debugger ddd:
----------------------------

g++ -I/usr/lib/R/include  -I/usr/local/include  -c combn.cpp -o combn.o -g
g++ -shared -L/usr/local/lib -o combn.so combn.o 
R -d ddd

Menu Program: cocher Run in Execution Window
Menu Program: Run puis Run

source("combn.R")

Menu File/Open source...
Cliquer sur Load Shared Object Library Symbols
Sélectionner combn.cpp
Cliquer sur Open
Mettre des breakpoints
Dans la fenêtre Execution Window de R, taper: 

combn(5,3)

Fin des commentaires */


// Inclusion de librairies et de fonctions extérieures
//----------------------------------------------------

#include <iostream>
using namespace std;
#include <math.h>


// Utilisation dans une fonction main:
// -----------------------------------

/*

#include "Cnp.cpp"

extern "C" {


   int main()

  {

    void combn(int *combmat, int *n, int *m);

    int *n, *m, *combmat, i, j;

    double Cnm;

    n = new int[1];
    m = new int[1];

    *(n+0)=5;
    *(m+0)=3;

    Cnm=Cnp(*(n+0),*(m+0));

    combmat = new int[(int)Cnm**(m+0)];

    combn(combmat,n,m);

    for (j = 1; j <= Cnm**(m+0); j++)      cout << *(combmat+j-1) << " ";
    
    cout << "\n";
    cout << "\n";
    
    for (j = 1; j <= Cnm; j++) {
      
      for (i = 1; i <= *(m+0); i++) {
	
	cout << *(combmat+(j-1)**(m+0)+i-1) << " ";
      }
      
      cout << "\n";
    }
    
    cout << "\n";
    
    for (i = 1; i <= *(m+0); i++) {
      
      for (j = 1; j <= Cnm; j++) {
	
	
	cout << *(combmat+(j-1)**(m+0)+i-1) << " ";
	
      }
      cout << "\n";
      
    }

  }

} // extern C

*/


//-----------------------------------------------------------------
// DEBUT DE LA FONCTION
//-----------------------------------------------------------------


extern "C" {



void combn(int *combmat, int *n, int *m)
{

  
  int i, j, e, h, nmmp1, mp1;
	
  int *a;
  a=new int[*(m+0)];
  for (i=1;i<=*(m+0);i=i+1) *(a+i-1)=i;

  e=0;
  h=*(m+0);
  
  for (i=1;i<=*(m+0);i=i+1) *(combmat+i-1)=i; 
	
  i=2;
  nmmp1=*(n+0) - *(m+0) + 1;
  mp1=*(m+0) + 1;
	while(*(a+0) != nmmp1) {

		if(e < *(n+0) - h) {
		  h=1;
		  e=*(a+*(m+0)-1);
			

		  *(a+*(m+0) - h)=e + 1;

		  for (j=1;j<=*(m+0);j=j+1) *(combmat+(i-1)**(m+0)+j-1)=*(a+j-1);
		  
		i=i+1;

		}
		else { 
		  h=h + 1;
		  e=*(a+mp1 - h-1);
		  
		  for (j=1;j<=h;j=j+1) *(a+*(m+0) - h + j-1)=e + j;
		  
		for (j=1;j<=*(m+0);j=j+1) *(combmat+(i-1)**(m+0)+j-1)=*(a+j-1); 

		i=i + 1;

			  }

	}
	//On libère de la mémoire
	delete[] a;
}


} // extern C

