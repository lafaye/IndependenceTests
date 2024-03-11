/* Debut des commentaires

Nom de la fonction: normeRnAwR
-------------------

Auteur: Pierre Lafaye de Micheaux
-------

Date: 09/11/2005
-----

Entrées:
--------
    int *N: nombre de points de la discrétisation de (0,PI)
    int *vecd: (d_1,d_2,...,d_p) où les d_j sont des entiers
    int *p: nombre de vecteurs aléatoires dont on veut tester l'indépendance
    double *X: matrice des données: n lignes et sum(vecd) colonnes
    int *n: nombre d'individus
    int *q: nombre de colonnes de la matrice X
    double *res: contiendra les résultats
    double *compt: affichage d'un compteur - si 0 pas de compteur, sinon il faut mettre 
        \sum_{k=2}^p 2^k N^{-k} \sum_{i=1}^{C_p^k} N^{\sum_{l=1}^k vecd[a_{i,l}]} où A_i=(a_{i,1}, ... , a_{i,k}) est le sous-ensemble (de taille k) de {1,..,p} sélectionné à chaque pas i

Sorties:
--------

    Le pointeur modifié est *res. La fonction ne renvoie rien.

Fonctions extérieures appelées:
-------------------------------

    combn
    Cnp

Description:
------------

Cette fonction calcule une approximation (due fait des discrétisations données) de ||R_{n,A}|| pour tous les 2^p-p-1 choix de A (cardA>1) possibles.

Références:
-----------

Exemples:
---------

N<-10
vecd.or.p<-c(3,2)
X<-matrix(c(0.7,0.4,0.1,0.8,0.4,0.13,3.7,0.4,2.1,1.8,0.1,0.1,0.7,0.2,0.1),byrow=F,nrow=3)
x<-Sys.time();normeRnAwR(X,vecd.or.p,N);Sys.time()-x

N<-5
vecd.or.p<-c(3,2,6)
X<-matrix(c(0.7,0.4,0.1,0.8,0.4,0.13,3.7,0.4,2.1,1.8,0.1,0.1,0.7,0.2,0.1,0.45,0.2,0.1,2,3.4,5.6,1.2,0.8,2.6,6.4,1.1,2.2,3.3,4.1,2.5,1.9,2.7,3.5),byrow=F,ncol=11)
x<-Sys.time();normeRnAwR(X,vecd.or.p,N,compt=0,affiche=1,seriel=0);Sys.time()-x
#On peut aussi donner un décompte, mais cela rallonge beaucoup le temps de calcul ...
x<-Sys.time();normeRnAwR(X,vecd.or.p,N,compt=1,affiche=1,seriel=0);Sys.time()-x


source("normeRnAwR.R")
N<-5
vecd.or.p<-4
X<-matrix(runif(30),nrow=10,ncol=3)
normeRnAwR(X,vecd.or.p,N,compt=0,affiche=1)

Equivalent en R:
----------------

normeRnAwRpurR<-function(X,vecd.or.p,N,compt=0,affiche=1) {

if (length(vecd.or.p) > 1) {
#on fait le cas non sériel
vecd<-vecd.or.p

p<-length(vecd)
resultat<-rep(0,2^p-p-1)

require(combinat)
RES<-as.list(1:(p-1))
for (cardA in 2:p) {RES[[cardA]]<-as.matrix(combn(p,cardA))}
nb<-0
for (cardA in 2:p) {
for (j in 1:(choose(p,cardA))) {
nb<-nb+1
resultat[nb]<-normeRnApurR(X,vecd,N,A=RES[[cardA]][,j],compt)
if (affiche != 0) {cat(c("A:",RES[[cardA]][,j],"||R_{n,A}||:",round(resultat[nb],3),"\n"))}
}
}
return(resultat)
}


if (length(vecd.or.p)==1) {
#On fait le cas sériel
#Attention, on a considéré que n est grand par rapport à p et donc que l'on peut calculer R_{n,A} à la place de S_{n,A}
#mais sur la matrice X de taille nprime x nprime ci-dessous
p<-vecd.or.p
Y<-X
n<-nrow(Y)
q<-ncol(Y)
vecd<-rep(q,p)
nprime<-n-p+1
X<-matrix(0,nrow=nprime,ncol=p*q)
for (i in 1:nprime) {
for (j in 1:p) {

X[i,(((j-1)*q+1):(j*q))] <- Y[i+j-1,]}}



resultat<-rep(0,2^(p-1)-1)


require(combinat)
RES<-as.list(1:(p-1))
for (cardA in 2:p) {RES[[cardA]]<-as.matrix(rbind(rep(1,choose(p-1,cardA-1)),as.matrix(combn(p-1,cardA-1)+1)))}
nb<-0
for (cardA in 2:p) {
for (j in 1:(choose(p-1,cardA-1))) {
nb<-nb+1
resultat[nb]<-normeRnApurR(X,vecd,N,A=RES[[cardA]][,j],compt)
if (affiche != 0) {cat(c("A:",RES[[cardA]][,j],"||R_{n,A}||:",round(resultat[nb],3),"\n"))}

}

}

return(resultat)
}
} 


normeRnApurR<-function(X,vecd,N,A,compt){
#Entrees:
#X est la matrice des données: n lignes et sum(vecd) colonnes
#A est un sous ensemble de {1,...,p}
#vecd=(d_1,d_2,...,d_p) où les d_j sont des entiers
#Sorties:
#Cette fonction calcule, pour un A donné, une approximation (du fait des discrétisations données) de ||R_{n,A}||
p<-length(vecd)
n<-nrow(X)
#On calcule les points sur les sphères
lesdirections<-directionspurR(N,vecd)
matrice<-matrix(0,nrow=n,ncol=n)
#On commence par calculer toutes les permutations des directions possibles s^{(1)},...,s^{(p)}
resultat<-permutpurR(2*N^(vecd-1))
globmaxRnA<-c()
#Ensuite, pour chaque jeu de directions s^{(1)},...,s^{(p)}, on calcule la matrice Psi_A, le vecteur psi_A et le max de R_{n,A} en t comme indiqué par Martin
for (j in 1:nrow(resultat)) {
#Affichage d'un compteur à l'écran pour savoir où on en est
if (compt != 0) {cat(" ",nrow(resultat)-j," ")}
#Initialisation de la matrice PsiA
PsiA<-matrix(1,nrow=n,ncol=n)
#On se donne l'identificateur d'un jeu de directions s^{(1)},...,s^{(p)}
ident<-resultat[j,]
for (k in A) {
#C'est le point s^{(k)}
sk<-lesdirections[[k]][ident[k],]
debut<-(sum(vecd[1:(k-1)])+1)
if (k==1) debut<-1
fin<-(sum(vecd[1:k]))
for (i in 1:n) {
#On calcule t_i^{(k)} comme indiqué par Martin
tik<-sum(X[i,debut:fin]*sk)
for (l in 1:n) {
#On remplit la matrice n*n comme indiqué par Martin
matrice[l,i]<-indHpurR(sk,tik,X[l,debut:fin])
}
}
#On multiplie (au fur et à mesure) les |A| matrices n*n obtenues 
PsiA<-PsiA*scale(matrice,scale=F)
}
psiA<-apply(PsiA,FUN=sum,MARGIN=2)
maxRnA<-max(abs(psiA))/(sqrt(n))
globmaxRnA<-c(globmaxRnA,maxRnA)
}
#Le max global, c'est-à-dire ||R_{n,A}|| si la discretisation est assez fine
resultat<-max(globmaxRnA) 
return(resultat)
}

Instructions de compilation pour utilisation dans un terminal:
--------------------------------------------------------------

g++ -Wall normeRnAwR.cpp 
./a.out

g++ -Wall -fPIC  -O2 -march=i686 -fomit-frame-pointer normeRnAwR.cpp 
./a.out

Instructions de compilation pour utilisation depuis R:
------------------------------------------------------

g++ -I/usr/lib/R/include  -I/usr/local/include  -fPIC  -O2 -march=i686 -fomit-frame-pointer -c normeRnAwR.cpp -o normeRnAwR.o
g++ -shared -L/usr/local/lib -o normeRnAwR.so normeRnAwR.o 

Pour utiliser dans R, taper source("normeRnAwR.R") où le fichier normeRnAwR.R contient le code R suivant:
 
normeRnAwR <- function(X,vecd.or.p,N=10,compt=0,affiche=1) {

#si length(vecd.or.p)>1 alors cas non sériel sinon cas sériel



	     dyn.load(paste("normeRnAwR", .Platform$dynlib.ext, sep=""))
           if (!is.numeric(vecd.or.p))
                   stop("argument vecd must be numeric")



if (length(vecd.or.p) > 1) {
#on fait le cas non sériel
seriel<-0
vecd<-vecd.or.p
#Laisser compt=0 si il y a beaucoup de sous-vecteurs ou s'ils sont de trop grande taille
if (compt != 0) {

resultat<-permutR(2*N^(vecd-1))
compt<-nrow(resultat)

}

p<-length(vecd)
resultat<-rep(0,2^p-p-1)

	
		     out <- .C("normeRnAwR",
		     as.integer(N),
		     as.integer(vecd),
		     as.integer(length(vecd)),
		     as.integer(p),
		     as.numeric(X),
		     as.integer(nrow(X)),
		     as.integer(ncol(X)),
			resultat=as.numeric(resultat),
		as.numeric(compt),
		as.integer(seriel))

if (affiche != 0) {
require(combinat)
RES<-as.list(1:(p-1))
for (cardA in 2:p) {RES[[cardA]]<-as.matrix(combn(p,cardA))}
nb<-0
for (cardA in 2:p) {
for (j in 1:(choose(p,cardA))) {
nb<-nb+1
cat(c("A:",RES[[cardA]][,j],"||R_{n,A}||:",round(out$resultat[nb],3),"\n"))
}
}
}
#cela renvoie un vecteur de taille 2^p-p-1 contenant la norme de RnA pour chacun des 2^p-p-1 A différents de |A|>1
return(out$resultat)
}


if (length(vecd.or.p)==1) {
#On fait le cas sériel
#Attention, on a considéré que n est grand par rapport à p et donc que l'on peut calculer R_{n,A} à la place de S_{n,A}
#mais sur la matrice X de taille nprime x nprime ci-dessous
seriel<-1
p<-vecd.or.p
vecd<-rep(ncol(X),p)


#Laisser compt=0 si il y a beaucoup de sous-vecteurs ou s'ils sont de trop grande taille
if (compt != 0) {

compt<-nrow(permutR(2*N^(vecd-1)))

}

resultat<-rep(0,2^(p-1)-1)

	
		     out <- .C("normeRnAwR",
		     as.integer(N),
		     as.integer(vecd),
		     as.integer(length(vecd)),
		     as.integer(p),
		     as.numeric(X),
		     as.integer(nrow(X)),
		     as.integer(ncol(X)),
			resultat=as.numeric(resultat),
		as.numeric(compt),
		as.integer(seriel))

if (affiche != 0) {
require(combinat)
RES<-as.list(1:(p-1))
for (cardA in 2:p) {RES[[cardA]]<-as.matrix(rbind(rep(1,choose(p-1,cardA-1)),as.matrix(combn(p-1,cardA-1)+1)))}
nb<-0
for (cardA in 2:p) {
for (j in 1:(choose(p-1,cardA-1))) {
nb<-nb+1
cat(c("A:",RES[[cardA]][,j],"||R_{n,A}||:",round(out$resultat[nb],3),"\n"))
}
}
}
#cela renvoie un vecteur de taille 2^(p-1)-1 contenant la norme de RnA pour chacun des 2^p-p-1 A différents de |A|>1 et qui contiennent 1
return(out$resultat)
}




		     dyn.unload(paste("normeRnAwR", .Platform$dynlib.ext, sep=""))
 
}



permutR<-function(k) {

#Ce programme calcule toutes les permutations du vecteur k=(k_1,...,k_p) où p=length(k)
#C'est-à-dire tous les vecteurs différents de longueur p depuis (1,...,1) jusqu'à (k_1,...,k_p)

  x<-matrix(1:k[1],ncol=1)
  for(i in k[-1]) x<-cbind(x[rep(1:nrow(x),rep(i,nrow(x))),],rep(1:i,nrow(x)))
  x 
}



Utilisation du debugger gdb avec electric fence:
------------------------------------------------

g++ -g normeRnAwR.cpp -lm -u malloc -lefence
gdb ./a.out
run

Utilisation du debugger ddd:
----------------------------

g++ -I/usr/lib/R/include  -I/usr/local/include  -c normeRnAwR.cpp -o normeRnAwR.o -g
g++ -shared -L/usr/local/lib -o normeRnAwR.so normeRnAwR.o 
R -d ddd

Menu Program: cocher Run in Execution Window
Menu Program: Run puis Run


dyn.load(paste("normeRnAwR", .Platform$dynlib.ext, sep=""))

Menu File/Open source...
Cliquer sur Load Shared Object Library Symbols
Sélectionner normeRnAwR.cpp
Cliquer sur Open
Mettre des breakpoints
Dans la fenêtre Execution Window de R, taper: 


source("normeRnAwR.R")
N<-10
vecd.or.p<-c(3,2)
X<-matrix(c(0.7,0.4,0.1,0.8,0.4,0.13,3.7,0.4,2.1,1.8,0.1,0.1,0.7,0.2,0.1),byrow=F,nrow=3)
normeRnAwR(X,vecd.or.p,N)

source("normeRnAwR.R")
N<-5
vecd.or.p<-3
X<-matrix(runif(30),nrow=10,ncol=3)
normeRnAwR(X,vecd.or.p=4,N,compt=0,affiche=1)


Fin des commentaires */



// Inclusion de librairies et de fonctions extérieures
//----------------------------------------------------

#include <iostream>
using namespace std;
#include <math.h>

#include "combn.cpp"
#include "Cnp.cpp"


// Utilisation dans une fonction main:
// -----------------------------------

/*

extern "C" {


   int main()

  {

    void normeRnAwR(int *N, int *vecd, int *lenvecd, int *p, double *X, int *n, int *q, double *res, double *compt, int *seriel);

    //Déclaration des variables
    int *N;
    int *vecd, *p;
    int *n, nbcol;
    double *X;
    double *res;
    double *compt;
    int i;
    int cardA, CpcardA;
    int numero;
    int *cardApoint;
    int j;
    int *combmat;
    int *seriel;

    //Initialisation des variables
    seriel=new int[1];
    *(seriel+0)=0;
    N=new int[1];
    *(N+0)=5;     //N: nombre de points de la discrétisation
    p=new int[1];
    *(p+0)=3;     //longueur de vecd
    vecd = new int[*(p+0)];
    *(vecd+0)=1;     //vecd=(d_1,d_2,...,d_p) où les d_j sont des entiers
    *(vecd+1)=1;   
    *(vecd+2)=1;   
    n=new int[1];
    *(n+0)=5;
    nbcol=3;
    X = new double[*(n+0)*nbcol];     //X est la matrice des données: n lignes et sum(vecd)=nbcol colonnes, remplie par colonnes
    *(X+0)=0.7;*(X+1)=0.4;*(X+2)=0.1;*(X+3)=0.8;*(X+4)=0.4;*(X+5)=0.13;
    *(X+6)=3.7;*(X+7)=0.4;*(X+8)=2.1;*(X+9)=1.8;
    *(X+10)=0.1;*(X+11)=0.1;*(X+12)=0.7;*(X+13)=0.2;*(X+14)=0.1;    
    res=new double[(int)pow(2.0,*(p+0))-*(p+0)-1];
    for (i = 1; i <= ((int)pow(2.0,*(p+0))-*(p+0)-1); i ++) *(res+i-1)=999999;
    compt=new double[1];
    //Si on met compt à 0, alors le décompte n'est pas affiché, sinon il faut mettre la valeur suivante:
    // *(compt+0)=\sum_{k=2}^p 2^k N^{-k} \sum_{i=1}^{C_p^k} N^{\sum_{l=1}^k vecd[a_{i,l}]} où A_i=(a_{i,1}, ... , a_{i,k}) est le sous-ensemble (de taille k) de {1,..,p} sélectionné 
    // à chaque pas i.
    *(compt+0)=0;
    numero=0;
    cardApoint = new int[1];
    *(cardApoint+0)=999999;

    //Appel de la fonction
    normeRnAwR(N, vecd, p, X, n, res, compt,seriel);

    //Affichage des résultats
    cout << "\nVecteur des maximum pour les différents A (pour |A| croissant et à |A| fixé pour la permutation croissante (par ordre lexicographique)):\n";
    cout << "\n";

    

    // On fixe |A|
    for (cardA=2; cardA<=*(p+0); cardA++) { 

      // On calcule Cp,cardA
      CpcardA=(int)Cnp(*(p+0),cardA);

      // On initialise combmat
      combmat=new int[CpcardA*cardA];
      for (i = 1; i <= CpcardA*cardA; i++) *(combmat+i-1)=999999;

      *(cardApoint+0)=cardA;

      // combmat va contenir toutes les permutations de 1:p de taille cardA, c'est-à-dire tous les ensembles A de taille |A|
      combn(combmat, p, cardApoint);
      

      // Pour chaque valeur de cardA, il y a Cp,cardA façons de choisir un ensemble de taille cardA parmi p valeurs
      // On va fixer A, pour |A| donné
            for (i = 1; i<=CpcardA; i++) {	

	numero=numero+1;
	

	// A est un sous ensemble (de taille cardA) qui contient la ième permutation de 1:p de taille cardA. C'est le i-ème ensemble A de taille |A|
	cout << "A= ";

			for (j=1; j<=cardA; j++) {
		  		  
		  cout << *(combmat+(i-1)*cardA+j-1) << " ";
		  }
	
	

	cout << "||R_{n,A}||= " << *(res+numero-1);

	
	 
	
	

	      } // fin de for i = 1 to CpcardA
      
      

	//On libère de la mémoire
      delete[] combmat;
      

      
    } // fin de for cardA = 2 to p

     cout << "\n";
   
    return(0);

  }

} // extern C

*/

//-----------------------------------------------------------------
// DEBUT DE LA FONCTION
//-----------------------------------------------------------------

extern "C" {


  void normeRnAwR(int *N, int *vecd, int *lenvecd, int *p, double *X, int *n, int *q, double *res, double *compt, int *seriel) 

  {


 
 
    //Déclaration des variables
    int taille;
    if (*(seriel+0) == 0) {
      taille=(int)pow(2.0,(double)*(p+0))-*(p+0)-1;
    }

    else {

      if (*(seriel+0) != 1) {return;}
      taille=(int)pow(2.0,(double)*(p+0)-1.0)-1;
    }
    

    int i, g;
    const double pi=3.1415926535897931160;
    int maxvecd, sommevecd;
    int cardA, CpcardA;
    double *PsiA;
    double *matrice;
    double x, produit, reste;
    int quotient;
    double prod;
    double *sk;
    int ii,k, j, jj, l, debut, fin;
    double tik, moyennecoli;
    int ind;
    double sommePsiA, maxpsiA;
    int *combmat;
    int *A, *cardApoint, *ppoint;
    double *permvecd, *limit;
    double max;
    double *maxRnA;
    int *decomp;
    int numero;
    double *discretization1, *discretization2;
    double somme;
    //On calcule max(vecd) et sum(vecd)
    maxvecd=*(vecd+0);
    sommevecd=*(vecd+0);
    for (i=1; i<=*(lenvecd+0)-1;i++) {
      if (maxvecd<*(vecd+i)) maxvecd=*(vecd+i);
      sommevecd=sommevecd+*(vecd+i);
    }


    double *newX;
    int nprime; 
    int *nbis;
    nbis = new int[1];

    if (*(seriel+0) == 0) {
      newX = new double[*(n+0)*sommevecd];
      for (i = 1; i <= *(n+0)*sommevecd; i++) *(newX+i-1)=*(X+i-1);
      *(nbis+0)=*(n+0);
    }

    else {
      if (*(seriel+0) != 1) {return;}
      nprime=*(n+0)-*(p+0)+1;
      newX = new double[nprime*sommevecd];


      for (j = 1; j <= *(p+0); j++) {
	  for (k = 1; k <= *(q+0); k++) {
	    for (i = 1; i <= nprime; i++) {

	      	      *(newX+(j-1)*nprime**(q+0)+(k-1)*nprime+i-1)=*(X+(k-1)*(*n+0)+i+j-1-1);
		      //a bien vérifier: on stocke une matrice à n lignes dans une matrice à nprime lignes: c'est bon??
	      
	  }
	}
      }

      *(nbis+0)=nprime;

      //Y<-X
      //n<-nrow(Y)
      //q<-ncol(Y)

      //nprime<-n-p+1
      //X<-matrix(0,nrow=nprime,ncol=p*q)
      //for (i in 1:nprime) {
      //for (j in 1:p) {

      //X[i,(((j-1)*q+1):(j*q))] <- Y[i+j-1,]}}

      //    }


    

    }

    //Initialisation des variables
    reste = 0.0;
    discretization1=new double[*(N+0)];
    discretization2=new double[2**(N+0)];
    PsiA = new double[*(nbis+0)**(nbis+0)];
    for (i = 1; i<= *(nbis+0)**(nbis+0); i++) *(PsiA+i-1)=999999;
    matrice = new double[*(nbis+0)**(nbis+0)];
    for (i = 1; i<= *(nbis+0)**(nbis+0); i++) *(matrice+i-1)=999999;
    sk = new double[maxvecd];
    for (i = 1; i <= maxvecd; i++) *(sk+i-1)=999999;
    cardApoint = new int[1];
    ppoint = new int[1];
    *(cardApoint+0)=999999;
    *(ppoint+0)=999999;
    numero=0;


    maxRnA=new double[taille];
    for (i=1; i<=taille; i++) *(maxRnA+i-1)=0;



    // Discrétisations de l'espace paramétrique (0,pi)^{d-2} x (0,2pi) de la sphère R^d: (O,pi) en N et (0,2*pi) 2*N points respectivement
    for (i = 1; i <= *(N+0); i++) {
      *(discretization1+i-1)=i*pi/(*(N+0)+1);
      *(discretization2+i-1)=i*pi/(*(N+0)+1);
    }
    for (i = *(N+0)+1; i <= 2**(N+0); i++) { 
      *(discretization2+i-1)=i*pi/(*(N+0)+1);
    }


    

    // On va calculer ||R_{n,A}|| pour chaque sous-ensemble A de {1,...,p} de |A|>1, et ceci par valeurs de |A| croissantes
    // Il faudra distinguer le cas sériel du cas non sériel. Dans le cas sériel, on ne prendra que les A qui contiennent 1.   
    // On fixe |A|
    for (cardA=2; cardA<=*(p+0); cardA++) { 
      

      *(cardApoint+0)=cardA;
      *(ppoint+0)=*(p+0);

      if (*(seriel+0) == 1) {
      *(cardApoint+0)=cardA-1;
      *(ppoint+0)=*(p+0)-1;
      }

      // On calcule Cppoint,cardApoint
      CpcardA=(int)Cnp(*(ppoint+0),*(cardApoint+0));

      // On initialise combmat
      combmat=new int[*(cardApoint+0)*CpcardA];
      for (i = 1; i<= *(cardApoint+0)*CpcardA; i++) *(combmat+i-1)=999999;

      //combmat va contenir toutes les permutations de 1:ppoint de taille cardApoint (matrice de taille cardApoint x Cppoint,cardApoint, c'est-à-dire tous les ensembles A de taille cardApoint
      combn(combmat, ppoint, cardApoint);

      // On réinitialise A
      A = new int[cardA];
      for (i = 1; i <= cardA; i++) *(A+i-1)=999999;	
      // On réinitialise permvecd et limit
      permvecd = new double[cardA];
      for (i = 1; i <= cardA; i++) *(permvecd+i-1)=999999;	
      limit = new double[cardA];
      for (i = 1; i <= cardA; i++) *(permvecd+i-1)=999999;	
	    
      // Pour chaque valeur de cardA, il y a Cp,cardA façons de choisir un ensemble de taille cardA parmi p valeurs
      // On va fixer A, pour |A| donné. i est donc l'indice de colonne de la matrice combmat
      for (i = 1; i<=CpcardA; i++) {	

	numero=numero+1;

	if (*(seriel+0) == 0) {
	// A est un sous ensemble (de taille cardA) qui contient la ième permutation de 1:p de taille cardA. C'est le i-ème ensemble A de taille |A|
	for (j=1; j<=cardA; j++) *(A+j-1)=*(combmat+(i-1)**(cardApoint+0)+j-1); // j est l'indice de ligne de la matrice combmat
	}
	if (*(seriel+0) == 1) {
	// A est un sous ensemble (de taille cardA) qui contient la ième permutation de 1:p de taille cardA qui contient 1. C'est le i-ème ensemble A de taille |A|
	  *(A+0)=1;
	  for (j=2; j<=cardA; j++) *(A+j-1)=*(combmat+(i-1)**(cardApoint+0)+j-2)+1; // j est l'indice de ligne de la matrice combmat
	}
	for (j = 1; j <= cardA; j++)  {
	  *(permvecd+j-1)=1;
	  *(limit+j-1)=2*(int)pow((double)*(N+0),*(vecd+*(A+j-1)-1)-1.0);
	}
	
	// A chaque passage dans la boucle while, *(permvecd+0), ..., *(permvecd+cardA-1) correspondra à une nouvelle permutation du vecteur (2*N^(vecd[A]-1))
	// Par exemple, si A={1,3,6} alors les permutations iront de (1,1,1) jusquà (2*N^{d_1-1},2*N^{d_3-1},2*N^{d_6-1})
	// Ainsi, si la permutation est (1,54,7) cela signifie que l'on prend le point 1 sur la sphère R^{d_1}, le point 54 sur la sphère R^{d_3} et le point 7 sur la sphère R^{d_6}
	while (1) {
	  

	  for (j = 1; j<= *(nbis+0)**(nbis+0); j++) *(PsiA+j-1)=1.0;
	  
	 
	  //Affichage d'un compteur à l'écran
	  //	  if (*(compt+0) != 0) {*(compt+0)=*(compt+0)-1; cout << "\n" << *(compt+0);}

	  
	  // Pour chaque k dans A (for k = 1; k<= cardA; k++), on calcule la permutation associée à x=*(permvecd+cardA-k) en utilisant l'équivalent C++ de 
	  // mon programme decomposition.R
	  // puis, pour cette permutation, on calcule le vecteur sk=s^{(k)} associé comme c'est fait dans le programme directions, c'est-à-dire qu'on pioche dans discretization1 et 2, 
	  // et on transforme en cordonnées cartésiennes. Ensuite, on calcule RnA.
	  // Par exemple, si discretization1=(alpha_1,...,alpha_N) et discretization2=(beta_1,...,beta_{2N})
	  // alors la ligne (1,2,1,...,3,4) (permutation associée à x) de longueur d_k-2+1=d_k-1 correspondra au choix des angles (alpha_1,alpha_2,alpha_1,...,alpha_3,beta_4)
	  for (k = 1; k<= cardA; k++) {

	    decomp = new int[*(vecd+*(A+k-1)-1)-1];
	    for ( jj = 1;  jj <= (*(vecd+*(A+k-1)-1)-1); jj++) *(decomp+jj-1)=999999;

	    //x est le numéro de la permutation pour s^{(k)}. La valeur maximale de x est 2*N^(vecd[A[k]]-1). Note: 2^32-2*50^(7-1) < 0.
	    // Par exemple, le point x=54 sur la sphère R^{d_3} correspondra (si d_3=4 et N=50) au choix (alpha=1,alpha=2,beta=4)
	    x = *(permvecd+k-1);


	    //Cas où d_i=1, il n'y a que deux points possibles sur la sphère dans R^1: matrice à 2 lignes et une colonne
	    if (*(vecd+*(A+k-1)-1)== 1) {
	      if ( (int)x == 1)   *(sk+0) = -1.0;
	      if ( (int)x == 2)   *(sk+0) = 1.0;
	    }
	
	    //Cas où d_i=2, c'est-à-dire sur un cercle: matrice à 2*N lignes et 2 colonnes
	    if (*(vecd+*(A+k-1)-1) == 2) {
	      *(sk+0) = cos(*(discretization2+(int)x-1));
	      *(sk+1) = sin(*(discretization2+(int)x-1));
	    }
	    
	    //Cas où d_i>2
	    if (*(vecd+*(A+k-1)-1) > 2) {
	      
	      //on calcule la permutation associée à x=*(permvecd+cardA-k) en utilisant l'équivalent C++ de mon programme decomposition.R, Elle doit être de longueur d_k-1
	      for (j = 1; j <= *(vecd+*(A+k-1)-1)-1; j++) {
		produit = pow((double)*(N+0),(double)*(vecd+*(A+k-1)-1)-1-j)*2**(N+0);
		quotient = (int)floor(x/produit);
		reste = fmod(x,produit);
		if (reste == 0) {
		  quotient = quotient-1;
		  reste = produit;
		}
		*(decomp+j-1)=quotient+1;
		x = reste;
	      }
	      // *(decomp+0), ..., *(decomp+*(vecd+*(A+k-1)-1)-2) correspond à la permutation (de longueur d_k-1) associée à x. Cela permettra de choisir  
	      // le x-ème vecteur s^{(k)} sur la sphère R^{d_k} en coordonnées sphériques
	      for (j = 1; j <= (*(vecd+*(A+k-1)-1)-2); j++) {
		*(decomp+j-1) = *(decomp+j);
	      }
	      *(decomp+*(vecd+*(A+k-1)-1)-2) = (int)reste;
	  
	      // On passe maintenant en coordonnées cartésiennes (s^{(k)}=(x_1,...,x_{d_A[k]}). Voir le livre de Martin p.32
	      *(sk+0)=cos(*(discretization1+*(decomp+0)-1)); // c'est x_1 (le x_n du livre de Martin)
	      
	      for (l = 2; l <= (*(vecd+*(A+k-1)-1)-1); l++) {
		prod=1;
		for (j=1; j<=l-1;j++) {
		  prod=prod*sin(*(discretization1+*(decomp+j-1)-1));     
		}
		if ((*(vecd+*(A+k-1)-1)-1) == 2) {*(sk+l-1) = prod*(cos(*(discretization2+*(decomp+l-1)-1))); }
		else {*(sk+l-1) = prod*(cos(*(discretization1+*(decomp+l-1)-1))); }      
	      }

	      // On calcule x_{d_A[k]} (le x_1 du livre de Martin qui est un produit de sinus)
	      prod=1;
	      for (j=1; j<=(*(vecd+*(A+k-1)-1)-2);j++) {
		prod=prod*sin(*(discretization1+*(decomp+j-1)-1));      
	      }
	      *(sk+*(vecd+*(A+k-1)-1)-1) = prod*sin(*(discretization2+*(decomp+*(vecd+*(A+k-1)-1)-2)-1)); // c'est x_{d_A[k]} (le x_1 du livre de Martin)
	      
	      
	    } // fin du Cas où d_i>2
	    

	    // Création des indices de début et de fin de colonnes pour extraire les données dans X 
	debut=0;
	if (*(A+k-1)==1) {debut=1;
	fin=*(vecd+0);
	}
	else {
	  fin=0;
	  for (j=1; j <= (*(A+k-1)-1); j++) {
	    debut=debut+*(vecd+j-1);
	    fin=fin+*(vecd+j-1);
	  }
	  debut=debut+1;
	  fin=fin+*(vecd+*(A+k-1)-1);
	}
	
	
	//indice de ligne	
	for (ii=1; ii<=*(nbis+0); ii++) {
	  
	  //On calcule t_i^{(k)} comme indiqué par Martin "For s^{(k)}, evaluate the n values <X_i^k,s^k>=t_i^k"
	  tik=0;
	  for (j=1; j<=(fin-debut+1); j++) tik=tik+*(newX+(debut+j-2)**(nbis+0)+ii-1)**(sk+j-1);
	  
	  
	  moyennecoli=0;
	  
	  for (l=1; l<=*(nbis+0); l++) {
	    
	    // The n x n matrix for s^k has an element in position (i,j) given by Ind{<X_i^k,s^k> \leq t_j^k}
	    somme=0.0;
	    ind=0;
	    for (j=1; j<=(fin-debut+1) ; j++) {
	      somme=somme+*(sk+j-1)**(newX+(debut+j-2)**(nbis+0)+l-1);
	    }
	    somme=somme-tik;
	    if (somme <= 0) ind=1;
	    
	    
	    // On remplit la matrice n*n comme indiqué par Martin, ligne l, colonne ii
	    // The n x n matrix for s^k has an element in position (i,j) given by Ind{<X_i^k,s^k> \leq t_j^k}
	    *(matrice+(ii-1)**(nbis+0)+l-1)=ind;
	    
	    moyennecoli=moyennecoli+ind;
	    
	    
	  }
	  
	  
	  
	  //On multiplie terme à terme (au fur et à mesure) les |A| matrices n*n obtenues. Il ne s'agit pas d'un produit matriciel mais d'un produit terme à terme.
	  //For a given subset A, one multiplies together the appropriate |A| such matrices to obtain an n x n matrix Psi (say)
	  for (l=1; l<=*(nbis+0); l++) *(PsiA+(ii-1)**(nbis+0)+l-1)=*(PsiA+(ii-1)**(nbis+0)+l-1)*(*(matrice+(ii-1)**(nbis+0)+l-1)-moyennecoli/(*(nbis+0)));
	  
	} //fin de for ii = 1 to n
	
	
	delete[] decomp;
	
      } // fin de for (k=1; k<=cardA; k++)
	  //A ce stade, la matrice PsiA est créée
      

	  // A vector psiA is then obtained by adding the rows of PsiA (i.e. psiA<-apply(PsiA,FUN=sum,MARGIN=2)). 
	  // The maximum value of the process R_{n,A} (for the given choice of directions) is then
	  // the max of the ABSOLUTE VALUES OF THE components of psiA divided by \sqrt{n} (i.e. maxRnA<-max(abs(psiA))/(sqrt(n)))
	  maxpsiA=0;
	  for (ii=1; ii<=*(nbis+0); ii++) {
	    sommePsiA=0;
	    for (l=1; l<=*(nbis+0); l++) sommePsiA=sommePsiA+*(PsiA+(ii-1)**(nbis+0)+l-1);
	    if (fabs(sommePsiA) > maxpsiA) maxpsiA=fabs(sommePsiA);
	  }
	  max=maxpsiA/sqrt((double)*(nbis+0));
	  
	  // The global max is then derived by varying the choice of directions
	  if (max > *(maxRnA+numero-1)) *(maxRnA+numero-1) = max;



	  for (ii = 1;  ii <= cardA;  ++ii) {
	    *(permvecd+cardA-ii) += 1;
	    if (*(permvecd+cardA-ii) <= *(limit+cardA-ii))
	      break;
	    *(permvecd+cardA-ii) = 1;
	  }
	  if (ii == cardA+1)
	    break;  // finish the cycle
	} // fin du while(1)
	
	
      } // fin de for i = 1 to CpcardA
      
	//On libère de la mémoire
      delete[] limit;
      delete[] permvecd;
      delete[] A;
      delete[] combmat;
      
    } // fin de for cardA = 2 to p
    


    //On modifie le pointeur d'entrée pour renvoyer les résultats
    for (g=1; g<=taille; g++)   *(res+g-1)=*(maxRnA+g-1);

    
    delete[] discretization1;
    delete[] discretization2;
    delete[] PsiA;
    delete[] matrice;
    delete[] sk;
    delete[] cardApoint;
    delete[] ppoint;
    delete[] maxRnA;
    delete[] newX;
    delete[] nbis;
    
  } // fin de la fonction normeRnAwR
  
  
  
} // extern C  


