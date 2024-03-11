/* Debut des commentaires

Nom de la fonction: dependogram
-------------------

Auteur: Pierre Lafaye de Micheaux
-------

Date: 21/11/2005
-----

Entrées:
--------

    int *N: nombre de points de la discrétisation de (0,PI)
    int *vecd: (d_1,d_2,...,d_p) où les d_j sont des entiers
    int *p: nombre de vecteurs aléatoires dont on veut tester l'indépendance
    double *X: matrice des données: n lignes et sum(vecd) colonnes
    int *n: nombre d'individus
    int *B: nombre de boucles bootstrap
    double *alpha: niveau du test
    double *RnAs: contiendra les résultats, vecteur des normes des RnA pour les différents A
    double *RnAsetoile:  contiendra les résultats, matrice des différents RnA pour tous les Xetoiles créées
    double *Rn:  contiendra les résultats, valeur de Rn

Sorties:
--------

    Les pointeurs modifiés sont *RnAs, *RnAsetoile et *Rn. La fonction ne renvoie rien.

Fonctions extérieures appelées:
-------------------------------

    normeRnAwR

Description:
------------

Cette fonction calcule une approximation (due au fait des discrétisations données) de ||R_{n,A}|| pour tous les 2^p-p-1 choix de A (cardA>1) possibles.
Elle calcule aussi  une approximation (due au fait des discrétisations données) de ||R_{n}||=\max_A ||R_{n,A}||.
Elle calcule enfin une approximation (due au fait des discrétisations données) de ||R_{n,A}^*|| pour tous les 2^p-p-1 choix de A (cardA>1) possibles, et pour 
tous les B échantillons (i.e. les matrices X^*) bootstrap générés.
L'appel depuis la fonction R trace aussi le dependogram.

Références:
-----------

Exemples:
---------

N<-10
vecd<-c(3,2)
X<-matrix(c(0.7,0.4,0.1,0.8,0.4,0.13,3.7,0.4,2.1,1.8,0.1,0.1,0.7,0.2,0.1),byrow=F,nrow=3)
x<-Sys.time();dependogram(X,vecd,N=10,B=200,alpha=0.05);Sys.time()-x

X<-matrix(c(runif(100),runif(100),rbinom(100,1,0.5)),nrow=100,byrow=F);x<-Sys.time();dependogram(X,vecd=c(1,1,1),N=10,B=7000,alpha=0.05);Sys.time()-x

Equivalent en R:
----------------

dependogram<-function(X,vecd.or.p,N=10,B=200,alpha=0.05,compt=1) {



if (length(vecd.or.p) > 1) {
#on fait le cas non sériel
vecd<-vecd.or.p
p<-length(vecd)
taille<-2^p-p-1


n<-nrow(X)
RnAs<-normeRnAwR(X,vecd.or.p,N,0,1) #vecteur de taille 'taille' contenant la norme de RnA pour chacun des 'taille' A différents avec |A|>1
Rn<-max(RnAs)
#Chaque ligne de la matrice RnAsetoile contiendra un vecteur de taille 'taille' contenant la norme de RnA pour chacun des 'taille' A différents (pour la matrice Xetoile considérée)
RnAsetoile<-matrix(0,nrow=B,ncol=taille)
Rnetoile<-rep(0,B)
debut<-c(1,cumsum(vecd)+1)[-length(vecd)-1]
fin<-cumsum(vecd)
for (b in 1:B) {
if (compt == 1) {cat(B-b);cat(" ")}
Xetoile<-c()
for (j in 1:p) {
aj<-sample(1:n,n)
Xetoile<-cbind(Xetoile,X[aj,(debut[j]:fin[j])])
}
RnAsetoile[b,]<-normeRnAwR(Xetoile,vecd,N,0,0)
Rnetoile[b]<-max(RnAsetoile[b,])
}
cat("\n")
#Ordonne les éléments de chaque colonne
matseuils<-apply(RnAsetoile,FUN=sort,MARGIN=2)


beta<-(1-alpha)^(1/taille)
#Contient les beta-quantiles des stats RnA pour chacun des des 2^p-p-1 ensembles A
seuilscomplets<-matseuils[round((1-beta)*B),]
#On trace une barre verticale pour chaque A de hauteur ||RnA||
plot(RnAs,type="h",ylim=c(0,max(max(seuilscomplets),max(RnAs))+0.1),xlim=c(0,2^p-p),main="Dependogram",xlab="Subsets",ylab="||RnA||")
#On met une étoile pour chaque beta-quantile de ||R_A||
points((1:(2^p-p-1)),seuilscomplets,pch="*")
res<-list(RnAs=RnAs,Rn=Rn,seuilscomplets=seuilscomplets)
return(res)


}


if (length(vecd.or.p)==1) {
#On fait le cas sériel
#Attention, on a considéré que n est grand par rapport à p et donc que l'on peut calculer R_{n,A} à la place de S_{n,A}
#mais sur la matrice X de taille nprime x nprime 
p<-vecd.or.p
n<-nrow(X)
q<-ncol(X)
vecd<-q
taille<-2^(p-1)-1


RnAs<-normeRnAwR(X,vecd.or.p,N,0,1) #vecteur de taille 'taille' contenant la norme de RnA pour chacun des 'taille' A différents avec |A|>1
Rn<-max(RnAs)
#Chaque ligne de la matrice RnAsetoile contiendra un vecteur de taille 'taille' contenant la norme de RnA pour chacun des 'taille' A différents (pour la matrice Xetoile considérée)
RnAsetoile<-matrix(0,nrow=B,ncol=taille)
Rnetoile<-rep(0,B)
for (b in 1:B) {
if (compt == 1) {cat(B-b);cat(" ")}
Xetoile<-c()
aj<-sample(1:n,n)
Xetoile<-X[aj,]
RnAsetoile[b,]<-normeRnAwR(Xetoile,vecd.or.p,N,0,0)
Rnetoile[b]<-max(RnAsetoile[b,])
}
cat("\n")

beta<-(1-alpha)^(1/taille)


#Le beta-quantile de S_n,A est calculé en amalgamant toutes les
#valeurs  S_n,A^* avec |A|=k comme dans l'article
#Il y a choose(p-1,|A|-1) ensembles A de taille |A| qui contiennent 1
#Il faut donc prendre dans la matrice RnAsetoile des paquets de choose(p-1,|A|-1) colonnes, pour |A|=2 to p.
#Il y aura p-1 paquets.
#Pour chacun de ces paquets, on crée un vecteur vecA en prenant tous les éléments du paquet, 
#puis on calcule seuilscomplets[|A|]<-vecA[round((1-beta)*B*choose(p-1,|A|-1))] pour |A|=2 to p
#Ce vecteur seuilscomplets contiendra donc les hauteurs des barres horizontales à placer sur le dependogram (une barre horizontale 
#pour chaque |A|)

seuilscomplets<-rep(0,p-1)
begin<-1
end<-0
for (cardA in 2:p) {
end<-end+choose(p-1,cardA-1)
vecA<-as.vector(RnAsetoile[,begin:end])
seuilscomplets[cardA-1]<-vecA[round((1-beta)*B*choose(p-1,cardA-1))]
begin<-end+1
}



#On trace une barre verticale pour chaque A de hauteur ||RnA||
plot(RnAs,type="h",ylim=c(0,max(c(max(seuilscomplets),max(RnAs),max(Rnetoile)))+0.1),xlim=c(0,2^(p-1)),main="Dependogram",xlab="Subsets",ylab="||RnA||")

Rnetoile<-sort(Rnetoile)
Leseuil<-Rnetoile[round((1-alpha)*B)]
abline(h=Leseuil,col="red")


##Il reste à placer les lignes horizontales des seuils critiques pour chaque |A|


begin<-1
end<-0
for (cardA in 2:p) {
end<-end+choose(p-1,cardA-1)

segments(begin-0.5,seuilscomplets[cardA-1],end+0.5,seuilscomplets[cardA-1],lty=4)

begin<-end+1
}






res<-list(RnAs=RnAs,Rn=Rn,seuilscomplets=seuilscomplets)
return(res)
}



}




Instructions de compilation pour utilisation dans un terminal:
--------------------------------------------------------------

g++ -Wall dependogram.cpp 
./a.out

g++ -Wall -fPIC  -O2 -march=i686 -fomit-frame-pointer dependogram.cpp 
./a.out

Instructions de compilation pour utilisation depuis R:
------------------------------------------------------

g++ -I/usr/lib/R/include  -I/usr/local/include  -fPIC  -O2 -march=i686 -fomit-frame-pointer -c dependogram.cpp -o dependogram.o
g++ -shared -L/usr/local/lib -o dependogram.so dependogram.o 

Pour utiliser dans R, taper source("dependogram.R") où le fichier dependogram.R contient le code R suivant:
 

dependogram <- function(X,vecd.or.p,N=10,B=200,alpha=0.05,affiche=1) {

X<-as.matrix(X)

#si length(vecd.or.p)>1 alors cas non sériel sinon cas sériel

if (length(vecd.or.p) > 1) {
#on fait le cas non sériel
seriel<-0
vecd<-vecd.or.p
p<-length(vecd)
taille<-2^p-p-1




		RnAs<-rep(0,2^p-p-1)
		RnAsetoile<-matrix(0,nrow=(2^p-p-1),ncol=B)
		Rn<-0
		Rnetoile<-rep(0,B)

#On charge la fonction C dans la mémoire
	     dyn.load(paste("dependogram", .Platform$dynlib.ext, sep=""))

#Remarque: quand on passe une matrice dans la fonction .C elle est reçue comme un vecteur obtenu en concaténant les colonnes de cette matrice
#et le résultat est lui aussi renvoyé sous la forme d'un tel vecteur

#On appelle la fonction C dependogram
		out <- .C("dependogram",
		as.integer(N),
		as.integer(vecd),
		as.integer(length(vecd)),
		as.integer(p),
		as.numeric(X),
		as.integer(nrow(X)),
		as.integer(ncol(X)),
		as.integer(B),
		as.numeric(alpha),
		RnAs=as.numeric(RnAs),
		RnAsetoile=as.numeric(RnAsetoile),
		Rn=as.numeric(Rn),
		Rnetoile=as.numeric(Rnetoile),
		as.integer(seriel))


#On décharge la fonction C de la mémoire
		dyn.unload(paste("dependogram", .Platform$dynlib.ext, sep=""))



if (affiche != 0) {
require(combinat)
RES<-as.list(1:(p-1))
for (cardA in 2:p) {RES[[cardA]]<-as.matrix(combn(p,cardA))}
nb<-0
for (cardA in 2:p) {
for (j in 1:(choose(p,cardA))) {
nb<-nb+1

cat(c("A:",RES[[cardA]][,j],"||RnA||:",round(out$RnAs[nb],3),"\n"))


}
}
}


#Ordonne les éléments de chaque ligne
RnAsetoile<-matrix(out$RnAsetoile,nrow=(2^p-p-1),ncol=B,byrow=F)
matseuils<-t(apply(RnAsetoile,FUN=sort,MARGIN=1))


beta<-(1-alpha)^(1/taille)
#Contient les beta-quantiles des stats RnA pour chacun des des 2^p-p-1 ensembles A
seuilscomplets<-matseuils[,round(beta*B)]
#On trace une barre verticale pour chaque A de hauteur ||RnA||
plot(out$RnAs,type="h",ylim=c(0,max(c(max(seuilscomplets),max(out$RnAs),max(out$Rnetoile)))+0.1),xlim=c(0,2^p-p),main="Dependogram",xlab="Subsets",ylab="||RnA||")

Rnetoile<-sort(out$Rnetoile)
Leseuil<-Rnetoile[round((1-alpha)*B)]
abline(h=Leseuil,col="red")


#On met une étoile pour chaque beta-quantile de ||R_A||
points((1:(2^p-p-1)),seuilscomplets,pch="*")
res<-list(RnAs=out$RnAs,Rn=out$Rn,seuilscomplets=seuilscomplets,Leseuil=Leseuil)
return(res)


}



if (length(vecd.or.p) == 1) {
#on fait le cas non sériel
seriel<-1
p<-vecd.or.p
vecd<-rep(ncol(X),p)
taille<-2^(p-1)-1




		RnAs<-rep(0,taille)
		RnAsetoile<-matrix(0,nrow=taille,ncol=B)
		Rn<-0
		Rnetoile<-rep(0,B)

#On charge la fonction C dans la mémoire
	     dyn.load(paste("dependogram", .Platform$dynlib.ext, sep=""))

#Remarque: quand on passe une matrice dans la fonction .C elle est reçue comme un vecteur obtenu en concaténant les colonnes de cette matrice
#et le résultat est lui aussi renvoyé sous la forme d'un tel vecteur

#On appelle la fonction C dependogram
		out <- .C("dependogram",
		as.integer(N),
		as.integer(vecd),
		as.integer(length(vecd)),
		as.integer(p),
		as.numeric(X),
		as.integer(nrow(X)),
		as.integer(ncol(X)),
		as.integer(B),
		as.numeric(alpha),
		RnAs=as.numeric(RnAs),
		RnAsetoile=as.numeric(RnAsetoile),
		Rn=as.numeric(Rn),
		Rnetoile=as.numeric(Rnetoile),
		as.integer(seriel))


#On décharge la fonction C de la mémoire
		dyn.unload(paste("dependogram", .Platform$dynlib.ext, sep=""))


if (affiche != 0) {

require(combinat)
RES<-as.list(1:(p-1))
for (cardA in 2:p) {RES[[cardA]]<-as.matrix(rbind(rep(1,choose(p-1,cardA-1)),as.matrix(combn(p-1,cardA-1)+1)))}
nb<-0
for (cardA in 2:p) {
for (j in 1:(choose(p-1,cardA-1))) {
nb<-nb+1

cat(c("A:",RES[[cardA]][,j],"||RnA||:",round(out$RnAs[nb],3),"\n"))

}
}
}


RnAsetoile<-matrix(out$RnAsetoile,nrow=taille,ncol=B,byrow=F)
Rn<-max(out$RnAs)


beta<-(1-alpha)^(1/taille)


#Le beta-quantile de S_n,A est calculé en amalgamant toutes les
#valeurs  S_n,A^* avec |A|=k comme dans l'article
#Il y a choose(p-1,|A|-1) ensembles A de taille |A| qui contiennent 1
#Il faut donc prendre dans la matrice RnAsetoile des paquets de choose(p-1,|A|-1) lignes, pour |A|=2 to p.
#Il y aura p-1 paquets.
#Pour chacun de ces paquets, on crée un vecteur vecA en prenant tous les éléments du paquet, 
#puis on calcule seuilscomplets[|A|]<-vecA[round(beta*B*choose(p-1,|A|-1))] pour |A|=2 to p
#Ce vecteur seuilscomplets contiendra donc les hauteurs des barres horizontales à placer sur le dependogram (une barre horizontale 
#pour chaque |A|)

seuilscomplets<-rep(0,p-1)
begin<-1
end<-0
for (cardA in 2:p) {
end<-end+choose(p-1,cardA-1)
vecA<-as.vector(RnAsetoile[begin:end],)
vecA<-sort(vecA)
seuilscomplets[cardA-1]<-vecA[round(beta*B*choose(p-1,cardA-1))]
begin<-end+1
}



#On trace une barre verticale pour chaque A de hauteur ||RnA||
plot(out$RnAs,type="h",ylim=c(0,max(c(max(seuilscomplets),max(out$RnAs),max(out$Rnetoile)))+0.1),xlim=c(0,2^(p-1)),main="Dependogram",xlab="Subsets",ylab="||RnA||")

Rnetoile<-sort(out$Rnetoile)
Leseuil<-Rnetoile[round((1-alpha)*B)]
abline(h=Leseuil,col="red")


##Il reste à placer les lignes horizontales des seuils critiques pour chaque |A|


begin<-1
end<-0
for (cardA in 2:p) {
end<-end+choose(p-1,cardA-1)

segments(begin-0.5,seuilscomplets[cardA-1],end+0.5,seuilscomplets[cardA-1],lty=4)

begin<-end+1
}






res<-list(RnAs=out$RnAs,Rn=out$Rn,seuilscomplets=seuilscomplets,Leseuil=Leseuil)
return(res)
}



}






Utilisation du debugger gdb avec electric fence:
------------------------------------------------

g++ -g dependogram.cpp -lm -u malloc -lefence
gdb ./a.out
run

Utilisation du debugger ddd avec R:
-----------------------------------

g++ -I/usr/lib/R/include  -I/usr/local/include  -fPIC -g -c dependogram.cpp -o dependogram.o 
g++ -shared -L/usr/local/lib -o dependogram.so dependogram.o 
R -d ddd

Menu Program: cocher Run in Execution Window
Menu Program: Run puis Run

dyn.load(paste("dependogram", .Platform$dynlib.ext, sep=""))

Menu File/Open source...
Cliquer sur Load Shared Object Library Symbols
Sélectionner dependogram.cpp
Cliquer sur Open
Mettre des breakpoints
Dans la fenêtre Execution Window de R, taper: 

source("dependogram.R")
N<-10
vecd<-c(3,2)
X<-matrix(c(0.7,0.4,0.1,0.8,0.4,0.13,3.7,0.4,2.1,1.8,0.1,0.1,0.7,0.2,0.1),byrow=F,nrow=3)
x<-Sys.time();dependogram(X,vecd,N=10,B=200,alpha=0.05);Sys.time()-x

Fin des commentaires */



// Inclusion de librairies et de fonctions extérieures
//----------------------------------------------------

#include <iostream>
using namespace std;
#include <math.h>
#include <ctime>    // For time()
//#include <cstdlib>  // For srand() and rand()

#include "others/normeRnAwR.cpp"
//#include "others/ran1.cpp"

// Pour debuger avec Electric Fence
//#include "efencepp.h"

// Utilisation dans une fonction main:
// -----------------------------------


/*

extern "C" {


   int main()

   {

     void dependogram(int *N, int *vecd, int *p, double *X, int *n, int *q, int *B, double *alpha, double *RnAs, double *RnAsetoile, double *Rn);
 
    //Déclaration des variables
     int i, j;
    int *N;
    int *vecd, *p;
    int *n, nbcol;
    double *X;
    int *B;
    double *alpha;
    double *Rn, *RnAs, *RnAsetoile;  

    //Initialisation des variables
    N=new int[1];
    *(N+0)=10;     //N: nombre de points de la discrétisation
    p=new int[1];
    *(p+0)=2;     //longueur de vecd
    vecd = new int[*(p+0)];
    *(vecd+0)=3;     //vecd=(d_1,d_2,...,d_p) où les d_j sont des entiers
    *(vecd+1)=2;   
    n=new int[1];
    *(n+0)=3;
    nbcol=5;
    X = new double[*(n+0)*nbcol];     //X est la matrice des données: n lignes et sum(vecd)=nbcol colonnes, remplie par colonnes
    *(X+0)=0.7;*(X+1)=0.4;*(X+2)=0.1;*(X+3)=0.8;*(X+4)=0.4;*(X+5)=0.13;
    *(X+6)=3.7;*(X+7)=0.4;*(X+8)=2.1;*(X+9)=1.8;
    *(X+10)=0.1;*(X+11)=0.1;*(X+12)=0.7;*(X+13)=0.2;*(X+14)=0.1;    
    B = new int[1];
    *(B+0)=7;
    alpha = new double[1];
    *(alpha+0)=0.05;
    Rn = new double[1];
    *(Rn+0)=999999;
    RnAs = new double[(int)pow(2.0,*(p+0))-*(p+0)-1];
    for (i = 1; i <= ((int)pow(2.0,*(p+0))-*(p+0)-1); i++) *(RnAs+i-1)=999999;
    RnAsetoile = new double[*(B+0)*((int)pow(2.0,*(p+0))-*(p+0)-1)];
    for (i = 1; i <= *(B+0)*((int)pow(2.0,*(p+0))-*(p+0)-1); i++) *(RnAsetoile+i-1)=999999;


    //Appel de la fonction
    dependogram(N, vecd, p, X, n, B, alpha, RnAs, RnAsetoile, Rn);


    //Affichage des résultats
    cout << "RnAs:\n";
    for (i = 1; i <= (int)pow(2.0,*(p+0))-*(p+0)-1; i++) {
      cout << *(RnAs+i-1) << " ";
    }
    cout << "\n";
    cout << "Rn:\n";
    cout << *(Rn+0);
    cout << "\n";
    cout << "RnAsetoile:\n";
    for (i = 1; i <= *(B+0); i++) {
      for (j = 1; j <= (int)pow(2.0,*(p+0))-*(p+0)-1; j++) {

	cout << *(RnAsetoile+(i-1)*((int)pow(2.0,*(p+0))-*(p+0)-1)+j-1) << " ";
      }
    cout << "\n";

}

    delete[] N;
    delete[] p;
    delete[] vecd;
    delete[] n;
    delete[] X;
    delete[] B;
    delete[] alpha;
    delete[] Rn;
    delete[] RnAs;
    delete[] RnAsetoile;


    return 0;

   }

} // extern C

*/

//-----------------------------------------------------------------
// DEBUT DE LA FONCTION
//-----------------------------------------------------------------

#include <R.h>
#include "Rmath.h"
#include <R_ext/Rdynload.h>

extern "C" {


  void dependogramC(int *N, int *vecd, int *lenvecd, int *p, double *X, int *n, int *q, int *B, double *alpha, double *RnAs, double *RnAsetoile, double *Rn, double *Rnetoile, int *seriel) 

  {
 
    double runif(double a, double b);  
    //    int randint(int debut, int fin, long *idum);

    GetRNGstate();

    //Déclaration des variables
    double *compt;
    int i, b, j, l;
    int *debut, *fin;
    double *Xetoile;
    int *aj;
    int lg,col;
    double *resul;
    double max;
    long *idum; //seed pour le générateur

   //Initialisation des variables
    int taille;
    if (*(seriel+0) == 0) {
      taille=(int)pow(2.0,(double)*(p+0))-*(p+0)-1;
    }
    else {
      if (*(seriel+0) != 1) {return;}
      taille=(int)pow(2.0,(double)*(p+0)-1.0)-1;
    }
    compt=new double[1];
    *(compt+0)=0;
    debut = new int[*(p+0)];
    for (i = 1; i <= *(p+0); i++) *(debut+i-1)=999999;
    fin = new int[*(p+0)];
    for (i = 1; i <= *(p+0); i++) *(fin+i-1)=999999;
    resul = new double[taille];
    for (i = 1; i <= taille; i++) *(resul+i-1)=999999;
    Xetoile = new double[*(n+0)**(q+0)];
    for (i = 1; i <= *(n+0)**(q+0); i++) *(Xetoile+i-1)=999999;
    aj = new int[*(n+0)];
    for (i = 1; i <= *(n+0); i++) *(aj+i-1)=999999;
    idum=new long[1];
    *(idum+0)=-1;
    //    *(idum+0)=-time(NULL);

    // RnAs: vecteur de taille 2^p-p-1 contenant la norme de RnA pour chacun des 2^p-p-1 A différents avec |A|>1
    normeRnAwR(N, vecd, lenvecd, p, X, n, q, RnAs, compt, seriel);

    //Rn=max(RnAs)
    *(Rn+0)=*(RnAs+0);
    for (i = 1; i <= taille; i++) {
      if (*(RnAs+i-1) > *(Rn+0)) {*(Rn+0)=*(RnAs+i-1);}
    }



    /*

#Chaque colonne de la matrice RnAsetoile contiendra un vecteur de taille 2^p-p-1 contenant la norme de RnA pour chacun des 2^p-p-1 A différents (pour la matrice Xetoile considérée)
RnAsetoile<-matrix(0,nrow=(2^p-p-1),ncol=B)
Rnetoile<-rep(0,B)

    */

    //fin=cumsum(vecd)
    *(fin+0)=*(vecd+0);
    for (i = 2; i <= *(p+0); i++) {*(fin+i-1)=*(fin+i-2)+*(vecd+i-1);}

    //debut=c(1,cumsum(vecd)+1)[-length(vecd)-1]
    *(debut+0)=1;
    for (i = 2; i <= *(p+0); i++) {*(debut+i-1)=*(fin+i-2)+1;}



    //On commence la boucle bootstrap
    for (b = 1; b<= *(B+0); b++) {

      if (*(seriel+0) == 0) {

      for (j = 1; j <= *(p+0); j++) {


	//aj=sample(1:n,n)
	for (l = 1; l<= *(n+0); l++) {
                                   	
	  *(aj+l-1) = 1 + (int)(runif(0.0,1.0)*n[0]);                  //randint(1,*(n+0),idum);

	} //aj est maintenant défini



	// Xetoile=cbind(Xetoile,X[aj,(debut[j]:fin[j])])
	// A VERIFIER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//	for (lg = 1; lg<= *(n+0); lg++) {
	// for (col = *(debut+j-1); col <= *(fin+j-1); col ++) {
	//   *(Xetoile+(lg-1)**(q+0)+col-1)=*(X+(*(aj+lg-1)-1)**(q+0)+col-1);
	// }}
	for (col = *(debut+j-1); col <= *(fin+j-1); col ++) {
	  for (lg = 1; lg<= *(n+0); lg++) {
	    *(Xetoile+(col-1)**(n+0)+lg-1)=*(X+(col-1)**(n+0)+*(aj+lg-1)-1);
	  }}


      } //fin de for j=1 to p

      }

      if (*(seriel+0) == 1) {

	//aj=sample(1:n,n)
	for (l = 1; l<= *(n+0); l++) {
                                   	
	  *(aj+l-1) = 1 + (int)(runif(0.0,1.0)*n[0]);                  //randint(1,*(n+0),idum);
	  //	*(aj+l-1)=randint(1,*(n+0),idum);

	} //aj est maintenant défini

	// Xetoile<-X[aj,]
	  for (col = 1; col <= *(q+0); col++) {
	    for (lg = 1; lg <= *(n+0); lg++) {
	      *(Xetoile+(col-1)**(n+0)+lg-1)=*(X+(col-1)**(n+0)+*(aj+lg-1)-1);
	  }}


      }

      //RnAsetoile[b,]<-normeRnAwR(Xetoile,vecd,N,0,0)

      //Appel de la fonction
      normeRnAwR(N, vecd, lenvecd, p, Xetoile, n, q, resul, compt, seriel);

      max=*(resul+0);
      for (col = 1; col <= taille; col++) {
	*(RnAsetoile+(b-1)*taille+col-1)=*(resul+col-1);
	if (*(resul+col-1) > max ) {max=*(resul+col-1);}
      }


      //Rnetoile[b]<-max(RnAsetoile[,b])
      *(Rnetoile+b-1)=max;



    } //fin de la boucle bootstrap


    //On libère de la mémoire
    delete[] aj;
    delete[] resul;
    delete[] Xetoile;
    delete[] compt;
    delete[] debut;
    delete[] fin;

    PutRNGstate();
    return;

  }

} //extern C

