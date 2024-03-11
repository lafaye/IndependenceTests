/* Begining comments

Function name: cubature1

Entries: b double smoothing factor, q integer dimension of space

Output: D DiagonalMatrix

Description:
Cubature formula:
E_n^{r^2}: 7-2 (n>2) Degree 7, 2^{n+1}+4n^2 Points.


Author: Pierre Lafaye de Micheaux

Date: 30/05/2014

End comments */

#include <R.h>
#include "Rmath.h"
#include <complex>
//#include <iostream>
using namespace std;

extern "C" {
#include "others/gauher.c"
#include "others/binarycode.c"


  void cubature1(double *a, int *q, int *N, double *prec, int *MAXIT, int * cuba, double *resweights, double *respoints) {

  //Declaration des variables
    double gammafn (double x);
    int i, j, k;
    double *weights;
    weights = new double[N[0]];
    double *points;
    points = new double[N[0]*q[0]];
   
    // double PI;

    //    PI = 3.1415926535897931160;

  //----------------------------------------------------------------------------------------
  //CUBATURE
  //----------------------------------------------------------------------------------------


    if (q[0] == 1) {
      //Declaration des fonctions
      void gauher(double a, int N, double prec, int MAXIT, double *weightsandpoints);

      double *weightsandpoints;
      weightsandpoints = new double[N[0]*2];
      gauher(a[0],N[0],prec[0],MAXIT[0],weightsandpoints);
      for (i=1; i<=N[0]; i++) {
	points[i-1] = weightsandpoints[i - 1];
	weights[i-1] = weightsandpoints[N[0] + i - 1];
      }
      delete[] weightsandpoints;
    }
  //End of if q == 1



  if (q[0] == 2) {
	
    //cuba[0]: Choice of cubature of the plane
    
    if (cuba[0] == 1) {
      //E_2^(r^2) : 15-1. Degree 15, 44 Points in Stroud (1971) p.326
      N[0] = 44;
      
      double s1,s2,s3,s4,s5,s6,s7,s8,s9;
      double r1,r2,r3,r4,r5,r6,r7,r8,r9;
      double B1,B2,B3,B4,B5,B6,B7,B8,B9;

      r1=3.53838872812180699759844816420977;
      s1=0.0;
      B1=8.00648356965962968878302657382053 * pow((double)10,(double)(-6));

      r2=2.35967641687792858456109163436214;
      s2=0.0;
      B2=3.60457742083826400410253786082198 * pow((double)10,(double)(-3));

      r3=1.31280184462092663946003191058349;
      s3=0.0;
      B3=0.118760933075913674443102663431639;

      r4=0.538955948211420514863252356676544;
      s4=0.0;
      B4=0.437248854379140237467918223689153; 

      r5=2.30027994980565789465239975311549;
      s5=2.30027994980565789465239975311549;
      B5=3.67173507583298936096131754975575 * pow((double)10,(double)(-5));

      r6=1.58113883008418966599944677221635;
      s6=1.58113883008418966599944677221635;
      B6=5.65486677646162782923275808990310 * pow((double)10,(double)(-3)); 

      r7=0.841850433581927898923665650469697;
      s7=0.841850433581927898923665650469697;
      B7=0.177777426842423967403376002318122;

      r8=2.68553358175534090063094214167163;
      s8=1.11238443177145724971821342047472;
      B8=2.73544964785329001953807301017241 * pow((double)10,(double)(-4)); 

      r9=1.74084751439740260707930075663572;
      s9=0.721082650486896005766801022499833;
      B9=0.0208798455693859454703613248130647; 

      points[0] = r1;
      points[1] = s1;
      weights[0] = B1;

      points[q[0]] = r2;
      points[q[0] + 1] = s2;
      weights[1] = B2;

      points[2 * q[0]] =r3;
      points[2 * q[0] + 1] =s3;
      weights[2] = B3;
	
      points[3 * q[0]] = r4;
      points[3 * q[0] + 1] = s4;
      weights[3] = B4;

      points[4 * q[0]] = r5;
      points[4 * q[0] + 1] = s5;
      weights[4] = B5;

      points[5 * q[0]] = r6;
      points[5 * q[0] + 1] = s6;
      weights[5] = B6;

      points[6 * q[0]] = r7;
      points[6 * q[0] + 1] = s7;
      weights[6] = B7;

      points[7 * q[0]] = r8;
      points[7 * q[0] + 1] = s8;
      weights[7] = B8;

      points[8 * q[0]] = r9;
      points[8 * q[0] + 1] = s9;
      weights[8] = B9;

      points[9 * q[0]] = r5;
      points[9 * q[0] + 1] = -s5;
      weights[9] = B5;

      points[10 * q[0]] = r6;
      points[10 * q[0] + 1] = -s6;
      weights[10] = B6;

      points[11 * q[0]] = r7;
      points[11 * q[0] + 1] = -s7;
      weights[11] = B7;

      points[12 * q[0]] = r8;
      points[12 * q[0] + 1] = -s8;
      weights[12] = B8;

      points[13 * q[0]] = r9;
      points[13 * q[0] + 1] = -s9;
      weights[13] = B9;

      points[14 * q[0]] = -r1;
      points[14 * q[0] + 1] = s1;
      weights[14] = B1;

      points[15 * q[0]] = -r2;
      points[15 * q[0] + 1] = s2;
      weights[15] = B2;

      points[16 * q[0]] = -r3;
      points[16 * q[0] + 1] = s3;
      weights[16] = B3;

      points[17 * q[0]] = -r4;
      points[17 * q[0] + 1] = s4;
      weights[17] = B4;

      points[18 * q[0]] = -r5;
      points[18 * q[0] + 1] = s5;
      weights[18] = B5;

      points[19 * q[0]] = -r6;
      points[19 * q[0] + 1] = s6;
      weights[19] = B6;

      points[20 * q[0]] = -r7;
      points[20 * q[0] + 1] = s7;
      weights[20] = B7;

      points[21 * q[0]] = -r8;
      points[21 * q[0] + 1] = s8;
      weights[21] = B8;

      points[22 * q[0]] = -r9;
      points[22 * q[0] + 1] = s9;
      weights[22] = B9;

      points[23 * q[0]] = -r5;
      points[23 * q[0] + 1] = -s5;
      weights[23] = B5;

      points[24 * q[0]] = -r6;
      points[24 * q[0] + 1] = -s6;
      weights[24] = B6;

      points[25 * q[0]] = -r7;
      points[25 * q[0] + 1] = -s7;
      weights[25] = B7;

      points[26 * q[0]] = -r8;
      points[26 * q[0] + 1] = -s8;
      weights[26] = B8;

      points[27 * q[0]] = -r9;
      points[27 * q[0] + 1] = -s9;
      weights[27] = B9;

      points[28 * q[0]] = s1;
      points[28 * q[0] + 1] = r1;
      weights[28] = B1;

      points[29 * q[0]] = s2;
      points[29 * q[0] + 1] = r2;
      weights[29] = B2;

      points[30 * q[0]] = s3;
      points[30 * q[0] + 1] = r3;
      weights[30] = B3;

      points[31 * q[0]] = s4;
      points[31 * q[0] + 1] = r4;
      weights[31] = B4;

      points[32 * q[0]] = s8;
      points[32 * q[0] + 1] = r8;
      weights[32] = B8;

      points[33 * q[0]] = s9;
      points[33 * q[0] + 1] = r9;
      weights[33] = B9;

      points[34 * q[0]] = s1;
      points[34 * q[0] + 1] = -r1;
      weights[34] = B1;

      points[35 * q[0]] = s2;
      points[35 * q[0] + 1] = -r2;
      weights[35] = B2;

      points[36 * q[0]] = s3;
      points[36 * q[0] + 1] = -r3;
      weights[36] = B3;

      points[37 * q[0] + (1- 1)] = s4;
      points[37 * q[0] + (2- 1)] = -r4;
      weights[37] = B4;

      points[38 * q[0]] = s8;
      points[38 * q[0] + 1] = -r8;
      weights[38] = B8;

      points[39 * q[0]] = s9;
      points[39 * q[0] + 1] = -r9;
      weights[39] = B9;

      points[40 * q[0]] = -s8;
      points[40 * q[0] + 1] = -r8;
      weights[40] = B8;

      points[41 * q[0]] = -s9;
      points[41 * q[0] + 1] = -r9;
      weights[41] = B9;

      points[42 * q[0]] = -s8;
      points[42 * q[0] + 1] = r8;
      weights[42] = B8;

      points[43 * q[0]] = -s9;
      points[43 * q[0] + 1] = r9;
      weights[43] = B9;
        
    }
    
    if (cuba[0]==2) {
      
      //Degree 15, 44 Points
      //Hae76 A. Haegemans, Circularly symmetrical integration formulas 
      //for two-dimensional circularly symmetrical regions, BIT 16 (1976), 52--59. 
      N[0] = 44;

      static double M[9][4] = {
	{1 , 8 , 1.88427977148214959597160789463215 * pow((double)10,(double)0) , 2.08798455693859454703613248130647 * pow((double)10,(double)(-2))} , 
	{2 , 8 , 2.90680060251527712149754356676149 * pow((double)10,(double)0) , 2.73544964785329001953807301017241 * pow((double)10,(double)(-4))} ,
	{3 , 4 , 1.19055630066123290136498536710293 * pow((double)10,(double)0) , 1.77777426842423967403376002318122 * pow((double)10,(double)(-1))} ,
	{4 , 4 , 2.23606797749978969640917366873127 * pow((double)10,(double)0) , 5.65486677646162782923275808990310 * pow((double)10,(double)(-3))} ,
	{5 , 4 , 3.25308710227006371908007902786315 * pow((double)10,(double)0) , 3.67173507583298936096131754975575 * pow((double)10,(double)(-5))} ,
	{6 , 4 , 5.38955948211420514863252356676544 * pow((double)10,(double)(-1)) , 4.37248854379140237467918223689153 * pow((double)10,(double)(-1))} ,
	{7 , 4 , 1.31280184462092663946003191058349 * pow((double)10,(double)0) , 1.18760933075913674443102663431639 * pow((double)10,(double)(-1))} ,
	{8 , 4 , 2.35967641687792858456109163436214 * pow((double)10,(double)0) , 3.60457742083826400410253786082198 * pow((double)10,(double)(-3))} ,
	{9 , 4 , 3.53838872812180699759844816420977 * pow((double)10,(double)0) , 8.00648356965962968878302657382053 * pow((double)10,(double)(-6))}
      };
      
      int somme, somme2;
      somme = 0;
      somme2 = 0;
      
      for (i=1; i<=9; i=i+1) {
	for (j=1; j<=((int)M[i-1][1]); j=j+1) {
	  weights[somme+j-1] = M[i-1][3];
	}
	somme2 = 0;
	for (j=1; j<=((int)M[i-1][1]/4); j=j+1) {
	  for (k=1; k<=4; k=k+1) {
	    points[(somme+somme2+k-1) * q[0]] = M[i-1][2] * cos((k-1)*PI/2.0);
	    points[(somme+somme2+k-1) * q[0] + 1] = M[i-1][2] * sin((k-1)*PI/2.0);
	  }
	  somme2 = somme2+((int)M[i-1][1])/4;
	}	
	somme = somme + (int)M[i-1][1];
      }
      
    }
    
    
    if (cuba[0]==3) {

      // Degree 31, 172 Points
      //Hae75 A. Haegemans, Tables of circularly symmetrical integration formulas of degree 2d-1 
      //for two-dimensional circularly symmetrical regions, Report TW 27, K.U. Leuven Applied Mathematics
      //and Programming Division, 1975. 
      N[0] = 172;

      static double M[25][4] = {
	{1 , 16 , 2.19342435218571098000852165314022 * pow((double)10,(double)0) , 2.73236456089821369938846899455853 * pow((double)10,(double)(-3))} ,
	{2 , 16 , 2.94143072785882658443890464121376 * pow((double)10,(double)0) , 7.45882708744579780334967451875029 * pow((double)10,(double)(-5))} ,
	{3 , 16 , 3.70011740521483590909984391822818 * pow((double)10,(double)0) , 6.51579202758359447548582435330121 * pow((double)10,(double)(-7))} ,
	{4 , 16 , 4.56574266380320843274636697274591 * pow((double)10,(double)0) , 7.75679616707262954339902390779625 * pow((double)10,(double)(-10))} ,
	{5 , 8 , 1.21580764835581239743824121978989 * pow((double)10,(double)0) , 6.29510894883832189721267177298424 * pow((double)10,(double)(-2))} ,
	{6 , 8 , 1.75022957356728645712052732323077 * pow((double)10,(double)0) , 1.58692609223783075894586018940927 * pow((double)10,(double)(-2))} ,
	{7 , 8 , 2.56265642157773665089210281907111 * pow((double)10,(double)0) , 5.27660220867058647807407843938869 * pow((double)10,(double)(-4))} ,
	{8 , 8 , 3.30618113478644619137730123417901 * pow((double)10,(double)0) , 8.77872589293851750256952094712742 * pow((double)10,(double)(-6))} ,
	{9 , 8 , 4.10268615562676937418264214216800 * pow((double)10,(double)0) , 3.38247075123041995374022839198793 * pow((double)10,(double)(-8))} ,
	{10 , 8 , 5.13066763122356142415670471249892 * pow((double)10,(double)0) , 5.26685826507637745552635596537658 * pow((double)10,(double)(-12))} ,
	{11 , 4 , 7.81924734415062654910759895512688 * pow((double)10,(double)(-1)) , 1.67233060864508092720103720312852 * pow((double)10,(double)(-1))} ,
	{12 , 4 , 1.39054046557164805349275297468774 * pow((double)10,(double)0) , 3.85569294188414023000355252036718 * pow((double)10,(double)(-2))} ,
	{13 , 4 , 1.81698526640070493276888680298209 * pow((double)10,(double)0) , 9.93275640431363439042042900758230 * pow((double)10,(double)(-3))} ,
	{14 , 4 , 2.58066849114911903633979005815686 * pow((double)10,(double)0) , 4.80414432703307303553096755498915 * pow((double)10,(double)(-4))} ,
	{15 , 4 , 3.32321802494503855192283519847459 * pow((double)10,(double)0) , 7.90743971209622142932287551395032 * pow((double)10,(double)(-6))} ,
	{16 , 4 , 4.12589037924347407669696623400778 * pow((double)10,(double)0) , 2.84841137511625856222914116013963 * pow((double)10,(double)(-8))} ,
	{17 , 4 , 5.20116644740776271749685173656130 * pow((double)10,(double)0) , 2.97220998361754914429195687074397 * pow((double)10,(double)(-12))} ,
	{18 , 4 , 3.57047855200640258244677801924088 * pow((double)10,(double)(-1)) , 2.24741451714247573527439486994526 * pow((double)10,(double)(-1))} ,
	{19 , 4 , 8.57582269837600053440432674664275 * pow((double)10,(double)(-1)) , 1.32619067539941897572390100531881 * pow((double)10,(double)(-1))} ,
	{20 , 4 , 1.45627192772164776850783593749915 * pow((double)10,(double)0) , 3.41635149494490740862897232913737 * pow((double)10,(double)(-2))} ,
	{21 , 4 , 1.85908354572076815175818819518961 * pow((double)10,(double)0) , 7.25277483230795428348171496790736 * pow((double)10,(double)(-3))} ,
	{22 , 4 , 2.58904047158998381595919190652148 * pow((double)10,(double)0) , 4.58588190682845115433705554969906 * pow((double)10,(double)(-4))} ,
	{23 , 4 , 3.33004698770489421170373178710890 * pow((double)10,(double)0) , 7.57531613358637802375888601661294 * pow((double)10,(double)(-6))} ,
	{24 , 4 , 4.13460787814490154951137004441549 * pow((double)10,(double)0) , 2.66835327349558608272806587429416 * pow((double)10,(double)(-8))} ,
	{25 , 4 , 5.22910467169634386126886183446137 * pow((double)10,(double)0) , 2.37617404612396438012114434000406 * pow((double)10,(double)(-12))}
      };
  

      int somme, somme2;
      somme = 0;
      somme2 = 0;
      
      for (i=1; i<=25; i=i+1) {
	for (j=1; j<=((int)M[i-1][1]); j=j+1) {
	  weights[somme+j-1] = M[i-1][3];
	}
	
	somme2 = 0;
	for (j=1; j<=((int)M[i-1][1]/4); j=j+1)
	  {
	    for (k=1; k<=4; k=k+1)
	      {
		points[(somme+somme2+k-1) * q[0]] = M[i-1][2] * cos((k-1)*PI/2.0);
		points[(somme+somme2+k-1) * q[0] + 1] = M[i-1][2] * sin((k-1)*PI/2.0);
	      }
	    somme2 = somme2 + ((int)M[i-1][1])/4;
	  }
	somme = somme + (int)M[i-1][1];
      }
      
    }


    for (j=1; j<=N[0]; j=j+1) weights[j-1] = pow(PI,(-q[0]/2.0)) * weights[j-1];
  
    for (i=0; i<(N[0]*q[0]); i++) points[i] = sqrt(2.0) * a[0] * points[i];

  }
  // End of if q==2


	


    if (q[0] >= 3) {

      //E_q^(r^2): 7-2 (q>=3) Degree 7, 2^(q+1)+4q^2 Points dans Stroud (1971) p319
      //Attention il y a une erreur dans le livre, voir plutot Stroud (1967) : Some seventh degree
      //integration formulas for symmetric regions, SIAM J. Numer. Anal. pp. 37-44


      N[0] = (int)pow(2,q[0]+1) + 4*(int)pow(q[0],2);

  void binarycode(int i,int q, int *y);
  double Absolute(double x);
  double A1, A2, r1, r2, sval, tval;

  A1 = (q[0] + 2 + sqrt((double)(2 * (q[0] + 2)))) / (4 * (q[0] + 2)) * (gammafn(q[0] / 2.0));
  A2 = (q[0] + 2 - sqrt((double)(2 * (q[0] + 2)))) / (4 * (q[0] + 2)) * (gammafn(q[0] / 2.0));
  r1 = sqrt((q[0] + 2 - sqrt((double)(2 * (q[0] + 2)))) / 2);
  r2 = sqrt((q[0] + 2 + sqrt((double)(2 * (q[0] + 2)))) / 2);
  sval = 1.0 / sqrt((double)q[0]);
  tval = 1.0 / (sqrt((double)2));

  int *y;
  y = new int[q[0]];

  double *B;
  B = new double[N[0]/2]; // 2^q+2q^2

  double *pointsB;
  pointsB = new double[(2*q[0]) * q[0]]; 
  for (i=0; i<(2*q[0])*q[0]; i++) pointsB[i] = 0.0;
  double *pointsC;
  pointsC = new double[((int)pow((double)2,(double)q[0])) * q[0]]; 
  for (i=0; i<(((int)pow((double)2,(double)q[0])) * q[0]); i++) pointsC[i] = 0.0;
  double *pointsDtemp;
  pointsDtemp = new double[((q[0]*(q[0]-1))/2) * q[0]];
  for (i=0; i<(((q[0]*(q[0]-1))/2) * q[0]); i++) pointsDtemp[i] = 0.0;
  double *AbspointsDtemp;
  AbspointsDtemp = new double[((q[0]*(q[0]-1))/2) * q[0]];
  double *pointsD;
  pointsD = new double[(4*(q[0]*(q[0]-1))/2) * q[0]];
  double *uj;
  uj = new double[( N[0]/2 ) * q[0]];

  

  for (j=1;j<=q[0];j=j+1) {
    pointsB[(j-1)*(2*q[0])+j-1] = 1;
    pointsB[(j-1)*(2*q[0])+q[0]+j-1] = -1;
  }


  for (i=0;i<=(int)pow((double)2,(double)q[0])-1;i=i+1) {
      binarycode(i,q[0],y);
      for (j=1;j<=q[0];j=j+1) pointsC[(j-1)*((int)pow((double)2,(double)q[0]))+i] = sval * pow((double)(-1),y[j-1]);
    }
  delete y;

      
for (i=1;i<=q[0]-1;i=i+1) {
        
  for (j=1;j<=q[0]-i;j=j+1) {
                
    pointsDtemp[(i-1)*((q[0]*(q[0]-1))/2) + (q[0]*(i-1)-i*(i-1)/2)+j-1] = -tval;
    
    if (i+j<q[0]+1) {
        
      pointsDtemp[(i+j-1)*((q[0]*(q[0]-1))/2) + (q[0]*(i-1)-i*(i-1)/2)+j-1] = tval; }
    
  }
 }

 

 for (i=1;i<=(q[0]*(q[0]-1))/2;i=i+1){
   for (j=1;j<=q[0];j=j+1){
     AbspointsDtemp[(j-1)*((q[0]*(q[0]-1))/2)+i-1] = Absolute(pointsDtemp[(j-1)*((q[0]*(q[0]-1))/2)+i-1]);
   }
 }

 

 for (j=1; j<=q[0]; j++) {
   for (i=1; i<=((q[0]*(q[0]-1))/2); i++) {
     pointsD[(j-1)*(4*(q[0]*(q[0]-1))/2) +i-1] = pointsDtemp[(j-1)*((q[0]*(q[0]-1))/2)+i-1];
     pointsD[(j-1)*(4*(q[0]*(q[0]-1))/2) +i-1 + ((q[0]*(q[0]-1))/2)] = -pointsDtemp[(j-1)*((q[0]*(q[0]-1))/2)+i-1];
     pointsD[(j-1)*(4*(q[0]*(q[0]-1))/2) +i-1 + 2*((q[0]*(q[0]-1))/2)] = AbspointsDtemp[(j-1)*((q[0]*(q[0]-1))/2)+i-1];
     pointsD[(j-1)*(4*(q[0]*(q[0]-1))/2) +i-1 + 3*((q[0]*(q[0]-1))/2)] = -AbspointsDtemp[(j-1)*((q[0]*(q[0]-1))/2)+i-1];
   }
   for (i=1; i<=2*q[0]; i++) {
     uj[(j-1)*N[0]/2+i-1] = pointsB[(j-1)*2*q[0]+i-1];
   }
   for (i=1; i<=((int)pow((double)2,(double)q[0])); i++) {
     uj[(j-1)*N[0]/2+i-1+2*q[0]] = pointsC[(j-1)*((int)pow((double)2,(double)q[0]))+i-1];
   }
   for (i=1; i<=(4*(q[0]*(q[0]-1))/2); i++) {
     uj[(j-1)*N[0]/2+i-1+2*q[0]+((int)pow((double)2,(double)q[0]))] = pointsD[(j-1)*(4*(q[0]*(q[0]-1))/2)+i-1];
   }

 }

 



 for (j=0; j<q[0]; j++) {
   for (i=0; i<(N[0]/2); i++) {
     points[j * N[0] + i] = r1 * sqrt(2.0) * a[0] * uj[j * N[0]/2 + i];
     points[j * N[0] + i + N[0]/2] = r2 * sqrt(2.0) * a[0] * uj[j * N[0]/2 + i];
   }
 }

 


 for (j=1;j<=2*q[0];j=j+1) B[j-1] = (8.0-q[0])/(q[0]*(q[0]+2)*(q[0]+4));
 for (j=1;j<=(int)pow((double)2,(double)q[0]);j=j+1) B[2*q[0]+j-1] = (pow(2.0,-q[0]) * pow((double)q[0],(double)3))/(q[0]*(q[0]+2)*(q[0]+4));
 for (j=1;j<=2*(int)pow((double)q[0],(double)2)-2*q[0];j=j+1) B[2*q[0]+(int)pow((double)2,(double)q[0])+j-1] = 4.0/(q[0]*(q[0]+2)*(q[0]+4));
 for (j=1; j<=(N[0]/2); j=j+1) {
   weights[j-1] = A1 * 2 * B[j-1] / ((gammafn(q[0]/2.0)));
   weights[N[0]/2+j-1] = A2 * 2 * B[j-1] / ((gammafn(q[0]/2.0)));
 }

    }
    // End of if q>=3


  //We output the results
    for (i=1; i<=N[0]; i=i+1) {
      resweights[i-1] = weights[i-1];
      for (j=1; j<=q[0]; j=j+1) { 
	respoints[(j-1)*N[0]+(i-1)] = points[(j-1)*N[0]+(i-1)];
      }
    }


    delete[] weights;
    delete[] points;

   return;

}



double Absolute(double x)

{

  
  if (x>0) return(x);
  else return(-x);

}

  }
