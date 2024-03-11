#define USE_FC_LEN_T
#include <R.h>
#include "Rmath.h"


#ifdef FC_LEN_T

extern "C" {
  
  void zhpevxC(char **jobz, char **range, char **uplo, int *n, Rcomplex *ap,
	       double *vl, double *vu, int *il, int *iu, double *abstol, int *m,
	       double *w, Rcomplex *z, int *ldz, Rcomplex *work, double *rwork,
	       int *iwork, int *ifail, int *info) {
    
    
    /*
    char cjobz[2];
    strncpy(cjobz, jobz[0], 1);
    cjobz[1] = '\0';
    char crange[2];
    strncpy(crange, range[0], 1);
    crange[1] = '\0';
    char cuplo[2];
    strncpy(cuplo, uplo[0], 1);
    cuplo[1] = '\0';
    */

    char *cjobz;
    cjobz = new char[2];
    strncpy(cjobz, jobz[0], 1);
    cjobz[1] = '\0';
    char *crange;
    crange = new char[2];
    strncpy(crange, range[0], 1);
    crange[1] = '\0';
    char *cuplo;
    cuplo = new char[2];
    strncpy(cuplo, uplo[0], 1);
    cuplo[1] = '\0';

       
    double _Complex *Cap;
    Cap = new double _Complex[sizeof(ap)];
    memcpy(&Cap, &ap, sizeof(ap));
    double _Complex *Cz;
    Cz = new double _Complex[sizeof(z)];
    memcpy(&Cz, &z, sizeof(z));
    double _Complex *Cwork;
    Cwork = new double _Complex[sizeof(work)];
    memcpy(&Cwork, &work, sizeof(work));

    
    void F77_NAME(zhpevx)(const char *cjobz, const char *crange, const char *cuplo,
			  const int *n, __complex__ double *Cap, const double *vl,
			  const double *vu, const int *il, const int *iu,
			  const double *abstol, int *m, double *w,
			  __complex__ double *Cz, const int *ldz, __complex__ double *Cwork, double *rwork,
			  int *iwork, int *ifail, int *info,
			  FC_LEN_T jobz_len,  FC_LEN_T range_len,  FC_LEN_T uplo_len);
    
    
    
    F77_CALL(zhpevx)(cjobz, crange, cuplo, &n[0], Cap, &vl[0], &vu[0], &il[0], &iu[0], &abstol[0], &m[0],
		     w, Cz, &ldz[0], Cwork, rwork, iwork, ifail, &info[0], strlen(cjobz), strlen(crange), strlen(cuplo));
    

    delete[] cjobz;
    delete[] crange;
    delete[] cuplo;
    
    delete[] Cap;
    delete[] Cz;
    delete[] Cwork;	  

  }
  
}
#else
extern "C" {
  
  void zhpevxC(char **jobz, char **range, char **uplo, int *n, Rcomplex *ap,
	       double *vl, double *vu, int *il, int *iu, double *abstol, int *m,
	       double *w, Rcomplex *z, int *ldz, Rcomplex *work, double *rwork,
	       int *iwork, int *ifail, int *info) {
    
    extern void F77_NAME(zhpevx)(const char *cjobz, const char *crange, const char *cuplo,
				 const int *n, __complex__ double *Cap, const double *vl,
				 const double *vu, const int *il, const int *iu,
				 const double *abstol, int *m, double *w,
				 __complex__ double *Cz, const int *ldz, __complex__ double *Cwork, double *rwork,
				 int *iwork, int *ifail, int *info);
    
    const char *CJOBZ = jobz[0];
    const char *CRANGE = range[0];
    const char *CUPLO = uplo[0];

    double _Complex *Cap;
    Cap = new double _Complex[sizeof(ap)];
    memcpy(&Cap, &ap, sizeof(ap));
    double _Complex *Cz;
    Cz = new double _Complex[sizeof(z)];
    memcpy(&Cz, &z, sizeof(z));
    double _Complex *Cwork;
    Cwork = new double _Complex[sizeof(work)];
    memcpy(&Cwork, &work, sizeof(work));
    
    F77_CALL(zhpevx)(CJOBZ, CRANGE, CUPLO, &n[0], Cap, &vl[0], &vu[0], &il[0], &iu[0], &abstol[0], &m[0],
		     w, Cz, &ldz[0], Cwork, rwork, iwork, ifail, &info[0]);

    
    delete[] Cap;
    delete[] Cz;
    delete[] Cwork;
  }
  
}
#endif

