//#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "mycomplex.h"

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void CnhatmatC(void *, void *, void *, void *, void *, void *, void *, void *);
extern void CnhatmatClower(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void cubature1(void *, void *, void *, void *, void *, void *, void *, void *);
extern void Dcov1Cnormed(void *, void *, void *, void *, void *, void *, void *, void *);
extern void Dcov2Cnormed(void *, void *, void *, void *, void *, void *, void *, void *);
extern void dependogramC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void zhpevxC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Fortran calls */
// extern void F77_NAME(zhpevx)(const void *, const void *, const void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"CnhatmatC",      (DL_FUNC) &CnhatmatC,       8},
    {"CnhatmatClower", (DL_FUNC) &CnhatmatClower,  9},
    {"cubature1",      (DL_FUNC) &cubature1,       8},
    {"Dcov1Cnormed",   (DL_FUNC) &Dcov1Cnormed,    8},
    {"Dcov2Cnormed",   (DL_FUNC) &Dcov2Cnormed,    8},
    {"dependogramC",    (DL_FUNC) &dependogramC,    14},
    {"zhpevxC",        (DL_FUNC) &zhpevxC,        19},
    {NULL, NULL, 0}
};

/*static const R_FortranMethodDef FortranEntries[] = {
    {"zhpevx", (DL_FUNC) &F77_NAME(zhpevx), 19},
    {NULL, NULL, 0}
    };*/

void R_init_IndependenceTests(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
