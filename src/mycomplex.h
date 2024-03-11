/////////////////////////////////////////////////////////////////////////
// macros to allow us to use either C++ or C (with C99 features)

#ifdef __cplusplus

# include <complex>
#  include <cfloat>
#  include <cmath>
#  include <limits>
using namespace std;

// use std::numeric_limits, since 1./0. and 0./0. fail with some compilers (MS)
#  define Inf numeric_limits<double>::infinity()
#  define NaN numeric_limits<double>::quiet_NaN()

typedef complex<double> cmplx;

// Use C-like complex syntax, since the C syntax is more restrictive
#  define cexp(z) exp((std::complex<double>)z)
#  define creal(z) real(z)
#  define cimag(z) imag(z)
#  define cpolar(r,t) polar(r,t)

#  define C(a,b) cmplx(a,b)

// isnan/isinf were introduced in C++11
#  if (__cplusplus < 201103L) && (!defined(HAVE_ISNAN) || !defined(HAVE_ISINF))
static inline bool my_isnan(double x) { return x != x; }
#    define isnan my_isnan
static inline bool my_isinf(double x) { return 1/x == 0.; }
#    define isinf my_isinf
#  elif (__cplusplus >= 201103L)
// g++ gets confused between the C and C++ isnan/isinf functions
#    define isnan std::isnan
#    define isinf std::isinf
#  endif

// copysign was introduced in C++11 (and is also in POSIX and C99)
#  if defined(_WIN32) || defined(__WIN32__)
#    define copysign _copysign // of course MS had to be different
#  elif defined(GNULIB_NAMESPACE) // we are using using gnulib <cmath>
#    define copysign GNULIB_NAMESPACE::copysign
#  elif (__cplusplus < 201103L) && !defined(HAVE_COPYSIGN) && !defined(__linux__) && !(defined(__APPLE__) && defined(__MACH__)) && !defined(_AIX)
static inline double my_copysign(double x, double y) { return x<0 != y<0 ? -x : x; }
#    define copysign my_copysign
#  endif

// If we are using the gnulib <cmath> (e.g. in the GNU Octave sources),
// gnulib generates a link warning if we use ::floor instead of gnulib::floor.
// This warning is completely innocuous because the only difference between
// gnulib::floor and the system ::floor (and only on ancient OSF systems)
// has to do with floor(-0), which doesn't occur in the usage below, but
// the Octave developers prefer that we silence the warning.
#  ifdef GNULIB_NAMESPACE
#    define floor GNULIB_NAMESPACE::floor
#  endif

#else // !__cplusplus, i.e. pure C (requires C99 features)

#  define _GNU_SOURCE // enable GNU libc NAN extension if possible

# include <complex.h>
#  include <float.h>
#  include <math.h>

typedef double complex cmplx;

/* Constructing complex numbers like 0+i*NaN is problematic in C99
   without the C11 CMPLX macro, because 0.+I*NAN may give NaN+i*NAN if
   I is a complex (rather than imaginary) constant.  For some reason,
   however, it works fine in (pre-4.7) gcc if I define Inf and NaN as
   1/0 and 0/0 (and only if I compile with optimization -O1 or more),
   but not if I use the INFINITY or NAN macros. */

/* __builtin_complex was introduced in gcc 4.7, but the C11 CMPLX macro
   may not be defined unless we are using a recent (2012) version of
   glibc and compile with -std=c11... note that icc lies about being
   gcc and probably doesn't have this builtin(?), so exclude icc explicitly */
#  if !defined(CMPLX) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7)) && !(defined(__ICC) || defined(__INTEL_COMPILER))
#    define CMPLX(a,b) __builtin_complex((double) (a), (double) (b))
#  endif

#  ifdef CMPLX // C11
#    define C(a,b) CMPLX(a,b)
#    define Inf INFINITY // C99 infinity
#    ifdef NAN // GNU libc extension
#      define NaN NAN
#    else
#      define NaN (0./0.) // NaN
#    endif
#  else
#    define C(a,b) ((a) + I*(b))
#    define Inf (1./0.) 
#    define NaN (0./0.) 
#  endif

static inline cmplx cpolar(double r, double t)
{
  if (r == 0.0 && !isnan(t))
    return 0.0;
  else
    return C(r * cos(t), r * sin(t));
}

#endif // !__cplusplus, i.e. pure C (requires C99 features)

