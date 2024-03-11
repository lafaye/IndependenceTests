*     > \brief <b> ZHPEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
*     
*     =========== DOCUMENTATION ===========
*     
*     Online html documentation available at
*     http://www.netlib.org/lapack/explore-html/
*     
*     > \htmlonly
*     > Download ZHPEVX + dependencies
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpevx.f">
*     > [TGZ]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpevx.f">
*     > [ZIP]</a>
*     > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpevx.f">
*     > [TXT]</a>
*     > \endhtmlonly
*     
*     Definition:
*     ===========
*     
*     SUBROUTINE ZHPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,
*     ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK,
*     IFAIL, INFO )
*     
*     .. Scalar Arguments ..
*     CHARACTER JOBZ, RANGE, UPLO
*     INTEGER IL, INFO, IU, LDZ, M, N
*     DOUBLE PRECISION ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
*     INTEGER IFAIL( * ), IWORK( * )
*     DOUBLE PRECISION RWORK( * ), W( * )
*     COMPLEX*16 AP( * ), WORK( * ), Z( LDZ, * )
*     ..
*     
*     
*     > \par Purpose:
*     =============
*     >
*     > \verbatim
*     >
*     > ZHPEVX computes selected eigenvalues and, optionally, eigenvectors
*     > of a complex Hermitian matrix A in packed storage.
*     > Eigenvalues/vectors can be selected by specifying either a range of
*     > values or a range of indices for the desired eigenvalues.
*     > \endverbatim
*     
*     Arguments:
*     ==========
*     
*     > \param[in] JOBZ
*     > \verbatim
*     > JOBZ is CHARACTER*1
*     > = 'N': Compute eigenvalues only;
*     > = 'V': Compute eigenvalues and eigenvectors.
*     > \endverbatim
*     >
*     > \param[in] RANGE
*     > \verbatim
*     > RANGE is CHARACTER*1
*     > = 'A': all eigenvalues will be found;
*     > = 'V': all eigenvalues in the half-open interval (VL,VU]
*     > will be found;
*     > = 'I': the IL-th through IU-th eigenvalues will be found.
*     > \endverbatim
*     >
*     > \param[in] UPLO
*     > \verbatim
*     > UPLO is CHARACTER*1
*     > = 'U': Upper triangle of A is stored;
*     > = 'L': Lower triangle of A is stored.
*     > \endverbatim
*     >
*     > \param[in] N
*     > \verbatim
*     > N is INTEGER
*     > The order of the matrix A. N >= 0.
*     > \endverbatim
*     >
*     > \param[in,out] AP
*     > \verbatim
*     > AP is COMPLEX*16 array, dimension (N*(N+1)/2)
*     > On entry, the upper or lower triangle of the Hermitian matrix
*     > A, packed columnwise in a linear array. The j-th column of A
*     > is stored in the array AP as follows:
*     > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*     > if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
*     >
*     > On exit, AP is overwritten by values generated during the
*     > reduction to tridiagonal form. If UPLO = 'U', the diagonal
*     > and first superdiagonal of the tridiagonal matrix T overwrite
*     > the corresponding elements of A, and if UPLO = 'L', the
*     > diagonal and first subdiagonal of T overwrite the
*     > corresponding elements of A.
*     > \endverbatim
*     >
*     > \param[in] VL
*     > \verbatim
*     > VL is DOUBLE PRECISION
*     > \endverbatim
*     >
*     > \param[in] VU
*     > \verbatim
*     > VU is DOUBLE PRECISION
*     > If RANGE='V', the lower and upper bounds of the interval to
*     > be searched for eigenvalues. VL < VU.
*     > Not referenced if RANGE = 'A' or 'I'.
*     > \endverbatim
*     >
*     > \param[in] IL
*     > \verbatim
*     > IL is INTEGER
*     > \endverbatim
*     >
*     > \param[in] IU
*     > \verbatim
*     > IU is INTEGER
*     > If RANGE='I', the indices (in ascending order) of the
*     > smallest and largest eigenvalues to be returned.
*     > 1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
*     > Not referenced if RANGE = 'A' or 'V'.
*     > \endverbatim
*     >
*     > \param[in] ABSTOL
*     > \verbatim
*     > ABSTOL is DOUBLE PRECISION
*     > The absolute error tolerance for the eigenvalues.
*     > An approximate eigenvalue is accepted as converged
*     > when it is determined to lie in an interval [a,b]
*     > of width less than or equal to
*     >
*     > ABSTOL + EPS * max( |a|,|b| ) ,
*     >
*     > where EPS is the machine precision. If ABSTOL is less than
*     > or equal to zero, then EPS*|T| will be used in its place,
*     > where |T| is the 1-norm of the tridiagonal matrix obtained
*     > by reducing AP to tridiagonal form.
*     >
*     > Eigenvalues will be computed most accurately when ABSTOL is
*     > set to twice the underflow threshold 2*DLAMCH('S'), not zero.
*     > If this routine returns with INFO>0, indicating that some
*     > eigenvectors did not converge, try setting ABSTOL to
*     > 2*DLAMCH('S').
*     >
*     > See "Computing Small Singular Values of Bidiagonal Matrices
*     > with Guaranteed High Relative Accuracy," by Demmel and
*     > Kahan, LAPACK Working Note #3.
*     > \endverbatim
*     >
*     > \param[out] M
*     > \verbatim
*     > M is INTEGER
*     > The total number of eigenvalues found. 0 <= M <= N.
*     > If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*     > \endverbatim
*     >
*     > \param[out] W
*     > \verbatim
*     > W is DOUBLE PRECISION array, dimension (N)
*     > If INFO = 0, the selected eigenvalues in ascending order.
*     > \endverbatim
*     >
*     > \param[out] Z
*     > \verbatim
*     > Z is COMPLEX*16 array, dimension (LDZ, max(1,M))
*     > If JOBZ = 'V', then if INFO = 0, the first M columns of Z
*     > contain the orthonormal eigenvectors of the matrix A
*     > corresponding to the selected eigenvalues, with the i-th
*     > column of Z holding the eigenvector associated with W(i).
*     > If an eigenvector fails to converge, then that column of Z
*     > contains the latest approximation to the eigenvector, and
*     > the index of the eigenvector is returned in IFAIL.
*     > If JOBZ = 'N', then Z is not referenced.
*     > Note: the user must ensure that at least max(1,M) columns are
*     > supplied in the array Z; if RANGE = 'V', the exact value of M
*     > is not known in advance and an upper bound must be used.
*     > \endverbatim
*     >
*     > \param[in] LDZ
*     > \verbatim
*     > LDZ is INTEGER
*     > The leading dimension of the array Z. LDZ >= 1, and if
*     > JOBZ = 'V', LDZ >= max(1,N).
*     > \endverbatim
*     >
*     > \param[out] WORK
*     > \verbatim
*     > WORK is COMPLEX*16 array, dimension (2*N)
*     > \endverbatim
*     >
*     > \param[out] RWORK
*     > \verbatim
*     > RWORK is DOUBLE PRECISION array, dimension (7*N)
*     > \endverbatim
*     >
*     > \param[out] IWORK
*     > \verbatim
*     > IWORK is INTEGER array, dimension (5*N)
*     > \endverbatim
*     >
*     > \param[out] IFAIL
*     > \verbatim
*     > IFAIL is INTEGER array, dimension (N)
*     > If JOBZ = 'V', then if INFO = 0, the first M elements of
*     > IFAIL are zero. If INFO > 0, then IFAIL contains the
*     > indices of the eigenvectors that failed to converge.
*     > If JOBZ = 'N', then IFAIL is not referenced.
*     > \endverbatim
*     >
*     > \param[out] INFO
*     > \verbatim
*     > INFO is INTEGER
*     > = 0: successful exit
*     > < 0: if INFO = -i, the i-th argument had an illegal value
*     > > 0: if INFO = i, then i eigenvectors failed to converge.
*     > Their indices are stored in array IFAIL.
*     > \endverbatim
*     
*     Authors:
*     ========
*     
*     > \author Univ. of Tennessee
*     > \author Univ. of California Berkeley
*     > \author Univ. of Colorado Denver
*     > \author NAG Ltd.
*     
*     > \date November 2011
*     
*     > \ingroup complex16OTHEReigen
*     
*     =====================================================================
      SUBROUTINE zhpevx(JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,
     $     abstol, m, w, z, ldz, work, rwork, iwork,
     $     ifail, info )
*     
*     -- LAPACK driver routine (version 3.4.0) --
*     -- LAPACK is a software package provided by Univ. of Tennessee, --
*     -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*     
*     .. Scalar Arguments ..
      CHARACTER jobz, range, uplo
      INTEGER il, info, iu, ldz, m, n
      DOUBLE PRECISION abstol, vl, vu
*     ..
*     .. Array Arguments ..
      INTEGER ifail( * ), iwork( * )
      DOUBLE PRECISION rwork( * ), w( * )
      COMPLEX*16 ap( * ), work( * ), z( ldz, * )
*     ..
*     
*     =====================================================================
*     
*     .. Parameters ..
      DOUBLE PRECISION zero, one
      parameter( zero = 0.0d0, one = 1.0d0 )
      COMPLEX*16 cone
      parameter( cone = ( 1.0d0, 0.0d0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL alleig, indeig, test, valeig, wantz
      CHARACTER order
      INTEGER i, iinfo, imax, indd, inde, indee, indibl,
     $     indisp, indiwk, indrwk, indtau, indwrk, iscale,
     $     itmp1, j, jj, nsplit
      DOUBLE PRECISION abstll, anrm, bignum, eps, rmax, rmin, safmin,
     $     sigma, smlnum, tmp1, vll, vuu
*     ..
*     .. External Functions ..
      LOGICAL lsame
      DOUBLE PRECISION dlamch, zlanhp
      EXTERNAL lsame, dlamch, zlanhp
*     ..
*     .. External Subroutines ..
      EXTERNAL dcopy, dscal, dstebz, dsterf, xerbla, zdscal,
     $     zhptrd, zstein, zsteqr, zswap, zupgtr, zupmtr
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dble, max, min, sqrt
*     ..
*     .. Executable Statements ..
*     
*     Test the input parameters.
*     
      wantz = lsame( jobz, 'V' )
      alleig = lsame( range, 'A' )
      valeig = lsame( range, 'V' )
      indeig = lsame( range, 'I' )
*     
      info = 0
      IF( .NOT.( wantz .OR. lsame( jobz, 'N' ) ) ) THEN
         info = -1
      ELSE IF( .NOT.( alleig .OR. valeig .OR. indeig ) ) THEN
         info = -2
      ELSE IF( .NOT.( lsame( uplo, 'L' ) .OR. lsame( uplo, 'U' ) ) )
     $        THEN
         info = -3
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE
         IF( valeig ) THEN
            IF( n.GT.0 .AND. vu.LE.vl )
     $           info = -7
         ELSE IF( indeig ) THEN
            IF( il.LT.1 .OR. il.GT.max( 1, n ) ) THEN
               info = -8
            ELSE IF( iu.LT.min( n, il ) .OR. iu.GT.n ) THEN
               info = -9
            END IF
         END IF
      END IF
      IF( info.EQ.0 ) THEN
         IF( ldz.LT.1 .OR. ( wantz .AND. ldz.LT.n ) )
     $        info = -14
      END IF
*     
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZHPEVX', -info )
         RETURN
      END IF
*     
*     Quick return if possible
*     
      m = 0
      IF( n.EQ.0 )
     $     RETURN
*     
      IF( n.EQ.1 ) THEN
         IF( alleig .OR. indeig ) THEN
            m = 1
            w( 1 ) = DBLE(ap( 1 ))
         ELSE
            IF( vl.LT.dble( ap( 1 ) ) .AND. vu.GE.dble( ap( 1 ) ) ) THEN
               m = 1
               w( 1 ) = DBLE(ap( 1 ))
            END IF
         END IF
         IF( wantz )
     $        z( 1, 1 ) = cone
         RETURN
      END IF
*     
*     Get machine constants.
*     
      safmin = dlamch( 'Safe minimum' )
      eps = dlamch( 'Precision' )
      smlnum = safmin / eps
      bignum = one / smlnum
      rmin = sqrt( smlnum )
      rmax = min( sqrt( bignum ), one / sqrt( sqrt( safmin ) ) )
*     
*     Scale matrix to allowable range, if necessary.
*     
      iscale = 0
      abstll = abstol
      IF( valeig ) THEN
         vll = vl
         vuu = vu
      ELSE
         vll = zero
         vuu = zero
      END IF
      anrm = zlanhp( 'M', uplo, n, ap, rwork )
      IF( anrm.GT.zero .AND. anrm.LT.rmin ) THEN
         iscale = 1
         sigma = rmin / anrm
      ELSE IF( anrm.GT.rmax ) THEN
         iscale = 1
         sigma = rmax / anrm
      END IF
      IF( iscale.EQ.1 ) THEN
         CALL zdscal( ( n*( n+1 ) ) / 2, sigma, ap, 1 )
         IF( abstol.GT.0 )
     $        abstll = abstol*sigma
         IF( valeig ) THEN
            vll = vl*sigma
            vuu = vu*sigma
         END IF
      END IF
*     
*     Call ZHPTRD to reduce Hermitian packed matrix to tridiagonal form.
*     
      indd = 1
      inde = indd + n
      indrwk = inde + n
      indtau = 1
      indwrk = indtau + n
      CALL zhptrd( uplo, n, ap, rwork( indd ), rwork( inde ),
     $     work( indtau ), iinfo )
*     
*     If all eigenvalues are desired and ABSTOL is less than or equal
*     to zero, then call DSTERF or ZUPGTR and ZSTEQR. If this fails
*     for some eigenvalue, then try DSTEBZ.
*     
      test = .false.
      IF (indeig) THEN
         IF (il.EQ.1 .AND. iu.EQ.n) THEN
            test = .true.
         END IF
      END IF
      IF ((alleig .OR. test) .AND. (abstol.LE.zero)) THEN
         CALL dcopy( n, rwork( indd ), 1, w, 1 )
         indee = indrwk + 2*n
         IF( .NOT.wantz ) THEN
            CALL dcopy( n-1, rwork( inde ), 1, rwork( indee ), 1 )
            CALL dsterf( n, w, rwork( indee ), info )
         ELSE
            CALL zupgtr( uplo, n, ap, work( indtau ), z, ldz,
     $           work( indwrk ), iinfo )
            CALL dcopy( n-1, rwork( inde ), 1, rwork( indee ), 1 )
            CALL zsteqr( jobz, n, w, rwork( indee ), z, ldz,
     $           rwork( indrwk ), info )
            IF( info.EQ.0 ) THEN
               DO 10 i = 1, n
                  ifail( i ) = 0
 10            CONTINUE
            END IF
         END IF
         IF( info.EQ.0 ) THEN
            m = n
            go to 20
         END IF
         info = 0
      END IF
*     
*     Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN.
*     
      IF( wantz ) THEN
         order = 'B'
      ELSE
         order = 'E'
      END IF
      indibl = 1
      indisp = indibl + n
      indiwk = indisp + n
      CALL dstebz( range, order, n, vll, vuu, il, iu, abstll,
     $     rwork( indd ), rwork( inde ), m, nsplit, w,
     $     iwork( indibl ), iwork( indisp ), rwork( indrwk ),
     $     iwork( indiwk ), info )
*     
      IF( wantz ) THEN
         CALL zstein( n, rwork( indd ), rwork( inde ), m, w,
     $        iwork( indibl ), iwork( indisp ), z, ldz,
     $        rwork( indrwk ), iwork( indiwk ), ifail, info )
*     
*     Apply unitary matrix used in reduction to tridiagonal
*     form to eigenvectors returned by ZSTEIN.
*     
         indwrk = indtau + n
         CALL zupmtr( 'L', uplo, 'N', n, m, ap, work( indtau ), z, ldz,
     $        work( indwrk ), iinfo )
      END IF
*     
*     If matrix was scaled, then rescale eigenvalues appropriately.
*     
 20   CONTINUE
      IF( iscale.EQ.1 ) THEN
         IF( info.EQ.0 ) THEN
            imax = m
         ELSE
            imax = info - 1
         END IF
         CALL dscal( imax, one / sigma, w, 1 )
      END IF
*     
*     If eigenvalues are not in order, then sort them, along with
*     eigenvectors.
*     
      IF( wantz ) THEN
         DO 40 j = 1, m - 1
            i = 0
            tmp1 = w( j )
            DO 30 jj = j + 1, m
               IF( w( jj ).LT.tmp1 ) THEN
                  i = jj
                  tmp1 = w( jj )
               END IF
 30         CONTINUE
*     
            IF( i.NE.0 ) THEN
               itmp1 = iwork( indibl+i-1 )
               w( i ) = w( j )
               iwork( indibl+i-1 ) = iwork( indibl+j-1 )
               w( j ) = tmp1
               iwork( indibl+j-1 ) = itmp1
               CALL zswap( n, z( 1, i ), 1, z( 1, j ), 1 )
               IF( info.NE.0 ) THEN
                  itmp1 = ifail( i )
                  ifail( i ) = ifail( j )
                  ifail( j ) = itmp1
               END IF
            END IF
 40      CONTINUE
      END IF
*     
      RETURN
*     
*     End of ZHPEVX
*     
      END
