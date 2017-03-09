      SUBROUTINE PDA_LSQR  ( M, N, DAMP, LENIW, LENJW, LENRW, IW,
     :                       JW, RW,
     :                       U, V, W, X, SE, ATOL, BTOL, CONLIM, ITNLIM,
     :                       ISTOP, ITN, ANORM, ACOND, RNORM, ARNORM,
     :                       XNORM )

      Implicit NONE
      INTEGER        M, N, LENIW,LENJW, LENRW, ITNLIM, ISTOP, ITN
      INTEGER        IW(LENIW), JW(LENJW)
      DOUBLE PRECISION    RW(LENRW), U(M), V(N), W(N), X(N), SE(N),
     :                   ATOL, BTOL, CONLIM, DAMP,
     :                   ANORM, ACOND, RNORM, ARNORM, XNORM
*-----------------------------------------------------------------------
*
*     PDA_LSQR  finds a solution x to the following problems:
*
*     1. Unsymmetric equations --    solve  A*x = b
*
*     2. Linear least squares  --    solve  A*x = b
*                                    in the least-squares sense
*
*     3. Damped least squares  --    solve  (   A    )*x = ( b )
*                                           ( damp*I )     ( 0 )
*                                    in the least-squares sense
*
*     where A is a matrix with m rows and n columns, b is an
*     m-vector, and damp is a scalar.  (All quantities are real.)
*     The matrix A is intended to be large and sparse.  It is accessed
*     by means of subroutine calls of the form
*
*              CALL APROD ( mode, m, n, x, y, LENIW, LENJW LENRW, IW, JW, RW )
*
*     which must perform the following functions:
*
*                If MODE = 1, compute  y = y + A*x.
*                If MODE = 2, compute  x = x + A(transpose)*y.
*
*     The vectors x and y are input parameters in both cases.
*     If  mode = 1,  y should be altered without changing x.
*     If  mode = 2,  x should be altered without changing y.
*     The parameters LENIW, LENRW, IW, RW may be used for workspace
*     as described below.
*
*     The rhs vector b is input via U, and subsequently overwritten.
*
*
*     Note:  PDA_LSQR uses an iterative method to approximate the solution.
*     The number of iterations required to reach a certain accuracy
*     depends strongly on the scaling of the problem.  Poor scaling of
*     the rows or columns of A should therefore be avoided where
*     possible.
*
*     For example, in problem 1 the solution is unaltered by
*     row-scaling.  If a row of A is very small or large compared to
*     the other rows of A, the corresponding row of ( A  b ) should be
*     scaled up or down.
*
*     In problems 1 and 2, the solution x is easily recovered
*     following column-scaling.  Unless better information is known,
*     the nonzero columns of A should be scaled so that they all have
*     the same Euclidean norm (e.g., 1.0).
*
*     In problem 3, there is no freedom to re-scale if damp is
*     nonzero.  However, the value of damp should be assigned only
*     after attention has been paid to the scaling of A.
*
*     The parameter damp is intended to help regularize
*     ill-conditioned systems, by preventing the true solution from
*     being very large.  Another aid to regularization is provided by
*     the parameter ACOND, which may be used to terminate iterations
*     before the computed solution becomes very large.
*
*
*     Notation
*     --------
*
*     The following quantities are used in discussing the subroutine
*     parameters:
*
*     Abar   =  (   A    ),          bbar  =  ( b )
*               ( damp*I )                    ( 0 )
*
*     r      =  b  -  A*x,           rbar  =  bbar  -  Abar*x
*
*     rnorm  =  sqrt( norm(r)**2  +  damp**2 * norm(x)**2 )
*            =  norm( rbar )
*
*     RELPR  =  the relative precision of floating-point arithmetic
*               on the machine being used.  For example, on the IBM 370,
*               RELPR is about 1.0E-6 and 1.0D-16 in single and double
*               precision respectively.
*
*     PDA_LSQR  minimizes the function rnorm with respect to x.
*
*
*     Parameters
*     ----------
*
*     M       input      m, the number of rows in A.
*
*     N       input      n, the number of columns in A.
*
*     APROD   external   See above.
*
*     DAMP    input      The damping parameter for problem 3 above.
*                        (DAMP should be 0.0 for problems 1 and 2.)
*                        If the system A*x = b is incompatible, values
*                        of DAMP in the range 0 to sqrt(RELPR)*norm(A)
*                        will probably have a negligible effect.
*                        Larger values of DAMP will tend to decrease
*                        the norm of x and reduce the number of
*                        iterations required by PDA_LSQR.
*
*                        The work per iteration and the storage needed
*                        by PDA_LSQR are the same for all values of DAMP.
*
*     LENIW   input      The length of the workspace array IW.
*     LENRW   input      The length of the workspace array RW.
*     IW      workspace  An integer array of length LENIW.
*     RW      workspace  A real array of length LENRW.
*
*             Note:  PDA_LSQR  does not explicitly use the previous four
*             parameters, but passes them to subroutine APROD for
*             possible use as workspace.  If APROD does not need
*             IW or RW, the values LENIW = 1 or LENRW = 1 should
*             be used, and the actual parameters corresponding to
*             IW or RW  may be any convenient array of suitable type.
*
*     U(M)    input      The rhs vector b.  Beware that U is
*                        over-written by PDA_LSQR.
*
*     V(N)    workspace
*     W(N)    workspace
*
*     X(N)    output     Returns the computed solution x.
*
*     SE(N)   output     Returns standard error estimates for the
*                        components of X.  For each i, SE(i) is set
*                        to the value  rnorm * sqrt( sigma(i,i) / T ),
*                        where sigma(i,i) is an estimate of the i-th
*                        diagonal of the inverse of Abar(transpose)*Abar
*                        and  T = 1      if  m .le. n,
*                             T = m - n  if  m .gt. n  and  damp = 0,
*                             T = m      if  damp .ne. 0.
*
*     ATOL    input      An estimate of the relative error in the data
*                        defining the matrix A.  For example,
*                        if A is accurate to about 6 digits, set
*                        ATOL = 1.0E-6 .
*
*     BTOL    input      An extimate of the relative error in the data
*                        defining the rhs vector b.  For example,
*                        if b is accurate to about 6 digits, set
*                        BTOL = 1.0E-6 .
*
*     CONLIM  input      An upper limit on cond(Abar), the apparent
*                        condition number of the matrix Abar.
*                        Iterations will be terminated if a computed
*                        estimate of cond(Abar) exceeds CONLIM.
*                        This is intended to prevent certain small or
*                        zero singular values of A or Abar from
*                        coming into effect and causing unwanted growth
*                        in the computed solution.
*
*                        CONLIM and DAMP may be used separately or
*                        together to regularize ill-conditioned systems.
*
*                        Normally, CONLIM should be in the range
*                        1000 to 1/RELPR.
*                        Suggested value:
*                        CONLIM = 1/(100*RELPR)  for compatible systems,
*                        CONLIM = 1/(10*sqrt(RELPR)) for least squares.
*
*             Note:  If the user is not concerned about the parameters
*             ATOL, BTOL and CONLIM, any or all of them may be set
*             to zero.  The effect will be the same as the values
*             RELPR, RELPR and 1/RELPR respectively.
*
*     ITNLIM  input      An upper limit on the number of iterations.
*                        Suggested value:
*                        ITNLIM = n/2   for well-conditioned systems
*                                       with clustered singular values,
*                        ITNLIM = 4*n   otherwise.
*
*     ISTOP   output     An integer giving the reason for termination:
*
*                0       x = 0  is the exact solution.
*                        No iterations were performed.
*
*                1       The equations A*x = b are probably
*                        compatible.  Norm(A*x - b) is sufficiently
*                        small, given the values of ATOL and BTOL.
*
*                2       The system A*x = b is probably not
*                        compatible.  A least-squares solution has
*                        been obtained that is sufficiently accurate,
*                        given the value of ATOL.
*
*                3       An estimate of cond(Abar) has exceeded
*                        CONLIM.  The system A*x = b appears to be
*                        ill-conditioned.  Otherwise, there could be an
*                        error in subroutine APROD.
*
*                4       The equations A*x = b are probably
*                        compatible.  Norm(A*x - b) is as small as
*                        seems reasonable on this machine.
*
*                5       The system A*x = b is probably not
*                        compatible.  A least-squares solution has
*                        been obtained that is as accurate as seems
*                        reasonable on this machine.
*
*                6       Cond(Abar) seems to be so large that there is
*                        no point in doing further iterations,
*                        given the precision of this machine.
*                        There could be an error in subroutine APROD.
*
*                7       The iteration limit ITNLIM was reached.
*
*     ITN     output     The number of iterations performed.
*
*     ANORM   output     An estimate of the Frobenius norm of  Abar.
*                        This is the square-root of the sum of squares
*                        of the elements of Abar.
*                        If DAMP is small and if the columns of A
*                        have all been scaled to have length 1.0,
*                        ANORM should increase to roughly sqrt(n).
*                        A radically different value for ANORM may
*                        indicate an error in subroutine APROD (there
*                        may be an inconsistency between modes 1 and 2).
*
*     ACOND   output     An estimate of cond(Abar), the condition
*                        number of Abar.  A very high value of ACOND
*                        may again indicate an error in APROD.
*
*     RNORM   output     An estimate of the final value of norm(rbar),
*                        the function being minimized (see notation
*                        above).  This will be small if A*x = b has
*                        a solution.
*
*     ARNORM  output     An estimate of the final value of
*                        norm( Abar(transpose)*rbar ), the norm of
*                        the residual for the usual normal equations.
*                        This should be small in all cases.  (ARNORM
*                        will often be smaller than the true value
*                        computed from the output vector X.)
*
*     XNORM   output     An estimate of the norm of the final
*                        solution vector X.
*
*
*     Precision
*     ---------
*
*     The number of iterations required by PDA_LSQR will usually decrease
*     if the computation is performed in higher precision.  To convert
*     PDA_LSQR between single and double precision, change the words
*                        DOUBLE PRECISION
*                        DCOPY, DNRM2, DSCAL
*     to the appropriate FORTRAN and BLAS equivalents.
*     Also change 'D+' or 'E+' in the PARAMETER statement.
*
*
*     References
*     ----------
*
*     C.C. Paige and M.A. Saunders,  PDA_LSQR: An algorithm for sparse
*          linear equations and sparse least squares,
*          ACM Transactions on Mathematical Software 8, 1 (March 1982),
*          pp. 43-71.
*
*     C.C. Paige and M.A. Saunders,  Algorithm 583, PDA_LSQR: Sparse
*          linear equations and least-squares problems,
*          ACM Transactions on Mathematical Software 8, 2 (June 1982),
*          pp. 195-209.
*
*     C.L. Lawson, R.J. Hanson, D.R. Kincaid and F.T. Krogh,
*          Basic linear algebra subprograms for Fortran usage,
*          ACM Transactions on Mathematical Software 5, 3 (Sept 1979),
*          pp. 308-323 and 324-325.
*-----------------------------------------------------------------------
*
*
*     PDA_LSQR development:
*     22 Feb 1982: LSQR sent to ACM TOMS to become Algorithm 583.
*     15 Sep 1985: Final F66 version.  LSQR sent to "misc" in netlib.
*     13 Oct 1987: Bug (Robert Davies, DSIR).  Have to delete
*                     IF ( (ONE + DDABS(T)) .LE. ONE ) GO TO 200
*                  from loop 200.  The test was an attempt to reduce
*                  underflows, but caused W(I) not to be updated.
*     17 Mar 1989: First F77 version.
*     04 May 1989: Bug (David Gay, AT&T).  When the second BETA is zero,
*                  RNORM = 0 and
*                  TEST2 = ARNORM / (ANORM * RNORM) overflows.
*                  Fixed by testing for RNORM = 0.
*     05 May 1989: Sent to "misc" in netlib.
*     Michael A. Saunders            (na.saunders @ NA-net.stanford.edu)
*     Department of Operations Research
*     Stanford University
*     Stanford, CA 94305-4022.
*
*     19 Sep 1996: Peter W. Draper. Removed NOUT argument and renamed
*                  PDA_LSQR. PDA routines may not write output.
*-----------------------------------------------------------------------

*     Intrinsics and local variables

      INTRINSIC          DABS, MOD, DSQRT
      INTEGER      ::I, NCONV, NSTOP 
      DOUBLE PRECISION         ::PDA_DNRM2
      DOUBLE PRECISION        ::ALFA, BBNORM, BETA, BNORM,
     :                   CS, CS1, CS2, CTOL, DAMPSQ, DDNORM, DELTA,
     :                   GAMMA, GAMBAR, PHI, PHIBAR, PSI,
     :                   RES1, RES2, RHO, RHOBAR, RHBAR1, RHBAR2,
     :                   RHS, RTOL, SN, SN1, SN2,
     :                   T, TAU, TEST1, TEST2, TEST3,
     :                   THETA, T1, T2, T3, XXNORM, Z, ZBAR

      INTEGER        ::NOUT
      DOUBLE PRECISION ::ZERO,  ONE

      CHARACTER(len=16) :: ENTER, EXIT
      CHARACTER(len=60) :: MSG(0:7)

      DATA               ENTER /' Enter LSQR.    '/,
     :                   EXIT  /' Exit  LSQR.    '/

      DATA               MSG
     : / 'The exact solution is  X = 0',
     :   'Ax - b is small enough, given ATOL, BTOL',
     :   'The least-squares solution is good enough, given ATOL',
     :   'The estimate of cond(Abar) has exceeded CONLIM',
     :   'Ax - b is small enough for this machine',
     :   'The least-squares solution is good enough for this machine',
     :   'Cond(Abar) seems to be too large for this machine',
     :   'The iteration limit has been reached' /
*-----------------------------------------------------------------------

*     Initialize.
      ZERO = 0.0
      ONE = 1.0
      NOUT = -1
      IF (NOUT .GT. 0)
     :   WRITE(NOUT, 1000) ENTER, M, N, DAMP, ATOL, CONLIM, BTOL, ITNLIM
      ITN    =   0
      ISTOP  =   0
      NSTOP  =   0
      CTOL   =   ZERO
      IF (CONLIM .GT. ZERO) CTOL = ONE / CONLIM
      ANORM  =   ZERO
      ACOND  =   ZERO
      BBNORM =   ZERO
      DAMPSQ =   DAMP**2
      DDNORM =   ZERO
      RES2   =   ZERO
      XNORM  =   ZERO
      XXNORM =   ZERO
      CS2    = - ONE
      SN2    =   ZERO
      Z      =   ZERO

      DO 10  I = 1, N
         V(I)  =  ZERO
         X(I)  =  ZERO
        SE(I)  =  ZERO
   10 CONTINUE

*     Set up the first vectors U and V for the bidiagonalization.
*     These satisfy  BETA*U = b,  ALFA*V = A(transpose)*U.

      ALFA   =   ZERO
      BETA   =   ZERO
      BETA   =   PDA_DNRM2 ( M, U, 1 )
      IF (BETA .GT. ZERO) THEN
         CALL PDA_DSCAL ( M, (ONE / BETA), U, 1 )
         CALL APROD ( 2, M, N, V, U, LENIW,LENJW, LENRW, IW,JW, RW )
         ALFA   =   PDA_DNRM2 ( N, V, 1 )
      END IF

      IF (ALFA .GT. ZERO) THEN
         CALL PDA_DSCAL ( N, (ONE / ALFA), V, 1 )
         CALL PDA_DCOPY ( N, V, 1, W, 1 )
      END IF

      ARNORM =   ALFA * BETA
      IF (ARNORM .EQ. ZERO) GO TO 800

      RHOBAR =   ALFA
      PHIBAR =   BETA
      BNORM  =   BETA
      RNORM  =   BETA

      IF (NOUT   .GT.  0  ) THEN
         IF (DAMPSQ .EQ. ZERO) THEN
             WRITE(NOUT, 1200)
         ELSE
             WRITE(NOUT, 1300)
         END IF
         TEST1  = ONE
         TEST2  = ALFA / BETA
         WRITE(NOUT, 1500) ITN, X(1), RNORM, TEST1, TEST2
         WRITE(NOUT, 1600)
      END IF

*     ------------------------------------------------------------------
*     Main iteration loop.
*     ------------------------------------------------------------------
  100 ITN    = ITN + 1

*     Perform the next step of the bidiagonalization to obtain the
*     next  BETA, U, ALFA, V.  These satisfy the relations
*                BETA*U  =  A*V  -  ALFA*U,
*                ALFA*V  =  A(transpose)*U  -  BETA*V.
      CALL PDA_DSCAL ( M, (- ALFA), U, 1 )
      CALL APROD ( 1, M, N, V, U, LENIW,LENJW, LENRW, IW,JW, RW )
      BETA   =   PDA_DNRM2 ( M, U, 1 )
      BBNORM =   BBNORM  +  ALFA**2  +  BETA**2  +  DAMPSQ
      IF (BETA .GT. ZERO) THEN
         CALL PDA_DSCAL ( M, (ONE / BETA), U, 1 )
         CALL PDA_DSCAL ( N, (- BETA), V, 1 )
         CALL APROD ( 2, M, N, V, U, LENIW,LENJW, LENRW, IW,JW, RW )
         ALFA   =   PDA_DNRM2 ( N, V, 1 )
         IF (ALFA .GT. ZERO) THEN
            CALL PDA_DSCAL ( N, (ONE / ALFA), V, 1 )
         END IF
      END IF
*     Use a plane rotation to eliminate the damping parameter.
*     This alters the diagonal (RHOBAR) of the lower-bidiagonal matrix.

      RHBAR2 = RHOBAR**2  +  DAMPSQ
      RHBAR1 = DSQRT( RHBAR2 )
      CS1    = RHOBAR / RHBAR1
      SN1    = DAMP   / RHBAR1
      PSI    = SN1 * PHIBAR
      PHIBAR = CS1 * PHIBAR

*     Use a plane rotation to eliminate the subdiagonal element (BETA)
*     of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.

      RHO    =   DSQRT( RHBAR2  +  BETA**2 )
      CS     =   RHBAR1 / RHO
      SN     =   BETA   / RHO
      THETA  =   SN * ALFA
      RHOBAR = - CS * ALFA
      PHI    =   CS * PHIBAR
      PHIBAR =   SN * PHIBAR
      TAU    =   SN * PHI

*     Update  X, W  and the standard error estimates.

      T1     =   PHI   / RHO
      T2     = - THETA / RHO
      T3     =   ONE   / RHO

      DO 200  I =  1, N
         T      =  W(I)
         X(I)   =  T1*T  +  X(I)
         W(I)   =  T2*T  +  V(I)
         T      = (T3*T)**2
         SE(I)  =  T     +  SE(I)
         DDNORM =  T     +  DDNORM
  200 CONTINUE

*     Use a plane rotation on the right to eliminate the
*     super-diagonal element (THETA) of the upper-bidiagonal matrix.
*     Then use the result to estimate  norm(X).

      DELTA  =   SN2 * RHO
      GAMBAR = - CS2 * RHO
      RHS    =   PHI    - DELTA * Z
      ZBAR   =   RHS    / GAMBAR
      XNORM  =   DSQRT( XXNORM    + ZBAR **2 )
      GAMMA  =   DSQRT( GAMBAR**2 + THETA**2 )
      CS2    =   GAMBAR / GAMMA
      SN2    =   THETA  / GAMMA
      Z      =   RHS    / GAMMA
      XXNORM =   XXNORM + Z**2

*     Test for convergence.
*     First, estimate the norm and condition of the matrix  Abar,
*     and the norms of  rbar  and  Abar(transpose)*rbar.

      ANORM  =   DSQRT( BBNORM )
      ACOND  =   ANORM * DSQRT( DDNORM )
      RES1   =   PHIBAR**2
      RES2   =   RES2  +  PSI**2
      RNORM  =   DSQRT( RES1 + RES2 )
      ARNORM =   ALFA  * DABS( TAU )

*     Now use these norms to estimate certain other quantities,
*     some of which will be small near a solution.

      TEST1  =   RNORM /  BNORM
      TEST2  =   ZERO
      IF (RNORM .GT. ZERO) TEST2 = ARNORM / (ANORM * RNORM)
      TEST3  =   ONE   /  ACOND
      T1     =   TEST1 / (ONE  +  ANORM * XNORM / BNORM)
      RTOL   =   BTOL  +  ATOL *  ANORM * XNORM / BNORM

*     The following tests guard against extremely small values of
*     ATOL, BTOL  or  CTOL.  (The user may have set any or all of
*     the parameters  ATOL, BTOL, CONLIM  to zero.)
*     The effect is equivalent to the normal tests using
*     ATOL = RELPR,  BTOL = RELPR,  CONLIM = 1/RELPR.

      T3     =   ONE + TEST3
      T2     =   ONE + TEST2
      T1     =   ONE + T1
      IF (ITN .GE. ITNLIM) ISTOP = 7
      IF (T3  .LE. ONE   ) ISTOP = 6
      IF (T2  .LE. ONE   ) ISTOP = 5
      IF (T1  .LE. ONE   ) ISTOP = 4

*     Allow for tolerances set by the user.

      IF (TEST3 .LE. CTOL) ISTOP = 3
      IF (TEST2 .LE. ATOL) ISTOP = 2
      IF (TEST1 .LE. RTOL) ISTOP = 1
*     ==================================================================

*     See if it is time to print something.

      IF (NOUT  .LE.  0       ) GO TO 600
      IF (N     .LE. 40       ) GO TO 400
      IF (ITN   .LE. 100       ) GO TO 400
      IF (ITN   .GE. ITNLIM-10) GO TO 400
      IF (MOD(ITN,10) .EQ. 0  ) GO TO 400
      IF (TEST3 .LE.  2.0*CTOL) GO TO 400
      IF (TEST2 .LE. 10.0*ATOL) GO TO 400
      IF (TEST1 .LE. 10.0*RTOL) GO TO 400
      IF (ISTOP .NE.  0       ) GO TO 400
      GO TO 600

*     Print a line for this iteration.

  400 WRITE(NOUT, 1500) ITN, X(1), RNORM, TEST1, TEST2, ANORM, ACOND
      IF (MOD(ITN,10) .EQ. 0) WRITE(NOUT, 1600)
*     ==================================================================

*     Stop if appropriate.
*     The convergence criteria are required to be met on  NCONV
*     consecutive iterations, where  NCONV  is set below.
*     Suggested value:  NCONV = 1, 2  or  3.

  600 IF (ISTOP .EQ. 0) NSTOP = 0
      IF (ISTOP .EQ. 0) GO TO 100
      NCONV  =   1
      NSTOP  =   NSTOP + 1
      IF (NSTOP .LT. NCONV  .AND.  ITN .LT. ITNLIM) ISTOP = 0
      IF (ISTOP .EQ. 0) GO TO 100
*     ------------------------------------------------------------------
*     End of iteration loop.
*     ------------------------------------------------------------------

*     Finish off the standard error estimates.

      T    =   ONE
      IF (M      .GT.   N )  T = M - N
      IF (DAMPSQ .GT. ZERO)  T = M
      T    =   RNORM / DSQRT( T )

      DO 700  I = 1, N
         SE(I)  = T * DSQRT( SE(I) )
  700 CONTINUE

*     Print the stopping condition.

  800 IF (NOUT .GT. 0) THEN
         WRITE(NOUT, 2000) EXIT, ISTOP, ITN,
     :                     EXIT, ANORM, ACOND,
     :                     EXIT, RNORM, ARNORM,
     :                     EXIT, BNORM, XNORM
         WRITE(NOUT, 3000) EXIT, MSG(ISTOP)
      END IF

  900 RETURN

*     ------------------------------------------------------------------
 1000 FORMAT(// 1P, A, '  Least-squares solution of  A*x = b'
     :    / ' The matrix  A  has', I7, ' rows   and', I7, ' columns'
     :    / ' The damping parameter is         DAMP   =', E10.2
     :    / ' ATOL   =', E10.2, 15X,        'CONLIM =', E10.2
     :    / ' BTOL   =', E10.2, 15X,        'ITNLIM =', I10)
 1200 FORMAT(// '   Itn       x(1)           Function',
     :   '     Compatible   LS        Norm A    Cond A' /)
 1300 FORMAT(// '   Itn       x(1)           Function',
     :   '     Compatible   LS     Norm Abar Cond Abar' /)
 1500 FORMAT(1P, I6, 2E17.9, 4E10.2)
 1600 FORMAT(1X)
 2000 FORMAT(/ 1P, A, 6X, 'ISTOP =', I3,   16X, 'ITN    =', I9
     :       /     A, 6X, 'ANORM =', E13.5, 6X, 'ACOND  =', E13.5
     :       /     A, 6X, 'RNORM =', E13.5, 6X, 'ARNORM =', E13.5,
     :       /     A, 6X, 'BNORM =', E13.5, 6X, 'XNORM  =', E13.5)
 3000 FORMAT( A, 6X, A )
*     ------------------------------------------------------------------
*     End of LSQR
      END

*DECK PDA_DSCAL
      SUBROUTINE PDA_DSCAL (N, DA, DX, INCX)
      Implicit NONE
C***BEGIN PROLOGUE  PDA_DSCAL
C***PURPOSE  Multiply a vector by a constant.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A6
C***TYPE      DOUBLE PRECISION (SSCAL-S, PDA_DSCAL-D, CSCAL-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scale factor
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C       DX  double precision result (unchanged if N.LE.0)
C
C     Replace double precision DX by double precision DA*DX.
C     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  PDA_DSCAL
      DOUBLE PRECISION       ::DA, DX(*)
      INTEGER    ::I, INCX, IX, M, MP1, N
C***FIRST EXECUTABLE STATEMENT  PDA_DSCAL
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DX(IX) = DA*DX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GOTO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I+1) = DA*DX(I+1)
        DX(I+2) = DA*DX(I+2)
        DX(I+3) = DA*DX(I+3)
        DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END


*DECK PDA_DCOPY
      SUBROUTINE PDA_DCOPY (N, DX, INCX, DY, INCY)
      Implicit NONE
C***BEGIN PROLOGUE  PDA_DCOPY
C***PURPOSE  Copy a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A5
C***TYPE      DOUBLE PRECISION (SCOPY-S, PDA_DCOPY-D, CCOPY-C, ICOPY-I)
C***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  copy of vector DX (unchanged if N .LE. 0)
C
C     Copy double precision DX to double precision DY.
C     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  PDA_DCOPY
      DOUBLE PRECISION       ::DX(*), DY(*)
      Integer    :: IX, IY, INCX, INCY
      Integer    :: M, N, MP1, I, NS
C***FIRST EXECUTABLE STATEMENT  PDA_DCOPY
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 7.
C
   20 M = MOD(N,7)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF (N .LT. 7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I+1) = DX(I+1)
        DY(I+2) = DX(I+2)
        DY(I+3) = DX(I+3)
        DY(I+4) = DX(I+4)
        DY(I+5) = DX(I+5)
        DY(I+6) = DX(I+6)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DX(I)
   70 CONTINUE
      RETURN
      END
C----------
*DECK PDA_DNRM2
      FUNCTION PDA_DNRM2 (N, DX, INCX)
      Implicit NONE
C***BEGIN PROLOGUE  PDA_DNRM2
C***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A3B
C***TYPE      DOUBLE PRECISION (SNRM2-S, PDA_DNRM2-D, SCNRM2-C)
C***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
C             LINEAR ALGEBRA, UNITARY, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C    PDA_DNRM2  double precision result (zero if N .LE. 0)
C
C     Euclidean norm of the N-vector stored in DX with storage
C     increment INCX.
C     If N .LE. 0, return with result = 0.
C     If N .GE. 1, then INCX must be .GE. 1
C
C     Four phase method using two built-in constants that are
C     hopefully applicable to all machines.
C         CUTLO = maximum of  DSQRT(U/EPS)  over all known machines.
C         CUTHI = minimum of  DSQRT(V)      over all known machines.
C     where
C         EPS = smallest no. such that EPS + 1. .GT. 1.
C         U   = smallest positive no.   (underflow limit)
C         V   = largest  no.            (overflow  limit)
C
C     Brief outline of algorithm.
C
C     Phase 1 scans zero components.
C     move to phase 2 when a component is nonzero and .LE. CUTLO
C     move to phase 3 when a component is .GT. CUTLO
C     move to phase 4 when a component is .GE. CUTHI/M
C     where M = N for X() real and M = 2*N for complex.
C
C     Values for CUTLO and CUTHI.
C     From the environmental parameters listed in the IMSL converter
C     document the limiting values are as follows:
C     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
C                   Univac and DEC at 2**(-103)
C                   Thus CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
C                   Thus CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
C                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
C     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  PDA_DNRM2
      INTEGER    ::N, I, J, NN, INCX
      INTEGER    ::NEXT
      DOUBLE PRECISION       ::PDA_DNRM2 
      DOUBLE PRECISION       ::DX(*), CUTLO, CUTHI, HITEST, SUM, XMAX,
     +                 ZERO, ONE
      SAVE CUTLO, CUTHI, ZERO, ONE
      DATA ZERO, ONE /0.0, 1.0/

      DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
C***FIRST EXECUTABLE STATEMENT  PDA_DNRM2
C      CUTLO = 8.232D-11
C      CUTHI = 1.304D19
C      ZERO = 0.0
C      ONE = 1.0
      IF (N .GT. 0) GO TO 10
         PDA_DNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C
C                                                 BEGIN MAIN LOOP
C
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF (DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF (DX(I) .EQ. ZERO) GO TO 200
      IF (DABS(DX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
C
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF (DABS(DX(I)) .GT. CUTLO) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF (DABS(DX(I)) .LE. XMAX) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI / N
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J = I,NN,INCX
      IF (DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      PDA_DNRM2 = DSQRT(SUM)
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF (I .LE. NN) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      PDA_DNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END

C  The code to calculate the 
*
*              CALL APROD ( mode, m, n, x, y, LENIW, LENJW LENRW, IW, JW, RW )
*
*     which must perform the following functions:
*
*                If MODE = 1, compute  y = y + A*x.
*                If MODE = 2, compute  x = x + A(transpose)*y.
*
*     The vectors x and y are input parameters in both cases.
*     If  mode = 1,  y should be altered without changing x.
*     If  mode = 2,  x should be altered without changing y.
*     The parameters LENIW, LENRW, IW, RW may be used for workspace
*     as described below.

      SUBROUTINE APROD (MODE,M,N,X,Y,LENIW,LENJW,LENRW,IW,JW,RW)
      Implicit NONE
      INTEGER       :: M, N, LENIW,LENJW, LENRW, MODE,I,J
      INTEGER       :: IW(LENIW), JW(LENJW)
      INTEGER       :: RSTART, REND, K
      DOUBLE PRECISION          :: RW(LENRW), X(N), Y(M)
      
C     
C      write(*,*) "leniw=", leniw," lenjw=",lenjw," lenrw =",lenjw
C      write(*,*) IW(1), IW(2), JW(1),JW(2)
      if(MODE .EQ. 1) then 
         DO 111 I = 1, LENIW - 1
            RSTART = IW(I)
            REND = IW(I + 1)
            DO 112 J = RSTART + 1, REND 
               K = JW(J) + 1
C               write(*,*) I, J, K," mode = 1"
               Y(I) = Y(I) + RW(J) * X(K)
 112        CONTINUE 
 111     CONTINUE        
      else if(MODE .eq. 2) then
         DO 113 I = 1, LENIW - 1
            RSTART = IW(I)
            REND = IW(I + 1)
            DO 114 J = RSTART + 1, REND
               K = JW(J) + 1
C               write(*,*) I, J, K, RSTART+1,REND, RW(J)," mode = 2"
               X(K) = X(K) + RW(J) * Y(I)
 114        CONTINUE 
 113     CONTINUE       
      END IF

      RETURN
      END

