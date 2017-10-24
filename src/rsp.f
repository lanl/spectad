c **********************************************************************
c This SOFTWARE  has been authored  by an employee   or employees of the
c University  of  California, operator    of the  Los   Alamos  National
c Laboratory  under Contract No. W-7405-ENG-36  with the U.S. Department
c of  Energy.   The U.S. Government  has  rights to use,  reproduce, and
c distribute this SOFTWARE.   The  public may copy,  prepare  derivative
c works and publicly display this SOFTWARE without charge, provided that
c this  Notice  and any statement  of  authorship are  reproduced on all
c copies.  Neither the Government nor the University makes any warranty,
c express or implied, or assumes any liability or responsibility for the
c use  of this SOFTWARE.  If SOFTWARE  is modified to produce derivative
c works, such  modified SOFTWARE should be clearly  marked, so as not to
c confuse it with the version available from LANL.
c
c It is expressly forbidden to distribute this code, in whole or in part,
c either as is or modified.
c **********************************************************************


  
      SUBROUTINE RSP(NM,N,NV,A,W,MATZ,Z,FV1,FV2,IERR)
C***BEGIN PROLOGUE  RSP
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D4A1
C***KEYWORDS  LIBRARY=SLATEC(EISPACK),TYPE=SINGLE PRECISION(RSP-S),
C             EIGENVALUES,EIGENVECTORS
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Compute eigenvalues and, optionally, eigenvectors of
C            real symmetric matrix packed into a one dimensional
C***DESCRIPTION
C
C     This subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (EISPACK)
C     to find the eigenvalues and eigenvectors (if desired)
C     of a REAL SYMMETRIC PACKED matrix.
C
C     On Input
C
C        NM  must be set to the row dimension of the two-dimensional
C        array parameters as declared in the calling program
C        dimension statement.
C
C        N  is the order of the matrix  A.
C
C        NV  is an integer variable set equal to the
C        dimension of the array  A  as specified for
C        A  in the calling program.  NV  must not be
C        less than  N*(N+1)/2.
C
C        A  contains the lower triangle of the real symmetric
C        packed matrix stored row-wise.
C
C        MATZ  is an integer variable set equal to zero if
C        only eigenvalues are desired.  Otherwise it is set to
C        any non-zero integer for both eigenvalues and eigenvectors.
C
C     On Output
C
C        W  contains the eigenvalues in ascending order.
C
C        Z  contains the eigenvectors if MATZ is not zero.
C
C        IERR  is an integer output variable set equal to an
C        error completion code described in section 2B OF the
C        documentation.  The normal completion code is zero.
C
C        FV1  and  FV2  are temporary storage arrays.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  TQL2,TQLRAT,TRBAK3,TRED3
C***END PROLOGUE  RSP
C
      implicit real*8(a-h,o-z)
      INTEGER I,J,N,NM,NV,IERR,MATZ
      REAL*8 A(NV),W(N),Z(NM,N),FV1(N),FV2(N)

C***FIRST EXECUTABLE STATEMENT  RSP
      IF (N .LE. NM) GO TO 5
      IERR = 10 * N
      GO TO 50
    5 IF (NV .GE. (N * (N + 1)) / 2) GO TO 10
      IERR = 20 * N
      GO TO 50

   10 CALL  TRED3(N,NV,A,W,FV1,FV2)
      IF (MATZ .gt. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      if(matz.eq.0) then
        CALL  TQLRAT(N,W,FV2,IERR)
      else if(matz.eq.-1) then
        CALL  TQL2x(NM,N,W,FV1,IERR)
      end if

      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 DO 40 I = 1, N

         DO 30 J = 1, N
            Z(J,I) = 0.0d0
   30    CONTINUE

         Z(I,I) = 1.0d0
   40 CONTINUE

      CALL  TQL2(NM,N,W,FV1,Z,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  TRBAK3(NM,N,NV,A,N,Z)
   50 RETURN
      END
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C***BEGIN PROLOGUE  TQL2
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D4A5,D4C2A
C***KEYWORDS  LIBRARY=SLATEC(EISPACK),TYPE=SINGLE PRECISION(TQL2-S),
C             EIGENVALUES,EIGENVECTORS
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Compute eigenvalues and eigenvectors of symmetric
C            tridiagonal matrix.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TQL2,
C     NUM. MATH. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
C     Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     This subroutine finds the eigenvalues and eigenvectors
C     of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
C     The eigenvectors of a FULL SYMMETRIC matrix can also
C     be found if  TRED2  has been used to reduce this
C     full matrix to tridiagonal form.
C
C     On Input
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        D contains the diagonal elements of the input matrix.
C
C        E contains the subdiagonal elements of the input matrix
C          in its last N-1 positions.  E(1) is arbitrary.
C
C        Z contains the transformation matrix produced in the
C          reduction by  TRED2, if performed.  If the eigenvectors
C          of the tridiagonal matrix are desired, Z must contain
C          the identity matrix.
C
C      On Output
C
C        D contains the eigenvalues in ascending order.  If an
C          error exit is made, the eigenvalues are correct but
C          unordered for indices 1,2,...,IERR-1.
C
C        E has been destroyed.
C
C        Z contains orthonormal eigenvectors of the symmetric
C          tridiagonal (or full) matrix.  If an error exit is made,
C          Z contains the eigenvectors associated with the stored
C          eigenvalues.
C
C        IERR is set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C
C     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  PYTHAG
C***END PROLOGUE  TQL2
C
      implicit real*8(a-h,o-z)
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      REAL*8 D(N),E(N),Z(NM,N)
      REAL*8 B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
      REAL*8 PYTHAG
C
C***FIRST EXECUTABLE STATEMENT  TQL2
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0d0
      B = 0.0d0
      E(N) = 0.0d0
C
      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (B .LT. H) B = H
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            IF (B + ABS(E(M)) .EQ. B) GO TO 120
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0d0 * E(L))
         R = PYTHAG(P,1.0d0)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
C
         DO 140 I = L2, N
  140    D(I) = D(I) - H
C
  145    F = F + H
C     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0d0
         C2 = C
         EL1 = E(L1)
         S = 0.0d0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0d0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0d0 / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0d0)
            E(I+1) = S * E(I) * R
            S = 1.0d0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
C
  200    CONTINUE
C
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         IF (B + ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      SUBROUTINE TQLRAT(N,D,E2,IERR)
C***BEGIN PROLOGUE  TQLRAT
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D4A5,D4C2A
C***KEYWORDS  LIBRARY=SLATEC(EISPACK),TYPE=SINGLE PRECISION(TQLRAT-S),
C             EIGENVALUES,EIGENVECTORS
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Computes eigenvalues of symmetric tridiagonal matrix
C            a rational variant of the QL method.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TQLRAT,
C     ALGORITHM 464, COMM. ACM 16, 689(1973) by Reinsch.
C
C     This subroutine finds the eigenvalues of a SYMMETRIC
C     TRIDIAGONAL matrix by the rational QL method.
C
C     On Input
C
C        N is the order of the matrix.
C
C        D contains the diagonal elements of the input matrix.
C
C        E2 contains the squares of the subdiagonal elements of the
C          input matrix in its last N-1 positions.  E2(1) is arbitrary.
C
C      On Output
C
C        D contains the eigenvalues in ascending order.  If an
C          error exit is made, the eigenvalues are correct and
C          ordered for indices 1,2,...IERR-1, but may not be
C          the smallest eigenvalues.
C
C        E2 has been destroyed.
C
C        IERR is set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C
C     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  PYTHAG
C***END PROLOGUE  TQLRAT
C
      implicit real*8(a-h,o-z)
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      REAL*8 D(N),E2(N)
      REAL*8 B,C,F,G,H,P,R,S,MACHEP
      REAL*8 PYTHAG
C
      SAVE MACHEP
      DATA MACHEP/1.0d0/
C***FIRST EXECUTABLE STATEMENT  TQLRAT
      IF (MACHEP .NE. 1.0d0) GO TO 10
   05 MACHEP = 0.5d0*MACHEP
      IF (1.0d0 + MACHEP .GT. 1.0E0) GO TO 05
      MACHEP = 2.0d0*MACHEP
C
   10 IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C
      F = 0.0d0
      B = 0.0d0
      E2(N) = 0.0d0
C
      DO 290 L = 1, N
         J = 0
         H = MACHEP * (ABS(D(L)) + SQRT(E2(L)))
         IF (B .GT. H) GO TO 105
         B = H
         C = B * B
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         S = SQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0d0 * S)
         R = PYTHAG(P,1.0d0)
         D(L) = S / (P + SIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     .......... RATIONAL QL TRANSFORMATION ..........
         G = D(M)
         IF (G .EQ. 0.0d0) G = B
         H = G
         S = 0.0d0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. 0.0d0) G = B
            H = G * P / R
  200    CONTINUE
C
         E2(L) = S * G
         D(L) = H
C     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
         IF (H .EQ. 0.0d0) GO TO 210
         IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0d0) GO TO 130
  210    P = D(L) + F
C     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      SUBROUTINE TRBAK3(NM,N,NV,A,M,Z)
C***BEGIN PROLOGUE  TRBAK3
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D4C4
C***KEYWORDS  LIBRARY=SLATEC(EISPACK),TYPE=SINGLE PRECISION(TRBAK3-S),
C             EIGENVALUES,EIGENVECTORS
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Forms eigenvectors of real symmetric matrix from the
C            eigenvectors of symmetric tridiagonal matrix formed
C            TRED3.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TRBAK3,
C     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     This subroutine forms the eigenvectors of a REAL SYMMETRIC
C     matrix by back transforming those of the corresponding
C     symmetric tridiagonal matrix determined by  TRED3.
C
C     On Input
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        NV must be set to the dimension of the array parameter A
C          as declared in the calling program dimension statement.
C
C        A contains information about the orthogonal transformations
C          used in the reduction by  TRED3  in its first
C          N*(N+1)/2 positions.
C
C        M is the number of eigenvectors to be back transformed.
C
C        Z contains the eigenvectors to be back transformed
C          in its first M columns.
C
C     On Output
C
C        Z contains the transformed eigenvectors
C          in its first M columns.
C
C     Note that TRBAK3 preserves vector Euclidean norms.
C
C     Questions and comments should be directed to b. s. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  TRBAK3
C
      implicit real*8(a-h,o-z)
      INTEGER I,J,K,L,M,N,IK,IZ,NM,NV
      REAL*8 A(NV),Z(NM,M)
      REAL*8 H,S
C
C***FIRST EXECUTABLE STATEMENT  TRBAK3
      IF (M .EQ. 0) GO TO 200
      IF (N .EQ. 1) GO TO 200
C
      DO 140 I = 2, N
         L = I - 1
         IZ = (I * L) / 2
         IK = IZ + I
         H = A(IK)
         IF (H .EQ. 0.0d0) GO TO 140
C
         DO 130 J = 1, M
            S = 0.0d0
            IK = IZ
C
            DO 110 K = 1, L
               IK = IK + 1
               S = S + A(IK) * Z(K,J)
  110       CONTINUE
C     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            S = (S / H) / H
            IK = IZ
C
            DO 120 K = 1, L
               IK = IK + 1
               Z(K,J) = Z(K,J) - S * A(IK)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE TRED3(N,NV,A,D,E,E2)
C***BEGIN PROLOGUE  TRED3
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D4C1B1
C***KEYWORDS  LIBRARY=SLATEC(EISPACK),TYPE=SINGLE PRECISION(TRED3-S),
C             EIGENVALUES,EIGENVECTORS
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Reduce real symmetric matrix stored in packed form to
C            symmetric tridiagonal matrix using orthogonal
C            transformations.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TRED3,
C     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     This subroutine reduces a REAL SYMMETRIC matrix, stored as
C     a one-dimensional array, to a symmetric tridiagonal matrix
C     using orthogonal similarity transformations.
C
C     On Input
C
C        n is the order of the matrix.
C
C        NV must be set to the dimension of the array parameter A
C          as declared in the calling program dimension statement.
C
C        A contains the lower triangle of the real symmetric
C          input matrix, stored row-wise as a one-dimensional
C          array, in its first N*(N+1)/2 positions.
C
C     On Output
C
C        A contains information about the orthogonal
C          transformations used in the reduction.
C
C        D contains the diagonal elements of the tridiagonal matrix.
C
C        E contains the subdiagonal elements of the tridiagonal
C          matrix in its last N-1 positions.  E(1) is set to zero.
C
C        E2 contains the squares of the corresponding elements of E.
C          E2 may coincide with E if the squares are not needed.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  TRED3
C
      implicit real*8(a-h,o-z)
      INTEGER I,J,K,L,N,II,IZ,JK,NV
      REAL*8 A(NV),D(N),E(N),E2(N)
      REAL*8 F,G,H,HH,SCALE
C
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
C***FIRST EXECUTABLE STATEMENT  TRED3
      DO  300 II = 1, N
         I = N + 1 - II
         L = I - 1
         IZ = (I * L) / 2
         H = 0.0d0
         SCALE = 0.0d0
         IF (L .LT. 1) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
            IZ = IZ + 1
            D(K) = A(IZ)
            SCALE = SCALE + ABS(D(K))
  120    CONTINUE
C
         IF (SCALE .NE. 0.0d0) GO TO 140
  130    E(I) = 0.0d0
         E2(I) = 0.0d0
         GO TO 290
C
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         F = D(L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         A(IZ) = SCALE * D(L)
         IF (L .EQ. 1) GO TO 290
         F = 0.0d0
C
         DO 240 J = 1, L
            G = 0.0d0
            JK = (J * (J-1)) / 2
C     .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, L
               JK = JK + 1
               IF (K .GT. J) JK = JK + K - 2
               G = G + A(JK) * D(K)
  180       CONTINUE
C     .......... FORM ELEMENT OF P ..........
            E(J) = G / H
            F = F + E(J) * D(J)
  240    CONTINUE
C
         HH = F / (H + H)
         JK = 0
C     .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = D(J)
            G = E(J) - HH * F
            E(J) = G
C
            DO 260 K = 1, J
               JK = JK + 1
               A(JK) = A(JK) - F * E(K) - G * D(K)
  260    CONTINUE
C
  290    D(I) = A(IZ+1)
         A(IZ+1) = SCALE * SQRT(H)
  300 CONTINUE
C
      RETURN
      END
      REAL*8 FUNCTION PYTHAG(A,B)
C***BEGIN PROLOGUE  PYTHAG
C***REFER TO  EISDOC
C
C     Finds sqrt(A**2+B**2) without overflow or destructive underflow
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  PYTHAG
      implicit real*8(a-h,o-z)
      REAL*8 A,B
C
      REAL*8 P,Q,R,S,T
C***FIRST EXECUTABLE STATEMENT  PYTHAG
      P = MAX(ABS(A),ABS(B))
      Q = MIN(ABS(A),ABS(B))
      IF (Q .EQ. 0.0d0) GO TO 20
   10 CONTINUE
         R = (Q/P)**2
         T = 4.0d0 + R
         IF (T .EQ. 4.0d0) GO TO 20
         S = R/T
         P = P + 2.0d0*P*S
         Q = Q*S
      GO TO 10
   20 PYTHAG = P
      RETURN
      END

cc**** beginning here are the nonredundant subroutines necessary for RS ****

      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)                       rs     3
C***BEGIN PROLOGUE  RS                                                  rs     4
C***DATE WRITTEN   760101   (YYMMDD)                                    rs     5
C***REVISION DATE  861211   (YYMMDD)                                    rs     6
C***CATEGORY NO.  D4A1                                                  rs     7
C***KEYWORDS  LIBRARY=SLATEC(EISPACK),TYPE=SINGLE PRECISION(RS-S CH-C), rs     8
C             EIGENVALUES,EIGENVECTORS                                  rs     9
C***AUTHOR  SMITH, B. T., ET AL.                                        rs    10
C***PURPOSE  Computes eigenvalues and, optionally, eigenvectors of      rs    11
C            real symmetric matrix.                                     rs    12
C***DESCRIPTION                                                         rs    13
C                                                                       rs    14
C     This subroutine calls the recommended sequence of                 rs    15
C     subroutines from the eigensystem subroutine package (EISPACK)     rs    16
C     to find the eigenvalues and eigenvectors (if desired)             rs    17
C     of a REAL SYMMETRIC matrix.                                       rs    18
C                                                                       rs    19
C     On Input                                                          rs    20
C                                                                       rs    21
C        NM  must be set to the row dimension of the two-dimensional    rs    22
C        array parameters as declared in the calling program            rs    23
C        dimension statement.                                           rs    24
C                                                                       rs    25
C        N  is the order of the matrix  A.                              rs    26
C                                                                       rs    27
C        A  contains the real symmetric matrix.                         rs    28
C                                                                       rs    29
C        MATZ  is an integer variable set equal to zero if              rs    30
C        only eigenvalues are desired.  Otherwise it is set to          rs    31
C        any non-zero integer for both eigenvalues and eigenvectors.    rs    32
C                                                                       rs    33
C     On Output                                                         rs    34
C                                                                       rs    35
C        W  contains the eigenvalues in ascending order.                rs    36
C                                                                       rs    37
C        Z  contains the eigenvectors if MATZ is not zero.              rs    38
C                                                                       rs    39
C        IERR  is an integer output variable set equal to an            rs    40
C        error completion code described in section 2B of the           rs    41
C        documentation.  The normal completion code is zero.            rs    42
C                                                                       rs    43
C        FV1  and  FV2  are temporary storage arrays.                   rs    44
C                                                                       rs    45
C     Questions and comments should be directed to B. S. Garbow,        rs    46
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         rs    47
C     ------------------------------------------------------------------rs    48
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW, rs    49
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-    rs    50
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,    rs    51
C                 1976.                                                 rs    52
C***ROUTINES CALLED  TQL2,TQLRAT,TRED1,TRED2                            rs    53
C***END PROLOGUE  RS                                                    rs    54
C     implicit real*8(a-h,o-z)                                                                  rs    55
      INTEGER N,NM,IERR,MATZ                                            rs    56
      REAL*8 A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)                         rs    57
C                                                                       rs    58
C***FIRST EXECUTABLE STATEMENT  RS                                      rs    59
      IF (N .LE. NM) GO TO 10                                           rs    60
      IERR = 10 * N                                                     rs    61
      GO TO 50                                                          rs    62
C                                                                       rs    63
   10 IF (MATZ .NE. 0) GO TO 20                                         rs    64
C     .......... FIND EIGENVALUES ONLY ..........                       rs    65
      CALL  TRED1(NM,N,A,W,FV1,FV2)                                     rs    66
      CALL  TQLRAT(N,W,FV2,IERR)                                        rs    67
      GO TO 50                                                          rs    68
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........      rs    69
   20 CALL  TRED2(NM,N,A,W,FV1,Z)                                       rs    70
      CALL  TQL2(NM,N,W,FV1,Z,IERR)                                     rs    71
   50 RETURN                                                            rs    72
      END                                                               rs    73
      SUBROUTINE TRED1(NM,N,A,D,E,E2)                                   tred1  3
C***BEGIN PROLOGUE  TRED1                                               tred1  4
C***DATE WRITTEN   760101   (YYMMDD)                                    tred1  5
C***REVISION DATE  861211   (YYMMDD)                                    tred1  6
C***CATEGORY NO.  D4C1B1                                                tred1  7
C***KEYWORDS  LIBRARY=SLATEC(EISPACK),TYPE=SINGLE PRECISION(TRED1-S),   tred1  8
C             EIGENVALUES,EIGENVECTORS                                  tred1  9
C***AUTHOR  SMITH, B. T., ET AL.                                        tred1 10
C***PURPOSE  Reduce real symmetric matrix to symmetric tridiagonal      tred1 11
C            matrix using orthogonal similarity transformations.        tred1 12
C***DESCRIPTION                                                         tred1 13
C                                                                       tred1 14
C     This subroutine is a translation of the ALGOL procedure TRED1,    tred1 15
C     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.   tred1 16
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   tred1 17
C                                                                       tred1 18
C     This subroutine reduces a REAL SYMMETRIC matrix                   tred1 19
C     to a symmetric tridiagonal matrix using                           tred1 20
C     orthogonal similarity transformations.                            tred1 21
C                                                                       tred1 22
C     On Input                                                          tred1 23
C                                                                       tred1 24
C        NM must be set to the row dimension of two-dimensional         tred1 25
C          array parameters as declared in the calling program          tred1 26
C          dimension statement.                                         tred1 27
C                                                                       tred1 28
C        N is the order of the matrix.                                  tred1 29
C                                                                       tred1 30
C        A contains the real symmetric input matrix.  Only the          tred1 31
C          lower triangle of the matrix need be supplied.               tred1 32
C                                                                       tred1 33
C     On Output                                                         tred1 34
C                                                                       tred1 35
C        A contains information about the orthogonal trans-             tred1 36
C          formations used in the reduction in its strict lower         tred1 37
C          triangle.  The full upper triangle of A is unaltered.        tred1 38
C                                                                       tred1 39
C        D contains the diagonal elements of the tridiagonal matrix.    tred1 40
C                                                                       tred1 41
C        E contains the subdiagonal elements of the tridiagonal         tred1 42
C          matrix in its last N-1 positions.  E(1) is set to zero.      tred1 43
C                                                                       tred1 44
C        E2 contains the squares of the corresponding elements of E.    tred1 45
C          E2 may coincide with E if the squares are not needed.        tred1 46
C                                                                       tred1 47
C     Questions and comments should be directed to B. S. Garbow,        tred1 48
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         tred1 49
C     ------------------------------------------------------------------tred1 50
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW, tred1 51
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-    tred1 52
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,    tred1 53
C                 1976.                                                 tred1 54
C***ROUTINES CALLED  (NONE)                                             tred1 55
C***END PROLOGUE  TRED1                                                 tred1 56
C                                                                       tred1 57
      implicit real*8 (a-h,o-z)
      INTEGER I,J,K,L,N,II,NM,JP1                                       tred1 58
      REAL*8 A(NM,N),D(N),E(N),E2(N)                                      tred1 59
      REAL*8 F,G,H,SCALE                                                  tred1 60
C                                                                       tred1 61
C***FIRST EXECUTABLE STATEMENT  TRED1                                   tred1 62
      write(6,*) 'hello from tred1'
      DO 100 I = 1, N                                                   tred1 63
  100 D(I) = A(I,I)                                                     tred1 64
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........               tred1 65
      DO 300 II = 1, N                                                  tred1 66
         I = N + 1 - II                                                 tred1 67
         L = I - 1                                                      tred1 68
         H = 0.0D0                                                      tred1 69
         SCALE = 0.0D0                                                  tred1 70
         IF (L .LT. 1) GO TO 130                                        tred1 71
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........       tred1 72
         DO 120 K = 1, L                                                tred1 73
  120    SCALE = SCALE + ABS(A(I,K))                                    tred1 74
C                                                                       tred1 75
         IF (SCALE .NE. 0.0D0) GO TO 140                                tred1 76
  130    E(I) = 0.0D0                                                   tred1 77
         E2(I) = 0.0D0                                                  tred1 78
         GO TO 290                                                      tred1 79
C                                                                       tred1 80
  140    DO 150 K = 1, L                                                tred1 81
            A(I,K) = A(I,K) / SCALE                                     tred1 82
            H = H + A(I,K) * A(I,K)                                     tred1 83
  150    CONTINUE                                                       tred1 84
C                                                                       tred1 85
         E2(I) = SCALE * SCALE * H                                      tred1 86
         F = A(I,L)                                                     tred1 87
         G = -SIGN(SQRT(H),F)                                           tred1 88
         E(I) = SCALE * G                                               tred1 89
         H = H - F * G                                                  tred1 90
         A(I,L) = F - G                                                 tred1 91
         IF (L .EQ. 1) GO TO 270                                        tred1 92
         F = 0.0D0                                                      tred1 93
C                                                                       tred1 94
         DO 240 J = 1, L                                                tred1 95
            G = 0.0D0                                                   tred1 96
C     .......... FORM ELEMENT OF A*U ..........                         tred1 97
            DO 180 K = 1, J                                             tred1 98
  180       G = G + A(J,K) * A(I,K)                                     tred1 99
C                                                                       tred1100
            JP1 = J + 1                                                 tred1101
            IF (L .LT. JP1) GO TO 220                                   tred1102
C                                                                       tred1103
            DO 200 K = JP1, L                                           tred1104
  200       G = G + A(K,J) * A(I,K)                                     tred1105
C     .......... FORM ELEMENT OF P ..........                           tred1106
  220       E(J) = G / H                                                tred1107
            F = F + E(J) * A(I,J)                                       tred1108
  240    CONTINUE                                                       tred1109
C                                                                       tred1110
         H = F / (H + H)                                                tred1111
C     .......... FORM REDUCED A ..........                              tred1112
         DO 260 J = 1, L                                                tred1113
            F = A(I,J)                                                  tred1114
            G = E(J) - H * F                                            tred1115
            E(J) = G                                                    tred1116
C                                                                       tred1117
            DO 260 K = 1, J                                             tred1118
               A(J,K) = A(J,K) - F * E(K) - G * A(I,K)                  tred1119
  260    CONTINUE                                                       tred1120
C                                                                       tred1121
  270    DO 280 K = 1, L                                                tred1122
  280    A(I,K) = SCALE * A(I,K)                                        tred1123
C                                                                       tred1124
  290    H = D(I)                                                       tred1125
         D(I) = A(I,I)                                                  tred1126
         A(I,I) = H                                                     tred1127
  300 CONTINUE                                                          tred1128
C                                                                       tred1129
      RETURN                                                            tred1130
      END                                                               tred1131
      SUBROUTINE TRED2(NM,N,A,D,E,Z)                                    tred2  3
C***BEGIN PROLOGUE  TRED2                                               tred2  4
C***DATE WRITTEN   760101   (YYMMDD)                                    tred2  5
C***REVISION DATE  861211   (YYMMDD)                                    tred2  6
C***CATEGORY NO.  D4C1B1                                                tred2  7
C***KEYWORDS  LIBRARY=SLATEC(EISPACK),TYPE=SINGLE PRECISION(TRED2-S),   tred2  8
C             EIGENVALUES,EIGENVECTORS                                  tred2  9
C***AUTHOR  SMITH, B. T., ET AL.                                        tred2 10
C***PURPOSE  Reduce real symmetric matrix to symmetric tridiagonal      tred2 11
C            matrix using and accumulating orthogonal transformation    tred2 12
C***DESCRIPTION                                                         tred2 13
C                                                                       tred2 14
C     This subroutine is a translation of the ALGOL procedure TRED2,    tred2 15
C     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.   tred2 16
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   tred2 17
C                                                                       tred2 18
C     This subroutine reduces a REAL SYMMETRIC matrix to a              tred2 19
C     symmetric tridiagonal matrix using and accumulating               tred2 20
C     orthogonal similarity transformations.                            tred2 21
C                                                                       tred2 22
C     On Input                                                          tred2 23
C                                                                       tred2 24
C        NM must be set to the row dimension of two-dimensional         tred2 25
C          array parameters as declared in the calling program          tred2 26
C          dimension statement.                                         tred2 27
C                                                                       tred2 28
C        N is the order of the matrix.                                  tred2 29
C                                                                       tred2 30
C        A contains the real symmetric input matrix.  Only the          tred2 31
C          lower triangle of the matrix need be supplied.               tred2 32
C                                                                       tred2 33
C     On Output                                                         tred2 34
C                                                                       tred2 35
C        D contains the diagonal elements of the tridiagonal matrix.    tred2 36
C                                                                       tred2 37
C        E contains the subdiagonal elements of the tridiagonal         tred2 38
C          matrix in its last N-1 positions.  E(1) is set to zero.      tred2 39
C                                                                       tred2 40
C        Z contains the orthogonal transformation matrix                tred2 41
C          produced in the reduction.                                   tred2 42
C                                                                       tred2 43
C        A and Z may coincide.  If distinct, A is unaltered.            tred2 44
C                                                                       tred2 45
C     Questions and comments should be directed to B. S. Garbow,        tred2 46
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         tred2 47
C     ------------------------------------------------------------------tred2 48
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW, tred2 49
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-    tred2 50
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,    tred2 51
C                 1976.                                                 tred2 52
C***ROUTINES CALLED  (NONE)                                             tred2 53
C***END PROLOGUE  TRED2                                                 tred2 54
C                                                                       tred2 55
      implicit real*8(a-h,o-z)
      INTEGER I,J,K,L,N,II,NM,JP1                                       tred2 56
      REAL*8 A(NM,N),D(N),E(N),Z(NM,N)                                  tred2 57
      REAL*8 F,G,H,HH,SCALE                                             tred2 58
C                                                                       tred2 59
C***FIRST EXECUTABLE STATEMENT  TRED2                                   tred2 60
      DO 100 I = 1, N                                                   tred2 61
C                                                                       tred2 62
         DO 100 J = 1, I                                                tred2 63
            Z(I,J) = A(I,J)                                             tred2 64
  100 CONTINUE                                                          tred2 65
C                                                                       tred2 66
      IF (N .EQ. 1) GO TO 320                                           tred2 67
C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........               tred2 68
      DO 300 II = 2, N                                                  tred2 69
         I = N + 2 - II                                                 tred2 70
         L = I - 1                                                      tred2 71
         H = 0.0D0                                                      tred2 72
         SCALE = 0.0D0                                                  tred2 73
         IF (L .LT. 2) GO TO 130                                        tred2 74
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........       tred2 75
         DO 120 K = 1, L                                                tred2 76
  120    SCALE = SCALE + ABS(Z(I,K))                                    tred2 77
C                                                                       tred2 78
         IF (SCALE .NE. 0.0D0) GO TO 140                                tred2 79
  130    E(I) = Z(I,L)                                                  tred2 80
         GO TO 290                                                      tred2 81
C                                                                       tred2 82
  140    DO 150 K = 1, L                                                tred2 83
            Z(I,K) = Z(I,K) / SCALE                                     tred2 84
            H = H + Z(I,K) * Z(I,K)                                     tred2 85
  150    CONTINUE                                                       tred2 86
C                                                                       tred2 87
         F = Z(I,L)                                                     tred2 88
         G = -SIGN(SQRT(H),F)                                           tred2 89
         E(I) = SCALE * G                                               tred2 90
         H = H - F * G                                                  tred2 91
         Z(I,L) = F - G                                                 tred2 92
         F = 0.0D0                                                      tred2 93
C                                                                       tred2 94
         DO 240 J = 1, L                                                tred2 95
            Z(J,I) = Z(I,J) / H                                         tred2 96
            G = 0.0D0                                                   tred2 97
C     .......... FORM ELEMENT OF A*U ..........                         tred2 98
            DO 180 K = 1, J                                             tred2 99
  180       G = G + Z(J,K) * Z(I,K)                                     tred2100
C                                                                       tred2101
            JP1 = J + 1                                                 tred2102
            IF (L .LT. JP1) GO TO 220                                   tred2103
C                                                                       tred2104
            DO 200 K = JP1, L                                           tred2105
  200       G = G + Z(K,J) * Z(I,K)                                     tred2106
C     .......... FORM ELEMENT OF P ..........                           tred2107
  220       E(J) = G / H                                                tred2108
            F = F + E(J) * Z(I,J)                                       tred2109
  240    CONTINUE                                                       tred2110
C                                                                       tred2111
         HH = F / (H + H)                                               tred2112
C     .......... FORM REDUCED A ..........                              tred2113
         DO 260 J = 1, L                                                tred2114
            F = Z(I,J)                                                  tred2115
            G = E(J) - HH * F                                           tred2116
            E(J) = G                                                    tred2117
C                                                                       tred2118
            DO 260 K = 1, J                                             tred2119
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)                  tred2120
  260    CONTINUE                                                       tred2121
C                                                                       tred2122
  290    D(I) = H                                                       tred2123
  300 CONTINUE                                                          tred2124
C                                                                       tred2125
  320 D(1) = 0.0D0                                                      tred2126
      E(1) = 0.0D0                                                      tred2127
C     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........     tred2128
      DO 500 I = 1, N                                                   tred2129
         L = I - 1                                                      tred2130
         IF (D(I) .EQ. 0.0D0) GO TO 380                                 tred2131
C                                                                       tred2132
         DO 360 J = 1, L                                                tred2133
            G = 0.0D0                                                   tred2134
C                                                                       tred2135
            DO 340 K = 1, L                                             tred2136
  340       G = G + Z(I,K) * Z(K,J)                                     tred2137
C                                                                       tred2138
            DO 360 K = 1, L                                             tred2139
               Z(K,J) = Z(K,J) - G * Z(K,I)                             tred2140
  360    CONTINUE                                                       tred2141
C                                                                       tred2142
  380    D(I) = Z(I,I)                                                  tred2143
         Z(I,I) = 1.0D0                                                 tred2144
         IF (L .LT. 1) GO TO 500                                        tred2145
C                                                                       tred2146
         DO 400 J = 1, L                                                tred2147
            Z(I,J) = 0.0D0                                              tred2148
            Z(J,I) = 0.0D0                                              tred2149
  400    CONTINUE                                                       tred2150
C                                                                       tred2151
  500 CONTINUE                                                          tred2152
C                                                                       tred2153
      RETURN                                                            tred2154
      END                                                               tred2155

      SUBROUTINE TQL2x(NM,N,D,E,IERR)
c  this is tql2 with the Z stuff pulled out -- i.e., it finds values, but not vectors
C***BEGIN PROLOGUE  TQL2
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D4A5,D4C2A
C***KEYWORDS  LIBRARY=SLATEC(EISPACK),TYPE=SINGLE PRECISION(TQL2-S),
C             EIGENVALUES,EIGENVECTORS
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Compute eigenvalues and eigenvectors of symmetric
C            tridiagonal matrix.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TQL2,
C     NUM. MATH. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
C     Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     This subroutine finds the eigenvalues and eigenvectors
C     of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
C     The eigenvectors of a FULL SYMMETRIC matrix can also
C     be found if  TRED2  has been used to reduce this
C     full matrix to tridiagonal form.
C
C     On Input
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        D contains the diagonal elements of the input matrix.
C
C        E contains the subdiagonal elements of the input matrix
C          in its last N-1 positions.  E(1) is arbitrary.
C
c  pulled:
C        Z contains the transformation matrix produced in the
C          reduction by  TRED2, if performed.  If the eigenvectors
C          of the tridiagonal matrix are desired, Z must contain
C          the identity matrix.
C
C      On Output
C
C        D contains the eigenvalues in ascending order.  If an
C          error exit is made, the eigenvalues are correct but
C          unordered for indices 1,2,...,IERR-1.
C
C        E has been destroyed.
C
c  pulled:
C        Z contains orthonormal eigenvectors of the symmetric
C          tridiagonal (or full) matrix.  If an error exit is made,
C          Z contains the eigenvectors associated with the stored
C          eigenvalues.
C
C        IERR is set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C
C     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  PYTHAG
C***END PROLOGUE  TQL2
C
      implicit real*8(a-h,o-z)
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      REAL*8 D(N),E(N)
      REAL*8 B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
      REAL*8 PYTHAG
C
C***FIRST EXECUTABLE STATEMENT  TQL2
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0d0
      B = 0.0d0
      E(N) = 0.0d0
C
      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (B .LT. H) B = H
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            IF (B + ABS(E(M)) .EQ. B) GO TO 120
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0d0 * E(L))
         R = PYTHAG(P,1.0d0)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
C
         DO 140 I = L2, N
  140    D(I) = D(I) - H
C
  145    F = F + H
C     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0d0
         C2 = C
         EL1 = E(L1)
         S = 0.0d0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0d0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0d0 / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0d0)
            E(I+1) = S * E(I) * R
            S = 1.0d0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
ccC     .......... FORM VECTOR ..........
cc            DO 180 K = 1, N
cc               H = Z(K,I+1)
cc               Z(K,I+1) = S * Z(K,I) + C * H
cc               Z(K,I) = C * Z(K,I) - S * H
cc  180       CONTINUE
C
  200    CONTINUE
C
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         IF (B + ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
cc         DO 280 J = 1, N
cc            P = Z(J,I)
cc            Z(J,I) = Z(J,K)
cc            Z(J,K) = P
cc  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
