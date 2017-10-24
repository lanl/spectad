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


  
c RCS keywords for sdiag.f:
c     $Id: sdiag.f,v 2.0 1998/12/18 21:31:04 afv Exp $
c     $Revision: 2.0 $
c     $Author: afv $
c     $Locker:  $
c     $Date: 1998/12/18 21:31:04 $
c     $Header: /n/home3/afv/source/RCS/sdiag.f,v 2.0 1998/12/18 21:31:04 afv Exp $
  
      subroutine sdiag(n,a,eigvec,eigval,work)
c  simply drives rsp routine or eisgiv routine
      implicit real*8(a-h,o-z)
      dimension a(1),eigvec(n,n),eigval(n),work(n)
      common/sdicom/work2(8000),iwork(1000)

      if(n.gt.8000) stop 'abort - 8000 work space exceeded in sdiag/rsp'

c  following turns on eigenvectors (zero for eigenvalues only)
      imat=1
      lena=n*(n+1)/2
      itest=1
      if(itest.eq.1) then
ccc        write(6,*) 'new eisgiv test'
        if(8*n.gt.8000) stop 'abort - workspace exceeded in sdiag/rsp'
        call eisgiv(6,n,n,lena,n,a,work2,iwork,eigval,eigvec,ierr)
        if(ierr.ne.0) stop 'abort - sdiag/eisgiv failed'
      else
        if(n.gt.8000) stop 'abort - workspace exceeded in sdiag/rsp'
        call rsp(n,n,lena,a,eigval,imat,eigvec,work,work2,ierr)
        if(ierr.ne.0) stop 'abort - sdiag/rsp failed'
      end if
      return
      end

      subroutine sdiagf(n,a,eigval,work)
c  just like sdiag, but does not find eigenvectors (saving space and time)
c  simply drives rsp routine
      implicit real*8(a-h,o-z)
      dimension a(1),eigval(n),work(n)
      dimension dum(1)
      common/sdicom/work2(8000),iwork(1000)

      if(n.gt.8000)stop 'abort - 8000 work space exceeded in sdiagf/rsp'

c  turn off eigenvectors
c  imat=-1 is to try out using tql2 instead of tqlrat for eigvecs only
c  tqlrat seems to have a bug at +OP2 compile option on HP735
      imat=-1
      lena=n*(n+1)/2
      call rsp(n,n,lena,a,eigval, imat,dum, work,work2,ierr)
      if(ierr.ne.0) stop 'abort - sdiagf/rsp failed'
      return
      end

c*************************************
      subroutine etrbk3(nm,n,nv,a,m,z)
      implicit real*8(a-h,o-z)
c*************************************
c
      double precision a,h,s,z,sdot,zero
c
      dimension a(nv),z(nm,m)
c     level 2,a,z
c
      data zero/0.0d0/
csng  data zero/0.0e0/
c
c     ------------------------------------------------------------------
c
c     this subroutine is a modification of eispack subroutine trbak3
c     which is a translation of the algol procedure trbak3,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine forms the eigenvectors of a real symmetric
c     matrix by back transforming those of the corresponding
c     symmetric tridiagonal matrix determined by  etred3.
c
c     the calculation is carried out by forming the matrix product qz,
c     where  q  is a product of the orthogonal symmetric matrices
c     encoded in  a  and  z  is the set of eigenvectors of the tri-
c     diagonal matrix  f  which was formed from the original symmetric
c     matrix  c  by the similarity transformation
c              f = q(transpose) c q
c
c
c     on input-
c
c        nm must be set to the row dimension of two-dimensional
c           array parameters as declared in the calling program
c           dimension statement.
c
c        n  is the order of the matrix.
c
c        nv must be set to the dimension of the array parameter a
c           as declared in the calling program dimension statement.
c
c        a  contains information about the orthogonal transformations
c           used in the reduction by  etred3 in its first
c           n*(n+1)/2 positions.
c
c        m  is the number of eigenvectors to be back transformed.
c
c        z  contains the eigenvectors to be back transformed
c           in its first m columns.
c
c     on output-
c
c        z  contains the transformed eigenvectors
c           in its first m columns.
c
c     note that etrbk3 preserves vector euclidean norms.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
      if (m .eq. 0) go to 200
      if (n .eq. 1) go to 200
c
      do 140 i = 2, n
         l = i - 1
         iz = (i * l) / 2
         ik = iz + i
         h = a(ik)
         if (h .eq. zero) go to 140
c
            do 130 j = 1, m
               s = -sdot(l,a(iz+1),1,z(1,j),1)
c
c              ********** double division avoids possible underflow ****
               s = (s / h) / h
c
               call saxpy(l,s,a(iz+1),1,z(1,j),1)
c
  130       continue
c
  140 continue
c
  200 return
      end

c 27 oct 1980 ste
c182********************************************
      subroutine eimqlv(n,d,e,e2,w,ind,ierr,rv1)
      implicit real*8(a-h,o-z)
c***********************************************
c
      double precision abs,anorm,b,c,d,dabs,dsign,dsqrt,e,e2,f,g,one
     +                ,p,parma,parmb,parm1,parm2,r,rv1,s,sign,sqrt,test
     +                ,two,w,zero
c
      dimension d(n),e(n),e2(n),w(n),rv1(n),ind(n)
c
      data zero/0.0d0/, one/1.0d0/, two/2.0d0/
csng  data zero/0.0e0/, one/1.0e0/, two/2.0e0/
c
      sqrt(parm1)=dsqrt(parm1)
      abs(parm2)=dabs(parm2)
      sign(parma,parmb)=dsign(parma,parmb)
c
c     ------------------------------------------------------------------
c
c     this subroutine is a modification of eispack subroutine imtqlv
c     which is a variant of  imtql1  which is a translation of
c     algol procedure imtql1, num. math. 12, 377-383(1968) by martin and
c     wilkinson, as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c
c     this subroutine finds the eigenvalues of a symmetric tridiagonal
c     matrix by the implicit ql method and associates with them
c     their corresponding submatrix indices.
c
c     the essence of this method is a process whereby a sequence of
c     symmetric tridiagonal matrices, unitarily similar to the
c     original symmetric tridiagonal matrix, is formed which
c     converges to a diagonal matrix.  the rate of convergence of
c     this sequence is improved by implicitly shifting the origin
c     at each iteration.  before each iteration, the symmetric tri-
c     diagonal matrix is checked for a possible splitting into
c     submatrices.  if a splitting occurs, only the uppermost sub-
c     matrix participates in the next iteration.  the eigenvalues
c     are ordered in ascending order as they are found.
c
c     the origin shift at each iteration is the eigenvalue of the
c     current uppermost 2x2 principal minor closer to the first
c     diagonal element of this minor.  whenever the uppermost 1x1
c     principal submatrix finally splits from the rest of the
c     matrix, its element is taken to be an eigenvalue of the
c     oritinal matrix and the algorithm proceeds with the remaining
c     submatrix.  this process is continued until the matrix has
c     split completely into submatrices of order 1.  the
c     tolerances in the splitting tests are proportional to the
c     relative machine precision, although the test is carried out
c     in a completely portable manner.
c
c
c     on input-
c
c        n is the order of the matrix,
c
c        d contains the diagonal elements of the input matrix,
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary,
c
c        e2 contains the squares of the corresponding elements of e.
c          e2(1) is arbitrary.
c
c     on output-
c
c        d and e are unaltered,
c
c        elements of e2, corresponding to elements of e regarded
c          as negligible, have been replaced by zero causing the
c          matrix to split into a direct sum of submatrices.
c          e2(1) is also set to zero,
c
c        w contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues,
c
c        ind contains the submatrix indices associated with the
c          corresponding eigenvalues in w -- 1 for eigenvalues
c          belonging to the first submatrix from the top,
c          2 for those belonging to the second submatrix, etc.,
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations,
c
c        rv1 is a temporary storage array.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
      ierr = 0
      k = 0
      itag = 0
c
      do 100 i = 1, n
         w(i) = d(i)
         if (i .ne. 1) rv1(i-1) = e(i)
  100 continue
c
      e2(1) = zero
      rv1(n) = zero
c
      do 290 l = 1, n
         j = 0
c
c        ********** look for small sub-diagonal element **********
c
  105    continue
         do 110 m = l, n
            if (m .eq. n) go to 120
            anorm = abs(w(m)) + abs(w(m+1))
            test  = anorm + abs(rv1(m))
            if(test .eq. anorm) go to 120
c
c
c        ********** guard against underflowed element of e2 **********
            if (e2(m+1) .eq. zero) go to 125
  110    continue
c
  120    continue
         if (m .le. k) go to 130
            if (m .ne. n) e2(m+1) = zero
  125       continue
            k = m
            itag = itag + 1
  130    continue
         p = w(l)
         if (m .eq. l) go to 215
            if (j .eq. 30) go to 1000
            j = j + 1
c
c           ********** form shift **********
c
            g = (w(l+1) - p) / (two * rv1(l))
            r = sqrt(g*g+one)
            g = w(m) - p + rv1(l) / (g + sign(r,g))
            s = one
            c = one
            p = zero
            mml = m - l
c
c           ********** for i=m-1 step -1 until l do -- **********
c
            do 200 ii = 1, mml
               i = m - ii
               f = s * rv1(i)
               b = c * rv1(i)
               if (abs(f) .lt. abs(g)) go to 150
                  c = g / f
                  r = sqrt(c*c+one)
                  rv1(i+1) = f * r
                  s = one / r
                  c = c * s
                  go to 160
c
  150          continue
               s = f / g
               r = sqrt(s*s+one)
               rv1(i+1) = g * r
               c = one / r
               s = s * c
  160          continue
               g = w(i+1) - p
               r = (w(i) - g) * s + two * c * b
               p = s * r
               w(i+1) = g + p
               g = c * r - b
  200       continue
c
            w(l) = w(l) - p
            rv1(l) = g
            rv1(m) = zero
            go to 105
c
c        ********** order eigenvalues **********
c
  215    continue
         if (l .eq. 1) go to 250
c
c           ********** for i=l step -1 until 2 do -- **********
c
            do 230 ii = 2, l
               i = l + 2 - ii
               if (p .ge. w(i-1)) go to 270
               w(i) = w(i-1)
               ind(i) = ind(i-1)
  230       continue
c
  250    continue
         i = 1
  270    continue
         w(i) = p
         ind(i) = itag
  290 continue
c
      go to 1001
c
c        ********** set error -- no convergence to an
c                   eigenvalue after 30 iterations **********
c
 1000 ierr = l
 1001 return
      end

c 27 oct 1980 ste
c182**********************************
      subroutine etred3(n,nv,a,d,e,e2)
      implicit real*8(a-h,o-z)
c*************************************
c
      double precision a,abs,d,dabs,dsign,dsqrt,dt,e,e2,f,g,h,hh,h1
     +                ,one,parma,parmb,parm1,parm2,sign,sqrt,scale
     +                ,scale1,zero
c
      dimension a(nv),d(n),e(n),e2(n)
c     level 2,a
c
      data zero/0.0d0/, one/1.0d0/
csng  data zero/0.0e0/, one/1.0e0/
c
      sqrt(parm1)=dsqrt(parm1)
      abs(parm2)=dabs(parm2)
      sign(parma,parmb)=dsign(parma,parmb)
c
c     ------------------------------------------------------------------
c
c     this subroutine is a modification of eispack subroutine tred3
c     which is a translation of the algol procedure tred3,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix, stored as
c     a one-dimensional array, to a symmetric tridiagonal matrix
c     using orthogonal similarity transformations.
c
c     the tridiagonal reduction is performed in the following way.
c     starting with j=n, the elements in the j-th row to the
c     left of the diagonal are first scaled, to avoid possible
c     underflow in the transformation that might result in severe
c     departure from orthogonality.  the sum of squares  sigma  of
c     these scaled elements is next formed.  then, a vector  u  and
c     a scalar
c                    h = u(transpose) * u / 2
c     define an operator
c                    p = i - u * u(transpose) / h
c     which is orthogonal and symmetric and for which the
c     similiarity transformation  pap  eliminates the elements in
c     the j-th row of  a  to the left of the subdiagonal and the
c     symmetrical elements in the j-th column.
c
c     the non-zero components of  u  are the elements of the j-th
c     row to the left of the diagonal with the last of them
c     augmented by the square root of  sigma  prefixed by the sign
c     of the subdiagonal element.  by storing the transformed sub-
c     diagonal element in  e(j)  and not overwriting the row
c     elements eliminated in the transformation, full information
c     about  p  is save for later use in  etrbk3.
c
c     the transformation sets  e2(j)  equal to  sigma  and  e(j)
c     equal to the square root of  sigma  prefixed by sign opposite
c     to that of the replaced subdiagonal element.
c
c     the above steps are repeated on further rows of the
c     transformed  a  in reverse order until  a  is reduced to tri-
c     diagonal form, that is, repeated for  j = n-1,n-2,...,3.
c
c
c     on input-
c
c        n  is the order of the matrix.
c
c        nv must be set to the dimension of the array parameter a
c           as declared in the calling program dimension statement.
c
c        a  contains the lower triangle of the real symmetric
c           input matrix, stored row-wise as a one-dimensional
c           array, in its first n*(n+1)/2 positions.
c
c     on output-
c
c        a  contains information about the orthogonal
c           transformations used in the reduction.
c
c        d  contains the diagonal elements of the tridiagonal matrix.
c
c        e  contains the subdiagonal elements of the tridiagonal
c           matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c           e2 may coincide with e if the squares are not needed.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
c     ********** for i=n step -1 until 1 do -- **********
c
      np1 = n + 1
      do  300 ii = 1, n
         i = np1 - ii
         l = i - 1
         iz = (i * l) / 2
         h = zero
         scale = zero
         if (l .lt. 1) go to 130
c
c           ********** scale row (algol tol then not needed) **********
c
            do 120 k = 1, l
               iz = iz + 1
               d(k) = a(iz)
               scale = scale + abs(d(k))
  120       continue
c
            if (scale .ne. zero) go to 140
  130    continue
               e(i) = zero
               e2(i) = zero
               go to 290
c
  140    continue
         scale1 = one / scale
         do 150 k = 1, l
            d(k) = d(k) * scale1
            h = h + d(k) * d(k)
  150    continue
c
         e2(i) = scale * scale * h
         f = d(l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         a(iz) = scale * d(l)
         if (l .eq. 1) go to 290
            f = zero
c
            jk = 1
            do 240 j = 1, l
               jm1 = j - 1
               dt = d(j)
               g = zero
c
c              ********** form element of a*u **********
c
               if(jm1 .eq. 0) go to 190
                  do 180 k = 1, jm1
                     e(k) = e(k) + dt * a(jk)
                     g = g + d(k) * a(jk)
                     jk = jk + 1
  180             continue
  190          continue
               e(j) = g + a(jk) * dt
               jk = jk + 1
c
c              ********** form element of p **********
c
  240       continue
            f = zero
            h1 = one / h
            do 250 j = 1, l
               e(j) = e(j) * h1
               f = f + e(j) * d(j)
  250       continue
c
            hh = f / (h + h)
            jk = 0
c
c              ********** form reduced a **********
c
            do 270 j = 1, l
               f = d(j)
               g = e(j) - hh * f
               e(j) = g
c
               do 260 k = 1, j
                  jk = jk + 1
                  a(jk) = a(jk) - f * e(k) - g * d(k)
  260          continue
  270       continue
c
  290    continue
         d(i) = a(iz+1)
         a(iz+1) = scale * sqrt(h)
  300 continue
c
      return
      end

c 27 oct 1980 ste
c174*****************************************************************
      subroutine eisgiv(msgfl,n,nvect,lena,nv,a,b,ind,root,vect,ierr)
      implicit real*8(a-h,o-z)
c********************************************************************
c
      double precision a,b,root,vect
c
      dimension a(lena),b(n,8),ind(n),root(n),vect(nv,nvect)
c     level 2,a,vect
c
c     ------------------------------------------------------------------
c
c     eispack-based substitute for qcpe subroutine givens.
c
c     finds all eigenvalues and some eigenvectors of a real symmetric
c     matrix.   author.. c. moler and d. spangler, n.r.c.c., 4/1/79.
c     further modifications by s. t. elbert, ames laboratory-usdoe.
c
c      the original reference to the givens technique is in oak ridge
c      report number ornl 1574 (physics), by wallace givens.
c      the method as presented in this program consists of four steps,
c      all modifications of the original method...
c      first, the input matrix is reduced to tridiagonal form by the
c      householder technique (orthogonal similarity transformations).
c      second, the roots are located using the implicit ql method.
c      third, the vectors of the tridiagonal form are evaluated by the
c      inverse iteration technique.  fourth, the tridiagonal vectors
c      are rotated to vectors of the original array.
c      vectors for degenerate (or near-degenerate) roots are forced
c      to be orthogonal.
c
c
c     input-
c
c     msgfl = file where error messages will be printed.
c             if msgfl is 0, error messages will be printed on file 6.
c             if msgfl is negative, error messages will not be printed.
c     n     = order of matrix.
c     nvect = number of vectors desired.  0 .le. nvect .le. n.
c     lena  = dimension of  a  in calling routine.  must not be less
c             than (n*n+n)/2.
c     nv    = row dimension of vect.   n .le. nv.
c     a     = input matrix, columns of the upper triangle packed into
c             linear array of dimension n*(n+1)/2.
c     b     = scratch array, 8*n elements (note this is more than
c             previous versions of givens.)
c     ind   = integer scratch array of length n.
c
c     output-
c
c     a       destoryed.
c     root  = all eigenvalues, root(1) .le. ... .le. root(n).
c             (for other orderings, see below.)
c     vect  = eigenvectors for root(1),..., root(nvect).
c     ierr  = 0 if no error detected,
c           = k if iteration for k-th eigenvalue failed,
c           = -k if iteration for k-th eigenvector failed.
c             (failures should be very rare.  contact moler.)
c
c     calls modified eispack subroutines etred3, eimqlv, einvit, and
c     etrbk3.
c     the original eispack subroutines tred3, imtqlv, tinvit, and trbak3
c     were modified by the introduction of two subroutines from the
c     blas library - sdot and saxpy. internal machine dependent
c     constants have been removed from imtqlv and tinvit.  calling
c     sequences have been modified to conserve variable type and
c     provide improved diagnostics.  the code has been passed
c     by the bell labs pfort verifier.  for further details
c      see eispack users guide, b. t. smith et al, springer-verlag
c     lecture notes in computer science, vol. 6, 2-nd edition, 1976.
c     ------------------------------------------------------------------
c
      lmsgfl=msgfl
      if(msgfl.eq.0) lmsgfl=6
      ierr = n + 1
      if( (n*n+n)/2 .gt. lena) go to 810
c
c        reduce real symmetric matrix a to tridiagonal form
c
      call etred3(n,lena,a,b(1,1),b(1,2),b(1,3))
c
c        find all eigenvalues of tridiagonal matrix via implicit ql
c
      call eimqlv(n,b(1,1),b(1,2),b(1,3),root,ind,ierr,b(1,4))
      if (ierr .ne. 0) go to 820
c
c     code to order roots in descending order (largest first)...
c     k = n/2
c     b(1,3) = 2.0
c     do 150 i = 1, k
c        j = n+1-i
c        t = root(i)
c        root(i) = root(j)
c        root(j) = t
c150  continue
c
      if (nvect .le. 0) return
      if (nv .lt. n) go to 830
c
c        find eigenvectors of tridiagonal matrix via inverse iteration
c
      call einvit(nv,n,b(1,1),b(1,2),b(1,3),nvect,root,ind,
     +            vect,ierr,b(1,4),b(1,5),b(1,6),b(1,7),b(1,8))
      if (ierr .ne. 0) go to 840
c
c        find eigenvectors of symmetric matrix via back transformation
c
      call etrbk3(nv,n,lena,a,nvect,vect)
      return
c
c        error message section
c
  810 if(lmsgfl.lt.0) return
      write(lmsgfl,901)
      go to 890
c
  820 if(lmsgfl.lt.0) return
      write(lmsgfl,902)
      go to 890
c
  830 if(lmsgfl.lt.0) return
      write(lmsgfl,903)
      go to 890
c
  840 if(lmsgfl.lt.0) return
      write(lmsgfl,904)
c
 890  continue
      write(lmsgfl,900) n,nvect,lena,nv,ierr
      return
c
  900 format(26h0*** eisgiv paramaters ***/
     +       14h ***      n = ,i8,4h ***/
     +       14h ***  nvect = ,i8,4h ***/
     +       14h ***   lena = ,i8,4h ***/
     +       14h ***     nv = ,i8,4h ***/
     +       14h ***   ierr = ,i8,4h ***)
  901 format(37h value of lena is less than (n*n+n)/2)
  902 format(39h eimqlv has detected an error condition)
  903 format(18h nv is less than n)
  904 format(39h einvit has returned an error condition)
      end

c 19 nov 1979
c174****************************************
      subroutine saxpy(n,sa,sx,incx,sy,incy)
c*******************************************
      implicit real*8(a-h,o-z)
c
      double precision sx,sy,sa,zero
c
      dimension sx(1),sy(1)
c     level 2,sx,sy
c
      data zero/0.0d0/
csng  data zero/0.0e0/
c
c        this is a basic linear algegra subroutine (blas)
c
c        constant times a vector plus a vector.
c        uses unrolled loop for increments equal to one.
c
      if(n.le.0) return
      if(sa.eq.zero) return
      if(incx.eq.1 .and. incy.eq.1) go to 120
c
c           code for unequal increments or equal increments
c              not equal to 1
         ix=1
         iy=1
         if(incx.lt.0) ix=(-n+1)*incx + 1
         if(incy.lt.0) iy=(-n+1)*incy + 1
         do 110 i=1,n
            sy(iy)=sy(iy) + sa*sx(ix)
            ix=ix+incx
            iy=iy+incy
  110    continue
         return
c
c              code for both increments equal to 1
c
  120 continue
      m=mod(n,4)
      if(m.eq.0) go to 140
c
c                    clean up loop
         do 130 i=1,m
            sy(i  ) = sy(i  ) + sa*sx(i  )
  130    continue
         if(n.lt.4) return
c
  140 continue
      mp1=m+1
      do 150 i=mp1,n,4
         sy(i  ) = sy(i  ) + sa*sx(i  )
         sy(i+1) = sy(i+1) + sa*sx(i+1)
         sy(i+2) = sy(i+2) + sa*sx(i+2)
         sy(i+3) = sy(i+3) + sa*sx(i+3)
  150 continue
      return
      end

c 22 may 1980
c174**********************************
csng  function sdot(n,sx,incx,sy,incy)
      double precision
     *function sdot(n,sx,incx,sy,incy)
c*************************************
      implicit real*8(a-h,o-z)
c
      double precision stemp,sx,sy,zero
c
      dimension sx(1),sy(1)
c     level 2,sx,sy
c
      data zero/0.0d0/
csng  data zero/0.0e0/
c
c        this is a basic lenear algegra subroutine (blas)
c
c        form the dot product of two vectors.
c        use unrolled loops for increments equal to one.
c
      stemp=zero
      sdot=zero
      if(n.le.0) return
      if(incx.eq.1 .and. incy.eq.1) go to 120
c
c           code for unequal increments or equal increments
c              not equal to 1
c
         ix=1
         iy=1
         if(incx.lt.0) ix=(-n+1)*incx + 1
         if(incy.lt.0) iy=(-n+1)*incy + 1
         do 110 i=1,n
            stemp=stemp+sx(ix)*sy(iy)
            ix=ix+incx
            iy=iy+incy
  110    continue
         sdot=stemp
         return
c
c                 code for both increments equal to 1
c
  120 continue
      m=mod(n,5)
      if(m.eq.0) go to 140
c
c           clean-up loop
c
         do 130 i=1,m
            stemp=stemp+sx(i)*sy(i)
  130    continue
         if(n.lt.5) go to 160
c
  140 continue
      mp1=m+1
      do 150 i=mp1,n,5
         stemp=stemp + sx(i)*sy(i) + sx(i+1)*sy(i+1) + sx(i+2)*sy(i+2)
     *                             + sx(i+3)*sy(i+3) + sx(i+4)*sy(i+4)
  150 continue
  160 sdot=stemp
      return
      end


c 27 oct 1980 ste
c192*******************************************************************
      subroutine einvit(nm,n,d,e,e2,m,w,ind,z,ierr,rv1,rv2,rv3,rv4,rv6)
      implicit real*8(a-h,o-z)
c**********************************************************************
c
      double precision abs,anorm,anorm2,d,dabs,dsqrt,e,e2,epmach
     +                ,eps,eps2,eps3,eps4,half,one,order,parm1,parm2
     +                ,pt02,rv1,rv2,rv3,rv4,rv6,sdot,sqrt,test,u,uk,v,w
     +                ,xu,x0,x1,z,zero
c
      dimension d(n),e(n),e2(n),w(m),z(nm,m),rv1(n),rv2(n),rv3(n),rv4(n)
     +         ,rv6(n),ind(m)
c
      data zero/0.0d0/, pt02/0.02d0/, half/0.5d0/, one/1.0d0/
csng  data zero/0.0e0/, pt02/0.02e0/, half/0.5e0/, one/1.0e0/
c
      sqrt(parm1)=dsqrt(parm1)
      abs(parm2)=dabs(parm2)
c
c     ------------------------------------------------------------------
c
c     this subroutine is a modification of eispack subroutine tinvit
c     which is a translation of the inverse iteration technique
c     in the algol procedure tristurm by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
c
c     this subroutine finds those eigenvectors of a tridiagonal
c     symmetric matrix corresponding to specified eigenvalues,
c     using inverse iteration.
c
c     the calculations proceed as follows.  first, the  e2  array
c     is inspected for the presence of a zero element defining a
c     submatrix.  the eigenvalues belonging to this submatrix are
c     identified by their common submatrix index in  ind.
c
c     the eigenvectors of the submatrix are then computed by
c     inverse iteration.  first, the  lu  decomposition of the sub-
c     matrix with an approximate eigenvalue subtracted from its
c     diagonal elements is achieved by gaussian elimination using
c     partial pivoting.  the multipliers defining the lower
c     triangular matrix  l  are stored in the temporary array  rv4
c     and the upper triangular matrix  u  is stored in the three
c     temporary arrays  rv1,  rv2, and  rv3.  saving these
c     quantities in rv1,  rv2,  rv3, and  rv4 avoids repeating
c     the  lu  decomposition if further iterations are required.
c     an approximate vector, stored in  rv6, is computed starting
c     from an initial vector, and the norm of the approximate
c     vector is compared with a norm of the submatrix to determine
c     whether the growth is sufficient to accept it as an
c     eigenvector.  if this vector is accepted, its euclidian norm
c     is made 1.  if the growth is not sufficient, this vector is
c     used as the initial vector in computing the next approximate
c     vector.  this iteration process is repeated at most  5  times.
c
c     eigenvectors computed in the above way corresponding to well-
c     separated eigenvalues of this submatrix will be orthogonal.
c     however, eigenvectors corresponding to close eigenvalues of
c     this submatrix may not be satisfactorily orthogonal.  hence,
c     to insure orthogonal eigenvectors, each approximate vector is
c     made orthogonal to those previously computed eigenvectors
c     whose eigenvalues are close to the current eigenvalue.  if
c     the orthogonalization process produces a zero vector, a
c     column of the identity matrix is used as an initial vector
c     for the next iteration.
c
c     identical eigenvalues are perturbed slightly in an attempt to
c     obtain independent eigenvectors.  these perturbations are not
c     recorded in the eigenvalue array  w.
c
c     the above steps are repeated on each submatrix until all the
c     eigenvectors are computed.
c
c
c     on input-
c
c        nm must be set to the row dimension of two-dimensional
c           array parameters as declared in the calling program
c           dimension statement,
c
c        n  is the order of the matrix,
c
c        d  contains the diagonal elements of the input matrix,
c
c        e  contains the subdiagonal elements of the input matrix
c           in its last n-1 positions.  e(1) is arbitrary,
c
c        e2 contains the squares of the corresponding elements of e,
c           with zeros corresponding to negligible elements of e.
c           e(i) is considered negligible if it is not larger than
c           the product of the relative machine precision and the sum
c           of the magnitudes of d(i) and d(i-1).  e2(1) must contain
c           0.0 if the eigenvalues are in ascending order, or 2.0
c           if the eigenvalues are in descending order.  if  bisect,
c           tridib, or  eimqlv  has been used to find the eigenvalues,
c           their output e2 array is exactly what is expected here,
c
c        m  is the number of specified eigenvalues,
c
c        w  contains the m eigenvalues in ascending or descending order,
c
c        ind contains in its first m positions the submatrix indices
c           associated with the corresponding eigenvalues in w --
c           1 for eigenvalues belonging to the first submatrix from
c           the top, 2 for those belonging to the second submatrix, etc.
c
c     on output-
c
c        all input arrays are unaltered,
c
c        z  contains the associated set of orthonormal eigenvectors.
c           any vector which fails to converge is set to zero,
c
c        ierr is set to
c          zero       for normal return,
c          -r         if the eigenvector corresponding to the r-th
c                     eigenvalue fails to converge in 5 iterations,
c
c        rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
c     ********** epmach is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c
      eps = half**20
   50 continue
         epmach = eps
         eps = eps*half
         test = one + eps
      if(test .ne. one) go to 50
c
      ierr = 0
      if (m .eq. 0) return
      itag = 0
      order = one - e2(1)
      iq = 0
c
c        ********** establish and process next submatrix **********
c
  100 continue
      ip = iq + 1
c
      do 120 iq = ip, n
         if (iq .eq. n) go to 140
         if (e2(iq+1) .eq. zero) go to 140
  120 continue
c
c        ********** find vectors by inverse iteration **********
c
  140 continue
      itag = itag + 1
      is = 0
c
      do 920 ir = 1, m
         if (ind(ir) .ne. itag) go to 920
            its = 1
            x1 = w(ir)
            if (is .ne. 0) go to 510
c
c           ********** check for isolated root **********
c
               xu = one
               if (ip .ne. iq) go to 490
                  rv6(ip) = one
                  go to 870
c
  490          continue
               anorm = abs(d(ip))
               jp = ip + 1
c
               do 500 i = jp, iq
                  anorm = anorm + abs(d(i)) + abs(e(i))
  500          continue
c
c              ********** eps2 is the criterion for grouping,
c                         eps3 replaces zero pivots and equal
c                         roots are modified by eps3,
c                         eps4 is taken very small to avoid overflow ***
c
               anorm2 = anorm / float(iq - ip + 1)
               eps2 = pt02 * anorm2
               eps3 = epmach * anorm
               uk = float(iq-ip+1)
               eps4 = uk * eps3
               uk = eps4 / sqrt(uk)
               is = ip
c
  505          continue
               igroup = 0
               go to 520
c
c           ********** look for close or coincident roots **********
c
  510       continue
            if (abs(x1-x0) .ge. eps2) go to 505
            igroup = igroup + 1
            if (order * (x1 - x0) .le. zero) x1 = x0 + order * eps3
c
c           ********** elimination with interchanges and
c                      initialization of vector **********
c
  520       continue
            v = zero
c
            do 580 i = ip,iq
               rv6(i) = uk
               if (i .eq. ip) go to 560
                  if (abs(e(i)) .lt. abs(u)) go to 540
c              ********** warning -- a divide check may occur here if
c                         e2 array has not been specified correctly ****
                     xu = u / e(i)
                     rv4(i) = xu
                     rv1(i-1) = e(i)
                     rv2(i-1) = d(i) - x1
                     rv3(i-1) = zero
                     if (i .ne. iq) rv3(i-1) = e(i+1)
                     u = v - xu * rv2(i-1)
                     v = -xu * rv3(i-1)
                     go to 580
c
  540             continue
                  xu = e(i) / u
                  rv4(i) = xu
                  rv1(i-1) = u
                  rv2(i-1) = v
                  rv3(i-1) = zero
  560          continue
               u = d(i) - x1 - xu * v
               if (i .ne. iq) v = e(i+1)
  580       continue
c
            if (u .eq. zero) u = eps3
            rv1(iq) = u
            rv2(iq) = zero
            rv3(iq) = zero
c
c           ********** back substitution
c                      for i=iq step -1 until ip do -- **********
c
  600       continue
            do 620 ii = ip,iq
               i = ip + iq - ii
               rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
               v = u
               u = rv6(i)
  620       continue
c
c           ********** orthogonalize with respect to previous
c                      members of group **********
c
            if (igroup .eq. 0) go to 700
               j = ir
c
               do 680 jj = 1,igroup
  630             j = j - 1
                  if (ind(j) .ne. itag) go to 630
                  xu = sdot(iq-ip+1,rv6(ip),1,z(ip,j),1)
c
                  call saxpy(iq-ip+1,-xu,z(ip,j),1,rv6(ip),1)
c
  680          continue
c
  700       continue
            anorm = zero
c
            do 720 i = ip,iq
               anorm = anorm + abs(rv6(i))
  720       continue
c
            if (anorm .ge. one) go to 840
c
c              ********** forward substitution **********
c
               if (its .eq. 5) go to 830
                  if (anorm .ne. zero) go to 740
                     rv6(is) = eps4
                     is = is + 1
                     if (is .gt. iq) is = ip
                     go to 780
c
  740             continue
                  xu = eps4 / anorm
c
                  do 760 i =ip,iq
                     rv6(i) = rv6(i) * xu
  760             continue
c
c                 ********** elimination operations on next vector
c                            iterate **********
c
  780             continue
                  do 820 i = jp,iq
                     u = rv6(i)
c
c                 ********** if rv1(i-1) .eq. e(i), a row interchange
c                            was performed earlier in the
c                            triangularization process **********
c
                     if (rv1(i-1) .ne. e(i)) go to 800
                        u = rv6(i-1)
                        rv6(i-1) = rv6(i)
  800                continue
                     rv6(i) = u - rv4(i) * rv6(i-1)
  820             continue
c
                  its = its + 1
                  go to 600
c
c                 ********** set error -- non-converged eigenvector ****
c
  830          continue
               ierr = -ir
               xu = zero
               go to 870
c
c              ********** normalize so that sum of squares is
c                         1 and expand to full order **********
c
  840       continue
            u = zero
c
            do 860 i = ip,iq
               u = u + rv6(i)**2
  860       continue
c
            xu = one / sqrt(u)
c
  870       continue
            do 880 i = 1, n
               z(i,ir) = zero
  880       continue
c
            do 900 i = ip,iq
               z(i,ir) = rv6(i) * xu
  900       continue
c
            x0 = x1
  920 continue
c
      if (iq .lt. n) go to 100
      return
      end
