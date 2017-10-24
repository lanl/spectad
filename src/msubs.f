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


c: msubs - assorted utility routines of a matrix or array nature
c:

      subroutine dmnmx(x,natom,xmin,xmax,ixmin,ixmax)
c: finds min and max in array x, returns values(xmin,xmax) and pointers(ixmin,ixmax)
c:
      implicit real*8(a-h,o-z)
      dimension x(natom)
      ixmin=1
      ixmax=1
      do 10 i=1,natom
      if(x(i).lt.x(ixmin)) ixmin=i
      if(x(i).gt.x(ixmax)) ixmax=i
   10 continue
      xmin=x(ixmin)
      xmax=x(ixmax)
      return
      end


      subroutine dmnmxstride(x,natom,xmin,xmax,ixmin,ixmax,istride)
c: finds min and max in array x, returns values(xmin,xmax) and pointers(ixmin,ixmax)
c: runs through array with stride istride, looking at natom elements
c:
      implicit real*8(a-h,o-z)
      dimension x(natom*istride)
      ixmin=1
      ixmax=1
      do 10 i=1,natom*istride,istride
      if(x(i).lt.x(ixmin)) ixmin=i
      if(x(i).gt.x(ixmax)) ixmax=i
   10 continue
      xmin=x(ixmin)
      xmax=x(ixmax)
      return
      end

      subroutine rmnmx(x,natom,xmin,xmax,ixmin,ixmax)
c: finds min and max in array x, returns values(xmin,xmax) and pointers(ixmin,ixmax)
c:
      implicit real*4(a-h,o-z)
      dimension x(natom)
      ixmin=1
      ixmax=1
      do 10 i=1,natom
      if(x(i).lt.x(ixmin)) ixmin=i
      if(x(i).gt.x(ixmax)) ixmax=i
   10 continue
      xmin=x(ixmin)
      xmax=x(ixmax)
      return
      end

      function imax(ia,n)
c: finds max integer in array ia
c:
      dimension ia(n)
      if(n.eq.0) then
        imax=0
        return
      else
        imax=ia(1)
      end if
      do 10 i=1,n
      if(ia(i).gt.imax) imax=ia(i)
   10 continue
      return
      end

      function dmax(a,n)
c:  finds max value in real*8 array a
c:
      implicit real*8(a-h,o-z)
      dimension a(n)
      dmax=a(1)
      do 10 i=1,n
      if(a(i).gt.dmax) dmax=a(i)
   10 continue
      return
      end

      function dmin(a,n)
c:  finds min value in real*8 array a
c:
      implicit real*8(a-h,o-z)
      dimension a(n)
      dmin=a(1)
      do 10 i=1,n
      if(a(i).lt.dmin) dmin=a(i)
   10 continue
      return
      end


      subroutine dmfac(a,n,fac)
c:  multiplies each of the n elements in real*8 array a by fac
c:
      implicit real*8(a-h,o-z)
      dimension a(n)
      do 10 i=1,n
      a(i)=a(i)*fac
   10 continue
      return
      end


      subroutine dmave(a,n,ave)
c:  finds average of the n elements in real*8 array a
c:
      implicit real*8(a-h,o-z)
      dimension a(n)
      sum=0.0d0
      do 10 i=1,n
      sum=sum+a(i)
   10 continue
      ave=sum/float(n)
      return
      end

      subroutine dmzer(a,n)
c:  zeros each of the n elements in real*8 array a
c:
      implicit real*8(a-h,o-z)
      dimension a(n)
      do 10 i=1,n
      a(i)=0.0d0
   10 continue
      return
      end

      subroutine imzer(ia,n)
c:  zeros each of the n elements in integer array ia
c:
      implicit real*8(a-h,o-z)
      dimension ia(n)
      do 10 i=1,n
      ia(i)=0
   10 continue
      return
      end

      subroutine nanchk(x,n)
c:  checks whether x is a bad number (either inf or NaN)
c:
      implicit real*8(a-h,o-z)
      dimension x(n)
      do 10 i=1,n
      if(x(i).ne.x(i)) then
        write(6,*) 'bad number found:',x(i)
         stop 'nanchk - bad number found'
      end if
   10 continue
      return
      end

      subroutine lag4pt(x,xn,fn,f)
c:  4 point lagrangian interpolation (Abrom. and Stegun, pg 878)
c:
      implicit real*8(a-h,o-z)
      real*8 l(4)
      dimension xn(4),fn(4)

      l(1)=(x-xn(2))*(x-xn(3))*(x-xn(4)) /
     >     ( (xn(1)-xn(2))*(xn(1)-xn(3))*(xn(1)-xn(4)) )
      l(2)=(x-xn(1))*(x-xn(3))*(x-xn(4)) /
     >     ( (xn(2)-xn(1))*(xn(2)-xn(3))*(xn(2)-xn(4)) )
      l(3)=(x-xn(1))*(x-xn(2))*(x-xn(4)) /
     >     ( (xn(3)-xn(1))*(xn(3)-xn(2))*(xn(3)-xn(4)) )
      l(4)=(x-xn(1))*(x-xn(2))*(x-xn(3)) /
     >     ( (xn(4)-xn(1))*(xn(4)-xn(2))*(xn(4)-xn(3)) )

      f=0.0d0
      do 10 i=1,4
        f=f+l(i)*fn(i)
   10 continue

      return
      end


      subroutine mwrit(iunit,a,n)
c:  writes square matrix out to unit iunit
c:
      implicit real*8(a-h,o-z)
      dimension a(n,n)

      do 10 i=1,n
      write(iunit,20) (a(i,j),j=1,n)
   20 format(1x,1p10d10.2)
   10 continue

      return
      end

      subroutine mwritg(iunit,title,a,nrdim,nr,nc)
c:  writes rectangular matrix out to unit iunit
c:
      implicit real*8(a-h,o-z)
      character*(*) title
      dimension a(nrdim,nc)

c   allow up to 3 blank lines before matrix print starts
c   signaled by '/' in title card
      iskip=0
      do 2 i=1,3
      if(title(i:i).eq.'/') then
        write(iunit,*) ' '
        iskip=iskip+1
      end if
    2 continue

      write(iunit,5) title(1+iskip:len(title))
    5 format(1x,a)
      do 10 i=1,nr
      write(iunit,20) (a(i,j),j=1,nc)
   20 format(1x,1p10d10.2)
   10 continue

      return
      end

      subroutine mwritt(iunit,a,n)
c:  writes trianular matrix out to unit iunit
c:
      implicit real*8(a-h,o-z)
      dimension a(n)

      ii=0
      do 10 i=1,n
      write(iunit,20) (a(ii+j),j=1,i)
      ii=ii+i
   20 format(1x,1p10d10.2)
   10 continue
      return
      end

      subroutine moscan(a,n)
c:  scans square matrix for largest off-diagonal element
c:
      implicit real*8(a-h,o-z)
      dimension a(n,n)

      amax=0.0d0
      do 20 i=1,n
      do 10 j=1,n
      if(abs(a(i,j)).gt.amax .and. i.ne.j) then
        amax=abs(a(i,j))
        iimax=i
        jjmax=j
      end if
   10 continue
   20 continue

      write(6,50) iimax,jjmax,a(iimax,jjmax)
   50 format(' largest off-diagonal element:',i4,i4,1pd12.3)

      return
      end

      subroutine mvmult(a,v,c,n)
c:  performs matrix-vector multiply:  c=A*v
c:
      implicit real*8(a-h,o-z)
      dimension a(n,n),v(n),c(n)

      do 10 i=1,n
      do 10 j=1,n
      c(i)=c(i) + a(i,j)*v(j)
   10 continue

      return
      end

      subroutine mvmulth(a,v,c,n)
c:  performs matrix-vector multiply:  c=A*v, where A is upper symmetric half matrix
c:
      implicit real*8(a-h,o-z)
      dimension a(n*(n+1)/2),v(n),c(n)

c      do 30 i=1,n
c      v(i)=0.0d0
c      ij=i*(i-1)/2
c      do 10 j=1,i
c      ij=ij+1
c      c(i)=c(i) + a(ij)*v(j)
c   10 continue
c      do 20 j=i+1,n
c      ij=ij+j-1
c      c(i)=c(i) + a(ij)*v(j)
c   20 continue
c   30 continue

      do 10 i=1,n
      c(i)=0.0d0
   10 continue

      ij=0
      do 30 i=1,n
      do 20 j=1,i
      ij=ij+1
      c(i)=c(i) + a(ij)*v(j)
      c(j)=c(j) + a(ij)*v(i)
   20 continue
      c(i)=c(i) - a(ij)*v(i)
   30 continue

      return
      end



      subroutine vecaxpy(n,a,x,y)
c: performs y(i) = a*x(i) + y(i)  , i=1,n
c: like daxpy, but without increment option
c:
      implicit real*8(a-h,o-z)
      dimension x(n),y(n)

      do 10 i=1,n
      y(i)=y(i) + a*x(i)
   10 continue

      return
      end

      subroutine vecypax(y,a,x,n)
c: performs y(i) = y(i) + a*x(i) , i=1,n
c: like daxpy, but reordered input, and without increment option
c:
      implicit real*8(a-h,o-z)
      dimension x(n),y(n)

      do 10 i=1,n
      y(i)=y(i) + a*x(i)
   10 continue

      return
      end

      subroutine vecswap(a,b,n,work)
c:  swaps n=long vector a with n-long vector b using work space in work(n)
c:
      implicit real*8(a-h,o-z)
      dimension a(n),b(n),work(n)

      do 10 i=1,n
      work(i)=a(i)
      a(i)=b(i)
      b(i)=work(i)
   10 continue

      return
      end

      subroutine scalswap(a,b)
c:  swaps scalar a with scalar b 
c:
      implicit real*8(a-h,o-z)

      temp=a
      a=b
      b=a

      return
      end

      subroutine mchebrecur(t3,x,t2,t1,n)
c:  performs Chebyshev recursion matrix operation:  T3 = 2*X*T2 - T1
      implicit real*8(a-h,o-z)
      dimension t1(n,n),t2(n,n),t3(n,n)
      dimension x(n,n) !! rjz -- wrong?
      call mmmult(x,t2,t3,n)
      do 10 i=1,n
      do 10 j=1,n
      t3(i,j)=2.0d0*t3(i,j)-t1(i,j)
   10 continue

      return
      end

      subroutine mmmult(a,b,c,n)
c:  performs matrix multiply:  C=A*B  (simply calls mmab)
      implicit real*8(a-h,o-z)
      dimension a(n,n),b(n,n),c(n,n)

      call mmab(a,b,c,n)

      return
      end


      subroutine mmab(a,b,c,n)
c:  performs matrix multiply:  C=A*B
c:  optimized for HP 735 -- compile with +OP2, 3 or 4 (4 is probably best)
c:
      implicit real*8(a-h,o-z)
      dimension a(n,n),b(n,n),c(n,n)

      if(n.le.300) then

        do 10 i=1,n
        do 10 j=1,n
        c(i,j)=0.0d0
        do 10 k=1,n
        c(i,j)=c(i,j) + a(i,k)*b(k,j)
   10   continue

      else

c  transpose A
        do 20 i=1,n
        do 20 j=1,n
        temp=a(i,j)
        a(i,j)=a(j,i)
        a(j,i)=temp
   20   continue
c  multiply ATT*B
        do 30 i=1,n
        do 30 j=1,n
        c(i,j)=0.0d0
        do 30 k=1,n
        c(i,j)=c(i,j) + a(k,i)*b(k,j)
   30   continue
c  untranspose A
        do 40 i=1,n
        do 40 j=1,n
        temp=a(i,j)
        a(i,j)=a(j,i)
        a(j,i)=temp
   40   continue
      end if

      return
      end

      subroutine vvadd(a,b,c,n)
c:  performs vector-vector add:  c=a+b
c:
      implicit real*8(a-h,o-z)
      dimension a(n),b(n),c(n)

        do 10 i=1,n
        c(i) = a(i) + b(i)
   10   continue

      return
      end

      subroutine xxprint(a,b,c,n)
      implicit real*8(a-h,o-z)
      dimension a(n,n),b(n,n),c(n,n)
      sum=0.0d0
      do 10 i=1,n
      do 10 j=1,n
      sum=sum+c(i,j)
   10 continue
      write(6,*) 'sum= ', sum
      return
      end

      real*8 function vecdot(a,b,n)
c: vector dot product of a and b, length n
c:
      implicit real*8(a-h,o-z)
      dimension a(n),b(n)
      dot=0.0d0
      do 10 i=1,n
      dot=dot+a(i)*b(i)
   10 continue
      vecdot=dot
      return
      end

      subroutine vecrenorm(a,n)
c: renormalize vector a of length n
c:
      implicit real*8(a-h,o-z)
      dimension a(n)

      dot=0.0d0
      do 10 i=1,n
      dot=dot+a(i)*a(i)
   10 continue

      fix=1.0d0/sqrt(dot)
      do 20 i=1,n
      a(i)=fix*a(i)
   20 continue

      return
      end

      subroutine vecmsk(a,bmask,n)
c: replace each element of vector a with a(i)=a(i)*bmask(i)
c:
      implicit real*8(a-h,o-z)
      dimension a(n),bmask(n)
      dot=0.0d0
      do 10 i=1,n
      a(i)=a(i)*bmask(i)
   10 continue
      return
      end

      subroutine vecred(v,n,crit)
c: zero each element of v whose abs value is less than crit
c: red== "reduce"
c:
      implicit real*8(a-h,o-z)
      dimension v(n)
      do 10 i=1,n
      if(abs(v(i)).lt.crit) v(i)=0.0d0
   10 continue
      return
      end

      subroutine vecmov(a,b,n)
c: copy vector a into b, length n
c:
      implicit real*8(a-h,o-z)
      dimension a(n),b(n)
      do 10 i=1,n
      b(i)=a(i)
   10 continue
      return
      end

      subroutine ivecmov(ia,ib,n)
c: copy integer vector a into b, length n
c:
      implicit real*8(a-h,o-z)
      dimension ia(n),ib(n)
      do 10 i=1,n
      ib(i)=ia(i)
   10 continue
      return
      end
