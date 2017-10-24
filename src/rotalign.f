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



cc *********** this version uses xyz convention ************
cc xxx put x-convention code in here and toss old routines
c      subroutine rotalignxxxxxx(x,y,z,n,xb,yb,zb,itermax,iwrt,a,b,c)
c      implicit real*8(a-h,o-z)
c      dimension x(n),y(n),z(n)
c      dimension xb(n),yb(n),zb(n)
c
cc xxx put iatommatch into this call
c
cc      write(84,*) 'entered rotalign'
c      
c      call rotalign2(x,y,z,n,xb,yb,zb,itermax,iwrt,iconverged,a,b,c)
c
c      return
c      end

      subroutine rotalign2(x,n,xb,itermax,iwrt,iconverged,phi
     $     ,theta,psi,iz)
      implicit real*8(a-h,o-z)
c shift and rotate A cluster (x,y,z) to match B cluster (xb,yb,zb) as closely as possible  - afv 9/2001
c      dimension x(n),y(n),z(n)
      dimension x(3*n)
c      dimension xb(n),yb(n),zb(n)
      dimension xb(3*n)
c      parameter (nmax=10000)
      include 'parameters.h'
c      dimension xt(nmax),yt(nmax),zt(nmax)    ! xxx move to a work common  xxxx
      dimension xt(3*nmaxrotalign)
      dimension a(3,3)
      dimension av(3)
      dimension v(3)
      dimension vb(3)
      dimension dadphi(3,3)
      dimension dadtheta(3,3)
      dimension dadpsi(3,3)
      dimension xsum(3),xcmb(3),xcma(3),xshift(3),xbshift(3)

      if(n.gt.nmaxrotalign) stop 'increase nmax in rotalign2'

c      write(6,690) itermax,iwrt,phi,theta,psi
c 690  format('RA: DB: p: ',2i5,3f12.5)
      
      pi=2.0d0*acos(0.0d0)
      twopi=2.0d0*pi

      gfac=0.1d0
      gfac=1.0d0*60.0d0/float(n)
c      iwrt=-1

c      write(6,*) 'ROTALIGN_DEBUG: iwrt: ',iwrt
      
c find center of "mass" of B
      do k=1,3
         xsum(k)=0.0d0
      enddo
      do i=1,n
         do k=1,3
            xsum(k) = xsum(k) + xb((i-1)*3+k)
         enddo
      end do
      do k=1,3
         xcmb(k)=xsum(k)/float(n)
      enddo

c find center of mass of A
      do k=1,3
         xsum(k)=0.0d0
      enddo
      do i=1,n
         do k=1,3
            xsum(k) = xsum(k) + x((i-1)*3+k)
         enddo
      end do
      do k=1,3
         xcma(k)=xsum(k)/float(n)
      enddo

c shift A to have zero center of mass, or match an atom in b

      iatommatch=0              ! =0 means use center of mass

      if(iatommatch.eq.0) then
         do k=1,3
            xshift(k)=xcma(k)
            xbshift(k)=xcmb(k)
         enddo
      else
         do k=1,3
            xshift(k)=x((iatommatch-1)*3+k)
            xbshift(k)=xb((iatommatch-1)*3+k)
         enddo
      end if

      do i=1,n
         do k=1,3
            x((i-1)*3+k)=x((i-1)*3+k)-xshift(k)
         enddo
      end do


c now do euler rotations to bring into alignment by minimizing squared displacement

c      phi=0.0d0    ! in degrees
c      theta=0.0d0  ! in degrees
c      psi=0.0d0    ! in degrees

c      write(6,*) 'Rotalign: initial angles: ',phi,theta,psi
      
      do 100 iter=1,itermax

c compute derivatives of squared displacement w.r.t. euler angles

c rotation matrix info
         a1=phi*pi/180.0d0
         a2=theta*pi/180.0d0
         a3=psi*pi/180.0d0
         c1=cos(a1)
         s1=sin(a1)
         c2=cos(a2)
         s2=sin(a2)
         c3=cos(a3)
         s3=sin(a3)
         
         a(1,1)=c1*c2
         a(1,2)=s1*c2
         a(1,3)=-s2
         a(2,1)=s3*s2*c1-c3*s1
         a(2,2)=s3*s2*s1+c3*c1
         a(2,3)=c2*s3
         a(3,1)=c3*s2*c1+s3*s1
         a(3,2)=c3*s2*s1-s3*c1
         a(3,3)=c2*c3


c notes on how to do this
c xrot(i) =  a11*x(i)+a12*y(i)+a13*z(i) 
c yrot(i) =  a21*x(i)+a22*y(i)+a23*z(i) 
c zrot(i) =  a31*x(i)+a32*y(i)+a33*z(i) 

c minimize sum over i of (xb(i)-xrot(i))**2 + (yb(i)-yrot(i))**2 + (zb(i)-zrot(i))**2
c by expressing each aij in terms of euler angles and differentiating w.r.t them

c ds/dphi = dTOTxDIST/dphi = sum over atoms i of 2*(xrot(i)-xb(i))*dxrot(i)/dphi
c  where dxrot(i)/dphi = x element of [[dAdphi]]*[xyz(i)]     (matrix times column)


c dA/dphi terms  - each times pi/180.0d0
         dadphi(1,1) = -s1*c2
         dadphi(1,2) = c1*c2
         dadphi(1,3) = 0.0d0
         dadphi(2,1) = -s1*s2*s3 - c1*c3
         dadphi(2,2) = c1*s2*s3 - s1*c3
         dadphi(2,3) = 0.0d0 
         dadphi(3,1) = -s1*s2*c3 + c1*s3
         dadphi(3,2) = c1*s2*c3 + s1*s3
         dadphi(3,3) = 0.0d0


c dA/dtheta terms  - each times pi/180.0d0
         dadtheta(1,1)=-c1*s2
         dadtheta(1,2)=-s1*s2
         dadtheta(1,3)=-c2
         dadtheta(2,1)=s3*c2*c1
         dadtheta(2,2)=s3*c2*s1
         dadtheta(2,3)=-s2*s3
         dadtheta(3,1)=c3*c2*c1
         dadtheta(3,2)=c3*c2*s1
         dadtheta(3,3)=-s2*c3

c dA/dpsi terms  - each times pi/180.0d0
         dadpsi(1,1)=0.0d0
         dadpsi(1,2)=0.0d0
         dadpsi(1,3)=0.0d0
         dadpsi(2,1)=c3*s2*c1+s3*s1
         dadpsi(2,2)=c3*s2*s1-s3*c1
         dadpsi(2,3)=c2*c3
         dadpsi(3,1)=-s3*s2*c1+c3*s1
         dadpsi(3,2)=-s3*s2*s1-c3*c1
         dadpsi(3,3)=-c2*s3

c compute A'x to get derivs

         dsdphi=0.0d0
         dsdtheta=0.0d0
         dsdpsi=0.0d0
         do 30 ii=1,n
            do k=1,3
               v(k)=x((ii-1)*3+k)
               vb(k)=xb((ii-1)*3+k)-xbshift(k)
            enddo
            do 20 i=1,3
               avi=0.0d0
               aviphi=0.0d0
               avitheta=0.0d0
               avipsi=0.0d0
               do 10 j=1,3
                  avi=avi+a(i,j)*v(j)
                  aviphi=aviphi+dadphi(i,j)*v(j)
                  avitheta=avitheta+dadtheta(i,j)*v(j)
                  avipsi=avipsi+dadpsi(i,j)*v(j)
 10            continue
               vroti= avi
               dsdphi=dsdphi + 2.d0*(vroti-vb(i))*aviphi*twopi/360.d0
               dsdtheta=dsdtheta + 2.d0*(vroti-vb(i))*avitheta*twopi/360
     +             .d0
               dsdpsi=dsdpsi + 2.d0*(vroti-vb(i))*avipsi*twopi/360.0d0
 20         continue
 30      continue

         if (iz.eq.1) then ! rotate along z axis only
            dsdpsi=0.0d0
            dsdtheta=0.0d0
         endif
         
         do 130 ii=1,n
            do k=1,3
               v(k)=x((ii-1)*3+k)
            enddo
            do 120 i=1,3
               avi=0.0d0
               do 110 j=1,3
                  avi=avi+a(i,j)*v(j)
 110           continue
               av(i)=avi
 120        continue
            do k=1,3
               xt((ii-1)*3+k)=av(k)
            enddo
 130     continue
         
         error=0.0d0
         do i=1,n
            do k=1,3
               error=error+(xt((i-1)*3+k)-xb((i-1)*3+k)+xbshift(k))**2
            enddo
         end do
         rmserror=sqrt(error/float(n))

         if (iz.eq.0) then
            gradlen=sqrt(dsdphi**2 + dsdtheta**2 + dsdpsi**2)
         else
            gradlen=sqrt(dsdphi**2)
         endif

c         write(6,*) 'RA: DB: rms: ',iter,rmserror
         
c         if (iter.eq.1.or.iter.eq.itermax.or.iter.gt.0) then
c            do i=1,2
c               write(6,65) i,iter,xt(i),xb(i),xbshift,xshift,x(i),xt(i)
c     +             -xb(i)+xbshift
c               write(6,66) i,iter,yt(i),yb(i),ybshift,yshift,y(i),yt(i)
c     +             -yb(i)+ybshift
c               write(6,67) i,iter,zt(i),zb(i),zbshift,zshift,z(i),zt(i)
c     +             -zb(i)+zbshift
c 65            format('RA: DB: x: ',2i6,6f12.5)
c 66            format('RA: DB: y: ',2i6,6f12.5)
c 67            format('RA: DB: z: ',2i6,6f12.5)
c            enddo
c         endif

         gradlen=sqrt(dsdphi**2 + dsdtheta**2 + dsdpsi**2)

         if(iwrt.ge.1) then
c      if(iter.eq.1) write(6,*)'**note Euler angles are xyz-convention**'  
c      write(6,50) iter,phi,theta,psi, dsdphi,dsdtheta,dsdpsi,rmserror
 50         format(i6,1x,3f9.3, 2x,1p3d12.2,2x,1pd12.3)
         end if

         if(rmserror.lt.1.d-5 .or. gradlen.lt.1.d-5) then
            if(iwrt.ge.10) then
               write(6,*) 'rotalign converged: iter=',iter,' 
     x             rms mismatch =',rmserror,' gradlen=',gradlen
            end if
            iconverged=1
            go to 150
         end if

c update euler angles

         phi=phi - dsdphi*gfac
         if (iz.eq.0) then
            theta=theta - dsdtheta*gfac
            psi=psi - dsdpsi*gfac
         endif

 100  continue
      if(iwrt.ge.6) write(6,*) 'rotalign not converged yet'
      if(iwrt.ge.4 .and. itermax.gt.100) then
         write(6,*) '*** rotalign did not converge: iter: ',iter
     +       ,', rms: ',rmserror
      end if
      iconverged=0
 150  continue

      do i=1,n
         do k=1,3
            x((i-1)*3+k)=xt((i-1)*3+k)+xbshift(k)
         enddo
      end do

      nupost=nupost+1
      iupost=iconverged
      dupost=rmserror
      
c      write(6,691) phi,theta,psi
c 691  format('RA: DB: a: ',3f12.5)
      
c      call flushtad(6)

      return
      end


