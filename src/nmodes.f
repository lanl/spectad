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


      subroutine nmodes(ip,ivec,
     x                     natom,nmove,xyz,itype,taxes,ietype,
     x                     maxtyp,rcut,gdelt, amu,  e,grad,curv, 
     x                     nneg,nprod,freqlogs)

c except for the name, this should be identical to nmodesxyz in clsman  2/02

c xyz-style nmodes routine, adapted from nmodesbig, 2/15/02
c now curv will be passed in from the very top
c eigvectors, if needed (if ivec.eq.1) will be held in work common here
c feynman/wigner stuff removed to make it look simpler and cleaner 
c even amu is now passed in
c calling routine should take responsibility for curv being large enough

c comment cards from nmodesbig:
cc this version has a separate, large, new work common (wrkcm21) to hold hessian
cc so that it can call the new numerical grad-based hessian routine without
cc causing any conflict for wrkcm1
cc  if desired: just make this new work common big when we want to make a 
cc  special version of clsman, and then make it small again.    afv 8/3/01
c
cc  I made as few changes as possible from nmodesbig.  Basically, I just changed 
cc  wrkcm1 to wrkcm21 and commented out the claims to wrkcm1, and then changed the
cc  call to cgcalc to a call to cgcalcnew
c
cc  I then had to reinsert wrkcm1 to hold the eigenvectors, and discovered that
cc  wrkm1 is not large enough to hold them unless movmax is reduced to about 365.
cc  Sigh - still need to think this through some more.
c
cc  compute normal mode frequencies via diagonalization
cc  sets up and diagonalizes full curvature matrix/sqrt(mi*mj)
cc  this routine uses cgsubs routines, which use trpfun
cc  if ivec=0 - do not bother calculating eigenvectors (much faster)

      use mod_mpi
      implicit real*8(a-h,o-z)

      dimension xyz(3,natom),itype(natom)
      dimension grad(*),curv(*)
      dimension taxes(3)
      dimension rcut(maxtyp,maxtyp)
      dimension amu(maxtyp)

c      parameter (movmax=10000)
      include 'parameters.h'
      parameter (movmaxvec=365)
      parameter (m3=movmax*3)
      parameter (m3v=movmaxvec*3)

c      parameter(lwrk1=1200 000)
c      parameter(lwrk2=100 000)
      include 'wrksubs_tad2.h'
      common/wrkcm1/iwrk1,jwrk1,eigvec(m3v*m3v),
     x              fill1(lwrk1 - m3v*m3v)
      common/wrkcm2/iwrk2,jwrk2, work2(m3),work3(m3),work4(m3),
     x              fill(lwrk2-3*m3)

      common/modcom/flogs,flogsl
      character*20 filnam

c size checks
      if(nmove.gt.movmax) then
        if(rank_f.eq.0)then
         write(6,*) 'nmodesxyz: nmove too large. max=',movmax
        endif
        return
      end if
      if(ivec.gt.0 .and. nmove.gt.movmaxvec) then
        if(rank_f.eq.0)then
         write(6,*) 'nmodesxyz: nmove w/vec too large. max=',movmax
        endif
        return
      end if
      if(nmove.le.0) then
       if(rank_f.eq.0)then
        write(6,*) 'nmodesxyz wont work till optset has been called'
       endif
       return
      end if

      call wrkclaimv(2,1,'nmodesxyz')

c  compute gradient and curvature
      call cgcalc(natom,nmove,xyz,itype,taxes,ietype,
     x                     maxtyp,rcut,gdelt,e,grad,curv)

c  divide curvature elements by sqrt(mi*mj)
c  units:  hartree/(bohr*atomic units of mass) (i.e. - all atomic units)
      k=0
      trace=0.0d0
      twopi=2.0d0*3.1415926d0
      do 40 iat=1,nmove
      rmi=1822.83d0*amu(itype(iat))
      do 30 ix=1,3
      i=3*(iat-1)+ix
      do 20 jat=1,nmove
      rmj=1822.83d0*amu(itype(jat))
      rmij=sqrt(rmi*rmj)
      do 10 jx=1,3
      j=3*(jat-1)+jx
      if(j.gt.i) go to 10
        k=k+1
        curv(k)=(curv(k)/rmij)*0.5291d0**2
        if(i.eq.j) trace=trace+curv(k)/(2.418d-17*twopi)**2
   10 continue
   20 continue
   30 continue
   40 continue

c  diagonalize mass-weighted curvature matrix
      n=nmove*3
      nm=n
      lencur=m3v*(m3v+1)/2
      if(ivec.gt.0) then
        call wrkclaimv(1,1,'nmodesxyz')
        ierr=0
        call rsp(nm,n,lencur,curv,work2,1,eigvec,work3,work4,ierr)
        if(rank_f.eq.0)then
         if(ierr.ne.0)write(6,*)'nmodesxyz - warning: rsp failed***'
        endif
      else if(ivec.eq.0) then
        call sdiagf(n,curv,work2,work3)
      end if

c  convert diagonal elements to frequencies in sec**-1
      twopi=2.0d0*3.1415926d0
      flogs=0.0d0
      nprod=0
      nneg=0
      do 60 i=1,nmove*3
      dag=work2(i)
      if(dag.le.0.0d0) nneg=nneg+1
      dag=sqrt(abs(dag))*sign(1.0d0,dag)/twopi
      work2(i)=dag/2.418d-17
      if(work2(i).gt.0.0d0) then     !!!! xxx note .ge. is now .gt. 
        flogs=flogs + log(work2(i))
        nprod=nprod+1
      end if
   60 continue

      freqsquaresum=0.0d0
      do 65 i=1,nmove*3
      freqsquaresum = freqsquaresum + work2(i)*abs(work2(i)) 
   65 continue


c compute gradient projection along lowest mode
      if(ivec.gt.0) then
        proj1=0.0d0
        do 70 i=1,nmove*3
        proj1=proj1+grad(i)*eigvec(i)
   70   continue
        proj1=proj1*27.21d0
      end if

c  print out frequencies
      if(ip.ge.1) then
       if(rank_f.eq.0) write(6,80) (work2(i),i=1,nmove*3)
   80   format(/'  Normal mode frequencies (Hertz):'/
     x   (1p5d13.3))
      end if

c  print out eigenvectors (normal modes)
c      do 100 i=1,nmove*3
      do 100 i=1,99
      j=nmove*3*(i-1)
      if(rank_f.eq.0)then
       if(ip.ge.2) write(7,90) i,work2(i),(eigvec(j+k),k=1,nmove*3)
       if(ip.ge.0) write(filnam,92) i
       if(ip.ge.0)
     +    open (11,file=filnam,form='formatted',status='unknown')
       if(ip.ge.0) write(11,91) nmove,(eigvec(j+k),k=1,nmove*3)  ! xxx blas - I made this contingent on ip
       if(ip.ge.0) close(11)
      endif
   90 format(' normal mode',i4,' frequency =',1pd13.3,' Hz'
     x            / (1x,0p3f10.4,2x,3f10.4,2x,3f10.4))
 91    format(i10/(0p3f20.10))
 92    format('eig.',i2)
  100 continue
      if(rank_f.eq.0)then
       if(ip.ge.0)write(7,105) nneg,flogs,exp(flogs/float(nprod))
      endif
  105 format(' neg. curvs =',i3,'  ln(pos freqs) =',1pd20.10/
     x       '        geometric ave of pos freqs =',1pd15.5,' Hz')
      if(rank_f.eq.0)then
       if(ip.ge.0.and.nneg.eq.0 .and. ivec.gt.0)write(6,107) 'pos',proj1
       if(ip.ge.0.and.nneg.ge.1 .and. ivec.gt.0)write(6,107) 'neg',proj1
      endif
  107 format(' grad projection onto lowest mode ('
     x                    ,a3,' curv) =',1pd12.3,' eV/A')
      if(rank_f.eq.0)then
       if(flogsl.ne.0.0d0 .and. abs(nprod-nprodl).eq.1 .and. ip.ge.0)
     x     write(7,120) exp(flogs-flogsl),nprod-nprodl
      endif
  120 format(' Vineyard ratio of freqs product with previous product =',
     x     1pd13.3,' Hz**',i2)
      if(rank_f.eq.0)then
       if(flogsl.ne.0.0d0 .and. abs(nprod-nprodl).eq.0 .and. ip.ge.0)
     x     write(7,121) exp(-(flogs-flogsl))
      endif
  121 format(' ratio of partition function with previous case =',
     x     1pd13.3)

      nprodl=nprod
      flogsl=flogs
      freqlogs=flogs
      freqsquaresuml = freqsquaresum
      tracel= trace

      call wrkrelease(2)
      if(ivec.ge.1) then
        call wrkrelease(1)
        iwrk1=-22126
      end if

      if(rank_f.eq.0) call flushtad(7)
      
      return
      end


      subroutine cgcalc(natom,nmove,xyz,itype,taxes,ietype,
     x                     maxtyp,rcut,gdelt,  e,grad,curv)
c cgcalcxyz - afv 2/15/02 - a cgcalc routine compatible with tad2 
c (e.g., xyz arrays, no or few commons)  Based on cgcalcn in testsubs.f
c comments from cgcalcnew:
c  8/3/01 - w,lenw removed, since they were not being used
c  this routine returns the curvature and gradient over the first nmove atoms
c  centered numerical derivatives using forces   6/2001
c  note that x,y,z are assumed to be 19000 atoms apart!
c  note that this routine does all one-atom moves (see ea routine)

      implicit real*8(a-h,o-z)
      dimension xyz(3,natom),itype(natom)
      dimension grad(*),curv(*)
      dimension taxes(3)
      dimension rcut(maxtyp,maxtyp)
c      parameter (maxmove=10000) ! consistent with elsewhere
      include 'parameters.h'

      include 'wrksubs_tad2.h'
c      parameter(lwrk8=100 000, maxmove=15000)
      common/wrkcm8/iwrk8,jwrk8,gp(3,maxmove),gm(3,maxmove),
     x                     wrk8(lwrk8-6*maxmove)

      call wrkclaimv(8,1,'cgcalcn')
      if(nmove.gt.maxmove) stop 'cgcaln - maxmove exceeded'

c compute gradient
      call gcalc(natom,nmove,xyz,itype, taxes,ietype,maxtyp,rcut,
     x                                          e,grad, hyperratio)

c zero hessian
      n=nmove*3
      do 20 ij=1,n*(n+1)/2
      curv(ij)=0.0d0
   20 continue

c compute hessian numerically using repeated gradient calls
      delt=gdelt

      do 60 i=1,nmove
      do 60 ix=1,3
      idx=3*(i-1)+ix

      xi0=xyz(ix,i)
      xyz(ix,i)=xi0+delt
      call gcalc(natom,nmove,xyz,itype, taxes,ietype,maxtyp,rcut,
     x                                      etemp,gp, hyperratio)
      xyz(ix,i)=xi0-delt
      call gcalc(natom,nmove,xyz,itype, taxes,ietype,maxtyp,rcut,
     x                                      etemp,gm, hyperratio)
      xyz(ix,i)=xi0


      do 50 j=1,nmove
      do 50 jx=1,3
      jdx=3*(j-1)+jx
      ijx=idx*(idx-1)/2+jdx
      if(jdx.gt.idx) ijx=jdx*(jdx-1)/2+idx
      
      hess = (gp(jx,j)-gm(jx,j))/(2.d0*delt)
      if(jdx.eq.idx) then
        curv(ijx)=hess
      else
        curv(ijx)=curv(ijx) + 0.5d0*hess
      end if

   50 continue
   60 continue

c for debugging
c      nn=n*(n+1)/2
c      if(nmove.le.10) write(6,70) (curv(i),i=1,nn)
c   70 format('full hessian: ', / (1p5d13.3))
c      call wrkrelease(8)

c may have deleted some code from this routine by mistake
c as I made the nmodesxyz routine (?) xxx
      call wrkrelease(8)

      return
      end


