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

      subroutine vineyardsubset(natom,nmove,xyz,xyz2,itype,taxes,ietype
     +    ,maxtyp,rcut,amu,rdynmatcrit,e,grad,prefac,nelements
     +    ,ierr,ipr)
      
      implicit real*8(a-h,o-z)
      
      dimension xyz(3,natom),xyz2(3,natom),itype(natom)
      dimension grad(*)
      dimension taxes(3)
      dimension rcut(maxtyp,maxtyp)
      dimension amu(maxtyp)
      
c      parameter (movmax=10000)
      include 'parameters.h'
      dimension iinclude(movmax)
      dimension rrmin(3),iimin(3)

      include 'wrksubs_tad2.h'
      common/wrkcm13/iwrk13,jwrk13,work13(lwrk13)
c      write(6,*) 'vineyardsubset: entered'
      
c find atoms to include in hessian
      nelements=0
      do jx=1,3
         rrmin(jx)=1.0d99
      enddo
      do ix=1,nmove
         iinclude(ix)=0
         rr=0.0d0
         do jx=1,3
            rr=rr+(xyz(jx,ix)-xyz2(jx,ix))**2
         enddo
         ixtest=ix
         rrtest=sqrt(rr)
         do jx=1,3
            if (rrtest.lt.rrmin(jx)) then
               rrhold=rrtest
               ixhold=ixtest
               rrtest=rrmin(jx)
               ixtest=iimin(jx)
               iimin(jx)=ixhold
               rrmin(jx)=rrhold
            endif
         enddo
         if (sqrt(rr).ge.rdynmatcrit) then
            iinclude(ix)=1
            nelements=nelements+1
c            write(6,*) 'vineyardsubset: atom ',ix,' meets crit'
         endif
      enddo
c freeze three atoms to remove translational and rotational modes
      nnotincluded=nmove-nelements
c      write(6,*) 'vineyardsubset: nmove,nelements,nnotincluded: ',nmove
c     +    ,nelements,nnotincluded,natom
      if (nnotincluded.lt.3.and.nmove.eq.natom) then
c         write(6,*) 'vineyardsubset: entered loop'
         do jx=nnotincluded+1,3
            iinclude(iimin(jx))=0
            nelements=nelements-1
            if (ipr.ge.6) write(6,*) 'vineyardsubset: removing atom '
     +          ,iimin(jx),' from hessian'
         enddo
      endif

      ip=-1
      ivec=0
      gdelt=1.d-4

      if (ipr.ge.4) write(6,*) 'vineyardsubset: rdynmatcrit: '
     +    ,rdynmatcrit,', nelements: ',nelements

      nn=nelements*3*(nelements*3+1)/2
      if(nn.gt.lwrk13) then
         write(6,*)
     +       'vineyardsubset: hessian array too small to call nmodes'
         write(6,*)
     +       'vineyardsubset: reducing size of hessian'
         call removeelements(natom,nmove,xyz,xyz2,nelements,lwrk13
     +       ,iinclude)
      end if

      if (nelements.gt.0) then

         call wrkclaimv(13,1,'vineyardsubset') ! for hessian curv array (work13)
         
         call nmodessubset(ip,ivec,natom,nmove,xyz,itype,taxes,ietype
     +       ,maxtyp,rcut,gdelt,amu,iinclude,nelements,e,grad,work13
     +       ,nnegmin,nprodmin,freqlogsmin)

         call nmodessubset(ip,ivec,natom,nmove,xyz2,itype,taxes,ietype
     +       ,maxtyp,rcut,gdelt,amu,iinclude,nelements,e,grad,work13
     +       ,nnegsad,nprodsad,freqlogssad)
         
         call wrkrelease(13)
         
         if (nnegmin.eq.0.and.nnegsad.eq.1) then
            prefac=exp(freqlogsmin-freqlogssad)
            ierr=0
         else
            write(6,*)
     +          'vineyard: something wrong, setting standard prefactor'
            write(6,*) 'vineyard: nnegmin, nnegsad: ',nnegmin,nnegsad
            prefac=1.0d12
            ierr=1
         endif
         
      else
         write(6,*) 'vineyardsubset: no elements! setting standard'
         prefac=1.0d12
      endif


      return
      end


      subroutine nmodessubset(ip,ivec,natom,nmove,xyz,itype,taxes,
     +    ietype,maxtyp,rcut,gdelt, amu,iinclude,nelements,e,grad,curv,
     +    nneg,nprod,freqlogs)

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

      include 'wrksubs_tad2.h'
c      parameter(lwrk1=1200 000)
c      parameter(lwrk2=100 000)
      common/wrkcm1/iwrk1,jwrk1,eigvec(m3v*m3v),
     x              fill1(lwrk1 - m3v*m3v)
      common/wrkcm2/iwrk2,jwrk2, work2(m3),work3(m3),work4(m3),
     x              fill(lwrk2-3*m3)

      common/modcom/flogs,flogsl
      character*20 filnam
      dimension iinclude(nmove)

c size checks
      if(nelements.gt.movmax) then
         write(6,*) 'nmodesxyz: nmove too large. max=',movmax
         return
      end if
      if(ivec.gt.0 .and. nmove.gt.movmaxvec) then
         write(6,*) 'nmodesxyz: nmove w/vec too large. max=',movmax
         return
      end if
      if(nmove.le.0) then
        write(6,*) 'nmodesxyz wont work till optset has been called'
        return
      end if

      call wrkclaimv(2,1,'nmodessubset')

c  compute gradient and curvature
      call cgcalcsubset(natom,nmove,xyz,itype,taxes,ietype,
     x    maxtyp,rcut,gdelt,e,grad,curv,nelements,iinclude)

c      write(6,*) 'vineyard: built hessian'

c  divide curvature elements by sqrt(mi*mj)
c  units:  hartree/(bohr*atomic units of mass) (i.e. - all atomic units)
      k=0
      trace=0.0d0
      twopi=2.0d0*3.1415926d0
      ii=0
      do iat=1,nmove
         if (iinclude(iat).eq.1) then
            ii=ii+1
            rmi=1822.83d0*amu(itype(iat))
            do ix=1,3
c               i=3*(iat-1)+ix
               i=3*(ii-1)+ix
               jj=0
               do jat=1,nmove
                  if (iinclude(jat).eq.1) then
                     jj=jj+1
                     rmj=1822.83d0*amu(itype(jat))
                     rmij=sqrt(rmi*rmj)
                     do jx=1,3
c                        j=3*(jat-1)+jx
                        j=3*(jj-1)+jx
                        if(j.gt.i) go to 10
                        k=k+1
                        curv(k)=(curv(k)/rmij)*0.5291d0**2
                        if(i.eq.j)
     +                      trace=trace+curv(k)/(2.418d-17*twopi)**2
 10                     continue
                     enddo
                  endif
               enddo
            enddo
         endif
      enddo

c      if (k.ne.nelements) write(6,*) 'nelements, k: ',nelements, k

c  diagonalize mass-weighted curvature matrix
c      n=nmove*3
      n=nelements*3
      nm=n
      lencur=m3v*(m3v+1)/2
      if(ivec.gt.0) then
        call wrkclaimv(1,1,'nmodesxyz')
        ierr=0
        call rsp(nm,n,lencur,curv,work2,1,eigvec,work3,work4,ierr)
        if(ierr.ne.0) write(6,*)'nmodesxyz - warning: rsp failed***'
      else if(ivec.eq.0) then
        call sdiagf(n,curv,work2,work3)
      end if

c  convert diagonal elements to frequencies in sec**-1
      twopi=2.0d0*3.1415926d0
      flogs=0.0d0
      nprod=0
      nneg=0
      do i=1,nelements*3
         dag=work2(i)
c         write(6,*) 'debug: dag(i): ',work2(i),i
         if(dag.le.0.0d0) nneg=nneg+1
         dag=sqrt(abs(dag))*sign(1.0d0,dag)/twopi
         work2(i)=dag/2.418d-17
         if(work2(i).gt.0.0d0) then !!!! xxx note .ge. is now .gt. 
            flogs=flogs + log(work2(i))
            nprod=nprod+1
         end if
      enddo

      freqsquaresum=0.0d0
      do i=1,nelements*3
         freqsquaresum = freqsquaresum + work2(i)*abs(work2(i)) 
      enddo


c compute gradient projection along lowest mode
c      if(ivec.gt.0) then
c         proj1=0.0d0
c         do i=1,nmove*3
c            proj1=proj1+grad(i)*eigvec(i)
c         enddo
c         proj1=proj1*27.21d0
c      end if

c  print out frequencies
      if(ip.ge.1) then
         write(6,80) (work2(i),i=1,nmove*3)
 80      format(/'  Normal mode frequencies (Hertz):'/
     x       (1p5d13.3))
      end if

c  print out eigenvectors (normal modes)
c      do 100 i=1,nmove*3
      do 100 i=1,99
      j=nelements*3*(i-1)
      if(ip.ge.2) write(7,90) i,work2(i),(eigvec(j+k),k=1,nmove*3)
      if(ip.ge.0) write(filnam,92) i
      if(ip.ge.0)
     +    open (11,file=filnam,form='formatted',status='unknown')
      if(ip.ge.0) write(11,91) nmove,(eigvec(j+k),k=1,nmove*3)  ! xxx blas - I made this contingent on ip
      if(ip.ge.0) close(11)
   90 format(' normal mode',i4,' frequency =',1pd13.3,' Hz'
     x            / (1x,0p3f10.4,2x,3f10.4,2x,3f10.4))
 91    format(i10/(0p3f20.10))
 92    format('eig.',i2)
  100 continue
      if(ip.ge.0) write(7,105) nneg,flogs,exp(flogs/float(nprod))
  105 format(' neg. curvs =',i3,'  ln(pos freqs) =',1pd20.10/
     x       '        geometric ave of pos freqs =',1pd15.5,' Hz')
      if(ip.ge.0.and.nneg.eq.0 .and. ivec.gt.0) write(6,107) 'pos',proj1
      if(ip.ge.0.and.nneg.ge.1 .and. ivec.gt.0) write(6,107) 'neg',proj1
  107 format(' grad projection onto lowest mode ('
     x                    ,a3,' curv) =',1pd12.3,' eV/A')
      if(flogsl.ne.0.0d0 .and. abs(nprod-nprodl).eq.1 .and. ip.ge.0)
     x     write(7,120) exp(flogs-flogsl),nprod-nprodl
  120 format(' Vineyard ratio of freqs product with previous product =',
     x     1pd13.3,' Hz**',i2)
      if(flogsl.ne.0.0d0 .and. abs(nprod-nprodl).eq.0 .and. ip.ge.0)
     x     write(7,121) exp(-(flogs-flogsl))
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

      call flushtad(7)

      return
      end


      subroutine cgcalcsubset(natom,nmove,xyz,itype,taxes,ietype,
     x    maxtyp,rcut,gdelt,e,grad,curv,nelements,iinclude)
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
c      parameter (maxmove=10000)
      include 'parameters.h'
      
      include 'wrksubs_tad2.h'
c      parameter(lwrk8=100 000, maxmove=10000) ! maxmove consistent with movmax above
      common/wrkcm8/iwrk8,jwrk8,gp(3,maxmove),gm(3,maxmove),
     x                     wrk8(lwrk8-6*maxmove)
      dimension iinclude(nmove)

      call wrkclaimv(8,1,'cgcalcn')
      if(natom.gt.maxmove) stop 'cgcaln - maxmove exceeded'

c compute gradient
c      write(6,*) 'debugA'
      call gcalc(natom,nmove,xyz,itype,taxes,ietype,maxtyp,rcut,e,grad
     +    ,hyperratio)

c zero hessian
      n=nelements*3
      do ij=1,n*(n+1)/2
         curv(ij)=0.0d0
      enddo

c compute hessian numerically using repeated gradient calls
      delt=gdelt

      ii=0
      do i=1,nmove
         if (iinclude(i).eq.1) then
            ii=ii+1
            do ix=1,3
c               idx=3*(i-1)+ix
               idx=3*(ii-1)+ix
               
               xi0=xyz(ix,i)

               xyz(ix,i)=xi0+delt
c               write(6,*) 'debugB'
               call gcalc(natom,nmove,xyz,itype, taxes,ietype,maxtyp
     +             ,rcut,etemp,gp, hyperratio)

               xyz(ix,i)=xi0-delt
c               write(6,*) 'debugC'
               call gcalc(natom,nmove,xyz,itype, taxes,ietype,maxtyp
     +             ,rcut,etemp,gm, hyperratio)

               xyz(ix,i)=xi0
               
               jj=0
               do j=1,nmove
                  if (iinclude(j).eq.1) then
                     jj=jj+1
                     do jx=1,3
c                        jdx=3*(j-1)+jx
                        jdx=3*(jj-1)+jx
                        ijx=idx*(idx-1)/2+jdx
                        if(jdx.gt.idx) ijx=jdx*(jdx-1)/2+idx
                        
                        hess = (gp(jx,j)-gm(jx,j))/(2.d0*delt)
                        if(jdx.eq.idx) then
                           curv(ijx)=hess
                        else
                           curv(ijx)=curv(ijx) + 0.5d0*hess
                        end if
c                        write(6,*) 'curv(ijx): ',ijx,curv(ijx),idx,jdx
                     enddo
                  endif
               enddo
            enddo
         endif
      enddo

      call wrkrelease(8)

      return
      end


      subroutine removeelements(natom,nmove,xyz,xyz2,nelements,lwrk13,
     +    iinclude)

      implicit real*8(a-h,o-z)
      
      dimension xyz(3,natom),xyz2(3,natom),iinclude(nmove)

      ndoable=int(1.0d0/6.0d0*(-1.0d0+sqrt(1.0d0+8.0d0*lwrk13)))-1

      nremove=nelements-ndoable
      
      do i=1,nremove
         rrmin=10d99
         iimin=0
         do j=1,nmove
            if (iinclude(j).eq.1) then
               rr=0.0d0
               do k=1,3
                  rr=rr+(xyz(k,j)-xyz2(k,j))**2
               enddo
               if (sqrt(rr).lt.rrmin) then
                  rrmin=sqrt(rr)
                  iimin=j
               endif
            endif
         enddo
         iinclude(j)=0
         nelements=nelements-1
      enddo

      return
      end
