c --------------------nik's conj. grad. code---------------

      subroutine acgstep(natom,nmove,xyz,itype,taxes,ietype,maxtyp,rcut
     $     ,hyperratio,ein,eout,dE,gmag,initcg,cgfac,a,g,h)
c     if iop=0 conjugate gradient minimization
c     if iop=1 gradient minimization
      implicit real*8(a-h,o-z)
c      parameter (mxatom= 10000)
      include'parameters.h'
      dimension rcut(maxtyp,maxtyp),itype(natom),xyz(3*natom),taxes(3)

c g1 is a work array for cg, but is needed to save info between calls
c am using the v array otherwise used for quickmin in this case
      dimension a(3*natom),g(3*natom),h(3*natom)

c xyzsave is to store a copy of xyz incase we need it
      dimension xyzsave(3*mxatom)
      save fsum,iter
c NOTE: a,g,h are based on gradient. In Kai's code, on the force.
c       careful with signs!

c      write(6,*) 'acg enter: a1: ',a(1),a(2),a(3)

      if (initcg.eq.1) then
         call gcalc(natom,nmove,xyz,itype, taxes,ietype,maxtyp,rcut,
     x       ein,a, hyperratio)
         do i=1,3*nmove
            a(i)=-a(i)
            g(i)=a(i)
            h(i)=a(i)
         enddo
         eout=ein
         fsum=0.0d0
         iter=0
         return
      endif

      nredo=0
      fprev=cgfac
      isteepest=0
      do i=1,3*nmove
         xyzsave(i)=xyz(i)
      enddo
 
 100  continue
      iter=iter+1
      do i=1,3*nmove
         a(i)=a(i)*cgfac
         xyz(i)=xyz(i)+a(i)
      enddo

      call gcalc(natom,nmove,xyz,itype, taxes,ietype,maxtyp,rcut,
     x    e,a,hyperratio)

      gmag=0.0d0
      do i=1,3*nmove
         a(i)=-a(i) ! force, not gradient
         gmag=gmag+a(i)*a(i)
      enddo
      gmag=sqrt(gmag)

c      write(6,*) 'acg: energy,eprev,de,isteepest: ',e,ein,e-ein
c     +    ,isteepest,gmag

      if (e.gt.ein) then
         nredo=nredo+1
         do i=1,3*nmove
            xyz(i)=xyzsave(i)
            if (nredo.lt.3) a(i)=h(i)
            if (nredo.ge.3) a(i)=g(i)
         enddo
         if (nredo.ge.10) then
            cgfac=0.0d0
            isteepest=1
c            if (nredo.eq.10) cgfac=fprev
c            cgfac=-cgfac
c            if (nredo.gt.50) cgfac=0.0d0
c            write(6,*) 'Emergency cgfac reselection: cgfac,prev,nredo: '
c     +          ,cgfac,fprev,nredo
         endif
         cgfac=cgfac*0.5d0
c         write(6,*) 'acg: energy increased, redoing: ',e-ein,cgfac
         goto 100
      endif
      cgfac=abs(cgfac)

      eout=e

      fsum=fsum+cgfac
      fmean=fsum/float(iter)
      fact=0.1d0+fmean/cgfac
      if (fact.lt.1.05d0) fact=1.05d0
      if (fact.gt.1.2d0) fact=1.2d0
      cgfac=cgfac*fact
      
c      write(6,*) 'acg: e-ein,cgfac,fact: ',e-ein,cgfac,fact

      gg=0.0d0
      dgg=0.0d0
      do i=1,3*nmove
         gg=gg+g(i)*g(i)
         dgg=dgg+(-a(i)+g(i))*(-a(i))
      enddo
c      write(6,*) 'acg: gmag: ',gmag,sqrt(gg)
      gmag=sqrt(gg)
      if (gg.eq.0.0d0) stop 'problem with gg in acg'
      gam=dgg/gg
c      write(6,*) 'acg: gam,dgg,gg: ',gam,dgg,gg
      do i=1,3*nmove
         g(i)=a(i)
         h(i)=g(i)+gam*h(i)
         if (isteepest.eq.0) a(i)=h(i)
      enddo
      if (isteepest.eq.1) cgfac=fprev*0.5d0
c      if (isteepest.eq.1) write(6,*) 'isteepest=1'

c      write(6,*) 'acg leave: a1: ',a(1),a(2),a(3)

      return
      end



c --------------------graeme's conj. grad. code---------------

      subroutine gcgstep(natom,nmove,xyz,itype,taxes,ietype,maxtyp,rcut
     $     ,hyperratio,ein,eout,dE,gmag,initcg,cgfac,a,g,h)
c     if iop=0 conjugate gradient minimization
c     if iop=1 gradient minimization
      implicit real*8(a-h,o-z)
c      parameter (mxatom= 10000)
      include 'parameters.h'
      dimension rcut(maxtyp,maxtyp),itype(natom),xyz(3*natom),taxes(3)

      dimension a(3*natom),g(3*natom),h(3*natom)

c xyzsave is to store a copy of xyz incase we need it
      dimension xyzsave(3*mxatom)
      save fsum,iter

      if (initcg.eq.1) then
         call gcalc(natom,nmove,xyz,itype, taxes,ietype,maxtyp,rcut,
     x       eout,a, hyperratio)
         do i=1,3*nmove
            a(i)=-a(i)          ! force, not gradient
            g(i)=a(i)
         enddo
         ein=eout
         fsum=0.0d0
         iter=0
         return
      endif

c1. Evaluate the force F0 at P0, the initial point.
c2. Evaluate the force F1 at P1=P0+UnitVec(F0)*dx.  
c  (In subsequent iterations, use the conjugate direction instead of F0).
c3. Move DeltaX=((F0+F1)/2)/((F0-F1)/dx)+dx/2 from the point P0.
c  (or move DeltaX=((F0+F1)/2)/((F0-F1)/dx)-dx/2 from point P1).
c4. Evaluate the conjugate direction and repeat.
      
      iter=iter+1
      amag=0.0d0
      gam=0.0d0
      ein=eout
      do i=1,3*nmove
         amag=amag+a(i)**2 ! a should be conjugate direction
      enddo
      amag=sqrt(amag)

      frp1=0.0d0
      do i=1,3*nmove
         frp1=frp1+g(i)*a(i)/amag
         xyzsave(i)=xyz(i)
         xyz(i)=xyz(i)+cgfac*a(i)/amag ! finite dif step to point P1
      enddo

      call gcalc(natom,nmove,xyz,itype, taxes,ietype,maxtyp,rcut,
     x    eout,h,hyperratio)       ! h=F1=F(P1)

      write(6,*) 'gcg: iter,energy: ',iter,e

      gmag=0.0d0
      frp2=0.0d0
      do i=1,3*nmove
         h(i)=-h(i) ! force, not gradient
         gmag=gmag+h(i)**2
         frp2=frp2+h(i)*a(i)/amag
      enddo
      gmag=sqrt(gmag)

      cr=(frp1-frp2)/cgfac
      frp=(frp1+frp2)/2.0d0
      rmaxmove=0.1d0
      if (cr.lt.0) then
         rdr=rmaxmove !maximum step size
      else
         rdr=frp/cr
         if (abs(rdr).gt.rmaxmove) rdr=sign(rmaxmove,rdr)
         rdr=rdr-cgfac/2.0d0
      endif

      do i=1,3*nmove
         xyz(i)=xyz(i)+a(i)*rdr/amag ! line minimization step
      enddo

      gg=0.0d0
      dgg=0.0d0
      h1=0.0d0
      do i=1,3*nmove
         gg=gg+g(i)**2          ! old f (F0)
         dgg=dgg+(h(i)-g(i))*(h(i)) ! h=F1, g=F0
         h1=h1+g(i)*h(i)
      enddo
      if (gg.eq.0.0d0) stop 'problem with gg in acg'
      gam=dgg/gg
      if (h1.le.0.5d0*gg) then  ! if force reduced, take CG step
         gam=dgg/gg
      else
         gam=0.0d0              ! otherwise, SD step
      endif
      write(6,*) 'gcg: gam,gg,dgg,gmag: ',gam,sqrt(gg),dgg,gmag
      do i=1,3*nmove
         g(i)=h(i)              ! save force of P1 for next step
         a(i)=g(i)+gam*a(i)     ! conjugate direction
      enddo

      return
      end

