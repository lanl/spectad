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

c 12/15/03: small bug in quickmin, wasn't doing vdotf over all components

      subroutine descent_check(natom,nmove,nstate,xin,xstate,ityp,
     x    taxes,ietype,maxtyp,rcut,hyperratio,barrevev,ereverse,isame,
     x     irotalign,itermax,gfac,transcrit,itranscrit,gcrit,dvcrit
     +     ,drcrit,ipr,imode,ifun)

c we need to update the comment cards in here

c performs a steepest descent on x (without changing it)
c to decide (as quickly as possible) which, if any, of
c the nstate different minima stored in xstate we are falling toward.
c Does a maximum of itermax steps in trying to reach convergence
c or a decision.

c xxx perhaps change "isame" to "imatch" since we are comparing to many states now

c output:
c isame=1 - it is the same
c isame=0 - it is not the same
c  x = replaced with partially (or completely) minimized version

c spring quench:
c the spring quenching is used to get a floppy system with long
c wavelength modes into alignment with a reference system.  atoms not
c involved in the transition have springs pulling them towards the
c reference state.  findzone is a routine to tell what atoms to not
c apply springs to.
c imode is to tell if a quench using springs to the reference state
c should be done.  if so, a second quench with the springs off is done
c immediately afterwards. added 020417
c imode=2 meands that dimer called descentcheck, so don't call dimer here
c Oct. 2008, cwp added ifunc = 7 using LBFGS minimization
      
      use mod_mpi
      implicit real*8(a-h,o-z)

c      parameter(mxatom=10000,nstatemax=1500,nneighmax=1500)
c      parameter(msave=7,nwork=3*mxatom*(2*msave+1)+2*msave) ! cwp for LBFGS
      include 'parameters.h'
      
      dimension xin(3*natom),xstate(3*natom*nstate),ityp(natom)
      dimension xtemp(3*mxatom),xtemp2(3*mxatom),rdimer(3*mxatom)
      dimension rcut(maxtyp,maxtyp)
c v is a work array for either quickmin or conjugate gradients
c gacg and hacg are work arrays for Kai's adaptive cg
      dimension grad(3*mxatom),v(3*mxatom),gacg(3*mxatom),hacg(3*mxatom)
      dimension gradold(3*mxatom)
      dimension xstatecm(3*nstatemax),xincm(3)
      dimension taxes(3)
      dimension icoordstate(mxatom*nstatemax),icoord(mxatom)
     +    ,izone(mxatom)
      dimension ineigh(nneighmax),rneigh(nneighmax)
     +    ,ineighstate(nneighmax*nstatemax),rneighstate(nneighmax
     +    *nstatemax),noffset(nstatemax),barrevev(nneighmax)
      dimension xspringave(3)
      dimension step(3*natom) !cwp for FIRE
      dimension nbd(mxatom),iwa(3*mxatom),isave(44),dsave(29) ! cwp for LBFGS-B
      dimension wa(2*24*3*mxatom+4*3*mxatom+12*24*24+12*24) ! cwp for LBFGS-B, 24 is mmax
      dimension xl(3*mxatom),xu(3*mxatom) ! cwp for LBFGS-B, stands for l,u
      dimension xtemp3(3*mxatom) ! cwp for LBFGS-B debug
      logical xlsave(4) ! cwp for LBFGS-B
      character*180 filnam,ctitle(2)
      character*60  task,csave !cwp for LBFGS-B
      

      common/callcom/ngradcall
      drmaxlbfgsb=0.5d0 ! cwp for LBFGS-B to limit step size
      sec2au=1.0d0/(2.418d-17)  ! au per second
      bohr2ang=0.529d0          ! angstrom per bohr

      rmaxstep=0.1d0            ! largest step to take in steepest descent or quickmin

      if (itranscrit.eq.3) then
         rcrit=rcut(1,1)
         if(ietype.eq.68) then ! cwp for meam
             rcrit=rcrit*0.8
         endif
         if(rank_f.eq.0) write(*,*) 'rcrit=',rcrit
         delta=0.01*rcrit
         delta=transcrit      ! xxxx how about this?  8/21/04 -art
      endif

      tadcharge(:) = 0.d0
      if(parallelmode.ne.'spectad')then !! cant write this many files when using spectad on a large system

          filnam=trim(filstrtLP)//'descentdoublecheck.dat'
c cwp added the following to write config before each descentcheck
          iunit=25
          if(rank_f.eq.0)then
           open(unit=iunit,file=filnam,form='formatted',
     x              status='unknown')
           rewind iunit
          endif
          iwrt=-1
          ntitle=1
          ctitle(1)='check descentcheck everytime'
          iop=1  ! 1 gives a compact, no-momentum file   - new as of 9/24/02
          if(rank_f.eq.0) call putclsnew(natom,ntitle,
     x                         ctitle,xin,xin,ityp,taxes,iunit,iwrt,iop)
          if(rank_f.eq.0) close(iunit)

      endif

c ifun defines minimzation routine:
c 0=steepest descent
c 1=quickmin
c 2=conjugate gradients (from clsman)
c 3=adaptive conjugate gradients from Kai Nordlund
c 4=graeme's conjugate gradients scheme
c 5=thermal quench
c      ifun=1
c 7=LBFGS-B minimization cwp
c 6=FIRE minimization cwp
      ifunhold=ifun
      if (ifun.lt.0.or.ifun.gt.7) stop 'ifun not defined' !cwp was .gt.5

      if (nstate.gt.nstatemax.and.itranscrit.eq.3) 
     +    stop 'Increase nstatemax!'
      
c      if(irotalign.ne.0.and.itranscrit.eq.0) stop
c     +    'need rotalign in descent_check'

      iz=0
      if (irotalign.eq.2) iz=1
      
      iusezone=0
      if (irotalign.ge.1.and.nstate.gt.0.and.imode.eq.1) iusezone=1
      if (iusezone.eq.1) ifun=0

      if ((ipr.ge.6).and.(rank_f.eq.0)) then
         if (ifun.eq.0) write(6,*) 'groupID ',groupID,
     +       ' descentcheck: using steepest descent'
         if (ifun.eq.1) write(6,*) 'descentcheck: using quickmin'
         if (ifun.eq.2) write(6,*)
     +       'descentcheck: using conjugate gradients'
         if (ifun.eq.3) write(6,*)
     +       'descentcheck: using adaptive conjugate gradients'
         if (ifun.eq.4) write(6,*)
     +       'descentcheck: using graemes conjugate gradient scheme'
         if (ifun.eq.5) write(6,*)
     +       'descentcheck: simulated annealing quench'
         if (ifun.eq.6) write(6,*)
     +       'descentcheck: FIRE quench'
         if (ifun.eq.7) write(6,*)
     +       'descentcheck: LBFGS quench'
         if (ifun.lt.0.or.ifun.gt.7) write(6,*)
     +       'descentcheck: unknown minimization function,'
     +       ,' expect crash!'
      endif

c make copy of input coordinates
      call vecmov(xin,xtemp,3*natom)
      
 1000 continue ! to do secondary quench if necessary
      if((ipr.ge.5).and.(rank_f.eq.0))
     +  write(6,*) 'groupID ',groupID,
     +  ' descentcheck: iusezone=',iusezone
      gfaclocal=gfac
      cgfac=gfac

      if (irotalign.eq.0.and.natom.eq.nmove.and.nstate.gt.0) then
         call shiftalign(xstate(1),xtemp,natom)
      endif
      if (itranscrit.eq.3) noffsetwork=1 ! maybe should be 0?
      do i=1,nstate
         call centerofmass_interface(xstate((i-1)*3*natom+1),natom
     +       ,xstatecm((i-1)*3+1))
         if (natom.ne.nmove) then
            do j=1,3
               xstatecm((i-1)*3+j)=0.0d0
            enddo
         endif
         if((ipr.ge.10).and.(rank_f.eq.0))
     x      write(6,*) 'CM, state ',i,': ',(xstatecm((i-1)*3+j),j=1,3)
         if (itranscrit.eq.1.or.itranscrit.eq.2) then
            call coordination(natom,xstate((i-1)*3*natom+1)
     +          ,icoordstate((i-1)*natom+1),transcrit,taxes,rcut
     +          ,ityp,itranscrit,maxtyp)
         endif
         if (itranscrit.eq.3) then
c            write(6,*) 'FB: noffset: ',noffsetwork
            noffset(i)=noffsetwork
            call fullbonding(natom,nmove,ineighstate
     +          ,rneighstate,rcrit,delta,nneighmax*nstatemax
     +          ,noffsetwork,xstate((i-1)*3*natom+1),taxes)
            if (noffsetwork.gt.nstatemax*nneighmax) then
              if(rank_f.eq.0)
     +          write(6,*) 'STOP: noffset.gt.nstatemax*nneighmax'
              stop 'noffset.gt.nstatemax*nneighmax'
            endif
         endif
      enddo
c      if (irotalign.eq.0.and.natom.eq.nmove.and.nstate.gt.0) then
c         call shiftalign(xstate(1),xtemp,natom)
c      endif
      if (itranscrit.eq.1.or.itranscrit.eq.2) then
         call coordination(natom,xtemp,icoord,transcrit,taxes
     +       ,rcut,ityp,itranscrit,maxtyp)
         if (iusezone.eq.1) then
            a=0.0d0
            b=0.0d0
            c=0.0d0
            call rotalign2(xtemp,natom,xstate(1),itermax,ipr
     +          ,iconverged,a,b,c,iz)
            call findzone(nmove,xtemp,xstate,icoord,icoordstate,izone
     +          ,xspringave,taxes)
         endif
      endif
      if (itranscrit.eq.3) then
         noffsetwork=1
c         write(6,*) 'FB here: ',noffsetwork
         call fullbonding(natom,nmove,ineigh,rneigh,rcrit,delta
     +       ,nneighmax,noffsetwork,xtemp,taxes)
      endif
      call centerofmass_interface(xtemp,natom,xincm)
      if (natom.ne.nmove) then
         do j=1,3
            xincm(j)=0.0d0
         enddo
      endif
      if(rank_f.eq.0)then
       if(ipr.ge.10) write(6,*) 'CM, xin',i,': ',(xincm(j),j=1,3)
      endif
      
      isame=0
      energy=0.0d0
      eprev=energy
  
c first set up call for each minimization method
      if (ifun.eq.0.or.ifun.eq.1.or.ifun.eq.6) then
c ifun=0: steepest decent
c ifun=1: quickmin (global)     
c ifun=6: FIRE    
         call gcalc(natom,nmove,xtemp,ityp,
     $       taxes,ietype, maxtyp,rcut, energy,grad, hyperratio)
         if (ifun.eq.1.or.ifun.eq.6) call dmzer(v,3*nmove)
         if (ifun.eq.6) then ! cwp for FIRE
             nmin2=5
             finc=1.1
             fdec=0.5
             alphastart=0.1
             alpha=alphastart
             falpha=0.99
             deltat=0.1
             deltatmax=deltat*10
             rmaxjump=0.2
         endif
      else if (ifun.eq.2) then
c ifun=2: conjugate gradient (per NR)
         initcg=1
         call cgstep(natom,nmove,xtemp,ityp,taxes,ietype,maxtyp
     $       ,rcut,hyperratio,ein,energy,dE,gmag,initcg,cgfac,v)
      else if (ifun.eq.3) then
c ifun=3: Kai Nordlund's adaptive conjugate gradients
         initcg=1
         call acgstep(natom,nmove,xtemp,ityp,taxes,ietype,maxtyp
     $       ,rcut,hyperratio,ein,energy,dE,gmag,initcg,cgfac,v,gacg
     +       ,hacg)
      else if (ifun.eq.4) then
c ifun=4: Graeme Henkelman's conjugate gradients
         initcg=1
         call gcgstep(natom,nmove,xtemp,ityp,taxes,ietype,maxtyp,rcut
     +       ,hyperratio,ein,energy,dE,gmag,initcg,cgfac,v,gacg,hacg)
      else if (ifun.eq.5) then
c ifun=5: damped dynamics minimzation
         call mdreset(1)
         call gcalc(natom,nmove,xtemp,ityp,
     $       taxes,ietype, maxtyp,rcut, energy,grad, hyperratio)
         call dmzer(v,3*nmove)
         icold=itermax/2.5
         tqtempnow=1d10
c ifun=7: LBFGS-B minimization cwp
      else if (ifun.eq.7) then
         call gcalc(natom,nmove,xtemp,ityp,
     $       taxes,ietype, maxtyp,rcut, energy,grad, hyperratio)
         
         mk=10 
         iprint2=0 ! <0 no output, =0 only first/last iter
         do i=1,3*nmove
            nbd(i)=2
         enddo
         do i=3*nmove+1,3*natom !cwp all atoms frozen during LBFGS-B
            nbd(i)=2
            xl(i)=xin(i)
            xu(i)=xin(i)
         enddo
         task= 'START'
         factr=0
         pgtol=0
      endif

      iconverged=0
c loop over iterations
      do i=1,itermax
         iter=i

         call vecmov(xtemp,xtemp2,3*nmove)
         call vecmov(grad,gradold,3*nmove)
         if (ifun.eq.0.or.ifun.eq.1.or.
     +       (ifun.eq.5.and.tqtempnow.lt.10.0d0)) then
            gmag=sqrt(vecdot(grad,grad,3*nmove))
            gmaxstepfact=1.0d0
            if (gfaclocal*gmag.gt.rmaxstep) then
               gmaxstepfact=1.0d0/gmag/gfaclocal*rmaxstep
            endif
            do j=1,3*nmove
               if (iusezone.eq.1) then
                  springk=0.005d0
                  if (izone(int(j/3)+1).eq.1) then
                     l=mod(j,3)+1
                     grad(j)=grad(j)-springk*xspringave(l)
                  else
                     grad(j)=grad(j)-springk*(xstate(j)-xtemp(j))
                  endif
               endif
               if (ifun.eq.0.or.ifun.eq.5) then ! steepest descent step
                  xtemp(j)=xtemp(j)-gfaclocal*grad(j)*gmaxstepfact
               endif
            enddo
            if (ifun.eq.1) then ! quickmin minimization
               vdotf=0.0d0
               fdotf=0.0d0
               do j=1,3*nmove
                  if (v(j)*grad(j).gt.0) v(j)=0.0d0
                  vdotf=vdotf-v(j)*grad(j)
                  fdotf=fdotf+grad(j)*grad(j)
               enddo
               vmag=0.0d0
               do j=1,nmove*3
                  v(j)=-grad(j)*(1.0d0+vdotf/fdotf)
c THIS IS WHAT GRAEME DID
c                  v(j)=-grad(j)*gfaclocal*(1.0d0+vdotf/fdotf/gfaclocal)
                  vmag=vmag+v(j)**2
               enddo
               vmag=sqrt(vmag)
               do j=1,nmove*3
                  if (vmag*gfaclocal.gt.rmaxstep) !max step size of rmaxstep
     +                v(j)=v(j)/vmag/gfaclocal*rmaxstep 
                  xtemp(j)=xtemp(j)+gfaclocal*v(j)
               enddo
            endif
            eprev=energy
            call gcalc(natom,nmove,xtemp,ityp,
     $          taxes,ietype, maxtyp,rcut, energy,grad, hyperratio)
c            write(6,*) 'energy report 1: ',energy,eprev,energy-eprev
         else if (ifun.eq.2) then
            initcg=0
            eprev=energy
            call cgstep(natom,nmove,xtemp,ityp,taxes,ietype,maxtyp
     $           ,rcut,hyperratio,ein,energy,dE,gmag,initcg,cgfac,v)
            gfaclocal=cgfac
         else if (ifun.eq.3) then
            initcg=0
            eprev=energy
            call acgstep(natom,nmove,xtemp,ityp,taxes,ietype,maxtyp,rcut
     +          ,hyperratio,eprev,energy,dE,gmag,initcg,cgfac,v,gacg
     +          ,hacg)
            gfaclocal=cgfac
         else if (ifun.eq.4) then
            initcg=0
            call gcgstep(natom,nmove,xtemp,ityp,taxes,ietype,maxtyp,rcut
     +          ,hyperratio,eprev,energy,dE,gmag,initcg,cgfac,v,gacg
     +          ,hacg)
            gfaclocal=cgfac
         else if (ifun.eq.5.and.tqtempnow.ge.10.0d0) then
            ifirst=0
            if (i.eq.1) ifirst=1
c            tempnow=(float(icold)-float(i))/float(icold)
            eprev=energy
            call dampedstep(natom,nmove,tqtempnow,xtemp,v,energy
     +          ,grad,rcut,maxtyp,taxes,ietype,ityp,ifirst)
         else if (ifun.eq.6) then ! cwp FIRE
            vdotf=0.0
            vmag=0.0
            fdotf=0.0
            do j=1,3*nmove
                  vdotf=vdotf-v(j)*grad(j) ! cwp get V dot F
                  vmag=vmag+v(j)*v(j) ! get |V|
                  fdotf=fdotf+grad(j)*grad(j) ! get |F|
            enddo
            
            vmag=sqrt(vmag)
            fdotf=sqrt(fdotf)
            if(vdotf.gt.0.0) then 
               do j=1,3*nmove
                   v(j)=(1-alpha)*v(j)-alpha*vmag*grad(j)/fdotf 
               enddo
               if(nvdotfgt0.ge.nmin2) then
                  deltat=min(deltat*finc,deltatmax)
                  alpha=alpha*fadec
               endif
               nvdotfgt0=nvdotfgt0+1
            else
               do j=1,3*nmove
                  v(j)=0.0
               enddo
               alpha=alphastart
               deltat=deltat*fdec
               nvdotfgt0=0
            endif

            stepmag=0.0
            do j=1,nmove*3
            v(j)=v(j)-deltat*grad(j)
            step(j)=deltat*v(j)
            stepmag=stepmag+step(j)*step(j)
            enddo
            
            stepmag=sqrt(stepmag)
            if(stepmag.gt.rmaxjump) then
               do j=1,3*nmove
                  step(j)=step(j)*rmaxjump/stepmag
               enddo
            endif
            do j=1,3*nmove !cwp update position and velocity
                  xtemp(j)=xtemp(j)+step(j)
            enddo
            eprev=energy
            call gcalc(natom,nmove,xtemp,ityp,
     $          taxes,ietype, maxtyp,rcut, energy,grad, hyperratio)

         else if (ifun.eq.7) then ! cwp LBFGS-B
 268     continue 
             call vecmov(xtemp,xtemp3,3*natom) ! cwp store current config in xtmp3
             do j=1,3*nmove
                xl(j)=xtemp(j)-0.2
                xu(j)=xtemp(j)+0.2
             enddo
             call setulb(3*natom,mk,xtemp,xl,xu,nbd,energy,grad,
     +            factr,pgtol,wa,iwa,task,iprint2,
     +            csave,lsave,isave,dsave)
             call maxdr(nmove,xtemp,xtemp3,xincm,xincm,drmax)
             if (drmax.ge.drmaxlbfgsb) then !cwp if a wild jump happens
                  ifun=0 !cwp switch to FIRE   
                  if(rank_f.eq.0) write(*,*) 'Wild jump happens!'
                  call vecmov(xtemp3,xtemp,3*natom) !cwp go back to old config
                  goto 1000
             endif
             if (task(1:8) .eq. 'ABNORMAL') then
                  if(rank_f.eq.0)
     +             write(*,*) 'Abnormal LBFGS minimization! use ifun=6'
                  call vecmov(xtemp3,xtemp,3*natom)
                  ifun=0
                  goto 1000
             endif
             if (task(1:11) .eq. 'CONVERGENCE') then
                  if(rank_f.eq.0) write(*,*)
     +             'Abnormal LBFGS minimization! use ifun=6'
                  call vecmov(xtemp3,xtemp,3*natom)
                  ifun=0
                  goto 1000
             endif
             if (task(1:2) .eq. 'FG') then
             eprev=energy ! perhaps add a displacement check here
             call gcalc(natom,nmove,xtemp,ityp,
     $       taxes,ietype, maxtyp,rcut, energy,grad, hyperratio)
             do j=3*nmove+1,3*natom
                 grad(j)=0
             enddo
             goto 268
             endif
             
         else
            if(rank_f.eq.0) write(6,*) 'ifun not recognized'
            stop 'ifun not recognized'
         endif
         dv=abs(energy-eprev)
         dvuse=energy-eprev
         if (ifun.ne.2.and.ifun.ne.3.and.ifun.ne.4.and.ifun.ne.7) !cwp
     +       gmag=sqrt(vecdot(grad,grad,3*nmove))
         if (ifun.ne.7) then !cwp
             call maxdr(nmove,xtemp,xtemp2,xincm,xincm,drmax)
         endif

         if (ifun.eq.7.and.task(1:5).eq.'NEW_X') then !cwp
             gmag=sqrt(vecdot(grad,grad,3*nmove))
             call maxdr(nmove,xtemp,xtemp2,xincm,xincm,drmax)
         endif

         if (ifun.eq.6) then !cwp
             gmag=sqrt(vecdot(grad,grad,3*nmove))
             call maxdr(nmove,xtemp,xtemp2,xincm,xincm,drmax)
         endif
c check if energy increased during step (moved here 021206 to fix
c convergence problems)
         ienergyhigher=0
c         write(6,*) 'energy report 2: ',energy,eprev,energy-eprev
c         write(6,269) xtemp(1),xtemp2(1),grad(1),gradold(1)
 269     format('check: ',4d20.5)
         if(ifun.ne.7.and.ifun.ne.6) then !cwp
         iadjustgfac=0
         if (iadjustgfac.eq.1) then
            if (energy.gt.eprev.and.abs(energy-eprev).gt.1.0e-7.and ! added 020620 by BPU to fix convergence problems
     +          .iusezone.eq.0.and.(ifun.ne.5.or.tqtempnow.lt.10.0d0)
     +          .and.iincrease.ne.1) then
               if (ipr.ge.10) write(6,259) i,energy,eprev,energy-eprev
 259           format('descentcheck: energy increased, iter: ',i5,1p3d30
     +             .16)
               ienergyhigher=1  ! do not break out of iteration loop if energy increased
               if (ifun.eq.0.or.ifun.eq.1.or.(ifun.eq.5.and.tqtempnow.lt
     +             .10.0d0)) then
                  gfaclocal=gfaclocal/2.0d0 
                  if (itranscrit.ne.0) then
c                  gfaclocal=gfaclocal/2.0d0 
                     energy=eprev
                     call vecmov(xtemp2,xtemp,3*nmove) ! energy went up, restore previous position (uuuu do for all ifun?)
                     call vecmov(gradold,grad,3*nmove)
                     gmag=sqrt(vecdot(grad,grad,3*nmove))
                  endif
               endif
               if (ifun.eq.1) then
                  call dmzer(v,3*nmove) ! zero v
               endif
               iincrease=1
            else if (energy.lt.eprev.and.
     +             (ifun.ne.5.or.tqtempnow.lt.10.0d0)) then
c            if (itranscrit.ne.0) 
               gfaclocal=gfaclocal*1.01d0
c            if (gfaclocal.gt.1.0d0) gfaclocal=1.0d0 ! xxxx temp ceiling
               iincrease=0
            else if (energy.eq.eprev) then
               iincrease=0
            endif
         endif
         endif !cwp
         isame=0
         isamerev=0
         istride = nstate/10 + 1   ! guessing that checking 10 states is roughly a force call's work
         if (mod(i,istride).eq.0.or.i.eq.1.or.iconverged.eq.1) then
            if (irotalign.eq.0.and.natom.eq.nmove.and.nstate.gt.0) then
               call shiftalign(xstate(1),xtemp,natom)
            endif
            if (irotalign.ne.0.and.itranscrit.eq.0.and.natom.eq.nmove
     +          .and.nstate.gt.0) then
               a=0.0d0
               b=0.0d0
               c=0.0d0
               call rotalign2(xtemp,natom,xstate(1),itermax,ipr
     +             ,iconverged,a,b,c,iz)
            endif
            if (itranscrit.eq.1.or.itranscrit.eq.2) then
               call coordination(natom,xtemp,icoord,transcrit,taxes
     +             ,rcut,ityp,itranscrit,maxtyp)
               if (iusezone.eq.1) then
                  call findzone(nmove,xtemp,xstate,icoord,icoordstate
     +                ,izone,xspringave,taxes)
               endif
            endif
            if (itranscrit.eq.3.and.nstate.gt.0) then
               noffsetwork=1
               call fullbonding(natom,nmove,ineigh,rneigh,rcrit,delta
     +             ,nneighmax,noffsetwork,xtemp,taxes)
            endif
            itwomatches=0
            do j=1,nstate
               if (itranscrit.eq.0) then
                  call maxdr(nmove,xtemp,xstate((j-1)*3*natom+1),xincm
     +                ,xstatecm((j-1)*3+1),deltarmax)
                  if ((mod(i,250)).eq.0.or.i.eq.1) then
                     if ((ipr.ge.10).and.(rank_f.eq.0)) write(6,*)
     +                   'descentcheck: maxdr: ',j,deltarmax,transcrit
                  endif
                  if (deltarmax.lt.transcrit) then ! potential match
c                     write(6,*) 'dc debug: ',j,deltarmax
                     if (j.gt.1) then
                        if (barrevev(j-1).ge.ereverse.and.isame.ne.0)
     +                      then
c go ahead and keep minimizing if more than one
c match, just stop if converged and still more than one
c                        if(isame.ne.0) then ! already have a match to something
c                           write(6,*)
c     +                         'descentcheck: state matches twostates '
c     +                         ,isame,' and ',j
c                           write(6,*)
c     +                         'descentcheck: state in twomatches.dat'
c                           filnam='twomatches.dat'
c                           call storefile(natom,xtemp,xtemp,ityp,taxes
c     +                         ,filnam)
c                           stop 'descentcheck: two state matches'
                           itwomatches=1
                        elseif (barrevev(j-1).ge.ereverse.and.isame.eq.0
     +                         ) then
                           isame=j
                        else
c                           if (ipr.ge.4) write(6,*)
c     +                         'dc debug: state matched ',j
c     +                      ,', but barrev small (',barrevev(j-1),')'
                           isamerev=j ! xxx overwrites previous matches to reverse barrier states, may not want
                        endif
                     elseif (j.eq.1) then
                        isame=j
                     endif
                  endif
               else if (itranscrit.eq.1.or.itranscrit.eq.2) then
                  icompare=0
                  do k=1,nmove
                     if (icoord(k).ne.icoordstate((j-1)*natom+k))
     +                   icompare=k
                  enddo
                  if (icompare.eq.0) isame=j
               else if (itranscrit.eq.3) then
                  icompare=0
                  noffsetwork=noffset(j)
c                  write(6,*) 'FBC: noffset: ',j,noffsetwork
                  call fullbondcheck(nmove,j,ineigh,rneigh
     +                ,ineighstate,rneighstate
     +                ,nstatemax,nneighmax,rcrit,delta,noffsetwork
     +                ,icompare)
 678              if (icompare.eq.0) isame=j
               else
                  stop 'itranscrit not recognized'
               endif
            enddo
         endif
         
         if(rank_f.eq.0)then
          if (ipr.ge.10) write(6,250) dvuse
     +          ,gmag,drmax,iter,isame,gfaclocal,rank_l
          if ((mod(i,250).eq.0.or.i.eq.1).and.ipr.ge.5) then
            write(6,250) dvuse
     +          ,gmag,drmax,iter,isame,gfaclocal,rank_l
            if ((ipr.ge.10).and.(rank_f.eq.0))
     +          write(6,*) 'ifun=5: tqtempnow: ',tqtempnow
          endif
         endif
         
         if (isame.ne.0) goto 210
         if(dv.lt.dvcrit.and.gmag.lt.gcrit.and.drmax.lt.drcrit
     +       .and.ienergyhigher.eq.0) then
            iconverged=iconverged+1
            if (iconverged.eq.2) goto 220
         endif
         
      enddo

      if(rank_f.eq.0) write(6,*)
     x           'descent_check warning - did not converge or find',
     x           'match before itermax'
      if (ipr.ge.1) then
       if(rank_f.eq.0) write(6,250)
     x                   dv,gmag,drmax,iter,isame,gfaclocal,rank_l
       call flushtad(6)
      endif

 220  continue                  ! convergence reached, xin doesn't match any xstate with reverse barrier > ereverse

      if (itwomatches.eq.1) then
         if(rank_f.eq.0) write(6,*)
     +       'descentcheck: state matches twostates '
     +       ,isame,' and ',j
         if(rank_f.eq.0) write(6,*)
     +       'descentcheck: state in twomatches.dat'
         filnam='twomatches.dat'
         call storefile(natom,xtemp,xtemp,ityp,taxes,filnam)
         stop 'descentcheck: two state matches'
      endif

      if (isamerev.ne.0) then   ! state matches another state with reverse barrier < ereverse
         if((ipr.ge.6).and.(rank_f.eq.0))
     +       write(6,*) 'dcdebug: isamerev=',isamerev
         isame=isamerev         ! TAD will know to ignore this as it will check for reverse barrier size
         goto 210
      endif

c following turned off for irotalign=0  9/1/04.  May want to also or instead
c key this on natom=nmove.  blas and art   xxx
c xxxx this dimer search is called when in a new basin and haven't
C matched anything yet.  this is a problem when using irecognize=2 and
C synthetic mode, as it does this before checking reordering, but wastes
C 100 force calls which dominate the synthetic mode immediate accept.
      if(parallelmode.ne.'spectad')then
      if ((imode.ne.2).and.(irotalign.eq.0).and.(natom.ne.nmove)) then ! not called by dimer
         idimerbail=0
         imodeds=3
         iranx=0
         irandir=1
         ngradcallhold=ngradcall
         idimerdof=3
c         filnam='before.dat'
c         call storefile(natom,xtemp,xtemp,ityp,taxes,filnam)
c xxxx ndimerdof, etc, not actually defined here, need to either set
c dummy or pass in
         call dimer_search(natom,nmove,idimerbail,emin,eminbar,idimerdof
     +       ,ndimerdof,rdimerdist,ndimercoord,idimeroverunder
     +       ,rdimerdisp,xtemp2,xtemp,xtemp2,xtemp2,energy,rdimer,ityp
     +       ,taxes,ietype,maxtyp,rcut,hyperratio,barrevev,ereverse,gfac
     +       ,irotalign,itermax,ierr,iranx,irandir,transcrit,itranscrit
     +       ,gcrit,dvcrit,drcrit,ipr,imodeds,ifun,eig)
         ngradcallds=ngradcall-ngradcallhold
         if ((ipr.ge.4).and.(rank_f.eq.0))
     +      write(6,*) 'descentcheck: dimer: ',ngradcallds,eig
         if (eig.lt.0.0d0) then
            if ((ipr.ge.4).and.(rank_f.eq.0))
     +       write(6,*) 'descentcheck: dimer found neg eig'
            do ids=1,3*nmove
               !xtemp(i)=xtemp(i)+0.01d0*(prngen(0)*2.0d0-1.0d0)
               if(rank_f.eq.0)then
                 xtemp(ids)=xtemp(ids)+0.01d0*(prngen(0)*2.0d0-1.0d0) ! RJZ -- BUG FOUND?
               endif
            enddo
            call MPI_BCAST(xtemp,3*nmove,MPI_REAL8,0,force_comm,ier)
            if(ifun.eq.7) ifun=0 !cwp switch to steepest descent
            goto 1000
         endif
      endif
      endif

c secondary quench: need to do if using springs to quench long
c -wavelength modes so that atoms are put in their actual minimized
c positions 
      if (iusezone.ne.0.and.(itranscrit.eq.1.or.itranscrit.eq.2).and
     +    .nstate.gt.0) then
         if((ipr.ge.5).and.(rank_f.eq.0)) write(6,250)
     +       dv,gmag,drmax,iter,isame,gfaclocal,rank_l
         if(rank_f.eq.0)
     +      write(6,*) 'descentcheck: doing secondary quench'
         iusezone=0
c         if(ifun.eq.7) ifun=0 ! cwp if need to do secondary quench use steepest descent
         goto 1000
      endif

      call vecmov(xtemp,xin,3*natom)
      if (itranscrit.eq.1.or.itranscrit.eq.2) then
         do j=1,natom
           if (icoord(j).ne.icoordstate(j).and.icoordstate(j).ne.0)then
             if(rank_f.eq.0) write(6,251) j,icoordstate(j),icoord(j)
           endif
         enddo
      endif
 251  format('descentcheck: atom ',i5,' changed coordination from'
     +    ,i5,' to ',i5)

 210  continue ! xin matches xstate(isame)
      if ((ipr.ge.5).and.(rank_f.eq.0))
     +       write(6,250) dv,gmag,drmax,iter,isame,gfaclocal,rank_l
      if (irotalign.eq.0.and.natom.eq.nmove.and.nstate.gt.0) then
         call shiftalign(xstate(1),xin,natom)
      endif
      if (irotalign.gt.0.and.nstate.gt.0) then
         a=0.0d0
         b=0.0d0
         c=0.0d0
         call rotalign2(xin,natom,xstate(1),itermax,ipr
     +       ,iconverged,a,b,c,iz)
      endif
 250  format('descentcheck: dv,gmag,drmax,iter,isame,rank_l='
     +    ,1p3d11.2,2i6,1p1d11.2,1i6)
   
      if (ifun.eq.5) call mdreset(1)

      ifun=ifunhold

      return
      end
c**********************************************************************

      subroutine maxdr(nmove,x1,x2,x1cm,x2cm,drmax)
      implicit real*8(a-h,o-z)
      dimension x1(3,nmove),x2(3,nmove)
      dimension x1cm(3),x2cm(3)

      drmax=0.0d0
      do j=1,nmove
         dr=0.0d0
         do k=1,3
            dr=dr+(x1(k,j)-x2(k,j))**2
         enddo
         if (drmax.lt.dr) jmax=j
         drmax=max(drmax,dr)
      enddo

      drmax=sqrt(drmax)
      return
      end
c**********************************************************************
      subroutine coordination(natom,x,icoord,transcrit,taxes,
     +    rtransarray,itype,itranscrit,maxtyp)

      implicit real*8(a-h,o-z)

      dimension icoord(natom),x(3*natom),taxes(3),itype(natom)
     +    ,rtransarray(maxtyp,maxtyp)

      do j=1,natom
         icoord(j)=0
         do k=1,natom
            rjk=0.0d0
            do l=1,3
               rr=x((j-1)*3+l)-x((k-1)*3+l)
               if (abs(rr).gt.taxes(l)/2.0d0) rr=rr-sign(taxes(l)
     +             ,rr)
               rjk=rjk+rr**2
            enddo
            if (itranscrit.eq.1) then
               if (k.ne.j.and.sqrt(rjk).lt.transcrit)
     +             icoord(j)=icoord(j)+1
            else if (itranscrit.eq.2) then
               if (k.ne.j.and.sqrt(rjk).lt.
     +             rtransarray(itype(k),itype(j))) icoord(j)=icoord(j)+1
            endif
         enddo
      enddo

      return
      end
c**********************************************************************
      subroutine findzone(n,x,xstate,icoord,icoordstate,izone,
     +    xspringave,taxes)
c find atoms that have changed coordination plus those neighbors within
c a radius
      
      implicit real*8(a-h,o-z)

      dimension x(3*n),xstate(3*n),icoord(n),icoordstate(n),izone(n)
      dimension xspringave(3),taxes(3)
      dimension ilist(100)

      nzone=0
      nlist=0
      do j=1,n
         izone(j)=0
         if (icoord(j).ne.icoordstate(j)) then
            nlist=nlist+1
            ilist(nlist)=j
         endif
      enddo

      nn=0
      do j=1,n
         rmin=1.0d10
         do k=1,nlist
            r=0.0d0
            do l=1,3
               rr=x((j-1)*3+l)-x((ilist(k)-1)*3+l)
               if (abs(rr).gt.taxes(l)/2.0d0) rr=rr-sign(taxes(l),rr)
               r=r+rr**2
            enddo
            rmin=min(rmin,r)
         enddo
         if (rmin.le.4.0d0) then ! radius is 2 (4=2**2)
            izone(j)=1
            nzone=nzone+1
c            write(6,*) 'fz4 debug: ',j,(x((j-1)*3+l),l=1,3)
c            write(6,*) 'fz4 debug: ',j,(xspringave(l),l=1,3),float(nn)
         else if (rmin.gt.4.0d0.and.rmin.le.9.0d0) then
            nn=nn+1
            do l=1,3
               if (nn.eq.1) xspringave(l)=0.0d0
               xspringave(l)=xspringave(l)+
     +             xstate((j-1)*3+l)-x((j-1)*3+l)
            enddo
c            write(6,*) 'fz debug: ',j,(xstate((j-1)*3+l)-x((j-1)*3+l),l
c     +          =1,3),nn
c            write(6,*) 'fz debug: ',j,(xstate((j-1)*3+l),l=1,3)
c            write(6,*) 'fz debug: ',j,(x((j-1)*3+l),l=1,3)
c            write(6,*) 'fz debug: ',j,(xspringave(l),l=1,3),float(nn)
c            write(6,*) 'fz debug: ',j,(xspringave(l)/float(nn),l=1,3)
         endif
      enddo

      do l=1,3
         if (nn.eq.0) then
            xspringave(l)=0.0d0
         else
            xspringave(l)=xspringave(l)/float(nn)
            xspringave(l)=0.0d0
         endif
      enddo
      
c      write(6,*) 'findzone: ',nlist,' atoms moved, ',nzone
c     +    ,' atoms in zone, ',nn,' in second shell'
c      write(6,*) 'findzone: xspringave: ',(xspringave(l),l=1,3)

      return
      end

c======================================================================
      subroutine shiftalign(x1,x2,natom)
c xxxx shiftalign does not use masses!
      implicit real*8(a-h,o-z)
      dimension x1(natom*3),x2(natom*3)

      dimension xcm1(3),xcm2(3)

      do j=1,3
         xcm1(j)=0.0d0
         xcm2(j)=0.0d0
      enddo
      
      do i=1,natom
         do j=1,3
            xcm1(j)=xcm1(j)+x1(3*(i-1)+j)
            xcm2(j)=xcm2(j)+x2(3*(i-1)+j)
         enddo
      enddo

      do j=1,3
         xcm1(j)=xcm1(j)/float(natom)
         xcm2(j)=xcm2(j)/float(natom)
      enddo

      do i=1,natom
         do j=1,3
            x2(3*(i-1)+j)=x2(3*(i-1)+j)+(xcm1(j)-xcm2(j))
         enddo
      enddo

      return
      end
      
c**********************************************************************
      subroutine dampedstep(natom,nmove,tqtempnow,xyz,pxyz
     +    ,energy,grad,rcut,maxtyp,taxes,ietype,itype,ifirst)

c this is a wrapper for mdstep for the damped step case
      use mod_mpi
      implicit real*8(a-h,o-z)
c      parameter (mxatom=10000)
      include 'parameters.h'
      dimension xyz(3*natom),pxyz(3*natom),grad(3*natom),itype(natom)
      dimension rcut(maxtyp,maxtyp),taxes(3)
      dimension moving(mxatom),mx(mxatom),my(mxatom),mz(mxatom)

      parameter(maxtherm=10)
      dimension itherm(maxtherm),jtherm(maxtherm)
      dimension ttherm(maxtherm),ctherm(maxtherm)
      dimension xtherm(maxtherm),mtherm(maxtherm),ktherm(maxtherm)

      common/bpuds/dt,amu(5),tqtempmax,tqrate,ntype ! ideally passed, but common for now (BPU type here)
      
c      write(6,*) 'dampedstep: ',natom,nmove,ntype
c      write(6,*) 'dampedstep: ',rcut(1,1),maxtyp,taxes(1),ietype

c specify moving atoms (ideally passed, but hardwire for now)
      do i=1,nmove
         moving(i)=1
         mx(i)=1
         my(i)=1
         mz(i)=1
      enddo
      do i=nmove+1,natom
         moving(i)=0 
         mx(i)=0
         my(i)=0
         mz(i)=0
      enddo
      time=0.0d0
      hyperratio=0.0d0

c thermostat parameters (ideally passed (except for tempnow), but
C hardwire for now)
c      tqtempmax=300.d0
      if (ifirst.eq.1) then
         tqtempnow=tqtempmax
      else
         tqtempnow=tqtempnow-tqrate
      endif
c      write(6,*) 'ifun5 debug: ',ifirst,tqtempnow,tqtempmax

      ntherm=1
      itherm(1)=1
      jtherm(1)=nmove
      ctherm(1)=1.0d12
      ttherm(1)=tqtempnow
      mtherm(1)=0
      ktherm(1)=4

      temptemp=ttherm(1)

      if (ifirst.eq.1) then
         iop=0
         iwrt=1
         call gaussp(natom,nmove,pxyz,itype,maxtyp,amu,ntype,moving,mx
     +       ,my,mz,tqtempnow,iop,iwrt)
         call MPI_BCAST(pxyz,3*nmove,MPI_REAL8,0,force_comm,ier)
      endif

      call mdreset(1)

c      write(6,*) 'dampedstep: temperature: ',ttherm(1),temptemp
c      write(6,*) 'dampedstep: amu: ',amu(1),amu(2),amu(3)

      call mdstep(natom,nmove,xyz,itype,pxyz, taxes,
     x    ietype, maxtyp,rcut,amu,
     x    maxtherm,itherm,jtherm,ttherm,ctherm,
     x    xtherm,mtherm,ktherm,ntherm,
     x    moving,mx,my,mz, dt,time,thyper,hyperratio, energy,grad)

      return
      end

c**********************************************************************
      subroutine fullbonding(natom,nmove,ineigh,rneigh,rcrit,delta
     +    ,nneighmax,noffset,x,taxes)

      implicit real*8(a-h,o-z)

      dimension ineigh(nneighmax),rneigh(nneighmax)
      dimension x(3*natom),taxes(3)

c noffset is a parameter to locate information in ineigh,rneigh arrays
c when there are many states in the array

c      write(6,*) 'entered fullbonding: noffset: ', noffset

      do j=1,nmove
         nneigh=0
         do k=1,natom
c            write(6,*) 'FB: ',natom,nmove,nneigh,noffset
            if (k.ne.j) then
               rjk=0.0d0
               do l=1,3
                  rr=x((j-1)*3+l)-x((k-1)*3+l)
                  if (abs(rr).gt.taxes(l)/2.0d0) rr=rr-sign(taxes(l)
     +                ,rr)
                  rjk=rjk+rr**2
               enddo
               rjk=sqrt(rjk)
               if (rjk.le.rcrit+delta) then
c                  write(6,*) 'FB: ',j,k,rjk,nneigh+1+noffset
                  nneigh=nneigh+1
                  ineigh(nneigh+noffset)=k
                  rneigh(nneigh+noffset)=rjk
               endif
            endif
         enddo
         ineigh(noffset)=nneigh
         rneigh(noffset)=0.0d0
         noffset=noffset+nneigh+1
      enddo
c      write(6,*) 'FB end: noffset: ',noffset

      return
      end

cc**********************************************************************
c      subroutine fullbondcheckbad(nmove,j,ineigh,rneigh,ineighstate
c     +    ,rneighstate,nstatemax,nneighmax,rcrit,delta,noffset,icompare)
c
c      implicit real*8(a-h,o-z)
c
c      dimension ineigh(nneighmax),rneigh(nneighmax)
c     +    ,ineighstate(nstatemax*nneighmax),rneighstate(nstatemax
c     +    *nneighmax)
c
cc      write(6,*) 'entered fullbondcheck'
cc      write(6,*) 'FBC: ',noffset
c
c      nindex=1
c      nindexj=noffset
cc      write(6,*) 'FBC: j: ',j,nindexj
c      do k=1,nmove
c         nneigh=ineigh(nindex)
c         nneighj=ineighstate(nindexj)
cc         write(6,*) 'FBC: nneigh,nneighj: ',nneigh,nneighj
c         nmin=min(nneigh,nneighj)
c         nmax=max(nneigh,nneighj)
c         do l=1,nmin
cc            write(6,*) 'FBC: ',j,k,ineigh(nindex+l),ineighstate(nindexj
cc     +          +l),rneigh(nindex+l),rneighstate(nindexj+l)
c            if (ineigh(nindex+l).eq.ineighstate(nindexj+l))then
c               rdelbond=rneigh(nindex+l)-rneighstate(nindexj+l)
c               if (abs(rdelbond).gt.delta) then
c                  icompare=k
c                  return
c               endif
c            else
c               icompare=k
c               return
c            endif
c         enddo
c         if (nneigh.ne.nneighj) then
c            iwhich=0
c            if (nneigh.gt.nneighj) iwhich=1
c            do l=nmin+1,nmax
c               if (iwhich.eq.1) then
c                  if (rneigh(nindex+l).le.rcrit) then
c                     icompare=k
c                     return
c                  endif
c               else
c                  if (rneighstate(nindexj+l).le.rcrit) then
c                     icompare=k
c                     return
c                  endif
c               endif
c            enddo
c         endif
c         nindex=nindex+nneigh+1
c         nindexj=nindexj+nneighj+1
c      enddo
c      noffset=nindexj
c      
c      return
c      end
c**********************************************************************
      subroutine fullbondcheck(nmove,j,ineigh,rneigh,ineighstate
     +    ,rneighstate,nstatemax,nneighmax,rcrit,delta,noffset,icompare)

c bpu: 08/05/03: found bug (only checked first bond for an atom if the
C neighbor was the same for both states), fixed

      implicit real*8(a-h,o-z)

      dimension ineigh(nneighmax),rneigh(nneighmax)
     +    ,ineighstate(nstatemax*nneighmax),rneighstate(nstatemax
     +    *nneighmax)

      nindex=1
      nindexj=noffset
      do k=1,nmove
         nneigh=ineigh(nindex)
         nneighj=ineighstate(nindexj)
         l=1
         lj=1
         nmin=min(nneigh,nneighj)
         nmax=max(nneigh,nneighj)
 567     continue
         iineigh=ineigh(nindex+l)
         iineighj=ineighstate(nindexj+lj)

c         write(6,*) 'FBC: j,k,in,inj,l,lj,nneigh,nneighj: ',j,k,iineigh
c     +       ,iineighj,l,lj,nneigh,nneighj

         if (iineigh.eq.iineighj)then
            rdelbond=rneigh(nindex+l)-rneighstate(nindexj+lj)
            if (abs(rdelbond).gt.delta) then
c               write(6,*) 'FBC: rdelbond: ',j,k,iineigh,iineighj
c     +             ,rdelbond
               icompare=k
               return
            else
               l=l+1
               lj=lj+1
               if (l.le.nneigh.and.lj.le.nneighj) goto 567
            endif
         else
            if (iineigh.lt.iineighj) then
               if (rneigh(nindex+l).le.rcrit) then
c                  write(6,*) 'FBC: iineigh bond small: ',j,k,iineigh
c     +                ,rneigh(nindex+l)
                  icompare=k
                  return
               else
                  l=l+1
                  if (l.le.nneigh) then
                     goto 567
                  endif
               endif
            else if (iineigh.gt.iineighj) then
               if (rneighstate(nindexj+lj).le.rcrit) then
c                  write(6,*) 'FBC: iineighj bond small: ',j,k,iineighj
c     +                ,rneighstate(nindexj+lj)
                  icompare=k
                  return
               else
                  lj=lj+1
                  if (lj.le.nneighj) then
                     goto 567
                  endif
               endif
            endif
         endif
         nindex=nindex+nneigh+1 ! xxxx need to check if too big (vs nneighmax)
         nindexj=nindexj+nneighj+1
      enddo
      
      return
      end
