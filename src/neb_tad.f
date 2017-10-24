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

c HISTORY:
c v1: * basic neb with various checks to make sure we get the saddle
c     connecting the initial state with the very next basin
c v2: * adds an option to do NEB for some number of steps and then launch
c     a dimer from the highest energy state to more quickly find the
c     saddle
c     * also: add in hot trajectory point to weight initial band
c     towards hot trajectory

c afv - changed imixhot integer to "hotmix", allowing control of fraction
c bpu (020827): added inebirc parameter to not do roll check (not fully
c               implemented yet)
c cwp (11052008): add ifunneb = 2 (FIRE) 
c cwp (11062008): add ifunneb = 3 (LBFGS)

      subroutine saddle_find(natom,x1,x2,xhot,hotmix,xsneb,es,nmove
     +   ,ityp,taxes,ietype,maxtyp,rcut,amu,hyperratio,gfac,irotalign
     +   ,itermax,ierr,ididintermediate,transcrit,itranscrit,gcritin
     +   ,dvcritin,drcritin,ipr,nimage,eshallow,xstates,nstate,barrevev
     +   ,ereverse,itan,iclimb,springk,lstore_neb,filnam1,filnam2
     +   ,idimerneb,inebirc,ifundc,ifunneb,imodecheck,rmodemag
     +   ,nnegsad,nprod,freqlogssad,freqlogsmin,nnegmin,prefacsad
     +   ,ivineyard,ntype,potnam,e0)
      
c ierr codes: ierr=0: no error
c             ierr=1: roll check doesn't connect correct saddles
c             ierr=2: no real saddle in chain, further convergence
c                     reveals no real transition      
c             ierr=3: shallow (neglected) minimum found, rollcheck will
c                     fail
c             ierr=5: NEB didn't converge
      use mod_mpi
      implicit real*8(a-h,o-z)
c      parameter(mxatom=10000)
      include 'parameters.h'
      dimension ityp(natom)
      dimension rcut(maxtyp,maxtyp)
      dimension amu(maxtyp)
      dimension taxes(3)
      character*180 filnam1,filnam2,filnam,filnamsad
      character*10 chri
      character*80 potnam(maxtyp)
      logical lstore_neb,lnoise

      dimension x1(3*natom),x2(3*natom),xsneb(3*natom),xhot(3*natom)
      dimension grad(3*natom), q(natom)
      dimension tan_neb(3*mxatom),xwork(3*mxatom),xleft(3*mxatom),
     +    xright(3*mxatom)      ! left,right for roll-checks

      dimension barrevev(nstate-1),xstates(nstate*natom*3)

c      parameter (movmax=10000)
c      include 'parameters.h'
      parameter (movmaxvec=365)
      parameter (m3=movmax*3)
      parameter (m3v=movmaxvec*3)
      include 'wrksubs_tad2.h'
      common/wrkcm13/iwrk13,jwrk13,work13(lwrk13)
      common/wrkcm2/iwrk2,jwrk2, work2(m3),work3(m3),work4(m3),
     x              fill(lwrk2-3*m3)

      save isaddlecount

      common/bpuneb/nebbailcount

      nebbailcount=0
      lnoise=.false.
      icnt = 0
      gcrit = gcritin
      dvcrit = dvcritin
      drcrit = drcritin

c needed quantities
      if((ipr.ge.4).and.(rank_f.eq.0)) then
        write(6,*) 'SF: natom: ',natom,nmove,gcrit,dvcrit,drcrit,ipr
     +    ,nimage,eshallow,itan,iclimb,springk,lstore_neb
        write(6,*) 'SF: filnam: ',filnam1,filnam2
        write(6,*) 'SF: hotmix: ',hotmix
c        write(6,*) 'SF: Debug: ifundc=',ifundc
      end if
      inebcount=0
      ibadmodecount=0
      imodeflag=0
      drroll=transcrit/10.0d0   ! initial displacement from saddle for roll check
      if (drroll.lt.0.05d0) drroll=0.05d0
      ierr=0
      ididintermediate=0
      iresolve=1 ! to increase resolution of band if endpoint is highest

      idimer=idimerneb
      if((ipr.ge.6).and.(rank_f.eq.0)) write(6,*)'SF: inebirc=',inebirc
      if((ipr.ge.6).and.(rank_f.eq.0)) write(6,*)'SF: ifunneb=',ifunneb

      if((ipr.ge.5.and.idimer.eq.1).and.(rank_f.eq.0)) write(6,*)
     $     'SF: using dimer to refine saddle'

c----------------------------------------------------------------------
 102  call neb_tad(natom,nimage,x1,x2,xhot,hotmix,xsneb,isad,es,tan_neb
     +    ,xleft,xright,nmove,ityp,taxes,ietype,maxtyp,rcut,hyperratio
     +    ,gfac,irotalign,itermax,ierr,intermediate,transcrit,itranscrit
     +    ,gcrit,dvcrit,drcrit,ipr,eshallow,xstates,nstate,barrevev
     +    ,ereverse,itan,iclimb,springk,lstore_neb,filnam1,filnam2
     +    ,idimer,ifundc,ifunneb,imodeflag,lnoise,ntype,potnam,amu)
      es = es + e0
      iresolvehold=iresolve ! value of iresolve actually used in NEB
c      write(6,*) 'SF: Debug, after neb_tad 1: ifundc=',ifundc
      inebcount=inebcount+1

      if (ierr.eq.2) then
c no saddle found, try a denser neb to find saddle
c xxx at some point, may add force interpolation
c         ididintermediate=1
         if((ipr.ge.0).and.(rank_f.eq.0)) write(6,*)
     +       ' saddle_find WARNING! NEB returned without saddle'
         iresolve=2
c         ierr=0 ! xxxxxxxx skip denser neb for now
c do a denser neb in next block to see if we can resolve a barrier
      endif

c xxxx this seems never to be met (ierr=5)      
      if (ierr.eq.5.and.intermediate.eq.0) then
         if((ipr.ge.0).and.(rank_f.eq.0)) write(6,*)
     +       ' SF: WARNING! Convergence not reached!'
      endif
      
c----------------------------------------------------------------------
 100  continue
      if (intermediate.eq.1.or.ierr.eq.2) then
         if((ipr.ge.4.and.ierr.ne.2).and.(rank_f.eq.0)) write(6,*)
     +       ' SF: Intermediate minimum found by NEB'
         if((ipr.ge.4.and.ierr.eq.2).and.(rank_f.eq.0)) write(6,*)
     +       ' SF: NEB returned without saddle, doing denser NEB'
         if (ierr.ne.2) ididintermediate=1
         if (intermediate.eq.1) then
            ifundcuse=ifundc
            if (ifundcuse.eq.5) ifundcuse=0
            imode=0
            ereverseuse=0.0d0 ! don't do reverse barrier stuff here
            call  descent_check(natom,nmove,1,x2,x1,
     +          ityp,taxes,ietype,maxtyp,rcut,hyperratio,barrevev
     +          ,ereverseuse,isame1,irotalign,itermax,gfac,transcrit
     +          ,itranscrit,gcrit,dvcrit,drcrit,ipr,imode,ifundcuse)
c            write(6,*) 'SF: Debug, after dc 1, ifundc=',ifundc
            if (isame1.eq.1) then
c uuuu what to do now?               
c should be taken care of in iexperimental block below (second
C iexperimental block), way down in neb_tad
            endif
         endif
         if((ipr.ge.4.and.ierr.ne.2).and.(rank_f.eq.0)) write(6,*)
     +       ' SF: Intermediate minimum put in x2'
c second NEB, in case no saddle, or intermediate minimum, etc.
         call neb_tad(natom,nimage*iresolve,x1,x2,xhot,0.0d0,xsneb,isad
     +       ,es,tan_neb,xleft,xright,nmove,ityp,taxes,ietype,maxtyp
     +       ,rcut,hyperratio,gfac,irotalign,itermax,ierr,intermediate
     +       ,transcrit,itranscrit,gcrit,dvcrit,drcrit,ipr,eshallow
     +       ,xstates,nstate,barrevev,ereverse,itan,iclimb,springk
     +       ,lstore_neb,filnam1,filnam2,idimer,ifundc,ifunneb
     +       ,imodeflag,lnoise,ntype,potnam,amu)
         es = es + e0
         iresolvehold=iresolve ! value of iresolve actually used in NEB
c         write(6,*) 'SF: Debug, after neb_tad 2, ifundc=',ifundc
         inebcount=inebcount+1
         iexperimental=1
         if (inebcount.gt.3.and.intermediate.eq.1.and.
     +       iexperimental.eq.1) then
            if (iresolve.ne.2) then
               iresolve=2
            else if (inebcount.gt.10) then
              if (rank_f.eq.0)then
               write(6,*) 'SF: WARNING! NEB seems stuck in loop!'
               write(6,*) 'SF: WARNING!   did NEB ',inebcount,' times,'
               write(6,*) 'SF: WARNING!   even did denser NEB,'
               write(6,*) 'SF: WARNING!   taking max energy as saddle,'
               write(6,*) 'SF: WARNING!   hoping for best.'
              endif
              return
            endif
         endif
         if (intermediate.eq.1) goto 100
         if (ierr.eq.1.and.intermediate.eq.0.and.ipr.ge.0) then
            if(rank_f.eq.0) write(6,*)
     +       ' SF: WARNING! Convergence not reached!'
         endif
         if (ierr.eq.2) then
c did a denser neb and still could find no saddle...
c passing back saddle geometry and energy (which should be same as an
C end point) and hoping for the best
c            ididintermediate=1 ! didn't do intermediate, just no saddle
            if ((ipr.ge.0).and.(rank_f.eq.0)) write(6,*)
     +       ' SF: denser NEB still returned with no saddle!'
            return
         endif
      endif

c Calculate saddle modes -- need to try again if nneg not 1:
      if((.false.).and.(ierr.ne.2).and.(ivineyard.gt.0))then
 117  continue
        if(nmove.eq.natom) nmove_tmp = nmove-1
        nnegsad = 0
        nnegmin = 0
        call dovineyard(natom,nmove_tmp,xsneb,ityp,taxes,ietype
     +         ,maxtyp,rcut,amu,e,grad,nnegsad,nprod,freqlogssad
     +         ,1,freqlogsmin,nnegmin,prefacdum,ierrv,ipr)
        if(nnegsad.ne.1)then
          lnoise=.true.
          icnt = icnt + 1
          if(icnt.gt.25)then
            if(rank_f.eq.0)then
              write(*,*) 'spawnID ',spawnID,' cannot find good saddle!!'
            endif
            STOP 'THIS SADDLE IS NO GOOD !!!'
          endif
          !dvcrit - energy (h) conv. crit; 5.d-8 is good for EAM
          !gcrit  - total grad length (h/A) conv. crit; 1.d-4 is good for EAM
          !drcrit - single-atom displacement (A) conv. crit; 1.d-5 is good for EAM
          gcrit = gcritin*(0.0001)
          dvcrit = dvcritin*(0.001)
          drcrit = drcritin*(0.001)
          if(rank_f.eq.0) write(*,*) 'Found saddle with nnegsad '
     +         ,nnegsad,' -- adding noise -- icnt ',icnt

          call neb_tad(natom,nimage*iresolve,x1,x2,xhot,0.0d0,xsneb,isad
     +       ,es,tan_neb,xleft,xright,nmove,ityp,taxes,ietype,maxtyp
     +       ,rcut,hyperratio,gfac,irotalign,itermax,ierr,intermediate
     +       ,transcrit,itranscrit,gcrit,dvcrit,drcrit,ipr,eshallow
     +       ,xstates,nstate,barrevev,ereverse,itan,iclimb,springk
     +       ,lstore_neb,filnam1,filnam2,idimer,ifundc,ifunneb
     +       ,imodeflag,lnoise,ntype,potnam,amu)
          es = es + e0

          goto 117
        endif
      endif
      lnoise=.false.

c----------------------------------------------------------------------
c roll check from saddle to see if connects minima as thought

      inebrefined=0
 300  continue
      if(rank_f.eq.0)then
       if(ipr.ge.10) write(6,*) 'SF: Debug 1: ierr,ifundc: ',ierr,ifundc
      endif
      if (ierr.eq.3.or.ifundc.eq.5) goto 400 ! shallow minimum (or thermal quench), rollcheck will fail so skip
c 300  continue !xxx maybe in wrong place, denser neb went here so skipped above check
c      write(6,*) 'SF: Debug 2: ierr,ifundc: ',ierr,ifundc
c check xa
c      if (ipr.ge.6) write(6,*) 'isad=',isad,': nimage*iresolve=',nimage
c     +    *iresolvehold
      if (itranscrit.eq.0.or.isad.eq.nimage*iresolvehold) then
         do j=1,3*natom
            xwork(j)=xsneb(j)-drroll*tan_neb(j) ! xxxx if descent_check only takes 1 iter, increase drroll
         enddo
      else
         call vecmov(xleft,xwork,3*natom)
         if((ipr.ge.6).and.(rank_f.eq.0))
     +       write(6,*) 'using xleft for rollcheck'
      endif
      imode=0
      ifundchold=ifundc
c      ifundc=0 cwp use default minimizer 
      ereverseuse=0.0d0 ! don't do reverse barrier stuff
      call descent_check(natom,nmove,1,xwork,x1,
     +    ityp,taxes,ietype,maxtyp,rcut,hyperratio,barrevev,ereverseuse
     +    ,isame1,irotalign,itermax,gfac,transcrit,itranscrit,gcrit
     +    ,dvcrit,drcrit,ipr,imode,ifundc)
      if (isame1.eq.1) then
         if ((ipr.ge.5).and.(rank_f.eq.0)) write(6,*)
     +   'saddle_find: note: Saddle connects to x1.'
      else if(isame1.eq.0) then
         if ((ipr.ge.4).and.(rank_f.eq.0)) write(6,*)
     +       'SF: WARNING! Saddle does not connect to x1!'
         filnam='rcfail.x1.dat'
         if ((ipr.ge.10).and.(rank_f.eq.0))
     +       call storefile(natom,xwork,xwork,ityp,taxes,filnam)
      else 
         if ((ipr.ge.0).and.(rank_f.eq.0)) write(6,*) 'isame1 =',isame1
         stop 'saddle_find: isame1 value unrecognized'
      endif

c check xb
      if (itranscrit.eq.0.or.isad.eq.nimage*iresolvehold) then
         do j=1,3*natom
            xwork(j)=xsneb(j)+drroll*tan_neb(j)
         enddo
      else
         call vecmov(xright,xwork,3*natom)
         if ((ipr.ge.6).and.(rank_f.eq.0))
     +        write(6,*) 'using xright for rollcheck'
      endif
      imode=0
      ereverseuse=0.0d0 ! don't do reverse barrier stuff
      call descent_check(natom,nmove,1,xwork,x2,
     +    ityp,taxes,ietype,maxtyp,rcut,hyperratio,barrevev,ereverseuse
     +    ,isame2,irotalign,itermax,gfac,transcrit,itranscrit,gcrit
     +    ,dvcrit,drcrit,ipr,imode,ifundc)
      if (isame2.eq.1) then
         if ((ipr.ge.5).and.(rank_f.eq.0)) write(6,*)
     +   'saddle_find: note: Saddle connects to x2.'
      else if(isame2.eq.0) then
         if ((ipr.ge.0).and.(rank_f.eq.0)) write(6,*)
     +       ' SF: WARNING! Saddle does not connect to x2!'
         filnam='rcfail.x2.dat'
         if ((ipr.ge.10).and.(rank_f.eq.0))
     +       call storefile(natom,xwork,xwork,ityp,taxes,filnam)
      else 
         if ((ipr.ge.0).and.(rank_f.eq.0)) write(6,*) 'isame2 =',isame2
         stop 'saddle_find: isame2 value unrecognized'
      endif
      ifundc=ifundchold

      if (idimer.eq.1.and.isame1.eq.1.and.isame2.eq.1) then
        if(rank_f.eq.0)then
         if (ipr.ge.5) write(6,*) 'SF: dimer connects both x1, x2'
        endif
        goto 400
      elseif (idimer.eq.1.and.isame1.eq.1.and.isame2.eq.0) then
        if(rank_f.eq.0)then
         if (ipr.ge.4) write(6,*) 'SF: dimer does not connect to x2'
         if (ipr.ge.4) write(6,*) 'SF: doing new NEB with new x2'
        endif
        do i=1,3*natom
            x2(i)=xwork(i)
        enddo
        intermediate=1
        goto 100
      elseif (idimer.eq.1.and.isame1.eq.0) then
        if(rank_f.eq.0)then
         if (ipr.ge.4) write(6,*) 'SF: dimer does not connect to x1'
        endif
        idimer=0
        goto 102
      endif

      if ((isame1.eq.0.or.isame2.eq.0).and.inebrefined.eq.0.and.
     +    inebirc.eq.0) then
         if((ipr.ge.4).and.(rank_f.eq.0))
     +       write(6,*) 'SF: doing denser NEB'
         call neb_tad(natom,nimage*2,x1,x2,xhot,hotmix,xsneb,isad,es
     +       ,tan_neb,xleft,xright,nmove,ityp,taxes,ietype,maxtyp,rcut
     +       ,hyperratio,gfac,irotalign,itermax,ierr,intermediate
     +       ,transcrit,itranscrit,gcrit,dvcrit,drcrit,ipr,eshallow
     +       ,xstates,nstate,barrevev,ereverse,itan,iclimb,springk
     +       ,lstore_neb,filnam1,filnam2,idimer,ifundc,ifunneb
     +       ,imodeflag,lnoise,ntype,potnam,amu)
         es = es + e0
         iresolvehold=2
c         write(6,*) 'SF: Debug, after neb_tad 3, ifundc=',ifundc
         inebcount=inebcount+1
         inebrefined=1
         goto 300
      else if ((isame1.eq.0.or.isame2.eq.0).and.inebrefined.eq.1.and
     +       .intermediate.eq.0) then
         ierr=1
         if ((ipr.ge.0).and.(rank_f.eq.0)) write(6,*)
     +       ' SF: WARNING! refined NEB finds no intermediate min,'
     +       ,' but roll check does; ignoring roll check.'
      else if (inebrefined.eq.1.and.intermediate.eq.1) then
         if ((ipr.ge.4).and.(rank_f.eq.0)) write(6,*)
     +       'SF: refined NEB found intermediate min,'
     +       ,' redoing NEB with new min.'
         goto 100
      else if (inebirc.eq.1) then
         if ((ipr.ge.4).and.(rank_f.eq.0)) 
     +       write(6,*) 'SF: inebirc on, ignoring roll check'
      endif

      if(rank_f.eq.0) write(6,*)'SF: about to do modecheck: ',imodecheck

      if (imodecheck.eq.1) then
         call wrkclaimv(13,1,'neb') ! for hessian curv array (work13)
         
         ip=-1
         ivec=0
         gdelt=1.d-4
         call nmodes(ip,ivec,natom,nmove,xsneb,ityp,taxes,ietype,maxtyp
     +       ,rcut,gdelt,amu,e,xwork,work13,nneg,nprod,freqlog)
         call wrkrelease(13)
         if(rank_f.eq.0) write(6,*)
     +      'SF: first two negative modes of saddle: ',work2(1),work2(2)
         if (work2(2).lt.rmodemag) then
            
            isaddlecount=isaddlecount+1
            
            filnamsad='badsaddle.'//chri(isaddlecount)//'.dat'
            if(rank_f.eq.0)then
              write(6,*) 'SF: writing badsaddle file: ',filnamsad
            endif
            call storefile(natom,xsneb,xsneb,ityp,taxes,filnamsad)

            ibadmodecount=ibadmodecount+1
            if (ibadmodecount.le.5) then
               if(rank_f.eq.0)then
                 write(6,*) 'SF: saddle has too many negative modes!'
                 write(6,*) 'SF: second mode: ',work2(2),rmodemag
                 write(6,*) 'SF: ibadmodecount: ',ibadmodecount
               endif
               dvcrituse=dvcrit/(10.0d0*real(ibadmodecount))
               drcrituse=drcrit/(10.0d0*real(ibadmodecount))
               gcrituse=gcrit/(10.0d0*real(ibadmodecount))
c            do i=2,nimage-1
c               do j=1,nmove*3
c                  (j+(i-1)*3*natom)=(j+(i-1)*3*natom)+rand
c               enddo
c            enddo
               imodeflag=1
               call neb_tad(natom,nimage,x1,x2,xhot,hotmix,xsneb,isad,es
     +             ,tan_neb,xleft,xright,nmove,ityp,taxes,ietype,maxtyp
     +             ,rcut,hyperratio,gfac,irotalign,itermax,ierr
     +             ,intermediate,transcrit,itranscrit,gcrituse,dvcrituse
     +             ,drcrituse,ipr,eshallow,xstates,nstate,barrevev
     +             ,ereverse,itan,iclimb,springk,lstore_neb,filnam1
     +             ,filnam2,idimer,ifundc,ifunneb,imodeflag,lnoise
     +             ,ntype,potnam,amu)
               es = es + e0
               imodeflag=0
               goto 100
            else 
               if(rank_f.eq.0)then
                 write(6,*) 'SF: saddle has too many negative modes!'
                 write(6,*) 'SF: second mode: ',work2(2),rmodemag
                 write(6,*) 'SF: ibadmodecount: ',ibadmodecount
                 write(6,*) 'SF: reached max ibadmodecount, continuing'
               endif
            endif
         endif
      endif

 400  continue

      return
      end
      
c**********************************************************************
      subroutine neb_tad(natom,nimage,x1,x2,xhot,hotmix,xs,imax,es
     +    ,tan_neb,xleft,xright,nmove,ityp,taxes,ietype,maxtyp,rcut
     +    ,hyperratio,gfac,irotalign,itermax,ierr,intermediate,transcrit
     +    ,itranscrit,gcrit,dvcrit,drcrit,ipr,eshallow,xstates,nstate
     +    ,barrevev,ereverse,itan,iclimb,springk,lstore_neb,filnam1
     +    ,filnam2,idimer,ifundc,ifunneb,imodeflag,lnoise
     +    ,ntype,potnam,amu)

c this is sort of a wrapper to NEB for TAD.  it expects only the two
c minima coordinates (plus some potential parameters and cell size),
c and returns the saddle point coordinates and energy.
c 
c input variables: natom: number of atoms
c                  xa,ya,za: coordinates of minimum 1
c                  xb,yb,zb: coordinates of minimum 2
c                  nmove: number of moving atoms
c                  ityp: array of atom types
c                  tax,tay,taz: cell size
c                  ietype: potential type
c                  maxtyp: maximum number of types (rcut size)
c                  rcut: matrix of potential cutoffs
c                  hyperratio: legacy from parrep code
c                  gfac: steepest decent step size
c                  itermax: maximum number of iterations to perform
c                  lpr: print lots of info?
c
c output variables: xs,ys,zs: saddle point coordinates
c                   es: saddle point energy (Hartree)
c                   ierr: error code
c                   xb,yb,zb are changed to intermediate minimum if there is one
c                   intermediate: 1 if an intermediate minimum found, zero otherwise

c ierr codes: ierr=0: no error
c             ierr=1: convergence not reached
c             ierr=2: NEB finds no saddle
      
      use mod_lammps
      implicit real*8(a-h,o-z)
c      parameter(mxatom=10000,mximage=42,
      include 'parameters.h'
      parameter(length=3*mximage*mxatom)
      dimension ityp(natom)

      dimension x1(3*natom),x2(3*natom),xs(3*natom),tan_neb(3*natom)
      dimension xleft(3*natom),xright(3*natom)
      dimension xhot(3*natom)
      dimension rcut(maxtyp,maxtyp)
      dimension q(length)

      dimension x(length),f(length),e(mximage),finfo(mximage*4)
      dimension v(length) ! for quickmin
      logical lstore_neb,ltemp,lnoise
      character*180 filnam1,filnam2
      dimension taxes(3)
      dimension x12(6*mxatom)

      !! LOCAL STUFF !!!
      logical :: use_local_neb = .false.
      real*8  :: freebuf   = 12.0 !!
      real*8  :: fixbuf    = 12.0 !!
      dimension x1_g(3*natom),x2_g(3*natom)
      dimension x1_save(3*natom),x2_save(3*natom)
      dimension xs_save(3*natom),tan_neb_save(3*natom)
      dimension x1lmp(3*mxatom)
      integer   l_to_g_map(3*natom)
      dimension ityp_g(natom)
      character*80 filn,filnb
      character*80 potnam(maxtyp)
      dimension p1(3*natom)
      dimension amu(maxtyp)
      !!!!!!!!!!!!!!!!!!

      dimension barrevev(nstate-1),xstates(nstate*natom*3)
      dimension nstep(mximage),deltat(mximage),alpha(mximage) ! cwp for FIRE
      common/bpuneb/nebbailcount

      if (nimage.gt.mximage) stop 'Increase mximage!'

c      write(6,*) 'NEB_TAD: Debug, ifundc=',ifundc
      
c minimization routines
c ifunneb defines which minimization routine is used:
c     ifunneb=0: steepest descent
c     ifunneb=1: quickmin
c      ifunneb=2: FIRE cwp
c     ifunneb=3: LBFGS-B minimization 

      if (ifunneb.lt.0.or.ifunneb.ge.4) stop 'ifunneb not defined'

      ev=27.21d0
      intermediate=0
      ierr=0
      imodedone=0

      !! LOCAL STUFF !!!
      ! Don't use local nebs for ionic systems yet...
      ! (might need to add/subtract from this list)
      if((ietype.eq.72).or.(ietype.eq.103).or.(ietype.eq.172))then
          use_local_neb = .false.
      endif
      !! Cant do Local NEBS with LAMMPS yet... :/
      !if((ietype.ge.100).and.(ietype.le.200))then
      !    use_local_neb = .false.
      !endif
      p1(:) = 0.d0
      if(use_local_neb)then

          ! Determin moving atoms:
          x1_save(:) = x1(:)
          x2_save(:) = x2(:)
          xs_save(:) = x1(:)
          tan_neb_save(:) = 0.d0
          x1_g(:) = x1(:)
          x2_g(:) = x2(:)
          ityp_g(:) = ityp(:)
          nmove_g = nmove
          natom_g = natom
          natom = 0
          nmove = 0
          x1(:) = 0.d0
          x2(:) = 0.d0
          ityp(:) = 0
          l_to_g_map(:) = 0
          xmax = -1.d10
          xmin = 1.d10
          ymax = -1.d10
          ymin = 1.d10
          zmax = -1.d10
          zmin = 1.d10
          xref = 0.d0
          yref = 0.d0
          zref = 0.d0
          ninv = 0
          do i=1, nmove_g

             ! Get change in atom coords:
             ddx = x2_g((i-1)*3+1) - x1_g((i-1)*3+1)
             do while (ddx.gt.(0.5*taxes(1)))
                 x2_g((i-1)*3+1) = x2_g((i-1)*3+1) - taxes(1)
                 ddx = x2_g((i-1)*3+1) - x1_g((i-1)*3+1)
             enddo
             do while (ddx.lt.(-0.5*taxes(1)))
                 x2_g((i-1)*3+1) = x2_g((i-1)*3+1) + taxes(1)
                 ddx = x2_g((i-1)*3+1) - x1_g((i-1)*3+1)
             enddo
             ddy = x2_g((i-1)*3+2) - x1_g((i-1)*3+2)
             do while (ddy.gt.(0.5*taxes(2)))
                 x2_g((i-1)*3+2) = x2_g((i-1)*3+2) - taxes(2)
                 ddy = x2_g((i-1)*3+2) - x1_g((i-1)*3+2)
             enddo
             do while (ddy.lt.(-0.5*taxes(2)))
                 x2_g((i-1)*3+2) = x2_g((i-1)*3+2) + taxes(2)
                 ddy = x2_g((i-1)*3+2) - x1_g((i-1)*3+2)
             enddo
             ddz = x2_g((i-1)*3+3) - x1_g((i-1)*3+3)
             do while (ddz.gt.(0.5*taxes(3)))
                 x2_g((i-1)*3+3) = x2_g((i-1)*3+3) - taxes(3)
                 ddz = x2_g((i-1)*3+3) - x1_g((i-1)*3+3)
             enddo
             do while (ddz.lt.(-0.5*taxes(3)))
                 x2_g((i-1)*3+3) = x2_g((i-1)*3+3) + taxes(3)
                 ddz = x2_g((i-1)*3+3) - x1_g((i-1)*3+3)
             enddo

             if (sqrt(ddx**2+ddy**2+ddz**2).gt.transcrit) then
                 ninv = ninv + 1

!                 if(rank_f.eq.0)
!     +             write(*,*) 'SpawnID ',spawnid,' adding atom in NEB!'

                 if (ninv.eq.1) then
                     xref = x1_g((i-1)*3+1)
                     yref = x1_g((i-1)*3+2)
                     zref = x1_g((i-1)*3+3)
                 else
                     ! Move coord to be close as possible to the referece atom
                     do while((x1_g((i-1)*3+1)-xref).gt.(0.5*taxes(1)))
                         x1_g((i-1)*3+1) = x1_g((i-1)*3+1) - taxes(1)
                     enddo
                     do while((x1_g((i-1)*3+1)-xref).lt.(-0.5*taxes(1)))
                         x1_g((i-1)*3+1) = x1_g((i-1)*3+1) + taxes(1)
                     enddo
                     do while((x1_g((i-1)*3+2)-yref).gt.(0.5*taxes(2)))
                         x1_g((i-1)*3+2) = x1_g((i-1)*3+2) - taxes(2)
                     enddo
                     do while((x1_g((i-1)*3+2)-yref).lt.(-0.5*taxes(2)))
                         x1_g((i-1)*3+2) = x1_g((i-1)*3+2) + taxes(2)
                     enddo
                     do while((x1_g((i-1)*3+3)-zref).gt.(0.5*taxes(3)))
                         x1_g((i-1)*3+3) = x1_g((i-1)*3+3) - taxes(3)
                     enddo
                     do while((x1_g((i-1)*3+3)-zref).lt.(-0.5*taxes(3)))
                         x1_g((i-1)*3+3) = x1_g((i-1)*3+3) + taxes(3)
                     enddo
                 endif

                 if (x1_g((i-1)*3+1).gt.xmax) then
                     xmax = x1_g((i-1)*3+1)
                 endif
                 if (x1_g((i-1)*3+1).lt.xmin) then
                     xmin = x1_g((i-1)*3+1)
                 endif

                 if (x1_g((i-1)*3+2).gt.ymax) then
                     ymax = x1_g((i-1)*3+2)
                 endif
                 if (x1_g((i-1)*3+2).lt.ymin) then
                     ymin = x1_g((i-1)*3+2)
                 endif

                 if (x1_g((i-1)*3+3).gt.zmax) then
                     zmax = x1_g((i-1)*3+3)
                 endif
                 if (x1_g((i-1)*3+3).lt.zmin) then
                     zmin = x1_g((i-1)*3+3)
                 endif

             endif
          enddo
          if (ninv <= 0) then
              !stop 'ninv <= 0 !!!'
              if(rank_f.eq.0) 
     +          write(*,*) 'SpawnID ',spawnid,' NOT DOING LOCAL NEB!'
     +            ,' ninv <= 0 !!!'
              use_local_neb = .false.
              x1(:) = x1_save(:)
              x2(:) = x2_save(:)
              ityp(:) = ityp_g(:)
              natom = natom_g
              nmove = nmove_g
              goto 217
          endif

          !! Move atoms into "fixed" cell boundary:
          xminf = xmin-freebuf-fixbuf
          yminf = ymin-freebuf-fixbuf
          zminf = zmin-freebuf-fixbuf
          xmaxf = xmax+freebuf+fixbuf
          ymaxf = ymax+freebuf+fixbuf
          zmaxf = zmax+freebuf+fixbuf
          xminfr = xmin-freebuf
          yminfr = ymin-freebuf
          zminfr = zmin-freebuf
          xmaxfr = xmax+freebuf
          ymaxfr = ymax+freebuf
          zmaxfr = zmax+freebuf

          xallfree = 0
          yallfree = 0
          zallfree = 0
          if((xmaxf-xminf).gt.taxes(1)) then
              xallfree = 1
              xminf = 0.d0
              xmaxf = taxes(1)
              xminfr = 0.d0 - 1.d0
              xmaxfr = taxes(1) + 1.d0
          endif
          if((ymaxf-yminf).gt.taxes(2)) then
              yallfree = 1
              yminf = 0.d0
              ymaxf = taxes(2)
              yminfr = 0.d0 - 1.d0
              ymaxfr = taxes(2) + 1.d0
          endif
          if((zmaxf-zminf).gt.taxes(3)) then
              zallfree = 1
              zminf = 0.d0
              zmaxf = taxes(3)
              zminfr = 0.d0 - 1.d0
              zmaxfr = taxes(3) + 1.d0
          endif


          if ((xallfree+yallfree+zallfree).lt.3) then

              !write(*,*)'SpawnID ',spawnid,' still using atom in NEB!'
              !write(*,*)'SpawnID ',spawnid,' xmax,xmin = ',xmax,xmin
              !write(*,*)'SpawnID ',spawnid,' ymax,ymin = ',ymax,ymin
              !write(*,*)'SpawnID ',spawnid,' zmax,zmin = ',zmax,zmin
              !write(*,*)'SpawnID ',spawnid,' xmaxf,xminf = ',xmaxf,xminf
              !write(*,*)'SpawnID ',spawnid,' ymaxf,yminf = ',ymaxf,yminf
              !write(*,*)'SpawnID ',spawnid,' zmaxf,zminf = ',zmaxf,zminf

              do i = 1,natom_g
                 do while (x1_g((i-1)*3+1).gt.xmaxf)
                     x1_g((i-1)*3+1) = x1_g((i-1)*3+1) - taxes(1)
                 enddo
                 do while (x1_g((i-1)*3+1).lt.xminf)
                     x1_g((i-1)*3+1) = x1_g((i-1)*3+1) + taxes(1)
                 enddo
                 do while (x1_g((i-1)*3+2).gt.ymaxf)
                     x1_g((i-1)*3+2) = x1_g((i-1)*3+2) - taxes(2)
                 enddo
                 do while (x1_g((i-1)*3+2).lt.yminf)
                     x1_g((i-1)*3+2) = x1_g((i-1)*3+2) + taxes(2)
                 enddo
                 do while (x1_g((i-1)*3+3).gt.zmaxf)
                     x1_g((i-1)*3+3) = x1_g((i-1)*3+3) - taxes(3)
                 enddo
                 do while (x1_g((i-1)*3+3).lt.zminf)
                     x1_g((i-1)*3+3) = x1_g((i-1)*3+3) + taxes(3)
                 enddo
              enddo

              !write(*,*) 'SpawnID ',spawnid,' put NEB atoms in pbcell'

              !! Define Local Arrays:
              do i = 1, nmove_g
                  xx = x1_g((i-1)*3+1)
                  yy = x1_g((i-1)*3+2)
                  zz = x1_g((i-1)*3+3)
                  ! If Inside Free Box:
                  if(
     +               (xx.lt.(xmaxfr)).and.(xx.gt.(xminfr)).and.
     +               (yy.lt.(ymaxfr)).and.(yy.gt.(yminfr)).and.
     +               (zz.lt.(zmaxfr)).and.(zz.gt.(zminfr))
     +              )then
                      ! Add LOCAL Free Atom
                      natom = natom + 1
                      nmove = nmove + 1
                      x1((natom-1)*3+1) = x1_save((i-1)*3+1)
                      x1((natom-1)*3+2) = x1_save((i-1)*3+2)
                      x1((natom-1)*3+3) = x1_save((i-1)*3+3)
                      x2((natom-1)*3+1) = x2_save((i-1)*3+1)
                      x2((natom-1)*3+2) = x2_save((i-1)*3+2)
                      x2((natom-1)*3+3) = x2_save((i-1)*3+3)
                      ityp(natom) = ityp_g(i)
                      l_to_g_map(natom) = i
                  endif
              enddo
              do i = 1, natom_g
                  xx = x1_g((i-1)*3+1)
                  yy = x1_g((i-1)*3+2)
                  zz = x1_g((i-1)*3+3)
                  ! If Inside Fixed Box:
                  if(
     +               (xx.lt.xmaxf).and.(xx.gt.xminf).and.
     +               (yy.lt.ymaxf).and.(yy.gt.yminf).and.
     +               (zz.lt.zmaxf).and.(zz.gt.zminf)
     +              )then
                      ! But Outside Free Box
                      if(
     +                   (xx.ge.xmaxfr).or.(xx.le.xminfr).or.
     +                   (yy.ge.ymaxfr).or.(yy.le.yminfr).or.
     +                   (zz.ge.zmaxfr).or.(zz.le.zminfr).or.
     +                   (i.gt.nmove_g) ! OR a Global Fixed Atom
     +                  )then
                          ! Add LOCAL Fixed Atom
                          natom = natom + 1
                          x1((natom-1)*3+1) = x1_save((i-1)*3+1)
                          x1((natom-1)*3+2) = x1_save((i-1)*3+2)
                          x1((natom-1)*3+3) = x1_save((i-1)*3+3)
                          x2((natom-1)*3+1) = x2_save((i-1)*3+1)
                          x2((natom-1)*3+2) = x2_save((i-1)*3+2)
                          x2((natom-1)*3+3) = x2_save((i-1)*3+3)
                          ityp(natom) = ityp_g(i)
                          l_to_g_map(natom) = i
                      endif
                  endif
              enddo
          else ! Not worth doing local stuff
              if(rank_f.eq.0) 
     +          write(*,*) 'SpawnID ',spawnid,' NOT DOING LOCAL NEB!'
              use_local_neb = .false.
              x1(:) = x1_save(:)
              x2(:) = x2_save(:)
              ityp(:) = ityp_g(:)
              natom = natom_g
              nmove = nmove_g
          endif
      endif

 217  continue

!       write(*,*) 'SpawnID ',spawnid
!     +   ,' doing local neb with natom,nmove = ',natom,nmove

       if(use_local_neb)then
           if(rank_f.eq.0) write(*,*) 'SpawnID ',spawnid
     +       ,' doing LOCAL neb with natom,nmove = ',natom,nmove

!           call flush(6)
!           filn=''
!           write(filn,*) spawnID
!           filnb=''
!           write(filnb,*) spawnID
!           filn='nebtest.'//trim(adjustl(filn))//
!     +               '.'//trim(adjustl(filnb))//'.dat'
!           call storefile(natom,x1,p1,ityp,taxes,filn)
!           call flush(6)
!           call MPI_BARRIER(force_comm,ier)

           ! Startup fresh lammps run
           if((ietype.ge.100).and.(ietype.le.200))then
             x1lmp(:) = 0.d0
             x1lmp(1:3*natom) = x1(1:3*natom)
             call MPI_BARRIER(force_comm,ier)
             call lmp_setup(natom,nmove,ntype,maxtyp,potnam,ietype
     +           ,x1lmp,taxes,ityp,rcut,amu,templow,temphigh
     +           ,.true.,dt)
             call MPI_BARRIER(force_comm,ier)
           endif

       else
           if(rank_f.eq.0) write(*,'(A,I9,A,I9,1X,I9)') 
     +        'SpawnID ',spawnid
     +       ,' doing GLOBAL neb with natom,nmove = ',natom,nmove
       endif
      !!!!!!!!!!!!!!!!!!


      if (3*natom*nimage.gt.length) then
        if(rank_f.eq.0)then
          if (ipr.ge.0) write(6,*) '#NEB: too many atoms*images!'
        endif
        stop 'aborting - natom*nimage space exceeded in neb_tad'
      endif

      if(rank_f.eq.0) write(6,*) 'neb_tad: debug: imodeflag: ',imodeflag

c cwp initialize FIRE ifunneb=2
      if (ifunneb.eq.2) then ! FIRE
         deltatmax=1.0
         rmaxjump=0.2
         finc=1.1
         fdec=0.5
         fadec=0.99
         nmin2=5
         alphaini=0.1
         do i=1,nimage
            deltat(i)=0.1
            alpha(i)=alphaini
            nstep(i)=0
         enddo
      endif

c----------------------------------------------------------------------
c prep chain
 100  continue
      !if(.not.lqeqfix) call refix_qeq()
      call neb_init_linear(nimage,natom,x1,x2,x,xhot,hotmix
     +                                                 ,q,nmove,lnoise)
      if((rank_f.eq.0).and.(ipr.ge.4))
     +   write(6,*)'#NEB: initialized chain'
      if (ietype.eq.0) then
         call neb_widen_chain(nimage,natom,x,f,e,nmove,ityp,taxes,ietype
     +       ,maxtyp,rcut,hyperratio,q)
c         write(6,*) '#NEB: widened chain'
      endif
      call neb_gcalc(nimage,natom,x,f,e,nmove,ityp,taxes,ietype
     +    ,maxtyp,rcut,hyperratio,1,dvmax,q)
      !if(lqeqfix) call unfix_qeq()
      call neb_nudge(nimage,natom,nmove,x,f,e,tan_neb,springk,itan
     $     ,iclimb,finfo,fmax,idimer)
      if ((ipr.ge.5).and.(rank_f.eq.0)) write(6,*) '#NEB: initial chain'
      if (ipr.ge.5) call neb_write_chain(nimage,natom,x,f,e,finfo
     +    ,lstore_neb,filnam1,taxes,ityp,springk,ipr)

      do i=1,3*nimage*natom
         v(i)=0.0d0
      enddo

c----------------------------------------------------------------------
c main loop
      i=0
      dvmax=1.0d0
      dxmax=1.0d0
      dxmaxold=1.0d10
      ndxmax=0
      ndxmaxlow=0
      fhighest=gcrit*1.1d0
      itermaxneb=itermax
      if (idimer.eq.1) itermaxneb=itermaxneb/20
      do while(i.lt.itermaxneb.and.(dvmax.gt.dvcrit.or.dxmax.gt.drcrit
     $     .or.fhighest.gt.gcrit))
         call neb_gcalc(nimage,natom,x,f,e,nmove,ityp,taxes,ietype
     +       ,maxtyp,rcut,hyperratio,0,dvmax,q)
         call neb_nudge(nimage,natom,nmove,x,f,e,tan_neb,springk,itan
     $        ,iclimb,finfo,fmax,idimer)
         
c adust gfac if force too big or step size not reducing
         gfacnow=gfac
         fbig=5.0d0
         if (fmax*ev.gt.fbig) gfacnow=gfacnow*fbig/fmax/ev
         if (abs(dxmax-dxmaxold).lt.drcrit.and.dxmax.gt.drcrit*10.0d0) 
     +       ndxmax=ndxmax+1
         if (ndxmax.gt.100) then
            gfacnow=gfacnow/10.0d0
            ndxmaxlow=ndxmaxlow+1
            if((ndxmaxlow.eq.1.and.ipr.ge.6).and.(rank_f.eq.0))
     +         write(6,*) 'NEB: turning down gfac=',gfacnow,dxmax
         endif
         if (ndxmaxlow.gt.100) then
            ndxmaxlow=0
            ndxmax=0
         endif
         
         dxmaxold=dxmax
         if (ifunneb.eq.0) then
            call neb_steepest(nimage,natom,nmove,x,f,gfacnow,dxmax)
         else if (ifunneb.eq.1) then
            call neb_quickmin(nimage,natom,nmove,x,f,v,gfacnow,dxmax)
         else if (ifunneb.eq.2) then ! cwp FIRE
            call neb_fire(nimage,natom,nmove,x,f,v,dxmax,rmaxjump
     +           ,alpha,deltat,deltatmax,finc,fadec,nstep
     +           ,nmin2,alphaini,fdec)
         else
            if (rank_f.eq.0) write(6,*) 'ifunneb not recognized'
            stop 'ifunneb not recognized'
         endif
         i=i+1
         ehighest=-1.0d10
         jmin=0
         jmax=0
         do j=2,nimage-1
            if (e(j).gt.ehighest) then
               fhighest=finfo((j-1)*4+1)
               ehighest=e(j)
               jmax=j
            endif
c xxxx bail on neb if intermediate minimum, now (050823) checking if
c belongs to basin or not
            if (e(j).lt.e(j-1).and.e(j).lt.e(j+1).and.jmin.eq.0) then
               imode=0
               niter=i
               call neb_inter_min(e,niter,j,jmin,nimage,eshallow,ierr
     +             ,xstates,nstate,barrevev,ereverse,natom,nmove,x(1
     +             +(j-1)*3*natom),ityp,taxes,ietype,maxtyp,rcut
     +             ,hyperratio,isame,irotalign,itermax,gfac,transcrit
     +             ,itranscrit,gcrit,dvcrit,drcrit,ipr,imode,ifundc)

               if((ipr.ge.6.and.imin.ne.0).and.(rank_f.eq.0)) then
                  write(6,*) '#NEB: image ',i,' is an intermediate min'
               endif

               iexperimental=1  ! this block copied from below because want to make sure intermediate minimum is not an end point
               if (iexperimental.eq.1.and.jmin.ne.0.and.j.eq.jmin.and
     +             .nebbailcount.lt.10) then
c verify that intermediate minimum really doesn't belong to either end
c point basin (for some transition criteria, some intermediate minima
c might be within transcrit of another basin)
                  do k=1,3*natom
                     x12(k)=x1(k)
                     x12(k+3*natom)=x2(k)
                     xs(k)=x(k+(jmin-1)*3*natom) ! use xs as temporary work array
                  enddo
                  imode=0
                  ereverseuse=0.0d0 ! don't do reverse barrier stuff
                  call descent_check(natom,nmove,2,xs,x12,
     +                ityp,taxes,ietype,maxtyp,rcut,hyperratio,barrevev
     +                ,ereverseuse,isame,irotalign,itermax,gfac
     +                ,transcrit,itranscrit,gcrit,dvcrit,drcrit,ipr
     +                ,imode,ifundc)
c                  write(6,*) 'NEB_TAD: Debug, after dc 1,ifundc=',ifundc
                  if (isame.eq.1.or.isame.eq.2) then
                     if((ipr.ge.4).and.(rank_f.eq.0)) write(6,*)
     +                   '#NEB: intermediate minimum belongs to basin '
     +                   ,isame,', neglecting'
                     jmin=0  
                  endif
               endif

               if (jmin.ne.0) then
                  intermediate=1
                  call vecmov(x(1+(jmin-1)*3*natom),x2,3*natom)
                  if (rank_f.eq.0) write(6,*)
     +                '#NEB: intermediate min during NEB,'
     +                ,' bailing after ',i,' iterations'
                  nebbailcount=nebbailcount+1
                  if (ipr.ge.4) call neb_write_chain(nimage,natom,x,f,e
     +                ,finfo,lstore_neb,filnam2,taxes,ityp,springk,ipr)
                  es=ehighest
                  do k=1,3*natom
                     xs(k)=x(k+(jmax-1)*3*natom)
                  enddo
                  if (use_local_neb) then
                     ityp(:) = ityp_g(:)
                     do i = 1, natom
                         x1_save((l_to_g_map(i)-1)*3+1) = x1((i-1)*3+1)
                         x1_save((l_to_g_map(i)-1)*3+2) = x1((i-1)*3+2)
                         x1_save((l_to_g_map(i)-1)*3+3) = x1((i-1)*3+3)
                         x2_save((l_to_g_map(i)-1)*3+1) = x2((i-1)*3+1)
                         x2_save((l_to_g_map(i)-1)*3+2) = x2((i-1)*3+2)
                         x2_save((l_to_g_map(i)-1)*3+3) = x2((i-1)*3+3)
                         xs_save((l_to_g_map(i)-1)*3+1) = xs((i-1)*3+1)
                         xs_save((l_to_g_map(i)-1)*3+2) = xs((i-1)*3+2)
                         xs_save((l_to_g_map(i)-1)*3+3) = xs((i-1)*3+3)
                         tan_neb_save((l_to_g_map(i)-1)*3+1) = 
     +                       tan_neb((i-1)*3+1)
                         tan_neb_save((l_to_g_map(i)-1)*3+2) = 
     +                       tan_neb((i-1)*3+2)
                         tan_neb_save((l_to_g_map(i)-1)*3+3) = 
     +                       tan_neb((i-1)*3+3)
                     enddo
                     x1(:) = x1_save(:)
                     x2(:) = x2_save(:)
                     xs(:) = xs_save(:)
                     natom = natom_g
                     nmove = nmove_g
                     ! Startup fresh lammps run
                     if((ietype.ge.100).and.(ietype.le.200))then
                       x1lmp(:) = 0.d0
                       x1lmp(1:3*natom) = x1(1:3*natom)
                       call MPI_BARRIER(force_comm,ier)
                       call lmp_setup(natom,nmove,ntype,maxtyp,potnam
     +                   ,ietype,x1lmp,taxes,ityp,rcut,amu
     +                   ,templow,temphigh,.true.,dt)
                       call MPI_BARRIER(force_comm,ier)
                     endif
                  endif
                  es = es - e(1)
                  return
               endif
            endif
         enddo
         if (iclimb.eq.1) then  ! if iclimb, converge only force of highest image
            dvmax=0.0d0
            dxmax=0.0d0
         endif
         ltemp=.false.
         if (mod(i,250).eq.0.and.ipr.ge.6) call neb_write_chain(nimage
     +       ,natom,x,f,e,finfo,ltemp,filnam1,taxes,ityp,springk,ipr)
         if ((mod(i,25).eq.0.and.ipr.ge.6).and.(rank_f.eq.0))
     +     write(6,6661) '#NEB: dvmax,dxmax,gmax: ',dvmax,dxmax,fhighest

         if (((dvmax.lt.dvcrit.and.dxmax.lt.drcrit.and.fhighest.lt.gcrit
     +       ).or.i.eq.int(itermaxneb/2)).and.imodeflag.eq.1.and
     +       .imodedone.eq.0) then
            dvmax=1.0d0
            dxmax=1.0d0
            fhighest=1.0d0
            imodedone=1
            if(rank_f.eq.0)
     +        write(6,*) 'neb_tad: adding noise to chain, step: ',i
            do j=2,nimage-1
               do k=1,nmove*3
                if(rank_f.eq.0)then
                  x((j-1)*3*natom+k)=x((j-1)*3*natom+k)+0.01d0*(0.5d0
     +                -prngen(0))
                endif
               enddo
               call MPI_BCAST
     +           (x((j-1)*3*natom+1),nmove*3,MPI_REAL8,0,force_comm,ier)
            enddo
         endif
      enddo ! cwp End of the main iteration loop

      if (i.ge.itermaxneb) then
         if(rank_f.eq.0)then
           if (ipr.ge.4) write(6,*) '#NEB: convergence not reached!'
           if (ipr.ge.6) write(6,6661) '#NEB: dvmax,dxmax,gmax: ',dvmax
     +       ,dxmax,fhighest
         endif
 6661    format(a30,3d20.10)
         ierr=1
      else
         if (rank_f.eq.0) then
           if (ipr.ge.4) write(6,*) '#NEB: convergence reached in ',i
     +       ,' iterations'
           if(ipr.ge.4)
     x        write(6,6663) dvmax,dvcrit,dxmax,drcrit,fhighest,gcrit
         endif
 6663    format('#NEB: converge crit: ',1p6d20.10)
      endif

      if ((nebbailcount.ge.10).and.(rank_f.eq.0)) write(6,*)
     +    '#NEB: WARNING: nebbailcount max reached'
      if ((ipr.ge.4.or.nebbailcount.ge.10).and.(rank_f.eq.0))
     +    write(6,*) '#NEB: final chain'
      if (ipr.ge.4.or.nebbailcout.ge.10) call neb_write_chain(nimage
     +    ,natom,x,f,e,finfo,lstore_neb,filnam2,taxes,ityp,springk,ipr)

      if (idimer.eq.1) then
         ni1max=0
         do i=2,nimage-1
            if (e(i).gt.e(i-1).and.e(i).gt.e(i+1).and.ni1max.eq.0) then ! find first maximum
               ni1max=i
            endif
         enddo
         if (ni1max.eq.0) then
            if ((ipr.ge.4).and.(rank_f.eq.0))
     $           write(6,*)'#NEB: Dimer: Max image is an end point'
            if (e(1).gt.e(nimage)) ni1max=2
            if (e(nimage).gt.e(1)) ni1max=nimage-1
         endif
         if ((ipr.ge.4).and.(rank_f.eq.0))
     +      write(6,*) '#NEB: first maximum is image ',ni1max
         do j=1,3*natom
            xs(j)=x(j+(ni1max-1)*3*natom)
            tan_neb(j)=x(j+(ni1max)*3*natom)-x(j+(ni1max-2)*3*natom)
         enddo
c needs version 2 of dimer_search
         if ((ipr.ge.6).and.(rank_f.eq.0))
     +      write(6,*) '#NEB: launching dimer'
         idimer=0
         ndimerdof=0
         iranx=0
         irandir=0
         imode=1
         rdimerdist=0.0d0
         ndimercoord=0
         idimeroverunder=0
         rdimerdisp=0.0d0
         idimerbail=0
         emin=10000d0
         eminbar=10000d0
         call dimer_search(natom,nmove,idimerbail,emin,eminbar,idimer
     +       ,ndimerdof,rdimerdist,ndimercoord,idimeroverunder
     +       ,rdimerdisp,x1,xs,x2,x2,es,tan_neb,ityp,taxes,ietype,maxtyp
     +       ,rcut,hyperratio,barrevev,ereverse,gfac,irotalign,itermax
     +       ,idimererr,iranx,irandir,transcrit,itranscrit,gcrit,dvcrit
     +       ,drcrit,ipr,imode,ifundc,eig)
         if(rank_f.eq.0)then
           if (ipr.ge.6) write(6,*) '#NEB: returned from dimer search'
           if (ipr.ge.4) write(6,*) '#NEB: energy of dimer: '
     +       ,es*27.21d0,' eV'
         endif
c no need to check for intermediate minima in the NEB as the NEB is not
c     too well converged, so return.  the roll check should catch these
c     things.
c     may need a new option that, if idimer is 1, then redo new NEB
c     not more densely but rather with the failed roll check point
c     (unless none connect to initial point)
         if (use_local_neb) then
             ityp(:) = ityp_g(:)
             do i = 1, natom
                 x1_save((l_to_g_map(i)-1)*3+1) = x1((i-1)*3+1)
                 x1_save((l_to_g_map(i)-1)*3+2) = x1((i-1)*3+2)
                 x1_save((l_to_g_map(i)-1)*3+3) = x1((i-1)*3+3)
                 x2_save((l_to_g_map(i)-1)*3+1) = x2((i-1)*3+1)
                 x2_save((l_to_g_map(i)-1)*3+2) = x2((i-1)*3+2)
                 x2_save((l_to_g_map(i)-1)*3+3) = x2((i-1)*3+3)
                 xs_save((l_to_g_map(i)-1)*3+1) = xs((i-1)*3+1)
                 xs_save((l_to_g_map(i)-1)*3+2) = xs((i-1)*3+2)
                 xs_save((l_to_g_map(i)-1)*3+3) = xs((i-1)*3+3)
                 tan_neb_save((l_to_g_map(i)-1)*3+1) =
     +               tan_neb((i-1)*3+1)
                 tan_neb_save((l_to_g_map(i)-1)*3+2) =
     +               tan_neb((i-1)*3+2)
                 tan_neb_save((l_to_g_map(i)-1)*3+3) =
     +               tan_neb((i-1)*3+3)
             enddo
             x1(:) = x1_save(:)
             x2(:) = x2_save(:)
             xs(:) = xs_save(:)
             natom = natom_g
             nmove = nmove_g
             ! Startup fresh lammps run
             if((ietype.ge.100).and.(ietype.le.200))then
               x1lmp(:) = 0.d0
               x1lmp(1:3*natom) = x1(1:3*natom)
               call MPI_BARRIER(force_comm,ier)
               call lmp_setup(natom,nmove,ntype,maxtyp,potnam
     +          ,ietype,x1lmp,taxes,ityp,rcut,amu
     +          ,templow,temphigh,.true.,dt)
               call MPI_BARRIER(force_comm,ier)
             endif
         endif
         es = es - e(1)
         return
      endif

c----------------------------------------------------------------------
c get saddle geometry, check for intermediate mins
      emax=-1.0d10
      imin=0
      imax=1
      do i=1,nimage
         if (e(i).gt.emax) then ! find saddle
            emax=e(i)
            imax=i
         endif
         if (i.gt.1.and.i.lt.nimage) then
            if (e(i).lt.e(i-1).and.e(i).lt.e(i+1).and.imin.eq.0) then ! find intermediate minimum
               imode=0
               niter=-1
               call neb_inter_min(e,niter,i,imin,nimage,eshallow,ierr
     +             ,xstates,nstate,barrevev,ereverse,natom,nmove,x(1
     +             +(i-1)*3*natom),ityp,taxes,ietype,maxtyp,rcut
     +             ,hyperratio,isame,irotalign,itermax,gfac,transcrit
     +             ,itranscrit,gcrit,dvcrit,drcrit,ipr,imode,ifundc)
c               imin=i
c               if (ipr.ge.4) then
c                  write(6,*) 'NEB note: Intermediate minimum found!'
c                  write(6,*) 'NEB note: image: ',imin,' with energy ',
c     $                e(i)*ev
c               endif
c               eleft=0.0d0
c               eright=0.0d0
c               do j=2,imin-1
c                  eleft=max(eleft,(e(j)-e(imin))*ev)
c               enddo
c               do j=imin+1,nimage-1
c                  eright=max(eright,(e(j)-e(imin))*ev)
c               enddo
c               if (eright.lt.eshallow.or.eleft.lt.eshallow) then
c                  imin=0
c                  if (ipr.ge.4) write(6,111) eleft,eright
c 111              format('#NEB: note: shallow minimum; barriers: '
c     +                ,2f12.5,' (eV); neglecting')
c                  ierr=3
c               endif
               if (ipr.ge.6.and.imin.ne.0) then
                  write(6,*) '#NEB: image ',i,' is an intermediate min'
               endif
               iexperimental=1
               if (iexperimental.eq.1.and.imin.ne.0.and.i.eq.imin) then
c verify that intermediate minimum really doesn't belong to either end
c point basin (for some transition criteria, some intermediate minima
c might be within transcrit of another basin)
                  do j=1,3*natom
                     x12(j)=x1(j)
                     x12(j+3*natom)=x2(j)
                     xs(j)=x(j+(imin-1)*3*natom) ! use xs as temporary work array
                  enddo
                  imode=0
                  ereverseuse=0.0d0 ! don't do reverse barrier stuff
                  call descent_check(natom,nmove,2,xs,x12,
     +                ityp,taxes,ietype,maxtyp,rcut,hyperratio,barrevev
     +                ,ereverseuse,isame,irotalign,itermax,gfac
     +                ,transcrit,itranscrit,gcrit,dvcrit,drcrit,ipr
     +                ,imode,ifundc)
c                  write(6,*) 'NEB_TAD: Debug, after dc 1,ifundc=',ifundc
                  if (isame.eq.1.or.isame.eq.2) then
                     if (ipr.ge.4) write(6,*)
     +                   '#NEB: intermediate minimum belongs to basin '
     +                   ,isame,', neglecting'
                     imin=0  
                  endif
               endif

            endif
         endif
      enddo
      call vecmov(x(1+(imax-1)*3*natom),xs,3*natom)
      if (imax.gt.1.and.imax.lt.nimage) then
         call vecmov(x(1+(imax-2)*3*natom),xleft,3*natom)
         call vecmov(x(1+(imax)*3*natom),xright,3*natom)
      endif
      if (imin.ne.0) call vecmov(x(1+(imin-1)*3*natom),x2,3*natom)
c      do j=1,3*natom
c         xs(j)=x(j+(imax-1)*3*natom)
c         if (imin.ne.0) x2(j)=x(j+(imin-1)*3*natom)
c      enddo
      es=emax
c----------------------------------------------------------------------
c no real saddle, "saddle" is one of the endpoints
      if (imax.eq.1.or.imax.eq.nimage.and.imin.eq.0) then
         if (ipr.ge.4) write(6,*)
     +       'NEB: note: no saddle in chain, endpoint is highest image'
c         if (imax.eq.1) stop 'NEB: imax.eq.1... not expected'
c xxx descent_check might needed tighter tolerances which need to be
c xxx passed, so means changes in call throughout code
c         imode=0
c xxx no saddle was found.  we could do a descent check here with
C tighter tolerance to see if it was a real transition
c         call descent_check(natom,nmove,1,x2,x1,
c     +       ityp,taxes,ietype,maxtyp,rcut,hyperratio,isame
c     +       ,irotalign,itermax*10,gfac,transcrit,itranscrit,gcrit*.1
c     +       ,dvcrit*.1,drcrit*.1,ipr,imode)
c         if (isame.eq.1) then
c xxx it was not a real transition, should have new error to send back
C to above routine to send to TAD, to reset ijump in mid stream.
C however, leak check will fail, so need to skip leak check as well.
c            do j=1,3*natom
c               x2(j)=x1(j)
c            enddo
c xxx if isame.eq.0, should set ierr=2, return, and neb will do denser
C neb to find saddle.
         if (imax.eq.nimage) ierr=2 ! only complain if end point is highest image as otherwise could get stuck
         if (imax.eq.1) write(6,*)
     $        'NEB: initial point highest, moving on'
c     intermediate=1
         if (use_local_neb) then
             ityp(:) = ityp_g(:)
             do i = 1, natom
                 x1_save((l_to_g_map(i)-1)*3+1) = x1((i-1)*3+1)
                 x1_save((l_to_g_map(i)-1)*3+2) = x1((i-1)*3+2)
                 x1_save((l_to_g_map(i)-1)*3+3) = x1((i-1)*3+3)
                 x2_save((l_to_g_map(i)-1)*3+1) = x2((i-1)*3+1)
                 x2_save((l_to_g_map(i)-1)*3+2) = x2((i-1)*3+2)
                 x2_save((l_to_g_map(i)-1)*3+3) = x2((i-1)*3+3)
                 xs_save((l_to_g_map(i)-1)*3+1) = xs((i-1)*3+1)
                 xs_save((l_to_g_map(i)-1)*3+2) = xs((i-1)*3+2)
                 xs_save((l_to_g_map(i)-1)*3+3) = xs((i-1)*3+3)
                 tan_neb_save((l_to_g_map(i)-1)*3+1) =
     +               tan_neb((i-1)*3+1)
                 tan_neb_save((l_to_g_map(i)-1)*3+2) =
     +               tan_neb((i-1)*3+2)
                 tan_neb_save((l_to_g_map(i)-1)*3+3) =
     +               tan_neb((i-1)*3+3)
             enddo
             x1(:) = x1_save(:)
             x2(:) = x2_save(:)
             xs(:) = xs_save(:)
             natom = natom_g
             nmove = nmove_g
             ! Startup fresh lammps run
             if((ietype.ge.100).and.(ietype.le.200))then
               x1lmp(:) = 0.d0
               x1lmp(1:3*natom) = x1(1:3*natom)
               call MPI_BARRIER(force_comm,ier)
               call lmp_setup(natom,nmove,ntype,maxtyp,potnam
     +           ,ietype,x1lmp,taxes,ityp,rcut,amu
     +           ,templow,temphigh,.true.,dt)
               call MPI_BARRIER(force_comm,ier)
             endif
         endif
         es = es - e(1)
         return
      endif
c descent_check was done with more iteration, hopefully x2 is more
c refined and a new NEB will give a barrier         
c         if (isame.eq.0) then
c            goto 100
c         endif
c      endif
c----------------------------------------------------------------------

c      do j=1,3*natom
c         xs(j)=x(j+(imax-1)*3*natom)
c         if (imin.ne.0) x2(j)=x(j+(imin-1)*3*natom)
c      enddo
      if (imin.ne.0) intermediate=1
c      es=emax
      if(rank_f.eq.0)then
        if (imin.eq.0.and.ipr.ge.4) write(6,*) '#NEB: saddle: image ',imax
     +    ,' with energy ',es*ev,' eV'
      endif
      
      if (use_local_neb) then
         ityp(:) = ityp_g(:)
         do i = 1, natom
             x1_save((l_to_g_map(i)-1)*3+1) = x1((i-1)*3+1)
             x1_save((l_to_g_map(i)-1)*3+2) = x1((i-1)*3+2)
             x1_save((l_to_g_map(i)-1)*3+3) = x1((i-1)*3+3)
             x2_save((l_to_g_map(i)-1)*3+1) = x2((i-1)*3+1)
             x2_save((l_to_g_map(i)-1)*3+2) = x2((i-1)*3+2)
             x2_save((l_to_g_map(i)-1)*3+3) = x2((i-1)*3+3)
             xs_save((l_to_g_map(i)-1)*3+1) = xs((i-1)*3+1)
             xs_save((l_to_g_map(i)-1)*3+2) = xs((i-1)*3+2)
             xs_save((l_to_g_map(i)-1)*3+3) = xs((i-1)*3+3)
             tan_neb_save((l_to_g_map(i)-1)*3+1) =
     +               tan_neb((i-1)*3+1)
             tan_neb_save((l_to_g_map(i)-1)*3+2) =
     +               tan_neb((i-1)*3+2)
             tan_neb_save((l_to_g_map(i)-1)*3+3) =
     +               tan_neb((i-1)*3+3)
         enddo
         x1(:) = x1_save(:)
         x2(:) = x2_save(:)
         xs(:) = xs_save(:)
         natom = natom_g
         nmove = nmove_g
         ! Startup fresh lammps run
         if((ietype.ge.100).and.(ietype.le.200))then
           x1lmp(:) = 0.d0
           x1lmp(1:3*natom) = x1(1:3*natom)
           call MPI_BARRIER(force_comm,ier)
           call lmp_setup(natom,nmove,ntype,maxtyp,potnam
     +       ,ietype,x1lmp,taxes,ityp,rcut,amu
     +       ,templow,temphigh,.true.,dt)
           call MPI_BARRIER(force_comm,ier)
         endif
      endif
      es = es - e(1)
      return
      end
c**********************************************************************
      subroutine neb_write_chain(nimage,natom,x,f,e,finfo,lstore_neb
     +    ,filnam,taxes,itype,springk,ipr)

      use mod_mpi
      implicit real*8(a-h,o-z)
      dimension x(3*nimage*natom),f(3*nimage*natom),e(nimage)
     +    ,finfo(nimage*4),itype(natom),taxes(3)
      logical lstore_neb
      character*180 filnam
      
      ev=27.21d0
      emax=-1.0d10
      do i=1,nimage
         if (e(i).gt.emax) then
            emax=e(i)
            imax=i
         endif
      enddo

      if(rank_f.eq.0)then

      if (lstore_neb) then
         if (ipr.ge.6) write(6,31) 'writing NEB file: ',filnam
 31      format(2a30)
         iunit=27
         lf=len(filnam)
         call chrpak(filnam,lf,lfil)
         open(unit=iunit,file=filnam(1:lfil),form='formatted',status
     +       ='unknown')
         rewind iunit
         write(iunit,30) nimage,natom,0,1
         write(iunit,40) 'NEB from TAD'
         write(iunit,50) (taxes(k),k=1,3)
         write(iunit,51) springk
         write(iunit,52) (itype(i),0,i=1,natom)
         
 30      format(i4,i5,i5,i3)
 40      format(a50)
 50      format(3f20.8)
 51      format(f20.8)
 52      format(2i2)

         do i=1,nimage
            do j=1,natom
               write(iunit,50) (x(k+(j-1)*3+(i-1)*3*natom),k=1,3)
            enddo
         enddo

         rewind iunit
         close(unit=iunit)
      endif

      rcum=0.0d0
      write(6,*) '#NEB: All Energies in eV'
      write(6,3) 'n','delr','rcum','E','delE','|F|','Fr','Fs','angle'
      do i=1,nimage
         ismin=0
         if (i.gt.1.and.i.lt.nimage) then
            if (e(i).lt.e(i-1).and.e(i).lt.e(i+1)) ismin=1
         endif
         r=0.0d0
         if (i.eq.1) then
            write(6,1) i,0.0d0,0.0d0,e(1)*ev,0.0d0,0.0d0,0.0d0,0.0d0,0
     $           .0d0
         else
            do j=1,3*natom
               r=r+(x(j+(i-1)*3*natom)-x(j+(i-2)*3*natom))**2
            enddo
            r=sqrt(r)
            rcum=rcum+r
            if (i.ne.imax.and.ismin.ne.1) write(6,1) i,r,rcum,e(i)*ev
     +          ,(e(i)-e(1))*ev,finfo((i-1)*4+1)*ev,finfo((i-1)*4+2)*ev
     +          ,finfo((i-1)*4+3)*ev,finfo((i-1)*4+4)
            if (i.eq.imax) write(6,2) i,r,rcum,e(i)*ev,(e(i
     +          )-e(1))*ev,finfo((i-1)*4+1)*ev,finfo((i-1)*4+2)*ev
     +          ,finfo((i-1)*4+3)*ev,finfo((i-1)*4+4)
            if (ismin.eq.1) write(6,4) i,r,rcum,e(i)*ev,(e(i
     +          )-e(1))*ev,finfo((i-1)*4+1)*ev,finfo((i-1)*4+2)*ev
     +          ,finfo((i-1)*4+3)*ev,finfo((i-1)*4+4)
         endif
      enddo
      
      endif
      
 1    format('NEB: ',i3,2f8.3,f10.3,f8.3,3f12.5,f8.1)
 2    format('NEB:>',i3,2f8.3,f10.3,f8.3,3f12.5,f8.1,' <')
 4    format('NEB:-',i3,2f8.3,f10.3,f8.3,3f12.5,f8.1,' -')
 3    format('#NEB:',a3,2a8,a10,a8,3a12,a8)
      
      return
      end
c**********************************************************************
      subroutine neb_init_linear(nimage,natom,x1,x2,x,xhot,hotmix
     +                                                 ,q,nmove,lnoise)

      use mod_mpi
      implicit real*8(a-h,o-z)
      logical lnoise
      dimension x1(3*natom),x2(3*natom),x(3*natom*nimage),xhot(3*natom)
      dimension xnoise(3*natom)
      dimension q(natom*nimage)

      if(lnoise)then
         do j=1,natom*3
           if(j.le.nmove*3)then
             if(rank_f.eq.0) xnoise(j)=x1(j)+(0.1d0)*gasdev(0) !xnoise(j)=x1(j)+(0.5d0)*(prngen(0)*2.0d0-1.0d0)
           endif
         enddo
         call MPI_BCAST(xnoise,natom*3,MPI_REAL8,0,force_comm,ier)
         d1=0.0d0
         d2=0.0d0
         do i=1,natom*3
            d1=d1+(xnoise(i)-x1(i))**2
            d2=d2+(x2(i)-xnoise(i))**2
         enddo
         d1=sqrt(d1)
         d2=sqrt(d2)
      elseif (hotmix.gt.0.0d0) then
         d1=0.0d0
         d2=0.0d0
         do i=1,natom*3
            d1=d1+(xhot(i)-x1(i))**2
            d2=d2+(x2(i)-xhot(i))**2
         enddo
         d1=sqrt(d1)
         d2=sqrt(d2)
      endif

       do i=1,nimage
         if((hotmix.gt.0.0d0).or.(lnoise))then
           di=(d1+d2)*float(i-1)/float(nimage-1)
         endif
         do j=1,natom*3
            x(j+(i-1)*3*natom)=x1(j)
     +          +(x2(j)-x1(j))*float(i-1)/float(nimage-1)
            if((hotmix.gt.0.0d0).or.(lnoise))then
               if (di.le.d1) then
                  xx=x1(j)+(xhot(j)-x1(j))*di/d1
               else
                  xx=xhot(j)+(x2(j)-xhot(j))*(di-d1)/d2
               endif
               x(j+(i-1)*3*natom)=x(j+(i-1)*3*natom)*(1-hotmix)+xx
     +             *hotmix
            endif
         enddo
       enddo

       do i = 1,nimage
         do j=1,natom
               q(j+(i-1)*natom) = 0.d0
         enddo
       enddo

      return
      end
c**********************************************************************
      subroutine neb_widen_chain(nimage,natom,x,grad,e,nmove,ityp,taxes
     +    ,ietype,maxtyp,rcut,hyperratio,q)

      use mod_mpi
      implicit real*8(a-h,o-z)
c      parameter(mxatom=10000)
      include 'parameters.h'

      dimension rcut(maxtyp,maxtyp),ityp(natom),grad(3*natom*nimage)
     +    ,e(nimage),x(3*natom*nimage),q(natom*nimage)
      dimension taxes(3)
      dimension rcutuse(3,3)

      do i=1,3
         do j=1,3
crz            rcutuse(i,j)=1.05d0 ! ag
            rcutuse(i,j)=1.5d0 ! ag -- rz test
c            rcutuse(i,j)=2.0d0  ! alnewsandia
         enddo
      enddo
      
      ietypeuse=999
      imode=0
      gfac=-1.0d0 ! neb_steepest assumes force, not gradient

      do i=1,10000
         call neb_gcalc(nimage,natom,x,grad,e,nmove,ityp,taxes,ietypeuse
     +       ,maxtyp,rcutuse,hyperratio,imode,dvmax,q)
c         write(6,*) 'widen: grada: ',grad(1),grad(2),grad(3)
         call neb_steepest(nimage,natom,nmove,x,grad,gfac,dxmax)
c         write(6,*) 'widen: gradb: ',grad(1),grad(2),grad(3)
         inonzero=0
         do j=3*natom+1,3*natom*(nimage-1)
            if (grad(j).ne.0.0d0) then
               inonzero=inonzero+1
c               write(6,*) 'neb_widen_chain: grad: ',j,grad(j),gfac
            endif
         enddo
         if (inonzero.eq.0) goto 111
c         write(6,*) 'neb_widen_chain: ',i,inonzero,dxmax
      enddo

 111  continue
      if(rank_f.eq.0) 
     +    write(6,*) 'SpawnID ',SpawnID,
     +        ' #NEB: widened chain, took ',i,' iterations'

      return
      end

c**********************************************************************
      subroutine neb_steepest(nimage,natom,nmove,x,f,gfac,dxmax)

      implicit real*8(a-h,o-z)
      dimension x(3*nimage*natom),f(3*nimage*natom)

      rmaxstep=0.1d0

c      write(6,*) 'widen: gfac: ',gfac

      dxmax=0.0d0
      dxold=0.0d0
      do i=2,nimage-1
         k=3
         dx=0.0d0
         gmag=sqrt(vecdot(f(1+(i-1)*3*natom),f(1+(i-1)*3*natom),
     $        3*nmove))
         gmaxstepfact=1.0d0
         if (gfac*gmag.gt.rmaxstep) then
            gmaxstepfact=1.0d0/gmag/gfac*rmaxstep
         endif
c         if (gfac.lt.0.0d0)
c     +       write(6,*) 'neb_steepest: ',i,gmaxstepfact,gfac*gmag,gmag
c     +       ,dxmax
         do j=1,nmove*3
            k=k+1
            if (k.gt.3) then
               dxmax=max(dxmax,dx)
               k=1
               dxold=dx
               dx=0.0d0
            endif
            xx=x(j+(i-1)*3*natom)
            x(j+(i-1)*3*natom)=x(j+(i-1)*3*natom)+
     +          gfac*f(j+(i-1)*3*natom)*gmaxstepfact
            dx=dx+(x(j+(i-1)*3*natom)-xx)**2
         enddo
      enddo
      dxmax=sqrt(dxmax)
      
      return
      end
c**********************************************************************
      subroutine neb_quickmin(nimage,natom,nmove,x,f,v,gfac,dxmax)

      implicit real*8(a-h,o-z)
      dimension x(3*nimage*natom),f(3*nimage*natom),v(3*nimage*natom)
      
      rmaxstep=0.1d0

c xxx may need masses here

      dxmax=0.0d0
      do i=2,nimage-1
         vdotf=0.0d0
         fdotf=0.0d0
         do j=1,nmove*3
            if (v(j+(i-1)*3*natom)*f(j+(i-1)*3*natom).lt.0) v(j+(i-1)*3
     $           *natom)=0.0d0
            vdotf=vdotf+v(j+(i-1)*3*natom)*f(j+(i-1)*3*natom)
            fdotf=fdotf+f(j+(i-1)*3*natom)*f(j+(i-1)*3*natom)
         enddo
         vmag=0.0d0
         do j=1,nmove*3
            v(j+(i-1)*3*natom)=f(j+(i-1)*3*natom)*(1.0d0+vdotf/fdotf)
            vmag=vmag+v(j+(i-1)*3*natom)**2
         enddo
         vmag=sqrt(vmag)
         if (vmag*gfac.gt.rmaxstep) then ! max step size of rmaxstep
            do j=1,nmove*3
               v(j+(i-1)*3*natom)=v(j+(i-1)*3*natom)/vmag/gfac*rmaxstep 
            enddo
         endif
         k=3
         dx=0.0d0
         do j=1,nmove*3
            k=k+1
            if (k.gt.3) then
               dxmax=max(dxmax,dx)
               k=1
               dx=0.0d0
            endif
            xx=x(j+(i-1)*3*natom)
            x(j+(i-1)*3*natom)=x(j+(i-1)*3*natom)+
     +           gfac*v(j+(i-1)*3*natom)
            dx=dx+(x(j+(i-1)*3*natom)-xx)**2
         enddo
      enddo
      dxmax=sqrt(dxmax)
      
      return
      end

c**********************************************************************
      subroutine neb_fire(nimage,natom,nmove,x,f,v,dxmax,rmaxjump
     +           ,alpha,deltat,deltatmax,finc,fadec,nstep
     +           ,nmin2,alphaini,fdec)
      implicit real*8(a-h,o-z)
      dimension x(3*nimage*natom),f(3*nimage*natom),v(3*nimage*natom)
      dimension step(3*natom)
      dimension alpha(nimage),nstep(nimage),deltat(nimage)

c xxx may need masses here
     
      dxmax=0.0d0
      do i=2,nimage-1
         power=0.0d0
         fdotf=0.0d0
         vdotv=0.0d0
         do j=1,nmove*3
            power=power+v(j+(i-1)*3*natom)*f(j+(i-1)*3*natom)
            fdotf=fdotf+f(j+(i-1)*3*natom)*f(j+(i-1)*3*natom)
            vdotv=vdotv+v(j+(i-1)*3*natom)*v(j+(i-1)*3*natom)
         enddo
         vdotv=sqrt(vdotv)
         fdotf=sqrt(fdotf)
         if(power.gt.0.0) then
            do j=1,nmove*3
               v(j+(i-1)*3*natom)=(1-alpha(i))*v(j+(i-1)*3*natom)+
     +         alpha(i)*vdotv*f(j+(i-1)*3*natom)/fdotf
            enddo
            if(nstep(i).ge.nmin2) then
              deltat(i)=min(deltat(i)*finc,deltatmax)
              alpha(i)=alpha(i)*fadec
            endif
            nstep(i)=nstep(i)+1
         else
            do j=1,3*nmove
              v(j+(i-1)*3*natom)=0.0
            enddo
            alpha(i)=alphaini
            deltat(i)=deltat(i)*fdec
            nstep(i)=0
         endif

         stepmag=0.0
         do j=1,nmove*3
            v(j+(i-1)*3*natom)=v(j+(i-1)*3*natom)+deltat(i)*
     +           f(j+(i-1)*3*natom)
            step(j)=deltat(i)*v(j+(i-1)*3*natom)
            stepmag=stepmag+step(j)*step(j)
         enddo

         stepmag=sqrt(stepmag)
         if(stepmag.gt.rmaxjump) then
            do j=1,3*nmove
               step(j)=step(j)*rmaxjump/stepmag
               
            enddo
         endif
         
         k=3
         dx=0.0d0
         do j=1,nmove*3
            k=k+1
            if (k.gt.3) then
               dxmax=max(dxmax,dx)
               k=1
               dx=0.0d0
            endif
            xx=x(j+(i-1)*3*natom)
            x(j+(i-1)*3*natom)=x(j+(i-1)*3*natom)+step(j)
    
            dx=dx+step(j)**2
            
         enddo
      enddo
      dxmax=sqrt(dxmax)
      
      return
      end
      
c**********************************************************************
      subroutine neb_gcalc(nimage,natom,x,f,e,nmove,ityp,taxes
     +    ,ietype,maxtyp,rcut,hyperratio,imode,dvmax,q)

      use mod_mpi
      implicit real*8(a-h,o-z)
c      parameter(mxatom=10000)
      include 'parameters.h'

      dimension rcut(maxtyp,maxtyp),ityp(natom),f(3*natom*nimage)
     +    ,e(nimage),x(3*natom*nimage),q(natom*nimage)
      dimension xgrad(3*mxatom)
      dimension taxes(3)

      ifirst=2
      ilast=nimage-1
      if (imode.eq.1) then
         ifirst=1
         ilast=nimage
      endif

      dvmax=0.0d0
      
      do i=ifirst,ilast
         do j=1,natom
            tadcharge(j)=q(j+(i-1)*natom)
         enddo
         call gcalc(natom,nmove,x(1+(i-1)*3*natom),ityp,taxes
     +       ,ietype,maxtyp,rcut,energy,xgrad,hyperratio)
         ! If charge was allowed to change, update it...
         do j=1,natom
           q(j+(i-1)*natom) = tadcharge(j)
         enddo
c         if (rcut(1,1).eq.1.0d0) write(6,*) 'widen: gcalca: ',i,xgrad(1)
c     +       ,xgrad(2),xgrad(3)
         do j=1,nmove*3
            f(j+(i-1)*3*natom)=-xgrad(j)
         enddo
c         if (rcut(1,1).eq.1.0d0) write(6,*) 'widen: gcalcb: ',i,f(1+(i-1
c     +       )*3*natom),f(2+(i-1)*3*natom),f(3+(i-1)*3*natom)
         do j=nmove*3+1,natom*3
            f(j+(i-1)*3*natom)=0.0d0
         enddo
         if (imode.ne.1) then
            dv=abs(e(i)-energy)
            dvmax=max(dvmax,dv)
         endif
         e(i)=energy
      enddo

c      if (rcut(1,1).eq.1.0d0) write(6,*) 'widen: gcalcf: ',f(1),f(2),f(3
c     +    )
      
      return
      end
c**********************************************************************
      subroutine neb_nudge(nimage,natom,nmove,x,f,e,tanmax,springk,itan
     +    ,iclimb,finfo,fmax,idimer)

c input: nimage,natom,x,f,springk,itan,iclimb
c output: f,e
      
      implicit real*8(a-h,o-z)
c      parameter(mxatom=10000,mximage=42)
      include 'parameters.h'
      dimension x(3*natom*nimage),f(3*natom*nimage),e(nimage),
     +    tanmax(3*natom)
      dimension dprev(3*mxatom),dnext(3*mxatom),tan_neb(3*mxatom)
      dimension finfo(nimage*4),fcm(3)
      dimension angles(mximage)
      logical eprev,enext

      pi=2.0d0*acos(0.0d0)

      do i=1,natom*3
         tanmax(i)=0.0d0
      enddo
      
      eimax=-1.0d10
      nimax=0
      ni1max=0
      do i=1,nimage
         angles(i)=pi
         if (i.gt.1.and.i.lt.nimage) then
            if (e(i).gt.eimax) then
               eimax=e(i)
               nimax=i
            endif
            if (e(i).gt.e(i-1).and.e(i).gt.e(i+1).and.ni1max.eq.0)
     $          ni1max=i
            pdotn=0.0d0
            rpmag=0.0d0
            rnmag=0.0d0
            if (i.gt.1.and.i.lt.nimage) then
               do j=1,nmove
                  do l=1,3
                     jj=(j-1)*3+l
                     rtemp1=x(jj+(i-2)*3*natom)-x(jj+(i-1)*3*natom)
                     rtemp2=x(jj+(i)*3*natom)-x(jj+(i-1)*3*natom)
                     pdotn=pdotn+rtemp1*rtemp2
                     rpmag=rpmag+rtemp1**2
                     rnmag=rnmag+rtemp2**2
                  enddo
               enddo
               rads=pdotn/sqrt(rpmag)/sqrt(rnmag)
               if (rads.lt.-1.0d0) rads=-1.0d0
               if (rads.gt.1.0d0) rads=1.0d0
               angles(i)=acos(rads)
c                  write(6,632) im,rpmag,rnmag,pdotn,angles(im)
 632           format('NEBD: rpmag,rnmag,pdotn: ',i5,4f12.6)
               if (e(i).lt.e(i-1).and.e(i).lt.e(i+1)) angles(i)=pi
            endif
         endif
         do j=1,4
            finfo((i-1)*4+j)=0.0d0
         enddo
      enddo
      fmax=0.0d0
      
      do i=2,nimage-1
         rprev=0.0d0
         rnext=0.0d0
         tmag=0.0d0
         fdott=0.0d0
         pdotn=0.0d0
         eprev=e(i-1).gt.e(i)
         enext=e(i+1).gt.e(i)

c----------------------------------------------------------------------
c     next, previous distance
         do j=1,nmove*3
            dprev(j)=x(j+(i-2)*3*natom)-x(j+(i-1)*3*natom)
            dnext(j)=x(j+(i)*3*natom)-x(j+(i-1)*3*natom)
            rprev=rprev+dprev(j)**2
            rnext=rnext+dnext(j)**2
            pdotn=pdotn+dprev(j)*dnext(j)
         enddo
         rprev=1.0d0/sqrt(rprev)
         rnext=1.0d0/sqrt(rnext)
         pdotn=pdotn*rprev*rnext
         if (pdotn.lt.-1.0d0.and.pdotn.gt.-1.001d0) pdotn=-1.0d0
         pdotn=acos(pdotn)*180d0/3.14159d0

c----------------------------------------------------------------------
c     tangent
         e1=e(i-1)-e(i)
         e2=e(i+1)-e(i)
         ismin=0
         if (e(i).lt.e(i-1).and.e(i).lt.e(i+1)) ismin=1
         emin=min(abs(e1),abs(e2))
         emax=max(abs(e1),abs(e2))
         if (itan.eq.1.or.itan.eq.2) then
            if (eprev.neqv.enext) then
               if (eprev) then
                  do j=1,natom*3
                     tan_neb(j)=-dprev(j)
                  enddo
               else
                  do j=1,natom*3
                     tan_neb(j)=dnext(j)
                  enddo
               endif
            else
cc               if ((abs(e1).le.abs(e2).and.(.not.eprev)).or.(abs(e1).gt
cc     $              .abs(e2).and.eprev)) then
               if (itan.eq.1) then ! regular new tangent
                  if (e1.gt.e2) then
                     do j=1,natom*3
                        tan_neb(j)=dnext(j)*emin-dprev(j)*emax
                     enddo
                  else
                     do j=1,natom*3
                        tan_neb(j)=dnext(j)*emax-dprev(j)*emin
                     enddo
                  endif
               else if (itan.eq.2) then ! old tangent for extrema (angle bisection)
                  do j=1,nmove*3
                     tan_neb(j)=dnext(j)*rnext-dprev(j)*rprev
                  enddo
               endif
c               do j=1,nmove*3
c                  tan_neb(j)=dnext(j)*rnext-dprev(j)*rprev
c               enddo
            endif
         else
            do j=1,nmove*3
               tan_neb(j)=dnext(j)*rnext-dprev(j)*rprev
            enddo
         endif
      
c     dot products
         do j=1,nmove*3
            fdott=fdott+f(j+(i-1)*3*natom)*tan_neb(j)
            tmag=tmag+tan_neb(j)**2
         enddo

c if climbing image, strengthen springs on each side
         sright=1.0d0
         sleft=1.0d0
c         if (iclimb.eq.1.and.i+1.eq.nimax) sright=3.0d0
c         if (iclimb.eq.1.and.i-1.eq.nimax) sleft=3.0d0
         
c----------------------------------------------------------------------               
c experimental angle dependence of spring constant
c         iangleneb=1
         iangleneb=0
         if (iangleneb.eq.1) then
            springmax=4.0d0 ! using for nanodiamond
            springmax=10.0d0 ! using for Pu
            springmax=10.0d0 ! fragments?
            if (i.eq.2) then
               angaveminus=angles(i)
               angaveplus=(angles(i)+angles(i+1))/2.0d0
            else if (i.eq.nimage-1) then
               angaveminus=(angles(i-1)+angles(i))/2.0d0
               angaveplus=angles(i)
            else 
               angaveminus=(angles(i-1)+angles(i))/2.0d0
               angaveplus=(angles(i)+angles(i+1))/2.0d0
            endif
            
            sleft=springmax+(1.0d0-springmax)*angaveminus/pi
            sright=springmax+(1.0d0-springmax)*angaveplus/pi
c            write(6,631) i,sleft,sright,angles(i-1),angles(i),
c     +          angles(i+1)
 631        format('NEBD: sleft,sright: ',i5,5f12.6)
         endif
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c     modify force
         tmag=1.0d0/sqrt(tmag)
         fspdott=-springk*sleft/rprev+springk*sright/rnext
         fdott=fdott*tmag
         frealperp=0.0d0
         fspar=0.0d0
         fspdotreal=0.0d0
         fmag=0.0d0

         fcm(1)=0.0d0
         fcm(2)=0.0d0
         fcm(3)=0.0d0

         iclimbuse=0
         iquenchuse=0
         if (iclimb.eq.1.or.iclimb.eq.11) iclimbuse=1
         if (iclimb.eq.10.or.iclimb.eq.11) iquenchuse=1

c         if (ipr.ge.10)
c            write(6,*)'NEB: iclimb,iclimbuse,iquenchuse,ismin: ',iclimb
c     +       ,iclimbuse,iquenchuse,ismin

         do j=1,nmove*3
            tan_neb(j)=tan_neb(j)*tmag
            fspar=fspar+(fspdott*tan_neb(j))**2
            frealperp=frealperp+(f(j+(i-1)*3*natom)-fdott*tan_neb(j))**2
            fspdotreal=fspdotreal+(fspdott*tan_neb(j))*(f(j+(i-1)*3
     $           *natom)-fdott*tan_neb(j))

            if (i.eq.nimax.and.idimer.eq.0) tanmax(j)=tan_neb(j)
            if (i.eq.ni1max.and.idimer.eq.1) tanmax(j)=tan_neb(j)
            if (i.eq.nimax.and.iclimbuse.eq.1) then
               f(j+(i-1)*3*natom)=f(j+(i-1)*3*natom)-2.0d0*fdott
     $             *tan_neb(j)
            else if (ismin.ne.0.and.iquenchuse.eq.1) then
               f(j+(i-1)*3*natom)=f(j+(i-1)*3*natom)
            else if ((i.ne.nimax.or.iclimbuse.eq.0).and.
     +             (ismin.eq.0.or.iquenchuse.eq.0)) then   
               f(j+(i-1)*3*natom)=f(j+(i-1)*3*natom)+fspdott*tan_neb(j)
     +             -fdott*tan_neb(j)
            else
               write(6,*) 'STOP: NEB: no condition met!'
               stop 'NEB: no condition met!'
            endif

            fmag=fmag+f(j+(i-1)*3*natom)**2
            index=mod(j-1,3)+1
c            fcm(index)=fcm(index)+f(j+(i-1)*3*natom)
         enddo

         if (nmove.eq.natom) then
            do j=1,nmove*3
               index=mod(j-1,3)+1
               f(j+(i-1)*3*natom)=f(j+(i-1)*3*natom)-fcm(index)
            enddo
         endif

         finfo((i-1)*4+1)=sqrt(fmag)
         finfo((i-1)*4+2)=sqrt(frealperp)
         finfo((i-1)*4+3)=sqrt(fspar)
         finfo((i-1)*4+4)=pdotn
         fmax=max(sqrt(fmag),fmax)
      enddo

      return
      end
c**********************************************************************
      subroutine map3to1(n,a1,a2,a3,b)

      implicit real*8(a-h,o-z)

      dimension a1(n),a2(n),a3(n),b(3*n)

      do i=1,n
         b(1+(i-1)*3)=a1(i)
         b(2+(i-1)*3)=a2(i)
         b(3+(i-1)*3)=a3(i)
      enddo
      
      return
      end
c**********************************************************************
      subroutine map1to3(n,a1,a2,a3,b)

      implicit real*8(a-h,o-z)

      dimension a1(n),a2(n),a3(n),b(3*n)

      do i=1,n
         a1(i)=b(1+(i-1)*3)
         a2(i)=b(2+(i-1)*3)
         a3(i)=b(3+(i-1)*3)
      enddo
      
      return
      end
c**********************************************************************
      subroutine centerofmass_interface(x,natom,xcm)

      implicit real*8(a-h,o-z)
      
c      parameter(mxatom=10000)
      include 'parameters.h'
      dimension x(3*natom),xcm(3)

      dimension xx(mxatom),yy(mxatom),zz(mxatom)

      call map1to3(natom,xx,yy,zz,x)

      call centerofmass(xx,yy,zz,natom,xc,yc,zc)

      xcm(1)=xc
      xcm(2)=yc
      xcm(3)=zc

      return
      end
c**********************************************************************      
      subroutine neb_inter_min(e,niter,i,imin,nimage,eshallow,ierr
     +    ,xstates,nstate,barrevev,ereverse,natom,nmove,xmin,ityp,taxes
     +    ,ietype,maxtyp,rcut,hyperratio,isame,irotalign,itermax,gfac
     +    ,transcrit,itranscrit,gcrit,dvcrit,drcrit,ipr,imode,ifundc)

      use mod_mpi
      implicit real*8(a-h,o-z)

      dimension e(nimage)
      dimension xstates(nstate*natom*3),barrevev(nstate-1)
      dimension ityp(natom),taxes(3),rcut(maxtyp,maxtyp)

      ev=27.21d0

      ishallow=0
      imin=i
      iminhold=imin

      eleft=0.0d0
      eright=0.0d0
      do j=2,imin-1
         eleft=max(eleft,(e(j)-e(imin))*ev)
      enddo
      do j=imin+1,nimage-1
         eright=max(eright,(e(j)-e(imin))*ev)
      enddo
      if (eright.lt.eshallow.or.eleft.lt.eshallow) then
         imin=0
         ishallow=1
      endif

      if (imin.ne.0) then
c         write(6,*) 'debug: niter=',niter
         if (mod(niter,100).eq.0.or.niter.lt.0) then
            ifundcuse=ifundc
            if (ifundcuse.eq.5) ifundcuse=0
            call descent_check(natom,nmove,nstate,xmin,xstates,
     +          ityp,taxes,ietype,maxtyp,rcut,hyperratio,barrevev
     +          ,ereverse,isame,irotalign,itermax,gfac,transcrit
     +          ,itranscrit,gcrit,dvcrit,drcrit,ipr,imode,ifundcuse)
            if (isame.ne.0.and.isame.ne.1) then
               if (barrevev(isame-1).lt.ereverse) then
                  ishallow=2
                  imin=0
                  if (1.eq.0) then
c xxx for free systems, before vecmov, need rotalign or something...
                     call vecmov(xstates((isame-1)*3*natom+1),xmin,3
     +                   *natom)
                     if((ipr.ge.6).and.(rank_f.eq.0)) write(6,*)
     +                   'replaced NEB local "minimum" with state '
     +                   ,isame
                  endif
               endif
            endif
         else
            imin=0              ! don't worry about intermediate minimum unless did above descent_check
         endif
      endif

      if ((ipr.ge.6.and.(ishallow.eq.1.or.ishallow.eq.2)).or.(ipr.ge.6
     +    .and.ishallow.eq.0))then
       if(rank_f.eq.0)then
         write(6,*) 'NEB note: Intermediate minimum found!'
         write(6,*) 'NEB note: image: ',iminhold,' with energy ',
     $       e(i)*ev
         if (ishallow.eq.1) write(6,111) eleft,eright
         if (ishallow.eq.2) write(6,112) barrevev(isame-1)
       endif
 111   format('#NEB: note: shallow minimum; barriers: '
     +       ,2f12.5,' (eV); neglecting')
 112   format('#NEB: note: minimum recognized, reverse barrier '
     +       ,'too small: ',f12.5,' (eV); neglecting')
c         ierr=3
      endif

      if (ishallow.ne.0) ierr=3 ! intermediate minimum that should be ignored
      
      return
      end
