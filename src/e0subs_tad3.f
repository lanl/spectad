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


c e0subs_tad2.f - EAM energy and gradient routines     12/6/01
c  like ones in verlet_alone.f (parmdg7), but now split off from md routines



      subroutine eam_setup(ntype,maxtyp,potnam,iwrt,
     x                        ietype,rcut,amu,azero)

      implicit real*8(a-h,o-z)
c sets up eam potential - reads in interp files
 
c passed in:
c  ntype   = number of types expected  (passed in)
c  maxtyp  = max. no. of types - (important: sets 2-D dimension rcut array)
c returned:
c  ietype  = energy type index.  For EAM, ietype=0
c  rcut    = 2-D array of cutoffs
c  amu     = 1-D array of masses
c  azero   = 1-D array of equilibrium T=0 lattice constants

      dimension rcut(maxtyp,maxtyp),amu(maxtyp),azero(maxtyp)
      character*80 dirnam,potnam(maxtyp),filnam,filnam0,cmdstr !,str1

c get dirnam from an environment variable  (added 9/24/02)
      write(6,*) 'getting potsdir name from POTSDIR environment var'
      dirnam='  '
      call getenv('POTSDIR',dirnam)
      !dirnam='../src/POTENTIALS/'
      if(dirnam.eq.'  ') then
        !write(6,*)'ERROR: you must set the environment variable POTSDIR'
        !write(6,*)'       to the location of the EAM potential files'
        !stop 'POTSDIR environment variable not set'
        write(6,*)'Assuming potential is in WORKING directory.'
        dirnam='./'
      end if
      write(6,*) 'dirnam =',dirnam

      do 8 i=1,ntype
      call chrpak(dirnam,80,ld)
      call chrpak(potnam(i),80,lp)

      filnam=dirnam(1:ld)//potnam(i)(1:lp)//'.doc'
!      filnam0=dirnam(1:ld)//potnam(i)(1:lp)//'.doc'
!      filnam=trim(ADJUSTL(filstrtLP))//potnam(i)(1:lp)//'.doc'
!      cmdstr = 'cp '//trim(ADJUSTL(filnam0))//' '//
!     x  trim(ADJUSTL(filnam))

      open(unit=60+i,file=filnam,status='old',
     x  form='formatted',err=99)
     
      filnam=dirnam(1:ld)//potnam(i)(1:lp)//'.rho'
!      filnam0=dirnam(1:ld)//potnam(i)(1:lp)//'.rho'
!      filnam=trim(ADJUSTL(filstrtLP))//potnam(i)(1:lp)//'.rho'
!      cmdstr = 'cp '//trim(ADJUSTL(filnam0))//' '//
!     x  trim(ADJUSTL(filnam))
      
      write(6,*) 'filnam: ',filnam
      open(unit=70+i,file=filnam,status='old',
     x  form='unformatted',err=99)

      filnam=dirnam(1:ld)//potnam(i)(1:lp)//'.f'
!      filnam0=dirnam(1:ld)//potnam(i)(1:lp)//'.f'
!      filnam=trim(ADJUSTL(filstrtLP))//potnam(i)(1:lp)//'.f'
!      cmdstr = 'cp '//trim(ADJUSTL(filnam0))//' '//
!     x  trim(ADJUSTL(filnam))
      
      write(6,*) 'filnam: ',filnam
      open(unit=90+i,file=filnam,status='old',
     x  form='unformatted',err=99)
     
      filnam=dirnam(1:ld)//potnam(i)(1:lp)//'.phi'
!      filnam0=dirnam(1:ld)//potnam(i)(1:lp)//'.phi'
!      filnam=trim(ADJUSTL(filstrtLP))//potnam(i)(1:lp)//'.phi'
!      cmdstr = 'cp '//trim(ADJUSTL(filnam0))//' '//
!     x  trim(ADJUSTL(filnam))
      
      write(6,*) 'filnam: ',filnam
      open(unit=80+i*i,file=filnam,status='old',
     x  form='unformatted',err=99)
      do 6 j=i+1,ntype
      lpi=lp
      call chrpak(potnam(j),80,lpj)
      
      filnam=dirnam(1:ld)//potnam(i)(1:lpi)//potnam(j)(1:lpj)//'.phi'
!      filnam0=dirnam(1:ld)//potnam(i)(1:lpi)//potnam(j)(1:lpj)//'.phi'
!      filnam=trim(ADJUSTL(filstrtLP))//potnam(i)(1:lp)//'.phi'
!      cmdstr = 'cp '//trim(ADJUSTL(filnam0))//' '//
!     x  trim(ADJUSTL(filnam))
      
      open(unit=80+i*j,file=filnam,status='old',
     x  form='unformatted',err=99)
    6 continue
    8 continue

c read nec. stuff from files
c  read in title cards and atomic masses from unit 60+itype
      if(iwrt.ge.0) write(6,2) ntype
    2 format(' here are the title cards from the',i2,' potentials:')
      do 5 i=1,ntype
      iu=60+i
      rewind iu
      read(iu,3) title
      read(iu,*) amu(i)
      read(iu,*) azero(i)
    3 format(10a8)
      if(iwrt.ge.0) write(6,4) title
    4 format(1x,10a8)
    5 continue


c procure cutoff values from trp arrays
c these calls automatically read in the trp arrays if necessary
      do 10 i=1,ntype
      call trpdat(70+i,npt,x0,x1,rhoci,y1,yn,iflo,ifhi)
      call trpdat(90+i,npt,x0,x1,xn,y1,yn,iflo,ifhi)
      do 10 j=1,ntype
      call trpdat(70+j,npt,x0,x1,rhocj,y1,yn,iflo,ifhi)
      call trpdat(80+i*j,npt,x0,x1,phicij,y1,yn,iflo,ifhi)
      rcut(i,j)=max(rhoci,rhocj,phicij)
   10 continue

      if(iwrt.ge.0) then
        write(6,20) ntype
   20   format(' trp arrays for ntype =',i2,
     x         ' have been read from units 70+i,80+i*j,90+i')
      end if

      ietype=0

      return

   99 stop 'open error in eam_setup'
      end

      subroutine g0(natom,nmove,x,y,z,itype, 
     x                          tax,tay,taz,
     x                          ietype, maxtyp,rcut,
     x                          energy,grad, hyperratio)
c  grad routine specific for EAM (ietype=0)
c  this is not a super fast routine
      implicit real*8(a-h,o-z)
      save
      dimension x(natom),y(natom),z(natom),itype(natom)
      dimension rcut(maxtyp,maxtyp)
      dimension grad(*)
      include 'wrksubs_tad2.h' 
ccc      parameter(natmax=10 000,lmax=242 000)  ! changed 3/30/98 to allow shorter lwrk1
c      parameter(natmax=10 000)
      include 'parameters.h'
      parameter(lmax1=2*(lwrk1-7*natmax/2)/9-1 )
      parameter(lmaxh=lmax1/2)
      parameter(lmax=lmaxh*2)  ! this makes it an even number
      common/wrkcm1/iwrk1,jwrk1,rhoj(natmax),jtype(natmax),
     x                 rhotot(natmax),nneigh(natmax),markn(natmax),
     x         jatomn(lmax),dxn(lmax),dyn(lmax),dzn(lmax),rn(lmax),
     x         fill1(lwrk1-3*natmax-natmax/2-4*lmax-lmax/2)
      common/potcom/rct(5,5),tx,ty,tz,ax,ay,az        !kludge for now

      ! RZ -- Uncomment bellow to allow binned neigh listing:
      common/nlists/ilist

c      write(6,*) 'lmax: ',lmax

c kludge for now - move these into call
c ilist setted by C.W.Pao 02222008
      ilist=1
      tx=tax
      ty=tay
      tz=taz
      ax=tax/2.0d0
      ay=tay/2.0d0
      az=taz/2.0d0
      do i=1,5
      do j=1,5
        rct(i,j)=rcut(i,j)
      end do
      end do

      ntype=imax(itype,natom)

      hyperratio=1.0d0

      if(ntype.eq.1) then
        call g0faster(natom,nmove,x,y,z,itype, 
     x                          tax,tay,taz,
     x                          maxtyp,rcut,
     x                            energy,grad, hyperratio)
        return
      else if(ntype.ge.2) then
        call g0faster2(natom,nmove,x,y,z,itype,      ! xxx installed 10/25/02
     x                          tax,tay,taz,
     x                          maxtyp,rcut,
     x                            energy,grad, hyperratio)
        return
      end if

      if(natom.gt.natmax) stop 'g0: natom exceeds natmax'
      if(iwrk1.gt.0)  stop 'g0 - wrkcm1 conflict'

c  make up full set of neighbor lists

      pad=0.0d0
      call rlistr(natom,1,natom,x,y,z,itype, pad,
     x         nneigh,markn, jatomn,dxn,dyn,dzn,rn,lmax)

c  evaulate rhotot for all atoms, and sum up energy

      energy=0.0d0
      do 20 i=1,natom
      it=itype(i)
      nn=nneigh(i)
      m=markn(i)+1
      call igthr(itype,jatomn(m),nn,jtype) ! gather jtype array
      call rhatm3(it,nn,rn(m),jtype, rho,ei)
      rhotot(i)=rho
      energy=energy+ei
   20 continue

c  calculate gradient for all moving atoms

      do 50 i=1,nmove
      ix=3*(i-1)+1
      iy=ix+1
      iz=ix+2
      it=itype(i)
      nn=nneigh(i)
      m=markn(i)+1
      call igthr(itype,jatomn(m),nn,jtype) ! gather jtype array
      call rgthr(rhotot,jatomn(m),nn,rhoj) ! gather rhoj array
      rhoi=rhotot(i)
      call gatom(it,rhoi,nn,dxn(m),dyn(m),dzn(m),rn(m),jtype,rhoj,
     x                                                   gx,gy,gz)
      grad(ix)=gx
      grad(iy)=gy
      grad(iz)=gz
   50 continue

      iwrk1=0
      return
      end


      subroutine rlistr(natom,iatom1,iatom2,x,y,z,itype, pad,
c                              ! supplied
     x          nneigh,markn,jatomn,dxn,dyn,dzn,rn,lmax)
c                              ! returned
c  like rlist_, but only computes neighbor lists for range of atoms iatom1,iatom2
c  like rlistn, but uses z-boxing info to improve speed
c  general neighbor-listing routine
c  sets up a full neighbor list for each of the natom atoms
c  user supplies:  natom, x(),y(),z(),itype(),
c                  lmax=maximum length for jatomn,dxn,dyn,dzn,rn arrays
c                  pad = extra distance (padding) to add to rcut
c                        (necessary if using later calls to rlistu)
c  rlistn returns:
c                  nneigh(natom) = no. of neighbors for each atom
c                  markn(natom) = pointer array for jatomn,dxn,dyn,dzn,rn
c                  jatomn(),dxn(),dyn(),dzn(),rn() - length<=lmax
c           define: m=markn(i)
c              jatomn(m+j) = j'th neighbor of atom i
c              dxn(m+j) = x(i)-x(j)  (dyn,dzn analagous)
c              rn(m+j) = distance from atom i to j'th neighbor of atom i
c
      implicit real*8(a-h,o-z)
      save
      dimension x(natom),y(natom),z(natom),itype(natom)
      dimension nneigh(natom),markn(natom)
      dimension jatomn(lmax),dxn(lmax),dyn(lmax),dzn(lmax),rn(lmax)
      common/potcom/rcut(5,5),tax,tay,taz,ax,ay,az
      common/nlists/ilist

      if(ilist.gt.0) then
        call rlistb(natom,iatom1,iatom2,x,y,z,itype, pad,
     x            nneigh,markn,jatomn,dxn,dyn,dzn,rn,lmax)
        return
      end if

        k=0
        do 190 i=iatom1,iatom2
        markn(i)=k
        it=itype(i)
        do 170 j=1,natom
        if(j.eq.i) go to 170
        dx=x(i)-x(j)
        if(dx.gt.ax) dx=dx-tax
        if(dx.lt.-ax) dx=dx+tax
        dy=y(i)-y(j)
        if(dy.gt.ay) dy=dy-tay
        if(dy.lt.-ay) dy=dy+tay
        dz=z(i)-z(j)
        if(dz.gt.az) dz=dz-taz
        if(dz.lt.-az) dz=dz+taz
        r2=dx*dx+dy*dy+dz*dz
        if(r2.gt.(rcut(it,itype(j))+pad)**2) go to 170
          k=k+1
          jatomn(k)=j
          dxn(k)=dx
          dyn(k)=dy
          dzn(k)=dz
          rn(k)=sqrt(r2)
  170   continue
        if(k.gt.lmax) then
           write(6,*) 'STOP: out of space in rlistr: ',k,lmax
           stop 'abort - out of space in rlistr'
        endif
        nneigh(i)=k-markn(i)
  190   continue

      return
      end

      subroutine rlistb(natom,iatom1,iatom2,x,y,z,itype, pad,
c                              ! supplied
     x          nneigh,markn,jatomn,dxn,dyn,dzn,rn,lmax)
c                              ! returned
c  computes full neighbor list each atom in range iatom1,iatom2
c  uses 3-D boxing to increase speed
c  user supplies:  natom, x(),y(),z(),itype(),
c                  lmax=maximum length for jatomn,dxn,dyn,dzn,rn arrays
c                  pad = extra distance (padding) to add to rcut
c                        (necessary if using later calls to rlistu)
c  returned:
c                  nneigh(natom) = no. of neighbors for each atom
c                  markn(natom) = pointer array for jatomn,dxn,dyn,dzn,rn
c                  jatomn(),dxn(),dyn(),dzn(),rn() - length<=lmax
c           define: m=markn(i)
c              jatomn(m+j) = j'th neighbor of atom i
c              dxn(m+j) = x(i)-x(j)  (dyn,dzn analagous)
c              rn(m+j) = distance from atom i to j'th neighbor of atom i
c
      implicit real*8(a-h,o-z)
      save
      dimension x(natom),y(natom),z(natom),itype(natom)
      dimension nneigh(natom),markn(natom)
      dimension jatomn(lmax),dxn(lmax),dyn(lmax),dzn(lmax),rn(lmax)
      common/potcom/rcut(5,5),tax,tay,taz,ax,ay,az
      include 'wrksubs_tad2.h' 
      parameter(mxijk=10000,mxpont=1000000)
      include 'parameters.h'
      common/wrkcm10/iwrk10,jwrk10,
     x            ixbox(mxatom),iybox(mxatom),izbox(mxatom),
     x            nijk(mxijk),ipoint(mxpont),
     x            fill10(lwrk10 - 3*mxatom/2 - mxijk/2 - mxpont/2)
      common/nlists/ilist
c  ilist=0 - rlist calls stay there
c  ilist=1 - rlist calls diverted here (rlistb)
c  ilist=2 - rlist calls diverted here (rlistb) with no boxing

      if(iwrk10.le.0 .and. ilist.ne.2) then
        iwrk10=1132
        ifbox=1
      else
        ifbox=0
        go to 234
      end if

c begin 3-D boxing scheme

c  determine number of boxes in each of x,y,z directions

c   make each box a little larger, as nec., to exactly fill tax length
      call dmnmx(rcut,9,rdum1,rcmax,idum1,idum2)
      nreach=1
      nxbox =  int(tax/(rcmax/float(nreach)))
      nybox =  int(tay/(rcmax/float(nreach)))
      nzbox =  int(taz/(rcmax/float(nreach)))

      ! RZ -- Min number of boxes in Z-Dir:
      ! USE FOR DEPOSITION
      nzbox =  min(1+2*nreach,int(taz/(rcmax/float(nreach))))

      !! RZ -- Use 1 box in X & Y Dirs:
      !! USE FOR SINGLE NANOWIRE
      !nxbox =  1
      !nybox =  1

      !! RZ -- Use 1 box in X-Dir:
      !! USE FOR CROSSED NANOWIRES
      !nxbox =  1

      if(min(nxbox,nybox,nzbox).lt.1+2*nreach) then
        ifbox=0
        go to 234
      end if
c    xbox = x length of the box  (a little larger than rcmax/nreach)
      xbox = tax/float(nxbox)
      ybox = tay/float(nybox)
      zbox = taz/float(nzbox)
      itest=nxbox*nybox*nzbox
      if(itest.gt.mxijk .or. itest.lt.0) then
          ifbox = 0
          write(6,*)  'note =- not able to box in rlistb -itest= ',itest
          go to 234
c         stop 'increase mxijk in rlistb'
      end if

c  determine which box each atom is in to find max in any one box

      do 10 ijk=1,nxbox*nybox*nzbox
   10 nijk(ijk)=0

c    following scheme good unless an atom coord is more negative than 10*tax
      do 20 i=1,natom
      ixbox(i)=mod(x(i)+100.0d0*tax,tax)/xbox + 1
      iybox(i)=mod(y(i)+100.0d0*tay,tay)/ybox + 1
      izbox(i)=mod(z(i)+100.0d0*taz,taz)/zbox + 1
      if(ixbox(i).lt.1 .or. ixbox(i).gt.nxbox) stop 'fix x boxing'
      if(iybox(i).lt.1 .or. iybox(i).gt.nybox) stop 'fix y boxing'
      if(izbox(i).lt.1 .or. izbox(i).gt.nzbox) stop 'fix z boxing'
      ijk = (izbox(i)-1)*nxbox*nybox + (iybox(i)-1)*nxbox + ixbox(i)
      nijk(ijk)=nijk(ijk)+1
   20 continue
      nbmax = imax(nijk,nxbox*nybox*nzbox)
      if(nbmax*nxbox*nybox*nzbox.gt.mxpont)stop'rlistb: increase mxpont'

      do 25 ijk=1,nxbox*nybox*nzbox
      nijk(ijk)=0
   25 continue

c  fill ipoint array so atoms in that box can be found

c     ipoint array is ordered as follows:
c      fastest index:  icell      max=nbmax
c                       x         max=nxbox
c                       y         max=nybox
c      slowest index:   z         max=nzbox

      do 30 i=1,natom
      ix = ixbox(i)
      iy = iybox(i)
      iz = izbox(i)
      ijk = (iz-1)*nxbox*nybox + (iy-1)*nxbox + ix
      nijk(ijk)=nijk(ijk)+1
      icell=nijk(ijk)
      ijkm = (ijk-1)*nbmax + icell
      ipoint(ijkm)=i
   30 continue

c now have the following arrays of interest:

c     ipoint(...)  given box, can look up each atom in that box
c     nijk(ijk) = number of atoms in box ijk
c     ixbox(i)  =  x box number for atom i
c     iybox(i)  =  y box number for atom i
c     izbox(i)  =  z box number for atom i

c make neighbor list
c  note that this code is only good if atoms are within roughly one tax of origin

 234  continue
      if(ifbox.ne.0) then
c        write(6,*) 'will perform 3-D boxing for fast neighbor list'
        k=0
        do 90 i=iatom1,iatom2
        markn(i)=k
        it=itype(i)
ccc       do 70 j=1,natom

        do 80 ixraw=ixbox(i)-nreach,ixbox(i)+nreach
        do 70 iyraw=iybox(i)-nreach,iybox(i)+nreach
        do 60 izraw=izbox(i)-nreach,izbox(i)+nreach
        ix =   mod(ixraw + nxbox-1,nxbox) + 1
        iy =   mod(iyraw + nybox-1,nybox) + 1
        iz =   mod(izraw + nzbox-1,nzbox) + 1
        ijk = (iz-1)*nxbox*nybox + (iy-1)*nxbox + ix
        do 50 icell=1,nijk(ijk)
        ijkm = (ijk-1)*nbmax + icell
        j=ipoint(ijkm)
        if(j.eq.i) go to 50
        dx=x(i)-x(j)
        if(dx.gt.ax) then
           dx=dx-tax
        else if(dx.lt.-ax) then
           dx=dx+tax
        end if
        dy=y(i)-y(j)
        if(dy.gt.ay) then
          dy=dy-tay
        else if(dy.lt.-ay) then
          dy=dy+tay
        end if
        dz=z(i)-z(j)
        if(dz.gt.az) then
          dz=dz-taz
        else if(dz.lt.-az) then
          dz=dz+taz
        end if
        r2=dx*dx+dy*dy+dz*dz
        if(r2.gt.(rcut(it,itype(j))+pad)**2) go to 50
          k=k+1
          jatomn(k)=j
          dxn(k)=dx
          dyn(k)=dy
          dzn(k)=dz
          rn(k)=sqrt(r2)
   50   continue
   60   continue
   70   continue
   80   continue
        if(k.gt.lmax) stop 'abort - out of space in rlistb'
        nneigh(i)=k-markn(i)
   90   continue

      else if(ifbox.eq.0) then
c        write(6,*) 'NOTE - skipping 3-D boxing: full neighbor list work'
        k=0
        do 190 i=iatom1,iatom2
        markn(i)=k
        it=itype(i)
        do 170 j=1,natom
        if(j.eq.i) go to 170
        dx=x(i)-x(j)
        if(dx.gt.ax) dx=dx-tax
        if(dx.lt.-ax) dx=dx+tax
        dy=y(i)-y(j)
        if(dy.gt.ay) dy=dy-tay
        if(dy.lt.-ay) dy=dy+tay
        dz=z(i)-z(j)
        if(dz.gt.az) dz=dz-taz
        if(dz.lt.-az) dz=dz+taz
        r2=dx*dx+dy*dy+dz*dz
        if(r2.gt.(rcut(it,itype(j))+pad)**2) go to 170
          k=k+1
          jatomn(k)=j
          dxn(k)=dx
          dyn(k)=dy
          dzn(k)=dz
          rn(k)=sqrt(r2)
  170   continue
        if(k.gt.lmax) stop 'abort - out of space in nonboxing rlistb'
        nneigh(i)=k-markn(i)
  190   continue
      end if

      if(iwrk10.eq.1132) iwrk10=0

      return
      end

      subroutine gatom(it,rhoi,nneigh,dx,dy,dz,r,itype,rhotot,
c                              ! supplied
     x                                           gx,gy,gz)
c                              ! returned
c  this routine computes the gradient (gx,gy,gz) on *one* atom.
c  calling routine supplies the type of this atom (ityp), the
c  electron density at this atom (rhoi),the number of
c  neighbors (nneigh), the relative displacements to each neighbor
c  (dx(nneigh),dy(nneigh),dz(nneigh),r(nneigh)), the type of each
c  neighbor (itype(nneigh)), and the electron density at each
c  neighbor atom (rhotot(nneigh)).
c
      implicit real*8(a-h,o-z)
      save
c     dimension dx(nneigh),dy(nneigh),dz(nneigh),r(nneigh)
c     dimension itype(nneigh),rhotot(nneigh)
      dimension dx(*),dy(*),dz(*),r(*)
      dimension itype(*),rhotot(*)
c
c  statement functions for interpolation fits
c  note:  these assume that sgcon has been performed to eliminate s() and g()
      rho(rr,it)=trpfun(rr,0,70+it)
      tgrho(rr,it)=trpfun(rr,11,70+it)
      phi(rr,it,jt)=trpfun(rr,0,80+it*jt)
      tgphi(rr,it,jt)=trpfun(rr,1,80+it*jt)/rr
      f(x,it)=trpfun(x,0,90+it)
      fp(x,it)=trpfun(x,1,90+it)
c
c  compute gradient on this atom
c
      gx=0.0d0
      gy=0.0d0
      gz=0.0d0
      fpi=fp(rhoi,it)
      do 20 j=1,nneigh
      jt=itype(j)
      tgrad=tgphi(r(j),it,jt) + fpi*tgrho(r(j),jt)
     x                         + fp(rhotot(j),jt)*tgrho(r(j),it)
      gx=gx+tgrad*dx(j)
      gy=gy+tgrad*dy(j)
      gz=gz+tgrad*dz(j)
   20 continue

      return
      end

      subroutine rgthr(arin,ip,n,arout)
c  gather n elements from arin, using pointer array ip, & put them in arout
      implicit real*8(a-h,o-z)
c     dimension arin(*),arout(n),ip(n)
      dimension arin(*),arout(*),ip(*)
c
      do 10 i=1,n
   10 arout(i)=arin(ip(i))
      return
      end

      subroutine igthr(irin,ip,n,irout)
c  gather n elements from irin, using pointer array ip, & put them in irout
      implicit real*8(a-h,o-z)
c     dimension irin(*),irout(n),ip(n)
      dimension irin(*),irout(*),ip(*)
c
      do 10 i=1,n
   10 irout(i)=irin(ip(i))
      return
      end

      subroutine rhatm3(ityp,nneigh,r,itype,
c                              ! supplied
     x                           rhotot,energy)
c                              ! returned
c like rhatom, but uses 3-point interpolation for higher accuracy
c  this routine computes the electron density at *one* atom.
c  the calling routine supplies the type of this atom (ityp), the number
c  of neighbors (nneigh), the distance to each neighbor (r(nneigh)),
c  and the type of each neighbor (itype(nneigh)).
c  this routine returns the energy of this atom (energy),
c  and the electron density (rhotot).
c
      implicit real*8(a-h,o-z)
      save
c     dimension r(nneigh),itype(nneigh)
      dimension r(*),itype(*)

      rho(rr,it)=trpfn3(rr,0,70+it)
      phi(rr,it,jt)=trpfn3(rr,0,80+it*jt)
      f(x,it)=trpfn3(x,0,90+it)


      epair=0.0d0
      it=ityp
      rhotot=0.0d0
      do 40 j=1,nneigh
      jt=itype(j)
      epair=epair + phi(r(j),it,jt)
      rhotot=rhotot+rho(r(j),jt)
   40 continue
      energy = f(rhotot,it) + 0.5d0*epair

      return
      end

      subroutine g0faster(natom,nmove,x,y,z,itype, 
     x                          tax,tay,taz,
     x                          maxtyp,rcut,
     x                               energy,grad, hyperratio)
c like g0fast, but does simple update if possible (if atoms have not moved much)
      implicit real*8(a-h,o-z)
      save
c      save padhold,taxhold,tayhold,tazhold,rcuthold
      dimension x(natom),y(natom),z(natom),itype(natom)
      dimension rcut(maxtyp,maxtyp)
      dimension grad(*)
      include 'wrksubs_tad2.h'
ccc      parameter(natmax=10 000,lmax=242 000)
c      parameter(natmax=10 000)
      include 'parameters.h'
      parameter(lmax1=2*(lwrk1-11*natmax)/9-1 )
      parameter(lmaxh=lmax1/2)
      parameter(lmax=lmaxh*2)  ! this makes it an even number
      common/wrkcm1/iwrk1,jwrk1,rhoj(natmax),jtype(natmax),fpi(natmax),
     x     xhold(natmax),yhold(natmax),zhold(natmax),itypehold(natmax),
     x                 dxnow(natmax),dynow(natmax),dznow(natmax),
     x                 rhotot(natmax),nneigh(natmax),markn(natmax),
     x         jatomn(lmax),dxn(lmax),dyn(lmax),dzn(lmax),rn(lmax),
     x         fill1(lwrk1 - 11*natmax - 9*lmax/2)
      common/mdcom/dt,tnow,mdwrt,mdfast,mdcool,mddum,padmd
      common/uflags/iuflag(100)
ccccccccc      common/potcom/rcut(3,3),tax,tay,taz,ax,ay,az
        common/holdtr/ a(40000),dum(80000)
        common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                         mark1(99),mark2(99),
     x                         iflo(99),ifhi(99),last
      dimension rcuthold(5,5)
      rho(rr,it)=trpfn3(rr,0,70+it)
      phi(rr,it,jt)=trpfn3(rr,0,80+it*jt)
      f(xx,it)=trpfn3(xx,0,90+it)
      tgrho(rr,it)=trpfun(rr,11,70+it)
      tgphi(rr,it,jt)=trpfun(rr,1,80+it*jt)/rr
      fp(xx,it)=trpfun(xx,1,90+it)
      
      if(natom.gt.natmax) stop 'g0faster abort: natom.gt.natmax'
c check that there is only one atom type
      ntype=imax(itype,natom)
      if(ntype.gt.1) stop 'g0faster abort: ntype.gt.1 not allowed'
c check that rho(r) and phi(r) are tabled the same
      if(x0(71).ne.x0(81)) stop 'g0faster abort - x071.ne.x081'
      if(dxi(71).ne.dxi(81)) stop 'g0faster abort - dxi71.ne.dxi81'
      if(npt(71).ne.npt(81)) stop 'g0faster abort - npt71.ne.npt81'
      
      pad=float(iuflag(92))/100.0d0

c check that cutoff+2*pad is less than half of tax,tay or taz.
c (only 1*pad necessary?)
      rcmax=0.0d0
      do 5 it=1,ntype
      do 5 jt=1,ntype
      rcmax=max(rcmax,rcut(it,jt))
    5 continue
      if(rcmax+2.d0*pad.gt.tax/2.0d0) stop 'g0faster - tax too small'
      if(rcmax+2.d0*pad.gt.tay/2.0d0) stop 'g0faster - tay too small'
      if(rcmax+2.d0*pad.gt.taz/2.0d0) stop 'g0faster - taz too small'
      
c check maximum atom move since last full update
      inew=0
      dr2max=0.0d0
      mxdeltyp=0
      do 10 i=1,natom
      dxnow(i)=x(i)-xhold(i)
      dynow(i)=y(i)-yhold(i)
      dznow(i)=z(i)-zhold(i)
      dr2=dxnow(i)**2+dynow(i)**2+dznow(i)**2
      dr2max=max(dr2max,dr2)
      mxdeltyp=max(mxdeltyp,abs(itype(i)-itypehold(i)))
   10 continue
      
      if(sqrt(dr2max).gt.pad/2.0d0 .or. mxdeltyp.ne.0) inew=1
      if(pad.ne.padhold) inew=1
      if(taxhold.ne.tax.or.tayhold.ne.tay.or.tazhold.ne.taz) inew=1
      do i=1,ntype
      do j=1,ntype
        if(rcut(i,j).ne.rcuthold(i,j)) inew=1
        rcuthold(i,j)=rcut(i,j)
      end do
      end do
      
      if(iwrk1.ne.-7778) inew=1
      call wrkclaimv(1,1,'g0faster')

c make fresh set of neighbor lists over all atoms if necessary
      if(inew.eq.1) then
        call rlistr_half(natom,1,natom,x,y,z,itype, pad,
     x            nneigh,markn, jatomn,dxn,dyn,dzn,rn, lmax)
        do 12 i=1,natom
        xhold(i)=x(i)
        yhold(i)=y(i)
        zhold(i)=z(i)
        itypehold(i)=itype(i)
        dxnow(i)=0.0d0
        dynow(i)=0.0d0
        dznow(i)=0.0d0
   12   continue
        padhold=pad
        taxhold=tax
        tayhold=tay
        tazhold=taz
      end if
      
c  evaulate rhotot and fpi for all atoms and sum up energy
c  works on a half list

      call dmzer(rhotot,natom)  
      epair=0.0d0

      do 40 i=1,natom
      rhoi=0.0d0
      do 20 jj=1,nneigh(i)
      mj=markn(i)+jj
      j=jatomn(mj)
        dxij=dxn(mj)+dxnow(i)-dxnow(j)     !update
        dyij=dyn(mj)+dynow(i)-dynow(j)     !update
        dzij=dzn(mj)+dznow(i)-dznow(j)     !update
        rij=sqrt(dxij**2+dyij**2+dzij**2)     !update
        rn(mj)=rij     !update

        rm=dxi(71)*(rij-x0(71))
        m0=rm
        m0=min(m0,npt(71)-1)
        p=min(rm-float(m0),1.0d0)

        m=mark(71)+m0

        rho3=a(m)+0.5d0*p*(p*(a(m+1)+a(m-1)-a(m)-a(m)) +a(m+1)-a(m-1))
cc          rho3=rho(rij,1)  ! for testing
        rhotot(i)=rhotot(i)+rho3
        rhotot(j)=rhotot(j)+rho3

        m=mark(81)+m0
cc          phi3=phi(r,1,1)  ! for testing
        phi3=a(m)+0.5d0*p*(p*(a(m+1)+a(m-1)-a(m)-a(m)) +a(m+1)-a(m-1))
        epair = epair + phi3

   20 continue
   40 continue
      
      embed=0.0d0
      do 60 i=1,natom
      embed = embed + f(rhotot(i),1)
      fpi(i)=fp(rhotot(i),1)
   60 continue
ccc      energy = embed + 0.5d0*epair    ! full list
      energy = embed + epair       ! half-neigh-list code
      
c calculate gradient over all moving atoms

      call dmzer(grad,3*natom) 
      do 80 i=1,nmove   
      ix=3*(i-1)+1
      iy=ix+1
      iz=ix+2
      do 70 jj=1,nneigh(i)
      mj=markn(i)+jj
      rij=rn(mj)    ! update (rn(mj) filled above to save the sqrt)
      if(rij.gt.rcut(1,1)) go to 70   ! 2/23

      j=jatomn(mj)
        jx=3*(j-1)+1   ! half-neigh-list code
        jy=jx+1        ! half-neigh-list code
        jz=jx+2        ! half-neigh-list code

        dxij=dxn(mj)+dxnow(i)-dxnow(j)     !update
        dyij=dyn(mj)+dynow(i)-dynow(j)     !update
        dzij=dzn(mj)+dznow(i)-dznow(j)     !update

c  2-pt interp for tgphi and tgrho
        rm=dxi(81)*(rij-x0(81))
        m0=rm
        m0=min(m0,npt(81)-1)
        p=min(rm-float(m0),1.0d0)

        m=mark(81)+m0
        gi=a(m+1)-a(m-1)
        gi1=a(m+2)-a(m)
        slope = ( gi + p*(gi1-gi) )*dxi(81)*0.5d0
        tgphi3=slope/rij

        m=mark(71)+m0
        gi=a(m+1)-a(m-1)
        gi1=a(m+2)-a(m)
        slope = ( gi + p*(gi1-gi) )*dxi(71)*0.5d0
        tgrho3=slope/rij

ccc        tgphi3=tgphi(rij,1,1)     ! for testing
ccc        tgrho3=tgrho(rij,1)       ! for testing
        tgrad = tgphi3 + (fpi(i)+fpi(j))*tgrho3
        grad(ix)=grad(ix)+tgrad*dxij
        grad(iy)=grad(iy)+tgrad*dyij
        grad(iz)=grad(iz)+tgrad*dzij
        grad(jx)=grad(jx)-tgrad*dxij  ! half-neigh-list code
        grad(jy)=grad(jy)-tgrad*dyij  ! half-neigh-list code
        grad(jz)=grad(jz)-tgrad*dzij  ! half-neigh-list code
   70 continue
   80 continue
      
      call wrkrelease(1)
      iwrk1=-7778
      return
      end


      subroutine rlistr_half(natom,iatom1,iatom2,x,y,z,itype, pad,
     x          nneigh,markn,jatomn,dxn,dyn,dzn,rn,lmax)
c like rlistr, but returns a half list (j.gt.i)

      implicit real*8(a-h,o-z)
      save
      dimension x(natom),y(natom),z(natom),itype(natom)
      dimension nneigh(natom),markn(natom)
      dimension jatomn(lmax),dxn(lmax),dyn(lmax),dzn(lmax),rn(lmax)
      common/potcom/rcut(5,5),tax,tay,taz,ax,ay,az
      common/nlists/ilist


      if(ilist.gt.0) then
        call rlistb_half(natom,iatom1,iatom2,x,y,z,itype, pad,
     x            nneigh,markn,jatomn,dxn,dyn,dzn,rn,lmax)
        return
      end if

        k=0
        do 190 i=iatom1,iatom2
        markn(i)=k
        it=itype(i)
        do 170 j=1,natom
        if(j.eq.i) go to 170
             if(j.le.i) go to 170  ! this makes it a half list
        dx=x(i)-x(j)
        if(dx.gt.ax) dx=dx-tax
        if(dx.lt.-ax) dx=dx+tax
        dy=y(i)-y(j)
        if(dy.gt.ay) dy=dy-tay
        if(dy.lt.-ay) dy=dy+tay
        dz=z(i)-z(j)
        if(dz.gt.az) dz=dz-taz
        if(dz.lt.-az) dz=dz+taz
        r2=dx*dx+dy*dy+dz*dz
        if(r2.gt.(rcut(it,itype(j))+pad)**2) go to 170
          k=k+1
          jatomn(k)=j
          dxn(k)=dx
          dyn(k)=dy
          dzn(k)=dz
          rn(k)=sqrt(r2)
  170   continue
        if(k.gt.lmax) stop 'abort - out of space in rlistr_half'
        nneigh(i)=k-markn(i)
  190   continue

      return
      end

      subroutine rlistb_half(natom,iatom1,iatom2,x,y,z,itype, pad,
c                              ! supplied
     x          nneigh,markn,jatomn,dxn,dyn,dzn,rn,lmax)
c                              ! returned
c  computes full neighbor list each atom in range iatom1,iatom2
c  uses 3-D boxing to increase speed
c  user supplies:  natom, x(),y(),z(),itype(),
c                  lmax=maximum length for jatomn,dxn,dyn,dzn,rn arrays
c                  pad = extra distance (padding) to add to rcut
c                        (necessary if using later calls to rlistu)
c  returned:
c                  nneigh(natom) = no. of neighbors for each atom
c                  markn(natom) = pointer array for jatomn,dxn,dyn,dzn,rn
c                  jatomn(),dxn(),dyn(),dzn(),rn() - length<=lmax
c           define: m=markn(i)
c              jatomn(m+j) = j'th neighbor of atom i
c              dxn(m+j) = x(i)-x(j)  (dyn,dzn analagous)
c              rn(m+j) = distance from atom i to j'th neighbor of atom i
c
      implicit real*8(a-h,o-z)
      save
      dimension x(natom),y(natom),z(natom),itype(natom)
      dimension nneigh(natom),markn(natom)
      dimension jatomn(lmax),dxn(lmax),dyn(lmax),dzn(lmax),rn(lmax)
      common/potcom/rcut(5,5),tax,tay,taz,ax,ay,az
      include 'wrksubs_tad2.h'
      parameter(mxijk=10000,mxpont=1000000)
c,mxatom=10000)
      include 'parameters.h'
      common/wrkcm10/iwrk10,jwrk10,
     x            ixbox(mxatom),iybox(mxatom),izbox(mxatom),
     x            nijk(mxijk),ipoint(mxpont),
     x            fill10(lwrk10 - 3*mxatom/2 - mxijk/2 - mxpont/2)
      common/nlists/ilist
c  ilist=0 - rlist calls stay there
c  ilist=1 - rlist calls diverted here (rlistb)
c  ilist=2 - rlist calls diverted here (rlistb) with no boxing

      if(iwrk10.le.0 .and. ilist.ne.2) then
        iwrk10=1132
        ifbox=1
      else
        ifbox=0
        go to 234
      end if

c begin 3-D boxing scheme

c  determine number of boxes in each of x,y,z directions

c   make each box a little larger, as nec., to exactly fill tax length
      call dmnmx(rcut,9,rdum1,rcmax,idum1,idum2)
      nreach=1
      nxbox =  int(tax/(rcmax/float(nreach)))
      nybox =  int(tay/(rcmax/float(nreach)))
      nzbox =  int(taz/(rcmax/float(nreach)))
      if(min(nxbox,nybox,nzbox).lt.1+2*nreach) then
        ifbox=0
        go to 234
      end if
c    xbox = x length of the box  (a little larger than rcmax/nreach)
      xbox = tax/float(nxbox)
      ybox = tay/float(nybox)
      zbox = taz/float(nzbox)
      itest=nxbox*nybox*nzbox
      if(itest.gt.mxijk .or. itest.lt.0) then
          ifbox = 0
          write(6,*)  'note =- not able to box in rlistb_half'
          go to 234
c         stop 'increase mxijk in rlistb_half'
      end if

c  determine which box each atom is in to find max in any one box

      do 10 ijk=1,nxbox*nybox*nzbox
   10 nijk(ijk)=0

c    following scheme good unless an atom coord is more negative than 10*tax
      do 20 i=1,natom
      ixbox(i)=mod(x(i)+100.0d0*tax,tax)/xbox + 1
      iybox(i)=mod(y(i)+100.0d0*tay,tay)/ybox + 1
      izbox(i)=mod(z(i)+100.0d0*taz,taz)/zbox + 1
      if(ixbox(i).lt.1 .or. ixbox(i).gt.nxbox) stop 'fix x boxing'
      if(iybox(i).lt.1 .or. iybox(i).gt.nybox) stop 'fix y boxing'
      if(izbox(i).lt.1 .or. izbox(i).gt.nzbox) stop 'fix z boxing'
      ijk = (izbox(i)-1)*nxbox*nybox + (iybox(i)-1)*nxbox + ixbox(i)
      nijk(ijk)=nijk(ijk)+1
   20 continue
      nbmax = imax(nijk,nxbox*nybox*nzbox)
      if(nbmax*nxbox*nybox*nzbox.gt.mxpont)stop'rlistb_half:  mxpont'

      do 25 ijk=1,nxbox*nybox*nzbox
      nijk(ijk)=0
   25 continue

c  fill ipoint array so atoms in that box can be found

c     ipoint array is ordered as follows:
c      fastest index:  icell      max=nbmax
c                       x         max=nxbox
c                       y         max=nybox
c      slowest index:   z         max=nzbox

      do 30 i=1,natom
      ix = ixbox(i)
      iy = iybox(i)
      iz = izbox(i)
      ijk = (iz-1)*nxbox*nybox + (iy-1)*nxbox + ix
      nijk(ijk)=nijk(ijk)+1
      icell=nijk(ijk)
      ijkm = (ijk-1)*nbmax + icell
      ipoint(ijkm)=i
   30 continue

c now have the following arrays of interest:

c     ipoint(...)  given box, can look up each atom in that box
c     nijk(ijk) = number of atoms in box ijk
c     ixbox(i)  =  x box number for atom i
c     iybox(i)  =  y box number for atom i
c     izbox(i)  =  z box number for atom i

c make neighbor list
c  note that this code is only good if atoms are within roughly one tax of origin

 234  continue
      if(ifbox.ne.0) then
c        write(6,*) 'will perform 3-D boxing for fast neighbor list'
        k=0
        do 90 i=iatom1,iatom2
        markn(i)=k
        it=itype(i)
ccc       do 70 j=1,natom

        do 80 ixraw=ixbox(i)-nreach,ixbox(i)+nreach
        do 70 iyraw=iybox(i)-nreach,iybox(i)+nreach
        do 60 izraw=izbox(i)-nreach,izbox(i)+nreach
        ix =   mod(ixraw + nxbox-1,nxbox) + 1
        iy =   mod(iyraw + nybox-1,nybox) + 1
        iz =   mod(izraw + nzbox-1,nzbox) + 1
        ijk = (iz-1)*nxbox*nybox + (iy-1)*nxbox + ix
        do 50 icell=1,nijk(ijk)
        ijkm = (ijk-1)*nbmax + icell
        j=ipoint(ijkm)
        if(j.eq.i) go to 50
             if(j.le.i) go to 50  ! this makes it a half list
        dx=x(i)-x(j)
        if(dx.gt.ax) then
           dx=dx-tax
        else if(dx.lt.-ax) then
           dx=dx+tax
        end if
        dy=y(i)-y(j)
        if(dy.gt.ay) then
          dy=dy-tay
        else if(dy.lt.-ay) then
          dy=dy+tay
        end if
        dz=z(i)-z(j)
        if(dz.gt.az) then
          dz=dz-taz
        else if(dz.lt.-az) then
          dz=dz+taz
        end if
        r2=dx*dx+dy*dy+dz*dz
        if(r2.gt.(rcut(it,itype(j))+pad)**2) go to 50
          k=k+1
          jatomn(k)=j
          dxn(k)=dx
          dyn(k)=dy
          dzn(k)=dz
          rn(k)=sqrt(r2)
   50   continue
   60   continue
   70   continue
   80   continue
        if(k.gt.lmax) stop 'abort - out of space in rlistb_half'
        nneigh(i)=k-markn(i)
   90   continue

      else if(ifbox.eq.0) then
c        write(6,*) 'NOTE - skipping 3-D boxing: full neighbor list work'
        k=0
        do 190 i=iatom1,iatom2
        markn(i)=k
        it=itype(i)
        do 170 j=1,natom
        if(j.eq.i) go to 170
             if(j.le.i) go to 170  ! this makes it a half list
        dx=x(i)-x(j)
        if(dx.gt.ax) dx=dx-tax
        if(dx.lt.-ax) dx=dx+tax
        dy=y(i)-y(j)
        if(dy.gt.ay) dy=dy-tay
        if(dy.lt.-ay) dy=dy+tay
        dz=z(i)-z(j)
        if(dz.gt.az) dz=dz-taz
        if(dz.lt.-az) dz=dz+taz
        r2=dx*dx+dy*dy+dz*dz
        if(r2.gt.(rcut(it,itype(j))+pad)**2) go to 170
          k=k+1
          jatomn(k)=j
          dxn(k)=dx
          dyn(k)=dy
          dzn(k)=dz
          rn(k)=sqrt(r2)
  170   continue
        if(k.gt.lmax) stop 'abort - out of space in nonbox rlistb_half'
        nneigh(i)=k-markn(i)
  190   continue
      end if

      if(iwrk10.eq.1132) iwrk10=0

      return
      end

c xxx put in place 10/25/02 - seems to work right.
      subroutine g0faster2(natom,nmove,x,y,z,itype, 
     x                                tax,tay,taz,
     x                                maxtyp,rcut,
     x                                energy,grad,hyperratio)
c like g0faster, but for two or more atom types, and interp not inline yet
c (i.e., this does halflists, only when necessary)

c taken from parmdg7.2, put in rcut check in loop 70, put in "save"  afv  10/24/02
c   taken from cgsubs.f, modified for verlet_alone style, and installed into
c   parmdg7 on 12/3/01
c   (added tax,tay,taz,maxtyp,rcut, and hyperratio to call, and
c    took lwrk1 out of parm stmt)
c    

      implicit real*8(a-h,o-z)
      dimension x(natom),y(natom),z(natom),itype(natom)
      dimension grad(*)
ccc      include 'wrksubs_alone.h'  ! this was in parmdg7
      include 'wrksubs_tad2.h' 
      dimension rcut(maxtyp,maxtyp)
ccc      parameter(lwrk1=1200000,natmax=19000,lmax=185000)
c      parameter(natmax=10 000)
      include 'parameters.h'
      parameter(lmax1=2*(lwrk1-11*natmax)/9-1 )
      parameter(lmaxh=lmax1/2)
      parameter(lmax=lmaxh*2)  ! this makes it an even number
      common/wrkcm1/iwrk1,jwrk1,rhoj(natmax),jtype(natmax),fpi(natmax),
     x     xhold(natmax),yhold(natmax),zhold(natmax),itypehold(natmax),
     x                 dxnow(natmax),dynow(natmax),dznow(natmax),
     x                 rhotot(natmax),nneigh(natmax),markn(natmax),
     x         jatomn(lmax),dxn(lmax),dyn(lmax),dzn(lmax),rn(lmax),
     x         fill1(lwrk1 - 22*natmax/2 - 9*lmax/2)
      common/mdcom/dt,tnow,mdwrt,mdfast,mdcool,mddum,padmd
      common/uflags/iuflag(100)
ccccccccc      common/potcom/rcut(3,3),tax,tay,taz,ax,ay,az
        common/holdtr/ a(40000),dum(80000)
        common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                         mark1(99),mark2(99),
     x                         iflo(99),ifhi(99),last
      dimension rcuthold(5,5)
c      save padhold,taxhold,tayhold,tazhold,rcuthold

      save

      rho(rr,it)=trpfn3(rr,0,70+it)
      phi(rr,it,jt)=trpfn3(rr,0,80+it*jt)
      f(xx,it)=trpfn3(xx,0,90+it)
      tgrho(rr,it)=trpfun(rr,11,70+it)
      tgphi(rr,it,jt)=trpfun(rr,1,80+it*jt)/rr
      fp(xx,it)=trpfun(xx,1,90+it)


      if(natom.gt.natmax) stop 'g0faster2 abort: natom.gt.natmax'

cc      pad=float(iuflag(92))/100.0d0
      pad=max(0.0d0,padmd)

c check that cutoff+2*pad is less than half of tax,tay or taz.
c (only 1*pad necessary?)
      rcmax=0.0d0
      do 5 it=1,ntype
      do 5 jt=1,ntype
      rcmax=max(rcmax,rcut(it,jt))
    5 continue
      if(rcmax+2.d0*pad.gt.tax/2.0d0) stop 'g0faster2 - tax too small'
      if(rcmax+2.d0*pad.gt.tay/2.0d0) stop 'g0faster2 - tay too small'
      if(rcmax+2.d0*pad.gt.taz/2.0d0) stop 'g0faster2 - taz too small'

c check maximum atom move since last full update, and anything else that might matter
c  8/08/01 - also checking tax,tay,taz, rcut now  (after Jim Sprague noticed problem)
      inew=0
      dr2max=0.0d0
      mxdeltyp=0
      do 10 i=1,natom
      dxnow(i)=x(i)-xhold(i)
      dynow(i)=y(i)-yhold(i)
      dznow(i)=z(i)-zhold(i)
      dr2=dxnow(i)**2+dynow(i)**2+dznow(i)**2
      dr2max=max(dr2max,dr2)
      mxdeltyp=max(mxdeltyp,abs(itype(i)-itypehold(i)))
   10 continue
      if(sqrt(dr2max).gt.pad/2.0d0 .or. mxdeltyp.ne.0) inew=1
      if(pad.ne.padhold) inew=1
      if(taxhold.ne.tax.or.tayhold.ne.tay.or.tazhold.ne.taz) inew=1
      do i=1,ntype
      do j=1,ntype
        if(rcut(i,j).ne.rcuthold(i,j)) inew=1
      end do
      end do

      if(iwrk1.ne.-7779) inew=1
      call wrkclaimv(1,1,'g0faster2')

c make fresh set of neighbor lists over all atoms if necessary
      if(inew.eq.1) then
        call rlistr_half(natom,1,natom,x,y,z,itype, pad,
     x            nneigh,markn, jatomn,dxn,dyn,dzn,rn, lmax)
        do 12 i=1,natom
        xhold(i)=x(i)
        yhold(i)=y(i)
        zhold(i)=z(i)
        itypehold(i)=itype(i)
        dxnow(i)=0.0d0
        dynow(i)=0.0d0
        dznow(i)=0.0d0
   12   continue

        padhold=pad
        taxhold=tax
        tayhold=tay
        tazhold=taz
        do i=1,ntype
        do j=1,ntype
          rcuthold(i,j)=rcut(i,j)
        end do
        end do
      end if


c  evaulate rhotot and fpi for all atoms and sum up energy
c  works on a half list

      call dmzer(rhotot,natom)
      epair=0.0d0

      do 40 i=1,natom
      rhoi=0.0d0
      it=itype(i)
      do 20 jj=1,nneigh(i)
      mj=markn(i)+jj
      j=jatomn(mj)
      jt=itype(j)
        dxij=dxn(mj)+dxnow(i)-dxnow(j)     !update
        dyij=dyn(mj)+dynow(i)-dynow(j)     !update
        dzij=dzn(mj)+dznow(i)-dznow(j)     !update
        rij=sqrt(dxij**2+dyij**2+dzij**2)     !update
        rn(mj)=rij     !update

        rho3j=rho(rij,jt)
        if(it.eq.jt) then
          rho3i=rho3j
        else
          rho3i=rho(rij,it)
        end if
        rhotot(i)=rhotot(i) + rho3j
        rhotot(j)=rhotot(j) + rho3i

        phi3=phi(rij,it,jt)
        epair = epair + phi3

   20 continue
   40 continue

      embed=0.0d0
      do 60 i=1,natom
      it=itype(i)
      embed = embed + f(rhotot(i),it)
      fpi(i)=fp(rhotot(i),it)
   60 continue

ccc      energy = embed + 0.5d0*epair    ! full list
      energy = embed + epair       ! half-neigh-list code


c calculate gradient over all moving atoms

      call dmzer(grad,3*natom)
      do 80 i=1,nmove
      it=itype(i)
      ix=3*(i-1)+1
      iy=ix+1
      iz=ix+2
      do 70 jj=1,nneigh(i)
      mj=markn(i)+jj
      j=jatomn(mj)
      jt=itype(j)
      rij=rn(mj)    ! (rn(mj) filled above to save the sqrt)
      if(rij.gt.rcut(it,jt)) go to 70
        jx=3*(j-1)+1   ! half-neigh-list code
        jy=jx+1        ! half-neigh-list code
        jz=jx+2        ! half-neigh-list code

        dxij=dxn(mj)+dxnow(i)-dxnow(j)     !update
        dyij=dyn(mj)+dynow(i)-dynow(j)     !update
        dzij=dzn(mj)+dznow(i)-dznow(j)     !update

        tgphi3=tgphi(rij,it,jt)
        tgrho3i=tgrho(rij,it)
        tgrho3j=tgrho(rij,jt)

        tgrad = tgphi3 + fpi(i)*tgrho3j + fpi(j)*tgrho3i
        grad(ix)=grad(ix)+tgrad*dxij
        grad(iy)=grad(iy)+tgrad*dyij
        grad(iz)=grad(iz)+tgrad*dzij
        grad(jx)=grad(jx)-tgrad*dxij  ! half-neigh-list code
        grad(jy)=grad(jy)-tgrad*dyij  ! half-neigh-list code
        grad(jz)=grad(jz)-tgrad*dzij  ! half-neigh-list code
   70 continue
   80 continue

      call wrkrelease(1)
      iwrk1=-7779
      return
      end
