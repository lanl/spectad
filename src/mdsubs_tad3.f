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



c mdsubs_tad2.f     12/6/01

c This is a set of MD-related routines (and various misc. routines)
c for the tad2 program
c still a bit rough
c This is a subset of the "verlet_alone.f" routine package in parmdg7.
c The rest of those routines are in "e0subs_tad2.f

c Langevin stuff was put in in a hurry.
c For safety, just use one or zero thermostats (ntherm), and make sure
c ctherm(1)=0.0 if ntherm=0


      subroutine gcalc(natom,nmove,xyz,itype,taxes,ietype,maxtyp,
     x                                  rcut,energy,grad, hyperratio)
c converted to single-array format  12/20/01
c gradient routine - compute gradient for first nmove atoms (or perhaps all natom)
c call is redirected to proper potential-specific routine
      use mod_lammps !crz
      
      implicit real*8(a-h,o-z)
      save
c      parameter(mxatom=10000)
      include 'parameters.h'
      dimension x(mxatom),y(mxatom),z(mxatom),itype(natom)
      dimension xyz(3,natom)
      dimension taxes(3)
      dimension rcut(maxtyp,maxtyp)
      dimension grad(3,natom)
      common/callcom/ngradcall
      
      
      call map1to3(natom,x,y,z,xyz)
      call map1to3(1,tax,tay,taz,taxes)

      if(ietype.eq.0) then
        natom_tmp = natom !- 3*98
        nmove_tmp = nmove !- 3*98
        call g0(natom_tmp,nmove_tmp,x,y,z,itype,
     x                          tax,tay,taz,
     x                          ietype, maxtyp,rcut,
     x                          energy,grad, hyperratio)
crz -- Lammps Force Call
      else if((ietype.ge.100).and.(ietype.lt.200)) then
        !call MPI_BARRIER(force_comm,ier)
        !if(rank_f.eq.0)write(217+spawnID,*)'spawnID ',spawnID,'IN g100'
        !if(rank_f.eq.0)call flushtad(217+spawnID)
         if(spawn_lammps)then
           call g100spawn(natom,nmove,xyz,ietype,itype,energy,grad)
           call map1to3(natom,x,y,z,xyz)
         else
           call g100(natom,nmove,x,y,z,itype,ietype,maxtyp,energy,grad
     x          ,taxes)
         endif
        !call MPI_BARRIER(force_comm,ier)
        !if(rank_f.eq.0)write(217+spawnID,*)'spawnID ',spawnID,'OUT g100'
        !if(rank_f.eq.0)call flushtad(217+spawnID)
      else if (ietype.eq.27) then
         call g27(natom,nmove,x,y,z,itype, 
     x                          tax,tay,taz,
     x                          ietype, maxtyp,rcut,
     x                          energy,grad, hyperratio)
      else if (ietype.eq.67) then
         call g67(natom,nmove,x,y,z,itype,tax,tay,taz,ietype,maxtyp,rcut
     +       ,energy,grad,hyperratio)
      else if (ietype.eq.68) then
         call g68(natom,nmove,x,y,z,itype,tax,tay,taz,ietype,maxtyp,rcut
     +       ,energy,grad,hyperratio)
      else if (ietype.eq.69) then
         call g69(natom,nmove,x,y,z,itype,tax,tay,taz,ietype,maxtyp,rcut
     +       ,energy,grad,hyperratio)
      else if (ietype.eq.70) then
         call g70(natom,nmove,x,y,z,itype,tax,tay,taz,ietype,maxtyp,rcut
     +       ,energy,grad,hyperratio)
      else if (ietype.eq.72) then
         call  g72(natom,nmove,x,y,z,itype,maxtyp,tax,tay,taz,rcut
     +       ,energy,grad)
      else if (ietype.eq.73) then
         call  g73(natom,nmove,x,y,z,itype,maxtyp,tax,tay,taz,rcut
     +       ,energy,grad)
      else if (ietype.eq.74) then
         call g74(natom,nmove,x,y,z,itype,tax,tay,taz,ietype,maxtyp,rcut
     +       ,energy,grad)
      else if (ietype.eq.75) then
         call g75(natom,nmove,x,y,z,itype,maxtyp,tax,tay,taz,rcut
     +       ,energy,grad)
      else if (ietype.eq.76) then
         call g76(natom,nmove,x,y,z,tax,tay,taz,energy,grad)
      else if (ietype.eq.78) then
         call g78(natom,nmove,x,y,z,itype,tax,tay,taz,ietype,maxtyp,rcut
     +       ,energy,grad)
      else if (ietype.eq.80) then
         call  g80(natom,nmove,x,y,z,itype,maxtyp,tax,tay,taz,rcut
     +       ,energy,grad)
      else if (ietype.eq.999) then
         call g999(natom,nmove,x,y,z,itype,maxtyp,tax,tay,taz,rcut
     +       ,energy,grad)
      else
        write(6,*) 'ietype',ietype
        stop 'gcalc: unrecognized ietype'
      end if

      call map3to1(natom,x,y,z,xyz)
      call map3to1(1,tax,tay,taz,taxes)

      ngradcall=ngradcall+1

c       write(6,*) ngradcall

      ! Do blocking here...
      if(allow_blocking.and.blocking_on)then
      if(statedata(istatenow)%nneigh_block.gt.0)then
        do ineigh = 1, 
     +         INT(SIZE(statedata(istatenow)%ibondneigh1))
          i1 = statedata(istatenow)%ibondneigh1(ineigh)
          i2 = statedata(istatenow)%ibondneigh2(ineigh)
          i3 = statedata(istatenow)%ibondneigh3(ineigh)
          dn12 = statedata(istatenow)%dbondneigh12(ineigh)
          dn13 = statedata(istatenow)%dbondneigh13(ineigh)
          dn12i = statedata(istatenow)%dbondneigh12i(ineigh)
          dn13i = statedata(istatenow)%dbondneigh13i(ineigh)
          if((i1.gt.0).and.(i2.gt.0).and.(i3.gt.0))then
              ! Ratio of "saddle" stretch to start adding blocking force:
              stol = 0.8
              ! Get simple bond distance 1-2:
              dx12 = (x(i1) - x(i2))
              dy12 = (y(i1) - y(i2))
              dz12 = (z(i1) - z(i2))
              dist12 = sqrt(dx12**2+dy12**2+dz12**2)
              st12 = dist12 / dn12i     ! <- Current stretch of 1-2 bond
              st12sd = dn12 / dn12i     ! <- Stretch of 1-2 bond at saddle pt
              st12mx = stol * st12sd    ! <- Maximum (unblocked) stretch to allow for 1-2 bond
              rat12  = st12 / (st12sd * 0.99)
              ! Get simple bond distance 1-3:
              dx13 = (x(i1) - x(i3))
              dy13 = (y(i1) - y(i3))
              dz13 = (z(i1) - z(i3))
              dist13 = sqrt(dx13**2+dy13**2+dz13**2)
              st13 = dist13 / dn13i     ! <- Current stretch of 1-3 bond
              st13sd = dn13 / dn13i     ! <- Stretch of 1-3 bond at saddle pt
              st13mx = stol * st13sd    ! <- Maximum (unblocked) stretch to allow for 1-3 bond
              rat13  = st13 / (st13sd * 0.99)
              if((st12.ge.st12mx).and.(st13.ge.st13mx))then
!                  gradmag12 = (st12 - st12mx) * 1.0d0
!                  gradmag13 = (st13 - st13mx) * 1.0d0
!                  dgrad1 =   (dx12/dist12) * gradmag12
!     +                     + (dx13/dist13) * gradmag13
!                  dgrad2 =   (dy12/dist12) * gradmag12
!     +                     + (dy13/dist13) * gradmag13
!                  dgrad3 =   (dz12/dist12) * gradmag12
!     +                     + (dz13/dist13) * gradmag13
!                  grad(1,i1) = dgrad1 + grad(1,i1)
!                  grad(2,i1) = dgrad2 + grad(2,i1)
!                  grad(3,i1) = dgrad3 + grad(3,i1)
                  ! Get normals in 'escape' direction
                  gnxe12 = (dx12/dist12)
                  gnye12 = (dy12/dist12)
                  gnze12 = (dz12/dist12)
                  gnxe13 = (dx13/dist13)
                  gnye13 = (dy13/dist13)
                  gnze13 = (dz13/dist13)
                  ! Get gradient in escape directions
                  gesc12 = gnxe12*grad(1,i1)+
     +                     gnye12*grad(2,i1)+gnze12*grad(3,i1)
                  gesc13 = gnxe13*grad(1,i1)+
     +                     gnye13*grad(2,i1)+gnze13*grad(3,i1)
                  ! Subtract off escape gradient
                  if ((gesc12).gt.0.0) gesc12 = 0.0
                  if ((gesc13).gt.0.0) gesc13 = 0.0
                  grad(1,i1) = grad(1,i1) 
     +                       - gesc12*gnxe12*(1.0*rat12)
     +                       - gesc13*gnxe13*(1.0*rat13)
                  grad(2,i1) = grad(2,i1)
     +                       - gesc12*pnye12*(1.0*rat12)
     +                       - gesc13*pnye13*(1.0*rat13)
                  grad(3,i1) = grad(3,i1)
     +                       - gesc12*pnze12*(1.0*rat12)
     +                       - gesc13*pnze13*(1.0*rat13)
                  !if((rat12.gt.(1.0)).or.(rat13.gt.(1.0)))then
                  if(.false.)then
                    ! Write to *.lis:
                    write(6,*)'SpawnID',SpawnID,
     +               'BLOCKING ineigh',ineigh,
     +               '- rat12',rat12,'gesc12',gesc12,
     +               '- rat13',rat13,'gesc13',gesc13
                  endif
              endif
          endif
        enddo
      endif
      endif


      return
      end

      subroutine getclsnew(natom,ntitle,ctitle,xyz,pxyz,itype,
     x                           taxes,iunit,iwrt,ierr)
c 12/20/01 - changed over to single-array format (xyz,pxyz,taxes)
c getclsnew is just like getcls, but uses a character title string instead of real*8 afv 12/15/01
c ctitle must be dimensioned character*80 and as an array at least 20 long
      use mod_mpi
      implicit real*8(a-h,o-z)
      !include 'mpif.h'
      dimension xyz(3,*),itype(*),pxyz(3,*),taxes(3)
      character*80 ctitle(*)
      character*30 form

c  read in cluster file

      ierr=0
      if(rank_w.eq.0) read(unit=iunit,fmt=30,end=766,err=666)
     + natom,iop,ntitle
      call MPI_BCAST(natom,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(iop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(ntitle,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      
      if(ntitle.gt.20) stop 'ntitle too large in getclsnew'
      
      if(rank_w.eq.0) write(*,*) 'naton,iop,ntitle = ',natom,iop,ntitle

c  original standard cluster file up to 7/91
      if(iop.eq.0) then
       if(rank_w.eq.0) read(unit=iunit,fmt=41,end=766,err=666)
     x   (ctitle(i),i=1,ntitle)
       if(rank_w.eq.0) read(unit=iunit,fmt=50,end=766,err=666)
     x   (taxes(j),j=1,3)
       if(rank_w.eq.0) read(unit=iunit,fmt=50,end=766,err=666)
     x      ((xyz(j,i),j=1,3),itype(i),
     x       (pxyz(j,i),j=1,3),i=1,natom)

      else if(iop.eq.1) then
c  x,y,z-only cluster file with format
        if(rank_w.eq.0) read(unit=iunit,fmt=41,end=766,err=666)
     x   (ctitle(i),i=1,ntitle)
        if(rank_w.eq.0) read(unit=iunit,fmt=50,end=766,err=666)
     x   (taxes(j),j=1,3)
        if(rank_w.eq.0) read(unit=iunit,fmt=55,end=766,err=666) form
        if(rank_w.eq.0) read(unit=iunit,fmt=form,end=766,err=666)
     x      ((xyz(j,i),j=1,3),itype(i),i=1,natom)

      else if(iop.eq.2) then
c  x,y,z,px,py,pz cluster file with format
        if(rank_w.eq.0) read(unit=iunit,fmt=41,end=766,err=666)
     x   (ctitle(i),i=1,ntitle)
        if(rank_w.eq.0) read(unit=iunit,fmt=50,end=766,err=666)
     x   (taxes(j),j=1,3)
        if(rank_w.eq.0) read(unit=iunit,fmt=55,end=766,err=666) form
         call MPI_BCAST(form,30,MPI_CHARACTER,0,
     x   MPI_COMM_WORLD,ier)
        if(rank_w.eq.0) read(unit=iunit,fmt=form,end=766,err=666)
     x      ((xyz(j,i),j=1,3),itype(i),
     x       (pxyz(j,i),j=1,3),i=1,natom)

      else if(iop.eq.3) then
c  simplest possible cluster file format
c  card 1: just put a 3 in column 7  (natom,ntitle ignored)
c  card 2: title card (a80)
c  card 3: tax,tay,taz  (free format)
c  cards 4-n: x,y,z  (free format)

c  eof signals end of atoms
        ntitle=1
        if(rank_w.eq.0) read(iunit,41,end=21) (ctitle(i),i=1,ntitle)
        if(rank_w.eq.0) read(iunit,*,end=21,err=21) (taxes(j),j=1,3)
        do 20 i=1,10000
        if(rank_w.eq.0) read(iunit,*,end=21,err=21) (xyz(j,i),j=1,3)
        itype(i)=1
   20   continue
         stop 'getcls: iop3 input file too long'
   21   natom=i-1

      else if(iop.eq.4) then
c     Import xcls96 output file
c     card2: natom
c     card3: taxes
c     card4-natom+4: x,y,z,type
         if(rank_w.eq.0) read(iunit,*,end=766,err=666) natom
         call MPI_BCAST(natom,1,MPI_INTEGER,0,
     x   MPI_COMM_WORLD,ier)
         if(rank_w.eq.0) read(iunit,50,end=766,err=666) (taxes(i),i=1,3)
         do 22 i = 1,natom
          if(rank_w.eq.0) read(iunit,50,end=766,err=666)
     x                           (xyz(j,i),j=1,3),itype(i)
 22      continue
      else if(iop.eq.5) then
c     Import xcls96 output file
c     card2: natom
c     card3: taxes
c     card4-natom+4: x,y,z,type,px,py,pz
         if(rank_w.eq.0) read(iunit,*,end=766,err=666) natom
         call MPI_BCAST(natom,1,MPI_INTEGER,0,
     x   MPI_COMM_WORLD,ier)
         if(rank_w.eq.0) read(iunit,50,end=766,err=666) (taxes(j),j=1,3)
         do 23 i = 1,natom
          if(rank_w.eq.0) read(iunit,50,end=766,err=666)
     x           (xyz(j,i),j=1,3),itype(i),(pxyz(j,i),j=1,3)
 23      continue

      else
        if(rank_w.eq.0) write(6,*) 'iop=',iop,' in getcls'
        stop 'getcls iop not understood'
      end if
   30 format(i5,i2,i3)
   41 format(a80)
   50 format(3f20.8,i5,3f20.8)
   55 format(a)
   
crz -- Bcast Info From File...
      call MPI_BCAST(taxes,3,MPI_REAL8,0,
     x   MPI_COMM_WORLD,ier)
      call MPI_BCAST(itype,natom,MPI_INTEGER,0,
     x   MPI_COMM_WORLD,ier)
      call MPI_BCAST(xyz,3*natom,MPI_REAL8,0,
     x   MPI_COMM_WORLD,ier)
      if(iop.lt.4)then
       call MPI_BCAST(ctitle,80*ntitle,MPI_CHARACTER,0,
     x   MPI_COMM_WORLD,ier)
      endif
      if((iop.ne.1).and.(iop.ne.3).and.(iop.ne.4))then
       call MPI_BCAST(pxyz,3*natom,MPI_REAL8,0,
     x   MPI_COMM_WORLD,ier)
      endif

      if (iwrt.ge.1) then
        ntype=imax(itype,natom)
        if(rank_w.eq.0)then
         write(6,*)
         write(6,*)'>>>> formatted cluster file read in <<<<'
         write(6,*)'________________________________________'
         write(6,77) natom,ntype,(taxes(j),j=1,3)
        endif
   77   format(' natom=',i6,'  ntype=',i2,
     x           '  tax=',f16.8,'  tay=',f16.8,'  taz=',f16.8)
c        write(6,70) natom
c        write(6,71) ntype
c        write(6,72) (taxes(j),j=1,3)
c   70   format(' number of atoms ',i4)
c   71   format(' number of types ',i2)
c   72   format(' tax=',f18.8,'  tay=',f18.8,'  taz=',f18.8)
        if (ntitle.gt.0) then
          if(rank_w.eq.0) write(6,73) ntitle
          if(rank_w.eq.0) write(6,74) (ctitle(i),i=1,ntitle)
   73     format(1x,i2,' title cards:')
   74     format(a80)
        endif
        !if(iop.ne.0) write(6,90) form
   90   format(' coords used format ',a)
      end if
      return

  666 continue
      ierr=666
      return

  766 ierr=766
      return
      end

      subroutine mdstep(natom,nmove,xyz,itype,pxyz,taxes,ietype, 
     x                  maxtyp,rcut,amu,
     x maxtherm,itherm,jtherm,ttherm,ctherm,xtherm,mtherm,ktherm,ntherm,
     x           moving,mx,my,mz, dt,tnow,thyper,hyperratio, e,grad)

c standalone routine that takes 1 Verlet step
c uses new langevin integration scheme from A & T 

      use mod_mpi
      implicit real*8(a-h,o-z)
c      parameter(mxatom=10000)
      include 'parameters.h'
      dimension xyz(3,natom)
      dimension pxyz(3,natom)
      dimension taxes(3)
      dimension x(mxatom),y(mxatom),z(mxatom),itype(natom)
      dimension px(mxatom),py(mxatom),pz(mxatom)
      dimension rcut(maxtyp,maxtyp),amu(maxtyp)
      dimension moving(natom),mx(natom),my(natom),mz(natom)
      dimension grad(3,natom)

      dimension itherm(maxtherm),jtherm(maxtherm)
      dimension ttherm(maxtherm),ctherm(maxtherm)
      dimension xtherm(maxtherm),mtherm(maxtherm),ktherm(maxtherm)

      real*8 dx,dy,dz,one

      include 'wrksubs_tad2.h'
c      parameter (movmax=10000)
c      include 'parameters.h'
      parameter (m3=movmax*3)
c      parameter (natmax=10000)
      common/dercom/vd,istatd
      common/wrkcm2/iwrk2,jwrk2,s(6*natmax),dsdt(6*natmax),
     x              fill2(lwrk2-2*6*natmax)
      common/browncom/sigr(5),sigv(5),crv(5),ibrown ! BPU type here
      common/print/ip1,ip2,ip3,ip4,ip5,ip6
c      parameter(nhold=10000) ! now in include
      dimension dsdtold(6*nholdmdsubs)
      dimension rmass(5) ! BPU type here
      save
c      save dsdtold,ntypebrown,rmass,sigr,sigv,crv
       
      call map1to3(natom,x,y,z,xyz)
      call map1to3(natom,px,py,pz,pxyz)
      call map1to3(1,tax,tay,taz,taxes)

c  ntherm = number of thermostats
c  for the i'th thermostat:
c    itherm(i)= start of atom range for this thermostat
c    jtherm(i)= end of atom range(-1=natom)
c    ttherm(i)= desired temperature(K)
c    ctherm(i)= response rate(/sec)
c    mtherm(i)= atom no. defining moving ref. frame
c    ktherm(i)= thermo type(0=off,1=Nose,2=Beren,3=frict,4=Langevin,5=Gauss)

      call wrkclaimv(2,1,'wrkcm2 conflict in md45step')

c kludge warning:
c        assuming 1 thermostat
      if(natom.gt.nholdmdsubs) stop 'redim dsdtold in mdstep'

      ntype=imax(itype,natom)
      if((ibrown.ne.1 .and. ibrown.ne.2) .or. ntype.gt.ntypebrown) then
        ibrown=1
        ntypebrown=ntype
        beta=ctherm(1)
        if(ntherm.eq.0 .or. ktherm(1).eq.0) beta=0.0d0
        boltz=3.167d-6
        tk=ttherm(1)*boltz
        do 20 it=1,ntype
        rmass(it)=1822.83d0*amu(it)
   20   continue

c set up coeficients for Brownian integration (Allen and Tildesley, p. 261)
        bdt=beta*dt
        if(beta.gt.0.0d0) then
          c0=exp(-bdt)         ! roughly 1 - bdt
          c1=(1.0d0-c0)/bdt    ! roughly 1 - bdt/2
          c2=(1.0d0-c1)/bdt    ! roughly 1/2 - bdt/6
        else
          c0=1.0d0
          c1=1.0d0
          c2=0.5d0
        end if

        if(ip1.ge.1) then
          write(6,*) 'Allen and Tildesley Brownian Verlet integration'
          write(6,*) 'bdt,c0,c1,c2:',bdt,c0,c1,c2
        end if

c set up variances and correlation coef. for x and v displacements (A&T, p. 262)
        bdt=beta*dt
        if(beta.gt.0.0d0) then
c calculate sigr and sigv in a.u.
          betaau=beta*2.418d-17
          do 30 it=1,ntype
          sig2r=tk/rmass(it)/betaau/betaau *
     x        (2.d0*bdt-3.d0+4.d0*exp(-bdt)-exp(-2.d0*bdt))
          sig2v=tk/rmass(it)*(1.0d0-exp(-2.d0*bdt))
          sigrsigvcrv=tk/rmass(it)/betaau*(1.0d0-exp(-bdt))**2
          sigr(it)=sqrt(sig2r)
          sigv(it)=sqrt(sig2v)
          crv(it)=sigrsigvcrv/(sigr(it)*sigv(it))
          if(ip1.ge.1) then
            write(6,*) it,'sigr,sigv,crv',sigr(it),sigv(it),crv(it)
          end if
   30     continue
        else
          do 35 it=1,ntype
          sigr(it)=0.0d0
          sigv(it)=0.0d0
          crv(it)=0.0d0
          if(ip1.ge.1) then
            write(6,*) it,'sigr,sigv,crv',sigr(it),sigv(it),crv(it)
          end if
   35     continue
        end if
      end if

      nthermhold=ntherm
      kthermhold=ktherm(1)
      ntherm=0
      ktherm(1)=0

      cfac=0.529d0/2.418d-17
c.....................
c  one time:
c      call vp(x,y,a,c,d,dvdx,dvdy)
c      accelx = -(dvdx)/rmass
c........
cc  Brownian integration step using Allen and Tildesley method, p 263.
c      x = x + c1*dt*velx + c2*dt**2*accelx + randx
c      call vp(x,y,a,c,d,dvdx,dvdy)
c      accelnewx = -(dvdx)/rmass
c      velx = c0*velx + (c1-c2)*dt*accelx + c2*dt*accelnewx + randvx
c      accelx=accelnewx   ! ready for next time
c...............................

      !call MPI_BARRIER(force_comm,ier)
      !if(rank_f.eq.0) write(217+spawnID,*)'spawnID ',spawnID,' IN bb..'
      !if(rank_f.eq.0) call flushtad(217+spawnID)
      im=0
      do 140 i=1,natom
      if(moving(i).eq.0) go to 140
      it=itype(i)
      im=im+1
      rminv=1.0d0/rmass(it)
c      x = x + c1*dt*velx + c2*dt**2*accelx + randx
      if(mx(i).ne.0) then
        dpdtoldx=0.0d0
        if(ibrown.eq.2) dpdtoldx=dsdtold(i+3*natom)
        call bivariate(sigr(it),sigv(it),crv(it),randx,randvx)
        x(i)=x(i)+c1*px(i)*rminv*cfac*dt +
     x    c2*dt**2*dpdtoldx*rminv*0.529d0/(2.418d-17) + randx*0.529d0
        px(i)=px(i)+randvx*rmass(it)
      end if
      if(my(i).ne.0) then
        dpdtoldy=0.0d0
        if(ibrown.eq.2) dpdtoldy=dsdtold(i+4*natom)
        call bivariate(sigr(it),sigv(it),crv(it),randy,randvy)
        y(i)=y(i)+c1*py(i)*rminv*cfac*dt +
     x    c2*dt**2*dpdtoldy*rminv*0.529d0/(2.418d-17) + randy*0.529d0
        py(i)=py(i)+randvy*rmass(it)
      end if
      if(mz(i).ne.0) then
        dpdtoldz=0.0d0
        if(ibrown.eq.2) dpdtoldz=dsdtold(i+5*natom)
        call bivariate(sigr(it),sigv(it),crv(it),randz,randvz)
        z(i)=z(i)+c1*pz(i)*rminv*cfac*dt +
     x    c2*dt**2*dpdtoldz*rminv*0.529d0/(2.418d-17) + randz*0.529d0
        pz(i)=pz(i)+randvz*rmass(it)
      end if
  140 continue
      ! Below BCASTs are now done after 's' is made
      !call MPI_BCAST(x,natom,MPI_REAL8,0,force_comm,ier)
      !call MPI_BCAST(y,natom,MPI_REAL8,0,force_comm,ier)
      !call MPI_BCAST(z,natom,MPI_REAL8,0,force_comm,ier)
      !call MPI_BCAST(px,natom,MPI_REAL8,0,force_comm,ier)
      !call MPI_BCAST(py,natom,MPI_REAL8,0,force_comm,ier)
      !call MPI_BCAST(pz,natom,MPI_REAL8,0,force_comm,ier)

      !call MPI_BARRIER(force_comm,ier)
      !if(rank_f.eq.0) write(217+spawnID,*)'spawnID ',spawnID,' OUT bb..'
      !if(rank_f.eq.0) call flushtad(217+spawnID)


      ! Do blocking here...
      if(allow_blocking.and.blocking_on)then
      if(statedata(istatenow)%nneigh_block.gt.0)then
        do ineigh = 1, 
     +         INT(SIZE(statedata(istatenow)%ibondneigh1))
          i1 = statedata(istatenow)%ibondneigh1(ineigh)
          i2 = statedata(istatenow)%ibondneigh2(ineigh)
          i3 = statedata(istatenow)%ibondneigh3(ineigh)
          dn12 = statedata(istatenow)%dbondneigh12(ineigh)
          dn13 = statedata(istatenow)%dbondneigh13(ineigh)
          dn12i = statedata(istatenow)%dbondneigh12i(ineigh)
          dn13i = statedata(istatenow)%dbondneigh13i(ineigh)
          if((i1.gt.0).and.(i2.gt.0).and.(i3.gt.0))then

              ! Ratio of "saddle" stretch to start adding blocking force:
              stol = 0.8

              ! Get simple bond distance 1-2:
              dx12 = (x(i1) - x(i2))
              dy12 = (y(i1) - y(i2))
              dz12 = (z(i1) - z(i2))
              dist12 = sqrt(dx12**2+dy12**2+dz12**2)
              st12 = dist12 / dn12i     ! <- Current stretch of 1-2 bond
              st12sd = dn12 / dn12i     ! <- Stretch of 1-2 bond at saddle pt
              st12mx = stol * st12sd    ! <- Maximum (unblocked) stretch to allow for 1-2 bond
              rat12  = st12 / (st12sd * (0.99))

              ! Get simple bond distance 1-3:
              dx13 = (x(i1) - x(i3))
              dy13 = (y(i1) - y(i3))
              dz13 = (z(i1) - z(i3))
              dist13 = sqrt(dx13**2+dy13**2+dz13**2)
              st13 = dist13 / dn13i     ! <- Current stretch of 1-3 bond
              st13sd = dn13 / dn13i     ! <- Stretch of 1-3 bond at saddle pt
              st13mx = stol * st13sd    ! <- Maximum (unblocked) stretch to allow for 1-3 bond
              rat13  = st13 / (st13sd * (0.99))

              if((st12.ge.st12mx).and.(st13.ge.st13mx))then
                  ! Get normals in 'escape' direction
                  pnxe12 = (dx12/dist12)
                  pnye12 = (dy12/dist12)
                  pnze12 = (dz12/dist12)
                  pnxe13 = (dx13/dist13)
                  pnye13 = (dy13/dist13)
                  pnze13 = (dz13/dist13)
                  ! Get escape momentum
                  pesc12 = pnxe12*px(i1)+pnye12*py(i1)+pnze12*pz(i1)
                  pesc13 = pnxe13*px(i1)+pnye13*py(i1)+pnze13*pz(i1)
                  ! Subtract off escape momentum
                  if ((pesc12).lt.0.0) pesc12 = 0.0
                  if ((pesc13).lt.0.0) pesc13 = 0.0
                  px(i1) = px(i1) - pesc12*pnxe12*(rat12*1.0)
     +                            - pesc13*pnxe13*(rat13*1.0)
                  py(i1) = py(i1) - pesc12*pnye12*(rat12*1.0)
     +                            - pesc13*pnye13*(rat13*1.0)
                  pz(i1) = pz(i1) - pesc12*pnze12*(rat12*1.0)
     +                            - pesc13*pnze13*(rat13*1.0)
                  !if((rat12.gt.(1.0)).or.(rat13.gt.(1.0)))then
                  if(.false.)then
                    ! Write to *.lis:
                    write(6,*)'SpawnID',SpawnID,
     +               'BLOCKING ineigh',ineigh,
     +               '- rat12',rat12,'- pesc12',pesc12,
     +               '- rat13',rat13,'- pesc13',pesc13
                  endif
              endif
          endif
        enddo
      endif
      endif


c  put x,y,z,px,py,pz into s
      do 10 i=1,natom
      s(i)=x(i)
      s(i+natom)=y(i)
      s(i+2*natom)=z(i)
      s(i+3*natom)=px(i)
      s(i+4*natom)=py(i)
      s(i+5*natom)=pz(i)
   10 continue
      neq=6*natom

      ! Match values to rank_f=0
      !call MPI_BCAST(s,6*natmax,MPI_REAL8,0,force_comm,ier)
      call MPI_BCAST(s,natom*6,MPI_REAL8,0,force_comm,ier)
      do i=1,natom
        if(rank_f.gt.0)then
          x(i)=s(i)
          y(i)=s(i+natom)
          z(i)=s(i+2*natom)
          px(i)=s(i+3*natom)
          py(i)=s(i+4*natom)
          pz(i)=s(i+5*natom)
        endif
      enddo


c  put the Nose-Hoover thermostatting variables into s
c  (note: this is done even if all the thermostats are of type 2)
      do 15 i=1,ntherm
      s(neq+i)=xtherm(i)
   15 continue
      neq=neq+ntherm

c compute derivs
      !call MPI_BARRIER(force_comm,ier)
      !if(rank_f.eq.0) write(217+spawnID,*)'spawnID ',spawnID,' IN D..'
      !if(rank_f.eq.0) call flushtad(217+spawnID)
      call derivs(s,dsdt,neq,natom,nmove,itype,moving,mx,my,mz,
     x                          tax,tay,taz,
     x                          ietype, maxtyp,rcut,amu,
     x maxtherm,itherm,jtherm,ttherm,ctherm,xtherm,mtherm,ktherm,ntherm,
     x                               dt,e,grad,hyperratio)
      !call MPI_BARRIER(force_comm,ier)
      !if(rank_f.eq.0) write(217+spawnID,*)'spawnID ',spawnID,' OUT D..'
      !if(rank_f.eq.0) call flushtad(217+spawnID)
      istat=istatd
      if(istat.ne.0) then
        call map3to1(natom,x,y,z,xyz)
        call map3to1(natom,px,py,pz,pxyz)
        call map3to1(1,tax,tay,taz,taxes)
        ntherm=nthermhold
        ktherm(1)=kthermhold
        iwrk2=0
        return
      end if
      v=vd

c ---------Brownian VERLET leapfrog  --------

c  move p forward 
      do 40 i=natom*3+1,neq
cc      velx = c0*velx + (c1-c2)*dt*accelx + c2*dt*accelnewx + randvx
      if(ibrown.eq.1) then
        accelx=dsdt(i)     ! first time through
      else if(ibrown.eq.2) then
        accelx=dsdtold(i)  ! normal case
      end if
      s(i) = c0*s(i) + (c1-c2)*dt*accelx + c2*dt*dsdt(i)
      dsdtold(i)=dsdt(i)
   40 continue
c -------------------------------------
      ibrown=2

c  put s back into x,y,z,px,py,pz
      do 60 i=1,natom
      x(i)=s(i)
      y(i)=s(i+natom)
      z(i)=s(i+2*natom)
      px(i)=s(i+3*natom)
      py(i)=s(i+4*natom)
      pz(i)=s(i+5*natom)
   60 continue

c  put the Nose-Hoover thermostatting variables back
      do 70 i=1,ntherm
      xtherm(i)=s(6*natom+i)
   70 continue

      ntherm=nthermhold
      ktherm(1)=kthermhold
      tnow=tnow+dt
      thyper=thyper+dt*hyperratio
      istat=istatd
      if(istat.ne.0) write(6,*) 'warning - coords have changed'
      iwrk2=0

      call map3to1(natom,x,y,z,xyz)
      call map3to1(natom,px,py,pz,pxyz)
      call map3to1(1,tax,tay,taz,taxes)

      return
      end


      subroutine mdreset(iop)
c reset Langevin integrator so that it doesn't try to use previously save derivs
c   iop=0 basic reset
c   iop=1 complete reset - recomputes integration parameters at next mdstep call
      implicit real*8(a-h,o-z)
      common/browncom/sigr(5),sigv(5),crv(5),ibrown ! BPU type here

      if(iop.eq.1) then
        ibrown=0
      else if(iop.eq.0) then
        if(ibrown.ge.1) ibrown=1
      else
        stop 'unrecognized iop in mdreset'
      end if

      return
      end


      subroutine derivs(s,dsdt,n, natom,nmove,itype,moving,mx,my,mz,
     x                          tax,tay,taz,
     x                          ietype, maxtyp,rcut,amu,
     x maxtherm,itherm,jtherm,ttherm,ctherm,xtherm,mtherm,ktherm,ntherm,
     x                                   dt,e,grad,hyperratio)
c  compute derivatives in form suitable for Runge Kutta driver
c  note that the coordinates (x,y,z,p,px,py,pz) are in s()
      use mod_mpi
      implicit real*8(a-h,o-z)
      include 'wrksubs_tad2.h'
c      parameter (movmax=10000)
      include 'parameters.h'
      parameter (m3=movmax*3)
c      parameter (natmax=10000)
      parameter (n3=natmax*3)
      dimension rcut(maxtyp,maxtyp),amu(maxtyp)
      common/wrkcm8/iwrk8,jwrk8,work(4*natmax),
     x                   fill8(lwrk8-4*natmax)
      common/wrkcm9/iwrk9,jwrk9,flast(n3),xlast(n3),
     x                   fill9(lwrk9-2*n3)
      dimension itherm(maxtherm),jtherm(maxtherm)
      dimension ttherm(maxtherm),ctherm(maxtherm)
      dimension xtherm(maxtherm),mtherm(maxtherm),ktherm(maxtherm)
      common/print/ip1,ip2,ip3,ip4,ip5,ip6
      common/dercom/vd,istatd
      dimension s(n),dsdt(n)
      dimension itype(natom)
      dimension moving(natom),mx(natom),my(natom),mz(natom)
      dimension grad(3,natom)
      dimension taxes(3)
      save
c      save ireport
      data ireport/0/

      if(nmove.gt.movmax) then
        write(6,*) ' - nmove too large in derivs - max=',movmax
        istatd=-5
        return
      end if
      if(natom.gt.natmax) then
        write(6,*) ' - natom too large in derivs - max=',natmax
        istatd=-6
        return
      end if

c  quick, one-time fluctuating force report
      if(ktherm(1).eq.4 .and. ireport.ne.1) then
        ireport=1
        boltz=3.167d-6
        alpha=ctherm(1)*1822.83d0*amu(itype(1))
        fmag=sqrt(alpha*2.0d0*boltz*ttherm(1)/dt)
c  fmag units should be a.u./sec   (here a.u. means momentum: aMu*bohr/atu)
c  converting to all a.u. should give hartree/bohr, then we go to per angstr.
        sigmaha=fmag*2.418d-17/0.529d0   ! this should be sigma in h/angstr.
        write(6,*)' Langevin thermostat info:'
        write(6,*)'   dp/dt = true force + Gaussian noise force -beta*p'
        write(6,*)'   variance of Gauss Force = 2kT*mass*beta/dt'
        write(6,*)'   Langevin alpha = beta*mass'
        write(6,*)'   temperature =',ttherm(1),'   dt =',dt
        write(6,*)'   beta (friction multiplying p) =',ctherm(1),' /sec'
        write(6,*)'   sigma of Gaussian force =',sigmaha,' (h/A)'
        write(6,*)'                           =',fmag, ' (a.u./sec)'
      end if

      if(n.ne.6*natom+ntherm) stop 'derivs: neq.ne.6*natom+ntherm'

      if(iwrk8.gt.0) then
        write(6,*) 'wrkcm8 conflict in derivs'
      else if(iwrk9.gt.0) then
        write(6,*) 'wrkcm9 conflict in derivs'
      else
        iwrk8=1
        iwrk9=1
      end if

c put moving atoms first, compute grad, order them back again
      call cpak(natom,s(1),s(1+natom),s(1+2*natom),itype,moving,work)
      lenw=lenw1
      call map3to1(natom,s(1),s(1+natom),s(1+2*natom),work)
      call map3to1(1,tax,tay,taz,taxes)
      call MPI_BARRIER(force_comm,ier)
      !if(rank_f.eq.0) write(217+spawnID,*)'spawnID ',spawnID,'IN GCALC'
      !if(rank_f.eq.0) call flushtad(217+spawnID)
      call gcalc(natom,nmove,work,itype,taxes,ietype,maxtyp,rcut,
     x                                        e,grad, hyperratio)
      call MPI_BARRIER(force_comm,ier)
      !if(rank_f.eq.0) write(217+spawnID,*)'spawnID ',spawnID,'OUT GCALC'
      !if(rank_f.eq.0) call flushtad(217+spawnID)
      vd=e
      call map1to3(natom,s(1),s(1+natom),s(1+2*natom),work)
      call map1to3(1,tax,tay,taz,taxes)
      call cunpak(natom,s(1),s(1+natom),s(1+2*natom),itype,moving,work)

      if(ip1.ge.2) write(6,20) ((grad(j,i),j=1,3),i=1,nmove)
   20 format('  gradient vector for moving atoms:'/ (3x,1p5d11.3) )

c  send back derivatives in the appropriate units
c  presently:
c       momentum is in a.u.
c       distance is in angstroms
c       gradient is in hartree/angstrom
c  units of dsdt();  appropriate for multiplication by dt in seconds:
c       distance deriv in angstroms/sec
c       momentum deriv in a.u./sec

c  all these gyrations with mx,my,mz allow user to turn on and off
c  atom motion in exactly the same way as for the molecular statics
c  note that now any atom with nonzero momentum will drift along,
c  even if it is not officially moving.  (7/18/89)

      cfac=0.529d0/2.418d-17
      im=0
      do 40 i=1,natom
      if(moving(i).ne.0) im=im+1
      rminv=1.0d0/(1822.83d0*amu(itype(i)))
      dsdt(i        )=s(i+3*natom)*rminv*cfac
      dsdt(i+1*natom)=s(i+4*natom)*rminv*cfac
      dsdt(i+2*natom)=s(i+5*natom)*rminv*cfac
      dsdt(i+3*natom)=0.0d0
      dsdt(i+4*natom)=0.0d0
      dsdt(i+5*natom)=0.0d0
cxxxxx check that following grad index conversion works OK
      if(mx(i).ne.0) dsdt(i+3*natom)=-grad(1,im)*cfac
      if(my(i).ne.0) dsdt(i+4*natom)=-grad(2,im)*cfac
      if(mz(i).ne.0) dsdt(i+5*natom)=-grad(3,im)*cfac
   40 continue
      istatd=0

c  calculate derivatives for thermostat(s)   (type stored in ktherm)
c      type 1 = Nose-Hoover
c      type 2 = Berendson (sp?)
c      type 3 = simple friction   dp/dt = F - beta*p
c      type 4 = friction + fluctuating force  dp/dt = F - beta*p + Fgauss
c               (Langevin dynamics)
c               Fgauss calc'd from fluct-diss. relation: beta=<F**2>*dt/2mkT
c      type 5 = Gauss (force added to achieve velocity rescaling - constant Tk)
c  and add in the drag term to the momentum derivatives
      boltz=3.167d-6
      do 70 i=1,ntherm
      if(ktherm(i).eq.0) go to 70
      ek=0.0d0
      if(mtherm(i).ne.0) then
        pxref = s(3*natom+mtherm(i))
        pyref = s(4*natom+mtherm(i))
        pzref = s(5*natom+mtherm(i))
      else
        pxref=0.0d0
        pyref=0.0d0
        pzref=0.0d0
      end if
c  loop over atoms in this thermostatted block
      sum5=0.0d0
      do 60 j=itherm(i),jtherm(i)
c   compute kinetic temperature in the moving reference frame of atom mtherm(i)
      rminv=1.0d0/(1822.83d0*amu(itype(i)))
      ek=ek + 0.5d0*(
     x  (s(3*natom+j) - pxref)**2 +
     x  (s(4*natom+j) - pyref)**2 +
     x  (s(5*natom+j) - pzref)**2
     x                               )*rminv
      sum5=sum5+               ! sum5 = sum of m*v*dvdt or p*dpdt/m  
     x ((s(3*natom+j) - pxref)*dsdt(3*natom+j) +
     x  (s(4*natom+j) - pyref)*dsdt(4*natom+j) +
     x  (s(5*natom+j) - pzref)*dsdt(5*natom+j)
     x                               )*rminv
   60 continue
      tk=ek*2.0d0/(3.0d0*boltz*(jtherm(i)-itherm(i)+1))

c set up drag coef.
      if(ktherm(i).eq.1) then
        eta=s(natom*6+i)
        drag=eta
        dsdt(natom*6+i)=ctherm(i)**2*(tk-ttherm(i))/ttherm(i)
      else if(ktherm(i).eq.2) then
        drag=ctherm(i)*(tk-ttherm(i))/ttherm(i)
        dsdt(natom*6+i)=0.0d0
      else if(ktherm(i).eq.3) then
        drag=ctherm(i)
        dsdt(natom*6+i)=0.0d0
      else if(ktherm(i).eq.4) then
        drag=ctherm(i)
        dsdt(natom*6+i)=0.0d0
      else if(ktherm(i).eq.5) then
        drag=sum5/(2.0d0*ek)    ! = - -sum[m*v*vdot]/sum[m*v*v]
        dsdt(natom*6+i)=0.0d0
      else
        write(6,*) i,ktherm(i)
        stop 'derivs: unrecognized thermostat type'
      end if

c  add drag to momentum derivatives in this block
      do 65 j=itherm(i),jtherm(i)
      if(mx(j).ne.0) dsdt(j+3*natom)=dsdt(j+3*natom)-drag*s(j+3*natom)
      if(my(j).ne.0) dsdt(j+4*natom)=dsdt(j+4*natom)-drag*s(j+4*natom)
      if(mz(j).ne.0) dsdt(j+5*natom)=dsdt(j+5*natom)-drag*s(j+5*natom)
   65 continue

c  for Langevin, add fluctuating force to momentum derivatives in this block
c    recall that beta is in inverse seconds, so must mult by mass to get alpha
      if(ktherm(i).eq.4) then
        iwarn=0
        do 66 j=itherm(i),jtherm(i)
c with luck, fmag units should come out to a.u./sec , so no cfac needed
        alpha=ctherm(i)*1822.83d0*amu(itype(j))
        if(itype(j).ne.1) iwarn=1
        fmag=sqrt(alpha*2.0d0*boltz*ttherm(i)/dt)
        call MPI_BARRIER(force_comm,ier)
!        if(rank_f.eq.0) write(217+spawnID,*)
!     x                   'spawnID ',spawnID,' IN aa.. '
        !if(rank_f.eq.0)then
          if(mx(j).ne.0) dsdt(j+3*natom)=dsdt(j+3*natom)+fmag*gasdev(0)
          if(my(j).ne.0) dsdt(j+4*natom)=dsdt(j+4*natom)+fmag*gasdev(0)
          if(mz(j).ne.0) dsdt(j+5*natom)=dsdt(j+5*natom)+fmag*gasdev(0)
        !endif
        !call MPI_BCAST(dsdt(j+3*natom),1,MPI_REAL8,0,force_comm,ier)
        !call MPI_BCAST(dsdt(j+4*natom),1,MPI_REAL8,0,force_comm,ier)
        !call MPI_BCAST(dsdt(j+5*natom),1,MPI_REAL8,0,force_comm,ier)
!        if(rank_f.eq.0) write(217+spawnID,*)
!     x                  'spawnID ',spawnID,' OUT aa.. '
   66   continue
        ! BCAST all of dsdt at once:
        call MPI_BCAST(dsdt(1),n,MPI_REAL8,0,force_comm,ier)
        if(iwarn.ne.0) 
     x    write(6,*) 'note - nonstandard Langevin dynamics for ntype>1'
      end if

   70 continue

      iwrk8=0
      iwrk9=0
      return
      end


      function gasdev(idum)
c returns Gaussian-distributed random number, unit variance
c uses Box-Muller transformation;  see p. 203 of Num. Rec., 1st ed.
      implicit real*8(a-h,o-z)
      save
      common/gasdevcom/iset
      if (iset.ne.1) then
1       v1=2.0d0*prngen(0)-1.0d0   ! pick coords in square (-1,1)x(-1,1)
        v2=2.0d0*prngen(0)-1.0d0
        r=v1**2+v2**2
        if(r.ge.1.0d0) go to 1     ! discard pairs that are not in unit circle
        fac=sqrt(-2.0d0*log(r)/r)  ! Box-Muller transformation
        gset=v1*fac                ! save this one for next time
        gasdev=v2*fac              ! return this one
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      end  

      subroutine bivariate(sig1,sig2,c12,g1,g2)
c returns 2 random numbers from bivariate Gaussian dist. with specified
c  std's and correlation coef.
c  i.e., 
c p(g1,g2)/norm 
c = exp(-(sig2**2*g1**2+sig1**2*g2**2)/(2[sig1**2+sig2**2-(c12*sig1*sig2)**2])
c = exp(-(var2*g1**2+var1*g2**2)/(2(var1+var2-c12**2*var1*var2))
c e.g.,  see p 28 of Wax (Chandrasekhar article, Eq.(178) )
c uses Box-Muller transformation;  see p. 203 of Num. Rec., 1st ed.
c and then covariance formula from App. G of Allen and Tildesley, p. 348.
      implicit real*8(a-h,o-z)
      save
1     v1=2.0d0*prngen(0)-1.0d0   ! pick coords in square (-1,1)x(-1,1)
      v2=2.0d0*prngen(0)-1.0d0
      r=v1**2+v2**2
      if(r.ge.1.0d0) go to 1     ! discard pairs that are not in unit circle
      fac=sqrt(-2.0d0*log(r)/r)  ! Box-Muller transformation
      gset1=v1*fac      
      gset2=v2*fac     

      g1=sig1*gset1
      g2=sig2*(c12*gset1 + sqrt(1.0d0-c12**2)*gset2)
      return
      end  

      SUBROUTINE UABORT(string)
      character*(*) string
      write(6,*) string
      if(i.eq.i) stop 'stopping in uabort'
      return
      end


      function prngen ( init )
c
c: >>> prngen -- portable random number generator
c
c            Note:  this function is lifted from the routine RAN3
c                   described by Press et. al. in _Numerical Recipes_
c                   I chose this function because of its portability and
c                   because it only requires one floating point multiply
c                   after the initialization call.
c
c     The routine is initialized (or reinitialized) if it is called with
c         a negative value for init.  Unlike the routine in
c         _Numerical Recipes_, the value of init is NOT altered by the
c         function.  Similarly, prngen returns ZERO on an initialization
c         call, and only returns random numbers when init.ge.0.
c
c     implicit none
      real*8 prngen
c ---                            !       avoids trouble
c
      integer   magic
c ---                            !       Knuth's magic number
      parameter (magic=55)
      integer   magoff
c ---                            !       another magic number
      parameter (magoff=31)
      integer   zero
c ---                            !       constant = 0
      parameter (zero=0)
      integer   BIGINT
c ---                            !       any large integer
      parameter (BIGINT=1000000000)
      integer   SEED
c ---                            !       any seed integer (.lt.BIGINT)
      parameter (SEED  = 161803398)
c
      logical   FIRST
c ---                            !       flag for first call
      integer   Next, NextP
c ---                            !       pointers to table slots
      integer   Table(magic)
c ---                            !       table of seed derivatives
      integer   i,j,k,ii
c ---                            !       miscellaneous counters
c
      real*8      scale
c ---                            !       scale factor to get reals
      parameter (scale=1e0/BIGINT)
c
      save FIRST, Table, Next, NextP
c
      data FIRST / .true. /
c
c--------------   start initialization section  here  --------------------
c
      IF ( FIRST .or. init.lt.zero ) THEN
        j = iabs(mod(SEED-iabs(init), BIGINT))   ! outer iabs added by afv 6/9/98
        Table(magic) = j
        k = 1
c
        do 10 i=1,magic-1
          ii        = mod(21*i,magic)
          Table(ii) = k
          k         = j - k
          if ( k.lt.zero ) k = k + BIGINT
          j = Table(ii)
   10   continue
c
        do 20 k = 1, 4
          do 30 i = 1, magic
            Table(i) = Table(i) - Table(1+mod(i+magoff-1,magic))
            if ( Table(i).lt.zero ) Table(i) = Table(i) + BIGINT
   30     continue
   20   continue
c
        Next  = 0
        NextP = magoff
c
        IF ( init.lt.zero ) THEN
          prngen = 0.0d0
          FIRST  = .false.
          return
        END IF
c
        FIRST = .false.
      END IF
c
c-------------  here's where we go if we are returning a value  -------
c
      Next  = Next  + 1
      if ( Next  .gt. magic ) Next  = 1
c
      NextP = NextP + 1
      if ( NextP .gt. magic ) NextP = 1
c
      j = Table(Next) - Table(NextP)
      if ( j.lt.zero ) j = j + BIGINT
c
      Table(Next) = j
      prngen = scale*j
c
      return
      end

      subroutine cpak(natom,x,y,z,itype, moving,w)
c  reorder x,y,z,itype arrays to put moving atoms first
c  work array (w) must be 4*natom real*8 long
      implicit real*8(a-h,o-z)
      dimension x(natom),y(natom),z(natom),itype(natom)
      dimension moving(natom)
      dimension w(*)
c                              ! work array
      nmove=0
      do 10 i=1,natom
      w(i)=x(i)
      w(i+natom)=y(i)
      w(i+natom+natom)=z(i)
      w(i+natom+natom+natom)=float(itype(i))
      if(moving(i).ne.0) nmove=nmove+1
   10 continue

      m=0
      nm=0
      do 20 i=1,natom
      if(moving(i).ne.0) then
        m=m+1
        idx=m
      else
        nm=nm+1
        idx=nmove+nm
      end if
      x(idx)=w(i)
      y(idx)=w(i+natom)
      z(idx)=w(i+natom+natom)
      itype(idx)=int(w(i+natom+natom+natom))
   20 continue
      if(m+nm.ne.natom) stop 'abort - screwed up in cpak'
      return
c
      entry cunpak(natom,x,y,z,itype, moving,w)
c  undoes what cpak does.
c
      nmove=0
      do 30 i=1,natom
      w(i)=x(i)
      w(i+natom)=y(i)
      w(i+natom+natom)=z(i)
      w(i+natom+natom+natom)=float(itype(i))
      if(moving(i).ne.0) nmove=nmove+1
   30 continue

      m=0
      nm=0
      do 40 i=1,natom
      if(moving(i).ne.0) then
        m=m+1
        idx=m
      else
        nm=nm+1
        idx=nmove+nm
      end if
      x(i)=w(idx)
      y(i)=w(idx+natom)
      z(i)=w(idx+natom+natom)
      itype(i)=int(w(idx+natom+natom+natom))
   40 continue
      if(m+nm.ne.natom) stop 'abort - screwed up in cunpak'
      return
      end

      subroutine putclsnew(natom,ntitle,ctitle,xyz,pxyz,itype,
     x                           taxes,iunit,iwrt,iop)
c 9/24/02 - changed to allow iop.  E.g., iop=1 does not write momenta (clsman 'outs' style)
c 12/20/01 - changed over to single-array format (xyz,pxyz,taxes)
c putclsnew is just like putcls, but uses a character title string instead of real*8 afv 12/15/01
c ctitle must be dimensioned character*80 and as an array ntitle long
      implicit real*8(a-h,o-z)
      dimension xyz(3,*),itype(*),pxyz(3,*),taxes(3)
      character*80 ctitle(*),fmt

      if(ntitle.gt.20) then
        write(6,10) ntitle
   10   format('  warning -- file not written, ntitle too large:',i5)
        return
      end if

c  write out a cluster file

      write(iunit,30) natom,iop,ntitle
      write(iunit,41) (ctitle(i),i=1,ntitle)
      write(iunit,50) (taxes(j),j=1,3)
      if(iop.eq.0) then  
        write(iunit,50) ((xyz(j,i),j=1,3),itype(i),
     x                  (pxyz(j,i),j=1,3),i=1,natom)
      else if(iop.eq.1) then  
        fmt='(3f12.6,i2)' 
        write(iunit,55) fmt
        write(iunit,fmt)
     x                  ((xyz(j,i),j=1,3),itype(i),i=1,natom)
      else if(iop.eq.2) then
        fmt='(3f9.3,i2,3f9.3)'
        write(iunit,55) fmt
        write(iunit,fmt) ((xyz(j,i),j=1,3),itype(i),
     x                    (pxyz(j,i),j=1,3),i=1,natom)
      else
        write(6,*) 'putclsnew - unrecognized iop - stopping',iop
        stop 'putclsnew - unrecognized iop - stopping'
      end if

   30 format(i5,i2,i3)
   40 format(10a8)
   41 format(a80)
   50 format(3f20.8,i5,3f20.8)
   55 format(a)


      if(iwrt.ge.1) then
        ntype=imax(itype,natom)
        write(6,70) natom,ntype,(ctitle(i),i=1,ntitle)
   70   format(' cluster file written out:  natom =',i6,'  ntype =',i2/
     x         (1x,a80))
        if(iwrt.ge.2) write(6,80) (taxes(j),j=1,3)
   80   format(/' tax =',f18.8,'  tay =',f18.8,'  taz =',f18.8)
      end if
      return
      end

      subroutine seedit(iseed)
c:  Seeds the random number generator with iseed.
c:  If iseed=0, iseed is set in here and returned for reference.
c:  WARNING - iseed gets changed if iseed=0  - - do not pass a constant integer!!!
c:
      implicit real*8(a-h,o-z)
      save
      common/gasdevcom/iset
      integer time
c
      if(iseed.eq.0) then
        iseed = time()*2 - 1
      end if

      iseed=mod(abs(iseed),100000000)  ! added 10/23/01 to fix seeds too large for prngen

c  after this first call to prngen, all future calls should pass zero arg.
      iseed=-abs(iseed)
      xdum=prngen(iseed)
      iset=0
      return
      end


      subroutine gaussp(natom,nmove,pxyz,itype,maxtyp,amu,ntype,
     x                   moving,mx,my,mz,temperature,iop,iwrt)
c  modular version   8/97
c  Set up gaussian random momenta for moving atoms
c  if iop.eq.0 - zero out the center of mass motion
c  if iop.eq.1 - do not zero out the center of mass motion
      use mod_mpi
      implicit real*8 (a-h,o-z)
      save
c      parameter(mxatom=10000)
      include 'parameters.h'
      dimension pxyz(3,natom)
      dimension px(mxatom),py(mxatom),pz(mxatom),itype(natom)
      dimension amu(maxtyp)
      dimension moving(natom),mx(natom),my(natom),mz(natom)
      dimension psigma(10),rmass(10)

      call map1to3(natom,px,py,pz,pxyz)

      if(ntype.gt.10) stop 'abort - redimension psigma in gaussp'
      do 5 it=1,ntype
      rmass(it)=amu(it)*1822.83d0
    5 psigma(it)=sqrt(rmass(it)*temperature*3.167d-6)

      ek = 0.0d0
      do 10 i=1,natom
      if(moving(i).eq.1) then
        px(i)=gasdev(0)*psigma(itype(i))
        if(mx(i).ne.1) px(i)=0.0d0
        py(i)=gasdev(0)*psigma(itype(i))
        if(my(i).ne.1) py(i)=0.0d0
        pz(i)=gasdev(0)*psigma(itype(i))
        if(mz(i).ne.1) pz(i)=0.0d0
        ek = ek + 0.5d0*(px(i)**2+py(i)**2+pz(i)**2)/rmass(itype(i))
      end if
   10 continue

c  Center of mass correction

      vcmx = 0.0d0
      vcmy = 0.0d0
      vcmz = 0.0d0
      tmassx=0.0d0
      tmassy=0.0d0
      tmassz=0.0d0
      do 40 i=1,natom
      if(moving(i).eq.1) then
        vcmx = vcmx + px(i)
        vcmy = vcmy + py(i)
        vcmz = vcmz + pz(i)
        if(mx(i).eq.1) tmassx=tmassx+rmass(itype(i))
        if(my(i).eq.1) tmassy=tmassy+rmass(itype(i))
        if(mz(i).eq.1) tmassz=tmassz+rmass(itype(i))
      end if
   40 continue
      tmass=tmass*1822.83d0
      vcmx = vcmx/tmassx
      vcmy = vcmy/tmassy
      vcmz = vcmz/tmassz

      ekp = 0.0d0
      do 50 i = 1,natom
      if(moving(i).eq.1) then
        if(iop.eq.0) then
          if(mx(i).eq.1) px(i)=px(i) - rmass(itype(i))*vcmx
          if(my(i).eq.1) py(i)=py(i) - rmass(itype(i))*vcmy
          if(mz(i).eq.1) pz(i)=pz(i) - rmass(itype(i))*vcmz
        end if
        ekp = ekp + (px(i)**2+py(i)**2+pz(i)**2)/rmass(itype(i))/2.0d0
      end if
   50 continue
      tfac=1.0d0/(3.0d0/2.0d0*3.167d-6)/float(nmove)
      if(iop.eq.0) then
        if(iwrt.ge.2) then
          write(6,70) nmove,ekp*tfac,ek*tfac
   70     format('  kin. temp of ',i4,' moving atoms =',f12.3,
     x             ' K (before c.m. correction =',f12.3,' K)')
        end if
      else if(iop.eq.1) then
        if(iwrt.ge.2) then
          write(6,80) nmove,ekp*tfac
   80     format('  kin. temp of ',i4,' moving atoms =',f12.3,
     x                              ' K (no c.m. correction)')
        end if
      else
        write(6,*) 'gaussp - iop =',iop,'  ???'
      end if

      call map3to1(natom,px,py,pz,pxyz)

      return
      end

      subroutine pbccum(natom,xyz,taxes,cdxyz)
c: move atoms into primary period, accumulating change in cdx,cdy,cdz
c: also see pbcuncum
c: user should zero cdx,cdy,cdz before first call here
      implicit real*8(a-h,o-z)
      dimension xyz(3,*),cdxyz(3,*)
      dimension taxes(3)
      dimension axes(3)

      do j=1,3
        axes(j)=taxes(j)/2.0d0
      end do

      do 20 i=1,natom
      do 10 j=1,3
      xold=xyz(j,i)
      xyz(j,i)=mod(xyz(j,i),taxes(j)) + axes(j) - sign(axes(j),xyz(j,i))
      cdxyz(j,i)=cdxyz(j,i)+xold-xyz(j,i)
   10 continue
   20 continue
      return
      end

      subroutine pbcuncum(natom,xyz,taxes,cdxyz)
c: move atoms back to their absolute positions (opposite of pbccum)
c: this uses cdx,cdy,cdz, and then sets those to zero
      implicit real*8(a-h,o-z)
      dimension xyz(3,*),cdxyz(3,*)

      do 20 i=1,natom
      do 10 j=1,3
      xyz(j,i)=xyz(j,i)+cdxyz(j,i)
      cdxyz(j,i)=0.0d0
   10 continue
   20 continue

      return
      end

      subroutine pbccumfake(natom,xyz,taxes,cdxyz)
c: move atoms into primary period, using -cdx,-cdy,-cdz
c: "fake" means that no check is made that this actually moves
c: the atoms into the primary period, and cdx,cdy,cdz are left UNCHANGED
c: Also see pbcuncumfake, which does the exact opposite.
c: The point of these routines is to apply the accumulated displacements
c: to what may be a slightly different cluster, without disrupting the
c: displacements.
c: tax,tay,taz are not used, but are passed so the call is identical to pbccum
      implicit real*8(a-h,o-z)
      dimension xyz(3,*),cdxyz(3,*)

      do 20 i=1,natom
      do 10 j=1,3
      xyz(j,i)=xyz(j,i)-cdxyz(j,i)
   10 continue
   20 continue
      return
      end

      subroutine pbcuncumfake(natom,xyz,taxes,cdxyz)
c: move atoms out of primary period, using +cdx,+cdy,+cdz
c: "fake" means that no check is made that this actually moves
c: the atoms into the primary period, and cdx,cdy,cdz are left UNCHANGED.
c: Also see pbccumfake, which does the exact opposite.
c: The point of these routines is to apply the accumulated displacements
c: to what may be a slightly different cluster, without disrupting the
c: displacements.
c: tax,tay,taz are not used, but are passed so the call is identical to pbccum
      implicit real*8(a-h,o-z)
      dimension xyz(3,*),cdxyz(3,*)

      do 20 i=1,natom
      do 10 j=1,3
      xyz(j,i)=xyz(j,i)+cdxyz(j,i)
   10 continue
   20 continue
      return
      end
      
