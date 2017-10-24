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
c     v1: this does regular dimer search, always randomizing dimer direction
c     v2: allows initial dimer direction to be passed in, also returns it
c     v3: returns new basin in x1
c**********************************************************************
c meaning of imode variable:
c imode=1: do regular dimer search, redoing dimers that fail criteria
c imode=2: do not redo dimers that fail criteria, return instead, but do
c full dimer search
c imode=3: only minimize direction, returning with eig, do not move
c dimer 
c**********************************************************************


      subroutine dimer_search(natom,nmove,idimerbail,emin,eminbar
     +    ,idimerdof,ndimerdof,rdimerdist,ndimercoord,idimeroverunder
     +    ,rdimerdisp,x1,x,xmin,xyzprev,energy,rdimer,ityp,taxes,ietype
     +    ,maxtyp,rcut,hyperratio,barrevev,ereverse,gfac,irotalign
     +    ,itermax,ierr,iranx,irandir,transcrit,itranscrit,gcritin
     +    ,dvcritin,drcritin,ipr,imode,ifundc,eig)

      use mod_mpi
      implicit real*8(a-h,o-z)
c      parameter(mxatom=10000,nneighmax=1500)
      include 'parameters.h'
      dimension ityp(natom)
      dimension rcut(maxtyp,maxtyp)
      dimension x1(3*natom),x(3*natom),xmin(3*natom),xyzprev(3*natom)
      dimension barrevev(nneighmax)

c x1 is original minimum, x is working saddle array, xmin is minimum
c again (maybe redundant?), xyzprev is previous minimum

      dimension grad(3*mxatom),rdimer(3*natom),xwork(3*mxatom),
     +    xwork2(3*mxatom)
      dimension F(3*natom),Fold(3*natom),F2(3*natom)
      dimension G(3*natom),Gold(3*natom),Gu(3*natom)
      dimension R(3*natom),Rold(3*natom),R2(3*natom)
      dimension dF(3*natom)
!      dimension dF(3*natom),Fp1(3*natom),Fp2(3*natom)
      dimension taxes(3)
      dimension icoord(mxatom)
      character*10 chri
      character*80 filnam ! uuuu unused

      common/callcom/ngradcall

      data itag/0/

      save

      nitermaxfail = 0
      nitermaxfailmax = 1

      ierr=0
      nredocount=0

      ev=27.21d0
c----------------------------------------------------------------------
c dimer parameters that may be passed later
      drdimer=0.0001d0
      drdimer=0.001d0
c      rcoordcrit=2.3d0      ! cutoff for coordination and neighbors of atoms
c      rranx=0.6d0           ! distance to displace atoms (+/- this number)
      rmaxmove=0.2d0         ! maximum distance allowed to move
      dR=0.001d0            ! finite difference translation step size
c----------------------------------------------------------------------
      rcoordcrit=rdimerdist
      rranx=rdimerdisp
      factdimer=10.0d0
      factdimer=1.0d0

      gcrit=gcritin*factdimer
      dvcrit=dvcritin*factdimer
      drcrit=drcritin*factdimer
      
      if(rank_f.eq.0)then
       if(ipr.ge.4) write(6,*) 'dimer_report: drdimer,rmaxmove,dR: ',
     x                         drdimer,rmaxmove,dR
      endif
     
 111  continue

      e100=0.0d0

      call flushtad(6)
      if (iranx.eq.1) then
         if (idimerdof.eq.2) then
            transcrittemp=rcoordcrit
            itranscrittemp=1
            call coordination(natom,x1,icoord,transcrittemp,taxes,rcut
     +          ,ityp,itranscrittemp,maxtyp)
         endif
         call dimer_ranx(natom,nmove,x1,x,ipr,idimerdof
     +       ,ndimerdof,ndimercoord,idimeroverunder,taxes,icoord
     +       ,rcoordcrit,rranx)
         if(rank_f.eq.0)then
           if (ipr.ge.4) write(6,*) 'dimer_search: randomized x'
         endif
      endif

      if (iranx.eq.2) then
         call dimer_ranx_hole(natom,nmove,x,xmin,xyzprev,ipr,idimerdof
     +       ,ndimerdof,ndimercoord,idimeroverunder,taxes,rranx
     +       ,transcrit,rcoordcrit)
         if((ipr.ge.4).and.(rank_f.eq.0))
     +       write(6,*) 'dimer_search: randomized x near hole'
         rdimermag=vecdot(rdimer,rdimer,3*natom)
         if((ipr.ge.4).and.(rank_f.eq.0))
     +       write(6,*) 'dimer_search: rdimermag: ',rdimermag
      endif

      if (iranx.eq.3) then
         do i=1,natom
            icoord(i)=0
         enddo
         iunit=25
         if(rank_f.eq.0)
     +       open(unit=iunit,file='dimer.activeatoms',status='old')
         if(rank_f.eq.0) read(iunit,*) nactive
         call MPI_BCAST(nactive,1,MPI_INTEGER,0,force_comm,ier)
         do i=1,nactive
            if(rank_f.eq.0) read(iunit,*) iactiveatom
            call MPI_BCAST(iactiveatom,1,MPI_INTEGER,0,force_comm,ier)
            icoord(iactiveatom)=1
         enddo
         if(rank_f.eq.0) close(unit=iunit)
         if(rank_f.eq.0)
     +       write(6,*) 'iranx.eq.3: read ',nactive,' active atoms'
         idimerdoftemp=2 ! based on atom number, not dof number
         idimeroverundertemp=0  ! over coordinated atoms 
         ndimercoordtemp=0 ! coordination of 0 (use as flag rather than coordination)
         call dimer_ranx(natom,nmove,x1,x,ipr,idimerdoftemp
     +       ,ndimerdof,ndimercoordtemp,idimeroverundertemp,taxes,icoord
     +       ,rcoordcrit,rranx)
      endif

      itag=itag+1
      filnam='dimerinit.'//chri(itag)//'.dat'
      if (ipr.ge.10.or.(ipr.ge.10.and.iranx.eq.2))then
        call storefile(natom,x,x,ityp,taxes,filnam)
      endif
      if (irandir.eq.1) then
         call dimer_initialize(natom,nmove,x1,x,idimerdof,ndimerdof
     +       ,rdimer)
         if(rank_f.eq.0)then
          if (ipr.ge.4) write(6,*) 'dimer_search: randomized dimer'
         endif
      else
c     make sure that last elements are zero
         do i=3*nmove+1,3*natom
            rdimer(i)=0.0d0
         enddo
         call vecrenorm(rdimer,3*natom)
         if(rank_f.eq.0)then
          if(ipr.ge.4) write(6,*) 'dimer_search: using passed direction'
         endif
      endif

! Dimer convergence loop      
      eprev=0.0d0
      ngradcallprev=ngradcall
      ngradcallorig=ngradcall
      icginit=1
      ngraddimer=0
      ngraddimerprev=ngradcall
      do i=1,itermax
         
! Dimer rotation
         call gcalc(natom,nmove,x,ityp,taxes,
     x       ietype,maxtyp,rcut,energy,grad,hyperratio)
     
!       Steepest descent rotation
c         call dimer_min_torque(natom,nmove,x,rdimer,drdimer,eig,ityp
c     +       ,taxes,ietype,maxtyp,rcut,energy,grad,hyperratio)
!       Quasi-Newton rotation
         call dimer_min_torque_2(imode,natom,nmove,x,rdimer,drdimer,eig
     +       ,ityp,taxes,ietype,maxtyp,rcut,energy,grad,hyperratio)
     
         if (imode.eq.3) return

! Steepest descent dimer translation
!         gdotdr=-vecdot(grad,rdimer,3*nmove)
!         call vecypax(grad,gdotdr,rdimer,3*nmove) ! zero grad along rdimer
!         gfactor=gdotdr/(gfac*(abs(eig)+0.05d0)) ! uphill move, tempered newt-raph-style along neg. curv
!         call vecypax(grad,gfactor,rdimer,3*nmove)
!         gmag=sqrt(vecdot(grad,grad,3*nmove))
!         dxmax=0.0d0
!         do j=1,nmove
!            dx=0.0d0
!            do k=1,3
!               xx=x((j-1)*3+k)
!               x((j-1)*3+k)=xx-gfac*grad((j-1)*3+k)
!               dx=dx+(x((j-1)*3+k)-xx)**2
!            enddo
!            dxmax=max(dxmax,sqrt(dx))
!         enddo

! Conjugate gradient dimer translation
         do j=nmove*3+1,natom*3
           rdimer(j)=0.0d0
         enddo
         call vecrenorm(rdimer,3*nmove)
         gmag=sqrt(vecdot(grad,grad,3*nmove))
         gdotdr=-vecdot(grad,rdimer,3*nmove)
         do j=1,3*natom
           R(j)=x(j)
           F(j)=-grad(j)
         enddo
         do j=nmove*3+1,natom*3
           F(j)=0.0d0
         enddo
         fdotdr=vecdot(F,rdimer,3*nmove)
         if(eig.lt.0) then
           do j=1,3*nmove
             F(j)=F(j)-rdimer(j)*fdotdr*2.0d0
           enddo
         else
            do j=1,3*natom
             F(j)=-rdimer(j)*fdotdr
           enddo
         endif
         if(icginit.eq.1) then
           icginit=0
!           write(6,*) 'dimer_search: CG Init'
           Gam=0
           do j=1,3*natom
             G(j)=F(j)
             Gold(j)=F(j)
             Fold(j)=F(j)
           enddo
         endif
         do j=nmove*3+1,natom*3
           G(j)=0.0d0
           Gu(j)=0.0d0
           Gold(j)=0.0d0
           Fold(j)=0.0d0
         enddo
         a1=abs(vecdot(F,Fold,3*nmove))
         a2=vecdot(Fold,Fold,3*nmove)
         if(a1.le.(0.5d0*a2)) then
           do j=1,3*nmove
             dF(j)=F(j)-Fold(j)
           enddo
           Gam=vecdot(F,dF,3*nmove)/a2
         else
!           write(6,*) 'dimer_search: Resetting CG in translation'
           Gam=0
         endif
!         write(6,*) 'dimer_search: a1 ',a1,' a2 ',a2,' Gam ',Gam
         do j=1,3*nmove
           G(j)=F(j)+Gold(j)*Gam
           Gu(j)=G(j)
         enddo
         call vecrenorm(Gu,3*nmove)
         do j=1,3*natom
           Fold(j)=F(j)
           Gold(j)=G(j)
           Rold(j)=R(j)
           R2(j)=R(j)+Gu(j)*dR
         enddo
         call gcalc(natom,nmove,R2,ityp,taxes,
     x       ietype,maxtyp,rcut,energy2,grad,hyperratio)
         do j=1,3*natom
           F2(j)=-grad(j)
         enddo
         do j=nmove*3+1,natom*3
           F2(j)=0.0d0
         enddo
         fdotdr=vecdot(F2,rdimer,3*nmove)         
         if(eig.lt.0) then
           do j=1,3*nmove
             F2(j)=F2(j)-rdimer(j)*fdotdr*2.0d0
           enddo
         else
            do j=1,3*nmove
             F2(j)=-rdimer(j)*fdotdr
           enddo
         endif
         Frp1=vecdot(F,Gu,3*nmove)
         Frp2=vecdot(F2,Gu,3*nmove)
!         do j=1,3*nmove
!           Fp1(j)=Gu(j)*Frp1
!           Fp2(j)=Gu(j)*Frp2
!         enddo
         Frp=(Frp1+Frp2)/2d0
         CR=(Frp1-Frp2)/dR
!         write(6,*) 'dimer_search: CR ',CR,' Frp1 ',Frp1,' Frp2 ',Frp2         
         if(CR.lt.0d0) then
           RdR=rMaxMove
!           write(6,*) 'dimer_search: CR<0, RdR = ',RdR        
         else
           RdR=Frp/CR
!           write(6,*) 'dimer_search: CR>0, RdR = ',RdR        
           if(abs(RdR).ge.rMaxMove) then
             RdR=sign(rMaxMove,RdR)
           else
             RdR=RdR+dR/2d0
           endif
         endif
!         write(6,*) 'dimer_search: Finally RdR = ',RdR        
         do j=1,3*natom
           R(j)=R(j)+Gu(j)*RdR
           x(j)=R(j)
         enddo
         dxmax=0
         do j=1,nmove
            dx=0.0d0
            do k=1,3
               dx=dx+(R((j-1)*3+k)-Rold((j-1)*3+k))**2
            enddo
            dxmax=max(dxmax,sqrt(dx))
         enddo
!         write(6,*) 'dimer_search: dxmax ',dxmax,' CR ',CR,' RdR ',RdR         
! Done dimer translation         
         
         ngraddimer=ngraddimer+ngradcall-ngraddimerprev
         ngraddimerprev=ngradcall

         dv=abs(energy-eprev)
         eprev=energy

         if (i.eq.100) e100=energy

         if (mod(i-1,10).eq.0.and.ipr.ge.4) then
            ndelgrad=ngradcall-ngradcallprev
            ngradcallprev=ngradcall
            eratio=(energy-emin)*27.21d0/eminbar
            if(rank_f.eq.0) then
             write(6,62) i,energy*27.21d0,dv,dxmax,gmag,gdotdr,eig
     $          ,ndelgrad,ngraddimer,eratio,e100
            endif
c            write(6,*) 'dimer_debug: ',energy,emin,eminbar,energy-emin
c     +          ,(energy-emin)*27.21d0/eminbar
         endif

 62      format('dimer_search: ',i5,f12.5,1p3d12.4,0p2f12.5,i12,i12,
     +       0pf12.5,f12.5)
         if (dv.lt.dvcrit.and.dxmax.lt.drcrit.and.gmag.lt.gcrit) goto
     +       100

         if (abs(eig).lt.0.0001d0) then
            nzero=nzero+1
            if (nzero.gt.0) then
               if(rank_f.eq.0) then
                write(6,*) 'dimer_search: seem to be on zero-mode: ',eig
               endif
               if (imode.eq.2) then
                  ierr=10
                  if(rank_f.eq.0) then
                   write(6,*)'dimer_search: going to next reused saddle'
                  endif
                  return
               endif
               if(rank_f.eq.0) then
                write(6,*) 'dimer_search: redoing dimer'
               endif
               nzero=0
               nredocount=nredocount+1
               if (nredocount.gt.100) then
                  if(rank_f.eq.0) then
                   write(6,*) 'STOP: nredocount.gt.100 in dimer!'
                  endif
                  stop 'nredocount.gt.100 in dimer!'
               endif
               goto 111
            endif
         else
            nzero=0
         endif
         
         if (ngraddimer.gt.itermax) then
            nitermaxfail = nitermaxfail + 1
            if(rank_f.eq.0) then
             write(6,*) 'SpawnID ',spawnID
     +          ,' - dimer_search: ngraddimer.gt.itermax: '
     +          ,ngraddimer,itermax
            endif
            if (imode.eq.2) then
               ierr=10
               if(rank_f.eq.0) then
                 write(6,*) 'dimer_search: going to next reused saddle'
               endif
               return
            endif
            if (ipr.ge.2) then 
               if(rank_f.eq.0) then
                write(6,61)'dimer_search: Search Did Not Converge: ',dv
     +             ,dxmax,gmag
               endif
               if (imode.eq.2) then
                  ierr=10
                  if(rank_f.eq.0) then
                   write(6,*)'dimer_search: going to next reused saddle'
                  endif
                  return
               endif
               if(nitermaxfail.lt.nitermaxfailmax)then
                if(rank_f.eq.0) write(6,*) 'dimer_search: redoing dimer'
               endif
            endif
            if(nitermaxfail.lt.nitermaxfailmax) goto 111
            ierr=10
            return
         endif
         
         if (idimerbail.eq.1) then ! stop dimer search if certain criteria are met
            if (i.gt.100.and.eig.gt.0.0d0) then
               if(rank_f.eq.0) write(6,*)
     +             'dimer_search: did 100 steps, eig still positive'
     +             ,' bailing on this search'
               if (imode.eq.2) then
                  ierr=10
                  if(rank_f.eq.0) write(6,*) 'dimer_search: returning'
                  return
               endif
               if(rank_f.eq.0) write(6,*) 'dimer_search: redoing dimer'
               goto 111
            endif
            ebar=(energy-emin)*27.21d0
            eratio=(energy-emin)*27.21d0/eminbar
            ebardif=ebar-eminbar
c            if (eig.lt.0.0d0.and.eratio.gt.2.0d0.and.gmag.lt.0.01d0)
            if (eig.lt.0.0d0.and.ebar.gt.0.5d0.and.ebardif.gt.5.0d0.and
     +          .gmag.lt.0.01d0)then
               if(rank_f.eq.0) write(6,*)
     +             'dimer_search: hit bailing condition, '
     +             ,'bailing on this dimer: ',ebar,ebardif
               if (imode.eq.2) then
                  if(rank_f.eq.0) write(6,*) 'dimer_search: returning'
                  ierr=10
                  return
               endif
               if(rank_f.eq.0) write(6,*) 'dimer_search: redoing dimer'
               goto 111
            endif
         endif
      enddo

 100  continue
      if (i.lt.itermax) then
         if(rank_f.eq.0) then
          if (ipr.ge.4) write(6,*) 'dimer_search: Search Converged'
         endif
         if (ipr.ge.4) then
           eratio=(energy-emin)*27.21d0/eminbar
           if(rank_f.eq.0) then
            write(6,62) i,energy*27.21d0,dv,dxmax,gmag,gdotdr,eig
     $          ,ndelgrad,ngraddimer,eratio,e100
           endif
         endif
      else
         ierr=1
         if (ipr.ge.2) then
            if(rank_f.eq.0) then
             write(6,61) 'dimer_search: Search Did Not Converge: ',dv
     +          ,dxmax,gmag
            endif
            if (imode.eq.2) then
              ierr=10
              if(rank_f.eq.0) then
               write(6,*) 'dimer_search: going to next reused saddle'
              endif
              return
            endif
            if(rank_f.eq.0) write(6,*) 'dimer_search: redoing dimer'
         endif
         
 61      format(a35,3d12.5)
         goto 111
      endif

      if(rank_f.eq.0) write(6,*) 'dimer_search: Total Force Calls '
     &            ,ngradcall-ngradcallorig
! Dimer search done

      filnam='dimerfinal.'//chri(itag)//'.dat'
      if (ipr.ge.10) call storefile(natom,x,x,ityp,taxes,filnam)

c if doing from NEB, don't want to do roll check, as NEB will do that
      if (imode.eq.1) return

c if nmove=natom, shift saddle to remove any translation
      if (nmove.eq.natom) then
         call shiftalign(x1,x,natom)
      endif

c do roll check to make sure dimer connects original minimum
      ngradrchold=ngradcall
      do j=1,3*natom
         xwork(j)=x(j)-x1(j)
      enddo
      xdotdimer=vecdot(xwork,rdimer,3*natom)
      isignd=-1
      if (xdotdimer.lt.0) isignd=1
      do j=1,3*natom
         xwork(j)=x(j)+0.1d0*rdimer(j)*isignd
      enddo
      imodedc=2
      ifundcuse=ifundc
      if (ifundcuse.eq.5) ifundcuse=0
      ereverseuse=0.0d0         ! not comparing to known neighbors, so don't use reverse barrier stuff
      if (imode.eq.0) then
         call descent_check(natom,nmove,1,xwork,x1,ityp,taxes,ietype
     +       ,maxtyp,rcut,hyperratio,barrevev,ereverseuse,isame1
     +       ,irotalign,itermax,gfac,transcrit,itranscrit,gcrit,dvcrit
     +       ,drcrit,ipr,imodedc,ifundcuse)
      else if (imode.eq.2) then
         call descent_check(natom,nmove,1,xwork,xmin,ityp,taxes,ietype
     +       ,maxtyp,rcut,hyperratio,barrevev,ereverseuse,isame1
     +       ,irotalign,itermax,gfac,transcrit,itranscrit,gcrit,dvcrit
     +       ,drcrit,ipr,imodedc,ifundcuse)
      endif

      filnam='dimerrc1.'//chri(itag)//'.dat'
      if(ipr.ge.10) call storefile(natom,xwork,xwork,ityp,taxes,filnam)
       
      if (isame1.eq.0) then
         call vecmov(xwork,xwork2,3*natom)
      endif
      do j=1,3*natom
         xwork(j)=x(j)-0.1d0*rdimer(j)*isignd
      enddo
      imodedc=2
      ifundcuse=ifundc
      if (ifundcuse.eq.5) ifundcuse=0
      if (imode.eq.0) then
         call descent_check(natom,nmove,1,xwork,x1,ityp,taxes,ietype
     +       ,maxtyp,rcut,hyperratio,barrevev,ereverseuse,isame2
     +       ,irotalign,itermax,gfac,transcrit,itranscrit,gcrit,dvcrit
     +       ,drcrit,ipr,imodedc,ifundcuse)
      else if (imode.eq.2) then
         call descent_check(natom,nmove,1,xwork,xmin,ityp,taxes,ietype
     +       ,maxtyp,rcut,hyperratio,barrevev,ereverseuse,isame2
     +       ,irotalign,itermax,gfac,transcrit,itranscrit,gcrit,dvcrit
     +       ,drcrit,ipr,imodedc,ifundcuse)
      endif

      ngradrc=ngradcall-ngradrchold

      filnam='dimerrc2.'//chri(itag)//'.dat'
      if(ipr.ge.10) call storefile(natom,xwork,xwork,ityp,taxes,filnam)

      if (isame2.eq.0) then
         call vecmov(xwork,xwork2,3*natom)
      endif

      if(rank_f.eq.0) write(6,666) isame1,isame2,ngradrc
 666  format ('dimersearch: after rollcheck: isame1=',i5,', isame2=',i5,
     +    ',ngrads=',i10)

      if (isame1.eq.0.and.isame2.eq.0.and.imode.eq.2) then
         if ((ipr.ge.4).and.(rank_f.eq.0)) write(6,*)
     +       'dimer_search: this reused saddle does not connect min'
         ierr=10
         return
      endif

      if (isame1.eq.0.and.isame2.eq.0) then
         if ((ipr.ge.4).and.(rank_f.eq.0)) write(6,*)
     +       'dimer_search: roll check does not connect to min;'
         if (imode.eq.2) then
            ierr=10
            if(rank_f.eq.0)
     +       write(6,*) 'dimer_search: going to next reused saddle'
            return
         endif
         if ((ipr.ge.4).and.(rank_f.eq.0)) write(6,*)
     +       'dimer_search: redoing this dimer from scratch'
         goto 111
      else if (isame1.eq.1.and.isame2.eq.1) then
         if ((ipr.ge.4).and.(rank_f.eq.0)) write(6,*)
     +       'dimer_search: both roll checks connect min 1;'
         if (imode.eq.2) then
            ierr=10
            if(rank_f.eq.0)
     +       write(6,*) 'dimer_search: going to next reused saddle'
            return
         endif
         if ((ipr.ge.4).and.(rank_f.eq.0)) write(6,*)
     +       'dimer_search: redoing this dimer from scratch'
         goto 111
      else
         call vecmov(xwork2,x1,3*natom)
      endif
      
      return
      end
      
c**********************************************************************
      subroutine dimer_min_torque(natom,nmove,x,rdimer,drdimer,eig,ityp
     +    ,taxes,ietype,maxtyp,rcut,ecenter,gradcenter,hyperratio)

      implicit real*8(a-h,o-z)
c      parameter(mxatom=10000)
      include 'parameters.h'

      dimension ityp(natom)
      dimension rcut(maxtyp,maxtyp)

      dimension x(3*natom),rdimer(3*natom)
      dimension gplus(3*mxatom),gminus(3*mxatom)
      dimension xtemp(3*mxatom),grad(3*mxatom)
      dimension taxes(3),gradcenter(3*mxatom)

c parameters that might need to be passed later
      gcrit=1.0d-4
      itermin=2
      dimergfac=1.0d0

      dimergfacorig=dimergfac

      iconverge=0
      eigprev=0.0d0
      i=0
 100  continue
      call vecrenorm(rdimer,3*natom)

      do j=1,3*natom
         xtemp(j)=x(j)
      enddo
      
      do j=1,3*nmove
         xtemp(j)=x(j)+drdimer*rdimer(j)
      enddo
      call gcalc(natom,nmove,xtemp,ityp, taxes,
     x    ietype, maxtyp,rcut, eplus,gplus, hyperratio)

      do j=1,3*nmove
         xtemp(j)=x(j)-drdimer*rdimer(j)
      enddo
      call gcalc(natom,nmove,xtemp,ityp, taxes,
     x    ietype, maxtyp,rcut, eminus,gminus, hyperratio)
      
      eig=(eplus+eminus-2.0d0*ecenter)/drdimer**2.0d0
      
      do j=1,3*nmove
         grad(j)=drdimer*(gplus(j)-gminus(j))/drdimer**2.0d0
      enddo

      if (i.gt.1.and.(eig-eigprev).gt.1.d-6) then ! eigenvalue increasing
         dimergfac=dimergfac/2.0d0
c         write(6,62) 'dimer_min_torque: eig increasing:',i,eig,eigprev
c 62      format(a35,i6,2d30.20)
         if (dimergfac/dimergfacorig.lt.1.d-3) then
            iconverge=-1
         endif
      endif

      gdotdr=vecdot(grad,rdimer,3*nmove)
      call vecypax(grad,-gdotdr,rdimer,3*nmove) ! get purely rotational grad
      call vecypax(rdimer,-dimergfac,grad,3*nmove) ! rotate rdimer
      gradmag=sqrt(vecdot(grad,grad,3*nmove))
      
c      write(6,63) eplus,eminus,eig,gradmag
c 63   format('eplus,eminus,eig,gradmag=',2p3d27.17,0p1d12.5)
c      if (mod(i-1,1).eq.0) write(6,61) i,gdotdr,gradmag
c 61   format('dimer_min_torque: i,gdotdr,gradmag=',i5,1pd10.2,1pd10.2)

      eigprev=eig
      i=i+1
      if ((gradmag.gt.gcrit.or.i.lt.itermin).and.iconverge.ge.0) goto
     +    100
c      write(6,*) 'dimer_min_torque: returning, iter,iconverge=',i,iconverge

      return
      end
c**********************************************************************
      subroutine dimer_initialize(natom,nmove,x1,x,idimerdof,ndimerdof,
     +    rdimer)

      use mod_mpi
      implicit real*8(a-h,o-z)
      dimension rdimer(3*natom),x1(3*natom),x(3*natom)

      if (ndimerdof.gt.nmove*3) then
         if (ipr.ge.4) then
            if(rank_f.eq.0) write(6,*)
     +          'dimer_init: ndimerdof gt nmove*3, setting to nmove*3',
     +          ndimerdof,nmove*3
            ndimerdof=nmove*3
         endif
      endif

      if (idimerdof.eq.1) then
         do j=1,ndimerdof
           if(rank_f.eq.0) rdimer(j)=prngen(0)-0.5d0
         enddo
         call MPI_BCAST(rdimer,ndimerdof,MPI_REAL8,0,force_comm,ier)
         do j=ndimerdof+1,3*natom
            rdimer(j)=0.0d0
         enddo
         call vecrenorm(rdimer,3*natom)
      else if (idimerdof.eq.2) then
         do j=1,nmove
            ii=(j-1)*3+1
            if (x1(ii).ne.x(ii)) then
               do k=1,3
                 if(rank_f.eq.0) rdimer(ii-1+k)=prngen(0)-0.5d0
               enddo
            endif
         enddo
         call MPI_BCAST(rdimer,3*natom,MPI_REAL8,0,force_comm,ier)
         call vecrenorm(rdimer,3*natom)
      else if (idimerdof.eq.3) then
         do j=1,3*nmove
           if(rank_f.eq.0) rdimer(j)=prngen(0)-0.5d0
         enddo
         call MPI_BCAST(rdimer,3*natom,MPI_REAL8,0,force_comm,ier)
         call vecrenorm(rdimer,3*natom)
      else
        if(rank_f.eq.0)
     +   write(6,*) 'dimer_init: idimerdof not recognized'
        stop 'dimer_init: idimerdof not recognized'
      endif

      return
      end
c**********************************************************************
      subroutine dimer_ranx(natom,nmove,x1,x,ipr,idimerdof
     +    ,ndimerdof,ndimercoord,idimeroverunder,taxes,icoord,rcoordcrit
     +    ,rranx)

      use mod_mpi
      implicit real*8(a-h,o-z)
      dimension x(3*natom),x1(3*natom)
c      parameter(mxatom=10000)
      include 'parameters.h'
      dimension icoord(nmove)
      dimension taxes(3)

      if (ndimerdof.gt.nmove*3) then
         if (ipr.ge.4) then
            if(rank_f.eq.0) write(6,*)
     +          'dimer_ranx: ndimerdof gt nmove*3, setting to nmove*3'
            ndimerdof=nmove*3
         endif
      endif

      if (ipr.ge.6.and.idimerdof.eq.1)then
        if(rank_f.eq.0)
     +    write(6,*)'dimer_ranx: randomizing along ',ndimerdof
     +    ,' degrees of freedom'
      endif
     
      if (idimerdof.eq.1) then
         xdist=0.0d0
         do j=1,ndimerdof
          if(rank_f.eq.0)then
            x(j)=x1(j)+(prngen(0)-0.5d0)*0.75d0
            xdist=xdist+(x(j)-x1(j))**2
          endif
         enddo
         call MPI_BCAST(xdist,ndimerdof,MPI_REAL8,0,force_comm,ier)
         do j=ndimerdof+1,3*natom
            x(j)=x1(j)
         enddo
         xdist=sqrt(xdist)
         if(rank_f.eq.0)then
          if(ipr.ge.6) write(6,*)'dimer_ranx: starting point is ',xdist
     +       ,' from minimum'
         endif
      else if (idimerdof.eq.2) then
         if(rank_f.eq.0) write(6,*) 'dimer_ranx: idimerdof=',idimerdof
         call vecmov(x1,x,3*natom)
         nmaxcoord=ndimercoord
         iunder=idimeroverunder
         nunder=0
         do i=1,nmove
c            write(6,*) 'atom: ',i,', coord: ',icoord(i)
            if (icoord(i).gt.nmaxcoord.and.iunder.eq.0) nunder=nunder+1 
            if (icoord(i).lt.nmaxcoord.and.iunder.eq.1) nunder=nunder+1 
            if (iunder.eq.2.and.(icoord(i).gt.nmaxcoord.or.
     +          icoord(i).lt.nmaxcoord)) nunder=nunder+1
         enddo
         if(rank_f.eq.0)then
           if(ipr.ge.6) write(6,*)'dimer_ranx: nunder=',nunder
         endif
c         nrand=int(prngen(0)*(float(nunder)-1)+1)
         if(rank_f.eq.0) nrand=prngen(0)*float(nunder)+1.0d0
         call MPI_BCAST(nrand,1,MPI_REAL8,0,force_comm,ier)
         if(rank_f.eq.0)then
           if(ipr.ge.6) write(6,*) 'nrand, nunder: ',nrand,nunder
         endif
         nunder=0
         iselect=0
         do i=1,nmove
            if (icoord(i).gt.nmaxcoord.and.iunder.eq.0) nunder=nunder+1
            if (icoord(i).lt.nmaxcoord.and.iunder.eq.1) nunder=nunder+1 
            if (iunder.eq.2.and.(icoord(i).gt.nmaxcoord.or.
     +          icoord(i).lt.nmaxcoord)) nunder=nunder+1
            if (nunder.eq.nrand.and.iselect.eq.0) iselect=i
         enddo
         if(rank_f.eq.0)then
           if(ipr.ge.6) write(6,*) 'dimer_ranx: iselect=',iselect
         endif
         if (iselect.eq.0) then
           if(rank_f.eq.0) write(6,*)
     +          'STOP: no atom met dimer displacement criterion'
           stop 'no atom met dimer displacement criterion'
         endif
         do i=1,nmove
            rij=0.0d0
            do j=1,3
               rr=x1((iselect-1)*3+j)-x1((i-1)*3+j)
               if (abs(rr).gt.taxes(j)/2.0d0) rr=rr-sign(taxes(j),rr)
               rij=rij+rr**2
            enddo
            if (sqrt(rij).lt.rcoordcrit) then
               do j=1,3
                  ii=(i-1)*3+j
c                  x(ii)=x1(ii)+(prngen(0)-0.5d0)*rranx*2.0d0
                  if(rank_f.eq.0) x(ii)=x1(ii)+rranx*gasdev(0)
               enddo
               if((ipr.ge.6).and.(rank_f.eq.0))
     +             write(6,*) 'dimer_ranx: atom ',
     +             i,' displaced'
            endif
         enddo
         call MPI_BCAST(x,nmove*3,MPI_REAL8,0,force_comm,ier)

      else
        if(rank_f.eq.0) write(6,*)'dimer_ranx: idimerdof not recognized'
        stop 'dimer_ranx: idimerdof not recognized'
      endif

      return
      end


c**********************************************************************
      subroutine dimer_min_torque_2(imode,natom,nmove,x,rdimer,drdimer
     +    ,eig,ityp,taxes,ietype,maxtyp,rcut,ecenter,gradcenter
     +    ,hyperratio)

c uses modified newton method described by Henkelman and Jonsson
c instead of steepest descent of dimer_min_torque
      use mod_mpi
c      implicit real*8(a-h,o-z)
      implicit none
c      save
c      integer mxatom,
      integer imode
c      parameter(mxatom=10000)
      include 'parameters.h'


      integer natom,nmove,maxtyp,ityp,ietype
      real*8 eig,ecenter,hyperratio,drdimer
      real*8 x(3*natom),rdimer(3*natom)
      real*8 taxes(3),rcut(maxtyp,maxtyp)
      real*8 vecdot,gradcenter(3*natom)
      dimension ityp(natom)

      integer i,j
      real*8 gcrit,itermin,itermax,PI
      real*8 dimergfac,dimergfacorig,gradmag
      real*8 NdR,FD1,FD2,GamN,U0,U1,U2,Cth,CN
      real*8 dTh,Th,FNr,FNrp1,FNrp2,FNrp,FNp
      real*8 R(3*mxatom),N(3*mxatom),dF(3*mxatom)
      real*8 F0(3*mxatom),FN(3*mxatom)
      real*8 R1(3*mxatom),R2(3*mxatom)
      real*8 F1(3*mxatom),F2(3*mxatom)
      real*8 FN1(3*mxatom),FN2(3*mxatom)
      real*8 FNp1(3*mxatom),FNp2(3*mxatom)
      real*8 GN(3*mxatom),GNu(3*mxatom)
      
      
c parameters that might need to be passed later
      gcrit=1.0d-2
      if (imode.eq.3) gcrit=1.0d-3
      itermin=1
      itermax=50
      dimergfac=1.0d0
      dimergfacorig=dimergfac
      U0=ecenter
      PI=3.1415926535d0

      NdR=drdimer
      GamN=0
      dTh=1.0d-3
      i=0
      
      do j=1,3*natom
         R(j)=x(j)
         N(j)=rdimer(j)
         F0(j)=-gradcenter(j)
      enddo
      
 100  continue
c      write(*,*) 'Iteration: ',i+1
      do j=nmove*3+1,natom*3
        N(j)=0.0d0
      enddo
      call vecrenorm(N,3*nmove)
      
      do j=1,3*natom
         R1(j)=R(j)+NdR*N(j) ! R1
         R2(j)=R(j)-NdR*N(j) ! R2
      enddo
      
      call gcalc(natom,nmove,R1,ityp, taxes,
     x    ietype, maxtyp,rcut, U1,F1, hyperratio)

      do j=1,3*nmove
        F1(j)=-F1(j)
      enddo
      do j=nmove*3+1,natom*3
        F1(j)=0.0d0
      enddo

! Central Difference
!      call gcalc(natom,nmove,R2,ityp, taxes,
!     x    ietype, maxtyp,rcut, U2,F2, hyperratio)
!      do j=1,3*nmove
!        F2(j)=-F2(j)
!      enddo
c      write(*,*) 'U2 (from gcalc): ',U2
! Forward Difference
      do j=1,3*nmove
        dF(j)=F1(j)-F0(j)
      enddo
      U2=U0-2.0d0*U1+NdR*vecdot(dF,N,3*nmove)
      do j=1,3*nmove
        F2(j)=F0(j)*2.0d0-F1(j)
      enddo
c      write(*,*) 'U2 (from frdif): ',U2

      do j=nmove*3+1,natom*3
        F2(j)=0.0d0
      enddo

      FD1=vecdot(F1,N,3*nmove) ! F1.N
      FD2=vecdot(F2,N,3*nmove) ! F2.N      
      do j=1,3*nmove
         FN1(j)=F1(j)-FD1*N(j)
         FN2(j)=F2(j)-FD2*N(j)
         FN(j)=(FN1(j)-FN2(j))/NdR
      enddo

      FNr=sqrt(vecdot(FN,FN,3*nmove))
      U0=(U1+U2)/2.0d0+(FD1-FD2)*(NdR/4.0d0)
      CN=(FD2-FD1)/(2.0d0*NdR)
      
      ! This is a steepest descent using the modified Newton's rotation
      do j=1,3*nmove
        GN(j)=FN(j)
        GNu(j)=GN(j)
      enddo
      do j=nmove*3+1,natom*3
        GN(j)=0.0d0
        GNu(j)=0.0d0
      enddo
      call vecrenorm(GNu,3*nmove)
      
      FNp=vecdot(FN,GNu,3*nmove)
      do j=1,3*nmove
        FNp1(j)=GNu(j)*FNp
      enddo
      FNrp1=vecdot(FNp1,GNu,3*nmove)

c      write(*,*) 'FNrp1',FNrp1,' FNr',FNr,' U0',U0,' CN',CN
      
      call Rotate(N,GNu,dTh,3*nmove)
      call vecrenorm(N,3*nmove)
      call vecrenorm(GNu,3*nmove)

      do j=1,3*natom
         R1(j)=R(j)+NdR*N(j) ! R1
         R2(j)=R(j)-NdR*N(j) ! R2
      enddo
      
      call gcalc(natom,nmove,R1,ityp, taxes,
     x    ietype, maxtyp,rcut, U1,F1, hyperratio)

      do j=1,3*nmove
        F1(j)=-F1(j)
      enddo

      do j=nmove*3+1,natom*3
        F1(j)=0.0d0
      enddo
! Central Difference
!      call gcalc(natom,nmove,R2,ityp, taxes,
!     x    ietype, maxtyp,rcut, U2,F2, hyperratio)
!      do j=1,3*nmove
!        F2(j)=-F2(j)
!      enddo

! Forward Difference
      do j=1,3*nmove
        dF(j)=F1(j)-F0(j)
      enddo

      U2=U0-2.0d0*U1+NdR*vecdot(dF,N,3*nmove)
      do j=1,3*nmove
        F2(j)=F0(j)*2.0d0-F1(j)
      enddo

      do j=nmove*3+1,natom*3
        F2(j)=0.0d0
      enddo

      FD1=vecdot(F1,N,3*nmove) ! F1.N
      FD2=vecdot(F2,N,3*nmove) ! F2.N
      do j=1,3*nmove
         FN1(j)=F1(j)-FD1*N(j)
         FN2(j)=F2(j)-FD2*N(j)
         FN(j)=(FN1(j)-FN2(j))/NdR
      enddo

      FNr=sqrt(vecdot(FN,FN,3*nmove))
      U0=(U1+U2)/2.0d0+(FD1-FD2)*(NdR/4.0d0)
      CN=(FD2-FD1)/(2.0d0*NdR)

      FNp=vecdot(FN,GNu,3*nmove)
      do j=1,3*nmove
        FNp2(j)=GNu(j)*FNp
      enddo
      FNrp2=vecdot(FNp2,GNu,3*nmove)
      Cth=(FNrp1-FNrp2)/dTh
      FNrP=(FNrp1+FNrp2)/2.0d0
      Th=atan((FNrp/Cth)*2.0d0)/2.0d0-dTh/2.0d0
      if(Cth.lt.0) Th=Th+PI/2.0d0
      
c      write(*,*) 'FNrp2',FNrp2,' FNr',FNr,' U0',U0,' CN',CN
c      write(*,*) 'Cth',Cth,' FNrP',FNrP,' Th',Th
c      write(*,*) 'Itr ',i+1,' Cth',Cth,' FNrP',FNrP,' Th',Th

      call Rotate(N,GNu,Th,3*nmove)
      call vecrenorm(N,3*nmove)
      call vecrenorm(GNu,3*nmove)

c      do j=1,3*natom
c         R1(j)=R(j)+NdR*N(j) ! R1
c         R2(j)=R(j)-NdR*N(j) ! R2
c      enddo

      i=i+1
c      write(*,*) 'gradmag',gradmag
      gradmag=FNrp1
c      if ((gradmag.gt.gcrit.or.i.lt.itermin).and.iconverge.ge.0) goto
      if ((gradmag.gt.gcrit.or.i.lt.itermin).and.(i.lt.itermax)) goto
     +    100
    
c return values
      eig=CN
!      ecenter=U0
      do j=1,3*natom
        rdimer(j)=N(j)
      enddo
     
      return
      end
      
!**********************************************************************
      subroutine Rotate(V1,V2,Th,num)
! Rotate the vectors V1 and V2
      integer num
      real*8 V1(num),V2(num),Th
      real*8 V1tmp(num),cth,sth
      cth=cos(Th)
      sth=sin(Th)
      do j=1,num
        V1tmp(j)=V1(j)
        V1(j)=V1(j)*cth+V2(j)*sth
        V2(j)=V2(j)*cth-V1tmp(j)*sth
      enddo
      end !subroutine Rotate


c**********************************************************************
      subroutine shiftsaddle(natom,nmove,rdimer,xyzsad,xyzmin1,xyzmin2
     +    ,transcrit,ierr)
      
      implicit real*8(a-h,o-z)

      dimension xyzsad(3*natom),xyzmin1(3*natom),xyzmin2(3*natom)
     +    ,rdimer(3*natom)
      
      ierr=0

      do i=1,nmove
         rij=0.0d0
         rik=0.0d0
         do j=1,3
            rij=rij+(xyzmin2((i-1)*3+j)-xyzmin1((i-1)*3+j))**2
            rik=rik+(xyzsad((i-1)*3+j)-xyzmin1((i-1)*3+j))**2
         enddo
         rij=sqrt(rij)
         rik=sqrt(rik)
         if (rij.ge.transcrit) then
            do j=1,3
               xyzsad((i-1)*3+j)=xyzmin2((i-1)*3+j)
               rdimer((i-1)*3+j)=0.0d0
            enddo
         endif
      enddo

      do i=1,nmove
         do k=1,nmove
            rik=0.0d0
            if (k.ne.i) then
               do j=1,3
                  rik=rik+(xyzsad((i-1)*3+j)-xyzsad((k-1)*3+j))**2 ! xxxx no pbcs here, need to add
               enddo
               if (sqrt(rik).le.1.05d0) then
                  ierr=1
                  write(6,*) 'shiftsaddle: ierr=1, atoms ',i,k,sqrt(rik)
                  return
               endif
            endif
         enddo
      enddo
      
      return
      end

         
c**********************************************************************
      subroutine dimer_ranx_hole(natom,nmove,x,xmin,xyzprev,ipr
     +    ,idimerdof,ndimerdof,ndimercoord,idimeroverunder,taxes
     +    ,rranx,transcrit,rcoordcrit)

      use mod_mpi
      implicit real*8(a-h,o-z)
      dimension xmin(3*natom),x(3*natom),xyzprev(3*natom)
c      parameter(mxatom=10000)
      include 'parameters.h'
      dimension iindex(100)
      dimension taxes(3)

      inumber=0
      call vecmov(xmin,x,3*natom)
      do i=1,nmove
         rij=0.0d0
         do j=1,3
            rij=rij+(xmin((i-1)*3+j)-xyzprev((i-1)*3+j))**2
         enddo
         rij=sqrt(rij)
         if (rij.ge.transcrit) then
            inumber=inumber+1
            if (inumber.gt.100) then
               write(6,*) 'STOP: inumber too big'
               stop 'inumber too big'
            endif
            iindex(inumber)=i
         endif
      enddo

      if(rank_f.eq.0) irand=prngen(0)*float(inumber)+1.0d0
      call MPI_BCAST(irand,1,MPI_REAL8,0,force_comm,ier)
      iselect=iindex(irand)
      if (ipr.ge.4) write(6,*) 'dimer_ranx_hole: iselect=',iselect,
     +      rcoordcrit,inumber,irand

      if (iselect.eq.0) then
         write(6,*) 'STOP: no atom met dimer displacement criterion'
         stop 'no atom met dimer displacement criterion'
      endif
      do i=1,nmove
         rij=0.0d0
         do j=1,3
            rr=xmin((iselect-1)*3+j)-xmin((i-1)*3+j)
            if (abs(rr).gt.taxes(j)/2.0d0) rr=rr-sign(taxes(j),rr)
            rij=rij+rr**2
         enddo
         if (sqrt(rij).lt.rcoordcrit) then
            do j=1,3
               ii=(i-1)*3+j
               if(rank_f.eq.0) x(ii)=xmin(ii)+rranx*gasdev(0)
            enddo
            if (ipr.ge.4) write(6,*) 'dimer_ranx_hole: atom ',
     +          i,' displaced'
         endif
      enddo
      call MPI_BCAST(x,nmove*3,MPI_REAL8,0,force_comm,ier)

      return
      end

