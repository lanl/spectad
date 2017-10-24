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


c**********************************************************************
      subroutine refine_transition(natom,nmove,nsnap,x1,x2,xsnap
     +    ,xyzneighs,nneighs,barrevev,ereverse,ityp,taxes,ietype,maxtyp
     +    ,rcut,hyperratio,irotalign,itermax,gfac,transcrit,itranscrit
     +    ,intermediate,gcrit,dvcrit,drcrit,ipr,imode,ileft,iright
     +    ,ifundc,rank_f)

      implicit real*8(a-h,o-z)

c      parameter (mxatom=10000)
      include 'parameters.h'
      
      dimension x1(3*natom),x2(3*natom),xsnap(3*natom*nsnap)
      dimension xyzneighs(3*natom*nneighs),barrevev(nneighs-1)
      dimension rcut(maxtyp,maxtyp),ityp(natom)
      dimension x12(3*mxatom*2)
      dimension taxes(3)
      character*180 filnam
      character*10 chri
      integer rank_f

      do j=1,3*natom            ! setup x12 for passing to descent_check
         x12(j)=x1(j)
         x12(j+3*natom)=x2(j)
      enddo

      intermediate=0
      ileft=1
      iright=nsnap

      isnap=(iright-ileft)/2+ileft
      if((ipr.ge.4).and.(rank_f.eq.0)) 
     +    write(6,*) 'RT: new isnap: ',ileft,isnap,iright

 200  continue                  !loop for checking new snap
      indexr=(isnap-1)*3*natom+1
      ereverseuse=0.0d0         ! don't do reverse barrier stuff as here comparing to x12
      call descent_check(natom,nmove,2,xsnap(indexr),x12,
     +    ityp,taxes,ietype,maxtyp,rcut,hyperratio,barrevev,ereverseuse
     +    ,isame,irotalign,itermax,gfac,transcrit,itranscrit,gcrit
     +    ,dvcrit,drcrit,ipr,imode,ifundc)
      if (isame.ne.0) goto 100

      if (ereverse.gt.0.0d0) then
         indexr=(isnap-1)*3*natom+1
         call descent_check(natom,nmove,nneighs,xsnap(indexr),xyzneighs
     +       ,ityp,taxes,ietype,maxtyp,rcut,hyperratio,barrevev,ereverse
     +       ,isamerev,irotalign,itermax,gfac,transcrit,itranscrit,gcrit
     +       ,dvcrit,drcrit,ipr,imode,ifundc)
         if (isamerev.gt.1) then
            if (barrevev(isamerev-1).lt.ereverse) then
               if ((ipr.ge.4).and.(rank_f.eq.0)) write(6,*)
     +             'refinetransition: intermediate min found '
     +             ,'but reverse barrier (',barrevev(isamerev-1)
     +             ,') less than ereverse=',ereverse
               isame=1
               goto 100
            endif
         endif
      endif

      if((ipr.ge.4).and.(rank_f.eq.0))
     +  write(6,*) 'refinetransition note: snap=',isnap
     +    ,' is in neither basin x1 nor x2!'
      if (ipr.ge.10) then
         filnam='snap.'//chri(isnap)//'.dat'
         call storefile(natom,xsnap(indexr),xsnap(indexr),ityp,taxes,
     +       filnam)
         filnam='snap.'//chri(isnap)//'.1'//'.dat'
         call storefile(natom,x12(1),x12(1),ityp,taxes,
     +       filnam)
         filnam='snap.'//chri(isnap)//'.2'//'.dat'
         call storefile(natom,x12(3*natom+1),x12(3*natom+1),ityp,taxes,
     +       filnam)
      endif

      intermediate=1
      do j=1,3*natom
         x2(j)=xsnap((isnap-1)*3*natom+j)
         x12(j+3*natom)=x2(j)
      enddo
      iright=isnap

 100  continue
      if (isame.eq.2) then
         if ((ipr.ge.4).and.(rank_f.eq.0))
     +     write(6,*) 'refinetransition: snap ',isnap,' is in basin 2'
         iright=isnap
      else if (isame.eq.1) then
         if ((ipr.ge.4).and.(rank_f.eq.0))
     +     write(6,*) 'refinetransition: snap ',isnap,' is in basin 1'
         ileft=isnap
      endif
      isnap=(iright-ileft)/2+ileft
      if ((ipr.ge.4).and.(rank_f.eq.0))
     +  write(6,*) 'refinetransition: new isnap: ',ileft,isnap,iright
      if (isnap.eq.ileft.or.isnap.eq.iright) return
      goto 200      
      
      return
      end
c**********************************************************************
