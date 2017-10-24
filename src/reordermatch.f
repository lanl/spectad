      subroutine reordermatch(natom,nmove,xyz,itype,xyz2,crit,
     x                                              imatch,kreorder)
c see if reording the atoms in xyz allows a match with those in xyz2 - afv 4/15/03
c If so, return imatch=1, and korder array filled with this transformation
c We assume that crit is small enough that there would never be two atoms "matching"
c one in the other cluster.
c For now (for simplicity), we do not make use of nmove.  All atoms are treated.
      use mod_mpi

      implicit real*8(a-h,o-z)
      dimension xyz(3,natom)
      dimension itype(natom)
      dimension xyz2(3,natom)
      dimension kreorder(natom)
c quick scratch array - later put this into work common 
      dimension iused(10000)

      if(natom.gt.10000) stop 'redimension iused in reordermatch'

c zero out arrays

      do i=1,natom
        kreorder(i)=0
        iused(i)=0
      end do

c match them up until failure occurs

      do 100 i=1,natom
      if(kreorder(i).ne.0) go to 100

      do 80 jj=1,natom
        if(jj.eq.1) then     ! try i=j first, since most likely to match
          j=i
        else if(jj.eq.i) then
          j=1
        else
          j=jj
        end if
      if(iused(j).ne.0) go to 80
      if(itype(j).ne.itype(i)) go to 80 ! Make sure types are the same!
      sum=0.0d0
      do jx=1,3
        sum = sum + (xyz(jx,i)-xyz2(jx,j))**2
      end do
      if(sum.lt.crit**2) then
        kreorder(i)=j
        iused(j)=i
        go to 100
      end if
   80 continue

      if(rank_f.eq.0)
     +  write(6,*) 'reordermatch: no match for this pair'
      imatch=0
      return

  100 continue

      imatch=1

c write out compact report of reordering to unit 35, as sort of a test
c reuse the iused array here, since we are done with it above

      nentry=0
      isum=0
      do i=1,nmove
        if(kreorder(i).eq.i) then
          isum=isum+1
        else
          if(isum.ne.0) then
            nentry=nentry+1
            iused(nentry)=-isum
          end if
          nentry=nentry+1
          iused(nentry)=kreorder(i)
          isum=0
        end if
      end do
      if(isum.ne.0) then
        nentry=nentry+1
        iused(nentry)=-isum
      end if
      if(rank_f.eq.0) write(6,110)
  110 format(' compact reorder array:')
      if(rank_f.eq.0) write(6,120) (iused(i),i=1,nentry)
  120 format(20i5)

c verify (for safety) that i=j for all kreorder values above nmove 

      do i=nmove+1,natom
        if(kreorder(i).ne.i)stop'moving nonmoving atoms in reordermatch'
      end do


      return

      end


      subroutine reorder(natom,nmove,xyz,kreorder,xyz2)
c fill xyz2 by reordering the atoms in xyz via the kreorder array
c (see reordermatch routine, which fills kreorder)

      implicit real*8(a-h,o-z)
      dimension xyz(3,natom)
      dimension xyz2(3,natom)
      dimension kreorder(natom)

      do i=1,natom
        do ix=1,3 
          xyz2(ix,kreorder(i)) = xyz(ix,i)
        end do
      end do

      end
