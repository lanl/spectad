      subroutine putclsxyz(natom,xyz,iunit)

      implicit real*8(a-h,o-z)
      dimension xyz(3,natom)

      write(iunit,50) natom
      write(iunit,*) 0.00d0

      do i=1,natom
         write(iunit,50) i,(xyz(j,i),j=1,3)
      enddo

 50   format(i10,3f20.8)

      return
      end
