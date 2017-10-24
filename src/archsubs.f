c machine dependent subroutines

      subroutine flushtad(iunit)
      implicit none
      integer iunit

      call flush(iunit)

      return
      end

      real*4 function etimetad(tarray)

      implicit none
      real*4 tarray(2),etime

      etimetad=etime(tarray)

      return
      end



      
