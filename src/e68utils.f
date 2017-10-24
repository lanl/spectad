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


*-----------------------------------------------------------------------
      subroutine dummye68(ntype,maxtyp,rcut,amu,azero)
*-----------------------------------------------------------------------

      implicit real*8(a-h,o-z)

      dimension rcut(maxtyp,maxtyp),amu(maxtyp)

      stop 'TB not implemented in this version of TAD'     
      
      return
      end
      
