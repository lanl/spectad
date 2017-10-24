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
      subroutine prmst70(natom,ntype,itype,maxtyp,rcut,amu,azero,
     +    element)
*-----------------------------------------------------------------------

      implicit real*8(a-h,o-z)

      dimension rcut(maxtyp,maxtyp),amu(maxtyp)
      character*80 element(maxtyp)
      dimension itype(natom),azero(maxtyp)

      stop 'etype=70 not implemented in this version of TAD'     
      
      return
      end
      
c**********************************************************************
      subroutine g70(natom,nmove,x,y,z,itype,tax,tay,taz,ietype,maxtyp
     +    ,rcut,energy,grad,hyperratio)
c      implicit real*8 (a-h,o-z)

      dimension grad(3*natom)
      dimension x(natom),y(natom),z(natom),itype(natom)
      dimension rcut(maxtyp,maxtyp)

      stop 'etype=70 not implemented in this version of TAD'     
      
      return
      end
