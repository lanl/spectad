       !---------------------------------------------------------------------
       !
       !  Dummy Module containing Lammps-Related Info/Subroutines
       !        
       !---------------------------------------------------------------------

       MODULE mod_lammps

          use mod_mpi
          logical :: doblocks = .false.

       CONTAINS
       
!       subroutine lmp_setup(natom,nmove,ntype,maxtyp,potnam,ietype,
!     +     xyz,taxes,itype,rcut,amu,azero,iseed,templow,temphigh,dt)
       subroutine lmp_setup(natom,nmove,ntype,maxtyp,potnam,ietype
     +     ,xyz,taxes,itype,rcut,amu,templow,temphigh
     +     ,lmp_end,dt)

         implicit real*8(a-h,o-z)
         include 'parameters.h'
         character*80 potnam(maxtyp)
         dimension taxes(3)
         dimension itype(mxatom)
         dimension xyz(3*mxatom)
         dimension pxyz(3*mxatom)
         logical lmp_end
         dimension rcut(maxtyp,maxtyp),amu(maxtyp),azero(maxtyp)
         write(*,*) ' ERROR.. This setup is not currently configured!!'
         STOP ' ERROR.. This setup is not currently configured!!'
        return
       end subroutine lmp_setup

      subroutine g100(natom,nmove,x,y,z,itype,ietype,maxtyp,energy,grad
     +     ,taxes)

         implicit real*8(a-h,o-z)
         save
         include 'parameters.h'
         dimension x(mxatom),y(mxatom),z(mxatom),itype(natom)
         dimension taxes(3)
         dimension grad(3,natom)
         write(*,*) ' ERROR.. g100 is not currently configured!!'
         STOP ' ERROR.. g100 is not currently configured!!'
        return
       end  subroutine g100

      subroutine g100spawn(natom,nmove,xyz,ietype,itype,energy,grad)

         implicit real*8(a-h,o-z)
         include 'parameters.h'
         dimension xyz(mxatom*3),itype(natom)
         dimension grad(3,natom)
         dimension taxes(3)
         write(*,*) ' ERROR.. g100spawn is not currently configured!!'
         STOP ' ERROR.. g100spawn is not currently configured!!'
        return
       end  subroutine g100spawn

      subroutine unfix_qeq()
         implicit real*8(a-h,o-z)
         include 'parameters.h'
         write(*,*) ' ERROR.. unfix_qeq is not currently configured!!'
         STOP ' ERROR.. unfix_qeq is not currently configured!!'
        return
      end  subroutine unfix_qeq

      subroutine refix_qeq()
         implicit real*8(a-h,o-z)
         include 'parameters.h'
         write(*,*) ' ERROR.. refix_qeq is not currently configured!!'
         STOP ' ERROR.. refix_qeq is not currently configured!!'
        return
       end  subroutine refix_qeq

      subroutine lmp_minimize(natom,nmove,xyz,itype,ietype,energy,taxes)
         implicit real*8(a-h,o-z)
         include 'parameters.h'
         !real (C_double), pointer :: energyl => NULL()
         dimension xyz(natom*3),itype(natom)
         dimension taxes(3)
         write(*,*) ' ERROR.. refix_qeq is not currently configured!!'
         STOP ' ERROR.. refix_qeq is not currently configured!!'
        return
       end  subroutine lmp_minimize

       subroutine lammps_command(lmp,str)
         implicit real*8(a-h,o-z)
         character*80 str
         include 'parameters.h'
         write(*,*) ' ERROR.. lammps_command is not configured!!'
         STOP ' ERROR.. lammps_command is not currently configured!!'
        return
       end  subroutine lammps_command

      subroutine vlmptopxyz(natom,nmove,pxyz,itype,amu,ntype,maxtyp)
         implicit real*8(a-h,o-z)
         include 'parameters.h'
         write(*,*) ' ERROR.. vlmptopxyz is not currently configured!!'
         STOP ' ERROR.. vlmptopxyz is not currently configured!!'
        return
        return
       end  subroutine vlmptopxyz

      subroutine pxyztovlmp(natom,nmove,pxyz,itype,amu,ntype,maxtyp)
         implicit real*8(a-h,o-z)
         include 'parameters.h'
         write(*,*) ' ERROR.. pxyztovlmp is not currently configured!!'
         STOP ' ERROR.. pxyztovlmp is not currently configured!!'
        return
       end  subroutine pxyztovlmp
       
      subroutine mdblock100(natom,nmove,xyz,pxyz,itype,maxtyp
     x   ,energy,grad,taxes,ietype,amu,temp,dt
     x   ,thermrate,nmd,iseed,ntype)
         implicit real*8(a-h,o-z)
         include 'parameters.h'
         dimension taxes(3)
         dimension itype(mxatom)
         dimension xyz(3*mxatom)
         dimension pxyz(3*mxatom)
         dimension grad(3,natom)
         dimension amu(maxtyp)
         write(*,*) ' ERROR.. mdblock100 is not currently configured!!'
         STOP ' ERROR.. mdblock100 is not currently configured!!'
        return
       end  subroutine mdblock100


       END MODULE mod_lammps
       
       
