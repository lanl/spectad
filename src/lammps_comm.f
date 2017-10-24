c  SpecTAD LAMMPS MPI Wrapper from Los Alamos National Laboratory;  2015
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

      use omp_lib
      use LAMMPS
      use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
      include 'mpif.h'
      include 'parameters.h'
      type (C_ptr) :: lmp
      real (C_double), dimension(:), pointer :: masslmp => NULL()
      real (C_double), pointer :: energyl => NULL()

      integer,   dimension(:), allocatable ::  itypelmp
      real*8,    dimension(:), allocatable ::  F_lmp
      integer,   dimension(:), allocatable ::  integerbuffer
      real*8,    dimension(:), allocatable ::  mlmp
      real*8,    dimension(:), allocatable ::  xyzlmp
      real*8,    dimension(:), allocatable ::  itype
      character, dimension(:), allocatable ::  characterbuffer
      real*8 realbuffer(3+mxatom*4)

      logical lgettypes

      integer :: groupkill = 100
      integer parent, ix, itypeadd, natomlmp
      integer msgstatus(MPI_STATUS_SIZE)
      integer ibuffsize, msgsource, msgdest, ntype, minimize
      integer ndof, size, ier, rank_lmp, nprocs_lmp
      character(LEN=80) str1,str2,str3,str4,comstr,command
      character(LEN=1) strspace

      integer, PARAMETER :: lmp_tag1 = 1
      integer, PARAMETER :: lmp_tag2 = 2
      integer, PARAMETER :: lmp_tag3 = 3
      integer, PARAMETER :: lmp_tag4 = 4
      integer, PARAMETER :: lmp_tag5 = 5
      integer, PARAMETER :: lmp_tag6 = 6

      ! Initialize MPI
      call MPI_Init (ier)

      write(*,*) ' lammps_comm.exe running... '

      ! Get info about parent
      call MPI_Comm_get_parent (parent, ier)
      if (parent.eq.MPI_COMM_NULL) STOP 'No parent!'
      call MPI_Comm_remote_size(parent, size, ier)
      write(*,*) 'parent size = ',size

      ! Get info about self
      call MPI_Comm_rank (MPI_COMM_WORLD, rank_lmp, ier)
      call MPI_Comm_size (MPI_COMM_WORLD, nprocs_lmp, ier)
      write(*,*) 'rank_lmp - ',rank_lmp
     +		    ,', nprocs_lmp - ',nprocs_lmp
      call MPI_BARRIER(MPI_COMM_WORLD,ier)

      ! Parallel code below:
      !   The manager is represented as the process with rank 0 in (the remote
      !   group of) parent.  If the workers need to communicate among
      !   themselves, they can use MPI_COMM_WORLD...

c ========= Wait for message to start lammps.
c ========= .. should contain char string to use in 'lammps_open'
      if(rank_lmp.eq.0)then
        call MPI_Probe(0,lmp_tag1,parent,msgstatus,ier)
        call MPI_Get_count(msgstatus,MPI_CHARACTER,ibuffsize,ier)
        command = ''
        call MPI_Recv(command,ibuffsize,MPI_CHARACTER,
     +      0,lmp_tag1,parent,msgstatus,ier)
      endif
      call MPI_Bcast(ibuffsize,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_Bcast(command,ibuffsize,MPI_CHARACTER,0,
     +      MPI_COMM_WORLD,ier)

c ========= Startup lammps:
c =========
      if(rank_lmp.eq.0) write(*,*) command
      call lammps_open (command, MPI_COMM_WORLD, lmp)
      call lammps_file (lmp, 'lammps.input')
      if(rank_lmp.eq.0) write(*,*) 'LAMMPS STARTED -- lmp'
          
c ========= Get Masses From Lammps...
c =========
      allocate(integerbuffer(4))
      if(rank_lmp.eq.0) call MPI_Recv(integerbuffer,4,MPI_INTEGER,
     +    0,lmp_tag2,parent,msgstatus,ier)
      call MPI_Bcast(integerbuffer,4,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      ntype = integerbuffer(1)
      natom = integerbuffer(2)
      itypeadd = integerbuffer(3)
      minimize = integerbuffer(4)
      deallocate(integerbuffer)
      call lammps_extract_atom (masslmp, lmp, 'mass')
      allocate(mlmp(ntype))
      !do ix = 1, ntype
      !  amu(ix) = masslmp(ix+1)*1.0
      !  write(*,*) 'mass -- ',amu(ix)
      !enddo
      do ix = 1, ntype
        mlmp(ix) = masslmp(ix+1)*1.0
      enddo
      if (rank_lmp.eq.0) then
        call MPI_Send(mlmp,ntype,MPI_REAL8,0,lmp_tag3,parent,ier)
      endif
      deallocate(mlmp)

c ========= Minimize Geometry in Lammps ...
c =========
      if(minimize.eq.1)then
        call lammps_command (lmp, 'minimize 0.0 0.005 1000 10000')
        call lammps_gather_atoms (lmp, 'x', 3, xyzlmp)
        if (rank_lmp.eq.0) then
          call MPI_Send(xyzlmp,natom*3,MPI_REAL8,0,lmp_tag4,parent,ier)
        endif
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ier)

      natomlmp = natom

c ========= Force-call Loop here ...
c =========
      do while(.true.)

        if (rank_lmp.eq.0) then
          call MPI_Probe(0,lmp_tag5,parent,msgstatus,ier)
          call MPI_Get_count(msgstatus,MPI_REAL8,ibuffsize,ier)
          call MPI_Recv(realbuffer,ibuffsize,MPI_REAL8,
     +        0,lmp_tag5,parent,msgstatus,ier)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,ier)
        call MPI_Bcast(ibuffsize,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
        call MPI_Bcast(realbuffer,ibuffsize,MPI_REAL8,0,
     +      MPI_COMM_WORLD,ier)

        ! Now do force-call work..
        natom =  NINT(realbuffer(1))
        nmove =  NINT(realbuffer(2))
        ietype = NINT(realbuffer(3))

         if((.not.allocated(xyzlmp)).and.(natom.eq.natomlmp))then
           call lammps_gather_atoms (lmp, 'type', 1, itypelmp)
           call lammps_gather_atoms (lmp, 'x', 3, xyzlmp)
         endif
         strspace = ' '

         ! Check that atom count in lammps is correct
         lgettypes = .false.
         if(natom.ne.natomlmp)then
           lgettypes = .true.
           write(*,*) 'natom, natomlmp - ',natom, natomlmp
           if(natom.gt.natomlmp)then
             nadd = natom-natomlmp
             do i = 1,nadd
               str1 = 'create_atoms'
               str2 = ''
               write(str2,*) itypeadd
               str3 = 'single 0.0 0.0 0.0'
               comstr =   trim(ADJUSTL(str1))//strspace
     +                  //trim(ADJUSTL(str2))//strspace
     +                  //trim(ADJUSTL(str3))
               write(*,*) comstr
               call lammps_command (lmp,comstr)
               natomlmp = natomlmp + 1
             enddo
           else
             ndel = natomlmp-natom
             do i = 1,ndel
               str1 = 'group'
               str2 = ''
               write(str2,*) groupkill
               str3 = 'id'
               str4 = ''
               write(str4,*) natomlmp-(i-1)
               comstr =   trim(ADJUSTL(str1))//strspace
     +                  //trim(ADJUSTL(str2))//strspace
     +                  //trim(ADJUSTL(str3))//strspace
     +                  //trim(ADJUSTL(str4))
               write(*,*) comstr
               call lammps_command (lmp,comstr)
               str1 = 'delete_atoms group'
               str2 = ''
               write(str2,*) groupkill
               groupkill = groupkill + 1
               comstr =   trim(ADJUSTL(str1))//strspace
     +                  //trim(ADJUSTL(str2))
               write(*,*) comstr
               call lammps_command (lmp,comstr)
               natomlmp = natomlmp - 1
             enddo
           endif
           !if(allocated(xyzlmp)) deallocate(xyzlmp)
           !if(allocated(F_lmp)) deallocate(F_lmp)
           !if(allocated(itypelmp)) deallocate(itypelmp)
           call lammps_gather_atoms (lmp, 'type', 1, itypelmp)
           call lammps_gather_atoms (lmp, 'x', 3, xyzlmp)
         endif

         do i = 1, nmove*3
           xyzlmp(i) = realbuffer(3+i)
           if(lgettypes.and.(i.le.nmove))then
             itypelmp(i) = realbuffer(3+nmove*3+i)
           endif
         enddo
         if(lgettypes) call lammps_scatter_atoms (lmp, 'type', itypelmp)
         call lammps_scatter_atoms (lmp, 'x', xyzlmp)

         call lammps_command (lmp, 'run 0')

         call lammps_gather_atoms (lmp, 'f', 3, F_lmp)
         call lammps_extract_compute (energyl,lmp,'1',0,0)
         energy = energyl * (1.0/27.21138506)

         imess_size = 1+nmove*3
         !deallocate(realbuffer)
         !allocate(realbuffer(imess_size))
         ! Construct message for force-call
         realbuffer(1) = energy
         do i = 1, nmove*3
           realbuffer(1+i) = -F_lmp(i) * (1.0/27.21138506)
         enddo
         if(rank_lmp.eq.0) call MPI_Send(realbuffer,imess_size,
     +       MPI_REAL8,0,lmp_tag6,parent,ier)
         call MPI_BARRIER(MPI_COMM_WORLD,ier)
         !deallocate(realbuffer)

         !if(allocated(xyzlmp)) deallocate(xyzlmp)
         !if(allocated(F_lmp)) deallocate(F_lmp)
         !if(allocated(itypelmp)) deallocate(itypelmp)

      enddo

      ! Finalize MPI
      call MPI_Finalize (ier)

      stop
      end
