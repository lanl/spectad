       !---------------------------------------------------------------------
       !
       !  Module containing Lammps-Related Info/Subroutines
       !        
       !---------------------------------------------------------------------

       MODULE mod_lammps

          use LAMMPS
          use mod_mpi
          use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
          !include 'mpif.h'
          type (C_ptr) :: lmp
          real*8, dimension(:), allocatable :: taxeslmp
          integer :: groupkill = 100
          logical :: lqeqfix = .true.
          integer natomlmp
          integer natomlmp_i
          integer :: lmpunits = 1
          logical :: doblocks = .false. ! Turn ON to do MD stepping in LAMMPS
          logical :: firstmdblock = .true.
          
         !real (C_double), pointer :: energyl => NULL()
         !integer, dimension(:), allocatable :: itypelmp
         !double precision, dimension(:), allocatable :: qlmp
         double precision, dimension(:), allocatable :: xyzlmp
         double precision, dimension(:), allocatable :: vlmp
         !double precision, dimension(:), allocatable :: F_lmp

       CONTAINS
       
!       subroutine lmp_setup(natom,nmove,ntype,maxtyp,potnam,ietype
!     +     ,xyz,taxes,itype,rcut,amu,azero,iseed,templow,temphigh
!     +     ,dt,lmp_end)
       subroutine lmp_setup(natom,nmove,ntype,maxtyp,potnam,ietype
     +     ,xyz,taxes,itype,rcut,amu,templow,temphigh
     +     ,lmp_end,dt)
       
         !use LAMMPS
         !use mod_lammps
         !use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
         implicit real*8(a-h,o-z)
         include 'parameters.h'
         !include 'mpif.h'
         integer, dimension(:), allocatable :: imessage
         dimension taxes(3)
         dimension itype(mxatom)
         dimension xyz(3*mxatom)
         dimension pxyz(3*mxatom)
         dimension rcut(maxtyp,maxtyp),amu(maxtyp),azero(maxtyp)
         character*80 str1,str2,str3,str4,hostname
         character*80 startnam
         character*80 inputnam
         character*80 filnam
         character*80 potnam(maxtyp)
         character*80 potfilnam,potfilnam1,potfilnam2
         character*80 potdirnam,potcomb
         character*80 filn,filnb
         character*160 lmpline
         real*8 mass_lmp(10)
         real*8 charge_lmp(10)
         !double precision, dimension(:), allocatable :: xyzlmp
         double precision, dimension(:), allocatable :: qlmp
         real (C_double), dimension(:), pointer :: masslmp => NULL()
         logical lmp_end, l_minimize_cell
       
        ! EAM_SETUP() Usually called here when NOT Lammps...
        !   call eam_setup(ntype,maxtyp,potnam,iwrt,ietype,rcut,amu,azero)
        ! Usually passed-in:
        !   ntype   = number of types expected  (passed in)
        !   maxtyp  = max. no. of types - (important: sets 2-D dimension rcut array)
        ! Usually returned:
        !   ietype  = energy type index.  For EAM, ietype=0
        !   rcut    = 2-D array of cutoffs
        !   amu     = 1-D array of masses
        !   azero   = 1-D array of equilibrium T=0 lattice constants

          if(lmp_end)then
              call lammps_command (lmp, 'clear')
              !call lammps_close(lmp)
          endif
          call MPI_BARRIER(force_comm,ier)

          if(allocated(xyzlmp)) deallocate(xyzlmp)
          if(allocated(vlmp)) deallocate(vlmp)

c ========= get dirnam from an environment variable
c          write(6,*) 'getting potsdir name from POTSDIR environment var'
          potdirnam='  '
          call getenv('POTSDIR',potdirnam)
          !potdirnam='../src/POTENTIALS/'
          if(potdirnam.eq.'  ') then
            if((ietype.ne.198).and.(ietype.ne.199))then
              !write(6,*) 'ERROR: must set environment variable POTSDIR'
              !write(6,*) '       to the location of the EAM pot files'
              !stop 'POTSDIR environment variable not set'
              write(6,*)'Assuming potential is in WORKING directory.'
              potdirnam='./'
            endif
          end if
          if(rank_w.eq.0) write(6,*) 'potdirnam =',potdirnam

          ! In LAMMPS-just declare the entire potential name for itype(1)
          if((ntype.eq.2).and.(ietype.eq.101))then
            ! Note.. This means that eam/alloy potential name MUST be
            ! concatinated element names.. ex. FeCr.eam.alloy for Fe-Cr system
            ! ...  Also, only 2 element types allowed
            ! ...  This should be modified in the future
            potcomb = trim(potnam(1))//trim(potnam(2))
          elseif((ntype.eq.2).and.(ietype.eq.102))then
            potcomb = trim(potnam(1))//trim(potnam(2))
          elseif((ntype.eq.2).and.(ietype.eq.103))then
            potcomb = trim(potnam(1))//trim(potnam(2))
          elseif((ntype.eq.2).and.(ietype.eq.104))then
            potcomb = trim(potnam(1))//trim(potnam(2))

          ! Special Case: FeC
          elseif ((ntype.eq.1).and.(trim(potnam(1)).eq.'Fe')) then
            potcomb = trim(potnam(1))//trim(potnam(2))

          else
            ! Note.. This means that eam potential name MUST be 
            ! element name.. ex. Ag.eam for Ag system
            potcomb = potnam(1)
          endif


          
          if(rank_w.eq.0) write(*,*) 'ntype,potcomb = ',ntype,potcomb
          potfilnam = trim(potdirnam)//trim(potcomb)
          
          ! ietype=100 -> EAM
          ! ietype=101 -> EAM/ALLOY
          if(ietype.eq.100)then
            potfilnam = trim(potfilnam)//'.eam'
            if(rank_w.eq.0) write(*,*) 'potfilnam = ',potfilnam
          else if(ietype.eq.101)then
            potfilnam = trim(potfilnam)//'.eam.alloy'
            if(rank_w.eq.0) write(*,*) 'potfilnam = ',potfilnam
          else if(ietype.eq.102)then
            potfilnam = trim(potfilnam)//'.adp'
            if(rank_w.eq.0) write(*,*) 'potfilnam = ',potfilnam
          else if(ietype.eq.103)then
            potfilnam2 = trim(potfilnam)//'.eam.alloy'
            if(rank_w.eq.0) write(*,*) 'potfilnam2 = ',potfilnam2
            potfilnam1 = trim(potfilnam)//'.streitz'
            if(rank_w.eq.0) write(*,*) 'potfilnam1 = ',potfilnam1
          else if(ietype.eq.104)then
            potfilnam = trim(potfilnam)//'.cdeam'
            if(rank_w.eq.0) write(*,*) 'potfilnam = ',potfilnam

          ! Blas' UO2
          else if(ietype.eq.105)then
            potfilnam = trim(potdirnam)//'CeThUNpPuAmCmO.eam.alloy'
            if(rank_w.eq.0) write(*,*) 'potfilnam = ',potfilnam

          else if(ietype.eq.167)then
            if(ntype.gt.1)then
              potfilnam = trim(potfilnam)//'.airebo'
            else
              potfilnam = trim(potfilnam)//'H.airebo'
            endif
            if(rank_w.eq.0) write(*,*) 'potfilnam = ',potfilnam
          else if((ietype.eq.172).or.(ietype.eq.198)
     x                           .or.(ietype.eq.199))then
            potfilnam = 'none'
            if(rank_w.eq.0) write(*,*) 'potfilnam = ',potfilnam
          else
            write(*,*) ' ERROR.. 1xx ietype wrong!! '
          endif
          
c ========= write LAMMPS geometry and input files (rz) ========
          str2 = ''
          write(str2,*) groupID
          str1='lammps.start.'
          startnam = trim(str1)//trim(ADJUSTL(str2))  
          if(rank_l.eq.0)then !if(rank_f.eq.0)then

          ! open *.buck.lmp file for reading if needed:
              ! cutoff
              ! ntype
              ! 1,charge-1,mass-1
              !  ... repeat for other types
              ! npairs
              ! 1,2,1279.69,0.29969,0.0
              ! 2,2,9547.96,0.21916,32.0
              !  ... repeat for all other pairs
          mass_lmp(:) = 0.d0
          charge_lmp(:) = 0.d0
          if((ietype.eq.172))then
            OPEN(UNIT=614,FILE='buck.lmp',STATUS='UNKNOWN')
            READ(614,*) cutoff_lmp
            READ(614,*) ntype_lmp
            do ix = 1, ntype_lmp
              READ(614,*) idummy,charge_lmp(ix),mass_lmp(ix)
            enddo
          elseif((ietype.eq.105))then
            ! NOTE: HARD-CODED for now, but could be read
            !       from file (like ietype=172)
            cutoff_lmp    = 11.0
            ntype_lmp     = 2
            charge_lmp(2) = -1.1104 ! O
            charge_lmp(1) =  2.2208 ! U
            mass_lmp(2)   =  15.99940
            mass_lmp(1)   = 238.02891
          endif

          open(unit=215,file=startnam,status='replace')
          write(215,*) 'AtomisticGeometry'
          write(215,*) ''
          if((ietype.eq.172).or.(ietype.eq.105))then
            write(215,*) natom,' atoms'
            write(215,*) '0 bonds'
            write(215,*) '0 angles'
            write(215,*) '0 dihedrals'
            write(215,*) '0 impropers'
            write(215,*) ntype,' atom types'
            write(215,*) '0 bond types'
            write(215,*) '0 angle types'
            write(215,*) '0 dihedral types'
            write(215,*) '0 improper types'
          else
            write(215,*) natom,' atoms'
            write(215,*) ntype,' atom types'
          endif
          write(215,*) ''
          write(215,*) 0.0,taxes(1),' xlo xhi'
          write(215,*) 0.0,taxes(2),' ylo yhi'
          write(215,*) 0.0,taxes(3),' zlo zhi'
          if((ietype.eq.172).or.(ietype.eq.105))then
            write(215,*) ''
            write(215,*) 'Masses'
            write(215,*) ''
            do ix = 1, ntype_lmp
              write(215,*) ix,mass_lmp(ix)
            enddo
          elseif(ietype.eq.167)then
            write(215,*) ''
            write(215,*) 'Masses'
            write(215,*) ''
            write(215,*) '1 12.01'
            if(ntype.gt.1) write(215,*) '2 1.00'
          endif
          write(215,*) ''
          write(215,*) 'Atoms'
          write(215,*) ''
          do ix = 1, natom
            if((ietype.eq.103).or.(ietype.eq.172)
     x     .or.(ietype.eq.105).or.(ietype.eq.198))then
              q = 0.0
              if((ietype.eq.172).or.(ietype.eq.105))then
                q = charge_lmp(itype(ix))
              endif
              write(215,'(I9,1X,I9,1X,F9.5,1X,F16.8,1X,F16.8,1X,F16.8)')
     x          ix,itype(ix),q,xyz((ix-1)*3+1)
     x          ,xyz((ix-1)*3+2)
     x          ,xyz((ix-1)*3+3)
            else
              write(215,'(I9,1X,I9,1X,F16.8,1X,F16.8,1X,F16.8)')
     x          ix,itype(ix),xyz((ix-1)*3+1)
     x          ,xyz((ix-1)*3+2)
     x          ,xyz((ix-1)*3+3)
            endif
          enddo
          close(215)
          endif
          call MPI_BARRIER(force_comm,ier)
      
          if(.not.allocated(taxeslmp)) allocate(taxeslmp(3))
          taxeslmp(1) = taxes(1)
          taxeslmp(2) = taxes(2)
          taxeslmp(3) = taxes(3)

          str2 = ''
          write(str2,*) groupID
          str1='lammps.input.'
          inputnam = trim(str1)//trim(ADJUSTL(str2))
          if(rank_l.eq.0) then !if(rank_f.eq.0) then

          open(unit=215,file=inputnam,status='replace')

          write(215,*) '#AtomisticSettings'
          write(215,*) '#.....'
          write(215,*) ''
          write(215,*) 'units metal'
          lmpunits = 1
          write(215,*) 'atom_modify map array sort 0 0.0'
          write(215,*) 'boundary p p p'
          if((ietype.eq.103).or.(ietype.eq.172)
     x       .or.(ietype.eq.105).or.(ietype.eq.198))then
            write(215,*) 'atom_style charge'
          else
            write(215,*) 'atom_style atomic'
          endif

          ! This code can be changed for each system...
          nzprocs = 1 
          if((nprocs_f.gt.1).and.(ietype.eq.172))then
            nzprocs = 1
            write(215,*) 'processors * * ',nzprocs
          endif
          write(215,'(A,A)',advance='no') 'read_data ',startnam
          write(215,*) ''
          if((nprocs_f.gt.1).and.(nzprocs.gt.1))then
            write(215,*) 'balance 0.9 shift z 20 1.0' !! Try to get perfect balance..
          endif

          if((ietype.eq.198).or.(ietype.eq.199)) then

            open(100, file='script.lmp')
            ilnmax=250
            iln=0
            do
              iln=iln+1
              lmpline =''
              read(100,'(A)',IOSTAT=io) lmpline
              write(215,'(A)') lmpline
              if (io<0) exit
              if (iln.ge.ilnmax) then
                write(*,*) 'WARNING check script.lmp !!!'
                exit
              endif
            end do
            close(100)

          else ! Closing 'endif' is far down
               ! look for 'if(ietype.eq.199)'

          if(ietype.eq.100)then
            write(215,*) 'pair_style eam'
            write(215,*) 'pair_coeff * * ',trim(potfilnam)
          else if(ietype.eq.101)then
            write(215,*) 'pair_style eam/alloy'
            if(ntype == 2) then
              write(215,*) 'pair_coeff * * ',trim(potfilnam),' ',
     x          trim(potnam(1)),' ',trim(potnam(2))
            else
              write(215,*) 'pair_coeff * * ',trim(potfilnam),' ',
     x          trim(potnam(1))
            endif
          else if(ietype.eq.102)then
            write(215,*) 'pair_style adp'
            if(ntype == 2) then
              write(215,*) 'pair_coeff * * ',trim(potfilnam),' ',
     x          trim(potnam(1)),' ',trim(potnam(2))
            else
              write(215,*) 'pair_coeff * * ',trim(potfilnam),' ',
     x          trim(potnam(1))
            endif
          else if(ietype.eq.104)then
            write(215,*) 'pair_style eam/cd'
            if(ntype == 2) then
              write(215,*) 'pair_coeff * * ',trim(potfilnam),' ',
     x          trim(potnam(1)),' ',trim(potnam(2))
            else
              write(215,*) 'pair_coeff * * ',trim(potfilnam),' ',
     x          trim(potnam(1))
            endif
          else if(ietype.eq.167)then
            write(215,*) 'pair_style airebo 3.0'
            if(ntype == 2) then
              write(215,*) 'pair_coeff * * ',trim(potfilnam),' ',
     x          trim(potnam(1)),' ',trim(potnam(2))
            else
              write(215,*) 'pair_coeff * * ',trim(potfilnam),' ',
     x          trim(potnam(1))
            endif
          elseif(ietype.eq.172)then
            write(215,*)''
            write(215,'(A,F16.8)',advance='no')
     x            'pair_style buck/coul/long ',cutoff_lmp
            write(215,*)''
            write(215,*) 'kspace_style ewald 1.0e-5'
            write(215,*) 'pair_modify table 0'
            READ(614,*) npairs_lmp
            do ix = 1, npairs_lmp
              READ(614,*) idum,jdum,par1,par2,par3
              write(215,'(A,I9,1X,I9,1X,F16.8,1X,F16.8,1X,F16.8)')
     x            'pair_coeff ',idum,jdum,par1,par2,par3
            enddo
            CLOSE(614)
          elseif(ietype.eq.105)then

            write(215,*) 'kspace_style pppm 1.0e-5'
            write(215,*)''
            write(215,'(A,F16.8,A)',advance='no')
     x            'pair_style hybrid/overlay coul/long '
     x            ,cutoff_lmp,' eam/alloy'
            write(215,*)''
            write(215,*) 'pair_coeff * * coul/long'
            write(215,'(A,A,A)',advance='no')
     x            'pair_coeff * * eam/alloy ',trim(potfilnam)
     x            ,' U O'

          elseif(ietype.eq.103)then
            if(ntype.ne.2) then
              STOP 'NOT set up for one atom type here!'
            endif
            write(215,*) 'group    	type1 type 1'
            write(215,*) 'compute   	charge1 type1 property/atom q'
            write(215,*) 'compute   	q1 type1 reduce ave c_charge1'
            write(215,*) 'group    	type2 type 2'
            write(215,*) 'compute   	charge2 type2 property/atom q'
            write(215,*) 'compute   	q2 type2 reduce ave c_charge2'
            write(215,*) ''
            write(215,*) 'variable   	qcat equal 2.8'
            write(215,*) 'variable  ',
     x          'qani equal -${qcat}*count(type1)/count(type2)'
            write(215,*) 'set   		group type1 charge ${qcat}'
            write(215,*) 'set   		group type2 charge ${qani}'
            write(215,*) 'variable  ',
     x          'qsum equal count(type1)*c_q1+count(type2)*c_q2'
            write(215,*) ''
            write(215,*) 'pair_style  	',
     x          'hybrid/overlay coul/streitz 12.0 ewald eam/alloy'
            write(215,*) 'kspace_style  	ewald 1e-7'
            write(215,*) 'pair_modify table 0'
            write(215,*) ''
            write(215,*) 'pair_coeff * * coul/streitz '
     x                                ,trim(potfilnam1),' Al O'
            write(215,*) 'pair_coeff * * eam/alloy '
     x                                ,trim(potfilnam2),' Al O'
          else
            write(*,*) ' ERROR.. 1xx ietype wrong!! '
          endif


          endif ! if(ietype.eq.199)


          write(215,*) ''
          if((ietype.eq.172))then
            write(215,*) 'neighbor 8.0 bin'
            write(215,*) 'neigh_modify every 1 delay 0 check no'
          else
            write(215,*) 'neighbor 2.0 bin'
            write(215,*) 'neigh_modify every 1 delay 0 check yes'
          endif
          write(215,*) ''

          !write(215,*) 'group frozen id > ',nmove
          !write(215,*) 'group free id < ',nmove+1
          !write(215,*) ''

          if(ietype.eq.103)then
!            write(215,*) 'thermo_style custom step temp etotal',
!     x          ' pe evdwl ecoul elong & '
!            write(215,*) ' c_q1 c_q2 v_qsum press spcpu '
            write(215,*)
     x        'thermo_style custom step pe etotal v_qsum c_q1 c_q2 '
            write(215,*) 'thermo_modify  	norm yes'
            write(215,*) 'thermo   	10'
            write(215,*) 'fix ',
     x       '1 all qeq/slater 1 12.0 1.0e-7 10000 coul/streitz'
          endif
          !write(215,*) 'fix 2 frozen setforce 0.0 0.0 0.0'
          !if(natom.eq.nmove)then
          !  write(215,*) 'fix 3 all recenter INIT INIT INIT'
          !endif

          !! LANGEVIN MD SETTINGS:
          write(215,*) ''
!          write(215,*) 'fix 5 all nve'
!          write(215,*) 'fix 6 all langevin '
!     x        ,tmephigh,temphigh,(1.d0/thermrate),' 48279'
          !!write(215,*) 'velocity all create 500.0 4928459 rot no dist gaussian'
          !!write(215,*) 'thermo 100'
          if(lmpunits.eq.1) then
              write(215,*) 'timestep ',(dt/(1.d-12)) ! ONLY GOOD FOR UNITS=METAL
          else
              ! ERROR - these units are not allowed yet
          endif
          !write(215,*) 'thermo 10'
          write(215,*) 'group frozen id > ',nmove

          write(215,*) ''
          write(215,*) 'compute 1 all pe'
          !if((.not.lmp_end).and.(rank_log.eq.1))then
            write(215,*) ''
            write(215,*) 'min_style fire #sd'
            write(215,*) '#min_modify dmax 0.05'
          !endif
          write(215,*) ''
          !write(215,*) 'thermo 100000'
!          write(215,*) 'timestep 0.002'
          write(215,*) '' 
          !str1 = ''
          !write(str1,*) groupID
          !str2 = 'lmpatoms.'//trim(ADJUSTL(str1))
          !write(215,*) 'dump 1 all custom 1 ',str2,' id x y z'
          !write(215,*) ''

          write(215,*) 'run 0'
          close(215)

          endif
          call MPI_BARRIER(force_comm,ier)
          !call sleep(10)

          ! Initialize natomlmp:
          natomlmp = natom
          natomlmp_i = natom

c ========= Start-up Lammps (lmp) Instance and read input files...

          if(.not.spawn_lammps)then !! Do things in mod_lammps?

           if(.not.lmp_end)then
            if(rank_f.eq.0) write(*,*) 'STARTING LAMMPS -- lmp'
            if(rank_log.eq.0)then
              call MPI_BARRIER(force_comm,ier)
              if(.true.)then
                call lammps_open ('lmp -log log.lammps -screen none',
     x              force_comm, lmp)
              else
                call lammps_open ('lmp -log none -screen none',
     x                                    force_comm, lmp)
              endif
            else
              call MPI_BARRIER(force_comm,ier)
              if(.false.)then
                str1 = ''
                write(str1,*) rank_log
                str2 = 'lmp -log log.lammps.'//trim(ADJUSTL(str1))//
     x             ' -screen none'
                call lammps_open (str2, force_comm, lmp)
              else
                call lammps_open ('lmp -log none -screen none',
     x                                 force_comm, lmp)
              endif
            endif
           endif

           call lammps_file (lmp,inputnam)
           write(*,*) 'rank_log ',rank_log,' good.. '
          
           if(rank_f.eq.0) write(*,*) 'LAMMPS STARTED -- lmp'


          ! HACK to get T-dependent lattice constant:
          ! (Set up for MgAl2O4 at the moment)
          if(.FALSE.)then

            ! 500K
            call lammps_command (lmp,
     +  'fix 5 all npt temp 500.0 500.0 0.1 iso 0.0 0.0 1.0')
            call lammps_command (lmp,
     +  'velocity all create 500.0 4928459 rot no dist gaussian')
            call lammps_command (lmp, 'thermo 100')
            call lammps_command (lmp, 'timestep 0.001')
            call lammps_command (lmp, 'run 10000')

            ! 1000K
            call lammps_command (lmp,'unfix 5')
            call lammps_command (lmp,
     +  'fix 5 all npt temp 1000.0 1000.0 0.1 iso 0.0 0.0 1.0')
            call lammps_command (lmp,
     +  'velocity all create 1000.0 4928460 rot no dist gaussian')
            call lammps_command (lmp, 'run 10000')

            ! 1500K
            call lammps_command (lmp,'unfix 5')
            call lammps_command (lmp,
     +  'fix 5 all npt temp 1500.0 1500.0 0.1 iso 0.0 0.0 1.0')
            call lammps_command (lmp,
     +  'velocity all create 1500.0 4928461 rot no dist gaussian')
            call lammps_command (lmp, 'run 10000')

            ! 2000K
            call lammps_command (lmp,'unfix 5')
            call lammps_command (lmp,
     +  'fix 5 all npt temp 2000.0 2000.0 0.1 iso 0.0 0.0 1.0')
            call lammps_command (lmp,
     +  'velocity all create 2000.0 4928462 rot no dist gaussian')
            call lammps_command (lmp, 'run 10000')
            
            stop 'THIS WAS AN MD HACK TO GET ao(T) !!' 

          
          elseif(.FALSE.)then

            ! 2000K
            call lammps_command (lmp,'fix 5 all nve')
            call lammps_command (lmp,
     +                    'fix 6 all langevin 2000.0 2000.0 0.1 48279')
            call lammps_command (lmp,
     +        'velocity all create 500.0 4928459 rot no dist gaussian')
            call lammps_command (lmp, 'thermo 100')
            call lammps_command (lmp, 'timestep 0.002')

            ii = 0
            pxyz(:) =0.d0
            do ii = 1, 1000
              call lammps_command (lmp, 'run 1000')
              call lammps_gather_atoms (lmp, 'x', 3, xyzlmp)
                call flush(6)
                filn=''
                write(filn,*) spawnID
                filnb=''
                write(filnb,*) ii
                filn='mdtest.'//trim(adjustl(filn))//
     +               '.'//trim(adjustl(filnb))//'.dat'
                call storefile(natom,xyzlmp,pxyz,itype,taxes,filn)
                call flush(6)
                call MPI_BARRIER(force_comm,ier)
            enddo
            
            stop 'THIS WAS A HACK TO RUN MD IN LAMMPS !!' 

          endif


c ========= Get Masses (and ao?) From Lammps...
           if(.not.lmp_end)then
            call lammps_extract_atom (masslmp, lmp, 'mass')
            if(rank_f.eq.0) write(*,*) 'mass masslmp -- ',masslmp
            do ix = 1, ntype
              amu(ix) = masslmp(ix+1)*1.0 !! Assumes point bound remapping is NOT used in LAMMPS.F90
              if(rank_f.eq.0) write(*,*) 'mass -- ',amu(ix)
              ! crz
              ! Lattice constants are currently hard-coded here..
              ! Not sure this matters if these values are never used.. 
              ! (need to look into this)
              if(potnam(ix) == 'Ag') then
                azero(ix) = 4.09
              else
                azero(ix) = 4.09 !crz <- Temporary fix
              endif
            enddo
            masslmp => NULL()
           endif

c ========= Minimize Geometry in Lammps (If this is SpawnID = 1)...
            if((.not.lmp_end).and.(rank_log.eq.1))then ! can only minimize with proper fix
              call
     +           lammps_command(lmp,'minimize 0.0 0.001 100000 1000000')

              l_minimize_cell = .false.
              if(l_minimize_cell)then
                call lammps_command(lmp,
     +            'fix 7 all box/relax iso 0.0 vmax 0.001')
                call lammps_command(lmp,
     +            'dump 1 all cfg 10 dump.*.out mass type xs ys zs')
                call lammps_command(lmp,'min_style cg')
                call
     +           lammps_command(lmp,'minimize 0.0 0.001 100000 1000000')
              endif

              call lammps_gather_atoms (lmp, 'x', 3, xyzlmp)
              call lammps_gather_atoms (lmp, 'v', 3, vlmp)
              do i = 1, natom
                do j = 1, 3
                  xyz((i-1)*3+j) = xyzlmp((i-1)*3+j)
                enddo
              enddo
            else
              call lammps_gather_atoms (lmp, 'x', 3, xyzlmp)
              call lammps_gather_atoms (lmp, 'v', 3, vlmp)
              do i = 1, natom
                do j = 1, 3
                  xyzlmp((i-1)*3+j) = xyz((i-1)*3+j)
                  vlmp((i-1)*3+j) = 0.0
                enddo
              enddo
            endif
            !call lammps_command (lmp, 'min_modify dmax 0.0')
            !if(allocated(xyzlmp)) deallocate(xyzlmp)

            ! Lets not update charges unless we are starting in a new state...
            ! (May cause issues for state recognition??)
            if((ietype.eq.103).or.(ietype.eq.198))then
               call lammps_gather_atoms (lmp, 'q', 1, qlmp)
               do i = 1, natom
                 tadcharge(i) = qlmp(i)
               enddo
               deallocate(qlmp)
               call lammps_command (lmp, 'unfix 1')
            endif

          else ! TRY USING MPI_SPAWN INSTEAD...

             call HOSTNM(hostname)
             write(*,*) 'hostname = ',hostname
             call MPI_Info_Create(host_info, ier)
             call MPI_Info_Set(host_info, 'host', hostname, ier)
             str1 = 'lammps_comm.exe'
             write(*,*) 'spawning ',str1
             call MPI_COMM_SPAWN(str1, MPI_ARGV_NULL, nforcecoresSpawn,
     +           host_info, 0, MPI_COMM_SELF, comm_everyone,
     +           MPI_ERRCODES_IGNORE,ier)
             call MPI_Comm_size (comm_everyone, nprocs_ev, ier)
             write(*,*) 'nprocs_ev ',nprocs_ev

             if(rank_log.eq.0)then

              if(.false.)then
                str1 = 'lmp -log log.lammps -screen none'
              else
                str1 = 'lmp -log none -screen none'
              endif
              imess_size = LEN(str1)
              lmp_tag1 = 1
              call MPI_Send(str1,imess_size,MPI_CHARACTER
     +               ,0,lmp_tag1,comm_everyone,ier)
              allocate(imessage(4))
              imessage(1) = ntype
              imessage(2) = natom
              imessage(3) = itypestart(1)
              imessage(4) = 0 ! Dont bother minimizing
              lmp_tag2 = 2
              call MPI_Send(imessage,4,MPI_INTEGER
     +               ,0,lmp_tag2,comm_everyone,ier)
              deallocate(imessage)
              lmp_tag3 = 3
              call MPI_Recv(amu,ntype,MPI_REAL8
     +               ,0,lmp_tag3,comm_everyone,msgstatus,ier)
              do ix = 1, ntype
                write(*,*) 'mass -- ',amu(ix)
                ! crz
                ! Lattice constants are currently hard-coded here..
                ! Not sure this matters if these values are never used..
                ! (need to look into this)
                if(potnam(ix) == 'Ag') then
                  azero(ix) = 4.09
                else
                  azero(ix) = 4.09 !crz <- Temporary fix
                endif
              enddo

             else

              if(.true.)then
                str1 = ''
                write(str1,*) rank_log
                str2 = 'lmp -log log.lammps.'//trim(ADJUSTL(str1))//
     x             ' -screen none'
                str1 = str2
              else
                str1 = 'lmp -log none -screen none'
              endif
              imess_size = LEN(str1)
              lmp_tag1 = 1
              call MPI_Send(str1,imess_size,MPI_CHARACTER
     +               ,0,lmp_tag1,comm_everyone,ier)
              allocate(imessage(4))
              imessage(1) = ntype
              imessage(2) = natom
              imessage(3) = itypestart(1)
              imessage(4) = 0 ! Dont bother minimizing
              if(rank_log.eq.1) imessage(4) = 1
              lmp_tag2 = 2
              call MPI_Send(imessage,4,MPI_INTEGER
     +               ,0,lmp_tag2,comm_everyone,ier)
              deallocate(imessage)
              lmp_tag3 = 3
              call MPI_Recv(amu,ntype,MPI_REAL8
     +               ,0,lmp_tag3,comm_everyone,msgstatus,ier)
              do ix = 1, ntype
                write(*,*) 'mass -- ',amu(ix)
                ! crz
                ! Lattice constants are currently hard-coded here..
                ! Not sure this matters if these values are never used..
                ! (need to look into this)
                if(potnam(ix) == 'Ag') then
                  azero(ix) = 4.09
                else
                  azero(ix) = 4.09 !crz <- Temporary fix
                endif
              enddo
              if (rank_log.eq.1) then
                lmp_tag4 = 4
                call MPI_Recv(xyz,3*natom,MPI_REAL8
     +               ,0,lmp_tag4,comm_everyone,msgstatus,ier)
              endif

             endif
          endif

        return
       end subroutine lmp_setup
       
       
      subroutine g100(natom,nmove,x,y,z,itype,ietype,maxtyp,energy,grad
     x               ,taxes)

         implicit real*8(a-h,o-z)
         save
         include 'parameters.h'
         real (C_double), pointer :: energyl => NULL()
         integer, dimension(:), allocatable :: itypelmp
         double precision, dimension(:), allocatable :: qlmp
         !double precision, dimension(:), allocatable :: xyzlmp
         double precision, dimension(:), allocatable :: F_lmp
         integer i,j
         dimension x(mxatom),y(mxatom),z(mxatom),itype(natom)
         dimension taxes(3)
         dimension grad(3,natom)
         character(LEN=80) str1,str2,str3,str4,comstr
         character(LEN=1) strspace
         logical lgettypes

         common/callcom/ngradcall

         call MPI_BARRIER(force_comm,ier)

         ! To be safe... deallocate everything...
         !if(allocated(xyzlmp)) deallocate(xyzlmp)
         if(allocated(F_lmp)) deallocate(F_lmp)
         if(allocated(itypelmp)) deallocate(itypelmp)
         if(allocated(qlmp)) deallocate(qlmp)

         ! Allow lammps to allocate arrays before modification
         if(natom.eq.natomlmp)then
           !call lammps_gather_atoms (lmp, 'type', 1, itypelmp)
           !call lammps_gather_atoms (lmp, 'x', 3, xyzlmp)
           if((ietype.eq.103).or.(ietype.eq.198))then
             call lammps_gather_atoms (lmp, 'q', 1, qlmp)
           endif
         endif

         ! Check that atom count in lammps is correct
         !natomlmp = SIZE(itypelmp)
         lgettypes = .false.
         if(natom.ne.natomlmp)then
           strspace = ' '
           lgettypes = .true.
           if(natom.gt.natomlmp)then
             nadd = natom-natomlmp
             do i = 1,nadd
               iadd = (natomlmp+i)-natomlmp_i
               str1 = 'create_atoms'
               str2 = ''
               write(str2,*) 1 !itypestart(1)
               !write(str2,*) itypedeposit(iadd)
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
           call lammps_gather_atoms (lmp, 'type', 1, itypelmp)
           if(allocated(xyzlmp)) deallocate(xyzlmp)
           call lammps_gather_atoms (lmp, 'x', 3, xyzlmp)
           if((ietype.eq.103).or.(ietype.eq.198))then
             call lammps_gather_atoms (lmp, 'q', 1, qlmp)
           endif
         endif

         do i = 1, natom
           xyzlmp((i-1)*3+1) = x(i)
           xyzlmp((i-1)*3+2) = y(i)
           xyzlmp((i-1)*3+3) = z(i)
           if(lgettypes) itypelmp(i) = itype(i)
           if((ietype.eq.103).or.(ietype.eq.198))then
              qlmp(i) = tadcharge(i)
           endif
         enddo
         
         if(lgettypes)call lammps_scatter_atoms (lmp, 'type', itypelmp)
         call lammps_scatter_atoms (lmp, 'x', xyzlmp)
         if((ietype.eq.103).or.(ietype.eq.198))then
           call lammps_scatter_atoms (lmp, 'q', qlmp)
         endif

         call lammps_command (lmp, 'run 0')

         call lammps_gather_atoms (lmp, 'f', 3, F_lmp)
         do i = 1, natom
           do j = 1, 3
             grad(j,i) = -F_lmp((i-1)*3+j) * (1.0/27.21138506)
           enddo
         enddo
      
         call lammps_extract_compute (energyl,lmp,'1',0,0)
         energy = energyl * (1.0/27.21138506)

         ! Get updated charges
         if((ietype.eq.103).or.(ietype.eq.198))then
           call lammps_gather_atoms (lmp, 'q', 1, qlmp)
           do i = 1, natom
             tadcharge(i) = qlmp(i)
           enddo
         endif

         !if(allocated(xyzlmp)) deallocate(xyzlmp)
         if(allocated(F_lmp)) deallocate(F_lmp)
         if(allocated(itypelmp)) deallocate(itypelmp)
         if(allocated(qlmp)) deallocate(qlmp)

        return
       end  subroutine g100
       
       
       
      subroutine g100spawn(natom,nmove,xyz,ietype,itype,energy,grad)

         implicit real*8(a-h,o-z)
         save
         include 'parameters.h'
         real*8 fcallbuffer(3+mxatom*4)
         integer i,j
         dimension xyz(mxatom*3),itype(natom)
         dimension grad(3,natom)

         ! Send Force-Call Info
         if(natomlmp.ne.natom)then
           imess_size = 3+nmove*3+nmove
           fcallbuffer(1) = REAL(natom,8)
           fcallbuffer(2) = REAL(nmove,8)
           fcallbuffer(3) = REAL(ietype,8)
           do i = 1, nmove*3
             fcallbuffer(3+i) = xyz(i)
             if(i.le.natom) fcallbuffer(3+natom*3+i) = REAL(itype(i),8)
           enddo
           natomlmp = natom
         else
           imess_size = 3+nmove*3
           fcallbuffer(1) = REAL(natom,8)
           fcallbuffer(2) = REAL(nmove,8)
           fcallbuffer(3) = REAL(ietype,8)
           do i = 1, nmove*3
             fcallbuffer(3+i) = xyz(i)
           enddo
         endif
         lmp_tag5 = 5
         call MPI_Send(fcallbuffer,imess_size,MPI_REAL8
     +       ,0,lmp_tag5,comm_everyone,ier)

         ! Get Forces
         imess_size = 1+nmove*3
         lmp_tag6 = 6
         call MPI_Recv(fcallbuffer,imess_size,MPI_REAL8
     +       ,0,lmp_tag6,comm_everyone,msgstatus,ier)
         energy = fcallbuffer(1)
         do i = 1, nmove
           do j = 1, 3
             grad(j,i) = fcallbuffer(1+(i-1)*3+j)
           enddo
         enddo

        return
       end  subroutine g100spawn


      subroutine unfix_qeq()
         implicit real*8(a-h,o-z)
         save
         include 'parameters.h'
         call lammps_command (lmp,'unfix 1')
         call lammps_command (lmp,'run 0')
         lqeqfix = .false.
        return
      end  subroutine unfix_qeq


      subroutine refix_qeq()
         implicit real*8(a-h,o-z)
         save
         include 'parameters.h'
         character(LEN=80) str1
         str1 = ''
         write(str1,*) 'fix ',
     x       '1 all qeq/slater 1 12.0 1.0e-7 10000 coul/streitz'
         write(*,*) str1
         call lammps_command (lmp,str1)
         call lammps_command (lmp,'run 0')
         lqeqfix = .true.
        return
       end  subroutine refix_qeq
       
       
      subroutine lmp_minimize(natom,nmove,xyz,itype,ietype
     x       ,energy,taxes)

         implicit real*8(a-h,o-z)
         save
         include 'parameters.h'
         real (C_double), pointer :: energyl => NULL()
         integer, dimension(:), allocatable :: itypelmp
         double precision, dimension(:), allocatable :: qlmp
         !double precision, dimension(:), allocatable :: xyzlmp
         integer i,j
         dimension xyz(natom*3),itype(natom)
         dimension taxes(3)

         common/callcom/ngradcall

         call MPI_BARRIER(force_comm,ier)

         ! To be safe... deallocate everything...
         !if(allocated(xyzlmp)) deallocate(xyzlmp)
         if(allocated(itypelmp)) deallocate(itypelmp)
         if(allocated(qlmp)) deallocate(qlmp)

         ! Allow lammps to allocate arrays before modification
         call lammps_gather_atoms (lmp, 'type', 1, itypelmp)
         call lammps_gather_atoms (lmp, 'x', 3, xyzlmp)
         if((ietype.eq.103).or.(ietype.eq.198))then
           call lammps_gather_atoms (lmp, 'q', 1, qlmp)
         endif

         do i = 1, natom
           xyzlmp((i-1)*3+1) = xyz((i-1)*3+1)
           xyzlmp((i-1)*3+2) = xyz((i-1)*3+2)
           xyzlmp((i-1)*3+3) = xyz((i-1)*3+3)
           itypelmp(i) = itype(i)
           if((ietype.eq.103).or.(ietype.eq.198))then
             qlmp(i) = tadcharge(i)
           endif
         enddo
         
         call lammps_scatter_atoms (lmp, 'type', itypelmp)
         call lammps_scatter_atoms (lmp, 'x', xyzlmp)
         if((ietype.eq.103).or.(ietype.eq.198))then
           call lammps_scatter_atoms (lmp, 'q', qlmp)
         endif

         call lammps_command (lmp,'minimize 0.0 0.001 1000 100000')
      
         call lammps_extract_compute (energyl,lmp,'1',0,0)
         energy = energyl * (1.0/27.21138506)

         ! Get updated charges
         if((ietype.eq.103).or.(ietype.eq.198))then
           call lammps_gather_atoms (lmp, 'q', 1, qlmp)
           do i = 1, natom
             tadcharge(i) = qlmp(i)
           enddo
         endif

         ! Get updated coords
         call lammps_gather_atoms (lmp, 'x', 3, xyzlmp)
         do i = 1, natom
           xyz((i-1)*3+1) = xyzlmp((i-1)*3+1)
           xyz((i-1)*3+2) = xyzlmp((i-1)*3+2)
           xyz((i-1)*3+3) = xyzlmp((i-1)*3+3)
         enddo

         !if(allocated(xyzlmp)) deallocate(xyzlmp)
         if(allocated(itypelmp)) deallocate(itypelmp)
         if(allocated(qlmp)) deallocate(qlmp)

        return
       end  subroutine lmp_minimize

      subroutine vlmptopxyz(natom,nmove,pxyz,itype,amu,ntype,maxtyp)
         implicit real*8(a-h,o-z)
         dimension amu(maxtyp)
         dimension rmass(10)
         dimension pxyz(3*natom)
         dimension itype(natom)
         do it = 1, ntype
             rmass(it) = amu(it) * 1822.83d0
         enddo
         do i = 1, natom
             do j = 1, 3
               if(lmpunits.eq.1)then
                 ! Convert Ang/ps velosity [lammps - metal] to a.u. velocity
                 v_conv = vlmp((i-1)*3+j) / (2.18769126e4)
               else
                 !! ERROR - this is not an option yet
               endif
               pxyz((i-1)*3+j) = rmass(itype(i))*v_conv !/10.0
             enddo
         enddo
        return
       end  subroutine vlmptopxyz

      subroutine pxyztovlmp(natom,nmove,pxyz,itype,amu,ntype,maxtyp)
         implicit real*8(a-h,o-z)
         dimension amu(maxtyp)
         dimension rmass(10)
         dimension pxyz(3*natom)
         dimension itype(natom)
         do it = 1, ntype
             rmass(it) = amu(it) * 1822.83d0
             !write(*,*) ' rmass, it -- ',rmass, it
             !call flush(6)
         enddo
         !if(.false.)then
         do i = 1, natom
             do j = 1, 3
               ! Get a.u. velocity from a.u. momentum
               vxyz = pxyz((i-1)*3+j) / rmass(itype(i))
               if(lmpunits.eq.1)then
                 ! Convert  a.u. velocity to Ang/ps velosity [lammps - metal]
                 vlmp((i-1)*3+j) = vxyz*(2.18769126e4)
               else
                 !! ERROR - this is not an option yet
               endif
             enddo
         enddo
         !endif
        return
       end  subroutine pxyztovlmp
       
      subroutine mdblock100(natom,nmove,xyz,pxyz,itype,maxtyp
     x   ,energy,grad,taxes,ietype,amu,temp,dt
     x   ,thermrate,nmd,iseed,ntype)

         implicit real*8(a-h,o-z)
         save
         include 'parameters.h'
         real (C_double), pointer :: energyl => NULL()
         integer, dimension(:), allocatable :: itypelmp
         double precision, dimension(:), allocatable :: qlmp
         !double precision, dimension(:), allocatable :: xyzlmp
         double precision, dimension(:), allocatable :: F_lmp
         integer i,j
         dimension xyz(3*mxatom),pxyz(3*mxatom),itype(natom)
         dimension taxes(3)
         dimension amu(maxtyp)
         dimension grad(3,natom)
         character(LEN=80) str0,str1,str2,str3,str4
         character(LEN=120) comstr
         character(LEN=1) strspace
         logical lagain

         common/callcom/ngradcall

         !write(*,*) 'rank_f ',rank_f,' entering mdblock100'
         call MPI_BARRIER(force_comm,ier)

         ! To be safe... deallocate everything...
         if(allocated(F_lmp)) deallocate(F_lmp)
         if(allocated(itypelmp)) deallocate(itypelmp)
         if(allocated(qlmp)) deallocate(qlmp)

         do i = 1, natom
           xyzlmp((i-1)*3+1) = xyz((i-1)*3+1)
           xyzlmp((i-1)*3+2) = xyz((i-1)*3+2)
           xyzlmp((i-1)*3+3) = xyz((i-1)*3+3)
           if((ietype.eq.103).or.(ietype.eq.198))then
             qlmp(i) = tadcharge(i)
           endif
         enddo
         call pxyztovlmp(natom,nmove,pxyz,itype,amu,ntype
     x       ,maxtyp)

         call lammps_scatter_atoms (lmp, 'x', xyzlmp)
         call lammps_scatter_atoms (lmp, 'v', vlmp)
         if((ietype.eq.103).or.(ietype.eq.198))then
           call lammps_scatter_atoms (lmp, 'q', qlmp)
         endif

         !! Add MD FIXES:
         call lammps_command (lmp,'fix 5 all nve')
         strspace = ' '
         str0 = ''
         write(str0,*) 'fix 6 all langevin'
         str1 = ''
         write(str1,*) (temp)
         str2 = ''
         write(str2,*) (temp)
         str3 = ''
         if(lmpunits.eq.1)then
           write(str3,*) ((1.d0/thermrate)*(1.d12)) ! for units=metal ONLY
         else
           !! ERROR - this is not an option yet
         endif
         str4 = ''
         if(rank_f.eq.0)
     +       iseedrand = abs(CEILING(prngen(0)*1000000))
         call MPI_BCAST(iseedrand,1,MPI_INTEGER,0,force_comm,ier)
         write(str4,*) iseedrand
         comstr =   trim(ADJUSTL(str0))//strspace
     +            //trim(ADJUSTL(str1))//strspace
     +            //trim(ADJUSTL(str2))//strspace
     +            //trim(ADJUSTL(str3))//strspace
     +            //trim(ADJUSTL(str4))
         !write(*,*) ' RICKTEST 1 -- ',comstr
         call lammps_command (lmp,comstr)
         call lammps_command (lmp,'fix 7 frozen setforce 0.0 0.0 0.0')
!         if (natom.eq.nmove) then
!           call lammps_command
!     +            (lmp,'fix 8 all recenter INIT INIT INIT')
!         endif
!         if (firstmdblock) then
!           if (rank_log.eq.1) then
!             call lammps_command
!     +            (lmp,'dump myDump all atom 10 dump.atom.*')
!           endif
!         endif

         !! Run MD Block:
         str1 = 'run '
         str2 = ''
         write(str2,*) nmd
         comstr =   trim(ADJUSTL(str1))//strspace
     +                  //trim(ADJUSTL(str2))
         !write(*,*) ' RICKTEST 2 -- ',comstr
         call lammps_command (lmp,comstr)

         ! Move back any atoms that may have crossed pbc's:
         call lammps_gather_atoms (lmp, 'x', 3, xyzlmp)
         do i = 1, natom
           lagain = .true.
           do while (lagain)
             dx = - xyz((i-1)*3+1) + xyzlmp((i-1)*3+1)
             dy = - xyz((i-1)*3+2) + xyzlmp((i-1)*3+2)
             dz = - xyz((i-1)*3+3) + xyzlmp((i-1)*3+3)
             lagain = .false.
             if (dx.gt.(taxes(1)*0.5))then
               xyzlmp((i-1)*3+1) = xyzlmp((i-1)*3+1)
     +                               - taxes(1)
               lagain = .true.
             elseif (dx.lt.(-taxes(1)*0.5))then
               xyzlmp((i-1)*3+1) = xyzlmp((i-1)*3+1)
     +                               + taxes(1)
               lagain = .true.
             endif
             if (dy.gt.(taxes(2)*0.5))then
               xyzlmp((i-1)*3+2) = xyzlmp((i-1)*3+2)
     +                               - taxes(2)
               lagain = .true.
             elseif (dy.lt.(-taxes(2)*0.5))then
               xyzlmp((i-1)*3+2) = xyzlmp((i-1)*3+2)
     +                               + taxes(2)
               lagain = .true.
             endif
             if (dz.gt.(taxes(3)*0.5))then
               xyzlmp((i-1)*3+3) = xyzlmp((i-1)*3+3)
     +                               - taxes(3)
               lagain = .true.
             elseif (dz.lt.(-taxes(3)*0.5))then
               xyzlmp((i-1)*3+3) = xyzlmp((i-1)*3+3)
     +                               + taxes(3)
               lagain = .true.
             endif
           enddo
           xyz((i-1)*3+1) = xyzlmp((i-1)*3+1)
           xyz((i-1)*3+2) = xyzlmp((i-1)*3+2)
           xyz((i-1)*3+3) = xyzlmp((i-1)*3+3)
         enddo

         !! Update momentum in tadcall
         call lammps_gather_atoms (lmp, 'v', 3, vlmp)
         call vlmptopxyz(natom,nmove,pxyz,itype,amu,ntype
     x       ,maxtyp)

         !! REMOVE MD FIXES and execute clean force-call:
         call lammps_command (lmp,'unfix 5')
         call lammps_command (lmp,'unfix 6')
         call lammps_command (lmp,'unfix 7')
!         if (natom.eq.nmove) then
!           call lammps_command (lmp,'unfix 8')
!         endif
!         if (firstmdblock) then
!           firstmdblock = .false.
!           if (rank_log.eq.1) then
!             !call lammps_command (lmp,'undump myDump')
!           endif
!         endif

         call lammps_scatter_atoms (lmp, 'x', xyzlmp)
         call lammps_command (lmp,'run 0')

         call lammps_gather_atoms (lmp, 'f', 3, F_lmp)
         do i = 1, natom
           do j = 1, 3
             grad(j,i)=-F_lmp((i-1)*3+j)*(1.0/27.21138506)
           enddo
         enddo
      
         call lammps_extract_compute (energyl,lmp,'1',0,0)
         energy = energyl * (1.0/27.21138506)

         ! Get updated charges
         if((ietype.eq.103).or.(ietype.eq.198))then
           call lammps_gather_atoms (lmp, 'q', 1, qlmp)
           do i = 1, natom
             tadcharge(i) = qlmp(i)
           enddo
         endif

         !if(allocated(xyzlmp)) deallocate(xyzlmp)
         if(allocated(F_lmp)) deallocate(F_lmp)
         if(allocated(itypelmp)) deallocate(itypelmp)
         if(allocated(qlmp)) deallocate(qlmp)

        return
       end  subroutine mdblock100

       END MODULE mod_lammps
       
       
