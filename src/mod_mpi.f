       !---------------------------------------------------------------------
       !
       !  Module containing mpi/spectad-related info
       !
       !---------------------------------------------------------------------

       MODULE mod_mpi

          use omp_lib
          include 'mpif.h'

          integer ier ! mpi error code
          
          character*80 filstrt,filstrt2,filstrtLP
          character*20 parallelmode
          integer nsrunning,nswaiting,lastspawn,natomfixed
          integer, dimension(:), allocatable :: corelist ! used to draw random core for activation
          integer, dimension(:), allocatable :: iofbuffer ! used for fma to update ma's ma_spawnlist_os
          logical spawnbuf_buf ! If false -- must wait for message to clear
          logical killbuf_buf ! If false -- must wait for message to clear
          logical, dimension(:), allocatable ::  mkillbuf_buf ! If false -- must wait for message to clear
          logical, dimension(:), allocatable ::  startbuf_buf
          logical execute_serial ! Will turn off parallel spectad behavior if true
          logical luselammps
          logical lsynthattempts

          integer groupID, spawnID, spawnIDlocal, packet_size
          integer world_group,  world_comm,  nprocs_w, rank_w !! All cores
          !integer global_group, global_comm, nprocs_g, rank_g !! ma, fma, and one core for each tad-lp
          integer local_group,  local_comm,  nprocs_l, rank_l !! once core for each replica within a tad-lp
          integer force_group,  force_comm,  nprocs_f, rank_f !! force-call cores used within each replica
          integer comm_everyone, host_info
          integer rank_log ! lables all replicas in order
          !integer current_comm, nprocs, rank
          integer ipr_exit
          real*8 cput_ref ! Reference, cputime

          integer, dimension(:), allocatable :: local_groups
          integer, dimension(:), allocatable :: local_comms
          integer, dimension(:), allocatable :: force_groups
          integer, dimension(:), allocatable :: force_comms
          real*8, dimension(:), allocatable :: prmdtarray ! Used by ParRep Cores
          real*8, dimension(:), allocatable :: prwctarray ! Used by ParRep Cores
          integer done_cnt
          integer :: icurrenteventmax = 0
          integer killcnt
          integer spawnDepth

          integer, dimension(:), allocatable :: spawnrunning !<- index = LP(groupID), value = spawnID
          integer, dimension(:), allocatable :: spawnwaiting !<-value = spawn ID waiting to be  
          integer, dimension(:), allocatable :: spawnwaiting_nmove
          real*8, dimension(:), allocatable :: spawnwaiting_t !<-Tlow waiting to be activated
          real*8, dimension(:), allocatable :: spawnwaitingxyz !<-xyz of waiting spawns
          integer, dimension(:), allocatable :: spawnstatelabel

!          integer, dimension(:), allocatable :: malp_spawnkey !<- translation of spawnIDglobal to spawnIDlocal
          integer, dimension(:), allocatable :: lpma_spawnkey !<- translation of spawnIDlocal to spawnIDglobal
          integer, dimension(:), allocatable :: ma_depthlist
          integer, dimension(:), allocatable :: ma_spawnlist_i
          integer, dimension(:), allocatable :: ma_spawnlist_os
          integer, dimension(:), allocatable :: ma_nstatesskipped
          real*8, dimension(:), allocatable :: ma_spawnlist_t
          real*8, dimension(:), allocatable :: ma_spawnlist_e
          real*8,dimension(:),allocatable :: ma_spawnlist_wct
          real*8, dimension(:), allocatable :: ma_boostdenom_t
          real*8, dimension(:), allocatable :: ma_winningbarrier
          integer, dimension(:), allocatable :: ma_winners
          real*8, dimension(:), allocatable :: ma_probs
          real*8, dimension(:), allocatable :: lp_spawnlist
          logical, dimension(:), allocatable :: lp_spawnlist_l
          integer, dimension(:), allocatable :: isdepspawn
          integer lp_winner, masterID, managerID
          integer fmanagerID, fmanagerCR

          ! ReplicaTAD testing variables...
          real*8 :: t_rep_tot = 0.0  ! Total MD time accumulated by this replica (all states)
          real*8 :: t_mas_tot = 0.0  ! Total MD time accumulated by this lp replica-master (all states)
          integer :: n_rep_trans = 0 ! Number of tranitions made by this replica
          integer :: n_mas_trans = 0 ! Number of transitions made by this lp replica-master

          integer :: kill1buf
          integer, dimension(:), allocatable :: kill2buf
          logical, dimension(:), allocatable :: donebuf

          real*8  previousbest, timechecklast
          logical spawnwrite
          logical activated ! flag to start main tad loop
          logical lfinished, neigh_known

          integer fma_nofficial

          integer trans_req, kill_req
          integer put_req, done_req, iof_req
          integer prtrans_req
          integer, dimension(:), allocatable :: share_req
          integer, dimension(:), allocatable :: mkill_req
          integer, dimension(:), allocatable :: start_req
          integer, dimension(:), allocatable :: deac_req

          integer, PARAMETER :: trans_tag1 = 1

          integer, PARAMETER :: start_tag1 = 2

          integer, PARAMETER :: kill_tag1  = 3
          integer, PARAMETER :: kill_tag2  = 4

          integer, PARAMETER :: done_tag1  = 5
          integer, PARAMETER :: done_tag2  = 6

          integer, PARAMETER :: pstate_tag1  = 7
          integer, PARAMETER :: iof_tag1  = 8

          integer, PARAMETER :: deac_tag1 = 9
          integer, PARAMETER :: deac_tag2 = 10
          integer, PARAMETER :: prcoords_tag1 = 11
          integer, PARAMETER :: prtrans_tag1 = 12
          integer, PARAMETER :: prtrans_tag2 = 13
          integer, PARAMETER :: prtrans_tag3 = 14
          integer, PARAMETER :: prtime_tag1 = 15
          integer, PARAMETER :: prtime_tag2 = 16

          integer :: pstate_tag2  = 100

          integer :: ss_ind = 0
          integer :: bcast_id = 0
          integer :: sharedstate_buf_size
          integer :: bcast_id_max
          real*8, dimension(:), allocatable :: sharedstate_buf
          real*8, dimension(:), allocatable :: putstate_buf
          real*8, dimension(:), allocatable :: spawnmes_buf
          real*8, dimension(:), allocatable :: startmes_buf
          real*8, dimension(:), allocatable :: officialmes_buf
          real*8, dimension(:), allocatable :: prtrans_buf

        ! NOTE That mpirun command should use [packet_size*lp_avail+1] ... (rank 0 is manager)

        !! Define Structure for holding official state information:
          integer, PARAMETER ::  mxatom2 = 2000    ! <- Increase for larger systems
          integer, PARAMETER ::  nneighmax2 = 200  ! <- Increase if more neighbors expected per state
          integer, PARAMETER ::  ntickmax2 = 50
          integer, PARAMETER ::  nattemptkMC = 10  !10
          real*8,  PARAMETER ::  minbarrier_hack = -1000.0d0
          !!integer, PARAMETER ::  nParRep = 1 ! <- # of TAD replica processes doing ParRep-MD
          !!integer, PARAMETER ::  nforcecores = 1 ! <- # of cores to use for force calls
          integer  ::  nParRep = 1 ! <- # of TAD replica processes doing ParRep-MD
          integer  ::  nforcecores = 1 ! <- # of cores to use for force calls
          integer, PARAMETER ::  nforcecoresSpawn = 1 ! <- # of cores to use for force calls
          logical, PARAMETER ::  spawn_lammps = .false. ! Run lammps force call on spawned cores (currently slower than inline cores - if improved, may simplify the code a bit)
          integer :: nstatedata = 1 ! Assume first is current state of LP (no matter what)
          integer :: n_states_skipped
          integer :: istatenow = 1
          TYPE :: StateInfo
            integer       :: nmove !moving atoms can change in a single run (deposition)
            integer       :: istate
            real*8        :: emin
            integer       :: nneigh
            integer       :: nnegmin
            integer       :: nprod
            integer, dimension(:), allocatable :: nattempted
            real*8        :: freqlogsmin
            real*8        :: barrierminev
            real*8        :: xyzmin(3*mxatom2)
            real*8, dimension(:), allocatable :: xyzneighs
            real*8, dimension(:), allocatable :: eneigh
            real*8, dimension(:), allocatable :: barrierev
            real*8, dimension(:), allocatable :: barrevev
            real*8, dimension(:), allocatable :: prefac
            real*8        :: timehighnow
            real*8        :: timehighprev
            real*8        :: timelownow ! Low temperature time accumulated in state (for )

            ! Record the most-stretched bond at the saddle for each neighbor:
            !
            ! -> nneigh_block  = the # of transitions to block during MD
            ! -> ibondneigh1   = the 1st atom ID in the tri-bond (the atom that moved most during the trans)
            ! -> ibondneigh2   = the 2nd atom ID in the tri-bond
            ! -> ibondneigh3   = the 3rd atom ID in the tri-bond
            ! -> dbondneigh12  = the dist between 1st and 2nd atoms at the saddle
            ! -> dbondneigh13  = the dist between 1st and 3rd atoms at the saddle
            ! -> dbondneigh12i = the dist between 1st and 2nd atoms at the min
            ! -> dbondneigh13i = the dist between 1st and 3rd atoms at the min
            !
            ! This information will not be shared with other workers (for now)
            !
            ! These bonds will be checked at each md step if they correspond to
            ! a "low-enough" barrier transition (which is likely to happen during warmup/blackout).
            !
            ! When these transitions are "reflected," they need to be counted as "attempts"
            ! to retain good statistics in synthetic mode...
            !
            integer                            :: nneigh_block
            integer, dimension(:), allocatable :: ibondneigh1
            integer, dimension(:), allocatable :: ibondneigh2
            integer, dimension(:), allocatable :: ibondneigh3
            real*8,  dimension(:), allocatable :: dbondneigh12
            real*8,  dimension(:), allocatable :: dbondneigh13
            real*8,  dimension(:), allocatable :: dbondneigh12i
            real*8,  dimension(:), allocatable :: dbondneigh13i

            integer, dimension(:), allocatable :: neighlabel
          END TYPE StateInfo
          type (StateInfo), dimension(:), allocatable :: statedata

          real*8, dimension(:), allocatable :: ma_timehighnow
          real*8, dimension(:), allocatable :: ma_timelownow

        !%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
          integer, PARAMETER ::  max_spawn = 100000 ! Cant have more than ~10 Milion Spawns in single run
          integer max_depth 
          integer lp_avail
          integer :: nmaxskip = 100000000
          logical :: blocking_on    = .false. ! <- Don't change this
          logical :: allow_blocking = .false.
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
          logical :: dodepositionrun
          integer :: natomstart          ! Set equal to natom at start of run
          integer :: itypestart(mxatom2) ! Set equal to itype at start of run
          real*8  :: tadcharge(mxatom2)
          integer, PARAMETER :: iblockdeposit = 25 !1 ! Set to 1 if inter-deposition event rates are rare
          integer, PARAMETER :: ndepositions = 100 ! Only does depositions if total-time is neg. in input
          integer, PARAMETER :: ndepMDsteps = 2000 !1000 <- Sprague
          real*8 tadtimestopspawn(ndepositions)
          integer itypedeposit(ndepositions)
          real*8 xA_deposit
          real*8 xB_deposit
          real*8 xC_deposit
          real*8 currenttadtimestop
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        !%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!

          !real*8 temphigh_global

       CONTAINS

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine startup_mpi(nofficialmax,irecognize)
         implicit real*8(a-h,o-z)  
         include 'parameters.h'
         integer, dimension(:), allocatable :: ranks
         character*80 str1,str2,str3
        
          ! Start-up MPI...
          ! Construct default WORLD communicator:
          call MPI_Init (ier)
          call MPI_Comm_rank (MPI_COMM_WORLD, rank_w, ier)
          call MPI_Comm_size (MPI_COMM_WORLD, nprocs_w, ier)
          call MPI_Comm_group(MPI_COMM_WORLD, world_group, ier)

          cput_ref = mpi_wtime()-cput_ref
          if(rank_w.eq.0)then
              write(*,*) 'cput_ref:   ',cput_ref
              write(*,*) 'mpi_wtick():',mpi_wtick()
          endif
          write(*,*) 'rank_w - ',rank_w
     +		    ,', nprocs_w - ',nprocs_w
          call MPI_BARRIER(MPI_COMM_WORLD,ier)

          done_cnt = 0
          execute_serial = .false.
          ! If fewer than 3 cores available, cannot do spectad...
          if(nprocs_w.lt.3)then
            stop 'Must use >=3 cores for now!!'
          endif

          !! First 5 lines of input are for SpecTAD (following vars):
          !!     nforcecores, 
          !!     nParRep,
          !!     xA_deposit,
          !!     xB_deposit,
          !!     xC_deposit
          if(rank_w.eq.0) read(5,*) nforcecores
          call MPI_BCAST(nforcecores,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
          if(rank_w.eq.0) read(5,*) nParRep
          call MPI_BCAST(nParRep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
          if(rank_w.eq.0) read(5,*) xA_deposit
          call MPI_BCAST(xA_deposit,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
          if(rank_w.eq.0) read(5,*) xB_deposit
          call MPI_BCAST(xB_deposit,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
          if(rank_w.eq.0) read(5,*) xC_deposit
          call MPI_BCAST(xC_deposit,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)

          lp_avail = (nprocs_w - 2) / (nParRep*nforcecores)
          packet_size = (nprocs_w - 2) / lp_avail
          allocate(ranks(packet_size))
          allocate(local_groups(lp_avail+2))
          allocate(local_comms(lp_avail+2))

          ! Construct local (ParRep) communicators
          do i = 1,lp_avail+2
            if(i.eq.1)then
              ifirst = 0
              ilen = 1
            else
              ifirst = (i-2)*packet_size+1
              ilen = packet_size
              if(i.eq.(lp_avail+2)) ilen = 1
            endif
            ranks(:) = 0
            do j = 1,ilen
              ranks(j) = ifirst + (j-1)
              ilast = ranks(j)
            enddo
            call MPI_Group_Incl(world_group,ilen,ranks,
     +		       local_groups(i),ier)
            call MPI_Comm_create(MPI_COMM_WORLD,local_groups(i),
     +		       local_comms(i),ier)
            if ((rank_w.ge.ifirst).and.(rank_w.le.ilast)) then
              local_comm = local_comms(i)
              local_group = local_groups(i)
              call MPI_Comm_rank(local_comms(i), rank_l, ier)
              call MPI_Comm_size(local_comms(i), nprocs_l, ier)
            endif
            call MPI_BARRIER(MPI_COMM_WORLD,ier)
          enddo

          write(*,*) 'rank_w - ',rank_w
     +		    ,', nprocs_w - ',nprocs_w
     +		    ,', rank_l - ',rank_l
     +		    ,', nprocs_l - ',nprocs_l
     +		    ,', MPI_COMM_WORLD - ',MPI_COMM_WORLD
     +		    ,', local_comm - ',local_comm
          call flushtad(6)
          call MPI_BARRIER(MPI_COMM_WORLD,ier)

          allocate(force_groups(lp_avail*nParRep+2))
          allocate(force_comms(lp_avail*nParRep+2))
          ! Construct force-call communicators
          do i = 1,lp_avail*nParRep+2
            if(i.eq.1)then
              ifirst = 0
              ilen = 1
            else
              ifirst = (i-2)*nforcecores+1
              ilen = nforcecores
              if(i.eq.(lp_avail*nParRep+2)) ilen = 1
            endif
            ranks(:) = 0
            do j = 1,ilen
              ranks(j) = ifirst + (j-1)
              ilast = ranks(j)
            enddo
            call MPI_Group_Incl(world_group,ilen,ranks,
     +		       force_groups(i),ier)
            call MPI_Comm_create(MPI_COMM_WORLD,force_groups(i),
     +		       force_comms(i),ier)
            if ((rank_w.ge.ifirst).and.(rank_w.le.ilast)) then
              force_comm = force_comms(i)
              force_group = force_groups(i)
              call MPI_Comm_rank(force_comms(i), rank_f, ier)
              call MPI_Comm_size(force_comms(i), nprocs_f, ier)
            endif
            call MPI_BARRIER(MPI_COMM_WORLD,ier)
          enddo

          ! Create group labels:
          if(rank_w.eq.(nprocs_w-1))then ! File Manager Core
            groupID = lp_avail + 1 ! very last core
          elseif(rank_w.gt.0)then ! MD Cores
            groupID = (rank_w - 1) / packet_size  + 1
          else ! Spawn Manager Core
            groupID = 0  ! very first core
          endif

          allocate(deac_req(nprocs_l))

          managerID = 0 ! use groupID=0 as spawn-Manager & Last group as state-file Manager
          fmanagerID = lp_avail + 1
          fmanagerCR = nprocs_w - 1

          if(groupID.eq.managerID)then

            allocate(spawnrunning(lp_avail+2)) !<- index = LP(groupID), value = spawnID
            allocate(spawnwaiting(lp_avail+2)) !<-value = spawnID waiting to be activated
            allocate(spawnwaiting_t(lp_avail+2))
            allocate(spawnwaiting_nmove(lp_avail+2)) !<-value = spawnID waiting to be activated

            allocate(spawnstatelabel(max_spawn))
            allocate(kill2buf(lp_avail))
            allocate(mkill_req(lp_avail))
            allocate(start_req(lp_avail))
            allocate(donebuf(lp_avail+1)) ! +1 is due to fmanager being an extra core to kill
            kill2buf(:) = 0
            spawnrunning(:) = 0
            spawnwaiting(:) = 0
            spawnwaiting_t(:) = 0.d0
            spawnstatelabel(:) = 0
            nsrunning = 0
            nswaiting = 0
            allocate(ma_timehighnow(nofficialmax))
            allocate(ma_timelownow(nofficialmax))
            ma_timehighnow(:) = 0.d0
            ma_timelownow(:) = 0.d0

          elseif(groupID.eq.fmanagerID)then

            sharedstate_buf_size =
     +		  1*(7+mxatom2*3+5*nneighmax2+(nneighmax2+1)*mxatom2*3+1)
            bcast_id_max = 100*lp_avail
            if(irecognize.gt.0)then
              allocate(sharedstate_buf(sharedstate_buf_size))
              allocate(share_req(bcast_id_max))
              allocate(statedata(nofficialmax))
              do istate = 1,nofficialmax
                statedata(istate)%timehighnow = 0.d0
                statedata(istate)%timehighprev = 0.d0
                statedata(istate)%timelownow = 0.d0
                call check_statedata(istate,mxatom2,1) ! Allocates Arrays
              enddo
            endif
            fma_nofficial = 1

          else

            if(irecognize.gt.0)then
              allocate(statedata(nofficialmax))
              do istate = 1,nofficialmax
                statedata(istate)%timehighnow = 0.d0
                statedata(istate)%timehighprev = 0.d0
                statedata(istate)%timelownow = 0.d0
                statedata(istate)%nneigh_block = 0
                call check_statedata(istate,mxatom2,1) ! Allocates Arrays
              enddo
            endif

          endif


          write(*,*) 'rank_w - ',rank_w
     +		    ,', nprocs_w - ',nprocs_w
     +		    ,', rank_l - ',rank_l
     +		    ,', nprocs_l - ',nprocs_l
     +		    ,', rank_f - ',rank_f
     +		    ,', nprocs_f - ',nprocs_f
     +		    ,', MPI_COMM_WORLD - ',MPI_COMM_WORLD
     +		    ,', local_comm - ',local_comm
     +		    ,', force_comm - ',force_comm
          call flushtad(6)
          call MPI_BARRIER(MPI_COMM_WORLD,ier)

          filstrtLP = ''
          str2 = ''
          do i = 1, lp_avail
           if(groupID.eq.i)then
             str1 = 'lp.'
             write(str2,*) groupID
             str3 = '.'
             filstrtLP = trim(ADJUSTL(str1))//trim(ADJUSTL(str2))
     +              //trim(ADJUSTL(str3))
           endif
           call MPI_BARRIER(MPI_COMM_WORLD,ier)
          enddo
          filstrt2 = 'sp.' !<- Used To Create filstrt for spawns

          deallocate(ranks)
          return
        end subroutine startup_mpi
        
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !subroutine mpiset_local()
        ! implicit none
        !   if(execute_serial) return
        !   rank = rank_l
        !   nprocs = nprocs_l
        !   current_comm = local_comm
        ! return
        !end subroutine mpiset_local
        
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !subroutine mpiset_global()
        ! implicit none
        !   if(execute_serial) return
        !   rank = rank_g
        !   nprocs = nprocs_g
        !   current_comm = MPI_COMM_WORLD
        ! return
        !end subroutine mpiset_global
        
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine setup_mpi2(nmove,natom,nofficialmax)
         implicit real*8(a-h,o-z)  
         include 'parameters.h'

          if(execute_serial) return

          natomfixed = natom - nmove

          lfinished = .false. ! ALL TADs not done yet...
          activated = .false. ! No processes ready for tad loop yet...
          spawnwrite = .false.
          spawnDepth = 0
          spawnbuf_buf = .false.
          killbuf_buf = .false.

          if (groupID.eq.managerID) then
            allocate(lpma_spawnkey(ntickmax2*lp_avail))
            allocate(ma_depthlist(max_spawn))
            allocate(ma_spawnlist_i(max_spawn))
            allocate(ma_spawnlist_os(max_spawn))
            allocate(ma_nstatesskipped(max_spawn))
            allocate(ma_spawnlist_t(max_spawn))
            allocate(ma_spawnlist_e(max_spawn))
            allocate(ma_spawnlist_wct(max_spawn))
            allocate(ma_boostdenom_t(max_spawn))
            allocate(ma_winningbarrier(max_spawn))
            allocate(ma_winners(max_spawn))
            allocate(spawnwaitingxyz(
     +          (lp_avail+2)*(mxatom2-natomfixed)*3))
            allocate(mkillbuf_buf(lp_avail))
            allocate(startbuf_buf(lp_avail))
            allocate(isdepspawn(max_spawn))
            allocate(startmes_buf(
     +          ((mxatom2-natomfixed)*3+8+2*nofficialmax)*lp_avail))
            allocate(corelist(lp_avail))
            do i = 1, lp_avail
               corelist(i) = i
            enddo

!            malp_spawnkey(:) = 0
            lpma_spawnkey(:) = 0
            ma_depthlist(:) = 0
            ma_spawnlist_i(:) = 0
            ma_spawnlist_os(:) = 0
            ma_nstatesskipped(:) = 0
            ma_spawnlist_t(:) = 0.d0
            ma_spawnlist_e(:) = 0.d0
            ma_spawnlist_wct(:) = 0.d0
            ma_boostdenom_t(:) = 0.d0
            ma_winningbarrier(:) = 0.d0
            ma_winners(:) = 0
            spawnwaitingxyz(:) = 0.0
            lastspawn = 1
            nsrunning = 1
            nwaiting = 0
            spawnrunning(1) = 1
            mkillbuf_buf(:) = .false.
            startbuf_buf(:) = .false.
            startmes_buf(:) = 0.d0
            isdepspawn(:) = 0

          elseif (groupID.eq.fmanagerID) then

            ! Want fmanager to know spawn-state relations:
            allocate(ma_spawnlist_os(max_spawn))
            allocate(officialmes_buf(nofficialmax*2+1))
            ma_spawnlist_os(:) = 0

          else

            allocate(spawnmes_buf(1+3*(mxatom2-natomfixed)+6))
            allocate(lp_spawnlist(ntickmax2))
            allocate(lp_spawnlist_l(ntickmax2))
            lp_spawnlist(:) = 0.0
            lp_spawnlist_l(:) = .false.
            if(groupID.eq.1)then
              activated = .true.
              spawnID = 1
              masterID = 1
            endif

          endif
          
         return
        end subroutine setup_mpi2
        
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine ma_check_spawn(uncertainty)
        
          implicit real*8(a-h,o-z)
          character*80 filn
          character*80 st1,st2,st3,st4,st5,comstr
          integer ii,jj,kk,spawnIDglobal,msgsource,nmove
          integer msgstatus(MPI_STATUS_SIZE)
          logical looking,msgflag,ldontignore
          real*8, dimension(:), allocatable :: spawn_message
          integer messagespawn,winningspawn_g,winningspawn_l

          if(execute_serial) return

          ij = 0
          msgflag = .true.

          ! Check For Spawns to Initiate
          do while (msgflag)

              ij = ij + 1

              ! Check For Spawn Message:
              call MPI_IProbe(MPI_ANY_SOURCE,trans_tag1,MPI_COMM_WORLD,
     +          msgflag,msgstatus,ier)

              if(msgflag)then ! If there is a message waiting...

                !write(*,*) 'Master detects a SPAWN message waiting..'
                call flushtad(6)

                call MPI_Get_count(msgstatus,MPI_REAL8,ibuffsize,ier)
                allocate(spawn_message(ibuffsize))
                msgsource = msgstatus(MPI_SOURCE) ! get the messages in order!
                i = ((msgsource-1)/packet_size)+1

!                write(*,*)
!     +          'Master getting SPAWN message with ibuffsize = '
!     +          ,ibuffsize
!                call flushtad(6)

                ! Recieve Message
                call MPI_recv(spawn_message,ibuffsize,MPI_REAL8,
     +           msgsource,trans_tag1,MPI_COMM_WORLD,msgstatus,ier)

                nmove = NINT(spawn_message(1))
                messagespawn = NINT(spawn_message(1+3*nmove+4))
                istatespawn = NINT(spawn_message(1+3*nmove+5))
                ! Do this correction somewhere later...
                !if(istatespawn.lt.0)then
                !  isneighshare = -istatespawn
                !  istatespawn = 0
                !endif

!                write(*,*) 'spawnID ',messagespawn
!     +           ,' sending BRANCH request!'
!                call flushtad(6)

                !look through ma_winners for parent
                ldontignore = .true.

                if(.true.)then ! DEBUG
                if(ldontignore)then
                 looking = .true.
                 ldontignore = .false.
                 ii = 1
                 jdepth = 0
                 do while (looking)
                  kk = ma_winners(ii)
                  if(ii.eq.messagespawn)then
                    ! Parent is a winner.. proceed
                    ldontignore = .true.
                    exit
                  endif
                  if(kk.gt.0) then
                    ii = kk
                    jdepth = jdepth + 1
                  else
                    looking = .false.
                  endif
                 enddo
                 if(.not.ldontignore)then
                  write(6,'(A,I9,A)')
     +                'Ignoring spawn from parent = '
     +                ,messagespawn
     +                ,' an ancestor has been killed!! '
                  call flushtad(6)
                 endif
                endif
                endif !if(.false.)then

                if(lastspawn.ge.max_spawn) ldontignore=.false.

                ! Done Add to wainting list if cant be a winner
                if(ldontignore)then

                nswaiting = nswaiting + 1
                spawnwaiting(nswaiting) = lastspawn + 1
                spawnwaiting_nmove(nswaiting) = NINT(spawn_message(1))
                spawnwaiting_t(nswaiting) = spawn_message(1+3*nmove+1)
                spawnIDglobal = lastspawn + 1
                spawnstatelabel(spawnIDglobal) = istatespawn
                isdepspawn(spawnIDglobal) = 
     +                                NINT(spawn_message(1+3*nmove+6))

                ! Define Depth of new spawnID
                ma_depthlist(spawnIDglobal) = 
     +           ma_depthlist(messagespawn)+1
                spawnIDlocal = NINT(spawn_message(1+3*nmove+2))
                !if(spawnIDlocal.eq.1)then
                !  do iii = 1, ntickmax2
                !   lpma_spawnkey((iii-1)*lp_avail+i) = 0
                !  enddo
                !endif
                lpma_spawnkey((spawnIDlocal-1)*lp_avail+i)
     +           = spawnIDglobal

                winningspawn_l = NINT(spawn_message(1+3*nmove+3))
                winningspawn_g = lpma_spawnkey((winningspawn_l-1)
     +           *lp_avail+i)

                if(winningspawn_g.le.messagespawn) then
                  write(6,*) '.. spawnIDglobal = ',spawnIDglobal
                  write(6,*) '^Something wrong with this spawn!!'
                  winningspawn_g = 0
                endif

                ! Add new state to the waiting list:
                if(nswaiting.lt.1)then
                  write(*,*) 'SpawnID ',messagespawn,
     +           ' ERROR, nswaiting = ',nswaiting
                endif
                if(mxatom2-natomfixed.lt.1)then
                  write(*,*) 'SpawnID ',messagespawn,
     +           ' ERROR, mxatom2,natomfixed = ',mxatom2,natomfixed
                endif
                call flushtad(6)
                do j = 1, nmove*3
                   spawnwaitingxyz((nswaiting-1)*(mxatom2-natomfixed)
     +               *3+j)=spawn_message(1+j)
                enddo

                write(*,'(A,I9,A,I9,A,I9)')
     +            'SpawnID ',messagespawn,
     +           ' Spawning ',lastspawn + 1,' in istate ',istatespawn
                call flushtad(6)

                lastspawn = lastspawn + 1

                ! Record current winning spawnID of spawnID source
                ma_winners(messagespawn) = winningspawn_g

                endif ! end if(ldontignore)

                deallocate(spawn_message)

              endif
          enddo
          
         return
        end subroutine ma_check_spawn
        
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine ma_check_done()
        
          implicit real*8(a-h,o-z)
          integer lwinner,gwinner,msgsource,icurrentevent
          integer msgstatus(MPI_STATUS_SIZE)
          logical msgflag,update_winlist,looking
          logical update_winfils
          character*80 st1,st2,st3,st4,st5,comstr
          real*8, dimension(:), allocatable :: messagebuffer
          integer*8 jlong

          if(execute_serial) return

          update_winlist = .false.
          update_winfils = .false.

          ij = 0
          msgflag = .true.

          ! Check For Spawns to Initiate
          !do while (msgflag.and.(ij.lt.lp_avail))
          do while (msgflag)

            ij = ij + 1

            if(.true.) then

              ! Check For Done Message:
              call MPI_IProbe(MPI_ANY_SOURCE,done_tag1,MPI_COMM_WORLD,
     +          msgflag,msgstatus,ier)

              if(msgflag)then ! If there is a message waiting...

!               write(*,*) 'Master detects a DONE message waiting..'
!               call flushtad(6)

               call MPI_Get_count
     +           (msgstatus,MPI_REAL8,ibuffsize,ier)
               allocate(messagebuffer(ibuffsize))

                msgsource = msgstatus(MPI_SOURCE) ! get the messages in order!
                i = ((msgsource-1)/packet_size)+1

!               write(*,*) 'spawnID',spawnrunning(i)
!     +           ,' on groupID ',i,' sending DONE message!'
!     +           ,' with ibuffsize = ',ibuffsize
!                call flushtad(6)

               ! Recieve Message
               call MPI_recv(messagebuffer,ibuffsize,
     +          MPI_REAL8,msgsource,
     +          done_tag1,MPI_COMM_WORLD,msgstatus,ier)
               lwinner = NINT(messagebuffer(1))
               ma_spawnlist_t(spawnrunning(i)) = messagebuffer(2)
               ma_boostdenom_t(spawnrunning(i)) = messagebuffer(3)
               ma_spawnlist_wct(spawnrunning(i)) = messagebuffer(4)
               icurrentevent = NINT(messagebuffer(5))
               instatesskipped = NINT(messagebuffer(6))
               ma_nstatesskipped(spawnrunning(i)) = instatesskipped
               istate = NINT(messagebuffer(7))
               ma_winningbarrier(spawnrunning(i)) = messagebuffer(8)
               if(istate.ne.0)then
                 ma_spawnlist_os(spawnrunning(i)) = istate ! official state of spawn
               endif
               ! Update summary files at least every maximum deposition event
               if(icurrentevent.ge.(icurrenteventmax+1))then
                 icurrenteventmax = icurrentevent
                 update_winlist = .true.
                 update_winfils = .true.
               endif
               ma_spawnlist_e(spawnrunning(i)) = messagebuffer(9)

               ! If there is tick info, record it now:
               if(ibuffsize.gt.9)then
                 if(istate.gt.0)then
                  timehighnowold = ma_timehighnow(istate)
                  timelownowold = ma_timelownow(istate)
                  if(messagebuffer(10).ge.(0.d0))then
                    timehighnownew = timehighnowold + messagebuffer(10)
                    timelownownew  = timelownowold  + messagebuffer(11)
                    if(messagebuffer(11).le.(1.d-15))then
                      write(*,*) 'WARNING!!! SpawnID ',spawnrunning(i)
     +               ,' HAS dtimehighnow = ',messagebuffer(10)
     +               ,' paired with dtimelownow ',messagebuffer(11)
                    endif
                  else
                    timehighnownew = 0.d0
                    timelownownew  = 0.d0
                    write(*,*) 'SpawnID ',spawnrunning(i)
     +               ,' RESETTING timehighnow = ',timehighnownew
     +               ,' to MANAGER for istate ',istate
     +               ,' (timehighold = ',timehighnowold,')'
                     call flushtad(6)
                  endif
                  ma_timehighnow(istate) = timehighnownew
                  ma_timelownow(istate)  = timelownownew
                 endif
               endif
               deallocate(messagebuffer)

               done_cnt = done_cnt + 1

               ! Should allow spawnlist.dat update every few mins or so:
               !if((MOD(done_cnt,1000).eq.0).or.(done_cnt.eq.100))then
               if(done_cnt.eq.100)then
                 update_winlist = .true.
                 update_winfils = .true.
               endif
               timecheck = mpi_wtime()-cput_ref
               ! Allow update every minute at most
               if((timecheck-timechecklast).ge.(60.0*1))then
                 update_winlist = .true.
                 update_winfils = .true.
               endif

               write(*,'(A,I9,A,I9)')
     +           'spawnID',spawnrunning(i)
     +           ,'FINISHING with spawnIDlocal winner '
     +           , lwinner
                call flushtad(6)

               nsrunning = nsrunning - 1
               
               write(*,'(A,I9,A,I9,A,I9,A,I9)')
     +           'SpawnID ',spawnrunning(i),
     +           ' Releasing GroupID ',i,
     +           '.. nsrunning = ',nsrunning,
     +           '.. nswaiting = ',nswaiting
                call flushtad(6)
               ma_spawnlist_i(spawnrunning(i)) = -1
               spawnrunning(i) = 0

              endif 
            endif
          enddo

           if(update_winlist)then
             timechecklast = timecheck
             ! Record Curent Winning Spawn List:
             open(UNIT=611,FILE='spawnlist.dat',STATUS='replace')
             looking = .true.
             ii = 1
             jlong = 1
             curprob = 1.d0
             acctime = 0.d0
             wcttime = 0.d0
             accdenom = 0.d0
             do while (looking)
              kk = ma_winners(ii)
              !acctime = acctime + ma_spawnlist_t(ii)
              acctime = ma_spawnlist_t(ii)
              wcttime = max(wcttime,ma_spawnlist_wct(ii))
              accdenom = wcttime*ma_boostdenom_t(ii)
              winbarrier = ma_winningbarrier(ii)
              energymin = ma_spawnlist_e(ii)*27.21d0
              if(accdenom.gt.0d0)then
                  cpuboost = acctime / accdenom
              else
                  cpuboost = 0.d0
              endif
              ispawnstate = ma_spawnlist_os(ii)
              jlong = jlong + ma_nstatesskipped(ii) ! Increase depth by number of "skipped" states
              write(611,'(I12,1X,I10,1X,ES16.6,1X,ES16.6,
     +          1X,ES16.6,1X,I9,1X,ES16.6,1X,ES16.6)')
     +          jlong,ii,acctime,wcttime,cpuboost
     +          ,ispawnstate,energymin,winbarrier
              update_winfils = .false. !! DEBUG -- Always set to FALSE
              if(update_winfils)then
                  st1 = 'cp sp.init.'
                  st2 = ''
                  write(st2,*) ii
                  st3 = '.dat spectad.'
                  st4 = ''
                  write(st4,*) jlong
                  st5 = '.dat'
                  comstr=trim(ADJUSTL(st1))//trim(ADJUSTL(st2))
     +              //trim(ADJUSTL(st3))//trim(ADJUSTL(st4))
     +              //trim(ADJUSTL(st5))
                  call system(comstr)
              endif
              if(kk.gt.0) then
                jlong = jlong + 1
                ii = kk
              else
                looking = .false.
              endif
             enddo
             close(611)

            if(ii.ne.masterID)then
              write(*,'(A,I9,A,I9)')
     +                   'Farthest SpawnID = ',ii,
     +                   ' at depth = ',jlong
              call flushtad(6)
              masterID = ii
            endif

           endif
          
         return
        end subroutine ma_check_done

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine ma_activate_new(nofficialmax)
        
          implicit real*8(a-h,o-z)
          include 'parameters.h'
          character*80 filn
          integer usej,msgdest,lp_ready
          integer msgstatus(MPI_STATUS_SIZE)
          logical looking,ldontignore,testflag

          if(execute_serial) return

          ! Send Activate Messages
          if((lp_avail.gt.nsrunning).and.(nswaiting.gt.0))then

          lp_ready = lp_avail - nsrunning   ! Number of cores that we can activate
          ix = min ( nswaiting , lp_ready ) ! Number of cores that we will activate

          if(ix.eq.0) return

          if(lp_ready.gt.0)then
            ! Shuffle cores:
            do j = 1, lp_avail
               jj = FLOOR(prngen(0)*real(lp_avail,8))+1
               tmp_val = corelist(jj)
               corelist(jj) = corelist(j)
               corelist(j) = tmp_val
            enddo
          endif

          do i = 1, min(ix,1) ! lets only activate 1 spawns between kill-checks.. (wastes fewer cores)
            if(lp_avail.gt.nsrunning)then

              do j = 1, lp_avail
                jj = corelist(j)
                if(spawnrunning(jj).eq.0)then
                  usej = jj
                  exit
                endif
              enddo
              msgdest = (usej-1)*packet_size+1

              if(startbuf_buf(usej))then
                !call MPI_Wait(start_req(usej),msgstatus,ier)
                testflag = .false.
                timenowi = mpi_wtime()-cput_ref
                do while(.not.testflag)
                  call MPI_Test(start_req(usej),testflag,msgstatus,ier)
                  if(.not.testflag) then
                    timenow = mpi_wtime()-cput_ref
                    if((timenow-timenowi).gt.(60.0))then
                      write(*,*) 
     +                  'Manager still waiting on spawn start_req !!!'
                      call flush(6)
                      timenowi = mpi_wtime()-cput_ref
                    endif
                  endif
                enddo
              endif
              startbuf_buf(usej) = .false.

                !look through ma_winners for parent
 512  continue

                 looking = .true.
                 ldontignore = .false.
                 ii = 1
                 do while (looking)
                  kk = ma_winners(ii)
                  if(ii.eq.spawnwaiting(1))then
                    ! Spawn is a current winner.. proceed
                    ldontignore = .true.
                    exit
                  endif
                  if(kk.gt.0) then
                    ii = kk
                  else
                    looking = .false.
                  endif
                 enddo
                 if(.not.ldontignore)then
                  write(6,'(A,I9,A,I9)') 'Ignoring spawn = '
     +                ,spawnwaiting(1)
     +                ,' ... it cannot win !! nswaiting = ',nswaiting
                  call flushtad(6)
                  goto 712
                 endif

              ! Update running state buffers...
              nmove = spawnwaiting_nmove(1)
              nmovemax = (mxatom2-natomfixed)
              startmes_buf((usej-1)*(nmovemax*3+8+2*nofficialmax)+1)
     +                                        = spawnwaiting_nmove(1)
              do k = 1, nmove*3
                startmes_buf((usej-1)*(nmovemax*3+8+2*nofficialmax)
     +                                 +1+k) = spawnwaitingxyz(k)
              enddo
              startmes_buf((usej-1)*(nmovemax*3+8+2*nofficialmax)
     +                         +1+nmove*3+1) = real(spawnwaiting(1),8)
              startmes_buf((usej-1)*(nmovemax*3+8+2*nofficialmax)
     +                         +1+nmove*3+2) = spawnwaiting_t(1)
              startmes_buf((usej-1)*(nmovemax*3+8+2*nofficialmax)
     +                         +1+nmove*3+3) 
     +                         = real(ma_depthlist(spawnwaiting(1)),8)
              startmes_buf((usej-1)*(nmovemax*3+8+2*nofficialmax)
     +                         +1+nmove*3+4) 
     +                         = real(isdepspawn(spawnwaiting(1)),8)

              istatespawn = spawnstatelabel(spawnwaiting(1))
              startmes_buf((usej-1)*(nmovemax*3+8+2*nofficialmax)
     +                         +nmove*3+6) 
     +                         = istatespawn
              imess_size = nmove*3+6
              if(istatespawn.gt.0)then
                imess_size = nmove*3+8
                ma_spawnlist_os(spawnwaiting(1)) = istatespawn
                startmes_buf
     +            ((usej-1)*(nmovemax*3+8+2*nofficialmax)+nmove*3+7) =
     +                                      ma_timehighnow(istatespawn)
                startmes_buf
     +            ((usej-1)*(nmovemax*3+8+2*nofficialmax)+nmove*3+8) =
     +                                      ma_timelownow(istatespawn)
                if(ma_timehighnow(istatespawn).gt.(1.d-11))then
                 do idum = 1,nofficialmax
                  startmes_buf((usej-1)*(nmovemax*3+8+2*nofficialmax)
     +                  +nmove*3+8+idum) = ma_timehighnow(idum)
                  startmes_buf((usej-1)*(nmovemax*3+8+2*nofficialmax)
     +                  +nmove*3+8+idum+nofficialmax) 
     +                  = ma_timelownow(idum)
                 enddo
                 imess_size = imess_size + 2*nofficialmax
                endif
              endif

              ! Send Messages...
              call MPI_ISend(
     +        startmes_buf((usej-1)*(nmovemax*3+8+2*nofficialmax)+1),
     +        imess_size,MPI_REAL8,msgdest,start_tag1,
     +        MPI_COMM_WORLD,start_req(usej),ier)
              startbuf_buf(usej) = .true.

              !Find spawn
              write(*,'(A,I9,A,I9,A,I9,A,I9,A,I9)')
     +           'MANAGER ACTIVATING SPAWN ',spawnwaiting(1),
     +           ' on GROUPID ',usej,
     +           '.. nsrunning = ',nsrunning+1,
     +           '.. nswaiting = ',nswaiting-1,
     +           '.. nmove = ',nmove
              call flushtad(6)

              nsrunning = nsrunning + 1
              ma_spawnlist_i(spawnwaiting(1)) = usej
              spawnrunning(usej) = spawnwaiting(1)

              ! Zero out spawnkey for new spawn...
              do idum = 1, ntickmax2
                lpma_spawnkey((idum-1)*lp_avail+usej) = 0
              enddo

 712  continue ! If sent directly here, we are ignoring something in the waiting list
              do j = 1, nswaiting-1
               do l = 1, (mxatom2-natomfixed)*3
                spawnwaitingxyz((j-1)*(mxatom2-natomfixed)*3+l) =
     +            spawnwaitingxyz(j*(mxatom2-natomfixed)*3+l)
               enddo
               spawnwaiting(j) = spawnwaiting(j+1)
               spawnwaiting_t(j) = spawnwaiting_t(j+1)
               spawnwaiting_nmove(j) = spawnwaiting_nmove(j+1)
              enddo
              do l = 1, (mxatom2-natomfixed)*3
               spawnwaitingxyz((nswaiting-1)*(mxatom2-natomfixed)*3+l)
     +                                                           = 0.d0
              enddo
              spawnwaiting(nswaiting) = 0
              spawnwaiting_t(nswaiting) = 0.d0
              spawnwaiting_nmove(nswaiting) = 0
              nswaiting = nswaiting - 1
              if((.not.ldontignore).and.(nswaiting.ge.1))then
                goto 512
              endif

            endif
          enddo

          endif
          
         return
        end subroutine ma_activate_new
        
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine ma_check_kill()
        
          implicit real*8(a-h,o-z)
          logical inwaiting
          integer ll,jj,ilocalID,msgdest,msgsource
          integer msgstatus(MPI_STATUS_SIZE)
          integer kill_list(max_spawn)
          logical msgflag,looking,testflag

          if(execute_serial) return

          ij = 0
          msgflag = .true.

          ! Check For Spawns to Initiate
          do while (msgflag.and.(ij.lt.lp_avail))

            ij = ij + 1

            if(.true.)then

              ! Check For Kill Message:
              call MPI_IProbe(MPI_ANY_SOURCE,kill_tag1,MPI_COMM_WORLD,
     +          msgflag,msgstatus,ier)

              if(msgflag)then ! If there is a message waiting...

!                write(*,*) 'Master detects a KILL message waiting..'
!                call flushtad(6)

                msgsource = msgstatus(MPI_SOURCE) ! get the messages in order!
                i = (msgsource-1)/packet_size+1

                ! Recieve Message
                call MPI_recv(kill1buf,1,MPI_INTEGER,
     +           msgsource,kill_tag1,MPI_COMM_WORLD,msgstatus,ier)

                if(kill1buf.ge.0)then
                   ilocalID = kill1buf
                   kill1buf = lpma_spawnkey((ilocalID-1)
     +                                       *lp_avail+i)
                   write(*,'(A,I9,A,I9,A,I9)')
     +               ' SpawnID ',spawnrunning(i),
     +               ' Killing spawnIDlocal ',ilocalID,
     +               ' aka spawnID ',kill1buf
                   call flushtad(6)
                endif

                if(    (kill1buf.le.spawnrunning(i))
     +            .and.(kill1buf.gt.0)   )then
                   write(*,*) 'STOP !!! SpawnID ',spawnrunning(i),
     +               ' cannot be killing ',kill1buf,'!!'
                   STOP 'clearly killing wrong spawnID !!!!'
                endif

                ! Find Kill :
                if(kill1buf.ne.0)then

                 !generate list of spawns in branch to kill:
                 kill_list(:) = 0
                 looking = .true.
                 ii = kill1buf
                 nkill = 1
                 jdepth = 0
                 do while (looking)
                  kill_list(nkill) = ii
                  kk = ma_winners(ii)
                  if(kk.gt.0) then
                    ii = kk
                    nkill = nkill + 1
                  else
                    looking = .false.
                  endif
                 enddo

                 do ikill = 1, nkill

                  kill1buf = kill_list(ikill)
                  l = ma_spawnlist_i(kill1buf)
                 
                  if(l.gt.0)then !IF RUNNING
                  
                    msgdest = (l-1)*packet_size+1

                    !write(*,*) ' Manager Killing SpawnID ',kill1buf
                    !call flushtad(6)

                     ! Check That Buffer is Ready for the target LP:
                     if(mkillbuf_buf(i))then
                       !call MPI_Wait(mkill_req(i),msgstatus,ier)
                       testflag = .false.
                       timenowi = mpi_wtime()-cput_ref
                       do while(.not.testflag)
                         call 
     +                     MPI_Test(mkill_req(i),testflag,msgstatus,ier)
                         if(.not.testflag) then
                           timenow = mpi_wtime()-cput_ref
                           if((timenow-timenowi).gt.(60.0))then
                             write(*,*)
     +                    'Manager still waiting on spawn mkill_req !!!'
                             call flush(6)
                             timenowi = mpi_wtime()-cput_ref
                           endif
                        endif
                       enddo
                     endif
                     mkillbuf_buf(i) = .false.

                    kill2buf(i) = kill1buf
                    call MPI_Isend(kill2buf(i),1,MPI_INTEGER,
     +               msgdest,kill_tag2,MPI_COMM_WORLD,mkill_req(i),ier)
                     mkillbuf_buf(i) = .true.
     
                    ! Remove spawn from managers list..
                    ! If there is a new spawn at that LP now, it will resend the kill
                    ma_spawnlist_i(kill1buf) = 0
                   
                  elseif(l.eq.0)then !IF WAITING

                    inwaiting = .false.
!                    write(*,*)'Manager Looking for SpawnID ',kill1buf
!                    call flushtad(6)
                    ! Remove spawn from waitning list if necessary..
                    ix = nswaiting
                    do j = 1, ix
                      if(spawnwaiting(j).eq.kill1buf) then
                       inwaiting = .true.
                       nmove = spawnwaiting_nmove(j)
!                       write(*,*)'Manager Removing SpawnID ',kill1buf
!                       call flushtad(6)
                       do jj = j, nswaiting-1
                        nmove_next = spawnwaiting_nmove(jj+1)
                        do ll = 1, (mxatom2-natomfixed)*3
                         spawnwaitingxyz
     +                     ((jj-1)*(mxatom2-natomfixed)*3+ll) =
     +                     spawnwaitingxyz(jj*(mxatom2-natomfixed)*3+ll)
                        enddo
                        spawnwaiting(jj) = spawnwaiting(jj+1)
                        spawnwaiting_nmove(jj) = nmove_next
                        spawnwaiting_t(jj) = spawnwaiting_t(jj+1)
                       enddo
                       do ll = 1,(mxatom2-natomfixed)*3
                        spawnwaitingxyz
     +                  ((nswaiting-1)*(mxatom2-natomfixed)*3+ll) = 0.d0
                       enddo
                       spawnwaiting(nswaiting) = 0
                       spawnwaiting_nmove(nswaiting) = 0
                       spawnwaiting_t(nswaiting) = 0.d0
                       nswaiting = nswaiting - 1
                      endif
                    enddo
                   
                  else

                    write(*,'(A,I9,A)') 'SpawnID ',kill1buf
     +                           ,' already finished (killed or done)'
                    call flushtad(6)
                   
                  endif ! end waiting/done test

                 enddo ! ikill = 1, nkill

                else ! executed if kill1buf == 0

                  write(*,'(A,I9,A)') 'WARNING.. groupid ',i
     +       ,' told to kill a spawn before it has reached master!!'
                  call flushtad(6)

                endif ! End kill1buf 1,0,(-) check

              endif ! msgflag
            endif ! if(.true.)then
          enddo

          
         return
        end subroutine ma_check_kill

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine ma_check_finish(isleep_ma)
        
          implicit real*8(a-h,o-z)
          integer, dimension(:), allocatable :: winlist
          integer finaldepth, ii, workingdepth, isleep_ma
          integer currentdepth, currentmaster, msgdest
          character*80 st1,st2,st3,st4,st5,comstr
          logical looking
          integer*8 jlong

          if(execute_serial) return
          if(isleep_ma.lt.10)then
            call sleep(1)
            if(isleep_ma.gt.0)then
               write(*,'(A,I9)') 'Manager has isleep_ma = ',isleep_ma
            endif
            return
          endif

          if((nsrunning.le.0).and.(nswaiting.le.0))then

             lfinished = .true.
             write(*,*) 'Manager Exiting!'

            do i = 1, lp_avail
              donebuf(i) = .true.
              msgdest = (i-1)*packet_size+1
              call MPI_ISend(donebuf(i),1,MPI_LOGICAL,msgdest,
     +            done_tag2,MPI_COMM_WORLD,done_req,ier)
            enddo

            open(UNIT=611,FILE='spawnlist.dat',STATUS='replace')
            acctime = 0.d0
            wcttime = 0.d0
            accdenom = 0.d0
            looking = .true.
            !stillgood = .true.
            i = 1
            l = 1
            jlong = 1
            do while (looking)
              k = ma_winners(i)
              st1 = 'cp sp.init.'
              st2 = ''
              write(st2,*) i
              st3 = '.dat spectad.'
              st4 = ''
              write(st4,*) jlong
              st5 = '.dat'
! DEBUG -- Dont make spectad.dat files -- too much memory ?
!              comstr=trim(ADJUSTL(st1))//trim(ADJUSTL(st2))
!     +            //trim(ADJUSTL(st3))//trim(ADJUSTL(st4))
!     +            //trim(ADJUSTL(st5))
!              call system(comstr)
              !acctime = acctime + ma_spawnlist_t(i)
              acctime = ma_spawnlist_t(i)
              wcttime = max(wcttime,ma_spawnlist_wct(i))
              accdenom = wcttime*ma_boostdenom_t(i)
              winbarrier = ma_winningbarrier(i)
              energymin = ma_spawnlist_e(i)*27.21d0
              if(accdenom.gt.0d0)then
                  cpuboost = acctime / accdenom
              else
                  cpuboost = 0.d0
              endif
              ispawnstate = ma_spawnlist_os(i)
              jlong = jlong + ma_nstatesskipped(i) ! Increase depth by number of "skipped" states
              write(611,
     +          '(I12,1X,I10,1X,ES16.6,1X,ES16.6,1X
     +          ,ES16.6,1X,I9,1X,ES16.6,1X,ES16.6)')
     +          jlong,i,acctime,wcttime,cpuboost
     +          ,ispawnstate,energymin,winbarrier
              l = i
              if(k.gt.0) then
                i = k
                jlong = jlong + 1
              else
                looking = .false.
              endif
            enddo

!            ! Make final.dat file:
!            st1 = 'cp spectad.'
!            st2 = ''
!            write(st2,*) jlong
!            st3 = '.dat final.dat'
!            comstr=trim(ADJUSTL(st1))//trim(ADJUSTL(st2))
!     +        //trim(ADJUSTL(st3))
!            call system(comstr)

            write(*,*) 'spawnlist file finalized.'

            i = lp_avail + 1 ! Need To kill file-manager
            donebuf(i) = .true.
            msgdest = (i-1)*packet_size+1
            call MPI_ISend(donebuf(i),1,MPI_LOGICAL,msgdest,
     +           done_tag2,MPI_COMM_WORLD,done_req,ier)

          endif
          
         return
        end subroutine ma_check_finish

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine ma_check_iof()
        
          implicit real*8(a-h,o-z)
          integer msgstatus(MPI_STATUS_SIZE)
          integer iofmessage(2)
          logical msgflag

          if(execute_serial) return

          msgflag = .true.

          do while (msgflag)

            ! Check For Done Message:
            call MPI_IProbe(fmanagerCR,iof_tag1,MPI_COMM_WORLD,
     +          msgflag,msgstatus,ier)

            if(msgflag)then ! If there is a message waiting...
              write(*,*) 'Master detects an IOF message waiting..'
              call flushtad(6)
                ! Recieve Message
              call MPI_recv(iofmessage,2,MPI_INTEGER,fmanagerCR
     +              ,iof_tag1,MPI_COMM_WORLD,msgstatus,ier)
              if(ma_spawnlist_os(iofmessage(1)).ne.0)then
                if(ma_spawnlist_os(iofmessage(1)).ne.iofmessage(2))then
                  write(*,*) 'ERROR!! fmanager says spawnID '
     +              ,iofmessage(1),' was in official state '
     +              ,iofmessage(2),' but manager knows it was in state '
     +              ,ma_spawnlist_os(iofmessage(1))
                  call flushtad(6)
                endif
              else
                ma_spawnlist_os(iofmessage(1)) = iofmessage(2)
                write(*,'(A,I9,A,I9)')
     +               'fmanager says spawnID ',iofmessage(1)
     +              ,' was in official state ',iofmessage(2)
                call flushtad(6)
              endif
            endif

          enddo

         return
        end subroutine ma_check_iof

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine lp_check_start(natom,nmove,xyz,tadtime,itype
     +             ,ldepspawn,iofficialmatch,irecognize,nofficialmax
     +             ,taxes,isneighshare)
        
          implicit real*8(a-h,o-z)  
          include 'parameters.h'
          integer msgstatus(MPI_STATUS_SIZE)
          dimension xyzi(mxatom2*3)
          dimension xyz(mxatom2*3)
          dimension pxyz(mxatom2*3)
          dimension itype(mxatom2)
          dimension taxes(3)
          character*80 filn
          logical msgflag
          integer ldepspawn
          real*8, dimension(:), allocatable :: message_buf

          if(execute_serial) return

          if(rank_l.eq.0)then
            call MPI_IProbe(managerID,start_tag1,MPI_COMM_WORLD,
     +        msgflag,msgstatus,ier)
          endif
          call MPI_BCAST(msgflag,1,MPI_LOGICAL,0,local_comm,ier)
          if(msgflag)then ! If there is a start message waiting...
            ! Recieve Message

             if(rank_l.eq.0)then
               call MPI_Get_count
     +             (msgstatus,MPI_REAL8,ibuffsize,ier)
             endif
             call MPI_BCAST(ibuffsize,1,MPI_INTEGER,0,local_comm,ier)
             allocate(message_buf(ibuffsize))
             if(rank_l.eq.0) call MPI_recv(message_buf,ibuffsize,
     +           MPI_REAL8,managerID,start_tag1,MPI_COMM_WORLD,
     +           msgstatus,ier)
             call MPI_BCAST
     +           (message_buf,ibuffsize,MPI_REAL8,0,local_comm,ier)
             call vecmov(xyz(nmove*3+1),xyzi,natomfixed*3)
             nmove = message_buf(1)
             do j = 1, nmove*3
               xyz(j) = message_buf(1+j)
             enddo
             do j = 1, natomfixed*3
               xyz(nmove*3+j) = xyzi(j)
             enddo
             spawnID = NINT(message_buf(1+nmove*3+1))
             tadtime = message_buf(1+nmove*3+2)
             spawnDepth = NINT(message_buf(1+nmove*3+3))
             ldepspawn = NINT(message_buf(1+nmove*3+4))
             istate = NINT(message_buf(nmove*3+6))
             isneighshare = 0
             if(istate.lt.0)then
               isneighshare = -istate
               istate = 0
             endif

             if(rank_l.lt.nforcecores)then
               if((ibuffsize.gt.(nmove*3+6)).and.(irecognize.gt.0))then
                 iofficialmatch = istate
                 statedata(istate)%timehighnow
     +                                   = message_buf(nmove*3+7)
                 statedata(istate)%timelownow
     +                                   = message_buf(nmove*3+8)

                 if((message_buf(nmove*3+7).gt.(1.d-15)).and.
     +                        (message_buf(nmove*3+8).le.(1.d-15)))then
                   if(rank_l.eq.0) write(*,*)
     +                 'WARNING-A!!! SpawnID ',spawnrunning(i)
     +                 ,' GETTING Thigh = ',message_buf(nmove*3+7)
     +                 ,' paired with Tlow ',message_buf(nmove*3+8)
                 endif

                 if(ibuffsize.gt.(nmove*3+8))then
                   nstaterd = (ibuffsize-(nmove*3+8))/2
                   if(nstaterd.ne.nofficialmax)then
                     if(rank_l.eq.0) write(*,*)
     +                 'spawnID ',spawnID
     +                ,' PROBLEM!! nstaterd.ne.nofficialmax -- '
     +                ,nstaterd,nofficialmax
                   endif
                   do idum = 1, nofficialmax !nstaterd
                     statedata(idum)%timehighnow = 
     +                            message_buf(nmove*3+8+idum+0*nstaterd)
                     statedata(idum)%timelownow =
     +                            message_buf(nmove*3+8+idum+1*nstaterd)

                     if((statedata(idum)%timehighnow.gt.(1.d-15)).and.
     +                     (statedata(idum)%timelownow.le.(1.d-15)))then
                      if(rank_l.eq.0) write(*,*)
     +                  'WARNING-B!!! SpawnID ',spawnrunning(i)
     +                 ,' GETTING Thigh = ',statedata(idum)%timehighnow
     +                 ,' paired with Tlow ',statedata(idum)%timelownow
                     endif
                   enddo
                 endif
               endif
             endif
             call MPI_BARRIER(local_comm,ier)
            if(rank_l.eq.0) deallocate(message_buf)

            natom = nmove + natomfixed
            activated = .true. 
            lp_spawnlist(:) = 0.0
            lp_spawnlist_l(:) = .false.

!            if(rank_l.eq.0)
!     +           write(*,'(A,I9,A,I9,A,I9,A,ES16.6,A,I9,1X,I9,A,I9)')
!     +            'spawnID ',spawnID
!     +           ,' now running at Depth = ',spawnDepth
!     +           ,' on groupID ',groupID
!     +           ,' at t = ',tadtime
!     +           ,' with nmove,natom = ',nmove,natom
!     +           ,' and ldepspawn = ',ldepspawn
!            if(rank_l.eq.0) call flushtad(6)

            ! Make sure the atom types are set correctly:
            natomdep = natom - natomstart
            icnt = 0
            do idum = 1,natom
              if(idum.le.natomdep)then
                !itype(idum) = itypestart(1) ! For now, can only add 1 type of atom during deposition
                itype(idum) = itypedeposit(natomdep-idum+1)
              else
                icnt = icnt + 1
                itype(idum) = itypestart(icnt)
              endif
            enddo

!            write(*,*) 'spawnID ',spawnID
!     +           ,' - natom, natomdep, natomstart = '
!     +           ,natom,natomdep,natomstart

          endif
         call MPI_BARRIER(local_comm,ier)
         return
        end subroutine lp_check_start

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine lp_check_finish()
        
          implicit real*8(a-h,o-z)
          logical msgflag
          integer msgstatus(MPI_STATUS_SIZE)

          if(execute_serial) return
          lfinished = .false.
          if(rank_l.gt.0) return

          call MPI_IProbe(0,done_tag2,MPI_COMM_WORLD,
     +      msgflag,msgstatus,ier)
          if(msgflag)then 
            ! Recieve Message
            call MPI_recv(lfinished,1,MPI_LOGICAL,
     +           0,done_tag2,MPI_COMM_WORLD,msgstatus,ier)
            write(*,*) 'GROUPID',groupID,'RECIEVING EXIT MESSAGE!'
            write(*,*) 'GROUPID',groupID,'lfinished = ', lfinished
          endif
          
         return
        end subroutine lp_check_finish 

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine lp_spec_init()
        
          implicit real*8(a-h,o-z)
          character*80 str2,str3

          if(execute_serial) return
          
          ! Update filstrt to list spawnID
          write(str2,*) spawnID
          str3 = '.'
          filstrt = trim(ADJUSTL(filstrt2))//trim(ADJUSTL(str2))
     +              //trim(ADJUSTL(str3))
      
          lp_winner = 0
          previousbest = 1.d99
          if(rank_l.eq.0) spawnwrite = .true.
          
         return
        end subroutine lp_spec_init

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine lp_spawn_update(ijump,timelowevent,
     x          tadtime,xyz2,nmove,istatespawn,ldepspawn,isneighshare)
        
          implicit real*8(a-h,o-z)
          dimension xyz2(nmove*3)
          integer msgstatus(MPI_STATUS_SIZE)
          logical msgflag, testflag
          integer ldepspawn,isneighshare

          if(execute_serial) return
          
          if(ijump.eq.1)then

            if(timelowevent.lt.previousbest) then
             previousbest = timelowevent

             spawnIDlocal = spawnIDlocal + 1

             ! If NOT first spawn, Check That Buffer is Ready:
             call MPI_BCAST(spawnbuf_buf,1,MPI_LOGICAL,0,force_comm,ier)
             if(spawnbuf_buf)then
              !if(rank_f.eq.0) call MPI_Wait(trans_req,msgstatus,ier)
              testflag = .false.
              timenowi = mpi_wtime()-cput_ref
              do while(.not.testflag)
                if(rank_f.eq.0)
     +              call MPI_Test(trans_req,testflag,msgstatus,ier)
                call MPI_BCAST(testflag,1,MPI_LOGICAL,0,force_comm,ier)
                if(.not.testflag) then
                  timenow = mpi_wtime()-cput_ref
                  if((timenow-timenowi).gt.(15.0))then
                    if(rank_f.eq.0) write(*,'(A,I9,A)')
     +              'spawnID ',spawnID,' still waiting on trans_req !!!'
                    if(rank_f.eq.0) call flush(6)
                    timenowi = mpi_wtime()-cput_ref
                  endif
                  call MPI_BARRIER(force_comm,ier)
                endif
              enddo
             endif
             call MPI_BARRIER(force_comm,ier)
             spawnbuf_buf = .false.

             spawnmes_buf(1) = real(nmove,8)
             do i = 1, nmove*3
               spawnmes_buf(1+i) = xyz2(i)
             enddo
             spawnmes_buf(1+nmove*3+1) = timelowevent+tadtime
             spawnmes_buf(1+nmove*3+2) = real(spawnIDlocal,8)

             lp_spawnlist(spawnIDlocal) = timelowevent+tadtime
             lp_spawnlist_l(spawnIDlocal) = .true.
             lp_winner = spawnIDlocal
             spawnmes_buf(1+nmove*3+3) = real(lp_winner,8)
             spawnmes_buf(1+nmove*3+4) = real(spawnID,8)
             if((istatespawn.eq.0).and.(isneighshare.gt.0))then
               ! If we are spawning to an unknown state..
               !  lets use locality to construct neighbors in new state
               !  (Only if the new spawn will be able to read info for this parent state)
               istatespawn = -isneighshare
               istatespawn = 0 !! RJZ DEBUG HACK
             endif
             spawnmes_buf(1+nmove*3+5) = real(istatespawn,8)
             spawnmes_buf(1+nmove*3+6) = real(ldepspawn,8)

             if(spawnIDlocal.gt.ntickmax2)then
               if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
     x          ,' -- Error spawnIDlocal,ntickmax2 = '
     x          ,spawnIDlocal,ntickmax2
               STOP 'Error spawnIDlocal,ntickmax2'
             endif

             ! SPAWN new branch !!
              if(rank_f.eq.0) write(*,'(A,I9,A,ES16.6,A,I9)')
     x           'spawnID ',spawnID
     x          ,' is branching! timelowevent = '
     x          ,timelowevent,' nmove = ',nmove
              if(rank_f.eq.0) call flushtad(6)
              if(rank_f.eq.0)
     x          call MPI_ISend(spawnmes_buf,1+nmove*3+6,MPI_REAL8,
     x          managerID,trans_tag1,MPI_COMM_WORLD,
     x          trans_req,ier)
              call MPI_BARRIER(force_comm,ier)
              spawnbuf_buf = .true.

             ! KILL any older branches spawned !!
             if(spawnIDlocal.gt.1) then
               ! If NOT first time killing, Check That Buffer is Ready:
               if(killbuf_buf)then
                 !if(rank_l.eq.0) call MPI_Wait(kill_req,msgstatus,ier)
                 testflag = .false.
                 timenowi = mpi_wtime()-cput_ref
                 do while(.not.testflag)
                   if(rank_f.eq.0)
     +                 call MPI_Test(kill_req,testflag,msgstatus,ier)
                   call MPI_BCAST
     +                 (testflag,1,MPI_LOGICAL,0,force_comm,ier)
                   if(.not.testflag) then
                     timenow = mpi_wtime()-cput_ref
                     if((timenow-timenowi).gt.(15.0))then
                       if(rank_f.eq.0) write(*,'(A,I9,A)')
     +                     'spawnID ',spawnID,
     +                     ' still waiting on kill_req !!!'
                       if(rank_f.eq.0) call flush(6)
                       timenowi = mpi_wtime()-cput_ref
                     endif
                     call MPI_BARRIER(force_comm,ier)
                   endif
                 enddo
               endif
               call MPI_BARRIER(force_comm,ier)
               killbuf_buf = .false.
               kill1buf = spawnIDlocal-1
               lp_spawnlist(spawnIDlocal-1) = 0.0
               lp_spawnlist_l(spawnIDlocal-1) = .false.
               if(rank_f.eq.0) call MPI_ISend(kill1buf,1,MPI_INTEGER,
     x            managerID,kill_tag1,MPI_COMM_WORLD,kill_req,ier)
               call MPI_BARRIER(force_comm,ier)
               killbuf_buf = .true.
               killcnt = killcnt + 1
             endif

           endif  
          endif
          
         return
        end subroutine lp_spawn_update

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine lp_done_mess(tadtime,boost_denom,cputimesave
     +          ,istate,lfirstinstate,lstillsynth,dtimelow
     +          ,dThigh,icurrentevent,nneigh,lsendtime,dt,nmd,winboost
     +          ,emin)
        
         implicit real*8(a-h,o-z)
         integer msgstatus(MPI_STATUS_SIZE)
         real*8 messagebuffer(11),emin
         logical lfirstinstate,lstillsynth,lsendtime,testflag

         if(execute_serial) return
         if(rank_l.ge.nforcecores) return

         ! If We have spawned, free buffer for next life:
         if(spawnbuf_buf)then
           !if(rank_f.eq.0) call MPI_Wait(trans_req,msgstatus,ier)
           testflag = .false.
           timenowi = mpi_wtime()-cput_ref
           do while(.not.testflag)
             if(rank_f.eq.0)
     +           call MPI_Test(trans_req,testflag,msgstatus,ier)
             call MPI_BCAST(testflag,1,MPI_LOGICAL,0,force_comm,ier)
             if(.not.testflag) then
               timenow = mpi_wtime()-cput_ref
               if((timenow-timenowi).gt.(15.0))then
                 if(rank_f.eq.0) write(*,'(A,I9,A)')
     +                     'spawnID ',spawnID,
     +                     ' still waiting on trans_req 2 !!!'
                 if(rank_f.eq.0) call flush(6)
                 timenowi = mpi_wtime()-cput_ref
               endif
               call MPI_BARRIER(force_comm,ier)
             endif
           enddo
         endif
         call MPI_BARRIER(force_comm,ier)
         spawnbuf_buf = .false.

         ! If We have killed, free buffer for next life:
         if(killbuf_buf)then
           !if(rank_f.eq.0) call MPI_Wait(kill_req,msgstatus,ier)
           testflag = .false.
           timenowi = mpi_wtime()-cput_ref
           do while(.not.testflag)
             if(rank_f.eq.0)
     +           call MPI_Test(kill_req,testflag,msgstatus,ier)
             call MPI_BCAST(testflag,1,MPI_LOGICAL,0,force_comm,ier)
             if(.not.testflag) then
               timenow = mpi_wtime()-cput_ref
               if((timenow-timenowi).gt.(15.0))then
                 if(rank_f.eq.0) write(*,'(A,I9,A)')
     +                     'spawnID ',spawnID,
     +                     ' still waiting on kill_req 2 !!!'
                 if(rank_f.eq.0) call flush(6)
                 timenowi = mpi_wtime()-cput_ref
               endif
               call MPI_BARRIER(force_comm,ier)
             endif
           enddo
         endif
         call MPI_BARRIER(force_comm,ier)
         killbuf_buf = .false.

         messagebuffer(1) = real(lp_winner,8)
         messagebuffer(2) = tadtime
         messagebuffer(3) = boost_denom
         messagebuffer(4) = cputimesave
         messagebuffer(5) = real(icurrentevent,8)
         messagebuffer(6) = real(n_states_skipped,8)
         if(lfirstinstate)then
           messagebuffer(7) = real(0,8)
         else
           messagebuffer(7) = real(istate,8)
         endif
         messagebuffer(8) = winboost
         messagebuffer(9) = emin
         imess_size = 9

         if(lsendtime)then
           if((dThigh.gt.(dt*real(nmd,8))).and.lstillsynth)then
             imess_size = 11
             messagebuffer(10) = dThigh
             messagebuffer(11) = dtimelow
             if(messagebuffer(11).le.(1.d-15))then
               if(rank_f.eq.0) write(*,*) 'WARNING!!! SpawnID ',spawnID
     +           ,' SENDING dtimehighnow = ',messagebuffer(10)
     +           ,' paired with dtimelownow ',messagebuffer(11)
             endif
           elseif(.not.lstillsynth)then
             imess_size = 10   !! This was set to 8 before??
             messagebuffer(10) = -1.d0
             if(rank_f.eq.0)
     +        write(*,*) 'SpawnID ',spawnID,' resetting thigh to ZERO!!'
             if(rank_f.eq.0) call flushtad(6)
           endif
         endif

         if(rank_f.eq.0) call MPI_Send(messagebuffer,imess_size,
     +      MPI_REAL8,managerID,done_tag1,MPI_COMM_WORLD,ier)
         call MPI_BARRIER(force_comm,ier)
          
         return
        end subroutine lp_done_mess
        
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine lp_check_kill()
        
         implicit real*8(a-h,o-z)
         logical msgflag
         integer msgstatus(MPI_STATUS_SIZE)
         integer spawn_to_kill

         if(execute_serial) return
         if(rank_l.gt.0) return

         call MPI_IProbe(managerID,kill_tag2,MPI_COMM_WORLD,
     +          msgflag,msgstatus,ier)
         if(msgflag)then
           spawnwrite = .false.
           call MPI_recv(spawn_to_kill,1,MPI_INTEGER
     +        ,managerID,kill_tag2,MPI_COMM_WORLD,msgstatus,ier)
           if(SpawnID.eq.spawn_to_kill)then
               activated = .false.
           else
               write(*,*) 'SpawnID ',SpawnID
     +        ,' recieving request to kill spawnid ',spawn_to_kill
               call flushtad(6)
           endif
         endif
          
         return
        end subroutine lp_check_kill

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%   PUT/GET STATE INFO   %%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine fma_state_update(taxes,ntype
     +     ,itype,erecognizecrit,lenxyzneighs
     +     ,ietype,maxtyp,rcut,nofficialmax
     +     ,irotalign,nsteep1,gfac,transcrit
     +     ,itranscrit,gcrit,dvcrit,drcrit,ipr,ifundcuse
     +     ,betalow,avgtadtimestop,lnewneigh,lnewstate)
        
          implicit real*8(a-h,o-z)  
          include 'parameters.h'
          dimension rcut(maxtyp,maxtyp)
          dimension xyzmin(mxatom2*3)
          dimension xyzneighs(lenxyzneighs)
          dimension eneigh(nneighmax)
          dimension barrierev(nneighmax)
          dimension barrevev(nneighmax)
          dimension prefac(nneighmax)
          dimension nattempted(nneighmax)
          dimension itype(mxatom2)
          dimension taxes(3)
          dimension pxyz(mxatom2*3)
          integer gstateid, matchid, msgsource
          integer msgstatus(MPI_STATUS_SIZE)
          logical ismatch, ismatch2, msgflag, sendstate, lsendiof
          logical lnewneigh,lnewstate
          real*8, dimension(:), allocatable :: statebuf
          logical, dimension(:), allocatable :: cpneighs
          character*80 st1,st2,st3,st4,st5,comstr
          character*80 filn

          ij = 0
          msgflag = .true.

          ! Check For Spawns to Initiate
          do while (msgflag.and.(ij.lt.lp_avail))

              ij = ij + 1

              ! Check For Spawn Message:
              call MPI_IProbe(MPI_ANY_SOURCE,pstate_tag1,MPI_COMM_WORLD
     +            ,msgflag,msgstatus,ier)

              if(msgflag)then ! If there is a message waiting...

               allocate(cpneighs(nneighmax))

               msgsource = msgstatus(MPI_SOURCE) ! get the messages in order!
               call MPI_Get_count
     +           (msgstatus,MPI_REAL8,ibuffsize,ier)
               igroup = ((msgsource-1)/packet_size)+1

               write(6,*) 'groupID ',igroup,
     +           ' putting state with ibuffsize = ', ibuffsize
               call flushtad(6)

               allocate(statebuf(ibuffsize))

               ! Recieve Message
               call MPI_recv(statebuf,ibuffsize,MPI_REAL8,
     +          msgsource,pstate_tag1,MPI_COMM_WORLD,msgstatus,ier)

               ispawnid = NINT(statebuf(1))
               lsendiof = .true.
               if(ispawnid.lt.0)then
                 lsendiof = .false.
                 ispawnid = -ispawnid
               endif
               nmove = NINT(statebuf(2))
               natom = nmove + natomfixed
               emin = statebuf(3)
               nneigh = NINT(statebuf(4))
               nnegmin = NINT(statebuf(5))
               nprod = NINT(statebuf(6))
               freqlogsmin = statebuf(7)
               do j = 1, natom*3
                 xyzmin(j) = statebuf(7+j)
               enddo
               do j = 1, nneigh
                 eneigh(j) = statebuf(7+j+natom*3+0*nneigh)
                 barrierev(j) = statebuf(7+j+natom*3+1*nneigh)
                 barrevev(j) = statebuf(7+j+natom*3+2*nneigh)
                 prefac(j) = statebuf(7+j+natom*3+3*nneigh)
                 nattempted(j) = statebuf(7+j+natom*3+4*nneigh)
               enddo
               j_start = (7+natom*3+5*nneigh)
               do j = 1, (nneigh+1)*natom*3
                 xyzneighs(j) = statebuf(j_start+j)
               enddo
               timehighprev = statebuf(j_start+(nneigh+1)*natom*3+1)
               timehighnow = 0.d0 !timehighnow = statebuf(j_start+(nneigh+1)*natom*3+2)
               !ieventlast = statebuf(j_start+(nneigh+1)*natom*3+2)

               deallocate(statebuf)

               ismatch = .false.
               drmax_min = 1.d99
               do j=2,fma_nofficial
                if(nmove.eq.statedata(j)%nmove)then !Only compare if same number of moving atoms...
                  call maxdr_rz(nmove,xyzmin
     +              ,statedata(j)%xyzmin,deltarmax)
                  !write(*,*)'deltarmax,transcrit: ',deltarmax,transcrit
                  if (deltarmax.lt.transcrit) then ! potential match
                       if(deltarmax.lt.drmax_min)then
                          drmax_min = deltarmax
                          matchid = j
                          ismatch = .true.
                       endif
                  endif
                endif
               enddo

                sendstate = .true.
                nneigh_prev = 0
                if(ismatch)then
                 ! We have a match
                 ! Should Not Create New Global State ... But may want to update global state file

                 write(*,*) 'groupID ',igroup,
     +           ' YES putstate match at Manager'
                 call flushtad(6)
                 istate = matchid

                 cpneighs(:) = .false.
                 ncopy = 0
                 sendstate = .false. ! Dont bother resending..

                 ! Update statedata for new neighbors (1):
                 call check_statedata(istate,natom,nneigh)

                 !look for new neighbors:
                 nneigh_prev = statedata(istate)%nneigh
                 do j = 1,nneigh
                   ismatch2 = .false.
                   drmax_min = 1.d99
                   idrmax_min = 0
                   do k = 1,statedata(istate)%nneigh
                     call maxdr_rz(nmove
     +                 ,xyzneighs(j*natom*3+1)
     +                 ,statedata(istate)%xyzneighs(k*natom*3+1)
     +                 ,deltarmax)
                     if (deltarmax.lt.transcrit) then ! potential match
                       if(deltarmax.lt.drmax_min)then
                          drmax_min = deltarmax
                          idrmax_min = k
                          ismatch2 = .true.
                          exit
                       endif
                     endif
                   enddo
                   if(.not.ismatch2)then
                       ncopy = ncopy + 1
                       cpneighs(j) = .true.
                   else
                     ! Still Need To Add Attempt Counts
                     statedata(istate)%nattempted(j) =
     +                   statedata(istate)%nattempted(j)
     +                   + nattempted(k)
                     if(statedata(istate)%nattempted(j)
     +                        .ge.(nattemptkMC))then
                       if((nattempted(k).gt.0).and.(lsynthattempts))then
                         sendstate = .true.
                       endif
                     endif
                   endif
                 enddo

                 statedata(istate)%istate = istate
                 inneigh = statedata(istate)%nneigh
                 statedata(istate)%nneigh = inneigh + ncopy
                 statedata(istate)%timehighprev =
     +               statedata(istate)%timehighprev + timehighprev

                 ! Update statedata for new neighbors (2):
                 call check_statedata(istate,natom,inneigh+ncopy)

                 if(timehighprev.gt.(1.d-14))then
                   if(lsynthattempts) sendstate = .true.
                 endif
                 icnt = 0
                 do j = 1, nneigh
                   if(cpneighs(j))then
                     icnt = icnt + 1
                     statedata(istate)%eneigh(inneigh+icnt) =
     +                    eneigh(j)
                     statedata(istate)%barrierev(inneigh+icnt) =
     +                    barrierev(j)
                     statedata(istate)%barrevev(inneigh+icnt) =
     +                    barrevev(j)
                     statedata(istate)%prefac(inneigh+icnt) =
     +                    prefac(j)
                     statedata(istate)%nattempted(inneigh+icnt) =
     +                    nattempted(j)
                   endif
                 enddo
                 icnt = 0
                 do j = 1, nneigh
                   if(cpneighs(j))then
                   icnt = icnt + 1
                   do k = 1, natom*3
                     statedata(istate)
     +                    %xyzneighs((inneigh+icnt)*natom*3+k) =
     +                    xyzneighs(j*natom*3+k)
                   enddo
                   endif
                 enddo

                 ! Now update neighbor count
                 nneigh = statedata(istate)%nneigh
                 if(ncopy.gt.0)then
                    sendstate = .true.
                 endif

                else

                 ! No match
                 ! Should Create New Global State ...
                 if(fma_nofficial.lt.nofficialmax)then
                   fma_nofficial = fma_nofficial + 1
                   istate = fma_nofficial
                   lnewstate = .true.
                 else
                   write(*,*) 'ERROR -- EXCEEDING nofficialmax !!!'
     +              ,' cant share future state information... '
                   deallocate(cpneighs)
                   CYCLE !! Dont do anything else for this message
                 endif

                 ! Prepare statedata for arrays:
                 call check_statedata(istate,natom,nneigh)

                 write(*,*) 'groupID ',igroup,
     +           ' NO putstate match at Manager.. ',
     +              ' Writing state ',istate,' info with nneigh = '
     +              ,nneigh
                 call flushtad(6)
                 statedata(istate)%nmove = nmove
                 statedata(istate)%istate = istate
                 statedata(istate)%emin = emin
                 statedata(istate)%nneigh = nneigh
                 statedata(istate)%nnegmin = nnegmin
                 statedata(istate)%nprod = nprod
                 statedata(istate)%freqlogsmin = freqlogsmin

                 ! Official Accumulated Simulation Time:
                 statedata(istate)%timehighprev = timehighprev

                 do j = 1, natom*3
                   statedata(istate)%xyzmin(j) =
     +                    xyzmin(j)
                 enddo
                 do j = 1, nneigh
                   statedata(istate)%eneigh(j) =
     +                    eneigh(j)
                   statedata(istate)%barrierev(j) =
     +                    barrierev(j)
                   statedata(istate)%barrevev(j) =
     +                    barrevev(j)
                   statedata(istate)%prefac(j) =
     +                    prefac(j)
                   statedata(istate)%nattempted(j) =
     +                    nattempted(j)
                 enddo
                 do j = 1, (nneigh+1)*natom*3
                   statedata(istate)%xyzneighs(j) =
     +                    xyzneighs(j)
                 enddo

                endif

                if(lsendiof)then ! ispawnid=0 if no message should be sent
                 if(ma_spawnlist_os(ispawnid).le.1)then
                  ma_spawnlist_os(ispawnid) = istate
                  msgdest = 0
                  if(.not.allocated(iofbuffer))then
                    allocate(iofbuffer(2*lp_avail))
                  endif
                  iofbuffer((igroup-1)*2+1) = ispawnid
                  iofbuffer((igroup-1)*2+2) = istate
                  write(*,*) 'spawnID ',ispawnid
     +              ,' will have iof message sent, istate = ',istate
                  call flushtad(6)
                  ! Let Manager know about this relationship
                  call MPI_ISend(iofbuffer((igroup-1)*2+1),2
     +                ,MPI_INTEGER,msgdest,iof_tag1,MPI_COMM_WORLD
     +                ,iof_req,ier)
                 endif
                endif

                !! Send State to LP's...
                if(sendstate)then

                 imessize = 7+natom*3+5*nneigh+(nneigh+1)*natom*3+1

                 ! Dont Bcast the message unless there is more space in the buffer
                 ! -> move back to the start of buffer if possible
                 !if(ss_ind.gt.0)then ! wait every time for now
                 if(((ss_ind+imessize).gt.sharedstate_buf_size)
     +             .or.((bcast_id+lp_avail).ge.(bcast_id_max))
     +             )then
                   write(*,*) 'spawnID ',ispawnid
     +              ,' WARNING: '
     +              ,'ss_ind,imessize,sharedstate_buf_size,bcast_id: '
     +              ,ss_ind,imessize,sharedstate_buf_size,bcast_id
                   call flushtad(6)
                   ! Make sure past messages are up-to-date
                   ! ...(to avoid overwriting info)
                   do j = 1, bcast_id
                     call MPI_Wait(share_req(j),msgstatus,ier)
                   enddo
                   write(*,*) 'spawnID ',ispawnid
     +              ,' messages all up-to-date: continuing... '
                   call flushtad(6)
                   ss_ind = 0
                   bcast_id = 0
                   sharedstate_buf(:) = 0.d0
                 endif

                 !update shared state buffer:
                 sharedstate_buf(ss_ind+1) = real(istate,8)
                 sharedstate_buf(ss_ind+2) = real(nmove,8)
                 sharedstate_buf(ss_ind+3) = emin
                 sharedstate_buf(ss_ind+4) = real(nneigh,8)
                 sharedstate_buf(ss_ind+5) = real(nnegmin,8)
                 sharedstate_buf(ss_ind+6) = real(nprod,8)
                 sharedstate_buf(ss_ind+7) = freqlogsmin

                 do j = 1, natom*3
                   sharedstate_buf(ss_ind+7+j) = 
     +                     statedata(istate)%xyzmin(j)
                 enddo
                 do j = 1, nneigh
                   sharedstate_buf
     +                     (ss_ind+7+natom*3+j+0*nneigh) =
     +                     statedata(istate)%eneigh(j)
                   sharedstate_buf
     +                     (ss_ind+7+natom*3+j+1*nneigh) =
     +                     statedata(istate)%barrierev(j)+bar_add
                   sharedstate_buf
     +                     (ss_ind+7+natom*3+j+2*nneigh) =
     +                     statedata(istate)%barrevev(j)
                   sharedstate_buf
     +                     (ss_ind+7+natom*3+j+3*nneigh) =
     +                     statedata(istate)%prefac(j)
                   sharedstate_buf
     +                     (ss_ind+7+natom*3+j+4*nneigh) =
     +                     statedata(istate)%nattempted(j)
                 enddo
                 j_start = (ss_ind+7+natom*3+5*nneigh)
                 do j = 1, (nneigh+1)*natom*3
                   sharedstate_buf(j_start+j) =
     +                     statedata(istate)%xyzneighs(j)
                 enddo
                 j_start = (j_start+(nneigh+1)*natom*3)
                 sharedstate_buf(j_start+1) =
     +                     statedata(istate)%timehighprev
!                 sharedstate_buf(j_start+2) =
!     +                     statedata(istate)%timehighnow

                 write(6,*) 'FMANAGER bcasting state '
     +             ,istate,' with pstate_tag2 = '
     +             ,pstate_tag2,' and imessize ',imessize,
     +                    ' -- natomfixed = ',natomfixed,
     +                    ' -- nmove,natom = ',nmove,natom,
     +                    ' -- nneigh = ',nneigh
                 call flushtad(6)
                 do j = 1, lp_avail
                   bcast_id = bcast_id + 1
                   msgdest = (j-1)*packet_size+1
                   call MPI_ISend(
     +             sharedstate_buf(ss_ind+1),
     +             imessize,MPI_REAL8,msgdest,pstate_tag2,
     +             MPI_COMM_WORLD,share_req(bcast_id),ier)
                 enddo
                 ss_ind = ss_ind + imessize
                 pstate_tag2 = pstate_tag2 + 1

!                else
!                 write(*,*)'ss_ind,imessize,sharedstate_buf_size: '
!     +             ,ss_ind,imessize,sharedstate_buf_size
!                 call flushtad(6)

               endif !if(sendstate)then

               if(sendstate) lnewneigh = .true.

               ! Prin new neighbor info:
               nneigh = statedata(istate)%nneigh
               nneigh_new = statedata(istate)%nneigh-nneigh_prev
               if(nneigh_new.gt.0)then
                 do idum = nneigh-(nneigh_new-1), nneigh
                   write(6,*) 'FMA: Recording istate',istate,
     +                 ' - ineigh',idum,
     +                 ' - Ea = ',statedata(istate)%barrierev(idum)
                   call flushtad(6)
                 enddo
               endif

               if(.not.ismatch)then
                 !CRZ -- Store New Official state
                 filn=''
                 write(filn,*) istate
                 filn='state.'//trim(adjustl(filn))//'.dat'
                 pxyz(:) = 0.d0
                 ! Make sure the atom types are set correctly:
                 natomdep = natom - natomstart
                 icnt = 0
                 do idum = 1,natom
                   if(idum.le.natomdep)then
                     !itype(idum) = itypestart(1) ! For now, can only add 1 type of atom during deposition
                     itype(idum) = itypedeposit(natomdep-idum+1)
                   else
                     icnt = icnt + 1
                     itype(idum) = itypestart(icnt)
                   endif
                 enddo
!                 write(6,*) 'MANAGER writing state file -- (state.'
!     +             ,istate,'.dat)'
                 call storefile(natom,xyzmin,pxyz,itype,taxes,filn)
               endif

               deallocate(cpneighs)
              endif ! msgflag
         enddo

         return
        end subroutine fma_state_update

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine lp_putstate(
     x  istate,
     x  natom,
     x  emin,
     x  nneigh,
     x  nnegmin,
     x  nprod,
     x  nneighmax,
     x  mxatom,
     x  lenxyzneighs,
     x  xyzmin,
     x  xyzneighs,
     x  eneigh,
     x  barrierev,
     x  barrevev,
     x  freqlogsmin,
     x  prefac,
     x  timehighnow,
     x  timehighprev,
     x  timelowprev,
     x  tickprev,
     x  ieventlast,
     x  lfirstinstate,
     x  lputyes,
     x  lendput,
     x  nattempted,
     x  nofficialmax
     x  )

        implicit real*8(a-h,o-z)
        dimension xyzmin(3*mxatom) ! called xyz1 in tad2.f
        dimension xyzneighs(lenxyzneighs)
        dimension eneigh(nneighmax)
        dimension barrierev(nneighmax)
        dimension barrevev(nneighmax)
        dimension prefac(nneighmax)
        dimension nattempted(nneighmax)
        integer msgdest

        integer msgstatus(MPI_STATUS_SIZE)
        logical lfirstinstate,lputyes,lendput,testflag

          if(execute_serial) return

          if(.not.lputyes)then
            return ! Done if lputyes = .false. Unless we have accumulated time and lsynthattempts
          endif

            msgdest = fmanagerCR

            if(allocated(putstate_buf)) then
              !if(rank_f.eq.0) call MPI_Wait(put_req,msgstatus,ier)
              testflag = .false.
              timenowi = mpi_wtime()-cput_ref
              do while(.not.testflag)
                if(rank_f.eq.0)
     +              call MPI_Test(put_req,testflag,msgstatus,ier)
                call MPI_BCAST(testflag,1,MPI_LOGICAL,0,force_comm,ier)
                if(.not.testflag) then
                  timenow = mpi_wtime()-cput_ref
                  call lp_update_official(nofficialmax)
                  if((timenow-timenowi).gt.(60.0))then
                    if(rank_f.eq.0) write(*,*)
     +                     'spawnID ',spawnID,
     +                     ' still waiting on put_req !!!'
                    if(rank_f.eq.0) call flush(6)
                    timenowi = mpi_wtime()-cput_ref
                  endif
                  call MPI_BARRIER(force_comm,ier)
                endif
              enddo
              call MPI_BARRIER(force_comm,ier)
              deallocate(putstate_buf)
            endif
            allocate(putstate_buf(7+natom*3+5*nneigh+
     +                                (nneigh+1)*natom*3+2))

            imessize = 7+natom*3+5*nneigh+(nneigh+1)*natom*3+1

            if(lfirstinstate)then
              putstate_buf(1) = real(spawnID,8)
            else
              putstate_buf(1) = real(-spawnID,8) ! manager wont need a message from fmanager..
            endif
            putstate_buf(2) = real(natom-natomfixed,8) ! Gives nmove
            putstate_buf(3) = emin
            putstate_buf(4) = real(nneigh,8)
            putstate_buf(5) = real(nnegmin,8)
            putstate_buf(6) = real(nprod,8)
            putstate_buf(7) = freqlogsmin
            do j = 1, natom*3
              putstate_buf(7+j) = xyzmin(j)
            enddo
            do j = 1, nneigh
              putstate_buf(7+j+natom*3+0*nneigh)
     +           = eneigh(j)
              putstate_buf(7+j+natom*3+1*nneigh)
     +           = barrierev(j)
              putstate_buf(7+j+natom*3+2*nneigh)
     +           = barrevev(j)
              putstate_buf(7+j+natom*3+3*nneigh)
     +           = prefac(j)
              putstate_buf(7+j+natom*3+4*nneigh)
     +           = nattempted(j)
            enddo
            j_start = (7+natom*3+5*nneigh)
            do j = 1, (nneigh+1)*natom*3
              putstate_buf(j_start+j) = xyzneighs(j)
            enddo

            ! DEBUG -- this might make it tough to zero-out thigh from 'ma'
            putstate_buf(j_start+(nneigh+1)*natom*3
     +           +1) = timehighprev
!            putstate_buf(j_start+(nneigh+1)*natom*3
!     +           +2) = timehighnow
!            putstate_buf(j_start+(nneigh+1)*natom*3
!     +           +2) = ieventlast

            if(rank_f.eq.0) call MPI_ISend(putstate_buf,imessize,
     +           MPI_REAL8,msgdest,pstate_tag1,MPI_COMM_WORLD,
     +           put_req,ier)
            call MPI_BARRIER(force_comm,ier)

        return
        end subroutine lp_putstate

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine lp_getstate(istate,natom,emin,xyzmin
     +       ,nneigh,nneighmax,mxatom,lenxyzneighs,xyzneighs
     +       ,eneigh,barrierev,barrevev,prefac,barrierminev
     +       ,nofficialmax,freqlogsmin,nnegmin,nprod
     +       ,timehighnow,timehighprev,timelowprev,tickprev
     +       ,ieventlast,betalow,betahigh,fstar,dt,nmd,itype,taxes)

        implicit real*8(a-h,o-z)

        dimension xyzmin(3*mxatom) ! called xyz1 in tad2.f
        dimension xyzneighs(lenxyzneighs)
        dimension eneigh(nneighmax)
        dimension barrierev(nneighmax)
        dimension barrevev(nneighmax)
        dimension prefac(nneighmax)
        dimension pxyz(mxatom2*3)
        dimension itype(mxatom2) ! not dimensioned previously: 01/21/04
        dimension taxes(3)
        character*80 filnam

        if(execute_serial) return

        if(rank_f.eq.0) write(*,*) 'SpawnID ',spawnID,
     +             ' recognizes it is at state ',istate
        if(rank_f.eq.0) call flushtad(6)

        emin = statedata(istate)%emin
        nneigh = statedata(istate)%nneigh
        nnegmin = statedata(istate)%nnegmin
        nprod = statedata(istate)%nprod
        freqlogsmin = statedata(istate)%freqlogsmin
        timehighnow = statedata(istate)%timehighnow
        timelowprev = statedata(istate)%timelownow
        tickprev = timelowprev
        ieventlast = 0
        do j = 1, natom*3
          xyzmin(j) = statedata(istate)%xyzmin(j)
        enddo

        if(nneigh.lt.1)then
          if(rank_f.eq.0) write(*,*) 'SpawnID ',spawnID
     +     ,' ERROR !! in lp_getstate nneigh = ',nneigh
     +     ,' for istate = ',istate
          if(rank_f.eq.0) call flushtad(6)
           filnam=trim(filstrtLP)//'.geterror.dat'
           pxyz(:) = 0.d0
           call storefile(natom,xyzmin,pxyz,itype,taxes,filnam)
        endif

        do j = 1, nneigh
          eneigh(j) = statedata(istate)%eneigh(j)
          barrierev(j) = statedata(istate)%barrierev(j)
          barrevev(j) = statedata(istate)%barrevev(j)
          prefac(j) = statedata(istate)%prefac(j)
        enddo
        barrierminev = statedata(istate)%barrierminev
        do j = 1, (nneigh+1)*natom*3
          xyzneighs(j) = statedata(istate)%xyzneighs(j)
        enddo

        if(timehighnow.ge.(dt*real(nmd,8)))then
          if(rank_f.eq.0) write(*,*) 'SpawnID ',spawnID,
     +             ' has thigh time to use!! timehighnow = '
     +             ,timehighnow,' and tickprev = ',tickprev
          if(rank_f.eq.0) call flushtad(6)
        endif

         call MPI_BARRIER(force_comm,ier)
         return
        end subroutine lp_getstate

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine lp_synthesize_neighbors(natom,emin,xyzmin
     +       ,nneigh,nneighmax,mxatom,lenxyzneighs,xyzneighs
     +       ,eneigh,barrierev,barrevev,prefac,barrierminev
     +       ,nofficialmax,freqlogsmin,nnegmin,nprod
     +       ,timehighnow,betalow,betahigh
     +       ,ineighdep,currenttadtimestop,transcrit,tadtime
     +       ,ntick,ticknext,xmintimelowdrawn
     +       ,jneighlow,isynthclass,isneighshare)

        implicit real*8(a-h,o-z)

        dimension xyzmin(3*mxatom) ! called xyz1 in tad2.f
        dimension xyzneighs(lenxyzneighs)
        dimension eneigh(nneighmax)
        dimension barrierev(nneighmax)
        dimension barrevev(nneighmax)
        dimension prefac(nneighmax)
        dimension ntick(nneighmax)
        dimension isynthclass(nneighmax)
        dimension ticknext(ntickmax2,nneighmax)
        dimension nmovedlist(mxatom)

        tooclose = 5.0 !
        nneigh_init = nneigh
        xmintimelowdrawn = 1.d99
        xminRATElowdrawn = 1.d99
        xminPREFlowdrawn = 1.d99
        xminBARRlowdrawn = 1.d99
        ntick(:) = 0
        ticknext(1,:) = 0.d0
        isynthclass(:) = 0
        jneighlow = 0

         ! Look through isneighshare information for reasonable neighbors
         nmove_i = statedata(isneighshare)%nmove
         nneigh_i = statedata(isneighshare)%nneigh

         if(nmove_i.ne.nmove) goto 711

         if(rank_f.eq.0) write(*,*) 'SpawnID ',spawnID,
     +     ' will try to synthesize new neighbors from ',isneighshare

         ! Determine moving atoms for accepted transition:
         nmoved = 0
         nmovedlist(:) = 0
         do i = 1, nmove
           dx = statedata(isneighshare)%xyzmin((i-1)*3+1)
     +         - xyzmin((i-1)*3+1)
           dy = statedata(isneighshare)%xyzmin((i-1)*3+2)
     +         - xyzmin((i-1)*3+2)
           dz = statedata(isneighshare)%xyzmin((i-1)*3+3)
     +         - xyzmin((i-1)*3+3)
           dist = sqrt(dx**2+dy**2+dz**2)
           if(dist.gt.transcrit)then
             nmoved = nmoved + 1
             nmovedlist(nmoved) = i
           endif
         enddo

         ! Find other neighbors with spatially distinct moving atoms:
         do 599 j=(1+ineighdep),nneigh_i
           do i = 1, nmove
             dx = statedata(isneighshare)%xyzneighs(j*natom*3+(i-1)*3+1)
     +           - statedata(isneighshare)%xyzmin((i-1)*3+1)
             dy = statedata(isneighshare)%xyzneighs(j*natom*3+(i-1)*3+2)
     +           - statedata(isneighshare)%xyzmin((i-1)*3+2)
             dz = statedata(isneighshare)%xyzneighs(j*natom*3+(i-1)*3+3)
     +           - statedata(isneighshare)%xyzmin((i-1)*3+3)
             dist = sqrt(dx**2+dy**2+dz**2)
             if(dist.gt.transcrit)then
               ! Check if starting position of moved atom in neighbor is too close to
               ! the final position of any moved atoms in the accepted neighbor
               do k = 1, nmoved
                 dx = xyzmin((i-1)*3+1) - xyzmin((nmovedlist(k)-1)*3+1)
                 dy = xyzmin((i-1)*3+2) - xyzmin((nmovedlist(k)-1)*3+2)
                 dz = xyzmin((i-1)*3+3) - xyzmin((nmovedlist(k)-1)*3+3)
                 dist = sqrt(dx**2+dy**2+dz**2)
                 if(dist.le.tooclose)then
                   goto 599
                 endif
               enddo
             endif
           enddo

           ! If you get here, you can add synthetic neighbor to our list
           nneigh = nneigh + 1
           eneigh(nneigh) = statedata(isneighshare)%eneigh(nneigh)  !! THIS IS WRONG.
           barrierev(nneigh) = statedata(isneighshare)%barrierev(nneigh)
           barrevev(nneigh) = statedata(isneighshare)%barrevev(nneigh)
           prefac(nneigh) = statedata(isneighshare)%prefac(nneigh)
           call vecmov(statedata(isneighshare)%xyzneighs(j*natom*3+1)
     +         ,xyzneighs(nneigh*natom*3+1),3*natom)
           do k = 1, nmoved
             xyzneighs(nneigh*natom*3+(nmovedlist(k)-1)*3+1) 
     +           = xyzmin((nmovedlist(k)-1)*3+1)
             xyzneighs(nneigh*natom*3+(nmovedlist(k)-1)*3+2)
     +           = xyzmin((nmovedlist(k)-1)*3+2)
             xyzneighs(nneigh*natom*3+(nmovedlist(k)-1)*3+3)
     +           = xyzmin((nmovedlist(k)-1)*3+3)
           enddo

 599  continue

         if(nneigh-nneigh_init.gt.0)then
           if(rank_f.eq.0) write(*,*) 'SpawnID ',spawnID,
     +        ' just synthesized ',nneigh-nneigh_init,' new neighbors'
         endif

         ! Draw times from synthesized neighbors
         xmintimelowdrawn = 1.d99
         xminRATElowdrawn = 1.d99
         xminPREFlowdrawn = 1.d99
         xminBARRlowdrawn = 1.d99
         ntick(:) = 0
         ticknext(1,:) = 0.d0
         isynthclass(:) = 0
         jneighlow = 0
         do jneigh=(1+ineighdep),nneigh
            isynthclass(jneigh) = 1
            ratelow = prefac(jneigh)*exp(-barrierev(jneigh)*betalow)
            if(rank_f.eq.0)then
              timelowdrawn = -(1.0d0/ratelow)*log(1.0d0-prngen(0))
            endif
            call MPI_Bcast(timelowdrawn,1,MPI_REAL8,0,force_comm,ier)
            ticknext(1,jneigh) = timelowdrawn
            ntick(jneigh) = 1
            if(xmintimelowdrawn.gt.timelowdrawn)then
              xmintimelowdrawn = timelowdrawn
              xminRATElowdrawn = ratelow
              xminPREFlowdrawn = prefac(jneigh)
              xminBARRlowdrawn = barrierev(jneigh)
              jneighlow = jneigh
            endif
         enddo

 711  continue

         if(ineighdep.eq.1) then
           ntick(ineighdep) = 1
           ticknext(1,ineighdep) = 1.d30
           prefac(ineighdep) = 1.d-30
           barrierev(ineighdep) = 1.d0
           barrevev(ineighdep) = 1.d0
           if(nneigh.eq.0)then
             nneigh = 1
             do j = 1, natom*3
               if(j.le.nmove)then
                 xyzneighs(ineighdep*natom*3+j) = 0.d0
               else
                 xyzneighs(ineighdep*natom*3+j) = xyzmin(j)
               endif
             enddo
           endif
         endif

         return
        end subroutine lp_synthesize_neighbors

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine lp_kMC(istate,natom,emin,xyzmin
     +       ,nneigh,nneighmax,mxatom,lenxyzneighs,xyzneighs
     +       ,eneigh,barrierev,barrevev,prefac,barrierminev
     +       ,nofficialmax,freqlogsmin,nnegmin,nprod
     +       ,timehighnow,timehighprev,timelowprev,tickprev
     +       ,betalow,betahigh,fstar,dt,nmd
     +       ,ineighdep,currenttadtimestop,transcrit,tadtime
     +       ,ntick,ticknext,xmintimelowdrawn
     +       ,jneighlow,isynthclass,nattempted,uncertainty)

        implicit real*8(a-h,o-z)

        dimension xyzmin(3*mxatom) ! called xyz1 in tad2.f
        dimension xyzneighs(lenxyzneighs)
        dimension eneigh(nneighmax)
        dimension barrierev(nneighmax)
        dimension barrevev(nneighmax)
        dimension prefac(nneighmax)
        dimension nattempted(nneighmax)
        dimension ntick(nneighmax)
        dimension isynthclass(nneighmax)
        dimension ticknext(ntickmax2,nneighmax)
        logical   lkMC,lzerowork

        lkMC = .true.
        istate_initial = istate
        tadtime_initial = tadtime
        n_states_skipped = 0
        jneighlow = 0
        nmove = natom - natomfixed
        xmintimelowdrawn = 1.d99
        xminRATElowdrawn = 1.d99
        xminPREFlowdrawn = 1.d99
        xminBARRlowdrawn = 1.d99
        ntick(:) = 0
        ticknext(:,:) = 0.d0
        isynthclass(:) = 0
        do while ((lkMC).and.(istate.ge.2))

             ! Load State Information
             nmove = statedata(istate)%nmove
             nneigh = statedata(istate)%nneigh

             ! Draw times from neighbors
             xmintimelowdrawn = 1.d99
             xminRATElowdrawn = 1.d99
             xminPREFlowdrawn = 1.d99
             xminBARRlowdrawn = 1.d99
             ntick(:) = 0
             ticknext(1,:) = 0.d0
             isynthclass(:) = 0
             jneighlow = 0
             timehighnow = statedata(istate)%timehighnow
             timehighprev = statedata(istate)%timehighprev
             tickprev = statedata(istate)%timelownow
             timelowprev = tickprev

             ! If we are already sweeping...
             ! Determine if we need more synthetic tick marks:
             nattemptkMC_use = nattemptkMC
             if(tickprev.gt.dt)then
               do jneigh=(1+ineighdep),nneigh
                 nattemptedx = statedata(istate)%nattempted(jneigh)
                 if((lsynthattempts).and.
     +                           (nattemptedx.ge.nattemptkMC))then
                    ! If we have started to sweep in this state
                    ! -> Need to draw synthetic ticks for ALL known states
                    nattemptkMC_use = 1
                 endif
               enddo
             endif

             ! Loop through neighbors again to get tick marks:
             do jneigh=(1+ineighdep),nneigh
                nattemptedx = statedata(istate)%nattempted(jneigh)
                barrier = statedata(istate)%barrierev(jneigh)
                if(lsynthattempts)then
                  if(nattemptedx.ge.nattemptkMC_use)then
                    ! rjz 6/10/2016 -- This line is "timetotx=timehighprev" .. correct?
                    ! --- It seems that this is correct for SpecTAD, bc inside statedata
                    ! --- timehighprev includes all previous time
                    timetotx=timehighprev !timehighnow+timehighprev
                    if(timetotx.lt.(dt*nmd))then !1.d-14))then
                      if(rank_f.eq.0) 
     +                  write(*,'(A,I9,A,A,ES16.6,A,ES16.6,A,I9,1X,I9)')
     +                   'spawnID ',spawnID,' ERROR -- '
     +                  ,' timehighnow = ',timehighnow
     +                  ,' timehighprev = ',timehighprev
     +                  ,' nattempted,jneigh = ',nattemptedx,jneigh
                      call MPI_BARRIER(force_comm,ier)
                      cycle
                    endif
                    ratehigh=float(nattemptedx)/timetotx
                    statedata(istate)%prefac(jneigh)
     +                       =ratehigh/exp(-barrier*betahigh)
                    if((n_states_skipped.eq.0).and.(rank_f.eq.0))
     +                  write(*,
     +                   '(A,I9,A,A,I9,A,A,ES16.6,A,ES16.6,A,I9,1X,I9)')
     +                   'spawnID ',spawnID,' UPDATING PREFACTOR '
     +                  ,' for istate = ',istate,' -- '
     +                  ,' timehighnow = ',timehighnow
     +                  ,' timehighprev = ',timehighprev
     +                  ,' nattempted,jneigh = ',nattemptedx,jneigh
                    call MPI_BARRIER(force_comm,ier)
                  else
                    if(statedata(istate)%prefac(jneigh).lt.(1.0))then
                      cycle ! If we have still a reasonable prefactor,
                            ! we must have synthesized this neighbor earlier
                    endif
                    cycle
                  endif
                endif
                isynthclass(jneigh) = 1
                pre = statedata(istate)%prefac(jneigh)
                ratelow = pre * exp(-barrier*betalow)
                if(lsynthattempts) ratelowatt = ratelow
                if(rank_f.eq.0)then
                  timelowdrawn = -(1.0d0/ratelow)*log(1.0d0-prngen(0))
                endif
                call 
     +            MPI_Bcast(timelowdrawn,1,MPI_REAL8,0,force_comm,ier)
                ticknext(1,jneigh) = timelowdrawn + tickprev
                ntick(jneigh) = 1
                if(xmintimelowdrawn.gt.timelowdrawn)then
                  xmintimelowdrawn = timelowdrawn
                  xminRATElowdrawn = ratelow
                  xminPREFlowdrawn = pre
                  xminBARRlowdrawn = barrier
                  jneighlow = jneigh
                endif
             enddo

             ! Find label of winning neighbor:
             if(jneighlow.gt.0) then
                 istate_next = statedata(istate)%neighlabel(jneighlow)
                 if(istate_next.eq.0)then 
                   neighlabel = 0
                   drmax_min = 1.d99
                   do jstate = 2, nstatedata
                   if ((nmove.eq.statedata(jstate)%nmove)
     +               .and.(istate.ne.jstate))then
                     call maxdr_rz(nmove
     +                 ,statedata(istate)%xyzneighs(3*natom*jneighlow+1)
     +                 ,statedata(jstate)%xyzmin,drmax)
                     if (drmax.lt.transcrit) then ! potential match
                       if(drmax.lt.drmax_min)then
                         drmax_min = drmax
                         neighlabel = jstate
                       endif
                     endif
                   endif
                   enddo
                   if(neighlabel.gt.0)then
                     statedata(istate)%neighlabel(jneighlow)
     +                                                 = neighlabel
                     istate_next = neighlabel
                   endif
                 endif

               ! Does the neighbor label make sense?...
               ! -> If this does happen...
               ! -> Perhaps your 'transcrit' is too small?
               if(istate_next.eq.istate)then
                 if(rank_f.eq.0)then
                   write(*,*) 'B - ERROR at SpawnID',SpawnID,
     x             'trying to spawn same state: istate=',istate,
     x             '- istate_next=',istate_next
                   call flush(6)
                 endif
                 call MPI_BARRIER(force_comm,ier)
                 istate_next = 0
                 statedata(istate)%neighlabel(jneighlow) = 0
               endif

             endif

             if((n_states_skipped.eq.0).and.(istate.ge.2)
     x                                      .and.(jneighlow.gt.0))then
               if(rank_f.eq.0)  
     x         write(6,
     x          '(A,I9,A,ES16.6,A,I9,A,I9,A,I9,A,ES16.6,A
     x          ,ES16.6,A,ES16.6,A,ES16.6,A)')
     x          'SpawnID ',spawnID,' -- xmintimelowdrawn = '
     x         ,xmintimelowdrawn,' in iofficial ',istate
     x         ,' and jneighlow = ',jneighlow,' of ',nneigh
     x         ,' ( prefactor: ',xminPREFlowdrawn
     x         ,' -- rate: ',xminRATElowdrawn
     x         ,' -- barrier: ',xminBARRlowdrawn,' ) - A'
               call MPI_BARRIER(force_comm,ier)
             endif

             ! Skip to Next State?
             lkMC = .false.
             if((jneighlow.gt.0).and.(timehighnow.gt.dt))then
               lzerowork = .false.
               timehighstop=(1.d0/fstar)*(fstar*ticknext(1,jneighlow))
     x                                 **(betahigh/betalow)
               dthigh_tmp = timehighstop-timehighnow
               if(dthigh_tmp.le.(0.d0))then
                 lzerowork = .true.
               endif

               !call MPI_Bcast(lzerowork,1,MPI_LOGICAL,0,force_comm,ier)
               call MPI_BARRIER(force_comm,ier)
               if((lzerowork).and.((tadtime+xmintimelowdrawn)
     x                                     .lt.currenttadtimestop))then
                 if((istate_next.gt.0)
     x             .and.(n_states_skipped.lt.nmaxskip))then
                     lkMC = .true.
                     istate = istate_next
                     tadtime = tadtime + xmintimelowdrawn
                     n_states_skipped = n_states_skipped + 1
                 else if (istate_next.gt.0) then
                   nmaxskip = nmaxskip ! <- Can increase here if you want.
                 endif
               elseif(lzerowork)then
                  if(rank_f.eq.0) write(*,'(A,I9,A)')
     x                 'spawnID ', spawnID
     x                ,' hitting currenttadtimestop in kMC mode '
                  call MPI_BARRIER(force_comm,ier)
               endif
             endif

        enddo !lkMC

         if((istate.ge.2).and.(rank_l.eq.0))then
           ntick(nneigh+1:nneighmax2) = 0
         endif

         if((istate.ge.2).and.(n_states_skipped.gt.0)
     x                   .and.(jneighlow.gt.0))then
           if(rank_f.eq.0) 
     x         write(6,
     x          '(A,I9,A,ES16.6,A,I9,A,I9,A,I9,A,ES16.6,A
     x          ,ES16.6,A,ES16.6,A,ES16.6,A)')
     x          'SpawnID ',spawnID,' -- xmintimelowdrawn = '
     x         ,xmintimelowdrawn,' in iofficial ',istate
     x         ,' and jneighlow = ',jneighlow,' of ',nneigh
     x         ,' ( prefactor: ',xminPREFlowdrawn
     x         ,' -- rate: ',xminRATElowdrawn
     x         ,' -- barrier: ',xminBARRlowdrawn,' ) - B '
c           if(ipr.gt.0)then
c               if(rank_f.eq.0) 
c     x         write(6,'(A,I9,A,I9,A,I9,A,ES16.6,A,ES16.6,A)')
c     x          'SpawnID ',spawnID,' -- jneighlow = ',jneighlow
c     x         ,' of ',nneigh
c     x         ,' ( prefactor: ',statedata(istate)%prefac(jneighlow)
c     x         ,' - barrier: ',statedata(istate)%barrierev(jneighlow),' )'
c               call MPI_BARRIER(force_comm,ier)
c           endif
           call MPI_BARRIER(force_comm,ier)
         elseif(n_states_skipped.gt.0)then
           !! We likely recognize the state but cant put any neighbors into synthetic mode yet...
           if(rank_f.eq.0) write(6,'(A,I9,A,I9,A,I9,A,I9,A)')
     x         'WARNING... SpawnID ',spawnID,' in iofficial ',istate
     x         ,' and jneighlow = ',jneighlow,' of ',nneigh,' - B '
           call MPI_BARRIER(force_comm,ier)
         endif

         if(istate.ne.istate_initial)then
           emin = statedata(istate)%emin
           nnegmin = statedata(istate)%nnegmin
           nprod = statedata(istate)%nprod
           freqlogsmin = statedata(istate)%freqlogsmin
           do j = 1, natom*3
              xyzmin(j) = statedata(istate)%xyzmin(j)
           enddo
           do j = 1, nneigh
              eneigh(j) = statedata(istate)%eneigh(j)
              barrierev(j) = statedata(istate)%barrierev(j)
              barrevev(j) = statedata(istate)%barrevev(j)
              prefac(j) = statedata(istate)%prefac(j)
              nattempted(j) = statedata(istate)%nattempted(j)
           enddo
           do j = 1, (nneigh+1)*natom*3
              xyzneighs(j) = statedata(istate)%xyzneighs(j)
           enddo
         endif

         if(ineighdep.eq.1) then
           ntick(ineighdep) = 1
           ticknext(1,ineighdep) = 1.d30
           prefac(ineighdep) = 1.d-30
           barrierev(ineighdep) = 1.d0
           barrevev(ineighdep) = 1.d0
           if(nneigh.eq.0)then
             nneigh = 1
             do j = 1, natom*3
               if(j.le.nmove)then
                 xyzneighs(ineighdep*natom*3+j) = 0.d0
               else
                 xyzneighs(ineighdep*natom*3+j) = xyzmin(j)
               endif
             enddo
           endif
         endif

        call MPI_BARRIER(force_comm,ier)
        if(n_states_skipped.gt.0)then
          if(rank_f.eq.0) write(*,'(A,I9,A,I10,A,ES16.6,A,ES16.6)')
     +             'SpawnID ',spawnID,
     +             ' just skipped over '
     +             ,n_states_skipped,' states and '
     +             ,tadtime-tadtime_initial,' s to tadtime = '
     +             ,tadtime
          if(rank_f.eq.0) call flushtad(6)
          call MPI_BARRIER(force_comm,ier)
        endif

         return
        end subroutine lp_kMC

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine lp_update_official(nofficialmax)

        implicit real*8(a-h,o-z)  
          include 'parameters.h'
        dimension eofficial(nofficialmax)
        character*80 st1,st2,st3,st4,st5,comstr

        integer msgstatus(MPI_STATUS_SIZE)
        logical msgflag
        real*8, dimension(:), allocatable :: statebuf

        if(execute_serial) return

        msgsource = fmanagerCR
        msgflag = .true.
        iupdated = 0

        do while(msgflag)
        !do while(msgflag.and.(iupdated.lt.lp_avail)) ! Only allow 10 updated before checking if its time to run...

           ! Check For PUT Message:
           if(rank_f.eq.0) call MPI_IProbe(msgsource,pstate_tag2
     +          ,MPI_COMM_WORLD,msgflag,msgstatus,ier)
            call MPI_BCAST(msgflag,1,MPI_LOGICAL,0,force_comm,ier)

           if(msgflag)then

               iupdated = iupdated + 1

               if(rank_f.eq.0) call MPI_Get_count
     +             (msgstatus,MPI_REAL8,ibuffsize,ier)
               call MPI_BCAST(ibuffsize,1,MPI_INTEGER,0,force_comm,ier)
               allocate(statebuf(ibuffsize))
               if(rank_f.eq.0) write(*,*) 'groupID ',groupID
     +           ,' getting state message with ibuffsize ',ibuffsize
               call flush(6)
               if(rank_f.eq.0) call MPI_recv(statebuf,ibuffsize,
     +           MPI_REAL8,msgsource,pstate_tag2,MPI_COMM_WORLD,
     +           msgstatus,ier)
               call MPI_BCAST(statebuf,ibuffsize,MPI_REAL8,0,
     +           force_comm,ier)

               istate = NINT(statebuf(1))

               if(istate.gt.nstatedata) then
                 nstatedata = istate
                 ! If new state... set synthetic stuff to zero
                 statedata(istate)%timehighnow = 0.d0
                 statedata(istate)%timelownow  = 0.d0
               endif

               statedata(istate)%istate = istate
               statedata(istate)%nmove = NINT(statebuf(2))
               nmove = statedata(istate)%nmove
               natom = statedata(istate)%nmove + natomfixed
               statedata(istate)%emin = statebuf(3)
               eofficial(istate) = statebuf(3)
               nneigh = NINT(statebuf(4))

               call check_statedata(istate,natom,nneigh)

               if(ibuffsize.ne.
     +            (7+natom*3+5*nneigh+(nneigh+1)*natom*3+1))then
                  if(rank_f.eq.0) write(*,*) 'groupID ',groupID
     +             ,' ERROR in lp_update !!'
     +             ,' natom = ',natom
     +             ,' nneigh = ',nneigh
     +             ,' ibuffsize = ',ibuffsize
     +             ,' expected-ibuffsize = '
     +             ,(7+natom*3+5*nneigh+(nneigh+1)*natom*3+1)
                  if(rank_f.eq.0) call flushtad(6)
               endif

                if(nneigh.lt.1)then
                  if(rank_f.eq.0) write(*,*) 'groupID ',groupID
     +             ,' ERROR in lp_update !!  nneigh = ',nneigh
     +             ,' for istate = ',istate
                  if(rank_f.eq.0) call flushtad(6)
                endif

               statedata(istate)%nneigh = nneigh
               statedata(istate)%nnegmin = NINT(statebuf(5))
               statedata(istate)%nprod = NINT(statebuf(6))
               statedata(istate)%freqlogsmin = statebuf(7)

               do j = 1, natom*3
                   statedata(istate)%xyzmin(j) =
     +                    statebuf(7+j)
               enddo
               xminbar = 1.d99
               do j = 1, nneigh
                   statedata(istate)%eneigh(j) =
     +                    statebuf(7+j+natom*3+0*nneigh)
                   statedata(istate)%barrierev(j) =
     +                    statebuf(7+j+natom*3+1*nneigh)
                   statedata(istate)%barrevev(j) =
     +                    statebuf(7+j+natom*3+2*nneigh)
                   statedata(istate)%prefac(j) =
     +                    statebuf(7+j+natom*3+3*nneigh)
                   statedata(istate)%nattempted(j) =
     +                    statebuf(7+j+natom*3+4*nneigh)
                   if(statedata(istate)%barrevev(j).lt.xminbar)then
                     xminbar = statedata(istate)%barrevev(j)
                   endif
               enddo
               statedata(istate)%barrierminev = xminbar
               j_start = (7+natom*3+5*nneigh)
               do j = 1, (nneigh+1)*natom*3
                 statedata(istate)%xyzneighs(j) = statebuf(j_start+j)
               enddo

               j_start = j_start + (nneigh+1)*natom*3
               timehighprev = statebuf(j_start+1)
               if(lsynthattempts)then
                 statedata(istate)%timehighprev = timehighprev
               endif

               deallocate(statebuf)
               pstate_tag2 = pstate_tag2 + 1

            endif

         enddo

         return
        end subroutine lp_update_official

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine maxdr_rz(nmove,x1,x2,drmax)
        implicit real*8(a-h,o-z)
        dimension x1(3,nmove),x2(3,nmove)

         drmax=0.0d0
         do j=1,nmove
           dr=0.0d0
           do k=1,3
             dr=dr+(x1(k,j)-x2(k,j))**2
           enddo
           !if (drmax.lt.dr) jmax=j
           drmax=max(drmax,dr)
         enddo

         drmax=sqrt(drmax)
         return
        end subroutine maxdr_rz

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine maxdr_rz_2(nmove,x1,x2,drmaxij,drmaxik,
     +                              drmaxiji,drmaxiki,imax,jmax,kmax)
        implicit real*8(a-h,o-z)
        dimension x1(3,nmove),x2(3,nmove)

         ! First determine which atom moved the most:
         drmax=0.0
         imax = 0
         do i = 1, nmove
             dr = 0.0d0
             do m=1,3
                 dr=dr+(x2(m,i)-x1(m,i))**2
             enddo
             if (dr.gt.drmax) then
                 drmax = dr
                 imax = i
             endif
         enddo

         ! Now find the 1st most-stretch bond including atom imax:
         drmaxij =0.0
         drmaxiji=0.0
         dratmax =0.0
         jmax = 0
         do j = 1, nmove
             if (imax.ne.j) then
                 dr1 = 0.0d0
                 dr2 = 0.0d0
                 do m=1,3
                     dr1=dr1+(x1(m,imax)-x1(m,j))**2
                     dr2=dr2+(x2(m,imax)-x2(m,j))**2
                 enddo
                 if ((dr2/dr1).gt.dratmax) then
                     drmaxij  = dr2
                     drmaxiji = dr1
                     dratmax  = dr2/dr1
                     jmax = j
                 endif
             endif
         enddo
         drmaxij=sqrt(drmaxij)
         drmaxiji=sqrt(drmaxiji)

         ! Now find the 2nd most-stretch bond including atom imax:
         drmaxik =0.0
         drmaxiki=0.0
         dratmax =0.0
         kmax = 0
         do k = 1, nmove
             if ((imax.ne.k).and.(jmax.ne.k))then
                 dr1 = 0.0d0
                 dr2 = 0.0d0
                 do m=1,3
                     dr1=dr1+(x1(m,imax)-x1(m,k))**2
                     dr2=dr2+(x2(m,imax)-x2(m,k))**2
                 enddo
                 if ((dr2/dr1).gt.dratmax) then
                     drmaxik  = dr2
                     drmaxiki = dr1
                     dratmax  = dr2/dr1
                     kmax = k
                 endif
             endif
         enddo
         drmaxik=sqrt(drmaxik)
         drmaxiki=sqrt(drmaxiki)

         return
        end subroutine maxdr_rz_2

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine pr_check_deactivate(transent)
          implicit real*8(a-h,o-z)
          logical msgflag, inbuffer(2), transent
          integer msgstatus(MPI_STATUS_SIZE)
          !activated = .true.
          !lfinished = .false.
!          if(rank_f.eq.0) write(6,*) 'spawnID ',spawnID
!     +             ,' rank_l ',rank_l,' probing deac_tag1..'
          call MPI_IProbe(0,deac_tag1,local_comm,
     +                                msgflag,msgstatus,ier)
          if(msgflag.and.(rank_f.gt.0))then
            write(*,*) 'spawnID,rank_l ',spawnID,rank_l
     +             ,' getting deac_tag1 message'
            call flush(6)
            STOP 'wrong rank recieving deac_tag1 message'
          endif
          call MPI_BCAST(msgflag,1,MPI_LOGICAL,0,force_comm,ier)
          if(msgflag)then 
            if(rank_f.eq.0) call MPI_recv(inbuffer,2,MPI_LOGICAL,
     +                     0,deac_tag1,local_comm,msgstatus,ier)
            call MPI_BCAST(inbuffer,2,MPI_LOGICAL,0,force_comm,ier)
            activated = inbuffer(1)
            lfinished = inbuffer(2)
!            if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
!     +       ,' rank_l ',rank_l
!     +       ,' - deactivate message: activated = ',activated
!     +       ,' - transent = ',transent
!            if(rank_f.eq.0) call flush(6)
!            if(rank_f.eq.0) call MPI_Send(.true.,1,MPI_LOGICAL
!     +                                 ,0,deac_tag2,local_comm,ier)
          endif
          call MPI_BARRIER(force_comm,ier)
         return
        end subroutine pr_check_deactivate

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine pr_send_deactivate()
          implicit real*8(a-h,o-z)
          logical outbuffer(2)
          
          outbuffer(1) = activated
          outbuffer(2) = lfinished
          if(nParRep.gt.1)then
          ireplast = nParRep-1
          do i = 1, ireplast !nprocs_l-1
            msgdest = i*nforcecores
            if(rank_f.eq.0) call MPI_ISend(outbuffer,2,MPI_LOGICAL,
     +           msgdest,deac_tag1,local_comm,deac_req(msgdest),ier)
            call MPI_BARRIER(force_comm,ier)
          enddo
          endif
         return
        end subroutine pr_send_deactivate

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine pr_get_coords(natom,xyz,lprcoords)
          implicit real*8(a-h,o-z)
          logical msgflag,lprcoords
          real*8, dimension(:), allocatable :: inbuffer
          integer msgstatus(MPI_STATUS_SIZE),imessage_size
          dimension xyz(mxatom2*3)
          call MPI_IProbe(0,prcoords_tag1,local_comm,
     +                                          msgflag,msgstatus,ier)
          if(msgflag.and.(rank_f.gt.0))then
            write(*,*) 'spawnID,rank_l ',spawnID,rank_l
     +             ,' getting prcoords_tag1 message'
            STOP 'wrong rank recieving prcoords_tag1 message'
          endif
          call MPI_BCAST(msgflag,1,MPI_LOGICAL,0,force_comm,ier)
          if(msgflag)then
            if(rank_f.eq.0)then
              call MPI_Get_count(msgstatus,MPI_REAL8,imessage_size,ier)
              allocate(inbuffer(imessage_size))
              call MPI_recv(inbuffer,imessage_size,MPI_REAL8,
     +           0,prcoords_tag1,local_comm,msgstatus,ier)
              natom = NINT(inbuffer(1))
              do i = 1,natom*3
                xyz(i) = inbuffer(1+i)
              enddo
              lprcoords = .true.
              deallocate(inbuffer)
            endif
!            if(rank_f.eq.0) write(*,*) 'spawnid,rank_l '
!     x                               ,spawnid,rank_l,' - has coords'
            call MPI_BCAST(natom,1,MPI_INTEGER,0,force_comm,ier)
            call MPI_BCAST(xyz,natom*3,MPI_REAL8,0,force_comm,ier)
            call MPI_BCAST(lprcoords,1,MPI_LOGICAL,0,force_comm,ier)
          endif
         return
        end subroutine pr_get_coords

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine pr_send_coords(natom,xyz)
          implicit real*8(a-h,o-z)
          real*8, dimension(:), allocatable :: outbuffer
          integer imessage_size
          dimension xyz(mxatom2*3)
          imessage_size = 1+natom*3
          allocate(outbuffer(imessage_size))
          outbuffer(1) = REAL(natom,8)
          do i = 1, natom*3
            outbuffer(i+1) = xyz(i)
          enddo
          if(nParRep.gt.1)then
            ireplast = nParRep-1
            do i = 1, ireplast !nprocs_l-1
              msgdest = i*nforcecores
              if(rank_f.eq.0)
     +          call MPI_Send(outbuffer,imessage_size,MPI_REAL8,
     +            msgdest,prcoords_tag1,local_comm,ier)
              call MPI_BARRIER(force_comm,ier)
            enddo
          endif
          deallocate(outbuffer)
         return
        end subroutine pr_send_coords

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine pr_send_trans(nmove,natom,ihold,xyz,xyzhold
     +             ,nholdmaxworstcase,iblock,nmd,dt,transent,hightime
     +             ,rat_lower,d_rat)
          implicit real*8(a-h,o-z)
          real*8 dt,rat_lower,d_rat
          integer msgstatus(MPI_STATUS_SIZE),imessage_size
          integer iblock,nmd,nmove
          dimension xyz(mxatom2*3)
          dimension xyzhold(3*mxatom2*nholdmaxworstcase)
          logical transent

          ! Find random time near "refined" transition time:
          if(rank_f.eq.0) tratio = prngen(0)*d_rat + rat_lower
          call MPI_BCAST(tratio,1,MPI_REAL8,0,force_comm,ier)

          wctblock = prwctarray(iblock) - prwctarray(iblock-1)
          transwct = prwctarray(iblock-1) + wctblock*tratio
          transmdt = prmdtarray(iblock-1) + nmd*dt*tratio
          prwctarray(iblock) = transwct
          prmdtarray(iblock) = transmdt
          ! MUST CORRECT THE "HIGHTIME":
          hightime = hightime - nmd*dt*(1.0-tratio)
          imessage_size = 1+nmove*3+2+(1+ihold*3*natom)
          if(allocated(prtrans_buf))then
            if(rank_f.eq.0)then
!              write(*,*) 'spawnID ',spawnID,' rank_l ',rank_l
!     +             ,' ERROR -- waiting to send transition info... '
              call MPI_Wait(prtrans_req,msgstatus,ier)
            endif
            call MPI_BARRIER(force_comm,ier)
            deallocate(prtrans_buf)
          endif
          allocate(prtrans_buf(imessage_size))
          prtrans_buf(1) = REAL(spawnID,8)
          do i = 1, nmove*3
            prtrans_buf(i+1) = xyz(i)
          enddo
          prtrans_buf(1+nmove*3+1) = transmdt
          prtrans_buf(1+nmove*3+2) = transwct
          prtrans_buf(1+nmove*3+2+1) = REAL(ihold,8)
          do i = 1, ihold*3*natom
            prtrans_buf(i+1+nmove*3+2+1) = xyzhold(i)
          enddo
!          if(rank_f.eq.0)
!     +        write(*,*) 'spawnID ',spawnID,' rank_l ',rank_l
!     +              ,' transmdt = ',transmdt,' transwct = ',transwct
!     +              ,' sending transition to master - tratio = ',tratio
!     +              ,' prwctarray(iblock-1) = ',prwctarray(iblock-1)
!     +              ,' prwctarray(iblock) = ',prwctarray(iblock)
!     +              ,' wctblock = ',wctblock,' iblock = ',iblock
          if(rank_f.eq.0)
     +        call MPI_ISend(prtrans_buf,imessage_size,MPI_REAL8,0
     +                       ,prtrans_tag1,local_comm,prtrans_req,ier)
          transent = .false. ! Have not confirmed that the send is complete yet


!          write(*,*) 'C rank_l= ',rank_l,' - iblock= ',
!     +      iblock,' - prwct, prmdt: ',
!     +      prwctarray(iblock), prmdtarray(iblock),
!     +      ' ROLLBACK:',(1.0-tratio)


          call MPI_BARRIER(force_comm,ier)
         return
        end subroutine pr_send_trans

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine pr_test_trans(transent)
          implicit real*8(a-h,o-z)
          integer msgstatus(MPI_STATUS_SIZE)
          logical transent
          if(.not.transent)then
            if(rank_f.eq.0)then
!              write(6,*) 'spawnID ',spawnID,' rank_l ',rank_l
!     +             ,' waiting for old transition messages..'
              call MPI_Test(prtrans_req,transent,msgstatus,ier)
            endif
            call MPI_Bcast(transent,1,MPI_LOGICAL,0,force_comm,ier)
            if(transent) deallocate(prtrans_buf)
          endif
         return
        end subroutine pr_test_trans

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine pr_clean_trans()
          implicit real*8(a-h,o-z)
          integer msgstatus1(MPI_STATUS_SIZE),imessage_size,msgsource
          integer msgstatus2(MPI_STATUS_SIZE),icnt
          logical msgflag,msgflag2,needclean,ldummy
          real*8, dimension(:), allocatable :: inbuffer
          integer, dimension(:), allocatable :: itest

          if(rank_l.gt.0) return ! Only need one core here
          iloop = 0
          needclean = .true.
          allocate(itest(nParRep-1))
          itest(:) = 0
          do while(needclean)

            ! Check For Late Transition Messages:
            call MPI_IProbe(MPI_ANY_SOURCE,prtrans_tag1,
     +          local_comm,msgflag,msgstatus1,ier)
            if(msgflag)then ! If there is a message waiting...
              ! Get Information about transition
              msgsource = msgstatus1(MPI_SOURCE) ! local rank sending message
              call MPI_Get_count(msgstatus1,MPI_REAL8,imessage_size,ier)
              allocate(inbuffer(imessage_size))
              call MPI_Recv(inbuffer,imessage_size,MPI_REAL8,
     +               msgsource,prtrans_tag1,local_comm,msgstatus1,ier)
              deallocate(inbuffer)
            endif

            ! Need more loops?
            icnt = 0
            if(nParRep.gt.1)then
            ireplast = nParRep-1
            do i = 1, ireplast !nprocs_l-1
              msgsource = i*nforcecores
              if(itest(i).eq.0)then
                call MPI_IProbe(msgsource,deac_tag2,local_comm,
     +                                    msgflag2,msgstatus2,ier)
                if(msgflag2)then
                  call MPI_Recv(ldummy,1,MPI_LOGICAL,
     +               msgsource,deac_tag2,local_comm,msgstatus2,ier)
                  itest(i) = 1
                  icnt = icnt + 1
!                  write(*,*) 'spawnID ',spawnID
!     +               ,' got deact confirmation from rank_l ',msgsource
!                  call flush(6)
                endif
              else
                icnt = icnt + 1
              endif
            enddo
            endif
            iloop = iloop + 1
            if(icnt.eq.(nParRep-1))then
              needclean = .false.
            elseif(MOD(iloop,10000000).eq.0)then
              write(*,*) 'spawnID ',spawnID
     +           ,' at pr_clean_trans(), iteration ',iloop
              call flush(6)
            endif

          enddo
          deallocate(itest)
         return
        end subroutine pr_clean_trans

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine pr_check_trans(nmove,natom,ihold,xyz,xyzhold
     +          ,nholdmaxworstcase,totaltime,lprtrans
     +          ,iprsafe,transwct,wctimehighstop,lrepstold)
          implicit real*8(a-h,o-z)
          integer msgstatus(MPI_STATUS_SIZE),imessage_size,msgsource
          integer msgstatus2(MPI_STATUS_SIZE),msgsource2
          integer transync_req(nParRep)
          logical msgflag,msgflag2,lprtrans,iprsafe,lrepstold
          real*8, dimension(:), allocatable :: inbuffer
          real*8  totaltime,nexttime,transwct,wctimehighstop
          real*8  recvbuf(2),transwctbuf(2)
          dimension xyz(mxatom2*3)
          dimension xyzhold(3*mxatom2*nholdmaxworstcase)
          ! Check For Spawn Message:
          if(rank_f.eq.0)
     +      call MPI_IProbe(MPI_ANY_SOURCE,prtrans_tag1,local_comm,
     +          msgflag,msgstatus,ier)
          call MPI_BCAST(msgflag,1,MPI_LOGICAL,0,force_comm,ier)
          if(msgflag)then ! If there is a message waiting...

             lprtrans = .true.
             
             if(rank_f.eq.0)then
!               write(*,*) 'spawnID ',spawnID
!     +              ,' master getting replica transition info '
               ! Get Information about transition
               msgsource = msgstatus(MPI_SOURCE) ! local rank sending message
               call MPI_Get_count(msgstatus,MPI_REAL8,imessage_size,ier)
               allocate(inbuffer(imessage_size))
               call MPI_Recv(inbuffer,imessage_size,MPI_REAL8,
     +               msgsource,prtrans_tag1,local_comm,msgstatus,ier)
               ispawnid = NINT(inbuffer(1))
             endif
             call MPI_BCAST(ispawnid,1,MPI_INTEGER,0,force_comm,ier)
             call MPI_BCAST(msgsource,1,MPI_INTEGER,0,force_comm,ier)
             if(spawnID.ne.ispawnid)then
               if(rank_f.eq.0) deallocate(inbuffer)
               write(*,*) 'spawnID ',spawnID
     +              ,' not matching ispawnid ',ispawnid
               return ! Lagging message is not for this state
             endif
             if(rank_f.eq.0)then
               do i = 1,nmove*3
                 xyz(i) = inbuffer(1+i)
               enddo
               transmdt = inbuffer(1+nmove*3+1)
               transwct = inbuffer(1+nmove*3+2)
               ihold = NINT(inbuffer(1+nmove*3+2+1))
               do i = 1,natom*3*ihold
                 xyzhold(i) = inbuffer(1+nmove*3+2+1+i)
               enddo
               deallocate(inbuffer)
             endif
             call MPI_BCAST(xyz,nmove*3,MPI_REAL8,0,force_comm,ier)
             call MPI_BCAST(transmdt,1,MPI_REAL8,0,force_comm,ier)
             call MPI_BCAST(transwct,1,MPI_REAL8,0,force_comm,ier)
             call MPI_BCAST(ihold,1,MPI_INTEGER,0,force_comm,ier)
             call MPI_BCAST
     +            (xyzhold,natom*3*ihold,MPI_REAL8,0,force_comm,ier)

             ! Don't bother with this transition if it can't win!
             ! ... at first -- send messages to tell replicas to stop MD
             if(lrepstold.and.(transwct.ge.wctimehighstop))then
               msgflag = .false.
               lprtrans = .false.
               if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
     +              ,' is ignoring a transition at wct = ',transwct
               if(rank_f.eq.0) call flush(6)
               call MPI_BARRIER(force_comm,ier)
               return
             endif

             ! Now get associated times from other ParRep Cores:
             totaltime = 0.d0
             if(nParRep.gt.1)then

             ireplast = nParRep-1
             transwctbuf(1) = transwct
             transwctbuf(2) = wctimehighstop

             ! Bcast transition time to replicas:
             do i = 1, ireplast !nprocs_l-1
               msgdest = i*nforcecores
               if(msgdest.ne.msgsource)then
                 if(rank_f.eq.0)then
                   call MPI_ISend(transwctbuf,2,MPI_REAL8,msgdest
     +                  ,prtrans_tag2,local_comm,transync_req(i),ier)
                 endif
                 call MPI_BARRIER(force_comm,ier)
               endif
             enddo

             ! Read in times for all replicas:
             idone = 0
             do while (idone.lt.(ireplast-1))
               nexttime = 0.d0
               if(rank_f.eq.0)
     +           call MPI_IProbe(MPI_ANY_SOURCE,prtrans_tag3,
     +               local_comm,msgflag2,msgstatus2,ier)
               call MPI_BCAST(msgflag2,1,MPI_LOGICAL,0,force_comm,ier)
               if(msgflag2)then ! If there is a message waiting...
                 if(rank_f.eq.0)then
                   msgsource2 = msgstatus2(MPI_SOURCE)
!                   write(*,*) 'spawnID ',spawnID
!     +                ,' master syncronizing trans time with rank_l '
!     +                ,msgsource2
                   call MPI_Recv(recvbuf,2,MPI_REAL8,msgsource2
     +                         ,prtrans_tag3,local_comm,msgstatus2,ier)
                   nexttime = recvbuf(1)
                   if(NINT(recvbuf(2)).ne.1) iprsafe = .false.
                   totaltime = totaltime + nexttime
                 endif
                 idone = idone + 1
                 call MPI_BARRIER(force_comm,ier)
               endif
             enddo
             endif
             call MPI_BCAST(totaltime,1,MPI_REAL8,0,force_comm,ier)
             call MPI_BCAST(iprsafe,1,MPI_LOGICAL,0,force_comm,ier)
             totaltime = totaltime + transmdt
!             write(*,*) 'spawnID ',spawnID
!     +                ,' master DONE syncronizing trans time'

             ! Don't bother with this transition if it can't win!
             if(transwct.ge.wctimehighstop)then
               totaltime = 0.d0
               lprtrans = .false.
               lrepstold = .true.
               if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
     +           ,' is FIRST ignoring a transition at wct = ',transwct
               if(rank_f.eq.0) call flush(6)
               call MPI_BARRIER(force_comm,ier)
             endif

          endif
         return
        end subroutine pr_check_trans

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine pr_sync_trans(iblock,dt,transent,lprneeded)
          implicit real*8(a-h,o-z)
          integer msgstatus(MPI_STATUS_SIZE),iblock
          logical msgflag
          real*8 sendbuf(2),transwctbuf(2),wctimehighstop
          real*8 lowbound,lowboundmdt
          logical transent,lprneeded

          ! Check For Spawn Message:
          if(rank_f.eq.0)then
!            write(*,*) 'spawnID,rank_l ',spawnID,rank_l
!     +             ,' probing in pr_sync_trans()'
!            call flush(6)
            call MPI_IProbe(0,prtrans_tag2,local_comm,
     +          msgflag,msgstatus,ier)
          endif
          call MPI_BCAST(msgflag,1,MPI_LOGICAL,0,force_comm,ier)
          if(msgflag)then ! If there is a message waiting...
             call pr_test_trans(transent)
             if(rank_f.eq.0)then
!               write(*,*) 'spawnID,rank_l ',spawnID,rank_l
!     +             ,' recieving pr_sync_trans() message'
!               call flush(6)
             endif
             if(rank_f.eq.0) call MPI_Recv(transwctbuf,2,MPI_REAL8,
     +               0,prtrans_tag2,local_comm,msgstatus,ier)
             call MPI_BCAST(transwctbuf,2,MPI_REAL8,0,force_comm,ier)
             transwct = transwctbuf(1)
             wctimehighstop = transwctbuf(2)
             if(lprneeded.and.wctimehighstop.le.(1.d30))then
               call MPI_BARRIER(force_comm,ier)
               lprneeded = .false.
!               if(rank_f.eq.0) write(*,*) 'spawnID,rank_l ',spawnID
!     +           ,rank_l,' is setting lprneeded to ',lprneeded
!               if(rank_f.eq.0) call flush(6)
             endif
             ifound = 0
             tratio = 0.d0
             dtime = 0.d0
             lowbound = 0.d0
             uppbound = 0.d0
             if(iblock.le.1)then
              timesend = 0.d0
             else
               do i = 1,iblock-1
                 uppbound = prwctarray(i+1)
                 lowbound = prwctarray(i)
                 if((transwct.gt.lowbound)
     +                          .and.(transwct.le.uppbound))then
                   ifound = i
                   Exit
                 endif
               enddo
               timesend = 0.d0
               if(ifound.eq.0)then
                 if(transwct.gt.prwctarray(iblock))then
                   timesend = prmdtarray(iblock)
                 else
                   timesend = 0.d0
                 endif
               else
                 uppbound = prwctarray(ifound+1)
                 lowbound = prwctarray(ifound)
                 uppboundmdt = prmdtarray(ifound+1)
                 lowboundmdt = prmdtarray(ifound)
                 if(uppboundmdt.ge.(lowboundmdt+dt))then
                   tratio = (transwct-lowbound)/(uppbound-lowbound)
                   dtime = uppboundmdt - lowboundmdt
                   timesend = lowboundmdt + tratio*dtime
                 else
                   timesend = prmdtarray(ifound)
                 endif
               endif
             endif
!             write(*,*) 'spawnID ',spawnID,' rank_l ',rank_l
!     +                  ,' - timesend = ',timesend
!     +                  ,' - tratio = ',tratio
!     +                  ,' - dtime = ',dtime
!     +                  ,' - transwct = ',transwct
!     +                  ,' - lowbound = ',lowbound
!     +                  ,' - uppbound = ',uppbound
             sendbuf(1) = timesend
             sendbuf(2) = 1.0
             if(.not.transent)then
               sendbuf(2) = 0.0
               ! If our waiting transition occured after the one we are
               ! syncing with -- don't need to warn the master:
               if(prwctarray(iblock).ge.transwct) sendbuf(2) = 1.0
               ! If our waiting transition occured after one that was beyond the
               ! stop time -- don't need to warn the master:
               if(prwctarray(iblock).ge.wctimehighstop) sendbuf(2) = 1.0
             endif
             if(rank_f.eq.0) call MPI_Send(sendbuf,2,MPI_REAL8,0
     +                                    ,prtrans_tag3,local_comm,ier)
             call MPI_BARRIER(force_comm,ier)
          endif
         return
        end subroutine pr_sync_trans

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine pr_check_time(totaltime,iprsafe,wctimehighstop)
          implicit real*8(a-h,o-z)
          integer msgstatus(MPI_STATUS_SIZE)
          logical msgflag,iprsafe
          real*8  totaltime,nexttime,wctimehighstop
          real*8  recvbuf(2)
          totaltime = 0.d0

!          call MPI_BCAST(wctimehighstop,1,MPI_REAL8,0,local_comm,ier)
!          call MPI_Reduce()
! int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
!               MPI_Op op, int root, MPI_Comm comm)

          if(rank_f.eq.0)then
!            write(*,*) 'spawnID ',spawnID,' gathering ParRep time.'
!            call flushtad(6)
            if(nParRep.gt.1)then
            ireplast = nParRep-1
            do i = 1, ireplast !nprocs_l-1
              nexttime = 0.d0
              msgdest = i*nforcecores
!              write(*,*) 'spawnID,rank_l ',spawnID,rank_l
!     +              ,' getting time from rank_l = ',msgdest
!              call flushtad(6)
              call MPI_Send(wctimehighstop,1,MPI_REAL8,msgdest
     +                                    ,prtime_tag1,local_comm,ier)
              call MPI_Recv(recvbuf,2,MPI_REAL8,msgdest
     +                          ,prtime_tag2,local_comm,msgstatus,ier)
              nexttime = recvbuf(1)
              if(NINT(recvbuf(2)).ne.1) iprsafe = .false.
              totaltime = totaltime + nexttime
              if(.not.iprsafe) EXIT !! Assumes time will NOT be used anyway
                                    !! Should also do nonblocking in future
            enddo
            endif
!            write(*,*) 'spawnID ',spawnID,' DONE gathering ParRep time.'
!            call flushtad(6)
          endif
          call MPI_BCAST(totaltime,1,MPI_REAL8,0,force_comm,ier)
          call MPI_BCAST(iprsafe,1,MPI_LOGICAL,0,force_comm,ier)
         return
        end subroutine pr_check_time

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine pr_share_time(iblock,transent,lprneeded)
          implicit real*8(a-h,o-z)
          integer msgstatus(MPI_STATUS_SIZE),msgsource
          logical msgflag,ldummy
          real*8  nexttime,wctimehighstop
          real*8  sendbuf(2)
          logical transent,lprneeded
          ! Check For Spawn Message:
!          if(rank_f.eq.0) write(*,*) 'spawnID,rank_l ',spawnID,rank_l
!     +             ,' probing in pr_share_time()'
          call MPI_IProbe(MPI_ANY_SOURCE
     +                   ,prtime_tag1,local_comm,msgflag,msgstatus,ier)
          if(msgflag.and.(rank_f.gt.0))then
                 write(*,*) 'spawnID,rank_l ',spawnID,rank_l
     +             ,' getting prtime_tag1 message'
                 STOP 'wrong rank recieving prtime_tag1 message'
          endif
          call MPI_BCAST(msgflag,1,MPI_LOGICAL,0,force_comm,ier)
          if(msgflag)then ! If there is a message waiting...
             call pr_test_trans(transent)
!            if(rank_f.eq.0) write(*,*) 'spawnID,rank_l ',spawnID,rank_l
!     +             ,' sharing time.'
!            if(rank_f.eq.0) call flushtad(6)
            if(rank_f.eq.0)then
               msgsource = msgstatus(MPI_SOURCE)
               if(msgsource.ne.0)then
                 write(*,*) 'spawnID,rank_l ',spawnID,rank_l
     +             ,' getting message from rank_l = ',msgsource
                 STOP 'getting prtime_tag1 mesasge from bad source'
               endif
            endif
            if(rank_f.eq.0) call MPI_Recv(wctimehighstop,1,MPI_REAL8,0
     +                          ,prtime_tag1,local_comm,msgstatus,ier)
            call MPI_BCAST(wctimehighstop,1,MPI_REAL8,0,force_comm,ier)
             if(lprneeded.and.wctimehighstop.le.(1.d30))then
               call MPI_BARRIER(force_comm,ier)
               lprneeded = .false.
!               if(rank_f.eq.0)
!     +           write(*,*) 'spawnID,rank_l ',spawnID,rank_l
!     +           ,' is setting lprneeded to ',lprneeded,' in time-share'
!               if(rank_f.eq.0) call flush(6)
             endif
            currenttime = prmdtarray(iblock)
            sendbuf(1) = currenttime
            sendbuf(2) = 1.0
            if(.not.transent)then
              sendbuf(2) = 0.0
              ! If our waiting transition occured after one that was at/beyond the
              ! stop time -- don't need to warn the master:
              if(prwctarray(iblock).ge.wctimehighstop) sendbuf(2) = 1.0
            endif
            if(rank_f.eq.0) call MPI_Send(sendbuf,2,MPI_REAL8,0
     +                                    ,prtime_tag2,local_comm,ier)
!            if(rank_f.eq.0) write(*,*) 'spawnID,rank_l ',spawnID,rank_l
!     +             ,' DONE sharing time.'
!            if(rank_f.eq.0) call flushtad(6)
          endif
          call MPI_BARRIER(force_comm,ier)
         return
        end subroutine pr_share_time

        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        !
        !
        !
        !=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!!=!=!=!=!
        subroutine check_statedata(istate,natom,nneigh)
          implicit real*8(a-h,o-z)
          integer istate,natom,nneigh,ninc,ninci,oldsize,newsize
          integer oldlen,newlen,ninc_initial

          ! Use these tmp arrays to reajust allocations:
          integer, dimension(:), allocatable :: intbuffer
          real*8, dimension(:), allocatable :: realbuffer

          ! Use this worker routine to make sure the statedata arrays
          ! are allocated correctly for this state

          ninc = 25
          ninc_initial = 1   ! First allocation can be small
                             ! (In case the state is never visited)

          !! Allocate arrays if they havent been before:
          if(.not.allocated(statedata(istate)%nattempted))then
          
              ninci = ninc_initial
              if(nneigh.gt.ninc_initial)then
                  ninci = (INT(nneigh/ninc_initial)+1)*ninc_initial
              endif

!              write(*,*)
!     +          'GroupID ',GroupID,' allocating statedata for istate '
!     +          ,istate,' with ninci = ',ninci

              allocate(statedata(istate)%nattempted(ninci))
              statedata(istate)%nattempted(:) = 0
              allocate(statedata(istate)%xyzneighs((ninci+1)*natom*3))
              statedata(istate)%xyzneighs(:) = 0.d0
              allocate(statedata(istate)%eneigh(ninci))
              statedata(istate)%eneigh(:) = 0.d0
              allocate(statedata(istate)%barrierev(ninci))
              statedata(istate)%barrierev(:) = 1.d0
              allocate(statedata(istate)%barrevev(ninci))
              statedata(istate)%barrevev(:) = 0.d0
              allocate(statedata(istate)%prefac(ninci))
              statedata(istate)%prefac(:) = 0.d0
              allocate(statedata(istate)%neighlabel(ninci))
              statedata(istate)%neighlabel(:) = 0

              ! For now, only "block" first 10 neighbors:
              allocate(statedata(istate)%ibondneigh1(10))
              statedata(istate)%ibondneigh1(:) = 0
              allocate(statedata(istate)%ibondneigh2(10))
              statedata(istate)%ibondneigh2(:) = 0
              allocate(statedata(istate)%ibondneigh3(10))
              statedata(istate)%ibondneigh3(:) = 0
              allocate(statedata(istate)%dbondneigh12(10))
              statedata(istate)%dbondneigh12(:) = 0.d0
              allocate(statedata(istate)%dbondneigh13(10))
              statedata(istate)%dbondneigh13(:) = 0.d0
              allocate(statedata(istate)%dbondneigh12i(10))
              statedata(istate)%dbondneigh12i(:) = 0.d0
              allocate(statedata(istate)%dbondneigh13i(10))
              statedata(istate)%dbondneigh13i(:) = 0.d0

          !! Otherwise, make arrays larger if there are too many neighbors:
          elseif(nneigh.gt.INT(SIZE(statedata(istate)%nattempted)))then

             oldsize = INT(SIZE(statedata(istate)%nattempted))
             oldlen = INT(SIZE(statedata(istate)%xyzneighs))

             newsize = (INT(nneigh/ninc)+1)*ninc
             newlen = (newsize+1)*natom*3

             write(*,*)
     +        'GroupID ',GroupID,' updating statedata to nneigh_max = ',
     +         newsize,' from ',oldsize,' for istate ',istate
             call flushtad(6)

             ! Allocate tmp arrays to max size needed:
             allocate(intbuffer( oldsize ))
             allocate(realbuffer( oldlen ))

             ! Re-allocate %nattempted
             intbuffer(1:oldsize) 
     +           = statedata(istate)%nattempted(1:oldsize)
             deallocate(statedata(istate)%nattempted)
             allocate(statedata(istate)%nattempted(newsize))
             statedata(istate)%nattempted(:) = 0
             statedata(istate)%nattempted(1:oldsize) 
     +           = intbuffer(1:oldsize)

             ! Re-allocate %xyzneighs (This is a bigger array)
             realbuffer(1:oldlen)=statedata(istate)%xyzneighs(1:oldlen)
             deallocate(statedata(istate)%xyzneighs)
             allocate(statedata(istate)%xyzneighs(newlen))
             statedata(istate)%xyzneighs(:) = 0.d0
             statedata(istate)%xyzneighs(1:min(oldlen,newlen))
     +           = realbuffer(1:min(oldlen,newlen))

             ! Re-allocate %eneigh
             realbuffer(1:oldsize) = statedata(istate)%eneigh(1:oldsize)
             deallocate(statedata(istate)%eneigh)
             allocate(statedata(istate)%eneigh(newsize))
             statedata(istate)%eneigh(:) = 0.d0
             statedata(istate)%eneigh(1:oldsize) = realbuffer(1:oldsize)

             ! Re-allocate %barrierev
             realbuffer(1:oldsize) 
     +           = statedata(istate)%barrierev(1:oldsize)
             deallocate(statedata(istate)%barrierev)
             allocate(statedata(istate)%barrierev(newsize))
             statedata(istate)%barrierev(:) = 1.d0
             statedata(istate)%barrierev(1:oldsize) 
     +           = realbuffer(1:oldsize)

             ! Re-allocate %barrevev
             realbuffer(1:oldsize) 
     +           = statedata(istate)%barrevev(1:oldsize)
             deallocate(statedata(istate)%barrevev)
             allocate(statedata(istate)%barrevev(newsize))
             statedata(istate)%barrevev(:) = 1.d0
             statedata(istate)%barrevev(1:oldsize) 
     +           = realbuffer(1:oldsize)

             ! Re-allocate %prefac
             realbuffer(1:oldsize) = statedata(istate)%prefac(1:oldsize)
             deallocate(statedata(istate)%prefac)
             allocate(statedata(istate)%prefac(newsize))
             statedata(istate)%prefac(:) = 0.d0
             statedata(istate)%prefac(1:oldsize) = realbuffer(1:oldsize)

             ! Re-allocate %neighlabel
             intbuffer(1:oldsize) 
     +           = statedata(istate)%neighlabel(1:oldsize)
             deallocate(statedata(istate)%neighlabel)
             allocate(statedata(istate)%neighlabel(newsize))
             statedata(istate)%neighlabel(:) = 0
             statedata(istate)%neighlabel(1:oldsize) 
     +           = intbuffer(1:oldsize)

             ! Deallocate tmp arrays:
             deallocate(intbuffer)
             deallocate(realbuffer)

          endif

         return
        end subroutine check_statedata

       END MODULE mod_mpi
