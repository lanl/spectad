c
c  SpecTAD - Principal Author(s):
c            Richard J. Zamora (2015-2017)
c            Theoretical Division,
c            Los Alamos National Laboratory, Los Alamos, NM
c
c  SpecTAD is the state-wise parallel version of the temperature accelerated
c  dynamics method. Built upon the TAD2 program from Los Alamos National
c  Laboratory.
c
c  TAD2 - Original Authors: Blas P. Uberuaga and Arthur F. Voter, 2001-2007
c
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

c T A D 2 - tad3.f - implements temperature accelerated dynamics method
c Principle authors:  Art Voter and Blas Uberuaga  December, 2001 - 2007

c======================================================================
c TAD2 CHANGE LOG
c
c tad-3.43  - cwp Oct 20, 2008
c           - added ifun = 7, LBFGS minimizer, can be 100 times faster
c             than steepest descent
c tad-3.42  - bpu Nov 06, 2007
c           - fixed problem in Verlet algorithm (random numbers were
c             scaled incorrectly, leading to incorrect temperatures)
c           - added ability to stop on given state, keyed on state
c             energy, for LJ38 only
c tad-3.41  - bpu Oct 12, 2007
c           - merged fork of Art's 3.39 and Blas' 3.40
c tad-3.40  - bpu Aug 07, 2007
c           - added new tad mode, itad=9 (userdimeretad), for which the
c             user can supply an initial list of active atoms for
c             use with the dimer-TAD mode (for complex materials
c             where identifying such regions would be tedious to code
c             up); list of atoms supplied in dimer.activeatoms file
c           - fixed NEB so zero barrier to leave state is accepted and
c             not ignored
c tad-3.39  - afv Sept 07, 2007
c           - made some changes to the README file
c tad-3.38  - bpu Jul 09, 2007
c           - added energy type 80 (Buckingham+Morse, LJ, Born-Meyer)
c tad-3.37  - bpu Feb 05, 2007
c           - added in extra pairsum check for energy-equiv structures when irecognize=9
c
c tad_3.36  - bpu Oct 31, 2006
c           - added stuff Art had originally put in 3.34 on dist4
c
c tad_3.35  - bpu Aug 16, 2006
c           - moved some global parameters to parameters.h include file
c           - added nebbailcount to neb_tad.f to make sure no infinite
c             loops on intermediate minimum
c           - fixed ninvolved calc
c           - changed rrecognizecrit from .05 to .004  (for LJ38 runs)
c           - fixed bug in e999subs.f with parameters.h file
c
c tad_3.34  - bpu Jan 24, 2006
c           - fixed a minor bug in the nforce report
c           - fixed a bug in NEB dealing with intermediate states
c           - other minor fixes
c tad_3.33  - bpu Sep 22, 2005
c           - fixed timelowshift for synthtad - 09/29/2005
c           - added fix so that times of events are distributed randomly
c             throughout blok, not always halfway through
c tad_3.32  - bpu Jun 28, 2005
c           - added checks to find bad saddles with LJ case (imodecheck code)
c           - added .boost file
c           - added rdimerredocrit for not restarting dimer recount
c           - redefined rrecognizecrit from 1d-2 to 5d-2
c           - fixed prngen -> 1-prngen to avoid log(0)
c tad_3.31s - afv Jun 23, 2005
c           - added hard-coded minimum barrier specific to LJ38 case
c tad_3.19  - Jun 19, 2003
c           - cleaned up the example .input files
c           - reports TAD confidence for events accepted in dimer modes
c           - cleaned up documentation a bit
c tad_3.18  - May 5, 2003
c           - check on NEB to make sure atoms do not pass closer than 1A, fixes if necessary
c           - dimer parameters added to input file
c           - reverse barrier more fully implemented, ignores transitions if
c             reverse barrier less than some value
c           - super-synthetic mode partially implemented, not finished
c           - if templow is different than what is in state files, variables are
c             reset
c           - maktad modified so that potentials are compiled by setting a flag
c tad_3.14  - Apr 18, 2003
c           - improved dimer routine, by GH
c           - vineyard now done on subsets of system, controlled by rdynmatcrit parameter
c tad_3.13  - Apr 16, 2003
c           - installed reordering matchup capability, accessed by irecognize=2
c           - put limit on runningstate file writes (nrunninglimit)
c           - writing out tick marks
c           - added some items to blok file and times file
c tad_3.12  - Apr 14, 2003
c           - reordered blocks of code to streamline the main logical flow
c           - now the acceptance check is performed in only one place
c           - similarly for the the tick-mark update and blok write
c tad_3.11  - Apr 10, 2003
c           - synthtad mode installed (itad=6)
c           - dimeretad mode installed (itad=2)
c           - synthdimeretad mode installed (itad=7)
c tad_3.08b - Oct 11, 2002
c           - P-TAD has now been installed.  Access it using itad=4 (vineyard) or itad=5 (fixed prefac)
c           - installed faster force-call routine (g0faster2) for EAM with multiple types
c           - put in a timing call to show how fast the forces run, based on the first warmup steps
c           - increased the max number of states (limited by .times file arrays) to 10000.
c           - transferred one or two fixes that Blas had put into the 3.09 version of tad3.f
c           - has Blas's updated descent-check and neb routines
c           - fixed up a few other minor things
c tad_3.08  - Oct 11, 2002
c           - got iequivcheck working (probably), although it is not yet an official input parameter
c           - changed the blok print slightly
c tad_3.07  - Oct 10, 2002
c           - newtad mode (itad=1) now appears to be working, saving time for revisited states
c           - reorganized the print levels.  ipr=4 is now the standard novice print level, and ipr=3 is normal level
c           - put in a "blok" print after every high-T MD block, which greps nicely for a summary of the whole run
c tad_3.03  - Sep 25, 2002
c           - improving use of tick marks, correcting errors there
c           - more precision in compact output files
c tad_3.02  - Sep 24, 2002
c           - state recognition is now really working, and you all the output .dat files,
c             except the "runningstate" files,  are numbered by the
c             official state number (iofficial),
c             rather than the istate number.  The mapping information shows up in the .times file.
c           - "state" files are now written out, containing basically all information to start up in a state.
c             Whenever tad3 realizes that it is visiting a state it has seen before, it reads in the appropriate
c             state file, and therefore doesn't have to refind the saddles, etc.
c           - NOTE:   state numbering now starts at 1, not zero.
c           - hopefully the framework is now ready for the more sophisticated tad modes.
c tad_2.02  - Mar 22, 2002
c           - changed to more sophisticated dimension check for nhold.
c             Now the code specifies "nholdmaxworstcase", which applies to
c             when natom=mxatom.  For fewer atoms, nhold can be greater.
c           - added "nmoved" (ninvolved) to both the .barriers file and the .times file
c           - created case "tetfast", which is tetramer with numin and delta set more
c             aggressively to make it run faster
c
c======================================================================

c For information on how to use this code, see the README file,
c the README.input file, and see the paper
c  M.R. Sorensen and A.F. Voter, J. Chem. Phys.  112, 9599 (2000).
c  "Temperature-accelerated dynamics for simulation of infrequent events"
c Also see
c   F. Montalenti and A.F. Voter, JCP 116, 4819 (2002).


c  subroutine packages needed:
c               tad2.f     (this one)
c               mdsubs_tad2.f
c               e0subs_tad2.f
c               descentcheck.f
c               refinetransition.f
ccc               computemoment.f
c               rotalign.f
c               neb_tad.f
c               nmodes.f
c               wrksubs_tad2.f
c               trpsubs.f
c               msubs.f
c               e67subs.f
c               tersubs_tad2.f
c               rsp.f
c               sdiag.f

c  include files needed:
c               wrksubs_tad2.h

c  input files needed at run time
c    <case>.input  - input file containing parameters (supply to standard input)
c    <casename>.start  - geometry file (clsman format)
c  if using etype=67 (AIREBO), also
c    <case>.coord - AIREBO format coordinate file
c    parmd.in - additional AIREBO parameters

c output files (some optional)
c       <casename>.lis   - output listing (redirect the output to make this)
c       <casename>.times  - list of transition times for accepted transitions
c       <casename>.barriers  - list of all saddles found out of each state
c   clsman-format configuration files ("cluster files")
c       <casename>.runningstate.<istate>.dat files  - every state, in sequence, during tad run
c       <casename>.min.<iofficial>.dat  - minimum for each official state
c       <casename>.min.unpbc.<iofficial>.dat  - unpbc'd minimum for each basin
c       <casename>.<iofficial>.attempt.<iattempt>.dat files - trajectory position, partly minimized,  when transition detected
c       <casename>.<iofficial>.sad.<isad>.dat files - saddle point for this pathway
c       <casename>.<iofficial>.end.<isad>.dat files - end state for this pathway


c Proposed new convention for xxx's  (4/4/03):
c    xxx = marks something we are working on, or that needs to be fixed or noted
c    xxxx = new or urgent item
c    xx   = downgraded, low-urgency status (good for keeping a comment around)
c   Feel free to put dates on the comments.
c   If it seems helpful to indicate who put in the xxx tag, then:
c    xxxa (e.g., xxa, xxxa, or  xxxxa) - entered by art
c    xxxb - entered by blas
c    xxxf - entered by francesco

      ! NOTE.. If not using mod_lammps... must use mod_mpi
      use mod_lammps

      implicit real*8(a-h,o-z)
      character*80 casename
      character*80 lmpcmd
      character*180 filnam,filnam1,filnam2
      character*80 filn,filnb
      character*10 chri
      character*2 ctype
      character*6 cdum6
c      parameter(mxatom=10000,maxtrans=10000)
      include 'parameters.h'
      !include 'mpif.h'
      parameter (maxtrans=10000)
c NOTE: if you change mxatom, you also need to change it in other
c locations:
c mdsubs_tad3: movmax,natmax,mxatom,nhold

      logical lstore_min
      logical lstore_minunpbc
      logical lstore_attempt
      logical lstore_attempt_chs
      logical lstore_sad
      logical lstore_end
      logical lstore_dimer
      logical lstore_neb
      logical lexist
      logical lputyes,lendput,lputyessynthkilled
      logical lsynthspawned
      logical lfirstinstate
      logical lstillsynth,lsendtime
      logical skipdeposition
      logical lprloop
      logical lprcoords
      logical lprwarm
      logical lprblack
      logical lprtrans
      logical lprtime
      logical lprneeded
      logical ldprcheckdone
      logical lsenddone
      logical transent,iprsafe,lrepstold
      logical cleanwarmup, cleanblackout
      logical lnewstate,lnewneigh

      !parameter (lenxyzneighs=1500000)
      parameter (lenxyzneighs=(nneighmax+1)*3*mxatom)

      parameter (nholdmaxworstcase=1000,nofficialmax=750,ninfo=9)
      dimension xyzhold(3*mxatom*nholdmaxworstcase)
      dimension xyz(3*mxatom)
      dimension xyzofficial(3*mxatom*nofficialmax)
      dimension eofficial(nofficialmax)
      dimension pxyz(3*mxatom)
      dimension xyzdep(3*mxatom)
      dimension xyz1(3*mxatom)
      dimension xyz2(3*mxatom)
      dimension xyz3(3*mxatom)
      dimension xyz4(3*mxatom)
      dimension xyz5(3*mxatom)
      dimension xyzs(3*mxatom)
      dimension xyzhot(3*mxatom)
      dimension xyzprev(3*mxatom)
      dimension rdimer(3*mxatom)
      dimension kreorder(mxatom)
      dimension taxes(3)
      dimension cdxyz(3*mxatom)
      dimension xyzneighs(lenxyzneighs)
      dimension xyzneighsads(lenxyzneighs)
      dimension rdimermodes(lenxyzneighs)
      dimension itype(mxatom)
      dimension moving(mxatom),mx(mxatom),my(mxatom),mz(mxatom)
      dimension grad(3,mxatom)

      integer ldepspawn
      !real*8 maxsynthrate

      character*20 ratemode
      character*20 tadmode
      character*20 synthrequirement

c  arrays over neighbors (i.e., neighboring states)
c      parameter (nneighmax=1500) ! now in include
      parameter (ntickmax=50)
      dimension esaddle(nneighmax)  ! hartree
      dimension eneigh(nneighmax)  ! hartree   ! new 9/3/04
      dimension barrierev(nneighmax)
      dimension barrevev(nneighmax)
      dimension prefac(nneighmax)
      dimension nattempted(nneighmax)  ! xx careful here - nattempt is different
      dimension nnewattempted(nneighmax)
      dimension naccepted(nneighmax)
      dimension ninvolved(nneighmax)
      dimension nebfound(nneighmax)
      dimension ntick(nneighmax)
      dimension ticknext(ntickmax,nneighmax)
      dimension isynthclass(nneighmax)

c equivalent-energy arrays
      parameter(nequivmax=1000)
      dimension eequiv(nequivmax)
      dimension bminequiv(nequivmax)
      dimension dimbminevequiv(nequivmax)
      dimension tprevequiv(nequivmax)

      character*80 ctitle(20)

      parameter(maxtyp=5)
      dimension rcut(maxtyp,maxtyp),amu(maxtyp),azero(maxtyp)
      character*80 potnam(maxtyp)
      character*80 fmt
      !character*80 potdirnam,potfilnam,potcomb

      parameter(maxtherm=10)
      dimension itherm(maxtherm),jtherm(maxtherm)
      dimension ttherm(maxtherm),ctherm(maxtherm)
      dimension xtherm(maxtherm),mtherm(maxtherm),ktherm(maxtherm)

      common/uflags/iuflag(100)
      common/trpcom/mutetrp
      common/callcom/ngradcall

      common/bpuds/dt,amu,tqtempmax,tqrate,ntype ! xxxb i want to pass these, but want to wait until i have everything i need before doing so

      ! Setup MPI/SpecTAD
      irecognize = 1 !! Need to change to zero to NOT allocate statedata memory
      call startup_mpi(nofficialmax,irecognize)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)

      !! Set rank,nprocs,current_comm to global settings for setup...
      !! ... Will switch to local settings when starting tad
      !call mpiset_global()
      !call MPI_BARRIER(MPI_COMM_WORLD,ier)
      xyzofficial(:) = 0.d0 ! Zero-out official state coordinates

      if(rank_w.eq.0)then
      write(6,*) ' '
      write(6,*) '--------------------------------------------'
      write(6,*) '                begin FORTRAN TAD run       '
      write(6,*) '          (version TAD 3.43 - Apr 10, 2009) '
      write(6,*) '--------------------------------------------'

      write(6,*)'**************************************************'
      write(6,*)'This SOFTWARE has been authored  by an employee or'
      write(6,*)'employees  of    the  University of    California,'
      write(6,*)'operator   of the  Los  Alamos National Laboratory'
      write(6,*)'under  Contract  No.  W-7405-ENG-36 with  the U.S.'
      write(6,*)'Department  of  Energy.   The U.S. Government  has'
      write(6,*)'rights  to  use,  reproduce,  and  distribute this'
      write(6,*)'SOFTWARE.  The public may copy, prepare derivative'
      write(6,*)'works  and publicly  display this SOFTWARE without'
      write(6,*)'charge,  provided   that   this  Notice   and  any'
      write(6,*)'statement  of  authorship  are  reproduced  on all'
      write(6,*)'copies.  Neither the Government nor the University'
      write(6,*)'makes any warranty, express or implied, or assumes'
      write(6,*)'any liability  or responsibility  for  the  use of'
      write(6,*)'this SOFTWARE.  If SOFTWARE is modified to produce'
      write(6,*)'derivative works, such modified SOFTWARE should be'
      write(6,*)'clearly marked, so as not  to confuse it with  the'
      write(6,*)'version available from LANL.'
      write(6,*) ' '
      write(6,*)'It is expressly forbidden to distribute this code,'
      write(6,*)'in whole or in part, either as is or modified.'
      write(6,*)'**************************************************'
      write(6,*) ' '
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ier)

c open a standard file for parameter data input
      fmt='(a80)'

      if(ntickmax.ne.ntickmax2)then
        STOP 'ntickmax.ne.ntickmax2 !!'
      endif
      if(nneighmax.ne.nneighmax2)then
        STOP 'nneighmax.ne.nneighmax2 !!'
      endif
      if(mxatom.ne.mxatom2)then
        STOP 'mxatom.ne.mxatom2 !!'
      endif

      if(rank_w.eq.0) read(5,fmt) casename   ! 'te st1' !  Ag/Ag(100), T=400K expansion,
      call MPI_BCAST(casename,80,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,fmt) potnam(1)  ! 'ag1'
      call MPI_BCAST(potnam(1),80,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,fmt) potnam(2)  ! 'xx'
      call MPI_BCAST(potnam(2),80,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) naccept         ! 10  ! number of accepted transitions at which run stops
      call MPI_BCAST(naccept,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) tadtimestop  ! total time after which to stop TAD
      call MPI_BCAST(tadtimestop,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) templow      !  300.
      call MPI_BCAST(templow,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) temphigh     ! 500.
      call MPI_BCAST(temphigh,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      !temphigh_global = temphigh
      if(rank_w.eq.0) read(5,*) ivariabletemp   ! adjust temphigh based on minimum barrier
      call MPI_BCAST(ivariabletemp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) temphighmax     ! maximum high temperature if in variable high temperature mode
      call MPI_BCAST(temphighmax,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) prefacmin    ! 1.d11
      call MPI_BCAST(prefacmin,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) uncertainty  ! 0.01
      call MPI_BCAST(uncertainty,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) ietype       ! 0
      call MPI_BCAST(ietype,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) natom        ! 91
      call MPI_BCAST(natom,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      natomstart = natom
      if(rank_w.eq.0) read(5,*) nmove        ! 55
      call MPI_BCAST(nmove,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) nmd          ! 1000 ! number of MD steps per block
      call MPI_BCAST(nmd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) nmdsub       ! 100  ! number of sub-divisions of nmd block
      call MPI_BCAST(nmdsub,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) nwarm        ! 1000 ! number of MD steps whenever warming from a minimum
      call MPI_BCAST(nwarm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) nthermalize  ! 1000 ! number of thermalization steps in a new state
      call MPI_BCAST(nthermalize,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) nblack          ! 1000 ! number of blackout steps after rejecting attempted transition
      call MPI_BCAST(nblack,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) ibadignore      ! ignore failed blackout or warmups
      call MPI_BCAST(ibadignore,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) dt           ! 2.d-15   ! MD time step (sec)
      call MPI_BCAST(dt,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) thermrate    ! 1.d12    ! Langevin coupling constant (sec^-1)
      call MPI_BCAST(thermrate,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) iseed        ! 0
      call MPI_BCAST(iseed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) nsteep1      ! 500  ! number of SD steps at each new basin
      call MPI_BCAST(nsteep1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) nsteep2         ! 500  ! max no. of SD steps to check for transition
      call MPI_BCAST(nsteep2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) ifundc          ! minimization method for descentchecks (0=sd, 1=quickmin, 2=cg, 3=nordlund cg, 4=henkelman cg, 5=thermal quench)
      call MPI_BCAST(ifundc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) tqtempmax       ! maximum temperature in thermal quench (ifundc=5 only)
      call MPI_BCAST(tqtempmax,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) tqrate          ! temperature quench rate (ifundc=5 only)
      call MPI_BCAST(tqrate,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) ivineyard       ! flag to calculate vineyard prefactors for transitions
      call MPI_BCAST(ivineyard,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(ivineyard.eq.0)then
        lsynthattempts = .true.
      else
        lsynthattempts = .false.
      endif
      if(rank_w.eq.0) read(5,*) irotalign    ! 0    ! turn this on to rotate system into alignment
      call MPI_BCAST(irotalign,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) gfac         ! 1.0d0  ! steepest-desc. factor
      call MPI_BCAST(gfac,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) transcrit       ! 1.0d0  ! transition displacement criterion (angstroms)
      call MPI_BCAST(transcrit,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) itranscrit      ! type of transition detection to do (see notes below)
      call MPI_BCAST(itranscrit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) ndimer          ! number of dimer searches to do in all unseen states to generate neighbor list
      call MPI_BCAST(ndimer,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) ndimeretad      ! number of dimers to do to find minimum barrier for dimeretad
      call MPI_BCAST(ndimeretad,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) idimerdof       ! degrees of freedom selection method
      call MPI_BCAST(idimerdof,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) ndimerdof       ! number of degrees of freedom to include in dimer searches
      call MPI_BCAST(ndimerdof,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) rdimerdist      ! distance to count neighbors in for dimer
      call MPI_BCAST(rdimerdist,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) ndimercoord     ! coordination to define over or under coordination
      call MPI_BCAST(ndimercoord,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) idimeroverunder ! if dimer should displace over (0), under (1), or both (2)
      call MPI_BCAST(idimeroverunder,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     +  ier)
      if(rank_w.eq.0) read(5,*) rdimerdisp      ! size of displacements in dimer
      call MPI_BCAST(rdimerdisp,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) lstore_min      ! store <casename>.min.<istate>.dat files
      call MPI_BCAST(lstore_min,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) lstore_minunpbc ! store <casename>.min.unpbc.<istate>.dat files
      call MPI_BCAST(lstore_minunpbc,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     +  ier)
      if(rank_w.eq.0) read(5,*) lstore_attempt  ! store <casename>.<istate>.attempt.<iattempt>.dat files
      call MPI_BCAST(lstore_attempt,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) lstore_sad      ! store <casename>.<istate>.sad.<isad>.dat files
      call MPI_BCAST(lstore_sad,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) lstore_end      ! store <casename>.<istate>.end.<isad>.dat files
      call MPI_BCAST(lstore_end ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) lstore_dimer    ! store <casename>.<istate>.dimer.<idimer>.dat files
      call MPI_BCAST(lstore_dimer,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) lstore_neb      ! store chs files for all new transitions
      call MPI_BCAST(lstore_neb,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) lstore_attempt_chs ! store snapshots of trajectories for attempted transitions in chs format
      call MPI_BCAST(lstore_attempt_chs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     +  ier)
      if(rank_w.eq.0) read(5,*) ipr          ! print level (1 or 2 is a good choice)
      call MPI_BCAST(ipr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) dvcrit       ! 5.d-8 is good for eam - keep these tight to prevent balancing on saddle
      call MPI_BCAST(dvcrit,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) gcrit        ! 1.d-4 is good for eam
      call MPI_BCAST(gcrit,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) drcrit       ! 1.d-5 is good for eam
      call MPI_BCAST(drcrit,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) itad         ! what kind of tad to do
      call MPI_BCAST(itad,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      parallelmode = 'none'
      if(itad.eq.100)then
          itad = 0
          if(.not.execute_serial)then
            parallelmode = 'spectad'
            tadtimestopspawn(:) = tadtimestop
            avgtadtimestop = abs(tadtimestop)
            tadtimestop = 1.d30
          endif
      else
          if(.not.execute_serial) then
              STOP 'itad=100 is only parallel setting supported!'
          endif
      endif
      if(rank_w.eq.0) read(5,*) irecognize   ! what kind of state recognition to do
      call MPI_BCAST(irecognize,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      !if((parallelmode.eq.'spectad').and.(irecognize.gt.1))then
      ! STOP 'irecognize must be <= 1 for spectad !!!!!'
      !endif
      if(rank_w.eq.0) read(5,*) ireadstates     ! read stored state files
      call MPI_BCAST(ireadstates,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) nimage       ! number of NEB images
      call MPI_BCAST(nimage,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) ifunneb         ! minimizing method for NEB (0=sd, 1=quickmin)
      call MPI_BCAST(ifunneb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) itan         ! tangent method
      call MPI_BCAST(itan,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) iclimb       ! use climbing method
      call MPI_BCAST(iclimb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) springk      ! spring constant for NEB
      call MPI_BCAST(springk,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) eshallow     ! maximum barrier to not resolve
      call MPI_BCAST(eshallow,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) hotmix       ! fraction of hot trajectory to mix into nudged elastic band guess
      call MPI_BCAST(hotmix,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) inebdimer    ! launch dimer from top of NEB?
      call MPI_BCAST(inebdimer,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) inebirc         ! ignore roll checks in NEB (if going to fail often)
      call MPI_BCAST(inebirc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(rank_w.eq.0) read(5,*) ereverse        ! reverse barriers to ignore transitions to
      call MPI_BCAST(ereverse,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      call chrparse1(casename,80,lcase)
      call chrparse1(potnam(1),80,lpot1)
      call chrparse1(potnam(2),80,lpot2)

      if(rank_w.eq.0) then
      write(6,*) '-----------------'
      write(6,*) 'input parameters:'
      write(6,*) '-----------------'
      write(6,*) 'casename =',casename(1:lcase)
      write(6,*) 'POTENTIAL 1 =',potnam(1)(1:lpot1)
      write(6,*) 'POTENTIAL 2 =',potnam(2)(1:lpot2)
      write(6,*) 'naccept =',naccept,' (no. of transitions to accept)'
      write(6,*) 'tadtimestop = ',tadtimestop,' (time to stop TAD)'
      write(6,*) 'templow =',templow,' (low T for TAD)'
      write(6,*) 'temphigh =',temphigh,' (high T for TAD)'
      write(6,*) 'ivariabletemp =',ivariabletemp,' (variable high T?)'
      write(6,*) 'temphighmax =',temphighmax,' (max high T)'
      write(6,*) 'prefacmin =',prefacmin,' (numin)'
      write(6,*) 'uncertainty =',uncertainty,' (delta (e.g., 0.01))'
crz -> Adding type 1xx potentials -- Lammps Potentials (100=EAM/Alloy)...
      write(6,*) 'ietype=',ietype
     +    ,' (0=eam,27=Tersoff Si,C,Ge,67=AIREBO, 100=Lammps-EAM/ALLOY)'
      write(6,*) 'natom =',natom, '(number of atoms)'
      write(6,*) 'nmove =',nmove, '(number of moving atoms)'
      write(6,*) 'nmd =',nmd,' (number of MD steps per block)'
      write(6,*) 'nmdsub =',nmdsub,' (number of MD steps per subblock)'
      write(6,*) 'nwarm =',nwarm,' '
      write(6,*) 'nthermalize =',nthermalize,' '
      write(6,*) 'nblack =',nblack,' '
      write(6,*) 'ibadignore =',ibadignore
     +    ,' (ignore failed blackout or warmup)'
      write(6,*) 'dt =',dt,' (MD time step (seconds))'
      write(6,*) 'thermrate =',thermrate,' (Langevin coupling rate)'
      write(6,*) 'iseed =',iseed,' (zero means pick random values)'
      write(6,*) 'nsteep1 =',nsteep1,' (at each fresh start)'
      write(6,*) 'nsteep2 =',nsteep2,' (max during event checks)'
      write(6,*) 'ifundc =',ifundc,
     +    ' (minimization method fordescentchecks)'
      write(6,*) 'tqtempmax =',tqtempmax
     +    ,' (max temp for thermal quench (ifundc=5 only))'
      write(6,*) 'tqrate =',tqrate,' (temp quench rate (ifundc=5 only))'
      write(6,*) 'ivineyard =',ivineyard
     +    ,' (calculate vineyard prefactors)'
      write(6,*) 'irotalign =',irotalign,' (rotate into alignment)'
      write(6,*) 'gfac =',gfac,' (steepest descent factor)'
      write(6,*) 'transcrit=',transcrit,
     x                     ' (transition distance crit. per atom (A))'
      write(6,*) 'itranscrit=',itranscrit,
     x                     ' (how to detect transitions)'
      write(6,*) 'ndimer =',ndimer
     +    ,' (number of dimer searches to generate neighbor list)'
      write(6,*) 'ndimeretad = ',ndimeretad
     +    ,' (number of dimers for minimum barrier for dimeretad)'
      write(6,*) 'idimerdof =',idimerdof,
     +    ' (degrees of freedom selection method)'
      write(6,*) 'ndimerdof =',ndimerdof,
     +    ' (number of degrees of freedom to include in dimer searches)'
      write(6,*) 'rdimerdist =',rdimerdist
     +    ,' (distance to define coordination in dimer)'
      write(6,*) 'ndimercoord =',ndimercoord,
     +    ' (coordination definingover or under coordination)'
      write(6,*) 'idimeroverunder =',idimeroverunder
     +    ,' (0=over, 1=under, 2=both coordination for dimer)'
      write(6,*) 'rdimerdisp =',rdimerdisp
     +    ,' (size of displacements in dimer)'
      write(6,*) 'lstore_min=',lstore_min,
     x    ' (store <case>.min.<istate>.dat files)'
      write(6,*) 'lstore_minunpbc=',lstore_minunpbc,
     x             ' (store <case>.min.unpbc.<istate>.dat files)'
      write(6,*) 'lstore_attempt=',lstore_attempt,
     x           ' (store <case>.<istate>.attempt.<iattempt>.dat files)'
      write(6,*) 'lstore_sad=',lstore_sad,
     x             ' (store <case>.<istate>.sad.<ipath>.dat files)'
      write(6,*) 'lstore_end=',lstore_end,
     x             ' (store <case>.<istate>.end.<ipath>.dat files)'
      write(6,*) 'lstore_dimer=',lstore_dimer,
     x             ' (store <case>.<istate>.dimer.<idimer>.dat files)'
      write(6,*) 'lstore_neb=',lstore_neb
     +    ,' (store chs files for new transitions)'
      write(6,*) 'lstore_attempt_chs=',lstore_attempt_chs,
     +    ' (store trajsnapshots for all attempted transitions)'
      write(6,*) 'print level= ',ipr
      write(6,*) 'dvcrit =',dvcrit,
     x      ' (delta energy conv. crit. (hartree) for descent checks)'
      write(6,*) 'gcrit =',gcrit,'  (total gradient',
     x 'conv. crit. (hartree/Angstrom) for descent checks)'
      write(6,*) 'drcrit =',drcrit,
     x ' (atom move conv. crit. (Angstrom) for descent checks)'
      write(6,*) 'itad =',itad, ' (what kind of tad to do)'
      write(6,*) 'irecognize =',irecognize,
     +    ' (what kind of state recognition to do)'
      write(6,*) 'ireadstates =',ireadstates,' (read state files)'
      write(6,*) 'nimage =',nimage, ' (NEB: number of images in chain)'
      write(6,*) 'ifunneb =',ifunneb,' (NEB: minimization method)'
      write(6,*) 'itan =',itan, ' (NEB: tangent method to use)'
      write(6,*) 'iclimb =',iclimb, ' (NEB: use climbing image, 1=yes)'
      write(6,*) 'springk =',springk,' (NEB: spring constant)'
      write(6,*) 'eshallow =',eshallow,
     +    ' (NEB: maximum barrier to not resolve)'
      write(6,*) 'hotmix =',hotmix,' (NEB: fract. of hot traj. in NEB)'
      write(6,*) 'inebdimer =',inebdimer,' (NEB: dimer-NEB)'
      write(6,*) 'inebirc =',inebirc,' (ignore roll checks in NEB)'
      write(6,*) 'ereverse =',ereverse,' (reverse barriers to ignore)'
      endif !(rank_w.eq.0)

c----------------------------------------------------------------------
c  parameters to later add to input file
      maxwarmreject=100
      maxblackreject=100
      if (ibadignore.eq.1) then
         maxwarmreject=10 !4   !!!!! HACK !!!! Should be larger to be safer...
         maxblackreject=10 !4
         if (ietype.eq.72.or.ietype.eq.73.or.ietype.eq.75.or.
     +       ietype.eq.80.or.ietype.eq.172) then
            maxwarmreject=5 !0 ! 0 is too dangerous
            maxblackreject=5 !0 ! 0 is too dangerous
         endif
      endif
      maxwarmrejectsave=maxwarmreject
      maxblackrejectsave=maxblackreject
      fixedprefac=1.d12  ! used if itad=5
      fixedbarriermin=0.0d0  ! used if itad=3
      iequivcheck=0             ! use energy equivalence for recognizing states
      eequivcrit=1.d-5   ! hartrees -  put into data?
      ratioet=3.0d0             ! scale factor for variable temp high mode
      if((parallelmode.eq.'spectad').and.(ivineyard.eq.1))then
        synthrequirement='always'
      else
cc      synthrequirement='barrierheight'
        synthrequirement='nattempts'
cc      synthrequirement='firstchosen'
        if(ivineyard.gt.0) synthrequirement='always'
      endif
      synthbarrier= 0.3d0 ! used if synthrequirement='barrierheight'
      nattemptsynth=nattemptkMC ! used if synthrequirement='nattempts'
      reordercrit=1.d-2  ! sameness criterion in check for atom reordering matchup
      idepflag=1                ! flag for deposition, to modify natom based on cluster
      nrunninglimit = 1000 ! max number of running state files to write (good for synth mode)
      rdynmatcrit=0.01d0        ! criteria for selecting atoms for vineyardsubset
      idynmatloop=0
      idosupersynth=0           ! do super synthetic mode
      nacceptsuper=2            ! number of events to accept before putting state into super-synthetic mode
      ratiosuper=50.0d0         ! number of events before accepting first non-synth event for super-synth mode
      erecognizecrit=1.d-6      ! used with irecognize=9, for official states, and for neighbors
      rrecognizecrit=4.d-3      ! used with irecognize=9, for neighbors - for hyperdistance match
      r4recognizecrit=4.d-3     ! used with irecognize=9, for neighbors - for hyperdistance match
      rdimerredocrit=0.01d0     ! criteria for restarting dimer count
      imodecheck=0              ! check neb saddles for negative modes
      if (ietype.eq.78) then
         imodecheck=1
         rmodemag=-1e10
      endif
      istorefileopt=0           ! optional file output in other formats (1=xyz)

      if(rank_w.eq.0) then
      write(6,*) ' '
      write(6,*) '--- hidden input parameters: ---'
      write(6,*) 'ibadignore =',ibadignore
      write(6,*) 'maxwarmreject =',maxwarmreject
      write(6,*) 'maxblackreject =',maxblackreject
      write(6,*) 'lstore_neb =',lstore_neb
      write(6,*) 'iequivcheck =',iequivcheck,' (energy-equiv check)'
      write(6,*) 'using eequivcrit =',eequivcrit,' h'
      write(6,*) 'ratioet =',ratioet
      write(6,*) 'synthrequirement =',synthrequirement
      write(6,*) 'synthbarrier =',synthbarrier
      write(6,*) 'nattemptsynth =',nattemptsynth
      write(6,*) 'reordercrit =',reordercrit
      write(6,*) 'idepflag = ',idepflag
      write(6,*) 'nrunninglimit = ',nrunninglimit
      write(6,*) 'rdynmatcrit = ',rdynmatcrit
      write(6,*) 'erecognizecrit = ',erecognizecrit,
     x                      ' (used with irecognize=9)'
      write(6,*) 'rrecognizecrit = ',rrecognizecrit,
     x                      ' (used with irecognize=9)'
      write(6,*) 'r4recognizecrit = ',r4recognizecrit,
     x                      ' (used with irecognize=9)'
      write(6,*) ' '
      endif !(rank_w.eq.0)
c----------------------------------------------------------------------

      if(rank_w.eq.0) then
      if(iequivcheck.ne.0) then
c        stop 'this needs to be debugged!'
        do i=1,10
        write(6,*) '********* NOTE: iequivcheck turned on!!!!********'
        end do
      end if
      if(ireadstates.ne.0) then
        do i=1,10
        write(6,*) '********* NOTE: ireadstates turned on!!!!********'
        end do
      end if
      if(itad.eq.4 .or. itad.eq.5) then
        do i=1,10
        write(6,*) '********* NOTE: trying out P-TAD!!!!********'
        end do
      end if
      if(itad.eq.6) then
        do i=1,10
        write(6,*) '******** NOTE: trying out synthtad mode!!!!*******'
        end do
      end if
      endif !(rank_w.eq.0)

ccc----------------------------------------------------------------------
ccc check for lock file
cc      filnam=casename(1:lcase)//'.lock'
cc      call chrpak(filnam,180,lfil)
cc      inquire(file=filnam(1:lfil),exist=lexist)
cc      if (lexist) then
cc         write(6,*) 'STOP: lockfile ',filnam(1:lfil),' exists!'
cc         write(6,*)
cc     +       '      if you want to run this casename, remove lockfile'
cc         stop 'lockfile <casename>.lock exists, stopping run'
cc      else
cc         open(unit=12,file=filnam(1:lfil),form='formatted',
cc     +       status='new')
cc         write(12,*) 'lockfile: ',filnam
cc         close(12)
cc      endif
ccc----------------------------------------------------------------------

c notes on some parameters:   NOTE NOTE 9/24/02 - some of this is still in progress or even wrong...

c     itranscrit: determines how transitions are detected.
c     itranscrit=0 - use an absolute difference between configurations, 1
c     itranscrit=1 - use coordination, one cutoff radius
c     itranscrit=2 - use coordination, but for multicomponent system,
c                    use rcut array for cutoff radius
c     itranscrit=3 - use fullbonding for state definition

c     irecognize:   determines how states are recognized.
c     irecognize=0 - no recognition is performed
c     irecognize=1 - compare absolute positions to see if states are identical
c     irecognize=2 - also try reordering atoms to see if there is a match - new 4/15/03 xxxa
c     irecognize=9 - simply match up to any official state with same energy.
c                    NOTE:  this will give a screwy looking state trajectory and will not
c                    work if properties like a diffusion constant are desired
c                    NOTE:  this setting does not use any of the iequivcheck stuff,
c                    which is meant to be less brutal.
c                    NOTE:  irecognize=9 option is not very graceful right now.  It simply
c                    pretends there is a nice official state match, dumps the system into
c                    that state, and proceeds.

c     iequivcheck
c     iequivcheck=0 - no energy equivalence recognition
c     iequivcheck=1 - compare energy to see if states are equivalent
c                    (or equivalent enough).  This could be very powerful within
c                    itad=1 or itad=2 method, since the same minimum barrier can be used for all
c                    energy-equivalent states.  Just starting to try this out.  Use with caution.

c     itad:  this determines what mode tad is running in.
c     itad=0 - Original TAD, as described in Sorensen and Voter, JCP, 2000.
c     itad=1 - TAD with revisit enhancement, as described in Montalenti and Voter, JCP 2002,
c              This ultimately converges on minimum-barrier TAD for states visited many times.
c              After this minimum-barrier convergence, even energy-equivalent states can benefit
c              from the ETAD speedup, if iequivcheck=1 (not yet fully implemented).
c     itad=2 - Dimer-based minimum-barrier method.  Note that this is approximate
c              since the random dimer searches are not guaranteed to find the lowest barrier.
c              See Henkelman and Jonsson, JCP 1999, JCP 2001,  and Montalenti and Voter, JCP 2002.
c              - requires ndimer>0 (working name: dimeretad)
c     itad=3 - Use minimum-barrier method with a fixed value for the Emin
c              - requires nonzero value for fixedbarriermin
c              This can be powerful for certain types of studies in which the
c              lowest possible barrier is already known for all states.  However,
c              if this not actually true, then the dynamics will be wrong.
c     itad=4 - Prefactor-TAD (working name: ptad) using Vineyard prefactors
c              - requires ivineyard=1
c              - much more powerful if states are revisited, and irecognize is turned on
c     itad=5 - Prefactor-TAD using a fixed prefactor
c              - requires a nonzero value for fixedprefac
c              - much more powerful if states are revisited, and irecognize is turned on
c              - should also work with iequivcheck=1
c     itad=6 - synthetic mode tad "synthtad".  This combines regular tad with synthetic
c              mode.  There will be other synthetic based modes as well, but this is the
c              first try.
c     itad=7 - synthetic mode tad "synthdimeretad".  This combines dimer-based minimum-barrier tad
c              with synthetic mode.
c     itad=8 - do regular TAD for first state, then dimer-TAD reusing saddles from first state

cc      if(naccept.gt.maxtrans) ! xxxx think we got rid of need for this, but save just in case
cc     x                   stop 'naccept.gt.maxtrans'
c xxx ivineyard=2 for       always, 1 for only when needed? - proposal of 4/03
      if(ivineyard.lt.0 .or. ivineyard.gt.3)
     x                   stop 'unallowed ivineyard'
      if(itad.gt.1)
     x             write(6,*) 'CAUTION - itad>1 is under construction'
      if(itad.eq.1 .and. irecognize.le.0)
     x                   stop 'itad=1,irecognize<=0'
      if(itad.eq.2 .and. ndimeretad.le.0)
     x                   stop 'itad=2 and ndimeretad.le.0'
      if(itad.eq.3.and.fixedbarriermin.le.0.0d0)
     x                   stop 'itad=3 and fixedbarriermin.le.0.0d0'
      if(itad.eq.4 .and. ivineyard.le.0)
     x                   stop 'itad=4 and ivineyard.le.0'
      if(itad.eq.5 .and. fixedprefac.le.0.0d0)
     x                   stop'itad=5,fixedprefac=0'
      if (itad.eq.8 .and. ndimeretad.le.0)
     x                   stop'itad=8 and ndimeretad.le.0'
      if (itad.eq.9 .and. ndimeretad.le.0)
     x                   stop'itad=9 and ndimeretad.le.0'
      if(itad.gt.9) stop 'unallowed itad value'
      if(itad.lt.0) stop 'unallowed itad value'

      tadmode='unknown'
      ratemode='unknown'
      if(itad.eq.0) tadmode='tad'
      if(itad.eq.1) tadmode='newtad'
      if(itad.eq.2) tadmode='dimeretad'
      if(itad.eq.3) tadmode='etad'
      if(itad.eq.4) tadmode='ptad'
      if(itad.eq.5) tadmode='ptad'
      if(itad.eq.6) tadmode='synthtad'
      if(itad.eq.7) tadmode='synthdimeretad'
      if(itad.eq.8) tadmode='protecteddimeretad'
      if(itad.eq.9) tadmode='userdimeretad'

c -------notes on new synthetic mode stuff, with some notes and comments--------

c We should probably later change ratemode to be it's own flag at input time,
c rather than one that gets set by what value itad has,
c and have as a new possible value for ratemode:  'attemptrate'.

c This will be like sorensen's synthetic mode if it is used in conjunction with a synthetic mode
c like synthtad, although at first we don't have the attempt-counting option to determine the
c prefactor.   FOR NOW, ratemode will be vineyard for this option.

cc============= xxxx synth do list:=============
c [done] put in basic synthtad mode
c [done] pass isynthclass to state file
c [done] figure out how to decide if last transition was synthetic or not
c [done] think about how to do the synthdimeretad
c [partially done] limit printout for synthetic visits
c [done] put in nattempts option for promotion
c [] put in nattempts-based option for prefactor determination
c [] debug all this!!!

      if(itad.eq.0 .or. itad.eq.1 .or. itad.eq.2 .or. itad.eq.4 .or.
     +    itad.eq.5 .or. itad.eq.6 .or. itad.eq.7 .or. itad.eq.8 .or.
     +    itad.eq.9 ) then
        continue
      else
        write(6,*) 'itad =',itad
        write(6,*) 'this value of itad not yet implemented'
        stop 'this value of itad not yet implemented'
      end if

      if(itad.eq.4) ratemode='vineyard'
      if(itad.eq.5) ratemode='fixedprefac'
      if(itad.eq.6) ratemode='nattempts'    ! xxx for now
      if(itad.eq.7) ratemode='vineyard'    ! xxx for now
      if(ivineyard.gt.0) ratemode='vineyard'

      if(rank_w.eq.0) then
        write(6,*) 'tadmode is being set to ',tadmode
        write(6,*) 'ratemode is being set to ',ratemode
      endif

      fstar=prefacmin/log(1.0d0/uncertainty)
      boltzev=8.617d-5  ! Boltzmann constant in eV/K
      betalow=1.0d0/(boltzev*templow)
      betahigh=1.0d0/(boltzev*temphigh)

      if(tadmode.eq.'ptad') then
        call ptadsetup(uncertainty,alpha,pmissing,pmatters)
        fpstar=prefacmin/log(1.0d0/pmissing)
      end if

c ========== check parameters for consistency ==============
      if (hotmix.ne.0.0d0 .and.irotalign.ne.0)
     x    stop 'cannot do hotmix.ne.0.0 and irotalign.ne.0 yet'

      if(ratemode.eq.'vineyard'.and.ivineyard.eq.0) then
         write(6,*) 'STOP: ratemode=vineyard but ivineyard=0'
         stop 'ratemode=vineyard but ivineyard=0'
      endif

      iok=0
      if(irecognize.eq.0) iok=1
      if(irecognize.eq.1) iok=1
      if(irecognize.eq.2) iok=1
      if(irecognize.eq.9) iok=1
      if(iok.eq.0) then
         write(6,*) 'stopping - bad value for irecognize:',irecognize
         stop 'bad value for irecognize:'
       end if

      if(ivariabletemp.ne.0) then
         if(tadmode.eq.'dimeretad'.or.tadmode.eq.'synthdimeretad'.or
     +       .tadmode.eq.'protecteddimeretad'.or.
     +       tadmode.eq.'userdimeretad') then
          continue
        else
          write(6,*)
     x   'STOP: variable temphigh only for dimeretad or synthdimeretad'
          stop 'variable temphigh only for dimeretad, synthdimeretad'
        end if
      end if

      if (idosupersynth.eq.1) then
         write(6,*) 'STOP: idosupersynth=1 not fully debugged'
         stop 'idosupersynth=1 not fully debugged'
      endif

      if (idosupersynth.eq.1.and.itad.ne.7) then
         write(6,*) 'STOP: to do super-synth mode, need itad=7'
         stop 'to do super-synth mode, need itad=7'
      endif

      if (idosupersynth.eq.1.and.iequivcheck.eq.1) then
         write(6,*)
     +       'STOP: super-synth mode and equivcheck not yet compatible'
         stop 'super-synth mode and equivcheck not yet compatible'
      endif

      if (iequivcheck.eq.1.and.(itad.ne.1.and.itad.ne.2))
     +    then
          write(6,*)
     +       'STOP: iequivcheck compatible only with itad=1'
          stop 'STOP: iequivcheck compatible only with itad=1'
c turned off working for itad=8, as there might be other things that
c need considered (want dimers to fill in anyways, even if match? also,
c need to know more than if iequiv.ne.0 to determine if need to do
c dimers, also need to know if this is a new equiv state or not)
      endif

      if (inebdimer.eq.1) then
         write(6,*) 'STOP: inebdimer not debugged'
         stop 'STOP: inebdimer not debugged'
      endif

      ibminspecialflag=0

      nhold=nmd/nmdsub
      if(nmdsub*nhold.ne.nmd) stop 'nmd/nmdsub is not an integer'
      nholdmax=nholdmaxworstcase*mxatom/natom   ! xx later, we could store just the moving atoms in this array
      if(rank_w.eq.0) write(6,*) 'nhold = nmd/nmdsub =',nhold
      if(rank_w.eq.0) write(6,*) 'for this natom, max nhold=',nholdmax
      if(nhold.gt.nholdmax) then
         write(6,*) 'STOP: increase nholdmaxworstcase'
         stop 'increase nholdmaxworstcase'
      endif

c      if (naccept+2.gt.maxtrans) then
c         write(6,*) 'STOP: naccept+2.gt.maxtrans, increase maxtrans'
c         stop 'naccept+2.gt.maxtrans, increase maxtrans'
c      endif

      if(irotalign.ne.0 .and. nmove.ne.natom) then
        write(6,*) 'DANGER - irotalign.ne.0 and nmove.ne.natom'
        stop 'aborting due to dangerous irotalign'
      end if

      iuflag(92)=10  ! iuflag(92)/100 is the extra neighbor list padding (angstroms)
      mutetrp=1

!c === seed the random number generator ===

      if(iseed.eq.0)then

        iseedlast = iseed
        do i = 0, nprocs_w-1
         if(rank_w.eq.i)then
           if(rank_w.eq.0)then
             call seedit(iseed)
             iseedlast = iseed
           else
             iseed = iseedlast+(rank_w+1)
             call seedit(iseed)
             iseedlast = iseed
           endif
           write(6,*) 'rank_w ',rank_w,' -- random number seed =',iseed
         endif
         call MPI_BCAST(iseedlast,1,MPI_INTEGER,i,MPI_COMM_WORLD,ier)
        enddo

      else

        iseed = iseed + 7*rank_w
        call seedit(iseed)
        write(6,*) 'rank_w ',rank_w,' -- random number seed =',iseed
        call MPI_BARRIER(MPI_COMM_WORLD,ier)

      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ier)

c === get Poisson stop time? ===    ! fixed by afv 8/04 (needed neg. sign)
      if(tadtimestop.lt.0) then
       tadtimestop=-tadtimestop*log(1.0d0-prngen(0))
       if(rank_w.eq.0) write(6,*)'Drew Poisson tadtimestop=',tadtimestop
      endif
      dodepositionrun = .false.
      if(tadtimestopspawn(1).lt.0.0)then
        dodepositionrun = .true.
        taureference = tadtimestopspawn(1) ! Note: taureference is negative here
        waittimelast = 0.d0
        do i = 1, ndepositions
         if(rank_w.eq.0)then
          tadtimestopspawn(i)
     +                 = taureference*log(1.0d0-prngen(0))+waittimelast
          waittimelast = tadtimestopspawn(i)
         endif
         call MPI_BCAST(tadtimestopspawn,ndepositions
     +                                ,MPI_REAL8,0,MPI_COMM_WORLD,ier)
        enddo

        ! Make shure the first deposition happens quickly (for
        ! restarting)
        shiftdeptime = tadtimestopspawn(1)
        do i = 1, ndepositions
          tadtimestopspawn(i) = tadtimestopspawn(i) - shiftdeptime
          tadtimestopspawn(i) = tadtimestopspawn(i) + nmd*dt !1.0d-12
        enddo

      endif
      currenttadtimestop = tadtimestopspawn(1)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)

c ===== read in cluster file =====
      iunit=25
      iwrt=1
      filnam=casename(1:lcase)//'.start'   ! note: rm case.*.dat will not delete this
      call chrpak(filnam,180,lfil)

      if(rank_w.eq.0) open(unit=iunit,file=filnam(1:lfil),
     +                                form='formatted',status='old')
      call getclsnew(natomx,ntitle,ctitle,xyz,pxyz,itype,
     x                           taxes,iunit,iwrt,ierr)
      if(rank_w.eq.0) close(unit=25)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
      call MPI_BCAST(itype,natom,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      if(natom.ne.natomx.and.idepflag.eq.0) then
        write(6,*) 'aborting - cluster natom not right',natom,natomx
        stop 'natom not right'
      else if (natom.ne.natomx.and.idepflag.eq.1) then
         if(rank_w.eq.0) write(6,*)
     +       'cluster natom does not match input, modifying input'
         nmove=nmove+(natomx-natom)
         natom=natomx
      end if
      ntype=imax(itype,natom)
      do i = 1, natom
        itypestart(i)=itype(i)
      enddo

c move all atoms into primary period (0<x<tax,0<y<tay,0<z<taz)
      call dmzer(cdxyz,3*mxatom)
      call pbccum(natom,xyz,taxes,cdxyz)

      ! Choose Deposition types now
      itypedeposit(:) = itype(1)
      if(dodepositionrun.and.(ntype.ge.2).and.(.not.spawn_lammps))then
        !xA_deposit = 0.00d0 ! MUST BE 0.00d0 if you want to deposit type-2 ! 1.00d0 for type-1
        !xB_deposit = 1.00d0 - xA_deposit
        !xC_deposit = 0.00d0
        if((ietype.eq.172).or.(ietype.eq.72))then
          !! Pyrochlore:
          !xA_deposit = 2.0/10.0
          !xB_deposit = 6.0/10.0
          !xC_deposit = 2.0/10.0
          !! Pyrochlore Choose 1 type:
          !!xA_deposit = 1.0
          !!xB_deposit = 0.0
          !!xC_deposit = 0.0
        endif
        if(rank_w.eq.0)then
          do i = 1, ndepositions
            xdraw = prngen(0)
            if(xdraw.lt.xA_deposit)then
              itypedeposit(i) = 1
            elseif(xdraw.lt.(xA_deposit+xB_deposit))then
              itypedeposit(i) = 2
            else
              itypedeposit(i) = 3
            endif
            write(6,*) 'itypedeposit(',i,') = ',itypedeposit(i)
            call flush(6)
          enddo
        endif
        call MPI_BCAST
     +    (itypedeposit,ndepositions,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
        !do i = 1, ndepositions
        !  if((itypedeposit(i).eq.2).and.(xA.gt.(0.1d0)))then
        !    itypedeposit(i) = 1
        !    exit ! Change 1st type-2 deposition since we are starting with a type-2 adatom
        !  endif
        !enddo
      endif

c ========= potential-specific setup calls =========
c note that amu array must be filled here

      luselammps = .false.
      call MPI_BARRIER(MPI_COMM_WORLD,ier)

      rank_log = rank_w
      if((ietype.ge.100).and.(ietype.lt.200))then

        ! First let master in
        if(rank_w.eq.0)then
           rank_log = 0
           call MPI_BARRIER(force_comm,ier)
           call lmp_setup(natom,nmove,ntype,maxtyp,potnam,ietype,xyz,
     +       taxes,itype,rcut,amu,templow,temphigh,.false.,dt)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,ier)
        !write(*,*) 'rank_w ',rank_w,' Through A'
        ! Now do each replica individually
        do i = 1,lp_avail
         if(groupID.eq.i)then
          do j = 1, nParRep
            if((rank_l.ge.(j-1)*nforcecores)
     +                              .and.(rank_l.lt.j*nforcecores))then
              ! lmp setup must synchronize across force_comm...
              call MPI_BARRIER(force_comm,ier)
              rank_log = (groupID-1)*nParRep + j
              call lmp_setup(natom,nmove,ntype,maxtyp,potnam,ietype,xyz,
     +          taxes,itype,rcut,amu,templow,temphigh,.false.,dt)
            endif
            call MPI_BARRIER(local_comm,ier)
          enddo
         endif
         call MPI_BARRIER(MPI_COMM_WORLD,ier)
        enddo
        !call MPI_BARRIER(MPI_COMM_WORLD,ier)
        !write(*,*) 'rank_w ',rank_w,' Through B'
        ! Now state (file) master
        if(rank_w.eq.(nprocs_w-1))then
           rank_log = 1 + lp_avail*nParRep
           call MPI_BARRIER(force_comm,ier)
           call lmp_setup(natom,nmove,ntype,maxtyp,potnam,ietype,xyz,
     +       taxes,itype,rcut,amu,templow,temphigh,.false.,dt)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,ier)
        luselammps = .true.
      elseif(nprocs_f.gt.1) then
          write(*,*)
     +        'ERROR -- only lammps can handle a parallel force call!'
          stop 'Only lammps can handle a parallel force call!'
      endif

      do i = 0,nprocs_w-1
      if(rank_w.eq.i)then

      if(ietype.eq.0) then
c       set up EAM potential
c       (note that currently ietype is set by user, and also by setup routine)
        iwrt=0
        call eam_setup(ntype,maxtyp,potnam,iwrt,ietype,rcut,amu,azero)
        if(ietype.ne.ietype) stop 'parmd confused by ietype'
!crz -- Lammps eam/alloy
!      else if((ietype.eq.100).or.(ietype.eq.101).or.(ietype.eq.102)
!     +                       .or.(ietype.eq.103).or.(ietype.eq.172))then
!      ! ietype=100 -> EAM
!      ! ietype=101 -> EAM/ALLOY
!         call lmp_setup(natom,nmove,ntype,maxtyp,potnam,ietype,xyz,
!     +     taxes,itype,rcut,amu,azero,iseed,templow,temphigh,dt)
!         luselammps = .true.
      else if((ietype.ge.100).and.(ietype.lt.200))then
         if(rank_w.eq.0) write(*,*) 'USING LAMMPS!'
      else if (ietype.eq.27) then
         call prmst27(39,maxtyp,rcut,amu,azero)
      else if (ietype.eq.67) then
         write(6,*)
     +       'WARNING! Make sure <casename>.in is copied to parmd.in!'
         call prmst67(ntype,maxtyp,rcut,amu,azero)
      else if (ietype.eq.68) then
         call prmst68(ntype,maxtyp,rcut,amu,azero)
      else if (ietype.eq.69) then
         call prmst69(natom,ntype,itype,maxtyp,rcut,amu,azero,taxes)
      else if (ietype.eq.70) then
         call prmst70(natom,ntype,itype,maxtyp,rcut,amu,azero,potnam)
      else if (ietype.eq.72) then
         call prmst72(ntype,maxtyp,rcut,amu,azero,taxes)
      else if (ietype.eq.73) then ! TBMD (Marco Cogoni)
         itype73=2
         call prmst73(ntype,maxtyp,rcut,amu,itype73)
      else if (ietype.eq.74) then ! Stillinger Weber (Marco Cogoni)
         call prmst74(ntype,maxtyp,rcut,amu,azero)
      else if (ietype.eq.75) then ! Lenosky TB  (Marco Cogoni)
         itype75=2
         call prmst75(ntype,maxtyp,rcut,amu,itype75)
      else if (ietype.eq.76) then ! Si EDIP (Marco Cogoni)
         itype76=1
         amu(1) = 28.0855d0
         azero(1) = 5.42981d0
         rcut (1,1) = 1.0d0
         call prmst76(ntype,maxtyp,rcut,amu,azero)
      else if (ietype.eq.77) then ! ReaxFF
         call prmst77(ntype,natom,maxtyp,rcut,amu,azero,itype,taxes)
      else if (ietype.eq.78) then ! Siegel LJ
         call prmst78(ntype,maxtyp,rcut,amu,azero)
      else if (ietype.eq.80) then ! Buckingham + Morse
         call prmst80(ntype,maxtyp,rcut,amu,azero,taxes)
      else
         stop 'tad2 main: ietype not recognized'
      end if

      ! End of mpi-loop for setup:
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
      enddo

      !! NOW SWITCH TO LOCAL MPI SETTINGS...
      !call mpiset_local()
      !call MPI_BARRIER(MPI_COMM_WORLD,ier)

c check that masses are not zero
      do i=1,ntype
        if(amu(i).le.0.0d0) stop 'zero masses not allowed'
        if(rank_w.eq.0) write(6,*) 'Mass of type: ',i,' is ',amu(i)
      end do

c specify moving atoms
      do i=1,natom
         moving(i)=0
         mx(i)=0
         my(i)=0
         mz(i)=0
      enddo
      do i=1,nmove
         moving(i)=1
         mx(i)=1
         my(i)=1
         mz(i)=1
      enddo
      tadcharge(:) = 0.d0

c specify thermostat properties
c  for the i'th thermostat:
c    itherm(i)= start of atom range for this thermostat
c    jtherm(i)= end of atom range(-1=natom)
c    ttherm(i)= desired temperature(K)
c    ctherm(i)= response rate(/sec)
c    mtherm(i)= atom no. defining moving ref. frame
c    ktherm(i)= thermo type(0=off,1=Nose,2=Beren,3=frict,4=Langevin,5=Gauss)

      ntherm=1
      itherm(1)=1
      jtherm(1)=nmove
      ctherm(1)=thermrate
      ttherm(1)=temphigh
      mtherm(1)=0
      ktherm(1)=4

      call dmzer(pxyz,3*natom)

c quick gradient check xxx worth relegating to a subroutine?
      if(ipr.ge.4) then ! xxxx these gcalc not accounted for in summary
       if((parallelmode.ne.'spectad').or.(groupID.eq.1)) then
       if(rank_l.lt.nforcecores)then ! First replica only
        call MPI_BARRIER(force_comm,ier)
        call vecmov(xyz,xyz2,3*natom)
        i=3 ! atom 3
        j=2 ! y direction
        idx=j+3*(i-1)
        xyz2(idx)=xyz2(idx)-0.001d0
        call gcalc(natom,nmove,xyz2,itype, taxes,ietype,maxtyp,rcut,
     x                                          em,grad, hyperratio)
        xyz2(idx)=xyz2(idx)+0.001d0
        call gcalc(natom,nmove,xyz2,itype,taxes,ietype,maxtyp,rcut,
     x                                          e0,grad, hyperratio)
        xyz2(idx)=xyz2(idx)+0.001d0
        xyz2(idx)=xyz2(idx)+0.001d0
        call gcalc(natom,nmove,xyz2,itype,taxes,ietype,maxtyp,rcut,
     x                                          epp,grad, hyperratio)
        xyz2(idx)=xyz2(idx)-0.001d0
        call gcalc(natom,nmove,xyz2,itype, taxes,ietype,maxtyp,rcut,
     x                                          ep,grad, hyperratio)
        if(rank_f.eq.0) then
         write(6,*) 'grad check y3 a,n: ',grad(j,i),(epp-e0)/0.002d0
         write(6,*) '(energies: em,e0,ep,epp: ',em,e0,ep,epp,')'
         if (abs(grad(j,i)-(epp-e0)/0.002d0).gt.5e-5) then
           write(6,*)
     +         'STOP: analytical and numerical grads do not match'
c           stop 'analytical and numerical grads do not match'
         endif
        endif
       endif
       endif
      endif

      ngradcall=0  ! this is incremented in gcalc
      nsteep2sum=0
      nsteep2tries=0

      tadtime=0.0d0
      tadtimeprevious=0.0d0

      ! IF SpecTAD... set max_depth to naccept and naccept to 1:
      if(parallelmode.eq.'spectad')then
        max_depth = naccept
        naccept = 1000000 !! nofficialmax !! max_spawn !! Allow Loop Forever
      endif

      !CRZ -- MPI SETUP FOR SPECTAD -- # 2
      call setup_mpi2(nmove,natom,nofficialmax)
      call MPI_BCAST(activated,1,MPI_LOGICAL,0,local_comm,ier)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)


c open the .blok file and write header   - blok file new as of 9/1/04
c close and reopen to get rid of any earlier version
      filnam=trim(filstrtLP)//casename(1:lcase)//'.blok'
      call chrpak(filnam,180,lfil)
      if(rank_l.eq.0) then
       open(unit=42,file=filnam(1:lfil),form='formatted',
     x                                   status='unknown')
       rewind 42
       write(42,121) casename(1:lcase)
       close(42)

       open(unit=42,file=filnam(1:lfil),form='formatted',
     x                                   status='unknown')
       rewind 42
       write(42,121) casename(1:lcase)
      endif
  121 format('blok file - info about every MD block - case = ',a,/)


c open the .times file and write header
c close and reopen to get rid of any earlier version   xxxa   4/16/03
      filnam=trim(filstrtLP)//casename(1:lcase)//'.times'
      call chrpak(filnam,180,lfil)
      if(rank_l.eq.0) then
       open(unit=22,file=filnam(1:lfil),form='formatted',
     x                                   status='unknown')
       rewind 22
       write(22,123) casename(1:lcase)
       close(22)

       open(unit=22,file=filnam(1:lfil),form='formatted',
     x                                   status='unknown')
       rewind 22
       write(22,123) casename(1:lcase)
      endif

c open the .boost file and write header
c close and reopen to get rid of any earlier version   xxxa   4/16/03
      filnam=trim(filstrtLP)//casename(1:lcase)//'.boost'
      call chrpak(filnam,180,lfil)
      if(rank_l.eq.0) then
       open(unit=23,file=filnam(1:lfil),form='formatted',
     x                                   status='unknown')
       rewind 23
       write(23,1233) casename(1:lcase)
       close(23)

       open(unit=23,file=filnam(1:lfil),form='formatted',
     x                                   status='unknown')
       rewind 23
       write(23,1233) casename(1:lcase)
      endif

 1233 format('boost file:  estimated boost factor - case = ',a,/
     x '     i  istate jstate ioffcl Erel  ',
     x '    hightime  ',
     x '       nforce  TAD time         state boost   run boost',
     x '     cpu boost')

c open the .barriers file
      filnam=trim(filstrtLP)//casename(1:lcase)//'.barriers'
      call chrpak(filnam,180,lfil)
      if(rank_l.eq.0) then
       open(unit=12,file=filnam(1:lfil),form='formatted',
     x                                   status='unknown')
       rewind 12
       write(12,2343) casename(1:lcase)
       write(12,2344) 'state','iofficial','saddle','Ea(eV)','Erev(eV)',
     x           'Eend(eV)','hyperdist','4-dist',
     x           'prefactor(Hz)','timehigh','timelow','nmoved'
       rewind 12
       close(unit=12)
      endif
 2343 format('Barrier information - case = ',a)
 2344 format(a10,a10,a9,a9,a9, a12,a12,a12, a16,a14,a14,a12)

      cumtime=0.0d0

      blocktime=dt*float(nmd)
      if(rank_w.eq.0) write(6,125) blocktime
  125 format(/' length of one integration block (s) =',1pd16.8,/)

      if(rank_w.eq.0) cputime = mpi_wtime()-cput_ref
      call MPI_BCAST(cputime,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      cputimelast=cputime

      !ihackwait = 0
      !if(rank_w.eq.1) ihackwait = 1
      !do while(ihackwait.eq.1)
      !  !
      !enddo
      !call MPI_BARRIER(MPI_COMM_WORLD,ier)

      nofficial=0
      nequiv=0
      ieventlast=0
      nneighlast=0

c total gradcalls in each section
      ngradcallwarmtot=0           ! force calls for warming
      ngradcallblacktot=0          ! blackout
      ngradcalldimertot=0          ! dimers
      ngradcallnmdtot=0            ! MD loop
      ngradcallrefinetot=0         ! refining
      ngradcallminimizetot=0       ! minimizing basin
      ngradcalltranstot=0          ! detecting transition
      ngradcallnebtot=0            ! doing NEB stuff
      ngradcallleaktot=0           ! doing leak check
      ngradcallvineyardtot=0       ! vineyard - note that we also have cost of diagonalization

      ispecialstop=0
      iofficialmatch = 0
      isneighshare = 0
      ldprcheckdone = .false.
      dprcheckti = 10
      dprcheckt = dprcheckti !Use until we know the ~length of an MD block
      ldepspawn = 0

      if(rank_w.eq.0) timechecklast = mpi_wtime()-cput_ref
      call MPI_BCAST(timechecklast,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      if(groupID.eq.fmanagerID) cputime_ref = mpi_wtime()-cput_ref
      lnewstate = .false.
      lnewneigh = .false.
      ifmastart = 2

c ==================== Begin master loop over transitions that will be accepted ====================
      !itrans = 0
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
      do 500 itrans=1,naccept+1
        !itrans = itrans + 1

C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C
C~C~C~C                              START SPECTAD IF-STATEMENT                             C~C~C~C
C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C~C

        if(execute_serial) activated = .true.

        isleep_ma = 0
        killcnt = 0
        spawnIDlocal = 0 ! Current LOCAL spawn number
        isneighshare = 0
        iprsafe = .true.

        do while (.not.activated)

         killcnt = 0
         spawnIDlocal = 0 ! Current LOCAL spawn number
         iprsafe = .true.

         ! Need to restart the time for every spawn
         tadtime=0.0d0
         tadtimeprevious=0.0d0

         ! Spawn Manager (ma) Work:
         if(groupID.eq.managerID)then

           ! Check For Spawns to Initiate
           call ma_check_spawn(uncertainty)

           ! Check For Spawns to Kill (1)
           call ma_check_kill()

           ! Send Activate Messages for New Spawns
           call ma_activate_new(nofficialmax)

           ! Check Exit Condition:
           if((nsrunning.le.0).and.(nswaiting.le.0))then
             call ma_check_finish(isleep_ma)
           endif

           ! Update ma_spawnlist_os()
           call ma_check_iof()

           ! Check For Spawns that have finished
           call ma_check_done()

           ! If the run is set to end, increment isleep_ma:
           if((nsrunning.le.0).and.(nswaiting.le.0))then
             isleep_ma = isleep_ma + 1
           else
             isleep_ma = 0
           endif

         ! File Manager (fma) Work:
         elseif(groupID.eq.fmanagerID)then

           ! Update Global States:
           if(irecognize.gt.0)then
             ! Update Global States:
             call fma_state_update(taxes,ntype,itype
     +         ,erecognizecrit,lenxyzneighs,ietype
     +         ,maxtyp,rcut,nofficialmax
     +         ,irotalign,nsteep1,gfac,transcrit
     +         ,itranscrit,gcrit,dvcrit,drcrit,ipr,ifundc
     +         ,betalow,avgtadtimestop,lnewneigh,lnewstate)
           endif

           ! Check if SpecTAD is finished
           call lp_check_finish() ! Can just use LP's function

           ! Write out statelist.dat:
           cputime = mpi_wtime()-cput_ref
           fma_cputime = cputime - cputime_ref
           if(irecognize.gt.0)then
           if((lfinished).or.
     +               ((fma_cputime.gt.60.0).and.(lnewneigh)))then
             lnewstate = .false.
             open(UNIT=614,FILE='statelist.dat',STATUS='replace')
             write(614,*)
     +       'StateID E_Rel Neigh# NeighID Ea(Forward) Ea(Rev)'
             eminrel = 0.0
             if(fma_nofficial.ge.2)then
                 eminrel = statedata(2)%emin*27.21d0
             endif
             do j = 2, fma_nofficial
               jstate = j                             ! Global state label "j"
               knmove = statedata(j)%nmove
               knatom = knmove + natomfixed
               nneigh = statedata(j)%nneigh           ! Number of known neighbors to state "j"
               emin   = statedata(j)%emin*27.21d0     ! Minimum energy of state
               do k = 1, nneigh
                 kneigh = k                          ! Neighbor Number "k"
                 klabel = statedata(j)%neighlabel(k) ! State label for neighbor k
                 EbarF  = statedata(j)%barrierev(k)  ! Forward barrier to neighbor (eV)
                 EbarR  = statedata(j)%barrevev(k)   ! Reverse barrier from neighbor (eV)

                 !if((klabel.eq.0).and.(lnewstate.or.lfinished))then
                 if(klabel.lt.1)then

                   drmax_min = 1.d99
                   isinitial = 2
                   if(klabel.eq.(-1)) isinitial = ifmastart

                   do istate = isinitial, fma_nofficial
                     if ((knmove.eq.statedata(istate)%nmove)
     +                 .and.(istate.ne.jstate))then
                       call maxdr_rz(nmove
     +                   ,statedata(jstate)%xyzneighs(3*knatom*k+1)
     +                   ,statedata(istate)%xyzmin,drmax)
                       if ((drmax.lt.transcrit)
     +                               .and.(drmax.lt.drmax_min))then
                         drmax_min = drmax
                         klabel = istate
                       endif
                     endif
                   enddo

                   ! If this neigh is still '0', only compare to new states in future:
                   if(klabel.eq.0)then
                       statedata(j)%neighlabel(k) = -1
                       if(ipr.ge.1)then
                         ! Write out this neighbor -
                         filn=''
                         write(filn,*) jstate
                         filnb=''
                         write(filnb,*) kneigh
                         filn='stateneigh.'//trim(adjustl(filn))//
     +                        '.'//trim(adjustl(filnb))//'.dat'
                         pxyz(:) = 0.d0
                         ! Make sure the atom types are set correctly:
                         knatomdep = knatom - natomstart
                         icnt = 0
                         do idum = 1,knatom
                           if(idum.le.knatomdep)then
                             itype(idum) =
     +                         itypedeposit(knatomdep-idum+1)
                           else
                             icnt = icnt + 1
                             itype(idum) = itypestart(icnt)
                           endif
                         enddo
                         call storefile(knatom,
     +                        statedata(j)%xyzneighs(k*knatom*3+1),
     +                        pxyz,itype,taxes,filn)
                       endif
                   else
                     if((klabel.gt.0).and.
     +                    (statedata(j)%neighlabel(k).eq.(-1)))then
                      if(ipr.ge.1)then
                       ! Remove this file - It's not needed...
                       ! -> Only needed if a state.*.dat file doesn't
                       ! -> already exist for this geometry
                       filn=''
                       write(filn,*) jstate
                       filnb=''
                       write(filnb,*) kneigh
                       filn='rm stateneigh.'//trim(adjustl(filn))//
     +                      '.'//trim(adjustl(filnb))//'.dat'
                       call system(filn)
                      endif
                     endif
                     ! Update neighlabel:
                     statedata(j)%neighlabel(k) = klabel
                   endif

                 endif
                 !if(klabel.eq.(-1)) klabel = 0
                 write(614,
     +        '(I10,1X,ES16.6,1X,I10,1X,I10,1X,ES16.6,1X,ES16.6)')
     +           jstate,(emin-eminrel),kneigh,klabel,EbarF,EbarR
               enddo
             enddo
             close(614)
             cputime_ref = mpi_wtime()-cput_ref
             ifmastart = fma_nofficial
           endif
           endif

         ! Compute Node (lp) Work:
         else

          ! Check for new states to hold:
          if(irecognize.gt.0)then
            if(rank_l.lt.nforcecores)then
              call lp_update_official(nofficialmax)
            endif
            ! should create xyzofficial if not deposition run:
            call MPI_BARRIER(local_comm,ier)
          endif

          ! Check if SpecTAD is finished
          lfinished = .false.
          if(rank_l.eq.0) call lp_check_finish()
          call MPI_BCAST(lfinished,1,MPI_LOGICAL,0,local_comm,ier)
          if(lfinished)then
            write(*,*) 'rank_w ',rank_w,', rank_l ',rank_l
     +            ,', nprocs_w = ',nprocs_w,', nprocs_l = ',nprocs_l
     +            ,' - lfinished = ',lfinished,' - ier = ',ier
          endif

          ! Check if you are activated with a Spec-Spawn
          if(.not.lfinished)then
            iofficialmatch = 0
            isneighshare = 0 ! previous state to synthesize neighbors from
            call lp_check_start(natom,nmove,xyz,tadtime,itype
     +            ,ldepspawn,iofficialmatch,irecognize,nofficialmax
     +            ,taxes,isneighshare)
            tadtimeprevious = tadtime
            ! specify moving atoms
            do i=1,natom
                 moving(i)=0
                 mx(i)=0
                 my(i)=0
                 mz(i)=0
            enddo
            do i=1,nmove
                 moving(i)=1
                 mx(i)=1
                 my(i)=1
                 mz(i)=1
            enddo
            ntherm=1
            itherm(1)=1
            jtherm(1)=nmove
            ctherm(1)=thermrate
            ttherm(1)=temphigh
            mtherm(1)=0
            ktherm(1)=4
          endif

         endif

         tadcharge(:) = 0.d0
         call MPI_BCAST(lfinished,1,MPI_LOGICAL,0,local_comm,ier)
         if(lfinished)then
           write(6,*) ' GROUPID ',groupID,' IS EXITING!'
           goto 600
         endif

        enddo ! END of do while (.not.activated) Loop

        if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
        call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
        cputimelast=cputime

        ! Initialize some stuff for this spawn
        call lp_spec_init()

        ! Don't spawn if you are at/near max_spawn
        if(spawnID.ge.(max_spawn-1))then
          max_depth = spawnDepth
        endif

      !!
      !!
      !!
      !!   THIS SPAWN NOW STARTING TO RUN TAD...
      !!
      !!
      !!

c ==== initializations for this state ====

        ngradcallwarm=0        ! force calls for warming
        ngradcallblack=0       ! blackout
        ngradcalldimer=0       ! dimers
        ngradcallnmd=0         ! MD loop
        ngradcallrefine=0      ! refining
        ngradcallminimize=0    ! minimizing basin
        ngradcalltrans=0       ! detecting transition
        ngradcallneb=0         ! doing NEB stuff
        ngradcallleak=0        ! doing leak check
        ngradcallvineyard=0     ! vineyard - note that we also have cost of diagonalization

        inojump=0
        nneighfast=0
        lsenddone = .false.
        lputyes = .false.
        lputyessynthkilled = .false.
        lfirstinstate = .true.
        lsynthspawned = .false.
        nofficial = nstatedata
        wctimehighstop = 1.d99
        timehighstop = 1.d99
        lrepstold = .false.
        cleanwarmup = .true.
        cleanblackout = .true.

        ! Make sure the high temperature starts at the
        ! target high temperature
        !temphigh = temphigh_global

        idebug = 0
        do while(idebug.eq.1)
          call sleep(1)
        enddo

        ! Find currenttadtimestop
        icurrentevent = 0
        currenttadtimestop = tadtimestopspawn(1)
        do idum = 1, ndepositions-1
          if(tadtime.ge.tadtimestopspawn(idum))then
            icurrentevent = idum
            currenttadtimestop = tadtimestopspawn(idum+1)
          endif
        enddo
        skipdeposition = .false.
        if(icurrentevent.ge.(ndepositions-1))then
          skipdeposition = .true. ! Dont bother spawning the last event (wont run it anyway)
        endif

        ! Be careful not to overshoot the next deposition time...
        if(ldepspawn.eq.1)then
          if(tadtime+(ndepMDsteps*dt).ge.currenttadtimestop)then
            currenttadtimestop = tadtime + (ndepMDsteps*dt) + dt
          endif
        endif

        if(rank_l.eq.0) write(*,'(A,I10,1X,I9,1X,ES16.6,1X,ES16.6)')
     +    'spawnID,deposition#,tadtime,currenttadtimestop = '
     +    ,spawnID,icurrentevent,tadtime,currenttadtimestop

        ngradcallstate=ngradcall
        ngradcallhold=ngradcall

        cputimevineyard=0.0d0
        cputimereorder=0.0d0
        if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
        call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
        cputimestate=cputime
        istopping=0

        if(rank_l.eq.0)then
         if(ipr.ge.1) write(6,*) '====================================='
	    endif

        istate=itrans

c NOTE - most (all?) of the following will be overwritten if getstate is called, but
c if not, we need these zeroed out
        do idum=1,nneighmax
          esaddle(idum)=0.0d0
          eneigh(idum)=0.0d0
          barrierev(idum)=0.0d0
          prefac(idum)=0.0d0
          nattempted(idum)=0
          nnewattempted(idum)=0
          naccepted(idum)=0
          ninvolved(idum)=0
          nebfound(idum)=0
          ntick(idum)=0
          isynthclass(idum)=0
          do itick=1,ntickmax
            ticknext(itick,idum)=0.0d0
          end do
        end do
        nneigh=0

        ieventlast=0
        barrierminev=1.0d10
        dimerbarrierminev=1.d10
        barrierminevss=0.0d0    ! for super-synthetic mode
        barrierminevns=0.0d0    ! for super-synthetic mode
        isupersynth=0           ! for super-synthetic mode
        bminlowerbound=0.d0     ! for newtad mode
        timehighnow=0.0d0
        timehighnowI=0.d0
        timehighprev=0.0d0
        timelowprev=0.0d0
        tickprev=0.0d0
        freqlogsmin=0.0d0
        nnegmin = 1 ! Need to do vineyard in minimum if not a recognized state
        ineedblack=0
        ineedwarm=nwarm+nthermalize
        ineedgauss=1
        nattempt=0  ! total number of attempts to any neigh for all visits
        irevisited=0  ! as far as we know right here, this state has not been seen before
        ireordermatch=0  ! needed?   xxxxa 4/15/03
        ibadtransyes=0 ! clear any previous setting

c === if necessary, move all atoms back into primary period ===
c (neighbor lists should be safe if total span is .le. 1.5*tax)

        do ix=1,3
          call  dmnmxstride(xyz(ix),natom,xmin,xmax,ixmin,ixmax,3)
          if(xmax-xmin.gt.1.2d0*taxes(ix)) then
            call pbccum(natom,xyz,taxes,cdxyz)
          end if
        end do

        ! Get reference WCT time for Replicas
        call MPI_BARRIER(local_comm,ier)
        if(rank_f.eq.0) parRepRefTime = mpi_wtime()-cput_ref
        call MPI_BCAST(parRepRefTime,1,MPI_REAL8,0,force_comm,ier)
        call MPI_BARRIER(local_comm,ier)
        t_mas_tot_i = t_mas_tot

CCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCcc
CC
CC     ParRep Cores Can Loop HERE...
CC
CCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCcc
        if (rank_l.ge.nforcecores) then
          lprloop = .true.    ! ParRep core should stay in loop
          lprcoords = .false. ! ParRep core has official coords
          lprwarm = .true.    ! ParRep core needs to warm up
          lprblack = .false.  ! ParRep core needs to thermalize/blackout
          transent = .true.   ! Last transition message was completed




          lprneeded = .true.  !!! Set to .false.
                               !!!  .. and write "masters" code to have
                               !!!  .. adaptive "nParRep"
          !! Start with 1 replica allowed to accumulate MD:
          !!   Others will be invited by MASTER
          !if(rank_l.lt.(nforcecores*2)) lprneeded = .true.




!          if(rank_f.eq.0) write(*,*) 'spawnid,rank_l '
!     x          ,spawnid,rank_l,' - entering parRep block'

          ! Allocate array for holding blocks of MD time
          if(.not.allocated(prmdtarray))then
            allocate(prmdtarray(100000000))
          else
            prmdtarray(:) = 0.d0
          endif
          ! Allocate array for holding blocks of Wall Clock time
          if(.not.allocated(prwctarray))then
            allocate(prwctarray(100000000))
          else
            prwctarray(:) = 0.d0
          endif
          hightime = 0.d0
          t_rep_tot_i = t_rep_tot

          do while (lprloop)

            if(activated.and.(.not.lfinished))then
              ! Check if we are acpr_send_coordstive
              call pr_check_deactivate(transent)
              call MPI_BARRIER(force_comm,ier)
            endif
              if(lfinished.and.transent)then
                if(rank_f.eq.0)
     x            write(6,*) ' GROUPID ',groupID,' ParRep Core EXITING.'
                if(rank_f.eq.0) call flush(6)
                if(rank_f.eq.0) call MPI_Send(.true.,1,MPI_LOGICAL
     +                                 ,0,deac_tag2,local_comm,ier)
                call MPI_BARRIER(force_comm,ier)
                goto 600
              endif
              if((.not.activated).and.transent)then
                if(rank_f.eq.0) call MPI_Send(.true.,1,MPI_LOGICAL
     +                                 ,0,deac_tag2,local_comm,ier)
                call MPI_BARRIER(force_comm,ier)
                goto 601
              endif

            ! Get Initial Coords
            if(.not.lprcoords) then
              call pr_get_coords(natom,xyz,lprcoords)
              if(lprcoords) then
                tadcharge(:) = 0.d0
                call vecmov(xyz,xyz1,3*natom)
                xyzneighs(:)=0
                call vecmov(xyz1,xyzneighs(1),3*natom)
                nneigh = 0
                iblock = 0
                hightime = 0.d0

                if((spawnid.gt.1).and.
     +             (ietype.ge.100).and.(ietype.lt.200))then
                  ! Startup fresh lammps run
                  call MPI_BARRIER(force_comm,ier)
                  call lmp_setup(natom,nmove,ntype,maxtyp,potnam
     +              ,ietype,xyz1,taxes,itype,rcut,amu
     +              ,templow,temphigh,.true.,dt)
                  call MPI_BARRIER(force_comm,ier)
                  call vecmov(xyz1,xyz,3*natom)
                endif

              endif
            endif

            ! ParRep Steps
            ! ... Wait to have official coords before doing any MD
            if(lprcoords)then

              if(activated.and.(.not.lfinished).and.(iblock.eq.0))then
                iblock = iblock + 1
                if(rank_f.eq.0) checktime = mpi_wtime()-cput_ref
                call MPI_BCAST(checktime,1,MPI_REAL8,0,force_comm,ier)
                prwctarray(iblock) = checktime-parRepRefTime
                prmdtarray(iblock) = hightime
!                write(*,*) '0 rank_l= ',rank_l,' - iblock= ',
!     +                     iblock,' - prwct, prmdt: ',
!     +                     prwctarray(iblock), prmdtarray(iblock)
              endif

              ! Check if master needs to synchronize for transition
              call pr_sync_trans(iblock,dt,transent,lprneeded)
              ! Check if master needs to synchronize for time
              call pr_share_time(iblock,transent,lprneeded)

              ! Warmup/Blackout (if necessary)
              if(lprwarm.or.lprblack) then

               call MPI_BARRIER(force_comm,ier)
!              if(rank_f.eq.0) write(*,*) 'spawnid,rank_l '
!     x                    ,spawnid,rank_l,' - entering warmup/blackout'
!              call flush(6)

               ineedwarm = 0
               if(lprwarm)then
                 ineedwarm = nwarm+nthermalize
                 call gaussp(natom,nmove,pxyz,itype,maxtyp,amu,ntype,
     x                moving,mx,my,mz,2.0*temphigh,0,-1)
                 call MPI_BCAST(pxyz,3*nmove,MPI_REAL8,0,force_comm,ier)
                 call mdreset(0)
               elseif(lprblack)then
                 ineedwarm = nblack
               endif
               ibadtrans=0
 101           continue                  ! redo warming loop
               call vecmov(xyz,xyz3,3*natom)
               if((ietype.ge.100).and.(ietype.lt.200).and.doblocks)then
                 call mdblock100(natom,nmove,xyz,pxyz,itype,
     x              maxtyp,e,grad,taxes,ietype,amu,temphigh,dt,
     x              thermrate,ineedwarm,iseed,ntype)
!                 if(rank_f.eq.0) write(*,*) 'spawnid,rank_l '
!     x                    ,spawnid,rank_l,' - did lmp warmup block '
                 call flush(6)
               else
                 blocking_on = .true.
                 do i=1,ineedwarm
                   call mdstep(natom,nmove,xyz,itype,pxyz, taxes,
     x                ietype, maxtyp,rcut,amu,
     x                maxtherm,itherm,jtherm,ttherm,ctherm,
     x                xtherm,mtherm,ktherm,ntherm,
     x                moving,mx,my,mz,dt,time,thyper,hyperratio,e,grad)
                 enddo
                 blocking_on = .false.
               endif
               !! RJZ - shiftTest
               if (natom.eq.nmove) then
                  call shiftalign(xyz3,xyz,natom)
               endif
               call vecmov(xyz,xyz2,3*natom)
               call MPI_BARRIER(force_comm,ier)

               ! Check if master needs to synchronize for transition
               call pr_sync_trans(iblock,dt,transent,lprneeded)
               ! Check if master needs to synchronize for time
               call pr_share_time(iblock,transent,lprneeded)

!              write(*,*) 'spawnid,rank_l '
!     x        ,spawnid,rank_l,' - doing warmup/blackout descent_check()'
!              call flushtad(6)
               imode=0 ! do not do spring quench
               imatch=0
               call descent_check(natom,nmove,1+nneigh,xyz2,xyzneighs
     +              ,itype,taxes,ietype,maxtyp,rcut,hyperratio
     +              ,barrevev,ereverse,imatch,irotalign,nsteep2
     +              ,gfac,transcrit,itranscrit,gcrit,dvcrit,drcrit
     +              ,ipr,imode,ifundc)
               call MPI_BARRIER(force_comm,ier)
!               write(*,*) 'spawnid,rank_l '
!     x     ,spawnid,rank_l,' - finished warmup/blackout descent_check()'
!     x     ,' - imatch = ',imatch
!               call flushtad(6)
                if ((ereverse.gt.0.0d0).and.(imatch.gt.1)) then
                   if (barrevev(imatch-1).lt.ereverse) then
                      imatch=1
                   endif
                endif
!              write(*,*) 'spawnid,rank_l '
!     x     ,spawnid,rank_l,' - descent_check() imatch = ',imatch
!              call MPI_BCAST(imatch,1,MPI_INTEGER,0,force_comm,ier)
!              call flushtad(6)

                if(imatch.ne.1) then
                   ibadtrans=ibadtrans+1
                   if(ibadtrans.lt.(maxwarmreject-1)) then
                     call vecmov(xyz3,xyz,3*natom)
                     call gaussp(natom,nmove,pxyz,itype,maxtyp,amu
     x                    ,ntype,moving,mx,my,mz,temphigh,iop,iwrt)
                     call MPI_BCAST
     x                    (pxyz,3*nmove,MPI_REAL8,0,force_comm,ier)
                     call mdreset(0)
                     if(rank_f.eq.0) write(*,*) 'spawnid,rank_l '
     x                    ,spawnid,rank_l,' - redoing warmup/blackout'
                     if(rank_f.eq.0) call flushtad(6)
                     ! Transition durring warmup, REDO
                     goto 101
                   elseif(ibadtrans.lt.maxwarmreject) then
                     call vecmov(xyz1,xyz,3*natom)
                     ineedwarm = nwarm+nthermalize
                     call gaussp(natom,nmove,pxyz,itype,maxtyp,amu
     x                    ,ntype,moving,mx,my,mz,2.0*temphigh,iop,iwrt)
                     call MPI_BCAST
     x                    (pxyz,3*nmove,MPI_REAL8,0,force_comm,ier)
                     call mdreset(0)
                     if(rank_f.eq.0) write(*,*) 'spawnid,rank_l '
     x                    ,spawnid,rank_l,' - redoing warmup/blackout'
                     if(rank_f.eq.0) call flushtad(6)
                     ! Transition durring warmup, REDO
                     goto 101
                   else
                     ! Too many warmup/blackout failures... Ignore (for now)
                     call vecmov(xyz3,xyz,3*natom)
                     if(rank_f.eq.0) write(*,*) 'spawnid,rank_l '
     x                    ,spawnid,rank_l
     x                    ,' - NOT bothering to redo warmup/blackout'
                     if(rank_f.eq.0) call flushtad(6)
                   end if
                else
                   call MPI_BARRIER(force_comm,ier)
!                   if(rank_f.eq.0) write(*,*) 'spawnid,rank_l '
!     x                    ,spawnid,rank_l,' - finished warmup/blackout'
!                   if(rank_f.eq.0) call flushtad(6)
                endif
              endif

              timenowi = mpi_wtime()-cput_ref
 102          continue
              ! Check if master needs to synchronize for transition
              call pr_sync_trans(iblock,dt,transent,lprneeded)
              ! Check if master needs to synchronize for time
              call pr_share_time(iblock,transent,lprneeded)
              ! must loop here until previous transition messages are recieved
              call pr_test_trans(transent)
              if(.not.transent) then
                timenow = mpi_wtime()-cput_ref
                if((timenow-timenowi).gt.(10*60.0))then
                  if(rank_f.eq.0) write(*,*)
     +                'SpawnID ',SpawnID,' rank_l ',rank_l
     +                ,' still waiting on pr_test_trans() -_-'
                  if(rank_f.eq.0) call flush(6)
                  timenowi = mpi_wtime()-cput_ref
                endif
                call MPI_BARRIER(force_comm,ier)
              endif
              if(.not.transent) goto 102

              call MPI_BARRIER(force_comm,ier)
              if(activated.and.(.not.lfinished))then
                ! Check if we are active
                call pr_check_deactivate(transent)
                call MPI_BARRIER(force_comm,ier)
              endif
              if(lfinished.and.transent)then
                if(rank_f.eq.0)
     x            write(6,*) ' GROUPID ',groupID,' ParRep Core EXITING.'
                if(rank_f.eq.0) call flush(6)
                if(rank_f.eq.0) call MPI_Send(.true.,1,MPI_LOGICAL
     +                                 ,0,deac_tag2,local_comm,ier)
                call MPI_BARRIER(force_comm,ier)
                goto 600
              endif
              if((.not.activated).and.transent)then
                if(rank_f.eq.0) call MPI_Send(.true.,1,MPI_LOGICAL
     +                                 ,0,deac_tag2,local_comm,ier)
                call MPI_BARRIER(force_comm,ier)
                goto 601
              endif

              if(.not.lprneeded)then
                if(iblock.eq.1)then
                  iblock = iblock + 1
                  if(rank_f.eq.0) checktime = mpi_wtime()-cput_ref
                  call MPI_BCAST(checktime,1,MPI_REAL8,0,force_comm,ier)
                  prwctarray(iblock) = checktime-parRepRefTime
                  prmdtarray(iblock) = hightime
!                  write(*,*) 'N rank_l= ',rank_l,' - iblock= ',
!     +                     iblock,' - prwct, prmdt: ',
!     +                     prwctarray(iblock), prmdtarray(iblock)
                endif
                goto 102 ! Loop back to 102 if not needed for MD yet
              endif

              if(activated.and.(.not.lfinished))then

                ! Run MD block
                call vecmov(xyz,xyz3,3*natom)
                ihold=0
                iblock = iblock + 1
                if(rank_f.eq.0) checktime = mpi_wtime()-cput_ref
                call MPI_BCAST(checktime,1,MPI_REAL8,0,force_comm,ier)
                prwctarray(iblock) = checktime-parRepRefTime
                prmdtarray(iblock) = hightime
!                write(*,*) 'A rank_l= ',rank_l,' - iblock= ',
!     +                     iblock,' - prwct, prmdt: ',
!     +                     prwctarray(iblock), prmdtarray(iblock)

                ! Do nhold loops if using lammps mdblock routine:
                nmd_use = nmd
                if((ietype.ge.100).and.(ietype.lt.200)
     x                                              .and.doblocks)then
                    nmd_use = nhold !nmd / nmdsub
                endif
                do i=1,nmd_use
                  if((ietype.ge.100).and.(ietype.lt.200)
     x                                              .and.doblocks)then
                     call mdblock100(natom,nmove,xyz,pxyz,itype,
     x                    maxtyp,e,grad,taxes,ietype,amu,temphigh,dt,
     x                    thermrate,nmdsub,iseed,ntype)
                     ngradcall = ngradcall + nmdsub
                     time = time + dt * nmdsub
                  else
                     !blocking_on = .true.
                     call mdstep(natom,nmove,xyz,itype,pxyz,
     x                 taxes,ietype, maxtyp,rcut,amu,
     x                 maxtherm,itherm,jtherm,ttherm,ctherm,
     x                 xtherm,mtherm,ktherm,ntherm,moving,
     x                 mx,my,mz,dt,time,thyper,hyperratio,e,grad)
                     !blocking_on = .false.
                  endif
                  ek=0.0d0
                  do ibpu=1,3*nmove
                    it=itype(int((ibpu-1)/3)+1)
                    ek=ek+(pxyz(ibpu)**2)/(amu(it)*1822.83d0)/2.0d0
                  enddo
                  tempnow=ek*1.0d0/(3.0d0/2.0d0*3.167d-6)/float(nmove)
!                  if((mod(i,nmdsub).eq.0).or.
!     +                 ((ietype.ge.100).and.(ietype.lt.200))) then
                  if(mod(i,nmdsub).eq.0) then
                    ihold=ihold+1
                    !! RJZ - shiftTest
                    if (natom.eq.nmove) then
                      call shiftalign(xyz3,xyz,natom)
                    endif
                    call vecmov
     x                  (xyz,xyzhold((ihold-1)*3*natom+1),3*natom)
                  endif
                enddo
                !! RJZ - shiftTest
                if (natom.eq.nmove) then
                  call shiftalign(xyz3,xyz,natom)
                endif

                ! Quickly check if master needs some time...
                call pr_share_time(iblock,transent,lprneeded)

                if(rank_f.eq.0) checktime = mpi_wtime()-cput_ref
                call MPI_BCAST(checktime,1,MPI_REAL8,0,force_comm,ier)
                hightime = hightime + nmd*dt
                iblock = iblock + 1
                prwctarray(iblock) = checktime-parRepRefTime
                prmdtarray(iblock) = hightime
!                write(*,*) 'B rank_l= ',rank_l,' - iblock= ',
!     +                     iblock,' - prwct, prmdt: ',
!     +                     prwctarray(iblock), prmdtarray(iblock)
                call vecmov(xyz,xyz2,3*natom)
                imode=0 ! do not do spring quench
                ijump=0
                call descent_check(natom,nmove,1+nneigh,xyz2,xyzneighs
     +             ,itype,taxes,ietype,maxtyp,rcut,hyperratio,barrevev
     +             ,ereverse,isame,irotalign,nsteep2,gfac,transcrit
     +             ,itranscrit,gcrit,dvcrit,drcrit,ipr,imode,ifundc)
                if(isame.ne.1) ijump=1
              else
                ijump = 0
              endif

              rat_upper = 1.0
              rat_lower = 0.0
              d_rat     = 1.0

              ! Now refine the transition:
              if((ijump.eq.1).and.(nhold.gt.1)) then

                !  refine the transition point, finding the first out-of-basin geometry
                imode=0 ! do not do spring quench
                call vecmov(xyz2,xyz5,3*natom)
                call refine_transition(natom,nmove,ihold,xyz1,xyz5
     +              ,xyzhold,xyzneighs,1+nneigh,barrevev,ereverse
     +              ,itype,taxes,ietype,maxtyp,rcut,hyperratio
     +              ,irotalign,nsteep1,gfac,transcrit,itranscrit
     +              ,intermediate,gcrit,dvcrit,drcrit,ipr,imode
     +              ,ileft,iright,ifundc,rank_f)
                if(ileft.eq.1) ileft = 0
                rat_upper = real(iright) / real(ihold)
                rat_lower = real(ileft)  / real(ihold)
                d_rat     = rat_upper - rat_lower
              endif


              ! Send dtime and transition result
              if(ijump.eq.1)then
                ijump = 0
!                if(rank_f.eq.0)
!     +            write(6,*) 'spawnID ',spawnID,' rank_l ',rank_l
!     +             ,' rank_w ',rank_w
!     +             ,' has ParRep Transition'
!                if(rank_f.eq.0) call flushtad(6)
                call pr_send_trans(nmove,natom,ihold,xyz2,xyzhold
     +              ,nholdmaxworstcase,iblock,nmd,dt,transent,hightime
     +              ,rat_lower,d_rat)

                n_rep_trans = n_rep_trans + 1
                t_rep_tot = t_rep_tot_i + hightime
                transtimerep = t_rep_tot / n_rep_trans
!                if(rank_f.eq.0)
!     +            write(6,*) 'spawnID ',spawnID,' rank_l ',rank_l
!     +             ,' has AVG trans time:',transtimerep
!     +                   ,'- n_rep_trans:',n_rep_trans
!     +                   ,'- t_rep_tot:',t_rep_tot
!                if(rank_f.eq.0) call flushtad(6)

                ! Prepare for next MD loop
                call vecmov(xyz3,xyz,3*natom)
                call gaussp(natom,nmove,pxyz,itype,maxtyp,amu,ntype
     x                               ,moving,mx,my,mz,temphigh,0,-1)
                call MPI_BCAST
     x                    (pxyz,3*nmove,MPI_REAL8,0,force_comm,ier)
                call mdreset(0)
                lprblack = .true.
                lprwarm  = .false.
              else
!                if(rank_f.eq.0)
!     +            write(6,*) 'spawnID ',spawnID,' rank_l ',rank_l
!     +             ,' has no ParRep Transition..'
!                if(rank_f.eq.0) call flushtad(6)
                lprblack = .false.
                lprwarm  = .false.

                t_rep_tot = t_rep_tot_i + hightime

              endif

!              write(6,*) 'spawnID ',spawnID,' rank_l ',rank_l
!     +             ,' has transent = ',transent

              ! Check if master needs to synchronize for transition
              call pr_sync_trans(iblock,dt,transent,lprneeded)
              ! Check if master needs to synchronize for time
              call pr_share_time(iblock,transent,lprneeded)

              if(activated.and.(.not.lfinished))then
                ! Check if we are active
                call pr_check_deactivate(transent)
                call MPI_BARRIER(force_comm,ier)
!                if(rank_f.eq.0) write(*,*) 'spawnid,rank_l '
!     x                    ,spawnid,rank_l,' - activated,lfinished C = '
!     x                    ,activated,lfinished
!                if(rank_f.eq.0) call flush(6)
              endif
              if(lfinished.and.transent)then
                if(rank_f.eq.0)
     +            write(6,*) ' GROUPID ',groupID,' ParRep Core EXITING.'
                if(rank_f.eq.0) call flush(6)
                if(rank_f.eq.0) call MPI_Send(.true.,1,MPI_LOGICAL
     +                                 ,0,deac_tag2,local_comm,ier)
                call MPI_BARRIER(force_comm,ier)
                goto 600
              endif
              if((.not.activated).and.transent)then
!                if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
!     +           ,' rank_l ',rank_l
!     +           ,' - Sending deac_tag2 message.'
!                if(rank_f.eq.0) call flush(6)
                if(rank_f.eq.0) call MPI_Send(.true.,1,MPI_LOGICAL
     +                                 ,0,deac_tag2,local_comm,ier)
                call MPI_BARRIER(force_comm,ier)
                goto 601
              endif

            endif

          enddo
          goto 601
        endif
CCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCcc
CC
CC     ...End ParRep Core Loop
CC
CCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCccCCcc

        ! Startup fresh lammps run
        if((spawnid.gt.1).and.(ietype.ge.100).and.(ietype.lt.200))then
          call MPI_BARRIER(force_comm,ier)
          call lmp_setup(natom,nmove,ntype,maxtyp,potnam
     +        ,ietype,xyz,taxes,itype,rcut,amu
     +        ,templow,temphigh,.true.,dt)
          call MPI_BARRIER(force_comm,ier)
        endif

c If ldepspawn... need to run a bit of MD:
        if(ldepspawn.eq.1)then

          if(ietype.eq.103)then
            lmpcmd =
     x        'fix 1 all qeq/slater 1 12.0 1.0e-7 10000 coul/streitz'
            call lammps_command (lmp,lmpcmd)
          endif

          iop=0
          iwrt=-1
          call gaussp(natom,nmove,pxyz,itype,maxtyp,amu,ntype,
     x                moving,mx,my,mz,templow,iop,iwrt)
          call MPI_BCAST(pxyz,3*nmove,MPI_REAL8,0,force_comm,ier)
          call mdreset(1)

          xyz(3) = xyz(3) + 4.d0 !Give deposited atom a bit more room
          timehold = time
          ii = 0
          call vecmov(xyz,xyz1,3*natom)
 215  continue
          call vecmov(xyz1,xyz,3*natom)

          !tempdep = 1356 ! Cu metling temperature
          !tempdep = 0.5*tempdep ! Half 0.2eV -> 0.1 eV ? As in Sprague?
          !ttherm(1)=templow
          !pxyz(1)=0.d0
          !pxyz(2)=0.d0
          !pxyz(3)=-sqrt(4.d0*amu(1)*1822.83d0*tempdep*3.167d-6)

          ! EAM Growth...
          if((etype.eq.0).or.(etype.eq.100))then
            pxyz(1)=0.d0
            pxyz(2)=0.d0
            massdep = amu(itype(1))*1822.83d0
            engdep = 0.1 ! DEPOSITION ENERGY
            pxyz(3)=-sqrt(2.0*massdep*engdep)
          endif

          ! Ceramic Surface Growth...
          if((etype.eq.172).or.(etype.eq.72))then
            pxyz(1)=0.d0
            pxyz(2)=0.d0
            massdep = amu(itype(1))*1822.83d0
            if(itype(1).eq.1)then
                engdep = 0.025
            elseif(itype(1).eq.2)then
                engdep = 0.025
            else
                engdep = 0.025
            endif
            pxyz(3)=-sqrt(2.0*massdep*engdep)
          endif

          if(rank_f.eq.0) xyz(1) = prngen(0)*taxes(1)
          if(rank_f.eq.0) xyz(2) = prngen(0)*taxes(2)
          call MPI_BCAST(xyz(1),2,MPI_REAL8,0,force_comm,ier)
          do i=1,ndepMDsteps
            if(((MOD(i,1000).eq.0).or.(i.eq.1)).and.(ipr.ge.4))then
              ii = ii + 1
              filn=''
              write(filn,*) spawnID
              filnb=''
              write(filnb,*) ii
              filn='sp.'//trim(adjustl(filn))//
     +               '.'//trim(adjustl(filnb))//'.dat'
              call storefile(natom,xyz,pxyz,itype,taxes,filn)
            endif
            !blocking_on = .true.
            call mdstep(natom,nmove,xyz,itype,pxyz,taxes,
     +             ietype,maxtyp,rcut,amu,
     +             maxtherm,itherm,jtherm,ttherm,ctherm,
     +             xtherm,mtherm,ktherm,ntherm,
     +             moving,mx,my,mz,dt,time,thyper,hyperratio,e,grad) ! xxx time
            !blocking_on = .false.
            !write(*,*) 'spawnID ',spawnID,' did deposition step ',i
          enddo

          !! RJZ - shiftTest
          if (natom.eq.nmove) then
              call shiftalign(xyz1,xyz,natom)
          endif

          ! Move stray atoms into PBCs
          do jj = 1, nmove
            if(xyz((jj-1)*3+1).gt.taxes(1))
     +        xyz((jj-1)*3+1) = xyz((jj-1)*3+1) - taxes(1)

            if(xyz((jj-1)*3+2).gt.taxes(2))
     +        xyz((jj-1)*3+2) = xyz((jj-1)*3+2) - taxes(2)

            if(xyz((jj-1)*3+1).lt.0.d0)
     +        xyz((jj-1)*3+1) = xyz((jj-1)*3+1) + taxes(1)

            if(xyz((jj-1)*3+2).lt.0.d0)
     +        xyz((jj-1)*3+2) = xyz((jj-1)*3+2) + taxes(2)
          enddo
          do idum = 1, nmove
              distmin = 1.d30
              do jj = 1, nmove
                if(idum.ne.jj)then
                  distx = xyz((jj-1)*3+1) - xyz((idum-1)*3+1)
                  disty = xyz((jj-1)*3+2) - xyz((idum-1)*3+2)
                  distz = xyz((jj-1)*3+3) - xyz((idum-1)*3+3)
                  if(sqrt(distx**2+disty**2+distz**2).lt.distmin)then
                    distmin = sqrt(distx**2+disty**2+distz**2)
                  endif
                endif
              enddo
              if(distmin.gt.(3.d0))then
               if(rank_f.eq.0)
     x          write(*,*) 'SpawnID ',spawnID,' - moving atom ',distmin
     +             ,' from nearest neighbor)'
               goto 215
              endif
          enddo
          tadtime = tadtime + dt*ndepMDsteps ! we have directly accumulated more tadtime now
          time = timehold
          ntherm=1
          itherm(1)=1
          jtherm(1)=nmove
          ctherm(1)=thermrate
          ttherm(1)=temphigh
          mtherm(1)=0
          ktherm(1)=4
          ineedgauss=1
          call mdreset(1)

          if(ietype.eq.103)then
            call lmp_minimize(natom,nmove,xyz,itype,ietype
     x        ,energy,taxes)
            lmpcmd = 'unfix 1'
            call lammps_command (lmp,lmpcmd)
          endif

        else

          if((spawnid.gt.1).and.(ietype.eq.103))then
             lmpcmd =
     x        'fix 1 all qeq/slater 1 12.0 1.0e-7 10000 coul/streitz'
            call lammps_command (lmp,lmpcmd)
            call lmp_minimize(natom,nmove,xyz,itype,ietype
     x        ,energy,taxes)
            lmpcmd = 'unfix 1'
            call lammps_command (lmp,lmpcmd)
          endif

        endif

        if(.FALSE.)then !! Run MD for debugging
          nneigh = 0
          xyzneighs(:) = 0
          call vecmov(xyz,xyzneighs(1),3*natom)
          ii = 0
          call vecmov(xyz,xyz3,3*natom)
          do i=1,10000000
            if((MOD(i,10000).eq.0).or.(i.eq.-1))then
              call flush(6)
              call vecmov(xyz,xyz2,3*natom)
              imode=0 ! do not do spring quench
              !call descent_check(natom,nmove,nneigh+1,xyz2,xyzneighs
              call descent_check(natom,nmove,0,xyz2,xyzneighs
     +          ,itype,taxes,ietype,maxtyp,rcut,hyperratio,barrevev
     +          ,ereverse,imatch,irotalign,nsteep2,gfac,transcrit
     +          ,itranscrit,gcrit,dvcrit,drcrit,ipr,imode,ifundc)
              !if(imatch.ne.1)then
              if(.false.)then
                iop=0
                iwrt=-1
                call gaussp(natom,nmove,pxyz,itype,maxtyp,amu,ntype,
     x                moving,mx,my,mz,temphigh,iop,iwrt)
                call MPI_BCAST(pxyz,3*nmove,MPI_REAL8,0,force_comm,ier)
                call mdreset(0)
                if(imatch.eq.0)then
                  nneigh = nneigh + 1
                  call vecmov(xyz2,xyzneighs(3*natom*nneigh+1),3*natom)
                endif
              endif
              if((i.eq.1).or.(imatch.eq.0).or.(.true.))then
                call flush(6)
                ii = ii + 1
                filn=''
                write(filn,*) spawnID
                filnb=''
                write(filnb,*) ii
                filn='mdtest.'//trim(adjustl(filn))//
     +               '.'//trim(adjustl(filnb))//'.dat'
                call storefile(natom,xyz2,pxyz,itype,taxes,filn)
                call flush(6)
                call MPI_BARRIER(force_comm,ier)
              endif
            endif
            call mdstep(natom,nmove,xyz,itype,pxyz,taxes,
     +             ietype,maxtyp,rcut,amu,
     +             maxtherm,itherm,jtherm,ttherm,ctherm,
     +             xtherm,mtherm,ktherm,ntherm,
     +             moving,mx,my,mz,dt,time,thyper,hyperratio,e,grad) ! xxx time
            !if((MOD(i,1000).eq.0).or.(i.eq.1))then
            if(.false.)then
                !! RJZ - shiftTest
                if (natom.eq.nmove) then
                  call shiftalign(xyz3,xyz,natom)
                endif
                call flush(6)
                ii = ii + 1
                filn=''
                write(filn,*) spawnID
                filnb=''
                write(filnb,*) ii
                filn='mdtest.'//trim(adjustl(filn))//
     +               '.'//trim(adjustl(filnb))//'.dat'
                call storefile(natom,xyz,pxyz,itype,taxes,filn)
                call flush(6)
                call MPI_BARRIER(force_comm,ier)
            endif
          enddo
          STOP ' THIS WAS A FAKE MD-DEBUG RUN!! '
        elseif(.false.)then ! Simple Hack to swap axes...
          call vecmov(xyz,xyz1,3*natom)
          xtax__tmp = taxes(1)
          ytax__tmp = taxes(2)
          ztax__tmp = taxes(3)
          taxes(1) = xtax__tmp
          taxes(2) = ztax__tmp
          taxes(3) = ytax__tmp
          do i=1,natom
            x__tmp = xyz1((i-1)*3+1)
            y__tmp = xyz1((i-1)*3+2)
            z__tmp = xyz1((i-1)*3+3)
            xyz1((i-1)*3+1) = x__tmp
            xyz1((i-1)*3+2) = z__tmp
            xyz1((i-1)*3+3) = y__tmp
            pxyz((i-1)*3+1) = 0.d0
            pxyz((i-1)*3+2) = 0.d0
            pxyz((i-1)*3+3) = 0.d0
          enddo
          filn=''
          write(filn,*) spawnID
          filn='rotate.'//trim(adjustl(filn))//'.dat'
          call storefile(natom,xyz1,pxyz,itype,taxes,filn)
          call MPI_BARRIER(force_comm,ier)
          STOP ' THIS WAS A HACK AXIS-SWAP RUN!! '
        endif

c move config into xyz1, which will become the minimum
        call vecmov(xyz,xyz1,3*natom)
        call vecmov(xyz,xyz4,3*natom) ! save in case we don't make it through warming first time

        ngradcallhold=ngradcall
        imode=0
        if (itranscrit.eq.1.and.itrans.gt.1.and.irotalign.ge.1)
     +      imode=1 ! spring quench mode
        if (itranscrit.eq.1.and.itrans.gt.1) iflag=1

        if(irecognize.eq.0) then
           iflag=0
           ifundcuse=ifundc
           if (ifundcuse.eq.5) ifundcuse=0
           ereverseuse=0.0d0    ! xyzprev has nothing to do with reverse barrier stuff, don't use it
             call descent_check(natom,nmove,iflag,xyz1,xyzprev
     +         ,itype,taxes,ietype,maxtyp,rcut,hyperratio,barrevev
     +         ,ereverseuse,iofficialmatch,irotalign,nsteep1,gfac
     +         ,transcrit,itranscrit,gcrit,dvcrit,drcrit,ipr,0,ifundcuse
     +         )
             if (imode.ne.0) then
              if (itranscrit.eq.1.and.itrans.gt.1) iflag=1
              call descent_check(natom,nmove,iflag,xyz1,xyzprev
     +            ,itype,taxes,ietype,maxtyp,rcut,hyperratio,barrevev
     +            ,ereverseuse,iofficialmatch,irotalign,nsteep1,gfac
     +            ,transcrit,itranscrit,gcrit,dvcrit,drcrit,ipr,imode
     +            ,ifundc)
             endif
             call gcalc(natom,nmove,xyz1,itype, taxes,
     x         ietype, maxtyp,rcut, energy,grad, hyperratio)

           emin=energy
           if(ipr.ge.1) then
              if(rank_f.eq.0)
     x          write(6,141) istate,emin*27.21d0
 143          format(' minimum no.',i5,' stored; energy =',f16.4
     +            ,' eV')
           end if
           iofficial=0
           ilabel=istate    ! xxx created ilabel to label output files since iofficial not defined for irecognize=0
c xx perhaps rename xyz1 to xyzmin

c ==== if irecognize=1, look for match with an official state while minimizing ====
        else if(irecognize.ge.1) then
c           imode=0 ! do not do spring quench if recognizing states
           ifundcuse=ifundc
           if (ifundcuse.eq.5) ifundcuse=0
           ereverseuse=0.0d0

          ! First try recognizing without any force calls:
          !iofficialmatch = 0
          drmax_min = 1.d99
          if(iofficialmatch.eq.0)then ! May already know what state you are in...
            if((parallelmode.eq.'spectad').and.(ldepspawn.ne.1))then
               !drmax_min = 1.d99
               if(nofficial.ge.2)then
               do j=2,nofficial
                if(nmove.eq.statedata(j)%nmove)then
                  call maxdr_rz(nmove,xyz1
     +              ,statedata(j)%xyzmin,deltarmax)
                  if (deltarmax.lt.transcrit) then ! potential match
                      if(deltarmax.lt.drmax_min)then
                        drmax_min = deltarmax
                        iofficialmatch = j
                        energy = statedata(j)%emin
                      endif
                  endif
                endif
               enddo
               endif
            endif
            if(iofficialmatch.ne.0)then
!               if(rank_f.eq.0)then
!                  write(*,*) 'WARNING -> SpawnID',SpawnID,
!     x             'decided that it is in istate=',iofficialmatch,
!     x             '- with drmax_min=',drmax_min
!                  call flush(6)
!               endif
!               call MPI_BARRIER(force_comm,ier)
               if(drmax_min.gt.MAX((drmax_min/2.0),0.1))then
                 if(rank_f.eq.0)then
                   write(*,*) 'ERROR -> SpawnID',SpawnID,
     x               'matched state with drmax_min=',drmax_min,
     x               ' - THIS IS PROBABLY WRONG !!!!'
                   call flush(6)
                 endif
               endif
            endif
          endif
          if((ipr.ge.4).and.(rank_f.eq.0))
     x      write(6,*) 'SPECTAD: irecognize=',irecognize,
     x                  '  iofficialmatch=',iofficialmatch,
     x                  '  nofficial=',nofficial,
     x                  '  drmax_min=',drmax_min

          if(parallelmode.eq.'spectad')then
            ! No match... now use force calls (make sure we are minimized):
            if(iofficialmatch.eq.0)then !.and.(.false.))then
              call descent_check(natom,nmove,0,xyz1,xyzofficial
     +         ,itype,taxes,ietype,maxtyp,rcut,hyperratio,barrevev
     +         ,ereverseuse,iofficialmatch,irotalign,nsteep1,gfac
     +         ,transcrit,itranscrit,gcrit,dvcrit,drcrit,ipr,0,ifundcuse
     +         )
            endif
          else
            call descent_check(natom,nmove,nofficial,xyz1,xyzofficial
     +         ,itype,taxes,ietype,maxtyp,rcut,hyperratio,barrevev
     +         ,ereverseuse,iofficialmatch,irotalign,nsteep1,gfac
     +         ,transcrit,itranscrit,gcrit,dvcrit,drcrit,ipr,0,ifundcuse
     +         )
            if (imode.ne.0) then
             itemp=1
             call descent_check(natom,nmove,itemp,xyz1,xyzofficial
     +           ,itype,taxes,ietype,maxtyp,rcut,hyperratio,barrevev
     +           ,ereverseuse,idum,irotalign,nsteep1,gfac,transcrit
     +           ,itranscrit,gcrit,dvcrit,drcrit,ipr,imode,ifundc)
            endif
          endif

          if((ipr.ge.4).and.(rank_f.eq.0))
     x      write(6,*) 'irecognize=',irecognize,
     x                  '  iofficialmatch=',iofficialmatch

          if (iofficialmatch.eq.0) then
             call gcalc(natom,nmove,xyz1,itype, taxes,
     x           ietype, maxtyp,rcut, energy,grad, hyperratio)
             emin=energy
          endif

      emin0 = emin

c ====energy is now minimized, so try simple energy comparison if irec..=9     9/1/04
c   Note: the indentation is goofed up in this general area,
c   but I will leave it alone for now so that it is easier to tell what
c   changes I made to put in the irecognize=9 option.  xxxafv 9/1/04
      if((iofficialmatch.eq.0).and.(irecognize.eq.9)
     x                        .and.(parallelmode.eq.'spectad')) then
        do i=2,nofficial
          if(iofficialmatch.eq.0) then
            de=abs(statedata(i)%emin-emin)
            if((de.lt.demin).or.(i.eq.2)) demin=de
            if(de.le.erecognizecrit) then    ! hartree
              callpairsum(natom,xyz1,r2sum)
              callpairsum(natom,statedata(i)%xyzmin,r2sumo)
              if(abs(r2sum-r2sumo).le.1.0d0) then   ! using 1.0 for the criterion for now
                iofficialmatch=i
                if((ipr.ge.4).and.(rank_l.eq.0))
     x            write(6,*)'irecognize=9 matchup:',istate,i
              else
                if(rank_l.eq.0) write(6,*)
     x           'energy match rejected based on pairsum',istate,i,
     x           r2sum,r2sumo
              end if
            end if
          endif
         end do
         if((ipr.ge.3).and.(rank_l.eq.0)) write(6,*)
     x         'irecognize=9: istate=',istate,'  iofficialmatch =',
     x         iofficialmatch,'  demin =',demin
         if((ipr.ge.5).and.(rank_l.eq.0))
     x         write(6,*) emin,(eofficial(ix),ix=1,nofficial)
      elseif((iofficialmatch.eq.0).and.(irecognize.eq.9)) then
        do i=1,nofficial
          if(iofficialmatch.eq.0) then
            de=abs(eofficial(i)-emin)
            if(de.lt.demin .or. i.eq.1) demin=de
            if(de.le.erecognizecrit) then    ! hartree
              callpairsum(natom,xyz1,r2sum)
              callpairsum(natom,xyzofficial((i-1)*3*natom+1),r2sumo)
              if(abs(r2sum-r2sumo).le.1.0d0) then   ! using 1.0 for the criterion for now
                iofficialmatch=i
                if((ipr.ge.4).and.(rank_l.eq.0))
     x            write(6,*)'irecognize=9 matchup:',istate,i
              else
                if(rank_l.eq.0) write(6,*)
     x           'energy match rejected based on pairsum',istate,i,
     x           r2sum,r2sumo
              end if
            end if
          endif
         end do
         if((ipr.ge.3).and.(rank_l.eq.0)) write(6,*)
     x         'irecognize=9: istate=',istate,'  iofficialmatch =',
     x         iofficialmatch,'  demin =',demin
         if((ipr.ge.5).and.(rank_l.eq.0))
     x         write(6,*) emin,(eofficial(ix),ix=1,nofficial)
      endif

c ==== if direct match failed, try to match with atom reordering  ==== xxxa new 4/15/03
      ireordermatch=0
      if((iofficialmatch.eq.0).and.(irecognize.eq.2)
     x                        .and.(parallelmode.eq.'spectad')) then
         if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
         call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
         cputimehold=cputime
         do i=2,nofficial
            if(ireordermatch.eq.0) then
               if(abs(statedata(i)%emin-emin).le.(1.d-4)) then    ! hartree
                 call reordermatch(natom,nmove,xyz1,itype,
     x                 statedata(i)%xyzmin,reordercrit,imatch,kreorder)    ! these two are returned
                 if(imatch.ne.0) then
                   ireordermatch=i
                   if((ipr.ge.5).and.(rank_l.eq.0))
     x                  write(6,*)'reordering matchup:',istate,i
                 endif
               endif
            endif
         enddo
         if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
         call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
         cputimereorder=cputimereorder+cputime-cputimehold
      elseif((iofficialmatch.eq.0).and.(irecognize.eq.2)) then
         if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
         call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
         cputimehold=cputime
         do i=1,nofficial
            if(ireordermatch.eq.0) then
               if(abs(eofficial(i)-emin).le.1.d-4) then    ! hartree
                 call reordermatch(natom,nmove,xyz1,itype,
     x                    xyzofficial((i-1)*3*natom+1),reordercrit,
     x                                             imatch,kreorder)    ! these two are returned
                 if(imatch.ne.0) then
                   ireordermatch=i
                   if((ipr.ge.5).and.(rank_l.eq.0))
     x                  write(6,*)'reordering matchup:',istate,i
                 endif
               endif
            endif
         enddo
         if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
         call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
         cputimereorder=cputimereorder+cputime-cputimehold
      endif

      if(ipr.ge.4)  then
       if(rank_l.eq.0) write(6,*)'istate,ioffmatch,ireordmatch:',istate
     x              ,iofficialmatch,ireordermatch
      end if

        ! Set current state to iofficial = 1
        if(parallelmode.eq.'spectad')then
          iofficialtmp=1
          call vecmov(xyz1,xyzofficial((iofficialtmp-1)*3*mxatom+1),3
     +          *natom)
          eofficial(iofficialtmp)=emin
          if(irecognize.gt.0)then
            call vecmov(xyz1,statedata(iofficialtmp)%xyzmin(1),3
     +          *natom)
            statedata(iofficialtmp)%emin = emin
            if(iofficialmatch.eq.0)then
              iofficial=iofficialtmp
            endif
          endif
        endif

c ==== still no matchup - create a new official state ====
          if(iofficialmatch.eq.0 .and. ireordermatch.eq.0) then
c            call gcalc(natom,nmove,xyz1,itype, taxes,
c     x          ietype, maxtyp,rcut, energy,grad, hyperratio)

           if(parallelmode.ne.'spectad')then
            nofficial=nofficial+1
            if (nofficial.gt.nofficialmax) then
              write(6,*) 'STOP: nofficial.gt.nofficialmax'
              stop 'nofficial.gt.nofficialmax'
            endif
            iofficial=nofficial
            call vecmov(xyz1,xyzofficial((iofficial-1)*3*mxatom+1),3
     +          *natom)
c            emin=energy
            eofficial(iofficial)=emin
            if(ipr.ge.1) then
              if(rank_l.eq.0) write(6,141) istate,emin*27.21d0
 141          format(' minimum no.',i5,' stored; energy =',f16.4,' eV')
            end if
            if((ipr.ge.3).and.(rank_l.eq.0))
     x          write(6,*) 'no match, nofficial increased to ',nofficial
           endif

c ==== match found to an official state (possibly with reordering) - use it =====
          else
            if(iofficialmatch.ne.0) then
              iofficial=iofficialmatch
            else if(ireordermatch.ne.0) then
              iofficial=ireordermatch    ! xxxa 4/15/03 - for this case, we will reorder atoms for runningstate
              if((ipr.ge.2).and.(rank_l.eq.0))
     x          write(6,*) 'reordering matchup found',istate, iofficial
            else
              stop 'problem - iofficial matching goofed up'
            end if
            irevisited=1
            if(parallelmode.ne.'spectad')then
             if(isneighshare.gt.0) isneighshare = 0
             call getstate(casename,iofficial,natomx,nmovex,taxes,ntype,
     x          itype,barrierminev,barrierminevss,dimerbarrierminev
     +          ,bminlowerbound,emin,templowtemp,temphightemp
     +          ,timehighprev,timehighnow,timelowprev,tickprev
     +          ,ieventlast,ineedblack,ineedwarm,ineedgauss,ibadtransyes
     +          ,freqlogsmin,nnegmin,nprod,nattempt,nneigh,ntickmax
     +          ,nneighmax,mxatom,lenxyzneighs,xyz,xyz1,xyz4,xyzneighs
     +     ,xyzneighsads,rdimermodes,esaddle,eneigh,nattempted,naccepted
     +          ,ninvolved,nebfound,ticknext,ntick,barrierev,barrevev
     +          ,prefac,isynthclass)
              if(natomx.ne.natom) stop 'natomx.ne.natom after getstate'
              if(nmovex.ne.nmove) then
                stop 'nmovex.ne.nmovex after getstate'
              endif
            else
               !call check_statedata(iofficial,natom,nneigh)
               call lp_getstate(iofficial,natom,emin,xyz1
     +          ,nneigh,nneighmax,mxatom,lenxyzneighs
     +          ,xyzneighs,eneigh,barrierev,barrevev,prefac
     +          ,barrierminev,nofficialmax,freqlogsmin,nnegmin,nprod
     +          ,timehighnow,timehighprev,timelowprev,tickprev
     +          ,ieventlast,betalow,betahigh,fstar,dt,nmd,itype,taxes)
               lfirstinstate = .false.
               timehighnowI=timehighnow
               call vecmov(xyz1,xyz,3*natom)
               templowtemp = templow
               temphightemp = temphigh
               ineedblack=0
               ineedwarm=nwarm+nthermalize
               ineedgauss=1
               bminlowerbound=0.0d0
               barrierminevss=0.0d0
               ibadtransyes=0

            endif
            if (templowtemp.ne.templow) then
               if((ipr.ge.4).and.(rank_l.eq.0)) write(6,*)
     +             'templow different than previously, resetting state'
cc               timehighprev=0.0d0   ! commented out 9/22/04 - see new hightemp ifblock below
cc               timehighnow=0.0d0   ! commented out 9/22/04 - see new hightemp ifblock below
               timelowprev=0.0d0
               tickprev=0.0d0
               ieventlast=0
               timehighprev=timehighprev+timehighnow ! xxx need this to keep real total time high
               timehighnow=0.0d0 ! xxx this might be too drastic
               ineedblack=0
               ineedwarm=nwarm+nthermalize
               ineedgauss=1
               bminlowerbound=0.0d0
               barrierminevss=0.0d0
               ibadtransyes=0
               do izero=1,nneighmax
                  naccepted(izero)=0
cc                  nattempted(izero)=0  ! commented out 9/22/04 - see new hightemp ifblock below
                  ntick(izero)=0
                  isynthclass(izero)=0
                  do idum=1,ntickmax
                     ticknext(idum,izero)=0.0d0
                  enddo
               enddo
               call vecmov(xyz1,xyz,3*natom) ! start from minimum of state
            endif
            if (temphightemp.ne.temphigh .and. ivariabletemp.eq.0) then
               timehighprev=0.0d0
               timehighnow=0.0d0
               do izero=1,nneighmax
                  nattempted(izero)=0
               enddo
            else if(ivariabletemp.eq.1) then
               temphigh=temphightemp
            endif

            if(ipr.ge.4) then
              if(rank_l.eq.0) write(6,142) istate,iofficial,emin*27.21d0
 142          format(' minimum no.',i5,' same as official ',i5,
     x               '; energy =',f16.4,' eV')
             if(rank_l.eq.0) write(6,144)
     x                     emin*27.21d0,barrierminev,timehighnow,
     x                     timehighprev,timelowprev,ibadtransyes
 144          format(' previous info: emin= ',f16.4,
     x             ', barrierminev=',f16.4,
     x             ', timehighnow=',1pd12.4,
     x             ', timehighprev=',1pd12.4,
     x             ', timelowprev=',1pd12.4,
     x             ', ibadtransyes=',i5)
            end if
c === fix betahigh in case temphigh different than last state ===
            if (ivariabletemp.eq.1) then
c xxxb parameter ratioet hardwired for now
               if (tadmode.eq.'dimeretad'.or.
     +             tadmode.eq.'synthdimeretad'.or.
     +             (tadmode.eq.'protecteddimeretad'.and.itrans.gt.1).or.
     +             tadmode.eq.'userdimeretad'
     +             ) then
                  ttherm(1)=temphigh
                  betahigh=1.0d0/(boltzev*temphigh)
                  call mdreset(1)
                  if((ipr.ge.4).and.(rank_l.eq.0))
     +              write(6,*) 'temphigh set to ',temphigh,
     +              ' for state ',istate
               endif
            endif
         endif
          ilabel=iofficial
       else
          stop 'irecognize not recognized'
       end if
       ngradcallminimize=ngradcallminimize+ngradcall-ngradcallhold
       call flushtad(6)

       if((parallelmode.eq.'spectad').and.(irecognize.gt.0))then
            if((ipr.ge.2).and.(rank_l.eq.0))then
              filn=''
              write(filn,*) spawnID
              filn='sp.init.'//trim(adjustl(filn))//'.dat'
              call storefile(natom,xyz1,pxyz,itype,taxes,filn)
            endif
       endif

      !! If this is state #1: make sure we arn't blocking anything yet:
      istatenow = iofficial
      statedata(1)%nneigh_block   = 0
      statedata(1)%ibondneigh1(:) = 0
      statedata(1)%ibondneigh2(:) = 0
      statedata(1)%ibondneigh3(:) = 0
      statedata(1)%dbondneigh12(:)  = 0.d0
      statedata(1)%dbondneigh13(:)  = 0.d0
      statedata(1)%dbondneigh12i(:) = 0.d0
      statedata(1)%dbondneigh13i(:) = 0.d0

CCCCCCC If deposition events are possible, make fake neighbor for this:
      ineighdep = 0
      if((parallelmode.eq.'spectad').and.(dodepositionrun)) then
       ineighdep = 1
       if(nneigh.eq.0) nneigh = 1
      endif

C ===== ++++++ ----- ++++++ ----- ++++++ ----- ++++++ ----- ++++++ ----- ++++++ =====
C ===== ++++++ ----- ++++++ ----- ++++++ ----- ++++++ ----- ++++++ ----- ++++++ =====
C =====                RJZ ---- Should Do SpecSynthTAD here if we can:
C ===== ++++++ ----- ++++++ ----- ++++++ ----- ++++++ ----- ++++++ ----- ++++++ =====
C ===== ++++++ ----- ++++++ ----- ++++++ ----- ++++++ ----- ++++++ ----- ++++++ =====

      time_to_stop = 1.d30
      n_states_skipped = 0
      !maxsynthrate = 0.d0
      if((parallelmode.eq.'spectad').and.
     x                  (irecognize.gt.0)) then

         if((irecognize.gt.1).and.(isneighshare.gt.0))then
           isneighshare = 0
           ! other recognition routines might not work correctly...
         endif

         if(.true.) isneighshare = 0 !! Turned off to debug other things...

         ! Skip states with zero work if possible
         if(isneighshare.gt.0)then

           ! Synthesize neighbors and spawn if possible:
           call lp_synthesize_neighbors(natom,emin,xyz1
     +          ,nneigh,nneighmax,mxatom,lenxyzneighs
     +          ,xyzneighs,eneigh,barrierev,barrevev,prefac,barrierminev
     +          ,nofficialmax,freqlogsmin,nnegmin,nprod
     +          ,timehighnow,betalow,betahigh
     +          ,ineighdep,currenttadtimestop,transcrit,tadtime
     +          ,ntick,ticknext,xmintimelowdrawn
     +          ,jneighlow,isynthclass,isneighshare)

         else
           call lp_kMC(iofficial,natom,emin,xyz1
     +          ,nneigh,nneighmax,mxatom,lenxyzneighs
     +          ,xyzneighs,eneigh,barrierev,barrevev,prefac,barrierminev
     +          ,nofficialmax,freqlogsmin,nnegmin,nprod
     +          ,timehighnow,timehighprev,timelowprev,tickprev
     +          ,betalow,betahigh,fstar,dt,nmd
     +          ,ineighdep,currenttadtimestop,transcrit,tadtime
     +          ,ntick,ticknext,xmintimelowdrawn
     +          ,jneighlow,isynthclass,nattempted,uncertainty)
         endif
         istatenow = iofficial
         call MPI_BARRIER(force_comm,ier)
         if(n_states_skipped.gt.0) then
           timehighnowI=timehighnow
           call vecmov(xyz1,xyz,3*natom)
           call vecmov(xyz,xyz4,3*natom)
           ieventlast = 0
         endif

         if(jneighlow.gt.0)then
           if((spawnDepth.lt.max_depth).and.
     x        ((tadtime+xmintimelowdrawn)
     x        .lt.currenttadtimestop))then

             ldepspawn = 0
             if(isneighshare.gt.0)then
               istatespawn = 0
             else
               istatespawn = statedata(iofficial)%neighlabel(jneighlow)
               if(istatespawn.eq.iofficial)then
                 if(rank_f.eq.0)then
                   write(*,*) 'A - ERROR at SpawnID',SpawnID,
     x             'trying to spawn same state: iofficial=',iofficial,
     x             '- istatespawn=',istatespawn
                   call flush(6)
                 endif
                 call MPI_BARRIER(force_comm,ier)
                 istatespawn = 0
                 statedata(iofficial)%neighlabel(jneighlow) = 0
               endif
             endif
             isneighshare = 0
             if((nneigh-ineighdep-1).gt.0)then
               isneighshare = iofficial
             endif
             ! SpecTAD -- Create New Spawn For Synthetic Event!
             call lp_spawn_update(1,xmintimelowdrawn
     x          ,tadtime,xyzneighs(3*natom*jneighlow+1)
     x          ,nmove,istatespawn,ldepspawn,isneighshare)

             lsynthspawned = .true.
             if(rank_f.eq.0)
     x          write(6,*) 'SpawnID ',spawnID,
     x         ' SPAWNING SYNTHETIC EVENT!! with',
     x         ' xmintimelowdrawn = ',xmintimelowdrawn,' in iofficial ',
     x         istatespawn
             if(rank_f.eq.0) call flush(6)

           elseif((spawnDepth.lt.max_depth).and.(dodepositionrun)
     x             .and.(.not.skipdeposition))then

             ! Time to spawn new deposition
             if(rank_f.eq.0) xyzdep(1) = prngen(0)*taxes(1)
             if(rank_f.eq.0) xyzdep(2) = prngen(0)*taxes(2)
             call MPI_BCAST(xyzdep(1),2,MPI_REAL8,0,force_comm,ier)
             highestatom = 0.d0
             do idum = 1, nmove
               if(xyz1((idum-1)*3+3).gt.highestatom)then
                 dxrz = abs( xyz1((idum-1)*3+1) - xyzdep(1) )
                 dyrz = abs( xyz1((idum-1)*3+2) - xyzdep(2) )
                 !if((dxrz.lt.(3.0)).and.(dyrz.lt.(3.0)))then
                   highestatom = xyz1((idum-1)*3+3)
                 !endif
               endif
             enddo
             xyzdep(3) = highestatom
             call vecmov(xyz1(1),xyzdep(4),3*nmove)
             if(rank_f.eq.0)
     x          write(6,*) 'SpawnID ',spawnID,
     x         ' SPAWNING DEPOSITION EVENT!! at',
     x         ' currenttadtimestop = ',currenttadtimestop
             ldepspawn = 1
             isneighshare = 0
             call lp_spawn_update(1,currenttadtimestop-tadtime
     x          ,tadtime,xyzdep,nmove+1,0,ldepspawn,isneighshare)

             ntick(ineighdep) = 1
             ticknext(1,ineighdep) =
     x         (currenttadtimestop-tadtime)+tickprev
             xmintimelowdrawn = (currenttadtimestop-tadtime)

           else
             previousbest = xmintimelowdrawn
           endif

!           time_to_stop =
!     x      (1.d0/fstar)*(fstar*(xmintimelowdrawn+tickprev))
!     x                        **(betahigh/betalow) - timehighnow
!           nblocksneeded = NINT(time_to_stop/(real(nmd,8)*dt))
!           write(6,*)'SpawnID ',spawnID,' needs ',nblocksneeded
!     x                ,' MD blocks to finish... '

         endif
      endif

      ! Check if we have been killed
      if(rank_f.eq.0) call lp_check_kill()
      call MPI_BCAST(activated,1,MPI_LOGICAL,0,force_comm,ier)
      call MPI_BCAST(spawnwrite,1,MPI_LOGICAL,0,force_comm,ier)
      if(.not.activated) goto 611 !CRZ

C ===== ++++++ ----- ++++++ ----- ++++++ ----- ++++++ ----- ++++++ ----- ++++++ =====
C ===== ++++++ ----- ++++++ ----- ++++++ ----- ++++++ ----- ++++++ ----- ++++++ =====
C =====                RJZ ---- END OF SpecSynthTAD Spawn
C ===== ++++++ ----- ++++++ ----- ++++++ ----- ++++++ ----- ++++++ ----- ++++++ =====
C ===== ++++++ ----- ++++++ ----- ++++++ ----- ++++++ ----- ++++++ ----- ++++++ =====

        ! If we have ParRep cores, give them the coords
        if(nprocs_l.gt.nforcecores)then
          call pr_send_coords(natom,xyz1)
        endif

!        ! Startup fresh lammps run
!        if((spawnid.gt.1).and.(ietype.ge.100).and.(ietype.lt.200))then
!          call MPI_BARRIER(force_comm,ier)
!          call lmp_setup(natom,nmove,ntype,maxtyp,potnam
!     +        ,ietype,xyz1,taxes,itype,rcut,amu
!     +        ,templow,temphigh,.true.,dt)
!          call MPI_BARRIER(force_comm,ier)
!          call vecmov(xyz1,xyz,3*natom)
!          call vecmov(xyz,xyz4,3*natom)
!        endif

        if(rank_f.eq.0) write(*,'(A,I9,A,ES16.6,A)')
     +     'spawnID ',spawnID,' has emin= ',emin*27.21d0,' (eV)'

c ========  try to find an energy-equivalent state ========

c see if there is an earlier state with an equivalent energy - "equiv" state
c NOTE - the point of this is that we can use the emin (in ETAD mode) from
c        an energy-equivalent state
c  the variables with "equiv" in them refer to the equiv-energy checking/matching
c       iequivcheck - controls whether to do these checks
c       eequivcrit = threshhold for considering two states different
c       iequiv = equiv energy case that this state matches
c       nequiv = number of equiv energies found/stored so far
c       eequiv() = array of equiv energies
c       bminequiv() = array of barrierminev values for equiv states
c       tprevequiv() = array of timehighprev values for equiv states
c       nequivmax = dimension of the above arrays
c in the future, we may store iequiv with the official state, but for now,
c let's keep them decoupled

        if(iequivcheck.ge.1) then
          iequivmatch=0
          iumatch=0
          do i=1,nequiv
            if(abs(emin-eequiv(i)).le.eequivcrit) then
              iumatch=iumatch+1
              iequiv=i
            endif
          enddo
          if(iumatch.gt.1) stop 'more than one equiv energy match!'
          if(iumatch.eq.1) then
            dimerbarrierminev=dimbminevequiv(iequiv) ! xxx maybe not ideal place, but better than before
c        xxxxa    could this ever cause a problem?
            if((ipr.ge.4).and.(rank_f.eq.0))
     x          write(6,171) istate,iequiv
 171        format('state ',i5,' matches equiv energy state ',i5)
            iequivmatch=1
          else
            nequiv=nequiv+1
            iequiv=nequiv
            if(nequiv.gt.nequivmax) stop 'increase nequivmax'
            eequiv(iequiv)=emin
            bminequiv(iequiv)=barrierminev
            tprevequiv(iequiv)=timehighprev
            dimbminevequiv(iequiv)=dimerbarrierminev
            if((ipr.ge.2).and.(rank_l.eq.0))
     x        write(6,172) istate,iequiv,emin
 172        format('state ',i5,' is a new equiv energy state', i5,
     x             '   e =',f16.6,' h')
          endif
        else if(iequivcheck.eq.0) then
          iequiv=0
        else
          if(rank_l.eq.0)
     x      write(6,*) 'stopping - iequivcheck must be 0 or 1'
          stop 'stopping - iequivcheck must be 0 or 1'
        endif

c === store this official state minimum as a cluster file called <case>.min.<iofficial>.dat ===
      !if (irevisited.eq.0) then ! xxx is this just how we want this?
        if(lstore_min) then
         filnam=trim(filstrt)//casename(1:lcase)//'.min.'//
     x          chri(ilabel)//'.dat'
         call storefile(natom,xyz1,pxyz,itype,taxes,filnam)
        endif
      !end if

      !CRZ -- Store Running state with spawnID label
      if((parallelmode.eq.'spectad').and.(irecognize.eq.0))then
         filn=''
         write(filn,*) spawnID
         filn='sp.init.'//trim(adjustl(filn))//'.dat'
         call storefile(natom,xyz1,pxyz,itype,
     +           taxes,filn)
      endif

c === store this minimum as <case>.runningstate.<istate>.dat ===
      if(lstore_min) then
        if(istate.le.nrunninglimit) then
         if(ireordermatch.ne.0) then
           call reorder(natom,nmove,xyz1,kreorder, xyz5)
         else
           call vecmov(xyz1,xyz5,3*natom)
         end if
         lc=lcase
         filnam=trim(filstrtLP)//casename(1:lc)//'.runningstate.'//
     x          chri(istate)//'.dat'
         call storefile(natom,xyz5,pxyz,itype,taxes,filnam)
        else
         if(rank_l.eq.0)then
          write(6,*) '********* running file limit reached ********'
         endif
        end if
      end if

c === also store this minimum as <case>.latest.dat ===
      if(lstore_min) then
        filnam=trim(filstrtLP)//casename(1:lcase)//'.latest.dat'
        call storefile(natom,xyz1,pxyz,itype,taxes,filnam)
      end if

c ===== zero out timehighnow and ticknext array, if appropriate =====

      if((parallelmode.ne.'spectad').or.(irecognize.eq.0))then
      if(   tadmode.eq.'tad'
     x .or. tadmode.eq.'newtad'
     x .or. tadmode.eq.'dimeretad'
     x .or. tadmode.eq.'protecteddimeretad'
     x .or. tadmode.eq.'userdimeretad'
     x    ) then
         timehighnow=0.0d0
         do i=1,nneigh
            ntick(i)=0
            ticknext(1,i)=0.0d0
         enddo
      endif
      endif

      if(tadmode.eq.'synthtad'.or.tadmode.eq.'synthdimeretad') then
        isynlast=0
        if(ieventlast.gt.0) isynlast=isynthclass(ieventlast)
        if(isynlast.eq.0) then
          timehighnow=0.0d0
          do i=1,nneigh
             if (isynthclass(i).eq.0) then
                ntick(i)=0
                ticknext(1,i)=0.0d0
             endif
          enddo
        end if
      endif

      iblock=0
      ijump=0
      ineigh=0
      newneighbor=0
c === jump below to do a zero-block acceptance check and then come right back to 222 ===

      timelowevent = 1.d99 ! crz -- make sure unseen 1st event took 'very' long
      if (iblock.eq.0) goto 333

c === jump back up to here to continue integrating in this basin ===
 222  continue

c ======= begin iblock.eq.0 if-block to do diag, dimer and warming on entering state =======
      if (iblock.eq.0) then

c ===== diagonalize hessian to obtain normal modes at this minimum =====

      ngradcallhold=ngradcall
      if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
      call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
      cputimehold=cputime
      icnt = 0
      if((ivineyard.eq.1).and.(nnegmin.gt.0)) then ! xxxx skip
         tmgcrit = gcrit
         tmdvcrit = dvcrit
         tmdrcrit = drcrit
 116  continue
         igeometry = 0 ! minimum
         nnegmin = 0
         freqlogsmin = 0.d0
         nprod = 0
         nmove_tmp = nmove
         if(nmove.eq.natom) nmove_tmp = nmove-1
         call dovineyard(natom,nmove_tmp,xyz1,itype,taxes,ietype,maxtyp
     +      ,rcut,amu,e,grad,nneg,nprod,freqlogsmin,igeometry,fdum,ndum
     +      ,fdum2,ierr,ipr)   ! xxx fdum,fdum2,ndum dummy variables

         nnegmin = nneg
         emin = e
         if(nnegmin.gt.0)then
           if(rank_f.eq.0)
     x      write(*,*) 'SpawnID ',spawnID,' is in state with nnegmin = '
     +       ,nnegmin,' -- adding noise and minimizing again! icnt - '
     +       ,icnt
           call flushtad(6)
           do idum=1,3*nmove
             if(rank_f.eq.0)then
               xyz1(idum)=xyz1(idum)+(0.01d0)*gasdev(0)
             endif
           enddo
           call MPI_BCAST(xyz1,3*nmove,MPI_REAL8,0,force_comm,ier)
           icnt = icnt + 1
           tmgcrit = tmgcrit*0.1
           tmdvcrit = tmdvcrit*0.1
           tmdrcrit = tmdrcrit*0.1
           call descent_check(natom,nmove,0,xyz1,xyzofficial
     +       ,itype,taxes,ietype,maxtyp,rcut,hyperratio,barrevev
     +       ,ereverseuse,iofficialmatch,irotalign,nsteep1,gfac
     +       ,transcrit,itranscrit,tmgcrit,tmdvcrit,tmdrcrit
     +       ,ipr,0,ifundcuse)
           call vecmov(xyz1,xyz,3*natom)
           if(icnt.ge.25) then
             filn=''
             write(filn,*) spawnID
             filn='sp.init.'//trim(adjustl(filn))//'.badmin.dat'
             call storefile(natom,xyz1,pxyz,itype,
     +           taxes,filn)
             STOP 'THIS MINIMUM IS NO GOOD !!!'
           endif
           goto 116
         endif
      endif
      ngradcallvineyard=ngradcallvineyard+ngradcall-ngradcallhold
      if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
      call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
      cputimevineyard=cputimevineyard+cputime-cputimehold
      call flushtad(6)

c move xyz1 into first position in xyzneighs super-array

      call vecmov(xyz1,xyzneighs(1),3*natom)

c == store unpbc'd coordinates of this minimum as a cluster file ==

      if(lstore_minunpbc.and.irevisited.eq.0) then
         filnam=trim(filstrtLP)//casename(1:lcase)//'.min.unpbc.'
     +       //chri(ilabel)//'.dat'
        call pbcuncumfake(natom,xyz1,taxes,cdxyz)
        call storefile(natom,xyz1,pxyz,itype,taxes,filnam)
        call pbccumfake(natom,xyz1,taxes,cdxyz)
      end if

!      !CRZ -- Store UPBC file with initial Positions...
!      if((.false.).and.(parallelmode.eq.'spectad'))then
!         filn=''
!         write(filn,*) spawnID
!         filn='sp.init.'//trim(adjustl(filn))//'.unpbc.dat'
!         call storefile(natom,xyz1,pxyz,itype,
!     +           taxes,filn)
!      end if

c this is here so that we get to see the final minimum, etc., before quitting
      if(parallelmode.ne.'spectad')then
      if(itrans.eq.naccept+1) then    ! xxx make this say istate.eq.naccept   xxxxa - OK here?
        if(rank_f.eq.0) write(6,*) ' '
        if(rank_f.eq.0) write(6,*) 'all requested transitions found'
        istopping=1
        go to 600
      else
        istopping=0
      end if

      if(ispecialstop.eq.1) then  ! certain final state reached (see below)
        if(rank_f.eq.0) write(6,*) ' '
        if(rank_f.eq.0) write(6,*) 'all requested transitions found'
        istopping=1
        go to 600
      else
        istopping=0
      end if
      endif

c ======= check if we are going to go into super-synth mode =======
      if (tadmode.eq.'synthdimeretad'.and.idosupersynth.eq.1.and.
     +    iblock.eq.0) then
         if(rank_f.eq.0)
     +     write(6,*)'LIST OF BARRIERS FOR STATE',iofficial
         if(rank_f.eq.0) write(6,*)'E   #attempted  #accepted'
         isupersynthdimers=0
         ndimerusesuper=0
         isupersynth=0
         icountsynth=0
         ienoughtosuper=0
         do ifm=1,nneigh
            if(isynthclass(ifm).eq.0) then
               if(rank_f.eq.0)
     +           write(6,3334)barrierev(ifm),nattempted(ifm)
     +             ,naccepted(ifm)
            else
               icountsynth=icountsynth+1
               if(rank_f.eq.0)
     +           write(6,3335)barrierev(ifm),nattempted(ifm),
     &             naccepted(ifm)
               if(naccepted(ifm).gt.nacceptsuper)
     +             ienoughtosuper=ienoughtosuper+1
            endif
         end do
         if(ienoughtosuper.ge.icountsynth.and.icountsynth.gt.0)
     +       then
            if(rank_f.eq.0) write(6,*)
     +          'READY TO GO SUPER, icountsynth,ienoughtosuper'
     +          ,icountsynth,ienoughtosuper
            isupersynth=1
         else
            if(rank_f.eq.0) write(6,*)'NOT READY TO GO SUPER,',
     +          ' icountsynth,ienoughtosuper',icountsynth
     +          ,ienoughtosuper
            isupersynth=0
         endif

 3334    format(f12.5,1x,'FRANZ Normal',1x,i6,1x,i6)
 3335    format(f12.5,1x,'FRANZ Synthetic',1x,i5,1x,i6)
         call findmaxbarriers(nneigh,barrierev,isynthclass
     +       ,barriermaxevs)
         call findminbarrierns(nneigh,barrierev,isynthclass
     +       ,naccepted,barrierminevns)
         call findminbarrier(nneigh,barrierev,barrierminev,barrevev
     +       ,ereverse)
         ratiohowmany=
     &       exp((barrierminevns-barrierminev)/(boltzev*templow))
         if(rank_f.eq.0) write(6,*)'lowest barrier:',barrierminev
         if(barriermaxevs.gt.1d-5.and.barrierminevns.gt.1d-5) then
          if(rank_f.eq.0)then
            write(6,*)'largest synthetic barrier:',barriermaxevs
            write(6,*)'smallest non-synthetic',
     +          ' or non-accepted barrier:',barrierminevns
            write(6,*)'how many transitions before',
     +          ' accepting the I N-S ?',ratiohowmany
          endif
         endif
         if(isupersynth.gt.0) then
            if(ratiohowmany.gt.ratiosuper) then
               if(rank_f.eq.0) write(6,*)'SWITCHING TO SUPER WITH E='
     +             ,barrierminevns,', iofficial=',iofficial
               if (abs(barrierminevns-barrierminevss).gt.1e-4) then
                  ndimerusesuper=10
                  isupersynthdimers=1
               endif
            else
               if(rank_f.eq.0)
     +           write(6,*)'I CAN LIVE WITHOUT GETTING SUPER'
               isupersynth=0
            endif
         endif
      endif                     ! end of super synthetic if-block

c ============== do some dimer searches ============

      ndimeretaduse=0
      ndimeruse=0
      ebarrierhole=10.0d10
c xxxb clean up this if-block
      if((tadmode.eq.'dimeretad'.or.tadmode.eq.'synthdimeretad'.or
     +    .tadmode.eq.'userdimeretad'.or
     +    .(tadmode.eq.'protecteddimeretad'.and.itrans.gt.1)).and
     +    .((iequivcheck.eq.1.and.iequivmatch.eq.0).or.(iequivcheck.eq.0
     +    .and.irevisited.eq.0)))ndimeretaduse=ndimeretad
      if(irevisited.eq.0) ndimeruse=ndimer
      ndimeruse=max(ndimeruse,ndimeretaduse)
      if (itrans.eq.1) ndimeruse=ndimeruse*2
      if (ipr.ge.10) then
        if(rank_f.eq.0)
     +   write(6,*) 'dimer search debug: tadmode: ',tadmode,', itrans: '
     +       ,itrans,', iequivcheck: ',iequivcheck,', iequiv: ',iequiv,
     +       ', irevisited: ',irevisted,', ndimeretaduse: '
     +       ,ndimeretaduse,', ndimeruse: ',ndimeruse
      endif
c      ndimerfull=int(float(ndimeruse)*0.2d0) ! xxxxx want to set as input parameter
c      ndimerfull=ndimeruse ! xxxx do all dimers fully
c      ndimernew=ndimeruse
      if (isupersynthdimers.eq.1) ! super-synth dimers (fixed number for now)
     +    ndimeruse=max(ndimeruse,ndimerusesuper)

      ndimerdoneold=0
c      iflag=0
      iflag2=0
      id=0
      idimerfocused=0
      idimertotal=0
      inewholebarrier=0
      inewminbarrier=0
      idimerstop=0
      idonewdimers=0            ! do new dimers on all atoms meeting criteria (not typically necessary)
      if (nneighlast.eq.0) idonewdimers=1 ! do new dimers if no neighbors last state (ie first state)
      if(ndimeruse.gt.0) then

c         ngradcallhold=ngradcall
         if((ipr.ge.4).and.(rank_f.eq.0))
     +        write(6,*) 'beginning dimer searches',ngradcall
         if(rank_f.eq.0) write(6,*) 'doing ',ndimeruse,' dimer searches'
c         do id=1,ndimeruse
         do while (idimerstop.eq.0)
            id=id+1
c            if (id.le.ndimerfull.or.(id.gt.nneighlast.and.id.le
c     +          .nneighlast+ndimerfull)) then
c               idimerbail=0
c            else
c               idimerbail=1
c            endif
            idimerbail=1 ! xxxx bailing (screening) on everything for now
            if ((ipr.ge.10).and.(rank_f.eq.0))
     +          write(6,*) 'ndimerfull: ',ndimerfull
     +          ,'; idimerbail: ',idimerbail,'; dimerbarrierminev: '
     +          ,dimerbarrierminev,'; nneighlast: ',nneighlast
     +          ,'; ndimerdoneold: ',ndimerdoneold
            imode=0
            ngradcallhold=ngradcall
            ndelgradvine=0

 411        if (itrans.ne.1.and.ndimerdoneold.lt.nneighlast.and
     +          .isupersynthdimers.ne.1) then ! reuse dimers from previous state
               ndimerdoneold=ndimerdoneold+1
               call vecmov(xyzneighsads(3*natom*ndimerdoneold+1)
     +             ,xyz2,3*natom) ! put previous saddle in xyz2
               call vecmov(rdimermodes(3*natom*ndimerdoneold+1)
     +             ,rdimer,3*natom) ! put previous dimer direction in rdimer
               call shiftsaddle(natom,nmove,rdimer,xyz2,xyzprev,xyz1
     +             ,transcrit,ierr) ! shift old saddle to new geometry
               rdimermag=vecdot(rdimer,rdimer,3*natom)
               if (ierr.eq.1.and.ietype.eq.0) then ! atoms too close for EAM
                  if((ipr.ge.6).and.(rank_f.eq.0)) write(6,*)
     +                'shiftsaddle caused atoms to be too close,',
     +                ' redoing'
                  goto 411
               endif
               iranx=0       ! do not randomize displacement
               if (rdimermag.lt.0.001d0) then ! no previous dimer direction, so do random direction
                  irandir=1     ! randomize direction
               else
                  irandir=0     ! do not randomize direction
               endif
               imode=2          ! do not try to redo dimer internally, abort back to here
               idimertype=2
            else                ! not using previous saddles
               if (iflag2.eq.0) then
                  if((ipr.ge.4).and.(rank_f.eq.0)) write(6,*)
     +                'done with reused dimers: ',id
                  iflag2=1
               endif
               call vecmov(xyz1,xyz3,3*natom) ! move minimum into xyz3
               if (nneighlast.gt.0.and.idimerfocused.lt.ndimeruse) then
                  iranx=2       ! focus new dimers on "hole" region; in protecteddimeretad mode, always focus on hole
                  if (inewholebarrier.eq.0) then
                     idimerfocused=idimerfocused+1
                  else
                     idimerfocused=0
                  endif
                  inewholebarrier=0
                  idimertype=3
               else if ((nneighlast.eq.0.or.idimerfocused.ge.10).and
     +                .idimertotal.lt.ndimeruse.and.idonewdimers.eq.1)
     +                then
                  iranx=1 ! randomize new dimers completely
                  if (tadmode.eq.'userdimeretad'.and.itrans.eq.1)
     +                iranx=3
                  if (inewminbarrier.eq.0) then
                     idimertotal=idimertotal+1
                  else
                     idimertotal=0
                  endif
                  inewminbarrier=0
                  idimertype=1
               else
                  idimerstop=1
               endif
               irandir=1
            endif
c            if (id.eq.ndimeruse.and.ndimerdoneold.lt.nneighlast)
c     +          ndimeruse=ndimeruse+1 ! make sure we don't stop too early if have more than ndimeruse old neighbors
c            if (iflag.eq.0.and.ndimerdoneold.ge.nneighlast) then
c               if (ndimeruse-id.lt.ndimernew-1) then
c                  ndimeruse=id+ndimernew ! make sure do at least ndimernew new dimers
c                  if (ipr.ge.4) write(6,*) 'dimers: for this state,',
c     +                ' ndimeruse increased to ',ndimeruse
c               endif
c               iflag=1
c            endif
            if((ipr.ge.4).and.(rank_f.eq.0))
     +          write(6,3339) istate,iranx,id,nneighlast
     +          ,idimertype,idimertotal,idimerfocused,ebarrierhole
     +          ,dimerbarrierminev
 3339       format('dlok: ',7i6,2d20.10)
            call dimer_search(natom,nmove,idimerbail,emin
     +          ,dimerbarrierminev,idimerdof,ndimerdof,rdimerdist
     +          ,ndimercoord,idimeroverunder,rdimerdisp,xyz3,xyz2,xyz1
     +          ,xyzprev,edimer,rdimer,itype,taxes,ietype,maxtyp,rcut
     +          ,hyperratio,barrevev,ereverse,gfac,irotalign,nsteep1
     +          ,ierr,iranx,irandir,transcrit,itranscrit,gcrit,dvcrit
     +          ,drcrit,ipr,imode,ifundc,eig)
            if (ierr.ne.10) then
               dimerbarrierev=(edimer-emin)*27.21d0
               if((ipr.ge.3).and.(rank_f.eq.0))
     +           write(6,51) id,dimerbarrierev
 51            format('dimer_search ',i5,' finds barrier of ',f12.5
     +             ,' eV')
c               if(dimerbarrierev-dimerbarrierminev.lt.-1.d-5) then
               rdimerdiff=(dimerbarrierev-dimerbarrierminev)
     +             /dimerbarrierminev
               if (rdimerdiff.lt.0.0d0.and.abs(rdimerdiff).gt
     +             .rdimerredocrit) then
                  dimerbarrierminev=dimerbarrierev
                  idlowest=id
                  inewminbarrier=1
               else if (rdimerdiff.lt.0.0d0) then
                  dimerbarrierminev=dimerbarrierev
                  idlowest=id
               end if
               if (iranx.eq.2) then ! dimer was focused on hole
                  inewholebarrier=0
                  if(dimerbarrierev-ebarrierhole.lt.-1.d-5) then
                     ebarrierhole=dimerbarrierev
                     inewholebarrier=1
                  endif
               endif
c check if endpoint is in neighborlist, if not, add it.
               ifundcuse=ifundc
               if (ifundcuse.eq.5) ifundcuse=0
               ereverseuse=0.0d0 ! here we are just trying to see if neighbor is in list, don't do reverse barrier stuff
               call descent_check(natom,nmove,1+nneigh,xyz3,xyzneighs
     +             ,itype,taxes,ietype,maxtyp,rcut,hyperratio,barrevev
     +             ,ereverseuse,isame,irotalign,nsteep2,gfac,transcrit
     +             ,itranscrit,gcrit,dvcrit,drcrit,ipr,0,ifundcuse)
               ngradcalldimer=ngradcalldimer+ngradcall-ngradcallhold ! xxxx wrong place?
               if (isame.eq.0) then
                  call gcalc(natom,nmove,xyz3,itype, taxes, ! get energy of endpoint
     x                ietype, maxtyp,rcut, energy2,grad, hyperratio)

                  nneigh=nneigh+1
                  if (nneigh.gt.nneighmax) then
                     if(rank_f.eq.0) write(6,*)
     +                   'STOP: nneigh.gt.nneighmax, too many neighs'
                     stop 'nneigh.gt.nneighmax, too many neighs'
                  endif
                  esaddle(nneigh)=edimer
                  eneigh(nneigh)=energy2
                  barrierev(nneigh)=(edimer-emin)*27.21d0
                  barrevev(nneigh)=(edimer-energy2)*27.21d0
                  if ((ipr.ge.6).and.(rank_f.eq.0)) write(6,*)
     +                'reverse debug: reverse barrier: ',nneigh
     +                ,barrevev(nneigh)
c xxx check for minimum barrier moved to loop above and below
c xxx we have a lot of these checks of array lengths versus parameters,
c should we consildate them into a call to a subroutine?
                  if(3*natom*(1+nneigh).gt.lenxyzneighs) then
                    if(rank_f.eq.0)
     +               write(6,*) 'STOP: nneigh*3*natom.gt.lenxyzneighs,',
     +                   ' too many neighs'
                    stop 'nneigh*3*natom>lenxyzneighs, too many neigh'
                  endif
                  call vecmov(xyz3,xyzneighs(3*natom*nneigh+1),3*natom)
                  call vecmov(xyz2,xyzneighsads(3*natom*nneigh+1),
     +                3*natom)
                  call vecmov(rdimer,rdimermodes(3*natom*nneigh+1),
     +                3*natom)
                  if(lstore_sad) then
                     filnam=trim(filstrt)//casename(1:lcase)//
     x                   '.'//chri(ilabel)//'.sad.'//chri(nneigh)//
     +                   '.dat'
                     call storefile(natom,xyz2,pxyz,itype,taxes,filnam)
                  end if
                  if (lstore_end) then
                     filnam=trim(filstrt)//casename(1:lcase)//
     x                   '.'//chri(ilabel)//'.end.'//chri(nneigh)//
     +                   '.dat'
                     call storefile(natom,xyz3,pxyz,itype,taxes,filnam)
                  endif
c xxx we should calculate vineyard and ninvolved if set for this
C neighbor, or should we check for events actually seen if we have these
C and calculate them then?
                  if (ivineyard.gt.0) then
                     ngradcallholdvine=ngradcall
                    if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
                    call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
                     cputimehold=cputime
                     if (ivineyard.eq.1) then
                        igeometry=1 ! saddle
                        call dovineyard(natom,nmove,xyz2,itype,taxes
     +                      ,ietype,maxtyp,rcut,amu,e,grad,nnegsad,nprod
     +                      ,freqlogssad,igeometry,freqlogsmin,nnegmin
     +                      ,prefac(nneigh),ierr,ipr)
c xxxx should we print warning if prefac < prefacmin?
                        if(ierr.eq.0) then
                           if((ipr.ge.2).and.(rank_f.eq.0)) write(6,*)
     +                         'Vineyard prefac for state ',istate
     +                         ,' to neigh(i.e., saddle) ',nneigh,'  ='
     +                         ,prefac(nneigh),' Hz'
                        end if
                     endif
c xxxx new subset vineyard
                     if (ivineyard.eq.2) then
                        iupper=4
                        if (idynmatloop.eq.0) iupper=0
                        do ibpu=iupper,0,-1
                           rdynmatuse=rdynmatcrit/(3**float(ibpu))
                           call vineyardsubset(natom,nmove,xyz1,xyz2
     +                         ,itype,taxes,ietype,maxtyp,rcut,amu
     +                         ,rdynmatuse,e,grad,prefacx,nelements,ierr
     +                         ,ipr)
                           if (ibpu.eq.iupper) prefacxx=prefacx
                           if(rank_f.eq.0)then
                             write(6,6543) istate,nneigh,rdynmatuse
     +                         ,prefacx,nelements,prefacxx
                           endif
 6543                      format('vineyardsubset: state: ',i4
     +                         ,', neigh: ',i4,', crit: ',f9.4
     +                         ,', prefac: ',1pd16.3,', nelements: ',i5
     +                         ,',prefacall: ',1pd16.3)
                        enddo
                        prefac(nneigh)=prefacxx
                     endif

                     ngradcallvineyard=ngradcallvineyard+ngradcall
     +                   -ngradcallholdvine
c                     ndelgradvine=ndelgradvine+ngradcall
c     +                   -ngradcallholdvine
                     if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
                    call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
                     cputimevineyard=cputimevineyard+cputime-cputimehold
                  endif
                  filnam=trim(filstrtLP)//casename(1:lcase)//'.barriers'
                  ctype='DI'
                  hyperdist=0.0d0 ! for now
                  call writebarrierentry(filnam,timehighnow
     +                ,barrierev(nneigh),betalow,betahigh,istate
     +                ,iofficial,nneigh,prefac(nneigh),ninvolved(nneigh)
     +                ,barrevev(nneigh),ctype,hyperdist,dist4)
               else if (isame.gt.1) then ! check energy of saddle to see if maybe states connected by more than one saddle
                  barriercurrent=(edimer-emin)*27.21d0
                  if (abs(barrierev(isame-1)-barriercurrent).gt.1e-4)
     +                then
                     if(rank_f.eq.0) write(6,*)
     +                   'WARNING: Dimer saddle connecting state '
     +                   ,istate,' with neighbor ',(isame-1)
     +                   ,' seems to have two saddles of energies: '
     +                   ,barrierev(isame-1),barriercurrent
     +                   ,abs(barrierev(isame-1)-barriercurrent)
                  endif
               else
                  if(rank_f.eq.0) write(6,*)
     +                'WARNING: Dimer end point is same as'
     +                ,' initial point'
                  if(rank_f.eq.0) write(6,*)
     +                'WARNING: Ignoring and moving on'
               endif
               if(lstore_dimer) then
                  filnam=trim(filstrtLP)//casename(1:lcase)//
     x                '.'//chri(ilabel)//'.dimer.'//chri(id)//'.dat'
                  call storefile(natom,xyz2,pxyz,itype,taxes,filnam)
               endif
            endif
         enddo

         if(ipr.ge.2) then
          if(rank_f.eq.0) write(6,55) idlowest,istate,dimerbarrierminev
         endif
 55      format('minimum dimer barrier (#',i4,') for state',i5,
     x       ' =',f12.5,' eV')
         if(ipr.ge.10) then
           if(rank_f.eq.0) write(6,*) 'done dimer searches',ngradcall
         endif
c         ngradcalldimer=ngradcalldimer+ngradcall-ngradcallhold ! xxxx not necessary?
c     +       -ndelgradvine
      endif
      call flushtad(6)
c------------------end of dimer section ------------------------

c ====== super-synth dimer check to see if we found any new ones ======
      if (tadmode.eq.'synthdimeretad'.and.idosupersynth.eq.1.and.
     +    iblock.eq.0.and.isupersynth.gt.0) then
         isupersynthdimers=0
         barrierminevnshold=barrierminevns
         call findminbarrierns(nneigh,barrierev
     +       ,isynthclass,naccepted,barrierminevns)
         if (abs(barrierminevns-barrierminevnshold).gt.1e-2) then
            if(rank_f.eq.0) write(6,*) 'supersynth: first DIMER did not'
     +          ,' find a barrier:'
     +          ,barrierminevns,barrierminevnshold,', iofficial='
     +          ,iofficial
            isupersynth=0
            if(rank_f.eq.0)
     +          write(6,*) 'supersynth: SWITCHING to super aborted'
            do ifm=1,nneigh
               if(isynthclass(ifm).eq.0) then
                 if(rank_f.eq.0)
     +            write(6,3334) barrierev(ifm),nattempted(ifm)
     +                ,naccepted(ifm)
               else
                 icountsynth=icountsynth+1
                 if(rank_f.eq.0)
     +            write(6,3335) barrierev(ifm),nattempted(ifm)
     +                ,naccepted(ifm)
               endif
            enddo
         endif
      endif

c === adjust high temperature if requested ===
      if (ivariabletemp.eq.1.and.irevisited.eq.0) then
c xxxb parameter ratioet hardwired for now
         if (tadmode.eq.'dimeretad'.or.tadmode.eq.'synthdimeretad'.or
     +       .tadmode.eq.'userdimeretad'.or
     +       .(tadmode.eq.'protecteddimeretad'.and.itrans.gt.1)) then
            dimerbarrierminevuse=1.0d10
c search for minimum barrier that has a reverse barrier bigger than
c ereverse, so we can set the temperature to a good value to escape from
c the "super" basin
            do i=1,nneigh
               if (barrevev(i).gt.ereverse.and.barrierev(i).lt
     +             .dimerbarrierminevuse) then
                  dimerbarrierminevuse=barrierev(i)
               endif
            enddo
            if (ipr.ge.4) then
               if(rank_f.eq.0) write(6,*)
     +             'Minimum Barrier with Reverse Barrier>Ereverse='
     +             ,dimerbarrierminevuse
            endif
            temphighx=dimerbarrierminevuse/boltzev/ratioet
c xxxxf FRANZ if dimerbarrierminev is too low, no thermalization will be required
            if((dimerbarrierminevuse/(boltzev*templow)).lt.0.5) then
               maxwarmreject=0
               maxblackreject=0
               if(rank_f.eq.0)
     +           write(6,*)'BARRIER VERY LOW: NO THERMALIZATION',
     +             dimerbarrierminev
            else
               maxwarmreject=maxwarmrejectsave
               maxblackreject=maxblackrejectsave
            endif

            if(temphighx.gt.temphighmax)temphighx=temphighmax
            if(temphighx.lt.templow)temphighx=templow+1.d0
c if state is revisited, temphigh should match what was found before, for now
            if (irevisited.eq.1.and.temphighx.ne.temphigh) then ! xxxb do we need to check difference
             if(rank_f.eq.0)then
               write(6,*) 'STOP: temphigh not same as previous visit'
               write(6,*) 'previous temphigh=',temphigh,'; new='
     $              ,temphighx
             endif
             stop 'temphigh not same as previous visit'
            endif
            temphigh=temphighx
            ttherm(1)=temphigh
            betahigh=1.0d0/(boltzev*temphigh)
            call mdreset(1)
            if(rank_f.eq.0)then
              if (ipr.ge.4) write(6,*) 'temphigh set to ',temphigh,
     +           ' for state ',istate
              if(ipr.ge.4) write(6,*)'maxwarm,maxblack',maxwarmreject
     +          ,maxblackreject
            endif
         endif
      endif

      if((ipr.ge.4).and.(rank_f.eq.0))
     +       write(6,*)'begin basin-constrained MD in state',istate

c === assign random thermal momenta  ===

      if(ineedgauss.ne.0) then
        iop=0
        iwrt=-1
        call gaussp(natom,nmove,pxyz,itype,maxtyp,amu,ntype,
     x                moving,mx,my,mz,2.0*temphigh,iop,iwrt)
        call MPI_BCAST(pxyz,3*nmove,MPI_REAL8,0,force_comm,ier)
        call mdreset(0)
      end if


c ======== take some warmup steps, rejecting transitions ========

      ! If we are doing ParRep... Don't 'usually' warm up Master..
      if((ldprcheckdone).and.(nprocs_l.gt.nforcecores)
     x                        .and.(rank_l.lt.nforcecores))then
        ineedwarm = 0
      endif

      call MPI_BARRIER(force_comm,ier)
!      if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID,
!     +      'warmup starting - ineedwarm = ',ineedwarm

c xxx this and blackout code very similar, can combine into subroutine?
      ngradcallhold=ngradcall
      if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
      call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
      cputimehold=cputime
      if(iblock.gt.1) ibadtransyes=0 ! do not reset if just entering as may be from before
      if(ineedwarm.gt.0) then
        ibadtrans=0
 200    continue                  ! redo warming loop

!        ! Block to reduce T if we are in shallow state:
!        if (.false.) then
! 614        continue ! Reduce T_high, and redo warming loop
!            ibadtrans = 0
!            temphigh = max((templow+10.0),(temphigh-100.0))
!            ttherm(1)=temphigh
!            betahigh=1.0d0/(boltzev*temphigh)
!            call mdreset(1)
!            if(rank_f.eq.0)then
!              if (ipr.ge.1)
!     +           write(6,*) 'SpawnID',SpawnID,'temphigh set to ',
!     +           temphigh,' for istate ',istate
!            endif
!            iop=0
!            iwrt=-1
!            call gaussp(natom,nmove,pxyz,itype,maxtyp,amu,ntype,
!     x                moving,mx,my,mz,temphigh,iop,iwrt)
!            call MPI_BCAST(pxyz,3*nmove,MPI_REAL8,0,force_comm,ier)
!            call mdreset(0)
!        endif

        call vecmov(xyz,xyz3,3*natom)

        if((ietype.ge.100).and.(ietype.lt.200).and.doblocks)then
            call mdblock100(natom,nmove,xyz,pxyz,itype,
     x              maxtyp,e,grad,taxes,ietype,amu,temphigh,dt,
     x              thermrate,ineedwarm,iseed,ntype)
            ngradcall = ngradcall + ineedwarm
        else
            blocking_on = .true.
            do i=1,ineedwarm
                call mdstep(natom,nmove,xyz,itype,pxyz,taxes,
     x              ietype, maxtyp,rcut,amu,
     x              maxtherm,itherm,jtherm,ttherm,ctherm,
     x              xtherm,mtherm,ktherm,ntherm,
     x              moving,mx,my,mz,dt,time,thyper,hyperratio,
     x              e,grad) ! xxx time
cdebug          filnam='badwarm.snap.'//chri(i)//'.dat'
cdebug          call storefile(natom,xyz,pxyz,itype,taxes,filnam)
            enddo
            blocking_on = .false.
        endif

        !! RJZ - shiftTest
        if (natom.eq.nmove) then
           call shiftalign(xyz3,xyz,natom)
        endif
        call vecmov(xyz,xyz2,3*natom)

        ! Check if we have been killed
        if(rank_f.eq.0) call lp_check_kill()
        call MPI_BCAST(spawnwrite,1,MPI_LOGICAL,0,force_comm,ier)
        call MPI_BCAST(activated,1,MPI_LOGICAL,0,force_comm,ier)
        if(.not.activated) goto 611 !CRZ

c   check against all known neighbors
c   note that xyzneighs has the original minimum in the first position,
c   the first neighbor in the 2nd position and so on.
        imode=0 ! do not do spring quench
!        if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID,
!     +      'warmup descent_check()'
        call descent_check(natom,nmove,1+nneigh,xyz2,xyzneighs,
     +      itype,taxes,ietype,maxtyp,rcut,hyperratio,barrevev,ereverse
     +      ,imatch,irotalign,nsteep2,gfac,transcrit,itranscrit,gcrit
     +      ,dvcrit,drcrit,ipr,imode,ifundc)
        if (ereverse.gt.0.0d0.and.imatch.gt.1) then ! bpu: added 09/09/2004
           if (barrevev(imatch-1).lt.ereverse) then
              imatch=1
              if (ipr.ge.6) write(6,*)
     +            'reverse debug: warm detected ignored state'
           endif
        endif

        if(imatch.ne.1) then
           ibadtrans=ibadtrans+1

           filnam=trim(filstrtLP)//casename(1:lcase)//'.'//
     +         chri(istate)//'.warmtrans.'//chri(ibadtrans)//'.dat'
           if (ipr.ge.10)
     +         call storefile(natom,xyz2,pxyz,itype,taxes,filnam)
           if(ibadtrans.lt.maxwarmreject) then
             if((ipr.ge.2).and.(rank_f.eq.0))
     x        write(6,*) 'SpawnID',SpawnID,
     x                   '- Transition during warming, redoing!'
             call vecmov(xyz3,xyz,3*natom)
             call gaussp(natom,nmove,pxyz,itype,maxtyp,amu,ntype,
     x           moving,mx,my,mz,2.0*temphigh,iop,iwrt)
             call MPI_BCAST(pxyz,3*nmove,MPI_REAL8,0,force_comm,ier)
             call mdreset(0)
!             ! Redo warmup with LOWER High-Temperature??
!             if ((ibadtrans.lt.2).or.(temphigh.le.templow+10.0))then
!                 goto 200
!             else
!                 goto 614
!             endif
             goto 200
           else
             if (ibadignore.eq.0) then
               write(6,*)'Waiting: Transition during warmup, aborting!'
               is=istate
               filnam=trim(filstrtLP)//casename(1:lcase)//'.'//chri(is)
     +                 //'.badwarm.dat'
               call storefile(natom,xyz2,pxyz,itype,taxes,filnam)
               write(6,*) 'abort-too many warming rejects',ibadtrans
              write(6,*)'see ...badwarm.dat file, named on state number'
               write(6,*) 'not official number'
               stop 'too many warming transitions found'
             else if (ibadignore.eq.1) then
                if(rank_f.eq.0) write(6,*) 'WARNING: Did '
     +              ,ibadtrans
     +              ,' warming rejects, continuing anyway!'
     +              ,'- SpawnID ',SpawnID
                ibadtransyes = 1 ! so we skip the leak check
                cleanwarmup = .false.
                call vecmov(xyz3,xyz,3*natom)
             endif
           end if
        endif
        call MPI_BARRIER(force_comm,ier)
        if((ipr.ge.4).and.(rank_f.eq.0))then
          write(6,*) 'spawnID ',spawnID,' finished warmup'
        endif
      end if
      ngradcallwarm=ngradcallwarm+ngradcall-ngradcallhold
      call flushtad(6)
c------------end of transition-rejecting warmup stage----------------

c quick, one-time force-call timing report

      if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
      call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
      if((.not.ldprcheckdone).and.(ipr.ge.0))then
        callspersec=float(ngradcallwarm)/(cputime-cputimehold)
        callspersecmd=callspersec
        dprcheckt = (1.0d0/callspersec)*min(nmd*nParRep,1000)
        ldprcheckdone = .true.
        !write(*,*) 'spawnID ',spawnID,' dprcheckt = ',dprcheckt
        if(rank_f.eq.0) write(6,110) callspersec,1.0d0/callspersec
  110   format('force-call timing from warmup steps:',
     x  ' speed = ',f10.1,' calls/sec,or timing =',1pd12.3,' secs/call')

        ! Also Report the temperature of the High-T trajectory
        ! after the initial warmup:
        ek=0.0d0
        do ibpu=1,3*nmove
           it=itype(int((ibpu-1)/3)+1)
           ek=ek+(pxyz(ibpu)**2)/(amu(it)*1822.83d0)/2.0d0
        enddo
        tempnow=ek*1.0d0/(3.0d0/2.0d0*3.167d-6)/float(nmove)
        write(6,*) 'TEMPERATURE: t= ',time,', K= ',ek
     +      ,', V= ',e,', T= ',tempnow

      end if
      call flushtad(6)

      endif
c ======== end of iblock.eq.0 if-block for diag,dimer and warming protection ========


c === zero the time for this basin ===

      time=timehighnow   ! xxx OK now?
      thyper=0.0d0 ! not used, but unsafe to not set it

c ============= superblock to integrate trajectory for a block of steps =============

c ======== transition-rejecting blackout stage ========

      ! If we are doing ParRep... Don't bother blacking out on Master..
      if((nprocs_l.gt.nforcecores).and.(rank_l.lt.nforcecores)
     x                                         .and.(iblock.gt.0))then
        ineedblack = 0
        ibadtransyes = 1
      endif

      ngradcallhold=ngradcall
!      wctime1 = mpi_wtime()-cput_ref
      if(ineedblack.gt.0.and.ibadtransyes.eq.0) then
        timehold=time   ! xxx time
        call vecmov(xyz,xyz3,3*natom)
        ibadtrans=0
 210    continue  ! jump back to here if transition detected
        if((ietype.ge.100).and.(ietype.lt.200).and.doblocks)then
            call mdblock100(natom,nmove,xyz,pxyz,itype,
     x              maxtyp,e,grad,taxes,ietype,amu,temphigh,dt,
     x              thermrate,ineedblack,iseed,ntype)
            ngradcall = ngradcall + ineedblack
        else
            blocking_on = .true.
            do i=1,ineedblack
                call mdstep(natom,nmove,xyz,itype,pxyz,taxes,
     x              ietype, maxtyp,rcut,amu,
     x              maxtherm,itherm,jtherm,ttherm,ctherm,
     x              xtherm,mtherm,ktherm,ntherm,
     x              moving,mx,my,mz,dt,time,thyper,hyperratio,
     x              e,grad)   ! xxx time
            enddo
            blocking_on = .false.
        endif

        !! RJZ - shiftTest
        if (natom.eq.nmove) then
           call shiftalign(xyz3,xyz,natom)
        endif

       ! Check if we have been killed
       if(rank_f.eq.0) call lp_check_kill()
       call MPI_BCAST(spawnwrite,1,MPI_LOGICAL,0,force_comm,ier)
       call MPI_BCAST(activated,1,MPI_LOGICAL,0,force_comm,ier)
       if(.not.activated) goto 611 !CRZ

c        write(*,*) 'I am here after blackout!'
        call vecmov(xyz,xyz2,3*natom)
        imode=0 ! do not do spring quench
        call descent_check(natom,nmove,1+nneigh,xyz2,xyzneighs,
     +      itype,taxes,ietype,maxtyp,rcut,hyperratio,barrevev,ereverse
     +      ,imatch,irotalign,nsteep2,gfac,transcrit,itranscrit,gcrit
     +      ,dvcrit,drcrit,ipr,imode,ifundc)
        if (ereverse.gt.0.0d0.and.imatch.gt.1) then ! bpu: added 09/09/2004
           if (barrevev(imatch-1).lt.ereverse) then
              imatch=1
              if((ipr.ge.6).and.(rank_f.eq.0)) write(6,*)
     +            'reverse debug: black detected ignored state'
           endif
        endif
        if (imatch.ne.1) then
           ibadtrans=ibadtrans+1
           if(ibadtrans.lt.(maxblackreject-1)) then
             ineedblack = nblack
             if((ipr.ge.3).and.(rank_f.eq.0))
     x        write(6,*) 'SpawnID',SpawnID,
     x                   'Waiting: Transition during blackout, redoing!'
             call vecmov(xyz3,xyz,3*natom)
             call gaussp(natom,nmove,pxyz,itype,maxtyp,amu,ntype,
     x           moving,mx,my,mz,temphigh,iop,iwrt)
             call MPI_BCAST(pxyz,3*nmove,MPI_REAL8,0,force_comm,ier)
             call mdreset(0)
             filnam=trim(filstrtLP)//casename(1:lcase)//'.'//
     +           chri(istate)//'.blacktrans.'//chri(iblock)//'.'//
     +           chri(ibadtrans)//'.dat'
             if (ipr.ge.10)
     +           call storefile(natom,xyz2,pxyz,itype,taxes,filnam)
             goto 210
           elseif(ibadtrans.lt.maxblackreject) then
             ineedblack = nwarm+nthermalize
             if((ipr.ge.3).and.(rank_f.eq.0))
     x        write(6,*) 'SpawnID',SpawnID,
     x                   'Waiting: Transition during blackout, redoing!'
             call vecmov(xyz1,xyz,3*natom)
             call gaussp(natom,nmove,pxyz,itype,maxtyp,amu,ntype,
     x           moving,mx,my,mz,temphigh,iop,iwrt)
             call MPI_BCAST(pxyz,3*nmove,MPI_REAL8,0,force_comm,ier)
             call mdreset(0)
             filnam=trim(filstrtLP)//casename(1:lcase)//'.'//
     +           chri(istate)//'.blacktrans.'//chri(iblock)//'.'//
     +           chri(ibadtrans)//'.dat'
             if (ipr.ge.10)
     +           call storefile(natom,xyz2,pxyz,itype,taxes,filnam)
             goto 210
           else
              if (ibadignore.eq.0) then
             write(6,*) 'SpawnID',SpawnID,
     x                  'Waiting: Transition during blackout, aborting!'
             is=istate
             filnam=trim(filstrtLP)//casename(1:lcase)//'.'//chri(is)//
     +           '.badblack.dat'
             call storefile(natom,xyz3,pxyz,itype,taxes,filnam)
             filnam=trim(filstrtLP)//casename(1:lcase)//'.'//chri(is)//
     +           '.badblackref.dat'
             call storefile(natom,xyzneighs(1),pxyz,itype,taxes,filnam)
             filnam=trim(filstrtLP)//casename(1:lcase)//'.'//chri(is)//
     +           '.badblackhot.dat'
             call storefile(natom,xyz,pxyz,itype,taxes,filnam)
             filnam=trim(filstrtLP)//casename(1:lcase)//'.'//chri(is)//
     +           '.badblackmin.dat'
             call storefile(natom,xyz2,pxyz,itype,taxes,filnam)

             write(6,*) 'abort-too many blackouts rejects',ibadtrans
             write(6,*) 'see ...badblack.dat file'
             write(6,*)
     +           'you can set ibadignore=1 to continue in spite of this'
             stop 'too many blackout transitions found'
             else if (ibadignore.eq.1) then
                if(rank_f.eq.0) write(6,*) 'WARNING: Did '
     +              ,ibadtrans
     +              ,' blackout rejects, continuing anyway!'
     +              ,'- SpawnID ',SpawnID
                cleanblackout = .false.
                call vecmov(xyz3,xyz,3*natom)
                ibadtransyes=1 ! so we can skip the leak check
             endif
           end if
        endif
        time=timehold   ! xxx time
        call MPI_BARRIER(force_comm,ier)
        if((ipr.ge.4).and.(rank_f.eq.0))then
          write(6,*) 'spawnID ',spawnID,' finished blackout'
        endif
      end if
      ngradcallblack=ngradcallblack+ngradcall-ngradcallhold

      call flushtad(6)
c------------end of transition-rejecting blackout stage----------------

      if (ibadtransyes.eq.0.and.inojump.eq.0.and.ifundc.ne.5) then
         ngradcallhold=ngradcall
c  store a temporary copy of the beginning of this block
         call vecmov(xyz,xyz4,3*natom) ! xxxx may want before if or store a copy above for safety

c ===== catastrophic check for leakage out of the basin =====

         call vecmov(xyz,xyz2,3*natom)
         imode=0                ! do not do spring quench
         call descent_check(natom,nmove,1+nneigh,xyz2,xyzneighs, ! xxxx these grads not in summary
     +       itype,taxes,ietype,maxtyp,rcut,hyperratio,barrevev,ereverse
     +       ,imatch,irotalign,nsteep2,gfac,transcrit,itranscrit,gcrit
     +       ,dvcrit,drcrit,ipr,imode,ifundc)

         if (ereverse.gt.0.0d0.and.imatch.gt.1) then
            if (barrevev(imatch-1).lt.ereverse) then
               imatch=1
               if (ipr.ge.6) write(6,*)
     +             'reverse debug: leak detected ignored state'
            endif
         endif
         if(imatch.ne.1) then
            filnam=trim(filstrtLP)//casename(1:lcase)//'.leak.dat'
            call storefile(natom,xyz2,pxyz,itype,taxes,filnam)
            filnam=trim(filstrtLP)//casename(1:lcase)//'.xyz.dat'
            call storefile(natom,xyz,pxyz,itype,taxes,filnam)
            filnam=trim(filstrtLP)//casename(1:lcase)//'.xyz1.dat'
            call storefile(natom,xyz1,pxyz,itype,taxes,filnam)
            filnam=trim(filstrtLP)//casename(1:lcase)//'.xyz2.dat'
            call storefile(natom,xyz2,pxyz,itype,taxes,filnam)
            filnam=trim(filstrtLP)//casename(1:lcase)//'.xyz3.dat'
            call storefile(natom,xyz3,pxyz,itype,taxes,filnam)
            filnam=trim(filstrtLP)//casename(1:lcase)//'.xyz4.dat'
            call storefile(natom,xyz4,pxyz,itype,taxes,filnam)
            write(6,*) 'STOP: abort-trajectory leaked out somehow'
            write(6,*) 'STOP: see <casename>.leak.dat file'
            write(6,*) 'spawnID ',spawnID,' (groupID ',groupID,'):'
     +       ,' trajectory has leaked out of istate - ignoring!'
            !stop 'bad news: trajectory has leaked out of istate'
         end if
         ngradcallleak=ngradcallleak+ngradcall-ngradcallhold
      endif

c ======== do the real MD for one block ========

      ihold=0
      iblock=iblock+1
      ngradcallhold=ngradcall
      ineberr=0                 ! reset neb error code
      ibadtransyes=0            ! reset if had trans during blackout or warming

      !!RJZ-HERE
      !if(iblock.eq.1) wctime1 = mpi_wtime()-cput_ref

      if(nprocs_l.eq.nforcecores)then

!        if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
!     x                 ,' entering MD block.. '
!        if(rank_f.eq.0) flush(6)
        call vecmov(xyz,xyz2,3*natom)

        ! Do nhold loops if using lammps mdblock routine:
        nmd_use = nmd
        if((ietype.ge.100).and.(ietype.lt.200).and.doblocks)then
            nmd_use = nhold !nmd / nmdsub
        endif
        do 160 i=1,nmd_use
          if((ietype.ge.100).and.(ietype.lt.200).and.doblocks)then
             call mdblock100(natom,nmove,xyz,pxyz,itype,
     x              maxtyp,e,grad,taxes,ietype,amu,temphigh,dt,
     x              thermrate,nmdsub,iseed,ntype)
             ngradcall = ngradcall + nmdsub
             time = time + dt * nmdsub
          else
            blocking_on = .true.
            call mdstep(natom,nmove,xyz,itype,pxyz, taxes,
     x       ietype, maxtyp,rcut,amu,
     x       maxtherm,itherm,jtherm,ttherm,ctherm,xtherm,mtherm,ktherm ! xxx time
     x       ,ntherm,moving,mx,my,mz,dt,time,thyper,hyperratio,e,grad)
            blocking_on = .false.
          endif
          ek=0.0d0
          do ibpu=1,3*nmove
           it=itype(int((ibpu-1)/3)+1)
c           write(6,*) 'blah: ',ibpu,it
           ek=ek+(pxyz(ibpu)**2)/(amu(it)*1822.83d0)/2.0d0
          enddo
          tempnow=ek*1.0d0/(3.0d0/2.0d0*3.167d-6)/float(nmove)
          if (ipr.ge.20) write(6,*) 'debug: t=',time,', K=',ek
     +      ,', V=',e,', T=',tempnow

          !rjz! write(6,*) 'debug: t=',time,', K=',ek,', V=',e,', T=',tempnow

c         hold regular snapshots for refining transition point
!          if((mod(i,nmdsub).eq.0).or.
!     +                 ((ietype.ge.100).and.(ietype.lt.200))) then
          if(mod(i,nmdsub).eq.0) then
             ihold=ihold+1
             !! RJZ - shiftTest
             if (natom.eq.nmove) then
                 call shiftalign(xyz2,xyz,natom)
             endif
             call vecmov(xyz,xyzhold((ihold-1)*3*natom+1),3*natom)
          endif
  160   continue
        !! RJZ - shiftTest
        if (natom.eq.nmove) then
           call shiftalign(xyz2,xyz,natom)
        endif
        ngradcallnmd=ngradcallnmd+ngradcall-ngradcallhold
c        write(6,*) 'spawnID',spawnID,'debug: t=',time,', K=',ek
c     +      ,', V=',e,', T=',tempnow

      endif


      ! Wait for ParRep transitions/time
      if(nprocs_l.gt.nforcecores)then
        iprsafe = .true.
        call vecmov(xyz1,xyz,3*natom) ! Make sure we are in minimum (so a timecheck wont find a transition)
!        if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
!     x                 ,' master entering ParRep loop '
!        if(rank_f.eq.0) flush(6)
        if(rank_f.eq.0) wctrz1 = mpi_wtime()-cput_ref
        prtranstime = 0.d0
        lprtrans = .false. ! Transition or time-check
        lprtime  = .false.
        do while ((.not.lprtrans).and.(.not.lprtime))
          transwct = -1.d0
          prtranstime = 0.d0
          lprtrans = .false.
          lprtime  = .false.
          iprsafe = .true.
          call vecmov(xyz1,xyz,3*natom)
          call pr_check_trans(nmove,natom,ihold,xyz,xyzhold
     +        ,nholdmaxworstcase,prtranstime,lprtrans
     +        ,iprsafe,transwct,wctimehighstop,lrepstold)
          !if(lprtrans.and.(spawnIDlocal.eq.0))then
          if(lprtrans)then
              if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
     x                 ,' new transition at time: ',prtranstime
     x                 ,' and wct: ',transwct
              n_mas_trans = n_mas_trans + 1
          endif
          if((parallelmode.eq.'spectad')
     x                 .and.(irecognize.ge.1).and.(.not.lprtrans))then
            call lp_update_official(nofficialmax)
            call MPI_BARRIER(force_comm,ier)
          endif
          if(rank_f.eq.0) call lp_check_kill()
          call MPI_BCAST(spawnwrite,1,MPI_LOGICAL,0,force_comm,ier)
          call MPI_BCAST(activated,1,MPI_LOGICAL,0,force_comm,ier)
          call MPI_BCAST(transwct,1,MPI_REAL8,0,force_comm,ier)
          if(.not.activated) goto 611 !CRZ
          ! Every so often get ParRep times
          if(.not.lprtrans)then
            if(rank_f.eq.0) wctrz2 = mpi_wtime()-cput_ref
            if(rank_f.eq.0) tlapsed = (wctrz2-wctrz1)
            call MPI_BCAST(tlapsed,1,MPI_REAL8,0,force_comm,ier)
            if((tlapsed.ge.dprcheckt).or.(ldepspawn.eq.1))then
!              if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
!     x                 ,' master checking for time...'
!              if(rank_f.eq.0) flush(6)
              lprtime = .true.
              call pr_check_time(prtranstime,iprsafe,wctimehighstop)
              if(.not.iprsafe)then ! Don't bother, lets wait for the trans message
                if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
     x                 ,' master bailing on time-check...'
                if(rank_f.eq.0) flush(6)
                lprtime = .false.
                prtranstime = 0.d0
              else
!                if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
!     x                 ,' master done with time-check...'
!                if(rank_f.eq.0) flush(6)
              endif
            endif
          endif
          ! Set wctimehighstop:
          if((timehighnowI+prtranstime).ge.timehighstop)then
            if(transwct.gt.0.d0)then
              if(transwct.lt.wctimehighstop) wctimehighstop = transwct
              if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
     +           ,' setting wctimehighstop to ',wctimehighstop
              if(rank_f.eq.0) flush(6)
            endif
          endif
          ! Don't deal with a transition if it happened after the stoptime:
          if(lprtrans.and.
     x                ((timehighnowI+prtranstime).ge.timehighstop))then
            ! Unless it will end TAD in this state... then let it through as a time-check:
            if(iprsafe)then
              lprtime = .true.
              lprtrans = .false.
              call vecmov(xyz1,xyz,3*natom)
            else
              if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
     x           ,' ignoring post-stop time transition at time: '
     x           ,prtranstime,' and wct: ',transwct
              if(rank_f.eq.0) flush(6)
              lprtrans = .false.
              call MPI_BARRIER(force_comm,ier)
            endif
          endif
        enddo
        time = timehighnowI + prtranstime

        t_mas_tot = t_mas_tot_i + prtranstime
        if(lprtrans)then
            transtimerep = t_mas_tot / n_mas_trans
            if(rank_f.eq.0)
     +         write(6,*) 'spawnID ',spawnID,' rank_l ',rank_l
     +                   ,' has AVG trans time:',transtimerep
     +                   ,'- n_mas_trans:',n_mas_trans
     +                   ,'- t_mas_tot:',t_mas_tot
            if(rank_f.eq.0) call flushtad(6)
        endif

        ! Make sure a transition will not be accidentally detected
        if(.not.lprtrans)then
          call vecmov(xyz1,xyz,3*natom)
        endif
      endif

      ! Check if we have been killed
      if(rank_f.eq.0) call lp_check_kill()
      call MPI_BCAST(spawnwrite,1,MPI_LOGICAL,0,force_comm,ier)
      call MPI_BCAST(activated,1,MPI_LOGICAL,0,force_comm,ier)
      if(.not.activated) goto 611 !CRZ

c      if(iblock.eq.112) then !cwp for LBFGS debugging!
c          open(unit=iunit,file='LBFGSdebug.dat',form='formatted',
c     x              status='unknown')
c      rewind iunit
c      iwrt=-1
c      ntitle=1
c      ctitle(1)='cluster file from TAD2; meaning encoded in file name'
c      iop=1  ! 1 gives a compact, no-momentum file   - new as of 9/24/02
c      call putclsnew(natom,ntitle,ctitle,xyz,pxyz,itype,
c     x                                 taxes,iunit,iwrt,iop)

c      endif

c ================ check for a transition ================

      ngradcallhold=ngradcall
      if((ipr.ge.4).and.(rank_f.eq.0))
     +  write(6,*) 'spID ',spawnID,' do a transition check',ngradcall
      call vecmov(xyz,xyz2,3*natom)
      imode=0                   ! do not do spring quench
      call descent_check(natom,nmove,1+nneigh,xyz2,xyzneighs,
     +    itype,taxes,ietype,maxtyp,rcut,hyperratio,barrevev,ereverse
     +    ,isame,irotalign,nsteep2,gfac,transcrit,itranscrit,gcrit
     +    ,dvcrit,drcrit,ipr,imode,ifundc)
      if(isame.eq.0) then
        ijump=1
        ineigh=0  ! neighbor number - zero if not recognized
        newneighbor=-1
        inojump=0
      else if(isame.ge.2) then
        ijump=1
        ineigh=isame-1  ! neighbor number - zero if not recognized
        newneighbor=0
        inojump=0
      else
        ijump=0
        ineigh=0  ! not needed if ijump=0, but less confusing    xxx - youch - does this goof up the logic below????
        newneighbor=0
        inojump=1               ! skip leak check
      end if
      call MPI_BARRIER(force_comm,ier)
      if(rank_f.eq.0)then
        !write(*,*) 'spID ',spawnID,' did MD to thigh = ',time
        if(ipr.ge.4) write(6,*) 'done with transition check',ngradcall
        if(ipr.ge.4) write(6,*) 'ijump =',ijump,' ineigh =',ineigh
      endif
      ngradcalltrans=ngradcalltrans+ngradcall-ngradcallhold

c =========== at this point, ijump is defined, but it may not be a real transition ============

c update the total number of attempts
      if(ijump.eq.1) then
        nattempt=nattempt+1
c       time=time-blocktime/2.0d0 ! subtract off out-of-state time   ! xxx time
        if(nprocs_l.gt.nforcecores)then
          ! ParRep time used here
          timehighnow=time
        else
          if(rank_f.eq.0)then
            time=time-blocktime*prngen(0) ! fix for distributing time throughout blok (BPU: 9/22/05)
          endif
          call MPI_BCAST(time,1,MPI_REAL8,0,force_comm,ier)
          timehighnow=time
        endif
        if(rank_f.eq.0)then
          if(ipr.ge.4) write(6,6665) nattempt,istate,ilabel,timehighnow
        endif
 6665   format('------- transition attempt no. ',i5,
     x      ' from official state ',i5,
     x      ' which is currently state ',i5,
     x       '; timehigh: =',1pd12.4)
        if((ipr.ge.4).and.(rank_f.eq.0)) write(6,*) 'ineigh=',ineigh
      else
        timehighnow=time
      end if
      call flushtad(6)

c store the after-transition minimized structure ("attempt") before refinement
c note that this may include attempts that later turn out to be no transition

      if(ijump.eq.1) then
        if(lstore_attempt) then
          filnam=trim(filstrt)//casename(1:lcase)//
     x         '.'//chri(ilabel)//'.attempt.'//chri(nattempt)//'.dat'
             call storefile(natom,xyz2,pxyz,itype,taxes,filnam)
        end if
c store chs of attempt for analysis xxx replaced with call
        if (lstore_attempt_chs) then
           filnam=trim(filstrtLP)//casename(1:lcase)//
     x         '.'//chri(ilabel)//'.attempt.'//chri(nattempt)//'.chs'
           call storechs(filnam,nhold,natom,taxes,springk,itype,xyzhold)
        endif
      end if

c ========= try a simple energy and distance match if irecognize=9 ==========

      if(ijump.eq.1 .and. ineigh.eq.0 .and. irecognize.eq.9) then
         call gcalc(natom,nmove,xyz2,itype, taxes, ! get energy of endpoint
     x       ietype, maxtyp,rcut, energy2,grad, hyperratio)

        call eneighmatch(natom,nmove,xyz2,energy2,taxes,eneigh,
     x       xyzneighs, transcrit, ninvolved,
     x        nneigh,erecognizecrit,rrecognizecrit,
     x                             r4recognizecrit,ipr,
     x       ineigh)
        if(ineigh.ne.0) then
          if(ipr.ge.2) write(6,*) 'irec9 energy match for neigh:',ineigh
        end if
      end if

c ==============  begin unrecognized-neighbor block ============

      if(ijump.eq.1 .and. ineigh.eq.0) then

         ngradcallhold=ngradcall
c  refine the transition point, finding the first out-of-basin geometry
         imode=0 ! do not do spring quench
         call vecmov(xyz2,xyz5,3*natom)
         call refine_transition(natom,nmove,nhold,xyz1,xyz2,xyzhold
     +       ,xyzneighs,1+nneigh,barrevev,ereverse,itype,taxes,ietype
     +       ,maxtyp,rcut,hyperratio,irotalign,nsteep1,gfac,transcrit
     +       ,itranscrit,intermediate,gcrit,dvcrit,drcrit,ipr,imode
     +       ,ileft,iright,ifundc,rank_f)
        if (intermediate.ne.0) then
         if(rank_f.eq.0)then
          if(ipr.ge.4)write(6,*)'TAD: Replacing x2 with refined minimum'
         endif
         filnam=trim(filstrtLP)//casename(1:lcase)//
     x        '.'//chri(ilabel)//'.refined.'//chri(nattempt)//'.dat'
         if (ipr.ge.10) call storefile(natom,xyz2,pxyz,itype,taxes
     +        ,filnam)
        endif
        call flushtad(6)

c get intermediate hot trajectory point to mix into NEB
        hotmixuse=hotmix
      if (ijump.eq.1.and.ineigh.eq.0.and.hotmix.gt.0.0d0) then
         if (irotalign.ne.0) write(6,*)
     +       'WARNING! May not work correctly with free systems!'
         do i=1,3*natom
            xyzhot(i)=xyzhold((ileft-1)*3*natom+i)/2.0d0+xyzhold((iright
     +          -1)*3*natom+i)/2.0d0
         enddo
         if (ipr.ge.5) then
            filnam=trim(filstrtLP)//casename(1:lcase)//
     x          '.'//chri(istate)//'.xhot.'//chri(nneigh)//'.dat'
            call storefile(natom,xyzhot,pxyz,itype,taxes,filnam)
         endif
      endif

c ========= again try a simple energy and distance match if irecognize=9 ==========

      if(ijump.eq.1 .and. ineigh.eq.0 .and. irecognize.eq.9) then
         call gcalc(natom,nmove,xyz2,itype, taxes, ! get energy of endpoint
     x       ietype, maxtyp,rcut, energy2,grad, hyperratio)

        call eneighmatch(natom,nmove,xyz2,energy2,taxes,eneigh,
     x         xyzneighs, transcrit,ninvolved,
     x         nneigh,erecognizecrit,rrecognizecrit,
     x                               r4recognizecrit,ipr,
     x         ineigh)
        if(ineigh.ne.0) then
          if(ipr.ge.2) write(6,*)'irec9 energy match2 for neigh:',ineigh
        end if
      end if


c try again to see if we recognize an existing neighbor
c note that the descent check should go very fast when already minimized

      if(intermediate.ne.0.or.irotalign.ge.1) then
         imode=0
         if (irotalign.ge.1.and.itranscrit.eq.1) imode=1
         call descent_check(natom,nmove,1+nneigh,xyz2,xyzneighs,
     +       itype,taxes,ietype,maxtyp,rcut,hyperratio,barrevev,ereverse
     +       ,isame,irotalign,nsteep2,gfac,transcrit,itranscrit,gcrit
     +       ,dvcrit,drcrit,ipr,imode,ifundc)
         if(isame.eq.1) then
            write(6,*) 'tad: confused: new neigh is state A'
            filnam=trim(filstrtLP)//casename(1:lcase)//
     x          '.'//chri(istate)//'.notnew.'//chri(nneigh)//'.dat'
            call storefile(natom,xyz2,pxyz,itype,taxes,filnam)
            if (iright+1.gt.nhold) then
               write(6,*)
     +             'do not have enough snaps to refine to right'
               write(6,*) 'pretending no jump occured'
               ijump=0
               ibadtransyes=1   ! don't do blackout or leak check as it will fail
            else
               write(6,*) 'trying iright+1 state: ',iright+1
               call vecmov(xyzhold((iright+1-1)*3*natom+1),xyz2,3
     +             *natom)
               call descent_check(natom,nmove,1+nneigh,xyz2,xyzneighs,
     +             itype,taxes,ietype,maxtyp,rcut,hyperratio,barrevev
     +             ,ereverse,isame,irotalign,nsteep2,gfac,transcrit
     +             ,itranscrit,gcrit,dvcrit,drcrit,ipr,imode,ifundc)
               if (isame.eq.1) then
                  write(6,*) 'tad: confused: iright+1 is also state A'
                  write(6,*)
     +                'moving xyz5 back to xyz2 and passing to NEB'
                  call vecmov(xyz5,xyz2,3*natom)
               endif
               hotmixuse=0.0d0
            endif
         endif
         if(isame.ne.0) then
            ineigh=isame-1
         end if
      end if
      ngradcallrefine=ngradcallrefine+ngradcall-ngradcallhold

c =========== if a transition was detected, find the saddle using NEB, and store saddle =============

c xxx find minimum for iright and ileft to pass to neb, to bracket saddle
C xxx better and minimize the amount of NEB to do
c xxx bracketing code pulled for now, until more fully decide if needed
      prefacsad = -1.d0
      if(ineigh.eq.0) then

         ! Check if we have been killed
         if(rank_f.eq.0) call lp_check_kill()
         call MPI_BCAST(spawnwrite,1,MPI_LOGICAL,0,force_comm,ier)
         call MPI_BCAST(activated,1,MPI_LOGICAL,0,force_comm,ier)
         if(.not.activated) goto 611 !CRZ

         ngradcallhold=ngradcall
         filnam1=trim(filstrtLP)//casename(1:lcase)//
     x       '.'//chri(ilabel)//'.nebi.'//chri(nneigh)//'.chs'
         filnam2=trim(filstrtLP)//casename(1:lcase)//
     x       '.'//chri(ilabel)//'.nebf.'//chri(nneigh)//'.chs'
         if (ipr.gt.4) write(6,*) 'hotmixuse: ',hotmixuse

         if(ipr.gt.1)then
           if(rank_f.eq.0) cputimeNEBi = mpi_wtime()-cput_ref
           call MPI_BCAST(cputimeNEBi,1,MPI_REAL8,0,force_comm,ier)
         endif

         call saddle_find(natom,xyz1,xyz2,xyzhot
     +       ,hotmixuse,xyzs,esad,nmove,itype,taxes,ietype,maxtyp,rcut
     +       ,amu,hyperratio,gfac,irotalign,nsteep1,ineberr
     +       ,intermediate2,transcrit,itranscrit,gcrit,dvcrit,drcrit,ipr
     +       ,nimage,eshallow,xyzneighs,nneigh+1,barrevev,ereverse,itan
     +       ,iclimb,springk,lstore_neb,filnam1,filnam2,inebdimer
     +       ,inebirc,ifundc,ifunneb,imodecheck,rmodemag
     +       ,nnegsad,nprod,freqlogssad,freqlogsmin,nnegmin,prefacsad
     +       ,ivineyard,ntype,potnam,emin)
         !if(.not.lqeqfix) call refix_qeq()

         if(ipr.gt.1)then
           if(rank_f.eq.0) cputimeNEBf = mpi_wtime()-cput_ref
           call MPI_BCAST(cputimeNEBf,1,MPI_REAL8,0,force_comm,ier)
           if(rank_f.eq.0) write(*,'(A,I9,A,ES16.6,A,ES16.6,A)')
     x        'SpawnID ',SpawnID
     x       ,' NEB TIME = ',cputimeNEBf-cputimeNEBi
     x       ,' s For Ea = ',(esad-emin)*27.21d0,' eV '
           call MPI_BARRIER(force_comm,ier)
         endif

c       note:  intermediate2=1 signals that an intermediate minimum was found
c       and that it has been passed back in x2.  If so, we will
c       want to use this (and warm it up) if we end up accepting this
c       transition.  It also means that there was no transition in the
c       chain, that one of the end points was really a saddle.  If so,
c       x2 should be made the same as x1.
         if (intermediate2.ne.0) then
            if(ipr.ge.4) write(6,*) 'TAD: Replacing x2 with NEB minimum'
            filnam=trim(filstrtLP)//casename(1:lcase)//
     x          '.'//chri(ilabel)//'.replaced.'//chri(nattempt)//'.dat'
            if (ipr.ge.4) call storefile(natom,xyz2,pxyz,itype,taxes
     +          ,filnam)
         endif
         if (ineberr.eq.2) then
            if((ipr.ge.0).and.(rank_f.eq.0)) write(6,*)
     +          'TAD: saddle_find returned without finding a saddle'
            ibadtransyes=1      ! xxx don't do blackout or leak check as it will fail
            ijump=0
         endif
         if((ijump.eq.1).and.((esad-emin).le.(0.d0)))then
            if(rank_f.eq.0) write(6,*)
     x          'spawnID,intermediate2:  ',spawnID,intermediate2,
     x          ' - ignoring transition with NEGATIVE barrier !!',
     x          ' - nnegmin = ',nnegmin,
     x          ' - emin, esad = ',emin, esad
            ibadtransyes=1      ! xxx don't do blackout or leak check as it will fail
            ijump=0
            filnam=trim(filstrtLP)//casename(1:lcase)//
     x          '.negmin.'//chri(spawnID)//'.dat'
            call storefile(natom,xyz1,pxyz,itype,taxes,filnam)
            filnam=trim(filstrtLP)//casename(1:lcase)//
     x          '.negsad.'//chri(spawnID)//'.dat'
            call storefile(natom,xyzs,pxyz,itype,taxes,filnam)
            filnam=trim(filstrtLP)//casename(1:lcase)//
     x          '.negend.'//chri(spawnID)//'.dat'
            call storefile(natom,xyz2,pxyz,itype,taxes,filnam)
            if(abs(esad-emin).ge.(0.001d0))
     x          STOP ' ERROR -- transition with NEGATIVE barrier '
         endif
      else
         ineberr=0
      end if
      call flushtad(6)

c try one last time to recognize this neighbor (since saddle search may have changed it)
c if this is a new neighbor, put the end state into the xyzneighs array.

      if(ineigh.eq.0 .and. intermediate2.eq.1) then
         imode=0 ! do not do spring quench
         ifundcuse=ifundc
         if (ifundcuse.eq.5) ifundcuse=0
        call descent_check(natom,nmove,1+nneigh,xyz2,xyzneighs,
     +       itype,taxes,ietype,maxtyp,rcut,hyperratio,barrevev,ereverse
     +       ,isame,irotalign,nsteep2,gfac,transcrit,itranscrit,gcrit
     +       ,dvcrit,drcrit,ipr,imode,ifundcuse)
        if(isame.eq.1) then
          if(ipr.ge.2) write(6,*) 'interesting - final state is initial'
          if(ipr.ge.2) write(6,*) 'resetting ijump to zero in midstream'
          filnam=trim(filstrtLP)//casename(1:lcase)//
     x        '.'//chri(ilabel)//'.same2.'//chri(nattempt)//'.dat'
          if (ipr.ge.10) call storefile(natom,xyz2,pxyz,itype,taxes
     +        ,filnam)
          filnam=trim(filstrtLP)//casename(1:lcase)//
     x        '.'//chri(ilabel)//'.same1.'//chri(nattempt)//'.dat'
          if (ipr.ge.10) call storefile(natom,xyzneighs(1),pxyz,itype
     +        ,taxes,filnam)
          ineigh=1
          ijump=0
          call vecmov(xyz4,xyz,3*natom) ! move hot point at beginning of block into xyz
        end if
        if(isame.ge.2) then
          ineigh=isame-1
        end if
      end if
      ngradcallneb=ngradcallneb+ngradcall-ngradcallhold
      call flushtad(6)

c ========= and again try a simple energy and distance match if irecognize=9 ==========

      if(ijump.eq.1 .and. ineigh.eq.0 .and. irecognize.eq.9) then
c     x          .and. intermediate2.eq.1) then ! bpu: commented out 09/09/2004 to force another check, just in case
        call gcalc(natom,nmove,xyz2,itype, taxes, ! get energy of endpoint
     x       ietype, maxtyp,rcut, energy2,grad, hyperratio)

        call eneighmatch(natom,nmove,xyz2,energy2,taxes,eneigh,
     x         xyzneighs, transcrit,ninvolved,
     x         nneigh,erecognizecrit,rrecognizecrit,
     x                               r4recognizecrit,ipr,
     x         ineigh)
        if(ineigh.ne.0) then
          if(ipr.ge.2) write(6,*)'irec9 energy match3 for neigh:',ineigh
        end if
      end if

c ======= OK - now we are sure whether or not a transition occurred ===========

c ======== If transition occurred,  add this unrecognized new neighbor to the list =======

      if(ijump.eq.1 .and. ineigh.eq.0) then    ! now need to be checking ijump, since it might have been zeroed

        ! CHECK if transcrit is small enough:
        call maxdr_rz(nmove,xyz1,xyz2,deltarmax)
        if (deltarmax.lt.transcrit) then ! potential match
          write(*,*) 'ERROR -- SpawnID',SpawnID,
     x               'detected that transcrit is too low!!!',
     x               '- transition has deltarmax =',deltarmax
          stop 'STOP -- transcrit too low'
        endif

        ! get energy of endpoint xxxx unaccounted for in ngracall summary
        call gcalc(natom,nmove,xyz2,itype, taxes,
     x       ietype, maxtyp,rcut, energy2,grad, hyperratio)

        nneigh=nneigh+1
        ineigh=nneigh
        newneighbor=1
        if (nneigh.gt.nneighmax) then
           write(6,*) 'STOP: nneigh.gt.nneighmax, too many neighs'
           stop 'nneigh.gt.nneighmax, too many neighs'
        endif
        esaddle(ineigh)=esad
!        write(*,*) 'spawnID ',spawnID,' FINDING emin,esad = '
!     +      ,emin,esad,' for ineigh ',ineigh,' in istate ',iofficial
        eneigh(nneigh)=energy2
        barrierev(ineigh)=(esad-emin)*27.21d0
        barrevev(ineigh)=(esad-energy2)*27.21d0
        ! RJZ - Hack to raise low barriers (gives wong dynamics, but allows faster exploration of system)
        if(barrierev(ineigh).lt.(minbarrier_hack))then
          barrierev(ineigh) = minbarrier_hack
        endif
        if(rank_f.eq.0)then
          if(ipr.ge.6) write(6,*) 'reverse debug: reverse barrier: '
     +      ,ineigh,barrevev(ineigh)
        endif
        if(3*natom*(1+nneigh).gt.lenxyzneighs) then
           write(6,*)
     +      'STOP: nneigh*3*natom.gt.lenxyzneighs, too many neighs'
           stop
     +      'nneigh*3*natom.gt.lenxyzneighs, too many neighs'
        endif
        call vecmov(xyz2,xyzneighs(3*natom*nneigh+1),3*natom)
        call vecmov(xyzs,xyzneighsads(3*natom*nneigh+1),3*natom)
        call dmzer(rdimermodes(3*natom*nneigh+1),3*natom)

c       store a cluster file of the saddle point
        if(lstore_sad) then
          filnam=trim(filstrt)//casename(1:lcase)//
     x         '.'//chri(ilabel)//'.sad.'//chri(ineigh)//'.dat'
             call storefile(natom,xyzs,pxyz,itype,taxes,filnam)
        end if

c ===  store a cluster file of the end (product) state ===
        if(lstore_end) then
           filnam=trim(filstrt)//casename(1:lcase)//
     x         '.'//chri(ilabel)//'.end.'//chri(ineigh)//'.dat'
           call storefile(natom,xyz2,pxyz,itype,taxes,filnam)
        end if

c   quick computation of hyperdistance and 4-dist just for later write to barriers file
        r2x=0.0d0
        do j=1,3*nmove
          r2x = r2x + (xyz1(j)-xyz2(j))**2
        end do
        hyperdist=sqrt(r2x)
c  sum up 4-distance
        r4=0.0d0
        jjj=0
        do j=1,nmove
        sumx=0.0d0
        do jj=1,3
        jjj=jjj+1
        sumx=sumx + (xyz1(jjj)-xyz2(jjj))**2
        end do
        r4=r4+sumx**2
        end do
        dist4=sqrt(sqrt(r4))

      else
        newneighbor=0
      end if

      end if   !----------------- end of unrecognized-neighbor block
      call flushtad(6)

c at this point, we have ijump, ineigh and nneigh set correctly

c =========== Find the reverse barrier.  If lower than ereverse, pretend no transition ============
c xxx might want to store states to ignore so don't have to refind them all the time
      if ((ijump.eq.1.and.ipr.ge.6).and.(rank_f.eq.0)) write(6,*)
     +    'reverse debug: nn,erev,barrevev: ',newneighbor,ereverse
     +    ,barrevev(ineigh)
      if (ereverse.gt.0.0d0.and.ijump.eq.1) then
         if (barrevev(ineigh).lt.ereverse) then
            ibadtransyes=1      ! don't do blackout or leak check
            ijump=0
            if((ipr.ge.2).and.(rank_f.eq.0))
     +        write(6,*)'TAD2: reverse barrier small ('
     +          ,barrevev(ineigh),' eV), ignoring'
            if (newneighbor.eq.1) then
               filnam=trim(filstrtLP)//casename(1:lcase)//'.barriers'
               ctype='MD'
c write barrier entry anyways so we have the info handy (added 29Mar04)
               call numberinvolved(xyz1,xyz2,nmove,transcrit,
     x                                     ninvolved(ineigh))
               call writebarrierentry(filnam,timehighnow
     +             ,barrierev(ineigh),betalow,betahigh,istate,iofficial
     +             ,ineigh,prefac(ineigh),ninvolved(ineigh)
     +             ,barrevev(ineigh),ctype,hyperdist,dist4)
            endif
            newneighbor=0
         endif
      endif
      call flushtad(6)

      ! For now - add blocking info for first 10 neighbors seen
      ! -> (This should really help with warmups/blackouts after many revisits)
      if((ijump.eq.1).and.(allow_blocking).and.(ineigh.le.10))then
        if((newneighbor.eq.1).or.
     +    (statedata(iofficial)%ibondneigh1(ineigh).eq.0)) then
          !! Record the most-stretched bond at the saddle for "blocking" purposes:
          ibmax = 0
          jbmax = 0
          kbmax = 0
          dbmax12 = 0.d0
          dbmax13 = 0.d0
          dbmax12i = 1.d0
          dbmax13i = 1.d0
          !call maxdr_rz_2(nmove,xyz1,xyzs,dbmax,ibmax,jbmax)
          call maxdr_rz_2
     +         (nmove,xyz1,xyzs,dbmax12,dbmax13,dbmax12i,
     +                               dbmax13i,ibmax,jbmax,kbmax)
          if((ibmax.gt.0).and.(jbmax.gt.0).and.(kbmax.gt.0))then
            statedata(iofficial)%nneigh_block =
     +                 statedata(iofficial)%nneigh_block + 1
            statedata(iofficial)%ibondneigh1(ineigh) = ibmax
            statedata(iofficial)%ibondneigh2(ineigh) = jbmax
            statedata(iofficial)%ibondneigh3(ineigh) = kbmax
            statedata(iofficial)%dbondneigh12(ineigh)=dbmax12
            statedata(iofficial)%dbondneigh13(ineigh)=dbmax13
            statedata(iofficial)%dbondneigh12i(ineigh)=dbmax12i
            statedata(iofficial)%dbondneigh13i(ineigh)=dbmax13i
            if(rank_f.eq.0)then
              write(*,*) 'SpawnID',SpawnID,' - for ineigh',ineigh
     +           ,'- ibmax:',ibmax,'jbmax:',jbmax,'kbmax:',kbmax
     +           ,'- dbmax12:',dbmax12,'ibmax13:',dbmax13
     +           ,'- dbmax12i:',dbmax12i,'ibmax13i:',dbmax13i
              call flushtad(6)
            endif
          endif
        endif
      endif


c ============== begin block of work for newly discovered transitions ============

      if(ijump.eq.1 .and. newneighbor.eq.1) then   !------+------+------+ new-jump block

        ! this indicates NEB found this one, rather than dimer:
        nebfound(ineigh)=1

c determine how many atoms were involved in the transition
c xxx this only works for idetect=0, for others is meaningless
c xxxxa after we are sure nebfound printout is working right, put this call into the dimer region also
        call numberinvolved(xyz1,xyz2,nmove,transcrit,ninvolved(ineigh))
        if(rank_f.eq.0)then
         if(ipr.ge.4) write(6,*) 'this transition out of state',istate,
     x    ' to neighbor ',ineigh,' involves ',ninvolved(ineigh),' atoms'
        endif

c compute vineyard prefactor for this new escape path
c xxx relegated to subroutine
        ngradcallhold=ngradcall
        if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
        call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
        cputimehold=cputime
        ! If spectad do this part of the prefactor after the spawn:
        if((ivineyard.gt.0).and.(parallelmode.ne.'spectad')) then
          if (ivineyard.eq.1) then
             igeometry=1        ! saddle
             freqlogssad=0.d0
             call dovineyard(natom,nmove,xyzs,itype,taxes,ietype,maxtyp
     +           ,rcut,amu,e,grad,nnegsad,nprod,freqlogssad,igeometry
     +           ,freqlogsmin,nnegmin,prefac(ineigh),ierr,ipr)
c xxx should we print warning if prefac < prefacmin?
             if(ierr.eq.0) then
                if(ipr.ge.2) write(6,*) 'Vineyard prefac for state '
     +              ,istate,' to neigh(i.e., saddle) ',ineigh,'  ='
     +              ,prefac(ineigh),' Hz'
             end if
          endif
c xxxx new subset vineyard
          if (ivineyard.eq.2) then
             iupper=4
             if (idynmatloop.eq.0) iupper=0
             do ibpu=iupper,0,-1
                rdynmatuse=rdynmatcrit/(3**float(ibpu))
                call vineyardsubset(natom,nmove,xyz1,xyzs,itype
     +              ,taxes,ietype,maxtyp,rcut,amu,rdynmatuse,
     +              e,grad,prefacx,nelements,ierr,ipr)
                if (ibpu.eq.iupper) prefacxx=prefacx
                write(6,6543) istate,nneigh,rdynmatuse
     +              ,prefacx,nelements,prefacxx
             enddo
             prefac(ineigh)=prefacxx
          endif

        end if
        ngradcallvineyard=ngradcallvineyard+ngradcall-ngradcallhold
        if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
        call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
        cputimevineyard=cputimevineyard+cputime-cputimehold
        call flushtad(6)

c write entry to barriers file xxx consolidated to a call
        filnam=trim(filstrtLP)//casename(1:lcase)//'.barriers'
        ctype='MD'
        call writebarrierentry(filnam,timehighnow,barrierev(ineigh)
     +      ,betalow,betahigh,istate,iofficial,ineigh,prefac(ineigh)
     +      ,ninvolved(ineigh),barrevev(ineigh),ctype,hyperdist,dist4)

c see if this becomes the lowest barrier xxx moved to loop below

c check whether dimer barrier still holds up as lowest barrier xxx moved
C to loop below

      end if                    !------+------+------+  end of new-jump block
      call flushtad(6)



      if(ijump.eq.1) then
        if(ineigh.le.0) stop 'code screwed up 556'
        nattempted(ineigh)=nattempted(ineigh)+1
        nnewattempted(ineigh)=nnewattempted(ineigh)+1 !! RJZ
        if((ipr.ge.4).and.(rank_f.eq.0))
     x   write(6,*) 'attempted transition was to neighbor',ineigh
        if(ineberr.ge.0) then
          if((ipr.ge.4).and.(rank_f.eq.0))
     x      write(6,*) 'note: saddle_find returned ineberr=',ineberr
        end if
        if(ineberr.gt.5) stop 'crashing on NEB error code'
      end if

c ===== compute the extraplolated low-temperature time associated with this event =====

      if(ijump.eq.1) then

        if(barrierev(ineigh).le.(0.d0))then
           write(*,*) 'SpawnID ',spawnID
     x      ,' detecting negative barrier ('
     x      ,barrierev(ineigh),') for istate,ineigh = '
     x      ,iofficial,ineigh
            barrierev(ineigh) = 100.0 !RJZ Set the barrier too high (band-aid)
        endif
        timelowevent =
     x      timehighnow*exp(barrierev(ineigh)*(betalow-betahigh))
        if(timelowevent-tickprev.le.0.d0)then
          timeloweventN = timelowevent
          timelowevent = (timehighnow-timehighnowI)
     x        *exp(barrierev(ineigh)*(betalow-betahigh)) + tickprev
          tickprev_check = (1.d0/fstar)
     +              *(timehighnowI*fstar)**(betalow/betahigh)
          if(newneighbor.eq.1)then
            write(*,*) 'spawnid ',spawnID
     x        ,' has NEGATIVE timeloweventN = ',timeloweventN
     x        ,' -- barrier = ',barrierev(ineigh)
     x        ,' -- ineigh = ',ineigh
     x        ,' -- timehighnow = ',timehighnow
     x        ,' -- tickprev = ',tickprev
     x        ,' -- tickprev_check = ',tickprev_check
     x        ,' -- CHANGED to timelowevent = ',timelowevent
            call flushtad(6)
          endif
        endif
      end if

c ===== update the nattempts-based prefactor if appropriate  (9/04)
c note change to make sure the prefac is ready before converting to synth mode  9/20/04

      if(ijump.eq.1  .and. ratemode.eq.'nattempts') then
         if (nattempted(ineigh).ge.nattemptsynth) then
c         Note that the prefac evolves (improves) as more attempts occur
            timetotx=timehighnow+timehighprev
            ratehigh=float(nattempted(ineigh))/timetotx
            preold=prefac(ineigh)
            prefac(ineigh)=ratehigh/exp(-barrierev(ineigh)*betahigh)
            if(ipr.ge.6) write(6,2234) iofficial,
     x          ineigh,nattempted(ineigh),preold,prefac(ineigh)
 2234       format('prefactor update: iofl,ineigh,nattempts,old,new:',
     x          i6,i5,i6,1p2d15.3)
         endif
      end if
c ====== place a tick mark on the low-temperature time line for this event ======

c (for ptad, only put a tick mark if it is a new type of event)
c (for synthtad, put an extrapolated tick mark for first occurence, and
c  let code below add synthetic tick marks for it after that)
c xxx think about this some more - is this what we want?

      if((lfirstinstate).and.(spawnIDlocal.eq.0)) lputyes = .true. ! Should "put" at end..

      if(ijump.eq.1) then

       !! Right Now spectad == specsynthtad:
       if(parallelmode.eq.'spectad')then
         if(ntick(ineigh).eq.1)then
           if(isynthclass(ineigh).eq.0)then
             if(ticknext(1,ineigh).gt.timelowevent)then
               ntick(ineigh) = 1
               ticknext(1,ineigh) = timelowevent
               !previousbest = timelowevent - tickprev
               !timelowevent = timelowevent + dt
               lputyes = .true.
             endif
           endif
         elseif((isynthclass(ineigh).eq.0).or.(newneighbor.eq.1)) then
           ntick(ineigh) = 1
           ticknext(1,ineigh) = timelowevent
           lputyes = .true.
         endif
       else if(tadmode.eq.'tad'.or.tadmode.eq.'newtad'.or.
     +       tadmode.eq.'dimeretad'.or.tadmode.eq.'protecteddimeretad'
     +       .or.tadmode.eq.'userdimeretad') then
          ntick(ineigh)=ntick(ineigh)+1
          if(ntick(ineigh).gt.ntickmax) then
            write(6,*) 'note - hit ntickmax for ineigh =',ineigh
            write(6,*) 'danger - expect unexpected consequences'
            ntick(ineigh)=ntick(ineigh)-1
          else
            ticknext(ntick(ineigh),ineigh)=timelowevent
          end if
       else if(tadmode.eq.'synthtad'.or.tadmode.eq.'synthdimeretad')
     x         then
          if(isynthclass(ineigh).eq.0 .or. newneighbor.eq.1) then
             ntick(ineigh)=ntick(ineigh)+1
             if(ntick(ineigh).gt.ntickmax) then
                write(6,*) 'note - hit ntickmax for ineigh =',ineigh
                write(6,*) 'danger - expect unexpected consequences'
                ntick(ineigh)=ntick(ineigh)-1
             else
                ticknext(ntick(ineigh),ineigh)=timelowevent
             end if
          end if
        else if(tadmode.eq.'ptad') then
          if(newneighbor.eq.1) then
            if(ratemode.eq.'vineyard') pre=prefac(ineigh)
            if(ratemode.eq.'fixedprefac') pre=fixedprefac
            ratelow=pre*exp(-barrierev(ineigh)*betalow)
            if(rank_f.eq.0)
     x        timelowdrawn=-(1.0d0/ratelow)*log(1.0d0-prngen(0))
            call MPI_BCAST(timelowdrawn,1,MPI_REAL8,0,force_comm,ier)
            ntick(ineigh)=ntick(ineigh)+1
            if(ipr.ge.3) then    ! xxx make this ipr=4 later
              write(6,*) 'drawing ptad tick mark; ineigh,tau,tick:',
     x        ineigh,1.0d0/ratelow,timelowdrawn
            end if
            if(ntick(ineigh).gt.ntickmax) then
              write(6,*) 'note - hit ntickmax for ineigh =',ineigh
              write(6,*) 'youch - this should never happen for ptad'
              stop 'ntickmax hit for ptad - something goofed up'
            else
              ticknext(ntick(ineigh),ineigh)=timelowdrawn
            end if
          else
            continue
          end if
        end if
      end if
      call flushtad(6)

c ======== prepare some things to be ready to continue the trajectory ========
c if rejecting an event, reposition correctly

      if ((ijump.eq.1).or.(ibadtransyes.eq.1)) then
        !if(nprocs_l.eq.nforcecores)then
        if(nParRep.eq.1)then ! Don't need to prepare new coords if you are master



          call vecmov(xyz4,xyz,3*natom)
          !call vecmov(xyz1,xyz,3*natom)


          iop=0
          iwrt=-1
          call gaussp(natom,nmove,pxyz,itype,maxtyp,amu,ntype,
     x                moving,mx,my,mz,temphigh,iop,iwrt)
          call MPI_BCAST(pxyz,3*nmove,MPI_REAL8,0,force_comm,ier)
          call mdreset(0)
        endif
        ineedblack=nblack
      else
        ineedblack=0
      end if
      if(rank_f.eq.0) call flushtad(6)

c ======== major block to set up tickmarks and accept transitions  ========

 333  continue

c === check whether dimer barrier still holds up as lowest barrier ===
c xxxxb might want to move this later

      call findminbarrier(nneigh,barrierev,barrierminev,barrevev ! xxx now checks for reverse barrier, not count barriers with reverse barrier too small
     +    ,ereverse)
      if(itad.eq.1 .and. barrierminev.lt.bminuse) then
         write(6,*) 'caution/note: barrierminev<bminuse',istate,iblock
      end if
      if(tadmode.eq.'dimeretad'.or.tadmode.eq.'synthdimeretad'.or
     +    .tadmode.eq.'userdimeretad'.or
     +    .(tadmode.eq.'protecteddimeretad'.and.itrans.gt.1)) then
         if(barrierminev-dimerbarrierminev.lt.-1.d-2) then ! xxx used to be 1d-4, franz changed it, should it be a parameter?
            if(ipr.ge.0) then
               write(6,*)
     +             'note - DIMER did not find lowest barrier,',
     +             ' iofficial=',iofficial,', irevisited=',irevisited
               write(6,*) barrierminev,dimerbarrierminev
               dimerbarrierminev=barrierminev
            end if
         end if
      end if


c =========== decide whether to promote any neighbors to synthetic class ============

      if((tadmode.eq.'synthtad'.or.tadmode.eq.'synthdimeretad')
     +          .and.(parallelmode.ne.'spectad')) then ! Dont promote in spectad for now
         do jneigh=1,nneigh
            if(isynthclass(jneigh).eq.0) then
               if(synthrequirement.eq.'always') then
                  isynthclass(jneigh)=1
               else if(synthrequirement.eq.'barrierheight') then
                  if(barrierev(jneigh).le.synthbarrier)
     +                isynthclass(jneigh)=1
               else if(synthrequirement.eq.'nattempts') then
                  if(nattempted(jneigh).ge.nattemptsynth)
     +                isynthclass(jneigh)=1
               else if(synthrequirement.eq.'firstchosen') then
               else if(synthrequirement.eq.'nattempts') then
                  if(nattempted(jneigh).ge.nattemptsynth)
     +                isynthclass(jneigh)=1
               else if(synthrequirement.eq.'firstchosen') then
                  if(naccepted(jneigh).ge.1)
     +                isynthclass(jneigh)=1
               else
                  write(6,*) 'synthrequirement =',synthrequirement
                  stop 'synthrequirement not recognized'
               end if
               if (isynthclass(jneigh).eq.1.and.ivineyard.eq.3) then
                  ngradcallhold=ngradcall
                  cputimehold=cputime
                  call vineyardsubset(natom,nmove,xyz1,xyzneighsads(3
     +                *natom*jneigh+1),itype,taxes,ietype,maxtyp,rcut
     +                ,amu,rdynmatcrit,e,grad,prefacx,nelements,ierr
     +                ,ipr)
                  prefac(jneigh)=prefacx

                  ngradcallvineyard=ngradcallvineyard+ngradcall
     +                -ngradcallhold
                  if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
                  call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
                  cputimevineyard=cputimevineyard+cputime-cputimehold
                  write(6,6543) istate,ineigh,rdynmatcrit
     +                ,prefacx,nelements,prefacx
                  call flushtad(6)
               endif
            end if
            if (isynthclass(jneigh).eq.1.and.ntick(jneigh).gt.1) then ! xxx erase all but most recent extrapolated ticks
               ntick(jneigh)=1
            end if
         enddo
      end if

c ======= assign any synthetic tickmarks not already assigned =======
      if(tadmode.eq.'synthtad'.or.tadmode.eq.'synthdimeretad'
     +  .and.(parallelmode.ne.'spectad')) then
         do jneigh=1,nneigh
            if (isynthclass(jneigh).eq.1.and.ntick(jneigh).eq.0) then
               if(ratemode.eq.'vineyard') pre=prefac(jneigh)
               if(ratemode.eq.'fixedprefac') pre=fixedprefac
               if(ratemode.eq.'nattempts') pre=prefac(jneigh)
               ratelow=pre*exp(-barrierev(jneigh)*betalow)
               if(rank_f.eq.0)
     +           timelowdrawn=-(1.0d0/ratelow)*log(1.0d0-prngen(0))
               call MPI_BCAST(timelowdrawn,1,MPI_REAL8,0,force_comm,ier)
               ntick(jneigh)=1
               if(ipr.ge.4) then ! xxx make this ipr=4 later
                  write(6,*)
     +                'drawing synthtad tick mark; jneigh,tau,tick:'
     +                ,jneigh,1.0d0/ratelow,timelowdrawn
               end if
               ticknext(1,jneigh)=timelowdrawn
            endif
         enddo
      endif

c ============ update ineighshortest ============

      timelowshortest=1.d99
      ineighshortest=0
      do i=1,nneigh
         if(ntick(i).gt.0 .and. ticknext(1,i).lt.timelowshortest) then
            timelowshortest=ticknext(1,i) ! xxxxbu in synth mode, this isn't calculated right to be delta t
c                         xxxafv - I think the problem is not here, but in the times  file
c                                  write.  See fix below (8/26/04)
            ineighshortest=i
         end if
      end do

      isynjump = 0
      if(ineigh.gt.0) isynjump = isynthclass(ineigh)

      if((spawnDepth.lt.max_depth).and.
     x   (parallelmode.eq.'spectad').and.
     x   ((tadtime+timelowevent-tickprev).lt.currenttadtimestop).and.
     x   (ijump.eq.1).and.
     x   ((newneighbor.eq.1)
     x   .or.((nprocs_l.gt.nforcecores).and.(isynjump.eq.0))
     x   .or.(lsynthattempts.and.(isynjump.eq.0)))
     x   )then

        ! Check if we have been killed
        if(rank_f.eq.0) call lp_check_kill()
        call MPI_BCAST(spawnwrite,1,MPI_LOGICAL,0,force_comm,ier)
        call MPI_BCAST(activated,1,MPI_LOGICAL,0,force_comm,ier)
        if(.not.activated) goto 611 !CRZ

        isneighshare = 0
        if((timelowevent-tickprev.lt.previousbest))then ! A real event will change the t_low
          if(lsynthspawned)then
             if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
     x          ,' (groupID',groupID,') KILLING SYNTH-SPAWN @ t='
     x          ,previousbest,' for t=',timelowevent-tickprev
     x          ,' with bar = ',barrierev(ineigh)
             lsynthspawned = .false.
             lputyessynthkilled = .true.
             if((nneigh-ineighdep-1).gt.0)then
               isneighshare = iofficial
             endif
          elseif(spawnIDlocal.gt.0)then
             if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
     x          ,' (groupID',groupID,') KILLING NORMAL-SPAWN @ t='
     x          ,previousbest,' for t=',timelowevent-tickprev
     x          ,' with bar = ',barrierev(ineigh)
          endif
          if(nneigh.ne.ineighshortest)then
            if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
     x        ,' ineighshortest = ',ineighshortest
     x        ,' ntick(ineighshortest) = ',ntick(ineighshortest)
     x        ,' ticknext(1,ineighshortest) = '
     x        ,ticknext(1,ineighshortest)
     x        ,' tadtime = ',tadtime
     x        ,' tickprev = ',tickprev
     x        ,' nneigh = ',nneigh
     x        ,' ntick(nneigh) = ',ntick(nneigh)
     x        ,' ticknext(1,nneigh) = ',ticknext(1,nneigh)
     x        ,' timelowevent = ',timelowevent
          endif
        endif

        ! SpecTAD -- Create New Spawns and Kill Old Spawns
        istatespawn = 0
        ldepspawn = 0
        !isneighshare = 0 !<- Set above (incase we are killing a sythetic spawn)
        call lp_spawn_update(ijump,timelowevent-tickprev,
     x          tadtime,xyz2,nmove,istatespawn,ldepspawn,isneighshare)

      elseif((ijump.eq.1).and.
     x   (spawnDepth.lt.max_depth).and.
     x   (parallelmode.eq.'spectad').and.
     x   (spawnIDlocal.eq.0).and.(dodepositionrun)
     x             .and.(.not.skipdeposition)
     x   )then

             ! Check if we have been killed
             if(rank_f.eq.0) call lp_check_kill()
             call MPI_BCAST(spawnwrite,1,MPI_LOGICAL,0,force_comm,ier)
             call MPI_BCAST(activated,1,MPI_LOGICAL,0,force_comm,ier)
             if(.not.activated) goto 611 !CRZ

             ! Time to spawn new deposition
             if(rank_f.eq.0) xyzdep(1) = prngen(0)*taxes(1)
             if(rank_f.eq.0) xyzdep(2) = prngen(0)*taxes(2)
             call MPI_BCAST(xyzdep(1),2,MPI_REAL8,0,force_comm,ier)

             highestatom = 0.d0
             do idum = 1, nmove
               if(xyz1((idum-1)*3+3).gt.highestatom)then
                 dxrz = abs( xyz1((idum-1)*3+1) - xyzdep(1) )
                 dyrz = abs( xyz1((idum-1)*3+2) - xyzdep(2) )
                 !if((dxrz.lt.(3.0)).and.(dyrz.lt.(3.0)))then
                   highestatom = xyz1((idum-1)*3+3)
                 !endif
               endif
             enddo
             xyzdep(3) = highestatom !+ 2.0
             call vecmov(xyz1(1),xyzdep(4),3*nmove)
             if(rank_f.eq.0) write(6,*) 'SpawnID ',spawnID,
     x         ' SPAWNING DEPOSITION EVENT!! (outside synth block) at',
     x         ' currenttadtimestop = ',currenttadtimestop
             call flushtad(6)

             ldepspawn = 1
             isneighshare = 0
             call lp_spawn_update(1,currenttadtimestop-tadtime
     x          ,tadtime,xyzdep,nmove+1,0,ldepspawn,isneighshare)

             ntick(ineighdep) = 1
             ticknext(1,ineighdep) =
     x         (currenttadtimestop-tadtime)+tickprev
             ineighshortest = ineighdep

      elseif((ijump.eq.0).and.(iblock.eq.iblockdeposit).and. ! Spawn next dep if we have reached 10 MD blocks with no transition
     x   (spawnDepth.lt.max_depth).and.
     x   (parallelmode.eq.'spectad').and.
     x   (spawnIDlocal.eq.0).and.(dodepositionrun)
     x             .and.(.not.skipdeposition)
     x   )then

             ! Check if we have been killed
             if(rank_f.eq.0) call lp_check_kill()
             call MPI_BCAST(spawnwrite,1,MPI_LOGICAL,0,force_comm,ier)
             call MPI_BCAST(activated,1,MPI_LOGICAL,0,force_comm,ier)
             if(.not.activated) goto 611 !CRZ

             ! Time to spawn new deposition
             if(rank_f.eq.0) xyzdep(1) = prngen(0)*taxes(1)
             if(rank_f.eq.0) xyzdep(2) = prngen(0)*taxes(2)
             call MPI_BCAST(xyzdep(1),2,MPI_REAL8,0,force_comm,ier)
             !write(6,*) 'x,y deposition coords: ',xyzdep(1),xyzdep(2)
             highestatom = 0.d0
             do idum = 1, nmove
               if(xyz1((idum-1)*3+3).gt.highestatom)then
                 dxrz = abs( xyz1((idum-1)*3+1) - xyzdep(1) )
                 dyrz = abs( xyz1((idum-1)*3+2) - xyzdep(2) )
                 !if((dxrz.lt.(3.0)).and.(dyrz.lt.(3.0)))then
                   highestatom = xyz1((idum-1)*3+3)
                 !endif
               endif
             enddo
             xyzdep(3) = highestatom !+ 2.0
             !write(6,*) 'z deposition coord: ',xyzdep(3)
             call vecmov(xyz1(1),xyzdep(4),3*nmove)

             ldepspawn = 1
             isneighshare = 0
             if(rank_f.eq.0) write(6,*) 'SpawnID ',spawnID,
     x         ' SPAWNING DEPOSITION EVENT!! (at iblockdeposit) at',
     x         ' currenttadtimestop = ',currenttadtimestop
             call lp_spawn_update(1,currenttadtimestop-tadtime
     x          ,tadtime,xyzdep,nmove+1,0,ldepspawn,isneighshare)

             ntick(ineighdep) = 1
             ticknext(1,ineighdep) =
     x         (currenttadtimestop-tadtime)+tickprev
             ineighshortest = ineighdep

      elseif((ijump.eq.1).and.
     x       (parallelmode.eq.'spectad').and.
     x       (timelowevent-tickprev.lt.previousbest))then

          previousbest = timelowevent-tickprev

      endif

        if((ijump.eq.1).and.(newneighbor.eq.1).and.
     +     (parallelmode.eq.'spectad').and.
     +     (irecognize.ge.1).and.(ivineyard.gt.0))then

           if(prefacsad.lt.0)then
             icnt = 0
             dimerbarrierev=barrierev(ineigh)
 217  continue
             igeometry=1        ! saddle
             nmove_tmp = nmove
             if(nmove.eq.natom) nmove_tmp = nmove-1
             nnegsad = 0
             freqlogssad=0.d0
             call dovineyard(natom,nmove_tmp,xyzs,itype,taxes,ietype
     +         ,maxtyp,rcut,amu,es,grad,nnegsad,nprod,freqlogssad
     +         ,igeometry,freqlogsmin,nnegmin,prefac(ineigh),ierr,ipr)

             ebarriertmp = dimerbarrierev
             if(ierr.gt.0)then
               if(rank_f.eq.0) write(*,*) 'SpawnID ',spawnID
     +           ,' (iofficial = ',iofficial
     +           ,') has problem with barrier to neighbor ',ineigh
             elseif(
     +           (abs(ebarriertmp-barrierev(ineigh))/barrierev(ineigh))
     +           .gt.0.025)then
               ! Should also initiate a roll check -- ensure barrier connects same minima
               if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
     +           ,' changing barrier to ',ebarriertmp
     +           ,'  - from - ',barrierev(ineigh)
     +           ,' for ineigh ',ineigh,' -- icnt = ',icnt
               barrierev(ineigh) = ebarriertmp

                 ! For now, lest ignore barriers that need to be changed this much..
                 filn=''
                 write(filn,*) spawnID
                 filn='sp.init.'//trim(adjustl(filn))//'.badsad.dat'
                 call storefile(natom,xyz1,pxyz,itype,
     +               taxes,filn)
                 if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
     +               , ': WARNING.. IGNORING THIS SADDLE !!!'
                 barrierev(ineigh) = 1000.0

             endif

             if(nnegsad.ne.1)then
               icnt = icnt + 1
               if(rank_f.eq.0) write(*,*) 'spawnID ',spawnid
     +           ,' Found saddle with nnegsad ',nnegsad
     +           ,' -- adding noise -- icnt ',icnt
               if(icnt.gt.4)then
                 filn=''
                 write(filn,*) spawnID
                 filn='sp.init.'//trim(adjustl(filn))//'.badsad.dat'
                 call storefile(natom,xyz1,pxyz,itype,
     +               taxes,filn)
                 !STOP 'THIS SADDLE IS NO GOOD !!!'
                 if(rank_f.eq.0) write(*,*) 'THIS SADDLE IS NO GOOD !!!'
                 if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
     +               , ': WARNING.. IGNORING THIS SADDLE !!!'
                 barrierev(ineigh) = 1000.0
               endif
               tgcrit = 1.d-10 !gcrit*0.00001
               tdvcrit = 1.d-10 !dvcrit*0.0001
               tdrcrit = 1.d-10 !drcrit*0.0001
               nmove_tmp = nmove
               if(nmove.eq.natom) nmove_tmp = nmove-1
               do idum=1,3*nmove_tmp
                if(rank_f.eq.0)then
                  xyzs(idum)=xyzs(idum)+(0.005d0)*gasdev(0)
                endif
               enddo
               call MPI_BCAST
     +               (xyzs(1),3*nmove_tmp,MPI_REAL8,0,force_comm,ier)
               idimerbail=0
               imodeds=1
               iranx=0
               irandir=1
               ngradcallhold=ngradcall
               idimerdof=3
               ierrdimer=0
               call dimer_search(natom,nmove_tmp,idimerbail,emin
     +          ,dimerbarrierminev,idimerdof,ndimerdof,rdimerdist
     +          ,ndimercoord,idimeroverunder,rdimerdisp,xyz1,xyzs,xyz1
     +          ,xyzprev,edimer,rdimer,itype,taxes,ietype,maxtyp,rcut
     +          ,hyperratio,barrevev,ereverse,gfac,irotalign,nsteep1
     +          ,ierrdimer,iranx,irandir,transcrit,itranscrit,tgcrit
     +          ,tdvcrit,tdrcrit,ipr,imodeds,ifundc,eig)
                ngradcallds=ngradcall-ngradcallhold
                dimerbarrierev=(edimer-emin)*27.21d0
               goto 217
             elseif(icnt.gt.0)then
               if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
     +           ,' FIXED nnegsad issue with noise/dimer -- '
     +           ,' dimerbarrierev = ',dimerbarrierev
     +           ,' -- OR emin,edimer = ',emin,edimer
     +           ,' -- ierrdimer = ',ierrdimer
             endif

           endif

           tmpvalue = prefac(ineigh)*exp(-barrierev(ineigh)*betalow)
           if(tmpvalue.gt.(lp_avail/(avgtadtimestop)))then ! Must be really fast to bother "putting" early
             nneighfast = nneighfast + 1
!             if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
!     +       ,' found neigh ',ineigh,' with transrate: ',tmpvalue
!     +       ,' prefac: ',prefac(ineigh)
!     +       ,' and barrierev: ',barrierev(ineigh)
!             if(rank_f.eq.0) call flushtad(6)
           endif



           ! IN ORDER TO CALL PUTSTATE BEFORE END ...
           !
           ! Need 2+ new/fast neighbors, after 3+ attempted transitions
           ! --- OR ---
           ! Need 1+ new/fast neighbors, after 4+ attempted transitions
           ! --- OR ---
           ! Found 1+ new/fast neighbors, and not first state visit
           ! --- OR ---
           ! Need to have just killed a Synthetic Spawn
           !
           ! ..Otherwise we are prone to "waste" many cores going into synthetic mode prematurely
           if(((nneighfast.ge.2).and.(nattempt.ge.3)).or.
     +        ((nneighfast.ge.1).and.(nattempt.ge.4)).or.
     +        ((nneighfast.ge.1).and.(.not.lfirstinstate)).or.
     +        lputyessynthkilled)then
             nneighfast = 0 ! Must find another fast neighbor to put again
             lputyes = .true. ! Should "put".. new neighbor found
             ineedgauss=1       ! required since we don't save pxyz
             ineedwarm=0
             lendput=.false.
            call putstate(casename,iofficial,natom,nmove,taxes,ntype,
     x       itype,barrierminev,barrierminevss,dimerbarrierminev
     +       ,bminlowerbound,emin,templow,temphigh,timehighprev
     +       ,timehighnow,timelowprev,tickprev,ieventlast,ineedblack
     +       ,ineedwarm,ineedgauss,ibadtransyes,freqlogsmin,nnegmin
     +       ,nprod,nattempt,nneigh,ntickmax,nneighmax,mxatom
     +       ,lenxyzneighs,xyz,xyz1,xyz4,xyzneighs,xyzneighsads
     +       ,rdimermodes,esaddle,eneigh,nattempted,naccepted,ninvolved
     +       ,nebfound,ticknext,ntick,barrierev,barrevev,prefac
     +       ,isynthclass,lfirstinstate,lputyes,lendput,nofficialmax)
             lputyes=.false.
             lputyessynthkilled=.false.
           endif
        endif

        ! TEMPORARY BAND-AID:
        ! Prevent large backlog of unread messages by checking often
        if((parallelmode.eq.'spectad').and.(irecognize.ge.1))then
          !if((ijump.eq.1).or.(MOD(iblock,2).eq.0))then
            call lp_update_official(nofficialmax)
            call MPI_BARRIER(force_comm,ier)
          !endif
        endif

      if((ijump.eq.1).and.(newneighbor.eq.1).and.
     +                    (parallelmode.eq.'spectad'))then
         if((irecognize.ge.1).and.(ivineyard.gt.0))then
             if(rank_f.eq.0)
     +        write(*,'(A,I9,A,I9,A,ES16.6,A,ES16.6,A,ES16.6)')
     +        'spawnID ',spawnID
     +       ,' found neigh ',ineigh,' with transrate: ',tmpvalue
     +       ,' prefac: ',prefac(ineigh)
     +       ,' and barrierev: ',barrierev(ineigh)
             if(rank_f.eq.0) call flushtad(6)
         else
             if(rank_f.eq.0) write(*,'(A,I9,A,I9,A,ES16.6)')
     +        'spawnID ',spawnID
     +       ,' found neigh ',ineigh
     +       ,' with barrierev: ',barrierev(ineigh)
             if(rank_f.eq.0) call flushtad(6)
         endif
      endif

c ==== for super-synth, find minimum non-synthetic barrier ====
      if (idosupersynth.eq.1) then
         call findminbarrierns(nneigh,barrierev,isynthclass
     +       ,naccepted,barrierminevns)
      endif

c ======= set special minimum barrier if appropriate (kludge of 6/23/05 afv/buber) =======

cc      if(ietype.eq.78 .and. iofficial.eq.1
      if(ietype.eq.78 .and. dabs(emin+2.192542).le.1.d-4) then
        bminspecialev=1.176
      else
        bminspecialev=0.0d0
      end if

      if(ibminspecialflag.ne.1 .and. bminspecialev.ne.0.0d0) then
c       xxx note that following write occurs only once, and cannot currently handle case of
c       xxx multiple states having bminspecials, since it only writes once
        write(6,456) bminspecialev, iofficial
  456   format(' ***** IMPORTANT NOTE:  bminspecialev =',f8.5,
     x        ' for iofficial =',i5)
        ibminspecialflag=1
      end if

c ======= check if we are ready to accept a transition =======

      iaccept=0
      call MPI_BARRIER(force_comm,ier)
      call MPI_BCAST(ineighshortest,1,MPI_INTEGER,0,force_comm,ier)
      if(ineighshortest.ne.0) then
        ipde=0 ! protecteddimeretad flag
        if(tadmode.eq.'protecteddimeretad') then
           if (itrans.eq.1) ipde=1
           if (itrans.gt.1) ipde=2
        endif
        if(parallelmode.eq.'spectad') then
          timelowx = ticknext(1,ineighshortest)
          ! Update stoptime if this event is better than any other events
          ! (note that depositions are not recorded as neighbors)
          timehighstop=(1.d0/fstar)*(fstar*timelowx)**(betahigh/betalow)
          if(.false.)then
            timelowsweptx=(timehighnow*fstar)**(betalow/betahigh)/fstar
            effslope=log(timehighnow/timelowsweptx)/(betalow-betahigh)
            bminuse=-effslope
            timetotx=timehighnow+timehighprev
            bminlowerbound=log((timetotx)*fstar)/betahigh  ! Eq. (6)
            bminusex=min(bminlowerbound,barrierminev)  ! Eq. (10)
            if(bminusex.gt.bminuse) bminuse=bminusex
            if(bminspecialev.gt.bminuse) bminuse=bminspecialev    !  6/05 - As in MVF02
            thsnewtad=timelowx*exp(bminuse*(betahigh-betalow))
            timehighstop=min(timehighstop,thsnewtad)
          endif
          if(timehighnow.ge.timehighstop)then
            iaccept=1
            call MPI_BARRIER(force_comm,ier)
            if(rank_f.eq.0) write(*,*)
     x        'SpawnID ',spawnID,' is accepting ineigh ',ineighshortest
     x        ,' at iblock = ',iblock
            if(rank_f.eq.0) call flushtad(6)
          endif
        else if(tadmode.eq.'tad'.or.ipde.eq.1) then
          timelowx = ticknext(1,ineighshortest)
          timehighstop=(1.d0/fstar)*(fstar*timelowx)**(betahigh/betalow)
          if(timehighnow.ge.timehighstop) iaccept=1
        else if(tadmode.eq.'synthtad') then
          timelowx = ticknext(1,ineighshortest)
          timehighstop=(1.d0/fstar)*(fstar*timelowx)**(betahigh/betalow)
c          if(timehighnow.ge.timehighstop) iaccept=1
          timelowsweptx=(timehighnow*fstar)**(betalow/betahigh)/fstar
          effslope=log(timehighnow/timelowsweptx)/(betalow-betahigh)
          bminuse=-effslope   ! for this mode, bminuse will be the neg. slope of sweeping line
c use whichever is better - newtad or synthtad    9/01/04
            timetotx=timehighnow+timehighprev
            bminlowerbound=log((timetotx)*fstar)/betahigh  ! Eq. (6)
            bminusex=min(bminlowerbound,barrierminev)  ! Eq. (10)
            if(bminusex.gt.bminuse) bminuse=bminusex
            if(bminspecialev.gt.bminuse) bminuse=bminspecialev    !  6/05 - As in MVF02
            thsnewtad=timelowx*exp(bminuse*(betahigh-betalow))
            timehighstop=min(timehighstop,thsnewtad)
            if(timehighnow.ge.timehighstop) iaccept=1
        else if(tadmode.eq.'newtad') then
          if(iequivcheck.eq.0) then
c           For following, see Montalenti and Voter, JCP 116, 4819 (2002).
            timetotx=timehighnow+timehighprev
            bminlowerbound=log((timetotx)*fstar)/betahigh  ! Eq. (6)
            bminuse=min(bminlowerbound,barrierminev)  ! Eq. (10)
            timelowx = ticknext(1,ineighshortest)
            timehighstop = timelowx*exp(bminuse*(betahigh-betalow))
            thstad=(1.d0/fstar)*(fstar*timelowx)**(betahigh/betalow) ! regular tad way
c            write(6,*) 'debug: newtadstop,oldtadstop: ',timehighstop
c     +          ,thstad
            timehighstop=min(timehighstop,thstad)                     ! take better of two
            if(timehighnow.ge.timehighstop) iaccept=1
            if(bminuse.ge.barrierev(ineighshortest)) iaccept=1
          else if(iequivcheck.eq.1) then
              tpuse=max(timehighprev,tprevequiv(iequiv))
              timetotx=timehighnow+tpuse
            bminlowerbound=log((timetotx)*fstar)/betahigh
            bminuse=min(bminlowerbound,barrierminev)
              bminuse=min(bminuse,bminequiv(iequiv))
            timelowx = ticknext(1,ineighshortest)
            timehighstop = timelowx*exp(bminuse*(betahigh-betalow))
            thstad=(1.d0/fstar)*(fstar*timelowx)**(betahigh/betalow)  ! regular tad way
            timehighstop=min(timehighstop,thstad)                     ! take better of two
            if(timehighnow.ge.timehighstop) iaccept=1
            if(bminuse.ge.barrierev(ineighshortest)) iaccept=1
          else
            write(6,*) 'unrecognized value for iequivcheck'
            write(6,*) 'iequivcheck =',iequivcheck
            stop 'unrecognized value for iequivcheck'
          end if
       else if(tadmode.eq.'dimeretad'.or.tadmode.eq.'userdimeretad'
     $         .or.ipde.eq.2) then
          if(iequivcheck.eq.0) then
c           For following, see Montalenti and Voter, JCP 116, 4819 (2002).
             timetotx=timehighnow+timehighprev
             bminlowerbound=barrierminev*2.d0 ! Eq. (6)
             bminuse=min(bminlowerbound,barrierminev) ! Eq. (10)
             timelowx = ticknext(1,ineighshortest)
             timehighstop = timelowx*exp(bminuse*(betahigh-betalow))
             if(timehighnow.ge.timehighstop) iaccept=1
             if(bminuse.ge.barrierev(ineighshortest)) iaccept=1
          else if(iequivcheck.eq.1) then
             tpuse=max(timehighprev,tprevequiv(iequiv))
             timetotx=timehighnow+tpuse
             bminlowerbound=barrierminev*2.d0
             bminuse=min(bminlowerbound,barrierminev)
             bminuse=min(bminuse,bminequiv(iequiv))
             timelowx = ticknext(1,ineighshortest)
             timehighstop = timelowx*exp(bminuse*(betahigh-betalow))
             if(timehighnow.ge.timehighstop) iaccept=1
             if(bminuse.ge.barrierev(ineighshortest)) iaccept=1
          else
             write(6,*) 'unrecognized value for iequivcheck'
             write(6,*) 'iequivcheck =',iequivcheck
             stop 'unrecognized value for iequivcheck'
          end if
       else if(tadmode.eq.'synthdimeretad') then
          if(iequivcheck.eq.0) then
c         For following, see Montalenti and Voter, JCP 116, 4819 (2002).
             timetotx=timehighnow+timehighprev
             bminlowerbound=2.d0*barrierminev
             bminuse=min(bminlowerbound,barrierminev) ! Eq. (10)
             timelowx = ticknext(1,ineighshortest)
             timehighstop = timelowx*exp(bminuse*(betahigh-betalow))
             if(timehighnow.ge.timehighstop) iaccept=1
c             if(bminuse.ge.barrierev(ineighshortest)) iaccept=1
c = super-synthetic test
             if (idosupersynth.eq.1.and. ! xxxx iblock.eq.0.and.isupersynth.gt.0) then
     +           (barrierminevss.gt.0.0d0.or.isupersynth.eq.1)) then
                if (barrierminevns.lt.barrierminevss)
     +              barrierminevss=barrierminevns
                bminuse=barrierminevns
c                if (bminuse.lt.barrierminevss) bminuse ! redundant?
c     +              =barrierminevss                    ! redundant?
                timelowx = ticknext(1,ineighshortest)
                timehighstop = timelowx*exp(bminuse*(betahigh
     +              -betalow))
                if (timehighnow.ge.timehighstop) iaccept=1
                barrierminevss=bminuse
             endif
          else if(iequivcheck.eq.1) then
c xxxxb need to put supersynthetic block here too, or in a separate place
             tpuse=max(timehighprev,tprevequiv(iequiv))
             timetotx=timehighnow+tpuse
             bminlowerbound=2.d0*barrierminev
             bminuse=min(bminlowerbound,barrierminev)
             bminuse=min(bminuse,bminequiv(iequiv))
             timelowx = ticknext(1,ineighshortest)
             timehighstop = timelowx*exp(bminuse*(betahigh-betalow))
             if(timehighnow.ge.timehighstop) iaccept=1
c             if(bminuse.ge.barrierev(ineighshortest)) iaccept=1
          else
             write(6,*) 'unrecognized value for iequivcheck'
             write(6,*) 'iequivcheck =',iequivcheck
             stop 'unrecognized value for iequivcheck'
          end if
          tdumxx=timehighstop
          if (ipr.ge.4.and.iaccept.eq.1) then
             write(6,*) 'bmin accept in synthdimeretad'
          endif
c     == check if we can accept based on sweeping ==
          timelowx = ticknext(1,ineighshortest)
          tdumx=(1.d0/fstar)*(fstar*timelowx)**(betahigh/betalow)
          if (tdumx.lt.timehighstop) then
             timehighstop=tdumx
             if(timehighnow.ge.timehighstop) iaccept=1
             if (ipr.ge.4.and.iaccept.eq.1) then
                write(6,*) 'sweeping accept in synthdimeretad'
             endif
          endif
          if(timehighnow.ge.timehighstop) iaccept=1
          if (ipr.ge.6.and.iaccept.eq.1)
     +        write(6,*) 'bmin stoptime=',tdumxx,'; pivoting stoptime='
     +        ,tdumx
       else if(tadmode.eq.'ptad') then
          timelowx = ticknext(1,ineighshortest)
          timehighstop =
     x     (1.0d0/fpstar)*(alpha*timelowx*fpstar)**(betahigh/betalow)
c           xxx can later put in bminuse extrapolating up from tick*alpha
c           xxx can also use iequivcheck in this
          bminuse=0.0d0
          if(timehighnow.ge.timehighstop) iaccept=1
       end if


      else ! if ineighshortest.eq.0
c   xxx do I still need something in here for synthtad?
!        if(.not.lfirstinstate) then
!          write(*,*) 'spawnID ',spawnID,' has ineighshortest = '
!     +        ,ineighshortest
!        endif
        timehighstop=1.d20
        if(timehighnow+timehighprev.gt.0.0d0)  then
           bminlowerbound=log((timehighnow+timehighprev)*fstar)/betahigh ! Eq. (6)
           if(tadmode.eq.'dimeretad') bminlowerbound=2.d0*barrierminev
           if(tadmode.eq.'synthdimeretad')bminlowerbound=2.d0
     +         *barrierminev
           if(tadmode.eq.'protecteddimeretad'.and.itrans.gt.1)
     +         bminlowerbound=2.d0*barrierminev
           if(tadmode.eq.'userdimeretad') bminlowerbound=2.d0
     +         *barrierminev
        else
           bminlowerbound=0.0d0
        end if
        bminuse=min(bminlowerbound,barrierminev) ! Eq. (10)
        if((tadmode.eq.'newtad' .or.tadmode.eq.'dimeretad'.or.
     $       tadmode.eq.'synthdimeretad'.or.tadmode.eq.'userdimeretad'
     $       .or.(tadmode.eq.'protecteddimeretad'.and.itrans.gt.1)).and
     $       . iequivcheck.eq.1) then
          tpuse=max(timehighprev,tprevequiv(iequiv))
          timetotx=timehighnow+tpuse
          if(timetotx.gt.0.0d0) then
             bminlowerbound=log((timetotx)*fstar)/betahigh
             if(tadmode.eq.'dimeretad'.or.tadmode.eq.'synthdimeretad'
     +           .or.tadmode.eq.'userdimeretad'
     +           .or.(tadmode.eq.'protecteddimeretad'.and.itrans.gt.1)
     +           )bminlowerbound=2.d0*barrierminev
          else
             bminlowerbound=0.0d0
          end if
          bminuse=min(bminlowerbound,barrierminev)
          bminuse=min(bminuse,bminequiv(iequiv))
       end if
      end if

      if(iaccept.eq.1)then
        if(.not.iprsafe)then
          if((.false.).and. ! Should we ignore waiting messages if we are far beyond t_stop??
     +           (timehighnow-timehighstop).ge.(dt*nmd*(nParRep-1)))then
            iaccept = 1
          else
            iaccept = 0 ! Not safe to accept the shortest trans.. delayed message might change things...
          endif
        endif
      endif

      if ((ipr.ge.10).and.(rank_f.eq.0)) then
         write(6,*) 'debug: tadmode: ',tadmode
         write(6,*) 'debug: bminuse: ',bminuse
         write(6,*) 'debug: timehighnow: ',timehighnow
         write(6,*) 'debug: timehighprev: ',timehighprev
         write(6,*) 'debug: bminlowerbound: ',bminlowerbound
         write(6,*) 'debug: barrierminev: ',barrierminev
         write(6,*) 'debug: tpuse: ',tpuse
         write(6,*) 'debug: timetotx: ',timetotx
c         write(6,*) 'debug: bminequiv: ',bminequiv(iequiv)
      endif

c ===== calculate confidence for accepted transitions for dimer modes =====

c      if (iaccept.eq.1.and.index(tadmode,'dimer').ne.0) then
      if (iaccept.eq.1) then
         effslope=log(timehighnow/timelowx)/(betalow-betahigh)
         deltaeff=exp(-timehighnow*prefacmin*exp(effslope*betahigh)) ! effslope < 0
         call MPI_BARRIER(force_comm,ier)
         if ((ipr.ge.2).and.(rank_f.eq.0)) write(6,*)
     +      'Given prefacmin=',prefacmin,', effective uncertainty for',
     +       ' this event=',deltaeff,
     +       ' (confidence=',(1.0d0-deltaeff)*100.0d0,'%)'
         if(rank_f.eq.0) call flushtad(6)
      endif

c =========== check if it is time to stop the whole run ===========

      if(((tadmode.eq.'tad')
     +  .and.((parallelmode.ne.'spectad').or.(irecognize.eq.0))).or.
     +    (tadmode.eq.'protecteddimeretad'.and.itrans.eq.1)) then
        timelowswept=(timehighnow*fstar)**(betalow/betahigh)/fstar
      else if(tadmode.eq.'synthtad'.or.tadmode.eq.'dimeretad'
     +                        .or.(parallelmode.eq.'spectad')) then ! xxxx dimeretad in two blocks!
        timelowsweptx=(timehighnow*fstar)**(betalow/betahigh)/fstar ! synsweep
        timelowswept = timehighnow*exp(-bminuse*(betahigh-betalow)) ! emin
        timelowswept=max(timelowswept,timelowsweptx)
      else if(tadmode.eq.'newtad'.or.tadmode.eq.'dimeretad'.or.
     +      tadmode.eq.'userdimeretad'.or.
     +      (tadmode.eq.'protecteddimeretad'.and.itrans.gt.1)) then
        timelowswept = timehighnow*exp(-bminuse*(betahigh-betalow))
      else if(tadmode.eq.'synthdimeretad') then
         timelowswept = timehighnow*exp(-bminuse*(betahigh-betalow)) ! emin
         timelowsweptx=(timehighnow*fstar)**(betalow/betahigh)/fstar ! synsweep
         timelowswept=max(timelowswept,timelowsweptx)
      else if(tadmode.eq.'ptad') then ! xxx this calculated above
        timelowswept =
     x   (1.0d0/alpha/fpstar)*(timehighnow*fpstar)**(betalow/betahigh)
      end if

      if((tadtime+timelowswept-tickprev.gt.tadtimestop)
     x      .and.(parallelmode.ne.'spectad')) then
        if(rank_f.eq.0)then
         write(6,*) 'terminating run based on tadtimestop'
         write(6,*) 'tadtime =',tadtime
         write(6,*) 'timelowswept =',timelowswept
         write(6,*) 'tadtime+timelowswept =',tadtime+timelowswept
         write(6,*) 'tadtimestop =',tadtimestop
        endif
        istopping=2
      else
        istopping=0
      end if

c ========= info write after every block =========

      if(ipr.ge.1) then
        if((istate.eq.1 .and. iblock.eq.0)  .or.
     x        (ipr.ge.3 .and. iblock.eq.0)) then
          if(rank_f.eq.0)then
            write(6,1121)
            write(42,1121)
            call flushtad(6)
            call flushtad(42)
          endif
        end if
 1121   format(
     x    'blokh  istate  iofl reord iqiv   iblk ijmp ingh ninv ineb ',   ! blokh, h=headerline
     x    'isyn    bar insh ssyn ',    ! xxxxa iaccept changed to ssyn - synthetic class of ineighshortest
     x    'bminuse   timehigh timehistop   tlowswpt  tickshort    ',
     x    'tadtime     nforce newforce')
c       xxxxa sorry blas, I pulled out the column numbers since I was afraid they would end up wrong
        if(istate.eq.1 .and. iblock.eq.0) ngradlast=0
        if(ijump.eq.1) then
          bar=barrierev(ineigh)
          ninv=ninvolved(ineigh)
          isyn=isynthclass(ineigh)
          ineb=nebfound(ineigh)
        else if(ijump.eq.0) then
          bar=0.0d0
          ninv=0
          isyn=0
          ineb=0
        else
          stop 'illegal value for ijump'
        end if
        issyn=0
        if(ineighshortest.ne.0) issyn=isynthclass(ineighshortest)

!        if(ipr.ge.3 .or.
!     x    (ipr.ge.2 .and. (ijump.eq.1 .or.iaccept.eq.1)) .or.
!     x    (ipr.ge.1 .and. newneighbor.eq.1) ) then
        if(.true.)then !! RJZ -- Temporary to see more data with a low ipr level
              cdum6='blok  '
              if(ijump.eq.1) cdum6='blokj '         ! j=jump
              if(newneighbor.eq.1) cdum6='blokn '   ! n=new jump
              if(iaccept.eq.1) cdum6='bloka '       ! a=accept
              if(iaccept.eq.1.and.iblock.eq.0) cdum6='blokai'
              if(istopping.eq.2) cdum6='bloks '     ! s=stopping whole run
            do iwx=1,2
            if(iwx.eq.1) iunit=6
            if(iwx.eq.2) iunit=42
            if(rank_f.eq.0)then
            write(iunit,1122) cdum6,
     x      istate,iofficial,ireordermatch,iequiv,iblock,ijump,ineigh,
     x      ninv,ineb,isyn,
     x      bar,ineighshortest,issyn,bminuse,timehighnow,timehighstop,
     x      timelowswept,timelowshortest,tadtime,ngradcall,
     x      ngradcall-ngradlast
            endif
            end do
 1122     format(a6,i7,2i6,i5,i7,5i5,f7.3,2i5,f8.3,1p5d11.3,i11,i9) ! afv 9/13/11

          ngradlast=ngradcall
          call flushtad(6)
          call flushtad(42)
        end if
      end if

c === write out tickmarks just for user reference ===

      if(rank_f.eq.0)then
      if(iaccept.eq.1 .or. istopping.eq.2) then
        if(ipr.ge.6) then
          write(6,*) 'final tickmarks for state ',istate,iofficial
          write(6,*)
     x         '  istate iofficial ineigh     tick      isyn shortest'
          do i=1,nneigh
            cdum6='      '
            if(i.eq.ineighshortest) cdum6='     *'
            write(6,132)
     x          istate,iofficial, i,ticknext(1,i),isynthclass(i),cdum6
  132       format('tickmarks: ',2i8, i8,1pd15.5,i5,a6)
          end do
        end if
      end if
      endif
      call MPI_BARRIER(force_comm,ier)

c ==== jump to 600 now if stopping based on total TAD time
      call MPI_BCAST(istopping,1,MPI_INTEGER,0,force_comm,ier)
      if(istopping.eq.2) goto 600

c =========== if we are not accepting some transition, jump back up to continue integrating ===========
      call MPI_BCAST(iaccept,1,MPI_INTEGER,0,force_comm,ier)
      if(iaccept.eq.0) go to 222

c ================== IF WE GET TO HERE, WE ARE ACCEPTING A TRANSITION ================

  444 tadtime=tadtime+ticknext(1,ineighshortest)-tickprev
      call MPI_BCAST(tadtime,1,MPI_REAL8,0,force_comm,ier)
      if((ticknext(1,ineighshortest)-tickprev).lt.0.d0)then
        if(rank_f.eq.0) write(*,*) 'spawnID ',spawnID
     x            ,' ERROR -- tadtime going NEGATIVE!'
        STOP ' ERROR -- tadtime going NEGATIVE!'
      endif
      if(tadtime.gt.currenttadtimestop)then
        ! There is no spawn after this one, adjust the actual tadtime
        tadtime = currenttadtimestop ! Real boost should only consider the tadtime<= currenttadtimestop
      endif
      naccepted(ineighshortest)=naccepted(ineighshortest)+1
      ieventlast=ineighshortest
      nneighlast=nneigh

c ========= set tickprev for next time =========
      dtimelow=ticknext(1,ineighshortest)-tickprev

      if(((tadmode.eq.'tad')
     +    .and.((parallelmode.ne.'spectad').or.(irecognize.eq.0)))
     +    .or.(tadmode.eq.'protecteddimeretad'.and.itrans.eq.1)) then
        tickprev=0.0d0
      else if(tadmode.eq.'synthtad'.or.tadmode.eq.'synthdimeretad'.or.
     +                          (parallelmode.eq.'spectad')) then
         if (isynthclass(ineighshortest).eq.0) then
            tickprev=0.0d0
         else
            tickprev=ticknext(1,ineighshortest)
         endif
      else if(tadmode.eq.'newtad') then
        tickprev=0.0d0
      else if(tadmode.eq.'dimeretad'.or.tadmode.eq.'userdimeretad'.or ! xxx this should be here, right?
     +      .(tadmode.eq.'protecteddimeretad'.and.itrans.gt.1)) then
         tickprev=0.0d0
      else if(tadmode.eq.'ptad') then
        tickprev=ticknext(1,ineighshortest)
      end if

c ==== save min of this basin for spring quench of next basin if itranscrit=1 ====
c xxxxa is this impacted by the reorder stuff?
      if (itranscrit.eq.1.or.ndimer.gt.0.or.ndimeretad.gt.0) then
         call vecmov(xyz1,xyzprev,3*natom)
      endif

c ====== print out some info for this state we are leaving ======

      if((ipr.ge.1).and.(rank_f.eq.0))then
        write(6,149) ineighshortest,timehighnow,
     +    ticknext(1,ineighshortest),tadtime
      endif
  149 format(' accepting trans ',i3,' at timehighnow=',1pd11.4,
     x       ' tlow=',1pd11.4,' TAD time =',1pd11.4)

      if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
      call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
      cpudelta=cputime-cputimelast
      cputimelast=cputime
      cputimesave=cputime
      cpudeltasave=cpudelta

c xxx will need to limit this summary stuff if in synthetic mode
      if((ipr.ge.1).and.(spawnwrite).and.(rank_l.eq.0))
     x    write(6,162) iofficial, iequiv,istate, emin*27.21d0,
     x    barrierminev,timehighprev,timehighnow,timehighnow+timehighprev
  162 format('------------------------------',/
     x    ' information on state we are leaving:',/
     x    ' info: iofficial =',i5,', iequiv =',i5,
     +                                ', istate =',i6,/
     +    ' info: basin minimum energy =',f16.5,' eV',/
     +    ' info: lowest barrier so far =',f12.5,' eV',/
     +    ' info: timehighprev =',1pd15.6,' timehighnow =',1pd15.6,/
     +    ' info: timehightot =',1pd15.6)

c write out force call information for this state
c xxx will need to limit this summary stuff if in synthetic mode

      ngradcallwarmtot=ngradcallwarmtot+ngradcallwarm
      ngradcallblacktot=ngradcallblacktot+ngradcallblack
      ngradcalldimertot=ngradcalldimertot+ngradcalldimer
      ngradcallnmdtot=ngradcallnmdtot+ngradcallnmd
      ngradcallrefinetot=ngradcallrefinetot+ngradcallrefine
      ngradcallminimizetot=ngradcallminimizetot+ngradcallminimize
      ngradcalltranstot=ngradcalltranstot+ngradcalltrans
      ngradcallnebtot=ngradcallnebtot+ngradcallneb
      ngradcallleaktot=ngradcallleaktot+ngradcallleak
      ngradcallvineyardtot=ngradcallvineyardtot+ngradcallvineyard

      if((ipr.ge.1).and.(spawnwrite).and.(rank_l.eq.0))
     x        write(6,163) istate,
     x        ngradcallminimize,ngradcallminimizetot,
     x        ngradcallwarm, ngradcallwarmtot,
     x        ngradcallblack, ngradcallblacktot,
     x        ngradcallleak, ngradcallleaktot,
     x        ngradcalldimer, ngradcalldimertot,
     x        ngradcallnmd, ngradcallnmdtot,
     x        ngradcalltrans, ngradcalltranstot,
     x        ngradcallrefine, ngradcallrefinetot,
     x        ngradcallneb, ngradcallnebtot,
     x        ngradcallvineyard, ngradcallvineyardtot,
     x        ngradcall-ngradcallstate, ngradcall
 163  format('--------------------',/
     x    ' ngradcall summary for state:          ',i15,
     x                                         '      total so far'/
     x    ' ngradcall: minimizing basin:          ',i15,i20/
     x    ' ngradcall: warming basin:             ',i15,i20/
     x    ' ngradcall: blacking out transitions:  ',i15,i20/
     x    ' ngradcall: checking for leaks:        ',i15,i20/
     x    ' ngradcall: doing dimer searches:      ',i15,i20/
     x    ' ngradcall: doing MD loop:             ',i15,i20/
     x    ' ngradcall: detecting transition:      ',i15,i20/
     x    ' ngradcall: refining transition:       ',i15,i20/
     x    ' ngradcall: finding saddle (via NEB):  ',i15,i20/
     x    ' ngradcall: Hessian setup for vineyard:',i15,i20/
     x    ' ngradcall: total for state:           ',i15,i20/
     x    '------------------------------')
      if(rank_f.eq.0) cputime = mpi_wtime()-cput_ref
      call MPI_BCAST(cputime,1,MPI_REAL8,0,force_comm,ier)
      if((ipr.ge.1).and.(spawnwrite).and.(rank_l.eq.0))
     x    write(6,164) istate,cputimereorder,cputimevineyard,
     +    cputime-cputimestate,  cputime
 164  format('--------------------',/,
     x       ' cpu time summary for state:',i5,/
     x       ' cpu time for reorder check          =',f10.2,/,
     x       ' cpu time for Hessian setup and diag =',f10.2,/,
     x       ' cpu time total for state            =',f10.2,/,
     x       ' cpu time total for run so far       =',f10.2,/,
     x    '------------------------------')
      if(rank_l.eq.0) call flushtad(6)

c ======== update the .times file ========

c the one-time write statement using the following format is up above
  123 format('times file: times between transitions - case = ',a,/
     x '     i  istate jstate ioffcl iequiv Erel    Ea(eV) neigh',
     x ' nmoved synth    hightime  ',
     x '     nforce    TAD time      delta t      cpu        dcpu')


      iofficialx=iofficial
      if(ireordermatch.gt.0) iofficialx=-iofficial
      deltatadtime=tadtime-tadtimeprevious
      if(rank_l.eq.0)
     x    write(22,225) istate,istate,istate+1,iofficialx,iequiv,
     x    (emin-eofficial(1))*27.21d0,
     x    barrierev(ineighshortest),ineighshortest,
     x    ninvolved(ineighshortest),isynthclass(ineighshortest),
     x    timehighnow,ngradcall,tadtime,
     x    deltatadtime,cputimesave,cpudeltasave  ! after fix - xxxafv
ccccccccc     +    ,ticknext(1,ineighshortest),cputimesave,cpudeltasave  ! before fix
  225 format(i6,i6,i7,i6,i7,0pf9.3,0pf8.3,i5,i7,i6,1pd17.6,i10,1pd15.6,
     x                              1pd11.3,0pf10.1,0pf11.2)
      tadtimepreviousx=tadtimeprevious
      tadtimeprevious=tadtime
      if(rank_l.eq.0) call flushtad(22)

c ======== check if we have reached a special state that stops the run =======

      if(ietype.eq.78
     x   .and. dabs(27.21*(emin+2.192542)-0.272).le.1.d-3) then
        ispecialstop=1
      end if

c ======== update the .boost file ========

      iofficialx=iofficial
      if(ireordermatch.gt.0) iofficialx=-iofficial
      runboost=tadtime/dt/real(ngradcall)
      rstateboost=(tadtime-tadtimepreviousx)/dt/real(ngradcall
     +    -ngradcallstate)
      rcpuboost=tadtime/(cputimesave*callspersec*dt)
      boost_denom = callspersecmd*dt
      if(rank_l.eq.0)
     x    write(23,2253) istate,istate,istate+1,iofficialx,
     x    (emin-eofficial(1))*27.21d0,
     x    timehighnow,ngradcall,tadtime,
     x    rstateboost,runboost,rcpuboost  ! after fix - xxxafv
 2253 format(i6,i6,i7,i6,0pf9.3,1pd17.6,i10,1pd15.6,1pd14.3,1pd14.3
     +    ,1pd14.3)
      if(rank_l.eq.0) call flushtad(23)


c ====== update equivalent-energy-state info ======

      if(iequivcheck.ge.1) then
        bminequiv(iequiv)=min(bminequiv(iequiv),barrierminev)
        tprevequiv(iequiv)=tprevequiv(iequiv)+timehighnow
        if(tadmode.eq.'dimeretad'.or.tadmode.eq.'synthdimeretad'.or
     +      .tadmode.eq.'userdimeretad'.or
     +      .(tadmode.eq.'protecteddimeretad'.and.itrans.gt.1)
     +      )dimbminevequiv(iequiv)=barrierminev
        if(ipr.ge.3) then   ! xxx change this to 4 when tested/debugged
          write(6,405)
     x    istate,iequiv,bminequiv(iequiv),tprevequiv(iequiv)
  405     format('leaving istate',i6,
     x    '  putting back into iequiv =',i6,
     x    ' bminequiv=',f10.4,'  tprevequiv =',1pd10.3)
        end if
      endif

c ======= shift the used up stop time ("timehighstop") into previous time =======
c xxxx this could be rewritten in terms of synclass

      if(((tadmode.eq.'tad')
     +  .and.((parallelmode.ne.'spectad').or.(irecognize.eq.0)))
     +  .or.(tadmode.eq.'protecteddimeretad'.and.itrans.eq.1)) then
        timehighshift=timehighstop
        timelowshift=(timehighshift*fstar)**(betalow/betahigh)/fstar
      else if(tadmode.eq.'newtad'.or.tadmode.eq.'dimeretad'.or.tadmode
     +      .eq.'userdimeretad'.or.(tadmode.eq.'protecteddimeretad'.and
     +      .itrans.gt.1)) then
        timehighshift=timehighstop    !!! xxxxa should we also zero out timehigh now after this?
        timelowshift=timehighshift*exp(bminuse*(betalow-betahigh))
      else if((tadmode.eq.'synthtad').or.
     +                              (parallelmode.eq.'spectad')) then
         if (isynthclass(ineighshortest).eq.1) then
            timehighshift=0.0d0
            timelowshift=0.0d0
         else
            timehighshift=timehighstop
c            timelowshift=(timehighshift*fstar)**(betalow/betahigh)/fstar ! xxxx not correct when using bmin! 09/29/2005
            timelowshift=ticknext(1,ineighshortest)
         endif
      else if(tadmode.eq.'synthdimeretad') then
         if (isynthclass(ineighshortest).eq.1) then
            timehighshift=0.0d0
            timelowshift=0.0d0
         else
            timehighshift=timehighstop !!! xxxxa should we also zero out timehigh now after this?
c            timelowshift=timehighshift*exp(bminuse*(betalow-betahigh)) ! xxxx not correct for synth mode, fixed 03/09/24
            timelowshift=ticknext(1,ineighshortest)
         endif
      else if(tadmode.eq.'ptad') then
        timehighshift=0.0d0
        timelowshift=0.0d0
      else
        write(6,*) 'uh oh - not clear how to compute timelowshift'
      end if


      dThigh = max(0.d0,timehighnow-timehighnowI)
      timehighnow=timehighnow-timehighshift
      timehighprev=timehighprev+timehighshift
      timelowprev=timelowprev+timelowshift

c ========= time-shift all the tick marks =========
      do i=1,nneigh
        do j=1,ntick(i)
          ticknext(j,i)=ticknext(j,i)-timelowshift
        end do
      end do
c
c ========= create new synthetic tick for one about to be erased  ==========
       isyn=isynthclass(ineighshortest)
       if(tadmode.eq.'ptad'
     x     .or. (tadmode.eq.'synthtad' .and. isyn.eq.1)
     x     .or. (tadmode.eq.'synthdimeretad' .and. isyn.eq.1)
     x     ) then
!     x     .or. (parallelmode.eq.'spectad' .and. isyn.eq.1)
        if(ntick(ineighshortest).ne.1) then
          write(6,*) 'uh oh - ntickmax.ne.1 in ptad'
          stop       'uh oh - ntickmax.ne.1 in ptad'
        end if
        if(ratemode.eq.'vineyard') pre=prefac(ineighshortest)
        if(ratemode.eq.'fixedprefac') pre=fixedprefac
        if(ratemode.eq.'nattempts') pre=prefac(ineighshortest)
        ratelow=pre*exp(-barrierev(ineighshortest)*betalow)
        if(rank_f.eq.0)
     x   timelowdrawn=-(1.0d0/ratelow)*log(1.0d0-prngen(0))
        call MPI_BCAST(timelowdrawn,1,MPI_REAL8,0,force_comm,ier)
        ntick(ineighshortest)=ntick(ineighshortest)+1
        if (ipr.ge.4) write(6,*) 'timelowdrawn=',timelowdrawn
        if(ntick(ineighshortest).gt.ntickmax) then
          write(6,*) ntick(ineighshortest)
          write(6,*) 'yow - ntickmax exceeded in ptad or synthtad'
          stop 'yow - ntickmax exceeded in ptad or synthtad'
        end if
        ticknext(2,ineighshortest) =
     x                    ticknext(1,ineighshortest) + timelowdrawn
       end if

c ======= erase the used-up tick mark =======
       if(parallelmode.ne.'spectad')then
         ntick(ineighshortest)=ntick(ineighshortest)-1
         if(ipr.ge.2) write(6,*) 'erasing tickmark for neighbor ',
     x                     ineighshortest,ntick(ineighshortest)
         if(ntick(ineighshortest).lt.0)stop 'ntick(ineighshortest).lt.0'
         do j=1,ntick(ineighshortest)
          ticknext(j,ineighshortest)=ticknext(j+1,ineighshortest)
         end do
       else
         if(isyn.eq.1)then
           ticknext(1,ineighshortest)=ticknext(2,ineighshortest)
           do j=2,ntick(ineighshortest)
             ticknext(j,ineighshortest)=0.d0
           end do
           ntick(ineighshortest) = 1
         else
           do j=1,ntick(ineighshortest)
             ticknext(j,ineighshortest)=0.d0
           end do
           ntick(ineighshortest) = 0
         endif
       endif

c ========= store state file for this state ==========
c xxx perhaps also put an interim putstate call above, in case the system
c xxx spends a long time before escaping from this state, and then crashes
c xxx before exiting.  This would make restarting more efficient
      if(irecognize.ge.1) then
        ineedgauss=1       ! required since we don't save pxyz
        ineedwarm=0
        lendput=.true.
        if(parallelmode.eq.'spectad')then
          if(isynthclass(ieventlast).eq.0) then
            if(lsynthattempts) lputyes = .true.
            timehighprev_send = dThigh
            timehighnow_send = 0.d0
          else
            if(dThigh.gt.(1.d-14))then
              if(lsynthattempts) lputyes = .true.
            endif
            ! Still Send Thigh to timehighprev... fine if kMC uses timehighprev only for rates
            timehighprev_send = dThigh
            timehighnow_send = 0.d0
          endif
          call putstate(casename,iofficial,natom,nmove,taxes,ntype
     +     ,itype,barrierminev,barrierminevss,dimerbarrierminev
     +     ,bminlowerbound,emin,templow,temphigh
     +     ,timehighprev_send,timehighnow_send
     +     ,timelowprev,tickprev,ieventlast,ineedblack
     +     ,ineedwarm,ineedgauss,ibadtransyes,freqlogsmin,nnegmin
     +     ,nprod,nattempt,nneigh,ntickmax,nneighmax,mxatom
     +     ,lenxyzneighs,xyz,xyz1,xyz4,xyzneighs,xyzneighsads
     +     ,rdimermodes,esaddle,eneigh,nnewattempted,naccepted,ninvolved
     +     ,nebfound,ticknext,ntick,barrierev,barrevev,prefac
     +     ,isynthclass,lfirstinstate,lputyes,lendput,nofficialmax)
        else
          call putstate(casename,iofficial,natom,nmove,taxes,ntype
     +     ,itype,barrierminev,barrierminevss,dimerbarrierminev
     +     ,bminlowerbound,emin,templow,temphigh,timehighprev
     +     ,timehighnow,timelowprev,tickprev,ieventlast,ineedblack
     +     ,ineedwarm,ineedgauss,ibadtransyes,freqlogsmin,nnegmin
     +     ,nprod,nattempt,nneigh,ntickmax,nneighmax,mxatom
     +     ,lenxyzneighs,xyz,xyz1,xyz4,xyzneighs,xyzneighsads
     +     ,rdimermodes,esaddle,eneigh,nattempted,naccepted,ninvolved
     +     ,nebfound,ticknext,ntick,barrierev,barrevev,prefac
     +     ,isynthclass,lfirstinstate,lputyes,lendput,nofficialmax)
        endif
      endif

c ======== move the next state into position ========
c  xxxa - now with reordering
      if(parallelmode.ne.'spectad')then
       if(ireordermatch.ne.0) then
        call reorder(natom,nmove,
     x               xyzneighs((ineighshortest-1+1)*3*natom+1),
     x               kreorder, xyz)
       else
        call vecmov(
     x      xyzneighs((ineighshortest-1+1)*3*natom+1),xyz,3*natom)
       endif
      endif

  611 continue

         isynlast=0
         lstillsynth=.true.
         if(ieventlast.gt.0) isynlast=isynthclass(ieventlast)
         ! Did a neighbor with (nattempts < nattemptkMC) win?
         ! Reset the sweeping line to zero!
         if(isynlast.eq.1)then
           n_attempt = statedata(iofficial)%nattempted(ieventlast)
           if(n_attempt.lt.nattemptkMC)then
               if(ivineyard.eq.0) then
                   isynlast=0
                   write(*,*)'SpawnID ',SpawnID
     x                 ,' -> RESET SWEEP!! n_attempt=',n_attempt
               endif
           endif
!           else if(.not.cleanwarmup)then
!               isynlast=0
!               write(*,*)'SpawnID ',SpawnID
!     x             ,' -> RESET SWEEP!! cleanwarmup=',cleanwarmup
!           endif
         endif
         if((isynlast.eq.0).and.(ivineyard.eq.0)) then
           lstillsynth=.false.
           timehighnow=0.d0
           timehighprev=0.d0
           timelowprev=0.d0
           tickprev=0.d0
           timehighprev=0.d0
           do i=1,nneigh
             ntick(i)=0
             ticknext(1,i)=0.0d0
           enddo
         endif
         lsendtime = .false.
         if(.not.activated)then
           lstillsynth = .false.
         else
           if(.not.lfirstinstate)then
             lsendtime = .true.
           endif
         endif

         ! Don't use this time to sweep if we werent properly
         ! warmed up and/or blacked out:
!         if((lstillsynth).and.
!     x              ((.not.cleanwarmup).or.(.not.cleanblackout)))then
         !!if((.not.cleanwarmup).or.(.not.cleanblackout))then
         if(.not.cleanwarmup)then
             lsendtime = .false.
             write(*,*)'SpawnID ',SpawnID
     x           ,' -> DONT SEND TIME!! cleanwarmup=',cleanwarmup
     x           ,' - cleanblackout=',cleanblackout
         endif

      activated = .false. !crz -- tad instance done
      ! crz -- Some SpecTAD Stuff
      dThigh = dThigh - timehighshift !Gives Change in timehighnow
      lsenddone = .true.

  601 continue
      if((nprocs_l.gt.nforcecores).and.(rank_l.lt.nforcecores))
     +                                        call pr_send_deactivate()
      if((nprocs_l.gt.nforcecores).and.(rank_l.lt.nforcecores))
     +                                        call pr_clean_trans()
      if(lsenddone.and.(rank_l.eq.0))then
        write(*,*) 'spawnID ',spawnID,' MASTER is DONE '
      endif

!      if(rank_l.gt.0) write(*,*) 'spawnID ',spawnID,' rank_l ',rank_l
!     +               ,' is DONE '
!      if(rank_l.gt.0) call flush(6)
      call MPI_BCAST(activated,1,MPI_LOGICAL,0,local_comm,ier)
      if(lsenddone.and.(rank_l.lt.nforcecores))then
        if((ineighshortest.le.nneigh).and.(ineighshortest.ge.1))then
          winningbarrier = barrierev(ineighshortest)
        else
          winningbarrier = 0.d0
        endif
        call lp_done_mess(tadtime,boost_denom,cputimesave
     +               ,iofficial,lfirstinstate,lstillsynth,dtimelow
     +               ,dThigh,icurrentevent,nneigh,lsendtime,dt,nmd
     +               ,winningbarrier,emin) !crz -- send done message
      endif

  500 continue
c =============== end of master loop over accepted transitions ==============
c ================ jump to here if ending run ===============
  600 continue
      call MPI_BCAST(lfinished,1,MPI_LOGICAL,0,local_comm,ier)

      if(.not.lfinished)then
        write(*,*) 'ERROR... groupID ', groupID,
     +    ' has been ousted form the activation loop!!'
        STOP ' LP ousted from the activation loop!!'
      endif

c ====== one last write to times file if have reached the TAD stop time
      if((istopping.eq.2).and.(rank_f.eq.0)) then
        write(22,225) istate,istate,-999,iofficial,iequiv
     +    ,0.0d0,0.0d0,0
     +    ,0,0,timehighnow,ngradcall,tadtime+timelowswept
     +    ,0.0d0,cputimesave,cpudeltasave
        call flushtad(22)
      end if

      if(rank_l.eq.0) close(22) ! close .times file
      if(rank_l.eq.0) close(23) ! close .boost file

c ========== print out summary of run ===========

      if ((iequivcheck.eq.1).and.(spawnwrite).and.(rank_l.eq.0)) then
         write(6,*) '========================================'
         write(6,*) 'Summary of equiv Energy States Found'
         do i=1,nequiv
            write(6,148) i,eequiv(i),bminequiv(i),tprevequiv(i)
         enddo
 148     format('state=',i5,', energy=',f12.5,' h, min barrier=',
     x           f12.5,' (eV), total time high=',d12.5)
         write(6,*) '========================================'
      endif

      if ((irecognize.eq.1).and.(spawnwrite).and.(rank_l.eq.0)) then
         write(6,*) '========================================'
         write(6,*) 'Summary of recognized states'
         write(6,*) 'Total number of states visited: ',istate
         write(6,*) 'Total number of official states: ',nofficial
         write(6,*) '========================================'
      endif

      fracmin=  float(ngradcallminimizetot)/float(ngradcall)*100.0d0
      fracwarm= float(ngradcallwarmtot)/float(ngradcall)*100.0d0
      fracblack=float(ngradcallblacktot)/float(ngradcall)*100.0d0
      fracdimer=float(ngradcalldimertot)/float(ngradcall)*100.0d0
      fracnmd=  float(ngradcallnmdtot)/float(ngradcall)*100.0d0
      fractrans=float(ngradcalltranstot)/float(ngradcall)*100.0d0
      fracref=  float(ngradcallrefinetot)/float(ngradcall)*100.0d0
      fracneb=  float(ngradcallnebtot)/float(ngradcall)*100.0d0
      fracleak= float(ngradcallleaktot)/float(ngradcall)*100.0d0
      fracvin=  float(ngradcallvineyardtot)/float(ngradcall)*100.0d0

      ngradsum=ngradcallminimizetot+ngradcallwarmtot
     +    +ngradcallblacktot+ngradcallleaktot+ngradcalldimertot
     +    +ngradcallnmdtot+ngradcalltranstot+ngradcallrefinetot
     +    +ngradcallnebtot+ngradcallvineyardtot

      fracsum=float(ngradsum)/float(ngradcall)*100.0d0

      if ((spawnwrite).and.(rank_l.eq.0))
     +    write(6,165) ngradcallminimizetot,fracmin,ngradcallwarmtot
     +    ,fracwarm,ngradcallblacktot,fracblack,ngradcallleaktot
     +    ,fracleak,ngradcalldimertot,fracdimer,ngradcallnmdtot,fracnmd
     +    ,ngradcalltranstot,fractrans,ngradcallrefinetot,fracref
     +    ,ngradcallnebtot,fracneb,ngradcallvineyardtot,fracvin,ngradsum
     +    ,fracsum,ngradcall
 165  format('=======================================',/
     x    ' ngradcall summary for run             ',/
     x    ' ngradcall: minimizing basin:          ',i15,f10.2,/
     x    ' ngradcall: warming basin:             ',i15,f10.2,/
     x    ' ngradcall: blacking out transitions:  ',i15,f10.2,/
     x    ' ngradcall: checking for leaks:        ',i15,f10.2,/
     x    ' ngradcall: doing dimer searches:      ',i15,f10.2,/
     x    ' ngradcall: doing MD loop:             ',i15,f10.2,/
     x    ' ngradcall: detecting transition:      ',i15,f10.2,/
     x    ' ngradcall: refining transition:       ',i15,f10.2,/
     x    ' ngradcall: finding saddle (via NEB):  ',i15,f10.2,/
     x    ' ngradcall: Hessian setup for vineyard:',i15,f10.2,/
     x    '---------------------------------------',/
     x    ' ngradcall: sum of above:              ',i15,f10.2,/
     x    ' ngradcall: total for run:             ',i15,/
     x    '=======================================')

      if((spawnwrite).and.(rank_l.eq.0))
     x  write(6,*) 'TAD2 terminating normally - all done'

      call MPI_BARRIER(MPI_COMM_WORLD,ier)
      call MPI_Finalize (ier)

      stop 'TAD2 terminating normally - all done'
      end


      subroutine centerofmass(x,y,z,natom,xcm,ycm,zcm)
      implicit real*8(a-h,o-z)
      dimension x(natom),y(natom),z(natom)
      xcm=0.0d0
      ycm=0.0d0
      zcm=0.0d0
      do 10 i=1,natom
      xcm=xcm+x(i)
      ycm=ycm+y(i)
      zcm=zcm+z(i)
   10 continue
      xcm=xcm/float(natom)
      ycm=ycm/float(natom)
      zcm=zcm/float(natom)
      return
      end


      subroutine numberinvolved(xyz1,xyz2,nmove,transcrit,ninvolved)
c determine the number of atoms that moved more than transcrit between
c the two files - i.e., the number of atoms involved in a transition
c Assumes nothing funny has been done, like a periodic wraparound
      implicit real*8(a-h,o-z)
      dimension xyz1(3,nmove)
      dimension xyz2(3,nmove)

      ninvolved=0
      do i=1,nmove
        dist2=0.0d0
        do ix=1,3
          dist2=dist2+(xyz1(ix,i)-xyz2(ix,i))**2
        end do
        if(dist2.gt.transcrit**2) ninvolved=ninvolved+1
      end do

      return
      end

      character*(*) function chri(i)
c:  returns a character version of the integer i
c:  chri must be declared as character*n in calling routine
c:
      implicit real*8(a-h,o-z)
      character*5 fmt
      character*40 hold
      character*755 cdum
      common/chrcm/fmt,hold,cdum
      fmt='(i40)'
      write(hold,fmt) i
      call chrpak(hold,40,l)
      chri=hold(1:l)
      return
      end


      subroutine chrparse1(c,lmax,lc)
c: picks out the first word from a character string
c: lc is returned as the packed length.
c:
      character*(*) c
      character*800 ct
      common/chrcm/ct
      k=0

c find the beginning of the word
      do 10 i=1,lmax
      if(c(i:i).ne.' ') then
        ibegin=i
        go to 20
      end if
   10 continue
      stop 'parse1 failed'

c find the end of the word
   20 continue
      do 30 i=ibegin,lmax
      if(c(i:i).eq.' ') then
        iend=i
        go to 40
      end if
   30 continue
      stop 'parse1 failed 2'

c return the word in c, and its length in lc
   40 continue
      lc=iend-ibegin+1
      if(lc.gt.800) stop 'parse1 out of space'
      ct=c(ibegin:iend)
      c=ct

      return
      end


      subroutine chrpak(c,lmax,lc)
c: packs out blanks in a character string whose maximum length is lmax.
c: lc is returned as the packed length.
c:
      character*(*) c
      character*800 ct
      common/chrcm/ct
      k=0
      do 10 i=1,lmax
      if(c(i:i).ne.' ') then
        k=k+1
        if(k.gt.800) stop 'abort - maxed out in chrpak'
        ct(k:k)=c(i:i)
      end if
   10 continue
      lc=k
      c=ct(1:lc)
      return
      end

      subroutine clock(time)
c: Obtains CPU clock time in seconds, measured since beginning of run
      implicit real*8(a-h,o-z)
      real*4 etimetad,tarray(2)
      time=etimetad(tarray)
      return
      end

      subroutine storefile(natom,xyz,pxyz,itype,taxes,filnam)
c stores a cluster file

      use mod_mpi
      implicit real*8(a-h,o-z)
      dimension xyz(3,natom)
      dimension pxyz(3,natom)
      dimension taxes(3)
      dimension itype(natom)
      character*(*) filnam
      character*80 ctitle(1)

      common/storefilecom/istorefileopt

      iunit=27
      lf=len(filnam)

      call chrpak(filnam,lf,lfil)

      if(rank_f.eq.0)then ! Only first core of a given replica can write files (make only rank_l=0 ?)

      write(6,10) 'storefile: ',filnam
c add ipr on this later ,'   atom 1: ',(xyz(k,1),k=1,3)
 10   format(a11,a35,a11,3f12.5)

      open(unit=iunit,file=filnam(1:lfil),form='formatted',
     x              status='unknown')
      rewind iunit
      iwrt=-1
      ntitle=1
      ctitle(1)='cluster file from TAD2; meaning encoded in file name'
      iop=1  ! 1 gives a compact, no-momentum file   - new as of 9/24/02
      call putclsnew(natom,ntitle,ctitle,xyz,pxyz,itype,
     x                                 taxes,iunit,iwrt,iop)
      rewind iunit
      close(unit=iunit)

      if (istorefileopt.eq.1) then
         filnam=filnam//'.xyz'
         lf=len(filnam)
         call chrpak(filnam,lf,lfil)
         open(unit=iunit,file=filnam(1:lfil),form='formatted',
     x       status='unknown')
         rewind iunit
c         call putclsxyz(natom,xyz)
         rewind iunit
         close(unit=iunit)
      endif

      endif

      return
      end

      subroutine putstate(
     x casename,
     x istate,
     x natom,
     x nmove,
     x taxes,
     x ntype,
     x itype,
     x barrierminev,
     x barrierminevss,
     x dimerbarrierminev,
     x bminlowerbound,
     x emin,
     x templow,
     x temphigh,
     x timehighprev,
     x timehighnow,
     x timelowprev,
     x tickprev,
     x ieventlast,
     x ineedblack,
     x ineedwarm,
     x ineedgauss,
     x ibadtransyes,
     x freqlogsmin,
     x nnegmin,
     x nprod,
     x nattempt,
     x nneigh,
     x ntickmax,
     x nneighmax,
     x mxatom,
     x lenxyzneighs,
     x xyz,    ! BEGIN ARRAYS HERE
     x xyzmin,
     x xyz4,
     x xyzneighs,
     x xyzneighsads,
     x rdimermodes,
     x esaddle,
     x eneigh,
     x nattempted,
     x naccepted,
     x ninvolved,
     x nebfound,
     x ticknext,
     x ntick,
     x barrierev,
     x barrevev,
     x prefac,
     x isynthclass,
     x lfirstinstate,
     x lputyes,
     x lendput,
     x nofficialmax
     x )

c things that might trip us up
c     cdxyz - should this be stored, or allowed to differ and controlled by main prog?

c     iofficial - for now, we omit this altogether.  In the present
c                 version of the code, we make no attempt to determine equivalency
c                 through reording atoms or symmetry operations.  Thus, iofficial
c                 will always be the same as istate, and to include it here might be
c                 misleading.  Later, we will want to have both iofficial and information
c                 about how to transform to that other equivalent state.
c     temphigh - needed to complete the understanding of nattempted
c     templow - needed to complete the understanding of ticknext
c
c     iequiv - for now, this has been pulled out to prevent any confusion between
c               equiv states and official states. If there is an equivalent-energy
c               equiv state, it will be determined each time.
c     ntickmax - don' write and read this, as it is set as a parameter above
c
c
c    also - changed name of following vars:
c             before              now
c            timehigh           timehighnow
c            xyz1                xyzmin
c            nattempt           nattempted - now an array
c             eaev               barrierev - now an array

c begin real code
      use mod_mpi

      implicit real*8(a-h,o-z)
      character*(*) casename
      dimension xyz(3*mxatom) ! the hot config
      dimension xyzmin(3*mxatom) ! called xyz1 in tad2.f
      dimension xyz4(3*mxatom)
      dimension xyzneighs(lenxyzneighs)
      dimension xyzneighsads(lenxyzneighs)
      dimension rdimermodes(lenxyzneighs)
      dimension esaddle(nneighmax)
      dimension eneigh(nneighmax)
      dimension nattempted(nneighmax)
      dimension naccepted(nneighmax)
      dimension ninvolved(nneighmax)
      dimension nebfound(nneighmax)
      dimension ticknext(ntickmax*nneighmax)
      dimension ntick(nneighmax)
      dimension barrierev(nneighmax)
      dimension barrevev(nneighmax)
      dimension prefac(nneighmax)
      dimension isynthclass(nneighmax)

      dimension itype(mxatom) ! not dimensioned previously: 01/21/04
      dimension taxes(3)

      logical lfirstinstate, lputyes, lendput

      character*10 chri
      character*250 filnam

      iversion=20040248  ! 9/3/04 version, afv

      ! IF fma, Putstate by writing file...
      if((groupID.eq.fmanagerID).or.(execute_serial))then

c open unformatted file
      filnam=casename//'.'//chri(istate)//'.state'
      call chrpak(filnam,250,lfil)
      open(unit=19,file=filnam(1:lfil),
     x  form='unformatted',status='replace')
      rewind 19

      write(19)
     x istate,
     x natom,
     x nmove,
     x taxes(3),
     x ntype,
     x itype(natom),
     x barrierminev,
     x barrierminevss,
     x dimerbarrierminev,
     x bminlowerbound,
     x emin,
     x templow,
     x temphigh,
     x timehighprev,
     x timehighnow,
     x timelowprev,
     x tickprev,
     x ieventlast,
     x ineedblack,
     x ineedwarm,
     x ineedgauss,
     x ibadtransyes,
     x freqlogsmin,
     x nnegmin,
     x nprod,
     x nattempt,
     x nneigh,
     x (xyz(i),i=1,3*natom),    ! BEGIN ARRAYS HERE
     x (xyzmin(i),i=1,3*natom),
     x (xyz4(i),i=1,3*natom),
     x (((xyzneighs(3*natom*(k-1)+3*(i-1)+j),
     x    j=1,3),i=1,natom),k=1,nneigh+1),
     x (((xyzneighsads(3*natom*(k-1)+3*(i-1)+j),
     x    j=1,3),i=1,natom),k=1,nneigh+1),
     x (((rdimermodes(3*natom*(k-1)+3*(i-1)+j),
     x    j=1,3),i=1,natom),k=1,nneigh+1),
     x (esaddle(i),i=1,nneigh),
     x (eneigh(i),i=1,nneigh),
     x (nattempted(i),i=1,nneigh),
     x (naccepted(i),i=1,nneigh),
     x (ninvolved(i),i=1,nneigh),
     x (nebfound(i),i=1,nneigh),
     x (ticknext(i),i=1,ntickmax*nneigh),
     x (ntick(i),i=1,nneigh),
     x (barrierev(i),i=1,nneigh),
     x (barrevev(i),i=1,nneigh),
     x (prefac(i),i=1,nneigh),
     x (isynthclass(i),i=1,nneigh),
     x iversion

      close(unit=19)

      write(6,*) 'putstate - full state info stored: ',filnam(1:lfil)

      else ! Putstate by Sending all info in a message...

       call lp_putstate(
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

      endif !end if(groupID.eq.fmanagerID)then

      return
      end


      subroutine getstate(
     x casename,
     x istate,
     x natom,
     x nmove,
     x taxes,
     x ntype,
     x itype,
     x barrierminev,
     x barrierminevss,
     x dimerbarrierminev,
     x bminlowerbound,
     x emin,
     x templow,
     x temphigh,
     x timehighprev,
     x timehighnow,
     x timelowprev,
     x tickprev,
     x ieventlast,
     x ineedblack,
     x ineedwarm,
     x ineedgauss,
     x ibadtransyes,
     x freqlogsmin,
     x nnegmin,
     x nprod,
     x nattempt,
     x nneigh,
     x ntickmax,
     x nneighmax,
     x mxatom,
     x lenxyzneighs,
     x xyz,    ! BEGIN ARRAYS HERE
     x xyzmin,
     x xyz4,
     x xyzneighs,
     x xyzneighsads,
     x rdimermodes,
     x esaddle,
     x eneigh,
     x nattempted,
     x naccepted,
     x ninvolved,
     x nebfound,
     x ticknext,
     x ntick,
     x barrierev,
     x barrevev,
     x prefac,
     x isynthclass
     x )

      use mod_mpi

      implicit real*8(a-h,o-z)
      character*(*) casename
      dimension xyz(3,mxatom)
      dimension xyzmin(3,mxatom) ! called xyz1 in tad2.f
      dimension xyz4(3,mxatom)
      dimension xyzneighs(lenxyzneighs)
      dimension xyzneighsads(lenxyzneighs)
      dimension rdimermodes(lenxyzneighs)
      dimension esaddle(nneighmax)
      dimension eneigh(nneighmax)
      dimension nattempted(nneighmax)
      dimension naccepted(nneighmax)
      dimension ninvolved(nneighmax)
      dimension nebfound(nneighmax)
      dimension ticknext(ntickmax,nneighmax)
      dimension ntick(nneighmax)
      dimension barrierev(nneighmax)
      dimension barrevev(nneighmax)
      dimension prefac(nneighmax)
      dimension isynthclass(nneighmax)

      dimension itype(mxatom) ! not dimensioned previously: 01/21/04
      dimension taxes(3)

      character*10 chri
      character*250 filnam

      iversion=20040248  ! 9/3/04 version, afv

c open unformatted file
      filnam=casename//'.'//chri(istate)//'.state'
      call chrpak(filnam,250,lfil)

      open(unit=19,file=filnam(1:lfil),form='unformatted',
     x                                 status='unknown')
      rewind 19

      read(19)
     x istate,
     x natom,
     x nmove,
     x taxes(3),
     x ntype,
     x itype(natom),
     x barrierminev,
     x barrierminevss,
     x dimerbarrierminev,
     x bminlowerbound,
     x emin,
     x templow,
     x temphigh,
     x timehighprev,
     x timehighnow,
     x timelowprev,
     x tickprev,
     x ieventlast,
     x ineedblack,
     x ineedwarm,
     x ineedgauss,
     x ibadtransyes,
     x freqlogsmin,
     x nnegmin,
     x nprod,
     x nattempt,
     x nneigh,
     x ((xyz(j,i),j=1,3),i=1,natom),
     x ((xyzmin(j,i),j=1,3),i=1,natom),
     x ((xyz4(j,i),j=1,3),i=1,natom),
     x (((xyzneighs(3*natom*(k-1)+3*(i-1)+j),
     x    j=1,3),i=1,natom),k=1,nneigh+1),
     x (((xyzneighsads(3*natom*(k-1)+3*(i-1)+j),
     x    j=1,3),i=1,natom),k=1,nneigh+1),
     x (((rdimermodes(3*natom*(k-1)+3*(i-1)+j),
     x    j=1,3),i=1,natom),k=1,nneigh+1),
     x (esaddle(i),i=1,nneigh),
     x (eneigh(i),i=1,nneigh),
     x (nattempted(i),i=1,nneigh),
     x (naccepted(i),i=1,nneigh),
     x (ninvolved(i),i=1,nneigh),
     x (nebfound(i),i=1,nneigh),
     x ((ticknext(i,j),i=1,ntickmax),j=1,nneigh),
     x (ntick(i),i=1,nneigh),
     x (barrierev(i),i=1,nneigh),
     x (barrevev(i),i=1,nneigh),
     x (prefac(i),i=1,nneigh),
     x (isynthclass(i),i=1,nneigh),
     x iversionread

      close(unit=19)

      if(rank_l.eq.0)then
       write(6,*) 'getstate -full state info recovered: '!,filnam(1:lfil)
       write(6,*) 'getstate -state has ',nneigh,' neighbors'
      endif

      if((iversionread.ne.iversion).and.(rank_l.eq.0)) then
        write(6,*) 'ERROR - stopping, state file is wrong version'
        write(6,*) 'version found =',iversionread
        write(6,*) 'version expected =',iversion
        stop 'wrong version'
      end if

      return
      end

      subroutine ptadsetup(delta,alpha,pmissing,pmatters)
c     given delta, ptadsetup calculates alpha, pmissing, and pmatters for P-TAD
c      afv 10/24/02
c     delta is the uncertainty value
c     alpha is the ratio of low-T stop time to the low-T tick mark in P-TAD
c     finding alpha requires solving a transcendental equation, so we
c     approximate it for simplicity.
c     pmissing is the probability that an event is still missing
c     pmatters is the probability that it would matter
c     pmissing*pmatters=delta

c     given these quantities, the high-temperature stop time is given by
c       timehighstop = taupstar*(alpha*ticklow/taupstar)**(betahigh/betalow)
c       timelowswept = (1/alpha)*taupstar*(timehigh/taupstar)**(betalow/betahigh)
c     where taupstar is like taustart, but uses pmissing instead of the whole delta
c         taupstar = ln(1/pmissing)/prefacmin

c      Useful to note: for any delta smaller than about 1/100, pmissing is basically 0.37,
c      and, as a result, taupstar = 1/prefacmin  - i.e., the stop line pivots about prefacmin
c      [in regular tad, it pivots around prefamin/ln(1/delta) ],
c      and intersects the low-Temperature time line at the tickmark*alpha,

      implicit real*8(a-h,o-z)

c See afv notes of 10/24/02 to derive this approximation for alpha.
c It comes from noticing that ln(1/delta) ~= (d/dx)[x*ln(x)] evaluated at x=alpha+0.5

      alpha = 1.0d0/(exp(1.0d0)*delta) - 0.5d0
      pmissing = (alpha/(1.0d0+alpha))**alpha
      pmatters = 1.0d0/(1.0d0+alpha)
      deltacheck = pmissing*pmatters
      if(abs(deltacheck-delta)/delta.gt.0.1d0) then
        write(6,*) delta,alpha,pmissing,pmatters,deltacheck
        write(6,*) 'ptadsetup: deltacheck failed'
        stop 'ptadsetup: deltacheck failed'
      end if

c  correct for the slight error in alpha, so that pmissing*pmatters=delta exactly
      pmatters = delta/pmissing


      write(6,10) delta,alpha,deltacheck,pmissing,pmatters,
     x            pmissing*pmatters
   10 format(' ptadsetup:  ',/
     x       '    delta (supplied) =',1pd12.2,/,
     x       '    alpha (approximated) =',1pd12.2,/,
     x       '    deltacheck  =',1pd12.2,/,
     x       '    pmissing  =',1pd12.2,/,
     x       '    pmatters  =',1pd12.2,/,
     x       '    pmissing*pmatters  =',1pd12.2)

      return
      end

c======================================================================
c xxx new subroutines by bpu to clean up code above somewhat
c----------------------------------------------------------------------
      subroutine findminbarrier(nneigh,barrierev,barrierminev,barrevev
     +    ,ereverse)

      implicit real*8(a-h,o-z)
      dimension barrierev(nneigh)
      dimension barrevev(nneigh)

      barrierminev=1.d99

      do i=1,nneigh
         if(barrierev(i)-barrierminev.lt.-1.d-4) then
            if (barrevev(i).gt.ereverse)
     +          barrierminev=barrierev(i)
         endif
      enddo

      return
      end

c--------------------------------------------------------------------
      subroutine findminbarrierns(nneigh,barrierev,isynthclass,
     &    naccepted,barrierminevns)

      implicit real*8(a-h,o-z)
      dimension barrierev(nneigh)
      dimension isynthclass(nneigh)
      dimension naccepted(nneigh)

      barrierminevns=1.d99
      do i=1,nneigh
         if(isynthclass(i).eq.0.or.naccepted(i).eq.0) then
            if(barrierev(i)-barrierminevns.lt.-1.d-4)
     +          barrierminevns=barrierev(i)
         endif
      enddo

      if(barrierminevns.gt.1.d98)barrierminevns=0.d0

      return
      end
c-------------------------------------------------------------------------------
      subroutine findmaxbarriers(nneigh,barrierev,isynthclass,
     &    barriermaxevs)

      implicit real*8(a-h,o-z)
      dimension barrierev(nneigh)
      dimension isynthclass(nneigh)

      barriermaxevs=-1.d99
      do i=1,nneigh
         if(isynthclass(i).eq.1) then
            if(barrierev(i).gt.(barriermaxevs+1.d-4))
     +          barriermaxevs=barrierev(i)
         endif
      enddo

      if(barriermaxevs.lt.-1.d98)barriermaxevs=0.d0

      return
      end

c----------------------------------------------------------------------
      subroutine writebarrierentry(filnam,timehighnow,barrierev
     +    ,betalow,betahigh,istate,iofficial,ineigh,prefac,ninvolved
     +    ,barrevev,ctype,hyperdist,dist4)

      use mod_mpi
      implicit real*8(a-h,o-z)
      character*1 cdum
      character*2 ctype
      character*180 filnam

      if(rank_l.eq.0)then

      call chrpak(filnam,180,lfil)
      open(unit=12,file=filnam(1:lfil),form='formatted',status
     +    ='old')
      do ibpu=1,1000000         ! xxx note that this file could become huge
         read(12,232,end=233) cdum
      enddo
      stop 'uh oh - didnt read far enough in barriers file'
 232  format(a1)
 233  backspace 12
      timeloweventx=timehighnow*
     x    exp(barrierev*(betalow-betahigh))
      eendev=barrierev-barrevev
      write(12,234) istate,iofficial,ineigh,barrierev,barrevev,
     x    eendev,hyperdist,dist4,
     x    prefac,
     x    timehighnow,timeloweventx,ninvolved,ctype
      close(12)
 234  format(i10,i10,i9,0pf9.4,0pf9.4, f12.6,f12.5, f12.5,
     x    1pd16.3,1pd14.3,1pd14.3,2x,i7,a4)

      endif
      return
      end

c----------------------------------------------------------------------
      subroutine storechs(filnam,nhold,natom,taxes,springk,itype,
     +    xyzhold)

      use mod_mpi
      implicit real*8(a-h,o-z)

      dimension taxes(3),itype(natom),xyzhold(3*natom*nhold)
      character*180 filnam

      iunit=27
      lf=len(filnam)
      call chrpak(filnam,lf,lfil)

      if(rank_l.eq.0)then

      open(unit=iunit,file=filnam(1:lfil),form='formatted',status
     +         ='unknown')
      rewind iunit
      write(iunit,330) nhold,natom,0,1
      write(iunit,340) 'attempt chs from TAD'
      write(iunit,350) (taxes(k),k=1,3)
      write(iunit,351) springk
      write(iunit,352) (itype(i),0,i=1,natom)
 330  format(i4,i5,i5,i3)
 340  format(a50)
 350  format(3f20.8)
 351  format(f20.8)
 352  format(2i2)
      do i=1,nhold
         do j=1,natom
            write(iunit,350)
     +          (xyzhold(k+(j-1)*3+(i-1)*3*natom),k=1,3)
         enddo
      enddo
      rewind iunit
      close(unit=iunit)

      endif

      return
      end
c----------------------------------------------------------------------
      subroutine dovineyard(natom,nmove,xyz,itype,taxes,ietype,maxtyp
     +    ,rcut,amu,e,grad,nneg,nprod,freqlog,igeometry,freqlogmin
     +    ,nnegmin,prefac,ierr,ipr)

      use mod_mpi
      implicit real*8(a-h,o-z)

c      parameter(lwrk13=375 000)
      include 'wrksubs_tad2.h'
      common/wrkcm13/iwrk13,jwrk13,work13(lwrk13)

      dimension itype(natom),taxes(3),rcut(maxtyp,maxtyp),amu(maxtyp)
     +    ,grad(3*natom),xyz(3*natom)

      ierr=0

      nn=nmove*3*(nmove*3+1)/2
      if(nn.gt.lwrk13) then
        if(rank_l.eq.0)
     +    write(6,*)
     +    'STOP: hessian array too small to call nmodes '
     +    ,'- nmove,lwrk13 = '
     +    ,nmove,lwrk13
         stop 'hessian array too small for nmodes'
      end if

      call wrkclaimv(13,1,'tad2 main') ! for hessian curv array (work13)

      ip=-1
      ivec=0

      !gdelt=1.d-6
      !gdelt=1.d-5
      gdelt=1.d-4
      !gdelt=1.d-3
      !gdelt=1.d-2
      call nmodes(ip,ivec,natom,nmove,xyz,itype,taxes,ietype,maxtyp
     +    ,rcut,gdelt,amu,e,grad,work13,nneg,nprod,freqlog)
      call wrkrelease(13)
      if (igeometry.eq.0) then
         if(nneg.ne.0) then
           if(rank_l.eq.0)
     +       write(6,*) 'WARNING: nneg=',nneg,' for this minimum'
           ierr=1
         end if
         if ((ipr.ge.4).and.(rank_l.eq.0)) write(6,*)
     +       'vineyard: igeometry,nneg,freq: '
     +       ,igeometry,nneg,freqlog
      else if (igeometry.eq.1) then
         prefac=1d-99
         if(nneg.ne.1) then
           if(rank_l.eq.0)
     +        write(6,*) 'WARNING: nneg=',nneg,' for this saddle'
           ierr=1
         endif
         if ((ipr.ge.4).and.(rank_l.eq.0))
     +       write(6,*) 'vineyard: igeometry,nneg,freq: '
     +       ,igeometry,nneg,freqlog
         if(nneg.eq.1 .and. nnegmin.eq.0) then
            prefac=exp(freqlogmin-freqlog)
          else
           if(rank_l.eq.0)then
            if(ipr.ge.0) write(6,*) 'WARNING: problem with vineyard'
            if(ipr.ge.0) write(6,*) 'nneg for the saddle =',nneg
            if(ipr.ge.0) write(6,*) 'nneg for the minimum =',nnegmin
           endif
           ierr=1
          end if
      end if

      return
      end

      subroutine eneighmatch(natom,nmove,xyz,enow,taxes,eneigh,
     x                       xyzneighs,
     x                       transcrit, ninvolved,
     x                       nneigh,ecrit,rcrit,r4crit,ipr,     ineigh)
c compare to known neighbors with simple energy match and distance
c to final state.  The idea is that (if irecognize=9 above) we can
c save the trouble of finding the saddle, and keep the number of neighbors
c down, and get into synthetic mode quicker also.
c  August 30, 2006 -- added ninvolved into the check
c  August 31, 2006 -- added 4-dist comparison

c set ineigh to the neighbor that matches, if any

      implicit real*8(a-h,o-z)
      dimension xyz(3*natom),xyzneighs(3*natom*(nneigh+1))
      dimension eneigh(nneigh)
      dimension taxes(3)
      dimension ninvolved(nneigh)

      call numberinvolved(xyzneighs(1),xyz(1),nmove,transcrit,ninv)

c Note also:  could precompute these neighbor distances and keep in state file

c First compute hyperdistance to our unknown neighbor in xyz:
c Note: taking advantage of fact that neigh number 1 is the original minimum

      r2=0.0d0
      do j=1,3*nmove
        r2 = r2 + (xyz(j)-xyzneighs(j))**2
      end do
      dist=sqrt(r2)

c  sum up 4-distance
      r4=0.0d0
      jjj=0
      do j=1,nmove
      sumx=0.0d0
      do jj=1,3
      jjj=jjj+1
      sumx=sumx + (xyz(jjj)-xyzneighs(jjj))**2
      end do
      r4=r4+sumx**2
      end do
      dist4=sqrt(sqrt(r4))

c Now loop over neighbors

      nmatch=0
      do i=1,nneigh
        imatch=0
        if(ninvolved(i).eq.ninv .and. abs(eneigh(i)-enow).le.ecrit) then
          imatch=i
ccc          idx=3*natom*(i-1)
          idx=3*natom*i   ! this skips over the self-neighbor
          r2=0.0d0
          do j=1,3*nmove
            r2 = r2 + (xyzneighs(idx+j)-xyzneighs(j))**2
          end do
          disti=sqrt(r2)
          if(abs(disti-dist).gt.rcrit) imatch=0
          if(imatch.eq.i) then
c  sum up 4-distance
             r4=0.0d0
             jjj=0
             do j=1,nmove
                sumx=0.0d0
                do jj=1,3
                   jjj=jjj+1
                   sumx = sumx + (xyzneighs(idx+jjj)-xyzneighs(jjj))**2
                end do
                r4=r4+sumx**2
             end do
             dist4i=sqrt(sqrt(r4))
             if(abs(dist4i-dist4).gt.r4crit) imatch=0
          end if
        end if
c        write(6,2222) 'matchdebug', i,ninvolved(i),ninv,eneigh(i),enow,
c     x                      dist4,dist4i,imatch
c 2222   format(a10,i4,i4,i4,f12.5,f12.5,2x,2f12.5,i5)
c        write(6,*) 'matchdebugger', i,ecrit,rcrit,r4crit
        if(imatch.ne.0) then
          nmatch=nmatch+1
          if(nmatch.gt.1) then
            write(6,*) 'ERROR - eneighmatch vars: ',enow,ninv,dist
            write(6,*) 'ERROR - eneighmatch matched two neighbors:',
     x                  ineigh,imatch
            stop 'ERROR - eneighmatch matched two neighbors'
          end if
          ineigh=imatch
          distihold=disti
        end if
      end do

      if(ipr.ge.5 .and. ineigh.gt.0) then
          write(6,123) ineigh,enow,eneigh(ineigh),eneigh(ineigh)-enow,
     x                            dist,distihold,distihold-dist
  123     format(' eneighmatch: i,e,ei,de,r,ri,dr',i5,
     x     3f13.7,1x,3f10.4)
        end if

      if(ipr.ge.2) write(6,*) 'eneighmatch - ineigh =',ineigh

      return
      end

      subroutine pairsum(natom,xyz,r2sum)
c simple sum of all pair interactions, ignoring any possible periodic b.c.'s   1/8/07
      implicit real*8(a-h,o-z)
      dimension xyz(3,natom)

      r2sum=0.0d0
      do i=1,natom
      do j=1,i-1
      do k=1,3
c do we want the square root here?
      r2sum=r2sum + (xyz(k,i)-xyz(k,j))**2
      end do
      end do
      end do

      return
      end
