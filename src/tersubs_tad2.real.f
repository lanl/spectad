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


  
c RCS keywords for tersubs.f:
c     $Id: tersubs.f,v 2.0 1998/12/18 21:31:48 afv Exp $
c     $Revision: 2.0 $
c     $Author: afv $
c     $Locker:  $
c     $Date: 1998/12/18 21:31:48 $
c     $Header: /n/home3/afv/source/RCS/tersubs.f,v 2.0 1998/12/18 21:31:48 afv Exp $
*-----------------------------------------------------------------------
      subroutine prmst27(nvar,maxtyp,rcut,amu,azero)
*-----------------------------------------------------------------------

c  binary tersoff potental
      implicit real*8(a-h,o-z)
      dimension var(100)
      dimension rcut(maxtyp,maxtyp),amu(maxtyp),azero(maxtyp)
      character*40 varnam
      common/varcm2/varnam(100)
c  Following is the common for tersoff binary compounds.
c  Those arrays dimensioned (3) have ii,jj,ij values,
c  while those dimensioned (2) have i,j (no cross term).
c  The last line (c2,d2...) holds secondary parms calculated from primaries.
c  These are precomputed in prmst27 routine.
      common/ters2/biga(6),bigb(6),bigr(6),bigs(6),
     x                           dlam1(6),dlam2(6),
     x             beta(3),en(3),c(3),d(3),h(3),
     x             chi(6),omega(6),
     x             c2(3),d2(3),c2od2(3),betaen(3),hon(3)

c  note that cutoff distance is different for binary pot:
c     rcut = bigs

       n1=11
       n2=22
       ntot=n2+n1+6

       open(21,file='tersubs.parms',status='old')
       do i=1,39
          read(21,*) var(i)
       enddo
       close(21)

       write(6,*)
     +     'prmst27: si=1,c=2,ge=3; azero(ge) might not be 100% correct'

      if(nvar.eq.ntot) then
c  type 1 parms
        biga(1)=var(1)
        bigb(1)=var(2)
        bigr(1)=var(3)
        bigs(1)=var(4)
        dlam1(1)=var(5)
        dlam2(1)=var(6)
        beta(1)=var(7)
        en(1)=var(8)
        c(1)=var(9)
        d(1)=var(10)
        h(1)=var(11)
c  type 2 parms
        biga(3)=var(n1+1)
        bigb(3)=var(n1+2)
        bigr(3)=var(n1+3)
        bigs(3)=var(n1+4)
        dlam1(3)=var(n1+5)
        dlam2(3)=var(n1+6)
        beta(2)=var(n1+7)
        en(2)=var(n1+8)
        c(2)=var(n1+9)
        d(2)=var(n1+10)
        h(2)=var(n1+11)
c  type 3 parms
        biga(6)=var(n2+1)
        bigb(6)=var(n2+2)
        bigr(6)=var(n2+3)
        bigs(6)=var(n2+4)
        dlam1(6)=var(n2+5)
        dlam2(6)=var(n2+6)
        beta(3)=var(n2+7)
        en(3)=var(n2+8)
        c(3)=var(n2+9)
        d(3)=var(n2+10)
        h(3)=var(n2+11)
c  cross parms
        chi(2)=var(n2+n1+1)
        omega(2)=var(n2+n1+2)
        chi(4)=var(n2+n1+3)
        omega(4)=var(n2+n1+4)
        chi(5)=var(n2+n1+5)
        omega(5)=var(n2+n1+6)
        chi(1)=1.0d0
        chi(3)=1.0d0
        chi(6)=1.0d0
        omega(1)=1.0d0
        omega(3)=1.0d0
        omega(6)=1.0d0

c set up derived cross parms

        biga(2) = sqrt(biga(1)*biga(3))
        bigb(2) = sqrt(bigb(1)*bigb(3))
        bigr(2) = sqrt(bigr(1)*bigr(3))
        bigs(2) = sqrt(bigs(1)*bigs(3))
        dlam1(2) = 0.5d0*(dlam1(1)+dlam1(3))
        dlam2(2) = 0.5d0*(dlam2(1)+dlam2(3))

        biga(4) = sqrt(biga(1)*biga(6))
        bigb(4) = sqrt(bigb(1)*bigb(6))
        bigr(4) = sqrt(bigr(1)*bigr(6))
        bigs(4) = sqrt(bigs(1)*bigs(6))
        dlam1(4) = 0.5d0*(dlam1(1)+dlam1(6))
        dlam2(4) = 0.5d0*(dlam2(1)+dlam2(6))

        biga(5) = sqrt(biga(6)*biga(3))
        bigb(5) = sqrt(bigb(6)*bigb(3))
        bigr(5) = sqrt(bigr(6)*bigr(3))
        bigs(5) = sqrt(bigs(6)*bigs(3))
        dlam1(5) = 0.5d0*(dlam1(6)+dlam1(3))
        dlam2(5) = 0.5d0*(dlam2(6)+dlam2(3))

c  precompute some useful scalars also kept in tersoff common
        do 10 i=1,3
        c2(i) = c(i)*c(i)
        d2(i) = d(i)*d(i)
        c2od2(i)=c2(i)/d2(i)
        betaen(i)=beta(i)**en(i)
        hon(i)=0.5d0/en(i)
   10   continue

c  assign cutoffs
        rcut(1,1)=bigs(1)
        rcut(1,2)=bigs(2)
        rcut(2,1)=bigs(2)
        rcut(2,2)=bigs(3)
        rcut(1,3)=bigs(4)
        rcut(3,1)=bigs(4)
        rcut(2,3)=bigs(5)
        rcut(3,2)=bigs(5)
        rcut(3,3)=bigs(6)

        amu(1)=28.086d0
        amu(2)=12.011d0
        amu(3)=72.61d0
        azero(1)=5.431976d0
        azero(2)=3.5648d0
        azero(3)=5.54d0

        varnam(1)='A - type 1$'
        varnam(2)='B - type 1$'
        varnam(3)='R - type 1$'
        varnam(4)='S - type 1$'
        varnam(5)='lambda1 - type 1$'
        varnam(6)='lambda2 - type 1$'
        varnam(7)='beta$'
        varnam(8)='n - type 1$'
        varnam(9)='c - type 1$'
        varnam(10)='d - type 1$'
        varnam(11)='h - type 1$'

        varnam(n1+1)='A - type 2$'
        varnam(n1+2)='B - type 2$'
        varnam(n1+3)='R - type 2$'
        varnam(n1+4)='S - type 2$'
        varnam(n1+5)='lambda1 - type 2$'
        varnam(n1+6)='lambda2 - type 2$'
        varnam(n1+7)='beta - type 2$'
        varnam(n1+8)='n - type 2$'
        varnam(n1+9)='c - type 2$'
        varnam(n1+10)='d - type 2$'
        varnam(n1+11)='h - type 2$'

        varnam(n2+1)='A - type 3$'
        varnam(n2+2)='B - type 3$'
        varnam(n2+3)='R - type 3$'
        varnam(n2+4)='S - type 3$'
        varnam(n2+5)='lambda1 - type 3$'
        varnam(n2+6)='lambda2 - type 3$'
        varnam(n2+7)='beta - type 3$'
        varnam(n2+8)='n - type 3$'
        varnam(n2+9)='c - type 3$'
        varnam(n2+10)='d - type 3$'
        varnam(n2+11)='h - type 3$'

        varnam(n2+n1+1)='chi 1-2$'
        varnam(n2+n1+2)='omega 1-2$'
        varnam(n2+n1+3)='chi 1-3$'
        varnam(n2+n1+4)='omega 1-3$'
        varnam(n2+n1+5)='chi 2-3$'
        varnam(n2+n1+6)='omega 2-3$'
      else
        stop 'prmst27: bogus tersoff parms'
      end if

      return
      end

*-----------------------------------------------------------------------
      function tercut(r)
*-----------------------------------------------------------------------

c  tersoff cutoff, sin factor starts at bigr-bigd and goes to bir+bigd
      implicit real*8(a-h,o-z)
      common/tersoff/biga,bigb,bigr,bigd,
     x               dlam1,dlam2,dlam3,alpha,beta,en,c,d,h,
     x               c2,d2,c2od2,alpen,betaen,hon
      data pihalf/1.5707963267948965579989817d0/
      if(r.lt.bigr-bigd) then
        tercut=1.0d0
      else if(r.gt.bigr+bigd) then
        tercut=0.0d0
      else
        tercut=0.5d0*(1.0d0-sin(pihalf*(r-bigr)/bigd))
      end if
      return
      end

*-----------------------------------------------------------------------
      function trcut(r,bigr,bigd)
*-----------------------------------------------------------------------

c  tersoff cutoff, sin factor starts at bigr-bigd and goes to bir+bigd

      implicit real*8(a-h,o-z)
      data pihalf/1.5707963267948965579989817d0/
      if(r.lt.bigr-bigd) then
        trcut=1.0d0
      else if(r.gt.bigr+bigd) then
        trcut=0.0d0
      else
        trcut=0.5d0*(1.0d0-sin(pihalf*(r-bigr)/bigd))
      end if
      return
      end

*-----------------------------------------------------------------------
      function trcut2(r,bigr,bigs)
*-----------------------------------------------------------------------

c  new binary-style tersoff cutoff - binary form - bigr and bigs 
c  passed in.
c  note that bigr has different meaning than in single-component parms

      implicit real*8(a-h,o-z)
      data pi/3.1415926535897931159979635d0/
      if(r.lt.bigr) then
        trcut2=1.0d0
      else if(r.gt.bigs) then
        trcut2=0.0d0
      else
        trcut2=0.5d0*(1.0d0+cos(pi*(r-bigr)/(bigs-bigr)))
      end if
      return
      end

*-----------------------------------------------------------------------
      function trcutg(r)
*-----------------------------------------------------------------------

c  gradient of the tersoff cutoff, sin factor starts at bigr-bigd and
c   goes to bir+bigd

      implicit real*8(a-h,o-z)
c  Tersoff parameters
      common/tersoff/biga,bigb,bigr,bigd,
     x               dlam1,dlam2,dlam3,alpha,beta,en,c,d,h,
     x               c2,d2,c2od2,alpen,betaen,hon
      data pihalf/1.57079632679489655800d0/
      data pifort/0.78539816339744827900d0/
      if(r.lt.bigr-bigd) then
        trcutg=0.0d0
      else if(r.gt.bigr+bigd) then
        trcutg=0.0d0
      else
        trcutg=-pifort*cos(pihalf*(r-bigr)/bigd)/bigd
      end if
      return
      end

  
c RCS keywords for tersubs2.f:
c     $Id: tersubs2.f,v 2.0 1998/12/18 21:31:48 afv Exp $
c     $Revision: 2.0 $
c     $Author: afv $
c     $Locker:  $
c     $Date: 1998/12/18 21:31:48 $
c     $Header: /n/home3/afv/source/RCS/tersubs2.f,v 2.0 1998/12/18 21:31:48 afv Exp $
  
*-----------------------------------------------------------------------
      subroutine g27(natom, nmove, x, y, z, itype, tax,tay,taz,ietype
     +    ,maxtyp,rcut,energy, grad,hyperratio)
*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

      dimension x(natom), y(natom), z(natom), itype(natom)
      dimension grad(*)
      dimension rcut(maxtyp,maxtyp)

c  Voter work common
c   (includes Tersoff precomputed values)
      integer lwrk1,natmax,nnmax,lmax,iwrk1,jwrk1,nneigh,markn,jatomn
      real*8 dxn,dyn,dzn,rn,fill1
      real*8 b, bd, fr, frd, fa, fad, fc, fcd
      include "wrksubs_tad2.h"
      parameter (natmax=1000,nnmax=13,lmax=40000)
      common/wrkcm1/iwrk1,jwrk1,nneigh(natmax),markn(natmax),
     x         jatomn(lmax),dxn(lmax),dyn(lmax),dzn(lmax),rn(lmax),
     x         b(natmax,nnmax), bd(natmax,nnmax),
     x         fr(natmax,nnmax), frd(natmax,nnmax),
     x         fa(natmax,nnmax), fad(natmax,nnmax),
     x         fc(natmax,nnmax), fcd(natmax,nnmax),
     x         fill1(lwrk1-natmax-9*lmax/2-8*natmax*nnmax)


      common/print/iprint(6)
      common/potcom/rct(5,5),tx,ty,tz,ax,ay,az

c kludge for now - move these into call
      tx=tax
      ty=tay
      tz=taz
      ax=tax/2.0d0
      ay=tay/2.0d0
      az=taz/2.0d0
      do i=1,3
      do j=1,3
        rct(i,j)=rcut(i,j)
      end do
      end do
      
      ntype=imax(itype,natom)

      hyperratio=1.0d0
      
      if(iwrk1.gt.0) then
        stop 'g27: wrkcm1 conflict'
        return
      else
        iwrk1=1
      end if

      pad = 0.0d0

* rlistr does the neighbor stuff

      call rlistr (natom, 1, natom, x, y, z, itype, pad,
     x             nneigh, markn, jatomn, dxn, dyn, dzn, rn, lmax)

* etrpr2: load parameters and precompute some Tersoff 2-body-like funcs

      call etrpr2 (natom, x, y, z, itype)

      energy = 0.0d0
      do 10 ia=1,nmove

* etrs2g: Tersoff derivatives

        call etrs2g (ia, natom, itype, ei, eidx, eidy, eidz)

        grad(3*(ia-1)+1) = eidx
        grad(3*(ia-1)+2) = eidy
        grad(3*(ia-1)+3) = eidz
        energy = energy + ei
        if(iprint(1).ge.2) then
          write(6,20) ia,ei
   20     format(' energy of atom',i3,' from g27 routine =',1pd15.6)
        end if
   10 continue

      if(iprint(1).ge.1) then
        write(6,30) energy
   30   format(' total energy from g27 routine =',1pd15.6)
      end if

      iwrk1=0
      return
      end

*-----------------------------------------------------------------------
      subroutine etrpr2 (natom, x, y, z, itype)
*-----------------------------------------------------------------------

c  load parameters and precompute some Tersoff two-body-like funcs

      implicit none
      real*8 x, y, z, rj, xj, yj, zj, fcj, trcut2, trcu2g, xiij, 
     &       xien, xien1, bijarg, xk, yk, zk, rk, rdotr, rxr, costh, 
     &       hmcos, gden, g, fck, omik
      integer natom, itype, i, it, j, jt, ijpt, k, kt, ikpt, marki, jj,
     &        ijdx, kk, ikdx

c  input
c    natom         number of atoms
c    x(natom)      x coordinate of the of atoms
c    y(natom)      y coordinate of the of atoms
c    z(natom)      z coordinate of the of atoms
c    itype(natom)  type of the atoms
c  output (through common tersfun)
c    b(i,j)    Tersoff b func
c    bd(i,j)   non-Cartesian part of Tersoff b func derivative
c    fr(i,j)   Tersoff fr func
c    frd(i,j)  non-Cartesian part of Tersoff fr func derivative
c    fa(i,j)   Tersoff fa func
c    frd(i,j)  non-Cartesian part of Tersoff fa func derivative
c    fc(i,j)   Tersoff fc (cutoff) func
c    fcd(i,j)  non-Cartesian part of Tersoff fc func derivative

      dimension x(natom), y(natom), z(natom), itype(natom)

c  Following is the common for tersoff binary compounds.
c  Those arrays dimensioned (3) have ii,jj,ij values,
c  while those dimensioned (2) have i,j (no cross term).
c  The last line (c2,d2...) holds secondary parms calculated from 
c  primaries.
c  These are precomputed in prmst27 routine.

      real*8 biga, bigb, bigr, bigs, dlam1, dlam2, beta, en, c, d, h,
     x       chi, omega, c2, d2, c2od2, betaen, hon
      common/ters2/biga(6),bigb(6),bigr(6),bigs(6),
     x                           dlam1(6),dlam2(6),
     x             beta(3),en(3),c(3),d(3),h(3),
     x             chi(6),omega(6),
     x             c2(3),d2(3),c2od2(3),betaen(3),hon(3)

c ipoint array is for defining ij index given i and j
c     i    j     ij
c    ---  ---    ---
c     1    1      1
c     2    1      2
c     1    2      2
c     2    2      3
c     1    3      4
c     3    1      4
c     2    3      5
c     3    2      5
c     3    3      6

      integer ipoint
      dimension ipoint(3,3)
      data ipoint /1,2,4,2,3,5,4,5,6/

      save

c  Voter work common
c   (includes Tersoff precomputed values)

      integer lwrk1,natmax,nnmax,lmax,iwrk1,jwrk1,nneigh,markn,jatomn
      integer lwrk2,lwrk3,lwrk4,lwrk5,lwrk6,lwrk7,lwrk8,lwrk9,lwrk10
     +    ,lwrk11,lwrk12,lwrk13,lwrk14,lwrk15
      real*8 dxn,dyn,dzn,rn,fill1
      real*8 b, bd, fr, frd, fa, fad, fc, fcd
      include "wrksubs_tad2.h"
      parameter (natmax=1000,nnmax=13,lmax=40000)
      common/wrkcm1/iwrk1,jwrk1,nneigh(natmax),markn(natmax),
     x         jatomn(lmax),dxn(lmax),dyn(lmax),dzn(lmax),rn(lmax),
     x         b(natmax,nnmax), bd(natmax,nnmax),
     x         fr(natmax,nnmax), frd(natmax,nnmax),
     x         fa(natmax,nnmax), fad(natmax,nnmax),
     x         fc(natmax,nnmax), fcd(natmax,nnmax),
     x         fill1(lwrk1-natmax-9*lmax/2-8*natmax*nnmax)

        do 10 i=1,natom
        marki = markn(i)
        it = itype(i)
          do 10 jj=1,nneigh(i)
          if (nneigh(i) .gt. nnmax) stop 'etrpr2: increase nnmax'

          ijdx = marki + jj
          j = jatomn(ijdx)
          jt = itype(j)
          ijpt = ipoint(it,jt)
          xj = dxn(ijdx)
          yj = dyn(ijdx)
          zj = dzn(ijdx)
          rj = rn(ijdx)
          fcj = trcut2( rj, bigr(ijpt), bigs(ijpt) )
          xiij = 0.0d0

          do 20 kk=1,nneigh(i)

          ikdx = marki + kk
          k = jatomn(ikdx)
          kt = itype(k)
          ikpt = ipoint(it,kt)

          if (k .ne. j) then

            xk = dxn(ikdx)
            yk = dyn(ikdx)
            zk = dzn(ikdx)
            rk = rn(ikdx)

            rdotr = ( xj*xk + yj*yk +  zj*zk )
            rxr = ( rj*rk )
            costh = rdotr / rxr
            hmcos = h(it) - costh
            gden = d2(it) + hmcos**2
            g = 1.0d0 + c2od2(it) - c2(it)/gden
            fck = trcut2( rk, bigr(ikpt), bigs(ikpt) )
            omik = omega(ikpt)
            xiij = xiij + fck*g*omik

         end if

   20    continue

cjdk  Activate for pair potential only contribution
c          xiij = 0.0d0
cjdk
           if (xiij .eq. 0.0d0) then
             xien =  0.0d0
             xien1 =  0.0d0
           else
             xien =  xiij**en(it)
             xien1 =  xien/xiij
           end if
           bijarg =  1.0d0 + betaen(it) * xien
           b(i,jj) = chi(ijpt) / ( bijarg ** hon(it))
           bd(i,jj) = -( b(i,jj) / bijarg )*betaen(it)*0.5d0*xien1

           fr(i,jj)= biga(ijpt)*exp(-dlam1(ijpt)*rj)
           fa(i,jj)= -bigb(ijpt)*exp(-dlam2(ijpt)*rj)
           frd(i,jj)= -dlam1(ijpt)*fr(i,jj)
           fad(i,jj)= -dlam2(ijpt)*fa(i,jj)
           fc(i,jj) = trcut2( rj, bigr(ijpt), bigs(ijpt) )
           fcd(i,jj) = trcu2g( rj, bigr(ijpt), bigs(ijpt) )

   10    continue

      return
      end

*-----------------------------------------------------------------------
      subroutine etrs2g (ii, natom, itype, ei, eidx, eidy, eidz)
*-----------------------------------------------------------------------

c  Tersoff derivatives

      implicit none
      real*8
     > ei, eidx, eidy, eidz, etot, etijdx, etjidx, etjdx,
     > etijdy, etjidy, etjdy, etijdz, etjidz, etjdz, rj, xj, yj, zj,
     > fcj, rijdx, rjidx, fcjdx, rijdy, rjidy, fcjdy, rijdz, rjidz,
     > fcjdz, xiijdx, xijidx, etjkdx, xiijdy, xijidy,
     > etjkdy, xiijdz, xijidz, etjkdz,
     > xk, yk, zk, rk, rikdx, rikdy, rikdz
       real*8
     > rdotr, rxr, costh, hmcos, gden, g, fck, omik, trcg,
     > cosdx, gdx, fckdx, cosdy, gdy, fckdy,
     > cosdz, gdz, fckdz, xjk, yjk, zjk, rjk, fcjk, omjk, omij,
     > xijkdx, xijkdy, xijkdz, fajk, bij, frij, faij, bji,
     > e, eji
       real*8
     > bjkd, bijd, bjid, frijd, faijd, gd,
     > bjkdx, ejkdx, bijdx, bjidx, frijdx, faijdx,
     > bjkdy, ejkdy, bijdy, bjidy, frijdy, faijdy,
     > bjkdz, ejkdz, bijdz, bjidz, frijdz, faijdz,
     > eijdx, ejidx, eijdy, ejidy, eijdz, ejidz
      integer natom, itype, ii, it, j, jt, ijpt, k, kt, ikpt, jj, kk,
     > jkpt, marki, ijdx, ikdx, markj, jkdx, jfind, iofjj

      dimension itype(natom)

c  input
c    ii        atom number of atom ii
c    natom     number of atoms
c    itype     type of atoms
c  output
c    ei        energy of atom ii
c    eidx      x gradient of energy of atom ii
c    eidy      y gradient of energy of atom ii
c    eidz      z gradient of energy of atom ii

c  Voter work common
c   (includes Tersoff precomputed values)
      integer lwrk1,natmax,nnmax,lmax,iwrk1,jwrk1,nneigh,markn,jatomn
      integer lwrk2,lwrk3,lwrk4,lwrk5,lwrk6,lwrk7,lwrk8,lwrk9,lwrk10
     +    ,lwrk11,lwrk12,lwrk13,lwrk14,lwrk15
      real*8 dxn,dyn,dzn,rn,fill1
      real*8 b, bd, fr, frd, fa, fad, fc, fcd
      include "wrksubs_tad2.h"
      parameter (natmax=1000,nnmax=13,lmax=40000)
      common/wrkcm1/iwrk1,jwrk1,nneigh(natmax),markn(natmax),
     x         jatomn(lmax),dxn(lmax),dyn(lmax),dzn(lmax),rn(lmax),
     x         b(natmax,nnmax), bd(natmax,nnmax),
     x         fr(natmax,nnmax), frd(natmax,nnmax),
     x         fa(natmax,nnmax), fad(natmax,nnmax),
     x         fc(natmax,nnmax), fcd(natmax,nnmax),
     x         fill1(lwrk1-natmax-9*lmax/2-8*natmax*nnmax)

c  Following is the common for tersoff binary compounds.
c  Those arrays dimensioned (3) have ii,jj,ij values,
c  while those dimensioned (2) have i,j (no cross term).
c  The last line (c2,d2...) holds secondary parms calculated from primaries.
c  These are precomputed in prmst27 routine.
      real*8 biga, bigb, bigr, bigs, dlam1, dlam2, beta, en, c, d, h,
     x       chi, omega, c2, d2, c2od2, betaen, hon
      common/ters2/biga(6),bigb(6),bigr(6),bigs(6),
     x                           dlam1(6),dlam2(6),
     x             beta(3),en(3),c(3),d(3),h(3),
     x             chi(6),omega(6),
     x             c2(3),d2(3),c2od2(3),betaen(3),hon(3)

c ipoint array is for defining ij index given i and j
c     i    j     ij
c    ---  ---    ---
c     1    1      1
c     2    1      2
c     1    2      2
c     2    2      3
c     1    3      4
c     3    1      4
c     2    3      5
c     3    2      5
c     3    3      6

      integer ipoint
      dimension ipoint(3,3)
      data ipoint /1,2,4,2,3,5,4,5,6/

      save

        etot=0.0d0
        etijdx=0.0d0
        etijdy=0.0d0
        etijdz=0.0d0
        etjidx=0.0d0
        etjidy=0.0d0
        etjidz=0.0d0
        etjdx=0.0d0
        etjdy=0.0d0
        etjdz=0.0d0

        it = itype(ii)
        marki = markn(ii)
        do 10 j=1,nneigh(ii)
        ijdx = marki + j
        jj = jatomn(ijdx)
        jt = itype(jj)
        ijpt = ipoint(it,jt)
        xj = dxn(ijdx)
        yj = dyn(ijdx)
        zj = dzn(ijdx)
        rj = rn(ijdx)

        fcj = fc(ii,j)
        rijdx = xj/rj
        rijdy = yj/rj
        rijdz = zj/rj
        rjidx = xj/rj
        rjidy = yj/rj
        rjidz = zj/rj
        trcg = fcd(ii,j)
        fcjdx = rijdx*trcg
        fcjdy = rijdy*trcg
        fcjdz = rijdz*trcg

        xiijdx = 0.0d0
        xiijdy = 0.0d0
        xiijdz = 0.0d0
        xijidx = 0.0d0
        xijidy = 0.0d0
        xijidz = 0.0d0
        etjkdx = 0.0d0
        etjkdy = 0.0d0
        etjkdz = 0.0d0

        do 20 k=1,nneigh(ii)
          ikdx = marki + k
          kk = jatomn(ikdx)
          kt = itype(kk)
          ikpt = ipoint(it,kt)
          if (jj .ne. kk) then

c  first, calculate dvij/dxi, dvij/dyi, dvij/dzi

            xk = dxn(ikdx)
            yk = dyn(ikdx)
            zk = dzn(ikdx)
            rk = rn(ikdx)
            rikdx = xk/rk
            rikdy = yk/rk
            rikdz = zk/rk

            rdotr = ( xj*xk + yj*yk +  zj*zk )
            rxr = ( rj*rk )
            costh = rdotr / rxr
            hmcos = h(it) - costh
            gden = d2(it) + hmcos**2
            g = 1.0d0 + c2od2(it) - c2(it)/gden
            cosdx = (xk + xj) / rxr -
     >              costh * ( (rijdx/rj) + (rikdx/rk) )
            cosdy = (yk + yj) / rxr -
     >              costh * ( (rijdy/rj) + (rikdy/rk) )
            cosdz = (zk + zj) / rxr -
     >              costh * ( (rijdz/rj) + (rikdz/rk) )
            gd = (-2.0d0*c2(it)*hmcos/gden**2)
            gdx = gd*cosdx
            gdy = gd*cosdy
            gdz = gd*cosdz
            fck = fc(ii,k)
            trcg = fcd(ii,k)
            fckdx = rikdx*trcg
            fckdy = rikdy*trcg
            fckdz = rikdz*trcg
            omik = omega(ikpt)
            xiijdx = xiijdx + (fckdx*g + fck*gdx)*omik
            xiijdy = xiijdy + (fckdy*g + fck*gdy)*omik
            xiijdz = xiijdz + (fckdz*g + fck*gdz)*omik

        end if
   20   continue

c  next, calculate dvji/dxi, dvji/dyi, dvji/dzi

            markj = markn(jj)
            do 30 k=1,nneigh(jj)
              jkdx = markj + k
              kk = jatomn(jkdx)
              kt = itype(kk)
              jkpt = ipoint(jt,kt)
              if (kk .ne. ii) then
              xjk = -dxn(jkdx)
              yjk = -dyn(jkdx)
              zjk = -dzn(jkdx)
              rjk = rn(jkdx)

            rdotr = ( xj*xjk + yj*yjk +  zj*zjk )
            rxr = rj*rjk
            costh = rdotr / rxr
            hmcos = h(jt) - costh
            gden = d2(jt) + hmcos**2
            g = 1.0d0 + c2od2(jt) - c2(jt)/gden
            cosdx = ( xjk ) / rxr -
     >              costh * ( (rjidx/rj) )
            cosdy = ( yjk ) / rxr -
     >              costh * ( (rjidy/rj) )
            cosdz = ( zjk ) / rxr -
     >              costh * ( (rjidz/rj) )
            gd = (-2.0d0*c2(jt)*hmcos/gden**2)
            gdx = gd*cosdx
            gdy = gd*cosdy
            gdz = gd*cosdz
            fcjk = fc(jj,k)
            trcg = fcd(jj,k)
            omjk = omega(jkpt)

            xijidx = xijidx + fcjk*gdx*omjk
            xijidy = xijidy + fcjk*gdy*omjk
            xijidz = xijidz + fcjk*gdz*omjk

c  finally, calculate dvjk/dxi, dvjk/dyi, dvjk/dzi

c  Compute l=i term
            omij = omega(ijpt)
            xijkdx = (fcjdx*g + fcj*gdx)*omij
            xijkdy = (fcjdy*g + fcj*gdy)*omij
            xijkdz = (fcjdz*g + fcj*gdz)*omij

c  Compute l not= i terms

            bjkd = bd(jj,k)
            bjkdx = bjkd*xijkdx
            bjkdy = bjkd*xijkdy
            bjkdz = bjkd*xijkdz
            fajk= fa(jj,k)
            ejkdx= bjkdx*fajk
            ejkdy= bjkdy*fajk
            ejkdz= bjkdz*fajk
            etjkdx = etjkdx + fcjk*ejkdx
            etjkdy = etjkdy + fcjk*ejkdy
            etjkdz = etjkdz + fcjk*ejkdz

c  (kk .ne. ii) end if
          end if

   30   continue

        bij = b(ii,j)
        bijd = bd(ii,j)
        bijdx = bijd*xiijdx
        bijdy = bijd*xiijdy
        bijdz = bijd*xiijdz

c  Find the neighbor index, iofjj, of atom ii relative to atom jj
        do 35 jfind=1,nneigh(jj)
          if ( ii .eq. jatomn( markj + jfind ) ) then
            iofjj = jfind
            goto 36
          end if
   35   continue
        stop 'etrsfg: neighbor ii to jj not found'
   36   continue

        bji = b(jj,iofjj)
        bjid = bd(jj,iofjj)
        bjidx = bjid*xijidx
        bjidy = bjid*xijidy
        bjidz = bjid*xijidz

        frij= fr(ii,j)
        faij= fa(ii,j)
        frijd= frd(ii,j)
        frijdx= frijd*rijdx
        frijdy= frijd*rijdy
        frijdz= frijd*rijdz
        faijd= fad(ii,j)
        faijdx= faijd*rijdx
        faijdy= faijd*rijdy
        faijdz= faijd*rijdz
        e= frij + bij*faij
        eji= frij + bji*faij
        eijdx= bijdx*faij + frijdx + bij*faijdx
        eijdy= bijdy*faij + frijdy + bij*faijdy
        eijdz= bijdz*faij + frijdz + bij*faijdz
        ejidx= bjidx*faij + frijdx + bji*faijdx
        ejidy= bjidy*faij + frijdy + bji*faijdy
        ejidz= bjidz*faij + frijdz + bji*faijdz

        etot = etot + fcj*e
        etijdx = etijdx + fcjdx*e + fcj*eijdx
        etijdy = etijdy + fcjdy*e + fcj*eijdy
        etijdz = etijdz + fcjdz*e + fcj*eijdz
        etjidx = etjidx + fcjdx*eji + fcj*ejidx
        etjidy = etjidy + fcjdy*eji + fcj*ejidy
        etjidz = etjidz + fcjdz*eji + fcj*ejidz
        etjdx = etjdx + etjkdx
        etjdy = etjdy + etjkdy
        etjdz = etjdz + etjkdz
   10   continue

        ei = 0.5d0*etot/27.21d0
        etijdx = 0.5d0*etijdx/27.21d0
        etijdy = 0.5d0*etijdy/27.21d0
        etijdz = 0.5d0*etijdz/27.21d0
        etjidx = 0.5d0*etjidx/27.21d0
        etjidy = 0.5d0*etjidy/27.21d0
        etjidz = 0.5d0*etjidz/27.21d0
        etjdx = 0.5d0*etjdx/27.21d0
        etjdy = 0.5d0*etjdy/27.21d0
        etjdz = 0.5d0*etjdz/27.21d0
        eidx = etijdx + etjidx + etjdx
        eidy = etijdy + etjidy + etjdy
        eidz = etijdz + etjidz + etjdz
cc      eidx = 0.5d0*(etijdx + etjidx + etjdx)/27.21d0

      return
      end

*-----------------------------------------------------------------------
      function trcu2g(r,bigr,bigs)
*-----------------------------------------------------------------------

c  radial derivative of the
c  new binary-style tersoff cutoff - binary form - bigr and bigs 
c  passed in.
c  note that bigr has different meaning than in single-component parms

      implicit real*8(a-h,o-z)
      data pi/3.1415926535897931159979635d0/
      if(r.lt.bigr) then
        trcu2g=0.0d0
      else if(r.gt.bigs) then
        trcu2g=0.0d0
      else
        trcu2g=-0.5d0*pi*sin(pi*(r-bigr)/(bigs-bigr))/(bigs-bigr)
      end if
      return
      end
