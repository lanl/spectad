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


c: recovered August 1996 from a March 10, 1996 clsman_tot.f file -- may need to 
c: get rbw to recover from a tape for me.
c:  This is the trpsubs system of subroutines written by art voter 3/86.
c:                last changed by afv, 11/16/89
c:  please do not distribute, copy or use without permission from A. F. Voter
c:  The routines store, retrieve and use interpolation arrays, with
c:  minimal c:  effort required by the user.  The interpolation arrays
c:  are stored initially using either the maktrp, mkstrp or puttrp
c:  subroutines (see descriptions below), and usually accessed using
c:  subroutine trpfun.  A variety of other subroutines are also available
c:  for handling situations that may arise for the user.
c:
c:  To obtain a summary of this subroutine package:
c:    on a VMS VAX:
c:      $search  trpsubs.for/nohead  subroutine,function,entry,"c:","ci:"
c:    on a UNIX machine:
c:      user routines: %egrep    'SUBROUTINE|FUNCTION|ENTRY|c:'  trpsubs.f
c:      all routines:  %egrep -i 'subroutine|function|entry|c:|ci:'  trpsubs.f
c:
c:
      FUNCTION TRPFUN(x,ideriv,iu)
c:  Returns ideriv'th deriv. of f(x) using 2-point interp on array on unit iu.
c:  ideriv=0,1,2 returns appropriate *continuous* derivative.
c:  ideriv=11 returns f'(x)/x  ;    ideriv=12 returns  f''(x)/x**2 - f'(x)/x**3
c:  Automatically calls trpin at first occurence of this unit.
c:  Multiple interpolation arrays can be stored, referenced by unit number (iu)
c:  NOTE: currently, no other trpfn_ routine can handle ifhi=2 or iflo=2 cases.
c:  (ifhi=2 means linear extrapolation if x.gt.xn)
c:  (iflo=2 means linear extrapolation if x.lt.x1)
c:
      implicit real*8(a-h,o-z)
      common/holdtr/ a(40000),dum(80000)
      common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                       mark1(99),mark2(99),
     x                       iflo(99),ifhi(99),last
      common/trplvl/lvltrp

      if(lvltrp.eq.4) then
        trpfun=trpfn4(x,ideriv,iu)
        return
      else if(lvltrp.eq.3) then
        trpfun=trpfn3(x,ideriv,iu)
        return
      else if(lvltrp.eq.10) then
        trpfun=trpfns(x,ideriv,iu)
        return
      end if

      if(npt(iu).le.0) call trpin(iu)
      rm=dxi(iu)*(x-x0(iu))
      m=rm
c
c  ideriv=0 cases
c
      if(ideriv.gt.0) go to 20
      if(m.ge.npt(iu)) then
        if(ifhi(iu).eq.0) then
          trpfun=0.0d0
        else if(ifhi(iu).eq.1) then
          go to 990
        else if(ifhi(iu).eq.2) then
          mn=npt(iu)+mark(iu)
          slope=(a(mn)-a(mn-1))*dxi(iu)
          trpfun = a(mn) + slope*(x-xn(iu))
        end if
      else if(m.lt.1) then
        if(iflo(iu).eq.0) then
          if(m.eq.0) then
            trpfun = rm*a(1+mark(iu))
          else
            trpfun=0.0d0
          end if
        else if(iflo(iu).eq.1) then
          go to 990
        else if(iflo(iu).eq.2) then
          mn=2+mark(iu)
          slope=(a(mn)-a(mn-1))*dxi(iu)
          trpfun = a(mn) + slope*(x-x0(iu))
        end if
      else
        trpfun =a(m+mark(iu)) +
     x             (rm-float(m))*(a(m+1+mark(iu))-a(m+mark(iu)))
      end if
      return
c
c  ideriv>0 cases
c
   20 mm=mark(iu)+m
      if(m.ge.npt(iu)-1) then
        if(ifhi(iu).eq.0) then
          trpfun=0.0d0
        else if(ifhi(iu).eq.1) then
          go to 990
        else if(ifhi(iu).eq.2) then
          mn=npt(iu)+mark(iu)
          slope=(a(mn)-a(mn-1))*dxi(iu)
          if(ideriv.eq.1) trpfun=slope
          if(ideriv.eq.2) trpfun=0.0d0
          if(ideriv.eq.11) trpfun=slope/x
          if(ideriv.eq.12) trpfun=-slope/x**3
        end if
      else if(m.lt.2) then
        if(iflo(iu).eq.0) then
          trpfun=0.0d0
        else if(iflo(iu).eq.1) then
          go to 990
        else if(iflo(iu).eq.2) then
          mn=2+mark(iu)
          slope=(a(mn)-a(mn-1))*dxi(iu)
          if(ideriv.eq.1) trpfun=slope
          if(ideriv.eq.2) trpfun=0.0d0
          if(ideriv.eq.11) trpfun=slope/x
          if(ideriv.eq.12) trpfun=-slope/x**3
        end if
      else if(ideriv.eq.1) then
        gi=a(mm+1)-a(mm-1)
        gi1=a(mm+2)-a(mm)
        trpfun = ( gi + (rm-float(m))*(gi1-gi) )*dxi(iu)*0.5d0
      else if(ideriv.eq.2) then
        ci=a(mm-1)+a(mm+1)-2.d0*a(mm)
        ci1=a(mm)+a(mm+2)-2.d0*a(mm+1)
        trpfun = ( ci + (rm-float(m))*(ci1-ci) )*dxi(iu)*dxi(iu)
      else if(ideriv.eq.11) then
        gi=a(mm+1)-a(mm-1)
        gi1=a(mm+2)-a(mm)
        grad = ( gi + (rm-float(m))*(gi1-gi) )*dxi(iu)*0.5d0
        trpfun = grad/x
      else if(ideriv.eq.12) then
        gi=a(mm+1)-a(mm-1)
        gi1=a(mm+2)-a(mm)
        grad  = ( gi + (rm-float(m))*(gi1-gi) )*dxi(iu)*0.5d0
        ci=a(mm-1)+a(mm+1)-2.d0*a(mm)
        ci1=a(mm)+a(mm+2)-2.d0*a(mm+1)
        curv = ( ci + (rm-float(m))*(ci1-ci) )*dxi(iu)*dxi(iu)
        trpfun = (curv - grad/x)/(x*x)
      end if
      return
c
  990 write(6,991) x,iu,x0(iu),xn(iu)
  991 format(/'  abort -  r =',1pd12.4,' passed to trpfun(',i2,');'/
     x     '  this is out of range: start =',1pd10.2,' stop =',1pd10.2)
      callabortx('abort - out of range in trpfn')
      end

      FUNCTION TRPFNS(x,ideriv,iu)
c:  returns f(x) using smooth 4-point interp on array on unit iu
c:
      implicit real*8(a-h,o-z)
      common/holdtr/ a(40000),fp(40000),fpp(40000)
      common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                       mark1(99),mark2(99),
     x                       iflo(99),ifhi(99),last
      common/trplvl/lvltrp

      if(npt(iu).le.0) call trpin(iu)
      rm=dxi(iu)*(x-x0(iu))
      m=rm
      m=min(npt(iu),max(0,m))
      q=rm-float(m)
      q2=q*q
      q3=q*q2
      mm=mark(iu)+m
      if(ideriv.eq.0) then
        trpfns=a(mm) + q*fp(mm) + 0.5d0*q2*fpp(mm)
     x           + (fpp(mm+1)-fpp(mm))*(-1.5d0*q3+2.5d0*q2*q2-q2*q3)
      else if(ideriv.eq.1) then
        dydq=           fp(mm) + q*fpp(mm)
     x        + (fpp(mm+1)-fpp(mm))*(-4.5d0*q2+10.0d0*q3-5.0d0*q2*q2)
        trpfns=dydq*dxi(iu)
      else if(ideriv.eq.2) then
        d2ydq2=                  fpp(mm)
     x           + (fpp(mm+1)-fpp(mm))*(-9.0d0*q+30.0d0*q2-20.0d0*q3)
        trpfns=d2ydq2*dxi(iu)**2
      else
        write(6,*) 'trpfns - this ideriv not yet programmed'
      end if

      return
      end

      FUNCTION TRPFN3(x,ideriv,iu)
c:  returns f(x) using 3-point interp on array on unit iu (only ideriv=0).
c:  this is like TRPFUN, except that a more accurate 3-point interp is used
c:
      use mod_mpi
      implicit real*8(a-h,o-z)
      common/holdtr/ a(40000),dum(80000)
      common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                       mark1(99),mark2(99),
     x                       iflo(99),ifhi(99),last
      common/trplvl/lvltrp

      if(lvltrp.eq.4) then
        trpfn3=trpfn4(x,ideriv,iu)
        return
      else if(lvltrp.eq.10) then
        trpfun=trpfns(x,ideriv,iu)
        return
      end if
      if(npt(iu).le.0) call trpin(iu)
      rm=dxi(iu)*(x-x0(iu))
      m=rm
      if(ideriv.eq.1 .or. ideriv.eq.11) go to 20
      if(ideriv.gt.0) stop 'abort - this ideriv not allowed in trpfn3'
      if(m.ge.npt(iu)) then
        if(ifhi(iu).eq.0) then
          trpfn3=0.0d0
        else if(ifhi(iu).eq.1) then
          go to 990
        else if(ifhi(iu).eq.2) then
          mn=npt(iu)+mark(iu)
          slope=(a(mn)-a(mn-1))*dxi(iu)
          trpfn3 = a(mn) + slope*(x-xn(iu))
        end if
      else if(m.eq.1) then
c        q=(rm-float(m))
c        am=a(mark(iu)+m-1)
c        a0=a(mark(iu)+m)
c        a1=a(mark(iu)+m+1)
c        trpfn3=a0 + 0.5d0*q*( (a1+am-a0-a0)*q + (a1-am) )
        if(iflo(iu).eq.0) then
          am=0.0d0
          a0=a(mark(iu)+1)
          a1=a(mark(iu)+2)
          q=rm-1.0d0
          trpfn3=a0 + 0.5d0*q*( (a1+am-a0-a0)*q + (a1-am) )
        else if(iflo(iu).eq.1) then
          am=a(mark(iu)+1)
          a0=a(mark(iu)+2)
          a1=a(mark(iu)+3)
          q=rm-2.0d0
          trpfn3=a0 + 0.5d0*q*( (a1+am-a0-a0)*q + (a1-am) )
        else if(iflo(iu).eq.2) then
          trpfn3=a(mark(iu)+1)+(rm-1.0d0)*(a(mark(iu)+2)-a(mark(iu)+1))
        end if
      else if(m.eq.0) then
        if(iflo(iu).eq.0) then
          am=0.0d0
          a0=a(mark(iu)+1)
          a1=a(mark(iu)+2)
          q=rm-1.0d0
          trpfn3=a0 + 0.5d0*q*( (a1+am-a0-a0)*q + (a1-am) )
        else if(iflo(iu).eq.1) then
          go to 990
        else if(iflo(iu).eq.2) then
          trpfn3=a(mark(iu)+1)+(rm-1.0d0)*(a(mark(iu)+2)-a(mark(iu)+1))
        end if

ccold      else if(m.lt.2) then
ccold        if(m.eq.1) then
ccold       trpfn3 =a(m+mark(iu)) +
ccold  x             (rm-float(m))*(a(m+1+mark(iu))-a(m+mark(iu)))
ccold     else if(m.eq.0 .and. iflo(iu).eq.0) then
ccold       trpfn3 = rm*a(1+mark(iu))
ccold     else if(iflo(iu).eq.2) then
ccold       mn=2+mark(iu)
ccold       slope=(a(mn)-a(mn-1))*dxi(iu)
ccold       trpfn3 = a(mn) + slope*(x-x0(iu))
ccold     else
ccold       trpfn3=0.0d0
ccold       if(iflo(iu).ne.0) go to 990
ccold    end if

      else
        q=(rm-float(m))
        if((mark(iu)+m-1).le.0)then
          write(*,*) 'spawnID ',spawnID,', rank_l ',rank_l
     x      ,', about to access negative index in a: ',mark(iu)+m-1
     x      ,' with mark(iu),m = ',mark(iu),m
          call flushtad(6)
          call sleep(5)
        endif
        am=a(mark(iu)+m-1)
        a0=a(mark(iu)+m)
        a1=a(mark(iu)+m+1)
        trpfn3=a0 + 0.5d0*q*( (a1+am-a0-a0)*q + (a1-am) )
      end if
      return
c
   20 mm=mark(iu)+m
      if(ideriv.eq.2 .or. ideriv.eq.12) stop 'trpfn3: ideriv no good'
      if(m.ge.npt(iu)-1) then
        if(ifhi(iu).eq.0) then
          trpfn3=0.0d0
        else if(ifhi(iu).eq.1) then
          go to 990
        else if(ifhi(iu).eq.2) then
          a0=a(mm)
          a1=a(mm+1)
          slope=(a(mark(iu)+npt(iu))-a(mark(iu)+npt(iu)-1))*dxi(iu)
          if(ideriv.eq.1) trpfn3=slope
          if(ideriv.eq.11) trpfn3=slope/x
        end if
      else if(m.eq.1) then
        if(iflo(iu).eq.0) then
          ioff=0
cc          am=a(mm-1+ioff)
          am=0.0d0
          a0=a(mm+ioff)
          a1=a(mm+1+ioff)
          a2=a(mm+2+ioff)
          p = rm-float(m+ioff)
          gi=(p+0.5)*a1 - 2.0*p*a0 + (p-0.5)*am
          q=p-1.0d0
          gi1=(q+0.5)*a2 - 2.0*q*a1 + (q-0.5)*a0
          slope = ( gi + (rm-float(m))*(gi1-gi) )*dxi(iu)
          if(ideriv.eq.1) trpfn3=slope
          if(ideriv.eq.11) trpfn3=slope/x
        else if(iflo(iu).eq.1) then
          ioff=1
          am=a(mm-1+ioff)
          a0=a(mm+ioff)
          a1=a(mm+1+ioff)
          a2=a(mm+2+ioff)
          p = rm-float(m+ioff)
          gi=(p+0.5)*a1 - 2.0*p*a0 + (p-0.5)*am
          q=p-1.0d0
          gi1=(q+0.5)*a2 - 2.0*q*a1 + (q-0.5)*a0
          slope = ( gi + (rm-float(m))*(gi1-gi) )*dxi(iu)
          if(ideriv.eq.1) trpfn3=slope
          if(ideriv.eq.11) trpfn3=slope/x
        else if(iflo(iu).eq.2) then
          a0=a(mark(iu)+1)
          a1=a(mark(iu)+2)
          slope=(a1-a0)*dxi(iu)
          if(ideriv.eq.1) trpfn3=slope
          if(ideriv.eq.11) trpfn3=slope/x
        end if
      else if(m.eq.0) then
        if(iflo(iu).eq.0) then
          ioff=1
cc          am=a(mm-1+ioff)
          am=0.0d0
          a0=a(mm+ioff)
          a1=a(mm+1+ioff)
          a2=a(mm+2+ioff)
          p = rm-float(m+ioff)
          gi=(p+0.5)*a1 - 2.0*p*a0 + (p-0.5)*am
          q=p-1.0d0
          gi1=(q+0.5)*a2 - 2.0*q*a1 + (q-0.5)*a0
          slope = ( gi + (rm-float(m))*(gi1-gi) )*dxi(iu)
          if(ideriv.eq.1) trpfn3=slope
          if(ideriv.eq.11) trpfn3=slope/x
        else if(iflo(iu).eq.1) then
          go to 990
        else if(iflo(iu).eq.2) then
          a0=a(mark(iu)+1)
          a1=a(mark(iu)+2)
          slope=(a1-a0)*dxi(iu)
          if(ideriv.eq.1) trpfn3=slope
          if(ideriv.eq.11) trpfn3=slope/x
        end if
      else if(ideriv.eq.1) then
        p = rm-float(m)
        gi=(p+0.5)*a(mm+1) - 2.0*p*a(mm) + (p-0.5)*a(mm-1)
        q=p-1.0d0
        gi1=(q+0.5)*a(mm+2) - 2.0*q*a(mm+1) + (q-0.5)*a(mm)
        slope = ( gi + (rm-float(m))*(gi1-gi) )*dxi(iu)
        trpfn3=slope
      else if(ideriv.eq.11) then
        p = rm-float(m)
        gi=(p+0.5)*a(mm+1) - 2.0*p*a(mm) + (p-0.5)*a(mm-1)
        q=p-1.0d0
        gi1=(q+0.5)*a(mm+2) - 2.0*q*a(mm+1) + (q-0.5)*a(mm)
        slope = ( gi + (rm-float(m))*(gi1-gi) )*dxi(iu)
        trpfn3=slope/x
      else
        stop 'this ideriv not yet implemented in trpfn3'
      end if
      return
c
  990 write(6,991) x,iu,x0(iu),xn(iu)
  991 format(/'  abort -  r =',1pd12.4,' passed to trpfn3(',i2,');'/
     x     '  this is out of range: start =',1pd10.2,' stop =',1pd10.2)
      callabortx('abort - out of range in trpfn')
      end

      FUNCTION TRPFN4(x,ideriv,iu)
c:  returns f(x) using 4-point interp on array on unit iu (only ideriv=0).
c:  this is like TRPFUN, except that a more accurate 4-point interp is used
c:
      implicit real*8(a-h,o-z)
      common/holdtr/ a(40000),dum(80000)
      common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                       mark1(99),mark2(99),
     x                       iflo(99),ifhi(99),last
      if(npt(iu).le.0) call trpin(iu)
      rm=dxi(iu)*(x-x0(iu))
      m=rm
      if(m.ge.npt(iu)) then
        if(ifhi(iu).eq.0) then
          trpfn4=0.0d0
        else if(ifhi(iu).eq.1) then
          go to 990
        else if(ifhi(iu).eq.2) then
          mn=npt(iu)+mark(iu)
          slope=(a(mn)-a(mn-1))*dxi(iu)
          if(ideriv.eq.0) then
            trpfn4 = a(mn) + slope*(x-xn(iu))
          else if(ideriv.eq.1) then
            trpfn4 = slope
          else if(ideriv.eq.11) then
            trpfn4 = slope/x
          else if(ideriv.eq.2) then
            trpfn4 = 0.0d0
          else if(ideriv.eq.12) then
            trpfn4 = 0.0d0
          end if
        end if
      else if(m.lt.2) then
        if(m.eq.1) then
          trpfn4 =a(m+mark(iu)) +
     x             (rm-float(m))*(a(m+1+mark(iu))-a(m+mark(iu)))
        else if(m.eq.0 .and. iflo(iu).eq.0) then
          trpfn4 = rm*a(1+mark(iu))
        else if(iflo(iu).eq.2) then
          mn=2+mark(iu)
          slope=(a(mn)-a(mn-1))*dxi(iu)
          if(ideriv.eq.0) then
            trpfn4 = a(mn) + slope*(x-x0(iu))
          else if(ideriv.eq.1) then
            trpfn4 = slope
          else if(ideriv.eq.11) then
            trpfn4 = slope/x
          else if(ideriv.eq.2) then
            trpfn4 = 0.0d0
          else if(ideriv.eq.12) then
            trpfn4 = 0.0d0
          end if
        else
          trpfn4=0.0d0
          if(iflo(iu).ne.0) go to 990
        end if
      else
        q=(rm-float(m))
        am=a(mark(iu)+m-1)
        a0=a(mark(iu)+m)
        a1=a(mark(iu)+m+1)
        a2=a(mark(iu)+m+2)
        if(m.eq.npt(iu)-1) then
c         3 pt interp
          if(ideriv.gt.0) stop 'trpfn4 abort -some ideriv>0 dont work'
          trpfn4=a0 + 0.5d0*q*( (a1+am-a0-a0)*q + (a1-am) )
        else
c         4 pt interp
          if(ideriv.eq.0) then
            trpfn4 = (-q*(q-1.0d0)*(q-2.0d0)*am
     x              +(q**2-1.0d0)*(q-2.0d0)*a0*3.0d0
     x              -q*(q+1.0d0)*(q-2.0d0)*a1*3.0d0
     x              +q*(q**2-1.0d0)*a2  )/6.0d0
          else if(ideriv.eq.1) then
            q3 = 3.0d0*q**2
            trpfn4 = (-(q3 - 6.0d0*q + 2.0d0)*am
     x                +(q3 - 4.0d0*q -1.0d0)*a0*3.0d0
     x                -(q3 - 2.0d0*q -2.0d0)*a1*3.0d0
     x                +(q3 - 1.0d0)*a2  )*dxi(iu)/6.0d0
          else if(ideriv.eq.11) then
            q3 = 3.0d0*q**2
            slope = (-(q3 - 6.0d0*q + 2.0d0)*am
     x                +(q3 - 4.0d0*q -1.0d0)*a0*3.0d0
     x                -(q3 - 2.0d0*q -2.0d0)*a1*3.0d0
     x                +(q3 - 1.0d0)*a2  )*dxi(iu)/6.0d0
            trpfn4 = slope/x
          else if(ideriv.eq.2) then
            q6 = 6.0d0*q
            trpfn4 = ( -(q6 - 6.0d0)*am
     x                +(q6 - 4.0d0)*a0*3.0d0
     x                -(q6 - 2.0d0)*a1*3.0d0
     x                +q6*a2  )*dxi(iu)**2/6.0d0
          else if(ideriv.eq.12) then
            q3 = 3.0d0*q**2
            slope = (-(q3 - 6.0d0*q + 2.0d0)*am
     x                +(q3 - 4.0d0*q -1.0d0)*a0*3.0d0
     x                -(q3 - 2.0d0*q -2.0d0)*a1*3.0d0
     x                +(q3 - 1.0d0)*a2  )*dxi(iu)/6.0d0
            q6 = 6.0d0*q
            curv = ( -(q6 - 6.0d0)*am
     x                +(q6 - 4.0d0)*a0*3.0d0
     x                -(q6 - 2.0d0)*a1*3.0d0
     x                +q6*a2  )*dxi(iu)**2/6.0d0
            trpfn4 = (curv - slope/x)/x**2
          end if
        end if


      end if
      return


  990 write(6,991) x,iu,x0(iu),xn(iu)
  991 format(/'  abort -  r =',1pd12.4,' passed to trpfn4(',i2,');'/
     x     '  this is out of range: start =',1pd10.2,' stop =',1pd10.2)
      callabortx('abort - out of range in trpfn')
      end

      SUBROUTINE TRPVEC(n,x,y, ideriv,iu)
c:  like trpfun, but works on list of values in x, returns them in y
c:  x and y are n-long vectors
c:  ideriv,iu are scalars
c:
c  note that we require an extra zero point at a(mark(iu)+npt(iu)+1) for iflhi=0
c  trpin is responsible for this

      implicit real*8(a-h,o-z)
      dimension x(n),y(n)
      common/holdtr/ a(40000),dum(80000)
      common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                       mark1(99),mark2(99),
     x                       iflo(99),ifhi(99),last
      dimension mvec(100),delvec(100)

      if(npt(iu).le.0) call trpin(iu)
      if(ideriv.ge.1 .and. mark1(iu).le.0) call trpin1(iu)

      do 10 i=1,n
      rm=dxi(iu)*(x(i)-x0(iu))
      m=rm
      mvec(i)=max(min(m,npt(iu)),0)
      delvec(i)=rm-float(mvec(i))
   10 continue

      if(ifhi(iu).ne.0) then
        mmax=0
        do 20 i=1,n
        mmax=max(mvec(i),mmax)
   20   continue
        if(mmax.ge.npt(iu)) then
          write(6,*) 'trpvec range violation, here are the x values:'
          write(6,*) (x(i),i=1,n)
          write(6,*) ' maximum allowed for unit',iu,' =',xn(iu)
          stop 'trpvec - max x range exceeded'
        end if
      end if

      if(iflo(iu).ne.0) then
        mmin=npt(iu)
        do 30 i=1,n
        mmin=min(mvec(i),mmin)
   30   continue
        if(mmin.le.0) then
          write(6,*) 'trpvec range violation, here are the x values:'
          write(6,*) (x(i),i=1,n)
          write(6,*) ' x0 allowed for unit',iu,' =',x0(iu)
          stop 'trpvec - min x range exceeded'
        end if
      end if

      if(ideriv.eq.0) then
        do 100 i=1,n
        m=mvec(i)+mark(iu)
        y(i)=a(m)+delvec(i)*(a(m+1)-a(m))
  100   continue
      else if(ideriv.eq.1) then
        do 101 i=1,n
        m=mvec(i)+mark1(iu)
        y(i)=a(m)+delvec(i)*(a(m+1)-a(m))
  101   continue
      else if(ideriv.eq.11) then
        do 111 i=1,n
        m=mvec(i)+mark1(iu)
        y(i)=a(m)+delvec(i)*(a(m+1)-a(m))
        y(i)=y(i)/x(i)
  111   continue
      else
        write(6,*) 'trpvec - illegal ideriv passed in:',ideriv
        stop 'trpvec - illegal ideriv'
      end if

      return

      end


      SUBROUTINE PUTTRP(a,npt,x0,xn,iflo,ifhi,iunit)
c:  stores the interp array in a(npt) on unit iunit; see below for iflo,ifhi.
c:  stores a .trp file.  note convention for x0,xn,a:
c:  internal storage area is not used
c:    a(i)=f(x0+dx*i), a(npt)=f(xn), where dx=(xn-x0)/npt
c:    iflo=0 means f(x)=0 for x<x0+dx ;  set iflo=1 if not allowed below x0+dx
c:    ifhi=0 means f(x)=0 for x>xn ;    set ifhi=1  if not allowed above xn
c:    ifhi=2 means use linear extrapolation for x>xn
c:
      implicit real*8(a-h,o-z)
      dimension a(npt)
      common/trpcom/mute
      dxi=float(npt)/(xn-x0)
c
      rewind iunit
      write(iunit) npt,x0,xn,dxi,iflo,ifhi
      write(iunit) a
      rewind iunit
c
      if(mute.eq.0) write(6,60) iunit,npt
   60 format(' a .trp file has been stored on unit',i3,', npt =',i6)
      return
      end

      subroutine GETTRP(a,npt,x0,xn,iflo,ifhi,iunit)
ci:  Retrieves the contents of the .trp file from unit iunit.
ci:  NOTE: User normally will NOT need to use this routine.
ci:  The array does not get added to internal storage area.
ci:  Calling arguments are identical to those of PUTTRP.
ci:  User must allocate enough space for real*8 array a(npt)
ci:  This routine is useful if one wishes to change ifhi or ihlo.
ci:  e.g. call gettrp(...); ifhi=2; call puttrp(...)
ci:
      implicit real*8(a-h,o-z)
      dimension a(*)
c
      rewind iunit
      read(iunit) npt,x0,xn,dxi,iflo,ifhi
      read(iunit) (a(i),i=1,npt)
      rewind iunit
c
      return
      end

      subroutine TRPIN(iu)
ci:  Internal routine. Reads an interpolation array from .trp file
ci:  and stores it in memory, with appropriate pointers so that the
ci:  other TRP routines can access it.
ci:  Called automatically by trpfun at first occurence of this unit.
ci:  Multiple trp arrays can be read in, referenced by unit number.
ci:  Storage is completely independent of the trp-generating routines,
ci:  since they write directly to disk, without storing internally.
ci:

      implicit real*8(a-h,o-z)
      common/holdtr/ a(40000),fp(40000),fpp(40000)
      common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                       mark1(99),mark2(99),
     x                       iflo(99),ifhi(99),last
      common/trpcom/mute
c stmt funcs:
c  e.g.,  pslope(i,j)=slope, at point j, of parabola through
c                                     points i,i+1,i+2 in array a
      pval(i,j)=(a(i+2)+a(i)-2.0d0*a(i+1))*float(j-i-1)**2/2.0d0
     x                  + 0.5d0*(a(i+2)-a(i))*float(j-i-1) + a(i+1)
      pslope(i,j)=(a(i+2)+a(i)-2.0d0*a(i+1))*float(j-i-1)
     x                  + 0.5d0*(a(i+2)-a(i))
      pcurv(i)=(a(i+2)+a(i)-2.0d0*a(i+1))

      if(mark(iu).ne.0) callabortx('abort - tried to read in trp twice')
      rewind iu
      read(iu) n,x0(iu),xn(iu),dxi(iu),iflo(iu),ifhi(iu)
      
      if(n+last.gt.40000) stop 'abort - redimension in trpin'
      npt(iu)=n
      mark(iu)=last+1
      last=last+n+1
      m=mark(iu)
      read(iu) (a(i),i=m+1,m+n)
      rewind iu

c  note points put in a() array off each end of interp array - see trpvec
      a(m)=0.0d0
      a(m+n+1)=0.0d0
      if(ifhi(iu).eq.2) a(m+n+1) = a(m+n) + a(m+n)-a(m+n-1)
c  this following stmt may be slightly wrong - see trpvec
      if(iflo(iu).eq.2) a(m) = a(m+1) + a(m+1)-a(m+2)

      x1=x0(iu)+1.0d0/dxi(iu)

c  setup for newer trpfns routine - smooth, differentiable interp

c fill fp,fpp
      do  10 ii=2,npt(iu)-1
      i=mark(iu)+ii
      fp(i)=(a(i+1)-a(i-1))/2.0d0
      fpp(i)=a(i+1)+a(i-1)-2.0d0*a(i)
   10 continue

      mrk=mark(iu)
c low end
      if(iflo(iu).eq.0) then
        a(mrk)=0.0d0
        fp(mrk)=pslope(mrk+0,mrk+0)
        fpp(mrk)=pcurv(mrk+0)
        fp(mrk+1)=(a(mrk+2)-a(mrk+0))/2.0d0
        fpp(mrk+1)=(a(mrk+2)+a(mrk+0)-2.0d0*a(mrk+1))
      else if(iflo(iu).eq.1) then
        a(mrk)=1.d50
        fp(mrk)=0.0d0
        fpp(mrk)=0.0d0
        fp(mrk+1)=pslope(mrk+1,mrk+1)
        fpp(mrk+1)=pcurv(mrk+1)
      else if(iflo(iu).eq.2) then
        fp(mrk+1)=pslope(mrk+1,mrk+1)
        fpp(mrk+1)=pcurv(mrk+1)
        a(mrk)=a(mrk+1)-1.0d0*fp(mrk+1)
        fp(mrk)=fp(mrk+1)
      end if
c high end
      nn=mark(iu)+npt(iu)
      if(ifhi(iu).eq.0) then
        a(nn)=0.0d0
        fp(nn)=0.0d0
        fpp(nn)=0.0d0
        fpp(nn+1)=0.0d0
      else if(ifhi(iu).eq.1) then
        a(nn)=1.d60
        fp(nn)=0.0d0
        fpp(nn)=0.0d0
        fpp(nn+1)=0.0d0
      else if(ifhi(iu).eq.2) then
        fp(nn)=pslope(nn-2,nn)
        fpp(nn)=0.0d0
        fpp(nn+1)=0.0d0
      end if
      if(mute.eq.0) write(6,50) iu,n,iflo(iu),ifhi(iu),x0(iu),xn(iu),
     x                                   x1,a(m+1),xn(iu),a(m+n)
   50 format(' trp in: unit=',i2,' npt=',i5,' iflo=',i2,' ifhi=',i2,
     x                                 '  x0,xn:',1p2d12.4/
     x       '         x1,y1,xn,yn:',1p4d10.2)
      return
      end

      subroutine TRPIN1(iu)
ci:  Internal routine. Generates the 1st deriv array marked by mark1
ci:  see trpvec and trpin
ci:  Called automatically at first occurence derivs for this unit.
ci:
      implicit real*8(a-h,o-z)
      common/holdtr/ a(40000),dum(80000)
      common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                       mark1(99),mark2(99),
     x                       iflo(99),ifhi(99),last
      common/trpcom/mute
      if(mark1(iu).ne.0) callabortx('trpin1 - tried to do twice')
      n=npt(iu)

      mark1(iu)=last+1
      last=last+n+1
      m1=mark1(iu)

      do 10 i=1,n
      a(m1+i) = (a(m1+i+1)-a(m1+i-1))*dxi(iu)*0.5d0
   10 continue
      if(iflo(iu).eq.0) a(m1+1)=0.0d0
      if(ifhi(iu).eq.0) a(m1+n)=0.0d0

c  note elements in a() array off each end of interp array - see trpvec
c   this should make iflo,ifhi = 0 or 2 work correctly
      a(m1)=a(m1+1)
      a(m1+n+1)=a(m1+n)

      write(6,*) 'trpin1 done for unit', iu

      return
      end

      subroutine TRPOUT(iu)
ci:  Expert/internal routine. Writes an interpolation array from
ci:  the internal trpsubs storage area to disk.  This is normally
ci:  not used, because any array in the internal storage area was
ci:  read from disk.  However, it is useful if TRPMOD has been used
ci:  to modify the internal storage area.
ci:  Be sure to open the file and close it afterwards.
ci:
      implicit real*8(a-h,o-z)
      common/holdtr/ a(40000),dum(80000)
      common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                       mark1(99),mark2(99),
     x                       iflo(99),ifhi(99),last
      common/trpcom/mute

      m=mark(iu)+1
      n=npt(iu)
      call puttrp(a(m),n,x0(iu),xn(iu),iflo(iu),ifhi(iu),iu)
      end

      SUBROUTINE NEWTRP
c:  clears all interp arrays from memory. use before changing an interp array.
c:  this routine is called to read purge all current .trp arrays
c:  so that new (different) ones may be read in by trpin
c:  user should call this routine if puttrp is called again w/ same unit
c:
      implicit real*8(a-h,o-z)
      common/holdtr/ a(40000),dum(80000)
      common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                       mark1(99),mark2(99),
     x                       iflo(99),ifhi(99),last
      last=0
      do 10 i=1,99
      npt(i)=0
      mark(i)=0
      mark1(i)=0
      mark2(i)=0
      x0(i)=0.0d0
      xn(i)=0.0d0
      dxi(i)=0.0d0
      iflo(i)=0
      ifhi(i)=0
   10 continue
      return
      end

      SUBROUTINE TRPDAT(iu,npoint,xx0,xx1,xxn,yy1,yyn,ifxlo,ifxhi)
c:  this routine can be called by user to get specs on this trp functi0n
c:  xx1 is the smallest usable x value; yy1 = trpfun(iu;xx1;0)
c:  (x can be smaller than xx1 if iflo=0)
c:  xxn = end of range (x must be <xxn unless iflo=0); ynn=trpfun(ynn)
c:
      implicit real*8(a-h,o-z)
      common/holdtr/ a(40000),dum(80000)
      common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                       mark1(99),mark2(99),
     x                       iflo(99),ifhi(99),last
      if(npt(iu).le.0) call trpin(iu)
      npoint=npt(iu)
      xx0=x0(iu)
      xx1=x0(iu) + 1.0d0/dxi(iu)
      xxn=xn(iu)
      yy1=a(mark(iu)+1)
      yyn=a(mark(iu)+npoint)
      ifxlo=iflo(iu)
      ifxhi=ifhi(iu)
      return
      end

      SUBROUTINE MAKTRP(a,npt,xstart,xstop,extfn,iflo,ifhi,iunit)
c:  makes & stores trp array over external func extfn. a(npt) is work space.
c:  makes and stores a real*8 interpolation array as a .trp file
c:  using functi0n extfn.
c:       (often it is easier for user to set up array directly,
c:        and call puttrp.  see below for demonstration of how.)
c:  user passes in npt,xstart,xstop,extfn,ifhi,iflo,iunit & work array a
c:  see code below for how to use routine trpfun to access the functi0n
c:  NOTE: interpolation is only valid in the range  xstart+dx < x < xstop
c:
      implicit real*8(a-h,o-z)
      dimension a(npt)
      common/trpcom/mute
      external extfn
      if(xstart.ge.xstop) callabortx('abort - invalid range in maktrp')
c
c  fill a(npt)
c
      dx=(xstop-xstart)/float(npt)
      do 20 i=1,npt
      x= xstart + float(i)*dx
      a(i)=extfn(x)
   20 continue
c
c  store this as a .trp file
c
      call puttrp(a,npt,xstart,xstop,iflo,ifhi,iunit)
c
c  test it out; find maximum deviation
c
      devmx=0.0d0
      do 40 i=1,npt-1
      x=xstart + float(i)*dx
      do 30 j=1,3
      x=x+dx*0.25d0
      atrue=extfn(x)
      atrp=trpfun(x,0,iunit)
      if(abs(atrue-atrp).gt.abs(devmx)) then
        xdevmx=x
        devmx=atrp-atrue
        idevmx=i
      end if
   30 continue
   40 continue
c
      if(mute.eq.0) write(6,50) devmx,xdevmx,idevmx,idevmx+1
   50 format('  max deviation =',1pd10.2,' at x =',1pd10.2,
     x                  ' (between points',i6,' and',i6,')')
      return
      end

      SUBROUTINE MKSTRP(a,npt,xstart,xstop,extfn,iflo,npow,iunit,cmax)
c:  smoothing version of routine maktrp. sets f(x),f'(x)=0 at xstop endpoint.
c:  the difference between this routine and maktrp is the smoothing done here
c:  makes and stores a real*8 interpolation array as a .trp file
c:  using functi0n extfn.
c:  smoothing:   zero functi0n and 1st derivative at x=xstop
c:               npow controls shape of smoothing
c:               npow.ge.2 recommended - see notes of 2/26/86 for details
c:
      implicit real*8(a-h,o-z)
      dimension a(npt)
      common/trpcom/mute
      external extfn
      if(npow.le.0) callabortx('abort - invalid npow in mkstrp')
      if(xstart.ge.xstop) callabortx('abort - invalid range in maktrp')
c
c  fill a(npt)
c
      dx=(xstop-xstart)/float(npt)
      yn=extfn(xstart+float(npt)*dx)
      ynm1=extfn(xstart+float(npt-1)*dx)
      ypn=(yn-ynm1)/dx
      if(yn*ypn.gt.0.0d0) then
        if(cmax.lt.0.0d0) then
           cmax=1000.0d0
           return
        else
          callabortx('abort- fn not decaying in mkstrp')
        end if
      end if
      temp=xstop*ypn/float(npow)
      do 20 i=1,npt
      x= xstart + float(i)*dx
      a(i)=extfn(x) - yn + temp*(1.0d0-(x/xstop)**npow)
   20 continue
c
c  return maximum deviation for user reference
c
      cmax= -yn + temp
c
c  store this as a .trp file
c
      ifhi=0
      call puttrp(a,npt,xstart,xstop,iflo,ifhi,iunit)
      return
      end

      SUBROUTINE MKSTRP2(a,npt,xstart,xstop,extfn,iflo,npow,iunit,cmax)
c:  smoothing version of routine maktrp. sets f(x),f'(x)=0 at xstop endpoint.
c:  this version zeros f''(xstop) also
c:  the difference between this routine and maktrp is the smoothing done here
c:  makes and stores a real*8 interpolation array as a .trp file
c:  using functi0n extfn.
c:  smoothing:   zero functi0n and 1st derivative at x=xstop
c:               npow controls shape of smoothing
c:               npow.ge.2 recommended - see notes of 2/26/86 for details
c:
      implicit real*8(a-h,o-z)
      dimension a(npt)
      common/trpcom/mute
      external extfn
      if(npow.le.0) callabortx('abort - invalid npow in mkstrp')
      if(xstart.ge.xstop) callabortx('abort - invalid range in maktrp')
c
c  fill a(npt)
c
      dx=(xstop-xstart)/float(npt)
      yn=extfn(xstart+float(npt)*dx)
      ynm1=extfn(xstart+float(npt-1)*dx)
      ynp1=extfn(xstart+float(npt+1)*dx)
cc      ypn=(yn-ynm1)/dx
cc use centered derivs instead
      ypn=(ynp1-ynm1)/(2.0d0*dx)
      yppn=(ynp1+ynm1-yn-yn)/dx**2
      if(yn*ypn.gt.0.0d0) then
        if(cmax.lt.0.0d0) then
           cmax=1000.0d0
           return
        else
          callabortx('abort- fn not decaying in mkstrp')
        end if
      end if
cc      temp=xstop*ypn/float(npow)
      pow=float(npow)
      do 20 i=1,npt
      x= xstart + float(i)*dx
      bigx=(xstop/pow)*(1.0d0-(x/xstop)**npow)
      yraw=extfn(x)
      a(i)=yraw - yn +
     x        bigx*ypn - 0.5d0*bigx**2*(yppn-ypn*(pow-1.0d0)/xstop)
   20 continue

c  return maximum deviation for user reference
c  (assumes last point)

      cmax= a(npt) - yraw

c  store this as a .trp file

      ifhi=0
      call puttrp(a,npt,xstart,xstop,iflo,ifhi,iunit)
      return
      end

      SUBROUTINE MKSTRP3(a,npt,xstart,xstop,extfn,iflo,npow,iunit,cmax)
c:  smoothing version of routine maktrp. sets f(x),f'(x)=0 at xstop endpoint.
c:  this version zeros f''(xstop) also, and f'''  (8/96)
c:  the difference between this routine and maktrp is the smoothing done here
c:  makes and stores a real*8 interpolation array as a .trp file
c:  using functi0n extfn.
c:  smoothing:   zero functi0n and 1st derivative at x=xstop
c:               npow controls shape of smoothing
c:               npow.ge.2 recommended - see notes of 2/26/86 for details
c:
      implicit real*8(a-h,o-z)
      dimension a(npt)
      common/trpcom/mute
      external extfn
      if(npow.le.0) stop 'abort - invalid npow in mkstrp3'
      if(xstart.ge.xstop) stop 'abort- invalid range in mkstrp3'
c
c  fill a(npt)
c
      dx=(xstop-xstart)/float(npt)
      yn=extfn(xstart+float(npt)*dx)
      ynm2=extfn(xstart+float(npt-2)*dx)
      ynm1=extfn(xstart+float(npt-1)*dx)
      ynp1=extfn(xstart+float(npt+1)*dx)
      ynp2=extfn(xstart+float(npt+2)*dx)
      ypn=(ynp1-ynm1)/(2.0d0*dx)
      yppn=(ynp1+ynm1-yn-yn)/dx**2
      ypppn=(ynp2-2.d0*ynp1+2.d0*ynm1-ynm2)/(2.d0*dx**3)   ! 5 pt 3rd deriv.
      if(yn*ypn.gt.0.0d0) then
        if(cmax.lt.0.0d0) then
           cmax=1000.0d0
           return
        else
          stop 'abort- fn not decaying in mkstrp3'
        end if
      end if
cc      temp=xstop*ypn/float(npow)
      pow=float(npow)
        bigxn=0.0d0
        bigxpn=-1.0d0
        bigxppn=bigxpn*(pow-1.d0)/xstop
        bigxpppn=bigxppn*(pow-2.d0)/xstop
        biga=-(ypn)/bigxpn
        bigb=-(yppn + bigxppn*biga)/bigxpn**2
        bigc=-(ypppn+bigxpppn*biga + 3.d0*bigxppn*bigxpn*bigb)/bigxpn**3
      do 20 i=1,npt
      x= xstart + float(i)*dx
      bigx=(xstop/pow)*(1.0d0-(x/xstop)**npow)
      yraw=extfn(x)
      a(i)=yraw - yn + biga*bigx + bigb*bigx**2/2.d0 + bigc*bigx**3/6.d0

c      a(i)=yraw - yn + bigx*ypn 
c     x  - 0.5d0*bigx**2*(yppn-ypn*(pow-1.0d0)/xstop)
c     x  + (1.d0/6.d0)*bigx**3*( ypppn 
c     x              - 3.d0*(pow-1.d0)/xstop*(yppn+ypn*(pow-1.d0)/xstop)
c     x                       - ypn*(pow-1.d0)*(pow-2.d0)/xstop**2 )
   20 continue

c  return maximum deviation for user reference
c  (assumes last point)

      cmax= a(npt) - yraw

c  store this as a .trp file

      ifhi=0
      call puttrp(a,npt,xstart,xstop,iflo,ifhi,iunit)
      return
      end

      SUBROUTINE TRPSPK(iu,npoint,iunit)
c:  using trp array on iu, writes npoint points to iunit (x,y,y',y'')(1p4d15.6).
c:  store the .trp functi0n iu over npoint points on iunit for speakeasy plots
c:
      implicit real*8(a-h,o-z)
      common/trpcom/mute
      call trpdat(iu,npt,x0,x1,xn,y1,yn,iflo,ifhi)
      dx=(xn-x1)/float(npoint+2)
      write(iunit,10) iu,npt,x1,xn,y1,yn,iflo,ihi
   10 format(' iu,npt,x1,xn,y1,yn,iflo,ifhi:',2i6,1p4d10.2,2i3)
c
      do 30 i=1,npoint
      x=x1+float(i)*dx
      if(i.eq.1) x=x+1.d-10
      y=trpfun(x,0,iu)
      yp=trpfun(x,1,iu)
      ypp=trpfun(x,2,iu)
      write(iunit,20) x,y,yp,ypp
   20 format(1p4d15.6)
   30 continue
      if(ifhi.eq.0) then
        r=x1+float(npoint+1)*dx
        write(iunit,20) r,trpfun(r,0,iu),trpfun(r,1,iu)
        r=x1+float(npoint+2)*dx
        write(iunit,20) r,trpfun(r,0,iu)
      end if
c
      if(mute.eq.0) write(6,60) iunit,npt
   60 format(' a file has been stored on unit',i3,', npoint =',i6)
      return
      end

      SUBROUTINE TRPTXO(iu,itxt)
c:  converts trp file on iu to a special (portable) text file on unit itxt.
c:  this routine is useful for converting from machine to machine.
c:  also see:  trptxi(itxt,iu) - converts text file back to standard
c:  unformatted file that the trp routines expect.
c:
      implicit real*8(a-h,o-z)
      common/holdtr/ a(40000),dum(80000)
      common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                       mark1(99),mark2(99),
     x                       iflo(99),ifhi(99),last
c
      if(npt(iu).le.0) call trpin(iu)
c
      write(itxt,10) npt(iu),x0(iu),xn(iu),dxi(iu),iflo(iu),ifhi(iu)
   10 format(i10,3d20.12,2i5)
      write(itxt,20) (a(i),i=mark(iu)+1,mark(iu)+npt(iu))
   20 format(4d20.12)
c
      write(6,30) iu,itxt,npt(iu),x0(iu),xn(iu),iflo(iu),ifhi(iu)
   30 format(' trptxo:  trp file from unit ',i3,' is now in text format'
     x      ,' on unit',i3/
     x       ' (npt,x0,xn,iflo,ifhi: ',i6,1p2d12.4,2i4,' )')
      return
      end

      SUBROUTINE TRPTXI(itxt,iu)
c:  converts special text file on unit itxt (see trptxo) back to trp file on iu.
c:  this is the complement of routine trptxi(iu,itxt)
c:
      implicit real*8(a-h,o-z)
      common/holdtr/ a(40000),dum(80000)
      common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                       mark1(99),mark2(99),
     x                       iflo(99),ifhi(99),last
c
      if(mark(iu).ne.0) callabortx('abort - file conflict in trptxi')
      rewind(itxt)
      read(itxt,10) n,x0(iu),xn(iu),dxi(iu),iflo(iu),ifhi(iu)
     
   10 format(i10,3d20.12,2i5)
      if(n+last.gt.20000) callabortx('abort - redimension in trptxi')
      npt(iu)=n
      mark(iu)=last
      last=last+n
      m=mark(iu)
      read(itxt,20) (a(i),i=m+1,m+n)
   20 format(4d20.12)
c
      rewind itxt
      rewind iu
      write(iu) n,x0(iu),xn(iu),dxi(iu),iflo(iu),ifhi(iu)
      write(iu) (a(i),i=m+1,m+n)
      rewind iu
c
      write(6,30) itxt,iu,npt(iu),x0(iu),xn(iu),iflo(iu),ifhi(iu)
   30 format(' trptxi: text format file from unit ',i3,
     x       ' is now standard on unit',i3/
     x       ' (npt,x0,xn,iflo,ifhi: ',i6,1p2d12.4,2i4,' )')
      return
      end

      subroutine TRPMOD(ipt,yi,iu)
ci:  Expert-level routine. Allows user to modify one element of
ci:  the interpolation array for unit iu, setting y(ipt)=yi.
ci:  Note that the the array is not stored on disk, so this
ci:  is a fast operation.  Use TRPOUT to make the modified
ci:  array go to disk. (also see TRPZER)
ci:
      implicit real*8(a-h,o-z)
      common/holdtr/ a(40000),dum(80000)
      common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                       mark1(99),mark2(99),
     x                       iflo(99),ifhi(99),last
c
      if(npt(iu).le.0) call trpin(iu)
      if(ipt.lt.1 .or. ipt.gt.npt(iu)) stop 'trpmod abort - illegal ipt'
      a(mark(iu)+ipt)=yi
      return
      end

      subroutine TRPZER(n,x0x,xnx,iflox,ifhix,iu)
ci:  Expert-level routine. Creates, in INTERNAL AREA ONLY,
ci:  a zeroed-out interpolation array.  Use this in
ci:  conjunction with TRPMOD to create a trpfile in place
ci:  without using storage in the calling routine.
ci:  To write it to disk, use TRPOUT.
ci:
      implicit real*8(a-h,o-z)
      common/holdtr/ a(40000),dum(80000)
      common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                       mark1(99),mark2(99),
     x                       iflo(99),ifhi(99),last
      common/trpcom/mute

c  put array at end if iu is new or if npt is greater than before
      if(n.gt.npt(iu)) then
        mark(iu)=last
        last=last+n
        if(last.gt.40000) stop 'abort - redimension in trpin'
      end if

      npt(iu)=n
      x0(iu)=x0x
      xn(iu)=xnx
      dxi(iu)=float(n)/(xnx-x0x)
      iflo(iu)=iflox
      ifhi(iu)=ifhix

      m=mark(iu)
      do 10 i=1,n
   10 a(m+i)=0.0d0

      return
      end

      subroutine TRPSET(array,n,x0x,xnx,iflox,ifhix,iu)
ci:  Expert-level routine. Sets, in INTERNAL AREA ONLY,
ci:  all values interpolation array.
ci:  Routines TRPZER and TRPMOD also
ci:  modify the interp array in memory only, but are
ci:  less general.
ci:  To write it to disk, use TRPOUT.
ci:
      implicit real*8(a-h,o-z)
      dimension array(n)
      common/holdtr/ a(40000),dum(80000)
      common/umarkt/ npt(99),mark(99),x0(99),xn(99),dxi(99),
     x                       mark1(99),mark2(99),
     x                       iflo(99),ifhi(99),last
      common/trpcom/mute


c  put array at end if iu is new or if npt is greater than before
      if(n.gt.npt(iu)) then
        mark(iu)=last
        last=last+n
        if(last.gt.40000) stop 'abort - redimension in trpin'
      end if

      npt(iu)=n
      x0(iu)=x0x
      xn(iu)=xnx
      dxi(iu)=float(n)/(xnx-x0x)
      iflo(iu)=iflox
      ifhi(iu)=ifhix

      m=mark(iu)
      do 10 i=1,n
   10 a(m+i)=array(i)

      return
      end
      
      subroutine abortx(string)
      character*(*) string
      write(6,*) string
      stop 'aborting'
      end
