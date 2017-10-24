c potcong - general, self-contained potential converter - afv 10/98
c  the idea is to make a simple program for converting pots for new machines
      implicit real*8(a-h,o-z)
      character*50 command,filnam,filnamt

c the input file should look like this:

c binary_to_text
c ni.rho
c ni.phi
c ni.f
c ...
c 
c or
c 
c text_to_binary
c ni.rho
c ni.phi
c ni.f
c ...
c
c the "t" will be automatically appended for the text versions
c this way, only the first line of input needs to be altered to
c reverse the conversion
 

      read(5,10,end=200) command
   10 format(2a20)
c     allowed commands:  binary_to_text  or text_to_binary

      write(6,20) command
   20 format(' potcong program, action requested:',a20)

      do 100 i=1,10000
      read(5,10,end=200) filnam
      lf=lstchr(filnam,50)
      filnamt=filnam(1:lf)//'t'

      open(unit=51,file=filnam,form='unformatted',status='unknown')
      open(unit=52,file=filnamt,form='formatted',status='unknown')

      if(command.eq.'binary_to_text') then
        write(6,30) filnam(1:lf),filnamt(1:lf+1)
   30   format(' converting ',a,' to ',a)
        call newtrp
        call trptxo(51,52)
      else if(command.eq.'text_to_binary') then
        write(6,30) filnamt(1:lf+1),filnam(1:lf)
        call newtrp
        call trptxi(52,51)
      else
        stop 'command not recognized'
      end if

  100 continue

  200 continue
      end


      integer function lstchr(char,lenc)
c: returns the location of the last nonblank character in char(1:lenc)
c:
      character*(*) char
      integer lenc
      do 10 i=lenc,1,-1
      if(char(i:i).ne.' ') then
        lstchr=i
        return
      end if
   10 continue
      lstchr=0
      return
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
      if(mark(iu).ne.0) stop 'abort - file conflict in trptxi'
      rewind(itxt)
      read(itxt,10) n,x0(iu),xn(iu),dxi(iu),iflo(iu),ifhi(iu)
   10 format(i10,3d20.12,2i5)
      if(n+last.gt.20000) stop 'abort - redimension in trptxi'
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

      if(mark(iu).ne.0) stop 'abort - tried to read in trp twice'
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

