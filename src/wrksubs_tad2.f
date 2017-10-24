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


c:  wrksub_alone - for use with parmdg program, 3/98  -  scaled down
c: wrksubs - routines for interfacing to "wrkcom" work commons
c: Note that these calls are optional.  wrkcms can be claimed,
c: used, and released without having to call any of these routines.
c: However, monitoring with wrkmon will only work correctly if all
c: claims are made using wrkclaim and all releases are done with wrkrelease.
c:
      subroutine wrkinf(iwrt)
c: iwrt=0 = quick summary of iwrk values in work commons
c: iwrt=1 = verbose summary of work commons
c:
      implicit real*8(a-h,o-z)
      include 'wrksubs_tad2.h'
      common/wrkcm1/iwrk1,jwrk1,work(lwrk1)
      common/wrkcm2/iwrk2,jwrk2,work2(lwrk2)
      common/wrkcm3/iwrk3,jwrk3,work3(lwrk3)
      common/wrkcm4/iwrk4,jwrk4,work4(lwrk4)
      common/wrkcm5/iwrk5,jwrk5,work5(lwrk5)
      common/wrkcm6/iwrk6,jwrk6,work6(lwrk6)
      common/wrkcm7/iwrk7,jwrk7,work7(lwrk7)
      common/wrkcm8/iwrk8,jwrk8,work8(lwrk8)
      common/wrkcm9/iwrk9,jwrk9,work9(lwrk9)
      common/wrkcm10/iwrk10,jwrk10,work10(lwrk10)
      common/wrkcm11/iwrk11,jwrk11,work11(lwrk11)
      common/wrkcm12/iwrk12,jwrk12,work12(lwrk12)
      common/wrkcm13/iwrk13,jwrk13,work13(lwrk13)
      common/wrkcm14/iwrk14,jwrk14,work14(lwrk14)
      common/wrkcm15/iwrk15,jwrk15,work15(lwrk15)

      nwrk=15

      if(iwrt.eq.0) then
        write(6,10) iwrk1,iwrk2,iwrk3,iwrk4,iwrk5,iwrk6,iwrk7,iwrk8,
     x              iwrk9,iwrk10,iwrk11,iwrk12,iwrk13,iwrk14,iwrk15
   10   format('wrkinf:',15i5)
      else if(iwrt.ge.1) then
        write(6,*) 'n        length      iwrkn'
        write(6,20) 1,lwrk1,iwrk1
        write(6,20) 2,lwrk2,iwrk2
        write(6,20) 3,lwrk3,iwrk3
        write(6,20) 4,lwrk4,iwrk4
        write(6,20) 5,lwrk5,iwrk5
        write(6,20) 6,lwrk6,iwrk6
        write(6,20) 7,lwrk7,iwrk7
        write(6,20) 8,lwrk8,iwrk8
        write(6,20) 9,lwrk9,iwrk9
        write(6,20) 10,lwrk10,iwrk10
        write(6,20) 11,lwrk11,iwrk11
        write(6,20) 12,lwrk12,iwrk12
        write(6,20) 13,lwrk13,iwrk13
        write(6,20) 14,lwrk14,iwrk14
        write(6,20) 15,lwrk15,iwrk15
   20   format(3i10)
      end if

      return
      end

      subroutine wrkclaim(n,iset)
c: This routine claims  work common number n (wrkcm<n>)
c: If the common is available, it sets iwrk<n>=iset and returns.
c: If the common is already in use (iwrk<n>.ne.0), it complains and crashes.
c: The traceback should tell which routine was calling.
c: Note that calls to this routine can be mixed with direct wrkcm claims.
c:
      implicit real*8(a-h,o-z)
      dimension ihold(20)
      common/wrkcmm/imon
      include 'wrksubs_tad2.h'
      common/wrkcm1/iwrk1,jwrk1,work(lwrk1)
      common/wrkcm2/iwrk2,jwrk2,work2(lwrk2)
      common/wrkcm3/iwrk3,jwrk3,work3(lwrk3)
      common/wrkcm4/iwrk4,jwrk4,work4(lwrk4)
      common/wrkcm5/iwrk5,jwrk5,work5(lwrk5)
      common/wrkcm6/iwrk6,jwrk6,work6(lwrk6)
      common/wrkcm7/iwrk7,jwrk7,work7(lwrk7)
      common/wrkcm8/iwrk8,jwrk8,work8(lwrk8)
      common/wrkcm9/iwrk9,jwrk9,work9(lwrk9)
      common/wrkcm10/iwrk10,jwrk10,work10(lwrk10)
      common/wrkcm11/iwrk11,jwrk11,work11(lwrk11)
      common/wrkcm12/iwrk12,jwrk12,work12(lwrk12)
      common/wrkcm13/iwrk13,jwrk13,work13(lwrk13)
      common/wrkcm14/iwrk14,jwrk14,work14(lwrk14)
      common/wrkcm15/iwrk15,jwrk15,work15(lwrk15)

      nwrk=15

      ihold(1)=iwrk1
      ihold(2)=iwrk2
      ihold(3)=iwrk3
      ihold(4)=iwrk4
      ihold(5)=iwrk5
      ihold(6)=iwrk6
      ihold(7)=iwrk7
      ihold(8)=iwrk8
      ihold(9)=iwrk9
      ihold(10)=iwrk10
      ihold(11)=iwrk11
      ihold(12)=iwrk12
      ihold(13)=iwrk13
      ihold(14)=iwrk14
      ihold(15)=iwrk15

      if(ihold(n).le.0) then
        ihold(n)=iset
        if(imon.ne.0) call wrkinf(0)
        iwrk1=ihold(1)
        iwrk2=ihold(2)
        iwrk3=ihold(3)
        iwrk4=ihold(4)
        iwrk5=ihold(5)
        iwrk6=ihold(6)
        iwrk7=ihold(7)
        iwrk8=ihold(8)
        iwrk9=ihold(9)
        iwrk10=ihold(10)
        iwrk11=ihold(11)
        iwrk12=ihold(12)
        iwrk13=ihold(13)
        iwrk14=ihold(14)
        iwrk15=ihold(15)
      else
       call wrkinf(1)
       write(6,*) 'wrkclaim: refusing claim for work common number',n
       call uabort('wrkclaim: wrkcom claim refused$')
      end if

      return
      end

      subroutine wrkclaimv(n,iset,string)
c: verbose version of wrkclaim
c: This routine claims  work common number n (wrkcm<n>)
c: If the common is available, it sets iwrk<n>=iset and returns.
c: If the common is already in use (iwrk<n>.ne.0), it complains and crashes.
c: The traceback should tell which routine was calling.
c: Note that calls to this routine can be mixed with direct wrkcm claims.
c:
      implicit real*8(a-h,o-z)
      dimension ihold(20)
      common/wrkcmm/imon
      character*(*) string
    
      include 'wrksubs_tad2.h'
      common/wrkcm1/iwrk1,jwrk1,work(lwrk1)
      common/wrkcm2/iwrk2,jwrk2,work2(lwrk2)
      common/wrkcm3/iwrk3,jwrk3,work3(lwrk3)
      common/wrkcm4/iwrk4,jwrk4,work4(lwrk4)
      common/wrkcm5/iwrk5,jwrk5,work5(lwrk5)
      common/wrkcm6/iwrk6,jwrk6,work6(lwrk6)
      common/wrkcm7/iwrk7,jwrk7,work7(lwrk7)
      common/wrkcm8/iwrk8,jwrk8,work8(lwrk8)
      common/wrkcm9/iwrk9,jwrk9,work9(lwrk9)
      common/wrkcm10/iwrk10,jwrk10,work10(lwrk10)
      common/wrkcm11/iwrk11,jwrk11,work11(lwrk11)
      common/wrkcm12/iwrk12,jwrk12,work12(lwrk12)
      common/wrkcm13/iwrk13,jwrk13,work13(lwrk13)
      common/wrkcm14/iwrk14,jwrk14,work14(lwrk14)
      common/wrkcm15/iwrk15,jwrk15,work15(lwrk15)

      nwrk=15

      ihold(1)=iwrk1
      ihold(2)=iwrk2
      ihold(3)=iwrk3
      ihold(4)=iwrk4
      ihold(5)=iwrk5
      ihold(6)=iwrk6
      ihold(7)=iwrk7
      ihold(8)=iwrk8
      ihold(9)=iwrk9
      ihold(10)=iwrk10
      ihold(11)=iwrk11
      ihold(12)=iwrk12
      ihold(13)=iwrk13
      ihold(14)=iwrk14
      ihold(15)=iwrk15

      if(ihold(n).le.0) then
        ihold(n)=iset
        if(imon.ne.0) call wrkinf(0)
        iwrk1=ihold(1)
        iwrk2=ihold(2)
        iwrk3=ihold(3)
        iwrk4=ihold(4)
        iwrk5=ihold(5)
        iwrk6=ihold(6)
        iwrk7=ihold(7)
        iwrk8=ihold(8)
        iwrk9=ihold(9)
        iwrk10=ihold(10)
        iwrk11=ihold(11)
        iwrk12=ihold(12)
        iwrk13=ihold(13)
        iwrk14=ihold(14)
        iwrk15=ihold(15)
      else
       call wrkinf(1)
       write(6,*) 'wrkclaim: refusing claim for work common number',n
       write(6,*) string
       call uabort('wrkclaim: wrkcom claim refused$')
      end if

      return
      end

      subroutine wrkrelease(n)
c: This routine releases work common number n (wrkcm<n>).
c:
      implicit real*8(a-h,o-z)
      dimension ihold(20)
      common/wrkcmm/imon

      include 'wrksubs_tad2.h'
      common/wrkcm1/iwrk1,jwrk1,work(lwrk1)
      common/wrkcm2/iwrk2,jwrk2,work2(lwrk2)
      common/wrkcm3/iwrk3,jwrk3,work3(lwrk3)
      common/wrkcm4/iwrk4,jwrk4,work4(lwrk4)
      common/wrkcm5/iwrk5,jwrk5,work5(lwrk5)
      common/wrkcm6/iwrk6,jwrk6,work6(lwrk6)
      common/wrkcm7/iwrk7,jwrk7,work7(lwrk7)
      common/wrkcm8/iwrk8,jwrk8,work8(lwrk8)
      common/wrkcm9/iwrk9,jwrk9,work9(lwrk9)
      common/wrkcm10/iwrk10,jwrk10,work10(lwrk10)
      common/wrkcm11/iwrk11,jwrk11,work11(lwrk11)
      common/wrkcm12/iwrk12,jwrk12,work12(lwrk12)
      common/wrkcm13/iwrk13,jwrk13,work13(lwrk13)
      common/wrkcm14/iwrk14,jwrk14,work14(lwrk14)
      common/wrkcm15/iwrk15,jwrk15,work15(lwrk15)

      nwrk=15

      ihold(1)=iwrk1
      ihold(2)=iwrk2
      ihold(3)=iwrk3
      ihold(4)=iwrk4
      ihold(5)=iwrk5
      ihold(6)=iwrk6
      ihold(7)=iwrk7
      ihold(8)=iwrk8
      ihold(9)=iwrk9
      ihold(10)=iwrk10
      ihold(11)=iwrk11
      ihold(12)=iwrk12
      ihold(13)=iwrk13
      ihold(14)=iwrk14
      ihold(15)=iwrk15

      if(ihold(n).gt.0) then
        ihold(n)=0
        if(imon.ne.0) call wrkinf(0)
        iwrk1=ihold(1)
        iwrk2=ihold(2)
        iwrk3=ihold(3)
        iwrk4=ihold(4)
        iwrk5=ihold(5)
        iwrk6=ihold(6)
        iwrk7=ihold(7)
        iwrk8=ihold(8)
        iwrk9=ihold(9)
        iwrk10=ihold(10)
        iwrk11=ihold(11)
        iwrk12=ihold(12)
        iwrk13=ihold(13)
        iwrk14=ihold(14)
        iwrk15=ihold(15)
      else
        write(6,*) 'confused in wrkrelease - iwrk',n,'.le.0'
        call wrkinf(1)
        call uabort('wrkrelease: confused$')
      end if

      return
      end

      subroutine wrkmonon(iop)
c: This routine turns on monitoring.  Only works if all work 
c: common acquisitions and releases during monitoring
c: are accomplished using routines wrkclaim and wrkrelease
c: iop is not presently used.
c:
      implicit real*8(a-h,o-z)
      common/wrkcmm/imon

      imon=1
      return
      end

      subroutine wrkmonoff
c: This routine turns off monitoring. see routine workmonon
c:
      implicit real*8(a-h,o-z)
      common/wrkcmm/imon

      imon=0
      return
      end
