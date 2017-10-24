c --------------------nik's conj. grad. code---------------

      subroutine cgstep(natom,nmove,xyz,itype,taxes,ietype,maxtyp,rcut
     $     ,hyperratio,ein,eout,dE,gmag,initcg,cgfac,g1)
c     if iop=0 conjugate gradient minimization
c     if iop=1 gradient minimization
      implicit real*8(a-h,o-z)
c      parameter (mxatom= 10000)
      include 'parameters.h'
      dimension rcut(maxtyp,maxtyp),itype(natom),xyz(3*natom),taxes(3)

c g1 is a work array for cg, but is needed to save info between calls
c am using the v array otherwise used for quickmin in this case
      dimension g1(3*natom)

c     work1 is some work array for storing things that some of these
C     routines need.
      dimension work1(4*mxatom),p1(mxatom,3),h2(3*mxatom),g2(3*mxatom)

c      write(6,*) 'entered cg routine'

      iop=0
      ifun=0
      tol=0.000001d0
      tol=-1.0d0
      pmin=0.1d0
      pmin=0.0d0
      pmax=2.0d0

c      pmin=0.0d0
c      pmax=0.5d0

c     ----------

      istat=0
      etmp=0.0d0
      nstp1=0
      nstp2=0

      do i=1,natom
         p1(i,1)=xyz((i-1)*3+1)
         p1(i,2)=xyz((i-1)*3+2)
         p1(i,3)=xyz((i-1)*3+3)
      enddo

      if (iop.eq.0.and.initcg.eq.1) then
         call gcalc(natom,nmove,xyz,itype, taxes,ietype,maxtyp,rcut,
     x        em,g1, hyperratio)

         sum2=0.0d0
         do i=1,3*nmove
            g1(i)=-g1(i)
            sum2=sum2+g1(i)**2
         enddo
         do i=3*nmove+1,3*natom
            g1(i)=0.0d0
         enddo

         g1len=sqrt(sum2)

c unit length prevents errors with boxing when called twice on same cluster
         do i=1,3*nmove
            g1(i)=g1(i)/g1len
         enddo

         if (ifun.eq.0) then
            call linmind(mxatom,natom,nmove,itype,taxes,ietype,maxtyp
     $           ,rcut,hyperratio,p1,g1,work1,drmax,nstp1,etmp,eout,dE
     $           ,cgfac,pmin,pmax,tol)
         else if (ifun.eq.1) then
            call linmine(mxatom,natom,nmove,itype,taxes,ietype,maxtyp
     $           ,rcut,hyperratio,p1,g1,drmax,nstp1,-g1len,work1,etmp
     $           ,eout,dE,cgfac,pmin,pmax,tol)
         endif
         
         do i=1,3*nmove
            g1(i)=g1(i)*g1len
         enddo
      endif

      call gcalc(natom,nmove,xyz,itype, taxes,ietype,maxtyp,rcut,
     x     em,g2, hyperratio)

      gmag=sqrt(vecdot(g2,g2,3*nmove))

      do i=1,3*nmove
         g2(i)=-g2(i)
      enddo
      do i=3*nmove+1,3*natom
         g2(i)=0.0d0
      enddo
      
      if (iop.eq.0) then
         ggn=0
         gn=0
         do i=1,3*nmove
            ggn=ggn + (g2(i)-g1(i))*g2(i)
            gn=gn + g1(i)*g1(i)
         enddo
         gamma=ggn/gn
      else if (iop.eq.1) then
         gamma=0.0d0
      endif

      sum2=0.0d0
      do i=1,3*nmove
         h2(i)=g2(i) + gamma*g1(i)
         sum2=sum2+h2(i)*h2(i)
      enddo
      h2len=sqrt(sum2)

      do i=1,3*nmove
         h2(i)=h2(i)/h2len
      enddo

c     minimization is carried out in the direction of vector h2.
      if (ifun.eq.0) then
         call linmind(mxatom,natom,nmove,itype,taxes,ietype,maxtyp,rcut
     $        ,hyperratio,p1,h2,work1,drmax,nstp2,ein,eout,dE,cgfac,pmin
     $        ,pmax,tol)
      else if(ifun.eq.1) then
         call linmine(mxatom,natom,nmove,itype,taxes,ietype,maxtyp,rcut
     $        ,hyperratio,p1,h2,drmax,nstp2,-h2len,work1,ein,eout,dE
     $        ,cgfac,pmin,pmax,tol)
      endif

      do i=1,natom
         xyz((i-1)*3+1)=p1(i,1)
         xyz((i-1)*3+2)=p1(i,2)
         xyz((i-1)*3+3)=p1(i,3)
      enddo
      
      do i=1,3*nmove
         g1(i)=g2(i)
      enddo
      
      if (etmp.ne.0.0d0) ein=etmp
      nstep=nstp1+nstp2
      
      return
      end

      subroutine linmind(nmax,natom,nmove,itype,taxes,ietype,maxtyp,rcut
     $     ,hyperratio,p0,h1,work1,drmax,nstep,einit,emin,dE,cgfac,pmin
     $     ,pmax,tol)
c     Implemented from within opto10...
c     Excess arguments for gcalc call...
      implicit real*8(a-h,o-z)
      parameter (gold=2.618)
      parameter (iter=500)

      dimension dr(5)
      dimension work1(4*nmax),h1(*),p0(nmax,3)

      dimension itype(natom),taxes(3),rcut(maxtyp,maxtyp)

      save icnt,dr

      if (pmin.le.0) pmin=0.1d0
      if (pmax.le.0) pmax=2.0d0
      pinit=0
      imin=0
      i=0

      if (tol.le.0.0d0) tol=1d-4
c      if (cgfac.le.0) cgfac=gfac
      if (drmax.le.0) then
         icnt1=0
         xtol=1.0d-2
      else
         dr(5)=dr(4)
         dr(4)=dr(3)
         dr(3)=dr(2)
         dr(2)=dr(1)
         dr(1)=drmax
         
         drsum=0.0d0
         icnt=0
         do 8 i=1,5
            if (dr(i).gt.0) then
               drsum=drsum+dr(i)
               icnt=icnt+1
            endif
 8       continue

         xtol=(drsum/icnt)*tol
      endif
      if (xtol.lt.1.0d-16) xtol=1.0d-16
      
      call df(nmax,natom,nmove,itype,taxes,ietype,maxtyp,rcut,hyperratio
     $     ,p0,pinit,h1,dfinit,edum,work1)
      nstep=1
      if(dfinit.eq.0.0d0) then
         call xmin=xtol
         call rasmbl(nmax,natom,nmove,p0,h1,xmin,work1)
         call energysub(emin,natom,natom,work1(1),work1(natom+1),
     x        work1(2*natom+1),itype,taxes,ietype,maxtyp,rcut,hyperratio
     $        )
         goto 210
      endif
      
      p1=pinit
      df1=dfinit
      e1=einit
      p2=p1 - df1*cgfac
      
    9 continue
      d21=abs(p2-p1)
      if (d21.gt.pmax) p2=sign(pmax,p2)
      if (d21.lt.pmin) then
         p2=sign(pmin,p2)
         imin=1
      endif
      call df(nmax,natom,nmove,itype,taxes,ietype,maxtyp,rcut,hyperratio
     $     ,p0,p2,h1,df2,e2,work1)
      nstep=nstep+1
      
      if (df1*df2.lt.0) then
         if (imin.eq.0) then
            cgfac=cgfac*0.5d0
            p2=(p2-p1)*0.1d0
            goto 9
         else
            goto 100
         endif
      else if (df2.eq.0) then
         xmin=p2
         emin=e2
         goto 210
      endif
      
      do 20 i=1,iter

c     do parabolic fit to get next point...
c     xcurv=(df1-df2)/(2*(p1-p2))
         pnum=p1*df2-p2*df1
         pden=df2-df1
         convex=(df2-df1)/(p2-p1)

         if(convex.lt.0.0d0) then
c     if geometry is convex then step towards minimum...
            p3=p2 - df2*cgfac
         else if (pden.ne.0.0d0 .and. pnum.ne.0.0d0) then
            p3=(pnum/pden)
         else
            p3=p1 + gold*(p2-p1)
         endif
         
         d32=abs(p3-p2)
         if (d32.ge.pmax) p3=p2 + sign(pmax,p3)
         if (d32.lt.pmin) p3=p2 + sign(pmin,p3)
         
         call df(nmax,natom,nmove,itype,taxes,ietype,maxtyp,rcut
     $        ,hyperratio,p0,p3,h1,df3,e3,work1)
         nstep=nstep+1
         if (df3.eq.0) then
            xmin=p3
            emin=e3
            goto 210
         endif
         
         p1=p2
         df1=df2
         e1=e2
         p2=p3
         df2=df3
         e2=e3
         if (df1*df2.lt.0) then
            goto 100
         endif
         
 20   continue

      stop 'max number of iters exceeded before bracketing'

 100  continue
c     check to make sure routine didn't step uphill...
      if (abs(p2).lt.abs(p1)) then
         xtmp=p1
         p1=p2
         p2=xtmp
         xtmp=e1
         e1=e2
         e2=xtmp
      endif

      eb1=e1
      eb2=e2
      b1=p1
      b2=p2
      
      if (dfinit*p2.gt.0) stop 'linmin: routine stepped uphill...'
c     zoom in on minimum, using the same derivative method....

      do 200 i=1,iter
         
         d21=abs(p2-p1)
         if (d21.le.tol) then
            if (abs(df1).le.abs(df2)) then
               xmin=p1
               emin=e1
            else
               xmin=p2
               emin=e2
            endif
            goto 210
         endif

         pnum=p1*df2-p2*df1
         pden=df2-df1
         
         if (pden.ne.0) then
            p3=pnum/pden
            d31=abs(p3-p1)
            d32=abs(p3-p2)
            if ((d31.gt.d21 .or. d32.gt.d21) .and. df2*df1.lt.0)
     x           p3=p2-(p2-p1)*0.5d0
         else
            p3=p2 - (p2-p1)*0.5d0
         endif

         call df(nmax,natom,nmove,itype,taxes,ietype,maxtyp,rcut
     $        ,hyperratio,p0,p3,h1,df3,e3,work1)
         nstep=nstep+1

         if (df3.eq.0) then
            xmin=p3
            emin=e3
            goto 210
         else if (abs(df1).lt.1d-1 .or. abs(df2).lt.1d-1) then
c Take two points with lowest derivative, namely p3 and...
            if (abs(df1).le.abs(df2)) then
               p2=p1
               df2=df1
               e2=e1
               p1=p3
               df1=df3
               e1=e3
            else
               p1=p3
               df1=df3
               e1=e3
            endif
         else
c If derivatives still large, then use p2 and p3 as newest points...
            if (df3.ne.0) then
               p1=p2
               df1=df2
               e1=e2
               p2=p3
               df2=df3
               e2=e3
            else
c     This means that df3=0, since df1.ne.0
               xmin=p3
               emin=e3
               goto 210
            endif
         endif

 200  continue

      stop 'linmin: minimum not found'

 210  continue
      if (xmin.eq.0) then
         if (dfinit.gt.0) then
            xmin=-xtol
            call rasmbl(nmax,natom,nmove,p0,h1,xmin,work1)
         else
            xmin=xtol
            call rasmbl(nmax,natom,nmove,p0,h1,xmin,work1)
         endif
         call energysub(emin,natom,natom,work1(1),work1(natom+1),
     x        work1(2*natom+1),itype,taxes,ietype,maxtyp,rcut,hyperratio
     $        )
      endif
      dE=emin-einit
      
      if (abs(xmin + dfinit*cgfac).lt.5.0d0) then
         icnt1=icnt1+1
         cgfac=icnt1*(cgfac - xmin/dfinit)/(icnt1+1)
      endif
      
      call rasmbl(nmax,natom,nmove,p0,h1,xmin,work1)
      do 220 i=1,natom
         p0(i,1)=work1(i)
         p0(i,2)=work1(natom+i)
         p0(i,3)=work1(2*natom+i)
 220  continue
      
 500  continue
      return
      end

      subroutine linmine(nmax,natom,nmove,itype,taxes,ietype,maxtyp,rcut
     $     ,hyperratio,p0,h1,drmax,nstep,dfinit,r,einit,emin,dE,cgfac
     $     ,pmin,pmax,tol)
c  Implemented from within opto10...
c  Excess arguments for gcalc call...
c  routine returns the number of energy calls necessary, nstep...
      implicit real*8(a-h,o-z)
      parameter (iter=100)
      
      dimension dr(5)
      dimension h1(*),p0(nmax,3),r(4*nmax)
      
      dimension itype(natom),taxes(3),rcut(maxtyp,maxtyp)

      save icnt,dr
      
      if (pmin.le.0.0d0) pmin=0.1d0
      if (pmax.le.0.0d0) pmax=2.0d0
      pinit=0.0d0
      nstep=0
      icnv=0
      i=0
      
      if (tol.le.0.0d0) tol=1d-4
c      if (cgfac.le.0.0d0) cgfac=gfac
      if (drmax.le.0.0d0) then
         icnt1=0
         xtol=1.0d-2
      else
         dr(5)=dr(4)
         dr(4)=dr(3)
         dr(3)=dr(2)
         dr(2)=dr(1)
         dr(1)=drmax
         
         drsum=0.0d0
         icnt=0
         do 7 i=1,5
            if (dr(i).gt.0.0d0) then
               drsum=drsum+dr(i)
               icnt=icnt+1
            endif
 7       continue
         
         xtol=(drsum/icnt)*tol
      endif
      if (xtol.lt.1.0d-7) xtol=1.0d-8
      
      if (dfinit.eq.0.0d0) then
         xmin=xtol
         emin=einit
         goto 210
      endif
      
      p1=pinit
      e1=einit
      
      p2=p1 - cgfac*dfinit
      d21=abs(p2-p1)
      if (d21.gt.pmax) p2=p1 + sign(pmax,p2)
      if (d21.lt.pmin) p2=p1 + sign(pmin,p2)
      
      call rasmbl(nmax,natom,nmove,p0,h1,p2,r)
      call energysub(e2,natom,natom,r(1),r(natom+1),r(2*natom+1),itype
     $     ,taxes,ietype,maxtyp,rcut,hyperratio)
      nstep=nstep+1
      
 12   continue
      if ((e2-e1).gt.0.0d0) then
         p3=p2
         e3=e2
         p2=(p2-p1)*0.5d0
         call rasmbl(nmax,natom,nmove,p0,h1,p2,r)
         call energysub(e2,natom,natom,r(1),r(natom+1),r(2*natom+1)
     $        ,itype,taxes,ietype,maxtyp,rcut,hyperratio)
         nstep=nstep+1
         
         d21=abs(p2-p1)
         if (d21.gt.pmin) goto 12
         goto 102
      endif
      
      pnum=2.0d0*(p1*(e2-e1)) - (p2**2-p1**2)*dfinit
      pden=-2.0d0*dfinit*(p2-p1) + 2.0d0*(e2-e1)
      if (pnum.ne.0.0d0 .or. pden.ne.0.0d0) then
         p3=pnum/pden
      else
         p3=p2 + sign(pmin,p1)
      endif
      if (p3*p2.lt.0.0d0) p3=p2 + sign(pmin,p2)
c     since search starts at zero, p3 and p2 must have same sign, else
c     in convex region...
      
 15   continue
      d32=abs(p3-p2)
      if (d32.gt.pmax) p3=p2 + sign(pmax,p3)
      if (d32.lt.pmin) then
         p3=p2 + sign(pmin,p3)
      endif
      
      call rasmbl(nmax,natom,nmove,p0,h1,p3,r)
      call energysub(e3,natom,natom,r(1),r(natom+1),r(2*natom+1),itype
     $     ,taxes,ietype,maxtyp,rcut,hyperratio)
      nstep=nstep+1
      
      d12=abs(p1-p2)
      d13=abs(p1-p3)
      d23=abs(p2-p3)
c     find middle state and check to see if it has lowest energy...
      xmax1=max(d12,d13)
      xmax2=max(d12,d23)
      xmax3=max(d13,d23)
      imin=1
      xlmin=xmax1
      if (xmax2.lt.xlmin) then
         imin=2
         xlmin=xmax2
      endif
      if (xmax3.lt.xlmin) then
         imin=3
         xlmin=xmax3
      endif
      
      if (d32.le.pmin) then
c     check to make sure imin is bound...
         if (imin.eq.1 .and. ( e3.gt.e1 .and. e2.gt.e1 )) goto 100
         if (imin.eq.2 .and. ( e3.gt.e2 .and. e1.gt.e2 )) goto 100
         if (imin.eq.3 .and. ( e2.gt.e3 .and. e1.gt.e3 )) goto 100
      else
c     check to make sure imin is bound...
         ibck=0
         if (imin.eq.1 .and. ( e3.gt.e1 .and. e2.gt.e1 )) ibck=1
         if (imin.eq.2 .and. ( e3.gt.e2 .and. e1.gt.e2 )) ibck=1
         if (imin.eq.3 .and. ( e2.gt.e3 .and. e1.gt.e3 )) then
            ibck=1
         else
            icnv=1
         endif
         if (ibck.eq.1) then
            p3=p2 + (p3-p2)*0.5d0
            goto 15
         endif
         
      endif
      
      do 20 i=1,iter

c     do parabolic fit to get next point...
c     check to see if colinear...
         pnum=(e2-e3)*(p1**2-p2**2) - (e1-e2)*(p2**2-p3**2)
         pden=2*( (e2-e3)*(p1-p2) - (e1-e2)*(p2-p3) )

         if (pden.ne.0.0d0 .and. pnum.ne.0.0d0) then
            p4=(pnum/pden)
         else
            p4=p3 + sign(pmin,p3)
         endif

         if ( (p4-p3)*(p3-p2).le.0.0d0 .or. icnv.eq.1) then
c     if geometry is convex then step towards minimum...
            p4=p3 + sign(pmin,p3)
         endif

c     make sure step is within minimum/maximum bounds...
         d43=abs(p4-p3)
         if (d43.ge.pmax) then
            p4=p3 + sign(pmax,p4)
         else if (d43.lt.pmin) then
            p4=p3 + sign(pmin,p4)
         endif
         
         call rasmbl(nmax,natom,nmove,p0,h1,p4,r)
         call energysub(e4,natom,natom,r(1),r(natom+1),r(2*natom+1)
     $        ,itype,taxes,ietype,maxtyp,rcut,hyperratio)
         nstep=nstep+1
         
         p1=p2
         e1=e2
         p2=p3
         e2=e3
         p3=p4
         e3=e4
         
         d12=abs(p1-p2)
         d13=abs(p1-p3)
         d23=abs(p2-p3)
c     find middle state and check to see if it has lowest energy...
         xmax1=max(d12,d13)
         xmax2=max(d12,d23)
         xmax3=max(d13,d23)
         imin=1
         xlmin=xmax1
         if (xmax2.lt.xlmin) then
            imin=2
            xlmin=xmax2
         endif
         if (xmax3.lt.xlmin) then
            imin=3
            xlmin=xmax3
         endif

c     check to make sure imin is bound...
         if (imin.eq.1 .and. ( e3.gt.e1 .and. e2.gt.e1 )) goto 100
         if (imin.eq.2 .and. ( e3.gt.e2 .and. e1.gt.e2 )) goto 100
         if (imin.eq.3 .and. ( e2.gt.e3 .and. e1.gt.e3 )) goto 100
         
 20   continue
      
      write(6,30)
 30   format(2x,'linmine: error --- loop finished before bracketing')
      stop 'linmine: error'
      
 100  continue
c     check to make sure routine didn't step uphill...
      if (dfinit*p3.gt.0) then
         write(6,101)
 101     format(2x,'linmine: routine stepped uphill...')
         stop'linmine: error'
      endif
 102  continue
      if (p1.gt.p3) then
         xtmp=e1
         e1=e3
         e3=xtmp
         xtmp=p1
         p1=p3
         p3=xtmp
      endif
      
      eb1=e1
      eb2=e3
      b1=p1
      b2=p3
      
      do 200 i=1,iter
         
         d31=abs(p3-p1)
         d32=abs(p3-p2)
         d21=abs(p2-p1)
         dmax=max(d31,d32)
         dmax=max(dmax,d21)
         if (dmax.le.xtol) then
            if (e2.le.eb1 .and. e2.le.eb2) then
               xmin=p2
               emin=e2
            else
               xmin=p2
               emin=e2
            endif
            goto 210
         endif

         pnum=(e2-e3)*(p1**2-p2**2) - (e1-e2)*(p2**2-p3**2)
         pden=2*( (e2-e3)*(p1-p2) - (e1-e2)*(p2-p3) )

         if (pden.ne.0.0d0 .and. pnum.ne.0.0d0) then
            p4=(pnum/pden)
            d41=abs(p4-p1)
            if (p4.lt.b1 .or. p4.gt.b2) then
               p4=p3 - (p3-p1)*0.5
            endif
         else
            p4=p3 - (p3-p1)*0.5
         endif

         call rasmbl(nmax,natom,nmove,p0,h1,p4,r)
         call energysub(e4,natom,natom,r(1),r(natom+1),r(2*natom+1)
     $        ,itype,taxes,ietype,maxtyp,rcut,hyperratio)
         nstep=nstep+1

         if((e1-e2).eq.0.0d0 .or. (e1-e3).eq.0.0d0 .or.
     x        (e1-e4).eq.0.0d0 .or. (e2-e3).eq.0.0d0 .or.
     x        (e2-e4).eq.0.0d0 .or. (e3-e4).eq.0.0d0) then
c     limited by roundoff-error... find lowest energy point and return
            emin=e1
            xmin=p1
            if (e2.lt.emin) then
               emin=e2
               xmin=p2
            endif
            if (e3.lt.emin) then
               emin=e3
               xmin=p3
            endif
            if (e4.lt.emin) then
               emin=e4
               xmin=p4
            endif
            goto 210
         endif

c     take three points with lowest energy...
         imax=1
         emax=e1
         if (e2.gt.emax) then
            imax=2
            emax=e2
         endif
         if (e3.gt.emax) then
            imax=3
            emax=e3
         endif
         if (e4.gt.emax) then
            imax=4
            emax=e4
         endif
         
         if (imax.eq.1) then
            p1=p2
            p2=p3
            p3=p4
            e1=e2
            e2=e3
            e3=e4
         else if (imax.eq.2) then
            p1=p1
            p2=p3
            p3=p4
            e1=e1
            e2=e3
            e3=e4
         else if (imax.eq.3) then
            p1=p1
            p2=p2
            p3=p4
            e1=e1
            e2=e2
            e3=e4
         else if (imax.eq.4) then
c     if new point is highest in energy, take it anyhow...
            p4=p3 - (p3-p1)*0.5
            call rasmbl(nmax,natom,nmove,p0,h1,p4,r)
            call energysub(e4,natom,natom,r(1),r(natom+1),r(2*natom+1)
     $           ,itype,taxes,ietype,maxtyp,rcut,hyperratio)
            p1=p1
            p2=p2
            p3=p4
            e1=e1
            e2=e2
            e3=e4
         endif
         
 200  continue
      
      stop 'linmin:  minimum not found'

 210  continue
      if (xmin.eq.0.0d0) then
         if (dfinit.gt.0.0d0) then
            xmin=-xtol
         else
            xmin=xtol
         endif
         call rasmbl(nmax,natom,nmove,p0,h1,xmin,r)
         call energysub(emin,natom,natom,r(1),r(natom+1),r(2*natom+1)
     $        ,itype,taxes,ietype,maxtyp,rcut,hyperratio)
      endif
      dE=emin-einit
      
      if (abs(xmin + dfinit*cgfac).lt.5.0d0) then
         icnt1=icnt1+1
         cgfac=icnt1*(cgfac - xmin/dfinit)/(icnt1+1)
      endif
      
      call rasmbl(nmax,natom,nmove,p0,h1,xmin,r)
      do 220 i=1,natom
         p0(i,1)=r(i)
         p0(i,2)=r(natom+i)
         p0(i,3)=r(2*natom+i)
 220  continue
      
 500  continue
      return
      end
      
      subroutine rasmbl(nmax,natom,nmove,p0,h1,p,pnew)
c     this routine does the following.... p0 + p*h1,
c     where p is a scalar, p0(1...3*natom), h1(1...3*nmove) are vectors
      implicit real*8(a-h,o-z)
      dimension p0(nmax,3),h1(*),pnew(4*nmax)

c      write(6,*) 'rasmbl: ',nmax,natom,p

      j=-3
      do i=1,nmove
         j=j+3
         pnew(i) = p0(i,1) + p*h1(j+1)
         pnew(natom+i) = p0(i,2) + p*h1(j+2)
         pnew(2*natom+i) = p0(i,3) + p*h1(j+3)
      enddo
      do i=nmove+1,natom
         pnew(i) = p0(i,1)
         pnew(natom+i) = p0(i,2)
         pnew(2*natom+i) = p0(i,3)
      enddo

c      write(6,*) 'rasmbl: leaving'

      return
      end

      subroutine df(nmax,natom,nmove,itype,taxes,ietype,maxtyp,rcut
     $     ,hyperratio,p0,p,dvec,der,e,workarray)
c  returns 1-dimension derivative at point p0 + p*dvec along dvec...
      implicit real*8(a-h,o-z)
c      parameter (mxatom=10000)
      include 'parameters.h'
      dimension p0(nmax,3),dvec(*),workarray(4*nmax)
      dimension itype(natom),taxes(3),rcut(maxtyp,maxtyp)
      
      igrad=3*natom

      j=-3
      do i=1,nmove
         j=j+3
         workarray((i-1)*3+1)=p0(i,1) + p*dvec(j+1)
         workarray((i-1)*3+2)=p0(i,2) + p*dvec(j+2)
         workarray((i-1)*3+3)=p0(i,3) + p*dvec(j+3)
      enddo
      do i=nmove+1,natom
         workarray((i-1)*3+1)=p0(i,1)
         workarray((i-1)*3+2)=p0(i,2)
         workarray((i-1)*3+3)=p0(i,3)
      enddo

      call gcalc(natom,nmove,workarray,itype, taxes,ietype,maxtyp,rcut,
     x     e,workarray(igrad+1),hyperratio)
      
      der=0
      do i=1,3*nmove
         der=der + workarray(igrad+i)*dvec(i)
      enddo
      
 100  continue
      return
      end

      subroutine energysub(energy,natom,nmove,x,y,z,itype,taxes,ietype
     $     ,maxtyp,rcut,hyperratio)
c wrapper for all of the energy calls in above stuff
      implicit real*8(a-h,o-z)
      dimension x(natom),y(natom),z(natom),itype(natom)
      dimension taxes(3),rcut(maxtyp,maxtyp)

c      parameter (mxatom=10000)
      include 'parameters.h'
      dimension workarray(4*mxatom)

      do i=1,natom
         workarray((i-1)*3+1)=x(i)
         workarray((i-1)*3+2)=y(i)
         workarray((i-1)*3+3)=z(i)
      enddo

      call gcalc(natom,nmove,workarray,itype, taxes,ietype,maxtyp,rcut,
     x     energy,workarray(igrad+1),hyperratio)
      
      return
      end
