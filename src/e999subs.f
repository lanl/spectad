c-----------------------------------------------------------------------
      subroutine g999
     x    (natom,nmove,x,y,z,itype,maxtyp,tax,tay,taz,rcut,energy,grad)
c-----------------------------------------------------------------------

      implicit real*8(a-h,o-z)
      real*8 x, y, z, xj, yj, zj, rj
      integer natom, itype, i, it, j, jt, marki, jj
c      integer lwrk1,llmax,nnmax,iwrk1,jwrk1,nneigh,markn,jatomn
c      integer maxtyp,nmove,lmax,mmax,nmax,ijdx
      real*8 tax,tay,taz,tx,ty,tz,ax,ay,az,energy,rcut,rct
c      real*8 dxn,dyn,dzn,rn,fill1
      real*8 grad, gx, gy, gz
      real*8 rk, pad

      dimension x(natom), y(natom), z(natom), itype(natom)
      dimension grad(3,natom)
      dimension rcut(maxtyp,maxtyp)
      common/potcom/rct(5,5),tx,ty,tz,ax,ay,az
      include 'wrksubs_tad2.h'
      include 'parameters.h'
      parameter(llmax1=2*(lwrk1-7*natmax/2)/9-1 )
      parameter(llmaxh=llmax1/2)
      parameter(llmax=llmaxh*2)   ! this makes it an even number
      parameter(nnmax=600) ! xxx verify what this is for

c  Voter work common
      common/wrkcm1/iwrk1,jwrk1,nneigh(natmax),markn(natmax),
     x         jatomn(llmax),dxn(llmax),dyn(llmax),dzn(llmax),rn(llmax),
     x         fill1(lwrk1-natmax-9*llmax/2)

      rk=0.1d0

c kludge for now - move these into call
      tx=tax
      ty=tay
      tz=taz
      ax=tax/2.0d0
      ay=tay/2.0d0
      az=taz/2.0d0
      do i=1,maxtyp
         do j=1,maxtyp
            rct(i,j)=rcut(i,j)
         end do
      end do

      if(iwrk1.gt.0) then
         stop 'g999: wrkcm1 conflict'
         return
      else
         iwrk1=1
      end if

      pad = 0.0d0

c rlistr does the neighbor stuff

      call rlistr (natom, 1, natom, x, y, z, itype, pad,
     x             nneigh, markn, jatomn, dxn, dyn, dzn, rn, llmax)

      do i=1,natom
         grad(1,i) = 0.0d0
         grad(2,i) = 0.0d0
         grad(3,i) = 0.0d0
      enddo
      do i=1,nmove
         gx = 0.0d0
         gy = 0.0d0
         gz = 0.0d0
         it=itype(i)
         marki = markn(i)
         
         do jj=1,nneigh(i)
            if (nneigh(i) .gt. nnmax) then
               write(6,'("nneigh(",i3,")=",i4," nnmax=",i4)') 
     x             i,nneigh(i),nnmax
               stop 'g999: increase nnmax'
            endif
            ijdx = marki + jj
            j = jatomn(ijdx)
            xj = dxn(ijdx)
            yj = dyn(ijdx)
            zj = dzn(ijdx)
            rj = rn(ijdx)
            jt = itype(j)

            if (xj.lt.rcut(1,1)) gx=gx+rk*xj/rj
            if (yj.lt.rcut(1,1)) gy=gy+rk*yj/rj
            if (zj.lt.rcut(1,1)) gz=gz+rk*zj/rj

c            if (rj.gt.rcut(1,1)) then
c               write(6,*) 'g999: i,j,rj: ',rj,xj,yj,zj
c               write(6,*) 'g999: i,j,gx: ',gx,gy,gz
c            endif

         enddo
         
c    Nearest Neighbour loop ended here. 
c    Total the forces on atom so far. 
c    These include short range and real space coloumbic
         
         grad(1,i)=gx
         grad(2,i)=gy
         grad(3,i)=gz

c         if (i.eq.1) write(6,*) 'g999: i: ',i,grad(1,i),grad(2,i),grad(3
c     +       ,i)
         
      enddo
      
      iwrk1=0
      
      return
      end

