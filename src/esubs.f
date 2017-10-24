      function energy(natom,nmove,x,y,z,itype)
c  this routine returns total energy over the first nmove atoms
      implicit real*8(a-h,o-z)
      dimension x(natom),y(natom),z(natom),itype(natom)
      parameter (nnmax=200)
      dimension jat1(nnmax),rn1(nnmax),jtyp1(nnmax)
c
      maxn=nnmax
      e=0.0d0
      do 20 i=1,nmove
      call rlist1(natom,x,y,z,itype,
     x                          i,  nn,jat1,jtyp1,rn1,maxn)
      it=itype(i)
      call rhatm3(it,nn,rn1,jtyp1, rho,ei)   
c                              ! 3-pt interp for higher accuracy
      e=e+ei
   20 continue
      energy=e
      return
      end

      subroutine gatom(it,rhoi,nneigh,dx,dy,dz,r,itype,rhotot, 
     x                                           gx,gy,gz)     
c  this routine computes the gradient (gx,gy,gz) on *one* atom.
c  calling routine supplies the type of this atom (ityp), the
c  electron density at this atom (rhoi),the number of
c  neighbors (nneigh), the relative displacements to each neighbor
c  (dx(nneigh),dy(nneigh),dz(nneigh),r(nneigh)), the type of each
c  neighbor (itype(nneigh)), and the electron density at each
c  neighbor atom (rhotot(nneigh)).
c
      implicit real*8(a-h,o-z)
c     dimension dx(nneigh),dy(nneigh),dz(nneigh),r(nneigh)
c     dimension itype(nneigh),rhotot(nneigh)
      dimension dx(*),dy(*),dz(*),r(*)
      dimension itype(*),rhotot(*)
c
c  statement functions for interpolation fits
c  note:  these assume that sgcon has been performed to eliminate s() and g()
      rho(rr,it)=trpfun(rr,0,70+it)
      tgrho(rr,it)=trpfun(rr,11,70+it)
      phi(rr,it,jt)=trpfun(rr,0,80+it*jt)
      tgphi(rr,it,jt)=trpfun(rr,1,80+it*jt)/rr
      f(x,it)=trpfun(x,0,90+it)
      fp(x,it)=trpfun(x,1,90+it)
c
c  compute gradient on this atom
c
      gx=0.0d0
      gy=0.0d0
      gz=0.0d0
      fpi=fp(rhoi,it)
      do 20 j=1,nneigh
      jt=itype(j)
      tgrad=tgphi(r(j),it,jt) + fpi*tgrho(r(j),jt)
     x                           + fp(rhotot(j),jt)*tgrho(r(j),it)
      gx=gx+tgrad*dx(j)
      gy=gy+tgrad*dy(j)
      gz=gz+tgrad*dz(j)
   20 continue
c
      return
      end

      subroutine rhatom(ityp,nneigh,r,itype,    
c                              ! supplied
     x                           rhotot,energy)   
c                              ! returned
c  this routine computes the electron density at *one* atom.
c  the calling routine supplies the type of this atom (ityp), the number
c  of neighbors (nneigh), the distance to each neighbor (r(nneigh)),
c  and the type of each neighbor (itype(nneigh)).
c  this routine returns the energy of this atom (energy),
c  and the electron density (rhotot).
c
      implicit real*8(a-h,o-z)
c     dimension r(nneigh),itype(nneigh)
      dimension r(*),itype(*)
      common/abparm/s(2),sinv(2),g(2)
c
c  statement functions for interpolation fits
c  note:  these assume that sgcon has been performed to eliminate s() and g()
      rho(rr,it)=trpfun(rr,0,70+it)
      phi(rr,it,jt)=trpfun(rr,0,80+it*jt)
      f(x,it)=trpfun(x,0,90+it)
c
c   compute energy and rhotot
c
      epair=0.0d0
      it=ityp
      rhotot=0.0d0
      do 40 j=1,nneigh
      jt=itype(j)
      epair=epair + phi(r(j),it,jt)
      rhotot=rhotot+rho(r(j),jt)
   40 continue
c
      energy = f(rhotot,it) + 0.5d0*epair
c
      return
      end

      subroutine rhatm3(ityp,nneigh,r,itype,    
c                              ! supplied
     x                           rhotot,energy)   
c                              ! returned
c like rhatom, but uses 3-point interpolation for higher accuracy
c  this routine computes the electron density at *one* atom.
c  the calling routine supplies the type of this atom (ityp), the number
c  of neighbors (nneigh), the distance to each neighbor (r(nneigh)),
c  and the type of each neighbor (itype(nneigh)).
c  this routine returns the energy of this atom (energy),
c  and the electron density (rhotot).
c
      implicit real*8(a-h,o-z)
c     dimension r(nneigh),itype(nneigh)
      dimension r(*),itype(*)
      common/abparm/s(2),sinv(2),g(2)
c  test stuff:
c      save sig,eps
c      data sig/2.282d0/,eps/1.9097d-2/  
c                              ! h&p lj parms for ni
c      phi(r,it,jt)=4.0d0*eps*((sig/r)**12-(sig/r)**6)
c
c ***  3-point interpolation for higher accuracy ***
      rho(rr,it)=trpfn3(rr,0,70+it)
      phi(rr,it,jt)=trpfn3(rr,0,80+it*jt)
      f(x,it)=trpfn3(x,0,90+it)
c
c   compute energy and rhotot
c
      epair=0.0d0
      it=ityp
      rhotot=0.0d0
      do 40 j=1,nneigh
      jt=itype(j)
      epair=epair + phi(r(j),it,jt)
      rhotot=rhotot+rho(r(j),jt)
   40 continue
c
      energy = f(rhotot,it) + 0.5d0*epair
c
      return
      end

      subroutine rlist1(natom,x,y,z,itype,iatom,      
c                              ! supplied
     x                  nn,jat1,jtyp1,rn1,maxn)
c  general neighbor-listing routine for *one* atom (iatom)
c  user supplies:  natom, x(),y(),z(),itype(),
c                  maxn=maximum length for jat1,jtyp1,rn1 arrays
c                  iatom = atom for which neighbor list is to be constructed
c  rlist1 returns:
c                  nn= no. of neighbors for this atom
c                  jat1(),jtyp1(),rn1() - length<=maxn
c              jat1(j) = j'th neighbor of atom iatom
c              jtyp1(j) = type of j'th neighbor
c              rn(j) = distance from atom iatom to j'th neighbor
c
      implicit real*8(a-h,o-z)
      dimension x(natom),y(natom),z(natom),itype(natom)
      dimension jat1(maxn),jtyp1(maxn),rn1(maxn)
      common/potcom/rcut(3,3),tax,tay,taz,ax,ay,az
c
      k=0
      i=iatom
      it=itype(i)
      do 10 j=1,natom
      if(j.eq.i) go to 10
      dx=x(i)-x(j)
      if(dx.gt.ax) dx=dx-tax
      if(dx.lt.-ax) dx=dx+tax
      dy=y(i)-y(j)
      if(dy.gt.ay) dy=dy-tay
      if(dy.lt.-ay) dy=dy+tay
      dz=z(i)-z(j)
      if(dz.gt.az) dz=dz-taz
      if(dz.lt.-az) dz=dz+taz
      r=sqrt(dx*dx+dy*dy+dz*dz)
      if(r.gt.rcut(it,itype(j))) go to 10
        k=k+1
        jat1(k)=j
        jtyp1(k)=itype(j)
        rn1(k)=r
   10 continue
      if(k.gt.maxn) stop 'abort - out of space in routine rlist1'
      nn=k
      return
      end

      subroutine potset(ntype,tax,tay,taz,iwrt)
c  initializes /potcom/ for use in cgsubs routines
c  user supplies ntype,tax,tay,taz
c  this routine fills rcut(i,j), checks that ntype is ok,
c  reads in trp arrays, etc.
      implicit real*8(a-h,o-z)
      dimension title(10)
      common/potcom/rcut(3,3),tx,ty,tz,ax,ay,az
      common/azeros/azero(3)
      common/masses/amu(3)
      common/trpcom/mute
c
      if(ntype.gt.3) stop 'abort - ntype.gt.3 in potset'
      if(iwrt.ge.1) write(6,1) tax,tay,taz
    1 format(' tax =',f18.8,'  tay =',f18.8,'  taz =',f18.8)
c
c  read in title cards and atomic masses from unit 60+itype
c
      if(iwrt.ge.0) write(6,2) ntype
    2 format(' here are the title cards from the',i2,' potentials:')
      do 5 i=1,ntype
      iu=60+i
      rewind iu
      read(iu,3) title
      read(iu,*) amu(i)
      read(iu,*) azero(i)
    3 format(10a8)
      if(iwrt.ge.0) write(6,4) title
    4 format(1x,10a8)
    5 continue

      write(*,*) 'IN ESUBS -- POTSET !!!' 


      if(iwrt.ge.1) mute=0
      if(iwrt.le.0) mute=1
      call taxset(tax,tay,taz)
c
c procure cutoff values from trp arrays
c note that the following calls are very cheap
c note that these calls automatically read in the trp arrays if necessary
      do 10 i=1,ntype
      call trpdat(70+i,npt,x0,x1,rhoci,y1,yn,iflo,ifhi)
      call trpdat(90+i,npt,x0,x1,xn,y1,yn,iflo,ifhi)
      do 10 j=1,ntype
      call trpdat(70+j,npt,x0,x1,rhocj,y1,yn,iflo,ifhi)
      call trpdat(80+i*j,npt,x0,x1,phicij,y1,yn,iflo,ifhi)
      rcut(i,j)=max(rhoci,rhocj,phicij)
   10 continue
c
      if(iwrt.eq.0) then
        write(6,20) ntype
   20   format(' trp arrays for ntype =',i2,
     x         ' have been read from units 70+i,80+i*j,90+i')
      end if
c
      nbad=0
      do 30 i=1,ntype
      do 30 j=1,i
      rc=rcut(i,j)
      if(rc.gt.ax.or.rc.gt.ay.or.rc.gt.az) then
        nbad=nbad+1
      end if
   30 continue
      if(nbad.gt.0) write(6,40) nbad
   40 format('  warning -',i3,' cutoffs are .gt. 1/2(tax,y, or z)')
c
      return
      end

      subroutine taxset(tax,tay,taz)
      implicit real*8(a-h,o-z)
      common/potcom/rcut(3,3),tx,ty,tz,ax,ay,az
      tx=tax
      ty=tay
      tz=taz
      ax=tax/2.0d0
      ay=tay/2.0d0
      az=taz/2.0d0
      return
      end
