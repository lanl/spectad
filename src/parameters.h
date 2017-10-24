      integer natmax,mxatom,nmaxrotalign,nholdmdsubs,movmax,maxmove
     +    ,mximage,nneighmax,nstatemax
            
     
c parameter natmax, mxatom, nmax, nhold used by
c acg.f
c conjgrad.f
c descentcheck.f
c dimer_tad.f
c e0subs_tad3.f
c e999subs.f
c mdsubs_tad3.f
c neb_tad.f
c refinetransition.f
c rotalign.f
c tad3.f
c tersubs_tad2.real.f

      parameter (natmax=2000)
      parameter (mxatom=natmax)
      parameter (nmaxrotalign=natmax) ! name changes from nmax since nmax is used elsewhere
      parameter (nholdmdsubs=natmax) ! name changed from nhold since that is used elsewhere for other things

c parameter movmax, maxmove used by
c mdsubs_tad3.f
c neb_tad.f
c nmodes.f
c nmodes_subset.f

      parameter(movmax=natmax)
      parameter(maxmove=movmax)

c parameter mximage used by
c neb_tad.f

      parameter(mximage=42)

c parameter nneighmax used by
c descentcheck.f
c dimer_tad.f
c tad3.f

      parameter(nneighmax=200)

c parameter nstatemax used by
c descentcheck.f (but only needed for itranscrit.eq.3)

      parameter(nstatemax=1500)

c other parameters to consider
c--
ce0subs_tad3.f:      parameter(mxijk=10000,mxpont=100000)
c--
ce0subs_tad3.f:      parameter(mxijk=10000,mxpont=100000)
c--
ce999subs.f:      parameter (lwrk1=15000000,nnmax=600,llmax=500000)
c--
ctersubs_tad2.real.f:      parameter (nnmax=13,lmax=40000)
c--
ctersubs_tad2.real.f:      parameter (nnmax=13,lmax=40000)
c--
ctersubs_tad2.real.f:      parameter (nnmax=13,lmax=40000)

