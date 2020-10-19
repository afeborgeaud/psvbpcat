      program psvbpcat
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  ************** psvbp.f ****************
c Computation of back-propagated wavefield 
c including displacements and their locally Cartesian derivatives
c for PSV synthetic seismograms in transversely isotropic media 
c                                                     2004.5 K.Kawai
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c ---------------------------<< constants >>---------------------------
      real*8 pi
      integer maxnlay,maxnslay,maxnllay,ilog,iup
      integer maxnzone,maxnr,maxnstack,maxlmax
      parameter ( pi=3.1415926535897932d0 )
      parameter ( maxnlay = 20220 )
      parameter ( maxnslay = 12210 )
      parameter ( maxnllay = 8010 )
      parameter ( maxnzone = 15 )
      parameter ( maxlmax = 80000 )
      parameter ( maxnr = 400 )    ! number of xy-perturbations
      parameter ( maxnstack = 20 ) ! number of z-perturbations
      parameter ( ilog = 0 )
      parameter ( iup = 1)      ! 0-both; 1-onlyLCD
c ---------------------------<< variables >>---------------------------
c variable for the trial function
      integer nnlayer,nlayer(maxnzone)
      integer nslay,nllay
      integer inlayer,jnlayer,jnslay,jnllay
      integer l,m
      real*8 ra(maxnlay+maxnzone+1),plm(3,0:3,maxnr)
      complex*16 bvec(3,-2:2,maxnr)
      complex*16 bvecdt(3,-2:2,maxnr),bvecdp(3,-2:2,maxnr)
c variable for the structure
      integer nzone,isl,ill,nsl,nll
      integer iphase(maxnzone),ndc,vnp
      real*8 rmin,rmax
      real*8 vrmin(maxnzone),vrmax(maxnzone)
      real*8 rrho(4,maxnzone)
      real*8 vpv(4,maxnzone),vph(4,maxnzone)
      real*8 vsv(4,maxnzone),vsh(4,maxnzone),eta(4,maxnzone)
      real*8 qmu(maxnzone),qkappa(maxnzone)
      real*8 vra(maxnlay+2*maxnzone+1)
      real*8 rho(maxnlay+2*maxnzone+1)
      real*8 kappa(maxnlay+2*maxnzone+1) 
      real*8 ecKx(maxnlay+2*maxnzone+1) !3*Kx=3A-4N
      real*8 ecKy(maxnlay+2*maxnzone+1) !3*Ky=3F+2N
      real*8 ecKz(maxnlay+2*maxnzone+1) !3*Kz=2F+C
      real*8 mu(maxnlay+2*maxnzone+1)
      real*8 ecL(maxnlay+2*maxnzone+1)
      real*8 ecN(maxnlay+2*maxnzone+1)
      real*8 rhoinv(maxnlay+2*maxnzone+1)
      real*8 kappainv(maxnlay+2*maxnzone+1)
      complex*16 coef1(maxnzone),coef2(maxnzone),coef(maxnzone)
c variable for the periodic range
      integer np,imin,imax
      real*8 tlen,omega,omegai
      complex*16 u(9,maxnstack,maxnr)
      complex*16 ulcd(27,maxnstack,maxnr) ! locally Cartesian derivatives
c variable for the station
      integer nr,ir
      real*8 theta(maxnr),phi(maxnr)
      real*8 lat(maxnr),lon(maxnr)
c variable for the matrix elements
      complex*16 a0(4,2*(2*(maxnslay+1)+(maxnllay+1)+maxnzone))
      complex*16 a1(4,2*(2*(maxnslay+1)+(maxnllay+1)+maxnzone))
      complex*16 a2(4,2*(2*(maxnslay+1)+(maxnllay+1)+maxnzone))
      complex*16 a( 4, 2*(maxnslay+1) + (maxnllay+1) )
      complex*16 c( 2, (maxnslay+1) + (maxnllay+1) )
      real*8  t( 8*maxnslay )
      real*8 h1x( 8*maxnslay ),h1y( 8*maxnslay ),h1z( 8*maxnslay )
      real*8 h2L( 8*maxnslay ),h2N( 8*maxnslay )
      real*8 h3ax( 8*maxnslay )
      real*8 h3ay( 8*maxnslay ),h3az( 8*maxnslay )
      real*8 h4aL( 8*maxnslay ),h4aN( 8*maxnslay )
      real*8 h5ax( 8*maxnslay ),h5ay( 8*maxnslay ),h5az( 8*maxnslay )
      real*8 h6aL( 8*maxnslay ),h6aN( 8*maxnslay )
      real*8 h3x( 8*maxnslay ),h3y( 8*maxnslay ),h3z( 8*maxnslay )
      real*8 h4L( 8*maxnslay ),h4N( 8*maxnslay )
      real*8 h5x( 8*maxnslay ),h5y( 8*maxnslay ),h5z( 8*maxnslay )
      real*8 h6L( 8*maxnslay ),h6N( 8*maxnslay )
      real*8 h7x( 8*maxnslay ),h7y( 8*maxnslay ),h7z( 8*maxnslay )
      real*8 h8L( 8*maxnslay ),h8N( 8*maxnslay )
      real*8  h3mx( -2:1,2*(maxnslay+maxnzone) )
      real*8  h3my( -2:1,2*(maxnslay+maxnzone) )
      real*8  h3mz( -2:1,2*(maxnslay+maxnzone) )
      real*8  h5mx( -1:2,2*(maxnslay+maxnzone) )
      real*8  h5my( -1:2,2*(maxnslay+maxnzone) )
      real*8  h5mz( -1:2,2*(maxnslay+maxnzone) )
      real*8 h4m1L( -1:2,2*(maxnslay+maxnzone) )
      real*8 h4m1N( -1:2,2*(maxnslay+maxnzone) )
      real*8 h4m2L( -2:1,2*(maxnslay+maxnzone) )
      real*8 h4m2N( -2:1,2*(maxnslay+maxnzone) )
      real*8 h6m1L( -1:2,2*(maxnslay+maxnzone) )
      real*8 h6m1N( -1:2,2*(maxnslay+maxnzone) )
      real*8 h6m2L( -2:1,2*(maxnslay+maxnzone) )
      real*8 h6m2N( -2:1,2*(maxnslay+maxnzone) )
      real*8 p1( 8*maxnllay ),p2( 8*maxnllay ),p3( 8*maxnllay )
      complex*16 g( 2*(maxnslay+1) + (maxnllay+1) )
      complex*16 d( (maxnslay+1) + (maxnllay+1) )
c variable for the file
      character*80 outu(maxnr),outlcd(maxnr)
      character*80 dir
      character*80 outfile
c variable for the stack point
      integer isp(maxnzone)
      integer issp(maxnzone)
      integer ilsp(maxnzone),jssp(maxnzone)
      integer jsp(maxnzone)
      integer ksp(maxnzone),lsp(maxnzone)
      integer isdr,jsdr,ildr,cista,cksta
c variables for the output stack point
      integer nsta
      real*8 rsta(maxnstack),rrsta(3,maxnstack)
      integer istazone(maxnstack)
      integer iista(3,maxnstack) ! output stack point for layer
      integer ksta(maxnstack)   ! output stack point for g
      integer jsta(maxnstack)   ! output stack point for d
      complex*16 gtmp(2),gdertmp(2),gtmp0(1),gdertmp0(1)
      integer jjj
c variable for ths observer
      real*8 obslat,obslon
      character*80 obs
c variables for the gridding
      integer jjdr(maxnzone),kkdr(maxnzone)
      integer jdr,kdr
      real*8 vmin(maxnzone),gridpar(maxnzone),dzpar(maxnzone)
      real*8 re,ratc,ratl
c variables for l cut off
      integer kc,lsuf,sufzone,ismall,llog
      real*8 maxamp
c other variables
      integer i,j,jj,nn,lda,ier,itmp,jtmp,mtmp,nn0,ista
      integer ll(12),lli(12),llj(12)
      real*8 eps,work(8*maxnslay),l2,lsq
      complex*16 z( 2*(maxnslay+1) + (maxnllay+1) )
      complex*16 w( 2*(maxnslay+1) + (maxnllay+1) )
      complex*16 cwork( 2*(16*maxnslay + 4*maxnllay) )
      complex*16 ctmp( 2, (maxnslay+1) + (maxnllay+1) )
c     
      data lda/ 4 /
      data eps/ -1.d0 /
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c *************** Inputting and computing the parameters ***************
c --- inputting parameter ---
      call pinput(maxnlay,maxnzone,maxnr,maxnstack,
     &     re,ratc,ratl,
     &     tlen,np,omegai,imin,imax,
     &     nzone,vrmin,vrmax,rrho,
     &     vpv,vph,vsv,vsh,eta,qmu,qkappa,
     &     nr,theta,phi,lat,lon,dir,
     &     obslat,obslon,obs, nsta,rsta ) 
c --- computing the required parameters ---
c counting of the nsl and nll
      call calnl( nzone,vsv,iphase,nsl,nll )
c computing and checking the parameters
      rmin = vrmin(1)
      rmax = vrmax(nzone)
      ndc = nzone - 1
      do ir=1,nr
         theta(ir)= theta(ir) / 1.8d2 * pi
         phi(ir)= phi(ir)   / 1.8d2 * pi
      enddo
c
c ************************** Files Handling **************************
      call makeoutputfile(dir,obs,nr,outu,outlcd)
      do ir=1,nr
         if(iup.eq.0) then
            open( unit=11,file=outu(ir) )
            write(11,*) tlen
            write(11,*) np,nsta
            write(11,*) omegai,lat(ir),lon(ir)
c     write(11,*) theta(ir)*1.8d2/pi,phi(ir)*1.8d2/pi
            write(11,*) obslat,obslon
            do ista=1,nsta
               write(11,*) rsta(ista)
            enddo
            close(11)
         endif
         open( unit=11,file=outlcd(ir) )
         write(11,*) tlen
         write(11,*) np,nsta
         write(11,*) omegai,lat(ir),lon(ir)
c         write(11,*) theta(ir)*1.8d2/pi,phi(ir)*1.8d2/pi
         write(11,*) obslat,obslon
         do ista=1,nsta
            write(11,*) rsta(ista)
         enddo
         close(11)
      enddo
      if(ilog.eq.1) then
         open(unit=11,file='llog.log',status='unknown')
         close(11)
      endif
c ************************** Files Handling end **********************
c
c computing of the number and the location of grid points
      call calgrid( nzone,vrmin,vrmax,vpv,vsv,rmin,rmax,
     &     imax,1,tlen,
     &     vmin,gridpar,dzpar )
      call calra( maxnlay,maxnslay,maxnllay,maxnzone,maxnstack,
     &     nnlayer,inlayer,jnlayer,jnslay,jnllay,
     &     gridpar,dzpar,nzone,vrmin,vrmax,iphase,
     &     rmin,rmax,nslay,nllay,nlayer,ra,re,
     &     nsta,rsta,rrsta,istazone,iista,rmax,cista)
c
c --- checking the parameter
      if ( inlayer.gt.maxnlay )
     &     pause 'The number of grid points is too large.'
      if ( nslay.gt.maxnslay )
     &     pause
     &     'The number of the grid points in the solid is too large.'
      if ( nllay.gt.maxnllay )
     &     pause
     &     'The number of the grid points in the liquid is too large.'
      if ( jnlayer.gt.2*maxnlay )
     &     pause 'The number of the total grid points is too large.'
      if ( ( jnslay.gt.2*maxnslay ).or.( jnllay.gt.2*maxnllay ) )
     &     pause 'The number of the total grid points is too large.'
c computing the stack points
      call calsp( maxnzone,ndc,nsl,nll,
     &     iphase,nlayer,nslay,nllay,
     &     isp,jsp,ksp,issp,ilsp,lsp,jssp,
     &        isdr,jsdr,ildr,jdr,kdr )
c ******************* Computing the matrix elements *******************
c data initialization
      call xmatinit( maxnslay,maxnllay,maxnzone,
     &     a,t,h1x,h1y,h1z,h2L,h2N,
     &     h3ax,h3ay,h3az,h4aL,h4aN,
     &     h5ax,h5ay,h5az,h6aL,h6aN,
     &     h3x,h3y,h3z,h4L,h4N,
     &     h5x,h5y,h5z,h6L,h6N,
     &     h7x,h7y,h7z,h8L,h8N,
     &     h3mx,h3my,h3mz,h5mx,h5my,h5mz,
     &     h4m1L,h4m1N,h4m2L,h4m2N,
     &     h6m1L,h6m1N,h6m2L,h6m2N,
     &     p1,p2,p3 )
c computing the structure grid points
      call calstg( maxnlay,maxnzone,
     &     nzone,iphase,rrho,
     &     vpv,vph,vsv,vsh,eta,nlayer,ra,rmax,
     &     vnp,vra,rho,kappa,ecKx,ecKy,ecKz,
     &     mu,ecL,ecN )
      call calinv( vnp,rho,kappa,rhoinv,kappainv )
      isl = 0
      ill = 0
      do i=1,ndc+1
         if ( iphase(i).eq.1 ) then
            isl = isl + 1
            itmp = isdr+issp(isl)
            call calmatc( nlayer(i),vnp,vra,
     &           rho ,2,0,0,ra(isp(i)), t(itmp) )
            call caltl( nlayer(i),vnp,vra,rho,
     &           ra(isp(i)),work(itmp) )
            call calt( nlayer(i), t(itmp),work(itmp), t(itmp) )
            call calmatc( nlayer(i),vnp,vra,
     &           ecKx,0,0,0,ra(isp(i)),h1x(itmp) )
            call calhl( nlayer(i),vnp,vra,ecKx,
     &           ra(isp(i)),work(itmp) )
            call calt( nlayer(i), h1x(itmp),work(itmp),h1x(itmp))
            call calmatc( nlayer(i),vnp,vra,
     &           ecKy,0,0,0,ra(isp(i)),h1y(itmp) )
            call calhl( nlayer(i),vnp,vra,ecKy,
     &           ra(isp(i)),work(itmp) )
            call calt( nlayer(i), h1y(itmp),work(itmp),h1y(itmp))
            call calmatc( nlayer(i),vnp,vra,
     &           ecKz,0,0,0,ra(isp(i)),h1z(itmp) )
            call calhl( nlayer(i),vnp,vra,ecKz,
     &           ra(isp(i)),work(itmp) )
            call calt( nlayer(i), h1z(itmp),work(itmp),h1z(itmp))
            call calmatc( nlayer(i),vnp,vra,
     &           ecL ,0,0,0,ra(isp(i)),h2L(itmp) )
            call calhl( nlayer(i),vnp,vra,ecL,
     &           ra(isp(i)),work(itmp) )
            call calt( nlayer(i), h2L(itmp),work(itmp),h2L(itmp))
            call calmatc( nlayer(i),vnp,vra,
     &           ecN ,0,0,0,ra(isp(i)),h2N(itmp) )
            call calhl( nlayer(i),vnp,vra,ecN,
     &           ra(isp(i)),work(itmp) )
            call calt( nlayer(i), h2N(itmp),work(itmp),h2N(itmp))
            call calmatc( nlayer(i),vnp,vra,
     &           ecKx,1,0,1,ra(isp(i)),h5ax(itmp) )
            call calmatc( nlayer(i),vnp,vra,
     &           ecKy,1,0,1,ra(isp(i)),h5ay(itmp) )
            call calmatc( nlayer(i),vnp,vra,
     &           ecKz,1,0,1,ra(isp(i)),h5az(itmp) )
            call calmatc( nlayer(i),vnp,vra,
     &           ecL ,1,0,1,ra(isp(i)),h6aL(itmp) )
            call calmatc( nlayer(i),vnp,vra,
     &           ecN ,1,0,1,ra(isp(i)),h6aN(itmp) )
            call mtrnp( nlayer(i),h5ax(itmp),h3ax(itmp) )
            call mtrnp( nlayer(i),h5ay(itmp),h3ay(itmp) )
            call mtrnp( nlayer(i),h5az(itmp),h3az(itmp) )
            call mtrnp( nlayer(i),h6aL(itmp),h4aL(itmp) )
            call mtrnp( nlayer(i),h6aN(itmp),h4aN(itmp) )
            call calmatc( nlayer(i),vnp,vra,
     &           ecKx,2,1,1,ra(isp(i)), h7x(itmp) )
            call calmatc( nlayer(i),vnp,vra,
     &           ecKy,2,1,1,ra(isp(i)), h7y(itmp) )
            call calmatc( nlayer(i),vnp,vra,
     &           ecKz,2,1,1,ra(isp(i)), h7z(itmp) )
            call calmatc( nlayer(i),vnp,vra,
     &           ecL ,2,1,1,ra(isp(i)), h8L(itmp) )
            call calmatc( nlayer(i),vnp,vra,
     &           ecN ,2,1,1,ra(isp(i)), h8N(itmp) )
         else
            ill = ill + 1
            itmp = ildr+ilsp(ill)
            call calmatc( nlayer(i),vnp,vra,
     &           rhoinv,2,1,1,ra(isp(i)),p1(itmp) )
            call calmatc( nlayer(i),vnp,vra,
     &           rhoinv,0,0,0,ra(isp(i)),p2(itmp) )
            call calhl( nlayer(i),vnp,vra,
     &           rhoinv,ra(isp(i)),work(itmp) )
            call calt( nlayer(i),p2(itmp),work(itmp),p2(itmp) )
            call calmatc( nlayer(i),vnp,vra,
     &           kappainv,2,0,0,ra(isp(i)),p3(itmp) )
            call caltl( nlayer(i),vnp,vra,
     &           kappainv,ra(isp(i)),work(itmp) )
            call calt( nlayer(i),p3(itmp),work(itmp),p3(itmp) )
         endif
      enddo
c Computing the modified operator of the 1st derivative
      call caltstg( maxnlay,maxnzone,
     &     nzone,rrho,vpv,vph,vsv,vsh,eta,
     &     nlayer,ra,rmax,
     &     vra,kappa,ecKx,ecKy,ecKz,mu,ecL,ecN )
      isl = 0
      do i=1,ndc+1
         if ( iphase(i).eq.1 ) then
            isl = isl + 1
            itmp = isdr+issp(isl)
            jtmp = isp(i)+i-1
            call calh5( nlayer(i),vra(jtmp),
     &           ecKx(jtmp),ra(isp(i)),work(itmp) )
            call submat( nlayer(i),h5ax(itmp),
     &           work(itmp),h5x(itmp) )
            call calh5( nlayer(i),vra(jtmp),
     &           ecKy(jtmp),ra(isp(i)),work(itmp) )
            call submat( nlayer(i),h5ay(itmp),
     &           work(itmp),h5y(itmp) )
            call calh5( nlayer(i),vra(jtmp),
     &           ecKz(jtmp),ra(isp(i)),work(itmp) )
            call submat( nlayer(i),h5az(itmp),
     &           work(itmp),h5z(itmp) )
            call calh5( nlayer(i),vra(jtmp),
     &           ecL(jtmp),ra(isp(i)),work(itmp) )
            call submat( nlayer(i),h6aL(itmp),
     &           work(itmp),h6L(itmp) )
            call calh5( nlayer(i),vra(jtmp),
     &           ecN(jtmp),ra(isp(i)),work(itmp) )
            call submat( nlayer(i),h6aN(itmp),
     &           work(itmp),h6N(itmp) )
            call mtrnp( nlayer(i),h5x(itmp),h3x(itmp) )
            call mtrnp( nlayer(i),h5y(itmp),h3y(itmp) )
            call mtrnp( nlayer(i),h5z(itmp),h3z(itmp) )
            call mtrnp( nlayer(i),h6L(itmp),h4L(itmp) )
            call mtrnp( nlayer(i),h6N(itmp),h4N(itmp) )
            itmp = jsdr+jssp(isl)
            call calhm1( nlayer(i),vra(jtmp),
     &           ecKx(jtmp),ra(isp(i)),h5mx(-1,itmp) )
            call calhm1( nlayer(i),vra(jtmp),
     &           ecKy(jtmp),ra(isp(i)),h5my(-1,itmp) )
            call calhm1( nlayer(i),vra(jtmp),
     &           ecKz(jtmp),ra(isp(i)),h5mz(-1,itmp) )
            call calhm1( nlayer(i),vra(jtmp),
     &           ecL(jtmp),ra(isp(i)),h6m1L(-1,itmp) )
            call calhm1( nlayer(i),vra(jtmp),
     &           ecN(jtmp),ra(isp(i)),h6m1N(-1,itmp) )
            call calhm2( nlayer(i),vra(jtmp),
     &           ecL(jtmp),ra(isp(i)),h6m2L(-2,itmp) )
            call calhm2( nlayer(i),vra(jtmp),
     &           ecN(jtmp),ra(isp(i)),h6m2N(-2,itmp) )
            call mtrnp2( nlayer(i),1,2,h5mx(-1,itmp),h3mx(-2,itmp) )
            call mtrnp2( nlayer(i),1,2,h5my(-1,itmp),h3my(-2,itmp) )
            call mtrnp2( nlayer(i),1,2,h5mz(-1,itmp),h3mz(-2,itmp) )
            call mtrnp2( nlayer(i),1,2,h6m1L(-1,itmp),h4m2L(-2,itmp) )
            call mtrnp2( nlayer(i),1,2,h6m1N(-1,itmp),h4m2N(-2,itmp) )
            call mtrnp2( nlayer(i),2,1,h6m2L(-2,itmp),h4m1L(-1,itmp) )
            call mtrnp2( nlayer(i),2,1,h6m2N(-2,itmp),h4m1N(-1,itmp) )
         endif
      enddo
c     
c ******************** Computing the displacement *********************
c
      llog = 0
      do i=imin,imax            ! f-loop start
         do ir=1,nr
            call cmatinit( 9,maxnstack,u(1,1,ir) )
            call cmatinit( 27,maxnstack,ulcd(1,1,ir) )
         enddo
         if ( i.ne.0 ) then
	    omega = 2.d0 * pi * dble(i) / tlen
	    call callsuf(omega,nzone,vrmax,vsv,lsuf)
	    do ir=1,nr
               call matinit( 3,4,plm(1,0,ir) )
            enddo
            call calcoef( nzone,omega,qmu,qkappa,coef1,coef2,coef )
c computing the matrix elements independent of l
            isl = 0
            ill = 0
            do j=1,ndc+1
               if ( iphase(j).eq.1 ) then
                  isl = isl + 1
                  itmp = isdr+issp(isl)
                  jtmp = jdr+jsp(j)
                  mtmp = kdr+ksp(j)
                  call cala0( nlayer(j),omega,omegai,
     &                 t(itmp), 
     &                 h1x(itmp), h1y(itmp),h1z(itmp),
     &                 h2L(itmp), h2N(itmp),
     &                 h3ax(itmp),h3ay(itmp),h3az(itmp),
     &                 h4aL(itmp),h4aN(itmp),
     &                 h5ax(itmp),h5ay(itmp),h5az(itmp),
     &                 h6aL(itmp),h6aN(itmp),
     &                 h7x(itmp),h7y(itmp),h7z(itmp),
     &                 h8L(itmp), h8N(itmp),
     &                 coef1(j),coef2(j),cwork(jtmp) )
                  call overlapa( nlayer(j),
     &                 cwork(jtmp),a0(1,mtmp))
                  call cala1( nlayer(j),
     &                 h1x(itmp),h1y(itmp),h1z(itmp),
     &                 h2L(itmp),h2N(itmp),
     &                 h3x(itmp), h3y(itmp), h3z(itmp), 
     &                 h4L(itmp), h4N(itmp), 
     &                 h5x(itmp), h5y(itmp), h5z(itmp), 
     &                 h6L(itmp), h6N(itmp),
     &                 coef1(j),coef2(j),cwork(jtmp) )
                  call overlapa( nlayer(j),
     &                 cwork(jtmp),a1(1,mtmp))
                  call cala2( nlayer(j), 
     &                 h1x(itmp), h1y(itmp),h1z(itmp),
     &                 h2L(itmp),h2N(itmp),
     &                 coef1(j),coef2(j),cwork(jtmp) )
                  call overlapa( nlayer(j),
     &                 cwork(jtmp),a2(1,mtmp))
                  jtmp = jsdr+jssp(isl)
                  call calhml( nlayer(j),coef1(j),coef2(j),
     &                 h3mx(-2,jtmp),h3my(-2,jtmp),h3mz(-2,jtmp),
     &                 h5mx(-1,jtmp),h5my(-1,jtmp),h5mz(-1,jtmp),
     &                 h4m1L(-1,jtmp),h4m1N(-1,jtmp),
     &                 h4m2L(-2,jtmp),h4m2N(-2,jtmp),
     &                 h6m1L(-1,jtmp),h6m1N(-1,jtmp),
     &                 h6m2L(-2,jtmp),h6m2N(-2,jtmp),
     &                 a1(1,mtmp) )
               else
                  ill = ill + 1
                  itmp = ildr+ilsp(ill)
                  jtmp = jdr+jsp(j)
                  mtmp = kdr+ksp(j)
                  call calb0( nlayer(j),omega,omegai,
     &                 p1(itmp),p3(itmp),coef(j),cwork(jtmp) )
                  call overlapb( nlayer(j),
     &                 cwork(jtmp),a0(1,mtmp))
                  call calb2( nlayer(j),omega,omegai,
     &                 p2(itmp),coef(j),cwork(jtmp) )
                  call overlapb( nlayer(j),
     &                 cwork(jtmp),a2(1,mtmp))
               endif
            enddo
c     
            kc = 1
            ismall = 0
            maxamp = -1.d0
            llog = maxlmax
            do l=0,maxlmax      ! l-loop start
               if( ismall.gt.20 ) then
                  if(llog.gt.l) llog = l
                  exit
               endif
               l2 = dble(l)*dble(l+1)
               lsq = dsqrt( l2 )
c computing the coefficient matrix elements
c --- renewing  mdr
               if ( mod(l,50).eq.0  ) then
                  call calmdr( omega,l,nzone,vrmin,vrmax,vmin,dzpar,
     &                 rmax,sufzone )
                  call calspdr( maxnzone,nzone,iphase,nlayer,
     &                 jjdr,kkdr )
                  do ista=1,nsta
                     ksta(ista) = kkdr(istazone(ista)) 
     &                    + 2 * iista(1,ista) - 1
                  enddo
                  cksta = kkdr(istazone(cista)) 
     &                 + 2 * iista(1,cista) - 1
                  nn = kkdr(nzone) + 2 * nlayer(nzone) + 1
               endif
c ***** Computing the trial function *****
               do ir=1,nr
                  call calbvec( l,theta(ir),phi(ir),
     &                 plm(1,0,ir),bvec(1,-2,ir),
     &                 bvecdt(1,-2,ir),bvecdp(1,-2,ir))
               enddo
c     computing the matrix elements
               call cala( maxnzone,ndc,iphase,nlayer,kkdr,kdr,ksp,
     &              l2,lsq,nn,a0,a1,a2,a )
c computing the boundary condition elements
               call calbc( maxnzone,ndc,vrmax,iphase,kkdr,a )
c     
               do m=-2,2        ! m-loop start---------------
                  if ( iabs(m).le.iabs(l) ) then
                     if ( l.eq.0 ) then
c     ---f_r---
                        call cvecinit( nn,g )
c computing the excitation vector
                        call calgpfpsv(l,m,1,g(nn-1))
c rearranging the matrix elements for l=0
                        call rea2( nn,a,g,c,d,
     &                       nzone,iphase,kkdr,nn0,
     &                       maxnstack,nsta,istazone,iista,jsta )
c computing the expansion coefficient
                        itmp=1
                        if ( rmin.eq.0.d0 ) itmp=2
                        call dcsymbdl0( c(1,itmp),1,nn0-itmp+1,1,eps,
     &                       z(itmp),w(itmp),ll,lli,llj,ier )
                        call dcsbdlv0( c(1,itmp),d(itmp),1,
     &                       nn0-itmp+1,eps,z(itmp),ier )
c computing the displacement
                        do ir=1,nr
                           do ista=1,nsta
                              call interpolate( 1,0,rsta(ista),
     &                             rrsta(1,ista),
     &                             d(jsta(ista)),gtmp0 )
                              call interpolate( 1,1,rsta(ista),
     &                             rrsta(1,ista),
     &                             d(jsta(ista)),gdertmp0 )
                              call calu0( gtmp0(1),bvec(1,m,ir),
     &                             u(1,ista,ir))
                              call calulcd0(gtmp0(1),gdertmp0(1),
     &                             rsta(ista),theta(ir),
     &                             bvec(1,m,ir),
     &                             bvecdt(1,m,ir),bvecdp(1,m,ir),
     &                             ulcd(1,ista,ir) )
                           enddo
                        enddo
c     ---f_\theta---
                        call cvecinit( nn,g )
c computing the excitation vector
                        call calgpfpsv(l,m,2,g(nn-1))
c rearranging the matrix elements for l=0
                        call rea2( nn,a,g,ctmp,d,
     &                       nzone,iphase,kkdr,nn0,
     &                       maxnstack,nsta,istazone,iista,jsta )
c computing the expansion coefficient
                        itmp=1
                        if ( rmin.eq.0.d0 ) itmp=2
                        call dcsbdlv0( c(1,itmp),d(itmp),1,
     &                       nn0-itmp+1,eps,z(itmp),ier )
c computing the displacement
                        do ir=1,nr
                           do ista=1,nsta
                              call interpolate( 1,0,rsta(ista),
     &                             rrsta(1,ista),
     &                             d(jsta(ista)),gtmp0 )
                              call interpolate( 1,1,rsta(ista),
     &                             rrsta(1,ista),
     &                             d(jsta(ista)),gdertmp0 )
                              call calu0( gtmp0(1),bvec(1,m,ir),
     &                             u(4,ista,ir))
                              call calulcd0(gtmp0(1),gdertmp0(1),
     &                             rsta(ista),theta(ir),
     &                             bvec(1,m,ir),
     &                             bvecdt(1,m,ir),bvecdp(1,m,ir),
     &                             ulcd(10,ista,ir) )
                           enddo
                        enddo
c     ---f_\phi---
                        call cvecinit( nn,g )
c computing the excitation vector
                        call calgpfpsv(l,m,3,g(nn-1))
c rearranging the matrix elements for l=0
                        call rea2( nn,a,g,ctmp,d,
     &                       nzone,iphase,kkdr,nn0,
     &                       maxnstack,nsta,istazone,iista,jsta )
c computing the expansion coefficient
                        itmp=1
                        if ( rmin.eq.0.d0 ) itmp=2
                        call dcsbdlv0( c(1,itmp),d(itmp),1,
     &                       nn0-itmp+1,eps,z(itmp),ier )
c computing the displacement
                        do ir=1,nr
                           do ista=1,nsta
                              call interpolate( 1,0,rsta(ista),
     &                             rrsta(1,ista),
     &                             d(jsta(ista)),gtmp0 )
                              call interpolate( 1,1,rsta(ista),
     &                             rrsta(1,ista),
     &                             d(jsta(ista)),gdertmp0 )
                              call calu0( gtmp0(1),bvec(1,m,ir),
     &                             u(7,ista,ir))
                              call calulcd0(gtmp0(1),gdertmp0(1),
     &                             rsta(ista),theta(ir),
     &                             bvec(1,m,ir),
     &                             bvecdt(1,m,ir),bvecdp(1,m,ir),
     &                             ulcd(19,ista,ir) )
                           enddo
                        enddo
                     else       ! l-branch for cal u (l!=0)
c     ---f_r---
                        call cvecinit( nn,g )
c computing the excitation vector
                        call calgpfpsv(l,m,1,g(nn-1))
c computing the expansion coefficient
                        itmp=1
                        if ( rmin.eq.0.d0 ) itmp=3
                        if ( ( m.eq.-2 ).or.( m.eq.-l ) ) then
                           call dcsymbdl0( a(1,itmp),3,nn-itmp+1,6,eps,
     &                          z(itmp),w(itmp),ll,lli,llj,ier )
                        endif
                        call dcsbdlv0( a(1,itmp),g(itmp),3,
     &                       nn-itmp+1,eps,z(itmp),ier )
c computing ampratio
                        call calamp( g(cksta),l,lsuf,maxamp,ismall,ratl)
c     
c computing the displacement
                        do ir=1,nr
                           do ista=1,nsta
                              call interpolate( 2,0,rsta(ista),
     &                             rrsta(1,ista),
     &                             g(ksta(ista)-1),gtmp )
                              call interpolate( 2,1,rsta(ista),
     &                             rrsta(1,ista),
     &                             g(ksta(ista)-1),gdertmp )
                              call calu( gtmp,lsq,
     &                             bvec(1,m,ir),u(1,ista,ir) )
                              call calulcd( gtmp,gdertmp,lsq,
     &                             rsta(ista),theta(ir),
     &                             bvec(1,m,ir),
     &                             bvecdt(1,m,ir),bvecdp(1,m,ir),
     &                             ulcd(1,ista,ir) )
                           enddo
                        enddo
c     ---f_\theta---
                        call cvecinit( nn,g )
c computing the excitation vector
                        call calgpfpsv(l,m,2,g(nn-1))
c computing the expansion coefficient
                        itmp=1
                        if ( rmin.eq.0.d0 ) itmp=3
                        call dcsbdlv0( a(1,itmp),g(itmp),3,
     &                       nn-itmp+1,eps,z(itmp),ier )
c computing the displacement
                        do ir=1,nr
                           do ista=1,nsta
                              call interpolate( 2,0,rsta(ista),
     &                             rrsta(1,ista),
     &                             g(ksta(ista)-1),gtmp )
                              call interpolate( 2,1,rsta(ista),
     &                             rrsta(1,ista),
     &                             g(ksta(ista)-1),gdertmp )
                              call calu( gtmp,lsq,
     &                             bvec(1,m,ir),u(4,ista,ir) )
                              call calulcd( gtmp,gdertmp,lsq,
     &                             rsta(ista),theta(ir),
     &                             bvec(1,m,ir),
     &                             bvecdt(1,m,ir),bvecdp(1,m,ir),
     &                             ulcd(10,ista,ir) )
                           enddo
                        enddo
c     ---f_\phi---
                        call cvecinit( nn,g )
c computing the excitation vector
                        call calgpfpsv(l,m,3,g(nn-1))
c computing the expansion coefficient
                        itmp=1
                        if ( rmin.eq.0.d0 ) itmp=3
                        call dcsbdlv0( a(1,itmp),g(itmp),3,
     &                       nn-itmp+1,eps,z(itmp),ier )
c computing the displacement
                        do ir=1,nr
                           do ista=1,nsta
                              call interpolate( 2,0,rsta(ista),
     &                             rrsta(1,ista),
     &                             g(ksta(ista)-1),gtmp )
                              call interpolate( 2,1,rsta(ista),
     &                             rrsta(1,ista),
     &                             g(ksta(ista)-1),gdertmp )
                              call calu( gtmp,lsq,
     &                             bvec(1,m,ir),u(7,ista,ir) )
                              call calulcd( gtmp,gdertmp,lsq,
     &                             rsta(ista),theta(ir),
     &                             bvec(1,m,ir),
     &                             bvecdt(1,m,ir),bvecdp(1,m,ir),
     &                             ulcd(19,ista,ir) )
                           enddo
                        enddo
                    endif      ! l-branch for calu
                  endif
               enddo            ! m-loop end--------------------
            enddo               ! l-loop end
         endif
c ************************** Files Handling **************************
         do ir=1,nr
            if(iup.eq.0) then
               open( unit=11,file=outu(ir),
     &              position='append',status='old')
            endif
            open( unit=12,file=outlcd(ir),
     &           position='append',status='old')
            do ista=1,nsta
               if(iup.eq.0) then
                  write(11,*) i,
     &                 dble(u(1,ista,ir)),dimag(u(1,ista,ir))
                  do jj=2,9
                     write(11,*)  dble(u(jj,ista,ir)),
     &                    dimag(u(jj,ista,ir))
                  enddo
               endif
               write(12,*) i,
     &              dble(ulcd(1,ista,ir)),dimag(ulcd(1,ista,ir))
               do jj=2,27
                  write(12,*) dble(ulcd(jj,ista,ir)),
     &                 dimag(ulcd(jj,ista,ir))
               enddo
            enddo
            if(iup.eq.0) close(11)
            close(12)
         enddo
         if(ilog.eq.1) then
            open(unit=11,file='llog.log',access='append',status='old')
            write(11,*) i,llog,inlayer
            close(11)
         endif
c**************************Files Handling end **********************
      enddo                     ! f-loop end
      end
