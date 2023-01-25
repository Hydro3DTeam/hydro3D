!#######################################################################
      DOUBLE PRECISION FUNCTION dh(dx,dy,dz,xij,yij,zij,Xl,Yl,Zl,order)
!-----------------------------------------------------------------------
!     The delta interpolation function menu to select the appriopriate
!     interpolation order. This case number is prescribed in in_geom.cin
!
!     PAPER: Yang,X., Zhang,X., Li,Z., & He, G.W.(2009). A smoothing
!     technique for discrete delta functions with application to 
!     immersed boundary method in moving boundary simulations. JCP
!#######################################################################
      implicit none 

      INTEGER,intent(in) :: order
      DOUBLE PRECISION,intent(in) :: dx,dy,dz,xij,yij,zij,Xl,Yl,Zl
      DOUBLE PRECISION :: phi_r2smth,phi_r3smth,phi_r4
      DOUBLE PRECISION :: phi_r3,phi_r1smth,phi_r4smth

      SELECT CASE (order)

         CASE (1)
      dh =  phi_r1smth((xij-Xl)/dx) 
     & * phi_r1smth((yij-Yl)/dy) * phi_r1smth((zij-Zl)/dz)
         CASE (2)
      dh =   phi_r2smth((xij-Xl)/dx) 
     & * phi_r2smth((yij-Yl)/dy) * phi_r2smth((zij-Zl)/dz)
         CASE (3)
       dh =   phi_r3smth((xij-Xl)/dx) 
     &   * phi_r3smth((yij-Yl)/dy) * phi_r3smth((zij-Zl)/dz)
         CASE (4)
       dh =   phi_r4smth((xij-Xl)/dx) 
     &   * phi_r4smth((yij-Yl)/dy) * phi_r4smth((zij-Zl)/dz)
         CASE (5)
       dh =   phi_r3((xij-Xl)/dx) 
     &   * phi_r3((yij-Yl)/dy) * phi_r3((zij-Zl)/dz)
         CASE (6)
      dh =   phi_r4((xij-Xl)/dx)
     &   * phi_r4((yij-Yl)/dy) * phi_r4((zij-Zl)/dz)
      Case default

        WRITE(6,'(a)')  '===ERROR==='
        WRITE(6,'(a)')' order of delta function is not selected '
        STOP

      END SELECT

      RETURN
      END FUNCTION dh
!#######################################################################
      DOUBLE PRECISION FUNCTION phi_r1smth(r)
!-----------------------------------------------------------------------
!     Case (1)
!     Smooth 2-point function (phi1* - Yang et al, 2009 - Eq.18)
!#######################################################################
      implicit none 
      
      DOUBLE PRECISION,intent(in) :: r
      DOUBLE PRECISION :: PI,abr
      
      PI = 4.D0*DATAN(1.D0)
      abr=SQRT(r*r)
      
      IF (abr.ge.1.5) THEN
        phi_r1smth = 0.0
      ELSE IF ((abr.lt.1.5).and.(abr.ge.0.5)) THEN
        phi_r1smth = 9./8.-3.*abr/2+abr**2/2
      ELSE IF ((abr.lt.0.5).and.(abr.ge.0.0)) THEN
        phi_r1smth = 3./4.-abr**2
      END IF
      
      RETURN
      END FUNCTION phi_r1smth
!#######################################################################
      DOUBLE PRECISION FUNCTION phi_r2smth(r)
!-----------------------------------------------------------------------
!     Case (2)
!     Smooth 4-point cosine function (phi2* - Yang et al, 2009 - Eq.19)
!#######################################################################
      implicit none 
      
      DOUBLE PRECISION, intent(in) :: r
      DOUBLE PRECISION :: PI
      
      PI = 4.D0*DATAN(1.D0)

      IF (r.le.-2.5) THEN
        phi_r2smth = 0.0
      ELSE IF ((r.ge.-2.5).and.(r.le.-1.5)) THEN
        phi_r2smth= -1./8./PI*(-5.*PI-2.*PI*r+4.*sin(PI/4.*(-2.*r-1.)))
      ELSE IF ((r.ge.-1.5).and.(r.le.0.0)) THEN
        phi_r2smth = 1./4./PI*(PI+2.*sin(PI/4.*(-2.*r+1.))
     &                           -2.*sin(PI/4.*(-2.*r-1.)))
      ELSE IF ((r.ge.0.0).and.(r.le.1.5)) THEN
        phi_r2smth = 1./4./PI*(PI+2.*sin(PI/4.*(2.*r+1.))
     &                           -2.*sin(PI/4.*(2.*r-1.)))
      ELSE IF ((r.ge.1.5).and.(r.le.2.5)) THEN
        phi_r2smth= -1./8./PI*(-5.*PI+2.*PI*r+4.*sin(PI/4.*(2.*r-1.)))
      ELSE IF (r.ge.2.5) THEN
        phi_r2smth = 0.0
      END IF

      RETURN
      END FUNCTION phi_r2smth
!#######################################################################
      DOUBLE PRECISION FUNCTION phi_r3smth(r)
!-----------------------------------------------------------------------
!     Case (3)
!     Smooth 2-point function (phi3* - Yang et al, 2009 - Eq.20)
!#######################################################################
      implicit none 

      DOUBLE PRECISION, intent(in) :: r
      DOUBLE PRECISION :: scal1,scal2,scal3,scal4,scal5,scal6

      scal1 = 1.095450017660160615261  ! 55.0/48.0 - sqrt(3.0)*pi/108.0
      scal2 = 1.083333333333333333333  ! 13.0/12.0
      scal3 = 0.4045499823398393847387 ! 17.0/48.0 + sqrt(3.0)*pi/108.0
      scal4 = 0.0481125224324688137091 ! sqrt(3.0)/36.0
      scal5 = 0.1443375672974064411273 ! sqrt(3.0)/12.0
      scal6 = 0.8660254037844386467637 ! sqrt(3.0)/2.0

      IF (r.le.-2.d0) THEN
        phi_r3smth = 0.0
      ELSE IF ((r.ge.-2.0).and.(r.le.-1.0)) THEN
        phi_r3smth = scal1 + scal2*r + scal4*ASIN(scal6*(-2.0*r-3.0))
     &+ 0.25*r**2 + (-2.0*r-3.0)/48.0*sqrt(-12.0*r**2-36.0*r-23.0)
      ELSE IF ((r.ge.-1.0).and.(r.le.0.0)) THEN
        phi_r3smth = scal3 - r/4.0 - scal5*ASIN(scal6*(-2.0*r-1.0))
     &- 0.25*r**2 + (2.0*r+1.0)/16.0*sqrt(-12.0*r**2-12.0*r+1.0)
      ELSE IF ((r.ge.0.0).and.(r.le.1.0)) THEN
        phi_r3smth = scal3 + r/4.0 - scal5*ASIN(scal6*(2.0*r-1.0))
     &- 0.25*r**2 + (-2.0*r+1.0)/16.0*sqrt(-12.0*r**2+12.0*r+1.0)
      ELSE IF ((r.ge.1.0).and.(r.le.2.0)) THEN
        phi_r3smth = scal1 - scal2*r + scal4*ASIN(scal6*(2.0*r-3.0))
     &+ 0.25*r**2 + (2.0*r-3.0)/48.0*sqrt(-12.0*r**2+36.0*r-23.0)
      ELSE IF (r.ge.2.0) THEN
        phi_r3smth = 0.0
      END IF

      RETURN
      END FUNCTION phi_r3smth
!#######################################################################
      Double precision function phi_r4smth(r)
!-----------------------------------------------------------------------
!     Case (4)
!     Smooth 4-point function (phi4 - Yang et al, 2009 - Eq.21)
!#######################################################################
      implicit none

      DOUBLE PRECISION, intent(in) :: r
      DOUBLE PRECISION :: PI,ar
      
      PI = 4.D0*DATAN(1.D0)
      ar=abs(r)
      
      IF ((ar.ge.0.0).and.(ar.le.0.5)) THEN
	phi_r4smth = 3.0/8.0+ PI/32.0 - r*r/4.0
      ELSE IF ((ar.ge.0.5).and.(ar.le.1.5)) THEN
	phi_r4smth = 1.0/4.0 + (1.0-ar)/8.0*SQRT(-2.0+8.0*ar -4.0*r*r)
     & -1.0/8.0*ASIN(sqrt(2.0)*(ar-1.0))
      ELSE IF ((ar.ge.1.5).and.(ar.le.2.5)) THEN
	phi_r4smth = 17.0/16.0 - PI/64.0 - 3.0*ar/4.0+ r*r/8.0
     & + (ar-2.0)/16.0 * SQRT(-14.0+ 16.0*ar - 4.0*r*r)
     & + 1.0/16.0*ASIN(sqrt(2.0)*(ar-2.0))
      ELSE IF (ar.ge.2.5) THEN
        phi_r4smth = 0.0
      END IF

      RETURN
      END FUNCTION phi_r4smth
!#######################################################################
      DOUBLE PRECISION FUNCTION phi_r3(r)
!-----------------------------------------------------------------------
!     Case (5)
!     3-point function (phi1* - Yang et al, 2009 - Eq.16)
!#######################################################################
      implicit none 
      
      DOUBLE PRECISION, intent(in) :: r

      IF (r.le.-1.5) THEN
        phi_r3 = 0.0
      ELSE IF ((r.ge.-1.5).and.(r.le.-0.5)) THEN
        phi_r3 = 1.0/6.0*(5.0+3.0*r-sqrt(-3.0*(1.0+r)**2+1.0))
      ELSE IF ((r.ge.-0.5).and.(r.le.0.0)) THEN
        phi_r3 = 1.0/3.0*(1.0+sqrt(-3.0*r**2+1.0))
      ELSE IF ((r.ge.0.0).and.(r.le.0.5)) THEN
        phi_r3 = 1.0/3.0*(1.0+sqrt(-3.0*r**2+1.0))
      ELSE IF ((r.ge.0.5).and.(r.le.1.5)) THEN
        phi_r3 = 1.0/6.0*(5.0-3.0*r-sqrt(-3.0*(1.0-r)**2+1.0))
      ELSE IF (r.ge.1.5) THEN
        phi_r3 = 0.0
      END IF

      RETURN
      END FUNCTION  phi_r3
!#######################################################################
      DOUBLE PRECISION FUNCTION phi_r4(r)
!-----------------------------------------------------------------------
!     Case (6)
!     4-point function (phi4 - Yang et al, 2009 - Eq.17)
!#######################################################################
      implicit none 
      
      DOUBLE PRECISION, intent(in) :: r
      
      IF (r.le.-2.0) THEN
        phi_r4 = 0.0
      ELSE IF ((r.ge.-2.0).and.(r.le.-1.0)) THEN
        phi_r4 = 1.0/8.0*(5.0+2.0*r-sqrt(-7.0-12.0*r-4.0*r**2))
      ELSE IF ((r.ge.-1.0).and.(r.le.0.0)) THEN
        phi_r4 = 1.0/8.0*(3.0+2.0*r+sqrt(1.0-4.0*r-4.0*r**2))
      ELSE IF ((r.ge.0.0).and.(r.le.1.0)) THEN
        phi_r4 = 1.0/8.0*(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*r**2))
      ELSE IF ((r.ge.1.0).and.(r.le.2.0)) THEN
        phi_r4 = 1.0/8.0*(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*r**2))
      ELSE IF (r.ge.2.0) THEN
        phi_r4 = 0.0
      END IF

      RETURN
      END FUNCTION phi_r4
!#######################################################################
      SUBROUTINE ShapeFunction_MLS(UorV,L,ib)
!-----------------------------------------------------------------------        
!     MOVING LEAST-SQUARE methodology needs to be selected in ibm.for
!     Implemented by Pablo (2015)
!     Algorithm from Ramirez and Nogueira (University of A Coruna)
!     
!     PAPER: Ouro,P., Cea,L., Ramírez,L., & Nogueira,X. (2016). An 
!     immersed boundary method for unstructured meshes in depth averaged
!     shallow water models. International JNMF.
!#######################################################################
      use vars
      use multidata
      use mpi
      use imb
      implicit none

      INTEGER, intent(in) :: UorV,L,ib
      INTEGER :: order_mls,nbase,MAXCELL,MAXNEIG
      INTEGER :: I,J,K,ipface,inodehalo,nn,neignum
      DOUBLE PRECISION :: nmls,k_mls
      DOUBLE PRECISION :: vmaxdist,xN,yN,zN,h0,xv2,yv2,zv2,hx,hy,hz
      DOUBLE PRECISION :: offsety,offsetx,offsetz,ratioHKx1,ratioHKx2
      DOUBLE PRECISION :: ratioHKy1,ratioHKy2,ratioHKz1,ratioHKz2
      DOUBLE PRECISION :: kwx1,kwx2,kwy1,kwy2,kwz1,kwz2,dist_ff
      LOGICAL :: vanella_mls
      INTEGER :: MLS_neignum(maxnodeIBS)
      INTEGER,dimension(9) :: MVefx,MVefy,MVefz
      DOUBLE PRECISION,dimension(9) :: parametros
      DOUBLE PRECISION,dimension(126) :: Xf,Yf,Zf
      DOUBLE PRECISION,dimension(126) :: difx,dify,difz
      DOUBLE PRECISION,dimension(126) :: dX2,dY2,dZ2
      DOUBLE PRECISION,dimension(126) :: W,FF,dFFX,dFFY,dFFZ
      DOUBLE PRECISION,dimension(126) :: ddFFXX,ddFFXY,ddFFYY
      DOUBLE PRECISION,dimension(126) :: ddFFXZ,ddFFZZ,ddFFYZ

!"""  Data for MLS

	MAXCELL=10000 ; MAXNEIG=65 
	MLS_neignum=MAXNEIG
	order_mls=1
	k_mls=5.0d+00

        if (order_mls .eq. 1) nbase=4
        if (order_mls .eq. 2) nbase=10
        if (order_mls .eq. 3) nbase=20
!""""

	Do i=1,MAXNEIG
 	  MVefx(i)=0 ; MVefy(i)=0 ; MVefz(i)=0 
	Enddo

	nn=0  ; nxl=1.9999999

	IF (UorV.eq.1) then
       DO I = 1, dom(ib)%ttc_i 
       IF ( dom(ib)%x(i).gt.(nodex_loc(L)+nxl*dom(ib)%dx) .or.
     &      dom(ib)%x(i).lt.(nodex_loc(L)-nxl*dom(ib)%dx))  GOTO 100
        DO J = 1, dom(ib)%ttc_j
       IF (dom(ib)%yc(j).gt.(nodey_loc(L)+nxl*dom(ib)%dy) .or.
     &     dom(ib)%yc(j).lt.(nodey_loc(L)-nxl*dom(ib)%dy))  GOTO 101
          DO K = 1, dom(ib)%ttc_k
       IF (dom(ib)%zc(k).gt.(nodez_loc(L)+nxl*dom(ib)%dz) .or.
     &     dom(ib)%zc(k).lt.(nodez_loc(L)-nxl*dom(ib)%dz) )  GOTO 102
		nn=nn+1
	IF (nn.gt.MAXNEIG) write(6,*)'More Neighbours than MAXNEIG!',nn
	    MVefx(nn)=I   ;  MVefy(nn)=J	;  MVefz(nn)=K
           Xf(nn)=dom(ib)%x(i)  ; Yf(nn)=dom(ib)%yc(j)
	    Zf(nn)=dom(ib)%zc(k)
102 	CONTINUE
         END DO
101 	CONTINUE
        END DO
100 	CONTINUE
       END DO
	ENDIF

	IF (UorV.eq.2) then
       DO I = 1, dom(ib)%ttc_i
       IF ( dom(ib)%xc(i).gt.(nodex_loc(L)+nxl*dom(ib)%dx) .or.
     &      dom(ib)%xc(i).lt.(nodex_loc(L)-nxl*dom(ib)%dx) )  GOTO 200 
        DO J = 1, dom(ib)%ttc_j
       IF (dom(ib)%y(j) .gt.(nodey_loc(L)+nxl*dom(ib)%dy) .or.
     &     dom(ib)%y(j) .lt.(nodey_loc(L)-nxl*dom(ib)%dy) )  GOTO 201
          DO K = 1, dom(ib)%ttc_k
	IF( dom(ib)%zc(k).gt.(nodez_loc(L)+nxl*dom(ib)%dz) .or.
     &     dom(ib)%zc(k).lt.(nodez_loc(L)-nxl*dom(ib)%dz) )  GOTO 202
		nn=nn+1
	IF (nn.gt.MAXNEIG)
     &       write(6,*)'More Neighbours than MAXNEIG!',nn,nxl
		MVefx(nn)=I   ;  MVefy(nn)=J	;  MVefz(nn)=K
              Xf(nn)=dom(ib)%xc(i)  
	       Yf(nn)=dom(ib)%y(j)
	       Zf(nn)=dom(ib)%zc(k)
202 	CONTINUE
         END DO
201 	CONTINUE
        END DO
200 	CONTINUE
       END DO
	ENDIF

	IF (UorV.eq.3) then
       DO I = 1, dom(ib)%ttc_i 
       IF ( dom(ib)%xc(i).gt.(nodex_loc(L)+nxl*dom(ib)%dx) .or.
     &      dom(ib)%xc(i).lt.(nodex_loc(L)-nxl*dom(ib)%dx)) GOTO 300
        DO J = 1, dom(ib)%ttc_j
       IF ( dom(ib)%yc(j).gt.(nodey_loc(L)+nxl*dom(ib)%dy) .or.
     &      dom(ib)%yc(j).lt.(nodey_loc(L)-nxl*dom(ib)%dy)) GOTO 301
          DO K = 1, dom(ib)%ttc_k
       IF ( dom(ib)%z(k) .gt.(nodez_loc(L)+nxl*dom(ib)%dz) .or.
     &      dom(ib)%z(k) .lt.(nodez_loc(L)-nxl*dom(ib)%dz) ) GOTO 302
		nn=nn+1
	IF (nn.gt.MAXNEIG)
     &       write(6,*)'More Neighbours than MAXNEIG!',nn,nxl
		MVefx(nn)=I   ;  MVefy(nn)=J	;  MVefz(nn)=K
              Xf(nn)=dom(ib)%xc(i)  
	       Yf(nn)=dom(ib)%yc(j)
	       Zf(nn)=dom(ib)%z(k)
302 	CONTINUE
         END DO
301 	CONTINUE
        END DO
300 	CONTINUE
       END DO
	ENDIF

	neignum=nn
	MLS_neignum(L)=neignum	
!--------------------------------------------------------
        vmaxdist= 0.d+00   
    	 do j=1, neignum
	  difx(j)= 0.d+0 ;    dify(j)= 0.d+0;    difz(j)= 0.d+00
	 enddo

    	  do j=1, neignum
	 dist_ff=DSQRT((nodex_loc(L)-Xf(j))**2+(nodey_loc(L)-Yf(j))**2
     &		+(nodez_loc(L)-Zf(j))**2)
	    if(dist_ff.gt.vmaxdist) vmaxdist= dist_ff
	     	difx(j)=Xf(j)-nodex_loc(L)
		dify(j)=Yf(j)-nodey_loc(L)
		difz(j)=Zf(j)-nodez_loc(L)
	  enddo

	h0=1.0d+00
	hx= h0*maxval(DABS(difx(1:neignum)))
	hy= h0*maxval(DABS(dify(1:neignum)))
	hz= h0*maxval(DABS(difz(1:neignum)))

	parametros(1)= k_mls
	parametros(2)= k_mls
	parametros(3)= k_mls
	parametros(4)= k_mls
	parametros(5)= k_mls
	parametros(6)= k_mls
	parametros(7)= 0.0d+00
	parametros(8)= 0.0d+00
	parametros(9)= 0.0d+00

	ratioHKx1= parametros(1)
	ratioHKx2= parametros(2)
	ratioHKy1= parametros(3)
	ratioHKy2= parametros(4)
	ratioHKz1= parametros(5)
	ratioHKz2= parametros(6)
	offsetx  = parametros(7)
	offsety  = parametros(8)
	offsetz  = parametros(9)

	h0= 1.0d+00
	kwx1= h0*ratioHKx1
	kwx2= h0*ratioHKx2
	kwy1= h0*ratioHKy1
	kwy2= h0*ratioHKy2
	kwz1= h0*ratioHKz1
	kwz2= h0*ratioHKz2

	call  kerneltododer(W,dX2,dY2,dZ2,difx,dify,difz,hx,hy,hz,
     &  neignum,kwx1,kwx2,kwy1,kwy2,kwz1,kwz2,MAXNEIG)

        call forma3Dder(FF,W,dX2,dY2,dZ2,neignum,Xf,Yf,Zf
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),hx,hy,hz,
     & dFFX,dFFY,dFFZ,ddFFXX,ddFFXY,ddFFXZ,ddFFYY,ddFFYZ,ddFFZZ,
     &  MAXNEIG,nbase)

	  do j=1, MLS_neignum(L)
	   IF (UorV.eq.1) then
	   	dh1_loc(L,j)=FF(j)  	;I_nr_U(L,j)=MVefx(j)
	    	J_nr_U(L,j)=MVefy(j)	;K_nr_U(L,j)=MVefz(j)
	   ELSEIF(UorV.eq.2) then
	    	dh2_loc(L,j)=FF(j)	;I_nr_V(L,j)=MVefx(j)
		J_nr_V(L,j)=MVefy(j)	;K_nr_V(L,j)=MVefz(j)
	   ELSEIF(UorV.eq.3) then
	   	dh3_loc(L,j)=FF(j)	;I_nr_W(L,j)=MVefx(j)
		J_nr_W(L,j)=MVefy(j)	;K_nr_W(L,j)=MVefz(j)
	   ENDIF
          enddo

	  IF(UorV.eq.1) KmaxU(L)=MLS_neignum(L)
	  IF(UorV.eq.2) KmaxV(L)=MLS_neignum(L)
	  IF(UorV.eq.3) KmaxW(L)=MLS_neignum(L)

      RETURN
      END SUBROUTINE ShapeFunction_MLS
!#######################################################################
      SUBROUTINE kerneltododer(W,dX1,dY1,dZ1,difx,dify,difz,hx,hy,hz,
     &  neignum,kwx1,kwx2,kwy1,kwy2,kwz1,kwz2,MAXNEIG)
!-----------------------------------------------------------------------
!     This kernel function provides a neightbour weight distribution for
!     the calculation of the MLS shape functions vector.
!#######################################################################
      use mpi
      implicit none

      INTEGER :: I,J
      INTEGER, intent(in) :: neignum,MAXNEIG
      DOUBLE PRECISION,intent(in) :: hx,hy,hz,kwx1,kwx2,kwy1
      DOUBLE PRECISION,intent(in) :: kwy2,kwz1,kwz2
      DOUBLE PRECISION,intent(in) :: difx(MAXNEIG)
      DOUBLE PRECISION,intent(in) :: dify(MAXNEIG)
      DOUBLE PRECISION,intent(in) :: difz(MAXNEIG)
      DOUBLE PRECISION :: cero, coef,cx1,cx2,cy1,cy2,cz1,cz2,PI
      DOUBLE PRECISION :: dmx,dmy,dmz,Wx,Wy,Wz,dWx,dWy,dWz
      DOUBLE PRECISION :: coef1x,coef2x,coefd1x,coefd2x,coefd3x
      DOUBLE PRECISION :: coef1y,coef2y,coefd1y,coefd2y,coefd3y
      DOUBLE PRECISION :: coef1z,coef2z,coefd1z,coefd2z,coefd3z
      DOUBLE PRECISION,dimension(MAXNEIG) :: sx,sy,sz
      DOUBLE PRECISION,dimension(MAXNEIG) :: W,dX1,dY1,dZ1

      cero=1.d-08  ;   coef=1.d+00  ; PI = 4.D0*DATAN(1.D0)

      do j=1, neignum
	 sx(j)=0.0d+00   ;  sy(j)=0.0d+00;  sz(j)=0.0d+00
      enddo

! First, we eveluate the smoothing length:
      do j=1, neignum
	sx(j)=dabs(difx(j)) ; sy(j)=dabs(dify(j)); sz(j)=dabs(difz(j))
	dmx=2.0d+00*hx 	; dmy=2.0d+00*hy 	;  dmz=2.0d+00*hz
	cx1= dmx/kwx1  	;  cx2= dmx/kwx2	
	cy1= dmy/kwy1	;  cy2= dmy/kwy2	
	cz1= dmz/kwz1	;  cz2= dmz/kwz2	

	if (difx(j).ge.0.0d+00)then
	  coef1x= exp(-(sx(j)/cx1)**2) - exp(-(dmx/cx1)**2)
	  coef2x= 1.0d+00 - exp(-(dmx/cx1)**2)
	  coefd1x=exp(-(sx(j)/cx1)**2)
	  coefd2x=coef2x
	  coefd3x=-2.0d+00*sx(j)/(cx1*cx1)
	else
	  coef1x= exp(-(sx(j)/cx2)**2) - exp(-(dmx/cx2)**2)
	  coef2x= 1.0d+00 - exp(-(dmx/cx2)**2)
	  coefd1x=exp(-(sx(j)/cx2)**2)
	  coefd2x=coef2x
	  coefd3x=-2.0d+00*sx(j)/(cx2*cx2)
	endif
	if (dify(j).ge.0.0d+00)then
	  coef1y= exp(-(sy(j)/cy1)**2) - exp(-(dmy/cy1)**2)
	  coef2y= 1.0d+00 - exp(-(dmy/cy1)**2)	
	  coefd1y=exp(-(sy(j)/cy1)**2)
	  coefd2y=coef2y
	  coefd3y=-2.0d+00*sy(j)/(cy1*cy1)
	else
	  coef1y= exp(-(sy(j)/cy2)**2) - exp(-(dmy/cy2)**2)
	  coef2y= 1.0d+00 - exp(-(dmy/cy2)**2)
	  coefd1y=exp(-(sy(j)/cy2)**2)
	  coefd2y=coef2y
	  coefd3y=-2.0d+00*sy(j)/(cy2*cy2)
	endif
	if (difz(j).ge.0.0d+00)then
	  coef1z= exp(-(sz(j)/cz1)**2) - exp(-(dmz/cz1)**2)
	  coef2z= 1.0d+00 - exp(-(dmz/cz1)**2)	
	  coefd1z=exp(-(sz(j)/cz1)**2)
	  coefd2z=coef2z
	  coefd3z=-2.0d+00*sz(j)/(cz1*cz1)
	else
	  coef1z= exp(-(sz(j)/cz2)**2) - exp(-(dmz/cz2)**2)
	  coef2z= 1.0d+00 - exp(-(dmz/cz2)**2)
	  coefd1z=exp(-(sz(j)/cz2)**2)
	  coefd2z=coef2z
	  coefd3z=-2.0d+00*sz(j)/(cz2*cz2)
	endif

	Wx= coef1x/coef2x ; Wy= coef1y/coef2y ;	Wz= coef1z/coef2z

	dWx= coefd3x*coefd1x/coefd2x
	dWy= coefd3y*coefd1y/coefd2y
	dWz= coefd3z*coefd1z/coefd2z

	W(j)= Wx*Wy*Wz 

	if(sx(j).gt.cero)then
	  dX1(j)= Wy*Wz*dWx*(difx(j)/sx(j))
	 else
	  dX1(j)= Wz*Wy*dWx*(difx(j))  
	endif
	if(sy(j).gt.cero)then
	  dY1(j)= Wx*Wz*dWy*(dify(j)/sy(j))
	 else
	  dY1(j)= Wx*Wz*dWy*(dify(j))
	endif
	if(sz(j).gt.cero)then
	  dZ1(j)= Wx*Wy*dWz*(difz(j)/sz(j))
	 else
	  dZ1(j)= Wx*Wy*dWz*(difz(j)) 
	endif
      enddo !j-loop

      RETURN
      END SUBROUTINE kerneltododer
!#######################################################################
      SUBROUTINE forma3Dder(FF,W2,dX2,dY2,dZ2,neignum,Xmls,Ymls,Zmls
     &  ,xgau,ygau,zgau,hx,hy,hz,
     &   dFFX,dFFY,dFFZ,ddFFXX,ddFFXY,ddFFXZ,ddFFYY,ddFFYZ,ddFFZZ,
     &   MAXNEIG,nbase)
!-----------------------------------------------------------------------
!     Solves the moment matrix for the MLS interpolation functions.
!#######################################################################
      implicit none
      
      INTEGER, intent(in) :: neignum,MAXNEIG,nbase
      DOUBLE PRECISION,intent(in) :: hx,hy,hz,xgau,ygau,zgau
      DOUBLE PRECISION,intent(in) :: W2(MAXNEIG)
      DOUBLE PRECISION,intent(in) :: Xmls(MAXNEIG),Ymls(MAXNEIG)
      DOUBLE PRECISION,intent(in) :: Zmls(MAXNEIG)
      DOUBLE PRECISION,intent(in) :: dX2(MAXNEIG),dY2(MAXNEIG)
      DOUBLE PRECISION,intent(in) :: dZ2(MAXNEIG)
      INTEGER :: I,J,K,ifila,icolumna,ibase,nelim,icol
      DOUBLE PRECISION :: xg,yg,zg,wmax,wmin,valor
      DOUBLE PRECISION,dimension(MAXNEIG) :: FF,dFFX,dFFY,dFFZ
      DOUBLE PRECISION,dimension(MAXNEIG) :: ddFFXX,ddFFXY,ddFFYY
      DOUBLE PRECISION,dimension(MAXNEIG) :: ddFFXZ,ddFFYZ,ddFFZZ
      DOUBLE PRECISION :: A(nbase,nbase),vIA(nbase,nbase)
      DOUBLE PRECISION :: eye(nbase,nbase)
      DOUBLE PRECISION :: C(nbase,MAXNEIG),PtC(MAXNEIG,MAXNEIG)
      DOUBLE PRECISION :: auxiliar(MAXNEIG,nbase)
      DOUBLE PRECISION :: auxiliar1(nbase),auxiliar2(nbase)
      DOUBLE PRECISION :: sol(nbase,nbase),wsv(10),vsv(10,10),coef

! Inicialisation of the matrix C
        do i=1,nbase
	 do j=1, MAXNEIG
	  C(i,j)=0.d+00
	 enddo
	enddo

! Creation of the matrix A
      do i=1,nbase
	do j=1,nbase
	 A(i,j)= 0.d+00
	 vIA(i,j)=0.d+00
	 if(i.eq.j)then
	  eye(i,j)= 1.d+00
	 else
	  eye(i,j)= 0.d+00
	 endif
	enddo
      enddo

! Creation of a auxilary vector to do the calculs:
      do k=1,neignum
       xg=(Xmls(k)-xgau)/hx ; yg=(Ymls(k)-ygau)/hy ;zg=(Zmls(k)-zgau)/hz

      if(nbase.eq.4) then
      	auxiliar(k, 1)= 1.d+00
	auxiliar(k, 2)= xg
	auxiliar(k, 3)= yg
	auxiliar(k, 4)= zg
      elseif(nbase.eq.10) then
	auxiliar(k, 1)= 1.d+00
	auxiliar(k, 2)= xg
	auxiliar(k, 3)= yg
	auxiliar(k, 4)= zg
	auxiliar(k, 5)= xg*xg
	auxiliar(k, 6)= xg*yg
	auxiliar(k, 7)= xg*zg
	auxiliar(k, 8)= yg*yg
	auxiliar(k, 9)= yg*zg
	auxiliar(k,10)= zg*zg
      elseif(nbase.eq.20) then
	auxiliar(k, 1)= 1.d+00
	auxiliar(k, 2)= xg
	auxiliar(k, 3)= yg
	auxiliar(k, 4)= zg
	auxiliar(k, 5)= xg*xg
	auxiliar(k, 6)= xg*yg
	auxiliar(k, 7)= xg*zg
	auxiliar(k, 8)= yg*yg
	auxiliar(k, 9)= yg*zg
	auxiliar(k,10)= zg*zg
	auxiliar(k,11)= xg*xg*xg
	auxiliar(k,12)= xg*xg*yg
	auxiliar(k,13)= xg*xg*zg
	auxiliar(k,14)= xg*yg*zg
	auxiliar(k,15)= yg*yg*yg
	auxiliar(k,16)= yg*yg*xg
	auxiliar(k,17)= yg*yg*zg
	auxiliar(k,18)= zg*zg*zg
	auxiliar(k,19)= zg*zg*xg
	auxiliar(k,20)= zg*zg*yg
       endif

       do ifila= 1,nbase  
	do icolumna= 1,nbase
     	 A(ifila,icolumna)= A(ifila,icolumna)+
     &         W2(k)*auxiliar(k,ifila)*auxiliar(k,icolumna)
	enddo
       enddo
      enddo

! Calculation of the inverse matrix of A
! For this calculation a system of 6 equations is resolved
! First, triangulation of the matrix A and perform the operation per
! rows upon the matrix identity.
 
      do i=1,nbase-1
! Check that the pivo is not null
        if(A(i,i).lt.1.d-06)then
	  write(6,*)' Pivote nulo',hx,hy,hz,i
	  pause
	endif
       do j=i+1,nbase
        coef=-A(j,i)/A(i,i)
          do k=i,nbase
            A(j,k)=A(j,k)+coef*A(i,k)
          enddo
! Apply the operacion on the identity columns
	    do icol=1,nbase
            eye(j,icol)= eye(j,icol)+coef*eye(i,icol)
	    enddo
        enddo
      enddo
      do icol=1,nbase
	do i=1,nbase
	 auxiliar1(i)= eye(i,icol)
	enddo
         call Back_Substitution(nbase,A,auxiliar1,auxiliar2)
	do i=1,nbase
	 vIA(i,icol)= auxiliar2(i)
	enddo
       enddo


! We calculate C, which is A^-1*B
      do i=1,nbase
	do k=1, neignum
 	 C(i,k)= 0.d+00
	  do j=1,nbase  
	  C(i,k)= C(i,k) + vIA(i,j)*auxiliar(k,j)
	  enddo
	 C(i,k)= C(i,k)*W2(k)
	enddo
      enddo

! We calculate I-Pt*C
      do i=1,neignum
	do j=1,neignum

	valor= 0.d+00
	do ibase=1,nbase
		valor= valor + C(ibase,j)*auxiliar(i,ibase)
	enddo

		PtC(i,j)=-valor
		if(i.eq.j)PtC(i,j)=PtC(i,j)+1.d+00
	  enddo
	enddo

      if(nbase.eq.4) then
        do k=1,neignum
	   FF(k)= C(1,k)
	enddo
        do j=1,neignum
	 dFFX(j)=C(2,j)/hx
	 dFFY(j)=C(3,j)/hy
	 dFFZ(j)=C(4,j)/hz
	 do k=1,neignum
           dFFX(j)=dFFX(j)+dX2(k)*C(1,k)*PtC(k,j)/W2(k)
           dFFY(j)=dFFY(j)+dY2(k)*C(1,k)*PtC(k,j)/W2(k)
           dFFZ(j)=dFFZ(j)+dZ2(k)*C(1,k)*PtC(k,j)/W2(k)
	 enddo
	enddo
      elseif(nbase.eq.10) then
        do k=1,neignum
 	  FF(k)= C(1,k)
	enddo
        do j=1,neignum
	 dFFX(j)=C(2,j)/hx
	 dFFY(j)=C(3,j)/hy
	 dFFZ(j)=C(4,j)/hz
            do k=1,neignum
             dFFX(j)=dFFX(j)+dX2(k)*C(1,k)*PtC(k,j)/W2(k)
             dFFY(j)=dFFY(j)+dY2(k)*C(1,k)*PtC(k,j)/W2(k)
             dFFZ(j)=dFFZ(j)+dZ2(k)*C(1,k)*PtC(k,j)/W2(k)
            enddo
	enddo
        do j=1,neignum
	  ddFFXY(j)=        C(5,j)/(hx*hy)
	  ddFFXZ(j)=        C(6,j)/(hx*hz)
	  ddFFYZ(j)=        C(7,j)/(hy*hz)
	  ddFFXX(j)= 2.d+00*C(8,j)/(hx*hx)  
	  ddFFYY(j)= 2.d+00*C(9,j)/(hy*hy)
	  ddFFZZ(j)=2.d+00*C(10,j)/(hz*hz)
        enddo
      elseif(nbase.eq.20) then
        do k=1,neignum
	  FF  (k)= C(1,k)
	enddo
        do j=1,neignum
	  dFFX(j)=C(2,j)/hx
	  dFFY(j)=C(3,j)/hy
	  dFFZ(j)=C(4,j)/hZ
	 do k=1,neignum
           dFFX(j)=dFFX(j)+dX2(k)*C(1,k)*PtC(k,j)/W2(k)
           dFFY(j)=dFFY(j)+dY2(k)*C(1,k)*PtC(k,j)/W2(k)
           dFFZ(j)=dFFZ(j)+dZ2(k)*C(1,k)*PtC(k,j)/W2(k)
	 enddo
	enddo
         do j=1,neignum
	  ddFFXX(j)= 2.d+00*C( 5,j)/(hx*hx)
	  ddFFXY(j)=        C( 6,j)/(hx*hy)
	  ddFFXZ(j)=        C( 7,j)/(hx*hz)
	  ddFFYY(j)= 2.d+00*C( 8,j)/(hy*hy)
	  ddFFYZ(j)=        C( 9,j)/(hy*hz)
	  ddFFZZ(j)= 2.d+00*C(10,j)/(hz*hz)
         enddo
	endif

	RETURN
	END SUBROUTINE forma3Dder
!#######################################################################
      SUBROUTINE Back_Substitution(ndim,A,b,x)
!#######################################################################
      implicit none

      INTEGER :: i,j
      INTEGER,intent(in) :: ndim
      DOUBLE PRECISION :: suma
      DOUBLE PRECISION :: A(ndim,ndim),b(ndim),x(ndim)

      x(ndim)=b(ndim)/A(ndim,ndim)
      do i=ndim-1,1,-1
        suma=0
        do j=i+1,ndim
          suma=suma+A(i,j)*x(j)
        enddo
        x(i)=(b(i)-suma)/A(i,i)
      enddo
      
      RETURN
      END SUBROUTINE Back_Substitution

