# ShapeFunction\_MLS

**File:** DeltaF\_MLS.for\
\
**Subroutine Called in:** \
deltah (ibm.for) -> IB_previous (ibm.for) -> flosol (flosol.for)_\
_imb\__openmpi (ibm.for) -> IBM (ibm.for) -> flosol (flosol.for)\
\
**Purpose:**\
****\
****\
******User Likeness to alter the file:** \
**Rarely**&#x20;

```
!######################################################################
	subroutine ShapeFunction_MLS(UorV,L,ib)
!######################################################################
         use vars
	 use multidata
         use mpi
	 use imb
         implicit none
         INTEGER :: order_mls,nbase,MAXCELL,MAXNEIG
	 LOGICAL :: vanella_mls
         INTEGER, intent(in) :: UorV,L,ib
         INTEGER :: I,J,K,ipface,inodehalo,nn,neignum
         INTEGER:: MVefx(126),MVefy(126),MVefz(126)
         Double precision :: vmaxdist,xN,yN,zN,h0,xv2,yv2,zv2
    	 Double precision::  Xf(126),Yf(126),Zf(126),hx,hy,hz
         Double precision :: offsety,offsetx,offsetz,ratioHKx1,ratioHKx2
         Double precision :: ratioHKy1,ratioHKy2,ratioHKz1,ratioHKz2
         Double precision :: kwx1,kwx2,kwy1,kwy2,kwz1,kwz2,dist_ff
	 Double precision :: parametros(9),W(126),difx(126)
	 Double precision :: dify(126),difz(126)
         Double precision :: dX2(126),dY2(126),dZ2(126)
         Double precision :: FF(126),dFFX(126),dFFY(126),dFFZ(126)
         Double precision :: ddFFXX(126),ddFFXY(126),ddFFYY(126)
         Double precision :: ddFFXZ(126),ddFFZZ(126),ddFFYZ(126)
	 Double precision :: nmls,k_mls
         INTEGER :: MLS_neignum(maxnodeIBS)

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

	return
	end
```
