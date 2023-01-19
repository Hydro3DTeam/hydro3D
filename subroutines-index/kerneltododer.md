# kerneltododer

**File:** DeltaF\__MLS.for_\
__\
__**Subroutine Called in:** __ \
_ShapeFunction (DeltaF_\_MLS.for) -> deltah (ibm.for) -> IB_previous (ibm.for) -> flosol (flosol.for)_\
_ShapeFunction (DeltaF_\_MLS.for_) -> imb\__openmpi (ibm.for) -> IBM (ibm.for) -> flosol (flosol.for)\
\
**Purpose:**\
****\
****\
******User Likeness to alter the file:** \
**Rarely**&#x20;

```
!######################################################################
      subroutine kerneltododer(W,dX1,dY1,dZ1,difx,dify,difz,hx,hy,hz,
     &  neignum,kwx1,kwx2,kwy1,kwy2,kwz1,kwz2,MAXNEIG)
!######################################################################
      use mpi
      implicit none
      Double precision, intent(in) :: hx,hy,hz,kwx1,kwx2,kwy1
      Double precision, intent(in) :: kwy2,kwz1,kwz2	
      Double precision, intent(in) :: difx(MAXNEIG)
      Double precision, intent(in) ::dify(MAXNEIG)
      Double precision, intent(in) :: difz(MAXNEIG)
      INTEGER, intent(in):: neignum,MAXNEIG
      INTEGER :: I,J
      Double precision:: cero, coef,cx1,cx2,cy1,cy2,cz1,cz2,PI
      Double precision:: dmx,dmy,dmz,Wx,Wy,Wz,dWx,dWy,dWz
      Double precision:: coef1x,coef2x,coefd1x,coefd2x,coefd3x
      Double precision:: coef1y,coef2y,coefd1y,coefd2y,coefd3y
      Double precision:: coef1z,coef2z,coefd1z,coefd2z,coefd3z
      Double precision:: sx(MAXNEIG),sy(MAXNEIG),sz(MAXNEIG)
      Double precision:: W(MAXNEIG)
      Double precision:: dX1(MAXNEIG),dY1(MAXNEIG),dZ1(MAXNEIG)

        cero=1.d-08  ;   coef=1.d+00  ;    PI = 4.D0*DATAN(1.D0)

        do j=1, neignum
	  sx(j)=0.0d+00   ;  sy(j)=0.0d+00;  sz(j)=0.0d+00
	enddo

! Primero cogemos las distancias que nos interesan
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

	return
	end
```
