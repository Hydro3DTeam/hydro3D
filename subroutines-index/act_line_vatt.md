# Act\_line\_VATT

File: actuator.for

```
!######################################################################
      SUBROUTINE Act_line_VATT(M,L,ib)
!######################################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: M,L,ib
	integer :: fileact
      DOUBLE PRECISION :: Utheta,Vrel,Vrot,phi,alpeff,Lift,Drag,c_in
	DOUBLE PRECISION :: betaeff,Vperp,Vtan,radlocal,alphalocal
!actuator line - VATT: forced velocities based on rotational speed
		c_in=0.14d0 !Change chord length accordingly
!Local velocity at the marker: Vlocal=Vperp,Vtan
		Vtan = U_Beta1_loc(L)*DCOS(rads(M)+alpha0_loc(L))
     &  		+U_Beta2_loc(L)*DSIN(rads(M)+alpha0_loc(L))
		Vperp= U_Beta1_loc(L)*DSIN(rads(M)+alpha0_loc(L))
     &  		-U_Beta2_loc(L)*DCOS(rads(M)+alpha0_loc(L))
!Effective angle of attack	
	alpeff=ATAN(Vperp/(Vtan+radsin(M)*R0_loc(L)))
!Module of the effective velocity
	Vrel=DSQRT(Vperp**2+(Vtan+radsin(M)*R0_loc(L))**2)

	radlocal = INT((rads(M)+alpha0_loc(L))/(2.d0*3.1416))       !Number of revs performed
	radlocal = rads(M)+alpha0_loc(L)-radlocal*2.d0*3.1416d0	!Angle rotated in current rev
!Angle of attack
	radlocal = radlocal *180.d0/3.1416d0				!Rads to degrees
!LIFT COEFFICIENT
	Cl_act= -0.470319202693645d0
     & -(9.83109403445012/10D0**30)*radlocal**14
     & +(2.35506917977834/10D0**26)*radlocal**13
     & -(2.49085514318290/10D0**23)*radlocal**12
     & +(1.52861741834597/10D0**20)*radlocal**11
     & -(6.01348602214706/10D0**18)*radlocal**10
     & +(1.57992317569433/10D0**15)*radlocal**9
     & -(2.79676685658241/10D0**13)*radlocal**8
     & +(3.27393735916174/10D0**11)*radlocal**7
     & -(2.39204712729446/10D0**9)*radlocal**6
     & +(9.28650148346635/10D0**8)*radlocal**5
     & -(5.91984268776261/10D0**7)*radlocal**4
     & -(9.47384394282790/10D0**5)*radlocal**3
     & +(3.14028704853696/10D0**3)*radlocal**2
     & +(2.19265827600333/10D0**2)*radlocal
	
	if(radlocal.gt.300 .and. Cl_act.lt.-0.5) Cl_act=-0.5

!DRAG COEFFICIENT
	if(radlocal.le.35.d0) then
	Cd_act = -0.140115557534005d0 
     & +(2.96862876671061/10D0**4)*radlocal**2
     & -(6.39641994261401/10D0**3)*radlocal
	elseif(radlocal.ge.352.d0) then
	Cd_act = 1.61032763683459 
     & -(4.87756197364469/10D0**3)*radlocal
	else
	Cd_act= -0.194577074144036d0
     & -(7.81195089028066/10D0**34)*radlocal**16
     & +(2.28728899458862/10D0**30)*radlocal**15
     & -(3.03482604943372/10D0**27)*radlocal**14
     & +(2.41321554737758/10D0**24)*radlocal**13
     & -(1.28157461368188/10D0**21)*radlocal**12
     & +(4.79472187285026/10D0**19)*radlocal**11
     & -(1.29899815564960/10D0**16)*radlocal**10
     & +(2.57937811840519/10D0**14)*radlocal**9
     & -(3.75662265252506/10D0**12)*radlocal**8
     & +(3.97565875154755/10D0**10)*radlocal**7
     & -(2.99730852027196/10D0**8)*radlocal**6
     & +(1.55841955992838/10D0**6)*radlocal**5
     & -(5.31865863851434/10D0**5)*radlocal**4
     & +(1.10154286382894/10D0**3)*radlocal**3
     & -(1.18896450651628/10D0**2)*radlocal**2
     & +(5.01729242562662/10D0**2)*radlocal
	end if

	radlocal = radlocal /(180.d0/3.1416d0)			!Deg to rads
!Compute equivalent hydrodynamic forces:
      Lift=0.5d0*Cl_act*c_in*dom(ib)%dz*(Vrel**2)
     &	   /(radsin(M)*R0_loc(L))**2  !   /2.5
      Drag=0.5d0*Cd_act*c_in*dom(ib)%dz*(Vrel**2)
     &	   /(radsin(M)*R0_loc(L))**2  !*Vrel**2   /2.5	
!Corrections ??
!Compute the projection of the forces over X-Y-Z coordinates
	betaeff = radlocal - alpeff
	FX1_loc(L)=-Drag*DCOS(betaeff)-Lift*DSIN(betaeff)
	FX2_loc(L)= Drag*DSIN(betaeff)-Lift*DCOS(betaeff)
	FX3_loc(L)= 0.d0

!Normalise by the area so the ibm subroutine can be used (Wu and Porte-Agel '11): - CHECK THIS STEP - 
	FX1_loc(L)= FX1_loc(L)/(dom(ib)%dz*c_in)
	FX2_loc(L)= FX2_loc(L)/(dom(ib)%dz*c_in)

   88 FORMAT (1f12.6,5e15.7)
      RETURN
      END SUBROUTINE 
```
