# ActuatorLine

File: actuator.for

```
!######################################################################
      SUBROUTINE ActuatorLine(M,L,ib)
!######################################################################
      use vars
      use multidata 
      use imb
      use mpi
      implicit none
      integer, intent (in) :: M,L,ib
      DOUBLE PRECISION::Utheta,Vrel,Vrot,phi,alpeff,Lift,Drag
			DOUBLE PRECISION::F1,GVAL,expval
			DOUBLE PRECISION,PARAMETER:: PI = 4.D0*DATAN(1.D0)
			real :: Radius,NumBlades
			Radius = 0.135
			NumBlades = 3.0
			
!actuator line - HATT: forced velocities based on rotational speed
!Relative velocity in the plane of rotation, i.e. y-z for HATT
	   Utheta= radsin(M)*R0_loc(L) - (U_Beta2_loc(L)*DCOS(rads(M)
     &    +alpha0_loc(L)) - U_Beta3_loc(L)*DSIN(rads(M)+alpha0_loc(L)))  

!Module of the effective velocity
		Vrel=DSQRT(U_Beta1_loc(L)**2+Utheta**2)

!Effective angle of attack	
		alpeff= (ATAN((U_Beta1_loc(L))/Utheta+1.D-10))

!ABS value is used because with negative alpeff the correction factor gives error
!Angle of attack, Pit_act= local pitch angle. Both in DEGREES
		phi = alpeff*180.D0/3.1416D0 - Pit_act(L)
		
!Calculate the lift and drag coefficients depending on the angle and Re
		call Stallard(phi)
		
!Compute equivalent hydrodynamic forces:
		Lift=0.5d0*1.d0*Cl_act*c_act(L)*dom(ib)%dz*Vrel**2
		Drag=0.5d0*1.d0*Cd_act*c_act(L)*dom(ib)%dz*Vrel**2		
		
!Corrections from Shen et al.
	   gval=dexp(-0.125d0*(NumBlades*radsin(M)*Radius/ubulk-21.d0))
     &   +0.1d0
	   expval=(-gval*(NumBlades*(Radius-r_act(L))+1.d-12)
     &       /(2.d0*r_act(L)*DABS(dsin(alpeff))))
		F1=2.d0/PI*DACOS(dexp(expval))
		Lift=Lift*F1 ; Drag=Drag*F1
		
!Compute the projection of the forces over X-Y-Z coordinates
		FX1_loc(L)= ( Drag*dsin(alpeff)+Lift*dcos(alpeff))
		FX2_loc(L)= (-Drag*dcos(alpeff)+Lift*dsin(alpeff))
     &  *Dcos(rads(M)+alpha0_loc(L))
		FX3_loc(L)=-(-Drag*dcos(alpeff)+Lift*dsin(alpeff))
     &  *Dsin(rads(M)+alpha0_loc(L))
     
!Normalise by the area so the ibm subroutine can be used (Wu and Porte-Agel '11): - CHECK THIS STEP - 
		FX1_loc(L)=-FX1_loc(L)/(dom(ib)%dz*c_act(L))
		FX2_loc(L)=-FX2_loc(L)/(dom(ib)%dz*c_act(L))
		FX3_loc(L)=-FX3_loc(L)/(dom(ib)%dz*c_act(L))

   88 FORMAT (1f16.8,5e15.7)
      RETURN
      END SUBROUTINE 
```
