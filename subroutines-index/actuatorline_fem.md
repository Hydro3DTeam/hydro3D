# ActuatorLine\_FEM

File: actuator.for

```
######################################################################
      SUBROUTINE ActuatorLine_FEM(numIB)
!######################################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: numIB
      integer ::M,L,fileact,fileact2
      DOUBLE PRECISION :: PI,alpharads,RtoHub
      DOUBLE PRECISION :: T1_ACT,T2_ACT,T3_ACT,Q1_ACT,Q2_ACT,Q3_ACT
      DOUBLE PRECISION :: BR1_ACT,BR2_ACT,BR3_ACT
      DOUBLE PRECISION :: BT1_ACT,BT2_ACT,BT3_ACT
      DOUBLE PRECISION :: Myaw,MtwrX,MtwrY,rho,FQ
      CHARACTER*8  :: char_block
      CHARACTER*31 :: gridfile
      INTEGER:: strlen

 	IF(itime.eq. itime_start) then
	  F_X_MEAN=0.D0  ;F_Q_MEAN=0.D0 ; F_Y_MEAN=0.D0 ; F_Z_MEAN=0.D0 
	ENDIF

	PI = 4.D0*DATAN(1.D0) ; M=numIB
	alpharads=rads(M)*180.d0/PI     ; rho=1000.d0		
	T1_ACT=0.D0 	; T2_ACT=0.D0 	; T3_ACT=0.D0 
	Q1_ACT=0.D0 	; Q2_ACT=0.D0 	; Q3_ACT=0.D0 
	BR1_ACT=0.D0	; BR2_ACT=0.D0 	; BR3_ACT=0.D0 
	BT1_ACT=0.D0 	; BT2_ACT=0.D0 	; BT3_ACT=0.D0 
	Myaw=0.D0 	; MtwrY=0.D0 	; MtwrX=0.D0 

! Time-averaged coefficients

!Open file to output mean coefficients:
	   write(char_block,'(i3)') numIB
         strlen=LEN(TRIM(ADJUSTL(char_block)))
         char_block=REPEAT('0',(3-strlen))//TRIM(ADJUSTL(char_block))
         gridfile='MeanCoefs_ActL_'//TRIM(ADJUSTL(char_block))//'.dat'
	   fileact=1400+numIB
         open (unit=fileact, file=gridfile)
	   write(fileact,*)'Variables="r/R","Fx","Fq","Fy","Fz"'
	 do L=1,nodes(M)
	   if(L.le.nodes(M)/3) THEN

	 F_X_MEAN(M,L)=(F_X_MEAN(M,L)*(ITIME-1)+FX1(M,L))/ITIME
	 F_Y_MEAN(M,L)=(F_Y_MEAN(M,L)*(ITIME-1)+FX2(M,L))/ITIME
	 F_Z_MEAN(M,L)=(F_Z_MEAN(M,L)*(ITIME-1)+FX3(M,L))/ITIME
	 FQ=(FX2(M,L)*dcos(rads(M)+alpha0(M,L))
     &    -FX3(M,L)*dsin(rads(M)+alpha0(M,L)))*R0(M,L)*radsin(M)
	 F_Q_MEAN(M,L)=(F_Q_MEAN(M,L)*(ITIME-1)+FQ)/ITIME

	write(fileact,'(10E18.9)')r_act(L)/r_act(nodes(M)/3)
     &    ,F_X_MEAN(M,L)*1000.d0,F_Q_MEAN(M,L)*1000.d0
     &    ,F_Y_MEAN(M,L)*1000.d0,F_Z_MEAN(M,L)*1000.d0

	endif;enddo
	close(fileact)

! Structural forces:
	 do L=1,nodes(M)
		RtoHub=R0(M,L)-r_act(1)
	   if(L.le.nodes(M)/3) THEN
		T1_ACT=T1_ACT+FX1(M,L)*dzm*c_act(L)*rho
		Q1_ACT=Q1_ACT+(FX2(M,L)*dcos(rads(M)+alpha0(M,L))
     &    	-FX3(M,L)*dsin(rads(M)+alpha0(M,L)))*R0(M,L)*dzm*c_act(L)*rho
		BR1_ACT=BR1_ACT+FX1(M,L)*dzm*c_act(L)*RtoHub*rho
		BT1_ACT=BT1_ACT+FX1(M,L)*dzm*c_act(L)
     &		*R0(M,L)*rho*dcos(rads(M)+alpha0(M,L))

	   ELSE IF(L.GT.nodes(M)/3 .and. L.le.nodes(M)*2/3 ) then
		T2_ACT=T2_ACT+FX1(M,L)*dzm*c_act(L)*rho
		Q2_ACT=Q2_ACT+(FX2(M,L)*dcos(rads(M)+alpha0(M,L))
     &    	-FX3(M,L)*dsin(rads(M)+alpha0(M,L)))*R0(M,L)*dzm*c_act(L)*rho	
		BR2_ACT=BR2_ACT+FX1(M,L)*dzm*c_act(L)*RtoHub*rho
		BT2_ACT=BT2_ACT+FX1(M,L)*dzm*c_act(L)
     &		*R0(M,L)*rho*dcos(rads(M)+alpha0(M,L))
	   ELSE IF(L.GT.nodes(M)*2/3 .and. L.le.nodes(M)) then
		T3_ACT=T3_ACT+FX1(M,L)*dzm*c_act(L)*rho
		Q3_ACT=Q3_ACT+(FX2(M,L)*dcos(rads(M)+alpha0(M,L))
     &    	-FX3(M,L)*dsin(rads(M)+alpha0(M,L)))*R0(M,L)*dzm*c_act(L)*rho
		BR3_ACT=BR3_ACT+FX1(M,L)*dzm*c_act(L)*RtoHub*rho
		BT3_ACT=BT3_ACT+FX1(M,L)*dzm*c_act(L)
     &		*R0(M,L)*rho*dcos(rads(M)+alpha0(M,L))
	   ENDIF	
		Myaw=Myaw+FX1(M,L)*dzm*c_act(L)
     &		*R0(M,L)*rho*dsin(rads(M)+alpha0(M,L))

		MtwrY=MtwrY+FX1(M,L)*dzm*c_act(L)*rho
     &		*(0.225d0-R0(M,L)*dcos(rads(M)+alpha0(M,L))) !Assume suport is 0.225 above torque shaft

		MtwrX=MtwrX+(FX2(M,L)*dcos(rads(M)+alpha0(M,L))
     &    	-FX3(M,L)*dsin(rads(M)+alpha0(M,L)))*dzm*c_act(L)*rho
     &		*(0.225d0-R0(M,L)*dcos(rads(M)+alpha0(M,L))) !Assume suport is 0.225 above torque shaft

	 end do
	fileact=1200+numIB ; fileact2=2200+numIB
	write(fileact,88)alpharads,T1_ACT,T2_ACT,T3_ACT,Q1_ACT,Q2_ACT
     &	,Q3_ACT,BR1_ACT,BR2_ACT,BR3_ACT,BT1_ACT,BT2_ACT,BT3_ACT
	write(fileact2,88)alpharads,Myaw,MtwrX,MtwrY

   88 FORMAT (1f18.7,18e18.9)
	END SUBROUTINE
```
