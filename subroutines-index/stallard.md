# Stallard

File: actuator.for

```
######################################################################
      SUBROUTINE Stallard(phi)
!######################################################################
      use imb
      Double precision, intent(in) :: phi
      Double precision :: Cl_1,Cl_2,d1,dtot,Cd_1,Cd_2,ang1,ang2

	ang1=0.d0 ;ang2=0.d0 ; Cl_act=0.d0 ; Cd_act=0.d0
	Cl_1=0.d0;Cd_1=0.d0;Cl_2=0.d0;Cd_2=0.d0

      IF ((phi.le.-30.d0)) THEN
	ang1=-100.d0 ;ang2=-30.d0
	Cl_1=-0.7d0;Cd_1=0.48d0;Cl_2=-0.7d0;Cd_2=0.48d0

      ELSE IF ((phi.gt.-30.d0).and.(phi.le.-20.d0)) THEN
	 ang1 = -30.d0 ; ang2 = -20.d0
	 Cl_1=-0.70 ; Cd_1= 0.48 ; Cl_2 = -0.50 ; Cd_2 = 0.25

      ELSE IF ((phi.gt.-20.d0).and.(phi.le.-10.d0)) THEN
	 ang1 = -20.d0 ; ang2 = -10.d0
	 Cl_1=-0.50 ; Cd_1= 0.25 ; Cl_2 = -0.30 ; Cd_2 = 0.13

      ELSE IF ((phi.gt.-10.d0).and.(phi.le.-7.d0)) THEN
	 ang1 = -10.d0 ; ang2 = -7.d0
	 Cl_1=-0.30 ; Cd_1= 0.13 ; Cl_2 = -0.40 ; Cd_2 = 0.12

      ELSE IF ((phi.gt.-7.d0).and.(phi.le.-6.d0)) THEN
	 ang1 = -7.d0 ; ang2 = -6.d0
	 Cl_1=-0.40 ; Cd_1= 0.12 ; Cl_2 = -0.51 ; Cd_2 = 0.11

      ELSE IF ((phi.gt.-6.d0).and.(phi.le.0.d0)) THEN
	 ang1 = -6.d0 ; ang2 = 0.d0
	 Cl_1=-0.51 ; Cd_1= 0.11 ; Cl_2 = 0.24 ; Cd_2 = 0.06

      ELSE IF ((phi.gt.0.d0).and.(phi.le.1.d0)) THEN
	 ang1 = 0.d0 ; ang2 = 1.d0
	 Cl_1= 0.24 ; Cd_1= 0.06 ; Cl_2 = 0.37 ; Cd_2 = 0.06

      ELSE IF ((phi.gt.1.d0).and.(phi.le.2.d0)) THEN
	 ang1 = 1.d0 ; ang2 = 2.d0
	 Cl_1= 0.37 ; Cd_1= 0.06 ; Cl_2 = 0.49 ; Cd_2 = 0.06

      ELSE IF ((phi.gt.2.d0).and.(phi.le.4.d0)) THEN
	 ang1 = 2.d0 ; ang2 = 4.d0
	 Cl_1= 0.49 ; Cd_1= 0.06 ; Cl_2 = 0.74 ; Cd_2 = 0.06

      ELSE IF ((phi.gt.4.d0).and.(phi.le.8.d0)) THEN
	 ang1 = 4.d0 ; ang2 = 8.d0
	 Cl_1= 0.74 ; Cd_1= 0.06 ; Cl_2 = 1.24 ; Cd_2 = 0.10

      ELSE IF ((phi.gt.8.d0).and.(phi.le.9.d0)) THEN
	 ang1 = 8.d0 ; ang2 = 9.d0
	 Cl_1= 1.24 ; Cd_1= 0.10 ; Cl_2 = 1.19 ; Cd_2 = 0.12

      ELSE IF ((phi.gt.9.d0).and.(phi.le.10.d0)) THEN
	 ang1 = 9.d0 ; ang2 = 10.d0
	 Cl_1= 1.19 ; Cd_1= 0.12 ; Cl_2 = 1.11 ; Cd_2 = 0.14

      ELSE IF ((phi.gt.10.d0).and.(phi.le.12.d0)) THEN
	 ang1 = 10.d0 ; ang2 = 12.d0
	 Cl_1= 1.11 ; Cd_1= 0.14 ; Cl_2 = 0.97 ; Cd_2 = 0.23

      ELSE IF ((phi.gt.12.d0).and.(phi.le.14.d0)) THEN
	 ang1 = 12.d0 ; ang2 = 14.d0
	 Cl_1= 0.97 ; Cd_1= 0.23 ; Cl_2 = 0.91; Cd_2 = 0.29

      ELSE IF ((phi.gt.14.d0).and.(phi.le.16.d0)) THEN
	 ang1 = 14.d0 ; ang2 = 16.d0
	 Cl_1= 0.91 ; Cd_1= 0.29 ; Cl_2 = 0.89 ; Cd_2 = 0.35

      ELSE IF ((phi.gt.16.d0).and.(phi.le.18.d0)) THEN
	 ang1 = 16.d0 ; ang2 = 18.d0
	 Cl_1= 0.89 ; Cd_1= 0.35 ; Cl_2 = 0.91 ; Cd_2 = 0.41

      ELSE IF ((phi.gt.18.d0).and.(phi.le.20.d0)) THEN
	 ang1 = 18.d0 ; ang2 = 20.d0
	 Cl_1= 0.91 ; Cd_1= 0.41 ; Cl_2 = 0.94 ; Cd_2 = 0.46

      ELSE IF ((phi.gt.20.d0).and.(phi.le.25.d0)) THEN
	 ang1 = 20.d0 ; ang2 = 25.d0
	 Cl_1= 0.94 ; Cd_1= 0.46 ; Cl_2 = 1.06 ; Cd_2 = 0.59

      ELSE IF ((phi.gt.25.d0).and.(phi.le.30.d0)) THEN
	 ang1 = 25.d0 ; ang2 = 30.d0
	 Cl_1= 1.06 ; Cd_1= 0.59 ; Cl_2 = 1.16 ; Cd_2 = 0.74

      ELSE IF ((phi.gt.30.d0).and.(phi.le.35.d0)) THEN
	 ang1 = 30.d0 ; ang2 = 35.d0
	 Cl_1= 1.16 ; Cd_1= 0.74 ; Cl_2 = 1.19 ; Cd_2 = 0.89

      ELSE IF ((phi.gt.35.d0).and.(phi.le.40.d0)) THEN
	 ang1 = 35.d0 ; ang2 = 40.d0
	 Cl_1= 1.19 ; Cd_1= 0.89 ; Cl_2 = 1.16 ; Cd_2 = 1.01

      ELSE IF ((phi.gt.40.d0).and.(phi.le.45.d0)) THEN
	 ang1 = 40.d0 ; ang2 = 45.d0
	 Cl_1= 1.16; Cd_1= 1.01 ; Cl_2 = 1.11 ; Cd_2 = 1.07

      ELSE IF ((phi.gt.45.d0).and.(phi.le.50.d0)) THEN
	 ang1 = 45.d0 ; ang2 = 50.d0
	 Cl_1= 1.11; Cd_1= 1.07 ; Cl_2 = 1.10 ; Cd_2 = 1.07

      ELSE IF ((phi.gt.50.d0)) THEN
	 ang1 = 50.d0 ; ang2 = 100.d0
	 Cl_1= 1.11; Cd_1= 1.07 ; Cl_2 = 1.11 ; Cd_2 = 1.07
      END IF

	 d1=DABS( ang1 - phi) ;  dtot = ang2 - ang1
	 Cl_act = Cl_1 + (Cl_2 - Cl_1)*d1/dtot
	 Cd_act = Cd_1 + (Cd_2 - Cd_1)*d1/dtot

      RETURN
      End

```

