# act\_line\_geom

File: actuator.for

```
!######################################################################
      SUBROUTINE  act_line_geom(numIB)
!######################################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: numIB
      DOUBLE PRECISION:: PI,an,d1,Dtot,alph,dummy
      INTEGER      :: L,nin,I,K,nlay,strlen,n_act,BladeN,fileact
      CHARACTER*8  :: char_block
      CHARACTER*31 :: gridfile,gridfile2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::r_in,c_in,Pit_in

       PI = 4.D0*DATAN(1.D0)

         write(char_block,'(i3)') numIB
         strlen=LEN(TRIM(ADJUSTL(char_block)))
         char_block=REPEAT('0',(3-strlen))//TRIM(ADJUSTL(char_block))
         gridfile='geom_ActL_'//TRIM(ADJUSTL(char_block))//'.dat'
         open (unit=2, file=gridfile)
!----      Load the file and proceed to interpolate     ----------
	 open(unit=1, file=filepoints(numIB))
	   read(1,*)n_act
	   read(1,*)
	   allocate(r_in(n_act),c_in(n_act),Pit_in(n_act))
	 DO L=1,n_act
	   read(1,*)r_in(L),c_in(L),Pit_in(L)
       ENDDO
      close (1)

	 BladeN=imbnumber(numIB)				!!!*****
	 nin=INT((r_in(n_act)-r_in(1))/dzm)+1

	 r_act(1+maxnodeIBS)  =r_in(1)      !+0.3999999d0*dzm
	 c_act(1+maxnodeIBS)  =c_in(1)
       Pit_act(1+maxnodeIBS)=Pit_in(1) 
	 nodex(numIB,1)=0.D0
	 nodey(numIB,1)=0.D0
	 nodez(numIB,1)=r_act(1+maxnodeIBS)

	 DO L = 2,nin
		r_act(L+maxnodeIBS)  = r_act(L-1+maxnodeIBS)+dzm
	  Do K =1,n_act-1
	   if(r_act(L+maxnodeIBS).le.r_act(1+maxnodeIBS)) then
	    	c_act(L+maxnodeIBS)  = c_in(K)
		Pit_act(L+maxnodeIBS)=Pit_in(K)
	    	nodex(numIB,L)=0.D0
	   	nodey(numIB,L)=0.D0
	   	nodez(numIB,L)=r_act(L+maxnodeIBS)
	   endif
	   if(r_act(L+maxnodeIBS).gt.r_in(K) 
     & .and.r_act(L+maxnodeIBS).le.r_in(K+1)) then
	  d1=DABS(r_in(K)-r_act(L+maxnodeIBS)) ;  dtot=r_in(K+1)-r_in(K) 
	  c_act(L+maxnodeIBS)  = c_in(K)  +(c_in(K+1)-c_in(K))*d1/dtot
	  Pit_act(L+maxnodeIBS)= Pit_in(K)+(Pit_in(K+1)-Pit_in(K))*d1/dtot
	  nodex(numIB,L)=0.D0
	  nodey(numIB,L)=0.D0
	  nodez(numIB,L)=r_act(L+maxnodeIBS)
		WRITE(6,'(i4,3f13.5)')L,r_act(L+maxnodeIBS)
     &   ,c_act(L+maxnodeIBS),Pit_act(L+maxnodeIBS)
	   endif

 	  enddo
	 ENDDO

!Multiple blades generation
	  write(6,*)'Total # points of one blade:',nin,BladeN*nin	
	  alph=2.D0*PI/BladeN  ; k=nin
	DO L=1,BladeN-1
	 DO I=1,nin
	  K=K+1
	  nodex(numIB,K)= nodex(numIB,nin*(L-1)+I)
	  nodey(numIB,K)= nodey(numIB,nin*(L-1)+I)*DCOS(alph)
     & 		     +nodez(numIB,nin*(L-1)+I)*DSIN(alph)
	  nodez(numIB,K)=-nodey(numIB,nin*(L-1)+I)*DSIN(alph)
     & 		     +nodez(numIB,nin*(L-1)+I)*DCOS(alph) 
	  r_act(K+maxnodeIBS)	= r_act(I+maxnodeIBS)
	  c_act(K+maxnodeIBS)   = c_act(I+maxnodeIBS)
	  Pit_act(K+maxnodeIBS) = Pit_act(I+maxnodeIBS)
	 ENDDO	
	ENDDO
	 maxnode=0  ; nodes(numIB)=nin*BladeN
	 Cxor(numIB)=Cx(numIB);Cyor(numIB)=Cy(numIB);Czor(numIB)=Cz(numIB)	
	 maxnode = max(maxnode,nodes(numIB))
   	   DO i=1,nodes(numIB)		
	    nodexlocal(numIB,i)=nodex(numIB,i)+1.d-7  !Set this as the local coordinates
	    nodeylocal(numIB,i)=nodey(numIB,i)+1.d-7 
	    nodezlocal(numIB,i)=nodez(numIB,i)+1.d-7 
	    nodex(numIB,i)=nodex(numIB,i)+Cxor(numIB)
	    nodey(numIB,i)=nodey(numIB,i)+Cyor(numIB)
	    nodez(numIB,i)=nodez(numIB,i)+Czor(numIB)
	   enddo   

	 write(2,*)'VARIABLES=x,y,z,Radius,Chord,Pitch'
	  do L=1,nodes(numIB)
	   write(2,89)nodex(numIB,L),nodey(numIB,L),nodez(numIB,L)
     &  ,r_act(L+maxnodeIBS),c_act(L+maxnodeIBS),Pit_act(L+maxnodeIBS)
	  enddo
	 close(2) 

	   deallocate(r_in,c_in,Pit_in)

!Open file to output structural forces:
	   write(char_block,'(i3)') numIB
         strlen=LEN(TRIM(ADJUSTL(char_block)))
         char_block=REPEAT('0',(3-strlen))//TRIM(ADJUSTL(char_block))
         gridfile='FEM_ActL_'//TRIM(ADJUSTL(char_block))//'.dat'
         gridfile2='FEM_ActT_'//TRIM(ADJUSTL(char_block))//'.dat'
	   fileact=1200+numIB
         open (unit=fileact, file=gridfile)
	   write(fileact,*)'Title = Structural forces Actuator blade'	
	   write(fileact,*)'Variables=Deg,T1,T2,T3,Q1,Q2,Q3'
     &  ,',BR1,BR2,BR3,BT1,BT2,BT3'

	   fileact=2200+numIB
         open (unit=fileact, file=gridfile2)
	   write(fileact,*)'Title = Structural forces Actuator turbine'	
	   write(fileact,*)'Variables=Deg,Myaw,MtwrX,MtwrY'
                           
   87 FORMAT (a,i2,a,i6)
   89 FORMAT (10f13.6)
      RETURN
      end
```
