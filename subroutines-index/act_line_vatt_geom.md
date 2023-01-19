# act\_line\_vatt\_geom

**File:** actuator.for\
\
**Subroutine Called in:** \


\


```
!######################################################################
      SUBROUTINE  act_line_vatt_geom(numIB)
!######################################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: numIB
      DOUBLE PRECISION:: PI,an,d1,Dtot,alph,dummy,c_in
      INTEGER      :: L,nin,I,K,nlay,strlen,n_act,BladeN,fileact
      CHARACTER*8  :: char_block
      CHARACTER*31 :: gridfile

       PI = 4.D0*DATAN(1.D0)

! Input parameters
	 c_in=0.14d0
	 BladeN=imbnumber(numIB)				!!!*****
	 nin=INT((zend(numIB)-zini(numIB))/dzm)+1

	 DO L = 1,nin
		 nodex(numIB,L)=0.D0		
		 nodey(numIB,L)=R(numIB)	
		 nodez(numIB,L)=zini(numIB)+(L-1)*dzm+0.09999d0*dzm	
	 ENDDO

!Multiple blades generation
	  write(6,*)'Total # points of one blade:',nin,BladeN*nin	
	  alph=2.D0*PI/BladeN  ; k=nin+maxnodeIBS
	DO L=1,BladeN-1
	 DO I=1,nin
	  K=K+1
	  nodex(numIB,K)= nodex(numIB,nin*(L-1)+I)*DCOS(alph)
     & 		     +nodey(numIB,nin*(L-1)+I)*DSIN(alph)
	  nodey(numIB,K)=-nodex(numIB,nin*(L-1)+I)*DSIN(alph)
     & 		     +nodey(numIB,nin*(L-1)+I)*DCOS(alph) 
	  nodez(numIB,K)= nodez(numIB,nin*(L-1)+I)
	 ENDDO	
	ENDDO
	 maxnode=0  ; nodes(numIB)=nin*BladeN
	 Cxor(numIB)=Cx(numIB);Cyor(numIB)=Cy(numIB);Czor(numIB)=Cz(numIB)	
	 maxnode = max(maxnode,nodes(numIB))
   	   DO L=1,nodes(numIB)		
	    nodexlocal(numIB,L)=nodex(numIB,L)+1.d-7  !Set this as the local coordinates
	    nodeylocal(numIB,L)=nodey(numIB,L)+1.d-7 
	    nodezlocal(numIB,L)=nodez(numIB,L)+1.d-7 
	    nodex(numIB,L)=nodex(numIB,L)+Cxor(numIB) !Global coordinates
	    nodey(numIB,L)=nodey(numIB,L)+Cyor(numIB)
	    nodez(numIB,L)=nodez(numIB,L)+Czor(numIB)
	   enddo   

	 write(2,*)'VARIABLES=x,y,z'
	  do L=1,nodes(numIB)
	   write(2,89)nodex(numIB,L),nodey(numIB,L),nodez(numIB,L)
	  enddo
	 close(2) 

           
   87 FORMAT (a,i2,a,i6)
   89 FORMAT (10f13.6)
      RETURN
      end
```
