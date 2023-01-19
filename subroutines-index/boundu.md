# boundu

**File:** bounds.for\
\
**Subroutine Called in:** \
**diffusion** (diffustion.for) -> **flosol** (flosol.for)\
**presure1sweep** (press.for) ->  **flosol** (flosol.for)\
**rungek\_conv4th** (rungek.for) -> **flosol** (flosol.for)\
**rungek\_conv2nd** (rungek.for) -> **flosol** (flosol.for)\
**rungek\_convWENO** (rungek.for) -> **flosol** (flosol.for)\
**newsolv\_mg** (newsolv\_mg.for) -> **flosol** (flosol.for)\
\
**Purpose:**\
****This subroutine allows you to prescribe u velocity boundary condition on all on the computational domain faces at each time step (West-East, South-North, Bottom-Top).\
\
**User Likeness to alter the file:** \
****Often

```
!##############################################################################     
        subroutine boundu
!	Note 1: prescribed inflow always comes from East
!	Note 2: wall functions apply in ALL directions/if node inside vicous layer->no slip
!	Note 3: if LSM/LSMbase are selected, inflow/outflow BCs only applied BELOW surface/lid
!##############################################################################
       use vars
       use multidata
	     use mpi
	     use imb
	      implicit none
	      integer :: i,j,k,ib,ly,im,kk,jj
        integer :: is,ie,js,je,ks,ke,ktop,itop,strlen,dummy
	  double precision ::Fwallu,dyy,dzz,dxx,up,ufric
        double precision ::n_r,U_corr,H,k_w,c,tau,dn,ddn,n_m,dddn,n_h,ep
	  character (LEN=29) :: filename,fileSEM
	  character (LEN=4) :: domain
	  character (LEN=5) :: name_end
        real, parameter :: PI = 4 * atan (1.0)	

        if (PERIODIC) call exchange_bc(1,pl_ex)
        do ly=0,pl_ex
        do ib=1,nbp
           dxx=dom(ib)%dx
           dyy=dom(ib)%dy
           dzz=dom(ib)%dz
           is=dom(ib)%isu; ie=dom(ib)%ieu
           js=dom(ib)%jsu; je=dom(ib)%jeu
           ks=dom(ib)%ksu; ke=dom(ib)%keu

! Boundary Conditions for u - i.e. shear at boundaries and/or diriclet conditions
!..............................................................................
!=== West ===> ..   4=wall  ..    1=Inflow  ..	7=read inflow   .. 31=Linear Waves
!..............................................................................
        if (dom(ib)%iprev.lt.0) then

           if ((dom(ib)%bc_west.eq.4).or.(dom(ib)%bc_east.ge.61)) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%u(is-1-ly,j,k)= 0.0	
              end do; end do

           else if (dom(ib)%bc_west.eq.1) then
              do k=ks-1,ke+1; do j=js-1,je+1
                if (L_LSM) then
			if (dom(ib)%phi(is,j,k) .ge. 0.d0) then
                    dom(ib)%u(is-1-ly,j,k)= ubulk
			else
                    dom(ib)%u(is-1-ly,j,k)= 0.d0
			end if
                elseif (L_LSMbase) then				
			if (dom(ib)%zc(k) .le. length) then
                    dom(ib)%u(is-1-ly,j,k)= ubulk
			else
                    dom(ib)%u(is-1-ly,j,k)= 0.d0
			end if
		    else									
                    dom(ib)%u(is-1-ly,j,k)= ubulk
		    end if
              end do; end do

           elseif (dom(ib)%bc_west.eq.31) then 			!===== Aristos 25/02/19 Solitary Wave=========
              do k=ks-1,ke+1; do j=js-1,je+1
                if (L_LSM) then
		H=0.02d0; k_w=sqrt(3.0d0*H/(4.0d0*length**3.0))     !specify some parameters for solitary wave
		c=sqrt(9.81*(H+length)); tau= 8.0d0; ep=H/length    !ep =epsilon, tau = phase to delay the generation of  waves 
!==============================================================
            n_h = 1.0d0/(cosh(k_w*(-c*ctime)+tau))**2.0
		dn = 2.0d0*c*k_w*tanh(k_w*(-c*ctime)+tau)/(cosh(k_w  ! dn,ddn,dddn = 1st 2nd 3rd derivatives d/dt of n_h
     & *(-c*ctime)+tau))**2.0
		ddn=2.0d0*c**2.0*k_w**2.0*(cosh(2.0d0*tau)-2.0d0)/
     & (cosh(tau))**4.0
		dddn = 4.0d0*c**3.0*k_w**3.0*(cosh(2.0*tau)-5.0d0)
     & *tanh(tau)/(cosh(tau))**4.0
		
		  if (dom(ib)%phi(is,j,k).ge.0.d0) then 
                  dom(ib)%u(is-1-ly,j,k)= sqrt(9.81*length)*ep*
     & (n_h-(1.0d0/4.0d0)*ep*n_h**2.0+(length**2.0d0/(3.0d0*
     & c**2.0))*(1.0d0-(3.0d0*dom(ib)%zc(k)**2.0/(2.0*
     & length**2.0)))*ddn) 
                  else
                  dom(ib)%u(is-1-ly,j,k)= 0.d0
		  end if

	        elseif (L_LSMbase) then
                        if (dom(ib)%zc(k) .le. length) then
                    dom(ib)%u(is-1-ly,j,k)= ubulk
                        else
                    dom(ib)%u(is-1-ly,j,k)= 0.d0
                        end if
                    else
                    dom(ib)%u(is-1-ly,j,k)= ubulk
                    end if
              end do; end do
!==============================================================================


          else if (dom(ib)%bc_west.eq.7) then					!Reading slices
        		write(name_end,'(I5)') ireadinlet
        		strlen=LEN(TRIM(ADJUSTL(name_end)))
       		name_end=REPEAT('0',(5-strlen))//
     &		TRIM(ADJUSTL(name_end)) 				
			write(domain,'(I4)') dom_id(ib)
		  	strlen=LEN(TRIM(ADJUSTL(domain)))
		 	domain=REPEAT('0',(4-strlen))//
     &		TRIM(ADJUSTL(domain)) 
		    filename='inflow/Inlet_'//domain//'_'//name_end//'.dat'
		      open (unit=405, file=filename)
		      do k=ks-1,ke+1; do j=js-1,je+1
				read(405,*)dom(ib)%u(is-1-ly,j,k)
		      end do; end do	
		      close (405)	

          else if (dom(ib)%bc_west.eq.77) then 					!Mapping inflow 
                write(domain,'(I4)') dom_id(ib)
                strlen=LEN(TRIM(ADJUSTL(domain)))
                domain=REPEAT('0',(4-strlen))//
     &                  TRIM(ADJUSTL(domain)) 
           filename='inflow/Mapping_'//domain//'.dat'
                open (unit=405, file=filename)
              do k=ks-1,ke+1; do j=js-1,je+1
                read(405,*)dom(ib)%u(is-1-ly,j,k)
              end do; end do                    
               close(405)

          else if (dom(ib)%bc_west.eq.8) then					!Reading SEM inlet
        		write(name_end,'(I5)') ireadinlet
        		strlen=LEN(TRIM(ADJUSTL(name_end)))
       		name_end=REPEAT('0',(5-strlen))//
     &			TRIM(ADJUSTL(name_end))
			write(domain,'(I4)') dom_id(ib)
        		strlen=LEN(TRIM(ADJUSTL(domain)))
       		domain=REPEAT('0',(4-strlen))//
     &			TRIM(ADJUSTL(domain)) 
		    fileSEM='inflow/Inlet_'//domain//'.dat'
            	open (unit=405, file=fileSEM)
			read(405,*)
			read(405,*)
              do k=ks-1,ke+1; do j=js-1,je+1
	   	 	read(405,*)up,dummy,dummy
			
	  IF (UPROF.eq.12) then 							!1/7th power law inlet condition
	   if (dom(ib)%yc(j).lt.((yen-yst)/2)) then
          dom(ib)%u(is-1-ly,j,k) = ubulk*(1.0d0+1.0d0/7.0d0)
     &	   *(DABS(2*dom(ib)%yc(j)/(yen-yst)))**(1./7.)
	   else
          dom(ib)%u(is-1-ly,j,k) =  ubulk*(1.0d0+1.0d0/7.0d0)
     &    *(DABS(2*((yen-yst)-dom(ib)%yc(j))/(yen-yst)))**(1.d0/7.d0)
	   endif
          dom(ib)%u(is-1-ly,j,k) = dom(ib)%u(is-1-ly,j,k)*(1.0d0+1./7.)
     &	   *(DABS(dom(ib)%zc(k)/(zen-zst)))**(1./7.) +up

	   ELSE IF (UPROF.eq.13) then 						!1/7th power law inlet condition. Only HORIZONTAL
	  if (dom(ib)%yc(j).lt.((yen-yst)/2)) then
          dom(ib)%u(is-1-ly,j,k) = ubulk*(1.0d0+1.0d0/7.0d0)
     &	   *(DABS(2*dom(ib)%yc(j)/(yen-yst)))**(1./7.)
	  else
          dom(ib)%u(is-1-ly,j,k) = ubulk*(1.0d0+1.0d0/7.0d0)
     &    *(DABS(2*((yen-yst)-dom(ib)%yc(j))/(yen-yst)))**(1.d0/7.d0)
	  endif
          dom(ib)%u(is-1-ly,j,k) = dom(ib)%u(is-1-ly,j,k)+up

	   ELSE IF (UPROF.eq.14) then 						!1/7th power law Only vertical
          dom(ib)%u(is-1-ly,j,k) = ubulk*(1.0d0+1./7.)
     &	   *(DABS(dom(ib)%zc(k)/(zen-zst)))**(1./7.)+up

	   ELSE IF (UPROF.eq.15) then 						!Logarithmic distribution
          dom(ib)%u(is-1-ly,j,k) = 0.0187d0*
     &  (1.d0/0.41d0*DLOG(DABS(dom(ib)%zc(k))*0.0187d0*10**6)+5.d0)+up

	  else 
          dom(ib)%u(is-1-ly,j,k) = ubulk+up
	  endif
               end do ; end do		
			close(405)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   else if (dom(ib)%bc_west.eq.12) then 					!1/7th power law inlet condition 
            do k=ks-1,ke+1; do j=js-1,je+1
	  if (dom(ib)%yc(j).lt.((yen-yst)/2)) then
          dom(ib)%u(is-1-ly,j,k) = ubulk*(1.0d0+1.0d0/7.0d0)
     &	   *(DABS(2*dom(ib)%yc(j)/(yen-yst)))**(1./7.)
	  else
          dom(ib)%u(is-1-ly,j,k) = ubulk*(1.0d0+1.0d0/7.0d0)
     &    *(DABS(2*((yen-yst)-dom(ib)%yc(j))/(yen-yst)))**(1.d0/7.d0)
	  endif
          dom(ib)%u(is-1-ly,j,k) = dom(ib)%u(is-1-ly,j,k)*(1.0d0+1./7.)
     &	   *(DABS(dom(ib)%zc(k)/(zen-zst)))**(1./7.)
	     enddo ; end do 

	   else if (dom(ib)%bc_west.eq.13) then 					!1/7th power law inlet condition. Only horizontal
            do k=ks-1,ke+1; do j=js-1,je+1
	  if (dom(ib)%yc(j).lt.((yen-yst)/2)) then
          dom(ib)%u(is-1-ly,j,k) = ubulk*(1.0d0+1.0d0/7.0d0)
     &	   *(DABS(2*dom(ib)%yc(j)/(yen-yst)))**(1./7.)
	  else
          dom(ib)%u(is-1-ly,j,k) = ubulk*(1.0d0+1.0d0/7.0d0)
     &    *(DABS(2*((yen-yst)-dom(ib)%yc(j))/(yen-yst)))**(1.d0/7.d0)
	  endif
	     enddo ; end do 

	   else if (dom(ib)%bc_west.eq.14) then 					!1/7th power law inlet condition. Only vertical
            do k=ks-1,ke+1; do j=js-1,je+1
          dom(ib)%u(is-1-ly,j,k) = ubulk*(1.0d0+1./7.)
     &	   *(DABS(dom(ib)%zc(k)/(zen-zst)))**(1./7.)
	     enddo ; end do 

          else if (dom(ib)%bc_west.eq.15) THEN					!Logarithmic distribution
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%u(is-1-ly,j,k)=0.0187d0*
     &  (1.d0/0.41d0*DLOG(DABS(dom(ib)%zc(k))*0.0187d0*10**6)+5.d0)
              end do; end do

          else if (dom(ib)%bc_west.eq.17) THEN			
              do k=ks-1,ke+1; do j=js-1,je+1
 			ufric=0.103594d0*ubulk+0.00568d0
                 dom(ib)%u(is-1-ly,j,k)=ufric*
     &  (1.d0/0.41d0*DLOG(DABS(dom(ib)%zc(k))*ufric*10**6))/LOG(10.d0)
              end do; end do

           end if
        end if
!...............................................................................
!=== East ===> ..  4=wall  ..   2=Outflow
!...............................................................................
        if (dom(ib)%inext.lt.0) then
           if ((dom(ib)%bc_east.eq.4).or.(dom(ib)%bc_east.ge.61)) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%u(ie+1+ly,j,k)= 0.0	
              end do; end do

           else if (dom(ib)%bc_east.eq.2 .or.dom(ib)%bc_east.eq.21) then
              if (alfabc.eq.1)then
                do k=ks-1,ke+1; do j=js-1,je+1

                if (L_LSM) then

                    if (dom(ib)%bc_east.eq.2) then
                       dom(ib)%u(ie+1+ly,j,k)=dom(ib)%u(ie,j,k)
                    else
                       dom(ib)%u(ie+1+ly,j,k)=dom(ib)%uoo(ie+1+ly,j,k) 
     & -dt*alfapr*(dom(ib)%uoo(ie+1+ly,j,k)-
     & dom(ib)%uoo(ie,j,k))/dom(ib)%dx
                    end if

		    else		!no LSM

                    if (dom(ib)%bc_east.eq.2) then
                       dom(ib)%u(ie+1+ly,j,k)=dom(ib)%u(ie,j,k)
                    else
                       dom(ib)%u(ie+1+ly,j,k)=dom(ib)%uoo(ie+1+ly,j,k) 
     & -dt*alfapr*(dom(ib)%uoo(ie+1+ly,j,k)-
     & dom(ib)%uoo(ie,j,k))/dom(ib)%dx
                    end if

		    end if
                    end do; end do
		  endif
           end if
        end if
!...............................................................................
!=== South ===> ..   4=wall ..   3=Symmetry .. 6=Wall function
!...............................................................................
        if (dom(ib)%jprev.lt.0) then

           if (dom(ib)%bc_south.eq.4) then 
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%u(i,js-1-ly,k)= -dom(ib)%u(i,js+ly,k)	
              end do; end do                      
           else if (dom(ib)%bc_south.eq.3) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%u(i,js-1-ly,k)= dom(ib)%u(i,js+ly,k)
              end do; end do
	     else if (dom(ib)%bc_south.ge.61) then				!Wall functions
		if (ly.eq.0) then
	   		if (dom(ib)%bc_south.lt.63) 
     &			call log_law(3,dom(ib)%bc_south)
			if (dom(ib)%bc_south.ge.63) 
     &			call wall_function(3,dom(ib)%bc_south)
		      do k=ks-1,ke+1; do i=is-1,ie+1
				Fwallu = dom(ib)%tauws2(i,k)
     &					*dom(ib)%u(i,js,k)/dyy			!*Acell/Vcell
				Fwallu=SIGN(Fwallu,-dom(ib)%u(i,js,k))
				dom(ib)%u(i,js,k)=dom(ib)%u(i,js,k) 
     &						+Fwallu*alfapr*dt
		           	dom(ib)%u(i,js-1,k)= dom(ib)%u(i,js,k)	
		      end do; end do
		else
		      do k=ks-1,ke+1; do i=is-1,ie+1
		           	dom(ib)%u(i,js-1-ly,k)= dom(ib)%u(i,js,k)	
		      end do; end do
		endif
	     endif
        end if
!.............................................................................
!=== North ===>  ..   4=wall ..   44=moving wall ..  3=Symmetry
!.............................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%bc_north.eq.4) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%u(i,je+1+ly,k)   = - dom(ib)%u(i,je-ly,k) 	
              end do; end do

           else if (dom(ib)%bc_north.eq.3) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%u(i,je+1+ly,k)   =  dom(ib)%u(i,je-ly,k) 
              end do; end do
	     else if (dom(ib)%bc_north.ge.61) then				!Wall functions 
		if (ly.eq.0) then
	   		if (dom(ib)%bc_north.lt.63) 
     &			call log_law(4,dom(ib)%bc_north)
			if (dom(ib)%bc_north.ge.63) 
     &			call wall_function(4,dom(ib)%bc_north)
		      do k=ks-1,ke+1; do i=is-1,ie+1
				Fwallu = dom(ib)%tauwn2(i,k)
     &					*dom(ib)%u(i,je,k)/dyy			!*Acell/Vcell
				Fwallu=sign(Fwallu,-dom(ib)%u(i,je,k))
				dom(ib)%u(i,je,k)=dom(ib)%u(i,je,k) 
     &						+Fwallu*alfapr*dt
		           	dom(ib)%u(i,je+1,k)= dom(ib)%u(i,je,k)	
		      end do; end do
		else
		      do k=ks-1,ke+1; do i=is-1,ie+1
		           	dom(ib)%u(i,je+1+ly,k)= dom(ib)%u(i,je,k)	
		      end do; end do
		endif
	     endif
        end if
!...............................................................................
!=== Bottom ===> ..   4=wall ..   3=Symmetry
!...............................................................................
        if (dom(ib)%kprev.lt.0) then
          if (dom(ib)%bc_bottom.eq.4) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%u(i,j,ks-1-ly)= -dom(ib)%u(i,j,ks+ly)	
              end do; end do

           else if (dom(ib)%bc_bottom.eq.1) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%u(i,j,ks-1-ly)= 0.0
              end do; end do

           else if (dom(ib)%bc_bottom.eq.3) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%u(i,j,ks-1-ly)= dom(ib)%u(i,j,ks+ly)
              end do; end do
	     else if (dom(ib)%bc_bottom.ge.61) then				!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_bottom.lt.63) then
     				call log_law(5,dom(ib)%bc_bottom)
			endif
			if (dom(ib)%bc_bottom.ge.63) 
     &			call wall_function(5,dom(ib)%bc_bottom)
		      do j=js-1,je+1; do i=is-1,ie+1
				Fwallu = dom(ib)%tauwb2(i,j)
     &					*dom(ib)%u(i,j,ks)/dzz			!*Acell/Vcell??
				Fwallu = sign(Fwallu,-dom(ib)%u(i,j,ks))

				dom(ib)%u(i,j,ks)=dom(ib)%u(i,j,ks) 
     &						+Fwallu*alfapr*dt

		           	dom(ib)%u(i,j,ks-1)= dom(ib)%u(i,j,ks)	
		      end do; end do
		else
		      do j=js-1,je+1; do i=is-1,ie+1
		           	dom(ib)%u(i,j,ks-1-ly)= dom(ib)%u(i,j,ks)
		      end do; end do
		endif
           end if
        end if
!.............................................................................
!=== Top ===>  ..   4=wall ..     3=Symmetry
!.............................................................................
	  if (L_LSMbase) then
        if (dom(ib)%z(1).le.length .and. dom(ib)%z(ke).ge.length) then		
	    k=1
          do while (dom(ib)%z(k).le.length)
	      k=k+1
	    end do
	    ktop=k-1
          do j=js-1,je+1; do i=is-1,ie+1
             dom(ib)%u(i,j,ktop) = dom(ib)%u(i,j,ktop-1) 	!THIS IS A SLIP CONDITION
             dom(ib)%u(i,j,ktop+1) = dom(ib)%u(i,j,ktop-2) 	!THIS IS A SLIP CONDITION
          enddo; enddo
	  !pablo:	03/19
	   	do k=ktop+2,ke	!Over the free-surface layer
		    do j=js-1,je+1; do i=is-1,ie+1
		       dom(ib)%u(i,j,k) = 0.0
		    enddo; enddo
		enddo
	  end if
        if (dom(ib)%z(1).ge.length) then		
	    do k=1,ke
           do j=js-1,je+1; do i=is-1,ie+1
             dom(ib)%u(i,j,k) = 0.d0 !ENFORCE VELOCITIES TO ZERO IN THE AIR
          enddo; enddo; enddo
	  end if
	  end if

        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%bc_top.eq.4) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%u(i,j,ke+1+ly)   = - dom(ib)%u(i,j,ke-ly)	
              end do; end do
           else if (dom(ib)%bc_top.eq.3) then
             do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%u(i,j,ke+1+ly)= dom(ib)%u(i,j,ke-ly)
             end do; end do
	     else if (dom(ib)%bc_top.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_top.lt.63) 
     &			call log_law(6,dom(ib)%bc_top)
			if (dom(ib)%bc_top.ge.63) 
     &			call wall_function(6,dom(ib)%bc_top)
		      do j=js-1,je+1; do i=is-1,ie+1
				Fwallu = dom(ib)%tauwt2(i,j)
     &					*dom(ib)%u(i,j,ke)/dzz			!*Acell/Vcell
				Fwallu = sign(Fwallu,-dom(ib)%u(i,j,ke))
				dom(ib)%u(i,j,ke)=dom(ib)%u(i,j,ke) 
     &						+Fwallu*alfapr*dt
                 		dom(ib)%u(i,j,ke+1) = dom(ib)%u(i,j,ke)	
		      end do; end do
		else
		      do j=js-1,je+1; do i=is-1,ie+1
                 		dom(ib)%u(i,j,ke+1+ly) = dom(ib)%u(i,j,ke)	
		      end do; end do
		endif
           end if

          if (trim(keyword).eq.'cavity') then	
		  do j=js-1,je+1
		   do i=is-1,ie+1
!                 dom(ib)%u(i,j,ke-ly)  = 1.d0	
                 dom(ib)%u(i,j,ke+1+ly)  = ubulk
              end do
		 	end do
	   	  endif
	  end if

        end do
        if(diff_sch.ne.3) call boundcoef(1)
        end do

	if (LHRS) call bc_human_u	  ! LHRS Alex & Yan

        end subroutine
```
