# boundw

**File:** bounds.for\
\
**Subroutine Called in:** \
diffusion (diffustion.for) -> flosol (flosol.for)\
presure1sweep (press.for) ->  flosol (flosol.for)\
rungek_conv4th (rungek.for) -> flosol (flosol.for)_\
_rungekconv2nd (rungek.for) -> flosol (flosol.for)_\
_rungekconvWENO (rungek.for) -> flosol (flosol.for)_\
_newsolv\_mg (newsolv\_mg.for) -> flosol (flosol.for)_\
__\
__**Purpose:**\
****This subroutine allows you to prescribe w velocity boundary condition on all on the computational domain faces (West-East, South-North, Bottom-Top).\
\
**User Likeness to alter the file:** \
****Often

```
!#############################################################################
        subroutine boundw
!#############################################################################
        use vars
        use multidata
	  use mpi
	  use imb
        implicit none
        integer :: i,j,k,ib,ly,kk,jj,is,ie,js,je,ks,ke,ktop,strlen
	  character (LEN=29):: filename,fileSEM
	  character (LEN=4) :: domain
	  character (LEN=5) :: name_end
	  double precision :: Fwallw,dxx,dyy,dummy,dzz
	  double precision :: H,k_w,c,tau,dn,ddn,dddn,n_h,ep

        if (PERIODIC) call exchange_bc(3,pl_ex)

        do ly=0,pl_ex
        do ib=1,nbp
	     dxx=dom(ib)%dx
	     dyy=dom(ib)%dy
	     dzz=dom(ib)%dz
           is=dom(ib)%isw; ie=dom(ib)%iew
           js=dom(ib)%jsw; je=dom(ib)%jew
           ks=dom(ib)%ksw; ke=dom(ib)%kew

! Boundary Conditions for W - i.e. shear at boundaries
!...............................................................................
!=== West ===>   4=wall   ..    1=Inflow
!...............................................................................
        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%bc_west.eq.4) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%w(is-1-ly,j,k)= -dom(ib)%w(is+ly,j,k)	
              end do; end do

           else if (dom(ib)%bc_west.eq.1 .or. dom(ib)%bc_west.eq.12
     &  .or. dom(ib)%bc_west.eq.14 .or. dom(ib)%bc_west.eq.15
     &  .or. dom(ib)%bc_west.eq.17)then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%w(is-1-ly,j,k)= 0.d0
              end do; end do
	     else if (dom(ib)%bc_west.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_west.lt.63) 
     &			call log_law(1,dom(ib)%bc_west)
			if (dom(ib)%bc_west.ge.63) 
     &			call wall_function(1,dom(ib)%bc_west)

		      do k=ks-1,ke+1; do j=js-1,je+1
				Fwallw = dom(ib)%tauww2(j,k)
     &					*dom(ib)%w(is,j,k)/dxx			!*Acell/Vcell
				Fwallw = sign(Fwallw,-dom(ib)%w(is,j,k))
				dom(ib)%w(is,j,k)=dom(ib)%w(is,j,k) 
     &						+Fwallw*alfapr*dt
                 		dom(ib)%w(is-1,j,k)= dom(ib)%w(is,j,k)	
			enddo; enddo
		else
		      do k=ks-1,ke+1; do j=js-1,je+1
                 		dom(ib)%w(is-1-ly,j,k)= dom(ib)%w(is,j,k)	
		      end do; end do
		endif

           elseif (dom(ib)%bc_west.eq.31) then    !===Aristos Christou Soitary Wave
              do k=ks-1,ke+1; do j=js-1,je+1

                H=0.02d0; k_w=sqrt(3.0d0*H/(4.0d0*length**3.0))         !specify some parameters for solitary wave
                c=sqrt(9.81*(H+length)); tau= 8.0d0; ep=H/length
!===============================================================
                n_h = 1.0d0/(cosh(k_w*(-c*ctime)+tau))**2.0
                dn = 2.0d0*c*k_w*tanh(k_w*(-c*ctime)+tau)/(cosh(k_w    ! dn,ddn,dddn = 1st 2nd 3rd derivatives d/dt of n_h
     & *(-c*ctime)+tau))**2.0
                ddn=2.0d0*c**2.0*k_w**2.0*(cosh(2.0d0*tau)-2.0d0)/
     & (cosh(tau))**4.0
                dddn = 4.0d0*c**3.0*k_w**3.0*(cosh(2.0*tau)-5.0d0)
     & *tanh(tau)/(cosh(tau))**4.0

                if (L_LSM) then
                  if (dom(ib)%phi(is,j,k) .ge. 0.d0) then
                  dom(ib)%w(is-1-ly,j,k)= dom(ib)%z(k)/c*
     & sqrt(9.81*length)*ep*((1.0d0-0.5d0*ep*n_h)*dn+
     & (length**2.0/(3.0d0*c**2.0))*(1.0d0-(dom(ib)%zc(k)**2.0
     & /(2.0d0*length**2.0)))*dddn)
		  else
                  dom(ib)%w(is-1-ly,j,k)= 0.d0
                  end if

                else if (L_LSMbase) then
                        if (dom(ib)%zc(k) .le. length) then
                    dom(ib)%u(is-1-ly,j,k)= ubulk
                        else
                    dom(ib)%u(is-1-ly,j,k)= 0.d0
                        end if
                    else
                    dom(ib)%u(is-1-ly,j,k)= ubulk
                    end if
              end do; end do
!===============================================================================

          else if (dom(ib)%bc_west.eq.7) then					!brunho2014 reading slices
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
		      do k=dom(ib)%ksu-1,dom(ib)%keu+1
			do j=dom(ib)%jsu-1,dom(ib)%jeu+1
				read(405,*)dummy
		      end do; end do	
		      do k=dom(ib)%ksv-1,dom(ib)%kev+1
			do j=dom(ib)%jsv-1,dom(ib)%jev+1
				read(405,*)dummy
		      end do; end do	
		      do k=ks-1,ke+1; do j=js-1,je+1
				read(405,*)dom(ib)%w(is-1-ly,j,k)
			end do; end do	
			close(405)	

          else if (dom(ib)%bc_west.eq.77) then 					!Mapping inflow
                write(domain,'(I4)') dom_id(ib)
                strlen=LEN(TRIM(ADJUSTL(domain)))
                domain=REPEAT('0',(4-strlen))//
     &                  TRIM(ADJUSTL(domain))
           filename='inflow/Mapping_'//domain//'.dat'
                open (unit=405, file=filename)
		      do k=dom(ib)%ksu-1,dom(ib)%keu+1
			do j=dom(ib)%jsu-1,dom(ib)%jeu+1
				read(405,*)dummy
		      end do; end do	
		      do k=dom(ib)%ksv-1,dom(ib)%kev+1
			do j=dom(ib)%jsv-1,dom(ib)%jev+1
				read(405,*)dummy
		      end do; end do	
              do k=ks-1,ke+1; do j=js-1,je+1
			read(405,*)dom(ib)%w(is-1-ly,j,k)
              end do; end do                    
               close(405)		
				
          else if (dom(ib)%bc_west.eq.8) then					!SEM 
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
	         read(405,*)dummy,dummy,dom(ib)%w(is-1-ly,j,k)
              end do; end do			
			close(405)

           end if
        end if
!...............................................................................
!=== East ===>   4=wall   ..    2=Outflow
!...............................................................................
        if (dom(ib)%inext.lt.0) then
           if (dom(ib)%bc_east.eq.4) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%w(ie+1+ly,j,k)= -dom(ib)%w(ie-ly,j,k)	
              end do; end do

           elseif (dom(ib)%bc_east.eq.2 .or.dom(ib)%bc_east.eq.21) then
              if (alfabc.eq.1) then
                 do k=ks-1,ke+1; do j=js-1,je+1

                   if (L_LSM) then

			   if (dom(ib)%phi(ie,j,k) .ge. 0.0) then
                       if (dom(ib)%bc_east.eq.2) then
                       dom(ib)%w(ie+1+ly,j,k)= dom(ib)%w(ie,j,k)
                       else
                       dom(ib)%w(ie+1+ly,j,k)= dom(ib)%woo(ie+1+ly,j,k)
     & -dt*alfapr*(dom(ib)%woo(ie+1+ly,j,k)-
     & dom(ib)%woo(ie,j,k))/dom(ib)%dx
                       end if
			   else
                       dom(ib)%w(ie+1+ly,j,k)= dom(ib)%w(ie,j,k) 
			   endif

			 else	!no LSM
                       if (dom(ib)%bc_east.eq.2) then
                       dom(ib)%w(ie+1+ly,j,k)= dom(ib)%w(ie,j,k)
                       else
                       dom(ib)%w(ie+1+ly,j,k)= dom(ib)%woo(ie+1+ly,j,k)
     & -dt*alfapr*(dom(ib)%woo(ie+1+ly,j,k)-
     & dom(ib)%woo(ie,j,k))/dom(ib)%dx
                       end if
			 endif

                 end do; end do
              end if
	     else if (dom(ib)%bc_east.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_east.lt.63) 
     &			call log_law(2,dom(ib)%bc_east)
			if (dom(ib)%bc_east.ge.63) 
     &			call wall_function(2,dom(ib)%bc_east)

		      do k=ks-1,ke+1; do j=js-1,je+1
				Fwallw = dom(ib)%tauwe2(j,k)
     &					*dom(ib)%w(ie,j,k)/dxx			!*Acell/Vcell
				Fwallw = sign(Fwallw,-dom(ib)%w(ie,j,k))
				dom(ib)%w(ie,j,k)=dom(ib)%w(ie,j,k) 
     &						+Fwallw*alfapr*dt
                 		dom(ib)%w(ie+1,j,k)= dom(ib)%w(ie,j,k)	
		      end do; end do
		else
		      do k=ks-1,ke+1; do j=js-1,je+1
                 		dom(ib)%w(ie+1+ly,j,k)= dom(ib)%w(ie,j,k)	
		      end do; end do
		endif
           end if
        end if
!...............................................................................
!=== South ===> ..   4=wall ..   3=Symmetry
!...............................................................................
        if (dom(ib)%jprev.lt.0) then
           if (dom(ib)%bc_south.eq.4) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%w(i,js-1-ly,k)= -dom(ib)%w(i,js+ly,k)	
              end do; end do

           else if (dom(ib)%bc_south.eq.3) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%w(i,js-1-ly,k)= dom(ib)%w(i,js+ly,k)
              end do; end do
	     else if (dom(ib)%bc_south.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_south.lt.63) 
     &			call log_law(3,dom(ib)%bc_south)
			if (dom(ib)%bc_south.ge.63) 
     &			call wall_function(3,dom(ib)%bc_south)
		      do k=ks-1,ke+1; do i=is-1,ie+1
				Fwallw = dom(ib)%tauws2(i,k)
     &					*dom(ib)%w(i,js,k)/dyy			!*Acell/Vcell
				Fwallw = sign(Fwallw,-dom(ib)%w(i,js,k))
				dom(ib)%w(i,js,k)=dom(ib)%w(i,js,k) 
     &						+Fwallw*alfapr*dt
		           	dom(ib)%w(i,js-1,k)= dom(ib)%w(i,js,k)	
		      end do; end do
		else
		      do k=ks-1,ke+1; do i=is-1,ie+1
		           	dom(ib)%w(i,js-1-ly,k)= dom(ib)%w(i,js,k)	
		      end do; end do
		endif
           end if
        end if
!.............................................................................
!=== North ===>  ..   4=wall ..   44=moving wall ..  3=Symmetry
!.............................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%bc_north.eq.4) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%w(i,je+1+ly,k)   = - dom(ib)%w(i,je-ly,k)	 
              end do; end do
           else if (dom(ib)%bc_north.eq.3) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%w(i,je+1+ly,k)   =  dom(ib)%w(i,je-ly,k) 
              end do; end do
	     else if (dom(ib)%bc_north.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_north.lt.63) 
     &			call log_law(4,dom(ib)%bc_north)
			if (dom(ib)%bc_north.ge.63) 
     &			call wall_function(4,dom(ib)%bc_north)
		      do k=ks-1,ke+1; do i=is-1,ie+1
				Fwallw = dom(ib)%tauwn2(i,k)
     &					*dom(ib)%w(i,je,k)/dyy			!*Acell/Vcell
				Fwallw = sign(Fwallw,-dom(ib)%w(i,je,k))
				dom(ib)%w(i,je,k)=dom(ib)%w(i,je,k) 
     &						+Fwallw*alfapr*dt
		           	dom(ib)%w(i,je+1,k)= dom(ib)%w(i,je,k)	
		      end do; end do
		else
		      do k=ks-1,ke+1; do i=is-1,ie+1
				dom(ib)%w(i,je+1+ly,k)= dom(ib)%w(i,je,k)	
		      end do; end do
		endif
           end if
        end if
!...............................................................................
!=== Bottom ===>   4=wall   ..    3=Symmetry
!...............................................................................
        if (dom(ib)%kprev.lt.0) then 
           if (dom(ib)%bc_bottom.eq.3) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%w(i,j,ks-1-ly)=0.0
              end do; end do

		else if ((dom(ib)%bc_bottom.eq.4).or.
     &		   (dom(ib)%bc_bottom.ge.61)) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%w(i,j,ks-1-ly)=0.0	
              end do; end do

           else if (dom(ib)%bc_bottom.eq.1) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%w(i,j,ks-1-ly)=1.1
              end do; end do
           end if
        end if
!...............................................................................
!=== Top ===>   4=wall   ..    3=Symmetry
!...............................................................................
        if (L_LSMbase) then
        if (dom(ib)%z(1).le.length .and. dom(ib)%z(ke).ge.length) then
          k=1
          do while (dom(ib)%z(k).le.length)
            k=k+1
	    enddo
	    ktop=k-1
          do j=js-1,je+1; do i=is-1,ie+1
             dom(ib)%w(i,j,ktop) = dom(ib)%w(i,j,ktop-1)
	       dom(ib)%w(i,j,ktop+1) = dom(ib)%w(i,j,ktop-2)
!	       dom(ib)%w(i,j,ktop+1) = 0.0				!Pablo 03/19
	    enddo; enddo
	  end if
        if (dom(ib)%z(1).ge.length) then
	   do k=1,ke
          do j=js-1,je+1; do i=is-1,ie+1
	       dom(ib)%w(i,j,k) = 0.0
	    enddo; enddo;enddo
	  end if
	  end if
        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%bc_top.eq.3) then
             do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%w(i,j,ke+1+ly)=0.0
             end do; end do
           else if ((dom(ib)%bc_top.eq.4).or.
     &				(dom(ib)%bc_top.ge.61)) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%w(i,j,ke+1+ly)=0.0	
              end do; end do
           end if

	  	  if (trim(keyword).eq.'cavity') then	
		  do j=js-1,je+1 ; do i=is-1,ie+1
                 dom(ib)%w(i,j,ke+1+ly)=0.0	
          end do ; end do
	   	 endif         
	  end if


        end do
        if(diff_sch.ne.3) call boundcoef(3)
        end do


	if (LHRS) call bc_human_w     ! LHRS Alex & Yan

        end subroutine
```
