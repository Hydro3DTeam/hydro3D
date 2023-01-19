# boundv

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
****This subroutine allows you to prescribe v velocity boundary condition on all on the computational domain faces (West-East, South-North, Bottom-Top).\
\
**User Likeness to alter the file:** \
****Often\


```
!#############################################################################
        subroutine boundv
!#############################################################################
        use vars
        use multidata
	  use mpi
	  use imb
        implicit none
        integer :: i,j,k,ib,ly,kk,jj
        integer :: is,ie,js,je,ks,ke,ktop
	  integer :: strlen
	  character (LEN=29):: filename,fileSEM
	  character (LEN=4) :: domain
	  character (LEN=5) :: name_end
	  double precision :: Fwallv,dxx,dzz,dummy,dyy

        if (PERIODIC) call exchange_bc(2,pl_ex)
 
        do ly=0,pl_ex
        do ib=1,nbp
	     dxx=dom(ib)%dx
	     dyy=dom(ib)%dy
	     dzz=dom(ib)%dz
           is=dom(ib)%isv; ie=dom(ib)%iev
           js=dom(ib)%jsv; je=dom(ib)%jev
           ks=dom(ib)%ksv; ke=dom(ib)%kev

! Boundary Conditions for V - i.e. shear at boundaries
!...............................................................................
!=== West ===>   4=wall   ..    1=Inflow
!...............................................................................
        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%bc_west.eq.4) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%v(is-1-ly,j,k)= -dom(ib)%v(is+ly,j,k)	
              end do; end do

           else if (dom(ib)%bc_west.eq.1 .or. dom(ib)%bc_west.eq.12
     & .or. dom(ib)%bc_west.eq.14 .or. dom(ib)%bc_west.eq.15
     & .or. dom(ib)%bc_west.eq.17 .or. dom(ib)%bc_west.eq.31)then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%v(is-1-ly,j,k)=0.0
              end do; end do

	     else if (dom(ib)%bc_west.ge.61) then					!Wall functions
		if (ly.eq.0) then
	   		if (dom(ib)%bc_west.lt.63) 
     &			call log_law(1,dom(ib)%bc_west)
			if (dom(ib)%bc_west.ge.63) 
     &			call wall_function(1,dom(ib)%bc_west)
		      do k=ks-1,ke+1; do j=js-1,je+1
				Fwallv = dom(ib)%tauww2(j,k)
     &					*dom(ib)%v(is,j,k)/dxx			!*Acell/Vcell
				Fwallv = sign(Fwallv,-dom(ib)%v(is,j,k))
				dom(ib)%v(is,j,k)=dom(ib)%v(is,j,k) 
     &						+Fwallv*alfapr*dt
                 		dom(ib)%v(is-1,j,k)= dom(ib)%v(is,j,k)	
		      end do; end do
		else
		      do k=ks-1,ke+1; do j=js-1,je+1
                 		dom(ib)%v(is-1-ly,j,k)= dom(ib)%v(is,j,k)	
		      end do; end do
		endif

          else if (dom(ib)%bc_west.eq.7) then					!Rading slices
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
		      do k=ks-1,ke+1; do j=js-1,je+1
				read(405,*)dom(ib)%v(is-1-ly,j,k)
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
              do k=ks-1,ke+1; do j=js-1,je+1
			read(405,*)dom(ib)%v(is-1-ly,j,k)
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
	   	   read(405,*)dummy,dom(ib)%v(is-1-ly,j,k),dummy
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
                 dom(ib)%v(ie+1+ly,j,k)= -dom(ib)%v(ie-ly,j,k)	
              end do; end do

           else if(dom(ib)%bc_east.eq.2 .or.dom(ib)%bc_east.eq.21) then
              if (alfabc.eq.1) then
                 do k=ks-1,ke+1; do j=js-1,je+1

                 if (L_LSM) then

			 if (dom(ib)%phi(ie,j,k) .ge. 0.0) then
                     if (dom(ib)%bc_east.eq.2) then
                       dom(ib)%v(ie+1+ly,j,k)= dom(ib)%v(ie,j,k)
                     else
                       dom(ib)%v(ie+1+ly,j,k)= dom(ib)%voo(ie+1+ly,j,k)
     & -dt*alfapr*(dom(ib)%voo(ie+1+ly,j,k)-
     & dom(ib)%voo(ie,j,k))/dom(ib)%dx
                     end if
			 else
			   dom(ib)%v(ie+1+ly,j,k)= dom(ib)%v(ie,j,k) 
			 endif

		     else									!no LSM

                     if (dom(ib)%bc_east.eq.2) then
                       dom(ib)%v(ie+1+ly,j,k)= dom(ib)%v(ie,j,k)
                     else
                       dom(ib)%v(ie+1+ly,j,k)= dom(ib)%voo(ie+1+ly,j,k)
     & -dt*alfapr*(dom(ib)%voo(ie+1+ly,j,k)-
     & dom(ib)%voo(ie,j,k))/dom(ib)%dx
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
				Fwallv = dom(ib)%tauwe2(j,k)
     &					*dom(ib)%v(ie,j,k)/dxx			!*Acell/Vcell
				Fwallv = sign(Fwallv,-dom(ib)%v(ie,j,k))
				dom(ib)%v(ie,j,k)=dom(ib)%v(ie,j,k) 
     &						+Fwallv*alfapr*dt
                 		dom(ib)%v(ie+1,j,k)= dom(ib)%v(ie,j,k)	
		      end do; end do
		else
		      do k=ks-1,ke+1; do j=js-1,je+1
                 		dom(ib)%v(ie+1+ly,j,k)= dom(ib)%v(ie,j,k)	
		      end do; end do
		endif
           end if
        end if
!...............................................................................
!=== South ===>   4=wall   ..    3=Symmetry
!...............................................................................
        if (dom(ib)%jprev.lt.0) then 
           if (dom(ib)%bc_south.eq.3) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%v(i,js-1-ly,k)=0.0
              end do; end do

	else if ((dom(ib)%bc_south.eq.4).or.(dom(ib)%bc_south.ge.61)) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%v(i,js-1-ly,k)=0.0	
              end do; end do
           end if
        end if
!...............................................................................
!=== North ===>   4=wall   ..    3=Symmetry
!...............................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%bc_north.eq.3) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%v(i,je+1+ly,k)=0.0
              end do; end do

	else if ((dom(ib)%bc_north.eq.4).or.(dom(ib)%bc_north.ge.61)) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%v(i,je+1+ly,k)=0.0	
              end do; end do
           end if
        end if
!...............................................................................
!=== Bottom ===> ..   4=wall ..   3=Symmetry
!...............................................................................
        if (dom(ib)%kprev.lt.0) then
          if (dom(ib)%bc_bottom.eq.4) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%v(i,j,ks-1-ly)= -dom(ib)%v(i,j,ks+ly)	
              end do; end do

           else if (dom(ib)%bc_bottom.eq.1) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%v(i,j,ks-1-ly)= 0.0
              end do; end do

           else if (dom(ib)%bc_bottom.eq.3) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%v(i,j,ks-1-ly)= dom(ib)%v(i,j,ks+ly)
              end do; end do
	     else if (dom(ib)%bc_bottom.ge.61) then				!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_bottom.lt.63) then
     				call log_law(5,dom(ib)%bc_bottom)
			endif
			if (dom(ib)%bc_bottom.ge.63) 
     &			call wall_function(5,dom(ib)%bc_bottom)
		      do j=js-1,je+1; do i=is-1,ie+1
				Fwallv = dom(ib)%tauwb2(i,j)
     &					*dom(ib)%v(i,j,ks)/dzz			!*Acell/Vcell
				Fwallv = sign(Fwallv,-dom(ib)%v(i,j,ks))
				dom(ib)%v(i,j,ks)=dom(ib)%v(i,j,ks) 
     &						+Fwallv*alfapr*dt
                 		dom(ib)%v(i,j,ks-1)= dom(ib)%v(i,j,ks)	
		      end do; end do
		else
		      do j=js-1,je+1; do i=is-1,ie+1
                 		dom(ib)%v(i,j,ks-1-ly)= dom(ib)%v(i,j,ks)	
		      end do; end do
		endif
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
	    end do
	    ktop=k-1
          do j=js-1,je+1; do i=is-1,ie+1
            dom(ib)%v(i,j,ktop) = dom(ib)%v(i,j,ktop-1)
            dom(ib)%v(i,j,ktop+1) = dom(ib)%v(i,j,ktop-2)
	    enddo; enddo
	  end if
        if (dom(ib)%z(1).ge.length) then				! PABLO 03/19
	    do k=1,ke
          do j=js-1,je+1; do i=is-1,ie+1
            dom(ib)%v(i,j,k) = 0.0
	    enddo; enddo;enddo
	  end if
	  end if

        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%bc_top.eq.4) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%v(i,j,ke+1+ly)= -dom(ib)%v(i,j,ke-ly)	
              end do; end do
           else if (dom(ib)%bc_top.eq.3) then
             do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%v(i,j,ke+1+ly)= dom(ib)%v(i,j,ke-ly)
             end do; end do
	     else if (dom(ib)%bc_top.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_top.lt.63) 
     &			call log_law(6,dom(ib)%bc_top)
			if (dom(ib)%bc_top.ge.63) 
     &			call wall_function(6,dom(ib)%bc_top)

		      do j=js-1,je+1; do i=is-1,ie+1
				Fwallv = dom(ib)%tauwt2(i,j)
     &					*dom(ib)%v(i,j,ke)/dzz			!*Acell/Vcell
				Fwallv = sign(Fwallv,-dom(ib)%v(i,j,ke))
				dom(ib)%v(i,j,ke)=dom(ib)%v(i,j,ke) 
     &						+Fwallv*alfapr*dt
                 		dom(ib)%v(i,j,ke+1)= dom(ib)%v(i,j,ke)	
		      end do; end do
		else
		      do j=js-1,je+1; do i=is-1,ie+1
                 		dom(ib)%v(i,j,ke+1+ly)= dom(ib)%v(i,j,ke)	
		      end do; end do
		endif
           end if !bc_top.eq.4

	     if (trim(keyword).eq.'cavity') then	
		  do j=js-1,je+1 ; do i=is-1,ie+1
                dom(ib)%v(i,j,ke+1+ly) = dom(ib)%v(i,j,ke-ly)
          end do ; end do
	   	 endif           
	  end if



        end do
        if(diff_sch.ne.3) call boundcoef(2)
        end do

	if (LHRS) call bc_human_v     ! LHRS Alex & Yan

        end subroutine
```
