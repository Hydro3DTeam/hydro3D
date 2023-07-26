!#######################################################################
! SUBROUTINES: - iniflux
!              - correctoutflux
!#######################################################################
      SUBROUTINE iniflux
!-----------------------------------------------------------------------
!     Calculate The mass flow entering the domain from the inlet faces
!     Note: This flowmass is only calculated once from the initial flowfield
!#######################################################################
      use vars
      use multidata
      use mpi
      implicit none

      INTEGER :: i,j,k,ib,ispr,iepr,jspr,jepr,kspr,kepr
      DOUBLE PRECISION :: buffer_flomas

      MPI_FLT = MPI_DOUBLE_PRECISION

      flomas=0.d0

      DO ib=1,nbp
        ispr=pl+1; iepr=dom(ib)%ttc_i-pl
        jspr=pl+1; jepr=dom(ib)%ttc_j-pl
        kspr=pl+1; kepr=dom(ib)%ttc_k-pl
        
        ! Mass inflow on West Face:
        if(dom(ib)%iprev.lt.0) then
	   if(dom(ib)%bc_west.lt.61 .and. dom(ib)%bc_west.ne.4) then
	     IF(L_LSM) THEN
              do j=jspr,jepr ; do k=kspr,kepr
		  if(dom(ib)%phi(dom(ib)%isu-1,j,k) .ge. 0.0) then
                  flomas = flomas + dom(ib)%u(dom(ib)%isu-1,j,k) * dom(ib)%dy*dom(ib)%dz
		  end if
              end do ; end do

	     ELSE
              do j=jspr,jepr ; do k=kspr,kepr
                flomas = flomas + dom(ib)%u(dom(ib)%isu-1,j,k) * dom(ib)%dy*dom(ib)%dz
              end do ; end do
	     END IF
	   end if
        end if


        ! Mass inflow on South Face:
        if(dom(ib)%jprev.lt.0) then
	   if(dom(ib)%bc_south.lt.61 .and. dom(ib)%bc_south.ne.4) then
	     IF(L_LSM) THEN
              do k=kspr,kepr ; do i=ispr,iepr
		  if(dom(ib)%phi(i,dom(ib)%jsv-1,k) .ge. 0.0) then
                  flomas = flomas + dom(ib)%v(i,dom(ib)%jsv-1,k) * dom(ib)%dx*dom(ib)%dz
		  end if
              end do ; end do
	     
            ELSE
              do k=kspr,kepr ; do i=ispr,iepr
                flomas = flomas + dom(ib)%v(i,dom(ib)%jsv-1,k) * dom(ib)%dx*dom(ib)%dz
              end do ; end do
	     END IF
	   end if
        end if


        ! Mass inflow on Bottom Face:
        if(dom(ib)%kprev.lt.0) then
	   if(dom(ib)%bc_bottom.lt.61 .and. dom(ib)%bc_bottom.ne.4) then
	     IF(L_LSM) THEN
              do j=jspr,jepr ; do i=ispr,iepr
	    	  if(dom(ib)%phi(i,j,dom(ib)%ksw-1) .ge. 0.0) then
                  flomas = flomas + dom(ib)%w(i,j,dom(ib)%ksw-1) * dom(ib)%dx*dom(ib)%dy
	    	  end if
              end do ; end do

	     ELSE
              do j=jspr,jepr ; do i=ispr,iepr
                flomas = flomas + dom(ib)%w(i,j,dom(ib)%ksw-1) * dom(ib)%dx*dom(ib)%dy
              end do ; end do
	     END IF
	   endif
        end if

      END DO ! ib
      
     
      buffer_flomas = flomas
      call MPI_ALLREDUCE(buffer_flomas,flomas,1,MPI_FLT,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE iniflux
!#######################################################################
      SUBROUTINE correctoutflux
!-----------------------------------------------------------------------
!     Calculates the mass exiting the domain face (east,north,top) when the
!     boundary condition is not no-slip (4,61,62,63). Compare the total exiting
!     mass with the initial flow mass entering the domain. Correct the mass
!     difference on the exiting face (east,north,top) with a factor (fct).
!#######################################################################
      use vars
      use multidata
      use mpi
      implicit none

      INTEGER :: i,j,k,ib,ispr,iepr,jspr,jepr,kspr,kepr
      DOUBLE PRECISION :: fmout,fct,buffer_fmout

      MPI_FLT = MPI_DOUBLE_PRECISION

      fmout=0.d0

      DO ib=1,nbp
        ispr=pl+1; iepr=dom(ib)%ttc_i-pl
        jspr=pl+1; jepr=dom(ib)%ttc_j-pl
        kspr=pl+1; kepr=dom(ib)%ttc_k-pl
        
        ! Mass Outflow on the East Face:
        if(dom(ib)%inext.lt.0) then
	   if(dom(ib)%bc_east.lt.61 .and. dom(ib)%bc_east.ne.4) then
	     IF(L_LSM) THEN
              do j=jspr,jepr ; do k=kspr,kepr
		  if(dom(ib)%phi(dom(ib)%ieu+1,j,k) .ge. 0.0) then
                  fmout = fmout + dom(ib)%u(dom(ib)%ieu+1,j,k) * dom(ib)%dy*dom(ib)%dz
		  end if
              end do ; end do 
	     ELSE
              do j=jspr,jepr ; do k=kspr,kepr
                fmout = fmout + dom(ib)%u(dom(ib)%ieu+1,j,k)* dom(ib)%dy*dom(ib)%dz
              end do ; end do 
	     END IF
	   endif
        end if
        
        ! Mass Outflow on the North Face:
        if(dom(ib)%jnext.lt.0) then
	   if(dom(ib)%bc_north.lt.61 .and. dom(ib)%bc_north.ne.4) then
	     if(L_LSM) then
              do k=kspr,kepr ; do i=ispr,iepr
	         if(dom(ib)%phi(i,dom(ib)%jev+1,k) .ge. 0.0) then
                  fmout = fmout + dom(ib)%v(i,dom(ib)%jev+1,k) * dom(ib)%dx*dom(ib)%dz
	         end if
              end do ; end do
	     else
              do k=kspr,kepr ; do i=ispr,iepr
                fmout = fmout + dom(ib)%v(i,dom(ib)%jev+1,k) * dom(ib)%dx*dom(ib)%dz
              end do ; end do
	     end if
	   endif
        end if
        
        ! Mass Ourflow on the Top Face:
        if(dom(ib)%knext.lt.0) then
	   if(dom(ib)%bc_top.lt.61 .and. dom(ib)%bc_top.ne.4) then
	     if(L_LSM) then
              do j=jspr,jepr ; do i=ispr,iepr
		  if(dom(ib)%phi(i,j,dom(ib)%kew+1) .ge. 0.0) then
                  fmout = fmout + dom(ib)%w(i,j,dom(ib)%kew+1) * dom(ib)%dx*dom(ib)%dy
		  end if
              end do ; end do
	     else
              do j=jspr,jepr ; do i=ispr,iepr
                fmout = fmout + dom(ib)%w(i,j,dom(ib)%kew+1) * dom(ib)%dx*dom(ib)%dy
              end do ; end do
	     end if
	   endif
        end if

      END DO

      buffer_fmout = fmout
      call MPI_ALLREDUCE(buffer_fmout,fmout,1,MPI_FLT,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!Ratio of inflow and outflow
      fct=flomas/(fmout+1.E-30)


!Mass deficit
      Mdef=flomas-fmout

!Update outflow regarding the mass deficit at the inlet
      DO ib=1,nbp
        
        ! Mass Correction on East Face:
        if(dom(ib)%inext.lt.0) then
          if(dom(ib)%bc_east.lt.61 .and. dom(ib)%bc_east.ne.4) then
            IF(L_LSM) THEN							!PABLO 04/18
              ! U velocity:
              do j=dom(ib)%jsu,dom(ib)%jeu ; do k=dom(ib)%ksu,dom(ib)%keu
            	  if(dom(ib)%phi(dom(ib)%ieu+1,j,k) .ge. 0.0) then
                  dom(ib)%u(dom(ib)%ieu+1,j,k) = dom(ib)%u(dom(ib)%ieu+1,j,k) * fct
            	  end if
              end do ; end do

              ! V velocity:
              do j=dom(ib)%jsv,dom(ib)%jev ; do k=dom(ib)%ksv,dom(ib)%kev
            	  if(dom(ib)%phi(dom(ib)%iev+1,j,k) .ge. 0.0) then
                  dom(ib)%v(dom(ib)%iev+1,j,k) = dom(ib)%v(dom(ib)%iev+1,j,k) * fct
                end if 
              end do ; end do

              ! W velocity:
              do j=dom(ib)%jsw,dom(ib)%jew ; do k=dom(ib)%ksw,dom(ib)%kew
            	  if(dom(ib)%phi(dom(ib)%iew+1,j,k) .ge. 0.0) then
                  dom(ib)%w(dom(ib)%iew+1,j,k) = dom(ib)%w(dom(ib)%iew+1,j,k) * fct
                end if
              end do ; end do

            ELSE ! .not.L_LSM
              ! U velocity:
              do j=dom(ib)%jsu,dom(ib)%jeu ; do k=dom(ib)%ksu,dom(ib)%keu
                dom(ib)%u(dom(ib)%ieu+1,j,k) = dom(ib)%u(dom(ib)%ieu+1,j,k) * fct
              end do ; end do

              ! V velocity:
              do j=dom(ib)%jsv,dom(ib)%jev ; do k=dom(ib)%ksv,dom(ib)%kev
                dom(ib)%v(dom(ib)%iev+1,j,k) = dom(ib)%v(dom(ib)%iev+1,j,k) * fct
              end do ; end do

              ! W velocity:
              do j=dom(ib)%jsw,dom(ib)%jew ; do k=dom(ib)%ksw,dom(ib)%kew
                dom(ib)%w(dom(ib)%iew+1,j,k) = dom(ib)%w(dom(ib)%iew+1,j,k) * fct
              end do ; end do

            END IF
          end if
        end if ! inext
        


        ! Mass Correction on the North Face:
        if(dom(ib)%jnext.lt.0) then
          if(dom(ib)%bc_north.lt.61 .and. dom(ib)%bc_north.ne.4) then
            IF(L_LSM) THEN							!PABLO 04/18
              ! U Velocity:
              do i=dom(ib)%isu,dom(ib)%ieu ; do k=dom(ib)%ksu,dom(ib)%keu
                if(dom(ib)%phi(i,dom(ib)%jeu+1,k) .ge. 0.0) then
                  dom(ib)%u(i,dom(ib)%jeu+1,k) = dom(ib)%u(i,dom(ib)%jeu+1,k) * fct
                end if
              end do ; end do 

              ! V Velocity:
              do i=dom(ib)%isv,dom(ib)%iev ; do k=dom(ib)%ksv,dom(ib)%kev
                if(dom(ib)%phi(i,dom(ib)%jev+1,k) .ge. 0.0) then
                  dom(ib)%v(i,dom(ib)%jev+1,k) = dom(ib)%v(i,dom(ib)%jev+1,k) * fct
                end if
              end do ; end do 

              ! W Velocity:
              do i=dom(ib)%isw,dom(ib)%iew ; do k=dom(ib)%ksw,dom(ib)%kew
                if(dom(ib)%phi(i,dom(ib)%jew+1,k) .ge. 0.0) then
                  dom(ib)%w(i,dom(ib)%jew+1,k) = dom(ib)%w(i,dom(ib)%jew+1,k) * fct
                end if
              end do ; end do 

            ELSE
              ! U Velocity:
              do i=dom(ib)%isu,dom(ib)%ieu ; do k=dom(ib)%ksu,dom(ib)%keu
                dom(ib)%u(i,dom(ib)%jeu+1,k) = dom(ib)%u(i,dom(ib)%jeu+1,k) * fct
              end do ; end do 

              ! V Velocity:
              do i=dom(ib)%isv,dom(ib)%iev ; do k=dom(ib)%ksv,dom(ib)%kev
                dom(ib)%v(i,dom(ib)%jev+1,k) = dom(ib)%v(i,dom(ib)%jev+1,k) * fct
              end do ; end do 

              ! W Velocity:
              do i=dom(ib)%isw,dom(ib)%iew ; do k=dom(ib)%ksw,dom(ib)%kew
                dom(ib)%w(i,dom(ib)%jew+1,k) = dom(ib)%w(i,dom(ib)%jew+1,k) * fct
              end do ; end do 
            END IF
          end if
        end if ! jnext



        ! Mass Correction on the Top Face:
        if(dom(ib)%knext.lt.0) then
          if(dom(ib)%bc_top.lt.61 .and. dom(ib)%bc_top.ne.4) then
            IF(L_LSM) THEN							!PABLO 04/18
              ! U Velocity:
              do i=dom(ib)%isu,dom(ib)%ieu ; do j=dom(ib)%jsu,dom(ib)%jeu
            	  if(dom(ib)%phi(i,j,dom(ib)%keu+1) .ge. 0.0) then
                  dom(ib)%u(i,j,dom(ib)%keu+1) = dom(ib)%u(i,j,dom(ib)%keu+1) * fct
                end if
              end do ; end do
              
              ! V Velocity:
              do i=dom(ib)%isv,dom(ib)%iev ; do j=dom(ib)%jsv,dom(ib)%jev
            	  if(dom(ib)%phi(i,j,dom(ib)%kev+1) .ge. 0.0) then
                  dom(ib)%v(i,j,dom(ib)%kev+1) = dom(ib)%v(i,j,dom(ib)%kev+1) * fct
            	  end if
              end do ; end do

              ! W Velocity:
              do i=dom(ib)%isw,dom(ib)%iew ; do j=dom(ib)%jsw,dom(ib)%jew
            	  if(dom(ib)%phi(i,j,dom(ib)%kew+1) .ge. 0.0) then
                  dom(ib)%w(i,j,dom(ib)%kew+1) = dom(ib)%w(i,j,dom(ib)%kew+1) * fct
            	  end if
              end do ; end do
              
            ELSE
              ! U Velocity:
              do i=dom(ib)%isu,dom(ib)%ieu ; do j=dom(ib)%jsu,dom(ib)%jeu
                dom(ib)%u(i,j,dom(ib)%keu+1) = dom(ib)%u(i,j,dom(ib)%keu+1) * fct
              end do ; end do
              
              ! V Velocity:
              do i=dom(ib)%isv,dom(ib)%iev ; do j=dom(ib)%jsv,dom(ib)%jev
                dom(ib)%v(i,j,dom(ib)%kev+1) = dom(ib)%v(i,j,dom(ib)%kev+1) * fct
              end do ; end do

              ! W Velocity:
              do i=dom(ib)%isw,dom(ib)%iew ; do j=dom(ib)%jsw,dom(ib)%jew
                dom(ib)%w(i,j,dom(ib)%kew+1) = dom(ib)%w(i,j,dom(ib)%kew+1) * fct
              end do ; end do

            END IF
          endif
        end if

      END DO ! ib

      RETURN
      END SUBROUTINE correctoutflux
!#######################################################################
      SUBROUTINE pressure_forcing
!-----------------------------------------------------------------------
!     Only applies to the periodic boundary conditions, calculates pressure forcing
!     required to maintain the flow mass conservation. The extra forcing is then 
!     distribute in the diffusion spatial term. 
!     Note: Using equation Q(t) = Q(t-1) + 0.3*(U(0)+U(t-1)-2*U(t))/dt
!#######################################################################
      use vars
      use mpi
      use multidata
      implicit none

      INTEGER :: i,j,k,ib,isp,jsp,ksp,jspr,jepr,kspr,kepr,ispr,iepr
      DOUBLE PRECISION ::  qstpp,fakfor,flwsum_loc,flwsum,A_loc,A
      DOUBLE PRECISION ::  qstpp_y,flwsum_loc_y,flwsum_y,A_loc_y,A_y
      DOUBLE PRECISION ::  qstpp_z,flwsum_loc_z,flwsum_z,A_loc_z,A_z
      DOUBLE PRECISION,dimension(2) :: sndbuffer,recbuffer
      DOUBLE PRECISION,dimension(2) :: sndbuffer_y,recbuffer_y
      DOUBLE PRECISION,dimension(2) :: sndbuffer_z,recbuffer_z

      MPI_FLT = MPI_DOUBLE_PRECISION

      flwsum_loc = 0.0    ;  A_loc = 0.0
      flwsum_loc_y = 0.0  ;  A_loc_y = 0.0
      flwsum_loc_z = 0.0  ;  A_loc_z = 0.0

      DO ib=1,nbp
        ispr=pl+1; iepr=dom(ib)%ttc_i-pl
        jspr=pl+1; jepr=dom(ib)%ttc_j-pl
        kspr=pl+1; kepr=dom(ib)%ttc_k-pl

	 if(pressureforce) then
          if(dom(ib)%iprev.lt.0) then						!west
            isp=dom(ib)%isu-1
	     IF(L_LSM) THEN
              do j=jspr,jepr ; do k=kspr,kepr
		  if(dom(ib)%phi(isp,j,k) .ge. 0.d0) then
                  flwsum_loc = flwsum_loc + dom(ib)%u(isp,j,k) * dom(ib)%dy*dom(ib)%dz &
     &                       * (dom(ib)%dens(isp,j,k)/densl)
                  A_loc = A_loc + dom(ib)%dy*dom(ib)%dz * (dom(ib)%dens(isp,j,k)/densl)
		  end if
              end do ; end do

	     ELSE IF(L_LSMbase) THEN
              do j=jspr,jepr ; do k=kspr,kepr
		  if(dom(ib)%z(k) .le. length) then
                  flwsum_loc = flwsum_loc + dom(ib)%u(isp,j,k) * dom(ib)%dy*dom(ib)%dz
                  A_loc = A_loc + dom(ib)%dy*dom(ib)%dz
		  end if
              end do ; end do

	     ELSE
              do j=jspr,jepr ; do k=kspr,kepr
                flwsum_loc = flwsum_loc + dom(ib)%u(isp,j,k) * dom(ib)%dy*dom(ib)%dz
                A_loc = A_loc + dom(ib)%dy*dom(ib)%dz
              end do ; end do
	     END IF
           end if ! iprev
	 end if ! pressureforce
          


        ! 
	 if(pressureforce_y) then
          if(dom(ib)%jprev.lt.0) then
            jsp=dom(ib)%jsu-1
	     IF(L_LSM) THEN
              do i=ispr,iepr ; do k=kspr,kepr
		  if(dom(ib)%phi(i,jsp,k) .ge. 0.d0) then
                  flwsum_loc_y = flwsum_loc_y + dom(ib)%v(i,jsp,k) * dom(ib)%dx*dom(ib)%dz &
     &                         * (dom(ib)%dens(i,jsp,k)/densl)
                  A_loc_y = A_loc_y + dom(ib)%dx*dom(ib)%dz * (dom(ib)%dens(i,jsp,k)/densl)
		  end if
              end do ; end do
	     
            ELSE IF(L_LSMbase) THEN
              do i=ispr,iepr ; do k=kspr,kepr
		  if(dom(ib)%z(k) .le. length) then
                  flwsum_loc_y = flwsum_loc_y + dom(ib)%v(i,jsp,k) * dom(ib)%dx*dom(ib)%dz
                  A_loc_y = A_loc_y + dom(ib)%dx*dom(ib)%dz
		  end if
              end do ; end do

	     ELSE
              do i=ispr,iepr ; do k=kspr,kepr
                flwsum_loc_y = flwsum_loc_y + dom(ib)%v(i,jsp,k) * dom(ib)%dx*dom(ib)%dz
                A_loc_y = A_loc_y + dom(ib)%dx*dom(ib)%dz
              end do ; end do
	     END IF
          end if ! jprev
	 end if ! pressureforce_y




	 if(pressureforce_z) then
          if(dom(ib)%kprev.lt.0) then
            ksp=dom(ib)%ksu-1
	     IF(L_LSM) THEN
              do i=ispr,iepr ; do j=jspr,jepr
		  if(dom(ib)%phi(i,j,ksp) .ge. 0.d0) then
                  flwsum_loc_z = flwsum_loc_z + dom(ib)%w(i,j,ksp) * dom(ib)%dx*dom(ib)%dy &
     &                         * (dom(ib)%dens(i,j,ksp)/densl)
                  A_loc_z = A_loc_z + dom(ib)%dx*dom(ib)%dy * (dom(ib)%dens(i,j,ksp)/densl)
		  end if
              end do ; end do

	     ELSE IF (L_LSMbase) THEN
              do i=ispr,iepr ; do j=jspr,jepr
		  if(dom(ib)%z(k) .le. length) then
                  flwsum_loc_z = flwsum_loc_z + dom(ib)%w(i,j,ksp) * dom(ib)%dx*dom(ib)%dy
                  A_loc_z = A_loc_z + dom(ib)%dx*dom(ib)%dy
		  end if
              end do ; end do
	     ELSE
              do i=ispr,iepr ; do j=jspr,jepr
                flwsum_loc_z = flwsum_loc_z + dom(ib)%w(i,j,ksp) * dom(ib)%dx*dom(ib)%dy
                A_loc_z = A_loc_z + dom(ib)%dx*dom(ib)%dy
              end do ; end do
	     END IF
          end if ! kprev
	 end if ! pressureforce_z

      END DO ! ib

!.....sum massflux and area over all processors
      if(pressureforce) then
        sndbuffer(1)   = flwsum_loc    ;  sndbuffer(2)   = A_loc
        call MPI_ALLREDUCE(sndbuffer,recbuffer,2,MPI_FLT,MPI_SUM,MPI_COMM_WORLD,ierr)
        flwsum   = recbuffer(1)  ;   A        = recbuffer(2)
      end if

      if(pressureforce_y) then
        sndbuffer_y(1) = flwsum_loc_y  ;  sndbuffer_y(2) = A_loc_y
        call MPI_ALLREDUCE(sndbuffer_y,recbuffer_y,2,MPI_FLT,MPI_SUM,MPI_COMM_WORLD,ierr)
        flwsum_y = recbuffer_y(1)  ;   A_y    = recbuffer_y(2)
      end if

      if(pressureforce_z) then
        sndbuffer_z(1) = flwsum_loc_z  ;  sndbuffer_z(2) = A_loc_z
        call MPI_ALLREDUCE(sndbuffer_z,recbuffer_z,2,MPI_FLT,MPI_SUM,MPI_COMM_WORLD,ierr)
        flwsum_z = recbuffer_z(1)  ;   A_z    = recbuffer_z(2)
      end if

! --- Calculate and store forcing term --------------------------------
      fakfor  = 0.3
      if(pressureforce)    qstpp   = flwsum/A
      if(pressureforce_y)  qstpp_y = flwsum_y/A_y
      if(pressureforce_z)  qstpp_z = flwsum_z/A_z

! Whatever the flow is lost, is added back again
! 	  Q(t)=Q(t-1)+0.3?  *(U(0) +U(t-1)-2*U(t))/dt
      if(pressureforce) then
        forcn = forcn + fakfor*(qzero+qstpn-2.0*qstpp)/dt
        qstpn = qstpp
      end if

      if(pressureforce_y) then
        forcn_y = forcn_y + fakfor*(0.0+qstpn_y-2.0*qstpp_y)/dt
        qstpn_y = qstpp_y
      end if

      if(pressureforce_z) then
        forcn_z = forcn_z + fakfor*(0.0+qstpn_z-2.0*qstpp_z)/dt
        qstpn_z = qstpp_z
      end if

!Write out values in the 3 directions:
      if(myrank.eq.0) then
        if(pressureforce) write (numfile1,'(4F20.12)') ctime,forcn,qstpn,flwsum
        if(pressureforce_y) write (numfile4,'(4F20.12)') ctime,forcn_y,qstpn_y,flwsum_y
        if(pressureforce_z) write (numfile5,'(4F20.12)') ctime,forcn_z,qstpn_z,flwsum_z
      end if
      
      RETURN
      END SUBROUTINE pressure_forcing
