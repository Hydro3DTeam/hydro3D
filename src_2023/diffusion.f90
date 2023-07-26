!#######################################################################
! SUBROUTINES: - diffusion
!              - rungek_diff2nd
!              - rungek_diff4th
!#######################################################################
      SUBROUTINE diffusion
!-----------------------------------------------------------------------
!     Numerical method that calculates the diffusion spatial term of
!     the Navier-Stokes equations.
!     Note 1: diffusion scheme for 1.Implicit Euler, 2.Crank-Nicholson
!     Note 2: LSM can not be used with this diffusion scheme
!     Note 3: This scheme only uses the SIP solver can not use MGsolver
!#######################################################################
      use vars
      use mpi
      use multidata
      implicit none

      INTEGER :: i,j,k,ib
      INTEGER :: is,ie,js,je,ks,ke
      DOUBLE PRECISION :: dxx,dyy,dzz

! ------ Compute coefficients for diffusion terms

      if(diff_sch.eq.1) then
         fac=dt
      else if(diff_sch.eq.2) then
         fac=dt/2.0
      end if


      do ib=1,nbp
         dom(ib)%ap = 0.0
         dom(ib)%su = 0.0
      end do

      call boundu

      do ib=1,nbp
        dxx=dom(ib)%dx*dom(ib)%dx
        dyy=dom(ib)%dy*dom(ib)%dy
        dzz=dom(ib)%dz*dom(ib)%dz

        is=dom(ib)%isu; ie=dom(ib)%ieu
        js=dom(ib)%jsu; je=dom(ib)%jeu
        ks=dom(ib)%ksu; ke=dom(ib)%keu

        ! Implicit Euler:
        IF(diff_sch.eq.1) THEN 
          do k=ks,ke ; do i=is,ie ; do j=js,je
          dom(ib)%su(i,j,k) = dom(ib)%ustar(i,j,k) * dom(ib)%su(i,j,k)
          end do ; end do ; end do

        ! Crank-Nicholson
        ELSE IF(diff_sch.eq.2) THEN
          do k=ks,ke ; do i=is,ie ; do j=js,je

            dom(ib)%aw(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dxx
            dom(ib)%ae(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dxx
            dom(ib)%an(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dyy
            dom(ib)%as(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dyy
            dom(ib)%at(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dzz
            dom(ib)%ab(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dzz

            if (dom(ib)%iprev.lt.0)  dom(ib)%aw(is,j,k)=0.0
            if (dom(ib)%inext.lt.0)  dom(ib)%ae(ie,j,k)=0.0
            if (dom(ib)%jprev.lt.0)  dom(ib)%as(i,js,k)=0.0
            if (dom(ib)%jnext.lt.0)  dom(ib)%an(i,je,k)=0.0
            if (dom(ib)%kprev.lt.0)  dom(ib)%ab(i,j,ks)=0.0
            if (dom(ib)%knext.lt.0)  dom(ib)%at(i,j,ke)=0.0

            dom(ib)%ap(i,j,k) = dom(ib)%ap(i,j,k)-1.0*( &
     &  dom(ib)%aw(i,j,k)+dom(ib)%ae(i,j,k)+ &
     &  dom(ib)%as(i,j,k)+dom(ib)%an(i,j,k)+ &
     &  dom(ib)%ab(i,j,k)+dom(ib)%at(i,j,k))+1.0

            dom(ib)%su(i,j,k)=((1.0- &
     &dom(ib)%ap(i,j,k))*dom(ib)%u(i,j,k)- &
     &dom(ib)%aw(i,j,k)*dom(ib)%u(i-1,j,k)- &
     &dom(ib)%ae(i,j,k)*dom(ib)%u(i+1,j,k)- &
     &dom(ib)%as(i,j,k)*dom(ib)%u(i,j-1,k)- &
     &dom(ib)%an(i,j,k)*dom(ib)%u(i,j+1,k)- &
     &dom(ib)%ab(i,j,k)*dom(ib)%u(i,j,k-1)- &
     &dom(ib)%at(i,j,k)*dom(ib)%u(i,j,k+1)+ &
     &dom(ib)%ustar(i,j,k))+dom(ib)%su(i,j,k)
        
          end do ;  end do ; end do
        END IF ! diff_sch

        dom(ib)%ustar(:,:,:)=0.0

      END DO ! ib

      if(pressureforce) dom(ib)%su(:,:,:) = dom(ib)%su(:,:,:) + forcn * dt
      
      
!                 if (LENERGY) dom(ib)%su(i,j,k)= &
!     & dom(ib)%su(i,j,k)+dt*grx*(1.d0-0.5*beta* &
!     & (dom(ib)%T(i+1,j,k)+dom(ib)%T(i,j,k)) )

      call sipsol(11)

      do ib=1,nbp
        dom(ib)%ap = 0.0
        dom(ib)%su = 0.0
      end do

      call boundv

      do ib=1,nbp
        dxx=dom(ib)%dx*dom(ib)%dx
        dyy=dom(ib)%dy*dom(ib)%dy
        dzz=dom(ib)%dz*dom(ib)%dz

        is=dom(ib)%isv; ie=dom(ib)%iev
        js=dom(ib)%jsv; je=dom(ib)%jev
        ks=dom(ib)%ksv; ke=dom(ib)%kev

        ! Implicit Euler:
        IF(diff_sch.eq.1) THEN
          do k=ks,ke ; do i=is,ie ; do j=js,je
            dom(ib)%su(i,j,k ) = dom(ib)%vstar(i,j,k) + dom(ib)%su(i,j,k)
          end do ; end do ; end do

        ! Crank-Nicholson
        ELSE IF(diff_sch.eq.2) THEN
          do k=ks,ke ; do i=is,ie ; do j=js,je

            dom(ib)%aw(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dxx
            dom(ib)%ae(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dxx
            dom(ib)%an(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dyy
            dom(ib)%as(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dyy
            dom(ib)%at(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dzz
            dom(ib)%ab(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dzz

            if (dom(ib)%iprev.lt.0)  dom(ib)%aw(is,j,k)=0.0
            if (dom(ib)%inext.lt.0)  dom(ib)%ae(ie,j,k)=0.0
            if (dom(ib)%jprev.lt.0)  dom(ib)%as(i,js,k)=0.0
            if (dom(ib)%jnext.lt.0)  dom(ib)%an(i,je,k)=0.0
            if (dom(ib)%kprev.lt.0)  dom(ib)%ab(i,j,ks)=0.0
            if (dom(ib)%knext.lt.0)  dom(ib)%at(i,j,ke)=0.0

            dom(ib)%ap(i,j,k) = dom(ib)%ap(i,j,k)-1.0*( &
     &  dom(ib)%aw(i,j,k)+dom(ib)%ae(i,j,k)+ &
     &  dom(ib)%as(i,j,k)+dom(ib)%an(i,j,k)+ &
     &  dom(ib)%ab(i,j,k)+dom(ib)%at(i,j,k))+1.0

            dom(ib)%su(i,j,k)=((1.0- &
     &dom(ib)%ap(i,j,k))*dom(ib)%v(i,j,k)- &
     &dom(ib)%aw(i,j,k)*dom(ib)%v(i-1,j,k)- &
     &dom(ib)%ae(i,j,k)*dom(ib)%v(i+1,j,k)- &
     &dom(ib)%as(i,j,k)*dom(ib)%v(i,j-1,k)- &
     &dom(ib)%an(i,j,k)*dom(ib)%v(i,j+1,k)- &
     &dom(ib)%ab(i,j,k)*dom(ib)%v(i,j,k-1)- &
     &dom(ib)%at(i,j,k)*dom(ib)%v(i,j,k+1)+ &
     &dom(ib)%vstar(i,j,k))+dom(ib)%su(i,j,k)

          end do ; end do ; end do
        END IF ! diff_sch

        dom(ib)%vstar(:,:,:)=0.d0
      END DO ! ib

      if(pressureforce_y) dom(ib)%su(:,:,:) = dom(ib)%su(:,:,:) + forcn_y * dt

!                 if (LENERGY) dom(ib)%su(i,j,k)= &
!     & dom(ib)%su(i,j,k)+dt*gry*(1.d0-0.5*beta* &
!     & (dom(ib)%T(i,j+1,k)+dom(ib)%T(i,j,k)) )

      call sipsol(22)

      do ib=1,nbp
        dom(ib)%ap = 0.0
        dom(ib)%su = 0.0
      end do

      call boundw

      do ib=1,nbp
        dxx=dom(ib)%dx*dom(ib)%dx
        dyy=dom(ib)%dy*dom(ib)%dy
        dzz=dom(ib)%dz*dom(ib)%dz

        is=dom(ib)%isw; ie=dom(ib)%iew
        js=dom(ib)%jsw; je=dom(ib)%jew
        ks=dom(ib)%ksw; ke=dom(ib)%kew
        
        ! Implicit Euler:
        if(diff_sch.eq.1) THEN
          do k=ks,ke ; do i=is,ie ; do j=js,je
            dom(ib)%su(i,j,k) = dom(ib)%wstar(i,j,k) + dom(ib)%su(i,j,k)
          end do ; end do ; end do 
        
        ! Crank-Nicholson
        ELSE IF(diff_sch.eq.2) THEN
          do k=ks,ke ; do i=is,ie ; do j=js,je

            dom(ib)%aw(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dxx
            dom(ib)%ae(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dxx
            dom(ib)%an(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dyy
            dom(ib)%as(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dyy
            dom(ib)%at(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dzz
            dom(ib)%ab(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dzz

            if (dom(ib)%iprev.lt.0)  dom(ib)%aw(is,j,k)=0.0
            if (dom(ib)%inext.lt.0)  dom(ib)%ae(ie,j,k)=0.0
            if (dom(ib)%jprev.lt.0)  dom(ib)%as(i,js,k)=0.0
            if (dom(ib)%jnext.lt.0)  dom(ib)%an(i,je,k)=0.0
            if (dom(ib)%kprev.lt.0)  dom(ib)%ab(i,j,ks)=0.0
            if (dom(ib)%knext.lt.0)  dom(ib)%at(i,j,ke)=0.0

            dom(ib)%ap(i,j,k) = dom(ib)%ap(i,j,k)-1.0*( &
     &  dom(ib)%aw(i,j,k)+dom(ib)%ae(i,j,k)+ &
     &  dom(ib)%as(i,j,k)+dom(ib)%an(i,j,k)+ &
     &  dom(ib)%ab(i,j,k)+dom(ib)%at(i,j,k))+1.0

            dom(ib)%su(i,j,k)=((1.0- &
     &dom(ib)%ap(i,j,k))*dom(ib)%w(i,j,k)- &
     &dom(ib)%aw(i,j,k)*dom(ib)%w(i-1,j,k)- &
     &dom(ib)%ae(i,j,k)*dom(ib)%w(i+1,j,k)- &
     &dom(ib)%as(i,j,k)*dom(ib)%w(i,j-1,k)- &
     &dom(ib)%an(i,j,k)*dom(ib)%w(i,j+1,k)- &
     &dom(ib)%ab(i,j,k)*dom(ib)%w(i,j,k-1)- &
     &dom(ib)%at(i,j,k)*dom(ib)%w(i,j,k+1)+ &
     &dom(ib)%wstar(i,j,k))+dom(ib)%su(i,j,k)
          end do ; end do ; end do
        END IF ! diff_sch

        dom(ib)%wstar(:,:,:)=0.d0
      END DO ! ib

!                 if (LENERGY) dom(ib)%su(i,j,k)= &
!     &dom(ib)%su(i,j,k)+dt*grz*(1.d0-0.5*beta* &
!     &(dom(ib)%T(i,j,k+1)+dom(ib)%T(i,j,k)) )

      if(pressureforce_z) dom(ib)%su(:,:,:) = dom(ib)%su(:,:,:) + forcn_z * dt

      call sipsol(33)

      if (LROUGH) call rough_velocity

      RETURN
      END SUBROUTINE diffusion
!#######################################################################
      SUBROUTINE rungek_diff2nd
!-----------------------------------------------------------------------
!     As part of the fractional step method we add the diffusion spatial term
!     onto the free-divergence field (ustar,vstar,wstar)
!     Note 1: LSM, Periodic, Scalar, Rough_bed compliant
!#######################################################################
      use vars
      use mpi
      use multidata
      implicit none

      INTEGER :: i,j,k,ib
      INTEGER :: is,ie,js,je,ks,ke
      DOUBLE PRECISION :: dxx,dyy,dzz,visc,diff

      fac=1.d0

      do ib=1,nbp
        dxx=dom(ib)%dx*dom(ib)%dx
        dyy=dom(ib)%dy*dom(ib)%dy
        dzz=dom(ib)%dz*dom(ib)%dz
        

        ! ADDING THE DIFFUSION TO USTAR:
        is=dom(ib)%isu; ie=dom(ib)%ieu
        js=dom(ib)%jsu; je=dom(ib)%jeu
        ks=dom(ib)%ksu; ke=dom(ib)%keu
        do k=ks,ke ; do i=is,ie ; do j=js,je
                 
          visc=dom(ib)%vis(i+1,j,k)

          dom(ib)%aw(i,j,k)=-visc/dxx
          dom(ib)%ae(i,j,k)=-visc/dxx
          dom(ib)%an(i,j,k)=-visc/dyy
          dom(ib)%as(i,j,k)=-visc/dyy
          dom(ib)%at(i,j,k)=-visc/dzz
          dom(ib)%ab(i,j,k)=-visc/dzz

          dom(ib)%ap(i,j,k) = -1.0*( &
     &dom(ib)%aw(i,j,k)+dom(ib)%ae(i,j,k)+ &
     &dom(ib)%as(i,j,k)+dom(ib)%an(i,j,k)+ &
     &dom(ib)%ab(i,j,k)+dom(ib)%at(i,j,k))

          diff=-(dom(ib)%ap(i,j,k)*dom(ib)%u(i,j,k)+ &
     & dom(ib)%an(i,j,k)*dom(ib)%u(i,j+1,k) + &
     & dom(ib)%as(i,j,k)*dom(ib)%u(i,j-1,k)+ &
     & dom(ib)%ae(i,j,k)*dom(ib)%u(i+1,j,k) + &
     & dom(ib)%aw(i,j,k)*dom(ib)%u(i-1,j,k)+ &
     & dom(ib)%at(i,j,k)*dom(ib)%u(i,j,k+1) + &
     & dom(ib)%ab(i,j,k)*dom(ib)%u(i,j,k-1))

          dom(ib)%ustar(i,j,k)=dom(ib)%ustar(i,j,k)+dt*alfapr*diff

        end do ; end do ; end do 

        if (L_LSM) then! .or. L_LSMbase) then
          dom(ib)%ustar(:,:,:) = dom(ib)%ustar(:,:,:) + dt*alfapr*(grx-grz*sin(atan(slope)))
        end if

        if(pressureforce) then
	   if(L_LSM) then
            do k=ks,ke ; do i=is,ie ; do j=js,je
	       if(dom(ib)%phi(i,j,k) .ge. 0.0) then
 	         dom(ib)%ustar(i,j,k) = dom(ib)%ustar(i,j,k) + dt*alfapr*forcn
	       endif
            end do ; end do ; end do 

	   else if (L_LSMbase) then
            do k=ks,ke ; do i=is,ie ; do j=js,je
		if(dom(ib)%zc(k).le.length) then
		  dom(ib)%ustar(i,j,k) = dom(ib)%ustar(i,j,k) + dt*alfapr*forcn
		endif
            end do ; end do ; end do 

	   else
!            do k=ks,ke ; do i=is,ie ; do j=js,je
 	       dom(ib)%ustar(:,:,:) = dom(ib)%ustar(:,:,:) + dt*alfapr*forcn
!            end do ; end do ; end do 
	   endif
	 endif ! pressureforce


        ! ADDING THE DIFFUSION TO VSTAR:
        is=dom(ib)%isv; ie=dom(ib)%iev
        js=dom(ib)%jsv; je=dom(ib)%jev
        ks=dom(ib)%ksv; ke=dom(ib)%kev
        do k=ks,ke ; do i=is,ie ; do j=js,je

          visc=dom(ib)%vis(i,j+1,k)

          dom(ib)%aw(i,j,k)=-visc/dxx
          dom(ib)%ae(i,j,k)=-visc/dxx
          dom(ib)%an(i,j,k)=-visc/dyy
          dom(ib)%as(i,j,k)=-visc/dyy
          dom(ib)%at(i,j,k)=-visc/dzz
          dom(ib)%ab(i,j,k)=-visc/dzz

          dom(ib)%ap(i,j,k) = -1.0*( &
     &dom(ib)%aw(i,j,k)+dom(ib)%ae(i,j,k)+ &
     &dom(ib)%as(i,j,k)+dom(ib)%an(i,j,k)+ &
     &dom(ib)%ab(i,j,k)+dom(ib)%at(i,j,k))

          diff=-(dom(ib)%ap(i,j,k)*dom(ib)%v(i,j,k)+ &
     & dom(ib)%an(i,j,k)*dom(ib)%v(i,j+1,k) + &
     & dom(ib)%as(i,j,k)*dom(ib)%v(i,j-1,k)+ &
     & dom(ib)%ae(i,j,k)*dom(ib)%v(i+1,j,k) + &
     & dom(ib)%aw(i,j,k)*dom(ib)%v(i-1,j,k)+ &
     & dom(ib)%at(i,j,k)*dom(ib)%v(i,j,k+1) + &
     & dom(ib)%ab(i,j,k)*dom(ib)%v(i,j,k-1))

           dom(ib)%vstar(i,j,k)=dom(ib)%vstar(i,j,k)+dt*alfapr*diff
        end do ; end do ; end do

        if(L_LSM) then! .or. L_LSMbase) then
          dom(ib)%vstar(:,:,:) = dom(ib)%vstar(i,j,k) + dt*alfapr*gry
        end if

        if(pressureforce_y) then				!Pablo 03/2018
	   if(L_LSM) then
            do k=ks,ke ; do i=is,ie ; do j=js,je
	       if(dom(ib)%phi(i,j,k) .ge. 0.0) then
 	         dom(ib)%vstar(i,j,k) = dom(ib)%vstar(i,j,k) + dt*alfapr*forcn_y
	       endif
            end do ; end do ; end do

	   else if(L_LSMbase) then
            do k=ks,ke ; do i=is,ie ; do j=js,je
		if (dom(ib)%zc(k).le.length) then
 	         dom(ib)%vstar(i,j,k) = dom(ib)%vstar(i,j,k) + dt*alfapr*forcn_y
		endif
            end do ; end do ; end do

	   else
!            do k=ks,ke ; do i=is,ie ; do j=js,je
 	       dom(ib)%vstar(:,:,:) = dom(ib)%vstar(:,:,:) + dt*alfapr*forcn_y
!            end do ; end do ; end do
	   endif
	 endif



        ! ADDING THE DIFFUSION TO WSTAR:
        is=dom(ib)%isw; ie=dom(ib)%iew
        js=dom(ib)%jsw; je=dom(ib)%jew
        ks=dom(ib)%ksw; ke=dom(ib)%kew
        do k=ks,ke ; do i=is,ie ; do j=js,je

          visc=dom(ib)%vis(i,j,k+1)

          dom(ib)%aw(i,j,k)=-visc/dxx
          dom(ib)%ae(i,j,k)=-visc/dxx
          dom(ib)%an(i,j,k)=-visc/dyy
          dom(ib)%as(i,j,k)=-visc/dyy
          dom(ib)%at(i,j,k)=-visc/dzz
          dom(ib)%ab(i,j,k)=-visc/dzz

          dom(ib)%ap(i,j,k) = -1.0*( &
     &dom(ib)%aw(i,j,k)+dom(ib)%ae(i,j,k)+ &
     &dom(ib)%as(i,j,k)+dom(ib)%an(i,j,k)+ &
     &dom(ib)%ab(i,j,k)+dom(ib)%at(i,j,k))

          diff=-(dom(ib)%ap(i,j,k)*dom(ib)%w(i,j,k)+ &
     & dom(ib)%an(i,j,k)*dom(ib)%w(i,j+1,k) + &
     & dom(ib)%as(i,j,k)*dom(ib)%w(i,j-1,k)+ &
     & dom(ib)%ae(i,j,k)*dom(ib)%w(i+1,j,k) + &
     & dom(ib)%aw(i,j,k)*dom(ib)%w(i-1,j,k)+ &
     & dom(ib)%at(i,j,k)*dom(ib)%w(i,j,k+1) + &
     & dom(ib)%ab(i,j,k)*dom(ib)%w(i,j,k-1))

          dom(ib)%wstar(i,j,k)=dom(ib)%wstar(i,j,k)+dt*alfapr*diff
          
!          if (L_LSM) then
          !Artificial damping zone Absorbing waves
!	   X_R_e = (dom(ib)%xc(i)-3.2d0)/0.8d0                              !sponge layer zone,artificial viscosity (2L)
	
!            if (X_R_e.ge.0.0)  then
!	     gamma = (exp(X_R_e**2.0)-1.0)/(exp(1.0)-1.0)
!	     dom(ib)%wstar(i,j,k)=dom(ib)%wstar(i,j,k)+dt*alfapr*grz*
!     & cos(atan(slope))-dt*alfapr*(100.0+0.0*abs(dom(ib)%wstar(i,j,k)))
!     & *gamma*dom(ib)%wstar(i,j,k)
!          end if

        end do ; end do ; end do 

        if(L_LSM) then
          dom(ib)%wstar(:,:,:) = dom(ib)%wstar(:,:,:) + dt*alfapr*grz*cos(atan(slope))
        end if

        if(pressureforce_z) then
	   if(L_LSM) then
            do k=ks,ke ; do i=is,ie ; do j=js,je
	       if(dom(ib)%phi(i,j,k) .ge. 0.0) then
 	         dom(ib)%wstar(i,j,k) = dom(ib)%wstar(i,j,k) + dt*alfapr*forcn_z
	       endif
            end do ; end do ; end do 
	   else if (L_LSMbase) then
            do k=ks,ke ; do i=is,ie ; do j=js,je
		if (dom(ib)%zc(k).le.length) then
 	         dom(ib)%wstar(i,j,k) = dom(ib)%wstar(i,j,k) + dt*alfapr*forcn_z
		endif
            end do ; end do ; end do 
	   else
!            do k=ks,ke ; do i=is,ie ; do j=js,je
 	       dom(ib)%wstar(:,:,:) = dom(ib)%wstar(:,:,:) + dt*alfapr*forcn_z
!            end do ; end do ; end do 
	   end if
	 end if ! pressureforce_z

      END DO ! ib

      if (LENERGY) call mom_buo
      if (LROUGH) call rough_velocity

      RETURN
      END SUBROUTINE rungek_diff2nd
!#######################################################################
      SUBROUTINE rungek_diff4th
!-----------------------------------------------------------------------
!     As part of the fractional step method we add the diffusion spatial term
!     onto the free-divergence field (ustar,vstar,wstar)
!     Note 1: LSM, Periodic, Scalar, Rough_bed compliant
!#######################################################################
      use vars
      use mpi
      use multidata
      implicit none

      INTEGER :: i,j,k,ib
      INTEGER :: is,ie,js,je,ks,ke
      DOUBLE PRECISION :: dxx,dyy,dzz,visc,diff

      fac=1.0

      do ib=1,nbp
        dxx=dom(ib)%dx*dom(ib)%dx
        dyy=dom(ib)%dy*dom(ib)%dy
        dzz=dom(ib)%dz*dom(ib)%dz


        ! ADDING THE DIFFUSION TO USTAR:
        is=dom(ib)%isu; ie=dom(ib)%ieu
        js=dom(ib)%jsu; je=dom(ib)%jeu
        ks=dom(ib)%ksu; ke=dom(ib)%keu
        do k=ks,ke ; do i=is,ie ; do j=js,je

          visc=dom(ib)%vis(i,j,k)

          diff=visc*( &
     &((-dom(ib)%u(i+2,j,k)+16.0*dom(ib)%u(i+1,j,k) &
     &  -30.0*dom(ib)%u(i,j,k)+16.0*dom(ib)%u(i-1,j,k) &
     &  -dom(ib)%u(i-2,j,k))/(12.0*dxx))+ &
     &((-dom(ib)%u(i,j+2,k)+16.0*dom(ib)%u(i,j+1,k) &
     &  -30.0*dom(ib)%u(i,j,k)+16.0*dom(ib)%u(i,j-1,k) &
     &  -dom(ib)%u(i,j-2,k))/(12.0*dyy))+ &
     &((-dom(ib)%u(i,j,k+2)+16.0*dom(ib)%u(i,j,k+1) &
     &  -30.0*dom(ib)%u(i,j,k)+16.0*dom(ib)%u(i,j,k-1) &
     &  -dom(ib)%u(i,j,k-2))/(12.0*dzz)))

          dom(ib)%ustar(i,j,k) = dom(ib)%ustar(i,j,k) + dt*alfapr*diff

        end do ; end do ; end do


        if (L_LSM) then! .or. L_LSMbase) then
          dom(ib)%ustar(:,:,:) = dom(ib)%ustar(:,:,:) + dt*alfapr*(grx-grz*sin(atan(slope)))
        end if

        if(pressureforce) then
	   if(L_LSM) then
            do k=ks,ke ; do i=is,ie ; do j=js,je
	       if(dom(ib)%phi(i,j,k) .ge. 0.0) then
 	         dom(ib)%ustar(i,j,k) = dom(ib)%ustar(i,j,k) + dt*alfapr*forcn
	       endif
            end do ; end do ; end do 

	   else if (L_LSMbase) then
            do k=ks,ke ; do i=is,ie ; do j=js,je
		if(dom(ib)%zc(k).le.length) then
		  dom(ib)%ustar(i,j,k) = dom(ib)%ustar(i,j,k) + dt*alfapr*forcn
		endif
            end do ; end do ; end do 

	   else
!            do k=ks,ke ; do i=is,ie ; do j=js,je
 	       dom(ib)%ustar(:,:,:) = dom(ib)%ustar(:,:,:) + dt*alfapr*forcn
!            end do ; end do ; end do 
	   endif
	 endif ! pressureforce



        ! ADDING THE DIFFUSION TO VSTAR:
        is=dom(ib)%isv; ie=dom(ib)%iev
        js=dom(ib)%jsv; je=dom(ib)%jev
        ks=dom(ib)%ksv; ke=dom(ib)%kev
        do k=ks,ke ; do i=is,ie ; do j=js,je

          visc=dom(ib)%vis(i,j+1,k)

          diff=visc*( &
     &((-dom(ib)%v(i+2,j,k)+16.0*dom(ib)%v(i+1,j,k) &
     &  -30.0*dom(ib)%v(i,j,k)+16.0*dom(ib)%v(i-1,j,k) &
     &  -dom(ib)%v(i-2,j,k))/(12.0*dxx))+ &
     &((-dom(ib)%v(i,j+2,k)+16.0*dom(ib)%v(i,j+1,k) &
     &  -30.0*dom(ib)%v(i,j,k)+16.0*dom(ib)%v(i,j-1,k) &
     &  -dom(ib)%v(i,j-2,k))/(12.0*dyy))+ &
     &((-dom(ib)%v(i,j,k+2)+16.0*dom(ib)%v(i,j,k+1) &
     &  -30.0*dom(ib)%v(i,j,k)+16.0*dom(ib)%v(i,j,k-1) &
     &  -dom(ib)%v(i,j,k-2))/(12.0*dzz)))

          dom(ib)%vstar(i,j,k) = dom(ib)%vstar(i,j,k) + dt*alfapr*diff
        end do ; end do ; end do 


        if(L_LSM) then! .or. L_LSMbase) then
          dom(ib)%vstar(:,:,:) = dom(ib)%vstar(:,:,:) + dt*alfapr*gry
        end if

        if(pressureforce_y) then				!Pablo 03/2018
	   if(L_LSM) then
            do k=ks,ke ; do i=is,ie ; do j=js,je
	       if(dom(ib)%phi(i,j,k) .ge. 0.0) then
 	         dom(ib)%vstar(i,j,k) = dom(ib)%vstar(i,j,k) + dt*alfapr*forcn_y
	       endif
            end do ; end do ; end do

	   else if(L_LSMbase) then
            do k=ks,ke ; do i=is,ie ; do j=js,je
		if (dom(ib)%zc(k).le.length) then
 	         dom(ib)%vstar(i,j,k) = dom(ib)%vstar(i,j,k) + dt*alfapr*forcn_y
		endif
            end do ; end do ; end do

	   else
!            do k=ks,ke ; do i=is,ie ; do j=js,je
 	       dom(ib)%vstar(:,:,:)=dom(ib)%vstar(:,:,:) + dt*alfapr*forcn_y
!            end do ; end do ; end do
	   endif
	 endif



        ! ADDING THE DIFFUSION TO VSTAR:
        is=dom(ib)%isw; ie=dom(ib)%iew
        js=dom(ib)%jsw; je=dom(ib)%jew
        ks=dom(ib)%ksw; ke=dom(ib)%kew
        do k=ks,ke ; do i=is,ie ; do j=js,je
          visc=dom(ib)%vis(i,j,k+1)

          diff=visc*( &
     &((-dom(ib)%w(i+2,j,k)+16.0*dom(ib)%w(i+1,j,k) &
     &  -30.0*dom(ib)%w(i,j,k)+16.0*dom(ib)%w(i-1,j,k) &
     &  -dom(ib)%w(i-2,j,k))/(12.0*dxx))+ &
     &((-dom(ib)%w(i,j+2,k)+16.0*dom(ib)%w(i,j+1,k) &
     &  -30.0*dom(ib)%w(i,j,k)+16.0*dom(ib)%w(i,j-1,k) &
     &  -dom(ib)%w(i,j-2,k))/(12.0*dyy))+ &
     &((-dom(ib)%w(i,j,k+2)+16.0*dom(ib)%w(i,j,k+1) &
     &  -30.0*dom(ib)%w(i,j,k)+16.0*dom(ib)%w(i,j,k-1) &
     &  -dom(ib)%w(i,j,k-2))/(12.0*dzz)))

          dom(ib)%wstar(i,j,k) = dom(ib)%wstar(i,j,k) + dt*alfapr*diff
        end do ; end do ; end do 

        if(L_LSM) then
          dom(ib)%wstar(:,:,:) = dom(ib)%wstar(:,:,:) + dt*alfapr*grz*cos(atan(slope))
        end if

        if(pressureforce_z) then
	   if(L_LSM) then
            do k=ks,ke ; do i=is,ie ; do j=js,je
	       if(dom(ib)%phi(i,j,k) .ge. 0.0) then
 	         dom(ib)%wstar(i,j,k) = dom(ib)%wstar(i,j,k) + dt*alfapr*forcn_z
	       endif
            end do ; end do ; end do 
	   else if (L_LSMbase) then
            do k=ks,ke ; do i=is,ie ; do j=js,je
		if (dom(ib)%zc(k).le.length) then
 	         dom(ib)%wstar(i,j,k) = dom(ib)%wstar(i,j,k) + dt*alfapr*forcn_z
		endif
            end do ; end do ; end do 
	   else
!            do k=ks,ke ; do i=is,ie ; do j=js,je
 	       dom(ib)%wstar(:,:,:) = dom(ib)%wstar(:,:,:) + dt*alfapr*forcn_z
!            end do ; end do ; end do 
	   end if
	 end if ! pressureforce_z

      END DO ! ib

      if (LENERGY) call mom_buo
      if (LROUGH) call rough_velocity

      RETURN
      END SUBROUTINE rungek_diff4th
!#######################################################################
