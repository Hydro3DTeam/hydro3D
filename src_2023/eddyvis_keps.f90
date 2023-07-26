!#######################################################################
! SUBROUTINES: - eddyv_keps
!              - eddyv_k
!              - eddyv_eps
!              - boundksgs
!              - boundeps
!#######################################################################
      SUBROUTINE eddyv_keps
!-----------------------------------------------------------------------
!     Calculates the subgrid-scale (SGS) viscosity for the unfiltered
!     scale. This model is known as the One equation k-eps model.
!     It includes low-Re near wall damping (Lam and Bremhorst, 1981).
!     Bruño Fraga Bugallo
!     Cardiff 2016
!#######################################################################
      use vars
      use multidata
      implicit none

      INTEGER :: i,j,k,rk
      INTEGER :: ib,is,ie,js,je,ks,ke
      DOUBLE PRECISION :: cmu
      DOUBLE PRECISION :: rrey,dx,dy,dz
      DOUBLE PRECISION :: alfark(3),strain

        alfark(1)=1./3.
        !alfark(1)=1./4.
        alfark(2)=0.5
       ! alfark(2)=0.4
        alfark(3)=1.0
       ! alfark(3)=0.8
	  !alfark(4)=1.0

	  !call exchange(11)
	  !call exchange(22)
	  !call exchange(33)

      do rk=1,3
	 call eddyv_k(alfark(rk))
	 call eddyv_eps(alfark(rk))
      enddo

       rrey=1.0/Re
       cmu=0.09

       do ib=1,nbp

       is=dom(ib)%isp; ie=dom(ib)%iep
       js=dom(ib)%jsp; je=dom(ib)%jep
       ks=dom(ib)%ksp; ke=dom(ib)%kep

	do k=ks,ke ; do j=js,je ; do i=is,ie

!	if (abs(dom(ib)%eps(i,j,k)).lt.1.0d-07) then
!		dom(ib)%eps(i,j,k)=1.0d-07
!	endif

!        if (L_LSM) rrey=dom(ib)%mu(i,j,k)/dom(ib)%dens(i,j,k)

        	dom(ib)%vis(i,j,k) = rrey + &
     &	cmu*dom(ib)%ksgs(i,j,k)**2.0/dom(ib)%eps(i,j,k)				!laminar + eddy viscosity

        	dom(ib)%vis(i,j,k) = min(dom(ib)%vis(i,j,k), &
     &	0.33*dom(ib)%ksgs(i,j,k)*sqrt(1.5/strain(i,j,k)))
	enddo ; enddo ; enddo

	     call exchange(7)

           dx=dom(ib)%dx
           dy=dom(ib)%dy
           dz=dom(ib)%dz

!..........................................................................
!=== West ===> ..  4=wall  ..   1=Inlet
!..........................................................................
        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%bc_west.eq.4) then
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    dom(ib)%vis(is,j,k)= rrey
                    dom(ib)%vis(is-1,j,k)= dom(ib)%vis(is,j,k)
                 end do
              end do
           else if (dom(ib)%bc_west.eq.1) then
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    dom(ib)%vis(is-1,j,k)= dom(ib)%vis(is,j,k)
                 end do
              end do
           else if (dom(ib)%bc_west.ge.63) then
	     call wall_function(1,dom(ib)%bc_west)
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    !dom(ib)%vis(is,j,k)= dom(ib)%tauww2(j,k)*dx		!m2/s
                    dom(ib)%vis(is-1,j,k)= dom(ib)%vis(is,j,k)
                 end do
              end do
           else if (dom(ib)%bc_west.eq.61.or.dom(ib)%bc_west.eq.62) then
	     call log_law(1,dom(ib)%bc_west)
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    !dom(ib)%vis(is,j,k)= dom(ib)%tauww2(j,k)*dx		!m2/s
                    dom(ib)%vis(is-1,j,k)= dom(ib)%vis(is,j,k)
                 end do
              end do
           else
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    dom(ib)%vis(is-1,j,k)= dom(ib)%vis(is,j,k)
                 end do
              end do
           end if
        else
           do k=ks-1,ke+1
              do j=js-1,je+1
                 dom(ib)%vis(is-1,j,k)= dom(ib)%vis(is,j,k)
              end do
           end do
        end if
!..........................................................................
!=== East ===> ..  4=wall  ..   2=Outflow
!..........................................................................
        if (dom(ib)%inext.lt.0) then
           if (dom(ib)%bc_east.eq.4) then
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    dom(ib)%vis(ie,j,k)= rrey
                    dom(ib)%vis(ie+1,j,k)= dom(ib)%vis(ie,j,k)
                 end do
              end do
           else if (dom(ib)%bc_east.eq.2) then
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    dom(ib)%vis(ie+1,j,k)= dom(ib)%vis(ie,j,k)
                 end do
              end do
           else if (dom(ib)%bc_east.ge.63) then
	     call wall_function(2,dom(ib)%bc_east)
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    !dom(ib)%vis(ie,j,k)= dom(ib)%tauwe2(j,k)*dx		!m2/s
                    dom(ib)%vis(ie+1,j,k)= dom(ib)%vis(ie,j,k)
                 end do
              end do
           else if (dom(ib)%bc_east.eq.61.or.dom(ib)%bc_east.eq.62) then
	     call log_law(2,dom(ib)%bc_east)
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    !dom(ib)%vis(ie,j,k)= dom(ib)%tauwe2(j,k)*dx		!m2/s
                    dom(ib)%vis(ie+1,j,k)= dom(ib)%vis(ie,j,k)
                 end do
              end do
           else
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    dom(ib)%vis(ie+1,j,k)= dom(ib)%vis(ie,j,k)
                 end do
              end do
           end if
        else
           do k=ks-1,ke+1
              do j=js-1,je+1
                 dom(ib)%vis(ie+1,j,k)= dom(ib)%vis(ie,j,k)
              end do
           end do
        end if
!..........................................................................
!=== South ===> ..   4=wall ..   3=Symmetry
!..........................................................................
        if (dom(ib)%jprev.lt.0) then
           if (dom(ib)%bc_south.eq.4) then
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,js,k)= rrey
                    dom(ib)%vis(i,js-1,k)= dom(ib)%vis(i,js,k)
                 end do
              end do
           else if (dom(ib)%bc_south.ge.63) then
	     call wall_function(3,dom(ib)%bc_south)
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                   ! dom(ib)%vis(i,js,k)= dom(ib)%tauws2(i,k)*dy		!m2/s
                    dom(ib)%vis(i,js-1,k)= dom(ib)%vis(i,js,k)
                 end do
              end do
           else if(dom(ib)%bc_south.eq.61.or.dom(ib)%bc_south.eq.62)then
	     call log_law(3,dom(ib)%bc_south)
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                    !dom(ib)%vis(i,js,k)= dom(ib)%tauws2(i,k)*dy		!m2/s
                    dom(ib)%vis(i,js-1,k)= dom(ib)%vis(i,js,k)
                 end do
              end do
           else
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,js-1,k)= dom(ib)%vis(i,js,k)
                 end do
              end do
           end if
        else
           do k=ks-1,ke+1
              do i=is-1,ie+1
                 dom(ib)%vis(i,js-1,k)= dom(ib)%vis(i,js,k)
              end do
           end do
        end if
!..........................................................................
!=== North ===>  ..   4=wall ..   44=moving wall ..  3=Symmetry
!..........................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%bc_north.eq.4) then
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,je,k) = rrey
                    dom(ib)%vis(i,je+1,k) =  dom(ib)%vis(i,je,k)
                 end do
              end do
           else if (dom(ib)%bc_north.eq.3) then
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,je+1,k) =  dom(ib)%vis(i,je,k)
                 end do
              end do
           else if (dom(ib)%bc_north.ge.63) then
	     call wall_function(4,dom(ib)%bc_north)
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                   ! dom(ib)%vis(i,je,k)= dom(ib)%tauwn2(i,k)*dy		!m2/s
                    dom(ib)%vis(i,je+1,k) =  dom(ib)%vis(i,je,k)
                 end do
              end do
           else if(dom(ib)%bc_north.eq.61.or.dom(ib)%bc_north.eq.62)then
	     call log_law(4,dom(ib)%bc_north)
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                  !  dom(ib)%vis(i,je,k)= dom(ib)%tauwn2(i,k)*dy		!m2/s
                    dom(ib)%vis(i,je+1,k) =  dom(ib)%vis(i,je,k)
                 end do
              end do
           else
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,je+1,k) =  dom(ib)%vis(i,je,k)
                 end do
              end do
           end if
        else
           do k=ks-1,ke+1
              do i=is-1,ie+1
                 dom(ib)%vis(i,je+1,k) =  dom(ib)%vis(i,je,k)
              end do
           end do
        end if
!..........................................................................
!=== Bottom ===> ..   4=wall ..   3=Symmetry
!..........................................................................
        if (dom(ib)%kprev.lt.0) then
           if (dom(ib)%bc_bottom.eq.4) then
              do j=js-1,je+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,j,ks)= rrey
                    dom(ib)%vis(i,j,ks-1)= dom(ib)%vis(i,j,ks)
                 end do
              end do
           else if (dom(ib)%bc_bottom.eq.3) then
              do j=js-1,je+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,j,ks-1)= dom(ib)%vis(i,j,ks)
                 end do
              end do
           else if (dom(ib)%bc_bottom.ge.63) then
	     call wall_function(5,dom(ib)%bc_bottom)
              do j=js-1,je+1
                 do i=is-1,ie+1
                   ! dom(ib)%vis(i,j,ks)= dom(ib)%tauwb2(i,j)*dz
                    dom(ib)%vis(i,j,ks-1)= dom(ib)%vis(i,j,ks)
                 end do
              end do
           else if(dom(ib)%bc_bottom.eq.61.or. &
     &					dom(ib)%bc_bottom.eq.62)then
	     call log_law(5,dom(ib)%bc_bottom)
              do j=js-1,je+1
                 do i=is-1,ie+1
                  !  dom(ib)%vis(i,j,ks)= dom(ib)%tauwb2(i,j)*dz
                    dom(ib)%vis(i,j,ks-1)= dom(ib)%vis(i,j,ks)
                 end do
              end do
           else
              do j=js-1,je+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,j,ks-1)= dom(ib)%vis(i,j,ks)
                 end do
              end do
           end if
        else
           do j=js-1,je+1
              do i=is-1,ie+1
                 dom(ib)%vis(i,j,ks-1)= dom(ib)%vis(i,j,ks)
              end do
           end do
        end if
!..........................................................................
!=== Top ===>  ..   4=wall ..     3=Symmetry
!..........................................................................
        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%bc_top.eq.4) then
              do j=js-1,je+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,j,ke)   = rrey
                    dom(ib)%vis(i,j,ke+1)   = dom(ib)%vis(i,j,ke)
                 end do
              end do
           else if (dom(ib)%bc_top.eq.3) then
              do j=js-1,je+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,j,ke+1)   = dom(ib)%vis(i,j,ke)
                 end do
              end do
           else if (dom(ib)%bc_top.ge.63) then
	     call wall_function(6,dom(ib)%bc_top)
              do j=js-1,je+1
                 do i=is-1,ie+1
                   ! dom(ib)%vis(i,j,ke)= dom(ib)%tauwt2(i,j)*dz
                    dom(ib)%vis(i,j,ke+1)   = dom(ib)%vis(i,j,ke)
                 end do
              end do
           else if(dom(ib)%bc_top.eq.61.or. &
     &					dom(ib)%bc_top.eq.62)then
	     call log_law(6,dom(ib)%bc_top)
              do j=js-1,je+1
                 do i=is-1,ie+1
                   ! dom(ib)%vis(i,j,ke)   = dom(ib)%tauwt2(i,j)*dz
                    dom(ib)%vis(i,j,ke+1)   = dom(ib)%vis(i,j,ke)
                 end do
              end do
           else
              do j=js-1,je+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,j,ke+1)   = dom(ib)%vis(i,j,ke)
                 end do
              end do
           end if
        else
           do j=js-1,je+1
              do i=is-1,ie+1
                 dom(ib)%vis(i,j,ke+1)   = dom(ib)%vis(i,j,ke)
              end do
           end do
        end if

        end do


      RETURN
      END SUBROUTINE eddyv_keps
!#######################################################################
      SUBROUTINE eddyv_k(alfark)
!-----------------------------------------------------------------------
!     Calculates the turbulent kinetic transport equation. It uses an
!     upwind scheme.
!#######################################################################
      use vars
      use multidata
      implicit none

      INTEGER :: i,j,k
      INTEGER :: ib,is,ie,js,je,ks,ke
      DOUBLE PRECISION :: dukdx,dvkdy,dwkdz,kp,km,ku,kc,kd,b_r
      DOUBLE PRECISION :: awT,aeT,asT,anT,abT,atT,apT
      DOUBLE PRECISION :: dxx,dyy,dzz,vsgs
      DOUBLE PRECISION :: visc_w,visc_e,visc_s,visc_n,visc_b,visc_t
      DOUBLE PRECISION :: sigmak,cmu
      DOUBLE PRECISION :: rrey,conv,diff,prod,other
      DOUBLE PRECISION :: alfark
      DOUBLE PRECISION :: strain

        rrey=1.0/Re
        sigmak=1.00
	 cmu = 0.09

        do ib=1,nbp

              do k=1,dom(ib)%ttc_k
              do i=1,dom(ib)%ttc_i
              do j=1,dom(ib)%ttc_j
                 dom(ib)%ksgso(i,j,k)=dom(ib)%ksgs(i,j,k)
                 dom(ib)%epso(i,j,k)=dom(ib)%eps(i,j,k)
			if (dom(ib)%epso(i,j,k).lt.1.0d-07) then
				dom(ib)%epso(i,j,k)=1.0d-07
			endif
			if (dom(ib)%ksgso(i,j,k).lt.1.0d-07) then
				dom(ib)%ksgso(i,j,k)=1.0d-07
			endif
              end do
              end do
              end do

	!if (dom_id(ib).eq.8) then
	!write(6,*)'k',dom(ib)%ksgso(20,20,20)
	!write(6,*)'eps',dom(ib)%epso(20,20,20)
	!write(6,*)'vis',dom(ib)%vis(20,20,20)
	!endif

        dxx=dom(ib)%dx*dom(ib)%dx
        dyy=dom(ib)%dy*dom(ib)%dy
        dzz=dom(ib)%dz*dom(ib)%dz

        is=dom(ib)%isp; ie=dom(ib)%iep
        js=dom(ib)%jsp; je=dom(ib)%jep
        ks=dom(ib)%ksp; ke=dom(ib)%kep

        do k=ks,ke
           do j=js,je
              do i=is,ie

        vsgs=cmu*dom(ib)%ksgso(i,j,k)**2.0/dom(ib)%epso(i,j,k)		!eddy viscosity based on t-1

        prod=0.8*2.0*vsgs*strain(i,j,k)

	 prod=min(prod,20.*dom(ib)%epso(i,j,k))

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if(dom(ib)%uoo(i-1,j,k).gt.0.0) then
           ku=dom(ib)%ksgso(i-2,j,k)
           kc=dom(ib)%ksgso(i-1,j,k)
           kd=dom(ib)%ksgso(i,j,k)
           b_r=max(0.0, &
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=kc+0.5*b_r*(kc-ku)
        else if(dom(ib)%uoo(i-1,j,k).lt.0.0) then
           ku=dom(ib)%ksgso(i+1,j,k)
           kc=dom(ib)%ksgso(i,j,k)
           kd=dom(ib)%ksgso(i-1,j,k)
           b_r=max(0.0, &
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=kc+0.5*b_r*(kc-ku)
        else
           km=0.5*(dom(ib)%ksgso(i,j,k)+dom(ib)%ksgso(i-1,j,k))
        end if
        if(dom(ib)%uoo(i,j,k).gt.0.0) then
           ku=dom(ib)%ksgso(i-1,j,k)
           kc=dom(ib)%ksgso(i,j,k)
           kd=dom(ib)%ksgso(i+1,j,k)
           b_r=max(0.0, &
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=kc+0.5*b_r*(kc-ku)
        else if(dom(ib)%uoo(i,j,k).lt.0.0) then
           ku=dom(ib)%ksgso(i+2,j,k)
           kc=dom(ib)%ksgso(i+1,j,k)
           kd=dom(ib)%ksgso(i,j,k)
           b_r=max(0.0, &
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=kc+0.5*b_r*(kc-ku)
        else
           kp=0.5*(dom(ib)%ksgso(i,j,k)+dom(ib)%ksgso(i+1,j,k))
        end if
        dukdx=(dom(ib)%uoo(i,j,k)*kp-dom(ib)%uoo(i-1,j,k)*km) &
     &	/dom(ib)%dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(dom(ib)%voo(i,j-1,k).gt.0.0) then
           ku=dom(ib)%ksgso(i,j-2,k)
           kc=dom(ib)%ksgso(i,j-1,k)
           kd=dom(ib)%ksgso(i,j,k)
           b_r=max(0.0, &
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=kc+0.5*b_r*(kc-ku)
        else if(dom(ib)%voo(i,j-1,k).lt.0.0) then
           ku=dom(ib)%ksgso(i,j+1,k)
           kc=dom(ib)%ksgso(i,j,k)
           kd=dom(ib)%ksgso(i,j-1,k)
           b_r=max(0.0, &
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=kc+0.5*b_r*(kc-ku)
        else
           km=0.5*(dom(ib)%ksgso(i,j,k)+dom(ib)%ksgso(i,j-1,k))
        end if
        if(dom(ib)%voo(i,j,k).gt.0.0) then
           ku=dom(ib)%ksgso(i,j-1,k)
           kc=dom(ib)%ksgso(i,j,k)
           kd=dom(ib)%ksgso(i,j+1,k)
           b_r=max(0.0, &
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=kc+0.5*b_r*(kc-ku)
        else if(dom(ib)%voo(i,j,k).lt.0.0) then
           ku=dom(ib)%ksgso(i,j+2,k)
           kc=dom(ib)%ksgso(i,j+1,k)
           kd=dom(ib)%ksgso(i,j,k)
           b_r=max(0.0, &
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=kc+0.5*b_r*(kc-ku)
        else
           kp=0.5*(dom(ib)%ksgso(i,j,k)+dom(ib)%ksgso(i,j+1,k))
        end if
        dvkdy=(dom(ib)%voo(i,j,k)*kp-dom(ib)%voo(i,j-1,k)*km) &
     &	/dom(ib)%dy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(dom(ib)%woo(i,j,k-1).gt.0.0) then
           ku=dom(ib)%ksgso(i,j,k-2)
           kc=dom(ib)%ksgso(i,j,k-1)
           kd=dom(ib)%ksgso(i,j,k)
           b_r=max(0.0, &
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=kc+0.5*b_r*(kc-ku)
        else if(dom(ib)%woo(i,j,k-1).lt.0.0) then
           ku=dom(ib)%ksgso(i,j,k+1)
           kc=dom(ib)%ksgso(i,j,k)
           kd=dom(ib)%ksgso(i,j,k-1)
           b_r=max(0.0, &
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=kc+0.5*b_r*(kc-ku)
        else
           km=0.5*(dom(ib)%ksgso(i,j,k)+dom(ib)%ksgso(i,j,k-1))
        end if
        if(dom(ib)%woo(i,j,k).gt.0.0) then
           ku=dom(ib)%ksgso(i,j,k-1)
           kc=dom(ib)%ksgso(i,j,k)
           kd=dom(ib)%ksgso(i,j,k+1)
           b_r=max(0.0, &
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=kc+0.5*b_r*(kc-ku)
        else if(dom(ib)%woo(i,j,k).lt.0.0) then
           ku=dom(ib)%ksgso(i,j,k+2)
           kc=dom(ib)%ksgso(i,j,k+1)
           kd=dom(ib)%ksgso(i,j,k)
           b_r=max(0.0, &
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=kc+0.5*b_r*(kc-ku)
        else
           kp=0.5*(dom(ib)%ksgso(i,j,k)+dom(ib)%ksgso(i,j,k+1))
        end if
        dwkdz=(dom(ib)%woo(i,j,k)*kp-dom(ib)%woo(i,j,k-1)*km) &
     &	/dom(ib)%dz
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        conv=(dukdx+dvkdy+dwkdz)


        visc_w=0.5*cmu*									 & !does rrey have to be here? According to Lars, yes
     & (dom(ib)%ksgso(i,j,k)**2.0/dom(ib)%epso(i,j,k) &
     & +dom(ib)%ksgso(i-1,j,k)**2.0/dom(ib)%epso(i-1,j,k))
        visc_e=0.5*cmu* &
     & (dom(ib)%ksgso(i,j,k)**2.0/dom(ib)%epso(i,j,k) &
     & +dom(ib)%ksgso(i+1,j,k)**2.0/dom(ib)%epso(i+1,j,k))
        visc_s=0.5*cmu* &
     & (dom(ib)%ksgso(i,j,k)**2.0/dom(ib)%epso(i,j,k) &
     & +dom(ib)%ksgso(i,j-1,k)**2.0/dom(ib)%epso(i,j-1,k))
        visc_n=0.5*cmu* &
     & (dom(ib)%ksgso(i,j,k)**2.0/dom(ib)%epso(i,j,k) &
     & +dom(ib)%ksgso(i,j+1,k)**2.0/dom(ib)%epso(i,j+1,k))
        visc_b=0.5*cmu* &
     & (dom(ib)%ksgso(i,j,k)**2.0/dom(ib)%epso(i,j,k) &
     & +dom(ib)%ksgso(i,j,k-1)**2.0/dom(ib)%epso(i,j,k-1))
        visc_t=0.5*cmu* &
     & (dom(ib)%ksgso(i,j,k)**2.0/dom(ib)%epso(i,j,k) &
     & +dom(ib)%ksgso(i,j,k+1)**2.0/dom(ib)%epso(i,j,k+1))

        awT=visc_w/(dxx*sigmak); aeT=visc_e/(dxx*sigmak)
        anT=visc_n/(dyy*sigmak); asT=visc_s/(dyy*sigmak)
        atT=visc_t/(dzz*sigmak); abT=visc_b/(dzz*sigmak)
        apT = -1.0*(awT+aeT+asT+anT+abT+atT)
        diff=(apT*dom(ib)%ksgso(i,j,k)+ &
     & anT*dom(ib)%ksgso(i,j+1,k) + asT*dom(ib)%ksgso(i,j-1,k)+ &
     & aeT*dom(ib)%ksgso(i+1,j,k) + awT*dom(ib)%ksgso(i-1,j,k)+ &
     & atT*dom(ib)%ksgso(i,j,k+1) + abT*dom(ib)%ksgso(i,j,k-1))

        other=dom(ib)%epso(i,j,k)

        dom(ib)%ksgs(i,j,k)=(dom(ib)%ksgso(i,j,k)+ &
     & alfark*dt*(diff-conv+prod-other))

        !if(dom(ib)%ksgs(i,j,k).lt.0.0) then
          ! print*,'ERRORRRR in eps eqn'
           !print*,dom_id(ib),i,j,k
		!write(6,*)'k',dom(ib)%ksgs(i,j,k)
		!write(6,*)'ko',dom(ib)%ksgso(i,j,k)		
		!write(6,*)'epso',dom(ib)%epso(i,j,k)
		!write(6,*)'u',dom(ib)%uoo(i,j,k)
		!write(6,*)'v',dom(ib)%voo(i,j,k)
		!write(6,*)'w',dom(ib)%woo(i,j,k)
		!write(6,*)'diff',diff
		!write(6,*)'conv',conv
		!write(6,*)'prod',prod
		!write(6,*)'dest',other
		!write(6,*)(diff-conv+prod-other)
 		!call tecplot_p(itime)
          !stop
        !end if

              end do
           end do
        end do

        end do

        call exchange(6)
        call boundksgs(cmu)

      RETURN
      END SUBROUTINE eddyv_k
!#######################################################################
      SUBROUTINE eddyv_eps(alfark)
!-----------------------------------------------------------------------
!     Calculates the dissipation transport equation using a upwind
!     scheme.
!#######################################################################
      use vars
      use multidata
      implicit none

      INTEGER :: i,j,k
      INTEGER :: ib,is,ie,js,je,ks,ke
      DOUBLE PRECISION :: epsdudx,epsdvdy,epsdwdz
      DOUBLE PRECISION :: epsp,epsm,epsu,epsc,epsd,b_r
      DOUBLE PRECISION :: awT,aeT,asT,anT,abT,atT,apT
      DOUBLE PRECISION :: dxx,dyy,dzz,vsgs
      DOUBLE PRECISION :: visc_w,visc_e,visc_s,visc_n,visc_b,visc_t
      DOUBLE PRECISION :: rrey
      DOUBLE PRECISION :: cmu,c1eps,c2eps,sigmaeps
      DOUBLE PRECISION :: conv,diff,prod,other
      DOUBLE PRECISION :: alfark,strain

        c1eps=1.44; c2eps=1.92 ; sigmaeps=1.31 ; cmu=0.09
        rrey=1.0/Re

        do ib=1,nbp				
			
              do k=1,dom(ib)%ttc_k
              do i=1,dom(ib)%ttc_i
              do j=1,dom(ib)%ttc_j
                 dom(ib)%ksgso(i,j,k)=dom(ib)%ksgs(i,j,k)
			if (dom(ib)%ksgso(i,j,k).lt.1.0d-07) then
				dom(ib)%ksgso(i,j,k)=1.0d-07
			endif
              end do
              end do
              end do

        dxx=dom(ib)%dx*dom(ib)%dx
        dyy=dom(ib)%dy*dom(ib)%dy
        dzz=dom(ib)%dz*dom(ib)%dz


        is=dom(ib)%isp; ie=dom(ib)%iep
        js=dom(ib)%jsp; je=dom(ib)%jep
        ks=dom(ib)%ksp; ke=dom(ib)%kep

        do k=ks,ke
           do j=js,je
              do i=is,ie

        vsgs=cmu*dom(ib)%ksgso(i,j,k)**2.0/dom(ib)%epso(i,j,k)		!eddy viscosity based on t-1

        prod=c1eps*(dom(ib)%epso(i,j,k)/dom(ib)%ksgso(i,j,k))		 & !production term based on t-1
     &	 *2.0*vsgs*strain(i,j,k)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Convection: Locating upstream, central and downstream nodes and interpolating k using previous time step
        if(dom(ib)%uoo(i-1,j,k).gt.0.0) then
           epsu=dom(ib)%epso(i-2,j,k)
           epsc=dom(ib)%epso(i-1,j,k)
           epsd=dom(ib)%epso(i,j,k)
           b_r=max(0.0, &
     & 	min(2.0*((epsd-epsc)/(epsc-epsu)), &
     &	0.75*((epsd-epsc)/(epsc-epsu))+0.25,4.0))
           epsm=epsc+0.5*b_r*(epsc-epsu)						!what's b_r?
        else if(dom(ib)%uoo(i-1,j,k).lt.0.0) then
           epsu=dom(ib)%epso(i+1,j,k)
           epsc=dom(ib)%epso(i,j,k)
           epsd=dom(ib)%epso(i-1,j,k)
           b_r=max(0.0, &
     &	min(2.0*((epsd-epsc)/(epsc-epsu)), &
     &	0.75*((epsd-epsc)/(epsc-epsu))+0.25,4.0))
           epsm=epsc+0.5*b_r*(epsc-epsu)
        else
           epsm=0.5*(dom(ib)%epso(i,j,k)+dom(ib)%epso(i-1,j,k))
        end if
        if(dom(ib)%uoo(i,j,k).gt.0.0) then
           epsu=dom(ib)%epso(i-1,j,k)
           epsc=dom(ib)%epso(i,j,k)
           epsd=dom(ib)%epso(i+1,j,k)
           b_r=max(0.0, &
     & 	min(2.0*((epsd-epsc)/(epsc-epsu)), &
     &	0.75*((epsd-epsc)/(epsc-epsu))+0.25,4.0))
           epsp=epsc+0.5*b_r*(epsc-epsu)
        else if(dom(ib)%uoo(i,j,k).lt.0.0) then
           epsu=dom(ib)%epso(i+2,j,k)
           epsc=dom(ib)%epso(i+1,j,k)
           epsd=dom(ib)%epso(i,j,k)
           b_r=max(0.0, &
     & 	min(2.0*((epsd-epsc)/(epsc-epsu)), &
     &	0.75*((epsd-epsc)/(epsc-epsu))+0.25,4.0))
           epsp=epsc+0.5*b_r*(epsc-epsu)
        else
           epsp=0.5*(dom(ib)%epso(i,j,k)+dom(ib)%epso(i+1,j,k))
        end if
        epsdudx=(dom(ib)%uoo(i,j,k)*epsp-dom(ib)%uoo(i-1,j,k)*epsm)		 & !calculating convective term
     &	    /dom(ib)%dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(dom(ib)%voo(i,j-1,k).gt.0.0) then
           epsu=dom(ib)%epso(i,j-2,k)
           epsc=dom(ib)%epso(i,j-1,k)
           epsd=dom(ib)%epso(i,j,k)
           b_r=max(0.0, &
     & 	min(2.0*((epsd-epsc)/(epsc-epsu)), &
     &	0.75*((epsd-epsc)/(epsc-epsu))+0.25,4.0))
           epsm=epsc+0.5*b_r*(epsc-epsu)
        else if(dom(ib)%voo(i,j-1,k).lt.0.0) then
           epsu=dom(ib)%epso(i,j+1,k)
           epsc=dom(ib)%epso(i,j,k)
           epsd=dom(ib)%epso(i,j-1,k)
           b_r=max(0.0, &
     & 	min(2.0*((epsd-epsc)/(epsc-epsu)), &
     &	0.75*((epsd-epsc)/(epsc-epsu))+0.25,4.0))
           epsm=epsc+0.5*b_r*(epsc-epsu)
        else
           epsm=0.5*(dom(ib)%epso(i,j,k)+dom(ib)%epso(i,j-1,k))
        end if
        if(dom(ib)%voo(i,j,k).gt.0.0) then
           epsu=dom(ib)%epso(i,j-1,k)
           epsc=dom(ib)%epso(i,j,k)
           epsd=dom(ib)%epso(i,j+1,k)
           b_r=max(0.0, &
     & 	min(2.0*((epsd-epsc)/(epsc-epsu)), &
     &	0.75*((epsd-epsc)/(epsc-epsu))+0.25,4.0))
           epsp=epsc+0.5*b_r*(epsc-epsu)
        else if(dom(ib)%voo(i,j,k).lt.0.0) then
           epsu=dom(ib)%epso(i,j+2,k)
           epsc=dom(ib)%epso(i,j+1,k)
           epsd=dom(ib)%epso(i,j,k)
           b_r=max(0.0, &
     & 	min(2.0*((epsd-epsc)/(epsc-epsu)), &
     &	0.75*((epsd-epsc)/(epsc-epsu))+0.25,4.0))
           epsp=epsc+0.5*b_r*(epsc-epsu)
        else
           epsp=0.5*(dom(ib)%epso(i,j,k)+dom(ib)%epso(i,j+1,k))
        end if
        epsdvdy=(dom(ib)%voo(i,j,k)*epsp-dom(ib)%voo(i,j-1,k)*epsm) &
     &	    /dom(ib)%dy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(dom(ib)%woo(i,j,k-1).gt.0.0) then
           epsu=dom(ib)%epso(i,j,k-2)
           epsc=dom(ib)%epso(i,j,k-1)
           epsd=dom(ib)%epso(i,j,k)
           b_r=max(0.0, &
     & 	min(2.0*((epsd-epsc)/(epsc-epsu)), &
     &	0.75*((epsd-epsc)/(epsc-epsu))+0.25,4.0))
           epsm=epsc+0.5*b_r*(epsc-epsu)
        else if(dom(ib)%woo(i,j,k-1).lt.0.0) then
           epsu=dom(ib)%epso(i,j,k+1)
           epsc=dom(ib)%epso(i,j,k)
           epsd=dom(ib)%epso(i,j,k-1)
           b_r=max(0.0, &
     & 	min(2.0*((epsd-epsc)/(epsc-epsu)), &
     &	0.75*((epsd-epsc)/(epsc-epsu))+0.25,4.0))
           epsm=epsc+0.5*b_r*(epsc-epsu)
        else
           epsm=0.5*(dom(ib)%epso(i,j,k)+dom(ib)%epso(i,j,k-1))
        end if
        if(dom(ib)%woo(i,j,k).gt.0.0) then
           epsu=dom(ib)%epso(i,j,k-1)
           epsc=dom(ib)%epso(i,j,k)
           epsd=dom(ib)%epso(i,j,k+1)
           b_r=max(0.0, &
     & 	min(2.0*((epsd-epsc)/(epsc-epsu)), &
     &	0.75*((epsd-epsc)/(epsc-epsu))+0.25,4.0))
           epsp=epsc+0.5*b_r*(epsc-epsu)
        else if(dom(ib)%woo(i,j,k).lt.0.0) then
           epsu=dom(ib)%epso(i,j,k+2)
           epsc=dom(ib)%epso(i,j,k+1)
           epsd=dom(ib)%epso(i,j,k)
           b_r=max(0.0, &
     & 	min(2.0*((epsd-epsc)/(epsc-epsu)), &
     &	0.75*((epsd-epsc)/(epsc-epsu))+0.25,4.0))
           epsp=epsc+0.5*b_r*(epsc-epsu)
        else
           epsp=0.5*(dom(ib)%epso(i,j,k)+dom(ib)%epso(i,j,k+1))
        end if
        epsdwdz=(dom(ib)%woo(i,j,k)*epsp-dom(ib)%woo(i,j,k-1)*epsm) &
     &	    /dom(ib)%dz

        conv=(epsdudx+epsdvdy+epsdwdz)						!convective term

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Difusion

        visc_w=0.5*cmu*									 & !does rrey have to be here? According to Lars, yes
     & (dom(ib)%ksgso(i,j,k)**2.0/dom(ib)%epso(i,j,k) &
     & +dom(ib)%ksgso(i-1,j,k)**2.0/dom(ib)%epso(i-1,j,k))
        visc_e=0.5*cmu* &
     & (dom(ib)%ksgso(i,j,k)**2.0/dom(ib)%epso(i,j,k) &
     & +dom(ib)%ksgso(i+1,j,k)**2.0/dom(ib)%epso(i+1,j,k))
        visc_s=0.5*cmu* &
     & (dom(ib)%ksgso(i,j,k)**2.0/dom(ib)%epso(i,j,k) &
     & +dom(ib)%ksgso(i,j-1,k)**2.0/dom(ib)%epso(i,j-1,k))
        visc_n=0.5*cmu* &
     & (dom(ib)%ksgso(i,j,k)**2.0/dom(ib)%epso(i,j,k) &
     & +dom(ib)%ksgso(i,j+1,k)**2.0/dom(ib)%epso(i,j+1,k))
        visc_b=0.5*cmu* &
     & (dom(ib)%ksgso(i,j,k)**2.0/dom(ib)%epso(i,j,k) &
     & +dom(ib)%ksgso(i,j,k-1)**2.0/dom(ib)%epso(i,j,k-1))
        visc_t=0.5*cmu* &
     & (dom(ib)%ksgso(i,j,k)**2.0/dom(ib)%epso(i,j,k) &
     & +dom(ib)%ksgso(i,j,k+1)**2.0/dom(ib)%epso(i,j,k+1))

        awT=visc_w/(dxx*sigmaeps); aeT=visc_e/(dxx*sigmaeps)
        anT=visc_n/(dyy*sigmaeps); asT=visc_s/(dyy*sigmaeps)
        atT=visc_t/(dzz*sigmaeps); abT=visc_b/(dzz*sigmaeps)
        apT = -1.0*(awT+aeT+asT+anT+abT+atT)			
        diff=(apT*dom(ib)%epso(i,j,k)+						 & !diffusive term
     & anT*dom(ib)%epso(i,j+1,k) + asT*dom(ib)%epso(i,j-1,k)+ &
     & aeT*dom(ib)%epso(i+1,j,k) + awT*dom(ib)%epso(i-1,j,k)+ &
     & atT*dom(ib)%epso(i,j,k+1) + abT*dom(ib)%epso(i,j,k-1))

        other=c2eps*dom(ib)%epso(i,j,k)**2.0/dom(ib)%ksgso(i,j,k)			!destruction

        dom(ib)%eps(i,j,k)=(dom(ib)%epso(i,j,k)+				 & !time derivative
     & alfark*dt*(diff-conv+prod-other))

        !if(dom(ib)%eps(i,j,k).lt.0.0) then
           !print*,'ERRORRRR in eps eqn'
           !print*,dom_id(ib),i,j,k
		!write(6,*)'eps',dom(ib)%eps(i,j,k)
		!write(6,*)'epso',dom(ib)%epso(i,j,k)		!why this doesnt print out?
		!write(6,*)'ko',dom(ib)%ksgso(i,j,k)
		!write(6,*)'u',dom(ib)%uoo(i,j,k)
		!write(6,*)'v',dom(ib)%voo(i,j,k)
		!write(6,*)'w',dom(ib)%woo(i,j,k)
		!write(6,*)'diff',diff
		!write(6,*)'conv',conv
		!write(6,*)'prod',prod
		!write(6,*)'dest',other
		!write(6,*)(diff-conv+prod-other)
 		!call tecplot_p(itime)
           !stop
       ! end if

              end do
           end do
        end do

        end do

        call exchange(9)
        call boundeps

      RETURN
      END SUBROUTINE eddyv_eps
!#######################################################################
      SUBROUTINE boundksgs(cmu)
!-----------------------------------------------------------------------
!     Precribes the kinetic energy boundary condition for the domain
!     limits (West,East,South,North,Bottom,Top).
!			Bruño Fraga Bugallo
!			Cardiff 2016
!#######################################################################
      use vars
      use multidata
      implicit none

      INTEGER :: i,j,k,ib,ly
      INTEGER :: is,ie,js,je,ks,ke
      DOUBLE PRECISION :: cmu

        if (PERIODIC) call exchange_bc(6,pl_ex)

        do ly=0,1 ! because eddy_keps needs to know within a buffer of 2 cells at each side
        do ib=1,nbp

           is=dom(ib)%isp; ie=dom(ib)%iep
           js=dom(ib)%jsp; je=dom(ib)%jep
           ks=dom(ib)%ksp; ke=dom(ib)%kep
	
! Boundary Conditions for ksgs
!.......................................................................
!=== West ===>
!.......................................................................
        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%bc_west.eq.4) then
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(is-1-ly,j,k)= 0.d0					!k=0 at the boundary?
           end do; end do
           else if (dom(ib)%bc_west.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(1,dom(ib)%bc_west)
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(is,j,k)=dom(ib)%tauww(j,k)/sqrt(cmu)		!k=ustar**2/cmu**0.5
              dom(ib)%ksgs(is-1-ly,j,k)=0.d0			
           end do; end do
           else if (dom(ib)%bc_west.eq.61 &
     &				.or.dom(ib)%bc_west.eq.62) then
	     if (ly.eq.0) call log_law(1,dom(ib)%bc_west)
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(is,j,k)=dom(ib)%tauww(j,k)/sqrt(cmu)
              dom(ib)%ksgs(is-1-ly,j,k)=0.d0
           end do; end do
           else if (dom(ib)%bc_west.eq.1) then	
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(is-1-ly,j,k)= &
     &	  (3.d0/2.d0)*(ubulk*0.1)**2.0
           end do; end do
           else if (dom(ib)%bc_west.ne.5) then					!if 5->exchange, if 2,3 -> dk/dn=0
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(is-1-ly,j,k)= dom(ib)%ksgs(is,j,k)
           end do; end do
           end if
        end if
!.......................................................................
!=== East ===>
!.......................................................................
        if (dom(ib)%inext.lt.0) then
           if (dom(ib)%bc_east.eq.4) then
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(ie+1+ly,j,k)= 0.d0
           end do; end do
           else if (dom(ib)%bc_east.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(2,dom(ib)%bc_east)
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(ie,j,k)=dom(ib)%tauwe(j,k)/sqrt(cmu)
              dom(ib)%ksgs(ie+1+ly,j,k)=0.d0
           end do; end do
           else if (dom(ib)%bc_east.eq.61 &
     &				.or.dom(ib)%bc_east.eq.62) then
		if (ly.eq.0) call log_law(2,dom(ib)%bc_east)
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(ie,j,k)=dom(ib)%tauwe(j,k)/sqrt(cmu)
              dom(ib)%ksgs(ie+1+ly,j,k)=0.d0
           end do; end do
           else if (dom(ib)%bc_east.ne.5) then
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(ie+1+ly,j,k)= dom(ib)%ksgs(ie,j,k)
           end do; end do
           end if
        end if
!.......................................................................
!=== South ===>
!.......................................................................
        if (dom(ib)%jprev.lt.0) then
           if (dom(ib)%bc_south.eq.4) then
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,js-1-ly,k)= 0.d0
           end do; end do
           else if (dom(ib)%bc_south.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(3,dom(ib)%bc_south)
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,js,k)=dom(ib)%tauws(i,k)/sqrt(cmu)
              dom(ib)%ksgs(i,js-1-ly,k)=0.d0
           end do; end do
           else if (dom(ib)%bc_south.eq.61 &
     &				.or.dom(ib)%bc_south.eq.62) then
		if (ly.eq.0) call log_law(3,dom(ib)%bc_south)
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,js,k)=dom(ib)%tauws(i,k)/sqrt(cmu)
              dom(ib)%ksgs(i,js-1-ly,k)=0.d0
           end do; end do
           else if (dom(ib)%bc_south.ne.5) then
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,js-1-ly,k)= dom(ib)%ksgs(i,js,k)
           end do; end do
           end if
        end if
!.......................................................................
!=== North ===>
!.......................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%bc_north.eq.4) then
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,je+1+ly,k) = 0.d0
           end do; end do
           else if (dom(ib)%bc_north.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(4,dom(ib)%bc_north)
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,je,k)=dom(ib)%tauwn(i,k)/sqrt(cmu)
              dom(ib)%ksgs(i,je+1+ly,k) =0.d0
           end do; end do
           else if (dom(ib)%bc_north.eq.61 &
     &				.or.dom(ib)%bc_north.eq.62) then
		if (ly.eq.0) call log_law(4,dom(ib)%bc_north)
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,je,k)=dom(ib)%tauwn(i,k)/sqrt(cmu)
              dom(ib)%ksgs(i,je+1+ly,k) =0.d0
           end do; end do
           else if (dom(ib)%bc_north.ne.5) then
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,je+1+ly,k) = dom(ib)%ksgs(i,je,k)	
           end do; end do
           end if
        end if
!.......................................................................
!=== Bottom ===>
!.......................................................................
        if (dom(ib)%kprev.lt.0) then
           if (dom(ib)%bc_bottom.eq.4) then
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,j,ks-1-ly)= 0.d0
           end do; end do
           else if (dom(ib)%bc_bottom.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(5,dom(ib)%bc_bottom)
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,j,ks)=dom(ib)%tauwb(i,j)/sqrt(cmu)
              dom(ib)%ksgs(i,j,ks-1-ly)=0.d0
           end do; end do
           else if (dom(ib)%bc_bottom.eq.61 &
     &				.or.dom(ib)%bc_bottom.eq.62) then
		if (ly.eq.0) then
			call log_law(5,dom(ib)%bc_bottom)
		endif
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,j,ks)=dom(ib)%tauwb(i,j)/sqrt(cmu)
              dom(ib)%ksgs(i,j,ks-1-ly)= 0.d0
           end do; end do
           else if (dom(ib)%bc_bottom.ne.5) then
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,j,ks-1-ly)= dom(ib)%ksgs(i,j,ks)
           end do; end do
           end if
        end if
!.......................................................................
!=== Top ===>
!.......................................................................
        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%bc_top.eq.4) then
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,j,ke+1+ly) = 0.d0
           end do; end do
           else if (dom(ib)%bc_top.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(6,dom(ib)%bc_top)
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,j,ke)=dom(ib)%tauwt(i,j)/sqrt(cmu)
              dom(ib)%ksgs(i,j,ke+1+ly) =0.d0
           end do; end do
           else if (dom(ib)%bc_top.eq.61 &
     &				.or.dom(ib)%bc_top.eq.62) then
		if (ly.eq.0) call log_law(6,dom(ib)%bc_top)
           do j=js-1,je+1; do i=is-1,ie+1
             dom(ib)%ksgs(i,j,ke)=dom(ib)%tauwt(i,j)/sqrt(cmu)
              dom(ib)%ksgs(i,j,ke+1+ly) =0.d0
           end do; end do
           else if (dom(ib)%bc_top.ne.5) then
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,j,ke+1+ly) = dom(ib)%ksgs(i,j,ke)
           end do; end do
           end if
        end if

        end do
        end do

      RETURN
      END SUBROUTINE boundksgs
!#######################################################################
      SUBROUTINE boundeps
!     Precribes the epsilon boundary condition for the domain
!     limits (West,East,South,North,Bottom,Top).
!#######################################################################
      use vars
      use multidata
      implicit none

      INTEGER :: i,j,k,ib,ly
      INTEGER :: is,ie,js,je,ks,ke
      DOUBLE PRECISION :: rrey,rk1,drkdy,lz
      DOUBLE PRECISION :: dx,dy,dz,kappa,delta

        rrey=1.0/Re
	  kappa=0.41
	  lz=zen-zst

        if (PERIODIC) call exchange_bc(9,pl_ex)

        do ly=0,2
        do ib=1,nbp

           dx=dom(ib)%dx
           dy=dom(ib)%dy
           dz=dom(ib)%dz

           is=dom(ib)%isp; ie=dom(ib)%iep
           js=dom(ib)%jsp; je=dom(ib)%jep
           ks=dom(ib)%ksp; ke=dom(ib)%kep
	
! Boundary Conditions for epsylon
!.......................................................................
!=== West ===>
!.......................................................................
        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%bc_west.eq.4) then
           do k=ks-1,ke+1; do j=js-1,je+1
              rk1 = sqrt(dom(ib)%ksgs(is,j,k))
              drkdy = rk1 /(0.5*dom(ib)%dx)
              dom(ib)%eps(is-1-ly,j,k)= 2.d0*rrey*drkdy**2.0
           end do; end do
           else if (dom(ib)%bc_west.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(1,dom(ib)%bc_west)
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%eps(is,j,k)= &
     &		dom(ib)%tauww(j,k)**1.5/(kappa*0.5*dx)
              dom(ib)%eps(is-1-ly,j,k)= dom(ib)%eps(is,j,k)
           end do; end do
           else if (dom(ib)%bc_west.eq.61 &
     &				.or.dom(ib)%bc_west.eq.62) then
		if (ly.eq.0) call log_law(1,dom(ib)%bc_west)
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%eps(is,j,k)= &
     &		dom(ib)%tauww(j,k)**1.5/(kappa*0.5*dx)
              dom(ib)%eps(is-1-ly,j,k)= dom(ib)%eps(is,j,k)
           end do; end do
           else if (dom(ib)%bc_west.eq.1) then	
           do k=ks-1,ke+1; do j=js-1,je+1			
              dom(ib)%eps(is-1-ly,j,k) = &
     &	  0.09**0.75*dom(ib)%ksgs(is-1-ly,j,k)**1.5/(0.07*lz)
           end do; end do
           else if (dom(ib)%bc_west.ne.5) then
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%eps(is-1-ly,j,k)= dom(ib)%eps(is,j,k)
           end do; end do
           end if
        end if
!.......................................................................
!=== East ===>
!.......................................................................
        if (dom(ib)%inext.lt.0) then
           if (dom(ib)%bc_east.eq.4) then
           do k=ks-1,ke+1; do j=js-1,je+1
              rk1 = sqrt(dom(ib)%ksgs(ie,j,k))
              drkdy = rk1 /(0.5*dom(ib)%dx)
              dom(ib)%eps(ie+1+ly,j,k)= 2.d0*rrey*drkdy**2.0
           end do; end do
           else if (dom(ib)%bc_east.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(2,dom(ib)%bc_east)
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%eps(ie,j,k)= &
     &		dom(ib)%tauwe(j,k)**1.5/(kappa*0.5*dx)
              dom(ib)%eps(ie+1+ly,j,k)= dom(ib)%eps(ie,j,k)
           end do; end do
           else if (dom(ib)%bc_east.eq.61 &
     &				.or.dom(ib)%bc_east.eq.62) then
		if (ly.eq.0) call log_law(2,dom(ib)%bc_east)
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%eps(ie,j,k)= &
     &		dom(ib)%tauwe(j,k)**1.5/(kappa*0.5*dx)
              dom(ib)%eps(ie+1+ly,j,k)= dom(ib)%eps(ie,j,k)
           end do; end do
           else if (dom(ib)%bc_east.ne.5) then
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%eps(ie+1+ly,j,k)= dom(ib)%eps(ie,j,k)
           end do; end do
           end if
        end if
!.......................................................................
!=== South ===>
!.......................................................................
        if (dom(ib)%jprev.lt.0) then
           if (dom(ib)%bc_south.eq.4) then
           do k=ks-1,ke+1; do i=is-1,ie+1
              rk1 = sqrt(dom(ib)%ksgs(i,js,k))
              drkdy = rk1 /(0.5*dom(ib)%dy)
              dom(ib)%eps(i,js-1-ly,k)= 2.d0*rrey*drkdy**2.0	
           end do; end do
           else if (dom(ib)%bc_south.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(3,dom(ib)%bc_south)
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%eps(i,js,k)= &
     &		dom(ib)%tauws(i,k)**1.5/(kappa*0.5*dy)
              dom(ib)%eps(i,js-1-ly,k)= dom(ib)%eps(i,js,k)
           end do; end do
           else if (dom(ib)%bc_south.eq.61 &
     &				.or.dom(ib)%bc_south.eq.62) then
		if (ly.eq.0) call log_law(3,dom(ib)%bc_south)
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%eps(i,js,k)= &
     &		dom(ib)%tauws(i,k)**1.5/(kappa*0.5*dy)
              dom(ib)%eps(i,js-1-ly,k)= dom(ib)%eps(i,js,k)
           end do; end do
           else if (dom(ib)%bc_south.ne.5) then
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%eps(i,js-1-ly,k)= dom(ib)%eps(i,js,k)
           end do; end do
           end if
        end if
!.......................................................................
!=== North ===>
!.......................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%bc_north.eq.4) then
           do k=ks-1,ke+1; do i=is-1,ie+1
              rk1 = sqrt(dom(ib)%ksgs(i,je,k))
              drkdy = rk1 /(0.5*dom(ib)%dy)
              dom(ib)%eps(i,je+1+ly,k)= 2.d0*rrey*drkdy**2.0
           end do; end do
           else if (dom(ib)%bc_north.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(4,dom(ib)%bc_north)
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%eps(i,je,k)= &
     &		dom(ib)%tauwn(i,k)**1.5/(kappa*0.5*dy)
              dom(ib)%eps(i,je+1+ly,k) = dom(ib)%eps(i,je,k)		
           end do; end do
           else if (dom(ib)%bc_north.eq.61 &
     &				.or.dom(ib)%bc_north.eq.62) then
		if (ly.eq.0) call log_law(4,dom(ib)%bc_north)
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%eps(i,je,k)= &
     &		dom(ib)%tauwn(i,k)**1.5/(kappa*0.5*dy)
              dom(ib)%eps(i,je+1+ly,k) = dom(ib)%eps(i,je,k)		
           end do; end do
           else if (dom(ib)%bc_north.ne.5) then
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%eps(i,je+1+ly,k) = dom(ib)%eps(i,je,k)	
           end do; end do
           end if
        end if
!.......................................................................
!=== Bottom ===>
!.......................................................................
        if (dom(ib)%kprev.lt.0) then
           if (dom(ib)%bc_bottom.eq.4) then
           do j=js-1,je+1; do i=is-1,ie+1
              rk1 = sqrt(dom(ib)%ksgs(i,j,ks))
              drkdy = rk1 /(0.5*dom(ib)%dz)
              dom(ib)%eps(i,j,ks-1-ly)= 2.d0*rrey*drkdy**2.0	
           end do; end do
           else if (dom(ib)%bc_bottom.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(5,dom(ib)%bc_bottom)
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%eps(i,j,ks)= &
     &		dom(ib)%tauwb(i,j)**1.5/(kappa*0.5*dz)
              dom(ib)%eps(i,j,ks-1-ly)= dom(ib)%eps(i,j,ks)
           end do; end do
           else if (dom(ib)%bc_bottom.eq.61 &
     &				.or.dom(ib)%bc_bottom.eq.62) then
		if (ly.eq.0) then
			call log_law(5,dom(ib)%bc_bottom)
		endif
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%eps(i,j,ks)= &
     &		dom(ib)%tauwb(i,j)**1.5/(kappa*0.5*dz)
              dom(ib)%eps(i,j,ks-1-ly)= dom(ib)%eps(i,j,ks)
           end do; end do
           else if (dom(ib)%bc_bottom.ne.5) then
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%eps(i,j,ks-1-ly)= dom(ib)%eps(i,j,ks)
           end do; end do
           end if
        end if
!.......................................................................
!=== Top ===>
!.......................................................................
        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%bc_top.eq.4) then
           do j=js-1,je+1; do i=is-1,ie+1
              rk1 = sqrt(dom(ib)%ksgs(i,j,ke))
              drkdy = rk1 /(0.5*dom(ib)%dz)
              dom(ib)%eps(i,j,ke+1+ly)= 2.d0*rrey*drkdy**2.0				
           end do; end do
           else if (dom(ib)%bc_top.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(6,dom(ib)%bc_top)
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%eps(i,j,ke)= &
     &		dom(ib)%tauwt(i,j)**1.5/(kappa*0.5*dz)
              dom(ib)%eps(i,j,ke+1+ly) = dom(ib)%eps(i,j,ke)
           end do; end do
           else if (dom(ib)%bc_top.eq.61 &
     &				.or.dom(ib)%bc_top.eq.62) then
		if (ly.eq.0) call log_law(6,dom(ib)%bc_top)
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%eps(i,j,ke)= &
     &		dom(ib)%tauwt(i,j)**1.5/(kappa*0.5*dz)
              dom(ib)%eps(i,j,ke+1+ly) = dom(ib)%eps(i,j,ke)
           end do; end do
           else if (dom(ib)%bc_top.ne.5) then
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%eps(i,j,ke+1+ly) = dom(ib)%eps(i,j,ke)
           end do; end do
           end if
        end if

        end do
        end do

      RETURN
      END SUBROUTINE boundeps
!#############################################################################
      DOUBLE PRECISION FUNCTION strain(i,j,k)
!-----------------------------------------------------------------------
!     Calculates the Strain tensor based on the velocity gradients.
!     This function is used in eddys_k and eddys_eps.
!#######################################################################
      use vars
      use multidata
      implicit none

      INTEGER :: i,j,k
      INTEGER :: ib,is,ie,js,je,ks,ke
      DOUBLE PRECISION :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
      DOUBLE PRECISION :: vr_a,vr_b,vr_c,vr_d
      DOUBLE PRECISION :: s12,s13,s23

        do ib=1,nbp				


!====================================================
        vr_a = 0.25*( dom(ib)%uoo(i,j,k)   + dom(ib)%uoo(i-1,j,k) + &
     &	 	  dom(ib)%uoo(i,j+1,k) + dom(ib)%uoo(i-1,j+1,k) )
        vr_b = 0.25*( dom(ib)%uoo(i,j,k)   + dom(ib)%uoo(i-1,j,k) + &
     &	 	   dom(ib)%uoo(i,j-1,k) + dom(ib)%uoo(i-1,j-1,k) )

        vr_c = 0.25*( dom(ib)%uoo(i,j,k)   + dom(ib)%uoo(i-1,j,k) + &
     &	 	  dom(ib)%uoo(i,j,k+1) + dom(ib)%uoo(i-1,j,k+1) )
        vr_d = 0.25*( dom(ib)%uoo(i,j,k)   + dom(ib)%uoo(i-1,j,k) + &
     &	 	  dom(ib)%uoo(i,j,k-1) + dom(ib)%uoo(i-1,j,k-1) )

        dudx = (dom(ib)%uoo(i,j,k)-dom(ib)%uoo(i-1,j,k) )/dom(ib)%dx
        dudy = ( vr_a - vr_b )/dom(ib)%dy
        dudz = ( vr_c - vr_d )/dom(ib)%dz
!====================================================
        vr_a = 0.25*( dom(ib)%voo(i,j,k)   + dom(ib)%voo(i,j-1,k) + &
     &	 	   dom(ib)%voo(i+1,j,k) + dom(ib)%voo(i+1,j-1,k) )
        vr_b = 0.25*( dom(ib)%voo(i,j,k)   + dom(ib)%voo(i,j-1,k) + &
     &	 	   dom(ib)%voo(i-1,j,k) + dom(ib)%voo(i-1,j-1,k) )

        vr_c = 0.25*( dom(ib)%voo(i,j,k)   + dom(ib)%voo(i,j-1,k) + &
     &	 	   dom(ib)%voo(i,j,k+1) + dom(ib)%voo(i,j-1,k+1) )
        vr_d = 0.25*( dom(ib)%voo(i,j,k)   + dom(ib)%voo(i,j-1,k) + &
     &	 	  dom(ib)%voo(i,j,k-1) + dom(ib)%voo(i,j-1,k-1) )

        dvdy = (dom(ib)%voo(i,j,k)-dom(ib)%voo(i,j-1,k) )/dom(ib)%dy
        dvdx = ( vr_a - vr_b )/dom(ib)%dx
        dvdz = ( vr_c - vr_d )/dom(ib)%dz
!====================================================
        vr_a = 0.25*( dom(ib)%woo(i,j,k)   + dom(ib)%woo(i,j,k-1) + &
     &	 	  dom(ib)%woo(i+1,j,k) + dom(ib)%woo(i+1,j,k-1) )
        vr_b = 0.25*( dom(ib)%woo(i,j,k)   + dom(ib)%woo(i,j,k-1) + &
     &	 	   dom(ib)%woo(i-1,j,k) + dom(ib)%woo(i-1,j,k-1) )

        vr_c = 0.25*( dom(ib)%woo(i,j,k)   + dom(ib)%woo(i,j,k-1) + &
     &	 	   dom(ib)%woo(i,j+1,k) + dom(ib)%woo(i,j+1,k-1) )
        vr_d = 0.25*( dom(ib)%woo(i,j,k)   + dom(ib)%woo(i,j,k-1) + &
     &	 	   dom(ib)%woo(i,j-1,k) + dom(ib)%woo(i,j-1,k-1) )

        dwdz = (dom(ib)%woo(i,j,k)-dom(ib)%woo(i,j,k-1) )/dom(ib)%dz
        dwdx = ( vr_a - vr_b )/dom(ib)%dx
        dwdy = ( vr_c - vr_d )/dom(ib)%dy


        s12 = 0.5 * (dudy + dvdx)
        s13 = 0.5 * (dudz + dwdx)
        s23 = 0.5 * (dvdz + dwdy)
        strain  = ( dudx*dudx + dvdy*dvdy   + dwdz*dwdz  + &
     &      2.0*s12*s12   + 2.0*s13*s13 + 2.0*s23*s23 )

	enddo

      RETURN
      END FUNCTION strain
