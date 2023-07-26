!#######################################################################
! SUBROUTINES: - convection  
!              - rungek_conv2nd
!              - rungek_conv4th
!              - rungek_convWENO
!#######################################################################
      SUBROUTINE convection
!-----------------------------------------------------------------------
!     Numerical method that calculates the convective spatial term of
!     the Navier-Stokes equations. Following the fractional step method
!     it adds it to the free-divergence field (ustar,vstar,wstar)
!     Note 1: convection scheme for 1.Explicit Euler, 2.Adams Bashfort
!     Note 2: each convection scheme can use the 3 differentiating scheme
!#######################################################################
      use vars
      use multidata
      implicit none

      INTEGER :: i,j,k,ib,ii,jj,kk,is,ie,js,je,ks,ke,tti,ttj,ttk
      DOUBLE PRECISION :: vel
      DOUBLE PRECISION :: up12,um12,vp12,vm12,wp12,wm12
      DOUBLE PRECISION :: uijk,vijk,wijk
      DOUBLE PRECISION :: dxx,dyy,dzz,dudx,dudy,dudz
      DOUBLE PRECISION :: dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
      DOUBLE PRECISION :: sp2,sp1,sn,sm1,sm2
      DOUBLE PRECISION,allocatable,dimension(:,:,:) :: du2dx,duvdy,duwdz
      DOUBLE PRECISION,allocatable,dimension(:,:,:) :: duvdx,dv2dy,dvwdz
      DOUBLE PRECISION,allocatable,dimension(:,:,:) :: duwdx,dvwdy,dw2dz

      do ib=1,nbp
        dom(ib)%ustar=dom(ib)%u
        dom(ib)%vstar=dom(ib)%v
        dom(ib)%wstar=dom(ib)%w
      end do

      if (differencing.eq.3) then
        call HJ_WENO_dx(1)
        call HJ_WENO_dy(1)
        call HJ_WENO_dz(1)
      end if
      
      ! U VELOCITY-STEAMWISE:
      DO ib=1,nbp
        is=dom(ib)%isu; ie=dom(ib)%ieu
        js=dom(ib)%jsu; je=dom(ib)%jeu
        ks=dom(ib)%ksu; ke=dom(ib)%keu

        tti=dom(ib)%ttc_i ; ttj=dom(ib)%ttc_j ; ttk=dom(ib)%ttc_k
        allocate(du2dx(tti,ttj,ttk),duvdy(tti,ttj,ttk),duwdz(tti,ttj,ttk))


        ! 2nd Order CDS:  
        IF(differencing.eq.1) THEN
          do k=ks,ke ; do i=is,ie ; do j=js,je
            up12=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i+1,j,k))
            um12=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
            du2dx(i,j,k)=(up12**2-um12**2)/dom(ib)%dx

            duvdy(i,j,k)=0.25*( (dom(ib)%v(i,j,k)+dom(ib)%v(i+1,j,k))* &
       &(dom(ib)%u(i,j,k)+dom(ib)%u(i,j+1,k)) - &
       &(dom(ib)%v(i,j-1,k)+dom(ib)%v(i+1,j-1,k))* &
       &(dom(ib)%u(i,j,k)+dom(ib)%u(i,j-1,k)) )/dom(ib)%dy

            duwdz(i,j,k)=0.25*( (dom(ib)%w(i,j,k)+dom(ib)%w(i+1,j,k))* &
       &(dom(ib)%u(i,j,k)+dom(ib)%u(i,j,k+1)) - &
       &(dom(ib)%w(i,j,k-1)+dom(ib)%w(i+1,j,k-1))* &
       &(dom(ib)%u(i,j,k)+dom(ib)%u(i,j,k-1)) )/dom(ib)%dz
          end do ; end do ; end do 

        ! 4th Order CDS:
        ELSE IF(differencing.eq.2) THEN 
          do k=ks,ke ; do i=is,ie ; do j=js,je
        ii=i+2; sp2=((-dom(ib)%u(ii+1,j,k)+9.0*dom(ib)%u(ii,j,k)+ &
     & 9.0*dom(ib)%u(ii-1,j,k)-dom(ib)%u(ii-2,j,k))/16.0)**2
        ii=i+1; sp1=((-dom(ib)%u(ii+1,j,k)+9.0*dom(ib)%u(ii,j,k)+ &
     & 9.0*dom(ib)%u(ii-1,j,k)-dom(ib)%u(ii-2,j,k))/16.0)**2
        ii=i;   sn=((-dom(ib)%u(ii+1,j,k)+9.0*dom(ib)%u(ii,j,k)+ &
     & 9.0*dom(ib)%u(ii-1,j,k)-dom(ib)%u(ii-2,j,k))/16.0)**2
        ii=i-1; sm1=((-dom(ib)%u(ii+1,j,k)+9.0*dom(ib)%u(ii,j,k)+ &
     & 9.0*dom(ib)%u(ii-1,j,k)-dom(ib)%u(ii-2,j,k))/16.0)**2
            du2dx(i,j,k)=(-sp2+27.0*sp1-27.0*sn+sm1)/(24.0*dom(ib)%dx)

        jj=j+1; sp1=((-dom(ib)%u(i,jj+2,k)+9.0*dom(ib)%u(i,jj+1,k)+ &
     & 9.0*dom(ib)%u(i,jj,k)-dom(ib)%u(i,jj-1,k))/16.0)* &
     &((-dom(ib)%v(i+2,jj,k)+9.0*dom(ib)%v(i+1,jj,k)+ &
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i-1,jj,k))/16.0)
        jj=j;   sn=((-dom(ib)%u(i,jj+2,k)+9.0*dom(ib)%u(i,jj+1,k)+ &
     & 9.0*dom(ib)%u(i,jj,k)-dom(ib)%u(i,jj-1,k))/16.0)* &
     &((-dom(ib)%v(i+2,jj,k)+9.0*dom(ib)%v(i+1,jj,k)+ &
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i-1,jj,k))/16.0)
        jj=j-1; sm1=((-dom(ib)%u(i,jj+2,k)+9.0*dom(ib)%u(i,jj+1,k)+ &
     & 9.0*dom(ib)%u(i,jj,k)-dom(ib)%u(i,jj-1,k))/16.0)* &
     &((-dom(ib)%v(i+2,jj,k)+9.0*dom(ib)%v(i+1,jj,k)+ &
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i-1,jj,k))/16.0)
        jj=j-2; sm2=((-dom(ib)%u(i,jj+2,k)+9.0*dom(ib)%u(i,jj+1,k)+ &
     & 9.0*dom(ib)%u(i,jj,k)-dom(ib)%u(i,jj-1,k))/16.0)* &
     &((-dom(ib)%v(i+2,jj,k)+9.0*dom(ib)%v(i+1,jj,k)+ &
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i-1,jj,k))/16.0)
            duvdy(i,j,k)=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dy)

        kk=k+1; sp1=((-dom(ib)%u(i,j,kk+2)+9.0*dom(ib)%u(i,j,kk+1)+ &
     & 9.0*dom(ib)%u(i,j,kk)-dom(ib)%u(i,j,kk-1))/16.0)* &
     &((-dom(ib)%w(i+2,j,kk)+9.0*dom(ib)%w(i+1,j,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i-1,j,kk))/16.0)
        kk=k;   sn=((-dom(ib)%u(i,j,kk+2)+9.0*dom(ib)%u(i,j,kk+1)+ &
     & 9.0*dom(ib)%u(i,j,kk)-dom(ib)%u(i,j,kk-1))/16.0)* &
     &((-dom(ib)%w(i+2,j,kk)+9.0*dom(ib)%w(i+1,j,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i-1,j,kk))/16.0)
        kk=k-1; sm1=((-dom(ib)%u(i,j,kk+2)+9.0*dom(ib)%u(i,j,kk+1)+ &
     & 9.0*dom(ib)%u(i,j,kk)-dom(ib)%u(i,j,kk-1))/16.0)* &
     &((-dom(ib)%w(i+2,j,kk)+9.0*dom(ib)%w(i+1,j,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i-1,j,kk))/16.0)
        kk=k-2; sm2=((-dom(ib)%u(i,j,kk+2)+9.0*dom(ib)%u(i,j,kk+1)+ &
     & 9.0*dom(ib)%u(i,j,kk)-dom(ib)%u(i,j,kk-1))/16.0)* &
     &((-dom(ib)%w(i+2,j,kk)+9.0*dom(ib)%w(i+1,j,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i-1,j,kk))/16.0)
              duwdz(i,j,k)=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dz)
            end do ; end do ; end do 
          
        ! WENO:
        ELSE IF(differencing.eq.3) THEN
          do k=ks,ke ; do i=is,ie ; do j=js,je
            if (dom(ib)%u(i,j,k).gt.0.0) then
               dudx = dom(ib)%dphi_dxminus(i,j,k)
            else if (dom(ib)%u(i,j,k).lt.0.0) then
               dudx = dom(ib)%dphi_dxplus(i,j,k)
            else
               dudx = 0.0
            end if

            vijk=0.25*(dom(ib)%v(i,j,k)+dom(ib)%v(i+1,j,k)+dom(ib)%v(i,j-1,k)+dom(ib)%v(i+1,j-1,k))
            if (vijk.gt.0.0) then
              dudy = dom(ib)%dphi_dyminus(i,j,k)
            else if (vijk.lt.0.0) then
              dudy = dom(ib)%dphi_dyplus(i,j,k)
            else
              dudy = 0.0
            end if

            wijk=0.25*(dom(ib)%w(i,j,k)+dom(ib)%w(i+1,j,k)+dom(ib)%w(i,j,k-1)+dom(ib)%w(i+1,j,k-1))
            if (wijk.gt.0.0) then
               dudz = dom(ib)%dphi_dzminus(i,j,k)
            else if (wijk.lt.0.0) then
               dudz = dom(ib)%dphi_dzplus(i,j,k)
            else
               dudz = 0.0
            end if

            du2dx(i,j,k)=dom(ib)%u(i,j,k)*dudx
            duvdy(i,j,k)=vijk*dudy
            duwdz(i,j,k)=wijk*dudz

          end do ; end do ; end do 
        END IF ! differecing


        ! IMPLICIT EULER:
        IF(conv_sch.eq.1) then
          do k=ks,ke ; do i=is,ie ; do j=js,je
            dom(ib)%ustar(i,j,k)=(dom(ib)%u(i,j,k) - dt*(du2dx(i,j,k)+duvdy(i,j,k)+duwdz(i,j,k)))
          end do ; end do ; end do 

        ! ADAM-BASHFORT
        ELSE IF(conv_sch.eq.2) then
          do k=ks,ke ; do i=is,ie ; do j=js,je
            vel=dom(ib)%uo(i,j,k)
            dom(ib)%uo(i,j,k)=(du2dx(i,j,k)+duvdy(i,j,k)+duwdz(i,j,k))
            dom(ib)%ustar(i,j,k)=(dom(ib)%u(i,j,k)-dt*(1.5*dom(ib)%uo(i,j,k)-0.5*vel))
          end do ; end do ; end do 
        END IF

        deallocate(du2dx,duvdy,duwdz)

      END DO ! ib


      if (differencing.eq.3) then
        call HJ_WENO_dx(2)
        call HJ_WENO_dy(2)
        call HJ_WENO_dz(2)
      end if

      ! V VELOCITY-TRANSVERSAL:
      DO ib=1,nbp
        is=dom(ib)%isv; ie=dom(ib)%iev
        js=dom(ib)%jsv; je=dom(ib)%jev
        ks=dom(ib)%ksv; ke=dom(ib)%kev

        tti=dom(ib)%ttc_i ; ttj=dom(ib)%ttc_j ; ttk=dom(ib)%ttc_k
        allocate(duvdx(tti,ttj,ttk),dv2dy(tti,ttj,ttk),dvwdz(tti,ttj,ttk))
        
        ! 2nd Order CDS:
        IF(differencing.eq.1) THEN
          do k=ks,ke ; do i=is,ie ; do j=js,je
            vp12=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j+1,k))
            vm12=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
            dv2dy(i,j,k)=(vp12**2-vm12**2)/dom(ib)%dy

            duvdx(i,j,k)=0.25*( (dom(ib)%u(i,j,k)+dom(ib)%u(i,j+1,k))* &
     &(dom(ib)%v(i,j,k)+dom(ib)%v(i+1,j,k)) - &
     &(dom(ib)%u(i-1,j,k)+dom(ib)%u(i-1,j+1,k))* &
     &(dom(ib)%v(i,j,k)+dom(ib)%v(i-1,j,k)) )/dom(ib)%dx

              dvwdz(i,j,k)=0.25*( (dom(ib)%w(i,j,k)+dom(ib)%w(i,j+1,k))* &
     &(dom(ib)%v(i,j,k)+dom(ib)%v(i,j,k+1)) - &
     &(dom(ib)%w(i,j,k-1)+dom(ib)%w(i,j+1,k-1))* &
     &(dom(ib)%v(i,j,k)+dom(ib)%v(i,j,k-1)) )/dom(ib)%dz
            end do ; end do ; end do 

        ! 4th Order CDS:
        ELSE IF(differencing.eq.2) THEN
          do k=ks,ke ; do i=is,ie ; do j=js,je
        jj=j+2; sp2=((-dom(ib)%v(i,jj+1,k)+9.0*dom(ib)%v(i,jj,k)+ &
     & 9.0*dom(ib)%v(i,jj-1,k)-dom(ib)%v(i,jj-2,k))/16.0)**2
        jj=j+1; sp1=((-dom(ib)%v(i,jj+1,k)+9.0*dom(ib)%v(i,jj,k)+ &
     & 9.0*dom(ib)%v(i,jj-1,k)-dom(ib)%v(i,jj-2,k))/16.0)**2
        jj=j  ; sn=((-dom(ib)%v(i,jj+1,k)+9.0*dom(ib)%v(i,jj,k)+ &
     & 9.0*dom(ib)%v(i,jj-1,k)-dom(ib)%v(i,jj-2,k))/16.0)**2
        jj=j-1; sm1=((-dom(ib)%v(i,jj+1,k)+9.0*dom(ib)%v(i,jj,k)+ &
     & 9.0*dom(ib)%v(i,jj-1,k)-dom(ib)%v(i,jj-2,k))/16.0)**2
            dv2dy(i,j,k)=(-sp2+27.0*sp1-27.0*sn+sm1)/(24.0*dom(ib)%dy)

        ii=i+1; sp1=((-dom(ib)%v(ii+2,j,k)+9.0*dom(ib)%v(ii+1,j,k)+ &
     & 9.0*dom(ib)%v(ii,j,k)-dom(ib)%v(ii-1,j,k))/16.0)* &
     &((-dom(ib)%u(ii,j+2,k)+9.0*dom(ib)%u(ii,j+1,k)+ &
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j-1,k))/16.0)
        ii=i;   sn=((-dom(ib)%v(ii+2,j,k)+9.0*dom(ib)%v(ii+1,j,k)+ &
     & 9.0*dom(ib)%v(ii,j,k)-dom(ib)%v(ii-1,j,k))/16.0)* &
     &((-dom(ib)%u(ii,j+2,k)+9.0*dom(ib)%u(ii,j+1,k)+ &
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j-1,k))/16.0)
        ii=i-1; sm1=((-dom(ib)%v(ii+2,j,k)+9.0*dom(ib)%v(ii+1,j,k)+ &
     & 9.0*dom(ib)%v(ii,j,k)-dom(ib)%v(ii-1,j,k))/16.0)* &
     &((-dom(ib)%u(ii,j+2,k)+9.0*dom(ib)%u(ii,j+1,k)+ &
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j-1,k))/16.0)
        ii=i-2; sm2=((-dom(ib)%v(ii+2,j,k)+9.0*dom(ib)%v(ii+1,j,k)+ &
     & 9.0*dom(ib)%v(ii,j,k)-dom(ib)%v(ii-1,j,k))/16.0)* &
     &((-dom(ib)%u(ii,j+2,k)+9.0*dom(ib)%u(ii,j+1,k)+ &
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j-1,k))/16.0)
            duvdx(i,j,k)=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dx)

        kk=k+1; sp1=((-dom(ib)%v(i,j,kk+2)+9.0*dom(ib)%v(i,j,kk+1)+ &
     & 9.0*dom(ib)%v(i,j,kk)-dom(ib)%v(i,j,kk-1))/16.0)* &
     &((-dom(ib)%w(i,j+2,kk)+9.0*dom(ib)%w(i,j+1,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i,j-1,kk))/16.0)
        kk=k;   sn=((-dom(ib)%v(i,j,kk+2)+9.0*dom(ib)%v(i,j,kk+1)+ &
     & 9.0*dom(ib)%v(i,j,kk)-dom(ib)%v(i,j,kk-1))/16.0)* &
     &((-dom(ib)%w(i,j+2,kk)+9.0*dom(ib)%w(i,j+1,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i,j-1,kk))/16.0)
        kk=k-1; sm1=((-dom(ib)%v(i,j,kk+2)+9.0*dom(ib)%v(i,j,kk+1)+ &
     & 9.0*dom(ib)%v(i,j,kk)-dom(ib)%v(i,j,kk-1))/16.0)* &
     &((-dom(ib)%w(i,j+2,kk)+9.0*dom(ib)%w(i,j+1,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i,j-1,kk))/16.0)
        kk=k-2; sm2=((-dom(ib)%v(i,j,kk+2)+9.0*dom(ib)%v(i,j,kk+1)+ &
     & 9.0*dom(ib)%v(i,j,kk)-dom(ib)%v(i,j,kk-1))/16.0)* &
     &((-dom(ib)%w(i,j+2,kk)+9.0*dom(ib)%w(i,j+1,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i,j-1,kk))/16.0)
          dvwdz(i,j,k)=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dz)
        end do ; end do ; end do 
          
        ! WENO:
        ELSE IF(differencing.eq.3) THEN
          do k=ks,ke ; do i=is,ie ; do j=js,je
            uijk=0.25*(dom(ib)%u(i,j,k)+dom(ib)%u(i,j+1,k)+dom(ib)%u(i-1,j,k)+dom(ib)%u(i-1,j+1,k))
            if (uijk.gt.0.0) then
              dvdx = dom(ib)%dphi_dxminus(i,j,k)
            else if (uijk.lt.0.0) then
              dvdx = dom(ib)%dphi_dxplus(i,j,k)
            else
              dvdx = 0.0
            end if

            if (dom(ib)%v(i,j,k).gt.0.0) then
              dvdy = dom(ib)%dphi_dyminus(i,j,k)
            else if (dom(ib)%v(i,j,k).lt.0.0) then
              dvdy = dom(ib)%dphi_dyplus(i,j,k)
            else
              dvdy = 0.0
            end if

            wijk=0.25*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j+1,k)+dom(ib)%w(i,j,k-1)+dom(ib)%w(i,j+1,k-1))
            if (wijk.gt.0.0) then
              dvdz = dom(ib)%dphi_dzminus(i,j,k)
            else if (wijk.lt.0.0) then
              dvdz = dom(ib)%dphi_dzplus(i,j,k)
            else
              dvdz = 0.0
            end if

            duvdx(i,j,k)=uijk*dvdx
            dv2dy(i,j,k)=dom(ib)%v(i,j,k)*dvdy
            dvwdz(i,j,k)=wijk*dvdz
          end do ; end do ; end do 
        END IF ! differencing 
        
        ! IMPLICIT EULER:
        IF(conv_sch.eq.1) THEN 
          do k=ks,ke ; do i=is,ie ; do j=js,je
            dom(ib)%vstar(i,j,k)=(dom(ib)%v(i,j,k) - dt*(duvdx(i,j,k)+dv2dy(i,j,k)+dvwdz(i,j,k)))
          end do ; end do ; end do 

        ! ADAMS-BASHFORD:
        ELSE IF(conv_sch.eq.2) THEN
          do k=ks,ke ; do i=is,ie ; do j=js,je
            vel=dom(ib)%vo(i,j,k)
            dom(ib)%vo(i,j,k)=(duvdx(i,j,k)+dv2dy(i,j,k)+dvwdz(i,j,k))
            dom(ib)%vstar(i,j,k)=(dom(ib)%v(i,j,k)-dt*(1.5*dom(ib)%vo(i,j,k)-0.5*vel))
          end do ; end do ; end do 
        END IF  ! conv_sch

        deallocate(duvdx,dv2dy,dvwdz)
      END DO ! ib

      if (differencing.eq.3) then
        call HJ_WENO_dx(3)
        call HJ_WENO_dy(3)
        call HJ_WENO_dz(3)
      end if

      DO ib=1,nbp
        is=dom(ib)%isw; ie=dom(ib)%iew
        js=dom(ib)%jsw; je=dom(ib)%jew
        ks=dom(ib)%ksw; ke=dom(ib)%kew

        tti=dom(ib)%ttc_i ; ttj=dom(ib)%ttc_j ; ttk=dom(ib)%ttc_k
        allocate(duwdx(tti,ttj,ttk),dvwdy(tti,ttj,ttk),dw2dz(tti,ttj,ttk))
        
        ! 2nd Order CDS:
        IF(differencing.eq.1) THEN
          do k=ks,ke ; do i=is,ie ; do j=js,je
            wp12=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k+1))
            wm12=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
            dw2dz(i,j,k)=(wp12**2-wm12**2)/dom(ib)%dz

            duwdx(i,j,k)=0.25*( (dom(ib)%u(i,j,k)+dom(ib)%u(i,j,k+1))* &
     &(dom(ib)%w(i,j,k)+dom(ib)%w(i+1,j,k)) - &
     &(dom(ib)%u(i-1,j,k)+dom(ib)%u(i-1,j,k+1))* &
     &(dom(ib)%w(i,j,k)+dom(ib)%w(i-1,j,k)) )/dom(ib)%dx

            dvwdy(i,j,k)=0.25*( (dom(ib)%v(i,j,k)+dom(ib)%v(i,j,k+1))* &
     &(dom(ib)%w(i,j,k)+dom(ib)%w(i,j+1,k)) - &
     &(dom(ib)%v(i,j-1,k)+dom(ib)%v(i,j-1,k+1))* &
     &(dom(ib)%w(i,j,k)+dom(ib)%w(i,j-1,k)) )/dom(ib)%dy
          end do ; end do ; end do
        
        ! 4th Order CDS:
        ELSE IF(differencing.eq.2) THEN
          do k=ks,ke ; do i=is,ie ; do j=js,je
        kk=k+2; sp2=((-dom(ib)%w(i,j,kk+1)+9.0*dom(ib)%w(i,j,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk-1)-dom(ib)%w(i,j,kk-2))/16.0)**2
        kk=k+1; sp1=((-dom(ib)%w(i,j,kk+1)+9.0*dom(ib)%w(i,j,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk-1)-dom(ib)%w(i,j,kk-2))/16.0)**2
        kk=k  ; sn=((-dom(ib)%w(i,j,kk+1)+9.0*dom(ib)%w(i,j,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk-1)-dom(ib)%w(i,j,kk-2))/16.0)**2
        kk=k-1; sm1=((-dom(ib)%w(i,j,kk+1)+9.0*dom(ib)%w(i,j,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk-1)-dom(ib)%w(i,j,kk-2))/16.0)**2
            dw2dz(i,j,k)=(-sp2+27.0*sp1-27.0*sn+sm1)/(24.0*dom(ib)%dz)

        ii=i+1; sp1=((-dom(ib)%w(ii+2,j,k)+9.0*dom(ib)%w(ii+1,j,k)+ &
     & 9.0*dom(ib)%w(ii,j,k)-dom(ib)%w(ii-1,j,k))/16.0)* &
     &((-dom(ib)%u(ii,j,k+2)+9.0*dom(ib)%u(ii,j,k+1)+ &
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j,k-1))/16.0)
        ii=i;   sn=((-dom(ib)%w(ii+2,j,k)+9.0*dom(ib)%w(ii+1,j,k)+ &
     & 9.0*dom(ib)%w(ii,j,k)-dom(ib)%w(ii-1,j,k))/16.0)* &
     &((-dom(ib)%u(ii,j,k+2)+9.0*dom(ib)%u(ii,j,k+1)+ &
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j,k-1))/16.0)
        ii=i-1; sm1=((-dom(ib)%w(ii+2,j,k)+9.0*dom(ib)%w(ii+1,j,k)+ &
     & 9.0*dom(ib)%w(ii,j,k)-dom(ib)%w(ii-1,j,k))/16.0)* &
     &((-dom(ib)%u(ii,j,k+2)+9.0*dom(ib)%u(ii,j,k+1)+ &
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j,k-1))/16.0)
        ii=i-2; sm2=((-dom(ib)%w(ii+2,j,k)+9.0*dom(ib)%w(ii+1,j,k)+ &
     & 9.0*dom(ib)%w(ii,j,k)-dom(ib)%w(ii-1,j,k))/16.0)* &
     &((-dom(ib)%u(ii,j,k+2)+9.0*dom(ib)%u(ii,j,k+1)+ &
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j,k-1))/16.0)
            duwdx(i,j,k)=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dx)

        jj=j+1; sp1=((-dom(ib)%w(i,jj+2,k)+9.0*dom(ib)%w(i,jj+1,k)+ &
     & 9.0*dom(ib)%w(i,jj,k)-dom(ib)%w(i,jj-1,k))/16.0)* &
     &((-dom(ib)%v(i,jj,k+2)+9.0*dom(ib)%v(i,jj,k+1)+ &
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i,jj,k-1))/16.0)
        jj=j;   sn=((-dom(ib)%w(i,jj+2,k)+9.0*dom(ib)%w(i,jj+1,k)+ &
     & 9.0*dom(ib)%w(i,jj,k)-dom(ib)%w(i,jj-1,k))/16.0)* &
     &((-dom(ib)%v(i,jj,k+2)+9.0*dom(ib)%v(i,jj,k+1)+ &
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i,jj,k-1))/16.0)
        jj=j-1; sm1=((-dom(ib)%w(i,jj+2,k)+9.0*dom(ib)%w(i,jj+1,k)+ &
     & 9.0*dom(ib)%w(i,jj,k)-dom(ib)%w(i,jj-1,k))/16.0)* &
     &((-dom(ib)%v(i,jj,k+2)+9.0*dom(ib)%v(i,jj,k+1)+ &
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i,jj,k-1))/16.0)
        jj=j-2; sm2=((-dom(ib)%w(i,jj+2,k)+9.0*dom(ib)%w(i,jj+1,k)+ &
     & 9.0*dom(ib)%w(i,jj,k)-dom(ib)%w(i,jj-1,k))/16.0)* &
     &((-dom(ib)%v(i,jj,k+2)+9.0*dom(ib)%v(i,jj,k+1)+ &
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i,jj,k-1))/16.0)
            dvwdy(i,j,k)=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dy)
          end do ; end do ; end do
        
        ! WENO:
        ELSE IF(differencing.eq.3) THEN
          do k=ks,ke ; do i=is,ie ; do j=js,je
            uijk=0.25*(dom(ib)%u(i,j,k)+dom(ib)%u(i,j,k+1)+dom(ib)%u(i-1,j,k)+dom(ib)%u(i-1,j,k+1))
            if (uijk.gt.0.0) then
              dwdx = dom(ib)%dphi_dxminus(i,j,k)
            else if (uijk.lt.0.0) then
              dwdx = dom(ib)%dphi_dxplus(i,j,k)
            else
              dwdx = 0.0
            end if

            vijk=0.25*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j,k+1)+dom(ib)%v(i,j-1,k)+dom(ib)%v(i,j-1,k+1))
            if(vijk.gt.0.0) then
              dwdy = dom(ib)%dphi_dyminus(i,j,k)
            else if (vijk.lt.0.0) then
              dwdy = dom(ib)%dphi_dyplus(i,j,k)
            else
              dwdy = 0.0
            end if

            if(dom(ib)%w(i,j,k).gt.0.0) then
              dwdz = dom(ib)%dphi_dzminus(i,j,k)
            else if (dom(ib)%w(i,j,k).lt.0.0) then
              dwdz = dom(ib)%dphi_dzplus(i,j,k)
            else
              dwdz = 0.0
            end if

            duwdx(i,j,k)=uijk*dwdx
            dvwdy(i,j,k)=vijk*dwdy
            dw2dz(i,j,k)=dom(ib)%w(i,j,k)*dwdz
          end do ; end do ; end do
        END IF ! differencing

        ! IMPLICIT EULER:
        IF(conv_sch.eq.1) THEN
          do k=ks,ke ; do i=is,ie ; do j=js,je
            dom(ib)%wstar(i,j,k)=(dom(ib)%w(i,j,k)-dt*(duwdx(i,j,k)+dvwdy(i,j,k)+dw2dz(i,j,k)))
          end do ; end do ; end do

        ! ADAM-BASHFORT: 
        ELSE IF(conv_sch.eq.2) THEN
          do k=ks,ke ; do i=is,ie ; do j=js,je
            vel=dom(ib)%wo(i,j,k)
            dom(ib)%wo(i,j,k)=(duwdx(i,j,k)+dvwdy(i,j,k)+dw2dz(i,j,k))
            dom(ib)%wstar(i,j,k)=(dom(ib)%w(i,j,k)-dt*(1.5*dom(ib)%wo(i,j,k)-0.5*vel))
          end do ; end do ; end do
        END IF

        deallocate(duwdx,dvwdy,dw2dz)

      END DO !ib

      RETURN
      END SUBROUTINE convection
!#######################################################################
      SUBROUTINE rungek_conv2nd
!#######################################################################
      use vars
      use mpi
      use multidata
      implicit none

      INTEGER :: i,j,k,ib
      DOUBLE PRECISION :: du2dx,dv2dy,dw2dz
      DOUBLE PRECISION :: up12,um12,vp12,vm12,wp12,wm12
      DOUBLE PRECISION :: duvdx,duvdy,duwdx,duwdz,dvwdy,dvwdz

      call boundu

      do ib=1,nbp

           do k=dom(ib)%ksu,dom(ib)%keu
              do i=dom(ib)%isu,dom(ib)%ieu
                 do j=dom(ib)%jsu,dom(ib)%jeu

                 up12=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i+1,j,k))
                 um12=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                 du2dx=(up12**2-um12**2)/dom(ib)%dx
                 duvdy=0.25*( (dom(ib)%v(i,j,k)+dom(ib)%v(i+1,j,k))* &
     &(dom(ib)%u(i,j,k)+dom(ib)%u(i,j+1,k)) - &
     &(dom(ib)%v(i,j-1,k)+dom(ib)%v(i+1,j-1,k))* &
     &(dom(ib)%u(i,j,k)+dom(ib)%u(i,j-1,k)) )/dom(ib)%dy
                 duwdz=0.25*( (dom(ib)%w(i,j,k)+dom(ib)%w(i+1,j,k))* &
     &(dom(ib)%u(i,j,k)+dom(ib)%u(i,j,k+1)) - &
     &(dom(ib)%w(i,j,k-1)+dom(ib)%w(i+1,j,k-1))* &
     &(dom(ib)%u(i,j,k)+dom(ib)%u(i,j,k-1)) )/dom(ib)%dz

                 dom(ib)%ustar(i,j,k)=(dom(ib)%uoo(i,j,k)- &
     & dt*alfapr*(du2dx+duvdy+duwdz))

                 end do
              end do
           end do
        end do

        call boundv

        do ib=1,nbp

           do k=dom(ib)%ksv,dom(ib)%kev
              do i=dom(ib)%isv,dom(ib)%iev
                 do j=dom(ib)%jsv,dom(ib)%jev

                 vp12=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j+1,k))
                 vm12=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                 dv2dy=(vp12**2-vm12**2)/dom(ib)%dy
                 duvdx=0.25*( (dom(ib)%u(i,j,k)+dom(ib)%u(i,j+1,k))* &
     &(dom(ib)%v(i,j,k)+dom(ib)%v(i+1,j,k)) - &
     &(dom(ib)%u(i-1,j,k)+dom(ib)%u(i-1,j+1,k))* &
     &(dom(ib)%v(i,j,k)+dom(ib)%v(i-1,j,k)) )/dom(ib)%dx
                 dvwdz=0.25*( (dom(ib)%w(i,j,k)+dom(ib)%w(i,j+1,k))* &
     &(dom(ib)%v(i,j,k)+dom(ib)%v(i,j,k+1)) - &
     &(dom(ib)%w(i,j,k-1)+dom(ib)%w(i,j+1,k-1))* &
     &(dom(ib)%v(i,j,k)+dom(ib)%v(i,j,k-1)) )/dom(ib)%dz

                 dom(ib)%vstar(i,j,k)=(dom(ib)%voo(i,j,k)- &
     & dt*alfapr*(duvdx+dv2dy+dvwdz))

                 end do
              end do
           end do
        end do


        call boundw

        do ib=1,nbp

           do k=dom(ib)%ksw,dom(ib)%kew
              do i=dom(ib)%isw,dom(ib)%iew
                 do j=dom(ib)%jsw,dom(ib)%jew

                 wp12=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k+1))
                 wm12=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                 dw2dz=(wp12**2-wm12**2)/dom(ib)%dz
                 duwdx=0.25*( (dom(ib)%u(i,j,k)+dom(ib)%u(i,j,k+1))* &
     &(dom(ib)%w(i,j,k)+dom(ib)%w(i+1,j,k)) - &
     &(dom(ib)%u(i-1,j,k)+dom(ib)%u(i-1,j,k+1))* &
     &(dom(ib)%w(i,j,k)+dom(ib)%w(i-1,j,k)) )/dom(ib)%dx
                 dvwdy=0.25*( (dom(ib)%v(i,j,k)+dom(ib)%v(i,j,k+1))* &
     &(dom(ib)%w(i,j,k)+dom(ib)%w(i,j+1,k)) - &
     &(dom(ib)%v(i,j-1,k)+dom(ib)%v(i,j-1,k+1))* &
     &(dom(ib)%w(i,j,k)+dom(ib)%w(i,j-1,k)) )/dom(ib)%dy

                 dom(ib)%wstar(i,j,k)=(dom(ib)%woo(i,j,k)- &
     & dt*alfapr*(duwdx+dvwdy+dw2dz))

                 end do
              end do
           end do

      end do

      RETURN
      END SUBROUTINE rungek_conv2nd
!#######################################################################
      SUBROUTINE rungek_conv4th
!#######################################################################
      use vars
      use mpi
      use multidata
      implicit none

      INTEGER :: i,j,k,ib,ii,jj,kk
      DOUBLE PRECISION :: du2dx,dv2dy,dw2dz
      DOUBLE PRECISION :: duvdx,duvdy,duwdx,duwdz,dvwdy,dvwdz
      DOUBLE PRECISION :: sp2,sp1,sn,sm1,sm2

        call boundu

        do ib=1,nbp

           do k=dom(ib)%ksu,dom(ib)%keu
              do i=dom(ib)%isu,dom(ib)%ieu
                 do j=dom(ib)%jsu,dom(ib)%jeu

        ii=i+2; sp2=((-dom(ib)%u(ii+1,j,k)+9.0*dom(ib)%u(ii,j,k)+ &
     & 9.0*dom(ib)%u(ii-1,j,k)-dom(ib)%u(ii-2,j,k))/16.0)**2
        ii=i+1; sp1=((-dom(ib)%u(ii+1,j,k)+9.0*dom(ib)%u(ii,j,k)+ &
     & 9.0*dom(ib)%u(ii-1,j,k)-dom(ib)%u(ii-2,j,k))/16.0)**2
        ii=i;   sn=((-dom(ib)%u(ii+1,j,k)+9.0*dom(ib)%u(ii,j,k)+ &
     & 9.0*dom(ib)%u(ii-1,j,k)-dom(ib)%u(ii-2,j,k))/16.0)**2
        ii=i-1; sm1=((-dom(ib)%u(ii+1,j,k)+9.0*dom(ib)%u(ii,j,k)+ &
     & 9.0*dom(ib)%u(ii-1,j,k)-dom(ib)%u(ii-2,j,k))/16.0)**2
                 du2dx=(-sp2+27.0*sp1-27.0*sn+sm1)/(24.0*dom(ib)%dx)

        jj=j+1; sp1=((-dom(ib)%u(i,jj+2,k)+9.0*dom(ib)%u(i,jj+1,k)+ &
     & 9.0*dom(ib)%u(i,jj,k)-dom(ib)%u(i,jj-1,k))/16.0)* &
     &((-dom(ib)%v(i+2,jj,k)+9.0*dom(ib)%v(i+1,jj,k)+ &
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i-1,jj,k))/16.0)
        jj=j;   sn=((-dom(ib)%u(i,jj+2,k)+9.0*dom(ib)%u(i,jj+1,k)+ &
     & 9.0*dom(ib)%u(i,jj,k)-dom(ib)%u(i,jj-1,k))/16.0)* &
     &((-dom(ib)%v(i+2,jj,k)+9.0*dom(ib)%v(i+1,jj,k)+ &
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i-1,jj,k))/16.0)
        jj=j-1; sm1=((-dom(ib)%u(i,jj+2,k)+9.0*dom(ib)%u(i,jj+1,k)+ &
     & 9.0*dom(ib)%u(i,jj,k)-dom(ib)%u(i,jj-1,k))/16.0)* &
     &((-dom(ib)%v(i+2,jj,k)+9.0*dom(ib)%v(i+1,jj,k)+ &
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i-1,jj,k))/16.0)
        jj=j-2; sm2=((-dom(ib)%u(i,jj+2,k)+9.0*dom(ib)%u(i,jj+1,k)+ &
     & 9.0*dom(ib)%u(i,jj,k)-dom(ib)%u(i,jj-1,k))/16.0)* &
     &((-dom(ib)%v(i+2,jj,k)+9.0*dom(ib)%v(i+1,jj,k)+ &
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i-1,jj,k))/16.0)
                 duvdy=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dy)

        kk=k+1; sp1=((-dom(ib)%u(i,j,kk+2)+9.0*dom(ib)%u(i,j,kk+1)+ &
     & 9.0*dom(ib)%u(i,j,kk)-dom(ib)%u(i,j,kk-1))/16.0)* &
     &((-dom(ib)%w(i+2,j,kk)+9.0*dom(ib)%w(i+1,j,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i-1,j,kk))/16.0)
        kk=k;   sn=((-dom(ib)%u(i,j,kk+2)+9.0*dom(ib)%u(i,j,kk+1)+ &
     & 9.0*dom(ib)%u(i,j,kk)-dom(ib)%u(i,j,kk-1))/16.0)* &
     &((-dom(ib)%w(i+2,j,kk)+9.0*dom(ib)%w(i+1,j,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i-1,j,kk))/16.0)
        kk=k-1; sm1=((-dom(ib)%u(i,j,kk+2)+9.0*dom(ib)%u(i,j,kk+1)+ &
     & 9.0*dom(ib)%u(i,j,kk)-dom(ib)%u(i,j,kk-1))/16.0)* &
     &((-dom(ib)%w(i+2,j,kk)+9.0*dom(ib)%w(i+1,j,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i-1,j,kk))/16.0)
        kk=k-2; sm2=((-dom(ib)%u(i,j,kk+2)+9.0*dom(ib)%u(i,j,kk+1)+ &
     & 9.0*dom(ib)%u(i,j,kk)-dom(ib)%u(i,j,kk-1))/16.0)* &
     &((-dom(ib)%w(i+2,j,kk)+9.0*dom(ib)%w(i+1,j,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i-1,j,kk))/16.0)
                 duwdz=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dz)

                 dom(ib)%ustar(i,j,k)=(dom(ib)%uoo(i,j,k)- &
     & dt*alfapr*(du2dx+duvdy+duwdz))

                 end do
              end do
           end do
        end do

        call boundv

        do ib=1,nbp

           do k=dom(ib)%ksv,dom(ib)%kev
              do i=dom(ib)%isv,dom(ib)%iev
                 do j=dom(ib)%jsv,dom(ib)%jev

        jj=j+2; sp2=((-dom(ib)%v(i,jj+1,k)+9.0*dom(ib)%v(i,jj,k)+ &
     & 9.0*dom(ib)%v(i,jj-1,k)-dom(ib)%v(i,jj-2,k))/16.0)**2
        jj=j+1; sp1=((-dom(ib)%v(i,jj+1,k)+9.0*dom(ib)%v(i,jj,k)+ &
     & 9.0*dom(ib)%v(i,jj-1,k)-dom(ib)%v(i,jj-2,k))/16.0)**2
        jj=j  ; sn=((-dom(ib)%v(i,jj+1,k)+9.0*dom(ib)%v(i,jj,k)+ &
     & 9.0*dom(ib)%v(i,jj-1,k)-dom(ib)%v(i,jj-2,k))/16.0)**2
        jj=j-1; sm1=((-dom(ib)%v(i,jj+1,k)+9.0*dom(ib)%v(i,jj,k)+ &
     & 9.0*dom(ib)%v(i,jj-1,k)-dom(ib)%v(i,jj-2,k))/16.0)**2
                 dv2dy=(-sp2+27.0*sp1-27.0*sn+sm1)/(24.0*dom(ib)%dy)

        ii=i+1; sp1=((-dom(ib)%v(ii+2,j,k)+9.0*dom(ib)%v(ii+1,j,k)+ &
     & 9.0*dom(ib)%v(ii,j,k)-dom(ib)%v(ii-1,j,k))/16.0)* &
     &((-dom(ib)%u(ii,j+2,k)+9.0*dom(ib)%u(ii,j+1,k)+ &
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j-1,k))/16.0)
        ii=i;   sn=((-dom(ib)%v(ii+2,j,k)+9.0*dom(ib)%v(ii+1,j,k)+ &
     & 9.0*dom(ib)%v(ii,j,k)-dom(ib)%v(ii-1,j,k))/16.0)* &
     &((-dom(ib)%u(ii,j+2,k)+9.0*dom(ib)%u(ii,j+1,k)+ &
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j-1,k))/16.0)
        ii=i-1; sm1=((-dom(ib)%v(ii+2,j,k)+9.0*dom(ib)%v(ii+1,j,k)+ &
     & 9.0*dom(ib)%v(ii,j,k)-dom(ib)%v(ii-1,j,k))/16.0)* &
     &((-dom(ib)%u(ii,j+2,k)+9.0*dom(ib)%u(ii,j+1,k)+ &
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j-1,k))/16.0)
        ii=i-2; sm2=((-dom(ib)%v(ii+2,j,k)+9.0*dom(ib)%v(ii+1,j,k)+ &
     & 9.0*dom(ib)%v(ii,j,k)-dom(ib)%v(ii-1,j,k))/16.0)* &
     &((-dom(ib)%u(ii,j+2,k)+9.0*dom(ib)%u(ii,j+1,k)+ &
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j-1,k))/16.0)
                 duvdx=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dx)

        kk=k+1; sp1=((-dom(ib)%v(i,j,kk+2)+9.0*dom(ib)%v(i,j,kk+1)+ &
     & 9.0*dom(ib)%v(i,j,kk)-dom(ib)%v(i,j,kk-1))/16.0)* &
     &((-dom(ib)%w(i,j+2,kk)+9.0*dom(ib)%w(i,j+1,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i,j-1,kk))/16.0)
        kk=k;   sn=((-dom(ib)%v(i,j,kk+2)+9.0*dom(ib)%v(i,j,kk+1)+ &
     & 9.0*dom(ib)%v(i,j,kk)-dom(ib)%v(i,j,kk-1))/16.0)* &
     &((-dom(ib)%w(i,j+2,kk)+9.0*dom(ib)%w(i,j+1,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i,j-1,kk))/16.0)
        kk=k-1; sm1=((-dom(ib)%v(i,j,kk+2)+9.0*dom(ib)%v(i,j,kk+1)+ &
     & 9.0*dom(ib)%v(i,j,kk)-dom(ib)%v(i,j,kk-1))/16.0)* &
     &((-dom(ib)%w(i,j+2,kk)+9.0*dom(ib)%w(i,j+1,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i,j-1,kk))/16.0)
        kk=k-2; sm2=((-dom(ib)%v(i,j,kk+2)+9.0*dom(ib)%v(i,j,kk+1)+ &
     & 9.0*dom(ib)%v(i,j,kk)-dom(ib)%v(i,j,kk-1))/16.0)* &
     &((-dom(ib)%w(i,j+2,kk)+9.0*dom(ib)%w(i,j+1,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i,j-1,kk))/16.0)
                 dvwdz=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dz)

                 dom(ib)%vstar(i,j,k)=(dom(ib)%voo(i,j,k)- &
     & dt*alfapr*(duvdx+dv2dy+dvwdz))

                 end do
              end do
           end do
        end do


        call boundw

        do ib=1,nbp

           do k=dom(ib)%ksw,dom(ib)%kew
              do i=dom(ib)%isw,dom(ib)%iew
                 do j=dom(ib)%jsw,dom(ib)%jew

        kk=k+2; sp2=((-dom(ib)%w(i,j,kk+1)+9.0*dom(ib)%w(i,j,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk-1)-dom(ib)%w(i,j,kk-2))/16.0)**2
        kk=k+1; sp1=((-dom(ib)%w(i,j,kk+1)+9.0*dom(ib)%w(i,j,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk-1)-dom(ib)%w(i,j,kk-2))/16.0)**2
        kk=k  ; sn=((-dom(ib)%w(i,j,kk+1)+9.0*dom(ib)%w(i,j,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk-1)-dom(ib)%w(i,j,kk-2))/16.0)**2
        kk=k-1; sm1=((-dom(ib)%w(i,j,kk+1)+9.0*dom(ib)%w(i,j,kk)+ &
     & 9.0*dom(ib)%w(i,j,kk-1)-dom(ib)%w(i,j,kk-2))/16.0)**2
                 dw2dz=(-sp2+27.0*sp1-27.0*sn+sm1)/(24.0*dom(ib)%dz)

        ii=i+1; sp1=((-dom(ib)%w(ii+2,j,k)+9.0*dom(ib)%w(ii+1,j,k)+ &
     & 9.0*dom(ib)%w(ii,j,k)-dom(ib)%w(ii-1,j,k))/16.0)* &
     &((-dom(ib)%u(ii,j,k+2)+9.0*dom(ib)%u(ii,j,k+1)+ &
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j,k-1))/16.0)
        ii=i;   sn=((-dom(ib)%w(ii+2,j,k)+9.0*dom(ib)%w(ii+1,j,k)+ &
     & 9.0*dom(ib)%w(ii,j,k)-dom(ib)%w(ii-1,j,k))/16.0)* &
     &((-dom(ib)%u(ii,j,k+2)+9.0*dom(ib)%u(ii,j,k+1)+ &
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j,k-1))/16.0)
        ii=i-1; sm1=((-dom(ib)%w(ii+2,j,k)+9.0*dom(ib)%w(ii+1,j,k)+ &
     & 9.0*dom(ib)%w(ii,j,k)-dom(ib)%w(ii-1,j,k))/16.0)* &
     &((-dom(ib)%u(ii,j,k+2)+9.0*dom(ib)%u(ii,j,k+1)+ &
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j,k-1))/16.0)
        ii=i-2; sm2=((-dom(ib)%w(ii+2,j,k)+9.0*dom(ib)%w(ii+1,j,k)+ &
     & 9.0*dom(ib)%w(ii,j,k)-dom(ib)%w(ii-1,j,k))/16.0)* &
     &((-dom(ib)%u(ii,j,k+2)+9.0*dom(ib)%u(ii,j,k+1)+ &
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j,k-1))/16.0)
                 duwdx=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dx)

        jj=j+1; sp1=((-dom(ib)%w(i,jj+2,k)+9.0*dom(ib)%w(i,jj+1,k)+ &
     & 9.0*dom(ib)%w(i,jj,k)-dom(ib)%w(i,jj-1,k))/16.0)* &
     &((-dom(ib)%v(i,jj,k+2)+9.0*dom(ib)%v(i,jj,k+1)+ &
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i,jj,k-1))/16.0)
        jj=j;   sn=((-dom(ib)%w(i,jj+2,k)+9.0*dom(ib)%w(i,jj+1,k)+ &
     & 9.0*dom(ib)%w(i,jj,k)-dom(ib)%w(i,jj-1,k))/16.0)* &
     &((-dom(ib)%v(i,jj,k+2)+9.0*dom(ib)%v(i,jj,k+1)+ &
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i,jj,k-1))/16.0)
        jj=j-1; sm1=((-dom(ib)%w(i,jj+2,k)+9.0*dom(ib)%w(i,jj+1,k)+ &
     & 9.0*dom(ib)%w(i,jj,k)-dom(ib)%w(i,jj-1,k))/16.0)* &
     &((-dom(ib)%v(i,jj,k+2)+9.0*dom(ib)%v(i,jj,k+1)+ &
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i,jj,k-1))/16.0)
        jj=j-2; sm2=((-dom(ib)%w(i,jj+2,k)+9.0*dom(ib)%w(i,jj+1,k)+ &
     & 9.0*dom(ib)%w(i,jj,k)-dom(ib)%w(i,jj-1,k))/16.0)* &
     &((-dom(ib)%v(i,jj,k+2)+9.0*dom(ib)%v(i,jj,k+1)+ &
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i,jj,k-1))/16.0)
                 dvwdy=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dy)

                 dom(ib)%wstar(i,j,k)=(dom(ib)%woo(i,j,k)- &
     & dt*alfapr*(duwdx+dvwdy+dw2dz))

                 end do
              end do
           end do

        end do


      RETURN
      END SUBROUTINE rungek_conv4th
!#######################################################################
      SUBROUTINE rungek_convWENO
!#######################################################################
      use vars
      use mpi
      use multidata
      implicit none

      INTEGER :: i,j,k,ib
      DOUBLE PRECISION :: du2dx,dv2dy,dw2dz
      DOUBLE PRECISION :: duvdx,duvdy,duwdx,duwdz,dvwdy,dvwdz
      DOUBLE PRECISION :: uijk,vijk,wijk,dudx,dudy,dudz
      DOUBLE PRECISION :: dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        call boundu
        call HJ_WENO_dx(1)
        call HJ_WENO_dy(1)
        call HJ_WENO_dz(1)

        do ib=1,nbp

           do k=dom(ib)%ksu,dom(ib)%keu
              do i=dom(ib)%isu,dom(ib)%ieu
                 do j=dom(ib)%jsu,dom(ib)%jeu

                    if (dom(ib)%u(i,j,k).gt.0.0) then
                       dudx = dom(ib)%dphi_dxminus(i,j,k)
                    else if (dom(ib)%u(i,j,k).lt.0.0) then
                       dudx = dom(ib)%dphi_dxplus(i,j,k)
                    else
                       dudx = 0.0
                    end if

                    vijk=0.25*(dom(ib)%v(i,j,k)+dom(ib)%v(i+1,j,k)+ &
     &     dom(ib)%v(i,j-1,k)+dom(ib)%v(i+1,j-1,k))
                    if (vijk.gt.0.0) then
                       dudy = dom(ib)%dphi_dyminus(i,j,k)
                    else if (vijk.lt.0.0) then
                       dudy = dom(ib)%dphi_dyplus(i,j,k)
                    else
                       dudy = 0.0
                    end if

                    wijk=0.25*(dom(ib)%w(i,j,k)+dom(ib)%w(i+1,j,k)+ &
     &     dom(ib)%w(i,j,k-1)+dom(ib)%w(i+1,j,k-1))
                    if (wijk.gt.0.0) then
                       dudz = dom(ib)%dphi_dzminus(i,j,k)
                    else if (wijk.lt.0.0) then
                       dudz = dom(ib)%dphi_dzplus(i,j,k)
                    else
                       dudz = 0.0
                    end if

                    du2dx=dom(ib)%u(i,j,k)*dudx
                    duvdy=vijk*dudy
                    duwdz=wijk*dudz

                 dom(ib)%ustar(i,j,k)=(dom(ib)%uoo(i,j,k)- &
     & dt*alfapr*(du2dx+duvdy+duwdz))

                 end do
              end do
           end do
        end do

        call boundv
        call HJ_WENO_dx(2)
        call HJ_WENO_dy(2)
        call HJ_WENO_dz(2)

        do ib=1,nbp

           do k=dom(ib)%ksv,dom(ib)%kev
              do i=dom(ib)%isv,dom(ib)%iev
                 do j=dom(ib)%jsv,dom(ib)%jev

                    uijk=0.25*(dom(ib)%u(i,j,k)+dom(ib)%u(i,j+1,k)+ &
     &     dom(ib)%u(i-1,j,k)+dom(ib)%u(i-1,j+1,k))
                    if (uijk.gt.0.0) then
                       dvdx = dom(ib)%dphi_dxminus(i,j,k)
                    else if (uijk.lt.0.0) then
                       dvdx = dom(ib)%dphi_dxplus(i,j,k)
                    else
                       dvdx = 0.0
                    end if

                    if (dom(ib)%v(i,j,k).gt.0.0) then
                       dvdy = dom(ib)%dphi_dyminus(i,j,k)
                    else if (dom(ib)%v(i,j,k).lt.0.0) then
                       dvdy = dom(ib)%dphi_dyplus(i,j,k)
                    else
                       dvdy = 0.0
                    end if

                    wijk=0.25*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j+1,k)+ &
     &     dom(ib)%w(i,j,k-1)+dom(ib)%w(i,j+1,k-1))
                    if (wijk.gt.0.0) then
                       dvdz = dom(ib)%dphi_dzminus(i,j,k)
                    else if (wijk.lt.0.0) then
                       dvdz = dom(ib)%dphi_dzplus(i,j,k)
                    else
                       dvdz = 0.0
                    end if

                    duvdx=uijk*dvdx
                    dv2dy=dom(ib)%v(i,j,k)*dvdy
                    dvwdz=wijk*dvdz

                 dom(ib)%vstar(i,j,k)=(dom(ib)%voo(i,j,k)- &
     & dt*alfapr*(duvdx+dv2dy+dvwdz))

                 end do
              end do
           end do
        end do


        call boundw
        call HJ_WENO_dx(3)
        call HJ_WENO_dy(3)
        call HJ_WENO_dz(3)

        do ib=1,nbp

           do k=dom(ib)%ksw,dom(ib)%kew
              do i=dom(ib)%isw,dom(ib)%iew
                 do j=dom(ib)%jsw,dom(ib)%jew

                    uijk=0.25*(dom(ib)%u(i,j,k)+dom(ib)%u(i,j,k+1)+ &
     &     dom(ib)%u(i-1,j,k)+dom(ib)%u(i-1,j,k+1))
                    if (uijk.gt.0.0) then
                       dwdx = dom(ib)%dphi_dxminus(i,j,k)
                    else if (uijk.lt.0.0) then
                       dwdx = dom(ib)%dphi_dxplus(i,j,k)
                    else
                       dwdx = 0.0
                    end if

                    vijk=0.25*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j,k+1)+ &
     &     dom(ib)%v(i,j-1,k)+dom(ib)%v(i,j-1,k+1))
                    if (vijk.gt.0.0) then
                       dwdy = dom(ib)%dphi_dyminus(i,j,k)
                    else if (vijk.lt.0.0) then
                       dwdy = dom(ib)%dphi_dyplus(i,j,k)
                    else
                       dwdy = 0.0
                    end if

                    if (dom(ib)%w(i,j,k).gt.0.0) then
                       dwdz = dom(ib)%dphi_dzminus(i,j,k)
                    else if (dom(ib)%w(i,j,k).lt.0.0) then
                       dwdz = dom(ib)%dphi_dzplus(i,j,k)
                    else
                       dwdz = 0.0
                    end if

                    duwdx=uijk*dwdx
                    dvwdy=vijk*dwdy
                    dw2dz=dom(ib)%w(i,j,k)*dwdz

                 dom(ib)%wstar(i,j,k)=(dom(ib)%woo(i,j,k)- &
     & dt*alfapr*(duwdx+dvwdy+dw2dz))

                 end do
              end do
           end do

        end do

      RETURN
      END SUBROUTINE rungek_convWENO
!##########################################################################
