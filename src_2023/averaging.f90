!#######################################################################
! SUBROUTINES: - update_mean
!              - add_noise
! FUNCTIONS:   - random_number_normal
!#######################################################################
      SUBROUTINE update_mean
!-----------------------------------------------------------------------
!     Calculates and stores the first and second time averaged variables
!     For the following field: 1.Velocity field: u,v,w
!                              2.Active and Passive Scalar field: S
!                              3.Viscosity field: vis
!#######################################################################
      use multidata
      use vars
      implicit none

      INTEGER :: i,j,k,ib,tti,ttj,ttk
!     DOUBLE PRECISION :: facp1,facm1,facp2,facm2
      DOUBLE PRECISION :: ufuf,vfvf,wfwf,ufvf,ufwf,vfwf
      DOUBLE PRECISION :: ucf,vcf,wcf,pfpf,TfTf,SfSf,SfUf,SfVf,SfWf
      DOUBLE PRECISION :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
      DOUBLE PRECISION :: uf_ip1,uf_im1,uf_jp1,uf_jm1,uf_kp1,uf_km1
      DOUBLE PRECISION :: vf_ip1,vf_im1,vf_jp1,vf_jm1,vf_kp1,vf_km1
      DOUBLE PRECISION :: wf_ip1,wf_im1,wf_jp1,wf_jm1,wf_kp1,wf_km1

      do ib=1,nbp
        tti=dom(ib)%ttc_i ; ttj=dom(ib)%ttc_j ; ttk=dom(ib)%ttc_k

!.....For first order moments
        if(ctime.ge.t_start_averaging1) then
	   IF(L_LSM) THEN
            do k=1,ttk ; do j=1,ttj ; do i=1,tti
              if(dom(ib)%phi(i,j,k).ge.0.0) then
		  dom(ib)%ntav1(i,j,k)=dom(ib)%ntav1(i,j,k)+1
                dom(ib)%facp1(i,j,k)=1./dom(ib)%ntav1(i,j,k)
        	  dom(ib)%facm1(i,j,k)=1.-dom(ib)%facp1(i,j,k)
		end if
            end do ; end do ; end do
	   ELSE
            do k=1,ttk ; do j=1,ttj ; do i=1,tti
		dom(ib)%ntav1(i,j,k)=dom(ib)%ntav1(i,j,k)+1
              dom(ib)%facp1(i,j,k)=1./dom(ib)%ntav1(i,j,k)
        	dom(ib)%facm1(i,j,k)=1.-dom(ib)%facp1(i,j,k)
            end do ; end do ; end do
	   END IF
        end if

!.....For second order moments
        if(ctime.ge.t_start_averaging2) then
	   IF(L_LSM) THEN
            do k=1,ttk ; do j=1,ttj ; do i=1,tti
		if(dom(ib)%phi(i,j,k).ge.0.0) then
		  dom(ib)%ntav2(i,j,k)=dom(ib)%ntav2(i,j,k)+1
                dom(ib)%facp2(i,j,k)=1./dom(ib)%ntav2(i,j,k)
        	  dom(ib)%facm2(i,j,k)=1.-dom(ib)%facp2(i,j,k)
		end if
            end do ; end do ; end do
	   ELSE
            do k=1,ttk ; do j=1,ttj ; do i=1,tti
		dom(ib)%ntav2(i,j,k)=dom(ib)%ntav2(i,j,k)+1
              dom(ib)%facp2(i,j,k)=1./dom(ib)%ntav2(i,j,k)
        	dom(ib)%facm2(i,j,k)=1.-dom(ib)%facp2(i,j,k)
            end do ; end do ; end do
	   END IF
        endif




        ! Mean U and Mean Variance UUM: Streamwise
	 IF(L_LSM) THEN
          do k=1,ttk ; do j=1,ttj ; do i=1,tti
	     if(dom(ib)%phi(i,j,k).ge.0.0) then
              dom(ib)%um(i,j,k)=dom(ib)%facm1(i,j,k)* &
     & dom(ib)%um(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%u(i,j,k)

              ufuf=(dom(ib)%u(i,j,k)-dom(ib)%um(i,j,k))* &
     & (dom(ib)%u(i,j,k)-dom(ib)%um(i,j,k))

              dom(ib)%uum(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%uum(i,j,k)+dom(ib)%facp2(i,j,k) * ufuf
	     end if
          end do ; end do ; end do

	 ELSE
          do k=1,ttk ; do j=1,ttj ; do i=1,tti
            dom(ib)%um(i,j,k)=dom(ib)%facm1(i,j,k)* &
     & dom(ib)%um(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%u(i,j,k)

            ufuf = (dom(ib)%u(i,j,k)-dom(ib)%um(i,j,k))* &
     & (dom(ib)%u(i,j,k)-dom(ib)%um(i,j,k))

            dom(ib)%uum(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%uum(i,j,k)+dom(ib)%facp2(i,j,k) * ufuf
          end do ; end do ; end do
	 END IF




        ! Mean V and Mean Variance VVM: Transversal
	 IF(L_LSM) THEN
          do k=1,ttk ; do j=1,ttj ; do i=1,tti
	     if(dom(ib)%phi(i,j,k).ge.0.0) then
              dom(ib)%vm(i,j,k)=dom(ib)%facm1(i,j,k)* &
     & dom(ib)%vm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%v(i,j,k)

              vfvf=(dom(ib)%v(i,j,k)-dom(ib)%vm(i,j,k))* &
     & (dom(ib)%v(i,j,k)-dom(ib)%vm(i,j,k))

              dom(ib)%vvm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%vvm(i,j,k)+dom(ib)%facp2(i,j,k)*vfvf
	     end if
          end do ; end do ; end do

	 ELSE
          do k=1,ttk ; do j=1,ttj ; do i=1,tti
            dom(ib)%vm(i,j,k)=dom(ib)%facm1(i,j,k)* &
     & dom(ib)%vm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%v(i,j,k)

            vfvf=(dom(ib)%v(i,j,k)-dom(ib)%vm(i,j,k))* &
     & (dom(ib)%v(i,j,k)-dom(ib)%vm(i,j,k))

            dom(ib)%vvm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%vvm(i,j,k)+dom(ib)%facp2(i,j,k)*vfvf
          end do ; end do ; end do
	 END IF




        ! Mean W and Mean Variance WWM: Vertical
	 IF(L_LSM) THEN
          do k=1,ttk ; do j=1,ttj ; do i=1,tti
	     if(dom(ib)%phi(i,j,k).ge.0.0) then
              dom(ib)%wm(i,j,k)=dom(ib)%facm1(i,j,k)* &
     & dom(ib)%wm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%w(i,j,k)

              wfwf=(dom(ib)%w(i,j,k)-dom(ib)%wm(i,j,k))* &
     & (dom(ib)%w(i,j,k)-dom(ib)%wm(i,j,k))

              dom(ib)%wwm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%wwm(i,j,k)+dom(ib)%facp2(i,j,k)*wfwf
	     end if
          end do ; end do ; end do

	 ELSE
          do k=1,ttk ; do j=1,ttj ; do i=1,tti
            dom(ib)%wm(i,j,k)=dom(ib)%facm1(i,j,k)* &
     & dom(ib)%wm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%w(i,j,k)

            wfwf=(dom(ib)%w(i,j,k)-dom(ib)%wm(i,j,k))* &
     & (dom(ib)%w(i,j,k)-dom(ib)%wm(i,j,k))

            dom(ib)%wwm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%wwm(i,j,k)+dom(ib)%facp2(i,j,k)*wfwf
          end do ; end do ; end do
	 END IF

        


        ! Mean P and Mean Variance PPM: Pressure
	 IF(L_LSM) THEN
          do k=1,ttk ; do j=1,ttj ; do i=1,tti
	     if(dom(ib)%phi(i,j,k).ge.0.0) then
              dom(ib)%pm(i,j,k)=dom(ib)%facm1(i,j,k)* &
     & dom(ib)%pm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%p(i,j,k)

              pfpf=(dom(ib)%p(i,j,k)-dom(ib)%pm(i,j,k))* &
     & (dom(ib)%p(i,j,k)-dom(ib)%pm(i,j,k))

              dom(ib)%ppm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%ppm(i,j,k)+dom(ib)%facp2(i,j,k) * pfpf
            end if
          end do ; end do ; end do

	 ELSE
          do k=1,ttk ; do j=1,ttj ; do i=1,tti
              dom(ib)%pm(i,j,k)=dom(ib)%facm1(i,j,k)* &
     & dom(ib)%pm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%p(i,j,k)

              pfpf=(dom(ib)%p(i,j,k)-dom(ib)%pm(i,j,k))* &
     & (dom(ib)%p(i,j,k)-dom(ib)%pm(i,j,k))

              dom(ib)%ppm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%ppm(i,j,k)+dom(ib)%facp2(i,j,k) * pfpf
          end do ; end do ; end do
        END IF



        ! Mean Vis: Kinetic viscosity
        do k=1,ttk ; do j=1,ttj ; do i=1,tti
           dom(ib)%vism(i,j,k)=dom(ib)%facm1(i,j,k)*dom(ib)%vism(i,j,k)+ &
     &  dom(ib)%facp1(i,j,k)*dom(ib)%vis(i,j,k)		
        end do ; end do ; end do



        ! Mean T and Mean Variance TTM: Temperature
        if(LENERGY) then
          IF(L_LSM) THEN
            do k=1,ttk ; do j=1,ttj ; do i=1,tti
              if(dom(ib)%phi(i,j,k).ge.0.0) then
                dom(ib)%Tm(i,j,k)=dom(ib)%facm1(i,j,k)* &
     & dom(ib)%Tm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%T(i,j,k)

                TfTf=(dom(ib)%T(i,j,k)-dom(ib)%Tm(i,j,k))* &
     & (dom(ib)%T(i,j,k)-dom(ib)%Tm(i,j,k))

                dom(ib)%Ttm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%Ttm(i,j,k)+dom(ib)%facp2(i,j,k)*TfTf
              end if
            end do ; end do ; end do

          ELSE
            do k=1,ttk ; do j=1,ttj ; do i=1,tti
              dom(ib)%Tm(i,j,k)=dom(ib)%facm1(i,j,k)* &
     & dom(ib)%Tm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%T(i,j,k)

              TfTf=(dom(ib)%T(i,j,k)-dom(ib)%Tm(i,j,k))* &
     & (dom(ib)%T(i,j,k)-dom(ib)%Tm(i,j,k))

              dom(ib)%Ttm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%Ttm(i,j,k)+dom(ib)%facp2(i,j,k)*TfTf
            end do ; end do ; end do
          END IF
        end if ! LENERGY


        ! Mean Phi: Free-Surface
	 if(L_LSM) then
          do k=1,ttk ; do j=1,ttj ; do i=1,tti
	     dom(ib)%phim(i,j,k)=dom(ib)%facm1(i,j,k)* &
     & dom(ib)%phim(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%phi(i,j,k)
          end do ; end do ; end do
        end if

        
        ! Mean S and Mean Variance SSM: Temperature
        if(LSCALAR) then
	   IF(L_LSM) THEN
            do k=1,ttk ; do j=1,ttj ; do i=1,tti
	       if(dom(ib)%phi(i,j,k).ge.0.0) then
                dom(ib)%Sm(i,j,k)=dom(ib)%facm1(i,j,k)* &
     & dom(ib)%Sm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%S(i,j,k)
                SfSf=(dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))* &
     & (dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))
                dom(ib)%Stm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%Stm(i,j,k)+dom(ib)%facp2(i,j,k)*SfSf

                SfUf=(dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))* &
     & (dom(ib)%U(i,j,k)-dom(ib)%Um(i,j,k))
                dom(ib)%SUtm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%SUtm(i,j,k)+dom(ib)%facp2(i,j,k)*SfUf

                SfVf=(dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))* &
     & (dom(ib)%V(i,j,k)-dom(ib)%Vm(i,j,k))
                dom(ib)%SVtm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%SVtm(i,j,k)+dom(ib)%facp2(i,j,k)*SfVf

                SfWf=(dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))* &
     & (dom(ib)%W(i,j,k)-dom(ib)%Wm(i,j,k))
                dom(ib)%SWtm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%SWtm(i,j,k)+dom(ib)%facp2(i,j,k)*SfWf
              end if
            end do ; end do ; end do

          ELSE
            do k=1,ttk ; do j=1,ttj ; do i=1,tti
              dom(ib)%Sm(i,j,k)=dom(ib)%facm1(i,j,k)* &
     & dom(ib)%Sm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%S(i,j,k)
              SfSf=(dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))* &
     & (dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))
              dom(ib)%Stm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%Stm(i,j,k)+dom(ib)%facp2(i,j,k)*SfSf

              SfUf=(dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))* &
     & (dom(ib)%U(i,j,k)-dom(ib)%Um(i,j,k))
              dom(ib)%SUtm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%SUtm(i,j,k)+dom(ib)%facp2(i,j,k)*SfUf

              SfVf=(dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))* &
     & (dom(ib)%V(i,j,k)-dom(ib)%Vm(i,j,k))
              dom(ib)%SVtm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%SVtm(i,j,k)+dom(ib)%facp2(i,j,k)*SfVf

              SfWf=(dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))* &
     & (dom(ib)%W(i,j,k)-dom(ib)%Wm(i,j,k))
              dom(ib)%SWtm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%SWtm(i,j,k)+dom(ib)%facp2(i,j,k)*SfWf
            end do ; end do ; end do
          END IF
        end if ! LSCALAR
        



        ! Eps and Mean Eps: 
        do k=2,ttk-1 ; do j=2,ttj-1 ; do i=2,tti-1
          uf_ip1 = dom(ib)%um(i+1,j,k)-dom(ib)%u(i+1,j,k)
          uf_im1 = dom(ib)%um(i-1,j,k)-dom(ib)%u(i-1,j,k)
          uf_jp1 = dom(ib)%um(i,j+1,k)-dom(ib)%u(i,j+1,k)
          uf_jm1 = dom(ib)%um(i,j-1,k)-dom(ib)%u(i,j-1,k)
          uf_kp1 = dom(ib)%um(i,j,k+1)-dom(ib)%u(i,j,k+1)
          uf_km1 = dom(ib)%um(i,j,k-1)-dom(ib)%u(i,j,k-1)

          vf_ip1 = dom(ib)%vm(i+1,j,k)-dom(ib)%v(i+1,j,k)
          vf_im1 = dom(ib)%vm(i-1,j,k)-dom(ib)%v(i-1,j,k)
          vf_jp1 = dom(ib)%vm(i,j+1,k)-dom(ib)%v(i,j+1,k)
          vf_jm1 = dom(ib)%vm(i,j-1,k)-dom(ib)%v(i,j-1,k)
          vf_kp1 = dom(ib)%vm(i,j,k+1)-dom(ib)%v(i,j,k+1)
          vf_km1 = dom(ib)%vm(i,j,k-1)-dom(ib)%v(i,j,k-1)

          wf_ip1 = dom(ib)%wm(i+1,j,k)-dom(ib)%w(i+1,j,k)
          wf_im1 = dom(ib)%wm(i-1,j,k)-dom(ib)%w(i-1,j,k)
          wf_jp1 = dom(ib)%wm(i,j+1,k)-dom(ib)%w(i,j+1,k)
          wf_jm1 = dom(ib)%wm(i,j-1,k)-dom(ib)%w(i,j-1,k)
          wf_kp1 = dom(ib)%wm(i,j,k+1)-dom(ib)%w(i,j,k+1)
          wf_km1 = dom(ib)%wm(i,j,k-1)-dom(ib)%w(i,j,k-1)

	   dudx = (uf_ip1-uf_im1)/(2*dom(ib)%dx)
	   dudy = (uf_jp1-uf_jm1)/(2*dom(ib)%dy)
	   dudz = (uf_kp1-uf_km1)/(2*dom(ib)%dz)

	   dvdx = (vf_ip1-vf_im1)/(2*dom(ib)%dx)
	   dvdy = (vf_jp1-vf_jm1)/(2*dom(ib)%dy)
	   dvdz = (vf_kp1-vf_km1)/(2*dom(ib)%dz)

	   dwdx = (wf_ip1-wf_im1)/(2*dom(ib)%dx)
	   dwdy = (wf_jp1-wf_jm1)/(2*dom(ib)%dy)
	   dwdz = (wf_kp1-wf_km1)/(2*dom(ib)%dz)

        dom(ib)%eps(i,j,k) = dom(ib)%vis(i,j,k) * &
     &  (2.d0*dudx**2+dvdx**2+dwdx**2 &
     &   +dudy**2+2.d0*dvdy**2+dwdy**2 &
     &   +dudz**2+dvdz**2+2.d0*dwdz**2 &
     &   +2.d0*dudy*dvdx+2.d0*dudz*dwdx+2.d0*dvdz*dwdy)

          dom(ib)%epsm(i,j,k)=dom(ib)%facm1(i,j,k)*dom(ib)%epsm(i,j,k)+ &
     &  dom(ib)%facp1(i,j,k)*dom(ib)%eps(i,j,k)					
        end do ; end do ; end do
        
        ! Mean UVM,UWM,VWM: Reynolds Shear Stress
        IF(L_LSM) THEN
          do k=2,ttk-1 ; do j=2,ttj-1 ; do i=2,tti-1
	     if(dom(ib)%phi(i,j,k).ge.0.0) then
              ucf=0.5*((dom(ib)%u(i-1,j,k)-dom(ib)%um(i-1,j,k))+ &
     & (dom(ib)%u(i,j,k)-dom(ib)%um(i,j,k)))
              vcf=0.5*((dom(ib)%v(i,j-1,k)-dom(ib)%vm(i,j-1,k))+ &
     & (dom(ib)%v(i,j,k)-dom(ib)%vm(i,j,k)))
              wcf=0.5*((dom(ib)%w(i,j,k-1)-dom(ib)%wm(i,j,k-1))+ &
     & (dom(ib)%w(i,j,k)-dom(ib)%wm(i,j,k)))

              ufvf = ucf * vcf
              ufwf = ucf * wcf
              vfwf = vcf * wcf

              dom(ib)%uvm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%uvm(i,j,k)+dom(ib)%facp2(i,j,k)*ufvf
              dom(ib)%uwm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%uwm(i,j,k)+dom(ib)%facp2(i,j,k)*ufwf
              dom(ib)%vwm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%vwm(i,j,k)+dom(ib)%facp2(i,j,k)*vfwf
	     end if
          end do ; end do ; end do
	 ELSE
          do k=2,ttk-1 ; do j=2,ttj-1 ; do i=2,tti-1
            ucf=0.5*((dom(ib)%u(i-1,j,k)-dom(ib)%um(i-1,j,k))+ &
     & (dom(ib)%u(i,j,k)-dom(ib)%um(i,j,k)))
            vcf=0.5*((dom(ib)%v(i,j-1,k)-dom(ib)%vm(i,j-1,k))+ &
     & (dom(ib)%v(i,j,k)-dom(ib)%vm(i,j,k)))
            wcf=0.5*((dom(ib)%w(i,j,k-1)-dom(ib)%wm(i,j,k-1))+ &
     & (dom(ib)%w(i,j,k)-dom(ib)%wm(i,j,k)))

            ufvf = ucf * vcf
            ufwf = ucf * wcf
            vfwf = vcf * wcf

            dom(ib)%uvm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%uvm(i,j,k)+dom(ib)%facp2(i,j,k)*ufvf
            dom(ib)%uwm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%uwm(i,j,k)+dom(ib)%facp2(i,j,k)*ufwf
            dom(ib)%vwm(i,j,k)=dom(ib)%facm2(i,j,k)* &
     & dom(ib)%vwm(i,j,k)+dom(ib)%facp2(i,j,k)*vfwf
          end do ; end do ; end do
        END IF

      END DO ! ib

      RETURN
      END SUBROUTINE update_mean
!#######################################################################
      SUBROUTINE add_noise(fnoise)
!-----------------------------------------------------------------------
!     When noise > 0.0 (control.cin)
!     Add noise in the INITIAL velocity field only at the beginning of
!     the simulation. It calls the FUNCTION random_number
!#######################################################################
      use vars
      use multidata
      implicit none

      INTEGER :: i,j,k,ib
      DOUBLE PRECISION :: fnoise
      DOUBLE PRECISION :: random_number_normal

!......add some gausian noise to flowfield

      call RANDOM_SEED

      DO ib=1,nbp
        
        ! Add noise to velocity field
        IF(L_LSM) THEN
          do k=1,dom(ib)%ttc_k ; do j=1,dom(ib)%ttc_j ; do i=1,dom(ib)%ttc_i
	     if(dom(ib)%phi(i,j,k) .ge. 0.0) then
              dom(ib)%u(i,j,k) = dom(ib)%u(i,j,k) + random_number_normal(0.0,fnoise)
              dom(ib)%v(i,j,k) = dom(ib)%v(i,j,k) + random_number_normal(0.0,fnoise)
              dom(ib)%w(i,j,k) = dom(ib)%w(i,j,k) + random_number_normal(0.0,fnoise)
	     endif
          end do ; end do ; end do

        ELSE IF(L_LSMbase) THEN
          do k=1,dom(ib)%ttc_k ; do j=1,dom(ib)%ttc_j ; do i=1,dom(ib)%ttc_i
            if(dom(ib)%zc(k).le.length) then
              dom(ib)%u(i,j,k) = dom(ib)%u(i,j,k) + random_number_normal(0.0,fnoise)
              dom(ib)%v(i,j,k) = dom(ib)%v(i,j,k) + random_number_normal(0.0,fnoise)
              dom(ib)%w(i,j,k) = dom(ib)%w(i,j,k) + random_number_normal(0.0,fnoise)
            end if
          end do ; end do ; end do

        ELSE 
          do k=1,dom(ib)%ttc_k ; do j=1,dom(ib)%ttc_j ; do i=1,dom(ib)%ttc_i
            dom(ib)%u(i,j,k) = dom(ib)%u(i,j,k) + random_number_normal(0.0,fnoise)
            dom(ib)%v(i,j,k) = dom(ib)%v(i,j,k) + random_number_normal(0.0,fnoise)
            dom(ib)%w(i,j,k) = dom(ib)%w(i,j,k) + random_number_normal(0.0,fnoise)
          end do ; end do ; end do
        END IF 

      END DO ! ib

      RETURN
      END SUBROUTINE add_noise
!#######################################################################
        FUNCTION random_number_normal(mean,sigma) result( fn_val )
!-----------------------------------------------------------------------
!       Generate random numbers
!       with a normal distribution with given mean and standard deviaton.
!
!       Generate a random normal deviate using the polar method.
!       Reference:
!         Marsaglia,G. & Bray,T.A. 'A convenient method for generating
!         normal variables', Siam Rev., vol.6, 260-264, 1964.
!#######################################################################
        implicit none

        DOUBLE PRECISION :: fn_val,mean,sigma
        DOUBLE PRECISION :: ull, sumall
        DOUBLE PRECISION, save :: vll, sln
        DOUBLE PRECISION, parameter :: one = 1.0, vsmall = tiny( one )
        LOGICAL, save   :: second = .false.

        if (second) then

!...... If second, use the second random number generated on last call
           second = .false.
           fn_val = vll*sln

        else
!...... First call; generate a pair of random normals
           second = .true.
           do
              call random_number(ull)
              call random_number(vll)

              ull = scale( ull, 1 ) - one
              vll = scale( vll, 1 ) - one

!.........vsmall added to prevent LOG(zero) / zero
              sumall = ull*ull + vll*vll + vsmall
              if(sumall < one) exit
           end do

           sln = sqrt(- scale( log(sumall), 1 ) / sumall)
           fn_val = ull*sln
        end if

!.....set mean and standart deviation
        fn_val = fn_val * sigma + mean

        RETURN
        END FUNCTION random_number_normal
!#######################################################################
