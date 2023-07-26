!#######################################################################
! SUBROUTINES: - calmas
!              - calvel
!              - pbound
!              - pressure_forcing
!#######################################################################
      SUBROUTINE calmas
!-----------------------------------------------------------------------
!     Calculates divergence free condition (su,sup) used in the poisson
!     pressure solver. Also calculates the residual error.
!#######################################################################
      use vars
      use mpi
      use multidata
      implicit none

      INTEGER :: i,j,k,ib
      DOUBLE PRECISION :: fact,fact1,buffer_rmax,ddt
      DOUBLE PRECISION,pointer,dimension(:,:,:) :: fi


      rmax=0.d0

      IF(differencing.eq.2) THEN !4th CDS
        DO ib=1,nbp
          if(L_LSM) then
            fi => dom(ib)%sup 
            ddt=dom(ib)%dx*dom(ib)%dy*dom(ib)%dz/dt
          else
            fi => dom(ib)%su
            ddt=1/dt
          end if

          fact=0.d0 ; fact1=0.d0
          do k=dom(ib)%ksp,dom(ib)%kep
            do i=dom(ib)%isp,dom(ib)%iep
              do j=dom(ib)%jsp,dom(ib)%jep

                fact=( &
     &(-dom(ib)%u(i+1,j,k)+27.0*dom(ib)%u(i,j,k)- &
     &27.0*dom(ib)%u(i-1,j,k)+dom(ib)%u(i-2,j,k))/(24.0*dom(ib)%dx)+ &
     &(-dom(ib)%v(i,j+1,k)+27.0*dom(ib)%v(i,j,k)- &
     &27.0*dom(ib)%v(i,j-1,k)+dom(ib)%v(i,j-2,k))/(24.0*dom(ib)%dy)+ &
     &(-dom(ib)%w(i,j,k+1)+27.0*dom(ib)%w(i,j,k)- &
     &27.0*dom(ib)%w(i,j,k-1)+dom(ib)%w(i,j,k-2))/(24.0*dom(ib)%dz))

                fi(i,j,k) = fact * ddt

                fact1=fact*dom(ib)%dx*dom(ib)%dy*dom(ib)%dz
                resor=max(rmax,abs(fact1))
                rmax=max(resor,abs(fact1))
              end do
            end do
          end do
        END DO

      ELSE !2nd CDS or WENO:
        DO ib=1,nbp
          if(L_LSM) then
            fi => dom(ib)%sup 
            ddt=dom(ib)%dx*dom(ib)%dy*dom(ib)%dz/dt
          else
            fi => dom(ib)%su
            ddt=1/dt
          end if

          fact=0.d0 ;  fact1=0.d0
          do k=dom(ib)%ksp,dom(ib)%kep
            do i=dom(ib)%isp,dom(ib)%iep
              do j=dom(ib)%jsp,dom(ib)%jep

                fact=( &
     & (dom(ib)%u(i,j,k)-dom(ib)%u(i-1,j,k))/dom(ib)%dx+ &
     & (dom(ib)%v(i,j,k)-dom(ib)%v(i,j-1,k))/dom(ib)%dy+ &
     & (dom(ib)%w(i,j,k)-dom(ib)%w(i,j,k-1))/dom(ib)%dz)

                fi(i,j,k) = fact * ddt

                fact1=fact*dom(ib)%dx*dom(ib)%dy*dom(ib)%dz
                resor=max(rmax,abs(fact1))
                rmax=max(resor,abs(fact1))
              end do
            end do
          end do
        END DO
      END IF

      buffer_rmax = rmax
      call MPI_ALLREDUCE(buffer_rmax,rmax,1,MPI_FLT,MPI_MAX,MPI_COMM_WORLD,ierr)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

      
      RETURN
      END SUBROUTINE calmas
!#######################################################################
      SUBROUTINE calvel
!-----------------------------------------------------------------------
!     As part of the fractional step method correct the free divergence 
!     with the pressure gradient.
!#######################################################################
      use vars
      use mpi
      use multidata
      implicit none

      INTEGER :: i,j,k,ib,is,ie,js,je,ks,ke


      IF(L_LSM) THEN

        DO ib=1,nbp
          ! U Velocity: Streamwise
          is=dom(ib)%isu ; ie=dom(ib)%ieu ; js=dom(ib)%jsu
          je=dom(ib)%jeu ; ks=dom(ib)%ksu ; ke=dom(ib)%keu

          do k=ks,ke ; do j=js,je ; do i=is,ie
            dom(ib)%u(i,j,k) = (dom(ib)%ustar(i,j,k) - dt * alfapr &
     &                       * (dom(ib)%p(i+1,j,k)-dom(ib)%p(i,j,k))/dom(ib)%dx &
     &                       / (0.5*(dom(ib)%dens(i,j,k)+dom(ib)%dens(i+1,j,k))))
	   end do ; end do ; end do 
  

          ! V Velocity: Transversal
          is=dom(ib)%isv ; ie=dom(ib)%iev ; js=dom(ib)%jsv
          je=dom(ib)%jev ; ks=dom(ib)%ksv ; ke=dom(ib)%kev

          do k=ks,ke ; do j=js,je ; do i=is,ie
            dom(ib)%v(i,j,k) = (dom(ib)%vstar(i,j,k) - dt * alfapr  &
     &                       * (dom(ib)%p(i,j+1,k)-dom(ib)%p(i,j,k))/dom(ib)%dy &
     &                       / (0.5*(dom(ib)%dens(i,j,k)+dom(ib)%dens(i,j+1,k))))
	   end do ; end do ; end do 

          
          ! W Velocity: Transversal
          is=dom(ib)%isw ; ie=dom(ib)%iew ; js=dom(ib)%jsw
          je=dom(ib)%jew ; ks=dom(ib)%ksw ; ke=dom(ib)%kew

          do k=ks,ke ; do j=js,je ; do i=is,ie
            dom(ib)%w(i,j,k) = (dom(ib)%wstar(i,j,k) - dt * alfapr  &
     &                       * (dom(ib)%p(i,j,k+1)-dom(ib)%p(i,j,k))/dom(ib)%dz &
     &                       / (0.5*(dom(ib)%dens(i,j,k)+dom(ib)%dens(i,j,k+1))))
	   end do ; end do ; end do 
        END DO ! ib


      ELSE !.not.L_LSM

        DO ib=1,nbp
          ! U Velocity: Streamwise
          is=dom(ib)%isu ; ie=dom(ib)%ieu ; js=dom(ib)%jsu
          je=dom(ib)%jeu ; ks=dom(ib)%ksu ; ke=dom(ib)%keu

          do k=ks,ke ; do j=js,je ; do i=is,ie
            dom(ib)%u(i,j,k) =  dom(ib)%ustar(i,j,k) - dt * alfapr &
     &                       * (dom(ib)%p(i+1,j,k)-dom(ib)%p(i,j,k))/dom(ib)%dx
	   end do ; end do ; end do 
  

          ! V Velocity: Transversal
          is=dom(ib)%isv ; ie=dom(ib)%iev ; js=dom(ib)%jsv
          je=dom(ib)%jev ; ks=dom(ib)%ksv ; ke=dom(ib)%kev

          do k=ks,ke ; do j=js,je ; do i=is,ie
            dom(ib)%v(i,j,k) =  dom(ib)%vstar(i,j,k) - dt * alfapr &
     &                       * (dom(ib)%p(i,j+1,k)-dom(ib)%p(i,j,k))/dom(ib)%dy
	   end do ; end do ; end do 


          ! W Velocity: Transversal
          is=dom(ib)%isw ; ie=dom(ib)%iew ; js=dom(ib)%jsw
          je=dom(ib)%jew ; ks=dom(ib)%ksw ; ke=dom(ib)%kew

          do k=ks,ke ; do j=js,je ; do i=is,ie
            dom(ib)%w(i,j,k) =  dom(ib)%wstar(i,j,k) - dt * alfapr  &
     &                       * (dom(ib)%p(i,j,k+1)-dom(ib)%p(i,j,k))/dom(ib)%dz
	   end do ; end do ; end do 
        END DO ! ib

      END IF ! L_LSM

      call exchange(1)
      call exchange(2)
      call exchange(3)
      
      RETURN
      END SUBROUTINE calvel
!#######################################################################
      SUBROUTINE pbound
!-----------------------------------------------------------------------
!     Transfer the pressure field to the ghost cells when the boundary condition
!     is not periodic. These values with then be exchanged between the 
!     neighbours of the subdomain. (exchange, exchangep) 
!#######################################################################
      use vars
      use multidata
      implicit none

      INTEGER :: i,j,k,ib,isp,iep,jsp,jep,ksp,kep

        do ib=1,nbp

           isp=dom(ib)%isp; iep=dom(ib)%iep
           jsp=dom(ib)%jsp; jep=dom(ib)%jep
           ksp=dom(ib)%ksp; kep=dom(ib)%kep
!...................... west and east............................
           if (dom(ib)%iprev.lt.0 .and. dom(ib)%bc_west.ne.5) then
              do k=ksp-1,kep+1
                 do j=jsp-1,jep+1
                    dom(ib)%p(isp-1,j,k)  =dom(ib)%p(isp,j,k)
                 end do
              end do
           end if
           if (dom(ib)%inext.lt.0 .and. dom(ib)%bc_east.ne.5) then
              do k=ksp-1,kep+1
                 do j=jsp-1,jep+1
                    dom(ib)%p(iep+1,j,k)  =dom(ib)%p(iep,j,k)
                 end do
              end do
           end if
!.....................south and north.........................
           if (dom(ib)%jprev.lt.0 .and. dom(ib)%bc_south.ne.5) then
              do k=ksp-1,kep+1
                 do i=isp-1,iep+1
                    dom(ib)%p(i,jsp-1,k)  =dom(ib)%p(i,jsp,k)
                 end do
              end do
           end if
           if (dom(ib)%jnext.lt.0 .and. dom(ib)%bc_north.ne.5) then
              do k=ksp-1,kep+1
                 do i=isp-1,iep+1
                    dom(ib)%p(i,jep+1,k)  =dom(ib)%p(i,jep,k)
                 end do
              end do
           end if
!.....................bottom and top.........................
           if (dom(ib)%kprev.lt.0 .and. dom(ib)%bc_bottom.ne.5) then
              do j=jsp-1,jep+1
                 do i=isp-1,iep+1
                    dom(ib)%p(i,j,ksp-1)  =dom(ib)%p(i,j,ksp)
                 end do
              end do
           end if
           if (dom(ib)%knext.lt.0 .and. dom(ib)%bc_top.ne.5) then
              do j=jsp-1,jep+1
                 do i=isp-1,iep+1
                    dom(ib)%p(i,j,kep+1)  =dom(ib)%p(i,j,kep)
                 end do
              end do
           end if

        end do
      
      RETURN
      END SUBROUTINE pbound
