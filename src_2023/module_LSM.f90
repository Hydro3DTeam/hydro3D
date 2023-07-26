
!#######################################################################
      module lsm
!#######################################################################
      SAVE

      DOUBLE PRECISION :: dt_reinit,initial_vol,total_vol
      DOUBLE PRECISION,allocatable,dimension(:,:):: f,g,hj_phi,hj_phi1

      end module lsm
!#######################################################################
