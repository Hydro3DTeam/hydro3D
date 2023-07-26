!#######################################################################
      module lpt
!#######################################################################
      SAVE

      INTEGER :: np,tsteps_pt,tsnr,cnt_pt,np_loc,ptnr
      DOUBLE PRECISION :: xp,yp,zp,uop,vop,wop,div,Dp,sigma
      LOGICAL :: DF,PSIcell,random,random_dp	
      INTEGER,allocatable,dimension(:)::	ptsinproc
      INTEGER,allocatable,dimension(:)::  id
      DOUBLE PRECISION,allocatable,dimension(:) :: xp_pt,yp_pt,zp_pt
      DOUBLE PRECISION,allocatable,dimension(:) :: xp_loc,yp_loc,zp_loc
      DOUBLE PRECISION,allocatable,dimension(:) :: uop_pt,vop_pt,wop_pt
      DOUBLE PRECISION,allocatable,dimension(:) :: uop_loc,vop_loc
      DOUBLE PRECISION,allocatable,dimension(:) :: wop_loc,dp_loc
      DOUBLE PRECISION,allocatable,dimension(:) :: Fu,Fv,Fw
      DOUBLE PRECISION,allocatable,dimension(:) :: Fpu,Fpv,Fpw
      DOUBLE PRECISION,allocatable,dimension(:) :: xpold,ypold,zpold
      DOUBLE PRECISION,allocatable,dimension(:) :: uopold,vopold,wopold
      DOUBLE PRECISION,allocatable,dimension(:) :: dp_pt,dp_old

      end module
!#######################################################################
