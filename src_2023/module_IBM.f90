!######################################################################
      module imb
!######################################################################
      SAVE

      INTEGER :: imp_proc_master,imb_block_master,forcefilej
      INTEGER :: bodynum,maxnode,master,maxnodeIBS,mdfsteps,yangcase
      INTEGER :: ibm_out_forces,nIBslv,nl

	
      DOUBLE PRECISION :: lambda,sigma,dxm,dym,dzm,nxl
      DOUBLE PRECISION,allocatable :: Cx(:),Cxor(:),Cy(:),Cyor(:)
      DOUBLE PRECISION,allocatable :: Cz(:),Czor(:),l2norm(:)			
      DOUBLE PRECISION,allocatable,dimension(:) :: xaero,yaero,zaero
      DOUBLE PRECISION,allocatable,dimension(:) :: radsin,rads,pitch
      DOUBLE PRECISION,allocatable,dimension(:) :: R,reddelta
      DOUBLE PRECISION,allocatable,dimension(:) :: nodex_loc,nodey_loc
      DOUBLE PRECISION,allocatable,dimension(:) :: nodez_loc
      DOUBLE PRECISION,allocatable,dimension(:) :: FX1_loc,FX2_loc
      DOUBLE PRECISION,allocatable,dimension(:) :: FX3_loc
      DOUBLE PRECISION,allocatable,dimension(:) :: R0_loc,alpha0_loc
      DOUBLE PRECISION,allocatable,dimension(:) :: U_Beta1_loc
      DOUBLE PRECISION,allocatable,dimension(:) :: U_Beta2_loc
      DOUBLE PRECISION,allocatable,dimension(:) :: U_Beta3_loc,zini,zend
      DOUBLE PRECISION,allocatable,dimension(:,:) :: dh1_loc,dh2_loc
      DOUBLE PRECISION,allocatable,dimension(:,:) :: dh3_loc
      DOUBLE PRECISION,allocatable,dimension(:,:) :: nodexlocal
      DOUBLE PRECISION,allocatable,dimension(:,:) :: nodeylocal
      DOUBLE PRECISION,allocatable,dimension(:,:) :: nodezlocal
      DOUBLE PRECISION,allocatable,dimension(:) :: FX1NF,FX2NF,FX3NF
      DOUBLE PRECISION,allocatable,dimension(:,:) :: FX1,FX2,FX3
      DOUBLE PRECISION,allocatable,dimension(:,:) :: FX1M,FX2M,FX3M
      DOUBLE PRECISION,allocatable,dimension(:,:) :: nodex,nodey,nodez
      DOUBLE PRECISION,allocatable,dimension(:,:) :: U_Beta1,U_Beta2
      DOUBLE PRECISION,allocatable,dimension(:,:) :: U_Beta3
      DOUBLE PRECISION,allocatable,dimension(:,:) :: alpha0,R0

      INTEGER,allocatable,dimension(:,:) :: I_nr_V,J_nr_V,K_nr_V
      INTEGER,allocatable,dimension(:,:) :: I_nr_U,J_nr_U,K_nr_U
      INTEGER,allocatable,dimension(:,:) :: I_nr_W,J_nr_W,K_nr_W
      INTEGER,allocatable,dimension(:) :: imbnumber,imb_shape
      INTEGER,allocatable,dimension(:) :: kmaxU,kmaxV,kmaxW,nodes
      INTEGER,allocatable,dimension(:) :: imb_block,turax	
      INTEGER,allocatable,dimension(:) :: lag_bod_loc,cmax,linfin
      INTEGER,allocatable,dimension(:) :: imb_block_loc,axis
      INTEGER,allocatable,dimension(:) :: imbinblock_loc,rott_loc
      INTEGER,allocatable,dimension(:) :: nscatter	

      CHARACTER(32),allocatable,dimension (:) :: filepoints
      LOGICAL,allocatable,dimension(:) :: rotating,ibturbine

      INTEGER,allocatable,dimension(:) :: numIBslave,Lslave,Lslv,dsplc
      INTEGER,allocatable,dimension(:) :: imb_proc_loc

      DOUBLE PRECISION,allocatable,dimension(:) :: xdh_U,ydh_U,zdh_U
      DOUBLE PRECISION,allocatable,dimension(:) :: xdh_V,ydh_V,zdh_V
      DOUBLE PRECISION,allocatable,dimension(:) :: xdh_W,ydh_W,zdh_W

      INTEGER,allocatable,dimension(:,:) :: Irec_U,Jrec_U,Krec_U
      INTEGER,allocatable,dimension(:,:) :: Irec_V,Jrec_V,Krec_V
      INTEGER,allocatable,dimension(:,:) :: Irec_W,Jrec_W,Krec_W

      DOUBLE PRECISION,allocatable,dimension(:,:) :: USTAR_mas,VSTAR_mas
      DOUBLE PRECISION,allocatable,dimension(:,:) :: WSTAR_mas

      DOUBLE PRECISION :: Cl_act,Cd_act
      INTEGER,allocatable,dimension(:) :: ppc_act
      DOUBLE PRECISION,allocatable,dimension(:) :: r_act,c_act,Pit_act
      DOUBLE PRECISION,allocatable,dimension(:,:) :: F_X_MEAN,F_Q_MEAN
      DOUBLE PRECISION,allocatable,dimension(:,:) :: F_Y_MEAN,F_Z_MEAN

      end module imb
