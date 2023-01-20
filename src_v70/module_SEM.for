!#############################################################################
      module module_SEM
!#############################################################################
      SAVE
	
      CHARACTER(44) :: FILEGLOBAL
      INTEGER,allocatable,dimension(:) :: elemyst,elemyen
      INTEGER,allocatable,dimension(:) :: elemzst,elemzen
      INTEGER,allocatable,dimension(:) :: iddom,ljdom,lkdom
      DOUBLE PRECISION,allocatable,dimension(:) ::  Ksem
      DOUBLE PRECISION,allocatable,dimension(:,:,:) :: Vsem ,Usem
      DOUBLE PRECISION,allocatable,dimension(:,:) :: X_EDDY,EPSILO,MOLT
      DOUBLE PRECISION,allocatable,dimension(:,:,:) :: SIGMA
      
      end 
