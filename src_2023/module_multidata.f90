!##########################################################################
      module multidata
!##########################################################################
      SAVE

      INTEGER :: nbp,nbpmax,num_domains
      INTEGER :: rdivmax,idom,jdom,kdom
      INTEGER,allocatable,dimension(:) :: rdiv
      INTEGER,allocatable,dimension(:,:) :: rdv
      INTEGER,allocatable,dimension(:) :: dom_id,dom_ad,dom_indid
      INTEGER,allocatable,dimension(:) :: i_unst,j_unst,k_unst
      INTEGER,allocatable,dimension(:) :: id_unst
      INTEGER,allocatable,dimension(:) :: imbinblk	
      DOUBLE PRECISION,allocatable,dimension(:,:) :: xcor,ycor,zcor

      type multidom
        INTEGER :: inext,iprev,jnext,jprev,knext,kprev
        INTEGER :: corprev1,corprev2,corprev3,corprev4
        INTEGER :: cornext1,cornext2,cornext3,cornext4
        INTEGER :: edgprev1,edgprev2,edgprev3
        INTEGER :: edgprev4,edgprev5,edgprev6
        INTEGER :: edgnext1,edgnext2,edgnext3
        INTEGER :: edgnext4,edgnext5,edgnext6
        INTEGER :: per_in,per_ip,per_jn,per_jp,per_kn,per_kp
        INTEGER :: ttc_i,ttc_j,ttc_k,ttc_ijk,ngrid
        INTEGER :: isp,iep,jsp,jep,ksp,kep
        INTEGER :: isu,ieu,jsu,jeu,ksu,keu
	 INTEGER :: isv,iev,jsv,jev,ksv,kev
        INTEGER :: isw,iew,jsw,jew,ksw,kew
        INTEGER :: nwork,nvars
        INTEGER :: mximb
        INTEGER :: rq_m1,rq_p1,rq_m2,rq_p2,rq_m3,rq_p3
        INTEGER :: rq_c1m,rq_c2m,rq_c3m,rq_c4m
        INTEGER :: rq_c1p,rq_c2p,rq_c3p,rq_c4p
        INTEGER :: rq_e1m,rq_e2m,rq_e3m,rq_e4m,rq_e5m,rq_e6m
        INTEGER :: rq_e1p,rq_e2p,rq_e3p,rq_e4p,rq_e5p,rq_e6p
        INTEGER :: bc_west,bc_east,bc_south,bc_north,bc_bottom,bc_top
        INTEGER :: Tbc_west,Tbc_east,Tbc_south,Tbc_north
        INTEGER :: Tbc_bottom,Tbc_top
        DOUBLE PRECISION :: xsl,ysl,zsl,xel,yel,zel,dx,dy,dz
	 LOGICAL :: coarse_ng,fine_ng
        INTEGER,dimension(26) :: tg
        INTEGER,allocatable,dimension(:) :: faz,cntp
        INTEGER,allocatable,dimension(:,:,:) :: ntav1,ntav2
	 INTEGER,allocatable,dimension(:,:,:) :: ibfactor
        INTEGER,allocatable,dimension(:,:,:,:) :: ndimb,imbbodynum
	 DOUBLE PRECISION,allocatable,dimension(:) :: tauw
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: S,So,Sm,Stm
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: SUtm,SVtm,SWtm
	 DOUBLE PRECISION,allocatable,dimension(:,:,:) :: sfactor
        DOUBLE PRECISION,allocatable,dimension(:) :: x,y,z,xc,yc,zc
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: u,v,w,p,pp
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: ksgs,ksgso
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: eps,epso
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: T,To,Tm,Ttm
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: su,ap,sup
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: ae,aw,as,an
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: at,ab
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: um,vm,wm
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: pm,ppm
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: vis,vism,epsm
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: uum,vvm,wwm
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: uvm,uwm,vwm
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: ustar,vstar
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: wstar
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: uo,uoo,vo
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: voo,wo,woo
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: stfcinf
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: facp1,facp2
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: facm1,facm2
        DOUBLE PRECISION,allocatable,dimension(:) :: dh1,dh2,dh3
        DOUBLE PRECISION,allocatable,dimension(:) :: cof
        DOUBLE PRECISION,allocatable,dimension(:,:) :: tauwe,tauww
        DOUBLE PRECISION,allocatable,dimension(:,:) :: tauws,tauwn
        DOUBLE PRECISION,allocatable,dimension(:,:) :: tauwt,tauwb
        DOUBLE PRECISION,allocatable,dimension(:,:) :: tauwe2,tauww2
        DOUBLE PRECISION,allocatable,dimension(:,:) :: tauws2,tauwn2
        DOUBLE PRECISION,allocatable,dimension(:,:) :: tauwt2,tauwb2
        DOUBLE PRECISION,allocatable,dimension(:) :: sendb_m1,sendb_p1
        DOUBLE PRECISION,allocatable,dimension(:) :: recvb_m1,recvb_p1
        DOUBLE PRECISION,allocatable,dimension(:) :: sendb_m2,sendb_p2
        DOUBLE PRECISION,allocatable,dimension(:) :: recvb_m2,recvb_p2
        DOUBLE PRECISION,allocatable,dimension(:) :: sendb_m3,sendb_p3
        DOUBLE PRECISION,allocatable,dimension(:) :: recvb_m3,recvb_p3
        DOUBLE PRECISION,allocatable,dimension(:) :: sc1m,sc1p,rc1m,rc1p
        DOUBLE PRECISION,allocatable,dimension(:) :: sc2m,sc2p,rc2m,rc2p
        DOUBLE PRECISION,allocatable,dimension(:) :: sc3m,sc3p,rc3m,rc3p
        DOUBLE PRECISION,allocatable,dimension(:) :: sc4m,sc4p,rc4m,rc4p
        DOUBLE PRECISION,allocatable,dimension(:) :: se1m,se1p,re1m,re1p
        DOUBLE PRECISION,allocatable,dimension(:) :: se2m,se2p,re2m,re2p
        DOUBLE PRECISION,allocatable,dimension(:) :: se3m,se3p,re3m,re3p
        DOUBLE PRECISION,allocatable,dimension(:) :: se4m,se4p,re4m,re4p
        DOUBLE PRECISION,allocatable,dimension(:) :: se5m,se5p,re5m,re5p
        DOUBLE PRECISION,allocatable,dimension(:) :: se6m,se6p,re6m,re6p
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: d1
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: dphi_dxplus
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: dphi_dyplus
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: dphi_dzplus
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: dphi_dxminus
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: dphi_dyminus
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: dphi_dzminus
!============================== LSM VARIABLES =============================
	 INTEGER :: tot
        INTEGER :: niul,njul,nkul,nivl,njvl,nkvl,niwl,njwl,nkwl
        INTEGER :: nipl,njpl,nkpl,nigl,njgl,nkgl
        INTEGER :: nipl2,njpl2,nkpl2
        INTEGER,allocatable,dimension(:) :: ijkp_lsm
        REAL,allocatable,dimension(:) :: sendb_m,sendb_p
        REAL,allocatable,dimension(:) :: recvb_m,recvb_p
        DOUBLE PRECISION,allocatable,dimension(:) :: dens_mg
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: phi_init,phi_new
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: phi_reinit,phi
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: dphi_dx,dphi_dy
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: dphi_dz,s_phi0
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: h_phi,dens,mu,phim
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: abs_dphi_check
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: resmax
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: resfact
        DOUBLE PRECISION,allocatable,dimension(:,:,:) :: resfact1
      end type multidom

      type (multidom), pointer, dimension(:) :: dom

      end module multidata
!##########################################################################
