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
        INTEGER,pointer,dimension(:) :: faz,cntp
        INTEGER,pointer,dimension(:,:,:) :: ntav1,ntav2
	 INTEGER,pointer,dimension(:,:,:) :: ibfactor
        INTEGER,pointer,dimension(:,:,:,:) :: ndimb,imbbodynum
	 DOUBLE PRECISION,pointer,dimension(:) :: tauw
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: S,So,Sm,Stm
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: SUtm,SVtm,SWtm
	 DOUBLE PRECISION,pointer,dimension(:,:,:) :: sfactor
        DOUBLE PRECISION,pointer,dimension(:) :: x,y,z,xc,yc,zc
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: u,v,w,p,pp
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: ksgs,ksgso
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: eps,epso
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: T,To,Tm,Ttm
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: su,ap,sup
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: ae,aw,as,an
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: at,ab
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: um,vm,wm
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: pm,ppm
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: vis,vism,epsm
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: uum,vvm,wwm
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: uvm,uwm,vwm
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: ustar,vstar
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: wstar
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: uo,uoo,vo
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: voo,wo,woo
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: stfcinf
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: facp1,facp2
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: facm1,facm2
        DOUBLE PRECISION,pointer,dimension(:) :: dh1,dh2,dh3
        DOUBLE PRECISION,pointer,dimension(:) :: cof
        DOUBLE PRECISION,pointer,dimension(:,:) :: tauwe,tauww
        DOUBLE PRECISION,pointer,dimension(:,:) :: tauws,tauwn
        DOUBLE PRECISION,pointer,dimension(:,:) :: tauwt,tauwb
        DOUBLE PRECISION,pointer,dimension(:,:) :: tauwe2,tauww2
        DOUBLE PRECISION,pointer,dimension(:,:) :: tauws2,tauwn2
        DOUBLE PRECISION,pointer,dimension(:,:) :: tauwt2,tauwb2
        DOUBLE PRECISION,pointer,dimension(:) :: sendb_m1,sendb_p1
        DOUBLE PRECISION,pointer,dimension(:) :: recvb_m1,recvb_p1
        DOUBLE PRECISION,pointer,dimension(:) :: sendb_m2,sendb_p2
        DOUBLE PRECISION,pointer,dimension(:) :: recvb_m2,recvb_p2
        DOUBLE PRECISION,pointer,dimension(:) :: sendb_m3,sendb_p3
        DOUBLE PRECISION,pointer,dimension(:) :: recvb_m3,recvb_p3
        DOUBLE PRECISION,pointer,dimension(:) :: sc1m,sc1p,rc1m,rc1p
        DOUBLE PRECISION,pointer,dimension(:) :: sc2m,sc2p,rc2m,rc2p
        DOUBLE PRECISION,pointer,dimension(:) :: sc3m,sc3p,rc3m,rc3p
        DOUBLE PRECISION,pointer,dimension(:) :: sc4m,sc4p,rc4m,rc4p
        DOUBLE PRECISION,pointer,dimension(:) :: se1m,se1p,re1m,re1p
        DOUBLE PRECISION,pointer,dimension(:) :: se2m,se2p,re2m,re2p
        DOUBLE PRECISION,pointer,dimension(:) :: se3m,se3p,re3m,re3p
        DOUBLE PRECISION,pointer,dimension(:) :: se4m,se4p,re4m,re4p
        DOUBLE PRECISION,pointer,dimension(:) :: se5m,se5p,re5m,re5p
        DOUBLE PRECISION,pointer,dimension(:) :: se6m,se6p,re6m,re6p
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: d1
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: dphi_dxplus
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: dphi_dyplus
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: dphi_dzplus
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: dphi_dxminus
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: dphi_dyminus
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: dphi_dzminus
!============================== LSM VARIABLES =============================
	 INTEGER :: tot
        INTEGER :: niul,njul,nkul,nivl,njvl,nkvl,niwl,njwl,nkwl
        INTEGER :: nipl,njpl,nkpl,nigl,njgl,nkgl
        INTEGER :: nipl2,njpl2,nkpl2
        INTEGER,pointer,dimension(:) :: ijkp_lsm
        REAL,pointer,dimension(:) :: sendb_m,sendb_p
        REAL,pointer,dimension(:) :: recvb_m,recvb_p
        DOUBLE PRECISION,pointer,dimension(:) :: dens_mg
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: phi_init,phi_new
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: phi_reinit,phi
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: dphi_dx,dphi_dy
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: dphi_dz,s_phi0
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: h_phi,dens,mu,phim
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: abs_dphi_check
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: resmax
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: resfact
        DOUBLE PRECISION,pointer,dimension(:,:,:) :: resfact1
      end type multidom

      type (multidom), pointer, dimension(:) :: dom

      end module multidata
!##########################################################################
