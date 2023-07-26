!#######################################################################
      module vars
!#######################################################################
      SAVE

      INTEGER :: alfabc,niter,nswp(4),nsweep,iter,nsweeps,ntime
      INTEGER :: ngg,nge,ngc,n_unstpt,OMP_threads,iaddinlet,ireadinlet
      INTEGER :: ngrid_input,ngrd_gl,maxcy,iproln,irestr
      INTEGER :: itmax,sweeps
      INTEGER :: numfile,numfile1,numfile2,numfile4,numfile5
      INTEGER :: conv_sch,diff_sch,sgs_model,solver,mg_itrsch
      INTEGER :: bc_w,bc_e,bc_s,bc_n,bc_b,bc_t
      INTEGER :: Tbc_w,Tbc_e,Tbc_s,Tbc_n,Tbc_b,Tbc_t
      INTEGER :: itime,itime_start,itime_end,n_out,ITMAX_PI
      INTEGER :: ipref,jpref,kpref,prefdom
      INTEGER :: pl,pl_ex,differencing,LMR,normal_inter,order
      DOUBLE PRECISION :: g_dx,g_dy,g_dz,dens,Re,eps,fac,Pr,beta,Th,Tc
      DOUBLE PRECISION :: qzero,qstpn,qstpn_y,qstpn_z
      DOUBLE PRECISION :: Sc_t,forcn,forcn_y,forcn_z
      DOUBLE PRECISION :: flomas,rmax,alfapr,resor,fric
      DOUBLE PRECISION :: ctime,dt,dtavg,dtsum,safety_factor,noise,Mdef
      DOUBLE PRECISION :: t_start_averaging1,t_start_averaging2,ubulk
      DOUBLE PRECISION :: xst,xen,yst,yen,zst,zen,wtime_cfl,dtinit
      CHARACTER(80) :: keyword,L_n
      LOGICAL :: LRESTART,LIMB,SGS,PERIODIC,LENERGY,LROUGH
      LOGICAL :: pressureforce,pressureforce_y,pressureforce_z
      LOGICAL :: time_averaging,reinitmean
      LOGICAL :: variTS,save_inflow,read_inflow,L_dt
      LOGICAL :: LTECP,LTECBIN,LTURB,LINST,LPLAN,LPAR,LLSM,LLPT,LSCALAR

!---- SEM VARIABLES ----------------------------------------------------
      INTEGER:: UPROF
      DOUBLE PRECISION :: TI_SEM

!---- LSM VARIABLES ----------------------------------------------------
      INTEGER :: ntime_reinit,ngrid,numfile3
      DOUBLE PRECISION :: reldif_LSM,length,densl,densg,nul,nug,mul,mug
      DOUBLE PRECISION :: cfl_lsm,uprev,grx,gry,grz,slope
      DOUBLE PRECISION :: densip12,densjp12,denskp12
      DOUBLE PRECISION :: densim12,densjm12,denskm12
      DOUBLE PRECISION :: muip12,mujp12,mukp12,muim12,mujm12,mukm12
      LOGICAL :: L_LSM,REINIT,L_LSMbase,L_LSMinit,LENDS
      LOGICAL :: L_anim_phi,L_anim_grd

      end module
!#######################################################################
