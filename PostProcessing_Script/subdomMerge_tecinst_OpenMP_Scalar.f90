!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     PURPOSE: Import all the tecbin files, merge their respective subdomains
!	       into one. The fulldomain can be written out as ASCII (.dat)
!	       or tecplot binary (.plt). Particularly usefull for tecplot
!	       macros and to monitor simulation from the supercomputer.

!     HYDRO3D: Works with V6.0 and V7.0

!     LIMITS: 
!     This software is not compatible with LMR simulations.
!     Check in (post.for) the tecbin variables and make sure that its export 
!     structure match the way it is read by this software. Anchor: !% TECBIN INPUT.

!     TUTO: WILL BE COMING SOON

!     AUTHORS: Arthur Hajaali, hajaalia@cardiff.ac.uk
!     Do not hesitate to contact the author if you have any question.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      program subdomMerge
      implicit none

      integer :: ID,num_threads,fileunit,n,Smax
      integer :: ierr,sn,ipts,l,dif,ie,itime,etime,dtime,itmax,it,t
      integer :: i,j,k,ii,jj,kk,iii,jjj,kkk,div,idiv,jdiv,kdiv,ni,nj,nk
      integer :: dom,ndom,npl,di,dj,dk,dT,fi,fj,fk,maxdom
      double precision :: ctime,dx,dy,dz,start,finish,delta
      double precision :: ccx,ccy,ccz,xp,yp,zp,xmag,ymag,zmag,thetaSw
      double precision,parameter :: pi=4.d0*atan(1.d0)
      integer,pointer,dimension(:) :: nd,dxe,dye,dze,rdiv
      integer,pointer,dimension(:) :: idom,jdom,kdom,idex,jdex,kdex
      integer,pointer,dimension(:) :: fi_inc,fj_inc,fk_inc
      integer,pointer,dimension(:) :: domtti,domttj,domttk
      character(200):: fn
      character(80) :: command,chID
      character(100) :: gf,chb,chb2,chdi,chdj,chdk,chctime,path,chb3
      character(200) :: fvariables,fs,chS
      logical :: binary

!%  FULLDOMAIN VARIABLES:
      double precision,pointer,dimension(:,:,:) :: xd,yd,zd
      double precision,pointer,dimension(:,:,:) :: ud,vd,wd,pd
      double precision,pointer,dimension(:,:,:,:) :: Sd,pthetad
      double precision,pointer,dimension(:,:,:) :: resd,celltyped
      double precision,pointer,dimension(:,:,:) :: rad

!% TECBIN VARIABLES:
      double precision,pointer,dimension(:,:,:) :: x,y,z
      double precision,pointer,dimension(:,:,:) :: u,v,w,p
      double precision,pointer,dimension(:,:,:) :: res,celltype
      double precision,pointer,dimension(:,:,:,:) :: S

!% TEMPORARILY VARIABLES TO INTERPOLATE AND AVERAGE THE FLOW FIELD:
      double precision,pointer,dimension(:,:,:) :: u_cn,v_cn,w_cn
      double precision,pointer,dimension(:,:,:) :: um_cn,vm_cn,wm_cn
      double precision,pointer,dimension(:,:,:) :: p_cn
      double precision,pointer,dimension(:,:,:,:) :: S_cn


      double precision,pointer,dimension(:,:) :: xcor,ycor,zcor
      integer,pointer,dimension(:) :: tti,ttj,ttk
      integer :: tti2,ttj2,ttk2,version



! *****  TO compile***
! ifort -qopenmp -parallel -fpp subdomMerge_tecinst_OpenMP_Scalar.f90 -o subdomMerge_tecinst_OpenMP_Scalar.exe


!% USER INPUT:
      version=7

      open(20,file='input_file.txt')
        read(20,*) delta,Smax
        read(20,*) itime,etime,dtime
        read(20,*) num_threads
      close(20)

      path='preplot'
      dx=delta ; dy=delta ; dz=delta
      itmax=nint((etime-itime)/dtime*1.d0)+1

      write(6,'(a,3E14.4)') 'dx,dy,dz:',dx,dy,dz
      write(6,'(a,4I9)') 'itime,etime,dtime,itmax:',itime,etime,dtime,itmax
      write(6,'(a,I9)') 'num_threads:',num_threads

      
      !$ call OMP_SET_NUM_THREADS(num_threads)

      call cpu_time(start)

!% OPEN THE INFODOM TO GATHER THE INFORMATION ABOUT THE DOMAIN:
      open(unit=22, file='infodom.cin',iostat=ierr) 
	if(ierr.ne.0) then
	  WRITE(6,*) 'The infodom.cin file is missing'
	  STOP
	end if
	
	if(version.eq.7) read(22,*) 
	read(22,*) ndom
	read(22,*) 

	if(version.eq.7) read(22,*)
	if(version.eq.7) read(22,*)
	if(version.eq.7) read(22,*)
	if(version.eq.7) read(22,*) idiv
	if(version.eq.7) read(22,*) jdiv
	if(version.eq.7) read(22,*) kdiv
	if(version.eq.7) read(22,*)
	
	  
  
	allocate(tti(ndom),ttj(ndom),ttk(ndom))
	allocate(xcor(2,ndom),ycor(2,ndom),zcor(2,ndom))
       allocate(nd(ndom),rdiv(ndom),dxe(ndom),dye(ndom),dze(ndom))
	allocate(idom(ndom),jdom(ndom),kdom(ndom))
	allocate(idex(ndom),jdex(ndom),kdex(ndom))
	allocate(domtti(ndom),domttj(ndom),domttk(ndom))
	allocate(fi_inc(ndom),fj_inc(ndom),fk_inc(ndom))


	do l=0,ndom-1
	  read(22,*) dom,rdiv(l),xcor(1,l),xcor(2,l),ycor(1,l),ycor(2,l),zcor(1,l),zcor(2,l)
	end do
	read(22,*)
	if(version.eq.6) read(22,*) idiv
	if(version.eq.6) read(22,*) jdiv
	if(version.eq.6) read(22,*) kdiv
      close(22)
      
      WRITE(6,*) 'INFODOM INFORMATION:'
      WRITE(6,*) 'ndom,idiv,jdiv,kdiv:',ndom,idiv,jdiv,kdiv
      WRITE(6,*) '                    '

      do l=0,ndom-1
	tti(l)=nint((xcor(2,l)-xcor(1,l))/dx)+2
	ttj(l)=nint((ycor(2,l)-ycor(1,l))/dy)+2
	ttk(l)=nint((zcor(2,l)-zcor(1,l))/dz)+2
      end do
      

      ni=0 ; nj=0 ; nk=0 ; fk=0 ; fj=0 ; fi=0  ; ipts=0 ; l=0
      do k=0,kdiv-1 
	 do j=0,jdiv-1 
	    do i=0,idiv-1
	      l=l+1
		ipts=ipts+1
		nd(ipts)=l-1				!% Store domain number
		fi_inc(ipts)=fi				!% Store istart div for each domain
		fj_inc(ipts)=fj				!% Store jstart div for each domain
		fk_inc(ipts)=fk				!% Store kstart div for each domain
		domtti(ipts)=1				!% Store domain i division minus the edge
		domttj(ipts)=1				!% Store domain k division minus the edge
		domttk(ipts)=1				!% Store domain k division minus the edge
		if(i.eq.idiv-1) domtti(ipts)=0		!% Add edge if domain shares with fulldom i edge
		if(j.eq.jdiv-1) domttj(ipts)=0		!% Add edge id domain shares with fulldom k edge
		if(k.eq.kdiv-1) domttk(ipts)=0		!% Add edge id domain shares with fulldom k edge
!	  WRITE(6,'(7I8)') nd(ipts),fi_inc(ipts),fj_inc(ipts)  &
!     &,fk_inc(ipts),domtti(ipts),domttj(ipts),domttk(ipts)
	      fi=fi+tti(l-1)-2                    !% Find the domain istart
	    end do
	    fj=fj+ttj(l-idiv)-2 ; fi=0		  !% Find the domain jstart
	 end do
	 fk=fk+ttk(l-idiv*jdiv)-2 ; fj=0	  !% Find the domain kstart
      end do 
      maxdom=ipts
      di=fi_inc(maxdom)+tti(nd(maxdom))-1
      dj=fj_inc(maxdom)+ttj(nd(maxdom))-1
      dk=fk_inc(maxdom)+ttk(nd(maxdom))-1
      dT=di*dj*dk
      WRITE(6,*)  '               '
      WRITE(6,*) 'TOTAL DIVISIONS:'
      WRITE(6,*) 'di,dj,dk,dT:',di,dj,dk,dT
      WRITE(6,*) 'ndom,maxdom:',ndom,maxdom
      WRITE(6,*) 


      WRITE(6,*) 'DOMAIN PARALLEL MERGING BEGIN!:'      

!$OMP PARALLEL DEFAULT(PRIVATE),shared(di,dj,dk,nd,itmax,itime,dtime,ndom,rdiv,fk_inc,fi_inc,fj_inc,domtti,domttj,domttk,Smax,ccx,ccy,ccz)
!$OMP DO
      DO it=1,itmax
	t=itime+(it-1)*dtime
!	WRITE(6,*) t
!% IMPORT TECBIN VARIABLES -> INTERPOLATE THEM -> STORE THEM IN GLOBAL DOMAIN
!============================================================================     
      allocate(xd(di,dj,dk),yd(di,dj,dk),zd(di,dj,dk))
      allocate(ud(di,dj,dk),vd(di,dj,dk),wd(di,dj,dk))
      allocate(pd(di,dj,dk),Sd(di,dj,dk,Smax),pthetad(di,dj,dk,2))
      allocate(rad(di,dj,dk))
      allocate(resd(di,dj,dk),celltyped(di,dj,dk))

      
      write(chb2,'(i6)') t
      fn='fullinst_'//trim(adjustl(chb2))//'.dat'
      WRITE(6,'(a)') 'Importing data for '//trim(adjustl(fn))//' ...' 

!      WRITE(6,*) 'IMPORTATION OF THE DATA INTO THE FULLDOMAIN:'
      fi=0 ; fj=0 ; fk=0
      do ipts=1,ndom
	
	  write(chb,'(i4)') nd(ipts)
	  sn=len(trim(adjustl(chb)))
	  chb=repeat('0',(4-sn))//trim(adjustl(chb))
	  
	  write(chb2,'(i6)') t
!	  WRITE(6,*) t,nd(ipts)
	  sn=len(trim(adjustl(chb2)))
	  chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
	  gf='tecinst'//trim(adjustl(chb2))//'_'//trim(adjustl(chb))//'.dat'
!	  WRITE(6,*) trim(adjustl(gf))
	  
	  fileunit=it*ndom
	  open (unit=fileunit, file=gf)
	    read(fileunit,*) tti2,ttj2,ttk2
	    read(fileunit,*) npl
          read(fileunit,*)     !inind,jnind,knind CHANGE/remove based on the version 
	
      allocate(x(tti2,ttj2,ttk2),y(tti2,ttj2,ttk2),z(tti2,ttj2,ttk2))
      allocate(u(tti2,ttj2,ttk2),v(tti2,ttj2,ttk2),w(tti2,ttj2,ttk2))
      allocate(p(tti2,ttj2,ttk2),S(tti2,ttj2,ttk2,Smax))
      allocate(res(tti2,ttj2,ttk2),celltype(tti2,ttj2,ttk2))
      
      allocate(u_cn(tti2,ttj2,ttk2),v_cn(tti2,ttj2,ttk2))
      allocate(w_cn(tti2,ttj2,ttk2))
      allocate(p_cn(tti2,ttj2,ttk2),S_cn(tti2,ttj2,ttk2,Smax))

	  ctime=0.d0
	  x(:,:,:)=0.d0	; y(:,:,:)=0.d0	; z(:,:,:)=0.d0
	  u(:,:,:)=0.d0	; v(:,:,:)=0.d0	; w(:,:,:)=0.d0
	  p(:,:,:)=0.d0	; S(:,:,:,:)=0.d0	

	  u_cn(:,:,:)=0.d0   ; v_cn(:,:,:)=0.d0	; w_cn(:,:,:)=0.d0
	  p_cn(:,:,:)=0.d0   ; S_cn(:,:,:,:)=0.d0	
	  res(:,:,:)=0.d0    ; celltype(:,:,:)=0.d0

	  do k=1,ttk2 ; do j=1,ttj2 ; do i=1,tti2
	  read(fileunit,*)  ctime,x(i,j,k),y(i,j,k),z(i,j,k) & !% CHANGE TECBIN INPUTa
     &	  ,u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k),S(i,j,k,1:Smax)
	  end do ; end do ; end do

	  close(fileunit)

!% INTERPOLATE AND AVERAGE THE NODE VELOCITY TO GENERATE CELLS VALUE:
!%====================================================================

	do k=npl,ttk2-npl
	 do j=npl,ttj2-npl
	  do i=npl,tti2-npl

		 u_cn(i,j,k)  =0.25*(u(i,j,k)+ &
     &u(i+1,j,k)+u(i,j,k+1)+		       &
     &u(i+1,j,k+1))
		 v_cn(i,j,k)  =0.25*(v(i,j,k)+ &
     &v(i+1,j,k)+v(i,j,k+1)+                   &
     &v(i+1,j,k+1))
		 w_cn(i,j,k)  =0.25*(w(i,j,k)+ &
     &w(i+1,j,k)+w(i,j,k+1)+                   &
     &w(i+1,j,k+1))
                 p_cn(i,j,k)  =0.125*(p(i,j,k)+ &
     &p(i+1,j,k)    +p(i,j+1,k)+                &
     &p(i+1,j+1,k)  +p(i,j,k+1)+                &
     &p(i+1,j,k+1)  +p(i,j+1,k+1)+              &
     &p(i+1,j+1,k+1))
		
	    do n=1,Smax
              S_cn(i,j,k,n)  =0.125*(S(i,j,k,n)+    &
     &S(i+1,j,k,n)    +S(i,j+1,k,n)+                &
     &S(i+1,j+1,k,n)  +S(i,j,k+1,n)+                &
     &S(i+1,j,k+1,n)  +S(i,j+1,k+1,n)+              &
     &S(i+1,j+1,k+1,n))
	    end do

          end do
	 end do
	end do

      u(:,:,:)=u_cn(:,:,:)
      v(:,:,:)=v_cn(:,:,:)
      w(:,:,:)=w_cn(:,:,:)
      p(:,:,:)=p_cn(:,:,:)
      if(Smax.ge.1) S(:,:,:,:)=S_cn(:,:,:,:)
!% INTEGRATE THE TECBIN SUBDOMAIN INTO THE GLOBAL DOMAIN:
!%=======================================================
      fk=fk_inc(ipts)
      do k=npl,ttk2-npl-domttk(ipts),rdiv(ipts-1)
	fk=fk+1 ; fj=fj_inc(ipts)
	do j=npl,ttj2-npl-domttj(ipts),rdiv(ipts-1)
	  fj=fj+1 ; fi=fi_inc(ipts)
	  do i=npl,tti2-npl-domtti(ipts),rdiv(ipts-1)
	    fi=fi+1
	    xd(fi,fj,fk)=x(i,j,k)
	    yd(fi,fj,fk)=y(i,j,k)
	    zd(fi,fj,fk)=z(i,j,k)
	    ud(fi,fj,fk)=u(i,j,k)
	    vd(fi,fj,fk)=v(i,j,k)
	    wd(fi,fj,fk)=w(i,j,k)
	    pd(fi,fj,fk)=p(i,j,k)
	    do n=1,Smax
	    Sd(fi,fj,fk,n)=S_cn(i,j,k,n)
	    end do
	  end do
	end do
      end do

33    format(a,3I4,3F10.2)  
    
      deallocate(x,y,z)
      deallocate(u,v,w,p,S)
      deallocate(celltype,res)
      
      deallocate(u_cn,v_cn,w_cn)
      deallocate(p_cn,S_cn)

      end do  !ipts
    

!% Write out the full domain, try to write as binary szplt (Check on tecplot)
!%========================================================================================
!      WRITE(6,*) '                         '
!      WRITE(6,*) 'FULLDOMAIN DATA FILE WRITE:'
      fvariables='Variables=x,y,z,u,v,w,p'!,celltype,res'
      do n=1,Smax
	write(chS,'(I4)') n 
	fs=',S'//trim(adjustl(chS))
	fvariables=trim(adjustl(fvariables))//trim(adjustl(fs))
      end do

      write(chb3,'(I10)') it*dtime
      write(chdi,'(I10)') di
      write(chdj,'(I10)') dj
      write(chdk,'(I10)') dk
      write(chctime,'(F10.5)') ctime
      fn='fullinst_'//trim(adjustl(chb2))//'.dat'
      fileunit=44+it*1000
      
      open(unit=fileunit,file=fn)
      write(fileunit,'(A)') trim(adjustl(fvariables))
      write(fileunit,'(A)') 'Zone T="it='//trim(adjustl(chb3))//'", i='           &
     & //trim(adjustl(chdi))//', j='                                              &
     & //trim(adjustl(chdj))//' ,k='//trim(adjustl(chdk))//', STRANDID='          &
     & //trim(adjustl(chb3))//', SOLUTIONTIME= '                                  &
     & //trim(adjustl(chctime))//', f=point'

      do fk=1,dk ; do fj=1,dj ; do fi=1,di
	write(fileunit,44) xd(fi,fj,fk),yd(fi,fj,fk),zd(fi,fj,fk)              &
     &		           ,ud(fi,fj,fk),vd(fi,fj,fk),wd(fi,fj,fk)              &
     &                   ,pd(fi,fj,fk),Sd(fi,fj,fk,1:Smax)
      end do ; end do  
      end do

44    format(24E25.8)
      close(fileunit)

!========================================================================================
      deallocate(xd,yd,zd)
      deallocate(ud,vd,wd)
      deallocate(pd,Sd,pthetad)
      deallocate(rad)
      deallocate(resd,celltyped)
      
!      !$ ID = omp_get_thread_num()
!      WRITE(chID,'(I4)') ID 
      WRITE(6,'(a)') trim(adjustl(fn))//' Completed'!, THREAD:'//trim(adjustl(chID))

      END DO
!$OMP END DO
!$OMP END PARALLEL
      
      
      !WRITE(6,*) 'Do you want it in binary format fulltecturb.plt? T/F'
      !read(*,*) binary
      !if(binary.eq..TRUE.) then
!	do it=1,itmax
!	command=trim(adjustl(path))//' '//trim(adjustl(fn(it)))
!	CALL SYSTEM (command)
!	end do
      !end if
      
!      WRITE(6,*) 'PROGRAM FINISHED!'
      call cpu_time(finish)
      WRITE(6,'(a,F15.3)') 'PROGRAM EXECUTION TIME:',(finish-start)
      
      end program
