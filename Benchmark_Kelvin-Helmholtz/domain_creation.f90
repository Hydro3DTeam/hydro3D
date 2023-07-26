	program domainsetup
	implicit none
	
!Domain resolution variables
	real :: xi,xe,dx
	real :: yi,ye,dy
	real :: zi,ze,dz
	double precision :: totcells
	integer :: i,j,k,xdiv,ydiv,zdiv,xcells,ycells,zcells
	integer :: np,npe,nd,ndpp,ndo,dom_number,rest_do,nl,l,d,n

!Interation variable and array declartion--------------------------------------
	integer :: mglevelxmin,mglevelymin,mglevelzmin,mglevelmin,divi,divj,divk
	integer :: jmax,status,i_div,j_div,k_div,element,version,ipts
	integer,allocatable,dimension(:,:) :: main
	integer,dimension(100) :: mglevelx,mglevely,mglevelz
        double precision,dimension(100) :: x,y,z
	real,dimension(100) :: idiv,jdiv,kdiv
	integer,allocatable,dimension(:,:,:) :: subdomatrix
	integer,pointer,dimension(:) :: division

!User input--------------------------------------------------------------------	
	x=0; y=0 ; z=0
	open(20,file='input_dom.cin')
	  do d=1,3
	    read(20,*)
	  end do
	  read(20,*) xdiv
	  read(20,*)
	  do i=1,xdiv
	    read(20,*) x(i)
	  end do
	  do d=1,3
	    read(20,*) 
	  end do 
          read(20,*) ydiv
	  read(20,*)
	  do j=1,ydiv
	    read(20,*) y(j)
	  end do
          do d=1,3
	    read(20,*) 
	  end do 
          read(20,*) zdiv
	  read(20,*)
	  do k=1,zdiv
	    read(20,*) z(k)
	  end do
	close(20)
	
	WRITE(6,'(a)') 'Enter (1) for v6 (2) for v7'
        read(*,*)  version
	WRITE(6,'(a)') 'Enter the coarser dx,dy,dz:'
	read(*,*) dx,dy,dz

	mglevelx(:)=0
	mglevelxmin=1000
	do i=1,xdiv-1
	  idiv(i)=(x(i+1)-x(i))/dx
	    mglevelx(i)=1
	    divi=nint(idiv(i))
	    do while (mod(divi,2).eq.0) 
	      divi=divi/2
	      mglevelx(i)=mglevelx(i)+1
	    end do
	  mglevelxmin=min(mglevelx(i),mglevelxmin)
	end do

	mglevely(:)=0
	mglevelymin=1000
	do j=1,ydiv-1
	  jdiv(j)=(y(j+1)-y(j))/dy
	    mglevely(j)=1
	    divj=nint(jdiv(j))
	    do while (mod(divj,2).eq.0) 
	      divj=divj/2
	      mglevely(j)=mglevely(j)+1
	    end do
	  mglevelymin=min(mglevely(j),mglevelymin)
	end do
	
	mglevelz(:)=0
	mglevelzmin=1000
	do k=1,zdiv-1
	  kdiv(k)=(z(k+1)-z(k))/dz
	    mglevelz(k)=1
	    divk=nint(kdiv(k))
	    do while (mod(divk,2).eq.0) 
	      divk=divk/2
	      mglevelz(k)=mglevelz(k)+1
	    end do
	  mglevelzmin=min(mglevelz(k),mglevelzmin)
	end do

	mglevelmin=min(mglevelxmin,mglevelymin,mglevelzmin)


	nd=(xdiv-1)*(ydiv-1)*(zdiv-1)
	i_div=xdiv-1
	j_div=ydiv-1
	k_div=zdiv-1

	xcells=int(x(xdiv)/dx)-1
	ycells=int(y(ydiv)/dy)-1
	zcells=int(z(zdiv)/dz)-1

	totcells=(xcells*ycells*zcells*1.d0)/10**(6)

	
	allocate(division(nd))
	ipts=0
	do n=1,nd
	  if(mod(nd,n).eq.0) then
	    ipts=ipts+1
	    division(ipts)=n
	  end if
	end do
	

	WRITE(6,'(a)') '================================================'
	WRITE(6,'(a)') 'Mglevel distribution (based on coarser mesh)' 
	WRITE(6,'(a)') '================================================'
	WRITE(6,'(a,I4)') 'x-direction: Mglevel x=',mglevelxmin
	WRITE(6,'(a)') '--------------------------'
	WRITE(6,*) '       x(i)      x(i+1)      idiv     Mglvl'
	do i=1,xdiv-1
	  WRITE(6,'(2F12.6,F10.2,I10)') x(i),x(i+1),idiv(i),mglevelx(i)
	end do
	WRITE(6,'(a)') '--------------------------'
	WRITE(6,'(a,I4)') 'y-direction: Mglevel y=',mglevelymin
	WRITE(6,'(a)') '--------------------------'
	WRITE(6,*) '       y(j)      y(j+1)      jdiv     Mglvl'
	do j=1,ydiv-1
	  WRITE(6,'(2F12.6,F10.2,I10)') y(j),y(j+1),jdiv(j),mglevely(j)
	end do
	WRITE(6,'(a)') '--------------------------'
	WRITE(6,'(a,I4)') 'z-direction: Mglevel z=',mglevelzmin
	WRITE(6,'(a)') '--------------------------'
	WRITE(6,*) '       z(k)      z(k+1)      kdiv     Mglvl'
	do k=1,zdiv-1
	  WRITE(6,'(2F12.6,F10.2,I10)') z(k),z(k+1),kdiv(k),mglevelz(k)
	end do
	WRITE(6,'(a)') '--------------------------'
	WRITE(6,'(a,I4)') 'Computational multigrid:',mglevelmin
	WRITE(6,'(a)') '================================================'
	WRITE(6,'(a,I10)') 'Total number of subdomain:',nd
	WRITE(6,'(a,F14.3)') 'Total number of cells:',totcells
	WRITE(6,'(a)') '================================================'
	WRITE(6,'(a)') 'Number of processor for balance load:'
	WRITE(6,'(a)') '-----------------------------------------'
	do n=1,ipts
	WRITE(6,'(I10)') division(n)
	end do
	WRITE(6,'(a)') '================================================'
	WRITE(6,'(a)') 'How many processor do you want to use ?'
	read *, np


	ndo=1
	dom_number=0
	
!INFODOM FILE
!------------------------------------------------------------------------	
	open (unit=30, file='infodom2.cin')
	
	write(*,*) 'Writing infodom.cin'
	
	if(version.eq.2) write(30,*) '1  (1)sub-domains description; (2)entire domain description (Lx,Ly,Lz)'
        write(30,*) nd, 'number of domains'
	write(30,'(a)') '========== Domain description =============='
	if(version.eq.2) write(30,'(a)') '0.0	      Lx: streamwise length'
	if(version.eq.2) write(30,'(a)') '0.0	      Ly: streamwise length'
	if(version.eq.2) write(30,'(a)') '0.0	      Lz: streamwise length'
	if(version.eq.2) write(30,'(I3,30a)') i_div,'           nx: subdomains in x-direction'
	if(version.eq.2) write(30,'(I3,30a)') j_div,'           ny: subdomains in x-direction'
	if(version.eq.2) write(30,'(I3,30a)') k_div,'           nz: subdomains in x-direction'
	if(version.eq.2) then
	  write(30,'(a)') '== CPU number, Local mesh refinement level,xini,xend,yini,yend,zini,zend=='
	end if
	
	do k=1,zdiv-1 ; do j=1,ydiv-1 ; do i=1,xdiv-1
	  write(30,20) dom_number, ndo, x(i),x(i+1), y(j), y(j+1), z(k), z(k+1) 
	  dom_number=dom_number+1
	enddo; enddo ;enddo
	write(30,*) '==============================================='
	write(30,*) i_div, 'number of divisions in i'
	write(30,*) j_div, 'number of divisions in j'
	write(30,*) k_div, 'number of divisions in k'
	
	do i=1,5
	  write(30,*) 
	end do 
        
      
        write(30,*) '################################'	    
	write(30,*) 'Subdomains disposition in Planes'
        write(30,*) '################################'	    

      element=0
      allocate(subdomatrix(0:i_div,0:j_div,0:k_div))
      do k=0,k_div-1 ; do j=0,j_div-1 ; do i=0,i_div-1
	subdomatrix(i,j,k)=element
	element=element+1
      end do ; end do ; end do
      
      write(30,*) '--------------------------------------------'
      write(30,*) 'XY-Planes'
      write(30,*) '--------------------------------------------'
      do k=0,k_div-1
	write(30,'(a,I3)') 'Level', k
        do j=j_div-1,0,-1
	  write(30,'(*(I5))') (subdomatrix(i,j,k), i=0,i_div-1)
        end do 
	write(30,*) 
      end do 
      
      write(30,*) '--------------------------------------------'
      write(30,*) 'XZ-Planes'
      write(30,*) '--------------------------------------------'
      do j=0,j_div-1
	write(30,'(a,I3)') 'Level', j
        do k=k_div-1,0,-1
	  write(30,'(*(I5))') (subdomatrix(i,j,k), i=0,i_div-1)
        end do 
	write(30,*) 
      end do 
 
      write(30,*) '--------------------------------------------'
      write(30,*) 'YZ-Planes'
      write(30,*) '--------------------------------------------'
      do i=0,i_div-1
	write(30,'(a,I3)') 'Level', j
        do k=k_div-1,0,-1
	  write(30,'(*(I5))') (subdomatrix(i,j,k), j=0,j_div-1)
        end do 
	write(30,*) 
      end do 

20	format(2I4, 6F15.5)
	close(unit=30)
!------------------------------------------------------------------------

!MDMAP FILE
!------------------------------------------------------------------------
	
	ndpp=nd/np
	npe=np-1
	rest_do=nd-ndpp*np
	write(*,*) ndpp,rest_do
	allocate(main(0:np,0:ndpp+2),stat=status)

!ASSIGN VALUES TO ARRAY
!----------------------
	do i=0,npe,1
	  main(i,0)=i
	  main(i,1)=ndpp
	  main(i,2)=i
!	  write(*,*) main(i,0), main(i,1), main(i,2)
	enddo
	
	do i=0,npe,1
	  do j=3,ndpp+1
	  main(i,j)= main(i,j-1)+np
!	  write(*,*) main(i,j-1), main(i,j)
	  enddo
	enddo
	
	if(rest_do.gt.0) then
	  do i=0,rest_do-1,1
	    main(i,1)=main(i,1)+1
	    main(i,ndpp+2)=main(i,ndpp+1)+np
	  end do
	end if

!WRITE ARRAY IN .CIN FILE
!------------------------
	open(unit=10, file='mdmap2.cin')
	write(*,*) 'Writing mdmap.cin'
	if(version.eq.2) write(10,*) '1  (1)sub-domains description; (2)entire domain description (Lx,Ly,Lz)'
        write(10,*) nd, 'number of domain'
        write(10,*) np, 'number of processor'
        write(10,*) '==============================================='
	do i=0,npe,1
	  if (i.le.rest_do-1) then
	     jmax=ndpp+2
	  else
	     jmax=ndpp+1
	  end if
        write(10,'(*(I5))') (main(i,j),j=0,jmax)    
	enddo
        write(10,*) '==============================================='
	close(unit=10)
!------------------------------------------------------------------------

	write(*,*) 'Both of your files have been created! Go check them!'
	end program
	
