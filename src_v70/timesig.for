!##########################################################################
      subroutine timesig
!##########################################################################
      use multidata
      use vars
      use mpi
      implicit none
        
      INTEGER :: i,j,idfile,ib,jtime
      CHARACTER(8) :: numpt,x,y,z
      CHARACTER(25) :: unst

	do i=1,n_unstpt
	  do ib=1,nbp
	    if (dom_id(ib).eq.id_unst(i)) then
			idfile=499+i
		 write(idfile,500)ctime,
     &	dom(ib)%u(i_unst(i),j_unst(i),k_unst(i)),
     &	dom(ib)%v(i_unst(i),j_unst(i),k_unst(i)),
     &	dom(ib)%w(i_unst(i),j_unst(i),k_unst(i))
 !Add any other variable, indicate also in the heading in initial.for
!     &	dom(ib)%p(i_unst(i),j_unst(i),k_unst(i))
		endif
	   enddo
	enddo

  500  format(F16.7,6e18.8)

      end subroutine
!##########################################################################
      subroutine write_inflow
!##########################################################################
      use vars
      use mpi
      use multidata
      implicit none
	 
      INTEGER :: i,j,k,ib,strlen,imapp,nimapp,nx
      DOUBLE PRECISION :: xmapp,lxdom(idom+1)
      CHARACTER(5) :: name_end
      CHARACTER(4) :: domain
      CHARACTER(38) :: filename

	 xmapp=2.5d0		!location where inflow planes are taken
        do ib=1,nbp
         if (dom(ib)%xc(1).le.xmapp .and. 
     &        dom(ib)%xc(dom(ib)%ttc_i).gt.xmapp) then

	  do i=1,dom(ib)%ttc_i-1
           if (dom(ib)%xc(i).le.xmapp .and. 
     &        dom(ib)%xc(i+1).ge.xmapp) then
		 imapp=i

		lxdom=0 ; lxdom(1)=xst 
		do nx=2,idom+1
		  lxdom(nx)=(xcor(nx-2,2)-xcor(nx-2,1))+lxdom(nx-1)
		enddo
		Do nx=1,idom 
		 if(xmapp.gt.lxdom(nx) .and.xmapp.le.lxdom(nx+1))THEN
		  nimapp=nx-1
		 endif 	
		Enddo

        		write(name_end,'(I5)') ireadinlet
        		strlen=LEN(TRIM(ADJUSTL(name_end)))
       		name_end=REPEAT('0',(5-strlen))//
     &			TRIM(ADJUSTL(name_end)) 
			write(domain,'(I4)') (dom_id(ib)-nimapp)  !Identification of targetted domain
        		strlen=LEN(TRIM(ADJUSTL(domain)))
       		domain=REPEAT('0',(4-strlen))//
     &			TRIM(ADJUSTL(domain)) 

			filename='inflow/Inlet_'//domain//'_'//name_end//'.dat'

            	open (unit=ireadinlet+500, file=filename)
		  do k=dom(ib)%ksu-1,dom(ib)%keu+1
		  	do j=dom(ib)%jsu-1,dom(ib)%jeu+1
				write(ireadinlet+500,*)dom(ib)%u(imapp,j,k)
			enddo
			enddo		
		  do k=dom(ib)%ksu-1,dom(ib)%keu+1
		  	do j=dom(ib)%jsu-1,dom(ib)%jeu+1
				write(ireadinlet+500,*)dom(ib)%v(imapp,j,k)
			enddo
			enddo	
		  do k=dom(ib)%ksu-1,dom(ib)%keu+1
		  	do j=dom(ib)%jsu-1,dom(ib)%jeu+1
				write(ireadinlet+500,*)dom(ib)%w(imapp,j,k)
			enddo
			enddo		
			close(ireadinlet+500)

	     endif
	   enddo
          endif !xmapp
         enddo !ib
      end subroutine
!##########################################################################
      subroutine write_mapping
!##########################################################################
      use vars
      use mpi
      use multidata
      implicit none
        
      INTEGER :: i,j,k,ib,strlen,imapp,nimapp,nx
      DOUBLE PRECISION :: xmapp,lxdom(idom+1)
      CHARACTER(29) :: filename
      CHARACTER(4) :: domain

	  xmapp=2.5d0		!location where inflow planes are taken
        do ib=1,nbp
         if (dom(ib)%xc(1).le.xmapp .and. 
     &        dom(ib)%xc(dom(ib)%ttc_i).gt.xmapp) then

	  do i=1,dom(ib)%ttc_i-1
           if (dom(ib)%xc(i).le.xmapp .and. 
     &        dom(ib)%xc(i+1).ge.xmapp) then
		 imapp=i

		lxdom=0 ; lxdom(1)=xst 
		do nx=2,idom+1
		  lxdom(nx)=(xcor(nx-2,2)-xcor(nx-2,1))+lxdom(nx-1)
		enddo
		Do nx=1,idom 
		 if(xmapp.gt.lxdom(nx) .and.xmapp.le.lxdom(nx+1))THEN
		  nimapp=nx-1
		 endif 	
		Enddo

		write(domain,'(I4)') (dom_id(ib)-nimapp)
  		strlen=LEN(TRIM(ADJUSTL(domain)))
 		domain=REPEAT('0',(4-strlen))//
     &	TRIM(ADJUSTL(domain)) ! e.g. "0001"
	       filename='inflow/Mapping_'//domain//'.dat'
      	open (unit=ireadinlet+500, file=filename)
		  do k=dom(ib)%ksu-1,dom(ib)%keu+1
		  	do j=dom(ib)%jsu-1,dom(ib)%jeu+1
			write(ireadinlet+500,*)dom(ib)%u(imapp,j,k)
		  enddo ;enddo		
		  do k=dom(ib)%ksu-1,dom(ib)%keu+1
		  	do j=dom(ib)%jsu-1,dom(ib)%jeu+1
		  	 write(ireadinlet+500,*)dom(ib)%v(imapp,j,k)
		  enddo ;enddo	
		  do k=dom(ib)%ksu-1,dom(ib)%keu+1
		  	do j=dom(ib)%jsu-1,dom(ib)%jeu+1
		  	 write(ireadinlet+500,*)dom(ib)%w(imapp,j,k)
		  enddo ;enddo	
		close(ireadinlet+500)
	     endif
	   enddo
          endif !xmapp
         enddo !ib
      end subroutine
