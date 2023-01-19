# read\_mdmap

**File:** alloc\_dom.for\
\
**Subroutine Called in:** \
fdstag (fdstag.for)\
\
**Purpose:**\
****This subroutine reads the mdmap.cin input file to store the parallisation set up of Hydro3D.\
\
**User Likeness to alter the file:** \
****Rarely (Only if you change the layout of mdmap.cin)

```
##########################################################################
        subroutine read_mdmap
!##########################################################################
        use multidata
        use vars
        use mpi
        implicit none
        integer i,j,k,ib,nptemp,myranktemp,nbtemp,cpu_no,type_mdmp
        integer :: subdom_cpu,remain_dom
        integer,allocatable,dimension(:) :: domtemp,buf_domindid
        character*80 :: dummyline

	if (myrank.eq.0) then
        numfile=1001
        open (unit=numfile, file='output.dat')
        end if

        open (unit=12, file='mdmap.cin')
        read (12,*) type_mdmp
        read (12,*) num_domains !number of domains
        read (12,*) nptemp !number of processors

        if (nprocs .ne. nptemp) then
           print*, '=====ERROR====='
           print*, 'number of cpus do not match map file'
           stop
        end if

        if (num_domains .gt. 9999) then
           print*, '=====ERROR====='
           print*, 'number of domains are exceeding',
     & ' the limit in exchange subroutine'
           stop
        end if

        allocate (dom_id(num_domains),domtemp(num_domains))
        allocate (dom_indid(0:num_domains-1),dom_ad(0:num_domains-1))
        allocate (buf_domindid(0:num_domains-1))
       	allocate (imbinblk(num_domains))	!Pablo

        dom_ad=-1
        read (12,*) dummyline

        IF(type_mdmp.eq.1) then
        nbpmax=0
        do i=1,nprocs
           read(12,*) myranktemp,nbtemp,domtemp(1:nbtemp)
           nbpmax=max(nbpmax,nbtemp)

           do ib=1,nbtemp
              dom_ad(domtemp(ib))=myranktemp
           end do

           if(myranktemp.eq.myrank) then
              nbp=nbtemp                                    !number of domains for this processor
              dom_id(1:nbp)=domtemp(1:nbp)
              do ib=0,num_domains-1
                 dom_indid(ib)=-1
              end do
              do ib=1,nbp
                 dom_indid(dom_id(ib))=ib
              end do
           end if
        end do

        close(12)

        ELSE
        
        close(12)

          dom_indid=-1
          subdom_cpu = int(num_domains/nptemp)
          remain_dom = num_domains-subdom_cpu*nptemp

	 if(mod(num_domains,nptemp).eq.0)remain_dom=subdom_cpu
          IF(myrank.eq.0) then
           write(6,*)'Mdmap generation "2"'
           write(6,*)'Sub-domains per CPU    :',subdom_cpu
           write(6,*)'Sub-domains to last CPU:',remain_dom
          ENDif

         do j=1,nptemp-1     !myranktemp
          do i=1,subdom_cpu          !nptemp
            domtemp(i)=(i-1)+(j-1)*subdom_cpu
            dom_ad(domtemp(i))=j-1
         enddo
          if(myrank.eq.j) then
           do i=1,subdom_cpu
             dom_id(i)=(i-1)+(j-1)*subdom_cpu
             dom_indid(dom_id(i))=i
          enddo
         endif
        enddo

          do i=1,remain_dom          !nptemp
            domtemp(i)=(i-1)+(nptemp-1)*subdom_cpu
            dom_ad(domtemp(i))=nptemp-1
         enddo
          if(myrank.eq.nptemp-1) then
           do i=1,remain_dom
             dom_id(i)=(i-1)+(nptemp-1)*subdom_cpu
             dom_indid(dom_id(i))=i
           enddo
          endif

        ENDIF

        buf_domindid = dom_indid

        call MPI_ALLREDUCE(buf_domindid,dom_indid,num_domains,
     &MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)

        if(myrank.eq.0) then
           do i=0,num_domains-1
              if(dom_ad(i).eq.-1) then
                 print*,'unidentified domain in mdmap no:',i
                 print*,'ERROR! check mdmap.cin'
                 stop
              end if
              print*, 'dom:',i,' cpu:',dom_ad(i),' ib:',dom_indid(i) 
        write(numfile,*) 'dom:',i,' cpu:',dom_ad(i),' ib:',dom_indid(i)
           end do
        end if

        allocate (dom(nbp))

        call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        return
        end subroutine read_mdmap
```
