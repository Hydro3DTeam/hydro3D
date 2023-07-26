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
      program subdomMerge_preplot
      use omp_lib
      implicit none
 
      integer :: ndom,itime,etime,dtime,itmax,it,jtime,sn,ID,num_threads
      double precision :: finish,start,dum
      character(80) :: command
      character(100) :: chb2,path,chID

! ifort -qopenmp -parallel -fpp subdomMerge_preplot_OpenMP.f90 -o subdomMerge_preplot_OpenMP.exe
	

      
!% USER INPUT:
      open(20,file='input_file.txt')
        read(20,*) dum,dum
        read(20,*) itime,etime,dtime
        read(20,*) num_threads
      close(20)
      itmax=nint((etime-itime)/dtime*1.d0)+1
      
      !$ call OMP_SET_NUM_THREADS(num_threads)
      
      call cpu_time(start)

!$OMP PARALLEL DEFAULT(PRIVATE),shared(itime,etime,dtime,itmax)
!$OMP DO
      do it=1,itmax
       ID = omp_get_thread_num()
       WRITE(chID,'(I4)') ID 
       jtime = itime+((it-1)*dtime)
!       WRITE(6,*) jtime,'/',etime,'     |    ',it,'/',itmax
       write(chb2,'(i6)') jtime
       sn=len(trim(adjustl(chb2)))
       chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
	
       WRITE(6,*) 'fullinst_'//trim(adjustl(chb2))//'.dat ASCII2 -> Binary'

       path='preplot'//trim(adjustl(chID))    
       command=trim(adjustl(path))//' '//trim(adjustl('fullinst_'//trim(adjustl(chb2))//'.dat > /dev/null'))
       CALL SYSTEM (command)
       
       WRITE(6,*) 'fullinst_'//trim(adjustl(chb2))//'.plt created'
      end do
!$OMP END DO
!$OMP END PARALLEL
      
      call cpu_time(finish)
      WRITE(6,'(a,F15.3)') 'PROGRAM EXECUTION TIME:',(finish-start)/num_threads
      end program
