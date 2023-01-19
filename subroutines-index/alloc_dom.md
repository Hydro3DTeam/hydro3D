# alloc\_dom

**File:** alloc\_dom.for\
\
**Subroutine Called in:** \
fdstag (fdstag.for)\
\
**Purpose:**\
****This subroutine defines the spatial location and neighbour of each subdomain, it also checks the integrity of the Local Mesh Refinement (LMR) to ensure that the eight neighbours around the subdomain possess an LMR difference of 2 maximum.\
\
**User Likeness to alter the file:** \
****Never

```
##########################################################################
        subroutine alloc_dom
!##########################################################################
        use multidata
        use vars
        use mpi
        implicit none
        integer i,j,k,ib,dir
        integer sync_dir,cpu_next,cpu_prev,rdivmy,rdivng,my_cor
        integer ng_p,ng_n

        if(rdivmax.gt.1) then

!===============================================================
        do sync_dir = 1,3
        do ib=1,nbp

           if (sync_dir.eq.1)  then
              cpu_prev=dom(ib)%iprev
              cpu_next=dom(ib)%inext
           else if (sync_dir.eq.2)  then
              cpu_prev=dom(ib)%jprev
              cpu_next=dom(ib)%jnext
           else
              cpu_prev=dom(ib)%kprev
              cpu_next=dom(ib)%knext
           end if

           if (cpu_prev.ge.0) then
              rdivmy=rdiv(dom_id(ib))
              rdivng=rdiv(cpu_prev)

              if (rdivmy.ne.rdivng) then
                 if (rdivmy.gt.rdivng) then
                    if (rdivmy .ne. 2*rdivng) then
                       print*, '====ERROR1===>  wrong rdiv',dom_id(ib)
                       stop
                    end if
                 else if (rdivmy.lt.rdivng) then
                    if (2*rdivmy .ne. rdivng) then
                       print*, '====ERROR2===>  wrong rdiv',dom_id(ib)
                       stop
                    end if
                 else
                    print*, '===ERROR==='
                 end if
              end if
           end if

           if (cpu_next.ge.0) then
              rdivmy=rdiv(dom_id(ib))
              rdivng=rdiv(cpu_next)

              if (rdivmy.ne.rdivng) then
                 if (rdivmy.gt.rdivng) then
                    if (rdivmy .ne. 2*rdivng) then
                       print*, '====ERROR3===>  wrong rdiv',dom_id(ib)
                       stop
                    end if
                 else if (rdivmy.lt.rdivng) then
                    if (2*rdivmy .ne. rdivng) then
                       print*, '====ERROR4===>  wrong rdiv',dom_id(ib)
                       stop
                    end if
                 else
                    print*, '===ERROR==='
                 end if
              end if
           end if


           if (sync_dir.eq.1)  then
           
!======================================================================
!******* Corner-Neighbors
!======================================================================
           do dir=1,4

           select case (dir)
              Case (1)  
                 ng_p=dom(ib)%corprev1
                 ng_n=dom(ib)%cornext1
              Case (2)
                 ng_p=dom(ib)%corprev2
                 ng_n=dom(ib)%cornext2
              Case (3)  
                 ng_p=dom(ib)%corprev3
                 ng_n=dom(ib)%cornext3
              Case (4)
                 ng_p=dom(ib)%corprev4
                 ng_n=dom(ib)%cornext4
           end select

           if (ng_p.ge.0) then
              rdivmy=rdiv(dom_id(ib))
              rdivng=rdiv(ng_p)

              if (rdivmy.ne.rdivng) then
                 if (rdivmy.gt.rdivng) then
                    if (rdivmy .ne. 2*rdivng) then
        print*, '====ERRORcor-p===>  wrong rdiv',dom_id(ib),dir
                       stop
                    end if
                 else !if (rdivmy.lt.rdivng) then
                    if (2*rdivmy .ne. rdivng) then
        print*, '====ERRORcor-p===>  wrong rdiv',dom_id(ib),dir
                       stop
                    end if
                 end if
              end if
           end if
           if (ng_n.ge.0) then
              rdivmy=rdiv(dom_id(ib))
              rdivng=rdiv(ng_n)

              if (rdivmy.ne.rdivng) then
                 if (rdivmy.gt.rdivng) then
                    if (rdivmy .ne. 2*rdivng) then
        print*, '====ERRORcor-n===>  wrong rdiv',dom_id(ib),dir
                       stop
                    end if
                 else !if (rdivmy.lt.rdivng) then
                    if (2*rdivmy .ne. rdivng) then
        print*, '====ERRORcor-n===>  wrong rdiv',dom_id(ib),dir
                       stop
                    end if
                 end if
              end if
           end if

           end do
!======================================================================
!******* Edge-Neighbors
!======================================================================
           do dir=1,6

           select case (dir)
              Case (1)  
                 ng_p=dom(ib)%edgprev1
                 ng_n=dom(ib)%edgnext1
              Case (2)
                 ng_p=dom(ib)%edgprev2
                 ng_n=dom(ib)%edgnext2
              Case (3)  
                 ng_p=dom(ib)%edgprev3
                 ng_n=dom(ib)%edgnext3
              Case (4)
                 ng_p=dom(ib)%edgprev4
                 ng_n=dom(ib)%edgnext4
              Case (5)  
                 ng_p=dom(ib)%edgprev5
                 ng_n=dom(ib)%edgnext5
              Case (6)
                 ng_p=dom(ib)%edgprev6
                 ng_n=dom(ib)%edgnext6
           end select

           if (ng_p.ge.0) then
              rdivmy=rdiv(dom_id(ib))
              rdivng=rdiv(ng_p)

              if (rdivmy.ne.rdivng) then
                 if (rdivmy.gt.rdivng) then
                    if (rdivmy .ne. 2*rdivng) then
        print*, '====ERRORedg-p===>  wrong rdiv',dom_id(ib),dir
                       stop
                    end if
                 else !if (rdivmy.lt.rdivng) then
                    if (2*rdivmy .ne. rdivng) then
        print*, '====ERRORedg-p===>  wrong rdiv',dom_id(ib),dir
                       stop
                    end if
                 end if
              end if
           end if
           if (ng_n.ge.0) then
              rdivmy=rdiv(dom_id(ib))
              rdivng=rdiv(ng_n)

              if (rdivmy.ne.rdivng) then
                 if (rdivmy.gt.rdivng) then
                    if (rdivmy .ne. 2*rdivng) then
        print*, '====ERRORedg-n===>  wrong rdiv',dom_id(ib),dir
                       stop
                    end if
                 else !if (rdivmy.lt.rdivng) then
                    if (2*rdivmy .ne. rdivng) then
        print*, '====ERRORedg-n===>  wrong rdiv',dom_id(ib),dir
                       stop
                    end if
                 end if
              end if
           end if

           end do
!======================================================================
!======================================================================
           end if !if (sync_dir.eq.1)  then
        end do
        end do
!===============================================================

        do ib=1,nbp

        if(dom(ib)%bc_west.eq.5 .or. dom(ib)%bc_east.eq.5) then
!=== West ===> 
           if (dom(ib)%iprev.lt.0) then
              if (dom(ib)%inext.ge.0) then
                 my_cor=dom(ib)%per_ip
                 if(rdiv(dom_id(ib)).ne.rdiv(my_cor)) then
                    print*,'==ERROR==>wrong rdiv for perBCw',dom_id(ib)
                    stop
                 end if
              end if
           end if
!=== East ===> 
           if (dom(ib)%inext.lt.0) then
              if (dom(ib)%iprev.ge.0) then
                 my_cor=dom(ib)%per_in
                 if(rdiv(dom_id(ib)).ne.rdiv(my_cor)) then
                    print*,'==ERROR==>wrong rdiv for perBCe',dom_id(ib)
                    stop
                 end if
              end if
           end if
        end if

        if(dom(ib)%bc_south.eq.5 .or. dom(ib)%bc_north.eq.5) then
!=== South ===> 
           if (dom(ib)%jprev.lt.0) then
              if (dom(ib)%jnext.ge.0) then
                 my_cor=dom(ib)%per_jp  
                 if(rdiv(dom_id(ib)).ne.rdiv(my_cor)) then
                    print*,'==ERROR==>wrong rdiv for perBCs',dom_id(ib)
                    stop
                 end if
              end if
           end if
!=== North ===> 
           if (dom(ib)%jnext.lt.0) then
              if (dom(ib)%jprev.ge.0) then
                 my_cor=dom(ib)%per_jn
                 if(rdiv(dom_id(ib)).ne.rdiv(my_cor)) then
                    print*,'==ERROR==>wrong rdiv for perBCn',dom_id(ib)
                    stop
                 end if
              end if
           end if
        end if

        if(dom(ib)%bc_bottom.eq.5 .or. dom(ib)%bc_top.eq.5) then
!=== Bottom ===> 
           if (dom(ib)%kprev.lt.0) then
              if (dom(ib)%knext.ge.0) then
                 my_cor=dom(ib)%per_kp
                 if(rdiv(dom_id(ib)).ne.rdiv(my_cor)) then
                    print*,'==ERROR==>wrong rdiv for perBCb',dom_id(ib)
                    stop
                 end if
              end if
           end if
!=== Top ===> 
           if (dom(ib)%knext.lt.0) then
              if (dom(ib)%kprev.ge.0) then
                 my_cor=dom(ib)%per_kn
                 if(rdiv(dom_id(ib)).ne.rdiv(my_cor)) then
                    print*,'==ERROR==>wrong rdiv for perBCt',dom_id(ib)
                    stop
                 end if
              end if
           end if
        end if

        end do
!===============================================================
        end if

!        call datainfo

        return
        end 
```
