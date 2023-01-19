# tauw\_noslip

**File:** eddyvis\_smag.for\
\
**Subroutine Called in:** \
**Not called anywhere**\
\
**Purpose:**\
****\
****\
******User Likeness to alter the file:** \
****Never (Because it has already been widely validated)

```
!##########################################################################
        subroutine tauw_noslip
!##########################################################################
        use vars
        use multidata
        implicit none
        integer i,j,k,ib
        double precision delta,n_x,n_y,n_z,vnor,vtan
        double precision uc,vc,wc,rrey,small


        rrey=1.0/Re
        small = 1.e-30

        do ib=1,nbp

        dom(ib)%tauww=0.0; dom(ib)%tauwe=0.0
        dom(ib)%tauwn=0.0; dom(ib)%tauws=0.0
        dom(ib)%tauwt=0.0; dom(ib)%tauwb=0.0

        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%bc_west.eq.4) then						!if it's 6 it was calculated before
              delta=0.5*dom(ib)%dx
              n_x=1.0
              n_y=0.0
              n_z=0.0
              i=dom(ib)%isp
              do k=dom(ib)%ksp,dom(ib)%kep
                 do j=dom(ib)%jsp,dom(ib)%jep
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                    dom(ib)%tauww(j,k) = rrey * vtan / delta
                 end do
              end do
           end if
        end if

        if (dom(ib)%inext.lt.0) then
           if (dom(ib)%bc_east.eq.4) then
              delta=0.5*dom(ib)%dx
              n_x=-1.0
              n_y=0.0
              n_z=0.0
              i=dom(ib)%iep
              do k=dom(ib)%ksp,dom(ib)%kep
                 do j=dom(ib)%jsp,dom(ib)%jep
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                    dom(ib)%tauwe(j,k) = rrey * vtan / delta
                 end do
              end do
           end if
        end if

        if (dom(ib)%jprev.lt.0) then
           if (dom(ib)%bc_south.eq.4) then
              delta=0.5*dom(ib)%dy
              n_x=0.0
              n_y=1.0
              n_z=0.0
              j=dom(ib)%jsp
              do k=dom(ib)%ksp,dom(ib)%kep
                 do i=dom(ib)%isp,dom(ib)%iep
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                    dom(ib)%tauws(i,k) = rrey * vtan / delta
                 end do
              end do
           end if
        end if

        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%bc_north.eq.4) then
              delta=0.5*dom(ib)%dy
              n_x=0.0
              n_y=-1.0
              n_z=0.0
              j=dom(ib)%jep
              do k=dom(ib)%ksp,dom(ib)%kep
                 do i=dom(ib)%isp,dom(ib)%iep
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                    dom(ib)%tauwn(i,k) = rrey * vtan / delta
                 end do
              end do
           end if
        end if

        if (dom(ib)%kprev.lt.0) then
           if (dom(ib)%bc_bottom.eq.4) then
              delta=0.5*dom(ib)%dz
              n_x=0.0
              n_y=0.0
              n_z=1.0
              k=dom(ib)%ksp
              do j=dom(ib)%jsp,dom(ib)%jep
                 do i=dom(ib)%isp,dom(ib)%iep
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                    dom(ib)%tauwb(i,j) = rrey * vtan / delta
                 end do
              end do
           end if
        end if

        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%bc_top.eq.4) then
              delta=0.5*dom(ib)%dz
              n_x=0.0
              n_y=0.0
              n_z=-1.0
              k=dom(ib)%kep
              do j=dom(ib)%jsp,dom(ib)%jep
                 do i=dom(ib)%isp,dom(ib)%iep
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                    dom(ib)%tauwt(i,j) = rrey * vtan / delta
                 end do
              end do
           end if
        end if

        end do

        return
        end subroutine tauw_noslip
!##########################################################################
```
