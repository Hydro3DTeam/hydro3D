# exchange

**File:** exchange.for\
\
**Subroutine Called in:** \
**flosol** (flosol.for)\
**calvel** (press.for) -> **pressure\_1sweep** (press.for) -> **flosol** (flosol.for)\
**calvel** (press.for) -> **newsolv\_mg** (newsolvmg.for) -> **flosol** (flosol.for)\
**newsolv\_mg** (newsolvmg.for) -> **flosol** (flosol.for)\
**IBM** (ibm.for) -> **flosol** (flosol.for)\
**sipsol** (sipsol.for) -> **flosol** (flosol.for)\
**eddyv\_smag** (eddyvis\_smag.for) -> **flosol** (flosol.for)\
**eddyv\_wale** (eddyvis\_wale.for) -> **flosol** (flosol.for)\
**eddyv\_1eqn** (eddyvis\_1eqn.for) -> **flosol** (flosol.for)\
**eddyv\_keps** (eddyvis\_keps.for) -> **flosol** (flosol.for)\
**exchange\_Scal** (scalar.for) ->  **scalar\_cal** (scalar.for) -> **flosol** (flosol.for)\
**bound\_LSM** (bounds\_lsm.for) ->  **lsm\_3d** (lsm.for) ->  **flosol** (flosol.for)\
**bound\_LSM** (bounds\_lsm.for) ->  **tvd\_rk\_reinit** (lsm.for) -> **initial\_lsm3d\_channel** (lsm.for) -> **fdstag** (fdstag.for)\
**bound\_LSM** (bounds\_lsm.for) -> **tvd\_rk\_renit** (lsm.for) -> **lsm\_3d** (lsm.for) -> **flosol** (flosol.for)\
**bound\_LSM** (bounds\_lsm.for) -> **heaviside** (lsm.for) -> **flosol** (flosol.for)\
\
**Purpose:**\
****This subroutine ensures that the boundary values of most variables are transferred to the neighbours' ghost cells.\
\
**User Likeness to alter the file:** \
****Rarely

```
##########################################################################
        subroutine  exchange(op)
!##########################################################################
        use vars
        use multidata
        implicit none
        integer :: op,ly,ib,ni,nj,nk,ijk,nly
        integer :: i,j,k,is,ie,js,je,ks,ke
        integer :: ispr,iepr,jspr,jepr,kspr,kepr
        double precision, pointer, dimension(:,:,:) :: fi

! stfcinf: 1=iprev, 3=jprev, 5=kprev, 2=inext, 4=jnext, 6=knext

        if (op.eq.4 .or. op.eq.44 .or. op.eq.7) then! .or. op.eq.14 .or.
!     & op.eq.15 .or. op.eq.16 .or. op.eq.17 .or. op.eq.18 .or. 
!     & op.eq.19) then
           nly=0
        else 
           nly=pl_ex
        end if

        if (PERIODIC) call exchange_bc(op,nly)
        call exchange_smlvl(op,nly)

        if(rdivmax.gt.1) then

        select case (op)
           case (1)  
              call exchangeu(op,nly)
           case (2)  
              call exchangev(op,nly)
           case (3)  
              call exchangew(op,nly)
           case (11)  
              call exchangeu(op,nly)
           case (22)  
              call exchangev(op,nly)
           case (33)  
              call exchangew(op,nly)
           case (4)  
              call exchangep(op)
           case (44)  
              call exchangep(op)
           case (5)  
              call exchangesca(op,nly)
           case (6)  
              call exchangesca(op,nly)
           case (7)  
              call exchangesca(op,nly)
           case (8)
              call exchangesca(op,nly)
           case (9)  
              call exchangesca(op,nly)
	     case (14)
		  call exchange_phi(op,nly) !phi
	     case (15)
		  call exchange_phi(op,nly) !phi_reinit
	     case (16)
		  call exchange_phi(op,nly) !phi_new
	     case (17)
		  call exchange_phi(op,nly) !phi_init
	     case (18)
		  call exchange_phi(op,nly)	!dens
	     case (19)
		  call exchange_phi(op,nly)	!mu
        end select

        if(LMR.eq.2) then 

        do ib=1,nbp

           ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k  

           select case (op)
              case (1)  
                 fi => dom(ib)%u
                 is=dom(ib)%isu; ie=dom(ib)%ieu
                 js=dom(ib)%jsu; je=dom(ib)%jeu
                 ks=dom(ib)%ksu; ke=dom(ib)%keu
              case (2)  
                 fi => dom(ib)%v
                 is=dom(ib)%isv; ie=dom(ib)%iev
                 js=dom(ib)%jsv; je=dom(ib)%jev
                 ks=dom(ib)%ksv; ke=dom(ib)%kev
              case (3)  
                 fi => dom(ib)%w
                 is=dom(ib)%isw; ie=dom(ib)%iew
                 js=dom(ib)%jsw; je=dom(ib)%jew
                 ks=dom(ib)%ksw; ke=dom(ib)%kew
              case (11)  
                 fi => dom(ib)%ustar
                 is=dom(ib)%isu; ie=dom(ib)%ieu
                 js=dom(ib)%jsu; je=dom(ib)%jeu
                 ks=dom(ib)%ksu; ke=dom(ib)%keu
              case (22)  
                 fi => dom(ib)%vstar
                 is=dom(ib)%isv; ie=dom(ib)%iev
                 js=dom(ib)%jsv; je=dom(ib)%jev
                 ks=dom(ib)%ksv; ke=dom(ib)%kev
              case (33)  
                 fi => dom(ib)%wstar
                 is=dom(ib)%isw; ie=dom(ib)%iew
                 js=dom(ib)%jsw; je=dom(ib)%jew
                 ks=dom(ib)%ksw; ke=dom(ib)%kew
              case (4) 
                 fi => dom(ib)%p
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
              case (44) 
                 fi => dom(ib)%pp
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
              case (5) 
                 fi => dom(ib)%T
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
              case (6) 
                 fi => dom(ib)%ksgs
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
              case (7) 
                 fi => dom(ib)%vis
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
              case (8)
                 fi => dom(ib)%Scal
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep  
              case (9)
                 fi => dom(ib)%eps
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep	
	        case (14)
		     fi => dom(ib)%phi
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
	        case (15)
		     fi => dom(ib)%phi_reinit
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
	        case (16)
		     fi => dom(ib)%phi_new
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
	        case (17)
		     fi => dom(ib)%phi_init
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
	        case (18)
		     fi => dom(ib)%dens
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
	        case (19)
		     fi => dom(ib)%mu
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
           end select

           do ly=0,nly

!..............................................................................
!=== Previous Neighbor  ===> 
!..............................................................................
              if (dom(ib)%iprev.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%iprev)) then
!===>> iprev !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ispr=pl+1;	jspr=pl+1;	kspr=pl+1
                 iepr=ni-pl;	jepr=nj-pl;	kepr=nk-pl
                 if (op.eq.2 .or. op.eq.22) then
        if(dom(ib)%edgprev5.ge.0 .and. 
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev5)) jspr=pl
                 else if (op.eq.3 .or. op.eq.33) then
        if(dom(ib)%edgprev1.ge.0 .and. 
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev1)) kspr=pl
                 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 do k=kspr,kepr; do j=jspr,jepr; ijk=(k-1)*nj+j
                    fi(is-1-ly,j,k)=dom(ib)%stfcinf(1,ly+1,ijk)
                 end do; end do
              end if
              end if

              if (dom(ib)%jprev.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%jprev)) then
!===>> jprev !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ispr=pl+1;	jspr=pl+1;	kspr=pl+1
                 iepr=ni-pl;	jepr=nj-pl;	kepr=nk-pl
                 if (op.eq.1 .or. op.eq.11) then
        if(dom(ib)%edgprev5.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev5)) ispr=pl
                 else if (op.eq.3 .or. op.eq.33) then
        if(dom(ib)%edgprev4.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev4)) kspr=pl
                 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 do k=kspr,kepr; do i=ispr,iepr; ijk=(k-1)*ni+i
                    fi(i,js-1-ly,k)=dom(ib)%stfcinf(3,ly+1,ijk)
                 end do; end do
              end if
              end if

              if (dom(ib)%kprev.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%kprev)) then
!===>> kprev !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ispr=pl+1;	jspr=pl+1;	kspr=pl+1
                 iepr=ni-pl;	jepr=nj-pl;	kepr=nk-pl
                 if (op.eq.1 .or. op.eq.11) then
        if(dom(ib)%edgprev1.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev1)) ispr=pl
                 else if (op.eq.2 .or. op.eq.22) then
        if(dom(ib)%edgprev4.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev4)) jspr=pl
                 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 do i=ispr,iepr; do j=jspr,jepr; ijk=(i-1)*nj+j
                    fi(i,j,ks-1-ly)=dom(ib)%stfcinf(5,ly+1,ijk)
                 end do; end do
              end if
              end if
!..............................................................................
!=== Next Neighbor  ===> 
!..............................................................................
              if (dom(ib)%inext.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%inext)) then
!===>> inext !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ispr=pl+1;	jspr=pl+1;	kspr=pl+1
                 iepr=ni-pl;	jepr=nj-pl;	kepr=nk-pl
                 if (op.eq.2 .or. op.eq.22) then
        if(dom(ib)%edgnext6.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext6)) jspr=pl
                 else if (op.eq.3 .or. op.eq.33) then
        if(dom(ib)%edgprev3.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev3)) kspr=pl
                 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 do k=kspr,kepr; do j=jspr,jepr; ijk=(k-1)*nj+j
                    fi(ie+1+ly,j,k)=dom(ib)%stfcinf(2,ly+1,ijk)
                 end do; end do
              end if
              end if

              if (dom(ib)%jnext.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%jnext)) then
!===>> jnext !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ispr=pl+1;	jspr=pl+1;	kspr=pl+1
                 iepr=ni-pl;	jepr=nj-pl;	kepr=nk-pl
                 if (op.eq.1 .or. op.eq.11) then
        if(dom(ib)%edgprev6.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev6)) ispr=pl
                 else if (op.eq.3 .or. op.eq.33) then
        if(dom(ib)%edgprev2.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev2)) kspr=pl
                 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 do k=kspr,kepr; do i=ispr,iepr; ijk=(k-1)*ni+i
                    fi(i,je+1+ly,k)=dom(ib)%stfcinf(4,ly+1,ijk)
                 end do; end do
              end if
              end if

              if (dom(ib)%knext.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%knext)) then
!===>> knext !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ispr=pl+1;	jspr=pl+1;	kspr=pl+1
                 iepr=ni-pl;	jepr=nj-pl;	kepr=nk-pl
                 if (op.eq.1 .or. op.eq.11) then
        if(dom(ib)%edgnext3.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext3)) ispr=pl
                 else if (op.eq.2 .or. op.eq.22) then
        if(dom(ib)%edgnext2.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext2)) jspr=pl
                 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 do i=ispr,iepr; do j=jspr,jepr; ijk=(i-1)*nj+j
                    fi(i,j,ke+1+ly)=dom(ib)%stfcinf(6,ly+1,ijk)
                 end do; end do
              end if
              end if

           end do
        end do
        end if
        end if

        return
        end subroutine exchange
```
