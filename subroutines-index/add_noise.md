# add\_noise

**File:** averaging.for\
\
**Subroutine Called in:** \
fdstag (fdstag.for - Only at start or restart)\
\
**Purpose:**\
****This subroutine adds random noise based on the initial computational field before the simulation starts.\
\
**User Likeness to alter the file:** \
****Never\


```
!##########################################################################
        subroutine add_noise(fnoise)
!##########################################################################
        use vars
        use multidata
        implicit none
        double precision    :: fnoise
        double precision    :: random_number_normal
        integer :: i,j,k,ib

!......add some gausian noise to flowfield

        call RANDOM_SEED

        do ib=1,nbp

        do k=1,dom(ib)%ttc_k
           do j=1,dom(ib)%ttc_j
              do i=1,dom(ib)%ttc_i
                 if (L_LSMbase) then
		      if (dom(ib)%zc(k).le.length) then
			 dom(ib)%u(i,j,k) = dom(ib)%u(i,j,k) +
     & random_number_normal(0.0,fnoise)
	 	      endif
		     else if (L_LSM) then
			if (dom(ib)%phi(i,j,k) .ge. 0.0) then
                  dom(ib)%u(i,j,k) = dom(ib)%u(i,j,k) +
     & random_number_normal(0.0,fnoise)
			endif
		     else
                  dom(ib)%u(i,j,k) = dom(ib)%u(i,j,k) +
     & random_number_normal(0.0,fnoise)
		     endif
              end do
           end do
        end do

        do k=1,dom(ib)%ttc_k
           do j=1,dom(ib)%ttc_j
              do i=1,dom(ib)%ttc_i
                 if (L_LSMbase) then
		      if (dom(ib)%zc(k).le.length) then
                   dom(ib)%v(i,j,k) = dom(ib)%v(i,j,k) +
     & random_number_normal(0.0,fnoise)
	 	      endif
		     else if (L_LSM) then
			if (dom(ib)%phi(i,j,k) .ge. 0.0) then
                    dom(ib)%v(i,j,k) = dom(ib)%v(i,j,k) +
     & random_number_normal(0.0,fnoise)
			endif
		     else
                  dom(ib)%v(i,j,k) = dom(ib)%v(i,j,k) +
     & random_number_normal(0.0,fnoise)
                 endif
              end do
           end do
        end do

        do k=1,dom(ib)%ttc_k
           do j=1,dom(ib)%ttc_j
              do i=1,dom(ib)%ttc_i
                 if (L_LSMbase) then
		      if (dom(ib)%zc(k).le.length) then
                   dom(ib)%w(i,j,k) = dom(ib)%w(i,j,k) +
     & random_number_normal(0.0,fnoise)
	 	      endif
		     else if (L_LSM) then
			if (dom(ib)%phi(i,j,k) .ge. 0.0) then
                    dom(ib)%w(i,j,k) = dom(ib)%w(i,j,k) +
     & random_number_normal(0.0,fnoise)
			endif
		     else
                  dom(ib)%w(i,j,k) = dom(ib)%w(i,j,k) +
     & random_number_normal(0.0,fnoise)
                 endif
              end do
           end do
        end do

        end do

        end subroutine add_noise
!##########################################################################
        function random_number_normal(mean,sigma) result( fn_val )
!##########################################################################
!       Generate random numbers
!       with a normal distribution with given mean and standard deviaton.
!
!       Generate a random normal deviate using the polar method.
!       Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
!                  normal variables', Siam Rev., vol.6, 260-264, 1964.
!       (source from internet)

        implicit none
        double precision       :: fn_val,mean,sigma
        double precision       :: ull, sumall
        double precision, save :: vll, sln
        logical, save   :: second = .false.
        double precision, parameter :: one = 1.0, vsmall = tiny( one )

        if (second) then

!...... If second, use the second random number generated on last call
           second = .false.
           fn_val = vll*sln

        else
!...... First call; generate a pair of random normals
           second = .true.
           do
              call random_number(ull)
              call random_number(vll)

              ull = scale( ull, 1 ) - one
              vll = scale( vll, 1 ) - one

!.........vsmall added to prevent LOG(zero) / zero
              sumall = ull*ull + vll*vll + vsmall
              if(sumall < one) exit
           end do

           sln = sqrt(- scale( log(sumall), 1 ) / sumall)
           fn_val = ull*sln
        end if

!.....set mean and standart deviation
        fn_val = fn_val * sigma + mean

        return

        end function random_number_normal
!##########################################################################

```



