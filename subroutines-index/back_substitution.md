# Back\_Substitution

**File:** DeltaF\__MLS.for_\
__\
__**Subroutine Called in:** __ \
_forma3Dder (DeltaF\_MLS.for) -> ShapeFunction (DeltaF_\_MLS.for) -> deltah (ibm.for) -> IB_previous (ibm.for) -> flosol (flosol.for)_\
__\
_forma3Dder (DeltaF\_MLS.for) ->  ShapeFunction (DeltaF_\_MLS.for_) -> imb\__openmpi (ibm.for) -> IBM (ibm.for) -> flosol (flosol.for)\
\
**Purpose:**\
****\
****\
******User Likeness to alter the file:** \
**Rarely**&#x20;

```
!######################################################################
	subroutine Back_Substitution(ndim,A,b,x)
!######################################################################
	INTEGER :: i
        INTEGER, intent(in) :: ndim
        Double precision :: A(ndim,ndim),b(ndim),x(ndim)

      x(ndim)=b(ndim)/A(ndim,ndim)
      do i=ndim-1,1,-1
        suma=0
          do j=i+1,ndim
            suma=suma+A(i,j)*x(j)
          enddo
        x(i)=(b(i)-suma)/A(i,i)
      enddo
	return
	end
```
