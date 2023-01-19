# format3Dder

**File:** DeltaF\__MLS.for_\
__\
__**Subroutine Called in:** __ \
_ShapeFunction (DeltaF_\_MLS.for) -> deltah (ibm.for) -> IB_previous (ibm.for) -> flosol (flosol.for)_\
_ShapeFunction (DeltaF_\_MLS.for_) -> imb\__openmpi (ibm.for) -> IBM (ibm.for) -> flosol (flosol.for)\
\
**Purpose:**\
****\
****\
******User Likeness to alter the file:** \
**Rarely**&#x20;

```
!######################################################################
      subroutine forma3Dder(FF,W2,dX2,dY2,dZ2,neignum,Xmls,Ymls,Zmls
     &  ,xgau,ygau,zgau,hx,hy,hz,
     &   dFFX,dFFY,dFFZ,ddFFXX,ddFFXY,ddFFXZ,ddFFYY,ddFFYZ,ddFFZZ,
     &   MAXNEIG,nbase)
!######################################################################
        implicit none
        INTEGER, intent(in) :: neignum,MAXNEIG,nbase
        Double precision, intent(in) :: hx,hy,hz,xgau,ygau,zgau
        Double precision, intent(in) :: W2(MAXNEIG)
	Double precision, intent(in) :: Xmls(MAXNEIG),Ymls(MAXNEIG)
	Double precision, intent(in) :: Zmls(MAXNEIG)
	Double precision, intent(in) :: dX2(MAXNEIG),dY2(MAXNEIG)
	Double precision, intent(in) :: dZ2(MAXNEIG)
        INTEGER :: I,J,K,ifila,icolumna,ibase,nelim,icol
        Double precision :: xg,yg,zg,wmax,wmin,valor
        Double precision :: FF(MAXNEIG),dFFX(MAXNEIG),dFFY(MAXNEIG)
        Double precision :: dFFZ(MAXNEIG)
        Double precision :: ddFFXX(MAXNEIG),ddFFXY(MAXNEIG)
        Double precision :: ddFFYY(MAXNEIG)
        Double precision :: ddFFXZ(MAXNEIG),ddFFYZ(MAXNEIG)
        Double precision :: ddFFZZ(MAXNEIG)
	Double precision :: A(nbase,nbase),vIA(nbase,nbase)
        Double precision :: eye(nbase,nbase)
        Double precision :: C(nbase,MAXNEIG),PtC(MAXNEIG,MAXNEIG)
        Double precision :: auxiliar(MAXNEIG,nbase)
	Double precision :: auxiliar1(nbase),auxiliar2(nbase)
        Double precision :: sol(nbase,nbase),wsv(10),vsv(10,10),coef

c Inicializamos la matriz C
        do i=1,nbase
	 do j=1, MAXNEIG
	  C(i,j)=0.d+00
	 enddo
	enddo
c Construimos la matriz A

      do i=1,nbase
	do j=1,nbase
	 A(i,j)= 0.d+00
	 vIA(i,j)=0.d+00
	 if(i.eq.j)then
	  eye(i,j)= 1.d+00
	 else
	  eye(i,j)= 0.d+00
	 endif
	enddo
      enddo

c Construimos un vector auxiliar para hacer los calculos
      do k=1,neignum
       xg=(Xmls(k)-xgau)/hx ; yg=(Ymls(k)-ygau)/hy ;zg=(Zmls(k)-zgau)/hz

      if(nbase.eq.4) then
      	auxiliar(k, 1)= 1.d+00
	auxiliar(k, 2)= xg
	auxiliar(k, 3)= yg
	auxiliar(k, 4)= zg
      elseif(nbase.eq.10) then
	auxiliar(k, 1)= 1.d+00
	auxiliar(k, 2)= xg
	auxiliar(k, 3)= yg
	auxiliar(k, 4)= zg
	auxiliar(k, 5)= xg*xg
	auxiliar(k, 6)= xg*yg
	auxiliar(k, 7)= xg*zg
	auxiliar(k, 8)= yg*yg
	auxiliar(k, 9)= yg*zg
	auxiliar(k,10)= zg*zg
      elseif(nbase.eq.20) then
	auxiliar(k, 1)= 1.d+00
	auxiliar(k, 2)= xg
	auxiliar(k, 3)= yg
	auxiliar(k, 4)= zg
	auxiliar(k, 5)= xg*xg
	auxiliar(k, 6)= xg*yg
	auxiliar(k, 7)= xg*zg
	auxiliar(k, 8)= yg*yg
	auxiliar(k, 9)= yg*zg
	auxiliar(k,10)= zg*zg
	auxiliar(k,11)= xg*xg*xg
	auxiliar(k,12)= xg*xg*yg
	auxiliar(k,13)= xg*xg*zg
	auxiliar(k,14)= xg*yg*zg
	auxiliar(k,15)= yg*yg*yg
	auxiliar(k,16)= yg*yg*xg
	auxiliar(k,17)= yg*yg*zg
	auxiliar(k,18)= zg*zg*zg
	auxiliar(k,19)= zg*zg*xg
	auxiliar(k,20)= zg*zg*yg
       endif

       do ifila= 1,nbase  
	do icolumna= 1,nbase
     	 A(ifila,icolumna)= A(ifila,icolumna)+
     &         W2(k)*auxiliar(k,ifila)*auxiliar(k,icolumna)
	enddo
       enddo
      enddo



c Calculamos la inversa de A
c Para calcularla resolvemos 6 sistemas de ecuaciones
c Primero triangularizamos la matriz A y hacemos las operaciones por filas
c sobre la matriz identidad
      do i=1,nbase-1
c Comprobamos que el pivote no sea nulo
        if(A(i,i).lt.1.d-06)then
	  write(6,*)' Pivote nulo',hx,hy,hz,i
	  read(*,*)
	endif
       do j=i+1,nbase
        coef=-A(j,i)/A(i,i)
          do k=i,nbase
            A(j,k)=A(j,k)+coef*A(i,k)
          enddo
c Aplicamos la operación a las columnas de la identidad
	    do icol=1,nbase
            eye(j,icol)= eye(j,icol)+coef*eye(i,icol)
	    enddo
        enddo
      enddo
      do icol=1,nbase
	do i=1,nbase
	 auxiliar1(i)= eye(i,icol)
	enddo
         call Back_Substitution(nbase,A,auxiliar1,auxiliar2)
	do i=1,nbase
	 vIA(i,icol)= auxiliar2(i)
	enddo
       enddo


c We calculate C, which is A^-1*B
      do i=1,nbase
	do k=1, neignum
 	 C(i,k)= 0.d+00
	  do j=1,nbase  
	  C(i,k)= C(i,k) + vIA(i,j)*auxiliar(k,j)
	  enddo
	 C(i,k)= C(i,k)*W2(k)
	enddo
      enddo
c We calculate I-Pt*C
      do i=1,neignum
	do j=1,neignum

	valor= 0.d+00
	do ibase=1,nbase
		valor= valor + C(ibase,j)*auxiliar(i,ibase)
	enddo

		PtC(i,j)=-valor
		if(i.eq.j)PtC(i,j)=PtC(i,j)+1.d+00
	  enddo
	enddo

      if(nbase.eq.4) then
        do k=1,neignum
	   FF(k)= C(1,k)
	enddo
        do j=1,neignum
	 dFFX(j)=C(2,j)/hx
	 dFFY(j)=C(3,j)/hy
	 dFFZ(j)=C(4,j)/hz
	 do k=1,neignum
           dFFX(j)=dFFX(j)+dX2(k)*C(1,k)*PtC(k,j)/W2(k)
           dFFY(j)=dFFY(j)+dY2(k)*C(1,k)*PtC(k,j)/W2(k)
           dFFZ(j)=dFFZ(j)+dZ2(k)*C(1,k)*PtC(k,j)/W2(k)
	 enddo
	enddo
      elseif(nbase.eq.10) then
        do k=1,neignum
 	  FF(k)= C(1,k)
	enddo
        do j=1,neignum
	  dFFX(j)=C(2,j)/hx
	  dFFY(j)=C(3,j)/hy
	 dFFZ(j)=C(4,j)/hz
            do k=1,neignum
             dFFX(j)=dFFX(j)+dX2(k)*C(1,k)*PtC(k,j)/W2(k)
             dFFY(j)=dFFY(j)+dY2(k)*C(1,k)*PtC(k,j)/W2(k)
             dFFZ(j)=dFFZ(j)+dZ2(k)*C(1,k)*PtC(k,j)/W2(k)
            enddo
	enddo
        do j=1,neignum
	  ddFFXY(j)=        C(5,j)/(hx*hy)
	  ddFFXZ(j)=        C(6,j)/(hx*hz)
	  ddFFYZ(j)=        C(7,j)/(hy*hz)
	  ddFFXX(j)= 2.d+00*C(8,j)/(hx*hx)  
	  ddFFYY(j)= 2.d+00*C(9,j)/(hy*hy)
	  ddFFZZ(j)=2.d+00*C(10,j)/(hz*hz)
        enddo
      elseif(nbase.eq.20) then
        do k=1,neignum
	  FF  (k)= C(1,k)
	enddo
        do j=1,neignum
	  dFFX(j)=C(2,j)/hx
	  dFFY(j)=C(3,j)/hy
	  dFFZ(j)=C(4,j)/hZ
	 do k=1,neignum
           dFFX(j)=dFFX(j)+dX2(k)*C(1,k)*PtC(k,j)/W2(k)
           dFFY(j)=dFFY(j)+dY2(k)*C(1,k)*PtC(k,j)/W2(k)
           dFFZ(j)=dFFZ(j)+dZ2(k)*C(1,k)*PtC(k,j)/W2(k)
	 enddo
	enddo
         do j=1,neignum
	  ddFFXX(j)= 2.d+00*C( 5,j)/(hx*hx)
	  ddFFXY(j)=        C( 6,j)/(hx*hy)
	  ddFFXZ(j)=        C( 7,j)/(hx*hz)
	  ddFFYY(j)= 2.d+00*C( 8,j)/(hy*hy)
	  ddFFYZ(j)=        C( 9,j)/(hy*hz)
	  ddFFZZ(j)= 2.d+00*C(10,j)/(hz*hz)
         enddo
	endif

	return
	end
```
