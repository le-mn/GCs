subroutine spline(X,Y,N,YP1,YPN,Y2)
  implicit none
  !
  ! Array dimensions
  integer N,NMAX
  parameter (NMAX=4096)
  !
  ! Integers
  integer I,K
  !
  ! Floats
  real P,QN,SIG,UN,U(NMAX),X(N),Y2(N),Y(N),YP1,YPN
  !
  ! Code
  if (N.gt.NMAX) stop 'spline(): FATAL - NMAX too small'
  if (YP1.gt.0.99e30) then
     Y2(1)=0.0
     U(1)=0.0
  else
     Y2(1)=-0.5
     U(1)=(3.0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
  end if
  do I=2,N-1
     SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
     P=SIG*Y2(I-1)+2.0
     Y2(I)=(SIG-1.0)/P
     U(I)=(6.0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))/(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
  end do
  if (YPN.gt.0.99e30) then
     QN=0.0
     UN=0.0
  else
     QN=0.5
     UN=(3.0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
  end if
  Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.0)
  do K=N-1,1,-1
     Y2(K)=Y2(K)*Y2(K+1)+U(K)
  end do
  return
end subroutine spline
