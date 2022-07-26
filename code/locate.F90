subroutine locate(XX,N,X,J)
  implicit none
  !
  ! Integers
  integer J,JL,JM,JU,N
  !
  ! Floats
  real X,XX(N)
  !
  ! Code
  JL=0
  JU=N+1
  do while (JU-JL.gt.1)
     JM=(JU+JL)/2
     if ((XX(N).gt.XX(1)).eqv.(X.gt.XX(JM))) then
        JL=JM
     else
        JU=JM
     endif
  end do
  J=JL
  return
end subroutine locate
