!
! Hyperbolic trigonometric functions


pure real function acosh(y)
  ! Returns the positive branch    
  implicit none
  !
  ! Floats
  real, intent(in) :: y
  !
  ! Code
  acosh=log(y+sqrt(y*y-1.0))
  return
end function acosh

pure real function atanh(y)
  implicit none
  !
  ! Floats
  real, intent(in) ::  y
  !
  ! Code
  atanh=0.5*log((1.0+y)/(1.0-y))
  return
end function atanh

pure real function asinh(y)
  implicit none
  !
  ! Floats
  real, intent(in) ::  y
  !
  ! Code 
  asinh=log(y+sqrt(y*y+1.0))
  return
end function asinh

