module Modified_Merger_Tree
  implicit none
  save
  ! Parameters used to modify the merger rate used in split.F90
  ! to be slightly different to the standard Press-Schechter formula
  !
  ! The modify factor is
  ! G0 [sigma(m1)/sigma(m2)]^gamma_1 [w/sigma(m2)]^gamma_2
  !
  !real, parameter :: gamma_1=0.02,gamma_2=0.1,G0=0.82
  !real :: gamma_1,gamma_2,G0=0.0
  real :: G0=0.57, gamma_1=0.38, gamma_2=-0.01    !using the same values for G0, gamma1, gamma2 as in the tress.f90   

end module Modified_Merger_Tree
