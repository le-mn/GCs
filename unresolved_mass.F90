! Integral related to the mass fraction in halos
! below the imposed resolutiom limit m_min.
!
! Here gamma_1 (set in module Modified_Merger_Tree) is
! a parameter involved in modifying the merger rate
! expression away from the Press-Schechter form. 
!
! z = sqrt(sigma^2(m_2)/(sigma^2(m_min)-sigma^2(m_2)))
!
! For the Press-Schecter case gamma_1=0 and J(z)= z
!
! In general J(gamma_1,z)= int_0^z (1+1/u^2)^(gamma_1/2) du
!
! Note: zmax should be >>1 so that dJ/dz= 1 is accurate for z>zmax.
!       dz=zmax/NTAB should be small so that the trapezium rule is
!       accurate for z>dz.
!
real function J_UNRESOLVED(z)
use Modified_Merger_Tree
real, intent(in) :: z
real, parameter :: EPS=1.0e-05,zmax=10.0
real :: dz
integer :: i
integer, parameter :: NTAB=1000
integer, save :: ifirst=0
real, save :: J_tab(NTAB),z_tab(NTAB),inv_dz

if (abs(gamma_1).gt.EPS) then !gamma_1.ne.0
   if (ifirst.eq.0) then !on the first call tabulate the integral.
!     Assuming dz<<1 we can do the integral from 0 to dz analytically
      dz=zmax/real(NTAB)
      inv_dz=1.0/dz
      if (abs(1.0-gamma_1).gt.EPS) then !gamma_1.ne.1
         J_tab(1)= (dz**(1.0-gamma_1))/(1.0-gamma_1)
      else  !special case if gamma_1=1 
         J_tab(1)= log(dz)
      end if   
      z_tab(1)=dz
!     Now continue to higher z by direct summation using the trapezium rule
      do i=2,NTAB
         z_tab(i)=real(i)*dz
         J_tab(i) = J_tab(i-1)+(1.0+1.0/z_tab(i)**2)**(0.5*gamma_1) * 0.5*dz  &
&                             +(1.0+1.0/z_tab(i-1)**2)**(0.5*gamma_1) *0.5*dz 
      end do   
      ifirst=1
   end if   
!
!  Look up using tabulated values and analytic extensions
   i=int(z*inv_dz)
   if (i.lt.1) then !use analytic result assuming z<<1
      if (abs(1.0-gamma_1).gt.EPS) then !gamma_1.ne.1
         J_UNRESOLVED= (z**(1.0-gamma_1))/(1.0-gamma_1)
      else  !special case if gamma_1=1 
         J_UNRESOLVED= log(z)
      end if   
   else if (i.ge.NTAB) then !if beyond tabulated region assume the z>>1
!                            analytic approimation that the integrand is unity
!                            and hence J grows linearly with z
      J_UNRESOLVED=J_tab(NTAB)+z-z_tab(NTAB) 
   else !linearly interpolate between tabulated values
      h=(z-z_tab(i))*inv_dz
      J_UNRESOLVED=J_tab(i)*(1.0-h) + J_tab(i+1)*h
   end if   
else !gamma_1=0
   J_UNRESOLVED=z  !Press-Schechter case
end if
end function J_UNRESOLVED
