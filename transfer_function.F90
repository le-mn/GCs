real function transfer_function(k,q,Gamma_eff)
  !
  ! Uses
  use Cosmological_Parameters
  use Power_Spectrum_Parameters
  implicit none
  !
  ! Floats
  real alphav,b1,b2,betac,Bk,C,fb,fc,fcb,fv,fvb,Gammaeff,Gamma_eff,k,kRfw,L,Nv,pc,pcb,q,qeff,qEH,qv&
       &,s,Theta27,Tsup,yd,zd,zeq
  !
  ! Functions
  real transfer_function_BBKS_CDM
  !
  ! Code
  select case (itrans)
  case (1)
     ! BBKS transfer function
     transfer_function=transfer_function_BBKS_CDM(q)
  case (2)
     ! Bond & Efstathiou transfer function
     transfer_function=(1.0+(((6.4*q)+((3.0*q)**1.5)+((1.7*q)**2.0))**1.13))**(-1.0/1.13)
  case (3)
     ! Eisenstein & Hu (1999, ApJ, 511, 5)
     Theta27=CMB_T0/2.7 ! CMB temperature [in units of 2.7K]
     zeq=2.50e4*omega0*(h0**2)/(Theta27**4)
     b1=0.313*((omega0*(h0**2))**(-0.419))*(1.0+0.607*((omega0*(h0**2))**0.674))
     b2=0.238*((omega0*(h0**2))**0.223)
     zd=1291.0*((omega0*(h0**2))**0.251)*(1.0+b1*((omegab*(h0**2))**b2))/(1.0+0.659*((omega0*(h0**2))**0.828))
     yd=(1.0+zeq)/(1.0+zd)
     s=44.5*log(9.83/omega0/(h0**2))/sqrt(1.0+10.0*((omegab*(h0**2))**0.75))
     qEH=k*(Theta27**2)/Gamma_eff
     fv=0.0 ! No neutrinos.
     Nv=0.0
     fb=omegab/omega0
     fc=(omega0-omegab)/omega0
     fcb=fb+fc
     fvb=fv+fb
     pc=0.25*(5.0-sqrt(1.0+24.0*fc))
     pcb=0.25*(5.0-sqrt(1.0+24.0*fcb))
     alphav=(fc/fcb)*((5.0-2.0*(pc+pcb))/(5.0-4.0*pcb))*((1.0-0.533 *fvb+0.126*(fvb**3.))*((1.0+yd)**(pcb-pc))/(1.0-0.193&
          &*sqrt(fv*Nv)+0.169*fv*(Nv**0.2)))*(1.0+0.5*(pc-pcb)*(1.0+1.0/(3.0-4.0*pc)/(7.0-4.0*pcb))/(1.+yd))
     Gammaeff=Gamma_eff*(sqrt(alphav)+(1.0-sqrt(alphav))/(1.0+((0.43*(k*h0)*s)**4)))
     qeff=k*(Theta27**2)/Gammaeff
     betac=1.0/(1.0-0.949*fvb)
     L=log(exp(1.0)+1.84*betac*sqrt(alphav)*qeff)
     C=14.4+325.0/(1.0+60.5*(qeff**1.11))
     Tsup=L/(L+C*(qeff**2))
     if (fv.eq.0.0) then
        qv=0.0
        Bk=1.0
     else
        qv=3.92*qEH*sqrt(Nv)/fv
        Bk=1.0+(1.2*(fv**0.64)*(Nv**(0.3+0.6*fv)))/((qv**(-1.6))+(qv**0.8))
     endif
     transfer_function=Tsup*Bk
  case (10)
     ! BBKS WDM matter transfer function
     kRfw=0.2*q*(omega0*h0**2/mwdm)**1.3333
     transfer_function=exp(-(kRfw/2.0)-((kRfw**2)/2.0))*transfer_function_BBKS_CDM(q)
  case default
     write (0,*) 'transfer_function(): FATAL - transfer function not recognised!'
     write (0,*) '                     itrans = ',itrans
     stop
  end select
  return
end function transfer_function

real function transfer_function_BBKS_CDM(q)
  !
  ! Returns the BBKS CDM transfer function
  implicit none
  !
  ! Floats
  real q
  !
  ! Code
  transfer_function_BBKS_CDM=(log(1.0+2.34*q)/(2.34*q**2))/((1.0/q)**4+3.89/q**3+(16.1/q)**2+5.46**3/q+6.71**4)**0.25
  return
end function transfer_function_BBKS_CDM


subroutine cobe_sigma8
  !
  ! cobe_sigma8():
  !     Return the sigma_8 appropriate for COBE normalization given
  !     Gamma, Omega_0 and Lambda_0.
  !
  !     Uses new fitting formulae from Bunn & White 96 in preparation
  !     1-sigma error bar now quoted at 7.5%.
  !
  !
  ! Uses
  use Cosmological_Parameters
  use Power_Spectrum_Parameters
  implicit none
  !
  ! Integers
  integer i,N
  !
  ! Floats
  real COBE_4YR,delta_H,dlnk,EPS,lnk,lnkmax,lnkmin,pk,pw2k3,pw2k3_kmax,pw2k3_kmin,pwdwk3,rf,sum
  !
  ! Parameters
  parameter (EPS=1.0e-04,N=100000,COBE_4YR=3.1768)
  !
  ! Code
  rf=8.0 ! Top-Hat filter radius in Mpc/h.
  !
  ! Liddle et al 1996 fitting formulae to the Bunn & White 1997 COBE
  ! 4yr data analysis.
  ! The first formula is from Liddle et al., 1996, MNRAS, 282, 281
  ! equation (10) but the leading constant has a different value
  ! due to the choice of units for wavenumber k and the Fourier
  ! transform convention.
  if (abs(lambda0+omega0-1.0).le.EPS) then ! Flat models
     delta_H=COBE_4YR*omega0**(-0.785-0.05*log(omega0)) 
     if (abs(nspec-1.0).ge.EPS) then
        ! The following formulae are Liddle et al., 1996, MNRAS, 282, 281
        ! equations (11) and (12).
        if (igwave.eq.0) then 
           delta_H=delta_H*exp(-0.95*(nspec-1.0)-0.169*(nspec-1)**2)  
        else if (igwave.eq.1) then
           delta_H=delta_H*exp(1.0*(nspec-1.0)+1.97*(nspec-1)**2)  
#ifdef DEBUG
           write (0,*) 'cobe_sigma8(): DEBUG - assuming gravitational waves for power-law inflation'
#endif
           if (nspec.gt.1.0) stop 'cobe_sigma8(): FATAL - assuming gravitational waves for power-law inflation - which is not&
                & applicable for n>1'
        else
#ifdef DEBUG
           write (0,*) 'cobe_sigma8(): FATAL - igwave must be 1 or 0'
#endif
        end if
     end if
  else if (abs(lambda0).le.EPS) then ! Open models
     if (abs(nspec-1.0).ge.EPS) then
#ifdef DEBUG
        write (0,*) 'cobe_sigma8(): DEBUG - assuming nspec=1 see Liddle, Lyth, Roberts and Viana.'
#endif
        nspec=1.0
     end if
     delta_H=COBE_4YR*omega0**(-0.35-0.19*log(omega0)) 
  else
     stop 'cobe_sigma8(): FATAL - cannot cope with this cosmology'
  end if
  !
  ! Integrate k^2 P(k) W(k) and scale by delta_H.
  lnkmax=5.0-log(rf)
  lnkmin=-9.0-log(rf)
  dlnk=(lnkmax-lnkmin)/float(N-1)
  sum=0.0
  do i=1,N
     lnk=lnkmin+dlnk*real(i-1)
     call pkfacs(exp(lnk),rf,Gamma,pk,pw2k3,pwdwk3)
     sum=sum+pw2k3
  end do
  lnkmin=lnkmin-dlnk
  call pkfacs(exp(lnkmin),rf,Gamma,pk,pw2k3_kmin,pwdwk3)
  lnkmax=lnkmax+dlnk 
  call pkfacs(exp(lnkmax),rf,Gamma,pk,pw2k3_kmax,pwdwk3)
  sigma8=(sum+0.5*pw2k3_kmin+0.5*pw2k3_kmax)*dlnk
  sigma8=delta_H*sqrt(sigma8)
  return
end subroutine cobe_sigma8
