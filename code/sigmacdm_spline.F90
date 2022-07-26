!     This is the top-level routine for computing sigma(M).
!     It returns both sigma and alpha = dln(sigma)/dlnM.
!     
!     It handles various cases, depending on value of itrans:
!
!     itrans<0: uses tabulated P(k) read in from file
!     itrans=0: uses power-law P(k) propto k^nspec
!     1<= itrans <= 9: uses analytical CDM power spectrum 
!                      P(k) propto k^nspec T(k)^2
!     itrans=10: uses analytical WDM power spectrum 
!                      P(k) propto k^nspec T(k)^2
!     
!     For the analytical CDM & WDM cases (itrans>0) it uses the variable Gamma.
!     
!     In all cases:
!
!     The input mass M is assumed to be in units Msun/h
!     sigma(M) is normalized to the sigma8 read in by parameters()
!     (overriding any normalization e.g. in an input file)
!     
!     On the first call, it calculates coefficients for a spline fit to
!     sigma(M), and also scaling factors sclm and scla for the mass and
!     normalization.
!     
real function sigmacdm(m,alpha)
  !
  ! Uses
  use Cosmological_Parameters
  use Numerical_Parameters
  use Power_Spectrum_Parameters
  implicit none
  !
  ! Integers
  integer fileend,ifirst,i
  !
  ! Floats
  real alpha,ktab,m,ms,m8,pktab,sigma
  !
  ! Saves
  save ifirst
  !
  ! Data
  data ifirst /1/
  !
  ! Code
  !
  ! If first call compute the appropriate scaling factors sclm and
  ! scla calculate sigma in an 8 Mpc/h sphere for the arbitrarily
  ! normalized input P(k), and hence get the factor scla by which the
  ! un-normalized sigma(M) must be multiplied to get the desired
  ! sigma8 normalization for CDM or WDM with an analytical transfer
  ! function, also calculate the scaling factor in mass sclm that
  ! relates the true P(k) to that for Omega0=1, Gamma=1.
  if (ifirst.eq.1.or.ireset.eq.1) then
     m8=M8CRIT*omega0 ! The mass within an 8Mpc/h sphere.
     transfer_function: select case (itrans)
     case (:-1) ! Read P(k) from input file.
#ifdef INFO
        write (0,*) 'sigmacdm(): INFO - reading P(k) from file ',trim(pkinfile)
#endif
        open (77,file=pkinfile,status='old') 
        read (77,*) ! Skip 1 line header.
        i=0 
        fileend = 0
        do while (fileend.eq.0)
           ! File contains k in units h/Mpc & P(k) in units (Mpc/h)^3
           read (77,*,iostat=fileend) ktab,pktab
           if (fileend.eq.0) then
              i=i+1
              if (i.gt.NKTABMAX) stop 'sigamcdm(): FATAL - increase NKTABMAX'
              ! Store k & P(k) as logs for interpolation.
              lnktab(i)=log(ktab)
              lnpktab(i)=log(pktab)
           endif
        enddo
        close(77) 
        nktab=i
#ifdef INFO
        write (0,*) 'sigmacdm(): INFO - number of wavenumbers read in =',nktab
#endif
        ! Mass scaling not used in this case.
        sclm=1.0
        !     calc sigma8 for input spectrum
        call spline_interp(m8,sigma,alpha)
        !     
     case (0) ! Power-law P(k)
        m8=M8CRIT*omega0 ! The mass within an 8Mpc/h sphere.
        sclm=1.0/m8
        ms=m8*sclm         
        sigma=ms**(-(nspec+3.0)/6.0)
     case default ! Analytic CDM or WDM P(k)
        ! Compute the required scaling factors sclm and scla.
        sclm=gamma**3/omega0
        m8=M8CRIT*omega0 ! The mass within an 8Mpc/h sphere.
        ms=m8*sclm
        ! The spline fit to CDM for Gamma=1.
        call spline_interp(ms,sigma,alpha)
     end select transfer_function
     scla=sigma8/sigma      !scales sigma_8 to required value
     !     
     ifirst=0               !indicates first call complete and sclm and scla are set
     ireset=0
  end if
  !     ----------------------------------------------------
  ms=m*sclm
  select case (itrans)
  case (0) ! Power-law P(k)
     alpha=-(nspec+3.0)/6.0
     sigmacdm=scla*(ms**alpha)
  case default ! CDM or WDM or tabulated
     !     Use spline fit
     call spline_interp(ms,sigma,alpha)
     sigmacdm=sigma*scla
  end select
  !     
  return
end function sigmacdm

!
! The spline evaluation called by sigmacdm_spline()
!
! If this first time it is called, it checks to see if a file
! containing the spline coefficients for the current P(k) already
! exists. if it does, it reads this file, otherwise it generates the
! spline coefficients afresh.
!
! Saves last NMOD positions in look-up table along with
! corresponding values of h h2 and inv. This is optimal if every
! NMOD calls the value of ms is approximately repeated. NMOD=2 is
! optimal when used in conjunction with split() as most calls don't
! produce a split and so m changes slowly and one keeps evaluating
! sigma(m) and sigma(m/2).

subroutine spline_interp(ms,sigma,alpha)
  !
  ! Uses
  use Cosmological_Parameters
  use Numerical_Parameters
  use Power_Spectrum_Parameters
  implicit none
  !
  ! Integers
  integer i,ifirst,imod,io,k,khi,klo,NMOD,NSPLINE
  !
  ! Array dimensions
  parameter(NMOD=2)
  !
  ! Integer arrays
  integer kphi(NMOD),kplo(NMOD)
  !
  ! Floats
  real a2(NSPL),a3,aa,alpha,a(NSPL),b3,bb,h,h2,hp2(NMOD),hp(NMOD),invh,invhp(NMOD),m(NSPL),ms,s2(NSPL),sigma,s(NSPL)
  !
  ! Characters
  character sform*10
  !
  ! Saves
  save a,a2,hp,hp2,ifirst,imod,invhp,kphi,kplo,m,s,s2
  !
  ! Data
  data ifirst /1/
  !
  ! Code
  !
  ! On first call set initial values and read the spline fit
  if ( ifirst.eq.1) then
     kplo=NSPL
     kphi=1
     m(1)=0.0 ! Set to here to avoid compiler warning but read below.
     imod=0
     nspec=0.01*nint(100.0*nspec)
     write (sform,'(sp,f6.3,ss)') dndlnk
     select case (itrans)
     case (:-1)  ! Tabulated P(k)
        write(splinefile,'(a,a)') trim(pkinfile),'.spline'
     case (3) ! Eisenstein & Hu
        write(splinefile,'(a25,f4.2,a1,a,a1,f4.2,a8,i1,a1,f4.2,a1,f6.4)') 'Data/Power_Spec/sigmacdm_',nspec,'_'&
             &,trim(sform),'_',kref ,'.spline.',itrans,'_',omega0,'_',omegab
     case (10) ! WDM
        write(splinefile,'(a25,f4.2,a1,a,a1,f4.2,a8,i1,a1,f4.2)') 'Data/Power_Spec/sigmacdm_',nspec,'_',trim(sform)&
             &,'_',kref,'.spline.',itrans,'_',mwdm
     case default ! CDM
        write(splinefile,'(a25,f4.2,a1,a,a1,f4.2,a8,i1)') 'Data/Power_Spec/sigmacdm_',nspec,'_',trim(sform),'_',kref&
             &,'.spline.',itrans
     end select
     io=0
     open (10,file=splinefile,status='unknown') 
     read (10,*,iostat=io) NSPLINE
     if (io.ne.0) then
        close (10)
#ifdef INFO
        write (0,*) 'spline_interp(): INFO - spline file not present creating it with make_spline()'
#endif
        call make_spline
        io=0
        open (10,file=splinefile,status='unknown') 
        read (10,*,iostat=io) NSPLINE
     end if
     if (NSPLINE.ne.NSPL) stop 'spline_interp(): FATAL - mismatch of NSPL( != NSPLINE)'
     do i=1,NSPL
        read (10,*,iostat=io) m(i),s(i),s2(i),a(i),a2(i)
     end do
     if (io.ne.0) then
        write (0,*) 'spline_interp(): FATAL - error reading ',trim(splinefile)
        write (0,*) '                         Run make_sigma_spline() subroutine with nspec = ',nspec
        write (0,*) '                         to create the required file.'
        stop
     end if
     close (10)        
     ifirst=0
  end if
  !
  imod=1+mod(imod,NMOD)
  klo=kplo(imod)  ! Look at position NMOD calls ago
  khi=kphi(imod)
  if (ms.lt.m(1)) then
     write (0,*) 'spline_interp(): FATAL - mass out of range.'
     write (0,*) '                 sigmacdm() was called with mass = ',ms/sclm     
     write (0,*) '                 which after scaling according the to adopted Gamma'
     write (0,*) '                 and Omega_0 corresponds to a look up mass of'
     write (0,*) '                 ms = ',ms ,' compared to the lower limit of ',m(1)
     stop
  else if (ms.gt.m(NSPL)) then
     write (0,*) 'spline_interp(): FATAL - mass out of range.'
     write (0,*) '                 sigmacdm() was called with mass = ',ms/sclm     
     write (0,*) '                 which after scaling according the to adopted Gamma'
     write (0,*) '                 and Omega_0 corresponds to a look up mass of'
     write (0,*) '                 ms = ',ms ,' compared to the upper limit of ',m(NSPL)
     stop
  else 
     if (m(khi).lt.ms .or. m(klo).gt.ms ) then ! Short cut if ms close to last
        klo=1                                  ! call, otherwise do binary search
        khi=NSPL
        do while (khi-klo.gt.1) 
           k=(khi+klo)/2
           if(m(k).gt.ms)then
              khi=k
           else
              klo=k
           endif
        end do
        h=m(khi)-m(klo)        ! Compute h factors.
        h2=(h**2)*0.1666667 
        invh=1.0/h
        !
        hp(imod)=h             ! Store for later call.
        invhp(imod)=invh
        hp2(imod)=h2
        kplo(imod)=klo
        kphi(imod)=khi
     else
        h=hp(imod)             ! Look up corresponding h factors.
        invh=invhp(imod)
        h2=hp2(imod)
     end if
     !            
     aa=(m(khi)-ms)*invh
     bb=(ms-m(klo))*invh
     a3=(aa**3-aa)
     b3=(bb**3-bb)
     sigma=aa*s(klo)+bb*s(khi)+(a3*s2(klo)+b3*s2(khi))*h2 ! Linear interpolation + spline.
     alpha=aa*a(klo)+bb*a(khi)+(a3*a2(klo)+b3*a2(khi))*h2 ! Linear interpolation + spline.
  end if
  return
end subroutine spline_interp

!
!  Fit a spline to sigma versus M.
!
!
!  Compute sigma(M) and the logarithmic slope alpha(m)
!  both from the old trusty fit and directly by integrating P(k).
!
!  This version sets Omega_0=h=Gamma=sigma_8=1. sigma(m) for
!  other values of Gamma and sigma_8 can then be computed from
!  this fit by two simple scaling of the mass and amplitude.
!  See the implementation in subroutine sigmacdm_spline.f .
!
!  The numerical integration is more accurate than the trusty fit
!  which is good to a few percent for M>10^10 Msol.
! 
subroutine make_spline
  !
  ! Uses
  use Cosmological_Parameters
  use Numerical_Parameters
  use Power_Spectrum_Parameters
  implicit none
  integer NT,i,ik
  parameter (NT=100000)
  real x(NSPL),y(NSPL),y2(NSPL),yp1,ypn,logm,a(NSPL),a2(NSPL),m,sigma,lnk,lnkmin,lnkmax,dlnk,sum,pw2k3,rf,Gamma_eff,pwdwk3,alph&
       &,suma
  character command*1024
  real pk,pw2k3_kmin,pwdwk3_kmin,pw2k3_kmax,pwdwk3_kmax
!  integer system_call,irem
  !
  ! First lock the power_spectrum files
  !
  open (10,file='Data/Power_Spec/lock',form='formatted',status='unknown')
  write (10,*) 'locked'
  close (10)
  !
  nspec=0.01*nint(100.0*nspec) ! Round to 2 decimal places.
#ifdef INFO
  if (dndlnk.eq.0.0) then
     write(0,'(a,f4.2)') 'make_spline(): INFO - adopting primordial spectal index nspec = ',nspec
  else
     write (0,'(a,f4.2,a,f6.3,a,f4.2,a)') 'make_spline(): INFO - adopting primordial spectral index nspec = ',nspec,'+0.5*'&
          &,dndlnk,'*ln(k/',kref,'h Mpc^-1)'
  endif
  write (0,*) '                N.B. To cover a mass range M_1 to M_2 the actual'
  write (0,*) '                tabulated values must cover M_1*SCLM to M_2*SCLM'
  write (0,*) '                where SCLM= Gamma**3/Omega_0'
#endif
  !
  Gamma_eff=1.0
  !
  ! We need to tabulate to sufficiently low masses to avoid galform
  ! crashing. AJB's method for estimating the minimum mass to tabulate
  ! is as follows:
  !     
  !  1) Take the lowest mres (merger tree mass resolution) that you
  !     need
  !
  !  2) Multiply by 0.01 (this should be the lowest mass ever required,
  !     and occurs in the a_nfw() routine, due to the presence of the
  !     factor F=0.01).
  !
  !  3) Multiply by the lowest value of Gamma**3/Omega (the scaling
  !     used in sigmacdm_spline). This is the lowest mass you need to
  !     tabulate.
  !
  do i=1,NSPL
     logm=0.0+18.0*real(i-1)/real(NSPL-1)
     m=10.0**logm
     x(i)=m
     ! Integrate 4 pi k^3 P(k) W^2(k) and 4 pi k^4 P(k) W(k) dW/du(kr)  
     ! calc top hat filter radius corresponding to mass M.
     select case (itrans)
     case (:-1) ! Tabulated P(k).
        ! Make tabulation with true Omega0.
        rf=(3.0*m/(4.0*PI*RHOCRIT*omega0))**(1.0/3.0)
     case (1:) ! Analytic CDM or WDM
        ! Make tabulation for Omega0=Gamma=1.
        rf=(3.0*m/(4.0*PI*RHOCRIT))**(1.0/3.0)
     case default
        stop 'make_spline(): FATAL - this subroutine does not handle itrans=0, which is power-law P(k)'
     end select
     lnkmax=5.0-log(rf)
     lnkmin=-9.0-log(rf)
     dlnk=(lnkmax-lnkmin)/float(NT-1)
     sum=0.0
     suma=0.0
     do ik=1,NT
        lnk=lnkmin+dlnk*float(ik-1)
        call pkfacs(exp(lnk),rf,Gamma_eff,pk,pw2k3,pwdwk3)
        sum=sum+pw2k3
        suma=suma+pwdwk3
     end do
     lnkmin=lnkmin-dlnk
     call pkfacs(exp(lnkmin),rf,Gamma_eff,pk,pw2k3_kmin,pwdwk3_kmin)
     lnkmax=lnkmax+dlnk
     call pkfacs(exp(lnkmax),rf,Gamma_eff,pk,pw2k3_kmax,pwdwk3_kmax)
     sigma=(sum+0.5*pw2k3_kmin+0.5*pw2k3_kmax)*dlnk
     alph=(suma+0.5*pwdwk3_kmin+0.5*pwdwk3_kmax)*dlnk
     alph=alph/(3.0*sigma)
     sigma=sqrt(sigma)
     y(i)=sigma ! Adopt values from
     a(i)=alph  ! numerical integration.
  end do
  yp1=2.0e+30   ! Indicates natural spline with zero second
  ypn=2.0e+30   ! Derivative at the end points.
  call spline(x,y,NSPL,yp1,ypn,y2)
  call spline(x,a,NSPL,yp1,ypn,a2)
  !
  ! Write spline coefficients to a file that later is read by
  ! the subroutine in sigmacdm_spline.f
  !
  open (10,file=splinefile,status='unknown')
  write (10,*) NSPL ! Nspline
#ifdef INFO
  write (0,*) 'M s=sigma(M) d^2s/dm^2  a=dlns/dlnm  d^2a/dm^2'
#endif
  do i=1,NSPL
#ifdef INFO
     write (0,*) x(i),y(i),y2(i),a(i),a2(i)
#endif
     write (10,*) x(i),y(i),y2(i),a(i),a2(i)
  end do
  close (10)
  !
  ! Unlock the power spectra files
  !    
  command='rm -f Data/Power_Spec/lock'
!  irem=system_call(command)
!#ifdef WARN
!  if (irem.ne.0) write (0,*) 'make_spline(): WARNING - failed to remove lock file'
!#endif
  return
end subroutine make_spline

!     
! This routine calculates:
!     
! pk =  P(k), and also
! pw2k3 = k^3 W(u)^2 P(k)
! pwdwk3 = k^3 W(u) u.dW/du
! where u=kr
!     
! for the case of an analytic CDM or WDM P(k) (itrans>0) or for a 
! tabulated P(k) (itrans<0)
!     
! itrans>0: P(k) propto k^nspec T(k)^2, analytic T(k)
! itrans<0: P(k) interpolated directly from table
!     
subroutine pkfacs(k,rf,Gamma_eff,pk,pw2k3,pwdwk3)
  !
  ! Uses
  use Cosmological_Parameters
  use Numerical_Parameters
  use Power_Spectrum_Parameters
  implicit none
  !
  ! Integers
  integer ifirst
  !
  ! Floats
  real dwin,Gamma_eff,k,lnk,lnpk,neff,pk,pw2k3,pwdwk3,q,rf,trans,transfer_function,u,win
  !     
  ! Saves
  save ifirst
  !
  ! Data
  data ifirst /1/
  !
  ! Code
  !     
  ! Calculate P(k)
  select case (itrans)
  case (:-1) ! Tabulated P(k).
     lnk=log(k)
     if (lnk.lt.lnktab(1)) then
        ! Extrapolate assuming a power-law P(k).
#ifdef INFO
        write (0,*) 'pkfacs(): INFO - k is below tabulated range: extrapolating'
#endif
        pk=exp(lnpktab(1)+(lnk-lnktab(1))*(lnpktab(2)-lnpktab(1))/(lnktab(2)-lnktab(1)))
     else if (lnk.gt.lnktab(nktab)) then
        ! Extrapolate assuming a power-law P(k)
#ifdef INFO
        write (0,*) 'pkfacs(): INFO - k is above tabulated range: extrapolating'
#endif
        pk=exp(lnpktab(nktab)+(lnk-lnktab(nktab))*(lnpktab(nktab)-lnpktab(nktab-1))/(lnktab(nktab)-lnktab(nktab-1)))
     else
        ! Calc P(k) by linear interpolation in ln(P(k)) vs ln(k)
        call interp(nktab,lnktab,lnpktab,lnk,lnpk)
        pk=exp(lnpk)
     endif
  case (1:) ! Analytic transfer function for CDM or WDM.
     ! Calculate transfer function T(k).
     q=k/Gamma_eff
     trans=transfer_function(k,q,Gamma_eff)
     ifirst=0
     ! Multiply by primordial P(k) propto k^nspec to get final P(k)
     neff=nspec+0.5*dndlnk*log(k/kref)
     pk=(trans**2)*(k/KHORIZON)**neff ! P(k)
  case (0)
     stop 'pkfacs(): FATAL - this function does not handle itrans=0, which is power-law P(k)'
  end select
  !     
  ! Calc k-space window function for top hat in real space, 
  ! and its derivative.
  u=k*rf
  ! W(u)
  win=3.0*(sin(u)/u-cos(u))/u**2
  ! u*dW(u)/du
  dwin=3.0*(3.0*cos(u)-3.0*sin(u)/u+u*sin(u))/u**2
  ! Multiply by P(k)
  pw2k3=(k**3)*(win**2)*pk
  pwdwk3=(k**3)*win*dwin*pk
  return
end subroutine pkfacs
