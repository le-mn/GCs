!% Contains a module which implements calculations of the power spectrum and linear theory mass variance, $\sigma(M)$.

module Power_Spectrum
  !% Contains subroutines which implement calculations of the power spectrum and linear theory mass variance, $\sigma(M)$.
  private
  public :: sigma_Mass,Power_Spectrum_Read_Parameters,sigma_Mass_Spline_Interp, sigma_Mass_no_alpha

  integer                 :: itrans
  real                    :: gamma,sigma8
  character               :: splinefile*1024
  real                    :: mwdm,dndlnk,kref,nspec,scla,sclm
  real                    :: M_min,M_max
  real,      parameter    :: M_min_default=1.0,M_max_default=1.0e16 ! Default minimum and maximum masses for sigma(M) table.

  ! Arrays for tabulated sigma(m).
  real,      pointer      :: a_sigma_Mass(:),a2_sigma_Mass(:),m_sigma_Mass(:),s2_sigma_Mass(:),s_sigma_Mass(:)
  integer,   dimension(2) :: a_sigma_Mass_idx,a2_sigma_Mass_idx,m_sigma_Mass_idx,s2_sigma_Mass_idx,s_sigma_Mass_idx
  integer                 :: NSPLINE

  ! Arrays for tabulated power-spectrum.
  real,      pointer      :: lnktab(:),lnpktab(:)
  integer                 :: nktab
  character               :: pkinfile*1024

contains
  subroutine Power_Spectrum_Read_Parameters
    !% Reads the parameters which describe the power spectrum.
    use Read_Parameters
    use Report_Error_Module
    use Cosmology
    ! Value of itrans controls form used for input P(k):
    !     itrans < 0 : read P(k) from file;
    !     itrans = 0 : power-law P(k) propto k^nspec;
    !  1<=itrans<= 9 : analytical CDM transfer function, uses Gamma, nspec;
    !     itrans =10 : analytical WDM transfer function, uses Mwdm, Gamma, nspec.
    implicit none
    call get_parameter("itrans", itrans,   range = (/ -1, 10 /) )
    call get_parameter("nspec" , nspec,    required = (itrans.ge.0))
    call get_parameter("dndlnk", dndlnk,   required = (itrans.ge.0))
    call get_parameter("kref"  , kref,     required = (itrans.ge.0))
    call get_parameter("gamma" , gamma,    required = (itrans.ge.0))
    call get_parameter("PKfile", pkinfile, required = (itrans.eq.-1) ) ! Tabulated power spectrum file.
    call get_parameter("Mwdm"  ,mwdm,      required = (itrans.eq.10) ) ! WDM particle mass.

    select case (itrans)
    case (:-1)     ! Read P(k) from file.
#ifdef INFO
       write (0,*) 'Power_Spectrum_Read_Parameters(): INFO - reading P(k) from file'
#endif
       nspec=0.0   ! nspec & gamma not used in this case.
       gamma=0.0
    case (0)       ! Power-law P(k).
#ifdef INFO
       write (0,*) 'Power_Spectrum_Read_Parameters(): INFO - using power-law P(k) with index nspec=',nspec
#endif
       gamma=0.0   ! gamma not used.
    case (1:3,5:9) ! Analytical CDM P(k).
#ifdef INFO
       write (0,*) 'Power_Spectrum_Read_Parameters(): INFO - using analytical CDM transfer function'
#endif
       if (gamma.le.0.0) call Report_Fatal_Error('Power_Spectrum_Read_Parameters(): FATAL - need Gamma>0')
       if (itrans.eq.3.and.gamma.ne.omega0*h0) then
          write (0,*) 'Power_Spectrum_Read_Parameters(): FATAL - Eisenstein & Hu T(k) selected, but gamma!=omega0*h0'
          write (0,*) '                                      gamma=',gamma
          write (0,*) '                                  omega0*h0=',omega0*h0
          call Report_Fatal_Error
       end if
    case (10)      ! Analytical WDM P(k).
       write (0,*) 'Power_Spectrum_Read_Parameters(): INFO - using analytical WDM transfer function'
       if (gamma.le.0.0) call Report_Fatal_Error('parameters(): FATAL - need Gamma>0')
    case default
       write (0,*) 'Power_Spectrum_Read_Parameters() FATAL - unrecognized value itrans=',itrans
       call Report_Fatal_Error
    end select

    call get_parameter("sigma8",sigma8, range = (/0.0,10.0/)) ! sigma8.
    return
  end subroutine Power_Spectrum_Read_Parameters

  subroutine sigma_Mass_Deallocate
    !% Deallocates arrays which store the tabulated $\sigma(M)$.
    use Memory_Usage_Module
    implicit none
    if (associated(m_sigma_Mass))  call Dealloc_Array(m_sigma_Mass ,m_sigma_Mass_idx )
    if (associated(s_sigma_Mass))  call Dealloc_Array(s_sigma_Mass ,s_sigma_Mass_idx )
    if (associated(s2_sigma_Mass)) call Dealloc_Array(s2_sigma_Mass,s2_sigma_Mass_idx)
    if (associated(a_sigma_Mass))  call Dealloc_Array(a_sigma_Mass ,a_sigma_Mass_idx )
    if (associated(a2_sigma_Mass)) call Dealloc_Array(a2_sigma_Mass,a2_sigma_Mass_idx)
    return
  end subroutine Sigma_Mass_Deallocate

  real function sigma_Mass_no_alpha(m)
    !% Interface to \hyperlink{func:sigma_mass}{{\tt sigma\_Mass()}} when no $\alpha$ is required.
    implicit none
    real, intent(in) :: m
    real             :: alpha
    sigma_Mass_no_alpha=sigma_Mass(m,alpha)
    return
  end function sigma_Mass_no_alpha

  real function sigma_Mass(m,alpha)
    !% This routine computes $\sigma(M)$. It returns both $\sigma$ and $\alpha = \d\ln(\sigma)/\d\ln M$.
    !%
    !% It handles various cases, depending on value of {\tt itrans}:
    !%
    !% {\tt itrans}$<0$: uses tabulated $P(k)$ read in from file
    !% {\tt itrans}$=0$: uses power-law $P(k) \propto k^{\tt nspec}$
    !% $1\le$ {\tt itrans} $\le 9$: uses analytical CDM power spectrum
    !% \begin{equation}
    !% P(k) \propto k^{n_{\rm s}} T^2(k).
    !% \end{equation}
    !% {\tt itrans}$=10$: uses analytical WDM power spectrum
    !% \begin{equation}
    !% P(k) \propto k^{n_{\rm s}} T^2(k).
    !% \end{equation}
    !%
    !% For the analytical CDM and WDM cases ({\tt itrans}$>0$) it uses the variable $\Gamma$.
    !%
    !% In all cases:
    !%
    !% The input mass {\tt m} is assumed to be in units $M_\odot/h$.
    !% $\sigma(M)$ is normalized to the $\sigma_8$ specified in the input parameter file
    !% (overriding any normalization, e.g. in an input file).
    !%
    !% On the first call, it calculates coefficients for a spline fit to
    !% $\sigma(M)$, and also scaling factors {\tt sclm} and {\tt scla} for the mass and
    !% normalization.
    use Cosmology
    use Numerical_Parameters
    use File_Units_Module
    use Memory_Usage_Module
    implicit none
    real,    intent(in)  :: m
    real,    intent(out) :: alpha
    integer              :: fileend,i,iunit
    integer, save        :: ifirst=1
    real                 :: ktab,ms,m8,pktab,sigma

    ! If first call compute the appropriate scaling factors sclm and
    ! scla, calculate sigma in an 8 Mpc/h sphere for the arbitrarily
    ! normalized input P(k), and hence get the factor scla by which the
    ! un-normalized sigma(M) must be multiplied to get the desired
    ! sigma8 normalization for CDM or WDM with an analytical transfer
    ! function, also calculate the scaling factor in mass sclm that
    ! relates the true P(k) to that for Omega0=1, Gamma=1.
    if (ifirst.eq.1) then
       m8=M8CRIT*omega0                        ! The mass within an 8Mpc/h sphere.
       transfer_function: select case (itrans) ! Branch on transfer function.
       case (:-1)                              ! Read P(k) from input file.
#ifdef INFO
          write (0,*) 'sigma_Mass(): INFO - reading P(k) from file ',trim(pkinfile)
#endif
          nktab=Count_Lines_in_File(pkinfile)-1    ! Count number of lines in power spectrum file.
#ifdef INFO
          write (0,*) 'sigma_Mass(): INFO - number of wavenumbers read in =',nktab
#endif
          call Alloc_Array(lnktab ,nktab,'lnktab' )
          call Alloc_Array(lnpktab,nktab,'lnpktab')
          iunit=Get_Unit()                         ! Get a file unit to read input file from.
          open (iunit,file=pkinfile,status='old')
          read (iunit,*)                           ! Skip 1 line header.
          do i=1,nktab
             read (iunit,*,iostat=fileend) ktab,pktab ! File contains k in units h/Mpc & P(k) in units (Mpc/h)^3.
             lnktab(i)=log(ktab)                      ! Store k & P(k) as logs for interpolation.
             lnpktab(i)=log(pktab)
          enddo
          close(iunit)
          sclm=1.0                                    ! Mass scaling not used in this case.
          call sigma_Mass_Spline_Interp(m8,sigma,alpha) ! Calculate sigma8 for input spectrum.
       case (0)                                    ! Power-law P(k).
          m8=M8CRIT*omega0                         ! The mass within an 8Mpc/h sphere.
          sclm=1.0/m8
          ms=m8*sclm
          sigma=ms**(-(nspec+3.0)/6.0)
       case default                                ! Analytic CDM or WDM P(k).
          sclm=gamma**3/omega0                     ! Compute the required scaling factors sclm and scla.
          m8=M8CRIT*omega0                         ! The mass within an 8Mpc/h sphere.
          ms=m8*sclm
          call sigma_Mass_Spline_Interp(ms,sigma,alpha)   ! The spline fit to CDM for Gamma=1.
       end select transfer_function
       scla=sigma8/sigma                                  ! Scales sigma_8 to required value.
       ifirst=0                                           ! Indicates first call complete and sclm and scla are set.
    end if
    ms=m*sclm
    select case (itrans)
    case (0)                                              ! Power-law P(k).
       alpha=-(nspec+3.0)/6.0
       sigma_Mass=scla*(ms**alpha)
    case default                                          ! CDM or WDM or tabulated.
       call sigma_Mass_Spline_Interp(ms,sigma,alpha)      ! Use spline fit.
       sigma_Mass=sigma*scla
    end select
    return
  end function sigma_Mass

  subroutine sigma_Mass_Spline_Interp(ms,sigma,alpha)
    !% The spline interpolation routine called by \hyperlink{func:sigma_mass}{{\tt sigma\_Mass()}}.
    !%
    !% The first time it is called, it checks to see if a file containing the spline coefficients for the current $P(k)$ already
    !% exists. If it does, it reads this file, otherwise it generates the spline coefficients afresh.
    !%
    !% Saves last {\tt NMOD} positions in look-up table along with corresponding values of {\tt h}, {\tt h2} and {\tt inv}. This is optimal if every
    !% {\tt NMOD} calls the value of {\tt ms} is approximately repeated. {\tt NMOD}$=2$ is optimal when used in conjunction with \hyperlink{func:split}{{\tt split()}}
    !% as most calls don't produce a split and so $M$ changes slowly and one keeps evaluating $\sigma(M)$ and $\sigma(M/2)$.
    use Band_Memory
    use Cosmology
    use Numerical_Parameters
    use File_Units_Module
    use Memory_Usage_Module
    use Report_Error_Module
    implicit none

    real,    intent(in)            :: ms
    real,    intent(out)           :: sigma,alpha
    integer, parameter             :: NMOD=2
    integer, save                  :: ifirst=1, imod, kphi(NMOD), kplo(NMOD)
    real,    save, dimension(NMOD) :: invhp(NMOD),hp2(NMOD),hp(NMOD)
    integer                        :: i,io,k,khi,klo,iunit
    real                           :: a3,aa,b3,bb,h,h2,invh,M_min_in_file,M_max_in_file
    character                      :: sform*10
    logical                        :: file_exists,loadfile,remake

    ! Determine if we need to remake the file
    loadfile=.false.
    remake=.false.
    if (ifirst.eq.1) then
       loadfile=.true.                  ! Always load the file on the first call.
       M_min=min(M_min_default,0.1*ms)  ! Specify the minimum...
       M_max=max(M_max_default,10.0*ms) ! ...and maximum masses required.

       write (0,*) 'First call - will load spline file'

    else
       if (ms.lt.M_min.or.ms.gt.M_max) then
          loadfile=.true.                  ! Remake the file if ms is out of range.
          remake=.true.
          M_min=min(M_min_default,min(0.1*ms,M_min))  ! Specify the minimum...
          M_max=max(M_max_default,max(10.0*ms,M_max)) ! ...and maximum masses required.

          write (0,*) 'Mass out of range - will remake spline file ',M_min,M_max

       end if
    end if

    ! Remake and reload the file if necessary.
    if (loadfile) then
       nspec=0.01*nint(100.0*nspec) ! Round nspec to two decimal places.
       write (sform,'(sp,f6.3,ss)') dndlnk
       select case (itrans)         ! Construct a file name for the power spectrum file.
       case (:-1)                   ! Tabulated P(k).
          write(splinefile,'(a,a)') trim(pkinfile),'.spline'
       case (3)                     ! Eisenstein & Hu
          write(splinefile,'(a25,f4.2,a1,a,a1,f4.2,a8,i1,a1,f4.2,a1,f6.4,a1,f6.4)') 'Data/Power_Spec/sigmacdm_',nspec,'_' &
               &,trim(sform),'_',kref ,'.spline.',itrans,'_',omega0,'_',omegab,'_',CMB_T0
       case (10)                    ! WDM.
          write(splinefile,'(a25,f4.2,a1,a,a1,f4.2,a8,i1,a1,f4.2)') 'Data/Power_Spec/sigmacdm_',nspec,'_',trim(sform)&
               &,'_',kref,'.spline.',itrans,'_',mwdm
       case default                 ! CDM.
          write(splinefile,'(a25,f4.2,a1,a,a1,f4.2,a8,i1)') 'Data/Power_Spec/sigmacdm_',nspec,'_',trim(sform),'_',kref&
               &,'.spline.',itrans
       end select
       if (.not.remake) then
          inquire(file=splinefile,exist=file_exists)
          if (.not.file_exists) then
#ifdef INFO
             write (0,*) 'spline_interp(): INFO - spline file not present creating it with sigma_Mass_Make_Spline()'
#endif
             remake=.true.

             write (0,*) 'File does not exist - will remake'

          else
             iunit=Get_Unit()
             open (iunit,file=splinefile,status='unknown',iostat=io)
             read (iunit,*,iostat=io) NSPLINE,M_min_in_file,M_max_in_file ! Find number of tabulated points and tabulated mass
             ! range.
             close (iunit)
             if (M_min_in_file.gt.M_min.or.M_max_in_file.lt.M_max) then
                remake=.true.

                write (0,*) 'File does not have sufficient mass range - will remake'

             end if
          end if
       end if
       if (remake) call sigma_Mass_Make_Spline
       iunit=Get_Unit()
       open (iunit,file=splinefile,status='unknown',iostat=io)
       read (iunit,*,iostat=io) NSPLINE,M_min,M_max   ! Find number of tabulated points and tabulated mass range.
       call sigma_Mass_Deallocate                     ! Deallocate old arrays.
       call Alloc_Array(m_sigma_Mass ,NSPLINE,'m_sigma_Mass' ,idx=m_sigma_Mass_idx )
       call Alloc_Array(s_sigma_Mass ,NSPLINE,'s_sigma_Mass' ,idx=s_sigma_Mass_idx )
       call Alloc_Array(s2_sigma_Mass,NSPLINE,'s2_sigma_Mass',idx=s2_sigma_Mass_idx)
       call Alloc_Array(a_sigma_Mass ,NSPLINE,'a_sigma_Mass' ,idx=a_sigma_Mass_idx )
       call Alloc_Array(a2_sigma_Mass,NSPLINE,'a2_sigma_Mass',idx=a2_sigma_Mass_idx)
       do i=1,NSPLINE                                ! Read in tabulated sigma(m).
          read (iunit,*,iostat=io) m_sigma_Mass(i),s_sigma_Mass(i),s2_sigma_Mass(i),a_sigma_Mass(i),a2_sigma_Mass(i)
          ! Check for conditions on alpha that are assumed by merger tree code.
          if (i.gt.1) then
             if (a_sigma_Mass(i).gt.0.0.and.m_sigma_Mass(i).gt.1.0e4*m_Sigma_Mass(1)) then
                write (0,*) 'sigma_Mass_Spline_Interp(): WARNING - alpha > 0'
                write (0,*) '                            conflicts with assumptions used in building merger trees'
             end if
             if (a_Sigma_Mass(i).ge.a_Sigma_Mass(i-1).and.m_sigma_Mass(i).gt.1.0e4*m_Sigma_Mass(1)) then
                write (0,*) 'sigma_Mass_Spline_Interp(): WARNING - d(alpha)/dm >= 0'
                write (0,*) '                            conflicts with assumptions used in building merger trees'
             end if
          end if
       end do

       if (io.ne.0) then
          write (0,*) 'spline_interp(): FATAL - error reading ',trim(splinefile)
          write (0,*) '                         Run sigma_Mass_Make_Spline() subroutine with nspec = ',nspec
          write (0,*) '                         to create the required file.'
          call Report_Fatal_Error
       end if
       close (iunit)
       ifirst=0
       kplo=NSPLINE
       kphi=1
       imod=0
    end if

    imod=1+mod(imod,NMOD)
    klo=kplo(imod)                                   ! Look at position NMOD calls ago.
    khi=kphi(imod)
    if (ms.lt.m_sigma_Mass(1)) then
       write (0,*) 'sigma_Mass_Spline_Interp(): FATAL - mass out of range.'
       write (0,*) '                            sigma_Mass() was called with mass = ',ms/sclm
       write (0,*) '                            which after scaling according the to adopted Gamma'
       write (0,*) '                            and Omega_0 corresponds to a look up mass of'
       write (0,*) '                            ms = ',ms ,' compared to the lower limit of ',m_sigma_Mass(1)
       write (0,*) '                            This should not happen!'
       call Report_Fatal_Error
    else if (ms.gt.m_sigma_Mass(NSPLINE)) then
       write (0,*) 'sigma_Mass_Spline_Interp(): FATAL - mass out of range.'
       write (0,*) '                            sigma_Mass() was called with mass = ',ms/sclm
       write (0,*) '                            which after scaling according the to adopted Gamma'
       write (0,*) '                            and Omega_0 corresponds to a look up mass of'
       write (0,*) '                            ms = ',ms ,' compared to the upper limit of ',m_sigma_Mass(NSPLINE)
       write (0,*) '                            This should not happen!'
       call Report_Fatal_Error
    else
       if (m_sigma_Mass(khi).lt.ms.or.m_sigma_Mass(klo).gt.ms) then ! Short cut if ms close to last call, otherwise do binary
          ! search.
          klo=1
          khi=NSPLINE
          do while (khi-klo.gt.1)
             k=(khi+klo)/2
             if (m_sigma_Mass(k).gt.ms) then
                khi=k
             else
                klo=k
             endif
          end do
          h=m_sigma_Mass(khi)-m_sigma_Mass(klo)        ! Compute h factors.
          h2=(h**2)*0.1666667
          invh=1.0/h
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
       aa=(m_sigma_Mass(khi)-ms)*invh
       bb=(ms-m_sigma_Mass(klo))*invh
       a3=(aa**3-aa)
       b3=(bb**3-bb)
       sigma=aa*s_sigma_Mass(klo)+bb*s_sigma_Mass(khi)+(a3*s2_sigma_Mass(klo)+b3*s2_sigma_Mass(khi))*h2 ! Linear interpolation +
       ! spline.
       alpha=aa*a_sigma_Mass(klo)+bb*a_sigma_Mass(khi)+(a3*a2_sigma_Mass(klo)+b3*a2_sigma_Mass(khi))*h2 ! Linear interpolation +
       ! spline.
    end if
    return
  end subroutine sigma_Mass_Spline_Interp

  subroutine sigma_Mass_Make_Spline
    !% Fit a spline to $\sigma(M)$.
    !%
    !% Compute $\sigma(M)$ and the logarithmic slope $\alpha(M)$ directly by integrating $P(k)$.
    !
    !% This version sets $\Omega_0=h=\Gamma=\sigma_8=1$. $\sigma(M)$ for other values of $\Gamma$ and $\sigma_8$ can then be computed from
    !% this fit by two simple scaling of the mass and amplitude. See the implementation in \hyperlink{sub:sigma_mass_spline_interp}{{\tt sigma\_Mass\_Spline\_Interp()}}.
    use Cosmology
    use Numerical_Parameters
    use System_Interaction_Module
    use File_Units_Module
    use Spline_Interpolation_Module
    use Report_Error_Module
    implicit none
    integer              :: i,ik
    integer, parameter   :: NT=1e5, NSPL_MIN=200
    real,    parameter   :: NSPL_PER_DECADE=10.0
    real,    allocatable :: x(:),y(:),y2(:),a(:),a2(:)
    real                 :: yp1,ypn,logm,m,sigma,lnk,lnkmin,lnkmax,dlnk,sum,pw2k3,rf,Gamma_eff,pwdwk3,alph,suma
    character            :: command*1024
    real                 :: pk,pw2k3_kmin,pwdwk3_kmin,pw2k3_kmax,pwdwk3_kmax,lg_M_min,lg_M_max
    integer              :: irem,iunit,alloc_err,NSPL

    ! Determine how many points we want in the spline fit. (Choose to have at least 10 per decade of mass.)
    lg_M_min=log10(M_min)
    lg_M_max=log10(M_max)
    NSPL=max(int((lg_M_max-lg_M_min)*NSPL_PER_DECADE),NSPL_MIN)

    ! Allocate workspace arrays.
    allocate(x(NSPL),stat=alloc_err)
    if (alloc_err.ne.0) call Report_Fatal_Error('sigma_Mass_Make_Spline(): FATAL - failed to allocate memory [x]')
    allocate(y(NSPL),stat=alloc_err)
    if (alloc_err.ne.0) call Report_Fatal_Error('sigma_Mass_Make_Spline(): FATAL - failed to allocate memory [y]')
    allocate(y2(NSPL),stat=alloc_err)
    if (alloc_err.ne.0) call Report_Fatal_Error('sigma_Mass_Make_Spline(): FATAL - failed to allocate memory [y2]')
    allocate(a(NSPL),stat=alloc_err)
    if (alloc_err.ne.0) call Report_Fatal_Error('sigma_Mass_Make_Spline(): FATAL - failed to allocate memory [a]')
    allocate(a2(NSPL),stat=alloc_err)
    if (alloc_err.ne.0) call Report_Fatal_Error('sigma_Mass_Make_Spline(): FATAL - failed to allocate memory [a2]')

    ! First lock the power_spectrum files. Other executables which use the sigma(m) files should not run while this lock file exists.
    iunit=Get_Unit()
    open (iunit,file='Data/Power_Spec/lock',form='formatted',status='unknown')
    write (iunit,*) 'locked'
    close (iunit)

    nspec=0.01*nint(100.0*nspec) ! Round to 2 decimal places.
    Gamma_eff=1.0
#ifdef INFO
    if (dndlnk.eq.0.0) then
       write(0,'(a,f4.2)') 'sigma_Mass_Make_Spline(): INFO - adopting primordial spectal index nspec = ',nspec
    else
       write (0,'(a,f4.2,a,f6.3,a,f4.2,a)') 'sigma_Mass_Make_Spline(): INFO - adopting primordial spectral index nspec = ',nspec&
            &,'+0.5*' ,dndlnk,'*ln(k/',kref,'h Mpc^-1)'
    endif
    write (0,*) '                          N.B. To cover a mass range M_1 to M_2 the actual'
    write (0,*) '                          tabulated values must cover M_1*SCLM to M_2*SCLM'
    write (0,*) '                          where SCLM= Gamma**3/Omega_0.'
#endif

    do i=1,NSPL
       logm=lg_M_min+(lg_M_max-lg_M_min)*real(i-1)/real(NSPL-1)
       m=10.0**logm
       x(i)=m
       ! Integrate 4 pi k^3 P(k) W^2(k) and 4 pi k^4 P(k) W(k) dW/du(kr). Calculate top hat filter radius corresponding to mass M.
       select case (itrans)
       case (:-1)                                        ! Tabulated P(k).
          rf=(3.0*m/(4.0*PI*RHOCRIT*omega0))**(1.0/3.0)  ! Make tabulation with true Omega0.
       case (1:)                                         ! Analytic CDM or WDM
          rf=(3.0*m/(4.0*PI*RHOCRIT))**(1.0/3.0)         ! Make tabulation for Omega0=Gamma=1.
       case default
          call Report_Fatal_Error('sigma_Mass_Make_Spline(): FATAL - this subroutine does not handle itrans=0, which is a power&
               &-law P(k)')
       end select
       lnkmax=5.0-log(rf)
       lnkmin=-9.0-log(rf)
       dlnk=(lnkmax-lnkmin)/float(NT-1)
       sum=0.0
       suma=0.0
       do ik=1,NT
          lnk=lnkmin+dlnk*float(ik-1)
          call Pk_Factors(exp(lnk),rf,Gamma_eff,pk,pw2k3,pwdwk3)
          sum=sum+pw2k3
          suma=suma+pwdwk3
       end do
       lnkmin=lnkmin-dlnk
       call Pk_Factors(exp(lnkmin),rf,Gamma_eff,pk,pw2k3_kmin,pwdwk3_kmin)
       lnkmax=lnkmax+dlnk
       call Pk_Factors(exp(lnkmax),rf,Gamma_eff,pk,pw2k3_kmax,pwdwk3_kmax)
       sigma=(sum+0.5*pw2k3_kmin+0.5*pw2k3_kmax)*dlnk
       alph=(suma+0.5*pwdwk3_kmin+0.5*pwdwk3_kmax)*dlnk
       alph=alph/(3.0*sigma)
       sigma=sqrt(sigma)
       y(i)=sigma ! Adopt values from
       a(i)=alph  ! numerical integration.
    end do
    yp1=2.0e+30   ! Indicates natural spline with zero second
    ypn=2.0e+30   ! Derivative at the end points.
    call Spline_Interpolation_Make(x,y,yp1,ypn,y2)
    call Spline_Interpolation_Make(x,a,yp1,ypn,a2)

    ! Write spline coefficients to a file that later is read by sigma_Mass_Spline_Interp().
    iunit=Get_Unit()
    open (iunit,file=splinefile,status='unknown')
    write (iunit,*) NSPL,M_min,M_max ! Number of points in spline file.
#ifdef INFO
    write (0,*) 'sigma_Mass_Make_Spline(): INFO - M s=sigma(M) d^2s/dm^2  a=dlns/dlnm  d^2a/dm^2'
#endif
    do i=1,NSPL
#ifdef INFO
       write (0,*) '                                 ',x(i),y(i),y2(i),a(i),a2(i)
#endif
       write (iunit,*) x(i),y(i),y2(i),a(i),a2(i)
    end do
    close (iunit)

    ! Unlock the power spectra files
    command='rm -f Data/Power_Spec/lock'
    irem=system_call(command)
#ifdef WARN
    if (irem.ne.0) write (0,*) 'sigma_Mass_Make_Spline(): WARNING - failed to remove lock file'
#endif

    ! Deallocate workspace arrays.
    deallocate(x,y,y2,a,a2)
    return
  end subroutine sigma_Mass_Make_Spline

  subroutine Pk_Factors(k,rf,Gamma_eff,pk,pw2k3,pwdwk3)
    !% Calculates:
    !%
    !% {\tt pk} $=  P(k)$, and also
    !% {\tt pw2k3} $= k^3 W(u)^2 P(k)$
    !% {\tt pwdwk3} $= k^3 W(u) u \d W/\d u$,
    !% where $u=kr$, for the case of an analytic CDM or WDM P(k) ({\tt itrans}$>0$) or for a
    !% tabulated $P(k)$ ({\tt itrans}$<0$).
    !%
    !% {\tt itrans}$>0$: $P(k) \propto k^{n_{\rm s}} T(k)^2$, analytic $T(k)$;
    !% {\tt itrans}$<0$: $P(k)$ interpolated directly from table.
    use Linear_Interpolation_Module
    use Numerical_Parameters
    use Report_Error_Module
    implicit none
    real,    intent(in)  :: k,rf,Gamma_eff
    real,    intent(out) :: pk,pw2k3,pwdwk3
    real,    parameter   :: UEPS=0.0316228  ! Value of u below which to switch to series solution.
    real                 :: dwin,lnk,neff,trans,u,win,k3winpk,u2
    logical              :: extrapolated

    select case (itrans)                    ! Calculate P(k).
    case (:-1)                              ! Tabulated P(k).
       lnk=log(k)
       pk=exp(Lin_Interp(nktab,lnktab,lnpktab,lnk,extrapolated))
#ifdef INFO
       if (extrapolated) write (0,*) 'Pk_Factors(): INFO - k is above tabulated range: extrapolating'
#endif
    case (1:)                               ! Analytic transfer function for CDM or WDM.
       trans=transfer_function(k,Gamma_eff) ! Calculate transfer function T(k).
       neff=nspec+0.5*dndlnk*log(k/kref)    ! Multiply by primordial P(k) propto k^nspec to get final P(k).
       pk=(trans**2)*(k/KHORIZON)**neff     ! P(k).
    case (0)
       call Report_Fatal_Error('Pk_Factors(): FATAL - this function does not handle itrans=0, which is power-law P(k)')
    end select
    u=k*rf                                  ! Calc k-space window function for top hat in real space, and its derivative.
    if (u.gt.UEPS) then                     ! Use full solutions:
       win=3.0*(sin(u)/u-cos(u))/u**2       ! W(u).
       dwin=3.0*(-win+sin(u)/u)             ! u*dW(u)/du.
    else                                    ! Use series solutions for small u:
       u2=u**2
       win=1.0+u2*(-0.1+u2*(1.0/280.0-u2/15120.0)) ! W(u).
       dwin=u2*(-0.2+u2*(1.0/70.0-u2/2520.0))      ! u*dW(u)/du.       
    end if
    ! Multiply by P(k).
    k3winpk=(k**3)*win*pk
    pw2k3=win*k3winpk
    pwdwk3=dwin*k3winpk
    return
  end subroutine Pk_Factors

  real function Transfer_Function(k,Gamma_eff)
    !% Computes and returns the transfer function, $T(k)$, at wavenumber $k=${\tt k} for an effective value of the shape parameter $\Gamma$
    !% specified by {\tt Gamma\_eff}. The transfer function implemented depends on the value of the model parameter {\tt itrans}:
    !% \begin{itemize}
    !% \item {\tt itrans}$=1$ : \href{http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=1986ApJ...304...15B}{BBKS cold dark matter $T(k)$}.
    !% \item {\tt itrans}$=2$ : \href{http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=1984ApJ...285L..45B}{Bond \& Efstathiou (1987) cold dark matter $T(k)$.}
    !% \item {\tt itrans}$=3$ : \href{http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=1999ApJ...511....5E}{Eisenstein \& Hu (1999) cold dark matter $T(k)$.}
    !% \item {\tt itrans}$=10$ : \href{http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=1986ApJ...304...15B}{BBKS warm dark matter $T(k)$}.
    !% \end{itemize}
    use Cosmology
    use Report_Error_Module
    implicit none
    real, intent(in) :: k,Gamma_eff
    real             :: alphav,b1,b2,betac,Bk,C,fb,fc,fcb,fv,fvb,Gammaeff,kRfw,L,Nv,pc,pcb,qeff,qEH,qv,s,Theta27,Tsup,yd,zd,zeq,q

    q=k/Gamma_eff ! Compute the scaled wavenumber.
    select case (itrans)
    case (1)
       ! BBKS cold dark matter transfer function.
       transfer_function=transfer_function_BBKS_CDM(q)
    case (2)
       ! Bond & Efstathiou (1984; eqn. 6) cold dark matter transfer function. (This uses the value of the fitting parameters from
       ! line 1 of their Table 1, rescaled to h0=1.0).
       transfer_function=(1.0+(((6.4*q)+((3.0*q)**1.5)+((1.7*q)**2.0))**1.13))**(-1.0/1.13)
    case (3)
       ! Eisenstein & Hu (1999, ApJ, 511, 5) cold dark matter transfer function (eqn. 18 and associated equations).
       Theta27=CMB_T0/2.7 ! Present day CMB temperature [in units of 2.7K].
       zeq=2.50e4*omega0*(h0**2)/(Theta27**4) ! Redshift of matter-radiation equality.
       ! Compute redshift at which baryons are released from Compton drag of photons (eqn. 2)
       b1=0.313*((omega0*(h0**2))**(-0.419))*(1.0+0.607*((omega0*(h0**2))**0.674))
       b2=0.238*((omega0*(h0**2))**0.223)
       zd=1291.0*((omega0*(h0**2))**0.251)*(1.0+b1*((omegab*(h0**2))**b2))/(1.0+0.659*((omega0*(h0**2))**0.828))
       ! Relative expansion factor between previous two computed redshifts.
       yd=(1.0+zeq)/(1.0+zd)
       ! Compute the comoving distance that a sound wave can propagate prior to zd (i.e. sound horizon; eq. 4)
       s=44.5*log(9.83/omega0/(h0**2))/sqrt(1.0+10.0*((omegab*(h0**2))**0.75))
       ! Compute effective q.
       qEH=k*(Theta27**2)/Gamma_eff
       ! Specify properties of neutrinos.
       fv=0.0 ! No neutrinos.
       Nv=0.0
       ! Compute baryonic and cold dark matter fractions.
       fb=omegab/omega0
       fc=(omega0-omegab)/omega0
       ! Total matter fraction.
       fcb=fb+fc
       ! Baryonic + neutrino fraction.
       fvb=fv+fb
       ! Compute small scale suppression factor (eqn. 15).
       pc=0.25*(5.0-sqrt(1.0+24.0*fc))
       pcb=0.25*(5.0-sqrt(1.0+24.0*fcb))
       alphav=(fc/fcb)*((5.0-2.0*(pc+pcb))/(5.0-4.0*pcb))*((1.0-0.533*fvb+0.126*(fvb**3))*((1.0+yd)**(pcb-pc))/(1.0-0.193&
            &*sqrt(fv*Nv)+0.169*fv*(Nv**0.2)))*(1.0+0.5*(pc-pcb)*(1.0+1.0/(3.0-4.0*pc)/(7.0-4.0*pcb))/(1.0+yd))
       ! Compute rescaled shape parameter (eqn. 16)
       Gammaeff=Gamma_eff*(sqrt(alphav)+(1.0-sqrt(alphav))/(1.0+((0.43*(k*h0)*s)**4)))
       qeff=k*(Theta27**2)/Gammaeff
       betac=1.0/(1.0-0.949*fvb)                    ! Eqn. 21.
       L=log(exp(1.0)+1.84*betac*sqrt(alphav)*qeff) ! Eqn. 19.
       C=14.4+325.0/(1.0+60.5*(qeff**1.11))         ! Eqn. 20.
       Tsup=L/(L+C*(qeff**2))                       ! Zero baryon form of the transfer function (eqn. 18).
       ! Apply correction for scales close to horizon.
       if (fv.gt.0.0.and.fv.le.0.3) then
          qv=3.92*qEH*sqrt(Nv)/fv
          Bk=1.0+(1.2*(fv**0.64)*(Nv**(0.3+0.6*fv)))/((qv**(-1.6))+(qv**0.8))
       else
          qv=0.0
          Bk=1.0
       endif
       transfer_function=Tsup*Bk
    case (10)
       ! BBKS WDM matter transfer function (eqn. G6)
       kRfw=0.2*q*(omega0*h0**2/mwdm)**1.3333 ! Factor involving warm dark matter free streaming length.
       transfer_function=exp(-(kRfw/2.0)-((kRfw**2)/2.0))*transfer_function_BBKS_CDM(q)
    case default
       ! Unrecognized transfer function requested.
       write (0,*) 'transfer_function(): FATAL - transfer function not recognised!'
       write (0,*) '                     itrans = ',itrans
       call Report_Fatal_Error
    end select
    return
  end function Transfer_Function

  pure real function Transfer_Function_BBKS_CDM(q)
    !% Evaluates and returns the
    !% \href{http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=1986ApJ...304...15B&amp;db_key=AST&amp;data_type=HTML&amp;format=&amp;high=41aefb602231989}{BBKS
    !% cold dark matter transfer function} (their eqn. G3).
    implicit none
    real, intent(in) :: q

    transfer_function_BBKS_CDM=(log(1.0+2.34*q)/(2.34*q**2))/((1.0/q)**4+3.89/q**3+(16.1/q)**2+5.46**3/q+6.71**4)**0.25
    return
  end function Transfer_Function_BBKS_CDM

end module Power_Spectrum

