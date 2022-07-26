!
! split(m2,w,mmin,sigma,iseed,dwmax,dw,nprog,mprog) 
!
! This version has been modified to allow a factor
!   G0 (sigma(m1)/sigma(m2)^gamma_1 (delta_c2/sigma(m2))^gamma_2
! to be inserted into the normal Press-Schechter expression for dn/dq .
! The parameters G0, gamma_1 and gamma_2 are specified in the module
!  modified_merger_tree.F90 .
! The subroutine has w as an extra argument and so a corresponding change 
! has to be made to make_tree.F90
!
! Subroutine to do single binary split.
! 
! Input:
! ------
! m2:    input mass
! mmin:  minimum mass to resolve
! iseed: random number seed
! dwmax: max timstep dw (w=deltc(a))
! sigma(m,dlsigdlm): FUNCTION to compute sigma(m) & dln(sigma)/dln(m)
!
! Output:
! -------
! dw:   actual time step in units of change in delta_c  
! nprog: (=0,1,2) number of progenitors with m>mmin after split
! mprog(2): array containing masses of progenitors with m>mmin
! 
! ----------------------------
! Theory:
! 
! Defining q=m/m2, qmin=mmin/m2
! (1) generate 0 or 1 progenitors with qmin<q<1/2 according to the 
! probability distribn dn/dq given by PS
! (2) generate the complement q'=1-q
! (3) compute the mean mass fraction f in progenitors with q<qmin according
! to PS, and subtract this from the larger progenitor, so final q'=1-q-f
! (4) discard any progenitors with q<qmin
! 
! dn/dq  = (2/pi)^(1/2) alpha(q*m2)/q^2  dw
! * sigma(q*m2)^2/(sigma(q*m2)^2-sigma(m2)^2)^(3/2) 
! * G0 * (sigma(q*m2)/sigma(m2))^gamma_1 * (delta_c2/sigma(m2))^gamma_2
! 
! where alpha(m) = -dln(sigma)/dln(m)
!
! and can be rewrite as the multiplication  
! dn/dq = P(q) * R(q)  where
! 
! P(q) = (2/pi)^(1/2)*alpha(m2/2)*B*dw /  q^(2-beta+gamma_1*mu) 
!       * G0 * (delta_c2/sigma(m2))^gamma_2                
! 
! (NB we call the exponent beta-1-gamma_1*mu = eta )
! and 
! 
! R(q) = 1/(B*q^beta) *sigma(q*m2)^2/(sigma(q*m2)^2-sigma(m2)^2)^(3/2)
! * alpha(q*m2)/alpha(m2/2)
! * [q^mu*sigma(q*m2)/sigma(m2)]^gamma_1
! 
! dn/dq is written this way so that P(q) is a simple power-law
! while R(q)<1 for 0<q<1/2 . 
! 
! If alpha(m)>0 and d(alpha)/dm>0, then the factor
! v(q) = sigma(q*m2)^2/(sigma(q*m2)^2-sigma(m2)^2)^(3/2)
! is monotonically increasing and ln(v) vs ln(q) is concave upwards 
! for 0<q<1/2, so v(q) is bounded from above by a power law B*q^beta
! chosen to equal v(q) at q=qmin and q=1/2
! 
! alpha(q m2)/alpha(m2/2)  is also less than unity for q<1/2 .
! as is
! [q^alpha_prime*sigma(q*m2)/sigma(m2)]^gamma_1
!
! if mu=alpha(m2/2) for gamma_1 >0
! and mu=ln(sigma(qmin m2)/sigma(m2/2)) / ln(qmin) for gamma_1<0
!
! --------------------------------------
! 
! Given a halo mass m2 and minimum mass mmin:
! 
! 1) Compute some values that are used repeatedly
! 
! 2) Find the power-law B*q^beta used in R(q)
! 
! 3) Compute a "timestep" dw  (dw = deltac1-deltac2 )
! subject to two constraints specified by parameters eps1 and eps2
! such that in this interval there is small but finite probablity
! that the halo will have formed from the merger of two halos both 
! of mass greater than mmin, and also the constraint dw<dwmax.
! specifically:
! (i) dw < eps1 *sqrt(2*(sigma(m2/2)^2-sigma(m2)^2))
! to ensure dn/dq propto dw
! (ii) nprog_av(qmin<q<1/2)(before rejection) < eps2 
! (iii) dw <= dwmax to allow output at a specified w
! 
! 4) Compute f, the mean mass fraction in progenitors m<mmin
! 
! 5) Decide whether to create a fragment.
!
! If yes:
! 
!    6) Generate a q from the power-law P(q)
! 
!    7) Keep or reject this q stochastically with probability R(q)
! 
!    8) If have generated a fragment q, then other fragment is q' = 1-q-f
!
! If not, q' = 1-f:
! 
!    9) Discard q' if q'<qmin
! 
subroutine split(m2,w,mmin,sigma,iseed,dwmax,dw,nprog,mprog) 
  !
  ! Uses
  use Numerical_Parameters
  use Time_Parameters
  use Modified_Merger_Tree
#ifdef DEBUG
  use Run_Statistics
#endif
  implicit none
  !
  ! Integers
  integer iseed,nprog
  !
  ! Floats
  REAL m2,mmin,dwmax,dw,mprog(2)

  !
  ! Functions
  real ran3,sigma
  !
  ! Externals
  external sigma
  !
  ! Common blocks
  REAL :: dlsigdlm,sigsq_m2,sigsq_hf,qmin,diff_qmin,diff12_qmin
  REAL :: diff32_qmin,v_qmin,diff_hf,diff12_hf,diff32_hf,v_hf,sfac
  REAL :: two_pow_beta ,beta ,eta_inv,half_pow_eta,qminexp,ffac,dn_dw,n_av,f
  REAL :: random ,q,alpha_q,sigsq_q,r_q,v_q,mtemp 
  REAL :: sig_m2,sig_hf,sig_q,diff12_q,b,dw_eps2 ,alpha_hf
  REAL :: mu,w,gfac0,gfac1,eta,q_pow_eta,z,J_UNRESOLVED
  !
  ! Parameters
  REAL, parameter :: EPSETA=1e-6,EPSQ=6.0e-6
  ! 
  ! Saves
  REAL, save :: mminlast,sigsq_qmin,sig_qmin
  !
  ! Data
  data mminlast /0.0/
  !
  ! Code
#ifdef DEBUG
  write (0,*) 'split(): DEBUG - start'
#endif
  ! 0) Since mmin will typically remain fixed we need only compute
  ! sigma(mmin) and slope alpha on the first call of this routine.
  ! 
  if (mmin.ne.mminlast) then  
     dlsigdlm=0.0 ! Set to avoid compiler warning
     sig_qmin=sigma(mmin,dlsigdlm)
     sigsq_qmin=sig_qmin**2 
     mminlast=mmin
  end if
#ifdef DEBUG
  ncall=ncall+1 ! Count number of calls to this routine
#endif
  !
  ! 1) Compute some useful numbers.
  ! 
  sig_m2=sigma(m2,dlsigdlm) ! sigma(m2)**2 and slope alpha
  sigsq_m2=sig_m2**2 
  sig_hf=sigma(0.5*m2,dlsigdlm) ! sigma(m2/2) and slope alpha_hf
  sigsq_hf=sig_hf**2
  alpha_hf=-dlsigdlm
  qmin=mmin/m2 ! Minimum mass ratio qmin.
  gfac0=G0 * ((w/sig_m2)**gamma_2) 
  ! 
  if (qmin.lt.(0.5-EPSQ)) then ! Generate fragment with qmin<q<1/2
!    More useful numbers needed only when qmin>1/2
     if (gamma_1 .le. 0.0) then
        mu =-log(sig_qmin/sig_hf) / log(2.0*qmin) 
     else
        mu=alpha_hf
     end if
     gfac1= gfac0 * ((sig_hf/sig_m2)**gamma_1) / (2.0**(mu*gamma_1))
     ! 
     ! 2) Find B and beta such that B*q^beta = v(q) at q=qmin and q=1/2
     !    where v(q) = sigma(q*m2)^2/(sigma(q*m2)^2-sigma(m2)^2)^(3/2).
     ! 
     diff_qmin=sigsq_qmin-sigsq_m2
     diff12_qmin=sqrt(diff_qmin)
     diff32_qmin=diff_qmin*diff12_qmin
#ifdef DEBUG
     if (diff32_qmin.le.0.0) stop 'split(): DEBUG/FATAL - diff32_qmin<=0 !'
#endif
     v_qmin=sigsq_qmin/diff32_qmin
     ! 
     diff_hf=sigsq_hf-sigsq_m2
     diff12_hf=sqrt(diff_hf)
     diff32_hf=diff_hf*diff12_hf
#ifdef DEBUG
     if (diff32_hf.le.0) stop 'split(): DEBUG/FATAL - diff32_hf<=0 !'
#endif
     v_hf=sigsq_hf/diff32_hf
     ! sfac=sqrt(2*(sigma(m2/2)^2 - sigma(m2)^2)), used in (3) below
     sfac=SQRT2*diff12_hf   
     ! 
#ifdef DEBUG
     if (v_hf.le.0) stop 'split(): DEBUG/FATAL v_hf<=0 !'
     if (qmin.le.0) stop 'split(): DEBUG/FATAL qmin<=0 !'
#endif
     beta=log(v_qmin/v_hf)/log(2.0*qmin)
     if (beta.le.0.0) then
        write (0,*) 'split(): FATAL - beta<0'
        write (0,*) '                     qmin = ',qmin,', log(2*qmin) = ',log(2.0*qmin)
        write (0,*) '                   v_qmin = ',v_qmin
        write (0,*) '                     v_hf = ',v_hf
        write (0,*) '         log(v_qmin/v_hf) = ',log(v_qmin/v_hf)
        write (0,*) '                     beta = ',beta
        write (0,*) '         This can probably be avoided by increasing EPSQ.'
        stop
     end if
     two_pow_beta=2.0**beta
     b=v_hf*two_pow_beta
     eta=beta-1.0-gamma_1*mu
#ifdef DEBUG
     if (two_pow_beta.le.0) stop 'split(): DEBUG/FATAL - two_pow_beta<=0 !'
#endif
     half_pow_eta=2.0**(-eta)! 2^(-eta)
     ! 
     ! 3) Set time step dw and compute mean number of progeny before rejection.
     !    N.B.: special case eta=0
     if (abs(eta).gt.EPSETA) then ! eta != 0
        eta_inv=1.0/eta
        qminexp=qmin**eta         ! Used again in section 6)
        ffac=half_pow_eta-qminexp ! Used again in section 6)
        dn_dw=SQRT2OPI*alpha_hf*b*eta_inv*ffac*gfac1
     else ! beta = 1
        dn_dw=-SQRT2OPI*alpha_hf*b*log(2.0*qmin)*gfac1 
     end if
     ! 
     if (dn_dw.gt.0.0) then
        dw_eps2=eps2/dn_dw
     else
        dw_eps2=dwmax
     end if
     dw=min(eps1*sfac,dw_eps2,dwmax) 
     n_av=dn_dw*dw
     ! 
     ! 4) Compute the mean fraction of mass in objects of mass less than mmin.
#ifdef DEBUG
     if (diff12_qmin.le.0) stop 'split(): DEBUG/FATAL - diff12_qmin<=0 !'
#endif
     z=sig_m2/diff12_qmin    
     f=SQRT2OPI*dw * gfac0 * J_UNRESOLVED(z) /sig_m2
     ! 
     ! 5) Randomly choose one or zero progeny depending on n_av.
     random=ran3(iseed)
     if (random.le.n_av) then
#ifdef DEBUG
        nc1=nc1+1 
#endif
        ! 6) Select a value of q with a power-law q**(eta-1) distribution.
        !    N.B.: eta=0 is special case.
        random=ran3(iseed)
        ! 
        if (abs(eta).gt.EPSETA) then ! eta != 0
           q_pow_eta=qminexp+random*ffac ! q^eta
           q=q_pow_eta**eta_inv
        else ! eta = 0
           q=qmin*(2.0*qmin)**(-random)
        end if
        ! 7) Compute rejection probability.
        sig_q=sigma(q*m2,dlsigdlm) ! sigma(q m2), alpha_q
        sigsq_q=sig_q**2
        alpha_q=-dlsigdlm
        diff12_q=sqrt(sigsq_q-sigsq_m2)
        v_q=sigsq_q/diff12_q**3
        ! Acceptance probability=R(q)<1
        r_q=(alpha_q/alpha_hf)*((sig_q*(2.0*q)**mu/sig_hf)**gamma_1) * &
 &                v_q/(b*q**beta)
        ! 
#ifdef DEBUG
!      If the assumptions regarding the behaviour of sigma(m) that 
!      are stated in the PCH appendix are true and the code is implemented
!      correctly then r_q should always be <1 (to within rounding) error.
!
!      If this is not the case print out some warning and diagnostic information
       if (r_q.gt.1.00001) then
          write(0,*) 'R(q)=',r_q,'>1!!:'
          write(0,*) 'alpha_q/alpha_hf,v_q/(b*q**beta), &
 &        (sig_q*(2q)**mu/sig_hf)**gamma_1',alpha_q/alpha_hf,v_q/(b*q**beta),&
 &        (sig_q*(2.0*q)**mu/sig_hf)**gamma_1
          if (alpha_q/alpha_hf.gt.1.00001) then
             write(0,*) 'alpha_q,alpha_hf,q,m2=',alpha_q,alpha_hf,q,m2
          end if   
          if (v_q/(b*q**beta).gt.1.00001) then
             write(0,*) 'v_q,b,q,beta=',v_q,b,q,beta
          end if   
          if ((sig_q*(2.0*q)**mu/sig_hf)**gamma_1.gt.1.00001) then
             write(0,*) 'sig_q,q,mu,sig_m2,gamma_1=',sig_q,q,mu,sig_m2,gamma_1
          end if   
       end if
#endif
        random=ran3(iseed)
        if (random.ge.r_q) q=0.0 ! Reject fragment
#ifdef DEBUG
        if (random.lt.r_q) nc2=nc2+1 ! Count number of acceptances.
#endif
     else
        q=0.0 ! No fragment.
     end if
     ! 
     ! Calculate mass of other progenitor in pair, and count how many with m>mmin.
     mprog(2)=q*m2
     if (mprog(2).le.mmin) then
        nprog=1
     else
        nprog=2
     end if
     ! 
     mprog(1)=(1-q-f)*m2 ! Adjust for mass frac in progenitors m<mmin.
     if (mprog(1).le.mmin) then
        nprog=nprog-1
        mprog(1)=mprog(2)
        mprog(2)=0.0
     else if  (mprog(1).lt.mprog(2)) then!can happen if 1-q-f<q, i.e. q>(1-f)/2
        mtemp=mprog(1)
        mprog(1)=mprog(2)
        mprog(2)=mtemp
     end if
     ! 
  else ! qmin>1/2 and so only progenitors are unresolved
     diff_hf=sigsq_hf-sigsq_m2
     diff12_hf=sqrt(diff_hf)
     sfac=SQRT2*diff12_hf   
     dw=min(eps1*sfac,dwmax) !ensure linear dependence on dw and < dwmax
     diff_qmin=sigsq_qmin-sigsq_m2
     if (diff_qmin.lt.0.0) then
        ! diff_qmin<0 which should be impossible but may have
        ! occured as within rounding error m2=mmin.
        diff_qmin=0.0
     end if
     diff12_qmin=sqrt(diff_qmin)
     if (diff12_qmin.gt.SQRT2OPI*dw) then
#ifdef DEBUG
        if (diff12_qmin.le.0) stop 'split(): DEBUG/FATAL diff12_qmin<=0 !' 
#endif
        z=sig_m2/diff12_qmin    
       f=SQRT2OPI*dw * gfac0 * J_UNRESOLVED(z) /sig_m2
     else
        f=1.0
     end if
     mprog(1)=(1-f)*m2 ! Adjust for mass frac in progenitors m<mmin.
     if (mprog(1).gt.mmin) then
        nprog=1
     else
        nprog=0
        mprog(1)=0
     end if
     mprog(2)=0
  end if
#ifdef DEBUG
  write (0,*) 'split(): DEBUG - done'
#endif
  return
end subroutine split

