real function deltcrit(a)
  !     
  ! Subroutine to calculate critical overdensity for collapse at time t,
  ! for density field normalized at reference epoch a0=1.
  ! N.B. delta is value extrapolated from collapse epoch to a0=1
  !
  ! e.g. For Omega=1 deltcrit = 1.686 * (1+z)
  !     
  ! Notation:
  !  a = expansion factor rel to a0=1
  !  t = time relative to t0=1
  !  omega0 = omega at a0=1
  !  delc = critical delta at time t
  !  dldelcdlt = dln(delc)/dln(t)     
  !
  !
  ! Uses
  use Cosmological_Parameters
  use Numerical_Parameters
  implicit none
  !
  ! Array dimensions
  integer NTABLE,NV
  parameter (NTABLE=200,NV=1000)
  !
  ! Integers
  integer i,io,is,NSUM
  !
  ! Floats
  real a !,dldelcdlt
  REAL EPSOM,delc0,eta0,sh0,ch0,tomega,d0,ch,sh,eta,t,acosh
  REAL AMIN
  REAL omega0_save,lambda0_save
  REAL density,omflat(NV),delflat(NV),aflat(NTABLE),delta_flat(NTABLE),sum,dlin,dlin0,x,x0,xp,dxp,h,aa,lambda,omega
  !
  ! Parameters
  parameter (AMIN=0.1,EPSOM=1.0e-5,NSUM=2000)
  !
  ! Saves
  save aflat,delta_flat,lambda0_save,omega0_save
  !
  ! Data
  data omega0_save/0.0/
  data lambda0_save/0.0/
  !
  ! Code     
  if (abs(1.0-omega0).le.EPSOM) then ! Omega_0=1
     delc0=3.0*(12.0*PI)**(2.0/3.0)/20.0
     deltcrit=delc0/a
     !dldelcdlt=-2.0/3.0
  else if ((1.0-omega0).gt.EPSOM.and.lambda0.lt.EPSOM) then ! Omega_0<1 Lambda_0=0
     ! Calculate properties at t=t0.
     eta0=acosh(2.0/omega0-1.0)
     sh0=sinh(eta0)
     ch0=cosh(eta0)
     tomega=2.0*PI/(sh0-eta0)
     d0=3.0*sh0*(sh0-eta0)/(ch0-1.0)**2-2.0 ! Linear growth factor.
     ! Calculate properties at expansion factor a.
     ch=a*(ch0-1.0)+1.0
     eta=acosh(ch)
     sh=sinh(eta)
     t=(sh-eta)/(sh0-eta0)
     deltcrit=1.5*d0*(1.0+(tomega/t)**(2.0/3.0))
     !dldelcdlt=-2.0/3.0/(1.0+(t/tomega)**(2.0/3.0))
  else ! General Omega_0 and Lambda_0
     if (omega0.ne.omega0_save.or.lambda0.ne.lambda0_save) then
#ifdef INFO
        write (0,*) 'deltcrit(): INFO - constructing look-up table for deltacrit(a)'
        write (0,*) '            for a flat Omega+Lambda=1 cosmology.'
        write (0,*) '            Note: derivative  dldelcdlt not implemented.'
#endif
        omega0_save=omega0
        lambda0_save=lambda0
        if (abs(omega0+lambda0-1.0).gt.EPSOM) stop 'Error: Omega_0+Lambda_0.ne.1'
        !          
        ! Read Vince Eke's file that tabulates deltcrit0 against omega0.
        open (10,file='Data/flat.data',status='old')
        read (10,*) ! Skip header.
        do i=1,NV
           read (10,*) omflat(i),density,delflat(i)
#ifdef DEBUG
           write (0,*) i,omflat(i),density,delflat(i)
#endif
        end do
        close (10)
        !
        ! Evaluate constant required to normalize the linear growth
        ! factor.
        x0=(2.0*(1.0/omega0-1.0))**0.333333
        sum=0.0
        dxp=x0/float(NSUM)
        do is=1,NSUM
           xp=x0*(float(is)-0.5)/float(NSUM)
           sum=sum+((xp/(xp**3+2.0))**1.5)*dxp
        end do
        dlin0=sum*sqrt(x0**3+2.0)/sqrt(x0**3)
        ! Tabulate deltcrit versus a for the specified values of
        ! omega0 lambda0. Spacing in a is linear in order to enable
        ! quick look up.
        do i=1,NTABLE
           aa=AMIN+(1.0-AMIN)*float(i-1)/float(NTABLE-1)
           aflat(i)=aa
           lambda=lambda0/(lambda0+(1.0-omega0-lambda0)/aa**2+omega0/aa**3)
           omega=omega0*lambda/(lambda0*aa**3)
           x=x0*aa
           sum=0.0
           dxp=x/float(NSUM)
           do is=1,NSUM
              xp=x*(float(is)-0.5)/float(NSUM)
              sum=sum+((xp/(xp**3+2.0))**1.5)*dxp
           end do
           dlin=(sum*sqrt(x**3+2.0)/sqrt(x**3))/dlin0
           call locate(omflat,NV,omega,io)
           if (io.lt.NV) then
              h=(omflat(io+1)-omega)/(omflat(io+1)-omflat(io))
              delta_flat(i)=(delflat(io)*h+delflat(io+1)*(1.0-h))/dlin
           else
              delta_flat(i)=delflat(NV)/dlin
           end if
        end do
     end if
     !
     ! Evaluate deltcrit using look-up table
     ! allow for a slightly greater than 1.
     if (a.gt.AMIN.and.a.le.1.001) then
        i=1+int((a-AMIN)*float(NTABLE-1)/(1.0-AMIN))
        i=min(NTABLE-1,i)
        h=(aflat(i+1)-a)/(aflat(i+1)-aflat(i))
        deltcrit=delta_flat(i)*h+delta_flat(i+1)*(1.0-h)
        !dldelcdlt=1.0 ! Not implemented.
     else if (a.le.AMIN) then
        ! Extrapolate to higher redshift by approximating to
        ! Omega=1 at high redshift. 
        deltcrit=delta_flat(1)*AMIN/a
        !dldelcdlt=1.0 ! Not implemented.
     else
        write (0,*) 'deltcrit(): FATAL - look-up table only for a<1'
        write (0,*) '            a = ',a
        stop
     end if
  end if
  return
end function deltcrit



         
         
