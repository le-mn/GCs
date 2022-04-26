module Numerical_Parameters
  ! Some numerical parameters in various galform routines

  !
  ! DECA - standard multipler
  real DECA,logDECA
  parameter (DECA=10.0)
  parameter (logDECA=1.0)

  !
  ! KILO - standard multipler
  real KILO,logKILO
  parameter (KILO=1.0e3)
  parameter (logKILO=3.0)

  !
  ! MILLI - standard multipler
  real MILLI
  parameter (MILLI=1.0e-3)

  !
  ! HECTO - standard multipler
  real HECTO,logHECTO
  parameter (HECTO=1.0e2)
  parameter (logHECTO=2.0)

  !
  ! MEGA - standard multipler
  real MEGA,logMEGA
  parameter (MEGA=1.0e6)
  parameter (logMEGA=6.0)

  !
  ! EPS6 is a small number
  !
  real EPS6
  parameter (EPS6=1.0e-6)

  !
  ! EPS3 is a (less) small number
  !
  real EPS3
  parameter (EPS3=1.0e-3)

  !
  ! EPS4 is another small number
  !
  real EPS4
  parameter (EPS4=1.0e-4)

  !
  ! EPSM is the mass below which a galaxy is considered to be defunct
  !
  real EPSM
  parameter(EPSM=1.0e-4)
  !
  ! EPSR is the smallest galaxy radius that we consider
  real EPSR
  parameter (EPSR=1.0e-6)

  !
  !   J2ERG - number of ergs in a Joule
  real J2ERG,logJ2erg
  parameter (J2ERG=1.0e7)
  parameter (logJ2ERG=7.0)

  !
  !
  ! RDISK_HALF_SCALE is the half-mass radius of an exponential disk
  ! in units of the disk scale-length
  !
  real RDISK_HALF_SCALE,EM_RDISK_HALF_SCALE_O2
  parameter (RDISK_HALF_SCALE=1.678346990)
  parameter (EM_RDISK_HALF_SCALE_O2=0.432067481) ! exp(-RDISK_HALF_SCALE/2)
  !
  ! SURF_DENS_NORM_HALF is exp(-RDISK_HALF_SCALE)/2/PI
  !
  !
  real SURF_DENS_NORM_HALF
  parameter (SURF_DENS_NORM_HALF=2.97114e-2)
  !
  !
  ! PI - the ratio of a circle's circumference to its diameter, plus
  ! some multiples of it
  real PI,PIO4,PI4,PIO2,PISQ,SQRT2PI,SQRT2OPI,SQRTPI,logPI
  parameter (PI=3.1415926536)
  parameter (logPI=0.497149873)
  parameter (PIO2=PI/2.0)
  parameter (PIO4=PI/4.0)
  parameter (PI4=4.0*PI)
  parameter (PISQ=PI*PI)
  parameter (SQRTPI=1.7724538509)
  parameter (SQRT2PI=2.5066282746)
  parameter (SQRT2OPI=0.7978845608)

  !
  !
  ! Various numerical constants
  real SQRT2,LN2,LN10,ISO_FAC,log4
  parameter (SQRT2=1.4142135624,LN2=0.6931471806,LN10=2.3025850930,ISO_FAC=2.5/LN10,log4=0.602059991)

  !
  ! GYR2S - the number of seconds in a Gyr (based on the Julian year
  ! of exactly 365.25 days - Allen's Astrophysical Quantities, page 15) 
  !
  real GYR2S,logGYR2S
  parameter (GYR2S=3.15576e16)
  parameter (logGYR2S=16.499103967)

  !
  ! GYR2YR - the number of years in a Gyr
  !
  real GYR2YR
  parameter (GYR2YR=1.0e9)

  !
  !
  ! MSOLAR - the mass of the Sun in kg (Allen's Astrophysical
  ! Quantities, page 12)
  !
  real MSOLAR,sqrtMSOLAR,MSOLAR_1030kg
  parameter (MSOLAR=1.9891e30)
  parameter (MSOLAR_1030kg=MSOLAR/1.0e30)
  parameter (sqrtMSOLAR=1.4103545653e15)

  !
  !
  ! MSOLAR_g - the mass of the Sun in g
  !
  real MSOLAR_g
  parameter (MSOLAR_g=MSOLAR*KILO)

  !
  !
  ! Ang2M - the number of metres in an Angstrom
  real Ang2M,logAng2M
  parameter (Ang2M=1.0e-10)
  parameter (logAng2M=-10.0)

  !
  !
  ! MPC2M - the number of metres in a Mpc (Particle Data Book 2002, page 6)
  real MPC2M,sqrtMPC2M
  parameter (MPC2M=3.0856775807e22)
  parameter (sqrtMPC2M=1.75660968365e11)

  !
  !
  ! PC2M - the number of metres in a pc (Particle Data Book 2002, page 6)
  real PC2M,logPC2M
  parameter (PC2M=3.0856775807e16)
  parameter (logPC2M=16.489350545)

  ! MPC2CM - the number of centimetres in a Mpc (Particle Data Book 2002, page 6)
  real MPC2CM
  parameter (MPC2CM=MPC2M*HECTO)

  ! PC2CM - the number of centimetres in a pc (Particle Data Book 2002, page 6)
  real PC2CM
  parameter (PC2CM=PC2M*HECTO)

  !
  ! G_SI - the gravitational constant in units of m^3/kg/s^2
  ! (Allen's Astrophysical Quantities, page 8)
  !
  real G_SI,sqrtG_SI
  parameter (G_SI=6.67259e-11)
  parameter (sqrtG_SI=8.16859228998e-6)

  !
  ! KM2M - the number of metres in a kilometre
  !
  real KM2M
  parameter (KM2M=1.0e3)

  !
  ! KMS2MPCGYR - convert velocity in km/s to Mpc/Gyr
  !
  real KMS2MPCGYR
  parameter (KMS2MPCGYR=KM2M*GYR2S/MPC2M)

  !
  ! G - the gravitational constant in units of (km/s)^2 Mpc/Msun
  !
  real G,sqrtG
  parameter (G=G_SI*MSOLAR/MPC2M/(KM2M**2))
  parameter (sqrtG=sqrtG_SI*sqrtMSOLAR/sqrtMPC2M/KM2M)

  !
  ! G_MPCGYR - the gravitational constant in units of km/s Mpc^2 Msun^-1 Gyr^-1
  !
  real G_MPCGYR
  parameter (G_MPCGYR=G_SI*MSOLAR*(GYR2S/MPC2M)/KM2M/MPC2M)


  !
  ! G_MPCGYR2 - the gravitational constant in units of Mpc^3 Msun^-1 Gyr^-2
  !
  real G_MPCGYR2
  parameter (G_MPCGYR2=G_SI*MSOLAR*((GYR2S/MPC2M)**2)/MPC2M)

  !
  ! G_GYRKMS3 - the gravitational constant in units of Gyr Msol^-1 km^3 s^-3
  !
  real G_GYRKMS3
  parameter (G_GYRKMS3=G_SI*MSOLAR/GYR2S/(KM2M**3))

  !
  ! H0100 - the Hubble constant in units of h km/s/Mpc
  !
  real H0100
  parameter (H0100=100.0)

  !
  ! H0100PGYR - the Hubble constant in units of h/Gyr
  !
  real H0100PGYR
  parameter (H0100PGYR=H0100*KM2M*GYR2S/MPC2M)

  !
  ! RHOCRIT - critical density of the Universe (3*H0^2/8*PI*G)
  ! in h^2 Msun/Mpc^3 
  real RHOCRIT,sqrtRHOCRIT
  parameter (RHOCRIT=3.0*(H0100**2)/8.0/PI/G)
  parameter (sqrtRHOCRIT=0.61237243570*H0100/sqrtPI/sqrtG)

  !
  !	k_Boltzmann - Boltzmann's constant in J/K (Particle Data
  !   Book 2002, page 5)
  real k_Boltzmann,logk_Boltzmann,sqrtk_Boltzmann,k_Boltzmann_erg
  parameter (k_Boltzmann=1.3806503e-23)
  parameter (k_Boltzmann_erg=k_Boltzmann*J2ERG)
  parameter (logk_Boltzmann=-22.859916308)
  parameter (sqrtk_Boltzmann=3.715710295489e-12)

  !
  !   M2CM - number of cm in a m
  real M2CM,logM2CM
  parameter (M2CM=100.0)
  parameter (logM2CM=2.0)

  !
  !   M32CM3 - number of cm^3 in a m^3
  real M32CM3
  parameter (M32CM3=M2CM**3)

  !
  !   CM32M3 - number of m^3 in a cm^3
  real CM32M3
  parameter (CM32M3=1.0/M32CM3)

  !
  !	c_light - speed of light in m/s (Allen's Astrophysical Quantities,
  !   page 8)
  real c_light,logc_light
  parameter (c_light=2.99792458e8)
  parameter (logc_light=8.47682070)

  !
  !	c_light_Angstroms - speed of light in Angstroms/s (Allen's Astrophysical Quantities,
  !   page 8)
  real c_light_Angstroms,logc_light_Angstroms
  parameter (c_light_Angstroms=c_light/Ang2M)
  parameter (logc_light_Angstroms=logc_light-logAng2M)

  !
  !	c_light_cm - speed of light in cm/s (Allen's Astrophysical Quantities,
  !   page 8)
  real c_light_cm
  parameter (c_light_cm=c_light*M2CM)

  !
  !	eV2ERG - convert electron volts to ergs
  real eV2ERG,logeV2ERG
  parameter (eV2ERG=1.60217733e-12)
  parameter (logeV2ERG=-11.7952894176)

  !
  !	eV2J - convert electron volts to Joules
  real eV2J
  parameter (eV2J=1.60217733e-19)

  !
  !   ERG2J - number of Joules in an erg
  real ERG2J
  parameter (ERG2J=1.0e-7)

  !
  !	h_Planck - Planck's constant in J s (Allen's Astrophysical Quantities,
  !   page 8)
  real h_Planck
  parameter (h_Planck=6.6260755e-34)

  !
  !	h_Planck_erg - Planck's constant in ergs s (Allen's Astrophysical Quantities,
  !   page 8)
  real h_Planck_erg
  parameter (h_Planck_erg=h_Planck*J2ERG)

  !
  !   LSUN - bolometric luminosity of the Sun in W (Allen's Astrophysical Quantities, page 340)
  real LSUN,logLSUN
  parameter (LSUN=3.845e26)
  parameter (logLSUN=26.5848963)

  !
  !   LSUN_ERG - bolometric luminosity of the Sun in erg/s (Allen's Astrophysical Quantities, page 340)
  real LSUN_ERG
  parameter (LSUN_ERG=LSUN*J2ERG)

  !
  !   LOG_LSUN_ERG - bolometric luminosity of the Sun in erg/s (Allen's Astrophysical Quantities, page 340)
  real LOG_LSUN_ERG
  parameter (LOG_LSUN_ERG=logLSUN+logJ2ERG)

  !
  !   LSUN_BC - bolometric luminosity of the Sun as used by Bruzual & Charlot [W]
  real LSUN_BC,logLSUN_BC
  parameter (LSUN_BC=3.862e26)
  parameter (logLSUN_BC=26.5868122)

  !
  !   LSUN40_BC - bolometric luminosity of the Sun as used by Bruzual & Charlot [10^40 ergs/s]
  real LSUN40_BC
  parameter (LSUN40_BC=3.862e-7)



  !
  !	sigma_Thomson - Thomson cross section for zero energy photons
  !   in m^2 (Particle Data Book 2002, page 5)
  real sigma_Thomson
  parameter (sigma_Thomson=6.65245854e-29)

  ! 
  !	M_Electron - mass of electron in kg (Particle Data Book, page 4)
  real M_Electron
  parameter (M_Electron=9.10938188e-31)

  ! 
  !	M_Atomic - mass of unit atomic weight in kg (12C=12 scale;
  !   Particle Data Book 2002, page 4)
  real M_Atomic,sqrtM_Atomic
  parameter (M_Atomic=1.66053873e-27)
  parameter (sqrtM_Atomic=4.07497083425e-14)

  ! 
  !	M_Atomic_g - mass of unit atomic weight in g (12C=12 scale)
  real M_Atomic_g
  parameter (M_Atomic_g=M_Atomic*KILO)

  !
  ! a_Radiation - radiation constant (J m^-3 K^4)
  real a_Radiation
  parameter (a_Radiation=8.0*(PI**5)*k_Boltzmann*((k_Boltzmann/c_light/h_Planck)**3)/15.0)

  !
  !	v_Halo_to_T_Virial - conversion factor from virial velocity (km/s) to
  !	virial temperature (K)
  real v_Halo_to_T_Virial
  parameter (v_Halo_to_T_Virial=0.5e6*M_Atomic/k_Boltzmann)

  !
  !	X_Hydrogen - mass fraction of hydrogen in primordial plasma
  real X_Hydrogen
  parameter (X_Hydrogen=0.778)

  !
  !	Y_Helium - mass fraction of helium in primordial plasma
  real Y_Helium
  parameter (Y_Helium=0.222)

  !
  !	Z_Metals - mass fraction of metals in primordial plasma
  real Z_Metals
  parameter (Z_Metals=5.36e-10)

  !
  !	X_Hydrogen_Solar - mass fraction of hydrogen in Solar plasma
  !   (Allen's Atrophysical Quantities, page 28)
  real X_Hydrogen_Solar
  parameter (X_Hydrogen_Solar=0.707)

  !
  !	Y_Helium_Solar - mass fraction of helium in Solar plasma
  !   (Allen's Atrophysical Quantities, page 28)
  real Y_Helium_Solar
  parameter (Y_Helium_Solar=0.274)

  !
  !	Z_Metals_Solar - mass fraction of metals in Solar plasma
  !   (Allen's Atrophysical Quantities, page 28)
  real Z_Metals_Solar
  parameter (Z_Metals_Solar=0.0189)

  !
  !	Atomic_Mass(Z) - mass of atom of atomic number Z in units
  !   of M_Atomic (Particle Data Book 2002, page 283)
  !	integer Z_Atomic
  !	real Atomic_Mass(2)
  !	data (Atomic_Mass(Z_Atomic),Z_Atomic=1,2) /1.00794,4.002602/
  real Atomic_Mass_Hydrogen,Atomic_Mass_Helium,sqrtAtomic_Mass_Hydrogen
  parameter (Atomic_Mass_Hydrogen=1.00794)
  parameter (sqrtAtomic_Mass_Hydrogen=1.00396215)
  parameter (Atomic_Mass_Helium=4.002602)

  !
  !	mu_Primordial - mean atomic weight for fully ionized plasma
  !   of primordial composition. 
  real mu_Primordial
  parameter (mu_Primordial=1.0/(2.0*X_Hydrogen/Atomic_Mass_Hydrogen+3.0*Y_Helium/Atomic_Mass_Helium))

  !
  ! Delta200 - mean density contrast of a halo at radius r_200
  real Delta200,sqrtDelta200
  parameter (Delta200=200.0)
  parameter (sqrtDelta200=14.142135623)

  !
  ! DeltaEdS - density contrast at virial radius of a halo in an Einstein-de Sitter universe
  real DeltaEdS
  parameter (DeltaEdS=18.0*PI*PI)

  !
  ! M8CRIT - Mass in a sphere of 8Mpc/h radius in a critical density Universe
  real M8CRIT
  parameter (M8CRIT=RHOCRIT*4.0*PI*(8.0**3)/3.0)

  !
  ! CRIT_FREQ=sqrt(4*pi*G*rhocrit) - characteristic gravitational frequency at critical density in units of s^-1
  real CRIT_FREQ
  parameter (CRIT_FREQ=2.0*SQRTPI*sqrtG*sqrtRHOCRIT*KM2M/MPC2M)

  !
  ! KOM=sqrt(k_B/m_H) in units of Mpc/s K^-1/2
  real KOM
  parameter (KOM=sqrtk_Boltzmann/sqrtAtomic_Mass_Hydrogen/sqrtM_Atomic/MPC2M)

  !
  ! KHORIZON - Defined as H_0/c in h Mpc^-1
  real KHORIZON
  parameter (KHORIZON=H0100*KM2M/c_light)

  !
  !	L_AB0 - luminosity density of a zero-magnitude object in the AB system [W/Hz] (Oke J. B., 1974, ApJS, 27, 21)
  real AB_NORM,L_AB0,logL_AB0,AB_MAG_NORM
  parameter (AB_MAG_NORM=48.60)
  parameter (AB_NORM=3.6307805477e-20)
  parameter (L_AB0=4.0*PI*((10.0*PC2M)**2)*AB_NORM*(M2CM**2)/J2ERG)
  parameter (logL_AB0=log4+logPI+(2.0*(1.0+logPC2M))+(-0.4*48.6)+(2.0*logM2CM)-logJ2ERG)

  !
  ! MPCKMS2GYR - convert units of Mpc/(km/s) to Gyr
  !
  real MPCKM2GYR
  parameter (MPCKM2GYR=MPC2M/GYR2S/KM2M)

  ! 
  ! GPI - G*pi in units of Gyr (km/s)^3 Msun^-1
  real GPI
  parameter (GPI=G*PI*MPC2M/GYR2S/KM2M)

  !
  ! Parameters used in figuring disk structure
  !
  ! kstrc_2 - constant relating disk angular momentum to r*V_c
  !   derived assuming a flat rotation curve
  !   (in Cedric's notes of August 1997 this has become kdisk)
  ! kstrc_1 - constant relating V_c(rdisk)^2 to GMdisk/rdisk
  !   in the disk plane (in Cedric's notes of August 1997 this has become k_h/2)
  real kstrc_1,kstrc_2
  parameter (kstrc_1=0.6251543028)
  parameter (kstrc_2=1.191648617)

  !
  ! fx_peak_NFW - the peak value of the function [ln(1+x)-x/(1+x)]/x as appears in the rotation curve of the NFW halo
  real fx_peak_NFW
  parameter (fx_peak_NFW=0.2162165956)

  !
  ! RAD2AS - convert radians to arcseconds
  real RAD2AS,logRAD2AS
  parameter (RAD2AS=180.0*60.0*60.0/PI)
  parameter (logRAD2AS=5.81157500587-logPI)

  !
  ! RAD2DEGREES
  real RAD2DEGREES
  parameter (RAD2DEGREES=180.0/PI)

  !
  ! Right_Angle_Degrees - number of degrees in a right angle
  real Right_Angle_Degrees
  parameter (Right_Angle_Degrees=90.0)

  !
  ! Semi_Circle_Degrees - number of degrees in a semi-circle
  real Semi_Circle_Degrees
  parameter (Semi_Circle_Degrees=180.0)

  !
  ! MPC2ASat10pc - converts Mpc to arcseconds for objects placed at a distance of 10pc
  real MPC2ASat10pc,logMPC2ASat10pc
  parameter (MPC2ASat10pc=(MEGA/DECA)*RAD2AS)
  parameter (logMPC2ASat10pc=logMEGA-logDECA+logRAD2AS)

end module Numerical_Parameters

module Atomic_Data
  !     
  ! Parameters of certain atomic lines
  !
  ! Floats
  real E_Halpha,E_Hbeta,E_Hydrogen_n1,N_Halpha_N_Lyc,N_Hbeta_N_Lyc,wH,wHeI,wHeII
  real logE_Hydrogen_n1,logE_Halpha,logE_Hbeta,logN_Halpha_N_Lyc,logN_Hbeta_N_Lyc
  !     
  ! Parameters
  parameter (wH=911.862)    ! Wavelength of hydrogen ionizing photon [Angstroms] (Allen's Astrophysical Quantities, p. 36).
  parameter (wHeI=504.319)  ! Wavelength of helium 1 ionizing photon [Angstroms] (Allen's Astrophysical Quantities, p. 36).
  parameter (wHeII=227.865) ! Wavelength of helium 2 ionizing photon [Angstroms] (Allen's Astrophysical Quantities, p. 36).
  parameter (E_Hydrogen_n1=13.59844) ! Energy of the hydrogen atom ground state [eV].
  parameter (logE_Hydrogen_n1=1.133489089)
  parameter (E_Halpha=E_Hydrogen_n1*((1.0/2.0)**2-(1.0/3.0)**2)) ! Energy of the Halpha line of hydrogen [eV].
  parameter (logE_Halpha=logE_Hydrogen_n1-0.857332496)
  parameter (E_Hbeta=E_Hydrogen_n1*((1.0/2.0)**2-(1.0/4.0)**2)) ! Energy of the Hbeta line of hydrogen [eV].
  parameter (logE_Hbeta=logE_Hydrogen_n1-0.7269987279)
  ! These ratios should be calcuable in principle, but I haven't found
  ! the relevant atomic calculation anywhere. So, the current values
  ! are just based on the numbers in the code from Alfonso
  ! Aragon-Salamanca's earlier calculation. Not very satisfactory!
  parameter (N_Halpha_N_Lyc=0.45) ! Number of Halpha photons produced per Lyman continuum photon.
  parameter (logN_Halpha_N_Lyc=-0.346787486)
  parameter (N_Hbeta_N_Lyc=0.116896) ! Number of Hbeta photons produced per Lyman continuum photon.
  parameter (logN_Hbeta_N_Lyc=-0.932200349468)
end module Atomic_Data
