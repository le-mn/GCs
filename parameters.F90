! Prompt for and read in all model parameters.
! Many are then distributed to the subroutines that need to know them
! via a set of modules.

module PreEnrich_Parameters
  real :: PreEnrich_DeltaMMin
  real :: PreEnrich_Z
  logical :: PreEnrich_DeltaMMin_Set=.false.
  logical :: PreEnrich_Z_Set=.false.
end module PreEnrich_Parameters

module Halo_Mass_List
  integer nmass,ntreemin,ntreemax,Halo_Mass_List_N_Halo
  integer :: Halo_Mass_List_Method
  logical :: Halo_Mass_List_Method_Set=.false.
  real ahalo,mlow,mfac,volume,xseek,lmlow
  real, allocatable :: Halo_Mass_List_Halo_Masses(:),Halo_Mass_List_Halo_Weights(:)
  integer, allocatable :: Halo_Mass_List_N_Trees(:)
  integer NMASSMAX,NTREESMAX
  parameter (NMASSMAX=1000,NTREESMAX=1000)
end module Halo_Mass_List

module Mass_to_Light
  real upsilon,upsilon_burst
end module Mass_to_Light

module SPH_Mimic
  !
  ! SPH comparison model parameters
  !
  ! Integers
  integer Nsph
  !
  ! Floats
  real mgas_sph,zmethard
  !
  ! Logicals
  logical nofb,SPHrun
end module SPH_Mimic

module Emission_Lines
  ! Name of the emission lines file
  !
  ! Characters
  character emlinefile*80
end module Emission_Lines

module Morphology_Parameters
  ! Variables associated with morphology calculations
  !
  ! Floats
  real fellip,fburst,btburst,fgasburst
  !
  ! Logicals
  logical idisk
end module Morphology_Parameters

module Merging_Parameters
  !
  ! Variables holding parameters associated with galaxy merging
  !
  ! Floats
  real tau0mrg,alphamrg
  !
  ! Logicals
  logical dyn_fric
end module Merging_Parameters

module Halo_Parameters
  !
  ! Variables holding details of halo density profiles
  !
  ! Integers
  integer iprofile
  !
  ! Floats
  real aprofile,anfw,core,c0,frac_dm,contractatol
  real, parameter :: ALPHA_ROT=0.0
  !
  ! Logicals
  logical dynamic_halo,use_formn_halo_props
end module Halo_Parameters

module Cooling_Table_and_Parameters
  !
  ! Array dimensions
  integer nmetals,nvalues
  parameter (nmetals=8,nvalues=90)
  !
  ! Integers
  integer icool
  !
  ! Floats
  real alpha_cool,loglam(nvalues,nmetals),logtc(nvalues,nmetals),logtemp(nvalues),VCONDUCTION,xe(nvalues,nmetals),zmetal(nmetals)&
       & ,qcore
  !
  ! Logicals
  logical frfall,allcool,starvation,Output_Feedback
end module Cooling_Table_and_Parameters

module Power_Spectrum_Parameters
  ! Variables used to hold properties of the power spectrum
  !
  ! Array dimensions
  integer Transfer_Function_Table_N_Max
  parameter (Transfer_Function_Table_N_Max=1800)
  !
  ! Integers
  integer igwave,ireset,itrans,nktab,NKTABMAX,NSPL,Trans_Func_Table_N_Points
  !
  ! Array dimensions
  parameter(NKTABMAX=1000)
  !
  ! Floats
  real dndlnk,gamma,kref,lnktab(NKTABMAX),lnpktab(NKTABMAX),mwdm,nspec,sigma8,scla,sclm
  real Transfer_Function_Table_lnk(Transfer_Function_Table_N_Max),Transfer_Function_Table_lnTk(Transfer_Function_Table_N_Max)
  !
  ! Logicals
  logical WDMrun
  !
  ! Characters
  character pkinfile*1024,splinefile*220,tffile*1024
  !
  ! Parameters
  parameter (NSPL=200)
end module Power_Spectrum_Parameters

module Mass_Conservation
  !
  ! Floats
  real :: Mass_Conservation_Tolerance=1.0e-2
  !
  ! Logicals
  logical :: Mass_Conservation_Tolerance_Set=.false.
end module Mass_Conservation

module Sizes_Parameters
!
!     Floats
      real f_orbit,spin_disp,spin_med
!
!     Logicals
      logical selfgravity
end module Sizes_Parameters

module Star_Formation_Feedback_Prmtr
  !
  ! Integers
  integer Star_Formation_Model
  !
  ! Floats
  real alphahot,alphastar,efold,fdyn,fsw0_burst,fsw0_disk,pstar,SF_critical_density,tau_star_min,tburst_min,epsilon_Star,tau0star&
       &,vhot_burst ,vhot_disk,VSTAR,vsw_burst ,vsw_disk
  !
  ! Logicals
  logical feedback,instant,tdisk,vc_bulge,vc_disk
  logical dstab_sf_law
  !
  !RGB7: added alpha_reheat - controls reheat timescale. an optional parameter.
  real alpha_reheat
  logical alpha_reheat_set
  !RGB8: added to control abrubtness of AGN feedback cutoff.
  !
  real :: can_cool_efac=0.0
  logical :: can_cool_efac_set=.false.
  !
  !RGB11: added to prevent vc used in feedback becoming excessive.
  real :: vcirc_fac=10000.0
  logical :: vcirc_fac_set=.false.
  !
  !!! RGB12 : added mcrit_fac to implement Croton et al surface brightness
  !   threshold.
  !
  !real :: mcrit_fac = 0.0
  !logical ::  mcrit_fac_set = .false.
  !
  ! Output format for Grasil
  logical :: GRASIL2_Output
  !
end module Star_Formation_Feedback_Prmtr

module Disk_Stability_Parameters
  !
  ! Floats
  real stable_disk
  real :: Disk_Stab_Transfer_Frac
  !
  ! Logicals
  logical :: Disk_Stab_Transfer_Frac_Set=.false.
  logical :: Disk_Stability_Grow_SMBH
  logical :: Disk_Stability_Grow_SMBH_Set=.false.
  logical :: Disk_Stability_Large_Radii
end module Disk_Stability_Parameters

module Satellite_Orbit_Model_Prmtr
  !
  ! Logicals
  logical EXtras
  !
  ! Floats
  real epsilonh,fdiskfac,fhalofac,frheat,frmerge
end module Satellite_Orbit_Model_Prmtr

module Time_Parameters
  ! Parameters used in making binary splits in the merger tree
  real eps1,eps2
  integer istep
end module Time_Parameters





