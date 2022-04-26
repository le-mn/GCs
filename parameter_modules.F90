module Stellar_Evolution
  !
  ! Module to hold details of stellar evolution
  !
  ! Integers
  integer ntsed(2),nzmetal
  !
  ! Logicals
  logical :: Output_ASCII_SSPs
end module Stellar_Evolution

module Cosmological_Parameters
  !
  ! Floats
  real h0,lambda0,omega0,omegab,CMB_T0
end module Cosmological_Parameters

module Run_Statistics
  ! Variables holding information on run-time statistics
  integer NI
  parameter (NI=60) ! Maximum number of iterations.
  !
  ! Integers
  integer icall(4),icall_sf,max_inode,nc1,nc2,ncall,nfail,nhist(20),niter_tab(NI),n_nelder_mead
#ifdef DEBUG
  data icall /0,0,0,0/
  data ncall /0/
  data nc1 /0/
  data nc2 /0/
#endif 
end module Run_Statistics

module IGM_Parameter_Data
  ! Variables associated with the IGM
  !
  ! Floats
  real VCUT,ZCUT
  !
  ! Logicals
  logical absorb
end module IGM_Parameter_Data

module Halo_Mass_Function
  !
  ! Integers
  integer massfun
  integer, parameter :: Halo_Mass_Function_PS=0,Halo_Mass_Function_SMT=1,Halo_Mass_Function_J2000=2
end module Halo_Mass_Function

module NBody_Galform
  !
  ! Module used to indicate to galform subroutines that this
  ! is a run with N-body merger trees and possibly substructure,
  ! in cases where this makes a difference.
  !
  ! Also used to pass N-body merger times to subroutines that need them,
  ! in order to minimise number of routines which must be changed.
  !
  ! Integers
  integer nbodyrun
  integer subdynfric
  integer nstepadd,nhalotree

  ! Snapshot redshifts from the simulation
  integer :: snap_first, snap_last
  real,    allocatable, dimension(:) :: snapshot_z
  !
  ! Floats
  real tcoolcut
  real mpart, lbox
  !
  ! Logicals
  logical trace_particles,use_nbody
  !
  ! Characters
  character*80 parttracefile
end module NBody_Galform


module Milli_Tree_Arrays

  use Kind_Numbers

  ! Array of redshifts from the tree file
  real, allocatable :: atrees(:)
  ! Number of steps
  integer, allocatable :: nhalolev(:)
  ! Halo properties
  integer, allocatable :: id(:),np(:),idesc(:),nsphalo(:),firstsub(:),halodecstep(:)
  ! Subhalo properties
  integer, allocatable :: subid(:),subnp(:),subidesc(:),subdecstep(:)
  real, allocatable :: pos(:,:),vel(:,:),subspin(:,:)
  real, allocatable :: veldisp(:),vmax(:)
  integer (kind=int8byte), allocatable :: mostboundid(:)
  real, allocatable :: rhalf(:)
  ! First halo at each step
  integer, allocatable :: firsthalo(:)
  ! Mapping between input and output halo arrays
  integer, allocatable :: mapin2out(:)
  integer, allocatable :: mapout2in(:)
  ! Galaxy merger times stored to save time if doing lots of outputs
  real, allocatable :: merger_time(:)
  integer nlev_trees_max
end module Milli_Tree_Arrays

module RCombined_Share
  logical :: Share_Use_Stars,Share_Use_Gas
end module RCombined_Share
