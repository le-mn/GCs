!---------------------------------------------------------------------------------------!
!                                                                                       !
! To count number of GCs inside a DM halo of mass M at z = 0                            !
! Using the conditions as in Creasey18, count all the subhalos which can form GCs before!
! not only assign one at z_form_GC, M_form_GC                                           !               
!---------------------------------------------------------------------------------------!

program tree7

 use Defined_Types               ! defined_types.F90
 use Cosmological_Parameters     ! parameter_modules.F90
 use Power_Spectrum_Parameters   ! parameters.F90
 use Tree_Memory_Arrays          ! memory_modules.F90
 use Tree_Memory_Arrays_Passable ! memory_modules.F90
 use Time_Parameters             ! parameters.F90
 use Tree_Routines               ! tree_routines.F90
 use Modified_Merger_Tree        ! modified_merger_tree.F90
 use use_subroutines

implicit none

  type (TreeNode), pointer :: This_Node
  integer                  :: i,j,count,ntree
  integer, parameter       :: long = selected_real_kind(9,99)
  real, allocatable        :: wlev(:),alev(:), mphalo_array(:), GC_array(:)
  integer, allocatable     :: ifraglev(:)
  real                     :: mphalo,mres,ahalo,deltcrit,sigmacdm,zmax
  integer,parameter        :: nlev=401  !number of levels of the tree to store
  integer                  :: ierr,nhalomax,nhalo,nhalolev(nlev),jphalo(nlev),ilev
  integer                  :: iter,iseed0,iseed
  EXTERNAL deltcrit,sigmacdm,split
  real                     :: dc, M_form_GC, z_form_GC, z, z_max, z_step 
  integer                  :: jlevelmax, n_z_bin, n_mass_bin, k, i_count
  character(len=30)        :: f1, f2
  real                     :: mphalo_max, mphalo_min
  real                     :: omega0_Creasey, lambda0_Creasey, h0_Creasey, omegab_Creasey, ns_Creasey, sigma8_Creasey, &
                              & omega0_Mill, lambda0_Mill, h0_Mill, omegab_Mill, ns_Mill, sigma8_Mill, &
                              & omega0_COCO, lambda0_COCO, h0_COCO, omegab_COCO, ns_COCO, sigma8_COCO

!---------------------------------------------------!
! Cosmological parameters from different sources    !
!---------------------------------------------------!

! Creasey et al. 2018 (Planck)
omega0_Creasey  = 0.315
lambda0_Creasey = 0.685
h0_Creasey      = 0.673
omegab_Creasey  = 0.0
ns_Creasey      = 0.9603
sigma8_Creasey  = 0.829

! Millenium
omega0_Mill  = 0.25
lambda0_Mill = 0.75
h0_Mill      = 0.73
omegab_Mill  = 0.04
ns_Mill      = 1.0
sigma8_Mill  = 0.9

! COCO WDM
omega0_COCO  = 0.272
lambda0_COCO = 0.728
h0_COCO      = 0.704
omegab_COCO  = 0.04455
ns_COCO      = 0.967
sigma8_COCO  = 0.81


! Cosmological parameters taken from Creasey+BBKS, Mill or COCO
 omega0  = omega0_Creasey              
 lambda0 = lambda0_Creasey          
 h0      = h0_Creasey                
 omegab  = omegab_Creasey           
 Gamma   = omega0*h0  ! Omega_m.h  ignoring effect of baryons
 
  
! Mass of halo for which the tree is to be grown [h^-1.M_sun]. The mass resolution of the tree and the number of trees to grow. 
 mres   = 1.0e+07*h0   ! [M_sun/h] mass resolution
 ntree  = 1000       !number of trees

! generate an array of M_halo at z=0, [M_sun] to compare to Creasey, then take M_halo*h0 <-> [M_sun/h] to put into code 
 mphalo_min   = 2.0e+8
 mphalo_max   = 1.0e+13
 n_mass_bin   = 48
 allocate(mphalo_array(n_mass_bin), GC_array(n_mass_bin))
 call logspace(mphalo_min,mphalo_max,mphalo_array) 
 
 z_max     = 20.0
 z_step    = 0.05
 ! redshift and halo mass to form GC, use values in Creasey et. al. 2018  
 z_form_GC = 8.65
 M_form_GC = 1.0e+08*h0    !at unit [M_sun/h0] to compare to the code 



!---------------------------------------------------------------------------------------------------------------------------------------------!
! Parameters of the Merger Tree Algorithm as defined in Parkinson, Cole and Helly (2007arXiv0708.138 version 3 and in MNRAS paper)

! G0      = 0.57
! gamma_1 = 0.38
! gamma_2 = -0.01
 eps1    = 0.1        
 eps2    = 0.1        

! parameters calibrated in Benson 2017
 G0      = 0.6353
 gamma_1 = 0.1761
 gamma_2 = 0.0411

!----------------------------------------------------------------------------------------------------------------------------------------------!


! Power Spectrum parameters
! (passed in module  Cosmological_Parameters and Power_Spectrum_Parameters)
! pkinfile = 'powerspec/pk_WDM33_COCO_WMAP7.dat'
! pkinfile= 'pk_Mill.dat' !Tabulated Millennium Simulation linear P(k)
! pkinfile = 'Pk-Creasey18.txt'

! itrans= -1             !indicates use transfer function tabulated in file pkinfile
 itrans= 1              !indicates use BBKS CDM transfer function with specified Gamma and Omega0
! itrans= 2              !indicates use Bond & Efstathiou CDM transfer function with specified Gamma and Omega0
! itrans= 3               !indicates use Eisenstein and Hu CDM transfer function with specified Omega0, Omegab and h0
 CMB_T0 = 2.73          !For Eisenstein and Hu CDM transfer function one must specify the CMB temperature


! Set primordial P(k) parameters (ignored if itrans=-1)
 nspec  = ns_Creasey              !primoridial power spectrum spectral index
 dndlnk = 0.0                  !allow running spectral index by setting ne.0
 kref   = 1.0                  !pivot point for running index
 sigma8 = sigma8_Creasey          !power spectrum amplitude set regardless of other parameters, default = 0.9

! inital values for generating trees 
 ierr    = 1                   !initial error status to control make_tree()
 nhalomax= 0                   !initialise
 nhalo   = 0
 iseed0  = -8635               !random number seed
 iseed   = iseed0

! Set up the array of redshifts at which the tree is to be stored
  write(0,*) 'The redshifts at which the tree will be stored:'
  allocate(wlev(nlev),alev(nlev),ifraglev(nlev))

! Specify output/storage times of the merger tree
  ahalo = 1.0                 !expansion factor at base of tree
  zmax  = z_max               !maximum redshift of stored tree
  do ilev = 1,nlev            !tree levels uniform between z=0 and zmax
     alev(ilev) = 1.0/(1.0+zmax*real(ilev-1)/real(nlev-1))
     dc         = deltcrit(alev(ilev))
     write(0,'(a2,1x,f6.3,1x,a,f6.3)')'z= ',(1/alev(ilev)) -1.0,'at deltcrit =',dc
  end do

!--------------------------------------------!
! Start generating trees                     !
!--------------------------------------------!

do i = 1, ntree
    i_count  = 1 !to avoid M(z=0) = 1e8 M_sun at i_count = 1
    GC_array = 0
    do while (i_count .le. n_mass_bin) 
       iter = 1   !if we run out of allocated memory, which is flagged by ierr=1 or ierr=2 then we do another iteration with more allocated memory
       mphalo = mphalo_array(i_count)*h0  ! convert [M_sun] -> [M_sun/h] to put in the code

       do while (ierr.ne.0 .or. iter.eq.1) 
            if (iter.eq.1) iseed0 = iseed0-19  !advance seed for new tree
            iseed = iseed0
            !if needed increase the amount of memory allocated
            call Memory(nhalo,nhalomax,ierr,nlev,mphalo,mres)
            do j = 1, nhalomax, 1
               MergerTree_Aux(j)%index = j
            end do
            MergerTree => MergerTree_Aux  !Maps MergerTree to allocated 

       call make_tree(mphalo,ahalo,mres,alev,nlev,iseed,split,sigmacdm,deltcrit,&
                      & nhalomax,ierr,nhalo,nhalolev,jphalo,wlev,ifraglev)
       iter = iter+1
       end do

       write(0,*)'made tree',i, 'with M(z=0) = 1e',log10(mphalo_array(i_count)),'[M_sun]'

       
       k = 0
       z = 0.

       This_Node => MergerTree(1)
       do while (associated(This_Node))
          if (This_Node%mhalo .lt. M_form_GC .and. This_Node%parent%mhalo .ge. M_form_GC) then
              z = 1.0/alev(This_Node%parent%jlevel) - 1.0
              if (z .ge. z_form_GC) then
                  GC_array(i_count) = GC_array(i_count) + 1         
              endif
              This_Node => Walk_To_Next_Branch(This_Node)
          endif
          if (associated(This_Node)) then
              This_Node => Walk_Tree(This_Node)
          endif
       end do
       

       This_Node => MergerTree(1)
       do while (associated(This_Node)) 
          if (This_Node%mhalo .ge. M_form_GC .and. This_Node%nchild .eq. 0) then
              z = 1.0/alev(This_Node%jlevel) - 1.0
              if (z .ge. z_form_GC) then
                  GC_array(i_count) = GC_array(i_count) + 1         
              endif
              This_Node => Walk_To_Next_Branch(This_Node)
          endif
          if (associated(This_Node)) then
              This_Node => Walk_Tree(This_Node)
          endif
       end do

       i_count = i_count + 1
   enddo

   !save in files the number count and redshift of halos which form GC
   write(f1,*) i

   open(unit = 101,file = 'results_model/GCs-tree'//trim(adjustl(f1))//'-model-Creasey.txt',status='replace') 
   do k = 1, n_mass_bin 
       write (101,*) log10(mphalo_array(k)), GC_array(k)     !mhalo_array [M_sun] to compare to Creasey    
   enddo
   close(101)
 
   print *, GC_array 
end do

deallocate(wlev,alev,ifraglev)


end program tree7
