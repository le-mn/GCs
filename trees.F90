program tree

 use Defined_Types               ! defined_types.F90
 use Cosmological_Parameters     ! parameter_modules.F90
 use Power_Spectrum_Parameters   ! parameters.F90
 use Tree_Memory_Arrays          ! memory_modules.F90
 use Tree_Memory_Arrays_Passable ! memory_modules.F90
 use Time_Parameters             ! parameters.F90
 use Tree_Routines               ! tree_routines.F90
 use Modified_Merger_Tree        ! modified_merger_tree.F90

implicit none

  type (TreeNode), pointer :: This_Node
  integer                  :: i,j,count,ntree
  integer, parameter       :: long = selected_real_kind(9,99)
  real, allocatable        :: wlev(:),alev(:), z_bin(:)
  integer, allocatable     :: ifraglev(:)!, count_GC_array(:)
  real                     :: mphalo,mres,ahalo,deltcrit,sigmacdm,zmax
  integer, parameter       :: nlev=180 !number of levels of the tree to store
  integer                  :: ierr,nhalomax,nhalo,nhalolev(nlev),jphalo(nlev),ilev
  integer                  :: iter,iseed0,iseed
  EXTERNAL deltcrit,sigmacdm,split
  real                     :: dc, M_form_GC, z_form_GC, z
  integer                  :: count_branch, count_GC, jlevelmax, n_z_bin
!  character                :: path
!  path = '/Users/ngoc/Desktop/project/GCs/code/'

! These cosmological parameters are taken from Creasey paper
 omega0  = 0.315      !code default = 0.25 
 lambda0 = 0.685      !code default = 0.75
 h0      = 0.67       !code default = 0.73
 omegab  = 0.0        !code default = 0.04
 Gamma   = omega0*h0  ! Omega_m.h  ignoring effect of baryons

! mass <- mass.h0 since unit in this code is [h^-1.M_sun] while mass [M_sun] in Creasey et. al. 2018
! Mass of halo for which the tree is to be grown [h^-1.M_sun]. The mass resolution of the tree and the number of trees to grow. 
 mphalo = 1.0e+12*h0   !halo mass at base of tree, default = 10^14
 mres   = 1.0e+07*h0   !mass resolution, default = 10^8
 ntree  = 1            !number of trees



! redshift and halo mass to form GC, use values in Creasey et. al. 2018  
 z_form_GC = 8.65
 M_form_GC = 1.0e+08*h0     

! Parameters of the Merger Tree Algorithm as defined in 
! Parkinson, Cole and Helly (2007arXiv0708.138 version 3 and in MNRAS paper)
! These values supercede the values given in the
! original astro-ph posting due to a small error in the code being
! identified. In this version of the code the error has been rectified
! and the fits redone. Using this code and these new parameters will
! produce near identical results to the old code with the old parameters.
! (passed in module Modified_Merger_Tree and Time_Parameters)

! G0      = 0.57
! gamma_1 = 0.38
! gamma_2 = -0.01
 eps1    = 0.1        
 eps2    = 0.1        

! parameters calibrated in Benson 2017
 G0      = 0.6353
 gamma_1 = 0.1761
 gamma_2 = 0.0411


! Cosmological and Power Spectrum parameters
! (passed in module  Cosmological_Parameters and Power_Spectrum_Parameters)

 pkinfile= 'pk_Mill.dat' !Tabulated Millennium Simulation linear P(k)
! itrans= -1             !indicates use transfer function tabulated in file pkinfile
  itrans= 1              !indicates use BBKS CDM transfer function with specified Gamma and Omega0
!  itrans= 2             !indicates use Bond & Efstathiou CDM transfer function with specified Gamma and Omega0
! itrans= 3              !indicates use Eisenstein and Hu CDM transfer function with specified Omega0, Omegab and h0
! CMB_T0 = 2.73          !For Eisenstein and Hu CDM transfer function one must specify the CMB temperature


!Set primordial P(k) parameters (ignored if itrans=-1)
 nspec  = 1.0           !primoridial power spectrum spectral index
 dndlnk = 0.0           !allow running spectral index by setting ne.0
 kref   = 1.0           !pivot point for running index
 sigma8 = 0.829         !power spectrum amplitude set regardless of other parameters, default = 0.9

 ierr    = 1            !initial error status to control make_tree()
 nhalomax= 0            !initialise
 nhalo   = 0
 iseed0  = -8635        !random number seed
 iseed   = iseed0

! Set up the array of redshifts at which the tree is to be stored
  write(0,*) 'The redshifts at which the tree will be stored:'
  allocate(wlev(nlev),alev(nlev),ifraglev(nlev))

! Specify output/storage times of the merger tree
  ahalo = 1.0         !expansion factor at base of tree
  zmax  = 9.0         !maximum redshift of stored tree, = 4. by default
  do ilev = 1,nlev    !tree levels uniform between z=0 and zmax
     alev(ilev) = 1.0/(1.0+zmax*real(ilev-1)/real(nlev-1))
     dc         = deltcrit(alev(ilev))
     write(0,'(a2,1x,f6.3,1x,a,f6.3)')'z=',(1/alev(ilev)) -1.0,'at which deltcrit=',dc
  end do


!Start generating trees
do i = 1, ntree
    iter = 1   !if we run out of allocated memory, which is flagged by ierr=1 or ierr=2 then we do another iteration with more allocated memory
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

     write(0,*)'made a tree',i
!write(0,*) 'Omega_0=',omega0,'Lambda_0=',lambda0
!write(0,*) 'sigma_8=',sigma8,'Gamma=',Gamma

!    Counting the number of nodes in the tree
     This_Node => MergerTree(1)
     count = 0
     do while (associated(This_Node))
        count = count + 1
        This_Node => Walk_Tree(This_Node)
     end do
     write(0,'(a,i3,a,i8)') 'number of nodes in tree',i,' is',count


!   Write out the information for the first couple of  halos in the tree
!     write(0,*) 'Example information from the tree:'
!     This_Node => MergerTree(1)
!     write(0,*) 'Base node:'
!     write(0,*) '  mass=',This_node%mhalo,' z= ',1.0/alev(This_node%jlevel)-1.0,' number of progenitors ',This_node%nchild
!     This_Node => This_node%child !move to first progenitor
!     write(0,*) 'First progenitor:'
!     write(0,*) '  mass=',This_node%mhalo,' z= ',1.0/alev(This_node%jlevel)-1.0
!     This_Node => This_node%sibling !move to 2nd progenitor
!     write(0,*) '  mass=',This_node%mhalo,' z= ',1.0/alev(This_node%jlevel)-1.0



!---------------------------------------------------------------------------!
!    My own code starts from here                                           !
!---------------------------------------------------------------------------!

!    Count the number of branches of this tree
     This_Node => MergerTree(1)
     print *, 'mass of halo at base', This_Node%mhalo/h0, '[M_sun]'
     count_branch = 0
     count_GC     = 0
     jlevelmax    = 1
     do while (associated(This_Node))
!        if (.not. associated(This_Node%child)) then        ! conditions are the same
        if (This_Node%nchild .eq. 0) then
             count_branch = count_branch+1
             jlevelmax = max(jlevelmax,This_Node%jlevel)
        endif
        This_Node => Walk_Tree(This_Node)
     end do
!     print *, 'number of branches in tree',i,' is',count_branch
!     print *, 'nlevel', jlevelmax
   

!    Count the total number of GCs can be formed in the tree
!    using the conditions from Creasey 2018: GCs can only formed inside halo of M >= 10^8 M_sun; z = 8.65
!    save in file the number count and redshift of halos which form GC

     open(unit = 101,file = '/Users/ngoc/Desktop/project/GCs/code/GCs-in-tree(Benson2017).txt',status='replace') 
!     open(unit = 101,file = '/Users/ngoc/Desktop/project/GCs/code/GCs-in-tree.txt',status='replace')        
     z = 0. 
     This_Node => MergerTree(1)
     do while (associated(This_Node))
        if (This_Node%mhalo .lt. M_form_GC) then 
            z = 1.0/alev(This_node%parent%jlevel)-1.0                
            if (z .ge. z_form_GC) then
                count_GC = count_GC + 1            
                ! write redshift of halos into file
                write(101,*) z
                !This_Node => Walk_To_Uncle(This_Node)
                This_Node => Walk_To_Next_Branch(This_Node)
            endif
        endif
        if (associated(This_Node)) then
            This_Node => Walk_Tree(This_Node)
        endif
     end do

     close(101)     
     print *, 'total number of primordial GCs in tree',i ,'is', count_GC









!---------------------------------------------------------------------------!
! My own code ends here                                                     !
!---------------------------------------------------------------------------!

end do

deallocate(wlev,alev,ifraglev)

end program tree
