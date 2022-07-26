!-------------------------------------------------------------!
!                                                             !
! This is a copy of program tree (in file trees.F90)          !
! to test the code for counting number of GCs which is        !
! a function of halo mass at redshift 0                       !
! N_GCs(M(z=0))                                               !
!                                                             !
!-------------------------------------------------------------!

program tree2

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
  integer, allocatable     :: ifraglev(:), count_GC_array(:)
  real                     :: mphalo,mres,ahalo,deltcrit,sigmacdm,zmax
  integer, parameter       :: nlev = 400!180 !number of levels of the tree to store
  integer                  :: ierr,nhalomax,nhalo,nhalolev(nlev),jphalo(nlev),ilev
  integer                  :: iter,iseed0,iseed
  EXTERNAL deltcrit,sigmacdm,split
  real                     :: dc, M_form_GC, z_form_GC, z, z_max, z_step 
  integer                  :: count_branch, count_GC, jlevelmax, n_z_bin, kk, ii, mphalo_convert
  character(len=30)        :: f1, f2
  real                     :: mphalo_max, mphalo_min, mphalo_step

! Cosmological parameters taken from Mill
 omega0  = 0.25 
 lambda0 = 0.75
 h0      = 0.73
 omegab  = 0.04
 Gamma   = omega0*h0  ! Omega_m.h  ignoring effect of baryons

  
! mass <- mass.h0 since unit in this code is [h^-1.M_sun] while mass [M_sun] in Creasey et. al. 2018
! Mass of halo for which the tree is to be grown [h^-1.M_sun]. The mass resolution of the tree and the number of trees to grow. 
! mphalo = 1.0e+12*h0   !halo mass at base of tree, default = 10^14
 mres   = 1.0e+07*h0   !mass resolution, default = 10^8
 ntree  = 10            !number of trees

!------------------------------------------------------------------! 
! now try with the halo mass M(z=0) = [10^11; 10^12.6] M_sun       !
! as in Creasey paper
!------------------------------------------------------------------!
 mphalo_min   = 1.0e+11*h0
 mphalo_max   = 1.0e+13*h0
 mphalo_step  = mphalo_min 
 mphalo       = mphalo_min

  z_max     = 20.0
  z_step    = 0.05
! redshift and halo mass to form GC, use values in Creasey et. al. 2018  
  z_form_GC = 8.65
  M_form_GC = 1.0e+08*h0     
  n_z_bin   = ceiling((z_max - z_form_GC)/z_step)
  allocate(count_GC_array(n_z_bin))
  

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


! Power Spectrum parameters
! (passed in module  Cosmological_Parameters and Power_Spectrum_Parameters)
! pkinfile = 'pk_WDMDove.dat'
 pkinfile= 'pk_Mill.dat' !Tabulated Millennium Simulation linear P(k)
! pkinfile = 'Pk-Creasey18.txt'
 itrans= -1             !indicates use transfer function tabulated in file pkinfile
!  itrans= 1              !indicates use BBKS CDM transfer function with specified Gamma and Omega0
!  itrans= 2            !indicates use Bond & Efstathiou CDM transfer function with specified Gamma and Omega0
! itrans= 3              !indicates use Eisenstein and Hu CDM transfer function with specified Omega0, Omegab and h0
 CMB_T0 = 2.73          !For Eisenstein and Hu CDM transfer function one must specify the CMB temperature


!Set primordial P(k) parameters (ignored if itrans=-1)
 nspec  = 1.0           !primoridial power spectrum spectral index
 dndlnk = 0.0           !allow running spectral index by setting ne.0
 kref   = 1.0           !pivot point for running index
 sigma8 = 0.9         !power spectrum amplitude set regardless of other parameters, default = 0.9

 ierr    = 1            !initial error status to control make_tree()
 nhalomax= 0            !initialise
 nhalo   = 0
 iseed0  = -8635        !random number seed
 iseed   = iseed0

! Set up the array of redshifts at which the tree is to be stored
  write(0,*) 'The redshifts at which the tree will be stored:'
  allocate(wlev(nlev),alev(nlev),ifraglev(nlev))

! Specify output/storage times of the merger tree
  ahalo = 1.0           !expansion factor at base of tree
  zmax  = z_max         !maximum redshift of stored tree, = 4. by default
  do ilev = 1,nlev      !tree levels uniform between z=0 and zmax
     alev(ilev) = 1.0/(1.0+zmax*real(ilev-1)/real(nlev-1))
     dc         = deltcrit(alev(ilev))
     write(0,'(a2,1x,f6.3,1x,a,f6.3)')'z=',(1/alev(ilev)) -1.0,'at which deltcrit=',dc
  end do

!--------------------------------------------!
! Start generating trees                     !
! adding a loop to change M(z=0)             !
!--------------------------------------------!

do while (mphalo .le. mphalo_max)
    !print *, mphalo
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

!---------------------------------------------------------------------------!
!                                                                           !
!    My own code starts from here                                           !
!                                                                           !
!---------------------------------------------------------------------------!
     !--------------------------------------------------------------------------------------------------------------------!
     ! counting number of GCs in each redshift bin                                                                        !
     ! using the conditions of M >= M_form_GC and z >= z_form_GC                                                          !
     ! by walking through the tree and count GC formed whenever hit a node that has M < M_form_GC,                        !
     ! then check its parent mass and redshift, if they satisfy the conditions then count 1 GC formed inside that DM halo !
     ! However, this doesnot cover the node at the end of the tree that satisfy all conditions but have no kid,           !
     ! so there is an additional part below this part                                                                     !
     !--------------------------------------------------------------------------------------------------------------------!
        z = 0. 
        count_GC_array = 0
        kk = 0
        This_Node => MergerTree(1)
        do while (associated(This_Node))
           if (This_Node%mhalo .lt. M_form_GC .and. This_Node%parent%mhalo .ge. M_form_GC) then
               z = 1.0/alev(This_node%parent%jlevel)-1.0
               if (z .ge. z_form_GC) then
                   kk = ceiling((z-z_form_GC)/z_step)
                   count_GC_array(kk) = count_GC_array(kk) + 1                             
                   This_Node => Walk_To_Next_Branch(This_Node)
               endif
           endif
           if (associated(This_Node)) then
               This_Node => Walk_Tree(This_Node)
           endif
        end do
        !print *,count_GC_array
     !--------------------------------------------------------------------------------------------!
     ! This part is for completing the part above, which forgot to count the node at end branches !
     ! The nodes at the end of branches that meet the conditions for forming GC are counted her   !
     ! as they can form GCs inside                                                                !
     !--------------------------------------------------------------------------------------------!
        This_Node => MergerTree(1)
        do while (associated(This_Node))
           if (This_Node%nchild .eq. 0 .and. This_Node%mhalo .ge. M_form_GC) then 
               z = 1.0/alev(This_node%jlevel)-1.0
               if (z .ge. z_form_GC) then
                   kk = ceiling((z-z_form_GC)/z_step)
                   count_GC_array(kk) = count_GC_array(kk) + 1                             
                   This_Node => Walk_To_Next_Branch(This_Node)
               endif
           endif
           if (associated(This_Node)) then
               This_Node => Walk_Tree(This_Node)
           endif
        end do
        !print *, count_GC_array

        !save in file the number count and redshift of halos which form GC
        write(f1,*) i
        mphalo_convert = int(mphalo/(h0*1.e+11))  ! save filename with mass unit [10^11 M_sun]
        write(f2,*) mphalo_convert
        open(unit = 101,file = 'test_results/GCs-tree'//trim(adjustl(f1))//'-'&
                        !& //trim(adjustl(f2))//'e11M_sun.txt',status='replace') 
        do ii = 1, n_z_bin
           write (101,*) z_form_GC + (ii-1)*z_step, mphalo, count_GC_array(ii)
          
        enddo
        close(101)


    enddo

   ! go to the next halo mass
    mphalo = mphalo + mphalo_step
end do
!---------------------------------------------------------------------------!
!                                                                           !
!                    My own code ends here                                  !
!                                                                           !
!---------------------------------------------------------------------------!
 
deallocate(wlev,alev,ifraglev, count_GC_array)

end program tree2
