! Allocate or alter the allocation of memory to halo tree arrays
! and galaxy property arrays
! 
! On first call nhalomax=0 and the array size is set to a estimate of
! based on the parent halo mass, the mass resolution and the number of 
! time steps.
!
! On subsequent calls no change is made if nhalo falls in the range
!      nhalomax/SHRINK_IF > nhalo > nhalomax
! with the second condition being signalled by ierr=1 and not the
! value of nhalo.
!
!  If nhalo is too large nhalomax  is increased by a factor GROW 
!
!  If nhalo is small compared to nhalomax for TWO consecutive calls
!  then nhalomax is reduced by the factor SHRINK
! 
!  nhalomax is not reduced below the value MINHALOS as once the arrays
!  are this small there is nothing really to be gained by shrinking them
!  further.
!
!  If nhalomax gets as large as MAXHALOS then the program gives up.
!
!    
subroutine memory(nhalo,nhalomax,ierr,nlev,mphalo,mres)
  use Tree_Memory_Arrays
  implicit none
  integer nhalo,nhalomax,ierr,nlev,nhalo_prev,ifirst_local,n_allocate
  real mphalo,mres
  !real mbytes
  real GROW,SHRINK,SHRINK_IF
  parameter (GROW=1.414,SHRINK=2.0,SHRINK_IF=4.0)
  integer MAXHALOS,MINHALOS,alloc_err,nhalomax_min
  parameter(MAXHALOS=1e+07,MINHALOS=1e+04)
  save nhalo_prev,ifirst_local
  data ifirst_local /0/
  !
  ! For N-body runs, the minimum value that nhalomax can take must be nhalo (as here we know in advance the size of the tree). So,
  ! we specify this in advance.
     nhalomax_min=1
  if (ierr.eq.2) then
     ! Do nothing as it is the arrays internal to make_tree() that
     ! needed resizing and not these main arrays      
  else if (ierr.eq.1.or.nhalo.gt.nhalomax.or.nhalomax.eq.0) then
     if (nhalomax.ne.0) then
        ! Release currently allocated memory
        call Memory_Deallocate
        nhalomax=max(int(GROW*float(nhalomax)),nhalomax_min)
     else
        ! If first call to the subroutine make an initial estimate
        ! of the memory requirement. This is very rough. 
        nhalomax=max(min(MAXHALOS,max(int(0.015*(mphalo/mres))*nlev,MINHALOS)),nhalomax_min)
     end if
     if(nhalomax.le.MAXHALOS) then
        allocate(MergerTree_Aux(nhalomax),stat=alloc_err)
        if (alloc_err.ne.0) stop 'Memory(): FATAL - Failed to allocate memory [MergerTree_Aux]'
        n_allocate=1
        ierr=0
        !mbytes=real((sizeof(Node)*nhalomax))/real(1024**2)
        !write(0,'(a,f7.1,a)') 'main storage= ',mbytes,' Mbytes'
     else
        stop 'memory(): nhalomax>MAXHALOS'
     end if
  else if (max(nhalo,nhalo_prev)*SHRINK_IF.lt.nhalomax.and.nhalomax.gt.MINHALOS) then 
     ! If for the 2nd consecutive time nhalo is much smaller
     ! than nhalomax then reduce nhalomax by a factor SHRINK 
     call Memory_Deallocate
     nhalomax=max(int(float(nhalomax)/SHRINK),nhalomax_min)
     write(0,*) 'shrink: nhalomax=',nhalomax
     allocate(MergerTree_Aux(nhalomax),stat=alloc_err)
     if (alloc_err.ne.0) stop 'Memory(): FATAL - Failed to allocate memory [MergerTree_Aux]'
     n_allocate=1
     if (alloc_err.ne.0) stop 'Memory(): FATAL - Failed to allocate memory [NBodyMergerTree_Aux]'
     !mbytes=real(sizeof(Node)*nhalomax)/real(1024**2)
     !write(0,'(a,f7.1,a)') 'main storage= ',mbytes,' Mbytes'
  end if
  nhalo_prev=nhalo ! save this value for the next call
  return
end subroutine memory

subroutine Memory_Deallocate
  !
  ! Uses
  use Tree_Memory_Arrays
  !
  ! Code
  !
  ! Release currently allocated memory.
  if (allocated(MergerTree_Aux)) deallocate(MergerTree_Aux)
  return
end subroutine Memory_Deallocate
