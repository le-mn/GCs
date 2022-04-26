module Tree_Memory_Arrays
  !
  ! Uses
  use Defined_Types
  !
  ! Other types
  type (TreeNode), allocatable, target :: MergerTree_Aux(:)
end module Tree_Memory_Arrays

module Tree_Memory_Arrays_Passable
  use Defined_Types
  type (TreeNode), pointer :: MergerTree(:)
contains

  pure integer function Tree_Index(Node)
    implicit none
    !
    ! * Arguments *
    !
    ! Other types
    type (TreeNode), pointer :: Node
    !
    ! Code
    Tree_Index = Node%index
    return
  end function Tree_Index

  integer function Tree_Formation_Index(i_Halo)
    implicit none
    !
    ! * Arguments *
    !
    ! Integers 
    integer, intent(in) :: i_Halo
    !
    ! * Other variables *
    !
    !
    ! Code
#ifdef TRAPS
    if (.not.associated(MergerTree(i_Halo)%formation)) then
       write (0,*) 'Tree_Formation_Index(): FATAL - formation pointer is not associated'
       write (0,*) '                        i_Halo = ',i_Halo
       stop
    end if
#endif

    Tree_Formation_Index=MergerTree(i_Halo)%formation%index
    return
  end function Tree_Formation_Index
end module Tree_Memory_Arrays_Passable
