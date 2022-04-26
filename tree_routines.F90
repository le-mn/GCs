!% Contains a module of routines for manipulating and utilizing merger trees.

module Tree_Routines
  !% A set of subroutines and functions for manipulating and utilizing the merger trees.
  use Defined_Types
  private
  public :: Walk_Tree,Tree_Get_Hierarchy_Level,Tree_Get_Current_Parent,Node_Galaxy_Active,Walk_Tree_To_Level
  public :: Walk_To_Next_Branch

contains

  function Walk_Tree_To_Level(This_Node,ilev_max) result (Next_Node)
    !% This function provides a mechanism for walking through the branches of the merger tree, going no higher than a given
    !% level. Given a pointer {\tt This\_Node}
    !% to a branch of the tree, it will return the next node that should be visited in the tree. Thus, if {\tt This\_Node} is
    !% initially set to the base of the merger tree and {\tt Walk\_Tree()} is called repeatedly it will walk through every branch
    !% of the merger tree. Once the entire tree has been walked, a {\tt null()} pointer will be returned, indicating that there
    !% are no more branches to walk. Each node will be visited once and once only if the tree is walked in this way.
    integer, intent(in) :: ilev_max
    type (TreeNode), pointer :: This_Node,Next_Node
    Next_Node => This_Node
    if (associated(Next_Node%child).and.Next_Node%jlevel.le.ilev_max) then
       Next_Node => Next_Node%child ! Walk up the tree to next child
    else
       if (associated(Next_Node%sibling)) then
          Next_Node => Next_Node%sibling ! No children, so walk to sibling if one exists.
       else
          do while (.not.associated(Next_Node%sibling).and.associated(Next_Node%parent))
             Next_Node => Next_Node%parent ! No siblings either, so walk back to parents until we find more siblings.
          end do
          if (associated(Next_Node%sibling)) then
             Next_Node => Next_Node%sibling ! Move to the next sibling
          else
             ! Node has no parent - we're back to the base of the tree. Dissassociate the pointer to indicate the tree has been
             ! completely walked.
             Next_Node => null()
          end if
       end if
    end if
    return
  end function Walk_Tree_To_Level

  function Walk_Tree(This_Node) result (Next_Node)
    !% This function provides a mechanism for walking through the branches of the merger tree. Given a pointer {\tt This\_Node}
    !% to a branch of the tree, it will return the next node that should be visited in the tree. Thus, if {\tt This\_Node} is
    !% initially set to the base of the merger tree and {\tt Walk\_Tree()} is called repeatedly it will walk through every branch
    !% of the merger tree. Once the entire tree has been walked, a {\tt null()} pointer will be returned, indicating that there
    !% are no more branches to walk. Each node will be visited once and once only if the tree is walked in this way.
    type (TreeNode), pointer :: This_Node,Next_Node
    Next_Node => This_Node
    if (associated(Next_Node%child)) then
       Next_Node => Next_Node%child ! Walk up the tree to next child
    else
       if (associated(Next_Node%sibling)) then
          Next_Node => Next_Node%sibling ! No children, so walk to sibling if one exists.
       else
          do while (.not.associated(Next_Node%sibling).and.associated(Next_Node%parent))
             Next_Node => Next_Node%parent ! No siblings either, so walk back to parents until we find more siblings.
          end do
          if (associated(Next_Node%sibling)) then
             Next_Node => Next_Node%sibling ! Move to the next sibling
          else
             ! Node has no parent - we're back to the base of the tree. Dissassociate the pointer to indicate the tree has been
             ! completely walked.
             Next_Node => null()
          end if
       end if
    end if
    return
  end function Walk_Tree

  integer function Tree_Get_Hierarchy_Level(This_Node,ilev)
    !% Compute the level of halo {\tt This\_Node} in the merger tree hierarchy (i.e. is it a distinct halo at level {\tt ilev}
    !% [0], or is it a sub-halo [1], or sub-sub-halo [2] etc.) at level {\tt ilev}.
    type (TreeNode), pointer :: This_Node,Parent_Node
    integer, intent(in) :: ilev
    !
    Tree_Get_Hierarchy_Level=0 ! Initialize to zero.
    Parent_Node => This_Node
    do while (Parent_Node%jlevel.gt.ilev) ! Loop over nodes until current level is reached.
       ! Increase level if not most massive child of parent (i.e. is a sub-halo of the parent).
       if (.not.associated(Parent_Node,Parent_Node%parent%child)) Tree_Get_Hierarchy_Level=Tree_Get_Hierarchy_Level+1
       Parent_Node => Parent_Node%parent
    end do
    return
  end function Tree_Get_Hierarchy_Level

!----------------------------------------------------------!
! this is the function to walk to the node in next branch  !
!----------------------------------------------------------!
  function Walk_To_Next_Branch(This_Node) result (Next_Node)
    !% This function provides a mechanism for walking through the branches of the merger tree. Given a pointer {\tt This\_Node}
    !% to a branch of the tree, it will return the 1st node in next branch of that tree

    type (TreeNode), pointer :: This_Node, Next_Node
    Next_Node => This_Node
    if (associated(Next_Node%sibling)) then
          Next_Node => Next_Node%sibling ! walk to sibling if one exists.
       else
          do while (.not.associated(Next_Node%sibling).and.associated(Next_Node%parent))
             Next_Node => Next_Node%parent ! No siblings either, so walk back to parents until we find more siblings.
          end do
          if (associated(Next_Node%sibling)) then
             Next_Node => Next_Node%sibling ! Move to the next sibling
          else
             ! Node has no parent - we're back to the base of the tree. Dissassociate the pointer to indicate the tree has been
             ! completely walked.
             Next_Node => null()
          end if
    end if
    return
  end function Walk_To_Next_Branch


end module Tree_Routines


