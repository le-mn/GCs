module Defined_Types
  use Kind_Numbers

  type TreeNode
     real :: mhalo
     integer :: jlevel,nchild,index
     type (TreeNode), pointer :: child,parent,formation,sibling
  end type TreeNode

end module Defined_Types
