! Subroutine to build a merger tree, based on binary splitting.
! 
! This version has the memory allocation for the workspace arrays
! wlev,ifraglev,ipar,isib,ichild,mtr,jindex
! be done by the calling program, rather than by the subroutine itself.
! 
! A 2nd re-ordering of the indices of tree nodes is done, so that the
! children of any node are ordered in decreasing mass.
! 
! INPUTS:
! -------
! m0 = final halo mass at a=a0 (a0 is arbitary)
! mmin = min progenitor mass that follow
! alev(ilev),(ilev=1,nlev) = array of values of a at which output 
! details of tree, with alev(1)=a0
! iseed = seed for random number generator
! nfragmax = max size of arrays for output tree structure (see below)
! 
! 
! INPUT SUBROUTINES/FUNCTIONS: (declare EXTERNAL in calling program)
! ----------------------------
! Any additional parameters needed by these subroutines (e.g. Omega, 
! nspec,eps) assumed to be communicated from calling program by 
! modules.
! 
! split(m,mmin,sigma,iseed,dwmax,dw,nprog,mprog) 
! - subroutine does single binary split
! dwmax (input) is max step in w=deltc(a), dw (output) is actual step
! nprog(=0,1,2) is no. of progenitors, mprog(2) is array of their masses
! mprog(1) >= mprog(2)
! 
! sigma(m,dlsigdlm) 
! - function calculates linear sigma(m) and dln(sigma)/dln(m) at a=1
! 
! deltc(a)
! - function calcs threshold linear overdensity for collapse at epoch a
! 
! 
! OUTPUTS:
! --------
! N.B. calling program must check ierr on return!
! ierr = 0 (worked OK) 
! ierr = 1 (need to increase nfragmax) 
! NB  nfragtot = total no of fragments in tree
! ierr = 2 (need to rerun now that inodemax has been increased 
!           but no need to increase nfragmax)
! 
! ilev=1,nlev:
! nfraglev(ilev) = no. of fragments at time level ilev
! jpfrag(ilev)   = index jfrag for first fragment at time 
! level ilev
! 
! jfrag=1,nfragtot:
! jlevel(jfrag)            = level ilev for fragment jfrag
! mfrag(jfrag)             = mass for fragment jfrag
! MergerTree(jfrag)%parent = jfrag for parent
! nchild(jfrag) = no. of children for fragment jfrag
! MergerTree(jfrag)%child = pointer to first child (largest mass)
! children have indices jfrag = Tree_Index(MergerTree%child) + i, i=0,nchild-1, in order of 
! decreasing mass
! 
! WORKSPACE ARRAYS:
! -----------------
! wlev(ilev),ifraglev(ilev), ilev = 1,nlev
! mtr(ifrag),ipar(ifrag),isib(ifrag),ichild(ifrag),jindex(ifrag), 
! ifrag = 1,nfragmax
! 
!
! 
! WHAT IT DOES:
! -------------
! Construct a merger tree, based on a prescription where each halo splits
! in 2 pieces at each merger event (i.e. going "up" tree, back in time).
! the prescription for the binary split is encoded in the subroutine 
! split(m,mmin,iseed,dwmax,dw,nprog,mprog,output). In fact, split() may 
! return nprog=0,1 or 2, the former indicating 0 or 1 progenitors above 
! the resolution mass mmin. Using split(), we can construct the entire
! ELEMENTARY TREE specifying every binary merger event.
! 
! However, the elementary tree contains more information than we usually
! want, so we use it to construct a REDUCED TREE, which provides 
! information about halo progenitors at a user-specified set of expansion
! factors alev(ilev). The information recorded for each time-level ilev
! includes the mass of each progenitor, a pointer to its parent (on the 
! next time level down i.e. later) and pointers to all of its children
! (on the next time level up i.e. earlier). The REDUCED TREE is not binary,
! each parent can have many children.
! 
! We construct the REDUCED TREE by following the ELEMENTARY TREE one branch
! at a time, storing the necessary information as we go. In fact, there
! are 2 stages: the process of its construction defines a natural array
! structure for the reduced tree, but this is not the most convenient. 
! Therefore, the information is re-arranged into array structures that
! allow easy movement up, down and sideways (i.e. at constant time-level)
! in tree.
! 
! N.B.: the subroutine enforces alev(1)=a0, i.e. ilev=1 corresponds to
! the final halo of mass m0 at a=a0
! 
!
! ARRAYS USED TO STORE TREES:
! ---------------------------
! 
! ELEMENTARY (BINARY) TREE (temporary):
! -------------------------------------
! inode = 1,inodemax:
! wnode(inode) = value of w (=deltc(a)) at that node
! mlnode(inode),mrnode(inode) = masses of "left" and "right" fragments
! at each point where binary tree splits in 2
! (moving up the tree, the 2 fragments are considered to appear at the 
! value of w given by wnode)
! lnode(inode) = .true. if node is at a level of wlev(), .false. otherwise
! 
! INTERMEDIATE REDUCED TREE (temporary):
! --------------------------------------
! ilev = 1,nlev:
! wlev(ilev)     = w at that time-level
! ifraglev(ilev) = ifrag for "active" (most recently processed) fragment
! at that ilev
! 
! ifrag = 1,nfragmax:
! mtr(ifrag)    = mass for that fragment
! ipar(ifrag)   = ifrag for parent of that fragment
! isib(ifrag)   = ifrag for sibling (= -1 if no sibling)
! ichild(ifrag) = ifrag for first child (= -1 if no children)
! 
! jindex(ifrag) = new index jfrag (for final tree) for that fragment
! 
! 1st level (ilev=1), 1st fragment (ifrag=1) is final halo m=m0 at a=a0
! 
! in re-arranging tree for 2nd time, to put children of any node in order
! of decreasing mass, isib(), ichild() and mtr() are re-used as temporary 
! arrays:
! isib() stores nchild() as function of new fragment index
! ichild() stores index of child node as function of new fragment index
! mtr() stores mfrag() as function of new fragment index
! jindex(kfrag) is old index jfrag as function of new index kfrag
!
! child_ref() is used to temporarily store the index of the nodes child
! 
! FINAL REDUCED TREE (output):
! ----------------------------
! as explained above under OUTPUTS
! 1st level (ilev=1), 1st fragment (jfrag=1) is final halo m=m0 at a=a0
! 
! an example of the indexing of the output final tree structure is:
! 
! .
! .
! .
! 3rd level (ilev=3):   jfrag = 5,6,7,8,9,10
! 2nd level (ilev=2):   jfrag = 2,3,4
! 1st level (ilev=1):   jfrag = 1
! 
! where: (-1 is null pointer)
! jlevel(1)=1,MergerTree(1)%parent=-1,nchild(1)=3,Tree_Index(MergerTree(1)%child)=2
! jlevel(2)=2,MergerTree(2)%parent=1,nchild(2)=2,Tree_Index(MergerTree(2)%child)=5
! jlevel(3)=2,MergerTree(3)%parent=1,nchild(3)=4,Tree_Index(MergerTree(3)%child)=7
! jlevel(4)=2,MergerTree(4)%parent=1,nchild(4)=0,Tree_Index(MergerTree(4)%child)=-1
! jlevel(5)=3,MergerTree(5)%parent=2,...
! .
! .
! i.e. 
! (1) has children (2,3,4)
! (2) has children (5,6)
! (3) has children (7,8,9,10)
! (4) has 0 children
! 
! nfraglev(1)=1,jpfraglev(1)=1
! nfraglev(2)=3,jpfraglev(1)=2
! nfraglev(3)=6,jpfraglev(3)=5
! .
! .
! 
!
! UNITS:
! ------
! units for mass m can be anything, but should match what is assumed by
! sigma(m,dlsigdlm) routine
! 
!

module Make_Tree_Arrays
  !
  ! Floats
  real, allocatable :: ml(:),mr(:),wnode(:)  
  !
  ! Logicals
  logical(kind=1), allocatable :: lnode(:)
end module Make_Tree_Arrays

module Make_Tree_Module
contains

  function Next_Sibling(Child_Node,Parent_Node,Sibs_Left) result (Sibling)
    use Defined_Types
    logical Sibs_Left
    type (TreeNode), pointer :: Child_Node,Parent_Node,Sibling
    if (associated(Child_Node,Parent_Node)) then
       if (associated(Parent_Node%child)) then
          Sibling => Parent_Node%child
       else
          stop 'Next_Sibling(): FATAL - parent has no children'
       endif
    else
       if (associated(Child_Node%sibling)) then
          Sibling => Child_Node%sibling
       else
          stop 'Next_Sibling(): FATAL - no more siblings remain'
       endif
    endif
    if (associated(Sibling%sibling)) then
       Sibs_Left=.true.
    else
       Sibs_Left=.false.
    endif
    return
  end function Next_Sibling

  subroutine Associate_Siblings(This_Node)
    ! Associate sibling pointers for children of This_Node
    !
    ! Uses
    use Defined_Types
    use Tree_Memory_Arrays_Passable
    implicit none
    ! * Arguments *
    ! Other types
    type (TreeNode) :: This_Node
    ! * Other variables *
    ! Integers
    integer Child_Index,ifrag
    ! Code
    if (This_Node%nchild.gt.1) then
       Child_Index=Tree_Index(This_Node%child)
       do ifrag=Child_Index,Child_Index+This_Node%nchild-2
          MergerTree(ifrag)%sibling => MergerTree(ifrag+1)
       end do
    endif
    return
  end subroutine Associate_Siblings

  subroutine Build_Sibling_Pointers(nfragtot)
    use Defined_Types
    use Tree_Memory_Arrays_Passable
    implicit none
    integer :: nfragtot
    ! * Other variables *
    integer :: i
    !
    ! Code
    do i = 1, nfragtot, 1
       call Associate_Siblings(MergerTree(i))
    end do
       
  end subroutine Build_Sibling_Pointers
end module Make_Tree_Module

subroutine make_tree(m0,a0,mmin,alev,nlev,iseed,split,sigma,deltc,nfragmax,ierr,nfragtot,nfraglev,jpfrag&
     &,wlev,ifraglev)
  !
  ! Uses
  use Defined_Types
  use Make_Tree_Arrays
  use Make_Tree_Module
#ifdef DEBUG
  use Run_Statistics
#endif
  use Tree_Memory_Arrays_Passable
  implicit none
  !
  ! Array dimensions
  integer NCHMAX,nfragmax,nlev
  parameter (NCHMAX=1e5)
  !
  ! Integers
  integer alloc_err,ichild(nfragmax),ierr,ifrag,ifragc,ifraglev(nlev),ifrag_prev,ilev,ilevwk,indxch(NCHMAX),inode,inodemax&
       &,ipar(nfragmax),ipar_prev,iseed,isib(nfragmax),iw,jchild1,jfrag,jfragc,jfragp,jindex(nfragmax),child_ref(nfragmax)&
       &,jpfrag(nlev),kchild,kfrag,kfragp,MAXNODES,nch,nchild1,nfraglev(nlev)&
       &,nfragtot,nprog
  !
  ! Floats
  real a0,alev(nlev),dw,dwmax,m,m0,mbytes,mmin,mprog(2),mtr(nfragmax),w,wfin,wlev(nlev)
  !
  ! Parameters
  parameter (MAXNODES=1e7) ! Absolute limit on number of nodes.
  !
  ! Functions
  real sigma,deltc
  !
  ! Externals
  external split,sigma,deltc
  !
  ! Saves
  save inodemax
  !
  ! Data
  data inodemax /0/ ! Init value for max no of nodes
  ! 
  ! Code
#ifdef DEBUG
  write (0,*) 'make_tree(): DEBUG - starting'
#endif
  !
  ! Create storage arrays
  if (inodemax.eq.0) then
     inodemax=2000
     allocate(wnode(inodemax),stat=alloc_err)
if (alloc_err.ne.0) stop 'make_tree(): FATAL - failed to allocate memory [wnode]'
     allocate(ml(inodemax),stat=alloc_err)
if (alloc_err.ne.0) stop 'make_tree(): FATAL - failed to allocate memory [ml]'
     allocate(mr(inodemax),stat=alloc_err)
if (alloc_err.ne.0) stop 'make_tree(): FATAL - failed to allocate memory [mr]'
     allocate(lnode(inodemax),stat=alloc_err)
if (alloc_err.ne.0) stop 'make_tree(): FATAL - failed to allocate memory [lnode]'
     mbytes=real(inodemax*(4*3+1))/real(1024**2)
!     write(0,'(a,f7.3,a)') 'make_tree(): internal storage = ',mbytes,' Mbytes' 
  end if

!write(0,*)'created storage arrays'
  ! 
  ! Reset output arrays to 0 or -1 (=pointer to nothing)
  ! 
  nfragtot         =0
  nfraglev(:)      =0
  jpfrag           =-1
  MergerTree%mhalo =0.0
  MergerTree%jlevel=0
  MergerTree%nchild=0
  do ifrag=1,nfragmax
     MergerTree(ifrag)%parent => null() ! Disassociate all parent pointers initially.
     MergerTree(ifrag)%child => null() ! Disassociate all child pointers initially.
     MergerTree(ifrag)%sibling => null() ! Disassociate all sibling pointers initially.
  end do
  !
#ifdef DEBUG
  write (0,*) 'make_tree(): DEBUG - creating array of values of w=deltc(a)'
#endif
  !
  ! Create array of values of w=deltc(a) at which output tree info. 
  alev(1)=a0 ! Enforce this!
  do ilev=1,nlev
     wlev(ilev)=deltc(alev(ilev))
  enddo
  inode = 1
#ifdef DEBUG
  write (0,*) 'mke_tree(): DEBUG - done creating array of values of w=deltc(a)'
#endif
  !
  ! Make binary tree.
  ! 
  ! ml and mr are masses of "left" and "right" forks.
  ! ml=-1 or mr=-1 indicates that fork empty (m<mmin)
  ! branch ends because either (a) m<mmin or (b) w>wfin
  ! 
  wfin=wlev(nlev) ! Earliest output time.
  ! 
  ! Initialize array of "active" fragments.
  ifraglev(:)=-1
  ! 
  ! Base (root) of tree (a=a0)
  inode   = 1
  ml(1)   = m0
  mr(1)   = -1 ! No right fork
  m       = m0
  wnode(1)= wlev(1)
  w       = wlev(1)
  ilev    = 2 ! Next output level up tree
  ifraglev(1)= 1 ! Active fragment at ilev=1
  ! 
  ifrag      = 1
  mtr(1)     = m0
  ! 
  ipar(1)    = -1
  isib(1)    = -1
  ichild(1)  = -1
  ! 
#ifdef DEBUG
  max_inode = 0 ! Count max no of nodes
#endif
  ! 
  do while(ml(inode).gt.0.0) ! ml=-1 signals end of branch.
     dwmax = wlev(ilev)-w
     ! Take step up tree (back in time) by binary split.
     call split(m,w,mmin,sigma,iseed,dwmax,dw,nprog,mprog) 
     w=w+dw
     !
     select case (nprog)
     case (2) ! Create new node.
        inode=inode+1
#ifdef DEBUG
        max_inode=max(inode,max_inode)
#endif
        ! If inode > inodemax, exit this loop and restart with larger inodemax.
        if (inode.gt.inodemax) then
#ifdef INFO
           write(0,*) 'make_tree(): INFO - increasing size inodemax'
#endif
           deallocate(wnode)
           deallocate(ml)
           deallocate(mr)
           deallocate(lnode)
           inodemax=2*inodemax
           if (inodemax.le.MAXNODES) then
              allocate(wnode(inodemax),stat=alloc_err)
if (alloc_err.ne.0) stop 'make_tree(): FATAL - failed to allocate memory&
                   & [wnode]'
              allocate(ml(inodemax),stat=alloc_err)
if (alloc_err.ne.0) stop 'make_tree(): FATAL - failed to allocate memory&
                   & [ml]'
              allocate(mr(inodemax),stat=alloc_err)
if (alloc_err.ne.0) stop 'make_tree(): FATAL - failed to allocate memory&
                   & [mr]'
              allocate(lnode(inodemax),stat=alloc_err)
if (alloc_err.ne.0) stop 'make_tree(): FATAL - failed to allocate memory&
                   & [lnode]'
              mbytes=real(inodemax*(4*3+1))/1024**2
!              write(*,'(a,f7.3,a)') 'internal halo tree storage= ',mbytes,' Mbytes'
              ierr=2 
              return
           else
              stop 'make_tree(): FATAL - inodemax>MAXNODES'
           end if
        end if
        wnode(inode)=w 
        lnode(inode)=.false.
        ml(inode)   =mprog(1)
        mr(inode)   =mprog(2)
        m           =ml(inode) ! Follow left branch.
     case (1) ! Update mass on this branch.
        m=mprog(1)
     case (0) ! End of branch
        m=0
        ml(inode)=-1
     end select
     ! 
     ! Store info for reduced tree, if w=wlev and nprog>0 from latest split
     ! If nprog=1, store inf for that fragment
     ! If nprog=2, store inf for both
     if(w.eq.wlev(ilev).and.nprog.gt.0) then          
        if((ifrag+2).gt.nfragmax) then ! Up to 2 more fragments.
           ierr=1 ! Signals failure of routine to calling program.
           return
        end if
        ! 
        ! 1st (or only) fragment from latest split.
        ifrag        =ifrag+1
        mtr(ifrag)   =mprog(1)
        ipar(ifrag)  =ifraglev(ilev-1) ! Parent on previous level.
        isib(ifrag)  =-1    ! Initialize sibling pointer.
        ichild(ifrag)=-1  ! Initialize child pointer.
        ! 
        ! Find out whether previous fragment on same ilev has same or different
        ! parent.
        ifrag_prev=ifraglev(ilev) ! Previous frag at same ilev.
        if(ifrag_prev.gt.0) then ! There is active frag on that ilev.
           ipar_prev=ipar(ifrag_prev) ! Parent for prev frag same ilev.
           if(ipar_prev.eq.ipar(ifrag)) then ! Same parent.
              isib(ifrag_prev)=ifrag ! Current frag is sibling of prev frag.
           else ! Different parent 
              ichild(ipar(ifrag))=ifrag ! So curr frag is 1st child.
           end if
        else ! This is 1st frag on this ilev.
           ichild(ipar(ifrag))=ifrag ! So curr frag is 1st child.
        end if
        ! 
        ifraglev(ilev)=ifrag ! Update index for "active" fragment.
        ! 
        ! 2nd fragment from latest split (if non-zero).
        if (nprog.eq.2) then
           ifrag        =ifrag+1
           mtr(ifrag)   =mprog(2)
           ipar(ifrag)  =ifraglev(ilev-1) ! Parent on previous level.
           isib(ifrag)  =-1 ! Initialize sibling pointer.
           ichild(ifrag)=-1 ! Initialize child pointer.
           isib(ifrag-1)=ifrag ! Sibling of previous frag.
           lnode(inode) =.true. ! This node occurs at a level of wlev.
           ! If this is top level, must increment ifraglev() now.
           if (ilev.eq.nlev) ifraglev(ilev) = ifrag
        end if
        ! 
        ! Next ilev up tree has changed.
        if (ilev.lt.nlev) ilev=ilev+1
     end if
     ! 
     ! 
     ! Terminate branch if reach earliest output epoch.
     if(w.ge.wfin) then     
        if(nprog.eq.2) then ! Current node has w=wfin.
           ml(inode)=-1
           mr(inode)=-1
        else ! Current node has w<wfin.
           ml(inode)=-1
        end if
     end if
     ! 
     ! If this branch finished, do other branch or go back down to previous node.
     do while(ml(inode).lt.0.0.and.inode.gt.1)
        if (mr(inode).gt.0.0) then ! Do other branch.
           ml(inode)=mr(inode) ! Swap left and right.
           mr(inode)=-1
           w        =wnode(inode)
           m        =ml(inode)
           ! If next wlev up tree has changed, find new ilev.
           if (w.lt.wlev(ilev-1).or.w.ge.wlev(ilev)) then
              ! Locate finds iw such that wlev(iw) < w <= wlev(iw+1).
              call locate(wlev,nlev,w,iw) 
              ilev=iw+1 ! Level above w.
              ! Require w < wlev(ilev).
              if (wlev(ilev).eq.w) ilev=ilev+1
           end if
           ! 
           ! If this node is at a level of wlev, must update index of active fragment
           ! at the level of the node, which is ilev-1.
           if (lnode(inode)) ifraglev(ilev-1)=ifraglev(ilev-1)+1
           ! 
        else ! Down to previous node
           inode    =inode-1
           ml(inode)=-1
        end if
     end do
  end do

  nfragtot=ifrag ! Total number of fragments in tree.
  ! 
  ! Now re-arrange reduced tree.
  ! 
  ! For each level ilevwk separately, move up through tree and process
  ! all fragments at that level, writing arrays for re-ordered tree
  ! first sweep assigns new index jfrag in terms of old index ifrag.
  ! 
  jfrag=0
  do ilevwk=1,nlev
     ! Start from root
     ilev =1
     ifrag=1
     do while(ifrag.gt.0)
        ! Move up tree to first ("leftmost") fragment in that subtree on working 
        ! level ilevwk, starting from current node.
        do while (ilev.lt.ilevwk.and.ichild(ifrag).gt.0)
           ifrag=ichild(ifrag)
           ilev =ilev+1
        end do
        if (ilev.eq.ilevwk) then ! Process siblings at that ilev.
           if (jpfrag(ilevwk).lt.0) then ! 1st frag on that level.
              jpfrag(ilevwk)=jfrag+1
           end if
           do while (ifrag.gt.0) ! Move through siblings at that ilev.
              jfrag                   =jfrag+1
              jindex(ifrag)           =jfrag
              nfraglev(ilevwk)        =nfraglev(ilevwk)+1
              MergerTree(jfrag)%mhalo =mtr(ifrag)
              MergerTree(jfrag)%jlevel=ilevwk
              ! New index jfrag for parent.
              if (ilev.gt.1) then
                 MergerTree(jfrag)%parent => MergerTree(jindex(ipar(ifrag)))
              else
                 MergerTree(jfrag)%parent => null()
              end if
              ifrag_prev=ifrag
              ifrag     =isib(ifrag) ! isib=-1 if no sibling.
           end do
           ifrag=ifrag_prev ! Reset ifrag to last sibling.
           if (ilev.gt.1) then ! Go back to parent level.
              ilev =ilev-1 
              ifrag=ipar(ifrag)
           end if
        end if
        ! 
        ! Go to next sibling if there is one, or go down tree until reach level
        ! where there is one.
        do while (isib(ifrag).lt.0.and.ilev.gt.1)
           ilev =ilev-1
           ifrag=ipar(ifrag)
        end do
        if (ilev.gt.1) then
           ifrag=isib(ifrag)
        else ! Back at root.
           ifrag=-1 ! To signal finished with tree search.
        end if
     end do
  end do

  ! 
  ! Have now found new index jfrag for every element in tree
  ! go through tree again, for each fragment counting its children
  ! and finding index of first child.
#ifdef DEBUG
  write (0,*) 'make_tree(): DEBUG - counting children'
#endif
  do ifrag =1,nfragtot
     jfrag =jindex(ifrag)
     ifragc=ichild(ifrag) ! Old index of 1st child.
     if (ifragc.gt.0) then ! There are children.
        child_ref(jfrag)=jindex(ifragc) ! New index of 1st child.

        MergerTree(jfrag)%child => MergerTree(jindex(ifragc))

        ! Count children.
        nch=0
        do while (ifragc.gt.0)
           nch   =nch+1
           ifragc=isib(ifragc)
        end do
        MergerTree(jfrag)%nchild=nch
     end if
  end do
  ! 
  ! Re-order fragments in each level of tree so that children of any node
  ! are indexed in order of decreasing mass. note that the re-indexing 
  ! only swaps elements on the same level ilev.
  ! 
  ! Note that isib(), ichild() and mtr() are re-used as temporary arrays:
  ! isib() stores nchild() as function of new fragment index
  ! ichild() stores index of child node as function of new fragment index
  ! mtr() stores mfrag() as function of new fragment index
  ! jindex(kfrag) is old index jfrag as function of new index kfrag.
  ! 

  kfrag     = 1
  jindex(1) = 1
  mtr(1)    = MergerTree(1)%mhalo
  do ilev=1,nlev-1
     do kfragp      =jpfrag(ilev),jpfrag(ilev)+nfraglev(ilev)-1
        jfragp      =jindex(kfragp)
        nchild1     =MergerTree(jfragp)%nchild
        isib(kfragp)=nchild1 ! Store nchild(kfragp) in isib(kfragp).
        jchild1     =child_ref(jfragp)
        ! Sort children by mass.
        select case (nchild1)
        case (1)
           indxch(1)=1
        case (2:)
           if (nchild1.gt.NCHMAX) then
              write (0,*) 'make_tree(): FATAL - increase NCHMAX'
              write (0,*) '             nchild1 = ',nchild1
              stop
           endif
           call indexxx(nchild1,MergerTree(jchild1:jchild1+nchild1-1)%mhalo,indxch)
        end select
        ! Assign new index kfrag to each child in order of decreasing mass.
        if (nchild1.ge.1) then
           ichild(kfragp)=kfrag+1 ! Store index of child node of node kfragp in ichild(kfragp)
        else
           ichild(kfragp)=-1
        end if
        do kchild       =nchild1,1,-1
           jfragc       =indxch(kchild)+jchild1-1
           kfrag        =kfrag+1
           jindex(kfrag)=jfragc
           MergerTree(kfrag)%parent => MergerTree(kfragp) ! Reset parent as fn of new index kfrag.
           mtr(kfrag)=MergerTree(jfragc)%mhalo ! Store mfrag(kfrag) in mtr(kfrag).
        end do
     end do
  end do
  ! 
  ! For top level, have to set the new values of child index,nchild()
  ! separately.
  forall (kfrag=jpfrag(nlev):jpfrag(nlev)+nfraglev(nlev)-1)
     isib(kfrag)  =0 ! nchild(kfrag)=0
     ichild(kfrag)=-1 ! index of child of node kfrag=-1
  end forall
  ! 
  ! Now copy elements of nchild(), mfrag() from
  ! temporary arrays to final arrays.
  MergerTree(1:nfragtot)%mhalo=mtr(1:nfragtot)
  child_ref(1:nfragtot)       =ichild(1:nfragtot)
  do ifrag=1,nfragtot
     if (ichild(ifrag).gt.0) then
        MergerTree(ifrag)%child => MergerTree(ichild(ifrag))
     else
        MergerTree(ifrag)%child => null() ! Disassociate pointer if no children.
     endif
  end do
  MergerTree(1:nfragtot)%nchild=isib(1:nfragtot)
  !
  ! Build sibling pointers
  call Build_Sibling_Pointers(nfragtot)
  ! 
  ierr=0 ! Routine successful
  ! 
#ifdef DEBUG
  write (0,*) 'make_tree(): DEBUG - done'
#endif
  return
end subroutine make_tree

subroutine Make_Tree_Deallocate
  use Make_Tree_Arrays
  if (allocated(wnode)) deallocate(wnode)
  if (allocated(ml))    deallocate(ml)
  if (allocated(mr))    deallocate(mr)
  if (allocated(lnode)) deallocate(lnode)
  return
end subroutine Make_Tree_Deallocate
