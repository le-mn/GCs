! Produce an index indx to the array arr(n) so that when accessed via
! the index the elements are in ascending numerical order.
! 
subroutine indexxx(n,arr,indx)
  implicit none
  !
  ! Array dimensions
  integer n,NSWITCH
  parameter (NSWITCH=60)
  !
  ! Integers
  integer indx(n)
  !
  ! Floats
  real arr(n)
  ! 
  ! Code
  if (n.gt.NSWITCH) then     ! If n is large use the Numerical Recipes routine
     call indexx(n,arr,indx) ! heapsort.
  else                       ! Simple shell sort which is faster for small n but scales as n^3/2 in worst case
     call indexsh(n,arr,indx)! rather than the n log(n) of heapsort.
  end if
  return
end subroutine indexxx

!
! Indexing or array arr(n) using shell sort.
subroutine indexsh(n,arr,indx)
  !
  ! Uses
  use Numerical_Parameters
  implicit none
  !
  ! Array dimensions
  integer n
  !
  ! Integers
  integer i,indx(n),j,k,l,m,lognb2,nn,tindx
  !
  ! Floats
  real aln2i,arr(n),TINY
  !
  ! Logicals
  logical done3
  !
  ! Parameters
  parameter (TINY=1.0e-5)
  parameter (ALN2I=1.0/LN2) ! 1/ln(2).
  ! 
  ! Code
  forall (i=1:n)
     indx(i)=i
  end forall
  if (n.ge.2) then
     lognb2=int(log(float(n))*ALN2I+TINY)
     m=n
     do nn=1,lognb2
        m=m/2
        k=n-m
        i=1
        do j=1,k
           i=j
           done3=.false.
           do while (.not.done3)
              l=i+m
              if (arr(indx(l)).lt.arr(indx(i))) then
                 tindx=indx(i)
                 indx(i)=indx(l)
                 indx(l)=tindx
                 i=i-m
                 if (i.lt.1) done3=.true.
              else
                 done3=.true.
              end if
           end do
        end do
     end do
  end if
  return
end subroutine indexsh

!
! Numerical Recipes heapsort
subroutine indexx(N,ARRIN,INDX)
  implicit none
  !
  ! Array dimensions
  integer N
  !
  ! Integers
  integer I,INDX(N),INDXT,IR,J,L
  !
  ! Floats
  real ARRIN(N),Q
  !
  ! Code
  forall (J=1:N)
     INDX(J)=J
  end forall
  if (n.le.1) return
  L=N/2+1
  IR=N
  do while (.true.)
     if (L.gt.1) then
        L=L-1
        INDXT=INDX(L)
        Q=ARRIN(INDXT)
     else
        INDXT=INDX(IR)
        Q=ARRIN(INDXT)
        INDX(IR)=INDX(1)
        IR=IR-1
        if (IR.eq.1) then
           INDX(1)=INDXT
           RETURN
        end if
     end if
     I=L
     J=L+L
     do while (J.le.IR)
        if (J.lt.IR) then
           if (ARRIN(INDX(J)).lt.ARRIN(INDX(J+1))) J=J+1
        endif
        if (Q.lt.ARRIN(INDX(J))) then
           INDX(I)=INDX(J)
           I=J
           J=J+J
        else
           J=IR+1
        end if
     end do
     INDX(I)=INDXT
  end do
end subroutine indexx
