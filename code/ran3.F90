module ran3_Storage
 !
 ! Module used to store state of ran3() random number generator
 !
 ! Integers
 integer IFF,II,INEXT,INEXTP,jran3,K,MA(55),MJ,MK
 integer ran3rst(63) ! Storage block - not directly used in some functions.
 !
 ! Equivalences
 equivalence (IFF,ran3rst)
end module ran3_Storage

real function ran3(IDUM)
  !
  ! Uses
  use ran3_Storage
  implicit none
  !
  ! Integers
  integer IDUM,MBIG,MSEED,MZ
  !
  ! Floats
  real FAC
  !
  ! Parameters
  parameter (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.0/MBIG)
  if (IDUM.lt.0.or.IFF.eq.0) then
     IFF=1
     MJ=MSEED-iabs(IDUM)
     MJ=mod(MJ,MBIG)	
     if (MJ.lt.MZ) MJ=MJ+MBIG
     MA(55)=MJ
     MK=1
     do jran3=1,54
        II=mod(21*jran3,55)
        MA(II)=MK
        MK=MJ-MK
        if (MK.lt.MZ) MK=MK+MBIG
        MJ=MA(II)
     end do
     do K=1,4
        do jran3=1,55
	   MA(jran3)=MA(jran3)-MA(1+mod(jran3+30,55))
	   if (MA(jran3).lt.MZ) MA(jran3)=MA(jran3)+MBIG
        end do
     end do
     INEXT=0
     INEXTP=31
     IDUM=1
  end if
  INEXT=INEXT+1
  if (INEXT.eq.56) INEXT=1
  INEXTP=INEXTP+1
  if (INEXTP.eq.56) INEXTP=1
  MJ=MA(INEXT)-MA(INEXTP)
  if (MJ.lt.MZ) MJ=MJ+MBIG
  MA(INEXT)=MJ
  ran3=MJ*FAC
  return
end function ran3
