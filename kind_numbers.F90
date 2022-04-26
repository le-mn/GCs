module kind_numbers
!
! Define parameters with the 'kind' numbers for reals and integers -
! kind numbers are sometimes not the same as *4 *8 etc - they vary
! between compilers!
!
! Reals:
! ------
!
! Kind for real*4
  INTEGER, PARAMETER :: real4byte = SELECTED_REAL_KIND(6,37)
! Kind for real*8
  INTEGER, PARAMETER :: real8byte = SELECTED_REAL_KIND(15,307)

! Integers:
! ---------
!
! Kind for integer*1
! (note that this will return kind for integer*2 or integer*4
! if *1 / *2 are not available on this system)
  INTEGER, PARAMETER :: int1byte = SELECTED_INT_KIND(2)
! Kind for integer*2
  INTEGER, PARAMETER :: int2byte = SELECTED_INT_KIND(4)
! Kind for integer*4
  INTEGER, PARAMETER :: int4byte = SELECTED_INT_KIND(9)
! Kind for integer*8
  INTEGER, PARAMETER :: int8byte = SELECTED_INT_KIND(18)

CONTAINS

  SUBROUTINE check_kind_numbers(need_i1, need_i2, need_i4, need_i8)
!
! This checks which kind numbers are available for integers
!
    IMPLICIT NONE
    LOGICAL, OPTIONAL, INTENT(IN) :: need_i1, need_i2, need_i4, need_i8

    INTEGER(KIND=int1byte) :: i1
    INTEGER(KIND=int2byte) :: i2
    INTEGER(KIND=int4byte) :: i4
    INTEGER(KIND=int8byte) :: i8

    INTEGER :: isize

    IF(BIT_SIZE(i1).NE.8)THEN
       WRITE(*,*)'check_kind_numbers(): WARNING - INTEGER*1 not available'
       IF(PRESENT(need_i1))THEN
          IF(need_i1)STOP 'check_kind_numbers(): FATAL - Required datatype INTEGER*1 is not available'
       END IF
       isize = BIT_SIZE(i1) / 8
       WRITE(*,*)'Using INTEGER*',isize,' instead'
    END IF

    IF(BIT_SIZE(i2).NE.16)THEN
       WRITE(*,*)'check_kind_numbers(): WARNING - INTEGER*2 not available'
       IF(PRESENT(need_i2))THEN
          IF(need_i2)STOP 'check_kind_numbers(): FATAL - Required datatype INTEGER*2 is not available'
       END IF
       isize = BIT_SIZE(i2) / 8
       WRITE(*,*)'Using INTEGER*',isize,' instead'
    END IF

    IF(BIT_SIZE(i4).NE.32)THEN
       WRITE(*,*)'check_kind_numbers(): WARNING - INTEGER*4 not available'
       IF(PRESENT(need_i4))THEN
          IF(need_i4)STOP 'check_kind_numbers(): FATAL - Required datatype INTEGER*4 is not available'
       END IF
       isize = BIT_SIZE(i4) / 8
       WRITE(*,*)'Using INTEGER*',isize,' instead'
    END IF

    IF(BIT_SIZE(i8).NE.64)THEN
       WRITE(*,*)'check_kind_numbers(): WARNING - INTEGER*8 not available'
       IF(PRESENT(need_i8))THEN
          IF(need_i8)STOP 'check_kind_numbers(): FATAL - Required datatype INTEGER*8 is not available'
       END IF
       isize = BIT_SIZE(i8) / 8
       WRITE(*,*)'Using INTEGER*',isize,' instead'
    END IF
    
    IF(real8byte.EQ.real4byte)WRITE(*,*)'check_kind_numbers(): WARNING - REAL*4 not available'

    RETURN
  END SUBROUTINE check_kind_numbers

END MODULE kind_numbers
