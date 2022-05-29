!--------------------------------
! module for useful subroutines
!-------------------------------

module use_subroutines
    implicit none
    contains
        !create an array of bins in linear space
        subroutine linspace(from, to, array)
            implicit none
            real, intent(in)  :: from, to
            real, intent(out) :: array(:)
            real              :: range
            integer           :: n, i
            n     = size(array)
            range = to - from
 
            if (n .eq. 0) return            
            
            if (n .eq. 1) then
                array(1) = from
                return
            end if

            do i = 1, n
                array(i) = from + range * (i - 1) / (n - 1)
            end do
            
        end subroutine linspace 
        
        !create an array of bins in logspace
        subroutine logspace(from,to,array)
            implicit none
            real,intent(in)  :: from, to
            real,intent(out) :: array(:)
            if (size(array) .eq. 1 .and. from .ne. to) then
                write(0,'("ERROR: Cannot call logspace with n=1 and xmin /= xmax")')
                stop
            end if
            call linspace(log10(from),log10(to),array)
            array = 10.**array
        end subroutine logspace




end module use_subroutines
