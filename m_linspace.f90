module m_linspace
    use precision, only: dp
    implicit none 
    contains

subroutine linspace(from, to, array)

    ! --------------------------------------------------------------------------------
    ! Description: Generates an array of equally spaced values between 'from' and 'to'.
    ! A replacement for the python "linspace" function. 
    ! Inputs:
    !   - from: Starting value (real, intent(in))
    !   - to: Ending value (real, intent(in))
    ! Outputs:
    !   - array: Resulting array of equally spaced values (real, intent(out))
    ! Notes:
    !   - The array must be pre-allocated with the desired size.
    !   - If the array size is zero, no values are generated.
    !   - If the array size is one, the 'from' value is assigned directly.
    !   - Otherwise, the array is filled with equally spaced values.
    !     The step size is calculated as (to - from) / (n - 1), where 'n' is the array size.
    !     Each element is computed as: array(i) = from + (i - 1) * step.
    ! --------------------------------------------------------------------------------


    real(dp), intent(in) :: from, to
    real(dp), intent(out) :: array(:)
    real(dp) :: range
    integer :: n, i
    n = size(array)
    range = to - from

    if (n == 0) return

    if (n == 1) then
        array(1) = from
        return
    end if

    do i=1, n
        array(i) = from + range * (i - 1) / (n - 1)
    end do


end subroutine
end module m_linspace