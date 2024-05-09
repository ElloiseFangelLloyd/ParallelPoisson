MODULE m_alloc
    use m_gldat
    use precision
    implicit none
    !allocation of necessary fields
    contains
    subroutine alloc(N)

    ! ---------------------------------------------------------------------------------------
    ! Description: Allocates memory for arrays u, u_old, f, x, y, and z.
    ! Inputs:
    !   - N: Integer representing the size of the domain (number of grid points along each dimension)
    ! Outputs:
    !   - Allocates memory for the arrays u, u_old, f, x, y, and z.
    !   - Checks the allocation status using the 'iostat' variable.
    !   - If allocation fails, prints an error message and stops execution.
    ! Notes:
    !   - The arrays u, u_old, and f are three-dimensional.
    !   - The vectors x, y, and z are one-dimensional.
    ! ---------------------------------------------------------------------------------------

        integer :: N
            allocate(u(N,N,N), stat=iostat)
            call check_iostat(iostat, "Could not allocate 'u' matrix!")
            allocate(u_old(N,N,N), stat=iostat)
            call check_iostat(iostat, "Could not allocate 'u_old' matrix!")
            allocate(f(N,N,N), stat=iostat)
            call check_iostat(iostat, "Could not allocate 'f' matrix!")
            allocate(x(N), stat=iostat)
            call check_iostat(iostat, "Could not allocate 'x' vector!")
            allocate(y(N), stat=iostat)
            call check_iostat(iostat, "Could not allocate 'y' vector!")
            allocate(z(N), stat=iostat)
            call check_iostat(iostat, "Could not allocate 'z' vector!")
    end subroutine alloc

    !Move this function here to avoid circular dependencies
    subroutine check_iostat(iostat, msg)
    ! ---------------------------------------------------------------------------------------
    ! Description: Allocates memory for arrays u, u_old, f, x, y, and z.
    ! Inputs:
    !   - N: Integer representing the size of the domain (number of grid points along each dimension)
    ! Outputs:
    !   - Allocates memory for the arrays u, u_old, f, x, y, and z.
    !   - Checks the allocation status using the 'iostat' variable.
    !   - If allocation fails, prints an error message and stops execution.
    ! Notes:
    !   - The arrays u, u_old, and f are three-dimensional.
    !   - The vectors x, y, and z are one-dimensional.
    ! ---------------------------------------------------------------------------------------
        integer, intent(in) :: iostat
        character(len=*), intent(in) :: msg
    
        if ( iostat == 0 ) return
    
        write(*,'(a,i0,/,tr3,a)') 'ERROR = ', iostat, msg
        stop
    
      end subroutine check_iostat
END MODULE m_alloc