MODULE m_gldat
    use precision
    implicit none
    real(dp), allocatable, dimension(:,:,:) :: u, u_old, f
    real(dp), allocatable, dimension(:) :: x, y, z
    integer :: iostat, iter, n_points, n_algorithm, threads
    real(dp) :: h,L, residual
    real(dp) :: wall_clock_time
END MODULE m_gldat