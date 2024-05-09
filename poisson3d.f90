program poisson3d

    ! Import different methods
    use precision
    use poisson_methods_m
    use m_gldat
    use m_alloc
    use m_init
    use omp_lib
    use m_input_output
  
    integer, save :: N = 100
    real(dp), save :: T0 = 0._dp
    integer, save :: itermax = 1000
    real(dp), save :: tolerance = 1.e-4_dp
    real(dp) :: t_start, t_end

   
    call read_namelist()   ! uncomment to read vals from namelist

    call alloc(N)   !allocates u, u_old, f fields and x,y,z vectors

    call init_bounds(N)   !initialise domain and boundary values

    call init_radiator(N)   !initialise f array

    iter = 0
    residual = 10 !initialise it with a nonzero value
    n_points = N
    n_algorithm = algorithm

    !$OMP MASTER
    t_start = omp_get_wtime() 
    !$OMP END MASTER

    do while (residual > tolerance .and. iter < itermax)    !quit if one of these is untrue

      call iterator(algorithm, N)   
      call update_old_field(N)

      iter = iter + 1
    end do


    !$OMP MASTER
    t_end = omp_get_wtime()
    wall_clock_time = t_end - t_start
    call write_performance()
    !$OMP END MASTER
  
  
  end program poisson3d
  