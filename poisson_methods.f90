module poisson_methods_m

    ! Please do not put any computational variable in this module.
    ! Pass all variables as arguments to the subroutine.
    use precision, only: dp ! use dp as the data-type (and for defining constants!)
    use m_gldat
    use omp_lib
    implicit none
  
    private
    public :: jacobi
    public :: gauss_seidel
    public :: iterator 
    public :: calculate_residual
    public :: update_old_field
  
  contains

    
    subroutine jacobi(N)

      ! -----------------------------------------------------------------------------------------------------------
      ! Description: Performs one iteration of the Jacobi method for solving Poisson's equation.
      ! Inputs:
      !   - N: Integer representing the size of the domain (number of grid points along each dimension)
      ! Outputs:
      !   - Updates the values of the global array u based on neighboring values and the source term f.
      ! Notes:
      !   - The Jacobi method is an iterative numerical technique for solving linear systems of equations.
      !   - In this context, it is used to update the solution field u at each grid point.
      !   - The update formula involves contributions from neighboring grid points and the source term f.
      !   - The factor 1/6.0 is used to weight the contributions appropriately.
      !   - The loop indices i, j, and k iterate over the interior grid points (excluding boundaries).
      !   - The updated values are stored in the array u.
      ! -----------------------------------------------------------------------------------------------------------

    integer :: N
    integer :: i,j,k
    real(kind=dp) :: temp
    residual = 0.0
        
!$omp parallel do shared(u, u_old, f, N) private(i,j,k)  schedule(dynamic) reduction(+:residual)
        DO k = 2,N-1
          DO j = 2,N-1
              DO i = 2,N-1
                  ! Use temporary variable to reduce memory access
                  temp = u_old(i-1,j,k) + u_old(i+1,j,k) + u_old(i,j-1,k) + &
                         u_old(i,j+1,k) + u_old(i,j,k-1) + u_old(i,j,k+1)
                  u(i,j,k) = 1.0/6.0 * (h**2 * f(i,j,k) + temp)
                  residual = residual + (u(i, j, k) - u_old(i, j, k))**2
              END DO 
          END DO 
      END DO
      !$OMP END parallel do 
      residual = (1.0/real((N-2)**3)) * SQRT(residual)
      
    end subroutine jacobi

    subroutine gauss_seidel(N)
      integer :: N
      !for this routine please see the Gauss-Seidel implementation
    end subroutine gauss_seidel

    subroutine iterator(algorithm, N)
      
      integer :: N, algorithm
        
        ! No support for different algorithms at this time
      call jacobi(N);

    end subroutine iterator

    subroutine calculate_residual(N) 
      integer :: N
      !subroutine not in use for this implementation
    end subroutine calculate_residual

    subroutine update_old_field(N)
    ! -----------------------------------------------------------------------------------------------------------
    ! Description: Updates the previous solution field u_old with the current solution field u.
    ! Inputs:
    !   - None (since it operates on global variables u and u_old)
    ! Outputs:
    !   - Updates the values of u_old to match the current solution u.
    ! Notes:
    !   - This subroutine is called after each iteration to maintain the previous solution.
    ! -----------------------------------------------------------------------------------------------------------
      integer :: i,j,k,N
      !$omp parallel do shared(u, u_old, N) private(i, j, k)
      DO k = 1, N
        DO j = 1, N
            DO i = 1, N
                u_old(i, j, k) = u(i,j,k)
            END DO
        END DO
    END DO
    !$omp end parallel do
    end subroutine update_old_field
  
  end module poisson_methods_m
  