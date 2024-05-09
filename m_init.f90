module m_init 
    use precision, only: dp ! use dp as the data-type (and for defining constants!)
    use m_gldat
    use m_linspace
    use omp_lib
    implicit none

    public :: init_bounds 
    public :: init_radiator
    
    contains

    subroutine init_bounds(N) 
        ! ---------------------------------------------------------------------------------------
        ! Description: Initializes boundary conditions for a simulation domain.
        ! Inputs:
        !   - N: Integer representing the size of the domain (number of grid points along each dimension)
        ! Outputs:
        !   - None (since it modifies global variables u, u_old, x, y, z, L, and h)
        ! Notes:
        !   - Sets u and u_old to specific values at the domain boundaries.
        !   - Initializes x, y, and z as linearly spaced vectors within the domain.
        !   - Computes grid spacing h based on the domain length L.
        ! ---------------------------------------------------------------------------------------
        integer :: N
        integer :: i,j,k

        DO k = 1, N
            DO j = 1, N
                DO i = 1, N
                    u(i, j, k) = 0.0
                    u_old(i, j, k) = 0.0
                END DO
            END DO
        END DO
        

        !1 is -1 in the  box, N is 1
        u(:, 1, :) = 0.0d0
        u(:, N, :) = 20.0d0
        u(1, :, :) = 20.0d0
        u(N, :, :) = 20.0d0
        u(:, :, 1) = 20.0d0
        u(:, :, N) = 20.0d0

        u_old(:, 1, :) = 0.0d0
        u_old(:, N, :) = 20.0d0
        u_old(1, :, :) = 20.0d0
        u_old(N, :, :) = 20.0d0
        u_old(:, :, 1) = 20.0d0
        u_old(:, :, N) = 20.0d0

        !initialises x,y,z to be linearly spaced vectors in the domain
        call linspace(from=-1._dp, to=1._dp, array=x)
        call linspace(from=-1._dp, to=1._dp, array=y)
        call linspace(from=-1._dp, to=1._dp, array=z)

        L = 2.0     !length of the box
        h = L/real(N-1)     !grid spacing

    end subroutine init_bounds

    subroutine init_radiator(N)

        ! ---------------------------------------------------------------------------------------
        ! Description: Initializes a radiating source within the simulation domain.
        ! Inputs:
        !   - N: Integer representing the size of the domain (number of grid points along each dimension)
        ! Outputs:
        !   - None (since it modifies global variable f)
        ! Notes:
        !   - Sets f to a value of 200.0 within a specific region of the domain defined by conditions on x, y, and z.
        !   - The region is defined by inequalities involving x, y, and z.
        ! ---------------------------------------------------------------------------------------

        integer :: N
        integer :: i,j,k

        f = 1.0
    
        DO i = 1,N
            DO j = 1,N
                DO k = 1,N
                    IF (x(i) > -1.0 .AND. x(i) < -(3.0/8.0) .AND. &
                        y(j) > -1.0 .AND. y(j) < -(1.0/2.0) .AND. &
                        z(k) > -(2.0/3.0) .AND. z(k) < 0.0) THEN
                        f(i,j,k) = 200.0
                    END IF
                END DO
            END DO
        END DO
        
    end subroutine init_radiator
end module m_init