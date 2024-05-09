module m_input_output
    
    use precision, only: dp
    use m_gldat
    use omp_lib, only: omp_get_num_threads
    
    implicit none 
    
    contains

    subroutine write_performance()

        ! -----------------------------------------------------------------
        ! Produces a file "performance.csv" in the root folder
        ! Contains runtime information for testing of parallelism
        ! -----------------------------------------------------------------
        integer :: ierr
        character(8)  :: date
        character(10) :: time
        character(200) :: line

        ! Get time stamp
        call date_and_time(date,time)

         ! Write performance to file
        open(10, file='performance.csv', status='unknown', action='write', position='append', form='formatted', iostat=ierr)

        ! Format: 
        ! Timestamp, algorithm, iterations, time, n_threads, n_points
        
        ! Write the new line to the file
        write(line, '(A, I10, I10, F12.6, I10, I10)') date//time, n_algorithm, iter, wall_clock_time, threads, n_points
        
        write(10, *) TRIM(line)

        ! Close the file
        close(10)

    end subroutine
end module m_input_output
