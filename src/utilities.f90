module utilities 

    use functions
    implicit none

    contains

subroutine preprocess()
    ! Preprocessing routine to set up necessary parameters and variables
    use variables
    use parameters
    use input_variables
    use types
    use file_utils
    use io_routines
    use lattice
    use basis
    use hamiltonian
    use printing_routines
    use symmetries
    use utilities

    implicit none

    call read_input(par)
    call setup_output_directory()
    call create_output_subdirs(outdir)
    call nsteps(v1min, v1max, dv1, ndv1)
    call nsteps(v2min, v2max, dv2, ndv2)
    call characters(symm, irrep, geo_par) ! Set characters according to chosen irrep 
    call timing(outdir, 0)
    call setvars()  
    call check_parallel()
    call nsteps(par%v1min, par%v1max, par%dv1, ndv1)
    call nsteps(par%v2min, par%v2max, par%dv2, ndv2)
    call stepunits(1, ndv1, ndv2, units_2)
    call threadunits(ndv1, ndv2, units)
    

    end subroutine preprocess

    subroutine cleanup(inner)
        
        use variables

        implicit none
        logical, intent(in) :: inner ! If true, deallocate arrays used within parameter-loops, otherwise deallocate arrays outside of parameter-loops.
        ! This subroutine deallocates all dynamically allocated arrays used in the program.

        if(inner) then ! Deallocate all dynamically allocated arrays of inner loops
            if(allocated(rc))          deallocate(rc)
            if(allocated(ham))         deallocate(ham)
            if(allocated(ham_dc))      deallocate(ham_dc)
            if(allocated(evals))       deallocate(evals)
            if(allocated(ham_dp))      deallocate(ham_dp)
            if(allocated(norm))        deallocate(norm)
            if(allocated(hamDi))       deallocate(hamDi)
            if(allocated(rcOff))       deallocate(rcOff)    
            if(allocated(abasis))      deallocate(abasis)
            if(allocated(bbasis))      deallocate(bbasis)
            if(allocated(hamOff))      deallocate(hamOff)
            if(allocated(mombasis))    deallocate(mombasis)
            if(allocated(parities))    deallocate(parities)
            if(allocated(eigstate))    deallocate(eigstate)
            if(allocated(energies))    deallocate(energies)
            if(allocated(hamOff_dp))   deallocate(hamOff_dp)
            if(allocated(hamOff_dc))   deallocate(hamOff_dc)
            if(allocated(occ))         deallocate(occ)
            if(allocated(dplcts))      deallocate(dplcts)
            if(allocated(eigstate_dc)) deallocate(eigstate_dc)
            if(allocated(basis_states)) deallocate(basis_states)
        else ! Deallocate all dynamically allocated arrays of outer loops
            if(allocated(xy))          deallocate(xy)
            if(allocated(xtransl))     deallocate(xtransl)
            if(allocated(ytransl))     deallocate(ytransl)
            if(allocated(alattice))    deallocate(alattice)
            if(allocated(blattice))    deallocate(blattice)
            if(allocated(phases))      deallocate(phases)
            if(allocated(bsites))      deallocate(bsites)
            if(allocated(hexsites))    deallocate(hexsites)
            if(allocated(asitesbonds)) deallocate(asitesbonds)
            if(allocated(bsitesbonds)) deallocate(bsitesbonds)
            if(allocated(geopar%latticevecs)) deallocate(geopar%latticevecs)
        end if 


    end subroutine cleanup

    subroutine timing(dir, start_end)
        implicit none

        integer, intent(in) :: start_end ! 0 = start, 1 = end
        character(len=*), intent(in) :: dir ! Directory for output
        character(8), save  :: datei, datef
        character(10), save :: timei, timef
        integer, save :: values_i(8), values_f(8)
        real, save :: start, finish
        character :: file*512


        if(start_end == 0) then
            call cpu_time(start)
            call date_and_time(date=datei,time=timei,values=values_i)
            write(*,"('Calculation started at',x,a,' h',x,a,' min',x,a,' sec')") timei(1:2), timei(3:4), timei(5:6)
            print*, ''
            print*, 'Start date: ',datei(7:8), '.',datei(5:6), '.',datei(1:4)
            print*, ''
        else
            call cpu_time(finish)
            print*, 'Start = ', start
            print*, 'Finish = ', finish
            if (finish - start < 60) then
                write(*,"(' Elapsed CPU time = ',f12.3,' seconds.')") finish-start
                print*, ''
            else if (finish - start < 3600) then
                write(*,"(' Elapsed CPU time = ',f12.3,' minutes.')") (finish-start)/60
                print*, ''
            else
                write(*,"(' Elapsed CPU time = ',f12.3,' hours.')") (finish-start)/3600
                print*, ''
            end if
            call date_and_time(date=datef,time=timef,values=values_f)
            write(*,"(' Calculation started at',x,a,'h',x,a,'min',x,a,'sec')") timei(1:2), timei(3:4), timei(5:6)
            print*, ''
            write(*,"(' Calculation ended at',x,a,'h',x,a,'min',x,a,'sec')") timef(1:2), timef(3:4), timef(5:6)
            print*, ''
            print*, 'Start date: ',datei(7:8), '.',datei(5:6), '.',datei(1:4)
            print*, ''
            print*, 'End date: ',datef(7:8), '.',datef(5:6), '.',datef(1:4)

            file=trim(dir)//"timing.dat"
            open(77,file = trim(file))
            write(77,"(' Elapsed CPU time = ',f20.10,' seconds.')") finish-start
            write(77,"(' Calculation started at',x,a,'h',x,a,'min',x,a,'sec')") timei(1:2), timei(3:4), timei(5:6)
            write(77,"(' Start date:',x,a,'.',x,a,'.',x,a)") ,datei(7:8), datei(5:6), datei(1:4)
            write(77,"(' Calculation ended at',x,a,'h',x,a,'min',x,a,'sec')") timef(1:2), timef(3:4), timef(5:6)
            write(77,"(' End date:',x,a,'.',x,a,'.',x,a)") ,datef(7:8), datef(5:6), datef(1:4)
            close(77)
        end if
    end subroutine timing

    !------------------------------------------!
    !            Set initial variables         !
    !------------------------------------------!

    subroutine setvars()

        use variables
        use input_variables
        implicit none

        character*10:: clusterl
        integer :: threads

        !call KMP_SET_STACKSIZE_S(1000000000)

        !$ call omp_set_dynamic(dynamic)
        !$ call omp_set_nested(nested)
        !$ call omp_set_num_threads(othrds)
        call mkl_set_num_threads(mthrds)
        call mkl_get_max_threads(threads)

        print*,'MKL THREADS',threads
        !$ threads = omp_get_num_threads()
        print*,'OMP THREADS',threads

        if(tilted == 0) then 
            sites = 2*ucx*ucy 
        else if(cluster(1:2) == "16") then 
            sites = 16 
        else if(cluster(1:2) == "18") then 
            sites = 18 
        else if(cluster(1:2) == "20") then 
            sites = 20 
        else if(cluster(1:2) == "24") then 
            sites = 24
        else if(cluster(1:2) == "30") then 
            sites = 30 
        else if(cluster(1:2) == "32") then 
            sites = 32
        else if(cluster(1:1) == "L") then 
            clusterl = cluster
            read(clusterl(2:3),*) sites
        end if 
        particles = int(filling * sites)
        

        if(symm == 1) print('(1x, a,a)'), 'Irrep = ', irrep
        if(tilted == 1) then
            print('(1x, a,a)'), 'Cluster = ', cluster
        else 
            print('(1x, a,i0)'), 'UCX = ', ucx
            print('(1x, a,i0)'), 'UCY = ', ucy
        end if 
        print('(1x, a,i0)'), 'Sites = ', sites
        print('(1x, a,i0)'), 'Particles = ', particles
        print*, ''
        

        return

    end subroutine setvars


    !------------------------------------------!
    !            Check parallelization         !
    !------------------------------------------!

    subroutine check_parallel()

        implicit none
        integer :: num_threads_loc
        integer :: thread_num_loc
        integer :: max_threads_loc

        num_threads_loc = 0
        thread_num_loc = 0
        max_threads_loc = 0
        !Critical block is a lock, forcing the instructions to be run in serial
        print*, ''
        !$omp parallel
            !$omp critical
                !$ thread_num_loc = omp_get_thread_num()
                print('(1x,100(a,i0))'), 'Thread number ', thread_num_loc, ' is online.'
                !$ num_threads_loc = omp_get_num_threads()
                if (thread_num_loc == 0) then
                    print*, ''
                    print('(1x,100(a,i0))'), 'Available number of threads = ', num_threads_loc
                    !$ max_threads_loc = omp_get_max_threads()
                    print('(1x,100(a,i0))'), 'Maximum number of threads = ', max_threads_loc
                    print*, ''
                end if
            !$omp end critical
        !$omp end parallel
        print*, ''
        return 

    end subroutine

    !-------------------------------!
    !            Step units         !
    !-------------------------------!

    subroutine stepunits(ndu, ndv, ndv2, units)
        !Creates array with units for parallel I/O. Each entry is a different unit number.
        implicit none
        integer, intent(in) :: ndu, ndv, ndv2
        integer, allocatable, intent(out) ::  units(:,:,:,:)
        integer :: i, j, k, l, q

        if (allocated(units)) deallocate(units)
        allocate(units(ndv2,ndv,ndu,2))
        q = 0
        do i = 1, ndv2
            do j = 1, ndv
                do k = 1, ndu
                    do l = 1, 2
                        q = q + 1
                        units(i,j,k,l) = q
                    end do
                end do
            end do
        end do

    end subroutine

    !---------------------------------!
    !            Thread units         !
    !---------------------------------!

    subroutine threadunits(ndv, ndv2, units)

        implicit none
        integer, intent(in) :: ndv, ndv2
        integer, allocatable, intent(out) :: units(:, :)
        integer :: i, j, h
        integer :: n_threads_tot

        n_threads_tot = 1

        if (ndv > 0) n_threads_tot = n_threads_tot + ndv
        if (ndv2 > 0) n_threads_tot = n_threads_tot + ndv2

        if(allocated(units)) deallocate(units)
        allocate(units(n_threads_tot,n_threads_tot))
        h = 0
        do i = 1, n_threads_tot
            do j = 1, n_threads_tot
                h = h + 1
                units(i,j) = h
            end do
        end do

        return 

    end subroutine


    !---------------------------------------------------!
    !            Number of discretization steps         !
    !---------------------------------------------------!

    subroutine nsteps(min, max, delta, steps)
        use types
        implicit none
        ! type(geometry) :: params
        double precision, intent(in) :: min, max, delta
        integer, intent(out) :: steps

        if(delta <= 1.d0) then
            steps = int(abs(max-min)/delta + delta/2) + 1
        else
            steps = int(abs(max-min)/delta) + 1
        end if


    end subroutine nsteps

end module utilities