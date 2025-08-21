program main

    use types
    use params    
    use functions
    use io_utils
    use lattice 
    use basis
    use file_utils
    use hamiltonian
    use diagonalization
    use observables

    implicit none
    type(sim_params) :: par ! Simulation input parameters
    type(output) :: out ! Output parameters
    type(diag_params) :: diag ! Diagonalization parameters
    type(thread_params) :: thrd ! Thread parameters
    type(geometry) :: geo ! Geometry includes lattice-, basis- and symmetry related parameters 
    type(hamiltonian_params) :: ham ! Hamiltonian parameters
    type(system_state) :: st ! System state parameters
    
    call preprocess(par, geo, thrd, out, st)  
    call define_lattice(out%outdir, par, geo, out)
   
    do k1_= 0, geo%k1_max
        do k2_= 0, geo%k2_max
            ! Loop over k1 and k2 values
            ! k1 and k2 are the momentum values in the first Brillouin zone
            st%k1 = k1_
            st%k2 = k2_    
            if((st%k1 .ne. 0) .or.(st%k2 .ne. 0) .or. (geo%id == 2)) then 
                st%mat_type = "C"
            else 
                st%mat_type = "R"
            end if    
            call make_basis(par, geo, st)
            call printing(par, geo, thrd, diag, out)
            if(par%otf == 0) call generate_hamiltonian(par, geo, ham, out, thrd, st%k1, st%k2)
            if(geo%dim == 0) then
                diag%nev = int(geo%dim, kind=4)
                diag%ncv = int(geo%dim, kind=4)
                print*, 'No basis states found.'
                goto 116
            end if
        
            call set_thrds(par, thrd)
        
            !$omp parallel do default(firstprivate) shared(geo%sites, geo%particles, par%dim, diag%nev, diag%ncv, geo%occ, parities, dplcts, hamOff, hamDi, basis, hamOff_dp, hamDi_d, rcOff, rcDi, hamOff_dc, num_thrds) num_threads(v2_thrds)
            do nv2_ = 0, thrd%ndv2-1, 1 ! Loop over V2 values
                st%nv2 = nv2_
                st%v2  = par%v2min + st%nv2 * par%dv2
                !$ thrd_id = omp_get_thread_num()
        
                !$omp parallel do default(firstprivate) shared(geo%sites, geo%particles, par%dim, diag%nev, diag%ncv, geo%occ, parities, dplcts, hamOff, hamDi, basis, hamOff_dp, hamDi_d, rcOff, rcDi, hamOff_dc, num_thrds) num_threads(v1_thrds)
                do nv1_ = 0, thrd%ndv1-1, 1 ! Loop over V1 values
                    st%nv1 = nv1_
                    st%v1  = par%v1min + st%nv1 * par%dv1
                    !$ thrd_id_2 = omp_get_thread_num()
                    out%unit = 11 + thrd%units(thrd%thrd_id + 1, thrd%thrd_id_2 + 1)
                    !$omp critical 
                    do conf_ = 1, par%nDis !Disorder loop 
                        st%conf = conf_
                        call printing(st%conf, par%ti, st%k1, st%k2, st%v1, st%v2)
                        call diagonalize(par, geo, ham, diag, st, out)
                        if(curr == 1) call current_correlations(par, geo, diag, out, st, thrd)                                 

                    end do !Disorder loops
                    !$omp end critical 
                end do !vloop
                !$omp end parallel do
            end do !v2loop
            !$omp end parallel do
            call cleanup(.true., geo, ham, diag)   
            116 continue

            call timing(out%outdir, 1)
            call printing(par, geo, thrd, diag, out)
            call cleanup(.false., geo, ham, diag)
        end do 
    end do !Momentum loops 

end program main