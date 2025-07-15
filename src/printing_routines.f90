module printing_routines 
    implicit none
    

    interface printing
        module procedure print_configuration, &
                         print_parameters!, &
                         !print_results
    end interface printing 
    contains

    subroutine print_configuration(conf, ti, k1, k2, v1, v2)
        implicit none
        integer, intent(in)          :: conf, ti, k1, k2
        double precision, intent(in) :: v1, v2
        ! Print the current configuration of the system
        print*, '----------------------------------'
        print*, 'Current Configuration: ', conf 
        if(ti == 1) print*, 'k1, k2:', k1, k2
        print*, 'V1:', v1
        print*, 'V2:', v2
        print*, '----------------------------------'

    end subroutine print_configuration

    subroutine print_parameters()
        ! Print the parameters used in the simulation

        use input_variables
        use variables
        use parameters 
        use types
        implicit none
        type(sim_params) par
        type(geometry) geo
        type(thread_params) thr
        type(out_params) out  
        
        print*, '----------------------------------'
        print*, 'Simulation Parameters:'
        print*, 'Sites:', geo%sites
        print*, 'Particles:', geo%particles
        print*, 'Dimension:', geo%dim
        print*, 'V1-range:', par%v1min, '-', par%v1max, 'in steps of ', par%dv1
        print*, 'V2-range:', par%v2min, '-', par%v2max, 'in steps of ', par%dv2
        print*, 'Number of eigenvalues (nev):', par%nev
        print*, 'Number of convergence vectors (ncv):', par%ncv0
        print*, 'Number of disorder configurations:', par%nDis
        if(dis > 0.d0) print*, 'Disorder strength:', par%dis
        if(mass > 0.d0) print*, 'Sublattice imbalance:', par%mass
        print*, 'Filling fraction:', par%filling
        if(t /= 1.d0) print*, 'Hopping strength (t):', par%t
        if(g_fact /= 2) print*, 'Electron g-factor:', par%g_fact
        print*, 'Diagonalization method:'
        if (feast == 1) then
            print*, '  FEAST'
        else if (arpack == 1) then
            print*, '  Arpack'
        else if (mkl == 1) then
            print*, '  MKL'
        else if (exact == 1) then
            print*, '  Exact diagonalization'
        else
            print*, '  No diagonalization method selected'
        end if
        if(othrds > 1) print*, 'Number of OpenMP threads for Arpack:', thr%othrds
        if(mthrds > 1) print*, 'Number of OpenMP threads for MKL:', thr%mthrds
        if(thr%dis_thrds > 1) print*, 'Disorder threads:', thr%dis_thrds
        if(thr%v1_thrds > 1) print*, 'V1 threads:', thr%v1_thrds
        if(thr%v2_thrds > 1) print*, 'V2 threads:', thr%v2_thrds
        if(thr%num_thrds > 1) print*, 'Number of threads for current calculations:', thr%num_thrds
        if(corr == 1) print*, 'Correlation calculations enabled.'
        if(curr == 1) print*, 'Current-current calculations enabled:', par%curr
        if(curr == 1) print*, 'Reference bonds for current calculations:', par%refbonds
        print*, 'degflag threshold:', out%deg
        if(states == 1) print*, 'Saving of eigenstates enabled.'
        print*, '----------------------------------'
    end subroutine print_parameters

 
end module printing_routines