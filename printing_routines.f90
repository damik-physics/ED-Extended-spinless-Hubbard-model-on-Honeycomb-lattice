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
        use variables
        use input_variables
        use variables
        use parameters 
        implicit none
        
        
        print*, '----------------------------------'
        print*, 'Simulation Parameters:'
        print*, 'Sites:', sites
        print*, 'Particles:', particles
        print*, 'Dimension:', dim
        print*, 'V1-range:', v1min, '-', v1max, 'in steps of ', dv1
        print*, 'V2-range:', v2min, '-', v2max, 'in steps of ', dv2
        print*, 'Number of eigenvalues (nev):', nev
        print*, 'Number of convergence vectors (ncv):', ncv
        print*, 'Number of disorder configurations:', nDis
        if(dis > 0.d0) print*, 'Disorder strength:', dis
        if(mass > 0.d0) print*, 'Sublattice imbalance:', mass
        print*, 'Filling fraction:', filling
        if(t /= 1.d0) print*, 'Hopping strength (t):', t
        if(g_fact /= 2) print*, 'Electron g-factor:', g_fact
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
        if(othrds > 1) print*, 'Number of OpenMP threads for Arpack:', othrds
        if(mthrds > 1) print*, 'Number of OpenMP threads for MKL:', mthrds
        if(dis_thrds > 1) print*, 'Disorder threads:', dis_thrds
        if(v1_thrds > 1) print*, 'V1 threads:', v1_thrds
        if(v2_thrds > 1) print*, 'V2 threads:', v2_thrds
        if(num_thrds > 1) print*, 'Number of threads for current calculations:', num_thrds
        if(debug == 1) print*, 'Debug mode enabled.'
        if(corr == 1) print*, 'Correlation calculations enabled.'
        if(curr == 1) print*, 'Current-current calculations enabled:', curr
        if(curr == 1) print*, 'Reference bonds for current calculations:', refbonds
        print*, 'Degeneracy threshold:', deg
        if(states == 1) print*, 'Saving of eigenstates enabled.'
        print*, '----------------------------------'
    end subroutine print_parameters

    ! subroutine print_results()
    !     use variables
    !     implicit none
    !     ! Print the results of the simulation
    !     print*, 'Simulation Results:'
    !     print*, 'Energies:', energies
    !     print*, 'Eigenstates:', eigstate
    ! end subroutine print_results

end module printing_routines