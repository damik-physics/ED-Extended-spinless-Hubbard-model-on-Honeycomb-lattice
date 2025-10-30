module hubbard_utilities
    ! Utility module for honeycomb Hubbard model calculations
    implicit none
    
contains

    subroutine parameter_sweep()
        ! Example of how to sweep interaction parameters
        implicit none
        integer, parameter :: n_points = 10
        double precision :: V_values(n_points)
        double precision :: ground_energies(n_points)
        integer :: i
        
        ! Create parameter sweep
        do i = 1, n_points
            V_values(i) = 0.1d0 * (i - 1)
        end do
        
        write(*,*) '=== Parameter Sweep Example ==='
        write(*,*) 'V1 values:', V_values
        write(*,*) 'This would require modifying the main program to loop over V1'
        write(*,*)
    end subroutine parameter_sweep

    subroutine analyze_correlations(N_sites, wavefunction)
        ! Calculate density-density correlations
        implicit none
        integer, intent(in) :: N_sites
        double precision, intent(in) :: wavefunction(:)
        
        write(*,*) '=== Correlation Analysis ==='
        write(*,*) 'Ground state wavefunction norm:', sum(wavefunction**2)
        write(*,*) 'This is where you would calculate:'
        write(*,*) '  - Density correlations <n_i n_j>'
        write(*,*) '  - Current correlations'
        write(*,*) '  - Structure factors'
        write(*,*)
    end subroutine analyze_correlations

    subroutine save_results(filename, energies, wavefunction)
        ! Save results to file
        implicit none
        character(len=*), intent(in) :: filename
        double precision, intent(in) :: energies(:)
        double precision, intent(in) :: wavefunction(:)
        
        integer :: unit, i
        
        open(newunit=unit, file=filename, status='replace')
        write(unit,'(A)') '# Hubbard Model Results'
        write(unit,'(A,I0)') '# Basis dimension: ', size(wavefunction)
        write(unit,'(A)') '# Eigenvalues:'
        
        do i = 1, min(10, size(energies))
            write(unit,'(I3,F16.8)') i, energies(i)
        end do
        
        write(unit,'(A)') '# Ground state wavefunction (first 20 components):'
        do i = 1, min(20, size(wavefunction))
            write(unit,'(I5,F16.8)') i, wavefunction(i)
        end do
        
        close(unit)
        write(*,*) 'Results saved to: ', trim(filename)
    end subroutine save_results

    function calculate_filling(N_sites, N_particles) result(filling)
        ! Calculate electron filling fraction
        implicit none
        integer, intent(in) :: N_sites, N_particles
        double precision :: filling
        
        filling = dble(N_particles) / dble(N_sites)
    end function calculate_filling

    subroutine print_physics_summary(N_sites, N_particles, t, V1)
        ! Print physical interpretation
        implicit none
        integer, intent(in) :: N_sites, N_particles
        double precision, intent(in) :: t, V1
        
        double precision :: filling, interaction_ratio
        
        filling = calculate_filling(N_sites, N_particles)
        interaction_ratio = V1 / t
        
        write(*,*) '=== Physics Summary ==='
        write(*,'(A,F6.3)') 'Electron filling: ', filling
        write(*,'(A,F6.3)') 'Interaction ratio V1/t: ', interaction_ratio
        write(*,*)
        
        if (interaction_ratio < 0.5d0) then
            write(*,*) 'Regime: Weakly interacting (metallic behavior expected)'
        else if (interaction_ratio < 2.0d0) then
            write(*,*) 'Regime: Intermediate coupling (competing effects)'
        else
            write(*,*) 'Regime: Strongly interacting (insulating behavior expected)'
        end if
        write(*,*)
        
        if (abs(filling - 0.5d0) < 1e-6) then
            write(*,*) 'Special case: Half filling - may show Mott insulator physics'
        end if
        write(*,*)
    end subroutine print_physics_summary

end module hubbard_utilities