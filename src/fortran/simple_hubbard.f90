program simple_hubbard_model
    ! Simplified honeycomb Hubbard model solver
    ! Focuses on core physics: hopping + nearest-neighbor interactions
    
    implicit none
    
    ! Physical parameters
    integer, parameter :: Lx = 4, Ly = 4                    ! Lattice dimensions
    integer, parameter :: N_sites = 2 * Lx * Ly             ! Total sites (A + B sublattices)
    integer, parameter :: N_particles = N_sites / 2         ! Half filling
    double precision, parameter :: t = 1.0d0                ! Hopping strength
    double precision, parameter :: V1 = 0.5d0               ! Nearest-neighbor interaction
    
    ! Computational parameters (calculated at runtime)
    integer(8) :: N_basis                                    ! Hilbert space dimension
    integer :: N_evals                                       ! Number of eigenvalues to find
    
    ! Arrays
    integer(8), allocatable :: basis_states(:)              ! Many-body basis
    integer, allocatable :: bonds(:,:)                      ! Nearest-neighbor bonds
    double precision, allocatable :: hamiltonian(:,:)       ! Hamiltonian matrix
    double precision, allocatable :: eigenvalues(:)         ! Ground state energies
    double precision, allocatable :: eigenvectors(:,:)      ! Ground state wavefunctions
    
    integer :: N_bonds
    
    ! Calculate basis dimension
    N_basis = choose(N_sites, N_particles)
    N_evals = min(20, int(N_basis))
    
    write(*,*) '=== Simple Honeycomb Hubbard Model ==='
    write(*,*) 'Sites:', N_sites, 'Particles:', N_particles
    write(*,*) 'Basis dimension:', N_basis
    write(*,*) 'Hopping t =', t, 'Interaction V1 =', V1
    write(*,*)
    
    ! Build the physics
    call generate_basis()
    call build_lattice()
    call construct_hamiltonian()
    call diagonalize_hamiltonian()
    call analyze_results()
    
    ! Cleanup
    deallocate(basis_states, bonds, hamiltonian, eigenvalues, eigenvectors)
    
contains

    subroutine generate_basis()
        ! Generate all possible N_particles-electron configurations on N_sites
        implicit none
        integer(8) :: i, count, state
        integer :: j, n_occupied
        
        allocate(basis_states(N_basis))
        
        count = 0
        do i = 0, 2_8**N_sites - 1
            ! Count occupied sites in this configuration
            n_occupied = 0
            state = i
            do j = 0, N_sites - 1
                if (btest(state, j)) n_occupied = n_occupied + 1
            end do
            
            ! Keep configurations with exactly N_particles electrons
            if (n_occupied == N_particles) then
                count = count + 1
                if (count <= N_basis) basis_states(count) = i
            end if
        end do
        
        write(*,*) 'Generated', count, 'basis states'
    end subroutine generate_basis

    subroutine build_lattice()
        ! Create honeycomb lattice bonds: each A site connects to 3 B neighbors
        implicit none
        integer :: ix, iy, site_A, site_B, bond_count, j
        integer :: dx(3), dy(3)  ! Relative positions of B neighbors
        
        ! Honeycomb bond vectors (A to B connections)
        dx = [0, 1, 0]
        dy = [0, 0, 1]
        
        N_bonds = 3 * Lx * Ly  ! Each unit cell has 3 bonds
        allocate(bonds(2, N_bonds))
        
        bond_count = 0
        do iy = 0, Ly - 1
            do ix = 0, Lx - 1
                site_A = 2 * (iy * Lx + ix) + 1  ! A sublattice (odd sites)
                
                ! Connect to 3 B neighbors with periodic boundary conditions
                do j = 1, 3
                    site_B = 2 * (mod(iy + dy(j), Ly) * Lx + mod(ix + dx(j), Lx)) + 2
                    bond_count = bond_count + 1
                    bonds(1, bond_count) = site_A
                    bonds(2, bond_count) = site_B
                end do
            end do
        end do
        
        write(*,*) 'Built lattice with', N_bonds, 'bonds'
    end subroutine build_lattice

    subroutine construct_hamiltonian()
        ! Build the many-body Hamiltonian matrix
        implicit none
        integer(8) :: i, j, new_state, pos
        integer :: bond, site1, site2
        double precision :: matrix_element
        
        allocate(hamiltonian(N_basis, N_basis))
        hamiltonian = 0.0d0
        
        do i = 1, N_basis
            ! Diagonal terms: interaction energy
            do bond = 1, N_bonds
                site1 = bonds(1, bond)
                site2 = bonds(2, bond)
                
                ! V1 * n_i * n_j interaction between neighboring sites
                if (btest(basis_states(i), site1-1) .and. btest(basis_states(i), site2-1)) then
                    hamiltonian(i, i) = hamiltonian(i, i) + V1
                end if
            end do
            
            ! Off-diagonal terms: hopping
            do bond = 1, N_bonds
                site1 = bonds(1, bond)
                site2 = bonds(2, bond)
                
                ! Hopping: c†_i c_j + c†_j c_i
                call try_hopping(basis_states(i), site1, site2, new_state, matrix_element)
                if (matrix_element /= 0.0d0) then
                    call find_state_index(new_state, pos)
                    if (pos > 0) hamiltonian(i, pos) = hamiltonian(i, pos) + matrix_element
                end if
                
                call try_hopping(basis_states(i), site2, site1, new_state, matrix_element)
                if (matrix_element /= 0.0d0) then
                    call find_state_index(new_state, pos)
                    if (pos > 0) hamiltonian(i, pos) = hamiltonian(i, pos) + matrix_element
                end if
            end do
        end do
        
        write(*,*) 'Constructed Hamiltonian matrix'
    end subroutine construct_hamiltonian

    subroutine try_hopping(state, from, to, new_state, matrix_element)
        ! Attempt electron hopping from site 'from' to site 'to'
        implicit none
        integer(8), intent(in) :: state
        integer, intent(in) :: from, to
        integer(8), intent(out) :: new_state
        double precision, intent(out) :: matrix_element
        
        
        matrix_element = 0.0d0
        new_state = 0
        
        ! Check if hopping is possible: occupied 'from', empty 'to'
        if (btest(state, from-1) .and. .not. btest(state, to-1)) then
            ! Perform hopping
            new_state = ibset(ibclr(state, from-1), to-1)
            
            ! Calculate fermionic sign
            sign_factor = 1
            bit_count = 0
            do i = min(from-1, to-1), max(from-1, to-1) - 1
                if (btest(state, i)) bit_count = bit_count + 1
            end do
            if (mod(bit_count, 2) == 1) sign_factor = -sign_factor
            
            matrix_element = -t * sign_factor  ! Negative sign in tight-binding
        end if
    end subroutine try_hopping

    subroutine find_state_index(target_state, index)
        ! Binary search for state in basis
        implicit none
        integer(8), intent(in) :: target_state
        integer(8), intent(out) :: index
        integer(8) :: left, right, mid
        
        left = 1
        right = N_basis
        index = 0
        
        do while (left <= right)
            mid = (left + right) / 2
            if (basis_states(mid) == target_state) then
                index = mid
                return
            else if (basis_states(mid) < target_state) then
                left = mid + 1
            else
                right = mid - 1
            end if
        end do
    end subroutine find_state_index

    subroutine diagonalize_hamiltonian()
        ! Find ground state using LAPACK
        implicit none
        integer :: info, lwork
        double precision, allocatable :: work(:)
        
        allocate(eigenvalues(N_basis))
        allocate(eigenvectors(N_basis, N_basis))
        
        eigenvectors = hamiltonian  ! DSYEV destroys input matrix
        
        ! Query optimal workspace
        lwork = -1
        allocate(work(1))
        call dsyev('V', 'U', int(N_basis), eigenvectors, int(N_basis), eigenvalues, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        
        ! Actual diagonalization
        allocate(work(lwork))
        call dsyev('V', 'U', int(N_basis), eigenvectors, int(N_basis), eigenvalues, work, lwork, info)
        deallocate(work)
        
        if (info /= 0) then
            write(*,*) 'ERROR: Diagonalization failed with info =', info
            stop
        end if
        
        write(*,*) 'Diagonalization completed successfully'
    end subroutine diagonalize_hamiltonian

    subroutine analyze_results()
        ! Print key results
        implicit none
        integer :: i
        double precision :: energy_per_site
        
        write(*,*)
        write(*,*) '=== RESULTS ==='
        write(*,'(A,F12.6)') 'Ground state energy: ', eigenvalues(1)
        
        energy_per_site = eigenvalues(1) / N_sites
        write(*,'(A,F12.6)') 'Energy per site: ', energy_per_site
        
        write(*,*)
        write(*,*) 'Lowest eigenvalues:'
        do i = 1, min(5_8, N_basis)
            write(*,'(I3,F12.6)') i, eigenvalues(i)
        end do
        
        ! Check for degeneracies
        if (N_basis > 1) then
            write(*,*)
            write(*,'(A,F8.4)') 'Gap to first excited state: ', eigenvalues(2) - eigenvalues(1)
        end if
    end subroutine analyze_results

    function choose(n, k) result(result)
        ! Binomial coefficient n choose k
        implicit none
        integer, intent(in) :: n, k
        integer(8) :: result
        integer :: i
        
        if (k > n .or. k < 0) then
            result = 0
            return
        end if
        
        if (k == 0 .or. k == n) then
            result = 1
            return  
        end if
        
        result = 1
        do i = 1, min(k, n-k)
            result = result * (n - i + 1) / i
        end do
    end function choose

end program simple_hubbard_model