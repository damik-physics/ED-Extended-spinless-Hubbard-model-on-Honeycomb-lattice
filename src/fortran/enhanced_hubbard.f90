program enhanced_hubbard_model
    ! Enhanced honeycomb Hubbard model with advanced features:
    ! - Sparse matrix storage and Lanczos diagonalization
    ! - Arbitrary lattice sizes and tilted geometries  
    ! - Translation invariance and C6v point group symmetries
    
    use advanced_lattice
    use sparse_matrix
    use symmetry_operations
    use hubbard_utilities
    
    implicit none
    
    ! Physical parameters
    integer, parameter :: Lx_default = 3, Ly_default = 3
    double precision, parameter :: t = 1.0d0
    double precision, parameter :: V1 = 0.8d0
    double precision, parameter :: filling = 0.5d0
    
    ! Computational parameters
    logical, parameter :: use_sparse = .true.
    logical, parameter :: use_tilted = .false.
    logical, parameter :: use_symmetries = .true.
    integer, parameter :: max_lanczos_iter = 200
    double precision, parameter :: lanczos_tol = 1.0d-10
    
    ! System variables
    type(honeycomb_lattice) :: lattice
    integer :: Lx, Ly, n_sites, n_particles, n_evals
    
    ! Basis and Hamiltonian
    integer(8) :: n_basis_full, n_basis_sym
    integer(8), allocatable :: basis_full(:), basis_sym(:)
    double precision, allocatable :: hamiltonian_dense(:,:)
    type(sparse_csr) :: hamiltonian_sparse
    
    ! Symmetry
    type(momentum_sector), allocatable :: k_sectors(:)
    integer :: n_k_sectors, current_k_sector, irrep_label
    
    ! Results
    double precision, allocatable :: eigenvalues(:), eigenvectors(:,:)
    
    ! User input
    character(len=32) :: arg
    logical :: interactive_mode
    
    ! Initialize
    call print_header()
    
    ! Check for command line arguments
    interactive_mode = .true.
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg)
        if (trim(arg) == 'auto') interactive_mode = .false.
    end if
    
    ! Get parameters
    if (interactive_mode) then
        call get_user_parameters(Lx, Ly)
    else
        Lx = Lx_default
        Ly = Ly_default
        write(*,'(A,I0,A,I0,A)') 'Using default parameters: ', Lx, 'x', Ly, ' lattice'
    end if
    
    n_sites = 2 * Lx * Ly
    n_particles = nint(filling * n_sites)
    n_evals = min(10, int(choose_8(n_sites, n_particles)))
    
    ! Create lattice
    call create_honeycomb_lattice(lattice, Lx, Ly, use_tilted, .true., .true.)
    call print_lattice_info(lattice)
    call print_physics_summary(n_sites, n_particles, t, V1)
    
    ! Generate full many-body basis
    write(*,*) 'Generating many-body basis...'
    call generate_full_basis(n_sites, n_particles, basis_full, n_basis_full)
    
    if (use_symmetries) then
        ! Create momentum sectors
        call create_momentum_sectors(lattice, k_sectors)
        n_k_sectors = size(k_sectors)
        
        ! Loop over momentum sectors
        do current_k_sector = 1, min(3, n_k_sectors)  ! Limit for demonstration
            write(*,*)
            write(*,'(A,I0,A,I0,A,I0,A,I0,A)') 'Processing momentum sector ', &
                current_k_sector, '/', n_k_sectors, ' (k = ', &
                k_sectors(current_k_sector)%kx, ',', k_sectors(current_k_sector)%ky, ')'
            
            ! Construct symmetrized basis for this momentum sector
            irrep_label = 1  ! A1 irrep for simplicity
            call construct_symmetrized_basis(lattice, basis_full, n_particles, &
                k_sectors(current_k_sector), irrep_label, basis_sym, n_basis_sym)
            
            if (n_basis_sym == 0) then
                write(*,*) 'No states in this symmetry sector, skipping...'
                cycle
            end if
            
            ! Build and diagonalize Hamiltonian for this sector
            call solve_hamiltonian_sector(lattice, basis_sym, n_basis_sym, &
                k_sectors(current_k_sector), eigenvalues, eigenvectors)
            
            ! Analyze results
            call analyze_sector_results(k_sectors(current_k_sector), eigenvalues, eigenvectors)
            
            if (allocated(basis_sym)) deallocate(basis_sym)
        end do
        
    else
        ! No symmetries - work with full basis
        write(*,*) 'Working with full basis (no symmetries)'
        n_basis_sym = n_basis_full
        allocate(basis_sym(n_basis_sym))
        basis_sym = basis_full
        
        ! Create dummy momentum sector
        allocate(k_sectors(1))
        k_sectors(1)%kx = 0
        k_sectors(1)%ky = 0
        k_sectors(1)%k_vector = [0.0d0, 0.0d0]
        
        call solve_hamiltonian_sector(lattice, basis_sym, n_basis_sym, &
            k_sectors(1), eigenvalues, eigenvectors)
        call analyze_sector_results(k_sectors(1), eigenvalues, eigenvectors)
    end if
    
    ! Cleanup
    call destroy_lattice(lattice)
    if (allocated(basis_full)) deallocate(basis_full)
    if (allocated(basis_sym)) deallocate(basis_sym)
    if (allocated(k_sectors)) deallocate(k_sectors)
    if (allocated(eigenvalues)) deallocate(eigenvalues)
    if (allocated(eigenvectors)) deallocate(eigenvectors)
    
    write(*,*)
    write(*,*) 'Enhanced Hubbard model calculation completed successfully!'
    
contains

    subroutine print_header()
        write(*,*) '======================================================='
        write(*,*) '      Enhanced Honeycomb Hubbard Model Solver'
        write(*,*) '======================================================='
        write(*,*) 'Features:'
        write(*,*) '  ✓ Arbitrary lattice sizes and tilted geometries'
        if (use_sparse) then
            write(*,*) '  ✓ Sparse matrix storage with Lanczos diagonalization'
        else
            write(*,*) '  ✓ Dense matrix diagonalization'
        end if
        if (use_symmetries) then
            write(*,*) '  ✓ Translation invariance and C6v point group symmetries'
        else
            write(*,*) '  - No symmetries (full Hilbert space)'
        end if
        write(*,*) '======================================================='
        write(*,*)
    end subroutine print_header

    subroutine get_user_parameters(Lx_out, Ly_out)
        integer, intent(out) :: Lx_out, Ly_out
        
        write(*,*) 'Enter lattice parameters:'
        write(*,'(A)', advance='no') 'Lx (lattice width, default 3): '
        read(*,*) Lx_out
        write(*,'(A)', advance='no') 'Ly (lattice height, default 3): '
        read(*,*) Ly_out
        
        if (Lx_out <= 0) Lx_out = Lx_default
        if (Ly_out <= 0) Ly_out = Ly_default
        
        ! Warn about computational cost
        if (Lx_out * Ly_out > 16) then
            write(*,*)
            write(*,*) 'WARNING: Large lattice may require significant computation time!'
            write(*,'(A,I0,A)') 'Estimated basis dimension: ~', &
                int(choose_8(2*Lx_out*Ly_out, nint(filling*2*Lx_out*Ly_out))/1000000), ' million'
            write(*,*) 'Continue? (y/n)'
            read(*,*) arg
            if (trim(arg) /= 'y' .and. trim(arg) /= 'Y') then
                write(*,*) 'Exiting...'
                stop
            end if
        end if
    end subroutine get_user_parameters

    subroutine generate_full_basis(n_sites, n_particles, basis, n_basis)
        integer, intent(in) :: n_sites, n_particles
        integer(8), allocatable, intent(out) :: basis(:)
        integer(8), intent(out) :: n_basis
        
        integer(8) :: i, count, state
        integer :: j, n_occupied
        
        n_basis = choose_8(n_sites, n_particles)
        allocate(basis(n_basis))
        
        count = 0
        do i = 0, 2_8**n_sites - 1
            n_occupied = 0
            state = i
            do j = 0, n_sites - 1
                if (btest(state, j)) n_occupied = n_occupied + 1
            end do
            
            if (n_occupied == n_particles) then
                count = count + 1
                if (count <= n_basis) basis(count) = i
            end if
        end do
        
        write(*,'(A,I0,A)') 'Generated full basis with ', count, ' states'
    end subroutine generate_full_basis

    subroutine solve_hamiltonian_sector(lattice, basis, n_basis, k_sector, evals, evecs)
        type(honeycomb_lattice), intent(in) :: lattice
        integer(8), intent(in) :: basis(:)
        integer(8), intent(in) :: n_basis
        type(momentum_sector), intent(in) :: k_sector
        double precision, allocatable, intent(out) :: evals(:), evecs(:,:)
        
        if (use_sparse .and. n_basis > 100) then
            call solve_sparse_hamiltonian(lattice, basis, n_basis, k_sector, evals, evecs)
        else
            call solve_dense_hamiltonian(lattice, basis, n_basis, k_sector, evals, evecs)
        end if
    end subroutine solve_hamiltonian_sector

    subroutine solve_dense_hamiltonian(lattice, basis, n_basis, k_sector, evals, evecs)
        type(honeycomb_lattice), intent(in) :: lattice
        integer(8), intent(in) :: basis(:)
        integer(8), intent(in) :: n_basis
        type(momentum_sector), intent(in) :: k_sector
        double precision, allocatable, intent(out) :: evals(:), evecs(:,:)
        
        double precision, allocatable :: hamiltonian(:,:)
        integer :: info, lwork, n
        double precision, allocatable :: work(:)
        
        n = int(n_basis)
        allocate(hamiltonian(n, n))
        allocate(evals(n))
        allocate(evecs(n, n))
        
        write(*,*) 'Constructing dense Hamiltonian matrix...'
        call construct_dense_hamiltonian(lattice, basis, n_basis, hamiltonian)
        
        write(*,*) 'Diagonalizing using LAPACK...'
        evecs = hamiltonian
        
        ! LAPACK diagonalization
        lwork = -1
        allocate(work(1))
        call dsyev('V', 'U', n, evecs, n, evals, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V', 'U', n, evecs, n, evals, work, lwork, info)
        
        if (info /= 0) then
            write(*,*) 'ERROR: Dense diagonalization failed'
            stop
        end if
        
        write(*,'(A,I0,A)') 'Dense diagonalization completed for ', n, ' states'
        deallocate(hamiltonian, work)
    end subroutine solve_dense_hamiltonian

    subroutine solve_sparse_hamiltonian(lattice, basis, n_basis, k_sector, evals, evecs)
        type(honeycomb_lattice), intent(in) :: lattice
        integer(8), intent(in) :: basis(:)
        integer(8), intent(in) :: n_basis
        type(momentum_sector), intent(in) :: k_sector
        double precision, allocatable, intent(out) :: evals(:), evecs(:,:)
        
        double precision, allocatable :: hamiltonian_dense(:,:)
        type(sparse_csr) :: hamiltonian_sparse
        integer :: n_evals_lanczos
        
        n_evals_lanczos = min(n_evals, int(n_basis))
        
        write(*,*) 'Constructing sparse Hamiltonian matrix...'
        
        ! First construct dense matrix, then convert to sparse
        allocate(hamiltonian_dense(n_basis, n_basis))
        call construct_dense_hamiltonian(lattice, basis, n_basis, hamiltonian_dense)
        call dense_to_sparse_csr(hamiltonian_dense, hamiltonian_sparse)
        deallocate(hamiltonian_dense)
        
        write(*,'(A,I0,A)') 'Starting Lanczos diagonalization for ', n_evals_lanczos, ' eigenvalues...'
        call lanczos_diagonalize(hamiltonian_sparse, evals, evecs, n_evals_lanczos, &
            max_lanczos_iter, lanczos_tol)
        
        call sparse_csr_destroy(hamiltonian_sparse)
    end subroutine solve_sparse_hamiltonian

    subroutine construct_dense_hamiltonian(lattice, basis, n_basis, hamiltonian)
        type(honeycomb_lattice), intent(in) :: lattice
        integer(8), intent(in) :: basis(:)
        integer(8), intent(in) :: n_basis
        double precision, intent(out) :: hamiltonian(:,:)
        
        integer(8) :: i, j, new_state, pos
        integer :: bond, site1, site2
        double precision :: matrix_element
        
        hamiltonian = 0.0d0
        
        write(*,'(A)', advance='no') 'Constructing matrix'
        
        do i = 1, n_basis
            if (mod(i, max(1_8, n_basis/10)) == 0) then
                write(*,'(A)', advance='no') '.'
            end if
            
            ! Diagonal terms: interaction energy
            do bond = 1, lattice%n_bonds_nn
                site1 = lattice%bonds_nn(bond)%site1
                site2 = lattice%bonds_nn(bond)%site2
                
                if (btest(basis(i), site1-1) .and. btest(basis(i), site2-1)) then
                    hamiltonian(i, i) = hamiltonian(i, i) + V1
                end if
            end do
            
            ! Off-diagonal terms: hopping
            do bond = 1, lattice%n_bonds_nn
                site1 = lattice%bonds_nn(bond)%site1
                site2 = lattice%bonds_nn(bond)%site2
                
                ! Forward hopping
                call try_hopping_advanced(basis(i), site1, site2, new_state, matrix_element)
                if (matrix_element /= 0.0d0) then
                    call find_state_index_advanced(basis, new_state, pos)
                    if (pos > 0) hamiltonian(i, pos) = hamiltonian(i, pos) + matrix_element
                end if
                
                ! Backward hopping
                call try_hopping_advanced(basis(i), site2, site1, new_state, matrix_element)
                if (matrix_element /= 0.0d0) then
                    call find_state_index_advanced(basis, new_state, pos)
                    if (pos > 0) hamiltonian(i, pos) = hamiltonian(i, pos) + matrix_element
                end if
            end do
        end do
        
        write(*,*) ' done'
    end subroutine construct_dense_hamiltonian

    subroutine try_hopping_advanced(state, from, to, new_state, matrix_element)
        integer(8), intent(in) :: state
        integer, intent(in) :: from, to
        integer(8), intent(out) :: new_state
        double precision, intent(out) :: matrix_element
        
        integer :: sign_factor, bit_count, i
        
        matrix_element = 0.0d0
        new_state = 0
        
        if (btest(state, from-1) .and. .not. btest(state, to-1)) then
            new_state = ibset(ibclr(state, from-1), to-1)
            
            sign_factor = 1
            bit_count = 0
            do i = min(from-1, to-1), max(from-1, to-1) - 1
                if (btest(state, i)) bit_count = bit_count + 1
            end do
            if (mod(bit_count, 2) == 1) sign_factor = -sign_factor
            
            matrix_element = -t * sign_factor
        end if
    end subroutine try_hopping_advanced

    subroutine find_state_index_advanced(basis, target_state, index)
        integer(8), intent(in) :: basis(:), target_state
        integer(8), intent(out) :: index
        integer(8) :: left, right, mid
        
        left = 1
        right = size(basis, kind=8)
        index = 0
        
        do while (left <= right)
            mid = (left + right) / 2
            if (basis(mid) == target_state) then
                index = mid
                return
            else if (basis(mid) < target_state) then
                left = mid + 1
            else
                right = mid - 1
            end if
        end do
    end subroutine find_state_index_advanced

    subroutine analyze_sector_results(k_sector, evals, evecs)
        type(momentum_sector), intent(in) :: k_sector
        double precision, intent(in) :: evals(:), evecs(:,:)
        
        integer :: n_evals_show, i
        double precision :: energy_per_site, gap
        
        write(*,*)
        write(*,'(A,I0,A,I0,A)') '=== Results for k = (', k_sector%kx, ',', k_sector%ky, ') ==='
        write(*,'(A,F14.8)') 'Ground state energy: ', evals(1)
        
        energy_per_site = evals(1) / n_sites
        write(*,'(A,F14.8)') 'Energy per site: ', energy_per_site
        
        n_evals_show = min(5, size(evals))
        write(*,*)
        write(*,*) 'Lowest eigenvalues:'
        do i = 1, n_evals_show
            write(*,'(I3,F14.8)') i, evals(i)
        end do
        
        if (size(evals) > 1) then
            gap = evals(2) - evals(1)
            write(*,*)
            write(*,'(A,F12.8)') 'Energy gap: ', gap
        end if
        
        ! Save results for this sector
        if (k_sector%kx == 0 .and. k_sector%ky == 0) then
            call save_results('hubbard_results.dat', evals, evecs(:,1))
        end if
    end subroutine analyze_sector_results

    function choose_8(n, k) result(result)
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
    end function choose_8

end program enhanced_hubbard_model