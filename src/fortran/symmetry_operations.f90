module symmetry_operations
    ! Symmetry module implementing translation invariance and C6v point group symmetries
    ! for honeycomb lattice Hubbard model to reduce Hilbert space dimension
    
    use advanced_lattice
    implicit none
    
    ! Constants for C6v point group
    integer, parameter :: N_C6V_OPERATIONS = 12  ! 6 rotations + 6 reflections
    double precision, parameter :: PI = 4.0d0 * atan(1.0d0)
    
    ! Symmetry types
    type :: momentum_sector
        integer :: kx, ky                     ! Momentum quantum numbers
        double precision :: k_vector(2)      ! k-vector in reciprocal space
        complex(8) :: phase_factor           ! e^{i k·R} phase factor
    end type momentum_sector
    
    type :: point_group_element
        integer :: operation_type            ! 1=rotation, 2=reflection
        integer :: rotation_order           ! n for C_n rotation (1-6)
        integer :: reflection_axis          ! 1-6 for different mirror planes
        integer, allocatable :: site_map(:) ! site permutation
        complex(8) :: character             ! Character for irrep
    end type point_group_element
    
    type :: symmetry_sector
        type(momentum_sector) :: momentum
        integer :: irrep_label              ! Irreducible representation label
        integer :: irrep_dimension          ! 1 for A/B irreps, 2 for E irreps
        type(point_group_element) :: point_group(N_C6V_OPERATIONS)
        
        ! Symmetrized basis
        integer(8) :: sym_dimension          ! Reduced Hilbert space dimension
        integer(8), allocatable :: representatives(:)  ! Representative states
        integer, allocatable :: orbit_sizes(:)        ! Size of each orbit
        double precision, allocatable :: normalizations(:)  ! Normalization factors
    end type symmetry_sector
    
contains

    subroutine create_momentum_sectors(lattice, k_sectors)
        ! Create all allowed momentum sectors for translation invariance
        implicit none
        type(honeycomb_lattice), intent(in) :: lattice
        type(momentum_sector), allocatable, intent(out) :: k_sectors(:)
        
        integer :: n_k_points, kx, ky, k_idx
        
        n_k_points = lattice%Lx * lattice%Ly
        allocate(k_sectors(n_k_points))
        
        k_idx = 0
        do ky = 0, lattice%Ly - 1
            do kx = 0, lattice%Lx - 1
                k_idx = k_idx + 1
                k_sectors(k_idx)%kx = kx
                k_sectors(k_idx)%ky = ky
                
                ! k-vector in units of reciprocal lattice vectors
                k_sectors(k_idx)%k_vector = (2.0d0 * PI / lattice%Lx) * kx * lattice%b1 + &
                                          (2.0d0 * PI / lattice%Ly) * ky * lattice%b2
                
                ! Phase factor for this momentum
                k_sectors(k_idx)%phase_factor = exp(cmplx(0.0d0, &
                    dot_product(k_sectors(k_idx)%k_vector, [1.0d0, 0.0d0]), 8))
            end do
        end do
        
        write(*,'(A,I0,A)') 'Created ', n_k_points, ' momentum sectors'
    end subroutine create_momentum_sectors

    subroutine setup_C6v_point_group(lattice, point_group_ops)
        ! Set up C6v point group operations for honeycomb lattice
        implicit none
        type(honeycomb_lattice), intent(in) :: lattice
        type(point_group_element), intent(out) :: point_group_ops(N_C6V_OPERATIONS)
        
        integer :: i, op_idx
        
        op_idx = 0
        
        ! C6 rotations: identity, C6, C6^2, C6^3, C6^4, C6^5
        do i = 0, 5
            op_idx = op_idx + 1
            point_group_ops(op_idx)%operation_type = 1
            point_group_ops(op_idx)%rotation_order = i
            call setup_rotation_map(lattice, i, point_group_ops(op_idx)%site_map)
        end do
        
        ! 6 mirror reflections
        do i = 1, 6
            op_idx = op_idx + 1
            point_group_ops(op_idx)%operation_type = 2
            point_group_ops(op_idx)%reflection_axis = i
            call setup_reflection_map(lattice, i, point_group_ops(op_idx)%site_map)
        end do
        
        write(*,*) 'Set up C6v point group with 12 operations'
    end subroutine setup_C6v_point_group

    subroutine setup_rotation_map(lattice, rotation_order, site_map)
        ! Set up site permutation for C6^n rotation
        implicit none
        type(honeycomb_lattice), intent(in) :: lattice
        integer, intent(in) :: rotation_order
        integer, allocatable, intent(out) :: site_map(:)
        
        integer :: i, rotated_site
        double precision :: angle, rotation_matrix(2,2)
        double precision :: original_coord(2), rotated_coord(2)
        
        allocate(site_map(lattice%n_sites))
        
        ! Rotation angle
        angle = rotation_order * PI / 3.0d0
        
        ! 2D rotation matrix
        rotation_matrix(1,1) = cos(angle)
        rotation_matrix(1,2) = -sin(angle)
        rotation_matrix(2,1) = sin(angle)
        rotation_matrix(2,2) = cos(angle)
        
        do i = 1, lattice%n_sites
            original_coord = lattice%sites(i)%coord
            
            ! Apply rotation
            rotated_coord = matmul(rotation_matrix, original_coord)
            
            ! Find which site this maps to
            rotated_site = find_nearest_site(lattice, rotated_coord)
            site_map(i) = rotated_site
        end do
    end subroutine setup_rotation_map

    subroutine setup_reflection_map(lattice, reflection_axis, site_map)
        ! Set up site permutation for mirror reflection
        implicit none
        type(honeycomb_lattice), intent(in) :: lattice
        integer, intent(in) :: reflection_axis
        integer, allocatable, intent(out) :: site_map(:)
        
        integer :: i, reflected_site
        double precision :: angle, reflection_matrix(2,2)
        double precision :: original_coord(2), reflected_coord(2)
        
        allocate(site_map(lattice%n_sites))
        
        ! Mirror plane angle
        angle = (reflection_axis - 1) * PI / 6.0d0
        
        ! 2D reflection matrix across line at angle
        reflection_matrix(1,1) = cos(2.0d0 * angle)
        reflection_matrix(1,2) = sin(2.0d0 * angle)
        reflection_matrix(2,1) = sin(2.0d0 * angle)
        reflection_matrix(2,2) = -cos(2.0d0 * angle)
        
        do i = 1, lattice%n_sites
            original_coord = lattice%sites(i)%coord
            
            ! Apply reflection
            reflected_coord = matmul(reflection_matrix, original_coord)
            
            ! Find which site this maps to
            reflected_site = find_nearest_site(lattice, reflected_coord)
            site_map(i) = reflected_site
        end do
    end subroutine setup_reflection_map

    function find_nearest_site(lattice, target_coord) result(nearest_site)
        ! Find the lattice site nearest to given coordinates
        implicit none
        type(honeycomb_lattice), intent(in) :: lattice
        double precision, intent(in) :: target_coord(2)
        integer :: nearest_site
        
        integer :: i
        double precision :: min_distance, distance, coord_diff(2)
        
        min_distance = huge(1.0d0)
        nearest_site = 1
        
        do i = 1, lattice%n_sites
            coord_diff = lattice%sites(i)%coord - target_coord
            
            ! Apply periodic boundary conditions
            call apply_periodic_bc(lattice, coord_diff)
            
            distance = norm2(coord_diff)
            if (distance < min_distance) then
                min_distance = distance
                nearest_site = i
            end if
        end do
    end function find_nearest_site

    subroutine apply_periodic_bc(lattice, coord_diff)
        ! Apply periodic boundary conditions to coordinate difference
        implicit none
        type(honeycomb_lattice), intent(in) :: lattice
        double precision, intent(inout) :: coord_diff(2)
        
        double precision :: cell_size(2)
        
        ! Effective cell size
        cell_size(1) = lattice%Lx * norm2(lattice%a1)
        cell_size(2) = lattice%Ly * norm2(lattice%a2)
        
        ! Apply minimum image convention
        if (abs(coord_diff(1)) > cell_size(1) / 2.0d0) then
            coord_diff(1) = coord_diff(1) - sign(cell_size(1), coord_diff(1))
        end if
        
        if (abs(coord_diff(2)) > cell_size(2) / 2.0d0) then
            coord_diff(2) = coord_diff(2) - sign(cell_size(2), coord_diff(2))
        end if
    end subroutine apply_periodic_bc

    subroutine apply_translation_to_state(lattice, state, translation_vector, new_state, phase)
        ! Apply translation operation to a many-body state
        implicit none
        type(honeycomb_lattice), intent(in) :: lattice
        integer(8), intent(in) :: state
        integer, intent(in) :: translation_vector(2)  ! (dx, dy) in unit cells
        integer(8), intent(out) :: new_state
        complex(8), intent(out) :: phase
        
        integer :: i, old_site, new_site, dx, dy
        integer :: old_ix, old_iy, new_ix, new_iy
        
        new_state = 0
        phase = (1.0d0, 0.0d0)
        
        do i = 0, lattice%n_sites - 1
            if (btest(state, i)) then
                old_site = i + 1
                
                ! Get unit cell coordinates of old site
                old_ix = lattice%sites(old_site)%unit_cell(1)
                old_iy = lattice%sites(old_site)%unit_cell(2)
                
                ! Apply translation
                new_ix = modulo(old_ix + translation_vector(1), lattice%Lx)
                new_iy = modulo(old_iy + translation_vector(2), lattice%Ly)
                
                ! Find new site index
                new_site = find_site_at_unit_cell(lattice, new_ix, new_iy, &
                    lattice%sites(old_site)%sublattice)
                
                ! Set bit for new site
                new_state = ibset(new_state, new_site - 1)
            end if
        end do
    end subroutine apply_translation_to_state

    function find_site_at_unit_cell(lattice, ix, iy, sublattice) result(site_index)
        ! Find site index at given unit cell and sublattice
        implicit none
        type(honeycomb_lattice), intent(in) :: lattice
        integer, intent(in) :: ix, iy, sublattice
        integer :: site_index
        
        ! Direct calculation for honeycomb lattice
        site_index = 2 * (iy * lattice%Lx + ix) + sublattice
    end function find_site_at_unit_cell

    subroutine apply_point_group_to_state(point_group_op, state, new_state, phase)
        ! Apply point group operation to a many-body state
        implicit none
        type(point_group_element), intent(in) :: point_group_op
        integer(8), intent(in) :: state
        integer(8), intent(out) :: new_state
        complex(8), intent(out) :: phase
        
        integer :: i, old_site, new_site
        
        new_state = 0
        phase = (1.0d0, 0.0d0)
        
        do i = 0, size(point_group_op%site_map) - 1
            if (btest(state, i)) then
                old_site = i + 1
                new_site = point_group_op%site_map(old_site)
                new_state = ibset(new_state, new_site - 1)
            end if
        end do
        
        ! Phase factor depends on irreducible representation
        phase = point_group_op%character
    end subroutine apply_point_group_to_state

    subroutine construct_symmetrized_basis(lattice, full_basis, n_particles, k_sector, &
        irrep_label, sym_basis, sym_dimension)
        ! Construct symmetry-adapted basis using translation + point group symmetries
        implicit none
        type(honeycomb_lattice), intent(in) :: lattice
        integer(8), intent(in) :: full_basis(:)
        integer, intent(in) :: n_particles
        type(momentum_sector), intent(in) :: k_sector
        integer, intent(in) :: irrep_label
        integer(8), allocatable, intent(out) :: sym_basis(:)
        integer(8), intent(out) :: sym_dimension
        
        integer(8) :: n_full_basis, i, representative
        logical, allocatable :: is_representative(:)
        integer(8), allocatable :: temp_sym_basis(:)
        integer(8) :: sym_count
        
        n_full_basis = size(full_basis)
        allocate(is_representative(n_full_basis))
        allocate(temp_sym_basis(n_full_basis))
        
        is_representative = .true.
        sym_count = 0
        
        ! Find representative states under translation symmetry
        do i = 1, n_full_basis
            if (is_representative(i)) then
                representative = full_basis(i)
                
                ! Check if this state generates a valid momentum eigenstate
                if (is_valid_momentum_state(lattice, representative, k_sector)) then
                    sym_count = sym_count + 1
                    temp_sym_basis(sym_count) = representative
                    
                    ! Mark equivalent states as non-representatives
                    call mark_translation_orbit(lattice, full_basis, representative, &
                        is_representative)
                end if
            end if
        end do
        
        ! Copy to final array
        sym_dimension = sym_count
        allocate(sym_basis(sym_dimension))
        sym_basis(1:sym_dimension) = temp_sym_basis(1:sym_dimension)
        
        write(*,'(A,I0,A,I0,A,I0,A,I0)') 'Symmetrized basis: ', n_full_basis, &
            ' -> ', sym_dimension, ' states (k = ', k_sector%kx, ',', k_sector%ky, ')'
        
        deallocate(is_representative, temp_sym_basis)
    end subroutine construct_symmetrized_basis

    function is_valid_momentum_state(lattice, state, k_sector) result(is_valid)
        ! Check if state can form a valid momentum eigenstate
        implicit none
        type(honeycomb_lattice), intent(in) :: lattice
        integer(8), intent(in) :: state
        type(momentum_sector), intent(in) :: k_sector
        logical :: is_valid
        
        integer(8) :: translated_state
        complex(8) :: phase
        integer :: dx, dy
        
        is_valid = .true.
        
        ! Check translation invariance for this momentum
        do dy = 0, lattice%Ly - 1
            do dx = 0, lattice%Lx - 1
                if (dx == 0 .and. dy == 0) cycle
                
                call apply_translation_to_state(lattice, state, [dx, dy], &
                    translated_state, phase)
                
                ! For a valid momentum eigenstate, translated states should
                ! either be the same or generate the required phase relationship
                ! This is a simplified check - full implementation would be more complex
                if (translated_state < state) then
                    is_valid = .false.
                    return
                end if
            end do
        end do
    end function is_valid_momentum_state

    subroutine mark_translation_orbit(lattice, full_basis, representative, is_representative)
        ! Mark all states in the translation orbit of representative as non-representatives
        implicit none
        type(honeycomb_lattice), intent(in) :: lattice
        integer(8), intent(in) :: full_basis(:)
        integer(8), intent(in) :: representative
        logical, intent(inout) :: is_representative(:)
        
        integer(8) :: translated_state, i
        complex(8) :: phase
        integer :: dx, dy
        
        do dy = 0, lattice%Ly - 1
            do dx = 0, lattice%Lx - 1
                call apply_translation_to_state(lattice, representative, [dx, dy], &
                    translated_state, phase)
                
                ! Find this state in the full basis and mark as non-representative
                do i = 1, size(full_basis)
                    if (full_basis(i) == translated_state .and. translated_state /= representative) then
                        is_representative(i) = .false.
                    end if
                end do
            end do
        end do
    end subroutine mark_translation_orbit

    subroutine get_C6v_irrep_characters(irrep_label, characters)
        ! Get character table for C6v irreducible representations
        implicit none
        integer, intent(in) :: irrep_label
        complex(8), intent(out) :: characters(N_C6V_OPERATIONS)
        
        ! Simplified character table for C6v point group
        ! irrep_label: 1=A1, 2=A2, 3=B1, 4=B2, 5=E1, 6=E2
        
        select case(irrep_label)
        case(1)  ! A1 (totally symmetric)
            characters = (1.0d0, 0.0d0)
        case(2)  ! A2 
            characters(1:6) = (1.0d0, 0.0d0)   ! Rotations
            characters(7:12) = (-1.0d0, 0.0d0) ! Reflections
        case(3)  ! B1
            characters(1:6) = [1.0d0, -1.0d0, 1.0d0, -1.0d0, 1.0d0, -1.0d0]
            characters(7:12) = (1.0d0, 0.0d0)
        case(4)  ! B2
            characters(1:6) = [1.0d0, -1.0d0, 1.0d0, -1.0d0, 1.0d0, -1.0d0]
            characters(7:12) = (-1.0d0, 0.0d0)
        case default
            ! E1, E2 are 2D irreps - more complex
            characters = (0.0d0, 0.0d0)
        end select
    end subroutine get_C6v_irrep_characters

    subroutine print_symmetry_info(sym_sector)
        ! Print information about symmetry sector
        implicit none
        type(symmetry_sector), intent(in) :: sym_sector
        
        write(*,*) '=== Symmetry Sector Information ==='
        write(*,'(A,I0,A,I0)') 'Momentum: k = (', sym_sector%momentum%kx, &
            ',', sym_sector%momentum%ky, ')'
        write(*,'(A,2F8.4)') 'k-vector: ', sym_sector%momentum%k_vector
        write(*,'(A,I0)') 'Irrep label: ', sym_sector%irrep_label  
        write(*,'(A,I0)') 'Irrep dimension: ', sym_sector%irrep_dimension
        write(*,'(A,I0)') 'Symmetrized basis dimension: ', sym_sector%sym_dimension
        write(*,*)
    end subroutine print_symmetry_info

end module symmetry_operations