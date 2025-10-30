module advanced_lattice
    ! Advanced lattice module supporting arbitrary sizes and tilted honeycomb clusters
    ! Handles both rectangular and tilted lattice geometries with proper connectivity
    
    implicit none
    
    ! Lattice geometry types
    type :: lattice_site
        integer :: index                    ! Global site index
        integer :: sublattice              ! A=1, B=2
        double precision :: coord(2)       ! x,y coordinates
        integer :: unit_cell(2)           ! Unit cell indices (ix, iy)
    end type lattice_site
    
    type :: lattice_bond
        integer :: site1, site2           ! Connected site indices
        integer :: bond_type              ! 1=nearest neighbor, 2=next-nearest
        double precision :: vector(2)     ! Bond vector
        double precision :: length        ! Bond length
    end type lattice_bond
    
    type :: honeycomb_lattice
        ! Lattice parameters
        integer :: Lx, Ly                 ! Lattice dimensions
        logical :: tilted                 ! Rectangular (false) vs tilted (true)
        logical :: periodic_x, periodic_y ! Boundary conditions
        
        ! Geometry
        integer :: n_sites                ! Total number of sites
        integer :: n_unit_cells          ! Number of unit cells
        integer :: n_bonds_nn            ! Number of nearest-neighbor bonds
        integer :: n_bonds_nnn           ! Number of next-nearest-neighbor bonds
        
        ! Arrays
        type(lattice_site), allocatable :: sites(:)
        type(lattice_bond), allocatable :: bonds_nn(:)
        type(lattice_bond), allocatable :: bonds_nnn(:)
        
        ! Lattice vectors
        double precision :: a1(2), a2(2)  ! Primitive lattice vectors
        double precision :: b1(2), b2(2)  ! Reciprocal lattice vectors
        
        ! Site lookup tables
        integer, allocatable :: site_map(:,:)  ! site_map(ix,iy) -> site index
        integer, allocatable :: A_sites(:)     ! List of A sublattice sites
        integer, allocatable :: B_sites(:)     ! List of B sublattice sites
    end type honeycomb_lattice
    
contains

    subroutine create_honeycomb_lattice(lattice, Lx, Ly, tilted, periodic_x, periodic_y)
        ! Create honeycomb lattice with specified parameters
        implicit none
        type(honeycomb_lattice), intent(out) :: lattice
        integer, intent(in) :: Lx, Ly
        logical, intent(in) :: tilted
        logical, intent(in), optional :: periodic_x, periodic_y
        
        lattice%Lx = Lx
        lattice%Ly = Ly
        lattice%tilted = tilted
        lattice%periodic_x = .true.
        lattice%periodic_y = .true.
        
        if (present(periodic_x)) lattice%periodic_x = periodic_x
        if (present(periodic_y)) lattice%periodic_y = periodic_y
        
        lattice%n_unit_cells = Lx * Ly
        lattice%n_sites = 2 * lattice%n_unit_cells
        
        ! Set up lattice vectors
        if (tilted) then
            call setup_tilted_vectors(lattice)
        else
            call setup_rectangular_vectors(lattice)
        end if
        
        ! Create sites and bonds
        call create_lattice_sites(lattice)
        call create_lattice_bonds(lattice)
        call setup_site_lookup_tables(lattice)
        
        write(*,'(A,I0,A,I0,A,L1,A)') 'Created honeycomb lattice: ', Lx, 'x', Ly, &
            ', tilted=', tilted, ', with ', lattice%n_sites, ' sites'
        write(*,'(A,I0,A,I0,A)') 'Bonds: ', lattice%n_bonds_nn, ' NN, ', &
            lattice%n_bonds_nnn, ' NNN'
    end subroutine create_honeycomb_lattice

    subroutine setup_rectangular_vectors(lattice)
        ! Set up lattice vectors for rectangular honeycomb lattice
        implicit none
        type(honeycomb_lattice), intent(inout) :: lattice
        
        double precision, parameter :: sqrt3 = sqrt(3.0d0)
        
        ! Primitive lattice vectors for rectangular honeycomb
        lattice%a1 = [3.0d0, 0.0d0]
        lattice%a2 = [0.0d0, sqrt3]
        
        ! Reciprocal lattice vectors
        lattice%b1 = [2.0d0*acos(-1.0d0)/3.0d0, 0.0d0]  ! 2π/3
        lattice%b2 = [0.0d0, 2.0d0*acos(-1.0d0)/sqrt3]  ! 2π/√3
    end subroutine setup_rectangular_vectors

    subroutine setup_tilted_vectors(lattice)
        ! Set up lattice vectors for tilted honeycomb lattice
        implicit none
        type(honeycomb_lattice), intent(inout) :: lattice
        
        double precision, parameter :: sqrt3 = sqrt(3.0d0)
        
        ! Primitive lattice vectors for tilted honeycomb
        lattice%a1 = [1.5d0, sqrt3/2.0d0]
        lattice%a2 = [1.5d0, -sqrt3/2.0d0]
        
        ! Reciprocal lattice vectors
        lattice%b1 = [2.0d0*acos(-1.0d0)/3.0d0, 2.0d0*acos(-1.0d0)/sqrt3]
        lattice%b2 = [2.0d0*acos(-1.0d0)/3.0d0, -2.0d0*acos(-1.0d0)/sqrt3]
    end subroutine setup_tilted_vectors

    subroutine create_lattice_sites(lattice)
        ! Create all lattice sites with coordinates and properties
        implicit none
        type(honeycomb_lattice), intent(inout) :: lattice
        
        integer :: ix, iy, site_idx
        double precision :: r_cell(2), r_A(2), r_B(2)
        double precision, parameter :: sqrt3 = sqrt(3.0d0)
        
        allocate(lattice%sites(lattice%n_sites))
        allocate(lattice%site_map(lattice%Lx, lattice%Ly))
        
        site_idx = 0
        
        do iy = 0, lattice%Ly - 1
            do ix = 0, lattice%Lx - 1
                ! Unit cell position
                r_cell = ix * lattice%a1 + iy * lattice%a2
                
                if (lattice%tilted) then
                    ! Tilted lattice site positions within unit cell
                    r_A = r_cell + [0.0d0, 0.0d0]
                    r_B = r_cell + [1.0d0, 0.0d0]
                else
                    ! Rectangular lattice site positions within unit cell
                    r_A = r_cell + [0.0d0, 0.0d0]
                    r_B = r_cell + [1.0d0, 1.0d0/sqrt3]
                end if
                
                ! A sublattice site
                site_idx = site_idx + 1
                lattice%sites(site_idx)%index = site_idx
                lattice%sites(site_idx)%sublattice = 1
                lattice%sites(site_idx)%coord = r_A
                lattice%sites(site_idx)%unit_cell = [ix, iy]
                
                ! B sublattice site
                site_idx = site_idx + 1
                lattice%sites(site_idx)%index = site_idx
                lattice%sites(site_idx)%sublattice = 2
                lattice%sites(site_idx)%coord = r_B
                lattice%sites(site_idx)%unit_cell = [ix, iy]
                
                ! Store mapping from unit cell to site indices
                lattice%site_map(ix + 1, iy + 1) = site_idx - 1  ! A site index
            end do
        end do
    end subroutine create_lattice_sites

    subroutine create_lattice_bonds(lattice)
        ! Create nearest-neighbor and next-nearest-neighbor bonds
        implicit none
        type(honeycomb_lattice), intent(inout) :: lattice
        
        integer :: ix, iy, site_A, site_B, bond_idx
        integer :: ix_nn, iy_nn, neighbor_bonds
        double precision :: bond_vector(2), bond_length
        integer, allocatable :: temp_bonds(:,:)
        double precision, allocatable :: temp_vectors(:,:)
        
        ! Count bonds first
        lattice%n_bonds_nn = 0
        lattice%n_bonds_nnn = 0
        
        ! Temporary storage for bond enumeration
        allocate(temp_bonds(2, 10 * lattice%n_sites))  ! Overallocate
        allocate(temp_vectors(2, 10 * lattice%n_sites))
        
        bond_idx = 0
        
        do iy = 0, lattice%Ly - 1
            do ix = 0, lattice%Lx - 1
                site_A = 2 * (iy * lattice%Lx + ix) + 1  ! A site index
                site_B = site_A + 1                      ! B site index
                
                ! Nearest-neighbor bonds
                call add_nearest_neighbor_bonds(lattice, ix, iy, site_A, site_B, &
                    bond_idx, temp_bonds, temp_vectors)
            end do
        end do
        
        ! Allocate final arrays
        lattice%n_bonds_nn = bond_idx
        allocate(lattice%bonds_nn(lattice%n_bonds_nn))
        
        ! Copy to final arrays
        do bond_idx = 1, lattice%n_bonds_nn
            lattice%bonds_nn(bond_idx)%site1 = temp_bonds(1, bond_idx)
            lattice%bonds_nn(bond_idx)%site2 = temp_bonds(2, bond_idx)
            lattice%bonds_nn(bond_idx)%bond_type = 1
            lattice%bonds_nn(bond_idx)%vector = temp_vectors(:, bond_idx)
            lattice%bonds_nn(bond_idx)%length = norm2(temp_vectors(:, bond_idx))
        end do
        
        deallocate(temp_bonds, temp_vectors)
    end subroutine create_lattice_bonds

    subroutine add_nearest_neighbor_bonds(lattice, ix, iy, site_A, site_B, bond_count, &
        temp_bonds, temp_vectors)
        ! Add nearest-neighbor bonds for a unit cell
        implicit none
        type(honeycomb_lattice), intent(in) :: lattice
        integer, intent(in) :: ix, iy, site_A, site_B
        integer, intent(inout) :: bond_count
        integer, intent(inout) :: temp_bonds(:,:)
        double precision, intent(inout) :: temp_vectors(:,:)
        
        integer :: ix_neighbor, iy_neighbor, neighbor_site
        double precision :: bond_vector(2)
        
        if (lattice%tilted) then
            ! Tilted lattice connections
            call add_tilted_nn_bonds(lattice, ix, iy, site_A, site_B, bond_count, &
                temp_bonds, temp_vectors)
        else
            ! Rectangular lattice connections
            call add_rectangular_nn_bonds(lattice, ix, iy, site_A, site_B, bond_count, &
                temp_bonds, temp_vectors)
        end if
    end subroutine add_nearest_neighbor_bonds

    subroutine add_rectangular_nn_bonds(lattice, ix, iy, site_A, site_B, bond_count, &
        temp_bonds, temp_vectors)
        ! Add nearest-neighbor bonds for rectangular lattice
        implicit none
        type(honeycomb_lattice), intent(in) :: lattice
        integer, intent(in) :: ix, iy, site_A, site_B
        integer, intent(inout) :: bond_count
        integer, intent(inout) :: temp_bonds(:,:)
        double precision, intent(inout) :: temp_vectors(:,:)
        
        integer :: ix_nn, iy_nn, site_B_nn
        double precision :: bond_vector(2)
        
        ! Bond 1: A(ix,iy) -- B(ix,iy) (vertical bond within unit cell)
        bond_count = bond_count + 1
        temp_bonds(1, bond_count) = site_A
        temp_bonds(2, bond_count) = site_B
        temp_vectors(:, bond_count) = lattice%sites(site_B)%coord - lattice%sites(site_A)%coord
        
        ! Bond 2: A(ix,iy) -- B(ix-1,iy) (horizontal bond to left)
        ix_nn = modulo(ix - 1, lattice%Lx)
        iy_nn = iy
        if (ix > 0 .or. lattice%periodic_x) then
            site_B_nn = 2 * (iy_nn * lattice%Lx + ix_nn) + 2
            bond_count = bond_count + 1
            temp_bonds(1, bond_count) = site_A
            temp_bonds(2, bond_count) = site_B_nn
            bond_vector = lattice%sites(site_B_nn)%coord - lattice%sites(site_A)%coord
            if (ix == 0 .and. lattice%periodic_x) then
                bond_vector(1) = bond_vector(1) - lattice%Lx * lattice%a1(1)
            end if
            temp_vectors(:, bond_count) = bond_vector
        end if
        
        ! Bond 3: A(ix,iy) -- B(ix,iy-1) (diagonal bond to bottom)
        ix_nn = ix
        iy_nn = modulo(iy - 1, lattice%Ly)
        if (iy > 0 .or. lattice%periodic_y) then
            site_B_nn = 2 * (iy_nn * lattice%Lx + ix_nn) + 2
            bond_count = bond_count + 1
            temp_bonds(1, bond_count) = site_A
            temp_bonds(2, bond_count) = site_B_nn
            bond_vector = lattice%sites(site_B_nn)%coord - lattice%sites(site_A)%coord
            if (iy == 0 .and. lattice%periodic_y) then
                bond_vector(2) = bond_vector(2) - lattice%Ly * lattice%a2(2)
            end if
            temp_vectors(:, bond_count) = bond_vector
        end if
    end subroutine add_rectangular_nn_bonds

    subroutine add_tilted_nn_bonds(lattice, ix, iy, site_A, site_B, bond_count, &
        temp_bonds, temp_vectors)
        ! Add nearest-neighbor bonds for tilted lattice
        implicit none
        type(honeycomb_lattice), intent(in) :: lattice
        integer, intent(in) :: ix, iy, site_A, site_B
        integer, intent(inout) :: bond_count
        integer, intent(inout) :: temp_bonds(:,:)
        double precision, intent(inout) :: temp_vectors(:,:)
        
        integer :: ix_nn, iy_nn, site_B_nn
        double precision :: bond_vector(2)
        
        ! Bond 1: A(ix,iy) -- B(ix,iy) (horizontal bond within unit cell)
        bond_count = bond_count + 1
        temp_bonds(1, bond_count) = site_A
        temp_bonds(2, bond_count) = site_B
        temp_vectors(:, bond_count) = [1.0d0, 0.0d0]
        
        ! Bond 2: A(ix,iy) -- B(ix,iy-1) (bond to lower-right)
        ix_nn = ix
        iy_nn = modulo(iy - 1, lattice%Ly)
        if (iy > 0 .or. lattice%periodic_y) then
            site_B_nn = 2 * (iy_nn * lattice%Lx + ix_nn) + 2
            bond_count = bond_count + 1
            temp_bonds(1, bond_count) = site_A
            temp_bonds(2, bond_count) = site_B_nn
            temp_vectors(:, bond_count) = [0.5d0, sqrt(3.0d0)/2.0d0]
        end if
        
        ! Bond 3: A(ix,iy) -- B(ix-1,iy-1) (bond to lower-left)
        ix_nn = modulo(ix - 1, lattice%Lx)
        iy_nn = modulo(iy - 1, lattice%Ly)
        if ((ix > 0 .or. lattice%periodic_x) .and. (iy > 0 .or. lattice%periodic_y)) then
            site_B_nn = 2 * (iy_nn * lattice%Lx + ix_nn) + 2
            bond_count = bond_count + 1
            temp_bonds(1, bond_count) = site_A
            temp_bonds(2, bond_count) = site_B_nn
            temp_vectors(:, bond_count) = [-0.5d0, sqrt(3.0d0)/2.0d0]
        end if
    end subroutine add_tilted_nn_bonds

    subroutine setup_site_lookup_tables(lattice)
        ! Create lookup tables for fast site access
        implicit none
        type(honeycomb_lattice), intent(inout) :: lattice
        
        integer :: i, n_A, n_B
        
        ! Count A and B sites
        n_A = 0
        n_B = 0
        do i = 1, lattice%n_sites
            if (lattice%sites(i)%sublattice == 1) then
                n_A = n_A + 1
            else
                n_B = n_B + 1
            end if
        end do
        
        allocate(lattice%A_sites(n_A))
        allocate(lattice%B_sites(n_B))
        
        ! Fill lookup tables
        n_A = 0
        n_B = 0
        do i = 1, lattice%n_sites
            if (lattice%sites(i)%sublattice == 1) then
                n_A = n_A + 1
                lattice%A_sites(n_A) = i
            else
                n_B = n_B + 1
                lattice%B_sites(n_B) = i
            end if
        end do
    end subroutine setup_site_lookup_tables

    subroutine print_lattice_info(lattice)
        ! Print detailed lattice information
        implicit none
        type(honeycomb_lattice), intent(in) :: lattice
        
        integer :: i
        
        write(*,*) '=== Lattice Information ==='
        write(*,'(A,I0,A,I0)') 'Dimensions: ', lattice%Lx, ' x ', lattice%Ly
        write(*,'(A,L1)') 'Tilted: ', lattice%tilted
        write(*,'(A,I0)') 'Total sites: ', lattice%n_sites
        write(*,'(A,I0)') 'NN bonds: ', lattice%n_bonds_nn
        
        write(*,*) 'Lattice vectors:'
        write(*,'(A,2F8.4)') 'a1 = ', lattice%a1
        write(*,'(A,2F8.4)') 'a2 = ', lattice%a2
        
        write(*,*) 'Sample site coordinates:'
        do i = 1, min(6, lattice%n_sites)
            write(*,'(A,I2,A,I1,A,2F8.4)') 'Site ', lattice%sites(i)%index, &
                ' (sublattice ', lattice%sites(i)%sublattice, '): ', lattice%sites(i)%coord
        end do
        
        write(*,*) 'Sample bonds:'
        do i = 1, min(6, lattice%n_bonds_nn)
            write(*,'(A,I2,A,I0,A,I0,A,F6.3)') 'Bond ', i, ': sites ', &
                lattice%bonds_nn(i)%site1, '--', lattice%bonds_nn(i)%site2, &
                ', length ', lattice%bonds_nn(i)%length
        end do
        write(*,*)
    end subroutine print_lattice_info

    subroutine destroy_lattice(lattice)
        ! Clean up lattice data structures
        implicit none
        type(honeycomb_lattice), intent(inout) :: lattice
        
        if (allocated(lattice%sites)) deallocate(lattice%sites)
        if (allocated(lattice%bonds_nn)) deallocate(lattice%bonds_nn)
        if (allocated(lattice%bonds_nnn)) deallocate(lattice%bonds_nnn)
        if (allocated(lattice%site_map)) deallocate(lattice%site_map)
        if (allocated(lattice%A_sites)) deallocate(lattice%A_sites)
        if (allocated(lattice%B_sites)) deallocate(lattice%B_sites)
    end subroutine destroy_lattice

end module advanced_lattice