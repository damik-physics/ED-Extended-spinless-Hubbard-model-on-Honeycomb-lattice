module variables 

    implicit none 
    
    ! System size and geometry
    integer                      :: sites, particles, orbsize, l1, l2, k1, k2, k1_max, k2_max
    integer(kind=8)              :: dim

    ! Hamiltonian and basis
    integer(kind=8)              :: nOff, nDi, nDiOff, nDi_dp, nnz 
    integer(kind=8), allocatable :: basis_states(:), abasis(:), bbasis(:), mombasis(:), period(:), orbits2D(:,:,:)
    integer(kind=8), allocatable :: occ(:,:), hamOff(:,:), hamDi(:,:), rc(:,:), rcOff(:,:), rcDi(:)
    double precision, allocatable:: ham(:), ham_dp(:,:), hamOff_dp(:), hamDi_dp(:)
    double complex, allocatable  :: ham_dc(:), hamOff_dc(:), hamDi_dc(:,:), hamDi_dc_2(:), hamDi_off_dc(:)

    ! Eigenvalues and eigenvectors
    double precision, allocatable:: energies(:), eigstate(:,:), evals(:), gs(:), norm(:), norm2D(:,:)
    double complex, allocatable  :: eigstate_dc(:,:), gs_dc(:)

    ! Bonds, sites, and lattice
    integer                      :: nnBonds, nnnBonds
    integer, allocatable         :: bsites(:,:), hexsites(:,:), xy(:,:), xyA(:,:), xyB(:,:), latticeVecs(:)
    integer, allocatable         :: Alattice(:,:), Blattice(:,:), AsitesBonds(:,:), BsitesBonds(:,:)
    integer, allocatable         :: xTransl(:,:), yTransl(:,:)
    double precision, allocatable :: nnnVec(:,:)

    ! Symmetry and group theory
    integer        , allocatable :: parities(:), dplcts(:), phases(:), reflections(:), refl(:,:), c6(:)
    double precision             :: mir(6), rot(5), id
    double complex, allocatable  :: phases2D(:,:,:)

    ! Calculation control and threading
    integer                      :: nDeg, unit, refsite, conf, nest, cntrA, cntrB
    integer                      :: nev, ncv, ndv, ndv2, nv, nv2, full, nHel, tilt
    integer                      :: v1_thrds = 1, v2_thrds = 1, dis_thrds = 1, thread_num, thread_num_2, num_thrds
    integer, allocatable         :: units(:,:), units_2(:,:,:,:)

    ! File and type
    character                    :: type*1, dir*1000, param_list*1000

    ! Parameters
    double precision             :: v1, v2

end module variables