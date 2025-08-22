module types
    ! Module containing derived type definitions for simulation parameters, geometry,
    ! system state, threading, output control, diagonalization, and Hamiltonian data structures
    implicit none
    
    type :: sim_params
        ! All input variables that set the simulation parameters

        ! Lattice and symmetry
        integer :: ucx
        integer :: ucy      
        character(len=1) :: bc
        character(len=2) :: irrep
        integer :: p1, p2, p3 ! Reflection eigenvalues of the symmetry group
        character(len=3) :: cluster ! Cluster name, e.g. '18A' for a 18-site cluster
        
        ! Hamiltonian parameters
        integer :: nDis ! Number of disorder configurations
        double precision :: t ! Hopping strength
        double precision :: dis ! Disorder strength
        double precision :: mass ! Sublattice imbalance
        double precision :: filling ! Filling fraction
        double precision :: dv1 ! Step size for v1
        double precision :: v1min ! Minimum value for v1
        double precision :: v1max ! Maximum value for v1
        double precision :: dv2 ! Step size for v2
        double precision :: v2min ! Minimum value for v2
        double precision :: v2max ! Maximum value for v2

        ! Parallelization parameters
        integer :: othrds ! Number of OpenMP threads for Arpack diagonalization
        integer :: mthrds ! Number of MKL threads for MKL diagonalization
        ! Note: othrds and mthrds should be set to 1 if the respective diagonalization method is not used
        
        ! Flags/enumerators 
        integer :: tilted, ti, k0, symm, corr, curr, refbonds, states
        integer :: degflag ! Flag for including degenerate ground states
        double precision :: deg ! Threshold for considering eigenvalues degenerate, i.e. Abs(E_{i}-E_{i-1}) <= deg are degenerate
        integer :: otf ! Flag for generating Hamiltonian on-the-fly
        ! otf = 1: Flag to generate Hamiltonian on-the-fly during vector-matrix multiplication(does not store sparse Hamiltonian in arrays)
        ! degflag = 1: Flag for including degenerate ground states:
        !                0 = assume non-degenerate ground state, 1 = assume truly degenerate ground states, 2 = assume quasidegenerate ground states

        ! Flags for diagonalization methods
        integer :: feast ! Use FEAST algorithm for diagonalization
        integer :: arpack ! Use Sparse Arpack diagonalization
        integer :: mkl ! Use Sparse MKL diagonalization
        integer :: exact ! Use full exact diagonalization for dense matrices
        logical :: rvec ! Flag for calculating eigenvectors
        !    feast  = 1: Use FEAST algorithm for diagonalization
        !    arpack = 1: Use Sparse Arpack diagonalization
        !    mkl    = 1: Use Sparse MKL diagonalization
        !    exact  = 1: Use full exact diagonalization for dense matrices
        integer :: dimthresh ! Threshold for dimension of the Hamiltonian matrix to use sparse diagonalization methods
        integer :: nevext ! Requested number of eigenvalues to calculate 
        integer :: ncv0 ! Requested number of Lanczos vectors to use(see Arpack manual). Typical on the order of 100-200
        integer :: nev0 ! Initial FEAST subspace size
        integer :: nst  ! Requested number of eigenstates to calculate
        integer :: nevmax ! Sets maximum energy for FEAST from previous run(see FEAST manual in MKL documentation)
        integer :: g_fact ! Electron g-factor

    end type sim_params

    type :: diag_params
        
        ! Diagonalization parameters

        integer :: nev  ! Number of eigenvalues to calculate at runtime
        integer :: nest ! Number of eigenstates to calculate at runtime
        integer :: ncv  ! Number of Lanczos convergence vectors to use at runtime 
        integer :: full ! Flag indicating type of diagonalization. 0: sparse diagonalization, 1: dense diagonalization
        integer :: nDeg ! Number of degenerate ground states
        double precision, allocatable :: energies(:), gs(:), evals(:), norm(:), eigstate(:,:), norm2d(:,:)
        double complex, allocatable :: eigstate_dc(:,:), gs_dc(:)
        
    end type diag_params

    type :: thread_params
        ! Parallelization parameters that are not set by the user
        integer :: dis_thrds ! Number of disorder threads
        integer :: v1_thrds ! Number of V1 threads
        integer :: v2_thrds ! Number of V2 threads
        integer :: num_thrds ! Number of threads for current calculations
        integer, allocatable :: units(:,:), units_2(:,:,:,:) ! Unit numbers for output files in parallel calculations
        integer :: thrd_id ! Thread number for V2 calculations
        integer :: thrd_id_2 ! Thread number for V1 calculations
        integer :: ndv1 ! Number of V1 values
        integer :: ndv2 ! Number of V2 values
    end type thread_params

    type output
        ! Output parameters    
        character(len=1000) :: dir ! Output directory
        integer :: unit ! Unit number for output files
        character(len=:), allocatable :: outdir
        ! character(len=512) :: outdir
        integer :: corr ! Flag for calculating single-operator correlation functions
        integer :: curr ! Flag for calculating current-current correlation functions
        integer :: refbonds ! Number of reference bonds for which currents are calculated
        integer :: states ! Flag for saving eigenstates
 
        ! Eigenvalues and eigenvectors
        double precision, allocatable :: energies(:) ! Array of eigenvalues
        double precision, allocatable :: gs(:) ! Ground state energy
        double complex, allocatable :: eigstate(:,:) ! Eigenstates
        double complex, allocatable :: eigstate_dc(:,:) ! Complex eigenstates
        double complex, allocatable :: gs_dc(:) ! Complex ground state energy

        ! Parameters for correlation function calculations
        integer :: refsite ! Reference site for correlation functions


    end type output

    type geometry
        ! Parameters and arrays related to the lattice and basis and its construction, symmetrization and lattice geometry  

        integer(kind=8) :: dim ! Dimension of the basis
        ! Lattice parameters
        integer :: sites ! Number of sites in the lattice
        integer :: particles ! Number of particles in the system
        integer :: orbsize ! Orbital size
        integer :: l1 ! Lattice dimension 1
        integer :: l2 ! Lattice dimension 2
        integer :: k1_max ! Maximum value for k1
        integer :: k2_max ! Maximum value for k2
        integer :: nUC  ! Number of unit cells
        integer :: nHel ! Number of cluster's helices
        integer :: tilt ! Tilt angle for the cluster
        integer :: cntrA ! Counter for A sites
        integer :: cntrB ! Counter for B sites
        ! Basis arrays
        integer(kind=8), allocatable :: basis_states(:) ! Array of basis states
        integer(kind=8), allocatable :: aBasis(:) ! Array of A basis states
        integer(kind=8), allocatable :: bBasis(:) ! Array of B basis states
        integer(kind=8), allocatable :: momBasis(:) ! Array of momentum basis states
        integer(kind=8), allocatable :: period(:) ! Periodicity of the basis states
        integer(kind=8), allocatable :: orbits2D(:,:,:) ! Orbits of basis states in 2D irrep
        double precision, allocatable :: norm(:) ! Normalization factors
        double precision, allocatable :: norm2D(:,:) ! Normalization factors for 2D irrep
        ! Bonds, sites, and lattice
        integer                      :: nnBonds, nnnBonds
        integer, allocatable         :: bsites(:,:), hexsites(:,:), xy(:,:), xyA(:,:), xyB(:,:), latticevecs(:)
        integer, allocatable         :: Alattice(:,:), Blattice(:,:), AsitesBonds(:,:), BsitesBonds(:,:)
        integer, allocatable         :: xTransl(:,:), yTransl(:,:)
        double precision, allocatable :: nnnVec(:,:)

        ! Symmetry and group theory
        integer        , allocatable :: parities(:), dplcts(:), phases(:), reflections(:), refl(:,:), c6(:), sitecoord(:,:)

        double precision             :: mir(6), rot(5), id
        double complex, allocatable  :: phases2D(:,:,:)


    end type geometry

    type hamiltonian_params
        ! Hamiltonian arrays
        integer(kind=8) :: nDiOff ! Number of diagonal-off-diagonal elements in the Hamiltonian
        integer(kind=8) :: nDi_dp ! Number of diagonal elements for double precision
        integer(kind=8) :: nnz ! Number of non-zero elements in the Hamiltonian    
        integer(kind=8) :: nOff ! Number of off-diagonal elements in the Hamiltonian
        integer(kind=8) :: nDi ! Number of diagonal elements in the Hamiltonian
        integer(kind=8), allocatable :: occ(:,:) ! Occupation numbers for each site and state
        integer(kind=8), allocatable :: hamOff(:,:) ! Off-diagonal parts of Hamiltonian matrix(non-symmetrized)
        integer(kind=8), allocatable :: hamDi(:,:)  ! Diagonal parts of Hamiltonian matrix(non-symmetrized)
        integer(kind=8), allocatable :: rc(:,:), rcOff(:,:), rcDi(:) ! Row-column indices for the Hamiltonian matrix, its off-diagonal and diagonal parts, respectively(COO format)
        double precision, allocatable :: ham(:), ham_dp(:,:), hamOff_dp(:), hamDi_dp(:), hamDense_dp(:,:) ! Entry values of double precision Hamiltonian matrix
        ! ham: non-symmetrized 
        ! ham_dp: translation invariance with real momenta k1 and k2 = {0,Pi} 
        ! hamOff_dp: off-diagonal entries of translation invariant Hamiltonian matrix 
        ! hamDi_dp: diagonal entries of translation invariant Hamiltonian matrix 
        ! hamDense_dp: dense Hamiltonian matrix in double precision
        double complex, allocatable :: ham_dc(:) ! Complex Hamiltonian matrix(symmetrized)
        double complex, allocatable :: hamOff_dc(:) ! Complex off-diagonal Hamiltonian matrix(symmetrized)
        double complex, allocatable :: hamDi_dc(:,:) ! Complex diagonal Hamiltonian matrix(symmetrized)
        double complex, allocatable :: hamDi_off_dc(:) ! Complex diagonal entries of Hamiltonian matrix generated by off-diagonal hopping. 
        double complex, allocatable :: hamDense_dc(:,:) ! Complex dense Hamiltonian matrix. 
        
    end type hamiltonian_params


    type system_state 
        ! Parameters characterizing the current state of the system
        character(len=1) :: mat_type ! "R" for real matrices, "C" for complex matrices
        integer :: conf ! Current disorder configuration
        integer :: nv1, nv2 ! Current V1 and V2 iteration in the loop
        integer :: k1, k2  ! Current momentum components
        double precision :: v1, v2 ! Current V1 and V2 values
        logical :: first_call ! Flag indicating if this is the first call to I/O writing. 

    end type system_state
        


end module types
