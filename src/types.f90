module types
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
        integer :: tilted, ti, k0, symm, corr, curr, refbonds, states, otf, degflag
        integer :: feast, arpack, mkl, exact
        integer :: dimthresh
        integer :: g_fact ! Electron g-factor

        ! Diagonalization parameters
        integer :: nev, nev0, nst, ncv0, nevmax 
        logical :: rvec ! Flag for calculating eigenvectors
    end type sim_params

    type :: diag_params
        
        ! Flags for diagonalization methods
        integer :: feast ! Use FEAST algorithm for diagonalization
        integer :: arpack ! Use Sparse Arpack diagonalization
        integer :: mkl ! Use Sparse MKL diagonalization
        integer :: exact ! Use full exact diagonalization for dense matrices
        !    feast  = 1: Use FEAST algorithm for diagonalization
        !    arpack = 1: Use Sparse Arpack diagonalization
        !    mkl    = 1: Use Sparse MKL diagonalization
        !    exact  = 1: Use full exact diagonalization for dense matrices
        integer :: otf ! Flag for generating Hamiltonian on-the-fly
        integer :: degflag ! Flag for including degenerate ground states
        logical :: rvec ! Flag for calculating eigenvectors
        ! otf = 1: Flag to generate Hamiltonian on-the-fly during vector-matrix multiplication (does not store sparse Hamiltonian in arrays)
        ! degflag = 1: Flag for including degenerate ground states:
        !                0 = assume non-degenerate ground state, 1 = assume truly degenerate ground states, 2 = assume quasidegenerate ground states

        ! Diagonalization parameters
        character(len=1) :: type ! "R" for real matrices, "C" for complex matrices
        integer :: nev ! Number of eigenvalues to calculate
        integer :: nst ! Number of eigenstates to calculate
        integer :: ncv ! Number of Lanczos convergence vectors to use
        integer :: ncv0 ! Number of Lanczos vectors to use
        integer :: nev0 ! Initial FEAST subspace size
        integer :: nevmax ! Sets maximum energy for FEAST from previous run
        integer :: full ! Flag indicating type of diagonalization. 0: sparse diagonalization, 1: dense diagonalization
        integer :: dimthresh ! Threshold for minimum dimension of the Hamiltonian matrix to use sparse diagonalization methods
        ! ncv0 : Number of Lanczos vectors to use (see Arpack manual). Typical on the order of 100-200
        ! nev0 : Initial FEAST subspace size
        ! nevmax : Sets maximum energy for FEAST from previous run
        ! dimthresh : Threshold for dimension of the Hamiltonian matrix to use sparse diagonalization methods

        double precision, allocatable :: energies(:), gs(:), evals(:), norm(:), eigstate(:,:), norm2d(:,:)
        double complex, allocatable :: eigstate_dc(:,:), gs_dc(:)
    end type diag_params

    type :: thread_params
        ! Parallelization parameters
        integer :: othrds ! Number of OpenMP threads for Arpack diagonalization
        integer :: mthrds ! Number of MKL threads for MKL diagonalization
        integer :: dis_thrds ! Number of disorder threads
        integer :: v1_thrds ! Number of V1 threads
        integer :: v2_thrds ! Number of V2 threads
        integer :: num_thrds ! Number of threads for current calculations
        integer, allocatable :: units(:,:), units_2(:,:,:,:) ! Unit numbers for output files in parallel calculations
        integer :: thread_num ! Thread number for V2 calculations
        integer :: thread_num_2 ! Thread number for V1 calculations
        integer :: ndv1 ! Number of V1 values
        integer :: ndv2 ! Number of V2 values
    end type thread_params

    type out_params
        ! Output parameters    
        character(len=1000) :: dir ! Output directory
        integer :: unit ! Unit number for output files
        character(len=512) :: outdir
        integer :: corr ! Flag for calculating single-operator correlation functions
        integer :: curr ! Flag for calculating current-current correlation functions
        integer :: refbonds ! Number of reference bonds for which currents are calculated
        integer :: states ! Flag for saving eigenstates
        double precision :: deg ! Threshold for considering eigenvalues degenerate, i.e. Abs(E_{i}-E_{i-1}) <= deg are degenerate
        ! Eigenvalues and eigenvectors
        integer :: nDeg ! Number of degenerate ground states
        double precision, allocatable :: energies(:) ! Array of eigenvalues
        double precision, allocatable :: gs(:) ! Ground state energy
        double complex, allocatable :: eigstate(:,:) ! Eigenstates
        double complex, allocatable :: eigstate_dc(:,:) ! Complex eigenstates
        double complex, allocatable :: gs_dc(:) ! Complex ground state energy

        ! Parameters for correlation function calculations
        integer :: refsite ! Reference site for correlation functions


    end type out_params

    type geometry
        ! Parameters and arrays related to the basis and its construction, symmetrization and lattice geometry  

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
        integer(kind=8), allocatable :: hamOff(:,:) ! Off-diagonal parts of Hamiltonian matrix (non-symmetrized)
        integer(kind=8), allocatable :: hamDi(:,:)  ! Diagonal parts of Hamiltonian matrix (non-symmetrized)
        integer(kind=8), allocatable :: rc(:,:), rcOff(:,:), rcDi(:) ! Row-column indices for the Hamiltonian matrix, its off-diagonal and diagonal parts, respectively (COO format)
        double precision, allocatable :: ham(:), ham_dp(:,:), hamOff_dp(:), hamDi_dp(:) ! Entry values of double precision Hamiltonian matrix
        ! ham: non-symmetrized 
        ! ham_dp: translation invariance with real momenta k1 and k2 = {0,Pi} 
        ! hamOff_dp: off-diagonal entries of translation invariant Hamiltonian matrix 
        ! hamDi_dp: diagonal entries of translation invariant Hamiltonian matrix 

        double complex, allocatable :: ham_dc(:) ! Complex Hamiltonian matrix (symmetrized)
        double complex, allocatable :: hamOff_dc(:) ! Complex off-diagonal Hamiltonian matrix (symmetrized)
        double complex, allocatable :: hamDi_dc(:,:) ! Complex diagonal Hamiltonian matrix (symmetrized)
        double complex, allocatable :: hamDi_off_dc(:) ! Complex diagonal entries of Hamiltonian matrix generated by off-diagonal hopping. 
    end type hamiltonian_params

    contains 
        

    subroutine init_params(params)
        ! Initialize simulation input parameters with default values 
        type(sim_params), intent(out) :: params
        params%ucx = 3
        params%ucy = 3
        params%bc = 'p'
        params%irrep = 'A1'
        params%p1 = 1
        params%p2 = 1
        params%p3 = 1
        params%cluster = '18A'
        
        params%nDis = 1
        params%t = 1.0d0
        params%dis = 0.0d0
        params%mass = 0.0d0
        params%filling = 0.5
        params%dv1 = 0.1d0
        params%v1min = 0.0d0
        params%v1max = 1.0d0
        params%dv2 = 0.1d0
        params%v2min = 0.0d0
        params%v2max = 1.0d0

        params%othrds = 1
        params%mthrds = 1

        params%tilted = 1
        params%ti = 1
        params%k0 = 1
        params%symm = 0
        params%corr = 0
        params%curr = 0
        params%refbonds = 0
        params%states = 1
        params%feast = 0
        params%arpack = 1
        params%mkl = 0
        params%exact = 0
        params%dimthresh = 1000
        params%otf = 0
        params%degflag = 1
        params%g_fact = 2

        params%rvec = .true.
        params%nev = 10
        params%nev0 = 10
        params%nevmax = 10
        params%nst = 10
        params%ncv0 = 120

    end subroutine init_params
    

end module types
