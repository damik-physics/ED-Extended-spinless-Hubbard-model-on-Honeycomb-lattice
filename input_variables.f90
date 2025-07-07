module input_variables

implicit none

!Module declaring input variables defined in the file 'input.nml' used for the diagonalization of the extended Hubbard model on the honeycomb lattice.

    ! Lattice parameters
    integer   :: ucx = 3, ucy = 3, tilted = 1     
    character :: cluster*4 = '18A', bc*1 = 'p' 
        ! ucx, ucy: Number of unit cells in x (horizontal) and y (vertical) directions
        ! tilted: Flag to use tilted clusters (1 = yes, 0 = no)
        ! cluster: Cluster name, e.g. '18A' for a 18-site cluster
        ! bc: Boundary conDitions, 'p' for periodic, 'o' for open

    ! Symmetization parameters
    integer          :: symmetrize = 0, p1 = 1, p2 = 1, p3 = 1, ti = 0, k0 = 0  
    character(len=2) :: irrep = 'A1' 
        ! ti = 1: Include translation symmetry, k0 = 1: Only k=0 sector
        ! symmetrize = 1: Symmetrization flag; compute basis in irreducible representation
        ! irrep = 'B2': Compute basis in irreducible representation B2 (currently supported irreps = {A1, A2, B1, B2})
        ! p1, p2, p3: Reflection symmetry eigenvalues for sigma-1, sigma-2, and sigma-3 axes
        ! irrep: Irreducible representation for symmetrization 


    ! Output parameters
    integer          :: corr = 0, curr = 0, refbonds = 3, states = 0
    double precision :: deg = 10.**(-10)
        ! corr = 1: Calculate single-operator correlation functions
        ! curr = 1: Calculate current-current correlation functions
        ! deg: Threshold for considering eigenvalues degenerate, i.e. Abs(E_{i}-E_{i-1}) <= deg are degenerate
        ! refbonds : Number of reference bonds for which currents are calculated:
        !              0: all reference bonds, >0: first 'refbonds'
        ! states: Flag for saving eigenstates


    ! Diagonalization parameters
    integer :: feast = 0, arpack = 1, mkl = 0, exact = 0, nevext = 10, n_st = 10, ncv0 = 180, otf = 0, degeneracy = 1, nev0 = 200, nevmax = 230, dimthresh = 250 
    logical :: rvec = .true. 
        ! Diagonalization methods:
        !    feast  = 1: Use FEAST algorithm for diagonalization
        !    arpack = 1: Use Sparse Arpack diagonalization
        !    mkl    = 1: Use Sparse MKL diagonalization
        !    exact  = 1: Use full exact diagonalization for dense matrices
        ! rvec : Flag for calculating eigenvectors
        ! nevext : Number of eigenvalues to calculate
        ! n_st : Number of eigenstates to calculate
        ! ncv0 : Number of Lanczos vectors to use (see Arpack manual). Typical on the order of 100-200
        ! otf = 1: Flag to generate Hamiltonian on-the-fly (does not store sparse Hamiltonian in arrays)
        ! degeneracy = 1: Flag for including degenerate ground states:
        !                0 = non-degenerate ground state, 1 = include degeneracies, 2 = quasidegeneracy
        ! nev0 : Initial FEAST subspace size
        ! nevmax : Sets maximum energy for FEAST from previous run
        ! dimthresh : Threshold for dimension of the Hamiltonian matrix to use sparse diagonalization methods


    ! Parallelization parameters
    integer :: othrds = 1, mthrds = 1
        ! othrds : Number of parallel OpenMP threads for Arpack diagonalization
        ! mthrds : Number of parallel MKL threads for MKL diagonalization
        ! Note: othrds and mthrds are set to 1 if the respective diagonalization method is not used
    

    ! Hamiltonian parameters
    integer          :: nDis = 1, g_fact = 2
    double precision :: dis = 0.d0, mass = 0.d0, filling = 0.5, t = 1.0d0, dv1 = 1.d0, v1min = 0.d0, v1max =0.d0, dv2 = 0.025d0, v2min = 0.d0, v2max = 3.d0
    ! nDis : Number of disorder configurations to average over
    ! g_fact : Electron g-factor, typically 2 for spin-1/2 electrons
    ! dis  : Disorder strength
    ! mass : Onsite potential on sublattices, M_A = mass, M_B = -mass
    ! filling : Electron filling fraction (0: empty lattice, 1: fully occupied lattice)
    ! t : Hopping strength
    ! Note: The Hamiltonian is defined as H = -t * sum_{<i,j>} c_i^dagger c_j + sum_i M_A * n_{i,A} - M_B * n_{i,B} + V * sum_{<i,j>} n_i n_j + V2 * sum_{<<i,j>>} n_i n_j 
    ! where c_i^dagger and c_i are the creation and annihilation operators, n_{i,A} and n_{i,B} are the number operators for sublattices A and B, V is the nearest-neighbor interaction strength, and V2 is the next-nearest-neighbor interaction strength.
    ! Interaction strengths are defined as:
    ! V = vmin + dv * i, where i = 0, 1, ..., (vmax - vmin) / dv
    ! V2 = v2min + dv2 * j, where j = 0, 1, ..., (v2max - v2min) / dv2
     

    integer :: debug = 0 ! Flag for debug mode 


end module input_variables
