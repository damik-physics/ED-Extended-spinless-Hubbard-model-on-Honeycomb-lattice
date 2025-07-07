module variables_old

implicit none

! Uncomment on HPC cluster
    ! integer :: othreads
    ! integer :: mthreads
    ! integer :: ucx
    ! integer :: ucy
    ! integer :: ti !Include translation symmetry

integer, parameter :: spartan = 0 
integer, parameter :: othrds  = 1
integer, parameter :: mthrds  = 1
integer, parameter :: otf     = 0  !Generate Hamiltonian on-the-fly
integer, parameter :: ucx     = 3   
integer, parameter :: ucy     = 3

integer, parameter :: ti      = 1  !Include translation symmetry
integer, parameter :: k0      = 1  !Only k=0 sector 
integer, parameter :: symmetrize = 1  !Compute basis in irrep
character(len=2), parameter :: irrep = 'E1'  !Compute basis in irrep 

integer, parameter :: tilted = 1 !Whether to generate rectangular cluster ('0') or tilted cluster ('1') 
character, parameter :: cluster*4 = '18A' !Label for tilted cluster if 'tilted == 1'

integer, parameter :: corr = 0 !Calculate correlation functions 
integer, parameter :: curr = 1 !Calculate current-current correlation functions 
integer, parameter :: refbonds = 3 !Number of refbonds for which currents are calculated: 0: all refbonds, >0: first 'refbonds'
integer, parameter :: states   = 0 !Flag for saving states

integer, parameter :: degeneracy = 1  !Flag for including degenerate ground states: 0 = non-degenerate GS, 1 = include degeneracies, 2 = quasidegeneracy
double precision, parameter :: deg = 10.**(-10) !Threshold for considering eigenvalues degenerate, i.e. Abs(E_{i}-E_{i-1}) <= deg are degenerate 

integer, parameter :: feast  = 0  !Use Feast algorithm
integer, parameter :: arpack = 1  !Use Arpack diagonalization
integer, parameter :: mkl    = 0  !Use MKL diagonalization
integer, parameter :: exact  = 0  !Force exact diagonalization

integer, parameter :: nevext = 30
integer, parameter :: n_st   = 20
integer, parameter :: ncv0   = 180

! Uncomment on HPC cluster
! integer :: nDis = 1!50
! double precision :: dis = 0.d0!1.0d0
integer, parameter :: nDis = 1!250
double precision, parameter :: dis = 0.0d0!1.0d0

integer, parameter :: nev0      = 200 !Initial FEAST subspace size
integer, parameter :: nevmax    = 230 !Sets maximum energy for FEAST from previous run   

integer, parameter :: debug     = 0  !Flag for debug mode 
integer, parameter :: dimthresh = 1 !250
integer, parameter :: g_fact    = 2

double precision, parameter :: mass = 0.d0! 10.**(-10) !Onsite potential on sublattices, M_A = mass, M_B = -mass 

double precision, parameter :: dv   = 1.0d0
double precision, parameter :: vmin = 0.0d0
double precision, parameter :: vmax = 0.0d0

double precision, parameter :: dv2   = 0.025d0
double precision, parameter :: v2min = 0.0d0
double precision, parameter :: v2max = 3.0d0

! Uncomment on HPC cluster
    ! double precision :: vmin   = 0.0d0
    ! double precision :: vmax   = 0.0d0
    ! double precision :: v2min  = 0.0d0
    ! double precision :: v2max  = 0.0d0 

! Comment out on HPC cluster
character :: gclusterext*1, tiext*1, othreadsext*3, mthreadsext*3, ucxext*1, ucyext*1, vext*16, v2minext*16, v2maxext*16, disext*16, nDisext*16
character :: prms*200
double precision, parameter :: filling = 0.5
double precision, parameter :: t       = 1.0d0
double precision, parameter :: eps     = 1
double precision, parameter :: thresh  = 10.**(-14)
double precision, parameter :: pi      = 4*atan(1.d0)

integer, parameter :: p      = 1  !Inversion symmetry eigenvalue
integer, parameter :: p1     = 1  !Reflection symmetry eigenvalue
integer, parameter :: p2     = 1  !Reflection symmetry eigenvalue
integer, parameter :: p3     = 1  !Reflection symmetry eigenvalue


character, parameter :: bc*1      = 'p' !'p'
character, parameter :: pattern*2 = 'AB' !'BA'

logical, parameter :: dynamic = .False.
logical, parameter :: nested  = .True.
logical, parameter :: rvec    = .True.

double complex, parameter :: ii = (0, 1)

external :: dsyev, dsaupd, dseupd, dmout

integer, save :: sites, particles 
integer(kind=8), save :: dim


end module variables_old
