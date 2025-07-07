module parameters

    implicit none

    double precision, parameter :: eps     = 1
    double precision, parameter :: thresh  = 10.**(-14) ! Threshold for numerical precision
    double precision, parameter :: pi      = 4*atan(1.d0) ! Pi constant
    double complex, parameter   :: ii = (0, 1) ! Imaginary unit
    character, parameter        :: mode*2    = 'SA' ! Mode for ARPACK diagonalization (SA = calculate nev smallest algebraic eigenvalues, LA = largest algebraic eigenvalues) 
    character, parameter        :: pattern*2 = 'AB' ! AB = First site of A sublattice, second site of B sublattice, BA = First site of B sublattice, second site of A sublattice 
    logical, parameter          :: dynamic = .False. 
    logical, parameter          :: nested  = .True.
  

end module parameters