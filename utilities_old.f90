module utilities_old
    ! Module for utility functions and subroutines

    use functions 
    use parameters
    implicit none

    !Spartan changes
    ! replace all dir//"..." in file names by "..dir///...". 
    ! for Hamitlonian I/O saving replace "hamiltonians/..." in file names by "../hamiltonians/...". 



    include "omp_lib.h"
    
    !INTERFACES
 

    interface check
        module procedure check_d, check_c
    end interface  

    interface gsdeg
        module procedure gsdeg_d, gsdeg_c
    end interface 

    interface cdw
        module procedure cdw_d, cdw_c
    end interface   


    external :: dsyev, dsaupd, dseupd, dmout
    external :: mkl_set_num_threads, mkl_get_max_thrds
contains



! subroutine deallocate_all(ham, ham_dc, hamOff, hamOff_dp, hamOff_dc, hamDi,  arg2)
    !     type1, intent(in) :: arg1
    !     type2, intent(out) ::  arg2

    !     if(allocated(ham))          deallocate(ham)
    !     if(allocated(ham_dc))          deallocate(ham_dc)
    !     if(allocated(hamOff))      deallocate(hamOff)
    !     if(allocated(hamOff_dp))      deallocate(hamOff_dp)
    !     if(allocated(hamOff_dc))      deallocate(hamOff_dc)
    !     if(allocated(hamDi))      deallocate(hamDi)
    !     if(allocated(rc))           deallocate(rc)
    !     if(allocated(rcoff))           deallocate(rcoff)
    !     if(allocated(eigstate))     deallocate(eigstate)
    !     if(allocated(energies))     deallocate(energies)
    !     if(allocated(eigstate_dc))     deallocate(eigstate_dc)
    !     if(allocated(evals))        deallocate(evals)
    !     if(allocated(permutations)) deallocate(permutations)
    !     if(allocated(mombasis)) deallocate(mombasis)
    !     if(allocated(ham_d))        deallocate(ham_d)
    !     if(allocated(occupation))   deallocate(occupation)
    !     if(allocated(dplcts))   deallocate(dplcts)
    !     if(allocated(prts))   deallocate(prts)
    !     if(allocated(norm))   deallocate(norm)
    !     if(allocated(bsites))   deallocate(bsites)
    !     if(allocated(hexsites))   deallocate(hexsites)
    !     if(allocated(latticevecs))   deallocate(latticevecs)
    !     if(allocated(xtransl))   deallocate(xtransl)
    !     if(allocated(ytransl))   deallocate(ytransl)
    !     if(allocated(abasis))   deallocate(abasis)
    !     if(allocated(bbasis))   deallocate(bbasis)
    !     if(allocated(alattice))   deallocate(alattice)
    !     if(allocated(blattice))   deallocate(blattice)
    !     if(allocated(asitesbonds))   deallocate(asitesbonds)
    !     if(allocated(bsitesbonds))   deallocate(bsitesbonds)
    !     if(allocated(phases))   deallocate(phases)
    !     if(allocated(xy))   deallocate(xy)

    
! end subroutine deallocate_all


subroutine check_d(dim, nnz, nev, nest, energies, eigstate, ndeg, rc, mat)
    implicit none
    integer(kind=8), intent(in) :: dim, nnz  
    integer, intent(in) :: nev, nest
    integer, intent(in), optional :: ndeg
    integer, intent(in), optional :: rc(nnz, 2)
    double precision, intent(inout) :: energies(nev)
    double precision, intent(inout) :: eigstate(dim, nest)
    double precision, intent(in), optional :: mat(nnz)

    integer :: i = 0 
    
    ! call sortstates(nev, nest, dim, energies, eigstate)
    call checnorm(nev, nest, dim, energies, eigstate)
    call checkorthogonality(nest, dim, eigstate)
    ! if(present(mat) .and. present(rc)) then 
    !     ! print*,'Degeneracy', ndeg 
    !     do i = 1, max(1, ndeg)
    !         call expectval_d(dim, nnz, energies(i), eigstate(1:dim, i), rc, mat)
    !     end do 
    ! end if 
    return 
    
end subroutine check_d

subroutine check_c(dim, feast, nev, nest, energies, eigstate)
    implicit none
    integer(kind=8), intent(in) :: dim  
    integer, intent(in) :: nev, nest
    logical, intent(in) :: feast
    double precision, intent(inout) :: energies(nev)
    double complex, intent(inout) :: eigstate(dim, nest)

    call sortstates_comp(nev, nest, dim, energies, eigstate)
    call checnorm_comp(nev, nest, dim, energies, eigstate)
    if(feast) then 
        call checkorthogonality_comp(nest, dim, eigstate)   
    end if 

    return 

end subroutine check_c

subroutine expectval_d(dim, nnz, eval, evec, rc, mat)
    implicit none 
    integer(kind=8), intent(in) :: dim, nnz  
    integer, intent(in) :: rc(nnz, 2)
    double precision, intent(in) :: eval
    double precision, intent(in) :: evec(dim)
    double precision, intent(in) :: mat(nnz)
    
    integer :: ia(dim + 1), ja(nnz)
    double precision :: exv, val(nnz), ax(dim)

    ia = 0 
    ja = 0 
    exv = 0.d0 
    val = 0.d0 

    call coocsr(dim, nnz, mat, rc(1:nnz,1), rc(1:nnz,2), val, ja, ia) 
    call pamux (1, dim, evec, ax, val, ja, ia)
    exv = dot_product(evec, ax)
    ! print*,'Expectation value ',exv 
    ! print*,'exv-eval',exv-eval 

    return 
    
    ! call mkl_dcsrmultcsr('N', 0, 3, n, n, n, val, cols, rows, val_b, cols_b, rows_b, c, jc, ic, nnz, info)


end subroutine expectval_d

subroutine sortstates(nev, nest, dim, evals, states)
    implicit none 
    integer, intent(in) :: nev, nest 
    integer(kind=8), intent(in) :: dim 

    double precision, intent(inout) :: evals(nev), states(dim, nev) 

    integer, save :: i = 0, j = 0, cntr = 0 
    double precision, save :: temp = 0 
    double precision, allocatable, save :: tempvec(:)

    cntr = 0
    if(allocated(tempvec)) deallocate(tempvec)
    allocate(tempvec(dim))
    do i = 1, nest
        do j = 1, nest
            if (dble(evals(j)) < dble(evals(i)) .and. i < j) then 
                print*,'Swap eigenvalues i, j',i, j
                temp = evals(i)
                evals(i) = evals(j)
                evals(j) = temp 
                cntr = cntr + 1 
                if(dble(evals(i)) > dble(evals(j))) print*,'Error in resorting routine.'
                tempvec = states(i, 1:nest)
                states(i, 1:nest) = states(j, 1:nest)
                states(j, 1:nest) = tempvec 
            end if 
        end do 
    end do 
    deallocate(tempvec)
    print*,'Eigenstates sorted.', cntr, ' eigenvalues were in wrong order.'
    return 

end subroutine sortstates

subroutine sortstates_comp(nev, nest, dim, evals, states)
    implicit none 
    integer, intent(in) :: nev, nest 
    integer(kind=8), intent(in) :: dim 
    double precision, intent(inout) :: evals(nev)
    double complex, intent(inout) :: states(dim, nev) 

    integer, save :: i = 0, j = 0 
    double complex, save :: temp = 0 
    double complex, allocatable, save :: tempvec(:)

    if(allocated(tempvec)) deallocate(tempvec)
    allocate(tempvec(dim))
    do i = 1, nest
        do j = 1, nest
            if (dble(evals(j)) < dble(evals(i)) .and. i < j) then 
                print*,'Swap eigenvalues i, j',i, j
                temp     = evals(i)
                evals(i) = evals(j)
                evals(j) = temp  
                if(dble(evals(i)) > dble(evals(j))) print*,'Error in resorting routine.'
                tempvec = states(i, 1:nest)
                states(i, 1:nest) = states(j, 1:nest)
                states(j, 1:nest) = tempvec 
            end if 
        end do 
    end do 
    deallocate(tempvec)

    return 

end subroutine sortstates_comp

subroutine checnorm(nev, nest, dim, evals, states)
    implicit none 
    integer, intent(in) :: nev, nest 
    integer(kind=8), intent(in) :: dim 

    double precision, intent(in) :: evals(nev), states(dim, nest) 

    integer, save :: i = 0 
    double precision, save :: norm = 0 
    
    do i = 1, nev
        write (* ,"(x, i0, '.Eigenvalue = ',f18.13)") i, dble(evals(i))
        if(i .le. nest) then 
            norm = dot_product(states(1:dim,i),states(1:dim,i))
            if(norm < 0.9999999d0) then
                write (* ,"(x, i0, '. Norm = ',f16.12)") i, dble(norm)
                print*,'Error in subroutine checnorm.'
                error stop  
            end if 
        end if 
    end do     
    print*, 'All eigenstates have norm 1.'
    return 

end subroutine checnorm

subroutine checnorm_comp(nev, nest, dim, evals, states)
    implicit none 
    integer, intent(in) :: nev, nest 
    integer(kind=8), intent(in) :: dim 
    double precision, intent(in) :: evals(nev)
    double complex, intent(in) :: states(dim, nest) 

    integer, save :: i = 0 
    double precision, save :: norm = 0 
    double complex :: psi(dim)
    psi = 0.d0  

    do i = 1, nev
        write (* ,"(x, i0, '.Eigenvalue = ',f18.13)") i, dble(evals(i))
        if(i .le. nest) then 
            ! norm = dot_product(dconjg(states(1:dim,i)), states(1:dim,i))
            norm = dot_product(states(1:dim,i), states(1:dim,i))           
            
            if(norm < 0.9999999d0) then
                write (* ,"(x, i0, '. Norm = ',f16.12)") i, dble(norm)
                print*,'Error in subroutine checnorm.'
                error stop
            end if 
        end if 
    end do     
    print*,'All eigenstates are normalized.'
    return 
    
end subroutine checnorm_comp

subroutine checnorm_comp_feast(nev, nest, dim, evals, states)
    implicit none 
    integer, intent(in) :: nev, nest 
    integer(kind=8), intent(in) :: dim 
    double precision, intent(in) :: evals(nev)
    double complex, intent(in) :: states(dim, nest) 

    integer, save :: i = 0 
    double precision, save :: norm = 0 
    double complex :: psi(dim)
    psi = 0.d0  

    do i = 1, nev
        write (* ,"(i0, '.Eigenvalue = ',f18.13)") i, dble(evals(i))
        if(i .le. nest) then 
            ! norm = dot_product(dconjg(states(1:dim,i)), states(1:dim,i))
            norm = dot_product(states(1:dim,i), states(1:dim,i))           
            
            if(norm < 0.9999999d0) then
                write (* ,"(i0, '. Norm = ',f16.12)") i, dble(norm)
                print*,'Error in subroutine checnorm.'
                error stop
            end if 
        end if 
    end do     
    print*,'All eigenstates are normalized.'
    return 

end subroutine checnorm_comp_feast

subroutine checkorthogonality(nest, dim, states)
    implicit none 
    integer, intent(in) :: nest 
    integer(kind=8), intent(in) :: dim 

    double precision, intent(in) :: states(dim, nest) 

    integer, save :: i = 0, j = 0 
    
    do i = 1, nest 
        do j = 1, nest 
            if(dot_product(states(1:dim, i), states(1:dim, j)) < 0.9999 .and. dot_product(states(1:dim, i), states(1:dim, j)) > 0.000000001) then
                print*, 'Psi(i).Psi(j)', i, j, dot_product(states(1:dim, i), states(1:dim, j))
                error stop 
            end if 
        end do 
    end do 
    print*,'All eigenstates are orthorgonal.'
    return 

end subroutine checkorthogonality

subroutine checkorthogonality_comp(nest, dim, states)
    implicit none 
    integer, intent(in) :: nest 
    integer(kind=8), intent(in) :: dim 

    double complex, intent(in) :: states(dim, nest) 

    integer, save :: i = 0, j = 0, k = 0  
    double complex :: res = 0.d0 
    double complex, allocatable, save :: vec1(:), vec2(:)

    if(allocated(vec1)) deallocate(vec1)
    if(allocated(vec2)) deallocate(vec2)
    allocate(vec1(dim))
    allocate(vec2(dim))
    vec1 = 0.d0
    vec2 = 0.d0 
    do i = 1, nest 
        do j = 1, nest 
            res = 0.d0 
            do k = 1, dim 
                res = res + states(k, i) * dconjg(states(k, j)) 
            end do
            ! print*, 'Psi(i).Psi(j)', i, j, dot_product(states(1:dim, i), states(1:dim, j)) 
            if(abs(dot_product(states(1:dim, i), states(1:dim, j))) < 0.9999 .and. abs(dot_product(states(1:dim, i), states(1:dim, j))) > 0.000000001) then
            ! if(abs(dot_product(dconjg(states(1:dim, i)), states(1:dim, j))) < 0.9999 .and. abs(dot_product(dconjg(states(1:dim, i)), states(1:dim, j))) > 0.000000001) then
                print*, 'Psi(i).Psi(j)', i, j, dot_product(states(1:dim, i), states(1:dim, j))
                ! print*, 'Psi(i)*.Psi(j)', i, j, dot_product(dconjg(states(1:dim, i)), states(1:dim, j))
                print*,'Error in checkorthogonality_comp.'
                error stop 
            end if 
        end do 
    end do 
    print*,'Eigenstates are orthorgonal.'

    return 

end subroutine checkorthogonality_comp

subroutine combinations(m_max, n_max, allcombs, binom) !Calculates all combinations of size m_max drawn from numbers 1-n_max 

    implicit none
    integer, intent(in) :: m_max
    integer, intent(in) :: n_max
    integer, allocatable, intent(out) :: allcombs(:,:) 
    integer(kind=8), intent(out) :: binom
    integer, dimension (m_max) :: comb
    
    ! comb = 0 
    ! character (*), parameter :: fmt = '(i0' // repeat (', 1x, i0', m_max - 1) // ')'

    binom = int(fact(n_max) / (fact(m_max) * fact(max(1,n_max-m_max))),8) !Number of combinations (binomial n_max choose m_max)
    
    if(allocated(allcombs)) deallocate(allcombs)
    allocate(allcombs(binom, m_max))
    call gen (1)
  
  contains
  
    recursive subroutine gen (m)
  
      implicit none
      integer, intent (in) :: m
      integer :: n
      integer :: cntr = 0 
        
      if (m > m_max) then
        write (*, *) comb
        cntr = cntr + 1 
        allcombs(cntr, 1:m_max) = comb 
        if(cntr == binom) cntr = 0 
        ! write (*, fmt) comb
      else
        do n = 1, n_max
          if ((m == 1) .or. (n > comb (m - 1))) then
            comb (m) = n
            call gen (m + 1)
          end if
        end do
      end if
  
    end subroutine gen
  
end subroutine combinations

subroutine permute(value_min, value_max)

    implicit none
    integer, intent(in) :: value_min
    integer, intent(in) :: value_max
    integer :: position_min
    integer :: position_max
    ! integer, dimension (position_min : position_max) :: permutation
    integer, allocatable :: permutation(:)
    
    position_min = value_min
    position_max = value_max
    allocate(permutation(position_min : position_max))
    call generate (position_min)
  
  contains
  
    recursive subroutine generate (position)
  
      implicit none
      integer, intent (in) :: position
      integer :: value
  
      if (position > position_max) then
        write (*, *) permutation
      else
        do value = value_min, value_max
          if (.not. any (permutation (: position - 1) == value)) then
            permutation (position) = value
            call generate (position + 1)
          end if
        end do
      end if
  
    end subroutine generate
  
end subroutine permute










!------------------------------------------------!
!            Lookup of basis states              !
!------------------------------------------------!

subroutine findstate(dim, s, basis, loc)

    implicit none

    integer(kind=8), intent(in) :: dim, basis(dim)
    integer(kind=8), intent(in) :: s
    integer(kind=8), intent(out) :: loc

    integer(kind=8) :: left, right, mean

    left = 1
    right = dim
    do while (left <= right)
        mean = floor((left + right)/2.d0)
        if (basis(mean) < s) then
            left = mean + 1
        else if (basis(mean) > s) then
            right = mean - 1
        else
            loc = mean
            return
        end if
    end do
    loc = -1 !If no state is found
    return


end subroutine findstate

!---------------------------------------------!
!            Momentum state basis             !
!---------------------------------------------!



!-------------------------------------------------------------------!
!            Translation representatives and periodicity            !
!-------------------------------------------------------------------!
    


!Checkstate for rectangular lattices (without point group symmetries)
! subroutine checkstate_rect(s, sites, Lx, Ly, kx, ky, xtransl, ytransl, nnnVec, r, norm)

    !     implicit none
    
    !     ! Given momentum
    !     integer(kind=8), intent(in) :: s
    !     integer, intent(in) :: kx, ky
    !     integer, intent(in) :: sites 
    !     integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
    !     ! Define lattice parameters
    !     integer, intent(in) :: Lx   ! Number of unit cells in the x-direction
    !     integer, intent(in) :: Ly   ! Number of unit cells in the y-direction
    !     double precision, intent(in) :: nnnVec(2,3)
    !     integer, intent(out) :: r
    !     double precision, intent(out) :: norm
        
    !     integer :: nx = 0, ny = 0, r_p = 0
    !     integer :: ntot = 0, edgeA = 0, edgeB = 0, signt = 0, rowSt = 0   
    !     integer :: orbit = 1, t = 0, tx = 0, ty = 0, i = 0, flag = 0, rt = 0, itemp = 0 
    !     integer, allocatable :: orbits(:)
    !     double precision :: k1(2), k2(2), k(2)
    !     double precision :: a1(2), a2(2), shift(2)
    !     double precision, parameter :: tolerance = 1.0e-8
    !     double complex :: phase = 0.d0, phase2 = 0.d0  
    !     double complex :: phase_p = 0.d0
        

    !     ntot = popcnt(s)
    !     ! a1 = (/sqrt(3.d0), 0.d0/)
    !     ! a2 = 0.5*(/-1*sqrt(3.d0), 3.d0/)
    !     a1 = nnnVec(1:2, 1)
    !     a2 = nnnVec(1:2, 2)
    !     k1 = (/(2.0d0 * pi)/sqrt(3.d0), (2.d0 * pi)/3.d0/)
    !     k2 = (/0.d0, (4.d0 * pi)/3.d0/)
    !     k  = (dble(kx)/dble(Lx)) * k1 + (dble(ky)/dble(Ly)) * k2 
        
    !     signt = 1 
    !     flag  = 0 
    !     itemp  = 0 
    !     orbit = 1 
    !     edgeA = 0 
    !     edgeB = 0 
    !     t  = 0 
    !     tx = 0 
    !     ty = 0 
    !     norm    = 0.d0 
    !     phase   = 0.d0 
    !     phase_p = 0.d0 
    !     phase2  = 0.d0 
    !     r   = - 1 
    !     r_p = 0 
    !     rt  = 0 
    !     t   = 0 
    !     if(allocated(orbits)) deallocate(orbits)
    !     allocate(orbits(Lx*Ly))
        
    !     orbits    = 0 
    !     orbits(1) = s

    !     do nx = 1, Lx !x translations
    !         if(nx == 1) then !Start with initial state s
    !             tx = s 
    !         else if(nx > 1) then 
    !             tx = t 
    !         end if  

    !         do i = 1, Ly !Calculate anticommutation sign from x-translations 
    !             edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
    !             edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
    !             if(btest(tx, edgeB) .and. btest(tx, edgeA)) then !If both edges occupied, no sign
    !                 cycle 
    !             else if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
    !                 rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
    !                 signt = signt * (-1)**popcnt( ibits( tx, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'
    !             end if 
    !         end do 

    !         do i = 1, sites
    !             call mvbits(tx, xtransl(1,i)-1, 1, t, xtransl(2,i)-1) !Translate one site in x direction
    !         end do

    !         do ny = 1, Ly 
    !             ty = t
    !             signt = signt * ((-1)**(ntot-popcnt(ibits(ty, 2*Lx*(Ly-1), 2*Lx))))**popcnt(ibits(ty, 2*Lx*(Ly-1), 2*Lx)) 
    !             do i = 1, sites
    !                 call mvbits(ty, ytransl(1,i)-1, 1, t, ytransl(2,i)-1)
    !             end do
    !             orbit = orbit + 1 
    !             flag  = 1 
    !             do i = 1, orbit 
    !                 if(orbits(i) == t) then 
    !                     flag = 0 
    !                     exit 
    !                 end if 
    !             end do 
    !             if(flag == 1) orbits(orbit) = t 
    !             if(t == s) then 
    !                 shift = nx * a1 + ny * a2  
    !                 phase = phase + exp( - ii * dot_product(k, shift)) * signt 
    !             end if 
    !             if (t < s) then !Representative is already in the list                
    !                 return
    !             end if
    !         end do
    !     end do

    !     if (t == s) then !New potential representative found
    !         phase = phase + phase_p 
    !         if (abs(dble(phase)) < tolerance .and. abs(aimag(phase)) < tolerance) then 
    !             r = -1
    !             return  
    !         end if 
    !         do i = 1, size(orbits)
    !             if( orbits(i) > 0 ) rt = rt + 1 
    !         end do 
    !         r = rt 
    !         norm = r * abs(phase)**2
    !         ! norm = r * dconjg(phase)*phase 
    !         return
    !     end if
    
! end subroutine checkstate_rect


subroutine checkstate_rect(s, sites, Lx, Ly, kx, ky, irrep, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, r, norm)

    implicit none
  
    ! Given momentum
    integer(kind=8), intent(in) :: s
    integer, intent(in) :: sites, kx, ky, irrep
    integer, intent(in) :: xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    ! Define lattice parameters
    integer, intent(in) :: Lx   ! Number of unit cells in the x-direction
    integer, intent(in) :: Ly   ! Number of unit cells in the y-direction
    double precision, intent(in) :: nnnVec(2,3), id, par(6), rot(5)
    integer, intent(out) :: r
    double precision, intent(out) :: norm
    
    integer(kind=8) :: sr = 0, s0 = 0
    integer :: ntot = 0, orbsize = 0, info = 0
    integer :: orbit = 1, sign = 1, i = 0, n = 0, flag = 0, rt = 0
    integer, allocatable :: orbits(:)
    double precision :: k1(2), k2(2), k(2)
    double precision :: a1(2), a2(2)
    double precision, parameter :: tolerance = 1.0e-8
    double complex :: phase = 0.d0
    
    ntot = popcnt(s)
    ! a1 = (/sqrt(3.d0), 0.d0/)
    ! a2 = 0.5*(/-1*sqrt(3.d0), 3.d0/)
    a1 = nnnVec(1:2, 1)
    a2 = nnnVec(1:2, 2)
    k1 = (/(2.0d0 * pi)/sqrt(3.d0), (2.d0 * pi)/3.d0/)
    k2 = (/0.d0, (4.d0 * pi)/3.d0/)
    k  = (dble(kx)/dble(Lx)) * k1 + (dble(ky)/dble(Ly)) * k2 
     
    flag    = 0 
    orbit   = 0
    norm    = 0.d0 
    phase   = 0.d0 
    r       = - 1 
    rt      = 0 
    s0      = s
    orbsize = size(par) * size(rot) * Lx * Ly 
    if(allocated(orbits)) deallocate(orbits)
    allocate(orbits(orbsize))    
    orbits = 0 
    sign   = 1 
    !Identity x translations
    call translation(s0, s, sites, ntot, orbit, orbits, 0, Ly, Lx, Ly, id, sign, a2, a1, xtransl, ytransl, k, phase, info)
    if(info < 0) return 
    
    if(irrep == 0) goto 11
    !Reflections 
    do i = 1, 6
        sign = 1   
        call reflect(s0, s, sites, refl(i,1:sites), sign, info, sr) 
        if(info < 0) return        
        call translation(s0, sr, sites, ntot, orbit, orbits, 0, Ly, Lx, Ly, par(i), sign, a2, a1, xtransl, ytransl, k, phase, info)
        if(info < 0) return 
    end do 

    !Rotations  
    do n = 1, 5  
        sign = 1
        call c6n(s0, s, sites, n, c6, sign, info, sr)
        if(info < 0) return 
        call translation(s0, sr, sites, ntot, orbit, orbits, 0, Ly, Lx, Ly, rot(n), sign, a2, a1, xtransl, ytransl, k, phase, info)
        if(info < 0) return 
    end do 

    11 continue 
       
    if((abs(dble(phase)) < tolerance) .and. (abs(aimag(phase)) < tolerance)) return  

    do i = 1, size(orbits)
        if(orbits(i) > 0) rt = rt + 1 
    end do 
    r = rt 
    norm = r * abs(phase)**2

    return

  
end subroutine checkstate_rect

!Include 2D irreps: First project on IRREP subspace using 'checkstate' and then project onto IRREP basis states using Wigner projections
    ! subroutine checkstate2D(s, orbsize, sites, tilted, nHel, Lx, Ly, kx, ky, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, orbits, norm, phases, r)

    !         implicit none
        
    !         ! Given momentum
    !         integer(kind=8), intent(in) :: s
    !         integer, intent(in) :: orbsize, kx, ky, sites, tilted, nHel, xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    !         ! Define lattice parameters
    !         integer, intent(in) :: Lx   ! Number of unit cells in the a2-direction
    !         integer, intent(in) :: Ly   ! Number of unit cells in the a1-direction
    !         double precision, intent(in) :: nnnVec(2,3), id, par(6), rot(5)
            
    !         integer, intent(out) :: r  
    !         integer(kind=8), intent(out) :: orbits(orbsize,2)
    !         double precision, intent(out) :: norm(2)
    !         double complex, intent(out) :: phases(orbsize, 2)

    !         integer(kind=8) :: sr = 0, s0 = 0
    !         integer :: rep, sign = 1, ntot = 0, info = 0, layers
    !         integer :: orbit = 1, i = 0, j = 0, n = 0, rt = 0, c = 0
    !         ! integer, allocatable :: orbits(:)
    !         double precision :: k1(2), k2(2), k(2)
    !         double precision :: a1(2), a2(2)
    !         double precision, parameter :: tolerance = 1.0e-8
            
    !         double precision :: sigma(2,2), rho(2,2), reflection(2,2), rotation(2,2), identity(2,2) 
    !         double complex :: phase = 0.d0 

    !         if(rot(1) == 1) rep = 1 !IRREP E1
    !         if(rot(1) == -1) rep = 2 !IRREP E2
    !         rho(1,1)   =  cos(2*pi*rep/6) !Rotation matrix entry (1,1)
    !         rho(1,2)   = -sin(2*pi*rep/6) !Rotation matrix entry (1,2)
    !         rho(2,1)   =  sin(2*pi*rep/6) !Rotation matrix entry (2,1)
    !         rho(2,2)   =  cos(2*pi*rep/6) !Rotation matrix entry (2,2)
    !         sigma(1,1) =  1               !Reflection matrix entry (1,1)
    !         sigma(1,2) =  0               !Reflection matrix entry (1,2)
    !         sigma(2,1) =  0               !Reflection matrix entry (2,1)
    !         sigma(2,2) = -1               !Reflection matrix entry (2,2)
    !         identity = 0.d0 
    !         identity(1, 1) = 1.d0
    !         identity(2, 2) = 1.d0
    !         if(tilted == 1) then 
    !             a1     = nnnVec(1:2, 1)
    !             a2     = nnnVec(1:2, 2)
    !             layers = nHel 
    !             k1 = (/(-2.d0 * pi) / 3.d0, (-2.d0 * pi) / sqrt(3.d0) /)
    !             k2 = (/(-4.d0 * pi) / 3.d0, 0.d0 /)
    !             if(nHel == 1) then 
    !                 k  = (dble(kx)/dble(Lx)) * k2
    !             else if(nHel > 1) then 
    !                 k  = (dble(kx)/dble(Lx)) * k2 + (dble(ky)/dble(Ly)) * k1 
    !             end if   
    !         else 
    !             layers = Ly
    !             a1     = nnnVec(1:2, 2)
    !             a2     = nnnVec(1:2, 1)
    !             k1     = (/(2.0d0 * pi)/sqrt(3.d0), (2.d0 * pi)/3.d0/)
    !             k2     = (/0.d0, (4.d0 * pi)/3.d0/)
    !             k      = (dble(kx)/dble(Lx)) * k1 + (dble(ky)/dble(Ly)) * k2 
    !         end if 
    !         ntot = popcnt(s) !Number of particles      
    !         r       = -1 
    !         s0      = s
    !         norm    = 0.d0 
            
            
    !         ! if(allocated(orbits)) deallocate(orbits)
    !         ! allocate(orbits(orbsize,2))
    !         orbits = 0 
    !         do c = 1, 2 !IRREP basis states
    !             if(s0 == 48300 .or. s0 == 39326) then 
    !                 print* ,s, 's'
    !                 print* ,c, 'c'
    !             end if 
    !             orbit  = 0            
    !             sign   = 1 
    !             rt     = 0 
    !             phase  = 0.d0 
    !             phases = 0.d0 
    !             if(s0 == 48300 .or. s0 == 39326) print*, 'TRANSLATIONS ---- ', s0 
    !             call translation2D(s0, s, sites, ntot, orbsize, orbit, orbits(1:orbsize, c), tilted, layers, Lx, Ly, 1.0d0, sign, a1, a2, xtransl, ytransl, k, phases(1:orbsize, c), phase, info)
    !             if(s0 == 48300 .or. s0 == 39326) print*, 'TRANSLATIONS ---- ', s0
    !             ! if(info < 0) return 
    !             !Reflections 
    !             reflection = 0.d0 
    !             do i = 1, 6
    !                 sign = 1
    !                 ! reflection = rho 
    !                 reflection = identity 

    !                 ! do j = 1, max(modulo(2*(i-1), 6) - 1, 0)
    !                 !     reflection = matmul(reflection, rho)
    !                 ! end do 
    !                 do j = 1, i - 1
    !                     reflection = matmul(reflection, rho)
    !                 end do        
    !                 reflection = matmul(reflection, sigma)
    !                 if(s0 == 48300 .or. s0 == 39326) then 
    !                     if(reflection(c,c) == 0.d0 ) then 
    !                         print* ,'checkstate'
    !                         print* ,reflection(c,c), 'reflection(c,c)'
    !                         print* ,j, i, c, 'j, i, c'
    !                         print* ,rho(c,c), 'rho(c,c)'
    !                         print* ,sigma(c,c), 'sigma(c,c)'
    !                         print* ,identity(c,c), 'identity(c,c)'
    !                         pause 
    !                     end if 
    !                 end if 
    !                 ! if(i == 1) reflection = sigma 

    !                 call reflect(s0, s, sites, refl(i,1:sites), sign, info, sr) 
    !                 ! if(info < 0) return        
    !                 if(s0 == 48300 .or. s0 == 39326) print*, 'REFLECTIONS ---- ', s0
    !                 call translation2D(s0, sr, sites, ntot, orbsize, orbit, orbits(1:orbsize, c), tilted, layers, Lx, Ly, reflection(c,c), sign, a1, a2, xtransl, ytransl, k, phases(1:orbsize, c), phase, info)
    !                 if(s0 == 48300 .or. s0 == 39326) print*, 'REFLECTIONS -----', s0
    !                 ! if(info < 0) return 

    !                 ! if(s0 == 48300 .or. s0 == 39326) then 
    !                 !     print* ,i,  'i'
    !                 !     print* ,reflection(c,c), 'reflection(c,c)'
    !                 ! end if 
    !             end do 

    !             rotation = rho 
    !             !Rotations
    !             do n = 1, 5  
    !                 sign = 1
    !                 call c6n(s0, s, sites, n, c6, sign, info, sr)
    !                 ! if(info < 0) return 
    !                 if(s0 == 48300 .or. s0 == 39326) print*, 'ROTATIONS -----', s0 
    !                 call translation2D(s0, sr, sites, ntot, orbsize, orbit, orbits(1:orbsize,c), tilted, layers, Lx, Ly, rotation(c,c), sign, a1, a2, xtransl, ytransl, k, phases(1:orbsize, c), phase, info)
    !                 ! if(info < 0) return 
    !                 if(s0 == 48300 .or. s0 == 39326) then 
    !                     print* ,n, 'n'
    !                     print* ,rotation(c,c), 'rotation(c,c)'
    !                 end if 
    !                 if(s0 == 48300 .or. s0 == 39326) print*, 'ROTATIONS -----', s0
    !                 rotation = matmul(rotation, rho)
                    
    !             end do 
            
    !             ! if((abs(dble(phase)) < tolerance) .and. (abs(aimag(phase)) < tolerance)) return  
                
                
    !             do i = 1, orbsize
    !                 if(orbits(i,c) > 0) rt = rt + 1 
    !             end do 
    !             r       = rt 
    !             norm(c) = rt * abs(phase)**2 * 0.5 !0.5 = 1/(dim of irrep) 
    !             ! norm(c) = rt !* abs(phase)**2
                
    !             if(s0 == 48300 .or. s0 == 39326) then 
    !                 print*, 'checkstate'
    !                 print* ,rt, 'rt'
    !                 print* ,phase, 'phase'
    !                 print* ,abs(phase)**2, 'abs(phase)**2'
    !                 print* ,norm(c), 'norm(c)'
    !                 pause 
    !             end if 
    !         end do 


    !         return


    
! end subroutine checkstate2D



! subroutine reps2D()

    !     implicit none

    !     integer, intent(in) :: 
    !     integer(kind=8), intent(in) ::
    !     integer, allocatable, intent(in) :: 
    !     integer(kind=8), allocatable, intent(in) ::
    !     double precision, intent(in) :: 
    !     double precision, allocatable, intent(in) :: 
    !     character(len=*), intent(in) :: 

    !     integer, intent(out) :: 
    !     integer(kind=8), intent(out) ::
    !     integer, allocatable, intent(out) :: 
    !     integer(kind=8), allocatable, intent(out) ::
    !     double precision, intent(out) :: 
    !     double precision, allocatable, intent(out) :: 

    !     if(irrep == "E1") then 
    !         sigmav1 = [[1.0,0.0],[0.0,-1.0]]
    !         sigmav2 = [[0.5,sqrt(3.0)/2.0],[sqrt(3.0)/2.0,-0.5]]
    !         sigmav3 = [[-0.5,sqrt(3.0)/2.0],[sqrt(3.0)/2.0,0.5]]
    !         sigmad1 = [[0.0,1.0],[1.0,0.0]]
    !         sigmad2 = [[-0.5,-sqrt(3.0)/2.0],[-sqrt(3.0)/2.0,0.5]]
    !         sigmad3 = [[0.5,-sqrt(3.0)/2.0],[-sqrt(3.0)/2.0,-0.5]]
    !     else if(irrep == "E2") then 
    !         sigmav1 = [[1.0,0.0],[0.0,-1.0]]
    !         sigmav2 = [[0.5,sqrt(3.0)/2.0],[sqrt(3.0)/2.0,-0.5]]
    !         sigmav3 = [[-0.5,sqrt(3.0)/2.0],[sqrt(3.0)/2.0,0.5]]
    !         sigmad1 = [[0.0,1.0],[1.0,0.0]]
    !         sigmad2 = [[-0.5,-sqrt(3.0)/2.0],[-sqrt(3.0)/2.0,0.5]]
    !         sigmad3 = [[0.5,-sqrt(3.0)/2.0],[-sqrt(3.0)/2.0,-0.5]]
    !     end if 


! end subroutine reps2D

!Compares translated state to reflected or rotated state. Does check for phase = 0.
subroutine translation2(s0, s, sites, ntot, orbit, orbits, nHel, Lx, Ly, char, signpg, a1, a2, xtransl, ytransl, k, phase, info)
    
    implicit none 
    integer, intent(in) :: sites, ntot, nHel, Lx, Ly, char, signpg
    integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
    integer(kind=8), intent(in) :: s0, s
    double precision, intent(in) :: a1(2), a2(2), k(2)
    integer, intent(out) :: info  
    integer, intent(inout) :: orbit
    integer, allocatable, intent(inout) :: orbits(:)
    double complex, intent(inout) :: phase

    integer :: nx, ny, i, rowst, edgeA, edgeB, flag, orbit0
    integer(kind=8) :: t, tx, ty, sign
    double precision :: shift(2)
    double complex :: phaset
    double precision, parameter :: tolerance = 1.0e-8

    
    i = 0
    t = 0 
    tx = 0 
    ty = 0
    nx = 0
    ny = 0
    info = 0 !On output: InDicates whether new state was found (info = 0) or state is already contained in list (info < 1)
    sign = 1
    edgeA = 0
    edgeB = 0
    rowSt = 0
    !Include initial state (if not already in the list) to orbit 
    orbit = orbit + 1 
    orbit0 = orbit 
    flag  = 1 
    ! phaset = phase + signpg * char 
    phaset = 0.d0 
    ! phaset = 0.d0 


    do i = 1, orbit 
        if(orbits(i) == s) then !State already in orbit
            flag = 0 
            exit 
        end if 
    end do 
    if(flag == 1) orbits(orbit) = s !New state found, add to orbit

    do nx = 1, Lx !x translations
        
        if(nx == 1) then !Start with initial state s
            tx = s 
        else if(nx > 1) then 
            tx = t 
        end if  
        
        if(nhel == 1) then !Single helix, quasi-1D
            edgeB = sites-1 !Last B site in helix
            edgeA = sites-2 !Last A site in helix
            if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                rowSt = 0 !First site 
                sign = sign * (-1)**popcnt( ibits( tx, rowst, sites - 2 ) ) !Occupation of rest of helix 
            end if 
        else !More than one helix
            do i = 1, nHel!Ly !Calculate anticommutation sign from x-translations 
                edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                if(btest(tx, edgeB) .and. btest(tx, edgeA)) then !If both edges occupied, no sign
                    cycle 
                else if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                    rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
                    sign  = sign * (-1)**popcnt( ibits( tx, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'
                end if 
                ! edgeB = (i-1)*2*Lx + 1 !First B site on y-layer 'i'
                ! edgeA = (i-1)*2*Lx     !First A site on y-layer 'i'
                ! if(btest(tx, edgeB) .and. btest(tx, edgeA)) then !If both edges occupied, no sign
                !     cycle 
                ! else if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                !     rowSt = edgeB + 1!First site on y-layer 'i' 
                !     sign  = sign * (-1)**popcnt( ibits( tx, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'
                ! end if 
            end do 
        end if 
        
        do i = 1, sites
            call mvbits(tx, xtransl(1,i)-1, 1, t, xtransl(2,i)-1) !Translate one site in x direction
        end do
 

        if(nHel == 1) then !No y-translations for single helix
            orbit = orbit + 1 
            flag  = 1 
            do i = 1, orbit 
                if(orbits(i) == t) then 
                    flag = 0 
                    exit 
                end if 
            end do 
            if(flag == 1) orbits(orbit) = t 
            shift = nx * a2  
            phaset = phaset + exp( - ii * dot_product(k, shift)) * sign * char * signpg 

            if(t == s) then 
                ! shift = nx * a2  
                ! phase = phase + exp( - ii * dot_product(k, shift)) * sign * char * signpg 
                phase = phase + phaset 
                phaset = 0.d0 
                
            end if 
            if (t < s0) then !Representative is already in the list 
                info = -1 
                return
            end if

        ! if(nHel > 1) then 
        else if(nHel > 1) then 
            do ny = 1, Ly 
                ty = t
                ! sign = sign * ((-1)**(ntot-popcnt( ibits( ty, 2*Lx*(nHel-1), 2*Lx ) ) ))**popcnt( ibits( ty, 2*Lx*(nHel-1), 2*Lx ) ) 
           
                ! sign = sign * ((-1)**modulo(popcnt( ibits( ty, 2*Lx*(nHel-1), 2*Lx )), 2))**modulo(ntot-popcnt( ibits( ty, 2*Lx*(nHel-1), 2*Lx ) ) , 2)  

                sign = sign * ((-1)**modulo(popcnt( ibits( ty, 0, 2*Lx )), 2))**modulo(ntot-popcnt( ibits( ty, 0, 2*Lx ) ) , 2) 
                if( ((-1)**modulo(popcnt( ibits( ty, 0, 2*Lx )), 2))**modulo(ntot-popcnt( ibits( ty, 0, 2*Lx ) ) , 2) == -1) stop 'minus sign' 

                ! sign = sign * ((-1)**(ntot-popcnt( ibits( ty, 0, 2*Lx ) ) ))**popcnt( ibits( ty, 0, 2*Lx ) ) 
                do i = 1, nHel
                    ! edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                    ! edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                    edgeB = (i-1)*2*Lx + 1 !First B site on y-layer 'i'
                    edgeA = (i-1)*2*Lx     !First A site on y-layer 'i'
                    if(btest(ty, edgeB) .and. btest(ty, edgeA)) then !If both edges occupied, no sign
                        cycle 
                    else if(btest(ty, edgeB) .neqv. btest(ty, edgeA)) then !If only one edge occupied 
                        ! rowSt = (i-1)*2*Lx!First site on y-layer 'i' 
                        rowSt = edgeB + 1!First site on y-layer 'i' 
                        sign = sign * (-1)**popcnt( ibits( ty, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'

                        ! if(s0 == 39326 .or. s0 == 6559 ) then 
                    
                        !     print* ,s0, ty, 's0, ty'
                        !     write(*,'(X,B0)') ty
                        !     print* ,edgeA, 'edgeA'
                        !     print* ,edgeB, 'edgeB'
                        !     print* ,rowSt, 'rowSt'
                        !     print* ,2*Lx-2, '2*Lx-2'
                        !     print* ,popcnt( ibits( ty, rowst, 2*Lx-2 ) ), 'popcnt( ibits( ty, rowst, 2*Lx-2 ) )'
                        !     write(*,'(X,B0)') ibits( ty, rowst, 2*Lx-2 )
                        !     print* ,(-1)**popcnt( ibits( ty, rowst, 2*Lx-2 ) ), '(-1)**popcnt( ibits( ty, rowst, 2*Lx-2 ) )'
                        !     print* ,sign, 'sign'
                            
                        ! end if 
                    end if                 
                    
                end do     


                do i = 1, sites
                    call mvbits(ty, ytransl(1,i)-1, 1, t, ytransl(2,i)-1)
                end do


                orbit = orbit + 1 
                flag  = 1 
                do i = 1, orbit 
                    if(orbits(i) == t) then 
                        flag = 0 
                        exit 
                    end if 
                end do 
                if(flag == 1) orbits(orbit) = t 

                shift = ny * a1 + nx * a2  
                ! phaset = phaset + exp( - ii * dot_product(k, shift)) * sign * char * signpg 
                phaset = phaset + exp( - ii * dot_product(k, shift)) * sign * char * signpg 
                if(t == s) then 
                    ! shift = ny * a1 + nx * a2  
                    ! phase = phase + exp( - ii * dot_product(k, shift)) * sign * char * signpg 
                    ! phase  = phase + phaset 
                    ! phaset = 0.d0
                    ! sign   = 1
                    if((abs(dble(phase)) < tolerance) .and. (abs(aimag(phase)) < tolerance)) then
                        print* ,s0, 's0'
                        print*, 'zero phase'
                        pause 
                    end if 
                end if 
                if (t < s0) then !Representative is already in the list 
                    info = -1
                    return
                end if
                
            end do
            
        end if 

    end do

    
    if((abs(dble(phase)) < tolerance) .and. (abs(aimag(phase)) < tolerance)) then 
        orbits(orbit0:orbit) = 0 
        print*, 'FLUSHED ORBIT ---------------------------------------------'
    end if 
    if(t .ne. s ) print*, 'Translation orbit incomplete.'
    if(t .ne. s ) stop

end subroutine translation2



subroutine translation2D(s0, s, sites, ntot, orbsize, orbit, orbits, tilted, nHel, Lx, Ly, char, signpg, a1, a2, xtransl, ytransl, k, phases, phase, info)

    !If tilted == 0, swap input vectors a1 and a2 and set nHel = Ly 
    implicit none 
    integer, intent(in) :: orbsize, sites, ntot, tilted, nHel, Lx, Ly, signpg
    integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
    integer(kind=8), intent(in) :: s0, s
    double precision, intent(in) :: a1(2), a2(2), k(2), char
    integer, intent(out) :: info  
    integer, intent(inout) :: orbit
    integer(kind=8), intent(inout) :: orbits(orbsize)
    double complex, intent(inout) :: phases(orbsize), phase 

    integer :: nx, ny, i, rowst, edgeA, edgeB, flag
    integer(kind=8) :: t, tx, ty, sign 
    double precision :: shift(2)

        t = 0 
        tx = 0 
        ty = 0
        nx = 0
        ny = 0
        info = 0 !On output: InDicates whether new state was found (info = 0) or state is already contained in list (info < 1)
        sign = 1
        edgeA = 0
        edgeB = 0
        rowSt = 0
        flag  = 1 
        ! !Include initial state (if not already in the list) to orbit 
        ! orbit = orbit + 1 
        

        ! do i = 1, orbit 
        !     if(orbits(i) == s) then !State already in orbit
        !         flag = 0 
        !         exit 
        !     end if 
        ! end do 
        ! if(flag == 1) orbits(orbit) = s !New state found, add to orbit

        do nx = 1, Lx !x translations
            if(nx == 1) then !Start with initial state s
                tx = s 
            else if(nx > 1) then 
                tx = t 
            end if  
            
            if(nhel == 1) then !Single helix, quasi-1D
                edgeB = sites-1 !Last B site in helix
                edgeA = sites-2 !Last A site in helix
                if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                    rowSt = 0 !First site 
                    sign = sign * (-1)**popcnt( ibits( tx, rowst, sites - 2 ) ) !Occupation of rest of helix 
                end if 
            else !More than one helix
                do i = 1, nHel!Ly !Calculate anticommutation sign from x-translations 
                    edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                    edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                    if(btest(tx, edgeB) .and. btest(tx, edgeA)) then !If both edges occupied, no sign
                        cycle 
                    else if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                        rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
                        sign = sign * (-1)**popcnt(ibits(tx, rowst, 2*Lx-2)) !Occupation of y-layer 'i'
                    end if 
                end do 
            end if 

            do i = 1, sites
                call mvbits(tx, xtransl(1,i)-1, 1, t, xtransl(2,i)-1) !Translate one site in x direction
            end do

            if(nHel == 1) then !No y-translations for single helix
                orbit = orbit + 1 
                flag  = 1 
                do i = 1, orbit 
                    if(orbits(i) == t) then 
                        flag = 0 
                        orbit = orbit - 1
                        exit 
                    end if 
                end do 
                if(flag == 1) orbits(orbit) = t 
                
                shift = nx * a2  
                phases(orbit) = phases(orbit) + exp( - ii * dot_product(k, shift)) * sign * char * signpg 
                if(t == s0) phase = phase + phases(orbit)
                if(t < s0) then !Representative is already in the list 
                    info = -1 
                    return
                end if

            else if(nHel > 1) then 

                do ny = 1, Ly 
                    ty    = t
                    shift = ny * a1 + nx * a2  
                    sign  = sign * ((-1)**(ntot-popcnt(ibits(ty, 2*Lx*(nHel-1), 2*Lx))))**popcnt(ibits(ty, 2*Lx*(nHel-1), 2*Lx)) 
                    if(tilted == 1) then !Minus sign from re-ordering after y-translation. Only if cluster is tilted.  
                        do i = 1, nHel
                            edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                            edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                            if(btest(ty, edgeB) .and. btest(ty, edgeA)) then !If both edges occupied, no sign
                                cycle 
                            else if(btest(ty, edgeB) .neqv. btest(ty, edgeA)) then !If only one edge occupied 
                                rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
                                sign = sign * (-1)**popcnt( ibits( ty, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'
                            end if 
                        end do     
                    end if 
                    do i = 1, sites 
                        call mvbits(ty, ytransl(1,i)-1, 1, t, ytransl(2,i)-1)
                    end do
                    orbit = orbit + 1 
                    flag  = 1 
                    do i = 1, orbit 
                        if(orbits(i) == t) then 
                            flag = 0 
                            orbit = orbit - 1 
                            phases(i) = phases(i) + exp(-ii * dot_product(k, shift)) * sign * char * signpg 

                            if(t == s0) phase = phase + exp(-ii * dot_product(k, shift)) * sign * char * signpg 
                            if(s0 == 511 .and. t == s0) then 
                                print*, 'trans2D.1'

                                print* ,i, 'orbit position'
                                print* ,sign, 'sign'
                                print* ,signpg, 'signpg'    
                                print* ,char, 'character'
                                print* , exp(-ii * dot_product(k, shift)) * sign * char * signpg , ' phi+'
                                print* ,phase, 'phase'
                            end if
                            exit 
                        end if 
                    end do 
                    if(flag == 1) then 
                        orbits(orbit) = t 
                        phases(orbit) = phases(orbit) + exp(-ii * dot_product(k, shift)) * sign * char * signpg 
                        if(t == s0) phase = phase + exp(-ii * dot_product(k, shift)) * sign * char * signpg
                        if(s0 == 511 .and. t == s0) then 
                            print*, 'trans2D.2'
                                
                            print* ,orbit, 'orbit position'
                            print* ,sign, 'sign'
                            print* ,signpg, 'signpg'    
                            print* ,char, 'character'
                            print* , exp(-ii * dot_product(k, shift)) * sign * char * signpg , ' phi+'
                            print* ,phase, 'phase'
                        end if
                    end if 
            
      

                    if (t < s0) then !Representative is already in the list 
                        info = -1
                        return
                    end if
                end do
            end if 
                
        end do

!         if(s0 == 39326 .and. t == 87382) then 
! pause
!         end if

        ! if(s0 == 39326) then 
            
        !     print* ,phases(1:orbsize), 'phases(1:orbsize)'
        !     print* ,s0, 's0, trans2D'
        !     ! pause 
        ! end if

        if(t .ne. s ) print*, 'Translation orbit incomplete.'
        if(t .ne. s ) stop

end subroutine translation2D

subroutine xtranslate(s, nHel, n, Lx, xtransl, sign, sx)
    implicit none 
    integer(kind=8), intent(in) :: s 
    integer, intent(in) :: nHel, n, Lx
    integer, intent(in) :: xtransl(2, n)
    integer, intent(inout) :: sign 
    integer(kind=8), intent(out) :: sx 

    integer :: i = 0, edgeA = 0, edgeB = 0, rowSt = 0
    integer(kind = 8) :: t = 0


    t = 0
    if(nhel == 1) then !Single helix, quasi-1D
        edgeB = n-1 !Last B site in helix
        edgeA = n-2 !Last A site in helix
        if(btest(s, edgeB) .neqv. btest(s, edgeA)) then !If only one edge occupied 
            rowSt = 0 !First site 
            sign = sign * (-1)**popcnt( ibits( s, rowst, n-2 ) ) !Occupation of rest of helix 
        end if 
    else
        do i = 1, nHel !Calculate anticommutation sign from x-translations 
            edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
            edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
            if(btest(s, edgeB) .and. btest(s, edgeA)) then !If both edges occupied, no sign
                cycle 
            else if(btest(s, edgeB) .neqv. btest(s, edgeA)) then !If only one edge occupied 
                rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
                sign = sign * (-1)**popcnt(ibits(s, rowst, 2*Lx-2)) !Occupation of y-layer 'i'
            end if 
        end do 
    end if 

    do i = 1, n
        call mvbits(s, xtransl(1,i)-1, 1, t, xtransl(2,i)-1)      
    end do
    sx = t 
    

end subroutine xtranslate

subroutine ytranslate(s, nHel, tilt, n, Lx, ytransl, sign, sy)
    !Set tilt < 0 if cluster is not tilted, i.e. rectangular. 
    implicit none 
    integer(kind=8), intent(in) :: s 
    integer, intent(in) :: nHel, tilt, n, Lx
    integer, intent(in) :: ytransl(2, n)
    integer, intent(inout) :: sign 
    integer(kind=8), intent(out) :: sy

    integer :: i = 0, ntot = 0, edgeA = 0, edgeB = 0, rowSt = 0
    integer(kind = 8) :: t = 0
    
    t = 0
    ntot = popcnt(s)
    if(nhel == 1) then !Single helix
        sign = sign * ((-1)**(ntot-popcnt(ibits(s, n - tilt, tilt))))**popcnt(ibits(s, n - tilt, tilt)) 
    else 
        sign = sign * ((-1)**(ntot-popcnt(ibits(s, 2*Lx*(nHel-1), 2*Lx))))**popcnt(ibits(s, 2*Lx*(nHel-1), 2*Lx)) 
        if(tilt < 0) goto 11
        do i = 1, nHel
            edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
            edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
            if(btest(s, edgeB) .and. btest(s, edgeA)) then !If both edges occupied, no sign
                cycle 
            else if(btest(s, edgeB) .neqv. btest(s, edgeA)) then !If only one edge occupied 
                rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
                sign = sign * (-1)**popcnt( ibits( s, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'
            end if 
        end do     
        11 continue 
    end if 

    do i = 1, n
        call mvbits(s, ytransl(1,i)-1, 1, t, ytransl(2,i)-1)
    end do
    sy = t

end subroutine ytranslate

!Uses modulo in calculating sign changes
subroutine ytranslate2(s, nHel, tilt, n, Lx, ytransl, sign, sy)
    implicit none 
    integer(kind=8), intent(in) :: s 
    integer, intent(in) :: nHel, tilt, n, Lx
    integer, intent(in) :: ytransl(2, n)
    integer, intent(inout) :: sign 
    integer(kind=8), intent(out) :: sy

    integer :: i = 0, ntot = 0, edgeA = 0, edgeB = 0, rowSt = 0
    
    sy = 0
    ntot = popcnt(s)
    if(nhel == 1) then !Single helix
        ! sign = sign * ((-1)**(ntot-popcnt( ibits( s, n - tilt, tilt ) ) ))**popcnt( ibits( s, n - tilt, tilt ) ) 
        sign = sign * ((-1)**(modulo(ntot-popcnt( ibits( s, n - tilt, tilt)), 2)))**modulo(popcnt(ibits(s, n - tilt, tilt)), 2)
    else 
        sign = sign * ((-1)**modulo(ntot-popcnt(ibits(s, 0, 2*Lx)), 2))**modulo(popcnt(ibits(s, 0, 2*Lx)), 2)
        ! sign = sign * ((-1)**(ntot-popcnt(ibits(s, 0, 2*Lx))))**popcnt(ibits(s, 0, 2*Lx))
        do i = 1, nHel
            edgeA = (i-1)*2*Lx !First A site on y-layer 'i'
            edgeB = (i-1)*2*Lx + 1 !First B site on y-layer 'i'
            if(btest(s, edgeB) .and. btest(s, edgeA)) then !If both edges occupied, no sign
                cycle 
            else if(btest(s, edgeB) .neqv. btest(s, edgeA)) then !If only one edge occupied 
                rowSt = (i-1)*2*Lx + 2 !First site after edge A on y-layer 'i' 
                sign = sign * (-1)**popcnt( ibits( s, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'
            end if 
        end do     
    end if 

    do i = 1, n
        call mvbits(s, ytransl(1,i)-1, 1, sy, ytransl(2,i)-1)
    end do

end subroutine ytranslate2


!---------------------------------------------!
!            Find representative              !
!---------------------------------------------!

subroutine representative_irrep(s, n, nHel, tilt, Lx, Ly, symmetrize, id, par, rot, xtransl, ytransl, refl, c6, r, l1, l2, sign)
    !Finds the representative 'r' for state 's'. 'n' is the number of sites and 'l' the number of shifts needed to translate 's' to 'r'.
    implicit none
    integer(kind=8), intent(in) :: s
    integer, intent(in) :: n, nHel, tilt, symmetrize, Lx, Ly, xtransl(2, n), ytransl(2, n), refl(6, n), c6(n)
    double precision, intent(in) :: id, par(6), rot(5)
    integer(kind=8), intent(out) :: r
    integer, intent(out) :: l1, l2, sign 

    integer :: nx = 0, ny = 0, i = 0 
    integer :: edgeA = 0, edgeB = 0, rowSt = 0, signt = 1
    integer(kind=8) :: t = 0, tx = 0, ty = 0

    l2    = Lx 
    l1    = Ly
    sign  = 1
    r     = s
    t     = 0
    tx    = 0 
    ty    = 0 
    signt = 1
    rowSt = 0 
    edgeA = 0 
    edgeB = 0 

    do nx = 1, Lx !Identity x translations 
        if(nx == 1) tx = s 
        if(nx > 1) tx = t
        call xtranslate(tx, nHel, n, Lx, xtransl, signt, t)
        if (t < r) then !New rep found
            sign = signt !Update sign 
            r    = t
            l1   = 0
            l2   = nx
        end if
        if(nHel > 1) then 
            do ny = 1, Ly 
                ty = t                
                !Note, 'signt' and 't' are carried on through out all translations. 
                call ytranslate(ty, nHel, tilt, n, Lx, ytransl, signt, t)
                if (t < r) then
                    sign = signt 
                    r    = t
                    l1   = ny
                    l2   = nx
                end if
            end do
        end if 
    end do !Leaving this loop, r <= s 

    if(symmetrize == 0) go to 12
    do i = 1, 6 !Check for representatives among reflected states
        call mirror_rep(r, s, n, nHel, tilt, par(i), Lx, Ly, refl(i, 1:n), xtransl, ytransl, sign, l1, l2)
    end do 
    do i = 1, 5 !Check for representatives among rotated states
        call rotate_rep(r, s, n, nHel, tilt, i, rot(i), Lx, Ly, c6, xtransl, ytransl, sign, l1, l2)
    end do 
    12 continue

    return 

end subroutine representative_irrep

subroutine representative_irrep_rect(s, n, Lx, Ly, symmetrize, id, par, rot, xtransl, ytransl, refl, c6, r, l1, l2, sign)
    !Finds the representative 'r' for state 's'. 'n' is the number of sites and 'l' the number of shifts needed to translate 's' to 'r'.
    implicit none
    integer(kind=8), intent(in) :: s
    integer, intent(in) :: n, symmetrize, Lx, Ly, xtransl(2, n), ytransl(2, n), refl(6, n), c6(n)
    double precision, intent(in) :: id, par(6), rot(5)
    integer(kind=8), intent(out) :: r
    integer, intent(out) :: l1, l2, sign 

    integer :: nx = 0, ny = 0, i = 0 
    integer :: edgeA = 0, edgeB = 0, rowSt = 0, signt = 1
    integer(kind=8) :: t = 0, tx = 0, ty = 0

    l1    = Lx 
    l2    = Ly
    sign  = 1
    r     = s
    t     = 0
    tx    = 0 
    ty    = 0 
    signt = 1
    rowSt = 0 
    edgeA = 0 
    edgeB = 0 

    do nx = 1, Lx !Identity x translations 
        if(nx == 1) tx = s 
        if(nx > 1) tx = t
        call xtranslate(tx, Ly, n, Lx, xtransl, signt, t)
        if (t < r) then !New rep found
            sign = signt !Update sign 
            r    = t
            l1   = 0
            l2   = nx
        end if

        do ny = 1, Ly 
            ty = t                
            !Note, 'signt' and 't' are carried on through out all translations. 
            call ytranslate(ty, Ly, -1, n, Lx, ytransl, signt, t)
            if (t < r) then
                sign = signt 
                r    = t
                l1   = nx
                l2   = ny
            end if
        end do

    end do !Leaving this loop, r <= s 

    if(symmetrize == 0) go to 12
    do i = 1, 6 !Check for representatives among reflected states
        call mirror_rep(r, s, n, Ly, -1, par(i), Lx, Ly, refl(i, 1:n), xtransl, ytransl, sign, l2, l1)
    end do 
    do i = 1, 5 !Check for representatives among rotated states
        call rotate_rep(r, s, n, Ly, -1, i, rot(i), Lx, Ly, c6, xtransl, ytransl, sign, l2, l1)
    end do 
    12 continue

    return 

end subroutine representative_irrep_rect

!Uses ytranslate2 for checkstate5 
subroutine representative_irrep2(s, n, nHel, tilt, Lx, Ly, symmetrize, id, par, rot, xtransl, ytransl, refl, c6, r, l1, l2, sign)
    !Finds the representative 'r' for state 's'. 'n' is the number of sites and 'l' the number of shifts needed to translate 's' to 'r'.
    implicit none
    integer(kind=8), intent(in) :: s
    integer, intent(in) :: n, nHel, tilt, symmetrize, Lx, Ly, xtransl(2, n), ytransl(2, n), refl(6, n), c6(n)
    double precision, intent(in) :: id, par(6), rot(5)
    integer(kind=8), intent(out) :: r
    integer, intent(out) :: l1, l2, sign 

    integer :: nx = 0, ny = 0, i = 0 
    integer :: edgeA = 0, edgeB = 0, rowSt = 0, signt = 1
    integer(kind=8) :: t = 0, tx = 0, ty = 0

    l2    = Lx 
    l1    = Ly
    sign  = 1
    r     = s
    t     = 0
    tx    = 0 
    ty    = 0 
    signt = 1
    rowSt = 0 
    edgeA = 0 
    edgeB = 0 

    do nx = 1, Lx !Identity x translations 
        if(nx == 1) tx = s 
        if(nx > 1) tx = t
        call xtranslate(tx, nHel, n, Lx, xtransl, signt, t)
        if (t < r) then !New rep found
            sign = signt !Update sign 
            r    = t
            l1   = 0
            l2   = nx
        end if
        if(nHel == 1) then 
            if (t < r) then !Why twice?
                sign = signt 
                r    = t
                l1   = 0
                l2   = nx
            end if
        else if(nHel > 1) then 
            do ny = 1, Ly 
                ty = t                
                !Note, 'signt' and 't' are carried on through out all translations. 
                call ytranslate2(ty, nHel, tilt, n, Lx, ytransl, signt, t)
                if (t < r) then
                    sign = signt 
                    r    = t
                    l1   = ny
                    l2   = nx
                end if
            end do
        end if 
    end do !Leaving this loop, r <= s 

    if(symmetrize == 0) go to 12
    do i = 1, 6 !Check for representatives among reflected states
        call mirror_rep2(r, s, n, nHel, tilt, par(i), Lx, Ly, refl(i, 1:n), xtransl, ytransl, sign, l1, l2)
    end do 
    do i = 1, 5 !Check for representatives among rotated states
        call rotate_rep2(r, s, n, nHel, tilt, i, rot(i), Lx, Ly, c6, xtransl, ytransl, sign, l1, l2)
    end do 
    12 continue

    return 

end subroutine representative_irrep2

subroutine representative(s, n, Lx, Ly, xtransl, ytransl, r, l1, l2, sign)
    !Finds the representative 'r' state for state 's'. 'n' is the number of sites and 'l' the number of shifts needed to translate 's' to 'r'.
    implicit none
    integer, intent(in) :: n
    integer(kind=8), intent(in) :: s
    integer, intent(in) :: Lx, Ly
    integer, intent(in) :: xtransl(2, n), ytransl(2, n)
    integer(kind=8), intent(out) :: r
    integer, intent(out) :: l1, l2, sign 

    integer :: i, nx, ny, ntot, signt, itemp 
    integer :: edgeA, edgeB, rowSt
    integer(kind=8) :: t, tx, ty 

    ntot  = popcnt(s)
    itemp = 0 
    l1    = Lx 
    l2    = Ly
    sign  = 1
    r     = s
    t     = 0
    tx    = 0 
    ty    = 0 
    signt = 1
    rowSt = 0 
    edgeA = 0 
    edgeB = 0 

    do nx = 1, Lx
        if(nx == 1) then 
            tx = s 
        else if(nx > 1) then 
            tx = t 
        end if  

        do i = 1, Ly
            edgeB = i*2*Lx-1
            edgeA = i*2*Lx-2 
            if(btest(tx, edgeB) .and. btest(tx, edgeA)) then 
                cycle 
            else if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then 
                rowSt = (i-1)*2*Lx
                signt = signt * (-1)**popcnt( ibits( tx, rowst, 2*Lx-2 ) )  
            end if 
        end do      

        do i = 1, n
            call mvbits(tx, xtransl(1,i)-1, 1, t, xtransl(2,i)-1)      
        end do

        do ny = 1, Ly
            ty    = t
            signt = signt * ((-1)**(ntot-popcnt( ibits( ty, 2*Lx*(Ly-1), 2*Lx ) ) ))**popcnt( ibits( ty, 2*Lx*(Ly-1), 2*Lx ) ) 
            do i = 1, n
                call mvbits(ty, ytransl(1,i)-1, 1, t, ytransl(2,i)-1)       
            end do
            if (t < r) then
                sign = signt 
                r    = t
                l1   = nx
                l2   = ny
            end if
        end do 
    end do 

    return 

end subroutine representative

!For tilted clusters
subroutine representative2(s, n, nHel, tilt, Lx, Ly, xtransl, ytransl, r, l1, l2, sign)
    !Finds the representative 'r' state for state 's'. 'n' is the number of sites and 'l' the number of shifts needed to translate 's' to 'r'.
    implicit none
    integer, intent(in) :: n, nHel, tilt 
    integer(kind=8), intent(in) :: s
    integer, intent(in) :: Lx, Ly
    integer, intent(in) :: xtransl(2, n), ytransl(2, n)
    integer(kind=8), intent(out) :: r
    integer, intent(out) :: l1, l2, sign 

    integer :: i, nx, ny, ntot, signt, itemp 
    integer :: edgeA, edgeB, rowSt
    integer(kind=8) :: t, tx, ty 

    ntot  = popcnt(s)
    itemp = 0 
    l2    = Lx 
    l1    = Ly
    sign  = 1
    r     = s
    t     = 0
    tx    = 0 
    ty    = 0 
    signt = 1
    rowSt = 0 
    edgeA = 0 
    edgeB = 0 

    do nx = 1, Lx
        if(nx == 1) then 
            tx = s 
        else if(nx > 1) then 
            tx = t 
        end if  

        if(nhel == 1) then !Single helix, quasi-1D
            edgeB = n-1 !Last B site in helix
            edgeA = n-2 !Last A site in helix
            if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                rowSt = 0 !First site 
                signt = signt * (-1)**popcnt( ibits( tx, rowst, n-2 ) ) !Occupation of rest of helix 
            end if 
        else
            do i = 1, nHel !Calculate anticommutation sign from x-translations 
                edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                if(btest(tx, edgeB) .and. btest(tx, edgeA)) then !If both edges occupied, no sign
                    cycle 
                else if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                    rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
                    signt = signt * (-1)**popcnt( ibits( tx, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'
                end if 
            end do 
        end if 
        do i = 1, n
            call mvbits(tx, xtransl(1,i)-1, 1, t, xtransl(2,i)-1)      
        end do
        if(nHel == 1) then 
            if (t < r) then
                sign = signt 
                r    = t
                l1   = ny
                l2   = nx
            end if
        else if(nHel > 1) then 
            do ny = 1, Ly 
                ty = t
                if(nhel == 1) then !Single helix
                    signt = signt * ((-1)**(ntot-popcnt( ibits( ty, n - tilt, tilt ) ) ))**popcnt( ibits( ty, n - tilt, tilt ) ) 
                else 
                    signt = signt * ((-1)**(ntot-popcnt( ibits( ty, 2*Lx*(nHel-1), 2*Lx ) ) ))**popcnt( ibits( ty, 2*Lx*(nHel-1), 2*Lx ) ) 
                    do i = 1, nHel
                        edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                        edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                        if(btest(ty, edgeB) .and. btest(ty, edgeA)) then !If both edges occupied, no sign
                            cycle 
                        else if(btest(ty, edgeB) .neqv. btest(ty, edgeA)) then !If only one edge occupied 
                            rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
                            signt = signt * (-1)**popcnt( ibits( ty, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'
                        end if 
                    end do     
                end if 

                do i = 1, n
                    call mvbits(ty, ytransl(1,i)-1, 1, t, ytransl(2,i)-1)
                end do
                if (t < r) then
                    sign = signt 
                    r    = t
                    l1   = ny
                    l2   = nx
                end if
            end do
        end if 
    end do

    return 

end subroutine representative2

subroutine mirror_rep(r, s, n, nHel, tilt, p, Lx, Ly, refl, xtransl, ytransl, sign, l1, l2)
    !Swap l1 and l2, set nHel = Ly, and tilt < 0, if cluster is not tilted, i.e. rectangular. 
    implicit none
    integer(kind=8), intent(in) :: s !original state
    integer, intent(in) :: n, nHel, tilt, Lx, Ly, refl(n), xtransl(2,n), ytransl(2,n)     
    double precision, intent(in) :: p
    integer(kind=8), intent(inout) :: r !potential rep (r<=s) 
    integer, intent(inout) :: l1, l2 
    integer, intent(inout) :: sign !If 'r' remains unchanged, previous sign remains unchanged. Otherwise updated. 

    integer :: nx = 0, ny = 0, signp = 1, signt = 1, info = 0
    integer(kind=8) :: t = 0, tx = 0, ty = 0 

    signp = 1
    signt = 1
    call reflect(s, s, n, refl, signp, info, t) ! compare to s, reflect s onto t
    if (t < r) then !compare to previously identified potenital representative 'r' 
        sign = signp * p 
        r    = t
        l1   = 0 !No shifts for k=0
        l2   = 0 !No shifts for k=0
    end if
    do nx = 1, Lx
        tx = t 
        call xtranslate(tx, nHel, n, Lx, xtransl, signt, t)
        if(nHel == 1) then 
            if (t < r) then
                sign = signt * p * signp  
                r    = t
                l1   = 0
                l2   = nx
            end if
        else if(nHel > 1) then 
            do ny = 1, Ly 
                ty = t               
                call ytranslate(ty, nHel, tilt, n, Lx, ytransl, signt, t)                
                if (t < r) then
                    sign = signt * signp * p
                    r    = t
                    l1   = ny
                    l2   = nx
                end if
            end do
        end if 
    end do

end subroutine mirror_rep

subroutine mirror_rep2(r, s, n, nHel, tilt, p, Lx, Ly, refl, xtransl, ytransl, sign, l1, l2)
    implicit none
    integer(kind=8), intent(in) :: s !original state
    integer, intent(in) :: n, nHel, tilt, Lx, Ly, refl(n), xtransl(2,n), ytransl(2,n)     
    double precision, intent(in) :: p
    integer(kind=8), intent(inout) :: r !potential rep (r<=s) 
    integer, intent(inout) :: l1, l2 
    integer, intent(inout) :: sign !If 'r' remains unchanged, previous sign remains unchanged. Otherwise updated. 

    integer :: nx = 0, ny = 0, signp = 1, signt = 1, info = 0
    integer(kind=8) :: t = 0, tx = 0, ty = 0 

    signp = 1
    signt = 1
    call reflect(s, s, n, refl, signp, info, t) ! compare to s, reflect s onto t
    if (t < r) then !compare to previously identified potenital representative 'r' 
        sign = signp * p 
        r    = t
        l1   = 0 !No shifts for k=0
        l2   = 0 !No shifts for k=0
    end if
    do nx = 1, Lx
        tx = t 
        call xtranslate(tx, nHel, n, Lx, xtransl, signt, t)
        if(nHel == 1) then 
            if (t < r) then
                sign = signt * p * signp  
                r    = t
                l1   = 0
                l2   = nx
            end if
        else if(nHel > 1) then 
            do ny = 1, Ly 
                ty = t               
                call ytranslate2(ty, nHel, tilt, n, Lx, ytransl, signt, t)                
                if (t < r) then
                    sign = signt * signp * p
                    r    = t
                    l1   = ny
                    l2   = nx
                end if
            end do
        end if 
    end do

end subroutine mirror_rep2

subroutine rotate_rep(r, s, n, nHel, tilt, nrot, eval, Lx, Ly, c6, xtransl, ytransl, sign, l1, l2)
    implicit none
    !Swap l1 and l2, set nHel = Ly, and tilt < 0, if cluster is not tilted, i.e. rectangular. 
    integer(kind=8), intent(in) :: s !original state
    integer, intent(in) :: n, nHel, tilt, nrot, Lx, Ly, c6(n), xtransl(2,n), ytransl(2,n)     
    double precision, intent(in) :: eval
    integer(kind=8), intent(inout) :: r !potential rep (r<=s) 
    integer, intent(inout) :: l1, l2, sign !If 'r' remains unchanged, previous sign remains unchanged. Otherwise updated. 

    integer :: nx = 0, ny = 0, signrot = 1, signt = 1, info = 0
    integer(kind=8) :: t = 0, tx = 0, ty = 0 

    signrot = 1
    signt   = 1 
    call c6n(s, s, n, nrot, c6, signrot, info, t) ! compare to s, reflect s onto t
    if (t < r) then !compare to previously identified potenital representative 'r' 
        sign = signrot * eval 
        r    = t
        l1   = 0 !No shifts for k=0
        l2   = 0 !No shifts for k=0
    end if
    do nx = 1, Lx
        tx = t 
        call xtranslate(tx, nHel, n, Lx, xtransl, signt, t)
        if(nHel == 1) then 
            if (t < r) then
                sign = signt * eval * signrot  
                r    = t
                l1   = 0
                l2   = nx
            end if
        else if(nHel > 1) then 
            do ny = 1, Ly 
                ty = t               
                call ytranslate(ty, nHel, tilt, n, Lx, ytransl, signt, t)                
                if (t < r) then
                    sign = signt * eval * signrot
                    r    = t
                    l1   = ny
                    l2   = nx
                end if
            end do
        end if 
    end do

end subroutine rotate_rep

subroutine rotate_rep2(r, s, n, nHel, tilt, nrot, eval, Lx, Ly, c6, xtransl, ytransl, sign, l1, l2)
    implicit none
    integer(kind=8), intent(in) :: s !original state
    integer, intent(in) :: n, nHel, tilt, nrot, Lx, Ly, c6(n), xtransl(2,n), ytransl(2,n)     
    double precision, intent(in) :: eval
    integer(kind=8), intent(inout) :: r !potential rep (r<=s) 
    integer, intent(inout) :: l1, l2 
    integer, intent(inout) :: sign !If 'r' remains unchanged, previous sign remains unchanged. Otherwise updated. 

    integer :: nx = 0, ny = 0, signrot = 1, signt = 1, info = 0
    integer(kind=8) :: t = 0, tx = 0, ty = 0 

    signrot = 1
    signt   = 1 
    call c6n(s, s, n, nrot, c6, signrot, info, t) ! compare to s, reflect s onto t
    if (t < r) then !compare to previously identified potenital representative 'r' 
        sign = signrot * eval 
        r    = t
        l1   = 0 !No shifts for k=0
        l2   = 0 !No shifts for k=0
    end if
    do nx = 1, Lx
        tx = t 
        call xtranslate(tx, nHel, n, Lx, xtransl, signt, t)
        if(nHel == 1) then 
            if (t < r) then
                sign = signt * eval * signrot  
                r    = t
                l1   = 0
                l2   = nx
            end if
        else if(nHel > 1) then 
            do ny = 1, Ly 
                ty = t               
                call ytranslate2(ty, nHel, tilt, n, Lx, ytransl, signt, t)                
                if (t < r) then
                    sign = signt * eval * signrot
                    r    = t
                    l1   = ny
                    l2   = nx
                end if
            end do
        end if 
    end do

end subroutine rotate_rep2

subroutine reflect_rep(r, s, n, nHel, tilt, p, Lx, Ly, refl, xtransl, ytransl, sign, l1, l2)
    implicit none
    integer(kind=8), intent(in) :: s
    integer, intent(in) :: n, nHel, tilt, Lx, Ly, refl(n), xtransl(2,n), ytransl(2,n)     
    double precision, intent(in) :: p
    integer(kind=8), intent(inout) :: r
    integer, intent(inout) :: l1, l2, sign

    integer :: nx = 0, ny = 0, signp = 1, signt = 1, info = 0
    integer(kind=8) :: t = 0, tx = 0, ty = 0 

    signp = 1
    signt = 1
    call reflect(r, s, n, refl, signp, info, t) ! compare to r, reflect s onto t
    if (t < r) then
        sign = signp * p 
        r    = t
        l1   = 0
        l2   = 0
    end if
    do nx = 1, Lx
        tx = t 
        call xtranslate(tx, nHel, n, Lx, xtransl, signt, t)
        if(nHel == 1) then 
            if (t < r) then
                sign = signt * p * signp  
                r    = t
                l1   = 0
                l2   = nx
            end if
        else if(nHel > 1) then 
            do ny = 1, Ly 
                ty = t               
                call ytranslate(ty, nHel, tilt, n, Lx, ytransl, signt, t)                
                if (t < r) then
                    sign = signt * signp * p
                    r    = t
                    l1   = ny
                    l2   = nx
                end if
            end do
        end if 
    end do

end subroutine reflect_rep
    








!-------------------------------------------------------------------------!
!            Calculate number of eigenvalues and Lanczos vectors          !
!-------------------------------------------------------------------------!

subroutine nevncv(unit, parameters, thresh, exact, nevext, nestext, ncv0, dim, full, nev, ncv, nest)

    implicit none
    integer, intent(in) :: unit, thresh, exact, nevext, nestext, ncv0
    integer(kind=8), intent(in) :: dim
    character(len=*), intent(in) :: parameters
    integer, intent(out) :: full, nev, ncv, nest
    character*200 :: filename

    if (exact == 0) then
        if (dim < thresh) then
            full = 1
            nev = dim
            ncv = dim
            nest = dim
        else
            full = 0
            if (dim == 1) then
                nev = 1
                ncv = 1
                nest = 1
            else
                ncv = min(ncv0, dim - 1) !nev+1<=ncv<=dim and ncv maximal
                ! ncv = max(ncv0, 2*nev + 10) !nev+1<=ncv<=dim and ncv maximal
                nev = min(nevext, ncv - 3)
                ! nev = min(nevext, dim - 10)
                ! if (ncv > dim) ncv = dim
                nest = min(nev, nestext)
            end if
        end if
    else if (exact == 1) then
        full = 1
        nev = dim
        ncv = dim
        nest = min(dim, nestext)
    end if
    print*, 'Number of eigenvalues = ', nev
    print*, 'Number of eigenvectors = ', nest
    print*, ''

    ! filename=dir//"Hilbert_space_dimension_"//parameters
    ! open(unit, file=filename)
    ! write(unit,*) dim 
    ! close(unit)

end subroutine nevncv

!------------------------------------------!
!            Set initial variables         !
!------------------------------------------!

subroutine set_vars()

    use variables
    use input_variables
    implicit none

    character*10:: clusterl
    integer :: threads

    ! Set variables for parallelization. 
    !$ call omp_set_dynamic(dynamic)
    !$ call omp_set_nested(nested)
    !$ call omp_set_num_threads(othrds)

    call mkl_set_num_threads(mthrds)
    call mkl_get_max_threads(threads)
    print*,'Number of MKL THREADS = ',threads
    
    !$ threads = omp_get_num_threads()
    print*,'Number of OMP THREADS = ',threads

    ! Set number of sites and particles depenDing on input parameters. 
    if(tilted == 0) then 
        sites = 2*ucx*ucy 
    else if(cluster(1:2) == "16") then 
        sites = 16 
    else if(cluster(1:2) == "18") then 
        sites = 18 
    else if(cluster(1:2) == "20") then 
        sites = 20 
    else if(cluster(1:2) == "24") then 
        sites = 24
    else if(cluster(1:2) == "30") then 
        sites = 30 
    else if(cluster(1:2) == "32") then 
        sites = 32
    else if(cluster(1:1) == "L") then 
        clusterl = cluster
        read(clusterl(2:3),*) sites
    end if 
    particles = int(filling * sites)
    
    
    if(symmetrize == 1)print('(1x, a,a)'),  'Irrep     = ', irrep
    if(tilted == 1) then
        print('(1x, a,a)'), 'Cluster   = ', cluster
    else 
        print('(1x, a,i0)'), 'UCX       = ', ucx
        print('(1x, a,i0)'), 'UCY       = ', ucy
    end if 
    print('(1x, a,i0)'), 'Sites     = ', sites
    print('(1x, a,i0)'), 'Particles = ', particles
    print*, ''
    

    return

end subroutine set_vars




  



!----------------------------------!
!         Test symmetry            !
!----------------------------------!

subroutine findstate2(dim, target, array, loc)

    implicit none

    integer(kind=8), intent(in) :: target(2), dim
    integer, intent(in) :: array(dim, 2)
    integer(kind=8), intent(out) :: loc
    

    integer(kind=8) :: i, left, right, mean

    left = 1
    right = dim
    do while (left <= right)
        mean = floor((left + right)/2.d0)
        if (array(mean, 1) < target(1)) then
            left = mean + 1
        else if (array(mean, 1) > target(1)) then
            right = mean - 1
        else
            loc = mean
            exit
        end if
    end do
    
    do i = 0, dim - loc - 1
        if(array(loc + i, 2) == target(2) .and. array(loc + i, 1) == target(1)) then 
            loc = loc + i 
            return
        end if 
    end do 

    loc = -1 !If no entry is found
    return


end subroutine findstate2

subroutine testsymm_sparse2(symmetric, nz, rc, ham)

    implicit none

    integer(kind=8), intent(in) :: nz 
    integer, intent(in) :: rc(nz,2)
    double precision, intent(in) :: ham(nz)
    
    logical, intent(out) :: symmetric
    
    integer(kind=8) :: i = 0, loc = 0, target(2) = 0
    double precision :: eps = 10.**(-10)
    
    symmetric = .true. 
    do i = 1, nz 
        target = rc(i,1:2) 
        call findstate2(nz, target, rc, loc)
        if(abs(ham(i) - ham(loc)) > eps) then 
            print* ,ham(i), rc(i,1), rc(i,2), 'ham(i), rc(i,1), rc(i,2)' 
            print* ,ham(loc), rc(loc,1), rc(loc,2), 'ham(loc), rc(loc,1), rc(loc,2)' 
            symmetric = .false. 
            exit 
        end if 
    end do 

    if(.not. symmetric) then
        print*, 'MATRIX IS NON-SYMMETRIC!'
        print*, ''
        error stop 
    else
        print*, 'MATRIX IS SYMMETRIC'
        print*, ''
    end if
    return



end subroutine testsymm_sparse2

subroutine testsymm_sparse(symmetric, dim, nz, rc, ham, norm, prts, dplcts)

    implicit none

    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: nz 
    integer, intent(in) :: rc(nz,2)
    integer, intent(in) :: prts(nz), dplcts(nz)
    double precision, intent(in) :: ham(nz)
    double precision, intent(in) :: norm(dim)
    
    logical, intent(out) :: symmetric
    
    integer(kind=8) :: i = 0, j = 0
    double precision :: eps = 10.**(-10)
    integer, allocatable :: parmat(:,:), dupmat(:,:)
    double precision, allocatable :: mat(:,:)
    
    if(allocated(mat)) deallocate(mat)
    if(allocated(parmat)) deallocate(parmat)
    if(allocated(dupmat)) deallocate(dupmat)
    allocate(mat(dim, dim))
    allocate(parmat(dim, dim))
    allocate(dupmat(dim, dim))
    mat = 0
    do i = 1, nz 
        mat(rc(i,1),rc(i,2)) = ham(i) 
        parmat(rc(i,1),rc(i,2)) = prts(i)
        dupmat(rc(i,1),rc(i,2)) = dplcts(i)
        if(rc(i,1) == 0 .or. rc(i,1) > dim ) print*, rc(i,1), 'rc(i,1)'
        if(rc(i,2) == 0 .or. rc(i,2) > dim ) print*, rc(i,2), 'rc(i,2)'
    end do 
 
    symmetric=.true.
    do i = 1, dim
        do j = 1, dim
            if(abs(mat(i,j) - mat(j,i)) > eps) then
                print*, i, 'i'
                print*, j, 'j'
                print*,'mat(i, j)', mat(i, j)
                print*,'mat(j, i)', mat(j, i)
                print*,'Difference', abs(mat(i,j) - mat(j,i))
                print*,'i-th norm ', norm(i) 
                print*,'j-th norm ', norm(j) 
                print*,'Parity ij', parmat(i, j)
                print*,'Parity ji', parmat(j, i)
                print*,'dplcts ij', dupmat(i, j)
                print*,'dplcts ji', dupmat(j, i)
                symmetric=.false.
                exit
            enDif
        end do
    end do
    if(.not. symmetric) then
        print*, 'MATRIX IS NON-SYMMETRIC!'
        print*, ''
        error stop 
    else
        print*, 'MATRIX IS SYMMETRIC'
        print*, ''
    end if
    return



end subroutine testsymm_sparse

subroutine testhermitian_sparse(hermit, irrep, symmetrize, dim, nz, rc, ham, norm, prts, dplcts)!, norm2D, orbits, phases

    implicit none

    integer(kind=8), intent(in) :: dim
    !integer(kind=8), allocatable, intent(in) :: orbits(:,:,:)
    integer, intent(in) :: nz, symmetrize
    integer, intent(in) :: rc(nz,2)
    integer, intent(in) :: prts(nz)
    integer, intent(in) :: dplcts(nz)
    double precision, allocatable, intent(in) :: norm(:)!, norm2D(:,:)
    double complex, allocatable, intent(in) :: ham(:)!, phases(:,:,:) 
    character, intent(in) :: irrep*2
    
    logical, intent(out) :: hermit
    
    integer(kind=8) :: i, j, mdim
    double precision :: eps = 10.**(-10)
    integer, allocatable :: parmat(:,:), dupmat(:,:)
    double complex, allocatable :: mat(:,:)

    if(irrep(1:1) .ne. "E") then 
        mdim = dim 
    else 
        mdim = 2 * dim 
    end if
    if(allocated(parmat)) deallocate(parmat)
    if(allocated(dupmat)) deallocate(dupmat)
    if(allocated(mat))    deallocate(mat)
    allocate(parmat(mdim, mdim))
    allocate(dupmat(mdim, mdim))
    allocate(mat(mdim, mdim))
    parmat = 0 
    dupmat = 0 
    mat = 0


    do i = 1, nz 
        ! print*, i, 'i'
        if(i == 25426 .or. i == 25488) then 
            print* , 'nz' ,nz
            print* ,i, 'i'
            print* ,ham(i), rc(i,1), rc(i,2), 'ham(i), row, col'
            pause 
        end if 

        mat(rc(i,1),rc(i,2)) = ham(i)

        ! parmat(rc(i,1),rc(i,2)) = prts(i)
        ! dupmat(rc(i,1),rc(i,2)) = dplcts(i)
        if(rc(i,1) == 0 .or. rc(i,1) > mdim ) print*, rc(i,1),i, 'rc(i,1), i'
        if(rc(i,2) == 0 .or. rc(i,2) > mdim ) print*, rc(i,2),i, 'rc(i,2), i'
    end do 

    
    hermit=.true.
    do i = 1, mdim
        do j = 1, mdim
            if(dble(conjg(mat(i,j))) - dble(mat(j,i)) > eps .or. aimag(conjg(mat(i,j))) - aimag(mat(j,i)) > eps) then
                if(dble(conjg(mat(i,j))) - dble(mat(j,i)) > eps) then
                    print*, 'RE(M*(i,j))', dble(conjg(mat(i,j)))
                    print*, 'RE(M(j,i))', dble(mat(j,i))
                    print*, 'Difference', abs(dble(conjg(mat(i,j))) - dble(mat(j,i)))
                    if(irrep(1:1) .ne. "E") then 
                        print*, 'i-th norm ', norm(i) 
                        print*, 'j-th norm ', norm(j)
                    end if 
                    print*, ''
                end if

                if(aimag(conjg(mat(i,j))) - aimag(mat(j,i)) > eps) then
                    print*, 'IM(M*(i,j))', aimag(conjg(mat(i,j)))
                    print*, 'IM(M(j,i))', aimag(mat(j,i))
                    print*, 'Difference', abs(aimag(conjg(mat(i,j))) - aimag(mat(j,i)))
                    if(irrep(1:1) .ne. "E") then 
                        print*, 'i-th norm ', norm(i) 
                        print*, 'j-th norm ', norm(j)
                    end if 
                    print*, ''
                end if

                print*, 'i', i
                print*, 'j', j
                print*, 'M(j,i)', mat(j,i)
                print*, 'M(i,j)', mat(i,j)
                print*,  ''

                hermit=.false.
                
            enDif
        end do
    end do
    if(.not. hermit) then
        print*, 'MATRIX IS NON-HERMITIAN!'
        print*, ''
        error stop
    else
        print*, 'MATRIX IS HERMITIAN'
        print*, ''
    end if
    return

end subroutine testhermitian_sparse

subroutine testsymm(sym, dim, mat)

    implicit none

    integer(kind=8), intent(in) :: dim
    double precision, intent(in) :: mat(dim, dim)
    double precision :: eps = 0.000001 
    logical, intent(out) :: sym

    integer :: i, j


    sym=.true.
    do i = 1, dim
        do j = 1, dim
            if(mat(i,j) - mat(j,i) > eps) then
                sym=.false.

                print*, i,j, 'i,j'
                print*, mat(i,j), 'mat(i,j)'
                print*, mat(j,i), 'mat(j,i)'

                exit
            enDif
        end do
    end do
    if(.not. sym) then
        print*, 'MATRIX IS NON-SYMMETRIC!'
        print*,''
        print*,''
        error stop
    else
        print*, 'MATRIX IS SYMMETRIC'
    end if
    return

end subroutine testsymm














!-----------------------------------------!
!            Charge density wave          !
!-----------------------------------------!

subroutine cdw_d(dir, parameters, conf, unit, dim, ucx, ucy, bc, pat, nbonds, nDis, dis, basis, bsites, psi, rho, rhotot)

    implicit none

    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: conf, unit, ucx, ucy, nbonds, nDis
    integer, intent(in) :: bsites(2,nbonds)
    integer(kind=8), intent(in) :: basis(dim) 
    double precision, intent(in) :: psi(dim), dis
    character(len=*), intent(in) :: dir, parameters, bc, pat

    double precision, allocatable, intent(out) :: rho(:,:)  
    double precision, allocatable, intent(out) :: rhotot(:,:)  

    integer, save :: nbx  = 0 !Number of x bonds 
    integer, save :: i = 0, j = 0 
    integer, save :: cntr = 0 
    integer, allocatable, save :: cntr1(:) !Counter for rho just within unit cells for both sublattices
    integer, allocatable, save :: cntr2(:) !Counter for rho just within unit cells for both sublattices
    integer, allocatable, save :: cntrtot1(:) !Counter for rho on all NN bonds for both sublattices
    integer, allocatable, save :: cntrtot2(:) !Counter for rho on all NN bonds for both sublattices
    logical :: append 

    !$omp threadprivate(nbx, i, j, cntr, cntr1, cntr2, cntrtot1, cntrtot2)
 
    if((nDis > 1) .and. (dis .ne. 0.d0) .and. (conf > 1)) then 
        append = .true.
    else 
        append = .false.
    end if 

    if(bc == 'o') then 
        nbx = (2 * ucx - 1) * ucy 
    else if(bc == 'p') then
        nbx = (2 * ucx) * ucy !total number of x bonds
    end if
    
    if( allocated(cntr1) ) deallocate(cntr1)
    allocate( cntr1(ucx*ucy) )
    if( allocated(cntr2) ) deallocate(cntr2)
    allocate( cntr2(ucx*ucy) )
    if( allocated(rho) ) deallocate(rho)
    allocate( rho(ucx*ucy, 3) )

    !Calculate charge differences for bonds lying within a unit cell (to be adjusted for general patterns AB or BA as well as PBC)
    rho = 0 
    do j = 1, dim 
        cntr  = 0 
        cntr1 = 0 
        cntr2 = 0 
        do i = 1, nbx !run over all x bonds
            if ( modulo( bsites(1,i), 2 ) == 0 ) cycle !Don't count x bonds between unit cells 
            cntr = cntr + 1
            if ( btest( basis(j), bsites(1,i) - 1 ) ) cntr1(cntr) = cntr1(cntr) + 1 
            if ( btest( basis(j), bsites(2,i) - 1 ) ) cntr2(cntr) = cntr2(cntr) + 1 
            ! rho(cntr) = rho(cntr) + cntr(cntr) * psi(j)**2 
            if (pat == 'AB' ) then !A = cntr1, B = cntr2
                rho(cntr,1) = rho(cntr,1) + (cntr1(cntr) - cntr2(cntr)) * psi(j)**2 !<n_A> - <n_B>
                rho(cntr,2) = rho(cntr,2) + cntr1(cntr) * psi(j)**2 !<n_A>
                rho(cntr,3) = rho(cntr,3) + cntr2(cntr) * psi(j)**2 !<n_B>
            else if (pat == 'BA' ) then !A = cntr2, B = cntr1
                rho(cntr,1) = rho(cntr,1) + (cntr2(cntr) - cntr1(cntr)) * psi(j)**2 !<n_A> - <n_B>
                rho(cntr,2) = rho(cntr,2) + cntr2(cntr) * psi(j)**2 !<n_A>
                rho(cntr,3) = rho(cntr,3) + cntr1(cntr) * psi(j)**2 !<n_B>
            end if 
        end do 
    end do 

    if( allocated(cntrtot1) ) deallocate(cntrtot1)
    allocate( cntrtot1(nbonds) )
    if( allocated(cntrtot2) ) deallocate(cntrtot2)
    allocate( cntrtot2(nbonds) )
    if( allocated(rhotot) ) deallocate(rhotot)
    allocate( rhotot(nbonds, 3) )

    rhotot = 0 
    do j = 1, dim  
        cntrtot1 = 0 
        cntrtot2 = 0 
        do i = 1, nbonds 
            ! Determine sign for each sublattice (to be fixed to A=+1 and B=-1)
            ! if ( modulo( bsites(1,i), 2 ) == 1 ) then 
            !     sgn = 1
            ! else if ( modulo( bsites(1,i), 2 ) == 0 ) then 
            !     sgn = - 1
            ! end if   
            ! if( btest(basis(j), bsites(1,i) - 1) .and. (.not. btest(basis(j), bsites(2,i) - 1) ) ) cntrtot1(i) = cntrtot1(i) + 1 
            ! if( btest(basis(j), bsites(2,i) - 1) .and. (.not. btest(basis(j), bsites(1,i) - 1) ) ) cntrtot2(i) = cntrtot2(i) + 1 
            if( btest(basis(j), bsites(1,i) - 1) ) cntrtot1(i) = cntrtot1(i) + 1 
            if( btest(basis(j), bsites(2,i) - 1) ) cntrtot2(i) = cntrtot2(i) + 1 

            
            ! rhotot(i,1) = rhotot(i,1) + cntrtot(i,1) * psi(j)**2 
            if (pat == 'AB' .and. modulo( bsites(1,i), 2 ) == 1 ) then !A = cntrtot1, B = cntrtot2
                rhotot(i,1) = rhotot(i,1) + (cntrtot1(i) - cntrtot2(i)) * psi(j)**2 !<n_A> - <n_B>
                rhotot(i,2) = rhotot(i,2) + cntrtot1(i) * psi(j)**2 !<n_A>
                rhotot(i,3) = rhotot(i,3) + cntrtot2(i) * psi(j)**2 !<n_B>
            else if (pat == 'AB' .and. modulo( bsites(1,i), 2 ) == 0 ) then !A = cntrtot2 B = cntrtot1
                rhotot(i,1) = rhotot(i,1) + (cntrtot2(i) - cntrtot1(i)) * psi(j)**2 !<n_A> - <n_B>
                rhotot(i,2) = rhotot(i,2) + cntrtot2(i) * psi(j)**2 !<n_A>
                rhotot(i,3) = rhotot(i,3) + cntrtot1(i) * psi(j)**2 !<n_B>
            else if (pat == 'BA' .and. modulo( bsites(1,i), 2 ) == 1 ) then !A = cntrtot2 B = cntrtot1
                rhotot(i,1) = rhotot(i,1) + (cntrtot2(i) - cntrtot1(i)) * psi(j)**2 !<n_A> - <n_B>
                rhotot(i,2) = rhotot(i,2) + cntrtot2(i) * psi(j)**2 !<n_A>
                rhotot(i,3) = rhotot(i,3) + cntrtot1(i) * psi(j)**2 !<n_B>
            else if (pat == 'BA' .and. modulo( bsites(1,i), 2 ) == 0 ) then !A = cntrtot1, B = cntrtot2
                rhotot(i,1) = rhotot(i,1) + (cntrtot1(i) - cntrtot2(i)) * psi(j)**2 !<n_A> - <n_B>
                rhotot(i,2) = rhotot(i,2) + cntrtot1(i) * psi(j)**2 !<n_A>
                rhotot(i,3) = rhotot(i,3) + cntrtot2(i) * psi(j)**2 !<n_B>
            end if 
        end do 

    end do 

    if( allocated(cntr1) ) deallocate(cntr1)
    if( allocated(cntr2) ) deallocate(cntr2)
    if( allocated(cntrtot1) ) deallocate(cntrtot1)
    if( allocated(cntrtot2) ) deallocate(cntrtot2)

    call save(dir, parameters, append, unit, ucx, ucy, nbonds, rho, rhotot)
    return 

end subroutine cdw_d 

subroutine ddcf(dir, refsite, conf, unit, dim, sites, k1, k2, nDis, dis, v1, v2, basis, psi)

    implicit none

    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: refsite, conf, unit, sites, nDis, k1, k2
    integer(kind=8), intent(in) :: basis(dim) 
    double precision, intent(in) :: psi(dim), dis, v1, v2
    character(len=*), intent(in) :: dir

    double precision :: refrho
    double precision, allocatable :: rhocn(:), rho(:), rhoij(:)  

    integer :: i = 0, j = 0 
    integer :: flag = 0 
    logical :: append 
    character*400 :: params 

    !$omp threadprivate(i, j)
    call parddcf(unit, k1, k2, refsite, v1, v2, params)
    if((nDis > 1) .and. (dis .ne. 0.d0) .and. (conf > 1)) then 
        append = .true.
    else 
        append = .false.
    end if 
    
    if(allocated(rho)) deallocate(rho)
    if(allocated(rhoij)) deallocate(rhoij)
    if(allocated(rhocn)) deallocate(rhocn)
    allocate(rho(sites))
    allocate(rhoij(sites))
    allocate(rhocn(sites))

    rho = 0.d0 !density correlation function of site i 
    refrho = 0.d0 !density correlation function of refsite 
    rhoij = 0.d0 !density-density correlation function between sites i and refsite 
    rhocn = 0.d0 !Connected density-density correlation function
    do i = 1, sites !run over all sites
        do j = 1, dim         
            if(btest(basis(j), i - 1)) then 
                rho(i) = rho(i) + psi(j)**2 
                flag = 1
                if(btest(basis(j), refsite - 1)) then     
                    if(i == 1) refrho = refrho + psi(j)**2
                    rhoij(i) = rhoij(i) + psi(j)**2 
                end if 
            end if 
        end do 
        rhocn(i) = rhoij(i) - rho(i) * refrho
    end do 
    rhocn = rhocn / sites 
    
    call saveddcf(dir, params, append, unit, sites, rhocn, rho)
    if(allocated(rho)) deallocate(rho)
    if(allocated(rhoij)) deallocate(rhoij)
    if(allocated(rhocn)) deallocate(rhocn)
    return 

end subroutine ddcf 

subroutine cdw_c(dir, parameters, conf, unit, dim, ucx, ucy, bc, pat, nbonds, nDis, dis, basis, bsites, psi, rho, rhotot)

    implicit none

    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: conf, unit, ucx, ucy, nbonds, nDis
    integer, intent(in) :: bsites(2,nbonds)
    integer(kind=8), intent(in) :: basis(dim) 
    double complex, intent(in) :: psi(dim), dis 
    character(len=*), intent(in) :: dir, parameters, bc, pat

    double precision, allocatable, intent(out) :: rho(:,:)  
    double precision, allocatable, intent(out) :: rhotot(:,:)  

    integer, save :: nbx  = 0 !Number of x bonds 
    integer, save :: i = 0, j = 0 
    integer, save :: cntr = 0 
    integer, allocatable, save :: cntr1(:) !Counter for rho just within unit cells for both sublattices
    integer, allocatable, save :: cntr2(:) !Counter for rho just within unit cells for both sublattices
    integer, allocatable, save :: cntrtot1(:) !Counter for rho on all NN bonds for both sublattices
    integer, allocatable, save :: cntrtot2(:) !Counter for rho on all NN bonds for both sublattices
    logical :: append 

    !$omp threadprivate(nbx, i, j, cntr, cntr1, cntr2, cntrtot1, cntrtot2)
   
    if((nDis > 1) .and. (dis .ne. 0.d0) .and. (conf > 1)) then 
        append = .true.
    else 
        append = .false.
    end if 
    if(bc == 'o') then 
        nbx = (2 * ucx - 1) * ucy 
    else if(bc == 'p') then
        nbx = (2 * ucx) * ucy !total number of x bonds
    end if
    
    if( allocated(cntr1) ) deallocate(cntr1)
    allocate( cntr1(ucx*ucy) )
    if( allocated(cntr2) ) deallocate(cntr2)
    allocate( cntr2(ucx*ucy) )
    if( allocated(rho) ) deallocate(rho)
    allocate( rho(ucx*ucy, 3) )

    !Calculate charge differences for bonds lying within a unit cell (to be adjusted for general patterns AB or BA as well as PBC)
    rho = 0 
    do j = 1, dim 
        cntr  = 0 
        cntr1 = 0 
        cntr2 = 0 
        do i = 1, nbx !run over all x bonds
            if ( modulo( bsites(1,i), 2 ) == 0 ) cycle !Don't count x bonds between unit cells 
            cntr = cntr + 1
            if ( btest( basis(j), bsites(1,i) - 1 ) ) cntr1(cntr) = cntr1(cntr) + 1 
            if ( btest( basis(j), bsites(2,i) - 1 ) ) cntr2(cntr) = cntr2(cntr) + 1 
            ! rho(cntr) = rho(cntr) + cntr(cntr) * psi(j)**2 
            if (pat == 'AB' ) then !A = cntr1, B = cntr2
                rho(cntr,1) = rho(cntr,1) + (cntr1(cntr) - cntr2(cntr)) * abs(psi(j))**2 !<n_A> - <n_B>
                rho(cntr,2) = rho(cntr,2) + cntr1(cntr) * abs(psi(j))**2 !<n_A>
                rho(cntr,3) = rho(cntr,3) + cntr2(cntr) * abs(psi(j))**2 !<n_B>
            else if (pat == 'BA' ) then !A = cntr2, B = cntr1
                rho(cntr,1) = rho(cntr,1) + (cntr2(cntr) - cntr1(cntr)) * abs(psi(j))**2 !<n_A> - <n_B>
                rho(cntr,2) = rho(cntr,2) + cntr2(cntr) * abs(psi(j))**2 !<n_A>
                rho(cntr,3) = rho(cntr,3) + cntr1(cntr) * abs(psi(j))**2 !<n_B>
            end if 
        end do 
    end do 

    if( allocated(cntrtot1) ) deallocate(cntrtot1)
    allocate( cntrtot1(nbonds) )
    if( allocated(cntrtot2) ) deallocate(cntrtot2)
    allocate( cntrtot2(nbonds) )
    if( allocated(rhotot) ) deallocate(rhotot)
    allocate( rhotot(nbonds, 3) )

    rhotot = 0 
    do j = 1, dim  
        cntrtot1 = 0 
        cntrtot2 = 0 
        do i = 1, nbonds 
            ! Determine sign for each sublattice (to be fixed to A=+1 and B=-1)
            ! if ( modulo( bsites(1,i), 2 ) == 1 ) then 
            !     sgn = 1
            ! else if ( modulo( bsites(1,i), 2 ) == 0 ) then 
            !     sgn = - 1
            ! end if   
            ! if( btest(basis(j), bsites(1,i) - 1) .and. (.not. btest(basis(j), bsites(2,i) - 1) ) ) cntrtot1(i) = cntrtot1(i) + 1 
            ! if( btest(basis(j), bsites(2,i) - 1) .and. (.not. btest(basis(j), bsites(1,i) - 1) ) ) cntrtot2(i) = cntrtot2(i) + 1 
            if( btest(basis(j), bsites(1,i) - 1) ) cntrtot1(i) = cntrtot1(i) + 1 
            if( btest(basis(j), bsites(2,i) - 1) ) cntrtot2(i) = cntrtot2(i) + 1 

            
            ! rhotot(i,1) = rhotot(i,1) + cntrtot(i,1) * psi(j)**2 
            if (pat == 'AB' .and. modulo( bsites(1,i), 2 ) == 1 ) then !A = cntrtot1, B = cntrtot2
                rhotot(i,1) = rhotot(i,1) + (cntrtot1(i) - cntrtot2(i)) * abs(psi(j))**2 !<n_A> - <n_B>
                rhotot(i,2) = rhotot(i,2) + cntrtot1(i) * abs(psi(j))**2 !<n_A>
                rhotot(i,3) = rhotot(i,3) + cntrtot2(i) * abs(psi(j))**2 !<n_B>
            else if (pat == 'AB' .and. modulo( bsites(1,i), 2 ) == 0 ) then !A = cntrtot2 B = cntrtot1
                rhotot(i,1) = rhotot(i,1) + (cntrtot2(i) - cntrtot1(i)) * abs(psi(j))**2 !<n_A> - <n_B>
                rhotot(i,2) = rhotot(i,2) + cntrtot2(i) * abs(psi(j))**2 !<n_A>
                rhotot(i,3) = rhotot(i,3) + cntrtot1(i) * abs(psi(j))**2 !<n_B>
            else if (pat == 'BA' .and. modulo( bsites(1,i), 2 ) == 1 ) then !A = cntrtot2 B = cntrtot1
                rhotot(i,1) = rhotot(i,1) + (cntrtot2(i) - cntrtot1(i)) * abs(psi(j))**2 !<n_A> - <n_B>
                rhotot(i,2) = rhotot(i,2) + cntrtot2(i) * abs(psi(j))**2 !<n_A>
                rhotot(i,3) = rhotot(i,3) + cntrtot1(i) * abs(psi(j))**2 !<n_B>
            else if (pat == 'BA' .and. modulo( bsites(1,i), 2 ) == 0 ) then !A = cntrtot1, B = cntrtot2
                rhotot(i,1) = rhotot(i,1) + (cntrtot1(i) - cntrtot2(i)) * abs(psi(j))**2 !<n_A> - <n_B>
                rhotot(i,2) = rhotot(i,2) + cntrtot1(i) * abs(psi(j))**2 !<n_A>
                rhotot(i,3) = rhotot(i,3) + cntrtot2(i) * abs(psi(j))**2 !<n_B>
            end if 
        end do 

    end do 

    if( allocated(cntr1) ) deallocate(cntr1)
    if( allocated(cntr2) ) deallocate(cntr2)
    if( allocated(cntrtot1) ) deallocate(cntrtot1)
    if( allocated(cntrtot2) ) deallocate(cntrtot2)

    call save(dir, parameters, append, unit, ucx, ucy, nbonds, rho, rhotot)
    return 

end subroutine cdw_c 

subroutine gsdeg_d(dir, rvec, unit, parameters, dim, nev, nest, thresh, energies, states, ndeg, gs)
    implicit none
    
    logical, intent(in) :: rvec
    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: unit, nev, nest
    double precision, intent(in) :: thresh 
    double precision, intent(in) :: energies(nev)
    double precision, intent(in) :: states(dim, nest)
    character(len=*), intent(in) :: dir, parameters 
    integer, intent(out) :: ndeg 
    double precision, allocatable, intent(out) :: gs(:)
    
    integer :: j 
    character*300 :: file 

    ndeg = 1
            
    if(rvec) then 
        if(allocated(gs)) deallocate(gs)
        allocate(gs(dim))
        gs = states(1:dim,1)
    end if 
    do j = 2, nest 
        if ( abs(energies(j) - energies(1)) <= thresh ) then 
            ! print*, 'dE=', abs(energies(j) - energies(1))
            ! if (j == 2 ) print*, 'E(1)',energies(1) 
            ! print*, 'E(j)=',energies(j)
            ndeg = ndeg + 1
            if(rvec) gs = gs + states(1:dim,j)
        else 
            exit 
        end if 
    end do 
    if ( nest .le. nev .and. abs(energies(nest + 1) - energies(1)) <= thresh ) print*,'Potentially not full degeneracy has been captured. Increase "nest".'
    if ( nest .le. ndeg ) print*,'Potentially not full degeneracy has been captured. Increase "nest".'
    ndeg = min(ndeg, nest)
    if(rvec) gs = gs/sqrt(dble(ndeg))
    print*,'Ground state degeneracy = ',ndeg 
   
    file = dir//"gs_deg_"//parameters
    file = trim_name(file)

    open(unit,file = file)
    write(unit,*) ndeg 
    close(unit)

    return 
   
end subroutine gsdeg_d

subroutine gsdeg_c(dir, rvec, unit, parameters, dim, nev, nest, thresh, energies, states, ndeg, gs)
    implicit none

    logical, intent(in) :: rvec
    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: unit, nev, nest
    double precision, intent(in) :: thresh 
    double precision, intent(in) :: energies(nev)
    double complex, intent(in) :: states(dim, nest)
    character(len=*), intent(in) :: dir, parameters 
    integer, intent(out) :: ndeg 
    double complex, allocatable, intent(out) :: gs(:)
    
    integer :: j 
    character*300 :: file 
    
    ndeg = 1

    if(rvec) then 
        if(allocated(gs)) deallocate(gs)
        allocate(gs(dim))
        gs = states(1:dim,1)
    end if 
    do j = 2, nest 
        if ( abs(energies(j) - energies(1)) <= thresh ) then 
            ! print*, 'dE=', abs(energies(j) - energies(1))
            ! if (j == 2 ) print*, 'E(1)',energies(1) 
            ! print*, 'E(j)=',energies(j)
            ndeg = ndeg + 1
            if(rvec) gs = gs + states(1:dim,j)
        else 
            exit 
        end if 
    end do 
    if ( nest .le. nev .and. abs(energies(nest + 1) - energies(1)) <= thresh ) print*,'Potentially not full degeneracy has been captured. Increase "nest".'
    if ( nest .le. ndeg ) print*,'Potentially not full degeneracy has been captured. Increase "nest".'
    ndeg = min(ndeg, nest)
    if(rvec) gs = gs/sqrt(dble(ndeg))
    print*,'Ground state degeneracy = ',ndeg 
   
    file = dir//"gs_deg_"//parameters
    file = trim_name(file)

    open(unit,file = file)
    write(unit,*) ndeg 
    close(unit)

    return 
   
end subroutine gsdeg_c

subroutine qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, states, ndeg, gs)
    implicit none
    
    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: unit, nev, nest
    double precision, intent(in) :: energies(nev)
    double precision, intent(in) :: states(dim, nest)
    character(len=*), intent(in) :: dir, parameters 
    integer, intent(out) :: ndeg 
    double precision, allocatable, intent(out) :: gs(:)
    
    integer :: j
    double precision :: dE(nest-1)
    character*300 :: file 

    ndeg = 1
            
    if(allocated(gs)) deallocate(gs)
    allocate(gs(dim))
    gs = states(1:dim,1)
    
    do j = 2, nest 
        dE(j-1) = abs(energies(j) - energies(j-1)) 
        ! print*,'j, dE(j)',j-1, dE(j-1)
    end do
    ndeg = maxloc(real(dE), dim=1) 
    
    ! print*,'dE(ndeg-1)',dE(ndeg-1)
    ! print*,'dE(ndeg)',dE(ndeg)
    
    if(ndeg > 1) then 
        do j = 2, ndeg 
            gs = gs + states(1:dim,j)
        end do 
    end if 
    gs = gs/sqrt(dble(ndeg))
    print*,'Ground state quasi-degeneracy = ',ndeg 
   
    file =dir//"quasi_gs_deg_"//parameters
    file = trim_name(file)

    open(unit,file = file)
    write(unit,*) ndeg 
    close(unit)

    return 
   
end subroutine qgsdeg_d

subroutine qgsdeg_c(dir, unit, parameters, dim, nev, nest, energies, states, ndeg, gs)
    implicit none
    
    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: unit, nev, nest
    double precision, intent(in) :: energies(nev)
    double complex, intent(in) :: states(dim, nest)
    character(len=*), intent(in) :: dir, parameters 
    integer, intent(out) :: ndeg 
    double complex, allocatable, intent(out) :: gs(:)
    
    integer :: j
    double precision :: dE(nest-1)
    character*300 :: file 

    ndeg = 1
            
    if(allocated(gs)) deallocate(gs)
    allocate(gs(dim))
    gs = states(1:dim,1)
    
    do j = 2, nest 
        dE(j-1) = abs(energies(j) - energies(j-1)) 
        ! print*,'j, dE(j)',j-1, dE(j-1)
    end do
    ndeg = maxloc(real(dE), dim=1) 
    
    ! print*,'dE(ndeg-1)',dE(ndeg-1)
    ! print*,'dE(ndeg)',dE(ndeg)
    
    if(ndeg > 1) then 
        do j = 2, ndeg 
            gs = gs + states(1:dim,j)
        end do 
    end if 
    gs = gs/sqrt(dble(ndeg))
    print*,'Ground state quasi-degeneracy = ',ndeg 
   
    file =dir//"quasi_gs_deg_"//parameters
    file = trim_name(file)

    open(unit,file = file)
    write(unit,*) ndeg 
    close(unit)

    return 
   
end subroutine qgsdeg_c

subroutine corrmat(ndeg, dim, psi, evals, mat)

    implicit none

    integer, intent(in) :: ndeg
    integer(kind=8), intent(in) :: dim
    double precision, intent(in) :: psi(dim, ndeg)
    double precision, allocatable, intent(out) :: evals(:), mat(:,:)

    integer :: i, j 

    if(allocated(evals)) deallocate(evals)
    if(allocated(mat)) deallocate(mat)
    allocate(evals(ndeg))
    allocate(mat(ndeg, ndeg))
    
    do i = 1, ndeg 
        do j = 1, ndeg 
            mat(i, j) = dot_product(psi(1:dim, i), psi(1:dim, j))
            print* ,mat(i, j), i, j, 'mat(i, j), i, j'
        end do 
    end do 

    call exactdiag(.True., ndeg, mat, evals)
    do i = 1, ndeg
        print* ,evals(i), i, 'evals(i), i'
    end do 

end subroutine corrmat

!----------------------------------------------!
!            Sublattice current                !
!----------------------------------------------!

subroutine currentcorrelations(dir, ti, tilted, nHel, tilt, threads, conf, degeneracy, full, feast, mkl, arpack, symmetrize, sites, l2, l1, ucx, ucy, k1, k2, id, par, rot, nDis, ndeg, refbonds, cntrA, cntrB, unit, dim, alattice, blattice, xyA, xyB, xtransl, ytransl, refl, c6, basis, v1, v2, dis, nnnVec, norm, eigstate, eigstate_dc)

    implicit none 
    
    integer, intent(in) :: ti, tilted, nHel, tilt, threads, conf, degeneracy, full, feast, mkl, arpack, symmetrize, sites, l2, l1, ucx, ucy, k1, k2, nDis, ndeg, refbonds, cntrA, cntrB, unit
    integer(kind=8), intent(in) :: dim
    integer, allocatable, intent(in) :: alattice(:,:), blattice(:,:), xyA(:,:), xyB(:,:), xtransl(:,:), ytransl(:,:), refl(:,:), c6(:) 
    integer(kind=8), allocatable, intent(in) :: basis(:)
    double precision, intent(in) :: v1, v2, dis, nnnVec(2,3), id, par(6), rot(5)
    double precision, allocatable, intent(in) :: norm(:), eigstate(:,:)
    double complex, allocatable, intent(in) :: eigstate_dc(:,:)
    character(len=*), intent(in) :: dir

    integer :: numthrds = 0, thread = 0, j = 0, nrefb = 0, l11 = 0, l22 = 0
    integer, save :: refbond = 0
    double precision :: current = 0.d0 
    double precision, allocatable :: bondcurrent(:)
    character :: sl*1, dirq*200
    character, save :: current_parameters*200
    logical :: append 
    !!$omp threadprivate(refbond, j, nrefb, current, bondcurrent, current_parameters, sl)
    !!$omp threadprivate(refbond, current_parameters)
    
    if(ti == 1 .and. tilted == 1) then 
        l11 = l2
        l22 = l1 
    else 
        l11 = ucx 
        l22 = ucy 
    end if 


    if((nDis > 1) .and. (dis .ne. 0.d0) .and. (conf > 1)) then 
        append = .true.
    else 
        append = .false.
    end if 

    if(degeneracy < 2 ) dirq = dir  
    if(degeneracy == 2 ) dirq = trim_name(dir//"QD_")
    dirq = trim_name(dirq)

    print*,'Calculate current current correlations...'
    print*,''
    if(refbonds > 0) then 
        nrefb = refbonds 
    else if(refbonds == 0) then 
        nrefb = cntrA
    end if 
    sl = "A"
    !!$ call omp_set_num_threads(threads)

    !!$omp parallel do default(firstprivate) shared(sites, particles, dim, basis) num_threads(threads)
    do j = 1, nrefb
        refbond = Alattice(4,j)
        call par2(j + 11, k1, k2, refbond, v1, v2, current_parameters)
        
        !!$ numthrds = omp_get_num_threads()
        !!$ thread = omp_get_thread_num()

        if( ndeg == 1 ) then  
            if(ti == 0) then 
                if(tilted == 0) call slcurrent(1.0d0, dim, sites, basis, Alattice, xyA, cntrA, j, eigstate(1:dim,1), nnnVec, bondcurrent, current)
                if(tilted == 1) call current_irrep(nHel, tilt, 1.0d0, dim, sites, 1, 1, ti, symmetrize, id, par, rot, basis, Alattice, xtransl, ytransl, refl, c6, cntrA, j, eigstate(1:dim,1), nnnVec, norm, bondcurrent, current)
            else if(ti == 1 .and. k1 == 0 .and. k2 == 0) then 
                if(tilted == 1) call current_irrep(nHel, tilt, 1.0d0, dim, sites, l2, l1, ti, symmetrize, id, par, rot, basis, Alattice, xtransl, ytransl, refl, c6, cntrA, j, eigstate(1:dim,1), nnnVec, norm, bondcurrent, current)
            else if(ti == 1 .and. ((k1 .ne. 0) .or. (k2 .ne. 0))) then 
                call current_c(tilted, nHel, tilt, k1, k2, 1.0d0, dim, sites, l11, l22, basis, Alattice, xyA, xtransl, ytransl, cntrA, j, eigstate_dc(1:dim,1), nnnVec, norm, bondcurrent, current)
            end if 
        else if(ndeg > 1) then 
            if(ti == 0) then 
                if(tilted == 1) call current_irrep_deg(ndeg, nHel, tilt, 1.d0, dim, sites, 1, 1, ti, symmetrize, id, par, rot, basis, Alattice, xtransl, ytransl, refl, c6, cntrA, refbond, eigstate, nnnVec, norm, bondcurrent, current)
                if(tilted == 0) call current_irrep_deg_rect(ndeg, 1.0d0, dim, sites, 1, 1, ti, symmetrize, id, par, rot, basis, Alattice, xtransl, ytransl, refl, c6, cntrA, refbond, eigstate, nnnVec, norm, bondcurrent, current)
            else if(ti == 1 .and. k1 == 0 .and. k2 == 0) then 
                if(tilted == 1) then 
                    call current_irrep_deg(ndeg, nHel, tilt, 1.d0, dim, sites, l2, l1, ti, symmetrize, id, par, rot, basis, Alattice, xtransl, ytransl, refl, c6, cntrA, refbond, eigstate, nnnVec, norm, bondcurrent, current)
                    ! call jij_irrep_deg(ndeg, 1.0d0, dim, sites, l2, l1, nHel, symmetrize, id, par, rot, basis, Alattice, xtransl, ytransl, refl, c6, cntrA, eigstate, nnnVec, norm) !, bcurrent, current
                else if(tilted == 0) then 
                    call current_irrep_deg_rect(ndeg, 1.0d0, dim, sites, l2, l1, ti, symmetrize, id, par, rot, basis, Alattice, xtransl, ytransl, refl, c6, cntrA, refbond, eigstate, nnnVec, norm, bondcurrent, current)
                    call jij_irrep_deg_rect(ndeg, l2, l1, id, par, rot, basis, Alattice, xtransl, ytransl, refl, c6, cntrA, eigstate, nnnVec, norm)
                    
                end if 
                
            else if(ti == 1 .and. ((k1 .ne. 0) .or. (k2 .ne. 0))) then 
                call current_ics_c2(tilted, nHel, tilt, k1, k2, ndeg, 1.d0, dim, sites, l2, l1, basis, Alattice, xyA, xtransl, ytransl, cntrA, j, eigstate_dc(1:dim,1:ndeg), nnnVec, norm, bondcurrent, current)    
                call current_cs_c2(tilted, nHel, tilt, dirq, j + 11, current_parameters, sl, ndeg, k1, k2, 1.0d0, dim, sites, l2, l1, basis, Alattice, xyA, xtransl, ytransl, cntrA, j, eigstate_dc(1:dim,1:ndeg), nnnVec, norm)
                call currmat_c2(tilted, nHel, tilt, dirq, j + 11, current_parameters, sl, ndeg, k1, k2, 1.d0, dim, sites, l2, l1, basis, Alattice, xyA, xtransl, ytransl, cntrA, j, eigstate_dc(1:dim,1:ndeg), nnnVec, norm)
            end if 
        end if 
        
        write (* ,"(' A lattice, refbond ',i2,': Current=',f20.17)") refbond, current       
        call save(dir, sl, current_parameters, append, degeneracy, j + 11, full, feast, mkl, arpack, cntrA, current, bondcurrent)   
        
    end do !Loop over bonds 
    !!$omp end parallel do 
    
    
    if(refbonds > 0) then 
        nrefb = refbonds 
    else if(refbonds == 0) then 
        nrefb = cntrB
    end if 
    
    ! sl = "B"
    ! !!$omp parallel do default(firstprivate) shared(sites, particles, dim, basis) num_threads(threads)
    ! do j = 1, nrefb
        
    !     !!$ numthrds = omp_get_num_threads()
    !     !!$ thread = omp_get_thread_num()

    !     refbond = Blattice(4,j)
    !     call par2(unit, k1, k2, refbond, v1, v2, current_parameters)
    !     if( ndeg == 1 ) then  
    !         if(ti == 0) then 
    !             call slcurrent(1.0d0, dim, sites, basis, Blattice, xyB, cntrB, j, eigstate(1:dim,1), nnnVec, bondcurrent, current)
    !         else if(ti == 1 .and. k1 == 0 .and. k2 == 0) then 
    !             if(tilted == 1) call current_irrep(nHel, tilt, 1.0d0, dim, sites, l2, l1, ti, symmetrize, id, par, rot, basis, Blattice, xtransl, ytransl, refl, c6, cntrB, j, eigstate(1:dim,1), nnnVec, norm, bondcurrent, current)
    !         else if(ti == 1 .and. ((k1 .ne. 0) .or. (k2 .ne. 0))) then 
    !             call current_c(tilted, nHel, tilt, k1, k2, 1.0d0, dim, sites, l11, l22, basis, Blattice, xyB, xtransl, ytransl, cntrB, j, eigstate_dc(1:dim,1), nnnVec, norm, bondcurrent, current)
    !         end if 
    !     else if ( ndeg > 1 ) then 
    !         if(ti == 0) then 
    !             call current_ics(ndeg, 1.0d0, dim, sites, basis, Blattice, cntrB, j, eigstate(1:dim,1:ndeg), nnnVec, bondcurrent, current)
    !             call currmat(dirq, j + 11, current_parameters, sl, ndeg, 1.0d0, dim, sites, basis, Blattice, xyB, cntrB, j, eigstate(1:dim,1:ndeg), nnnVec)
    !             call current_cs(dirq, j + 11, current_parameters, sl, ndeg, 1.0d0, dim, sites, basis, Blattice, cntrB, j, eigstate(1:dim,1:ndeg), nnnVec)
    !         else if(ti == 1 .and. k1 == 0 .and. k2 == 0) then 
    !             if(tilted == 1) call current_irrep_deg(ndeg, nHel, tilt, 1.d0, dim, sites, l2, l1, ti, symmetrize, id, par, rot, basis, Blattice, xtransl, ytransl, refl, c6, cntrB, refbond, eigstate, nnnVec, norm, bondcurrent, current)
    !             ! if(tilted == 0) call jij_irrep_deg_rect(ndeg, 1.0d0, dim, sites, l2, l1, symmetrize, id, par, rot, basis, Blattice, xtransl, ytransl, refl, c6, cntrB, eigstate, nnnVec, norm)
    !             if(tilted == 0) call current_irrep_deg_rect(ndeg, 1.0d0, dim, sites, l2, l1, ti, symmetrize, id, par, rot, basis, Blattice, xtransl, ytransl, refl, c6, cntrB, refbond, eigstate, nnnVec, norm, bondcurrent, current)
    !             if(tilted == 1) call current_irrep_deg(ndeg, nHel, tilt, 1.d0, dim, sites, l2, l1, ti, symmetrize, id, par, rot, basis, Blattice, xtransl, ytransl, refl, c6, cntrB, refbond, eigstate, nnnVec, norm, bondcurrent, current)
    !             ! if(tilted == 0) call jij_irrep_deg_rect(ndeg, 1.0d0, dim, sites, l2, l1, symmetrize, id, par, rot, basis, Blattice, xtransl, ytransl, refl, c6, cntrB, eigstate, nnnVec, norm)
    !             if(tilted == 0) call current_irrep_deg_rect(ndeg, 1.0d0, dim, sites, l2, l1, ti, symmetrize, id, par, rot, basis, Blattice, xtransl, ytransl, refl, c6, cntrB, refbond, eigstate, nnnVec, norm, bondcurrent, current)
    !         else if(ti == 1 .and. ((k1 .ne. 0) .or. (k2 .ne. 0))) then 
    !             call current_ics_c2(tilted, nHel, tilt, k1, k2, ndeg, 1.d0, dim, sites, l2, l1, basis, Blattice, xyB, xtransl, ytransl, cntrB, j, eigstate_dc(1:dim,1:ndeg), nnnVec, norm, bondcurrent, current)    
    !             call current_cs_c2(tilted, nHel, tilt, dirq, j + 11, current_parameters, sl, ndeg, k1, k2, 1.0d0, dim, sites, l2, l1, basis, Blattice, xyB, xtransl, ytransl, cntrB, j, eigstate_dc(1:dim,1:ndeg), nnnVec, norm)
    !             call currmat_c2(tilted, nHel, tilt, dirq, j + 11, current_parameters, sl, ndeg, k1, k2, 1.d0, dim, sites, l2, l1, basis, Blattice, xyB, xtransl, ytransl, cntrB, j, eigstate_dc(1:dim,1:ndeg), nnnVec, norm)
    !         end if 
    !     end if 

    !     write (* ,"(' B lattice, refbond ',i2,': Current=',f20.17)") refbond, current                
    !     call save(dir, sl, current_parameters, append, degeneracy, unit, full, feast, mkl, arpack, cntrB, current, bondcurrent)

    ! end do !Loop over bonds 
    ! !!$omp end parallel do 

    return 
end subroutine currentcorrelations

subroutine slcurrent(t, dim, sites, basis, bsites, xy, nbonds, refbond, psi, nnnVec, bondcurrent, current)

    implicit none
    ! save 
    integer, intent(in) :: nbonds, refbond, sites
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    integer, intent(in) :: xy(4,nbonds)
    double precision, intent(in) :: t 
    double precision, intent(in) :: psi(dim), nnnVec(2, 3)
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bondcurrent(:)
    
    integer(kind=8) :: loc = 0, newst = 0
    integer :: i = 0, j = 0
    integer :: sign = 0
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec, coord1, coord2 
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bondcurrent)) deallocate(bondcurrent)
    if( allocated(psiprime) ) deallocate(psiprime)
    allocate(bondcurrent(nbonds))
    allocate(psiprime(dim))

    psiprime    = 0.d0
    refvec      = vecs(1:2,bsites(5,refbond))
    coord1(1:2) = int(0.5*xy(3,refbond))*vecs(1:2,1) + xy(4,refbond)*vecs(1:2,2) - (int(0.5*xy(1,refbond))*vecs(1:2,1) + xy(2,refbond)*vecs(1:2,2))
    refsite1    = bsites(1,refbond)
    refsite2    = bsites(2,refbond)      

    bondcurrent = 0.d0
    current = 0.d0     
    do i = 1, nbonds            
        bondcurrent(i) = 0 
        site1 = bsites(1,i)
        site2 = bsites(2,i)
        if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
            .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
    
        
        coord2(1:2) = int(0.5*xy(3,i))*vecs(1:2,1) + xy(4,i)*vecs(1:2,2) - (int(0.5*xy(1,i))*vecs(1:2,1) + xy(2,i)*vecs(1:2,2))

        phase    = bsites(3,i)   
        vec      = vecs(1:2,bsites(5,i))
        dist     = dot_product(vec, refvec)
        psiprime = 0.d0  
        do j = 1, dim
            if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                    ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                    newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1   
                    call findstate(dim, newst, basis, loc) 
                    if(loc <= 0) cycle   
                    parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                    parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                    psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j)
                else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                    ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                    newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2  
                    call findstate(dim, newst, basis, loc) 
                    if(loc <= 0) cycle   
                    parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                    parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                    psiprime(loc) = psiprime(loc) + (-1) * sign * (-(t)**2) * psi(j)
                end if     
            else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                    ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                    newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1   
                    call findstate(dim, newst, basis, loc) 
                    if(loc <= 0) cycle   
                    parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                    parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                    psiprime(loc) = psiprime(loc) + (-1) * sign * (-(t)**2) * psi(j)
                else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                    ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                    newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2  
                    call findstate(dim, newst, basis, loc) 
                    if(loc <= 0) cycle   
                    parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                    parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                    psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j)
                end if     
            end if 
        end do !Loop over basis states                  
        bondcurrent(i) = bondcurrent(i) + dist * dot_product(psi, psiprime)
        current = current + dble(phase) * bondcurrent(i)
    end do !Loop over bonds 
    
    if(allocated(psiprime)) deallocate(psiprime)

    return

end subroutine slcurrent

!LOC and J reversed (ti = 1, k1 = 0, k2 = 0, ndeg = 1) (for all clusters)
subroutine current_irrep(nHel, tilt, t, dim, sites, Lx, Ly, ti, symmetrize, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, refbond, psi, nnnVec, norm, bcurrent, current)

    implicit none
    integer, intent(in) :: nHel, tilt, nbonds, refbond, sites, Lx, Ly, ti, symmetrize, bsites(5,nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    integer(kind=8), intent(in) :: dim, basis(dim)
    double precision, intent(in) :: id, mir(6), rot(5), t, psi(dim), nnnVec(2, 3), norm(dim) 
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bcurrent(:)
    
    integer(kind=8) :: loc = 0, newst = 0, rep = 0, state = 0
    integer :: x, y, info, i, j, k, l1, l2
    integer :: sign = 1, signrep = 1, signstate = 1, signDir = 1, order = 0 
    integer :: site1, site2, refsite1, refsite2, phase, parity
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bcurrent)) deallocate(bcurrent)
    if(allocated(psiprime)) deallocate(psiprime)
    allocate(bcurrent(nbonds))
    allocate(psiprime(dim))

    refvec   = vecs(1:2,bsites(5,refbond))
    refsite1 = bsites(1,refbond)
    refsite2 = bsites(2,refbond)  
    current  = 0.d0
    bcurrent = 0.d0 
   
    if(symmetrize == 1 .and. ti == 1) then 
        order = 12
    else 
        order = 1
    end if 

    do i = 1, nbonds !Loop over bonds    
        psiprime = 0.0d0  
        site1    = bsites(1,i)
        site2    = bsites(2,i)
        phase    = bsites(3,i)   
        refsite1 = bsites(1,refbond)
        refsite2 = bsites(2,refbond)
        vec      = vecs(1:2,bsites(5,i))
        dist     = dot_product(vec, refvec)
        dist     = 1.d0

        if(site1 == refsite1 .or. site1 == refsite2) cycle
        if(site2 == refsite1 .or. site2 == refsite2) cycle !Exclude bonds which share a site with the reference bond

        do j = 1, dim !Loop over basis states 

            do k = 1, order 
                if(k == 1) then 
                    state = basis(j)
                else if(k <= 7) then 
                    call reflect(basis(j), basis(j), sites, refl(k - 1, 1:sites), signstate, info, state)
                    signstate = signstate * mir(k - 1)
                else if(8 <= k) then 
                    call c6n(basis(j), basis(j), sites, k - 7, c6, signstate, info, state)
                    signstate = signstate * rot(k - 7)
                end if 
                
                do x = 1, Lx
                    if(ti == 1) call xtranslate(state, nHel, sites, Lx, xtransl, signstate, state)

                    do y = 1, Ly
                        if(ti == 1) call ytranslate(state, nHel, tilt, sites, Lx, ytransl, signstate, state)
                        if(btest(state, refsite1-1) .and. .not.(btest(state, refsite2-1))) then !Forward hopping: 1 -> 2
                            newst = ibclr(ibset(state, refsite2-1), refsite1-1) !Create on refsite 2, annihilate on refsite 1   
                            parity = popcnt(ibits(state, refsite1, sites)) + popcnt(ibits(ibclr(state, refsite1-1), refsite2, sites))
                            signDir = 1 
                        else if(btest(state, refsite2-1) .and. .not.(btest(state, refsite1-1))) then !Backwards hopping: 2 -> 1 
                            newst = ibclr(ibset(state, refsite1-1), refsite2-1) !Create on refsite 1, annihilate on refsite 2
                            parity = popcnt(ibits(state, refsite2, sites)) + popcnt(ibits(ibclr(state, refsite2-1), refsite1, sites))
                            signDir = -1 
                        else 
                            cycle 
                        end if 

                        if(btest(state, site1-1) .and. .not.(btest(state, site2-1))) then !Forward forward hopping: 1 -> 2       
                            parity = parity + popcnt(ibits(newst, site1, sites)) + popcnt(ibits(ibclr(newst, site1-1), site2, sites))                     
                            newst = ibclr(ibset(newst, site2-1), site1-1) !Create on site 2, annihilate on site 1                             
                        else if(btest(state, site2-1) .and. .not.(btest(state, site1-1))) then !Backwards hopping: 2 -> 1
                            parity = parity + popcnt(ibits(newst, site2, sites)) + popcnt(ibits(ibclr(newst, site2-1), site1, sites)) 
                            newst = ibclr(ibset(newst, site1-1), site2-1) !Create on site 1, annihilate on site 2
                            signDir = signDir * (-1)
                        else 
                            cycle 
                        end if 

                        if(ti == 1) call representative_irrep(newst, sites, nHel, tilt, Lx, Ly, symmetrize, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                        if(ti == 0) rep = newst
                        call findstate(dim, rep, basis, loc) 
                        
                        if(loc <= 0) cycle   
                        if(abs(signstate) > 1) stop 'signstate > 1' 
                        if(abs(signrep) > 1)   stop 'signrep > 1' 
                        if(abs(signDir) > 1)   stop 'signDir > 1' 
                        
                        sign = signstate * signrep * signDir * (-1)**parity 
                        ! psiprime(j) = psiprime(j) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(loc)
                        psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(j)
                        
                    end do !y
                end do !x 
            end do !G 
            
        end do !j 
        bcurrent(i) = bcurrent(i) + dist * dot_product(psi, psiprime) / (Lx * Ly * order) !Bond current J_ij 
        current = current + dble(phase) * bcurrent(i)       
    end do !Loop over bonds    
    print*, current, 'current'
    deallocate(psiprime)

    return

end subroutine current_irrep

subroutine current_rect(t, dim, sites, Lx, Ly, ti, symmetrize, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, refbond, psi, nnnVec, norm, bcurrent, current)

    implicit none
    integer, intent(in) :: nbonds, refbond, sites, Lx, Ly, ti, symmetrize
    integer, intent(in) :: bsites(5,nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    integer(kind=8), intent(in) :: dim, basis(dim)
    double precision, intent(in) :: id, mir(6), rot(5), t, psi(dim), nnnVec(2, 3), norm(dim) 
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bcurrent(:)
    
    integer(kind=8) :: loc = 0, newst = 0, rep = 0, state = 0
    integer :: x, y, info, i, j, k, l1, l2
    integer :: sign = 1, signrep = 1, signstate = 1, signDir = 1, order = 0 
    integer :: site1, site2, refsite1, refsite2, phase, parity
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bcurrent)) deallocate(bcurrent)
    if(allocated(psiprime)) deallocate(psiprime)
    allocate(bcurrent(nbonds))
    allocate(psiprime(dim))

    refvec   = vecs(1:2,bsites(5,refbond))
    refsite1 = bsites(1,refbond)
    refsite2 = bsites(2,refbond)  
    current  = 0.d0
    bcurrent = 0.d0 
   
    if(symmetrize == 1 .and. ti == 1) then 
        order = 12
    else 
        order = 1
    end if 

    do i = 1, nbonds !Loop over bonds    
        psiprime = 0.0d0  
        site1    = bsites(1,i)
        site2    = bsites(2,i)
        phase    = bsites(3,i)   
        refsite1 = bsites(1,refbond)
        refsite2 = bsites(2,refbond)
        vec      = vecs(1:2,bsites(5,i))
        dist     = dot_product(vec, refvec)
        dist     = 1.d0

        if(site1 == refsite1 .or. site1 == refsite2) cycle
        if(site2 == refsite1 .or. site2 == refsite2) cycle !Exclude bonds which share a site with the reference bond

        do j = 1, dim !Loop over basis states 

            do k = 1, order 
                if(k == 1) then 
                    state = basis(j)
                else if(k <= 7) then 
                    call reflect(basis(j), basis(j), sites, refl(k - 1, 1:sites), signstate, info, state)
                    signstate = signstate * mir(k - 1)
                else if(8 <= k) then 
                    call c6n(basis(j), basis(j), sites, k - 7, c6, signstate, info, state)
                    signstate = signstate * rot(k - 7)
                end if 
                
                do x = 1, Lx
                    call xtranslate(state, Ly, sites, Lx, xtransl, signstate, state)

                    do y = 1, Ly
                        call ytranslate(state, Ly, 0, sites, Lx, ytransl, signstate, state)
                        if(btest(state, refsite1-1) .and. .not.(btest(state, refsite2-1))) then !Forward hopping: 1 -> 2
                            newst = ibclr(ibset(state, refsite2-1), refsite1-1) !Create on refsite 2, annihilate on refsite 1   
                            parity = popcnt(ibits(state, refsite1, sites)) + popcnt(ibits(ibclr(state, refsite1-1), refsite2, sites))
                            signDir = 1 
                        else if(btest(state, refsite2-1) .and. .not.(btest(state, refsite1-1))) then !Backwards hopping: 2 -> 1 
                            newst = ibclr(ibset(state, refsite1-1), refsite2-1) !Create on refsite 1, annihilate on refsite 2
                            parity = popcnt(ibits(state, refsite2, sites)) + popcnt(ibits(ibclr(state, refsite2-1), refsite1, sites))
                            signDir = -1 
                        else 
                            cycle 
                        end if 

                        if(btest(state, site1-1) .and. .not.(btest(state, site2-1))) then !Forward forward hopping: 1 -> 2       
                            parity = parity + popcnt(ibits(newst, site1, sites)) + popcnt(ibits(ibclr(newst, site1-1), site2, sites))                     
                            newst = ibclr(ibset(newst, site2-1), site1-1) !Create on site 2, annihilate on site 1                             
                        else if(btest(state, site2-1) .and. .not.(btest(state, site1-1))) then !Backwards hopping: 2 -> 1
                            parity = parity + popcnt(ibits(newst, site2, sites)) + popcnt(ibits(ibclr(newst, site2-1), site1, sites)) 
                            newst = ibclr(ibset(newst, site1-1), site2-1) !Create on site 1, annihilate on site 2
                            signDir = signDir * (-1)
                        else 
                            cycle 
                        end if 
                        if(ti == 1) call representative_irrep(newst, sites, Ly, 0, Lx, Ly, symmetrize, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                        if(ti == 0) rep = newst
                        call findstate(dim, rep, basis, loc) 
                        
                        if(loc <= 0) cycle   

                        sign = signstate * signrep * signDir * (-1)**parity 
                        ! psiprime(j) = psiprime(j) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(loc)
                        psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(j)
                        
                    end do !y
                end do !x 
            end do !G 
            
        end do !j 
        bcurrent(i) = bcurrent(i) + dist * dot_product(psi, psiprime) / (Lx * Ly * order) !Bond current J_ij 
        current = current + dble(phase) * bcurrent(i)       
    end do !Loop over bonds    

    deallocate(psiprime)

    return

end subroutine current_rect

subroutine current_irrep_deg(ndeg, nHel, tilt, t, dim, sites, Lx, Ly, ti, symmetrize, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, refbond, psi, nnnVec, norm, bcurrent, current)

    implicit none
    integer, intent(in) :: ndeg, nHel, tilt, nbonds, refbond, sites, Lx, Ly, ti, symmetrize, bsites(5,nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    integer(kind=8), intent(in) :: dim, basis(dim)
    double precision, intent(in) :: id, mir(6), rot(5), t, psi(dim, ndeg), nnnVec(2, 3), norm(dim) 
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bcurrent(:)
    integer(kind=8) :: loc = 0, newst = 0, rep = 0, state = 0
    integer :: x, y, info, i, j, k, l1, l2, nd1, nd2
    integer :: sign = 1, signrep = 1, signstate = 1, signDir = 1, order = 0 
    integer :: site1, site2, refsite1, refsite2, phase, parity
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:), currmat(:,:), bcurrmat(:,:,:), currents(:), evals(:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bcurrent)) deallocate(bcurrent)
    if(allocated(bcurrmat)) deallocate(bcurrmat)
    if(allocated(psiprime)) deallocate(psiprime)
    if(allocated(currents)) deallocate(currents)
    if(allocated(currmat)) deallocate(currmat)
    
    allocate(bcurrmat(nbonds, ndeg, ndeg))
    allocate(currmat(ndeg, ndeg))
    allocate(bcurrent(nbonds))
    allocate(currents(ndeg))
    allocate(psiprime(dim))

    refvec   = vecs(1:2,bsites(5,refbond))
    refsite1 = bsites(1,refbond)
    refsite2 = bsites(2,refbond)  
    bcurrent = 0.d0 
    currents = 0.d0 
    bcurrmat = 0.d0 
    current = 0.d0

    if(symmetrize == 1) then 
        order = 12
    else 
        order = 1
    end if 

    do nd1 = 1, ndeg 
        do nd2 = 1, ndeg 
            do i = 1, nbonds !Loop over bonds    
                psiprime = 0.0d0  
                site1    = bsites(1,i)
                site2    = bsites(2,i)
                phase    = bsites(3,i)   
                refsite1 = bsites(1,refbond)
                refsite2 = bsites(2,refbond)
                vec      = vecs(1:2,bsites(5,i))
                dist     = dot_product(vec, refvec)
                dist     = 1.d0

                if(site1 == refsite1 .or. site1 == refsite2) cycle
                if(site2 == refsite1 .or. site2 == refsite2) cycle !Exclude bonds which share a site with the reference bond

                do j = 1, dim !Loop over basis states 
                    do k = 1, order !Loop over point groups operators 
                        if(k == 1) then 
                            state = basis(j)
                        else if(k <= 7) then 
                            call reflect(basis(j), basis(j), sites, refl(k - 1, 1:sites), signstate, info, state)
                            signstate = signstate * mir(k - 1)
                        else if(8 <= k) then 
                            call c6n(basis(j), basis(j), sites, k - 7, c6, signstate, info, state)
                            signstate = signstate * rot(k - 7)
                        end if 
                        
                        do x = 1, Lx
                            if(ti == 1) call xtranslate(state, nHel, sites, Lx, xtransl, signstate, state)

                            do y = 1, Ly
                                if(ti == 1) call ytranslate(state, nHel, tilt, sites, Lx, ytransl, signstate, state)
                                if(btest(state, refsite1-1) .and. .not.(btest(state, refsite2-1))) then !Forward hopping: 1 -> 2
                                    newst = ibclr(ibset(state, refsite2-1), refsite1-1) !Create on refsite 2, annihilate on refsite 1   
                                    parity = popcnt(ibits(state, refsite1, sites)) + popcnt(ibits(ibclr(state, refsite1-1), refsite2, sites))
                                    signDir = 1 
                                else if(btest(state, refsite2-1) .and. .not.(btest(state, refsite1-1))) then !Backwards hopping: 2 -> 1 
                                    newst = ibclr(ibset(state, refsite1-1), refsite2-1) !Create on refsite 1, annihilate on refsite 2
                                    parity = popcnt(ibits(state, refsite2, sites)) + popcnt(ibits(ibclr(state, refsite2-1), refsite1, sites))
                                    signDir = -1 
                                else 
                                    cycle 
                                end if 

                                if(btest(state, site1-1) .and. .not.(btest(state, site2-1))) then !Forward forward hopping: 1 -> 2       
                                    parity = parity + popcnt(ibits(newst, site1, sites)) + popcnt(ibits(ibclr(newst, site1-1), site2, sites))                     
                                    newst = ibclr(ibset(newst, site2-1), site1-1) !Create on site 2, annihilate on site 1                             
                                else if(btest(state, site2-1) .and. .not.(btest(state, site1-1))) then !Backwards hopping: 2 -> 1
                                    parity = parity + popcnt(ibits(newst, site2, sites)) + popcnt(ibits(ibclr(newst, site2-1), site1, sites)) 
                                    newst = ibclr(ibset(newst, site1-1), site2-1) !Create on site 1, annihilate on site 2
                                    signDir = signDir * (-1)
                                else 
                                    cycle 
                                end if 

                                if(ti == 1) call representative_irrep(newst, sites, nHel, tilt, Lx, Ly, symmetrize, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                                if(ti == 0) rep = newst
                                call findstate(dim, rep, basis, loc) 
                                
                                if(loc <= 0) cycle   
                                
                                sign = signstate * signrep * signDir * (-1)**parity 
                                ! psiprime(j) = psiprime(j) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(loc)
                                if(ti == 1) psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(j, nd1)
                                if(ti == 0) psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j, nd1)
                            end do !y
                        end do !x 
                    end do !G 
                    
                end do !j 
                bcurrmat(i,nd2,nd1) = bcurrmat(i,nd2,nd1) + dble(phase) * dist * dot_product(psi(1:dim, nd2), psiprime)
                currmat(nd2,nd1) = currmat(nd2,nd1) + bcurrmat(i,nd2,nd1)
            end do !Loop over bonds    
        end do 
    end do !Degenerate states
    
    bcurrmat = bcurrmat / ndeg
    do i = 1, nbonds          
        call exactdiag(.True., ndeg, bcurrmat(i,1:ndeg,1:ndeg), evals)
        bcurrent(i) = sum(evals) 
    end do 

    currmat = currmat / ndeg 

    call exactdiag(.True., ndeg, currmat, currents)
    current = sum(currents)
    do i = 1, ndeg 
        print* ,currents(i), i, 'currents(i), i'
    end do
    deallocate(psiprime)

    return

end subroutine current_irrep_deg

subroutine current_irrep_deg_rect(ndeg, t, dim, sites, Lx, Ly, ti, symmetrize, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, refbond, psi, nnnVec, norm, bcurrent, current)

    implicit none
    integer, intent(in) :: ndeg, nbonds, refbond, sites, Lx, Ly, ti, symmetrize
    integer, intent(in) :: bsites(5,nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    integer(kind=8), intent(in) :: dim, basis(dim)
    double precision, intent(in) :: id, mir(6), rot(5), t, psi(dim, ndeg), nnnVec(2, 3), norm(dim) 
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bcurrent(:)
    integer(kind=8) :: loc = 0, newst = 0, rep = 0, state = 0
    integer :: x, y, info, i, j, k, l1, l2, nd1, nd2
    integer :: sign = 1, signrep = 1, signstate = 1, signDir = 1, order = 0 
    integer :: site1, site2, refsite1, refsite2, phase, parity
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:), currmat(:,:), bcurrmat(:,:,:), currents(:), evals(:), gs(:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bcurrent)) deallocate(bcurrent)
    if(allocated(bcurrmat)) deallocate(bcurrmat)
    if(allocated(psiprime)) deallocate(psiprime)
    if(allocated(currents)) deallocate(currents)
    if(allocated(currmat)) deallocate(currmat)
    
    allocate(bcurrmat(nbonds, ndeg, ndeg))
    allocate(currmat(ndeg, ndeg))
    allocate(bcurrent(nbonds))
    allocate(currents(ndeg))
    allocate(psiprime(dim))

    refvec   = vecs(1:2,bsites(5,refbond))
    refsite1 = bsites(1,refbond)
    refsite2 = bsites(2,refbond)  
    bcurrent = 0.d0 
    currents = 0.d0 
    bcurrmat = 0.d0 
    current = 0.d0

    if(symmetrize == 1) then 
        order = 12
    else 
        order = 1
    end if 

    do nd1 = 1, ndeg 
        do nd2 = 1, ndeg 
            do i = 1, nbonds !Loop over bonds    
                psiprime = 0.0d0  
                site1    = bsites(1,i)
                site2    = bsites(2,i)
                phase    = bsites(3,i)   
                refsite1 = bsites(1,refbond)
                refsite2 = bsites(2,refbond)
                vec      = vecs(1:2,bsites(5,i))
                dist     = dot_product(vec, refvec)
                dist     = 1.d0

                if(site1 == refsite1 .or. site1 == refsite2) cycle
                if(site2 == refsite1 .or. site2 == refsite2) cycle !Exclude bonds which share a site with the reference bond

                do j = 1, dim !Loop over basis states 
                    do k = 1, order !Loop over point groups operators 
                        if(k == 1) then 
                            state = basis(j)
                        else if(k <= 7) then 
                            call reflect(basis(j), basis(j), sites, refl(k - 1, 1:sites), signstate, info, state)
                            signstate = signstate * mir(k - 1)
                        else if(8 <= k) then 
                            call c6n(basis(j), basis(j), sites, k - 7, c6, signstate, info, state)
                            signstate = signstate * rot(k - 7)
                        end if 
                        
                        do x = 1, Lx
                            if(ti == 1) call xtranslate(state, Ly, sites, Lx, xtransl, signstate, state)

                            do y = 1, Ly
                                if(ti == 1) call ytranslate(state, Ly, 0, sites, Lx, ytransl, signstate, state)
                                if(btest(state, refsite1-1) .and. .not.(btest(state, refsite2-1))) then !Forward hopping: 1 -> 2
                                    newst   = ibclr(ibset(state, refsite2-1), refsite1-1) !Create on refsite 2, annihilate on refsite 1   
                                    parity  = popcnt(ibits(state, refsite1, sites)) + popcnt(ibits(ibclr(state, refsite1-1), refsite2, sites))
                                    signDir = 1 
                                else if(btest(state, refsite2-1) .and. .not.(btest(state, refsite1-1))) then !Backwards hopping: 2 -> 1 
                                    newst   = ibclr(ibset(state, refsite1-1), refsite2-1) !Create on refsite 1, annihilate on refsite 2
                                    parity  = popcnt(ibits(state, refsite2, sites)) + popcnt(ibits(ibclr(state, refsite2-1), refsite1, sites))
                                    signDir = -1 
                                else 
                                    cycle 
                                end if 

                                if(btest(state, site1-1) .and. .not.(btest(state, site2-1))) then !Forward forward hopping: 1 -> 2       
                                    parity = parity + popcnt(ibits(newst, site1, sites)) + popcnt(ibits(ibclr(newst, site1-1), site2, sites))                     
                                    newst  = ibclr(ibset(newst, site2-1), site1-1) !Create on site 2, annihilate on site 1                             
                                else if(btest(state, site2-1) .and. .not.(btest(state, site1-1))) then !Backwards hopping: 2 -> 1
                                    parity  = parity + popcnt(ibits(newst, site2, sites)) + popcnt(ibits(ibclr(newst, site2-1), site1, sites)) 
                                    newst   = ibclr(ibset(newst, site1-1), site2-1) !Create on site 1, annihilate on site 2
                                    signDir = signDir * (-1)
                                else 
                                    cycle 
                                end if 

                                if(ti == 1) call representative_irrep(newst, sites, Ly, 0, Lx, Ly, symmetrize, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                                if(ti == 0) rep = newst
                                call findstate(dim, rep, basis, loc) 
                                
                                if(loc <= 0) cycle   

                                sign = signstate * signrep * signDir * (-1)**parity 
                                if(ti == 1) psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(j, nd1)
                                if(ti == 0) psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j, nd1)
                                ! psiprime(j) = psiprime(j) + sign * (-(t)**2) * sqrt(norm(j)/norm(loc)) * psi(loc, nd1)
                                
                            end do !y
                        end do !x 
                    end do !G 
                    
                end do !j 
                bcurrmat(i,nd2,nd1) = bcurrmat(i,nd2,nd1) + dble(phase) * dist * dot_product(psi(1:dim, nd2), psiprime)
                currmat(nd2,nd1) = currmat(nd2,nd1) + bcurrmat(i,nd2,nd1)
            end do !Loop over bonds    
        end do 
    end do !Degenerate states
    
    bcurrmat = bcurrmat / ndeg
    do i = 1, nbonds          
        call exactdiag(.True., ndeg, bcurrmat(i,1:ndeg,1:ndeg), evals)
        bcurrent(i) = sum(evals) 
    end do 

    currmat = currmat / ndeg 
    ! do i = 1, ndeg 
    !     do j = 1, ndeg  
    !         print* ,currmat(i, j), i, j, '(JJ)_ij, i, j'      
    !     end do 
    !     print*, '------------------'
    ! end do 
    ! pause 
    call exactdiag(.True., ndeg, currmat, currents)
    current = sum(currents)
    if(allocated(gs)) deallocate(gs)
    allocate(gs(dim))
    gs = 0.d0 
    ! do i = 1, ndeg 
    !     print* ,currents(i), i, 'currents(i), i'      
    ! end do 
    
    
    ! do i = 1, ndeg 
        !     do j = 1, ndeg 
        !         print* ,norm2(currmat(1:ndeg,i)), 'norm2(currmat(j,i))'
        !         print* ,norm2(psi(1:dim,i)), 'norm2(psi(1:dim,i))'
        !         gs = gs + currmat(j, i) * psi(1:dim, j)
        !     end do 
        !     ! gs = gs/ norm2(gs)
        !     print* ,norm2(gs), 'norm2(gs)'
        !     call current_rect(t, dim, sites, Lx, Ly, ti, symmetrize, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, refbond, gs, nnnVec, norm, bcurrent, current)
        !     print* ,current, 'current'
        !     gs = 0.d0 
        !     pause     
    ! end do 
    deallocate(psiprime)

    return

end subroutine current_irrep_deg_rect

subroutine jij_irrep_deg_rect(ndeg, Lx, Ly, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, psi, nnnVec, norm)

    use variables_tmi

    implicit none
    
    integer, intent(in) :: ndeg, nbonds, Lx, Ly
    integer, intent(in) :: bsites(5,nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    integer(kind=8), intent(in) :: basis(dim)
    double precision, intent(in) :: id, mir(6), rot(5), psi(dim, ndeg), nnnVec(2, 3), norm(dim) 
    integer(kind=8) :: loc = 0, newst = 0, rep = 0, state = 0
    integer :: x, y, info, i, j, k, l1, l2, nd1, nd2
    integer :: sign = 1, signrep = 1, signstate = 1, signDir = 1, order = 0 
    integer :: site1, site2, phase, parity
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: evals(:), bcurrent(:,:)
    double complex, allocatable :: psiprime(:), bcurrmat(:,:,:), mat(:,:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

   
    if(allocated(bcurrmat)) deallocate(bcurrmat)
    if(allocated(bcurrent)) deallocate(bcurrent)
    if(allocated(psiprime)) deallocate(psiprime)
    
    allocate(bcurrmat(nbonds, ndeg, ndeg))
    allocate(bcurrent(nbonds, ndeg))
    allocate(psiprime(dim))
  
    bcurrmat = (0.d0,0.d0) 
    bcurrent = 0.d0 

    if(symmetrize == 1) then 
        order = 12
    else 
        order = 1
    end if 

    do nd1 = 1, ndeg 
        do nd2 = 1, ndeg 
            do i = 1, nbonds !Loop over bonds    
                psiprime = 0.0d0  
                site1    = bsites(1,i)
                site2    = bsites(2,i)
                phase    = bsites(3,i)   
                vec      = vecs(1:2,bsites(5,i))
                dist     = 1.d0

                do j = 1, dim !Loop over basis states 
                    do k = 1, order !Loop over point groups operators 
                        if(k == 1) then !Identity
                            state = basis(j)
                        else if(k <= 7) then !Reflections
                            call reflect(basis(j), basis(j), sites, refl(k - 1, 1:sites), signstate, info, state)
                            signstate = signstate * mir(k - 1)
                        else if(8 <= k) then !Rotations
                            call c6n(basis(j), basis(j), sites, k - 7, c6, signstate, info, state)
                            signstate = signstate * rot(k - 7)
                        end if 
                        
                        do x = 1, Lx !X-translations
                            call xtranslate(state, Ly, sites, Lx, xtransl, signstate, state)

                            do y = 1, Ly !Y-translations
                                call ytranslate(state, Ly, 0, sites, Lx, ytransl, signstate, state)

                                if(btest(state, site1-1) .and. .not.(btest(state, site2-1))) then !Forward forward hopping: 1 -> 2       
                                    parity  = popcnt(ibits(state, site1, sites)) + popcnt(ibits(ibclr(state, site1-1), site2, sites))                     
                                    newst   = ibclr(ibset(state, site2-1), site1-1) !Create on site 2, annihilate on site 1                             
                                    signDir = 1
                                else if(btest(state, site2-1) .and. .not.(btest(state, site1-1))) then !Backwards hopping: 2 -> 1
                                    parity  = popcnt(ibits(state, site2, sites)) + popcnt(ibits(ibclr(state, site2-1), site1, sites)) 
                                    newst   = ibclr(ibset(state, site1-1), site2-1) !Create on site 1, annihilate on site 2
                                    signDir = -1
                                else 
                                    cycle 
                                end if 

                                call representative_irrep(newst, sites, Ly, 0, Lx, Ly, symmetrize, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                                call findstate(dim, rep, basis, loc) 
                                
                                if(loc <= 0) cycle   
                                                                
                                sign = signstate * signrep * signDir * (-1)**parity 
                                
                                bcurrmat(i, nd2, nd1) = bcurrmat(i, nd2, nd1) + sign * ii * sqrt(norm(loc)/norm(j)) * psi(j, nd1) * psi(loc, nd2)
                                ! bcurrmat(i, nd2, nd1) = bcurrmat(i, nd2, nd1) + dble(phase) * dist * sign * t * ii * sqrt(norm(loc)/norm(j)) * psi(j, nd1) * psi(loc, nd2)
                                ! psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(j, nd1)
                                ! psiprime(j) = psiprime(j) + sign * (-(t)**2) * sqrt(norm(j)/norm(loc)) * psi(loc, nd1)
                                
                            end do !y
                        end do !x 
                    end do !G 
                    
                end do !j 
                ! bcurrmat(i,nd2,nd1) = bcurrmat(i,nd2,nd1) + dble(phase) * dist * dot_product(psi(1:dim, nd2), psiprime)
            end do !Loop over bonds    
        end do 
    end do !Degenerate states

    if(allocated(mat)) deallocate(mat)
    allocate(mat(ndeg, ndeg))
    ! bcurrmat = bcurrmat / ndeg
    do i = 1, nbonds          
        mat = bcurrmat(i,1:ndeg,1:ndeg)
        call exactdiag_c(.True., ndeg, mat, evals)
        bcurrent(i, 1:ndeg) = evals
        
        ! do j = 1, ndeg 
        !     print* ,evals(j), j, i, 'evals(state), state, bond'
        ! end do 
        
        ! do j = 1, ndeg 
        !     print*, 'state'
        !     do k = 1, ndeg 
        !         print* ,mat(k,j)
        !     end do 
        ! end do 
        ! print* ,sum(evals), 'sum(evals)'
        ! print* , ''    
        if(ndeg == 2) pause 
        mat = (0.d0, 0.d0)  
    end do
    ! print* ,prms, 'prms'
    call saveD("A", 11, "output/", "Jij_", prms, nbonds, ndeg, 1.d0, bcurrent, 1) 
    if(ndeg == 2) pause 

    deallocate(psiprime)
    
    
    return

end subroutine jij_irrep_deg_rect

! subroutine jij_irrep_deg_rect(ndeg, t, dim, sites, Lx, Ly, symmetrize, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, psi, nnnVec, norm) !, bcurrent, current

    !     implicit none
        
    !     integer, intent(in) :: ndeg, nbonds, sites, Lx, Ly, symmetrize
    !     integer, intent(in) :: bsites(5,nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    !     integer(kind=8), intent(in) :: dim, basis(dim)
    !     double precision, intent(in) :: id, mir(6), rot(5), t, psi(dim, ndeg), nnnVec(2, 3), norm(dim) 
    !     integer(kind=8) :: loc = 0, newst = 0, rep = 0, state = 0
    !     integer :: x, y, info, i, j, k, l1, l2, nd1, nd2
    !     integer :: sign = 1, signrep = 1, signstate = 1, signDir = 1, order = 0 
    !     integer :: site1, site2, phase, parity
    !     double precision :: dist = 0.d0
    !     double precision, dimension(2)   :: vec
    !     double precision, dimension(2,3) :: vecs
    !     double precision, allocatable :: evals(:), bcurrent(:,:)
    !     double complex, allocatable :: psiprime(:), bcurrmat(:,:,:), mat(:,:)

    !     !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bcurrent, psiprime, psiprimeij, psiprimeprime)

    !     vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    !     vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    !     vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    
    !     if(allocated(bcurrmat)) deallocate(bcurrmat)
    !     if(allocated(bcurrent)) deallocate(bcurrent)
    !     if(allocated(psiprime)) deallocate(psiprime)
        
    !     allocate(bcurrmat(nbonds, ndeg, ndeg))
    !     allocate(bcurrent(nbonds, ndeg))
    !     allocate(psiprime(dim))
    
    !     bcurrmat = (0.d0,0.d0) 
    !     bcurrent = 0.d0 

    !     if(symmetrize == 1) then 
    !         order = 12
    !     else 
    !         order = 1
    !     end if 

    !     do nd1 = 1, ndeg 
    !         do nd2 = 1, ndeg 
    !             do i = 1, nbonds !Loop over bonds    
    !                 psiprime = 0.0d0  
    !                 site1    = bsites(1,i)
    !                 site2    = bsites(2,i)
    !                 phase    = bsites(3,i)   
    !                 vec      = vecs(1:2,bsites(5,i))
    !                 dist     = 1.d0

    !                 do j = 1, dim !Loop over basis states 
    !                     do k = 1, order !Loop over point groups operators 
    !                         if(k == 1) then !Identity
    !                             state = basis(j)
    !                         else if(k <= 7) then !Reflections
    !                             call reflect(basis(j), basis(j), sites, refl(k - 1, 1:sites), signstate, info, state)
    !                             signstate = signstate * mir(k - 1)
    !                         else if(8 <= k) then !Rotations
    !                             call c6n(basis(j), basis(j), sites, k - 7, c6, signstate, info, state)
    !                             signstate = signstate * rot(k - 7)
    !                         end if 
                            
    !                         do x = 1, Lx !X-translations
    !                             call xtranslate(state, Ly, sites, Lx, xtransl, signstate, state)

    !                             do y = 1, Ly !Y-translations
    !                                 call ytranslate(state, Ly, 0, sites, Lx, ytransl, signstate, state)

    !                                 if(btest(state, site1-1) .and. .not.(btest(state, site2-1))) then !Forward forward hopping: 1 -> 2       
    !                                     parity  = popcnt(ibits(state, site1, sites)) + popcnt(ibits(ibclr(state, site1-1), site2, sites))                     
    !                                     newst   = ibclr(ibset(state, site2-1), site1-1) !Create on site 2, annihilate on site 1                             
    !                                     signDir = 1
    !                                 else if(btest(state, site2-1) .and. .not.(btest(state, site1-1))) then !Backwards hopping: 2 -> 1
    !                                     parity  = popcnt(ibits(state, site2, sites)) + popcnt(ibits(ibclr(state, site2-1), site1, sites)) 
    !                                     newst   = ibclr(ibset(state, site1-1), site2-1) !Create on site 1, annihilate on site 2
    !                                     signDir = -1
    !                                 else 
    !                                     cycle 
    !                                 end if 

    !                                 call representative_irrep(newst, sites, Ly, 0, Lx, Ly, symmetrize, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
    !                                 call findstate(dim, rep, basis, loc) 
                                    
    !                                 if(loc <= 0) cycle   
                                                                    
    !                                 sign = signstate * signrep * signDir * (-1)**parity 
                                    
    !                                 bcurrmat(i, nd2, nd1) = bcurrmat(i, nd2, nd1) + sign * ii * sqrt(norm(loc)/norm(j)) * psi(j, nd1) * psi(loc, nd2)
    !                                 ! bcurrmat(i, nd2, nd1) = bcurrmat(i, nd2, nd1) + dble(phase) * dist * sign * t * ii * sqrt(norm(loc)/norm(j)) * psi(j, nd1) * psi(loc, nd2)
    !                                 ! psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(j, nd1)
    !                                 ! psiprime(j) = psiprime(j) + sign * (-(t)**2) * sqrt(norm(j)/norm(loc)) * psi(loc, nd1)
                                    
    !                             end do !y
    !                         end do !x 
    !                     end do !G 
                        
    !                 end do !j 
    !                 ! bcurrmat(i,nd2,nd1) = bcurrmat(i,nd2,nd1) + dble(phase) * dist * dot_product(psi(1:dim, nd2), psiprime)
    !             end do !Loop over bonds    
    !         end do 
    !     end do !Degenerate states

    !     if(allocated(mat)) deallocate(mat)
    !     allocate(mat(ndeg, ndeg))
    !     ! bcurrmat = bcurrmat / ndeg
    !     do i = 1, nbonds          
    !         mat = bcurrmat(i,1:ndeg,1:ndeg)
    !         call exactdiag_c(.True., ndeg, mat, evals)
    !         bcurrent(i, 1:ndeg) = evals
            
    !         do j = 1, ndeg 
    !             print* ,evals(j), j, i, 'evals(state), state, bond'
    !         end do 
            
    !         do j = 1, ndeg 
    !             print*, 'state'
    !             do k = 1, ndeg 
    !                 print* ,mat(k,j)
    !             end do 
    !         end do 
    !         ! print* ,sum(evals), 'sum(evals)'
    !         print* , ''    
    !         if(ndeg == 2) pause 
    !         mat = (0.d0, 0.d0)  
    !     end do
    !     call saveZ("A", 11, dir, name, pars, nbonds, ndeg, 1, bcurrent, 1) 
    !     if(ndeg == 2) pause 

    !     deallocate(psiprime)
        
        
    !     return

! end subroutine jij_irrep_deg_rect

subroutine jij_irrep_deg(ndeg, t, dim, sites, Lx, Ly, nHel, symmetrize, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, psi, nnnVec, norm) !, bcurrent, current

    implicit none
    integer, intent(in) :: ndeg, nbonds, sites, Lx, Ly, nHel, symmetrize
    integer, intent(in) :: bsites(5,nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    integer(kind=8), intent(in) :: dim, basis(dim)
    double precision, intent(in) :: id, mir(6), rot(5), t, psi(dim, ndeg), nnnVec(2, 3), norm(dim) 
    integer(kind=8) :: loc = 0, newst = 0, rep = 0, state = 0
    integer :: x, y, info, i, j, k, l1, l2, nd1, nd2
    integer :: sign = 1, signrep = 1, signstate = 1, signDir = 1, order = 0 
    integer :: site1, site2, phase, parity
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: evals(:), bcurrent(:,:)
    double complex, allocatable :: psiprime(:), bcurrmat(:,:,:), mat(:,:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

   
    if(allocated(bcurrmat)) deallocate(bcurrmat)
    if(allocated(bcurrent)) deallocate(bcurrent)
    if(allocated(psiprime)) deallocate(psiprime)
    
    allocate(bcurrmat(nbonds, ndeg, ndeg))
    allocate(bcurrent(nbonds, ndeg))
    allocate(psiprime(dim))
  
    bcurrmat = (0.d0,0.d0) 
    bcurrent = 0.d0 

    if(symmetrize == 1) then 
        order = 12
    else 
        order = 1
    end if 

    do nd1 = 1, ndeg 
        do nd2 = 1, ndeg 
            do i = 1, nbonds !Loop over bonds    
                psiprime = 0.0d0  
                site1    = bsites(1,i)
                site2    = bsites(2,i)
                phase    = bsites(3,i)   
                vec      = vecs(1:2,bsites(5,i))
                dist     = 1.d0

                do j = 1, dim !Loop over basis states 
                    do k = 1, order !Loop over point groups operators 
                        if(k == 1) then !Identity
                            state = basis(j)
                        else if(k <= 7) then !Reflections
                            call reflect(basis(j), basis(j), sites, refl(k - 1, 1:sites), signstate, info, state)
                            signstate = signstate * mir(k - 1)
                        else if(8 <= k) then !Rotations
                            call c6n(basis(j), basis(j), sites, k - 7, c6, signstate, info, state)
                            signstate = signstate * rot(k - 7)
                        end if 
                        
                        do x = 1, Lx !X-translations
                            call xtranslate(state, nHel, sites, Lx, xtransl, signstate, state)

                            do y = 1, Ly !Y-translations
                                call ytranslate(state, nHel, 0, sites, Lx, ytransl, signstate, state)

                                if(btest(state, site1-1) .and. .not.(btest(state, site2-1))) then !Forward forward hopping: 1 -> 2       
                                    parity  = popcnt(ibits(state, site1, sites)) + popcnt(ibits(ibclr(state, site1-1), site2, sites))                     
                                    newst   = ibclr(ibset(state, site2-1), site1-1) !Create on site 2, annihilate on site 1                             
                                    signDir = 1
                                else if(btest(state, site2-1) .and. .not.(btest(state, site1-1))) then !Backwards hopping: 2 -> 1
                                    parity  = popcnt(ibits(state, site2, sites)) + popcnt(ibits(ibclr(state, site2-1), site1, sites)) 
                                    newst   = ibclr(ibset(state, site1-1), site2-1) !Create on site 1, annihilate on site 2
                                    signDir = -1
                                else 
                                    cycle 
                                end if 

                                call representative_irrep(newst, sites, nHel, 1, Lx, Ly, symmetrize, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                                call findstate(dim, rep, basis, loc) 
                                
                                if(loc <= 0) cycle   
                                                                
                                sign = signstate * signrep * signDir * (-1)**parity 
                                
                                bcurrmat(i, nd2, nd1) = bcurrmat(i, nd2, nd1) + sign * ii * sqrt(norm(loc)/norm(j)) * psi(j, nd1) * psi(loc, nd2)
                                ! bcurrmat(i, nd2, nd1) = bcurrmat(i, nd2, nd1) + dble(phase) * dist * sign * t * ii * sqrt(norm(loc)/norm(j)) * psi(j, nd1) * psi(loc, nd2)
                                ! psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(j, nd1)
                                ! psiprime(j) = psiprime(j) + sign * (-(t)**2) * sqrt(norm(j)/norm(loc)) * psi(loc, nd1)
                                
                            end do !y
                        end do !x 
                    end do !G 
                    
                end do !j 
                ! bcurrmat(i,nd2,nd1) = bcurrmat(i,nd2,nd1) + dble(phase) * dist * dot_product(psi(1:dim, nd2), psiprime)
            end do !Loop over bonds    
        end do 
    end do !Degenerate states

    do k = 1, nbonds
        ! print* , k, 'bond'
        ! do i = 1, ndeg
        !     do j = 1, ndeg               
        !         print* ,bcurrmat(k, i, j), i, j, '(J_k)_ij, i, j'
        !     end do 
        ! end do 
        ! print*, '---------------'
    end do 
    ! pause 
    if(allocated(mat)) deallocate(mat)
    allocate(mat(ndeg, ndeg))
    ! bcurrmat = bcurrmat / ndeg
    do i = 1, nbonds          
        mat = bcurrmat(i,1:ndeg,1:ndeg)
        call exactdiag_c(.True., ndeg, mat, evals)
        bcurrent(i, 1:ndeg) = evals
        
        ! do j = 1, ndeg 
        !     print* ,evals(j), j, i, 'evals(state), state, bond'
        ! end do 
        ! do j = 1, ndeg 
        !     print*, 'state'
        !     do k = 1, ndeg 
        !         print* ,mat(k,j)
        !     end do 
        ! end do 
        ! print* ,sum(evals), 'sum(evals)'
        ! print* , ''    
        ! if(ndeg == 2) pause
         
        mat = (0.d0, 0.d0)  
    end do 
    ! if(ndeg == 2) pause 

    deallocate(psiprime)
    
    
    return

end subroutine jij_irrep_deg






subroutine ccf(t, dim, sites, basis, bsites, phases, nbonds, refbond, psi, bondcurcurrent, current)

    implicit none
    integer, intent(in) :: nbonds, refbond, sites  
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(2,nbonds)
    integer, intent(in) :: phases(nbonds)
    double precision, intent(in) :: t 
    double precision, intent(in) :: psi(dim)

    double precision, intent(out) :: current 
    ! double precision, intent(out) :: bondcurrent(nbonds)
    double precision, intent(out) :: bondcurcurrent(nbonds)
    
    integer(kind=8) :: loc = 0, newst = 0
    integer :: i = 0, j = 0, k = 0 
    integer :: forward = 0, backwards = 0
    integer :: signf = 0, signb = 0, cntr = 0 
    integer :: site1 = 0, site2 = 0, phase = 0   

    ! double precision :: bondcurcurrent(nbonds)
    double precision :: current_temp
    double precision, allocatable :: bondcurrent(:)
    double precision, allocatable :: psiprime(:)
    double precision, allocatable :: psiprimeij(:)
    double precision, allocatable :: psiprimeprime(:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)

    !j_\alpha = -i \sum_{i,\sigma} t (c^\dagger_{i\sigma}c_{i+\delta_\alpha,\sigma} -c^\dagger_{i+\delta_\alpha,\sigma} c_{i\sigma})
    if( allocated(bondcurrent) ) deallocate(bondcurrent)
    allocate(bondcurrent(nbonds))
    current        = 0 
    bondcurrent    = 0 
    bondcurcurrent = 0 

    if( allocated(psiprime) ) deallocate(psiprime)
    if( allocated(psiprimeij) ) deallocate(psiprimeij)
    if( allocated(psiprimeprime) ) deallocate(psiprimeprime)
    allocate( psiprime(dim) )
    allocate( psiprimeij(dim) )
    allocate( psiprimeprime(dim) )
    psiprime      = 0 
    psiprimeij    = 0 
    psiprimeprime = 0 

    i = refbond 
    do j = 1, dim
        ! do k = 1, 2
        !     if(k == 1) then 
        ! if(bsites(1,i) < bsites(2,i)) then 
                site1 = bsites(1,i)
                site2 = bsites(2,i)                
        ! else 
        !         site2 = bsites(1,i)
        !         site1 = bsites(2,i)
        ! end if 
            ! else if(k == 2) then 
            !     site1 = bsites(2,i)
            !     site2 = bsites(1,i)            
            ! end if 
            current_temp   = 0 
            forward        = 0
            backwards      = 0 
            if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Hopping from 1->2
                forward = forward + 1 
                newst   = ibclr( ibset( basis(j), site2 - 1 ), site1 - 1 ) !Create site 2, annihilate site 1
                call findstate (dim, newst, basis, loc)     
                if( basis(loc) .ne. newst) error stop "Subroutine slcurrent: Wrong location of scattered basis state."
                signf = (-1)**(popcnt( ibits( basis(j), site1, sites - site1 ) ) + &
                                popcnt( ibits( ibclr(basis(j), site1 - 1), site2, sites - site2 ) ) ) !Count occupied sites over which operators have to hop. -1 for already occupied site. -1 for hopping over already annihilated site.
            else if( btest( basis(j), site2 - 1 ) .and. .not.( btest( basis(j), site1 - 1 ) ) ) then !Hopping from 2->1 
                backwards = backwards + 1 
                newst = ibclr( ibset( basis(j), site1 - 1 ), site2 - 1 ) 
                call findstate (dim, newst, basis, loc)
                if( basis(loc) .ne. newst) error stop "Subroutine slcurrent: Wrong location of scattered basis state."
                signb = (-1)**(popcnt( ibits( ibclr(basis(j), site2 - 1), site1, sites - site1 ) ) + &
                                popcnt( ibits( basis(j), site2, sites - site2 ) ) ) 
            else 
                cycle    
            end if 
            current_temp  = -t * ( signf * forward - signb * backwards )
            psiprime(loc) = psiprime(loc) + current_temp * psi(j)
            ! psiprime(loc) = psiprime(loc) + bsites(3,i) * current_temp * psi(j)
        ! end do 
    end do             

    do i = 1, nbonds 
    ! if(bsites(1,i) < bsites(2,i)) then 
            site1 = bsites(1,i)
            site2 = bsites(2,i)                     
    ! else 
    !         site2 = bsites(1,i)
    !         site1 = bsites(2,i)             
    ! end if 
            phase = phases(i)   
        if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
            .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
        psiprimeprime = 0 
        do j = 1, dim            
            current_temp   = 0 
            forward        = 0
            backwards      = 0 
            if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Hopping from 1->2                 
                forward = forward + 1 
                newst   = ibclr( ibset( basis(j), site2 - 1 ), site1 - 1 ) 
                call findstate (dim, newst, basis, loc)
                if( basis(loc) .ne. newst) error stop "Subroutine slcurrent: Wrong location of scattered basis state."
                signf = (-1)**(popcnt(ibits(basis(j), site1, sites - site1 ) ) + &
                            popcnt(ibits(ibclr(basis(j), site1 - 1), site2, sites - site2 ) ) ) !Count occupied sites over which operators have to hop. -1 for already occupied site. -1 for hopping over already annihilated site.
            else if(btest(basis(j), site2  - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then                 
                backwards = backwards + 1 
                newst     = ibclr( ibset( basis(j), site1 - 1 ), site2 - 1 ) 
                call findstate (dim, newst, basis, loc)
                if( basis(loc) .ne. newst) error stop "Subroutine slcurrent: Wrong location of scattered basis state."
                signb = (-1)**(popcnt(ibits(ibclr(basis(j), site2 - 1), site1, sites - site1 ) ) + &
                            popcnt(ibits(basis(j), site2, sites - site2 ) ) ) 
            else
                cycle 
            end if 
            current_temp       = - t * ( signf * forward - signb * backwards )
            psiprimeij(loc)    = psiprimeij(loc) + current_temp * psi(j)
            psiprimeprime(loc) = psiprimeprime(loc) + current_temp * psiprime(j)
        end do           
        bondcurrent(i)    = 0 
        bondcurrent(i)    = dot_product(psi, psiprimeij)
        bondcurcurrent(i) = 0 
        bondcurcurrent(i) = dot_product(psi, psiprimeprime)
        current           = current + phase * bondcurcurrent(i)
        ! current = current +  bondcurcurrent(i)
    end do

    deallocate(psiprime)
    deallocate(psiprimeij)
    deallocate(psiprimeprime)

    return

end subroutine ccf 




!---------------------------------!
!            Thread units         !
!---------------------------------!

subroutine threadunits(ndv, ndv2, units)

    implicit none
    integer, intent(in) :: ndv, ndv2
    integer, allocatable, intent(out) :: units(:, :)
    integer :: i, j, h
    integer :: n_threads_tot

    n_threads_tot = 1

    if (ndv > 0) n_threads_tot = n_threads_tot + ndv
    if (ndv2 > 0) n_threads_tot = n_threads_tot + ndv2

    if(allocated(units)) deallocate(units)
    allocate(units(n_threads_tot,n_threads_tot))
    h = 0
    do i = 1, n_threads_tot
        do j = 1, n_threads_tot
            h = h + 1
            units(i,j) = h
        end do
    end do

    return 

end subroutine

!-------------------------------!
!            Step units         !
!-------------------------------!

subroutine stepunits(ndu, ndv, ndv2, units)
    !Creates array with units for parallel I/O. Each entry is a different unit number.
    implicit none
    integer, intent(in) :: ndu, ndv, ndv2
    integer, allocatable, intent(out) ::  units(:,:,:,:)
    integer :: i, j, k, l, q

    if (allocated(units)) deallocate(units)
    allocate(units(ndv2,ndv,ndu,2))
    q = 0
    do i = 1, ndv2
        do j = 1, ndv
            do k = 1, ndu
                do l = 1, 2
                    q = q + 1
                    units(i,j,k,l) = q
                end do
            end do
        end do
    end do

end subroutine

!-------------------------------------------!
!            Date and Time                  !
!-------------------------------------------!

subroutine datetime(dir, startend, params)
    implicit none

    integer, intent(in) :: startend
    character(len=*), intent(in) :: dir
    character(len=*), intent(in), optional :: params
    character(8), save  :: datei, datef
    character(10), save :: timei, timef
    integer,dimension(8), save :: valuesi, valuesf
    real, save :: start, finish
    character :: file_name*400


    if(startend == 0) then
        call cpu_time(start)
        call date_and_time(date=datei,time=timei,values=valuesi)
        write(*,"('Calculation started at',x,a,' h',x,a,' min',x,a,' sec')") timei(1:2), timei(3:4), timei(5:6)
        print*, ''
        print*, 'Start date: ',datei(7:8), '.',datei(5:6), '.',datei(1:4)
        print*, ''
    else
        call cpu_time(finish)
        print*, 'Start = ', start
        print*, 'Finish = ', finish
        if (finish - start < 60) then
            write(*,"(' Elapsed CPU time = ',f12.3,' seconds.')") finish-start
            print*, ''
        else if (finish - start < 3600) then
            write(*,"(' Elapsed CPU time = ',f12.3,' minutes.')") (finish-start)/60
            print*, ''
        else
            write(*,"(' Elapsed CPU time = ',f12.3,' hours.')") (finish-start)/3600
            print*, ''
        end if
        call date_and_time(date=datef,time=timef,values=valuesf)
        write(*,"(' Calculation started at',x,a,'h',x,a,'min',x,a,'sec')") timei(1:2), timei(3:4), timei(5:6)
        print*, ''
        write(*,"(' Calculation ended at',x,a,'h',x,a,'min',x,a,'sec')") timef(1:2), timef(3:4), timef(5:6)
        print*, ''
        print*, 'Start date: ',datei(7:8), '.',datei(5:6), '.',datei(1:4)
        print*, ''
        print*, 'End date: ',datef(7:8), '.',datef(5:6), '.',datef(1:4)

        file_name=dir//"times_"//params
        file_name=trim_name(file_name)
        open(77,file = file_name)
        write(77,"(' Elapsed CPU time = ',f20.10,' seconds.')") finish-start
        write(77,"(' Calculation started at',x,a,'h',x,a,'min',x,a,'sec')") timei(1:2), timei(3:4), timei(5:6)
        write(77,"(' Start date:',x,a,'.',x,a,'.',x,a)") ,datei(7:8), datei(5:6), datei(1:4)
        write(77,"(' Calculation ended at',x,a,'h',x,a,'min',x,a,'sec')") timef(1:2), timef(3:4), timef(5:6)
        write(77,"(' End date:',x,a,'.',x,a,'.',x,a)") ,datef(7:8), datef(5:6), datef(1:4)
        close(77)
    end if
end subroutine datetime

!---------------------------------------------------!
!            Number of discretization steps         !
!---------------------------------------------------!

subroutine nsteps(min, max, delta, steps)

    implicit none
    double precision, intent(in) :: min, max, delta
    integer, intent(out) :: steps

    if(delta <= 1.d0) then
        steps = int(abs(max-min)/delta + delta/2) + 1
    else
        steps = int(abs(max-min)/delta) + 1
    end if


end subroutine nsteps

!------------------------------------------!
!            Check parallelization         !
!------------------------------------------!

subroutine parallelcheck()

    implicit none
    integer :: num_threads_loc
    integer :: thread_num_loc
    integer :: max_threads_loc

    num_threads_loc = 0
    thread_num_loc = 0
    max_threads_loc = 0
    !Critical block is a lock, forcing the instructions to be run in serial
    print*, ''
    !$omp parallel
        !$omp critical
            !$ thread_num_loc = omp_get_thread_num()
            print('(1x,100(a,i0))'), 'Thread number ', thread_num_loc, ' is online.'
            !$ num_threads_loc = omp_get_num_threads()
            if (thread_num_loc == 0) then
                print*, ''
                print('(1x,100(a,i0))'), 'Available number of threads = ', num_threads_loc
                !$ max_threads_loc = omp_get_max_threads()
                print('(1x,100(a,i0))'), 'Maximum number of threads = ', max_threads_loc
                print*, ''
            end if
        !$omp end critical
    !$omp end parallel
    print*, ''
    return 

end subroutine

!-------------------------------------!
!         Disorder average            !
!-------------------------------------!

!subroutine disav ()
    !
    !    implicit none
    !
    !    do i = 1, nDis
    !        av_evals =
    !    end do
    !
!end subroutine disav


!-------------------------------------!
!            Error messages           !
!-------------------------------------!

!Error messages for diagonalization routines
subroutine dsaupderrormessage(dsaupdinfo)
    implicit none
    integer :: dsaupdinfo

    if (dsaupdinfo .EQ. 0) THEN
    PRINT *, 'Normal Exit in dsaupd: info = 0.'
    elseif (dsaupdinfo .EQ. -1) THEN
    PRINT *, 'Error in dsaupd: info = -1.'
    PRINT *, 'N must be positive.'
    elseif (dsaupdinfo .EQ. -2) THEN
    PRINT *, 'Error in dsaupd: info = -2.'
    PRINT *, 'NEV must be positive.'
    elseif (dsaupdinfo .EQ. -3) THEN
    PRINT *, 'Error in dsaupd: info = -3.'
    PRINT *, 'NCV must be between NEV and N. '
    elseif (dsaupdinfo .EQ. -4) THEN
    PRINT *, 'Error in dsaupd: info = -4'
    PRINT *, 'The maximum number of Arnoldi update iterations allowed must be greater than zero.'
    elseif (dsaupdinfo .EQ. -5) THEN
    PRINT *, 'Error in dsaupd: info = -5'
    PRINT *, 'WHICH must be LM, SM, LA, SA, or BE. info = -5.'
    elseif (dsaupdinfo .EQ. -6) THEN
    PRINT *, 'Error in dsauupd: info = -6. '
    PRINT *, 'BMAT must be I or G. '
    elseif (dsaupdinfo .EQ. -7) THEN
    PRINT *, 'Error in dsaupd: info = -7.'
    PRINT *, 'Length of private work work WORKL array isnt sufficient.'
    elseif (dsaupdinfo .EQ. -8) THEN
    PRINT *, 'Error in dsaupd: info = -8.'
    PRINT *, 'Error in return from trid. eval calc. Error info from LAPACK dsteqr. info =-8'
    elseif (dsaupdinfo .EQ. -9) THEN
    PRINT *, 'Error in dsaupd: info = -9.'
    PRINT *, 'Starting vector is 0.'
    elseif (dsaupdinfo .EQ. -10) THEN
    PRINT *, 'Error in dsaupd: info = -10. '
    PRINT *, 'IPARAM(7) must be 1,2,3,4, or 5.'
    elseif (dsaupdinfo .EQ. -11) THEN
    PRINT *, 'Error in dsaupd: info = -11.'
    PRINT *, 'IPARAM(7)=1 and BMAT=G are incompatible.'
    elseif (dsaupdinfo .EQ. -12) THEN
    PRINT *, 'Error in dsaupd: info = -12'
    PRINT *, 'NEV and WHICH=BE are incompatible.'
    elseif (dsaupdinfo .EQ. -13) THEN
    PRINT *, 'Error in dsaupd: info = -13.'
    PRINT *, 'DSAUPD did find any eigenvalues to sufficient accuracy.'
    elseif (dsaupdinfo .EQ. -9999) THEN
    PRINT *, 'Error in dsaupd: info = -9999'
    PRINT *, 'Could not build an Arnoldi factorization. IPARAM(5) returns the size of the current Arnoldi factorization. &
    The user is advised to check that enough workspace and array storage     has been allocated. '
    elseif (dsaupdinfo .EQ. 1) THEN
    PRINT *, 'Error in dsaupd: info = 1'
    PRINT *, 'Maximum number of iterations taken. All possible eigenvalues of OP has been found. '
    PRINT *, 'IPARAM(5) returns the number of wanted converged Ritz values.'
    elseif (dsaupdinfo .EQ. 3) THEN
    PRINT *, 'Error in dsaupd: info =3'
    PRINT *, 'No shifts could be applied during a cycle of the Implicitly restarted Arnoldi iteration.'
    PRINT *, 'One possibility is to increase the size of NCV relative to NEV.'
    else
    PRINT *, 'Unknown error.  info =', dsaupdinfo
    END IF
end subroutine dsaupderrormessage

subroutine dseupderrormessage(dseupdinfo)
    implicit none
    integer :: dseupdinfo

    if (dseupdinfo .EQ. 0) THEN
    PRINT *, 'Normal Exit in dseupd: info = 0.'
    elseif (dseupdinfo .EQ. -1) THEN
    PRINT *, 'Error in deseupd: N must be positive. info =-1.'
    elseif (dseupdinfo .EQ. -2) THEN
    PRINT *, 'Error in deseupd: NEV must be positive. info = -2.'
    elseif (dseupdinfo .EQ. -3) THEN
    PRINT *, 'Error in deseupd: NCV must be between NEV and N. info = -3.'
    elseif (dseupdinfo .EQ. -5) THEN
    PRINT *, 'Error in deseupd: WHICH must be LM, SM, LA, SA, or BE info = -5.'
    elseif (dseupdinfo .EQ. -6) THEN
    PRINT *, 'Error in deseupd: BMAT must be I or G. info = -6.'
    elseif (dseupdinfo .EQ. -7) THEN
    PRINT *, 'Error in deseupd: N Length of private work work WORKL array isnt sufficient. info = -7.'
    elseif (dseupdinfo .EQ. -8) THEN
    PRINT *, 'Error in deseupd: Error in return from trid. eval calc. Error info from LAPACK dsteqr. info = -8.'
    elseif (dseupdinfo .EQ. -9) THEN
    PRINT *, 'Error in deseupd: Starting vector is 0. info = -9.'
    elseif (dseupdinfo .EQ. -10) THEN
    PRINT *, 'Error in deseupd: IPARAM(7) must be 1,2,3,4, or 5. info = -10.'
    elseif (dseupdinfo .EQ. -11) THEN
    PRINT *, 'Error in deseupd: IPARAM(7)=1 and BMAT=G are incompatible. info = -11.'
    elseif (dseupdinfo .EQ. -12) THEN
    PRINT *, 'Error in deseupd: NEV and WHICH=BE are incompatible. info = -12.'
    elseif (dseupdinfo .EQ. -14) THEN
    PRINT *, 'Error in deseupd: DSAUPD did find any eigenvalues to sufficient accuracy. info = -14.'
    elseif (dseupdinfo .EQ. -15) THEN
    PRINT *, 'Error in deseupd: HOWMANY must one A or S if RVEC=1. info = -15.'
    elseif (dseupdinfo .EQ. -16) THEN
    PRINT *, 'Error in deseupd: HOWMANY =S not yet implemented. info = -16.'
    elseif (dseupdinfo .EQ. -17) THEN
    PRINT *, 'Error in deseupd: info =-17.'
    PRINT *, 'DSEUPD got a different count of the number of converged Ritz values than DSAUPD.'
    PRINT *, 'User likely made an error in passing data DSAUPD -> DSEUPD. info = -17.'
    else
    PRINT *, 'Unknown error.  info =', dseupdinfo
    END IF

end subroutine dseupderrormessage

subroutine printparams(nev, ndeg, nHel, tilt, k1, k2)

    use variables_tmi
    implicit none
    integer, intent(in) :: nev, ndeg, nHel, tilt, k1, k2
    print*, ''
    print*, '!-------------------------------------------------------------!'
    print*, '!                                                             !'
    print*, '!                          Parameters                         !'
    print*, '!                                                             !'
    print*, '!-------------------------------------------------------------!'
    print*, ''
    print('(1x, a,i0)'), 'L = ', sites
    print('(1x, a,i0)'), 'N = ', particles
    if(tilted == 0) then 
    print('(1x, a,i0)'), 'Unit cells in x direction = ', ucx
    print('(1x, a,i0)'), 'Unit cells in y direction = ', ucy
    else if(tilted == 1) then 
        print('(1x, a, a)'), 'Cluster = ', cluster 
        print('(1x, a, i0)'),'Helices = ', nHel 
        print('(1x, a, i0)'),'Tilt    = ', tilt 
    end if 
    print*, ''
    print('(1x, a,i0)'), 'DIM = ', dim
    print('(1x, a,i0)'), 'NEV = ', nev

    if(ndeg > 0) then 
        if(degeneracy == 1) then 
            print*, ''
            ! print('(1x, a,f20.16)'), 'Degeneracy threshold = ', deg
            print('(1x, a,i0)'), 'NDEG = ', ndeg
            print*, ''
        else if(degeneracy == 2) then 
            print('(1x, a,f20.16)'), 'QNDEG = ', ndeg
            print*, ''
        end if 
    end if
    if(ti == 1) then 
        print('(1x, a,i0)'), 'K1   = ', k1
        print('(1x, a,i0)'), 'K2   = ', k2
        print('(1x, a,i0)'), 'Ktot = ', k1 + ucx * k2
        print*, ''    
    end if 
    if(symmetrize == 1) print('(1x, a,a)'),'IRREP = ', irrep 
    if(symmetrize == 1) print*, ''
    print*, '---------------------------------------------------------------'
    print*, ''    
    print('(1x, a,i0)'), 'OMP = ', othrds
    print('(1x, a,i0)'), 'MKL = ', mthrds
    if(exact == 1) then 
        print('(1x, a)'), 'EXACT '
        print*, ''
    end if 
    if(arpack == 1) then 
        print('(1x, a)'), 'ARPACK '
        print*, ''
    end if 
    if(mkl == 1) then 
        print('(1x, a)'), 'MKL '
        print*, ''
    end if 
    if(feast == 1) then 
        print('(1x, a)'), 'FEAST '
        print*, ''
    end if 
    if(curr == 1) then 
        print('(1x, a)'), 'CURRENTS '
        print*, ''
    end if 

    print*, '---------------------------------------------------------------'
    print*, ''
    print('(1x, a)'), 'DISORDER'
    print('(1x, a,f6.4)'), 'W    = ', dis
    print('(1x, a,i0)'),   'Ndis = ', nDis
    print*, ''
    print*, '---------------------------------------------------------------'
    print*, ''
    print('(1x, a,f6.4)'), 'dv   = ',dv
    print('(1x, a,f6.4)'), 'Vmin = ', vmin
    print('(1x, a,f6.4)'), 'Vmax = ', vmax
    print*, ''
    print*, '---------------------------------------------------------------'
    print*, ''
    print('(1x, a,f6.4)'), 'dv2   = ',dv2
    print('(1x, a,f6.4)'), 'V2min = ', v2min
    print('(1x, a,f6.4)'), 'V2max = ', v2max
    print*, ''
    print*, '---------------------------------------------------------------'    
    print*, ''
    if(mass > 0.d0) print('(1x, a,f6.4)'), 'Mass = ', mass
    if(bc == 'p') then
        print*, ''
        print('(1x, a)'), 'PBC'
        print*, ''
    elseif (bc == 'o') then
        print*, ''
        print('(1x, a)'), 'OBC'
        print*, ''
    elseif (bc == 't') then
        print*, ''
        print('(1x, a)'), 'TBC'
        print*, ''
    end if
    print*, '!-------------------------------------------------------------!'
    print*, '!-------------------------------------------------------------!'
    print*, ''

end subroutine printparams




















!--------------------------------------------!
!           Older current routines           !
!--------------------------------------------!



subroutine current_c(gencluster, nHel, tilt, k1, k2, t, dim, sites, Lx, Ly, basis, bsites, xy, xtransl, ytransl, nbonds, refbond, psi, nnnVec, norm, bondcurrent, current)

    implicit none
    integer, intent(in) :: gencluster, nHel, tilt, k1, k2, nbonds, refbond, sites, Lx, Ly
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    integer, intent(in) :: xy(4,nbonds)
    integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
    double precision, intent(in) :: t 
    double precision, intent(in) :: nnnVec(2, 3)
    double complex, intent(in) :: psi(dim)
    double precision, intent(in) :: norm(dim)
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bondcurrent(:)
    
    integer(kind=8) :: loc = 0, newst = 0, rep = 0, cntr = 0 
    integer :: i = 0, j = 0, l1 = 0, l2 = 0, l = 0, l11 = 0, l22 = 0 
    integer :: x = 0, y = 0
    integer :: sign = 0, signt = 0 
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double complex, allocatable :: psiprime(:)
    character*250 :: name, name2 

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bondcurrent)) deallocate(bondcurrent)
    allocate(bondcurrent(nbonds))
    if(allocated(psiprime)) deallocate(psiprime)
    allocate(psiprime(dim))

    psiprime = 0.0d0
    refvec   = vecs(1:2,bsites(5,refbond))
    refsite1 = bsites(1,refbond)
    refsite2 = bsites(2,refbond)       
    current = 0.d0
    bondcurrent = 0.d0 
    
    do i = 1, nbonds !Loop over bonds    
        bondcurrent(i) = 0.0d0
        psiprime = 0.0d0  
        site1 = bsites(1,i)
        site2 = bsites(2,i)
        refsite1 = bsites(1,refbond)
        refsite2 = bsites(2,refbond)    
        phase = bsites(3,i)   
        vec   = vecs(1:2,bsites(5,i))
        dist  = dot_product(vec, refvec)

        do x = 0, Lx - 1
            do y = 0, Ly - 1
                if (site1 == refsite1 .or. site2 == refsite2 &
                    .or. site1 == refsite2 .or. site2 == refsite1) go to 110 !Exclude bonds which share a site with the reference bond

                do j = 1, dim !Loop over basis states 
                    if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                        newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                        if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                            ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                            newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                            if(gencluster == 0) then 
                                call representative2(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                l11 = Ly
                                l22 = Lx
                            else if(gencluster == 1) then 
                                call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                l11 = Lx
                                l22 = Ly
                            end if
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(j) = psiprime(j) + sign * signt * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(loc) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                        else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                            ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                            newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                            if(gencluster == 0) then 
                                call representative2(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                l11 = Ly
                                l22 = Lx
                            else if(gencluster == 1) then 
                                call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                l11 = Lx
                                l22 = Ly
                            end if
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(j) = psiprime(j) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(loc) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                        else 
                            cycle 
                        end if     
                    else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                        newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                        if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                            ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                            newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                            if(gencluster == 0) then 
                                call representative2(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                l11 = Ly
                                l22 = Lx
                            else if(gencluster == 1) then 
                                call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                l11 = Lx
                                l22 = Ly
                            end if
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            
                            parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(j) = psiprime(j) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(loc) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                        else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                            ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                            newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                            if(gencluster == 0) then 
                                call representative2(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                l11 = Ly
                                l22 = Lx
                            else if(gencluster == 1) then 
                                call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                l11 = Lx
                                l22 = Ly
                            end if
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            
                            parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(j) = psiprime(j) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(loc) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                        else 
                            cycle 
                        end if     
                    else 
                        cycle 
                    end if 
                end do !Loop over basis states
                
                110 continue 
                
                !a1 translations 
                cntr = 0
                do l = 1, nbonds
                    if(bsites(1,l) == site1 .and. bsites(5, l) == 1) then 
                        site1 = bsites(2,l)
                        cntr = cntr + 1
                    else if(bsites(1,l) == site2 .and. bsites(5, l) == 1) then
                        site2 = bsites(2,l)
                        cntr = cntr + 1
                    else if(bsites(1,l) == refsite1 .and. bsites(5, l) == 1) then
                        refsite1 = bsites(2,l)
                        cntr = cntr + 1
                    else if(bsites(1,l) == refsite2 .and. bsites(5, l) == 1) then
                        refsite2 = bsites(2,l)
                        cntr = cntr + 1
                    else if(cntr == 4) then 
                        exit 
                    end if
                end do 

            end do 

            !a2 translations 
            cntr = 0
            do l = 1, nbonds
                if(bsites(1,l) == site1 .and. bsites(5, l) == 2) then 
                    site1 = bsites(2,l)
                    cntr = cntr + 1
                else if(bsites(1,l) == site2 .and. bsites(5, l) == 2) then
                    site2 = bsites(2,l)
                    cntr = cntr + 1
                else if(bsites(1,l) == refsite1 .and. bsites(5, l) == 2) then
                    refsite1 = bsites(2,l)
                    cntr = cntr + 1
                else if(bsites(1,l) == refsite2 .and. bsites(5, l) == 2) then
                    refsite2 = bsites(2,l)
                    cntr = cntr + 1
                else if(cntr == 4) then 
                    exit 
                end if
            end do 

        end do 
        bondcurrent(i) = bondcurrent(i) + dist * dot_product(psi, psiprime) / (Lx*Ly)
        current = current + dble(phase) * bondcurrent(i)
    end do !Loop over bonds 
    

    deallocate(psiprime)

    return

end subroutine current_c

subroutine currmat_c(gencluster, nHel, tilt, dir, unit, filename, sl, ndeg, k1, k2, t, dim, sites, Lx, Ly, basis, bsites, xy, xtransl, ytransl, nbonds, refbond, psi, nnnVec, norm)

    implicit none
    integer, intent(in) :: gencluster, nHel, tilt, unit, ndeg, k1, k2, nbonds, refbond, sites, Lx, Ly  
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    integer, intent(in) :: xy(4,nbonds)
    integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
    double precision, intent(in) :: t 
    double precision, intent(in) :: nnnVec(2, 3)
    double complex, intent(in) :: psi(dim, ndeg)
    double precision, intent(in) :: norm(dim)
    character(len=*), intent(in) :: dir, sl, filename
    
    integer(kind=8) :: loc = 0, newst = 0, rep = 0
    integer :: i = 0, j = 0, k = 0, l = 0, l1 = 0, l2 = 0
    integer :: sign = 0, signt = 0 
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0, current= 0.d0 
    double precision, dimension(2)   :: vec, refvec, coord1, coord2 
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: currents(:)
    double precision, allocatable :: currents2(:)
    double precision, allocatable :: evals(:)
    double precision, allocatable :: bondcurrents(:,:), bondcurrent(:)
    double precision, allocatable :: currentmat(:,:)
    double precision, allocatable :: bondcurrentmat(:,:,:)
    double complex, allocatable :: psiprime(:), gs(:)
    character*300 :: name 

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(psiprime)) deallocate(psiprime)
    if(allocated(bondcurrentmat)) deallocate(bondcurrentmat)
    if(allocated(currentmat)) deallocate(currentmat)
    if(allocated(currents2)) deallocate(currents2)
    if(allocated(bondcurrents)) deallocate(bondcurrents)
    allocate(psiprime(dim))
    allocate(bondcurrentmat(nbonds, ndeg, ndeg))
    allocate(currentmat(ndeg, ndeg))
    allocate(currents2(ndeg))
    allocate(bondcurrents(nbonds, ndeg))

    psiprime = 0.d0
    refvec   = vecs(1:2,bsites(5,refbond))
    coord1(1:2) = int(0.5*xy(3,refbond))*vecs(1:2,1) + xy(4,refbond)*vecs(1:2,2) - (int(0.5*xy(1,refbond))*vecs(1:2,1) + xy(2,refbond)*vecs(1:2,2))
    refsite1          = bsites(1,refbond)
    refsite2          = bsites(2,refbond)      

    ! do i = 1, ndeg
        !     do j = 1, ndeg
        !         if ( abs(dot_product(psi(1:dim,i),psi(1:dim,j))) > 10.**(-10) .and. abs(dot_product(psi(1:dim,i),psi(1:dim,j))) < 0.999999999d0 ) then 
        !             print*,'currmat_c: Groundstates ',i, 'and' ,j,' are not orthorgonal.', dot_product(psi(1:dim,i),psi(1:dim,j))               
        !             error stop 
        !         end if 
        !     end do 
        ! end do 
    
    currentmat = 0.d0     
    bondcurrentmat = 0.d0  
    do l = 1, ndeg 
        do k = 1, ndeg 
            do i = 1, nbonds                            
                site1 = bsites(1,i)
                site2 = bsites(2,i)
                if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
                    .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
                coord2(1:2) = int(0.5*xy(3,i))*vecs(1:2,1) + xy(4,i)*vecs(1:2,2) - (int(0.5*xy(1,i))*vecs(1:2,1) + xy(2,i)*vecs(1:2,2))
                ! coord2(1:2) = coord2(1:2)/norm2(coord2(1:2)) 
                phase = bsites(3,i)   
                vec   = vecs(1:2,bsites(5,i))
                dist = dot_product(vec, refvec)
                psiprime = 0 

                do j = 1, dim
                    if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                        newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                        if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                            ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                            newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                            call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)  
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(Lx)) + (dble(k2*l2)/dble(Ly))))
                        else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                            ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                            newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                            call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)  
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(Lx)) + (dble(k2*l2)/dble(Ly))))
                        end if     
                    else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                        newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                        if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                            ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                            newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                            call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)  
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(Lx)) + (dble(k2*l2)/dble(Ly))))
                        else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                            ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                            newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                            call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)  
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(Lx)) + (dble(k2*l2)/dble(Ly))))
                        end if     
                    end if 
                end do !Loop over basis states                  
                bondcurrentmat(i,k,l) = bondcurrentmat(i,k,l) + dble(phase) * dist * dot_product(psi(1:dim,k), psiprime)
                currentmat(k,l) = currentmat(k,l) + bondcurrentmat(i,k,l)
                ! bondcurrentmat(i,k,l) = bondcurrentmat(i,k,l) + dist * dot_product(psi(1:dim,k), psiprime)
                ! currentmat(k,l) = currentmat(k,l) + dble(phase) * bondcurrentmat(i,k,l)

            end do !Loop over bonds 
        end do !Loop over degenerate ground states
    end do !Loop over degenerate ground states
    deallocate(psiprime)

    call exactdiag(.True., ndeg, currentmat, currents)

    ! print*,'Current eigenvalues'
        ! do i = 1, ndeg 
        !     print*, currents(i)
        ! end do
        ! print*,'sum(currents)',sum(currents) 
    
        ! do i = 1, nbonds
        !     call exactdiag(.True., ndeg, bondcurrentmat(i, 1:ndeg, 1:ndeg), evals)
        !     bondcurrents(i, 1:ndeg) = evals
        !     if(allocated(evals)) deallocate(evals)
        ! end do 
        
        
        ! if(sl == "A") then 
        !     name = dir//"Acurrent_dens_mat_"//filename
        ! else if(sl == "B") then 
        !     name = dir//"Bcurrent_dens_mat_"//filename
        ! end if 
        ! name = trim_name(name)
        ! open(unit,file=name)
        ! do i = 1, ndeg 
        !     currents2(i) = sum(bondcurrents(1:nbonds, i))
        !     write(unit,*) currents2(i)
        !     ! print*,'Current of psi ',i, currents2(i)
        ! end do 
        ! close(unit)

    currentmat = currentmat / ndeg 


    call exactdiag(.True., ndeg, currentmat, currents)
    ! print*,'|===================|'
        ! print*,' Current eigenvalues'
        ! do i = 1, ndeg 
        !     print*, currents(i)
        ! end do
        ! print*,'|===================|'
        ! print*,' Total'
        ! print*,sum(currents) 
        ! currents = 0.d0 
        ! print*,'|===================|'
    if(allocated(gs)) deallocate(gs)
    allocate(gs(dim))

    do j = 1, ndeg 
        gs = 0.d0 
        do i = 1, ndeg 
            gs = gs + currentmat(i,j) * psi(1:dim, i)
        end do
        ! print*,'gs',gs
        call current_c(gencluster, nHel, tilt, k1, k2, t, dim, sites, Lx, Ly, basis, bsites, xy, xtransl, ytransl, nbonds, refbond, gs, nnnVec, norm, bondcurrent, current)
        bondcurrents(1:nbonds, j) = bondcurrent 
        currents(j) = current 
        ! print*,'|===================|'
            ! print*,' Current, eigenstate'
            ! print*, current, j
            ! print*,'|===================|'
        current = 0.d0 
    end do 
    ! print*,'|===================|'
        ! print*,' Sum'
        ! print*, sum(currents)
        ! print*,'|===================|'


    if(sl == "A") then 
        name = trim_name(dir//"Acurrent_matrix_"//filename)
    else if(sl == "B") then 
        name = trim_name(dir//"Bcurrent_matrix_"//filename)
    end if 
    name = trim_name(name)
    
    open(unit,file=name)
    do i = 1, ndeg 
        write(unit,*) currents(i)
    end do 
    close(unit)
    if(sl == "A") then 
        name = dir//"Abondcurrent_matrix_"//filename
    else if(sl == "B") then 
        name = dir//"Bbondcurrent_matrix_"//filename
    end if 
    name = trim_name(name)
    open(unit,file=name)
    do i = 1, ndeg 
        write(unit,*) bondcurrents(1:nbonds,i)
    end do 
    close(unit)

    ! bondcurrents = 0.d0 
        ! do i = 1, nbonds
        !     call exactdiag(.True., ndeg, bondcurrentmat(i, 1:ndeg, 1:ndeg), evals)
        !     bondcurrents(i, 1:ndeg) = evals
        !     if(allocated(evals)) deallocate(evals)
        ! end do 
        ! if(sl == "A") then 
        !     name = dir//"Acurrent_matrix2_"//filename
        ! else if(sl == "B") then 
        !     name = dir//"Bcurrent_matrix2_"//filename
        ! end if 
        ! name = trim_name(name)
        ! open(unit,file=name)
        ! do i = 1, ndeg 
        !     currents2(i) = sum(bondcurrents(1:nbonds, i))
        !     write(unit,*) currents2(i)
        !     print*,'Current of psi ',i, currents2(i)
        ! end do 
        ! close(unit)

    return

end subroutine currmat_c

subroutine currmat_c2(gencluster, nHel, tilt, dir, unit, filename, sl, ndeg, k1, k2, t, dim, sites, Lx, Ly, basis, bsites, xy, xtransl, ytransl, nbonds, refbond, psi, nnnVec, norm)

    implicit none
    integer, intent(in) :: gencluster, nHel, tilt, unit, ndeg, k1, k2, nbonds, refbond, sites, Lx, Ly  
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    integer, intent(in) :: xy(4,nbonds)
    integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
    double precision, intent(in) :: t 
    double precision, intent(in) :: nnnVec(2, 3)
    double complex, intent(in) :: psi(dim, ndeg)
    double precision, intent(in) :: norm(dim)
    character(len=*), intent(in) :: dir, sl, filename
    
    integer(kind=8) :: loc = 0, newst = 0, rep = 0
    integer :: i = 0, j = 0, k = 0, l = 0, l1 = 0, l2 = 0, x = 0, y = 0, cntr = 0, q = 0 
    integer :: sign = 0, signt = 0, l11 = 0, l22 = 0 
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0, current= 0.d0 
    double precision, dimension(2)   :: vec, refvec, coord1, coord2 
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: currents(:)
    double precision, allocatable :: currents2(:)
    double precision, allocatable :: bondcurrents(:,:), bondcurrent(:)
    double precision, allocatable :: currentmat(:,:)
    double precision, allocatable :: bondcurrentmat(:,:,:)
    double complex, allocatable :: psiprime(:), gs(:)
    character*300 :: name 

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(psiprime)) deallocate(psiprime)
    if(allocated(currents2)) deallocate(currents2)
    if(allocated(currentmat)) deallocate(currentmat)
    if(allocated(bondcurrents)) deallocate(bondcurrents)
    if(allocated(bondcurrentmat)) deallocate(bondcurrentmat)
    allocate(psiprime(dim))
    allocate(currents2(ndeg))
    allocate(currentmat(ndeg,ndeg))
    allocate(bondcurrents(nbonds,ndeg))
    allocate(bondcurrentmat(nbonds,ndeg,ndeg))

    psiprime = 0.d0
    refvec   = vecs(1:2,bsites(5,refbond))
    coord1(1:2) = int(0.5*xy(3,refbond))*vecs(1:2,1) + xy(4,refbond)*vecs(1:2,2) - (int(0.5*xy(1,refbond))*vecs(1:2,1) + xy(2,refbond)*vecs(1:2,2))
    refsite1          = bsites(1,refbond)
    refsite2          = bsites(2,refbond)      

    ! do i = 1, ndeg
        !     do j = 1, ndeg
        !         if ( abs(dot_product(psi(1:dim,i),psi(1:dim,j))) > 10.**(-10) .and. abs(dot_product(psi(1:dim,i),psi(1:dim,j))) < 0.999999999d0 ) then 
        !             print*,'currmat_c: Groundstates ',i, 'and' ,j,' are not orthorgonal.', dot_product(psi(1:dim,i),psi(1:dim,j))               
        !             error stop 
        !         end if 
        !     end do 
        ! end do 
        
    currentmat = 0.d0     
    bondcurrentmat = 0.d0  
    do l = 1, ndeg 
        do k = 1, ndeg 
            do i = 1, nbonds                            
                site1 = bsites(1,i)
                site2 = bsites(2,i)
                if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
                    .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
                coord2(1:2) = int(0.5*xy(3,i))*vecs(1:2,1) + xy(4,i)*vecs(1:2,2) - (int(0.5*xy(1,i))*vecs(1:2,1) + xy(2,i)*vecs(1:2,2))
                ! coord2(1:2) = coord2(1:2)/norm2(coord2(1:2)) 
                phase = bsites(3,i)   
                vec   = vecs(1:2,bsites(5,i))
                dist = dot_product(vec, refvec)
                psiprime = 0 

                do x = 0, Lx - 1
                    do y = 0, Ly - 1
                        do j = 1, dim
                            if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                                newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                                if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                                    ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                                    newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                                    if(gencluster == 0) then 
                                        call representative2(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                        l11 = Ly
                                        l22 = Lx
                                    else if(gencluster == 1) then 
                                        call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                        l11 = Lx
                                        l22 = Ly
                                    end if
                                    call findstate(dim, rep, basis, loc) 
                                    if(loc <= 0) cycle   
                                    parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                                    parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                    psiprime(loc) = psiprime(loc) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                                    ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                                    newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                                    if(gencluster == 0) then 
                                        call representative2(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                        l11 = Ly
                                        l22 = Lx
                                    else if(gencluster == 1) then 
                                        call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                        l11 = Lx
                                        l22 = Ly
                                    end if
                                    call findstate(dim, rep, basis, loc) 
                                    if(loc <= 0) cycle   
                                    parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                                    parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                    psiprime(loc) = psiprime(loc) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                end if     
                            else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                                newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                                if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                                    ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                                    newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                                    if(gencluster == 0) then 
                                        call representative2(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                        l11 = Ly
                                        l22 = Lx
                                    else if(gencluster == 1) then 
                                        call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                        l11 = Lx
                                        l22 = Ly
                                    end if
                                    call findstate(dim, rep, basis, loc) 
                                    if(loc <= 0) cycle   
                                    parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                                    parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                    psiprime(loc) = psiprime(loc) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                                    ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                                    newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                                    if(gencluster == 0) then 
                                        call representative2(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                        l11 = Ly
                                        l22 = Lx
                                    else if(gencluster == 1) then 
                                        call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                        l11 = Lx
                                        l22 = Ly
                                    end if
                                    call findstate(dim, rep, basis, loc) 
                                    if(loc <= 0) cycle   
                                    parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                                    parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                    psiprime(loc) = psiprime(loc) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                end if     
                            end if 
                        end do !Loop over basis states                  

                        !a1 translations 
                        cntr = 0
                        do q = 1, nbonds
                            if(bsites(1,q) == site1 .and. bsites(5,q) == 1) then 
                                site1 = bsites(2,q)
                                cntr = cntr + 1
                            else if(bsites(1,q) == site2 .and. bsites(5,q) == 1) then
                                site2 = bsites(2,q)
                                cntr = cntr + 1
                            else if(bsites(1,q) == refsite1 .and. bsites(5,q) == 1) then
                                refsite1 = bsites(2,q)
                                cntr = cntr + 1
                            else if(bsites(1,q) == refsite2 .and. bsites(5,q) == 1) then
                                refsite2 = bsites(2,q)
                                cntr = cntr + 1
                            else if(cntr == 4) then 
                                exit 
                            end if
                        end do 

                    end do 

                    !a2 translations 
                    cntr = 0
                    do q = 1, nbonds
                        if(bsites(1,q) == site1 .and. bsites(5,q) == 2) then 
                            site1 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == site2 .and. bsites(5,q) == 2) then
                            site2 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == refsite1 .and. bsites(5,q) == 2) then
                            refsite1 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == refsite2 .and. bsites(5,q) == 2) then
                            refsite2 = bsites(2,q)
                            cntr = cntr + 1
                        else if(cntr == 4) then 
                            exit 
                        end if
                    end do 
                end do 
                bondcurrentmat(i,k,l) = bondcurrentmat(i,k,l) + dble(phase) * dist * dot_product(psi(1:dim,k), psiprime)
                currentmat(k,l) = currentmat(k,l) + bondcurrentmat(i,k,l) / (Lx*Ly)
                ! bondcurrentmat(i,k,l) = bondcurrentmat(i,k,l) + dist * dot_product(psi(1:dim,k), psiprime)
                ! currentmat(k,l) = currentmat(k,l) + dble(phase) * bondcurrentmat(i,k,l)

            end do !Loop over bonds 
        end do !Loop over degenerate ground states
    end do !Loop over degenerate ground states
    deallocate(psiprime)

    call exactdiag(.True., ndeg, currentmat, currents)

    ! print*,'Current eigenvalues'
        ! do i = 1, ndeg 
        !     print*, currents(i)
        ! end do
        ! print*,'sum(currents)',sum(currents) 
        ! do i = 1, nbonds
        !     call exactdiag(.True., ndeg, bondcurrentmat(i, 1:ndeg, 1:ndeg), evals)
        !     bondcurrents(i, 1:ndeg) = evals
        !     if(allocated(evals)) deallocate(evals)
        ! end do     
        ! if(sl == "A") then 
        !     name = dir//"Acurrent_dens_mat_"//filename
        ! else if(sl == "B") then 
        !     name = dir//"Bcurrent_dens_mat_"//filename
        ! end if 
        ! name = trim_name(name)
        ! open(unit,file=name)
        ! do i = 1, ndeg 
        !     currents2(i) = sum(bondcurrents(1:nbonds, i))
        !     write(unit,*) currents2(i)
        !     ! print*,'Current of psi ',i, currents2(i)
        ! end do 
        ! close(unit)

    currentmat = currentmat / ndeg 


    call exactdiag(.True., ndeg, currentmat, currents)
    ! print*,'|===================|'
        ! print*,' Current eigenvalues'
        ! do i = 1, ndeg 
        !     print*, currents(i)
        ! end do
        ! print*,'|===================|'
        ! print*,' Total'
        ! print*,sum(currents) 
        ! currents = 0.d0 
        ! print*,'|===================|'
    if(allocated(gs)) deallocate(gs)
    allocate(gs(dim))

    do j = 1, ndeg 
        gs = 0.d0 
        do i = 1, ndeg            
            gs = gs + currentmat(i,j) * psi(1:dim, i)
        end do
        call current_c(gencluster, nHel, tilt, k1, k2, t, dim, sites, Lx, Ly, basis, bsites, xy, xtransl, ytransl, nbonds, refbond, gs, nnnVec, norm, bondcurrent, current)
        bondcurrents(1:nbonds, j) = bondcurrent 
        currents(j) = current 
        ! print*,'|===================|'
            ! print*,' Current, eigenstate'
            ! print*, current, j
            ! print*,'|===================|'
        current = 0.d0 
    end do 
    ! print*,'|===================|'
        ! print*,' Sum'
        ! print*, sum(currents)
        ! print*,'|===================|'


    if(sl == "A") then 
        name = trim_name(dir//"Acurrent_matrix_"//filename)
    else if(sl == "B") then 
        name = trim_name(dir//"Bcurrent_matrix_"//filename)
    end if 
    name = trim_name(name)
    print* ,name, 'name'
    open(unit,file=name)
    do i = 1, ndeg 
        write(unit,*) currents(i)
    end do 
    close(unit)
    if(sl == "A") then 
        name = trim_name(dir//"Abondcurrent_matrix_"//filename)
    else if(sl == "B") then 
        name = trim_name(dir//"Bbondcurrent_matrix_"//filename)
    end if 
    name = trim_name(name)
    open(unit,file=name)
    do i = 1, ndeg 
        write(unit,*) bondcurrents(1:nbonds,i)
    end do 
    close(unit)


    ! bondcurrents = 0.d0 
        ! do i = 1, nbonds
        !     call exactdiag(.True., ndeg, bondcurrentmat(i, 1:ndeg, 1:ndeg), evals)
        !     bondcurrents(i, 1:ndeg) = evals
        !     if(allocated(evals)) deallocate(evals)
        ! end do 
            
        ! if(sl == "A") then 
        !     name = dir//"Acurrent_matrix2_"//filename
        ! else if(sl == "B") then 
        !     name = dir//"Bcurrent_matrix2_"//filename
        ! end if 
        ! name = trim_name(name)
        ! open(unit,file=name)
        ! do i = 1, ndeg 
        !     currents2(i) = sum(bondcurrents(1:nbonds, i))
        !     write(unit,*) currents2(i)
        !     print*,'Current of psi ',i, currents2(i)
        ! end do 
        ! close(unit)

    return

end subroutine currmat_c2

subroutine current_ics_c2(gencluster, nHel, tilt, k1, k2, ndeg, t, dim, sites, Lx, Ly, basis, bsites, xy, xtransl, ytransl, nbonds, refbond, psi, nnnVec, norm, bondcurrent, current)

    implicit none
    integer, intent(in) :: gencluster, nHel, tilt, k1, k2, ndeg, nbonds, refbond, sites, Lx, Ly  
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    integer, intent(in) :: xy(4,nbonds)
    integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
    double precision, intent(in) :: t, nnnVec(2, 3) 
    double complex, intent(in) :: psi(dim, ndeg)
    double precision, intent(in) :: norm(dim)
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bondcurrent(:)
    
    integer(kind=8) :: loc = 0, newst = 0, rep = 0
    integer :: i = 0, j = 0, l = 0, l1 = 0, l2 = 0, x = 0, y = 0, q = 0, cntr = 0
    integer :: sign = 0, signt = 0, l11 = 0, l22 = 0  
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec, coord1, coord2 
    double precision, dimension(2,3) :: vecs
    double complex, allocatable :: psiprime(:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bondcurrent)) deallocate(bondcurrent)
    allocate(bondcurrent(nbonds))
    if( allocated(psiprime) ) deallocate(psiprime)
    allocate( psiprime(dim) )

    psiprime = 0.d0
    refvec   = vecs(1:2,bsites(5,refbond))
    coord1(1:2)    = int(0.5*xy(3,refbond))*vecs(1:2,1) + xy(4,refbond)*vecs(1:2,2) - (int(0.5*xy(1,refbond))*vecs(1:2,1) + xy(2,refbond)*vecs(1:2,2))
    refsite1          = bsites(1,refbond)
    refsite2          = bsites(2,refbond)      
   
    ! do i = 1, ndeg
        !     do j = 1, ndeg
        !         if ( abs(dot_product(psi(1:dim,i),psi(1:dim,j))) > 10.**(-10) .and. abs(dot_product(psi(1:dim,i),psi(1:dim,j))) < 0.999999999d0 ) then 
        !             print*,'current_ics_c: Groundstates ',i, 'and' ,j,' are not orthorgonal.', dot_product(psi(1:dim,i),psi(1:dim,j))               
        !             error stop 
        !         end if 
        !     end do 
        ! end do 

    bondcurrent = 0.d0 
    current = 0.d0     
    do i = 1, nbonds            
        bondcurrent(i) = 0 
        site1 = bsites(1,i)
        site2 = bsites(2,i)
        if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
            .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
        coord2(1:2) = int(0.5*xy(3,i))*vecs(1:2,1) + xy(4,i)*vecs(1:2,2) - (int(0.5*xy(1,i))*vecs(1:2,1) + xy(2,i)*vecs(1:2,2))
        phase = bsites(3,i)   
        vec   = vecs(1:2,bsites(5,i))
        dist = dot_product(vec, refvec)

        do l = 1, ndeg 
            psiprime = 0 
            do x = 0, Lx - 1
                do y = 0, Ly - 1
                    do j = 1, dim
                        if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                            newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                            if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                                ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                                newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                                if(gencluster == 0) then 
                                    call representative2(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                    l11 = Ly
                                    l22 = Lx
                                else if(gencluster == 1) then 
                                    call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                    l11 = Lx
                                    l22 = Ly
                                end if
                                call findstate(dim, rep, basis, loc) 
                                if(loc <= 0) cycle   
                                
                                parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                                parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                                parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                                parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                                sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                psiprime(loc) = psiprime(loc) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                            else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                                ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                                newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                                if(gencluster == 0) then 
                                    call representative2(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                    l11 = Ly
                                    l22 = Lx
                                else if(gencluster == 1) then 
                                    call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                    l11 = Lx
                                    l22 = Ly
                                end if
                                call findstate(dim, rep, basis, loc) 
                                if(loc <= 0) cycle   
                                
                                parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                                parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                                parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                                parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                                sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                psiprime(loc) = psiprime(loc) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                            end if     
                        else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                            newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                            if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                                ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                                newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                                if(gencluster == 0) then 
                                    call representative2(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                    l11 = Ly
                                    l22 = Lx
                                else if(gencluster == 1) then 
                                    call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                    l11 = Lx
                                    l22 = Ly
                                end if
                                call findstate(dim, rep, basis, loc) 
                                if(loc <= 0) cycle   
                                
                                parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                                parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                                parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                                parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                                sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                psiprime(loc) = psiprime(loc) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                            else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                                ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                                newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                                if(gencluster == 0) then 
                                    call representative2(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                    l11 = Ly
                                    l22 = Lx
                                else if(gencluster == 1) then 
                                    call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                    l11 = Lx
                                    l22 = Ly
                                end if
                                call findstate(dim, rep, basis, loc) 
                                if(loc <= 0) cycle   
                                
                                parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                                parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                                parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                                parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                                sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                psiprime(loc) = psiprime(loc) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                            end if     
                        end if 
                    end do !Loop over basis states      
                    
                    !a1 translations 
                    cntr = 0
                    do q = 1, nbonds
                        if(bsites(1,q) == site1 .and. bsites(5,q) == 1) then 
                            site1 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == site2 .and. bsites(5,q) == 1) then
                            site2 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == refsite1 .and. bsites(5,q) == 1) then
                            refsite1 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == refsite2 .and. bsites(5,q) == 1) then
                            refsite2 = bsites(2,q)
                            cntr = cntr + 1
                        else if(cntr == 4) then 
                            exit 
                        end if
                    end do 

                end do 

                !a2 translations 
                cntr = 0
                do q = 1, nbonds
                    if(bsites(1,q) == site1 .and. bsites(5,q) == 2) then 
                        site1 = bsites(2,q)
                        cntr = cntr + 1
                    else if(bsites(1,q) == site2 .and. bsites(5,q) == 2) then
                        site2 = bsites(2,q)
                        cntr = cntr + 1
                    else if(bsites(1,q) == refsite1 .and. bsites(5,q) == 2) then
                        refsite1 = bsites(2,q)
                        cntr = cntr + 1
                    else if(bsites(1,q) == refsite2 .and. bsites(5,q) == 2) then
                        refsite2 = bsites(2,q)
                        cntr = cntr + 1
                    else if(cntr == 4) then 
                        exit 
                    end if
                end do 
            end do 
            
            bondcurrent(i) = bondcurrent(i) + dist * dot_product(psi(1:dim,l), psiprime)
        end do !Loop over degenerate ground states
        current = current + dble(phase) * bondcurrent(i) / (Lx*Ly)
    end do !Loop over bonds 

    current = current/dble(ndeg)!sqrt(dble(ndeg))
    deallocate(psiprime)
 
    
    return

end subroutine current_ics_c2

subroutine current_cs_c2(gencluster, nHel, tilt, dir, unit, parameters, sl, ndeg, k1, k2, t, dim, sites, Lx, Ly, basis, bsites, xy, xtransl, ytransl, nbonds, refbond, psi, nnnVec, norm)

    implicit none
    integer, intent(in) :: gencluster, nHel, tilt, unit, ndeg, k1, k2, nbonds, refbond, sites, Lx, Ly  
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    integer, intent(in) :: xy(4,nbonds)
    integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
    double precision, intent(in) :: t, nnnVec(2, 3)
    double complex, intent(in) :: psi(dim, ndeg)
    double precision, intent(in) :: norm(dim)
    character(len=*), intent(in) :: dir, sl, parameters
    
    integer(kind=8) :: loc = 0, newst = 0, rep = 0
    integer :: i = 0, j = 0, k = 0, l = 0, l1 = 0, l2 = 0, x = 0, y = 0, q = 0, cntr = 0
    integer :: sign = 0, signt = 0, l11 = 0, l22 = 0  
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec, coord1, coord2 
    double precision, dimension(2,3) :: vecs
    double precision :: current
    double precision, allocatable :: bondcurrent(:)
    double complex, allocatable :: psiprime(:)
    character :: filename*500 
    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bondcurrent)) deallocate(bondcurrent)
    allocate(bondcurrent(nbonds))
    if( allocated(psiprime) ) deallocate(psiprime)
    allocate( psiprime(dim) )

    psiprime = 0.d0
    refvec   = vecs(1:2,bsites(5,refbond))
    coord1(1:2)    = int(0.5*xy(3,refbond))*vecs(1:2,1) + xy(4,refbond)*vecs(1:2,2) - (int(0.5*xy(1,refbond))*vecs(1:2,1) + xy(2,refbond)*vecs(1:2,2))
    refsite1          = bsites(1,refbond)
    refsite2          = bsites(2,refbond)      

    ! do i = 1, ndeg
        !     do j = 1, ndeg
        !         if ( abs(dot_product(psi(1:dim,i),psi(1:dim,j))) > 10.**(-10) .and. abs(dot_product(psi(1:dim,i),psi(1:dim,j))) < 0.999999999d0 ) then 
        !             print*,'current_cs_c: Groundstates ',i, 'and' ,j,' are not orthorgonal.', dot_product(psi(1:dim,i),psi(1:dim,j))               
        !             error stop 
        !         end if 
        !     end do 
        ! end do 
    bondcurrent = 0.d0
    current = 0.d0     
    do i = 1, nbonds            
        bondcurrent(i) = 0 
        site1 = bsites(1,i)
        site2 = bsites(2,i)
        if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
            .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
        coord2(1:2) = int(0.5*xy(3,i))*vecs(1:2,1) + xy(4,i)*vecs(1:2,2) - (int(0.5*xy(1,i))*vecs(1:2,1) + xy(2,i)*vecs(1:2,2))
        ! coord2(1:2) = coord2(1:2)/norm2(coord2(1:2)) 
        phase = bsites(3,i)   
        vec   = vecs(1:2,bsites(5,i))
        dist = dot_product(vec, refvec)

        do l = 1, ndeg 
            do k = 1, ndeg 
                psiprime = 0 
                do x = 0, Lx - 1
                    do y = 0, Ly - 1
                        do j = 1, dim
                            if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                                newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                                if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                                    ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                                    newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                                    if(gencluster == 0) then 
                                        call representative2(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                        l11 = Ly
                                        l22 = Lx
                                    else if(gencluster == 1) then 
                                        call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                        l11 = Lx
                                        l22 = Ly
                                    end if  
                                    call findstate(dim, rep, basis, loc) 
                                    if(loc <= 0) cycle   
                                    parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                                    parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                    psiprime(loc) = psiprime(loc) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                                    ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                                    newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                                    if(gencluster == 0) then 
                                            call representative2(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                            l11 = Ly
                                            l22 = Lx
                                        else if(gencluster == 1) then 
                                            call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                            l11 = Lx
                                            l22 = Ly
                                        end if  
                                    call findstate(dim, rep, basis, loc) 
                                    if(loc <= 0) cycle   
                                    parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                                    parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                    psiprime(loc) = psiprime(loc) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                end if     
                            else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                                newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                                if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                                    ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                                    newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                                    if(gencluster == 0) then 
                                        call representative2(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                        l11 = Ly
                                        l22 = Lx
                                    else if(gencluster == 1) then 
                                        call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                        l11 = Lx
                                        l22 = Ly
                                    end if  
                                    call findstate(dim, rep, basis, loc) 
                                    if(loc <= 0) cycle   
                                    parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                                    parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                    psiprime(loc) = psiprime(loc) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                                    ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                                    newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                                    if(gencluster == 0) then 
                                            call representative2(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                            l11 = Ly
                                            l22 = Lx
                                        else if(gencluster == 1) then 
                                            call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                            l11 = Lx
                                            l22 = Ly
                                        end if  
                                    call findstate(dim, rep, basis, loc) 
                                    if(loc <= 0) cycle   
                                    parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                                    parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                    psiprime(loc) = psiprime(loc) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                end if     
                            end if 
                        end do !Loop over basis states    
                        !a1 translations 
                        cntr = 0
                        do q = 1, nbonds
                            if(bsites(1,q) == site1 .and. bsites(5,q) == 1) then 
                                site1 = bsites(2,q)
                                cntr = cntr + 1
                            else if(bsites(1,q) == site2 .and. bsites(5,q) == 1) then
                                site2 = bsites(2,q)
                                cntr = cntr + 1
                            else if(bsites(1,q) == refsite1 .and. bsites(5,q) == 1) then
                                refsite1 = bsites(2,q)
                                cntr = cntr + 1
                            else if(bsites(1,q) == refsite2 .and. bsites(5,q) == 1) then
                                refsite2 = bsites(2,q)
                                cntr = cntr + 1
                            else if(cntr == 4) then 
                                exit 
                            end if
                        end do 

                    end do 

                    !a2 translations 
                    cntr = 0
                    do q = 1, nbonds
                        if(bsites(1,q) == site1 .and. bsites(5,q) == 2) then 
                            site1 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == site2 .and. bsites(5,q) == 2) then
                            site2 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == refsite1 .and. bsites(5,q) == 2) then
                            refsite1 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == refsite2 .and. bsites(5,q) == 2) then
                            refsite2 = bsites(2,q)
                            cntr = cntr + 1
                        else if(cntr == 4) then 
                            exit 
                        end if
                    end do 
                end do          
                bondcurrent(i) = bondcurrent(i) + dist * dot_product(psi(1:dim,k), psiprime)
            end do !Loop over degenerate ground states
        end do !Loop over degenerate ground states
        current = current + dble(phase) * bondcurrent(i) / (Lx*Ly)
    end do !Loop over bonds 

    deallocate(psiprime)

    if(sl == "A") then 
        filename = trim_name(dir//"Acurrent_coh_sum_"//parameters)
    else if(sl == "B") then 
        filename = trim_name(dir//"Bcurrent_coh_sum_"//parameters)
    end if 
    filename = trim_name(filename)
    print* ,filename, 'name'
    open(unit,file=filename)
    write(unit,*) current 
    close(unit)
    
    return

end subroutine current_cs_c2

subroutine currmat(dir, unit, filename, sl, ndeg, t, dim, sites, basis, bsites, xy, nbonds, refbond, psi, nnnVec)

    implicit none
    integer, intent(in) :: unit, ndeg, nbonds, refbond, sites
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    integer, intent(in) :: xy(4,nbonds)
    double precision, intent(in) :: t 
    double precision, intent(in) :: psi(dim, ndeg), nnnVec(2, 3)
    character(len=*), intent(in) :: dir, sl, filename
    
    integer(kind=8) :: loc = 0, newst = 0
    integer :: i = 0, j = 0, k = 0, l = 0
    integer :: sign = 0
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0, current = 0.d0 
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:), gs(:)
    double precision, allocatable :: currents(:)
    double precision, allocatable :: currents2(:)
    double precision, allocatable :: evals(:)
    double precision, allocatable :: bondcurrents(:,:), bondcurrent(:)
    double precision, allocatable :: currentmat(:,:)
    double precision, allocatable :: bondcurrentmat(:,:,:)
    character*500 :: name 

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)

    !Definition of lattice vectors: (5th row of bsites)
    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(psiprime)) deallocate(psiprime)
    if(allocated(bondcurrentmat)) deallocate(bondcurrentmat)
    if(allocated(currentmat)) deallocate(currentmat)
    if(allocated(currents2)) deallocate(currents2)
    if(allocated(bondcurrents)) deallocate(bondcurrents)
    allocate(psiprime(dim))
    allocate(bondcurrentmat(nbonds, ndeg, ndeg))
    allocate(currentmat(ndeg, ndeg))
    allocate(currents2(ndeg))
    allocate(bondcurrents(nbonds, ndeg))

    psiprime = 0.d0
    refvec   = vecs(1:2,bsites(5,refbond))
    refsite1 = bsites(1,refbond)
    refsite2 = bsites(2,refbond)      

    ! do i = 1, ndeg
        !     do j = 1, ndeg
        !         if ( abs(dot_product(psi(1:dim,i),psi(1:dim,j))) > 10.**(-10) .and. abs(dot_product(psi(1:dim,i),psi(1:dim,j))) < 0.999999999d0 ) then 
        !             print*,'Currmat: Groundstates ',i, 'and' ,j,' are not orthorgonal.', dot_product(psi(1:dim,i),psi(1:dim,j))               
        !             error stop 
        !         end if 
        !     end do 
        ! end do 
    
    currentmat = 0.d0     
    bondcurrentmat = 0.d0  
    do l = 1, ndeg 
        do k = 1, ndeg 
            do i = 1, nbonds                            
                site1 = bsites(1,i)
                site2 = bsites(2,i)
                if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
                    .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
                phase = bsites(3,i)   
                vec   = vecs(1:2,bsites(5,i))
                dist  = dot_product(vec, refvec)
                psiprime = 0.d0
                do j = 1, dim
                    if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                        newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                        if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                            ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                            newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                            call findstate(dim, newst, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j, l)
                        else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                            ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                            newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                            call findstate(dim, newst, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + (-1) * sign * (-(t)**2) * psi(j, l)
                        end if     
                    else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                        newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                        if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                            ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                            newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                            call findstate(dim, newst, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + (-1) * sign * (-(t)**2) * psi(j, l)
                        else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                            ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                            newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                            call findstate(dim, newst, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j, l)
                        end if     
                    end if 
                end do !Loop over basis states                  
                bondcurrentmat(i,k,l) = bondcurrentmat(i,k,l) + dble(phase) * dist * dot_product(psi(1:dim,k), psiprime)
                ! bondcurrentmat(i,k,l) = bondcurrentmat(i,k,l) + dist * dot_product(psi(1:dim,k), psiprime)
                currentmat(k,l) = currentmat(k,l) + bondcurrentmat(i,k,l)
                ! currentmat(k,l) = currentmat(k,l) + dble(phase) * bondcurrentmat(i,k,l)
            end do !Loop over bonds 
        end do !Loop over degenerate ground states
    end do !Loop over degenerate ground states
    deallocate(psiprime)

    currentmat = currentmat / ndeg 


    call exactdiag(.True., ndeg, currentmat, currents)
print*,'|===================|'
    print*,' Current eigenvalues'
    do i = 1, ndeg 
        print*, currents(i)
    end do
    print*,'|===================|'
    print*,' Total'
    print*,sum(currents) 
    ! currents = 0.d0 
    print*,'|===================|'
    if(allocated(gs)) deallocate(gs)
    allocate(gs(dim))

    gs = 0.d0 
    do i = 1, ndeg 
        gs = gs + currents(i) * psi(1:dim, i)
        print* ,norm2(psi(1:dim,i)), 'norm2(psi(1:dim,i))'
    end do
    print* ,norm2(gs), 'norm2(gs)'
    gs = gs / sqrt(dble(ndeg))
    print* ,norm2(gs), 'norm2(gs)'

    call slcurrent(1.0d0, dim, sites, basis, bsites, xy, nbonds, refbond, gs, nnnVec, bondcurrent, current)
    
    print*,'|===================|'
    print*,' Current, manifold'
    print*, current
    print*,'|===================|'
    ! do j = 1, ndeg 
    !     gs = 0.d0 
    !     do i = 1, ndeg 
    !         gs = gs + currentmat(i,j) * psi(1:dim, i)
    !     end do
    !     call slcurrent(1.0d0, dim, sites, basis, bsites, xy, nbonds, refbond, gs, nnnVec, bondcurrent, current)
    !     bondcurrents(1:nbonds, j) = bondcurrent 
    !     currents(j) = current 
    !     print*,'|===================|'
    !     print*,' Current, eigenstate'
    !     print*, current, j
    !     print*,'|===================|'
    !     current = 0.d0 
    ! end do 
    
    if(sl == "A") then 
        name = trim_name(dir//"Acurrent_matrix_"//filename)
    else if(sl == "B") then 
        name = trim_name(dir//"Bcurrent_matrix_"//filename)
    end if 
    name = trim_name(name)
    
    open(unit,file=name)
    do i = 1, ndeg 
        write(unit,*) currents(i)
    end do 
    close(unit)
    if(sl == "A") then 
        name = dir//"Abondcurrent_matrix_"//filename
    else if(sl == "B") then 
        name = dir//"Bbondcurrent_matrix_"//filename
    end if 
    name = trim_name(name)
    open(unit,file=name)
    do i = 1, ndeg 
        write(unit,*) bondcurrents(1:nbonds,i)
    end do 
    close(unit)

    ! bondcurrents = 0.d0 
        ! do i = 1, nbonds
        !     call exactdiag(.True., ndeg, bondcurrentmat(i, 1:ndeg, 1:ndeg), evals)
        !     bondcurrents(i, 1:ndeg) = evals
        !     if(allocated(evals)) deallocate(evals)
        ! end do 
            
        ! if(sl == "A") then 
        !     name = dir//"Acurrent_matrix2_"//filename
        ! else if(sl == "B") then 
        !     name = dir//"Bcurrent_matrix2_"//filename
        ! end if 
        ! name = trim_name(name)
        ! open(unit,file=name)
        ! do i = 1, ndeg 
        !     currents2(i) = sum(bondcurrents(1:nbonds, i))
        !     write(unit,*) currents2(i)
        !     print*,'Current of psi ',i, currents2(i)
        ! end do 
        ! close(unit)

    return

end subroutine currmat

subroutine current_ics(ndeg, t, dim, sites, basis, bsites, nbonds, refbond, psi, nnnVec, bondcurrent, current)

    implicit none
    integer, intent(in) :: ndeg, nbonds, refbond, sites
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    double precision, intent(in) :: t 
    double precision, intent(in) :: psi(dim, ndeg), nnnVec(2, 3)
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bondcurrent(:)
    
    integer(kind=8) :: loc = 0, newst = 0
    integer :: i = 0, j = 0, l = 0
    integer :: sign = 0
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)
    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bondcurrent)) deallocate(bondcurrent)
    if(allocated(psiprime)) deallocate(psiprime)
    allocate(bondcurrent(nbonds))
    allocate(psiprime(dim))

    psiprime = 0.d0
    refvec   = vecs(1:2,bsites(5,refbond))
    refsite1 = bsites(1,refbond)
    refsite2 = bsites(2,refbond)      

    current = 0     
    do i = 1, nbonds            
        bondcurrent(i) = 0 
        site1 = bsites(1,i)
        site2 = bsites(2,i)
        if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
            .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
        phase = bsites(3,i)   
        vec   = vecs(1:2,bsites(5,i))
        dist = dot_product(vec, refvec)

        do l = 1, ndeg  
            psiprime = 0 
            do j = 1, dim
                if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                    newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                    if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                        ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                        newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                        call findstate(dim, newst, basis, loc) 
                        if(loc <= 0) cycle   
                        
                        parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                        parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                        parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                        parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                        sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                        psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j, l)
                    else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                        ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                        newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                        call findstate(dim, newst, basis, loc) 
                        if(loc <= 0) cycle   
                        
                        parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                        parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                        parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                        parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                        sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                        psiprime(loc) = psiprime(loc) + (-1) * sign * (-(t)**2) * psi(j, l)
                    end if     
                else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                    newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                    if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                        ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                        newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                        call findstate(dim, newst, basis, loc) 
                        if(loc <= 0) cycle   
                        
                        parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                        parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                        parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                        parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                        sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                        psiprime(loc) = psiprime(loc) + (-1) * sign * (-(t)**2) * psi(j, l)
                    else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                        ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                        newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                        call findstate(dim, newst, basis, loc) 
                        if(loc <= 0) cycle   
                        
                        parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                        parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                        parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                        parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                        sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                        psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j, l)
                    end if     
                end if 
            end do !Loop over basis states                  
            bondcurrent(i) = bondcurrent(i) + dist * dot_product(psi(1:dim,l), psiprime)
        
        end do !Loop over degenerate ground states
        current = current + dble(phase) * bondcurrent(i)
    end do !Loop over bonds 

    deallocate(psiprime)
    
    return

end subroutine current_ics

subroutine current_cs(dir, unit, parameters, sl, ndeg, t, dim, sites, basis, bsites, nbonds, refbond, psi, nnnVec)

    implicit none
    integer, intent(in) :: unit, ndeg, nbonds, refbond, sites
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    double precision, intent(in) :: t 
    double precision, intent(in) :: psi(dim, ndeg), nnnVec(2, 3)
    character(len=*), intent(in) :: dir, sl, parameters
    
    integer(kind=8) :: loc = 0, newst = 0
    integer :: i = 0, j = 0, k = 0, l = 0
    integer :: sign = 0
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:)
    double precision :: current
    double precision, allocatable :: bondcurrent(:)
    character*400 :: name

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)
    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bondcurrent)) deallocate(bondcurrent)
    if(allocated(psiprime)) deallocate(psiprime)
    allocate(bondcurrent(nbonds))
    allocate(psiprime(dim))

    psiprime    = 0.d0
    refvec      = vecs(1:2,bsites(5,refbond))
    refsite1    = bsites(1,refbond)
    refsite2    = bsites(2,refbond)      

    do i = 1, ndeg
        do j = 1, ndeg
            if ( abs(dot_product(psi(1:dim,i),psi(1:dim,j))) > 10.**(-10) .and. abs(dot_product(psi(1:dim,i),psi(1:dim,j))) < 0.999999999d0 ) then 
                print*,'current_cs coherent sum: Groundstates ',i, 'and' ,j,' are not orthorgonal.', dot_product(psi(1:dim,i),psi(1:dim,j))               
                error stop 
            end if 
        end do 
    end do 
    current = 0     
    do i = 1, nbonds            
        bondcurrent(i) = 0 
        site1 = bsites(1,i)
        site2 = bsites(2,i)
        if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
            .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
        phase = bsites(3,i)   
        vec   = vecs(1:2,bsites(5,i))
        dist = dot_product(vec, refvec)

        do l = 1, ndeg 
            do k = 1, ndeg 
                psiprime = 0 
                do j = 1, dim
                    if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                        newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                        if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                            ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                            newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                            call findstate(dim, newst, basis, loc) 
                            if(loc <= 0) cycle   
                            if( basis(loc) .ne. newst) error stop "Subroutine current_cs coherent sum: Wrong location of scattered basis state."
                            parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j, l)
                        else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                            ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                            newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                            call findstate(dim, newst, basis, loc) 
                            if(loc <= 0) cycle   
                            if( basis(loc) .ne. newst) error stop "Subroutine current_cs coherent sum: Wrong location of scattered basis state."
                            parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + (-1) * sign * (-(t)**2) * psi(j, l)
                        end if     
                    else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                        newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                        if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                            ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                            newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                            call findstate(dim, newst, basis, loc) 
                            if(loc <= 0) cycle   
                            if( basis(loc) .ne. newst) error stop "Subroutine current_cs coherent sum: Wrong location of scattered basis state."
                            parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + (-1) * sign * (-(t)**2) * psi(j, l)
                        else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                            ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                            newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                            call findstate(dim, newst, basis, loc) 
                            if(loc <= 0) cycle   
                            if( basis(loc) .ne. newst) error stop "Subroutine current_cs coherent sum: Wrong location of scattered basis state."
                            parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j, l)
                        end if     
                    end if 
                end do !Loop over basis states                  
                bondcurrent(i) = bondcurrent(i) + dist * dot_product(psi(1:dim,k), psiprime)
            end do !Loop over degenerate ground states
        end do !Loop over degenerate ground states
        current = current + dble(phase) * bondcurrent(i)
    end do !Loop over bonds 
    print* ,current, 'coherent sum '
    

    if(sl == "A") then 
        name = trim_name(dir//"Acurrent_coh_sum_"//parameters)
    else if(sl == "B") then 
        name = trim_name(dir//"Bcurrent_coh_sum_"//parameters)
    end if 
    open(unit,file=name)
    write(unit,*) current
    close(unit)

    if(sl == "A") then 
        name = trim_name(dir//"Abondcurrent_coh_sum_"//parameters)
    else if(sl == "B") then 
        name = trim_name(dir//"Bbondcurrent_coh_sum_"//parameters)
    end if 
    name = trim_name(name)
    open(unit,file=name)
    write(unit,*) bondcurrent
    close(unit)

    deallocate(psiprime)
    
    return

end subroutine current_cs







end module utilities_old    ! Module for utility functions and subroutines
