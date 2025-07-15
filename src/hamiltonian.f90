module hamiltonian 

    use printing_routines
    use io_routines
    use basis

    implicit none 

    interface unify 
        module procedure unify, unify_dp, unify_dc, unify_dense, unify_dense_dp, unify_dense_dc
    end interface

    interface hopping 
        module procedure hopping_i, hopping_dc, hopping_irrep, hopping_irrep2D
    end interface

    interface diagonal 
        module procedure diagonal, diagonal_dp, diagonal_irrep2D, diagonal_p
    end interface

    contains 


    !------------------------------------------!
    !               Hamiltonian                !
    !------------------------------------------!

    subroutine generate_hamiltonian()
    ! subroutine generate_hamiltonian(thrds, ti, unit, prms, sites, nbb, trisites, dim, basis, hamOff, nOff, k1, k2, tilted, nHel, tilt, lx, ly, ucx, ucy, orbsize, norm, norm2D, orbits, phases, xtransl, ytransl, symm, id, par, rot, refl, c6, t, rcoff, rcdi, parities, dplcts, hamOff_dp, hamDi_dp, nDi_off, hamOff_dc, hamDi_dc, hamDi_c2, nnnbb, hexsites, hamDi, occ, nDi)
        use variables
        use input_variables
        implicit none 

        ! integer, intent(in)                          :: thrds, ti, ucx, ucy, tilted, nHel, tilt, lx, ly, k1, k2, symm, orbsize, refl(6, sites), c6(sites)
        ! double precision, intent(in)                 :: t, id, par(6), rot(5) 

        ! integer, intent(inout)                       :: unit, sites, nbb, nnnbb
        ! integer, allocatable, intent(inout)          :: parities(:), xtransl(:,:), ytransl(:,:), trisites(:,:), hexsites(:,:)
        ! integer(kind=8), intent(inout)               :: dim, nOff, nDi_off, nDi 
        ! integer(kind=8), allocatable, intent(inout)  :: basis(:), orbits(:,:,:), hamDi(:,:), occ(:,:), hamOff(:,:), rcoff(:,:), rcdi(:), dplcts(:)
        ! double precision, allocatable, intent(inout) :: norm(:), norm2D(:,:), hamDi_dp(:), hamOff_dp(:)
        ! double complex, allocatable, intent(inout)   :: hamOff_dc(:), hamDi_dc(:), hamDi_c2(:,:), phases(:,:,:)
        ! character(len=*), intent(inout)              ::         ! Local variables
        ! character :: type*1 
        integer :: flag_i, flag_2D
        double precision :: flag_dp
        logical :: exist1, exist2, exist3

        print*, 'Generating Hamiltonian ...'
        print*, ''

        call print_parameters()
        call loadham(type, unit, ti, sites, nOff, nDi, hamDi, hamOff, rcoff, hamOff_dp, hamOff_dc, occ, exist1, exist2, exist3) ! Loads the sub-Hamiltonians from file if they exists. If not, sets corresponDing exist1, exist2, exist3 to .false. and generates the Hamiltonian.

        if(.not.(exist1)) then ! Generates the off-diagonal part of the Hamiltonian if it does not exist yet.
            if(ti == 0) then ! If no translation symmetry is used, generate the Hamiltonian with the standard hopping procedure
                ! call hopping(flag_i)
                call hopping_i(othrds, unit, sites, nnBonds, bsites, dim, basis_states, hamOff, nOff)
            else if(ti == 1 .and. k1 == 0 .and. k2 == 0) then ! If translation symmetry is used and k1 = k2 = 0, generate the Hamiltonian with the hopping procedure for irreducible representations
                
                if(id == 1) then ! Hopping for irreducible representations in 1D 
                    ! call hopping(flag_dp)
                    call hopping_irrep(othrds, tilted, nHel, tilt, l2, l1, ucx, ucy, sites, nnBonds, dim, basis_states, bsites, norm, xtransl, ytransl, symm, id, mir, rot, refl, c6, t, rcoff, rcdi, parities, dplcts, hamOff_dp, hamDi_dp, nOff, nDiOff)
                else if(id == 2) then ! Hopping for irreducible representations in 2D (currently not supported)
                    ! call hopping(flag_dp, flag_2D)
                    call hopping_irrep2D(othrds, tilted, nHel, tilt, l2, l1, sites, nnBonds, dim, basis_states, orbsize, orbits2D, phases2D, norm2D, bsites, xtransl, ytransl, symm, id, mir, rot, refl, c6, t, rcoff, rcdi, parities, dplcts, hamOff_dc, hamDi_dc_2, nOff, nDiOff)
                    ! call diagonal(flag_dp, flag_2D)
                    call diagonal_irrep2D(sites, dim, bsites, hexsites, orbsize, orbits2D, phases2D, norm2D, hamDi_dc, occ)
                end if 
            else if(ti == 1 .and. ((k1 .ne. 0) .or. (k2 .ne. 0))) then ! If translation symmetry is used and k1 or k2 are not zero, generate the Hamiltonian with the hopping procedure for complex entries. Currently only supported for 1D irreps. 
                call hopping(unit, tilted, nHel, tilt, l2, l1, ucx, ucy, sites, nnBonds, dim, basis_states, bsites, norm, xtransl, ytransl, t, k1, k2, rcoff, hamOff_dc, nOff)
            end if 
        end if
        
        !Generates diagonal part of the Hamiltonian if it does not exist yet.
        if(id == 1 .and. (.not.(exist2) .or. .not.(exist3))) call diagonal(unit, sites, nnBonds, nnnBonds, dim, bsites, hexsites, basis_states, hamDi, occ, nDi) 
        
        ! Saves the Hamiltonian to file if it does not exist yet.
        if(id == 1 .and. (.not.(exist1) .or. .not.(exist2) .or. .not.(exist3))) call save_ham(type, unit, ti, sites, nOff, nDi, hamDi, hamOff, rcoff, hamOff_dp, hamOff_dc, occ)

        return 

    end subroutine generate_hamiltonian

    ! With variable arguments
    ! subroutine generate_hamiltonian(thrds, ti, unit, prms, sites, nbb, trisites, dim, basis, hamOff, nOff, k1, k2, tilted, nHel, tilt, lx, ly, ucx, ucy, orbsize, norm, norm2D, orbits, phases, xtransl, ytransl, symm, id, par, rot, refl, c6, t, rcoff, rcdi, parities, dplcts, hamOff_dp, hamDi_dp, nDi_off, hamOff_dc, hamDi_dc, hamDi_c2, nnnbb, hexsites, hamDi, occ, nDi)

        !     implicit none 

        !     integer, intent(in)                          :: thrds, ti, ucx, ucy, tilted, nHel, tilt, lx, ly, k1, k2, symm, orbsize, refl(6, sites), c6(sites)
        !     double precision, intent(in)                 :: t, id, par(6), rot(5) 

        !     integer, intent(inout)                       :: unit, sites, nbb, nnnbb
        !     integer, allocatable, intent(inout)          :: parities(:), xtransl(:,:), ytransl(:,:), trisites(:,:), hexsites(:,:)
        !     integer(kind=8), intent(inout)               :: dim, nOff, nDi_off, nDi 
        !     integer(kind=8), allocatable, intent(inout)  :: basis(:), orbits(:,:,:), hamDi(:,:), occ(:,:), hamOff(:,:), rcoff(:,:), rcdi(:), dplcts(:)
        !     double precision, allocatable, intent(inout) :: norm(:), norm2D(:,:), hamDi_dp(:), hamOff_dp(:)
        !     double complex, allocatable, intent(inout)   :: hamOff_dc(:), hamDi_dc(:), hamDi_c2(:,:), phases(:,:,:)
        !     character(len=*), intent(inout)              :: prms

        !     ! Local variables
        !     character :: type*1 
        !     logical   :: exist1, exist2, exist3

        !     print*, 'Generating Hamiltonian ...'
        !     print*, ''
        !      call print_parameters()
        !!     call par2(unit, k1, k2, -1, 0.d0, 0.!d0, prms, type)
        !     call loadham(type, prms, unit, ti, sites, nOff, nDi, hamDi, hamOff, rcoff, hamOff_dp, hamOff_dc, occ, exist1, exist2, exist3) ! Loads the sub-Hamiltonians from file if they exists. If not, sets corresponDing exist1, exist2, exist3 to .false. and generates the Hamiltonian.

        !     if(.not.(exist1)) then ! Generates the off-diagonal part of the Hamiltonian if it does not exist yet.
        !         if(ti == 0) then ! If no translation symmetry is used, generate the Hamiltonian with the standard hopping procedure
        !             call hopping(thrds, unit, prms, sites, nbb, trisites, dim, basis, hamOff, nOff)
        !         else if(ti == 1 .and. k1 == 0 .and. k2 == 0) then ! If translation symmetry is used and k1 = k2 = 0, generate the Hamiltonian with the hopping procedure for irreducible representations
                    
        !             if(id == 1) then ! Hopping for irreducible representations in 1D 
        !                 call hopping(thrds, tilted, nHel, tilt, lx, ly, ucx, ucy, sites, nbb, dim, basis, trisites, norm, xtransl, ytransl, symm, id, par, rot, refl, c6, t, rcoff, rcdi, parities, dplcts, hamOff_dp, hamDi_dp, nOff, nDi_off)
        !             else if(id == 2) then ! Hopping for irreducible representations in 2D (currently not supported)
        !                 call hopping(thrds, tilted, nHel, tilt, lx, ly, sites, nbb, dim, basis, orbsize, orbits, phases, norm2D, trisites, xtransl, ytransl, symm, id, par, rot, refl, c6, t, rcoff, rcdi, parities, dplcts, hamOff_dc, hamDi_dc, nOff, nDi_off)
        !                 call diagonal(sites, dim, trisites, hexsites, orbsize, orbits, phases, norm2D, hamDi_c2, occ)
        !             end if 
        !         else if(ti == 1 .and. ((k1 .ne. 0) .or. (k2 .ne. 0))) then ! If translation symmetry is used and k1 or k2 are not zero, generate the Hamiltonian with the hopping procedure for complex entries. Currently only supported for 1D irreps. 
        !             call hopping(unit, prms, tilted, nHel, tilt, lx, ly, ucx, ucy, sites, nbb, dim, basis, trisites, norm, xtransl, ytransl, t, k1, k2, rcoff, hamOff_dc, nOff)
        !         end if 
        !     end if
            
        !     !Generates diagonal part of the Hamiltonian if it does not exist yet.
        !     if(id == 1 .and. (.not.(exist2) .or. .not.(exist3))) call diagonal(unit, prms, sites, nbb, nnnbb, dim, trisites, hexsites, basis, hamDi, occ, nDi) 
            
        !     ! Saves the Hamiltonian to file if it does not exist yet.
        !     if(id == 1 .and. (.not.(exist1) .or. .not.(exist2) .or. .not.(exist3))) call save_ham(type, prms, unit, ti, sites, nOff, nDi, hamDi, hamOff, rcoff, hamOff_dp, hamOff_dc, occ)

        !     return 

        ! end subroutine generate_hamiltonian

        ! subroutine hopping_i(flag)
        !     ! subroutine hopping(othrds, unit, parameters, sites, nnBonds, bsites, dim, basis_list, hamOff, nnz)

        !     !-------------------------------------------------------------!
        !     !          Hopping procedure for real matrix entries          !
        !     !-------------------------------------------------------------!
        !     ! This subroutine generates the hopping Hamiltonian for a given basis and bond structure.
        !     use variables
        !     use input_variables
        !     implicit none
            
        !     integer, intent(in) :: flag ! Integer flag to indicate this routine 
        !     ! integer, intent(in)                       :: othrds, unit, sites, nnBonds, bsites(2, nnBonds) 
        !     ! integer(kind=8), intent(in)               :: dim, basis_list(dim)
        !     ! character(len=*), intent(in)              :: parameters
            
        !     ! integer(kind=8), intent(out)              :: nnz
        !     ! integer(kind=8), allocatable, intent(out) :: ham(:,:)

        !     ! Local variables
        !     integer                      :: parity1 = 0, parity2 = 0
        !     integer(kind=8)              :: i = 0, j = 0, s = 0, arrsize = 0, cntr2 = 0, loc = 0, pos = 0, cntrj = 0, newst = 0  
        !     integer(kind=8), allocatable :: cntr(:), ham_temp(:,:)
                
            
        !     arrsize = 2*nnBonds*dim ! Maximum number of non-zero elements in the Hamiltonian matrix. Each basis state can have at most 2*nnBonds non-zero elements (one for each bond).
            
        !     if(allocated(ham_temp)) deallocate(ham_temp) ! Temporary array to store non-zero elements of the Hamiltonian matrix before number of non-zero elements (nnz) is known.
        !     if(allocated(cntr))     deallocate(cntr)
        !     allocate(ham_temp(arrsize, 3))
        !     allocate(cntr(othrds))
            
        !     ham_temp = 0 
        !     cntr     = 0
        !     nnz      = 0

        !     !$omp parallel do default(firstprivate) shared(ham_temp, cntr, bsites, basis_states) num_threads(othrds)
        !     do j = 1, dim
        !         !$ id = omp_get_thread_num()
        !         cntrj = 0 
        !         do s = 1, nnBonds    

        !             if(btest(basis_states(j), bsites(1, s) - 1) .and. .not. btest(basis_states(j), bsites(2, s) - 1)) then 
        !                 newst = ibclr(ibset(basis_states(j), bsites(2, s) - 1), bsites(1, s) - 1) !Create on site 2, annihilate on site 1 
        !                 call findstate(dim, newst, basis_states, loc) 
        !                 if(loc > 0) then 
        !                     cntr(id+1)      = cntr(id+1) + 1 !Count non-zero elements in each thread 
        !                     cntrj           = cntrj + 1 !Count non-zero elements per basis state 
        !                     pos             = (j-1)*nnBonds + cntrj !Position of non-zero element in sparse array: Offset ((j-1)*nnBonds): Each basis state (j) has at most nnBonds non-zero elements. Cntrj gives the current increment.
        !                     parity1         = popcnt(ibits(basis_states(j), bsites(1, s), sites)) !Parity of site 1
        !                     parity2         = popcnt(ibits(ibclr(basis_states(j), bsites(1, s) - 1), bsites(2, s), sites)) !Parity of site 2
        !                     ham_temp(pos,1) = (-1)**(parity1 + parity2) 
        !                     ham_temp(pos,2) = j   !Row index of non-zero element in sparse matrix
        !                     ham_temp(pos,3) = loc !Column index of non-zero element in sparse matrix
        !                 end if 
        !             else if(btest(basis_states(j), bsites(2, s) - 1) .and. .not. btest(basis_states(j), bsites(1, s) - 1)) then 
        !                 newst = ibclr(ibset(basis_states(j), bsites(1, s) - 1), bsites(2, s) - 1) !Create on site 2, annihilate on site 1 
        !                 call findstate(dim, newst, basis_states, loc) 
        !                 if(loc > 0) then 
        !                     cntr(id+1)      = cntr(id+1) + 1
        !                     cntrj           = cntrj + 1
        !                     pos             = (j-1)*nnBonds + cntrj 
        !                     parity1         = popcnt(ibits(basis_states(j), bsites(2, s), sites)) !Parity of site 2
        !                     parity2         = popcnt(ibits(ibclr(basis_states(j), bsites(2,s) - 1), bsites(1, s), sites)) !Parity of site 1
        !                     ham_temp(pos,1) = (-1)**(parity1 + parity2)  
        !                     ham_temp(pos,2) = j
        !                     ham_temp(pos,3) = loc
        !                 end if                    
        !             end if 
                
        !         end do

        !     end do
        !     !$omp end parallel do

        !     nnz = sum(cntr)

        !     if(allocated(ham)) deallocate(ham) !Allocate the final Hamiltonian matrix with the correct number of non-zero elements.
        !     allocate(ham(nnz,3))
        !     ham   = 0 
        !     cntr2 = 0 
        !     do i = 1, arrsize
        !         if(ham_temp(i, 1) .ne. 0) then 
        !             cntr2        = cntr2 + 1 
        !             ham(cntr2,1) = ham_temp(i,1)
        !             ham(cntr2,2) = ham_temp(i,2)
        !             ham(cntr2,3) = ham_temp(i,3)    
        !         end if 
        !     end do
        !     if(cntr2 .ne. nnz) stop 'Subroutine hopping: Counter > NNZ'
        !     if(allocated(ham_temp)) deallocate(ham_temp)

        !     print*, 'Generated hopping Hamiltonian.'
        !     print*, 'Number of non-zero matrix elements: ', nnz
        !     print*, ''

    ! end subroutine hopping_i

    !With variable arguments
    subroutine hopping_i(threads, unit, sites, nbonds, bsites, dim, basis_list, ham, nnz)

        !-------------------------------------------------------------!
        !          Hopping procedure for real matrix entries          !
        !-------------------------------------------------------------!
        ! This subroutine generates the hopping Hamiltonian for a given basis and bond structure.
        
        implicit none

        integer, intent(in)                       :: threads, unit, sites, nbonds, bsites(2, nbonds) 
        integer(kind=8), intent(in)               :: dim, basis_list(dim)
                
        integer(kind=8), intent(out)              :: nnz
        integer(kind=8), allocatable, intent(out) :: ham(:,:)

        ! Local variables
        integer                      :: id = 0, parity1 = 0, parity2 = 0
        integer(kind=8)              :: i = 0, j = 0, s = 0, arrsize = 0, cntr2 = 0, loc = 0, pos = 0, cntrj = 0, newst = 0  
        integer(kind=8), allocatable :: cntr(:), ham_temp(:,:)
            
        
        arrsize = 2*nbonds*dim ! Maximum number of non-zero elements in the Hamiltonian matrix. Each basis state can have at most 2*nbonds non-zero elements (one for each bond).
        
        if(allocated(ham_temp)) deallocate(ham_temp) ! Temporary array to store non-zero elements of the Hamiltonian matrix before number of non-zero elements (nnz) is known.
        if(allocated(cntr))     deallocate(cntr)
        allocate(ham_temp(arrsize, 3))
        allocate(cntr(threads))
        
        ham_temp = 0 
        cntr     = 0
        nnz      = 0

        !$omp parallel do default(firstprivate) shared(ham_temp, cntr, bsites, basis_list) num_threads(threads)
        do j = 1, dim
            !$ id = omp_get_thread_num()
            cntrj = 0 
            do s = 1, nbonds    

                if(btest(basis_list(j), bsites(1, s) - 1) .and. .not. btest(basis_list(j), bsites(2, s) - 1)) then 
                    newst = ibclr(ibset(basis_list(j), bsites(2, s) - 1), bsites(1, s) - 1) !Create on site 2, annihilate on site 1 
                    call findstate(dim, newst, basis_list, loc) 
                    if(loc > 0) then 
                        cntr(id+1)      = cntr(id+1) + 1 !Count non-zero elements in each thread 
                        cntrj           = cntrj + 1 !Count non-zero elements per basis state 
                        pos             = (j-1)*nbonds + cntrj !Position of non-zero element in sparse array: Offset ((j-1)*nbonds): Each basis state (j) has at most nbonds non-zero elements. Cntrj gives the current increment.
                        parity1         = popcnt(ibits(basis_list(j), bsites(1, s), sites)) !Parity of site 1
                        parity2         = popcnt(ibits(ibclr(basis_list(j), bsites(1, s) - 1), bsites(2, s), sites)) !Parity of site 2
                        ham_temp(pos,1) = (-1)**(parity1 + parity2) 
                        ham_temp(pos,2) = j   !Row index of non-zero element in sparse matrix
                        ham_temp(pos,3) = loc !Column index of non-zero element in sparse matrix
                    end if 
                else if(btest(basis_list(j), bsites(2, s) - 1) .and. .not. btest(basis_list(j), bsites(1, s) - 1)) then 
                    newst = ibclr(ibset(basis_list(j), bsites(1, s) - 1), bsites(2, s) - 1) !Create on site 2, annihilate on site 1 
                    call findstate(dim, newst, basis_list, loc) 
                    if(loc > 0) then 
                        cntr(id+1)      = cntr(id+1) + 1
                        cntrj           = cntrj + 1
                        pos             = (j-1)*nbonds + cntrj 
                        parity1         = popcnt(ibits(basis_list(j), bsites(2, s), sites)) !Parity of site 2
                        parity2         = popcnt(ibits(ibclr(basis_list(j), bsites(2,s) - 1), bsites(1, s), sites)) !Parity of site 1
                        ham_temp(pos,1) = (-1)**(parity1 + parity2)  
                        ham_temp(pos,2) = j
                        ham_temp(pos,3) = loc
                    end if                    
                end if 
            
            end do

        end do
        !$omp end parallel do

        nnz = sum(cntr)

        if(allocated(ham)) deallocate(ham) !Allocate the final Hamiltonian matrix with the correct number of non-zero elements.
        allocate(ham(nnz,3))
        ham   = 0 
        cntr2 = 0 
        do i = 1, arrsize
            if(ham_temp(i, 1) .ne. 0) then 
                cntr2        = cntr2 + 1 
                ham(cntr2,1) = ham_temp(i,1)
                ham(cntr2,2) = ham_temp(i,2)
                ham(cntr2,3) = ham_temp(i,3)    
            end if 
        end do
        if(cntr2 .ne. nnz) stop 'Subroutine hopping: Counter > NNZ'
        if(allocated(ham_temp)) deallocate(ham_temp)

        print*, 'Generated hopping Hamiltonian.'
        print*, 'Number of non-zero matrix elements: ', nnz
        print*, ''

    end subroutine hopping_i

    subroutine hopping_dc(unit, gencluster, nHel, tilt, lx, ly, ucx, ucy, sites, nbonds, dim, basis_list, bsites, norm, xtransl, ytransl, t, k1, k2, rc, ham, nnz)
       
        !----------------------------------------------------------------!
        !          Hopping procedure for complex matrix entries          !
        !----------------------------------------------------------------!
        ! This subroutine generates the hopping Hamiltonian for a given basis and bond structure in momentum sector (k1, k2).
       
        use parameters

        implicit none

        integer, intent(in) :: unit, gencluster, nHel, tilt, lx, ly, ucx, ucy, sites, nbonds, k1, k2, bsites(2, nbonds), xtransl(2, sites), ytransl(2, sites)
        integer(kind=8), intent(in) :: dim, basis_list(dim)
        double precision, intent(in) :: t, norm(dim)

        integer(kind=8), intent(out) :: nnz
        integer(kind=8), allocatable, intent(out) :: rc(:,:)
        double complex, allocatable, intent(out) :: ham(:)
        
        ! Local variables
        integer :: dbl = 0, parity1 = 0, parity2 = 0, l1 = 0, l2 = 0, sign = 0, l11 = 0, l22 = 0 
        integer(kind=8) :: arrsize = 0, i = 0, j = 0, k = 0, l = 0, s = 0, loc = 0, rep = 0, mask = 0, newst = 0, n_temp = 0, cntr = 0, cntrj = 0
        integer(kind=8), allocatable :: rc_temp(:,:)
        double complex, allocatable :: ham_temp(:)

        arrsize = 6*sites*dim! Maximum number of non-zero elements in the Hamiltonian matrix. Each basis state can have at most 2*nbonds non-zero elements (one for each bond).
        
        if(allocated(ham_temp)) deallocate(ham_temp)
        if(allocated(rc_temp))  deallocate(rc_temp)
        allocate(ham_temp(arrsize))
        allocate(rc_temp(arrsize, 2))
        
        ham_temp = 0.d0  
        rc_temp  = 0 
        n_temp   = 0 
        nnz      = 0 

        do j = 1, dim
            cntrj = n_temp      !Starting point in temp_rc for matrix elements of row 'j' to look for double entries with the same column 'loc'
            cntr  = 0           !Counts the number of already calculated matrix elements of each row 'j'
            do s = 1, nbonds 
                ! !Create scattered basis state I'
                    ! mask = ibset(ibset(0,bsites(1,s)-1),bsites(2,s)-1) !Procedure suggested in Lin-paper. Sets a mask with only non-zero components on sites 1 and 2
                    ! k    = iand(mask, basis_list(j)) !K records the occupancy of those two sites.
                    ! l    = ieor(mask,k) !L records whether hopping is allowed or not. If it's 0 or the same as the mask, hopping is not allowed. If it is allowed the occupations of both sites (01 or 10) are swapped (10 or 01)
                    ! if(l == 0 .or. l == mask) then
                    !     cycle
                    ! end if
                    ! call representative(basis_list(j) - k + l, sites, ucx, ucy, xtransl, ytransl, rep, l1, l2, sign) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            

                if(btest(basis_list(j), bsites(1,s) - 1) .and. .not. (btest(basis_list(j), bsites(2,s) - 1))) then 
                    newst = ibclr(ibset(basis_list(j), bsites(2, s) - 1), bsites(1, s) - 1) 
                else if(btest(basis_list(j), bsites(2,s) - 1) .and. .not. (btest(basis_list(j), bsites(1,s) - 1))) then 
                    newst = ibclr(ibset(basis_list(j), bsites(1, s) - 1), bsites(2, s) - 1) 
                else 
                    cycle 
                end if 
                if(gencluster == 1) then 
                    call representative(newst, sites, ucx, ucy, xtransl, ytransl, rep, l1, l2, sign) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                    l11 = ucx 
                    l22 = ucy
                else if(gencluster == 0) then 
                    call representative_tilted(newst, sites, nHel, tilt, lx, ly, xtransl, ytransl, rep, l1, l2, sign)
                    l11 = ly 
                    l22 = lx
                end if 
                call findstate(dim, rep, basis_list, loc) !Finds the location of representative in basis_list
                if(loc > 0) then
                    if(btest(basis_list(j),bsites(1,s)-1) .and. .not.(btest(basis_list(j),bsites(2,s)-1))) then !Hopping from 1->2
                        parity1 = popcnt(ibits(basis_list(j), bsites(1,s), sites))
                        parity2 = popcnt(ibits(ibclr(basis_list(j), bsites(1,s) - 1), bsites(2,s), sites)) 
                    else if(btest(basis_list(j),bsites(2,s)-1) .and. .not.(btest(basis_list(j),bsites(1,s)-1))) then !Hopping from 2->1
                        parity1 = popcnt(ibits(ibclr(basis_list(j), bsites(2,s) - 1), bsites(1,s), sites ))
                        parity2 = popcnt(ibits(basis_list(j), bsites(2,s), sites)) 
                    end if 
                    dbl  = 0 !Flag for whether a new matrix element should be created (0) or an old matrix element was updated (1)
                    cntr = cntr + 1
                    if(cntr > 1) then
                        do i = max(cntrj,1), cntrj + cntr !Loop for updating existing elements: Goes through all matrix elements already calculated for row 'j' and checks whether column 'loc' already exists.
                            if(rc_temp(i,2) == loc .and. rc_temp(i,1) == j) then !If yes, update existing element.                            
                                ham_temp(i) = ham_temp(i) + sign * (-1) * t * (-1)**(parity1 + parity2) * sqrt(real(norm(loc))/real(norm(j))) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))       
                                dbl = 1 !Found and updated existing element.
                            end if
                        end do
                    end if
                    if(dbl == 0) then !If no existing element was found, create new matrix element.
                        n_temp = n_temp + 1 !Counter for non-zero matrix elements.
                        ham_temp(n_temp) = ham_temp(n_temp) + sign * (-1) * t * (-1)**(parity1 + parity2) * sqrt(real(norm(loc))/real(norm(j))) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                        ! ham_temp(n_temp) = ham_temp(n_temp) + (-1) * t * (-1)**(parity1 + parity2) * exp(-ii*2*pi*(k1*l1 + k2*l2)/sites)
                        ! ham_temp(n_temp) = ham_temp(n_temp) + (-1) * t * (-1)**(parity1 + parity2) * sqrt(real(norm(loc))/real(norm(j))) * exp(-ii*2*pi*(k1*dble(l1-Lx)/2.d0 + k2*dble(l2-Ly)/2.d0)/sites)
                        rc_temp(n_temp,1) = j !Row of matrix element
                        rc_temp(n_temp,2) = loc !Column of matrix element
                    end if
                end if
            end do 
        end do 
            
        nnz = n_temp 
        if(allocated(ham)) deallocate(ham)
        if(allocated(rc))  deallocate(rc)
        allocate(ham(nnz))
        allocate(rc(nnz,2))
        ham = 0.d0  
        rc  = 0

        do i = 1, nnz
            ham(i)  = ham_temp(i)
            rc(i,1) = rc_temp(i,1)
            rc(i,2) = rc_temp(i,2)
        end do

        if(allocated(ham_temp)) deallocate(ham_temp)
        if(allocated(rc_temp))  deallocate(rc_temp)

        print*,'Generated hopping Hamiltonian.'
        print*,'Number of non-zero matrix elements: ', nnz
        print*,''

    end subroutine hopping_dc

    
    ! subroutine hopping_irrep(flag)

        !     !---------------------------------------------------------------------!
        !     !          Hopping procedure for irreducible representations          !
        !     !---------------------------------------------------------------------!
        !     ! This subroutine generates the hopping Hamiltonian for irreducible representations in 1D.
        !     ! It uses the point group operators and translations to generate the Hamiltonian matrix.
            
        !     implicit none

        !     double precision, intent(in) :: flag ! Double precision flag to indicate this routine
        !     ! integer(kind=8), intent(in)  :: dim, basis(dim)
        !     ! integer, intent(in)          :: threads, tilted, nHel, tilt, lx, ly, ucx, ucy, sites, nbonds, symm, bsites(2, nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
        !     ! double precision, intent(in) :: t, id, par(6), rot(5), norm(dim)
            
        !     ! integer(kind=8), intent(out)               :: nnz, nDi 
        !     ! integer        , allocatable, intent(out)  :: parities(:), dplcts(:)
        !     ! integer(kind=8), allocatable, intent(out)  :: rc(:,:), rcdi(:)
        !     ! double precision, allocatable, intent(out) :: ham(:), hamDi(:)
            
        !     !Local variables
        !     integer                       :: signrep = 1, dbl = 0, order, info, s = 0, k = 0, l1 = 0, l2 = 0, parity = 0
        !     integer(kind=8)               :: i = 0, j = 0, loc = 0, rep = 0, pos = 0, rowst = 0, newst = 0, state = 0, cntr = 0, cntrj = 0, n_temp = 0, arrsize = 0
        !     integer(kind=8), allocatable  :: rc_temp(:,:), rcdi_temp(:), parities_temp(:), dplcts_temp(:), cntr_di(:,:)
        !     double precision              :: sign = 1.d0, signstate = 1.d0, h_add = 0.d0 
        !     double precision, allocatable :: ham_temp(:), hamDi_temp(:)
            
        !     if(symm == 1) then 
        !         order = 12
        !     else 
        !         order = 1
        !     end if 

        !     arrsize = dim * nbonds * order ! Maximum number of non-zero elements in the Hamiltonian matrix. Each basis state can have at most nbonds*order non-zero elements (one for each bond and point group operator).
        !     if(allocated(ham_temp))    deallocate(ham_temp)
        !     if(allocated(rc_temp))     deallocate(rc_temp)
        !     if(allocated(hamDi_temp))  deallocate(hamDi_temp)
        !     if(allocated(rcdi_temp))   deallocate(rcdi_temp)
        !     if(allocated(parities_temp))   deallocate(parities_temp)
        !     if(allocated(dplcts_temp)) deallocate(dplcts_temp)
        !     if(allocated(cntr_di))     deallocate(cntr_di)
            
        !     allocate(ham_temp(arrsize))
        !     allocate(rc_temp(arrsize, 2))
        !     allocate(hamDi_temp(arrsize))
        !     allocate(rcdi_temp(arrsize))
        !     allocate(parities_temp(arrsize))
        !     allocate(dplcts_temp(arrsize))
        !     allocate(cntr_di(dim, 2))
        !     h_add       = 0.d0  
        !     ham_temp    = 0.d0  
        !     hamDi_temp  = 0.d0  
        !     rc_temp     = 0 
        !     rcdi_temp   = 0 
        !     parities_temp   = 0 
        !     dplcts_temp = 0 
        !     n_temp      = 0 
        !     nDi         = 0 
        !     nnz         = 0 
        !     cntr_di     = 0 

        !     !$omp parallel do default(firstprivate) shared(ham_temp, bsites, basis) num_threads(threads)
        !     do j = 1, dim
        !         cntrj = 0 !Counts the number of already calculated matrix elements of each row 'j'
        !         rowst = (j-1) * nbonds * order 
        !         !$ id = omp_get_thread_num()
        !         do k = 1, order !Loop over point groups operators 
        !             h_add     = 0.d0 !Scattering contribution to matrix element
        !             signstate = 1 !Sign due to symmetry transformations and translations
        !             if(k == 1) then !Prepare initial state 
        !                 state = basis(j)
        !             else if(k <= 7) then 
        !                 call reflect(basis(j), basis(j), sites, refl(k - 1, 1:sites), signstate, info, state)
        !                 signstate = signstate * par(k - 1)
        !             else if(8 <= k) then 
        !                 call c6n(basis(j), basis(j), sites, k - 7, c6, signstate, info, state)
        !                 signstate = signstate * rot(k - 7)
        !             end if 

        !             do s = 1, nbonds !Takes care of translations
        !                 if(btest(state, bsites(1,s)-1) .and. .not.(btest(state, bsites(2,s)-1))) then !Hopping from 1->2
        !                     newst  = ibclr(ibset(state, bsites(2, s) - 1), bsites(1, s) - 1)            
        !                     parity = popcnt(ibits(state, bsites(1,s), sites)) + popcnt(ibits(ibclr(state, bsites(1,s) - 1), bsites(2,s), sites))      
        !                 else if(btest(state, bsites(2,s)-1) .and. .not.(btest(state, bsites(1,s)-1))) then !Hopping from 2->1
        !                     newst  = ibclr(ibset(state, bsites(1, s) - 1), bsites(2, s) - 1) 
        !                     parity = popcnt(ibits(state, bsites(2,s), sites)) + popcnt(ibits(ibclr(state, bsites(2,s) - 1), bsites(1,s), sites))
        !                 else 
        !                     cycle 
        !                 end if 
        !                 if(tilted == 0) then 
        !                     call representative_irrep_rect(newst, sites, Lx, Ly, symm, id, par, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
        !                 else if(tilted == 1) then 
        !                     call representative_irrep(newst, sites, nHel, tilt, Lx, Ly, symm, id, par, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
        !                 end if 
                        
        !                 call findstate(dim, rep, basis, loc) !Finds the location of representative in basis
        !                 if(loc <= 0) cycle 
        !                 sign  = signrep * signstate * (-1) * (-1)**parity 
        !                 h_add = sign * t * sqrt(dble(norm(loc))/dble(norm(j)))  

        !                 if(loc == j) then !Diagonal element contributes to diagonal Hamiltonian
        !                     if(cntr_di(j, 1) == 1) then !Not the first diagonal element of state j. Add to previous element inDicated by cntr_di(j, 2). 
        !                         hamDi_temp(cntr_di(j, 2)) = hamDi_temp(cntr_di(j, 2)) + h_add
        !                     else !First diagonal element of state j. Create new element and save position in cntr_di(j, 2).
        !                         nDi             = nDi + 1
        !                         hamDi_temp(nDi) = hamDi_temp(nDi) + h_add
        !                         rcdi_temp(nDi)  = j 
        !                         cntr_di(j, 1)   = 1
        !                         cntr_di(j, 2)   = nDi
        !                     end if 
        !                 else !Off-diagonal elements 
        !                     dbl = 0 !Flag for whether a new matrix element should be created (0) or an old matrix element was updated (1)
        !                     if(cntrj > 0) then !Not the first element in row j. Check elements for double entries. 
        !                         do i = rowst + 1, rowst + cntrj + 1 !Go from beginning of row j to current element. 
        !                             if(rc_temp(i,2) == loc .and. rc_temp(i,1) == j) then !If yes, update existing element.
        !                                 if(i > arrsize) stop 'hopping irrep: i > arrsize'
        !                                 ham_temp(i)    = ham_temp(i) + h_add 
        !                                 parities_temp(i)   = (-1)**parity
        !                                 dplcts_temp(i) = dplcts_temp(i) + 1 
        !                                 dbl            = 1 !Found and updated existing element.       
        !                             end if
        !                         end do
        !                     end if
        !                     if(dbl == 0) then !If no existing element was found, create new matrix element.
        !                         if(pos > arrsize) stop 'hopping irrep: pos > arrsize'
        !                         ! n_temp(id) = n_temp(id) + 1 !Counter for non-zero matrix elements.       
        !                         cntrj            = cntrj + 1 !One more element found for row j. 
        !                         n_temp           = n_temp + 1 !Counter for non-zero matrix elements. Needed later.
        !                         pos              = rowst + cntrj !Position corresponds to position in row j. 
        !                         ham_temp(pos)    = ham_temp(pos) + h_add !Value of matrix element. 
        !                         parities_temp(pos)   = (-1)**parity !Save parities for later. 
        !                         dplcts_temp(pos) = dplcts_temp(pos) + 1 !Save number of dplcts for later. 
        !                         rc_temp(pos,1)   = j !Row of matrix element. 
        !                         rc_temp(pos,2)   = loc !Column of matrix element.                          
        !                     end if
        !                 end if 
        !             end do 

        !         end do !point group operations
        !     end do 

        !     nnz = n_temp

        !     if(allocated(rc))     deallocate(rc)
        !     if(allocated(ham))    deallocate(ham)
        !     if(allocated(rcdi))   deallocate(rcdi)
        !     if(allocated(hamDi))  deallocate(hamDi)
        !     if(allocated(parities))   deallocate(parities)
        !     if(allocated(dplcts)) deallocate(dplcts)

        !     allocate(ham(nnz))
        !     allocate(rc(nnz,2))
        !     allocate(rcdi(nDi))
        !     allocate(hamDi(nDi))
        !     allocate(parities(nnz))
        !     allocate(dplcts(nnz))

        !     rc    = 0
        !     cntr  = 0
        !     rcdi  = 0
        !     ham   = 0.d0  
        !     hamDi = 0.d0  

        !     do i = 1, arrsize
        !         if(rc_temp(i,1) > 0) then 
        !             cntr         = cntr + 1
        !             ham(cntr)    = ham_temp(i)
        !             rc(cntr,1)   = rc_temp(i,1)
        !             rc(cntr,2)   = rc_temp(i,2)
        !             parities(cntr)   = parities_temp(i)
        !             dplcts(cntr) = dplcts_temp(i)
        !         end if 
        !     end do
        !     ham = ham / dble(order)

        !     do i = 1, nDi !Fill up diagonal Hamiltonian
        !         hamDi(i) = hamDi_temp(i)
        !         rcdi(i)  = rcdi_temp(i)
        !     end do
        !     hamDi = hamDi / order 

        !     print*, 'Generated hopping Hamiltonian.'
        !     print*, 'Number of non-zero matrix elements: ', nnz
        !     print*, 'Number of diagonal matrix elements: ', nDi
        !     print*, 'Number of representatives: ', dim
        !     print*, ''

    ! end subroutine hopping_irrep
    
    !With variable arguments
    subroutine hopping_irrep(threads, tilted, nHel, tilt, lx, ly, ucx, ucy, sites, nbonds, dim, basis, bsites, norm, xtransl, ytransl, symm, id, par, rot, refl, c6, t, rc, rcdi, parities, dplcts, ham, hamDi, nnz, nDi)

        !---------------------------------------------------------------------!
        !          Hopping procedure for irreducible representations          !
        !---------------------------------------------------------------------!
        ! This subroutine generates the hopping Hamiltonian for irreducible representations in 1D.
        ! It uses the point group operators and translations to generate the Hamiltonian matrix.
        
        implicit none
        integer(kind=8), intent(in)  :: dim, basis(dim)
        integer, intent(in)          :: threads, tilted, nHel, tilt, lx, ly, ucx, ucy, sites, nbonds, symm, bsites(2, nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
        double precision, intent(in) :: t, id, par(6), rot(5), norm(dim)
        
        integer(kind=8), intent(out)               :: nnz, nDi 
        integer        , allocatable, intent(out)  :: parities(:), dplcts(:)
        integer(kind=8), allocatable, intent(out)  :: rc(:,:), rcdi(:)
        double precision, allocatable, intent(out) :: ham(:), hamDi(:)
        
        !Local variables
        integer                       :: dbl = 0, order, info, s = 0, k = 0, l1 = 0, l2 = 0, parity = 0
        integer(kind=8)               :: i = 0, j = 0, loc = 0, rep = 0, pos = 0, rowst = 0, newst = 0, state = 0, cntr = 0, cntrj = 0, n_temp = 0, arrsize = 0
        integer(kind=8), allocatable  :: rc_temp(:,:), rcdi_temp(:), cntr_di(:,:)
        integer        , allocatable  :: parities_temp(:), dplcts_temp(:)
        double precision              :: sign = 1.d0, signstate = 1.d0, h_add = 0.d0, signrep = 1.d0 
        double precision, allocatable :: ham_temp(:), hamDi_temp(:)
        
        if(symm == 1) then 
            order = 12
        else 
            order = 1
        end if 

        arrsize = dim * nbonds * order ! Maximum number of non-zero elements in the Hamiltonian matrix. Each basis state can have at most nbonds*order non-zero elements (one for each bond and point group operator).
        if(allocated(ham_temp))    deallocate(ham_temp)
        if(allocated(rc_temp))     deallocate(rc_temp)
        if(allocated(hamDi_temp))  deallocate(hamDi_temp)
        if(allocated(rcdi_temp))   deallocate(rcdi_temp)
        if(allocated(parities_temp))   deallocate(parities_temp)
        if(allocated(dplcts_temp)) deallocate(dplcts_temp)
        if(allocated(cntr_di))     deallocate(cntr_di)
        
        allocate(ham_temp(arrsize))
        allocate(rc_temp(arrsize, 2))
        allocate(hamDi_temp(arrsize))
        allocate(rcdi_temp(arrsize))
        allocate(parities_temp(arrsize))
        allocate(dplcts_temp(arrsize))
        allocate(cntr_di(dim, 2))
        h_add       = 0.d0  
        ham_temp    = 0.d0  
        hamDi_temp  = 0.d0  
        rc_temp     = 0 
        rcdi_temp   = 0 
        parities_temp   = 0 
        dplcts_temp = 0 
        n_temp      = 0 
        nDi         = 0 
        nnz         = 0 
        cntr_di     = 0 

        !$omp parallel do default(firstprivate) shared(ham_temp, bsites, basis) num_threads(threads)
        do j = 1, dim
            cntrj = 0 !Counts the number of already calculated matrix elements of each row 'j'
            rowst = (j-1) * nbonds * order 
            !$ id = omp_get_thread_num()
            do k = 1, order !Loop over point groups operators 
                h_add     = 0.d0 !Scattering contribution to matrix element
                signstate = 1 !Sign due to symmetry transformations and translations
                if(k == 1) then !Prepare initial state 
                    state = basis(j)
                else if(k <= 7) then 
                    call reflect(basis(j), basis(j), sites, refl(k - 1, 1:sites), signstate, info, state)
                    signstate = signstate * par(k - 1)
                else if(8 <= k) then 
                    call c6n(basis(j), basis(j), sites, k - 7, c6, signstate, info, state)
                    signstate = signstate * rot(k - 7)
                end if 

                do s = 1, nbonds !Takes care of translations
                    if(btest(state, bsites(1,s)-1) .and. .not.(btest(state, bsites(2,s)-1))) then !Hopping from 1->2
                        newst  = ibclr(ibset(state, bsites(2, s) - 1), bsites(1, s) - 1)            
                        parity = popcnt(ibits(state, bsites(1,s), sites)) + popcnt(ibits(ibclr(state, bsites(1,s) - 1), bsites(2,s), sites))      
                    else if(btest(state, bsites(2,s)-1) .and. .not.(btest(state, bsites(1,s)-1))) then !Hopping from 2->1
                        newst  = ibclr(ibset(state, bsites(1, s) - 1), bsites(2, s) - 1) 
                        parity = popcnt(ibits(state, bsites(2,s), sites)) + popcnt(ibits(ibclr(state, bsites(2,s) - 1), bsites(1,s), sites))
                    else 
                        cycle 
                    end if 
                    if(tilted == 0) then 
                        call representative_irrep_rect(newst, sites, Lx, Ly, symm, id, par, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                    else if(tilted == 1) then 
                        call representative_irrep(newst, sites, nHel, tilt, Lx, Ly, symm, id, par, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                    end if 
                    
                    call findstate(dim, rep, basis, loc) !Finds the location of representative in basis
                    if(loc <= 0) cycle 
                    sign  = signrep * signstate * (-1) * (-1)**parity 
                    h_add = sign * t * sqrt(dble(norm(loc))/dble(norm(j)))  

                    if(loc == j) then !Diagonal element contributes to diagonal Hamiltonian
                        if(cntr_di(j, 1) == 1) then !Not the first diagonal element of state j. Add to previous element inDicated by cntr_di(j, 2). 
                            hamDi_temp(cntr_di(j, 2)) = hamDi_temp(cntr_di(j, 2)) + h_add
                        else !First diagonal element of state j. Create new element and save position in cntr_di(j, 2).
                            nDi             = nDi + 1
                            hamDi_temp(nDi) = hamDi_temp(nDi) + h_add
                            rcdi_temp(nDi)  = j 
                            cntr_di(j, 1)   = 1
                            cntr_di(j, 2)   = nDi
                        end if 
                    else !Off-diagonal elements 
                        dbl = 0 !Flag for whether a new matrix element should be created (0) or an old matrix element was updated (1)
                        if(cntrj > 0) then !Not the first element in row j. Check elements for double entries. 
                            do i = rowst + 1, rowst + cntrj + 1 !Go from beginning of row j to current element. 
                                if(rc_temp(i,2) == loc .and. rc_temp(i,1) == j) then !If yes, update existing element.
                                    if(i > arrsize) stop 'hopping irrep: i > arrsize'
                                    ham_temp(i)    = ham_temp(i) + h_add 
                                    parities_temp(i)   = (-1)**parity
                                    dplcts_temp(i) = dplcts_temp(i) + 1 
                                    dbl            = 1 !Found and updated existing element.       
                                end if
                            end do
                        end if
                        if(dbl == 0) then !If no existing element was found, create new matrix element.
                            if(pos > arrsize) stop 'hopping irrep: pos > arrsize'
                            ! n_temp(id) = n_temp(id) + 1 !Counter for non-zero matrix elements.       
                            cntrj            = cntrj + 1 !One more element found for row j. 
                            n_temp           = n_temp + 1 !Counter for non-zero matrix elements. Needed later.
                            pos              = rowst + cntrj !Position corresponds to position in row j. 
                            ham_temp(pos)    = ham_temp(pos) + h_add !Value of matrix element. 
                            parities_temp(pos)   = (-1)**parity !Save parities for later. 
                            dplcts_temp(pos) = dplcts_temp(pos) + 1 !Save number of dplcts for later. 
                            rc_temp(pos,1)   = j !Row of matrix element. 
                            rc_temp(pos,2)   = loc !Column of matrix element.                          
                        end if
                    end if 
                end do 

            end do !point group operations
        end do 

        nnz = n_temp

        if(allocated(rc))     deallocate(rc)
        if(allocated(ham))    deallocate(ham)
        if(allocated(rcdi))   deallocate(rcdi)
        if(allocated(hamDi))  deallocate(hamDi)
        if(allocated(parities))   deallocate(parities)
        if(allocated(dplcts)) deallocate(dplcts)

        allocate(ham(nnz))
        allocate(rc(nnz,2))
        allocate(rcdi(nDi))
        allocate(hamDi(nDi))
        allocate(parities(nnz))
        allocate(dplcts(nnz))

        rc    = 0
        cntr  = 0
        rcdi  = 0
        ham   = 0.d0  
        hamDi = 0.d0  

        do i = 1, arrsize
            if(rc_temp(i,1) > 0) then 
                cntr         = cntr + 1
                ham(cntr)    = ham_temp(i)
                rc(cntr,1)   = rc_temp(i,1)
                rc(cntr,2)   = rc_temp(i,2)
                parities(cntr) = parities_temp(i)
                dplcts(cntr) = dplcts_temp(i)
            end if 
        end do
        ham = ham / dble(order)

        do i = 1, nDi !Fill up diagonal Hamiltonian
            hamDi(i) = hamDi_temp(i)
            rcdi(i)  = rcdi_temp(i)
        end do
        hamDi = hamDi / order 

        print*, 'Generated hopping Hamiltonian.'
        print*, 'Number of non-zero matrix elements: ', nnz
        print*, 'Number of diagonal matrix elements: ', nDi
        print*, 'Number of representatives: ', dim
        print*, ''

    end subroutine hopping_irrep

    ! subroutine hopping_irrep2D(flag1, flag2)
        
        !     !------------------------------------------------------------------------!
        !     !          Hopping procedure for 2D irreducible representations          !
        !     !------------------------------------------------------------------------!
        !     ! This subroutine generates the hopping Hamiltonian for irreducible representations in 2D.
        !     ! Work in progress.

        !     use parameters
        !     implicit none
        !     integer, intent(in) :: flag2 ! Integer precision flag to indicate 2D irrep
        !     double precision, intent(in) :: flag1 ! Double precision flag to indicate this routine
        !     ! integer, intent(in):: orbsize, threads, tilted, nHel, tilt, lx, ly, sites, nbonds, symm, bsites(2, nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
        !     ! integer(kind=8), intent(in) :: dim, basis(dim)
        !     ! integer(kind=8), allocatable, intent(in) :: orbits(:,:,:)
        !     ! double precision, intent(in) :: t, id, par(6), rot(5), norm(dim,2)
        !     ! double complex, allocatable, intent(in) :: phases(:,:,:)
            
        !     ! integer(kind=8), intent(out) :: nnz, nDi 
        !     ! integer(kind=8), allocatable, intent(out) :: rc(:,:), rcdi(:), parities(:), dplcts(:)
        !     ! double complex, allocatable, intent(out) :: ham(:), hamDi(:)
            
        !     !Local variables
        !     integer :: s, c, d, o, l1, l2, parity, sign = 1, signrep, dbl, order, site1, site2
        !     integer(kind=8) :: i, j, loc, rep, pos, rowst, rowIndx, colIndx, newst, state, cntr, cntrjc, n_temp, arrsize 
        !     integer(kind=8), allocatable :: rc_temp(:,:), rcdi_temp(:), parities_temp(:), dplcts_temp(:), cntr_di(:,:)
        !     double precision, parameter :: tol = 1.0e-14
        !     double complex :: h_add, coeff, newcf 
        !     double complex, allocatable :: ham_temp(:), hamDi_temp(:)

            
        !     if(symm == 1) then 
        !         order = 12
        !     else 
        !         order = 1
        !     end if 

        !     arrsize = dim * orbsize * 2 !Each representative (dim) has two basis states (2), each of which is a linear combination of at most 'orbsize' states. 
        !     if(allocated(ham_temp))    deallocate(ham_temp)
        !     if(allocated(rc_temp))     deallocate(rc_temp)
        !     if(allocated(hamDi_temp))  deallocate(hamDi_temp)
        !     if(allocated(rcdi_temp))   deallocate(rcdi_temp)
        !     if(allocated(parities_temp))   deallocate(parities_temp)
        !     if(allocated(dplcts_temp)) deallocate(dplcts_temp)
        !     if(allocated(cntr_di))     deallocate(cntr_di)
            
        !     allocate(dplcts_temp(arrsize))
        !     allocate(hamDi_temp(arrsize))
        !     allocate(rc_temp(arrsize, 2))
        !     allocate(rcdi_temp(arrsize))
        !     allocate(parities_temp(arrsize))
        !     allocate(ham_temp(arrsize)) 
        !     allocate(cntr_di(dim, 2))

        !     h_add       = 0.d0
        !     ham_temp    = 0.d0  
        !     hamDi_temp  = 0.d0  
        !     rc_temp     = 0 
        !     rcdi_temp   = 0 
        !     parities_temp   = 0 
        !     dplcts_temp = 0 
        !     n_temp      = 0 
        !     nDi         = 0 
        !     nnz         = 0 
        !     cntr_di     = 0 
        !     site1       = bsites(1,1)
        !     site2       = bsites(2,1)
            
        !     !$omp parallel do default(firstprivate) shared(ham_temp, bsites, basis) num_threads(threads)
        !     do j = 1, dim !All representatives 
        !         !$ id = omp_get_thread_num()
                
        !         do c = 1, 2 !First and second basis state of 2D irrep 
                
        !             h_add = 0.d0 !Scattering contribution to matrix element
        !             cntrjc = 0 !Counts the number of already calculated matrix elements of each combination row 'j' and basist state 'c', i.e. tuple (j, c) 
        !             rowIndx = 2*(j-1) + c
        !             rowst = (rowIndx-1) * orbsize !Position-1 of beginning of tuple (j,c)   
        !             do s = 1, orbsize !Loop through orbit of state basis(j)
        !                 state = orbits(j, s, c)
        !                 coeff = phases(j, s, c)
        !                 if(state == 0) exit 
        !                 if(abs(dble(coeff)) <= tol .and. abs(aimag(coeff)) <= tol) cycle 
        !                 if(btest(state, site1-1) .and. .not.(btest(state, site2-1))) then !Hopping from 1->2
        !                     newst  = ibclr(ibset(state, site2-1), site1-1)            
        !                     parity = popcnt(ibits(state, site1, sites)) + popcnt(ibits(ibclr(state, site1 - 1), site2, sites))      
        !                 else if(btest(state, site2-1) .and. .not.(btest(state, site1-1))) then !Hopping from 2->1
        !                     newst  = ibclr(ibset(state, site1-1), site2-1) 
        !                     parity = popcnt(ibits(state, site2, sites)) + popcnt(ibits(ibclr(state, site2-1), site1, sites))
        !                 else 
        !                     cycle 
        !                 end if 
        !                 if(tilted == 0) then 
        !                     call representative_irrep_rect(newst, sites, Lx, Ly, symm, id, par, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
        !                 else if(tilted == 1) then 
        !                     call representative_irrep(newst, sites, nHel, tilt, Lx, Ly, symm, id, par, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
        !                 end if 
        !                 call findstate(dim, rep, basis, loc) !Finds the location of representative in basis
        !                 if(loc <= 0) cycle !New state is compatible with momentum and symmetry 
        !                 do d = 1, 2 !Irrep basis states of scattered representative
        !                     colIndx = 2*(loc-1) + d
        !                     newcf   = 0.d0 
        !                     do o = 1, orbsize !Search for position of scattered state in orbit of Irrep's basis state and assign corresponDing coefficient 'newcf'   
        !                         if(orbits(loc, o, d) == newst) then 
        !                             newcf = phases(loc, o, d)
        !                             exit 
        !                         end if 
        !                     end do 
        !                     if(abs(dble(newcf)) <= tol .and. abs(aimag(newcf)) <= tol) cycle 

        !                     sign  = (-1) * (-1)**parity     
        !                     h_add = sign * t * sqrt(dble(norm(loc, d))/dble(norm(j, c))) * coeff * newcf !Contribution to value of matrix element <j,c| H | loc, c>  
                            
        !                     if(colIndx == rowIndx) then !Matrix element contributes to diagonal Hamiltonian
        !                         if(cntr_di(rowIndx, 1) == 1) then !Not the first diagonal element of state j. Add to previous element inDicated by cntr_di(rowIndx, 2). 
        !                             hamDi_temp(cntr_di(rowIndx, 2)) = hamDi_temp(cntr_di(rowIndx, 2)) + h_add
        !                         else !First diagonal element of state j. Create new element and save position in cntr_di(j, 2).
        !                             nDi                 = nDi + 1 !Number of diagonal matrix elements
        !                             hamDi_temp(nDi)     = hamDi_temp(nDi) + h_add !Diagonal matrix element
        !                             rcdi_temp(nDi)      = rowIndx !Row and column index of matrix element 
        !                             cntr_di(rowIndx, 1) = 1 !Flag for existing diagonal matrix element in row 'rowIndx'
        !                             cntr_di(rowIndx, 2) = nDi !Position of existing element of 'rowIndx' in array 'hamDi_temp'
        !                         end if 
        !                     else !Off-diagonal elements 
        !                         dbl = 0 !Flag for whether a new matrix element should be created (0) or an old matrix element was updated (1)
        !                         if(cntrjc > 0) then !Not the first element in row j. Check elements for double entries. 
        !                             do i = rowst + 1, rowst + cntrjc + 1 !Go from beginning of row j to current element. 
        !                                 if(rc_temp(i,2) == colIndx .and. rc_temp(i,1) == rowIndx) then !If yes, update existing element.
        !                                     if(i > arrsize) stop 'Subroutine hopping_irrep2D: i > arrsize'
        !                                     ham_temp(i)    = ham_temp(i) + h_add 
        !                                     parities_temp(i)   = (-1)**parity
        !                                     dplcts_temp(i) = dplcts_temp(i) + 1 
        !                                     dbl            = 1 !Found and updated existing element.       
        !                                 end if
        !                             end do
        !                         end if
        !                         if(dbl == 0) then !If no existing element was found, create new matrix element.
        !                             ! n_temp(id) = n_temp(id) + 1 !Counter for non-zero matrix elements.       
        !                             cntrjc           = cntrjc + 1 !One more element found for row j, basis state c. 
        !                             n_temp           = n_temp + 1 !Counter for non-zero matrix elements. Needed later.
        !                             pos              = rowst + cntrjc !Position corresponds to position in row j. 
        !                             if(pos > arrsize) stop 'Subroutine hopping_irrep2D: pos > arrsize'
        !                             ham_temp(pos)    = ham_temp(pos) + h_add !Value of matrix element. 
        !                             parities_temp(pos)   = (-1)**parity !Save parities for later. 
        !                             dplcts_temp(pos) = dplcts_temp(pos) + 1 !Save number of dplcts for later. 
        !                             rc_temp(pos,1)   = rowIndx !Row of matrix element. 
        !                             rc_temp(pos,2)   = colIndx !Column of matrix element.                       
        !                         end if
        !                     end if 
        !                 end do !Loop over scattered irrep basis states (d) 
        !             end do !Loop over symmetry orbit (s)
        !         end do !Loop over irrep basis states (c)
        !     end do !Loop over representatives (j)

        !     nnz = n_temp
            

        !     if(allocated(rc))     deallocate(rc)
        !     if(allocated(ham))    deallocate(ham)
        !     if(allocated(rcdi))   deallocate(rcdi)
        !     if(allocated(hamDi))  deallocate(hamDi)
        !     if(allocated(parities))   deallocate(parities)
        !     if(allocated(dplcts)) deallocate(dplcts)

        !     allocate(ham(nnz))
        !     allocate(rc(nnz,2))
        !     allocate(rcdi(nDi))
        !     allocate(hamDi(nDi))
        !     allocate(parities(nnz))
        !     allocate(dplcts(nnz))

        !     rc    = 0
        !     cntr  = 0
        !     rcdi  = 0
        !     ham   = 0.d0  
        !     hamDi = 0.d0  

        !     do i = 1, arrsize!nnz
        !         if(rc_temp(i,1) > 0) then 
        !             cntr         = cntr + 1
        !             ham(cntr)    = ham_temp(i)
        !             rc(cntr,1)   = rc_temp(i,1)
        !             rc(cntr,2)   = rc_temp(i,2)
        !             parities(cntr)   = parities_temp(i)
        !             dplcts(cntr) = dplcts_temp(i)        
        !         end if 
        !     end do
        !     ham = ham / dble(order)

        !     do i = 1, nDi !Fill up diagonal Hamiltonian
        !         hamDi(i) = hamDi_temp(i)
        !         rcdi(i)  = rcdi_temp(i)
        !     end do
        !     hamDi = hamDi / order 

        !     print*, 'Generated hopping Hamiltonian'
        !     print*, 'Number of non-zero matrix elements: ', nnz
        !     print*, 'Number of diagonal matrix elements: ', nDi
        !     print*, 'Number of representatives: ', dim
        !     print*, ''

    ! end subroutine hopping_irrep2D
    
    !With variable arguments
    subroutine hopping_irrep2D(threads, tilted, nHel, tilt, lx, ly, sites, nbonds, dim, basis, orbsize, orbits, phases, norm, bsites, xtransl, ytransl, symm, id, par, rot, refl, c6, t, rc, rcdi, parities, dplcts, ham, hamDi, nnz, nDi)
    
        !------------------------------------------------------------------------!
        !          Hopping procedure for 2D irreducible representations          !
        !------------------------------------------------------------------------!
        ! This subroutine generates the hopping Hamiltonian for irreducible representations in 2D.
        ! Work in progress.

        use parameters
        implicit none
        integer, intent(in):: orbsize, threads, tilted, nHel, tilt, lx, ly, sites, nbonds, symm, bsites(2, nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
        integer(kind=8), intent(in) :: dim, basis(dim)
        integer(kind=8), allocatable, intent(in) :: orbits(:,:,:)
        double precision, intent(in) :: t, id, par(6), rot(5), norm(dim,2)
        double complex, allocatable, intent(in) :: phases(:,:,:)
        
        integer,         allocatable, intent(out) :: parities(:), dplcts(:)
        integer(kind=8),              intent(out) :: nnz, nDi 
        integer(kind=8), allocatable, intent(out) :: rc(:,:), rcdi(:)
        double complex,  allocatable, intent(out) :: ham(:), hamDi(:)
        
        !Local variables
        integer :: s, c, d, o, l1, l2, parity, dbl, order, site1, site2
        integer(kind=8) :: i, j, loc, rep, pos, rowst, rowIndx, colIndx, newst, state, cntr, cntrjc, n_temp, arrsize 
        integer(kind=8), allocatable :: rc_temp(:,:), rcdi_temp(:), cntr_di(:,:)
        integer        , allocatable :: parities_temp(:), dplcts_temp(:)
        double precision, parameter :: tol = 1.0e-14
        double precision :: sign = 1, signrep = 1 
        double complex :: h_add, coeff, newcf 
        double complex, allocatable :: ham_temp(:), hamDi_temp(:)

        
        if(symm == 1) then 
            order = 12
        else 
            order = 1
        end if 

        arrsize = dim * orbsize * 2 !Each representative (dim) has two basis states (2), each of which is a linear combination of at most 'orbsize' states. 
        if(allocated(ham_temp))    deallocate(ham_temp)
        if(allocated(rc_temp))     deallocate(rc_temp)
        if(allocated(hamDi_temp))  deallocate(hamDi_temp)
        if(allocated(rcdi_temp))   deallocate(rcdi_temp)
        if(allocated(parities_temp))   deallocate(parities_temp)
        if(allocated(dplcts_temp)) deallocate(dplcts_temp)
        if(allocated(cntr_di))     deallocate(cntr_di)
        
        allocate(dplcts_temp(arrsize))
        allocate(hamDi_temp(arrsize))
        allocate(rc_temp(arrsize, 2))
        allocate(rcdi_temp(arrsize))
        allocate(parities_temp(arrsize))
        allocate(ham_temp(arrsize)) 
        allocate(cntr_di(dim, 2))

        h_add       = 0.d0
        ham_temp    = 0.d0  
        hamDi_temp  = 0.d0  
        rc_temp     = 0 
        rcdi_temp   = 0 
        parities_temp   = 0 
        dplcts_temp = 0 
        n_temp      = 0 
        nDi         = 0 
        nnz         = 0 
        cntr_di     = 0 
        site1       = bsites(1,1)
        site2       = bsites(2,1)
        
        !$omp parallel do default(firstprivate) shared(ham_temp, bsites, basis) num_threads(threads)
        do j = 1, dim !All representatives 
            !$ id = omp_get_thread_num()
            
            do c = 1, 2 !First and second basis state of 2D irrep 
            
                h_add = 0.d0 !Scattering contribution to matrix element
                cntrjc = 0 !Counts the number of already calculated matrix elements of each combination row 'j' and basist state 'c', i.e. tuple (j, c) 
                rowIndx = 2*(j-1) + c
                rowst = (rowIndx-1) * orbsize !Position-1 of beginning of tuple (j,c)   
                do s = 1, orbsize !Loop through orbit of state basis(j)
                    state = orbits(j, s, c)
                    coeff = phases(j, s, c)
                    if(state == 0) exit 
                    if(abs(dble(coeff)) <= tol .and. abs(aimag(coeff)) <= tol) cycle 
                    if(btest(state, site1-1) .and. .not.(btest(state, site2-1))) then !Hopping from 1->2
                        newst  = ibclr(ibset(state, site2-1), site1-1)            
                        parity = popcnt(ibits(state, site1, sites)) + popcnt(ibits(ibclr(state, site1 - 1), site2, sites))      
                    else if(btest(state, site2-1) .and. .not.(btest(state, site1-1))) then !Hopping from 2->1
                        newst  = ibclr(ibset(state, site1-1), site2-1) 
                        parity = popcnt(ibits(state, site2, sites)) + popcnt(ibits(ibclr(state, site2-1), site1, sites))
                    else 
                        cycle 
                    end if 
                    if(tilted == 0) then 
                        call representative_irrep_rect(newst, sites, Lx, Ly, symm, id, par, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                    else if(tilted == 1) then 
                        call representative_irrep(newst, sites, nHel, tilt, Lx, Ly, symm, id, par, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                    end if 
                    call findstate(dim, rep, basis, loc) !Finds the location of representative in basis
                    if(loc <= 0) cycle !New state is compatible with momentum and symmetry 
                    do d = 1, 2 !Irrep basis states of scattered representative
                        colIndx = 2*(loc-1) + d
                        newcf   = 0.d0 
                        do o = 1, orbsize !Search for position of scattered state in orbit of Irrep's basis state and assign corresponDing coefficient 'newcf'   
                            if(orbits(loc, o, d) == newst) then 
                                newcf = phases(loc, o, d)
                                exit 
                            end if 
                        end do 
                        if(abs(dble(newcf)) <= tol .and. abs(aimag(newcf)) <= tol) cycle 

                        sign  = (-1) * (-1)**parity     
                        h_add = sign * t * sqrt(dble(norm(loc, d))/dble(norm(j, c))) * coeff * newcf !Contribution to value of matrix element <j,c| H | loc, c>  
                        
                        if(colIndx == rowIndx) then !Matrix element contributes to diagonal Hamiltonian
                            if(cntr_di(rowIndx, 1) == 1) then !Not the first diagonal element of state j. Add to previous element inDicated by cntr_di(rowIndx, 2). 
                                hamDi_temp(cntr_di(rowIndx, 2)) = hamDi_temp(cntr_di(rowIndx, 2)) + h_add
                            else !First diagonal element of state j. Create new element and save position in cntr_di(j, 2).
                                nDi                 = nDi + 1 !Number of diagonal matrix elements
                                hamDi_temp(nDi)     = hamDi_temp(nDi) + h_add !Diagonal matrix element
                                rcdi_temp(nDi)      = rowIndx !Row and column index of matrix element 
                                cntr_di(rowIndx, 1) = 1 !Flag for existing diagonal matrix element in row 'rowIndx'
                                cntr_di(rowIndx, 2) = nDi !Position of existing element of 'rowIndx' in array 'hamDi_temp'
                            end if 
                        else !Off-diagonal elements 
                            dbl = 0 !Flag for whether a new matrix element should be created (0) or an old matrix element was updated (1)
                            if(cntrjc > 0) then !Not the first element in row j. Check elements for double entries. 
                                do i = rowst + 1, rowst + cntrjc + 1 !Go from beginning of row j to current element. 
                                    if(rc_temp(i,2) == colIndx .and. rc_temp(i,1) == rowIndx) then !If yes, update existing element.
                                        if(i > arrsize) stop 'Subroutine hopping_irrep2D: i > arrsize'
                                        ham_temp(i)    = ham_temp(i) + h_add 
                                        parities_temp(i)   = (-1)**parity
                                        dplcts_temp(i) = dplcts_temp(i) + 1 
                                        dbl            = 1 !Found and updated existing element.       
                                    end if
                                end do
                            end if
                            if(dbl == 0) then !If no existing element was found, create new matrix element.
                                ! n_temp(id) = n_temp(id) + 1 !Counter for non-zero matrix elements.       
                                cntrjc           = cntrjc + 1 !One more element found for row j, basis state c. 
                                n_temp           = n_temp + 1 !Counter for non-zero matrix elements. Needed later.
                                pos              = rowst + cntrjc !Position corresponds to position in row j. 
                                if(pos > arrsize) stop 'Subroutine hopping_irrep2D: pos > arrsize'
                                ham_temp(pos)    = ham_temp(pos) + h_add !Value of matrix element. 
                                parities_temp(pos)   = (-1)**parity !Save parities for later. 
                                dplcts_temp(pos) = dplcts_temp(pos) + 1 !Save number of dplcts for later. 
                                rc_temp(pos,1)   = rowIndx !Row of matrix element. 
                                rc_temp(pos,2)   = colIndx !Column of matrix element.                       
                            end if
                        end if 
                    end do !Loop over scattered irrep basis states (d) 
                end do !Loop over symmetry orbit (s)
            end do !Loop over irrep basis states (c)
        end do !Loop over representatives (j)

        nnz = n_temp
        

        if(allocated(rc))     deallocate(rc)
        if(allocated(ham))    deallocate(ham)
        if(allocated(rcdi))   deallocate(rcdi)
        if(allocated(hamDi))  deallocate(hamDi)
        if(allocated(parities))   deallocate(parities)
        if(allocated(dplcts)) deallocate(dplcts)

        allocate(ham(nnz))
        allocate(rc(nnz,2))
        allocate(rcdi(nDi))
        allocate(hamDi(nDi))
        allocate(parities(nnz))
        allocate(dplcts(nnz))

        rc    = 0
        cntr  = 0
        rcdi  = 0
        ham   = 0.d0  
        hamDi = 0.d0  

        do i = 1, arrsize!nnz
            if(rc_temp(i,1) > 0) then 
                cntr         = cntr + 1
                ham(cntr)    = ham_temp(i)
                rc(cntr,1)   = rc_temp(i,1)
                rc(cntr,2)   = rc_temp(i,2)
                parities(cntr)   = parities_temp(i)
                dplcts(cntr) = dplcts_temp(i)        
            end if 
        end do
        ham = ham / dble(order)

        do i = 1, nDi !Fill up diagonal Hamiltonian
            hamDi(i) = hamDi_temp(i)
            rcdi(i)  = rcdi_temp(i)
        end do
        hamDi = hamDi / order 

        print*, 'Generated hopping Hamiltonian'
        print*, 'Number of non-zero matrix elements: ', nnz
        print*, 'Number of diagonal matrix elements: ', nDi
        print*, 'Number of representatives: ', dim
        print*, ''

    end subroutine hopping_irrep2D

    ! subroutine diagonal(flag)

        !     implicit none
        !     integer, intent(in) :: flag ! Integer flag to indicate this routine
        !     ! integer, intent(in)          :: unit, sites, nbonds, nbonds2, bsites(2, nbonds), bsites2(2, nbonds2)
        !     ! integer(kind=8), intent(in)  :: dim, basis(dim)
        !     ! character(len=*), intent(in) :: parameters

        !     ! integer(kind=8), intent(out)              :: nnz
        !     ! integer(kind=8), allocatable, intent(out) :: ham(:,:), occ(:,:)

        !     !Local variables
        !     integer(kind=8) :: j = 0, m = 0, s = 0, counter_v = 0, counter_v2 = 0

        !     if(allocated(occ))   deallocate(occ)
        !     if(allocated(ham))   deallocate(ham)
        !     allocate(occ(sites, dim))
        !     allocate(ham(dim, 2))
        !     occ = 0 
        !     ham = 0 
        !     do j = 1, dim
        !         counter_v  = 0
        !         counter_v2 = 0
        !         do m = 0, sites - 1 !m goes through all digits of each configuration
        !             if(btest(basis(j), m)) occ(m + 1, j) = occ(m + 1, j) + 1
        !         end do
        !         do s = 1, nbonds !Loop runs over all nn-bonds
        !             if(btest(basis(j), bsites(1, s) - 1) .and. &
        !                 btest(basis(j), bsites(2, s) - 1)) counter_v = counter_v + 1
        !             if(.not. (btest(basis(j), bsites(1, s) - 1)) .and. &
        !                 .not. (btest(basis(j), bsites(2, s) - 1))) counter_v = counter_v + 1
        !             if(btest(basis(j), bsites(1, s) - 1) .and. &
        !                 .not. (btest(basis(j), bsites(2, s) - 1))) counter_v = counter_v - 1
        !             if(.not. (btest(basis(j), bsites(1, s) - 1)) .and. &
        !                 btest(basis(j), bsites(2, s) - 1)) counter_v = counter_v - 1
        !         end do

        !         do s = 1, nbonds2 !Loop runs over all nnn-bonds
        !             if(btest(basis(j), bsites2(1, s) - 1) .and. &
        !                 btest(basis(j), bsites2(2, s) - 1)) counter_v2 = counter_v2 + 1
        !             if(.not. (btest(basis(j), bsites2(1, s) - 1)) .and. &
        !                 .not. (btest(basis(j), bsites2(2, s) - 1))) counter_v2 = counter_v2 + 1
        !             if(btest(basis(j), bsites2(1, s) - 1) .and. &
        !                 .not. (btest(basis(j), bsites2(2, s) - 1))) counter_v2 = counter_v2 - 1
        !             if(.not. (btest(basis(j), bsites2(1, s) - 1)) .and. &
        !                 btest(basis(j), bsites2(2, s) - 1)) counter_v2 = counter_v2 - 1
        !         end do
    
        !         ham(j, 1) = ham(j, 1) + counter_v
        !         ham(j, 2) = ham(j, 2) + counter_v2
        !     end do

        !     nnz = dim

        !     print*, 'Generated diagonal Hamiltonian.'
        !     print*, 'Number of non-zero matrix elements: ', nnz
        !     print*, ''

    ! end subroutine diagonal

    ! Version with variable arguments
    subroutine diagonal(unit, sites, nbonds, nbonds2, dim, bsites, bsites2, basis, ham, occ, nnz)

        implicit none

        integer, intent(in)          :: unit, sites, nbonds, nbonds2, bsites(2, nbonds), bsites2(2, nbonds2)
        integer(kind=8), intent(in)  :: dim, basis(dim)

        integer(kind=8), intent(out)              :: nnz
        integer(kind=8), allocatable, intent(out) :: ham(:,:), occ(:,:)

        !Local variables
        integer(kind=8) :: j = 0, m = 0, s = 0, counter_v = 0, counter_v2 = 0

        if(allocated(occ))   deallocate(occ)
        if(allocated(ham))   deallocate(ham)
        allocate(occ(sites, dim))
        allocate(ham(dim, 2))
        occ = 0 
        ham = 0 
        do j = 1, dim
            counter_v  = 0
            counter_v2 = 0
            do m = 0, sites - 1 !m goes through all digits of each configuration
                if(btest(basis(j), m)) occ(m + 1, j) = occ(m + 1, j) + 1
            end do
            do s = 1, nbonds !Loop runs over all nn-bonds
                if(btest(basis(j), bsites(1, s) - 1) .and. &
                    btest(basis(j), bsites(2, s) - 1)) counter_v = counter_v + 1
                if(.not. (btest(basis(j), bsites(1, s) - 1)) .and. &
                    .not. (btest(basis(j), bsites(2, s) - 1))) counter_v = counter_v + 1
                if(btest(basis(j), bsites(1, s) - 1) .and. &
                    .not. (btest(basis(j), bsites(2, s) - 1))) counter_v = counter_v - 1
                if(.not. (btest(basis(j), bsites(1, s) - 1)) .and. &
                    btest(basis(j), bsites(2, s) - 1)) counter_v = counter_v - 1
            end do

            do s = 1, nbonds2 !Loop runs over all nnn-bonds
                if(btest(basis(j), bsites2(1, s) - 1) .and. &
                    btest(basis(j), bsites2(2, s) - 1)) counter_v2 = counter_v2 + 1
                if(.not. (btest(basis(j), bsites2(1, s) - 1)) .and. &
                    .not. (btest(basis(j), bsites2(2, s) - 1))) counter_v2 = counter_v2 + 1
                if(btest(basis(j), bsites2(1, s) - 1) .and. &
                    .not. (btest(basis(j), bsites2(2, s) - 1))) counter_v2 = counter_v2 - 1
                if(.not. (btest(basis(j), bsites2(1, s) - 1)) .and. &
                    btest(basis(j), bsites2(2, s) - 1)) counter_v2 = counter_v2 - 1
            end do
 
            ham(j, 1) = ham(j, 1) + counter_v
            ham(j, 2) = ham(j, 2) + counter_v2
        end do

        nnz = dim

        print*, 'Generated diagonal Hamiltonian.'
        print*, 'Number of non-zero matrix elements: ', nnz
        print*, ''

    end subroutine diagonal

    subroutine diagonal_dp(unit, sites, nbonds, nbonds2, dim, bsites, bsites2, basis, ham, occ, nnz)

        implicit none

        integer, intent(in)          :: unit, sites, nbonds, nbonds2, bsites(2, nbonds), bsites2(2, nbonds2)
        integer(kind=8), intent(in)  :: dim, basis(dim)
        
        integer(kind=8), intent(out)               :: nnz
        integer(kind=8), allocatable, intent(out)  :: occ(:,:)
        double precision, allocatable, intent(out) :: ham(:,:)

        integer(kind=8)  :: j = 0, m = 0, s = 0
        double precision :: counter_v = 0.d0, counter_v2 = 0.d0
        

        if(allocated(occ)) deallocate(occ)
        if(allocated(ham)) deallocate(ham)
        allocate(occ(sites, dim))
        allocate(ham(dim, 2))
        occ = 0 
        ham = 0.d0  
        do j = 1, dim
            counter_v  = 0
            counter_v2 = 0
            do m = 0, sites - 1 !m goes through all digits of each configuration
                if(btest(basis(j), m)) occ(m + 1, j) = occ(m + 1, j) + 1
            end do
            do s = 1, nbonds !Loop runs over all nn-bonds
                if(btest(basis(j), bsites(1, s) - 1) .and. &
                    btest(basis(j), bsites(2, s) - 1)) counter_v = counter_v + 1
                if(.not. (btest(basis(j), bsites(1, s) - 1)) .and. &
                    .not. (btest(basis(j), bsites(2, s) - 1))) counter_v = counter_v + 1
                if(btest(basis(j), bsites(1, s) - 1) .and. &
                    .not. (btest(basis(j), bsites(2, s) - 1))) counter_v = counter_v - 1
                if(.not. (btest(basis(j), bsites(1, s) - 1)) .and. &
                    btest(basis(j), bsites(2, s) - 1)) counter_v = counter_v - 1
            end do
            do s = 1, nbonds2 !Loop runs over all nnn-bonds
                if(btest(basis(j), bsites2(1, s) - 1) .and. &
                    btest(basis(j), bsites2(2, s) - 1)) counter_v2 = counter_v2 + 1
                if(.not. (btest(basis(j), bsites2(1, s) - 1)) .and. &
                    .not. (btest(basis(j), bsites2(2, s) - 1))) counter_v2 = counter_v2 + 1
                if(btest(basis(j), bsites2(1, s) - 1) .and. &
                    .not. (btest(basis(j), bsites2(2, s) - 1))) counter_v2 = counter_v2 - 1
                if(.not. (btest(basis(j), bsites2(1, s) - 1)) .and. &
                    btest(basis(j), bsites2(2, s) - 1)) counter_v2 = counter_v2 - 1
            end do

            ham(j, 1) = ham(j, 1) + counter_v
            ham(j, 2) = ham(j, 2) + counter_v2
        end do

        nnz = dim

        print*,'Generated diagonal Hamiltonian.'
        print*,'Number of non-zero matrix elements: ', nnz
        print*, ''

    end subroutine diagonal_dp

    ! subroutine diagonal_irrep2D(flag1, flag2)

        !     implicit none
        !     integer, intent(in) :: flag2 ! Integer flag for 2D irreps
        !     double precision, intent(in) :: flag1 ! Double precision flag 
        !     ! integer, intent(in)                       :: sites, orbsize
        !     ! integer(kind=8), intent(in)               :: dim
        !     ! integer, allocatable, intent(in)          :: trisites(:,:), hexsites(:,:)
        !     ! integer(kind=8), allocatable, intent(in)  :: orbits(:,:,:)
        !     ! double precision, allocatable, intent(in) :: norm(:,:)
        !     ! double complex, allocatable, intent(in)   :: phases(:,:,:)
            
        !     ! integer(kind=8), allocatable, intent(out) :: occ(:,:)
        !     ! double complex, allocatable, intent(out)  :: ham(:,:)

        !     integer          :: m, s, c, o, site1, site2 
        !     integer(kind=8)  :: j = 0, arrsize = 0, state, rowIndx
        !     double precision :: cntr_v = 0.d0, cntr_v2 = 0.d0
            
        !     arrsize = 2 * dim !Each representative has two basis states in 2D irrep.
        !     if(allocated(occ)) deallocate(occ)
        !     if(allocated(ham)) deallocate(ham)
        !     allocate(occ(sites,arrsize))
        !     allocate(ham(arrsize,2))
        !     occ = 0 
        !     ham = 0.d0  
        !     do j = 1, dim !Representatives
        !         do c = 1, 2 !First and second basis state of 2D irrep 
        !             cntr_v  = 0
        !             cntr_v2 = 0
        !             rowIndx = 2*(j-1) + c
        !             do o = 1, orbsize
        !                 state = orbits(j, o, c)
        !                 do m = 0, sites - 1 !m goes through all digits of each configuration
        !                     if(btest(state, m)) occ(m + 1, rowIndx) = occ(m + 1, rowIndx) + 1
        !                 end do

        !                 do s = 1, 3 !Loop runs over the three different NN-bonds
        !                     site1 = trisites(1,s)
        !                     site2 = trisites(2,s)
                        
        !                     if(btest(state, site1-1) .eqv. btest(state, site2-1)) cntr_v = cntr_v + 1
        !                     if(btest(state, site1-1) .neqv. btest(state, site2-1)) cntr_v = cntr_v - 1
        !                 end do

        !                 site1 = hexsites(1, s)
        !                 site2 = hexsites(2, s)
                        
        !                 if(btest(state, site1-1) .eqv. btest(state, site2-1))  cntr_v2 = cntr_v2 + 1
        !                 if(btest(state, site1-1) .neqv. btest(state, site2-1)) cntr_v2 = cntr_v2 - 1
                        
        !                 ham(rowIndx, 1) = ham(rowIndx, 1) + phases(j, o, c) * cntr_v / sqrt(norm(j, c))
        !                 ham(rowIndx, 2) = ham(rowIndx, 2) + phases(j, o, c) * cntr_v2 / sqrt(norm(j, c))
        !             end do !Loop over orbit of irrep basis states (o)
        !         end do !Loop over irrep basis states (c)
        !     end do !Loop over representatives (j)

        !     print*,'Generated diagonal Hamiltonian.'
        !     print*, ''

    ! end subroutine diagonal_irrep2D
    
    ! Version with variable arguments
    subroutine diagonal_irrep2D(sites, dim, trisites, hexsites, orbsize, orbits, phases, norm, ham, occ)

        implicit none

        integer, intent(in)                       :: sites, orbsize
        integer(kind=8), intent(in)               :: dim
        integer, allocatable, intent(in)          :: trisites(:,:), hexsites(:,:)
        integer(kind=8), allocatable, intent(in)  :: orbits(:,:,:)
        double precision, allocatable, intent(in) :: norm(:,:)
        double complex, allocatable, intent(in)   :: phases(:,:,:)
        
        integer(kind=8), allocatable, intent(out) :: occ(:,:)
        double complex, allocatable, intent(out)  :: ham(:,:)

        integer          :: m, s, c, o, site1, site2 
        integer(kind=8)  :: j = 0, arrsize = 0, state, rowIndx
        double precision :: cntr_v = 0.d0, cntr_v2 = 0.d0
        
        arrsize = 2 * dim !Each representative has two basis states in 2D irrep.
        if(allocated(occ)) deallocate(occ)
        if(allocated(ham)) deallocate(ham)
        allocate(occ(sites,arrsize))
        allocate(ham(arrsize,2))
        occ = 0 
        ham = 0.d0  
        do j = 1, dim !Representatives
            do c = 1, 2 !First and second basis state of 2D irrep 
                cntr_v  = 0
                cntr_v2 = 0
                rowIndx = 2*(j-1) + c
                do o = 1, orbsize
                    state = orbits(j, o, c)
                    do m = 0, sites - 1 !m goes through all digits of each configuration
                        if(btest(state, m)) occ(m + 1, rowIndx) = occ(m + 1, rowIndx) + 1
                    end do

                    do s = 1, 3 !Loop runs over the three different NN-bonds
                        site1 = trisites(1,s)
                        site2 = trisites(2,s)
                    
                        if(btest(state, site1-1) .eqv. btest(state, site2-1)) cntr_v = cntr_v + 1
                        if(btest(state, site1-1) .neqv. btest(state, site2-1)) cntr_v = cntr_v - 1
                    end do

                    site1 = hexsites(1, s)
                    site2 = hexsites(2, s)
                    
                    if(btest(state, site1-1) .eqv. btest(state, site2-1))  cntr_v2 = cntr_v2 + 1
                    if(btest(state, site1-1) .neqv. btest(state, site2-1)) cntr_v2 = cntr_v2 - 1
                    
                    ham(rowIndx, 1) = ham(rowIndx, 1) + phases(j, o, c) * cntr_v / sqrt(norm(j, c))
                    ham(rowIndx, 2) = ham(rowIndx, 2) + phases(j, o, c) * cntr_v2 / sqrt(norm(j, c))
                end do !Loop over orbit of irrep basis states (o)
            end do !Loop over irrep basis states (c)
        end do !Loop over representatives (j)

        print*,'Generated diagonal Hamiltonian.'
        print*, ''

    end subroutine diagonal_irrep2D

    subroutine diagonal_p(parity, p1, p2, p3, refl1, refl2, refl3, sites, nbonds, nbonds2, dim, bsites, bsites2, basis, ham, occ, nnz)

        implicit none

        integer(kind=8), intent(in) :: dim, basis(dim)
        integer, intent(in) :: parity, p1, p2, p3, sites, nbonds, nbonds2, bsites(2, nbonds), bsites2(2, nbonds2), refl1(sites), refl2(sites), refl3(sites)
        
        integer(kind=8), intent(out) :: nnz
        integer(kind=8), allocatable, intent(out) :: occ(:,:), ham(:,:)

        !Local variables
        integer :: info = 0
        integer(kind=8) :: j = 0, m = 0, s = 0, state = 0, nrefl = 0, p = 0, counter_v = 0, counter_v2 = 0
        double precision :: sign = 1.d0

        if(allocated(occ)) deallocate(occ)
        if(allocated(ham)) deallocate(ham)
        allocate(occ(sites, dim))
        allocate(ham(dim, 2))
        occ = 0 
        ham = 0.d0  

        if(parity == 0) then 
            nrefl = 1
        else 
            nrefl = 1 + abs(p1) + abs(p2) + abs(p3) + abs(p1*p2) + abs(p1*p3)
        end if 

        do j = 1, dim
            counter_v  = 0
            counter_v2 = 0
            do p = 1, nrefl
                
                if(p == 2) call reflect(int(1,8), basis(j), sites, refl1, sign, info, state)
                if(p == 3) call reflect(int(1,8), basis(j), sites, refl2, sign, info, state)
                if(p == 4) call reflect(int(1,8), basis(j), sites, refl3, sign, info, state)
                if(p == 5) then 
                    call reflect(int(1,8), basis(j), sites, refl1, sign, info, state)
                    call reflect(int(1,8), state, sites, refl2, sign, info, state)
                else if(p == 6) then 
                    call reflect(int(1,8), basis(j), sites, refl1, sign, info, state)
                    call reflect(int(1,8), state, sites, refl3, sign, info, state)
                end if 
                if(state == basis(j)) cycle 
                if(p == 1) state = basis(j)

                do m = 0, sites - 1 !m goes through all digits of each configuration
                    if(btest(state, m)) occ(m + 1, j) = occ(m + 1, j) + 1
                end do
                do s = 1, nbonds !Loop runs over all nn-bonds
                    if(btest(state, bsites(1, s) - 1) .and. &
                        btest(state, bsites(2, s) - 1)) counter_v = counter_v + 1
                    if(.not.(btest(state, bsites(1, s) - 1)) .and. &
                        .not.(btest(state, bsites(2, s) - 1))) counter_v = counter_v + 1
                    if(btest(state, bsites(1, s) - 1) .and. &
                        .not.(btest(state, bsites(2, s) - 1))) counter_v = counter_v - 1
                    if(.not.(btest(state, bsites(1, s) - 1)) .and. &
                        btest(state, bsites(2, s) - 1)) counter_v = counter_v - 1
                end do
                do s = 1, nbonds2 !Loop runs over all nnn-bonds
                    if(btest(state, bsites2(1, s) - 1) .and. &
                        btest(state, bsites2(2, s) - 1)) counter_v2 = counter_v2 + 1
                    if(.not.(btest(state, bsites2(1, s) - 1)) .and. &
                        .not.(btest(state, bsites2(2, s) - 1))) counter_v2 = counter_v2 + 1
                    if(btest(state, bsites2(1, s) - 1) .and. &
                        .not.(btest(state, bsites2(2, s) - 1))) counter_v2 = counter_v2 - 1
                    if(.not.(btest(state, bsites2(1, s) - 1)) .and. &
                        btest(state, bsites2(2, s) - 1)) counter_v2 = counter_v2 - 1
                end do

                ham(j, 1) = ham(j, 1) + counter_v
                ham(j, 2) = ham(j, 2) + counter_v2

            end do 
        end do

        nnz = dim

        print*,'Generated diagonal Hamiltonian.'
        print*, ''

    end subroutine diagonal_p

    !---------------------------------------------!
    !            Unify sparse Hamiltonian         !
    !---------------------------------------------!

    subroutine unify(t, v1, v2, w, mass, pattern, sites, occ, nOff, nDi, hamOff, hamDi, ham, rc, nnz)

        implicit none
        ! save 
        integer, intent(in):: sites, nOff, nDi, occ(sites,nDi)
        integer, intent(in):: hamOff(nOff,3)
        integer, intent(in):: hamDi(nDi,2)
        double precision, intent(in) :: t, v1, v2, w, mass 
        character*2, intent(in) :: pattern

        integer, intent(out)::  nnz
        integer, allocatable, intent(out):: rc(:,:)
        double precision, allocatable, intent(out):: ham(:)

        integer :: j = 0, ab = 0 
        integer :: slstaggering(sites)
        double precision :: rand(sites)

        if(allocated(ham)) deallocate(ham)
        if(allocated(rc))  deallocate(rc)
        nnz = nOff + nDi !Number of non-zero elements
        allocate(ham(nnz))
        allocate(rc(nnz,2))
        
        ham = 0
        rc  = 0
        if(pattern == "AB") then 
            ab = 1 
        else if(pattern == "BA") then 
            ab = 0
        end if 

        !Staggered sublattice potential: +1 on A lattice, -1 on B lattice
        do j = 1, sites 
            slstaggering(j) = (-1)**(j-ab)
        end do 
        
        call random_number(rand)

        rand = 2 * (rand - 0.5)

        ham(1:nOff)  = - t * hamOff(1:nOff,1)
        rc(1:nOff,1) = hamOff(1:nOff,2)
        rc(1:nOff,2) = hamOff(1:nOff,3)
        do j = 1, nDi !Fill diagonal Hamiltonian
            ham(nOff+j)  = v1 * 0.25 * hamDi(j,1) + v2 * 0.25 * hamDi(j,2) + mass * dot_product(slstaggering, occ(1:sites,j)) + w * sum(rand * occ(1:sites,j))
            ! if(j == 1) ham(nOff+j)  = ham(nOff+j) + 0.01
            rc(nOff+j,1) = j
            rc(nOff+j,2) = j
        end do
        
        return
        
    end subroutine unify

    subroutine unify_dp(v1, v2, w, mass, pattern, sites, occ, nOff, nDi, nDi2, hamOff, rcoff, rcdi, hamDi2, hamDi, ham, rc, nnz)

        implicit none

        integer, intent(in):: sites, nOff, nDi, nDi2, occ(sites,nDi)
        double precision, intent(in):: hamOff(nOff), hamDi2(nDi2)
        integer, intent(in):: rcoff(nOff,2), rcdi(nDi2)
        integer, intent(in):: hamDi(nDi,2)
        double precision, intent(in) :: v1, v2, w, mass 
        character*2, intent(in) :: pattern

        integer, intent(out)::  nnz
        integer, allocatable, intent(out):: rc(:,:)
        double precision, allocatable, intent(out):: ham(:)

        integer :: j = 0, ab = 0 
        integer :: slstaggering(sites)
        double precision :: rand(sites)

        if(allocated(ham)) deallocate(ham)
        if(allocated(rc))  deallocate(rc)
        if(v1 == 0.d0 .and. v2 == 0.d0 .and. mass == 0.d0) then 
            nnz = nOff + nDi2
        else 
            nnz = nOff + max(nDi , nDi2) !Number of non-zero elements
        end if 
        allocate(ham(nnz))
        allocate(rc(nnz,2))   

        ham = 0.d0
        rc  = 0
        if(pattern == "AB") then 
            ab = 1 
        else if(pattern == "BA") then 
            ab = 0
        end if 

        !Staggered sublattice potential: +1 on A lattice, -1 on B lattice
        do j = 1, sites 
            slstaggering(j) = (-1)**(j-ab)
        end do   
        
        call random_number(rand)
        rand = 2 * (rand - 0.5)
        ham(1:nOff)  = hamOff(1:nOff)
        rc(1:nOff,1) = rcoff(1:nOff,1)
        rc(1:nOff,2) = rcoff(1:nOff,2)

        if(v1 == 0.d0 .and. v2 == 0.d0 .and. mass == 0.d0) go to 110 !Don't fill diagonal with zeroes
        do j = 1, nDi !Fill diagonal Hamiltonian        
            ham(nOff+j)  = v1 * 0.25 * hamDi(j,1) + v2 * 0.25 * hamDi(j,2) + mass * dot_product(slstaggering, occ(1:sites,j)) + w * sum(rand * occ(1:sites,j)) !- 3 * v1 * sum(occ(1:sites,j)) - 6 * v2 * sum(occ(1:sites,j))        
            rc(nOff+j,1) = j
            rc(nOff+j,2) = j
        end do
        
        
        if(nDi2 > 0) then 
            do j = 1, nDi2 
                ham(nOff+rcdi(j)) = ham(nOff+rcdi(j)) + hamDi2(j)
            end do 
        end if 
        110 continue 

        if(v1 == 0.d0 .and. v2 == 0.d0 .and. mass == 0.d0) then !Still need to fill diagonal entries from hopping 
            if(nDi2 > 0) then 
                do j = 1, nDi2 
                    ham(nOff+ j) = ham(nOff+ j) + hamDi2(j)
                    rc(nOff + j, 1) = rcdi(j) 
                    rc(nOff + j, 2) = rcdi(j)
                end do 
            end if 
        end if 


        return

    end subroutine unify_dp 

    subroutine unify_dc(v1, v2, w, mass, pattern, sites, occ, nOff, nDi, hamOff, rcoff, hamDi, ham, rc, nnz)

        implicit none

        integer, intent(in)          :: sites 
        integer(kind=8), intent(in)  :: nOff, nDi, occ(sites,nDi), rcoff(nOff,2), hamDi(nDi,2)
        double precision, intent(in) :: v1, v2, w, mass
        double complex, intent(in)   :: hamOff(nOff)
        character*2, intent(in)      :: pattern

        integer(kind=8), intent(out)              :: nnz
        integer(kind=8), allocatable, intent(out) :: rc(:,:)
        double complex, allocatable, intent(out)  :: ham(:)

        !Local variables
        integer          :: ab = 0, slstaggering(sites) 
        integer(kind=8)  :: j = 0
        double precision :: rand(sites)

        if(allocated(ham)) deallocate(ham)
        if(allocated(rc))  deallocate(rc)
        
        nnz = 0
        nnz = nOff + nDi !Number of non-zero elements
        allocate(ham(nnz))
        allocate(rc(nnz,2))    

        ham = 0
        rc  = 0
        if(pattern == "AB") then 
            ab = 1 
        else if(pattern == "BA") then 
            ab = 0
        end if 

        !Staggered sublattice potential: +1 on A lattice, -1 on B lattice
        do j = 1, sites 
            slstaggering(j) = int((-1)**(j-ab), kind=4)
        end do 
        
        call random_number(rand)
        rand         = 2 * (rand - 0.5)
        ham(1:nOff)  = hamOff(1:nOff)
        rc(1:nOff,1) = rcoff(1:nOff,1)
        rc(1:nOff,2) = rcoff(1:nOff,2)
        do j = 1, nDi !Fill diagonal Hamiltonian
            ham(nOff+j)  = v1 * 0.25 * hamDi(j,1) + v2 * 0.25 * hamDi(j,2) + mass * dot_product(slstaggering, occ(1:sites,j)) + w * sum(rand * occ(1:sites,j))
            rc(nOff+j,1) = j
            rc(nOff+j,2) = j
        end do

        return

    end subroutine unify_dc 

    subroutine unify_dc_2D(sites, nOff, nDi_off, dim, v1, v2, w, mass, pattern, occ, hamOff, rcoff, hamDi_off, rcdi, hamDi, ham, rc, nnz)

        implicit none

        integer, intent(in)                      :: sites
        integer(kind=8), intent(in)              :: dim, nOff, nDi_off
        integer(kind=8), allocatable, intent(in) :: rcoff(:,:), rcdi(:), occ(:,:)
        double precision, intent(in)             :: v1, v2, w, mass 
        double complex, allocatable, intent(in)  :: hamOff(:), hamDi_off(:), hamDi(:,:)
        character*2, intent(in)                  :: pattern

        integer(kind=8), intent(out)              :: nnz
        integer(kind=8), allocatable, intent(out) :: rc(:,:)
        double complex, allocatable, intent(out)  :: ham(:)

        integer          :: slstaggering(sites)
        integer(kind=8)  :: j, ab, nDi
        double precision :: rand(sites)

        if(allocated(ham)) deallocate(ham)
        if(allocated(rc))  deallocate(rc)
        nDi = 2*dim 

        if(v1 == 0.d0 .and. v2 == 0.d0 .and. mass == 0.d0) then 
            nnz = nOff + nDi_off
        else 
            nnz = nOff + max(nDi, nDi_off) !Number of non-zero elements
        end if 
        allocate(ham(nnz))
        allocate(rc(nnz,2))    

        ham = 0.d0
        rc  = 0
        if(pattern == "AB") then 
            ab = 1 
        else if(pattern == "BA") then 
            ab = 0
        end if 

        !Staggered sublattice potential: +1 on A lattice, -1 on B lattice
        do j = 1, sites 
            slstaggering(j) = int((-1)**(j-ab), kind=4)
        end do 
        
        call random_number(rand)
        rand         = 2 * (rand - 0.5)
        !Off diagonal matrix elements
        ham(1:nOff)  = hamOff(1:nOff)
        rc(1:nOff,1) = rcoff(1:nOff,1)
        rc(1:nOff,2) = rcoff(1:nOff,2)

        if(v1 == 0.d0 .and. v2 == 0.d0 .and. mass == 0.d0) go to 110 !Don't fill diagonal with zeroes
        do j = 1, nDi !Fill diagonal Hamiltonian
            ham(nOff+j)  = v1 * 0.25 * hamDi(j,1) + v2 * 0.25 * hamDi(j,2) + mass * dot_product(slstaggering, occ(1:sites,j)) + w * sum(rand * occ(1:sites,j))
            rc(nOff+j,1) = j
            rc(nOff+j,2) = j
        end do

        if(nDi_off > 0) then !Fill diagonal elements generated by hopping 
            do j = 1, nDi_off 
                ham(nOff+rcdi(j)) = ham(nOff+rcdi(j)) + hamDi_off(j)
            end do 
        end if 
        110 continue 

        if(v1 == 0.d0 .and. v2 == 0.d0 .and. mass == 0.d0) then !Still need to fill diagonal entries from hopping 
            if(nDi_off > 0) then 
                do j = 1, nDi_off 
                    ham(nOff+j) = ham(nOff+j) + hamDi_off(j)
                    rc(nOff+j, 1) = rcdi(j) 
                    rc(nOff+j, 2) = rcdi(j)
                end do 
            end if 
        end if 


        return

    end subroutine unify_dc_2D 

    subroutine unify_dense(dim, t, v1, v2, w, sites, occ, nOff, nDi, hamOff, hamDi, ham)

        implicit none

        integer, intent(in)          :: sites
        integer(kind=8), intent(in)  :: dim,  nOff, nDi, occ(sites,*), hamOff(nOff,3), hamDi(nDi,2)       
        double precision, intent(in) :: t, v1, v2, w

        double precision, allocatable, intent(out):: ham(:,:)

        integer(kind=8)  :: j = 0
        double precision :: rand(sites) 
        logical          :: symm 

        if(allocated(ham)) deallocate(ham)
        allocate(ham(dim, dim))
        ham = 0
        
        call random_number(rand)
        rand = 2 * (rand - 0.5)
        do j = 1, nOff
            ham(hamOff(j,2), hamOff(j,3))  = -t * hamOff(j,1)
        end do 
        do j = 1, nDi !Fill diagonal Hamiltonian
            ham(j,j) = v1 * hamDi(j,1) + v2 * hamDi(j,2) + w * sum(rand * occ(1:sites,j))
        end do
        call test_symm(symm, dim, ham)
        
        return

    end subroutine unify_dense

    subroutine unify_dense_dp(dim, v1, v2, w, sites, occ, nOff, nDi, nDi2, hamOff, rcoff, rcdi, hamDi2, hamDi, ham)
        
        implicit none

        integer, intent(in)          :: sites
        integer(kind=8), intent(in)  :: dim, nOff, nDi, nDi2, occ(sites,*), rcoff(nOff,2), rcdi(nDi2), hamDi(nDi,2)
        double precision, intent(in) :: v1, v2, w
        double precision, intent(in) :: hamOff(nOff), hamDi2(nDi2)
        
        double precision, allocatable, intent(out) :: ham(:,:)

        integer(kind=8)  :: j = 0
        double precision :: rand(sites)
        logical          :: symm 

        if(allocated(ham)) deallocate(ham)
        allocate(ham(dim, dim))
        ham = 0.d0 
        
        call random_number(rand)
        rand = 2 * (rand - 0.5)
        do j = 1, nOff
            ham(rcoff(j,1), rcoff(j,2)) = ham(rcoff(j,1), rcoff(j,2)) + hamOff(j)
        end do 
        do j = 1, nDi !Fill diagonal Hamiltonian
            ham(j,j) = v1 * hamDi(j,1) + v2 * hamDi(j,2) + w * sum(rand * occ(1:sites,j))
        end do
        call test_symm(symm, dim, ham)
        
        do j = 1, nDi2 
            ham(rcdi(j), rcdi(j)) = ham(rcdi(j), rcdi(j)) + hamDi2(j)
        end do 

        return

    end subroutine unify_dense_dp

    subroutine unify_dense_dc(dim, v1, v2, w, sites, occ, nOff, nDi, hamOff, rcoff, hamDi, ham)
        
        implicit none

        integer, intent(in)          :: sites
        integer(kind=8), intent(in)  :: dim, nOff, nDi, occ(sites,*), rcoff(nOff,2), hamDi(nDi,2)
        double precision, intent(in) :: v1, v2, w
        double complex, intent(in)   :: hamOff(nOff)
        
        double complex, allocatable, intent(out) :: ham(:,:)

        integer(kind=8)  :: j = 0
        double precision :: rand(sites)
        logical          :: symm 

        if(allocated(ham)) deallocate(ham)
        allocate(ham(dim, dim))
        ham = 0.d0 
        
        call random_number(rand)
        rand = 2 * (rand - 0.5)
        do j = 1, nOff
            ham(rcoff(j,1), rcoff(j,2)) = ham(rcoff(j,1), rcoff(j,2)) + hamOff(j)
        end do 
        do j = 1, nDi !Fill diagonal Hamiltonian
            ham(j,j) = v1 * hamDi(j,1) + v2 * hamDi(j,2) + w * sum(rand * occ(1:sites,j))
        end do

        return

    end subroutine unify_dense_dc



end module hamiltonian 