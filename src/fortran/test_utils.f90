module test_utils 
    ! Module for testing the functionality of the code
    use file_utils
    implicit none 

    interface check_spec 
        module procedure check_spectrum_dp
        module procedure check_spectrum_dc
    end interface check_spec

    interface check_matrix
        module procedure check_symmetry
        module procedure check_hermiticity
        module procedure check_symmetry_dense
    end interface check_matrix

    contains 

    
    subroutine check_spectrum_dp(dim, nnz, nev, nest, energies, eigstate, ndeg, rc, mat)
        ! Runs sanity checks on the spectrum of eigenvalues and real eigenstates. 
        ! Checks if the eigenstates are normalized, orthogonal, sorted in ascending order and compares the expectation value of the Hamiltonian with the eigenvalue.
        implicit none
        integer(kind=8), intent(in) :: dim, nnz  
        integer, intent(in) :: nev, nest
        integer, intent(in), optional :: ndeg
        integer(kind=8), intent(in), optional :: rc(nnz, 2)
        double precision, intent(inout) :: energies(nev)
        double precision, intent(inout) :: eigstate(dim, nest)
        double precision, intent(in), optional :: mat(nnz)

        integer :: i = 0 
        
        call sort_states_dp(nev, nest, dim, energies, eigstate)
        call check_norm_dp(nev, nest, dim, energies, eigstate)
        call check_orthogonal(nest, dim, eigstate)
        if(present(mat) .and. present(rc)) then 
            do i = 1, max(1, ndeg)
                call expectval_dp(dim, nnz, energies(i), eigstate(1:dim, i), rc, mat)
            end do 
        end if 
        return 
    
    end subroutine check_spectrum_dp

    subroutine check_spectrum_dc(dim, feast, nev, nest, energies, eigstate)
        ! Runs sanity checks on the spectrum of eigenvalues and complex eigenstates. 
        ! Checks if the eigenstates are normalized, orthogonal, and sorted in ascending order.
        implicit none
        integer(kind=8), intent(in) :: dim  
        integer, intent(in) :: nev, nest
        logical, intent(in) :: feast
        double precision, intent(inout) :: energies(nev)
        double complex, intent(inout) :: eigstate(dim, nest)

        call sort_states_dc(nev, nest, dim, energies, eigstate)
        call check_norm_dc(nev, nest, dim, energies, eigstate)
        if(feast) then 
            call check_orthogonal_dc(nest, dim, eigstate)   
        end if 

        return 

    end subroutine check_spectrum_dc

    subroutine sort_states_dp(nev, nest, dim, evals, states)
        ! Sorts the eigenvalues (and corresponding eigenstates) in ascending order, i.e., E_1 < E_2 < ... < E_nev
        implicit none 
        integer, intent(in) :: nev, nest 
        integer(kind=8), intent(in) :: dim 

        double precision, intent(inout) :: evals(nev), states(dim, nev) 

        integer :: i = 0, j = 0, cntr = 0 
        double precision :: temp = 0 
        double precision, allocatable :: temp_vec(:)

        cntr = 0
        if(allocated(temp_vec)) deallocate(temp_vec)
        allocate(temp_vec(dim))
        do i = 1, nest
            do j = 1, nest
                if(dble(evals(j)) < dble(evals(i)) .and. i < j) then 
                    print*,'Swap eigenvalues i, j',i, j
                    temp = evals(i)
                    evals(i) = evals(j)
                    evals(j) = temp 
                    cntr = cntr + 1 
                    if(dble(evals(i)) > dble(evals(j))) print*,'Error in resorting routine.'
                    temp_vec = states(i, 1:nest)
                    states(i, 1:nest) = states(j, 1:nest)
                    states(j, 1:nest) = temp_vec 
                end if 
            end do 
        end do 
        deallocate(temp_vec)
        print*,'Eigenstates sorted.', cntr, ' eigenvalues were in the wrong order.'
        return 

    end subroutine sort_states_dp

    subroutine sort_states_dc(nev, nest, dim, evals, states)
        ! Sorts the eigenvalues (and corresponding eigenstates) in ascending order, i.e., E_1 < E_2 < ... < E_nev
        
        implicit none 
        integer, intent(in) :: nev, nest 
        integer(kind=8), intent(in) :: dim 
        double precision, intent(inout) :: evals(nev)
        double complex, intent(inout) :: states(dim, nev) 

        integer :: i = 0, j = 0, info = 0 
        double complex :: temp = 0 
        double complex, allocatable :: temp_vec(:)

        
        if(allocated(temp_vec)) deallocate(temp_vec)
        allocate(temp_vec(dim))
        do i = 1, nest
            do j = 1, nest
                if(dble(evals(j)) < dble(evals(i)) .and. i < j) then 
                    print*,'Swap eigenvalues i = ', i, 'and j = ', j
                    temp     = evals(i)
                    evals(i) = evals(j)
                    evals(j) = temp  
                    if(dble(evals(i)) > dble(evals(j))) print*,'Error in sort_states_dc.'
                    temp_vec = states(i, 1:nest)
                    states(i, 1:nest) = states(j, 1:nest)
                    states(j, 1:nest) = temp_vec 
                end if 
            end do 
        end do 
        deallocate(temp_vec)

        return 

    end subroutine sort_states_dc

    subroutine expectval_dp(dim, nnz, eval, evec, rc, mat)
        ! Calculates the expectation value of the Hamiltonian for a given eigenstate and checks if the calculated expectation value matches the given eigenvalue.
        ! Uses the coo_to_csr and amux_dp routines to compute the matrix-vector product
        ! and then calculates the dot product with the eigenstate vector
        
        implicit none 
        integer(kind=8), intent(in) :: dim, nnz  
        integer(kind=8), intent(in) :: rc(nnz, 2)
        double precision, intent(in) :: eval
        double precision, intent(in) :: evec(dim)
        double precision, intent(in) :: mat(nnz)
        
        integer(kind=8)  :: ia(dim+1), ja(nnz)
        double precision :: exv, val(nnz), ax(dim)

        ia = 0 
        ja = 0 
        exv = 0.d0 
        val = 0.d0 

        call coo_to_csr(dim, nnz, mat, rc(1:nnz,1), rc(1:nnz,2), val, ja, ia) 
        call spmv(1, dim, evec, ax, val, ja, ia)
        exv = dot_product(evec, ax)
        if(abs(exv - eval) > 1.d-14) then 
            print*,'Error in expectation value calculation.'
            print*,'Expected value:', eval, 'Calculated value:', exv
            error stop 
        end if
        return 

    end subroutine expectval_dp

    subroutine check_norm_dp(nev, nest, dim, evals, states)
        ! Checks if the eigenstates are normalized, i.e., ||psi_i||^2 = 1 for i = 1, ..., nev
        
        implicit none 
        integer, intent(in) :: nev, nest 
        integer(kind=8), intent(in) :: dim 

        double precision, intent(in) :: evals(nev), states(dim, nest) 

        integer :: i = 0 
        double precision :: norm = 0 
        
        do i = 1, nev
            write(* ,"(x, i0, '.Eigenvalue = ',f18.13)") i, dble(evals(i))
            if(i .le. nest) then 
                norm = dot_product(states(1:dim,i),states(1:dim,i))
                if(norm < (1.0-1.d-10)) then
                    write(* ,"(x, i0, '. Norm = ',f16.12)") i, dble(norm)
                    print*,'Error in subroutine check_norm_dp.'
                    error stop  
                end if 
            end if 
        end do     
        print*, 'All eigenstates have norm 1.'
        return 

    end subroutine check_norm_dp

    subroutine check_norm_dc(nev, nest, dim, evals, states)
        ! Checks if the eigenstates are normalized, i.e., ||psi_i||^2 = 1 for i = 1, ..., nev

        implicit none 
        integer, intent(in) :: nev, nest 
        integer(kind=8), intent(in) :: dim 
        double precision, intent(in) :: evals(nev)
        double complex, intent(in) :: states(dim, nest) 

        integer :: i = 0 
        double precision :: norm = 0 
        double complex :: psi(dim)
        
        psi = 0.d0  

        do i = 1, nev
            write(* ,"(x, i0, '.Eigenvalue = ',f18.13)") i, dble(evals(i))
            if(i .le. nest) then 
                norm = dot_product(states(1:dim,i), states(1:dim,i))           
                
                if(norm < (1.0-1.d-10)) then
                    write(* ,"(x, i0, '. Norm = ',f16.12)") i, dble(norm)
                    print*,'Error in subroutine check_norm_dc.'
                    error stop
                end if 
            end if 
        end do     
        print*,'All eigenstates are normalized.'
        return 
        
    end subroutine check_norm_dc    

    subroutine check_orthogonal(nest, dim, states)
        ! Checks if the eigenstates are orthogonal, i.e., <psi_i|psi_j> = 0 for i != j
        implicit none 
        integer, intent(in) :: nest 
        integer(kind=8), intent(in) :: dim 

        double precision, intent(in) :: states(dim, nest) 

        integer :: i = 0, j = 0, info = 0
        
        
        do i = 1, nest 
            do j = 1, nest 
                if(abs(dot_product(states(1:dim, i), states(1:dim, j))) < (1.0-1.d-10) .and. dot_product(states(1:dim, i), states(1:dim, j)) > 1.d-10) then
                    print*, 'Error in check_orthogonal.'
                    print*, 'Eigenstates', i, j, 'are not orthogonal.'
                    print*, 'Psi(i).Psi(j)', dot_product(states(1:dim, i), states(1:dim, j))
                    info = 1
                end if 
            end do 
        end do 
        if(info == 1) then 
            error stop 
        else
            print*,'All eigenstates are orthorgonal.'
        end if 

        return 

    end subroutine check_orthogonal

    subroutine check_orthogonal_dc(nest, dim, states)
        ! Checks if the eigenstates are orthogonal, i.e., <psi_i|psi_j> = 0 for i != j
        implicit none 
        integer, intent(in) :: nest 
        integer(kind=8), intent(in) :: dim 

        double complex, intent(in) :: states(dim, nest) 

        integer :: i = 0, j = 0, info = 0
        integer(kind=8) :: k = 0
        double complex :: res = 0.d0 
        double complex, allocatable :: vec1(:), vec2(:)

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
                if(abs(dot_product(states(1:dim, i), states(1:dim, j))) < (1.0-1.d-10) .and. abs(dot_product(states(1:dim, i), states(1:dim, j))) > 1.d-10) then
                    print*, 'Error in check_orthogonal_dc.'
                    print*, 'Eigenstates', i, j, 'are not orthogonal.'
                    print*, 'Psi(i).Psi(j)', i, j, dot_product(states(1:dim, i), states(1:dim, j))
                    info = 1
                end if 
            end do 
        end do 
        if(info == 1) then 
            error stop 
        else
            print*,'All eigenstates are orthorgonal.'
        end if

        return 

    end subroutine check_orthogonal_dc


    subroutine check_symmetry(symmetric, dim, nz, rc, ham)

        implicit none

        integer(kind=8),  intent(in) :: dim
        integer(kind=8),  intent(in) :: nz 
        integer(kind=8),  intent(in) :: rc(nz,2)
        double precision, intent(in) :: ham(nz)
        
        logical, intent(out) :: symmetric
        
        integer(kind=8)  :: i = 0, j = 0
        double precision :: eps = 10.**(-10)
        double precision, allocatable :: mat(:,:)
        
        if(allocated(mat))    deallocate(mat)
        allocate(mat(dim, dim))

        mat = 0
        do i = 1, nz 
            mat(rc(i,1),rc(i,2)) = ham(i) 
        end do 
    
        symmetric=.true.
        do i = 1, dim
            do j = 1, dim
                if(abs(mat(i,j) - mat(j,i)) > eps) then
                    print*, i, j, 'i, j'
                    symmetric = .false.
                    exit
                endif
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

    end subroutine check_symmetry

    subroutine check_hermiticity(hermitian, irrep, dim, nz, rc, ham)

        implicit none

        integer(kind=8),             intent(in) :: dim, nz
        integer(kind=8),             intent(in) :: rc(nz,2)
        double complex, allocatable, intent(in) :: ham(:) 
        character,                   intent(in) :: irrep*2
        logical,                    intent(out) :: hermitian
        
        integer(kind=8)             :: i, j, mdim
        double precision            :: eps = 10.**(-10)
        double complex, allocatable :: mat(:,:)

        if(irrep(1:1) .ne. "E") then 
            mdim = dim 
        else 
            mdim = 2 * dim 
        end if
        if(allocated(mat)) deallocate(mat)
        allocate(mat(mdim, mdim))
        mat = 0


        do i = 1, nz 
            mat(rc(i,1),rc(i,2)) = ham(i)
        end do 

        
        hermitian = .true.
        do i = 1, mdim
            do j = 1, mdim
                if(dble(conjg(mat(i,j))) - dble(mat(j,i)) > eps .or. aimag(conjg(mat(i,j))) - aimag(mat(j,i)) > eps) then
                    if(dble(conjg(mat(i,j))) - dble(mat(j,i)) > eps) then
                        print*, 'RE(M*(i,j))', dble(conjg(mat(i,j)))
                        print*, 'RE(M(j,i))', dble(mat(j,i))
                        print*, 'Difference', abs(dble(conjg(mat(i,j))) - dble(mat(j,i)))
                        print*, ''
                    end if

                    if(aimag(conjg(mat(i,j))) - aimag(mat(j,i)) > eps) then
                        print*, 'IM(M*(i,j))', aimag(conjg(mat(i,j)))
                        print*, 'IM(M(j,i))', aimag(mat(j,i))
                        print*, 'Difference', abs(aimag(conjg(mat(i,j))) - aimag(mat(j,i)))
                        print*, ''
                    end if

                    print*, 'i, j', i, j
                    hermitian = .false.
                    
                endif
            end do
        end do
        if(.not. hermitian) then
            print*, 'MATRIX IS NON-HERMITIAN!'
            print*, ''
            error stop
        else
            print*, 'MATRIX IS HERMITIAN'
            print*, ''
        end if
        return

    end subroutine check_hermiticity

    subroutine check_symmetry_dense(sym, dim, mat)

        implicit none

        integer(kind=8),  intent(in) :: dim
        double precision, intent(in) :: mat(dim, dim)
        double precision :: eps = 0.000001 
        logical, intent(out) :: sym

        integer(kind=8) :: i, j


        sym=.true.
        do i = 1, dim
            do j = 1, dim
                if(abs(mat(i,j) - mat(j,i)) > eps) then
                    sym=.false.
                    print*, i,j, 'i,j'
                endif
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

    end subroutine check_symmetry_dense

end module test_utils 
