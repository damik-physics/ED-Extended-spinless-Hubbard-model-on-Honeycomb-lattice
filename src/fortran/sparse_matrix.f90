module sparse_matrix
    ! Sparse matrix module with CSR storage and Lanczos diagonalization
    ! Implements iterative eigenvalue solver for large sparse Hamiltonian matrices
    
    implicit none
    
    ! Sparse matrix type in Compressed Sparse Row (CSR) format
    type :: sparse_csr
        integer :: n                              ! Matrix dimension
        integer :: nnz                            ! Number of non-zero elements
        double precision, allocatable :: values(:) ! Non-zero values
        integer, allocatable :: col_indices(:)     ! Column indices
        integer, allocatable :: row_ptr(:)         ! Row pointers
    end type sparse_csr
    
    ! Lanczos iteration data
    type :: lanczos_data
        integer :: max_iter                       ! Maximum Lanczos iterations
        integer :: n_evals                        ! Number of eigenvalues requested
        double precision :: tolerance             ! Convergence tolerance
        double precision, allocatable :: alpha(:) ! Diagonal elements of tridiagonal matrix
        double precision, allocatable :: beta(:)  ! Off-diagonal elements
        double precision, allocatable :: q(:,:)   ! Lanczos vectors
        logical :: converged                      ! Convergence flag
    end type lanczos_data
    
contains

    subroutine sparse_csr_init(matrix, n, nnz)
        ! Initialize sparse CSR matrix structure
        implicit none
        type(sparse_csr), intent(out) :: matrix
        integer, intent(in) :: n, nnz
        
        matrix%n = n
        matrix%nnz = nnz
        
        if (allocated(matrix%values)) deallocate(matrix%values)
        if (allocated(matrix%col_indices)) deallocate(matrix%col_indices)
        if (allocated(matrix%row_ptr)) deallocate(matrix%row_ptr)
        
        allocate(matrix%values(nnz))
        allocate(matrix%col_indices(nnz))
        allocate(matrix%row_ptr(n + 1))
        
        matrix%values = 0.0d0
        matrix%col_indices = 0
        matrix%row_ptr = 0
    end subroutine sparse_csr_init

    subroutine dense_to_sparse_csr(dense_matrix, sparse_matrix, threshold)
        ! Convert dense matrix to sparse CSR format
        implicit none
        double precision, intent(in) :: dense_matrix(:,:)
        type(sparse_csr), intent(out) :: sparse_matrix
        double precision, intent(in), optional :: threshold
        
        integer :: n, i, j, nnz, idx
        double precision :: thresh
        
        n = size(dense_matrix, 1)
        thresh = 1.0d-12
        if (present(threshold)) thresh = threshold
        
        ! Count non-zero elements
        nnz = 0
        do i = 1, n
            do j = 1, n
                if (abs(dense_matrix(i, j)) > thresh) nnz = nnz + 1
            end do
        end do
        
        call sparse_csr_init(sparse_matrix, n, nnz)
        
        ! Fill CSR format
        idx = 1
        sparse_matrix%row_ptr(1) = 1
        
        do i = 1, n
            do j = 1, n
                if (abs(dense_matrix(i, j)) > thresh) then
                    sparse_matrix%values(idx) = dense_matrix(i, j)
                    sparse_matrix%col_indices(idx) = j
                    idx = idx + 1
                end if
            end do
            sparse_matrix%row_ptr(i + 1) = idx
        end do
        
        write(*,'(A,I0,A,I0,A,F6.2,A)') 'Converted to sparse CSR: ', n, 'x', n, &
            ' matrix with ', real(nnz)/real(n*n)*100.0, '% non-zeros'
    end subroutine dense_to_sparse_csr

    subroutine sparse_matvec(matrix, x, y)
        ! Sparse matrix-vector multiplication: y = A*x
        implicit none
        type(sparse_csr), intent(in) :: matrix
        double precision, intent(in) :: x(:)
        double precision, intent(out) :: y(:)
        
        integer :: i, j
        
        y = 0.0d0
        
        do i = 1, matrix%n
            do j = matrix%row_ptr(i), matrix%row_ptr(i + 1) - 1
                y(i) = y(i) + matrix%values(j) * x(matrix%col_indices(j))
            end do
        end do
    end subroutine sparse_matvec

    subroutine lanczos_init(lanczos, n, n_evals, max_iter, tolerance)
        ! Initialize Lanczos iteration data structure
        implicit none
        type(lanczos_data), intent(out) :: lanczos
        integer, intent(in) :: n, n_evals, max_iter
        double precision, intent(in) :: tolerance
        
        lanczos%max_iter = min(max_iter, n)  ! Cannot exceed matrix dimension
        lanczos%n_evals = n_evals
        lanczos%tolerance = tolerance
        lanczos%converged = .false.
        
        if (allocated(lanczos%alpha)) deallocate(lanczos%alpha)
        if (allocated(lanczos%beta)) deallocate(lanczos%beta)
        if (allocated(lanczos%q)) deallocate(lanczos%q)
        
        allocate(lanczos%alpha(lanczos%max_iter))
        allocate(lanczos%beta(lanczos%max_iter))
        allocate(lanczos%q(n, lanczos%max_iter))
        
        lanczos%alpha = 0.0d0
        lanczos%beta = 0.0d0
        lanczos%q = 0.0d0
    end subroutine lanczos_init

    subroutine lanczos_diagonalize(matrix, eigenvalues, eigenvectors, n_evals, max_iter, tolerance)
        ! Main Lanczos diagonalization routine
        implicit none
        type(sparse_csr), intent(in) :: matrix
        double precision, allocatable, intent(out) :: eigenvalues(:)
        double precision, allocatable, intent(out) :: eigenvectors(:,:)
        integer, intent(in) :: n_evals, max_iter
        double precision, intent(in) :: tolerance
        
        type(lanczos_data) :: lanczos
        double precision, allocatable :: r(:), w(:), tridiag(:,:), work(:)
        double precision, allocatable :: tridiag_evals(:), tridiag_evecs(:,:)
        integer :: n, iter, info, lwork, j
        double precision :: norm_r, beta_old
        
        n = matrix%n
        call lanczos_init(lanczos, n, n_evals, max_iter, tolerance)
        
        allocate(r(n), w(n))
        allocate(tridiag(lanczos%max_iter, lanczos%max_iter))
        allocate(tridiag_evals(lanczos%max_iter))
        allocate(tridiag_evecs(lanczos%max_iter, lanczos%max_iter))
        
        ! Initialize with random vector
        call random_seed()
        call random_number(lanczos%q(:, 1))
        lanczos%q(:, 1) = lanczos%q(:, 1) / norm2(lanczos%q(:, 1))
        
        write(*,*) 'Starting Lanczos iterations...'
        
        ! Lanczos iterations
        do iter = 1, lanczos%max_iter
            ! Matrix-vector multiplication
            call sparse_matvec(matrix, lanczos%q(:, iter), w)
            
            ! Compute alpha_i = q_i^T * A * q_i
            lanczos%alpha(iter) = dot_product(lanczos%q(:, iter), w)
            
            ! Compute residual: r = A*q_i - alpha_i*q_i - beta_{i-1}*q_{i-1}
            r = w - lanczos%alpha(iter) * lanczos%q(:, iter)
            if (iter > 1) then
                r = r - lanczos%beta(iter-1) * lanczos%q(:, iter-1)
            end if
            
            ! Compute beta_i = ||r||
            norm_r = norm2(r)
            
            if (iter < lanczos%max_iter) then
                lanczos%beta(iter) = norm_r
                
                ! Check for breakdown
                if (norm_r < 1.0d-14) then
                    write(*,'(A,I0)') 'Lanczos breakdown at iteration ', iter
                    exit
                end if
                
                ! Next Lanczos vector: q_{i+1} = r / beta_i
                lanczos%q(:, iter + 1) = r / norm_r
                
                ! Reorthogonalization (full reorthogonalization)
                call reorthogonalize(lanczos%q(:, 1:iter+1), iter + 1)
            end if
            
            ! Check convergence every few iterations
            if (mod(iter, 10) == 0 .or. iter == lanczos%max_iter) then
                call check_lanczos_convergence(lanczos, iter, n_evals, tolerance)
                if (lanczos%converged) then
                    write(*,'(A,I0,A)') 'Lanczos converged after ', iter, ' iterations'
                    exit
                end if
            end if
        end do
        
        ! Diagonalize tridiagonal matrix
        tridiag = 0.0d0
        do iter = 1, lanczos%max_iter
            if (iter > 0) tridiag(iter, iter) = lanczos%alpha(iter)
            if (iter < lanczos%max_iter .and. lanczos%beta(iter) > 0) then
                tridiag(iter, iter+1) = lanczos%beta(iter)
                tridiag(iter+1, iter) = lanczos%beta(iter)
            end if
        end do
        
        ! Use LAPACK to diagonalize tridiagonal matrix
        tridiag_evecs = tridiag
        lwork = -1
        allocate(work(1))
        call dsyev('V', 'U', lanczos%max_iter, tridiag_evecs, lanczos%max_iter, &
                   tridiag_evals, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        
        call dsyev('V', 'U', lanczos%max_iter, tridiag_evecs, lanczos%max_iter, &
                   tridiag_evals, work, lwork, info)
        
        if (info /= 0) then
            write(*,*) 'ERROR: Tridiagonal diagonalization failed'
            stop
        end if
        
        ! Extract requested eigenvalues and transform eigenvectors
        allocate(eigenvalues(n_evals))
        allocate(eigenvectors(n, n_evals))
        
        eigenvalues(1:n_evals) = tridiag_evals(1:n_evals)
        
        ! Transform eigenvectors: eigenvectors = Q * tridiag_evecs
        do iter = 1, n_evals
            eigenvectors(:, iter) = 0.0d0
            do j = 1, lanczos%max_iter
                eigenvectors(:, iter) = eigenvectors(:, iter) + &
                    lanczos%q(:, j) * tridiag_evecs(j, iter)
            end do
        end do
        
        write(*,'(A,I0,A)') 'Computed ', n_evals, ' eigenvalues using Lanczos'
        
        ! Cleanup
        deallocate(r, w, tridiag, tridiag_evals, tridiag_evecs, work)
    end subroutine lanczos_diagonalize

    subroutine reorthogonalize(q_vectors, n_vecs)
        ! Modified Gram-Schmidt reorthogonalization
        implicit none
        double precision, intent(inout) :: q_vectors(:,:)
        integer, intent(in) :: n_vecs
        
        integer :: i, j
        double precision :: projection
        
        do i = 1, n_vecs - 1
            projection = dot_product(q_vectors(:, n_vecs), q_vectors(:, i))
            q_vectors(:, n_vecs) = q_vectors(:, n_vecs) - projection * q_vectors(:, i)
        end do
        
        ! Normalize
        q_vectors(:, n_vecs) = q_vectors(:, n_vecs) / norm2(q_vectors(:, n_vecs))
    end subroutine reorthogonalize

    subroutine check_lanczos_convergence(lanczos, current_iter, n_evals, tolerance)
        ! Check convergence of Lanczos iteration
        implicit none
        type(lanczos_data), intent(inout) :: lanczos
        integer, intent(in) :: current_iter, n_evals
        double precision, intent(in) :: tolerance
        
        double precision, allocatable :: tridiag(:,:), evals(:), work(:)
        integer :: info, lwork, i
        
        if (current_iter < n_evals + 5) return  ! Need enough iterations
        
        ! Quick convergence check using tridiagonal eigenvalues
        allocate(tridiag(current_iter, current_iter))
        allocate(evals(current_iter))
        
        tridiag = 0.0d0
        do i = 1, current_iter
            tridiag(i, i) = lanczos%alpha(i)
            if (i < current_iter) then
                tridiag(i, i+1) = lanczos%beta(i)
                tridiag(i+1, i) = lanczos%beta(i)
            end if
        end do
        
        ! Diagonalize to check convergence
        lwork = -1
        allocate(work(1))
        call dsyev('N', 'U', current_iter, tridiag, current_iter, evals, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        call dsyev('N', 'U', current_iter, tridiag, current_iter, evals, work, lwork, info)
        
        ! Simple convergence criterion: beta should be small compared to eigenvalue gaps
        if (current_iter > 1 .and. lanczos%beta(current_iter-1) < tolerance * abs(evals(1))) then
            lanczos%converged = .true.
        end if
        
        deallocate(tridiag, evals, work)
    end subroutine check_lanczos_convergence

    subroutine sparse_csr_destroy(matrix)
        ! Cleanup sparse CSR matrix
        implicit none
        type(sparse_csr), intent(inout) :: matrix
        
        if (allocated(matrix%values)) deallocate(matrix%values)
        if (allocated(matrix%col_indices)) deallocate(matrix%col_indices)
        if (allocated(matrix%row_ptr)) deallocate(matrix%row_ptr)
        
        matrix%n = 0
        matrix%nnz = 0
    end subroutine sparse_csr_destroy

    subroutine lanczos_destroy(lanczos)
        ! Cleanup Lanczos data structure
        implicit none
        type(lanczos_data), intent(inout) :: lanczos
        
        if (allocated(lanczos%alpha)) deallocate(lanczos%alpha)
        if (allocated(lanczos%beta)) deallocate(lanczos%beta)
        if (allocated(lanczos%q)) deallocate(lanczos%q)
    end subroutine lanczos_destroy

end module sparse_matrix