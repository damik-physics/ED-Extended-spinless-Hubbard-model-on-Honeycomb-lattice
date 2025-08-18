module diagonalization
    use functions
    use test_utils
    use file_utils, only: spmv        
    use types
    use params
    use io_utils



    implicit none 


    interface exactdiag
        module procedure exactdiag_dp, exactdiag_dc
    end interface exactdiag

    interface unify 
            module procedure unify_i8, unify_dp, unify_dc, unify_dc_2D, unifyDense_i8, unifyDense_dp, unifyDense_dc 
    end interface unify

    interface gsdeg
        module procedure gsdeg_dp, gsdeg_dc, qgsdeg_dp, qgsdeg_dc
    end interface gsdeg
    

    contains 

    !----------------------------------------------!
    !               Diagonalization                !
    !----------------------------------------------!

    subroutine diagonalize(par, geo, ham, diag, st, out)
        
        implicit none 
        type(sim_params), intent(inout) :: par
        type(geometry), intent(inout) :: geo
        type(hamiltonian_params), intent(inout) :: ham
        type(diag_params), intent(inout) :: diag
        type(system_state), intent(inout) :: st
        type(output), intent(inout) :: out

        integer :: i = 0 
        double precision :: Emin = 0.d0, Emax = 0.d0
        logical :: symmetric, append 

        if((par%nDis > 1) .and.(par%dis .ne. 0.d0) .and.(st%conf > 1)) then 
            append = .true.
        else 
            append = .false.
        end if 

        100 format(1000(x,A,x,F6.4,x,A,x,F6.4))
        
        ! Determine number of eigenvalues (nev) and number of Lanczos vectors (ncv)
        call tune_lanczos_parameters(out%unit, par%dimthresh, par%exact, par%nevext, par%nst, par%ncv0, geo%dim, diag%full, diag%nev, diag%ncv, diag%nest)

        if(diag%full == 0) then
            if(st%mat_type== "R") then 
                if(par%ti == 0) then 
                    call unify(par%t, st%v1, st%v2, par%dis, par%mass, pattern, geo%sites, ham%occ, ham%noff, ham%ndi, ham%hamoff, ham%hamdi, ham%ham, ham%rc, ham%nnz)
                    if(par%arpack == 1) then 
                        if(par%otf == 0) then
                            call lanczos_dp(par%othrds, geo%dim, ham%nnz, diag%nev, diag%ncv, diag%nest, mode, par%rvec, ham%ham, ham%rc, diag%energies, diag%eigstate)          
                        else if(par%otf == 1) then 
                            call lanczos_dp_otf(out%unit, par%othrds, geo%dim, geo%sites, geo%sites/2*3, geo%sites*3, diag%nev, diag%ncv, diag%nest, par%t, st%v1, st%v2, mode, par%rvec, geo%basis_states, geo%bsites, geo%hexsites, diag%energies, diag%eigstate)
                        end if 
                        if(par%rvec) call check_spec(geo%dim, ham%nnz, diag%nev, diag%nest, diag%energies, diag%eigstate)        
                        ! call save(out%outdir, "R", append, st%conf, out%unit, geo%dim, par%states, diag%nev, diag%nest, par%rvec, diag%energies, diag%eigstate)
                    end if 
                    if(par%mkl == 1) then 
                        call lanczos_mkl(geo%dim, diag%nev, diag%ncv, diag%nest, ham%nnz, ham%ham, ham%rc, diag%energies, diag%eigstate) 
                        if(par%rvec) call check_spec(geo%dim, ham%nnz, diag%nev, diag%nest, diag%energies, diag%eigstate)                  
                        ! call save(out%outdir, "R", append, st%conf, out%unit, geo%dim, par%states, diag%nev, diag%nest, par%rvec, diag%energies, diag%eigstate)
                    end if 
                    if(par%degflag == 1) call gsdeg(out%outdir, par%rvec, out%unit, geo%dim, diag%nev, diag%nest, par%deg, diag%energies, diag%eigstate, diag%ndeg, diag%gs)     
                    if(par%degflag == 2 .and. par%rvec) call gsdeg(out%outdir, out%unit, geo%dim, diag%nev, diag%nest, diag%energies, diag%eigstate, diag%ndeg, diag%gs)
                else if(par%ti == 1) then 
                    call unify(st%v1, st%v2, par%dis, par%mass, pattern, geo%sites, ham%occ, ham%noff, ham%ndi, ham%ndioff, ham%hamoff_dp, ham%rcoff, ham%rcdi, ham%hamdi_dp, ham%hamdi, ham%ham, ham%rc, ham%nnz)
                    call check_matrix(symmetric, geo%dim, ham%nnz, ham%rc, ham%ham)
                    if(par%arpack == 1) then 
                        call lanczos_dp(par%othrds, geo%dim, ham%nnz, diag%nev, diag%ncv, diag%nest, mode, par%rvec, ham%ham, ham%rc, diag%energies, diag%eigstate)
                        if(par%degflag == 1) call gsdeg(out%outdir, par%rvec, out%unit, geo%dim, diag%nev, diag%nest, par%deg, diag%energies, diag%eigstate, diag%ndeg, diag%gs)
                        if(par%degflag == 2 .and. par%rvec) call gsdeg(out%outdir, out%unit, geo%dim, diag%nev, diag%nest, diag%energies, diag%eigstate, diag%ndeg, diag%gs)
                        if(par%rvec) call check_spec(geo%dim, ham%nnz, diag%nev, diag%nest, diag%energies, diag%eigstate, ndeg = diag%ndeg, rc = ham%rc, mat = ham%ham) 
                        ! call save(out%outdir, "R", append, st%conf, out%unit, geo%dim, par%states, diag%nev, diag%nest, par%rvec, diag%energies, diag%eigstate)    
                        if(par%degflag == 1) call gsdeg(out%outdir, par%rvec, out%unit, geo%dim, diag%nev, diag%nest, par%deg, diag%energies, diag%eigstate, diag%ndeg, diag%gs)
                        if(par%degflag == 2 .and. par%rvec) call gsdeg(out%outdir, out%unit, geo%dim, diag%nev, diag%nest, diag%energies, diag%eigstate, diag%ndeg, diag%gs)
                    end if 
                    if(par%mkl == 1) then 
                        call lanczos_mkl(geo%dim, diag%nev, diag%ncv, diag%nest, ham%nnz, ham%ham, ham%rc, diag%energies, diag%eigstate)                
                        if(par%degflag == 1) call gsdeg(out%outdir, par%rvec, out%unit, geo%dim, diag%nev, diag%nest, par%deg, diag%energies, diag%eigstate, diag%ndeg, diag%gs)
                        if(par%degflag == 2 .and. par%rvec) call gsdeg(out%outdir, out%unit, geo%dim, diag%nev, diag%nest, diag%energies, diag%eigstate, diag%ndeg, diag%gs)
                        if(par%rvec) call check_spec(geo%dim, ham%nnz, diag%nev, diag%nest, diag%energies, diag%eigstate, ndeg = diag%ndeg, rc = ham%rc, mat = ham%ham) 
                        ! call save(out%outdir, "R", append, st%conf, out%unit, geo%dim, par%states, diag%nev, diag%nest, par%rvec, diag%energies, diag%eigstate)
                        if(par%degflag == 1) call gsdeg(out%outdir, par%rvec, out%unit, geo%dim, diag%nev, diag%nest, par%deg, diag%energies, diag%eigstate, diag%ndeg, diag%gs)     
                        if(par%degflag == 2 .and. par%rvec) call gsdeg(out%outdir, out%unit, geo%dim, diag%nev, diag%nest, diag%energies, diag%eigstate, diag%ndeg, diag%gs)    
                    end if 
                end if
                if(par%feast == 1) then 
                    Emin = dble(diag%energies(1)) - 0.01
                    Emax = dble(diag%energies(min(diag%nev, par%nevmax)))+0.01
                    
                    call feast_dp(geo%dim, ham%nnz, par%nev0, diag%nest, ham%rc, ham%ham, Emin, Emax, diag%nev, diag%energies, diag%eigstate)

                    if(par%degflag == 1) call gsdeg(out%outdir, par%rvec, out%unit, geo%dim, diag%nev, diag%nest, par%deg, diag%energies, diag%eigstate, diag%ndeg, diag%gs)
                    if(par%degflag == 2 .and. par%rvec) call gsdeg(out%outdir, out%unit, geo%dim, diag%nev, diag%nest, diag%energies, diag%eigstate, diag%ndeg, diag%gs)
                    if(par%rvec) call check_spec(geo%dim, ham%nnz, diag%nev, diag%nest, diag%energies, diag%eigstate, ndeg = diag%ndeg, rc = ham%rc, mat = ham%ham)
                    ! call save(out%outdir, "R", append, st%conf, out%unit, geo%dim, par%states, diag%nev, diag%nest, par%rvec, diag%energies, diag%eigstate)
                    if(par%degflag == 1) call gsdeg(out%outdir, par%rvec, out%unit, geo%dim, diag%nev, diag%nest, par%deg, diag%energies, diag%eigstate, diag%ndeg, diag%gs)
                    if(par%degflag == 2 .and. par%rvec) call gsdeg(out%outdir, out%unit, geo%dim, diag%nev, diag%nest, diag%energies, diag%eigstate, diag%ndeg, diag%gs)
                end if                          
            else if(par%ti == 1 .and. st%mat_type== "C") then 
                if(par%irrep(1:1) == "E") then 
                    call unify(geo%sites, ham%noff, ham%ndioff, geo%dim, st%v1, st%v2, par%dis, par%mass, pattern, ham%occ, ham%hamoff_dc, ham%rcoff, ham%hamdi_off_dc, ham%rcdi, ham%hamDi_dc, ham%ham_dc, ham%rc, ham%nnz)
                else
                    call unify(st%v1, st%v2, par%dis, par%mass, pattern, geo%sites, ham%occ, ham%noff, ham%ndi, ham%hamoff_dc, ham%rcoff, ham%hamdi, ham%ham_dc, ham%rc, ham%nnz)   
                end if 
                call check_matrix(symmetric, par%irrep, geo%dim, ham%nnz, ham%rc, ham%ham_dc)
                call lanczos_dc(par%othrds, geo%dim, diag%nev, diag%ncv, diag%nest, mode, par%rvec, ham%nnz, ham%ham_dc, ham%rc, diag%energies, diag%eigstate_dc)
                if(par%rvec) call check_spec(geo%dim, .false., diag%nev, diag%nest, diag%energies, diag%eigstate_dc)
                ! call save(out%outdir, "C", append, st%conf, out%unit, geo%dim, par%states, diag%nev, diag%nest, par%rvec, diag%energies, diag%eigstate_dc)          
                if(par%degflag == 1) call gsdeg(out%outdir, par%rvec, out%unit, geo%dim, diag%nev, diag%nest, par%deg, diag%energies, diag%eigstate_dc, diag%ndeg, diag%gs_dc)      
                if(par%degflag == 2 .and. par%rvec) call qgsdeg_dc(out%outdir, out%unit, geo%dim, diag%nev, diag%nest, diag%energies, diag%eigstate_dc, diag%ndeg, diag%gs_dc)
                if(par%feast == 1) then 
                    Emin = dble(diag%energies(1))-0.01      
                    Emax = dble(diag%energies(min(diag%nev, par%nevmax)))+0.01
                    call feast_dc(geo%dim, ham%nnz, par%nev0, diag%nest, ham%rc, ham%ham_dc, Emin, Emax, diag%nev, diag%energies, diag%eigstate_dc)
                    if(par%rvec) call check_spec(geo%dim, .True., diag%nev, diag%nest, diag%energies, diag%eigstate_dc)
                    ! call save(out%outdir, "C", append, st%conf, out%unit, geo%dim, par%states, diag%nev, diag%nest, par%rvec, diag%energies, diag%eigstate_dc)           
                    if(par%degflag == 1) call gsdeg(out%outdir, par%rvec, out%unit, geo%dim, diag%nev, diag%nest, par%deg, diag%energies, diag%eigstate_dc, diag%ndeg, diag%gs_dc)            
                    if(par%degflag == 2 .and. par%rvec) call qgsdeg_dc(out%outdir, out%unit, geo%dim, diag%nev, diag%nest, diag%energies, diag%eigstate_dc, diag%ndeg, diag%gs_dc)                      
                end if 
            end if
        else if(diag%full == 1) then
            if(par%ti == 0 .and. st%mat_type== "R") then 
                call unify(geo%dim, par%t, st%v1, st%v2, par%dis, geo%sites, ham%occ, ham%noff, ham%ndi, ham%hamoff, ham%hamdi, ham%hamDense_dp)
                call check_matrix(symmetric, geo%dim, ham%hamDense_dp)
                call exactdiag(par%rvec, geo%dim, ham%hamDense_dp, diag%energies)    
            else if(par%ti == 1 .and. st%mat_type== "R") then 
                call unify(geo%dim, st%v1, st%v2, par%dis, geo%sites, ham%occ, ham%noff, ham%ndi, ham%ndioff, ham%hamoff_dp, ham%rcoff, ham%rcdi, ham%hamdi_dp, ham%hamdi, ham%hamDense_dp)
                call check_matrix(symmetric, geo%dim, ham%hamDense_dp)
                call exactdiag(par%rvec, geo%dim, ham%hamDense_dp, diag%energies)
            else if(par%ti == 1 .and. st%mat_type== "C") then 
                call unify(geo%dim, st%v1, st%v2, par%dis, geo%sites, ham%occ, ham%noff, ham%ndi, ham%hamoff_dc, ham%rcoff, ham%hamdi, ham%hamDense_dc)
                call exactdiag(par%rvec, geo%dim, ham%hamDense_dc, diag%energies)
            end if 
            if(st%mat_type== "R") then 
                if(par%rvec) then 
                    if(allocated(diag%eigstate)) deallocate(diag%eigstate)
                    allocate(diag%eigstate(geo%dim, diag%nest))
                    diag%eigstate = 0.d0 
                    do i = 1, diag%nest
                        diag%eigstate(1:geo%dim,i) = ham%hamDense_dp(1:geo%dim,i)
                    end do
                    if(par%degflag == 1) call gsdeg(out%outdir, par%rvec, out%unit, geo%dim, diag%nev, diag%nest, par%deg, diag%energies, diag%eigstate, diag%ndeg, diag%gs)
                    if(par%degflag == 2 .and. par%rvec) call gsdeg(out%outdir, out%unit, geo%dim, diag%nev, diag%nest, diag%energies, diag%eigstate, diag%ndeg, diag%gs)
                    if(par%rvec) call check_spec(geo%dim, ham%nnz, diag%nev, diag%nest, diag%energies, diag%eigstate, ndeg = diag%ndeg)
                end if
                ! call save(out%outdir, "R", append, st%conf, out%unit, geo%dim, par%states, diag%nev, diag%nest, par%rvec, diag%energies, diag%eigstate)
            else if(st%mat_type== "C") then 
                if(par%rvec) then 
                    if(allocated(diag%eigstate_dc)) deallocate(diag%eigstate_dc)
                    allocate(diag%eigstate_dc(geo%dim, diag%nest))
                    diag%eigstate_dc = 0.d0 
                    do i = 1, diag%nest
                        diag%eigstate_dc(1:geo%dim,i) = ham%hamDense_dc(1:geo%dim,i)
                    end do
                    call gsdeg(out%outdir, par%rvec, out%unit, geo%dim, diag%nev, diag%nest, par%deg, diag%energies, diag%eigstate_dc, diag%ndeg, diag%gs_dc)
                    if(par%rvec) call check_spec(geo%dim, .True., diag%nev, diag%nest, diag%energies, diag%eigstate_dc)
                end if
                ! call save(out%outdir, "C", append, st%conf, out%unit, geo%dim, par%states, diag%nev, diag%nest, par%rvec, diag%energies, diag%eigstate_dc)
            end if 
        end if
        
        call append_energy_csv(out%outdir, par%ti, par%ndis, st%k1, st%k2, st%conf, st%v1, st%v2, diag%energies, st%first_call)
        if(st%mat_type == "R" .and. par%rvec .and. par%states == 1) then 
            call append_states_csv_dp(out%outdir, par%ti, par%ndis, st%k1, st%k2, st%conf, st%v1, st%v2, diag%eigstate)
        else if(st%mat_type == "C" .and. par%rvec .and. par%states == 1) then
            call append_states_csv_dc(out%outdir, par%ti, par%ndis, st%k1, st%k2, st%conf, st%v1, st%v2, diag%eigstate_dc)
        end if 

    end subroutine diagonalize

    !-------------------------------------------!
    !            Exact diagonalization          !
    !-------------------------------------------!

    subroutine exactdiag_dp(eigvec, dim, mat, evals)
    
        implicit none
        logical,                       intent(in)    :: eigvec
        integer(kind=8),               intent(in)    :: dim
        double precision,              intent(inout) :: mat(dim, dim)
        double precision, allocatable, intent(out)   :: evals(:)

        integer, parameter :: lwmax = 100000
        integer :: lda = 0, info = 0, lwork = 0
        character :: jobz
        double precision, allocatable :: matloc(:,:), evalsloc(:), work(:)

        external :: dsyev, dsyev_2stage

        !$OMP THREADPRIVATE(matloc, evalsloc, lda, info, lwork, jobz, work)
        lda   = 0
        info  = 0
        lwork = 0

        if(eigvec) then
            jobz = 'V'
        else
            jobz = 'N'
        end if

        if(allocated(evals)) deallocate(evals)
        allocate(evals(dim))
        evals = 0.0d0

        if(allocated(evalsloc)) deallocate(evalsloc)
        if(allocated(matloc)) deallocate(matloc)
        allocate(evalsloc(dim))
        allocate(matloc(dim, dim))
        evalsloc = 0.0d0
        matloc = mat

        if(dim > 1) then
            lda = dim
            if(allocated(work)) deallocate(work)
            allocate(work(lwmax))
            lwork = -1
            if(jobz == 'V') then
                !$omp critical
                call dsyev(jobz, 'U', dim, matloc, lda, evalsloc, work, lwork, info)
                !$omp end critical
            else if(jobz == 'N') then
                call dsyev_2stage('N', 'U', dim, matloc, lda, evalsloc, work, lwork, info)
            end if
            lwork = int(work(1))
            if(lwork < max(1, 6*(dim-1) ) ) lwork = max(1, 6*(dim-1) )
            if(allocated(work)) deallocate(work)
            allocate(work(lwork))
            work = 0.0d0

            if(jobz == 'V') then
                !$omp critical
                call dsyev(jobz, 'U', dim, matloc, lda, evalsloc, work, lwork, info)
                !$omp end critical
            else if(jobz == 'N') then
                call dsyev_2stage('N', 'U', dim, matloc, lda, evalsloc, work, lwork, info)
            end if

            if(info .gt. 0) then
                write(*,*)'The algorithm failed to compute eigenvalues.'
                stop
            end if
            if(allocated(work)) deallocate(work)
        else
            evalsloc = matloc(1,1)
        end if

        evals = evalsloc
        mat   = matloc
        if(allocated(evalsloc)) deallocate(evalsloc)
        if(allocated(matloc)) deallocate(matloc)

        return

    end subroutine exactdiag_dp

    subroutine exactdiag_dc(eigvec, dim, mat, evals)
        
        implicit none
        logical,                       intent(in)    :: eigvec
        integer(kind=8),               intent(in)    :: dim
        double complex,                intent(inout) :: mat(dim, dim)
        double precision, allocatable, intent(out)   :: evals(:)

        double complex,   allocatable :: matloc(:,:)
        double precision, allocatable :: evalsloc(:)
        
        integer            :: lda = 0, info = 0, lwork = 0
        character          :: jobz
        integer, parameter :: lwmax = 100000
        double precision, allocatable :: rwork(:)
        double complex,   allocatable :: work(:)

        external :: zheev

        !$OMP THREADPRIVATE(matloc, evalsloc, lda, info, lwork, jobz, work)
        lda   = 0
        info  = 0
        lwork = 0

        if(eigvec) then
            jobz = 'V'
        else
            jobz = 'N'
        end if

        if(allocated(evals)) deallocate(evals)
        allocate(evals(dim))
        evals = 0.0d0

        if(allocated(evalsloc)) deallocate(evalsloc)
        if(allocated(matloc)) deallocate(matloc)
        allocate(evalsloc(dim))
        allocate(matloc(dim, dim))
        evalsloc = 0.0d0
        matloc = mat

        if(dim > 1) then
            lda = dim
            if(allocated(work)) deallocate(work)
            allocate(work(lwmax))
            if(allocated(rwork)) deallocate(rwork)
            allocate(rwork(max(1, 3*dim-2)))
            lwork = -1
            if(jobz == 'V') then
                !$omp critical
                call zheev(jobz, 'U', dim, matloc, lda, evalsloc, work, lwork, rwork, info)
                !$omp end critical
            end if
            lwork = int(work(1))
            if(lwork < max(1, 6*(dim-1) ) ) lwork = max(1, 6*(dim-1) )
            if(allocated(work)) deallocate(work)
            allocate(work(lwork))
            work  = 0.0d0
            rwork = 0.d0 

            if(jobz == 'V') then
                !$omp critical
                call zheev(jobz, 'U', dim, matloc, lda, evalsloc, work, lwork, rwork, info)
                !$omp end critical
            end if

            if(info .gt. 0) then
                write(*,*)'The algorithm failed to compute eigenvalues.'
                stop
            end if
            if(allocated(work))  deallocate(work)
            if(allocated(rwork)) deallocate(rwork)
        else
            evalsloc = matloc(1,1)
        end if
        evals = evalsloc
        mat   = matloc
        if(allocated(evalsloc)) deallocate(evalsloc)
        if(allocated(matloc)) deallocate(matloc)

        return

    end subroutine exactdiag_dc

    !---------------------------------------------------!
    !            Sparse diagonalization (ARPACK)       !
    !---------------------------------------------------!
    ! This subroutine uses ARPACK to perform Lanczos diagonalization.

    subroutine lanczos_dp(threads, dim, nnz, nev, ncv, nest, mode, rvec, ham, rc, evals, evecs)

        implicit none
        integer(kind=8), intent(in) :: dim, nnz
        integer, intent(in) :: threads, nev, ncv, nest
        integer(kind=8), intent(in) :: rc(nnz,2)
        double precision, intent(in) :: ham(nnz)
        character(len=*), intent(in) :: mode
        logical, intent(in)          :: rvec
        double precision, allocatable, intent(out) :: evals(:), evecs(:,:)

        integer(kind=8) :: maxn, n
        integer :: maxnev, maxncv, ldv, j
        integer :: ishfts, lworkl, maxitr, mode1, nconv, iparam(11), ipntr(11), ido, ierr, info
        integer(kind=8), allocatable :: jao(:), iao(:)
        double precision :: zero = 0.0D+00
        double precision :: sigma, tol
        double precision, allocatable :: ax(:), workl(:), workd(:), v(:,:), resid(:), d(:,:), ao(:)
        character :: bmat*1, which*2
        logical, allocatable :: selector(:)
        double precision, external :: dnrm2
        external :: daxpy, dsaupd, dseupd, dmout
        intrinsic :: abs

        integer :: cntr = 0

        !$omp threadprivate(maxn, n, maxnev, maxncv, ldv, j, ishfts, lworkl, maxitr, mode1, nconv, iparam, ipntr, ido, ierr, info, jao, iao, zero, sigma, tol, ax, workl, workd, v, resid, d, ao, bmat, which, selector)


        if(allocated(evals)) deallocate(evals)
        if(allocated(evecs)) deallocate(evecs)
        allocate(evals(nev))
        allocate(evecs(dim, nest))
        evals = 0.d0 
        evecs = 0.d0    

        if(dim == 1) then
            evecs(1,1) = 1
            evals(1) = ham(1)
            print*, ''
            write(*,"('Only eigenvalue is ',f8.4)") ham(1)
        else
            maxn   = 10 + dim
            maxnev = 10 + nev
            maxncv = 10 + ncv
            ldv    = maxn
        
            if(allocated(workd)) deallocate(workd)
            if(allocated(ax)) deallocate(ax)
            if(allocated(d)) deallocate(d)
            if(allocated(resid)) deallocate(resid)
            if(allocated(selector)) deallocate(selector)
            if(allocated(workl)) deallocate(workl)
            if(allocated(v)) deallocate(v)
            allocate(ax(maxn),d(maxncv,2),resid(maxn),selector(maxncv),v(ldv,maxncv),workl(maxncv*(maxncv+8)),workd(3*maxn))

            !SET ARPACK PARAMETERS: Make sure that 1) maxn >= n 2) maxnev >= nev 3) maxncv >= ncv
            !Set dimensions for this problem.
            n     = int(dim, 4)
            bmat  = 'I'
            which = mode

            if(n .gt. maxn) then
                print *, ' ERROR: N is greater than MAXN '
                go to 9870
            else if(nev .gt. maxnev) then
                print *, ' ERROR: NEV is greater than MAXNEV '
                go to 9870
            else if(ncv .gt. maxncv) then
                print *, ' ERROR: NCV is greater than MAXNCV '
                go to 9870
            end if

            lworkl = ncv *( ncv + 8 )
            tol    = zero
            info   = 0
            ido    = 0
            ishfts = 1
            maxitr = 300
            mode1  = 1
            iparam = 0
            cntr   = 0 

            iparam(1) = ishfts
            iparam(3) = maxitr
            iparam(7) = mode1
            iparam(4) = 1 
            !MAIN LOOP(Reverse communication loop)! Repeatedly call DSAUPD and take actions indicated by parameter !IDO until convergence is indicated or MAXITR is exceeded.
            ierr     = 0
            ax       = 0.d0
            v        = 0.d0
            workl    = 0.d0
            workd    = 0.d0
            sigma    = 0.d0
            d        = 0.d0
            resid    = 0.d0
            ipntr    = 0.d0
            selector = 0.d0

            !Create sparse matrix in CSR format
            if(allocated(jao)) deallocate(jao)
            if(allocated(iao)) deallocate(iao)
            if(allocated(ao)) deallocate(ao)
            allocate(jao(nnz))
            allocate(iao(dim+1))
            allocate(ao(nnz))
            jao = 0 
            iao = 0 
            ao  = 0 
            
            call coo_to_csr(dim, nnz, ham, rc(1:nnz,1), rc(1:nnz,2), ao, jao, iao) 
            
            !$omp critical
            do
            !Sparse eigensolver routine
            call dsaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info)
            if(ido /= -1 .and. ido /= 1) then
                exit
            end if
            !Vector-sparse matrix multiplication                 
            cntr = cntr + 1
            call amux_dp(threads, dim, workd(ipntr(1)), workd(ipntr(2)), ao, jao, iao)
            end do
            !$omp end critical
            print*,'Number of MV-multiplications',cntr
            if (info /= 0) then
                print *, 'DSSIMP - Fatal error! Error with DSAUPD, IERR = ', info, 'Check the documentation of DSAUPD.'
                go to 9870
            end if
            
            
            !Extract eigenvalues and -vectors
            !$omp critical
            call dseupd(rvec, 'A', selector, d, v, ldv, sigma, bmat, n, which, nev, tol, &
                        resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, ierr)
            !$omp end critical
            if(ierr /= 0) then
            print *, 'DSSIMP - Fatal error! Error with DSEUPD, IERR = ', ierr, 'Check the documentation of DSEUPD.'
            else
            nconv =  iparam(5)
            do j = 1, nconv
                call amux_dp(threads, dim, v(1,j), ax, ao, jao, iao)
                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1 )
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
            end do
            call dmout(6, nconv, 2, d, maxncv, -6, 'Ritz values are named relative residuals' ) !Display computed residuals. 6: Output to screen Write(6, #internalnumber)! nconv: number of rows in the matrix d! 2: Number of columns in matrix d! maxncv: Leading dimension of the matrix data! -6: print the matrix d with iabs(-6) decimal digits per number! Use formatting indexed by -6 to print A
        end if
            
        evals(1:nev) = d(1:nev,1)
        if(rvec) then
            do j= 1, min(iparam(5), nest) 
                evecs(1:dim, j) = v(1:dim,j)
            end do
        end if
        
        9870 continue

        end if


        100 format(1000(F40.30))

        return 

    end subroutine lanczos_dp

    subroutine lanczos_dp_otf(unit, threads, dim, sites, nbonds, nnnbonds, nev, ncv, nest, t1, v1, v2, mode, rvec, basis, bsites, hexsites, evals, evecs)

        implicit none

        integer(kind=8), intent(in) :: dim, basis(dim)
        integer, intent(in) :: unit, threads, sites, nbonds, nnnbonds, nev, ncv, nest
        integer, intent(in) :: bsites(2, nbonds), hexsites(2, nnnbonds)
        double precision, intent(in) :: t1, v1, v2
        character(len=*), intent(in) :: mode
        logical, intent(in)          :: rvec
        double precision, allocatable, intent(out) :: evals(:), evecs(:,:)

        integer(kind=8)  :: maxn, n
        integer  :: maxnev, maxncv, ldv, j, cntr = 0 
        integer  :: ishfts, lworkl, maxitr, mode1, nconv, iparam(11), ipntr(11), ido, ierr, info
        double precision  :: zero = 0.0D+00
        double precision  :: sigma, tol
        double precision, allocatable  :: ax(:), workl(:), workd(:), v(:,:), resid(:), d(:,:)
        character  :: bmat*1, which*2
        logical, allocatable  :: selector(:)
        double precision, external :: dnrm2
        external :: daxpy, dsaupd, dseupd, dmout
        intrinsic :: abs

        !$omp threadprivate(maxn, n, maxnev, maxncv, ldv, j, ishfts, lworkl, maxitr, mode1, nconv, iparam, ipntr, ido, ierr, info, zero, sigma, tol, ax, workl, workd, v, resid, d, bmat, which, selector)


        if(allocated(evals)) deallocate(evals)
        if(allocated(evecs)) deallocate(evecs)
        allocate(evals(nev))
        allocate(evecs(dim, nest))
        evals = 0.d0 
        evecs = 0.d0    

        if(dim == 1) then
            evecs(1,1) = 1
            evals(1) = 1
            print*, ''
            write(*,"('Only eigenvalue is ',f8.4)") evals(1)
        else
            maxn   = 10 + dim
            maxnev = 10 + nev
            maxncv = 10 + ncv
            ldv    = maxn
            
            if(allocated(v)) deallocate(v)
            if(allocated(d)) deallocate(d)
            if(allocated(ax)) deallocate(ax)
            if(allocated(resid)) deallocate(resid)
            if(allocated(workd)) deallocate(workd)
            if(allocated(workl)) deallocate(workl)
            if(allocated(selector)) deallocate(selector)
            allocate(ax(maxn),d(maxncv,2),resid(maxn),selector(maxncv),v(ldv,maxncv),workl(maxncv*(maxncv+8)),workd(3*maxn))

            !SET ARPACK PARAMETERS: Make sure that 1) maxn >= n 2) maxnev >= nev 3) maxncv >= ncv
            !Set dimensions for this problem.
            n     = int(dim, 4)
            bmat  = 'I'
            which = mode

            if(n .gt. maxn) then
                print *, ' ERROR: N is greater than MAXN '
                go to 9870
            else if(nev .gt. maxnev) then
                print *, ' ERROR: NEV is greater than MAXNEV '
                go to 9870
            else if(ncv .gt. maxncv) then
                print *, ' ERROR: NCV is greater than MAXNCV '
                go to 9870
            end if

            lworkl = ncv *( ncv + 8 )
            tol    = zero
            info   = 0
            ido    = 0
            ishfts = 1
            maxitr = 300! 5000! 1000 !300
            mode1  = 1
            iparam = 0

            iparam(1) = ishfts
            iparam(3) = maxitr
            iparam(7) = mode1
            iparam(4) = 1 
            !MAIN LOOP(Reverse communication loop)! Repeatedly call DSAUPD and take actions indicated by parameter !IDO until convergence is indicated or MAXITR is exceeded.
            ierr     = 0
            cntr     = 0
            ax       = 0.d0
            v        = 0.d0
            workl    = 0.d0
            workd    = 0.d0
            sigma    = 0.d0
            d        = 0.d0
            resid    = 0.d0
            ipntr    = 0.d0
            selector = 0.d0 
            !$omp critical
            do
            !Sparse eigensolver routine
            call dsaupd( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
            if(ido /= -1 .and. ido /= 1) then
                exit
            end if
            cntr = cntr + 1
            !Vector-sparse matrix multiplication                 
            call amux_dp_otf(threads, dim, sites, nbonds, nnnbonds, t1, v1, v2, basis, bsites, hexsites, workd(ipntr(1)), workd(ipntr(2)))
   
        end do
            !$omp end critical
            print*,'Number of MV-multiplications',cntr
            if (info /= 0) then
                print *, 'DSSIMP - Fatal error! Error with DSAUPD, IERR = ', info, 'Check the documentation of DSAUPD.'
                go to 9870
            end if
            
            !Extract eigenvalues and -vectors
            !$omp critical
            call dseupd( rvec, 'A', selector, d, v, ldv, sigma, bmat, n, which, nev, tol, &
                        resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, ierr )
            !$omp end critical
            if(ierr /= 0) then
            print *, 'DSSIMP - Fatal error! Error with DSEUPD, IERR = ', ierr, 'Check the documentation of DSEUPD.'
            else
            nconv =  iparam(5)
            do j = 1, nconv
                call amux_dp_otf(threads, dim, sites, nbonds, nnnbonds, t1, v1, v2, basis, bsites, hexsites, v(1,j), ax)
                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1 )
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
            end do
            call dmout(6, nconv, 2, d, maxncv, -6, 'Ritz values are named relative residuals' ) !Display computed residuals. 6: Output to screen Write(6, #internalnumber)! nconv: number of rows in the matrix d! 2: Number of columns in matrix d! maxncv: Leading dimension of the matrix data! -6: print the matrix d with iabs(-6) decimal digits per number! Use formatting indexed by -6 to print A
            end if

            evals(1:nev) = d(1:nev,1)

            if(rvec) then
                do j= 1, min(iparam(5), nest) 
                    evecs(1:dim, j) = v(1:dim,j)
                end do
            end if
            
            9870 continue

        end if


        100 format(1000(F40.30))
        return 

    end subroutine lanczos_dp_otf

    subroutine lanczos_dc( threads, dim, nev, ncv, nst, mode, rvec, nnz, ham, rc, evals, evecs )
            
            implicit none

            integer, intent(in) :: threads, nev, ncv, nst
            integer(kind=8), intent(in) :: dim, nnz
            integer(kind=8), intent(in) :: rc(nnz,2)
            double complex, intent(in) :: ham(nnz)   
            character*2, intent(in) :: mode
            logical, intent(in) :: rvec
            double precision, allocatable, intent(out) :: evals(:)
            double complex, allocatable, intent(out) :: evecs(:, :)
            integer :: printing = 1 
            integer(kind=8) :: maxn = 0, maxnev = 0, maxncv = 0, ldv = 0
            integer(kind=8) :: n = 0, nx = 0
            integer :: buffer = 0 
            integer :: j = 0, iparam(11), ipntr(14)
            integer :: ido = 0, ishfts = 0, lworkl = 0, info = 0, maxitr = 0, mode1 = 0, nconv = 0, ierr = 0    
            integer(kind=8), allocatable :: jao(:), iao(:), jao_ord(:)
            double precision :: tol = 0
            double precision, allocatable :: rwork(:), rd(:,:)
            double complex :: sigma = 0    
            double complex, allocatable :: ao(:), ao_ord(:)
            double complex, allocatable :: ax(:), d(:), v(:,:), workd(:), workev(:), resid(:), workl(:)
            character :: bmat*1, which*2
            logical, allocatable :: selector(:)

            intrinsic :: abs

            double precision :: dznrm2 , dlapy2
            external :: dznrm2 , zaxpy , dlapy2, znaupd, zneupd, dmout
            double precision, external :: dnrm2

            !$omp threadprivate(ao, jao, iao, bmat, which, j, iparam, sigma, ipntr, ido, ishfts, lworkl, info, maxitr, mode1, nconv, ierr, selector, ax, d, v, workd, workev, resid, workl, rwork, rd)
            
            
            if(dim == 1) then
                evecs(1,1) = 1
                evals(1) = ham(1)
                print*, ''
                write(*,"('Only eigenvalue is ',f8.4)") ham(1)
            else
            
            buffer = 50 
            maxn   = buffer + dim
            maxnev = buffer + nev
            maxncv = buffer + ncv
            ldv    = maxn

            if(allocated(evals))    deallocate(evals)
            if(allocated(evecs))    deallocate(evecs)
            if(allocated(workd))    deallocate(workd)
            if(allocated(workev))   deallocate(workev)
            if(allocated(rwork))    deallocate(rwork)
            if(allocated(rd))       deallocate(rd)
            if(allocated(ax))       deallocate(ax)
            if(allocated(d))        deallocate(d)
            if(allocated(resid))    deallocate(resid)
            if(allocated(selector)) deallocate(selector)
            if(allocated(workl))    deallocate(workl)
            if(allocated(v))        deallocate(v)
            allocate(evals(nev), evecs(dim, nst))
            allocate(ax(maxn), d(maxncv), resid(maxn), selector(maxncv), v(ldv,maxncv), &
                    workl(3*maxncv*maxncv+5*maxncv), workd(3*maxn), workev(3*maxncv), &
                    rwork(maxncv), rd(maxncv,3))
            ax     = 0.d0 
            d      = 0.d0 
            resid  = 0.d0 
            v      = 0.d0
            workl  = 0.d0
            workd  = 0.d0 
            workev = 0.d0
            rwork  = 0.d0 
            rd     = 0.d0 


            nx = dim
            n  = nx 

            if(n .gt. maxn) then
                print *, ' ERROR with _NDRV1: N is greater than MAXN '
                go to 9870
            else if(nev .gt. maxnev) then
                print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
                go to 9870
            else if(ncv .gt. maxncv) then
                print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
                go to 9870
            end if
            bmat   = 'I'
            which  = 'SR'
            lworkl = 3*ncv**2+5*ncv
            tol    = 0
            ido    = 0
            info   = 0
            nconv  = 0
            ishfts = 1
            maxitr = 300
            mode1  = 1

            iparam(1) = ishfts
            iparam(3) = maxitr
            iparam(7) = mode1

            

            !c
            !c Create sparse matrix in CSR format
            !c
            
            if(allocated(jao)) deallocate(jao)
            if(allocated(iao)) deallocate(iao)
            if(allocated(ao))  deallocate(ao)
            allocate(jao(nnz))
            allocate(iao(dim+1))
            allocate(ao(nnz))
            ao  = 0.d0
            iao = 0 
            jao = 0 
            
            call coo_to_csr(dim, nnz, ham, rc(1:nnz,1), rc(1:nnz,2), ao, jao, iao)
    
            !$omp critical
            do
                !Sparse eigensolver routine
                call znaupd( ido, bmat, n, which, nev, tol, resid, ncv, &
                            v, ldv, iparam, ipntr, workd, workl, lworkl, &
                            rwork, info )
                if(ido /= -1 .and. ido /= 1) then
                    exit
                end if
                !Vector-sparse matrix multiplication
                call amux_dc(threads, nx, workd(ipntr(1)), workd(ipntr(2)), ao, jao, iao)
            end do
            !$omp end critical
        
            if(info .lt. 0) then
                print *, ' '
                print *, ' Error with znaupd, info = ', info
                print *, ' Check the documentation of _naupd'
                print *, ' '
            else
                call zneupd(rvec, 'A', selector, d, v, ldv, sigma, &
                            workev, bmat, n, which, nev, tol, resid, ncv, &
                            v, ldv, iparam, ipntr, workd, workl, lworkl, &
                            rwork, ierr)
    

                if(ierr .ne. 0) then
                    print *, ' '
                    print *, ' Error with zneupd, info = ', ierr
                    print *, ' Check the documentation of _neupd. '
                    print *, ' '

                else

                    nconv = iparam(5)
                    do 20 j = 1, nconv
                        call amux_dc(threads, nx, v(1,j), ax, ao, jao, iao)
                        call zaxpy(n, -d(j), v(1,j), 1, ax, 1)
                        rd(j,1) = dble(d(j))
                        rd(j,2) = dimag(d(j))
                        rd(j,3) = dznrm2(n, ax, 1)
                        rd(j,3) = rd(j,3) / dlapy2(rd(j,1),rd(j,2))
                    20 continue
                    !c
                    !c            %-----------------------------%
                    !c            | Display computed residuals. |
                    !c            %-----------------------------%
                    !c
                    if(printing == 1) then 
                        call dmout(6, nconv, 3, rd, maxncv, -6, &
                                    'Ritz values(Real, Imag) and relative residuals')
                    end if 
                end if


                    ! c
                    ! c        %-------------------------------------------%
                    ! c        | Print additional convergence information. |
                    ! c        %-------------------------------------------%
                    ! c
                if(info .eq. 1) then
                    print *, ' '
                    print *, ' Maximum number of iterations reached.'
                    print *, ' '
                else if(info .eq. 3) then
                    print *, ' '
                    print *, ' No shifts could be applied during implicit', &
                            ' Arnoldi update, try increasing NCV.'
                    print *, ' '
                end if


            end if

            evals(1:nev) = dble(d(1:nev))
            evecs = 0 
            
            if(rvec) then
                do j = 1, nst
                    evecs(1:dim, j) = v(1:dim,j)
                end do
            end if

            
            end if

            9870 continue


            if(allocated(workd)) deallocate(workd)
            if(allocated(workev)) deallocate(workev)
            if(allocated(rwork)) deallocate(rwork)
            if(allocated(rd)) deallocate(rd)
            if(allocated(ax)) deallocate(ax)
            if(allocated(d)) deallocate(d)
            if(allocated(resid)) deallocate(resid)
            if(allocated(selector)) deallocate(selector)
            if(allocated(workl)) deallocate(workl)
            if(allocated(v)) deallocate(v)
            if(allocated(jao)) deallocate(jao)
            if(allocated(iao)) deallocate(iao)
            if(allocated(ao)) deallocate(ao)


    end subroutine lanczos_dc

    !-------------------------------------------------!
    !            Sparse diagonalization (MKL)        !
    !-------------------------------------------------!

    subroutine lanczos_mkl(dim, nev, ncv, nest, nnz, values, rc, evals, evecs)
        
        USE MKL_SPBLAS
        USE MKL_SOLVERS_EE    
        use iso_c_binding, only: c_double, c_int
        implicit none
        
        !!!!!!!!!!!!!!!!! Matrix declaration variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer(kind=8), intent(in) :: dim
        integer, intent(in) :: nev
        integer, intent(in) :: ncv
        integer, intent(in) :: nest
        integer(kind=8), intent(in) :: nnz
        integer(kind=8), intent(in) :: rc(nnz,2)
        double precision, intent(in) ::  values(nnz)
        
        double precision, allocatable, intent(out) :: evals(:)
        double precision, allocatable, intent(out) :: evecs(:,:)

        integer :: i, threads 
        !   Matrix descriptor
        TYPE(MATRIX_DESCR) :: descrA
        !   CSR matrix structure
        TYPE(SPARSE_MATRIX_T) ::  cooA, csrA

        !!!!!!!!!!!!!!!!! Declaration of MKL_SPARSE_S_EV variables !!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! E - eigenvalues, X - eigenvectors, res - residual !!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! XL - left singular vectors, XR - right singular vectors !!!!!!

        character*1 WHICH
        parameter (WHICH='S') !L - for largest eigenvalues to find
        integer, allocatable :: pm(:)
        integer :: n, K, info

        double precision, allocatable :: E(:)
        double precision, allocatable :: sigma(:)
        double precision, allocatable :: X(:,:)
        double precision, allocatable :: res(:)

    
        !$omp threadprivate(i, descrA, cooA, csrA, pm, n, k, info, E, sigma, X, res)

        
        if(allocated(E)) deallocate(E)
        if(allocated(X)) deallocate(X)
        if(allocated(pm)) deallocate(pm)
        if(allocated(res)) deallocate(res)
        if(allocated(evals)) deallocate(evals)
        if(allocated(evecs)) deallocate(evecs)
        if(allocated(sigma)) deallocate(sigma)
        
        allocate(E(nev))
        allocate(pm(128))
        allocate(res(dim))
        allocate(sigma(dim))
        allocate(evals(nev))
        allocate(X(dim, nev))
        allocate(evecs(dim, nest))    
        

        evals = 0.d0 
        evecs = 0.d0 
        n = int(dim, kind=4) 
        print *, 'Sparse matrix size', n
        !
        !        Task 0. Call MKL_SPARSE_C_CREATE_COO to create matrix handle 
        !      
        info = mkl_sparse_d_create_coo(cooA, SPARSE_INDEX_BASE_ONE, n, n, int(nnz, kind=4), int(rc(1:nnz,1), kind=4), int(rc(1:nnz,2), kind=4), values)

        !
        !        Task 1. Call MKL_SPARSE_C_CREATE_CSR to create matrix handle
        !      

        info = mkl_sparse_convert_csr(cooA, SPARSE_OPERATION_NON_TRANSPOSE, csrA)
        info = MKL_SPARSE_DESTROY(cooA)
        ! info = mkl_sparse_d_create_csr(csrA,SPARSE_INDEX_BASE_ZERO,N,N,rows,rows(2),cols,val)

        !         Create matrix descriptor
        descrA % TYPE = SPARSE_MATRIX_TYPE_GENERAL
        !
        !        Task 2. Call MKL_SPARSE_EE_INIT to define the default values for the input
        !        parameters.
        !
        info = mkl_sparse_ee_init(pm)
        pm(2) = 6 !Setting tolerance 10**(-pm(2) + 1)  
        pm(3) = 1
        !0 - Decided at runtime
        !1 - The Krylov-Schur method
        !2 - Subspace Iteration technique based on FEAST algorithm(for eigenvalues in specific interval)
        pm(4) = ncv 
        ! This parameter is referenced only for Krylov-Schur Method. It indicates the number of Lanczos/Arnoldi vectors(NCV) generated at each iteration.
        ! This parameter must be less than or equal to size of matrix and greater than number of eigenvalues(k0) to be computed. If unspecified, NCV is set to be at least 1.5 times larger than k0. 

        pm(8) = 1 ! Use absolute stopping criteria
        !
        !         Task 3. Solve the standard eigenvalue problem Ax=ex.
        ! 
        !

        !!$ threads = omp_get_num_threads()
        ! print*,'threads',threads
        !!$ threads = omp_get_thread_num()
        ! print*,'id',threads
        !!$omp critical 

        ! print* ,x, 'x'
        ! print* ,e, 'e'
        ! print* ,res, 'res'
        ! print* ,nest, 'nest'
        ! print* ,nev, 'nev'
        ! print* ,dim, 'dim'
        ! print* ,pm, 'pm'
        ! print* ,k, 'k'
        
        info = mkl_sparse_d_ev(WHICH,pm,csrA,descrA,nev,k,E,X,res)
        !!$omp end critical 
        print  *,' OUTPUT INFO ',info
        print  *,' ITERATION EXIT REASON(0 = converged) ', pm(10)
        
        if(info.ne.0) stop 1
        if(pm(10).ne.0) stop 2
        print *, 'Number of eigenvalues found ', k
        ! print *, ' Computed    |    Expected  '
        ! print *, ' Eigenvalues |    Eigenvalues '
        do i = 1, K
            ! write(*,'(f20.16)') E(i)
            ! print *, E(i)
            evals(i) = E(i)
        end do
        !   Release internal representation of CSR matrix
        info = MKL_SPARSE_DESTROY(csrA)

        ! print*, 'Number of eigenstates to be saved ', nest 
        do i = 1, nest 
            ! evecs(1:n,i) = X(i,1:n) 
            evecs(1:n,i) = X(1:n,i) 
        end do 
        
        return 

    end subroutine lanczos_mkl

    !------------------------------!
    !            FEAST             !
    !------------------------------!

    subroutine feast_dp(n, nnz, nev0, nest, rc, ham, Emin, Emax, nev, eigvals, eigvecs)
    
        !       1. The code calls  FEASTINIT  to define the default values for the input
        !          FEAST parameters.
        
        !       2. The  code  solves  the  standard eigenvalue  problem  Ax=ex   using 
        !          ZFEAST_HCSREV.
        
        !       3. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i) 
        !           are the expected eigenvalues  and E(i) are eigenvalues computed 
        !           with the help of ZFEAST_HCSREV().
        
        !       4. The code computes the maximum absolute value of the matrix  Y=(X')*X-I  
        !          where X is the matrix of eigenvectors computed with the help of 
        !          ZFEAST_HCSREV. ZGEMM(BLAS Level 3 Routine) is called  to compute(X')*X.
        ! *******************************************************************************
        USE MKL_SPBLAS   
        USE ISO_C_BINDING 
        implicit none
        !!!!!!!!!!!!!!!!! External function declarations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        external zgemm, zfeast_hcsrev, zfeast_hcsrgv, feastinit
        !!!!!!!!!!!!!!!!! Matrix declaration variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer, intent(in) :: nev0
        integer(kind=8),  intent(in) :: n, nnz
        integer(kind=8),  intent(in) :: rc(nnz,2)
        double precision, intent(in) :: ham(nnz)
        double precision, intent(in) :: Emin,Emax
        integer,          intent(inout) :: nest
        integer,          intent(out) :: nev 
        double precision, allocatable, intent(out) :: eigvals(:)
        double precision, allocatable, intent(out) :: eigvecs(:,:) 
        
        TYPE(MATRIX_DESCR) :: descrA     ! Sparse matrix descriptor
        !   CSR matrix representation
        TYPE(SPARSE_MATRIX_T) :: cooA, csrA    ! Structure with sparse matrix

        integer(kind=8)  :: rows(n+1), cols(nnz), rows_b(n+1), cols_b(n) 
        integer(C_INT)   :: indexing
        double precision :: val(nnz), val_b(n)
        double precision :: beta
        character*1 :: UPLO
        parameter  (UPLO='F')
        !!!!!!!!!!!!!!!!! Feast declaration variable !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! E - eigenvalues, X - eigenvectors, res - residual !!!!!!!!!!!!
        integer :: fpm(128), loop, L, M0, M, info  
        parameter  (L=10)
        double precision :: epsout, E(nev0)
        !On output, the first m columns of x contain the orthonormal eigenvectors corresponding 
        !to the computed eigenvalues e, with the i-th column of x holding the eigenvector associated with e(i).   
        double precision :: X(n,nev0), res(nev0)
        !!!!!!!!!!!!!!!!! Declaration of local variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! Eig - exact eigenvalues, R=|E-Eig| !!!!!!!!!!!!!!!
        double precision :: Eig(nev0), R(nev0)
        integer :: i,j, mkl, dim, UTnnz = 0, ldx, ldy
        double precision :: one, zero, smax, eigabs 

        integer(kind=8),  allocatable :: ib(:), jb(:), ic(:), jc(:)
        double precision, allocatable :: b(:), c(:)
        integer,          allocatable :: UTcols(:), UTrows(:)
        double precision, allocatable :: UTval(:)

        print*,'Sparse matrix size',n
        print*,'Use Feast algorithm'

        M0=nev0
        M=10
        print *,'Search interval ', Emin,' ', Emax
        ldx = n
        ldy = n

        if(allocated(b)) deallocate(b)
        if(allocated(jb)) deallocate(jb)
        if(allocated(ib)) deallocate(ib)
        allocate(b(n))
        allocate(jb(n))
        allocate(ib(n))
        val    = 0.d0
        val_b  = 0.d0
        cols   = 0 
        rows   = 0 
        cols_b = 0 
        rows_b = 0 

        do I = 1, N
            b(i)  = 1.d0
            ib(i) = i 
            jb(i) = i
        end do
        call coo_to_csr(n, nnz, ham, rc(1:nnz,1), rc(1:nnz,2), val, cols, rows)
        call coo_to_csr(n, n, b, ib, jb, val_b, cols_b, rows_b)

        deallocate(b, ib, jb)
        if(allocated(c)) deallocate(c)
        if(allocated(jc)) deallocate(jc)
        if(allocated(ic)) deallocate(ic)
        allocate(ic(n+1))
        ic = 0 
        allocate(c(nnz))
        allocate(jc(nnz))
        jc = 0 
        c  = 0.d0

        call mkl_dcsrmultcsr('N', 0, 3, n, n, n, val, cols, rows, val_b, cols_b, rows_b, c, jc, ic, nnz, info)

        deallocate(c, ic, jc)
    
        !
        !        Task 1. Call  FEASTINIT  to define the default values for the input
        !        FEAST parameters.
        !
        call feastinit(fpm)
        fpm(1)  = 1 !Specifies whether Extended Eigensolver routines print runtime status. 
        fpm(26) = 1 !Specifies whether Extended Eigensolver routines check input matrices(applies to CSR format only). 
        fpm(27) = 1 !Specifies whether Extended Eigensolver routines check input matrices(applies to CSR format only). 

        print *, ' Testing feast_dp_scsrev '
        !
        !         Task 2. Solve the standard eigenvalue problem Ax=ex.
        !
        
        call zfeast_hcsrev(UPLO,N,val,rows,cols,fpm,epsout,loop,&
        Emin,Emax,M0,E,X,M,res,info)
        
        print*,' FEAST OUTPUT INFO ', info
        if(info.ne.0) stop 1
        nev  = m 
        nest = min(nest, nev)
        if(allocated(eigvals)) deallocate(eigvals)
        if(allocated(eigvecs)) deallocate(eigvecs)
        allocate(eigvals(nev))
        allocate(eigvecs(n,nest))
        eigvals = 0.d0 
        eigvecs = 0.d0
        do i = 1, nev
            eigvals(i) = E(i)        
        end do
        do i = 1, nest
            eigvecs(1:n, i) = X(1:n, i)
        end do 

        print*, 'Number of eigenvalues found ', M
        print*, ' Computed    |    Expected  '
        print*, ' Eigenvalues |    Eigenvalues '
        eigabs = zero
        do i = 1, M
        R(i)   = dabs(E(i)-Eig(i))
        eigabs = max(eigabs, R(i))
        print*, E(i), Eig(i)
        enddo
        print*, ' Max value of '
        print*, ' | computed eigenvalue -expected eigenvalues | ', eigabs

        return 

    end subroutine feast_dp

    subroutine feast_dc(n, nnz, nev0, nest, rc, ham, Emin, Emax, nev, eigvals, eigvecs)
    
        !       1. The code calls  FEASTINIT  to define the default values for the input
        !          FEAST parameters.
        
        !       2. The  code  solves  the  standard eigenvalue  problem  Ax=ex   using 
        !          ZFEAST_HCSREV.
        
        !       3. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i) 
        !           are the expected eigenvalues  and E(i) are eigenvalues computed 
        !           with the help of ZFEAST_HCSREV().
        
        !       4. The code computes the maximum absolute value of the matrix  Y=(X')*X-I  
        !          where X is the matrix of eigenvectors computed with the help of 
        !          ZFEAST_HCSREV. ZGEMM(BLAS Level 3 Routine) is called  to compute(X')*X.
        ! *******************************************************************************
        use MKL_SPBLAS   
        use ISO_C_BINDING 
        implicit none
        !!!!!!!!!!!!!!!!! External function declarations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        external zgemm, zfeast_hcsrev, zfeast_hcsrgv, feastinit
        !!!!!!!!!!!!!!!!! Matrix declaration variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer(kind=8), intent(in) :: n
        integer(kind=8), intent(in) :: nnz
        integer, intent(in) :: nev0
        integer(kind=8), intent(in) :: rc(nnz,2)
        double complex, intent(in) :: ham(nnz)
        double precision, intent(in) :: Emin, Emax
        integer, intent(inout) :: nest
        integer, intent(out) :: nev 
        double precision, allocatable, intent(out) :: eigvals(:)
        double complex, allocatable, intent(out) :: eigvecs(:,:) 
        
        type(MATRIX_DESCR) :: descrA     ! Sparse matrix descriptor
        !   CSR matrix representation
        type(SPARSE_MATRIX_T) :: cooA, csrA    ! Structure with sparse matrix
        
        integer(kind=8) :: rows(n+1), cols(nnz), rows_b(n+1), cols_b(n)
        integer(C_INT)  :: indexing
        double complex  :: val(nnz), val_b(n), beta
        character*1 :: UPLO
        parameter  (UPLO='U')
        !!!!!!!!!!!!!!!!! Feast declaration variable !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! E - eigenvalues, X - eigenvectors, res - residual !!!!!!!!!!!!
        integer :: fpm(128), loop, L, M0, M, info  
        parameter  (L=10)
        double precision :: epsout, E(nev0)
        !On output, the first m columns of x contain the orthonormal eigenvectors corresponding 
        !to the computed eigenvalues e, with the i-th column of x holding the eigenvector associated with e(i).   
        double complex :: X(n,nev0), res(nev0) 
        !!!!!!!!!!!!!!!!! Declaration of local variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! Eig - exact eigenvalues, R=|E-Eig| !!!!!!!!!!!!!!!
        double precision :: Eig(nev0), R(nev0)

        integer          :: i, j, mkl, dim, UTnnz = 0, ldx, ldy 
        double complex   :: one, zero
        double precision :: smax, eigabs 

        integer(kind=8), allocatable :: ib(:), jb(:), ic(:), jc(:)
        integer,         allocatable :: UTcols(:), UTrows(:)
        double complex,  allocatable :: b(:), c(:)
        double complex,  allocatable :: UTval(:)

        print*,'Sparse matrix size',n
        print*,'Use Feast algorithm'
        M0 = nev0
        M  = int(0.5*M0)
        print *,'Search interval ', Emin,' ', Emax
        ldx = n
        ldy = n

        if(allocated(b)) deallocate(b)
        if(allocated(jb)) deallocate(jb)
        if(allocated(ib)) deallocate(ib)
        allocate(b(n))
        allocate(jb(n))
        allocate(ib(n))
        val =(0.d0,0.d0)
        val_b =(0.d0,0.d0)
        cols = 0 
        rows = 0 
        cols_b = 0 
        rows_b = 0 

        do I = 1, N
            b(i)  = 1.d0
            ib(i) = i 
            jb(i) = i
        end do

        call coo_to_csr(n, nnz, ham, rc(1:nnz,1), rc(1:nnz,2), val, cols, rows)
        call coo_to_csr(n, n, b, ib, jb, val_b, cols_b, rows_b)

        deallocate(b, ib, jb)
        if(allocated(c)) deallocate(c)
        if(allocated(jc)) deallocate(jc)
        if(allocated(ic)) deallocate(ic)
        allocate(ic(n+1))
        ic = 0 
        allocate(c(nnz))
        allocate(jc(nnz))
        jc = 0 
        c  = 0.d0 
        
        call mkl_zcsrmultcsr('N', 0, 3, n, n, n, val, cols, rows, val_b, cols_b, rows_b, c, jc, ic, nnz, info)
        deallocate(c, ic, jc)

        call HermitianUpperTriangularCSR(n, nnz, val, rows, cols, UTval, UTrows, UTcols, UTnnz)
   
        call feastinit(fpm)
        fpm(1)=1  !Specifies whether Extended Eigensolver routines print runtime status. 
        fpm(26)=1 !Specifies whether Extended Eigensolver routines check input matrices(applies to CSR format only). 
        fpm(27)=1 !Specifies whether Extended Eigensolver routines check input matrices(applies to CSR format only). 

        print *, ' Testing zfeast_hcsrev '

        call zfeast_hcsrev(UPLO,N,UTval,UTrows,UTcols,fpm,epsout,loop,&
        Emin,Emax,M0,E,X,M,res,info)

        print  *,' FEAST OUTPUT INFO ',info
        if(info.ne.0) stop 1
        nev = m 
        nest = min(nest, nev)
        if(allocated(eigvals)) deallocate(eigvals)
        if(allocated(eigvecs)) deallocate(eigvecs)
        allocate(eigvals(nev))
        allocate(eigvecs(n, nest))
        eigvals = 0.d0 
        eigvecs = 0.d0 
        do i = 1, nev
            eigvals(i) = E(i)        
        end do
        do i = 1, nest
            eigvecs(1:n, i) = X(1:n, i)
        end do 
  
        print *, 'Number of eigenvalues found ', M
        print *, ' Computed    |    Expected  '
        print *, ' Eigenvalues |    Eigenvalues '
        eigabs=zero
        do i=1,M
        R(i)=dabs(E(i)-Eig(i))
        eigabs=max(eigabs, R(i))
        print *, E(i), Eig(i)
        enddo
        print *, ' Max value of '
        print *, ' | computed eigenvalue -expected eigenvalues | ', eigabs
    
        
        return 

    end subroutine feast_dc

    subroutine HermitianUpperTriangularCSR(n, nnz, A, IA, JA, UTriangular, IA_UTriangular, JA_UTriangular, nnz_UpperTriangular)
        implicit none
        integer(kind=8), intent(in) :: n, nnz
        integer(kind=8), intent(in) :: IA(n+1), JA(nnz)
        double complex,  intent(in) :: A(nnz)

        integer,                     intent(out) :: nnz_UpperTriangular
        integer,        allocatable, intent(out) :: IA_UTriangular(:), JA_UTriangular(:)
        double complex, allocatable, intent(out) :: UTriangular(:)
    
        ! Calculate the number of non-zero elements in the upper triangular matrix
        integer :: i, j 

        nnz_UpperTriangular = 0
        do i = 1, n
        do j = IA(i), IA(i + 1) - 1
            if(JA(j) >= i) nnz_UpperTriangular = nnz_UpperTriangular + 1
        end do
        end do
        
        ! Allocate memory for the upper triangular matrix
        allocate(UTriangular(nnz_UpperTriangular))
        allocate(IA_UTriangular(n + 1))
        allocate(JA_UTriangular(nnz_UpperTriangular))
    
        ! Set the first element of IA_UTriangular to 1
        IA_UTriangular = 0 
        IA_UTriangular(1) = 1

        ! Calculate the CSR format for the upper triangular matrix
        nnz_UpperTriangular = 1
        IA_UTriangular(1) = nnz_UpperTriangular
        do i = 1, n
            do j = IA(i), IA(i + 1) - 1
            if(JA(j) >= i) then
                UTriangular(nnz_UpperTriangular) = A(j)
                JA_UTriangular(nnz_UpperTriangular) = JA(j)
                nnz_UpperTriangular = nnz_UpperTriangular + 1
            end if
            end do
            IA_UTriangular(i + 1) = nnz_UpperTriangular
        end do
        nnz_UpperTriangular = nnz_UpperTriangular - 1

    end subroutine HermitianUpperTriangularCSR




    subroutine gsdeg_dp(dir, rvec, unit, dim, nev, nest, thresh, energies, states, ndeg, gs)
        implicit none
        
        logical, intent(in)          :: rvec
        integer(kind=8), intent(in)  :: dim
        integer, intent(in)          :: unit, nev, nest
        double precision, intent(in) :: thresh 
        double precision, intent(in) :: energies(nev)
        double precision, intent(in) :: states(dim, nest)
        character(len=*), intent(in) :: dir 
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
            if(abs(energies(j) - energies(1)) <= thresh) then 
                ndeg = ndeg + 1
                if(rvec) gs = gs + states(1:dim,j)
            else 
                exit 
            end if 
        end do 
        
        if(nest < nev) then 
            if(abs(energies(nest + 1) - energies(1)) <= thresh) print*,'Potentially not full degeneracy has been captured. Increase "nest".'
        end if
        if(nest .le. ndeg) print*,'Potentially not full degeneracy has been captured. Increase "nest".'
        ndeg = min(ndeg, nest)
        if(rvec) gs = gs/sqrt(dble(ndeg))
        print*,'Ground state degeneracy = ',ndeg 
    

        return 
    
    end subroutine gsdeg_dp

    subroutine gsdeg_dc(dir, rvec, unit, dim, nev, nest, thresh, energies, states, ndeg, gs)
        implicit none

        logical,          intent(in) :: rvec
        integer(kind=8),  intent(in) :: dim
        integer,          intent(in) :: unit, nev, nest
        double precision, intent(in) :: thresh 
        double precision, intent(in) :: energies(nev)
        double complex,   intent(in) :: states(dim, nest)
        character(len=*), intent(in) :: dir 
        integer,         intent(out) :: ndeg 
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
            if(abs(energies(j) - energies(1)) <= thresh) then 
                ndeg = ndeg + 1
                if(rvec) gs = gs + states(1:dim,j)
            else 
                exit 
            end if 
        end do 
        if(nest < nev) then
            if(abs(energies(nest + 1) - energies(1)) <= thresh) print*,'Potentially not full degeneracy has been captured. Increase "nest".'
        end if
        if(nest .le. ndeg ) print*,'Potentially not full degeneracy has been captured. Increase "nest".'
        ndeg = min(ndeg, nest)
        if(rvec) gs = gs/sqrt(dble(ndeg))
        print*,'Ground state degeneracy = ',ndeg 

        return 
    
    end subroutine gsdeg_dc

    subroutine qgsdeg_dp(dir, unit, dim, nev, nest, energies, states, ndeg, gs)
        implicit none
        
        integer(kind=8), intent(in) :: dim
        integer, intent(in) :: unit, nev, nest
        double precision, intent(in) :: energies(nev)
        double precision, intent(in) :: states(dim, nest)
        character(len=*), intent(in) :: dir 
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
        end do
        ndeg = maxloc(real(dE), dim=1) 
                
        if(ndeg > 1) then 
            do j = 2, ndeg 
                gs = gs + states(1:dim,j)
            end do 
        end if 
        gs = gs/sqrt(dble(ndeg))
        print*,'Ground state quasi-degeneracy = ',ndeg 

        return 
    
    end subroutine qgsdeg_dp

    subroutine qgsdeg_dc(dir, unit, dim, nev, nest, energies, states, ndeg, gs)
        implicit none
        
        integer(kind=8), intent(in) :: dim
        integer, intent(in) :: unit, nev, nest
        double precision, intent(in) :: energies(nev)
        double complex, intent(in) :: states(dim, nest)
        character(len=*), intent(in) :: dir 
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
        end do
        ndeg = maxloc(real(dE), dim=1) 

        if(ndeg > 1) then 
            do j = 2, ndeg 
                gs = gs + states(1:dim,j)
            end do 
        end if 
        gs = gs/sqrt(dble(ndeg))
        print*,'Ground state quasi-degeneracy = ',ndeg 

        return 
    
    end subroutine qgsdeg_dc





    subroutine unify_i8(t, v1, v2, w, mass, pattern, sites, occ, nOff, nDi, hamOff, hamDi, ham, rc, nnz)
        !---------------------------------------------!
        !            Unify sparse Hamiltonian         !
        !---------------------------------------------!
        ! This subroutine unifies the sparse Hamiltonian parts into a single matrix
        ! in the compressed sparse row format.
        ! Use this routine if no basis symmetrization is used. 
        implicit none
        
        integer,          intent(in) :: sites
        integer(kind=8),  intent(in) :: nOff, nDi
        integer(kind=8),  intent(in) :: hamOff(nOff,3), hamDi(nDi,2) ,occ(sites,nDi)
        double precision, intent(in) :: t, v1, v2, w, mass 
        character*2,      intent(in) :: pattern

        integer(kind=8),               intent(out) ::  nnz
        integer(kind=8),  allocatable, intent(out) :: rc(:,:)
        double precision, allocatable, intent(out) :: ham(:)

        integer          :: j = 0, ab = 0, slstaggering(sites) 
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
            slstaggering(j) =(-1)**(j-ab)
        end do 
        
        call random_number(rand)

        rand = 2 *(rand - 0.5)

        ham(1:nOff)  = - t * hamOff(1:nOff,1)
        rc(1:nOff,1) = hamOff(1:nOff,2)
        rc(1:nOff,2) = hamOff(1:nOff,3)
        do j = 1, nDi !Fill diagonal Hamiltonian
            ham(nOff+j)  = v1 * 0.25 * hamDi(j,1) + v2 * 0.25 * hamDi(j,2) + mass * dot_product(slstaggering, occ(1:sites,j)) + w * sum(rand * occ(1:sites,j))
            rc(nOff+j,1) = j
            rc(nOff+j,2) = j
        end do
        
        return
        
    end subroutine unify_i8

    subroutine unify_dp(v1, v2, w, mass, pattern, sites, occ, nOff, nDi, nDi2, hamOff, rcoff, rcdi, hamDi2, hamDi, ham, rc, nnz)
        !---------------------------------------------!
        !            Unify sparse Hamiltonian         !
        !---------------------------------------------!
        ! This subroutine unifies the sparse Hamiltonian parts into a single matrix
        ! in the compressed sparse row format.
        ! Use this routine for purely real basis states.  
        implicit none

        integer,          intent(in) :: sites
        integer(kind=8),  intent(in) :: nOff, nDi, nDi2
        integer(kind=8),  intent(in) :: rcoff(nOff,2), rcdi(nDi2), occ(sites,nDi), hamDi(nDi,2)
        double precision, intent(in) :: v1, v2, w, mass 
        double precision, intent(in) :: hamOff(nOff), hamDi2(nDi2)
        character*2,      intent(in) :: pattern

        integer(kind=8),               intent(out) :: nnz
        integer(kind=8),  allocatable, intent(out) :: rc(:,:)
        double precision, allocatable, intent(out) :: ham(:)

        integer(kind=8)  :: j = 0 
        integer          :: slstaggering(sites), ab = 0
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
            slstaggering(j) =(-1)**(j-ab)
        end do   
        
        call random_number(rand)
        rand = 2 *(rand - 0.5)
        ham(1:nOff)  = hamOff(1:nOff)
        rc(1:nOff,1) = rcoff(1:nOff,1)
        rc(1:nOff,2) = rcoff(1:nOff,2)

        if(v1 == 0.d0 .and. v2 == 0.d0 .and. mass == 0.d0) go to 110 !Don't fill diagonal with zeroes
        do j = 1, nDi !Fill diagonal Hamiltonian        
            ham(nOff+j)  = v1 * 0.25 * hamDi(j,1) + v2 * 0.25 * hamDi(j,2) + mass * dot_product(slstaggering, occ(1:sites,j)) + w * sum(rand * occ(1:sites,j))        
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
        !---------------------------------------------!
        !            Unify sparse Hamiltonian         !
        !---------------------------------------------!
        ! This subroutine unifies the sparse Hamiltonian parts into a single matrix
        ! in the compressed sparse row format.
        ! Use this routine for complex basis states.  
        implicit none

        integer,          intent(in) :: sites 
        integer(kind=8),  intent(in) :: nOff, nDi, occ(sites,nDi), rcoff(nOff,2), hamDi(nDi,2)
        double precision, intent(in) :: v1, v2, w, mass
        double complex,   intent(in) :: hamOff(nOff)
        character*2,      intent(in) :: pattern

        integer(kind=8),              intent(out) :: nnz
        integer(kind=8), allocatable, intent(out) :: rc(:,:)
        double complex,  allocatable, intent(out) :: ham(:)

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
        rand         = 2 *(rand - 0.5)
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
        !---------------------------------------------!
        !            Unify sparse Hamiltonian         !
        !---------------------------------------------!
        ! This subroutine unifies the sparse Hamiltonian parts into a single matrix
        ! in the compressed sparse row format.
        ! Use this routine for symmetrization in 2D irreps with complex basis states.  
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
        rand         = 2 *(rand - 0.5)
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

    subroutine unifyDense_i8(dim, t, v1, v2, w, sites, occ, nOff, nDi, hamOff, hamDi, ham)
        !---------------------------------------------!
        !            Unify dense Hamiltonian          !
        !---------------------------------------------!
        ! This subroutine unifies the dense Hamiltonian parts into a single dense matrix.
        ! Use this routine if no basis symmetrization is used. 
        implicit none

        integer, intent(in)          :: sites
        integer(kind=8), intent(in)  :: dim,  nOff, nDi, occ(sites,*), hamOff(nOff,3), hamDi(nDi,2)       
        double precision, intent(in) :: t, v1, v2, w

        double precision, allocatable, intent(out) :: ham(:,:)

        integer(kind=8)  :: j = 0
        double precision :: rand(sites) 
        logical          :: symm 

        if(allocated(ham)) deallocate(ham)
        allocate(ham(dim, dim))
        ham = 0
        
        call random_number(rand)
        rand = 2 *(rand - 0.5)
        do j = 1, nOff
            ham(hamOff(j,2), hamOff(j,3))  = -t * hamOff(j,1)
        end do 
        do j = 1, nDi !Fill diagonal Hamiltonian
            ham(j,j) = v1 * hamDi(j,1) + v2 * hamDi(j,2) + w * sum(rand * occ(1:sites,j))
        end do
        call check_matrix(symm, dim, ham)
        
        return

    end subroutine unifyDense_i8

    subroutine unifyDense_dp(dim, v1, v2, w, sites, occ, nOff, nDi, nDi2, hamOff, rcoff, rcdi, hamDi2, hamDi, ham)
        !---------------------------------------------!
        !            Unify dense Hamiltonian          !
        !---------------------------------------------!
        ! This subroutine unifies the dense Hamiltonian parts into a single dense matrix.
        ! Use this routine for purely real basis states.  
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
        rand = 2 *(rand - 0.5)
        do j = 1, nOff
            ham(rcoff(j,1), rcoff(j,2)) = ham(rcoff(j,1), rcoff(j,2)) + hamOff(j)
        end do 
        do j = 1, nDi !Fill diagonal Hamiltonian
            ham(j,j) = v1 * hamDi(j,1) + v2 * hamDi(j,2) + w * sum(rand * occ(1:sites,j))
        end do
        call check_matrix(symm, dim, ham)
        
        do j = 1, nDi2 
            ham(rcdi(j), rcdi(j)) = ham(rcdi(j), rcdi(j)) + hamDi2(j)
        end do 

        return

    end subroutine unifyDense_dp

    subroutine unifyDense_dc(dim, v1, v2, w, sites, occ, nOff, nDi, hamOff, rcoff, hamDi, ham)
        !---------------------------------------------!
        !            Unify dense Hamiltonian          !
        !---------------------------------------------!
        ! This subroutine unifies the dense Hamiltonian parts into a single dense matrix.
        ! Use this routine for complex basis states.  
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
        rand = 2 *(rand - 0.5)
        do j = 1, nOff
            ham(rcoff(j,1), rcoff(j,2)) = ham(rcoff(j,1), rcoff(j,2)) + hamOff(j)
        end do 
        do j = 1, nDi !Fill diagonal Hamiltonian
            ham(j,j) = v1 * hamDi(j,1) + v2 * hamDi(j,2) + w * sum(rand * occ(1:sites,j))
        end do

        return

    end subroutine unifyDense_dc


end module diagonalization