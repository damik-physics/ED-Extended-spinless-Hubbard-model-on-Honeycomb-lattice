module diagonalization

    use test_utils
    implicit none 

    contains 

    !----------------------------------------------!
    !               Diagonalization                !
    !----------------------------------------------!

    subroutine diagonalize(nconf, nev, ncv, full, v1, v2, threads, parameters, type, basis, bsites, hexsites, occ, noff, ndi, ndi_off, hamoff, hamdi, hamoff_d, hamdi_d, hamoff_c, hamdi_c, hamdi_off_c, ham, ham_c, ham_d, ham_dc, norm, rcoff, rcdi, rc, prts, dplcts, nnz, ndeg, unit, nest, mode, energies, eigstate, eigstate_c, gs, gs_c)

    ! Old subroutine interface     
    ! subroutine diagonalize(dir, nconf, nev, ncv, full, v1, v2, threads, parameters, type, basis, bsites, hexsites, occ, noff, ndi, ndi_off, hamoff, hamdi, hamoff_d, hamdi_d, hamoff_c, hamdi_c, hamdi_off_c, ham, ham_c, ham_d, ham_dc, norm, rcoff, rcdi, rc, prts, dplcts, nnz, ndeg, unit, nest, mode, energies, eigstate, eigstate_c, gs, gs_c)

        use input_variables
        use variables

        implicit none 
        
        integer, intent(in) :: nconf, threads
        integer, intent(inout) :: nev, ncv, full, nest, unit, noff, ndi, ndi_off
        ! integer(kind=8), intent(in) :: dim
        double precision, intent(in) :: v1, v2
        character(len=*), intent(in) :: irrep, parameters, type, mode
        integer, allocatable, intent(inout) :: occ(:,:), hamoff(:,:), rcoff(:,:), rcdi(:), hamdi(:,:), prts(:), dplcts(:), bsites(:,:), hexsites(:,:)
        integer(kind=8), allocatable, intent(inout) :: basis(:)
        double precision, allocatable, intent(inout) :: hamoff_d(:), hamdi_d(:), norm(:)
        double complex, allocatable, intent(inout) :: hamoff_c(:), hamdi_c(:,:), hamdi_off_c(:)
        integer, intent(out) :: nnz, ndeg  
        integer, allocatable, intent(out) :: rc(:,:) 
        double precision, allocatable, intent(out) :: energies(:), eigstate(:,:), gs(:), ham(:), ham_d(:,:)
        double complex, allocatable, intent(out) :: eigstate_c(:,:), gs_c(:), ham_c(:), ham_dc(:,:)

        integer :: i = 0 
        double precision :: Emin = 0.d0, Emax = 0.d0
        logical :: symmetric, append 
        
        if((nDis > 1) .and. (dis .ne. 0.d0) .and. (nconf > 1)) then 
            append = .true.
        else 
            append = .false.
        end if 
        
        100 format(1000(x,A,x,F6.4,x,A,x,F6.4))

        if(full == 0) then
            if(type == "R") then 
                if(ti == 0) then 
                    call unify(t, v1, v2, dis, mass, pattern, sites, occ, noff, ndi, hamoff, hamdi, ham, rc, nnz)
                    if(arpack == 1) then 
                        if(otf == 0) then
                            call lanczos_d(threads, dim, nnz, nev, ncv, nest, mode, rvec, ham, rc, energies, eigstate)          
                        else if(otf == 1) then 
                            call lanczos_d_otf(parameters, unit, threads, dim, sites, sites/2*3, sites*3, nev, ncv, nest, t, v1, v2, mode, rvec, basis, bsites, hexsites, energies, eigstate)
                        end if 
                        if(rvec) call check_spectrum_dp(dim, nnz, nev, nest, energies, eigstate)                  
                        call save(dir, "R", parameters, append, nconf, unit, 3, dim, states, nev, nest, nDis, rvec, energies, eigstate)
                    end if 
                    if(mkl == 1) then 
                        call lanczos_mkl(dim, nev, ncv, nest, nnz, ham, rc, energies, eigstate) 
                        if(rvec) call check_spectrum_dp(dim, nnz, nev, nest, energies, eigstate)                  
                        call save(dir, "R", parameters, append, nconf, unit, 2, dim, states, nev, nest, nDis, rvec, energies, eigstate)
                    end if 
                    if(degflag == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate, ndeg, gs)     
                    if(degflag == 2 .and. rvec) call qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, eigstate, ndeg, gs)
                else if(ti == 1) then 
                    call unify_real(v1, v2, dis, mass, pattern, sites, occ, noff, ndi, ndi_off, hamoff_d, rcoff, rcdi, hamdi_d, hamdi, ham, rc, nnz)
                    call test_symm_sparse(symmetric, dim, nnz, rc, ham, norm, prts, dplcts)
                    if(arpack == 1) then 
                        call lanczos_d(threads, dim, nnz, nev, ncv, nest, mode, rvec, ham, rc, energies, eigstate)
                        if(degflag == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate, ndeg, gs)
                        if(degflag == 2 .and. rvec) call qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, eigstate, ndeg, gs)
                        if(rvec) call check_spectrum_dp(dim, nnz, nev, nest, energies, eigstate, ndeg = ndeg, rc = rc, mat = ham) 
                        call save(dir, "R", parameters, append, nconf, unit, 3, dim, states, nev, nest, nDis, rvec, energies, eigstate)    
                        if(degflag == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate, ndeg, gs)
                        if(degflag == 2 .and. rvec) call qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, eigstate, ndeg, gs)
                    end if 
                    if(mkl == 1) then 
                        call lanczos_mkl(dim, nev, ncv, nest, nnz, ham, rc, energies, eigstate)                
                        if(degflag == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate, ndeg, gs)
                        if(degflag == 2 .and. rvec) call qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, eigstate, ndeg, gs)
                        if(rvec) call check_spectrum_dp(dim, nnz, nev, nest, energies, eigstate, ndeg = ndeg, rc = rc, mat = ham) 
                        call save(dir, "R", parameters, append, nconf, unit, 2, dim, states, nev, nest, nDis, rvec, energies, eigstate)
                        if(degflag == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate, ndeg, gs)     
                        if(degflag == 2 .and. rvec) call qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, eigstate, ndeg, gs)    
                    end if 
                end if
                if(feast == 1) then 
                    Emin = dble(energies(1)) - 0.01
                    Emax = dble(energies(min(nev, nevmax)))+0.01
                    
                    call dfeast(dim, nnz, nev0, nest, rc, ham, Emin, Emax, nev, energies, eigstate)

                    if(degflag == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate, ndeg, gs)
                    if(degflag == 2 .and. rvec) call qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, eigstate, ndeg, gs)
                    if(rvec) call check_spectrum_dp(dim, nnz, nev, nest, energies, eigstate, ndeg = ndeg, rc = rc, mat = ham)
                    call save(dir, "R", parameters, append, conf, unit, 1, dim, states, nev, nest, nDis, rvec, energies, eigstate)
                    if(degflag == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate, ndeg, gs)
                    if(degflag == 2 .and. rvec) call qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, eigstate, ndeg, gs)
                end if                          
            else if(ti == 1 .and. type == "C") then 
                if(irrep(1:1) == "E") then 
                    call unify_c_2D(sites, noff, ndi_off, dim, v1, v2, dis, mass, pattern, occ, hamoff_c, rcoff, hamdi_off_c, rcdi, hamdi_c, ham_c, rc, nnz)
                else
                    call unify_comp(v1, v2, dis, mass, pattern, sites, occ, noff, ndi, hamoff_c, rcoff, hamdi, ham_c, rc, nnz)   
                end if 
        
                call test_hermitian_sparse(symmetric, irrep, symm, dim, nnz, rc, ham_c, norm, prts, dplcts)
                call lanczos_c(threads, dim, nev, ncv, nest, mode, rvec, nnz, ham_c, rc, energies, eigstate_c)
                if(rvec) call check_spectrum_dc(dim, .false., nev, nest, energies, eigstate_c)
                call save(dir, "C", parameters, append, conf, unit, 3, dim, states, nev, nest, nDis, rvec, energies, st_c=eigstate_c)          
                if(degflag == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate_c, ndeg, gs_c)      
                if(degflag == 2 .and. rvec) call qgsdeg_c(dir, unit, parameters, dim, nev, nest, energies, eigstate_c, ndeg, gs_c)
                if(feast == 1) then 
                    Emin = dble(energies(1))-0.01      
                    Emax = dble(energies(min(nev, nevmax)))+0.01
                    print*,'Emin',Emin
                    print*,'Emax',Emax
                    print*,'Emax-Emin',Emax-Emin
                    print*,''
                    print*,'---------------- START COMPLEX FEAST DIAGONALIZATION: V1 =', v1, ' V2 =', v2, ' -------------------'
                    call cfeast(dim, nnz, nev0, nest, rc, ham_c, Emin, Emax, nev, energies, eigstate_c)
                    if(rvec) call check_spectrum_dc(dim, .True., nev, nest, energies, eigstate_c)
                    call save(dir, "C", parameters, append, conf, unit, 1, dim, states, nev, nest, nDis, rvec, energies, st_c=eigstate_c)           
                    if(degflag == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate_c, ndeg, gs_c)            
                    if(degflag == 2 .and. rvec) call qgsdeg_c(dir, unit, parameters, dim, nev, nest, energies, eigstate_c, ndeg, gs_c)                      
                end if 
            end if
        else if(full == 1) then
            print*,''
            write(*,100) '---------------- START EXACT DIAGONALIZATION: V1 =', v1, ' V2 =', v2, ' -------------------'
            print*,''
            if(ti == 0 .and. type == "R") then 
                call unify_dense(dim, t, v1, v2, dis, sites, occ, noff, ndi, hamoff, hamdi, ham_d)
                call test_symm(symmetric, dim, ham_d)
                call exactdiag(rvec, dim, ham_d, energies)    
            else if(ti == 1 .and. type == "R") then 
                call unify_dense_d(dim, v1, v2, dis, sites, occ, noff, ndi, ndi_off, hamoff_d, rcoff, rcdi, hamdi_d, hamdi, ham_d)
                call test_symm(symmetric, dim, ham_d)
                call exactdiag(rvec, dim, ham_d, energies)
            else if(ti == 1 .and. type == "C") then 
                call unify_dense_c(dim, v1, v2, dis, sites, occ, noff, ndi, hamoff_c, rcoff, hamdi, ham_dc)
                call exactdiag_c(rvec, dim, ham_dc, energies)
            end if 
            if(type == "R") then 
                if(rvec) then 
                    if(allocated(eigstate)) deallocate(eigstate)
                    allocate(eigstate(dim, nest))
                    eigstate = 0.d0 
                    do i = 1, nest
                        eigstate(1:dim,i) = ham_d(1:dim,i)
                    end do
                    if(degflag == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate, ndeg, gs)
                    if(degflag == 2 .and. rvec) call qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, eigstate, ndeg, gs)
                    if(rvec) call check_spectrum_dp(dim, nnz, nev, nest, energies, eigstate, ndeg = ndeg)
                end if
                call save(dir, "R", parameters, append, conf, unit, 0, dim, states, nev, nest, nDis, rvec, energies, eigstate)
            else if(type == "C") then 
                if(rvec) then 
                    if(allocated(eigstate_c)) deallocate(eigstate_c)
                    allocate(eigstate_c(dim, nest))
                    eigstate_c = 0.d0 
                    do i = 1, nest
                        eigstate_c(1:dim,i) = ham_dc(1:dim,i)
                    end do
                    call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate_c, ndeg, gs_c)
                    if(rvec) call check_spectrum_dc(dim, .True., nev, nest, energies, eigstate_c)
                end if
                call save(dir, "C", parameters, append, conf, unit, 0, dim, states, nev, nest, nDis, rvec, energies, st_c = eigstate_c)
            end if 
        end if

    end subroutine diagonalize


    !-------------------------------------------------------------------------!
    !            Calculate number of eigenvalues and Lanczos vectors          !
    !-------------------------------------------------------------------------!

    subroutine ncv_from_nev(unit, parameters, thresh, exact, nevext, nestext, ncv0, dim, full, nev, ncv, nest)

        implicit none
        integer, intent(in) :: unit, thresh, exact, nevext, nestext, ncv0
        integer(kind=8), intent(in) :: dim
        character(len=*), intent(in) :: parameters
        integer, intent(out) :: full, nev, ncv, nest
        character*512 :: filename

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

    end subroutine ncv_from_nev

    !-------------------------------------------!
    !            Exact diagonalization          !
    !-------------------------------------------!

    subroutine exactdiag(eigvec, dim, mat, evals)
    
        implicit none
        logical, intent(in) :: eigvec
        integer(kind=8), intent(in) :: dim
        double precision, intent(inout) :: mat(dim, dim)
        double precision, allocatable, intent(out) :: evals(:)

        double precision, allocatable, save :: matloc(:,:)
        double precision, allocatable, save :: evalsloc(:)
        integer, parameter :: lwmax = 100000
        integer, save :: lda = 0, info = 0, lwork = 0
        character, save :: jobz
        double precision, allocatable, save :: work(:)

        external :: dsyev, dsyev_2stage

        !$OMP THREADPRIVATE (matloc, evalsloc, lda, info, lwork, jobz, work)
        lda   = 0
        info  = 0
        lwork = 0

        if (eigvec) then
            jobz = 'V'
        else
            jobz = 'N'
        end if

        if (allocated(evals)) deallocate(evals)
        allocate(evals(dim))
        evals = 0.0d0

        if (allocated(evalsloc)) deallocate(evalsloc)
        if (allocated(matloc)) deallocate(matloc)
        allocate(evalsloc(dim))
        allocate(matloc(dim, dim))
        evalsloc = 0.0d0
        matloc = mat

        if (dim > 1) then
            lda = dim
            if (allocated(work)) deallocate(work)
            allocate(work(lwmax))
            lwork = -1
            if (jobz == 'V') then
                !$omp critical
                call dsyev(jobz, 'U', dim, matloc, lda, evalsloc, work, lwork, info)
                !$omp end critical
            else if (jobz == 'N') then
                call dsyev_2stage ('N', 'U', dim, matloc, lda, evalsloc, work, lwork, info)
            end if
            !lwork = min(int(work(1)), lwmax)
            lwork = int(work(1))
            if (lwork < max(1, 6*(dim-1) ) ) lwork = max(1, 6*(dim-1) )
            if (allocated(work)) deallocate(work)
            allocate(work(lwork))
            work = 0.0d0

            !print*, 'Running diagoalization routine ...'
            if (jobz == 'V') then
                !$omp critical
                call dsyev(jobz, 'U', dim, matloc, lda, evalsloc, work, lwork, info)
                !$omp end critical
            else if (jobz == 'N') then
                call dsyev_2stage ('N', 'U', dim, matloc, lda, evalsloc, work, lwork, info)
            end if

            !print*, 'Finished diagoalization routine.'
            if (info .gt. 0 ) THEN
                write(*,*)'The algorithm failed to compute eigenvalues.'
                stop
            end if
            if (allocated(work)) deallocate(work)
        else
            evalsloc = matloc(1,1)
        end if
        ! print*, 'Eigenvalues '
        ! do i = 1, dim 
        !     print*, evalsloc(i)
        ! end do 
        evals = evalsloc
        mat   = matloc
        if(allocated(evalsloc)) deallocate(evalsloc)
        if(allocated(matloc)) deallocate(matloc)

        return

    end subroutine exactdiag

    subroutine exactdiag_c(eigvec, dim, mat, evals)
        
        implicit none
        logical, intent(in) :: eigvec
        integer(kind=8), intent(in) :: dim
        double complex, intent(inout) :: mat(dim, dim)
        double precision, allocatable, intent(out) :: evals(:)

        double complex, allocatable, save :: matloc(:,:)
        double precision, allocatable, save :: evalsloc(:)
        integer, parameter :: lwmax = 100000
        integer, save :: lda = 0, info = 0, lwork = 0
        character, save :: jobz
        double precision, allocatable, save :: rwork(:)
        double complex, allocatable, save :: work(:)

        external :: zheev

        !$OMP THREADPRIVATE (matloc, evalsloc, lda, info, lwork, jobz, work)
        lda   = 0
        info  = 0
        lwork = 0

        if (eigvec) then
            jobz = 'V'
        else
            jobz = 'N'
        end if

        if (allocated(evals)) deallocate(evals)
        allocate(evals(dim))
        evals = 0.0d0

        if (allocated(evalsloc)) deallocate(evalsloc)
        if (allocated(matloc)) deallocate(matloc)
        allocate(evalsloc(dim))
        allocate(matloc(dim, dim))
        evalsloc = 0.0d0
        matloc = mat

        if (dim > 1) then
            lda = dim
            if (allocated(work)) deallocate(work)
            allocate(work(lwmax))
            if (allocated(rwork)) deallocate(rwork)
            allocate(rwork(max(1, 3*dim-2)))
            lwork = -1
            if (jobz == 'V') then
                !$omp critical
                call zheev(jobz, 'U', dim, matloc, lda, evalsloc, work, lwork, rwork, info)
                !$omp end critical
            ! else if (jobz == 'N') then
            !     call dsyev_2stage ('N', 'U', dim, matloc, lda, evalsloc, work, lwork, info)
            end if
            !lwork = min(int(work(1)), lwmax)
            lwork = int(work(1))
            if (lwork < max(1, 6*(dim-1) ) ) lwork = max(1, 6*(dim-1) )
            if (allocated(work)) deallocate(work)
            allocate(work(lwork))
            work  = 0.0d0
            rwork = 0.d0 

            !print*, 'Running diagoalization routine ...'
            if (jobz == 'V') then
                !$omp critical
                call zheev(jobz, 'U', dim, matloc, lda, evalsloc, work, lwork, rwork, info)
                !$omp end critical
            ! else if (jobz == 'N') then
            !     call dsyev_2stage ('N', 'U', dim, matloc, lda, evalsloc, work, lwork, info)
            end if

            !print*, 'Finished diagoalization routine.'
            if (info .gt. 0 ) THEN
                write(*,*)'The algorithm failed to compute eigenvalues.'
                stop
            end if
            if (allocated(work)) deallocate(work)
            if (allocated(rwork)) deallocate(rwork)
        else
            evalsloc = matloc(1,1)
        end if
        ! print*, 'Eigenvalues '
        ! do i = 1, dim 
        !     print*, evalsloc(i)
        ! end do 
        evals = evalsloc
        mat   = matloc
        if(allocated(evalsloc)) deallocate(evalsloc)
        if(allocated(matloc)) deallocate(matloc)

        return

    end subroutine exactdiag_c


    !---------------------------------------------------!
    !            Sparse diagonalization  (ARPACK)       !
    !---------------------------------------------------!
    ! This subroutine uses ARPACK to perform Lanczos diagonalization.

    subroutine lanczos_d(threads, dim, nnz, nev, ncv, nest, mode, rvec, ham, rc, evals, evecs)

        implicit none
        ! save
        integer(kind=8), intent(in) :: dim
        integer, intent(in) :: threads, nev, ncv, nest, nnz
        integer, intent(in) :: rc(nnz,2)
        double precision, intent(in) :: ham(nnz)
        character(len=*), intent(in) :: mode
        logical, intent(in)          :: rvec
        double precision, allocatable, intent(out) :: evals(:), evecs(:,:)

        integer(kind=8), save :: maxn, n
        integer, save :: maxnev, maxncv, ldv, j
        integer, save :: ishfts, lworkl, maxitr, mode1, nconv, iparam(11), ipntr(11), ido, ierr, info
        integer, allocatable, save :: jao(:), iao(:)
        double precision, save :: zero = 0.0D+00
        double precision, save :: sigma, tol
        double precision, allocatable, save :: ax(:), workl(:), workd(:), v(:,:), resid(:), d(:,:), ao(:)
        character, save :: bmat*1, which*2
        logical, allocatable, save :: selector(:)
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

        if (dim == 1) then
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

            if ( n .gt. maxn ) then
                print *, ' ERROR: N is greater than MAXN '
                go to 9870
            else if ( nev .gt. maxnev ) then
                print *, ' ERROR: NEV is greater than MAXNEV '
                go to 9870
            else if ( ncv .gt. maxncv ) then
                print *, ' ERROR: NCV is greater than MAXNCV '
                go to 9870
            end if

            lworkl = ncv * ( ncv + 8 )
            tol    = zero
            info   = 0
            ido    = 0
            ishfts = 1
            maxitr = 300! 5000! 1000 !300
            mode1  = 1
            iparam = 0
            cntr   = 0 

            iparam(1) = ishfts
            iparam(3) = maxitr
            iparam(7) = mode1
            iparam(4) = 1 
            !MAIN LOOP (Reverse communication loop)! Repeatedly call DSAUPD and take actions indicated by parameter !IDO until convergence is indicated or MAXITR is exceeded.
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
            
            call coocsr(dim, nnz, ham, rc(1:nnz,1), rc(1:nnz,2), ao, jao, iao) 
            
        !!$omp critical
            do
            !Sparse eigensolver routine
            call dsaupd ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
            if ( ido /= -1 .and. ido /= 1 ) then
                exit
            end if
            !Vector-sparse matrix multiplication                 
            cntr = cntr + 1
            call pamux (threads, dim, workd(ipntr(1)), workd(ipntr(2)), ao, jao, iao)
            end do
        !!$omp end critical
            print*,'Number of MV-multiplications',cntr
            call dsaupderrormessage(info)
            
            !Extract eigenvalues and -vectors
            ! !$omp critical
            call dseupd ( rvec, 'A', selector, d, v, ldv, sigma, bmat, n, which, nev, tol, &
                        resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, ierr )
            call dseupderrormessage(ierr)
            ! !$omp end critical
            if (ierr /= 0) then
            print *, 'DSSIMP - Fatal error! Error with DSEUPD, IERR = ', ierr, 'Check the documentation of DSEUPD.'
            else
            nconv =  iparam(5)
            do j = 1, nconv
                call pamux (threads, dim, v(1,j), ax, ao, jao, iao)
                call daxpy (n, -d(j,1), v(1,j), 1, ax, 1 )
                d(j,2) = dnrm2 (n, ax, 1)
                d(j,2) = d(j,2) / abs (d(j,1))
            end do
            call dmout (6, nconv, 2, d, maxncv, -6, 'Ritz values are named relative residuals' ) !Display computed residuals. 6: Output to screen Write(6, #internalnumber)! nconv: number of rows in the matrix d! 2: Number of columns in matrix d! maxncv: Leading dimension of the matrix data! -6: print the matrix d with iabs(-6) decimal digits per number! Use formatting indexed by -6 to print A
            end if
            


    !     10 continue
        ! !c
        ! !c        %---------------------------------------------%
        ! !c        | Repeatedly call the routine DSAUPD and take |
        ! !c        | actions indicated by parameter IDO until    |
        ! !c        | either convergence is indicated or maxitr   |
        ! !c        | has been exceeded.                          |
        ! !c        %---------------------------------------------%
        ! !c

        !     call dsaupd ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
        !     if (ido .eq. -1 .or. ido .eq. 1) then
        ! !c
        ! !c           %--------------------------------------%
        ! !c           | Perform matrix vector multiplication |
        ! !c           |              y <--- OP*x             |
        ! !c           | The user should supply his/her own   |
        ! !c           | matrix vector multiplication routine |
        ! !c           | here that takes workd(ipntr(1)) as   |
        ! !c           | the input, and return the result to  |
        ! !c           | workd(ipntr(2)).                     |
        ! !c           %--------------------------------------%

        !         call pamux (threads, dim, workd(ipntr(1)), workd(ipntr(2)), ao, jao, iao)

        ! !c           %-----------------------------------------%
        ! !c           | L O O P   B A C K to call DSAUPD again. |
        ! !c           %-----------------------------------------%

        !         go to 10

        !     end if

        ! !c     %----------------------------------------%
        ! !c     | Either we have convergence or there is |
        ! !c     | an error.                              |
        ! !c     %----------------------------------------%

        !     if ( info .lt. 0 ) then
        ! !c
        ! !c        %--------------------------%
        ! !c        | Error message. Check the |
        ! !c        | documentation in DSAUPD. |
        ! !c        %--------------------------%
        ! !c
        !         print *, ' '
        !         print *, ' Error with dsaupd, info = ', info
        !         print *, ' Check documentation in dsaupd '
        !         print *, ' '
        !     else
        ! !c
        ! !c        %-------------------------------------------%
        ! !c        | No fatal errors occurred.                 |
        ! !c        | Post-Process using DSEUPD.                |
        ! !c        |                                           |
        ! !c        | Computed eigenvalues may be extracted.    |
        ! !c        |                                           |
        ! !c        | Eigenvectors may also be computed now if  |
        ! !c        | desired.  (indicated by rvec = .true.)    |
        ! !c        %-------------------------------------------%
        ! !c

        !         call dseupd ( rvec, 'All', selector, d, v, ldv, sigma, bmat, n, which, nev, tol, &
        !                       resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, ierr )

        ! !c        %----------------------------------------------%
        ! !c        | Eigenvalues are returned in the first column |
        ! !c        | of the two dimensional array D and the       |
        ! !c        | corresponding eigenvectors are returned in   |
        ! !c        | the first NEV columns of the two dimensional |
        ! !c        | array V if requested.  Otherwise, an         |
        ! !c        | orthogonal basis for the invariant subspace  |
        ! !c        | corresponding to the eigenvalues in D is     |
        ! !c        | returned in V.                               |
        ! !c        %----------------------------------------------%
        ! !c
    
        !         if ( ierr .ne. 0) then
        ! !c
        ! !c            %------------------------------------%
        ! !c            | Error condition:                   |
        ! !c            | Check the documentation of DSEUPD. |
        ! !c            %------------------------------------%
        ! !c
        !             print *, ' '
        !             print *, ' Error with dseupd, info = ', ierr
        !             print *, ' Check the documentation of dseupd. '
        !             print *, ' '

        !         else

        !             nconv =  iparam(5)
        !             do j = 1, nconv

        ! !c
        ! !c               %---------------------------%
        ! !c               | Compute the residual norm |
        ! !c               |                           |
        ! !c               |   ||  A*x - lambda*x ||   |
        ! !c               |                           |
        ! !c               | for the NCONV accurately  |
        ! !c               | computed eigenvalues and  |
        ! !c               | eigenvectors.  (iparam(5) |
        ! !c               | indicates how many are    |
        ! !c               | accurate to the requested |
        ! !c               | tolerance)                |
        ! !c               %---------------------------%
        ! !c
        !                 call pamux (threads, dim, v(1,j), ax, ao, jao, iao)
        !                 call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
        !                 d(j,2) = dnrm2(n, ax, 1)
        !                 d(j,2) = d(j,2) / abs(d(j,1))

        !             end do
        ! !c
        ! !c            %-------------------------------%
        ! !c            | Display computed residuals    |
        ! !c            %-------------------------------%
        ! !c
        !             call dmout(6, nconv, 2, d, maxncv, -6, &
        !                 'Ritz values and relative residuals')
        !             end if
    !         end if

                evals(1:nev) = d(1:nev,1)
                if (rvec) then
                    do j= 1, min(iparam(5), nest) 
                        evecs(1:dim, j) = v(1:dim,j)
                    end do
                end if
                
                9870 continue

            end if


            100 format(1000(F40.30))

            ! if(allocated(workd))    deallocate(workd)
            ! if(allocated(ax))       deallocate(ax)
            ! if(allocated(d))        deallocate(d)
            ! if(allocated(resid))    deallocate(resid)
            ! if(allocated(selector)) deallocate(selector)
            ! if(allocated(workl))    deallocate(workl)
            ! if(allocated(v))        deallocate(v)
            ! if(allocated(jao))      deallocate(jao)
            ! if(allocated(iao))      deallocate(iao)
            ! if(allocated(ao))       deallocate(ao)

            return 

    end subroutine lanczos_d

    subroutine lanczos_d_otf(dir, params, unit, threads, dim, sites, nbonds, nnnbonds, nev, ncv, nest, t1, v1, v2, mode, rvec, basis, bsites, hexsites, evals, evecs)

        implicit none

        integer(kind=8), intent(in) :: dim, basis(dim)
        integer, intent(in) :: unit, threads, sites, nbonds, nnnbonds, nev, ncv, nest
        integer, intent(in) :: bsites(2, nbonds), hexsites(2, nnnbonds)
        double precision, intent(in) :: t1, v1, v2
        character(len=*), intent(in) :: dir, params, mode
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

        !!$omp threadprivate(maxn, n, maxnev, maxncv, ldv, j, ishfts, lworkl, maxitr, mode1, nconv, iparam, ipntr, ido, ierr, info, zero, sigma, tol, ax, workl, workd, v, resid, d, bmat, which, selector)


        if(allocated(evals)) deallocate(evals)
        if(allocated(evecs)) deallocate(evecs)
        allocate(evals(nev))
        allocate(evecs(dim, nest))
        evals = 0.d0 
        evecs = 0.d0    

        if (dim == 1) then
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

            if ( n .gt. maxn ) then
                print *, ' ERROR: N is greater than MAXN '
                go to 9870
            else if ( nev .gt. maxnev ) then
                print *, ' ERROR: NEV is greater than MAXNEV '
                go to 9870
            else if ( ncv .gt. maxncv ) then
                print *, ' ERROR: NCV is greater than MAXNCV '
                go to 9870
            end if

            lworkl = ncv * ( ncv + 8 )
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
            !MAIN LOOP (Reverse communication loop)! Repeatedly call DSAUPD and take actions indicated by parameter !IDO until convergence is indicated or MAXITR is exceeded.
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
        !!$omp critical
            do
            !Sparse eigensolver routine
            call dsaupd( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
            if ( ido /= -1 .and. ido /= 1 ) then
                exit
            end if
            cntr = cntr + 1
            if(mod(cntr, 100) == 0) print*,'cntr',cntr
            !Vector-sparse matrix multiplication                 
            call pamux2(threads, dim, sites, nbonds, nnnbonds, t1, v1, v2, basis, bsites, hexsites, workd(ipntr(1)), workd(ipntr(2)))
            end do
        !!$omp end critical
            print*,'Number of MV-multiplications',cntr
            call dsaupderrormessage(info)
            
            !Extract eigenvalues and -vectors
            ! !$omp critical
            call dseupd( rvec, 'A', selector, d, v, ldv, sigma, bmat, n, which, nev, tol, &
                        resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, ierr )
            call dseupderrormessage(ierr)
            ! !$omp end critical
            if (ierr /= 0) then
            print *, 'DSSIMP - Fatal error! Error with DSEUPD, IERR = ', ierr, 'Check the documentation of DSEUPD.'
            else
            nconv =  iparam(5)
            do j = 1, nconv
                call pamux2(threads, dim, sites, nbonds, nnnbonds, t1, v1, v2, basis, bsites, hexsites, v(1,j), ax)
                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1 )
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
            end do
            call dmout(6, nconv, 2, d, maxncv, -6, 'Ritz values are named relative residuals' ) !Display computed residuals. 6: Output to screen Write(6, #internalnumber)! nconv: number of rows in the matrix d! 2: Number of columns in matrix d! maxncv: Leading dimension of the matrix data! -6: print the matrix d with iabs(-6) decimal digits per number! Use formatting indexed by -6 to print A
            end if

            evals(1:nev) = d(1:nev,1)

            if (rvec) then
                do j= 1, min(iparam(5), nest) 
                    evecs(1:dim, j) = v(1:dim,j)
                end do
            end if
            
            9870 continue

        end if


        100 format(1000(F40.30))

        ! if(allocated(workd))    deallocate(workd)
        ! if(allocated(ax))       deallocate(ax)
        ! if(allocated(d))        deallocate(d)
        ! if(allocated(resid))    deallocate(resid)
        ! if(allocated(selector)) deallocate(selector)
        ! if(allocated(workl))    deallocate(workl)
        ! if(allocated(v))        deallocate(v)
        ! if(allocated(jao))      deallocate(jao)
        ! if(allocated(iao))      deallocate(iao)
        ! if(allocated(ao))       deallocate(ao)


        return 

    end subroutine lanczos_d_otf

    !-------------------------------------------------!
    !            Sparse diagonalization  (MKL)        !
    !-------------------------------------------------!
    ! This subroutine uses Intel MKL Sparse BLAS to perform Lanczos diagonalization
    ! It is a wrapper for MKL_SPARSE_S_EV and MKL_SPARSE_EE routines.
    ! It is used to compute the eigenvalues and eigenvectors of a sparse matrix.

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
        integer, intent(in) :: nnz
        ! integer, intent(in) :: row_indx(nnz), col_indx(nnz)   
        integer, intent(in) :: rc(nnz,2)
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
        parameter  (WHICH='S') !L - for largest eigenvalues to find
        integer, allocatable :: pm(:)
        integer :: n, K, info

        double precision, allocatable :: E(:)
        double precision, allocatable :: sigma(:)
        double precision, allocatable :: X(:,:)
        double precision, allocatable :: res(:)

    
        !!$omp threadprivate(i, descrA, cooA, csrA, pm, n, k, info, E, sigma, X, res)

        
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
        n = int(dim, 4) 
        print *, 'Sparse matrix size', n
        !
        !        Task 0. Call MKL_SPARSE_C_CREATE_COO to create matrix handle 
        !      

        ! info = mkl_sparse_d_create_coo (cooA, SPARSE_INDEX_BASE_ZERO, n, n, nnz, row_indx, col_indx, values)
        info = mkl_sparse_d_create_coo (cooA, SPARSE_INDEX_BASE_ONE, n, n, nnz, rc(1:nnz,1), rc(1:nnz,2), values)

        !
        !        Task 1. Call MKL_SPARSE_C_CREATE_CSR to create matrix handle
        !      

        info = mkl_sparse_convert_csr (cooA, SPARSE_OPERATION_NON_TRANSPOSE, csrA)
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
        !2 - Subspace Iteration technique based on FEAST algorithm (for eigenvalues in specific interval)
        pm(4) = ncv 
        ! This parameter is referenced only for Krylov-Schur Method. It indicates the number of Lanczos/Arnoldi vectors (NCV) generated at each iteration.
        ! This parameter must be less than or equal to size of matrix and greater than number of eigenvalues (k0) to be computed. If unspecified, NCV is set to be at least 1.5 times larger than k0. 

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
        print  *,' ITERATION EXIT REASON (0 = converged) ', pm(10)
        
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

    subroutine lanczos_c( threads, dim, nev, ncv, nst, mode, rvec, nnz, ham, rc, evals, evecs )
        
        implicit none

        integer, intent(in) :: threads, nev, ncv, nst, nnz
        integer(kind=8), intent(in) :: dim
        integer, intent(in) :: rc(nnz,2)
        double complex, intent(in) :: ham(nnz)   
        character*2, intent(in) :: mode
        logical, intent(in) :: rvec
        double precision, allocatable, intent(out) :: evals(:)
        double complex, allocatable, intent(out) :: evecs(:, :)

        integer :: printing = 1 
        integer(kind=8), save :: maxn = 0, maxnev = 0, maxncv = 0, ldv = 0
        integer(kind=8), save :: n = 0, nx = 0
        integer, save :: buffer = 0 
        integer, save :: j = 0, iparam(11), ipntr(14)
        integer, save :: ido = 0, ishfts = 0, lworkl = 0, info = 0, maxitr = 0, mode1 = 0, nconv = 0, ierr = 0    !VARIABLES FOR DIAGONALIZATION ROUTINE
        integer, allocatable, save :: jao(:), iao(:), jao_ord(:)
        double precision, save :: tol = 0
        double precision, allocatable, save :: rwork(:), rd(:,:)
        double complex, save :: sigma = 0    
        double complex, allocatable, save :: ao(:), ao_ord(:)
        double complex, allocatable, save :: ax(:), d(:), v(:,:), workd(:), workev(:), resid(:), workl(:)!VARIABLES FOR DIAGONALIZATION ROUTINE
        character, save :: bmat*1, which*2
        logical, allocatable, save :: selector(:)

        intrinsic :: abs
        !c
        !c     %-----------------------------%
        !c     | BLAS & LAPACK routines used |
        !c     %-----------------------------%
        !c

        double precision :: dznrm2 , dlapy2
        external :: dznrm2 , zaxpy , dlapy2, znaupd, zneupd, dmout
        double precision, external :: dnrm2

        ! c  WHICH   Character*2.  (INPUT)
        !c
        !c     %-----------------------%
        !c     | Executable Statements |
        !c     %-----------------------%
        !c
        !c     %--------------------------------------------------%
        !c     | The number NX is the number of interior points   |
        !c     | in the discretization of the 2-dimensional       |
        !c     | convection-diffusion operator on the unit        |
        !c     | square with zero Dirichlet boundary condition.   |
        !c     | The number N(=NX*NX) is the dimension of the     |
        !c     | matrix.  A standard eigenvalue problem is        |
        !c     | solved (BMAT = 'I').  NEV is the number of       |
        !c     | eigenvalues to be approximated.  The user can    |
        !c     | modify NX, NEV, NCV, WHICH to solve problems of  |
        !c     | different sizes, and to get different parts of   |
        !c     | the spectrum.  However, The following            |
        !c     | conditions must be satisfied:                    |
        !c     |                   N <= MAXN                      |
        !c     |                 NEV <= MAXNEV                    |
        !c     |           NEV + 2 <= NCV <= MAXNCV               |
        !c     %--------------------------------------------------%
        !c

        !$omp threadprivate(ao, jao, iao, bmat, which, j, iparam, sigma, ipntr, ido, ishfts, lworkl, info, maxitr, mode1, nconv, ierr, selector, ax, d, v, workd, workev, resid, workl, rwork, rd)
        
        
        if (dim == 1) then
            evecs(1,1) = 1
            evals(1) = ham(1)
            print*, ''
            write(*,"('Only eigenvalue is ',f8.4)") ham(1)
        else
        
        buffer = 50 
        ! maxn = 10 + dim * dim !debug
        maxn   = buffer + dim !* dim !debug
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


        nx = dim !debug
        ! n  = nx * nx !debug
        n  = nx !* nx !debug

        if ( n .gt. maxn ) then
            print *, ' ERROR with _NDRV1: N is greater than MAXN '
            go to 9870
        else if ( nev .gt. maxnev ) then
            print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
            go to 9870
        else if ( ncv .gt. maxncv ) then
            print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
            go to 9870
        end if
        bmat  = 'I'
        which = 'SR'
        !c
        !c     %---------------------------------------------------%
        !c     | The work array WORKL is used in ZNAUPD  as         |
        !c     | workspace.  Its dimension LWORKL is set as        |
        !c     | illustrated below.  The parameter TOL determines  |
        !c     | the stopping criterion. If TOL<=0, machine        |
        !c     | precision is used.  The variable IDO is used for  |
        !c     | reverse communication, and is initially set to 0. |
        !c     | Setting INFO=0 indicates that a random vector is  |
        !c     | generated to start the ARNOLDI iteration.         |
        !c     %---------------------------------------------------%
        !c
        lworkl = 3*ncv**2+5*ncv
        tol    = 0
        ido    = 0
        info   = 0
        nconv  = 0


        !c
        !c     %---------------------------------------------------%
        !c     | This program uses exact shift with respect to     |
        !c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
        !c     | IPARAM(3) specifies the maximum number of Arnoldi |
        !c     | iterations allowed.  Mode 1 of ZNAUPD  is used     |
        !c     | (IPARAM(7) = 1). All these options can be changed |
        !c     | by the user. For details see the documentation in |
        !c     | ZNAUPD .                                           |
        !c     %---------------------------------------------------%
        !c
        ishfts = 1
        maxitr = 300 !5000
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
        
        call ccoocsr(dim, nnz, ham, rc(1:nnz,1), rc(1:nnz,2), ao, jao, iao)

        ! if(allocated(jao_ord)) deallocate(jao_ord)
        ! if(allocated(ao_ord)) deallocate(ao_ord)
        ! allocate(jao_ord(nnz))
        ! allocate(ao_ord(nnz))
        ! ao_ord  = 0.d0 
        ! jao_ord = 0
        ! call sortCols(dim, nnz, iao, jao, ao, jao_ord, ao_ord)
        ! ao = ao_ord 
        ! jao = jao_ord 
        !$omp critical
        do
            !Sparse eigensolver routine
            call znaupd ( ido, bmat, n, which, nev, tol, resid, ncv, &
                        v, ldv, iparam, ipntr, workd, workl, lworkl, &
                        rwork, info )
            if ( ido /= -1 .and. ido /= 1 ) then
                exit
            end if
            !Vector-sparse matrix multiplication
            ! call mkl_zcsrmv(transa, m, k, alpha, matdescra, val, indx, pntrb, pntre, x, beta, y)
            call cpamux (threads, nx, workd(ipntr(1)), workd(ipntr(2)), ao, jao, iao)
        end do
        !$omp end critical
    
        ! !c
        ! !c     %-------------------------------------------%
        ! !c     | M A I N   L O O P (Reverse communication) |
        ! !c     %-------------------------------------------%
        ! !c
        ! 10   continue
        ! !c
        ! !c        %---------------------------------------------%
        ! !c        | Repeatedly call the routine ZNAUPD  and take |
        ! !c        | actions indicated by parameter IDO until    |
        ! !c        | either convergence is indicated or maxitr   |
        ! !c        | has been exceeded.                          |
        ! !c        %---------------------------------------------%
        ! !c
        ! call znaupd ( ido, bmat, n, which, nev, tol, resid, ncv, &
        !             v, ldv, iparam, ipntr, workd, workl, lworkl, &
        !             rwork, info )
        ! !c
        ! !c           %-------------------------------------------%
        ! !c           | Perform matrix vector multiplication      |
        ! !c           |                y <--- OP*x                |
        ! !c           | The user should supply his/her own        |
        ! !c           | matrix vector multiplication routine here |
        ! !c           | that takes workd(ipntr(1)) as the input   |
        ! !c           | vector, and return the matrix vector      |
        ! !c           | product to workd(ipntr(2)).               |
        ! !c           %-------------------------------------------%
        ! !c

        !     call cpamux (threads, nx, workd(ipntr(1)), workd(ipntr(2)), ao, jao, iao)

        ! !c
        ! !c           %-----------------------------------------%
        ! !c           | L O O P   B A C K to call ZNAUPD  again. |
        ! !c           %-----------------------------------------%
        ! !c
        !     go to 10

        ! end if


        !c
        !c     %----------------------------------------%
        !c     | Either we have convergence or there is |
        !c     | an error.                              |
        !c     %----------------------------------------%
        !c
        if ( info .lt. 0 ) then
            !c
            !c        %--------------------------%
            !c        | Error message, check the |
            !c        | documentation in ZNAUPD   |
            !c        %--------------------------%
            !c
            print *, ' '
            print *, ' Error with znaupd, info = ', info
            print *, ' Check the documentation of _naupd'
            print *, ' '

        else
        !c
        !c        %-------------------------------------------%
        !c        | No fatal errors occurred.                 |
        !c        | Post-Process using ZNEUPD .                |
        !c        |                                           |
        !c        | Computed eigenvalues may be extracted.    |
        !c        |                                           |
        !c        | Eigenvectors may also be computed now if  |
        !c        | desired.  (indicated by rvec = .true.)    |
        !c        %-------------------------------------------%
        !c


            call zneupd (rvec, 'A', selector, d, v, ldv, sigma, &
                        workev, bmat, n, which, nev, tol, resid, ncv, &
                        v, ldv, iparam, ipntr, workd, workl, lworkl, &
                        rwork, ierr)
        !c
        !c        %----------------------------------------------%
        !c        | Eigenvalues are returned in the one          |
        !c        | dimensional array D.  The corresponding      |
        !c        | eigenvectors are returned in the first NCONV |
        !c        | (=IPARAM(5)) columns of the two dimensional  |
        !c        | array V if requested.  Otherwise, an         |
        !c        | orthogonal basis for the invariant subspace  |
        !c        | corresponding to the eigenvalues in D is     |
        !c        | returned in V.                               |
        !c        %----------------------------------------------%
        !c


            if ( ierr .ne. 0) then
                !c
                !c           %------------------------------------%
                !c           | Error condition:                   |
                !c           | Check the documentation of ZNEUPD  |
                !c           %------------------------------------%
                !c
                print *, ' '
                print *, ' Error with zneupd, info = ', ierr
                print *, ' Check the documentation of _neupd. '
                print *, ' '

            else

                nconv = iparam(5)
                do 20 j = 1, nconv
                    !c
                    !c               %---------------------------%
                    !c               | Compute the residual norm |
                    !c               |                           |
                    !c               |   ||  A*x - lambda*x ||   |
                    !c               |                           |
                    !c               | for the NCONV accurately  |
                    !c               | computed eigenvalues and  |
                    !c               | eigenvectors.  (iparam(5) |
                    !c               | indicates how many are    |
                    !c               | accurate to the requested |
                    !c               | tolerance)                |
                    !c               %---------------------------%
                    !c

                    call cpamux (threads, nx, v(1,j), ax, ao, jao, iao)
                    call zaxpy (n, -d(j), v(1,j), 1, ax, 1)
                    rd(j,1) = dble (d(j))
                    rd(j,2) = dimag (d(j))
                    rd(j,3) = dznrm2 (n, ax, 1)
                    rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
                20 continue
                !c
                !c            %-----------------------------%
                !c            | Display computed residuals. |
                !c            %-----------------------------%
                !c
                if(printing == 1) then 
                    call dmout (6, nconv, 3, rd, maxncv, -6, &
                                'Ritz values (Real, Imag) and relative residuals')
                end if 
            end if


                ! c
                ! c        %-------------------------------------------%
                ! c        | Print additional convergence information. |
                ! c        %-------------------------------------------%
                ! c
            if ( info .eq. 1) then
                print *, ' '
                print *, ' Maximum number of iterations reached.'
                print *, ' '
            else if ( info .eq. 3) then
                print *, ' '
                print *, ' No shifts could be applied during implicit', &
                        ' Arnoldi update, try increasing NCV.'
                print *, ' '
            end if

            ! print *, ' '
            ! print *, '_NDRV1'
            ! print *, '====== '
            ! print *, ' '
            ! print *, ' Size of the matrix is ', n
            ! print *, ' The number of Ritz values requested is ', nev
            ! print *, ' The number of Arnoldi vectors generated', &
            !          ' (NCV) is ', ncv
            ! print *, ' What portion of the spectrum: ', which
            ! print *, ' The number of converged Ritz values is ', &
            !            nconv
            ! print *, ' The number of Implicit Arnoldi update', &
            !          ' iterations taken is ', iparam(3)
            ! print *, ' The number of OP*x is ', iparam(9)
            ! print *, ' The convergence criterion is ', tol
            ! print *, ' '

        end if

        evals(1:nev) = dble(d(1:nev))
        evecs = 0 
        
        if (rvec) then
            do j = 1, nst
                evecs(1:dim, j) = v(1:dim,j)
            end do
        end if

        
        end if


        !c
        !c     %---------------------------%
        !c     | Done with program zndrv1 . |
        !c     %---------------------------%
        !c

        !c          Error flag for ZNAUPD on output.
        !c          =  0: Normal exit.
        !c          =  1: Maximum number of iterations taken.
        !c                All possible eigenvalues of OP has been found. IPARAM(5)
        !c                returns the number of wanted converged Ritz values.
        !c          =  2: No longer an informational error. Deprecated starting
        !c                with release 2 of ARPACK.
        !c          =  3: No shifts could be applied during a cycle of the
        !c                Implicitly restarted Arnoldi iteration. One possibility
        !c                is to increase the size of NCV relative to NEV.
        !c                See remark 4 below.
        !c          = -1: N must be positive.
        !c          = -2: NEV must be positive.
        !c          = -3: NCV-NEV >= 2 and less than or equal to N.
        !c          = -4: The maximum number of Arnoldi update iteration
        !c                must be greater than zero.
        !c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
        !c          = -6: BMAT must be one of 'I' or 'G'.
        !c          = -7: Length of private work array is not sufficient.
        !c          = -8: Error return from LAPACK eigenvalue calculation;
        !c          = -9: Starting vector is zero.
        !c          = -10: IPARAM(7) must be 1,2,3.
        !c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
        !c          = -12: IPARAM(1) must be equal to 0 or 1.
        !c          = -9999: Could not build an Arnoldi factorization.
        !c                   User input error highly likely.  Please
        !c                   check actual array dimensions and layout.
        !c                   IPARAM(5) returns the size of the current Arnoldi
        !c                   factorization.


        !    Error flag for ZNEUPD on output.
        !    0: Normal exit.
        !    1: The Schur form computed by LAPACK routine csheqr could not be reordered by LAPACK routine ztrsen.
        !       Re-enter subroutine zneupd with IPARAM(5) = NCV and increase the size of the array D to have dimension at least dimension NCV and allocate at least NCV columns for Z.
        !       NOTE: Not necessary if Z and V share the same space. Please notify the authors if this error occurs.
        !    -1: N must be positive.
        !    -2: NEV must be positive.
        !    -3: NCV-NEV >= 1 and less than or equal to N.
        !    -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'.
        !    -6: BMAT must be one of 'I' or 'G'.
        !    -7: Length of private work WORKL array is not sufficient.
        !    -8: Error return from LAPACK eigenvalue calculation. This should never happened.
        !    -9: Error return from calculation of eigenvectors. Informational error from LAPACK routine ztrevc.
        !    -10: IPARAM(7) must be 1, 2, 3.
        !    -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
        !    -12: HOWMANY = 'S' not yet implemented.
        !    -13: HOWMANY must be one of 'A' or 'P' if RVEC = .true.
        !    -14: ZNAUPD did not find any eigenvalues to sufficient accuracy.
        !    -15: ZNEUPD got a different count of the number of converged Ritz values than ZNAUPD got. This indicates the user probably made an error in passing data from ZNAUPD to ZNEUPD or that the data was modified before entering ZNEUPD.

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


    end subroutine lanczos_c


    !------------------------------!
    !            FEAST             !
    !------------------------------!

    subroutine cfeast(n, nnz, nev0, nest, rc, ham, Emin, Emax, nev, eigvals, eigvecs)
    
        !       1. The code calls  FEASTINIT  to define the default values for the input
        !          FEAST parameters.
        
        !       2. The  code  solves  the  standard eigenvalue  problem  Ax=ex   using 
        !          ZFEAST_HCSREV.
        
        !       3. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i) 
        !           are the expected eigenvalues  and E(i) are eigenvalues computed 
        !           with the help of ZFEAST_HCSREV().
        
        !       4. The code computes the maximum absolute value of the matrix  Y=(X')*X-I  
        !          where X is the matrix of eigenvectors computed with the help of 
        !          ZFEAST_HCSREV. ZGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
        ! *******************************************************************************
        USE MKL_SPBLAS   
        USE ISO_C_BINDING 
        implicit none
        !!!!!!!!!!!!!!!!! External function declarations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        external zgemm, zfeast_hcsrev, zfeast_hcsrgv, feastinit
        !!!!!!!!!!!!!!!!! Matrix declaration variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer(kind=8), intent(in) :: n
        integer, intent(in) :: nnz
        integer, intent(in) :: nev0
        integer, intent(in) :: rc(nnz,2)
        double complex, intent(in) :: ham(nnz)
        double precision, intent(in) :: Emin, Emax
        integer, intent(inout) :: nest
        integer, intent(out) :: nev 
        double precision, allocatable, intent(out) :: eigvals(:)
        double complex, allocatable, intent(out) :: eigvecs(:,:) 
        
        TYPE(MATRIX_DESCR) :: descrA     ! Sparse matrix descriptor
        !   CSR matrix representation
        TYPE(SPARSE_MATRIX_T) :: cooA, csrA    ! Structure with sparse matrix

        ! integer :: nrows, ncols 
        ! integer, allocatable :: csr_row_ptr(:), csr_col_indices(:)
        ! complex*16, allocatable :: csr_values(:)

        integer :: rows(n+1), cols(nnz)
        integer :: rows_b(n+1), cols_b(n)
        INTEGER(C_INT) :: indexing
        ! TYPE(C_PTR) :: rows_start_c, rows_end_c
        ! TYPE(C_PTR) :: val_c, col_c
        ! INTEGER, POINTER :: rows_start(:), rows_end(:)
        double complex :: val(nnz), val_b(n)
        double complex :: beta
        character*1 :: UPLO
        parameter   (UPLO='U')
        !!!!!!!!!!!!!!!!! Feast declaration variable !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! E - eigenvalues, X - eigenvectors, res - residual !!!!!!!!!!!!
        integer :: fpm(128)  
        double precision :: epsout
        integer :: loop
        integer :: L
        parameter   (L=10)
        integer :: M0,M,info
        double precision :: E(nev0)
        ! double precision :: E(n)
        !On output, the first m columns of x contain the orthonormal eigenvectors corresponding 
        !to the computed eigenvalues e, with the i-th column of x holding the eigenvector associated with e(i).   
        double complex :: X(n,nev0) 
        ! double complex :: X(n,n) 
        ! double precision :: res(n)
        double precision :: res(nev0)
        !!!!!!!!!!!!!!!!! Declaration of local variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! Eig - exact eigenvalues, R=|E-Eig|, Y=(X')*X-I !!!!!!!!!!!!!!!
        double precision :: Eig(nev0)
        double precision :: R(nev0)
        ! double precision      Eig(n)
        ! double precision      R(n)
        ! double complex  Y(n,n)
        integer :: i,j,n4
        integer ::  mkl, dim 
        integer ::  UTnnz = 0 
        integer ::  ldx, ldy
        double complex :: one, zero
        double precision :: smax, eigabs 

        integer, allocatable :: ib(:), jb(:), ic(:), jc(:)
        double complex, allocatable :: b(:), c(:)
        integer, allocatable :: UTcols(:), UTrows(:)
        double complex, allocatable :: UTval(:)

        print*,'Sparse matrix size',n
        print*,'Use Feast algorithm'
        M0=nev0
        M=int(0.5*M0)
        print *,'Search interval ', Emin,' ', Emax
        ldx=n
        ldy=n

        ! mkl = mkl_sparse_z_create_coo (cooA, SPARSE_INDEX_BASE_ONE, n, n, nnz, rc(1:nnz,1), rc(1:nnz,2), ham)
        ! mkl = mkl_sparse_convert_csr (cooA, SPARSE_OPERATION_NON_TRANSPOSE, csrA)
        ! mkl = mkl_sparse_z_export_csr (csrA, indexing, dim, dim, rows_start_c, rows_end_c, col_c, val_c)
        ! cols = rc(1:nnz,2) 
        ! val  = ham 
        ! call C_F_POINTER(rows_start_c, rows_start, [n])
        ! call C_F_POINTER(rows_end_c  , rows_end  , [n])
        ! rows(1:n) = rows_start(1:n)
        ! rows(n+1) = nnz + 1 
        if(allocated(b)) deallocate(b)
        if(allocated(jb)) deallocate(jb)
        if(allocated(ib)) deallocate(ib)
        allocate(b(n))
        allocate(jb(n))
        allocate(ib(n))
        val = (0.d0,0.d0)
        val_b = (0.d0,0.d0)
        cols = 0 
        rows = 0 
        cols_b = 0 
        rows_b = 0 

        do I = 1, N
            b(i)  = 1.d0
            ib(i) = i 
            jb(i) = i
        end do
        ! open(20,file='A.dat')
        !     do i = 1, nnz 
        ! write(20,'(f10.6,x,f10.8,x,i0,x,i0)') dble(ham(i)), aimag(ham(i)), rc(i,1), rc(i,2)
        !     end do 
        ! close(20)
        ! open(20,file='B.dat')
        ! do i = 1, n 
        ! write(20,'(f10.6,x,f10.8,x,i0,x,i0)') dble(b(i)), aimag(b(i)), ib(i), jb(i)
        !     end do 
        ! close(20)
        n4 = n 
        ! call coo_to_csr_hermitian(nnz, ham, rc(1:nnz,1), rc(1:nnz,2), csr_values, csr_row_ptr, csr_col_indices, nrows)
        call ccoocsr(n, nnz, ham, rc(1:nnz,1), rc(1:nnz,2), val, cols, rows)
        call ccoocsr(n, n4, b, ib, jb, val_b, cols_b, rows_b)

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
        
        ! call mkl_zcsradd('N', 0, 1, n, n, ham, cols, rows, beta, b, jb, ib, c, jc, ic, 1, info)
        call mkl_zcsrmultcsr('N', 0, 3, n, n, n, val, cols, rows, val_b, cols_b, rows_b, c, jc, ic, nnz, info)
        ! call mkl_zcsrmultcsr('N', 0, 1, n, n, n, csr_values, csr_col_indices, csr_row_ptr, val_b, cols_b, rows_b, c, jc, ic, nnz, info)
        deallocate(c, ic, jc)

        call HermitianUpperTriangularCSR(n, nnz, val, rows, cols, UTval, UTrows, UTcols, UTnnz)
            
        ! call mkl_zcsrmultcsr('T', 1, 0, n4, n4, n4, val, cols, rows, val, cols, rows, c, jc, ic, nnz, info)
        ! allocate(c(ic(n+1)-1))
        ! allocate(jc(ic(n+1)-1))
        ! jc = 0 
        ! c  = 0.d0 
        ! call mkl_zcsrmultcsr('T', 2, 0, n4, n4, n4, val, cols, rows, val, cols, rows, c, jc, ic, nnz, info)
        
        !
        !        Task 1. Call  FEASTINIT  to define the default values for the input
        !        FEAST parameters.
        !
        call feastinit(fpm)
        fpm(1)=1  !Specifies whether Extended Eigensolver routines print runtime status. 
        fpm(26)=1 !Specifies whether Extended Eigensolver routines check input matrices (applies to CSR format only). 
        fpm(27)=1 !Specifies whether Extended Eigensolver routines check input matrices (applies to CSR format only). 

        print *, ' Testing zfeast_hcsrev '
        !
        !         Task 2. Solve the standard eigenvalue problem Ax=ex.
        !
        ! call zfeast_hcsrev(UPLO,N,csr_values,csr_row_ptr,csr_col_indices,fpm,epsout,loop,&
        !    Emin,Emax,M0,E,X,M,res,info)

        call zfeast_hcsrev(UPLO,N,UTval,UTrows,UTcols,fpm,epsout,loop,&
        Emin,Emax,M0,E,X,M,res,info)
        ! call zfeast_hcsrev(UPLO,N,val,rows,cols,fpm,epsout,loop,&
        !    Emin,Emax,M0,E,X,M,res,info)
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
        !
        !         Task 3. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
        !         are the expected eigenvalues  and E(i) are eigenvalues computed
        !         with the help of ZFEAST_HCSREV().
        !
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

    end subroutine cfeast

    subroutine dfeast(n, nnz, nev0, nest, rc, ham, Emin, Emax, nev, eigvals, eigvecs)
    
        !       1. The code calls  FEASTINIT  to define the default values for the input
        !          FEAST parameters.
        
        !       2. The  code  solves  the  standard eigenvalue  problem  Ax=ex   using 
        !          ZFEAST_HCSREV.
        
        !       3. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i) 
        !           are the expected eigenvalues  and E(i) are eigenvalues computed 
        !           with the help of ZFEAST_HCSREV().
        
        !       4. The code computes the maximum absolute value of the matrix  Y=(X')*X-I  
        !          where X is the matrix of eigenvectors computed with the help of 
        !          ZFEAST_HCSREV. ZGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
        ! *******************************************************************************
        USE MKL_SPBLAS   
        USE ISO_C_BINDING 
        implicit none
        !!!!!!!!!!!!!!!!! External function declarations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        external zgemm, zfeast_hcsrev, zfeast_hcsrgv, feastinit
        !!!!!!!!!!!!!!!!! Matrix declaration variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer(kind=8), intent(in) :: n
        integer, intent(in) :: nnz
        integer, intent(in) :: nev0
        integer, intent(in) :: rc(nnz,2)
        double precision, intent(in) :: ham(nnz)
        double precision, intent(in) :: Emin,Emax
        integer, intent(inout) :: nest
        integer, intent(out) :: nev 
        double precision, allocatable, intent(out) :: eigvals(:)
        double precision, allocatable, intent(out) :: eigvecs(:,:) 
        
        TYPE(MATRIX_DESCR) :: descrA     ! Sparse matrix descriptor
        !   CSR matrix representation
        TYPE(SPARSE_MATRIX_T) :: cooA, csrA    ! Structure with sparse matrix

        integer :: rows(n+1), cols(nnz)
        integer :: rows_b(n+1), cols_b(n)
        INTEGER(C_INT) :: indexing
        double precision :: val(nnz), val_b(n)
        double precision :: beta
        character*1 :: UPLO
        parameter   (UPLO='F')
        !!!!!!!!!!!!!!!!! Feast declaration variable !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! E - eigenvalues, X - eigenvectors, res - residual !!!!!!!!!!!!
        integer :: fpm(128)  
        double precision :: epsout
        integer :: loop
        integer :: L
        parameter   (L=10)
        integer :: M0, M, info
        double precision :: E(nev0)
        !On output, the first m columns of x contain the orthonormal eigenvectors corresponding 
        !to the computed eigenvalues e, with the i-th column of x holding the eigenvector associated with e(i).   
        double precision :: X(n,nev0)
        double precision :: res(nev0)
        !!!!!!!!!!!!!!!!! Declaration of local variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! Eig - exact eigenvalues, R=|E-Eig|, Y=(X')*X-I !!!!!!!!!!!!!!!
        double precision :: Eig(nev0)
        double precision :: R(nev0)
        integer :: i,j,n4
        integer :: mkl, dim 
        integer :: UTnnz = 0 
        integer :: ldx, ldy
        double precision :: one, zero
        double precision :: smax, eigabs 

        integer, allocatable :: ib(:), jb(:), ic(:), jc(:)
        double precision, allocatable :: b(:), c(:)
        integer, allocatable :: UTcols(:), UTrows(:)
        double precision, allocatable :: UTval(:)

        print*,'Sparse matrix size',n
        print*,'Use Feast algorithm'

        M0=nev0
        M=10
        ! M=int(0.5*M0)
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

        n4 = n 

        call coocsr(n, nnz, ham, rc(1:nnz,1), rc(1:nnz,2), val, cols, rows)
        call coocsr(n, n4, b, ib, jb, val_b, cols_b, rows_b)

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

        ! call HermitianUpperTriangularCSR(n, nnz, val, rows, cols, UTval, UTrows, UTcols, UTnnz)
    
        !
        !        Task 1. Call  FEASTINIT  to define the default values for the input
        !        FEAST parameters.
        !
        call feastinit(fpm)
        fpm(1)  = 1 !Specifies whether Extended Eigensolver routines print runtime status. 
        fpm(26) = 1 !Specifies whether Extended Eigensolver routines check input matrices (applies to CSR format only). 
        fpm(27) = 1 !Specifies whether Extended Eigensolver routines check input matrices (applies to CSR format only). 

        print *, ' Testing dfeast_scsrev '
        !
        !         Task 2. Solve the standard eigenvalue problem Ax=ex.
        !
        
        call dfeast_scsrev(UPLO,N,val,rows,cols,fpm,epsout,loop,&
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
        !
        !         Task 3. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
        !         are the expected eigenvalues  and E(i) are eigenvalues computed
        !         with the help of ZFEAST_HCSREV().
        !
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

    end subroutine dfeast

    subroutine HermitianUpperTriangularCSR(n, nnz, A, IA, JA, UTriangular, IA_UTriangular, JA_UTriangular, nnz_UpperTriangular)
        implicit none
        integer(kind=8), intent(in) :: n
        integer, intent(in) :: nnz
        complex*16, intent(in) :: A(nnz)
        integer, intent(in) :: IA(n+1), JA(nnz)
        integer, intent(out) :: nnz_UpperTriangular
        complex*16, allocatable, intent(out) :: UTriangular(:)
        integer, allocatable, intent(out) :: IA_UTriangular(:), JA_UTriangular(:)
    
        ! Calculate the number of non-zero elements in the upper triangular matrix
        integer :: i, j 

        nnz_UpperTriangular = 0
        do i = 1, n
        do j = IA(i), IA(i + 1) - 1
            if (JA(j) >= i) nnz_UpperTriangular = nnz_UpperTriangular + 1
        end do
        end do
        
        ! Allocate memory for the upper triangular matrix
        allocate(UTriangular(nnz_UpperTriangular))
        allocate(IA_UTriangular(n + 1))
        allocate(JA_UTriangular(nnz_UpperTriangular))
    
        ! Set the first element of IA_UTriangular to 1
            IA_UTriangular = 0 
            IA_UTriangular(1) = 1

            ! ! Calculate the CSR format for the upper triangular matrix
            ! nnz_UpperTriangular = 0
            ! do i = 1, n
            !   IA_UTriangular(i) = nnz_UpperTriangular + 1
            !   do j = IA(i), IA(i + 1) - 1
            !     if (JA(j) >= i) then
            !       nnz_UpperTriangular = nnz_UpperTriangular + 1
            !       UTriangular(nnz_UpperTriangular) = A(j)
            !       JA_UTriangular(nnz_UpperTriangular) = JA(j)
            !     end if
            !   end do
            ! end do
            ! Calculate the CSR format for the upper triangular matrix
        nnz_UpperTriangular = 1
        IA_UTriangular(1) = nnz_UpperTriangular
        do i = 1, n
            do j = IA(i), IA(i + 1) - 1
            if (JA(j) >= i) then
                UTriangular(nnz_UpperTriangular) = A(j)
                JA_UTriangular(nnz_UpperTriangular) = JA(j)
                nnz_UpperTriangular = nnz_UpperTriangular + 1
            end if
            end do
            IA_UTriangular(i + 1) = nnz_UpperTriangular
        end do
        nnz_UpperTriangular = nnz_UpperTriangular - 1
        ! IA_UTriangular(n + 1) = nnz_UpperTriangular + 1

    end subroutine HermitianUpperTriangularCSR


    !----------------------------------!
    !            Utilities             !
    !----------------------------------!


    subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)

        !-----------------------------------------------------!
        !            COO to CSR for real sparse matrix        !
        !-----------------------------------------------------!
        implicit none

        integer(kind=8), intent(in) :: nrow
        integer, intent(in) :: nnz
        double precision, intent(in)  :: a(*)
        double precision, intent(out) :: ao(*)
        
        integer, intent(in) :: ir(*), jc(*)
        integer, intent(out) :: jao(*), iao(*)

        double precision, save :: x = 0 
        integer, save :: k = 0, i = 0, j = 0, k0 = 0, iad = 0

        !$omp threadprivate(x, k, i, j, k0, iad)

        ! c-----------------------------------------------------------------------
        ! c  Coordinate     to   Compressed Sparse Row 
        ! c----------------------------------------------------------------------- 
        ! c converts a matrix that is stored in coordinate format
        ! c  a, ir, jc into a row general sparse ao, jao, iao format.
        ! c
        ! c on entry:
        ! c--------- 
        ! c nrow  = dimension of the matrix 
        ! c nnz   = number of nonzero elements in matrix
        ! c a,
        ! c ir, 
        ! c jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
        ! c         nonzero elements of the matrix with a(k) = actual real value of
        ! c     the elements, ir(k) = its row number and jc(k) = its column 
        ! c     number. The order of the elements is arbitrary. 
        ! c
        ! c on return:
        ! c----------- 
        ! c ir    is destroyed
        ! c
        ! c ao, jao, iao = matrix in general sparse matrix format with ao 
        ! c   continung the real values, jao containing the column indices, 
        ! c   and iao being the pointer to the beginning of the row, 
        ! c   in arrays ao, jao.
        ! c
        ! c Notes:
        ! c------ This routine is NOT in place.  See coicsr
        ! c       On return the entries  of each row are NOT sorted by increasing 
        ! c       column number
        ! c
        ! c------------------------------------------------------------------------
        do k = 1, nrow + 1
            iao(k) = 0
        end do 
        ! c determine row-lengths.
        do k = 1, nnz
            iao(ir(k)) = iao(ir(k)) + 1
        end do 
        ! c starting position of each row..
        k = 1
        do j = 1, nrow + 1
            k0 = iao(j)
            iao(j) = k
            k = k + k0
        end do 
        ! c go through the structure  once more. Fill in output matrix.
        do k = 1, nnz
            i = ir(k)
            j = jc(k)
            x = a(k)
            iad = iao(i)
            ao(iad) =  x
            jao(iad) = j
            iao(i) = iad+1
        end do 
        ! c shift back iao
        do j = nrow, 1, -1
            iao(j+1) = iao(j)
        end do 
        iao(1) = 1    
        return
        ! c----¬--------- end of coocsr ------------------------------------------- 
        ! c----------------------------------------------------------------------- 
    end subroutine coocsr
            
    subroutine ccoocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
        implicit none

        !--------------------------------------------------------!
        !            COO to CSR for complex sparse matrix        !
        !--------------------------------------------------------!

        ! integer, intent(in) :: nrow, nnz
        integer(kind=8), intent(in) :: nrow
        integer, intent(in) :: nnz
        double complex, intent(in) :: a(*)
        integer, intent(in) :: ir(*),jc(*)    
        double complex, intent(out) :: ao(*)
        integer, intent(out) :: jao(*),iao(*)
        
        ! double complex, save :: x = 0 
        ! integer, save :: k = 0, k0 = 0, i = 0, j = 0, iad = 0
        double complex :: x = 0 
        integer :: k = 0, k0 = 0, i = 0, j = 0, iad = 0

        !$omp threadprivate(x, k, k0, i, j, iad)

        ! do 1 k = 1, nrow + 1
        do k = 1, nrow + 1
            iao(k) = 0
        end do 
        ! 1 continue

        ! do 2 k = 1, nnz
        do k = 1, nnz
            iao(ir(k)) = iao(ir(k)) + 1
        end do 
        ! 2 continue

        k = 1
        ! do 3 j = 1, nrow + 1
        do j = 1, nrow + 1
            k0 = iao(j)
            iao(j) = k
            k = k + k0
        end do 
        ! 3 continue

        ! do 4 k = 1, nnz
        do k = 1, nnz
            i = ir(k)
            j = jc(k)
            x = a(k)
            iad = iao(i)
            ao(iad) =  x
            jao(iad) = j
            iao(i) = iad + 1
        end do 
        ! 4 continue

        ! do 5 j = nrow, 1, -1
        do j = nrow, 1, -1
            iao(j+1) = iao(j)
        end do 
        ! 5 continue
        iao(1) = 1
        return

    end subroutine ccoocsr

    subroutine pamux(threads, n, x, y, a, ja, ia) !My parallelized sparse matrix-vector multiplication

        !--------------------------------------------!
        !            Parallel Spmmv: y = A*x         !
        !--------------------------------------------!


        implicit none

        integer(kind=8), intent(in) :: n
        integer, intent(in) :: threads
        integer, intent(in) :: ja(*), ia(*)
        double precision, intent(in)  :: x(n), a(*)

        double precision, intent(out) :: y(n)
        
        !-----------------------------------------------------------------------
        !         A times a vector
        !-----------------------------------------------------------------------
        ! multiplies a matrix by a vector using the dot product form
        ! Matrix A is stored in compressed sparse row storage.
        !
        ! on entry:
        !----------
        ! n     = row dimension of A
        ! x     = real array of length equal to the column dimension of
        !         the A matrix.
        ! a, ja,
        !    ia = input matrix in compressed sparse row format.
        !
        ! on return:
        !-----------
        ! y     = real array of length n, containing the product y=Ax
        !
        !-----------------------------------------------------------------------
        ! local variables
        !

        double precision, save :: t = 0 
        integer, save :: i = 0, k = 0
        integer :: cntr = 0  
        !-----------------------------------------------------------------------

        !$omp parallel do private(i, k, t) num_threads(threads)
        do i = 1, n

            !
            !     compute the inner product of row i with vector x
            !

            t = 0.0d0
            do k = ia(i), ia(i+1) - 1
                ! cntr = cntr + 1
                t = t + a(k) * x(ja(k))
            end do 

            !
            !     store result in y(i)
            !

            y(i) = t
        end do 
        !$omp end parallel do

        return
    end subroutine pamux 

    subroutine pamux2(threads, dim, sites, nbonds, nnnbonds, t1, v1, v2, basis, bsites, hexsites, x, y) !My parallelized sparse matrix-vector multiplication. 
        !Hamiltonian matrix entry is created on the fly 

        !--------------------------------------------!
        !            Parallel Spmmv: y = A*x         !
        !--------------------------------------------!


        implicit none

        integer(kind=8), intent(in) :: dim
        integer(kind=8), intent(in) :: basis(dim)
        integer, intent(in) :: threads, sites, nbonds, nnnbonds
        integer, intent(in) :: bsites(2, nbonds), hexsites(2, nnnbonds)
        double precision, intent(in)  :: t1, v1, v2
        double precision, intent(in)  :: x(dim)

        double precision, intent(out) :: y(dim)
        
        !-----------------------------------------------------------------------
        !         A times a vector
        !-----------------------------------------------------------------------
        ! multiplies a matrix by a vector using the dot product form
        ! Matrix A is stored in compressed sparse row storage.
        !
        ! on entry:
        !----------
        ! n     = row dimension of A
        ! x     = real array of length equal to the column dimension of
        !         the A matrix.
        ! a, ja,
        !    ia = input matrix in compressed sparse row format.
        !
        ! on return:
        !-----------
        ! y     = real array of length n, containing the product y=Ax
        !
        !-----------------------------------------------------------------------
        ! local variables
        !

        double precision :: a = 0.d0, t = 0.d0  
        integer(kind=8) :: i = 0, newst = 0, loc = 0
        integer :: s = 0, parity1 = 0, parity2 = 0
        ! integer :: id = 0

        !-----------------------------------------------------------------------
        
        
        !$omp parallel do private(i, loc, newst, a, t, s, parity1, parity2) num_threads(threads)
        do i = 1, dim
            loc = 0 
            newst = 0
            !
            !     compute the inner product of row i with vector x
            !
            a = 0.0d0
            t = 0.0d0
            !!$ id = omp_get_thread_num()
            
            do s = 1, nbonds  
                if (btest(basis(i), bsites(1, s) - 1) .and. .not. btest(basis(i), bsites(2, s) - 1)) then 
                    newst = ibclr( ibset( basis(i), bsites(2, s) - 1 ), bsites(1, s) - 1 ) !Create on site 2, annihilate on site 1 
                    call findstate(dim, newst, basis, loc) 
                    if(loc > 0) then 
                        parity1 = popcnt( ibits( basis(i), bsites(1, s), sites ) ) !Parity of site1
                        parity2 = popcnt( ibits( ibclr(basis(i), bsites(1, s) - 1), bsites(2, s), sites ) ) !Parity of site2
                        a = - t1 * (-1)**(parity1 + parity2) 
                        t = t + a * x(loc) 
                    end if 
                else if (btest(basis(i), bsites(2, s) - 1) .and. .not. btest(basis(i), bsites(1, s) - 1)) then 
                    newst = ibclr( ibset( basis(i), bsites(1, s) - 1 ), bsites(2, s) - 1 ) !Create on site 2, annihilate on site 1 
                    call findstate(dim, newst, basis, loc) 
                    if(loc > 0) then 
                        parity1 = popcnt( ibits( basis(i), bsites(2, s), sites ) ) !Parity of site2
                        parity2 = popcnt( ibits( ibclr(basis(i), bsites(2,s) - 1), bsites(1, s), sites ) ) !Parity of site1
                        a = - t1 * (-1)**(parity1 + parity2)  
                        t = t + a * x(loc) 
                    end if                    
                end if 
                if( btest(basis(i), bsites(1, s) - 1) .and. btest(basis(i), bsites(2, s) - 1) ) t = t + 0.25 * v1 * x(i)
                if( .not. (btest(basis(i), bsites(1, s) - 1)) .and. .not. (btest(basis(i), bsites(2, s) - 1)) ) t = t + 0.25 * v1 * x(i)
                if( btest(basis(i), bsites(1, s) - 1) .and. .not. (btest(basis(i), bsites(2, s) - 1)) ) t = t - 0.25 * v1 * x(i)
                if( .not. (btest(basis(i), bsites(1, s) - 1)) .and. btest(basis(i), bsites(2, s) - 1) ) t = t - 0.25 * v1 * x(i)
            end do !bond loop
            do s = 1, nnnbonds  
                if( btest(basis(i), hexsites(1, s) - 1) .and. btest(basis(i), hexsites(2, s) - 1) ) t = t + 0.25 * v2 * x(i)
                if( .not. (btest(basis(i), hexsites(1, s) - 1)) .and. .not. (btest(basis(i), hexsites(2, s) - 1)) ) t = t + 0.25 * v2 * x(i)
                if( btest(basis(i), hexsites(1, s) - 1) .and. .not. (btest(basis(i), hexsites(2, s) - 1)) ) t = t - 0.25 * v2 * x(i)
                if( .not. (btest(basis(i), hexsites(1, s) - 1)) .and. btest(basis(i), hexsites(2, s) - 1) ) t = t - 0.25 * v2 * x(i)
            end do !bond loop
            
            !
            !     store result in y(i)
            !

            y(i) = t
        end do 
        !$omp end parallel do

        return

    end subroutine pamux2 

    subroutine cpamux(threads, n, x, y, a, ja, ia) !My parallelized sparse matrix-vector multiplication
        implicit none

        !----------------------------------------------------!
        !            Parallel complex Spmmv: y = A*x         !
        !----------------------------------------------------!

        !include "omp_lib.h"
        ! complex*16, intent(in)  :: x(*), a(*)
        double complex, intent(in)  :: x(*), a(*)
        integer(kind=8), intent(in) :: n
        integer, intent(in) :: threads
        integer, intent(in) :: ja(*), ia(*)

        ! complex*16, intent(out) :: y(*)
        double complex, intent(out) :: y(*)
        
        
        !-----------------------------------------------------------------------
        !         A times a vector
        !-----------------------------------------------------------------------
        ! multiplies a matrix by a vector using the dot product form
        ! Matrix A is stored in compressed sparse row storage.
        !
        ! on entry:
        !----------
        ! n     = row dimension of A
        ! x     = real array of length equal to the column dimension of
        !         the A matrix.
        ! a, ja,
        !    ia = input matrix in compressed sparse row format.
        !
        ! on return:
        !-----------
        ! y     = real array of length n, containing the product y=Ax
        !
        !-----------------------------------------------------------------------
        ! local variables
        !
        ! complex*16 :: t
        double complex, save :: t
        integer, save :: i, k
        !-----------------------------------------------------------------------
        
        
            !$omp parallel do private(k,t,i) num_threads(threads)
            !   do 100 i = 1, n
            do i = 1, n
        !
        !     compute the inner product of row i with vector x
        !
        
                t = 0.0d0
                !  do 99 k = ia(i), ia(i+1)-1
                do k = ia(i), ia(i+1)-1
                    t = t + a(k)*x(ja(k))
        !  99      continue
                end do 
        !
        !     store result in y(i)
        !
                y(i) = t
        !  100  continue
                end do 
            !$omp end parallel do
        
        
            return
    end subroutine cpamux

end module diagonalization