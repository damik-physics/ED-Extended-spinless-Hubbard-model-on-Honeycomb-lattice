module io_routines 

    use functions 
    use file_utils
    
    implicit none 

    interface save
        module procedure save_momenta, save_spectrum, save_ham, save_currents, save_cdw
    end interface save

    interface params_string
        module procedure params_string, params_string_dd_cf
    end interface params_string

    contains 

    subroutine save_momenta(dir, cluster, unit, genCluster, sites, particles, bc, pattern, k1_max, k2_max)

        implicit none 
        
        integer,          intent(in), optional :: unit, genCluster, sites, particles, k1_max, k2_max
        character(len=*), intent(in), optional :: dir, cluster, bc, pattern

        integer   :: k1 = 0, k2 = 0 

        open(unit=unit, file=trim(outdir)//'momenta.dat', status='replace')
        do k1 = 0, k1_max
            do k2 = 0, k2_max
                write(unit,*) k1, k2 
            end do 
        end do
        close(unit)

        return 

    end subroutine save_momenta

    subroutine save_spectrum(type, par, append, conf, unit, dim, states, nev, nest, rvec, en, st_dp, st_c)
        
        implicit none 
        
        integer,           intent(in) :: conf, unit, states, nev, nest 
        integer(kind=8),   intent(in) :: dim  
        character(len=*),  intent(in) :: type, par
        logical,           intent(in) :: append, rvec 
        double precision,  intent(in) :: en(dim), st_dp(dim, nest) 
        double complex,    intent(in) :: st_c(dim, nest) 

        integer         :: j = 0
        character       :: file*256, appchar*20

        100 format(1000(F40.30))
        
        if(append) then 
            appchar = 'APPEND'
        else  
            appchar = 'SEQUENTIAL'
        end if 
        file = trim(outdir)//'energy_'//par//'.dat'
        open(unit=unit, file=trim(file), status='replace', access=appchar)

        ! file = dir // "energies_" // par
        ! file2 = dir // "states_" // par
        
        ! file = trim_name(file)
        ! file2 = trim_name(file2)

        ! open(unit, file = file, access = appchar)
        if(conf == 1) write(unit,*) nev

        do j = 1, nev
            write(unit, 100) en(j)
        end do
        close(unit)

        if(states == 1 .and. rvec) then
            file = trim(outdir)//'eigenstates_'//par//'.dat'
            if(type == "R") then 
                open(unit=unit, file=trim(file) , status='replace', access=appchar)
                ! open(unit, file = file2, access = appchar)
                write(unit,*) dim
                write(unit,*) nev
                do j = 1, nest
                    write(unit, 100) st_dp(1:dim, j)
                end do
            else if(type == "C") then
                open(unit=unit, file=trim(file) , status='replace', access=appchar) 
                ! open(unit, file = file2, access = appchar)
                write(unit,*) dim
                write(unit,*) nev
                do j = 1, nest               
                    write(unit, 100) dble(st_c(1:dim, j)), dimag(st_c(1:dim, j)) 
                end do
            end if
            close(unit)
        end if 


        return 

    end subroutine save_spectrum

    subroutine save_ham(type, unit, ti, sites, nnz, nDi, hamDi, ham_i, rc, ham_dp, ham_dc, occ)

        implicit none 

        integer,          intent(in)           :: ti, unit, sites  
        integer(kind=8),  intent(in)           :: nnz, nDi       
        integer(kind=8),  intent(in), optional :: hamDi(nDi, 2), ham_i(nnz, 3), rc(nnz, 2), occ(sites, nDi)
        double precision, intent(in), optional :: ham_dp(nnz)
        double complex,   intent(in), optional :: ham_dc(nnz) 
        character(len=*), intent(in)           :: type

        integer(kind=8) :: i = 0
        character       :: dir*1024, file*1024!, dir*256

        print*, 'Saving Hamiltonian to disk...'
        dir = trim(outdir) // "hamiltonians/"
        ! dir = trim_name(dir)
        file = trim(dir) // 'hamiltonian_off_diagonal.dat'
        ! file = trim_name(file)
        print*, 'Output directory: ', trim(outdir)
        print*, 'File: ', trim(file)
        open(unit, file=trim(file), status='replace')
        write(unit,*) nnz 
        
        if(type == "R" .and. ti == 0) then 
            do i = 1, nnz
                if(ti == 0) write(unit,*) ham_i(i,1), ham_i(i,2), ham_i(i,3)
            end do 
        else if(type == "R" .and. ti == 1) then         
            do i = 1, nnz
                write(unit,*) ham_dp(i), rc(i,1), rc(i,2)
            end do 
        else if(type == "C" .and. ti == 1) then             
            do i = 1, nnz
                write(unit,*) ham_dc(i), rc(i,1), rc(i,2)
            end do 
        end if 
        close(unit)
        
        file = trim(dir) // 'hamiltonian_diagonal' // '.dat'
        ! file = trim(dir // "diagonal_hamiltonian_" // params)
        open(unit, file=trim(file), status='replace')
        write(unit,*) nDi 
        do i = 1, nDi
            write(unit,*) hamDi(i,1), hamDi(i,2)
        end do 
        close(unit)

        file = trim(dir) // 'site_occupation' // '.dat'
        open(unit, file=trim(file), status='replace')
        write(unit,*) nDi 
        do i = 1, nDi
            write(unit,*) occ(1:sites,i)
        end do 
        close(unit)

        print*,'Hamiltonian and site-occupation saved to disk.'
        print*,''
        
        return 

    end subroutine save_ham

    subroutine loadham(type, unit, ti, sites, nnz, nDi, hamDi, ham_i, rc, ham_dp, ham_dc, occupation, exist1, exist2, exist3)
        
        implicit none 
        
        character(len=*), intent(in) :: type
        integer, intent(in) :: ti, unit, sites  
        integer(kind=8), intent(out) :: nnz, nDi       
        integer(kind=8), allocatable, intent(out), optional :: hamDi(:, :), ham_i(:, :), rc(:, :), occupation(:, :)
        double precision, allocatable, intent(out), optional :: ham_dp(:)
        double complex, allocatable, intent(out), optional :: ham_dc(:) 
        logical, intent(out) :: exist1, exist2, exist3

        integer(kind=8) :: i = 0
        character :: file*2024, dir*1024
        
        dir = trim("/hamiltonians/")
        file = trim(dir) // "hamiltonian_off_diagonal"
        inquire(file=trim(file), exist=exist1)

        if(exist1) then 
            open(unit=unit, file=trim(file))
            read(unit, *) nnz 

            if(type == "R" .and. ti == 0) then 
                if(allocated(ham_i)) deallocate(ham_i)
                allocate(ham_i(nnz, 3))
                ham_i = 0
                do i = 1, nnz
                    read(unit,*) ham_i(i,1), ham_i(i,2), ham_i(i,3)
                end do 
            else if(type == "R" .and. ti == 1) then 
                if(allocated(ham_dp)) deallocate(ham_dp)
                if(allocated(rc)) deallocate(rc)
                allocate(ham_dp(nnz))
                allocate(rc(nnz, 2))
                ham_dp = 0.d0
                rc = 0

                do i = 1, nnz
                    read(unit,*) ham_dp(i), rc(i,1), rc(i,2)
                end do 
            else if(type == "C" .and. ti == 1) then      
                if(allocated(ham_dc)) deallocate(ham_dc)
                if(allocated(rc)) deallocate(rc)
                allocate(ham_dc(nnz))
                allocate(rc(nnz, 2))       
                do i = 1, nnz
                    read(unit,*) ham_dc(i), rc(i,1), rc(i,2)
                end do 
            end if 
            close(unit)
            print*,'Off-diagonal Hamiltonian uploaded.'
            print*,''
        end if

        if(.not.(exist1)) print*,'Off-diagonal Hamiltonian does not yet exist.'

        file = trim(dir)//"hamiltonian_diagonal"
        inquire(file=trim(file), exist=exist2)

        if( exist2 ) then 
            open(unit, file=trim(file))
            read(unit, *) nDi 
            if(allocated(hamDi)) deallocate(hamDi)
            allocate(hamDi(nDi, 2))
            do i = 1, nDi
                read(unit,*) hamDi(i,1), hamDi(i,2)
            end do 
            close(unit)
            print*,'Diagonal Hamiltonian uploaded.'
            print*,''
        end if 
        
        if(.not. (exist2)) print*,'Diagonal Hamiltonian does not yet exist.'

        file = trim(dir)//"site_occupation"
        inquire(file=trim(file), exist=exist3)
        if( exist3 ) then 
            open(unit, file=trim(file))
            read(unit, *) nDi 
            if(allocated(occupation)) deallocate(occupation)
            allocate(occupation(sites, nDi))
            do i = 1, nDi
                read(unit,*) occupation(1:sites,i)
            end do 
            close(unit)
            print*,'Occupation uploaded.'
            print*,''
        end if 

        if(.not. (exist3)) print*,'Site occupation does not yet exist.'
        if(.not. (exist3)) print*,''

        return 

    end subroutine loadham

    subroutine save_currents(dir, sl, params, append, degeneracy, unit, nbonds, current, bondcurrent)
        
        implicit none 

        integer,          intent(in) :: degeneracy, unit, nbonds
        double precision, intent(in) :: current, bondcurrent(nbonds)
        character(len=*), intent(in) :: dir, sl, params
        logical,          intent(in) :: append 

        character :: file*512, file2*512, appchar*20, dirq*512

        100 format(1000(F30.20))

        if(append) then 
            appchar = 'APPEND'
        else  
            appchar = 'SEQUENTIAL'
        end if 

        if(degeneracy < 2) dirq = dir  
        if(degeneracy == 2) dirq = trim_name(dir // "QD_")
        dirq = trim_name(dirq)
        
        if(sl == "A") then 
            file = trim_name(dirq // "current_A_" // params)
            file = trim_name(file)
            file2 = trim_name(dirq // "bondcurrent_A_" // params)
            file2 = trim_name(file2)
        else if(sl == "B") then 
            file = dirq // "current_B_" // params
            file = trim_name(file)
            file2 = dirq // "bondcurrent_B_" // params
            file2 = trim_name(file2)
        end if 


        open(unit, file = file, access = appchar)
        write(unit, 100) current
        close(unit)
        open(unit, file = file2, access = appchar)
        write(unit, 100) bondcurrent
        close(unit)

        return 

    end subroutine save_currents

    subroutine save_cdw(dir, params, append, unit, ucx, ucy, nbonds, rho, rhotot)
        
        implicit none 

        integer,          intent(in) :: unit, ucx, ucy, nbonds
        double precision, intent(in) :: rho(ucx*ucy, 3), rhotot(nbonds, 3) 
        character(len=*), intent(in) :: dir, params
        logical,          intent(in) :: append 

        integer   :: i = 0
        character :: file*256, file2*256, appchar*20

        100 format(1000(F30.20))
        
        if(append) then 
            appchar = 'APPEND'
        else  
            appchar = 'SEQUENTIAL'
        end if 

        file = dir // "cdw_" // params
        file = trim_name(file)
        file2 = dir // "cdw_total_" // params
        file2 = trim_name(file2)
        
        open(unit, file = file, access = appchar)
        do i = 1, ucx*ucy
            write(unit,100) rho(i,1), rho(i,2), rho(i,3)    
        end do 
        close(unit)
        open(unit, file = file2, access = appchar)
        do i = 1, nbonds
            write(unit,100) rhotot(i,1), rhotot(i,2), rhotot(i,3)    
        end do 
        close(unit)

        return 

    end subroutine save_cdw

    subroutine save_dpd_cf(dir, params, append, unit, sites, rhocn, rho)
        ! Calculates the density-density correlation function (dd_cf) as well as density correlation function and saves them to disk.
        implicit none

        integer,          intent(in) :: unit, sites
        double precision, intent(in) :: rhocn(sites), rho(sites) 
        character(len=*), intent(in) :: dir, params
        logical,          intent(in) :: append 

        integer   :: i = 0
        character :: file*400, file2*400, appchar*20

        100 format(1000(F30.20))
        
        if(append) then 
            appchar = 'APPEND'
        else  
            appchar = 'SEQUENTIAL'
        end if 

        file = dir // "dd_cf_" // params 
        file = trim_name(file)
        file2 = dir // "density_cf_" // params
        file2 = trim_name(file2)

        open(unit, file = file, access = appchar)
        do i = 1, sites
            write(unit,100) rhocn(i)
        end do 
        close(unit)
        open(unit, file = file2, access = appchar)
        do i = 1, sites
            write(unit,100) rho(i)
        end do 
        close(unit)

        return 

    end subroutine save_dpd_cf

    subroutine save_i(sa, unit, dir, name, pars, rows, cols, scalar, arr)
        
        implicit none 
        
        integer,          intent(in)           :: unit, scalar
        integer(kind=8),  intent(in)           :: rows, cols
        integer,          intent(in), optional :: arr(rows, cols)
        character(len=*), intent(in)           :: sa, pars, name, dir

        integer(kind=8) :: j = 0
        character       :: file*512

        j = 0 
        
        file = trim_name(dir // name // pars)
        open(unit, file = file)

        if(sa == "S") then 
            write(unit,*) scalar
        end if 
        if(sa == "A") then  
            do j = 1, cols
                write(unit,*) arr(1:rows, j)
            end do
        end if 
        close(unit)

        return 

    end subroutine save_i

    subroutine save_i8(sa, unit, dir, name, pars, rows, cols, scalar, arr)
        
        implicit none 
        
        integer,          intent(in)           :: unit
        integer(kind=8),  intent(in)           :: scalar, rows, cols  
        integer(kind=8),  intent(in), optional :: arr(rows, cols) 
        character(len=*), intent(in)           :: sa, pars, name, dir

        integer(kind=8) :: j = 0
        character       :: file*512

        j = 0 
        file = trim_name(dir // name // pars)
        
        open(unit, file = file)

        if(sa == "S") then 
            write(unit,*) scalar
        end if 
        if(sa == "A") then 
            do j = 1, cols
                write(unit,*) arr(1:rows, j)
            end do
        end if 
        close(unit)

        return 

    end subroutine save_i8

    subroutine save_dp(sa, unit, dir, name, pars, rows, cols, scalar, arr, transp)
        
        implicit none 
        
        integer,          intent(in)           :: unit
        integer,          intent(in), optional :: transp
        integer(kind=8),  intent(in)           :: rows, cols  
        double precision, intent(in)           :: scalar
        double precision, intent(in), optional :: arr(rows, cols) 
        character(len=*), intent(in)           :: sa, pars, name, dir

        integer(kind=8) :: j = 0
        character       :: file*512

        j = 0 
        99 format(1000(F40.30))
        
        file = trim_name(dir // name // pars)
        open(unit, file = file)
        
        if(sa == "S") then 
            write(unit,99) scalar
        end if 
        if(sa == "A") then 
            if(rows < 1000) then 
                if(present(transp)) then 
                    do j = 1, rows
                        write(unit,100) arr(j, 1:cols)
                    end do
                else
                    do j = 1, cols
                        write(unit,100) arr(1:rows, j)
                    end do
                end if 
            else 
                100 format(1000000000(F40.30))
                do j = 1, cols
                    write(unit,100) arr(1:rows, j)
                end do
            end if
        end if 
        close(unit)

        return 

    end subroutine save_dp

    subroutine save_dc(sa, unit, dir, name, pars, rows, cols, scalar, arr, transp)
        
        implicit none 
        
        integer,           intent(in)           :: unit
        integer(kind=8),   intent(in)           :: rows, cols  
        integer,           intent(in), optional :: transp
        double complex,    intent(in)           :: scalar
        double complex,    intent(in), optional :: arr(rows, cols) 
        character(len=*),  intent(in)           :: sa, pars, name, dir

        integer(kind=8) :: j = 0
        character       :: file*512

        j = 0 
        100 format(1000(F40.30))
        
        file = trim_name(dir // name // pars)
    
        open(unit, file = file)

        if(sa == "S") then 
            write(unit,100) scalar
        else if(sa == "A") then 
            if(present(transp)) then 
                do j = 1, rows
                    write(unit,100) arr(j, 1:cols)
                end do
            else
                do j = 1, cols
                    write(unit,100) arr(1:rows, j)
                end do
            end if 
        end if 
        close(unit)

        return 

    end subroutine save_dc

    subroutine params_string(sites, particles, k1, k2, refb, v1, v2, pattern, params, type)

        ! Generates a string for the params of the Hamiltonian or energy and eigenstate files.
        ! The string is used to create the file names for the respective files.
        use input_variables

        implicit none 

        integer,          intent(in) :: sites, particles, k1, k2, refb
        double precision, intent(in) :: v1, v2 
        character(len=*), intent(in) :: pattern 

        character,        intent(out)           :: params*512 
        character,        intent(out), optional :: type*1 

        character :: string*256

        if(tilted == 1) then 
            write(string, "(a,'_')") cluster 
            string = trim_name(string)
        end if 
        params = ""    
        if(refb < 0) then !Name for Hamiltonian file
            type = 'R'
            if(((k1 .ne. 0) .or. (k2 .ne. 0))) type = 'C'
            if(ti == 0) write (params, "('L=',i0,'N=',i0,'BC=',a,'_pat=',a2'.dat')") sites, particles, bc, pattern
            if(ti == 1) write (params, "('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'BC=',a,'_pat=',a2'.dat')") sites, particles, k1, k2, bc, pattern
            if(symmetrize == 1) write (params, "('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'_',a,'_BC=',a,'_pat=',a2'.dat')") sites, particles, k1, k2, irrep, bc, pattern
        else if(refb == 0) then !Name for energy/states file
            if(ti == 0) write (params, "('L=',i0,'N=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2'.dat')") sites, particles, v1, v2, mass, dis, nDis, bc, pattern
            if(ti == 1) write (params, "('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2'.dat')") sites, particles, k1, k2, v1, v2, mass, dis, nDis, bc, pattern
            if(symmetrize == 1) write (params, "('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'_',a,'_V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2'.dat')") sites, particles, k1, k2, irrep, v1, v2, mass, dis, nDis, bc, pattern
        else if(refb > 0) then !Name for current files
            if(ti == 0) write (params, "('L=',i0,'N=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refb=',i0,'.dat')") sites, particles, v1, v2, mass, dis, nDis, bc, pattern, refb            
            if(ti == 1) write (params, "('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refb=',i0,'.dat')") sites, particles, k1, k2, v1, v2, mass, dis, nDis, bc, pattern, refb 
            if(symmetrize == 1) write (params, "('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'_',a,'_V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refb=',i0,'.dat')") sites, particles, k1, k2, irrep, v1, v2, mass, dis, nDis, bc, pattern, refb
        end if 
        params = trim_name(params)
        if(tilted == 1) params = trim_name(string // params)
        params = trim_name(params)    

    end subroutine params_string

    subroutine params_string_dd_cf(sites, particles, k1, k2, v1, v2, refsite, pattern, params)
        
        ! Generates a string for the params of the density-density correlation function (dd_cf) for the output file-name.
        use input_variables

        implicit none 

        integer,          intent(in)  :: sites, particles, k1, k2, refsite
        double precision, intent(in)  :: v1, v2 
        character(len=*), intent(in)  :: pattern 
 
        character,        intent(out) :: params*400 

        character :: string*256

        if(tilted == 1) then 
            write(string, "(a,'_')") cluster 
            string = trim_name(string)
        end if 
        params = ""    
        if(refsite > 0) then !Name for current files
            if(ti == 0) write (params, "('L=',i0,'N=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refsite=',i0,'.dat')") sites, particles, v1, v2, mass, dis, nDis, bc, pattern, refsite            
            if(ti == 1) write (params, "('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refsite=',i0,'.dat')") sites, particles, k1, k2, v1, v2, mass, dis, nDis, bc, pattern, refsite 
            if(symmetrize == 1) write (params, "('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,a,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refsite=',i0,'.dat')") sites, particles, k1, k2, irrep, v1, v2, mass, dis, nDis, bc, pattern, refsite
        end if 
        params = trim_name(params)
        if(tilted == 1) then 
            params = trim_name(string // params)
        end if 
        params = trim_name(params)    

    end subroutine params_string_dd_cf

    subroutine load_ham(type, params, unit, ti, sites, nnz, nDi, hamDi, ham_i, rc, ham_dp, ham_dc, occupation, exist1, exist2, exist3)
        ! Loads the Hamiltonian from disk.
        ! The Hamiltonian is stored in two files: one for the off-diagonal part and
        ! one for the diagonal part. The off-diagonal part is stored in a file named
        ! "hopping_hamiltonian_<params>" and the diagonal part is stored in a
        ! file named "diagonal_hamiltonian_<params>". The occupation is stored in a
        ! file named "occupation_<params>". 
        implicit none 
        
        integer,          intent(in) :: ti, unit, sites  
        character(len=*), intent(in) :: type, params

        integer,                       intent(out)           :: nnz, nDi       
        logical,                       intent(out)           :: exist1, exist2, exist3
        integer,          allocatable, intent(out), optional :: hamDi(:, :), ham_i(:, :), rc(:, :), occupation(:, :)
        double precision, allocatable, intent(out), optional :: ham_dp(:)
        double complex,   allocatable, intent(out), optional :: ham_dc(:) 

        integer   :: i = 0
        character :: file*256, dir*256
        
        dir = "hamiltonians/"
        dir = trim_name(dir)

        file = dir // "hopping_hamiltonian_" // params
        file = trim_name(file)
        inquire(file = file, exist = exist1)

        if( exist1) then 
            open(unit, file = file)
            read(unit, *) nnz 

            if(type == "R" .and. ti == 0) then 
                if(allocated(ham_i)) deallocate(ham_i)
                allocate(ham_i(nnz, 3))
                ham_i = 0
                do i = 1, nnz
                    read(unit,*) ham_i(i,1), ham_i(i,2), ham_i(i,3)
                end do 
            else if(type == "R" .and. ti == 1) then 
                if(allocated(ham_dp)) deallocate(ham_dp)
                if(allocated(rc)) deallocate(rc)
                allocate(ham_dp(nnz))
                allocate(rc(nnz, 2))
                ham_dp = 0.d0
                rc = 0

                do i = 1, nnz
                    read(unit,*) ham_dp(i), rc(i,1), rc(i,2)
                end do 
            else if(type == "C" .and. ti == 1) then      
                if(allocated(ham_dc)) deallocate(ham_dc)
                if(allocated(rc)) deallocate(rc)
                allocate(ham_dc(nnz))
                allocate(rc(nnz, 2))       
                do i = 1, nnz
                    read(unit,*) ham_dc(i), rc(i,1), rc(i,2)
                end do 
            end if 
            close(unit)
            print*,'Off-diagonal Hamiltonian uploaded.'
            print*,''
        end if

        if(.not. (exist1)) print*,'Off-diagonal Hamiltonian does not yet exist.'

        file = ""
        file = dir // "diagonal_hamiltonian_" // params
        file = trim_name(file)
        inquire(file = file, exist = exist2)

        if( exist2) then 
            open(unit, file = file)
            read(unit, *) nDi 
            if(allocated(hamDi)) deallocate(hamDi)
            allocate(hamDi(nDi, 2))
            do i = 1, nDi
                read(unit,*) hamDi(i,1), hamDi(i,2)
            end do 
            close(unit)
            print*,'Diagonal Hamiltonian uploaded.'
            print*,''
        end if 
        
        if(.not. (exist2)) print*,'Diagonal Hamiltonian does not yet exist.'

        file = ""
        file = "hamiltonians/occupation_" // params
        file = trim_name(file)
        inquire(file = file, exist = exist3)
        if( exist3) then 
            open(unit, file = file)
            read(unit, *) nDi 
            if(allocated(occupation)) deallocate(occupation)
            allocate(occupation(sites, nDi))
            do i = 1, nDi
                read(unit,*) occupation(1:sites,i)
            end do 
            close(unit)
            print*,'Occupation uploaded.'
            print*,''
        end if 

        if(.not. (exist3)) print*,'Occupation does not yet exist.'
        if(.not. (exist3)) print*,''

        return 

    end subroutine load_ham


end module io_routines