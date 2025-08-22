module io_utils 
    ! Module for input/output utilities including file I/O operations, parameter string generation,
    ! spectrum saving, Hamiltonian I/O, and configuration printing functionality
    use types
    use functions 
    use file_utils
    use params 
        
    implicit none 

    interface save
        module procedure save_momenta, save_spectrum_dp, save_spectrum_dc, save_ham, save_currents, save_cdw
    end interface save
    interface params_string
        module procedure params_string, params_string_dd_cf
    end interface params_string
    interface printing
        module procedure print_configuration, &
                         print_parameters
    end interface printing 

    contains 

    subroutine save_momenta(dir, cluster, unit, genCluster, sites, particles, bc, pattern, k1_max, k2_max)
        ! Saves discrete momentum values (k1, k2) to file for momentum space calculations
        
        implicit none 
        
        integer,          intent(in), optional :: unit, genCluster, sites, particles, k1_max, k2_max
        character(len=*), intent(in), optional :: dir, cluster, bc, pattern

        integer                                :: k1, k2

        open(unit=unit, file=trim(dir)//'momenta.dat', status='replace')
        do k1 = 0, k1_max
            do k2 = 0, k2_max
                write(unit,*) k1, k2 
            end do 
        end do
        close(unit)

        return 

    end subroutine save_momenta

    subroutine save_spectrum_dp(dir, par, append, conf, unit, dim, states, nev, nest, rvec, en, st)
        ! Saves the spectrum of eigenvalues and real eigenstates to disk.
        implicit none 
        character(len=*), intent(in) :: dir
        integer,          intent(in) :: conf, unit, states, nev, nest 
        integer(kind=8),  intent(in) :: dim  
        character(len=*), intent(in) :: par
        logical,          intent(in) :: append, rvec 
        double precision, intent(in) :: en(dim), st(dim, nest) 

        integer                      :: j
        character                    :: file*256, appchar*20

        100 format(1000(F40.30))
        
        if(append) then 
            appchar = 'APPEND'
        else  
            appchar = 'SEQUENTIAL'
        end if 
        file = trim(dir)//'energy_'//par//'.dat'
        open(unit=unit, file=trim(file), status='replace', access=appchar)

        if(conf == 1) write(unit,*) nev
        do j = 1, nev
            write(unit, 100) en(j)
        end do
        close(unit)

        if(states == 1 .and. rvec) then
            file = trim(dir)//'eigenstates_'//par//'.dat' 
            open(unit=unit, file=trim(file) , status='replace', access=appchar)
            write(unit,*) dim
            write(unit,*) nev
            do j = 1, nest
                write(unit, 100) st(1:dim, j)
            end do
            close(unit)
        end if 


        return 

    end subroutine save_spectrum_dp

    subroutine save_spectrum_dc(dir, par, append, conf, unit, dim, states, nev, nest, rvec, en, st)
        ! Saves the spectrum of eigenvalues and complex eigenstates to disk.
        implicit none 
        character(len=*), intent(in) :: dir
        integer,          intent(in) :: conf, unit, states, nev, nest 
        integer(kind=8),  intent(in) :: dim  
        character(len=*), intent(in) :: par
        logical,          intent(in) :: append, rvec 
        double precision, intent(in) :: en(dim) 
        double complex,   intent(in) :: st(dim, nest) 

        integer                      :: j
        character                    :: file*256, appchar*20

        100 format(1000(F40.30))
        
        if(append) then 
            appchar = 'APPEND'
        else  
            appchar = 'SEQUENTIAL'
        end if 
        file = trim(dir)//'energy_'//par//'.dat'
        open(unit=unit, file=trim(file), status='replace', access=appchar)

        if(conf == 1) write(unit,*) nev
        do j = 1, nev
            write(unit, 100) en(j)
        end do
        close(unit)

        if(states == 1 .and. rvec) then
            file = trim(dir)//'eigenstates_'//par//'.dat'
            open(unit=unit, file=trim(file) , status='replace', access=appchar) 
            write(unit,*) dim
            write(unit,*) nev
            do j = 1, nest               
                write(unit, 100) dble(st(1:dim, j)), dimag(st(1:dim, j)) 
            end do
            close(unit)
        end if 


        return 

    end subroutine save_spectrum_dc

    subroutine save_ham(outdir, mat_type, unit, ti, sites, nnz, nDi, hamDi, ham_i, rc, ham_dp, ham_dc, occ)
        ! Saves Hamiltonian matrix elements and site occupations to disk files in .dat format. To be replaced by I/O storage to binary HDF5 files in the future.

        implicit none

        integer,                       intent(in)           :: ti, unit, sites  
        integer(kind=8),               intent(in)           :: nnz, nDi       
        integer(kind=8),               intent(in), optional :: hamDi(nDi, 2), ham_i(nnz, 3), rc(nnz, 2), occ(sites, nDi)
        double precision,              intent(in), optional :: ham_dp(nnz)
        double complex,                intent(in), optional :: ham_dc(nnz) 
        character(len=:), allocatable, intent(in)           :: outdir, mat_type

        integer(kind=8)                                     :: i
        character(len=:), allocatable                       :: dir, file

        print*, 'Saving Hamiltonian to disk...'
        dir  = trim(outdir) // "hamiltonians/"
        file = trim(dir) // 'hamiltonian_off_diagonal.dat'
        print*, 'Output directory: ', trim(dir)
        print*, 'File: ', trim(file)
        open(unit, file=trim(file), status='replace')
        write(unit,*) nnz 
        
        if(mat_type == "R" .and. ti == 0) then 
            do i = 1, nnz
                if(ti == 0) write(unit,*) ham_i(i,1), ham_i(i,2), ham_i(i,3)
            end do 
        else if(mat_type == "R" .and. ti == 1) then         
            do i = 1, nnz
                write(unit,*) ham_dp(i), rc(i,1), rc(i,2)
            end do 
        else if(mat_type == "C" .and. ti == 1) then             
            do i = 1, nnz
                write(unit,*) ham_dc(i), rc(i,1), rc(i,2)
            end do 
        end if 
        close(unit)
        
        file = trim(dir) // 'hamiltonian_diagonal' // '.dat'
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

    subroutine loadham(mat_type, dir, unit, ti, sites, nnz, nDi, hamDi, ham_i, rc, ham_dp, ham_dc, occupation, exist1, exist2, exist3)
        ! Loads previously saved Hamiltonian matrix elements from disk files in dat format. To be replaced by I/O storage to binary HDF5 files in the future.
        
        implicit none 
        
        character(len=1),              intent(in)            :: mat_type
        character(len=:), allocatable, intent(in)            :: dir
        integer,                       intent(in)            :: ti, unit, sites  
        integer(kind=8),               intent(out)           :: nnz, nDi       
        integer(kind=8), allocatable,  intent(out), optional :: hamDi(:, :), ham_i(:, :), rc(:, :), occupation(:, :)
        double precision, allocatable, intent(out), optional :: ham_dp(:)
        double complex, allocatable,   intent(out), optional :: ham_dc(:) 
        logical,                       intent(out)           :: exist1, exist2, exist3

        integer(kind=8)                                      :: i
        character(len=:),    allocatable                     :: file, base

        base = trim("hamiltonians/")
        file = trim(dir) // trim(base) // "hamiltonian_off_diagonal" // '.dat'
        inquire(file=trim(file), exist=exist1)

        if(exist1) then 
            open(unit=unit, file=trim(file))
            read(unit, *) nnz 

            if(mat_type == "R" .and. ti == 0) then 
                if(allocated(ham_i)) deallocate(ham_i)
                allocate(ham_i(nnz, 3))
                ham_i = 0
                do i = 1, nnz
                    read(unit,*) ham_i(i,1), ham_i(i,2), ham_i(i,3)
                end do 
            else if(mat_type == "R" .and. ti == 1) then 
                if(allocated(ham_dp)) deallocate(ham_dp)
                if(allocated(rc)) deallocate(rc)
                allocate(ham_dp(nnz))
                allocate(rc(nnz, 2))
                ham_dp = 0.d0
                rc = 0

                do i = 1, nnz
                    read(unit,*) ham_dp(i), rc(i,1), rc(i,2)
                end do 
            else if(mat_type == "C" .and. ti == 1) then      
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

        file = trim(dir) // trim(base) // "hamiltonian_diagonal" // ".dat"
        inquire(file=trim(file), exist=exist2)

        if(exist2) then 
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
        
        if(.not.(exist2)) print*,'Diagonal Hamiltonian does not yet exist.'

        file = trim(dir) // trim(base) // "site_occupation" // ".dat"
        inquire(file=trim(file), exist=exist3)
        if(exist3) then 
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

        if(.not.(exist3)) print*,'Occupation does not yet exist.'
        if(.not.(exist3)) print*,''

        return 

    end subroutine loadham

    subroutine save_currents(dir, sl, append, unit, nbonds, current, bondcurrent)
        ! Saves current and bond current data to .dat files with appropriate sublattice labeling. Obsolete due to I/O routines to csv file format.
        
        implicit none 

        integer,                       intent(in) :: unit, nbonds
        double precision,              intent(in) :: current, bondcurrent(nbonds)
        character(len=:), allocatable, intent(in) :: dir, sl
        logical,                       intent(in) :: append

        character(len=:), allocatable :: file, file2, appchar, base

        100 format(1000(F30.20))
        base = "correlations/"
        if(append) then 
            appchar = 'APPEND'
        else  
            appchar = 'SEQUENTIAL'
        end if 

        
        if(sl == "A") then 
            file  = trim(dir) // trim(base) // "current_A_" // ".dat"
            file2 = trim(dir) // trim(base) // "bondcurrent_A_" // ".dat"
        else if(sl == "B") then 
            file  = trim(dir) // trim(base) // "current_B_" // ".dat"
            file2 = trim(dir) // trim(base) // "bondcurrent_B_" // ".dat"
        end if 

        open(unit, file = file, access = appchar)
        write(unit, 100) current
        close(unit)
        open(unit, file = file2, access = appchar)
        write(unit, 100) bondcurrent
        close(unit)

        return 

    end subroutine save_currents

    subroutine save_cdw(dir, append, unit, ucx, ucy, nbonds, rho, rhotot)
        ! Saves charge density wave (CDW) data for individual sites and total density to .dat files on disk. Obsolete due to I/O routines to csv file format.

        implicit none 

        integer,                       intent(in) :: unit, ucx, ucy, nbonds
        double precision,              intent(in) :: rho(ucx*ucy, 3), rhotot(nbonds, 3) 
        character(len=:), allocatable, intent(in) :: dir
        logical,                       intent(in) :: append 

        integer                                   :: i
        character(len=:), allocatable             :: file, file2, appchar, base

        100 format(1000(F30.20))
        base = "correlations/"
        
        if(append) then 
            appchar = 'APPEND'
        else  
            appchar = 'SEQUENTIAL'
        end if 

        file = trim(dir) // trim(base) // "cdw_" // ".dat"
        file2 = trim(dir) // trim(base) // "cdw_total_" // ".dat"

        open(unit, file = trim(file), access = appchar)
        do i = 1, ucx*ucy
            write(unit,100) rho(i,1), rho(i,2), rho(i,3)    
        end do 
        close(unit)
        open(unit, file = trim(file2), access = appchar)
        do i = 1, nbonds
            write(unit,100) rhotot(i,1), rhotot(i,2), rhotot(i,3)    
        end do 
        close(unit)

        return 

    end subroutine save_cdw

    subroutine save_dd_cf(dir, append, unit, sites, rhocn, rho)
        ! Calculates the density-density correlation function (dd_cf) as well as density correlation function and saves them to disk as .dat files. Obsolete due to I/O routines to csv file format.
        implicit none

        integer,                       intent(in) :: unit, sites
        double precision,              intent(in) :: rhocn(sites), rho(sites) 
        character(len=:), allocatable, intent(in) :: dir
        logical,                       intent(in) :: append 

        integer                                   :: i
        character(len=:), allocatable             :: file, file2, appchar, base

        100 format(1000(F30.20))
        base = "correlations/"
        if(append) then 
            appchar = 'APPEND'
        else  
            appchar = 'SEQUENTIAL'
        end if 

        file  = dir // base // "dd_cf_" // ".dat" 
        file2 = dir // base // "density_cf_" // ".dat"

        open(unit, file = trim(file), access = appchar)
        do i = 1, sites
            write(unit,100) rhocn(i)
        end do 
        close(unit)
        open(unit, file = trim(file2), access = appchar)
        do i = 1, sites
            write(unit,100) rho(i)
        end do 
        close(unit)

        return 

    end subroutine save_dd_cf

    subroutine save_i(sa, unit, dir, name, pars, rows, cols, scalar, arr)
        ! Generic subroutine to save integer scalars or arrays to .dat files. Obsolete due to I/O routines to csv file format.
        
        implicit none 
        
        integer,                       intent(in)           :: unit, scalar
        integer(kind=8),               intent(in)           :: rows, cols
        integer,                       intent(in), optional :: arr(rows, cols)
        character(len=*),              intent(in)           :: sa, pars, name, dir

        integer(kind=8)                                     :: j
        character                                           :: file*512

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
        ! Generic subroutine to save 8-byte integer scalars or arrays to .dat files. Obsolete due to I/O routines to csv file format.
        
        implicit none 
        
        integer,                       intent(in)           :: unit
        integer(kind=8),               intent(in)           :: scalar, rows, cols  
        integer(kind=8),               intent(in), optional :: arr(rows, cols) 
        character(len=*),              intent(in)           :: sa, pars, name, dir

        integer(kind=8)                                     :: j
        character                                           :: file*512

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
        ! Generic subroutine to save double precision scalars or arrays to .dat files. Obsolete due to I/O routines to csv file format.
        implicit none 
        
        integer,                       intent(in)           :: unit
        integer,                       intent(in), optional :: transp
        integer(kind=8),               intent(in)           :: rows, cols  
        double precision,              intent(in)           :: scalar
        double precision,              intent(in), optional :: arr(rows, cols) 
        character(len=*),              intent(in)           :: sa, pars, name, dir

        integer(kind=8)                                     :: j
        character                                           :: file*512

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
        ! Generic subroutine to save double complex scalars or arrays to .dat files. Obsolete due to I/O routines to csv file format.
        implicit none 
        
        integer,                       intent(in)           :: unit
        integer(kind=8),               intent(in)           :: rows, cols  
        integer,                       intent(in), optional :: transp
        double complex,                intent(in)           :: scalar
        double complex,                intent(in), optional :: arr(rows, cols) 
        character(len=*),              intent(in)           :: sa, pars, name, dir

        integer(kind=8)                                     :: j
        character                                           :: file*512

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

    subroutine params_string(tilted, cluster, bc, ti, symm, irrep, sites, particles, k1, k2, refb, ndis, v1, v2, mass, dis, pattern, params, mat_type)
        ! Generates a string for the parameters of the Hamiltonian or energy and eigenstate files.
        ! The string is used to create the file names for the respective files.

        implicit none 
    
        integer,          intent(in)  :: tilted, ti, symm, sites, particles, k1, k2, refb, ndis
        double precision, intent(in)  :: v1, v2, mass, dis 
        character(len=*), intent(in)  :: pattern, cluster, bc, irrep 
        character,        intent(out) :: params*512, mat_type*1 

        character                     :: string*256

        if(tilted == 1) then 
            write(string, "(a,'_')") cluster 
            string = trim_name(string)
        end if 
        params = ""    
        if(refb < 0) then !Name for Hamiltonian file
            mat_type = 'R'
            if(((k1 .ne. 0) .or.(k2 .ne. 0))) mat_type = 'C'
            if(ti == 0) write(params, "('L=',i0,'N=',i0,'BC=',a,'_pat=',a2'.dat')") sites, particles, bc, pattern
            if(ti == 1) write(params, "('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'BC=',a,'_pat=',a2'.dat')") sites, particles, k1, k2, bc, pattern
            if(symm == 1) write(params, "('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'_',a,'_BC=',a,'_pat=',a2'.dat')") sites, particles, k1, k2, irrep, bc, pattern
        else if(refb == 0) then !Name for energy/states file
            if(ti == 0) write(params, "('L=',i0,'N=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2'.dat')") sites, particles, v1, v2, mass, dis, nDis, bc, pattern
            if(ti == 1) write(params, "('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2'.dat')") sites, particles, k1, k2, v1, v2, mass, dis, nDis, bc, pattern
            if(symm == 1) write(params, "('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'_',a,'_V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2'.dat')") sites, particles, k1, k2, irrep, v1, v2, mass, dis, nDis, bc, pattern
        else if(refb > 0) then !Name for current files
            if(ti == 0) write(params, "('L=',i0,'N=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refb=',i0,'.dat')") sites, particles, v1, v2, mass, dis, nDis, bc, pattern, refb            
            if(ti == 1) write(params, "('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refb=',i0,'.dat')") sites, particles, k1, k2, v1, v2, mass, dis, nDis, bc, pattern, refb 
            if(symm == 1) write(params, "('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'_',a,'_V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refb=',i0,'.dat')") sites, particles, k1, k2, irrep, v1, v2, mass, dis, nDis, bc, pattern, refb
        end if 
        params = trim_name(params)
        if(tilted == 1) params = trim_name(string // params)
        params = trim_name(params)    

    end subroutine params_string

    subroutine params_string_dd_cf(tilted, cluster, bc, ti, symm, irrep, sites, particles, k1, k2, ndis, refsite, v1, v2, mass, dis, pattern, params)
        ! Generates a string for the params of the density-density correlation function(dd_cf) for the output file-name.

        implicit none 
        
        integer,          intent(in)  :: tilted, ti, symm, sites, particles, k1, k2, refsite, ndis
        double precision, intent(in)  :: v1, v2, mass, dis 
        character(len=*), intent(in)  :: pattern, cluster, bc, irrep 
        character,        intent(out) :: params*512 

        character                     :: string*256

        if(tilted == 1) then 
            write(string, "(a,'_')") cluster 
            string = trim_name(string)
        end if 
        params = ""    
        if(refsite > 0) then !Name for current files
            if(ti == 0) write(params, "('L=',i0,'N=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refsite=',i0,'.dat')") sites, particles, v1, v2, mass, dis, nDis, bc, pattern, refsite            
            if(ti == 1) write(params, "('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refsite=',i0,'.dat')") sites, particles, k1, k2, v1, v2, mass, dis, nDis, bc, pattern, refsite 
            if(symm == 1) write(params, "('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,a,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refsite=',i0,'.dat')") sites, particles, k1, k2, irrep, v1, v2, mass, dis, nDis, bc, pattern, refsite
        end if 
        params = trim_name(params)
        if(tilted == 1) then 
            params = trim_name(string // params)
        end if 
        params = trim_name(params)    

    end subroutine params_string_dd_cf

    subroutine load_ham(mat_type, params, unit, ti, sites, nnz, nDi, hamDi, ham_i, rc, ham_dp, ham_dc, occupation, exist1, exist2, exist3)
        ! Loads the Hamiltonian from disk.
        ! The Hamiltonian is stored in two files: one for the off-diagonal part and
        ! one for the diagonal part. The off-diagonal part is stored in a file named
        ! "hopping_hamiltonian_<params>" and the diagonal part is stored in a
        ! file named "diagonal_hamiltonian_<params>". The occupation is stored in a
        ! file named "occupation_<params>". 
        implicit none 
        
        integer,                       intent(in)            :: ti, unit, sites  
        character(len=*),              intent(in)            :: mat_type, params

        integer,                       intent(out)           :: nnz, nDi       
        logical,                       intent(out)           :: exist1, exist2, exist3
        integer,          allocatable, intent(out), optional :: hamDi(:, :), ham_i(:, :), rc(:, :), occupation(:, :)
        double precision, allocatable, intent(out), optional :: ham_dp(:)
        double complex,   allocatable, intent(out), optional :: ham_dc(:) 

        integer                                              :: i
        character                                            :: file*256, dir*256

        dir = "hamiltonians/"
        dir = trim_name(dir)

        file = dir // "hopping_hamiltonian_" // params
        file = trim_name(file)
        inquire(file = file, exist = exist1)

        if(exist1) then 
            open(unit, file = file)
            read(unit, *) nnz 

            if(mat_type == "R" .and. ti == 0) then 
                if(allocated(ham_i)) deallocate(ham_i)
                allocate(ham_i(nnz, 3))
                ham_i = 0
                do i = 1, nnz
                    read(unit,*) ham_i(i,1), ham_i(i,2), ham_i(i,3)
                end do 
            else if(mat_type == "R" .and. ti == 1) then 
                if(allocated(ham_dp)) deallocate(ham_dp)
                if(allocated(rc)) deallocate(rc)
                allocate(ham_dp(nnz))
                allocate(rc(nnz, 2))
                ham_dp = 0.d0
                rc = 0

                do i = 1, nnz
                    read(unit,*) ham_dp(i), rc(i,1), rc(i,2)
                end do 
            else if(mat_type == "C" .and. ti == 1) then      
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

        file = ""
        file = dir // "diagonal_hamiltonian_" // params
        file = trim_name(file)
        inquire(file = trim(file), exist = exist2)

        if(exist2) then 
            open(unit, file = trim(file))
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
        
        if(.not.(exist2)) print*,'Diagonal Hamiltonian does not yet exist.'

        file = ""
        file = "hamiltonians/occupation_" // params
        file = trim_name(file)
        inquire(file = file, exist = exist3)
        if(exist3) then 
            open(unit, file = trim(file))
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

        if(.not.(exist3)) print*,'Occupation does not yet exist.'
        if(.not.(exist3)) print*,''

        return 

    end subroutine load_ham

    subroutine print_configuration(conf, ti, k1, k2, v1, v2)
        ! Print the current configuration parameters for debugging and monitoring
        implicit none
        
        integer,          intent(in) :: conf, ti, k1, k2
        double precision, intent(in) :: v1, v2

        print*, '----------------------------------'
        print*, 'Current Configuration: ', conf 
        if(ti == 1) print*, 'k1, k2:', k1, k2
        print*, 'V1:', v1
        print*, 'V2:', v2
        print*, '----------------------------------'

    end subroutine print_configuration

    subroutine print_parameters(par, geo, thr, diag, out)
        ! Print comprehensive simulation parameters for verification and logging
        implicit none
        
        type(sim_params),    intent(in) :: par
        type(geometry),      intent(in) :: geo
        type(thread_params), intent(in) :: thr
        type(diag_params),   intent(in) :: diag
        type(output),        intent(in) :: out  
        
        print*, '----------------------------------'
        print*, 'Simulation Parameters:'
        print*, 'Sites:', geo%sites
        print*, 'Particles:', geo%particles
        print*, 'Dimension:', geo%dim
        print*, 'V1-range:', par%v1min, '-', par%v1max, 'in steps of ', par%dv1
        print*, 'V2-range:', par%v2min, '-', par%v2max, 'in steps of ', par%dv2
        print*, 'Number of eigenvalues(nev):', diag%nev
        print*, 'Number of convergence vectors(ncv):', par%ncv0
        print*, 'Number of disorder configurations:', par%nDis
        if(par%dis > 0.d0) print*, 'Disorder strength:', par%dis
        if(par%mass > 0.d0) print*, 'Sublattice imbalance:', par%mass
        print*, 'Filling fraction:', par%filling
        if(par%t /= 1.d0) print*, 'Hopping strength(t):', par%t
        if(par%g_fact /= 2) print*, 'Electron g-factor:', par%g_fact
        print*, 'Diagonalization method:'
        if(par%feast == 1) then
            print*, '  FEAST'
        else if(par%arpack == 1) then
            print*, '  Arpack'
        else if(par%mkl == 1) then
            print*, '  MKL'
        else if(par%exact == 1) then
            print*, '  Exact diagonalization'
        else
            print*, '  No diagonalization method selected'
        end if
        if(par%othrds > 1) print*, 'Number of OpenMP threads for Arpack:', par%othrds
        if(par%mthrds > 1) print*, 'Number of OpenMP threads for MKL:', par%mthrds
        if(thr%dis_thrds > 1) print*, 'Disorder threads:', thr%dis_thrds
        if(thr%v1_thrds > 1) print*, 'V1 threads:', thr%v1_thrds
        if(thr%v2_thrds > 1) print*, 'V2 threads:', thr%v2_thrds
        if(thr%num_thrds > 1) print*, 'Number of threads for current calculations:', thr%num_thrds
        if(par%corr == 1) print*, 'Correlation calculations enabled.'
        if(par%curr == 1) print*, 'Current-current calculations enabled:', par%curr
        if(par%curr == 1) print*, 'Reference bonds for current calculations:', par%refbonds
        print*, 'degflag threshold:', par%deg
        if(par%states == 1) print*, 'Saving of eigenstates enabled.'
        print*, '----------------------------------'
    
    end subroutine print_parameters

end module io_utils