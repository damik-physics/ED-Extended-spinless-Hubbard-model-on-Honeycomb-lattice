module core_utilities
    ! Module containing essential utilities for simulation setup, sparse matrix operations,
    ! file I/O, mathematical functions, performance monitoring, and resource management
    ! for the honeycomb Hubbard model simulation framework
    use types
    use params
    use symmetries
    implicit none 
    
    interface coo_to_csr
        module procedure coocsr_dp_i4, coocsr_dp_i8, coocsr_dc_i4, coocsr_dc_i8
    end interface coo_to_csr

    interface spmv
        module procedure amux_dp, amux_dc, amux_dp_otf
    end interface spmv

    contains

    subroutine preprocess(par, geo, thrd, out, st)
        ! Main preprocessing routine that initializes simulation parameters, creates output directories,
        ! sets up threading parameters, and prepares the computational environment
                
        implicit none
        type(sim_params),    intent(inout) :: par
        type(geometry),      intent(inout) :: geo
        type(thread_params), intent(inout) :: thrd
        type(output),        intent(inout) :: out
        type(system_state),  intent(inout) :: st

        call read_input(par)
        call setup_output_directory(out%outdir)
        call create_output_subdirs(out%outdir)
        call nsteps(par%v1min, par%v1max, par%dv1, thrd%ndv1)
        call nsteps(par%v2min, par%v2max, par%dv2, thrd%ndv2)
        call characters(par%symm, par%irrep, geo) ! Set characters according to chosen irrep 
        call timing(out%outdir, 0)
        call setvars(par, geo, thrd)  
        call nsteps(par%v1min, par%v1max, par%dv1, thrd%ndv1)
        call nsteps(par%v2min, par%v2max, par%dv2, thrd%ndv2)
        call stepunits(1, thrd%ndv1, thrd%ndv2, thrd%units_2)
        call threadunits(thrd%ndv1, thrd%ndv2, thrd%units)
        call write_parameters_json(out%outdir, par)
        st%first_call = .true. ! Flag to indicate if this is the first call to I/O writing.

    end subroutine preprocess

    subroutine setvars(par, geo, thrd)
        ! Sets basic simulation variables including lattice size based on cluster type
        ! and particle number based on filling fraction
        implicit none

        type(sim_params)    :: par
        type(geometry)      :: geo
        type(thread_params) :: thrd
        character*10        :: clusterl
        integer             :: threads

        !$ call omp_set_num_threads(thrd%othrds)
        !call mkl_set_num_threads(thrd%mthrds)

        if(par%tilted == 0) then 
            geo%sites = 2*par%ucx*par%ucy 
        else if(par%cluster(1:2) == "16") then 
            geo%sites = 16 
        else if(par%cluster(1:2) == "18") then 
            geo%sites = 18 
        else if(par%cluster(1:2) == "20") then 
            geo%sites = 20 
        else if(par%cluster(1:2) == "24") then 
            geo%sites = 24
        else if(par%cluster(1:2) == "30") then 
            geo%sites = 30 
        else if(par%cluster(1:2) == "32") then 
            geo%sites = 32
        else if(par%cluster(1:1) == "L") then 
            clusterl = par%cluster
            read(clusterl(2:3),*) geo%sites
        end if 
        geo%particles = int(par%filling * geo%sites)   

        return

    end subroutine setvars

    subroutine tune_lanczos_parameters(unit, thresh, exact, nevext, nestext, ncv0, dim, full, nev, ncv, nest)
        ! Automatically tunes Lanczos diagonalization parameters based on matrix dimension
        ! and chooses between sparse and dense diagonalization methods
        implicit none
        integer,          intent(in)  :: unit, thresh, exact, nevext, nestext, ncv0
        integer(kind=8),  intent(in)  :: dim
        integer,          intent(out) :: full, nev, ncv, nest
        character*512                 :: filename

        if(exact == 0) then
            if(dim < thresh) then
                full = 1
                nev  = dim
                ncv  = dim
                nest = dim
            else
                full = 0
                if(dim == 1) then
                    nev  = 1
                    ncv  = 1
                    nest = 1
                else
                    ncv  = min(ncv0, dim - 1) !nev+1<=ncv<=dim and ncv maximal
                    nev  = min(nevext, ncv - 3)
                    nest = min(nev, nestext)
                end if
            end if
        else if(exact == 1) then
            full = 1
            nev  = dim
            ncv  = dim
            nest = min(dim, nestext)
        end if


    end subroutine tune_lanczos_parameters

    subroutine stepunits(ndu, ndv, ndv2, units)
        ! Assigns unique unit numbers for file I/O operations in nested parameter loops
        implicit none
        integer,              intent(in)  :: ndu, ndv, ndv2
        integer, allocatable, intent(out) :: units(:,:,:,:)
        integer                           :: i, j, k, l, q

        if(allocated(units)) deallocate(units)
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

    end subroutine stepunits

    subroutine threadunits(ndv, ndv2, units)
        ! Assigns unit numbers for parallel file operations across different thread configurations
        implicit none
        integer,              intent(in)  :: ndv, ndv2
        integer, allocatable, intent(out) :: units(:, :)
        integer                           :: i, j, h
        integer                           :: n_threads_tot

        n_threads_tot = 1

        if(ndv > 0) n_threads_tot = n_threads_tot + ndv
        if(ndv2 > 0) n_threads_tot = n_threads_tot + ndv2

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

    end subroutine threadunits

    subroutine set_thrds(par, th)
        ! Calculates optimal thread distribution across different parallelization levels
        ! (V2, V1, disorder, and computational threads)
        implicit none 
        type(sim_params),    intent(inout) :: par
        type(thread_params), intent(inout) :: th

        ! Assign thread counts from th to local variables(assuming these are module variables)
        th%v2_thrds = max(min(th%ndv2, par%othrds), 1)
        th%v1_thrds = max(min(th%ndv1, int((par%othrds - th%v2_thrds) / th%v2_thrds)), 1) 
        th%dis_thrds = max(min(par%nDis, int((par%othrds - th%v2_thrds * th%v1_thrds) /(th%v2_thrds * th%v1_thrds))), 1)
        th%num_thrds = max(int((par%othrds - th%v2_thrds * th%v1_thrds * th%dis_thrds) /(th%v2_thrds * th%v1_thrds * th%dis_thrds)), 1)

        
        print*, 'Number of V2 threads = ', th%v2_thrds
        print*, 'Number of V1 threads = ', th%v1_thrds
        print*, 'Number of disorder threads = ', th%dis_thrds
        print*, 'Number of threads left = ', th%num_thrds
        if(th%num_thrds < 1) error stop "NO THREADS LEFT AVAILABLE!"
        print*, ''

    end subroutine set_thrds
    
    subroutine read_input(params)
        ! Reads simulation parameters from input.nml namelist file and initializes
        ! derived type components with default or user-specified values
        implicit none
        
        type(sim_params), intent(out) :: params

        integer                       :: ucx, ucy, tilted, ti, k0, symm, p1, p2, p3
        integer                       :: otf,corr, curr, refbonds
        integer                       :: states, deg, feast, arpack, mkl, exact, degflag, dimthresh, nevext, nst, ncv0, nev0, nevmax, othrds, mthrds, nDis, g_fact
        double precision              :: t, dis, mass, filling, dv1, v1min, v1max, dv2, v2min, v2max
        logical                       :: rvec
        character(len=1)              :: bc
        character(len=2)              :: irrep
        character(len=3)              :: cluster
        integer                       :: ios

        ! Set defaults first
        call init_params(params) ! Initializes all components of derived type sim_params to default values
        ! Now transfer default values to local variables 
        ucx = params%ucx 
        ucy = params%ucy 
        bc = params%bc 
        irrep = params%irrep
        p1 = params%p1 
        p2 = params%p2 
        p3 = params%p3 
        cluster = params%cluster
        
        nDis = params%nDis 
        t = params%t 
        dis = params%dis 
        mass = params%mass 
        filling = params%filling 
        dv1 = params%dv1 
        v1min = params%v1min 
        v1max = params%v1max 
        dv2 = params%dv2 
        v2min = params%v2min 
        v2max = params%v2max 

        othrds = params%othrds 
        mthrds = params%mthrds 

        tilted = params%tilted 
        ti = params%ti 
        k0 = params%k0 
        symm = params%symm 
        corr = params%corr 
        curr = params%curr 
        refbonds = params%refbonds 
        states = params%states 
        feast = params%feast 
        arpack = params%arpack 
        mkl = params%mkl 
        exact = params%exact 
        dimthresh = params%dimthresh 
        otf = params%otf 
        degflag = params%degflag 
        g_fact = params%g_fact 

        rvec = params%rvec
        nevext = params%nevext 
        nev0 = params%nev0 
        nevmax = params%nevmax 
        nst = params%nst 
        ncv0 = params%ncv0 

        ! Read input parameters from input.nml file into local variables, overriding the default values initialized above. 
        ! This file should contain the namelist /params_nml/ with all the variables that are subject to change.  
        namelist /params_nml/ ucx, ucy, tilted, cluster, bc, ti, k0, symm, irrep, p1, p2, p3, corr, curr, refbonds, states, deg, feast, arpack, mkl, exact, dimthresh, rvec, nevext, nst, ncv0, nev0, nevmax, otf, degflag, othrds, mthrds, nDis, dis, mass, filling, t, g_fact, dv1, v1min, v1max, dv2, v2min, v2max

        open(unit=10, file='input.nml', status='old', action='read', iostat=ios)
        if(ios == 0) then
            read(10, nml=params_nml)
            close(10)
        else
            print *, 'No input.nml found, using defaults.'
        end if

        ! transfer into derived type
        params%ucx = ucx
        params%ucy = ucy  
        params%bc = bc 
        params%irrep = irrep 
        params%p1 = p1  
        params%p2 = p2  
        params%p3 = p3  
        params%cluster = cluster 

        params%nDis = nDis  
        params%t = t  
        params%dis = dis  
        params%mass = mass  
        params%filling = filling  
        params%dv1 = dv1  
        params%v1min = v1min  
        params%v1max = v1max  
        params%dv2 = dv2  
        params%v2min = v2min  
        params%v2max = v2max  

        params%othrds = othrds  
        params%mthrds = mthrds  

        params%tilted = tilted  
        params%ti = ti  
        params%k0 = k0  
        params%symm = symm  
        params%corr = corr  
        params%curr = curr  
        params%refbonds = refbonds  
        params%states = states  
        params%feast = feast  
        params%arpack = arpack  
        params%mkl = mkl  
        params%exact = exact  
        params%dimthresh = dimthresh  
        params%otf = otf  
        params%degflag = degflag  
        params%g_fact = g_fact  

        params%rvec = rvec 
        params%nevext = nevext  
        params%nev0 = nev0  
        params%nevmax = nevmax  
        params%nst = nst  
        params%ncv0 = ncv0  

    end subroutine read_input

    subroutine init_params(params)
        ! Initializes all simulation parameters to their default values
        implicit none
        type(sim_params), intent(inout) :: params
        
        params%ucx = 3
        params%ucy = 3
        params%bc = 'p'
        params%irrep = 'A1'
        params%p1 = 1
        params%p2 = 1
        params%p3 = 1
        params%cluster = '18A'
        
        params%nDis = 1
        params%t = 1.0d0
        params%dis = 0.0d0
        params%mass = 0.0d0
        params%filling = 0.5
        params%dv1 = 0.1d0
        params%v1min = 0.0d0
        params%v1max = 1.0d0
        params%dv2 = 0.1d0
        params%v2min = 0.0d0
        params%v2max = 1.0d0

        params%othrds = 1
        params%mthrds = 1

        params%tilted = 1
        params%ti = 1
        params%k0 = 1
        params%symm = 0
        params%corr = 0
        params%curr = 0
        params%refbonds = 0
        params%states = 1
        params%feast = 0
        params%arpack = 1
        params%mkl = 0
        params%exact = 0
        params%dimthresh = 1000
        params%otf = 0
        params%degflag = 1
        params%g_fact = 2

        params%rvec = .true.
        params%nevext = 10
        params%nev0   = 200
        params%nevmax = 230
        params%nst = 10
        params%ncv0 = 120

    end subroutine init_params   

    subroutine setup_output_directory(outdir)
        ! Creates timestamped output directory and copies input file for record-keeping
        implicit none
        character(len=:), allocatable, intent(inout) :: outdir
        character(len=64)                             :: timestamp
        integer                                       :: values(8)
        character(len=512)                            :: cmd

        ! Create a timestamp for the output directory
        ! This will create a directory named 'run_YYYYMMDD_HHMMSS'
        ! where YYYYMMDD is the date and HHMMSS is the time
        call date_and_time(values = values)
        write(timestamp, '(I4.4,I2.2,I2.2,"_",I2.2,I2.2,I2.2)') &
             values(1), values(2), values(3), values(5), values(6), values(7)

        outdir = trim('output/run_' // trim(timestamp) // '/')
        cmd = 'mkdir -p ' // trim(outdir)
        call system(cmd)

        ! Also copy the input file for record-keeping
        cmd = 'cp input.nml ' // trim(outdir) // 'input.nml'
        call system(cmd)

    end subroutine setup_output_directory

    subroutine create_output_subdirs(out_dir)
        implicit none
        character(len=*), intent(in) :: out_dir
        character(len=512)            :: cmd
        integer                       :: exitstat

        ! Create subdirectories for different output types
        cmd = 'mkdir -p ' // trim(out_dir) // 'correlations ' // &
                        trim(out_dir) // 'hamiltonians ' // &
                        trim(out_dir) // 'lattice_data ' // &
                        trim(out_dir) // 'logs ' // &
                        trim(out_dir) // 'parameters ' // &
                        trim(out_dir) // 'plots ' // &
                        trim(out_dir) // 'spectra ' // &
                        trim(out_dir) // 'states'

        call execute_command_line(trim(cmd), exitstat=exitstat)
        if(exitstat /= 0) then
            print *, 'Error creating subdirectories. Exit status:', exitstat
            stop 1
        end if
    end subroutine create_output_subdirs

    subroutine write_parameters_json(dir, par)
        ! Writes all simulation parameters to a JSON file for analysis and reproducibility
        implicit none
        character(len=*),     intent(in) :: dir
        type(sim_params),     intent(in) :: par
        integer                          :: unit

        ! Open the JSON file in the specified directory
        open(newunit=unit, file=trim(dir)//'/parameters/parameters.json', status='replace', action='write')

        write(unit,*) '{'

        write(unit,*) '  "ucx": ', par%ucx, ','
        write(unit,*) '  "ucy": ', par%ucy, ','
        write(unit,*) '  "tilted": ', par%tilted, ','

        if(par%tilted == 1) then
            write(unit,*) '  "cluster": "'//trim(par%cluster)//'",'
        else
            write(unit,*) '  "cluster": "0",'
        end if

        write(unit,*) '  "bc": "'//trim(par%bc)//'",'
        write(unit,*) '  "ti": ', par%ti, ','
        write(unit,*) '  "k0": ', par%k0, ','
        write(unit,*) '  "symm": ', par%symm, ','
        write(unit,*) '  "irrep": "'//trim(par%irrep)//'",'
        write(unit,*) '  "p1": ', par%p1, ','
        write(unit,*) '  "p2": ', par%p2, ','
        write(unit,*) '  "p3": ', par%p3, ','
        write(unit,*) '  "corr": ', par%corr, ','
        write(unit,*) '  "curr": ', par%curr, ','
        write(unit,*) '  "refbonds": ', par%refbonds, ','
        write(unit,*) '  "states": ', par%states, ','
        write(unit,*) '  "deg": ', par%deg, ','
        write(unit,*) '  "feast": ', par%feast, ','
        write(unit,*) '  "arpack": ', par%arpack, ','
        write(unit,*) '  "mkl": ', par%mkl, ','
        write(unit,*) '  "exact": ', par%exact, ','
        write(unit,*) '  "dimthresh": ', par%dimthresh, ','
        write(unit,*) '  "rvec": ', merge('true ','false', par%rvec), ','
        write(unit,*) '  "nevext": ', par%nevext, ','
        write(unit,*) '  "nst": ', par%nst, ','
        write(unit,*) '  "ncv0": ', par%ncv0, ','
        write(unit,*) '  "otf": ', par%otf, ','
        write(unit,*) '  "degflag": ', par%degflag, ','
        write(unit,*) '  "nev0": ', par%nev0, ','
        write(unit,*) '  "nevmax": ', par%nevmax, ','
        write(unit,*) '  "othrds": ', par%othrds, ','
        write(unit,*) '  "mthrds": ', par%mthrds, ','
        write(unit,*) '  "ndis": ', par%ndis, ','
        write(unit,*) '  "dis": ', par%dis, ','
        write(unit,*) '  "mass": ', par%mass, ','
        write(unit,*) '  "filling": ', par%filling, ','
        write(unit,*) '  "t": ', par%t, ','
        write(unit,*) '  "g_fact": ', par%g_fact, ','
        write(unit,*) '  "dv1": ', par%dv1, ','
        write(unit,*) '  "v1min": ', par%v1min, ','
        write(unit,*) '  "v1max": ', par%v1max, ','
        write(unit,*) '  "dv2": ', par%dv2, ','
        write(unit,*) '  "v2min": ', par%v2min, ','
        write(unit,*) '  "v2max": ', par%v2max

        write(unit,*) '}'

        close(unit)
    end subroutine write_parameters_json
    
    subroutine cleanup(inner, geo, ham, diag)
        ! Deallocates dynamically allocated arrays to free memory
        ! inner=.true. deallocates inner loop arrays, inner=.false. deallocates outer loop arrays
        implicit none

        type(geometry),         intent(inout) :: geo
        type(hamiltonian_params), intent(inout) :: ham
        type(diag_params),      intent(inout) :: diag
        logical,                intent(in)    :: inner

        ! Deallocate all dynamically allocated arrays used in the program.
        if(inner) then ! Deallocate all dynamically allocated arrays of inner loops
            if(allocated(geo%dplcts))      deallocate(geo%dplcts)
            if(allocated(geo%norm))        deallocate(geo%norm)
            if(allocated(geo%abasis))      deallocate(geo%abasis)
            if(allocated(geo%bbasis))      deallocate(geo%bbasis)
            if(allocated(geo%mombasis))    deallocate(geo%mombasis)
            if(allocated(geo%parities))    deallocate(geo%parities)
            if(allocated(geo%basis_states)) deallocate(geo%basis_states)

            if(allocated(ham%occ))         deallocate(ham%occ)
            if(allocated(ham%rc))          deallocate(ham%rc)
            if(allocated(ham%ham))         deallocate(ham%ham)
            if(allocated(ham%hamDi))       deallocate(ham%hamDi)
            if(allocated(ham%rcOff))       deallocate(ham%rcOff)    
            if(allocated(ham%ham_dc))      deallocate(ham%ham_dc)
            if(allocated(ham%ham_dp))      deallocate(ham%ham_dp)
            if(allocated(ham%hamOff))      deallocate(ham%hamOff)
            if(allocated(ham%hamOff_dp))   deallocate(ham%hamOff_dp)
            if(allocated(ham%hamOff_dc))   deallocate(ham%hamOff_dc)

            if(allocated(diag%evals))       deallocate(diag%evals)
            if(allocated(diag%eigstate))    deallocate(diag%eigstate)
            if(allocated(diag%energies))    deallocate(diag%energies)
            if(allocated(diag%eigstate_dc)) deallocate(diag%eigstate_dc)

        else ! Deallocate all dynamically allocated arrays of outer loops
            if(allocated(geo%xy))          deallocate(geo%xy)
            if(allocated(geo%phases))      deallocate(geo%phases)
            if(allocated(geo%bsites))      deallocate(geo%bsites)
            if(allocated(geo%xtransl))     deallocate(geo%xtransl)
            if(allocated(geo%ytransl))     deallocate(geo%ytransl)
            if(allocated(geo%alattice))    deallocate(geo%alattice)
            if(allocated(geo%blattice))    deallocate(geo%blattice)
            if(allocated(geo%hexsites))    deallocate(geo%hexsites)
            if(allocated(geo%asitesbonds)) deallocate(geo%asitesbonds)
            if(allocated(geo%bsitesbonds)) deallocate(geo%bsitesbonds)
            if(allocated(geo%latticevecs)) deallocate(geo%latticevecs)
        end if 


    end subroutine cleanup

    subroutine timing(dir, start_end)
        ! Records and reports CPU timing information for performance monitoring
        ! start_end=0 starts timer, start_end=1 stops timer and writes results
        implicit none

        integer,          intent(in) :: start_end
        character(len=*), intent(in) :: dir
        character(8),     save       :: datei, datef
        character(10),    save       :: timei, timef
        integer,          save       :: values_i(8), values_f(8)
        real,             save       :: start, finish
        character                    :: file*512


        if(start_end == 0) then
            call cpu_time(start)
            call date_and_time(date=datei,time=timei,values=values_i)
            write(*,"('Calculation started at',x,a,' h',x,a,' min',x,a,' sec')") timei(1:2), timei(3:4), timei(5:6)
            print*, ''
            print*, 'Start date: ',datei(7:8), '.',datei(5:6), '.',datei(1:4)
            print*, ''
        else
            call cpu_time(finish)
            print*, 'Start = ', start
            print*, 'Finish = ', finish
            if(finish - start < 60) then
                write(*,"(' Elapsed CPU time = ',f12.3,' seconds.')") finish-start
                print*, ''
            else if(finish - start < 3600) then
                write(*,"(' Elapsed CPU time = ',f12.3,' minutes.')")(finish-start)/60
                print*, ''
            else
                write(*,"(' Elapsed CPU time = ',f12.3,' hours.')")(finish-start)/3600
                print*, ''
            end if
            call date_and_time(date=datef,time=timef,values=values_f)
            write(*,"(' Calculation started at',x,a,'h',x,a,'min',x,a,'sec')") timei(1:2), timei(3:4), timei(5:6)
            print*, ''
            write(*,"(' Calculation ended at',x,a,'h',x,a,'min',x,a,'sec')") timef(1:2), timef(3:4), timef(5:6)
            print*, ''
            print*, 'Start date: ',datei(7:8), '.',datei(5:6), '.',datei(1:4)
            print*, ''
            print*, 'End date: ',datef(7:8), '.',datef(5:6), '.',datef(1:4)

            file=trim(dir)//"timing.dat"
            open(77,file = trim(file))
            write(77,"(' Elapsed CPU time = ',f20.10,' seconds.')") finish-start
            write(77,"(' Calculation started at',x,a,'h',x,a,'min',x,a,'sec')") timei(1:2), timei(3:4), timei(5:6)
            write(77,"(' Start date:',x,a,'.',x,a,'.',x,a)") ,datei(7:8), datei(5:6), datei(1:4)
            write(77,"(' Calculation ended at',x,a,'h',x,a,'min',x,a,'sec')") timef(1:2), timef(3:4), timef(5:6)
            write(77,"(' End date:',x,a,'.',x,a,'.',x,a)") ,datef(7:8), datef(5:6), datef(1:4)
            close(77)
        end if
    end subroutine timing

    subroutine nsteps(min, max, delta, steps)
        ! Calculates the number of steps in a parameter range given min, max, and step size
        implicit none
        
        double precision, intent(in)  :: min, max, delta
        integer,          intent(out) :: steps

        if(delta <= 1.d0) then
            steps = int(abs(max-min)/delta + delta/2) + 1
        else
            steps = int(abs(max-min)/delta) + 1
        end if


    end subroutine nsteps
    
    subroutine binary_search(dim, s, basis, loc)
        ! Performs binary search to find the location of state s in the sorted basis array
        implicit none

        integer(kind=8), intent(in)  :: dim, basis(dim)
        integer(kind=8), intent(in)  :: s
        integer(kind=8), intent(out) :: loc

        integer(kind=8)              :: left, right, mean

        left = 1
        right = dim
        do while(left <= right)
            mean = floor((left + right)/2.d0)
            if(basis(mean) < s) then
                left = mean + 1
            else if(basis(mean) > s) then
                right = mean - 1
            else
                loc = mean
                return
            end if
        end do
        loc = -1 !If no state is found
        return

    end subroutine binary_search

    subroutine coocsr_dp_i4(nrow,nnz,a,ir,jc,ao,jao,iao)
        ! Converts double precision sparse matrix from COO to CSR format using 4-byte integers
        implicit none

        integer(kind=8),  intent(in)  :: nrow
        integer(kind=8),  intent(in)  :: nnz
        double precision, intent(in)  :: a(*)
        integer,          intent(in)  :: ir(*), jc(*)
        
        double precision, intent(out) :: ao(*)
        integer,          intent(out) :: jao(*), iao(*)

        double precision              :: x
        integer(kind=8)               :: k, i, j, k0
        integer                       :: iad

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

    end subroutine coocsr_dp_i4

    subroutine coocsr_dp_i8(nrow,nnz,a,ir,jc,ao,jao,iao)
        ! Converts double precision sparse matrix from COO to CSR format using 8-byte integers
        implicit none

        integer(kind=8),  intent(in)  :: nrow
        integer(kind=8),  intent(in)  :: nnz
        double precision, intent(in)  :: a(*)
        integer(kind=8),  intent(in)  :: ir(*), jc(*)
        
        double precision, intent(out) :: ao(*)
        integer(kind=8),  intent(out) :: jao(*), iao(*)

        double precision              :: x
        integer(kind=8)               :: k, i, j, k0, iad

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

    end subroutine coocsr_dp_i8
            
    subroutine coocsr_dc_i4(nrow,nnz,a,ir,jc,ao,jao,iao)
        ! Converts double complex sparse matrix from COO to CSR format using 4-byte integers
        implicit none

        integer(kind=8), intent(in)  :: nrow
        integer(kind=8), intent(in)  :: nnz
        integer,         intent(in)  :: ir(*), jc(*)
        double complex,  intent(in)  :: a(*)
    
        integer,         intent(out) :: jao(*),iao(*)
        double complex,  intent(out) :: ao(*)
        
        double complex               :: x
        integer(kind=8)              :: k, k0, i, j
        integer                      :: iad

        !$omp threadprivate(x, k, k0, i, j, iad)

        do k = 1, nrow + 1
            iao(k) = 0
        end do 

        do k = 1, nnz
            iao(ir(k)) = iao(ir(k)) + 1
        end do 

        k = 1
        do j = 1, nrow + 1
            k0 = iao(j)
            iao(j) = k
            k = k + k0
        end do 

        do k = 1, nnz
            i = ir(k)
            j = jc(k)
            x = a(k)
            iad = iao(i)
            ao(iad) =  x
            jao(iad) = j
            iao(i) = iad + 1
        end do 

        do j = nrow, 1, -1
            iao(j+1) = iao(j)
        end do 
        iao(1) = 1
        return

    end subroutine coocsr_dc_i4

    subroutine coocsr_dc_i8(nrow,nnz,a,ir,jc,ao,jao,iao)
        ! Converts double complex sparse matrix from COO to CSR format using 8-byte integers
        implicit none

        integer(kind=8), intent(in)  :: nrow
        integer(kind=8), intent(in)  :: nnz
        integer(kind=8), intent(in)  :: ir(*), jc(*)
        double complex,  intent(in)  :: a(*)
    
        double complex,  intent(out) :: ao(*)
        integer(kind=8), intent(out) :: jao(*),iao(*)
        
        double complex               :: x
        integer(kind=8)              :: k, k0, i, j, iad

        !$omp threadprivate(x, k, k0, i, j, iad)

        do k = 1, nrow + 1
            iao(k) = 0
        end do 

        do k = 1, nnz
            iao(ir(k)) = iao(ir(k)) + 1
        end do 

        k = 1
        do j = 1, nrow + 1
            k0 = iao(j)
            iao(j) = k
            k = k + k0
        end do 

        do k = 1, nnz
            i = ir(k)
            j = jc(k)
            x = a(k)
            iad = iao(i)
            ao(iad) =  x
            jao(iad) = j
            iao(i) = iad + 1
        end do 

        do j = nrow, 1, -1
            iao(j+1) = iao(j)
        end do 
        iao(1) = 1
        return

    end subroutine coocsr_dc_i8

    subroutine amux_dp(threads, n, x, y, a, ja, ia) 
        ! Performs sparse matrix-vector multiplication y = A*x for double precision matrices in CSR format
        implicit none

        integer,          intent(in)  :: threads
        integer(kind=8),  intent(in)  :: n
        integer(kind=8),  intent(in)  :: ja(*), ia(*)
        double precision, intent(in)  :: x(n), a(*)
        double precision, intent(out) :: y(n)
        
        double precision              :: t
        integer(kind=8)               :: i, k

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
        !-----------------------------------------------------------------------

        !$omp parallel do private(i, k, t) num_threads(threads)
        do i = 1, n

            !
            !     compute the inner product of row i with vector x
            !

            t = 0.0d0
            do k = ia(i), ia(i+1) - 1
                t = t + a(k) * x(ja(k))
            end do 

            !
            !     store result in y(i)
            !

            y(i) = t
        end do 
        !$omp end parallel do

        return
    end subroutine amux_dp 

    subroutine amux_dp_otf(threads, dim, sites, nbonds, nnnbonds, t1, v1, v2, basis, bsites, hexsites, x, y) 
        ! Performs matrix-vector multiplication on-the-fly without storing the full Hamiltonian matrix
        ! Generates Hamiltonian elements during the multiplication for memory efficiency
        implicit none

        integer(kind=8),  intent(in)  :: dim
        integer(kind=8),  intent(in)  :: basis(dim)
        integer,          intent(in)  :: threads, sites, nbonds, nnnbonds
        integer,          intent(in)  :: bsites(2, nbonds), hexsites(2, nnnbonds)
        double precision, intent(in)  :: t1, v1, v2
        double precision, intent(in)  :: x(dim)
        double precision, intent(out) :: y(dim)
        
        double precision              :: a, t
        integer(kind=8)               :: i, newst, loc
        integer                       :: s, parity1, parity2

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
            !$ id = omp_get_thread_num()
            
            do s = 1, nbonds  
                if(btest(basis(i), bsites(1, s) - 1) .and. .not. btest(basis(i), bsites(2, s) - 1)) then 
                    newst = ibclr(ibset(basis(i), bsites(2, s) - 1 ), bsites(1, s) - 1 ) !Create on site 2, annihilate on site 1 
                    call binary_search(dim, newst, basis, loc) 
                    if(loc > 0) then 
                        parity1 = popcnt(ibits(basis(i), bsites(1, s), sites) ) !Parity of site1
                        parity2 = popcnt(ibits(ibclr(basis(i), bsites(1, s) - 1), bsites(2, s), sites) ) !Parity of site2
                        a = - t1 *(-1)**(parity1 + parity2) 
                        t = t + a * x(loc) 
                    end if 
                else if(btest(basis(i), bsites(2, s) - 1) .and. .not. btest(basis(i), bsites(1, s) - 1)) then 
                    newst = ibclr(ibset(basis(i), bsites(1, s) - 1 ), bsites(2, s) - 1 ) !Create on site 2, annihilate on site 1 
                    call binary_search(dim, newst, basis, loc) 
                    if(loc > 0) then 
                        parity1 = popcnt(ibits(basis(i), bsites(2, s), sites) ) !Parity of site2
                        parity2 = popcnt(ibits(ibclr(basis(i), bsites(2,s) - 1), bsites(1, s), sites) ) !Parity of site1
                        a = - t1 *(-1)**(parity1 + parity2)  
                        t = t + a * x(loc) 
                    end if                    
                end if 
                if(btest(basis(i), bsites(1, s) - 1) .and. btest(basis(i), bsites(2, s) - 1) ) t = t + 0.25 * v1 * x(i)
                if(.not.(btest(basis(i), bsites(1, s) - 1)) .and. .not.(btest(basis(i), bsites(2, s) - 1)) ) t = t + 0.25 * v1 * x(i)
                if(btest(basis(i), bsites(1, s) - 1) .and. .not.(btest(basis(i), bsites(2, s) - 1)) ) t = t - 0.25 * v1 * x(i)
                if(.not.(btest(basis(i), bsites(1, s) - 1)) .and. btest(basis(i), bsites(2, s) - 1) ) t = t - 0.25 * v1 * x(i)
            end do !bond loop
            do s = 1, nnnbonds  
                if(btest(basis(i), hexsites(1, s) - 1) .and. btest(basis(i), hexsites(2, s) - 1) ) t = t + 0.25 * v2 * x(i)
                if(.not.(btest(basis(i), hexsites(1, s) - 1)) .and. .not.(btest(basis(i), hexsites(2, s) - 1)) ) t = t + 0.25 * v2 * x(i)
                if(btest(basis(i), hexsites(1, s) - 1) .and. .not.(btest(basis(i), hexsites(2, s) - 1)) ) t = t - 0.25 * v2 * x(i)
                if(.not.(btest(basis(i), hexsites(1, s) - 1)) .and. btest(basis(i), hexsites(2, s) - 1) ) t = t - 0.25 * v2 * x(i)
            end do !bond loop
            
            !
            !     store result in y(i)
            !

            y(i) = t
        end do 
        !$omp end parallel do

        return

    end subroutine amux_dp_otf 
        
    subroutine amux_dc(threads, n, x, y, a, ja, ia)
        ! Performs sparse matrix-vector multiplication y = A*x for double complex matrices in CSR format
        implicit none
        
        integer,         intent(in)  :: threads
        integer(kind=8), intent(in)  :: n, ja(*), ia(*)
        double complex,  intent(in)  :: x(*), a(*)
        double complex,  intent(out) :: y(*)
                
        double complex               :: t
        integer(kind=8)              :: i, k

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


        !$omp parallel do private(k,t,i) num_threads(threads)
        do i = 1, n
        !
        !     compute the inner product of row i with vector x
        !
        
            t = 0.0d0
        
            do k = ia(i), ia(i+1)-1
                t = t + a(k)*x(ja(k))
            end do 
        !
        !     store result in y(i)
        !
            y(i) = t

        end do 
        !$omp end parallel do        
        return

    end subroutine amux_dc

    
end module core_utilities
