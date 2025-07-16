module file_utils
    implicit none

contains


    subroutine preprocess(par, geo, thrd, out)
        ! Preprocessing routine to set up necessary parameters and variables
        use vars
        use params
        use input_vars
        use types
        
        implicit none
        type(sim_params), intent(inout) :: par
        type(geometry), intent(inout) :: geo
        type(thread_params), intent(inout) :: thrd
        type(out_params), intent(inout) :: out
        call read_input(par)
        call setup_output_directory(out%outdir)
        call create_output_subdirs(out%outdir)
        call nsteps(par%v1min, par%v1max, par%dv1, thrd%ndv1)
        call nsteps(par%v2min, par%v2max, par%dv2, thrd%ndv2)
        call characters(par%symm, par%irrep, geo) ! Set characters according to chosen irrep 
        call timing(out%outdir, 0)
        call setvars()  
        call check_parallel()
        call nsteps(par%v1min, par%v1max, par%dv1, thrd%ndv1)
        call nsteps(par%v2min, par%v2max, par%dv2, thrd%ndv2)
        call stepunits(1, thrd%ndv1, thrd%ndv2, thrd%units_2)
        call threadunits(thrd%ndv1, thrd%ndv2, thrd%units)
        

    end subroutine preprocess


    subroutine set_thrds()
        use types
        implicit none 
        type(sim_params) :: par
        type(thread_params) :: th
        ! Assign thread counts from th to local variables (assuming these are module variables)
        th%v2_thrds = max(min(th%ndv2, th%othrds), 1)
        th%v1_thrds = max(min(th%ndv1, int((th%othrds - th%v2_thrds) / th%v2_thrds)), 1) 
        th%dis_thrds = max(min(par%nDis, int((th%othrds - th%v2_thrds * th%v1_thrds) / (th%v2_thrds * th%v1_thrds))), 1)
        th%num_thrds = max(int((th%othrds - th%v2_thrds * th%v1_thrds * th%dis_thrds) / (th%v2_thrds * th%v1_thrds * th%dis_thrds)), 1)

        
        print*, 'Number of V2 threads = ', th%v2_thrds
        print*, 'Number of V1 threads = ', th%v1_thrds
        print*, 'Number of disorder threads = ', th%dis_thrds
        print*, 'Number of threads left = ', th%num_thrds
        if(th%num_thrds < 1) error stop "NO THREADS LEFT AVAILABLE!"
        print*, ''

    end subroutine set_thrds
    
    subroutine read_input(params)
        use types
        implicit none
        type(sim_params), intent(out) :: params

        ! local copies for reading from namelist
        integer :: ucx, ucy, tilted, ti, k0, symm, p1, p2, p3, corr, curr, refbonds, states, deg, feast, arpack, mkl, exact, dimthresh, nev, nst, ncv0, otf, degflag, nev0, nevmax, othrds, mthrds, nDis, g_fact
        double precision :: t, dis, mass, filling, dv1, v1min, v1max, dv2, v2min, v2max
        logical :: rvec
        character(len=1) :: bc
        character(len=2) :: irrep
        character(len=3) :: cluster

        integer :: ios

        ! Read input parameters from input.nml file
        ! This file should contain the namelist /params/ with all the variables that are subject to change.  
        namelist /params_nml/ ucx, ucy, tilted, cluster, bc, ti, k0, symm, irrep, p1, p2, p3, corr, curr, refbonds, states, deg, feast, arpack, mkl, exact, dimthresh, rvec, nev, nst, ncv0, otf, degflag, nev0, nevmax, othrds, mthrds, nDis, dis, mass, filling, t, g_fact, dv1, v1min, v1max, dv2, v2min, v2max
        open(unit=10, file='input.nml', status='old', action='read')
        read(10, nml=params_nml)
        close(10)


        ! set defaults first
        call init_params(params)
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
        nev = params%nev 
        nev0 = params%nev0 
        nevmax = params%nevmax 
        nst = params%nst 
        ncv0 = params%ncv0 

        ! attempt to read from file
        open(unit=10, file='input.nml', status='old', action='read', iostat=ios)
        if (ios == 0) then
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
        params%nev = nev  
        params%nev0 = nev0  
        params%nevmax = nevmax  
        params%nst = nst  
        params%ncv0 = ncv0  

    end subroutine read_input


    subroutine setup_output_directory(outdir)
        use types 
        character(len=*) :: outdir
        implicit none
        character(len=64) :: timestamp
        integer :: values(8)
        character(len=512) :: cmd

        ! Create a timestamp for the output directory
        ! This will create a directory named 'run_YYYYMMDD_HHMMSS'
        ! where YYYYMMDD is the date and HHMMSS is the time
        call date_and_time(values = values)
        write(timestamp, '(I4.4,I2.2,I2.2,"_",I2.2,I2.2,I2.2)') &
             values(1), values(2), values(3), values(5), values(6), values(7)

        outdir = 'output/run_' // trim(timestamp) // '/'
        cmd = 'mkdir -p ' // trim(outdir)
        call system(cmd)

        ! Also copy the input file for record-keeping
        cmd = 'cp input.nml ' // trim(outdir) // 'input.nml'
        call system(cmd)

    end subroutine setup_output_directory

    subroutine create_output_subdirs(out_dir)
        implicit none
        character(len=*), intent(in) :: out_dir
        character(len=512) :: cmd
        integer :: exitstat
        ! Create subdirectories for different output types
        cmd = 'mkdir -p ' // trim(out_dir) // '/correlations ' // &
                        trim(out_dir) // '/hamiltonians ' // &
                        trim(out_dir) // '/lattice_data ' // &
                        trim(out_dir) // '/logs ' // &
                        trim(out_dir) // '/parameters ' // &
                        trim(out_dir) // '/plots ' // &
                        trim(out_dir) // '/spectra ' // &
                        trim(out_dir) // '/states'

        call execute_command_line(trim(cmd), exitstat=exitstat)
        if (exitstat /= 0) then
            print *, 'Error creating subdirectories. Exit status:', exitstat
            stop 1
        end if
    end subroutine create_output_subdirs


    subroutine write_parameters_json(ucx, ucy, tilted, cluster, bc, V, V2, filling, irrep)
        use types
        implicit none
    
        integer, intent(in) :: ucx, ucy, tilted
        character(*), intent(in) :: cluster, bc, irrep
        real(8), intent(in) :: V, V2, filling
        integer :: unit
        type(out_params) :: out 
        
        open(newunit=unit, file=trim(out%outdir)//'parameters.json', status='replace')
        write(unit,'(A)') '{'
        write(unit,'(A,I0,A)') '"ucx": ', ucx, ','
        write(unit,'(A,I0,A)') '"ucy": ', ucy, ','
        write(unit,'(A,I0,A)') '"tilted": ', tilted, ','
        write(unit,'(A,"""",A,"""",A)') '"cluster": "', trim(cluster), '",'
        write(unit,'(A,"""",A,"""",A)') '"bc": "', trim(bc), '",'
        write(unit,'(A,F6.3,A)') '"V": ', V, ','
        write(unit,'(A,F6.3,A)') '"V2": ', V2, ','
        write(unit,'(A,F6.3,A)') '"filling": ', filling, ','
        write(unit,'(A,"""",A,"""")') '"irrep": "', trim(irrep), '"'
        write(unit,'(A)') '}'
        close(unit)
    end subroutine write_parameters_json

    
end module file_utils
