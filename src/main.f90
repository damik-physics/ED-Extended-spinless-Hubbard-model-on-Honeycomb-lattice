program main

    use types
    use parameters
    use variables
    use input_variables
    use functions
    use file_utils
    use io_routines
    use lattice 
    use basis
    use hamiltonian
    use printing_routines
    ! use symmetries
    ! use utilities
    ! use diagonalization
    ! use correlationfunctions

    implicit none
    type(sim_params) :: par
    type(out_params) :: out_par
    type(system_params) :: sys_par
    type(diag_params) :: diag_par
    type(thread_params) :: thread_par
    type(geometry) :: geo ! Geometry includes lattice-, basis- and symmetry related parameters 
    type(hamiltonian_params) :: ham_par

    
    ! Read input parameters from input.nml file
    ! This file should contain the namelist /params/ with all the variables that are subject to change.  
    ! namelist /params/ ucx, ucy, tilted, cluster, bc, ti, k0, symm, irrep, p1, p2, p3, corr, curr, refbonds, states, deg, feast, arpack, mkl, exact, dimthresh, rvec, nevext, nst, ncv0, otf, degflag, nev0, nevmax, othrds, mthrds, nDis, dis, mass, filling, t, g_fact, dv1, v1min, v1max, dv2, v2min, v2max
    ! open(unit=10, file='input.nml', status='old')
    ! read(10, nml=params)
    ! close(10)
    
    ! call get_environment_variable("pwd", cwd)
    ! call read_input(par)
    ! call setup_output_directory()
    ! call create_output_subdirs(outdir)
    

    ! Now variables are set: those in input.nml are overwritten, others keep defaults.

    ! call characters(symm, irrep, geo) ! Set characters according to chosen irrep 
    ! call characters(symm, irrep, mir, rot, geo%id) ! Set characters according to chosen irrep 
    




    ! call timing(outdir, 0)
    ! call setvars()  
    ! call check_parallel()
    ! call nsteps(par%v1min, par%v1max, par%dv1, ndv1)
    ! call nsteps(par%v2min, par%v2max, par%dv2, ndv2)
    ! call stepunits(1, ndv1, ndv2, units_2)
    ! call threadunits(ndv1, ndv2, units)
    
    
    call define_lattice(outdir, par, geo)
    ! call define_lattice(outdir, tilted, ucx, ucy, nnBonds, nnnBonds, bc, pattern, cluster, bsites, hexsites, geopar%latticevecs, alattice, blattice, xyA, xyB, asitesbonds, bsitesbonds, cntrA, cntrB, nHel, tilt, phases, xy, xtransl, ytransl, reflections, nnnVec)
    

    if(ti == 0 .or. symm == 1 .or. k0 == 1) then 
        geo%k1_max = 0
        geo%k2_max = 0 
    else if(ti == 1 .and. tilted == 0) then 
        geo%k1_max = ucx - 1
        geo%k2_max = ucy - 1
        call save(outdir, cluster, sys_par%unit, tilted, geo%sites, geo%particles, bc, pattern, geo%k1_max, geo%k2_max)
    else if(ti == 1 .and. tilted == 1) then 
        if(geo%nHel == 1) then 
            geo%k1_max = 0
        else if(geo%nHel > 1) then
            if(modulo(dble(geo%sites)/dble((geo%nHel*geo%tilt)), dble(geo%nHel)) == 0.d0) then 
                geo%k1_max = geo%sites/(geo%nHel * geo%tilt) - 1
            else if(modulo(dble(geo%sites)/dble((geo%nHel*geo%tilt)), dble(geo%nHel)) >= 1.d0) then 
                geo%k1_max = geo%sites/geo%tilt - 1
            else if(modulo(dble(geo%sites)/dble((geo%nHel*geo%tilt)), dble(geo%nHel)) < 1.d0) then 
                geo%k2_max = ceiling(geo%sites/(geo%nHel * geo%tilt * modulo(dble(geo%sites)/dble((geo%nHel*geo%tilt)), 1.d0))) - 1
            end if 
        end if 
        geo%k2_max = geo%sites/(2*geo%nHel) - 1
        call save(outdir, cluster, sys_par%unit, tilted, geo%sites, geo%particles, bc, pattern, geo%k1_max, geo%k2_max)
    end if


    do k1 = 0, geo%k1_max
    do k2 = 0, geo%k2_max
    if(ti == 1) print('(1x,3(a,i0))'), '---------------------- K1 = ', k1, ' K2 = ', k2, ' -------------------------'
    
    if((k1 .ne. 0) .or. (k2 .ne. 0) .or. (geo%id == 2)) then 
        diag_par%type = "C"
        print('(1x, a)'), 'Matrix type: Complex'
    else 
        diag_par%type = "R"
        print('(1x, a)'), 'Matrix type: Real'
    end if    

    call make_basis(ti, tilted, pattern, nnnVec, sites, geo%particles, dim, symm, ucx, ucy, l1, l2, basis_states, abasis, bbasis, tilt, nHel, k1, k2, xtransl, ytransl, geo%id, mir, rot, refl, c6, period, norm, orbsize, orbits2D, phases2D, norm2D)

    ! call print_params(nev, ndeg, nHel, tilt, k1, k2)
    call printing()
    ! if(otf == 0) call generate_hamiltonian(othrds, ti, sys_par%unit, param_list, sites, nnBonds, bsites, dim, basis_states, hamOff, nOff, k1, k2, tilted, nHel, tilt, l2, l1, ucx, ucy, orbsize, norm, norm2D, orbits2D, phases2D, xtransl, ytransl, symm, geo%id, mir, rot, refl, c6, t, rcOff, rcDi, parities, dplcts, hamOff_dp, hamDi_dp, nDi_dp, hamOff_dc, hamDi_off_dc, hamDi_dc, nnnBonds, hexsites, hamDi, occ, nDiOff)    
    if(otf == 0) call generate_hamiltonian()

    if (dim == 0) then
        nev = int(dim, kind=4)
        ncv = int(dim, kind=4)
        print*, 'No basis states found.'
        goto 116
    end if

    v2_thrds = max(min(ndv2, othrds), 1)
    v1_thrds = max(min(ndv, int((othrds - v2_thrds) / v2_thrds)), 1) 
    dis_thrds = max(min(nDis, int((othrds - v2_thrds * v1_thrds) / (v2_thrds * v1_thrds))), 1)
    num_thrds = max(int((othrds - v2_thrds * v1_thrds * dis_thrds) / (v2_thrds * v1_thrds * dis_thrds)), 1)

    print*, 'Number of V2 threads = ', v2_thrds
    print*, 'Number of V1 threads = ', v1_thrds
    print*, 'Number of disorder threads = ', dis_thrds
    print*, 'Number of threads left = ', num_thrds
    if(num_thrds < 1) error stop "NO THREADS LEFT AVAILABLE!"
    print*, ''
  
    !$omp parallel do default(firstprivate) shared(sites, geo%particles, dim, nev, ncv, occ, parities, dplcts, hamOff, hamDi, basis, hamOff_dp, hamDi_d, rcOff, rcDi, hamOff_dc, num_thrds) num_threads(v2_thrds)
    do nv2 = 0, ndv2-1, 1 ! Loop over V2 values
        v2 = v2min + nv2*dv2
        !$ thread_num = omp_get_thread_num()
            
        !$omp parallel do default(firstprivate) shared(sites, geo%particles, dim, nev, ncv, occ, parities, dplcts, hamOff, hamDi, basis, hamOff_dp, hamDi_d, rcOff, rcDi, hamOff_dc, num_thrds) num_threads(v1_thrds)
        do nv = 0, ndv-1, 1 ! Loop over V1 values
            v1 = v1min + nv*dv1

            !$ thread_num_2 = omp_get_thread_num()
            unit = 11 + units(thread_num + 1, thread_num_2 + 1)

            call ncv_from_nev(unit, param_list, dimthresh, exact, nevext, nst, ncv0, dim, full, nev, ncv, nest)

            ! !$omp critical 
            do conf = 1, nDis !Disorder loop 
                call print(conf, ti, k1, k2, v1, v2)
                ! call diagonalization(outdir, conf, nev, ncv, full, v1, v2, othrds, param_list, type, basis, bsites, hexsites, occ, nOff, nDi, nDi_d, hamOff, hamDi, hamOff_dp, hamDi_d, hamOff_dc, hamDi_c, hamDi_off_c, ham, ham_dc, ham_d, ham_dc, norm, rcOff, rcDi, rc, parities, dplcts, nnz, ndeg, sys_par%unit, nest, mode, energies, eigstate, eigstate_dc, gs, gs_c)
                call diagonalize()
    
                if(curr == 0) goto 15          
            
                
                if(nv2 > 0) ndeg = 2
                call currentcorrelations(outdir, ti, tilted, nHel, tilt, num_thrds, conf, degflag, full, feast, mkl, arpack, symm, sites, l2, l1, ucx, ucy, k1, k2, geo%id, mir, rot, nDis, ndeg, refbonds, cntrA, cntrB, sys_par%unit, dim, alattice, blattice, xyA, xyB, xtransl, ytransl, refl, c6, basis, v1, v2, dis, nnnVec, norm, eigstate, eigstate_dc)
                15 continue 

                if(corr == 0) goto 16
                do refsite = 1, sites
                    call ddcf(outdir, refsite, conf, sys_par%unit, dim, sites, k1, k2, nDis, dis, v1, v2, basis, eigstate(1:dim, 1))
                end do 
                16 continue

            end do !Disorder loops
            ! !$omp end critical 


            100 format(1000(F30.20))
            101 format(2000(F30.20))
        
    
        end do !vloop
        !!$omp end parallel do

    end do !v2loop
    !!$omp end parallel do
   
    call cleanup(inner=.true.)

    
    116 continue

    call timing(outdir, 1)
    ! call printparams(nev, ndeg, nHel, tilt, k1, k2)
    call printing()
    call cleanup(inner=.false.)
    end do 
    end do !Momentum loops 


    

end program main