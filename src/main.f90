program main

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
    use symmetries
    use utilities
    ! use diagonalization
    ! use correlationfunctions

    implicit none
    

    ! Read input parameters from input.nml file
    ! This file should contain the namelist /params/ with all the variables that are subject to change.  
    namelist /params/ ucx, ucy, tilted, cluster, bc, ti, k0, symmetrize, irrep, p1, p2, p3, corr, curr, refbonds, states, deg, feast, arpack, mkl, exact, dimthresh, rvec, nevext, n_st, ncv0, otf, degeneracy, nev0, nevmax, othrds, mthrds, nDis, dis, mass, filling, t, g_fact, dv1, v1min, v1max, dv2, v2min, v2max, debug
    open(unit=10, file='input.nml', status='old')
    read(10, nml=params)
    close(10)
    
    ! call get_environment_variable("pwd", cwd)

    call setup_output_directory()
    call create_output_subdirs(outdir)
    

    ! Now variables are set: those in input.nml are overwritten, others keep defaults.

    call characters(symmetrize, irrep, mir, rot, id) ! Set characters according to chosen irrep 
    




    call datetime(outdir, 0)
    call setvars()  
    call check_parallel()
    call nsteps(v1min, v1max, dv1, ndv)
    call nsteps(v2min, v2max, dv2, ndv2)
    call stepunits(1, ndv, ndv2, units_2)
    call threadunits(ndv, ndv2, units)

    call define_lattice(outdir, tilted, ucx, ucy, nnBonds, nnnBonds, bc, pattern, cluster, bsites, hexsites, latticevecs, alattice, blattice, xyA, xyB, asitesbonds, bsitesbonds, cntrA, cntrB, nHel, tilt, phases, xy, xtransl, ytransl, reflections, nnnVec)
    

    if(ti == 0 .or. symmetrize == 1 .or. k0 == 1) then 
        k1_max = 0
        k2_max = 0 
    else if(ti == 1 .and. tilted == 0) then 
        k1_max = ucx - 1
        k2_max = ucy - 1
        call save(outdir, cluster, unit, tilted, sites, particles, bc, pattern, k1_max, k2_max)
    else if(ti == 1 .and. tilted == 1) then 
        if(nHel == 1) then 
            k1_max = 0
        else if(nHel > 1) then
            if(modulo(dble(sites)/dble((nHel*tilt)), dble(nhel)) == 0.d0) then 
                k1_max = sites/(nHel * tilt) - 1
            else if(modulo(dble(sites)/dble((nHel*tilt)), dble(nhel)) >= 1.d0) then 
                k1_max = sites/tilt - 1
            else if(modulo(dble(sites)/dble((nHel*tilt)), dble(nhel)) < 1.d0) then 
                k2_max = ceiling(sites/(nHel * tilt * modulo(dble(sites)/dble((nHel*tilt)), 1.d0))) - 1
            end if 
        end if 
        k2_max = sites/(2*nHel) - 1
        call save(outdir, cluster, unit, tilted, sites, particles, bc, pattern, k1_max, k2_max)
    end if


    do k1 = 0, k1_max
    do k2 = 0, k2_max
    if(ti == 1) print('(1x,3(a,i0))'), '---------------------- K1 = ', k1, ' K2 = ', k2, ' -------------------------'
    
    if((k1 .ne. 0) .or. (k2 .ne. 0) .or. (id == 2)) then 
        type = "C"
        print('(1x, a)'), 'Matrix type: Complex'
    else 
        type = "R"
        print('(1x, a)'), 'Matrix type: Real'
    end if    

    call make_basis(ti, tilted, pattern, nnnVec, sites, particles, dim, symmetrize, ucx, ucy, l1, l2, basis_states, abasis, bbasis, tilt, nHel, k1, k2, xtransl, ytransl, id, mir, rot, refl, c6, period, norm, orbsize, orbits2D, phases2D, norm2D)

    ! call print_params(nev, ndeg, nHel, tilt, k1, k2)
    call printing()
    ! if(otf == 0) call generate_hamiltonian(othrds, ti, unit, param_list, sites, nnBonds, bsites, dim, basis_states, hamOff, nOff, k1, k2, tilted, nHel, tilt, l2, l1, ucx, ucy, orbsize, norm, norm2D, orbits2D, phases2D, xtransl, ytransl, symmetrize, id, mir, rot, refl, c6, t, rcOff, rcDi, parities, dplcts, hamOff_dp, hamDi_dp, nDi_dp, hamOff_dc, hamDi_off_dc, hamDi_dc, nnnBonds, hexsites, hamDi, occ, nDiOff)    
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
  
    ! !$omp parallel do default(firstprivate) shared(sites, particles, dim, nev, ncv, occ, parities, dplcts, hamOff, hamDi, basis, hamOff_dp, hamDi_d, rcOff, rcDi, hamOff_dc, num_thrds) num_threads(v2_thrds)
    ! do nv2 = 0, ndv2-1, 1 ! Loop over V2 values
    !     v2 = v2min + nv2*dv2
    !     !$ thread_num = omp_get_thread_num()
            
    !     !$omp parallel do default(firstprivate) shared(sites, particles, dim, nev, ncv, occ, parities, dplcts, hamOff, hamDi, basis, hamOff_dp, hamDi_d, rcOff, rcDi, hamOff_dc, num_thrds) num_threads(v1_thrds)
    !     do nv = 0, ndv-1, 1 ! Loop over V1 values
    !         v1 = v1min + nv*dv1

    !         !$ thread_num_2 = omp_get_thread_num()
    !         unit = 11 + units(thread_num + 1, thread_num_2 + 1)

    !         call nevncv(unit, param_list, dimthresh, exact, nevext, n_st, ncv0, dim, full, nev, ncv, nest)

    !         ! !$omp critical 
    !         do conf = 1, nDis !Disorder loop 
    !             call print(conf, ti, k1, k2, v1, v2)
    !             ! call diagonalization(outdir, conf, nev, ncv, full, v1, v2, othrds, param_list, type, basis, bsites, hexsites, occ, nOff, nDi, nDi_d, hamOff, hamDi, hamOff_dp, hamDi_d, hamOff_dc, hamDi_c, hamDi_off_c, ham, ham_dc, ham_d, ham_dc, norm, rcOff, rcDi, rc, parities, dplcts, nnz, ndeg, unit, nest, mode, energies, eigstate, eigstate_dc, gs, gs_c)
    !             call diagonalize()
    
    !             if(curr == 0) goto 15          
            
                
    !             if(nv2 > 0) ndeg = 2
    !             call currentcorrelations(outdir, ti, tilted, nHel, tilt, num_thrds, conf, degeneracy, full, feast, mkl, arpack, symmetrize, sites, l2, l1, ucx, ucy, k1, k2, id, mir, rot, nDis, ndeg, refbonds, cntrA, cntrB, unit, dim, alattice, blattice, xyA, xyB, xtransl, ytransl, refl, c6, basis, v1, v2, dis, nnnVec, norm, eigstate, eigstate_dc)
    !             15 continue 

    !             if(corr == 0) goto 16
    !             do refsite = 1, sites
    !                 call ddcf(outdir, refsite, conf, unit, dim, sites, k1, k2, nDis, dis, v1, v2, basis, eigstate(1:dim, 1))
    !             end do 
    !             16 continue

    !         end do !Disorder loops
    !         ! !$omp end critical 


    !         100 format(1000(F30.20))
    !         101 format(2000(F30.20))
        
    
    !     end do !vloop
    !     !!$omp end parallel do

    ! end do !v2loop
    ! !!$omp end parallel do
   


    
    116 continue

    call datetime(outdir, 1)
    ! call printparams(nev, ndeg, nHel, tilt, k1, k2)
    call printing()
    end do 
    end do !Momentum loops 


    

end program main