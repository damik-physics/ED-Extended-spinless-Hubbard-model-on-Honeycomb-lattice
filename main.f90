program main

    use routines
    use variables
    

    implicit none
    
    integer :: k1, k2, k1_max, k2_max, orbsize, ndeg, unit, refsite, nOff, nDi, conf, nest, nnz, nDi_d, cntrA, cntrB, l1, l2, i, j, v1_thrds = 1, v2_thrds = 1, dis_thrds = 1, thread_num, thread_num_2, num_threads, nev, ncv, nbb, nnnbb, ndv, ndv2, nv, nv2, full, nHel, tilt 
    integer, allocatable :: units(:,:), units2(:,:,:,:), bsites(:,:), hexsites(:,:), occ(:,:), hamOff(:,:), hamDi(:,:), rc(:,:), rcOff(:,:), rcDi(:), parities(:), dplcts(:), phases(:), xy(:,:), xyA(:,:), xyB(:,:), latticevecs(:), alattice(:,:), blattice(:,:), asitesbonds(:,:), bsitesbonds(:,:),  xtransl(:,:), ytransl(:,:), reflections(:), refl(:,:), c6(:)
    integer(kind=8), allocatable :: basis(:), abasis(:), bbasis(:), mombasis(:), period(:), orbits2D(:,:,:)

    double precision :: v1, v2, mir(6), rot(5), id 
    double precision, allocatable :: gs(:), nnnVec(:,:), energies(:), eigstate(:,:), evals(:), ham(:), ham_dp(:,:), hamOff_dp(:), hamDi_dp(:), norm(:), norm2D(:,:)
    double complex, allocatable :: eigstate_dc(:,:), hamOff_dc(:), ham_dc(:), gs_c(:), ham_dc(:,:), phases2D(:,:,:), hamDi_c(:,:), hamDi_off_c(:)
    character :: type*1, dir*100!, prms*200
    character, parameter :: mode*2 = 'SA'!'SA'
    ! Uncomment on HPC cluster
    ! character :: gclusterext*1, tiext*1, othrdsext*3, mthrdsext*3, ucxext*1, ucyext*1
    ! character :: vext*16, v2minext*16, v2maxext*16, disext*16, nDisext*16
    external :: mkl_set_num_threads, mkl_get_max_thrds
    
    i = 0
    j = 0
    call characters(symmetrize, irrep, mir, rot, id)
    if(spartan == 0) dir = "output/"
    if(spartan == 1) dir = "../output/"
    dir = trim_name(dir)

    call datetime(dir, 0)

    write(*,*) '-----------------------'
    write(*,*) ' Exact diagonalization '
    write(*,*) '-----------------------'


    ! Uncomment on HPC cluster
        ! if(spartan == 1) then 
        !     call get_environment_variable('tilted',gclusterext)
        !     call get_environment_variable('cluster',cluster)
        !     call get_environment_variable('ti',tiext)
        !     call get_environment_variable('othrds',othrdsext)
        !     call get_environment_variable('mthrds',mthrdsext)
        !     call get_environment_variable('ucx',ucxext)
        !     call get_environment_variable('ucy',ucyext)
        !     call get_environment_variable('v1min',vext)
        !     call get_environment_variable('v2min',v2minext)
        !     call get_environment_variable('v2max',v2maxext)
        !     call get_environment_variable('dis',disext)
        !     call get_environment_variable('nDis',nDisext)
        !     read(gclusterext,*) tilted
        !     read(tiext,*) ti
        !     read(othrdsext,*) othrds
        !     read(mthrdsext,*) mthrds
        !     read(ucxext,*) ucx
        !     read(ucyext,*) ucy
        !     read(vext,*) vmin
        !     read(vext,*) vmax
        !     read(v2minext,*) v2min
        !     read(v2maxext,*) v2max
        !     read(disext,*) dis
        !     read(nDisext,*) nDis 
    ! end if 

    call setvars()  
    call parallelcheck()
    call nsteps(vmin, vmax, dv, ndv)
    call nsteps(v2min, v2max, dv2, ndv2)
    call stepunits(1, ndv, ndv2, units2)
    call threadunits(ndv, ndv2, units)

    call lattice(dir, tilted, ucx, ucy, nbb, nnnbb, bc, pattern, cluster, bsites, hexsites, latticevecs, alattice, blattice, xyA, xyB, asitesbonds, bsitesbonds, cntrA, cntrB, nHel, tilt, phases, xy, xtransl, ytransl, reflections, nnnVec)
    

    if(ti == 0 .or. symmetrize == 1 .or. k0 == 1) then 
        k1_max = 0
        k2_max = 0 
    else if(ti == 1 .and. tilted == 0) then 
        k1_max = ucx - 1
        k2_max = ucy - 1
        call savemom(dir, cluster, unit, tilted, sites, particles, bc, pattern, k1_max, k2_max)
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
        call savemom(dir, cluster, unit, tilted, sites, particles, bc, pattern, k1_max, k2_max)
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

    call make_basis(ti, tilted, pattern, nnnVec, sites, particles, dim, symmetrize, ucx, ucy, l1, l2, basis, abasis, bbasis, tilt, nHel, k1, k2, xtransl, ytransl, id, mir, rot, refl, c6, period, norm, orbsize, orbits2D, phases2D, norm2D)   
    call printparams(nev, ndeg, nHel, tilt, k1, k2)

    if(otf == 0) call hamiltonian(spartan, othrds, ti, unit, prms, sites, nbb, bsites, dim, basis, hamOff, nOff, k1, k2, tilted, nHel, tilt, l2, l1, ucx, ucy, orbsize, norm, norm2D, orbits2D, phases2D, xtransl, ytransl, symmetrize, id, mir, rot, refl, c6, t, rcOff, rcDi, parities, dplcts, hamOff_dp, hamDi_d, nDi_d, hamOff_dc, hamDi_off_c, hamDi_c, nnnbb, hexsites, hamDi, occ, nDi)    


    if (dim == 0) then
        nev = dim
        ncv = dim
        goto 116
    end if

    num_threads = 1

    print*, 'Number of V2 thrds = ', v2_thrds
    print*, 'Number of V1 thrds = ', v1_thrds
    print*, 'Number of disorder thrds = ', dis_thrds
    print*, 'Number of thrds left = ', num_threads
    if(num_threads < 1) error stop "NO thrds LEFT AVAILABLE! NUMBER OF thrds"
    print*, ''

    !--------------------------!
    !         V2LOOP           !
    !--------------------------!
    
    ! !$omp parallel do default(firstprivate) shared(sites, particles, dim, nev, ncv, occ, parities, dplcts, hamOff, hamDi, basis, hamOff_dp, hamDi_d, rcOff, rcDi, hamOff_dc, num_threads) num_threads(v2_thrds)
    do nv2 = 0, ndv2-1, 1
        v2 = v2min + nv2*dv2
        !$ thread_num = omp_get_thread_num()
            
        !--------------------------!
        !          V1LOOP          !
        !--------------------------!

        ! !$omp parallel do default(firstprivate) shared(sites, particles, dim, nev, ncv, occ, parities, dplcts, hamOff, hamDi, basis, hamOff_dp, hamDi_d, rcOff, rcDi, hamOff_dc, num_threads) num_threads(v1_thrds)
        do nv = 0, ndv-1, 1
            v1 = vmin + nv*dv

            !$ thread_num_2 = omp_get_thread_num()
            unit = 11 + units(thread_num+1,thread_num_2+1)
            
            call par2(unit, k1, k2, 0, v1, v2, prms)
            call nevncv(unit, prms, dimthresh, exact, nevext, n_st, ncv0, dim, full, nev, ncv, nest)

            ! !$omp critical 
            do conf = 1, nDis !Disorder loop 
                print*,''
                if(ti == 1) print('(1x,3(a,i0),3(a,f5.3))'), '---------------------- Conf. = ', conf,' k1 = ',k1,' k2 = ',k2,' V2 = ', v2, ' V = ', v1, ' -------------------------'
                if(ti == 0) print('(1x,a,i0,3(a,f5.3))'), '---------------------- Conf. = ', conf,' V2 = ', v2, ' V = ', v1, ' -------------------------'
                print*, ''
                call diagonalization(dir, conf, nev, ncv, full, v1, v2, othrds, prms, type, basis, bsites, hexsites, occ, nOff, nDi, nDi_d, hamOff, hamDi, hamOff_dp, hamDi_d, hamOff_dc, hamDi_c, hamDi_off_c, ham, ham_dc, ham_d, ham_dc, norm, rcOff, rcDi, rc, parities, dplcts, nnz, ndeg, unit, nest, mode, energies, eigstate, eigstate_dc, gs, gs_c)


                if(curr == 0) goto 15          
                
                num_threads = 1
                
                if(nv2 > 0) ndeg = 2
                call currentcorrelations(dir, ti, tilted, nHel, tilt, num_threads, conf, degeneracy, full, feast, mkl, arpack, symmetrize, sites, l2, l1, ucx, ucy, k1, k2, id, mir, rot, nDis, ndeg, refbonds, cntrA, cntrB, unit, dim, alattice, blattice, xyA, xyB, xtransl, ytransl, refl, c6, basis, v1, v2, dis, nnnVec, norm, eigstate, eigstate_dc)
                15 continue 

                if(corr == 0) goto 16
                do refsite = 1, sites
                    call ddcf(dir, refsite, conf, unit, dim, sites, k1, k2, nDis, dis, v1, v2, basis, eigstate(1:dim, 1))
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
   

    if(allocated(rc))           deallocate(rc)
    if(allocated(ham))          deallocate(ham)
    if(allocated(ham_dc))        deallocate(ham_dc)
    if(allocated(evals))        deallocate(evals)
    if(allocated(ham_d))        deallocate(ham_d)
    if(allocated(norm))         deallocate(norm)
    if(allocated(hamDi))        deallocate(hamDi)
    if(allocated(rcOff))       deallocate(rcOff)    
    if(allocated(abasis))       deallocate(abasis)
    if(allocated(bbasis))       deallocate(bbasis)
    if(allocated(hamOff))       deallocate(hamOff)
    if(allocated(mombasis))     deallocate(mombasis)
    if(allocated(parities))     deallocate(parities)
    if(allocated(eigstate))     deallocate(eigstate)
    if(allocated(energies))     deallocate(energies)
    if(allocated(hamOff_dp))     deallocate(hamOff_dp)
    if(allocated(hamOff_dc))     deallocate(hamOff_dc)
    if(allocated(occ))          deallocate(occ)
    if(allocated(dplcts))       deallocate(dplcts)
    if(allocated(eigstate_dc))   deallocate(eigstate_dc)
    if(allocated(basis))        deallocate(basis)
    
    116 continue

    call datetime(dir, 1)
    call printparams(nev, ndeg, nHel, tilt, k1, k2)
    end do 
    end do !Momentum loops 

    if(allocated(xy))          deallocate(xy)
    if(allocated(xtransl))     deallocate(xtransl)
    if(allocated(ytransl))     deallocate(ytransl)
    if(allocated(alattice))    deallocate(alattice)
    if(allocated(blattice))    deallocate(blattice)
    if(allocated(phases))      deallocate(phases)
    if(allocated(bsites))      deallocate(bsites)
    if(allocated(hexsites))    deallocate(hexsites)
    if(allocated(asitesbonds)) deallocate(asitesbonds)
    if(allocated(bsitesbonds)) deallocate(bsitesbonds)
    if(allocated(latticevecs)) deallocate(latticevecs)
    
    

end program

!------------------------------!
!        END OF PROGRAM        !
!------------------------------!
      