module observables
    !>  This module contains routines for calculating current-current correlations 
    !>  in the honeycomb Hubbard model simulation. Density-density correlations are added soon. 
     
    use types
    use functions 
    use core_utilities
    use corr_writer
    use symmetries
    implicit none

    interface current_cf
        module procedure current_reg_dp
        module procedure current_tilted_dp
        module procedure current_dc
    end interface current_cf

    contains 
    
    subroutine current_correlations(par, geo, diag, out, st, thr)
        implicit none 
        type(sim_params),             intent(in) :: par
        type(geometry),               intent(in) :: geo
        type(diag_params),            intent(in) :: diag
        type(output),                 intent(in) :: out
        type(system_state),           intent(in) :: st
        type(thread_params),          intent(in) :: thr

        integer                                   :: refbond, case, nbonds, numthrds, thread, i, j, nrefb, l11, l22
        integer,          allocatable             :: site_coords(:,:), lattice(:,:)
        double precision                          :: current
        double precision, allocatable             :: bondcurrent(:)
        character                                 :: sl*1, fname1*512, fname2*512 
        character(len=:), allocatable             :: dir
        
        !$omp threadprivate(par%refbond, j, nrefb, current, bondcurrent, sl)

        dir = out%outdir
        if(par%ti == 1 .and. par%tilted == 1) then 
            l11 = geo%l2
            l22 = geo%l1 
        else 
            l11 = par%ucx 
            l22 = par%ucy 
        end if 

        print*,'Calculate current current correlations...'
        print*,''
        do i = 1, 2 ! Select sublattice A or B
            if(i == 1) then 
                sl = "A"
                fname1 = "current_A"
                fname2 = "bond_current_A"
                nbonds = geo%cntrA
                if(allocated(lattice)) deallocate(lattice)
                if(allocated(site_coords)) deallocate(site_coords)
                allocate(lattice(5, nbonds))
                allocate(site_coords(4, nbonds))
                lattice     = geo%alattice
                site_coords = geo%xyA
            else if(i == 2) then 
                sl = "B"
                fname1 = "current_B"
                fname2 = "bond_current_B"
                nbonds = geo%cntrB
                if(allocated(lattice)) deallocate(lattice)
                if(allocated(site_coords)) deallocate(site_coords)
                allocate(lattice(5, nbonds))
                allocate(site_coords(4, nbonds))
                lattice     = geo%blattice
                site_coords = geo%xyB
            end if 

            if(par%refbonds > 0)  nrefb = par%refbonds 
            if(par%refbonds == 0) nrefb = nbonds

            !$ call omp_set_num_threads(thrd%num_threads)
            !$omp parallel do default(firstprivate) shared(geo%sites, geo%particles, geo%dim, geo%basis_states) num_threads(thrd%num_threads)
            do j = 1, nrefb
                refbond = lattice(4,j)       
                !$ numthrds = omp_get_num_threads()
                !$ thread = omp_get_thread_num()
                if(par%ti == 1 .and. ((st%k1 .ne. 0) .or. (st%k2 .ne. 0))) then ! Complex eigenstates
                    ! Regular and tilted cluster
                    call current_cf(par%tilted, geo%nHel, geo%tilt, st%k1, st%k2, par%t, geo%dim, geo%sites, l11, l22, par%ti, par%symm, geo%mir, geo%rot, geo%basis_states, lattice, site_coords, geo%xtransl, geo%ytransl, geo%refl, geo%c6, nbonds, refbond, diag%eigstate_dc(1:geo%dim,1), geo%nnnVec, geo%norm, bondcurrent, current)
                else ! Real eigenstates
                    ! Regular (rectangular) cluster  
                    if(par%tilted==0) call current_reg_dp(par%t, geo%dim, geo%sites, l11, l22, par%ti, par%symm, geo%mir, geo%rot, geo%basis_states, lattice, geo%xtransl, geo%ytransl, geo%refl, geo%c6, nbonds, refbond, diag%eigstate(1:geo%dim,1), geo%nnnVec, geo%norm, bondcurrent, current)
                    ! Tilted cluster
                    if(par%tilted==1) call current_tilted_dp(geo%nHel, geo%tilt, par%t, geo%dim, geo%sites, l11, l22, par%ti, par%symm, geo%mir, geo%rot, geo%basis_states, lattice, geo%xtransl, geo%ytransl, geo%refl, geo%c6, nbonds, refbond, diag%eigstate(1:geo%dim,1), geo%nnnVec, geo%norm, bondcurrent, current)    
                end if 
                ! Store results in csv file  
                call append_corr_csv(fname1, dir, par%ti, par%ndis, st%k1, st%k2, st%conf, st%v1, st%v2, current)
                call append_corr_csv(fname2, dir, par%ti, par%ndis, st%k1, st%k2, st%conf, st%v1, st%v2, bondcurrent)
            end do !Loop over bonds 
            !$omp end parallel do 
        end do

        return 
    end subroutine current_correlations

    subroutine current_reg_dp(t, dim, sites, Lx, Ly, ti, symm, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, refbond, psi, nnnVec, norm, bcurrent, current)
        implicit none
        integer,                      intent(in)  :: nbonds, refbond, sites, Lx, Ly, ti, symm
        integer,          allocatable, intent(in) :: bsites(:,:), xtransl(:,:), ytransl(:,:), refl(:,:), c6(:)
        integer(kind=8),              intent(in)  :: dim, basis(dim)
        double precision,             intent(in)  :: t
        double precision,             intent(in)  :: mir(6), rot(5)
        double precision,             intent(in)  :: psi(:) 
        double precision, allocatable, intent(in) :: nnnVec(:, :), norm(:) 
        double precision,             intent(out) :: current
        double precision, allocatable, intent(out):: bcurrent(:)

        integer(kind=8)                           :: loc, newst, rep, state
        integer                                   :: order, x, y, info, i, j, k, l1, l2, nd1, nd2
        integer                                   :: site1, site2, refsite1, refsite2, phase, parity
        double precision                          :: dist
        double precision                          :: sign, signrep, signstate, signDir
        double precision, dimension(2)            :: vec, refvec
        double precision, dimension(2,3)          :: vecs
        double precision, allocatable             :: psiprime(:), evals(:), gs(:)

        !$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, sublattice_current_temp, bcurrent, psiprime, psiprimeij, psiprimeprime)

        vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
        vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
        vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

        if(allocated(bcurrent)) deallocate(bcurrent)
        if(allocated(psiprime)) deallocate(psiprime)
        
        allocate(bcurrent(nbonds))
        allocate(psiprime(dim))

        refvec   = vecs(1:2,bsites(5,refbond))
        refsite1 = bsites(1,refbond)
        refsite2 = bsites(2,refbond)  
        bcurrent = 0.d0 
        current  = 0.d0

        if(symm == 1) then 
            order = 12
        else 
            order = 1
        end if 

        do i = 1, nbonds !Loop over bonds    
            psiprime = 0.0d0  
            site1    = bsites(1,i)
            site2    = bsites(2,i)
            phase    = bsites(3,i)   
            refsite1 = bsites(1,refbond)
            refsite2 = bsites(2,refbond)
            vec      = vecs(1:2,bsites(5,i))
            dist     = dot_product(vec, refvec)
            dist     = 1.d0

            if(site1 == refsite1 .or. site1 == refsite2) cycle
            if(site2 == refsite1 .or. site2 == refsite2) cycle !Exclude bonds which share a site with the reference bond

            do j = 1, dim !Loop over basis states 
                do k = 1, order !Loop over point groups operators 
                    if(k == 1) then 
                        state = basis(j)
                    else if(k <= 7) then 
                        call reflect(basis(j), basis(j), sites, refl(k-1, 1:sites), signstate, info, state)
                        signstate = signstate * mir(k-1)
                    else if(8 <= k) then 
                        call c6n(basis(j), basis(j), sites, k-7, c6, signstate, info, state)
                        signstate = signstate * rot(k-7)
                    end if 
                    
                    do x = 1, Lx
                        if(ti == 1) call xtranslate(state, Ly, sites, Lx, xtransl, signstate, state)
                        do y = 1, Ly
                            if(ti == 1) call ytranslate(state, Ly, 0, sites, Lx, ytransl, signstate, state)
                            if(btest(state, refsite1-1) .and. .not.(btest(state, refsite2-1))) then !Forward hopping: 1 -> 2
                                newst   = hopping(basis(j), refsite1, refsite2) !Create on refsite 2, annihilate on refsite 1   
                                parity  = popcnt(ibits(state, refsite1, sites)) + popcnt(ibits(ibclr(state, refsite1-1), refsite2, sites))
                                signDir = 1 
                            else if(btest(state, refsite2-1) .and. .not.(btest(state, refsite1-1))) then !Backwards hopping: 2 -> 1 
                                newst   = hopping(basis(j), refsite2, refsite1) !Create on refsite 1, annihilate on refsite 2
                                parity  = popcnt(ibits(state, refsite2, sites)) + popcnt(ibits(ibclr(state, refsite2-1), refsite1, sites))
                                signDir = -1 
                            else 
                                cycle 
                            end if 
                            if(btest(state, site1-1) .and. .not.(btest(state, site2-1))) then !Forward forward hopping: 1 -> 2       
                                parity = parity + popcnt(ibits(newst, site1, sites)) + popcnt(ibits(ibclr(newst, site1-1), site2, sites))                     
                                newst  = hopping(basis(j), site1, site2) !Create on site 2, annihilate on site 1                             
                            else if(btest(state, site2-1) .and. .not.(btest(state, site1-1))) then !Backwards hopping: 2 -> 1
                                parity  = parity + popcnt(ibits(newst, site2, sites)) + popcnt(ibits(ibclr(newst, site2-1), site1, sites)) 
                                newst   = ibclr(ibset(newst, site1-1), site2-1) !Create on site 1, annihilate on site 2
                                newst   = hopping(basis(j), site1, site2) !Create on site 1, annihilate on site 2
                                signDir = signDir * (-1)
                            else 
                                cycle 
                            end if
                            if(ti == 1) call representative_rect(newst, sites, Lx, Ly, symm, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                            if(ti == 0) rep = newst
                            call binary_search(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            sign = signstate * signrep * signDir * (-1)**parity 
                            if(ti == 1) psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(j)
                            if(ti == 0) psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j)                               
                        end do !y
                    end do !x 
                end do !G 
            end do !j 
            bcurrent(i) = bcurrent(i) + dist * dot_product(psi, psiprime) / (Lx * Ly * order) !Bond current J_ij 
            current = current + dble(phase) * bcurrent(i)       
        end do !Loop over bonds    

        deallocate(psiprime)

        return

    end subroutine current_reg_dp

    subroutine current_tilted_dp(nHel, tilt, t, dim, sites, Lx, Ly, ti, symm, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, refbond, psi, nnnVec, norm, bcurrent, current)
        ! Calculates current-current structure factor (current) and local current current correlations (bondcurrent) on the sublattice for a given reference bond.
        ! This subroutine is used for the tilted case.
        implicit none
        integer,                      intent(in)  :: nHel, tilt, nbonds, refbond, sites, Lx, Ly, ti, symm
        integer,          allocatable, intent(in) :: bsites(:,:), xtransl(:,:), ytransl(:,:), refl(:,:), c6(:)
        integer(kind=8),              intent(in)  :: dim, basis(dim)
        double precision,             intent(in)  :: t    
        double precision,             intent(in)  :: psi(:) 
        double precision, allocatable, intent(in) :: nnnVec(:, :), norm(:) 
        double precision,             intent(in)  :: mir(6), rot(5) 
        double precision,             intent(out) :: current
        double precision, allocatable, intent(out):: bcurrent(:)
        
        integer(kind=8)                           :: loc, newst, rep, state
        integer                                   :: x, y, info, i, j, k, l1, l2
        integer                                   :: order, site1, site2, refsite1, refsite2, phase, parity
        double precision                          :: dist
        double precision                          :: sign, signrep, signstate, signDir 
        double precision, dimension(2)            :: vec, refvec
        double precision, dimension(2,3)          :: vecs
        double precision, allocatable             :: psiprime(:)

        !$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, sublattice_current_temp, bcurrent, psiprime, psiprimeij, psiprimeprime)

        vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
        vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
        vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

        if(allocated(bcurrent)) deallocate(bcurrent)
        if(allocated(psiprime)) deallocate(psiprime)
        allocate(bcurrent(nbonds))
        allocate(psiprime(dim))

        refvec    = vecs(1:2,bsites(5,refbond))
        refsite1  = bsites(1,refbond)
        refsite2  = bsites(2,refbond)  
        current   = 0.d0
        bcurrent  = 0.d0 
        sign      = 1.d0 
        signrep   = 1.d0
        signstate = 1.d0
        signDir   = 1.d0

        if(symm == 1 .and. ti == 1) then 
            order = 12
        else 
            order = 1
        end if 

        do i = 1, nbonds !Loop over bonds    
            psiprime = 0.0d0  
            site1    = bsites(1,i)
            site2    = bsites(2,i)
            phase    = bsites(3,i)   
            refsite1 = bsites(1,refbond)
            refsite2 = bsites(2,refbond)
            vec      = vecs(1:2,bsites(5,i))
            dist     = dot_product(vec, refvec)
            dist     = 1.d0

            if(site1 == refsite1 .or. site1 == refsite2) cycle
            if(site2 == refsite1 .or. site2 == refsite2) cycle !Exclude bonds which share a site with the reference bond

            do j = 1, dim !Loop over basis states 

                do k = 1, order 
                    if(k == 1) then 
                        state = basis(j)
                    else if(k <= 7) then 
                        call reflect(basis(j), basis(j), sites, refl(k-1, 1:sites), signstate, info, state)
                        signstate = signstate * mir(k-1)
                    else if(8 <= k) then 
                        call c6n(basis(j), basis(j), sites, k-7, c6, signstate, info, state)
                        signstate = signstate * rot(k-7)
                    end if 
                    
                    do x = 1, Lx ! Loop over x-layers 
                        if(ti == 1) call xtranslate(state, nHel, sites, Lx, xtransl, signstate, state) ! Translate state in x-direction

                        do y = 1, Ly ! Loop over y-layers
                            if(ti == 1) call ytranslate(state, nHel, tilt, sites, Lx, ytransl, signstate, state) ! Translate state in y-direction
                            if(btest(state, refsite1-1) .and. .not.(btest(state, refsite2-1))) then !Forward hopping: 1 -> 2
                                newst = ibclr(ibset(state, refsite2-1), refsite1-1) !Create on refsite 2, annihilate on refsite 1   
                                parity = popcnt(ibits(state, refsite1, sites)) + popcnt(ibits(ibclr(state, refsite1-1), refsite2, sites))
                                signDir = 1 
                            else if(btest(state, refsite2-1) .and. .not.(btest(state, refsite1-1))) then !Backwards hopping: 2 -> 1 
                                newst = ibclr(ibset(state, refsite1-1), refsite2-1) !Create on refsite 1, annihilate on refsite 2
                                parity = popcnt(ibits(state, refsite2, sites)) + popcnt(ibits(ibclr(state, refsite2-1), refsite1, sites))
                                signDir = -1 
                            else 
                                cycle 
                            end if 

                            if(btest(state, site1-1) .and. .not.(btest(state, site2-1))) then !Forward forward hopping: 1 -> 2       
                                parity = parity + popcnt(ibits(newst, site1, sites)) + popcnt(ibits(ibclr(newst, site1-1), site2, sites))                     
                                newst = ibclr(ibset(newst, site2-1), site1-1) !Create on site 2, annihilate on site 1                             
                            else if(btest(state, site2-1) .and. .not.(btest(state, site1-1))) then !Backwards hopping: 2 -> 1
                                parity = parity + popcnt(ibits(newst, site2, sites)) + popcnt(ibits(ibclr(newst, site2-1), site1, sites)) 
                                newst = ibclr(ibset(newst, site1-1), site2-1) !Create on site 1, annihilate on site 2
                                signDir = signDir * (-1)
                            else 
                                cycle 
                            end if 

                            if(ti == 1) call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, symm, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                            if(ti == 0) rep = newst
                            call binary_search(dim, rep, basis, loc) 
                            
                            if(loc <= 0) cycle   
                            if(abs(signstate) > 1) stop 'signstate > 1' 
                            if(abs(signrep) > 1)   stop 'signrep > 1' 
                            if(abs(signDir) > 1)   stop 'signDir > 1' 
                            
                            sign = signstate * signrep * signDir * (-1)**parity 
                            if(ti == 1) psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(j)
                            if(ti == 0) psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j)
                            
                        end do !y
                    end do !x 
                end do !G 
                
            end do !j 
            bcurrent(i) = bcurrent(i) + dist * dot_product(psi, psiprime) / (Lx * Ly * order) !Bond current J_ij 
            current = current + dble(phase) * bcurrent(i)       
        end do !Loop over bonds    

        deallocate(psiprime)

        return

    end subroutine current_tilted_dp

    subroutine current_dc(tilted, nHel, tilt, k1, k2, t, dim, sites, Lx, Ly, ti, symm, mir, rot, basis, bsites, xy, xtransl, ytransl, refl, c6, nbonds, refbond, psi, nnnVec, norm, bondcurrent, current)
        ! Based on complex eigenstates, this subroutine calculates current-current structure factor (current) and local current current correlations (bondcurrent) on the sublattice for a given reference bond.
        ! This subroutine is used for both regular and tilted clusters.
        implicit none
        integer,                      intent(in)  :: tilted, nHel, tilt, k1, k2, nbonds, refbond, sites, Lx, Ly, ti, symm
        integer(kind=8),              intent(in)  :: dim
        integer(kind=8),              intent(in)  :: basis(dim)
        integer,          allocatable, intent(in) :: bsites(:,:)
        integer,          allocatable, intent(in) :: xy(:,:)
        integer,          allocatable, intent(in) :: xtransl(:,:), ytransl(:,:)
        integer,          allocatable, intent(in) :: refl(:,:), c6(:)
        double precision,             intent(in)  :: t 
        double precision, allocatable, intent(in) :: nnnVec(:,:)
        double precision,             intent(in)  :: mir(6), rot(5)
        double complex,               intent(in)  :: psi(dim)
        double precision, allocatable, intent(in) :: norm(:)
        double precision,             intent(out) :: current
        double precision, allocatable, intent(out):: bondcurrent(:)
        
        integer(kind=8)                           :: loc, newst, rep, cntr
        integer                                   :: site1, site2, phase 
        integer                                   :: x, y, i, j, l1, l2, l, l11, l22 
        integer                                   :: refsite1, refsite2, parity1, parity2, parity3, parity4
        double precision                          :: dist
        double precision                          :: sign, signt 
        double precision, dimension(2)            :: vec, refvec
        double precision, dimension(2,3)          :: vecs
        double complex,   allocatable             :: psiprime(:)
        character*250                             :: name, name2 

        !$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, sublattice_current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)

        vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
        vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
        vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

        if(allocated(bondcurrent)) deallocate(bondcurrent)
        allocate(bondcurrent(nbonds))
        if(allocated(psiprime)) deallocate(psiprime)
        allocate(psiprime(dim))

        psiprime = 0.0d0
        refvec   = vecs(1:2,bsites(5,refbond))
        refsite1 = bsites(1,refbond)
        refsite2 = bsites(2,refbond)       
        current = 0.d0
        bondcurrent = 0.d0 
        
        do i = 1, nbonds !Loop over bonds    
            bondcurrent(i) = 0.0d0
            psiprime       = 0.0d0  
            site1          = bsites(1,i)
            site2          = bsites(2,i)
            refsite1       = bsites(1,refbond)
            refsite2       = bsites(2,refbond)    
            phase          = bsites(3,i)   
            vec            = vecs(1:2,bsites(5,i))
            dist           = dot_product(vec, refvec)
            do x = 0, Lx-1 !Loop over x-layers
                do y = 0, Ly-1 !Loop over y-layers
                    if (site1 == refsite1 .or. site2 == refsite2 &
                        .or. site1 == refsite2 .or. site2 == refsite1) go to 110 !Exclude bonds which share a site with the reference bond
                    do j = 1, dim !Loop over basis states 
                        if(btest(basis(j), refsite1-1) .and. .not.(btest(basis(j), refsite2-1))) then !Forward hopping from refsites 1 -> 2
                            newst = ibclr(ibset(basis(j), refsite2-1), refsite1-1) !Create on refsite 2, annihilate on refsite 1   
                            if(btest(basis(j), site1-1) .and. .not.(btest(basis(j), site2-1))) then !Forward forward hopping from sites 1 -> 2       
                                ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                                newst = ibclr(ibset(newst, site2-1), site1-1) !Create on site 2, annihilate on site 1 
                                if(tilted == 1) then 
                                    call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, symm, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signt)
                                    l11 = Ly
                                    l22 = Lx
                                else if(tilted == 0) then 
                                    call representative_rect(newst, sites, Lx, Ly, symm, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                    l11 = Lx
                                    l22 = Ly
                                end if
                                call binary_search(dim, rep, basis, loc) 
                                if(loc <= 0) cycle   
                                parity1 = popcnt(ibits(basis(j), refsite1, sites)) !Parity of c_refs1
                                parity2 = popcnt(ibits(ibclr(basis(j), refsite1-1), refsite2, sites)) !Parity of cd_refs2 
                                parity3 = popcnt(ibits(ibset( ibclr(basis(j), refsite1-1), refsite2-1), site1, sites)) !Parity of c_site1
                                parity4 = popcnt(ibits(ibclr(ibset(ibclr(basis(j), refsite1-1), refsite2-1), site1 - 1), site2, sites)) !Parity of cd_site2
                                sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                if(ti == 1) psiprime(j) = psiprime(j) + (-1) * sign * signt * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(loc) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                if(ti == 0) psiprime(j) = psiprime(j) + (-1) * sign * signt * (-(t)**2) * psi(loc) 
                            else if(btest(basis(j), site2-1) .and. .not.(btest(basis(j), site1-1))) then !Forward backwards hopping from sites 2 -> 1                   
                                ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                                newst = ibclr(ibset(newst, site1-1), site2-1) !Create on site 1, annihilate on site 2
                                if(tilted == 1) then 
                                    call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, symm, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signt)                                
                                    l11 = Ly
                                    l22 = Lx
                                else if(tilted == 0) then 
                                    call representative_rect(newst, sites, Lx, Ly, symm, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.                      
                                    l11 = Lx
                                    l22 = Ly
                                end if
                                call binary_search(dim, rep, basis, loc) 
                                if(loc <= 0) cycle   
                                parity1 = popcnt(ibits(basis(j), refsite1, sites)) !Parity of c_refs1
                                parity2 = popcnt(ibits(ibclr(basis(j), refsite1-1), refsite2, sites)) !Parity of cd_refs2 
                                parity3 = popcnt(ibits(ibset( ibclr(basis(j), refsite1-1), refsite2-1), site2, sites)) !Parity of c_site2
                                parity4 = popcnt(ibits(ibclr(ibset(ibclr(basis(j), refsite1-1), refsite2-1), site2-1), site1, sites)) !Parity of cd_site1
                                sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                if(ti == 1) psiprime(j) = psiprime(j) + (-1) * sign * signt * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(loc) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                if(ti == 0) psiprime(j) = psiprime(j) + (-1) * sign * signt * (-(t)**2) * psi(loc) 
                            else 
                                cycle 
                            end if     
                        else if( btest( basis(j), refsite2-1) .and. .not.( btest( basis(j), refsite1-1))) then !Backwards hopping from refsites 2 -> 1 
                            newst = ibclr(ibset(basis(j), refsite1-1), refsite2-1) !Create on refsite 1, annihilate on refsite 2
                            if(btest(basis(j), site1-1) .and. .not.(btest(basis(j), site2-1))) then !Backwards forward hopping from sites 1 -> 2       
                                ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                                newst = ibclr(ibset(newst, site2-1), site1-1) !Create on site 2, annihilate on site 1 
                                if(tilted == 1) then 
                                    call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, symm, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signt)                                
                                    l11 = Ly
                                    l22 = Lx
                                else if(tilted == 0) then 
                                    call representative_rect(newst, sites, Lx, Ly, symm, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.                   
                                    l11 = Lx
                                    l22 = Ly
                                end if
                                call binary_search(dim, rep, basis, loc) 
                                if(loc <= 0) cycle                              
                                parity1 = popcnt(ibits(basis(j), refsite2, sites)) !Parity of c_refs2
                                parity2 = popcnt(ibits(ibclr(basis(j), refsite2-1), refsite1, sites)) !Parity of cd_refs1 
                                parity3 = popcnt(ibits(ibset( ibclr(basis(j), refsite2-1), refsite1-1), site1, sites)) !Parity of c_site1
                                parity4 = popcnt(ibits(ibclr(ibset(ibclr(basis(j), refsite2-1), refsite1-1), site1 - 1), site2, sites)) !Parity of cd_site2
                                sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                if(ti == 1) psiprime(j) = psiprime(j) + (-1) * sign * signt * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(loc) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                if(ti == 0) psiprime(j) = psiprime(j) + (-1) * sign * signt * (-(t)**2) * psi(loc) 
                            else if(btest(basis(j), site2-1) .and. .not.(btest(basis(j), site1-1))) then !Backwards backwards hopping from sites 2 -> 1                   
                                ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                                newst = ibclr(ibset(newst, site1-1), site2-1) !Create on site 1, annihilate on site 2
                                if(tilted == 1) then 
                                    call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, symm, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signt)                                
                                    l11 = Ly
                                    l22 = Lx
                                else if(tilted == 0) then 
                                    call representative_rect(newst, sites, Lx, Ly, symm, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.   
                                    l11 = Lx
                                    l22 = Ly
                                end if
                                call binary_search(dim, rep, basis, loc) 
                                if(loc <= 0) cycle   
                                parity1 = popcnt(ibits(basis(j), refsite2, sites)) !Parity of c_refs2
                                parity2 = popcnt(ibits(ibclr(basis(j), refsite2-1), refsite1, sites)) !Parity of cd_refs1 
                                parity3 = popcnt(ibits(ibset( ibclr(basis(j), refsite1-1), refsite2-1), site2, sites)) !Parity of c_site2
                                parity4 = popcnt(ibits(ibclr(ibset(ibclr(basis(j), refsite1-1), refsite2-1), site2-1), site1, sites)) !Parity of cd_site1
                                sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                if(ti == 1) psiprime(j) = psiprime(j) + (-1) * sign * signt * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(loc) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                if(ti == 0) psiprime(j) = psiprime(j) + (-1) * sign * signt * (-(t)**2) * psi(loc) 
                            else 
                                cycle 
                            end if     
                        else 
                            cycle 
                        end if 
                    end do !Loop over basis states
                    110 continue 
                    !a1 translations 
                    cntr = 0
                    do l = 1, nbonds
                        if(bsites(1,l) == site1 .and. bsites(5, l) == 1) then 
                            site1 = bsites(2,l)
                            cntr = cntr + 1
                        else if(bsites(1,l) == site2 .and. bsites(5, l) == 1) then
                            site2 = bsites(2,l)
                            cntr = cntr + 1
                        else if(bsites(1,l) == refsite1 .and. bsites(5, l) == 1) then
                            refsite1 = bsites(2,l)
                            cntr = cntr + 1
                        else if(bsites(1,l) == refsite2 .and. bsites(5, l) == 1) then
                            refsite2 = bsites(2,l)
                            cntr = cntr + 1
                        else if(cntr == 4) then 
                            exit 
                        end if
                    end do 
                end do 
                !a2 translations 
                cntr = 0
                do l = 1, nbonds
                    if(bsites(1,l) == site1 .and. bsites(5, l) == 2) then 
                        site1 = bsites(2,l)
                        cntr = cntr + 1
                    else if(bsites(1,l) == site2 .and. bsites(5, l) == 2) then
                        site2 = bsites(2,l)
                        cntr = cntr + 1
                    else if(bsites(1,l) == refsite1 .and. bsites(5, l) == 2) then
                        refsite1 = bsites(2,l)
                        cntr = cntr + 1
                    else if(bsites(1,l) == refsite2 .and. bsites(5, l) == 2) then
                        refsite2 = bsites(2,l)
                        cntr = cntr + 1
                    else if(cntr == 4) then 
                        exit 
                    end if
                end do 
            end do 
            bondcurrent(i) = bondcurrent(i) + dist * dot_product(psi, psiprime) / (Lx*Ly)
            current = current + dble(phase) * bondcurrent(i)
        end do !Loop over bonds 
        
        deallocate(psiprime)

        return

    end subroutine current_dc



end module observables