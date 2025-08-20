module observables
    use types
    use functions 
    use file_utils
    use corr_writer
    implicit none

    interface current_cf
        module procedure current_reg_dp
        module procedure current_tilted_dp
        module procedure current_dc
    end interface current_cf

    interface representative
        module procedure representative_reg
        module procedure representative_tilted
    end interface representative

    contains 
    
subroutine current_correlations(par, geo, out, st, thr)
    implicit none 
    type(sim_params), intent(in) :: par
    type(geometry), intent(in) :: geo
    type(output), intent(in) :: out
    type(system_state), intent(in) :: st
    type(thread_params), intent(in) :: thr
    character(len=*), parameter :: dir = out%outdir

    integer                       :: refbond = 0, case = 0, nbonds = 0, numthrds = 0, thread = 0, i = 0, j = 0, nrefb = 0, l11 = 0, l22 = 0
    integer,          allocatable :: site_coords(:,:), lattice(:,:)
    double precision              :: current = 0.d0 
    double precision, allocatable :: bondcurrent(:)
    character                     :: sl*1, fname1*512, fname2*512 
    class(*)                      :: gs(:)
    
    !$omp threadprivate(par%refbond, j, nrefb, current, bondcurrent, sl)
    
    if(ti == 1 .and. par%tilted == 1) then 
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

        if(allocated(gs)) deallocate(gs)
        allocate(gs(geo%dim))
        if(ti == 1 .and. ((k1 .ne. 0) .or. (k2 .ne. 0))) then 
            gs   = diag%eigstate_dc(1:geo%dim, 1)
            case = 2
        else 
            gs   = diag%eigstate(1:geo%dim, 1)
            case = 1
        end if

        !$ call omp_set_num_threads(thrd%num_threads)
        !$omp parallel do default(firstprivate) shared(geo%sites, geo%particles, geo%dim, geo%basis_states) num_threads(thrd%num_threads)
        do j = 1, nrefb
            refbond = lattice(4,j)       
            !$ numthrds = omp_get_num_threads()
            !$ thread = omp_get_thread_num()
            select case(case)
            case(1) ! Real eigenstates
                ! Regular (rectangular) cluster  
                if(par%tilted==0) call current_cf(par%t, geo%dim, geo%sites, l11, l22, par%ti, par%symm, geo%id, geo%mir, geo%rot, geo%basis_states, lattice, geo%xtransl, geo%ytransl, geo%refl, geo%c6, nbonds, refbond, gs, geo%nnnVec, geo%norm, bondcurrent, current)
                ! Tilted cluster
                if(par%tilted==1) call current_cf(geo%nHel, geo%tilt, par%t, geo%dim, geo%sites, l11, l22, par%ti, par%symm, geo%id, geo%mir, geo%rot, geo%basis_states, lattice, geo%xtransl, geo%ytransl, geo%refl, geo%c6, nbonds, refbond, gs, geo%nnnVec, geo%norm, bondcurrent, current)
            case(2) ! Complex eigenstates
                ! Regular and tilted cluster
                call current_cf(par%tilted, geo%nHel, geo%tilt, k1, k2, t, geo%dim, geo%sites, l11, l22, geo%basis_states, lattice, site_coords, geo%xtransl, geo%ytransl, nbonds, refbond, gs, geo%nnnVec, geo%norm, bondcurrent, current)
            end select
            
            ! Store results in csv file  
            call append_corr_csv(fname1, dir, par%ti, par%ndis, st%k1, st%k2, st%conf, st%v1, st%v2, current)
            call append_corr_csv(fname2, dir, par%ti, par%ndis, st%k1, st%k2, st%conf, st%v1, st%v2, bondcurrent)
            
        end do !Loop over bonds 
        !$omp end parallel do 
    
    end do

    return 
end subroutine current_correlations


subroutine current_reg_dp(t, dim, sites, Lx, Ly, ti, symm, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, refbond, psi, nnnVec, norm, bcurrent, current)
    implicit none
    integer, intent(in) :: nbonds, refbond, sites, Lx, Ly, ti, symm
    integer, intent(in) :: bsites(5,nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(geo%sites)
    integer(kind=8), intent(in) :: dim, basis(dim)
    double precision, intent(in) :: id, mir(6), rot(5), t, psi(dim), nnnVec(2, 3), norm(dim) 
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bcurrent(:)
    

    integer(kind=8) :: loc = 0, newst = 0, rep = 0, state = 0
    integer :: x, y, info, i, j, k, l1, l2, nd1, nd2
    integer :: sign = 1, signrep = 1, signstate = 1, signDir = 1, order = 0 
    integer :: site1, site2, refsite1, refsite2, phase, parity
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:), evals(:), gs(:)

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
                        if(ti == 1) call representative_reg(newst, sites, Lx, Ly, symm, id, par, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                        if(ti == 0) rep = newst
                        call findstate(dim, rep, basis, loc) 
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

subroutine current_tilted_dp(nHel, tilt, t, dim, sites, Lx, Ly, ti, symm, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, refbond, psi, nnnVec, norm, bcurrent, current)
    ! Calculates current-current structure factor (current) and local current current correlations (bondcurrent) on the sublattice for a given reference bond.
    ! This subroutine is used for the tilted case.
    implicit none
    integer,                       intent(in)  :: nHel, tilt, nbonds, refbond, sites, Lx, Ly, ti, symm, bsites(5,nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(geo%sites)
    integer(kind=8),               intent(in)  :: dim, basis(dim)
    double precision,              intent(in)  :: id, mir(6), rot(5), t, psi(dim), nnnVec(2, 3), norm(dim) 
    double precision,              intent(out) :: current
    double precision, allocatable, intent(out) :: bcurrent(:)
    
    integer(kind=8)                  :: loc = 0, newst = 0, rep = 0, state = 0
    integer                          :: x, y, info, i, j, k, l1, l2
    integer                          :: sign = 1, signrep = 1, signstate = 1, signDir = 1, order = 0 
    integer                          :: site1, site2, refsite1, refsite2, phase, parity
    double precision                 :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:)

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
    current  = 0.d0
    bcurrent = 0.d0 
   
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

                        if(ti == 1) call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, symm, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                        if(ti == 0) rep = newst
                        call findstate(dim, rep, basis, loc) 
                        
                        if(loc <= 0) cycle   
                        if(abs(signstate) > 1) stop 'signstate > 1' 
                        if(abs(signrep) > 1)   stop 'signrep > 1' 
                        if(abs(signDir) > 1)   stop 'signDir > 1' 
                        
                        sign = signstate * signrep * signDir * (-1)**parity 
                        psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(j)
                        
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

subroutine current_dc(tilted, nHel, tilt, k1, k2, t, dim, sites, Lx, Ly, basis, bsites, xy, xtransl, ytransl, nbonds, refbond, psi, nnnVec, norm, bondcurrent, current)

    implicit none
    integer,                       intent(in)  :: tilted, nHel, tilt, k1, k2, nbonds, refbond, sites, Lx, Ly
    integer(kind=8),               intent(in)  :: dim
    integer(kind=8),               intent(in)  :: basis(dim)
    integer,                       intent(in)  :: bsites(5,nbonds)
    integer,                       intent(in)  :: xy(4,nbonds)
    integer,                       intent(in)  :: xtransl(2, sites), ytransl(2, sites)
    double precision,              intent(in)  :: t 
    double precision,              intent(in)  :: nnnVec(2, 3)
    double complex,                intent(in)  :: psi(dim)
    double precision,              intent(in)  :: norm(dim)
    double precision,              intent(out) :: current
    double precision, allocatable, intent(out) :: bondcurrent(:)
    
    integer(kind=8)                  :: loc = 0, newst = 0, rep = 0, cntr = 0 
    integer                          :: x = 0, y = 0, i = 0, j = 0, l1 = 0, l2 = 0, l = 0, l11 = 0, l22 = 0 
    integer                          :: sign = 0, signt = 0, site1 = 0, site2 = 0, phase = 0 
    integer                          :: refsite1 = 0, refsite2 = 0, parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision                 :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double complex, allocatable      :: psiprime(:)
    character*250                    :: name, name2 

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
                                call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, symm, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signt)
                                l11 = Ly
                                l22 = Lx
                            else if(tilted == 0) then 
                                call representative_reg(newst, sites, Lx, Ly, symm, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                l11 = Lx
                                l22 = Ly
                            end if
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt(ibits(basis(j), refsite1, sites)) !Parity of c_refs1
                            parity2 = popcnt(ibits(ibclr(basis(j), refsite1-1), refsite2, sites)) !Parity of cd_refs2 
                            parity3 = popcnt(ibits(ibset( ibclr(basis(j), refsite1-1), refsite2-1), site1, sites)) !Parity of c_site1
                            parity4 = popcnt(ibits(ibclr(ibset(ibclr(basis(j), refsite1-1), refsite2-1), site1 - 1), site2, sites)) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(j) = psiprime(j) + sign * signt * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(loc) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                        else if(btest(basis(j), site2-1) .and. .not.(btest(basis(j), site1-1))) then !Forward backwards hopping from sites 2 -> 1                   
                            ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                            newst = ibclr(ibset(newst, site1-1), site2-1) !Create on site 1, annihilate on site 2
                            if(tilted == 1) then 
                                call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, symm, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signt)                                
                                l11 = Ly
                                l22 = Lx
                            else if(tilted == 0) then 
                                call representative_reg(newst, sites, Lx, Ly, symm, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.                      
                                l11 = Lx
                                l22 = Ly
                            end if
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt(ibits(basis(j), refsite1, sites)) !Parity of c_refs1
                            parity2 = popcnt(ibits(ibclr(basis(j), refsite1-1), refsite2, sites)) !Parity of cd_refs2 
                            parity3 = popcnt(ibits(ibset( ibclr(basis(j), refsite1-1), refsite2-1), site2, sites)) !Parity of c_site2
                            parity4 = popcnt(ibits(ibclr(ibset(ibclr(basis(j), refsite1-1), refsite2-1), site2-1), site1, sites)) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(j) = psiprime(j) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(loc) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                        else 
                            cycle 
                        end if     
                    else if( btest( basis(j), refsite2-1) .and. .not.( btest( basis(j), refsite1-1))) then !Backwards hopping from refsites 2 -> 1 
                        newst = ibclr(ibset(basis(j), refsite1-1), refsite2-1) !Create on refsite 1, annihilate on refsite 2
                        if(btest(basis(j), site1-1) .and. .not.(btest(basis(j), site2-1))) then !Backwards forward hopping from sites 1 -> 2       
                            ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                            newst = ibclr(ibset(newst, site2-1), site1-1) !Create on site 2, annihilate on site 1 
                            if(tilted == 1) then 
                                call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, symm, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signt)                                
                                l11 = Ly
                                l22 = Lx
                            else if(tilted == 0) then 
                                call representative_reg(newst, sites, Lx, Ly, symm, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.                   
                                l11 = Lx
                                l22 = Ly
                            end if
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle                              
                            parity1 = popcnt(ibits(basis(j), refsite2, sites)) !Parity of c_refs2
                            parity2 = popcnt(ibits(ibclr(basis(j), refsite2-1), refsite1, sites)) !Parity of cd_refs1 
                            parity3 = popcnt(ibits(ibset( ibclr(basis(j), refsite2-1), refsite1-1), site1, sites)) !Parity of c_site1
                            parity4 = popcnt(ibits(ibclr(ibset(ibclr(basis(j), refsite2-1), refsite1-1), site1 - 1), site2, sites)) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(j) = psiprime(j) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(loc) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                        else if(btest(basis(j), site2-1) .and. .not.(btest(basis(j), site1-1))) then !Backwards backwards hopping from sites 2 -> 1                   
                            ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                            newst = ibclr(ibset(newst, site1-1), site2-1) !Create on site 1, annihilate on site 2
                            if(tilted == 1) then 
                                call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, symm, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signt)                                
                                l11 = Ly
                                l22 = Lx
                            else if(tilted == 0) then 
                                call representative_reg(newst, sites, Lx, Ly, symm, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.   
                                l11 = Lx
                                l22 = Ly
                            end if
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt(ibits(basis(j), refsite2, sites)) !Parity of c_refs2
                            parity2 = popcnt(ibits(ibclr(basis(j), refsite2-1), refsite1, sites)) !Parity of cd_refs1 
                            parity3 = popcnt(ibits(ibset( ibclr(basis(j), refsite1-1), refsite2-1), site2, sites)) !Parity of c_site2
                            parity4 = popcnt(ibits(ibclr(ibset(ibclr(basis(j), refsite1-1), refsite2-1), site2-1), site1, sites)) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(j) = psiprime(j) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(loc) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
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

subroutine representative_reg(s, n, Lx, Ly, symm, id, mir, rot, xtransl, ytransl, refl, c6, r, l1, l2, sign)
    !Finds the representative 'r' for state 's'. 'n' is the number of sites and 'l' the number of shifts needed to translate 's' to 'r'.
    implicit none
    integer(kind=8), intent(in) :: s
    integer, intent(in) :: n, symm, Lx, Ly, xtransl(2, n), ytransl(2, n), refl(6, n), c6(n)
    double precision, intent(in) :: id, mir(6), rot(5)
    integer(kind=8), intent(out) :: r
    integer, intent(out) :: l1, l2, sign 

    integer :: nx = 0, ny = 0, i = 0 
    integer :: edgeA = 0, edgeB = 0, rowSt = 0, signt = 1
    integer(kind=8) :: t = 0, tx = 0, ty = 0

    l1    = Lx 
    l2    = Ly
    sign  = 1
    r     = s
    t     = 0
    tx    = 0 
    ty    = 0 
    signt = 1
    rowSt = 0 
    edgeA = 0 
    edgeB = 0 

    do nx = 1, Lx !Identity x translations 
        if(nx == 1) tx = s 
        if(nx > 1) tx = t
        call xtranslate(tx, Ly, n, Lx, xtransl, signt, t)
        if (t < r) then !New rep found
            sign = signt !Update sign 
            r    = t
            l1   = 0
            l2   = nx
        end if

        do ny = 1, Ly 
            ty = t                
            !Note, 'signt' and 't' are carried on through out all translations. 
            call ytranslate(ty, Ly, -1, n, Lx, ytransl, signt, t)
            if (t < r) then
                sign = signt 
                r    = t
                l1   = nx
                l2   = ny
            end if
        end do

    end do !Leaving this loop, r <= s 

    if(symm == 0) go to 12
    do i = 1, 6 !Check for representatives among reflected states
        call mirror_rep(r, s, n, Ly, -1, mir(i), Lx, Ly, refl(i, 1:n), xtransl, ytransl, sign, l2, l1)
    end do 
    do i = 1, 5 !Check for representatives among rotated states
        call rotate_rep(r, s, n, Ly, -1, i, rot(i), Lx, Ly, c6, xtransl, ytransl, sign, l2, l1)
    end do 
    12 continue

    return 

end subroutine representative_reg

subroutine representative_tilted(s, n, nHel, tilt, Lx, Ly, symm, id, mir, rot, xtransl, ytransl, refl, c6, r, l1, l2, sign)
    !Finds the representative 'r' for state 's'. 'n' is the number of sites and 'l' the number of shifts needed to translate 's' to 'r'.
    implicit none
    integer(kind=8), intent(in) :: s
    integer, intent(in) :: n, nHel, tilt, symm, Lx, Ly, xtransl(2, n), ytransl(2, n), refl(6, n), c6(n)
    double precision, intent(in) :: id, mir(6), rot(5)
    integer(kind=8), intent(out) :: r
    integer, intent(out) :: l1, l2, sign 

    integer :: nx = 0, ny = 0, i = 0 
    integer :: edgeA = 0, edgeB = 0, rowSt = 0, signt = 1
    integer(kind=8) :: t = 0, tx = 0, ty = 0

    l2    = Lx 
    l1    = Ly
    sign  = 1
    r     = s
    t     = 0
    tx    = 0 
    ty    = 0 
    signt = 1
    rowSt = 0 
    edgeA = 0 
    edgeB = 0 

    do nx = 1, Lx !Identity x translations 
        if(nx == 1) tx = s 
        if(nx > 1) tx = t
        call xtranslate(tx, nHel, n, Lx, xtransl, signt, t)
        if (t < r) then !New rep found
            sign = signt !Update sign 
            r    = t
            l1   = 0
            l2   = nx
        end if
        if(nHel > 1) then 
            do ny = 1, Ly 
                ty = t                
                !Note, 'signt' and 't' are carried on through out all translations. 
                call ytranslate(ty, nHel, tilt, n, Lx, ytransl, signt, t)
                if (t < r) then
                    sign = signt 
                    r    = t
                    l1   = ny
                    l2   = nx
                end if
            end do
        end if 
    end do !Leaving this loop, r <= s 

    if(symm == 0) go to 12
    do i = 1, 6 !Check for representatives among reflected states
        call mirror_rep(r, s, n, nHel, tilt, mir(i), Lx, Ly, refl(i, 1:n), xtransl, ytransl, sign, l1, l2)
    end do 
    do i = 1, 5 !Check for representatives among rotated states
        call rotate_rep(r, s, n, nHel, tilt, i, rot(i), Lx, Ly, c6, xtransl, ytransl, sign, l1, l2)
    end do 
    12 continue

    return 

end subroutine representative_tilted


subroutine mirror_rep(r, s, n, nHel, tilt, p, Lx, Ly, refl, xtransl, ytransl, sign, l1, l2)
    !Swap l1 and l2, set nHel = Ly, and tilt < 0, if cluster is not tilted, i.e. rectangular. 
    implicit none
    integer(kind=8), intent(in) :: s !original state
    integer, intent(in) :: n, nHel, tilt, Lx, Ly, refl(n), xtransl(2,n), ytransl(2,n)     
    double precision, intent(in) :: p
    integer(kind=8), intent(inout) :: r !potential rep (r<=s) 
    integer, intent(inout) :: l1, l2 
    integer, intent(inout) :: sign !If 'r' remains unchanged, previous sign remains unchanged. Otherwise updated. 

    integer :: nx = 0, ny = 0, signp = 1, signt = 1, info = 0
    integer(kind=8) :: t = 0, tx = 0, ty = 0 

    signp = 1
    signt = 1
    call reflect(s, s, n, refl, signp, info, t) ! compare to s, reflect s onto t
    if (t < r) then !compare to previously identified potenital representative 'r' 
        sign = signp * p 
        r    = t
        l1   = 0 !No shifts for k=0
        l2   = 0 !No shifts for k=0
    end if
    do nx = 1, Lx
        tx = t 
        call xtranslate(tx, nHel, n, Lx, xtransl, signt, t)
        if(nHel == 1) then 
            if (t < r) then
                sign = signt * p * signp  
                r    = t
                l1   = 0
                l2   = nx
            end if
        else if(nHel > 1) then 
            do ny = 1, Ly 
                ty = t               
                call ytranslate(ty, nHel, tilt, n, Lx, ytransl, signt, t)                
                if (t < r) then
                    sign = signt * signp * p
                    r    = t
                    l1   = ny
                    l2   = nx
                end if
            end do
        end if 
    end do

end subroutine mirror_rep

subroutine rotate_rep(r, s, n, nHel, tilt, nrot, eval, Lx, Ly, c6, xtransl, ytransl, sign, l1, l2)
    implicit none
    !Swap l1 and l2, set nHel = Ly, and tilt < 0, if cluster is not tilted, i.e. rectangular. 
    integer(kind=8), intent(in) :: s !original state
    integer, intent(in) :: n, nHel, tilt, nrot, Lx, Ly, c6(n), xtransl(2,n), ytransl(2,n)     
    double precision, intent(in) :: eval
    integer(kind=8), intent(inout) :: r !potential rep (r<=s) 
    integer, intent(inout) :: l1, l2, sign !If 'r' remains unchanged, previous sign remains unchanged. Otherwise updated. 

    integer :: nx = 0, ny = 0, signrot = 1, signt = 1, info = 0
    integer(kind=8) :: t = 0, tx = 0, ty = 0 

    signrot = 1
    signt   = 1 
    call c6n(s, s, n, nrot, c6, signrot, info, t) ! compare to s, reflect s onto t
    if (t < r) then !compare to previously identified potenital representative 'r' 
        sign = signrot * eval 
        r    = t
        l1   = 0 !No shifts for k=0
        l2   = 0 !No shifts for k=0
    end if
    do nx = 1, Lx
        tx = t 
        call xtranslate(tx, nHel, n, Lx, xtransl, signt, t)
        if(nHel == 1) then 
            if (t < r) then
                sign = signt * eval * signrot  
                r    = t
                l1   = 0
                l2   = nx
            end if
        else if(nHel > 1) then 
            do ny = 1, Ly 
                ty = t               
                call ytranslate(ty, nHel, tilt, n, Lx, ytransl, signt, t)                
                if (t < r) then
                    sign = signt * eval * signrot
                    r    = t
                    l1   = ny
                    l2   = nx
                end if
            end do
        end if 
    end do

end subroutine rotate_rep

subroutine reflect(s0, s, sites, refl, sign, info, sr)

    implicit none
  
    integer(kind=8), intent(in) :: s0, s
    integer, intent(in) :: sites 
    integer, intent(in) :: refl(geo%sites)
    integer(kind=8), intent(out) :: sr !Reflected state 
    integer, intent(out) :: sign, info 

    integer(kind=8) :: t = 0
    integer :: i = 0, flag = 0

    i = 0
    t = 0
    sr = 0 
    sign = 1 
    flag = 0 
    info = 0
    
    do i = 1, sites  !Reflect state 
        call mvbits(s, i - 1, 1, sr, refl(i) - 1)
    end do

    t = sr 
    do i = 1, sites !Sign change due to reflection 
        if(btest(t, refl(i) - 1)) then 
            ! if(i < refl(i)) sign = sign * (-1)**popcnt(ibits(t, i-1, refl(i) - 1)) 
            sign = sign * (-1)**popcnt(ibits(t, 0, refl(i) - 1)) 
            t = ibclr(t, refl(i) - 1)
        end if 
    end do 

    if(sr < s0) info = - 1 !Reflected state is already contained
    
    return  
  
end subroutine reflect

subroutine c6n(s0, s, sites, n, c6, sign, info, sr)

    implicit none
  
    integer(kind=8), intent(in) :: s0, s
    integer, intent(in) :: sites, n  
    integer, intent(in) :: c6(geo%sites)
    integer(kind=8), intent(out) :: sr !Rotated state 
    integer, intent(out) :: sign, info 

    integer(kind=8) :: t = 0, si = 0
    integer :: i = 0, j = 0

    i = 0
    t = 0
    sr = 0 
    si = s
    sign = 1 
    info = 0

    do j = 1, n !Number of subsequent C6 rotations 
        if(j > 1) si = sr !If more than one C6 rotation, start with rotated state from previous rotation 
        do i = 1, sites  !Rotate state 
            call mvbits(si, i - 1, 1, sr, c6(i) - 1)
        end do
        t = sr 
        do i = 1, sites !Sign change due anticommutations resulting from rotation 
            if(btest(t, c6(i) - 1)) then 
                sign = sign * (-1)**popcnt(ibits(t, 0, c6(i) - 1)) 
                t = ibclr(t, c6(i) - 1)
            end if 
        end do 

    end do 

    if(sr < s0) info = - 1 !Rotated state is already contained
   
    return  
  
end subroutine c6n

subroutine xtranslate(s, nHel, n, Lx, xtransl, sign, sx)
    implicit none 
    integer(kind=8), intent(in)    :: s 
    integer,         intent(in)    :: nHel, n, Lx
    integer,         intent(in)    :: xtransl(2, n)
    integer,         intent(inout) :: sign 
    integer(kind=8), intent(out)   :: sx 

    integer           :: i = 0, edgeA = 0, edgeB = 0, rowSt = 0
    integer(kind = 8) :: t = 0


    t = 0
    if(nhel == 1) then !Single helix, quasi-1D
        edgeB = n-1 !Last B site in helix
        edgeA = n-2 !Last A site in helix
        if(btest(s, edgeB) .neqv. btest(s, edgeA)) then !If only one edge occupied 
            rowSt = 0 !First site 
            sign = sign * (-1)**popcnt(ibits(s, rowst, n-2 ) ) !Occupation of rest of helix 
        end if 
    else
        do i = 1, nHel !Calculate anticommutation sign from x-translations 
            edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
            edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
            if(btest(s, edgeB) .and. btest(s, edgeA)) then !If both edges occupied, no sign
                cycle 
            else if(btest(s, edgeB) .neqv. btest(s, edgeA)) then !If only one edge occupied 
                rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
                sign = sign * (-1)**popcnt(ibits(s, rowst, 2*Lx-2)) !Occupation of y-layer 'i'
            end if 
        end do 
    end if 

    do i = 1, n
        call mvbits(s, xtransl(1,i)-1, 1, t, xtransl(2,i)-1)      
    end do
    sx = t 
    

end subroutine xtranslate

subroutine ytranslate(s, nHel, tilt, n, Lx, ytransl, sign, sy)
    !Set tilt < 0 if cluster is not tilted, i.e. rectangular. 
    implicit none 
    integer(kind=8), intent(in)    :: s 
    integer,         intent(in)    :: nHel, tilt, n, Lx
    integer,         intent(in)    :: ytransl(2, n)
    integer,         intent(inout) :: sign 
    integer(kind=8), intent(out)   :: sy

    integer           :: i = 0, ntot = 0, edgeA = 0, edgeB = 0, rowSt = 0
    integer(kind = 8) :: t = 0
    
    t = 0
    ntot = popcnt(s)
    if(nhel == 1) then !Single helix
        sign = sign * ((-1)**(ntot-popcnt(ibits(s, n - tilt, tilt))))**popcnt(ibits(s, n - tilt, tilt)) 
    else 
        sign = sign * ((-1)**(ntot-popcnt(ibits(s, 2*Lx*(nHel-1), 2*Lx))))**popcnt(ibits(s, 2*Lx*(nHel-1), 2*Lx)) 
        if(tilt < 0) goto 11
        do i = 1, nHel
            edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
            edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
            if(btest(s, edgeB) .and. btest(s, edgeA)) then !If both edges occupied, no sign
                cycle 
            else if(btest(s, edgeB) .neqv. btest(s, edgeA)) then !If only one edge occupied 
                rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
                sign = sign * (-1)**popcnt(ibits(s, rowst, 2*Lx-2)) !Occupation of y-layer 'i'
            end if 
        end do     
        11 continue 
    end if 

    do i = 1, n
        call mvbits(s, ytransl(1,i)-1, 1, t, ytransl(2,i)-1)
    end do
    sy = t

end subroutine ytranslate

end module observables