module basis ! Module for defining the basis for the extended Hubbard model on the honeycomb lattice.
    use parameters
    use functions
    
    implicit none
    
    contains 

    subroutine make_basis(ti, tilted, pat, nnnVec, sites, particles, dim, symmetrize, ucx, ucy, l1, l2, basis_states, abasis, bbasis, tilt, nHel, k1, k2, xtransl, ytransl, id, par, rot, refl, c6, period, norm, orbsize, orbits2D, phases2D, norm2D)

        implicit none

        integer, intent(in) :: ti, tilted, sites, particles, tilt, nHel, ucx, ucy, k1, k2, symmetrize
        integer, allocatable, intent(in) :: xtransl(:,:), ytransl(:,:)
        double precision, intent(in) :: nnnVec(2,3), id, par(6), rot(5)
        character(len=*), intent(in) :: pat

        integer, intent(out) :: l1, l2, orbsize
        integer(kind=8), intent(out) :: dim
        integer(kind=8), allocatable, intent(out) :: basis_states(:), abasis(:), bbasis(:),period(:), orbits2D(:,:,:)
        integer, allocatable, intent(out) :: refl(:,:), c6(:)
        double precision, allocatable, intent(out) :: norm(:), norm2D(:,:)
        double complex, allocatable, intent(out) :: phases2D(:,:,:)
        
        integer :: a(sites), temp = 0, I_in = 0, I_tail = 0, i, l
        integer(kind=8) :: symdim, amask, bmask, j
        integer(kind=8), allocatable :: symbasis(:)    

        dim = int(fact(sites) / (fact(particles) * fact(max(1,sites-particles))),8) !Number of basis states


        if (dim == 0) then
            print*, 'No basis states available.'
            return 
        end if

        !Permutations contains integer values I of basis states, perm_up/dn contain the integer values I_up/dn and the corresponDing J_up/dn
        if (allocated(basis_states)) deallocate(basis_states)
        allocate(basis_states(dim))
        if (allocated(abasis)) deallocate(abasis)
        allocate(abasis(dim))
        if (allocated(bbasis)) deallocate(bbasis)
        allocate(bbasis(dim))
        basis_states  = 0
        abasis = 0
        bbasis = 0
        amask  = 0 
        bmask  = 0 
        if(pat=='AB') then 
            do i = 0, sites-2, 2
                amask = ibset(amask, i)
                bmask = ibset(bmask, i + 1)
            end do 
        else if(pat=='BA') then 
            do i = 0, sites-2, 2
                amask = ibset(amask, i + 1)
                bmask = ibset(bmask, i)
            end do 
        end if 

        a(1 : sites-particles) = 0 !'a' contains the initial basis state with all '1's to the right
        a(sites-particles+1 : sites) = 1
        I_in = 0
        do l = 1, sites !Binary representation of 'a'
            I_in = I_in + a(l) * 2**(sites-l)
        end do

        do j = 1, dim !Generates all possible configurations in 'basis'
            I_tail = 0
            basis_states(j)  = I_in
            abasis(j) = iand( basis_states(j), amask ) 
            bbasis(j) = iand( basis_states(j), bmask )
            temp = I_in
            do i = 0, 64
                if (btest(temp,i)) then
                    temp = ibclr(temp,i)
                    if (.not. btest(temp,i+1)) then   ! asks if pos i+1 is zero
                        temp = ibset(temp,i+1)
                        exit
                    end if
                I_tail = ibset(ishft(I_tail,1),0) ! only generated if first loop finds no '01' pair and has to be repeated
                end if
            end do
            I_in = temp + I_tail
        end do

        print*, 'Basis generated.'
        if(ti == 0) return 
        if(allocated(refl)) deallocate(refl)
        allocate(refl(6, sites))
        if(allocated(c6)) deallocate(c6)
        allocate(c6(sites))
        refl = 0
        c6   = 0    
        if(tilted == 1) then !18A
            refl(1, 1:sites) = [16,17,18,13,14,15,10,11,12,7,8,9,4,5,6,1,2,3] !Mirror axis through edges
            refl(2, 1:sites) = [8,7,4,3,18,17,2,1,16,15,12,11,14,13,10,9,6,5] !Mirror axis through edges
            refl(3, 1:sites) = [12,15,14,5,4,7,6,9,8,17,16,1,18,3,2,11,10,13] !Mirror axis through edges
            refl(4, 1:sites) = [11,2,3,18,13,10,17,8,9,6,1,16,5,14,15,12,7,4] !Mirror axis through sites
            refl(5, 1:sites) = [1,6,5,4,3,2,7,12,11,10,9,8,13,18,17,16,15,14] !Mirror axis through sites
            refl(6, 1:sites) = [9,10,13,14,5,6,15,16,1,2,11,12,3,4,7,8,17,18] !Mirror axis through sites
            c6 = [8,17,18,3,4,7,2,11,12,15,16,1,14,5,6,9,10,13] !C6 rotation
        else if(tilted == 0) then !Rectangular L=18
            !Odd index = axis through edges, even index = axis through sites
            !Index - 1 = power of C6 rotations
            refl(1, 1:sites) = [2,1,8,7,14,13,4,3,10,9,16,15,6,5,12,11,18,17] !Mirror axis through edges
            refl(2, 1:sites) = [3,2,1,6,5,4,11,10,9,8,7,12,13,18,17,16,15,14] !Mirror axis through sites
            refl(3, 1:sites) = [10,3,2,13,18,11,16,9,8,1,6,17,4,15,14,7,12,5] !Mirror axis through edges
            refl(4, 1:sites) = [9,10,3,4,15,16,7,8,1,2,13,14,11,12,5,6,17,18] !Mirror axis through sites
            refl(5, 1:sites) = [8,9,10,11,12,7,6,1,2,3,4,5,16,17,18,13,14,15] !Mirror axis through edges
            refl(6, 1:sites) = [1,8,9,16,17,6,13,2,3,10,11,18,7,14,15,4,5,12] !Mirror axis through sites
            c6               = [2,3,10,11,18,13,6,1,8,9,16,17,4,5,12,7,14,15] !C6 rotation
        end if 

        if(tilted == 0) then 
            l2 = ucx 
            l1 = ucy 
        else 
            if(nHel == 1) then
                    l1 = 1
            else if(nHel > 1) then    
                if(modulo(dble(sites)/dble((nHel*tilt)), dble(nhel)) == 0.d0) then 
                    l1 = sites/(nHel * tilt) 
                else if(modulo(dble(sites)/dble((nHel*tilt)), dble(nhel)) >= 1.d0) then 
                    l1 = sites/tilt 
                else if(modulo(dble(sites)/dble((nHel*tilt)), dble(nhel)) < 1.d0) then 
                    l1 = ceiling(sites/(nHel * tilt * modulo(dble(sites)/dble((nHel*tilt)), 1.d0)))
                end if 
            end if 
            l2 = sites/(2*nHel) 
        end if 
        call symm_basis(tilted, dim, sites, nHel, l2, l1, k2, k1, symmetrize, id, par, rot, nnnVec, basis_states, xtransl, ytransl, refl, c6, symdim, symbasis, period, norm, orbsize, orbits2D, phases2D, norm2D)
        
        dim = symdim
        deallocate(basis_states)
        allocate(basis_states(dim))
        basis_states = symbasis 
        print*, 'Basis symmetrized.'

    end subroutine make_basis

    !Rotations, reflections and translations
    subroutine symm_basis(tilted, dim, n, nHel, Lx, Ly, k1, k2, irrep, id, par, rot, nnnVec, basis, xtransl, ytransl, refl, c6, symdim, symbasis, period, norm, orbsize, orbits2D, phases2D, norm2D)
        
        implicit none
        integer(kind=8), intent(in) :: dim, basis(dim)
        integer, intent(in) :: tilted, n, nHel, k1, k2, irrep, Lx, Ly
        integer, intent(in) :: xtransl(2, n), ytransl(2, n), refl(6, n), c6(n)
        double precision, intent(in) :: nnnVec(2,3), id, par(6), rot(5)
        integer, intent(out) :: orbsize
        integer(kind=8), intent(out) :: symdim
        integer(kind=8), allocatable, intent(out) :: symbasis(:), orbits2D(:,:,:), period(:)
        double precision, allocatable, intent(out) :: norm(:), norm2D(:,:)
        double complex, allocatable, intent(out) :: phases2D(:,:,:)

        integer :: r = 0, temp(n), kx, ky, maxorb
        integer(kind=8) :: j = 0, cntr = 0, cntr2 = 0 
        integer(kind=8), allocatable :: symbasis_temp(:), period_temp(:), orbits2D_temp(:,:,:), orbarr(:,:)
        double precision, allocatable :: norm_temp(:), norm2D_temp(:,:)
        double precision :: normalization = 0.d0
        double complex, allocatable :: phases2D_temp(:,:,:)
        
        if(tilted == 1) then 
            kx = k2
            ky = k1 
        else 
            kx = k1 
            ky = k2 
        end if 
        temp    = 1
        
        !Define temporary arrays due to unknown final Hilbert space dimension of symmetrized block
        if(allocated(period_temp)) deallocate(period_temp)
        
    
        if(id == 1) then     
            if(allocated(norm_temp)) deallocate(norm_temp)
            allocate(norm_temp(dim))
            norm_temp = 0 
        else if(id == 2) then
            orbsize = size(par) * size(rot) * Lx * Ly
            if(allocated(norm2D_temp)) deallocate(norm2D_temp)
            if(allocated(orbits2D_temp)) deallocate(orbits2D_temp)
            if(allocated(phases2D_temp)) deallocate(phases2D_temp)
            allocate(norm2D_temp(dim, 2))
            allocate(orbits2D_temp(dim, orbsize, 2))
            allocate(orbarr(orbsize, 2))
            allocate(phases2D_temp(dim, orbsize, 2))
            norm2D_temp   = 0.d0 
            phases2D_temp = 0.d0
            orbits2D_temp = 0 
            orbarr = 0
        end if
        allocate(symbasis_temp(dim))
        allocate(period_temp(dim))
        normalization = 0.d0 
        symbasis_temp = 0
        period_temp   = 0
        cntr          = 0 
        cntr2         = 0 
        symdim        = 0 
        r             = 0 
        

        do j = 1, dim
            if(id == 1) then 
                if(tilted == 0) then 
                    call checkstate_rect(basis(j), n, Lx, Ly, kx, ky, irrep, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, r, normalization)
                else if(tilted == 1) then 
                    call checkstate(basis(j), n, nHel, Lx, Ly, kx, ky, irrep, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, r, normalization)
                end if 
                if(r >= 0) then !New representative state found 
                    cntr                = cntr + 1
                    symbasis_temp(cntr) = basis(j)
                    period_temp(cntr)   = r       
                    norm_temp(cntr)     = normalization   
                end if
            else if(id == 2) then 
                ! if(r >= 0) then !New representative state found 
                !     cntr                 = cntr + 1
                !     symbasis_temp(cntr)  = basis(j)
                !     period_temp(cntr)    = r       
                !     if(id == 1) norm_temp(cntr) = normalization   
                !     if(id == 2) then                    
                
                ! allocate(orbits2D(r, 2), phases2D(r, 2))
                cntr = cntr + 1 !Number of representatives (tentative)
                call checkstate2D(basis(j), orbsize, n, tilted, nHel, Lx, Ly, kx, ky, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, orbits2D_temp(cntr, 1:orbsize, 1:2), norm2D_temp(cntr, 1:2), phases2D_temp(cntr, 1:orbsize, 1:2), r)



                if(r >= 0) then !New representative state found 
                    ! cntr                 = cntr + 1 
                    symbasis_temp(cntr)  = basis(j) !Representative
                    period_temp(cntr)    = r        !Orbit size 
                    if(cntr == 1) maxorb = r        !
                    if(r > maxorb) maxorb = r       !Size of largest orbit for array allocation (later)
                else if(r < 0) then 
                    cntr = cntr - 1 !Reset counter 
                end if
                
            ! end if 

            end if


        end do

        symdim = cntr
        
        if(allocated(period)) deallocate(period)
        
        if(id == 1) then 
            if(allocated(norm)) deallocate(norm)
            allocate(norm(symdim))
            norm = 0.d0  
            norm(1:symdim) = norm_temp(1:symdim)
            if(allocated(norm_temp)) deallocate(norm_temp)
        else if(id == 2) then 
            if(allocated(norm2D)) deallocate(norm2D)
            if(allocated(orbits2D)) deallocate(orbits2D)
            if(allocated(phases2D)) deallocate(phases2D)
            allocate(norm2D(symdim, 2))
            allocate(orbits2D(symdim, maxorb, 2))
            allocate(phases2D(symdim, maxorb, 2))
            
            norm2D   = 0.d0 
            phases2D = 0.d0 
            orbits2D = 0
            norm2D(1:symdim, 1) = norm2D_temp(1:symdim, 1)
            norm2D(1:symdim, 2) = norm2D_temp(1:symdim, 2)
            ! do i = 1, orbsize            
            orbits2D(1:symdim, 1:maxorb, 1) = orbits2D_temp(1:symdim, 1:maxorb, 1)
            orbits2D(1:symdim, 1:maxorb, 2) = orbits2D_temp(1:symdim, 1:maxorb, 2)
            phases2D(1:symdim, 1:maxorb, 1) = phases2D_temp(1:symdim, 1:maxorb, 1)
            phases2D(1:symdim, 1:maxorb, 2) = phases2D_temp(1:symdim, 1:maxorb, 2)
            
            orbsize = maxorb 
            if(allocated(norm2D_temp)) deallocate(norm2D_temp)
            if(allocated(orbits2D_temp)) deallocate(orbits2D_temp)
            if(allocated(phases2D_temp)) deallocate(phases2D_temp)
        end if 
        if(allocated(symbasis)) deallocate(symbasis)
        allocate(symbasis(symdim), period(symdim))
        symbasis           = 0
        period             = 0 
        symbasis(1:symdim) = symbasis_temp(1:symdim)
        period(1:symdim)   = period_temp(1:symdim)
        
        if(allocated(symbasis_temp)) deallocate(symbasis_temp)
        if(allocated(period_temp))   deallocate(period_temp)

    end subroutine symm_basis

    !Checkstate for tilted lattices: For rotations, reflections and translations (single sum of phases for all orbits.)
    subroutine checkstate(s, sites, nHel, Lx, Ly, kx, ky, irrep, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, r, norm)

        implicit none
    
        ! Given momentum
        integer(kind=8), intent(in) :: s
        integer, intent(in) :: kx, ky
        integer, intent(in) :: sites, nHel, irrep
        integer, intent(in) :: xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)

        ! Define lattice parameters
        integer, intent(in) :: Lx   ! Number of unit cells in the a2-direction
        integer, intent(in) :: Ly   ! Number of unit cells in the a1-direction
        double precision, intent(in) :: nnnVec(2,3), id, par(6), rot(5)
        integer, intent(out) :: r
        double precision, intent(out) :: norm

        
        integer(kind=8) :: sr = 0, s0 = 0
        integer :: ntot = 0, info = 0, orbsize = 0
        integer :: sign = 1, orbit = 1, i = 0, n = 0, flag = 0, rt = 0
        integer, allocatable :: orbits(:)
        double precision :: k1(2), k2(2), k(2)
        double precision :: a1(2), a2(2)
        double precision, parameter :: tolerance = 1.0e-8
        double complex :: phase = 0.d0

        ntot = popcnt(s) !Number of particles 
        a1 = nnnVec(1:2, 1)
        a2 = nnnVec(1:2, 2)
        k1 = (/(-2.d0 * pi) / 3.d0, (-2.d0 * pi) / sqrt(3.d0) /)
        k2 = (/(-4.d0 * pi) / 3.d0, 0.d0 /)
        if(nHel == 1) then 
            k  = (dble(kx)/dble(Lx)) * k2
        else if(nHel > 1) then 
            k  = (dble(kx)/dble(Lx)) * k2 + (dble(ky)/dble(Ly)) * k1 
        end if

        r  = -1 !r is the orbit size. If r = -1, the state is not compatible with the momentum k or point group symmetry.
        rt = 0 
        s0 = s
        flag  = 0 
        norm  = 0.d0 
        orbit = 0 
        phase = 0.d0 
        orbsize = size(par) * size(rot) * Lx * Ly 

        if(allocated(orbits)) deallocate(orbits)
        allocate(orbits(orbsize))
        
        orbits = 0 
        ! orbits(1) = s
        sign = 1 
        call translation(s0, s, sites, ntot, orbit, orbits, 1, nHel, Lx, Ly, id, sign, a1, a2, xtransl, ytransl, k, phase, info)
        if(info < 0) return 

        if(irrep == 0) goto 11
        !Reflections 
        do i = 1, 6
            sign = 1   
            call reflect(s0, s, sites, refl(i,1:sites), sign, info, sr) 
            if(info < 0) return        
            call translation(s0, sr, sites, ntot, orbit, orbits, 1, nHel, Lx, Ly, par(i), sign, a1, a2, xtransl, ytransl, k, phase, info)
            if(info < 0) return 
        end do 

        !Rotations
        do n = 1, 5  
            sign = 1
            call c6n(s0, s, sites, n, c6, sign, info, sr)
            if(info < 0) return 
            call translation(s0, sr, sites, ntot, orbit, orbits, 1, nHel, Lx, Ly, rot(n), sign, a1, a2, xtransl, ytransl, k, phase, info)
            if(info < 0) return 
        end do 

        11 continue 
       
        if((abs(dble(phase)) < tolerance) .and. (abs(aimag(phase)) < tolerance)) return  

        do i = 1, size(orbits)
            if(orbits(i) > 0) rt = rt + 1 
        end do 
        r    = rt 
        norm = r * abs(phase)**2
        return


   
    end subroutine checkstate

    subroutine checkstate_rect(s, sites, Lx, Ly, kx, ky, irrep, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, r, norm)

        implicit none
    
        ! Given momentum
        integer(kind=8), intent(in) :: s
        integer, intent(in) :: sites, kx, ky, irrep
        integer, intent(in) :: xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
        ! Define lattice parameters
        integer, intent(in) :: Lx   ! Number of unit cells in the x-direction
        integer, intent(in) :: Ly   ! Number of unit cells in the y-direction
        double precision, intent(in) :: nnnVec(2,3), id, par(6), rot(5)
        integer, intent(out) :: r
        double precision, intent(out) :: norm
        
        integer(kind=8) :: sr = 0, s0 = 0
        integer :: ntot = 0, orbsize = 0, info = 0
        integer :: orbit = 1, sign = 1, i = 0, n = 0, flag = 0, rt = 0
        integer, allocatable :: orbits(:)
        double precision :: k1(2), k2(2), k(2)
        double precision :: a1(2), a2(2)
        double precision, parameter :: tolerance = 1.0e-8
        double complex :: phase = 0.d0
        
        ntot = popcnt(s)
        ! a1 = (/sqrt(3.d0), 0.d0/)
        ! a2 = 0.5*(/-1*sqrt(3.d0), 3.d0/)
        a1 = nnnVec(1:2, 1)
        a2 = nnnVec(1:2, 2)
        k1 = (/(2.0d0 * pi)/sqrt(3.d0), (2.d0 * pi)/3.d0/)
        k2 = (/0.d0, (4.d0 * pi)/3.d0/)
        k  = (dble(kx)/dble(Lx)) * k1 + (dble(ky)/dble(Ly)) * k2 
        
        flag    = 0 
        orbit   = 0
        norm    = 0.d0 
        phase   = 0.d0 
        r       = - 1 
        rt      = 0 
        s0      = s
        orbsize = size(par) * size(rot) * Lx * Ly 
        if(allocated(orbits)) deallocate(orbits)
        allocate(orbits(orbsize))    
        orbits = 0 
        sign   = 1 
        !Identity x translations
        call translation(s0, s, sites, ntot, orbit, orbits, 0, Ly, Lx, Ly, id, sign, a2, a1, xtransl, ytransl, k, phase, info)
        if(info < 0) return 
        
        if(irrep == 0) goto 11
        !Reflections 
        do i = 1, 6
            sign = 1   
            call reflect(s0, s, sites, refl(i,1:sites), sign, info, sr) 
            if(info < 0) return        
            call translation(s0, sr, sites, ntot, orbit, orbits, 0, Ly, Lx, Ly, par(i), sign, a2, a1, xtransl, ytransl, k, phase, info)
            if(info < 0) return 
        end do 

        !Rotations  
        do n = 1, 5  
            sign = 1
            call c6n(s0, s, sites, n, c6, sign, info, sr)
            if(info < 0) return 
            call translation(s0, sr, sites, ntot, orbit, orbits, 0, Ly, Lx, Ly, rot(n), sign, a2, a1, xtransl, ytransl, k, phase, info)
            if(info < 0) return 
        end do 

        11 continue 
        
        if((abs(dble(phase)) < tolerance) .and. (abs(aimag(phase)) < tolerance)) return  

        do i = 1, size(orbits)
            if(orbits(i) > 0) rt = rt + 1 
        end do 
        r = rt 
        norm = r * abs(phase)**2

        return

  
    end subroutine checkstate_rect

    subroutine checkstate2D(s, orbsize, sites, tilted, nHel, Lx, Ly, kx, ky, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, orbits, norm, phases, r)

        implicit none
    
        ! Given momentum
        integer(kind=8), intent(in) :: s
        integer, intent(in) :: orbsize, kx, ky, sites, tilted, nHel, xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
        ! Define lattice parameters
        integer, intent(in) :: Lx   ! Number of unit cells in the a2-direction
        integer, intent(in) :: Ly   ! Number of unit cells in the a1-direction
        double precision, intent(in) :: nnnVec(2,3), id, par(6), rot(5)
        
        integer, intent(out) :: r
        integer(kind=8), intent(out) :: orbits(orbsize,2)
        double precision, intent(out) :: norm(2)
        double complex, intent(out) :: phases(orbsize, 2)

        integer(kind=8) :: sr = 0, s0 = 0
        integer :: rep, sign = 1, ntot, info, layers, flag
        integer :: orbit = 1, i, n, rt, c
        ! integer, allocatable :: orbits(:)
        double precision :: k1(2), k2(2), k(2)
        double precision :: a1(2), a2(2), dcntr
        double precision, parameter :: tolerance = 1.0e-8
        
        double precision :: sigma(2,2), rho(2,2), reflection(2,2), rotation(2,2), identity(2,2) 
        double complex :: phase = 0.d0 

        if(rot(1) == 1)  rep = 1 !IRREP E1
        if(rot(1) == -1) rep = 2 !IRREP E2
        rho(1,1)   =  cos(2*pi*rep/6) !Rotation matrix entry (1,1)
        rho(1,2)   = -sin(2*pi*rep/6) !Rotation matrix entry (1,2)
        rho(2,1)   =  sin(2*pi*rep/6) !Rotation matrix entry (2,1)
        rho(2,2)   =  cos(2*pi*rep/6) !Rotation matrix entry (2,2)
        sigma(1,1) =  1               !Reflection matrix entry (1,1)
        sigma(1,2) =  0               !Reflection matrix entry (1,2)
        sigma(2,1) =  0               !Reflection matrix entry (2,1)
        sigma(2,2) = -1               !Reflection matrix entry (2,2)
        identity = 0.d0 
        identity(1, 1) = 1.d0
        identity(2, 2) = 1.d0
        if(tilted == 1) then 
            a1     = nnnVec(1:2, 1)
            a2     = nnnVec(1:2, 2)
            layers = nHel 
            k1 = (/(-2.d0 * pi) / 3.d0, (-2.d0 * pi) / sqrt(3.d0) /)
            k2 = (/(-4.d0 * pi) / 3.d0, 0.d0 /)
            if(nHel == 1) then 
                k  = (dble(kx)/dble(Lx)) * k2
            else if(nHel > 1) then 
                k  = (dble(kx)/dble(Lx)) * k2 + (dble(ky)/dble(Ly)) * k1 
            end if   
        else 
            layers = Ly
            a1     = nnnVec(1:2, 2)
            a2     = nnnVec(1:2, 1)
            k1     = (/(2.0d0 * pi)/sqrt(3.d0), (2.d0 * pi)/3.d0/)
            k2     = (/0.d0, (4.d0 * pi)/3.d0/)
            k      = (dble(kx)/dble(Lx)) * k1 + (dble(ky)/dble(Ly)) * k2 
        end if 
        ntot    = popcnt(s) !Number of particles      
        r       = -1 
        s0      = s
        norm    = 0.d0 
        flag    = 0
        ! if(allocated(orbits)) deallocate(orbits)
        ! allocate(orbits(orbsize,2))
        orbits = 0 
        phases = 0.d0
        do c = 1, 2 !IRREP basis states
            orbit  = 0            
            sign   = 1 
            rt     = 0 
            phase  = 0.d0 
             
            call translation2D(s0, s, sites, ntot, orbsize, orbit, orbits(1:orbsize, c), tilted, layers, Lx, Ly, 1.0d0, sign, a1, a2, xtransl, ytransl, k, phases(1:orbsize, c), phase, info)

            if(info < 0) phases(1:orbsize, c) = 0.d0 !All prefactors = 0.  
            if(info < 0) orbits(1:orbsize, c) = 0    !Remove orbit. 
            if(info < 0) cycle
            
            !Reflections 
            reflection = 0.d0 
            do i = 1, 6
                sign = 1
                reflection = identity 
                reflection(1,1) = cos(2*pi*(i-1)/3.0d0)
                reflection(2,2) = -cos(2*pi*(i-1)/3.0d0)
 
                call reflect(s0, s, sites, refl(i,1:sites), sign, info, sr) 
                if(info < 0) then 
                    r = -2        
                    exit 
                end if
                call translation2D(s0, sr, sites, ntot, orbsize, orbit, orbits(1:orbsize, c), tilted, layers, Lx, Ly, reflection(c,c), sign, a1, a2, xtransl, ytransl, k, phases(1:orbsize, c), phase, info)

                if(info < 0) then 
                    r = -2        
                    exit 
                end if

            end do 
            if(r == -2) then 
                r = -1  
                phases(1:orbsize, c) = 0.d0 !All prefactors = 0. 
                orbits(1:orbsize, c) = 0    !Remove orbit. 
                cycle 
            end if 

            rotation = rho 
            !Rotations
            do n = 1, 5  
                sign = 1
                call c6n(s0, s, sites, n, c6, sign, info, sr)
                if(info < 0) then 
                    r = -2        
                    exit 
                end if
                call translation2D(s0, sr, sites, ntot, orbsize, orbit, orbits(1:orbsize,c), tilted, layers, Lx, Ly, rotation(c,c), sign, a1, a2, xtransl, ytransl, k, phases(1:orbsize, c), phase, info)
                
                if(info < 0) then 
                    r = -2        
                    exit 
                end if
                rotation = matmul(rotation, rho)
                
            end do 

            if(r == -2) then 
                r = -1 
                phases(1:orbsize, c) = 0.d0 !All prefactors = 0. 
                orbits(1:orbsize, c) = 0    !Remove orbit. 
                cycle 
            end if 

            if((abs(dble(phase)) < tolerance) .and. (abs(aimag(phase)) < tolerance)) then !Basis state c has norm = 0, i.e. does not contribute to the symmetry state. Set all prefactors = 0. 
                phases(1:orbsize, c) = 0.d0 !All prefactors = 0. 
                orbits(1:orbsize, c) = 0    !Remove orbit. 
                cycle 
            end if 
            
            if(r == -2) cycle

                
            do i = 1, orbsize
                if(orbits(i,c) > 0) rt = rt + 1 
            end do 
            r       = rt 
            norm(c) = rt * abs(phase)**2 
            dcntr = 0.d0 
            do i = 1, orbsize
                dcntr = dcntr + abs(phases(i, c))**2
            end do 
        end do 



        if(r == -2) r = -1 
       
        return


   
    end subroutine checkstate2D

end module basis
