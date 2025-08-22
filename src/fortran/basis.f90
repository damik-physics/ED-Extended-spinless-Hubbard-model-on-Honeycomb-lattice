module basis 
    ! Module for constructing and managing many-body basis states including symmetry operations,
    ! momentum basis construction, and irreducible representation handling for the Hubbard model
    use types
    use params
    use functions
    use symmetries
    
    implicit none
    
    contains 

    subroutine make_basis(par, geo, st)
        ! Main basis construction routine that generates all many-body states and applies symmetry constraints
        
        implicit none
        type(sim_params),     intent(inout) :: par
        type(geometry),       intent(inout) :: geo
        type(system_state),   intent(inout) :: st
            
        integer                             :: a(geo%sites), temp, I_in, I_tail, i, l
        integer(kind=8)                     :: symdim, amask, bmask, j
        integer(kind=8), allocatable        :: momBasis(:)    

        geo%dim = int(fact(geo%sites) /(fact(geo%particles) * fact(max(1,geo%sites-geo%particles))),8) !Number of basis states


        if(geo%dim == 0) then
            print*, 'No basis states available.'
            return 
        end if

        !Permutations contains integer values I of basis states, perm_up/dn contain the integer values I_up/dn and the corresponDing J_up/dn
        if(allocated(geo%basis_states)) deallocate(geo%basis_states)
        allocate(geo%basis_states(geo%dim))
        if(allocated(geo%abasis)) deallocate(geo%abasis)
        allocate(geo%abasis(geo%dim))
        if(allocated(geo%bbasis)) deallocate(geo%bbasis)
        allocate(geo%bbasis(geo%dim))
        geo%basis_states  = 0
        geo%abasis = 0
        geo%bbasis = 0
        amask  = 0 
        bmask  = 0 
        if(pattern == 'AB') then 
            do i = 0, geo%sites-2, 2
                amask = ibset(amask, i)
                bmask = ibset(bmask, i + 1)
            end do 
        else if(pattern == 'BA') then 
            do i = 0, geo%sites-2, 2
                amask = ibset(amask, i + 1)
                bmask = ibset(bmask, i)
            end do 
        end if 

        a(1 : geo%sites-geo%particles) = 0 !'a' contains the initial basis state with all '1's to the right
        a(geo%sites-geo%particles+1 : geo%sites) = 1
        I_in = 0
        do l = 1, geo%sites !Binary representation of 'a'
            I_in = I_in + a(l) * 2**(geo%sites-l)
        end do

        do j = 1, geo%dim !Generates all possible configurations in 'basis'
            I_tail = 0
            geo%basis_states(j)  = I_in
            geo%abasis(j) = iand(int(geo%basis_states(j),8), amask) 
            geo%bbasis(j) = iand(int(geo%basis_states(j),8), bmask)
            temp = I_in
            do i = 0, 64
                if(btest(temp,i)) then
                    temp = ibclr(temp,i)
                    if(.not. btest(temp,i+1)) then   ! asks if pos i+1 is zero
                        temp = ibset(temp,i+1)
                        exit
                    end if
                I_tail = ibset(ishft(I_tail,1),0) ! only generated if first loop finds no '01' pair and has to be repeated
                end if
            end do
            I_in = temp + I_tail
        end do

        if(allocated(geo%refl)) deallocate(geo%refl)
        allocate(geo%refl(6, geo%sites))
        if(allocated(geo%c6)) deallocate(geo%c6)
        allocate(geo%c6(geo%sites))
        geo%refl = 0
        geo%c6   = 0    
        print*, 'Basis generated.'
        if(par%ti == 0) return 

        if(par%tilted == 1) then !18A
            geo%refl(1, 1:geo%sites) = [16,17,18,13,14,15,10,11,12,7,8,9,4,5,6,1,2,3] !Mirror axis through edges
            geo%refl(2, 1:geo%sites) = [8,7,4,3,18,17,2,1,16,15,12,11,14,13,10,9,6,5] !Mirror axis through edges
            geo%refl(3, 1:geo%sites) = [12,15,14,5,4,7,6,9,8,17,16,1,18,3,2,11,10,13] !Mirror axis through edges
            geo%refl(4, 1:geo%sites) = [11,2,3,18,13,10,17,8,9,6,1,16,5,14,15,12,7,4] !Mirror axis through geo%sites
            geo%refl(5, 1:geo%sites) = [1,6,5,4,3,2,7,12,11,10,9,8,13,18,17,16,15,14] !Mirror axis through geo%sites
            geo%refl(6, 1:geo%sites) = [9,10,13,14,5,6,15,16,1,2,11,12,3,4,7,8,17,18] !Mirror axis through geo%sites
            geo%c6 = [8,17,18,3,4,7,2,11,12,15,16,1,14,5,6,9,10,13] !C6 rotation
        else if(par%tilted == 0) then !Rectangular L=18
            !Odd index = axis through edges, even index = axis through geo%sites
            !Index - 1 = power of C6 rotations
            geo%refl(1, 1:geo%sites) = [2,1,8,7,14,13,4,3,10,9,16,15,6,5,12,11,18,17] !Mirror axis through edges
            geo%refl(2, 1:geo%sites) = [3,2,1,6,5,4,11,10,9,8,7,12,13,18,17,16,15,14] !Mirror axis through geo%sites
            geo%refl(3, 1:geo%sites) = [10,3,2,13,18,11,16,9,8,1,6,17,4,15,14,7,12,5] !Mirror axis through edges
            geo%refl(4, 1:geo%sites) = [9,10,3,4,15,16,7,8,1,2,13,14,11,12,5,6,17,18] !Mirror axis through geo%sites
            geo%refl(5, 1:geo%sites) = [8,9,10,11,12,7,6,1,2,3,4,5,16,17,18,13,14,15] !Mirror axis through edges
            geo%refl(6, 1:geo%sites) = [1,8,9,16,17,6,13,2,3,10,11,18,7,14,15,4,5,12] !Mirror axis through geo%sites
            geo%c6               = [2,3,10,11,18,13,6,1,8,9,16,17,4,5,12,7,14,15] !C6 rotation
        end if 


        if(par%tilted == 0) then 
            geo%l2 = par%ucx 
            geo%l1 = par%ucy 
        else 
            if(geo%nHel == 1) then
                    geo%l1 = 1
            else if(geo%nHel > 1) then    
                if(modulo(dble(geo%sites)/dble((geo%nHel* geo%tilt)), dble(geo%nhel)) == 0.d0) then 
                    geo%l1 = geo%sites/(geo%nHel * geo%tilt) 
                else if(modulo(dble(geo%sites)/dble((geo%nHel* geo%tilt)), dble(geo%nhel)) >= 1.d0) then 
                    geo%l1 = geo%sites/geo%tilt 
                else if(modulo(dble(geo%sites)/dble((geo%nHel* geo%tilt)), dble(geo%nhel)) < 1.d0) then 
                    geo%l1 = ceiling(geo%sites/(geo%nHel * geo%tilt * modulo(dble(geo%sites)/dble((geo%nHel* geo%tilt)), 1.d0)))
                end if 
            end if 
            geo%l2 = geo%sites/(2*geo%nHel) 
        end if 

        call symm_basis(par, geo, st)
        
        geo%dim = symdim
        deallocate(geo%basis_states)
        allocate(geo%basis_states(geo%dim))
        geo%basis_states = momBasis 

    end subroutine make_basis

    subroutine symm_basis(par, geo, st)
        ! Constructs symmetry-adapted basis by finding representative states under group operations
        
        implicit none

        type(sim_params),   intent(inout) :: par
        type(geometry),     intent(inout) :: geo
        type(system_state), intent(inout) :: st 
        
        integer                           :: r, period, temp(geo%sites), kx, ky, maxorb
        integer(kind=8)                   :: j, cntr, cntr2, symdim 
        integer(kind=8), allocatable      :: momBasis_temp(:), period_temp(:), orbits2D_temp(:,:,:), orbarr(:,:)
        double precision, allocatable     :: norm_temp(:), norm2D_temp(:,:)
        double precision                  :: normalization, normalization2D(2)
        double complex, allocatable       :: phases2D_temp(:,:,:)

        if(par%tilted == 1) then 
            kx = st%k2
            ky = st%k1 
        else 
            kx = st%k1 
            ky = st%k2 
        end if 
        temp = 1
        
        !Define temporary arrays due to unknown final Hilbert space dimension of symmd block
        if(allocated(period_temp)) deallocate(period_temp)
        
    
        if(geo%id == 1) then     
            if(allocated(norm_temp)) deallocate(norm_temp)
            allocate(norm_temp(geo%dim))
            norm_temp = 0 
        else if(geo%id == 2) then
            geo%orbsize = size(geo%mir) * size(geo%rot) * geo%l1 * geo%l2
            if(allocated(norm2D_temp)) deallocate(norm2D_temp)
            if(allocated(orbits2D_temp)) deallocate(orbits2D_temp)
            if(allocated(phases2D_temp)) deallocate(phases2D_temp)
            allocate(norm2D_temp(geo%dim, 2))
            allocate(orbits2D_temp(geo%dim, geo%orbsize, 2))
            allocate(orbarr(geo%orbsize, 2))
            allocate(phases2D_temp(geo%dim, geo%orbsize, 2))
            norm2D_temp   = 0.d0 
            phases2D_temp = 0.d0
            orbits2D_temp = 0 
            orbarr = 0
        end if
        allocate(momBasis_temp(geo%dim))
        allocate(period_temp(geo%dim))
        normalization = 0.d0 
        momBasis_temp = 0
        period_temp   = 0
        cntr          = 0 
        cntr2         = 0 
        symdim        = 0 
        period        = 0 
        

        do j = 1, geo%dim
            if(geo%id == 1) then 
                if(par%tilted == 0) then 
                    call checkstate_rect(par, geo, st, geo%basis_states(j), period, normalization)
                else if(par%tilted == 1) then 
                    call checkstate(par, geo, st, geo%basis_states(j), period, normalization)
                end if 
                if(period >= 0) then !New representative state found 
                    cntr                = cntr + 1
                    momBasis_temp(cntr) = geo%basis_states(j)
                    period_temp(cntr)   = period       
                    norm_temp(cntr)     = normalization   
                end if
            else if(geo%id == 2) then 
                cntr = cntr + 1 !Number of representatives(tentative)
                call checkstate2D(par, geo, st, geo%basis_states(j), period, orbits2D_temp(cntr, 1:geo%orbsize, 1:2), norm2D_temp(cntr, 1:2), phases2D_temp(cntr, 1:geo%orbsize, 1:2))

                if(period >= 0) then !New representative state found 
                    momBasis_temp(cntr)  = geo%basis_states(j) !Representative
                    period_temp(cntr)    = period       !Orbit size 
                    if(cntr == 1) maxorb = period       !
                    if(period > maxorb) maxorb = period      !Size of largest orbit for array allocation(later)
                else if(period < 0) then 
                    cntr = cntr - 1 !Reset counter 
                end if
                    
            end if 

        end do

        symdim = cntr
        
        if(allocated(geo%period)) deallocate(geo%period)
        
        if(geo%id == 1) then 
            if(allocated(geo%norm)) deallocate(geo%norm)
            allocate(geo%norm(symdim))
            geo%norm = 0.d0  
            geo%norm(1:symdim) = norm_temp(1:symdim)
            if(allocated(norm_temp)) deallocate(norm_temp)
        else if(geo%id == 2) then 
            if(allocated(geo%norm2D)) deallocate(geo%norm2D)
            if(allocated(geo%orbits2D)) deallocate(geo%orbits2D)
            if(allocated(geo%phases2D)) deallocate(geo%phases2D)
            allocate(geo%norm2D(symdim, 2))
            allocate(geo%orbits2D(symdim, maxorb, 2))
            allocate(geo%phases2D(symdim, maxorb, 2))
            
            geo%norm2D   = 0.d0 
            geo%phases2D = 0.d0 
            geo%orbits2D = 0
            geo%norm2D(1:symdim, 1) = norm2D_temp(1:symdim, 1)
            geo%norm2D(1:symdim, 2) = norm2D_temp(1:symdim, 2)
            geo%orbits2D(1:symdim, 1:maxorb, 1) = orbits2D_temp(1:symdim, 1:maxorb, 1)
            geo%orbits2D(1:symdim, 1:maxorb, 2) = orbits2D_temp(1:symdim, 1:maxorb, 2)
            geo%phases2D(1:symdim, 1:maxorb, 1) = phases2D_temp(1:symdim, 1:maxorb, 1)
            geo%phases2D(1:symdim, 1:maxorb, 2) = phases2D_temp(1:symdim, 1:maxorb, 2)
            
            geo%orbsize = maxorb 
            if(allocated(norm2D_temp)) deallocate(norm2D_temp)
            if(allocated(orbits2D_temp)) deallocate(orbits2D_temp)
            if(allocated(phases2D_temp)) deallocate(phases2D_temp)
        end if 
        if(allocated(geo%momBasis)) deallocate(geo%momBasis)
        allocate(geo%momBasis(symdim), geo%period(symdim))
        geo%momBasis           = 0
        geo%period             = 0 
        geo%momBasis(1:symdim) = momBasis_temp(1:symdim)
        geo%period(1:symdim)   = period_temp(1:symdim)
        geo%dim = symdim
        if(allocated(momBasis_temp)) deallocate(momBasis_temp)
        if(allocated(period_temp))   deallocate(period_temp)

    end subroutine symm_basis

    !Checkstate for tilted lattices: For rotations, reflections and translations(single sum of phases for all orbits.)
    subroutine checkstate(par, geo, st, state, period, norm_sym)
        ! Checks if a state is compatible with momentum and point group symmetry for tilted lattices
        implicit none

        type(sim_params),   intent(inout) :: par
        type(geometry),     intent(inout) :: geo
        type(system_state), intent(inout) :: st
        integer(kind=8),    intent(in)    :: state
        integer,            intent(out)   :: period
        double precision,   intent(out)   :: norm_sym

        integer(kind=8)                   :: sr, s0, orbit = 1
        integer                           :: ntot, info
        integer                           :: i, n, flag, rt
        integer(kind=8), allocatable      :: orbits(:)
        double precision                  :: sign, kx(2), ky(2), k(2)
        double precision                  :: a1(2), a2(2)
        double precision, parameter       :: tolerance = 1.0e-8
        double complex                    :: phase

        
        ntot = popcnt(state) !Number of geo%particles 
        a1 = geo%nnnVec(1:2, 1)
        a2 = geo%nnnVec(1:2, 2)
        kx =(/(-2.d0 * pi) / 3.d0,(-2.d0 * pi) / sqrt(3.d0) /)
        ky =(/(-4.d0 * pi) / 3.d0, 0.d0 /)
        if(geo%nHel == 1) then 
            k  =(dble(st%k1)/dble(geo%l1)) * ky
        else if(geo%nHel > 1) then 
            k  =(dble(st%k1)/dble(geo%l1)) * ky +(dble(st%k2)/dble(geo%l2)) * kx 
        end if

        sign = 1.d0
        info = 0
        period = -1 !r is the orbit size. If r = -1, the state is not compatible with the momentum k or point group symmetry.
        rt = 0 
        s0 = state
        sr = 0
        flag  = 0 
        norm_sym  = 0.d0 
        orbit = 0 
        phase = 0.d0 
        geo%orbsize = size(geo%mir) * size(geo%rot) * geo%l1 * geo%l2
 

        if(allocated(orbits)) deallocate(orbits)
        allocate(orbits(geo%orbsize))
        
        orbits = 0 
        sign = 1.d0
        call translation(s0, state, geo%sites, ntot, orbit, orbits, 1, geo%nHel, geo%l1, geo%l2, geo%id, sign, a1, a2, geo%xtransl, geo%ytransl, k, phase, info)

        if(info < 0) return 

        if(par%symm == 0) goto 11
        !Reflections 
        do i = 1, 6
            sign = 1.d0   
            call reflect(s0, state, geo%sites, geo%refl(i,1:geo%sites), sign, info, sr) 
            if(info < 0) return        
            call translation(s0, sr, geo%sites, ntot, orbit, orbits, 1, geo%nHel, geo%l1, geo%l2, geo%mir(i), sign, a1, a2, geo%xtransl, geo%ytransl, k, phase, info)
            if(info < 0) return 
        end do 

        !Rotations
        do n = 1, 5  
            sign = 1
            call c6n(s0, state, geo%sites, n, geo%c6, sign, info, sr)
            if(info < 0) return 
            call translation(s0, sr, geo%sites, ntot, orbit, orbits, 1, geo%nHel, geo%l1, geo%l2, geo%rot(n), sign, a1, a2, geo%xtransl, geo%ytransl, k, phase, info)
            if(info < 0) return 
        end do 

        11 continue 
       
        if((abs(dble(phase)) < tolerance) .and.(abs(aimag(phase)) < tolerance)) return  

        do i = 1, size(orbits)
            if(orbits(i) > 0) rt = rt + 1 
        end do 
        period   = rt 
        norm_sym = period * abs(phase)**2
        return
   
    end subroutine checkstate

    subroutine checkstate_rect(par, geo, st, state, period, norm_sym)
        ! Checks if a state is compatible with momentum and point group symmetry for rectangular lattices
        
        implicit none
    
        type(sim_params),   intent(inout) :: par
        type(geometry),     intent(inout) :: geo
        type(system_state), intent(inout) :: st
        integer(kind=8),    intent(in)    :: state
        integer,            intent(out)   :: period
        double precision,   intent(out)   :: norm_sym
        
        integer(kind=8)                   :: sr, s0, orbit = 1
        integer                           :: ntot, info
        integer                           :: i, n, flag, rt
        integer(kind=8), allocatable      :: orbits(:)
        double precision                  :: sign, kx(2), ky(2), k(2)
        double precision                  :: a1(2), a2(2)
        double precision, parameter       :: tolerance = 1.0e-8
        double complex                    :: phase

        ntot = popcnt(state)
        ! a1 =(/sqrt(3.d0), 0.d0/)
        ! a2 = 0.5*(/-1*sqrt(3.d0), 3.d0/)
        a1 = geo%nnnVec(1:2, 1)
        a2 = geo%nnnVec(1:2, 2)
        kx =(/(2.0d0 * pi)/sqrt(3.d0),(2.d0 * pi)/3.d0/)
        ky =(/0.d0,(4.d0 * pi)/3.d0/)
        k  =(dble(st%k1)/dble(geo%l1)) * kx +(dble(st%k2)/dble(geo%l2)) * ky
        
        flag    = 0 
        orbit   = 0
        norm_sym = 0.d0 
        phase   = 0.d0 
        period  = - 1 
        rt      = 0 
        s0      = state
        geo%orbsize = size(geo%mir) * size(geo%rot) * geo%l1 * geo%l2
 
        if(allocated(orbits)) deallocate(orbits)
        allocate(orbits(geo%orbsize))    
        orbits = 0 
        sign   = 1 
        !Identity x translations
        call translation(s0, state, geo%sites, ntot, orbit, orbits, 0, geo%l2, geo%l1, geo%l2, geo%id, sign, a2, a1, geo%xtransl, geo%ytransl, k, phase, info)
        if(info < 0) return 
        
        if(par%symm == 0) goto 11
        !Reflections 
        do i = 1, 6
            sign = 1   
            call reflect(s0, state, geo%sites, geo%refl(i,1:geo%sites), sign, info, sr) 
            if(info < 0) return        
            call translation(s0, sr, geo%sites, ntot, orbit, orbits, 0, geo%nHel, geo%l1, geo%l2, geo%mir(i), sign, a2, a1, geo%xtransl, geo%ytransl, k, phase, info)
            if(info < 0) return 
        end do 

        !Rotations  
        do n = 1, 5  
            sign = 1
            call c6n(s0, state, geo%sites, n, geo%c6, sign, info, sr)
            if(info < 0) return 
            call translation(s0, sr, geo%sites, ntot, orbit, orbits, 0, geo%nHel, geo%l1, geo%l2, geo%rot(n), sign, a2, a1, geo%xtransl, geo%ytransl, k, phase, info)
            if(info < 0) return 
        end do 

        11 continue 
        
        if((abs(dble(phase)) < tolerance) .and.(abs(aimag(phase)) < tolerance)) return  

        do i = 1, size(orbits)
            if(orbits(i) > 0) rt = rt + 1 
        end do 
        period = rt 
        norm_sym = period * abs(phase)**2

        return

  
    end subroutine checkstate_rect                   

    subroutine checkstate2D(par, geo, st, state, period, orbits, norm, phases)
        ! Checks state compatibility for 2D irreducible representations with complex phases
        
        implicit none

        type(sim_params),   intent(inout) :: par
        type(geometry),     intent(inout) :: geo
        type(system_state), intent(inout) :: st
        integer(kind=8),    intent(in)    :: state
        integer,            intent(out)   :: period
        integer(kind=8),    intent(out)   :: orbits(geo%orbsize,2)
        double precision,   intent(out)   :: norm(2)
        double complex,     intent(out)   :: phases(1:geo%orbsize, 2)

        integer(kind=8)                   :: sr, s0
        integer                           :: rep, ntot, info, layers, flag
        integer                           :: orbit = 1, i, n, rt, c
        double precision                  :: sign, kx(2), ky(2), k(2)
        double precision                  :: a1(2), a2(2), dcntr
        double precision, parameter       :: tolerance = 1.0e-8
        
        double precision                  :: sigma(2,2), rho(2,2), reflection(2,2), rotation(2,2), identity(2,2) 
        double complex                    :: phase

        if(geo%rot(1) == 1)  rep = 1 !IRREP E1
        if(geo%rot(1) == -1) rep = 2 !IRREP E2
        rho(1,1)   =  cos(2*pi*rep/6) !Rotation matrix entry(1,1)
        rho(1,2)   = -sin(2*pi*rep/6) !Rotation matrix entry(1,2)
        rho(2,1)   =  sin(2*pi*rep/6) !Rotation matrix entry(2,1)
        rho(2,2)   =  cos(2*pi*rep/6) !Rotation matrix entry(2,2)
        sigma(1,1) =  1               !Reflection matrix entry(1,1)
        sigma(1,2) =  0               !Reflection matrix entry(1,2)
        sigma(2,1) =  0               !Reflection matrix entry(2,1)
        sigma(2,2) = -1               !Reflection matrix entry(2,2)
        identity = 0.d0 
        identity(1, 1) = 1.d0
        identity(2, 2) = 1.d0
        if(par%tilted == 1) then 
            a1     = geo%nnnVec(1:2, 1)
            a2     = geo%nnnVec(1:2, 2)
            layers = geo%nHel 
            kx =(/(-2.d0 * pi) / 3.d0,(-2.d0 * pi) / sqrt(3.d0) /)
            ky =(/(-4.d0 * pi) / 3.d0, 0.d0 /)
            if(geo%nHel == 1) then 
                k  =(dble(st%k1)/dble(geo%l1)) * ky
            else if(geo%nHel > 1) then 
                k  =(dble(st%k1)/dble(geo%l1)) * ky +(dble(st%k2)/dble(geo%l2)) * kx 
                end if   
            else 
                layers = geo%l2
                a1     = geo%nnnVec(1:2, 2)
                a2     = geo%nnnVec(1:2, 1)
                kx     =(/(2.0d0 * pi)/sqrt(3.d0),(2.d0 * pi)/3.d0/)
                ky     =(/0.d0,(4.d0 * pi)/3.d0/)
                k      =(dble(st%k1)/dble(geo%l1)) * kx +(dble(st%k2)/dble(geo%l2)) * ky 
            end if 
            ntot    = popcnt(state) !Number of geo%particles      
            period  = -1 
            s0      = state
            norm    = 0.d0 
            flag    = 0
            orbits  = 0 
            geo%phases2D = 0.d0
            do c = 1, 2 !IRREP basis states
                orbit  = 0            
                sign   = 1 
                rt     = 0 
                phase  = 0.d0 

                call translation2D(par, geo, s0, state, orbit, orbits(1:geo%orbsize, c), ntot, layers, 1.0d0, sign, a1, a2, k, phases(1:geo%orbsize, c), phase, info)

                if(info < 0) phases(1:geo%orbsize, c) = 0.d0 !All prefactors = 0.  
                if(info < 0) orbits(1:geo%orbsize, c) = 0    !Remove orbit. 
                if(info < 0) cycle
                
                !Reflections 
                reflection = 0.d0 
                do i = 1, 6
                    sign = 1
                    reflection = identity 
                    reflection(1,1) = cos(2*pi*(i-1)/3.0d0)
                    reflection(2,2) = -cos(2*pi*(i-1)/3.0d0)
    
                    call reflect(s0, state, geo%sites, geo%refl(i,1:geo%sites), sign, info, sr) 
                    if(info < 0) then 
                        period = -2        
                        exit 
                    end if

                    call translation2D(par, geo, s0, sr, orbit, orbits(1:geo%orbsize, c), ntot, layers, reflection(c,c), sign, a1, a2, k,phases(1:geo%orbsize, c), phase, info)

                    if(info < 0) then 
                        period = -2        
                        exit 
                    end if

                end do 
                if(period == -2) then 
                    period = -1  
                    phases(1:geo%orbsize, c) = 0.d0 !All prefactors = 0. 
                    orbits(1:geo%orbsize, c) = 0    !Remove orbit. 
                    cycle 
                end if 

                rotation = rho 
                !Rotations
                do n = 1, 5  
                    sign = 1
                    call c6n(s0, state, geo%sites, n, geo%c6, sign, info, sr)
                    if(info < 0) then 
                        period = -2        
                        exit 
                    end if

                    call translation2D(par, geo, s0, sr, orbit, orbits(1:geo%orbsize,c), ntot, layers, rotation(c,c), sign, a1, a2, k, phases(1:geo%orbsize, c), phase, info)
                    
                    if(info < 0) then 
                        period = -2        
                        exit 
                    end if
                    rotation = matmul(rotation, rho)
                    
                end do 

                if(period == -2) then 
                    period = -1 
                    phases(1:geo%orbsize, c) = 0.d0 !All prefactors = 0. 
                    orbits(1:geo%orbsize, c) = 0    !Remove orbit. 
                    cycle 
                end if 

                if((abs(dble(phase)) < tolerance) .and.(abs(aimag(phase)) < tolerance)) then !Basis state c has norm = 0, i.e. does not contribute to the symmetry state. Set all prefactors = 0. 
                    phases(1:geo%orbsize, c) = 0.d0 !All prefactors = 0. 
                    orbits(1:geo%orbsize, c) = 0    !Remove orbit. 
                    cycle 
                end if 
                
                if(period == -2) cycle

                    
                do i = 1, geo%orbsize
                    if(orbits(i,c) > 0) rt = rt + 1 
                end do 
                period = rt 
                norm(c) = rt * abs(phase)**2 
                dcntr = 0.d0 
                do i = 1, geo%orbsize
                    dcntr = dcntr + abs(phases(i, c))**2
                end do 
            end do 

            if(period == -2) period = -1 
        
            return
   
    end subroutine checkstate2D

end module basis

