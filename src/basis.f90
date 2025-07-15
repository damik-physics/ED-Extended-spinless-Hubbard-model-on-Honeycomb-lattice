module basis ! Module for defining the basis for the extended Hubbard model on the honeycomb lattice.
    use parameters
    use functions
    use symmetries
    use variables
    
    implicit none
    
    contains 

    ! Looks for the state s in the basis and returns its location in loc.
    subroutine findstate(dim, s, basis, loc)

        implicit none

        integer(kind=8), intent(in) :: dim, basis(dim)
        integer(kind=8), intent(in) :: s
        integer(kind=8), intent(out) :: loc

        integer(kind=8) :: left, right, mean

        left = 1
        right = dim
        do while (left <= right)
            mean = floor((left + right)/2.d0)
            if (basis(mean) < s) then
                left = mean + 1
            else if (basis(mean) > s) then
                right = mean - 1
            else
                loc = mean
                return
            end if
        end do
        loc = -1 !If no state is found
        return


    end subroutine findstate

    subroutine make_basis(par, geo)
    ! subroutine make_basis(ti, par%tilted, pat, nnnVec, geo%sites, geo%particles, dim, symm, par%ucx, par%ucy, geo%l1, geo%l2, geo%basis_states, geo%abasis, geo%bbasis, tilt, nHel, k1, k2, xtransl, ytransl, id, par, rot, refl, c6, period, norm, orbsize, orbits2D, phases2D, norm2D)

        implicit none
        type(sim_params), intent(inout) :: par
        type(geometry), intent(inout) :: geo
            ! integer, intent(in) :: ti, par%tilted, geo%sites, geo%particles, tilt, nHel, par%ucx, par%ucy, k1, k2, symm
            ! integer, allocatable, intent(in) :: xtransl(:,:), ytransl(:,:)
            ! double precision, intent(in) :: nnnVec(2,3), id, par(6), rot(5)
            ! character(len=*), intent(in) :: pat

            ! integer, intent(out) :: geo%l1, geo%l2, orbsize
            ! integer(kind=8), intent(out) :: dim
            ! integer(kind=8), allocatable, intent(out) :: geo%basis_states(:), geo%abasis(:), geo%bbasis(:),period(:), orbits2D(:,:,:)
            ! integer, allocatable, intent(out) :: geo%refl(:,:), geo%c6(:)
            ! double precision, allocatable, intent(out) :: norm(:), geo%norm2D(:,:)
            ! double complex, allocatable, intent(out) :: phases2D(:,:,:)
            
            integer :: a(geo%sites), temp = 0, I_in = 0, I_tail = 0, i, l
            integer(kind=8) :: symdim, amask, bmask, j
            integer(kind=8), allocatable :: momBasis(:)    

        geo%dim = int(fact(geo%sites) / (fact(geo%particles) * fact(max(1,geo%sites-geo%particles))),8) !Number of basis states


        if (geo%dim == 0) then
            print*, 'No basis states available.'
            return 
        end if

        !Permutations contains integer values I of basis states, perm_up/dn contain the integer values I_up/dn and the corresponDing J_up/dn
        if (allocated(geo%basis_states)) deallocate(geo%basis_states)
        allocate(geo%basis_states(geo%dim))
        if (allocated(geo%abasis)) deallocate(geo%abasis)
        allocate(geo%abasis(geo%dim))
        if (allocated(geo%bbasis)) deallocate(geo%bbasis)
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
        if(par%ti == 0) return 
        if(allocated(geo%refl)) deallocate(geo%refl)
        allocate(geo%refl(6, geo%sites))
        if(allocated(geo%c6)) deallocate(geo%c6)
        allocate(geo%c6(geo%sites))
        geo%refl = 0
        geo%c6   = 0    
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
            if (geo%nHel == 1) then
                    geo%l1 = 1
            else if (geo%nHel > 1) then    
                if(modulo(dble(geo%sites)/dble( (geo%nHel* geo%tilt)), dble (geo%nhel)) == 0.d0) then 
                    geo%l1 = geo%sites/ (geo%nHel * geo%tilt) 
                else if(modulo(dble(geo%sites)/dble( (geo%nHel* geo%tilt)), dble (geo%nhel)) >= 1.d0) then 
                    geo%l1 = geo%sites/geo%tilt 
                else if(modulo(dble(geo%sites)/dble( (geo%nHel* geo%tilt)), dble(geo%nhel)) < 1.d0) then 
                    geo%l1 = ceiling(geo%sites/ (geo%nHel * geo%tilt * modulo(dble(geo%sites)/dble( (geo%nHel* geo%tilt)), 1.d0)))
                end if 
            end if 
            geo%l2 = geo%sites/(2*geo%nHel) 
        end if 
        call symm_basis(par, geo)
        ! call symm_basis(par%tilted, geo%dim, geo%sites, nHel, geo%l2, geo%l1, k2, k1, symm, id, par, rot, nnnVec, geo%basis_states, xtransl, ytransl, refl, geo%c6, symdim, momBasis, period, norm, orbsize, orbits2D, phases2D, norm2D)
        
        geo%dim = symdim
        deallocate(geo%basis_states)
        allocate(geo%basis_states(geo%dim))
        geo%basis_states = momBasis 
        print*, 'Basis symmd.'

    end subroutine make_basis

    !Rotations, reflections and translations
    subroutine symm_basis(par, geo)
    ! subroutine symm_basis(tilted, dim, n, nHel, Lx, geo%l2
! , k1, k2, irrep, id, par, rot, nnnVec, basis, xtransl, ytransl, refl, c6, symdim, momBasis, period, norm, orbsize, orbits2D, phases2D, norm2D)
        
        implicit none
        type(sim_params), intent(inout) :: par
        type(geometry), intent(inout) :: geo
  
        ! integer(kind=8), intent(in) :: dim, basis(dim)
        ! integer, intent(in) :: par%tilted, n, nHel, k1, k2, irrep, Lx, geo%l2

        ! integer, intent(in) :: xtransl(2, n), ytransl(2, n), geo%refl(6, n), geo%c6(n)
        ! double precision, intent(in) :: nnnVec(2,3), id, par(6), rot(5)
        ! integer, intent(out) :: orbsize
        ! integer(kind=8), intent(out) :: symdim
        ! integer(kind=8), allocatable, intent(out) :: momBasis(:), orbits2D(:,:,:), period(:)
        ! double precision, allocatable, intent(out) :: norm(:), norm2D(:,:)
        ! double complex, allocatable, intent(out) :: phases2D(:,:,:)

        integer :: r = 0, temp(geo%sites), kx, ky, maxorb
        integer(kind=8) :: j = 0, cntr = 0, cntr2 = 0, symdim = 0 
        integer(kind=8), allocatable :: momBasis_temp(:), period_temp(:), orbits2D_temp(:,:,:), orbarr(:,:)
        double precision, allocatable :: norm_temp(:), norm2D_temp(:,:)
        double precision :: normalization = 0.d0
        double complex, allocatable :: phases2D_temp(:,:,:)
        
        if(par%tilted == 1) then 
            kx = k2
            ky = k1 
        else 
            kx = k1 
            ky = k2 
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
        r             = 0 
        

        do j = 1, geo%dim
            if(geo%id == 1) then 
                if(par%tilted == 0) then 
                    ! call checkstate_rect(geo%basis_states(j), geo%sites, geo%l1, geo%l2, kx, ky, irrep, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, r, normalization)
                    call checkstate_rect(geo)
                else if(par%tilted == 1) then 
                    ! call checkstate(geo%basis_states(j), geo%sites, nHel, geo%l1, geo%l2, kx, ky, irrep, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, r, normalization)
                    call checkstate(geo, par)
                end if 
                if(r >= 0) then !New representative state found 
                    cntr                = cntr + 1
                    momBasis_temp(cntr) = geo%basis_states(j)
                    period_temp(cntr)   = r       
                    norm_temp(cntr)     = normalization   
                end if
            else if(geo%id == 2) then 
                ! if(r >= 0) then !New representative state found 
                !     cntr                 = cntr + 1
                !     momBasis_temp(cntr)  = geo%basis_states(j)
                !     period_temp(cntr)    = r       
                !     if(geo%id == 1) norm_temp(cntr) = normalization   
                !     if(geo%id == 2) then                    
                
                ! allocate(orbits2D(r, 2), phases2D(r, 2))
                cntr = cntr + 1 !Number of representatives (tentative)
                call checkstate2D(geo, par)
!                 call checkstate2D(geo%basis_states(j), geo%orbsize, geo%sites, par%tilted, nHel, geo%l1, geo%l2
! , kx, ky, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, orbits2D_temp(cntr, 1:geo%orbsize, 1:2), norm2D_temp(cntr, 1:2), phases2D_temp(cntr, 1:geo%orbsize, 1:2), r)



                if(r >= 0) then !New representative state found 
                    ! cntr                 = cntr + 1 
                    momBasis_temp(cntr)  = geo%basis_states(j) !Representative
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
            ! do i = 1, geo%orbsize            
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

    !Checkstate for par%tilted lattices: For rotations, reflections and translations (single sum of phases for all orbits.)
    subroutine checkstate(s, par, geo, r, norm)
!     subroutine checkstate(s, sites, nHel, geo%l1, geo%l2
! , kx, ky, irrep, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, r, norm)

        implicit none
    
        ! Given momentum
        integer(kind=8), intent(in) :: s
        integer, intent(in) :: kx, ky
        integer, intent(in) :: geo%sites, nHel, irrep
        integer, intent(in) :: xtransl(2, geo%sites), ytransl(2, geo%sites), geo%refl(6, geo%sites), geo%c6(geo%sites)

        ! Define lattice parameters
        integer, intent(in) :: Lx   ! Number of unit cells in the a2-direction
        integer, intent(in) :: geo%l2
   ! Number of unit cells in the a1-direction
        double precision, intent(in) :: nnnVec(2,3), id, par(6), rot(5)
        integer, intent(out) :: r
        double precision, intent(out) :: norm

        
        integer(kind=8) :: sr = 0, s0 = 0, orbit = 1
        integer :: ntot = 0, info = 0, geo%orbsize = 0
        integer :: i = 0, n = 0, flag = 0, rt = 0
        integer(kind=8), allocatable :: orbits(:)
        double precision :: sign = 1.d0, k1(2), k2(2), k(2)
        double precision :: a1(2), a2(2)
        double precision, parameter :: tolerance = 1.0e-8
        double complex :: phase = 0.d0

        ntot = popcnt(s) !Number of geo%particles 
        a1 = nnnVec(1:2, 1)
        a2 = nnnVec(1:2, 2)
        k1 = (/(-2.d0 * pi) / 3.d0, (-2.d0 * pi) / sqrt(3.d0) /)
        k2 = (/(-4.d0 * pi) / 3.d0, 0.d0 /)
        if (geo%nHel == 1) then 
            k  = (dble(kx)/dble(Lx)) * k2
        else if (geo%nHel > 1) then 
            k  = (dble(kx)/dble(Lx)) * k2 + (dble(ky)/dble(geo%l2
)) * k1 
        end if

        r  = -1 !r is the orbit size. If r = -1, the state is not compatible with the momentum k or point group symmetry.
        rt = 0 
        s0 = s
        flag  = 0 
        norm  = 0.d0 
        orbit = 0 
        phase = 0.d0 
        geo%orbsize = size(par) * size(rot) * Lx * geo%l2
 

        if(allocated(orbits)) deallocate(orbits)
        allocate(orbits(orbsize))
        
        orbits = 0 
        ! orbits(1) = s
        sign = 1.d0
        call translation(s0, s, geo%sites, ntot, orbit, orbits, 1, nHel, Lx, geo%l2
, id, sign, a1, a2, xtransl, ytransl, k, phase, info)

        if(info < 0) return 

        if(irrep == 0) goto 11
        !Reflections 
        do i = 1, 6
            sign = 1.d0   
            call reflect(s0, s, geo%sites, geo%refl(i,1:geo%sites), sign, info, sr) 
            if(info < 0) return        
            call translation(s0, sr, geo%sites, ntot, orbit, orbits, 1, nHel, Lx, geo%l2
, par(i), sign, a1, a2, xtransl, ytransl, k, phase, info)
            if(info < 0) return 
        end do 

        !Rotations
        do n = 1, 5  
            sign = 1
            call c6n(s0, s, geo%sites, n, c6, sign, info, sr)
            if(info < 0) return 
            call translation(s0, sr, geo%sites, ntot, orbit, orbits, 1, nHel, Lx, geo%l2
, rot(n), sign, a1, a2, xtransl, ytransl, k, phase, info)
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

    subroutine checkstate_rect(s, geo%sites, Lx, geo%l2
, kx, ky, irrep, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, r, norm)

        implicit none
    
        ! Given momentum
        integer(kind=8), intent(in) :: s
        integer, intent(in) :: geo%sites, kx, ky, irrep
        integer, intent(in) :: xtransl(2, geo%sites), ytransl(2, geo%sites), geo%refl(6, geo%sites), geo%c6(geo%sites)
        ! Define lattice parameters
        integer, intent(in) :: Lx   ! Number of unit cells in the x-direction
        integer, intent(in) :: geo%l2
   ! Number of unit cells in the y-direction
        double precision, intent(in) :: nnnVec(2,3), id, par(6), rot(5)
        integer, intent(out) :: r
        double precision, intent(out) :: norm
        
        integer(kind=8) :: sr = 0, s0 = 0, orbit = 1
        integer :: ntot = 0, orbsize = 0, info = 0
        integer :: i = 0, n = 0, flag = 0, rt = 0
        integer(kind=8), allocatable :: orbits(:)
        double precision :: sign = 1.d0, k1(2), k2(2), k(2)
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
        k  = (dble(kx)/dble(Lx)) * k1 + (dble(ky)/dble(geo%l2
)) * k2 
        
        flag    = 0 
        orbit   = 0
        norm    = 0.d0 
        phase   = 0.d0 
        r       = - 1 
        rt      = 0 
        s0      = s
        orbsize = size(par) * size(rot) * Lx * geo%l2
 
        if(allocated(orbits)) deallocate(orbits)
        allocate(orbits(orbsize))    
        orbits = 0 
        sign   = 1 
        !Identity x translations
        call translation(s0, s, geo%sites, ntot, orbit, orbits, 0, geo%l2
, Lx, geo%l2
, id, sign, a2, a1, xtransl, ytransl, k, phase, info)
        if(info < 0) return 
        
        if(irrep == 0) goto 11
        !Reflections 
        do i = 1, 6
            sign = 1   
            call reflect(s0, s, geo%sites, geo%refl(i,1:geo%sites), sign, info, sr) 
            if(info < 0) return        
            call translation(s0, sr, geo%sites, ntot, orbit, orbits, 0, geo%l2
, Lx, geo%l2
, par(i), sign, a2, a1, xtransl, ytransl, k, phase, info)
            if(info < 0) return 
        end do 

        !Rotations  
        do n = 1, 5  
            sign = 1
            call c6n(s0, s, geo%sites, n, c6, sign, info, sr)
            if(info < 0) return 
            call translation(s0, sr, geo%sites, ntot, orbit, orbits, 0, geo%l2
, Lx, geo%l2
, rot(n), sign, a2, a1, xtransl, ytransl, k, phase, info)
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

    subroutine checkstate2D(s, orbsize, geo%sites, par%tilted, nHel, Lx, geo%l2
, kx, ky, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, orbits, norm, phases, r)

        implicit none
    
        ! Given momentum
        integer(kind=8), intent(in) :: s
        integer, intent(in) :: orbsize, kx, ky, geo%sites, par%tilted, nHel, xtransl(2, geo%sites), ytransl(2, geo%sites), geo%refl(6, geo%sites), geo%c6(geo%sites)
        ! Define lattice parameters
        integer, intent(in) :: Lx   ! Number of unit cells in the a2-direction
        integer, intent(in) :: geo%l2
   ! Number of unit cells in the a1-direction
        double precision, intent(in) :: nnnVec(2,3), id, par(6), rot(5)
        
        integer, intent(out) :: r
        integer(kind=8), intent(out) :: orbits(orbsize,2)
        double precision, intent(out) :: norm(2)
        double complex, intent(out) :: phases(orbsize, 2)

        integer(kind=8) :: sr = 0, s0 = 0
        integer :: rep, ntot, info, layers, flag
        integer :: orbit = 1, i, n, rt, c
        ! integer, allocatable :: orbits(:)
        double precision :: sign = 1.d0, k1(2), k2(2), k(2)
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
        if(par%tilted == 1) then 
            a1     = nnnVec(1:2, 1)
            a2     = nnnVec(1:2, 2)
            layers = nHel 
            k1 = (/(-2.d0 * pi) / 3.d0, (-2.d0 * pi) / sqrt(3.d0) /)
            k2 = (/(-4.d0 * pi) / 3.d0, 0.d0 /)
            if (geo%nHel == 1) then 
                k  = (dble(kx)/dble(Lx)) * k2
            else if (geo%nHel > 1) then 
                k  = (dble(kx)/dble(Lx)) * k2 + (dble(ky)/dble(geo%l2
)) * k1 
            end if   
        else 
            layers = geo%l2

            a1     = nnnVec(1:2, 2)
            a2     = nnnVec(1:2, 1)
            k1     = (/(2.0d0 * pi)/sqrt(3.d0), (2.d0 * pi)/3.d0/)
            k2     = (/0.d0, (4.d0 * pi)/3.d0/)
            k      = (dble(kx)/dble(Lx)) * k1 + (dble(ky)/dble(geo%l2
)) * k2 
        end if 
        ntot    = popcnt(s) !Number of geo%particles      
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
             
            call translation2D(s0, s, geo%sites, ntot, orbsize, orbit, orbits(1:orbsize, c), par%tilted, layers, Lx, geo%l2
, 1.0d0, sign, a1, a2, xtransl, ytransl, k, phases(1:orbsize, c), phase, info)

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
 
                call reflect(s0, s, geo%sites, geo%refl(i,1:geo%sites), sign, info, sr) 
                if(info < 0) then 
                    r = -2        
                    exit 
                end if
                call translation2D(s0, sr, geo%sites, ntot, orbsize, orbit, orbits(1:orbsize, c), par%tilted, layers, Lx, geo%l2
, reflection(c,c), sign, a1, a2, xtransl, ytransl, k, phases(1:orbsize, c), phase, info)

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
                call c6n(s0, s, geo%sites, n, c6, sign, info, sr)
                if(info < 0) then 
                    r = -2        
                    exit 
                end if
                call translation2D(s0, sr, geo%sites, ntot, orbsize, orbit, orbits(1:orbsize,c), par%tilted, layers, Lx, geo%l2
, rotation(c,c), sign, a1, a2, xtransl, ytransl, k, phases(1:orbsize, c), phase, info)
                
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
