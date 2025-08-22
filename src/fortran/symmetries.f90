module symmetries
    ! Module for handling symmetries in the C6v point group and translations in a 2D lattice. 
    ! This module provides functions to calculate characters, perform translations, and find representatives of states under symmetries.
    use params
    use functions
    use types
    implicit none

    interface representative
        module procedure representative_rect
        module procedure representative_tilted
    end interface representative

    contains 
        
    subroutine characters(symm, irrep, geo)
        !---------------------------------------!
        !            Character table            !
        !---------------------------------------!
        ! This subroutine provides the character table for the C6v point group.
        ! symm = 0: No symmetrization, only identity
        ! symm = 1: Symmetrization, use character table
        ! irrep = "A1", "A2", "B1", "B2", "E1", "E2"
        ! mir: reflection(mirror) characters, rot: rotation characters: identity character
        ! For C6v group, the characters are defined as
        ! mir(1) = sigma_v_1, mir(2) = sigma_v_2, mir(3) = sigma_v_3, mir(4) = sigma_d_1, mir(5) = sigma_d_2, mir(6) = sigma_d_3
        ! rot(1) = c6, rot(2) =(c6)^2 = c3, rot(3) =(c6)^3 = c2, rot(4) =(c6)^4 = -c3, rot(5) =(c6)^5 = -c6

        implicit none
        
        integer,          intent(in)  :: symm
        character(len=2), intent(in)  :: irrep
        type(geometry),   intent(out) :: geo

        if(symm == 0) then 
            geo%id  = 1
            geo%mir = 1
            geo%rot = 1   
        else 
            if(irrep == "A1") then 
                geo%id  = 1
                geo%mir = 1
                geo%rot = 1
            else if(irrep == "A2") then 
                geo%id  = 1
                geo%mir = -1
                geo%rot = 1
            else if(irrep == "B1") then     
                geo%id     = 1
                geo%mir    = 1
                geo%mir(4) = -1
                geo%mir(5) = -1
                geo%mir(6) = -1
                geo%rot    = -1
                geo%rot(2) = 1
                geo%rot(4) = 1
            else if(irrep == "B2") then     
                geo%id     = 1
                geo%mir    = -1
                geo%mir(4) = 1
                geo%mir(5) = 1
                geo%mir(6) = 1
                geo%rot    = -1
                geo%rot(2) = 1
                geo%rot(4) = 1
            else if(irrep == "E1") then 
                geo%id     = 2
                geo%mir    = 0 
                geo%rot(1) = +1 
                geo%rot(5) = +1
                geo%rot(2) = -1
                geo%rot(4) = -1
                geo%rot(3) = -2
            else if(irrep == "E2") then 
            geo%id     = 2
            geo%mir    = 0 
            geo%rot(1) = -1 
            geo%rot(5) = -1
            geo%rot(2) = -1
            geo%rot(4) = -1
            geo%rot(3) = +2
            end if
        end if

    end subroutine characters

    subroutine translation(s0, s, sites, ntot, orbit, orbits, tilted, nHel, Lx, Ly, char, signpg, a1, a2, xtransl, ytransl, k, phase, info)
        !Compares translated state to initially tested basis state, NOT the reflected or rotated state. Does not check for phase = 0.
        !If tilted == 0, swap input vectors a1 and a2 
        implicit none 
        
        integer,                       intent(in)    :: sites, ntot, tilted, nHel, Lx, Ly
        integer,                       intent(in)    :: xtransl(2, sites), ytransl(2, sites)
        integer(kind=8),               intent(in)    :: s0, s
        double precision,              intent(in)    :: a1(2), a2(2), k(2), char, signpg
        integer,                       intent(out)   :: info  
        integer(kind=8),               intent(inout) :: orbit
        integer(kind=8), allocatable,  intent(inout) :: orbits(:)
        double complex,                intent(inout) :: phase

        integer          :: nx, ny, flag
        integer(kind=8)  :: rowst, edgeA, edgeB, i, t, tx, ty 
        double precision :: shift(2), sign
            
        i     = 0
        t     = 0 
        tx    = 0 
        ty    = 0
        nx    = 0
        ny    = 0
        info  = 0    !On output: InDicates whether new state was found(info = 0) or state is already contained in list(info < 1).
        sign  = 1
        edgeA = 0
        edgeB = 0
        rowSt = 0
        orbit = orbit + 1    !Include initial state(if not already in the list) to orbit.
        flag  = 1 

        do i = 1, orbit 
            if(orbits(i) == s) then !State already in orbit
                flag = 0 
                exit 
            end if 
        end do 
        if(flag == 1) orbits(orbit) = s !New state found, add to orbit

        do nx = 1, Lx !x translations
            if(nx == 1) then !Start with initial state s
                tx = s 
            else if(nx > 1) then 
                tx = t 
            end if  
            
            if(nhel == 1) then !Single helix, quasi-1D
                edgeB = sites-1 !Last B site in helix
                edgeA = sites-2 !Last A site in helix
                if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                    rowSt = 0 !First site 
                    sign = sign *(-1)**popcnt(ibits(tx, rowst, sites - 2 ) ) !Occupation of rest of helix 
                end if 
            else !More than one helix
                do i = 1, nHel !Calculate anticommutation sign from x-translations 
                    edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                    edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                    if(btest(tx, edgeB) .and. btest(tx, edgeA)) then !If both edges occupied, no sign
                        cycle 
                    else if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                        rowSt =(i-1)*2*Lx !First site on y-layer 'i' 
                        sign = sign *(-1)**popcnt(ibits(tx, rowst, 2*Lx-2)) !Occupation of y-layer 'i'
                    end if 
                end do 
            end if 

            do i = 1, sites
                call mvbits(tx, xtransl(1,i)-1, 1, t, xtransl(2,i)-1) !Translate one site in x direction
            end do

            if(nHel == 1) then !No y-translations for single helix
                orbit = orbit + 1 
                flag  = 1 
                do i = 1, orbit 
                    if(orbits(i) == t) then 
                        flag = 0 
                        exit 
                    end if 
                end do 
                if(flag == 1) orbits(orbit) = t 
                if(t == s0) then 
                    shift = nx * a2  
                    phase = phase + exp( - ii * dot_product(k, shift)) * sign * char * signpg 
                end if 
                if(t < s0) then !Representative is already in the list 
                    info = -1 
                    return
                end if

            else if(nHel > 1) then 

                do ny = 1, Ly 
                    ty = t
                    sign = sign *((-1)**(ntot-popcnt(ibits(ty, 2*Lx*(nHel-1), 2*Lx))))**popcnt(ibits(ty, 2*Lx*(nHel-1), 2*Lx)) 
                    if(tilted == 1) then !Minus sign from re-ordering after y-translation. Only if cluster is tilted.  
                        do i = 1, nHel
                            edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                            edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                            if(btest(ty, edgeB) .and. btest(ty, edgeA)) then !If both edges occupied, no sign
                                cycle 
                            else if(btest(ty, edgeB) .neqv. btest(ty, edgeA)) then !If only one edge occupied 
                                rowSt =(i-1)*2*Lx !First site on y-layer 'i' 
                                sign = sign *(-1)**popcnt(ibits(ty, rowst, 2*Lx-2)) !Occupation of y-layer 'i'
                            end if 
                        end do     
                    end if 
                    do i = 1, sites
                        call mvbits(ty, ytransl(1,i)-1, 1, t, ytransl(2,i)-1)
                    end do
                    orbit = orbit + 1 
                    flag  = 1 
                    do i = 1, orbit 
                        if(orbits(i) == t) then 
                            flag = 0 
                            exit 
                        end if 
                    end do 
                    if(flag == 1) orbits(orbit) = t 
                    if(t == s0) then 
                        shift = ny * a1 + nx * a2  
                        phase = phase + exp(-ii * dot_product(k, shift)) * sign * char * signpg 
                    end if 
                    if(t < s0) then !Representative is already in the list 
                        info = -1
                        return
                    end if
                end do
            end if 
            
        end do

        if(t .ne. s ) print*, 'Translation orbit incomplete.'
        if(t .ne. s ) stop

    end subroutine translation

    subroutine translation2D(par, geo, s0, s, orbit, orbits, ntot, nHel, char, signpg, a1, a2, k, phases, phase, info)
        !Compares translated state to initially tested basis state, NOT the reflected or rotated state. Does not check for phase=0.
        !If tilted==0, swap input vectors a1 and a2 and set nHel=Ly 
        implicit none 
        
        type(sim_params),              intent(in)    :: par
        type(geometry),                intent(in)    :: geo
        integer,                       intent(in)    :: ntot, nHel

        integer(kind=8),               intent(in)    :: s0, s
        double precision,              intent(in)    :: a1(2), a2(2), k(2), char, signpg
        integer,                       intent(out)   :: info  
        integer,                       intent(inout) :: orbit
        integer(kind=8),               intent(inout) :: orbits(geo%orbsize)
        double complex,                intent(inout) :: phases(geo%orbsize), phase 

        integer          :: nx, ny, i, rowst, edgeA, edgeB, flag, lx, ly
        integer(kind=8)  :: t, tx, ty
        double precision :: shift(2), sign

        lx = geo%l2
        ly = geo%l1
        t = 0 
        tx = 0 
        ty = 0
        nx = 0
        ny = 0
        info = 0 !On output: InDicates whether new state was found(info = 0) or state is already contained in list(info < 1)
        sign = 1.d0
        edgeA = 0
        edgeB = 0
        rowSt = 0
        flag  = 1 

        do nx = 1, Lx !x translations
            if(nx == 1) then !Start with initial state s
                tx = s 
            else if(nx > 1) then 
                tx = t 
            end if  
            
            if(nhel == 1) then !Single helix, quasi-1D
                edgeB = geo%sites-1 !Last B site in helix
                edgeA = geo%sites-2 !Last A site in helix
                if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                    rowSt = 0 !First site 
                    sign = sign *(-1)**popcnt(ibits(tx, rowst, geo%sites - 2 ) ) !Occupation of rest of helix 
                end if 
            else !More than one helix
                do i = 1, nHel!Ly !Calculate anticommutation sign from x-translations 
                    edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                    edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                    if(btest(tx, edgeB) .and. btest(tx, edgeA)) then !If both edges occupied, no sign
                        cycle 
                    else if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                        rowSt =(i-1)*2*Lx !First site on y-layer 'i' 
                        sign = sign *(-1)**popcnt(ibits(tx, rowst, 2*Lx-2)) !Occupation of y-layer 'i'
                    end if 
                end do 
            end if 

            do i = 1, geo%sites
                call mvbits(tx, geo%xtransl(1,i)-1, 1, t, geo%xtransl(2,i)-1) !Translate one site in x direction
            end do

            if(nHel == 1) then !No y-translations for single helix
                orbit = orbit + 1 
                flag  = 1 
                do i = 1, orbit 
                    if(orbits(i) == t) then 
                        flag = 0 
                        exit 
                    end if 
                end do 
                if(flag == 1) orbits(orbit) = t 
                
                shift = nx * a2  
                phases(orbit) = phases(orbit) + exp( - ii * dot_product(k, shift)) * sign * char * signpg 
                if(t == s0) phase = phase + phases(orbit)
                if(t < s0) then !Representative is already in the list 
                    info = -1 
                    return
                end if
            else if(nHel > 1) then 
                do ny = 1, Ly 
                    ty    = t
                    shift = ny * a1 + nx * a2  
                    sign  = sign *((-1)**(ntot-popcnt(ibits(ty, 2*Lx*(nHel-1), 2*Lx))))**popcnt(ibits(ty, 2*Lx*(nHel-1), 2*Lx)) 
                    if(par%tilted == 1) then !Minus sign from re-ordering after y-translation. Only if cluster is tilted.  
                        do i = 1, nHel
                            edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                            edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                            if(btest(ty, edgeB) .and. btest(ty, edgeA)) then !If both edges occupied, no sign
                                cycle 
                            else if(btest(ty, edgeB) .neqv. btest(ty, edgeA)) then !If only one edge occupied 
                                rowSt =(i-1)*2*Lx !First site on y-layer 'i' 
                                sign = sign *(-1)**popcnt(ibits(ty, rowst, 2*Lx-2)) !Occupation of y-layer 'i'
                            end if 
                        end do     
                    end if 
                    do i = 1, geo%sites 
                        call mvbits(ty, geo%ytransl(1,i)-1, 1, t, geo%ytransl(2,i)-1)
                    end do
                    orbit = orbit + 1 
                    flag  = 1 
                    do i = 1, orbit 
                        if(orbits(i) == t) then 
                            flag = 0 
                            orbit = orbit - 1 
                            phases(i) = phases(i) + exp(-ii * dot_product(k, shift)) * sign * char * signpg 

                            if(t == s0) phase = phase + exp(-ii * dot_product(k, shift)) * sign * char * signpg                             
                            exit 
                        end if 
                    end do 
                    if(flag == 1) then 
                        orbits(orbit) = t 
                        phases(orbit) = phases(orbit) + exp(-ii * dot_product(k, shift)) * sign * char * signpg 
                        if(t == s0) phase = phase + exp(-ii * dot_product(k, shift)) * sign * char * signpg
                    end if       

                    if(t < s0) then !Representative is already in the list 
                        info = -1
                        return
                    end if
                end do
            end if
        end do

    end subroutine translation2D

    subroutine representative_rect(s, n, Lx, Ly, symm, mir, rot, xtransl, ytransl, refl, c6, r, l1, l2, sign)
        !Finds the representative 'r' for state 's'. 'n' is the number of sites and 'l' the number of shifts needed to translate 's' to 'r'.
        implicit none
        
        integer(kind=8),  intent(in)  :: s
        integer,          intent(in)  :: n, symm, Lx, Ly, xtransl(2, n), ytransl(2, n), refl(6, n), c6(n)
        double precision, intent(in)  :: mir(6), rot(5)
        integer(kind=8),  intent(out) :: r
        integer,          intent(out) :: l1, l2
        double precision, intent(out) :: sign

        integer          :: nx, ny, i 
        integer          :: edgeA, edgeB, rowSt
        integer(kind=8)  :: t, tx, ty
        double precision :: signt

        l1    = Lx 
        l2    = Ly
        sign  = 1.d0
        signt = 1.d0
        r     = s
        t     = 0
        tx    = 0 
        ty    = 0 
        rowSt = 0 
        edgeA = 0 
        edgeB = 0 
        nx    = 0
        ny    = 0
        i     = 0


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

    end subroutine representative_rect

    subroutine representative_tilted(s, n, nHel, tilt, Lx, Ly, symm, mir, rot, xtransl, ytransl, refl, c6, r, l1, l2, sign)
        !Finds the representative 'r' for state 's'. 'n' is the number of sites and 'l' the number of shifts needed to translate 's' to 'r'.
        implicit none
        
        integer(kind=8),  intent(in)  :: s
        integer,          intent(in)  :: n, nHel, tilt, symm, Lx, Ly, xtransl(2, n), ytransl(2, n), refl(6, n), c6(n)
        double precision, intent(in)  :: mir(6), rot(5)
        integer(kind=8),  intent(out) :: r
        integer,          intent(out) :: l1, l2
        double precision, intent(out) :: sign

        integer          :: nx, ny, i 
        integer          :: edgeA, edgeB, rowSt
        integer(kind=8)  :: t, tx, ty
        double precision :: signt

        l2    = Lx 
        l1    = Ly
        sign  = 1
        r     = s
        t     = 0
        tx    = 0 
        ty    = 0 
        signt = 1.d0
        rowSt = 0 
        edgeA = 0 
        edgeB = 0 

        do nx = 1, Lx !Identity x translations 
            if(nx == 1) tx = s 
            if(nx > 1) tx = t
            call xtranslate(tx, nHel, n, Lx, xtransl, signt, t)
            if(t < r) then !New rep found
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
                    if(t < r) then
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
        ! Finds representative state among mirror-reflected configurations
        !Swap l1 and l2, set nHel = Ly, and tilt < 0, if cluster is not tilted, i.e. rectangular. 
        implicit none
        
        integer(kind=8),  intent(in)    :: s                                    ! original state
        integer,          intent(in)    :: n, nHel, tilt, Lx, Ly               
        integer,          intent(in)    :: refl(n), xtransl(2,n), ytransl(2,n)     
        double precision, intent(in)    :: p
        integer(kind=8),  intent(inout) :: r                                    ! potential rep (r<=s) 
        integer,          intent(inout) :: l1, l2 
        double precision, intent(inout) :: sign                                 ! If 'r' remains unchanged, previous sign remains unchanged. Otherwise updated. 

        integer          :: nx, ny, info
        integer(kind=8)  :: t, tx, ty 
        double precision :: signp, signt 

        nx    = 0
        ny    = 0
        info  = 0
        t     = 0
        tx    = 0
        ty    = 0
        signp = 1.0d0
        signt = 1.0d0

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
        ! Finds representative state among rotationally transformed configurations
        !Swap l1 and l2, set nHel = Ly, and tilt < 0, if cluster is not tilted, i.e. rectangular. 
        implicit none
        
        integer(kind=8),  intent(in)    :: s                                    ! original state
        integer,          intent(in)    :: n, nHel, tilt, nrot, Lx, Ly         
        integer,          intent(in)    :: c6(n), xtransl(2,n), ytransl(2,n)     
        double precision, intent(in)    :: eval
        integer(kind=8),  intent(inout) :: r                                    ! potential rep (r<=s) 
        integer,          intent(inout) :: l1, l2
        double precision, intent(inout) :: sign                                 ! If 'r' remains unchanged, previous sign remains unchanged. Otherwise updated. 

        integer          :: nx, ny, info
        integer(kind=8)  :: t, tx, ty 
        double precision :: signrot, signt 
        
        nx      = 0
        ny      = 0
        info    = 0
        t       = 0
        tx      = 0
        ty      = 0
        signrot = 1.0d0
        signt   = 1.0d0
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
        ! Applies mirror reflection operation to a many-body state with proper fermion sign
        implicit none
        
        integer(kind=8),  intent(in)  :: s0, s
        integer,          intent(in)  :: sites 
        integer,          intent(in)  :: refl(sites)
        integer(kind=8),  intent(out) :: sr ! Reflected state 
        integer,          intent(out) :: info 
        double precision, intent(out) :: sign
        
        integer(kind=8)  :: t
        integer          :: i, flag

        i    = 0
        t    = 0
        sr   = 0
        flag = 0
        info = 0
        sign = 1.0d0
        
        do i = 1, sites !Reflect state 
            call mvbits(s, i - 1, 1, sr, refl(i) - 1)
        end do

        t = sr 
        do i = 1, sites !Sign change due to reflection 
            if(btest(t, refl(i) - 1)) then 
                sign = sign * (-1)**popcnt(ibits(t, 0, refl(i) - 1)) 
                t = ibclr(t, refl(i) - 1)
            end if 
        end do 

        if(sr < s0) info = - 1 !Reflected state is already contained
    
        return  
    
    end subroutine reflect

    subroutine c6n(s0, s, sites, n, c6, sign, info, sr)
        ! Applies n successive C6 rotation operations to a many-body state with proper fermion sign
        implicit none
        
        integer(kind=8),  intent(in)  :: s0, s
        integer,          intent(in)  :: sites, n  
        integer,          intent(in)  :: c6(sites)
        integer(kind=8),  intent(out) :: sr ! Rotated state 
        integer,          intent(out) :: info 
        double precision, intent(out) :: sign

        integer(kind=8)  :: t, si
        integer          :: i, j

        i = 0
        t = 0
        sr = 0 
        si = s
        info = 0
        sign = 1.0d0

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
        ! Performs translation in x-direction with appropriate fermion anticommutation signs
        implicit none 
        
        integer(kind=8),  intent(in)    :: s 
        integer,          intent(in)    :: nHel, n, Lx
        integer,          intent(in)    :: xtransl(2, n)
        integer(kind=8),  intent(out)   :: sx 
        double precision, intent(inout) :: sign 

        integer           :: i, edgeA, edgeB, rowSt
        integer(kind = 8) :: t

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
        ! Performs translation in y-direction with appropriate fermion anticommutation signs
        !Set tilt < 0 if cluster is not tilted, i.e. rectangular. 
        implicit none 
        
        integer(kind=8),  intent(in)    :: s 
        integer,          intent(in)    :: nHel, tilt, n, Lx
        integer,          intent(in)    :: ytransl(2, n)
        integer(kind=8),  intent(out)   :: sy
        double precision, intent(inout) :: sign 

        integer           :: i, ntot, edgeA, edgeB, rowSt
        integer(kind = 8) :: t

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

    subroutine checkstate_tilted(par, geo, st, state, period, norm_sym)
        ! Checks if a state is compatible with momentum and point group symmetry for tilted lattices (single sum of phases for all orbits).
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

    end subroutine checkstate_tilted

    subroutine checkstate_rect(par, geo, st, state, period, norm_sym)
        ! Checks if a state is compatible with momentum and point group symmetry for rectangular lattices.
        
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
        ! a1 = (/sqrt(3.d0), 0.d0/)
        ! a2 = 0.5*(/-1*sqrt(3.d0), 3.d0/)
        a1 = geo%nnnVec(1:2, 1)
        a2 = geo%nnnVec(1:2, 2)
        kx = (/(2.0d0 * pi)/sqrt(3.d0),(2.d0 * pi)/3.d0/)
        ky = (/0.d0,(4.d0 * pi)/3.d0/)
        k  = (dble(st%k1)/dble(geo%l1)) * kx +(dble(st%k2)/dble(geo%l2)) * ky
        
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

    subroutine checkstate_irrep2d(par, geo, st, state, period, orbits, norm, phases)
        ! Checks state compatibility for 2D irreducible representations with complex phases. Works for both tilted and rectangular lattices.
        
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
   
    end subroutine checkstate_irrep2d
    
end module symmetries