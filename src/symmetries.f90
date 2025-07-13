module symmetries
    
    use parameters
    use functions

    implicit none

    contains 


    subroutine characters(symmetrize, irrep, mir, rot, id)
        !---------------------------------------!
        !            Character table            !
        !---------------------------------------!
        ! This subroutine provides the character table for the C6v point group.
        ! It returns the characters for the irreducible representations (irreps) of the group. 
        ! Currently only defined for the C6v point group.
        ! symmetrize = 0: No symmetrization, only identity
        ! symmetrize = 1: Symmetrization, use character table
        ! irrep = "A1", "A2", "B1", "B2", "E1", "E2"
        ! mir = reflection (mirror) characters, rot = rotation characters, id = identity character
        ! For C6v group, the characters are defined as
        ! mir(1) = sigma_v_1, mir(2) = sigma_v_2, mir(3) = sigma_v_3, mir(4) = sigma_d_1, mir(5) = sigma_d_2, mir(6) = sigma_d_3
        ! rot(1) = c6, rot(2) = (c6)^2 = c3, rot(3) = (c6)^3 = c2, rot(4) = (c6)^4 = -c3, rot(5) = (c6)^5 = -c6

        implicit none
        integer, intent(in) :: symmetrize
        character(len=2), intent(in) :: irrep
        double precision, intent(out) :: mir(6), rot(5), id

        if(symmetrize == 0) then 
            id  = 1
            mir = 1
            rot = 1    
        else if(irrep == "A1") then 
            id  = 1
            mir = 1
            rot = 1
        else if(irrep == "A2") then 
            id  = 1
            mir = -1
            rot = 1
        else if(irrep == "B1") then     
            id     = 1
            mir    = 1
            mir(4) = -1
            mir(5) = -1
            mir(6) = -1
            rot    = -1
            rot(2) = 1
            rot(4) = 1
        else if(irrep == "B2") then     
            id     = 1
            mir    = -1
            mir(4) = 1
            mir(5) = 1
            mir(6) = 1
            rot    = -1
            rot(2) = 1
            rot(4) = 1
        else if(irrep == "E1") then 
            id     = 2
            mir    = 0 
            rot(1) = +1 
            rot(5) = +1
            rot(2) = -1
            rot(4) = -1
            rot(3) = -2
        else if(irrep == "E2") then 
            id     = 2
            mir    = 0 
            rot(1) = -1 
            rot(5) = -1
            rot(2) = -1
            rot(4) = -1
            rot(3) = +2
        end if 

    end subroutine characters

    subroutine translation(s0, s, sites, ntot, orbit, orbits, tilted, nHel, Lx, Ly, char, signpg, a1, a2, xtransl, ytransl, k, phase, info)
            !Compares translated state to initially tested basis state, NOT the reflected or rotated state. Does not check for phase = 0.
            !If tilted == 0, swap input vectors a1 and a2 
            implicit none 
            integer, intent(in) :: sites, ntot, tilted, nHel, Lx, Ly
            integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
            integer(kind=8), intent(in) :: s0, s
            double precision, intent(in) :: a1(2), a2(2), k(2), char, signpg
            integer, intent(out) :: info  
            integer(kind=8), intent(inout) :: orbit
            integer(kind=8), allocatable, intent(inout) :: orbits(:)
            double complex, intent(inout) :: phase

            integer :: nx, ny, flag
            integer(kind=8) :: rowst, edgeA, edgeB, i, t, tx, ty 
            double precision :: shift(2), sign

            
            i     = 0
            t     = 0 
            tx    = 0 
            ty    = 0
            nx    = 0
            ny    = 0
            info  = 0    !On output: InDicates whether new state was found (info = 0) or state is already contained in list (info < 1).
            sign  = 1
            edgeA = 0
            edgeB = 0
            rowSt = 0
            orbit = orbit + 1    !Include initial state (if not already in the list) to orbit.
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
                        sign = sign * (-1)**popcnt( ibits( tx, rowst, sites - 2 ) ) !Occupation of rest of helix 
                    end if 
                else !More than one helix
                    do i = 1, nHel !Calculate anticommutation sign from x-translations 
                        edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                        edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                        if(btest(tx, edgeB) .and. btest(tx, edgeA)) then !If both edges occupied, no sign
                            cycle 
                        else if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                            rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
                            sign = sign * (-1)**popcnt(ibits(tx, rowst, 2*Lx-2)) !Occupation of y-layer 'i'
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
                    if (t < s0) then !Representative is already in the list 
                        info = -1 
                        return
                    end if

                else if(nHel > 1) then 

                    do ny = 1, Ly 
                        ty = t
                        sign = sign * ((-1)**(ntot-popcnt(ibits(ty, 2*Lx*(nHel-1), 2*Lx))))**popcnt(ibits(ty, 2*Lx*(nHel-1), 2*Lx)) 
                        if(tilted == 1) then !Minus sign from re-ordering after y-translation. Only if cluster is tilted.  
                            do i = 1, nHel
                                edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                                edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                                if(btest(ty, edgeB) .and. btest(ty, edgeA)) then !If both edges occupied, no sign
                                    cycle 
                                else if(btest(ty, edgeB) .neqv. btest(ty, edgeA)) then !If only one edge occupied 
                                    rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
                                    sign = sign * (-1)**popcnt( ibits( ty, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'
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
                        if (t < s0) then !Representative is already in the list 
                            info = -1
                            return
                        end if
                    end do
                end if 
                
            end do

            if(t .ne. s ) print*, 'Translation orbit incomplete.'
            if(t .ne. s ) stop

    end subroutine translation

    subroutine reflect(s0, s, sites, refl, sign, info, sr)

        implicit none

        integer(kind=8), intent(in) :: s0, s
        integer, intent(in) :: sites 
        integer, intent(in) :: refl(sites)
        integer(kind=8), intent(out) :: sr !Reflected state 
        double precision, intent(out) :: sign
        integer, intent(out) :: info 

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
        integer, intent(in) :: c6(sites)
        integer(kind=8), intent(out) :: sr !Rotated state 
        integer, intent(out) :: info 
        double precision, intent(out) :: sign 

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
                    ! if(i < c6(i)) sign = sign * (-1)**popcnt(ibits(t, i-1, c6(i) - 1)) 
                    sign = sign * (-1)**popcnt(ibits(t, 0, c6(i) - 1)) 

                    t = ibclr(t, c6(i) - 1)
                end if 
            end do 

        end do 

        if(sr < s0) info = - 1 !Rotated state is already contained
    
        return  
    
    end subroutine c6n


    subroutine translation2D(s0, s, sites, ntot, orbsize, orbit, orbits, tilted, nHel, Lx, Ly, char, signpg, a1, a2, xtransl, ytransl, k, phases, phase, info)

        !If tilted == 0, swap input vectors a1 and a2 and set nHel = Ly 
        implicit none 
        integer, intent(in) :: orbsize, sites, ntot, tilted, nHel, Lx, Ly
        integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
        integer(kind=8), intent(in) :: s0, s
        double precision, intent(in) :: a1(2), a2(2), k(2), char, signpg
        integer, intent(out) :: info  
        integer, intent(inout) :: orbit
        integer(kind=8), intent(inout) :: orbits(orbsize)
        double complex, intent(inout) :: phases(orbsize), phase 

        integer :: nx, ny, i, rowst, edgeA, edgeB, flag
        integer(kind=8) :: t, tx, ty 
        double precision :: shift(2), sign

            t = 0 
            tx = 0 
            ty = 0
            nx = 0
            ny = 0
            info = 0 !On output: InDicates whether new state was found (info = 0) or state is already contained in list (info < 1)
            sign = 1.d0
            edgeA = 0
            edgeB = 0
            rowSt = 0
            flag  = 1 
            ! !Include initial state (if not already in the list) to orbit 
            ! orbit = orbit + 1 
            

            ! do i = 1, orbit 
            !     if(orbits(i) == s) then !State already in orbit
            !         flag = 0 
            !         exit 
            !     end if 
            ! end do 
            ! if(flag == 1) orbits(orbit) = s !New state found, add to orbit

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
                        sign = sign * (-1)**popcnt( ibits( tx, rowst, sites - 2 ) ) !Occupation of rest of helix 
                    end if 
                else !More than one helix
                    do i = 1, nHel!Ly !Calculate anticommutation sign from x-translations 
                        edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                        edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                        if(btest(tx, edgeB) .and. btest(tx, edgeA)) then !If both edges occupied, no sign
                            cycle 
                        else if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                            rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
                            sign = sign * (-1)**popcnt(ibits(tx, rowst, 2*Lx-2)) !Occupation of y-layer 'i'
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
                            orbit = orbit - 1
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
                        sign  = sign * ((-1)**(ntot-popcnt(ibits(ty, 2*Lx*(nHel-1), 2*Lx))))**popcnt(ibits(ty, 2*Lx*(nHel-1), 2*Lx)) 
                        if(tilted == 1) then !Minus sign from re-ordering after y-translation. Only if cluster is tilted.  
                            do i = 1, nHel
                                edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                                edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                                if(btest(ty, edgeB) .and. btest(ty, edgeA)) then !If both edges occupied, no sign
                                    cycle 
                                else if(btest(ty, edgeB) .neqv. btest(ty, edgeA)) then !If only one edge occupied 
                                    rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
                                    sign = sign * (-1)**popcnt( ibits( ty, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'
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
                                orbit = orbit - 1 
                                phases(i) = phases(i) + exp(-ii * dot_product(k, shift)) * sign * char * signpg 

                                if(t == s0) phase = phase + exp(-ii * dot_product(k, shift)) * sign * char * signpg 
                                if(s0 == 511 .and. t == s0) then 
                                    print*, 'trans2D.1'

                                    print* ,i, 'orbit position'
                                    print* ,sign, 'sign'
                                    print* ,signpg, 'signpg'    
                                    print* ,char, 'character'
                                    print* , exp(-ii * dot_product(k, shift)) * sign * char * signpg , ' phi+'
                                    print* ,phase, 'phase'
                                end if
                                exit 
                            end if 
                        end do 
                        if(flag == 1) then 
                            orbits(orbit) = t 
                            phases(orbit) = phases(orbit) + exp(-ii * dot_product(k, shift)) * sign * char * signpg 
                            if(t == s0) phase = phase + exp(-ii * dot_product(k, shift)) * sign * char * signpg
                            if(s0 == 511 .and. t == s0) then 
                                print*, 'trans2D.2'
                                    
                                print* ,orbit, 'orbit position'
                                print* ,sign, 'sign'
                                print* ,signpg, 'signpg'    
                                print* ,char, 'character'
                                print* , exp(-ii * dot_product(k, shift)) * sign * char * signpg , ' phi+'
                                print* ,phase, 'phase'
                            end if
                        end if 
                
        

                        if (t < s0) then !Representative is already in the list 
                            info = -1
                            return
                        end if
                    end do
                end if 
                    
            end do

    end subroutine translation2D

    !For tilted clusters
    subroutine representative_tilted(s, n, nHel, tilt, Lx, Ly, xtransl, ytransl, r, l1, l2, sign)
        !Finds the representative 'r' state for state 's'. 'n' is the number of sites and 'l' the number of shifts needed to translate 's' to 'r'.
        implicit none
        integer, intent(in) :: n, nHel, tilt 
        integer(kind=8), intent(in) :: s
        integer, intent(in) :: Lx, Ly
        integer, intent(in) :: xtransl(2, n), ytransl(2, n)
        integer(kind=8), intent(out) :: r
        integer, intent(out) :: l1, l2, sign 

        integer :: i, nx, ny, ntot, signt, itemp 
        integer :: edgeA, edgeB, rowSt
        integer(kind=8) :: t, tx, ty 

        ntot  = popcnt(s)
        itemp = 0 
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

        do nx = 1, Lx
            if(nx == 1) then 
                tx = s 
            else if(nx > 1) then 
                tx = t 
            end if  

            if(nhel == 1) then !Single helix, quasi-1D
                edgeB = n-1 !Last B site in helix
                edgeA = n-2 !Last A site in helix
                if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                    rowSt = 0 !First site 
                    signt = signt * (-1)**popcnt( ibits( tx, rowst, n-2 ) ) !Occupation of rest of helix 
                end if 
            else
                do i = 1, nHel !Calculate anticommutation sign from x-translations 
                    edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                    edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                    if(btest(tx, edgeB) .and. btest(tx, edgeA)) then !If both edges occupied, no sign
                        cycle 
                    else if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                        rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
                        signt = signt * (-1)**popcnt( ibits( tx, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'
                    end if 
                end do 
            end if 
            do i = 1, n
                call mvbits(tx, xtransl(1,i)-1, 1, t, xtransl(2,i)-1)      
            end do
            if(nHel == 1) then 
                if (t < r) then
                    sign = signt 
                    r    = t
                    l1   = ny
                    l2   = nx
                end if
            else if(nHel > 1) then 
                do ny = 1, Ly 
                    ty = t
                    if(nhel == 1) then !Single helix
                        signt = signt * ((-1)**(ntot-popcnt( ibits( ty, n - tilt, tilt ) ) ))**popcnt( ibits( ty, n - tilt, tilt ) ) 
                    else 
                        signt = signt * ((-1)**(ntot-popcnt( ibits( ty, 2*Lx*(nHel-1), 2*Lx ) ) ))**popcnt( ibits( ty, 2*Lx*(nHel-1), 2*Lx ) ) 
                        do i = 1, nHel
                            edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                            edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                            if(btest(ty, edgeB) .and. btest(ty, edgeA)) then !If both edges occupied, no sign
                                cycle 
                            else if(btest(ty, edgeB) .neqv. btest(ty, edgeA)) then !If only one edge occupied 
                                rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
                                signt = signt * (-1)**popcnt( ibits( ty, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'
                            end if 
                        end do     
                    end if 

                    do i = 1, n
                        call mvbits(ty, ytransl(1,i)-1, 1, t, ytransl(2,i)-1)
                    end do
                    if (t < r) then
                        sign = signt 
                        r    = t
                        l1   = ny
                        l2   = nx
                    end if
                end do
            end if 
        end do

        return 

    end subroutine representative_tilted


    subroutine representative(s, n, Lx, Ly, xtransl, ytransl, r, l1, l2, sign)
        !Finds the representative 'r' state for state 's'. 'n' is the number of sites and 'l' the number of shifts needed to translate 's' to 'r'.
        implicit none
        integer, intent(in) :: n
        integer(kind=8), intent(in) :: s
        integer, intent(in) :: Lx, Ly
        integer, intent(in) :: xtransl(2, n), ytransl(2, n)
        integer(kind=8), intent(out) :: r
        integer, intent(out) :: l1, l2, sign 

        integer :: i, nx, ny, ntot, signt, itemp 
        integer :: edgeA, edgeB, rowSt
        integer(kind=8) :: t, tx, ty 

        ntot  = popcnt(s)
        itemp = 0 
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

        do nx = 1, Lx
            if(nx == 1) then 
                tx = s 
            else if(nx > 1) then 
                tx = t 
            end if  

            do i = 1, Ly
                edgeB = i*2*Lx-1
                edgeA = i*2*Lx-2 
                if(btest(tx, edgeB) .and. btest(tx, edgeA)) then 
                    cycle 
                else if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then 
                    rowSt = (i-1)*2*Lx
                    signt = signt * (-1)**popcnt( ibits( tx, rowst, 2*Lx-2 ) )  
                end if 
            end do      

            do i = 1, n
                call mvbits(tx, xtransl(1,i)-1, 1, t, xtransl(2,i)-1)      
            end do

            do ny = 1, Ly
                ty    = t
                signt = signt * ((-1)**(ntot-popcnt( ibits( ty, 2*Lx*(Ly-1), 2*Lx ) ) ))**popcnt( ibits( ty, 2*Lx*(Ly-1), 2*Lx ) ) 
                do i = 1, n
                    call mvbits(ty, ytransl(1,i)-1, 1, t, ytransl(2,i)-1)       
                end do
                if (t < r) then
                    sign = signt 
                    r    = t
                    l1   = nx
                    l2   = ny
                end if
            end do 
        end do 

        return 

    end subroutine representative

    !---------------------------------------------!
    !            Find representative              !
    !---------------------------------------------!

    subroutine representative_irrep(s, n, nHel, tilt, Lx, Ly, symmetrize, id, par, rot, xtransl, ytransl, refl, c6, r, l1, l2, sign)
        !Finds the representative 'r' for state 's'. 'n' is the number of sites and 'l' the number of shifts needed to translate 's' to 'r'.
        implicit none
        integer(kind=8), intent(in) :: s
        integer, intent(in) :: n, nHel, tilt, symmetrize, Lx, Ly, xtransl(2, n), ytransl(2, n), refl(6, n), c6(n)
        double precision, intent(in) :: id, par(6), rot(5)
        integer(kind=8), intent(out) :: r
        integer, intent(out) :: l1, l2
        double precision, intent(out) :: sign 

        integer :: nx = 0, ny = 0, i = 0 
        integer :: edgeA = 0, edgeB = 0, rowSt = 0
        double precision :: signt = 1
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

        if(symmetrize == 0) go to 12
        do i = 1, 6 !Check for representatives among reflected states
            call mirror_rep(r, s, n, nHel, tilt, par(i), Lx, Ly, refl(i, 1:n), xtransl, ytransl, sign, l1, l2)
        end do 
        do i = 1, 5 !Check for representatives among rotated states
            call rotate_rep(r, s, n, nHel, tilt, i, rot(i), Lx, Ly, c6, xtransl, ytransl, sign, l1, l2)
        end do 
        12 continue

        return 

    end subroutine representative_irrep


    subroutine representative_irrep_rect(s, n, Lx, Ly, symmetrize, id, par, rot, xtransl, ytransl, refl, c6, r, l1, l2, sign)
            !Finds the representative 'r' for state 's'. 'n' is the number of sites and 'l' the number of shifts needed to translate 's' to 'r'.
            implicit none
            integer(kind=8), intent(in) :: s
            integer, intent(in) :: n, symmetrize, Lx, Ly, xtransl(2, n), ytransl(2, n), refl(6, n), c6(n)
            double precision, intent(in) :: id, par(6), rot(5)
            integer(kind=8), intent(out) :: r
            integer, intent(out) :: l1, l2 
            double precision, intent(out) :: sign 

            integer :: nx = 0, ny = 0, i = 0 
            integer :: edgeA = 0, edgeB = 0, rowSt = 0
            integer(kind=8) :: t = 0, tx = 0, ty = 0
            double precision :: signt = 1.d0

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

            if(symmetrize == 0) go to 12
            do i = 1, 6 !Check for representatives among reflected states
                call mirror_rep(r, s, n, Ly, -1, par(i), Lx, Ly, refl(i, 1:n), xtransl, ytransl, sign, l2, l1)
            end do 
            do i = 1, 5 !Check for representatives among rotated states
                call rotate_rep(r, s, n, Ly, -1, i, rot(i), Lx, Ly, c6, xtransl, ytransl, sign, l2, l1)
            end do 
            12 continue

            return 

    end subroutine representative_irrep_rect

    subroutine testsymm(sym, dim, mat)

        implicit none

        integer(kind=8), intent(in) :: dim
        double precision, intent(in) :: mat(dim, dim)
        double precision :: eps = 0.000001 
        logical, intent(out) :: sym

        integer :: i, j


        sym=.true.
        do i = 1, dim
            do j = 1, dim
                if(mat(i,j) - mat(j,i) > eps) then
                    sym=.false.

                    print*, i,j, 'i,j'
                    print*, mat(i,j), 'mat(i,j)'
                    print*, mat(j,i), 'mat(j,i)'

                    exit
                enDif
            end do
        end do
        if(.not. sym) then
            print*, 'MATRIX IS NON-SYMMETRIC!'
            print*,''
            print*,''
            error stop
        else
            print*, 'MATRIX IS SYMMETRIC'
        end if
        return

    end subroutine testsymm

    subroutine xtranslate(s, nHel, n, Lx, xtransl, sign, sx)
        implicit none 
        integer(kind=8), intent(in) :: s 
        integer, intent(in) :: nHel, n, Lx
        integer, intent(in) :: xtransl(2, n)
        double precision, intent(inout) :: sign 
        integer(kind=8), intent(out) :: sx 

        integer :: i = 0, edgeA = 0, edgeB = 0, rowSt = 0
        integer(kind = 8) :: t = 0


        t = 0
        if(nhel == 1) then !Single helix, quasi-1D
            edgeB = n-1 !Last B site in helix
            edgeA = n-2 !Last A site in helix
            if(btest(s, edgeB) .neqv. btest(s, edgeA)) then !If only one edge occupied 
                rowSt = 0 !First site 
                sign = sign * (-1)**popcnt( ibits( s, rowst, n-2 ) ) !Occupation of rest of helix 
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
        integer(kind=8), intent(in) :: s 
        integer, intent(in) :: nHel, tilt, n, Lx
        integer, intent(in) :: ytransl(2, n)
        double precision, intent(inout) :: sign 
        integer(kind=8), intent(out) :: sy

        integer :: i = 0, ntot = 0, edgeA = 0, edgeB = 0, rowSt = 0
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
                    sign = sign * (-1)**popcnt( ibits( s, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'
                end if 
            end do     
            11 continue 
        end if 

        do i = 1, n
            call mvbits(s, ytransl(1,i)-1, 1, t, ytransl(2,i)-1)
        end do
        sy = t

    end subroutine ytranslate

    subroutine mirror_rep(r, s, n, nHel, tilt, p, Lx, Ly, refl, xtransl, ytransl, sign, l1, l2)
        !Swap l1 and l2, set nHel = Ly, and tilt < 0, if cluster is not tilted, i.e. rectangular. 
        implicit none
        integer(kind=8), intent(in) :: s !original state
        integer, intent(in) :: n, nHel, tilt, Lx, Ly, refl(n), xtransl(2,n), ytransl(2,n)     
        double precision, intent(in) :: p
        integer(kind=8), intent(inout) :: r !potential rep (r<=s) 
        integer, intent(inout) :: l1, l2 
        double precision, intent(inout) :: sign !If 'r' remains unchanged, previous sign remains unchanged. Otherwise updated. 

        integer :: nx = 0, ny = 0, info = 0
        integer(kind=8) :: t = 0, tx = 0, ty = 0 
        double precision :: signp = 1.0d0, signt = 1.d0

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
        integer, intent(inout) :: l1, l2
        double precision, intent(inout) :: sign !If 'r' remains unchanged, previous sign remains unchanged. Otherwise updated. 

        integer :: nx = 0, ny = 0, info = 0
        integer(kind=8) :: t = 0, tx = 0, ty = 0 
        double precision :: signrot = 1.0d0, signt = 1.d0

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


end module symmetries