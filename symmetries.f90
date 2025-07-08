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
            integer, intent(in) :: sites, ntot, tilted, nHel, Lx, Ly, signpg
            integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
            integer(kind=8), intent(in) :: s0, s
            double precision, intent(in) :: a1(2), a2(2), k(2), char
            integer, intent(out) :: info  
            integer(kind=8), intent(inout) :: orbit
            integer(kind=8), allocatable, intent(inout) :: orbits(:)
            double complex, intent(inout) :: phase

            integer :: nx, ny, flag
            integer(kind=8) :: rowst, edgeA, edgeB, i, t, tx, ty, sign 
            double precision :: shift(2)

            
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
        integer, intent(in) :: c6(sites)
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
                    ! if(i < c6(i)) sign = sign * (-1)**popcnt(ibits(t, i-1, c6(i) - 1)) 
                    sign = sign * (-1)**popcnt(ibits(t, 0, c6(i) - 1)) 

                    t = ibclr(t, c6(i) - 1)
                end if 
            end do 

        end do 

        if(sr < s0) info = - 1 !Rotated state is already contained
    
        return  
    
    end subroutine c6n

end module symmetries