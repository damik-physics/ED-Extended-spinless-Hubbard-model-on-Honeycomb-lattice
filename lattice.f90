module lattice

    use functions 
    implicit none
    

    contains 

    subroutine define_lattice(dir, tilted, ucx, ucy, nnbonds, nnnbonds, bc, pattern, cluster, bsites, hexsites, latticevecs, alattice, blattice, xyA, xyB, asitesbonds, bsitesbonds, cntrA, cntrB, nHel, tilt, phase, xy, xtransl, ytransl, reflections, nnnVec)

        ! use clusters 

        implicit none

        integer, intent(in) :: ucx, ucy, tilted
        character(len=*), intent(in) :: dir, bc, pattern, cluster
        integer, intent(out) :: nnbonds, nnnbonds, cntrA, cntrB, nHel, tilt 
        integer, allocatable, intent(out) :: bsites(:,:)
        integer, allocatable, intent(out) :: hexsites(:,:), phase(:), xy(:,:), latticevecs(:)
        integer, allocatable, intent(out) :: alattice(:,:), blattice(:,:), xyA(:,:), xyB(:,:)
        integer, allocatable, intent(out) :: asitesbonds(:,:), bsitesbonds(:,:)
        integer, allocatable, intent(out) :: xtransl(:,:), ytransl(:,:)
        integer, allocatable, intent(out) :: reflections(:)
        double precision, allocatable, intent(out) :: nnnVec(:,:)
        integer, allocatable :: sitecoord(:,:)
        integer :: i, nUC
        character :: name*200

        if(tilted == 0) then     
            call honeycomb(dir, ucx, ucy, bc, pattern, nnbonds, bsites)
            call honeycomb_nnn(dir, ucx, ucy, bc, pattern, nnnbonds, hexsites, latticevecs, alattice, blattice, asitesbonds, bsitesbonds, cntrA, cntrB, phase, xy, xtransl, ytransl, sitecoord, nnnVec)
            call coordinates(dir, 20, ucx, ucy, nnnbonds, cntrA, cntrB, xy, alattice, blattice, bc, pattern, xyA, xyB)
            call reflection(ucx, ucy, sitecoord, reflections)        
            print*, 'Honeycomb lattice created'
        else if(tilted == 1) then 
            if(cluster == '18A') then 
                nUC  = 9
                nHel = 3
                tilt = 2
            else if(cluster == '18B') then 
                nUC  = 9
                nHel = 1
                tilt = 4
            else if(cluster == '18C') then 
                nUC  = 9
                nHel = 1
                tilt = 8
            else if(cluster == '20A') then 
                nUC  = 10
                nHel = 1
                tilt = 6
            else if(cluster == '20B') then 
                nUC  = 10
                nHel = 5
                tilt = 2
            else if(cluster == '24C') then 
                nUC  = 12
                nHel = 1
                tilt = 8
            else if(cluster == '24D') then 
                nUC  = 12
                nHel = 2
                tilt = 2
            else if(cluster == '30A') then 
                nUC  = 15
                nHel = 1
                tilt = 20
            else if(cluster == '16A') then 
                nUC  = 8
                nHel = 2
                tilt = 2
            else if(cluster == '16B') then 
                nUC  = 8
                nHel = 1
                tilt = 4
            else if(cluster == 'L12') then 
                nUC  = 6
                nHel = 1
                tilt = 2
            else if(cluster == 'L18') then 
                nUC  = 9
                nHel = 3
                tilt = 2
            else if(cluster == 'L24') then 
                nUC  = 12
                nHel = 2
                tilt = 2
            else if(cluster == 'L24B') then 
                nUC  = 12
                nHel = 1
                tilt = 4
            else if(cluster == 'L26') then 
                nUC  = 13
                nHel = 1
                tilt = 6
            else if(cluster == 'L28') then 
                nUC  = 14
                nHel = 2
                tilt = 2
            else if(cluster == 'L32') then 
                nUC  = 16
                nHel = 4
                tilt = 2
            else if(cluster == 'L34') then 
                nUC  = 17
                nHel = 1
                tilt = 20
            else if(cluster == 'L36') then 
                nUC  = 18
                nHel = 3
                tilt = 2
            else if(cluster == 'L38') then 
                nUC  = 19
                nHel = 1
                tilt = 22
            else if(cluster == 'L42') then 
                nUC  = 21
                nHel = 1
                tilt = 8
            else if(cluster == 'L54') then 
                nUC  = 27
                nHel = 3
                tilt = 2
            end if
            
            write (name,"('Parameters_Cluster=',a,'_BC=',a,'.dat')") cluster,bc
            name = trim_name(dir//name)
            open(30,file=name)
            write(30,*) 'N_helices', nHel
            write(30,*) 'Tilt', tilt
            close(30)

            call build_cluster(nUC, nHel, tilt, cntrA, cntrB, nnbonds, nnnVec, alattice, blattice, hexsites, bsites, xtransl, ytransl)
            nnnbonds = cntrA + cntrB 
            
                 
            write(*,'(a,a,a)') ' Honeycomb lattice cluster ',cluster,' created' 
            if(allocated(xyA)) deallocate(xyA)
            if(allocated(xyB)) deallocate(xyB)
            allocate(xyA(4, cntrA))
            allocate(xyB(4, cntrB))
            xyA = 0 
            xyB = 0 

            write (name,"('Alattice_Cluster=',a,'_BC=',a,'.dat')") cluster,bc
            name = trim_name(dir//name)
            open(30,file=name)
            do i = 1, cntrA
                write(30,*) alattice(1, i), alattice(2, i), alattice(3, i), alattice(4, i), alattice(5, i)
            end do 
            close(30)
            write (name,"('Blattice_Cluster=',a,'_BC=',a,'.dat')") cluster,bc
            name = trim_name(dir//name)
            open(30,file=name)
            do i = 1, cntrB
                write(30,*) blattice(1, i), blattice(2, i), blattice(3, i), blattice(4, i), blattice(5, i)
            end do 
            close(30)
            write (name,"('NNlattice_Cluster=',a,'_BC=',a,'.dat')") cluster,bc
            name = trim_name(dir//name)
            open(30,file=name)
            do i = 1, nnbonds
                write(30,*) bsites(1, i), bsites(2, i)
            end do 
            close(30)
            write (name,"('NNNlattice_Cluster=',a,'_BC=',a,'.dat')") cluster,bc
            name = trim_name(dir//name)
            open(30,file=name)
            do i = 1, nnnbonds
                write(30,*) hexsites(1, i), hexsites(2, i)
            end do 
            close(30)
        end if 

    end subroutine define_lattice

    subroutine build_cluster(nUC, nHel, tilt, nabonds, nbbonds, nbonds, nnnVecs, sitesA, sitesB, nnnbonds, nnbonds, xtransl, ytransl)
        implicit none 
        integer, intent(in) :: nUC, nHel, tilt 

        integer, intent(out) :: nabonds, nbbonds, nbonds 
        integer, allocatable, intent(out) :: sitesA(:,:), sitesB(:,:), nnnbonds(:,:), nnbonds(:,:), xtransl(:,:), ytransl(:,:)
        double precision, allocatable, intent(out) :: nnnVecs(:,:)
        integer :: site1, site2, sites, nS
        integer :: cntrA, cntrB, cntr, xcntr, ycntr
        integer :: h, i

        if(allocated(nnnVecs)) deallocate(nnnVecs)
        allocate(nnnVecs(2,3))
        nnnVecs(1:2, 1) = (/0.d0, sqrt(3.d0)/) !a1 |
        nnnVecs(1:2, 2) = 0.5 * (/3.d0, -sqrt(3.d0)/) !a2 \
        nnnVecs(1:2, 3) = 0.5 * (/-3.d0, -sqrt(3.d0)/) !a3 /

        sites = nUC * 2
        nS = sites / nHel !Number of sites per helix in a2 direction 
        nabonds = nUC * 3 
        nbbonds = nUC * 3 
        nbonds  = nUC * 3 
        if(allocated(sitesA)) deallocate(sitesA)
        if(allocated(sitesB)) deallocate(sitesB)
        if(allocated(nnnbonds)) deallocate(nnnbonds)
        if(allocated(nnbonds)) deallocate(nnbonds)
        if(allocated(xtransl)) deallocate(xtransl)
        if(allocated(ytransl)) deallocate(ytransl)
        allocate(sitesA(5, nabonds))
        allocate(sitesB(5, nbbonds))
        allocate(nnnbonds(5, nabonds + nbbonds))
        allocate(nnbonds(2, nbonds))
        allocate(xtransl(2, sites))
        allocate(ytransl(2, sites))

        cntrA = 0 
        cntrB = 0
        xcntr = 0 
        ycntr = 0
        
        do h = 0, nHel - 1 !Loop over all helices
            do i = 0, nS - 1 !Loop over all sites per helix
                site1 = i + h * nS + 1
                if( modulo(i, 2) == 0 ) then 
                    
                    !a1 direction 
                    ycntr = ycntr + 1 
                    cntrA = cntrA + 1 
                    site2 = modulo(i + tilt, nS) + modulo(h+1, nHel) * nS + 1
                    sitesA(1, cntrA) = site1
                    sitesA(2, cntrA) = site2 
                    ytransl(1, ycntr) = site1
                    ytransl(2, ycntr) = site2
                    sitesA(3, cntrA) = 1
                    sitesA(4, cntrA) = cntrA + cntrB 
                    sitesA(5, cntrA) = 1
                        
                    !a2 direction 
                    xcntr = xcntr + 1
                    cntrA = cntrA + 1 
                    site2 = modulo(i + 2, nS) + h * nS + 1
                    sitesA(1, cntrA) = site1   
                    sitesA(2, cntrA) = site2 
                    sitesA(3, cntrA) = 1
                    sitesA(4, cntrA) = cntrA + cntrB 
                    sitesA(5, cntrA) = 2
                    xtransl(1, xcntr) = site1
                    xtransl(2, xcntr) = site2

                    !a3 direction 
                    cntrA = cntrA + 1 
                    site2 = modulo(modulo(i - tilt, nS) - 2, nS) + modulo(h-1, nHel) * nS + 1 
                    sitesA(1, cntrA) = site1 
                    sitesA(2, cntrA) = site2 
                    ! sitesA(2, cntrA) = modulo(site1 - 1 - nUp - 2, sites) + 1!? mod(..., sites)
                    sitesA(3, cntrA) = 1
                    sitesA(4, cntrA) = cntrA + cntrB 
                    sitesA(5, cntrA) = 3

                else if( modulo(i, 2) == 1 ) then  
                    
                    !a1 direction 
                    ycntr = ycntr + 1
                    cntrB = cntrB + 1 
                    site2 = modulo(i + tilt, nS) + modulo(h+1, nHel) * nS + 1
                    sitesB(1, cntrB) = site1
                    sitesB(2, cntrB) = site2 
                    ! if( i + nUp .le. nS - i) then 
                    !     ! site2 = modulo(site + nUp + nS, sites) !? mod(..., sites)
                    !     site2 = modulo(site1 - 1 + nUp + nS, sites) + 1 !? mod(..., sites)
                    !     ! site2 = modulo(site1 - 1 + nUp + h * nS, sites) + 1 !? mod(..., sites)
                    !     sitesB(2, cntrB) = site2 
                    ! else 
                    !     site2 = modulo(site1 - 1 + nUp, sites) + 1
                    !     sitesB(2, cntrB) = site2 
                    ! end if
                    sitesB(3, cntrB) = -1
                    sitesB(4, cntrB) = cntrA + cntrB 
                    sitesB(5, cntrB) = 1
                    ytransl(1, ycntr) = site1
                    ytransl(2, ycntr) = site2
                        
                    !a2 direction 
                    xcntr = xcntr + 1
                    cntrB = cntrB + 1 
                    site2 = modulo(i + 2, nS) + h * nS + 1
                    sitesB(1, cntrB) = site1   
                    sitesB(2, cntrB) = site2 
                    sitesB(3, cntrB) = -1
                    sitesB(4, cntrB) = cntrA + cntrB 
                    sitesB(5, cntrB) = 2
                    xtransl(1, xcntr) = site1
                    xtransl(2, xcntr) = site2

                    !a3 direction 
                    cntrB = cntrB + 1 
                    site2 = modulo(modulo(i - tilt, nS) - 2, nS) + modulo(h-1, nHel) * nS + 1 
                    ! site2 = modulo(site1 - 1 - nUp - 2, sites) + 1!? mod(..., sites)
                    sitesB(1, cntrB) = site1 
                    sitesB(2, cntrB) = site2 
                    sitesB(3, cntrB) = -1
                    sitesB(4, cntrB) = cntrA + cntrB 
                    sitesB(5, cntrB) = 3

                end if 
                if (cntrA > nabonds) print*,'Build cluster: Too many A bonds'
                if (cntrB > nbbonds) print*,'Build cluster: Too many B bonds'
            end do 
        end do 
        
        do i = 1, cntrA
            nnnbonds(1, 2*(i-1) + 1) = sitesA(1, i)
            nnnbonds(2, 2*(i-1) + 1) = sitesA(2, i)
            nnnbonds(1, 2*i) = sitesB(1, i)
            nnnbonds(2, 2*i) = sitesB(2, i)
        
            if(sitesA(1, i) > sites) print*,'Build cluster: A site ',sitesA(1, i), 'is larger than maximum number of sites.'
            if(sitesA(2, i) > sites) print*,'Build cluster: A site ',sitesA(2, i), 'is larger than maximum number of sites.'
            if(sitesB(1, i) > sites) print*,'Build cluster: B site ',sitesB(1, i), 'is larger than maximum number of sites.'
            if(sitesB(2, i) > sites) print*,'Build cluster: B site ',sitesB(2, i), 'is larger than maximum number of sites.'
        end do 

        cntr = 0
        do i = 1, cntrA, 3
            !e1
            cntr = cntr + 1
            nnbonds(1, cntr) = sitesA(1, i)
            nnbonds(2, cntr) = sitesB(1, i)
            !e2
            cntr = cntr + 1
            nnbonds(1, cntr) = sitesA(1, i)
            nnbonds(2, cntr) = sitesB(2, 3 * sitesB(2, i + 2) / 2 - 2)
            !e3
            cntr = cntr + 1
            nnbonds(1, cntr) = sitesA(1, i)
            nnbonds(2, cntr) = sitesB(2, i + 2)
        end do 
        
        if(cntr > nbonds) print*,'Build cluster: Too many NN bonds.'   

    end subroutine build_cluster

    subroutine reflection(ucx, ucy, sitecoord, reflections)
        implicit none 
        integer, intent(in) :: ucx, ucy
        integer, intent(in) :: sitecoord(2, *)
        
        integer, allocatable, intent(out) :: reflections(:)
        
        integer :: sites = 0 
        integer :: i = 0, j = 0, xr = 0, yr = 0 
        integer, allocatable :: uccoord(:,:)

        sites = 2 * ucx * ucy 
        if(allocated(uccoord)) deallocate(uccoord)
        allocate(uccoord(3, sites))
        if(allocated(reflections)) deallocate(reflections)
        allocate(reflections(sites))
        
        do j = 1, sites 
            uccoord(1, j) = sitecoord(1, j)/2 - 1
            uccoord(2, j) = sitecoord(2, j) - 1
            if(modulo(j,2) > 0) then 
                uccoord(3, j) = 0
            else if(modulo(j,2) == 0) then 
                uccoord(3, j) = 1
            end if 
        end do 
        i = 0
        j = 0
        xr = 0
        yr = 0 
        do j = 1, sites            

            xr = -uccoord(1,j)
            yr = -uccoord(2,j)
            do i = 1, sites
                if(uccoord(1,i) == xr .and. uccoord(2,i) == yr .and. uccoord(3,i) .ne. uccoord(3,j) ) then 
                    reflections(j) = i
                    exit 
                end if 
            end do 
        end do 
        return 

    end subroutine reflection

    subroutine coordinates(dir, unit, ucx, ucy, nbonds, abonds, bbonds, xy, alattice, blattice, bc, pattern, xyA, xyB)
        implicit none 
        integer, intent(in) :: unit, ucx, ucy, nbonds, abonds, bbonds
        integer, intent(in) :: xy(4, nbonds), alattice(5, abonds), blattice(5, bbonds)
        character(len=*), intent(in) :: dir, bc, pattern
        integer, allocatable, intent(out) :: xyA(:,:), xyB(:,:)

        integer   :: j = 0 
        character :: filename*300
        character :: parameters*300

        write (parameters,"('L=',i0,'ucx=',i0,'ucy=',i0,'BC=',a,'_pat=',a2'.dat')") 2*ucx*ucy,ucx,ucy,bc,pattern
        parameters = trim_name(parameters)
        if(allocated(xyA)) deallocate(xyA)
        if(allocated(xyB)) deallocate(xyB)
        allocate(xyA(4, abonds))
        allocate(xyB(4, bbonds))
        filename=dir//"A_coordinates_"//parameters
        filename=trim_name(filename)
        open(unit,file = filename)
        do j = 1, abonds             
            xyA(1,j) = xy(1,alattice(4,j))
            xyA(2,j) = xy(2,alattice(4,j))
            xyA(3,j) = xy(3,alattice(4,j))
            xyA(4,j) = xy(4,alattice(4,j))
            write(unit,'(4(I0,x))') xyA(1,j), xyA(2,j), xyA(3,j), xyA(4,j) 
        end do   
        close(unit)

        filename=dir//"B_coordinates_"//parameters
        filename=trim_name(filename)
        open(unit,file = filename)
        do j = 1, bbonds             
            xyB(1,j) = xy(1,blattice(4,j))
            xyB(2,j) = xy(2,blattice(4,j))
            xyB(3,j) = xy(3,blattice(4,j))
            xyB(4,j) = xy(4,blattice(4,j))
            write(unit,'(4(I0,x))') xyB(1,j), xyB(2,j), xyB(3,j), xyB(4,j) 
        end do 
        close(unit)

        return 
    end subroutine coordinates

    subroutine honeycomb(dir, ucx, ucy, bc, pattern, nbonds, bsites)
        !Subroutine obtained from: http://physics.bu.edu/~sandvik/vietri/sse/ssebasic.f90
        implicit none

        integer, intent(in) :: ucx, ucy
        character(len=*), intent(in) :: dir, bc, pattern

        integer, intent(out) :: nbonds
        integer, allocatable, intent(out) :: bsites(:,:)

        integer, allocatable :: bsites_temp(:,:)
        integer :: nn, nx, ny       ! number of sites
        integer :: y0
        integer :: xcounter, ycounter
        integer :: nbx, nby
        integer :: s, x1, x2, y1, y2
        integer, allocatable :: xy(:,:)    
        character :: name*200

        !d1 = (0,1)
        !d2 = (-Sqrt[3]/2,-1/2)
        !d3 = (+Sqrt[3]/2,-1/2)
        !First site on A lattice (pattern AB)
        xcounter = 0
        ycounter = 0

        nn = 2 * ucx * ucy !Number of lattice sites
        nx = 2 * ucx
        ny = ucy

        if(bc == 'o') then !nx = number of bonds per x-layer
            nbx = (nx - 1) * ucy !total number of x bonds
            nby = ((nx-1) + y0)/2 * (ny-1) !total number of y bonds
        else if(bc == 'p') then
            nbx = nx * ucy !total number of x bonds
            ! nby = int(nx/2) * ny !total number of y bonds
            nby = ucx * ny !total number of y bonds
        end if
        nbonds = nbx + nby
        ! print*, 'Number of NN bonds', nbonds 

        if (allocated(bsites_temp)) deallocate(bsites_temp)
        allocate(bsites_temp(2,nbonds))
        allocate(xy(4,nbonds))

        do y1 = 0, ny - 1
            do x1 = 0, nx - 1
                if (mod(x1,2) == 0 .and. nbx + ycounter <= nbonds) then !y-bonds d1 
                    if(y1 == ny - 1 .and. bc == 'o') goto 55
                    x2 = modulo(x1 + 1, nx) !x-coordinate of second site
                    y2 = modulo(y1 + 1, ny) !y-coordinate of second site
                    if (y2 == y1) goto 55 !debug For 1D case no y-hopping
                    ycounter = ycounter + 1 !Current y-bond = nbx + ycounter
                    bsites_temp(1,nbx + ycounter) = 1 + x1 + y1 * nx
                    bsites_temp(2,nbx + ycounter) = 1 + x2 + y2 * nx
                    ! print*, 'Ysites', bsites_temp(1,nbx + ycounter), bsites_temp(2,nbx + ycounter)
                    xy(1, nbx + ycounter) = x1
                    xy(2, nbx + ycounter) = y1
                    xy(3, nbx + ycounter) = x2
                    xy(4, nbx + ycounter) = y2                
                end if

                55 continue
                if (x1 == nx-1 .and. bc == 'o') cycle !x-bonds d3 
                xcounter = xcounter + 1
                x2 = modulo(x1 + 1, nx)
                y2 = y1
                bsites_temp(1,xcounter) = 1 + x1 + y1 * nx
                bsites_temp(2,xcounter) = 1 + x2 + y2 * nx
                ! print*, 'Xsites', bsites_temp(1,xcounter), bsites_temp(2,xcounter)                        
                xy(1, xcounter) = x1
                xy(2, xcounter) = y1
                xy(3, xcounter) = x2
                xy(4, xcounter) = y2            
            end do
        end do

        nbonds = xcounter + ycounter
        allocate(bsites(2,nbonds))
        bsites(1,1:xcounter) = bsites_temp(1,1:xcounter)
        bsites(2,1:xcounter) = bsites_temp(2,1:xcounter)
        bsites(1,xcounter + 1:xcounter  + ycounter) = bsites_temp(1,nbx + 1:nbx  + ycounter)
        bsites(2,xcounter + 1:xcounter  + ycounter) = bsites_temp(2,nbx + 1:nbx  + ycounter)

        deallocate(bsites_temp)
        write (name,"('L=',i0,'ucx=',i0,'ucy=',i0,'pat=',a2,'BC=',a,'.dat')") nn,ucx,ucy,pattern,bc
        name = trim_name(name)
        name=dir//"NN_lattice_"//name
        open(unit=31, file=name)
        do s = 1, nbonds
            write(31,*) bsites(1,s), bsites(2,s)
        end do
        close(31)

        write (name,"('L=',i0,'ucx=',i0,'ucy=',i0,'pat=',a2,'BC=',a,'.dat')") nn,ucx,ucy,pattern,bc
        name = trim_name(name)
        name=dir//"NN_xy_"//name
        open(unit=32, file=name)
        do s = 1, nbonds
            write(32,*) xy(1,s), xy(2,s), xy(3,s), xy(4,s)
        end do
        close(32)

        ! print*, 'Finished honeycomb lattice'

    end subroutine honeycomb

    subroutine honeycomb_nnn(dir, ucx, ucy, bc, pattern, nbonds, hexsites, latticevecs, alattice, blattice, asitesbonds, bsitesbonds, cntrA, cntrB, phase, xy, xtransl, ytransl, sitecoord, nnnVec)
        !Subroutine obtained from: http://physics.bu.edu/~sandvik/vietri/sse/ssebasic.f90
        implicit none

        integer, intent(in) :: ucx, ucy
        character(len=*), intent(in) :: dir, bc, pattern
        integer, intent(out) :: nbonds
        integer, intent(out) :: cntrA
        integer, intent(out) :: cntrB
        integer, allocatable, intent(out) :: hexsites(:,:)
        integer, allocatable, intent(out) :: latticevecs(:)
        integer, allocatable, intent(out) :: alattice(:,:)
        integer, allocatable, intent(out) :: blattice(:,:)
        integer, allocatable, intent(out) :: asitesbonds(:,:)
        integer, allocatable, intent(out) :: bsitesbonds(:,:)
        integer, allocatable, intent(out) :: phase(:)
        integer, allocatable, intent(out) :: xy(:,:)
        integer, allocatable, intent(out) :: xtransl(:,:)
        integer, allocatable, intent(out) :: ytransl(:,:)
        integer, allocatable, intent(out) :: sitecoord(:,:)
        double precision, allocatable, intent(out) :: nnnVec(:,:)

        ! integer, allocatable :: hexsites_temp(:,:), phase_temp(:)
        integer :: nn, nx, ny
        integer :: y0
        integer :: xcounter, ycounter, ytcounter
        integer :: nbx, nby
        integer :: s, x1, x2, y1, y2
        ! integer, allocatable :: xy(:,:)
        integer, allocatable :: alatt(:,:)
        integer, allocatable :: blatt(:,:)
        integer, allocatable :: cntrAsites(:)
        integer, allocatable :: cntrBsites(:)
        
        character :: name*200

        !Definition of lattice vectors: 
        ! v1 = {sqrt(3), 0}
        ! v2 = 1/2*{-sqrt(3), 3}
        ! v3 = 1/2*{-sqrt(3), -3}
        if(allocated(nnnVec)) deallocate(nnnVec)
        allocate(nnnVec(2, 3))
        nnnVec(1:2,1) = (/sqrt(3.d0), 0.d0/)
        nnnVec(1:2,2) = 0.5*(/-1*sqrt(3.d0), 3.d0/)
        nnnVec(1:2,3) = 0.5*(/-1*sqrt(3.d0), -3.d0/)

        ! print*, 'Start creating Honeycomb sublattices'
        xcounter = 0
        ycounter = 0
        ytcounter = 0
        cntrA    = 0 
        cntrB    = 0 
        if(pattern == 'AB') then !Determines pattern in x-direction: ABAB... = 1; BABA... = 0
            y0 = 1
        else if(pattern == 'BA') then
            y0 = 0
        end if
        nn = 2 * ucx * ucy !Number of lattice sites
        nx = 2 * ucx !nx = number of bonds per x-layer (PBC)
        ny = ucy !Number of y layers

        if(bc == 'o') then
            nbx = 2 * (ucx - 1) * ucy !total number of x bonds: 2 = # sublattices
            nby = 2 * (ucy - 1) * ( 2 * (ucx-1) + 1 ) !total number of y bonds: 2 = # sublattices; ucx-1 = unit cells with 2 bonds/site; 1 = unit cell with 1 bond per site (first uc)
        else if(bc == 'p') then
            nbx = nx * ny !total number of x bonds
            nby = 2 * nx * ny !total number of y bonds
        end if
        nbonds = nbx + nby
        ! print*, 'Number of NNN bonds', nbonds 
        

        if (allocated(hexsites)) deallocate(hexsites)
        if (allocated(latticevecs)) deallocate(latticevecs)
        if (allocated(phase)) deallocate(phase)
        if (allocated(xy)) deallocate(xy)
        if (allocated(asitesbonds)) deallocate(asitesbonds)
        if (allocated(bsitesbonds)) deallocate(bsitesbonds)
        if (allocated(cntrAsites)) deallocate(cntrAsites)
        if (allocated(cntrBsites)) deallocate(cntrBsites)
        if(allocated(xtransl)) deallocate(xtransl)
        if(allocated(ytransl)) deallocate(ytransl)
        if(allocated(sitecoord)) deallocate(sitecoord)
        
        allocate(xtransl(2, nn))
        allocate(ytransl(2, nn))
        allocate(hexsites(2,nbonds))
        allocate(latticevecs(nbonds))
        allocate(alatt(5,nbonds))
        allocate(blatt(5,nbonds))
        allocate(phase(nbonds))
        allocate(xy(4, nbonds))
        allocate(asitesbonds(6,nn))
        allocate(bsitesbonds(6,nn))
        allocate(cntrAsites(nn))
        allocate(cntrBsites(nn))
        allocate(sitecoord(2,nn))


        hexsites  = 0
        latticevecs = 0 
        phase       = 0
        cntrAsites  = 0
        cntrBsites  = 0
        asitesbonds = 0 
        bsitesbonds = 0 
        sitecoord   = 0 

        do y1 = 0, ny - 1
            do x1 = 0, nx - 1
                if(y1 == ny - 1 .and. bc == 'o') goto 65 !No y-bonds for last row at OBC 
                if (nbx + ycounter <= nbonds) then !y-bonds

                    !First vertical bond: Down left (lattice vector v3)
                    if((x1 == 0 .or. x1 == 1) .and. bc == 'o') goto 60 !No down-left NNN-bond for 1st and 2nd sites             
                    x2 = modulo(x1 - 2, nx) !x-coordinate of second site 
                    y2 = modulo(y1 - 1, ny) !y-coordinate of second site
                    if (y2 == y1) goto 60 !debug For 1D case no y-hopping
                    ycounter = ycounter + 1 !Current y-bond = nbx + ycounter
                    hexsites(1,nbx + ycounter) = 1 + x1 + y1 * nx
                    hexsites(2,nbx + ycounter) = 1 + x2 + y2 * nx     
                    latticevecs(nbx + ycounter)  = 3
                    ! print*, 'v3sites', hexsites(1,nbx + ycounter), hexsites(2,nbx + ycounter)
                    if( pattern == 'AB' .and. modulo(x1, 2) == 0 ) then 
                        cntrA = cntrA + 1
                        alatt(1, cntrA) = hexsites(1,nbx + ycounter)
                        alatt(2, cntrA) = hexsites(2,nbx + ycounter)
                        alatt(3, cntrA) = 1 !QAH phase 
                        alatt(4, cntrA) = nbx + ycounter !Bond number 
                        alatt(5, cntrA) = 3 !Lattice vector 
                        phase(nbx + ycounter) = 1
                        cntrAsites(hexsites(1,nbx + ycounter)) = cntrAsites(hexsites(1,nbx + ycounter)) + 1
                        asitesbonds(cntrAsites(hexsites(1,nbx + ycounter)),hexsites(1,nbx + ycounter)) = cntrA
                        cntrAsites(hexsites(2,nbx + ycounter)) = cntrAsites(hexsites(2,nbx + ycounter)) + 1
                        asitesbonds(cntrAsites(hexsites(2,nbx + ycounter)),hexsites(2,nbx + ycounter)) = -cntrA
                    else if( pattern == 'AB' .and. modulo(x1, 2) == 1 ) then 
                        cntrB = cntrB + 1
                        blatt(1, cntrB) = hexsites(1,nbx + ycounter)
                        blatt(2, cntrB) = hexsites(2,nbx + ycounter)                    
                        blatt(3, cntrB) = - 1 !QAH phase 
                        blatt(4, cntrB) = nbx + ycounter !Bond number
                        blatt(5, cntrB) = 3 !Lattice vector 
                        phase(nbx + ycounter) = - 1
                        cntrBsites(hexsites(1,nbx + ycounter)) = cntrBsites(hexsites(1,nbx + ycounter)) + 1
                        bsitesbonds(cntrBsites(hexsites(1,nbx + ycounter)),hexsites(1,nbx + ycounter)) = cntrB
                        cntrBsites(hexsites(2,nbx + ycounter)) = cntrBsites(hexsites(2,nbx + ycounter)) + 1
                        bsitesbonds(cntrBsites(hexsites(2,nbx + ycounter)),hexsites(2,nbx + ycounter)) = -cntrB
                    else if( pattern == 'BA' .and. modulo(x1, 2) == 0 ) then 
                        cntrB = cntrB + 1
                        blatt(1, cntrB) = hexsites(1,nbx + ycounter)
                        blatt(2, cntrB) = hexsites(2,nbx + ycounter)
                        blatt(3, cntrB) = - 1 !QAH phase 
                        blatt(4, cntrB) = nbx + ycounter !Bond number
                        blatt(5, cntrB) = 3 !Lattice vector 
                        phase(nbx + ycounter) = - 1
                        cntrBsites(hexsites(1,nbx + ycounter)) = cntrBsites(hexsites(1,nbx + ycounter)) + 1
                        bsitesbonds(cntrBsites(hexsites(1,nbx + ycounter)),hexsites(1,nbx + ycounter)) = cntrB
                        cntrBsites(hexsites(2,nbx + ycounter)) = cntrBsites(hexsites(2,nbx + ycounter)) + 1
                        bsitesbonds(cntrBsites(hexsites(2,nbx + ycounter)),hexsites(2,nbx + ycounter)) = -cntrB
                    else if( pattern == 'BA' .and. modulo(x1, 2) == 1 ) then 
                        cntrA = cntrA + 1              
                        alatt(1, cntrA) = hexsites(1,nbx + ycounter)
                        alatt(2, cntrA) = hexsites(2,nbx + ycounter)  
                        alatt(3, cntrA) = 1 !QAH phase       
                        alatt(4, cntrA) = nbx + ycounter !Bond number     
                        alatt(5, cntrA) = 3 !Lattice vector        
                        phase(nbx + ycounter) = 1
                        cntrAsites(hexsites(1,nbx + ycounter)) = cntrAsites(hexsites(1,nbx + ycounter)) + 1
                        asitesbonds(cntrAsites(hexsites(1,nbx + ycounter)),hexsites(1,nbx + ycounter)) = cntrA
                        cntrAsites(hexsites(2,nbx + ycounter)) = cntrAsites(hexsites(2,nbx + ycounter)) + 1
                        asitesbonds(cntrAsites(hexsites(2,nbx + ycounter)),hexsites(2,nbx + ycounter)) = -cntrA
                    end if 
                    xy(1, nbx + ycounter) = x1
                    xy(2, nbx + ycounter) = y1
                    xy(3, nbx + ycounter) = x2
                    xy(4, nbx + ycounter) = y2
                    
                    ! phase(nbx + ycounter) = - 1
                    60 continue
                    !if(y0 == 0 .and. x1 == nx - 1 .and. bc == 'o') goto 65

                    !Second vertical bond: Up left (lattice vector v2)
                    x2 = x1 !x-coordinate of second site
                    y2 = modulo(y1 + 1, ny) !y-coordinate of second site
                    if (y2 == y1) goto 65 !debug For 1D case no y-hopping
                    ycounter = ycounter + 1
                    hexsites(1,nbx + ycounter) = 1 + x1 + y1 * nx
                    hexsites(2,nbx + ycounter) = 1 + x2 + y2 * nx
                    latticevecs(nbx + ycounter)  = 2
                    ytcounter = ytcounter + 1 
                    ytransl(1, ytcounter) = hexsites(1,nbx + ycounter)
                    ytransl(2, ytcounter) = hexsites(2,nbx + ycounter)
                    ! print*, 'v2sites', hexsites(1,nbx + ycounter), hexsites(2,nbx + ycounter)                                                
                    if( pattern == 'AB' .and. modulo(x1, 2) == 0 ) then 
                        cntrA = cntrA + 1
                        alatt(1, cntrA) = hexsites(1,nbx + ycounter)
                        alatt(2, cntrA) = hexsites(2,nbx + ycounter)   
                        alatt(3, cntrA) = 1 !QAH phase     
                        alatt(4, cntrA) = nbx + ycounter !Bond number       
                        alatt(5, cntrA) = 2 !Lattice vector        
                        phase(nbx + ycounter) = 1      
                        cntrAsites(hexsites(1,nbx + ycounter)) = cntrAsites(hexsites(1,nbx + ycounter)) + 1
                        asitesbonds(cntrAsites(hexsites(1,nbx + ycounter)),hexsites(1,nbx + ycounter)) = cntrA
                        cntrAsites(hexsites(2,nbx + ycounter)) = cntrAsites(hexsites(2,nbx + ycounter)) + 1
                        asitesbonds(cntrAsites(hexsites(2,nbx + ycounter)),hexsites(2,nbx + ycounter)) = -cntrA
                    else if( pattern == 'AB' .and. modulo(x1, 2) == 1 ) then 
                        cntrB = cntrB + 1
                        blatt(1, cntrB) = hexsites(1,nbx + ycounter)
                        blatt(2, cntrB) = hexsites(2,nbx + ycounter)  
                        blatt(3, cntrB) = - 1 !QAH phase     
                        blatt(4, cntrB) = nbx + ycounter !Bond number    
                        blatt(5, cntrB) = 2 !Lattice vector            
                        phase(nbx + ycounter) = - 1      
                        cntrBsites(hexsites(1,nbx + ycounter)) = cntrBsites(hexsites(1,nbx + ycounter)) + 1
                        bsitesbonds(cntrBsites(hexsites(1,nbx + ycounter)),hexsites(1,nbx + ycounter)) = cntrB
                        cntrBsites(hexsites(2,nbx + ycounter)) = cntrBsites(hexsites(2,nbx + ycounter)) + 1
                        bsitesbonds(cntrBsites(hexsites(2,nbx + ycounter)),hexsites(2,nbx + ycounter)) = -cntrB
                    else if( pattern == 'BA' .and. modulo(x1, 2) == 0 ) then 
                        cntrB = cntrB + 1
                        blatt(1, cntrB) = hexsites(1,nbx + ycounter)
                        blatt(2, cntrB) = hexsites(2,nbx + ycounter)
                        blatt(3, cntrB) = 1 !QAH phase 
                        blatt(4, cntrB) = nbx + ycounter !Bond number
                        blatt(5, cntrB) = 2 !Lattice vector            
                        phase(nbx + ycounter) = - 1
                        cntrBsites(hexsites(1,nbx + ycounter)) = cntrBsites(hexsites(1,nbx + ycounter)) + 1
                        bsitesbonds(cntrBsites(hexsites(1,nbx + ycounter)),hexsites(1,nbx + ycounter)) = cntrB
                        cntrBsites(hexsites(2,nbx + ycounter)) = cntrBsites(hexsites(2,nbx + ycounter)) + 1
                        bsitesbonds(cntrBsites(hexsites(2,nbx + ycounter)),hexsites(2,nbx + ycounter)) = -cntrB
                    else if( pattern == 'BA' .and. modulo(x1, 2) == 1 ) then 
                        cntrA = cntrA + 1
                        alatt(1, cntrA) = hexsites(1,nbx + ycounter)
                        alatt(2, cntrA) = hexsites(2,nbx + ycounter)      
                        alatt(3, cntrA) = 1 !QAH phase        
                        alatt(4, cntrA) = nbx + ycounter !Bond number     
                        alatt(5, cntrA) = 2 !Lattice vector        
                        phase(nbx + ycounter) = 1     
                        cntrAsites(hexsites(1,nbx + ycounter)) = cntrAsites(hexsites(1,nbx + ycounter)) + 1
                        asitesbonds(cntrAsites(hexsites(1,nbx + ycounter)),hexsites(1,nbx + ycounter)) = cntrA
                        cntrAsites(hexsites(2,nbx + ycounter)) = cntrAsites(hexsites(2,nbx + ycounter)) + 1
                        asitesbonds(cntrAsites(hexsites(2,nbx + ycounter)),hexsites(2,nbx + ycounter)) = -cntrA
                    end if                 
                    xy(1, nbx + ycounter) = x1
                    xy(2, nbx + ycounter) = y1
                    xy(3, nbx + ycounter) = x2
                    xy(4, nbx + ycounter) = y2
                    ! phase(nbx + ycounter) = 1
                end if
                65 continue

                !Horizontal bond (lattice vector v1)
                if (x1 == nx-1 .and. bc == 'o') cycle !No NNN xbond to the right for last site
                if (x1 == nx-2 .and. bc == 'o') cycle !No NNN xbond to the right for second-to-last site
                xcounter = xcounter + 1
                x2 = modulo(x1 + 2, 2*ucx)
                y2 = y1
                hexsites(1,xcounter) = 1 + x1 + y1 * 2*ucx
                hexsites(2,xcounter) = 1 + x2 + y2 * 2*ucx
                latticevecs(xcounter)  = 1
                xtransl(1, xcounter) = hexsites(1,xcounter)
                xtransl(2, xcounter) = hexsites(2,xcounter)
                ! print*, 'v1sites', hexsites(1,xcounter), hexsites(2,xcounter)
                if( pattern == 'AB' .and. modulo(x1, 2) == 0 ) then 
                    cntrA = cntrA + 1
                    alatt(1, cntrA) = hexsites(1,xcounter)
                    alatt(2, cntrA) = hexsites(2,xcounter)      
                    alatt(3, cntrA) = 1 !QAH phase               
                    alatt(4, cntrA) = xcounter !Bond number   
                    alatt(5, cntrA) = 1 !Lattice vector                
                    phase(xcounter) = 1
                    cntrAsites(hexsites(1,xcounter)) = cntrAsites(hexsites(1,xcounter)) + 1
                    asitesbonds(cntrAsites(hexsites(1,xcounter)),hexsites(1,xcounter)) = cntrA
                    cntrAsites(hexsites(2,xcounter)) = cntrAsites(hexsites(2,xcounter)) + 1
                    asitesbonds(cntrAsites(hexsites(2,xcounter)),hexsites(2,xcounter)) = -cntrA
                else if( pattern == 'AB' .and. modulo(x1, 2) == 1 ) then 
                    cntrB = cntrB + 1
                    blatt(1, cntrB) = hexsites(1,xcounter)
                    blatt(2, cntrB) = hexsites(2,xcounter)
                    blatt(3, cntrB) = - 1 !QAH phase  
                    blatt(4, cntrB) = xcounter !Bond number    
                    blatt(5, cntrB) = 1 !Lattice vector            
                    phase(xcounter) = - 1              
                    cntrBsites(hexsites(1,xcounter)) = cntrBsites(hexsites(1,xcounter)) + 1
                    bsitesbonds(cntrBsites(hexsites(1,xcounter)),hexsites(1,xcounter)) = cntrB
                    cntrBsites(hexsites(2,xcounter)) = cntrBsites(hexsites(2,xcounter)) + 1
                    bsitesbonds(cntrBsites(hexsites(2,xcounter)),hexsites(2,xcounter)) = -cntrB
                else if( pattern == 'BA' .and. modulo(x1, 2) == 0 ) then 
                    cntrB = cntrB + 1
                    blatt(1, cntrB) = hexsites(1,xcounter)
                    blatt(2, cntrB) = hexsites(2,xcounter)
                    blatt(3, cntrB) = - 1 !QAH phase
                    blatt(4, cntrB) = xcounter !Bond number       
                    blatt(5, cntrB) = 1 !Lattice vector             
                    phase(xcounter) = - 1 
                    cntrBsites(hexsites(1,xcounter)) = cntrBsites(hexsites(1,xcounter)) + 1
                    bsitesbonds(cntrBsites(hexsites(1,xcounter)),hexsites(1,xcounter)) = cntrB
                    cntrBsites(hexsites(2,xcounter)) = cntrBsites(hexsites(2,xcounter)) + 1
                    bsitesbonds(cntrBsites(hexsites(2,xcounter)),hexsites(2,xcounter)) = -cntrB
                else if( pattern == 'BA' .and. modulo(x1, 2) == 1 ) then 
                    cntrA = cntrA + 1
                    alatt(1, cntrA) = hexsites(1,xcounter)
                    alatt(2, cntrA) = hexsites(2,xcounter)           
                    alatt(3, cntrA) = 1 !QAH phase           
                    alatt(4, cntrA) = xcounter !Bond number    
                    alatt(5, cntrA) = 1 !Lattice vector                
                    phase(xcounter) = 1
                    cntrAsites(hexsites(1,xcounter)) = cntrAsites(hexsites(1,xcounter)) + 1
                    asitesbonds(cntrAsites(hexsites(1,xcounter)),hexsites(1,xcounter)) = cntrA
                    cntrAsites(hexsites(2,xcounter)) = cntrAsites(hexsites(2,xcounter)) + 1
                    asitesbonds(cntrAsites(hexsites(2,xcounter)),hexsites(2,xcounter)) = -cntrA
                end if    
                xy(1, xcounter) = x1
                xy(2, xcounter) = y1
                xy(3, xcounter) = x2
                xy(4, xcounter) = y2
                sitecoord(1, hexsites(1,xcounter)) = x1 
                sitecoord(2, hexsites(1,xcounter)) = y1 
                ! phase(xcounter) = 1
            end do
        end do

        allocate( alattice(5,cntrA) )
        allocate( blattice(5,cntrB) )
        alattice(1, 1:cntrA) = alatt(1, 1:cntrA)
        alattice(2, 1:cntrA) = alatt(2, 1:cntrA)
        alattice(3, 1:cntrA) = alatt(3, 1:cntrA)
        alattice(4, 1:cntrA) = alatt(4, 1:cntrA)
        alattice(5, 1:cntrA) = alatt(5, 1:cntrA)
        blattice(1, 1:cntrB) = blatt(1, 1:cntrB)
        blattice(2, 1:cntrB) = blatt(2, 1:cntrB)
        blattice(3, 1:cntrB) = blatt(3, 1:cntrB)
        blattice(4, 1:cntrB) = blatt(4, 1:cntrB)
        blattice(5, 1:cntrB) = blatt(5, 1:cntrB)
        deallocate(alatt)
        deallocate(blatt)
        write (name,"('L=',i0,'ucx=',i0,'ucy=',i0,'pat=',a2,'BC=',a,'.dat')") nn,ucx,ucy,pattern,bc
        name = trim_name(name)
        name=dir//"A_lattice_"//name
        open(unit=32, file=name)
        do s = 1, cntrA
            write(32,'(i0,x,i0,x,i0,x,i0,x,i0)') alattice(1,s), alattice(2,s), alattice(3,s), alattice(4,s), alattice(5,s)
        end do
        close(32)


        write (name,"('L=',i0,'ucx=',i0,'ucy=',i0,'pat=',a2,'BC=',a,'.dat')") nn,ucx,ucy,pattern,bc
        name = trim_name(name)
        name=dir//"B_lattice_"//name
        open(unit=32, file=name)
        do s = 1, cntrB
            write(32,'(i0,x,i0,x,i0,x,i0,x,i0)') blattice(1,s), blattice(2,s), blattice(3,s), blattice(4,s), blattice(5,s)
        end do
        close(32)

        nbonds = xcounter + ycounter

        ! allocate(phase(nbonds))
        ! phase(1:xcounter) = phase(1:xcounter)
        ! phase(xcounter + 1:xcounter +  ycounter) = phase(nbx + 1:nbx +  ycounter)

        write (name,"('L=',i0,'ucx=',i0,'ucy=',i0,'pat=',a2,'BC=',a,'.dat')") nn,ucx,ucy,pattern,bc
        name = trim_name(name)
        name=dir//"NNN_lattice_"//name
        open(unit=32, file=name)
        do s = 1, nbonds
            write(32,*) hexsites(1,s), hexsites(2,s)      
        end do

        close(32)

        write (name,"('L=',i0,'ucx=',i0,'ucy=',i0,'pat=',a2,'BC=',a,'.dat')") nn,ucx,ucy,pattern,bc
        name = trim_name(name)
        name=dir//"NNN_xy_"//name
        open(unit=32, file=name)
        do s = 1, nbonds
            write(32,'(4(i0,x))') xy(1,s), xy(2,s), xy(3,s), xy(4,s)
        end do
        close(32)

        ! print*, 'Finished honeycomb nnn-lattice'

    end subroutine honeycomb_nnn




end module lattice 