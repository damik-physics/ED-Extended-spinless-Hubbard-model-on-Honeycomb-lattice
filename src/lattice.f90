module lattice
    use types
    use functions 
    use parameters
    implicit none
    

    contains 
    subroutine define_lattice(dir, par, geo)
    ! subroutine define_lattice(dir, tilted, ucx, ucy, nnnBonds, geo%nnBonds, par%bc, pattern, cluster, bsites, geo%hexsites, geo%latticevecs, alattice, geo%blattice, xyA, xyB, geo%BsitesBonds, geo%BsitesBonds, geo%cntrA,geo%cntrB, geo%nHel, geo%tilt, phase, xy, geo%xtransl, geo%ytransl, reflections, nnnVec)
    ! call define_lattice(outdir, tilted, ucx, ucy, nnnBonds, geo%nnBonds, bc, pattern, cluster, bsites, geo%hexsites, geo%latticevecs, alattice, geo%blattice, xyA, xyB, geo%BsitesBonds, BsitesBonds, geo%cntrA,geo%cntrB, geo%nHel, geo%tilt, phases, xy, geo%xtransl, geo%ytransl, reflections, nnnVec)
        implicit none
    
        ! integer, intent(in) :: ucx, ucy, tilted
        ! character(len=*), intent(in) :: dir, bc, pattern, cluster
        ! integer, intent(out) :: nnnBonds, geo%nnBonds, geo%cntrA,geo%cntrB, geo%nHel, geo%tilt 
        ! integer, allocatable, intent(out) :: bsites(:,:)
        ! integer, allocatable, intent(out) :: geo%hexsites(:,:), geo%phases(:), geo%xy(:,:), geo%latticevecs(:)
        ! integer, allocatable, intent(out) :: alattice(:,:), geo%blattice(:,:), xyA(:,:), xyB(:,:)
        ! integer, allocatable, intent(out) :: geo%BsitesBonds(:,:), BsitesBonds(:,:)
        ! integer, allocatable, intent(out) :: geo%xtransl(:,:), geo%ytransl(:,:)
        ! integer, allocatable, intent(out) :: reflections(:)
        ! double precision, allocatable, intent(out) :: nnnVec(:,:)
        character(len=*), intent(inout) :: dir
        type(sim_params), intent(inout) :: par
        type(geometry), intent(inout) :: geo
        ! integer, allocatable :: geo%sitecoord(:,:)
        integer :: i!, geo%nUC
        character :: name*512

        if(par%tilted == 0) then     

            call honeycomb(dir, par, geo)
            call honeycomb_nnn(dir, par, geo)
            call coordinates(dir, par, geo)
            call reflection(par, geo)        

            ! call honeycomb(dir, ucx, ucy, bc, pattern, nnnBonds, bsites)
            ! call honeycomb_nnn(dir, ucx, ucy, bc, pattern, geo%nnBonds, geo%hexsites, geo%latticevecs, alattice, geo%blattice, geo%BsitesBonds, BsitesBonds, geo%cntrA,geo%cntrB, phase, xy, geo%xtransl, geo%ytransl, geo%sitecoord, nnnVec)
            ! call coordinates(dir, 20, ucx, ucy, geo%nnBonds, geo%cntrA,geo%cntrB, xy, alattice, geo%blattice, bc, pattern, xyA, xyB)
            ! call reflection(par%ucx, par%ucy, sitecoord, geo%reflections)        
            print*, 'Honeycomb lattice created'
        else if(par%tilted == 1) then 
            if(par%cluster == '18A') then 
                geo%nUC  = 9
                geo%nHel = 3
                geo%tilt = 2
            else if(par%cluster == '18B') then 
                geo%nUC  = 9
                geo%nHel = 1
                geo%tilt = 4
            else if(par%cluster == '18C') then 
                geo%nUC  = 9
                geo%nHel = 1
                geo%tilt = 8
            else if(par%cluster == '20A') then 
                geo%nUC  = 10
                geo%nHel = 1
                geo%tilt = 6
            else if(par%cluster == '20B') then 
                geo%nUC  = 10
                geo%nHel = 5
                geo%tilt = 2
            else if(par%cluster == '24C') then 
                geo%nUC  = 12
                geo%nHel = 1
                geo%tilt = 8
            else if(par%cluster == '24D') then 
                geo%nUC  = 12
                geo%nHel = 2
                geo%tilt = 2
            else if(par%cluster == '30A') then 
                geo%nUC  = 15
                geo%nHel = 1
                geo%tilt = 20
            else if(par%cluster == '16A') then 
                geo%nUC  = 8
                geo%nHel = 2
                geo%tilt = 2
            else if(par%cluster == '16B') then 
                geo%nUC  = 8
                geo%nHel = 1
                geo%tilt = 4
            else if(par%cluster == 'L12') then 
                geo%nUC  = 6
                geo%nHel = 1
                geo%tilt = 2
            else if(par%cluster == 'L18') then 
                geo%nUC  = 9
                geo%nHel = 3
                geo%tilt = 2
            else if(par%cluster == 'L24') then 
                geo%nUC  = 12
                geo%nHel = 2
                geo%tilt = 2
            else if(par%cluster == 'L24B') then 
                geo%nUC  = 12
                geo%nHel = 1
                geo%tilt = 4
            else if(par%cluster == 'L26') then 
                geo%nUC  = 13
                geo%nHel = 1
                geo%tilt = 6
            else if(par%cluster == 'L28') then 
                geo%nUC  = 14
                geo%nHel = 2
                geo%tilt = 2
            else if(par%cluster == 'L32') then 
                geo%nUC  = 16
                geo%nHel = 4
                geo%tilt = 2
            else if(par%cluster == 'L34') then 
                geo%nUC  = 17
                geo%nHel = 1
                geo%tilt = 20
            else if(par%cluster == 'L36') then 
                geo%nUC  = 18
                geo%nHel = 3
                geo%tilt = 2
            else if(par%cluster == 'L38') then 
                geo%nUC  = 19
                geo%nHel = 1
                geo%tilt = 22
            else if(par%cluster == 'L42') then 
                geo%nUC  = 21
                geo%nHel = 1
                geo%tilt = 8
            else if(par%cluster == 'L54') then 
                geo%nUC  = 27
                geo%nHel = 3
                geo%tilt = 2
            end if
            
            write(name,"('Parameters_Cluster=',a,'_BC=',a,'.dat')") par%cluster, par%bc
            name = trim(dir//name)
            open(30,file=name)
            write(30,*) 'N_helices', geo%nHel
            write(30,*) 'Tilt', geo%tilt
            close(30)

            call build_cluster(geo)
            ! call build_cluster(nUC, geo%nHel, geo%tilt, geo%cntrA,geo%cntrB, nnnBonds, nnnVec, alattice, geo%blattice, geo%hexsites, bsites, geo%xtransl, geo%ytransl)
            geo%nnnBonds = geo%cntrA + geo%cntrB 
            
                 
            write(*,'(a,a,a)') ' Honeycomb lattice cluster ',par%cluster,' created' 
            if(allocated(geo%xyA)) deallocate(geo%xyA)
            if(allocated(geo%xyB)) deallocate(geo%xyB)
            allocate(geo%xyA(4, geo%cntrA))
            allocate(geo%xyB(4, geo%cntrB))
            geo%xyA = 0 
            geo%xyB = 0 

            write(name,"('Alattice_Cluster=',a,'_BC=',a,'.dat')") par%cluster, par%bc
            name = trim(dir//name)
            open(30,file=name)
            do i = 1, geo%cntrA
                write(30,*) geo%alattice(1, i), geo%alattice(2, i), geo%alattice(3, i), geo%alattice(4, i), geo%alattice(5, i)
            end do 
            close(30)
            write(name,"('geo%Blattice_Cluster=',a,'_BC=',a,'.dat')") par%cluster, par%bc
            name = trim(dir//name)
            open(30,file=name)
            do i = 1, geo%cntrB
                write(30,*) geo%blattice(1, i), geo%blattice(2, i), geo%blattice(3, i), geo%blattice(4, i), geo%blattice(5, i)
            end do 
            close(30)
            write(name,"('NNlattice_Cluster=',a,'_BC=',a,'.dat')") par%cluster, par%bc
            name = trim(dir//name)
            open(30,file=name)
            do i = 1, geo%nnnBonds
                write(30,*) geo%bsites(1, i), geo%bsites(2, i)
            end do 
            close(30)
            write(name,"('NNNlattice_Cluster=',a,'_BC=',a,'.dat')") par%cluster, par%bc
            name = trim(dir//name)
            open(30,file=name)
            do i = 1, geo%nnnBonds
                write(30,*) geo%hexsites(1, i), geo%hexsites(2, i)
            end do 
            close(30)
        end if 

    end subroutine define_lattice

    subroutine build_cluster(geo)
    ! subroutine build_cluster(nUC, geo%nHel, geo%tilt, nabonds, cntrA, geo%nnBonds, nnnVec, Alattice, geo%Blattice, geo%nnBonds, nnnBonds, geo%xtransl, geo%ytransl)
        implicit none 
        type(geometry), intent(inout) :: geo
        ! integer, intent(in) :: nUC, geo%nHel, geo%tilt 

        ! integer, intent(out) :: nabonds, cntrA, geo%nnBonds 
        ! integer, allocatable, intent(out) :: Alattice(:,:), geo%Blattice(:,:), geo%nnBonds(:,:), geo%bsites(:,:), geo%xtransl(:,:), geo%ytransl(:,:)
        ! double precision, allocatable, intent(out) :: nnnVec(:,:)

        integer :: nBondsA, nBondsB, site1, site2, sites, nS
        integer :: cntr, xcntr, ycntr
        integer :: h, i

         
        if(allocated(geo%nnnVec)) deallocate(geo%nnnVec)
        allocate(geo%nnnVec(2,3))
        geo%nnnVec(1:2, 1) = (/0.d0, sqrt(3.d0)/) !a1 |
        geo%nnnVec(1:2, 2) = 0.5 * (/3.d0, -sqrt(3.d0)/) !a2 \
        geo%nnnVec(1:2, 3) = 0.5 * (/-3.d0, -sqrt(3.d0)/) !a3 /

        sites = geo%nUC * 2
        nS = sites / geo%nHel !Number of sites per helix in a2 direction 
        nBondsA = geo%nUC * 3 
        geo%cntrB = geo%nUC * 3 
        geo%nnnBonds = geo%nUC * 3 
        if(allocated(geo%Alattice)) deallocate(geo%Alattice)
        if(allocated(geo%Blattice)) deallocate(geo%Blattice)
        if(allocated(geo%nnnVec)) deallocate(geo%nnnVec)
        if(allocated(geo%bsites)) deallocate(geo%bsites)
        if(allocated(geo%xtransl)) deallocate(geo%xtransl)
        if(allocated(geo%ytransl)) deallocate(geo%ytransl)
        allocate(geo%Alattice(5, nBondsA))
        allocate(geo%Blattice(5, nBondsB))
        allocate(geo%bsites(5, nBondsA + nBondsB))
        allocate(geo%hexsites(2, geo%nnnBonds))
        allocate(geo%xtransl(2, sites))
        allocate(geo%ytransl(2, sites))

        geo%cntrA = 0 
        geo%cntrB = 0
        xcntr = 0 
        ycntr = 0
        
        do h = 0, geo%nHel - 1 !Loop over all helices
            do i = 0, nS - 1 !Loop over all sites per helix
                site1 = i + h * nS + 1
                if(modulo(i, 2) == 0) then 
                    
                    !a1 direction 
                    ycntr = ycntr + 1 
                    geo%cntrA = geo%cntrA + 1 
                    site2 = modulo(i + geo%tilt, nS) + modulo(h+1, geo%nHel) * nS + 1
                    geo%Alattice(1, geo%cntrA) = site1
                    geo%Alattice(2, geo%cntrA) = site2 
                    geo%ytransl(1, ycntr) = site1
                    geo%ytransl(2, ycntr) = site2
                    geo%Alattice(3, geo%cntrA) = 1
                    geo%Alattice(4, geo%cntrA) = geo%cntrA + geo%cntrB 
                    geo%Alattice(5, geo%cntrA) = 1
                        
                    !a2 direction 
                    xcntr = xcntr + 1
                    geo%cntrA = geo%cntrA + 1 
                    site2 = modulo(i + 2, nS) + h * nS + 1
                    geo%Alattice(1, geo%cntrA) = site1   
                    geo%Alattice(2, geo%cntrA) = site2 
                    geo%Alattice(3, geo%cntrA) = 1
                    geo%Alattice(4, geo%cntrA) = geo%cntrA + geo%cntrB 
                    geo%Alattice(5, geo%cntrA) = 2
                    geo%xtransl(1, xcntr) = site1
                    geo%xtransl(2, xcntr) = site2

                    !a3 direction 
                    geo%cntrA = geo%cntrA + 1 
                    site2 = modulo(modulo(i - geo%tilt, nS) - 2, nS) + modulo(h-1, geo%nHel) * nS + 1 
                    geo%Alattice(1, geo%cntrA) = site1 
                    geo%Alattice(2, geo%cntrA) = site2 
                    ! geo%Alattice(2, geo%cntrA) = modulo(site1 - 1 - nUp - 2, sites) + 1!? mod(..., sites)
                    geo%Alattice(3, geo%cntrA) = 1
                    geo%Alattice(4, geo%cntrA) = geo%cntrA + geo%cntrB 
                    geo%Alattice(5, geo%cntrA) = 3

                else if( modulo(i, 2) == 1 ) then  
                    
                    !a1 direction 
                    ycntr = ycntr + 1
                   geo%cntrB =geo%cntrB + 1 
                    site2 = modulo(i + geo%tilt, nS) + modulo(h+1, geo%nHel) * nS + 1
                    geo%Blattice(1,geo%cntrB) = site1
                    geo%Blattice(2,geo%cntrB) = site2 
                    ! if( i + nUp .le. nS - i) then 
                    !     ! site2 = modulo(site + nUp + nS, sites) !? mod(..., sites)
                    !     site2 = modulo(site1 - 1 + nUp + nS, sites) + 1 !? mod(..., sites)
                    !     ! site2 = modulo(site1 - 1 + nUp + h * nS, sites) + 1 !? mod(..., sites)
                    !     geo%Blattice(2,geo%cntrB) = site2 
                    ! else 
                    !     site2 = modulo(site1 - 1 + nUp, sites) + 1
                    !     geo%Blattice(2,geo%cntrB) = site2 
                    ! end if
                    geo%Blattice(3,geo%cntrB) = -1
                    geo%Blattice(4,geo%cntrB) = geo%cntrA + geo%cntrB 
                    geo%Blattice(5,geo%cntrB) = 1
                    geo%ytransl(1, ycntr) = site1
                    geo%ytransl(2, ycntr) = site2
                        
                    !a2 direction 
                    xcntr = xcntr + 1
                   geo%cntrB =geo%cntrB + 1 
                    site2 = modulo(i + 2, nS) + h * nS + 1
                    geo%Blattice(1,geo%cntrB) = site1   
                    geo%Blattice(2,geo%cntrB) = site2 
                    geo%Blattice(3,geo%cntrB) = -1
                    geo%Blattice(4,geo%cntrB) = geo%cntrA + geo%cntrB 
                    geo%Blattice(5,geo%cntrB) = 2
                    geo%xtransl(1, xcntr) = site1
                    geo%xtransl(2, xcntr) = site2

                    !a3 direction 
                    geo%cntrB = geo%cntrB + 1 
                    site2 = modulo(modulo(i - geo%tilt, nS) - 2, nS) + modulo(h-1, geo%nHel) * nS + 1 
                    ! site2 = modulo(site1 - 1 - nUp - 2, sites) + 1!? mod(..., sites)
                    geo%Blattice(1,geo%cntrB) = site1 
                    geo%Blattice(2,geo%cntrB) = site2 
                    geo%Blattice(3,geo%cntrB) = -1
                    geo%Blattice(4,geo%cntrB) = geo%cntrA + geo%cntrB 
                    geo%Blattice(5,geo%cntrB) = 3

                end if 
                if (geo%cntrA > nBondsA) print*,'Build cluster: Too many A bonds'
                if (geo%cntrB > nBondsB) print*,'Build cluster: Too many B bonds'
            end do 
        end do 
        
        ! Store next-nearest neighbor bonds
        do i = 1, geo%cntrA
            geo%hexsites(1, 2 *(i-1) + 1) = geo%Alattice(1, i)
            geo%hexsites(2, 2 *(i-1) + 1) = geo%Alattice(2, i)
            geo%hexsites(1, 2 *i) = geo%Blattice(1, i)
            geo%hexsites(2, 2 *i) = geo%Blattice(2, i)
        
            if(geo%Alattice(1, i) > sites) print*,'Build cluster: A site ',geo%Alattice(1, i), 'is larger than maximum number of sites.'
            if(geo%Alattice(2, i) > sites) print*,'Build cluster: A site ',geo%Alattice(2, i), 'is larger than maximum number of sites.'
            if(geo%Blattice(1, i) > sites) print*,'Build cluster: B site ',geo%Blattice(1, i), 'is larger than maximum number of sites.'
            if(geo%Blattice(2, i) > sites) print*,'Build cluster: B site ',geo%Blattice(2, i), 'is larger than maximum number of sites.'
        end do 

        ! Store nearest neighbor bonds 
        cntr = 0
        do i = 1, geo%cntrA, 3
            !e1
            cntr = cntr + 1
            geo%bsites(1, cntr) = geo%Alattice(1, i)
            geo%bsites(2, cntr) = geo%Blattice(1, i)
            !e2
            cntr = cntr + 1
            geo%bsites(1, cntr) = geo%Alattice(1, i)
            geo%bsites(2, cntr) = geo%Blattice(2, 3 * geo%Blattice(2, i + 2) / 2 - 2)
            !e3
            cntr = cntr + 1
            geo%bsites(1, cntr) = geo%Alattice(1, i)
            geo%bsites(2, cntr) = geo%Blattice(2, i + 2)
        end do 
        
        if(cntr > geo%nnBonds) print*,'Build cluster: Too many NN bonds.'   

    end subroutine build_cluster

    subroutine reflection(par, geo)
    ! subroutine reflection(ucx, ucy, sitecoord, reflections)
        implicit none 
        type(sim_params), intent(inout) :: par
        type(geometry), intent(inout) :: geo
        ! integer, intent(in) :: ucx, ucy
        ! integer, intent(in) :: sitecoord(2, *)
        
        ! integer, allocatable, intent(out) :: reflections(:)
        
        integer :: sites = 0 
        integer :: i = 0, j = 0, xr = 0, yr = 0 
        integer, allocatable :: uccoord(:,:)

        sites = 2 * par%ucx * par%ucy 
        if(allocated(uccoord)) deallocate(uccoord)
        allocate(uccoord(3, sites))

        if(allocated(geo%reflections)) deallocate(geo%reflections)
        allocate(geo%reflections(sites))
        
        do j = 1, sites 
            uccoord(1, j) = geo%sitecoord(1, j)/2 - 1
            uccoord(2, j) = geo%sitecoord(2, j) - 1
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
                    geo%reflections(j) = i
                    exit 
                end if 
            end do 
        end do 
        return 

    end subroutine reflection

    subroutine coordinates(dir, par, geo)!unit, ucx, ucy, geo%nnBonds, abonds, bbonds, xy, geo%alattice, geo%blattice, bc, pattern, xyA, xyB)
        implicit none 
        ! integer, intent(in) :: unit, ucx, ucy, geo%nnBonds, abonds, bbonds
        ! integer, intent(in) :: geo%xy(4, geo%nnBonds), geo%alattice(5, abonds), geo%blattice(5, bbonds)
        character(len=*), intent(in) :: dir!, bc, pattern
        ! integer, allocatable, intent(out) :: xyA(:,:), xyB(:,:)
        type(sim_params), intent(inout) :: par
        type(geometry), intent(inout) :: geo
        type(system_params) :: sys

        integer   :: j = 0 
        character :: filename*256
        character :: par_char*256

        write(par_char,"('L=',i0,'ucx=',i0,'ucy=',i0,'BC=',a,'_pat=',a2'.dat')") 2 *par%ucx*par%ucy,par%ucx,par%ucy, par%bc,pattern
        par_char = trim(par_char)
        if(allocated(geo%xyA)) deallocate(geo%xyA)
        if(allocated(geo%xyB)) deallocate(geo%xyB)
        allocate(geo%xyA(4, geo%cntrA))
        allocate(geo%xyB(4, geo%cntrB))
        filename = dir//"A_coordinates_"//par_char
        filename = trim(filename)
        open(sys%unit,file = filename)
        do j = 1, geo%cntrA             
            geo%xyA(1,j) = geo%xy(1,geo%alattice(4,j))
            geo%xyA(2,j) = geo%xy(2,geo%alattice(4,j))
            geo%xyA(3,j) = geo%xy(3,geo%alattice(4,j))
            geo%xyA(4,j) = geo%xy(4,geo%alattice(4,j))
            write(sys%unit,'(4(I0,x))') geo%xyA(1,j), geo%xyA(2,j), geo%xyA(3,j), geo%xyA(4,j) 
        end do   
        close(sys%unit)

        filename = dir//"B_coordinates_"//par_char
        filename = trim(filename)
        open(sys%unit,file = filename)
        do j = 1, geo%cntrB             
            geo%xyB(1,j) = geo%xy(1,geo%blattice(4,j))
            geo%xyB(2,j) = geo%xy(2,geo%blattice(4,j))
            geo%xyB(3,j) = geo%xy(3,geo%blattice(4,j))
            geo%xyB(4,j) = geo%xy(4,geo%blattice(4,j))
            write(sys%unit,'(4(I0,x))') geo%xyB(1,j), geo%xyB(2,j), geo%xyB(3,j), geo%xyB(4,j) 
        end do 
        close(sys%unit)

        return 
    end subroutine coordinates

    subroutine honeycomb(dir, par, geo)
    ! subroutine honeycomb(dir, ucx, ucy, bc, pattern, geo%nnBonds, geo%bsites)
        use types
        use parameters, only: pattern
        !Subroutine obtained from: http://physics.bu.edu/~sandvik/vietri/sse/ssebasic.f90
        implicit none

        ! integer, intent(in) :: ucx, ucy
        character(len=*), intent(in) :: dir!, bc, pattern

        ! integer, intent(out) :: geo%nnBonds
        ! integer, allocatable, intent(out) :: geo%bsites(:,:)
        type(sim_params), intent(in) :: par ! Simulation parameters
        type(geometry), intent(inout) :: geo ! geometry parameters

        integer, allocatable :: bsites_temp(:,:)
        integer :: nn, nx, ny       ! number of sites
        integer :: y0
        integer :: xcounter, ycounter
        integer :: nbx, nby
        integer :: s, x1, x2, y1, y2
        integer, allocatable :: xy(:,:)    
        character :: name*512

        !d1 = (0,1)
        !d2 = (-Sqrt[3]/2,-1/2)
        !d3 = (+Sqrt[3]/2,-1/2)
        !First site on A lattice (pattern AB)
        xcounter = 0
        ycounter = 0

        nn = 2 * par%ucx * par%ucy !Number of lattice sites
        nx = 2 * par%ucx
        ny = par%ucy

        if(par%bc == 'o') then !nx = number of bonds per x-layer
            nbx = (nx - 1) * par%ucy !total number of x bonds
            nby = ((nx-1) + y0)/2 * (ny-1) !total number of y bonds
        else if(par%bc == 'p') then
            nbx = nx * par%ucy !total number of x bonds
            ! nby = int(nx/2) * ny !total number of y bonds
            nby = par%ucx * ny !total number of y bonds
        end if
        geo%nnnBonds = nbx + nby
        ! print*, 'Number of NN bonds', geo%nnBonds 

        if (allocated(bsites_temp)) deallocate(bsites_temp)
        allocate(bsites_temp(2,geo%nnnBonds))
        allocate(xy(4,geo%nnnBonds))

        do y1 = 0, ny - 1
            do x1 = 0, nx - 1
                if (mod(x1,2) == 0 .and. nbx + ycounter <= geo%nnnBonds) then !y-bonds d1 
                    if(y1 == ny - 1 .and. par%bc == 'o') goto 55
                    x2 = modulo(x1 + 1, nx) !x-coordinate of second site
                    y2 = modulo(y1 + 1, ny) !y-coordinate of second site
                    if (y2 == y1) goto 55 !debug For 1D case no y-hopping
                    ycounter = ycounter + 1 !Current y-bond = nbx + ycounter
                    bsites_temp(1,nbx + ycounter) = 1 + x1 + y1 * nx
                    bsites_temp(2,nbx + ycounter) = 1 + x2 + y2 * nx
                    ! print*, 'Ysites', bsites_temp(1,nbx + ycounter), bsites_temp(2,nbx + ycounter)
                    geo%xy(1, nbx + ycounter) = x1
                    geo%xy(2, nbx + ycounter) = y1
                    geo%xy(3, nbx + ycounter) = x2
                    geo%xy(4, nbx + ycounter) = y2                
                end if

                55 continue
                if (x1 == nx-1 .and. par%bc == 'o') cycle !x-bonds d3 
                xcounter = xcounter + 1
                x2 = modulo(x1 + 1, nx)
                y2 = y1
                bsites_temp(1,xcounter) = 1 + x1 + y1 * nx
                bsites_temp(2,xcounter) = 1 + x2 + y2 * nx
                ! print*, 'Xsites', bsites_temp(1,xcounter), bsites_temp(2,xcounter)                        
                geo%xy(1, xcounter) = x1
                geo%xy(2, xcounter) = y1
                geo%xy(3, xcounter) = x2
                geo%xy(4, xcounter) = y2            
            end do
        end do

        geo%nnnBonds = xcounter + ycounter
        allocate(geo%bsites(2,geo%nnnBonds))
        geo%bsites(1,1:xcounter) = bsites_temp(1,1:xcounter)
        geo%bsites(2,1:xcounter) = bsites_temp(2,1:xcounter)
        geo%bsites(1,xcounter + 1:xcounter  + ycounter) = bsites_temp(1,nbx + 1:nbx  + ycounter)
        geo%bsites(2,xcounter + 1:xcounter  + ycounter) = bsites_temp(2,nbx + 1:nbx  + ycounter)

        deallocate(bsites_temp)
        write(name,"('L=',i0,'ucx=',i0,'ucy=',i0,'pat=',a2,'BC=',a,'.dat')") geo%sites,par%ucx,par%ucy,pattern, par%bc
        name = trim(name)
        name=dir//"NN_lattice_"//name
        open(unit=31, file=name)
        do s = 1, geo%nnnBonds
            write(31,*) geo%bsites(1,s), geo%bsites(2,s)
        end do
        close(31)

        write(name,"('L=',i0,'ucx=',i0,'ucy=',i0,'pat=',a2,'BC=',a,'.dat')") nn,par%ucx,par%ucy,pattern, par%bc
        name = trim(name)
        name=dir//"NN_xy_"//name
        open(unit=32, file=name)
        do s = 1, geo%nnnBonds
            write(32,*) geo%xy(1,s), geo%xy(2,s), geo%xy(3,s), geo%xy(4,s)
        end do
        close(32)

        ! print*, 'Finished honeycomb lattice'

    end subroutine honeycomb

    subroutine honeycomb_nnn(dir, par, geo)
    ! subroutine honeycomb_nnn(dir, ucx, ucy, bc, pattern, geo%nnBonds, geo%hexsites, geo%latticevecs, alattice, geo%blattice, geo%BsitesBonds, BsitesBonds, geo%cntrA,geo%cntrB, phase, xy, geo%xtransl, geo%ytransl, sitecoord, nnnVec)
        !Subroutine obtained from: http://physics.bu.edu/~sandvik/vietri/sse/ssebasic.f90
        use parameters, only: pattern
        implicit none

        ! integer, intent(in) :: ucx, ucy
        character(len=*), intent(in) :: dir!, bc, pattern
        ! integer, intent(out) :: geo%nnBonds
        ! integer, intent(out) :: geo%cntrA
        ! integer, intent(out) ::geo%cntrB
        ! integer, allocatable, intent(out) :: geo%hexsites(:,:)
        ! integer, allocatable, intent(out) :: geo%latticevecs(:)
        ! integer, allocatable, intent(out) :: alattice(:,:)
        ! integer, allocatable, intent(out) :: geo%blattice(:,:)
        ! integer, allocatable, intent(out) :: geo%BsitesBonds(:,:)
        ! integer, allocatable, intent(out) :: BsitesBonds(:,:)
        ! integer, allocatable, intent(out) :: geo%phases(:)
        ! integer, allocatable, intent(out) :: geo%xy(:,:)
        ! integer, allocatable, intent(out) :: geo%xtransl(:,:)
        ! ! integer, allocatable, intent(out) :: geo%ytransl(:,:)
        ! integer, allocatable, intent(out) :: sitecoord(:,:)
        ! double precision, allocatable, intent(out) :: nnnVec(:,:)
        type(sim_params), intent(inout) :: par
        type(geometry), intent(inout) :: geo
        ! integer, allocatable :: geo%hexsites_temp(:,:), phase_temp(:)
        integer :: nn, nx, ny
        integer :: y0
        integer :: xcounter, ycounter, ytcounter
        integer :: nbx, nby
        integer :: s, x1, x2, y1, y2
        ! integer, allocatable :: geo%xy(:,:)
        integer, allocatable :: alatt(:,:)
        integer, allocatable :: blatt(:,:)
        integer, allocatable :: cntrAsites(:)
        integer, allocatable :: cntrBsites(:)
        
        character :: name*512

        !Definition of lattice vectors: 
        ! v1 = {sqrt(3), 0}
        ! v2 = 1/2 *{-sqrt(3), 3}
        ! v3 = 1/2 *{-sqrt(3), -3}
        if(allocated(geo%nnnVec)) deallocate(geo%nnnVec)
        allocate(geo%nnnVec(2, 3))
        geo%nnnVec(1:2,1) = (/sqrt(3.d0), 0.d0/)
        geo%nnnVec(1:2,2) = 0.5*(/-1*sqrt(3.d0), 3.d0/)
        geo%nnnVec(1:2,3) = 0.5*(/-1*sqrt(3.d0), -3.d0/)

        ! print*, 'Start creating Honeycomb sublattices'
        xcounter = 0
        ycounter = 0
        ytcounter = 0
        geo%cntrA = 0 
        geo%cntrB = 0 
        if(pattern == 'AB') then !Determines pattern in x-direction: ABAB... = 1; BABA... = 0
            y0 = 1
        else if(pattern == 'BA') then
            y0 = 0
        end if
        nn = 2 * par%ucx * par%ucy !Number of lattice sites
        nx = 2 * par%ucx !nx = number of bonds per x-layer (PBC)
        ny = par%ucy !Number of y layers

        if(par%bc == 'o') then
            nbx = 2 * (par%ucx - 1) * par%ucy !total number of x bonds: 2 = # sublattices
            nby = 2 * (par%ucy - 1) * ( 2 * (par%ucx-1) + 1 ) !total number of y bonds: 2 = # sublattices; ucx-1 = unit cells with 2 bonds/site; 1 = unit cell with 1 bond per site (first uc)
        else if(par%bc == 'p') then
            nbx = nx * ny !total number of x bonds
            nby = 2 * nx * ny !total number of y bonds
        end if
        geo%nnBonds = nbx + nby
        ! print*, 'Number of NNN bonds', geo%nnBonds 
        

        if(allocated(geo%hexsites)) deallocate(geo%hexsites)
        if(allocated(geo%latticevecs)) deallocate(geo%latticevecs)
        if(allocated(geo%phases)) deallocate(geo%phases)
        if(allocated(geo%xy)) deallocate(geo%xy)
        if(allocated(geo%AsitesBonds)) deallocate(geo%AsitesBonds)
        if(allocated(geo%BsitesBonds)) deallocate(geo%BsitesBonds)
        if(allocated(cntrAsites)) deallocate(cntrAsites)
        if(allocated(cntrBsites)) deallocate(cntrBsites)
        if(allocated(geo%xtransl)) deallocate(geo%xtransl)
        if(allocated(geo%ytransl)) deallocate(geo%ytransl)
        if(allocated(geo%nnnVec)) deallocate(geo%nnnVec)
        
        allocate(geo%xtransl(2, nn))
        allocate(geo%ytransl(2, nn))
        allocate(geo%hexsites(2,geo%nnBonds))
        allocate(geo%latticevecs(geo%nnBonds))
        allocate(alatt(5,geo%nnBonds))
        allocate(blatt(5,geo%nnBonds))
        allocate(geo%phases(geo%nnBonds))
        allocate(geo%xy(4, geo%nnBonds))
        allocate(geo%AsitesBonds(6, nn))
        allocate(geo%BsitesBonds(6, nn))
        allocate(cntrAsites(nn))
        allocate(cntrBsites(nn))
        allocate(geo%nnnVec(2, nn))


        geo%hexsites  = 0
        geo%latticevecs = 0 
        geo%phases       = 0
        cntrAsites  = 0
        cntrBsites  = 0
        geo%AsitesBonds = 0 
        geo%BsitesBonds = 0 
        geo%nnnVec   = 0 

        do y1 = 0, ny - 1
            do x1 = 0, nx - 1
                if(y1 == ny - 1 .and. par%bc == 'o') goto 65 !No y-bonds for last row at OBC 
                if (nbx + ycounter <= geo%nnBonds) then !y-bonds

                    !First vertical bond: Down left (lattice vector v3)
                    if((x1 == 0 .or. x1 == 1) .and. par%bc == 'o') goto 60 !No down-left NNN-bond for 1st and 2nd sites             
                    x2 = modulo(x1 - 2, nx) !x-coordinate of second site 
                    y2 = modulo(y1 - 1, ny) !y-coordinate of second site
                    if (y2 == y1) goto 60 !debug For 1D case no y-hopping
                    ycounter = ycounter + 1 !Current y-bond = nbx + ycounter
                    geo%hexsites(1,nbx + ycounter) = 1 + x1 + y1 * nx
                    geo%hexsites(2,nbx + ycounter) = 1 + x2 + y2 * nx     
                    geo%latticevecs(nbx + ycounter)  = 3
                    ! print*, 'v3sites', geo%hexsites(1,nbx + ycounter), geo%hexsites(2,nbx + ycounter)
                    if( pattern == 'AB' .and. modulo(x1, 2) == 0 ) then 
                        geo%cntrA = geo%cntrA + 1
                        alatt(1, geo%cntrA) = geo%hexsites(1,nbx + ycounter)
                        alatt(2, geo%cntrA) = geo%hexsites(2,nbx + ycounter)
                        alatt(3, geo%cntrA) = 1 !QAH phase 
                        alatt(4, geo%cntrA) = nbx + ycounter !Bond number 
                        alatt(5, geo%cntrA) = 3 !Lattice vector 
                        geo%phases(nbx + ycounter) = 1
                        cntrAsites(geo%hexsites(1,nbx + ycounter)) = cntrAsites(geo%hexsites(1,nbx + ycounter)) + 1
                        geo%AsitesBonds(cntrAsites(geo%hexsites(1,nbx + ycounter)),geo%hexsites(1,nbx + ycounter)) = geo%cntrA
                        cntrAsites(geo%hexsites(2,nbx + ycounter)) = cntrAsites(geo%hexsites(2,nbx + ycounter)) + 1
                        geo%AsitesBonds(cntrAsites(geo%hexsites(2,nbx + ycounter)),geo%hexsites(2,nbx + ycounter)) = -geo%cntrA
                    else if( pattern == 'AB' .and. modulo(x1, 2) == 1 ) then 
                        geo%cntrB = geo%cntrB + 1
                        blatt(1,geo%cntrB) = geo%hexsites(1,nbx + ycounter)
                        blatt(2,geo%cntrB) = geo%hexsites(2,nbx + ycounter)                    
                        blatt(3,geo%cntrB) = - 1 !QAH phase 
                        blatt(4,geo%cntrB) = nbx + ycounter !Bond number
                        blatt(5,geo%cntrB) = 3 !Lattice vector 
                        geo%phases(nbx + ycounter) = - 1
                        cntrBsites(geo%hexsites(1,nbx + ycounter)) = cntrBsites(geo%hexsites(1,nbx + ycounter)) + 1
                        geo%Bsitesbonds(cntrBsites(geo%hexsites(1,nbx + ycounter)),geo%hexsites(1,nbx + ycounter)) =geo%cntrB
                       cntrBsites(geo%hexsites(2,nbx + ycounter)) = cntrBsites(geo%hexsites(2,nbx + ycounter)) + 1
                        geo%Bsitesbonds(cntrBsites(geo%hexsites(2,nbx + ycounter)),geo%hexsites(2,nbx + ycounter)) = -geo%cntrB
                    else if( pattern == 'BA' .and. modulo(x1, 2) == 0 ) then 
                       geo%cntrB =geo%cntrB + 1
                        blatt(1,geo%cntrB) = geo%hexsites(1,nbx + ycounter)
                        blatt(2,geo%cntrB) = geo%hexsites(2,nbx + ycounter)
                        blatt(3,geo%cntrB) = - 1 !QAH phase 
                        blatt(4,geo%cntrB) = nbx + ycounter !Bond number
                        blatt(5,geo%cntrB) = 3 !Lattice vector 
                        geo%phases(nbx + ycounter) = - 1
                        cntrBsites(geo%hexsites(1,nbx + ycounter)) = cntrBsites(geo%hexsites(1,nbx + ycounter)) + 1
                        geo%Bsitesbonds(cntrBsites(geo%hexsites(1,nbx + ycounter)),geo%hexsites(1,nbx + ycounter)) =geo%cntrB
                        cntrBsites(geo%hexsites(2,nbx + ycounter)) = cntrBsites(geo%hexsites(2,nbx + ycounter)) + 1
                        geo%Bsitesbonds(cntrBsites(geo%hexsites(2,nbx + ycounter)),geo%hexsites(2,nbx + ycounter)) = -geo%cntrB
                    else if( pattern == 'BA' .and. modulo(x1, 2) == 1 ) then 
                        geo%cntrA = geo%cntrA + 1              
                        alatt(1, geo%cntrA) = geo%hexsites(1,nbx + ycounter)
                        alatt(2, geo%cntrA) = geo%hexsites(2,nbx + ycounter)  
                        alatt(3, geo%cntrA) = 1 !QAH phase       
                        alatt(4, geo%cntrA) = nbx + ycounter !Bond number     
                        alatt(5, geo%cntrA) = 3 !Lattice vector        
                        geo%phases(nbx + ycounter) = 1
                        cntrAsites(geo%hexsites(1,nbx + ycounter)) = cntrAsites(geo%hexsites(1,nbx + ycounter)) + 1
                        geo%BsitesBonds(cntrAsites(geo%hexsites(1,nbx + ycounter)),geo%hexsites(1,nbx + ycounter)) = geo%cntrA
                        cntrAsites(geo%hexsites(2,nbx + ycounter)) = cntrAsites(geo%hexsites(2,nbx + ycounter)) + 1
                        geo%BsitesBonds(cntrAsites(geo%hexsites(2,nbx + ycounter)),geo%hexsites(2,nbx + ycounter)) = -geo%cntrA
                    end if 
                    geo%xy(1, nbx + ycounter) = x1
                    geo%xy(2, nbx + ycounter) = y1
                    geo%xy(3, nbx + ycounter) = x2
                    geo%xy(4, nbx + ycounter) = y2
                    
                    ! geo%phases(nbx + ycounter) = - 1
                    60 continue
                    !if(y0 == 0 .and. x1 == nx - 1 .and. par%bc == 'o') goto 65

                    !Second vertical bond: Up left (lattice vector v2)
                    x2 = x1 !x-coordinate of second site
                    y2 = modulo(y1 + 1, ny) !y-coordinate of second site
                    if (y2 == y1) goto 65 !debug For 1D case no y-hopping
                    ycounter = ycounter + 1
                    geo%hexsites(1,nbx + ycounter) = 1 + x1 + y1 * nx
                    geo%hexsites(2,nbx + ycounter) = 1 + x2 + y2 * nx
                    geo%latticevecs(nbx + ycounter)  = 2
                    ytcounter = ytcounter + 1 
                    geo%ytransl(1, ytcounter) = geo%hexsites(1,nbx + ycounter)
                    geo%ytransl(2, ytcounter) = geo%hexsites(2,nbx + ycounter)
                                        
                    if(pattern == 'AB' .and. modulo(x1, 2) == 0 ) then 
                        geo%cntrA = geo%cntrA + 1
                        alatt(1, geo%cntrA) = geo%hexsites(1,nbx + ycounter)
                        alatt(2, geo%cntrA) = geo%hexsites(2,nbx + ycounter)   
                        alatt(3, geo%cntrA) = 1 !QAH phase     
                        alatt(4, geo%cntrA) = nbx + ycounter !Bond number       
                        alatt(5, geo%cntrA) = 2 !Lattice vector        
                        geo%phases(nbx + ycounter) = 1      
                        cntrAsites(geo%hexsites(1,nbx + ycounter)) = cntrAsites(geo%hexsites(1,nbx + ycounter)) + 1
                        geo%BsitesBonds(cntrAsites(geo%hexsites(1,nbx + ycounter)),geo%hexsites(1,nbx + ycounter)) = geo%cntrA
                        cntrAsites(geo%hexsites(2,nbx + ycounter)) = cntrAsites(geo%hexsites(2,nbx + ycounter)) + 1
                        geo%BsitesBonds(cntrAsites(geo%hexsites(2,nbx + ycounter)),geo%hexsites(2,nbx + ycounter)) = -geo%cntrA
                    else if( pattern == 'AB' .and. modulo(x1, 2) == 1 ) then 
                        geo%cntrB =geo%cntrB + 1
                        blatt(1,geo%cntrB) = geo%hexsites(1,nbx + ycounter)
                        blatt(2,geo%cntrB) = geo%hexsites(2,nbx + ycounter)  
                        blatt(3,geo%cntrB) = - 1 !QAH phase     
                        blatt(4,geo%cntrB) = nbx + ycounter !Bond number    
                        blatt(5,geo%cntrB) = 2 !Lattice vector            
                        geo%phases(nbx + ycounter) = - 1      
                        cntrBsites(geo%hexsites(1,nbx + ycounter)) =cntrBsites(geo%hexsites(1,nbx + ycounter)) + 1
                        geo%BsitesBonds(cntrBsites(geo%hexsites(1,nbx + ycounter)),geo%hexsites(1,nbx + ycounter)) =geo%cntrB
                        cntrBsites(geo%hexsites(2,nbx + ycounter)) =cntrBsites(geo%hexsites(2,nbx + ycounter)) + 1
                        geo%BsitesBonds(cntrBsites(geo%hexsites(2,nbx + ycounter)),geo%hexsites(2,nbx + ycounter)) = -geo%cntrB
                    else if( pattern == 'BA' .and. modulo(x1, 2) == 0 ) then 
                        geo%cntrB =geo%cntrB + 1
                        blatt(1,geo%cntrB) = geo%hexsites(1,nbx + ycounter)
                        blatt(2,geo%cntrB) = geo%hexsites(2,nbx + ycounter)
                        blatt(3,geo%cntrB) = 1 !QAH phase 
                        blatt(4,geo%cntrB) = nbx + ycounter !Bond number
                        blatt(5,geo%cntrB) = 2 !Lattice vector            
                        geo%phases(nbx + ycounter) = - 1
                        cntrBsites(geo%hexsites(1,nbx + ycounter)) =cntrBsites(geo%hexsites(1,nbx + ycounter)) + 1
                        geo%BsitesBonds(cntrBsites(geo%hexsites(1,nbx + ycounter)),geo%hexsites(1,nbx + ycounter)) =geo%cntrB
                        cntrBsites(geo%hexsites(2,nbx + ycounter)) =cntrBsites(geo%hexsites(2,nbx + ycounter)) + 1
                        geo%BsitesBonds(cntrBsites(geo%hexsites(2,nbx + ycounter)),geo%hexsites(2,nbx + ycounter)) = -geo%cntrB
                    else if( pattern == 'BA' .and. modulo(x1, 2) == 1 ) then 
                        geo%cntrA = geo%cntrA + 1
                        alatt(1, geo%cntrA) = geo%hexsites(1,nbx + ycounter)
                        alatt(2, geo%cntrA) = geo%hexsites(2,nbx + ycounter)      
                        alatt(3, geo%cntrA) = 1 !QAH phase        
                        alatt(4, geo%cntrA) = nbx + ycounter !Bond number     
                        alatt(5, geo%cntrA) = 2 !Lattice vector        
                        geo%phases(nbx + ycounter) = 1     
                        cntrAsites(geo%hexsites(1,nbx + ycounter)) = cntrAsites(geo%hexsites(1,nbx + ycounter)) + 1
                        geo%BsitesBonds(cntrAsites(geo%hexsites(1,nbx + ycounter)),geo%hexsites(1,nbx + ycounter)) = geo%cntrA
                        cntrAsites(geo%hexsites(2,nbx + ycounter)) = cntrAsites(geo%hexsites(2,nbx + ycounter)) + 1
                        geo%BsitesBonds(cntrAsites(geo%hexsites(2,nbx + ycounter)),geo%hexsites(2,nbx + ycounter)) = -geo%cntrA
                    end if                 
                    geo%xy(1, nbx + ycounter) = x1
                    geo%xy(2, nbx + ycounter) = y1
                    geo%xy(3, nbx + ycounter) = x2
                    geo%xy(4, nbx + ycounter) = y2
                    ! geo%phases(nbx + ycounter) = 1
                end if
                65 continue

                !Horizontal bond (lattice vector v1)
                if (x1 == nx-1 .and. par%bc == 'o') cycle !No NNN xbond to the right for last site
                if (x1 == nx-2 .and. par%bc == 'o') cycle !No NNN xbond to the right for second-to-last site
                xcounter = xcounter + 1
                x2 = modulo(x1 + 2, 2 * par%ucx)
                y2 = y1
                geo%hexsites(1,xcounter) = 1 + x1 + y1 * 2 * par%ucx
                geo%hexsites(2,xcounter) = 1 + x2 + y2 * 2 * par%ucx
                geo%latticevecs(xcounter)  = 1
                geo%xtransl(1, xcounter) = geo%hexsites(1,xcounter)
                geo%xtransl(2, xcounter) = geo%hexsites(2,xcounter)
                if( pattern == 'AB' .and. modulo(x1, 2) == 0 ) then 
                    geo%cntrA = geo%cntrA + 1
                    alatt(1, geo%cntrA) = geo%hexsites(1,xcounter)
                    alatt(2, geo%cntrA) = geo%hexsites(2,xcounter)      
                    alatt(3, geo%cntrA) = 1 !QAH phase               
                    alatt(4, geo%cntrA) = xcounter !Bond number   
                    alatt(5, geo%cntrA) = 1 !Lattice vector                
                    geo%phases(xcounter) = 1
                    cntrAsites(geo%hexsites(1,xcounter)) = cntrAsites(geo%hexsites(1,xcounter)) + 1
                    geo%BsitesBonds(cntrAsites(geo%hexsites(1,xcounter)),geo%hexsites(1,xcounter)) = geo%cntrA
                    cntrAsites(geo%hexsites(2,xcounter)) = cntrAsites(geo%hexsites(2,xcounter)) + 1
                    geo%BsitesBonds(cntrAsites(geo%hexsites(2,xcounter)),geo%hexsites(2,xcounter)) = -geo%cntrA
                else if( pattern == 'AB' .and. modulo(x1, 2) == 1 ) then 
                    geo%cntrB =geo%cntrB + 1
                    blatt(1,geo%cntrB) = geo%hexsites(1,xcounter)
                    blatt(2,geo%cntrB) = geo%hexsites(2,xcounter)
                    blatt(3,geo%cntrB) = - 1 !QAH phase  
                    blatt(4,geo%cntrB) = xcounter !Bond number    
                    blatt(5,geo%cntrB) = 1 !Lattice vector            
                    geo%phases(xcounter) = - 1              
                    cntrBsites(geo%hexsites(1,xcounter)) =cntrBsites(geo%hexsites(1,xcounter)) + 1
                    geo%BsitesBonds(cntrBsites(geo%hexsites(1,xcounter)),geo%hexsites(1,xcounter)) =geo%cntrB
                    cntrBsites(geo%hexsites(2,xcounter)) =cntrBsites(geo%hexsites(2,xcounter)) + 1
                    geo%BsitesBonds(cntrBsites(geo%hexsites(2,xcounter)),geo%hexsites(2,xcounter)) = -geo%cntrB
                else if( pattern == 'BA' .and. modulo(x1, 2) == 0 ) then 
                    geo%cntrB =geo%cntrB + 1
                    blatt(1,geo%cntrB) = geo%hexsites(1,xcounter)
                    blatt(2,geo%cntrB) = geo%hexsites(2,xcounter)
                    blatt(3,geo%cntrB) = - 1 !QAH phase
                    blatt(4,geo%cntrB) = xcounter !Bond number       
                    blatt(5,geo%cntrB) = 1 !Lattice vector             
                    geo%phases(xcounter) = - 1 
                    cntrBsites(geo%hexsites(1,xcounter)) =cntrBsites(geo%hexsites(1,xcounter)) + 1
                    geo%BsitesBonds(cntrBsites(geo%hexsites(1,xcounter)),geo%hexsites(1,xcounter)) =geo%cntrB
                    cntrBsites(geo%hexsites(2,xcounter)) =cntrBsites(geo%hexsites(2,xcounter)) + 1
                    geo%BsitesBonds(cntrBsites(geo%hexsites(2,xcounter)),geo%hexsites(2,xcounter)) = -geo%cntrB
                else if( pattern == 'BA' .and. modulo(x1, 2) == 1 ) then 
                    geo%cntrA = geo%cntrA + 1
                    alatt(1, geo%cntrA) = geo%hexsites(1,xcounter)
                    alatt(2, geo%cntrA) = geo%hexsites(2,xcounter)           
                    alatt(3, geo%cntrA) = 1 !QAH phase           
                    alatt(4, geo%cntrA) = xcounter !Bond number    
                    alatt(5, geo%cntrA) = 1 !Lattice vector                
                    geo%phases(xcounter) = 1
                    cntrAsites(geo%hexsites(1,xcounter)) = cntrAsites(geo%hexsites(1,xcounter)) + 1
                    geo%BsitesBonds(cntrAsites(geo%hexsites(1,xcounter)),geo%hexsites(1,xcounter)) = geo%cntrA
                    cntrAsites(geo%hexsites(2,xcounter)) = cntrAsites(geo%hexsites(2,xcounter)) + 1
                    geo%BsitesBonds(cntrAsites(geo%hexsites(2,xcounter)),geo%hexsites(2,xcounter)) = -geo%cntrA
                end if    
                geo%xy(1, xcounter) = x1
                geo%xy(2, xcounter) = y1
                geo%xy(3, xcounter) = x2
                geo%xy(4, xcounter) = y2
                geo%sitecoord(1, geo%hexsites(1,xcounter)) = x1 
                geo%sitecoord(2, geo%hexsites(1,xcounter)) = y1 
                ! geo%phases(xcounter) = 1
            end do
        end do

        allocate(geo%Alattice(5,geo%cntrA))
        allocate(geo%Blattice(5,geo%cntrB))
        geo%Alattice(1, 1:geo%cntrA) = alatt(1, 1:geo%cntrA)
        geo%Alattice(2, 1:geo%cntrA) = alatt(2, 1:geo%cntrA)
        geo%Alattice(3, 1:geo%cntrA) = alatt(3, 1:geo%cntrA)
        geo%Alattice(4, 1:geo%cntrA) = alatt(4, 1:geo%cntrA)
        geo%Alattice(5, 1:geo%cntrA) = alatt(5, 1:geo%cntrA)
        geo%Blattice(1, 1:geo%cntrB) = blatt(1, 1:geo%cntrB)
        geo%Blattice(2, 1:geo%cntrB) = blatt(2, 1:geo%cntrB)
        geo%Blattice(3, 1:geo%cntrB) = blatt(3, 1:geo%cntrB)
        geo%Blattice(4, 1:geo%cntrB) = blatt(4, 1:geo%cntrB)
        geo%Blattice(5, 1:geo%cntrB) = blatt(5, 1:geo%cntrB)
        deallocate(alatt)
        deallocate(blatt)
        write(name,"('L=',i0,'ucx=',i0,'ucy=',i0,'pat=',a2,'par%BC=',a,'.dat')") nn, par%ucx, par%ucy, pattern, par%bc
        name = trim(dir)//"A_lattice_"//trim(name)
        open(unit=32, file=name)
        do s = 1, geo%cntrA
            write(32,'(i0,x,i0,x,i0,x,i0,x,i0)') geo%alattice(1,s), geo%alattice(2,s), geo%alattice(3,s), geo%alattice(4,s), geo%alattice(5,s)
        end do
        close(32)


        write(name,"('L=',i0,'ucx=',i0,'ucy=',i0,'pat=',a2,'par%BC=',a,'.dat')") nn, par%ucx, par%ucy, pattern, par%bc
        name = trim(name)
        name=dir//"B_lattice_"//name
        open(unit=32, file=name)
        do s = 1,geo%cntrB
            write(32,'(i0,x,i0,x,i0,x,i0,x,i0)') geo%blattice(1,s), geo%blattice(2,s), geo%blattice(3,s), geo%blattice(4,s), geo%blattice(5,s)
        end do
        close(32)

        geo%nnBonds = xcounter + ycounter

        ! allocate(geo%phases(nnBonds))
        ! geo%phases(1:xcounter) = geo%phases(1:xcounter)
        ! geo%phases(xcounter + 1:xcounter +  ycounter) = geo%phases(nbx + 1:nbx +  ycounter)

        write(name,"('L=',i0,'ucx=',i0,'ucy=',i0,'pat=',a2,'par%BC=',a,'.dat')") nn, par%ucx, par%ucy, pattern, par%bc
        name = trim(name)
        name = dir//"NNN_lattice_"//name
        open(unit=32, file=name)
        do s = 1, geo%nnBonds
            write(32,*) geo%hexsites(1,s), geo%hexsites(2,s)      
        end do

        close(32)

        write(name,"('L=',i0,'ucx=',i0,'ucy=',i0,'pat=',a2,'par%BC=',a,'.dat')") nn, par%ucx, par%ucy,pattern, par%bc
        name = trim(name)
        name=dir//"NNN_xy_"//name
        open(unit=32, file=name)
        do s = 1, geo%nnBonds
            write(32,'(4(i0,x))') geo%xy(1,s), geo%xy(2,s), geo%xy(3,s), geo%xy(4,s)
        end do
        close(32)

        ! print*, 'Finished honeycomb nnn-lattice'

    end subroutine honeycomb_nnn




end module lattice 