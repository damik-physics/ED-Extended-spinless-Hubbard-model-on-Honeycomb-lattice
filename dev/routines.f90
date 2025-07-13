module routines

    implicit none

    !Spartan changes
    ! replace all dir//"..." in file names by "..dir///...". 
    ! for Hamitlonian I/O saving replace "hamiltonians/..." in file names by "../hamiltonians/...". 



    include "omp_lib.h"
    
    !INTERFACES
    interface save
        module procedure savespec, savecurr, savecdw, savemom
    end interface  

    interface check
        module procedure check_d, check_c
    end interface  

    interface gsdeg
        module procedure gsdeg_d, gsdeg_c
    end interface 

    interface cdw
        module procedure cdw_d, cdw_c
    end interface   

    interface unify 
        module procedure unify, unify_dp, unify_comp, unify_dense
    end interface    

    double precision, parameter :: pi = 4*atan(1.d0)
    double complex, parameter   :: ii = (0, 1)


contains

!---------------------------!
!         FUNCTIONS         !
!---------------------------!

function trim_name(file_name) !Function to trim file names

    implicit none

    character*512 :: trim_name
    character(len=*):: file_name
    integer :: ls1, ls2, i

    trim_name=''
    ls1 = len_trim(file_name)
    ls2 = 0
    do i = 1, ls1
        if(file_name(i:i).ne.' ') then
           ls2 = ls2 + 1
           trim_name(ls2:ls2) = file_name(i:i)
        end if
    end do
    return
end function

integer function combs (nsites, part, magn) !Calculates the number of possible combinations for a given magnetization and particle number!
    !-------------------------------------------!
    !            Number of basis states         !
    !-------------------------------------------!
    implicit none

    integer, intent(in) :: nsites, part
    double precision, intent(in) :: magn
    integer :: part_up, part_dn

    part_up=int((part+2*magn)/2)
    part_dn=part-part_up

    if (part_up < 0 .or. part_dn < 0) then
        combs = 0
    else
        combs = int(fact(nsites)/(fact(part_up)*fact(max(1,nsites-part_up))),8)*int(fact(nsites)/(fact(part_dn)*fact(max(1,nsites-part_dn))),8)
    end if

end function combs

double precision function fact(nn) !Calculates the factorial n!

    !----------------------------------!
    !            Factorial: n!         !
    !----------------------------------!

    implicit none

    integer, intent(in) :: nn
    integer :: i

    if (nn < 0) error stop 'factorial is singular for negative integers'

    if(nn==0) then
        fact=1
    else
        fact = 1
        do i = 2, nn
            fact = fact * i
        end do
    end if

end function fact

double precision function binomial(n, k)

    !---------------------------------------!
    !            Binomial coefficients      !
    !---------------------------------------!

    implicit none
    integer, intent(in) :: n, k

    binomial = fact(n) / (fact(k) * fact(n - k))

end function binomial

function getPosition(b1, b2, N, M)

    implicit none
    integer, intent(in) :: b1, b2, N, M
    integer :: getPosition
    integer :: i, offset_b1, offset_b2
  
    offset_b1 = 0
    do i = 0, popcnt(b1)-1!b1-1
       offset_b1 = offset_b1 + binomial(N/2, i)
    end do
  
    offset_b2 = 0
    do i = 0, popcnt(b2)-1! b2-1
       offset_b2 = offset_b2 + binomial(N/2, i)
    end do
  
    getPosition = offset_b1 + offset_b2 + 1
end function getPosition

double precision function log2(x)
  
    implicit none

    double precision, intent(in) :: x   

    log2 = log(x) / log(2.)

end function

!-----------------------------!
!         SUBROUTINES         !
!-----------------------------!

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

!-----------------------------!
!         Save data           !
!-----------------------------!

subroutine savemom(dir, cluster, unit, gencluster, sites, particles, bc, pattern, k1_max, k2_max)
    
    implicit none 
    
    character(len=*), intent(in) :: dir, cluster, bc, pattern
    integer, intent(in) :: unit, gencluster, sites, particles, k1_max, k2_max

    integer :: k1 = 0, k2 = 0 
    character :: file*256, prefix*20, params*256

    if(gencluster == 0) then 
        write(prefix, "(a,'_')"), cluster 
        prefix = trim_name(prefix)
    end if 
    
    write (params,"('L=',i0,'N=',i0,'BC=',a,'_pat=',a2'.dat')") sites,particles,bc,pattern
    params=trim_name(prefix//params)
    file=trim_name(dir//"momenta_"//params)
    
    open(unit,file=file)
    do k1 = 0, k1_max
        do k2 = 0, k2_max
            write(unit,*) k1, k2 
        end do 
    end do
    close(unit)

    return 

end subroutine savemom

subroutine savespec(dir, type, parameters, append, conf, unit, routine, dim, states, nev, nest, nDis, rvec, en, st_dp, st_c)
    
    implicit none 
    
    character(len=*), intent(in) :: dir, type, parameters
    integer, intent(in) :: conf, unit, routine 
    integer(kind=8), intent(in) :: dim  
    integer, intent(in) :: states, nev, nest, nDis
    logical, intent(in) :: append, rvec 
    double precision, intent(in), optional :: en(dim), st_dp(dim, nest) 
    double complex, intent(in), optional :: st_c(dim, nest) 

    integer :: j = 0
    character :: file*256, file2*256, appchar*20

    100 format(1000(F40.30))
    
    if(append) then 
        appchar = 'APPEND'
    else  
        appchar = 'SEQUENTIAL'
    end if 

    file=dir//"energies_"//parameters
    file2=dir//"states_"//parameters
    
    file=trim_name(file)
    file2=trim_name(file2)
    open(unit,file=file, access=appchar)
    if(conf == 1) write(unit,*) nev

    do j = 1, nev
        write(unit,100) en(j)
    end do
    close(unit)

    if(states == 1 .and. rvec) then
        if(type == "R") then 
            open(unit,file=file2, access=appchar)
            write(unit,*) dim
            write(unit,*) nev
            do j = 1, nest
                write(unit,100) st_dp(1:dim, j)
            end do
            close(unit)
        else if(type == "C") then 
            open(unit,file=file2, access=appchar)
            write(unit,*) dim
            write(unit,*) nev
            do j = 1, nest               
                write(unit,100) dble(st_c(1:dim, j)), dimag(st_c(1:dim, j)) 
            end do
            close(unit)
        end if
    end if 

    if(routine == 0) then 
        file=dir//"exact_energies_"//parameters
        file2=dir//"exact_states_"//parameters
    else if(routine == 1) then 
        file=dir//"feast_energies_"//parameters
        file2=dir//"feast_states_"//parameters
    else if(routine == 2) then 
        file=dir//"mkl_energies_"//parameters
        file2=dir//"mkl_states_"//parameters
    else if(routine == 3) then 
        file=dir//"arpack_energies_"//parameters
        file2=dir//"arpack_states_"//parameters
    end if 

    ! file=trim_name(file)
    ! file2=trim_name(file2)
    ! open(unit,file=file, access=appchar)
    ! if(conf == 1) write(unit,*) nev

    ! do j = 1, nev
    !     write(unit,100) en(j)
    ! end do
    ! close(unit)

    ! if(states == 1 .and. rvec) then
    !     if(type == "R") then 
    !         open(unit,file=file2, access=appchar)
    !         write(unit,*) dim
    !         write(unit,*) nev
    !         do j = 1, nest
    !             write(unit,100) st_dp(1:dim, j)
    !         end do
    !         close(unit)
    !     else if(type == "C") then 
    !         open(unit,file=file2, access=appchar)
    !         write(unit,*) dim
    !         write(unit,*) nev
    !         do j = 1, nest               
    !             write(unit,100) dble(st_c(1:dim, j)), dimag(st_c(1:dim, j)) 
    !         end do
    !         close(unit)
    !     end if
    ! end if 


    return 

end subroutine savespec

subroutine saveham(spartan, type, parameters, unit, ti, sites, nnz, nDi, hamDi, ham_i, rc, ham_dp, ham_dc, occ)
    implicit none 
    character(len=*), intent(in) :: type, parameters
    integer, intent(in) :: spartan, ti, unit, sites  
    integer, intent(in) :: nnz, nDi       
    integer, intent(in), optional :: hamDi(nDi, 2), ham_i(nnz, 3), rc(nnz, 2), occ(sites, nDi)
    double precision, intent(in), optional :: ham_dp(nnz)
    double complex, intent(in), optional :: ham_dc(nnz) 

    integer :: i = 0
    character :: file*256, dir*100

    print*, 'Saving Hamiltonian to disk...'
    if(spartan == 0) then 
        dir = "hamiltonians/"
    else if(spartan == 1) then 
        dir = "../hamiltonians/"
    end if 
    dir = trim_name(dir)
    file=dir//"hopping_hamiltonian_"//parameters
    file=trim_name(file)
    open(unit, file=file)
    write(unit,*) nnz 
    
    if(type == "R" .and. ti == 0) then 
        do i = 1, nnz
            if(ti == 0) write(unit,*) ham_i(i,1), ham_i(i,2), ham_i(i,3)
        end do 
    else if(type == "R" .and. ti == 1) then         
        do i = 1, nnz
            write(unit,*) ham_dp(i), rc(i,1), rc(i,2)
        end do 
    else if(type == "C" .and. ti == 1) then             
        do i = 1, nnz
            write(unit,*) ham_dc(i), rc(i,1), rc(i,2)
        end do 
    end if 
    close(unit)
    
    file = ""
    file = trim_name(dir//"diagonal_hamiltonian_"//parameters)
    open(unit, file=file)
    write(unit,*) nDi 
    do i = 1, nDi
        write(unit,*) hamDi(i,1), hamDi(i,2)
    end do 
    close(unit)

    file= trim_name(dir//"occupation_"//parameters)
    open(unit, file=file)
    write(unit,*) nDi 
    do i = 1, nDi
        write(unit,*) occ(1:sites,i)
    end do 
    close(unit)
    print*,'Hamiltonian saved to disk.'
    print*,''
    return 

end subroutine saveham

subroutine loadham(spartan, type, parameters, unit, ti, sites, nnz, nDi, hamDi, ham_i, rc, ham_dp, ham_dc, occupation, exist1, exist2, exist3)
    
    implicit none 
    
    character(len=*), intent(in) :: type, parameters
    integer, intent(in) :: spartan, ti, unit, sites  
    integer, intent(out) :: nnz, nDi       
    integer, allocatable, intent(out), optional :: hamDi(:, :), ham_i(:, :), rc(:, :), occupation(:, :)
    double precision, allocatable, intent(out), optional :: ham_dp(:)
    double complex, allocatable, intent(out), optional :: ham_dc(:) 
    logical, intent(out) :: exist1, exist2, exist3

    integer :: i = 0
    character :: file*256, dir*100
    
    if(spartan == 0) then 
        dir = "hamiltonians/"
    else if(spartan == 1) then 
        dir = "../hamiltonians/"
    end if 
    dir = trim_name(dir)

    file=dir//"hopping_hamiltonian_"//parameters
    file=trim_name(file)
    inquire(file=file, exist=exist1)

    if( exist1 ) then 
        open(unit, file=file)
        read(unit, *) nnz 

        if(type == "R" .and. ti == 0) then 
            if(allocated(ham_i)) deallocate(ham_i)
            allocate(ham_i(nnz, 3))
            ham_i = 0
            do i = 1, nnz
                read(unit,*) ham_i(i,1), ham_i(i,2), ham_i(i,3)
            end do 
        else if(type == "R" .and. ti == 1) then 
            if(allocated(ham_dp)) deallocate(ham_dp)
            if(allocated(rc)) deallocate(rc)
            allocate(ham_dp(nnz))
            allocate(rc(nnz, 2))
            ham_d = 0.d0
            rc = 0

            do i = 1, nnz
                read(unit,*) ham_dp(i), rc(i,1), rc(i,2)
            end do 
        else if(type == "C" .and. ti == 1) then      
            if(allocated(ham_dc)) deallocate(ham_dc)
            if(allocated(rc)) deallocate(rc)
            allocate(ham_dc(nnz))
            allocate(rc(nnz, 2))       
            do i = 1, nnz
                read(unit,*) ham_dc(i), rc(i,1), rc(i,2)
            end do 
        end if 
        close(unit)
        print*,'Off-diagonal Hamiltonian uploaded.'
        print*,''
    end if

    if(.not. (exist1)) print*,'Off-diagonal Hamiltonian does not yet exist.'

    file = ""
    file = dir//"diagonal_hamiltonian_"//parameters
    file = trim_name(file)
    inquire(file=file, exist=exist2)

    if( exist2 ) then 
        open(unit, file=file)
        read(unit, *) nDi 
        if(allocated(hamDi)) deallocate(hamDi)
        allocate(hamDi(nDi, 2))
        do i = 1, nDi
            read(unit,*) hamDi(i,1), hamDi(i,2)
        end do 
        close(unit)
        print*,'Diagonal Hamiltonian uploaded.'
        print*,''
    end if 
    
    if(.not. (exist2)) print*,'Diagonal Hamiltonian does not yet exist.'

    file = ""
    file = "hamiltonians/occupation_"//parameters
    file = trim_name(file)
    inquire(file=file, exist=exist3)
    if( exist3 ) then 
        open(unit, file=file)
        read(unit, *) nDi 
        if(allocated(occupation)) deallocate(occupation)
        allocate(occupation(sites, nDi))
        do i = 1, nDi
            read(unit,*) occupation(1:sites,i)
        end do 
        close(unit)
        print*,'Occupation uploaded.'
        print*,''
    end if 

    if(.not. (exist3)) print*,'Occupation does not yet exist.'
    if(.not. (exist3)) print*,''

    return 

end subroutine loadham

subroutine savecurr(dir, sl, parameters, append, degeneracy, unit, full, feast, mkl, arpack, nbonds, current, bondcurrent)
    
    implicit none 
    character(len=*), intent(in) :: dir, sl, parameters
    integer, intent(in) :: degeneracy, unit, nbonds
    integer, intent(in) :: full, feast, mkl, arpack  
    double precision, intent(in) :: current, bondcurrent(nbonds)
    logical, intent(in) :: append 

    character :: file*512, file2*512, file3*512, file4*512, appchar*20
    character :: dirq*512
    100 format(1000(F30.20))

    if(append) then 
        appchar = 'APPEND'
    else  
        appchar = 'SEQUENTIAL'
    end if 

    if(degeneracy < 2 ) dirq = dir  
    if(degeneracy == 2 ) dirq = trim_name(dir//"QD_")
    dirq = trim_name(dirq)
    
    if(sl == "A") then 
        file=dirq//"Acurrent_"//parameters
        file=trim_name(file)
        file2=dirq//"Abondcurrent_"//parameters
        file2=trim_name(file2)
        if(full == 1) then 
            file3=dirq//"exact_Acurrent_"//parameters
            file3=trim_name(file3)
            file4=dirq//"exact_Abondcurrent_"//parameters
            file4=trim_name(file4)
        else if(feast == 1) then 
            file3=dirq//"feast_Acurrent_"//parameters
            file3=trim_name(file3)
            file4=dirq//"feast_Abondcurrent_"//parameters
            file4=trim_name(file4)
        else if(mkl == 1) then 
            file3=dirq//"mkl_Acurrent_"//parameters
            file3=trim_name(file3)
            file4=dirq//"mkl_Abondcurrent_"//parameters
            file4=trim_name(file4)
        else if(arpack == 1) then 
            file3=dirq//"arpack_Acurrent_"//parameters
            file3=trim_name(file3)
            file4=dirq//"arpack_Abondcurrent_"//parameters
            file4=trim_name(file4)
        end if 
    else if(sl == "B") then 
        file=dirq//"Bcurrent_"//parameters
        file=trim_name(file)
        file2=dirq//"Bbondcurrent_"//parameters
        file2=trim_name(file2)
        if(full == 1) then 
            file3=dirq//"exact_Bcurrent_"//parameters
            file3=trim_name(file3)
            file4=dirq//"exact_Bbondcurrent_"//parameters
            file4=trim_name(file4)
        else if(feast == 1) then 
            file3=dirq//"feast_Bcurrent_"//parameters
            file3=trim_name(file3)
            file4=dirq//"feast_Bbondcurrent_"//parameters
            file4=trim_name(file4)
        else if(mkl == 1) then 
            file3=dirq//"mkl_Bcurrent_"//parameters
            file3=trim_name(file3)
            file4=dirq//"mkl_Bbondcurrent_"//parameters
            file4=trim_name(file4)
        else if(arpack == 1) then 
            file3=dirq//"arpack_Bcurrent_"//parameters
            file3=trim_name(file3)
            file4=dirq//"arpack_Bbondcurrent_"//parameters
            file4=trim_name(file4)
        end if 
    end if 

    open(unit,file=file, access=appchar)
    write(unit,100) current
    close(unit)
    open(unit,file=file2, access=appchar)
    write(unit,100) bondcurrent
    close(unit)
    ! open(unit,file=file3, access=appchar)
    ! write(unit,100) current
    ! close(unit)
    ! open(unit,file=file4, access=appchar)
    ! write(unit,100) bondcurrent
    ! close(unit)

    return 

end subroutine savecurr

subroutine savecdw(dir, parameters, append, unit, ucx, ucy, nbonds, rho, rhotot)
    implicit none 
    character(len=*), intent(in) :: dir, parameters
    integer, intent(in) :: unit 
    integer, intent(in) :: ucx, ucy     
    integer, intent(in) :: nbonds  
    double precision, intent(in) :: rho(ucx*ucy, 3) 
    double precision, intent(in) :: rhotot(nbonds, 3) 
    logical, intent(in) :: append 

    integer :: i = 0
    character :: file*256, file2*256, appchar*20

    100 format(1000(F30.20))
     
    if(append) then 
        appchar = 'APPEND'
    else  
        appchar = 'SEQUENTIAL'
    end if 

    file=dir//"cdw_"//parameters
    file=trim_name(file)
    file2=dir//"cdw_tot_"//parameters
    file2=trim_name(file2)
    
    open(unit,file=file, access=appchar)
    do i = 1, ucx*ucy
        write(unit,100) rho(i,1), rho(i,2), rho(i,3)    
    end do 
    close(unit)
    open(unit,file=file2, access=appchar)
    do i = 1, nbonds
        write(unit,100) rhotot(i,1), rhotot(i,2), rhotot(i,3)    
    end do 
    close(unit)

    return 

end subroutine savecdw

subroutine saveddcf(dir, parameters, append, unit, sites, rhocn, rho)
    implicit none 
    character(len=*), intent(in) :: dir, parameters
    integer, intent(in) :: unit, sites
    double precision, intent(in) :: rhocn(sites), rho(sites) 
    logical, intent(in) :: append 

    integer :: i = 0
    character :: file*400, file2*400, appchar*20

    100 format(1000(F30.20))
     
    if(append) then 
        appchar = 'APPEND'
    else  
        appchar = 'SEQUENTIAL'
    end if 

    file=dir//"ddcf_"//parameters
    file=trim_name(file)
    file2=dir//"density_"//parameters
    file2=trim_name(file2)

    open(unit,file=file, access=appchar)
    do i = 1, sites
        write(unit,100) rhocn(i)
    end do 
    close(unit)
    open(unit,file=file2, access=appchar)
    do i = 1, sites
        write(unit,100) rho(i)
    end do 
    close(unit)

    return 

end subroutine saveddcf

subroutine saveI(sa, unit, dir, name, pars, rows, cols, scalar, arr)
    
    implicit none 
    
    character(len=*), intent(in) :: sa, pars, name, dir
    integer, intent(in) :: unit
    integer(kind=8), intent(in) :: rows, cols  
    integer, intent(in) :: scalar
    integer, intent(in), optional :: arr(rows, cols) 
    
    integer(kind=8) :: j = 0
    character :: file*512

    j = 0 
    ! 100 format(1000(F40.30))
    
    file=trim_name(dir//name//pars)
    print* ,file, 'file'    
    open(unit,file=file)

    if(sa == "S") then 
        write(unit,*) scalar
    end if 
    if(sa == "A") then  
        do j = 1, cols
            write(unit,*) arr(1:rows, j)
        end do
    end if 
    close(unit)

    return 

end subroutine saveI

subroutine saveI8(sa, unit, dir, name, pars, rows, cols, scalar, arr)
    
    implicit none 
    
    character(len=*), intent(in) :: sa, pars, name, dir
    integer, intent(in) :: unit
    integer(kind=8), intent(in) :: rows, cols  
    integer(kind=8), intent(in) :: scalar
    integer(kind=8), intent(in), optional :: arr(rows, cols) 
    
    integer(kind=8) :: j = 0
    character :: file*512
    
    j = 0 
    file=trim_name(dir//name//pars)

    
    
    
    open(unit,file=file)

    if(sa == "S") then 
        write(unit,*) scalar
    end if 
    if(sa == "A") then 
        do j = 1, cols
            write(unit,*) arr(1:rows, j)
        end do
    end if 
    close(unit)

    return 

end subroutine saveI8

subroutine saveD(sa, unit, dir, name, pars, rows, cols, scalar, arr, transp)
    
    implicit none 
    
    character(len=*), intent(in) :: sa, pars, name, dir
    integer, intent(in) :: unit
    integer, intent(in), optional :: transp
    integer(kind=8), intent(in) :: rows, cols  
    double precision, intent(in) :: scalar
    double precision, intent(in), optional :: arr(rows, cols) 
    
    integer(kind=8) :: j = 0
    character :: file*512

    j = 0 
    99 format(1000(F40.30))
    
    file=trim_name(dir//name//pars)
    open(unit,file=file)
     
    if(sa == "S") then 
        write(unit,99) scalar
    end if 
    if(sa == "A") then 
        if(rows < 1000) then 
            if(present(transp)) then 
                do j = 1, rows
                    write(unit,100) arr(j, 1:cols)
                end do
            else
                do j = 1, cols
                    write(unit,100) arr(1:rows, j)
                end do
            end if 
        else 
            100 format(1000000000(F40.30))
            do j = 1, cols
                write(unit,100) arr(1:rows, j)
            end do
        end if
    end if 
    close(unit)

    return 

end subroutine saveD

subroutine saveZ(sa, unit, dir, name, pars, rows, cols, scalar, arr, transp)
    
    implicit none 
    
    character(len=*), intent(in) :: sa, pars, name, dir
    integer, intent(in) :: unit
    integer(kind=8), intent(in) :: rows, cols  
    integer, intent(in), optional :: transp
    double complex, intent(in) :: scalar
    double complex, intent(in), optional :: arr(rows, cols) 
    
    integer(kind=8) :: j = 0
    character :: file*512

    j = 0 
    100 format(1000(F40.30))
    
    file=trim_name(dir//name//pars)
   
    open(unit,file=file)

    if(sa == "S") then 
        write(unit,100) scalar
    else if(sa == "A") then 
        if(present(transp)) then 
            do j = 1, rows
                write(unit,100) arr(j, 1:cols)
            end do
        else
            do j = 1, cols
                write(unit,100) arr(1:rows, j)
            end do
        end if 
    end if 
    close(unit)

    return 

end subroutine saveZ

subroutine par2(unit, k1, k2, refb, v1, v22, params, type)
    
    use variables

    implicit none 
    integer, intent(in) :: unit, k1, k2, refb
    double precision, intent(in) :: v1, v22 
    character, intent(out) :: params*512 
    character, intent(out), optional :: type*1 
    

    character :: string*100

    if(tilted == 1) then 
        write(string, "(a,'_')"), cluster 
        string = trim_name(string)
    end if 
    params = ""    
    if(refb < 0) then !Name for Hamiltonian file
        type = 'R'
        if(((k1 .ne. 0) .or. (k2 .ne. 0))) type = 'C'
        if(ti == 0) write (params,"('L=',i0,'N=',i0,'BC=',a,'_pat=',a2'.dat')") sites,particles,bc,pattern
        if(ti == 1) write (params,"('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'BC=',a,'_pat=',a2'.dat')") sites,particles,k1,k2,bc,pattern
        if(symmetrize == 1) write (params,"('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'_',a,'_BC=',a,'_pat=',a2'.dat')") sites,particles,k1,k2,irrep,bc,pattern
    else if(refb == 0) then !Name for energy/states file
        if(ti == 0 ) write (params,"('L=',i0,'N=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2'.dat')") sites,particles,v1,v22,mass,dis,nDis,bc,pattern
        if(ti == 1) write (params,"('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2'.dat')") sites,particles,k1,k2,v1,v22,mass,dis,nDis,bc,pattern
        if(symmetrize == 1) write (params,"('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'_',a,'_V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2'.dat')") sites,particles,k1,k2,irrep,v1,v22,mass,dis,nDis,bc,pattern
    else if(refb > 0) then !Name for current files
        if(ti == 0 ) write (params,"('L=',i0,'N=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refb=',i0,'.dat')") sites,particles,v1,v22,mass,dis,nDis,bc,pattern,refb            
        if(ti == 1) write (params,"('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refb=',i0,'.dat')") sites,particles,k1,k2,v1,v22,mass,dis,nDis,bc,pattern,refb 
        if(symmetrize == 1) write (params,"('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'_',a,'_V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refb=',i0,'.dat')") sites,particles,k1,k2,irrep,v1,v22,mass,dis,nDis,bc,pattern,refb
    end if 
    params=trim_name(params)
    if(tilted == 1) params=trim_name(string//params)
    params=trim_name(params)    

end subroutine par2

subroutine parddcf(unit, k1, k2, refsite, v1, v22, params)
    
    use variables

    implicit none 
    integer, intent(in) :: unit, k1, k2, refsite
    double precision, intent(in) :: v1, v22 
    character, intent(out) :: params*400 

    character :: string*100

    if(tilted == 1) then 
        write(string, "(a,'_')"), cluster 
        string = trim_name(string)
    end if 
    params = ""    
    if(refsite > 0) then !Name for current files
        if(ti == 0 ) write (params,"('L=',i0,'N=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refsite=',i0,'.dat')") sites,particles,v1,v22,mass,dis,nDis,bc,pattern,refsite            
        if(ti == 1) write (params,"('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refsite=',i0,'.dat')") sites,particles,k1,k2,v1,v22,mass,dis,nDis,bc,pattern,refsite 
        if(symmetrize == 1) write (params,"('L=',i0,'N=',i0,'k1=',i0,'k2=',i0,a,'V=',f7.4,'V2=',f7.4,'M=',f10.7,'W=',f7.4,'nDis=',i0,'BC=',a,'_pat=',a2,'_refsite=',i0,'.dat')") sites,particles,k1,k2,irrep,v1,v22,mass,dis,nDis,bc,pattern,refsite
    end if 
    params=trim_name(params)
    if(tilted == 1) then 
        params=trim_name(string//params)
    end if 
    params=trim_name(params)    

end subroutine parddcf

subroutine check_d(dim, nnz, nev, nest, energies, eigstate, ndeg, rc, mat)
    implicit none
    integer(kind=8), intent(in) :: dim, nnz  
    integer, intent(in) :: nev, nest
    integer, intent(in), optional :: ndeg
    integer, intent(in), optional :: rc(nnz, 2)
    double precision, intent(inout) :: energies(nev)
    double precision, intent(inout) :: eigstate(dim, nest)
    double precision, intent(in), optional :: mat(nnz)

    integer :: i = 0 
    
    ! call sortstates(nev, nest, dim, energies, eigstate)
    call checnorm(nev, nest, dim, energies, eigstate)
    call checkorthogonality(nest, dim, eigstate)
    ! if(present(mat) .and. present(rc)) then 
    !     ! print*,'Degeneracy', ndeg 
    !     do i = 1, max(1, ndeg)
    !         call expectval_d(dim, nnz, energies(i), eigstate(1:dim, i), rc, mat)
    !     end do 
    ! end if 
    return 
    
end subroutine check_d

subroutine check_c(dim, feast, nev, nest, energies, eigstate)
    implicit none
    integer(kind=8), intent(in) :: dim  
    integer, intent(in) :: nev, nest
    logical, intent(in) :: feast
    double precision, intent(inout) :: energies(nev)
    double complex, intent(inout) :: eigstate(dim, nest)

    call sortstates_comp(nev, nest, dim, energies, eigstate)
    call checnorm_comp(nev, nest, dim, energies, eigstate)
    if(feast) then 
        call checkorthogonality_comp(nest, dim, eigstate)   
    end if 

    return 

end subroutine check_c

subroutine expectval_d(dim, nnz, eval, evec, rc, mat)
    implicit none 
    integer(kind=8), intent(in) :: dim, nnz  
    integer, intent(in) :: rc(nnz, 2)
    double precision, intent(in) :: eval
    double precision, intent(in) :: evec(dim)
    double precision, intent(in) :: mat(nnz)
    
    integer :: ia(dim + 1), ja(nnz)
    double precision :: exv, val(nnz), ax(dim)

    ia = 0 
    ja = 0 
    exv = 0.d0 
    val = 0.d0 

    call coocsr(dim, nnz, mat, rc(1:nnz,1), rc(1:nnz,2), val, ja, ia) 
    call pamux (1, dim, evec, ax, val, ja, ia)
    exv = dot_product(evec, ax)
    ! print*,'Expectation value ',exv 
    ! print*,'exv-eval',exv-eval 

    return 
    
    ! call mkl_dcsrmultcsr('N', 0, 3, n, n, n, val, cols, rows, val_b, cols_b, rows_b, c, jc, ic, nnz, info)


end subroutine expectval_d

subroutine sortstates(nev, nest, dim, evals, states)
    implicit none 
    integer, intent(in) :: nev, nest 
    integer(kind=8), intent(in) :: dim 

    double precision, intent(inout) :: evals(nev), states(dim, nev) 

    integer, save :: i = 0, j = 0, cntr = 0 
    double precision, save :: temp = 0 
    double precision, allocatable, save :: tempvec(:)

    cntr = 0
    if(allocated(tempvec)) deallocate(tempvec)
    allocate(tempvec(dim))
    do i = 1, nest
        do j = 1, nest
            if (dble(evals(j)) < dble(evals(i)) .and. i < j) then 
                print*,'Swap eigenvalues i, j',i, j
                temp = evals(i)
                evals(i) = evals(j)
                evals(j) = temp 
                cntr = cntr + 1 
                if(dble(evals(i)) > dble(evals(j))) print*,'Error in resorting routine.'
                tempvec = states(i, 1:nest)
                states(i, 1:nest) = states(j, 1:nest)
                states(j, 1:nest) = tempvec 
            end if 
        end do 
    end do 
    deallocate(tempvec)
    print*,'Eigenstates sorted.', cntr, ' eigenvalues were in wrong order.'
    return 

end subroutine sortstates

subroutine sortstates_comp(nev, nest, dim, evals, states)
    implicit none 
    integer, intent(in) :: nev, nest 
    integer(kind=8), intent(in) :: dim 
    double precision, intent(inout) :: evals(nev)
    double complex, intent(inout) :: states(dim, nev) 

    integer, save :: i = 0, j = 0 
    double complex, save :: temp = 0 
    double complex, allocatable, save :: tempvec(:)

    if(allocated(tempvec)) deallocate(tempvec)
    allocate(tempvec(dim))
    do i = 1, nest
        do j = 1, nest
            if (dble(evals(j)) < dble(evals(i)) .and. i < j) then 
                print*,'Swap eigenvalues i, j',i, j
                temp     = evals(i)
                evals(i) = evals(j)
                evals(j) = temp  
                if(dble(evals(i)) > dble(evals(j))) print*,'Error in resorting routine.'
                tempvec = states(i, 1:nest)
                states(i, 1:nest) = states(j, 1:nest)
                states(j, 1:nest) = tempvec 
            end if 
        end do 
    end do 
    deallocate(tempvec)

    return 

end subroutine sortstates_comp

subroutine checnorm(nev, nest, dim, evals, states)
    implicit none 
    integer, intent(in) :: nev, nest 
    integer(kind=8), intent(in) :: dim 

    double precision, intent(in) :: evals(nev), states(dim, nest) 

    integer, save :: i = 0 
    double precision, save :: norm = 0 
    
    do i = 1, nev
        write (* ,"(x, i0, '.Eigenvalue = ',f18.13)") i, dble(evals(i))
        if(i .le. nest) then 
            norm = dot_product(states(1:dim,i),states(1:dim,i))
            if(norm < 0.9999999d0) then
                write (* ,"(x, i0, '. Norm = ',f16.12)") i, dble(norm)
                print*,'Error in subroutine checnorm.'
                error stop  
            end if 
        end if 
    end do     
    print*, 'All eigenstates have norm 1.'
    return 

end subroutine checnorm

subroutine checnorm_comp(nev, nest, dim, evals, states)
    implicit none 
    integer, intent(in) :: nev, nest 
    integer(kind=8), intent(in) :: dim 
    double precision, intent(in) :: evals(nev)
    double complex, intent(in) :: states(dim, nest) 

    integer, save :: i = 0 
    double precision, save :: norm = 0 
    double complex :: psi(dim)
    psi = 0.d0  

    do i = 1, nev
        write (* ,"(x, i0, '.Eigenvalue = ',f18.13)") i, dble(evals(i))
        if(i .le. nest) then 
            ! norm = dot_product(dconjg(states(1:dim,i)), states(1:dim,i))
            norm = dot_product(states(1:dim,i), states(1:dim,i))           
            
            if(norm < 0.9999999d0) then
                write (* ,"(x, i0, '. Norm = ',f16.12)") i, dble(norm)
                print*,'Error in subroutine checnorm.'
                error stop
            end if 
        end if 
    end do     
    print*,'All eigenstates are normalized.'
    return 
    
end subroutine checnorm_comp

subroutine checnorm_comp_feast(nev, nest, dim, evals, states)
    implicit none 
    integer, intent(in) :: nev, nest 
    integer(kind=8), intent(in) :: dim 
    double precision, intent(in) :: evals(nev)
    double complex, intent(in) :: states(dim, nest) 

    integer, save :: i = 0 
    double precision, save :: norm = 0 
    double complex :: psi(dim)
    psi = 0.d0  

    do i = 1, nev
        write (* ,"(i0, '.Eigenvalue = ',f18.13)") i, dble(evals(i))
        if(i .le. nest) then 
            ! norm = dot_product(dconjg(states(1:dim,i)), states(1:dim,i))
            norm = dot_product(states(1:dim,i), states(1:dim,i))           
            
            if(norm < 0.9999999d0) then
                write (* ,"(i0, '. Norm = ',f16.12)") i, dble(norm)
                print*,'Error in subroutine checnorm.'
                error stop
            end if 
        end if 
    end do     
    print*,'All eigenstates are normalized.'
    return 

end subroutine checnorm_comp_feast

subroutine checkorthogonality(nest, dim, states)
    implicit none 
    integer, intent(in) :: nest 
    integer(kind=8), intent(in) :: dim 

    double precision, intent(in) :: states(dim, nest) 

    integer, save :: i = 0, j = 0 
    
    do i = 1, nest 
        do j = 1, nest 
            if(dot_product(states(1:dim, i), states(1:dim, j)) < 0.9999 .and. dot_product(states(1:dim, i), states(1:dim, j)) > 0.000000001) then
                print*, 'Psi(i).Psi(j)', i, j, dot_product(states(1:dim, i), states(1:dim, j))
                error stop 
            end if 
        end do 
    end do 
    print*,'All eigenstates are orthorgonal.'
    return 

end subroutine checkorthogonality

subroutine checkorthogonality_comp(nest, dim, states)
    implicit none 
    integer, intent(in) :: nest 
    integer(kind=8), intent(in) :: dim 

    double complex, intent(in) :: states(dim, nest) 

    integer, save :: i = 0, j = 0, k = 0  
    double complex :: res = 0.d0 
    double complex, allocatable, save :: vec1(:), vec2(:)

    if(allocated(vec1)) deallocate(vec1)
    if(allocated(vec2)) deallocate(vec2)
    allocate(vec1(dim))
    allocate(vec2(dim))
    vec1 = 0.d0
    vec2 = 0.d0 
    do i = 1, nest 
        do j = 1, nest 
            res = 0.d0 
            do k = 1, dim 
                res = res + states(k, i) * dconjg(states(k, j)) 
            end do
            ! print*, 'Psi(i).Psi(j)', i, j, dot_product(states(1:dim, i), states(1:dim, j)) 
            if(abs(dot_product(states(1:dim, i), states(1:dim, j))) < 0.9999 .and. abs(dot_product(states(1:dim, i), states(1:dim, j))) > 0.000000001) then
            ! if(abs(dot_product(dconjg(states(1:dim, i)), states(1:dim, j))) < 0.9999 .and. abs(dot_product(dconjg(states(1:dim, i)), states(1:dim, j))) > 0.000000001) then
                print*, 'Psi(i).Psi(j)', i, j, dot_product(states(1:dim, i), states(1:dim, j))
                ! print*, 'Psi(i)*.Psi(j)', i, j, dot_product(dconjg(states(1:dim, i)), states(1:dim, j))
                print*,'Error in checkorthogonality_comp.'
                error stop 
            end if 
        end do 
    end do 
    print*,'Eigenstates are orthorgonal.'

    return 

end subroutine checkorthogonality_comp

subroutine combinations(m_max, n_max, allcombs, binom) !Calculates all combinations of size m_max drawn from numbers 1-n_max 

    implicit none
    integer, intent(in) :: m_max
    integer, intent(in) :: n_max
    integer, allocatable, intent(out) :: allcombs(:,:) 
    integer(kind=8), intent(out) :: binom
    integer, dimension (m_max) :: comb
    
    ! comb = 0 
    ! character (*), parameter :: fmt = '(i0' // repeat (', 1x, i0', m_max - 1) // ')'

    binom = int(fact(n_max) / (fact(m_max) * fact(max(1,n_max-m_max))),8) !Number of combinations (binomial n_max choose m_max)
    
    if(allocated(allcombs)) deallocate(allcombs)
    allocate(allcombs(binom, m_max))
    call gen (1)
  
  contains
  
    recursive subroutine gen (m)
  
      implicit none
      integer, intent (in) :: m
      integer :: n
      integer :: cntr = 0 
        
      if (m > m_max) then
        write (*, *) comb
        cntr = cntr + 1 
        allcombs(cntr, 1:m_max) = comb 
        if(cntr == binom) cntr = 0 
        ! write (*, fmt) comb
      else
        do n = 1, n_max
          if ((m == 1) .or. (n > comb (m - 1))) then
            comb (m) = n
            call gen (m + 1)
          end if
        end do
      end if
  
    end subroutine gen
  
end subroutine combinations

subroutine permute(value_min, value_max)

    implicit none
    integer, intent(in) :: value_min
    integer, intent(in) :: value_max
    integer :: position_min
    integer :: position_max
    ! integer, dimension (position_min : position_max) :: permutation
    integer, allocatable :: permutation(:)
    
    position_min = value_min
    position_max = value_max
    allocate(permutation(position_min : position_max))
    call generate (position_min)
  
  contains
  
    recursive subroutine generate (position)
  
      implicit none
      integer, intent (in) :: position
      integer :: value
  
      if (position > position_max) then
        write (*, *) permutation
      else
        do value = value_min, value_max
          if (.not. any (permutation (: position - 1) == value)) then
            permutation (position) = value
            call generate (position + 1)
          end if
        end do
      end if
  
    end subroutine generate
  
end subroutine permute


!----------------------------------------!
!            Generate lattice            !
!----------------------------------------!

subroutine lattice(dir, tilted, ucx, ucy, nnbonds, nnnbonds, bc, pattern, cluster, bsites, hexsites, latticevecs, alattice, blattice, xyA, xyB, asitesbonds, bsitesbonds, cntrA, cntrB, nHel, tilt, phase, xy, xtransl, ytransl, reflections, nnnVec)

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
    character :: name*512

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
        
        ! if(allocated(xtransl)) deallocate(xtransl)
            ! if(allocated(ytransl)) deallocate(ytransl)
            ! allocate(xtransl(2, 2 * nUC))
            ! allocate(ytransl(2, 2 * nUC))
            ! xtransl = 0 
            ! ytransl = 0 
        
        !Old cluster
            ! if(cluster == '18A') then 
            !     call cluster_18A(nnbonds, nnnbonds, cntrA, cntrB, nnnVec, bsites, hexsites, alattice, blattice, xtransl, ytransl )
            ! else if(cluster == '18B') then 
            !     call cluster_18B(nnbonds, nnnbonds, cntrA, cntrB, nnnVec, bsites, hexsites, alattice, blattice, xtransl, ytransl )
            ! else if(cluster == '18C') then 
            !     call cluster_18C(nnbonds, nnnbonds, cntrA, cntrB, nnnVec, bsites, hexsites, alattice, blattice, xtransl, ytransl )
            ! else if(cluster == '20A') then 
            !     call cluster_20A(nnbonds, nnnbonds, cntrA, cntrB, nnnVec, bsites, hexsites, alattice, blattice, xtransl, ytransl )
            ! else if(cluster == '20B') then 
            !     call cluster_20B(nnbonds, nnnbonds, cntrA, cntrB, nnnVec, bsites, hexsites, alattice, blattice, xtransl, ytransl )
            ! else if(cluster == '24C') then 
            !     call cluster_24C(nnbonds, nnnbonds, cntrA, cntrB, nnnVec, bsites, hexsites, alattice, blattice, xtransl, ytransl )
            ! else if(cluster == '24D') then 
            !     call cluster_24D(nnbonds, nnnbonds, cntrA, cntrB, nnnVec, bsites, hexsites, alattice, blattice, xtransl, ytransl )
            ! else if(cluster == '30A') then 
            !     call cluster_30A(nnbonds, nnnbonds, cntrA, cntrB, nnnVec, bsites, hexsites, alattice, blattice, xtransl, ytransl )
            ! else if(cluster == '16A') then 
            !     call cluster_16A(nnbonds, nnnbonds, cntrA, cntrB, nnnVec, bsites, hexsites, alattice, blattice, xtransl, ytransl )
            ! else if(cluster == '16B') then 
            !     call cluster_16B(nnbonds, nnnbonds, cntrA, cntrB, nnnVec, bsites, hexsites, alattice, blattice, xtransl, ytransl )
            ! end if
        !
        
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

end subroutine lattice

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
    character :: filename*256
    character :: parameters*256

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
    character :: name*512

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
    write (name,"('L=',i0,'ucx=',i0, &
        'ucy=',i0,'pat=',a2 &
        ,'BC=',a,'.dat')") nn,ucx,ucy,pattern,&
        bc
    name = trim_name(name)
    name=dir//"NN_lattice_"//name
    open(unit=31, file=name)
    do s = 1, nbonds
        write(31,*) bsites(1,s), bsites(2,s)
    end do
    close(31)

    write (name,"('L=',i0,'ucx=',i0, &
        'ucy=',i0,'pat=',a2 &
        ,'BC=',a,'.dat')") nn,ucx,ucy,pattern,bc
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
    
    character :: name*512

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

    write (name,"('L=',i0,'ucx=',i0, &
        'ucy=',i0,'pat=',a2 &
        ,'BC=',a,'.dat')") nn,ucx,ucy,pattern,bc
    name = trim_name(name)
    name=dir//"NNN_lattice_"//name
    open(unit=32, file=name)
    do s = 1, nbonds
        write(32,*) hexsites(1,s), hexsites(2,s)      
    end do

    close(32)

    write (name,"('L=',i0,'ucx=',i0, &
        'ucy=',i0,'pat=',a2 &
        ,'BC=',a,'.dat')") nn,ucx,ucy,pattern,bc
    name = trim_name(name)
    name=dir//"NNN_xy_"//name
    open(unit=32, file=name)
    do s = 1, nbonds
        write(32,'(4(i0,x))') xy(1,s), xy(2,s), xy(3,s), xy(4,s)
    end do
    close(32)

    ! print*, 'Finished honeycomb nnn-lattice'

end subroutine honeycomb_nnn

subroutine honeycomb_nnn2(dir, ucx, ucy, bc, pattern, nbonds, hexsites, latticevecs, alattice, blattice, asitesbonds, bsitesbonds, cntrA, cntrB, phase, xy, xtransl, ytransl, sitecoord)
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
    
    character :: name*512

    !Definition of lattice vectors: 
    ! v1 = {sqrt(3), 0}
    ! v2 = 1/2*{-sqrt(3), 3}
    ! v3 = 1/2*{-sqrt(3), -3}

    print*, 'Start creating Honeycomb sublattices'
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

    if ( ucx == 2 ) then
        nx = 2 * ucx - 2!nx = number of bonds per x-layer (PBC)
        ny = ucy !Number of y layers
    else 
        nx = 2 * ucx !nx = number of bonds per x-layer (PBC)
        ny = ucy !Number of y layers
    end if


    if(bc == 'o') then
        nbx = 2 * (ucx - 1) * ucy !total number of x bonds: 2 = # sublattices
        nby = 2 * (ucy - 1) * ( 2 * (ucx-1) + 1 ) !total number of y bonds: 2 = # sublattices; ucx-1 = unit cells with 2 bonds/site; 1 = unit cell with 1 bond per site (first uc)
    else if(bc == 'p') then
            nbx = nx * ny !total number of x bonds
            nby = 2 * nx * ny !total number of y bonds
    end if
    nbonds = nbx + nby
    print*, 'Number of NNN bonds', nbonds 

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
    
    ! allocate(xtransl(2, nn))
    allocate(xtransl(2, nbx))
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

    write (name,"('L=',i0,'ucx=',i0, &
        'ucy=',i0,'pat=',a2 &
        ,'BC=',a,'.dat')") nn,ucx,ucy,pattern,bc
    name = trim_name(name)
    name=dir//"NNN_lattice_"//name
    open(unit=32, file=name)
    do s = 1, nbonds
        write(32,*) hexsites(1,s), hexsites(2,s)      
    end do

    close(32)

    write (name,"('L=',i0,'ucx=',i0, &
        'ucy=',i0,'pat=',a2 &
        ,'BC=',a,'.dat')") nn,ucx,ucy,pattern,bc
    name = trim_name(name)
    name=dir//"NNN_xy_"//name
    open(unit=32, file=name)
    do s = 1, nbonds
        write(32,'(4(i0,x))') xy(1,s), xy(2,s), xy(3,s), xy(4,s)
    end do
    close(32)



    print*, 'Finished honeycomb nnn-lattice'

end subroutine honeycomb_nnn2


!--------------------------!
!            BASIS         !
!--------------------------!

subroutine characters(symmetrize, irrep, mir, rot, id)
    implicit none
    integer, intent(in) :: symmetrize
    character(len=2), intent(in) :: irrep
    double precision, intent(out) :: mir(6), rot(5), id

    !rot(1) = c6 
    !rot(2) = (c6)^2 = c3  
    !rot(3) = (c6)^3 = c2
    !rot(4) = (c6)^4 = -c3 
    !rot(5) = (c6)^5 = -c6
    if(symmetrize == 0) then 
        id = 1
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
        id  = 1
        mir = 1
        mir(4) = -1
        mir(5) = -1
        mir(6) = -1
        rot = - 1
        rot(2) = 1
        rot(4) = 1
    else if(irrep == "B2") then     
        id  = 1
        mir = - 1
        mir(4) = 1
        mir(5) = 1
        mir(6) = 1
        rot = - 1
        rot(2) = 1
        rot(4) = 1
    else if(irrep == "E1") then 
        id  = 2
        mir = 0 
        rot(1) = +1 
        rot(5) = +1
        rot(2) = -1
        rot(4) = -1
        rot(3) = -2
    else if(irrep == "E2") then 
        id  = 2
        mir = 0 
        rot(1) = -1 
        rot(5) = -1
        rot(2) = -1
        rot(4) = -1
        rot(3) = +2
    end if 

end subroutine characters

subroutine make_basis(ti, tilted, pat, nnnVec, sites, particles, dim, symmetrize, ucx, ucy, l1, l2, basis, abasis, bbasis, tilt, nHel, k1, k2, xtransl, ytransl, id, par, rot, refl, c6, period, norm, orbsize, orbits2D, phases2D, norm2D)

    implicit none

    integer, intent(in) :: ti, tilted, sites, particles, tilt, nHel, ucx, ucy, k1, k2, symmetrize
    integer, allocatable, intent(in) :: xtransl(:,:), ytransl(:,:)
    double precision, intent(in) :: nnnVec(2,3), id, par(6), rot(5)
    character(len=*), intent(in) :: pat

    integer, intent(out) :: l1, l2, orbsize
    integer(kind=8), intent(out) :: dim
    integer(kind=8), allocatable, intent(out) :: basis(:), abasis(:), bbasis(:),period(:), orbits2D(:,:,:)
    integer, allocatable, intent(out) :: refl(:,:), c6(:)
    double precision, allocatable, intent(out) :: norm(:), norm2D(:,:)
    double complex, allocatable, intent(out) :: phases2D(:,:,:)
    
    integer :: a(sites), amask, bmask, temp = 0, I_in = 0, I_tail = 0, i, j, l
    integer(kind=8) :: symdim
    integer(kind=8), allocatable :: symbasis(:)    

    dim = int(fact(sites) / (fact(particles) * fact(max(1,sites-particles))),8) !Number of basis states


    if (dim == 0) then
        print*, 'No basis states available.'
        return 
    end if

    !Permutations contains integer values I of basis states, perm_up/dn contain the integer values I_up/dn and the corresponDing J_up/dn
    if (allocated(basis)) deallocate(basis)
    allocate(basis(dim))
    if (allocated(abasis)) deallocate(abasis)
    allocate(abasis(dim))
    if (allocated(bbasis)) deallocate(bbasis)
    allocate(bbasis(dim))
    basis  = 0
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
        basis(j)  = I_in
        abasis(j) = iand( basis(j), amask ) 
        bbasis(j) = iand( basis(j), bmask )
        temp = I_in
        do i = 0, 64
            if (btest(temp,i)) then
                temp = ibclr(temp,i)
                if (not(btest(temp,i+1))) then   ! asks if pos i+1 is zero
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
    call symm_basis(tilted, dim, sites, nHel, l2, l1, k2, k1, symmetrize, id, par, rot, nnnVec, basis, xtransl, ytransl, refl, c6, symdim, symbasis, period, norm, orbsize, orbits2D, phases2D, norm2D)
    
    dim = symdim
    deallocate(basis)
    allocate(basis(dim))
    basis = symbasis 
    print*, 'Basis symmetrized.'

end subroutine make_basis



!------------------------------------------------!
!            Lookup of basis states              !
!------------------------------------------------!

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

!---------------------------------------------!
!            Momentum state basis             !
!---------------------------------------------!

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
                ! call checkstate_rect(basis(j), n, Lx, Ly, kx, ky, xtransl, ytransl, nnnVec, r, normalization)
                call checkstate_rect(basis(j), n, Lx, Ly, kx, ky, irrep, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, r, normalization)
            else if(tilted == 1) then 
                call checkstate(basis(j), n, nHel, Lx, Ly, kx, ky, irrep, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, r, normalization)
                ! call checkstate2(basis(j), n, nHel, 1, Lx, Ly, kx, ky, 0, 0, refl, xtransl, ytransl, nnnVec, r, normalization)
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

        ! call checkstate5(basis(j), n, nHel, Lx, Ly, kx, ky, irrep, id, par, rot, refl, c6, xtransl, ytransl, nnnVec, r, normalization)

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

!-------------------------------------------------------------------!
!            Translation representatives and periodicity            !
!-------------------------------------------------------------------!
    
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
        integer :: orbit = 1, i, j, n, rt, c
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
             
            ! if(s0 == 39326 .and. c == 2) then 
                    
            !     print* ,phases(1:orbsize, 1), 'phases(1:orbsize, 1)'
            !     print*, 'before translation'
            !     print* ,s0, 's0, checkst2D'
            !     print* ,c, 'c'
            !     print* ,info, 'info'

            !     pause 
            ! end if
            call translation2D(s0, s, sites, ntot, orbsize, orbit, orbits(1:orbsize, c), tilted, layers, Lx, Ly, 1.0d0, sign, a1, a2, xtransl, ytransl, k, phases(1:orbsize, c), phase, info)

            ! if(s0 == 39326 .and. c == 2) then 
                    
            !     print* ,phases(1:orbsize, 1), 'phases(1:orbsize, 1)'
            !     print*, 'after translation'
            !     print* ,s0, 's0, checkst2D'
            !     print* ,c, 'c'
            !     print* ,info, 'info'

            !     pause 
            ! end if

            ! if(s0 == 39326) then 
                
            !     print* ,phases(1:orbsize, c), 'phases(1:orbsize, c)'
            !     print*, 'translation'
            !     print* ,s0, 's0, checkst2D'
            !     print* ,info, 'info'
            !     pause 
            ! end if

            if(info < 0) phases(1:orbsize, c) = 0.d0 !All prefactors = 0.  
            if(info < 0) orbits(1:orbsize, c) = 0    !Remove orbit. 
            ! if(s0 == 39326 .and. info<0) print* ,r, 'r, transl'
            ! if(s0 == 39326 .and. info<0) pause
            if(info < 0) cycle
            ! if(info < 0) r = -2
            
            
            ! if(info < 0 .and. (s0 == 174761 .or. s0 == 87382 .or. s0 == 39326 .or. s0 == 39326 .or. s0 == 174762) ) then 
            !     print* , 'translations'
            !     print* , 's0' ,s0
            !     print* , 'c' ,c
            ! end if 

            
            !Reflections 
            reflection = 0.d0 
            do i = 1, 6
                sign = 1
                reflection = identity 

                ! do j = 1, i - 1
                !     reflection = matmul(reflection, rho)
                ! end do        
                ! reflection = matmul(reflection, sigma)
                
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
                !debug: rotation(1, c) ??? 
                
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
            norm(c) = rt * abs(phase)**2! * 0.5 !0.5 = 1/(dim of irrep) 
            dcntr = 0.d0 
            do i = 1, orbsize
                dcntr = dcntr + abs(phases(i, c))**2
            end do 
            ! norm(c) = dcntr * 0.5
            print* ,s0, 's0'
            print* ,phase, 'phase'
            print* ,r, 'r'
            print* ,norm(c), 'norm(c)'
            print* ,dcntr, c, 'dcntr, c'
            pause 
        end do 



        if(r == -2) r = -1 
       
        return


   
end subroutine checkstate2D

!Compares translated state to reflected or rotated state. Does check for phase = 0.
subroutine translation2(s0, s, sites, ntot, orbit, orbits, nHel, Lx, Ly, char, signpg, a1, a2, xtransl, ytransl, k, phase, info)
    
    implicit none 
    integer, intent(in) :: sites, ntot, nHel, Lx, Ly, char, signpg
    integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
    integer(kind=8), intent(in) :: s0, s
    double precision, intent(in) :: a1(2), a2(2), k(2)
    integer, intent(out) :: info  
    integer, intent(inout) :: orbit
    integer, allocatable, intent(inout) :: orbits(:)
    double complex, intent(inout) :: phase

    integer :: nx, ny, i, rowst, edgeA, edgeB, flag, orbit0
    integer(kind=8) :: t, tx, ty, sign
    double precision :: shift(2)
    double complex :: phaset
    double precision, parameter :: tolerance = 1.0e-8

    
    i = 0
    t = 0 
    tx = 0 
    ty = 0
    nx = 0
    ny = 0
    info = 0 !On output: InDicates whether new state was found (info = 0) or state is already contained in list (info < 1)
    sign = 1
    edgeA = 0
    edgeB = 0
    rowSt = 0
    !Include initial state (if not already in the list) to orbit 
    orbit = orbit + 1 
    orbit0 = orbit 
    flag  = 1 
    ! phaset = phase + signpg * char 
    phaset = 0.d0 
    ! phaset = 0.d0 


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
            do i = 1, nHel!Ly !Calculate anticommutation sign from x-translations 
                edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                if(btest(tx, edgeB) .and. btest(tx, edgeA)) then !If both edges occupied, no sign
                    cycle 
                else if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                    rowSt = (i-1)*2*Lx !First site on y-layer 'i' 
                    sign  = sign * (-1)**popcnt( ibits( tx, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'
                end if 
                ! edgeB = (i-1)*2*Lx + 1 !First B site on y-layer 'i'
                ! edgeA = (i-1)*2*Lx     !First A site on y-layer 'i'
                ! if(btest(tx, edgeB) .and. btest(tx, edgeA)) then !If both edges occupied, no sign
                !     cycle 
                ! else if(btest(tx, edgeB) .neqv. btest(tx, edgeA)) then !If only one edge occupied 
                !     rowSt = edgeB + 1!First site on y-layer 'i' 
                !     sign  = sign * (-1)**popcnt( ibits( tx, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'
                ! end if 
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
            shift = nx * a2  
            phaset = phaset + exp( - ii * dot_product(k, shift)) * sign * char * signpg 

            if(t == s) then 
                ! shift = nx * a2  
                ! phase = phase + exp( - ii * dot_product(k, shift)) * sign * char * signpg 
                phase = phase + phaset 
                phaset = 0.d0 
                
            end if 
            if (t < s0) then !Representative is already in the list 
                info = -1 
                return
            end if

        ! if(nHel > 1) then 
        else if(nHel > 1) then 
            do ny = 1, Ly 
                ty = t
                ! sign = sign * ((-1)**(ntot-popcnt( ibits( ty, 2*Lx*(nHel-1), 2*Lx ) ) ))**popcnt( ibits( ty, 2*Lx*(nHel-1), 2*Lx ) ) 
           
                ! sign = sign * ((-1)**modulo(popcnt( ibits( ty, 2*Lx*(nHel-1), 2*Lx )), 2))**modulo(ntot-popcnt( ibits( ty, 2*Lx*(nHel-1), 2*Lx ) ) , 2)  

                sign = sign * ((-1)**modulo(popcnt( ibits( ty, 0, 2*Lx )), 2))**modulo(ntot-popcnt( ibits( ty, 0, 2*Lx ) ) , 2) 
                if( ((-1)**modulo(popcnt( ibits( ty, 0, 2*Lx )), 2))**modulo(ntot-popcnt( ibits( ty, 0, 2*Lx ) ) , 2) == -1) stop 'minus sign' 

                ! sign = sign * ((-1)**(ntot-popcnt( ibits( ty, 0, 2*Lx ) ) ))**popcnt( ibits( ty, 0, 2*Lx ) ) 
                do i = 1, nHel
                    ! edgeB = i*2*Lx-1 !Last B site on y-layer 'i'
                    ! edgeA = i*2*Lx-2 !Last A site on y-layer 'i'
                    edgeB = (i-1)*2*Lx + 1 !First B site on y-layer 'i'
                    edgeA = (i-1)*2*Lx     !First A site on y-layer 'i'
                    if(btest(ty, edgeB) .and. btest(ty, edgeA)) then !If both edges occupied, no sign
                        cycle 
                    else if(btest(ty, edgeB) .neqv. btest(ty, edgeA)) then !If only one edge occupied 
                        ! rowSt = (i-1)*2*Lx!First site on y-layer 'i' 
                        rowSt = edgeB + 1!First site on y-layer 'i' 
                        sign = sign * (-1)**popcnt( ibits( ty, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'

                        ! if(s0 == 39326 .or. s0 == 6559 ) then 
                    
                        !     print* ,s0, ty, 's0, ty'
                        !     write(*,'(X,B0)') ty
                        !     print* ,edgeA, 'edgeA'
                        !     print* ,edgeB, 'edgeB'
                        !     print* ,rowSt, 'rowSt'
                        !     print* ,2*Lx-2, '2*Lx-2'
                        !     print* ,popcnt( ibits( ty, rowst, 2*Lx-2 ) ), 'popcnt( ibits( ty, rowst, 2*Lx-2 ) )'
                        !     write(*,'(X,B0)') ibits( ty, rowst, 2*Lx-2 )
                        !     print* ,(-1)**popcnt( ibits( ty, rowst, 2*Lx-2 ) ), '(-1)**popcnt( ibits( ty, rowst, 2*Lx-2 ) )'
                        !     print* ,sign, 'sign'
                            
                        ! end if 
                    end if                 
                    
                end do     


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

                shift = ny * a1 + nx * a2  
                ! phaset = phaset + exp( - ii * dot_product(k, shift)) * sign * char * signpg 
                phaset = phaset + exp( - ii * dot_product(k, shift)) * sign * char * signpg 
                if(t == s) then 
                    ! shift = ny * a1 + nx * a2  
                    ! phase = phase + exp( - ii * dot_product(k, shift)) * sign * char * signpg 
                    ! phase  = phase + phaset 
                    ! phaset = 0.d0
                    ! sign   = 1
                    if((abs(dble(phase)) < tolerance) .and. (abs(aimag(phase)) < tolerance)) then
                        print* ,s0, 's0'
                        print*, 'zero phase'
                        pause 
                    end if 
                end if 
                if (t < s0) then !Representative is already in the list 
                    info = -1
                    return
                end if
                
            end do
            
        end if 

    end do

    
    if((abs(dble(phase)) < tolerance) .and. (abs(aimag(phase)) < tolerance)) then 
        orbits(orbit0:orbit) = 0 
        print*, 'FLUSHED ORBIT ---------------------------------------------'
    end if 
    if(t .ne. s ) print*, 'Translation orbit incomplete.'
    if(t .ne. s ) stop

end subroutine translation2

!Compares translated state to initially tested basis state, NOT the reflected or rotated state. Does not check for phase = 0.
subroutine translation(s0, s, sites, ntot, orbit, orbits, tilted, nHel, Lx, Ly, char, signpg, a1, a2, xtransl, ytransl, k, phase, info)
        !If tilted == 0, swap input vectors a1 and a2 
        implicit none 
        integer, intent(in) :: sites, ntot, tilted, nHel, Lx, Ly, signpg
        integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
        integer(kind=8), intent(in) :: s0, s
        double precision, intent(in) :: a1(2), a2(2), k(2), char
        integer, intent(out) :: info  
        integer, intent(inout) :: orbit
        integer, allocatable, intent(inout) :: orbits(:)
        double complex, intent(inout) :: phase

        integer :: nx, ny, i, rowst, edgeA, edgeB, flag
        integer(kind=8) :: t, tx, ty, sign 
        double precision :: shift(2)

        
        i = 0
        t = 0 
        tx = 0 
        ty = 0
        nx = 0
        ny = 0
        info = 0 !On output: InDicates whether new state was found (info = 0) or state is already contained in list (info < 1)
        sign = 1
        edgeA = 0
        edgeB = 0
        rowSt = 0
        !Include initial state (if not already in the list) to orbit 
        orbit = orbit + 1 
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

subroutine translation2D(s0, s, sites, ntot, orbsize, orbit, orbits, tilted, nHel, Lx, Ly, char, signpg, a1, a2, xtransl, ytransl, k, phases, phase, info)

    !If tilted == 0, swap input vectors a1 and a2 and set nHel = Ly 
    implicit none 
    integer, intent(in) :: orbsize, sites, ntot, tilted, nHel, Lx, Ly, signpg
    integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
    integer(kind=8), intent(in) :: s0, s
    double precision, intent(in) :: a1(2), a2(2), k(2), char
    integer, intent(out) :: info  
    integer, intent(inout) :: orbit
    integer(kind=8), intent(inout) :: orbits(orbsize)
    double complex, intent(inout) :: phases(orbsize), phase 

    integer :: nx, ny, i, rowst, edgeA, edgeB, flag
    integer(kind=8) :: t, tx, ty, sign 
    double precision :: shift(2)

        t = 0 
        tx = 0 
        ty = 0
        nx = 0
        ny = 0
        info = 0 !On output: InDicates whether new state was found (info = 0) or state is already contained in list (info < 1)
        sign = 1
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

subroutine xtranslate(s, nHel, n, Lx, xtransl, sign, sx)
    implicit none 
    integer(kind=8), intent(in) :: s 
    integer, intent(in) :: nHel, n, Lx
    integer, intent(in) :: xtransl(2, n)
    integer, intent(inout) :: sign 
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
    integer, intent(inout) :: sign 
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

!Uses modulo in calculating sign changes
subroutine ytranslate2(s, nHel, tilt, n, Lx, ytransl, sign, sy)
    implicit none 
    integer(kind=8), intent(in) :: s 
    integer, intent(in) :: nHel, tilt, n, Lx
    integer, intent(in) :: ytransl(2, n)
    integer, intent(inout) :: sign 
    integer(kind=8), intent(out) :: sy

    integer :: i = 0, ntot = 0, edgeA = 0, edgeB = 0, rowSt = 0
    
    sy = 0
    ntot = popcnt(s)
    if(nhel == 1) then !Single helix
        ! sign = sign * ((-1)**(ntot-popcnt( ibits( s, n - tilt, tilt ) ) ))**popcnt( ibits( s, n - tilt, tilt ) ) 
        sign = sign * ((-1)**(modulo(ntot-popcnt( ibits( s, n - tilt, tilt)), 2)))**modulo(popcnt(ibits(s, n - tilt, tilt)), 2)
    else 
        sign = sign * ((-1)**modulo(ntot-popcnt(ibits(s, 0, 2*Lx)), 2))**modulo(popcnt(ibits(s, 0, 2*Lx)), 2)
        ! sign = sign * ((-1)**(ntot-popcnt(ibits(s, 0, 2*Lx))))**popcnt(ibits(s, 0, 2*Lx))
        do i = 1, nHel
            edgeA = (i-1)*2*Lx !First A site on y-layer 'i'
            edgeB = (i-1)*2*Lx + 1 !First B site on y-layer 'i'
            if(btest(s, edgeB) .and. btest(s, edgeA)) then !If both edges occupied, no sign
                cycle 
            else if(btest(s, edgeB) .neqv. btest(s, edgeA)) then !If only one edge occupied 
                rowSt = (i-1)*2*Lx + 2 !First site after edge A on y-layer 'i' 
                sign = sign * (-1)**popcnt( ibits( s, rowst, 2*Lx-2 ) ) !Occupation of y-layer 'i'
            end if 
        end do     
    end if 

    do i = 1, n
        call mvbits(s, ytransl(1,i)-1, 1, sy, ytransl(2,i)-1)
    end do

end subroutine ytranslate2

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

!Uses ytranslate2 for checkstate5 
subroutine representative_irrep2(s, n, nHel, tilt, Lx, Ly, symmetrize, id, par, rot, xtransl, ytransl, refl, c6, r, l1, l2, sign)
    !Finds the representative 'r' for state 's'. 'n' is the number of sites and 'l' the number of shifts needed to translate 's' to 'r'.
    implicit none
    integer(kind=8), intent(in) :: s
    integer, intent(in) :: n, nHel, tilt, symmetrize, Lx, Ly, xtransl(2, n), ytransl(2, n), refl(6, n), c6(n)
    double precision, intent(in) :: id, par(6), rot(5)
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
        if(nHel == 1) then 
            if (t < r) then !Why twice?
                sign = signt 
                r    = t
                l1   = 0
                l2   = nx
            end if
        else if(nHel > 1) then 
            do ny = 1, Ly 
                ty = t                
                !Note, 'signt' and 't' are carried on through out all translations. 
                call ytranslate2(ty, nHel, tilt, n, Lx, ytransl, signt, t)
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
        call mirror_rep2(r, s, n, nHel, tilt, par(i), Lx, Ly, refl(i, 1:n), xtransl, ytransl, sign, l1, l2)
    end do 
    do i = 1, 5 !Check for representatives among rotated states
        call rotate_rep2(r, s, n, nHel, tilt, i, rot(i), Lx, Ly, c6, xtransl, ytransl, sign, l1, l2)
    end do 
    12 continue

    return 

end subroutine representative_irrep2

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

subroutine mirror_rep2(r, s, n, nHel, tilt, p, Lx, Ly, refl, xtransl, ytransl, sign, l1, l2)
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
                call ytranslate2(ty, nHel, tilt, n, Lx, ytransl, signt, t)                
                if (t < r) then
                    sign = signt * signp * p
                    r    = t
                    l1   = ny
                    l2   = nx
                end if
            end do
        end if 
    end do

end subroutine mirror_rep2

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

subroutine rotate_rep2(r, s, n, nHel, tilt, nrot, eval, Lx, Ly, c6, xtransl, ytransl, sign, l1, l2)
    implicit none
    integer(kind=8), intent(in) :: s !original state
    integer, intent(in) :: n, nHel, tilt, nrot, Lx, Ly, c6(n), xtransl(2,n), ytransl(2,n)     
    double precision, intent(in) :: eval
    integer(kind=8), intent(inout) :: r !potential rep (r<=s) 
    integer, intent(inout) :: l1, l2 
    integer, intent(inout) :: sign !If 'r' remains unchanged, previous sign remains unchanged. Otherwise updated. 

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
                call ytranslate2(ty, nHel, tilt, n, Lx, ytransl, signt, t)                
                if (t < r) then
                    sign = signt * eval * signrot
                    r    = t
                    l1   = ny
                    l2   = nx
                end if
            end do
        end if 
    end do

end subroutine rotate_rep2

subroutine reflect_rep(r, s, n, nHel, tilt, p, Lx, Ly, refl, xtransl, ytransl, sign, l1, l2)
    implicit none
    integer(kind=8), intent(in) :: s
    integer, intent(in) :: n, nHel, tilt, Lx, Ly, refl(n), xtransl(2,n), ytransl(2,n)     
    double precision, intent(in) :: p
    integer(kind=8), intent(inout) :: r
    integer, intent(inout) :: l1, l2, sign

    integer :: nx = 0, ny = 0, signp = 1, signt = 1, info = 0
    integer(kind=8) :: t = 0, tx = 0, ty = 0 

    signp = 1
    signt = 1
    call reflect(r, s, n, refl, signp, info, t) ! compare to r, reflect s onto t
    if (t < r) then
        sign = signp * p 
        r    = t
        l1   = 0
        l2   = 0
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

end subroutine reflect_rep
    








!-------------------------------------------------------------------------!
!            Calculate number of eigenvalues and Lanczos vectors          !
!-------------------------------------------------------------------------!

subroutine nevncv(unit, parameters, thresh, exact, nevext, nestext, ncv0, dim, full, nev, ncv, nest)

    implicit none
    integer, intent(in) :: unit, thresh, exact, nevext, nestext, ncv0
    integer(kind=8), intent(in) :: dim
    character(len=*), intent(in) :: parameters
    integer, intent(out) :: full, nev, ncv, nest
    character*512 :: filename

    if (exact == 0) then
        if (dim < thresh) then
            full = 1
            nev = dim
            ncv = dim
            nest = dim
        else
            full = 0
            if (dim == 1) then
                nev = 1
                ncv = 1
                nest = 1
            else
                ncv = min(ncv0, dim - 1) !nev+1<=ncv<=dim and ncv maximal
                ! ncv = max(ncv0, 2*nev + 10) !nev+1<=ncv<=dim and ncv maximal
                nev = min(nevext, ncv - 3)
                ! nev = min(nevext, dim - 10)
                ! if (ncv > dim) ncv = dim
                nest = min(nev, nestext)
            end if
        end if
    else if (exact == 1) then
        full = 1
        nev = dim
        ncv = dim
        nest = min(dim, nestext)
    end if
    print*, 'Number of eigenvalues = ', nev
    print*, 'Number of eigenvectors = ', nest
    print*, ''

    ! filename=dir//"Hilbert_space_dimension_"//parameters
    ! open(unit, file=filename)
    ! write(unit,*) dim 
    ! close(unit)

end subroutine nevncv

!------------------------------------------!
!            Set initial variables         !
!------------------------------------------!

subroutine setvars()

    use variables
    implicit none

    character*10:: clusterl
    integer :: threads

    !call KMP_SET_STACKSIZE_S(1000000000)

    !$ call omp_set_dynamic(dynamic)
    !$ call omp_set_nested(nested)
    !$ call omp_set_num_threads(othrds)
    call mkl_set_num_threads(mthrds)
    call mkl_get_max_threads(threads)

    print*,'MKL THREADS',threads
    !$ threads = omp_get_num_threads()
    print*,'OMP THREADS',threads

    if(tilted == 0) then 
        sites = 2*ucx*ucy 
    else if(cluster(1:2) == "16") then 
        sites = 16 
    else if(cluster(1:2) == "18") then 
        sites = 18 
    else if(cluster(1:2) == "20") then 
        sites = 20 
    else if(cluster(1:2) == "24") then 
        sites = 24
    else if(cluster(1:2) == "30") then 
        sites = 30 
    else if(cluster(1:2) == "32") then 
        sites = 32
    else if(cluster(1:1) == "L") then 
        clusterl = cluster
        read(clusterl(2:3),*) sites
    end if 
    particles = int(filling * sites)
    

    if(symmetrize == 1) print('(1x, a,a)'), 'Irrep = ', irrep
    if(tilted == 1) then
        print('(1x, a,a)'), 'Cluster = ', cluster
    else 
        print('(1x, a,i0)'), 'UCX = ', ucx
        print('(1x, a,i0)'), 'UCY = ', ucy
    end if 
    print('(1x, a,i0)'), 'Sites = ', sites
    print('(1x, a,i0)'), 'Particles = ', particles
    print*, ''
    

    return

end subroutine setvars

!------------------------------------------!
!               Hamiltonian                !
!------------------------------------------!

subroutine hamiltonian(spartan, thrds, ti, unit, prms, sites, nbb, trisites, dim, basis, hamOff, nOff, k1, k2, tilted, nHel, tilt, lx, ly, ucx, ucy, orbsize, norm, norm2D, orbits, phases, xtransl, ytransl, symmetrize, id, par, rot, refl, c6, t, rcoff, rcdi, prts, dplcts, hamOff_dp, hamDi_d, nDi_off, hamOff_dc, hamDi_c, hamDi_c2, nnnbb, hexsites, hamDi, occ, nDi)

    implicit none 

    integer, intent(in) :: spartan, thrds, ti, ucx, ucy, tilted, nHel, tilt, lx, ly, k1, k2, symmetrize, orbsize
    double precision, intent(in) :: t, id, par(6), rot(5) 
    integer, intent(inout) :: unit, sites, nbb, nnnbb, nOff, nDi_off, nDi
    integer(kind=8), intent(inout) :: dim 
    integer(kind=8), allocatable, intent(inout) :: basis(:), orbits(:,:,:) 
    integer, allocatable, intent(inout) :: hamDi(:,:), occ(:,:), prts(:), hamOff(:,:), xtransl(:,:), ytransl(:,:), dplcts(:), trisites(:,:), hexsites(:,:), rcoff(:,:), rcdi(:)
    integer, intent(in) :: refl(6, sites), c6(sites)
    double precision, allocatable, intent(inout) :: norm(:), norm2D(:,:), hamDi_dp(:), hamOff_dp(:)
    double complex, allocatable, intent(inout) :: hamOff_dc(:), hamDi_c(:), hamDi_c2(:,:), phases(:,:,:)
    character(len=*), intent(inout) :: prms

    character :: type*1 
    logical :: exist1, exist2, exist3

    print*, 'Creating Hamiltonian ...'
    print*, ''

    call par2(unit, k1, k2, -1, 0.d0, 0.d0, prms, type)
    ! call loadham(spartan, type, prms, unit, ti, sites, nOff, nDi, hamDi, hamOff, rcoff, hamOff_dp, hamOff_dc, occ, exist1, exist2, exist3)
    exist1 = .false.
    exist2 = .false.
    exist3 = .false.
    if(.not.(exist1)) then 
        if(ti == 0) then 
            call hopping_p2(thrds, unit, prms, sites, nbb, trisites, dim, basis, hamOff, nOff)
        else if(ti == 1 .and. k1 == 0 .and. k2 == 0) then 
            ! call hopping_symm(thrds, tilted, nHel, tilt, lx, ly, ucx, ucy, sites, nbb, dim, basis, trisites, norm, xtransl, ytransl, symmetrize, id, par, rot, refl, c6, t, rcoff, rcdi, ptrs, dplcts, hamOff_dp, hamDi_d, nOff, nDi_off)
            if(id == 1) call hopping_irrep(thrds, tilted, nHel, tilt, lx, ly, ucx, ucy, sites, nbb, dim, basis, trisites, norm, xtransl, ytransl, symmetrize, id, par, rot, refl, c6, t, rcoff, rcdi, prts, dplcts, hamOff_dp, hamDi_d, nOff, nDi_off)
            if(id == 2) then 
                call hopping_irrep2D(thrds, tilted, nHel, tilt, lx, ly, sites, nbb, dim, basis, orbsize, orbits, phases, norm2D, trisites, xtransl, ytransl, symmetrize, id, par, rot, refl, c6, t, rcoff, rcdi, prts, dplcts, hamOff_dc, hamDi_c, nOff, nDi_off)
                call diagonal_irrep2D(sites, dim, trisites, hexsites, orbsize, orbits, phases, norm2D, hamDi_c2, occ)
            end if 
        else if(ti == 1 .and. ((k1 .ne. 0) .or. (k2 .ne. 0))) then 
            call khopping_c(unit, prms, tilted, nHel, tilt, lx, ly, ucx, ucy, sites, nbb, dim, basis, trisites, norm, xtransl, ytransl, t, k1, k2, rcoff, hamOff_dc, nOff)
        end if 
    end if
    

    if(id == 1 .and. (.not.(exist2) .or. .not.(exist3))) call diagonal(unit, prms, sites, nbb, nnnbb, dim, trisites, hexsites, basis, hamDi, occ, nDi) 
    
    if(id == 1 .and. (.not.(exist1) .or. .not.(exist2) .or. .not.(exist3))) call saveham(spartan, type, prms, unit, ti, sites, nOff, nDi, hamDi, hamOff, rcoff, hamOff_dp, hamOff_dc, occ)

    
    ! if(parity == 1) then     
    !     open(unit, file=trim_name("hamiltonians/diagonal_hamiltonian_2_"//prms))
    !     write(unit,*) hamDi_d
    ! end if 

    return 

end subroutine hamiltonian

subroutine hopping_p2(threads, unit, parameters, sites, nbonds, bsites, dim, basis, ham, nnz)

    implicit none

    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: threads, unit, sites, nbonds
    integer, intent(in) :: bsites(2, nbonds)
    character(len=*), intent(in) :: parameters
    integer, intent(out) :: nnz
    integer, allocatable, intent(out) :: ham(:,:)

    integer(kind=8) :: arrsize = 0, cntr2 = 0, loc = 0, pos = 0 
    integer(kind=8) :: i = 0, j = 0, s = 0
    integer(kind=8) :: newst = 0 
    integer :: cntrj = 0, id = 0, parity1 = 0, parity2 = 0
    integer(kind=8), allocatable :: cntr(:)
    integer, allocatable :: ham_temp(:,:)
        
    
    arrsize = 2*nbonds*dim
    
    if (allocated(ham_temp)) deallocate(ham_temp)
    if (allocated(cntr)) deallocate(cntr)
    allocate(ham_temp(arrsize, 3))
    allocate(cntr(threads))
    
    ham_temp = 0 
    cntr     = 0
    nnz      = 0

    !$omp parallel do default(firstprivate) shared(ham_temp, cntr, bsites, basis) num_threads(threads)
    do j = 1, dim
        !$ id = omp_get_thread_num()
        cntrj = 0 
        do s = 1, nbonds    

            if (btest(basis(j), bsites(1, s) - 1) .and. .not. btest(basis(j), bsites(2, s) - 1)) then 
                newst = ibclr( ibset( basis(j), bsites(2, s) - 1 ), bsites(1, s) - 1 ) !Create on site 2, annihilate on site 1 
                call findstate(dim, newst, basis, loc) 
                if(loc > 0) then 
                    cntr(id+1) = cntr(id+1) + 1 !Count nnz elements in each thread 
                    cntrj = cntrj + 1 !Count nnz elements per basis state 
                    pos = (j-1)*nbonds + cntrj !Position of nnz element in sparse array: Offset ((j-1)*nbonds): Each basis state (j) has at most nbonds nnz elements. Cntrj gives the current increment.
                    parity1 = popcnt( ibits( basis(j), bsites(1, s), sites ) ) !Parity of site1
                    parity2 = popcnt( ibits( ibclr(basis(j), bsites(1, s) - 1), bsites(2, s), sites ) ) !Parity of site2
                    ham_temp(pos,1) = (-1)**(parity1 + parity2) 
                    ham_temp(pos,2) = j
                    ham_temp(pos,3) = loc
                end if 
            else if (btest(basis(j), bsites(2, s) - 1) .and. .not. btest(basis(j), bsites(1, s) - 1)) then 
                newst = ibclr( ibset( basis(j), bsites(1, s) - 1 ), bsites(2, s) - 1 ) !Create on site 2, annihilate on site 1 
                call findstate(dim, newst, basis, loc) 
                if(loc > 0) then 
                    cntr(id+1) = cntr(id+1) + 1
                    cntrj = cntrj + 1
                    pos = (j-1)*nbonds + cntrj 
                    parity1 = popcnt( ibits( basis(j), bsites(2, s), sites ) ) !Parity of site2
                    parity2 = popcnt( ibits( ibclr(basis(j), bsites(2,s) - 1), bsites(1, s), sites ) ) !Parity of site1
                    ham_temp(pos,1) = (-1)**(parity1 + parity2)  
                    ham_temp(pos,2) = j
                    ham_temp(pos,3) = loc
                end if                    
            end if 
        
        end do

    end do
    !$omp end parallel do

    nnz = sum(cntr)

    if(allocated(ham)) deallocate(ham)
    allocate(ham(nnz,3))
    ham = 0 
    cntr2 = 0 
    do i = 1, arrsize
        if(ham_temp(i, 1) .ne. 0 ) then 
            cntr2 = cntr2 + 1 
            ham(cntr2,1) = ham_temp(i,1)
            ham(cntr2,2) = ham_temp(i,2)
            ham(cntr2,3) = ham_temp(i,3)    
        end if 
    end do
    if(cntr2 .ne. nnz) stop 'SR hopping: Counter > NNZ'
    if (allocated(ham_temp)) deallocate(ham_temp)

    print*,'Generated hopping Hamiltonian.'
    print*, ''

end subroutine hopping_p2

subroutine khopping_c(unit, parameters, gencluster, nHel, tilt, lx, ly, ucx, ucy, sites, nbonds, dim, basis, bsites, norm, xtransl, ytransl, t, k1, k2, rc, ham, nnz)

    implicit none

    integer, intent(in) :: unit, gencluster, nHel, tilt, lx, ly, ucx, ucy, sites, nbonds
    integer, intent(in) :: k1
    integer, intent(in) :: k2
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(2, nbonds)
    double precision, intent(in) :: norm(dim)
    integer, intent(in) :: xtransl(2, sites)
    integer, intent(in) :: ytransl(2, sites)
    double precision, intent(in) :: t
    character(len=*), intent(in) :: parameters
    integer, intent(out) :: nnz
    integer, allocatable, intent(out) :: rc(:,:)
    double complex, allocatable, intent(out) :: ham(:)
    
    integer :: dbl = 0 
    integer :: parity1 = 0, parity2 = 0 
    integer :: l1 = 0, l2 = 0, sign = 0, l11 = 0, l22 = 0   
    integer(kind=8) :: loc = 0, rep = 0
    integer(kind=8) :: i = 0, j = 0, k = 0, l = 0, s = 0
    integer(kind=8) :: mask = 0, newst = 0 
    integer(kind=8) :: n_temp = 0, cntr = 0, cntrj = 0
    integer, allocatable :: rc_temp(:,:)
    double complex, allocatable :: ham_temp(:)
    
    if (allocated(ham_temp)) deallocate(ham_temp)
    if (allocated(rc_temp))  deallocate(rc_temp)
    allocate(ham_temp(6*sites*dim))
    allocate(rc_temp(6*sites*dim, 2))
    ham_temp = 0.d0  
    rc_temp = 0 
    n_temp  = 0 
    nnz     = 0 

    do j = 1, dim
        cntrj = n_temp      !Starting point in temp_rc for matrix elements of row 'j' to look for double entries with the same column 'loc'
        cntr  = 0           !Counts the number of already calculated matrix elements of each row 'j'
        do s = 1, nbonds 
            ! !Create scattered basis state I'
                ! mask = ibset(ibset(0,bsites(1,s)-1),bsites(2,s)-1) !Procedure suggested in Lin-paper. Sets a mask with only non-zero components on sites 1 and 2
                ! k    = iand(mask, basis(j)) !K records the occupancy of those two sites.
                ! l    = ieor(mask,k) !L records whether hopping is allowed or not. If it's 0 or the same as the mask, hopping is not allowed. If it is allowed the occupations of both sites (01 or 10) are swapped (10 or 01)
                ! if (l == 0 .or. l == mask) then
                !     cycle
                ! end if
                ! call representative(basis(j) - k + l, sites, ucx, ucy, xtransl, ytransl, rep, l1, l2, sign) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            


            if(btest(basis(j), bsites(1,s) - 1) .and. .not. (btest(basis(j), bsites(2,s) - 1)) ) then 
                newst = ibclr( ibset(basis(j), bsites(2, s) - 1), bsites(1, s) - 1) 
            else if(btest(basis(j), bsites(2,s) - 1) .and. .not. (btest(basis(j), bsites(1,s) - 1)) ) then 
                newst = ibclr( ibset(basis(j), bsites(1, s) - 1), bsites(2, s) - 1) 
            else 
                cycle 
            end if 
            if(gencluster == 1) then 
                call representative(newst, sites, ucx, ucy, xtransl, ytransl, rep, l1, l2, sign) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                l11 = ucx 
                l22 = ucy
            else if(gencluster == 0) then 
                call representative_tilted(newst, sites, nHel, tilt, lx, ly, xtransl, ytransl, rep, l1, l2, sign)
                l11 = ly 
                l22 = lx
            end if 
            call findstate(dim, rep, basis, loc) !Finds the location of representative in basis
            if (loc > 0) then
                if(btest(basis(j),bsites(1,s)-1) .and. .not.(btest(basis(j),bsites(2,s)-1))) then !Hopping from 1->2
                    parity1 = popcnt( ibits( basis(j), bsites(1,s), sites ) )
                    parity2 = popcnt( ibits( ibclr(basis(j), bsites(1,s) - 1), bsites(2,s), sites ) ) 
                else if(btest(basis(j),bsites(2,s)-1) .and. .not.(btest(basis(j),bsites(1,s)-1))) then !Hopping from 2->1
                    parity1 = popcnt( ibits( ibclr(basis(j), bsites(2,s) - 1), bsites(1,s), sites  ) )
                    parity2 = popcnt( ibits( basis(j), bsites(2,s), sites ) ) 
                end if 
                if(loc == j) then 
                    print*,'complex : j, loc',j, loc
                    pause 
                    print*,'--------------------------------------------------------------------------------------------------------------------------------------------------------------------------' 
                end if
                dbl = 0 !Flag for whether a new matrix element should be created (0) or an old matrix element was updated (1)
                cntr = cntr + 1
                if (cntr > 1) then
                    do i = max(cntrj,1), cntrj + cntr !Loop for updating existing elements: Goes through all matrix elements already calculated for row 'j' and checks whether column 'loc' already exists.
                        if(rc_temp(i,2) == loc .and. rc_temp(i,1) == j) then !If yes, update existing element.                            
                            ham_temp(i) = ham_temp(i) + sign * (-1) * t * (-1)**(parity1 + parity2) * sqrt(real(norm(loc))/real(norm(j))) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                            ! ham_temp(i) = ham_temp(i) + (-1) * t * (-1)**(parity1 + parity2) * exp(-ii*2*pi*(k1*l1 + k2*l2)/sites)        
                            ! ham_temp(i) = ham_temp(i) + (-1) * t * (-1)**(parity1 + parity2) * sqrt(real(norm(loc))/real(norm(j))) * exp(-ii*2*pi*(k1*dble(l1-ucx)/2.d0 + k2*dble(l2-l22)/2.d0)/sites)        
                            dbl = 1 !Found and updated existing element.
                        end if
                    end do
                end if
                if (dbl == 0) then !If no existing element was found, create new matrix element.
                    n_temp = n_temp + 1 !Counter for non-zero matrix elements.
                    ham_temp(n_temp) = ham_temp(n_temp) + sign * (-1) * t * (-1)**(parity1 + parity2) * sqrt(real(norm(loc))/real(norm(j))) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                    ! ham_temp(n_temp) = ham_temp(n_temp) + (-1) * t * (-1)**(parity1 + parity2) * exp(-ii*2*pi*(k1*l1 + k2*l2)/sites)
                    ! ham_temp(n_temp) = ham_temp(n_temp) + (-1) * t * (-1)**(parity1 + parity2) * sqrt(real(norm(loc))/real(norm(j))) * exp(-ii*2*pi*(k1*dble(l1-Lx)/2.d0 + k2*dble(l2-Ly)/2.d0)/sites)
                    rc_temp(n_temp,1) = j !Row of matrix element
                    rc_temp(n_temp,2) = loc !Column of matrix element
                end if
            end if
        end do 
    end do 
          
    nnz = n_temp 
    if(allocated(ham)) deallocate(ham)
    if(allocated(rc)) deallocate(rc)
    allocate(ham(nnz))
    allocate(rc(nnz,2))
    ham = 0.d0  
    rc = 0

    do i = 1, nnz
        ham(i) = ham_temp(i)
        rc(i,1) = rc_temp(i,1)
        rc(i,2) = rc_temp(i,2)
    end do

    if (allocated(ham_temp)) deallocate(ham_temp)
    if (allocated(rc_temp)) deallocate(rc_temp)

    print*, 'Finished hopping Hamiltonian'
    print*, ''

end subroutine khopping_c

subroutine hopping_symm(threads, tilted, nHel, tilt, lx, ly, ucx, ucy, sites, nbonds, dim, basis, bsites, norm, xtransl, ytransl, symmetrize, id, par, rot, refl, c6, t, rc, rcdi, prts, dplcts, ham, hamDi, nnz, nDi)

    implicit none

    integer, intent(in) :: threads, tilted, nHel, tilt, lx, ly, ucx, ucy, sites, nbonds, symmetrize
    integer(kind=8), intent(in) :: dim, basis(dim)
    integer, intent(in) :: bsites(2, nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    double precision, intent(in) :: t, id, par(6), rot(5), norm(dim)
    integer, intent(out) :: nnz, nDi 
    integer, allocatable, intent(out) :: rc(:,:), rcdi(:), prts(:), dplcts(:)
    double precision, allocatable, intent(out) :: ham(:), hamDi(:)
    
    integer :: sign = 0, dbl = 0
    integer :: s = 0, l1 = 0, l2 = 0, parity1 = 0, parity2 = 0 
    integer(kind=8) :: i = 0, j = 0, loc = 0, rep = 0, pos = 0
    integer(kind=8) :: newst = 0, cntr = 0, cntrj = 0, n_temp = 0, arrsize = 0
    ! integer(kind=8), allocatable :: n_temp(:)
    integer, allocatable :: rc_temp(:,:), rcdi_temp(:), prts_temp(:), dplcts_temp(:), cntr_di(:,:)
    double precision, allocatable :: ham_temp(:), hamDi_temp(:)
    
    arrsize = 2*dim*nbonds
    if (allocated(ham_temp)) deallocate(ham_temp)
    if (allocated(rc_temp)) deallocate(rc_temp)
    if (allocated(hamDi_temp)) deallocate(hamDi_temp)
    if (allocated(rcdi_temp)) deallocate(rcdi_temp)
    if (allocated(prts_temp)) deallocate(prts_temp)
    if (allocated(dplcts_temp)) deallocate(dplcts_temp)
    if (allocated(cntr_di)) deallocate(cntr_di)
    
    ! if (allocated(n_temp)) deallocate(n_temp)
    allocate(ham_temp(arrsize))
    allocate(rc_temp(arrsize, 2))
    allocate(hamDi_temp(arrsize))
    allocate(rcdi_temp(arrsize))
    allocate(prts_temp(arrsize))
    allocate(dplcts_temp(arrsize))
    allocate(cntr_di(dim, 2))
    ! allocate(n_temp(threads))
    ham_temp        = 0.d0  
    rc_temp         = 0 
    hamDi_temp     = 0.d0  
    rcdi_temp      = 0 
    prts_temp   = 0 
    dplcts_temp = 0 
    n_temp          = 0 
    nDi             = 0 
    nnz             = 0 
    cntr_di         = 0 

    !!$omp parallel do default(firstprivate) shared(ham_temp, bsites, basis) num_threads(threads)
    do j = 1, dim
        cntrj = 0 !Counts the number of already calculated matrix elements of each row 'j'
        !!$ id = omp_get_thread_num()
        do s = 1, nbonds 
            if(btest(basis(j), bsites(1,s) - 1) .and. .not.(btest(basis(j), bsites(2,s) - 1))) then 
                newst = ibclr( ibset(basis(j), bsites(2, s) - 1), bsites(1, s) - 1)                 
            else if(btest(basis(j), bsites(2,s) - 1) .and. .not.(btest(basis(j), bsites(1,s) - 1))) then 
                newst = ibclr( ibset(basis(j), bsites(1, s) - 1), bsites(2, s) - 1) 
            else 
                cycle 
            end if 
            if(tilted == 0) then 
                ! call representative(newst, sites, ucx, ucy, xtransl, ytransl, rep, l1, l2, sign) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                call representative_irrep_rect(newst, sites, Lx, Ly, symmetrize, id, par, rot, xtransl, ytransl, refl, c6, rep, l1, l2, sign)
            else if(tilted == 1) then 
                call representative_irrep(newst, sites, nHel, tilt, Lx, Ly, symmetrize, id, par, rot, xtransl, ytransl, refl, c6, rep, l1, l2, sign)
                ! call representative_irrep2(newst, sites, nHel, tilt, Lx, Ly, symmetrize, id, par, rot, xtransl, ytransl, refl, c6, rep, l1, l2, sign)
                ! call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, sign)
            end if 
            
            call findstate(dim, rep, basis, loc) !Finds the location of representative in basis
            if (loc <= 0) cycle 

            if(btest(basis(j),bsites(1,s)-1) .and. .not.(btest(basis(j),bsites(2,s)-1))) then !Hopping from 1->2
                parity1 = popcnt(ibits(basis(j), bsites(1,s), sites))
                parity2 = popcnt(ibits(ibclr(basis(j), bsites(1,s) - 1), bsites(2,s), sites)) 
            else if(btest(basis(j),bsites(2,s)-1) .and. .not.(btest(basis(j),bsites(1,s)-1))) then !Hopping from 2->1
                parity2 = popcnt(ibits(basis(j), bsites(2,s), sites)) 
                parity1 = popcnt(ibits(ibclr(basis(j), bsites(2,s) - 1), bsites(1,s), sites))
            end if 
            if(loc == j) then
                if(cntr_di(j, 1) == 1) then 
                    hamDi_temp(cntr_di(j, 2)) = hamDi_temp(cntr_di(j, 2)) + sign * (-1) * t * (-1)**(parity1 + parity2) * sqrt(dble(norm(loc))/dble(norm(j))) 
                else 
                    nDi = nDi + 1
                    hamDi_temp(nDi) = hamDi_temp(nDi) + sign * (-1) * t * (-1)**(parity1 + parity2) * sqrt(dble(norm(loc))/dble(norm(j))) 
                    rcdi_temp(nDi) = j 
                    cntr_di(j, 1) = 1
                    cntr_di(j, 2) = nDi
                end if 
            else    
                dbl  = 0 !Flag for whether a new matrix element should be created (0) or an old matrix element was updated (1)
                ! cntr = cntr + 1
                if (cntrj > 0) then 
                    do i =  1, (j)*nbonds + cntrj + 1 !Q 
                    ! do i =  1, (j-1)*nbonds + cntrj + 1 !Q 
                    ! do i = (j-1)*nbonds + 1, (j-1)*nbonds + cntrj + 1 !Q 
                        if(rc_temp(i,2) == loc .and. rc_temp(i,1) == j) then !If yes, update existing element. 
                            if(i > arrsize) stop 'hopping_symm: i > arrsize'
                            ham_temp(i) = ham_temp(i) + sign * (-1) * t * (-1)**(parity1 + parity2) * sqrt(dble(norm(loc))/dble(norm(j))) 
                            prts_temp(i) = (-1)**(parity1 + parity2)
                            dbl = 1 !Found and updated existing element.
                            dplcts_temp(i) = dplcts_temp(i) + 1 
                        end if
                    end do
                end if
                if (dbl == 0) then !If no existing element was found, create new matrix element.
                    cntrj = cntrj + 1
                    pos = (j-1) * nbonds + cntrj 
                    if(pos > arrsize) stop 'hopping_symm: pos > arrsize'
                    n_temp = n_temp + 1 !Counter for non-zero matrix elements.         
                    ! n_temp(id) = n_temp(id) + 1 !Counter for non-zero matrix elements.       
                    ham_temp(pos) = ham_temp(pos) + sign * (-1) * t * (-1)**(parity1 + parity2) * sqrt(dble(norm(loc))/dble(norm(j))) 
                    prts_temp(pos) = (-1)**(parity1 + parity2)
                    dplcts_temp(pos) = dplcts_temp(pos) + 1
                    rc_temp(pos,1) = j !Row of matrix element
                    rc_temp(pos,2) = loc !Column of matrix element                          
                end if
            end if 
        end do 
    end do 

    nnz = n_temp
    ! nnz = sum(n_temp)

    if(allocated(rc)) deallocate(rc)
    if(allocated(ham)) deallocate(ham)
    if(allocated(rcdi)) deallocate(rcdi)
    if(allocated(hamDi)) deallocate(hamDi)
    if(allocated(prts)) deallocate(prts)
    if(allocated(dplcts)) deallocate(dplcts)

    allocate(ham(nnz))
    allocate(rc(nnz,2))
    allocate(rcdi(nDi))
    allocate(hamDi(nDi))
    allocate(prts(nnz))
    allocate(dplcts(nnz))

    rc = 0
    cntr = 0
    rcdi = 0
    ham = 0.d0  
    hamDi = 0.d0  

    do i = 1, arrsize!nnz
        if(rc_temp(i,1) > 0) then 
        ! if(ham_temp(i) .ne. 0.d0) then 
            cntr = cntr + 1
            ham(cntr) = ham_temp(i)
            rc(cntr,1) = rc_temp(i,1)
            rc(cntr,2) = rc_temp(i,2)
            prts(cntr) = prts_temp(i)
            dplcts(cntr) = dplcts_temp(i)
        end if 
    end do
    

    ! if(cntr .ne. nnz) then !If matrix elements are zero due to dplcts 
    !     ham_temp = 0.d0 
    !     rc_temp = 0 
    !     ham_temp(1:cntr) = ham(1:cntr)
    !     rc_temp(1:cntr, 1) = rc(1:cntr, 1)
    !     rc_temp(1:cntr, 2) = rc(1:cntr, 2)
    !     if(allocated(ham)) deallocate(ham)
    !     if(allocated(rc)) deallocate(rc)
    !     allocate(ham(cntr))
    !     allocate(rc(cntr, 2))
    !     ham = ham_temp(1:cntr) 
    !     rc(1:cntr, 1) = rc_temp(1:cntr, 1) 
    !     rc(1:cntr, 2) = rc_temp(1:cntr, 2) 
    !     nnz = cntr 
    ! end if 

    do i = 1, nDi
        hamDi(i) = hamDi_temp(i)
        rcdi(i)  = rcdi_temp(i)
    end do

    ! if (allocated(ham_temp)) deallocate(ham_temp)
    ! if (allocated(rc_temp)) deallocate(rc_temp)
    ! if (allocated(hamDi_temp)) deallocate(hamDi_temp)
    ! if (allocated(rcdi_temp)) deallocate(rcdi_temp)
    ! if (allocated(prts_temp)) deallocate(prts_temp)
    ! if (allocated(dplcts_temp)) deallocate(dplcts_temp)

    print*, 'Finished hopping Hamiltonian'
    print*, ''

end subroutine hopping_symm

subroutine hopping_irrep(threads, tilted, nHel, tilt, lx, ly, ucx, ucy, sites, nbonds, dim, basis, bsites, norm, xtransl, ytransl, symmetrize, id, par, rot, refl, c6, t, rc, rcdi, prts, dplcts, ham, hamDi, nnz, nDi)

    implicit none

    integer(kind=8), intent(in) :: dim, basis(dim)
    integer, intent(in) :: threads, tilted, nHel, tilt, lx, ly, ucx, ucy, sites, nbonds, symmetrize
    integer, intent(in) :: bsites(2, nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    double precision, intent(in) :: t, id, par(6), rot(5), norm(dim)
    integer, intent(out) :: nnz, nDi 
    integer, allocatable, intent(out) :: rc(:,:), rcdi(:), prts(:), dplcts(:)
    double precision, allocatable, intent(out) :: ham(:), hamDi(:)
    
    integer :: sign = 1, signrep = 1, signstate = 1, dbl = 0, order, info
    integer :: s = 0, k = 0, l1 = 0, l2 = 0, parity = 0 
    integer(kind=8) :: i = 0, j = 0, loc = 0, rep = 0, pos = 0, rowst = 0
    integer(kind=8) :: newst = 0, state = 0, cntr = 0, cntrj = 0, n_temp = 0, arrsize = 0
    ! integer(kind=8), allocatable :: n_temp(:)
    integer, allocatable :: rc_temp(:,:), rcdi_temp(:), prts_temp(:), dplcts_temp(:), cntr_di(:,:)
    double precision :: h_add = 0.d0 
    double precision, allocatable :: ham_temp(:), hamDi_temp(:)
    
    if(symmetrize == 1) then 
        order = 12
    else 
        order = 1
    end if 

    arrsize = dim * nbonds * order 
    if (allocated(ham_temp)) deallocate(ham_temp)
    if (allocated(rc_temp)) deallocate(rc_temp)
    if (allocated(hamDi_temp)) deallocate(hamDi_temp)
    if (allocated(rcdi_temp)) deallocate(rcdi_temp)
    if (allocated(prts_temp)) deallocate(prts_temp)
    if (allocated(dplcts_temp)) deallocate(dplcts_temp)
    if (allocated(cntr_di)) deallocate(cntr_di)
    
    ! if (allocated(n_temp)) deallocate(n_temp)
    allocate(ham_temp(arrsize))
    allocate(rc_temp(arrsize, 2))
    allocate(hamDi_temp(arrsize))
    allocate(rcdi_temp(arrsize))
    allocate(prts_temp(arrsize))
    allocate(dplcts_temp(arrsize))
    allocate(cntr_di(dim, 2))
    ! allocate(n_temp(threads))
    h_add           = 0.d0  
    ham_temp        = 0.d0  
    hamDi_temp     = 0.d0  
    rc_temp         = 0 
    rcdi_temp      = 0 
    prts_temp   = 0 
    dplcts_temp = 0 
    n_temp          = 0 
    nDi             = 0 
    nnz             = 0 
    cntr_di         = 0 

    !!$omp parallel do default(firstprivate) shared(ham_temp, bsites, basis) num_threads(threads)
    do j = 1, dim
        cntrj = 0 !Counts the number of already calculated matrix elements of each row 'j'
        rowst = (j-1) * nbonds * order 
        !!$ id = omp_get_thread_num()
        do k = 1, order !Loop over point groups operators 
            h_add     = 0.d0 !Scattering contribution to matrix element
            signstate = 1 !Sign due to symmetry transformations and translations
            if(k == 1) then !Prepare initial state 
                state = basis(j)
            else if(k <= 7) then 
                call reflect(basis(j), basis(j), sites, refl(k - 1, 1:sites), signstate, info, state)
                signstate = signstate * par(k - 1)
            else if(8 <= k) then 
                call c6n(basis(j), basis(j), sites, k - 7, c6, signstate, info, state)
                signstate = signstate * rot(k - 7)
            end if 

            do s = 1, nbonds !Takes care of translations
                if(btest(state, bsites(1,s)-1) .and. .not.(btest(state, bsites(2,s)-1))) then !Hopping from 1->2
                    newst  = ibclr(ibset(state, bsites(2, s) - 1), bsites(1, s) - 1)            
                    parity = popcnt(ibits(state, bsites(1,s), sites)) + popcnt(ibits(ibclr(state, bsites(1,s) - 1), bsites(2,s), sites))      
                else if(btest(state, bsites(2,s)-1) .and. .not.(btest(state, bsites(1,s)-1))) then !Hopping from 2->1
                    newst  = ibclr(ibset(state, bsites(1, s) - 1), bsites(2, s) - 1) 
                    parity = popcnt(ibits(state, bsites(2,s), sites)) + popcnt(ibits(ibclr(state, bsites(2,s) - 1), bsites(1,s), sites))
                else 
                    cycle 
                end if 
                if(tilted == 0) then 
                    call representative_irrep_rect(newst, sites, Lx, Ly, symmetrize, id, par, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                else if(tilted == 1) then 
                    call representative_irrep(newst, sites, nHel, tilt, Lx, Ly, symmetrize, id, par, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                end if 
                
                call findstate(dim, rep, basis, loc) !Finds the location of representative in basis
                if(loc <= 0) cycle 
                sign  = signrep * signstate * (-1) * (-1)**parity 
                h_add = sign * t * sqrt(dble(norm(loc))/dble(norm(j)))  

                if(loc == j) then !Diagonal element contributes to diagonal Hamiltonian
                    if(cntr_di(j, 1) == 1) then !Not the first diagonal element of state j. Add to previous element inDicated by cntr_di(j, 2). 
                        hamDi_temp(cntr_di(j, 2)) = hamDi_temp(cntr_di(j, 2)) + h_add
                    else !First diagonal element of state j. Create new element and save position in cntr_di(j, 2).
                        nDi              = nDi + 1
                        hamDi_temp(nDi) = hamDi_temp(nDi) + h_add
                        rcdi_temp(nDi)  = j 
                        cntr_di(j, 1)    = 1
                        cntr_di(j, 2)    = nDi
                    end if 
                else !Off-diagonal elements 
                    dbl = 0 !Flag for whether a new matrix element should be created (0) or an old matrix element was updated (1)
                    if(cntrj > 0) then !Not the first element in row j. Check elements for double entries. 
                        do i = rowst + 1, rowst + cntrj + 1 !Go from beginning of row j to current element. 
                            if(rc_temp(i,2) == loc .and. rc_temp(i,1) == j) then !If yes, update existing element.
                                if(i > arrsize) stop 'hopping irrep: i > arrsize'
                                ham_temp(i)        = ham_temp(i) + h_add 
                                prts_temp(i)   = (-1)**parity
                                dplcts_temp(i) = dplcts_temp(i) + 1 
                                dbl                = 1 !Found and updated existing element.       
                            end if
                        end do
                    end if
                    if (dbl == 0) then !If no existing element was found, create new matrix element.
                        if(pos > arrsize) stop 'hopping irrep: pos > arrsize'
                        ! n_temp(id) = n_temp(id) + 1 !Counter for non-zero matrix elements.       
                        cntrj                = cntrj + 1 !One more element found for row j. 
                        n_temp               = n_temp + 1 !Counter for non-zero matrix elements. Needed later.
                        pos                  = rowst + cntrj !Position corresponds to position in row j. 
                        ham_temp(pos)        = ham_temp(pos) + h_add !Value of matrix element. 
                        prts_temp(pos)   = (-1)**parity !Save prts for later. 
                        dplcts_temp(pos) = dplcts_temp(pos) + 1 !Save number of dplcts for later. 
                        rc_temp(pos,1)       = j !Row of matrix element. 
                        rc_temp(pos,2)       = loc !Column of matrix element.                          
                    end if
                end if 
            end do 

        end do !point group operations
    end do 

    nnz = n_temp

    if(allocated(rc)) deallocate(rc)
    if(allocated(ham)) deallocate(ham)
    if(allocated(rcdi)) deallocate(rcdi)
    if(allocated(hamDi)) deallocate(hamDi)
    if(allocated(prts)) deallocate(prts)
    if(allocated(dplcts)) deallocate(dplcts)

    allocate(ham(nnz))
    allocate(rc(nnz,2))
    allocate(rcdi(nDi))
    allocate(hamDi(nDi))
    allocate(prts(nnz))
    allocate(dplcts(nnz))

    rc     = 0
    cntr   = 0
    rcdi  = 0
    ham    = 0.d0  
    hamDi = 0.d0  

    do i = 1, arrsize!nnz
        if(rc_temp(i,1) > 0) then 
        ! if(ham_temp(i) .ne. 0.d0) then 
            cntr       = cntr + 1
            ham(cntr)  = ham_temp(i)
            rc(cntr,1) = rc_temp(i,1)
            rc(cntr,2) = rc_temp(i,2)
            prts(cntr) = prts_temp(i)
            dplcts(cntr) = dplcts_temp(i)
        end if 
    end do
    ham = ham / dble(order)

    ! if(cntr .ne. nnz) then !If matrix elements are zero due to dplcts 
    !     ham_temp = 0.d0 
    !     rc_temp = 0 
    !     ham_temp(1:cntr) = ham(1:cntr)
    !     rc_temp(1:cntr, 1) = rc(1:cntr, 1)
    !     rc_temp(1:cntr, 2) = rc(1:cntr, 2)
    !     if(allocated(ham)) deallocate(ham)
    !     if(allocated(rc)) deallocate(rc)
    !     allocate(ham(cntr))
    !     allocate(rc(cntr, 2))
    !     ham = ham_temp(1:cntr) 
    !     rc(1:cntr, 1) = rc_temp(1:cntr, 1) 
    !     rc(1:cntr, 2) = rc_temp(1:cntr, 2) 
    !     nnz = cntr 
    ! end if 

    do i = 1, nDi !Fill up diagonal Hamiltonian
        hamDi(i) = hamDi_temp(i)
        rcdi(i)  = rcdi_temp(i)
    end do
    hamDi = hamDi / order 

    ! if (allocated(ham_temp)) deallocate(ham_temp)
    ! if (allocated(rc_temp)) deallocate(rc_temp)
    ! if (allocated(hamDi_temp)) deallocate(hamDi_temp)
    ! if (allocated(rcdi_temp)) deallocate(rcdi_temp)
    ! if (allocated(prts_temp)) deallocate(prts_temp)
    ! if (allocated(dplcts_temp)) deallocate(dplcts_temp)

    print*, 'Finished hopping Hamiltonian'
    print*, ''

end subroutine hopping_irrep

subroutine hopping_irrep2D(threads, tilted, nHel, tilt, lx, ly, sites, nbonds, dim, basis, orbsize, orbits, phases, norm, bsites, xtransl, ytransl, symmetrize, id, par, rot, refl, c6, t, rc, rcdi, prts, dplcts, ham, hamDi, nnz, nDi)
    
    implicit none
    integer, intent(in):: orbsize, threads, tilted, nHel, tilt, lx, ly, sites, nbonds, symmetrize
    integer(kind=8), intent(in) :: dim, basis(dim)
    integer, intent(in) :: bsites(2, nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    integer(kind=8), intent(in) :: orbits(:,:,:)
    double precision, intent(in) :: t, id, par(6), rot(5), norm(dim,2)
    double complex, allocatable, intent(in) :: phases(:,:,:)
    integer, intent(out) :: nnz, nDi 
    integer, allocatable, intent(out) :: rc(:,:), rcdi(:), prts(:), dplcts(:)
    double complex, allocatable, intent(out) :: ham(:), hamDi(:)
    
    integer :: sign = 1, signrep, dbl, order, site1, site2
    integer :: s, c, d, o, l1, l2, parity 
    integer(kind=8) :: i, j, loc, rep, pos, rowst, rowIndx, colIndx
    integer(kind=8) :: newst, state, cntr, cntrjc, n_temp, arrsize
    ! integer(kind=8), allocatable :: n_temp(:)
    integer, allocatable :: rc_temp(:,:), rcdi_temp(:), prts_temp(:), dplcts_temp(:), cntr_di(:,:)
    double precision, parameter :: tol = 1.0e-14
    double complex, allocatable :: ham_temp(:), hamDi_temp(:)
    double complex :: h_add, coeff, newcf 
    
    if(symmetrize == 1) then 
        order = 12
    else 
        order = 1
    end if 

    arrsize = dim * orbsize * 2 !Each representative (dim) has two basis states (2), each of which is a linear combination of at most 'orbsize' states. 
    if(allocated(ham_temp))    deallocate(ham_temp)
    if(allocated(rc_temp))     deallocate(rc_temp)
    if(allocated(hamDi_temp))  deallocate(hamDi_temp)
    if(allocated(rcdi_temp))   deallocate(rcdi_temp)
    if(allocated(prts_temp))   deallocate(prts_temp)
    if(allocated(dplcts_temp)) deallocate(dplcts_temp)
    if(allocated(cntr_di))     deallocate(cntr_di)
    
    ! if (allocated(n_temp)) deallocate(n_temp)
    allocate(dplcts_temp(arrsize))
    allocate(hamDi_temp(arrsize))
    allocate(rc_temp(arrsize, 2))
    allocate(rcdi_temp(arrsize))
    allocate(prts_temp(arrsize))
    allocate(ham_temp(arrsize)) 
    allocate(cntr_di(dim, 2))
    ! allocate(n_temp(threads))
    h_add       = 0.d0
    ham_temp    = 0.d0  
    hamDi_temp  = 0.d0  
    rc_temp     = 0 
    rcdi_temp   = 0 
    prts_temp   = 0 
    dplcts_temp = 0 
    n_temp      = 0 
    nDi         = 0 
    nnz         = 0 
    cntr_di     = 0 
    site1       = bsites(1,1)
    site2       = bsites(2,1)
    !!$omp parallel do default(firstprivate) shared(ham_temp, bsites, basis) num_threads(threads)
    
    do j = 1, dim !All representatives 
        !!$ id = omp_get_thread_num()
        
        do c = 1, 2 !First and second basis state of 2D irrep 
        
            h_add = 0.d0 !Scattering contribution to matrix element
            cntrjc = 0 !Counts the number of already calculated matrix elements of each combination row 'j' and basist state 'c', i.e. tuple (j, c) 
            rowIndx = 2*(j-1) + c
            rowst = (rowIndx-1) * orbsize !Position-1 of beginning of tuple (j,c)   
            do s = 1, orbsize !Loop through orbit of state basis(j)
                state = orbits(j, s, c)
                coeff = phases(j, s, c)
                if(state == 0) exit 
                if(abs(dble(coeff)) <= tol .and. abs(aimag(coeff)) <= tol) cycle 
                if(btest(state, site1-1) .and. .not.(btest(state, site2-1))) then !Hopping from 1->2
                    newst  = ibclr(ibset(state, site2-1), site1-1)            
                    parity = popcnt(ibits(state, site1, sites)) + popcnt(ibits(ibclr(state, site1 - 1), site2, sites))      
                else if(btest(state, site2-1) .and. .not.(btest(state, site1-1))) then !Hopping from 2->1
                    newst  = ibclr(ibset(state, site1-1), site2-1) 
                    parity = popcnt(ibits(state, site2, sites)) + popcnt(ibits(ibclr(state, site2-1), site1, sites))
                else 
                    cycle 
                end if 
                if(tilted == 0) then 
                    call representative_irrep_rect(newst, sites, Lx, Ly, symmetrize, id, par, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                else if(tilted == 1) then 
                    call representative_irrep(newst, sites, nHel, tilt, Lx, Ly, symmetrize, id, par, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                end if 
                call findstate(dim, rep, basis, loc) !Finds the location of representative in basis
                if(loc <= 0) cycle !New state is compatible with momentum and symmetry 
                do d = 1, 2 !Irrep basis states of scattered representative
                    colIndx = 2*(loc-1) + d
                    ! if(colIndx == 992 .and. rowIndx == 987) print* ,phases(j, 1:orbsize, c), 'phases(j, 1:orbsize, c)'
                    ! if(colIndx == 992 .and. rowIndx == 987) print* ,phases(j, s, c), 'phases(j, s, c)'
                    newcf = 0.d0 
                    do o = 1, orbsize !Search for position of scattered state in orbit of Irrep's basis state and assign corresponDing coefficient 'newcf'   
                        if(orbits(loc, o, d) == newst) then 
                            newcf = phases(loc, o, d)
                            ! if(orbits(loc, o, d) == 91750 .or. orbits(loc, o, d) == 91749 .or. orbits(loc, o, d) == 109157 .or. orbits(loc, o, d) == 109158 .or. orbits(loc, o, d) == 173417 .or. orbits(loc, o, d) == 173418) then 
                            !     print* ,int(state), 'state'
                            !     print* ,int(s),int(c), 's, c'
                            !     print* ,int(o),int(d), 'o, d '
                            !     print* ,int(orbits(loc, o, d)), 'orbits(loc, o, d)'
                            !     ! print* ,phases(loc, o, d), 'phases(loc, o, d)'
                            !     print*, ''
                            ! end if
                            exit 
                        end if 
                        ! if(orbits(loc, o, d) == newst .and. newst == 39326) then 
                        !     print* ,j, 'j'
                        !     print* ,loc, 'loc'
                        !     print* ,c, 'c'
                        !     print* ,d, 'd'
                        !     print* ,s, 's'
                        !     print* ,o, 'o'
                        !     print* ,rowindx, 'rowindx'
                        !     print* ,colindx, 'colindx'
                        !     print* ,phases(j, s, c), 'phases(j, s, c)'
                        !     print* ,phases(loc, o, d), 'phases(loc, o, d)'
                        !     print* ,state, 'state'
                        !     print* ,newst, 'newst'
                        !     print* ,rep, 'rep'
                        !     pause 
                        ! end if 
                    end do 
                    if(abs(dble(newcf)) <= tol .and. abs(aimag(newcf)) <= tol) cycle 

                    sign  = (-1) * (-1)**parity     
                    h_add = sign * t * sqrt(dble(norm(loc, d))/dble(norm(j, c))) * coeff * newcf !Contribution to value of matrix element <j,c| H | loc, c> 
                    ! h_add = sign * t * coeff * newcf !Contribution to value of matrix element <j,c| H | loc, c> 
                    
                    if(colIndx == rowIndx) then !Matrix element contributes to diagonal Hamiltonian
                        if(cntr_di(rowIndx, 1) == 1) then !Not the first diagonal element of state j. Add to previous element inDicated by cntr_di(rowIndx, 2). 
                            hamDi_temp(cntr_di(rowIndx, 2)) = hamDi_temp(cntr_di(rowIndx, 2)) + h_add
                        else !First diagonal element of state j. Create new element and save position in cntr_di(j, 2).
                            nDi                 = nDi + 1 !Number of diagonal matrix elements
                            hamDi_temp(nDi)     = hamDi_temp(nDi) + h_add !Diagonal matrix element
                            rcdi_temp(nDi)      = rowIndx !Row and column index of matrix element 
                            cntr_di(rowIndx, 1) = 1 !Flag for existing diagonal matrix element in row 'rowIndx'
                            cntr_di(rowIndx, 2) = nDi !Position of existing element of 'rowIndx' in array 'hamDi_temp'
                        end if 
                    else !Off-diagonal elements 
                        dbl = 0 !Flag for whether a new matrix element should be created (0) or an old matrix element was updated (1)
                        if(cntrjc > 0) then !Not the first element in row j. Check elements for double entries. 
                            do i = rowst + 1, rowst + cntrjc + 1 !Go from beginning of row j to current element. 
                                if(rc_temp(i,2) == colIndx .and. rc_temp(i,1) == rowIndx) then !If yes, update existing element.
                                    if(i > arrsize) stop 'hopping irrep2D: i > arrsize'
                                    ham_temp(i)    = ham_temp(i) + h_add 
                                    prts_temp(i)   = (-1)**parity
                                    dplcts_temp(i) = dplcts_temp(i) + 1 
                                    dbl            = 1 !Found and updated existing element.       
                                    if(rowIndx == 992 .and. colIndx == 987) then 
                                        print* ,int8(colIndx), 'colIndx'
                                        print* ,int8(rowIndx), 'rowIndx'
                                        ! print* ,int8(i), 'i'
                                        ! print* ,int8(j), 'j'
                                        ! print* ,int8(loc), 'loc'
                                        print* , int8(c),  ' c'
                                        print* , int8(d),  ' d'
                                        print* , int8(s),  ' s'
                                        print* , int8(o),  ' o'
                                        ! print* ,phases(loc, o, d), 'phases(loc, o, d)'
                                        print* ,int8(sign),  'sign'
                                        print* ,int8(prts_temp(pos)), 'parity'
                                        print* ,state, 'state'
                                        print* ,newst, 'newst'
                                        print* ,rep, 'rep'
                                        print*, norm(j, c), 'norm(j, c)'
                                        print*, norm(loc, d), 'norm(loc, d)'
                                        print* ,sqrt(dble(norm(loc, d))/dble(norm(j, c))), 'sqrt'
                                        print* ,coeff, 'coeff'
                                        print* ,newcf, 'newcf'
                                        print* ,h_add, 'h_add'
                                        print* ,ham_temp(i) , 'ham_temp(i) '
                                        
                                        pause
                                    else if(rowIndx == 987 .and. colIndx == 992) then 
                                        print* ,colIndx, 'colIndx'
                                        print* ,rowIndx, 'rowIndx'
                                        ! print* ,int8(i), 'i'
                                        ! print* ,int8(j), 'j'
                                        ! print* ,int8(loc), 'loc'
                                        print* , int8(c),  ' c'
                                        print* , int8(d),  ' d'
                                        print* , int8(s),  ' s'
                                        print* , int8(o),  ' o'
                                        ! print* ,phases(loc, o, d), 'phases(loc, o, d)'
                                        print* ,int8(sign),  'sign'
                                        print* ,int8(prts_temp(pos)), 'parity'
                                        print* ,state, 'state'
                                        print* ,newst, 'newst'
                                        print* ,rep, 'rep'
                                        print*, norm(j, c), 'norm(j, c)'
                                        print*, norm(loc, d), 'norm(loc, d)'
                                        print* ,sqrt(dble(norm(loc, d))/dble(norm(j, c))), 'sqrt'
                                        print* ,coeff, 'coeff'
                                        print* ,newcf, 'newcf'
                                        print* ,h_add, 'h_add'
                                        print* ,ham_temp(i) , 'ham_temp(i) '
                                        
                                        pause 
                                    end if 
                                end if
                            end do
                        end if
                        if (dbl == 0) then !If no existing element was found, create new matrix element.
                            ! n_temp(id) = n_temp(id) + 1 !Counter for non-zero matrix elements.       
                            cntrjc           = cntrjc + 1 !One more element found for row j, basis state c. 
                            n_temp           = n_temp + 1 !Counter for non-zero matrix elements. Needed later.
                            pos              = rowst + cntrjc !Position corresponds to position in row j. 
                            if(pos > arrsize) stop 'hopping irrep2D: pos > arrsize'
                            ham_temp(pos)    = ham_temp(pos) + h_add !Value of matrix element. 
                            prts_temp(pos)   = (-1)**parity !Save prts for later. 
                            dplcts_temp(pos) = dplcts_temp(pos) + 1 !Save number of dplcts for later. 
                            rc_temp(pos,1)   = rowIndx !Row of matrix element. 
                            rc_temp(pos,2)   = colIndx !Column of matrix element.                       
                            if(rowIndx == 992 .and. colIndx == 987) then 
                                print* ,colIndx, 'colIndx'
                                print* ,rowIndx, 'rowIndx'
                                print* ,int8(i), 'i'
                                print* ,int8(j), 'j'
                                print* ,int8(loc), 'loc'
                                print* , int8(c),  ' c'
                                print* , int8(d),  ' d'
                                print* , int8(s),  ' s'
                                print* , int8(o),  ' o'
                                print* , int8(sign),  'sign'

                                print* ,int8(prts_temp(pos)), 'parity'
                                print* ,state, 'state'
                                print* ,newst, 'newst'
                                print* ,rep, 'rep'
                                print*, norm(j, c), 'norm(j, c)'
                                print*, norm(loc, d), 'norm(loc, d)'
                                print* ,sqrt(dble(norm(loc, d))/dble(norm(j, c))), 'sqrt'
                                print* ,coeff, 'coeff'
                                print* ,newcf, 'newcf'
                                print* ,h_add, 'h_add'
                                print* ,ham_temp(pos) , 'ham_temp(pos) '

                                pause

                            else if(rowIndx == 987 .and. colIndx == 992) then 
                                print* ,colIndx, 'colIndx'
                                print* ,rowIndx, 'rowIndx'
                                print* ,int8(i), 'i'
                                print* ,int8(j), 'j'
                                print* ,int8(loc), 'loc'
                                print* , int8(c),  ' c'
                                print* , int8(d),  ' d'
                                print* , int8(s),  ' s'
                                print* , int8(o),  ' o'
                                print* ,int8(sign),  'sign'
                                ! print* ,phases(loc, o, d), 'phases(loc, o, d)'
                                print* ,int8(prts_temp(pos)), 'parity'
                                print* ,state, 'state'
                                print* ,newst, 'newst'
                                print* ,rep, 'rep'
                                print*, norm(j, c), 'norm(j, c)'
                                print*, norm(loc, d), 'norm(loc, d)'
                                print* ,sqrt(dble(norm(loc, d))/dble(norm(j, c))), 'sqrt'
                                print* ,coeff, 'coeff'
                                print* ,newcf, 'newcf'
                                print* ,h_add, 'h_add'
                                print* ,ham_temp(pos) , 'ham_temp(pos) '

                                pause

                            end if    
                        end if
                    end if 
                end do !Loop over scattered irrep basis states (d) 
            end do !Loop over symmetry orbit (s)
        end do !Loop over irrep basis states (c)
    end do !Loop over representatives (j)

    nnz = n_temp
    

    if(allocated(rc)) deallocate(rc)
    if(allocated(ham)) deallocate(ham)
    if(allocated(rcdi)) deallocate(rcdi)
    if(allocated(hamDi)) deallocate(hamDi)
    if(allocated(prts)) deallocate(prts)
    if(allocated(dplcts)) deallocate(dplcts)

    allocate(ham(nnz))
    allocate(rc(nnz,2))
    allocate(rcdi(nDi))
    allocate(hamDi(nDi))
    allocate(prts(nnz))
    allocate(dplcts(nnz))

    rc    = 0
    cntr  = 0
    rcdi  = 0
    ham   = 0.d0  
    hamDi = 0.d0  

    do i = 1, arrsize!nnz
        if(rc_temp(i,1) > 0) then 
        ! if(ham_temp(i) .ne. 0.d0) then 
            cntr         = cntr + 1
            ham(cntr)    = ham_temp(i)
            rc(cntr,1)   = rc_temp(i,1)
            rc(cntr,2)   = rc_temp(i,2)
            prts(cntr)   = prts_temp(i)
            dplcts(cntr) = dplcts_temp(i)        
        end if 
    end do
    ham = ham / dble(order)

    
    ! if(cntr .ne. nnz) then !If matrix elements are zero due to dplcts 
    !     ham_temp = 0.d0 
    !     rc_temp = 0 
    !     ham_temp(1:cntr) = ham(1:cntr)
    !     rc_temp(1:cntr, 1) = rc(1:cntr, 1)
    !     rc_temp(1:cntr, 2) = rc(1:cntr, 2)
    !     if(allocated(ham)) deallocate(ham)
    !     if(allocated(rc)) deallocate(rc)
    !     allocate(ham(cntr))
    !     allocate(rc(cntr, 2))
    !     ham = ham_temp(1:cntr) 
    !     rc(1:cntr, 1) = rc_temp(1:cntr, 1) 
    !     rc(1:cntr, 2) = rc_temp(1:cntr, 2) 
    !     nnz = cntr 
    ! end if 

    do i = 1, nDi !Fill up diagonal Hamiltonian
        hamDi(i) = hamDi_temp(i)
        rcdi(i)  = rcdi_temp(i)
    end do
    hamDi = hamDi / order 

    ! if (allocated(ham_temp)) deallocate(ham_temp)
    ! if (allocated(rc_temp)) deallocate(rc_temp)
    ! if (allocated(hamDi_temp)) deallocate(hamDi_temp)
    ! if (allocated(rcdi_temp)) deallocate(rcdi_temp)
    ! if (allocated(prts_temp)) deallocate(prts_temp)
    ! if (allocated(dplcts_temp)) deallocate(dplcts_temp)

    print*, 'Finished hopping Hamiltonian'
    print*, ''

end subroutine hopping_irrep2D

subroutine diagonal(unit, parameters, sites, nbonds, nbonds2, dim, bsites, bsites2, basis, &
                    ham, occ, nnz)

    implicit none

    integer(kind=8), intent(in) :: dim, basis(dim)
    integer, intent(in) :: unit, sites, nbonds, nbonds2
    integer, intent(in) :: bsites(2, nbonds), bsites2(2, nbonds2)
    character(len=*), intent(in) :: parameters
    integer, intent(out) :: nnz
    integer, allocatable, intent(out) :: occ(:,:)
    integer, allocatable, intent(out) :: ham(:,:)

    integer :: counter_v = 0, counter_v2 = 0
    integer :: j = 0, m = 0, s = 0

    if (allocated( occ ) )   deallocate( occ )
    if (allocated( ham ) )   deallocate( ham )
    allocate( occ( sites, dim ) )
    allocate( ham( dim, 2 ) )
    occ = 0 
    ham = 0 
    do j = 1, dim
        counter_v     = 0
        counter_v2    = 0
        do m = 0, sites - 1 !m goes through all digits of each configuration
            if ( btest( basis( j ), m ) ) occ( m + 1, j ) = occ( m + 1, j ) + 1
        end do
        do s = 1, nbonds !Loop runs over all nn-bonds
            ! if ( btest( basis( j ), bsites( 1, s ) - 1 ) .and. &
            !      btest( basis( j ), bsites( 2, s ) - 1 ) ) counter_v = counter_v + 1
            if ( btest( basis( j ), bsites( 1, s ) - 1 ) .and. &
                btest( basis( j ), bsites( 2, s ) - 1 ) ) counter_v = counter_v + 1
            if ( .not. (btest( basis( j ), bsites( 1, s ) - 1 )) .and. &
                .not. (btest( basis( j ), bsites( 2, s ) - 1 )) ) counter_v = counter_v + 1
            if ( btest( basis( j ), bsites( 1, s ) - 1 ) .and. &
                .not. (btest( basis( j ), bsites( 2, s ) - 1 )) ) counter_v = counter_v - 1
            if ( .not. (btest( basis( j ), bsites( 1, s ) - 1 )) .and. &
                btest( basis( j ), bsites( 2, s ) - 1 ) ) counter_v = counter_v - 1
        end do

        do s = 1, nbonds2 !Loop runs over all nnn-bonds

            
            ! if ( btest( basis( j ), bsites2( 1, s ) - 1 ) .and. &
            !     btest( basis( j ), bsites2( 2, s ) - 1 ) ) counter_v2 = counter_v2 + 1
            if ( btest( basis( j ), bsites2( 1, s ) - 1 ) .and. &
                btest( basis( j ), bsites2( 2, s ) - 1 ) ) counter_v2 = counter_v2 + 1
            if ( .not. (btest( basis( j ), bsites2( 1, s ) - 1 )) .and. &
                .not. (btest( basis( j ), bsites2( 2, s ) - 1 )) ) counter_v2 = counter_v2 + 1
            if ( btest( basis( j ), bsites2( 1, s ) - 1 ) .and. &
                .not. (btest( basis( j ), bsites2( 2, s ) - 1 )) ) counter_v2 = counter_v2 - 1
            if ( .not. (btest( basis( j ), bsites2( 1, s ) - 1 )) .and. &
                btest( basis( j ), bsites2( 2, s ) - 1 ) ) counter_v2 = counter_v2 - 1
        end do
        
        ! do s = 1, nbonds !Loop runs over all nn-bonds
        !     if ( .not. (btest( basis( j ), bsites( 1, s ) - 1 )) .and. &
        !          .not. (btest( basis( j ), bsites( 2, s ) - 1 )) ) counter_v = counter_v + 1
        ! end do
        ! do s = 1, nbonds2 !Loop runs over all nnn-bonds
        ! if ( .not. (btest( basis( j ), bsites2( 1, s ) - 1 )) .and. &
        !      .not.btest( basis( j ), bsites2( 2, s ) - 1 ) ) counter_v2 = counter_v2 + 1
        ! end do
        ham( j, 1 ) = ham( j, 1 ) + counter_v
        ham( j, 2 ) = ham( j, 2 ) + counter_v2
    end do

    nnz = dim

    print*,'Generated diagonal Hamiltonian.'
    print*, ''

end subroutine diagonal

subroutine diagonal_d(unit, parameters, sites, nbonds, nbonds2, dim, bsites, bsites2, basis, ham, occ, nnz)

    implicit none

    integer(kind=8), intent(in) :: dim, basis(dim)
    integer, intent(in) :: unit, sites, nbonds, nbonds2
    integer, intent(in) :: bsites(2, nbonds), bsites2(2, nbonds2)
    character(len=*), intent(in) :: parameters
    integer, intent(out) :: nnz
    integer, allocatable, intent(out) :: occ(:,:)
    double precision, allocatable, intent(out) :: ham(:,:)

    double precision :: counter_v = 0.d0, counter_v2 = 0.d0
    integer :: j = 0, m = 0, s = 0

    if (allocated( occ ) )   deallocate( occ )
    if (allocated( ham ) )   deallocate( ham )
    allocate( occ( sites, dim ) )
    allocate( ham( dim, 2 ) )
    occ = 0 
    ham = 0.d0  
    do j = 1, dim
        counter_v     = 0
        counter_v2    = 0
        do m = 0, sites - 1 !m goes through all digits of each configuration
            if ( btest( basis( j ), m ) ) occ( m + 1, j ) = occ( m + 1, j ) + 1
        end do
        do s = 1, nbonds !Loop runs over all nn-bonds
            if ( btest( basis( j ), bsites( 1, s ) - 1 ) .and. &
                 btest( basis( j ), bsites( 2, s ) - 1 ) ) counter_v = counter_v + 1
            if ( .not. (btest( basis( j ), bsites( 1, s ) - 1 )) .and. &
                 .not. (btest( basis( j ), bsites( 2, s ) - 1 )) ) counter_v = counter_v + 1
            if ( btest( basis( j ), bsites( 1, s ) - 1 ) .and. &
                 .not. (btest( basis( j ), bsites( 2, s ) - 1 )) ) counter_v = counter_v - 1
            if ( .not. (btest( basis( j ), bsites( 1, s ) - 1 )) .and. &
                 btest( basis( j ), bsites( 2, s ) - 1 ) ) counter_v = counter_v - 1
        end do
        do s = 1, nbonds2 !Loop runs over all nnn-bonds
            if ( btest( basis( j ), bsites2( 1, s ) - 1 ) .and. &
                 btest( basis( j ), bsites2( 2, s ) - 1 ) ) counter_v2 = counter_v2 + 1
            if ( .not. (btest( basis( j ), bsites2( 1, s ) - 1 )) .and. &
                 .not. (btest( basis( j ), bsites2( 2, s ) - 1 )) ) counter_v2 = counter_v2 + 1
            if ( btest( basis( j ), bsites2( 1, s ) - 1 ) .and. &
                 .not. (btest( basis( j ), bsites2( 2, s ) - 1 )) ) counter_v2 = counter_v2 - 1
            if ( .not. (btest( basis( j ), bsites2( 1, s ) - 1 )) .and. &
                 btest( basis( j ), bsites2( 2, s ) - 1 ) ) counter_v2 = counter_v2 - 1
        end do
        ! do s = 1, nbonds !Loop runs over all nn-bonds
            ! if ( .not. (btest( basis( j ), bsites( 1, s ) - 1 )) .and. &
            !      .not. (btest( basis( j ), bsites( 2, s ) - 1 )) ) counter_v = counter_v + 1
        ! end do
        ! do s = 1, nbonds2 !Loop runs over all nnn-bonds
        ! if ( .not. (btest( basis( j ), bsites2( 1, s ) - 1 )) .and. &
        !      .not.btest( basis( j ), bsites2( 2, s ) - 1 ) ) counter_v2 = counter_v2 + 1
        ! end do
        ham( j, 1 ) = ham( j, 1 ) + counter_v
        ham( j, 2 ) = ham( j, 2 ) + counter_v2
    end do

    nnz = dim

    print*,'Generated diagonal Hamiltonian.'
    print*, ''

end subroutine diagonal_d

subroutine diagonal_irrep2D(sites, dim, trisites, hexsites, orbsize, orbits, phases, norm, ham, occ)

    implicit none

    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: sites, orbsize
    integer, allocatable, intent(in) :: trisites(:,:), hexsites(:,:)
    integer(kind=8), allocatable, intent(in) :: orbits(:,:,:)
    double precision, allocatable, intent(in) :: norm(:,:)
    double complex, allocatable, intent(in) :: phases(:,:,:)
    
    integer, allocatable, intent(out) :: occ(:,:)
    double complex, allocatable, intent(out) :: ham(:,:)

    integer :: j, m, s, c, o, site1, site2 
    integer(kind=8) :: state, rowIndx
    double precision :: cntr_v = 0.d0, cntr_v2 = 0.d0
    

    if(allocated(occ)) deallocate(occ)
    if(allocated(ham)) deallocate(ham)
    allocate(occ(sites,2*dim))
    allocate(ham(2*dim,2))
    occ = 0 
    ham = 0.d0  
    do j = 1, dim !Representatives
        do c = 1, 2 !First and second basis state of 2D irrep 
            cntr_v  = 0
            cntr_v2 = 0
            rowIndx = 2*(j-1) + c
            do o = 1, orbsize
                state = orbits(j, o, c)
                do m = 0, sites - 1 !m goes through all digits of each configuration
                    if(btest(state, m)) occ(m + 1, rowIndx) = occ(m + 1, rowIndx) + 1
                end do

                do s = 1, 3 !Loop runs over the three different NN-bonds
                    site1 = trisites(1,s)
                    site2 = trisites(2,s)
                
                    if(btest(state, site1-1) .eqv. btest(state, site2-1)) cntr_v = cntr_v + 1
                    if(btest(state, site1-1) .neqv. btest(state, site2-1)) cntr_v = cntr_v - 1
                end do

                site1 = hexsites(1, s)
                site2 = hexsites(2, s)
                
                if(btest(state, site1-1) .eqv. btest(state, site2-1))  cntr_v2 = cntr_v2 + 1
                if(btest(state, site1-1) .neqv. btest(state, site2-1)) cntr_v2 = cntr_v2 - 1
                
                ham(rowIndx, 1) = ham(rowIndx, 1) + phases(j, o, c) * cntr_v / sqrt(norm(j, c))
                ham(rowIndx, 2) = ham(rowIndx, 2) + phases(j, o, c) * cntr_v2 / sqrt(norm(j, c))
            end do !Loop over orbit of irrep basis states (o)
        end do !Loop over irrep basis states (c)
    end do !Loop over representatives (j)

    print*,'Generated diagonal Hamiltonian.'
    print*, ''

end subroutine diagonal_irrep2D

subroutine diagonal_p(parity, p1, p2, p3, refl1, refl2, refl3, sites, nbonds, nbonds2, dim, bsites, bsites2, basis, ham, occ, nnz)

    implicit none

    integer(kind=8), intent(in) :: dim, basis(dim)
    integer, intent(in) :: parity, p1, p2, p3, sites, nbonds, nbonds2
    integer, intent(in) :: bsites(2, nbonds), bsites2(2, nbonds2), refl1(sites), refl2(sites), refl3(sites)
    integer, intent(out) :: nnz
    integer, allocatable, intent(out) :: occ(:,:)
    integer, allocatable, intent(out) :: ham(:,:)

    integer :: counter_v = 0.d0, counter_v2 = 0.d0
    integer(kind=8) :: j = 0, m = 0, s = 0, state = 0, nrefl = 0, p = 0
    integer :: sign = 1, info = 0
    

    if(allocated(occ)) deallocate(occ)
    if(allocated(ham)) deallocate(ham)
    allocate(occ(sites, dim))
    allocate(ham(dim, 2))
    occ = 0 
    ham = 0.d0  
    if(parity == 0) then 
        nrefl = 1
    else 
        nrefl = 1 + abs(p1) + abs(p2) + abs(p3) + abs(p1*p2) + abs(p1*p3)
    end if 

    do j = 1, dim
        counter_v  = 0
        counter_v2 = 0
        do p = 1, nrefl
            
            if(p == 2) call reflect(1, basis(j), sites, refl1, sign, info, state)
            if(p == 3) call reflect(1, basis(j), sites, refl2, sign, info, state)
            if(p == 4) call reflect(1, basis(j), sites, refl3, sign, info, state)
            if(p == 5) then 
                call reflect(1, basis(j), sites, refl1, sign, info, state)
                call reflect(1, state, sites, refl2, sign, info, state)
            else if(p == 6) then 
                call reflect(1, basis(j), sites, refl1, sign, info, state)
                call reflect(1, state, sites, refl3, sign, info, state)
            end if 
            if(state == basis(j)) cycle 
            if(p == 1) state = basis(j)

            do m = 0, sites - 1 !m goes through all digits of each configuration
                if(btest(state, m)) occ(m + 1, j) = occ(m + 1, j) + 1
            end do
            do s = 1, nbonds !Loop runs over all nn-bonds
                if(btest(state, bsites(1, s) - 1) .and. &
                    btest(state, bsites(2, s) - 1)) counter_v = counter_v + 1
                if(.not.(btest(state, bsites(1, s) - 1)) .and. &
                    .not.(btest(state, bsites(2, s) - 1))) counter_v = counter_v + 1
                if(btest(state, bsites(1, s) - 1) .and. &
                    .not.(btest(state, bsites(2, s) - 1))) counter_v = counter_v - 1
                if(.not.(btest(state, bsites(1, s) - 1)) .and. &
                    btest(state, bsites(2, s) - 1)) counter_v = counter_v - 1
            end do
            do s = 1, nbonds2 !Loop runs over all nnn-bonds
                if(btest(state, bsites2(1, s) - 1) .and. &
                    btest(state, bsites2(2, s) - 1)) counter_v2 = counter_v2 + 1
                if(.not.(btest(state, bsites2(1, s) - 1)) .and. &
                    .not.(btest(state, bsites2(2, s) - 1))) counter_v2 = counter_v2 + 1
                if(btest(state, bsites2(1, s) - 1) .and. &
                    .not.(btest(state, bsites2(2, s) - 1))) counter_v2 = counter_v2 - 1
                if(.not.(btest(state, bsites2(1, s) - 1)) .and. &
                    btest(state, bsites2(2, s) - 1)) counter_v2 = counter_v2 - 1
            end do

            ham(j, 1) = ham(j, 1) + counter_v
            ham(j, 2) = ham(j, 2) + counter_v2

        end do 
    end do

    nnz = dim

    print*,'Generated diagonal Hamiltonian.'
    print*, ''

end subroutine diagonal_p

  



!----------------------------------!
!         Test symmetry            !
!----------------------------------!

subroutine findstate2(dim, target, array, loc)

    implicit none

    integer(kind=8), intent(in) :: target(2), dim
    integer, intent(in) :: array(dim, 2)
    integer(kind=8), intent(out) :: loc
    

    integer(kind=8) :: i, left, right, mean

    left = 1
    right = dim
    do while (left <= right)
        mean = floor((left + right)/2.d0)
        if (array(mean, 1) < target(1)) then
            left = mean + 1
        else if (array(mean, 1) > target(1)) then
            right = mean - 1
        else
            loc = mean
            exit
        end if
    end do
    
    do i = 0, dim - loc - 1
        if(array(loc + i, 2) == target(2) .and. array(loc + i, 1) == target(1)) then 
            loc = loc + i 
            return
        end if 
    end do 

    loc = -1 !If no entry is found
    return


end subroutine findstate2

subroutine testsymm_sparse2(symmetric, nz, rc, ham)

    implicit none

    integer(kind=8), intent(in) :: nz 
    integer, intent(in) :: rc(nz,2)
    double precision, intent(in) :: ham(nz)
    
    logical, intent(out) :: symmetric
    
    integer(kind=8) :: i = 0, loc = 0, target(2) = 0
    double precision :: eps = 10.**(-10)
    
    symmetric = .true. 
    do i = 1, nz 
        target = rc(i,1:2) 
        call findstate2(nz, target, rc, loc)
        if(abs(ham(i) - ham(loc)) > eps) then 
            print* ,ham(i), rc(i,1), rc(i,2), 'ham(i), rc(i,1), rc(i,2)' 
            print* ,ham(loc), rc(loc,1), rc(loc,2), 'ham(loc), rc(loc,1), rc(loc,2)' 
            symmetric = .false. 
            exit 
        end if 
    end do 

    if(.not. symmetric) then
        print*, 'MATRIX IS NON-SYMMETRIC!'
        print*, ''
        error stop 
    else
        print*, 'MATRIX IS SYMMETRIC'
        print*, ''
    end if
    return



end subroutine testsymm_sparse2

subroutine testsymm_sparse(symmetric, dim, nz, rc, ham, norm, prts, dplcts)

    implicit none

    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: nz 
    integer, intent(in) :: rc(nz,2)
    integer, intent(in) :: prts(nz), dplcts(nz)
    double precision, intent(in) :: ham(nz)
    double precision, intent(in) :: norm(dim)
    
    logical, intent(out) :: symmetric
    
    integer(kind=8) :: i = 0, j = 0
    double precision :: eps = 10.**(-10)
    integer, allocatable :: parmat(:,:), dupmat(:,:)
    double precision, allocatable :: mat(:,:)
    
    if(allocated(mat)) deallocate(mat)
    if(allocated(parmat)) deallocate(parmat)
    if(allocated(dupmat)) deallocate(dupmat)
    allocate(mat(dim, dim))
    allocate(parmat(dim, dim))
    allocate(dupmat(dim, dim))
    mat = 0
    do i = 1, nz 
        mat(rc(i,1),rc(i,2)) = ham(i) 
        parmat(rc(i,1),rc(i,2)) = prts(i)
        dupmat(rc(i,1),rc(i,2)) = dplcts(i)
        if(rc(i,1) == 0 .or. rc(i,1) > dim ) print*, rc(i,1), 'rc(i,1)'
        if(rc(i,2) == 0 .or. rc(i,2) > dim ) print*, rc(i,2), 'rc(i,2)'
    end do 
 
    symmetric=.true.
    do i = 1, dim
        do j = 1, dim
            if(abs(mat(i,j) - mat(j,i)) > eps) then
                print*, i, 'i'
                print*, j, 'j'
                print*,'mat(i, j)', mat(i, j)
                print*,'mat(j, i)', mat(j, i)
                print*,'Difference', abs(mat(i,j) - mat(j,i))
                print*,'i-th norm ', norm(i) 
                print*,'j-th norm ', norm(j) 
                print*,'Parity ij', parmat(i, j)
                print*,'Parity ji', parmat(j, i)
                print*,'dplcts ij', dupmat(i, j)
                print*,'dplcts ji', dupmat(j, i)
                symmetric=.false.
                exit
            enDif
        end do
    end do
    if(.not. symmetric) then
        print*, 'MATRIX IS NON-SYMMETRIC!'
        print*, ''
        error stop 
    else
        print*, 'MATRIX IS SYMMETRIC'
        print*, ''
    end if
    return



end subroutine testsymm_sparse

subroutine testhermitian_sparse(hermit, irrep, symmetrize, dim, nz, rc, ham, norm, prts, dplcts)!, norm2D, orbits, phases

    implicit none

    integer(kind=8), intent(in) :: dim
    !integer(kind=8), allocatable, intent(in) :: orbits(:,:,:)
    integer, intent(in) :: nz, symmetrize
    integer, intent(in) :: rc(nz,2)
    integer, intent(in) :: prts(nz)
    integer, intent(in) :: dplcts(nz)
    double precision, allocatable, intent(in) :: norm(:)!, norm2D(:,:)
    double complex, allocatable, intent(in) :: ham(:)!, phases(:,:,:) 
    character, intent(in) :: irrep*2
    
    logical, intent(out) :: hermit
    
    integer(kind=8) :: i, j, mdim
    double precision :: eps = 10.**(-10)
    integer, allocatable :: parmat(:,:), dupmat(:,:)
    double complex, allocatable :: mat(:,:)

    if(irrep(1:1) .ne. "E") then 
        mdim = dim 
    else 
        mdim = 2 * dim 
    end if
    if(allocated(parmat)) deallocate(parmat)
    if(allocated(dupmat)) deallocate(dupmat)
    if(allocated(mat))    deallocate(mat)
    allocate(parmat(mdim, mdim))
    allocate(dupmat(mdim, mdim))
    allocate(mat(mdim, mdim))
    parmat = 0 
    dupmat = 0 
    mat = 0


    do i = 1, nz 
        ! print*, i, 'i'
        if(i == 25426 .or. i == 25488) then 
            print* , 'nz' ,nz
            print* ,i, 'i'
            print* ,ham(i), rc(i,1), rc(i,2), 'ham(i), row, col'
            pause 
        end if 

        mat(rc(i,1),rc(i,2)) = ham(i)

        ! parmat(rc(i,1),rc(i,2)) = prts(i)
        ! dupmat(rc(i,1),rc(i,2)) = dplcts(i)
        if(rc(i,1) == 0 .or. rc(i,1) > mdim ) print*, rc(i,1),i, 'rc(i,1), i'
        if(rc(i,2) == 0 .or. rc(i,2) > mdim ) print*, rc(i,2),i, 'rc(i,2), i'
    end do 

    
    hermit=.true.
    do i = 1, mdim
        do j = 1, mdim
            if(dble(conjg(mat(i,j))) - dble(mat(j,i)) > eps .or. aimag(conjg(mat(i,j))) - aimag(mat(j,i)) > eps) then
                if(dble(conjg(mat(i,j))) - dble(mat(j,i)) > eps) then
                    print*, 'RE(M*(i,j))', dble(conjg(mat(i,j)))
                    print*, 'RE(M(j,i))', dble(mat(j,i))
                    print*, 'Difference', abs(dble(conjg(mat(i,j))) - dble(mat(j,i)))
                    if(irrep(1:1) .ne. "E") then 
                        print*, 'i-th norm ', norm(i) 
                        print*, 'j-th norm ', norm(j)
                    end if 
                    print*, ''
                end if

                if(aimag(conjg(mat(i,j))) - aimag(mat(j,i)) > eps) then
                    print*, 'IM(M*(i,j))', aimag(conjg(mat(i,j)))
                    print*, 'IM(M(j,i))', aimag(mat(j,i))
                    print*, 'Difference', abs(aimag(conjg(mat(i,j))) - aimag(mat(j,i)))
                    if(irrep(1:1) .ne. "E") then 
                        print*, 'i-th norm ', norm(i) 
                        print*, 'j-th norm ', norm(j)
                    end if 
                    print*, ''
                end if

                print*, 'i', i
                print*, 'j', j
                print*, 'M(j,i)', mat(j,i)
                print*, 'M(i,j)', mat(i,j)
                print*,  ''

                hermit=.false.
                
            enDif
        end do
    end do
    if(.not. hermit) then
        print*, 'MATRIX IS NON-HERMITIAN!'
        print*, ''
        error stop
    else
        print*, 'MATRIX IS HERMITIAN'
        print*, ''
    end if
    return

end subroutine testhermitian_sparse

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

!----------------------------------------------!
!               Diagonalization                !
!----------------------------------------------!

subroutine diagonalization(dir, conf, nev, ncv, full, v1, v2, threads, parameters, type, basis, bsites, hexsites, occ, nOff, nDi, nDi_off, hamOff, hamDi, hamOff_dp, hamDi_d, hamOff_dc, hamDi_c, hamDi_off_c, ham, ham_dc, ham_d, ham_dc, norm, rcoff, rcdi, rc, prts, dplcts, nnz, ndeg, unit, nest, mode, energies, eigstate, eigstate_dc, gs, gs_c)
    
    use variables
    implicit none 
    
    integer, intent(in) :: conf, threads
    integer, intent(inout) :: nev, ncv, full, nest, unit, nOff, nDi, nDi_off
    ! integer(kind=8), intent(in) :: dim
    double precision, intent(in) :: v1, v2
    character(len=*), intent(in) :: dir, parameters, type, mode
    integer, allocatable, intent(inout) :: occ(:,:), hamOff(:,:), rcoff(:,:), rcdi(:), hamDi(:,:), prts(:), dplcts(:), bsites(:,:), hexsites(:,:)
    integer(kind=8), allocatable, intent(inout) :: basis(:)
    double precision, allocatable, intent(inout) :: hamOff_dp(:), hamDi_dp(:), norm(:)
    double complex, allocatable, intent(inout) :: hamOff_dc(:), hamDi_c(:,:), hamDi_off_c(:)
    integer, intent(out) :: nnz, ndeg  
    integer, allocatable, intent(out) :: rc(:,:) 
    double precision, allocatable, intent(out) :: energies(:), eigstate(:,:), gs(:), ham(:), ham_dp(:,:)
    double complex, allocatable, intent(out) :: eigstate_dc(:,:), gs_c(:), ham_dc(:), ham_dc(:,:)

    integer :: i = 0 
    double precision :: Emin = 0.d0, Emax = 0.d0
    logical :: symmetric, append 
    
    if((nDis > 1) .and. (dis .ne. 0.d0) .and. (conf > 1)) then 
        append = .true.
    else 
        append = .false.
    end if 
    
    100 format(1000(x,A,x,F6.4,x,A,x,F6.4))

    if(full == 0) then
        if(type == "R") then 
            if(ti == 0) then 
                call unify(t, v1, v2, dis, mass, pattern, sites, occ, nOff, nDi, hamOff, hamDi, ham, rc, nnz)
                if(arpack == 1) then 
                    if(otf == 0) then
                        print*,''
                        write(*,100) '---------------- START ARPACK DIAGONALIZATION: V1 =', v1, ' V2 =', v2, ' -------------------'
                        print*,''
                        call lanczos_d(threads, dim, nnz, nev, ncv, nest, mode, rvec, ham, rc, energies, eigstate)          
                    else if(otf == 1) then 
                        print*,''
                        write(*,100) '---------------- START ARPACK DIAGONALIZATION ON THE FLY: V1 =', v1, ' V2 =', v2, ' -------------------'
                        print*,''
                        call lanczos_d_otf(dir, parameters, unit, threads, dim, sites, sites/2*3, sites*3, nev, ncv, nest, t, v1, v2, mode, rvec, basis, bsites, hexsites, energies, eigstate)
                    end if 
                    if(rvec) call check_d(dim, nnz, nev, nest, energies, eigstate)                  
                    call save(dir, "R", parameters, append, conf, unit, 3, dim, states, nev, nest, nDis, rvec, energies, eigstate)
                end if 
                if(mkl == 1) then 
                    print*,''
                    write(*,100) '---------------- START MKL DIAGONALIZATION: V1 =', v1, ' V2 =', v2, ' -------------------'
                    call lanczos_mkl(dim, nev, ncv, nest, nnz, ham, rc, energies, eigstate) 
                    if(rvec) call check_d(dim, nnz, nev, nest, energies, eigstate)                  
                    call save(dir, "R", parameters, append, conf, unit, 2, dim, states, nev, nest, nDis, rvec, energies, eigstate)
                end if 
                if(degeneracy == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate, ndeg, gs)     
                if(degeneracy == 2 .and. rvec) call qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, eigstate, ndeg, gs)
            else if(ti == 1) then 
                call unify_dp(v1, v2, dis, mass, pattern, sites, occ, nOff, nDi, nDi_off, hamOff_dp, rcoff, rcdi, hamDi_d, hamDi, ham, rc, nnz)
                call testsymm_sparse(symmetric, dim, nnz, rc, ham, norm, prts, dplcts)
                ! call testsymm_sparse2(symmetric, nnz, rc, ham)
                if(arpack == 1) then 
                    print*,''
                    write(*,100) '---------------- START ARPACK DIAGONALIZATION: V1 =', v1, ' V2 =', v2, ' -------------------'
                    print*,''
                    call lanczos_d(threads, dim, nnz, nev, ncv, nest, mode, rvec, ham, rc, energies, eigstate)
                    if(degeneracy == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate, ndeg, gs)
                    if(degeneracy == 2 .and. rvec) call qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, eigstate, ndeg, gs)
                    if(rvec) call check_d(dim, nnz, nev, nest, energies, eigstate, ndeg = ndeg, rc = rc, mat = ham) 
                    call save(dir, "R", parameters, append, conf, unit, 3, dim, states, nev, nest, nDis, rvec, energies, eigstate)    
                    if(degeneracy == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate, ndeg, gs)
                    if(degeneracy == 2 .and. rvec) call qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, eigstate, ndeg, gs)
                end if 
                if(mkl == 1) then 
                    print*,''
                    write(*,100) '---------------- START MKL DIAGONALIZATION: V1 =', v1, ' V2 =', v2, ' -------------------'
                    print*,''
                    call lanczos_mkl(dim, nev, ncv, nest, nnz, ham, rc, energies, eigstate)                
                    if(degeneracy == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate, ndeg, gs)
                    if(degeneracy == 2 .and. rvec) call qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, eigstate, ndeg, gs)
                    if(rvec) call check_d(dim, nnz, nev, nest, energies, eigstate, ndeg = ndeg, rc = rc, mat = ham) 
                    call save(dir, "R", parameters, append, conf, unit, 2, dim, states, nev, nest, nDis, rvec, energies, eigstate)
                    if(degeneracy == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate, ndeg, gs)     
                    if(degeneracy == 2 .and. rvec) call qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, eigstate, ndeg, gs)    
                end if 
            end if
            if(feast == 1) then 
                Emin = dble(energies(1)) - 0.01
                Emax = dble(energies(min(nev, nevmax)))+0.01
                print*,'Emin',Emin
                print*,'Emax',Emax
                print*,'Emax-Emin',Emax-Emin                
                
                print*,''
                write(*,100) '---------------- START FEAST DIAGONALIZATION: V1 =', v1, ' V2 =', v2, ' -------------------'
                print*,''
                call dfeast(dim, nnz, nev0, nest, rc, ham, Emin, Emax, nev, energies, eigstate)

                if(degeneracy == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate, ndeg, gs)
                if(degeneracy == 2 .and. rvec) call qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, eigstate, ndeg, gs)
                if(rvec) call check_d(dim, nnz, nev, nest, energies, eigstate, ndeg = ndeg, rc = rc, mat = ham)
                call save(dir, "R", parameters, append, conf, unit, 1, dim, states, nev, nest, nDis, rvec, energies, eigstate)
                if(degeneracy == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate, ndeg, gs)
                if(degeneracy == 2 .and. rvec) call qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, eigstate, ndeg, gs)
            end if                          
        else if(ti == 1 .and. type == "C") then 
            if(irrep(1:1) == "E") then 
                call unify_c_2D(sites, nOff, nDi_off, dim, v1, v2, dis, mass, pattern, occ, hamOff_dc, rcoff, hamDi_off_c, rcdi, hamDi_c, ham_dc, rc, nnz)
            else
                call unify_comp(v1, v2, dis, mass, pattern, sites, occ, nOff, nDi, hamOff_dc, rcoff, hamDi, ham_dc, rc, nnz)   
            end if 
    
            call testhermitian_sparse(symmetric, irrep, symmetrize, dim, nnz, rc, ham_dc, norm, prts, dplcts)
            print*,''
            write(*,100) '---------------- START COMPLEX ARPACK DIAGONALIZATION: V1 =', v1, ' V2 =', v2, ' -------------------'
            print*,''

            call lanczos_c(threads, dim, nev, ncv, nest, mode, rvec, nnz, ham_dc, rc, energies, eigstate_dc)
            if(rvec) call check_c(dim, .False., nev, nest, energies, eigstate_dc)
            call save(dir, "C", parameters, append, conf, unit, 3, dim, states, nev, nest, nDis, rvec, energies, st_c=eigstate_dc)          
            if(degeneracy == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate_dc, ndeg, gs_c)      
            if(degeneracy == 2 .and. rvec) call qgsdeg_c(dir, unit, parameters, dim, nev, nest, energies, eigstate_dc, ndeg, gs_c)
            if(feast == 1) then 
                Emin = dble(energies(1))-0.01      
                Emax = dble(energies(min(nev, nevmax)))+0.01
                print*,'Emin',Emin
                print*,'Emax',Emax
                print*,'Emax-Emin',Emax-Emin
                print*,''
                print*,'---------------- START COMPLEX FEAST DIAGONALIZATION: V1 =', v1, ' V2 =', v2, ' -------------------'
                call cfeast(dim, nnz, nev0, nest, rc, ham_dc, Emin, Emax, nev, energies, eigstate_dc)
                if(rvec) call check_c(dim, .True., nev, nest, energies, eigstate_dc)
                call save(dir, "C", parameters, append, conf, unit, 1, dim, states, nev, nest, nDis, rvec, energies, st_c=eigstate_dc)           
                if(degeneracy == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate_dc, ndeg, gs_c)            
                if(degeneracy == 2 .and. rvec) call qgsdeg_c(dir, unit, parameters, dim, nev, nest, energies, eigstate_dc, ndeg, gs_c)                      
            end if 
        end if
    else if (full == 1) then
        print*,''
        write(*,100) '---------------- START EXACT DIAGONALIZATION: V1 =', v1, ' V2 =', v2, ' -------------------'
        print*,''
        if(ti == 0 .and. type == "R") then 
            call unify_dense(dim, t, v1, v2, dis, sites, occ, nOff, nDi, hamOff, hamDi, ham_dp)
            call testsymm(symmetric, dim, ham_dp)
            call exactdiag(rvec, dim, ham_d, energies)    
        else if(ti == 1 .and. type == "R") then 
            call unify_dense_d(dim, v1, v2, dis, sites, occ, nOff, nDi, nDi_off, hamOff_dp, rcoff, rcdi, hamDi_d, hamDi, ham_dp)
            call testsymm(symmetric, dim, ham_dp)
            call exactdiag(rvec, dim, ham_d, energies)
        else if(ti == 1 .and. type == "C") then 
            call unify_dense_c(dim, v1, v2, dis, sites, occ, nOff, nDi, hamOff_dc, rcoff, hamDi, ham_dc)
            call exactdiag_c(rvec, dim, ham_dc, energies)
        end if 
        if(type == "R") then 
            if(rvec) then 
                if(allocated(eigstate)) deallocate(eigstate)
                allocate(eigstate(dim, nest))
                eigstate = 0.d0 
                do i = 1, nest
                    eigstate(1:dim,i) = ham_dp(1:dim,i)
                end do
                if(degeneracy == 1) call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate, ndeg, gs)
                if(degeneracy == 2 .and. rvec) call qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, eigstate, ndeg, gs)
                if(rvec) call check_d(dim, nnz, nev, nest, energies, eigstate, ndeg = ndeg)
            end if
            call save(dir, "R", parameters, append, conf, unit, 0, dim, states, nev, nest, nDis, rvec, energies, eigstate)
        else if(type == "C") then 
            if(rvec) then 
                if(allocated(eigstate_dc)) deallocate(eigstate_dc)
                allocate(eigstate_dc(dim, nest))
                eigstate_dc = 0.d0 
                do i = 1, nest
                    eigstate_dc(1:dim,i) = ham_dc(1:dim,i)
                end do
                call gsdeg(dir, rvec, unit, parameters, dim, nev, nest, deg, energies, eigstate_dc, ndeg, gs_c)
                if(rvec) call check_c(dim, .True., nev, nest, energies, eigstate_dc)
            end if
            call save(dir, "C", parameters, append, conf, unit, 0, dim, states, nev, nest, nDis, rvec, energies, st_c = eigstate_dc)
        end if 
    end if

end subroutine diagonalization


!---------------------------------------------!
!            Unify sparse Hamiltonian         !
!---------------------------------------------!

subroutine unify(t, v1, v2, w, mass, pattern, sites, occ, nOff, nDi, hamOff, hamDi, ham, rc, nnz)

    implicit none
    ! save 
    integer, intent(in):: sites, nOff, nDi, occ(sites,nDi)
    integer, intent(in):: hamOff(nOff,3)
    integer, intent(in):: hamDi(nDi,2)
    double precision, intent(in) :: t, v1, v2, w, mass 
    character*2, intent(in) :: pattern

    integer, intent(out)::  nnz
    integer, allocatable, intent(out):: rc(:,:)
    double precision, allocatable, intent(out):: ham(:)

    integer :: j = 0, ab = 0 
    integer :: slstaggering(sites)
    double precision :: rand(sites)

    if( allocated(ham) ) deallocate(ham)
    if( allocated(rc) )  deallocate(rc)
    nnz = nOff + nDi !Number of non-zero elements
    allocate( ham(nnz) )
    allocate( rc(nnz,2) )
    
    ham = 0
    rc  = 0
    if (pattern == "AB") then 
        ab = 1 
    else if (pattern == "BA") then 
        ab = 0
    end if 

    !Staggered sublattice potential: +1 on A lattice, -1 on B lattice
    do j = 1, sites 
        slstaggering(j) = (-1)**(j-ab)
    end do 
    
    call random_number(rand)

    rand = 2 * (rand - 0.5)

    ham(1:nOff)  = - t * hamOff(1:nOff,1)
    rc(1:nOff,1) = hamOff(1:nOff,2)
    rc(1:nOff,2) = hamOff(1:nOff,3)
    do j = 1, nDi !Fill diagonal Hamiltonian
        ham(nOff+j)  = v1 * 0.25 * hamDi(j,1) + v2 * 0.25 * hamDi(j,2) + mass * dot_product(slstaggering, occ(1:sites,j)) + w * sum( rand * occ(1:sites,j) )
        ! if(j == 1) ham(nOff+j)  = ham(nOff+j) + 0.01
        rc(nOff+j,1) = j
        rc(nOff+j,2) = j
    end do
    
    return
    
end subroutine unify

subroutine unify_dp(v1, v2, w, mass, pattern, sites, occ, nOff, nDi, nDi2, hamOff, rcoff, rcdi, hamDi2, hamDi, ham, rc, nnz)

    implicit none

    integer, intent(in):: sites, nOff, nDi, nDi2, occ(sites,nDi)
    double precision, intent(in):: hamOff(nOff), hamDi2(nDi2)
    integer, intent(in):: rcoff(nOff,2), rcdi(nDi2)
    integer, intent(in):: hamDi(nDi,2)
    double precision, intent(in) :: v1, v2, w, mass 
    character*2, intent(in) :: pattern

    integer, intent(out)::  nnz
    integer, allocatable, intent(out):: rc(:,:)
    double precision, allocatable, intent(out):: ham(:)

    integer :: j = 0, ab = 0 
    integer :: slstaggering(sites)
    double precision :: rand(sites)

    if(allocated(ham)) deallocate(ham)
    if(allocated(rc))  deallocate(rc)
    if(v1 == 0.d0 .and. v2 == 0.d0 .and. mass == 0.d0) then 
        nnz = nOff + nDi2
    else 
        nnz = nOff + max(nDi , nDi2) !Number of non-zero elements
    end if 
    allocate(ham(nnz))
    allocate(rc(nnz,2))   

    ham = 0.d0
    rc  = 0
    if(pattern == "AB") then 
        ab = 1 
    else if(pattern == "BA") then 
        ab = 0
    end if 

    !Staggered sublattice potential: +1 on A lattice, -1 on B lattice
    do j = 1, sites 
        slstaggering(j) = (-1)**(j-ab)
        ! print*,'Sublattice potential on site', j,' is ', slstaggering(j) 
    end do   
    
    call random_number(rand)
    rand = 2 * (rand - 0.5)
    ham(1:nOff)  = hamOff(1:nOff)
    rc(1:nOff,1) = rcoff(1:nOff,1)
    rc(1:nOff,2) = rcoff(1:nOff,2)

    if(v1 == 0.d0 .and. v2 == 0.d0 .and. mass == 0.d0 ) go to 110 !Don't fill diagonal with zeroes
    do j = 1, nDi !Fill diagonal Hamiltonian        
        ham(nOff+j)  = v1 * 0.25 * hamDi(j,1) + v2 * 0.25 * hamDi(j,2) + mass * dot_product(slstaggering, occ(1:sites,j)) + w * sum( rand * occ(1:sites,j) ) !- 3 * v1 * sum(occ(1:sites,j)) - 6 * v2 * sum(occ(1:sites,j))        
        rc(nOff+j,1) = j
        rc(nOff+j,2) = j
    end do
    
    
    if(nDi2 > 0) then 
        do j = 1, nDi2 
            ham(nOff+rcdi(j)) = ham(nOff+rcdi(j)) + hamDi2(j)
        end do 
    end if 
    110 continue 

    if(v1 == 0.d0 .and. v2 == 0.d0 .and. mass == 0.d0 ) then !Still need to fill diagonal entries from hopping 
        if(nDi2 > 0) then 
            do j = 1, nDi2 
                ham(nOff+ j) = ham(nOff+ j) + hamDi2(j)
                rc(nOff + j, 1) = rcdi(j) 
                rc(nOff + j, 2) = rcdi(j)
            end do 
        end if 
    end if 

    ! open(20,file=dir//'ham.dat')
    ! do j = 1, size(ham)
    ! write(20,*) ham(j), rc(j,1), rc(j,2) 
    ! end do 
    ! close(20)


    return

end subroutine unify_dp 

subroutine unify_comp(v1, v2, w, mass, pattern, sites, occ, nOff, nDi, hamOff, rcoff, hamDi, ham, rc, nnz)

    implicit none

    integer, intent(in):: sites, nOff, nDi, occ(sites,nDi)
    double complex, intent(in):: hamOff(nOff)
    integer, intent(in):: rcoff(nOff,2)
    integer, intent(in):: hamDi(nDi,2)
    double precision, intent(in) :: v1, v2, w, mass 
    character*2, intent(in) :: pattern

    integer, intent(out)::  nnz
    integer, allocatable, intent(out):: rc(:,:)
    double complex, allocatable, intent(out):: ham(:)

    integer :: j = 0, ab = 0 
    integer :: slstaggering(sites)
    double precision :: rand(sites)

    if( allocated(ham) ) deallocate(ham)
    if( allocated(rc) )  deallocate(rc)
    nnz = 0
    nnz = nOff + nDi !Number of non-zero elements
    allocate( ham(nnz) )
    allocate( rc(nnz,2) )    

    ham = 0
    rc  = 0
    if (pattern == "AB") then 
        ab = 1 
    else if (pattern == "BA") then 
        ab = 0
    end if 

    !Staggered sublattice potential: +1 on A lattice, -1 on B lattice
    do j = 1, sites 
        slstaggering(j) = (-1)**(j-ab)
    end do 
    
    call random_number(rand)
    rand = 2 * (rand - 0.5)
    ham(1:nOff)  = hamOff(1:nOff)
    rc(1:nOff,1) = rcoff(1:nOff,1)
    rc(1:nOff,2) = rcoff(1:nOff,2)
    do j = 1, nDi !Fill diagonal Hamiltonian
        ham(nOff+j)  = v1 * 0.25 * hamDi(j,1) + v2 * 0.25 * hamDi(j,2) + mass * dot_product(slstaggering, occ(1:sites,j)) !- 3 * v1 * sum(occ(1:sites,j)) - 6 * v2 * sum(occ(1:sites,j)) !+ w * sum( rand * occ(1:sites,j) )
        rc(nOff+j,1) = j
        rc(nOff+j,2) = j
    end do

    return

end subroutine unify_comp 

subroutine unify_c_2D(sites, nOff, nDi_off, dim, v1, v2, w, mass, pattern, occ, hamOff, rcoff, hamDi_off, rcdi, hamDi, ham, rc, nnz)

    implicit none

    integer, intent(in):: sites, nOff, nDi_off
    integer(kind=8), intent(in):: dim
    integer, allocatable, intent(in):: rcoff(:,:), rcdi(:), occ(:,:)
    double precision, intent(in) :: v1, v2, w, mass 
    double complex, allocatable, intent(in):: hamOff(:), hamDi_off(:), hamDi(:,:)
    character*2, intent(in) :: pattern

    integer, intent(out)::  nnz
    integer, allocatable, intent(out):: rc(:,:)
    double complex, allocatable, intent(out):: ham(:)

    integer(kind=8) :: j, ab, nDi
    integer :: slstaggering(sites)
    double precision :: rand(sites)

    if(allocated(ham)) deallocate(ham)
    if(allocated(rc))  deallocate(rc)
    nDi = 2*dim
    if(v1 == 0.d0 .and. v2 == 0.d0 .and. mass == 0.d0) then 
        nnz = nOff + nDi_off
        print* ,nnz, nDi_off, 'nnz, nDi_off'
    else 
        nnz = nOff + max(nDi, nDi_off) !Number of non-zero elements
    end if 
    allocate(ham(nnz))
    allocate(rc(nnz,2))    

    ham = 0.d0
    rc  = 0
    if (pattern == "AB") then 
        ab = 1 
    else if (pattern == "BA") then 
        ab = 0
    end if 

    !Staggered sublattice potential: +1 on A lattice, -1 on B lattice
    do j = 1, sites 
        slstaggering(j) = (-1)**(j-ab)
    end do 
    
    call random_number(rand)
    rand = 2 * (rand - 0.5)
    !Off diagonal matrix elements
    ham(1:nOff)  = hamOff(1:nOff)
    rc(1:nOff,1) = rcoff(1:nOff,1)
    rc(1:nOff,2) = rcoff(1:nOff,2)

    if(v1 == 0.d0 .and. v2 == 0.d0 .and. mass == 0.d0 ) go to 110 !Don't fill diagonal with zeroes
    do j = 1, nDi !Fill diagonal Hamiltonian
        ham(nOff+j)  = v1 * 0.25 * hamDi(j,1) + v2 * 0.25 * hamDi(j,2) + mass * dot_product(slstaggering, occ(1:sites,j)) + w * sum(rand * occ(1:sites,j))
        rc(nOff+j,1) = j
        rc(nOff+j,2) = j
    end do

    if(nDi_off > 0) then !Fill diagonal elements generated by hopping 
        do j = 1, nDi_off 
            ham(nOff+rcdi(j)) = ham(nOff+rcdi(j)) + hamDi_off(j)
        end do 
    end if 
    110 continue 

    if(v1 == 0.d0 .and. v2 == 0.d0 .and. mass == 0.d0 ) then !Still need to fill diagonal entries from hopping 
        if(nDi_off > 0) then 
            do j = 1, nDi_off 
                ham(nOff+j) = ham(nOff+j) + hamDi_off(j)
                rc(nOff+j, 1) = rcdi(j) 
                rc(nOff+j, 2) = rcdi(j)
            end do 
        end if 
    end if 


    return

end subroutine unify_c_2D 

subroutine unify_dense(dim, t, v1, v2, w, sites, occ, nOff, nDi, hamOff, hamDi, ham)

    implicit none

    integer(kind=8), intent(in) :: dim 
    integer, intent(in):: sites, nOff, nDi, occ(sites,*)
    integer, intent(in):: hamOff(nOff,3)
    integer, intent(in):: hamDi(nDi,2)
    double precision, intent(in) :: t, v1, v2, w

    double precision, allocatable, intent(out):: ham(:,:)

    integer :: j = 0
    double precision :: rand(sites) 
    logical :: symm 

    if( allocated(ham) ) deallocate(ham)
    allocate( ham(dim, dim) )
    ham = 0
    
    call random_number(rand)
    rand = 2 * (rand - 0.5)
    do j = 1, nOff
        ham( hamOff(j,2), hamOff(j,3) )  = -t * hamOff(j,1)
    end do 
    do j = 1, nDi !Fill diagonal Hamiltonian
        ham(j,j) = v1 * hamDi(j,1) + v2 * hamDi(j,2) + w * sum( rand * occ(1:sites,j) )
    end do
    call testsymm(symm, dim, ham)
    
    return

end subroutine unify_dense

subroutine unify_dense_d(dim, v1, v2, w, sites, occ, nOff, nDi, nDi2, hamOff, rcoff, rcdi, hamDi2, hamDi, ham)
    
    implicit none

    integer(kind=8), intent(in) :: dim 
    integer, intent(in):: sites, nOff, nDi, nDi2, occ(sites,*)
    integer, intent(in):: rcoff(nOff,2), rcdi(nDi2)
    integer, intent(in):: hamDi(nDi,2)
    double precision, intent(in) :: v1, v2, w
    double precision, intent(in):: hamOff(nOff), hamDi2(nDi2)
    double precision, allocatable, intent(out):: ham(:,:)

    integer :: j = 0
    double precision :: rand(sites)
    logical :: symm 

    if(allocated(ham)) deallocate(ham)
    allocate(ham(dim, dim))
    ham = 0.d0 
    
    call random_number(rand)
    rand = 2 * (rand - 0.5)
    do j = 1, nOff
        ham(rcoff(j,1), rcoff(j,2)) = ham(rcoff(j,1), rcoff(j,2)) + hamOff(j)
    end do 
    do j = 1, nDi !Fill diagonal Hamiltonian
        ham(j,j) = v1 * hamDi(j,1) + v2 * hamDi(j,2) + w * sum( rand * occ(1:sites,j) )
    end do
    call testsymm(symm, dim, ham)
    
    do j = 1, nDi2 
        ham(rcdi(j), rcdi(j)) = ham(rcdi(j), rcdi(j)) + hamDi2(j)
    end do 

    return

end subroutine unify_dense_d

subroutine unify_dense_c(dim, v1, v2, w, sites, occ, nOff, nDi, hamOff, rcoff, hamDi, ham)
    
    implicit none

    integer(kind=8), intent(in) :: dim 
    integer, intent(in):: sites, nOff, nDi, occ(sites,*)
    integer, intent(in):: rcoff(nOff,2)
    integer, intent(in):: hamDi(nDi,2)
    double precision, intent(in) :: v1, v2, w
    double complex, intent(in):: hamOff(nOff)
    double complex, allocatable, intent(out):: ham(:,:)

    integer :: j = 0
    double precision :: rand(sites)
    logical :: symm 

    if(allocated(ham)) deallocate(ham)
    allocate(ham(dim, dim))
    ham = 0.d0 
    
    call random_number(rand)
    rand = 2 * (rand - 0.5)
    do j = 1, nOff
        ham(rcoff(j,1), rcoff(j,2)) = ham(rcoff(j,1), rcoff(j,2)) + hamOff(j)
    end do 
    do j = 1, nDi !Fill diagonal Hamiltonian
        ham(j,j) = v1 * hamDi(j,1) + v2 * hamDi(j,2) + w * sum( rand * occ(1:sites,j) )
    end do
    ! call testsymm(symm, dim, ham)
    
    ! do j = 1, nDi2 
    !     ham(rcdi(j), rcdi(j)) = ham(rcdi(j), rcdi(j)) + hamDi2(j)
    ! end do 

    return

end subroutine unify_dense_c

!-------------------------------------------!
!            Exact diagonalization          !
!-------------------------------------------!
    
subroutine exactdiag(eigvec, dim, mat, evals)
    
    implicit none
    logical, intent(in) :: eigvec
    integer(kind=8), intent(in) :: dim
    double precision, intent(inout) :: mat(dim, dim)
    double precision, allocatable, intent(out) :: evals(:)

    double precision, allocatable, save :: matloc(:,:)
    double precision, allocatable, save :: evalsloc(:)
    integer, parameter :: lwmax = 100000
    integer, save :: lda = 0, info = 0, lwork = 0
    character, save :: jobz
    double precision, allocatable, save :: work(:)

    external :: dsyev, dsyev_2stage

    !$OMP THREADPRIVATE (matloc, evalsloc, lda, info, lwork, jobz, work)
    lda   = 0
    info  = 0
    lwork = 0

    if (eigvec) then
        jobz = 'V'
    else
        jobz = 'N'
    end if

    if (allocated(evals)) deallocate(evals)
    allocate(evals(dim))
    evals = 0.0d0

    if (allocated(evalsloc)) deallocate(evalsloc)
    if (allocated(matloc)) deallocate(matloc)
    allocate(evalsloc(dim))
    allocate(matloc(dim, dim))
    evalsloc = 0.0d0
    matloc = mat

    if (dim > 1) then
        lda = dim
        if (allocated(work)) deallocate(work)
        allocate(work(lwmax))
        lwork = -1
        if (jobz == 'V') then
            !$omp critical
            call dsyev(jobz, 'U', dim, matloc, lda, evalsloc, work, lwork, info)
            !$omp end critical
        else if (jobz == 'N') then
            call dsyev_2stage ('N', 'U', dim, matloc, lda, evalsloc, work, lwork, info)
        end if
        !lwork = min(int(work(1)), lwmax)
        lwork = int(work(1))
        if (lwork < max(1, 6*(dim-1) ) ) lwork = max(1, 6*(dim-1) )
        if (allocated(work)) deallocate(work)
        allocate(work(lwork))
        work = 0.0d0

        !print*, 'Running diagoalization routine ...'
        if (jobz == 'V') then
            !$omp critical
            call dsyev(jobz, 'U', dim, matloc, lda, evalsloc, work, lwork, info)
            !$omp end critical
        else if (jobz == 'N') then
            call dsyev_2stage ('N', 'U', dim, matloc, lda, evalsloc, work, lwork, info)
        end if

        !print*, 'Finished diagoalization routine.'
        if (info .gt. 0 ) THEN
            write(*,*)'The algorithm failed to compute eigenvalues.'
            stop
        end if
        if (allocated(work)) deallocate(work)
    else
        evalsloc = matloc(1,1)
    end if
    ! print*, 'Eigenvalues '
    ! do i = 1, dim 
    !     print*, evalsloc(i)
    ! end do 
    evals = evalsloc
    mat   = matloc
    if(allocated(evalsloc)) deallocate(evalsloc)
    if(allocated(matloc)) deallocate(matloc)

    return

end subroutine exactdiag

subroutine exactdiag_c(eigvec, dim, mat, evals)
    
    implicit none
    logical, intent(in) :: eigvec
    integer(kind=8), intent(in) :: dim
    double complex, intent(inout) :: mat(dim, dim)
    double precision, allocatable, intent(out) :: evals(:)

    double complex, allocatable, save :: matloc(:,:)
    double precision, allocatable, save :: evalsloc(:)
    integer, parameter :: lwmax = 100000
    integer, save :: lda = 0, info = 0, lwork = 0
    character, save :: jobz
    double precision, allocatable, save :: rwork(:)
    double complex, allocatable, save :: work(:)

    external :: zheev

    !$OMP THREADPRIVATE (matloc, evalsloc, lda, info, lwork, jobz, work)
    lda   = 0
    info  = 0
    lwork = 0

    if (eigvec) then
        jobz = 'V'
    else
        jobz = 'N'
    end if

    if (allocated(evals)) deallocate(evals)
    allocate(evals(dim))
    evals = 0.0d0

    if (allocated(evalsloc)) deallocate(evalsloc)
    if (allocated(matloc)) deallocate(matloc)
    allocate(evalsloc(dim))
    allocate(matloc(dim, dim))
    evalsloc = 0.0d0
    matloc = mat

    if (dim > 1) then
        lda = dim
        if (allocated(work)) deallocate(work)
        allocate(work(lwmax))
        if (allocated(rwork)) deallocate(rwork)
        allocate(rwork(max(1, 3*dim-2)))
        lwork = -1
        if (jobz == 'V') then
            !$omp critical
            call zheev(jobz, 'U', dim, matloc, lda, evalsloc, work, lwork, rwork, info)
            !$omp end critical
        ! else if (jobz == 'N') then
        !     call dsyev_2stage ('N', 'U', dim, matloc, lda, evalsloc, work, lwork, info)
        end if
        !lwork = min(int(work(1)), lwmax)
        lwork = int(work(1))
        if (lwork < max(1, 6*(dim-1) ) ) lwork = max(1, 6*(dim-1) )
        if (allocated(work)) deallocate(work)
        allocate(work(lwork))
        work  = 0.0d0
        rwork = 0.d0 

        !print*, 'Running diagoalization routine ...'
        if (jobz == 'V') then
            !$omp critical
            call zheev(jobz, 'U', dim, matloc, lda, evalsloc, work, lwork, rwork, info)
            !$omp end critical
        ! else if (jobz == 'N') then
        !     call dsyev_2stage ('N', 'U', dim, matloc, lda, evalsloc, work, lwork, info)
        end if

        !print*, 'Finished diagoalization routine.'
        if (info .gt. 0 ) THEN
            write(*,*)'The algorithm failed to compute eigenvalues.'
            stop
        end if
        if (allocated(work)) deallocate(work)
        if (allocated(rwork)) deallocate(rwork)
    else
        evalsloc = matloc(1,1)
    end if
    ! print*, 'Eigenvalues '
    ! do i = 1, dim 
    !     print*, evalsloc(i)
    ! end do 
    evals = evalsloc
    mat   = matloc
    if(allocated(evalsloc)) deallocate(evalsloc)
    if(allocated(matloc)) deallocate(matloc)

    return

end subroutine exactdiag_c

!---------------------------------------------------!
!            Sparse diagonalization  (ARPACK)       !
!---------------------------------------------------!

subroutine lanczos_d(threads, dim, nnz, nev, ncv, nest, mode, rvec, ham, rc, evals, evecs)

    implicit none
    ! save
    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: threads, nev, ncv, nest, nnz
    integer, intent(in) :: rc(nnz,2)
    double precision, intent(in) :: ham(nnz)
    character(len=*), intent(in) :: mode
    logical, intent(in)          :: rvec
    double precision, allocatable, intent(out) :: evals(:), evecs(:,:)

    integer(kind=8), save :: maxn, n
    integer, save :: maxnev, maxncv, ldv, j
    integer, save :: ishfts, lworkl, maxitr, mode1, nconv, iparam(11), ipntr(11), ido, ierr, info
    integer, allocatable, save :: jao(:), iao(:)
    double precision, save :: zero = 0.0D+00
    double precision, save :: sigma, tol
    double precision, allocatable, save :: ax(:), workl(:), workd(:), v(:,:), resid(:), d(:,:), ao(:)
    character, save :: bmat*1, which*2
    logical, allocatable, save :: selector(:)
    double precision, external :: dnrm2
    external :: daxpy, dsaupd, dseupd, dmout
    intrinsic :: abs

    integer :: cntr = 0

    !$omp threadprivate(maxn, n, maxnev, maxncv, ldv, j, ishfts, lworkl, maxitr, mode1, nconv, iparam, ipntr, ido, ierr, info, jao, iao, zero, sigma, tol, ax, workl, workd, v, resid, d, ao, bmat, which, selector)


    if(allocated(evals)) deallocate(evals)
    if(allocated(evecs)) deallocate(evecs)
    allocate(evals(nev))
    allocate(evecs(dim, nest))
    evals = 0.d0 
    evecs = 0.d0    

    if (dim == 1) then
        evecs(1,1) = 1
        evals(1) = ham(1)
        print*, ''
        write(*,"('Only eigenvalue is ',f8.4)") ham(1)
    else
        maxn   = 10 + dim
        maxnev = 10 + nev
        maxncv = 10 + ncv
        ldv    = maxn
    
        if(allocated(workd)) deallocate(workd)
        if(allocated(ax)) deallocate(ax)
        if(allocated(d)) deallocate(d)
        if(allocated(resid)) deallocate(resid)
        if(allocated(selector)) deallocate(selector)
        if(allocated(workl)) deallocate(workl)
        if(allocated(v)) deallocate(v)
        allocate(ax(maxn),d(maxncv,2),resid(maxn),selector(maxncv),v(ldv,maxncv),workl(maxncv*(maxncv+8)),workd(3*maxn))

        !SET ARPACK PARAMETERS: Make sure that 1) maxn >= n 2) maxnev >= nev 3) maxncv >= ncv
        !Set dimensions for this problem.
        n     = int(dim, 4)
        bmat  = 'I'
        which = mode

        if ( n .gt. maxn ) then
            print *, ' ERROR: N is greater than MAXN '
            go to 9870
        else if ( nev .gt. maxnev ) then
            print *, ' ERROR: NEV is greater than MAXNEV '
            go to 9870
        else if ( ncv .gt. maxncv ) then
            print *, ' ERROR: NCV is greater than MAXNCV '
            go to 9870
        end if

        lworkl = ncv * ( ncv + 8 )
        tol    = zero
        info   = 0
        ido    = 0
        ishfts = 1
        maxitr = 300! 5000! 1000 !300
        mode1  = 1
        iparam = 0
        cntr   = 0 

        iparam(1) = ishfts
        iparam(3) = maxitr
        iparam(7) = mode1
        iparam(4) = 1 
        !MAIN LOOP (Reverse communication loop)! Repeatedly call DSAUPD and take actions inDicated by parameter !IDO until convergence is inDicated or MAXITR is exceeded.
        ierr     = 0
        ax       = 0.d0
        v        = 0.d0
        workl    = 0.d0
        workd    = 0.d0
        sigma    = 0.d0
        d        = 0.d0
        resid    = 0.d0
        ipntr    = 0.d0
        selector = 0.d0

        !Create sparse matrix in CSR format
        if(allocated(jao)) deallocate(jao)
        if(allocated(iao)) deallocate(iao)
        if(allocated(ao)) deallocate(ao)
        allocate(jao(nnz))
        allocate(iao(dim+1))
        allocate(ao(nnz))
        jao = 0 
        iao = 0 
        ao  = 0 
        
        call coocsr(dim, nnz, ham, rc(1:nnz,1), rc(1:nnz,2), ao, jao, iao) 
        
       !!$omp critical
        do
           !Sparse eigensolver routine
           call dsaupd ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
           if ( ido /= -1 .and. ido /= 1 ) then
               exit
           end if
           !Vector-sparse matrix multiplication                 
           cntr = cntr + 1
           call pamux (threads, dim, workd(ipntr(1)), workd(ipntr(2)), ao, jao, iao)
        end do
       !!$omp end critical
        print*,'Number of MV-multiplications',cntr
        call dsaupderrormessage(info)
        
        !Extract eigenvalues and -vectors
        ! !$omp critical
        call dseupd ( rvec, 'A', selector, d, v, ldv, sigma, bmat, n, which, nev, tol, &
                     resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, ierr )
        call dseupderrormessage(ierr)
        ! !$omp end critical
        if (ierr /= 0) then
           print *, 'DSSIMP - Fatal error! Error with DSEUPD, IERR = ', ierr, 'Check the documentation of DSEUPD.'
        else
           nconv =  iparam(5)
           do j = 1, nconv
               call pamux (threads, dim, v(1,j), ax, ao, jao, iao)
               call daxpy (n, -d(j,1), v(1,j), 1, ax, 1 )
               d(j,2) = dnrm2 (n, ax, 1)
               d(j,2) = d(j,2) / abs (d(j,1))
           end do
           call dmout (6, nconv, 2, d, maxncv, -6, 'Ritz values are named relative residuals' ) !Display computed residuals. 6: Output to screen Write(6, #internalnumber)! nconv: number of rows in the matrix d! 2: Number of columns in matrix d! maxncv: Leading dimension of the matrix data! -6: print the matrix d with iabs(-6) decimal digits per number! Use formatting indexed by -6 to print A
        end if
        


   !     10 continue
    ! !c
    ! !c        %---------------------------------------------%
    ! !c        | Repeatedly call the routine DSAUPD and take |
    ! !c        | actions inDicated by parameter IDO until    |
    ! !c        | either convergence is inDicated or maxitr   |
    ! !c        | has been exceeded.                          |
    ! !c        %---------------------------------------------%
    ! !c

    !     call dsaupd ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
    !     if (ido .eq. -1 .or. ido .eq. 1) then
    ! !c
    ! !c           %--------------------------------------%
    ! !c           | Perform matrix vector multiplication |
    ! !c           |              y <--- OP*x             |
    ! !c           | The user should supply his/her own   |
    ! !c           | matrix vector multiplication routine |
    ! !c           | here that takes workd(ipntr(1)) as   |
    ! !c           | the input, and return the result to  |
    ! !c           | workd(ipntr(2)).                     |
    ! !c           %--------------------------------------%

    !         call pamux (threads, dim, workd(ipntr(1)), workd(ipntr(2)), ao, jao, iao)

    ! !c           %-----------------------------------------%
    ! !c           | L O O P   B A C K to call DSAUPD again. |
    ! !c           %-----------------------------------------%

    !         go to 10

    !     end if

    ! !c     %----------------------------------------%
    ! !c     | Either we have convergence or there is |
    ! !c     | an error.                              |
    ! !c     %----------------------------------------%

    !     if ( info .lt. 0 ) then
    ! !c
    ! !c        %--------------------------%
    ! !c        | Error message. Check the |
    ! !c        | documentation in DSAUPD. |
    ! !c        %--------------------------%
    ! !c
    !         print *, ' '
    !         print *, ' Error with dsaupd, info = ', info
    !         print *, ' Check documentation in dsaupd '
    !         print *, ' '
    !     else
    ! !c
    ! !c        %-------------------------------------------%
    ! !c        | No fatal errors occurred.                 |
    ! !c        | Post-Process using DSEUPD.                |
    ! !c        |                                           |
    ! !c        | Computed eigenvalues may be extracted.    |
    ! !c        |                                           |
    ! !c        | Eigenvectors may also be computed now if  |
    ! !c        | desired.  (inDicated by rvec = .true.)    |
    ! !c        %-------------------------------------------%
    ! !c

    !         call dseupd ( rvec, 'All', selector, d, v, ldv, sigma, bmat, n, which, nev, tol, &
    !                       resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, ierr )

    ! !c        %----------------------------------------------%
    ! !c        | Eigenvalues are returned in the first column |
    ! !c        | of the two dimensional array D and the       |
    ! !c        | corresponDing eigenvectors are returned in   |
    ! !c        | the first NEV columns of the two dimensional |
    ! !c        | array V if requested.  Otherwise, an         |
    ! !c        | orthogonal basis for the invariant subspace  |
    ! !c        | corresponDing to the eigenvalues in D is     |
    ! !c        | returned in V.                               |
    ! !c        %----------------------------------------------%
    ! !c
 
    !         if ( ierr .ne. 0) then
    ! !c
    ! !c            %------------------------------------%
    ! !c            | Error conDition:                   |
    ! !c            | Check the documentation of DSEUPD. |
    ! !c            %------------------------------------%
    ! !c
    !             print *, ' '
    !             print *, ' Error with dseupd, info = ', ierr
    !             print *, ' Check the documentation of dseupd. '
    !             print *, ' '

    !         else

    !             nconv =  iparam(5)
    !             do j = 1, nconv

    ! !c
    ! !c               %---------------------------%
    ! !c               | Compute the residual norm |
    ! !c               |                           |
    ! !c               |   ||  A*x - lambda*x ||   |
    ! !c               |                           |
    ! !c               | for the NCONV accurately  |
    ! !c               | computed eigenvalues and  |
    ! !c               | eigenvectors.  (iparam(5) |
    ! !c               | inDicates how many are    |
    ! !c               | accurate to the requested |
    ! !c               | tolerance)                |
    ! !c               %---------------------------%
    ! !c
    !                 call pamux (threads, dim, v(1,j), ax, ao, jao, iao)
    !                 call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
    !                 d(j,2) = dnrm2(n, ax, 1)
    !                 d(j,2) = d(j,2) / abs(d(j,1))

    !             end do
    ! !c
    ! !c            %-------------------------------%
    ! !c            | Display computed residuals    |
    ! !c            %-------------------------------%
    ! !c
    !             call dmout(6, nconv, 2, d, maxncv, -6, &
    !                 'Ritz values and relative residuals')
    !             end if
   !         end if

            evals(1:nev) = d(1:nev,1)
            if (rvec) then
                do j= 1, min(iparam(5), nest) 
                    evecs(1:dim, j) = v(1:dim,j)
                end do
            end if
            
            9870 continue

        end if


        100 format(1000(F40.30))

        ! if(allocated(workd))    deallocate(workd)
        ! if(allocated(ax))       deallocate(ax)
        ! if(allocated(d))        deallocate(d)
        ! if(allocated(resid))    deallocate(resid)
        ! if(allocated(selector)) deallocate(selector)
        ! if(allocated(workl))    deallocate(workl)
        ! if(allocated(v))        deallocate(v)
        ! if(allocated(jao))      deallocate(jao)
        ! if(allocated(iao))      deallocate(iao)
        ! if(allocated(ao))       deallocate(ao)

        return 

end subroutine lanczos_d

subroutine lanczos_d_otf(dir, params, unit, threads, dim, sites, nbonds, nnnbonds, nev, ncv, nest, t1, v1, v2, mode, rvec, basis, bsites, hexsites, evals, evecs)

    implicit none

    integer(kind=8), intent(in) :: dim, basis(dim)
    integer, intent(in) :: unit, threads, sites, nbonds, nnnbonds, nev, ncv, nest
    integer, intent(in) :: bsites(2, nbonds), hexsites(2, nnnbonds)
    double precision, intent(in) :: t1, v1, v2
    character(len=*), intent(in) :: dir, params, mode
    logical, intent(in)          :: rvec
    double precision, allocatable, intent(out) :: evals(:), evecs(:,:)

    integer(kind=8)  :: maxn, n
    integer  :: maxnev, maxncv, ldv, j, cntr = 0 
    integer  :: ishfts, lworkl, maxitr, mode1, nconv, iparam(11), ipntr(11), ido, ierr, info
    double precision  :: zero = 0.0D+00
    double precision  :: sigma, tol
    double precision, allocatable  :: ax(:), workl(:), workd(:), v(:,:), resid(:), d(:,:)
    character  :: bmat*1, which*2
    logical, allocatable  :: selector(:)
    double precision, external :: dnrm2
    external :: daxpy, dsaupd, dseupd, dmout
    intrinsic :: abs

    !!$omp threadprivate(maxn, n, maxnev, maxncv, ldv, j, ishfts, lworkl, maxitr, mode1, nconv, iparam, ipntr, ido, ierr, info, zero, sigma, tol, ax, workl, workd, v, resid, d, bmat, which, selector)


    if(allocated(evals)) deallocate(evals)
    if(allocated(evecs)) deallocate(evecs)
    allocate(evals(nev))
    allocate(evecs(dim, nest))
    evals = 0.d0 
    evecs = 0.d0    

    if (dim == 1) then
        evecs(1,1) = 1
        evals(1) = 1
        print*, ''
        write(*,"('Only eigenvalue is ',f8.4)") evals(1)
    else
        maxn   = 10 + dim
        maxnev = 10 + nev
        maxncv = 10 + ncv
        ldv    = maxn
        
        if(allocated(v)) deallocate(v)
        if(allocated(d)) deallocate(d)
        if(allocated(ax)) deallocate(ax)
        if(allocated(resid)) deallocate(resid)
        if(allocated(workd)) deallocate(workd)
        if(allocated(workl)) deallocate(workl)
        if(allocated(selector)) deallocate(selector)
        allocate(ax(maxn),d(maxncv,2),resid(maxn),selector(maxncv),v(ldv,maxncv),workl(maxncv*(maxncv+8)),workd(3*maxn))

        !SET ARPACK PARAMETERS: Make sure that 1) maxn >= n 2) maxnev >= nev 3) maxncv >= ncv
        !Set dimensions for this problem.
        n     = int(dim, 4)
        bmat  = 'I'
        which = mode

        if ( n .gt. maxn ) then
            print *, ' ERROR: N is greater than MAXN '
            go to 9870
        else if ( nev .gt. maxnev ) then
            print *, ' ERROR: NEV is greater than MAXNEV '
            go to 9870
        else if ( ncv .gt. maxncv ) then
            print *, ' ERROR: NCV is greater than MAXNCV '
            go to 9870
        end if

        lworkl = ncv * ( ncv + 8 )
        tol    = zero
        info   = 0
        ido    = 0
        ishfts = 1
        maxitr = 300! 5000! 1000 !300
        mode1  = 1
        iparam = 0

        iparam(1) = ishfts
        iparam(3) = maxitr
        iparam(7) = mode1
        iparam(4) = 1 
        !MAIN LOOP (Reverse communication loop)! Repeatedly call DSAUPD and take actions inDicated by parameter !IDO until convergence is inDicated or MAXITR is exceeded.
        ierr     = 0
        cntr     = 0
        ax       = 0.d0
        v        = 0.d0
        workl    = 0.d0
        workd    = 0.d0
        sigma    = 0.d0
        d        = 0.d0
        resid    = 0.d0
        ipntr    = 0.d0
        selector = 0.d0 
       !!$omp critical
        do
           !Sparse eigensolver routine
           call dsaupd( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
           if ( ido /= -1 .and. ido /= 1 ) then
               exit
           end if
           cntr = cntr + 1
           if(mod(cntr, 100) == 0) print*,'cntr',cntr
           !Vector-sparse matrix multiplication                 
           call pamux2(threads, dim, sites, nbonds, nnnbonds, t1, v1, v2, basis, bsites, hexsites, workd(ipntr(1)), workd(ipntr(2)))
        end do
       !!$omp end critical
        print*,'Number of MV-multiplications',cntr
        call dsaupderrormessage(info)
        
        !Extract eigenvalues and -vectors
        ! !$omp critical
        call dseupd( rvec, 'A', selector, d, v, ldv, sigma, bmat, n, which, nev, tol, &
                     resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, ierr )
        call dseupderrormessage(ierr)
        ! !$omp end critical
        if (ierr /= 0) then
           print *, 'DSSIMP - Fatal error! Error with DSEUPD, IERR = ', ierr, 'Check the documentation of DSEUPD.'
        else
           nconv =  iparam(5)
           do j = 1, nconv
               call pamux2(threads, dim, sites, nbonds, nnnbonds, t1, v1, v2, basis, bsites, hexsites, v(1,j), ax)
               call daxpy(n, -d(j,1), v(1,j), 1, ax, 1 )
               d(j,2) = dnrm2(n, ax, 1)
               d(j,2) = d(j,2) / abs(d(j,1))
           end do
           call dmout(6, nconv, 2, d, maxncv, -6, 'Ritz values are named relative residuals' ) !Display computed residuals. 6: Output to screen Write(6, #internalnumber)! nconv: number of rows in the matrix d! 2: Number of columns in matrix d! maxncv: Leading dimension of the matrix data! -6: print the matrix d with iabs(-6) decimal digits per number! Use formatting indexed by -6 to print A
        end if

        evals(1:nev) = d(1:nev,1)

        if (rvec) then
            do j= 1, min(iparam(5), nest) 
                evecs(1:dim, j) = v(1:dim,j)
            end do
        end if
        
        9870 continue

    end if


    100 format(1000(F40.30))

    ! if(allocated(workd))    deallocate(workd)
    ! if(allocated(ax))       deallocate(ax)
    ! if(allocated(d))        deallocate(d)
    ! if(allocated(resid))    deallocate(resid)
    ! if(allocated(selector)) deallocate(selector)
    ! if(allocated(workl))    deallocate(workl)
    ! if(allocated(v))        deallocate(v)
    ! if(allocated(jao))      deallocate(jao)
    ! if(allocated(iao))      deallocate(iao)
    ! if(allocated(ao))       deallocate(ao)


    return 

end subroutine lanczos_d_otf

!-------------------------------------------------!
!            Sparse diagonalization  (MKL)        !
!-------------------------------------------------!

subroutine lanczos_mkl(dim, nev, ncv, nest, nnz, values, rc, evals, evecs)
    
    USE MKL_SPBLAS
    USE MKL_SOLVERS_EE    
    use iso_c_binDing, only: c_double, c_int
    implicit none
    
    !!!!!!!!!!!!!!!!! Matrix declaration variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: nev
    integer, intent(in) :: ncv
    integer, intent(in) :: nest
    integer, intent(in) :: nnz
    ! integer, intent(in) :: row_indx(nnz), col_indx(nnz)   
    integer, intent(in) :: rc(nnz,2)
    double precision, intent(in) ::  values(nnz)
    
    double precision, allocatable, intent(out) :: evals(:)
    double precision, allocatable, intent(out) :: evecs(:,:)

    integer :: i, threads 
    !   Matrix descriptor
    TYPE(MATRIX_DESCR) :: descrA
    !   CSR matrix structure
    TYPE(SPARSE_MATRIX_T) ::  cooA, csrA

    !!!!!!!!!!!!!!!!! Declaration of MKL_SPARSE_S_EV variables !!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! E - eigenvalues, X - eigenvectors, res - residual !!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! XL - left singular vectors, XR - right singular vectors !!!!!!

    character*1 WHICH
    parameter  (WHICH='S') !L - for largest eigenvalues to find
    integer, allocatable :: pm(:)
    integer :: n, K, info

    double precision, allocatable :: E(:)
    double precision, allocatable :: sigma(:)
    double precision, allocatable :: X(:,:)
    double precision, allocatable :: res(:)

   
    !!$omp threadprivate(i, descrA, cooA, csrA, pm, n, k, info, E, sigma, X, res)

    
    if(allocated(E)) deallocate(E)
    if(allocated(X)) deallocate(X)
    if(allocated(pm)) deallocate(pm)
    if(allocated(res)) deallocate(res)
    if(allocated(evals)) deallocate(evals)
    if(allocated(evecs)) deallocate(evecs)
    if(allocated(sigma)) deallocate(sigma)
    
    allocate(E(nev))
    allocate(pm(128))
    allocate(res(dim))
    allocate(sigma(dim))
    allocate(evals(nev))
    allocate(X(dim, nev))
    allocate(evecs(dim, nest))    
    

    evals = 0.d0 
    evecs = 0.d0 
    n = int(dim, 4) 
    print *, 'Sparse matrix size', n
    !
    !        Task 0. Call MKL_SPARSE_C_CREATE_COO to create matrix handle 
    !      

    ! info = mkl_sparse_d_create_coo (cooA, SPARSE_INDEX_BASE_ZERO, n, n, nnz, row_indx, col_indx, values)
    info = mkl_sparse_d_create_coo (cooA, SPARSE_INDEX_BASE_ONE, n, n, nnz, rc(1:nnz,1), rc(1:nnz,2), values)

    !
    !        Task 1. Call MKL_SPARSE_C_CREATE_CSR to create matrix handle
    !      

    info = mkl_sparse_convert_csr (cooA, SPARSE_OPERATION_NON_TRANSPOSE, csrA)
    info = MKL_SPARSE_DESTROY(cooA)
    ! info = mkl_sparse_d_create_csr(csrA,SPARSE_INDEX_BASE_ZERO,N,N,rows,rows(2),cols,val)

    !         Create matrix descriptor
    descrA % TYPE = SPARSE_MATRIX_TYPE_GENERAL
    !
    !        Task 2. Call MKL_SPARSE_EE_INIT to define the default values for the input
    !        parameters.
    !
    info = mkl_sparse_ee_init(pm)
    pm(2) = 6 !Setting tolerance 10**(-pm(2) + 1)  
    pm(3) = 1
    !0 - Decided at runtime
    !1 - The Krylov-Schur method
    !2 - Subspace Iteration technique based on FEAST algorithm (for eigenvalues in specific interval)
    pm(4) = ncv 
    ! This parameter is referenced only for Krylov-Schur Method. It inDicates the number of Lanczos/Arnoldi vectors (NCV) generated at each iteration.
    ! This parameter must be less than or equal to size of matrix and greater than number of eigenvalues (k0) to be computed. If unspecified, NCV is set to be at least 1.5 times larger than k0. 

    pm(8) = 1 ! Use absolute stopping criteria
    !
    !         Task 3. Solve the standard eigenvalue problem Ax=ex.
    ! 
    !

    !!$ threads = omp_get_num_threads()
    ! print*,'threads',threads
    !!$ threads = omp_get_thread_num()
    ! print*,'id',threads
    !!$omp critical 

    ! print* ,x, 'x'
    ! print* ,e, 'e'
    ! print* ,res, 'res'
    ! print* ,nest, 'nest'
    ! print* ,nev, 'nev'
    ! print* ,dim, 'dim'
    ! print* ,pm, 'pm'
    ! print* ,k, 'k'
    
    info = mkl_sparse_d_ev(WHICH,pm,csrA,descrA,nev,k,E,X,res)
    !!$omp end critical 
    print  *,' OUTPUT INFO ',info
    print  *,' ITERATION EXIT REASON (0 = converged) ', pm(10)
    
    if(info.ne.0) stop 1
    if(pm(10).ne.0) stop 2
    print *, 'Number of eigenvalues found ', k
    ! print *, ' Computed    |    Expected  '
    ! print *, ' Eigenvalues |    Eigenvalues '
    do i = 1, K
        ! write(*,'(f20.16)') E(i)
        ! print *, E(i)
        evals(i) = E(i)
    end do
    !   Release internal representation of CSR matrix
    info = MKL_SPARSE_DESTROY(csrA)

    ! print*, 'Number of eigenstates to be saved ', nest 
    do i = 1, nest 
        ! evecs(1:n,i) = X(i,1:n) 
        evecs(1:n,i) = X(1:n,i) 
    end do 
    
    return 

end subroutine lanczos_mkl

subroutine lanczos_c( threads, dim, nev, ncv, nst, mode, rvec, nnz, ham, rc, evals, evecs )
    
    implicit none

    integer, intent(in) :: threads, nev, ncv, nst, nnz
    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: rc(nnz,2)
    double complex, intent(in) :: ham(nnz)   
    character*2, intent(in) :: mode
    logical, intent(in) :: rvec
    double precision, allocatable, intent(out) :: evals(:)
    double complex, allocatable, intent(out) :: evecs(:, :)

    integer :: printing = 1 
    integer(kind=8), save :: maxn = 0, maxnev = 0, maxncv = 0, ldv = 0
    integer(kind=8), save :: n = 0, nx = 0
    integer, save :: buffer = 0 
    integer, save :: j = 0, iparam(11), ipntr(14)
    integer, save :: ido = 0, ishfts = 0, lworkl = 0, info = 0, maxitr = 0, mode1 = 0, nconv = 0, ierr = 0    !VARIABLES FOR DIAGONALIZATION ROUTINE
    integer, allocatable, save :: jao(:), iao(:), jao_ord(:)
    double precision, save :: tol = 0
    double precision, allocatable, save :: rwork(:), rd(:,:)
    double complex, save :: sigma = 0    
    double complex, allocatable, save :: ao(:), ao_ord(:)
    double complex, allocatable, save :: ax(:), d(:), v(:,:), workd(:), workev(:), resid(:), workl(:)!VARIABLES FOR DIAGONALIZATION ROUTINE
    character, save :: bmat*1, which*2
    logical, allocatable, save :: selector(:)

    intrinsic :: abs
    !c
    !c     %-----------------------------%
    !c     | BLAS & LAPACK routines used |
    !c     %-----------------------------%
    !c

    double precision :: dznrm2 , dlapy2
    external :: dznrm2 , zaxpy , dlapy2, znaupd, zneupd, dmout
    double precision, external :: dnrm2

    ! c  WHICH   Character*2.  (INPUT)
    !c
    !c     %-----------------------%
    !c     | Executable Statements |
    !c     %-----------------------%
    !c
    !c     %--------------------------------------------------%
    !c     | The number NX is the number of interior points   |
    !c     | in the discretization of the 2-dimensional       |
    !c     | convection-diffusion operator on the unit        |
    !c     | square with zero Dirichlet boundary conDition.   |
    !c     | The number N(=NX*NX) is the dimension of the     |
    !c     | matrix.  A standard eigenvalue problem is        |
    !c     | solved (BMAT = 'I').  NEV is the number of       |
    !c     | eigenvalues to be approximated.  The user can    |
    !c     | modify NX, NEV, NCV, WHICH to solve problems of  |
    !c     | different sizes, and to get different parts of   |
    !c     | the spectrum.  However, The following            |
    !c     | conDitions must be satisfied:                    |
    !c     |                   N <= MAXN                      |
    !c     |                 NEV <= MAXNEV                    |
    !c     |           NEV + 2 <= NCV <= MAXNCV               |
    !c     %--------------------------------------------------%
    !c

    !$omp threadprivate(ao, jao, iao, bmat, which, j, iparam, sigma, ipntr, ido, ishfts, lworkl, info, maxitr, mode1, nconv, ierr, selector, ax, d, v, workd, workev, resid, workl, rwork, rd)
    
    
    if (dim == 1) then
        evecs(1,1) = 1
        evals(1) = ham(1)
        print*, ''
        write(*,"('Only eigenvalue is ',f8.4)") ham(1)
    else
    
    buffer = 50 
    ! maxn = 10 + dim * dim !debug
    maxn   = buffer + dim !* dim !debug
    maxnev = buffer + nev
    maxncv = buffer + ncv
    ldv    = maxn

    if(allocated(evals))    deallocate(evals)
    if(allocated(evecs))    deallocate(evecs)
    if(allocated(workd))    deallocate(workd)
    if(allocated(workev))   deallocate(workev)
    if(allocated(rwork))    deallocate(rwork)
    if(allocated(rd))       deallocate(rd)
    if(allocated(ax))       deallocate(ax)
    if(allocated(d))        deallocate(d)
    if(allocated(resid))    deallocate(resid)
    if(allocated(selector)) deallocate(selector)
    if(allocated(workl))    deallocate(workl)
    if(allocated(v))        deallocate(v)
    allocate(evals(nev), evecs(dim, nst))
    allocate(ax(maxn), d(maxncv), resid(maxn), selector(maxncv), v(ldv,maxncv), &
            workl(3*maxncv*maxncv+5*maxncv), workd(3*maxn), workev(3*maxncv), &
            rwork(maxncv), rd(maxncv,3))
    ax     = 0.d0 
    d      = 0.d0 
    resid  = 0.d0 
    v      = 0.d0
    workl  = 0.d0
    workd  = 0.d0 
    workev = 0.d0
    rwork  = 0.d0 
    rd     = 0.d0 


    nx = dim !debug
    ! n  = nx * nx !debug
    n  = nx !* nx !debug

    if ( n .gt. maxn ) then
        print *, ' ERROR with _NDRV1: N is greater than MAXN '
        go to 9870
    else if ( nev .gt. maxnev ) then
        print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
        go to 9870
    else if ( ncv .gt. maxncv ) then
        print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
        go to 9870
    end if
    bmat  = 'I'
    which = 'SR'
    !c
    !c     %---------------------------------------------------%
    !c     | The work array WORKL is used in ZNAUPD  as         |
    !c     | workspace.  Its dimension LWORKL is set as        |
    !c     | illustrated below.  The parameter TOL determines  |
    !c     | the stopping criterion. If TOL<=0, machine        |
    !c     | precision is used.  The variable IDO is used for  |
    !c     | reverse communication, and is initially set to 0. |
    !c     | Setting INFO=0 inDicates that a random vector is  |
    !c     | generated to start the ARNOLDI iteration.         |
    !c     %---------------------------------------------------%
    !c
    lworkl = 3*ncv**2+5*ncv
    tol    = 0
    ido    = 0
    info   = 0
    nconv  = 0


    !c
    !c     %---------------------------------------------------%
    !c     | This program uses exact shift with respect to     |
    !c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
    !c     | IPARAM(3) specifies the maximum number of Arnoldi |
    !c     | iterations allowed.  Mode 1 of ZNAUPD  is used     |
    !c     | (IPARAM(7) = 1). All these options can be changed |
    !c     | by the user. For details see the documentation in |
    !c     | ZNAUPD .                                           |
    !c     %---------------------------------------------------%
    !c
    ishfts = 1
    maxitr = 300 !5000
    mode1  = 1

    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(7) = mode1

    

    !c
    !c Create sparse matrix in CSR format
    !c
    
    if(allocated(jao)) deallocate(jao)
    if(allocated(iao)) deallocate(iao)
    if(allocated(ao))  deallocate(ao)
    allocate(jao(nnz))
    allocate(iao(dim+1))
    allocate(ao(nnz))
    ao  = 0.d0
    iao = 0 
    jao = 0 
    
    call ccoocsr(dim, nnz, ham, rc(1:nnz,1), rc(1:nnz,2), ao, jao, iao)

    ! if(allocated(jao_ord)) deallocate(jao_ord)
    ! if(allocated(ao_ord)) deallocate(ao_ord)
    ! allocate(jao_ord(nnz))
    ! allocate(ao_ord(nnz))
    ! ao_ord  = 0.d0 
    ! jao_ord = 0
    ! call sortCols(dim, nnz, iao, jao, ao, jao_ord, ao_ord)
    ! ao = ao_ord 
    ! jao = jao_ord 
    !$omp critical
    do
        !Sparse eigensolver routine
        call znaupd ( ido, bmat, n, which, nev, tol, resid, ncv, &
                      v, ldv, iparam, ipntr, workd, workl, lworkl, &
                      rwork, info )
        if ( ido /= -1 .and. ido /= 1 ) then
            exit
        end if
        !Vector-sparse matrix multiplication
        ! call mkl_zcsrmv(transa, m, k, alpha, matdescra, val, indx, pntrb, pntre, x, beta, y)
        call cpamux (threads, nx, workd(ipntr(1)), workd(ipntr(2)), ao, jao, iao)
    end do
    !$omp end critical
 
    ! !c
    ! !c     %-------------------------------------------%
    ! !c     | M A I N   L O O P (Reverse communication) |
    ! !c     %-------------------------------------------%
    ! !c
    ! 10   continue
    ! !c
    ! !c        %---------------------------------------------%
    ! !c        | Repeatedly call the routine ZNAUPD  and take |
    ! !c        | actions inDicated by parameter IDO until    |
    ! !c        | either convergence is inDicated or maxitr   |
    ! !c        | has been exceeded.                          |
    ! !c        %---------------------------------------------%
    ! !c
    ! call znaupd ( ido, bmat, n, which, nev, tol, resid, ncv, &
    !             v, ldv, iparam, ipntr, workd, workl, lworkl, &
    !             rwork, info )
    ! !c
    ! !c           %-------------------------------------------%
    ! !c           | Perform matrix vector multiplication      |
    ! !c           |                y <--- OP*x                |
    ! !c           | The user should supply his/her own        |
    ! !c           | matrix vector multiplication routine here |
    ! !c           | that takes workd(ipntr(1)) as the input   |
    ! !c           | vector, and return the matrix vector      |
    ! !c           | product to workd(ipntr(2)).               |
    ! !c           %-------------------------------------------%
    ! !c

    !     call cpamux (threads, nx, workd(ipntr(1)), workd(ipntr(2)), ao, jao, iao)

    ! !c
    ! !c           %-----------------------------------------%
    ! !c           | L O O P   B A C K to call ZNAUPD  again. |
    ! !c           %-----------------------------------------%
    ! !c
    !     go to 10

    ! end if


    !c
    !c     %----------------------------------------%
    !c     | Either we have convergence or there is |
    !c     | an error.                              |
    !c     %----------------------------------------%
    !c
    if ( info .lt. 0 ) then
        !c
        !c        %--------------------------%
        !c        | Error message, check the |
        !c        | documentation in ZNAUPD   |
        !c        %--------------------------%
        !c
        print *, ' '
        print *, ' Error with znaupd, info = ', info
        print *, ' Check the documentation of _naupd'
        print *, ' '

    else
    !c
    !c        %-------------------------------------------%
    !c        | No fatal errors occurred.                 |
    !c        | Post-Process using ZNEUPD .                |
    !c        |                                           |
    !c        | Computed eigenvalues may be extracted.    |
    !c        |                                           |
    !c        | Eigenvectors may also be computed now if  |
    !c        | desired.  (inDicated by rvec = .true.)    |
    !c        %-------------------------------------------%
    !c


        call zneupd (rvec, 'A', selector, d, v, ldv, sigma, &
                    workev, bmat, n, which, nev, tol, resid, ncv, &
                    v, ldv, iparam, ipntr, workd, workl, lworkl, &
                    rwork, ierr)
    !c
    !c        %----------------------------------------------%
    !c        | Eigenvalues are returned in the one          |
    !c        | dimensional array D.  The corresponDing      |
    !c        | eigenvectors are returned in the first NCONV |
    !c        | (=IPARAM(5)) columns of the two dimensional  |
    !c        | array V if requested.  Otherwise, an         |
    !c        | orthogonal basis for the invariant subspace  |
    !c        | corresponDing to the eigenvalues in D is     |
    !c        | returned in V.                               |
    !c        %----------------------------------------------%
    !c


        if ( ierr .ne. 0) then
            !c
            !c           %------------------------------------%
            !c           | Error conDition:                   |
            !c           | Check the documentation of ZNEUPD  |
            !c           %------------------------------------%
            !c
            print *, ' '
            print *, ' Error with zneupd, info = ', ierr
            print *, ' Check the documentation of _neupd. '
            print *, ' '

        else

            nconv = iparam(5)
            do 20 j = 1, nconv
                !c
                !c               %---------------------------%
                !c               | Compute the residual norm |
                !c               |                           |
                !c               |   ||  A*x - lambda*x ||   |
                !c               |                           |
                !c               | for the NCONV accurately  |
                !c               | computed eigenvalues and  |
                !c               | eigenvectors.  (iparam(5) |
                !c               | inDicates how many are    |
                !c               | accurate to the requested |
                !c               | tolerance)                |
                !c               %---------------------------%
                !c

                call cpamux (threads, nx, v(1,j), ax, ao, jao, iao)
                call zaxpy (n, -d(j), v(1,j), 1, ax, 1)
                rd(j,1) = dble (d(j))
                rd(j,2) = dimag (d(j))
                rd(j,3) = dznrm2 (n, ax, 1)
                rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
            20 continue
            !c
            !c            %-----------------------------%
            !c            | Display computed residuals. |
            !c            %-----------------------------%
            !c
            if(printing == 1) then 
                call dmout (6, nconv, 3, rd, maxncv, -6, &
                            'Ritz values (Real, Imag) and relative residuals')
            end if 
        end if


            ! c
            ! c        %-------------------------------------------%
            ! c        | Print additional convergence information. |
            ! c        %-------------------------------------------%
            ! c
        if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
        else if ( info .eq. 3) then
            print *, ' '
            print *, ' No shifts could be applied during implicit', &
                     ' Arnoldi update, try increasing NCV.'
            print *, ' '
        end if

        ! print *, ' '
        ! print *, '_NDRV1'
        ! print *, '====== '
        ! print *, ' '
        ! print *, ' Size of the matrix is ', n
        ! print *, ' The number of Ritz values requested is ', nev
        ! print *, ' The number of Arnoldi vectors generated', &
        !          ' (NCV) is ', ncv
        ! print *, ' What portion of the spectrum: ', which
        ! print *, ' The number of converged Ritz values is ', &
        !            nconv
        ! print *, ' The number of Implicit Arnoldi update', &
        !          ' iterations taken is ', iparam(3)
        ! print *, ' The number of OP*x is ', iparam(9)
        ! print *, ' The convergence criterion is ', tol
        ! print *, ' '

    end if

    evals(1:nev) = dble(d(1:nev))
    evecs = 0 
    
    if (rvec) then
        do j = 1, nst
            evecs(1:dim, j) = v(1:dim,j)
        end do
    end if

    
    end if


    !c
    !c     %---------------------------%
    !c     | Done with program zndrv1 . |
    !c     %---------------------------%
    !c

    !c          Error flag for ZNAUPD on output.
    !c          =  0: Normal exit.
    !c          =  1: Maximum number of iterations taken.
    !c                All possible eigenvalues of OP has been found. IPARAM(5)
    !c                returns the number of wanted converged Ritz values.
    !c          =  2: No longer an informational error. Deprecated starting
    !c                with release 2 of ARPACK.
    !c          =  3: No shifts could be applied during a cycle of the
    !c                Implicitly restarted Arnoldi iteration. One possibility
    !c                is to increase the size of NCV relative to NEV.
    !c                See remark 4 below.
    !c          = -1: N must be positive.
    !c          = -2: NEV must be positive.
    !c          = -3: NCV-NEV >= 2 and less than or equal to N.
    !c          = -4: The maximum number of Arnoldi update iteration
    !c                must be greater than zero.
    !c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
    !c          = -6: BMAT must be one of 'I' or 'G'.
    !c          = -7: Length of private work array is not sufficient.
    !c          = -8: Error return from LAPACK eigenvalue calculation;
    !c          = -9: Starting vector is zero.
    !c          = -10: IPARAM(7) must be 1,2,3.
    !c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
    !c          = -12: IPARAM(1) must be equal to 0 or 1.
    !c          = -9999: Could not build an Arnoldi factorization.
    !c                   User input error highly likely.  Please
    !c                   check actual array dimensions and layout.
    !c                   IPARAM(5) returns the size of the current Arnoldi
    !c                   factorization.


    !    Error flag for ZNEUPD on output.
    !    0: Normal exit.
    !    1: The Schur form computed by LAPACK routine csheqr could not be reordered by LAPACK routine ztrsen.
    !       Re-enter subroutine zneupd with IPARAM(5) = NCV and increase the size of the array D to have dimension at least dimension NCV and allocate at least NCV columns for Z.
    !       NOTE: Not necessary if Z and V share the same space. Please notify the authors if this error occurs.
    !    -1: N must be positive.
    !    -2: NEV must be positive.
    !    -3: NCV-NEV >= 1 and less than or equal to N.
    !    -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'.
    !    -6: BMAT must be one of 'I' or 'G'.
    !    -7: Length of private work WORKL array is not sufficient.
    !    -8: Error return from LAPACK eigenvalue calculation. This should never happened.
    !    -9: Error return from calculation of eigenvectors. Informational error from LAPACK routine ztrevc.
    !    -10: IPARAM(7) must be 1, 2, 3.
    !    -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
    !    -12: HOWMANY = 'S' not yet implemented.
    !    -13: HOWMANY must be one of 'A' or 'P' if RVEC = .true.
    !    -14: ZNAUPD did not find any eigenvalues to sufficient accuracy.
    !    -15: ZNEUPD got a different count of the number of converged Ritz values than ZNAUPD got. This inDicates the user probably made an error in passing data from ZNAUPD to ZNEUPD or that the data was modified before entering ZNEUPD.

    9870 continue


    if(allocated(workd)) deallocate(workd)
    if(allocated(workev)) deallocate(workev)
    if(allocated(rwork)) deallocate(rwork)
    if(allocated(rd)) deallocate(rd)
    if(allocated(ax)) deallocate(ax)
    if(allocated(d)) deallocate(d)
    if(allocated(resid)) deallocate(resid)
    if(allocated(selector)) deallocate(selector)
    if(allocated(workl)) deallocate(workl)
    if(allocated(v)) deallocate(v)
    if(allocated(jao)) deallocate(jao)
    if(allocated(iao)) deallocate(iao)
    if(allocated(ao)) deallocate(ao)


end subroutine lanczos_c

!------------------------------!
!            FEAST             !
!------------------------------!

subroutine cfeast(n, nnz, nev0, nest, rc, ham, Emin, Emax, nev, eigvals, eigvecs)
   
    !       1. The code calls  FEASTINIT  to define the default values for the input
    !          FEAST parameters.
    
    !       2. The  code  solves  the  standard eigenvalue  problem  Ax=ex   using 
    !          ZFEAST_HCSREV.
    
    !       3. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i) 
    !           are the expected eigenvalues  and E(i) are eigenvalues computed 
    !           with the help of ZFEAST_HCSREV().
    
    !       4. The code computes the maximum absolute value of the matrix  Y=(X')*X-I  
    !          where X is the matrix of eigenvectors computed with the help of 
    !          ZFEAST_HCSREV. ZGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
    ! *******************************************************************************
    USE MKL_SPBLAS   
    USE ISO_C_BINDING 
    implicit none
    !!!!!!!!!!!!!!!!! External function declarations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    external zgemm, zfeast_hcsrev, zfeast_hcsrgv, feastinit
    !!!!!!!!!!!!!!!!! Matrix declaration variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(kind=8), intent(in) :: n
    integer, intent(in) :: nnz
    integer, intent(in) :: nev0
    integer, intent(in) :: rc(nnz,2)
    double complex, intent(in) :: ham(nnz)
    double precision, intent(in) :: Emin, Emax
    integer, intent(inout) :: nest
    integer, intent(out) :: nev 
    double precision, allocatable, intent(out) :: eigvals(:)
    double complex, allocatable, intent(out) :: eigvecs(:,:) 
    
    TYPE(MATRIX_DESCR) :: descrA     ! Sparse matrix descriptor
    !   CSR matrix representation
    TYPE(SPARSE_MATRIX_T) :: cooA, csrA    ! Structure with sparse matrix

    ! integer :: nrows, ncols 
    ! integer, allocatable :: csr_row_ptr(:), csr_col_inDices(:)
    ! complex*16, allocatable :: csr_values(:)

    integer :: rows(n+1), cols(nnz)
    integer :: rows_b(n+1), cols_b(n)
    INTEGER(C_INT) :: indexing
    ! TYPE(C_PTR) :: rows_start_c, rows_end_c
    ! TYPE(C_PTR) :: val_c, col_c
    ! INTEGER, POINTER :: rows_start(:), rows_end(:)
    double complex :: val(nnz), val_b(n)
    double complex :: beta
    character*1 :: UPLO
    parameter   (UPLO='U')
    !!!!!!!!!!!!!!!!! Feast declaration variable !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! E - eigenvalues, X - eigenvectors, res - residual !!!!!!!!!!!!
    integer :: fpm(128)  
    double precision :: epsout
    integer :: loop
    integer :: L
    parameter   (L=10)
    integer :: M0,M,info
    double precision :: E(nev0)
    ! double precision :: E(n)
    !On output, the first m columns of x contain the orthonormal eigenvectors corresponDing 
    !to the computed eigenvalues e, with the i-th column of x holding the eigenvector associated with e(i).   
    double complex :: X(n,nev0) 
    ! double complex :: X(n,n) 
    ! double precision :: res(n)
    double precision :: res(nev0)
    !!!!!!!!!!!!!!!!! Declaration of local variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! Eig - exact eigenvalues, R=|E-Eig|, Y=(X')*X-I !!!!!!!!!!!!!!!
    double precision :: Eig(nev0)
    double precision :: R(nev0)
    ! double precision      Eig(n)
    ! double precision      R(n)
    ! double complex  Y(n,n)
    integer :: i,j,n4
    integer ::  mkl, dim 
    integer ::  UTnnz = 0 
    integer ::  ldx, ldy
    double complex :: one, zero
    double precision :: smax, eigabs 

    integer, allocatable :: ib(:), jb(:), ic(:), jc(:)
    double complex, allocatable :: b(:), c(:)
    integer, allocatable :: UTcols(:), UTrows(:)
    double complex, allocatable :: UTval(:)

    print*,'Sparse matrix size',n
    print*,'Use Feast algorithm'
    M0=nev0
    M=int(0.5*M0)
    print *,'Search interval ', Emin,' ', Emax
    ldx=n
    ldy=n

    ! mkl = mkl_sparse_z_create_coo (cooA, SPARSE_INDEX_BASE_ONE, n, n, nnz, rc(1:nnz,1), rc(1:nnz,2), ham)
    ! mkl = mkl_sparse_convert_csr (cooA, SPARSE_OPERATION_NON_TRANSPOSE, csrA)
    ! mkl = mkl_sparse_z_export_csr (csrA, indexing, dim, dim, rows_start_c, rows_end_c, col_c, val_c)
    ! cols = rc(1:nnz,2) 
    ! val  = ham 
    ! call C_F_POINTER(rows_start_c, rows_start, [n])
    ! call C_F_POINTER(rows_end_c  , rows_end  , [n])
    ! rows(1:n) = rows_start(1:n)
    ! rows(n+1) = nnz + 1 
    if(allocated(b)) deallocate(b)
    if(allocated(jb)) deallocate(jb)
    if(allocated(ib)) deallocate(ib)
    allocate(b(n))
    allocate(jb(n))
    allocate(ib(n))
    val = (0.d0,0.d0)
    val_b = (0.d0,0.d0)
    cols = 0 
    rows = 0 
    cols_b = 0 
    rows_b = 0 

    do I = 1, N
        b(i)  = 1.d0
        ib(i) = i 
        jb(i) = i
    end do
    ! open(20,file='A.dat')
    !     do i = 1, nnz 
    ! write(20,'(f10.6,x,f10.8,x,i0,x,i0)') dble(ham(i)), aimag(ham(i)), rc(i,1), rc(i,2)
    !     end do 
    ! close(20)
    ! open(20,file='B.dat')
    ! do i = 1, n 
    ! write(20,'(f10.6,x,f10.8,x,i0,x,i0)') dble(b(i)), aimag(b(i)), ib(i), jb(i)
    !     end do 
    ! close(20)
    n4 = n 
    ! call coo_to_csr_hermitian(nnz, ham, rc(1:nnz,1), rc(1:nnz,2), csr_values, csr_row_ptr, csr_col_inDices, nrows)
    call ccoocsr(n, nnz, ham, rc(1:nnz,1), rc(1:nnz,2), val, cols, rows)
    call ccoocsr(n, n4, b, ib, jb, val_b, cols_b, rows_b)

    deallocate(b, ib, jb)
    if(allocated(c)) deallocate(c)
    if(allocated(jc)) deallocate(jc)
    if(allocated(ic)) deallocate(ic)
    allocate(ic(n+1))
    ic = 0 
    allocate(c(nnz))
    allocate(jc(nnz))
    jc = 0 
    c  = 0.d0 
    
    ! call mkl_zcsradd('N', 0, 1, n, n, ham, cols, rows, beta, b, jb, ib, c, jc, ic, 1, info)
    call mkl_zcsrmultcsr('N', 0, 3, n, n, n, val, cols, rows, val_b, cols_b, rows_b, c, jc, ic, nnz, info)
    ! call mkl_zcsrmultcsr('N', 0, 1, n, n, n, csr_values, csr_col_inDices, csr_row_ptr, val_b, cols_b, rows_b, c, jc, ic, nnz, info)
    deallocate(c, ic, jc)

    call HermitianUpperTriangularCSR(n, nnz, val, rows, cols, UTval, UTrows, UTcols, UTnnz)
        
    ! call mkl_zcsrmultcsr('T', 1, 0, n4, n4, n4, val, cols, rows, val, cols, rows, c, jc, ic, nnz, info)
    ! allocate(c(ic(n+1)-1))
    ! allocate(jc(ic(n+1)-1))
    ! jc = 0 
    ! c  = 0.d0 
    ! call mkl_zcsrmultcsr('T', 2, 0, n4, n4, n4, val, cols, rows, val, cols, rows, c, jc, ic, nnz, info)
    
    !
    !        Task 1. Call  FEASTINIT  to define the default values for the input
    !        FEAST parameters.
    !
    call feastinit(fpm)
    fpm(1)=1  !Specifies whether Extended Eigensolver routines print runtime status. 
    fpm(26)=1 !Specifies whether Extended Eigensolver routines check input matrices (applies to CSR format only). 
    fpm(27)=1 !Specifies whether Extended Eigensolver routines check input matrices (applies to CSR format only). 

    print *, ' Testing zfeast_hcsrev '
    !
    !         Task 2. Solve the standard eigenvalue problem Ax=ex.
    !
    ! call zfeast_hcsrev(UPLO,N,csr_values,csr_row_ptr,csr_col_inDices,fpm,epsout,loop,&
    !    Emin,Emax,M0,E,X,M,res,info)

    call zfeast_hcsrev(UPLO,N,UTval,UTrows,UTcols,fpm,epsout,loop,&
       Emin,Emax,M0,E,X,M,res,info)
    ! call zfeast_hcsrev(UPLO,N,val,rows,cols,fpm,epsout,loop,&
    !    Emin,Emax,M0,E,X,M,res,info)
    print  *,' FEAST OUTPUT INFO ',info
    if(info.ne.0) stop 1
    nev = m 
    nest = min(nest, nev)
    if(allocated(eigvals)) deallocate(eigvals)
    if(allocated(eigvecs)) deallocate(eigvecs)
    allocate(eigvals(nev))
    allocate(eigvecs(n, nest))
    eigvals = 0.d0 
    eigvecs = 0.d0 
    do i = 1, nev
        eigvals(i) = E(i)        
    end do
    do i = 1, nest
        eigvecs(1:n, i) = X(1:n, i)
    end do 
    !
    !         Task 3. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
    !         are the expected eigenvalues  and E(i) are eigenvalues computed
    !         with the help of ZFEAST_HCSREV().
    !
    print *, 'Number of eigenvalues found ', M
    print *, ' Computed    |    Expected  '
    print *, ' Eigenvalues |    Eigenvalues '
    eigabs=zero
    do i=1,M
       R(i)=dabs(E(i)-Eig(i))
       eigabs=max(eigabs, R(i))
       print *, E(i), Eig(i)
    enddo
    print *, ' Max value of '
    print *, ' | computed eigenvalue -expected eigenvalues | ', eigabs
  
    
    return 

end subroutine cfeast

subroutine dfeast(n, nnz, nev0, nest, rc, ham, Emin, Emax, nev, eigvals, eigvecs)
   
    !       1. The code calls  FEASTINIT  to define the default values for the input
    !          FEAST parameters.
    
    !       2. The  code  solves  the  standard eigenvalue  problem  Ax=ex   using 
    !          ZFEAST_HCSREV.
    
    !       3. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i) 
    !           are the expected eigenvalues  and E(i) are eigenvalues computed 
    !           with the help of ZFEAST_HCSREV().
    
    !       4. The code computes the maximum absolute value of the matrix  Y=(X')*X-I  
    !          where X is the matrix of eigenvectors computed with the help of 
    !          ZFEAST_HCSREV. ZGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
    ! *******************************************************************************
    USE MKL_SPBLAS   
    USE ISO_C_BINDING 
    implicit none
    !!!!!!!!!!!!!!!!! External function declarations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    external zgemm, zfeast_hcsrev, zfeast_hcsrgv, feastinit
    !!!!!!!!!!!!!!!!! Matrix declaration variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(kind=8), intent(in) :: n
    integer, intent(in) :: nnz
    integer, intent(in) :: nev0
    integer, intent(in) :: rc(nnz,2)
    double precision, intent(in) :: ham(nnz)
    double precision, intent(in) :: Emin,Emax
    integer, intent(inout) :: nest
    integer, intent(out) :: nev 
    double precision, allocatable, intent(out) :: eigvals(:)
    double precision, allocatable, intent(out) :: eigvecs(:,:) 
    
    TYPE(MATRIX_DESCR) :: descrA     ! Sparse matrix descriptor
    !   CSR matrix representation
    TYPE(SPARSE_MATRIX_T) :: cooA, csrA    ! Structure with sparse matrix

    integer :: rows(n+1), cols(nnz)
    integer :: rows_b(n+1), cols_b(n)
    INTEGER(C_INT) :: indexing
    double precision :: val(nnz), val_b(n)
    double precision :: beta
    character*1 :: UPLO
    parameter   (UPLO='F')
    !!!!!!!!!!!!!!!!! Feast declaration variable !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! E - eigenvalues, X - eigenvectors, res - residual !!!!!!!!!!!!
    integer :: fpm(128)  
    double precision :: epsout
    integer :: loop
    integer :: L
    parameter   (L=10)
    integer :: M0, M, info
    double precision :: E(nev0)
    !On output, the first m columns of x contain the orthonormal eigenvectors corresponDing 
    !to the computed eigenvalues e, with the i-th column of x holding the eigenvector associated with e(i).   
    double precision :: X(n,nev0)
    double precision :: res(nev0)
    !!!!!!!!!!!!!!!!! Declaration of local variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! Eig - exact eigenvalues, R=|E-Eig|, Y=(X')*X-I !!!!!!!!!!!!!!!
    double precision :: Eig(nev0)
    double precision :: R(nev0)
    integer :: i,j,n4
    integer :: mkl, dim 
    integer :: UTnnz = 0 
    integer :: ldx, ldy
    double precision :: one, zero
    double precision :: smax, eigabs 

    integer, allocatable :: ib(:), jb(:), ic(:), jc(:)
    double precision, allocatable :: b(:), c(:)
    integer, allocatable :: UTcols(:), UTrows(:)
    double precision, allocatable :: UTval(:)

    print*,'Sparse matrix size',n
    print*,'Use Feast algorithm'

    M0=nev0
    M=10
    ! M=int(0.5*M0)
    print *,'Search interval ', Emin,' ', Emax
    ldx = n
    ldy = n

    if(allocated(b)) deallocate(b)
    if(allocated(jb)) deallocate(jb)
    if(allocated(ib)) deallocate(ib)
    allocate(b(n))
    allocate(jb(n))
    allocate(ib(n))
    val    = 0.d0
    val_b  = 0.d0
    cols   = 0 
    rows   = 0 
    cols_b = 0 
    rows_b = 0 

    do I = 1, N
        b(i)  = 1.d0
        ib(i) = i 
        jb(i) = i
    end do

    n4 = n 

    call coocsr(n, nnz, ham, rc(1:nnz,1), rc(1:nnz,2), val, cols, rows)
    call coocsr(n, n4, b, ib, jb, val_b, cols_b, rows_b)

    deallocate(b, ib, jb)
    if(allocated(c)) deallocate(c)
    if(allocated(jc)) deallocate(jc)
    if(allocated(ic)) deallocate(ic)
    allocate(ic(n+1))
    ic = 0 
    allocate(c(nnz))
    allocate(jc(nnz))
    jc = 0 
    c  = 0.d0

    call mkl_dcsrmultcsr('N', 0, 3, n, n, n, val, cols, rows, val_b, cols_b, rows_b, c, jc, ic, nnz, info)

    deallocate(c, ic, jc)

    ! call HermitianUpperTriangularCSR(n, nnz, val, rows, cols, UTval, UTrows, UTcols, UTnnz)
  
    !
    !        Task 1. Call  FEASTINIT  to define the default values for the input
    !        FEAST parameters.
    !
    call feastinit(fpm)
    fpm(1)  = 1 !Specifies whether Extended Eigensolver routines print runtime status. 
    fpm(26) = 1 !Specifies whether Extended Eigensolver routines check input matrices (applies to CSR format only). 
    fpm(27) = 1 !Specifies whether Extended Eigensolver routines check input matrices (applies to CSR format only). 

    print *, ' Testing dfeast_scsrev '
    !
    !         Task 2. Solve the standard eigenvalue problem Ax=ex.
    !
    
    call dfeast_scsrev(UPLO,N,val,rows,cols,fpm,epsout,loop,&
       Emin,Emax,M0,E,X,M,res,info)
    
    print*,' FEAST OUTPUT INFO ', info
    if(info.ne.0) stop 1
    nev  = m 
    nest = min(nest, nev)
    if(allocated(eigvals)) deallocate(eigvals)
    if(allocated(eigvecs)) deallocate(eigvecs)
    allocate(eigvals(nev))
    allocate(eigvecs(n,nest))
    eigvals = 0.d0 
    eigvecs = 0.d0
    do i = 1, nev
        eigvals(i) = E(i)        
    end do
    do i = 1, nest
        eigvecs(1:n, i) = X(1:n, i)
    end do 
    !
    !         Task 3. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
    !         are the expected eigenvalues  and E(i) are eigenvalues computed
    !         with the help of ZFEAST_HCSREV().
    !
    print*, 'Number of eigenvalues found ', M
    print*, ' Computed    |    Expected  '
    print*, ' Eigenvalues |    Eigenvalues '
    eigabs = zero
    do i = 1, M
       R(i)   = dabs(E(i)-Eig(i))
       eigabs = max(eigabs, R(i))
       print*, E(i), Eig(i)
    enddo
    print*, ' Max value of '
    print*, ' | computed eigenvalue -expected eigenvalues | ', eigabs

    return 

end subroutine dfeast

subroutine HermitianUpperTriangularCSR(n, nnz, A, IA, JA, UTriangular, IA_UTriangular, JA_UTriangular, nnz_UpperTriangular)
    implicit none
    integer(kind=8), intent(in) :: n
    integer, intent(in) :: nnz
    complex*16, intent(in) :: A(nnz)
    integer, intent(in) :: IA(n+1), JA(nnz)
    integer, intent(out) :: nnz_UpperTriangular
    complex*16, allocatable, intent(out) :: UTriangular(:)
    integer, allocatable, intent(out) :: IA_UTriangular(:), JA_UTriangular(:)
  
    ! Calculate the number of non-zero elements in the upper triangular matrix
    integer :: i, j 

    nnz_UpperTriangular = 0
    do i = 1, n
      do j = IA(i), IA(i + 1) - 1
        if (JA(j) >= i) nnz_UpperTriangular = nnz_UpperTriangular + 1
      end do
    end do
    
    ! Allocate memory for the upper triangular matrix
    allocate(UTriangular(nnz_UpperTriangular))
    allocate(IA_UTriangular(n + 1))
    allocate(JA_UTriangular(nnz_UpperTriangular))
  
   ! Set the first element of IA_UTriangular to 1
    IA_UTriangular = 0 
    IA_UTriangular(1) = 1

    ! ! Calculate the CSR format for the upper triangular matrix
    ! nnz_UpperTriangular = 0
    ! do i = 1, n
    !   IA_UTriangular(i) = nnz_UpperTriangular + 1
    !   do j = IA(i), IA(i + 1) - 1
    !     if (JA(j) >= i) then
    !       nnz_UpperTriangular = nnz_UpperTriangular + 1
    !       UTriangular(nnz_UpperTriangular) = A(j)
    !       JA_UTriangular(nnz_UpperTriangular) = JA(j)
    !     end if
    !   end do
    ! end do
      ! Calculate the CSR format for the upper triangular matrix
  nnz_UpperTriangular = 1
  IA_UTriangular(1) = nnz_UpperTriangular
  do i = 1, n
    do j = IA(i), IA(i + 1) - 1
      if (JA(j) >= i) then
        UTriangular(nnz_UpperTriangular) = A(j)
        JA_UTriangular(nnz_UpperTriangular) = JA(j)
        nnz_UpperTriangular = nnz_UpperTriangular + 1
      end if
    end do
    IA_UTriangular(i + 1) = nnz_UpperTriangular
  end do
  nnz_UpperTriangular = nnz_UpperTriangular - 1
    ! IA_UTriangular(n + 1) = nnz_UpperTriangular + 1

end subroutine HermitianUpperTriangularCSR





!--------------------------------------------!
!            Parallel Spmmv: y = A*x         !
!--------------------------------------------!

subroutine pamux(threads, n, x, y, a, ja, ia) !My parallelized sparse matrix-vector multiplication

    implicit none

    integer(kind=8), intent(in) :: n
    integer, intent(in) :: threads
    integer, intent(in) :: ja(*), ia(*)
    double precision, intent(in)  :: x(n), a(*)

    double precision, intent(out) :: y(n)
    
    !-----------------------------------------------------------------------
    !         A times a vector
    !-----------------------------------------------------------------------
    ! multiplies a matrix by a vector using the dot product form
    ! Matrix A is stored in compressed sparse row storage.
    !
    ! on entry:
    !----------
    ! n     = row dimension of A
    ! x     = real array of length equal to the column dimension of
    !         the A matrix.
    ! a, ja,
    !    ia = input matrix in compressed sparse row format.
    !
    ! on return:
    !-----------
    ! y     = real array of length n, containing the product y=Ax
    !
    !-----------------------------------------------------------------------
    ! local variables
    !

    double precision, save :: t = 0 
    integer, save :: i = 0, k = 0
    integer :: cntr = 0  
    !-----------------------------------------------------------------------

    !$omp parallel do private(i, k, t) num_threads(threads)
    do i = 1, n

        !
        !     compute the inner product of row i with vector x
        !

        t = 0.0d0
        do k = ia(i), ia(i+1) - 1
            ! cntr = cntr + 1
            t = t + a(k) * x(ja(k))
        end do 

        !
        !     store result in y(i)
        !

        y(i) = t
    end do 
    !$omp end parallel do

    return
end subroutine pamux 

subroutine pamux2(threads, dim, sites, nbonds, nnnbonds, t1, v1, v2, basis, bsites, hexsites, x, y) !My parallelized sparse matrix-vector multiplication. 
    !Hamiltonian matrix entry is created on the fly 

    implicit none

    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: threads, sites, nbonds, nnnbonds
    integer, intent(in) :: bsites(2, nbonds), hexsites(2, nnnbonds)
    double precision, intent(in)  :: t1, v1, v2
    double precision, intent(in)  :: x(dim)

    double precision, intent(out) :: y(dim)
    
    !-----------------------------------------------------------------------
    !         A times a vector
    !-----------------------------------------------------------------------
    ! multiplies a matrix by a vector using the dot product form
    ! Matrix A is stored in compressed sparse row storage.
    !
    ! on entry:
    !----------
    ! n     = row dimension of A
    ! x     = real array of length equal to the column dimension of
    !         the A matrix.
    ! a, ja,
    !    ia = input matrix in compressed sparse row format.
    !
    ! on return:
    !-----------
    ! y     = real array of length n, containing the product y=Ax
    !
    !-----------------------------------------------------------------------
    ! local variables
    !

    double precision :: a = 0.d0, t = 0.d0  
    integer(kind=8) :: i = 0, newst = 0, loc = 0
    integer :: s = 0, parity1 = 0, parity2 = 0
    ! integer :: id = 0

    !-----------------------------------------------------------------------
    
    
    !$omp parallel do private(i, loc, newst, a, t, s, parity1, parity2) num_threads(threads)
    do i = 1, dim
        loc = 0 
        newst = 0
        !
        !     compute the inner product of row i with vector x
        !
        a = 0.0d0
        t = 0.0d0
        !!$ id = omp_get_thread_num()
        
        do s = 1, nbonds  
            if (btest(basis(i), bsites(1, s) - 1) .and. .not. btest(basis(i), bsites(2, s) - 1)) then 
                newst = ibclr( ibset( basis(i), bsites(2, s) - 1 ), bsites(1, s) - 1 ) !Create on site 2, annihilate on site 1 
                call findstate(dim, newst, basis, loc) 
                if(loc > 0) then 
                    parity1 = popcnt( ibits( basis(i), bsites(1, s), sites ) ) !Parity of site1
                    parity2 = popcnt( ibits( ibclr(basis(i), bsites(1, s) - 1), bsites(2, s), sites ) ) !Parity of site2
                    a = - t1 * (-1)**(parity1 + parity2) 
                    t = t + a * x(loc) 
                end if 
            else if (btest(basis(i), bsites(2, s) - 1) .and. .not. btest(basis(i), bsites(1, s) - 1)) then 
                newst = ibclr( ibset( basis(i), bsites(1, s) - 1 ), bsites(2, s) - 1 ) !Create on site 2, annihilate on site 1 
                call findstate(dim, newst, basis, loc) 
                if(loc > 0) then 
                    parity1 = popcnt( ibits( basis(i), bsites(2, s), sites ) ) !Parity of site2
                    parity2 = popcnt( ibits( ibclr(basis(i), bsites(2,s) - 1), bsites(1, s), sites ) ) !Parity of site1
                    a = - t1 * (-1)**(parity1 + parity2)  
                    t = t + a * x(loc) 
                end if                    
            end if 
            if( btest(basis(i), bsites(1, s) - 1) .and. btest(basis(i), bsites(2, s) - 1) ) t = t + 0.25 * v1 * x(i)
            if( .not. (btest(basis(i), bsites(1, s) - 1)) .and. .not. (btest(basis(i), bsites(2, s) - 1)) ) t = t + 0.25 * v1 * x(i)
            if( btest(basis(i), bsites(1, s) - 1) .and. .not. (btest(basis(i), bsites(2, s) - 1)) ) t = t - 0.25 * v1 * x(i)
            if( .not. (btest(basis(i), bsites(1, s) - 1)) .and. btest(basis(i), bsites(2, s) - 1) ) t = t - 0.25 * v1 * x(i)
        end do !bond loop
        do s = 1, nnnbonds  
            if( btest(basis(i), hexsites(1, s) - 1) .and. btest(basis(i), hexsites(2, s) - 1) ) t = t + 0.25 * v2 * x(i)
            if( .not. (btest(basis(i), hexsites(1, s) - 1)) .and. .not. (btest(basis(i), hexsites(2, s) - 1)) ) t = t + 0.25 * v2 * x(i)
            if( btest(basis(i), hexsites(1, s) - 1) .and. .not. (btest(basis(i), hexsites(2, s) - 1)) ) t = t - 0.25 * v2 * x(i)
            if( .not. (btest(basis(i), hexsites(1, s) - 1)) .and. btest(basis(i), hexsites(2, s) - 1) ) t = t - 0.25 * v2 * x(i)
        end do !bond loop
        
        !
        !     store result in y(i)
        !

        y(i) = t
    end do 
    !$omp end parallel do

    return

end subroutine pamux2 

!----------------------------------------------------!
!            Parallel complex Spmmv: y = A*x         !
!----------------------------------------------------!
    
subroutine cpamux(threads, n, x, y, a, ja, ia) !My parallelized sparse matrix-vector multiplication
    implicit none
    !include "omp_lib.h"
    ! complex*16, intent(in)  :: x(*), a(*)
    double complex, intent(in)  :: x(*), a(*)
    integer(kind=8), intent(in) :: n
    integer, intent(in) :: threads
    integer, intent(in) :: ja(*), ia(*)

    ! complex*16, intent(out) :: y(*)
    double complex, intent(out) :: y(*)
    
    
    !-----------------------------------------------------------------------
    !         A times a vector
    !-----------------------------------------------------------------------
    ! multiplies a matrix by a vector using the dot product form
    ! Matrix A is stored in compressed sparse row storage.
    !
    ! on entry:
    !----------
    ! n     = row dimension of A
    ! x     = real array of length equal to the column dimension of
    !         the A matrix.
    ! a, ja,
    !    ia = input matrix in compressed sparse row format.
    !
    ! on return:
    !-----------
    ! y     = real array of length n, containing the product y=Ax
    !
    !-----------------------------------------------------------------------
    ! local variables
    !
    ! complex*16 :: t
    double complex, save :: t
    integer, save :: i, k
    !-----------------------------------------------------------------------
    
    
          !$omp parallel do private(k,t,i) num_threads(threads)
        !   do 100 i = 1, n
          do i = 1, n
    !
    !     compute the inner product of row i with vector x
    !
    
             t = 0.0d0
            !  do 99 k = ia(i), ia(i+1)-1
             do k = ia(i), ia(i+1)-1
                t = t + a(k)*x(ja(k))
    !  99      continue
             end do 
    !
    !     store result in y(i)
    !
             y(i) = t
    !  100  continue
             end do 
          !$omp end parallel do
    
    
          return
end subroutine

!-----------------------------------------!
!            Charge density wave          !
!-----------------------------------------!

subroutine cdw_d(dir, parameters, conf, unit, dim, ucx, ucy, bc, pat, nbonds, nDis, dis, basis, bsites, psi, rho, rhotot)

    implicit none

    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: conf, unit, ucx, ucy, nbonds, nDis
    integer, intent(in) :: bsites(2,nbonds)
    integer(kind=8), intent(in) :: basis(dim) 
    double precision, intent(in) :: psi(dim), dis
    character(len=*), intent(in) :: dir, parameters, bc, pat

    double precision, allocatable, intent(out) :: rho(:,:)  
    double precision, allocatable, intent(out) :: rhotot(:,:)  

    integer, save :: nbx  = 0 !Number of x bonds 
    integer, save :: i = 0, j = 0 
    integer, save :: cntr = 0 
    integer, allocatable, save :: cntr1(:) !Counter for rho just within unit cells for both sublattices
    integer, allocatable, save :: cntr2(:) !Counter for rho just within unit cells for both sublattices
    integer, allocatable, save :: cntrtot1(:) !Counter for rho on all NN bonds for both sublattices
    integer, allocatable, save :: cntrtot2(:) !Counter for rho on all NN bonds for both sublattices
    logical :: append 

    !$omp threadprivate(nbx, i, j, cntr, cntr1, cntr2, cntrtot1, cntrtot2)
 
    if((nDis > 1) .and. (dis .ne. 0.d0) .and. (conf > 1)) then 
        append = .true.
    else 
        append = .false.
    end if 

    if(bc == 'o') then 
        nbx = (2 * ucx - 1) * ucy 
    else if(bc == 'p') then
        nbx = (2 * ucx) * ucy !total number of x bonds
    end if
    
    if( allocated(cntr1) ) deallocate(cntr1)
    allocate( cntr1(ucx*ucy) )
    if( allocated(cntr2) ) deallocate(cntr2)
    allocate( cntr2(ucx*ucy) )
    if( allocated(rho) ) deallocate(rho)
    allocate( rho(ucx*ucy, 3) )

    !Calculate charge differences for bonds lying within a unit cell (to be adjusted for general patterns AB or BA as well as PBC)
    rho = 0 
    do j = 1, dim 
        cntr  = 0 
        cntr1 = 0 
        cntr2 = 0 
        do i = 1, nbx !run over all x bonds
            if ( modulo( bsites(1,i), 2 ) == 0 ) cycle !Don't count x bonds between unit cells 
            cntr = cntr + 1
            if ( btest( basis(j), bsites(1,i) - 1 ) ) cntr1(cntr) = cntr1(cntr) + 1 
            if ( btest( basis(j), bsites(2,i) - 1 ) ) cntr2(cntr) = cntr2(cntr) + 1 
            ! rho(cntr) = rho(cntr) + cntr(cntr) * psi(j)**2 
            if (pat == 'AB' ) then !A = cntr1, B = cntr2
                rho(cntr,1) = rho(cntr,1) + (cntr1(cntr) - cntr2(cntr)) * psi(j)**2 !<n_A> - <n_B>
                rho(cntr,2) = rho(cntr,2) + cntr1(cntr) * psi(j)**2 !<n_A>
                rho(cntr,3) = rho(cntr,3) + cntr2(cntr) * psi(j)**2 !<n_B>
            else if (pat == 'BA' ) then !A = cntr2, B = cntr1
                rho(cntr,1) = rho(cntr,1) + (cntr2(cntr) - cntr1(cntr)) * psi(j)**2 !<n_A> - <n_B>
                rho(cntr,2) = rho(cntr,2) + cntr2(cntr) * psi(j)**2 !<n_A>
                rho(cntr,3) = rho(cntr,3) + cntr1(cntr) * psi(j)**2 !<n_B>
            end if 
        end do 
    end do 

    if( allocated(cntrtot1) ) deallocate(cntrtot1)
    allocate( cntrtot1(nbonds) )
    if( allocated(cntrtot2) ) deallocate(cntrtot2)
    allocate( cntrtot2(nbonds) )
    if( allocated(rhotot) ) deallocate(rhotot)
    allocate( rhotot(nbonds, 3) )

    rhotot = 0 
    do j = 1, dim  
        cntrtot1 = 0 
        cntrtot2 = 0 
        do i = 1, nbonds 
            ! Determine sign for each sublattice (to be fixed to A=+1 and B=-1)
            ! if ( modulo( bsites(1,i), 2 ) == 1 ) then 
            !     sgn = 1
            ! else if ( modulo( bsites(1,i), 2 ) == 0 ) then 
            !     sgn = - 1
            ! end if   
            ! if( btest(basis(j), bsites(1,i) - 1) .and. (.not. btest(basis(j), bsites(2,i) - 1) ) ) cntrtot1(i) = cntrtot1(i) + 1 
            ! if( btest(basis(j), bsites(2,i) - 1) .and. (.not. btest(basis(j), bsites(1,i) - 1) ) ) cntrtot2(i) = cntrtot2(i) + 1 
            if( btest(basis(j), bsites(1,i) - 1) ) cntrtot1(i) = cntrtot1(i) + 1 
            if( btest(basis(j), bsites(2,i) - 1) ) cntrtot2(i) = cntrtot2(i) + 1 

            
            ! rhotot(i,1) = rhotot(i,1) + cntrtot(i,1) * psi(j)**2 
            if (pat == 'AB' .and. modulo( bsites(1,i), 2 ) == 1 ) then !A = cntrtot1, B = cntrtot2
                rhotot(i,1) = rhotot(i,1) + (cntrtot1(i) - cntrtot2(i)) * psi(j)**2 !<n_A> - <n_B>
                rhotot(i,2) = rhotot(i,2) + cntrtot1(i) * psi(j)**2 !<n_A>
                rhotot(i,3) = rhotot(i,3) + cntrtot2(i) * psi(j)**2 !<n_B>
            else if (pat == 'AB' .and. modulo( bsites(1,i), 2 ) == 0 ) then !A = cntrtot2 B = cntrtot1
                rhotot(i,1) = rhotot(i,1) + (cntrtot2(i) - cntrtot1(i)) * psi(j)**2 !<n_A> - <n_B>
                rhotot(i,2) = rhotot(i,2) + cntrtot2(i) * psi(j)**2 !<n_A>
                rhotot(i,3) = rhotot(i,3) + cntrtot1(i) * psi(j)**2 !<n_B>
            else if (pat == 'BA' .and. modulo( bsites(1,i), 2 ) == 1 ) then !A = cntrtot2 B = cntrtot1
                rhotot(i,1) = rhotot(i,1) + (cntrtot2(i) - cntrtot1(i)) * psi(j)**2 !<n_A> - <n_B>
                rhotot(i,2) = rhotot(i,2) + cntrtot2(i) * psi(j)**2 !<n_A>
                rhotot(i,3) = rhotot(i,3) + cntrtot1(i) * psi(j)**2 !<n_B>
            else if (pat == 'BA' .and. modulo( bsites(1,i), 2 ) == 0 ) then !A = cntrtot1, B = cntrtot2
                rhotot(i,1) = rhotot(i,1) + (cntrtot1(i) - cntrtot2(i)) * psi(j)**2 !<n_A> - <n_B>
                rhotot(i,2) = rhotot(i,2) + cntrtot1(i) * psi(j)**2 !<n_A>
                rhotot(i,3) = rhotot(i,3) + cntrtot2(i) * psi(j)**2 !<n_B>
            end if 
        end do 

    end do 

    if( allocated(cntr1) ) deallocate(cntr1)
    if( allocated(cntr2) ) deallocate(cntr2)
    if( allocated(cntrtot1) ) deallocate(cntrtot1)
    if( allocated(cntrtot2) ) deallocate(cntrtot2)

    call save(dir, parameters, append, unit, ucx, ucy, nbonds, rho, rhotot)
    return 

end subroutine cdw_d 

subroutine ddcf(dir, refsite, conf, unit, dim, sites, k1, k2, nDis, dis, v1, v2, basis, psi)

    implicit none

    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: refsite, conf, unit, sites, nDis, k1, k2
    integer(kind=8), intent(in) :: basis(dim) 
    double precision, intent(in) :: psi(dim), dis, v1, v2
    character(len=*), intent(in) :: dir

    double precision :: refrho
    double precision, allocatable :: rhocn(:), rho(:), rhoij(:)  

    integer :: i = 0, j = 0 
    integer :: flag = 0 
    logical :: append 
    character*400 :: params 

    !$omp threadprivate(i, j)
    call parddcf(unit, k1, k2, refsite, v1, v2, params)
    if((nDis > 1) .and. (dis .ne. 0.d0) .and. (conf > 1)) then 
        append = .true.
    else 
        append = .false.
    end if 
    
    if(allocated(rho)) deallocate(rho)
    if(allocated(rhoij)) deallocate(rhoij)
    if(allocated(rhocn)) deallocate(rhocn)
    allocate(rho(sites))
    allocate(rhoij(sites))
    allocate(rhocn(sites))

    rho = 0.d0 !density correlation function of site i 
    refrho = 0.d0 !density correlation function of refsite 
    rhoij = 0.d0 !density-density correlation function between sites i and refsite 
    rhocn = 0.d0 !Connected density-density correlation function
    do i = 1, sites !run over all sites
        do j = 1, dim         
            if(btest(basis(j), i - 1)) then 
                rho(i) = rho(i) + psi(j)**2 
                flag = 1
                if(btest(basis(j), refsite - 1)) then     
                    if(i == 1) refrho = refrho + psi(j)**2
                    rhoij(i) = rhoij(i) + psi(j)**2 
                end if 
            end if 
        end do 
        rhocn(i) = rhoij(i) - rho(i) * refrho
    end do 
    rhocn = rhocn / sites 
    
    call saveddcf(dir, params, append, unit, sites, rhocn, rho)
    if(allocated(rho)) deallocate(rho)
    if(allocated(rhoij)) deallocate(rhoij)
    if(allocated(rhocn)) deallocate(rhocn)
    return 

end subroutine ddcf 

subroutine cdw_c(dir, parameters, conf, unit, dim, ucx, ucy, bc, pat, nbonds, nDis, dis, basis, bsites, psi, rho, rhotot)

    implicit none

    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: conf, unit, ucx, ucy, nbonds, nDis
    integer, intent(in) :: bsites(2,nbonds)
    integer(kind=8), intent(in) :: basis(dim) 
    double complex, intent(in) :: psi(dim), dis 
    character(len=*), intent(in) :: dir, parameters, bc, pat

    double precision, allocatable, intent(out) :: rho(:,:)  
    double precision, allocatable, intent(out) :: rhotot(:,:)  

    integer, save :: nbx  = 0 !Number of x bonds 
    integer, save :: i = 0, j = 0 
    integer, save :: cntr = 0 
    integer, allocatable, save :: cntr1(:) !Counter for rho just within unit cells for both sublattices
    integer, allocatable, save :: cntr2(:) !Counter for rho just within unit cells for both sublattices
    integer, allocatable, save :: cntrtot1(:) !Counter for rho on all NN bonds for both sublattices
    integer, allocatable, save :: cntrtot2(:) !Counter for rho on all NN bonds for both sublattices
    logical :: append 

    !$omp threadprivate(nbx, i, j, cntr, cntr1, cntr2, cntrtot1, cntrtot2)
   
    if((nDis > 1) .and. (dis .ne. 0.d0) .and. (conf > 1)) then 
        append = .true.
    else 
        append = .false.
    end if 
    if(bc == 'o') then 
        nbx = (2 * ucx - 1) * ucy 
    else if(bc == 'p') then
        nbx = (2 * ucx) * ucy !total number of x bonds
    end if
    
    if( allocated(cntr1) ) deallocate(cntr1)
    allocate( cntr1(ucx*ucy) )
    if( allocated(cntr2) ) deallocate(cntr2)
    allocate( cntr2(ucx*ucy) )
    if( allocated(rho) ) deallocate(rho)
    allocate( rho(ucx*ucy, 3) )

    !Calculate charge differences for bonds lying within a unit cell (to be adjusted for general patterns AB or BA as well as PBC)
    rho = 0 
    do j = 1, dim 
        cntr  = 0 
        cntr1 = 0 
        cntr2 = 0 
        do i = 1, nbx !run over all x bonds
            if ( modulo( bsites(1,i), 2 ) == 0 ) cycle !Don't count x bonds between unit cells 
            cntr = cntr + 1
            if ( btest( basis(j), bsites(1,i) - 1 ) ) cntr1(cntr) = cntr1(cntr) + 1 
            if ( btest( basis(j), bsites(2,i) - 1 ) ) cntr2(cntr) = cntr2(cntr) + 1 
            ! rho(cntr) = rho(cntr) + cntr(cntr) * psi(j)**2 
            if (pat == 'AB' ) then !A = cntr1, B = cntr2
                rho(cntr,1) = rho(cntr,1) + (cntr1(cntr) - cntr2(cntr)) * abs(psi(j))**2 !<n_A> - <n_B>
                rho(cntr,2) = rho(cntr,2) + cntr1(cntr) * abs(psi(j))**2 !<n_A>
                rho(cntr,3) = rho(cntr,3) + cntr2(cntr) * abs(psi(j))**2 !<n_B>
            else if (pat == 'BA' ) then !A = cntr2, B = cntr1
                rho(cntr,1) = rho(cntr,1) + (cntr2(cntr) - cntr1(cntr)) * abs(psi(j))**2 !<n_A> - <n_B>
                rho(cntr,2) = rho(cntr,2) + cntr2(cntr) * abs(psi(j))**2 !<n_A>
                rho(cntr,3) = rho(cntr,3) + cntr1(cntr) * abs(psi(j))**2 !<n_B>
            end if 
        end do 
    end do 

    if( allocated(cntrtot1) ) deallocate(cntrtot1)
    allocate( cntrtot1(nbonds) )
    if( allocated(cntrtot2) ) deallocate(cntrtot2)
    allocate( cntrtot2(nbonds) )
    if( allocated(rhotot) ) deallocate(rhotot)
    allocate( rhotot(nbonds, 3) )

    rhotot = 0 
    do j = 1, dim  
        cntrtot1 = 0 
        cntrtot2 = 0 
        do i = 1, nbonds 
            ! Determine sign for each sublattice (to be fixed to A=+1 and B=-1)
            ! if ( modulo( bsites(1,i), 2 ) == 1 ) then 
            !     sgn = 1
            ! else if ( modulo( bsites(1,i), 2 ) == 0 ) then 
            !     sgn = - 1
            ! end if   
            ! if( btest(basis(j), bsites(1,i) - 1) .and. (.not. btest(basis(j), bsites(2,i) - 1) ) ) cntrtot1(i) = cntrtot1(i) + 1 
            ! if( btest(basis(j), bsites(2,i) - 1) .and. (.not. btest(basis(j), bsites(1,i) - 1) ) ) cntrtot2(i) = cntrtot2(i) + 1 
            if( btest(basis(j), bsites(1,i) - 1) ) cntrtot1(i) = cntrtot1(i) + 1 
            if( btest(basis(j), bsites(2,i) - 1) ) cntrtot2(i) = cntrtot2(i) + 1 

            
            ! rhotot(i,1) = rhotot(i,1) + cntrtot(i,1) * psi(j)**2 
            if (pat == 'AB' .and. modulo( bsites(1,i), 2 ) == 1 ) then !A = cntrtot1, B = cntrtot2
                rhotot(i,1) = rhotot(i,1) + (cntrtot1(i) - cntrtot2(i)) * abs(psi(j))**2 !<n_A> - <n_B>
                rhotot(i,2) = rhotot(i,2) + cntrtot1(i) * abs(psi(j))**2 !<n_A>
                rhotot(i,3) = rhotot(i,3) + cntrtot2(i) * abs(psi(j))**2 !<n_B>
            else if (pat == 'AB' .and. modulo( bsites(1,i), 2 ) == 0 ) then !A = cntrtot2 B = cntrtot1
                rhotot(i,1) = rhotot(i,1) + (cntrtot2(i) - cntrtot1(i)) * abs(psi(j))**2 !<n_A> - <n_B>
                rhotot(i,2) = rhotot(i,2) + cntrtot2(i) * abs(psi(j))**2 !<n_A>
                rhotot(i,3) = rhotot(i,3) + cntrtot1(i) * abs(psi(j))**2 !<n_B>
            else if (pat == 'BA' .and. modulo( bsites(1,i), 2 ) == 1 ) then !A = cntrtot2 B = cntrtot1
                rhotot(i,1) = rhotot(i,1) + (cntrtot2(i) - cntrtot1(i)) * abs(psi(j))**2 !<n_A> - <n_B>
                rhotot(i,2) = rhotot(i,2) + cntrtot2(i) * abs(psi(j))**2 !<n_A>
                rhotot(i,3) = rhotot(i,3) + cntrtot1(i) * abs(psi(j))**2 !<n_B>
            else if (pat == 'BA' .and. modulo( bsites(1,i), 2 ) == 0 ) then !A = cntrtot1, B = cntrtot2
                rhotot(i,1) = rhotot(i,1) + (cntrtot1(i) - cntrtot2(i)) * abs(psi(j))**2 !<n_A> - <n_B>
                rhotot(i,2) = rhotot(i,2) + cntrtot1(i) * abs(psi(j))**2 !<n_A>
                rhotot(i,3) = rhotot(i,3) + cntrtot2(i) * abs(psi(j))**2 !<n_B>
            end if 
        end do 

    end do 

    if( allocated(cntr1) ) deallocate(cntr1)
    if( allocated(cntr2) ) deallocate(cntr2)
    if( allocated(cntrtot1) ) deallocate(cntrtot1)
    if( allocated(cntrtot2) ) deallocate(cntrtot2)

    call save(dir, parameters, append, unit, ucx, ucy, nbonds, rho, rhotot)
    return 

end subroutine cdw_c 

subroutine gsdeg_d(dir, rvec, unit, parameters, dim, nev, nest, thresh, energies, states, ndeg, gs)
    implicit none
    
    logical, intent(in) :: rvec
    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: unit, nev, nest
    double precision, intent(in) :: thresh 
    double precision, intent(in) :: energies(nev)
    double precision, intent(in) :: states(dim, nest)
    character(len=*), intent(in) :: dir, parameters 
    integer, intent(out) :: ndeg 
    double precision, allocatable, intent(out) :: gs(:)
    
    integer :: j 
    character*256 :: file 

    ndeg = 1
            
    if(rvec) then 
        if(allocated(gs)) deallocate(gs)
        allocate(gs(dim))
        gs = states(1:dim,1)
    end if 
    do j = 2, nest 
        if ( abs(energies(j) - energies(1)) <= thresh ) then 
            ! print*, 'dE=', abs(energies(j) - energies(1))
            ! if (j == 2 ) print*, 'E(1)',energies(1) 
            ! print*, 'E(j)=',energies(j)
            ndeg = ndeg + 1
            if(rvec) gs = gs + states(1:dim,j)
        else 
            exit 
        end if 
    end do 
    if ( nest .le. nev .and. abs(energies(nest + 1) - energies(1)) <= thresh ) print*,'Potentially not full degeneracy has been captured. Increase "nest".'
    if ( nest .le. ndeg ) print*,'Potentially not full degeneracy has been captured. Increase "nest".'
    ndeg = min(ndeg, nest)
    if(rvec) gs = gs/sqrt(dble(ndeg))
    print*,'Ground state degeneracy = ',ndeg 
   
    file = dir//"gs_deg_"//parameters
    file = trim_name(file)

    open(unit,file = file)
    write(unit,*) ndeg 
    close(unit)

    return 
   
end subroutine gsdeg_d

subroutine gsdeg_c(dir, rvec, unit, parameters, dim, nev, nest, thresh, energies, states, ndeg, gs)
    implicit none

    logical, intent(in) :: rvec
    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: unit, nev, nest
    double precision, intent(in) :: thresh 
    double precision, intent(in) :: energies(nev)
    double complex, intent(in) :: states(dim, nest)
    character(len=*), intent(in) :: dir, parameters 
    integer, intent(out) :: ndeg 
    double complex, allocatable, intent(out) :: gs(:)
    
    integer :: j 
    character*256 :: file 
    
    ndeg = 1

    if(rvec) then 
        if(allocated(gs)) deallocate(gs)
        allocate(gs(dim))
        gs = states(1:dim,1)
    end if 
    do j = 2, nest 
        if ( abs(energies(j) - energies(1)) <= thresh ) then 
            ! print*, 'dE=', abs(energies(j) - energies(1))
            ! if (j == 2 ) print*, 'E(1)',energies(1) 
            ! print*, 'E(j)=',energies(j)
            ndeg = ndeg + 1
            if(rvec) gs = gs + states(1:dim,j)
        else 
            exit 
        end if 
    end do 
    if ( nest .le. nev .and. abs(energies(nest + 1) - energies(1)) <= thresh ) print*,'Potentially not full degeneracy has been captured. Increase "nest".'
    if ( nest .le. ndeg ) print*,'Potentially not full degeneracy has been captured. Increase "nest".'
    ndeg = min(ndeg, nest)
    if(rvec) gs = gs/sqrt(dble(ndeg))
    print*,'Ground state degeneracy = ',ndeg 
   
    file = dir//"gs_deg_"//parameters
    file = trim_name(file)

    open(unit,file = file)
    write(unit,*) ndeg 
    close(unit)

    return 
   
end subroutine gsdeg_c

subroutine qgsdeg_d(dir, unit, parameters, dim, nev, nest, energies, states, ndeg, gs)
    implicit none
    
    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: unit, nev, nest
    double precision, intent(in) :: energies(nev)
    double precision, intent(in) :: states(dim, nest)
    character(len=*), intent(in) :: dir, parameters 
    integer, intent(out) :: ndeg 
    double precision, allocatable, intent(out) :: gs(:)
    
    integer :: j
    double precision :: dE(nest-1)
    character*256 :: file 

    ndeg = 1
            
    if(allocated(gs)) deallocate(gs)
    allocate(gs(dim))
    gs = states(1:dim,1)
    
    do j = 2, nest 
        dE(j-1) = abs(energies(j) - energies(j-1)) 
        ! print*,'j, dE(j)',j-1, dE(j-1)
    end do
    ndeg = maxloc(real(dE), dim=1) 
    
    ! print*,'dE(ndeg-1)',dE(ndeg-1)
    ! print*,'dE(ndeg)',dE(ndeg)
    
    if(ndeg > 1) then 
        do j = 2, ndeg 
            gs = gs + states(1:dim,j)
        end do 
    end if 
    gs = gs/sqrt(dble(ndeg))
    print*,'Ground state quasi-degeneracy = ',ndeg 
   
    file =dir//"quasi_gs_deg_"//parameters
    file = trim_name(file)

    open(unit,file = file)
    write(unit,*) ndeg 
    close(unit)

    return 
   
end subroutine qgsdeg_d

subroutine qgsdeg_c(dir, unit, parameters, dim, nev, nest, energies, states, ndeg, gs)
    implicit none
    
    integer(kind=8), intent(in) :: dim
    integer, intent(in) :: unit, nev, nest
    double precision, intent(in) :: energies(nev)
    double complex, intent(in) :: states(dim, nest)
    character(len=*), intent(in) :: dir, parameters 
    integer, intent(out) :: ndeg 
    double complex, allocatable, intent(out) :: gs(:)
    
    integer :: j
    double precision :: dE(nest-1)
    character*256 :: file 

    ndeg = 1
            
    if(allocated(gs)) deallocate(gs)
    allocate(gs(dim))
    gs = states(1:dim,1)
    
    do j = 2, nest 
        dE(j-1) = abs(energies(j) - energies(j-1)) 
        ! print*,'j, dE(j)',j-1, dE(j-1)
    end do
    ndeg = maxloc(real(dE), dim=1) 
    
    ! print*,'dE(ndeg-1)',dE(ndeg-1)
    ! print*,'dE(ndeg)',dE(ndeg)
    
    if(ndeg > 1) then 
        do j = 2, ndeg 
            gs = gs + states(1:dim,j)
        end do 
    end if 
    gs = gs/sqrt(dble(ndeg))
    print*,'Ground state quasi-degeneracy = ',ndeg 
   
    file =dir//"quasi_gs_deg_"//parameters
    file = trim_name(file)

    open(unit,file = file)
    write(unit,*) ndeg 
    close(unit)

    return 
   
end subroutine qgsdeg_c

subroutine corrmat(ndeg, dim, psi, evals, mat)

    implicit none

    integer, intent(in) :: ndeg
    integer(kind=8), intent(in) :: dim
    double precision, intent(in) :: psi(dim, ndeg)
    double precision, allocatable, intent(out) :: evals(:), mat(:,:)

    integer :: i, j 

    if(allocated(evals)) deallocate(evals)
    if(allocated(mat)) deallocate(mat)
    allocate(evals(ndeg))
    allocate(mat(ndeg, ndeg))
    
    do i = 1, ndeg 
        do j = 1, ndeg 
            mat(i, j) = dot_product(psi(1:dim, i), psi(1:dim, j))
            print* ,mat(i, j), i, j, 'mat(i, j), i, j'
        end do 
    end do 

    call exactdiag(.True., ndeg, mat, evals)
    do i = 1, ndeg
        print* ,evals(i), i, 'evals(i), i'
    end do 

end subroutine corrmat

!----------------------------------------------!
!            Sublattice current                !
!----------------------------------------------!

subroutine currentcorrelations(dir, ti, tilted, nHel, tilt, threads, conf, degeneracy, full, feast, mkl, arpack, symmetrize, sites, l2, l1, ucx, ucy, k1, k2, id, par, rot, nDis, ndeg, refbonds, cntrA, cntrB, unit, dim, alattice, blattice, xyA, xyB, xtransl, ytransl, refl, c6, basis, v1, v2, dis, nnnVec, norm, eigstate, eigstate_dc)

    implicit none 
    
    integer, intent(in) :: ti, tilted, nHel, tilt, threads, conf, degeneracy, full, feast, mkl, arpack, symmetrize, sites, l2, l1, ucx, ucy, k1, k2, nDis, ndeg, refbonds, cntrA, cntrB, unit
    integer(kind=8), intent(in) :: dim
    integer, allocatable, intent(in) :: alattice(:,:), blattice(:,:), xyA(:,:), xyB(:,:), xtransl(:,:), ytransl(:,:), refl(:,:), c6(:) 
    integer(kind=8), allocatable, intent(in) :: basis(:)
    double precision, intent(in) :: v1, v2, dis, nnnVec(2,3), id, par(6), rot(5)
    double precision, allocatable, intent(in) :: norm(:), eigstate(:,:)
    double complex, allocatable, intent(in) :: eigstate_dc(:,:)
    character(len=*), intent(in) :: dir

    integer :: numthrds = 0, thread = 0, j = 0, nrefb = 0, l11 = 0, l22 = 0
    integer, save :: refbond = 0
    double precision :: current = 0.d0 
    double precision, allocatable :: bondcurrent(:)
    character :: sl*1, dirq*512
    character, save :: current_parameters*512
    logical :: append 
    !!$omp threadprivate(refbond, j, nrefb, current, bondcurrent, current_parameters, sl)
    !!$omp threadprivate(refbond, current_parameters)
    
    if(ti == 1 .and. tilted == 1) then 
        l11 = l2
        l22 = l1 
    else 
        l11 = ucx 
        l22 = ucy 
    end if 


    if((nDis > 1) .and. (dis .ne. 0.d0) .and. (conf > 1)) then 
        append = .true.
    else 
        append = .false.
    end if 

    if(degeneracy < 2 ) dirq = dir  
    if(degeneracy == 2 ) dirq = trim_name(dir//"QD_")
    dirq = trim_name(dirq)

    print*,'Calculate current current correlations...'
    print*,''
    if(refbonds > 0) then 
        nrefb = refbonds 
    else if(refbonds == 0) then 
        nrefb = cntrA
    end if 
    sl = "A"
    !!$ call omp_set_num_threads(threads)

    !!$omp parallel do default(firstprivate) shared(sites, particles, dim, basis) num_threads(threads)
    do j = 1, nrefb
        refbond = Alattice(4,j)
        call par2(j + 11, k1, k2, refbond, v1, v2, current_parameters)
        
        !!$ numthrds = omp_get_num_threads()
        !!$ thread = omp_get_thread_num()
        
        if( ndeg == 1 ) then  
            if(ti == 0) then 
                call slcurrent(1.0d0, dim, sites, basis, Alattice, xyA, cntrA, j, eigstate(1:dim,1), nnnVec, bondcurrent, current)
                if(tilted == 1) call current_irrep(nHel, tilt, 1.0d0, dim, sites, 1, 1, ti, symmetrize, id, par, rot, basis, Alattice, xtransl, ytransl, refl, c6, cntrA, j, eigstate(1:dim,1), nnnVec, norm, bondcurrent, current)
            else if(ti == 1 .and. k1 == 0 .and. k2 == 0) then 
                if(tilted == 1) call current_irrep(nHel, tilt, 1.0d0, dim, sites, l2, l1, ti, symmetrize, id, par, rot, basis, Alattice, xtransl, ytransl, refl, c6, cntrA, j, eigstate(1:dim,1), nnnVec, norm, bondcurrent, current)
            else if(ti == 1 .and. ((k1 .ne. 0) .or. (k2 .ne. 0))) then 
                call current_c(tilted, nHel, tilt, k1, k2, 1.0d0, dim, sites, l11, l22, basis, Alattice, xyA, xtransl, ytransl, cntrA, j, eigstate_dc(1:dim,1), nnnVec, norm, bondcurrent, current)
            end if 
        else if(ndeg > 1) then 
            if(ti == 0) then 
                if(tilted == 1) call current_irrep_deg(ndeg, nHel, tilt, 1.d0, dim, sites, 1, 1, ti, symmetrize, id, par, rot, basis, Alattice, xtransl, ytransl, refl, c6, cntrA, refbond, eigstate, nnnVec, norm, bondcurrent, current)
                if(tilted == 0) call current_irrep_deg_rect(ndeg, 1.0d0, dim, sites, 1, 1, ti, symmetrize, id, par, rot, basis, Alattice, xtransl, ytransl, refl, c6, cntrA, refbond, eigstate, nnnVec, norm, bondcurrent, current)
            else if(ti == 1 .and. k1 == 0 .and. k2 == 0) then 
                if(tilted == 1) then 
                    ! call current_irrep_deg(ndeg, nHel, tilt, 1.d0, dim, sites, l2, l1, ti, symmetrize, id, par, rot, basis, Alattice, xtransl, ytransl, refl, c6, cntrA, refbond, eigstate, nnnVec, norm, bondcurrent, current)
                    call jij_irrep_deg(ndeg, 1.0d0, dim, sites, l2, l1, nHel, symmetrize, id, par, rot, basis, Alattice, xtransl, ytransl, refl, c6, cntrA, eigstate, nnnVec, norm) !, bcurrent, current
                else if(tilted == 0) then 
                    call current_irrep_deg_rect(ndeg, 1.0d0, dim, sites, l2, l1, ti, symmetrize, id, par, rot, basis, Alattice, xtransl, ytransl, refl, c6, cntrA, refbond, eigstate, nnnVec, norm, bondcurrent, current)
                    call jij_irrep_deg_rect(ndeg, l2, l1, id, par, rot, basis, Alattice, xtransl, ytransl, refl, c6, cntrA, eigstate, nnnVec, norm)
                    
                end if 
                
            else if(ti == 1 .and. ((k1 .ne. 0) .or. (k2 .ne. 0))) then 
                call current_ics_c2(tilted, nHel, tilt, k1, k2, ndeg, 1.d0, dim, sites, l2, l1, basis, Alattice, xyA, xtransl, ytransl, cntrA, j, eigstate_dc(1:dim,1:ndeg), nnnVec, norm, bondcurrent, current)    
                call current_cs_c2(tilted, nHel, tilt, dirq, j + 11, current_parameters, sl, ndeg, k1, k2, 1.0d0, dim, sites, l2, l1, basis, Alattice, xyA, xtransl, ytransl, cntrA, j, eigstate_dc(1:dim,1:ndeg), nnnVec, norm)
                call currmat_c2(tilted, nHel, tilt, dirq, j + 11, current_parameters, sl, ndeg, k1, k2, 1.d0, dim, sites, l2, l1, basis, Alattice, xyA, xtransl, ytransl, cntrA, j, eigstate_dc(1:dim,1:ndeg), nnnVec, norm)
            end if 
        end if 
        
        write (* ,"(' A lattice, refbond ',i2,': Current=',f20.17)") refbond, current       
        call save(dir, sl, current_parameters, append, degeneracy, j + 11, full, feast, mkl, arpack, cntrA, current, bondcurrent)   
        
    end do !Loop over bonds 
    !!$omp end parallel do 
    
    
    if(refbonds > 0) then 
        nrefb = refbonds 
    else if(refbonds == 0) then 
        nrefb = cntrB
    end if 
    
    ! sl = "B"
    ! !!$omp parallel do default(firstprivate) shared(sites, particles, dim, basis) num_threads(threads)
    ! do j = 1, nrefb
        
    !     !!$ numthrds = omp_get_num_threads()
    !     !!$ thread = omp_get_thread_num()

    !     refbond = Blattice(4,j)
    !     call par2(unit, k1, k2, refbond, v1, v2, current_parameters)
    !     if( ndeg == 1 ) then  
    !         if(ti == 0) then 
    !             call slcurrent(1.0d0, dim, sites, basis, Blattice, xyB, cntrB, j, eigstate(1:dim,1), nnnVec, bondcurrent, current)
    !         else if(ti == 1 .and. k1 == 0 .and. k2 == 0) then 
    !             if(tilted == 1) call current_irrep(nHel, tilt, 1.0d0, dim, sites, l2, l1, ti, symmetrize, id, par, rot, basis, Blattice, xtransl, ytransl, refl, c6, cntrB, j, eigstate(1:dim,1), nnnVec, norm, bondcurrent, current)
    !         else if(ti == 1 .and. ((k1 .ne. 0) .or. (k2 .ne. 0))) then 
    !             call current_c(tilted, nHel, tilt, k1, k2, 1.0d0, dim, sites, l11, l22, basis, Blattice, xyB, xtransl, ytransl, cntrB, j, eigstate_dc(1:dim,1), nnnVec, norm, bondcurrent, current)
    !         end if 
    !     else if ( ndeg > 1 ) then 
    !         if(ti == 0) then 
    !             call current_ics(ndeg, 1.0d0, dim, sites, basis, Blattice, cntrB, j, eigstate(1:dim,1:ndeg), nnnVec, bondcurrent, current)
    !             call currmat(dirq, j + 11, current_parameters, sl, ndeg, 1.0d0, dim, sites, basis, Blattice, xyB, cntrB, j, eigstate(1:dim,1:ndeg), nnnVec)
    !             call current_cs(dirq, j + 11, current_parameters, sl, ndeg, 1.0d0, dim, sites, basis, Blattice, cntrB, j, eigstate(1:dim,1:ndeg), nnnVec)
    !         else if(ti == 1 .and. k1 == 0 .and. k2 == 0) then 
    !             if(tilted == 1) call current_irrep_deg(ndeg, nHel, tilt, 1.d0, dim, sites, l2, l1, ti, symmetrize, id, par, rot, basis, Blattice, xtransl, ytransl, refl, c6, cntrB, refbond, eigstate, nnnVec, norm, bondcurrent, current)
    !             ! if(tilted == 0) call jij_irrep_deg_rect(ndeg, 1.0d0, dim, sites, l2, l1, symmetrize, id, par, rot, basis, Blattice, xtransl, ytransl, refl, c6, cntrB, eigstate, nnnVec, norm)
    !             if(tilted == 0) call current_irrep_deg_rect(ndeg, 1.0d0, dim, sites, l2, l1, ti, symmetrize, id, par, rot, basis, Blattice, xtransl, ytransl, refl, c6, cntrB, refbond, eigstate, nnnVec, norm, bondcurrent, current)
    !             if(tilted == 1) call current_irrep_deg(ndeg, nHel, tilt, 1.d0, dim, sites, l2, l1, ti, symmetrize, id, par, rot, basis, Blattice, xtransl, ytransl, refl, c6, cntrB, refbond, eigstate, nnnVec, norm, bondcurrent, current)
    !             ! if(tilted == 0) call jij_irrep_deg_rect(ndeg, 1.0d0, dim, sites, l2, l1, symmetrize, id, par, rot, basis, Blattice, xtransl, ytransl, refl, c6, cntrB, eigstate, nnnVec, norm)
    !             if(tilted == 0) call current_irrep_deg_rect(ndeg, 1.0d0, dim, sites, l2, l1, ti, symmetrize, id, par, rot, basis, Blattice, xtransl, ytransl, refl, c6, cntrB, refbond, eigstate, nnnVec, norm, bondcurrent, current)
    !         else if(ti == 1 .and. ((k1 .ne. 0) .or. (k2 .ne. 0))) then 
    !             call current_ics_c2(tilted, nHel, tilt, k1, k2, ndeg, 1.d0, dim, sites, l2, l1, basis, Blattice, xyB, xtransl, ytransl, cntrB, j, eigstate_dc(1:dim,1:ndeg), nnnVec, norm, bondcurrent, current)    
    !             call current_cs_c2(tilted, nHel, tilt, dirq, j + 11, current_parameters, sl, ndeg, k1, k2, 1.0d0, dim, sites, l2, l1, basis, Blattice, xyB, xtransl, ytransl, cntrB, j, eigstate_dc(1:dim,1:ndeg), nnnVec, norm)
    !             call currmat_c2(tilted, nHel, tilt, dirq, j + 11, current_parameters, sl, ndeg, k1, k2, 1.d0, dim, sites, l2, l1, basis, Blattice, xyB, xtransl, ytransl, cntrB, j, eigstate_dc(1:dim,1:ndeg), nnnVec, norm)
    !         end if 
    !     end if 

    !     write (* ,"(' B lattice, refbond ',i2,': Current=',f20.17)") refbond, current                
    !     call save(dir, sl, current_parameters, append, degeneracy, unit, full, feast, mkl, arpack, cntrB, current, bondcurrent)

    ! end do !Loop over bonds 
    ! !!$omp end parallel do 

    return 
end subroutine currentcorrelations

subroutine slcurrent(t, dim, sites, basis, bsites, xy, nbonds, refbond, psi, nnnVec, bondcurrent, current)

    implicit none
    ! save 
    integer, intent(in) :: nbonds, refbond, sites
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    integer, intent(in) :: xy(4,nbonds)
    double precision, intent(in) :: t 
    double precision, intent(in) :: psi(dim), nnnVec(2, 3)
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bondcurrent(:)
    
    integer(kind=8) :: loc = 0, newst = 0
    integer :: i = 0, j = 0
    integer :: sign = 0
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec, coord1, coord2 
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bondcurrent)) deallocate(bondcurrent)
    if( allocated(psiprime) ) deallocate(psiprime)
    allocate(bondcurrent(nbonds))
    allocate(psiprime(dim))

    psiprime    = 0.d0
    refvec      = vecs(1:2,bsites(5,refbond))
    coord1(1:2) = int(0.5*xy(3,refbond))*vecs(1:2,1) + xy(4,refbond)*vecs(1:2,2) - (int(0.5*xy(1,refbond))*vecs(1:2,1) + xy(2,refbond)*vecs(1:2,2))
    refsite1    = bsites(1,refbond)
    refsite2    = bsites(2,refbond)      

    bondcurrent = 0.d0
    current = 0.d0     
    do i = 1, nbonds            
        bondcurrent(i) = 0 
        site1 = bsites(1,i)
        site2 = bsites(2,i)
        if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
            .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
    
        
        coord2(1:2) = int(0.5*xy(3,i))*vecs(1:2,1) + xy(4,i)*vecs(1:2,2) - (int(0.5*xy(1,i))*vecs(1:2,1) + xy(2,i)*vecs(1:2,2))

        phase    = bsites(3,i)   
        vec      = vecs(1:2,bsites(5,i))
        dist     = dot_product(vec, refvec)
        psiprime = 0.d0  
        do j = 1, dim
            if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                    ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                    newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1   
                    call findstate(dim, newst, basis, loc) 
                    if(loc <= 0) cycle   
                    parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                    parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                    psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j)
                else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                    ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                    newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2  
                    call findstate(dim, newst, basis, loc) 
                    if(loc <= 0) cycle   
                    parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                    parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                    psiprime(loc) = psiprime(loc) + (-1) * sign * (-(t)**2) * psi(j)
                end if     
            else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                    ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                    newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1   
                    call findstate(dim, newst, basis, loc) 
                    if(loc <= 0) cycle   
                    parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                    parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                    psiprime(loc) = psiprime(loc) + (-1) * sign * (-(t)**2) * psi(j)
                else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                    ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                    newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2  
                    call findstate(dim, newst, basis, loc) 
                    if(loc <= 0) cycle   
                    parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                    parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                    psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j)
                end if     
            end if 
        end do !Loop over basis states                  
        bondcurrent(i) = bondcurrent(i) + dist * dot_product(psi, psiprime)
        current = current + dble(phase) * bondcurrent(i)
    end do !Loop over bonds 
    
    if(allocated(psiprime)) deallocate(psiprime)

    return

end subroutine slcurrent

!LOC and J reversed (ti = 1, k1 = 0, k2 = 0, ndeg = 1) (for all clusters)
subroutine current_irrep(nHel, tilt, t, dim, sites, Lx, Ly, ti, symmetrize, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, refbond, psi, nnnVec, norm, bcurrent, current)

    implicit none
    integer, intent(in) :: nHel, tilt, nbonds, refbond, sites, Lx, Ly, ti, symmetrize, bsites(5,nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    integer(kind=8), intent(in) :: dim, basis(dim)
    double precision, intent(in) :: id, mir(6), rot(5), t, psi(dim), nnnVec(2, 3), norm(dim) 
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bcurrent(:)
    
    integer(kind=8) :: loc = 0, newst = 0, rep = 0, state = 0
    integer :: x, y, info, i, j, k, l1, l2
    integer :: sign = 1, signrep = 1, signstate = 1, signDir = 1, order = 0 
    integer :: site1, site2, refsite1, refsite2, phase, parity
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bcurrent, psiprime, psiprimeij, psiprimeprime)

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
   
    if(symmetrize == 1 .and. ti == 1) then 
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
                    call reflect(basis(j), basis(j), sites, refl(k - 1, 1:sites), signstate, info, state)
                    signstate = signstate * mir(k - 1)
                else if(8 <= k) then 
                    call c6n(basis(j), basis(j), sites, k - 7, c6, signstate, info, state)
                    signstate = signstate * rot(k - 7)
                end if 
                
                do x = 1, Lx
                    if(ti == 1) call xtranslate(state, nHel, sites, Lx, xtransl, signstate, state)

                    do y = 1, Ly
                        if(ti == 1) call ytranslate(state, nHel, tilt, sites, Lx, ytransl, signstate, state)
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

                        if(ti == 1) call representative_irrep(newst, sites, nHel, tilt, Lx, Ly, symmetrize, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                        if(ti == 0) rep = newst
                        call findstate(dim, rep, basis, loc) 
                        
                        if(loc <= 0) cycle   
                        if(abs(signstate) > 1) stop 'signstate > 1' 
                        if(abs(signrep) > 1)   stop 'signrep > 1' 
                        if(abs(signDir) > 1)   stop 'signDir > 1' 
                        
                        sign = signstate * signrep * signDir * (-1)**parity 
                        ! psiprime(j) = psiprime(j) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(loc)
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

end subroutine current_irrep

subroutine current_rect(t, dim, sites, Lx, Ly, ti, symmetrize, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, refbond, psi, nnnVec, norm, bcurrent, current)

    implicit none
    integer, intent(in) :: nbonds, refbond, sites, Lx, Ly, ti, symmetrize
    integer, intent(in) :: bsites(5,nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    integer(kind=8), intent(in) :: dim, basis(dim)
    double precision, intent(in) :: id, mir(6), rot(5), t, psi(dim), nnnVec(2, 3), norm(dim) 
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bcurrent(:)
    
    integer(kind=8) :: loc = 0, newst = 0, rep = 0, state = 0
    integer :: x, y, info, i, j, k, l1, l2
    integer :: sign = 1, signrep = 1, signstate = 1, signDir = 1, order = 0 
    integer :: site1, site2, refsite1, refsite2, phase, parity
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bcurrent, psiprime, psiprimeij, psiprimeprime)

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
   
    if(symmetrize == 1 .and. ti == 1) then 
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
                    call reflect(basis(j), basis(j), sites, refl(k - 1, 1:sites), signstate, info, state)
                    signstate = signstate * mir(k - 1)
                else if(8 <= k) then 
                    call c6n(basis(j), basis(j), sites, k - 7, c6, signstate, info, state)
                    signstate = signstate * rot(k - 7)
                end if 
                
                do x = 1, Lx
                    call xtranslate(state, Ly, sites, Lx, xtransl, signstate, state)

                    do y = 1, Ly
                        call ytranslate(state, Ly, 0, sites, Lx, ytransl, signstate, state)
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
                        if(ti == 1) call representative_irrep(newst, sites, Ly, 0, Lx, Ly, symmetrize, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                        if(ti == 0) rep = newst
                        call findstate(dim, rep, basis, loc) 
                        
                        if(loc <= 0) cycle   

                        sign = signstate * signrep * signDir * (-1)**parity 
                        ! psiprime(j) = psiprime(j) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(loc)
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

end subroutine current_rect

subroutine current_irrep_deg(ndeg, nHel, tilt, t, dim, sites, Lx, Ly, ti, symmetrize, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, refbond, psi, nnnVec, norm, bcurrent, current)

    implicit none
    integer, intent(in) :: ndeg, nHel, tilt, nbonds, refbond, sites, Lx, Ly, ti, symmetrize, bsites(5,nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    integer(kind=8), intent(in) :: dim, basis(dim)
    double precision, intent(in) :: id, mir(6), rot(5), t, psi(dim, ndeg), nnnVec(2, 3), norm(dim) 
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bcurrent(:)
    integer(kind=8) :: loc = 0, newst = 0, rep = 0, state = 0
    integer :: x, y, info, i, j, k, l1, l2, nd1, nd2
    integer :: sign = 1, signrep = 1, signstate = 1, signDir = 1, order = 0 
    integer :: site1, site2, refsite1, refsite2, phase, parity
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:), currmat(:,:), bcurrmat(:,:,:), currents(:), evals(:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bcurrent)) deallocate(bcurrent)
    if(allocated(bcurrmat)) deallocate(bcurrmat)
    if(allocated(psiprime)) deallocate(psiprime)
    if(allocated(currents)) deallocate(currents)
    if(allocated(currmat)) deallocate(currmat)
    
    allocate(bcurrmat(nbonds, ndeg, ndeg))
    allocate(currmat(ndeg, ndeg))
    allocate(bcurrent(nbonds))
    allocate(currents(ndeg))
    allocate(psiprime(dim))

    refvec   = vecs(1:2,bsites(5,refbond))
    refsite1 = bsites(1,refbond)
    refsite2 = bsites(2,refbond)  
    bcurrent = 0.d0 
    currents = 0.d0 
    bcurrmat = 0.d0 
    current = 0.d0

    if(symmetrize == 1) then 
        order = 12
    else 
        order = 1
    end if 

    do nd1 = 1, ndeg 
        do nd2 = 1, ndeg 
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
                            call reflect(basis(j), basis(j), sites, refl(k - 1, 1:sites), signstate, info, state)
                            signstate = signstate * mir(k - 1)
                        else if(8 <= k) then 
                            call c6n(basis(j), basis(j), sites, k - 7, c6, signstate, info, state)
                            signstate = signstate * rot(k - 7)
                        end if 
                        
                        do x = 1, Lx
                            if(ti == 1) call xtranslate(state, nHel, sites, Lx, xtransl, signstate, state)

                            do y = 1, Ly
                                if(ti == 1) call ytranslate(state, nHel, tilt, sites, Lx, ytransl, signstate, state)
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

                                if(ti == 1) call representative_irrep(newst, sites, nHel, tilt, Lx, Ly, symmetrize, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                                if(ti == 0) rep = newst
                                call findstate(dim, rep, basis, loc) 
                                
                                if(loc <= 0) cycle   
                                
                                sign = signstate * signrep * signDir * (-1)**parity 
                                ! psiprime(j) = psiprime(j) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(loc)
                                if(ti == 1) psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(j, nd1)
                                if(ti == 0) psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j, nd1)
                            end do !y
                        end do !x 
                    end do !G 
                    
                end do !j 
                bcurrmat(i,nd2,nd1) = bcurrmat(i,nd2,nd1) + dble(phase) * dist * dot_product(psi(1:dim, nd2), psiprime)
                currmat(nd2,nd1) = currmat(nd2,nd1) + bcurrmat(i,nd2,nd1)
            end do !Loop over bonds    
        end do 
    end do !Degenerate states
    
    bcurrmat = bcurrmat / ndeg
    do i = 1, nbonds          
        call exactdiag(.True., ndeg, bcurrmat(i,1:ndeg,1:ndeg), evals)
        bcurrent(i) = sum(evals) 
    end do 

    currmat = currmat / ndeg 

    call exactdiag(.True., ndeg, currmat, currents)
    current = sum(currents)
    do i = 1, ndeg 
        print* ,currents(i), i, 'currents(i), i'
    end do
    deallocate(psiprime)

    return

end subroutine current_irrep_deg

subroutine current_irrep_deg_rect(ndeg, t, dim, sites, Lx, Ly, ti, symmetrize, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, refbond, psi, nnnVec, norm, bcurrent, current)

    implicit none
    integer, intent(in) :: ndeg, nbonds, refbond, sites, Lx, Ly, ti, symmetrize
    integer, intent(in) :: bsites(5,nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    integer(kind=8), intent(in) :: dim, basis(dim)
    double precision, intent(in) :: id, mir(6), rot(5), t, psi(dim, ndeg), nnnVec(2, 3), norm(dim) 
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bcurrent(:)
    integer(kind=8) :: loc = 0, newst = 0, rep = 0, state = 0
    integer :: x, y, info, i, j, k, l1, l2, nd1, nd2
    integer :: sign = 1, signrep = 1, signstate = 1, signDir = 1, order = 0 
    integer :: site1, site2, refsite1, refsite2, phase, parity
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:), currmat(:,:), bcurrmat(:,:,:), currents(:), evals(:), gs(:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bcurrent)) deallocate(bcurrent)
    if(allocated(bcurrmat)) deallocate(bcurrmat)
    if(allocated(psiprime)) deallocate(psiprime)
    if(allocated(currents)) deallocate(currents)
    if(allocated(currmat)) deallocate(currmat)
    
    allocate(bcurrmat(nbonds, ndeg, ndeg))
    allocate(currmat(ndeg, ndeg))
    allocate(bcurrent(nbonds))
    allocate(currents(ndeg))
    allocate(psiprime(dim))

    refvec   = vecs(1:2,bsites(5,refbond))
    refsite1 = bsites(1,refbond)
    refsite2 = bsites(2,refbond)  
    bcurrent = 0.d0 
    currents = 0.d0 
    bcurrmat = 0.d0 
    current = 0.d0

    if(symmetrize == 1) then 
        order = 12
    else 
        order = 1
    end if 

    do nd1 = 1, ndeg 
        do nd2 = 1, ndeg 
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
                            call reflect(basis(j), basis(j), sites, refl(k - 1, 1:sites), signstate, info, state)
                            signstate = signstate * mir(k - 1)
                        else if(8 <= k) then 
                            call c6n(basis(j), basis(j), sites, k - 7, c6, signstate, info, state)
                            signstate = signstate * rot(k - 7)
                        end if 
                        
                        do x = 1, Lx
                            if(ti == 1) call xtranslate(state, Ly, sites, Lx, xtransl, signstate, state)

                            do y = 1, Ly
                                if(ti == 1) call ytranslate(state, Ly, 0, sites, Lx, ytransl, signstate, state)
                                if(btest(state, refsite1-1) .and. .not.(btest(state, refsite2-1))) then !Forward hopping: 1 -> 2
                                    newst   = ibclr(ibset(state, refsite2-1), refsite1-1) !Create on refsite 2, annihilate on refsite 1   
                                    parity  = popcnt(ibits(state, refsite1, sites)) + popcnt(ibits(ibclr(state, refsite1-1), refsite2, sites))
                                    signDir = 1 
                                else if(btest(state, refsite2-1) .and. .not.(btest(state, refsite1-1))) then !Backwards hopping: 2 -> 1 
                                    newst   = ibclr(ibset(state, refsite1-1), refsite2-1) !Create on refsite 1, annihilate on refsite 2
                                    parity  = popcnt(ibits(state, refsite2, sites)) + popcnt(ibits(ibclr(state, refsite2-1), refsite1, sites))
                                    signDir = -1 
                                else 
                                    cycle 
                                end if 

                                if(btest(state, site1-1) .and. .not.(btest(state, site2-1))) then !Forward forward hopping: 1 -> 2       
                                    parity = parity + popcnt(ibits(newst, site1, sites)) + popcnt(ibits(ibclr(newst, site1-1), site2, sites))                     
                                    newst  = ibclr(ibset(newst, site2-1), site1-1) !Create on site 2, annihilate on site 1                             
                                else if(btest(state, site2-1) .and. .not.(btest(state, site1-1))) then !Backwards hopping: 2 -> 1
                                    parity  = parity + popcnt(ibits(newst, site2, sites)) + popcnt(ibits(ibclr(newst, site2-1), site1, sites)) 
                                    newst   = ibclr(ibset(newst, site1-1), site2-1) !Create on site 1, annihilate on site 2
                                    signDir = signDir * (-1)
                                else 
                                    cycle 
                                end if 

                                if(ti == 1) call representative_irrep(newst, sites, Ly, 0, Lx, Ly, symmetrize, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                                if(ti == 0) rep = newst
                                call findstate(dim, rep, basis, loc) 
                                
                                if(loc <= 0) cycle   

                                sign = signstate * signrep * signDir * (-1)**parity 
                                if(ti == 1) psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(j, nd1)
                                if(ti == 0) psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j, nd1)
                                ! psiprime(j) = psiprime(j) + sign * (-(t)**2) * sqrt(norm(j)/norm(loc)) * psi(loc, nd1)
                                
                            end do !y
                        end do !x 
                    end do !G 
                    
                end do !j 
                bcurrmat(i,nd2,nd1) = bcurrmat(i,nd2,nd1) + dble(phase) * dist * dot_product(psi(1:dim, nd2), psiprime)
                currmat(nd2,nd1) = currmat(nd2,nd1) + bcurrmat(i,nd2,nd1)
            end do !Loop over bonds    
        end do 
    end do !Degenerate states
    
    bcurrmat = bcurrmat / ndeg
    do i = 1, nbonds          
        call exactdiag(.True., ndeg, bcurrmat(i,1:ndeg,1:ndeg), evals)
        bcurrent(i) = sum(evals) 
    end do 

    currmat = currmat / ndeg 
    do i = 1, ndeg 
        do j = 1, ndeg  
            print* ,currmat(i, j), i, j, '(JJ)_ij, i, j'      
        end do 
        print*, '------------------'
    end do 
    pause 
    call exactdiag(.True., ndeg, currmat, currents)
    current = sum(currents)
    if(allocated(gs)) deallocate(gs)
    allocate(gs(dim))
    gs = 0.d0 
    do i = 1, ndeg 
        print* ,currents(i), i, 'currents(i), i'      
    end do 
    
    
    ! do i = 1, ndeg 
        !     do j = 1, ndeg 
        !         print* ,norm2(currmat(1:ndeg,i)), 'norm2(currmat(j,i))'
        !         print* ,norm2(psi(1:dim,i)), 'norm2(psi(1:dim,i))'
        !         gs = gs + currmat(j, i) * psi(1:dim, j)
        !     end do 
        !     ! gs = gs/ norm2(gs)
        !     print* ,norm2(gs), 'norm2(gs)'
        !     call current_rect(t, dim, sites, Lx, Ly, ti, symmetrize, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, refbond, gs, nnnVec, norm, bcurrent, current)
        !     print* ,current, 'current'
        !     gs = 0.d0 
        !     pause     
    ! end do 
    deallocate(psiprime)

    return

end subroutine current_irrep_deg_rect

subroutine jij_irrep_deg_rect(ndeg, Lx, Ly, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, psi, nnnVec, norm)

    use variables

    implicit none
    
    integer, intent(in) :: ndeg, nbonds, Lx, Ly
    integer, intent(in) :: bsites(5,nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    integer(kind=8), intent(in) :: basis(dim)
    double precision, intent(in) :: id, mir(6), rot(5), psi(dim, ndeg), nnnVec(2, 3), norm(dim) 
    integer(kind=8) :: loc = 0, newst = 0, rep = 0, state = 0
    integer :: x, y, info, i, j, k, l1, l2, nd1, nd2
    integer :: sign = 1, signrep = 1, signstate = 1, signDir = 1, order = 0 
    integer :: site1, site2, phase, parity
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: evals(:), bcurrent(:,:)
    double complex, allocatable :: psiprime(:), bcurrmat(:,:,:), mat(:,:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

   
    if(allocated(bcurrmat)) deallocate(bcurrmat)
    if(allocated(bcurrent)) deallocate(bcurrent)
    if(allocated(psiprime)) deallocate(psiprime)
    
    allocate(bcurrmat(nbonds, ndeg, ndeg))
    allocate(bcurrent(nbonds, ndeg))
    allocate(psiprime(dim))
  
    bcurrmat = (0.d0,0.d0) 
    bcurrent = 0.d0 

    if(symmetrize == 1) then 
        order = 12
    else 
        order = 1
    end if 

    do nd1 = 1, ndeg 
        do nd2 = 1, ndeg 
            do i = 1, nbonds !Loop over bonds    
                psiprime = 0.0d0  
                site1    = bsites(1,i)
                site2    = bsites(2,i)
                phase    = bsites(3,i)   
                vec      = vecs(1:2,bsites(5,i))
                dist     = 1.d0

                do j = 1, dim !Loop over basis states 
                    do k = 1, order !Loop over point groups operators 
                        if(k == 1) then !Identity
                            state = basis(j)
                        else if(k <= 7) then !Reflections
                            call reflect(basis(j), basis(j), sites, refl(k - 1, 1:sites), signstate, info, state)
                            signstate = signstate * mir(k - 1)
                        else if(8 <= k) then !Rotations
                            call c6n(basis(j), basis(j), sites, k - 7, c6, signstate, info, state)
                            signstate = signstate * rot(k - 7)
                        end if 
                        
                        do x = 1, Lx !X-translations
                            call xtranslate(state, Ly, sites, Lx, xtransl, signstate, state)

                            do y = 1, Ly !Y-translations
                                call ytranslate(state, Ly, 0, sites, Lx, ytransl, signstate, state)

                                if(btest(state, site1-1) .and. .not.(btest(state, site2-1))) then !Forward forward hopping: 1 -> 2       
                                    parity  = popcnt(ibits(state, site1, sites)) + popcnt(ibits(ibclr(state, site1-1), site2, sites))                     
                                    newst   = ibclr(ibset(state, site2-1), site1-1) !Create on site 2, annihilate on site 1                             
                                    signDir = 1
                                else if(btest(state, site2-1) .and. .not.(btest(state, site1-1))) then !Backwards hopping: 2 -> 1
                                    parity  = popcnt(ibits(state, site2, sites)) + popcnt(ibits(ibclr(state, site2-1), site1, sites)) 
                                    newst   = ibclr(ibset(state, site1-1), site2-1) !Create on site 1, annihilate on site 2
                                    signDir = -1
                                else 
                                    cycle 
                                end if 

                                call representative_irrep(newst, sites, Ly, 0, Lx, Ly, symmetrize, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                                call findstate(dim, rep, basis, loc) 
                                
                                if(loc <= 0) cycle   
                                                                
                                sign = signstate * signrep * signDir * (-1)**parity 
                                
                                bcurrmat(i, nd2, nd1) = bcurrmat(i, nd2, nd1) + sign * ii * sqrt(norm(loc)/norm(j)) * psi(j, nd1) * psi(loc, nd2)
                                ! bcurrmat(i, nd2, nd1) = bcurrmat(i, nd2, nd1) + dble(phase) * dist * sign * t * ii * sqrt(norm(loc)/norm(j)) * psi(j, nd1) * psi(loc, nd2)
                                ! psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(j, nd1)
                                ! psiprime(j) = psiprime(j) + sign * (-(t)**2) * sqrt(norm(j)/norm(loc)) * psi(loc, nd1)
                                
                            end do !y
                        end do !x 
                    end do !G 
                    
                end do !j 
                ! bcurrmat(i,nd2,nd1) = bcurrmat(i,nd2,nd1) + dble(phase) * dist * dot_product(psi(1:dim, nd2), psiprime)
            end do !Loop over bonds    
        end do 
    end do !Degenerate states

    if(allocated(mat)) deallocate(mat)
    allocate(mat(ndeg, ndeg))
    ! bcurrmat = bcurrmat / ndeg
    do i = 1, nbonds          
        mat = bcurrmat(i,1:ndeg,1:ndeg)
        call exactdiag_c(.True., ndeg, mat, evals)
        bcurrent(i, 1:ndeg) = evals
        
        do j = 1, ndeg 
            print* ,evals(j), j, i, 'evals(state), state, bond'
        end do 
        
        do j = 1, ndeg 
            print*, 'state'
            do k = 1, ndeg 
                print* ,mat(k,j)
            end do 
        end do 
        ! print* ,sum(evals), 'sum(evals)'
        print* , ''    
        if(ndeg == 2) pause 
        mat = (0.d0, 0.d0)  
    end do
    print* ,prms, 'prms'
    call saveD("A", 11, "output/", "Jij_", prms, nbonds, ndeg, 1.d0, bcurrent, 1) 
    if(ndeg == 2) pause 

    deallocate(psiprime)
    
    
    return

end subroutine jij_irrep_deg_rect

subroutine jij_irrep_deg(ndeg, t, dim, sites, Lx, Ly, nHel, symmetrize, id, mir, rot, basis, bsites, xtransl, ytransl, refl, c6, nbonds, psi, nnnVec, norm) !, bcurrent, current

    implicit none
    integer, intent(in) :: ndeg, nbonds, sites, Lx, Ly, nHel, symmetrize
    integer, intent(in) :: bsites(5,nbonds), xtransl(2, sites), ytransl(2, sites), refl(6, sites), c6(sites)
    integer(kind=8), intent(in) :: dim, basis(dim)
    double precision, intent(in) :: id, mir(6), rot(5), t, psi(dim, ndeg), nnnVec(2, 3), norm(dim) 
    integer(kind=8) :: loc = 0, newst = 0, rep = 0, state = 0
    integer :: x, y, info, i, j, k, l1, l2, nd1, nd2
    integer :: sign = 1, signrep = 1, signstate = 1, signDir = 1, order = 0 
    integer :: site1, site2, phase, parity
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: evals(:), bcurrent(:,:)
    double complex, allocatable :: psiprime(:), bcurrmat(:,:,:), mat(:,:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

   
    if(allocated(bcurrmat)) deallocate(bcurrmat)
    if(allocated(bcurrent)) deallocate(bcurrent)
    if(allocated(psiprime)) deallocate(psiprime)
    
    allocate(bcurrmat(nbonds, ndeg, ndeg))
    allocate(bcurrent(nbonds, ndeg))
    allocate(psiprime(dim))
  
    bcurrmat = (0.d0,0.d0) 
    bcurrent = 0.d0 

    if(symmetrize == 1) then 
        order = 12
    else 
        order = 1
    end if 

    do nd1 = 1, ndeg 
        do nd2 = 1, ndeg 
            do i = 1, nbonds !Loop over bonds    
                psiprime = 0.0d0  
                site1    = bsites(1,i)
                site2    = bsites(2,i)
                phase    = bsites(3,i)   
                vec      = vecs(1:2,bsites(5,i))
                dist     = 1.d0

                do j = 1, dim !Loop over basis states 
                    do k = 1, order !Loop over point groups operators 
                        if(k == 1) then !Identity
                            state = basis(j)
                        else if(k <= 7) then !Reflections
                            call reflect(basis(j), basis(j), sites, refl(k - 1, 1:sites), signstate, info, state)
                            signstate = signstate * mir(k - 1)
                        else if(8 <= k) then !Rotations
                            call c6n(basis(j), basis(j), sites, k - 7, c6, signstate, info, state)
                            signstate = signstate * rot(k - 7)
                        end if 
                        
                        do x = 1, Lx !X-translations
                            call xtranslate(state, nHel, sites, Lx, xtransl, signstate, state)

                            do y = 1, Ly !Y-translations
                                call ytranslate(state, nHel, 0, sites, Lx, ytransl, signstate, state)

                                if(btest(state, site1-1) .and. .not.(btest(state, site2-1))) then !Forward forward hopping: 1 -> 2       
                                    parity  = popcnt(ibits(state, site1, sites)) + popcnt(ibits(ibclr(state, site1-1), site2, sites))                     
                                    newst   = ibclr(ibset(state, site2-1), site1-1) !Create on site 2, annihilate on site 1                             
                                    signDir = 1
                                else if(btest(state, site2-1) .and. .not.(btest(state, site1-1))) then !Backwards hopping: 2 -> 1
                                    parity  = popcnt(ibits(state, site2, sites)) + popcnt(ibits(ibclr(state, site2-1), site1, sites)) 
                                    newst   = ibclr(ibset(state, site1-1), site2-1) !Create on site 1, annihilate on site 2
                                    signDir = -1
                                else 
                                    cycle 
                                end if 

                                call representative_irrep(newst, sites, nHel, 1, Lx, Ly, symmetrize, id, mir, rot, xtransl, ytransl, refl, c6, rep, l1, l2, signrep)
                                call findstate(dim, rep, basis, loc) 
                                
                                if(loc <= 0) cycle   
                                                                
                                sign = signstate * signrep * signDir * (-1)**parity 
                                
                                bcurrmat(i, nd2, nd1) = bcurrmat(i, nd2, nd1) + sign * ii * sqrt(norm(loc)/norm(j)) * psi(j, nd1) * psi(loc, nd2)
                                ! bcurrmat(i, nd2, nd1) = bcurrmat(i, nd2, nd1) + dble(phase) * dist * sign * t * ii * sqrt(norm(loc)/norm(j)) * psi(j, nd1) * psi(loc, nd2)
                                ! psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(j, nd1)
                                ! psiprime(j) = psiprime(j) + sign * (-(t)**2) * sqrt(norm(j)/norm(loc)) * psi(loc, nd1)
                                
                            end do !y
                        end do !x 
                    end do !G 
                    
                end do !j 
                ! bcurrmat(i,nd2,nd1) = bcurrmat(i,nd2,nd1) + dble(phase) * dist * dot_product(psi(1:dim, nd2), psiprime)
            end do !Loop over bonds    
        end do 
    end do !Degenerate states

    do k = 1, nbonds
        print* , k, 'bond'
        do i = 1, ndeg
            do j = 1, ndeg               
                print* ,bcurrmat(k, i, j), i, j, '(J_k)_ij, i, j'
            end do 
        end do 
        print*, '---------------'
    end do 
    pause 
    if(allocated(mat)) deallocate(mat)
    allocate(mat(ndeg, ndeg))
    ! bcurrmat = bcurrmat / ndeg
    do i = 1, nbonds          
        mat = bcurrmat(i,1:ndeg,1:ndeg)
        call exactdiag_c(.True., ndeg, mat, evals)
        bcurrent(i, 1:ndeg) = evals
        
        do j = 1, ndeg 
            print* ,evals(j), j, i, 'evals(state), state, bond'
        end do 
        
        do j = 1, ndeg 
            print*, 'state'
            do k = 1, ndeg 
                print* ,mat(k,j)
            end do 
        end do 
        ! print* ,sum(evals), 'sum(evals)'
        print* , ''    
        if(ndeg == 2) pause 
        mat = (0.d0, 0.d0)  
    end do 
    if(ndeg == 2) pause 

    deallocate(psiprime)
    
    
    return

end subroutine jij_irrep_deg

subroutine ccf(t, dim, sites, basis, bsites, phases, nbonds, refbond, psi, bondcurcurrent, current)

    implicit none
    integer, intent(in) :: nbonds, refbond, sites  
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(2,nbonds)
    integer, intent(in) :: phases(nbonds)
    double precision, intent(in) :: t 
    double precision, intent(in) :: psi(dim)

    double precision, intent(out) :: current 
    ! double precision, intent(out) :: bondcurrent(nbonds)
    double precision, intent(out) :: bondcurcurrent(nbonds)
    
    integer(kind=8) :: loc = 0, newst = 0
    integer :: i = 0, j = 0, k = 0 
    integer :: forward = 0, backwards = 0
    integer :: signf = 0, signb = 0, cntr = 0 
    integer :: site1 = 0, site2 = 0, phase = 0   

    ! double precision :: bondcurcurrent(nbonds)
    double precision :: current_temp
    double precision, allocatable :: bondcurrent(:)
    double precision, allocatable :: psiprime(:)
    double precision, allocatable :: psiprimeij(:)
    double precision, allocatable :: psiprimeprime(:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)

    !j_\alpha = -i \sum_{i,\sigma} t (c^\dagger_{i\sigma}c_{i+\delta_\alpha,\sigma} -c^\dagger_{i+\delta_\alpha,\sigma} c_{i\sigma})
    if( allocated(bondcurrent) ) deallocate(bondcurrent)
    allocate(bondcurrent(nbonds))
    current        = 0 
    bondcurrent    = 0 
    bondcurcurrent = 0 

    if( allocated(psiprime) ) deallocate(psiprime)
    if( allocated(psiprimeij) ) deallocate(psiprimeij)
    if( allocated(psiprimeprime) ) deallocate(psiprimeprime)
    allocate( psiprime(dim) )
    allocate( psiprimeij(dim) )
    allocate( psiprimeprime(dim) )
    psiprime      = 0 
    psiprimeij    = 0 
    psiprimeprime = 0 

    i = refbond 
    do j = 1, dim
        ! do k = 1, 2
        !     if(k == 1) then 
        ! if(bsites(1,i) < bsites(2,i)) then 
                site1 = bsites(1,i)
                site2 = bsites(2,i)                
        ! else 
        !         site2 = bsites(1,i)
        !         site1 = bsites(2,i)
        ! end if 
            ! else if(k == 2) then 
            !     site1 = bsites(2,i)
            !     site2 = bsites(1,i)            
            ! end if 
            current_temp   = 0 
            forward        = 0
            backwards      = 0 
            if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Hopping from 1->2
                forward = forward + 1 
                newst   = ibclr( ibset( basis(j), site2 - 1 ), site1 - 1 ) !Create site 2, annihilate site 1
                call findstate (dim, newst, basis, loc)     
                if( basis(loc) .ne. newst) error stop "Subroutine slcurrent: Wrong location of scattered basis state."
                signf = (-1)**(popcnt( ibits( basis(j), site1, sites - site1 ) ) + &
                                popcnt( ibits( ibclr(basis(j), site1 - 1), site2, sites - site2 ) ) ) !Count occupied sites over which operators have to hop. -1 for already occupied site. -1 for hopping over already annihilated site.
            else if( btest( basis(j), site2 - 1 ) .and. .not.( btest( basis(j), site1 - 1 ) ) ) then !Hopping from 2->1 
                backwards = backwards + 1 
                newst = ibclr( ibset( basis(j), site1 - 1 ), site2 - 1 ) 
                call findstate (dim, newst, basis, loc)
                if( basis(loc) .ne. newst) error stop "Subroutine slcurrent: Wrong location of scattered basis state."
                signb = (-1)**(popcnt( ibits( ibclr(basis(j), site2 - 1), site1, sites - site1 ) ) + &
                                popcnt( ibits( basis(j), site2, sites - site2 ) ) ) 
            else 
                cycle    
            end if 
            current_temp  = -t * ( signf * forward - signb * backwards )
            psiprime(loc) = psiprime(loc) + current_temp * psi(j)
            ! psiprime(loc) = psiprime(loc) + bsites(3,i) * current_temp * psi(j)
        ! end do 
    end do             

    do i = 1, nbonds 
    ! if(bsites(1,i) < bsites(2,i)) then 
            site1 = bsites(1,i)
            site2 = bsites(2,i)                     
    ! else 
    !         site2 = bsites(1,i)
    !         site1 = bsites(2,i)             
    ! end if 
            phase = phases(i)   
        if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
            .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
        psiprimeprime = 0 
        do j = 1, dim            
            current_temp   = 0 
            forward        = 0
            backwards      = 0 
            if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Hopping from 1->2                 
                forward = forward + 1 
                newst   = ibclr( ibset( basis(j), site2 - 1 ), site1 - 1 ) 
                call findstate (dim, newst, basis, loc)
                if( basis(loc) .ne. newst) error stop "Subroutine slcurrent: Wrong location of scattered basis state."
                signf = (-1)**(popcnt(ibits(basis(j), site1, sites - site1 ) ) + &
                            popcnt(ibits(ibclr(basis(j), site1 - 1), site2, sites - site2 ) ) ) !Count occupied sites over which operators have to hop. -1 for already occupied site. -1 for hopping over already annihilated site.
            else if(btest(basis(j), site2  - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then                 
                backwards = backwards + 1 
                newst     = ibclr( ibset( basis(j), site1 - 1 ), site2 - 1 ) 
                call findstate (dim, newst, basis, loc)
                if( basis(loc) .ne. newst) error stop "Subroutine slcurrent: Wrong location of scattered basis state."
                signb = (-1)**(popcnt(ibits(ibclr(basis(j), site2 - 1), site1, sites - site1 ) ) + &
                            popcnt(ibits(basis(j), site2, sites - site2 ) ) ) 
            else
                cycle 
            end if 
            current_temp       = - t * ( signf * forward - signb * backwards )
            psiprimeij(loc)    = psiprimeij(loc) + current_temp * psi(j)
            psiprimeprime(loc) = psiprimeprime(loc) + current_temp * psiprime(j)
        end do           
        bondcurrent(i)    = 0 
        bondcurrent(i)    = dot_product(psi, psiprimeij)
        bondcurcurrent(i) = 0 
        bondcurcurrent(i) = dot_product(psi, psiprimeprime)
        current           = current + phase * bondcurcurrent(i)
        ! current = current +  bondcurcurrent(i)
    end do

    deallocate(psiprime)
    deallocate(psiprimeij)
    deallocate(psiprimeprime)

    return

end subroutine ccf 




!---------------------------------!
!            Thread units         !
!---------------------------------!

subroutine threadunits(ndv, ndv2, units)

    implicit none
    integer, intent(in) :: ndv, ndv2
    integer, allocatable, intent(out) :: units(:, :)
    integer :: i, j, h
    integer :: n_threads_tot

    n_threads_tot = 1

    if (ndv > 0) n_threads_tot = n_threads_tot + ndv
    if (ndv2 > 0) n_threads_tot = n_threads_tot + ndv2

    if(allocated(units)) deallocate(units)
    allocate(units(n_threads_tot,n_threads_tot))
    h = 0
    do i = 1, n_threads_tot
        do j = 1, n_threads_tot
            h = h + 1
            units(i,j) = h
        end do
    end do

    return 

end subroutine

!-------------------------------!
!            Step units         !
!-------------------------------!

subroutine stepunits(ndu, ndv, ndv2, units)
    !Creates array with units for parallel I/O. Each entry is a different unit number.
    implicit none
    integer, intent(in) :: ndu, ndv, ndv2
    integer, allocatable, intent(out) ::  units(:,:,:,:)
    integer :: i, j, k, l, q

    if (allocated(units)) deallocate(units)
    allocate(units(ndv2,ndv,ndu,2))
    q = 0
    do i = 1, ndv2
        do j = 1, ndv
            do k = 1, ndu
                do l = 1, 2
                    q = q + 1
                    units(i,j,k,l) = q
                end do
            end do
        end do
    end do

end subroutine

!-------------------------------------------!
!            Date and Time                  !
!-------------------------------------------!

subroutine datetime(dir, startend, params)
    implicit none

    integer, intent(in) :: startend
    character(len=*), intent(in) :: dir
    character(len=*), intent(in), optional :: params
    character(8), save  :: datei, datef
    character(10), save :: timei, timef
    integer,dimension(8), save :: valuesi, valuesf
    real, save :: start, finish
    character :: file_name*400


    if(startend == 0) then
        call cpu_time(start)
        call date_and_time(date=datei,time=timei,values=valuesi)
        write(*,"('Calculation started at',x,a,' h',x,a,' min',x,a,' sec')") timei(1:2), timei(3:4), timei(5:6)
        print*, ''
        print*, 'Start date: ',datei(7:8), '.',datei(5:6), '.',datei(1:4)
        print*, ''
    else
        call cpu_time(finish)
        print*, 'Start = ', start
        print*, 'Finish = ', finish
        if (finish - start < 60) then
            write(*,"(' Elapsed CPU time = ',f12.3,' seconds.')") finish-start
            print*, ''
        else if (finish - start < 3600) then
            write(*,"(' Elapsed CPU time = ',f12.3,' minutes.')") (finish-start)/60
            print*, ''
        else
            write(*,"(' Elapsed CPU time = ',f12.3,' hours.')") (finish-start)/3600
            print*, ''
        end if
        call date_and_time(date=datef,time=timef,values=valuesf)
        write(*,"(' Calculation started at',x,a,'h',x,a,'min',x,a,'sec')") timei(1:2), timei(3:4), timei(5:6)
        print*, ''
        write(*,"(' Calculation ended at',x,a,'h',x,a,'min',x,a,'sec')") timef(1:2), timef(3:4), timef(5:6)
        print*, ''
        print*, 'Start date: ',datei(7:8), '.',datei(5:6), '.',datei(1:4)
        print*, ''
        print*, 'End date: ',datef(7:8), '.',datef(5:6), '.',datef(1:4)

        file_name=dir//"times_"//params
        file_name=trim_name(file_name)
        open(77,file = file_name)
        write(77,"(' Elapsed CPU time = ',f20.10,' seconds.')") finish-start
        write(77,"(' Calculation started at',x,a,'h',x,a,'min',x,a,'sec')") timei(1:2), timei(3:4), timei(5:6)
        write(77,"(' Start date:',x,a,'.',x,a,'.',x,a)") ,datei(7:8), datei(5:6), datei(1:4)
        write(77,"(' Calculation ended at',x,a,'h',x,a,'min',x,a,'sec')") timef(1:2), timef(3:4), timef(5:6)
        write(77,"(' End date:',x,a,'.',x,a,'.',x,a)") ,datef(7:8), datef(5:6), datef(1:4)
        close(77)
    end if
end subroutine datetime

!---------------------------------------------------!
!            Number of discretization steps         !
!---------------------------------------------------!

subroutine nsteps(min, max, delta, steps)

    implicit none
    double precision, intent(in) :: min, max, delta
    integer, intent(out) :: steps

    if(delta <= 1.d0) then
        steps = int(abs(max-min)/delta + delta/2) + 1
    else
        steps = int(abs(max-min)/delta) + 1
    end if


end subroutine nsteps

!------------------------------------------!
!            Check parallelization         !
!------------------------------------------!

subroutine check_parallel()

    implicit none
    integer :: num_threads_loc
    integer :: thread_num_loc
    integer :: max_threads_loc

    num_threads_loc = 0
    thread_num_loc = 0
    max_threads_loc = 0
    !Critical block is a lock, forcing the instructions to be run in serial
    print*, ''
    !$omp parallel
        !$omp critical
            !$ thread_num_loc = omp_get_thread_num()
            print('(1x,100(a,i0))'), 'Thread number ', thread_num_loc, ' is online.'
            !$ num_threads_loc = omp_get_num_threads()
            if (thread_num_loc == 0) then
                print*, ''
                print('(1x,100(a,i0))'), 'Available number of threads = ', num_threads_loc
                !$ max_threads_loc = omp_get_max_threads()
                print('(1x,100(a,i0))'), 'Maximum number of threads = ', max_threads_loc
                print*, ''
            end if
        !$omp end critical
    !$omp end parallel
    print*, ''
    return 

end subroutine

!-------------------------------------!
!            Error messages           !
!-------------------------------------!

!Error messages for diagonalization routines
subroutine dsaupderrormessage(dsaupdinfo)
    implicit none
    integer :: dsaupdinfo

    if (dsaupdinfo .EQ. 0) THEN
    PRINT *, 'Normal Exit in dsaupd: info = 0.'
    elseif (dsaupdinfo .EQ. -1) THEN
    PRINT *, 'Error in dsaupd: info = -1.'
    PRINT *, 'N must be positive.'
    elseif (dsaupdinfo .EQ. -2) THEN
    PRINT *, 'Error in dsaupd: info = -2.'
    PRINT *, 'NEV must be positive.'
    elseif (dsaupdinfo .EQ. -3) THEN
    PRINT *, 'Error in dsaupd: info = -3.'
    PRINT *, 'NCV must be between NEV and N. '
    elseif (dsaupdinfo .EQ. -4) THEN
    PRINT *, 'Error in dsaupd: info = -4'
    PRINT *, 'The maximum number of Arnoldi update iterations allowed must be greater than zero.'
    elseif (dsaupdinfo .EQ. -5) THEN
    PRINT *, 'Error in dsaupd: info = -5'
    PRINT *, 'WHICH must be LM, SM, LA, SA, or BE. info = -5.'
    elseif (dsaupdinfo .EQ. -6) THEN
    PRINT *, 'Error in dsauupd: info = -6. '
    PRINT *, 'BMAT must be I or G. '
    elseif (dsaupdinfo .EQ. -7) THEN
    PRINT *, 'Error in dsaupd: info = -7.'
    PRINT *, 'Length of private work work WORKL array isnt sufficient.'
    elseif (dsaupdinfo .EQ. -8) THEN
    PRINT *, 'Error in dsaupd: info = -8.'
    PRINT *, 'Error in return from trid. eval calc. Error info from LAPACK dsteqr. info =-8'
    elseif (dsaupdinfo .EQ. -9) THEN
    PRINT *, 'Error in dsaupd: info = -9.'
    PRINT *, 'Starting vector is 0.'
    elseif (dsaupdinfo .EQ. -10) THEN
    PRINT *, 'Error in dsaupd: info = -10. '
    PRINT *, 'IPARAM(7) must be 1,2,3,4, or 5.'
    elseif (dsaupdinfo .EQ. -11) THEN
    PRINT *, 'Error in dsaupd: info = -11.'
    PRINT *, 'IPARAM(7)=1 and BMAT=G are incompatible.'
    elseif (dsaupdinfo .EQ. -12) THEN
    PRINT *, 'Error in dsaupd: info = -12'
    PRINT *, 'NEV and WHICH=BE are incompatible.'
    elseif (dsaupdinfo .EQ. -13) THEN
    PRINT *, 'Error in dsaupd: info = -13.'
    PRINT *, 'DSAUPD did find any eigenvalues to sufficient accuracy.'
    elseif (dsaupdinfo .EQ. -9999) THEN
    PRINT *, 'Error in dsaupd: info = -9999'
    PRINT *, 'Could not build an Arnoldi factorization. IPARAM(5) returns the size of the current Arnoldi factorization. &
    The user is advised to check that enough workspace and array storage     has been allocated. '
    elseif (dsaupdinfo .EQ. 1) THEN
    PRINT *, 'Error in dsaupd: info = 1'
    PRINT *, 'Maximum number of iterations taken. All possible eigenvalues of OP has been found. '
    PRINT *, 'IPARAM(5) returns the number of wanted converged Ritz values.'
    elseif (dsaupdinfo .EQ. 3) THEN
    PRINT *, 'Error in dsaupd: info =3'
    PRINT *, 'No shifts could be applied during a cycle of the Implicitly restarted Arnoldi iteration.'
    PRINT *, 'One possibility is to increase the size of NCV relative to NEV.'
    else
    PRINT *, 'Unknown error.  info =', dsaupdinfo
    END IF
end subroutine dsaupderrormessage

subroutine dseupderrormessage(dseupdinfo)
    implicit none
    integer :: dseupdinfo

    if (dseupdinfo .EQ. 0) THEN
    PRINT *, 'Normal Exit in dseupd: info = 0.'
    elseif (dseupdinfo .EQ. -1) THEN
    PRINT *, 'Error in deseupd: N must be positive. info =-1.'
    elseif (dseupdinfo .EQ. -2) THEN
    PRINT *, 'Error in deseupd: NEV must be positive. info = -2.'
    elseif (dseupdinfo .EQ. -3) THEN
    PRINT *, 'Error in deseupd: NCV must be between NEV and N. info = -3.'
    elseif (dseupdinfo .EQ. -5) THEN
    PRINT *, 'Error in deseupd: WHICH must be LM, SM, LA, SA, or BE info = -5.'
    elseif (dseupdinfo .EQ. -6) THEN
    PRINT *, 'Error in deseupd: BMAT must be I or G. info = -6.'
    elseif (dseupdinfo .EQ. -7) THEN
    PRINT *, 'Error in deseupd: N Length of private work work WORKL array isnt sufficient. info = -7.'
    elseif (dseupdinfo .EQ. -8) THEN
    PRINT *, 'Error in deseupd: Error in return from trid. eval calc. Error info from LAPACK dsteqr. info = -8.'
    elseif (dseupdinfo .EQ. -9) THEN
    PRINT *, 'Error in deseupd: Starting vector is 0. info = -9.'
    elseif (dseupdinfo .EQ. -10) THEN
    PRINT *, 'Error in deseupd: IPARAM(7) must be 1,2,3,4, or 5. info = -10.'
    elseif (dseupdinfo .EQ. -11) THEN
    PRINT *, 'Error in deseupd: IPARAM(7)=1 and BMAT=G are incompatible. info = -11.'
    elseif (dseupdinfo .EQ. -12) THEN
    PRINT *, 'Error in deseupd: NEV and WHICH=BE are incompatible. info = -12.'
    elseif (dseupdinfo .EQ. -14) THEN
    PRINT *, 'Error in deseupd: DSAUPD did find any eigenvalues to sufficient accuracy. info = -14.'
    elseif (dseupdinfo .EQ. -15) THEN
    PRINT *, 'Error in deseupd: HOWMANY must one A or S if RVEC=1. info = -15.'
    elseif (dseupdinfo .EQ. -16) THEN
    PRINT *, 'Error in deseupd: HOWMANY =S not yet implemented. info = -16.'
    elseif (dseupdinfo .EQ. -17) THEN
    PRINT *, 'Error in deseupd: info =-17.'
    PRINT *, 'DSEUPD got a different count of the number of converged Ritz values than DSAUPD.'
    PRINT *, 'User likely made an error in passing data DSAUPD -> DSEUPD. info = -17.'
    else
    PRINT *, 'Unknown error.  info =', dseupdinfo
    END IF

end subroutine dseupderrormessage

subroutine printparams(nev, ndeg, nHel, tilt, k1, k2)

    use variables
    implicit none
    integer, intent(in) :: nev, ndeg, nHel, tilt, k1, k2
    print*, ''
    print*, '!-------------------------------------------------------------!'
    print*, '!                                                             !'
    print*, '!                          Parameters                         !'
    print*, '!                                                             !'
    print*, '!-------------------------------------------------------------!'
    print*, ''
    print('(1x, a,i0)'), 'L = ', sites
    print('(1x, a,i0)'), 'N = ', particles
    if(tilted == 0) then 
    print('(1x, a,i0)'), 'Unit cells in x direction = ', ucx
    print('(1x, a,i0)'), 'Unit cells in y direction = ', ucy
    else if(tilted == 1) then 
        print('(1x, a, a)'), 'Cluster = ', cluster 
        print('(1x, a, i0)'),'Helices = ', nHel 
        print('(1x, a, i0)'),'Tilt    = ', tilt 
    end if 
    print*, ''
    print('(1x, a,i0)'), 'DIM = ', dim
    print('(1x, a,i0)'), 'NEV = ', nev

    if(ndeg > 0) then 
        if(degeneracy == 1) then 
            print*, ''
            ! print('(1x, a,f20.16)'), 'Degeneracy threshold = ', deg
            print('(1x, a,i0)'), 'NDEG = ', ndeg
            print*, ''
        else if(degeneracy == 2) then 
            print('(1x, a,f20.16)'), 'QNDEG = ', ndeg
            print*, ''
        end if 
    end if
    if(ti == 1) then 
        print('(1x, a,i0)'), 'K1   = ', k1
        print('(1x, a,i0)'), 'K2   = ', k2
        print('(1x, a,i0)'), 'Ktot = ', k1 + ucx * k2
        print*, ''    
    end if 
    if(symmetrize == 1) print('(1x, a,a)'),'IRREP = ', irrep 
    if(symmetrize == 1) print*, ''
    print*, '---------------------------------------------------------------'
    print*, ''    
    print('(1x, a,i0)'), 'OMP = ', othrds
    print('(1x, a,i0)'), 'MKL = ', mthrds
    if(exact == 1) then 
        print('(1x, a)'), 'EXACT '
        print*, ''
    end if 
    if(arpack == 1) then 
        print('(1x, a)'), 'ARPACK '
        print*, ''
    end if 
    if(mkl == 1) then 
        print('(1x, a)'), 'MKL '
        print*, ''
    end if 
    if(feast == 1) then 
        print('(1x, a)'), 'FEAST '
        print*, ''
    end if 
    if(curr == 1) then 
        print('(1x, a)'), 'CURRENTS '
        print*, ''
    end if 

    print*, '---------------------------------------------------------------'
    print*, ''
    print('(1x, a)'), 'DISORDER'
    print('(1x, a,f6.4)'), 'W    = ', dis
    print('(1x, a,i0)'),   'Ndis = ', nDis
    print*, ''
    print*, '---------------------------------------------------------------'
    print*, ''
    print('(1x, a,f6.4)'), 'dv   = ',dv
    print('(1x, a,f6.4)'), 'Vmin = ', vmin
    print('(1x, a,f6.4)'), 'Vmax = ', vmax
    print*, ''
    print*, '---------------------------------------------------------------'
    print*, ''
    print('(1x, a,f6.4)'), 'dv2   = ',dv2
    print('(1x, a,f6.4)'), 'V2min = ', v2min
    print('(1x, a,f6.4)'), 'V2max = ', v2max
    print*, ''
    print*, '---------------------------------------------------------------'    
    print*, ''
    if(mass > 0.d0) print('(1x, a,f6.4)'), 'Mass = ', mass
    if(bc == 'p') then
        print*, ''
        print('(1x, a)'), 'PBC'
        print*, ''
    elseif (bc == 'o') then
        print*, ''
        print('(1x, a)'), 'OBC'
        print*, ''
    elseif (bc == 't') then
        print*, ''
        print('(1x, a)'), 'TBC'
        print*, ''
    end if
    print*, '!-------------------------------------------------------------!'
    print*, '!-------------------------------------------------------------!'
    print*, ''

end subroutine printparams




















!--------------------------------------------!
!           Older current routines           !
!--------------------------------------------!


subroutine current_c(gencluster, nHel, tilt, k1, k2, t, dim, sites, Lx, Ly, basis, bsites, xy, xtransl, ytransl, nbonds, refbond, psi, nnnVec, norm, bondcurrent, current)

    implicit none
    integer, intent(in) :: gencluster, nHel, tilt, k1, k2, nbonds, refbond, sites, Lx, Ly
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    integer, intent(in) :: xy(4,nbonds)
    integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
    double precision, intent(in) :: t 
    double precision, intent(in) :: nnnVec(2, 3)
    double complex, intent(in) :: psi(dim)
    double precision, intent(in) :: norm(dim)
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bondcurrent(:)
    
    integer(kind=8) :: loc = 0, newst = 0, rep = 0, cntr = 0 
    integer :: i = 0, j = 0, l1 = 0, l2 = 0, l = 0, l11 = 0, l22 = 0 
    integer :: x = 0, y = 0
    integer :: sign = 0, signt = 0 
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double complex, allocatable :: psiprime(:)
    character*250 :: name, name2 

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)

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
        psiprime = 0.0d0  
        site1 = bsites(1,i)
        site2 = bsites(2,i)
        refsite1 = bsites(1,refbond)
        refsite2 = bsites(2,refbond)    
        phase = bsites(3,i)   
        vec   = vecs(1:2,bsites(5,i))
        dist  = dot_product(vec, refvec)

        do x = 0, Lx - 1
            do y = 0, Ly - 1
                if (site1 == refsite1 .or. site2 == refsite2 &
                    .or. site1 == refsite2 .or. site2 == refsite1) go to 110 !Exclude bonds which share a site with the reference bond

                do j = 1, dim !Loop over basis states 
                    if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                        newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                        if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                            ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                            newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                            if(gencluster == 0) then 
                                call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                l11 = Ly
                                l22 = Lx
                            else if(gencluster == 1) then 
                                call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                l11 = Lx
                                l22 = Ly
                            end if
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(j) = psiprime(j) + sign * signt * (-(t)**2) * sqrt(norm(loc)/norm(j)) * psi(loc) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                        else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                            ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                            newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                            if(gencluster == 0) then 
                                call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                l11 = Ly
                                l22 = Lx
                            else if(gencluster == 1) then 
                                call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                l11 = Lx
                                l22 = Ly
                            end if
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(j) = psiprime(j) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(loc) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                        else 
                            cycle 
                        end if     
                    else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                        newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                        if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                            ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                            newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                            if(gencluster == 0) then 
                                call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                l11 = Ly
                                l22 = Lx
                            else if(gencluster == 1) then 
                                call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                l11 = Lx
                                l22 = Ly
                            end if
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            
                            parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(j) = psiprime(j) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(loc) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                        else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                            ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                            newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                            if(gencluster == 0) then 
                                call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                l11 = Ly
                                l22 = Lx
                            else if(gencluster == 1) then 
                                call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                l11 = Lx
                                l22 = Ly
                            end if
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            
                            parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
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

end subroutine current_c

subroutine currmat_c(gencluster, nHel, tilt, dir, unit, filename, sl, ndeg, k1, k2, t, dim, sites, Lx, Ly, basis, bsites, xy, xtransl, ytransl, nbonds, refbond, psi, nnnVec, norm)

    implicit none
    integer, intent(in) :: gencluster, nHel, tilt, unit, ndeg, k1, k2, nbonds, refbond, sites, Lx, Ly  
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    integer, intent(in) :: xy(4,nbonds)
    integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
    double precision, intent(in) :: t 
    double precision, intent(in) :: nnnVec(2, 3)
    double complex, intent(in) :: psi(dim, ndeg)
    double precision, intent(in) :: norm(dim)
    character(len=*), intent(in) :: dir, sl, filename
    
    integer(kind=8) :: loc = 0, newst = 0, rep = 0
    integer :: i = 0, j = 0, k = 0, l = 0, l1 = 0, l2 = 0
    integer :: sign = 0, signt = 0 
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0, current= 0.d0 
    double precision, dimension(2)   :: vec, refvec, coord1, coord2 
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: currents(:)
    double precision, allocatable :: currents2(:)
    double precision, allocatable :: evals(:)
    double precision, allocatable :: bondcurrents(:,:), bondcurrent(:)
    double precision, allocatable :: currentmat(:,:)
    double precision, allocatable :: bondcurrentmat(:,:,:)
    double complex, allocatable :: psiprime(:), gs(:)
    character*256 :: name 

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(psiprime)) deallocate(psiprime)
    if(allocated(bondcurrentmat)) deallocate(bondcurrentmat)
    if(allocated(currentmat)) deallocate(currentmat)
    if(allocated(currents2)) deallocate(currents2)
    if(allocated(bondcurrents)) deallocate(bondcurrents)
    allocate(psiprime(dim))
    allocate(bondcurrentmat(nbonds, ndeg, ndeg))
    allocate(currentmat(ndeg, ndeg))
    allocate(currents2(ndeg))
    allocate(bondcurrents(nbonds, ndeg))

    psiprime = 0.d0
    refvec   = vecs(1:2,bsites(5,refbond))
    coord1(1:2) = int(0.5*xy(3,refbond))*vecs(1:2,1) + xy(4,refbond)*vecs(1:2,2) - (int(0.5*xy(1,refbond))*vecs(1:2,1) + xy(2,refbond)*vecs(1:2,2))
    refsite1          = bsites(1,refbond)
    refsite2          = bsites(2,refbond)      

    ! do i = 1, ndeg
        !     do j = 1, ndeg
        !         if ( abs(dot_product(psi(1:dim,i),psi(1:dim,j))) > 10.**(-10) .and. abs(dot_product(psi(1:dim,i),psi(1:dim,j))) < 0.999999999d0 ) then 
        !             print*,'currmat_c: Groundstates ',i, 'and' ,j,' are not orthorgonal.', dot_product(psi(1:dim,i),psi(1:dim,j))               
        !             error stop 
        !         end if 
        !     end do 
        ! end do 
    
    currentmat = 0.d0     
    bondcurrentmat = 0.d0  
    do l = 1, ndeg 
        do k = 1, ndeg 
            do i = 1, nbonds                            
                site1 = bsites(1,i)
                site2 = bsites(2,i)
                if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
                    .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
                coord2(1:2) = int(0.5*xy(3,i))*vecs(1:2,1) + xy(4,i)*vecs(1:2,2) - (int(0.5*xy(1,i))*vecs(1:2,1) + xy(2,i)*vecs(1:2,2))
                ! coord2(1:2) = coord2(1:2)/norm2(coord2(1:2)) 
                phase = bsites(3,i)   
                vec   = vecs(1:2,bsites(5,i))
                dist = dot_product(vec, refvec)
                psiprime = 0 

                do j = 1, dim
                    if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                        newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                        if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                            ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                            newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                            call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)  
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(Lx)) + (dble(k2*l2)/dble(Ly))))
                        else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                            ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                            newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                            call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)  
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(Lx)) + (dble(k2*l2)/dble(Ly))))
                        end if     
                    else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                        newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                        if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                            ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                            newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                            call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)  
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(Lx)) + (dble(k2*l2)/dble(Ly))))
                        else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                            ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                            newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                            call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)  
                            call findstate(dim, rep, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(Lx)) + (dble(k2*l2)/dble(Ly))))
                        end if     
                    end if 
                end do !Loop over basis states                  
                bondcurrentmat(i,k,l) = bondcurrentmat(i,k,l) + dble(phase) * dist * dot_product(psi(1:dim,k), psiprime)
                currentmat(k,l) = currentmat(k,l) + bondcurrentmat(i,k,l)
                ! bondcurrentmat(i,k,l) = bondcurrentmat(i,k,l) + dist * dot_product(psi(1:dim,k), psiprime)
                ! currentmat(k,l) = currentmat(k,l) + dble(phase) * bondcurrentmat(i,k,l)

            end do !Loop over bonds 
        end do !Loop over degenerate ground states
    end do !Loop over degenerate ground states
    deallocate(psiprime)

    call exactdiag(.True., ndeg, currentmat, currents)

    ! print*,'Current eigenvalues'
        ! do i = 1, ndeg 
        !     print*, currents(i)
        ! end do
        ! print*,'sum(currents)',sum(currents) 
    
        ! do i = 1, nbonds
        !     call exactdiag(.True., ndeg, bondcurrentmat(i, 1:ndeg, 1:ndeg), evals)
        !     bondcurrents(i, 1:ndeg) = evals
        !     if(allocated(evals)) deallocate(evals)
        ! end do 
        
        
        ! if(sl == "A") then 
        !     name = dir//"Acurrent_dens_mat_"//filename
        ! else if(sl == "B") then 
        !     name = dir//"Bcurrent_dens_mat_"//filename
        ! end if 
        ! name = trim_name(name)
        ! open(unit,file=name)
        ! do i = 1, ndeg 
        !     currents2(i) = sum(bondcurrents(1:nbonds, i))
        !     write(unit,*) currents2(i)
        !     ! print*,'Current of psi ',i, currents2(i)
        ! end do 
        ! close(unit)

    currentmat = currentmat / ndeg 


    call exactdiag(.True., ndeg, currentmat, currents)
    ! print*,'|===================|'
        ! print*,' Current eigenvalues'
        ! do i = 1, ndeg 
        !     print*, currents(i)
        ! end do
        ! print*,'|===================|'
        ! print*,' Total'
        ! print*,sum(currents) 
        ! currents = 0.d0 
        ! print*,'|===================|'
    if(allocated(gs)) deallocate(gs)
    allocate(gs(dim))

    do j = 1, ndeg 
        gs = 0.d0 
        do i = 1, ndeg 
            gs = gs + currentmat(i,j) * psi(1:dim, i)
        end do
        ! print*,'gs',gs
        call current_c(gencluster, nHel, tilt, k1, k2, t, dim, sites, Lx, Ly, basis, bsites, xy, xtransl, ytransl, nbonds, refbond, gs, nnnVec, norm, bondcurrent, current)
        bondcurrents(1:nbonds, j) = bondcurrent 
        currents(j) = current 
        ! print*,'|===================|'
            ! print*,' Current, eigenstate'
            ! print*, current, j
            ! print*,'|===================|'
        current = 0.d0 
    end do 
    ! print*,'|===================|'
        ! print*,' Sum'
        ! print*, sum(currents)
        ! print*,'|===================|'


    if(sl == "A") then 
        name = trim_name(dir//"Acurrent_matrix_"//filename)
    else if(sl == "B") then 
        name = trim_name(dir//"Bcurrent_matrix_"//filename)
    end if 
    name = trim_name(name)
    
    open(unit,file=name)
    do i = 1, ndeg 
        write(unit,*) currents(i)
    end do 
    close(unit)
    if(sl == "A") then 
        name = dir//"Abondcurrent_matrix_"//filename
    else if(sl == "B") then 
        name = dir//"Bbondcurrent_matrix_"//filename
    end if 
    name = trim_name(name)
    open(unit,file=name)
    do i = 1, ndeg 
        write(unit,*) bondcurrents(1:nbonds,i)
    end do 
    close(unit)

    ! bondcurrents = 0.d0 
        ! do i = 1, nbonds
        !     call exactdiag(.True., ndeg, bondcurrentmat(i, 1:ndeg, 1:ndeg), evals)
        !     bondcurrents(i, 1:ndeg) = evals
        !     if(allocated(evals)) deallocate(evals)
        ! end do 
        ! if(sl == "A") then 
        !     name = dir//"Acurrent_matrix2_"//filename
        ! else if(sl == "B") then 
        !     name = dir//"Bcurrent_matrix2_"//filename
        ! end if 
        ! name = trim_name(name)
        ! open(unit,file=name)
        ! do i = 1, ndeg 
        !     currents2(i) = sum(bondcurrents(1:nbonds, i))
        !     write(unit,*) currents2(i)
        !     print*,'Current of psi ',i, currents2(i)
        ! end do 
        ! close(unit)

    return

end subroutine currmat_c

subroutine currmat_c2(gencluster, nHel, tilt, dir, unit, filename, sl, ndeg, k1, k2, t, dim, sites, Lx, Ly, basis, bsites, xy, xtransl, ytransl, nbonds, refbond, psi, nnnVec, norm)

    implicit none
    integer, intent(in) :: gencluster, nHel, tilt, unit, ndeg, k1, k2, nbonds, refbond, sites, Lx, Ly  
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    integer, intent(in) :: xy(4,nbonds)
    integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
    double precision, intent(in) :: t 
    double precision, intent(in) :: nnnVec(2, 3)
    double complex, intent(in) :: psi(dim, ndeg)
    double precision, intent(in) :: norm(dim)
    character(len=*), intent(in) :: dir, sl, filename
    
    integer(kind=8) :: loc = 0, newst = 0, rep = 0
    integer :: i = 0, j = 0, k = 0, l = 0, l1 = 0, l2 = 0, x = 0, y = 0, cntr = 0, q = 0 
    integer :: sign = 0, signt = 0, l11 = 0, l22 = 0 
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0, current= 0.d0 
    double precision, dimension(2)   :: vec, refvec, coord1, coord2 
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: currents(:)
    double precision, allocatable :: currents2(:)
    double precision, allocatable :: bondcurrents(:,:), bondcurrent(:)
    double precision, allocatable :: currentmat(:,:)
    double precision, allocatable :: bondcurrentmat(:,:,:)
    double complex, allocatable :: psiprime(:), gs(:)
    character*256 :: name 

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(psiprime)) deallocate(psiprime)
    if(allocated(currents2)) deallocate(currents2)
    if(allocated(currentmat)) deallocate(currentmat)
    if(allocated(bondcurrents)) deallocate(bondcurrents)
    if(allocated(bondcurrentmat)) deallocate(bondcurrentmat)
    allocate(psiprime(dim))
    allocate(currents2(ndeg))
    allocate(currentmat(ndeg,ndeg))
    allocate(bondcurrents(nbonds,ndeg))
    allocate(bondcurrentmat(nbonds,ndeg,ndeg))

    psiprime = 0.d0
    refvec   = vecs(1:2,bsites(5,refbond))
    coord1(1:2) = int(0.5*xy(3,refbond))*vecs(1:2,1) + xy(4,refbond)*vecs(1:2,2) - (int(0.5*xy(1,refbond))*vecs(1:2,1) + xy(2,refbond)*vecs(1:2,2))
    refsite1          = bsites(1,refbond)
    refsite2          = bsites(2,refbond)      

    ! do i = 1, ndeg
        !     do j = 1, ndeg
        !         if ( abs(dot_product(psi(1:dim,i),psi(1:dim,j))) > 10.**(-10) .and. abs(dot_product(psi(1:dim,i),psi(1:dim,j))) < 0.999999999d0 ) then 
        !             print*,'currmat_c: Groundstates ',i, 'and' ,j,' are not orthorgonal.', dot_product(psi(1:dim,i),psi(1:dim,j))               
        !             error stop 
        !         end if 
        !     end do 
        ! end do 
        
    currentmat = 0.d0     
    bondcurrentmat = 0.d0  
    do l = 1, ndeg 
        do k = 1, ndeg 
            do i = 1, nbonds                            
                site1 = bsites(1,i)
                site2 = bsites(2,i)
                if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
                    .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
                coord2(1:2) = int(0.5*xy(3,i))*vecs(1:2,1) + xy(4,i)*vecs(1:2,2) - (int(0.5*xy(1,i))*vecs(1:2,1) + xy(2,i)*vecs(1:2,2))
                ! coord2(1:2) = coord2(1:2)/norm2(coord2(1:2)) 
                phase = bsites(3,i)   
                vec   = vecs(1:2,bsites(5,i))
                dist = dot_product(vec, refvec)
                psiprime = 0 

                do x = 0, Lx - 1
                    do y = 0, Ly - 1
                        do j = 1, dim
                            if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                                newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                                if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                                    ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                                    newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                                    if(gencluster == 0) then 
                                        call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                        l11 = Ly
                                        l22 = Lx
                                    else if(gencluster == 1) then 
                                        call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                        l11 = Lx
                                        l22 = Ly
                                    end if
                                    call findstate(dim, rep, basis, loc) 
                                    if(loc <= 0) cycle   
                                    parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                                    parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                    psiprime(loc) = psiprime(loc) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                                    ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                                    newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                                    if(gencluster == 0) then 
                                        call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                        l11 = Ly
                                        l22 = Lx
                                    else if(gencluster == 1) then 
                                        call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                        l11 = Lx
                                        l22 = Ly
                                    end if
                                    call findstate(dim, rep, basis, loc) 
                                    if(loc <= 0) cycle   
                                    parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                                    parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                    psiprime(loc) = psiprime(loc) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                end if     
                            else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                                newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                                if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                                    ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                                    newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                                    if(gencluster == 0) then 
                                        call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                        l11 = Ly
                                        l22 = Lx
                                    else if(gencluster == 1) then 
                                        call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                        l11 = Lx
                                        l22 = Ly
                                    end if
                                    call findstate(dim, rep, basis, loc) 
                                    if(loc <= 0) cycle   
                                    parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                                    parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                    psiprime(loc) = psiprime(loc) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                                    ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                                    newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                                    if(gencluster == 0) then 
                                        call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                        l11 = Ly
                                        l22 = Lx
                                    else if(gencluster == 1) then 
                                        call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                        l11 = Lx
                                        l22 = Ly
                                    end if
                                    call findstate(dim, rep, basis, loc) 
                                    if(loc <= 0) cycle   
                                    parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                                    parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                    psiprime(loc) = psiprime(loc) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                end if     
                            end if 
                        end do !Loop over basis states                  

                        !a1 translations 
                        cntr = 0
                        do q = 1, nbonds
                            if(bsites(1,q) == site1 .and. bsites(5,q) == 1) then 
                                site1 = bsites(2,q)
                                cntr = cntr + 1
                            else if(bsites(1,q) == site2 .and. bsites(5,q) == 1) then
                                site2 = bsites(2,q)
                                cntr = cntr + 1
                            else if(bsites(1,q) == refsite1 .and. bsites(5,q) == 1) then
                                refsite1 = bsites(2,q)
                                cntr = cntr + 1
                            else if(bsites(1,q) == refsite2 .and. bsites(5,q) == 1) then
                                refsite2 = bsites(2,q)
                                cntr = cntr + 1
                            else if(cntr == 4) then 
                                exit 
                            end if
                        end do 

                    end do 

                    !a2 translations 
                    cntr = 0
                    do q = 1, nbonds
                        if(bsites(1,q) == site1 .and. bsites(5,q) == 2) then 
                            site1 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == site2 .and. bsites(5,q) == 2) then
                            site2 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == refsite1 .and. bsites(5,q) == 2) then
                            refsite1 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == refsite2 .and. bsites(5,q) == 2) then
                            refsite2 = bsites(2,q)
                            cntr = cntr + 1
                        else if(cntr == 4) then 
                            exit 
                        end if
                    end do 
                end do 
                bondcurrentmat(i,k,l) = bondcurrentmat(i,k,l) + dble(phase) * dist * dot_product(psi(1:dim,k), psiprime)
                currentmat(k,l) = currentmat(k,l) + bondcurrentmat(i,k,l) / (Lx*Ly)
                ! bondcurrentmat(i,k,l) = bondcurrentmat(i,k,l) + dist * dot_product(psi(1:dim,k), psiprime)
                ! currentmat(k,l) = currentmat(k,l) + dble(phase) * bondcurrentmat(i,k,l)

            end do !Loop over bonds 
        end do !Loop over degenerate ground states
    end do !Loop over degenerate ground states
    deallocate(psiprime)

    call exactdiag(.True., ndeg, currentmat, currents)

    ! print*,'Current eigenvalues'
        ! do i = 1, ndeg 
        !     print*, currents(i)
        ! end do
        ! print*,'sum(currents)',sum(currents) 
        ! do i = 1, nbonds
        !     call exactdiag(.True., ndeg, bondcurrentmat(i, 1:ndeg, 1:ndeg), evals)
        !     bondcurrents(i, 1:ndeg) = evals
        !     if(allocated(evals)) deallocate(evals)
        ! end do     
        ! if(sl == "A") then 
        !     name = dir//"Acurrent_dens_mat_"//filename
        ! else if(sl == "B") then 
        !     name = dir//"Bcurrent_dens_mat_"//filename
        ! end if 
        ! name = trim_name(name)
        ! open(unit,file=name)
        ! do i = 1, ndeg 
        !     currents2(i) = sum(bondcurrents(1:nbonds, i))
        !     write(unit,*) currents2(i)
        !     ! print*,'Current of psi ',i, currents2(i)
        ! end do 
        ! close(unit)

    currentmat = currentmat / ndeg 


    call exactdiag(.True., ndeg, currentmat, currents)
    ! print*,'|===================|'
        ! print*,' Current eigenvalues'
        ! do i = 1, ndeg 
        !     print*, currents(i)
        ! end do
        ! print*,'|===================|'
        ! print*,' Total'
        ! print*,sum(currents) 
        ! currents = 0.d0 
        ! print*,'|===================|'
    if(allocated(gs)) deallocate(gs)
    allocate(gs(dim))

    do j = 1, ndeg 
        gs = 0.d0 
        do i = 1, ndeg            
            gs = gs + currentmat(i,j) * psi(1:dim, i)
        end do
        call current_c(gencluster, nHel, tilt, k1, k2, t, dim, sites, Lx, Ly, basis, bsites, xy, xtransl, ytransl, nbonds, refbond, gs, nnnVec, norm, bondcurrent, current)
        bondcurrents(1:nbonds, j) = bondcurrent 
        currents(j) = current 
        ! print*,'|===================|'
            ! print*,' Current, eigenstate'
            ! print*, current, j
            ! print*,'|===================|'
        current = 0.d0 
    end do 
    ! print*,'|===================|'
        ! print*,' Sum'
        ! print*, sum(currents)
        ! print*,'|===================|'


    if(sl == "A") then 
        name = trim_name(dir//"Acurrent_matrix_"//filename)
    else if(sl == "B") then 
        name = trim_name(dir//"Bcurrent_matrix_"//filename)
    end if 
    name = trim_name(name)
    print* ,name, 'name'
    open(unit,file=name)
    do i = 1, ndeg 
        write(unit,*) currents(i)
    end do 
    close(unit)
    if(sl == "A") then 
        name = trim_name(dir//"Abondcurrent_matrix_"//filename)
    else if(sl == "B") then 
        name = trim_name(dir//"Bbondcurrent_matrix_"//filename)
    end if 
    name = trim_name(name)
    open(unit,file=name)
    do i = 1, ndeg 
        write(unit,*) bondcurrents(1:nbonds,i)
    end do 
    close(unit)


    ! bondcurrents = 0.d0 
        ! do i = 1, nbonds
        !     call exactdiag(.True., ndeg, bondcurrentmat(i, 1:ndeg, 1:ndeg), evals)
        !     bondcurrents(i, 1:ndeg) = evals
        !     if(allocated(evals)) deallocate(evals)
        ! end do 
            
        ! if(sl == "A") then 
        !     name = dir//"Acurrent_matrix2_"//filename
        ! else if(sl == "B") then 
        !     name = dir//"Bcurrent_matrix2_"//filename
        ! end if 
        ! name = trim_name(name)
        ! open(unit,file=name)
        ! do i = 1, ndeg 
        !     currents2(i) = sum(bondcurrents(1:nbonds, i))
        !     write(unit,*) currents2(i)
        !     print*,'Current of psi ',i, currents2(i)
        ! end do 
        ! close(unit)

    return

end subroutine currmat_c2

subroutine current_ics_c2(gencluster, nHel, tilt, k1, k2, ndeg, t, dim, sites, Lx, Ly, basis, bsites, xy, xtransl, ytransl, nbonds, refbond, psi, nnnVec, norm, bondcurrent, current)

    implicit none
    integer, intent(in) :: gencluster, nHel, tilt, k1, k2, ndeg, nbonds, refbond, sites, Lx, Ly  
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    integer, intent(in) :: xy(4,nbonds)
    integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
    double precision, intent(in) :: t, nnnVec(2, 3) 
    double complex, intent(in) :: psi(dim, ndeg)
    double precision, intent(in) :: norm(dim)
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bondcurrent(:)
    
    integer(kind=8) :: loc = 0, newst = 0, rep = 0
    integer :: i = 0, j = 0, l = 0, l1 = 0, l2 = 0, x = 0, y = 0, q = 0, cntr = 0
    integer :: sign = 0, signt = 0, l11 = 0, l22 = 0  
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec, coord1, coord2 
    double precision, dimension(2,3) :: vecs
    double complex, allocatable :: psiprime(:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bondcurrent)) deallocate(bondcurrent)
    allocate(bondcurrent(nbonds))
    if( allocated(psiprime) ) deallocate(psiprime)
    allocate( psiprime(dim) )

    psiprime = 0.d0
    refvec   = vecs(1:2,bsites(5,refbond))
    coord1(1:2)    = int(0.5*xy(3,refbond))*vecs(1:2,1) + xy(4,refbond)*vecs(1:2,2) - (int(0.5*xy(1,refbond))*vecs(1:2,1) + xy(2,refbond)*vecs(1:2,2))
    refsite1          = bsites(1,refbond)
    refsite2          = bsites(2,refbond)      
   
    ! do i = 1, ndeg
        !     do j = 1, ndeg
        !         if ( abs(dot_product(psi(1:dim,i),psi(1:dim,j))) > 10.**(-10) .and. abs(dot_product(psi(1:dim,i),psi(1:dim,j))) < 0.999999999d0 ) then 
        !             print*,'current_ics_c: Groundstates ',i, 'and' ,j,' are not orthorgonal.', dot_product(psi(1:dim,i),psi(1:dim,j))               
        !             error stop 
        !         end if 
        !     end do 
        ! end do 

    bondcurrent = 0.d0 
    current = 0.d0     
    do i = 1, nbonds            
        bondcurrent(i) = 0 
        site1 = bsites(1,i)
        site2 = bsites(2,i)
        if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
            .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
        coord2(1:2) = int(0.5*xy(3,i))*vecs(1:2,1) + xy(4,i)*vecs(1:2,2) - (int(0.5*xy(1,i))*vecs(1:2,1) + xy(2,i)*vecs(1:2,2))
        phase = bsites(3,i)   
        vec   = vecs(1:2,bsites(5,i))
        dist = dot_product(vec, refvec)

        do l = 1, ndeg 
            psiprime = 0 
            do x = 0, Lx - 1
                do y = 0, Ly - 1
                    do j = 1, dim
                        if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                            newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                            if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                                ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                                newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                                if(gencluster == 0) then 
                                    call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                    l11 = Ly
                                    l22 = Lx
                                else if(gencluster == 1) then 
                                    call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                    l11 = Lx
                                    l22 = Ly
                                end if
                                call findstate(dim, rep, basis, loc) 
                                if(loc <= 0) cycle   
                                
                                parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                                parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                                parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                                parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                                sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                psiprime(loc) = psiprime(loc) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                            else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                                ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                                newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                                if(gencluster == 0) then 
                                    call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                    l11 = Ly
                                    l22 = Lx
                                else if(gencluster == 1) then 
                                    call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                    l11 = Lx
                                    l22 = Ly
                                end if
                                call findstate(dim, rep, basis, loc) 
                                if(loc <= 0) cycle   
                                
                                parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                                parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                                parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                                parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                                sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                psiprime(loc) = psiprime(loc) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                            end if     
                        else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                            newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                            if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                                ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                                newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                                if(gencluster == 0) then 
                                    call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                    l11 = Ly
                                    l22 = Lx
                                else if(gencluster == 1) then 
                                    call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                    l11 = Lx
                                    l22 = Ly
                                end if
                                call findstate(dim, rep, basis, loc) 
                                if(loc <= 0) cycle   
                                
                                parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                                parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                                parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                                parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                                sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                psiprime(loc) = psiprime(loc) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                            else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                                ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                                newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                                if(gencluster == 0) then 
                                    call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                    l11 = Ly
                                    l22 = Lx
                                else if(gencluster == 1) then 
                                    call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                    l11 = Lx
                                    l22 = Ly
                                end if
                                call findstate(dim, rep, basis, loc) 
                                if(loc <= 0) cycle   
                                
                                parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                                parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                                parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                                parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                                sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                psiprime(loc) = psiprime(loc) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                            end if     
                        end if 
                    end do !Loop over basis states      
                    
                    !a1 translations 
                    cntr = 0
                    do q = 1, nbonds
                        if(bsites(1,q) == site1 .and. bsites(5,q) == 1) then 
                            site1 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == site2 .and. bsites(5,q) == 1) then
                            site2 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == refsite1 .and. bsites(5,q) == 1) then
                            refsite1 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == refsite2 .and. bsites(5,q) == 1) then
                            refsite2 = bsites(2,q)
                            cntr = cntr + 1
                        else if(cntr == 4) then 
                            exit 
                        end if
                    end do 

                end do 

                !a2 translations 
                cntr = 0
                do q = 1, nbonds
                    if(bsites(1,q) == site1 .and. bsites(5,q) == 2) then 
                        site1 = bsites(2,q)
                        cntr = cntr + 1
                    else if(bsites(1,q) == site2 .and. bsites(5,q) == 2) then
                        site2 = bsites(2,q)
                        cntr = cntr + 1
                    else if(bsites(1,q) == refsite1 .and. bsites(5,q) == 2) then
                        refsite1 = bsites(2,q)
                        cntr = cntr + 1
                    else if(bsites(1,q) == refsite2 .and. bsites(5,q) == 2) then
                        refsite2 = bsites(2,q)
                        cntr = cntr + 1
                    else if(cntr == 4) then 
                        exit 
                    end if
                end do 
            end do 
            
            bondcurrent(i) = bondcurrent(i) + dist * dot_product(psi(1:dim,l), psiprime)
        end do !Loop over degenerate ground states
        current = current + dble(phase) * bondcurrent(i) / (Lx*Ly)
    end do !Loop over bonds 

    current = current/dble(ndeg)!sqrt(dble(ndeg))
    deallocate(psiprime)
 
    
    return

end subroutine current_ics_c2

subroutine current_cs_c2(gencluster, nHel, tilt, dir, unit, parameters, sl, ndeg, k1, k2, t, dim, sites, Lx, Ly, basis, bsites, xy, xtransl, ytransl, nbonds, refbond, psi, nnnVec, norm)

    implicit none
    integer, intent(in) :: gencluster, nHel, tilt, unit, ndeg, k1, k2, nbonds, refbond, sites, Lx, Ly  
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    integer, intent(in) :: xy(4,nbonds)
    integer, intent(in) :: xtransl(2, sites), ytransl(2, sites)
    double precision, intent(in) :: t, nnnVec(2, 3)
    double complex, intent(in) :: psi(dim, ndeg)
    double precision, intent(in) :: norm(dim)
    character(len=*), intent(in) :: dir, sl, parameters
    
    integer(kind=8) :: loc = 0, newst = 0, rep = 0
    integer :: i = 0, j = 0, k = 0, l = 0, l1 = 0, l2 = 0, x = 0, y = 0, q = 0, cntr = 0
    integer :: sign = 0, signt = 0, l11 = 0, l22 = 0  
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec, coord1, coord2 
    double precision, dimension(2,3) :: vecs
    double precision :: current
    double precision, allocatable :: bondcurrent(:)
    double complex, allocatable :: psiprime(:)
    character :: filename*512 
    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)

    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bondcurrent)) deallocate(bondcurrent)
    allocate(bondcurrent(nbonds))
    if( allocated(psiprime) ) deallocate(psiprime)
    allocate( psiprime(dim) )

    psiprime = 0.d0
    refvec   = vecs(1:2,bsites(5,refbond))
    coord1(1:2)    = int(0.5*xy(3,refbond))*vecs(1:2,1) + xy(4,refbond)*vecs(1:2,2) - (int(0.5*xy(1,refbond))*vecs(1:2,1) + xy(2,refbond)*vecs(1:2,2))
    refsite1          = bsites(1,refbond)
    refsite2          = bsites(2,refbond)      

    ! do i = 1, ndeg
        !     do j = 1, ndeg
        !         if ( abs(dot_product(psi(1:dim,i),psi(1:dim,j))) > 10.**(-10) .and. abs(dot_product(psi(1:dim,i),psi(1:dim,j))) < 0.999999999d0 ) then 
        !             print*,'current_cs_c: Groundstates ',i, 'and' ,j,' are not orthorgonal.', dot_product(psi(1:dim,i),psi(1:dim,j))               
        !             error stop 
        !         end if 
        !     end do 
        ! end do 
    bondcurrent = 0.d0
    current = 0.d0     
    do i = 1, nbonds            
        bondcurrent(i) = 0 
        site1 = bsites(1,i)
        site2 = bsites(2,i)
        if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
            .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
        coord2(1:2) = int(0.5*xy(3,i))*vecs(1:2,1) + xy(4,i)*vecs(1:2,2) - (int(0.5*xy(1,i))*vecs(1:2,1) + xy(2,i)*vecs(1:2,2))
        ! coord2(1:2) = coord2(1:2)/norm2(coord2(1:2)) 
        phase = bsites(3,i)   
        vec   = vecs(1:2,bsites(5,i))
        dist = dot_product(vec, refvec)

        do l = 1, ndeg 
            do k = 1, ndeg 
                psiprime = 0 
                do x = 0, Lx - 1
                    do y = 0, Ly - 1
                        do j = 1, dim
                            if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                                newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                                if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                                    ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                                    newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                                    if(gencluster == 0) then 
                                        call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                        l11 = Ly
                                        l22 = Lx
                                    else if(gencluster == 1) then 
                                        call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                        l11 = Lx
                                        l22 = Ly
                                    end if  
                                    call findstate(dim, rep, basis, loc) 
                                    if(loc <= 0) cycle   
                                    parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                                    parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                    psiprime(loc) = psiprime(loc) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                                    ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                                    newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                                    if(gencluster == 0) then 
                                            call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                            l11 = Ly
                                            l22 = Lx
                                        else if(gencluster == 1) then 
                                            call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                            l11 = Lx
                                            l22 = Ly
                                        end if  
                                    call findstate(dim, rep, basis, loc) 
                                    if(loc <= 0) cycle   
                                    parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                                    parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                    psiprime(loc) = psiprime(loc) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                end if     
                            else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                                newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                                if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                                    ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                                    newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                                    if(gencluster == 0) then 
                                        call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                        l11 = Ly
                                        l22 = Lx
                                    else if(gencluster == 1) then 
                                        call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                        l11 = Lx
                                        l22 = Ly
                                    end if  
                                    call findstate(dim, rep, basis, loc) 
                                    if(loc <= 0) cycle   
                                    parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                                    parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                    psiprime(loc) = psiprime(loc) + (-1) * sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                                    ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                                    newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                                    if(gencluster == 0) then 
                                            call representative_tilted(newst, sites, nHel, tilt, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt)
                                            l11 = Ly
                                            l22 = Lx
                                        else if(gencluster == 1) then 
                                            call representative(newst, sites, Lx, Ly, xtransl, ytransl, rep, l1, l2, signt) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.            
                                            l11 = Lx
                                            l22 = Ly
                                        end if  
                                    call findstate(dim, rep, basis, loc) 
                                    if(loc <= 0) cycle   
                                    parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                                    parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                                    parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                                    parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                                    sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                                    psiprime(loc) = psiprime(loc) + sign * signt * (-(t)**2)*sqrt(norm(loc)/norm(j)) * psi(j, l) * exp(-ii*2*pi*((dble(k1*l1)/dble(l11)) + (dble(k2*l2)/dble(l22))))
                                end if     
                            end if 
                        end do !Loop over basis states    
                        !a1 translations 
                        cntr = 0
                        do q = 1, nbonds
                            if(bsites(1,q) == site1 .and. bsites(5,q) == 1) then 
                                site1 = bsites(2,q)
                                cntr = cntr + 1
                            else if(bsites(1,q) == site2 .and. bsites(5,q) == 1) then
                                site2 = bsites(2,q)
                                cntr = cntr + 1
                            else if(bsites(1,q) == refsite1 .and. bsites(5,q) == 1) then
                                refsite1 = bsites(2,q)
                                cntr = cntr + 1
                            else if(bsites(1,q) == refsite2 .and. bsites(5,q) == 1) then
                                refsite2 = bsites(2,q)
                                cntr = cntr + 1
                            else if(cntr == 4) then 
                                exit 
                            end if
                        end do 

                    end do 

                    !a2 translations 
                    cntr = 0
                    do q = 1, nbonds
                        if(bsites(1,q) == site1 .and. bsites(5,q) == 2) then 
                            site1 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == site2 .and. bsites(5,q) == 2) then
                            site2 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == refsite1 .and. bsites(5,q) == 2) then
                            refsite1 = bsites(2,q)
                            cntr = cntr + 1
                        else if(bsites(1,q) == refsite2 .and. bsites(5,q) == 2) then
                            refsite2 = bsites(2,q)
                            cntr = cntr + 1
                        else if(cntr == 4) then 
                            exit 
                        end if
                    end do 
                end do          
                bondcurrent(i) = bondcurrent(i) + dist * dot_product(psi(1:dim,k), psiprime)
            end do !Loop over degenerate ground states
        end do !Loop over degenerate ground states
        current = current + dble(phase) * bondcurrent(i) / (Lx*Ly)
    end do !Loop over bonds 

    deallocate(psiprime)

    if(sl == "A") then 
        filename = trim_name(dir//"Acurrent_coh_sum_"//parameters)
    else if(sl == "B") then 
        filename = trim_name(dir//"Bcurrent_coh_sum_"//parameters)
    end if 
    filename = trim_name(filename)
    print* ,filename, 'name'
    open(unit,file=filename)
    write(unit,*) current 
    close(unit)
    
    return

end subroutine current_cs_c2

subroutine currmat(dir, unit, filename, sl, ndeg, t, dim, sites, basis, bsites, xy, nbonds, refbond, psi, nnnVec)

    implicit none
    integer, intent(in) :: unit, ndeg, nbonds, refbond, sites
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    integer, intent(in) :: xy(4,nbonds)
    double precision, intent(in) :: t 
    double precision, intent(in) :: psi(dim, ndeg), nnnVec(2, 3)
    character(len=*), intent(in) :: dir, sl, filename
    
    integer(kind=8) :: loc = 0, newst = 0
    integer :: i = 0, j = 0, k = 0, l = 0
    integer :: sign = 0
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0, current = 0.d0 
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:), gs(:)
    double precision, allocatable :: currents(:)
    double precision, allocatable :: currents2(:)
    double precision, allocatable :: evals(:)
    double precision, allocatable :: bondcurrents(:,:), bondcurrent(:)
    double precision, allocatable :: currentmat(:,:)
    double precision, allocatable :: bondcurrentmat(:,:,:)
    character*512 :: name 

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)

    !Definition of lattice vectors: (5th row of bsites)
    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(psiprime)) deallocate(psiprime)
    if(allocated(bondcurrentmat)) deallocate(bondcurrentmat)
    if(allocated(currentmat)) deallocate(currentmat)
    if(allocated(currents2)) deallocate(currents2)
    if(allocated(bondcurrents)) deallocate(bondcurrents)
    allocate(psiprime(dim))
    allocate(bondcurrentmat(nbonds, ndeg, ndeg))
    allocate(currentmat(ndeg, ndeg))
    allocate(currents2(ndeg))
    allocate(bondcurrents(nbonds, ndeg))

    psiprime = 0.d0
    refvec   = vecs(1:2,bsites(5,refbond))
    refsite1 = bsites(1,refbond)
    refsite2 = bsites(2,refbond)      

    ! do i = 1, ndeg
        !     do j = 1, ndeg
        !         if ( abs(dot_product(psi(1:dim,i),psi(1:dim,j))) > 10.**(-10) .and. abs(dot_product(psi(1:dim,i),psi(1:dim,j))) < 0.999999999d0 ) then 
        !             print*,'Currmat: Groundstates ',i, 'and' ,j,' are not orthorgonal.', dot_product(psi(1:dim,i),psi(1:dim,j))               
        !             error stop 
        !         end if 
        !     end do 
        ! end do 
    
    currentmat = 0.d0     
    bondcurrentmat = 0.d0  
    do l = 1, ndeg 
        do k = 1, ndeg 
            do i = 1, nbonds                            
                site1 = bsites(1,i)
                site2 = bsites(2,i)
                if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
                    .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
                phase = bsites(3,i)   
                vec   = vecs(1:2,bsites(5,i))
                dist  = dot_product(vec, refvec)
                psiprime = 0.d0
                do j = 1, dim
                    if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                        newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                        if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                            ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                            newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                            call findstate(dim, newst, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j, l)
                        else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                            ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                            newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                            call findstate(dim, newst, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + (-1) * sign * (-(t)**2) * psi(j, l)
                        end if     
                    else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                        newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                        if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                            ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                            newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                            call findstate(dim, newst, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + (-1) * sign * (-(t)**2) * psi(j, l)
                        else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                            ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                            newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                            call findstate(dim, newst, basis, loc) 
                            if(loc <= 0) cycle   
                            parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j, l)
                        end if     
                    end if 
                end do !Loop over basis states                  
                bondcurrentmat(i,k,l) = bondcurrentmat(i,k,l) + dble(phase) * dist * dot_product(psi(1:dim,k), psiprime)
                ! bondcurrentmat(i,k,l) = bondcurrentmat(i,k,l) + dist * dot_product(psi(1:dim,k), psiprime)
                currentmat(k,l) = currentmat(k,l) + bondcurrentmat(i,k,l)
                ! currentmat(k,l) = currentmat(k,l) + dble(phase) * bondcurrentmat(i,k,l)
            end do !Loop over bonds 
        end do !Loop over degenerate ground states
    end do !Loop over degenerate ground states
    deallocate(psiprime)

    currentmat = currentmat / ndeg 


    call exactdiag(.True., ndeg, currentmat, currents)

    print*,' Current eigenvalues'
    do i = 1, ndeg 
        print*, currents(i)
    end do
    print*,'|===================|'
    print*,' Total'
    print*,sum(currents) 
    ! currents = 0.d0 
    print*,'|===================|'
    if(allocated(gs)) deallocate(gs)
    allocate(gs(dim))

    gs = 0.d0 
    do i = 1, ndeg 
        gs = gs + currents(i) * psi(1:dim, i)
        print* ,norm2(psi(1:dim,i)), 'norm2(psi(1:dim,i))'
    end do
    print* ,norm2(gs), 'norm2(gs)'
    gs = gs / sqrt(dble(ndeg))
    print* ,norm2(gs), 'norm2(gs)'

    call slcurrent(1.0d0, dim, sites, basis, bsites, xy, nbonds, refbond, gs, nnnVec, bondcurrent, current)
    
    print*,'|===================|'
    print*,' Current, manifold'
    print*, current
    print*,'|===================|'
    ! do j = 1, ndeg 
    !     gs = 0.d0 
    !     do i = 1, ndeg 
    !         gs = gs + currentmat(i,j) * psi(1:dim, i)
    !     end do
    !     call slcurrent(1.0d0, dim, sites, basis, bsites, xy, nbonds, refbond, gs, nnnVec, bondcurrent, current)
    !     bondcurrents(1:nbonds, j) = bondcurrent 
    !     currents(j) = current 
    !     print*,'|===================|'
    !     print*,' Current, eigenstate'
    !     print*, current, j
    !     print*,'|===================|'
    !     current = 0.d0 
    ! end do 
    
    if(sl == "A") then 
        name = trim_name(dir//"Acurrent_matrix_"//filename)
    else if(sl == "B") then 
        name = trim_name(dir//"Bcurrent_matrix_"//filename)
    end if 
    name = trim_name(name)
    
    open(unit,file=name)
    do i = 1, ndeg 
        write(unit,*) currents(i)
    end do 
    close(unit)
    if(sl == "A") then 
        name = dir//"Abondcurrent_matrix_"//filename
    else if(sl == "B") then 
        name = dir//"Bbondcurrent_matrix_"//filename
    end if 
    name = trim_name(name)
    open(unit,file=name)
    do i = 1, ndeg 
        write(unit,*) bondcurrents(1:nbonds,i)
    end do 
    close(unit)

    ! bondcurrents = 0.d0 
        ! do i = 1, nbonds
        !     call exactdiag(.True., ndeg, bondcurrentmat(i, 1:ndeg, 1:ndeg), evals)
        !     bondcurrents(i, 1:ndeg) = evals
        !     if(allocated(evals)) deallocate(evals)
        ! end do 
            
        ! if(sl == "A") then 
        !     name = dir//"Acurrent_matrix2_"//filename
        ! else if(sl == "B") then 
        !     name = dir//"Bcurrent_matrix2_"//filename
        ! end if 
        ! name = trim_name(name)
        ! open(unit,file=name)
        ! do i = 1, ndeg 
        !     currents2(i) = sum(bondcurrents(1:nbonds, i))
        !     write(unit,*) currents2(i)
        !     print*,'Current of psi ',i, currents2(i)
        ! end do 
        ! close(unit)

    return

end subroutine currmat

subroutine current_ics(ndeg, t, dim, sites, basis, bsites, nbonds, refbond, psi, nnnVec, bondcurrent, current)

    implicit none
    integer, intent(in) :: ndeg, nbonds, refbond, sites
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    double precision, intent(in) :: t 
    double precision, intent(in) :: psi(dim, ndeg), nnnVec(2, 3)
    double precision, intent(out) :: current
    double precision, allocatable, intent(out) :: bondcurrent(:)
    
    integer(kind=8) :: loc = 0, newst = 0
    integer :: i = 0, j = 0, l = 0
    integer :: sign = 0
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:)

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)
    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bondcurrent)) deallocate(bondcurrent)
    if(allocated(psiprime)) deallocate(psiprime)
    allocate(bondcurrent(nbonds))
    allocate(psiprime(dim))

    psiprime = 0.d0
    refvec   = vecs(1:2,bsites(5,refbond))
    refsite1 = bsites(1,refbond)
    refsite2 = bsites(2,refbond)      

    current = 0     
    do i = 1, nbonds            
        bondcurrent(i) = 0 
        site1 = bsites(1,i)
        site2 = bsites(2,i)
        if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
            .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
        phase = bsites(3,i)   
        vec   = vecs(1:2,bsites(5,i))
        dist = dot_product(vec, refvec)

        do l = 1, ndeg  
            psiprime = 0 
            do j = 1, dim
                if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                    newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                    if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                        ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                        newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                        call findstate(dim, newst, basis, loc) 
                        if(loc <= 0) cycle   
                        
                        parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                        parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                        parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                        parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                        sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                        psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j, l)
                    else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                        ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                        newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                        call findstate(dim, newst, basis, loc) 
                        if(loc <= 0) cycle   
                        
                        parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                        parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                        parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                        parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                        sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                        psiprime(loc) = psiprime(loc) + (-1) * sign * (-(t)**2) * psi(j, l)
                    end if     
                else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                    newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                    if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                        ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                        newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                        call findstate(dim, newst, basis, loc) 
                        if(loc <= 0) cycle   
                        
                        parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                        parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                        parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                        parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                        sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                        psiprime(loc) = psiprime(loc) + (-1) * sign * (-(t)**2) * psi(j, l)
                    else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                        ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                        newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                        call findstate(dim, newst, basis, loc) 
                        if(loc <= 0) cycle   
                        
                        parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                        parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                        parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                        parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                        sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                        psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j, l)
                    end if     
                end if 
            end do !Loop over basis states                  
            bondcurrent(i) = bondcurrent(i) + dist * dot_product(psi(1:dim,l), psiprime)
        
        end do !Loop over degenerate ground states
        current = current + dble(phase) * bondcurrent(i)
    end do !Loop over bonds 

    deallocate(psiprime)
    
    return

end subroutine current_ics

subroutine current_cs(dir, unit, parameters, sl, ndeg, t, dim, sites, basis, bsites, nbonds, refbond, psi, nnnVec)

    implicit none
    integer, intent(in) :: unit, ndeg, nbonds, refbond, sites
    integer(kind=8), intent(in) :: dim
    integer(kind=8), intent(in) :: basis(dim)
    integer, intent(in) :: bsites(5,nbonds)
    double precision, intent(in) :: t 
    double precision, intent(in) :: psi(dim, ndeg), nnnVec(2, 3)
    character(len=*), intent(in) :: dir, sl, parameters
    
    integer(kind=8) :: loc = 0, newst = 0
    integer :: i = 0, j = 0, k = 0, l = 0
    integer :: sign = 0
    integer :: site1 = 0, site2 = 0, phase = 0   
    integer :: refsite1 = 0, refsite2 = 0
    integer :: parity1 = 0, parity2 = 0, parity3 = 0, parity4 = 0
    double precision :: dist = 0.d0
    double precision, dimension(2)   :: vec, refvec
    double precision, dimension(2,3) :: vecs
    double precision, allocatable :: psiprime(:)
    double precision :: current
    double precision, allocatable :: bondcurrent(:)
    character*400 :: name

    !!$omp threadprivate(loc, newst, i, j, forward, backwards, signf, signb, current_temp, bondcurrent, psiprime, psiprimeij, psiprimeprime)
    vecs(1:2,1) = nnnVec(1:2,1)/norm2(nnnVec(1:2,1))    
    vecs(1:2,2) = nnnVec(1:2,2)/norm2(nnnVec(1:2,2))
    vecs(1:2,3) = nnnVec(1:2,3)/norm2(nnnVec(1:2,3))

    if(allocated(bondcurrent)) deallocate(bondcurrent)
    if(allocated(psiprime)) deallocate(psiprime)
    allocate(bondcurrent(nbonds))
    allocate(psiprime(dim))

    psiprime    = 0.d0
    refvec      = vecs(1:2,bsites(5,refbond))
    refsite1    = bsites(1,refbond)
    refsite2    = bsites(2,refbond)      

    do i = 1, ndeg
        do j = 1, ndeg
            if ( abs(dot_product(psi(1:dim,i),psi(1:dim,j))) > 10.**(-10) .and. abs(dot_product(psi(1:dim,i),psi(1:dim,j))) < 0.999999999d0 ) then 
                print*,'current_cs coherent sum: Groundstates ',i, 'and' ,j,' are not orthorgonal.', dot_product(psi(1:dim,i),psi(1:dim,j))               
                error stop 
            end if 
        end do 
    end do 
    current = 0     
    do i = 1, nbonds            
        bondcurrent(i) = 0 
        site1 = bsites(1,i)
        site2 = bsites(2,i)
        if (site1 == bsites(1,refbond) .or. site2 == bsites(1,refbond) &
            .or. site1 == bsites(2,refbond) .or. site2 == bsites(2,refbond)) cycle !Exclude bonds which share a site with the reference bond
        phase = bsites(3,i)   
        vec   = vecs(1:2,bsites(5,i))
        dist = dot_product(vec, refvec)

        do l = 1, ndeg 
            do k = 1, ndeg 
                psiprime = 0 
                do j = 1, dim
                    if(btest(basis(j), refsite1 - 1 ) .and. .not.(btest(basis(j), refsite2 - 1 ) ) ) then !Forward hopping from refsites 1->2
                        newst = ibclr( ibset( basis(j), refsite2 - 1 ), refsite1 - 1 ) !Create on refsite 2, annihilate on refsite 1   
                        if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Forward forward hopping from sites 1->2       
                            ! cd_site2 * c_site1 * cd_refs2 * c_refs1
                            newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                            call findstate(dim, newst, basis, loc) 
                            if(loc <= 0) cycle   
                            if( basis(loc) .ne. newst) error stop "Subroutine current_cs coherent sum: Wrong location of scattered basis state."
                            parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1, sites ) ) !Parity of c_site1
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j, l)
                        else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Forward backwards hopping from sites 2->1                   
                            ! cd_site1 * c_site2 * cd_refs2 * c_refs1
                            newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                            call findstate(dim, newst, basis, loc) 
                            if(loc <= 0) cycle   
                            if( basis(loc) .ne. newst) error stop "Subroutine current_cs coherent sum: Wrong location of scattered basis state."
                            parity1 = popcnt( ibits( basis(j), refsite1, sites ) ) !Parity of c_refs1
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite1 - 1), refsite2, sites ) ) !Parity of cd_refs2 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + (-1) * sign * (-(t)**2) * psi(j, l)
                        end if     
                    else if( btest( basis(j), refsite2 - 1 ) .and. .not.( btest( basis(j), refsite1 - 1 ) ) ) then !Backwards hopping from refsites 2->1 
                        newst = ibclr( ibset( basis(j), refsite1 - 1 ), refsite2 - 1 ) !Create on refsite 1, annihilate on refsite 2
                        if(btest(basis(j), site1 - 1 ) .and. .not.(btest(basis(j), site2 - 1 ) ) ) then !Backwards forward hopping from sites 1->2       
                            ! cd_site2 * c_site1 * cd_refs1 * c_refs2
                            newst = ibclr( ibset( newst, site2 - 1 ), site1 - 1 ) !Create on site 2, annihilate on site 1 
                            call findstate(dim, newst, basis, loc) 
                            if(loc <= 0) cycle   
                            if( basis(loc) .ne. newst) error stop "Subroutine current_cs coherent sum: Wrong location of scattered basis state."
                            parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1, sites ) ) !Parity of c_site1
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite2 - 1), refsite1 - 1), site1 - 1), site2, sites ) ) !Parity of cd_site2
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + (-1) * sign * (-(t)**2) * psi(j, l)
                        else if(btest(basis(j), site2 - 1 ) .and. .not.(btest(basis(j), site1 - 1 ) ) ) then !Backwards backwards hopping from sites 2->1                   
                            ! cd_site1 * c_site2 * cd_refs1 * c_refs2
                            newst = ibclr( ibset( newst, site1 - 1 ), site2 - 1 ) !Create on site 1, annihilate on site 2
                            call findstate(dim, newst, basis, loc) 
                            if(loc <= 0) cycle   
                            if( basis(loc) .ne. newst) error stop "Subroutine current_cs coherent sum: Wrong location of scattered basis state."
                            parity1 = popcnt( ibits( basis(j), refsite2, sites ) ) !Parity of c_refs2
                            parity2 = popcnt( ibits( ibclr(basis(j), refsite2 - 1), refsite1, sites ) ) !Parity of cd_refs1 
                            parity3 = popcnt( ibits( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2, sites ) ) !Parity of c_site2
                            parity4 = popcnt( ibits( ibclr( ibset( ibclr(basis(j), refsite1 - 1), refsite2 - 1), site2 - 1), site1, sites ) ) !Parity of cd_site1
                            sign = (-1)**(parity1 + parity2 + parity3 + parity4) 
                            psiprime(loc) = psiprime(loc) + sign * (-(t)**2) * psi(j, l)
                        end if     
                    end if 
                end do !Loop over basis states                  
                bondcurrent(i) = bondcurrent(i) + dist * dot_product(psi(1:dim,k), psiprime)
            end do !Loop over degenerate ground states
        end do !Loop over degenerate ground states
        current = current + dble(phase) * bondcurrent(i)
    end do !Loop over bonds 
    print* ,current, 'coherent sum '
    

    if(sl == "A") then 
        name = trim_name(dir//"Acurrent_coh_sum_"//parameters)
    else if(sl == "B") then 
        name = trim_name(dir//"Bcurrent_coh_sum_"//parameters)
    end if 
    open(unit,file=name)
    write(unit,*) current
    close(unit)

    if(sl == "A") then 
        name = trim_name(dir//"Abondcurrent_coh_sum_"//parameters)
    else if(sl == "B") then 
        name = trim_name(dir//"Bbondcurrent_coh_sum_"//parameters)
    end if 
    name = trim_name(name)
    open(unit,file=name)
    write(unit,*) bondcurrent
    close(unit)

    deallocate(psiprime)
    
    return

end subroutine current_cs







end module routines
