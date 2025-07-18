module test_utils 
    implicit none 

    ! Module for testing the functionality of the code
    contains 

    subroutine check_spectrum_dp(dim, nnz, nev, nest, energies, eigstate, ndeg, rc, mat)
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
    
    end subroutine check_spectrum_dp

subroutine check_spectrum_dc(dim, feast, nev, nest, energies, eigstate)
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

end subroutine check_spectrum_dc

end module test_utils 
