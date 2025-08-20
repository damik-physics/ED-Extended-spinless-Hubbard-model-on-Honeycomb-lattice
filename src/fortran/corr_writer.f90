module corr_writer
    implicit none
    private

    interface append_corr_csv
        module procedure append_corr_csv_dp
        module procedure append_corr_csv_dc
    end interface append_corr_csv
    
    public :: append_corr_csv_dp, append_corr_csv_dc, append_corr_csv


    
contains

    !=========================
    ! PUBLIC WRAPPERS
    !=========================

    subroutine append_corr_csv_dp(fname, dir, ti, ndis, k1, k2, conf, v1, v2, corr)
        character(len=*), intent(in) :: fname, dir
        integer, intent(in), optional :: ti, ndis, k1, k2, conf
        real(kind=8), intent(in) :: v1, v2
        real(kind=8), intent(in) :: corr(..)   ! assumed-rank

        call append_corr_csv_core(fname, dir, ti, ndis, k1, k2, conf, v1, v2, corr, '(",",F12.6)')
    end subroutine append_corr_csv_dp


    subroutine append_corr_csv_dc(fname, dir, ti, ndis, k1, k2, conf, v1, v2, corr)
        character(len=*), intent(in) :: fname, dir
        integer, intent(in), optional :: ti, ndis, k1, k2, conf
        real(kind=8), intent(in) :: v1, v2
        complex(kind=8), intent(in) :: corr(..)   ! assumed-rank

        call append_corr_csv_core(fname, dir, ti, ndis, k1, k2, conf, v1, v2, corr, '(",",(F12.6,F12.6))')
    end subroutine append_corr_csv_dc


    !=========================
    ! PRIVATE CORE
    !=========================

    subroutine append_corr_csv_core(fname, dir, ti, ndis, k1, k2, conf, v1, v2, corr, fmt)
        character(len=*), intent(in) :: fname, dir
        integer, intent(in), optional :: ti, ndis, k1, k2, conf
        real(kind=8), intent(in) :: v1, v2
        class(*), intent(in) :: corr(..)      ! works for both real and complex
        character(len=*), intent(in) :: fmt   ! format descriptor

        character(len=512) :: fpath, errmsg
        integer :: i, ios
        logical :: file_exists

        ! Construct full path
        fpath = trim(dir)//'correlations/'//fname//'.csv'
        errmsg = 'Error opening '//fname//'.csv'

        ! Check existence
        inquire(file=fpath, exist=file_exists)

        ! Create header if needed
        if (.not. file_exists) then
            open(unit=40, file=fpath, status='new', action='write', iostat=ios)
            if (ios /= 0) stop errmsg
            write(40,'(A)', advance='no') 'v1,v2'
            if (present(conf)) write(40,'(A)', advance='no') ',conf'

            select rank(corr)
            rank (0)  ! scalar
                write(40,'(A)', advance='no') ',c'
            rank (1)  ! 1D array
                do i=1, size(corr)
                    write(40,'(A,I0)', advance='no') ',c', i
                end do
            rank default
                stop 'Unsupported rank for corr'
            end select

            write(40,*)
            close(40)
        end if

        ! Append data
        open(unit=41, file=fpath, status='old', action='write', position='append', iostat=ios)
        if (ios /= 0) stop errmsg
        write(41,'(F8.4,",",F8.4)', advance='no') v1, v2
        if (present(conf)) write(41,'(",",I0)', advance='no') conf

        select rank(corr)
        rank (0)  ! scalar
            select type(corr)
            type is (real(kind=8))
                write(41, fmt, advance='no') corr
            type is (complex(kind=8))
                write(41, fmt, advance='no') corr
            class default
                stop 'Unsupported type for corr'
            end select
        rank (1)  ! 1D array
            do i=1, size(corr)
                select type(corr)
                type is (real(kind=8))
                    write(41, fmt, advance='no') corr(i)
                type is (complex(kind=8))
                    write(41, fmt, advance='no') corr(i)
                class default
                    stop 'Unsupported type for corr'
                end select
            end do
        rank default
            stop 'Unsupported rank for corr'
        end select

        write(41,*)
        close(41)
    end subroutine append_corr_csv_core

end module corr_writer
