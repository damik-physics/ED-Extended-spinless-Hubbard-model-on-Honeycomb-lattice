module ham_hdf5_io
    use hdf5
    implicit none
    private

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: i8 = selected_int_kind(18)

    public :: save_ham_hdf5_coo, load_hdf5_ham_coo

contains

    subroutine save_ham_hdf5_coo(h5file, group, nrows, ncols, rc, values, overwrite)
        !! Save COO (rc, values) under /ham/<group>/... . Handles real(dp) and complex(dp).
        character(len=*), intent(in)              :: h5file, group
        integer(i8),      intent(in)              :: nrows, ncols
        integer(i8),      intent(in)              :: rc(:, :)      ! (nnz,2) 1-based
        class(*),         intent(in)              :: values(:)     ! real(dp) or complex(dp)
        logical,          intent(in), optional    :: overwrite

        integer(HID_T) :: fid, gid, dspace, dset
        integer(HID_T) :: dtype 
        integer :: ierr
        integer(HSIZE_T) :: dims1(1), dims2(2)
        integer(i8), allocatable :: shape(:)
        real(dp), allocatable :: vr(:), vi(:)
        logical :: ow
        character(len=:), allocatable :: base

        ow = .false.; if (present(overwrite)) ow = overwrite
        base = "/ham/"//trim(group)

        ! Open or create file
        if (ow) then
            call h5fcreate_f(trim(h5file), H5F_ACC_TRUNC_F, fid, ierr)
        else
            call h5fopen_f(trim(h5file), H5F_ACC_RDWR_F, fid, ierr)
            if (ierr /= 0) call h5fcreate_f(trim(h5file), H5F_ACC_TRUNC_F, fid, ierr)
        end if
        if (ierr /= 0) stop "HDF5: cannot open/create file"

        ! Create /ham and /ham/<group> (idempotent: create then close; if exists, open)
        call ensure_group(fid, "/ham", gid, ierr)
        call h5gclose_f(gid, ierr)
        call ensure_group(fid, base, gid, ierr)

        ! Write shape = [nrows, ncols]
        allocate(shape(2)); shape = [nrows, ncols]
        dims1 = 2
        call h5screate_simple_f(1, dims1, dspace, ierr)
        dtype = H5T_NATIVE_INTEGER
        call create_or_replace_dataset(fid, trim(base)//"/shape", dtype, dspace, dset, ierr)
        call h5dwrite_f(dset, dtype, shape, dims1, ierr)
        call h5dclose_f(dset, ierr); call h5sclose_f(dspace, ierr)
        deallocate(shape)

        ! Write rc (nnz,2)
        dims2 = [int(size(rc,1),HSIZE_T), 2_HSIZE_T]
        call h5screate_simple_f(2, dims2, dspace, ierr)
        call create_or_replace_dataset(fid, trim(base)//"/rc", H5T_NATIVE_INTEGER, dspace, dset, ierr)
        call h5dwrite_f(dset, H5T_NATIVE_INTEGER, rc, dims2, ierr)
        call h5dclose_f(dset, ierr); call h5sclose_f(dspace, ierr)

        ! Write values (real or complex split)
        select type(values)
        type is (integer(i8))
            dims1 = int(size(values, kind=8), HSIZE_T)
            call h5screate_simple_f(1, dims1, dspace, ierr)
            call create_or_replace_dataset(fid, trim(base)//"/values", H5T_NATIVE_INTEGER, dspace, dset, ierr)
            call h5dwrite_f(dset, H5T_NATIVE_INTEGER, values, dims1, ierr)
            call h5dclose_f(dset, ierr)
            call h5sclose_f(dspace, ierr)

        type is (real(dp))
            dims1 = int(size(values),HSIZE_T)
            call h5screate_simple_f(1, dims1, dspace, ierr)
            call create_or_replace_dataset(fid, trim(base)//"/values", H5T_NATIVE_DOUBLE, dspace, dset, ierr)
            call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, values, dims1, ierr)
            call h5dclose_f(dset, ierr); call h5sclose_f(dspace, ierr)

        type is (complex(dp))
            allocate(vr(size(values)), vi(size(values)))
            vr = real(values, kind=dp); vi = aimag(values)

            dims1 = int(size(values),HSIZE_T)
            call h5screate_simple_f(1, dims1, dspace, ierr)
            call create_or_replace_dataset(fid, trim(base)//"/values_real", H5T_NATIVE_DOUBLE, dspace, dset, ierr)
            call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, vr, dims1, ierr)
            call h5dclose_f(dset, ierr); call h5sclose_f(dspace, ierr)

            call h5screate_simple_f(1, dims1, dspace, ierr)
            call create_or_replace_dataset(fid, trim(base)//"/values_imag", H5T_NATIVE_DOUBLE, dspace, dset, ierr)
            call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, vi, dims1, ierr)
            call h5dclose_f(dset, ierr); call h5sclose_f(dspace, ierr)

            deallocate(vr, vi)

        class default
            call h5gclose_f(gid, ierr); call h5fclose_f(fid, ierr)
            stop "save_ham_hdf5_coo: unsupported values type"
        end select

        call h5gclose_f(gid, ierr)
        call h5fclose_f(fid, ierr)
    end subroutine save_ham_hdf5_coo


    subroutine load_hdf5_ham_coo(h5file, group, nrows, ncols, rc, values_int, values_real, values_complex)
        !! Load COO from /ham/<group>/... . Exactly one of values_real or values_complex is allocated.
        character(len=*), intent(in)            :: h5file, group
        integer(i8),      intent(out)           :: nrows, ncols
        integer(i8),      allocatable, intent(out) :: rc(:, :)
        integer(i8), allocatable, intent(out), optional :: values_int(:)
        real(dp),         allocatable, intent(out), optional :: values_real(:)
        complex(dp),      allocatable, intent(out), optional :: values_complex(:)

        integer(HID_T) :: fid, dset, dspace
        integer :: ierr
        integer(HSIZE_T) :: dims1(1), dims2(2), maxdims1(1), maxdims2(2)
        integer(i8), allocatable :: shape(:)
        integer(HSIZE_T) :: nnz_h

        character(len=:), allocatable :: base
        logical :: has_re, has_im, has_val, has_val_int

        base = "/ham/"//trim(group)

        call h5fopen_f(trim(h5file), H5F_ACC_RDONLY_F, fid, ierr)
        if (ierr /= 0) stop "HDF5: cannot open file for reading"

        ! Read shape
        call h5dopen_f(fid, trim(base)//"/shape", dset, ierr); if (ierr /= 0) stop "shape missing"
        call h5dget_space_f(dset, dspace, ierr)
        call h5sget_simple_extent_dims_f(dspace, dims1, maxdims1, ierr)
        if (dims1(1) /= 2) stop "shape should have length 2"
        allocate(shape(2))
        call h5dread_f(dset, H5T_NATIVE_INTEGER, shape, dims1, ierr)
        call h5sclose_f(dspace, ierr); call h5dclose_f(dset, ierr)
        nrows = shape(1); ncols = shape(2)
        deallocate(shape)

        ! Read rc (nnz,2)
        call h5dopen_f(fid, trim(base)//"/rc", dset, ierr); if (ierr /= 0) stop "rc missing"
        call h5dget_space_f(dset, dspace, ierr)
        call h5sget_simple_extent_dims_f(dspace, dims2, maxdims2,ierr)
        if (dims2(2) /= 2) stop "rc second dim must be 2"
        nnz_h = dims2(1)
        allocate(rc(int(nnz_h,kind=i8),2))
        call h5dread_f(dset, H5T_NATIVE_INTEGER, rc, dims2, ierr)
        call h5sclose_f(dspace, ierr); call h5dclose_f(dset, ierr)

        ! Detect value layout
        has_val = dataset_exists(fid, trim(base)//"/values")
        has_re  = dataset_exists(fid, trim(base)//"/values_real")
        has_im  = dataset_exists(fid, trim(base)//"/values_imag")
        has_val_int = dataset_exists(fid, trim(base)//"/values_int")

        if (has_val .and. .not. has_re .and. .not. has_im) then
            if (.not. present(values_real)) stop "real values present but caller didn't provide values_real"
            allocate(values_real(int(nnz_h,kind=i8)))
            call read_1d_real(fid, trim(base)//"/values", values_real)

        else if (has_re .and. has_im .and. .not. has_val) then
            if (.not. present(values_complex)) stop "complex values present but caller didn't provide values_complex"
            call read_complex_pair(fid, trim(base)//"/values_real", trim(base)//"/values_imag", values_complex, nnz_h)
                
        else if (has_val_int .and. .not. has_val .and. .not. has_re .and. .not. has_im) then
            if (.not. present(values_int)) stop "integer(i8) values present but caller didn't provide values_int"
            allocate(values_int(int(nnz_h,kind=i8)))
            call read_1d_int8(fid, trim(base)//"/values_int", values_int)

        else
            stop "Invalid or mixed value datasets in HDF5"
        end if

        call h5fclose_f(fid, ierr)
    end subroutine load_hdf5_ham_coo


    !------------------ helpers ------------------

    subroutine ensure_group(fid, path, gid, ierr)
        integer(HID_T), intent(in)  :: fid
        character(len=*), intent(in):: path
        integer(HID_T), intent(out) :: gid
        integer, intent(out)        :: ierr
        ! Try open; if fail, create
        call h5gopen_f(fid, trim(path), gid, ierr)
        if (ierr /= 0) then
            call h5gcreate_f(fid, trim(path), gid, ierr)
        end if
    end subroutine ensure_group

    subroutine create_or_replace_dataset(fid, path, h5type, dspace, dset, ierr)
        integer(HID_T), intent(in)  :: fid, dspace
        integer(HID_T), intent(in)  :: h5type
        character(len=*), intent(in):: path
        integer(HID_T), intent(out) :: dset
        integer, intent(out)        :: ierr
        integer(HID_T) :: old
        call h5dopen_f(fid, trim(path), old, ierr)
        if (ierr == 0) then
            call h5dclose_f(old, ierr)
            call h5ldelete_f(fid, trim(path), ierr)
        end if
        call h5dcreate_f(fid, trim(path), h5type, dspace, dset, ierr)
    end subroutine create_or_replace_dataset

    logical function dataset_exists(loc_id, path)
        integer(HID_T), intent(in) :: loc_id
        character(len=*), intent(in):: path
        integer(HID_T) :: dset
        integer :: ierr
        call h5dopen_f(loc_id, trim(path), dset, ierr)
        if (ierr == 0) then
            call h5dclose_f(dset, ierr)
            dataset_exists = .true.
        else
            dataset_exists = .false.
        end if
    end function dataset_exists

    subroutine read_1d_int8(fid, name, arr)
        integer(HID_T), intent(in) :: fid
        character(len=*), intent(in) :: name
        integer(i8), allocatable, intent(out) :: arr(:)

        integer(HID_T) :: dset, dspace
        integer(HSIZE_T) :: dims(1), maxdims(1)
        integer :: ierr
        call h5dopen_f(fid, name, dset, ierr)
        call h5dget_space_f(dset, dspace, ierr)
        call h5sget_simple_extent_dims_f(dspace, dims, maxdims, ierr)
        allocate(arr(int(dims(1),kind=i8)))
        call h5dread_f(dset, H5T_NATIVE_INTEGER, arr, dims, ierr)
        call h5sclose_f(dspace, ierr)
        call h5dclose_f(dset, ierr)
    end subroutine read_1d_int8


    subroutine read_1d_real(fid, path, arr)
        integer(HID_T), intent(in) :: fid
        character(len=*), intent(in) :: path
        real(dp), allocatable, intent(inout) :: arr(:)
        integer(HID_T) :: dset, dspace
        integer :: ierr
        integer(HSIZE_T) :: dims(1), maxdims(1)
        call h5dopen_f(fid, trim(path), dset, ierr)
        call h5dget_space_f(dset, dspace, ierr)
        call h5sget_simple_extent_dims_f(dspace, dims, maxdims, ierr)
        if (.not. allocated(arr)) allocate(arr(int(dims(1),kind=i8)))
        call h5dread_f(dset, H5T_NATIVE_DOUBLE, arr, dims, ierr)
        call h5sclose_f(dspace, ierr); call h5dclose_f(dset, ierr)
    end subroutine read_1d_real

    subroutine read_complex_pair(fid, path_re, path_im, z, nnz_h)
        integer(HID_T), intent(in) :: fid
        character(len=*), intent(in) :: path_re, path_im
        complex(dp), allocatable, intent(out) :: z(:)
        integer(HSIZE_T), intent(in) :: nnz_h
        real(dp), allocatable :: vr(:), vi(:)
        allocate(z(int(nnz_h,kind=i8)), vr(int(nnz_h,kind=i8)), vi(int(nnz_h,kind=i8)))
        call read_1d_real(fid, path_re, vr)
        call read_1d_real(fid, path_im, vi)
        z = cmplx(vr, vi, kind=dp)
        deallocate(vr, vi)
    end subroutine read_complex_pair

end module ham_hdf5_io
