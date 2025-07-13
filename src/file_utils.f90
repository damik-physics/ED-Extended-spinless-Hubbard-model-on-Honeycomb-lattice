module file_utils
    implicit none
    character(len=256) :: outdir

contains

    subroutine setup_output_directory()
        implicit none
        character(len=64) :: timestamp
        integer :: values(8)
        character(len=512) :: cmd

        ! Create a timestamp for the output directory
        ! This will create a directory named 'run_YYYYMMDD_HHMMSS'
        ! where YYYYMMDD is the date and HHMMSS is the time
        call date_and_time(values = values)
        write(timestamp, '(I4.4,I2.2,I2.2,"_",I2.2,I2.2,I2.2)') &
             values(1), values(2), values(3), values(5), values(6), values(7)

        outdir = 'output/run_' // trim(timestamp) // '/'
        cmd = 'mkdir -p ' // trim(outdir)
        call system(cmd)

        ! Also copy the input file for record-keeping
        cmd = 'cp input.nml ' // trim(outdir) // 'input.nml'
        call system(cmd)
        
    end subroutine setup_output_directory

    subroutine create_output_subdirs(out_dir)
        implicit none
        character(len=*), intent(in) :: out_dir
        character(len=512) :: cmd
        integer :: exitstat
        ! Create subdirectories for different output types
        cmd = 'mkdir -p ' // trim(out_dir) // '/correlations ' // &
                        trim(out_dir) // '/hamiltonians ' // &
                        trim(out_dir) // '/lattice_data ' // &
                        trim(out_dir) // '/logs ' // &
                        trim(out_dir) // '/parameters ' // &
                        trim(out_dir) // '/plots ' // &
                        trim(out_dir) // '/spectra ' // &
                        trim(out_dir) // '/states'

        call execute_command_line(trim(cmd), exitstat=exitstat)
        if (exitstat /= 0) then
            print *, 'Error creating subdirectories. Exit status:', exitstat
            stop 1
        end if
    end subroutine create_output_subdirs


    subroutine write_parameters_json(ucx, ucy, tilted, cluster, bc, V, V2, filling, irrep)
        implicit none
        integer, intent(in) :: ucx, ucy, tilted
        character(*), intent(in) :: cluster, bc, irrep
        real(8), intent(in) :: V, V2, filling
        integer :: unit

        open(newunit=unit, file=trim(outdir)//'parameters.json', status='replace')
        write(unit,'(A)') '{'
        write(unit,'(A,I0,A)') '"ucx": ', ucx, ','
        write(unit,'(A,I0,A)') '"ucy": ', ucy, ','
        write(unit,'(A,I0,A)') '"tilted": ', tilted, ','
        write(unit,'(A,"""",A,"""",A)') '"cluster": "', trim(cluster), '",'
        write(unit,'(A,"""",A,"""",A)') '"bc": "', trim(bc), '",'
        write(unit,'(A,F6.3,A)') '"V": ', V, ','
        write(unit,'(A,F6.3,A)') '"V2": ', V2, ','
        write(unit,'(A,F6.3,A)') '"filling": ', filling, ','
        write(unit,'(A,"""",A,"""")') '"irrep": "', trim(irrep), '"'
        write(unit,'(A)') '}'
        close(unit)
    end subroutine write_parameters_json

    
end module file_utils
