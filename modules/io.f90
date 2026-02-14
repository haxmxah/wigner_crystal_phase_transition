module io
    use config

    implicit none

    contains
        ! Reads external input files that follow the structure of
        ! input_parameters.in
        subroutine read_input(input_file)
            implicit none

            character(*), intent(in) :: input_file

            integer :: ios, input_file_unit

            open(newunit = input_file_unit, file = input_file, &
                 status = 'old', action = 'read', iostat = ios)

            if (ios /= 0) then
                print *, 'Error reading input file ', input_file
                stop
            end if

            !
            ! Read input file contents,
            !
            ! Make sure these are in the same order that the ones defined in
            ! the input file.
            ! In case you wish to add an extra variable to be read from the
            ! input file, it should also be added here below.
            !
            read(unit = input_file_unit, fmt = *) melting
            read(unit = input_file_unit, fmt = *) input_positions_file
            read(unit = input_file_unit, fmt = *) N
            read(unit = input_file_unit, fmt = *) freeze_mc_steps_scale
            read(unit = input_file_unit, fmt = *) alpha
            read(unit = input_file_unit, fmt = *) C0
            read(unit = input_file_unit, fmt = *) density
            read(unit = input_file_unit, fmt = *) q
            read(unit = input_file_unit, fmt = *) T_i
            read(unit = input_file_unit, fmt = *) T_f
            read(unit = input_file_unit, fmt = *) T_step

            close(unit = input_file_unit)
        end subroutine read_input

    ! Reads the positions from a file .xyz 
    subroutine read_positions_xyz(filename, positions, N_check)
        character(*), intent(in) :: filename
        real(dp), intent(out) :: positions(:, :)
        integer, intent(in) :: N_check 
        
        integer :: unit_xyz, ios, i, N_file
        character(len=100) :: dummy_line
        character(len=2) :: dummy_char

        open(newunit=unit_xyz, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Cannot open the file ", filename
            stop
        end if

        ! Read header
        read(unit_xyz, *) N_file
        read(unit_xyz, '(A)') dummy_line 

        if (N_file /= N_check) then
            print *, "Error: N in file (", N_file, ") does not match the input (", N_check, ")"
            stop
        end if

        ! Read positions
        do i = 1, N_file
            read(unit_xyz, *) dummy_char, positions(1, i), positions(2, i), positions(3, i)
        end do

        close(unit_xyz)
        print *, "Positions loaded: ", filename
    end subroutine read_positions_xyz
end module io
