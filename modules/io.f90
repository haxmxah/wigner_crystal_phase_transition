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
            read(unit = input_file_unit, fmt = *) N
            read(unit = input_file_unit, fmt = *) freeze_mc_steps_scale
            read(unit = input_file_unit, fmt = *) alpha
            read(unit = input_file_unit, fmt = *) density
            read(unit = input_file_unit, fmt = *) q
            read(unit = input_file_unit, fmt = *) T_i
            read(unit = input_file_unit, fmt = *) T_f
            read(unit = input_file_unit, fmt = *) T_step

            close(unit = input_file_unit)
        end subroutine read_input
end module io
