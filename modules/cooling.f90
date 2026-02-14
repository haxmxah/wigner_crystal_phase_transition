module cooling
    use, intrinsic :: iso_fortran_env, only: error_unit
    use config

    contains

    subroutine get_temperatures_list(temperatures, cooling_method)
        implicit none

        integer :: idx, num_temperatures
        real(kind = dp), intent(inout), allocatable :: temperatures(:)
        character(len=256), intent(in) :: cooling_method

        select case (trim(cooling_method))
            case ('linear')
                
                num_temperatures = nint(abs((T_i - T_f))/T_step) + 1
                allocate(temperatures(num_temperatures))

                T_i = T_i + T_step
                do idx = 1, num_temperatures
                    temperatures(idx) = T_i + T_step*idx
                end do

            case ('annealing')
                num_temperatures = nint(log(T_f/T_i) &
                                        /log(T_step)) &
                                 + 1

                allocate(temperatures(num_temperatures))

                temperatures(1) = T_i
                do idx = 2, num_temperatures
                    temperatures(idx) = temperatures(idx - 1)*T_step
                end do

            case default
                write (unit = error_unit, fmt = '(2a/, a)') &
                    'Invalid value for "cooling_method": ', trim(cooling_method), &
                    'The cooling method must either be "linear" or "annealing".'
                stop
        end select
    end subroutine get_temperatures_list
end module cooling
