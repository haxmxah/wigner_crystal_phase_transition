module config
    use iso_fortran_env, only: dp => real64
    implicit none

    integer :: N, melting
    character(len=256) :: input_positions_file
    real(dp) :: T_i, T_f, T_step, alpha, C0, density, q, freeze_mc_steps_scale
    
    ! Mathematical constants
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

contains

    ! Uniform Random generator [min, max]
    function random_uniform(min_val, max_val) result(r)
        real(dp), intent(in) :: min_val, max_val
        real(dp) :: r, u

        call random_number(u)
        r = min_val + (max_val - min_val) * u
    end function random_uniform
end module config
