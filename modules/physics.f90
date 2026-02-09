module physics
    use config
    implicit none

contains

    ! Applies PBC to a vector
    pure function apply_pbc(r, L) result(r_pbc)
        real(dp), intent(in) :: r(3)
        real(dp), intent(in) :: L
        real(dp) :: r_pbc(3)

        r_pbc = r - L * dnint(r / L)
    end function apply_pbc

    ! PBC euclidean distance
    pure function get_pbc_dist(r1, r2, L) result(dist)
        real(dp), intent(in) :: r1(3), r2(3)
        real(dp), intent(in) :: L
        real(dp) :: dist, dr(3)

        dr = apply_pbc(r1 - r2, L)
        dist = sqrt(dot_product(dr, dr))
    end function get_pbc_dist

    ! Coulomb total potential
    function get_total_potential(positions, L) result(V)
        real(dp), intent(in) :: positions(3, N)
        real(dp), intent(in) :: L
        real(dp) :: V, r
        integer :: i, j

        V = 0.0_dp
        do i = 1, N-1
            do j = i+1, N
                r = get_pbc_dist(positions(:, i), positions(:, j), L)
                if (r > 1.0e-12_dp) then ! Evitar división por cero
                    V = V + (q * q) / r
                endif
            end do
        end do
    end function get_total_potential

    ! Coulomb potential for a single particle
    function get_single_particle_potential(pid, pos_vec, all_positions, L) result(V)
        integer, intent(in) :: pid          ! particle id
        real(dp), intent(in) :: pos_vec(3)  ! position of pid
        real(dp), intent(in) :: all_positions(3, N)
        real(dp), intent(in) :: L
        real(dp) :: V, r
        integer :: i

        V = 0.0_dp
        do i = 1, N
            if (i /= pid) then
                r = get_pbc_dist(pos_vec, all_positions(:, i), L)
                if (r > 1.0e-12_dp) then
                    V = V + (q * q) / r
                endif
            end if
        end do
    end function get_single_particle_potential
end module physics
