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

    subroutine init_bcc_lattice(positions, N, L)
        use iso_fortran_env, only: dp => real64
        implicit none
        real(dp), intent(out) :: positions(3, N)
        integer, intent(in) :: N
        real(dp), intent(in) :: L
        integer :: n_cells, ix, iy, iz, idx
        real(dp) :: a, offset

        ! Check for N=128 (4x4x4 BCC)
        n_cells = nint((real(N, dp)/2.0_dp)**(1.0_dp/3.0_dp))
        a = L / real(n_cells, dp)
        offset = -L / 2.0_dp
        idx = 1

        do ix = 0, n_cells - 1
            do iy = 0, n_cells - 1
                do iz = 0, n_cells - 1
                    ! Corner Atom
                    positions(1, idx) = offset + ix*a
                    positions(2, idx) = offset + iy*a
                    positions(3, idx) = offset + iz*a
                    idx = idx + 1
                    ! Body Center Atom
                    positions(1, idx) = offset + (ix + 0.5_dp)*a
                    positions(2, idx) = offset + (iy + 0.5_dp)*a
                    positions(3, idx) = offset + (iz + 0.5_dp)*a
                    idx = idx + 1
                end do
            end do
        end do
    end subroutine init_bcc_lattice
end module physics
