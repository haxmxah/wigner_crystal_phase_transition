module montecarlo
    use physics
    use config

contains

    subroutine montecarlo_step(L, T, dr, positions, total_energy, accepted_moves)
        implicit none
        real(dp), intent(in) :: L, T, dr
        real(dp), intent(inout) :: positions(3, N)
        real(dp), intent(inout) :: total_energy
        integer, intent(inout) :: accepted_moves

        real(dp) :: rnd_val, dE, current_E_part, new_E_part
        real(dp) :: rvec(3), old_pos_vec(3), new_pos_vec(3)
        integer :: particle_id
        logical :: move_accepted

        move_accepted = .false.

        ! Pick a particle randomly
        call random_number(rnd_val)
        particle_id = floor(rnd_val * N) + 1
        if (particle_id > N) particle_id = N

        ! Generate a random displacement
        call random_number(rvec) ! rvec en [0,1]
        rvec = (2.0_dp * rvec - 1.0_dp) * dr ! rvec en [-dr, dr]

        ! Propose new position applying pbc
        old_pos_vec = positions(:, particle_id)
        new_pos_vec = apply_pbc(old_pos_vec + rvec, L)

        ! Compute the energy difference due to that particle
        current_E_part = get_single_particle_potential(particle_id, old_pos_vec, positions, L)
        new_E_part     = get_single_particle_potential(particle_id, new_pos_vec, positions, L)

        dE = new_E_part - current_E_part

        ! Metropolis algorithm
        if (dE < 0.0_dp) then
            ! Accept
            move_accepted = .true.
        else
            call random_number(rnd_val)
            if (rnd_val < exp(-dE / T)) then
                ! Accept
                move_accepted = .true.
            end if
        end if

        if (move_accepted) then
            ! Update position and total energy
            positions(:, particle_id) = new_pos_vec
            total_energy = total_energy + dE
            accepted_moves = accepted_moves + 1
        end if
    end subroutine montecarlo_step
end module montecarlo
