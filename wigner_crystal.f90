module config
    use iso_fortran_env, only: dp => real64
    implicit none

    ! Parameters
    integer, parameter :: N = 128
    
    real(dp), parameter :: density = 1.0_dp
    real(dp), parameter :: q = 1.0_dp
    
    ! Mathematical constants
    real(dp), parameter :: PI = 4.0_dp * atan(1.0_dp)

contains

    ! Uniform Random generator [min, max]
    function random_uniform(min_val, max_val) result(r)
        real(dp), intent(in) :: min_val, max_val
        real(dp) :: r, u
        call random_number(u)
        r = min_val + (max_val - min_val) * u
    end function random_uniform

end module config

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

module montecarlo
    use physics
    use config

contains

    subroutine montecarlo_step(L, T, positions, total_energy)
        implicit none
        real(dp), intent(in) :: L, T
        real(dp), intent(inout), allocatable :: positions(:, :)
        real(dp), intent(inout) :: total_energy

        real(dp) :: rnd_val, dr, dE, boltzmann_factor, current_E_part, new_E_part
        real(dp) :: rvec(3), old_pos_vec(3), new_pos_vec(3)
        integer :: particle_id
        

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
        call random_number(rnd_val)
        boltzmann_factor = exp(-dE / T)

        if (dE < 0.0_dp .or. rnd_val < boltzmann_factor) then
            ! Movement accepted
            positions(:, particle_id) = new_pos_vec
            total_energy = total_energy + dE
        end if

    end subroutine montecarlo_step

end module montecarlo

program annealing
    use config
    use physics
    implicit none

    ! Simulation variables
    real(dp) :: L
    real(dp), allocatable :: positions(:,:) ! Shape (3, N)
    real(dp) :: total_energy, current_E_part, new_E_part, dE
    real(dp) :: T, T_i, T_f, T_step, C, C0, alpha, dr
    real(dp) :: acceptance_ratio
    real(dp) :: rnd_val, boltzmann_factor
    real(dp) :: rvec(3), old_pos_vec(3), new_pos_vec(3)
    
    integer :: mc_steps, warm_up_steps, step, i, particle_id
    integer :: num_temps, t_idx, positions_unit, energy_unit
    integer :: accepted_moves

    ! Parameter inicialization
    L = (real(N, dp) / density)**(1.0_dp/3.0_dp)
    
    mc_steps = 100 * N**2
    warm_up_steps = int(0.1_dp * mc_steps)
    
    T_i = 2.0_dp
    T_f = 0.0001_dp
    T_step = 0.98_dp
    
    C0 = 0.8_dp * L
    alpha = 0.5_dp

    ! Number of temperature steps
    num_temps = nint(log(T_f / T_i) / log(T_step)) + 1


    allocate(positions(3, N))
    call random_seed() 

    ! Random initialization of the particles
    do i = 1, N
        positions(1, i) = random_uniform(-L/2.0_dp, L/2.0_dp)
        positions(2, i) = random_uniform(-L/2.0_dp, L/2.0_dp)
        positions(3, i) = random_uniform(-L/2.0_dp, L/2.0_dp)
    end do

    ! Initial energy
    total_energy = get_total_potential(positions, L)

    print *, "========================================="
    print *, " Simulation Parameters"
    print *, " N =", N
    print *, " L =", L
    print *, " MC Steps per T =", mc_steps
    print *, "========================================="

    ! Energy unit
    open (newunit = energy_unit, file = 'energy.out', action = 'write', status = 'replace')
    

    ! Loop temperature (annealing)
    T = T_i
    do t_idx = 1, num_temps
        
        accepted_moves = 0
        C = C0 * (T / T_i)**alpha
        dr = C * sqrt(T)
        
        print '(A, F10.5, A, F10.5)', " * Temperature =", T, " | C =", C

        ! Warm up MC ---------------------------------------------------------------------------------------------------
        do step = 1, warm_up_steps
            
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
            call random_number(rnd_val)
            boltzmann_factor = exp(-dE / T)

            if (dE < 0.0_dp .or. rnd_val < boltzmann_factor) then
                ! Movement accepted
                positions(:, particle_id) = new_pos_vec
                total_energy = total_energy + dE
            end if

        end do

        ! Loop MC
        do step = 1, mc_steps
            
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
            call random_number(rnd_val)
            boltzmann_factor = exp(-dE / T)

            if (dE < 0.0_dp .or. rnd_val < boltzmann_factor) then
                ! Movement accepted
                positions(:, particle_id) = new_pos_vec
                total_energy = total_energy + dE
                accepted_moves = accepted_moves + 1
            end if

            if ( mod(step, 100) == 0 ) then
                write(unit = energy_unit, fmt = *) step, total_energy, T
            end if

        end do

        acceptance_ratio = real(accepted_moves, dp) / real(mc_steps, dp)
        
        print '(A, E14.6)', "   -> Final Energy: ", total_energy
        print '(A, F10.4)', "   -> Accept Ratio: ", acceptance_ratio
        print * 
        T = T * T_step ! Update temperature
    end do

    close(unit  = energy_unit)


    open (newunit = positions_unit, file = 'positions.xyz', action = 'write', status = 'replace')
    write(unit = positions_unit, fmt = '(i0/)') N
    do i = 1, N
        write(unit = positions_unit, fmt = *) 'A', positions(:, i)
    end do
    close(unit  = positions_unit)
    

    deallocate(positions)

end program annealing