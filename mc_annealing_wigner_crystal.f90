program annealing
    use config
    use physics
    use montecarlo
    use io
    use radial_distribution_function

    implicit none

    ! Simulation variables
    real(dp) :: L
    real(dp) :: total_energy, best_energy
    real(dp) :: T, C, C0, dr
    real(dp) :: acceptance_ratio
    real(dp) :: boltzmann_factor
    real(dp) :: old_pos_vec(3), new_pos_vec(3)
    real(dp), allocatable :: positions(:, :), best_positions(:, :), x(:, :), y(:, :), z(:, :), time_points(:) ! Shape (3, N)
    real(dp) :: sum_E, sum_E_sq, avg_E, heat_cap, E_current, gamma

    integer :: mc_steps, warm_up_steps, freeze_mc_steps, step, best_step, i, particle_id
    integer :: num_temps, t_idx, positions_unit, energy_unit, heat_capacity_unit
    integer :: accepted_moves, num_samples
    logical :: melting

    call read_input('input_parameters.in')

    ! Parameter initialization
    L = (real(N, dp) / density)**(1.0_dp/3.0_dp)

    mc_steps = 100 * N**2
    warm_up_steps = int(0.1_dp * mc_steps)
    freeze_mc_steps = int(freeze_mc_steps_scale*mc_steps)

    C0 = 0.8_dp * L

    ! Number of temperature steps
    num_temps = nint(log(T_f / T_i) / log(T_step)) + 1

    allocate(positions(3, N), best_positions(3, N))
    call random_seed()

    ! ! Random initialization of the particles
    ! do i = 1, N
    !     positions(1, i) = random_uniform(-L/2.0_dp, L/2.0_dp)
    !     positions(2, i) = random_uniform(-L/2.0_dp, L/2.0_dp)
    !     positions(3, i) = random_uniform(-L/2.0_dp, L/2.0_dp)
    ! end do

    call init_bcc_lattice(positions, N, L)

    ! print *, "Guardando estado inicial en check_bcc.xyz..."
    
    ! open(newunit=positions_unit, file='check_bcc.xyz', status='replace', action='write')
    ! write(positions_unit, '(I0)') N
    ! write(positions_unit, *) 'Comprobacion Estructura BCC'
    ! do i = 1, N
    !     ! Escribimos 'H' o 'A' como átomo dummy para visualización
    !     write(positions_unit, *) 'A', positions(1, i), positions(2, i), positions(3, i)
    ! end do
    ! close(positions_unit)

    ! Initial energy
    total_energy = get_total_potential(positions, L)
    best_energy = total_energy

    print *, "========================================="
    print *, " Simulation Parameters"
    print *, " N =", N
    print *, " L =", L
    print *, " MC Steps per T =", mc_steps
    print *, "========================================="

    ! Create output file to save results
    open (newunit = energy_unit, file = 'energy.out', action = 'write', status = 'replace')
    open (newunit = heat_capacity_unit, file = 'heat_capacity.out', action = 'write', status = 'replace')

    ! Loop temperature (annealing)
    T = T_i
    do t_idx = 1, num_temps
        accepted_moves = 0
        C = C0 * (T / T_i)**alpha
        dr = C * sqrt(T)
        
        sum_E = 0.0_dp
        sum_E_sq = 0.0_dp

        print '(A, F10.5, A, F10.5)', " * Temperature =", T, " | C =", C

        ! Warm up MC ---------------------------------------------------------------------------------------------------
        do step = 1, warm_up_steps
            call montecarlo_step(L, T, dr, positions, total_energy, accepted_moves)
        end do

        ! Loop MC ---------------------------------------------------------------------------------------------------
        accepted_moves = 0
        num_samples = 0
        do step = 1, mc_steps

            call montecarlo_step(L, T, dr, positions, total_energy, accepted_moves)

            if ( mod(step, 100) == 0 ) then
                write(unit = energy_unit, fmt = *) step, total_energy, T
                
                ! Accumulate sums for averages        
                sum_E = sum_E + total_energy
                sum_E_sq = sum_E_sq + (total_energy**2)
                num_samples = num_samples + 1
                
            end if

        end do

        avg_E = sum_E / real(num_samples, dp)
        heat_cap = ( (sum_E_sq / real(num_samples, dp)) - (avg_E**2) ) / (T**2)
        gamma = (q**2) * ( (4.0_dp/3.0_dp * pi * density)**(1.0_dp / 3.0_dp) ) / T
        acceptance_ratio = real(accepted_moves, dp) / real(mc_steps, dp)

        write(unit = heat_capacity_unit, fmt = *) T, avg_E, heat_cap, gamma, acceptance_ratio


        print '(A, E14.6)', "   -> Final Energy: ", total_energy
        print '(A, F10.4)', "   -> Accept Ratio: ", acceptance_ratio
        print *
        T = T * T_step ! Update temperature
    end do

    close(unit  = energy_unit)
    close(unit  = heat_capacity_unit)

    ! Freeze MC steps ---------------------------------------------------------------------------------------------------    
    open (newunit = positions_unit, file = 'positions_t0.xyz', action = 'write', status = 'replace')
    do step = 1, freeze_mc_steps
        call montecarlo_step(L, T, dr, positions, total_energy, accepted_moves)
        write(unit = positions_unit, fmt = '(i0)') N
        write(unit = positions_unit, fmt = *) step, T
        do i = 1, N
            write(unit = positions_unit, fmt = *) 'A', positions(:, i)
        end do

        ! Saves the energy, step and particle positions for the minimum energy configuration.
        if (total_energy < best_energy) then
            best_energy = total_energy
            best_positions(:, :) = positions(:, :)
            best_step = step
        end if
    end do
    close(unit  = positions_unit)

    open (newunit = positions_unit, file = 'final_position.xyz', action = 'write', status = 'replace')
    write(unit = positions_unit, fmt = '(i0)') N
    write(unit = positions_unit, fmt = '(i0)') best_step
    do i = 1, N
        write(unit = positions_unit, fmt = *) 'A', best_positions(:, i)
    end do
    close(unit  = positions_unit)

    deallocate(positions)

    ! RDF computation ---------------------------------------------------------------------------------------------------
    
    print *
    print *, 'Performing RDF analysis...'


    allocate( &
        x(N, freeze_mc_steps), &
        y(N, freeze_mc_steps), &
        z(N, freeze_mc_steps), &
        time_points(freeze_mc_steps) &
    )

    call read_trajectory('positions_t0.xyz', x, y, z, time_points, freeze_mc_steps)

    ! Compute and save RDF data to an output file.
    call compute_rdf(x, y, z, L, freeze_mc_steps, 'rdf.out')

    deallocate(x, y, z, time_points)

end program annealing
