program heating
    use config
    use physics
    use montecarlo
    use io

    implicit none

    real(dp) :: L
    real(dp) :: total_energy
    real(dp) :: T, C, C0, dr
    real(dp) :: acceptance_ratio
    real(dp) :: sum_E, sum_E_sq, avg_E, heat_cap, gamma
    real(dp), allocatable :: positions(:, :)

    integer :: mc_steps, warm_up_steps, step, i, t_idx, num_temps
    integer :: energy_unit, heat_capacity_unit, accepted_moves
    
    ! Input file containing the equilibrated crystal structure (from annealing)
    character(len=50) :: input_crystal_file = 'final_position.xyz'

    ! Read parameters. Note: 'input_heating.in' should define a heating schedule
    ! (T_i = low/solid, T_f = high/liquid, T_step > 1.0)
    call read_input('input_heating.in')

    ! Box length calculation based on density
    L = (real(N, dp) / density)**(1.0_dp/3.0_dp)

    mc_steps = 100 * N**2
    warm_up_steps = int(0.1_dp * mc_steps)
    
    ! Base scalar for displacement adjustments
    C0 = 0.8_dp * L

    ! Calculate total number of temperature steps for geometric heating
    num_temps = nint(log(T_f / T_i) / log(T_step)) + 1

    allocate(positions(3, N))

    ! ---------------------------------------------------------
    ! Load the pre-equilibrated crystal structure
    ! ---------------------------------------------------------
    print *, "Loading initial crystal structure..."
    call read_positions_xyz(input_crystal_file, positions, N)

    ! Calculate initial potential energy of the loaded crystal
    total_energy = get_total_potential(positions, L)

    print *, "========================================="
    print *, " HEATING RUN "
    print *, " N =", N
    print *, " T_initial =", T_i
    print *, " T_final   =", T_f
    print *, " T_step    =", T_step
    print *, "========================================="

    open (newunit = energy_unit, file = 'heating_energy.out', action = 'write', status = 'replace')
    open (newunit = heat_capacity_unit, file = 'heating_cv.out', action = 'write', status = 'replace')

    T = T_i ! Start from the cold (solid) state

    ! Main Heating Loop
    do t_idx = 1, num_temps
        accepted_moves = 0
        
        ! Scale max displacement (dr) with sqrt(T) to maintain reasonable acceptance 
        ! as the system melts and particles diffuse faster.
        C = C0 * (T / T_i)**alpha 
        dr = C * sqrt(T) 
        
        sum_E = 0.0_dp
        sum_E_sq = 0.0_dp

        print '(A, F10.5, A, F10.5)', " * Heating T =", T, " | dr =", dr

        ! 1. Equilibration Phase (Let system adjust to new T)
        do step = 1, warm_up_steps
            call montecarlo_step(L, T, dr, positions, total_energy, accepted_moves)
        end do

        ! 2. Production Phase (Data Collection)
        accepted_moves = 0
        do step = 1, mc_steps
            call montecarlo_step(L, T, dr, positions, total_energy, accepted_moves)

            ! Sample energy periodically to reduce correlation
            if ( mod(step, 100) == 0 ) then
                sum_E = sum_E + total_energy
                sum_E_sq = sum_E_sq + (total_energy**2)
            end if
        end do

        ! Statistics Calculation
        ! Note: Ensure consistent type conversion for the average
        avg_E = sum_E / real(mc_steps/100, dp)
        
        ! Heat Capacity (Cv) calculation via energy fluctuations: <E^2> - <E>^2 / T^2
        heat_cap = ( (sum_E_sq / real(mc_steps/100, dp)) - (avg_E**2) ) / (T**2)
        
        ! Coulomb coupling parameter
        gamma = (4.0_dp/3.0_dp * pi * density)**(1.0_dp / 3.0_dp) / T
        
        write(unit = heat_capacity_unit, fmt = *) T, avg_E, heat_cap, gamma

        acceptance_ratio = real(accepted_moves, dp) / real(mc_steps, dp)
        
        ! Increase temperature (Geometric heating)
        T = T * T_step 
    end do

    close(unit = energy_unit)
    close(unit = heat_capacity_unit)
    deallocate(positions)
    
end program heating