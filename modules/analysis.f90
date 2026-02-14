module radial_distribution_function
    use config

    implicit none

    contains
        ! Read data from the positions file and store it into arrays.
        ! Particle positions in time are stored in columns: time x y z
        ! Each row refers to a particle.
        subroutine read_trajectory(positions_file, x, y, z, time, freeze_mc_steps)
            implicit none

            character(*), intent(in) :: positions_file
            integer, intent(in) :: freeze_mc_steps
            real(kind = dp), allocatable, intent(inout) :: x(:, :), y(:, :), z(:, :), time(:)

            character(3) :: atom_type
            integer :: part, step, ios, n_
            real(kind = dp) :: t, temp

            print *, 'Reading trajectory data from ', positions_file

            open(4, file = positions_file, status = 'old', action = 'read')

            do step = 1, freeze_mc_steps
                read(4, *, iostat = ios) n_
                read(4, *, iostat = ios) t, temp

                do part = 1, N
                    read(4, *, iostat = ios) &
                        atom_type, x(part, step), y(part, step), z(part, step)

                    ! Store the different values of t (once for each different time).
                    if (part == 1) then
                        time(step) = t
                    end if
                end do
            end do

            close(4)
        end subroutine read_trajectory

        ! Compute RDF using the stored data.
        ! Must be executed after read_trajectory, as it depends on the arrays it
        ! creates.
        subroutine compute_rdf(x, y, z, L, freeze_mc_steps, rdf_file)
            implicit none

            real(kind = dp), allocatable, intent(in) :: x(:, :), y(:, :), z(:, :)
            real(kind = dp), intent(in) :: L
            integer, intent(in) :: freeze_mc_steps
            character(*), intent(in) :: rdf_file

            integer :: i, j, k, time_index, bins, bin_index
            real(kind = dp), parameter :: dr = 0.01_dp
            real(kind = dp), allocatable :: h(:), rdf(:), r_values(:)
            real(kind = dp) :: maximum_radius, r, r_sq, dx, &
                dy, dz, dv, r_lo, r_hi, const, nid

            ! Parameters
            maximum_radius = L/2.0_dp
            bins = int(maximum_radius / dr)

            ! Allocate arrays
            allocate(h(bins), rdf(bins), r_values(bins))
            h = 0.0_dp

            ! Compute histogram h(k)
            do time_index = 1, freeze_mc_steps
                do i = 1, N - 1
                    do j = i + 1, N
                        
                        ! Compute minimum image distance between particles
                        dx = x(j, time_index) - x(i, time_index)
                        dy = y(j, time_index) - y(i, time_index)
                        dz = z(j, time_index) - z(i, time_index)
                        
                        ! PBC
                        dx = dx - L * anint(dx / L)
                        dy = dy - L * anint(dy / L)
                        dz = dz - L * anint(dz / L)

                        r_sq = dx**2 + dy**2 + dz**2
                        r = sqrt(r_sq)

                        if (r < maximum_radius) then
                            bin_index = floor(r / dr) + 1

                            ! Boundary check to prevent out-of-bounds access
                            if (bin_index >= 1 .and. bin_index <= bins) then
                                h(bin_index) = h(bin_index) + 2  ! pairwise counting
                            else
                                print *, 'Warning: bin_index out of range:', bin_index, ' (max bins:', bins, ')'
                            endif
                        endif
                    end do
                end do
            end do

            ! Normalize RDF
            ! const = 4.0_dp*pi/3.0_dp
            ! do k = 1, bins
            !     r_lo = (k - 1) * dr
            !     r_hi = r_lo + dr
            !     dv = const * (r_hi**3 - r_lo**3)  ! Shell volume
            !     nid = density * dv
            !     rdf(k) = h(k) / (N * freeze_mc_steps * nid)
            !     r_values(k) = (k - 0.5_dp) * dr  ! Bin center
            ! end do
            const = 4.0_dp*pi/3.0_dp

            do k = 1, bins
                r_lo = (k - 1) * dr
                r_hi = r_lo + dr
                dv = const * (r_hi**3 - r_lo**3)
                
                ! USA LA DENSIDAD REAL:
                nid = real(N, dp) / (L**3)  * dv  
                
                ! Asegura que N y steps se traten como reales en la división
                rdf(k) = h(k) / (real(N, dp) * real(freeze_mc_steps, dp) * nid)
                
                r_values(k) = (k - 0.5_dp) * dr
            end do

            ! Save RDF results to file
            open(12, file = rdf_file, status = 'replace')
            do k = 1, bins
                write(12, *) r_values(k), rdf(k)
            end do
            close(12)

            print *, 'RDF calculation completed and saved to ', rdf_file

            deallocate(h, rdf, r_values)
        end subroutine compute_rdf

        ! subroutine compute_sq(r_values, rdf, L, density, sq_file)
        !     implicit none
        !     real(kind=dp), intent(in) :: r_values(:), rdf(:), L, density
        !     character(*), intent(in) :: sq_file
            
        !     integer :: i, j, n_q, bins
        !     real(kind=dp) :: q, dq, sum_int, r, gr, s_val
        !     real(kind=dp), parameter :: q_max = 20.0_dp
            
        !     n_q = 500
        !     dq = q_max / n_q
        !     bins = size(r_values)

        !     open(13, file=sq_file, status='replace')

        !     do i = 1, n_q
        !         q = i * dq
        !         sum_int = 0.0_dp
                
        !         do j = 1, bins
        !             r = r_values(j)
        !             gr = rdf(j)
                    
        !             if (r > L/2.0_dp) exit
        !             sum_int = sum_int + (r**2 * (gr - 1.0_dp) * (sin(q*r)/(q*r))) * dr
        !         end do
                
        !         s_val = 1.0_dp + 4.0_dp * pi * density * sum_int
        !         write(13, *) q, s_val
        !     end do
            
        !     close(13)
        !     print *, 'S(q) calculation completed and saved to ', sq_file
        ! end subroutine compute_sq

end module radial_distribution_function