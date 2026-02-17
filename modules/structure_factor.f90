module sq_module
    implicit none
    integer, parameter :: dp = kind(0.d0) ! Doble precisión (float64)

contains

    ! Subrutina optimizada con OpenMP
    subroutine calculate_frame_sq(x, y, z, qx, qy, qz, n_atoms, n_q, sq_out)
        ! Argumentos de entrada
        integer, intent(in) :: n_atoms, n_q
        real(dp), intent(in) :: x(n_atoms), y(n_atoms), z(n_atoms)
        real(dp), intent(in) :: qx(n_q), qy(n_q), qz(n_q)
        
        ! Argumento de salida
        real(dp), intent(out) :: sq_out(n_q)
        
        ! Variables locales
        integer :: i, j
        real(dp) :: fase, sum_cos, sum_sin
        
        ! Paralelización: Distribuimos los vectores q entre los núcleos de la CPU
        !$omp parallel do private(i, j, fase, sum_cos, sum_sin)
        do i = 1, n_q
            sum_cos = 0.0_dp
            sum_sin = 0.0_dp
            
            ! Bucle sobre átomos (vectorizado por el compilador)
            do j = 1, n_atoms
                fase = qx(i)*x(j) + qy(i)*y(j) + qz(i)*z(j)
                sum_cos = sum_cos + cos(fase)
                sum_sin = sum_sin + sin(fase)
            end do
            
            ! S(q) = |sum(exp(iqr))|^2 / N
            sq_out(i) = (sum_cos**2 + sum_sin**2) / real(n_atoms, dp)
        end do
        !$omp end parallel do

    end subroutine calculate_frame_sq

end module sq_module