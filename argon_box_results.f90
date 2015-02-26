module argon_box_results
	implicit none
  	private
	
	public write_energy_file, calc_specific_heat, write_pos_correlation, pos_correlation, write_histogram_file

contains

	subroutine write_energy_file(H, kin_energy, pot_energy, T, cnt)
		real(8), intent(in) :: H, kin_energy, T, pot_energy
		integer, intent(in) :: cnt
		open (unit=1,file="energy_matrix.dat",action="write")
		write (1,"(I6, 4F18.6)")  cnt, H, kin_energy, pot_energy, T
	end subroutine

	subroutine write_histogram_file(average_number, num_intervals, N_part, step )
		integer, intent(in) :: num_intervals, N_part, step
		integer, intent(in), dimension(1:num_intervals) :: average_number
		real(8) :: constant_factor, temp_factor
		integer :: i

		temp_factor = num_intervals**3/(N_part * (N_part - 1))
		constant_factor = temp_factor / ( step * 4 * abs(atan(1d0)) * 4)


		open (unit=6,file="histogram.dat",action="write")
		do i=1,num_intervals

			write (6,"(I3, 4F18.6)")  i, (constant_factor * average_number(i) )/ (i**2)
!			print *, num_intervals, N_part, step, average_number(i), i
!I3, 4F18.6

		end do 


	end subroutine


	subroutine calc_specific_heat(end_of_routine,N_part, kin_energy, sum_kin_energy, sum_kin_energy_sqr)
		logical, intent(in) :: end_of_routine
		real(8), intent(in) :: kin_energy
		real(8), intent(out) :: sum_kin_energy, sum_kin_energy_sqr
		integer, intent(in) :: N_part	
		real(8) :: specific_heat

		sum_kin_energy_sqr = sum_kin_energy_sqr +  kin_energy**2
		sum_kin_energy = sum_kin_energy + kin_energy

		if (end_of_routine .eqv. .true.) then
			specific_heat = ((2d0/(3d0*N_part)) - ((sum_kin_energy_sqr - sum_kin_energy**2) / (sum_kin_energy**2)))**(-1)

			print *, "The specific heat is ", specific_heat
		end if
		
	end subroutine

	subroutine write_pos_correlation(L_side, num_points, pos, N, step)
		real(8), intent(in) :: L_side
		integer, intent(in) :: num_points, N, step
		real(8), intent(in), dimension(1:3, 1:N) :: pos
		integer :: i, j, k
		real(8) :: delta_x, average_sum_correlation, absolute_delta_r
		real(8), dimension(1:3) :: delta_r
	
		delta_x = L_side/num_points
		average_sum_correlation = 0d0

		open (unit=2,file="histogram_correlation.dat",action="write")
		open (unit=3,file="delta_r_file.dat",action="write")
		open (unit=4,file="ijk.dat",action="write")
				

		do i=1,num_points
			do j=1,num_points
				do k=1, num_points
					delta_r(1) = i * delta_x
					delta_r(2) = j * delta_x
					delta_r(3) = k * delta_x
					absolute_delta_r = sqrt(delta_r(1)**2 + delta_r(2)**2 + delta_r(3)**2)
					call pos_correlation(pos, N, delta_r, L_side, average_sum_correlation, step)
					write (2,"(4F18.6)") absolute_delta_r, average_sum_correlation
					write (3,"(4F18.6)") delta_r(1), delta_r(2), delta_r(3), absolute_delta_r
					write (4,"(I6)") i, j, k
				end do 
			end do
		end do

	end subroutine


	subroutine pos_correlation(pos, N, delta_r, L_side, average_sum_correlation, step)
		real(8), intent(in) :: L_side
		real(8), intent(out) :: average_sum_correlation
		integer,intent(in) :: N, step
		real(8), intent(in), dimension(1:3, 1:N) :: pos
		real(8), intent(in), dimension(1:3) :: delta_r
		integer :: i,j
		real(8) :: small_epsilon, sum_correlation, difference_in_x, difference_in_y, difference_in_z

		small_epsilon =  1d-5* L_side/N
		sum_correlation = 0d0
		
		do i=1,N
			do j=1,N
				difference_in_x = pos(1,i) + delta_r(1) - pos(1,j)
				difference_in_y = pos(2,i) + delta_r(2) - pos(2,j)
				difference_in_z = pos(3,i) + delta_r(3) - pos(3,j)
				if (i /= j) then
				if  (difference_in_x .lt. small_epsilon) then
				if (difference_in_x .lt. small_epsilon) then
				if  (difference_in_z .lt. small_epsilon) then
					sum_correlation = sum_correlation + 1			
				end if
				end if
				end if
				end if
			end do
		end do

		sum_correlation = (sum_correlation * L_side**3)/(N*(N-1))

		average_sum_correlation = (average_sum_correlation * (step-1) + sum_correlation) / step

		

	end subroutine

end module
