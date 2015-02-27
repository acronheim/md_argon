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
		real(8) :: constant_factor, temp_factor, temp_factor2, temp_factor3
		integer :: i

		temp_factor = 1d0 * num_intervals / N_part
		temp_factor2 = 1d0 * num_intervals / (N_part - 1)
		temp_factor3 = 1d0 * num_intervals / step
		constant_factor = (temp_factor * temp_factor2 * temp_factor3) / ( 4 * abs(atan(1d0)) * 4)


		open (unit=6,file="histogram.dat",action="write")
		do i=1,num_intervals

			write (6,"(I3, 4F18.6)")  i, (constant_factor * average_number(i) )/ (i**2)

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


end module
