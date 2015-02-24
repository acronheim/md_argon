module argon_box_results
	implicit none
  	private
	
	public write_energy_file, calc_specific_heat

contains

	subroutine write_energy_file(H, kin_energy, pot_energy, T, cnt)
		real(8), intent(in) :: H, kin_energy, T, pot_energy
		integer, intent(in) :: cnt
		open (unit=1,file="energy_matrix.dat",action="write")
		write (1,"(I6, 4F18.6)")  cnt, H, kin_energy, pot_energy, T
	end subroutine

	subroutine calc_specific_heat(N_part, sum_kin_energy, sum_deltaK_sqr)
		real(8), intent(in) :: sum_kin_energy, sum_deltaK_sqr
		integer, intent(in) :: N_part	
		real(8) :: specific_heat

		specific_heat = ((2d0/(3d0*N_part)) - (sum_deltaK_sqr / (sum_kin_energy**2)))**(-1)
		
		print *, "The specific heat is ", specific_heat
		
	end subroutine

end module
