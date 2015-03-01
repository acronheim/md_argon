module argon_box_results
	implicit none
  	private
	
	public calc_specific_heat

contains

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
