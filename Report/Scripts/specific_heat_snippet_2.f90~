subroutine calc_specific_heat(end_of_routine,N_part, kin_energy, sum_kin_energy, sum_kin_energy_sqr, step)
	....

	if (end_of_routine .eqv. .true.) then
		specific_heat = ((2d0/(3d0*N_part)) - (((sum_kin_energy_sqr * step) - sum_kin_energy**2) / (sum_kin_energy**2)))**(-1)

		print *, "The specific heat is ", specific_heat
	end if		
end subroutine


