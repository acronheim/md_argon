subroutine calc_specific_heat(end_of_routine,N_part, kin_energy, sum_kin_energy, sum_kin_energy_sqr)
	....

	sum_kin_energy_sqr = sum_kin_energy_sqr +  kin_energy**2
	sum_kin_energy = sum_kin_energy + kin_energy

	....		
end subroutine


