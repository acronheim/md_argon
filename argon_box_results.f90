module argon_box_results
	implicit none
  	private
	
	public write_energy_file

contains

	subroutine write_energy_file(H, kin_energy, pot_energy, T, cnt)
		real(8), intent(in) :: H, kin_energy, T, pot_energy
		integer, intent(in) :: cnt
		open (unit=1,file="energy_matrix.dat",action="write")
		write (1,"(I6, 4F18.6)")  cnt, H, kin_energy, pot_energy, T
	end subroutine

	subroutine calc_specific_heat
		



	end subroutine
	
end module
