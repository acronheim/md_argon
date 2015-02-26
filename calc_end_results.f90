!Test file to make makefile compile two executables
		

program calc_end_results
!	use argon_box_results
	implicit none	

	integer :: N_part, step
	real(8) :: L_side
	

	print *, "Reading constants from argon_box program"
	open (unit =4, file = "constants.dat")
	read(4,*) step, N_part, L_side

	print *, "Number of steps are: ", step
	print *, "Number of particles are: ", N_part
	print *, "Length of box is: ", L_side
	


		print *,  "Calculating position correlation function"

!		call write_pos_correlation(L_side, 10, pos, N_part, step)
		



end program
