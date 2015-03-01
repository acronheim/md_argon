subroutine write_histogram_file(average_number, num_intervals, N_part, step )

	temp_factor = 2d0 * num_intervals / N_part
	temp_factor2 = 1d0 * num_intervals / (N_part - 1)
	temp_factor3 = 1d0 * num_intervals / step
	constant_factor = (temp_factor * temp_factor2 * temp_factor3) / ( 4 * abs(atan(1d0)) * 4)

	open (unit=6,file="histogram.dat",action="write")
	do i=1,num_intervals
		write (6,"(I3, 4F18.6)")  i, (constant_factor * average_number(i) )/ (i**2)
	end do 
end subroutine
