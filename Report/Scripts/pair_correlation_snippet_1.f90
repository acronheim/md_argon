subroutine calc_dynamics(calc_quant, N_part, L_side, time_step, m, e, s, r_cut, pos, kin_energy, pot_energy,&
                                & virial, vel, num_intervals, histogram_vector)
		
	histogram_vector = 0
			
	do n = 1,N_part 	
		do i = 	1,N_part !integrate over all particles inside box except i = n					
			do j = -1, 1 
			do k = -1, 1 !periodic boundary condition
			do l = -1, 1				
				if (n/=i) then 
					r_vec = .... !Calculation of the difference vector
					r = sqrt(dot_product(r_vec, r_vec))

					! histogram for the pair correlation function 
					if ((n > i) .and. (calc_quant .eqv. .true.)) then
						hist_i = 1 + floor(r/delta_r_hist)
						if (hist_i < num_intervals + 1) then ! defines a cut off distance
							histogram_vector(hist_i) = histogram_vector(hist_i) + 1
						end if
					end if
				end if		
			end do 
			end do 
			end do
		end do		
	end do	
		
end subroutine

