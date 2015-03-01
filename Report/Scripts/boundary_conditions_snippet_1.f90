subroutine calc_dynamics(calc_quant, N_part, L_side, time_step, m, e, s, r_cut, pos, kin_energy, pot_energy, virial, vel, num_intervals, histogram_vector)
	....

	do n = 1,N_part 	
		F = 0
		do i = 	1,N_part !integrate over all particles inside box except i = n					
			do j = -1, 1 
			do k = -1, 1 !periodic boundary condition
			do l = -1, 1				
				if (n/=i) then 
					r_vec = (/pos(1,n)-pos(1,i), pos(2,n)-pos(2,i), pos(3,n)-pos(3,i)/) + L_side*(/j,k,l/)
					r = sqrt(dot_product(r_vec, r_vec))
					!force calculation
					if (r<r_cut) then
						dF = e*(48*s**12/r**14 - 24*s**6/r**8) * r_vec
						F = F + dF
					end if
				end if		
			end do 
			end do 
			end do
		end do		
	end do	
end subroutine

