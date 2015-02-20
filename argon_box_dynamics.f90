module argon_box_dynamics
	implicit none
	private
	
	public calc_dynamics, rescale_vel, new_pos
	
contains
	
	subroutine calc_dynamics(N_part, L_side, time_step, Kb, m, e, s, r_cut, pos, kin_energy, pot_energy, virial, vel)
		integer, intent(in) :: N_part
		real(8), intent(in) :: e, s, r_cut !Lennard Jones
		real(8), intent(in) :: m, time_step, L_side, Kb
		real(8), intent(inout), dimension(1:3, 1:N_part) :: vel
		real(8), intent(in), dimension(1:3, 1:N_part) :: pos
		integer :: i,j,k,l,n		
		real(8), intent(out) :: kin_energy, pot_energy, virial
		real(8) :: sum_v_2, F(3), dF(3), r, r_vec(3)
		virial = 0
		pot_energy = 0 !potential		
		do n = 1,N_part 	
			F = 0
			sum_v_2 = 0
			do i = 	1,N_part !integrate over all particles inside box except i = n					
				do j = -1, 1 
				do k = -1, 1 !periodic boundary condition for potential forces..
				do l = -1, 1				
					if (n/=i) then !(.not.(n==i .and. ((/j,k,l/) == (/0,0,0/)))) then 
						r_vec = (/pos(1,n)-pos(1,i), pos(2,n)-pos(2,i), pos(3,n)-pos(3,i)/) + L_side*(/j,k,l/)
						r = sqrt(dot_product(r_vec, r_vec))  
						if (r<r_cut) then
							dF = e*(48*s**12/r**14 - 24*s**6/r**8) * r_vec
							F = F + dF
							if (n > i) then	
								pot_energy = pot_energy + 4*e*((s/r)**12-(s/r)**6)
							end if
							virial =  virial + dot_product(r_vec, dF)	 
						end if
					end if		
				end do 
				end do 
				end do
			end do			
			vel(:,n) = vel(:,n) + F/m*time_step	
			sum_v_2 = sum_V_2 + dot_product(vel(:,n),vel(:,n))
			!print *,"particle:", n, "/",  N_part, "velocity:", (VEL(i,n), i=1,3)
		end do		
		kin_energy = m/2*sum_v_2	
	end subroutine

	subroutine rescale_vel(T_intended, T_actual, time_step, N_part, Vel)
		!rescale velocities in order to keep temperature constant
		real(8), intent(in) :: T_intended, T_actual, time_step
		integer, intent(in) :: N_part
		real(8), intent(inout), dimension(1:3, 1:N_part) :: vel
		real(8), parameter :: tau = 1d-3
		real(8) :: scaling_factor
		! Berendsen thermostat
		scaling_factor = sqrt(1 + time_step/tau*(T_intended/T_actual -1))
		vel = scaling_factor*vel
	end subroutine
	
	subroutine new_pos(N_part, L_side, dt, vel, pos)
		! Postion calculation
		real(8), intent(in) :: L_side, dt
		integer, intent(in) :: N_part
		real(8), intent(in), dimension(1:3, 1:N_part) :: vel
		real(8), intent(inout), dimension(1:3, 1:N_part) :: pos
		integer :: n, i
		do n = 1,N_part 
			pos(:,n) = pos(:,n) + vel(:,n)*dt			
			do i = 1,3 !implements periodic boundary conditions
				if (pos(i,n) < 0d0) then 
					pos(i,n) = pos(i,n) + L_side
				else if (pos(i,n) > L_side) then
					pos(i,n) = pos(i,n) - L_side 
				end if
			end do		
	!		print *, "particle:", n, "/",  N_part, "position:", pos(:,n)
		end do	
	end subroutine
	
end module 