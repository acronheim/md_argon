module argon_box_dynamics
	implicit none
	private
	
	public calc_dynamics, rescale_vel, new_pos
	
contains
	
	subroutine calc_dynamics(calc_quant, N_part, L_side, time_step, m, e, s, r_cut, pos, kin_energy, pot_energy, virial, vel)
		
		logical, intent(in) :: calc_quant !for improving efficiency with velocity verlet method
		integer, intent(in) :: N_part
		real(8), intent(in) :: e, s, r_cut !Lennard Jones
		real(8), intent(in) :: m, time_step, L_side
		real(8), intent(inout), dimension(1:3, 1:N_part) :: vel
		real(8), intent(in), dimension(1:3, 1:N_part) :: pos
		integer :: i,j,k,l,n		
		real(8), intent(out) :: kin_energy, pot_energy, virial
		real(8) :: sum_v_2, F(3), dF(3), r, r_vec(3)
		
		virial = 0
		pot_energy = 0 
		sum_v_2 = 0
			
		do n = 1,N_part 	
			F = 0
			do i = 	1,N_part !integrate over all particles inside box except i = n					
				do j = -1, 1 
				do k = -1, 1 !periodic boundary condition
				do l = -1, 1				
					if (n/=i) then 
						r_vec = (/pos(1,n)-pos(1,i), pos(2,n)-pos(2,i), pos(3,n)-pos(3,i)/) + L_side*(/j,k,l/)
						r = sqrt(dot_product(r_vec, r_vec))  
						if (r<r_cut) then
							dF = e*(48*s**12/r**14 - 24*s**6/r**8) * r_vec
							F = F + dF
							if (calc_quant .eqv. .true.) then
								if (n > i) then	
									pot_energy = pot_energy + 4*e*((s/r)**12-(s/r)**6)
									virial =  virial + dot_product(r_vec, dF) 
								end if	
							end if
						end if
					end if		
				end do 
				end do 
				end do
			end do		
			vel(:,n) = vel(:,n) + F/m*time_step/2 ! velocity verlet method -> factor 1/2 !
			if (calc_quant .eqv. .true.) then	
				sum_v_2 = sum_V_2 + dot_product(vel(:,n),vel(:,n))
			end if
		end do	
		kin_energy = m/2*sum_v_2	
		
	end subroutine

	subroutine rescale_vel(T_intended, kin_energy, Kb, N_part, Vel)
		!rescale velocities in order to keep temperature constant
		
		real(8), intent(in) :: T_intended, kin_energy, kb
		integer, intent(in) :: N_part
		real(8), intent(inout), dimension(1:3, 1:N_part) :: vel
		real(8) :: scaling_factor
		
		scaling_factor = sqrt((N_part - 1)*3/2*kb*T_intended/kin_energy)
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
		end do	
		
	end subroutine
	
end module 
