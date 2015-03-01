subroutine init_vel(T, Kb, m, N_part, vel)
	.....	
	!Set center of mass velocity to zero
	do i = 1,3 
		vel(i,:) = vel(i,:) - sum(vel(i,:))/N_part
	end do
end subroutine

