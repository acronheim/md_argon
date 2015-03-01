subroutine init_vel(T, Kb, m, N_part, vel)
		
	do n = 1,N_part
		do i = 1,3
			CALL RANDOM_NUMBER(xs(1))
			CALL RANDOM_NUMBER(xs(2))						
			vel(i,n) = sqrt(Kb*T/m) * sqrt(-2d0*log(xs(1)))*cos(2*pi*xs(2)) !sigma * box_muller
		end do
	end do
	
	.....
end subroutine

