module argon_box_init
	implicit none
	private
	
	public cubic_fcc_lattice, init_vel, init_random_seed

contains

	function fcc_cell(i) result(output)
		implicit none
		integer, intent(in) :: i
		real(8) :: output(3)
		! face centered cubic unit cell with basis particle positions:
		real(8), dimension(1:3), parameter :: &
		fcc_part1 = (/0d0,  0d0,  0d0/),  &
		fcc_part2 = (/0d0,  5d-1, 5d-1/), &
		fcc_part3 = (/5d-1, 0d0,  5d-1/), &
		fcc_part4 = (/5d-1, 5d-1, 0d0/)		
		real(8), dimension(1:3,1:4), parameter :: &
		R_cell = reshape( (/fcc_part1, fcc_part2, fcc_part3, fcc_part4/), (/3,4/))
		output = R_cell(:,i)
	
		!	print *, "face centered cubic unit cell"
		!	do i = 1,3		!		print *, (R_cell(i,j), j=1,4)
		!	end do		
	end function fcc_cell

	subroutine cubic_fcc_lattice(N_cell_dim, L_side, pos)	
		!initial positions of all particles according to an fcc lattice structure
		integer :: i,j,k,l,n
		integer, intent(in) :: N_cell_dim
		real(8), intent(in) :: L_side
		real(8), intent(out), dimension(1:3, 1:(4*N_cell_dim**3)) :: pos

		n = 0
		do i = 1,N_cell_dim
			do j = 1,N_cell_dim
				do k = 1,N_cell_dim
					do l = 1,4
						n = n + 1
						pos(:,n) = L_side/N_cell_dim*((/ i-1, j-1, k-1 /) + fcc_cell(l))
						!			print *, "particle:", n, "/",  N_part, "position:", pos(:,n)
					end do
				end do
			end do
		end do
	end subroutine

	subroutine init_vel(T, Kb, m, N_part, vel)
		!initial particles velocities according to the maxwell distribution	
		! in the maxwell boltzman distribution each velocity component is normally distributed:
		! Box muller transform used for converting uniform dist to normal dist
		!
		! f(v) = sqrt(m/(2*PI*Kb*T)) * exp(-(v**2)/2 *m/(*Kb*T))
		! sigma**2 = Kb*T/m, and zero mean	
		real(8), intent(in)  :: T, Kb, m
		integer, intent(in) :: N_part
		real(8), intent(out), dimension(1:3, 1:N_part) :: vel
		real(8), parameter :: pi = 4*atan(1d0)
		real(8) :: xs(2) !two random numbers
		integer :: n, i
		
		do n = 1,N_part
			do i = 1,3
				CALL RANDOM_NUMBER(xs(1))
				CALL RANDOM_NUMBER(xs(2))						
				vel(i,n) = sqrt(Kb*T/m) * sqrt(-2d0*log(xs(1)))*cos(2*pi*xs(2)) !sigma * box_muller
			end do
		end do
		!Set center of mass velocity to zero
		do i = 1,3 
			vel(i,:) = vel(i,:) - sum(vel(i,:))/N_part
		end do
	end subroutine

	! copied from ICCP coding-notes
	subroutine init_random_seed()
		implicit none
		integer, allocatable :: seed(:)
		integer :: i, n, un, istat, dt(8), pid, t(2), s
		integer(8) :: count, tms

		call random_seed(size = n)
		allocate(seed(n))
		open(newunit=un, file="/dev/urandom", access="stream",&
		form="unformatted", action="read", status="old", &
		iostat=istat)
		if (istat == 0) then
			read(un) seed
			close(un)
		else
			call system_clock(count)
			if (count /= 0) then
				t = transfer(count, t)
			else
				call date_and_time(values=dt)
				tms = (dt(1) - 1970)*365_8 * 24 * 60 * 60 * 1000 &
				+ dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
				+ dt(3) * 24 * 60 * 60 * 60 * 1000 &
				+ dt(5) * 60 * 60 * 1000 &
				+ dt(6) * 60 * 1000 + dt(7) * 1000 &
				+ dt(8)
				t = transfer(tms, t)
			end if
			s = ieor(t(1), t(2))
			pid = getpid() + 1099279 ! Add a prime
			s = ieor(s, pid)
			if (n >= 3) then
				seed(1) = t(1) + 36269
				seed(2) = t(2) + 72551
				seed(3) = pid
				if (n > 3) then
					seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
				end if
			else
				seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
			end if
		end if
		call random_seed(put=seed)
	end subroutine init_random_seed

end module 