!argon gas in a box simulation, molecular dynamics
program argon_box

    implicit none

	!the cubic geometry sides of length = L
	!initial positions initialized according to fcc lattice structure	
	!number of fcc cells per cartesian dimension = Ncell, 
	!number of particles is N, (4 particles per cube)		
	integer, parameter :: N_cell_dim = 6, N_cell = N_cell_dim**3, N_part = N_cell*4
	real(8), parameter :: L_side = 1d0, T = 300, m = 39.948*1.660538921d-27
	!constants
	real(8), parameter :: PI = 4 * atan(1d0), Kb = 1.3806488d-23
	
	integer :: i,j,k,l,n !iteration variables
	real(8) :: xs(2) !random number variable

	! face centered cubic unit cell with basis particle positions
	real(8), dimension(1:3), parameter :: &
		fcc_part1 = (/0d0,  0d0,  0d0/),  &
		fcc_part2 = (/0d0,  5d-1, 5d-1/), &
		fcc_part3 = (/5d-1, 0d0,  5d-1/), &
		fcc_part4 = (/5d-1, 5d-1, 0d0/)	
	
	real(8), dimension(1:3,1:4), parameter :: &
	R_cell = reshape( (/fcc_part1, fcc_part2, fcc_part3, fcc_part4/), (/3,4/))
	
	! particle system arrays
	real(8), dimension(1:3, 1:N_part) :: POS, VEL
	
	
	
	print *, "face centered cubic unit cell"
	do i = 1,3
		print *, (R_cell(i,j), j=1,4)
	end do
	
	!initial positions of all particles according to an fcc lattice structure
	print *, "initializing initial particle positions"
	
	n = 0
	do i = 1,N_cell_dim
		do j = 1,N_cell_dim
			do k = 1,N_cell_dim
				do l = 1,4
					n = n + 1
					POS(:,n) = L_side/N_cell_dim*( (/ i-1, j-1, k-1 /) + R_cell(:,l) )
					print *, "particle:", n, "/",  N_part, "position:", POS(:,n)
				end do
			end do
		end do
	end do
	
	!initial particles velocities according to the maxwell distribution
	print *, "initializing initial particle velocities"
	call init_random_seed
	do n = 1,N_part
		do i = 1,3
			! each velocity component is normally distributed:
			! f(v) = sqrt(m/(2*PI*Kb*T))*exp(-m/(2*Kb*T) * (v**2))
			! -> box-muller transform from uniform to normal distribution
			CALL RANDOM_NUMBER(xs(1))
			CALL RANDOM_NUMBER(xs(2))
			! how to incorporate the standard deviation?
			VEL(i,n) = sqrt(m/(2*PI*Kb*T)) * box_muller(xs)
		end do
		print *,"particle:", n, "/",  N_part, "velocity:", VEL(:,n)
	end do
	
	!calculate time evolution
	
!	do n = 1,N_part		
!		U
!		F = -1d0*(/ , , /)
!		dv = F/m*dt		
!		v = v + dv
!		x = x + v*dt
!	end do
	

contains
	
	! adapted from ICCP coding-notes
	function box_muller(xs) result(zs)
	  implicit none
	  real(8) :: pi = 4*atan(1d0), zs !(2)
	  real(8), intent(in) :: xs(2)

	  zs = sqrt(-2d0*log(xs(1)))*cos(2*pi*xs(2))
	  !zs(2) = sqrt(-2d0*log(xs(1)))*sin(2*pi*xs(2))
	end function box_muller
	
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
	

end program

