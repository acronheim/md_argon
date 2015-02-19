<<<<<<< Local Changes
!argon gas in a box simulation, molecular dynamics.
!compile with: gfortran argon-box.f90 $(pkg-config --cflags --libs plplotd-f95)

!the cubic geometry sides of length = L
!initial positions initialized according to fcc lattice structure	
!number of fcc cells per cartesian dimension = Ncell, 
!number of particles is N, (4 particles per cube)	

!linear aproximation: x = x + v*dt, dv= F/m*dt: semi iplicit euler method
!time evolution for particles in lennard jones potential: U = 4*e*((s/r)**12-(s/r)**6),
!Fij = -du/dx = -du/dr*dr/dx = e*(48*s**12/r**13 - 6*s**6/r**7) * x/r, 
!r = sqrt(x**2+y**2+z**2)		

program argon_box
	use plplot
    implicit none
	
	integer, parameter :: N_cell_dim = 6
	real(8), parameter :: dt = 0.004_8, T_initial = 9d-1, rho = 0.85_8
	
	integer, parameter :: N_cell = N_cell_dim**3, N_part = N_cell*4
	real(8), parameter :: L_side = (N_part/rho)**(1._8/3), m = 1d0
	real(8), parameter :: s = 1d0, e = 1d0, r_cut = L_side ! lennard jones potential
	real(8), parameter :: t_stop = 1d0
	real(8), parameter :: Kb = 1d0 	!constants	
	!integer, parameter :: N_avSteps = 100 ! #steps used for ensemble average
	integer :: i,j,k,l,n, step !iteration variables	
	real(8) :: time, U, H, P, Temp
	real(8), dimension(1:3, 1:N_part) :: pos, vel 	! particle system arrays		
	
	call cubic_fcc_lattice(N_cell_dim, L_side, pos)
	call init_random_seed
	call init_vel(T_initial, Kb, m, vel)

	print *, "calculating time evolution"
	time = 0d0
	step = 0
	call plot_init(0d0, L_side,0d0, L_side,0d0, L_side) 
	
	do while (time < t_stop)
		call plot_points(pos)
		time = time + dt	
		step = step + 1	
		call calc_dynamics(N_part, L_side, dt, Kb, m, e, s, r_cut, pos, Temp, P, H, vel, step)
		call new_pos(N_part, L_side, vel, pos)
		print *, step,  "T =", Temp, "P =", P,  "vel1 =", vel(1,1), "pos1 =", pos(1,1)	
		!call rescale_vel(T_initial, Temp, dt, vel)
	end do	
	
	call plot_end

contains
	
	subroutine calc_dynamics(N_part, L_side, time_step, Kb, m, e, s, r_cut, pos, Temperature, Pressure, tot_energy, vel, cnt)
		! Force calculation	
		integer, intent(in) :: N_part, cnt
		real(8), intent(in) :: e, s, r_cut !Lennard Jones
		real(8), intent(in) :: m, time_step, L_side, Kb
		real(8), intent(inout), dimension(1:3, 1:N_part) :: vel
		real(8), intent(in), dimension(1:3, 1:N_part) :: pos
		integer :: i,j,k,l,n		
		real(8), intent(out) :: tot_energy, Temperature, Pressure
		real(8) :: v_2(N_part), F(3), dF(3), r, r_vec(3), virial, kin_energy, pot_energy
		virial = 0
		pot_energy = 0 !potential		
		do n = 1,N_part 	
			F = 0
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
							pot_energy = pot_energy + 5d-1*4*e*((s/r)**12-(s/r)**6)
							virial =  virial + dot_product(r_vec, dF)	 
						end if
					end if		
				end do 
				end do 
				end do
			end do			
			vel(:,n) = vel(:,n) + F/m*time_step	
			v_2(n) = dot_product(vel(:,n),vel(:,n))
			!print *,"particle:", n, "/",  N_part, "velocity:", (VEL(i,n), i=1,3)
		end do		
		kin_energy = sum(m/2*V_2)
		tot_energy = pot_energy + kin_energy
		Temperature = 2*kin_energy/(3*N_part*Kb)
		call write_energy_file(tot_energy, kin_energy, pot_energy, Temperature, cnt)
		Pressure = (1 + 1/(6*Kb*Temperature*N_part)* virial) !P/(Kb T rho) + correction cuttoff,.		
	end subroutine

	subroutine rescale_vel(T_intended, T_actual, time_step, Vel)
		!rescale velocities in order to keep temperature constant
		real(8), intent(in) :: T_intended, T_actual, time_step
		real(8), intent(inout), dimension(1:3, 1:N_part) :: vel
		real(8), parameter :: tau = 1d-3
		real(8) :: scaling_factor
		! Berendsen thermostat
		scaling_factor = sqrt(1 + time_step/tau*(T_intended/T_actual -1))
		vel = scaling_factor*vel
	end subroutine
	
	subroutine new_pos(N_part, L_side, vel, pos)
		! Postion calculation
		real(8), intent(in) :: L_side
		integer, intent(in) :: N_part
		real(8), intent(in), dimension(1:3, 1:N_part) :: vel
		real(8), intent(inout), dimension(1:3, 1:N_part) :: pos
		integer :: n
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
		real(8), intent(out), dimension(1:3, 1:N_part) :: pos
		print *, "initializing initial particle positions"
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
	
	subroutine init_vel(T, Kb, m, vel)
		!initial particles velocities according to the maxwell distribution	
		! in the maxwell boltzman distribution each velocity component is normally distributed:
		! Box muller transform used for converting uniform dist to normal dist
		!
		! f(v) = sqrt(m/(2*PI*Kb*T)) * exp(-(v**2)/2 *m/(*Kb*T))
		! sigma**2 = Kb*T/m, and zero mean	
		real(8), intent(in)  :: T, Kb, m
		real(8), intent(out), dimension(1:3, 1:N_part) :: vel
		real(8), parameter :: pi = 4*atan(1d0)
		real(8):: xs(2) !two random numbers
		print *, "initializing initial particle velocities"
		do n = 1,N_part
			do i = 1,3
				CALL RANDOM_NUMBER(xs(1))
				CALL RANDOM_NUMBER(xs(2))						
				vel(i,n) = sqrt(Kb*T/m) * sqrt(-2d0*log(xs(1)))*cos(2*pi*xs(2)) !sigma * box_muller
			end do
	!		print *,"particle:", n, "/",  N_part, "velocity:", VEL(:,n)
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
	
	!adapted from coding notes
	subroutine plot_init(xmin,xmax, ymin,ymax,zmin,zmax)
		implicit none
		real(8), intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax
		
	    call plsdev("xcairo")
		call plinit()
	    !call plparseopts(PL_PARSE_FULL)
		
		call pladv(0)
		call plvpor(0d0, 1d0, 0d0, 1d0)
		call plwind(-1d0, 1d0, -2d0 / 3, 4d0 / 3)
		call plw3d(1d0, 1d0, 1d0, xmin, xmax, ymin, ymax, &
		             zmin, zmax, 45d0, -45d0)
			
	end subroutine plot_init
	
	subroutine plot_end
		call plspause(.false.)
		call plend()
	end subroutine 
	
	!copied from coding notes
    subroutine plot_points(xyz)
	implicit none
		
      real(8), intent(in) :: xyz(:, :)

      call plclear()
      call plcol0(1)
      call plbox3("bnstu", "x", 0d0, 0, "bnstu", "y", &
                    0d0, 0, "bcnmstuv", "z", 0d0, 0)
      call plcol0(2)
      call plpoin3(xyz(1, :), xyz(2, :), xyz(3, :), 4)
      call plflush()
    end subroutine 

	subroutine write_energy_file(H, kin_energy, pot_energy, T, cnt)

		real(8), intent(in) :: H, kin_energy, T, pot_energy
		integer, intent(in) :: cnt

			

		open (unit=1,file="energy_matrix.dat",action="write")
  	
		write (1,"(I6, 4F18.6)")  cnt, H, kin_energy, pot_energy, T
  		

	end subroutine


end program
