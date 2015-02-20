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
	use argon_box_init 
	use argon_box_dynamics
	use md_plot
	implicit none
	
	integer, parameter :: N_cell_dim = 4
	real(8), parameter :: dt = 0.004_8, T_initial = 9d-1, rho = 0.85_8, t_stop = 1d0
	
	integer, parameter :: N_cell = N_cell_dim**3, N_part = N_cell*4
	real(8), parameter :: L_side = (N_part/rho)**(1._8/3), m = 1d0
	
	real(8), parameter :: s = 1d0, e = 1d0, r_cut = L_side ! lennard jones potential
	real(8), parameter :: Kb = 1d0 	!Boltzman constant
	
	!integer, parameter :: N_avSteps = 100 ! #steps used for ensemble average
	integer :: i,j,k,l,n, step !iteration variables	
	real(8), dimension(1:3, 1:N_part) :: pos, vel 	
	real(8) :: time, kin_energy, pot_energy, virial
	real(8) :: Pressure, Temperature, tot_energy
	
	
	
	call cubic_fcc_lattice(N_cell_dim, L_side, pos)
	call init_random_seed
	call init_vel(T_initial, Kb, m, N_part, vel)

	time = 0d0
	step = 0
	call plot_init(0d0, L_side,0d0, L_side,0d0, L_side) 	
	do while (time < t_stop)
		time = time + dt	
		step = step + 1	
		call calc_dynamics(N_part, L_side, dt, Kb, m, e, s, r_cut, pos, kin_energy, pot_energy, virial, vel)																																		
		call new_pos(N_part, L_side, dt, vel, pos)
		!call rescale_vel(T_initial, Temperature, dt, N_part, vel)			
		
		tot_energy = pot_energy + kin_energy
		Temperature = 2*kin_energy/(3*N_part*Kb)		
		Pressure = (1 + 1/(6*Kb*Temperature*N_part)* virial) !P/(Kb T rho) + correction cuttoff	
		
		call plot_points(pos)	
		call write_energy_file(tot_energy, kin_energy, pot_energy, Temperature, step)
		print *, step,  "t=", time, "H=", tot_energy, "K=", kin_energy, "U=", pot_energy, "T =", Temperature, "P =", Pressure
	end do	
	
	call plot_end

contains

	subroutine write_energy_file(H, kin_energy, pot_energy, T, cnt)
		real(8), intent(in) :: H, kin_energy, T, pot_energy
		integer, intent(in) :: cnt
		open (unit=1,file="energy_matrix.dat",action="write")
		write (1,"(I6, 4F18.6)")  cnt, H, kin_energy, pot_energy, T
	end subroutine

end program
