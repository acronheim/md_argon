!argon gas in a box simulation, molecular dynamics.

!the cubic geometry sides of length = L
!initial positions initialized according to fcc lattice structure	
!number of fcc cells per cartesian dimension = Ncell, 
!number of particles is N, (4 particles per cube)	

!velocity verlet method: v' = v+1/2*F(x)/m*dt; x = x+v'*dt; v = v'+1/2*F/m*dt. 
! (coverted from initially an implementation of the semi iplicit euler method)
!time evolution for particles in lennard jones potential: U = 4*e*((s/r)**12-(s/r)**6),
!Fij = -du/dx = -du/dr*dr/dx = e*(48*s**12/r**13 - 6*s**6/r**7) * x/r, 
!r = sqrt(x**2+y**2+z**2)		

program argon_box
	use argon_box_init 
	use argon_box_dynamics
	use argon_box_results
!	use md_plot
	implicit none
	
	integer, parameter :: N_cell_dim = 6, velocity_rescale_steps = 50
	real(8), parameter :: dt = 0.004_8, T_initial = 1d0, rho = 0.88_8, t_stop = 5d0
	
	integer, parameter :: N_cell = N_cell_dim**3, N_part = N_cell*4
	real(8), parameter :: L_side = (N_part/rho)**(1._8/3) 
	
	real(8), parameter :: s = 1d0, e = 1d0, r_cut = 5d-1*L_side ! lennard jones potential
	real(8), parameter :: m = 1d0, Kb = 1d0 	!mass and boltzman constant
	
	!integer, parameter :: N_avSteps = 100 ! #steps used for ensemble average
!	integer :: i,j,k,l,n, step !iteration variables	
	integer :: step 
	real(8), dimension(1:3, 1:N_part) :: pos, vel 	
	real(8) :: time, kin_energy, pot_energy, virial, prev_kin_energy, sum_deltaK_sqr, sum_kin_energy
	real(8) :: Pressure, Temperature, tot_energy

	! Create initial state
	call cubic_fcc_lattice(N_cell_dim, L_side, pos)
	call init_random_seed
	call init_vel(T_initial, Kb, m, N_part, vel)
	
!	call plot_init(0d0, L_side,0d0, L_side,0d0, L_side)

	time = 0d0
	step = 0
	prev_kin_energy = 0d0
	sum_deltaK_sqr = 0d0
	sum_kin_energy = 0d0

	do while (time < t_stop)
		time = time + dt	
		step = step + 1	
		
		!velocity verlet integration method, .true. triggers the calculation of thermodynamic quantities.
		call calc_dynamics(.false., N_part, L_side, dt, m, e, s, r_cut, pos, kin_energy, pot_energy, virial, vel) 
		call new_pos(N_part, L_side, dt, vel, pos)
		call calc_dynamics(.true., N_part, L_side, dt, m, e, s, r_cut, pos, kin_energy, pot_energy, virial, vel)
		
		!Temperature control
		if (step < velocity_rescale_steps) then
			call rescale_vel(T_initial, kin_energy, Kb, N_part, vel)			
		end if
		
		tot_energy = pot_energy + kin_energy
		Temperature = 2*kin_energy/(3* (N_part-1) *Kb)	!Center of mass degrees of freedom substracted..	
		Pressure = (1 + 1/(3*Kb*Temperature*N_part)* virial) !P/(Kb T rho) + correction cuttoff	
		
		!call plot_points(pos)	
		call write_energy_file(tot_energy, kin_energy, pot_energy, Temperature, step)
		call calc_specific_heat(.false., N_part, kin_energy, prev_kin_energy, sum_kin_energy, sum_deltaK_sqr)

		prev_kin_energy = kin_energy

		tot_energy = tot_energy/N_part
		pot_energy = pot_energy/N_part 
		kin_energy = kin_energy/N_part		
		print *, step,  "t=", time, "H=", tot_energy, "K=", kin_energy, "U=", pot_energy, "T =", Temperature, "P =", Pressure
		
	end do	
	
	call write_pos_correlation(L_side, 10, pos, N_part, step)
	call calc_specific_heat(.true., N_part, kin_energy, prev_kin_energy, sum_kin_energy, sum_deltaK_sqr)
!	call plot_end


end program
