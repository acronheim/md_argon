program argon_box

    implicit none

	! N = # fcc cells in the cubic geometry, 4 particles per fcc cube
	! b1, b2 are the boundaries of the cubic geometry with sides of of equal length, boundary in origin and at L in every spatial dimension
	integer, parameter :: N = 16
	real(8), parameter, dimension(3) :: bc = (/0,0,0/), b2 = (/1,1,1/)
	
	! Variables
	integer :: i, j, Ncells
	real(8), dimension(N,3) :: r(N,3), v(N,3), dr(N,3), dv(N,3)
	
	Ncells = N/4
	Ncells_per_dim = Ncelss**(1/3)
	
	
	

	
end program

