module md_plot
	use plplot	
	implicit none
  	private
	
	public plot_init, plot_end, plot_points

contains
	
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
	
end module
