subroutine cubic_fcc_lattice(N_cell_dim, L_side, pos)	
	....

	n = 0
	do i = 1,N_cell_dim
		do j = 1,N_cell_dim
			do k = 1,N_cell_dim
				do l = 1,4
					n = n + 1
					pos(:,n) = L_side/N_cell_dim*((/ i-1, j-1, k-1 /) + fcc_cell(l))
				end do
			end do
		end do
	end do
end subroutine

