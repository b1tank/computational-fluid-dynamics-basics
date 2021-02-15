program project1
!------------------------------------------------------------------------------------------------
! Fortran program to solve Laplace's equation on a plate of unit width and height 
! subject to the boundary conditions
!					T(0,y) = 0,
!					T(x,0) = 0,
!					T(1,y) = sin(pi*y),
!					T(x,1) = sin(pi*x)
!	The main program includes:
!	(a) array storage for the points of the grid;
!	(b) boundary conditions at the boundary points;
!	(c) some initial conditions at the interior points;
!	(d) T at each interior point by the average value of T for the nearest four points;
!	(e) the L(infty) norm of the difference between T(n+1) and T(n) to see if the scheme
!		is converged (Use 1E-5 as your primary convergence criteria, for some cases you
!		may need to consider going smaller);
!	(f) to step d. if convergence is not satisfied, else;
!	(g) write the solution to an output file;
!	(h) write the convergence history to another output file.

implicit none

! define variables
integer :: alstat_1, alstat_2, alstat_3, hori_grid, vert_grid, n1, n2, h, n, i, j, iter, 		&
			wropst_1, wropst_2, wropst_3, wropst_4, wropst_5
real, allocatable, dimension(:,:) :: T, errors, trun_err
real, parameter :: PI = 4. * atan(1.), criteria = 1E-5
real :: dx, dy, conver_norm, r, temp, L_1, L_2, L_infty
character(30) :: wrfile_convergence, wrfile_solution, wrfile_errors, wrfile_trunerrors,			&
					wrfile_norms


	write (*, '(/a/a/)') "------------------------------------------------------------",		&
				"Fortran solver to Laplace equation with boundary conditions."
	do
		write (*, '(2x, a)', advance = "no") "Enter grid sizes (horizontal and vertical): "
		
		! read grid size (horizontal and vertical)
		read (*,*) hori_grid, vert_grid
		
		! constrain the grid size between 1 and 1000000
		if (hori_grid >= 1 .and. hori_grid <= 1000000) then
			if (vert_grid >= 1 .and. vert_grid <= 1000000) then 
				exit
			else
				write (*, '(2x, a)') "Error, grid size must be integers between 1 and 1000000."
				write (*, '(2x, a)') "Please re-enter."
			end if
		else
			write (*, '(2x, a)') "Error, grid size must be integers between 1 and 1000000."
			write (*, '(2x, a)') "Please re-enter."
		end if
	end do
	
! allocate two dimension array for storing the points of the grid
	allocate (T(0:hori_grid, 0:vert_grid), stat = alstat_1)
	if (alstat_1 /= 0) then
		write (*, '(2x, a, /a)') "Error, unable to allocate memory. ", "Program terminated."
		stop
	end if
	
! initial conditions (all zero) at the interior points
	T = 0.
	
! allocate two dimension array for storing the errors
	allocate (errors(1:(hori_grid-1), 1:(vert_grid-1)), stat = alstat_2)
	if (alstat_2 /= 0) then
		write (*, '(2x, a, /a)') "Error, unable to allocate memory. ", "Program terminated."
		stop
	end if
	
! initial conditions for the errors
	errors = 0.	
	
! allocate two dimension array for storing the truncation errors
	allocate (trun_err(1:(hori_grid-1), 1:(vert_grid-1)), stat = alstat_3)
	if (alstat_3 /= 0) then
		write (*, '(2x, a, /a)') "Error, unable to allocate memory. ", "Program terminated."
		stop
	end if
	
! initial conditions for the truncation errors
	trun_err = 0.	
	
! prompt for the convergence history output file name
	do	
		write (*,'(/2x, a)', advance = "no") "Convergence history output File Name: "
		
		! read output file name
		read (*,*) wrfile_convergence
		
		! open output file (read access)
		! if open unsuccessful, display error message
		! 		otherwise, end loop
		open(12, file = wrfile_convergence, status="new", 				&
					action="write", position="rewind", 					&
					iostat=wropst_1)
		if (wropst_1 == 0) exit	
		write (*,'(2x, a,a/,a)') "Unable to create ",					&
				"convergence output file.", "Please re-enter..."
	end do
	write (12,'(a, T16, a)') "Iteration", "Convergence_norm"

! grid steps
	dx = 1.0 / hori_grid
	dy = 1.0 / vert_grid
	
! Boundary conditions
	do i = 0, hori_grid
		T(i,vert_grid) = sin(pi*(dx*i))	
		
		! eliminate the round error when T = sin(pi) = 0.
		if (abs(T(i,vert_grid) - 0.) < 1E-6) T(i,vert_grid) = 0.
	end do
	do j = 0, vert_grid
		T(hori_grid,j) = sin(pi*(dy*j))
		
		! eliminate the round error when T = sin(pi) = 0.
		if (abs(T(hori_grid,j) - 0.) < 1E-6) T(hori_grid,j) = 0.
	end do
	
! initialization including convergence norm and iteration time
	n1 = hori_grid - 1; n2 = vert_grid - 1; conver_norm = 0; iter = 0
	do 
		do i = 1, n1
			do j = 1, n2
				! five-point scheme
				r = 0.25 * (T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1))
				
				! to find maximum |T(n+1)-T(n)|, which is L_infinity
				temp = abs(T(i,j) - r)
				if (temp > conver_norm) conver_norm = temp
				
				! update temperature of the interior point
				T(i,j) = r
			end do
		end do
		
		iter = iter + 1
		
		! write the convergence history into the output file
		write (12, '(i5, T16, f9.7)') iter, conver_norm
		
		! convergence criteria for convergence norm
		if (conver_norm < criteria) then
			! if converged, output the convergence history to the file
			write (*, '(2x, a, a, i4, a)') "The solution has converged.",			 		&
							"Iteration time is ", iter, "."
			write (*, '(2x, a/)') "Convergence history output successfully done!"
			exit
		else
			! if not converged yet, set the L_infty value to be zero
			conver_norm = 0
		end if
		
		! if it has taken too long to converge, the program will be terminated
		if (iter > 1E5) then
			write (*, '(2x, a/, 2x, a)') "It has taken so long to converge.", 				&
										"Program terminated."
			stop
		end if
	end do

! initialize L_1, L_2 and L_infinity norm
	L_1 = 0.; L_2 = 0.; L_infty = 0.
	
! evaluate actual errors, truncation errors and three norms
	do i = 1, n1
		do j = 1, n2		
			errors(i,j) = abs((sin(PI*dx*i) * sinh(PI*dy*j) + sin(PI*dy*j) * sinh(PI*dx*i))		&
							/ (sinh(PI)) - T(i,j))
			trun_err(i,j) = (sin(PI*dx*i) * sinh(PI*dy*j) + sin(PI*dy*j) * sinh(PI*dx*i))		&
							/ (sinh(PI))*PI**4/48*(dx**4+dy**4)
			L_1 = L_1 + errors(i,j)
			L_2 = L_2 + errors(i,j)**2
			if (L_infty < errors(i,j)) L_infty = errors(i,j)
		end do
	end do
	
! calculate L_1 and L_2 norm
	L_1 = L_1 / (n1 * n2)
	L_2 = sqrt(L_2 / (n1 * n2))
	
! prompt for the solution output file name
	do	
		write (*, '(2x, a)', advance = "no") "Solution output File Name: "
		
		! read solution output file name
		read (*,*) wrfile_solution
		
		! open output file (read access)
		! if open successful, exit loop
		! 		otherwise, display error message
		open(14, file = wrfile_solution, status="new", 					&
					action="write", position="rewind", 					&
					iostat=wropst_2)
		if (wropst_2 == 0) exit		
		write (*, '(2x, a,a/,a)') "Unable to create ",					&
								"output file.", "Please re-enter..."								
	end do
	
! write the solution into the output file
	do j = 0, vert_grid
		do i = 0, hori_grid
			write (14, '(2x, f9.7)', advance = "no") T(i,j)
		end do
		write (14,*)
	end do	
	
	write (*, '(2x, a/)') "Solution output successfully done!"

	
! prompt for the errors output file name
	do	
		write (*, '(2x, a)', advance = "no") "Errors output File Name: "
		
		! read errors output file name
		read (*,*) wrfile_errors
		
		! open output file (read access)
		! if open successful, exit loop
		! 		otherwise, display error message
		open(16, file = wrfile_errors, status="new", 					&
					action="write", position="rewind", 					&
					iostat=wropst_3)
		if (wropst_3 == 0) exit		
		write (*, '(2x, a,a/,a)') "Unable to create ",					&
								"output file.", "Please re-enter..."								
	end do
	
! write the errors into the output file
	do j = 1, vert_grid-1
		do i = 1, hori_grid-1
			write (16, '(2x, f9.7)', advance = "no") errors(i,j)
		end do
		write (16,*)
	end do	
	
	write (*, '(2x, a/)') "Errors output successfully done!"	

! prompt for the truncation errors output file name
	do	
		write (*, '(2x, a)', advance = "no") "Truncation errors output File Name: "
		
		! read truncation errors output file name
		read (*,*) wrfile_trunerrors
		
		! open output file (read access)
		! if open successful, exit loop
		! 		otherwise, display error message
		open(18, file = wrfile_trunerrors, status="new", 					&
					action="write", position="rewind", 					&
					iostat=wropst_4)
		if (wropst_4 == 0) exit		
		write (*, '(2x, a,a/,a)') "Unable to create ",					&
								"output file.", "Please re-enter..."								
	end do
	
! write the truncation errors into the output file
	do j = 1, vert_grid-1
		do i = 1, hori_grid-1
			write (18, '(2x, f9.7)', advance = "no") trun_err(i,j)
		end do
		write (18,*)
	end do	
	
	write (*, '(2x, a/)') "Truncation errors output successfully done!"	


! prompt for the norms output file name
	do	
		write (*, '(2x, a)', advance = "no") "Norms output File Name: "
		
		! read norms output file name
		read (*,*) wrfile_norms
		
		! open output file (read access)
		! if open successful, exit loop
		! 		otherwise, display error message
		open(20, file = wrfile_norms, status="new", 					&
					action="write", position="rewind", 					&
					iostat=wropst_5)
		if (wropst_5 == 0) exit		
		write (*, '(2x, a,a/,a)') "Unable to create ",					&
								"output file.", "Please re-enter..."								
	end do
	
! write the norms into the output file
	write (20,'(a, T12, a, T24, a)') "L_1", "L_2", "L_infinity"
	write (20, '(f9.7, T12, f9.7, T24, f9.7)', advance = "no") L_1, L_2, L_infty
	
	write (*, '(2x, a/)') "Norms output successfully done!"	
	
	write (*, '(a)') "------------------------------------------------------------"
	
end program project1

