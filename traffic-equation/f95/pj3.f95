program project3
!------------------------------------------------------------------------------------------------
! Fortran program to solve the traffic problem using the original closure equation.
! v = v_max * (1 - rho**2/rho_max**2)
! Various flux functions:
! (1) Exact Riemann Solution
!                 F(rho_L,rho_R) = rho_max * v_max / 4  (a_L <= 0 <= a_R)
!                                = f(rho_L)             (a_L + a_R >= 0)
!                                = f(rho_R)             else
! (2) Lax-Wendroff Flux Formula
!                 F(rho_L,rho_R) = f(1/2*(rho_L+rho_R) - Delta_t/(2*Delta_x)*(f(rho_R)-f(rho_L)))
! (3) Roe's Approximate Riemann Solver
!                 F(rho_L,rho_R) = f(rho_L)             (a_L + a_R >= 0)
!                                = f(rho_R)             else
! (4) Roe's Approximate Riemann Solver with Entropy Fix
!                 F(rho_L,rho_R) = 1/2*(f(rho_L) + f(rho_R) - 1/2*(a_R-a_L)*(rho_R-rho_L) sonic case
!                                = f(rho_L)                                               (a_L + a_R >= 0)
!                                = f(rho_R)                                               else

implicit none

! define variables
integer :: nn, i
integer, parameter :: M = 151, N = 251 ! large enough mesh sizes
real :: rho_max, rho_L, rho_R, v_max, Dx, Dt, SafeFactor, temp, a_max
real, allocatable, dimension(:,:) :: rho1, rho2, rho3, rho4 
real, dimension(1:M) :: F1, F2, F3, F4, a, x
real, dimension(1:N) :: t

	write (*, '(/a/a/)') "------------------------------------------------------------",		&
				"Fortran solver to traffic equation using the original closure equation."
  write (*, '(a)') "Enter safe factor:"
  read (*,*) SafeFactor
  
  Dx = 50. ! space step Delta_x
  rho_max = 0.1 ! maximum density
  rho_L = rho_max ! initial left density
  rho_R = 0. ! initial right density
  v_max = 50. ! maximum velocity

  ! create arrays storing data for four schemes
  call create_array(M,N,rho1)
  call create_array(M,N,rho2)
  call create_array(M,N,rho3)
  call create_array(M,N,rho4)

  ! specify the range [-1000,6500] of space coordinate 
  do i = 1, M
    x(i) = 0.0 - 1000.0 + (real(i)-1.0)*Dx
  end do

  ! 1st scheme - Exact Riemann solver
  nn = 1
  t(1) = 0. ! initial time 
  do
    a_max = 0.
    do i = 1, M
      a(i) = v_max * (1.0 - 2.0*rho1(i,nn)/rho_max)
      temp = abs(a(i))
      if (temp > a_max) then
        a_max = temp
      end if
    end do
    Dt = SafeFactor * Dx/a_max ! the timestep has to allow the numerical scheme to catch all the characteristics

    if (t(nn)+Dt <= 60.) then
      F1(1) = 0. ! the very left flux is zero
      do i = 2, M-1
        F1(i) = F_ER(rho1(i,nn), rho1(i+1,nn)) ! calculate all fluxes
      end do
      F1(M) = f(rho1(M,nn)) ! determine the rightmost flux using its left neighbor

      ! calculate solution in next time step
      do i = 2, M
        rho1(i,nn+1) = rho1(i,nn) - Dt/Dx * (F1(i) - F1(i-1))
      end do
      
      t(nn+1) = t(nn) + Dt ! increment the total simulation time
      nn = nn + 1 ! increment the time step number
    else 
      exit ! if the total time exceeds 60 seconds, terminate the loop
    end if
  end do

  write (*,'(a)') "Output for Exact Riemann:"
  call output(rho1,M,N,x,t,17) ! store the solution in a single file


  ! 2nd scheme - Lax-Wendroff scheme
  nn = 1
  t(1) = 0. ! initial time 
  do
    a_max = 0.
    do i = 1, M
      a(i) = v_max * (1.0 - 2.0*rho1(i,nn)/rho_max)
      temp = abs(a(i))
      if (temp > a_max) then
        a_max = temp
      end if
    end do
    Dt = SafeFactor * Dx/a_max ! the timestep has to allow the numerical scheme to catch all the characteristics

    if (t(nn)+Dt <= 60.) then
      F2(1) = 0. ! the very left flux is zero
      do i = 2, M-1
        F2(i) = F_LW(rho2(i,nn),rho2(i+1,nn),Dt,Dx) ! calculate all fluxes
      end do
      F2(M) = f(rho2(M,nn)) ! determine the rightmost flux using its left neighbor

      ! calculate solution in next time step
      do i = 2, M
        rho2(i,nn+1) = rho2(i,nn) - Dt/Dx * (F2(i) - F2(i-1))
      end do
      
      t(nn+1) = t(nn) + Dt ! increment the total simulation time
      nn = nn + 1 ! increment the time step number
    else 
      exit ! if the total time exceeds 60 seconds, terminate the loop
    end if
  end do

  write (*,'(a)') "Output for Lax-Wendroff:"
  call output(rho2,M,N,x,t,21) ! store the solution in a single file

  ! 3rd scheme - Roe's Approximate Riemann Solver
  nn = 1
  t(1) = 0. ! initial time
  do
    a_max = 0.
    do i = 1, M
      a(i) = v_max * (1.0 - 2.0*rho1(i,nn)/rho_max)
      temp = abs(a(i))
      if (temp > a_max) then
        a_max = temp
      end if
    end do
    Dt = SafeFactor * Dx/a_max ! the timestep has to allow the numerical scheme to catch all the characteristics

    if (t(nn)+Dt <= 60.) then
      F3(1) = 0. ! the very left flux is zero
      do i = 2, M-1
        F3(i) = F_RA(rho3(i,nn), rho3(i+1,nn)) ! calculate all fluxes
      end do
      F3(M) = f(rho3(M,nn)) ! determine the rightmost flux using its left neighbor

      ! calculate solution in next time step
      do i = 2, M
        rho3(i,nn+1) = rho3(i,nn) - Dt/Dx * (F3(i) - F3(i-1))
      end do
      
      t(nn+1) = t(nn) + Dt ! increment the total simulation time
      nn = nn + 1 ! increment the time step number
    else 
      exit ! if the total time exceeds 60 seconds, terminate the loop
    end if
  end do

  write (*,'(a)') "Output for Roe Approximate:"
  call output(rho3,M,N,x,t,25) ! store the solution in a single file

  ! 4th scheme - Roe's Approximate Riemann Solver with Entropy Fix
  nn = 1
  t(1) = 0. ! initial time
  do 
    a_max = 0.
    do i = 1, M
      a(i) = v_max * (1.0 - 2.0*rho1(i,nn)/rho_max)
      temp = abs(a(i))
      if (temp > a_max) then
        a_max = temp
      end if
    end do
    Dt = SafeFactor * Dx/a_max ! the timestep has to allow the numerical scheme to catch all the characteristics

    if (t(nn)+Dt <= 60.) then
      F4(1) = 0. ! the very left flux is zero
      do i = 2, M-1
        F4(i) = F_EF(rho4(i,nn), rho4(i+1,nn)) ! calculate all fluxes
      end do
      F4(M) = f(rho4(M,nn)) ! determine the rightmost flux using its left neighbor

      ! calculate solution in next time step
      do i = 2, M
        rho4(i,nn+1) = rho4(i,nn) - Dt/Dx * (F4(i) - F4(i-1))
      end do
      
      t(nn+1) = t(nn) + Dt ! increment the total simulation time
      nn = nn + 1 ! increment the time step number
    else 
      exit ! if the total time exceeds 60 seconds, terminate the loop
    end if
  end do

  write (*,'(a)') "Output for Entropy Fix:"
  call output(rho4,M,N,x,t,29) ! store the solution in a single file


  ! output the exact solution in a separate file
  call output_exact(M,N,x,t,31) 


  ! end of the main program
  ! *****************************************************************
  ! beginning of the subroutines and functions


contains
  ! create solution array
  subroutine create_array(M,N,rho)
    implicit none
    integer :: alstat_1, i
    integer, intent(in) :: M, N
    real, allocatable, dimension(:,:), intent(out) :: rho

    ! allocate two dimension array for storing the points of the grid
    allocate (rho(1:M, 1:N), stat = alstat_1)
    if (alstat_1 /= 0) then
      write (*, '(2x, a, /a)') "Error, unable to allocate memory. ", "Program terminated."
      stop
    end if
    ! initialize the array
    rho = 0.
    ! initial condition
    do i = 2, 21 
      rho(i,1) = 0.1 ! left boundary densities are all zero, while the densities in the range of [-1000,0] are all 0.1
    end do
  end subroutine create_array 
  
  ! interface flux function for four schemes
  function F_ER(rho_l,rho_r) ! Exact Riemann
    real :: rho_l, rho_r, a_L, a_R, F_ER 
    a_L = v_max * (1 - 2*rho_l/rho_max)
    a_R = v_max * (1 - 2*rho_r/rho_max)
    if (a_L <= 0 .and. 0 <= a_R) then
      F_ER = rho_max*v_max/4.0
    else if (a_L + a_R >= 0) then
      F_ER = f(rho_l)
    else
      F_ER = f(rho_r)
    end if
  end function F_ER

  function F_LW(rho_l,rho_r,Delta_t,Delta_x) ! Lax-Wendroff
    real :: rho_l, rho_r, Delta_t, Delta_x, F_LW, ff
    ff = (rho_l+rho_r)/2.0 - Delta_t/2.0/Delta_x*(f(rho_r)-f(rho_l))
    F_LW = f(ff)
  end function F_LW

  function F_RA(rho_l,rho_r) ! Roe's Approximate Riemann Solver
    real :: rho_l, rho_r, a_L, a_R, F_RA
    a_L = v_max * (1 - 2*rho_l/rho_max)
    a_R = v_max * (1 - 2*rho_r/rho_max)
    if (a_L + a_R >= 0) then
      F_RA = f(rho_l)
    else
      F_RA = f(rho_r)
    end if
  end function F_RA

  function F_EF(rho_l,rho_r) ! Roe's Approximate Riemann Solver Entropy Fix
    real :: rho_l, rho_r, a_L, a_R, F_EF
    a_L = v_max * (1 - 2*rho_l/rho_max)
    a_R = v_max * (1 - 2*rho_r/rho_max)
    if (a_L <= 0 .and. a_R >= 0) then ! sonic case
      F_EF = 1.0/2.0*(f(rho_l)+f(rho_r)-1.0/2.0*(a_R-a_L)*(rho_r-rho_l))
    else if (a_L + a_R >= 0) then
      F_EF = f(rho_l)
    else 
      F_EF = f(rho_r)
    end if
  end function F_EF

  function f(x)
    real :: f, x
    f = x * v_max * (1.0 - x/rho_max)
  end function f

  ! solution output subroutine
  subroutine output(rho,M,N,x,t,file_num)
    implicit none
    character(30) :: wrfile_output
    integer :: wropst_1, i, nn
    integer, intent(in) :: M, N, file_num
    real, dimension(1:M), intent(in) :: x
    real, dimension(1:N), intent(in) :: t
    real, allocatable, dimension(:,:), intent(in) :: rho

    ! prompt for the errors output file name
    do	
      write (*, '(2x, a)', advance = "no") "Output File Name: "
      
      ! read errors output file name
      read (*,*) wrfile_output
      
      ! open output file (read access)
      ! if open successful, exit loop
      ! 		otherwise, display error message
      open(file_num, file = wrfile_output, status="new", 					&
            action="write", position="rewind", 					&
            iostat=wropst_1)
      if (wropst_1 == 0) exit		
      write (*, '(2x, a, a/, a)') "Unable to create ",					&
                  "output file.", "Please re-enter..."								
    end do

    ! store the values of space coordinates in the first row of the solution file 
    write (file_num, '(2x, f10.5)', advance = "no") 0.0 ! merely a placeholder without physical meaning
    do i = 1, M
      write (file_num, '(2x, f10.4)', advance = "no") x(i)
    end do
    write (file_num,*)

    do nn = 1, N
      write (file_num, '(2x, f10.7)', advance = "no") t(nn) ! store the values of time in the first column of the solution file
      do i = 1, M
        write (file_num, '(2x, f10.7)', advance = "no") rho(i,nn)
      end do
      write (file_num,*)
    end do	
    
    write (*, '(2x, a)') "Solution output successfully done!"
  end subroutine output

  ! exact solution output subroutine
  subroutine output_exact(M,N,x,t,file_num)
    implicit none
    character(30) :: wrfile_output
    integer :: wropst_1, i, nn
    integer, intent(in) :: M, N, file_num
    real :: a_L, a_R
    real, dimension(1:M), intent(in) :: x
    real, dimension(1:N), intent(in) :: t

    ! prompt for the errors output file name
    do	
      write (*, '(2x, a)', advance = "no") "Output File Name: "
      
      ! read errors output file name
      read (*,*) wrfile_output
      
      ! open output file (read access)
      ! if open successful, exit loop
      ! 		otherwise, display error message
      open(file_num, file = wrfile_output, status="new", 					&
            action="write", position="rewind", 					&
            iostat=wropst_1)
      if (wropst_1 == 0) exit		
      write (*, '(2x, a, a/, a)') "Unable to create ",					&
                  "output file.", "Please re-enter..."								
    end do

    ! store the values of space coordinates in the first row of the solution file 
    write (file_num, '(2x, f10.5)', advance = "no") 0.0 ! merely a placeholder without physical meaning
    do i = 1, M
      write (file_num, '(2x, f10.4)', advance = "no") x(i)
    end do
    write (file_num,*)

    ! write the initial condition into the second row of the solution file
    do i = 1, M 
      if (i >= 3 .and. i <= 22 ) then
        write (file_num, '(2x, f10.5)', advance = "no") 0.1  
      else
        write (file_num, '(2x, f10.5)', advance = "no") 0.0
      end if
    end do
    write (file_num,*)

    a_L = v_max * (1 - 2*0.1/rho_max)
    a_R = v_max * (1 - 2*0.0/rho_max)
    do nn = 2, N
      write (file_num, '(2x, f10.7)', advance = "no") t(nn) ! store the values of time in the first column of the solution file
      do i = 1, M
        if (x(i)/t(nn) < a_L) then
          write (file_num, '(2x, f10.7)', advance = "no") 0.1
        else if (x(i)/t(nn) > a_R) then
          write (file_num, '(2x, f10.7)', advance = "no") 0.0
        else 
          write (file_num, '(2x, f10.7)', advance = "no") rho_max/2.0*(1.0-x(i)/t(nn)/v_max)
        end if
      end do
      write (file_num,*)
    end do	
    
    write (*, '(2x, a)') "Exact solution output successfully done!"
  end subroutine output_exact

end program project3

