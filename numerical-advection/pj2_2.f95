program project2_2
!------------------------------------------------------------------------------------------------
! Fortran program to solve the problem for the following
! special choices of q.
!                 q(v) = 1               Lax-Friedrichs
!                 q(v) = |v|             First-order upwind
!                 q(v) = v^2             Lax-Wendroff
!                 q(v) = 1/3 + 2/3*v^2   Minimum dispersion
! initial condition: u(x,0) = sin(2*PI*x)
! for 0 <= x <= 1 and periodic boundary conditions, u(0,t) = u(1,t).
! uniform grid spacings: 32
! compute to a time T where a*T = 3/2
! using v = 1.0, 0.75, 0.50, 0.25
!
! figure out amplitude and phase errors for each case.
! compute L_1 norm where N = 3/(2*v*Delta x), the final time step.

implicit none

! define variables
integer :: N, nn, i, j, j1, j2, j3, j4, vv
integer, parameter :: M = 32
real :: a, v, Delta_x, Delta_t, &
  pi1, pi2, pi3, pi4, p2_1, p2_2, p2_3, p2_4, amp1, amp2, amp3, amp4, &
  ph1, ph2, ph3, ph4, pherr1, pherr2, pherr3, pherr4, &
  L1_1, L1_2, L1_3, L1_4, s1, s2, s3, s4
real, parameter :: PI = 4. * atan(1.)
real, allocatable, dimension(:,:) :: u1, u2, u3, u4
real, dimension(1:M) :: F1, F2, F3, F4

	write (*, '(/a/a/)') "------------------------------------------------------------",		&
				"Fortran solver to advection equation with 32 grid uniform spacing."
	write (*, '(a)') "Please enter the wave speed:"
  read (*,*) a
  
  write (*,'(/a)') "Output for amplitudes:"
  call create_file(12)
  write (12, '(2x, a, T20, a, T40, a, T60, a, T80, a)') "v", "Lax-Friedrichs", "First-orderUpwind", &
    "Lax-Wendroff", "MinimumDispersion"
  
  write (*,'(a)') "Output for phase errors:"
  call create_file(14)
  write (14, '(2x, a, T20, a, T40, a, T60, a, T80, a)') "v", "Lax-Friedrichs", "First-orderUpwind", &
    "Lax-Wendroff", "MinimumDispersion"

  write (*,'(a)') "Output for L_1 norms:"
  call create_file(16)
  write (16, '(2x, a, T20, a, T40, a, T60, a, T80, a)') "v", "Lax-Friedrichs", "First-orderUpwind", &
    "Lax-Wendroff", "MinimumDispersion"

  v = 1.0

  do vv = 1, 4
    write (*, '(/a, f4.2)') "CFL number v = ", v
    ! confine the x in the range of [0,1]
    Delta_x = 1.0/M 
    Delta_t = v * Delta_x / a
    ! total time steps
    N = 3.0/2.0/(a*Delta_t)

    ! create arrays storing data for four schemes
    call create_array(M,u1,a,v,N,pi1)
    call create_array(M,u2,a,v,N,pi2)
    call create_array(M,u3,a,v,N,pi3)
    call create_array(M,u4,a,v,N,pi4)

    do nn = 1, N
      ! calculate all fluxes
      do i = 1, M-1
        F1(i) = LaxFriedrichs(u1(i+1,nn), u1(i,nn))
        F2(i) = FirstOrderUpwind(u2(i+1,nn), u2(i,nn))
        F3(i) = LaxWendroff(u3(i+1,nn), u3(i,nn))
        F4(i) = MinimumDispersion(u4(i+1,nn), u4(i,nn))
      end do
      F1(M) = LaxFriedrichs(u1(1,nn), u1(M,nn))
      F2(M) = FirstOrderUpwind(u2(1,nn), u2(M,nn))
      F3(M) = LaxWendroff(u3(1,nn), u3(M,nn))
      F4(M) = MinimumDispersion(u4(1,nn), u4(M,nn))

      ! calculate solution in next time step
      u1(1,nn+1) = u1(1,nn) - v/a * (F1(1) - F1(M))
      u2(1,nn+1) = u2(1,nn) - v/a * (F2(1) - F2(M))
      u3(1,nn+1) = u3(1,nn) - v/a * (F3(1) - F3(M))
      u4(1,nn+1) = u4(1,nn) - v/a * (F4(1) - F4(M))
      ! to obtain the peak value of the solution in 16th time step, which can be used for calculation of amplitude
      if (nn == 15) then
        p2_1 = u1(1,nn+1)
        p2_2 = u2(1,nn+1)
        p2_3 = u3(1,nn+1)
        p2_4 = u4(1,nn+1)
      end if
      ! calculate solution in next time step
      do i = 2, M
        u1(i,nn+1) = u1(i,nn) - v/a * (F1(i) - F1(i-1))
        u2(i,nn+1) = u2(i,nn) - v/a * (F2(i) - F2(i-1))
        u3(i,nn+1) = u3(i,nn) - v/a * (F3(i) - F3(i-1))
        u4(i,nn+1) = u4(i,nn) - v/a * (F4(i) - F4(i-1))
        if (nn == 15) then
          ! to obtain the peak value of the solution in 16th time step
          if (p2_1 < u1(i,nn+1)) then
            p2_1 = u1(i,nn+1)
          end if
          if (p2_2 < u2(i,nn+1)) then
            p2_2 = u2(i,nn+1)
          end if
          if (p2_3 < u3(i,nn+1)) then
            p2_3 = u3(i,nn+1)
          end if
          if (p2_4 < u4(i,nn+1)) then
            p2_4 = u4(i,nn+1)
          end if
          ! to obtain the phase of the point of inflection (u=0)
          if (u1(i,nn+1)<0 .and. u1(i-1,nn+1)>0) then
            ph1 = (u1(i-1,nn+1)/(u1(i-1,nn+1)-u1(i,nn+1)) + (i-2))*Delta_x
          end if
          if (u2(i,nn+1)<0 .and. u2(i-1,nn+1)>0) then
            ph2 = (u2(i-1,nn+1)/(u2(i-1,nn+1)-u2(i,nn+1)) + (i-2))*Delta_x
          end if
          if (u3(i,nn+1)<0 .and. u3(i-1,nn+1)>0) then
            ph3 = (u3(i-1,nn+1)/(u3(i-1,nn+1)-u3(i,nn+1)) + (i-2))*Delta_x
          end if
          if (u4(i,nn+1)<0 .and. u4(i-1,nn+1)>0) then
            ph4 = (u4(i-1,nn+1)/(u4(i-1,nn+1)-u4(i,nn+1)) + (i-2))*Delta_x
          end if
        end if
      end do
      ! calculate solution at the right boundary in next time step
      u1(M+1,nn+1) = u1(1,nn+1)
      u2(M+1,nn+1) = u2(1,nn+1)
      u3(M+1,nn+1) = u3(1,nn+1)
      u4(M+1,nn+1) = u4(1,nn+1)
      ! to obtain the peak value of the solution in 16th time step
      if (nn == 15) then
        if (p2_1 < u1(M+1,nn+1)) then
          p2_1 = u1(M+1,nn+1)
        end if
        if (p2_2 < u2(M+1,nn+1)) then
          p2_2 = u2(M+1,nn+1)
        end if
        if (p2_3 < u3(M+1,nn+1)) then
          p2_3 = u3(M+1,nn+1)
        end if
        if (p2_4 < u4(M+1,nn+1)) then
          p2_4 = u4(M+1,nn+1)
        end if
      end if
    end do

    ! calculate the relative phase error
    pherr1 = (ph1/(0.5 + a * 15.0 * Delta_t) - 1.0) / 15.0
    pherr2 = (ph2/(0.5 + a * 15.0 * Delta_t) - 1.0) / 15.0
    pherr3 = (ph3/(0.5 + a * 15.0 * Delta_t) - 1.0) / 15.0
    pherr4 = (ph4/(0.5 + a * 15.0 * Delta_t) - 1.0) / 15.0

    ! calculate the amplitude using two peak values in two time steps
    amp1 = (p2_1/pi1)**(1.0/15.0)
    amp2 = (p2_2/pi2)**(1.0/15.0)
    amp3 = (p2_3/pi3)**(1.0/15.0)
    amp4 = (p2_4/pi4)**(1.0/15.0)

    ! calculate the L_1 error norms
    s1 = 0.; s2 = 0.; s3 = 0.; s4 = 0.
    do j = 1, M+1
      s1 = s1 + abs(u1(j,N+1) - sin(2*PI*((j-1)*Delta_x-1.5)))
      s2 = s2 + abs(u2(j,N+1) - sin(2*PI*((j-1)*Delta_x-1.5)))
      s3 = s3 + abs(u3(j,N+1) - sin(2*PI*((j-1)*Delta_x-1.5)))
      s4 = s4 + abs(u4(j,N+1) - sin(2*PI*((j-1)*Delta_x-1.5)))
    end do
    L1_1 = s1 / (M+1)
    L1_2 = s2 / (M+1)
    L1_3 = s3 / (M+1)
    L1_4 = s4 / (M+1)


    write (12, '(2x, f4.2, T20, f8.5, T40, f8.5, T60, f8.5, T80, f8.5)') v, amp1, amp2, amp3, amp4
    write (14, '(2x, f4.2, T20, f11.7, T40, f11.7, T60, f11.7, T80, f11.7)') v, pherr1, pherr2, pherr3, pherr4
    write (16, '(2x, f4.2, T20, f8.4, T40, f8.4, T60, f8.4, T80, f8.4)') v, L1_1, L1_2, L1_3, L1_4
    write (*,'(a)') "Output for Lax-Friedrichs:"
    call output(u1,M,N,17+vv)
    write (*,'(a)') "Output for First-order Upwind:"
    call output(u2,M,N,21+vv)
    write (*,'(a)') "Output for Lax-Wendroff:"
    call output(u3,M,N,25+vv)
    write (*,'(a)') "Output for Minimum Dispersion:"
    call output(u4,M,N,29+vv)
    
    v = v - 0.25

  end do

  write (*, '(/2x, a)') "Amplitudes output successfully done!"
  write (*, '(2x, a)') "Phase errors output successfully done!"
  write (*, '(2x, a)') "L_1 norms output successfully done!"

contains
  ! create solution array
  subroutine create_array(M,u,a,v,N,peak)
    implicit none
    integer :: alstat_1, i
    integer, intent(in) :: M
    integer, intent(out) :: N
    real, intent(in) :: a, v
    real, intent(out) :: peak
    real, allocatable, dimension(:,:), intent(out) :: u

    ! allocate two dimension array for storing the points of the grid
    allocate (u(1:(M+1), 1:N+1), stat = alstat_1)
    if (alstat_1 /= 0) then
      write (*, '(2x, a, /a)') "Error, unable to allocate memory. ", "Program terminated."
      stop
    end if
    ! initialize the array
    u = 0.
    ! initialize the peak value on the boundary curve
    peak = 0.
    ! initial condition
    do i = 1, M
      u(i,1) = sin(2*PI*(i-1)*Delta_x)
      if (peak < u(i,1)) then
        peak = u(i,1)
      end if
    end do
    ! periodic boundary condition
    u(M+1,1) = u(1,1)
  end subroutine create_array 
  
  ! interface flux function for four schemes
  function LaxFriedrichs(u11,u12)
    real :: u11, u12, LaxFriedrichs
    LaxFriedrichs = 1.0/2.0*(f(u11)+f(u12)) - 1.0/2.0*a/v*(u11-u12)
  end function LaxFriedrichs

  function FirstOrderUpwind(u21,u22)
    real :: u21, u22, FirstOrderUpwind
    FirstOrderUpwind = 1.0/2.0*(f(u21)+f(u22)) - abs(v)/2.0*a/v*(u21-u22)
  end function FirstOrderUpwind

  function LaxWendroff(u31,u32)
    real :: u31, u32, LaxWendroff
    LaxWendroff = 1.0/2.0*(f(u31)+f(u32)) - v**2.0/2.0*a/v*(u31-u32)
  end function LaxWendroff

  function MinimumDispersion(u41,u42)
    real :: u41, u42, MinimumDispersion
    MinimumDispersion = 1.0/2.0*(f(u41)+f(u42)) - (1.0/3.0+2.0/3.0*v**2.0)/2.0*a/v*(u41-u42)
  end function MinimumDispersion

  function f(x)
    real :: f, x
    f = a * x
  end function f

  ! solution output subroutine
  subroutine output(u,M,N,file_num)
    implicit none
    character(30) :: wrfile_output
    integer :: wropst_1, i, nn
    integer, intent(in) :: M, N, file_num
    real, allocatable, dimension(:,:), intent(in) :: u

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

    do nn = 1, N+1
      do i = 1, M+1
        write (file_num, '(2x, f10.7)', advance = "no") u(i,nn)
      end do
      write (file_num,*)
    end do	
    
    write (*, '(2x, a)') "Solution output successfully done!"
  end subroutine output

  ! create output file for amplitude, phase error and L_1 norm
  subroutine create_file(file_num)
    implicit none
    character(30) :: wrfile_output
    integer :: wropst_1, file_num
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
      write (*, '(2x, a,a/,a)') "Unable to create ",					&
                  "output file.", "Please re-enter..."								
    end do
  end subroutine create_file

end program project2_2

