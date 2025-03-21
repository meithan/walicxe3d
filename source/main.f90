!===============================================================================
!> @file main.f90
!> @brief Walixce3D main program unit
!> @author Juan C. Toledo
!> @date 2/Jun/2011

! Copyright (c) 2014 Juan C. Toledo and Alejandro Esquivel
!
! This file is part of Walicxe3D.
!
! Walicxe3D is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.

!===============================================================================

!> @brief Walicxe3D Main Program

!> This is the main program unit of the Walicxe3D code. It manages all the
!! high-level operations of the code, which include initializations, initial
!! flow conditions, adaptive mesh refinements, flow variables integration,
!! boundary conditions calculation, parallel load balancing.
program Walicxe3D

  use parameters   ! Code parameters
  use globals      ! Runtime global variables
  use tictoc       ! Time-measuring library
  implicit none    ! ALWAYS mandatory

  ! Timing mark
  call tic(start_mark)

  ! Initialize global arrays and variables
  call initmain ()

  ! Initialize base grid (cold start)
  call basegrid ()

  ! Impose initial conditions (cold start) or load them from file (warm start)
  if (dowarm) then
    call warmstart ()
  else
    call initflow ()
  end if

  ! Load balance initial state of simulation among all processors
  call loadBalance ()

  ! Update primitives
  call updatePrims ()

  ! Write initial condition to disk
  if (nextout.eq.0) then
    call output (0)
    nextout = 1
  end if

  ! Main Loop
!  do while(.false.)
  do while (time.le.tfin/t_sc)

    it = it + 1
    call tic(it_mark)

    write(logu,'(a)') "================================================================================"
    write(logu,'(1x,a,i0)') "Starting Iteration " , it
    write(logu,'(1x,a)') stamp()
    write(logu,'(a)') "================================================================================"

    ! Hydro Solver
    call hydroSolver ()

    ! Update primitives in all blocks
    call updatePrims ()

    ! Radiative cooling
    call cooling ()

    ! Update AMR grid
    call admesh ()

    ! Load balance
    call loadBalance ()

    ! Data output (if scheduled)
    if (time.ge.nextout*dtout/t_sc) then
      call output(nextout)
      nextout = nextout + 1
    end if

    ! Report progress
    call main_report ()

    ! Everyone stop here
    call mpi_barrier (mpi_comm_world, ierr)

  end do

  ! Deallocate globals and terminate execution
  write(logu,*) ""
  write(logu,'(a)') "================================================================================"  
  write(logu,'(a)') STAMP()
  write(logu,'(a)') 'Execution complete!'
  write(logu,'(a,a)') 'Total elapsed time: ', nicetoc(start_mark)
  call deinit()

end program Walicxe3D
