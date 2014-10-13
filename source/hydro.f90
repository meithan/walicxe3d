!===============================================================================
!> @file hydro.f90
!> @brief Hydrodynamical solver
!> @author Juan C. Toledo
!> @date 7/Mar/2012

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

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.

!===============================================================================

!> @brief High-level wrapper for the numerical integrator
subroutine hydroSolver ()

  use parameters
  use globals
  use tictoc
  implicit none

  integer :: mark, mark1
  real :: dt_new

  call tic(mark)
  write(logu,*) ""
  write(logu,'(1x,a)') "============================================"
  write(logu,'(1x,a)') " Integrating hydrodynamical equations ..."
  write(logu,'(1x,a)') "============================================"
  write(logu,*) ""

  ! Compute global timestep
  call tic(mark1)
  write(logu,'(1x,a)') "> Calculating timestep ..."
  call getTimestep (dt_new)
  if ((dt.gt.0.0).and.(abs(dt-dt_new)/dt.gt.50)) then
    write(logu,'(a)') "Abnormal timestep change!"
    write(logu,'(a,es12.5)') "Previous timestep:", dt
    write(logu,'(a,es12.5)') "New timestep:     ", dt_new
    write(logu,'(a)') "***ABORTING***"
    call clean_abort(ERROR_TIMESTEP)
  end if
  dt = dt_new
  write(logu,'(1x,a,es12.5,a)') "dt= ", dt, " code units"
  write(logu,'(1x,a,f12.5,a)') "dt= ", dt*t_sc/YR, " yr"
  write(logu,'(1x,a,a)') "Timestep calculated in ", nicetoc(mark1)

  ! Call appropriate numerical solver
  ! Once this routine returns, the UPs are expected to have correct
  ! timestepped values in all physical cells (ghost cells are unimportant)
  call tic(mark1)

  select case (solver_type)
  
    case (SOLVER_LAX)
      call LaxFriedrichs ()

    case (SOLVER_HLL1)
      call Godunov (1)

    case (SOLVER_HLL)
      call Godunov (2)

    case (SOLVER_HLLC)
      call Godunov (2)

  end select
  
  write(logu,*) ""
  write(logu,'(1x,a,a)') "Integrator (total)=", nicetoc(mark1)

  ! Step hydro variables
  call tic(mark1)
  call doStep()
  write(logu,*) "Stepping=", nicetoc(mark1)

  ! Done
  write(logu,*) ""
  write(logu,'(1x,a,a)') "> Numerical integration completed in ", nicetoc(mark)

end subroutine hydroSolver

!===============================================================================

!> @brief Obtains the simulation-wide numerical timestep
!> @details Will compute all the relevant timesteps (hydrodynamical, cooling,
!! etc) and return the largest, modified by the CFL parameter
!> @param dt_glob The global timestep
subroutine getTimestep (dt_glob)

  use parameters
  use globals
  implicit none

  real, intent(out) :: dt_glob

  integer :: nb, bID, i, j, k, ilev
  real :: cs, dt_hydro, dt_cool, dt_loc

  dt_loc = 1.0e30
  dt_glob = 1.0e30
  dt_hydro = 1.0e30
  dt_cool = 1.0e30

  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then

      call meshlevel(bID, ilev)

      do i=1,ncells_x
        do j=1,ncells_y
          do k=1,ncells_z

            ! Hydrodynamical timestep
            call sound (PRIM(nb,:,i,j,k),cs)
            dt_hydro = min( dt_hydro, &
              dx(ilev)/(abs(PRIM(nb,2,i,j,k))+cs), &
              dy(ilev)/(abs(PRIM(nb,3,i,j,k))+cs), &
              dz(ilev)/(abs(PRIM(nb,4,i,j,k))+cs) )

            ! Cooling timestep - not needed at the moment
      
          end do
        end do
      end do

    end if
  end do

  ! Limit the timestep (optional)
  if (cooling_type.eq.COOL_NONE) then
    dt_loc = dt_hydro*CFL
  else
    dt_loc = min( dt_hydro, dt_cool ) * CFL
  end if

#ifdef MPIP
  ! Obtain the global minimum timestep through MPI
  call MPI_ALLREDUCE(dt_loc, dt_glob, 1, mpi_real_kind, mpi_min, mpi_comm_world, ierr)
#else
  dt_glob = dt_loc
#endif

  return

end subroutine getTimestep

!===============================================================================

!> @brief Advances simulation by one timestep
subroutine doStep ()

  use parameters
  use globals
  implicit none

  integer :: nb, bID, i, j, k

  time = time + dt
 
  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then

      do i=1,ncells_x
        do j=1,ncells_y
          do k=1,ncells_z          
            U(nb,:,i,j,k) = UP(nb,:,i,j,k)
          end do
        end do
      end do

    end if
  end do
  
end subroutine doStep

!===============================================================================

!> @brief Applies artificial viscosity to a block's flow variables
!> @details Viscosity: up(n)=up(n)+eta*\nabla^2(u(n))
!> @param bIndx Block's local index
!> @param U Flow variables at t^n
!> @param UP Flow variables at t^(n+1)
subroutine viscosity (bIndx, U, UP)

  use parameters
  implicit none

  integer, intent(in) :: bIndx
  real, intent(in) :: U (nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)
  real, intent(inout) :: UP (nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)

  integer :: i, j, k

  do i=1,ncells_x
    do j=1,ncells_y
      do k=1,ncells_z
        UP(bIndx,:,i,j,k)   &
          = UP(bIndx,:,i,j,k) + visc_eta*(   &
          + U(bIndx,:,i+1,j,k) + U(bIndx,:,i-1,j,k)   &
          + U(bIndx,:,i,j+1,k) + U(bIndx,:,i,j-1,k)   &
          + U(bIndx,:,i,j,k+1) + U(bIndx,:,i,j,k-1)   &
          - 6.0*U(bIndx,:,i,j,k) )
      end do
    end do
  end do

end subroutine viscosity

!===============================================================================

!> @brief Calculates the soundspeed given a vector of primitives
!> @param pvars Vector of primitive variables
!> @param csound The hydrodynamical speed of sound
subroutine sound (pvars, csound)

  use parameters
  implicit none

  real, intent(in) :: pvars(neqtot)
  real, intent(out) :: csound

  csound = sqrt(gamma*pvars(5)/pvars(1))

  return

end subroutine sound

!===============================================================================
