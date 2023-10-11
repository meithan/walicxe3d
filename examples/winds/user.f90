!===============================================================================
!> @file user.f90
!> @brief User-specified initial and boundary conditions
!> @author Juan C. Toledo
!> @date 20/May/2013

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

! EXAMPLE: winds
! This example shows how to impose wind sources, both plane and spherical.

! In this example we simulate the bow shock produced when a spherical wind
! source moves through the ISM. This is achieved by imposing a spherical wind
! source at the center of a 1x1x1 pc box as well as a planar wind incoming
! from one of the edges of the box.

! More specifically, we do three things:
! 1) The initial ISM is set in motion initially.
! 2) As part of the initial conditions, we impose both the spherical wind
!    source and the planar wind source (coming in from the top of the
!    computational box moving towards -Z).
! 3) At each boundary condition exchange, these wind sources are refreshed.

! The wind sources are imposed by calling subroutines from the winds module.
! These routines must be given a wind parameters object (of the correspodning
! type) which holds all the parameters of the wind.

module userconds
! ============================================
!< @brief User-specified initial and boundary conditions
!< @details This is the file where the user sets the initial conditions of
!! the simulation, as well as any special boundary conditions or custom
!! modifications to the simulation.
!!
!! How to use this module:
!!   1) Add "use" statements to include any modules you require in the
!!      section marked [1].
!!   2) You can define aditional module-wide global variables that you
!!      may require in the section marked [2].
!!   3) Fill in the subroutine userInitialCondition(), marked [3], which
!!      is called at the beginning of the simulation.
!!   4) Optionally, fill in the subroutine userBoundary(), marked [4],
!!      which is called at the end of each boundary exchange operation.
!! All subroutines in this module automatically have access to the global
!! parameters and variables.
! ============================================

  use parameters
  use globals
  ! ============================================
  ! [1] Add HERE any aditional modules required by your simulation

  use winds     ! Here we import the winds module

  ! ============================================
  implicit none

  ! ============================================
  ! [2] Define HERE any additional parameters or variables required by
  ! your simulation, unless they are provided by an external module

  ! We declare the two wind parameters objects
  type(SphericalWindType) :: spherical_wind
  type(PlaneWindType) :: plane_wind

  ! ============================================

contains

  subroutine userInitialCondition (uvars)
  ! ============================================
  ! [3] USER-DEFINED INITIAL CONDITIONS
  !
  !< @brief User-defined Initial Conditions
  !< @details This subroutine is called at the beginning of the simulation,
  !! after the base grid is built and a basic uniform initial condition
  !! is imposed. It is to be modified by the user to define the problem-
  !! specific Initial Condition.
  !!
  !! IMPORTANT: This subroutine receives the FLOW variables array to be 
  !! modified as argument 'uvars'. The subroutine must modify this array,
  !! *NOT* the global arrays U, UP or PRIMS.
  !!
  !! The array has the following structure:
  !!   uvars ( block ID, equation number, cell_i, cell_j, cell_k )
  !! where the equation numbering is:
  !!   1: rho
  !!   2: rho*u
  !!   3: rho*v
  !!   4: rho*w
  !!   5: E (kinetic+thermal)
  !! If the passive magnetic field is enabled:
  !!   6: B_x
  !!   7: B_y
  !!   8: B_z
  !! If passive scalars are enabled, they begin after all other flow
  !! variables (i.e., 9+ if passive B-field enabled, 6+ if not).
  !!
  !! Note that the cell indices include ghost cells. For instance, if we
  !! had a 1-deep ghost cell layer then cells from 1 to ncells_x would be
  !! physical cells, while cells 0 and ncells_x+1  would be ghost cells.
  ! ============================================

    implicit none
    real, intent(inout) :: uvars (nbMaxProc, neqtot, &
                           nxmin:nxmax, nymin:nymax, nzmin:nzmax)

    ! ============================================

    ! First we fill in the stellar wind parameters
    spherical_wind%xc = 0.5 * PC
    spherical_wind%yc = 0.5 * PC
    spherical_wind%zc = 0.5 * PC
    spherical_wind%radius = 0.1 * PC
    spherical_wind%mdot = 3e-7 * MSUN/YR
    spherical_wind%vinf = 5 * KPS
    spherical_wind%temp = 100
    spherical_wind%mu = 1.3

    ! For the plane wind we just copy the ISM parameters
    plane_wind%plane = PLANE_TOP
    plane_wind%rho = ism_dens
    plane_wind%vel = abs(ism_vz)
    plane_wind%temp = ism_temp
    plane_wind%mu = ism_mu0

    ! Then we impose the wind sources
    call imposeSphericalWind (spherical_wind, uvars)
    call imposePlaneWind (plane_wind, uvars)

    ! ============================================
    
  end subroutine userInitialCondition

  !=============================================================================

  subroutine userBoundary (uvars)
  ! ============================================
  ! [4] USER-DEFINED BOUNDARY CONDITIONS
  !
  !< @brief User-defined Boundary Conditions
  !< @details This subroutine is called once per timestep *after* standard
  !! boundary have been applied to all blocks. It allows the user to
  !! to impose an arbitrary boundary condition on the simulation.
  !!
  !! IMPORTANT: This subroutine receives the FLOW variables array to be 
  !! modified as argument 'uvars'. The subroutine must modify this array,
  !! *NOT* the global arrays U and UP.
  !!
  !! The structure of this array is described in the userInitialConditions()
  !! subroutine documentation above.
  ! ============================================

    implicit none
    real, intent(inout) :: uvars (nbMaxProc, neqtot, &
                           nxmin:nxmax, nymin:nymax, nzmin:nzmax)
    ! ============================================
    
    ! After each boundary exchange, we refresh the wind sources
    call imposeSphericalWind (spherical_wind, uvars)
    call imposePlaneWind (plane_wind, uvars)
    
    ! ============================================

  end subroutine userBoundary

end module userconds
