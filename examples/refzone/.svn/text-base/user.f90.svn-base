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

! EXAMPLE: refzone
! This example shows how to ask the code to refine a specific zone of
! the mesh to a specified level.

! The zone refinement routine is called twice in order to refine two zones of
! the mesh to the maximum level: a 1 pc wide zone centered on (2.5,2.5,2.5)
! and another 1 pc wide zone centered on (8.5,8.5,8.5). Note how the mesh is
! also modified outside these zones to conform to the level proximity rules.

! The zone refinement routine is called here as part of the initial conditions
! (userInitialConditions subroutine). The 6-element array argument 'zone'
! contains the coordinates of the bounding box of the zone to be refined:
!  zone(1): left edge (low X)
!  zone(2): right edge (high X)
!  zone(3): front edge (low Y)
!  zone(4): back edge (high Y)
!  zone(5): bottom edge (low Z)
!  zone(6): top edge (high Z)
! These must be given in cgs (cm).
! The argument 'zlevel' is the level to which the zone is to be refined. In
! this case the global parameter maxlev has been used.

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



  ! ============================================
  implicit none

  ! ============================================
  ! [2] Define HERE any additional parameters or variables required by
  ! your simulation, unless they are provided by an external module



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
  !! Note that the cell indices include ghost cells. For instance, cells
  !! from 1 to ncells_x are physical cells, while cells 0 and ncells_x+1
  !! would be ghost cells.
  ! ============================================

    implicit none
    real, intent(inout) :: uvars (nbMaxProc, neqtot, &
                           nxmin:nxmax, nymin:nymax, nzmin:nzmax)

    ! ============================================

    real :: zone(6)

    zone(1) = 2.0*PC
    zone(2) = 3.0*PC
    zone(3) = 2.0*PC
    zone(4) = 3.0*PC
    zone(5) = 2.0*PC
    zone(6) = 3.0*PC
    call refineZone (zone, maxlev)

    zone(1) = 8.0*PC
    zone(2) = 9.0*PC
    zone(3) = 8.0*PC
    zone(4) = 9.0*PC
    zone(5) = 8.0*PC
    zone(6) = 9.0*PC
    call refineZone (zone, maxlev)

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
    
   
    
    ! ============================================

  end subroutine userBoundary

end module userconds
