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

! EXAMPLE FILE: snr

! This example shows how to impose supernova remnants. We simulate a 10^3 pc
! box and impose two snr at the start (one with a "simple" model and the
! other with a more realistic Type Ia ejecta model). At t=500 yr we detonate
! a third snr in the center. We then let the simulation evole until 1 kyr.
! This is done using an equivalent max-level resolution of 128^3 and using
! four parallel processes.

! In the userInitialConditions subroutine we use imposeSNR and imposeSNRIa
! subroutines to detonate the initial remnants. These subroutines must be
! passed two arguments: a snr_params_type data object which contains the
! parameters of the snr, and the flow variables vector (uvars), which are
! received by the user IC routine from its caller.

! Note that the snr module is imported at the top of the file.

! The snr_params data object contains several fields that must be filled
! for the snr to be properly defined:
!  xc: x-coordinate of the center of the remnant (cm)
!  yc: y-coordinate of the center of the remnant (cm)
!  zc: z-coordinate of the center of the remnant (cm)
!  radius: the radius of the remnant (cm)
!  mass: the total mass of the ejecta (g)
!  energy: the total energy (kinetic+thermal) of the explosion (erg)
!  chi: the fraction of the total energy deposited as kinetic energy (0.0-1.0)
!  bx, by, bz: magnetic field components inside the remnant (G)
!  time: the time at which the snr is to be detonated. This is only for the
!   user to use if a snr is to be detonated at a later point in time.
!  armed: a logical flag that is set to .false. after the snr has been
!   detonated. This can be used by the user to prevent multiple detonations
!   of the same remnant.

! We also use the userBoundary subroutine, which is called at the end of the
! boundary conditions at each timestep, to detonate a SNR at a poin in time
! after the start of the simulation. Note how the time and armed fields of
! the snr parameters data object are used.

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
!!   2) You can define aditional module-wide global variables and parameters
!!      in the section marked [2].
!!   3) Fill in the subroutine userInitialCondition(), marked [3], which
!!      is called at the beginning of the simulation.
!!   4) Optionally, fill in the subroutine userBoundary(), marked [4],
!!      which is called at the end of each boundary exchange operation.
!!
!! All subroutines in this module automatically have access to the global
!! parameters and variables.
! ============================================

  use parameters
  use globals
  ! ============================================
  ! [1] Add HERE any aditional modules required by your simulation
 
  ! It's important to import the snr module here
  use snr

  ! ============================================
  implicit none

  ! ============================================
  ! [2] Define HERE any additional parameters or variables required
  ! by the user subroutines below if they they are not provided by
  ! an external module

  ! We can declare the snr parameters objects here. We'll fill them out
  ! during the initial conditions.
  type(snr_params_type) :: snr1
  type(snr_params_type) :: snr2
  type(snr_params_type) :: snr3

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

    ! Parameters of the first SNR
    snr1%xc = 5.0*PC
    snr1%yc = 5.0*PC
    snr1%zc = 2.5*PC
    snr1%radius = 1.0*PC
    snr1%mass = 4*MSUN
    snr1%energy = 1E51
    snr1%chi = 0.5
    snr1%time = 0.0

    ! Parameters of the second SNR
    snr2%xc = 5.0*PC
    snr2%yc = 5.0*PC
    snr2%zc = 7.5*PC
    snr2%radius = 1.0*PC
    snr2%mass = 1.4*MSUN
    snr2%energy = 1E51
    snr2%chi = 0.5
    snr2%time = 0.0

    ! Parameters of a third SNR, to be detonated at 500 yr
    snr3%xc = 5.0*PC
    snr3%yc = 5.0*PC
    snr3%zc = 5.0*PC
    snr3%radius = 1.0*PC
    snr3%mass = 1.4*MSUN
    snr3%energy = 1E51
    snr3%chi = 0.5
    snr3%time = 500*YR

    ! The first two remnants are detonated as initial conditions
    call detonateSNR(snr1, uvars)
    call detonateSNRIa(snr2, uvars)

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
    
    ! We check both for time and armed status, and detonate when 
    ! both become true. The armed property is set to false automatically
    ! by the subroutine (but can be reset by the user, if he so wishes).
    ! Also note how we de-scale the time variable to physical units.
    if ((time*t_sc.ge.snr3%time).and.(snr3%armed)) then
      call detonateSNRIa(snr3,uvars)
    end if
    
    ! ============================================

  end subroutine userBoundary

end module userconds
