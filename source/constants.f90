!===============================================================================
!> @file constants.f90
!> @brief Generic constants
!> @author Juan C. Toledo
!> @date 10/Apr/2011

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

!> @brief Physical, astronomical and generic named constants
!> @details This file contains a list of constants, which include physical and
!! astrophysical constants, as well of an assortment of named constants used
!! throughout the code for easier reading
module constants

  implicit none

  ! ============================================

  ! Fundamental and astrophysical constants (cgs)
  real, parameter :: AMU  = 1.660538782e-24   !< Atomic Mass Unit
  real, parameter :: KB   = 1.380650400e-16   !< Boltzmann constant
  real, parameter :: PC   = 3.085677588e+18   !< Parsec
  real, parameter :: AU   = 1.495978707e+13   !< Astronomical unit
  real, parameter :: YR   = 3.155673600e+7    !< Year (Earth, sidereal)
  real, parameter :: KYR  = 3.155673600e+10   !< One thousand years
  real, parameter :: MSUN = 1.988920000e+33   !< Solar mass
  real, parameter :: KPS  = 1.0e5             !< km/s in cgs
  real, parameter :: PI   = 3.14159265358979  !< Ratio of perimeter to diameter

  ! Named Constants

  ! Neighbor directions
  integer, parameter :: LEFT   = 1       !< Direction LEFT (low x)
  integer, parameter :: RIGHT  = 2       !< Direction RIGHT (hi x)
  integer, parameter :: FRONT  = 3       !< Direction FRONT (low y)
  integer, parameter :: BACK   = 4       !< Direction BACK (hi y)
  integer, parameter :: BOTTOM = 5       !< Direction BOTTOM (low z)
  integer, parameter :: TOP    = 6       !< Direction TOP (hi z)

  ! Dimensions
  integer, parameter :: DIM_X = 1     !< X-dimension
  integer, parameter :: DIM_Y = 2     !< Y-dimension
  integer, parameter :: DIM_Z = 3     !< Z-dimension

  ! BC types
  integer, parameter :: BC_REFLECTIVE = 1     !< Reflection (velocity) BC
  integer, parameter :: BC_FREE_FLOW  = 2     !< Free flow (zero gradient) BC
  integer, parameter :: BC_PERIODIC   = 3     !< Periodic BC

  ! Neighbor Types
  integer, parameter :: NEIGH_SAME     = 1    !< Same-level neighbor
  integer, parameter :: NEIGH_COARSER  = 2    !< Coarser neighbor
  integer, parameter :: NEIGH_FINER    = 3    !< Finer neighbor
  integer, parameter :: NEIGH_BOUNDARY = 4    !< Computational domain boundary

  ! Refinement Flags
  integer, parameter :: FLAG_NONE   = 0    !< No refinement flag
  integer, parameter :: FLAG_REFINE = 1    !< Flagged for Refinement
  integer, parameter :: FLAG_COARSE = 2    !< Flagged for Coarsening

  ! Numerical solvers
  integer, parameter :: SOLVER_LAX  = 1    !< Lax-Friedrichs solver
  integer, parameter :: SOLVER_HLL1 = 2    !< First-order HLL Riemann solver
  integer, parameter :: SOLVER_HLL  = 3    !< Full (2nd-order) HLL Riemann solver
  integer, parameter :: SOLVER_HLLC = 4    !< HLLC Riemann solver

  ! Slope limiters
  integer, parameter :: LIMITER_NONE     = 0     !< No slope limiter
  integer, parameter :: LIMITER_MINMOD   = 1     !< Minmod limiter
  integer, parameter :: LIMITER_VANLEER  = 2     !< Falle limiter
  integer, parameter :: LIMITER_ALBADA   = 3     !< Van Albada limiter
  integer, parameter :: LIMITER_UMIST    = 4     !< UMIST limiter
  integer, parameter :: LIMITER_WOODWARD = 5   !< Woodoward limiter
  integer, parameter :: LIMITER_SUPERBEE = 6   !< Superbee limiter

  ! Cooling algorithms
  integer, parameter :: COOL_NONE = 0      !< No radiative cooling
  integer, parameter :: COOL_TABLE = 1     !< Tabulated cooling, Lambda(T)
  integer, parameter :: COOL_TABLE_METAL = 2  !< Tabulated cooling, Lambda(T,Z)

  ! Output concurrency types
  integer, parameter :: OUT_SIMULT = 0
  integer, parameter :: OUT_TURNS = 1

  ! Cell ranges on blocks
  integer, parameter :: CELLS_ALL   = 0
  integer, parameter :: CELLS_PHYS  = 1
  integer, parameter :: CELLS_GHOST = 2

  ! Output units types
  integer, parameter :: PHYS_UNITS = 1
  integer, parameter :: CODE_UNITS = 2

  ! Mesh creation methods
  integer, parameter :: MESH_AUTO   = 1
  integer, parameter :: MESH_MANUAL = 2

  ! Cartesian axes (for cut extractor)
  integer, parameter :: AXIS_X = 1
  integer, parameter :: AXIS_Y = 2
  integer, parameter :: AXIS_Z = 3

  ! Error codes
  integer, parameter :: ERROR_GENERIC = 255
  integer, parameter :: ERROR_WRONG_NPROCS = 1
  integer, parameter :: ERROR_NO_LOGFILE = 2
  integer, parameter :: ERROR_INSUFFICIENT_RAM = 3
  integer, parameter :: ERROR_NOALLOC_BIGARR = 4
  integer, parameter :: ERROR_LOCAL_BID_NOT_FOUND = 5
  integer, parameter :: ERROR_ALREADY_MAX_LEV = 6
  integer, parameter :: ERROR_INSUFICIENT_NBMAXPROC = 7
  integer, parameter :: ERROR_REGISTER_CHILD = 8
  integer, parameter :: ERROR_REFINING_PAST_MAX = 9
  integer, parameter :: ERROR_REGISTER_FATHER = 10
  integer, parameter :: ERROR_INVALID_LOC_INDEX = 11
  integer, parameter :: ERROR_TOO_MANY_LEVS = 12
  integer, parameter :: ERROR_BASEGRID_BAD_ASPECT = 13
  integer, parameter :: ERROR_COOLING_LOAD_COEFS = 14
  integer, parameter :: ERROR_HILBERT_ORDER = 15
  integer, parameter :: ERROR_LOADBAL_NO_LOCAL_BID = 16
  integer, parameter :: ERROR_LOADBAL_INSUFFICIENT_RAM = 17
  integer, parameter :: ERROR_OUTPUT_FILE_OPEN = 18
  integer, parameter :: ERROR_PUT = 19
  integer, parameter :: ERROR_WARM_FILE_OPEN = 20
  integer, parameter :: ERROR_WARM_READ_BLOCKS = 21
  integer, parameter :: ERROR_INVALID_CELL_RANGE = 22
  integer, parameter :: ERROR_DIVISION_BY_ZERO = 23
  integer, parameter :: ERROR_TIMESTEP = 24
  integer, parameter :: ERROR_BAD_OUTPUT_MODE = 25
  integer, parameter :: ERROR_NOT_ENOUGH_PASSIVES = 26

contains

  ! Boundary conditions names
  character(10) function bcname (bc_type)

    integer, intent(in) :: bc_type

    if (bc_type.eq.BC_REFLECTIVE) then
      bcname = 'Reflective'
    else if (bc_type.eq.BC_FREE_FLOW) then
      bcname = 'Free flow'
    else if (bc_type.eq.BC_PERIODIC) then
      bcname = 'Periodic'
    else
      bcname = 'Unknown'
    end if

  end function bcname

  ! Hydro solver names
  character(30) function solvername (solver)

    integer, intent(in) :: solver

    if (solver.eq.SOLVER_LAX) then
      solvername = 'Lax-Friedrichs'
    else if (solver.eq.SOLVER_HLL1) then
      solvername = '1st-order HLL Riemann solver'
    else if (solver.eq.SOLVER_HLL) then
      solvername = 'HLL Riemann solver'
    else if (solver.eq.SOLVER_HLLC) then
      solvername = 'HLLC Riemann solver'
    else
      solvername = 'Unknown'
    end if

  end function solvername

  ! Limiter names
  character(32) function limitername (lim)

    integer, intent(in) :: lim

    if (lim.eq.LIMITER_NONE) then
      limitername = 'Arithmetic average (no limiter)'
    else if (lim.eq.LIMITER_MINMOD) then
      limitername = 'Minmod (most diffusive)'
    else if (lim.eq.LIMITER_VANLEER) then
      limitername = 'van Leer'
    else if (lim.eq.LIMITER_ALBADA) then
      limitername = 'Van Albada'
    else if (lim.eq.LIMITER_UMIST) then
      limitername = 'UMIST'
    else if (lim.eq.LIMITER_WOODWARD) then
      limitername = 'Woodward'
    else if (lim.eq.LIMITER_SUPERBEE) then
      limitername = 'Superbee (least diffusive)'
    else
      limitername = 'Unknown'
    end if

  end function limitername

end module
