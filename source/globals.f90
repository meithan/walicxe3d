!===============================================================================
!> @file globals.f90
!> @brief Global variables module
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

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.

!===============================================================================

!> @brief Declaration of global variables
!> @details This module declares all the global variables used by the code
!! (as few as possible) and is included by some subroutines. These globals
!! include the integration variables for every cell and the adaptive mesh
!! information.
module globals

  implicit none

  !> Rank of this process
  integer :: rank

  ! Mesh globals
  integer :: nbmin         !< Global index of first block
  integer :: nbmax         !< Global index of last block
  integer :: nbrootx       !< Number of root blocks along x
  integer :: nbrooty       !< Number of root blocks along y
  integer :: nbrootz       !< Number of root blocks along z
  integer :: maxcells_x    !< Number of max-resolution cells along x
  integer :: maxcells_y    !< Number of max-resolution cells along y
  integer :: maxcells_z    !< Number of max-resolution cells along z
  integer :: maxlev        !< Number of refinement levels
  integer :: nbLocal       !< Number of *local* blocks (per process)
  integer :: nbActive      !< Number of *global* blocks (all processes)

  ! Simulation state globals
  integer :: it            !< Iteration number
  real :: time             !< Time at start of current timestep
  real :: dt               !< Current numerical timestep
  integer :: nextout       !< Number of next output

  ! Generic globals
  integer :: start_mark     !< Timing mark (start of simulation)
  integer :: logu           !< Logfile unit number
  character(80) :: logfile  !< Logfile name
  character(15) :: host     !< The host on which the code runs

  ! MPI globals
#ifdef MPIP 
  integer :: ierr           !< Standard MPI return error
#endif

  ! Big Data Arrays
  !> @brief Flow variables
  !> @details The big arrays (U, UP and P) have all the following structure.
  !! @li  U/UP ( blockID, ieq, i, j, k )
  !! @n For the flow variables (U/UP), this is:
  !! @li U(:,1,i,j,k) : @f$\rho@f$
  !! @li U(:,2,i,j,k) : @f$\rho u@f$
  !! @li U(:,3,i,j,k) : @f$\rho v@f$
  !! @li U(:,4,i,j,k) : @f$\rho w@f$
  !! @li U(:,5,i,j,k) : @f$E@f$ (kinetic+thermal)
  !! @li U(:,6,i,j,k) : @f$B_x@f$
  !! @li U(:,7,i,j,k) : @f$B_y@f$
  !! @li U(:,8,i,j,k) : @f$B_z@f$
  !! @li U(:,9+,i,j,k) : passive scalars
  !! @li -OR-
  !! @li U(:,6+,i,j,k) : passive scalars
  !! @n and the primitive variables (PRIM):
  !! @li PRIM(:,1,i,j,k) : @f$\rho@f$
  !! @li PRIM(:,2,i,j,k) : @f$u@f$
  !! @li PRIM(:,3,i,j,k) : @f$v@f$
  !! @li PRIM(:,4,i,j,k) : @f$w@f$
  !! @li PRIM(:,5,i,j,k) : @f$P@f$
  !! @li PRIM(:,6,i,j,k) : @f$B_x@f$
  !! @li PRIM(:,7,i,j,k) : @f$B_y@f$
  !! @li PRIM(:,8,i,j,k) : @f$B_z@f$
  !! @li PRIM(:,9+,i,j,k) : passive scalars
  !! @li -OR-
  !! @li PRIM(:,6+,i,j,k) : passive scalars
  real, allocatable :: U(:,:,:,:,:)
  real, allocatable :: UP(:,:,:,:,:)      !< Stepped flow variables
  real, allocatable :: PRIM(:,:,:,:,:)    !< Primitives

  ! Auxiliary data arrays
  real, allocatable :: F (:,:,:,:)      !< Physical flux along x
  real, allocatable :: G (:,:,:,:)      !< Physical flux along y
  real, allocatable :: H (:,:,:,:)      !< Physical flux along z
  real, allocatable :: FC (:,:,:,:)     !< Numerical flux along x
  real, allocatable :: GC (:,:,:,:)     !< Numerical flux along y
  real, allocatable :: HC (:,:,:,:)     !< Numerical flux along z

  !> @brief Global block registry: list of all blocks (bIDs)
  !> @details This list is shared and synchronized across all processes.
  !! It is never modified directly. The index of a bID in this list
  !! indicates its owner and is the same as the index in the refineFlags list.
  integer, allocatable :: globalBlocks(:)

  !> @brief Global list of refinement flags
  !! This list is shared and synchronized across all processes
  !! A process is only allowed to modify the sublist betwen nbmin and nbmax
  integer, allocatable :: refineFlags(:)

  !> @brief Local block registry: list of blocks (bIDs) on this process
  !> @details This list is modified by each process, and then combined
  !! into globalBlocks. The index of a bID in this list is the same as
  !! the index in the U, UP and PRIM arrays
  integer, allocatable :: localBlocks(:)

  !> @brief Number of block at each refinement level
  integer, allocatable :: nblockslev(:)

  !> @brief Number of blocks per dimension at each level
  integer, allocatable :: nbx(:)
  integer, allocatable :: nby(:)
  integer, allocatable :: nbz(:)

  !> @brief bIDs of the first and last blocks in each level
  integer, allocatable :: minID(:)
  integer, allocatable :: maxID(:)

  !> @brief Grid Spacings
  real, allocatable :: dx(:)     !< Grid spacings along X (for every level)
  real, allocatable :: dy(:)     !< Grid spacings along Y (for every level)
  real, allocatable :: dz(:)     !< Grid spacings along Z (for every level)  

  !> @brief Radiative Cooling
  real, allocatable :: cooldata(:,:)
  integer :: coolpts
  real :: cool_Tmin, cool_Tmax

  ! List of timing marks
  ! 1: Initialization: big array allocation and main initializations
  ! 2: Initialization: basegrid
  ! 3: Load balance
  ! 4: Boundary exchange
  ! 5: Solver: Timestep calculation
  ! 6: Solver: Primitives calculation
  ! 7: Solver: Integration
  ! 8: Disk output
  real :: timings(3)

end module globals

