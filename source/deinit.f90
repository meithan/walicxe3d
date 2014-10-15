!===============================================================================
!> @file deinit.f90
!> @brief Cleans up and deinitializes
!> @author Juan C. Toledo
!> @date 3/Jun/2011

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

!> @brief Cleans up before exiting program
!> @details Deinitializations
subroutine deinit ()

  use parameters
  use globals
  use tictoc
  implicit none

  ! Deallocate dynamic arrays
  write(logu,'(a)') "Deallocating memory ..."  
  deallocate (U)
  deallocate (UP)
  deallocate (PRIM)
  deallocate (globalBlocks)
  deallocate (localBlocks)
  deallocate (refineFlags)
  deallocate (nblockslev)
  deallocate (nbx)
  deallocate (nby)
  deallocate (nbz)
  deallocate (minID)
  deallocate (maxID)
  deallocate (dx)
  deallocate (dy)
  deallocate (dz)

#ifdef MPIP
  write(logu,'(a)') "Finalizing MPI ..."  
  call mpi_finalize (ierr)
#endif

  if (logged) then
    close(logu)
  end if

end subroutine deinit

!===============================================================================

!> @brief Aborts the execution cleanly, returning an error code
!> @details The error code will be one of the constants defined in
!! constants.f90
subroutine clean_abort (errcode)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: errcode

  ! Close log file to make sure all output is flushed to disk
  if (logged) then
    close(unit=logu)
  end if

#ifdef MPIP
  call mpi_abort (MPI_COMM_WORLD, errcode, ierr)
  call mpi_finalize (ierr)
#endif
  stop

end subroutine clean_abort
