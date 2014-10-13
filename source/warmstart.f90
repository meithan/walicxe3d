!===============================================================================
!> @file warmstart.f90
!> @brief Warm start module
!> @author Juan C. Toledo
!> @date 5/Jun/2012

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

!> @brief Performs a warm start. This requires reading the state file specified
!! in parameters.f90 and loading hydro data.
subroutine warmstart ()

  use parameters
  use globals
  implicit none

  integer :: unitin, istat, noutput, l, nb, nblocks
  character(256) :: datadir_old
  character(256) :: blocksfile
  character(1) :: slash
  character(4) :: noutstr
  character(3) :: rankstr

  write(logu,*) ""
  write(logu,*) "============================================"
  write(logu,'(1x,a)') " Performing warm start ..."
  write(logu,*) "============================================"
  write(logu,*) ""

  ! Open state file
  write(logu,'(1x,a,a,a)') "Reading state file '", trim(warm_file), "' ..."
  unitin = 10 + rank
  open (unit=unitin, file=warm_file, status='old', iostat=istat)
  if (istat.ne.0) then
    write(logu,'(a,a,a)') "Could not open the file '", trim(warm_file), "' !"
    close(unitin)
    call clean_abort (ERROR_WARM_FILE_OPEN)
  end if

  ! Read simulation state variables and datadir
  read(unitin,'(es22.15,i8,i5)') time, it, noutput
  read(unitin,'(a)') datadir_old

  write(logu,*) ""
  write(logu,'(1x,a,i0,a)') "> Restarting simulation from output ", noutput, " ..."
  write(logu,'(1x,a,i0)') "Last iteration = ", it
  write(logu,'(1x,a,es22.15)') "Current time = ", time
  write(logu,*) ""

  time = time/t_sc
  nextout = noutput + 1

  close(unitin)

  ! Generate filename based on templates for data file
  l = len_trim(datadir_old)
  if (datadir(l:l).ne.'/') then
    slash = '/'
  else
    slash = ''
  end if
  write(rankstr,'(I3.3)') rank
  write(noutstr,'(I4.4)') noutput
  blocksfile = blockstpl
  call replace (blocksfile, 'XXX', rankstr)
  call replace (blocksfile, 'YYYY', noutstr)
  write(blocksfile,'(a)') trim(datadir_old) // trim(slash) // trim(blocksfile) // ".bin"

  ! Open data file
  write(logu,'(1x,a,a,a)') "Reading data file '", trim(blocksfile), "' ..."
  unitin = 10 + rank
  open (unit=unitin, file=blocksfile, status='old', access='stream', iostat=istat)
  if (istat.ne.0) then
    write(logu,'(a,a,a)') "Could not open the file '", trim(blocksfile), "' !"
    write(logu,'(a,a,a)') "Does the datadir '", trim(datadir_old), "' exist?"
    close(unitin)
    call clean_abort (ERROR_WARM_FILE_OPEN)
  end if

  ! Read data
  localBlocks(:) = -1
  nblocks = 0
  read(unitin) nbLocal
  do nb=1,nbLocal
    read(unitin) localBlocks(nb)
    write(logu,'(1x,a,i0,a)') "Loading block ", localBlocks(nb), " ..."
    read(unitin) U(nb,:,1:ncells_x,1:ncells_y,1:ncells_z)
    nblocks = nblocks + 1
  end do

  if (nbLocal.ne.nblocks) then
    write(logu,'(a,i0,a,i0,a)') "Read ", nblocks, " blocks from file, but expected ", nbLocal, " !"
    write(logu,'(a)') "***Aborting***"
    call clean_abort (ERROR_WARM_READ_BLOCKS)
  else
    write(logu,'(1x,a,i0,a,i0,a)') "Loaded ", nblocks, " blocks; expected ", nbLocal, "."
  end if

  close (unitin)

end subroutine warmstart
