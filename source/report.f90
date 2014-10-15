!===============================================================================
!> @file report.f90
!> @brief Subroutines to print progress reports
!> @author Juan C. Toledo
!> @date 22/Apr/2013

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

!> @brief Prints out the "main" progress report (after each iteration)
subroutine main_report ()

  use parameters
  use globals
  use tictoc
  implicit none

  character(30) :: timestr, dtstr

  call nicetime(time*t_sc, timestr)
  call nicetime(dt*t_sc, dtstr)

  ! Report progress
  write(logu,*) "============================================"
  write(logu,'(1x,a,i0,a)') "Iteration " , it, " complete!"
  write(logu,'(1x,a)') stamp()
  write(logu,'(1x,a,a)') "Elapsed: ", nicetoc(start_mark)
  write(logu,'(1x,a,a)') "time = ", trim(timestr)
  write(logu,'(1x,a,a)') "dt = ", trim(dtstr)
  write(logu,'(1x,i0,a)') nbLocal, " local blocks"
  write(logu,*) "============================================"
  write(logu,'(1x,a,i0,a,a,a,a,a,i0,a,i0)') "it ", it, "; t= ", &
    trim(timestr), "; dt=", trim(dtstr), "; blocks=", nbLocal, "/", nbActive
  write(logu,*) "============================================"
  write(logu,*) ""

end subroutine main_report





