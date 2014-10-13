!===============================================================================
!> @file tictoc.f90
!> @brief Benchmarking module
!> @author Juan C. Toledo
!> @version 1.5
!> @date 14/Dec/2011

! Copyright (c) 2014 Juan C. Toledo
!
! This program is free software; you can redistribute it and/or modify
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

!> @brief Benchmarking functions
!> @details Contains functions and subroutines designed to measure
!! the execution time between points in the program execution.
!> @section changelog Changelog
!! @li @b 1.5
!! @c TOC now returns a real with the number of elapsed seconds, and two
!! new functions were added: @c FTIME accepts a real representing a
!! time span in seconds, and returns a character string with hours, minutes
!! and seconds, as needed. @c NICETOC receives an @c INTEGER time mark, passes it
!! to toc, and then returns the result of timeformat, i.e., a formatted
!! character string represented the time elapsed.
!! @li @b 1.4
!! Some subroutines were converted into functions that return
!! a character string instead of printing the result to screen.
!! @li @b 1.3
!! The subroutine @c ETA was added. It allows to estimate the remaining
!! execution time based on a user provided progress estimator.
!! @li @b 1.2
!! The intrinsic @c SYSTEM_CLOCK is now used instead of
!! @c DATE_AND_TIME to calculate time intervals. See Comments.
!! @li @b 1.1
!! @c STAMP was added to generate nicely formatted timestamps.
!! @li @b 1.0
!! First version.
!> @section comments Comments
!! Starting from version 1.2, the module utilizes @c SYSTEM_CLOCK
!! instead of @c DATE_AND_TIME because it provides higher time
!! resolution. The disadvantage of this is that @c SYSTEM_CLOCK
!! has a relatively small looparound period, which is compiler-dependent.
!! For the GCC 4.2.1, this value is 24.855 days, so trying to measure
!! a time interval larger than this will produce incorrect results.
!! @n @n This program measures wall time.
!! @section usage Module Usage
!! Time measurements are done with the function pair @c TIC and @c TOC.
!! One must first call @c TIC to mark a certain point in the execution
!! of the program, and then @c TOC can be called at any later time in
!! order to calculate the elapsed time since the moment marked by @c TIC.
!! @n @n
!! Subroutine @c TIC modifies the given @c INTEGER argument @c mark storing
!! information about the point in time when it was called. Multiple
!! time marks (one variable per mark) can be used in the program to mark
!! different points. Time resolution is system dependent but is generally
!! on the order of milliseconds. 
!! @n @n
!! The @c ETA function estimates the execution time remaining. It must
!! be given a real number between 0 and 1 that indicates the percentage
!! progress of the execution at the moment of the call, as well as a time
!! mark obtained at the program start.
!! @n @n
!! The @c STAMP function returns a string with the current date and time.
!! It is common to call @c STAMP together with a call to @c TOC.

! ==============================================================================

module tictoc

  implicit none

  contains

! ==============================================================================

  !> @brief Creates a time mark
  !> @details Used to create a time mark at certain point in the execution.
  !> @retval mark An @c INTEGER representing a time mark
  subroutine tic (mark)    
    integer, intent(out) :: mark
    
    call system_clock(mark)

  end subroutine

! ==============================================================================

  !> @brief Calculates the number of seconds elapsed since a given time mark
  !> @details Computes the elapsed number of seconds since a time mark
  !! represented as an @c INTEGER previously returned by @c TIC.
  !> @return A @C REAL value with the number of elapsed seconds
  real function toc (mark)
    integer, intent(in) :: mark

    integer :: now, RATE, CMAX

    CALL system_clock(now,RATE,CMAX)

    toc = (now-mark)/real(RATE)
    return
    
  end function

! ==============================================================================

  !> @brief A 'nice' version of toc which returns the time in hours, mins, secs
  !> @details Calls toc and then post-processes the result through timeformat.
  !> @return A character string with the elapsed hours, minutes and seconds.
  character(len=30) function nicetoc (mark)
    integer, intent(in) :: mark

    real :: secs
    character(len=30) :: buffer

    secs = toc(mark)
    buffer = ftime(secs)
    nicetoc = buffer
    return
    
  end function

! ==============================================================================

  !> @brief Returns a nice time representation in hours, minutes and seconds
  !> @details Accepts a @C REAL value with a number of seconds and returns
  !! a character string with the time formatted as hours, minutes, seconds.
  !> @return A character string representingt the time laps in h, m, s
  character(len=30) function ftime (seconds)
    real, intent(in) :: seconds

    integer :: days, hours, mins
    real :: secs
    character(len=30) :: buffer

    secs = seconds
    days = int(secs/86400)
    secs = seconds - days*86400
    hours = int(secs/3600)
    secs = secs - hours*3600
    mins = int(secs/60)
    secs = secs - mins*60

    if (seconds.lt.60.0) then
      write(buffer,'(f6.3,a)') secs, 's'
    elseif (seconds.lt.3600.0) then
      write(buffer,'(i2,a,f6.3,a)') mins, 'm ', secs, 's'
    elseif (seconds.lt.86400.0) then
      write(buffer,'(i2,a,i2.2,a,f6.3,a)') hours, 'h ', mins, 'm ', secs, 's'
    else
      write(buffer,'(i3,a,i2.2,a,i2.2,a,F6.3,a)') &
        days, 'd ', hours, 'h ', mins, 'm ', secs, 's'
    end if

    ftime = buffer
    return
    
  end function

! ==============================================================================

  !> @brief Computes the estimated remaining time
  !> @details Returns the estimated remaining time of execution.
  !! This requires a time mark previously returned by @c TIC marking
  !! the start of execution (or, alternatively, a zero value), as well
  !! as a progress value given as a real number between 0 and 1.
  !> @return A character string with the hours, minutes and seconds remaining.
  character(len=20) function eta (mark, progress)
    integer, intent(in) :: mark
    real, intent(in) :: progress

    integer :: now, elapsed, RATE, CMAX
    integer :: hours, mins
    real :: remaining, secs
    character(len=20) :: buffer

    call system_clock(now,RATE,CMAX)

    elapsed = now - mark
    remaining = elapsed/progress - elapsed

    hours = int(remaining/real(RATE)/3600.)
    mins = int(remaining/real(RATE)/60.) - hours*60
    secs = remaining/real(RATE) - hours*3600 - mins*60

    if ((hours==0).and.(mins==0)) then
      201 format(F6.3,'s')
      write(buffer,201) secs
    elseif (hours==0) then
      202 format(I2,'m ',F6.3,'s')
      write(buffer,202) mins, secs
    else
      203 format(I4,'h ',I2,'m ',F6.3,'s')
      write(buffer,203) hours, mins, secs
    end if

    eta = buffer
    return

    end function

! ==============================================================================

  !> @brief Creates a timestamp
  !> @details Returns the current date and time in a nice format.
  !> @return A character string with the current date and time.  
  character(len=22) function stamp ()
    
    ! Local variables
    character(len=22) :: buffer
    integer :: now(8)
    character(8) :: datestr
    character(10) :: time
    character(5) :: zone
    character(3) :: month

    call date_and_time(datestr,time,zone,now)

    select case (datestr(5:6))
      case ('01')
        month = 'Jan'
      case ('02')
        month = 'Feb'
      case ('03')
        month = 'Mar'
      case ('04')
        month = 'Apr'
      case ('05')
        month = 'May'
      case ('06')
        month = 'Jun'
      case ('07')
        month = 'Jul'
      case ('08')
        month = 'Aug'
      case ('09')
        month = 'Sep'
      case ('10')
        month = 'Oct'
      case ('11')
        month = 'Nov'
      case ('12')
        month = 'Dec'
    end select

    100 format(A,'/',A,'/',A,' @ ',A,':',A,':',A)
    write(buffer,100) datestr(7:8), month, datestr(1:4), time(1:2),&
                      time(3:4), time(5:6)

    stamp = buffer
    return
    
  end function


end module
