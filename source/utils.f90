!===============================================================================
!> @file utils.f90
!> @brief Generic subroutines
!> @author Juan C. Toledo
!> @date 9/Jun/2011

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

!> @brief Replaces substring a with subtring b in string str.
!> @details If a is not found within str, str is left unchanged.
!! It is up to the caller to make sure str is long
!! enough to acommodate b (when it is longer than a).
!> @param str The string in which to perform the substitution
!> @param a The substring to be replaced by 'b'
!> @param b The substring to replace 'a' with
subroutine replace (str,a,b)

  implicit none

  character(*), intent(inout) :: str
  character(*), intent(in) :: a
  character(*), intent(in) :: b

  integer :: i, s, e, l
  character(len(str)) :: buff

  l = len(a)
  s = index(str,a,.false.)
  e = s + l - 1

  if (s.ne.0) then
    buff = str(1:s-1) // b // str(e+1:len(str))
    do i=1,len(str)-1
      str(i:i) =  buff(i:i)
    end do
  end if

end subroutine

!===============================================================================

!> @brief Returns the position of the first appearance of an element in a list
!> @details 'indx' is returned with the element's index in the list.
!! If the element is not found, 'indx' returns as -1.
!! Warning: 'list' must be a 1D integer array, 'element' must be an integer.
!> @param element The element to be found
!> @param list The list (a 1D integer array)
!> @param length The length of the list
!> @param indx The returned position of the element within the list 
subroutine find (element, list, length, indx)

  implicit none

  integer, intent(in) :: element
  integer, intent(in) :: length
  integer, intent(in) :: list(length)
  integer, intent(out) :: indx

  integer :: nb

  indx = -1

  ! Linear search
  do nb=1,length
    if (list(nb).eq.element) then
      indx = nb
      return
    end if
  end do

end subroutine find

!===============================================================================

!> @brief Inserts an element in the first empty slot of a list
!> @details Empty slots are defined as having value -1. Parameter 'indx' is
!! returned with the position in which the item was inserted. If the item has
!! no empty slots, an indx of -1 is returned.
!! Warning: 'list' must be a 1D integer array, 'element' must be an integer.
!> @param element The element to be inserted
!> @param list The list (a 1D integer array)
!> @param length The length of the list
!> @param indx The returned position of the element within the list 
subroutine put (element, list, length, indx)

  implicit none

  integer, intent(in) :: element
  integer, intent(in) :: length
  integer, intent(inout) :: list(length)
  integer, intent(inout) :: indx

  indx = -1

  call find (-1, list, length, indx)
! DEBUG
  if (indx.eq.0) then
    print*, "PUT routine"
    print*, "index = 0 !!!!"
    stop
  end if
! DEBUG
  if (indx.ne.-1) then
    list(indx) = element
  else
! DEBUG
!  print*, "returning -1: ", indx
! DEBUG
    return
  end if

! DEBUG
!  print*, "returning: ", indx
! DEBUG

end subroutine put

!===============================================================================

!> @brief Deletes the first occurence of an element in a list
!> @details This will only delete the first appearance. The slot will be set
!! to -1. Parameter 'indx' will return the index where the element was found,
!! or -1 if the element was not found (and no deletion occured)
!! Warning: 'list' must be a 1D integer array, 'element' must be an integer.
!> @param element The element to be deleted
!> @param list The list (a 1D integer array)
!> @param length The length of the list
!> @param indx The position where the first occurence was found
subroutine pop (element, list, length, indx)

  implicit none

  integer, intent(in) :: element
  integer, intent(in) :: length
  integer, intent(inout) :: list(length)
  integer, intent(out) :: indx

  indx = -1

  call find (element, list, length, indx)

  if (indx.ne.-1) then
    list(indx) = -1
  else
    return
  end if

end subroutine pop

!===============================================================================

!> @brief Deletes all occurences of an element in a list
!> @details This will search and delete all occurrences of 'element'. The slot
!! will be set to -1. Parameter 'dels' will return with the number of occurences
!! that were deleted.
!! Warning: 'list' must be a 1D integer array, 'element' must be an integer.
!> @param element The element of which all occurrences are to be deleted
!> @param list The list (a 1D integer array)
!> @param length The length of the list
!> @param dels The numbers of time the element was found and deleted
subroutine popAll (element, list, length, dels)

  implicit none

  integer, intent(in) :: element
  integer, intent(in) :: length
  integer, intent(inout) :: list(length)
  integer, intent(out) :: dels

  integer :: indx

  dels = 0
  indx = -1

  ! Linear search
  do indx=1,length
    if (list(indx).eq.element) then
      list(indx) = -1
      dels = dels + 1
    end if
  end do

  return

end subroutine popAll

!===============================================================================

!> @brief Takes in a time in seconds and returns it as a nice string
!> @details It will scale the time unit used so as to produce a 
!! human-friendly output, in kyr, yr or days.
!> @param time A time to be converted to a nice string. MUST BE IN SECONDS.
!> @param timestr The returned string
subroutine nicetime (time, timestr)

  implicit none

  real, intent(in) :: time
  character(*), intent(out) :: timestr

  real, parameter :: YR = 3.15576e7

  if (time.ge.10000*YR) then
    write(timestr,'(f0.3,a)') time/(1000.0*YR), " kyr"
  else if (time.ge.1.0*YR) then
    write(timestr,'(f0.3,a)') time/YR, " yr"
  else if (time.ge.86400*1.0E-2) then
    write(timestr,'(f0.3,a)') time/86400.0, " d"
  else
    write(timestr,'(es9.2,a)') time, " s"
  end if

  return

end subroutine nicetime

!===============================================================================

