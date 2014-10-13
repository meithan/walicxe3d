!===============================================================================
!> @file cooling.f90
!> @brief Radiative cooling
!> @author Juan C. Toledo
!> @date 8/May/2012

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

!> @brief High-level wrapper subroutine for radiative cooling
subroutine cooling

  use parameters
  use globals
  use tictoc
  implicit none

  integer :: mark, nb, bID
  real :: maxcool
 
  if (cooling_type.ne.COOL_NONE) then 

    call tic(mark)
    write(logu,*) ""
    write(logu,'(1x,a)') "============================================"
    write(logu,'(1x,a)') " Applying Radiative Cooling ..."
    write(logu,'(1x,a)') "============================================"
    write(logu,*) ""

    select case (cooling_type)

    case (COOL_TABLE)

      ! Apply tabulated cooling to all local blocks
      do nb=1,nbMaxProc
        bID = localBlocks(nb)
        if (bID.ne.-1) then
          call apply_cooling (nb, maxcool)
!          if (maxcool.ge.cooling_warning) then
!            write(logu,'(1x,a,i8,a)') "Cooling warning for block ", bID, " !"
!            write(logu,'(1x,a,f6.3)') "Max thermal energy loss= ", maxcool
!          end if
        end if
      end do

    end select

    write(logu,'(1x,a,a)') "> Cooling applied in ", nicetoc(mark)

  end if

end subroutine cooling

!===============================================================================

!> @brief Loads radiative cooling coefficients from a data file
!> @details Loads radiative cooling coefficients from a data file into the
!! cooldata global array. The number of data points must defined in the first
!! line of the data file, and every subsequent line must have a pair of data
!! values: log10(temp) -1*log10(coef).
subroutine loadcooldata ()

  use parameters
  use globals
  implicit none

  integer :: i, istat
  real :: a, b

  ! Master loads cooling table from file
  if (rank.eq.master) then

    write(logu,'(2x,a,a,a)') "Loading cooling data from file ", trim(cooling_file), " ..."
    open (unit=99, file=cooling_file, status='old', iostat=istat)
    if (istat.ne.0) then
      write(logu,'(a,a,a)') "Could not open the file ", trim(cooling_file), " !"
      write(logu,*) "***ABORTING***"
      close(99)
      call clean_abort (ERROR_COOLING_LOAD_COEFS)
    end if

    read(99,*) coolpts
    allocate( cooldata(2,coolpts) )

    do i=1,coolpts
      read(99,*) a, b
      cooldata(1,i) = 10.0**a
      cooldata(2,i) = 10.0**(-b)
! DEBUG
write(logu,*) 10.0**a, 10.0**(-b)
! DEBUG          
    end do
    close (unit=99)

  end if

  ! Broadcast cooling data
  call mpi_bcast(coolpts, 1, mpi_integer, 0, mpi_comm_world, ierr)
  if (rank.ne.master) then
    write (logu,'(2x,a)') "Receiving cooling data from master process ..."
    allocate (cooldata(2,coolpts))
  end if
  call mpi_bcast(cooldata, coolpts*2, mpi_real_kind, 0, mpi_comm_world, ierr)
  write(logu,'(2x,a,i0,a)') "Loaded ", coolpts, " cooling coefficients."

  ! Set minimum and maximum temperatures
  cool_Tmin = cooldata(1,1)
  cool_Tmax = cooldata(1,coolpts)
  
end subroutine loadcooldata

!===============================================================================

!> @brief Applies radiative cooling to a block
!> @details This employs the previously loaded tabulated cooling function.
!! Cooling is turned off below for temperatures below the minimum value
!! of the table (global var cool_Tmin), while above the maximum value
!! (cool_Tmax) free-free emission, propto sqrt(T), is used.
!! Cooling is applied by decreasing the pressure of the cell by the factor
!!   exp(-n^2*Lambda(T)*dt/e0)
!! where lambda(T) is the temperature-dependent cooling function (in units of
!! erg cm^3 s^-1), n is the number density, dt is the timestep and e0 is the
!! initial thermal energy density of the cell. This comes from taking the
!! energy density loss rate, -n^2*Lambda(T), and assuming Lambda(T) is approx.
!! linear throughout the small interval dt.
!! The cooling factor is limited to being at least as large as the parameter
!! cooling_limit. Also, a locally-defined temperature floor is enforced.
!> @param bIndx The local index of the block
subroutine apply_cooling (bIndx, maxcool)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bIndx
  real, intent(out) :: maxcool
  
  real, parameter :: T_floor = 10.0

  real :: temp, new_temp, vel2, ETH, EK, aloss, numdens, ce, cool_factor
  integer :: i, j, k

  maxcool = 0.0

  do i=1,ncells_x
    do j=1,ncells_y
      do k=1,ncells_z

        ! Calculate temperature of this cell
        call calcTemp (PRIM(bIndx,:,i,j,k), temp)

        ! Cooling not applied below T_min
        if (temp.ge.cool_Tmin) then

          ! Code units
          ETH = CV * PRIM(bIndx,5,i,j,k)
          vel2 = PRIM(bIndx,2,i,j,k)**2 + PRIM(bIndx,3,i,j,k)**2 + PRIM(bIndx,4,i,j,k)**2
          EK  = 0.5 * PRIM(bIndx,1,i,j,k) * vel2

          ! Calculate radiative loss and new temperature
          call find_coef (temp, cooldata, coolpts, aloss)
          numdens = PRIM(bIndx,1,i,j,k)*d_sc/(mui*AMU)  ! cgs, gas ionized
          ce = aloss*(numdens**2)/(ETH*e_sc)  ! cgs
          cool_factor = exp(-ce*(dt*t_sc))

          ! Limit cool_factor directly, if needed
          if (cool_factor.lt.1.0-cooling_limit) then
            cool_factor = 1.0-cooling_limit
          end if

          ! Floor temperature by adjusting cool_factor, if needed
          new_temp = temp * cool_factor
          if (new_temp.lt.T_floor) then
            new_temp = T_floor
            cool_factor = T_floor / temp
          end if

          ! Update pressure and total energy          
          PRIM(bIndx,5,i,j,k) = PRIM(bIndx,5,i,j,k) * cool_factor
          ETH = CV * PRIM(bIndx,5,i,j,k)
          U(bIndx,5,i,j,k) = EK + ETH

          ! Maximum cooling for this block (% of thermal energy lost)
          maxcool = max(maxcool, (1.0-cool_factor))

        end if

      end do
    end do
  end do

end subroutine apply_cooling

!===============================================================================

!> @brief Interpolates a cooling cofficient from tabulated data
!> @param T Temperature (Kelvin)
!> @param cooldata Tabulated cooling coefficients
!> @param numdata Size of coefficients table
!> @param coef Interpolated cooling coefficient
subroutine find_coef (T, cooldata, numdata, coef)

    implicit none

    real, intent(in) :: T
    integer, intent(in) :: numdata
    real, intent(in) :: cooldata(2,numdata)
    real, intent(out) :: coef

    integer :: i
    real :: Tmin, Tmax, T1, T2, C1, C2

    Tmin = cooldata(1,1)
    Tmax = cooldata(1,numdata)

    if (T.lt.Tmin) then

      coef = 0.0
      return
    
    else if (T.ge.Tmax) then
    
!      coef = 0.21D-26*sqrt(dble(T))
      coef = cooldata(2,numdata)*sqrt(real(T)/Tmax)
      return
    
    else

      do i=2,numdata
        if (cooldata(1,i).gt.T) then
          T1 = cooldata(1,i-1)
          C1 = cooldata(2,i-1)
          T2 = cooldata(1,i)
          C2 = cooldata(2,i)
          coef = (C2-C1)/(T2-T1)*(real(T)-T1)+C1 
          return
        end if
      end do
      
    end if

  end subroutine find_coef
  
