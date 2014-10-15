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
  real :: maxloss
 
  if (cooling_type.ne.COOL_NONE) then 

    call tic(mark)
    write(logu,*) ""
    write(logu,'(1x,a)') "============================================"
    write(logu,'(1x,a)') " Applying Radiative Cooling ..."
    write(logu,'(1x,a)') "============================================"
    write(logu,*) ""

    select case (cooling_type)

    case (COOL_TABLE, COOL_TABLE_METAL)

      ! Apply tabulated cooling to all local blocks
      do nb=1,nbMaxProc
        bID = localBlocks(nb)
        if (bID.ne.-1) then

          call apply_cooling (nb, maxloss)
          if (maxloss.ge.cooling_limit) then
            write(logu,'(1x,a,i0,a,f6.3,a,f6.3)') &
              "Cooling warning for block ", bID, "! Max loss ", &
              maxloss, ", limited to", cooling_limit
          end if

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
!! values containing the logarithms of temperature and the cooling rate, with
!! the latter multiplied by -1 so it's positive.
!! Example table:
!! 5
!! 4.0  23.4
!! 4.5  22.9
!! 5.0  22.5
!! 5.5  21.7
!! 6.0  20.8
subroutine loadcooldata ()

  use parameters
  use globals
  implicit none

  integer :: i, istat
  real :: a, b

  ! Master loads cooling table from file
  if (rank.eq.master) then

    write(logu,'(2x,a)') "Loading cooling data from file ..."
    open (unit=99, file=cooling_file, status='old', iostat=istat)
    if (istat.ne.0) then
      write(logu,'(a,a,a)') "Could not open the file ", trim(cooling_file), " !"
      write(logu,*) "***ABORTING***"
      close(99)
      call clean_abort (ERROR_COOLING_LOAD_COEFS)
    end if

    read(99,*) nptsT
    allocate( cooltable(2,nptsT) )

    do i=1,nptsT
      read(99,*) a, b
      cooltable(1,i) = a
      cooltable(2,i) = -b
      write(logu,*) 10.0**a, 10.0**(-b)
    end do
    close (unit=99)

  end if

  ! Broadcast cooling data to all processes
  call mpi_bcast(nptsT, 1, mpi_integer, 0, mpi_comm_world, ierr)
  if (rank.ne.master) then
    write (logu,'(2x,a)') "Receiving cooling data from master process ..."
    allocate (cooltable(2,nptsT))
  end if
  call mpi_bcast(cooltable, nptsT*2, mpi_real_kind, 0, mpi_comm_world, ierr)
  write(logu,'(2x,a,i0,a)') "Loaded ", nptsT, " cooling coefficients."

  ! Set global vars for minimum and maximum temperatures
  ! Note that the temperatures are logarithmic
  cool_Tmin = cooltable(1,1)
  cool_Tmax = cooltable(1,nptsT)
 
end subroutine loadcooldata

!===============================================================================

!> @brief Loads radiative cooling coefficients for 2D (T,Z) data
!> @details Loads radiative cooling coefficients in the (T,Z) plane from a
!! data file into the cooldata global array. The number of temperature and
!! metallicity data points must defined in the first line of the data file,
!! separated by blank. The second line must contain a blank-separated list
!! of metallicity values. Then, every subsequent line must contain, first,
!! the temperature, and then the corresponding cooling coefficients for each
!! metallicity value. The temperatures and cooling coefficients must be given
!! logarithmically, and are kept that way.
!! Example table:
!! 3 2
!!        1.0   10.0
!! 4.0  -23.4  -22.3
!! 4.5  -22.7  -21.4
!! 5.0  -22.1  -20.8
subroutine loadcooldata_metal ()

  use parameters
  use globals
  implicit none

  integer :: i, istat

  ! Master loads cooling table from file
  if (rank.eq.master) then

    write(logu,'(2x,a)') "Loading cooling data from file ..."
    open (unit=99, file=cooling_file, status='old', iostat=istat)
    if (istat.ne.0) then
      write(logu,'(a,a,a)') "Could not open the file ", trim(cooling_file), " !"
      write(logu,*) "***ABORTING***"
      close(99)
      call clean_abort (ERROR_COOLING_LOAD_COEFS)
    end if

    read(99,*,iostat=istat) nptsT, nptsZ
    if (istat.ne.0) then
      write(logu,'(a,a,a)') "Expected two values in first line of cooling &
        &file ", trim(cooling_file), " !"
      write(logu,'(a)') "Are you using a 2D temperature-metallicity table?"
      write(logu,*) "***ABORTING***"
      close(99)
      call clean_abort (ERROR_COOLING_LOAD_COEFS)
    end if
    allocate( cooltable(nptsT+1, nptsZ+1) )
    
    cooltable(1,1) = -1.0
    read(99,*) cooltable(1,2:)
    do i=2,nptsT+1
      read(99,*) cooltable(i,:)
    end do
    
    close(unit=99)

  end if

  ! Broadcast cooling data too all processes
  call mpi_bcast(nptsT, 1, mpi_integer, 0, mpi_comm_world, ierr)
  call mpi_bcast(nptsZ, 1, mpi_integer, 0, mpi_comm_world, ierr)  
  if (rank.ne.master) then
    write (logu,'(2x,a)') "Receiving cooling data from master process ..."
    allocate ( cooltable(nptsT+1,nptsZ+1) )
  end if
  call mpi_bcast(cooltable, (nptsT+1)*(nptsZ+1), mpi_real_kind, 0, &
    mpi_comm_world, ierr)
  write(logu,'(2x,a,i0,a,i0,a,i0,a)') "Loaded ", nptsT*nptsZ, &
    " cooling coefficients for ", nptsT, " temperature values and ", &
    nptsZ, " metallicity values."

  ! Set global vars for minimum/maximum temperatures and metallicities
  ! Note that the temperatures are logarithmic
  cool_Tmin = cooltable(2,1)
  cool_Tmax = cooltable(nptsT+1,1)
  cool_Zmin = cooltable(1,2)
  cool_Zmax = cooltable(1,nptsZ+1)
  
end subroutine loadcooldata_metal

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
!> @param maxloss Maximum fractional thermal energy loss observed
subroutine apply_cooling (bIndx, maxloss)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bIndx
  real, intent(out) :: maxloss
  
  real, parameter :: T_floor = 10.0

  real :: temp, new_temp, logT, vel2, ETH, EK, aloss, numdens, ce, cool_factor
  real :: frac_loss, metal
  integer :: i, j, k

  maxloss = 0.0

  do i=1,ncells_x
    do j=1,ncells_y
      do k=1,ncells_z

        ! Calculate temperature of this cell
        call calcTemp (PRIM(bIndx,:,i,j,k), temp)
        logT = log10(temp)

        ! Obtain metallicity of this cell (for metal cooling)
        if (cooling_type.eq.COOL_TABLE_METAL) then
          metal = abs(PRIM(bIndx,metalpas,i,j,k)/PRIM(bIndx,1,i,j,k))
        end if

        ! Cooling not applied below cool_Tmin
        if (logT.ge.cool_Tmin) then

          ! Code units
          ETH = CV * PRIM(bIndx,5,i,j,k)
          vel2 = PRIM(bIndx,2,i,j,k)**2 + &
                 PRIM(bIndx,3,i,j,k)**2 + &
                 PRIM(bIndx,4,i,j,k)**2
          EK = 0.5 * PRIM(bIndx,1,i,j,k) * vel2

          ! Interpolate cooling coefficient from table
          if (cooling_type.eq.COOL_TABLE) then
            call find_coef (logT, aloss)
          else if (cooling_type.eq.COOL_TABLE_METAL) then
            call find_coef_metal (logT, metal, aloss)
          end if

          ! Calculate radiative loss and new temperature
          numdens = PRIM(bIndx,1,i,j,k)*d_sc/(mui*AMU)  ! cgs, gas ionized
          ce = aloss*(numdens**2)/(ETH*e_sc)  ! cgs
          cool_factor = exp(-ce*(dt*t_sc))
          frac_loss = 1.0-cool_factor

          ! DEBUG
!          write(logu,*) localBlocks(bIndx), i, j, k
!          write(logu,*) log10(temp), aloss
!          write(logu,*) numdens, frac_loss

          ! Record maximum cooling for this block before limiting
          maxloss = max(maxloss, frac_loss)

          ! Limit cool_factor directly, if needed
          if (cool_factor.lt.1.0-cooling_limit) then
            cool_factor = 1.0-cooling_limit
          end if

          ! Impose a temperature floor by adjusting cool_factor, if needed
          ! Set to 10 K by default
          new_temp = temp * cool_factor
          if (new_temp.lt.T_floor) then
            new_temp = T_floor
            cool_factor = T_floor / temp
          end if

          ! Update pressure and total energy          
          PRIM(bIndx,5,i,j,k) = PRIM(bIndx,5,i,j,k) * cool_factor
          ETH = CV * PRIM(bIndx,5,i,j,k)
          U(bIndx,5,i,j,k) = EK + ETH

        end if

      end do
    end do
  end do

end subroutine apply_cooling

!===============================================================================

!> @brief Interpolates a cooling cofficient from 1D temperature-tabulated data
!> @param logT The base-10 logarithm of the temperature (Kelvin)
!> @param coef Interpolated cooling coefficient
subroutine find_coef (logT, coef)

    use globals, only: cooltable, nptsT, cool_Tmin, cool_Tmax
    implicit none

    real, intent(in) :: logT
    real, intent(out) :: coef

    integer :: i
    real :: T0, T1, C0, C1

    ! Below Tmin, the returned coefficient is zero (no cooling)
    if (logT.lt.cool_Tmin) then

      coef = 0.0
      return

    ! Above Tmax, we assume the coefficient is proportional to sqrt(T)
    ! (corresponding to the free-free regime)
    else if (logT.ge.cool_Tmax) then
    
!      coef = cooltable(2,nptsT)*sqrt(real(T)/cool_Tmax)
      coef = 10**( cooltable(2,nptsT) + 0.5*(logT-cool_Tmax) )
      return

    ! For T in range of the table, we do linear interpolation    
    else

      do i=2,nptsT
        if (cooltable(1,i).gt.logT) then
          T0 = cooltable(1,i-1)
          C0 = cooltable(2,i-1)
          T1 = cooltable(1,i)
          C1 = cooltable(2,i)
          coef = 10**( (C1-C0)/(T1-T0)*(logT-T0) + C0 )
          return
        end if
      end do
      
    end if

  end subroutine find_coef

!===============================================================================

!> @brief Interpolates a cooling cofficient from 2D temperature-metallicity 
!! tabulated data
!> @param logT Base-10 logarithm of temperature (Kelvin)
!> @param Z Metallicity (fraction of metals)
!> @param coef Interpolated cooling coefficient
!> @details The subroutines assumes cooltable holds 2D tabulated coefficients
!! in the form cooltable(T,Z), where each dimension holds nptsT+1 and nptsZ+1
!! values. The first row and the first column state the temperature/metallicity
!! values, while the inner rows and columns hold the corresponding coefficients.
!! Below Tmin, cooling is set to zero. Above Tmax, the linear interpolation
!! with respect to temperature is changed to an extrapolation using a ~sqrt(T)
!! law. Extrapolation over metallicity is linear on both ends.
subroutine find_coef_metal (logT, Z, coef)

    use globals, only: cooltable, nptsT, nptsZ, cool_Tmin, cool_Tmax, cool_Zmin, cool_Zmax
    implicit none

    real, intent(in) :: logT
    real, intent(in) :: Z
    real, intent(out) :: coef

    integer :: Tlo, Thi, Tmid, Zlo, Zhi, Zmid
    real :: T0, T1, Z0, Z1, C0, C1, TC, ZC

    ! Below Tmin, the returned coefficient is zero (no cooling)
    if (logT.lt.cool_Tmin) then
      coef = 0.0
      return
    end if
    
    ! Find the T and Z brackets in which the given values are located
    if (logT.le.cool_Tmax) then
      ! Binary search
      Tlo = 2
      Thi = nptsT+1
      do while ((Thi-Tlo).gt.1)
        Tmid = (Tlo+Thi)/2
        if (logT <= cooltable(Tmid,1)) then
          Thi = Tmid
        else
          Tlo = Tmid
        end if
      end do
      T0 = cooltable(Tlo,1)
      T1 = cooltable(Thi,1)
    end if

    if (Z.lt.cool_Zmin) then
      Zlo = 2
      Zhi = 3
      Z0 = cooltable(1,Zlo)
      Z1 = cooltable(1,Zhi)
    else if (Z.gt.cool_Zmax) then
      Zlo = nptsZ
      Zhi = nptsZ+1
      Z0 = cooltable(1,Zlo)
      Z1 = cooltable(1,Zhi)
    else
      ! Binary search
      Zlo = 2
      Zhi = nptsZ+1
      do while ((Zhi-Zlo).gt.1)
        Zmid = (Zlo+Zhi)/2
        if (Z <= cooltable(1,Zmid)) then
          Zhi = Zmid
        else
          Zlo = Zmid
        end if
      end do
      Z0 = cooltable(1,Zlo)
      Z1 = cooltable(1,Zhi)
    end if
    
    ! Above Tmax, we assume the coefficient is proportional to sqrt(T),
    ! but still linearly interpolate/extrapolate over Z
    if (logT.ge.cool_Tmax) then
    
      C0 = cooltable(nptsT,Zlo) + 0.5*(logT-cool_Tmax)
      C1 = cooltable(nptsT,Zhi) + 0.5*(logT-cool_Tmax)
      coef = 10**( (C1-C0)/(Z1-Z0)*(Z-Z0) + C0 )
      return 
    
    ! For T in range, we do bilinear interpolation over T and Z
    else

      TC = (logT-T0)/(T1-T0)
      ZC = (Z-Z0)/(Z1-Z0)
      coef = (1.0-TC)*(1.0-ZC)*cooltable(Tlo,Zlo) + &
             (1.0-TC)*ZC*cooltable(Tlo,Zhi) + &
             TC*(1.0-ZC)*cooltable(Thi,Zlo) + &
             TC*ZC*cooltable(Thi,Zhi)
      coef = 10**coef
      return
      
    end if

  end subroutine find_coef_metal

