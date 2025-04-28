!===============================================================================
!> @file windsource.f90
!> @brief Boundary condition module: spherical wind source
!> @author Alejandro Esquivel and Juan C. Toledo
!> @date 17/Apr/2013

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

!> @brief Imposes a spherical wind source
!> @details This user-defined module establishes the parameters, 
!! variables and subroutines required to impose a spherical wind source.
!!
!! How to use this module.
!! The module currently supports
!! 1) The parameters of each wind must be defined by creating an object of
!!    the appropriate type. The currently supported wind types are:
!!      SphericalWindType : a spherical r^-2 wind
!!      PlaneWindType : a planar wind (to be used as in inflow condition)
!!    Then, the user must specify the corresponding wind parameters.
!! 2) The user must modify the userBC() subroutine in the user.f90 source
!!    file, adding a call to the corresponding subroutine: 

module winds

  implicit none

  !===============================
  ! Spherical wind source parameters
  !
  ! ALL VALUES IN CGS
  ! xc, yc, zc: position of the source center
  ! vx, vy, vz: systematic velocity components of the wind source
  ! radius: radius of the wind source
  ! mdot: mass loss rate
  ! vinf: terminal speed
  ! temp: temperature
  ! mu: gas mean atomic mass (in amus)
  ! -- Optional fields --
  ! bx, by, bz: magnetic field components (for runs with passive B field)
  ! metal: metallicity of the gas (for runs with metallicity-dependent cooling)
  ! -- Internal fields (not set by the user) --
  ! bbox: the physical coordinates of the wind source's bounding box
  type SphericalWindType
    real :: xc, yc, zc
    real :: vx, vy, vz
    real :: radius
    real :: mdot
    real :: vinf
    real :: temp
    real :: mu
    real :: bx = 0.0
    real :: by = 0.0
    real :: bz = 0.0
    real :: metal = 1.0
    real :: bbox(6) = 0.0
    logical :: bbox_init = .false.
  end type SphericalWindType
  !===============================

  !===============================
  ! Plane wind source definition
  !
  ! ALL VALUES IN CGS
  ! plane: the plane upon which the wind enters the simulation.
  ! Must be one of the following constants:
  !   PLANE_LEFT: the YZ plane crossing at X=0
  !   PLANE_RIGHT: the YZ plane crossing at X=size_x
  !   PLANE_FRONT: the XZ plane crossing at Y=0
  !   PLANE_BACK: the XZ plane crossing at Y=size_y
  !   PLANE_BOTTOM: the XY plane crossing at Z=0
  !   PLANE_TOP: the XY plane crossing at Z=size_z
  ! rho: density of the flow to be imposed on the plane
  ! vel: velocity *magnitude* of the flow (the sign will be calculated
  !      based on the plane of entry)
  ! temps: temperature of the flow
  ! mu: gas mean atomic mass (in amus)
  ! -- Optional fields --
  ! bx, by, bz: magnetic field components (for runs with passive B field)
  ! metal: metallicity of the gas (for runs with metallicity-dependent cooling)
  integer, parameter :: PLANE_LEFT   = 1
  integer, parameter :: PLANE_RIGHT  = 2
  integer, parameter :: PLANE_FRONT  = 3
  integer, parameter :: PLANE_BACK   = 4
  integer, parameter :: PLANE_BOTTOM = 5
  integer, parameter :: PLANE_TOP    = 6
  type PlaneWindType
    integer :: plane
    real :: rho
    real :: vel
    real :: temp
    real :: mu
    real :: bx = 0.0
    real :: by = 0.0
    real :: bz = 0.0
    real :: metal = 1.0
  end type PlaneWindType
  !===============================

contains
  
  ! ============================================
  !> @brief Imposes a spherical wind source in a specific simulation block
  !> @details Simulates a spherical wind source by setting the flow
  !! conditions in a region of given radius centered at (xc,yc,zc), 
  !! in which a steady-state wind solution is imposed with the given
  !! mass-loss rate, terminal speed and temperature. This is called by
  !! imposeSphericalWindsList().
  !> @param nb The index of the block in the local blocks registry
  !> @param wind A SphericalWindType variable containing the
  !! wind parameters. See the module documentation for further details.
  !> @param uvars Flow variables array to be modified
  subroutine imposeSphericalWindBlock (nb, wind, uvars)
 
    use parameters
    use globals
    implicit none

    integer, intent(in) :: nb
    type(SphericalWindType), intent(in) :: wind
    real, intent(inout) :: uvars(nbMaxProc, neqtot, &
                           nxmin:nxmax, nymin:nymax, nzmin:nzmax)

    integer :: i, j, k, bID
    real :: x, y, z, dist, dens, vx, vy, vz, pres
    real :: primit(neqtot)

    bID = localblocks(nb)

    ! Loop over all cells of the block
    do k=nzmin,nzmax
      do j=nymin,nymax
        do i=nxmin,nxmax

          ! Compute absolute cell position, de-scaled
          call cellPos (bID, i, j, k, x, y, z)
          x = x*l_sc;  y = y*l_sc;  z = z*l_sc

          ! Distance to wind source center
          dist = sqrt((x-wind%xc)**2+(y-wind%yc)**2+(z-wind%zc)**2)   ! phys units
          
          if (dist.lt.wind%radius) then

            ! Calculate wind density, velocity components and pressure in this cell
            dens = wind%mdot/wind%vinf/dist/dist/(4.0*PI)
            vx = wind%vinf*(x-wind%xc)/dist + wind%vx
            vy = wind%vinf*(y-wind%yc)/dist + wind%vy
            vz = wind%vinf*(z-wind%zc)/dist + wind%vz
            pres = dens/(wind%mu*AMU)*KB*wind%temp

            ! Scale primitives
            primit(1) = dens/d_sc
            primit(2) = vx/v_sc
            primit(3) = vy/v_sc
            primit(4) = vz/v_sc
            primit(5) = pres/p_sc
#ifdef PASBP
            ! Magnetic field
            primit(6) = wind%bx
            primit(7) = wind%by
            primit(8) = wind%bz
#endif
            ! Passive scalar for metalicity
            if (cooling_type.eq.COOL_TABLE_METAL) then
              primit(metalpas) = wind%metal*primit(1)
            end if

            ! Convert primitives and set flow vars for this cell
            call prim2flow( primit, uvars(nb,:,i,j,k) )

          end if

        end do ! i
      end do ! j
    end do ! k 

  end subroutine imposeSphericalWindBlock
  
  ! ============================================
  !> @brief Imposes a list of spherical wind sources on the simulation
  !> @details This is the main subroutine to call to impose multiple
  !! wind sources, as it is more efficient than calling imposeSphericalWind()
  !! for each wind in turn.
  !> @num_winds The number of wind sources contained in the array
  !> @winds_list An array of SphericalWindType variables with the properties
  !! of the winds to be imposed
  !> @param uvars Flow variables array to be modified
  subroutine imposeSphericalWindsList (num_winds, winds_list, uvars)
 
    use parameters
    use globals
    implicit none

    integer, intent(in) :: num_winds
    type(SphericalWindType), intent(inout) :: winds_list(:)
    real, intent(inout) :: uvars(nbMaxProc, neqtot, &
                           nxmin:nxmax, nymin:nymax, nzmin:nzmax)

    integer :: i, l, nb, bID
    type(SphericalWindType) :: wind
    logical :: intersects
    real :: bbox(6)

    write(logu,*) ""
    write(logu,'(1x,a)') "> Imposing spherical wind sources ..."
    write(logu,'(1x,a,i,a)') "There are ", num_winds, " wind sources"

    ! Only once: refine zones, compute bounding box and report
    ! parameters for all wind sources
    do i=1,num_winds
      wind = winds_list(i)
      if (.not.wind%bbox_init) then
        winds_list(i)%bbox(1) = wind%xc - wind%radius
        winds_list(i)%bbox(2) = wind%xc + wind%radius
        winds_list(i)%bbox(3) = wind%yc - wind%radius
        winds_list(i)%bbox(4) = wind%yc + wind%radius
        winds_list(i)%bbox(5) = wind%zc - wind%radius
        winds_list(i)%bbox(6) = wind%zc + wind%radius
        winds_list(i)%bbox_init = .true.
      end if
      if (it.eq.0) then
        write(logu,*)
        write(logu,*) "----------------------------"
        write(logu,'(1x,a,i)') "Source #", i
        write(logu,'(1x,a,es12.5,a)') "Mdot = ", wind%mdot / (MSUN/YR), " Msun/yr"
        write(logu,'(1x,a,es12.5,a)') "vinf = ", wind%vinf / KPS, " km/s"
        write(logu,'(1x,a,es12.5,a)') "radius = ", wind%radius / PC, " pc"
        write(logu,'(1x,a,es12.5,a)') "rho(R) = ", wind%mdot/wind%vinf/wind%radius/wind%radius/(4.0*PI), " g/cm^3"
        write(logu,'(1x,a,es12.5,a)') "temp = ", wind%temp, " K"
        write(logu,'(1x,a,es12.5,1x,es12.5,1x,es12.5,a)') "Location: ", wind%xc/PC, wind%yc/PC, wind%zc/PC, " pc"
        write(logu,'(1x,a,f8.5,1x,f8.5,1x,f8.5,1x,f8.5,1x,f8.5,1x,f8.5,a)') "Bbox: ", winds_list(i)%bbox(1)/PC, winds_list(i)%bbox(2)/PC, winds_list(i)%bbox(3)/PC, winds_list(i)%bbox(4)/PC, winds_list(i)%bbox(5)/PC, winds_list(i)%bbox(6)/PC, " pc"
        call refineZone (winds_list(i)%bbox, maxlev)
      end if
    end do

    ! Impose flow conditions, where applicable
    do nb=1,nbMaxProc
      bID = localBlocks(nb)
      if (bID.ne.-1) then

        ! Iterate over wind sources
        do l=1,num_winds
           
          wind = winds_list(l)

          ! Check whether the wind bounding box intersects with the block
          call getBoundingBox(bID, bbox)
          bbox(1) = bbox(1) * l_sc
          bbox(2) = bbox(2) * l_sc
          bbox(3) = bbox(3) * l_sc
          bbox(4) = bbox(4) * l_sc
          bbox(5) = bbox(5) * l_sc
          bbox(6) = bbox(6) * l_sc
          intersects = ((wind%bbox(1).le.bbox(2)) .and. (wind%bbox(2).ge.bbox(1)) &
          .and. (wind%bbox(3).le.bbox(4)) .and. (wind%bbox(4).ge.bbox(3)) &
          .and. (wind%bbox(5).le.bbox(6)) .and. (wind%bbox(6).ge.bbox(5)))

          ! Only impose wind if the wind bbox intersects
          if (intersects) then 
            call imposeSphericalWindBlock(nb, winds_list(l), uvars)
          end if

        end do

      end if
    end do

  end subroutine imposeSphericalWindsList

! =============================================
  
  !> @brief Impose a single spherical wind source in the simulation
  !> @details This is a wrapper that will create a wind_list with a single
  !! wind and call imposeSphericalWindsList() on it
  !> @params wind A SphericalWindType variable containing the wind
  !! parameters. See the module documentation for further details.
  !> @param uvars Flow variables array to be modified
  subroutine imposeSphericalWind (wind, uvars)
    
    use parameters
    use globals
    implicit none

    type(SphericalWindType), intent(inout) :: wind
    real, intent(inout) :: uvars(nbMaxProc, neqtot, &
                           nxmin:nxmax, nymin:nymax, nzmin:nzmax)

    type(SphericalWindType), dimension(1) :: winds_list_temp

    ! Compute wind bbox if needed
    if (.not.wind%bbox_init) then
      wind%bbox(1) = wind%xc - wind%radius
      wind%bbox(2) = wind%xc + wind%radius
      wind%bbox(3) = wind%yc - wind%radius
      wind%bbox(4) = wind%yc + wind%radius
      wind%bbox(5) = wind%zc - wind%radius
      wind%bbox(6) = wind%zc + wind%radius
      wind%bbox_init = .true.
    end if

    winds_list_temp(1) = wind
    call imposeSphericalWindsList(1, winds_list_temp, uvars)

  end subroutine imposeSphericalWind

! =============================================

  !> @brief Imposes a planar wind on one of the edges of the simulation box.
  !> @details Simulates a planar wind entering the simulation box
  !! by setting the flow conditions in the appropriate boundary layer. 
  !> @params wind_params A PlaneWindType variable containing the wind
  !! parameters. See the module documentation for further details.
  !> @param uvars Flow variables array to be modified
  subroutine imposePlaneWind (wind_params, uvars)
 
    use parameters
    use globals
    implicit none

    type(PlaneWindType), intent(in) :: wind_params
    real, intent(inout) :: uvars(nbMaxProc, neqtot, &
                           nxmin:nxmax, nymin:nymax, nzmin:nzmax)

    integer :: nb, bID, i, j, k, plane
    integer :: i1, i2, j1, j2, k1, k2, dir
    real :: dens, vx, vy, vz, pres, temp, vel, mu, metal
    real :: primit(neqtot)
    integer :: neighType
    integer :: neighList(4)

    write(logu,*) ""
    write(logu,'(1x,a)') "> Refreshing plane wind ..."

    ! Unpack wind source parameters
    plane = wind_params%plane
    dens  = wind_params%rho
    vel   = abs(wind_params%vel)
    temp  = wind_params%temp
    mu    = wind_params%mu
    metal = wind_params%metal

    ! Set appropriate variables depending on selected plane boundary
    i1 = nxmin
    i2 = nxmax
    j1 = nymin
    j2 = nymax
    k1 = nzmin
    k2 = nzmax
    select case(plane)
      case (PLANE_LEFT)
        i2 = 0
        dir = LEFT
      case (PLANE_RIGHT)
        i1 = ncells_x + 1
        dir = RIGHT
      case (PLANE_FRONT)
        j2 = 0
        dir = FRONT
      case (PLANE_BACK)
        j1 = ncells_y + 1
        dir = BACK
      case (PLANE_BOTTOM)
        k2 = 0
        dir = BOTTOM
      case (PLANE_TOP)
        k1 = ncells_z + 1
        dir = TOP
    end select

    ! Impose flow conditions on ghost cells
    do nb=1,nbMaxProc
      bID = localBlocks(nb)
      if (bID.ne.-1) then

        ! Determine if this block is at the selected boundary
        call neighbors (bID, dir, neighType, neighList)

        ! If so, modify the ghost cells along that boundary
        if (neighType.eq.NEIGH_BOUNDARY) then

          do i=i1,i2
            do j=j1,j2
              do k=k1,k2

                ! Compute velocity components and pressure
                vx = 0.0
                vy = 0.0
                vz = 0.0
                select case(plane)
                  case (PLANE_LEFT)
                    vx = vel
                  case (PLANE_RIGHT)
                    vx = -vel
                  case (PLANE_FRONT)
                    vy = vel
                  case (PLANE_BACK)
                    vy = -vel
                  case (PLANE_BOTTOM)
                    vz = vel
                  case (PLANE_TOP)
                    vz = -vel
                end select
                pres = dens/(mu*AMU)*KB*temp

                ! Scale primitives
                primit(1) = dens/d_sc
                primit(2) = vx/v_sc
                primit(3) = vy/v_sc
                primit(4) = vz/v_sc
                primit(5) = pres/p_sc

                ! Passive scalar for metalicity
                if (cooling_type.eq.COOL_TABLE_METAL) then
                  primit(metalpas) = metal*primit(1)
                end if

                ! Convert primitives and set flow vars for this cell
                call prim2flow( primit, uvars(nb,:,i,j,k) )

              end do
            end do
          end do

        end if

      end if
    end do

  end subroutine imposePlaneWind


! =============================================

end module winds
