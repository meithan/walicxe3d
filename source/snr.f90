!===============================================================================
!> @file snr.f90
!> @brief Boundary condition module: supernova remnant
!> @author Juan C. Toledo
!> @date 2/Feb/2012

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

!> @brief Imposes a supernova remnant
!> @details This module contains code to initialize supernova remnants
!!
!! NOTE: do not modify the code in this module! Import it in your user
!! module and utilize the provided code.
!!
!! How to use this module. To detonate a SNR:
!! 1) Import the 'snr' module in your user module
!! 2) Create an object of type snr_params_type to hold the parameters of
!!    the SNR you want to detonate. The descriptions of the object's fields
!!    are given below. All fields are mandatory.
!! 3) Whenever the SNR is to be detonated, call the detonateSNR() subroutine,
!!    passing as arguments the snrparams object and the flow variables array. 
module snr

  implicit none

  !===============================
  ! Please don't modify this section

  ! This custom type contains the data required by the detonateSNR to
  ! set off a supernova remnant.
  ! The variables have the following meanings:
  !  xc, yc, zc: physical coordinates of the SNR center
  !  radius: initial radius of the SNR
  !  mass: mass of SN ejecta
  !  energy: total energy (thermal+kinetic) of explosion
  !  chi: ratio of kinetic to total energy
  ! -- The following fields are optional --
  !  bx, by, bz: magnetic field strength in ejecta (may be ignored)
  !  time: time at which the SNR will be imposed
  !  metal: metallicity of the ejecta (may be ignored)
  !  armed: a logical flag used to disarm the SNR once it is imposed
  ! All physical values must be given in *CGS*.
  type snr_params_type
    real :: xc
    real :: yc
    real :: zc
    real :: radius
    real :: mass
    real :: energy
    real :: chi
    real :: bx = 0.0
    real :: by = 0.0
    real :: bz = 0.0
    real :: time = 0.0
    real :: metal = 1.0
    logical :: armed = .true.
  end type snr_params_type
 
  ! ============================================================================

contains

  ! ============================================
  !> @brief Simulates the remnant of a supernova explosion (simple ejecta)
  !> @details Simulates the detonation of a supernova by imposing simple
  !! flow properties for the ejecta in a region of radius RSN centered at
  !! (xc,yc,zc), having a total mass of MSN, a total (kinetic+thermal) energy
  !! of ESN, and a ratio of (kinetic energy)/(total energy) given by chi.
  !> @params snr_params A snr_params_type object containing the SNR parameters.
  !! See the snr module main documentation for further details.
  !> @param uvars Full flow variables array to be modified (U or UP)
  subroutine detonateSNR (snr_params, uvars)

    use parameters
    use globals
    implicit none
    type(snr_params_type), intent(inout) :: snr_params
    real, intent(inout) :: uvars(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)

    integer :: nb, bID, i, j, k
    real :: xc, yc, zc, RSN, MSN, ESN, chi, bx, by, bz
    real :: x, y, z, dist
    real(kind=8) :: Ek, Et, Pres, Dens, Vmax, metal
    real(kind=8) :: Vel, Vx, Vy, Vz
    real :: primit(neqtot), flowvars(neqtot)
    real :: zone(6)
    integer :: zlevel

    write(logu,*) ""
    write(logu,'(1x,a)') "> Imposing supernova remnant ..."

    ! Unpack SNR parameters
    xc = snr_params%xc
    yc = snr_params%yc
    zc = snr_params%zc
    RSN = snr_params%radius
    MSN = snr_params%mass
    ESN = snr_params%energy
    chi = snr_params%chi
    bx = snr_params%bx
    by = snr_params%by
    bz = snr_params%bz
    metal = snr_params%metal

    ! Calculate derived parameters in cgs
    Ek = chi*ESN
    Et = (1.0-chi)*ESN
    Dens = MSN/(4.0/3.0*PI*RSN**3)
!    Vmax = sqrt(2.0*Ek/mass)        ! for constant v(r) profile
    Vmax = sqrt(10.0/3.0*Ek/MSN)     ! for linear v(r) profile
    Pres = (1.0/CV)*Et/(4.0/3.0*PI*RSN**3)

    write(logu,'(1x,a)') "> SNR parameters (cgs):"
    write(logu,'(1x,a,es12.5)') "Ekin = ", Ek
    write(logu,'(1x,a,es12.5)') "Eter = ", Et
    write(logu,'(1x,a,es12.5)') "Etot = ", Ek+Et
    write(logu,'(1x,a,es12.5)') "Dens = ", Dens
    write(logu,'(1x,a,es12.5)') "Vmax = ", Vmax
    write(logu,'(1x,a,es12.5)') "Pres = ", Pres
    write(logu,'(1x,a,es12.5,1x,es12.5,1x,es12.5)') "Location: ", xc, yc, zc

    ! A zone around the SNR will be refined to provide adequate resolution
    zone(1) = xc - RSN
    zone(2) = xc + RSN
    zone(3) = yc - RSN
    zone(4) = yc + RSN
    zone(5) = zc - RSN
    zone(6) = zc + RSN
    zlevel = maxlev
    call refineZone (zone, zlevel)

    ! Impose flow conditions
    do nb=1,nbMaxProc
      bID = localBlocks(nb)
      if (bID.ne.-1) then

        do i=nxmin,nxmax
          do j=nymin,nymax
            do k=nzmin,nzmax

              call cellPos (bID, i, j, k, x, y, z)
              x = x*l_sc;  y = y*l_sc;  z = z*l_sc
              dist = sqrt((x-xc)**2+(y-yc)**2+(z-zc)**2)   ! phys units

              ! DEBUG
!              write(logu,*) "CELL coords:", x/PC, y/PC, z/PC
!              write(logu,*) "dist to center of SN=", dist/PC
!              write(logu,*) "RSN=", RSN/PC

              if (dist.lt.RSN) then

                ! DEBUG
!                write(logu,*) "Cell", x, y, z, "in SNR region"

                ! Calculate velocity components
                Vel = (dist/RSN)*Vmax
                Vx = Vel*(x-xc)/dist
                Vy = Vel*(y-yc)/dist
                Vz = Vel*(z-zc)/dist

                ! Scale primitives
                primit(1) = Dens/d_sc
                primit(2) = Vx/v_sc
                primit(3) = Vy/v_sc
                primit(4) = Vz/v_sc
                primit(5) = Pres/p_sc

                ! Magnetic field (if applicable)
#ifdef PASBP
                primit(6) = bx
                primit(7) = by
                primit(8) = bz
#endif
                ! Passive scalar for metalicity
                if (cooling_type.eq.COOL_TABLE_METAL) then
                  primit(metalpas) = metal*primit(1)
                end if

                ! Convert primitives and set flow vars for this cell
                call prim2flow (primit, flowvars)
                uvars(nb,:,i,j,k) = flowvars(:)

              end if

            end do
          end do
        end do

      end if
    end do

    ! Disarm SNR
    snr_params%armed = .false.

  end subroutine detonateSNR

  ! ============================================
  !> @brief Simulates the remnant of a supernova explosion (SN Ia ejecta model)
  !> @details Simulates the detonation of a supernova by imposing flow 
  !! properties following the SN Ia ejecta model of Jun & Norman (1996).
  !! The method is similar to the 'simple' version, except for the following:
  !! Density is taken as:
  !!   rho_c                 if   r < r_c
  !!   rho_c*(r/r_c)^(-n)    if   r > r_c
  !! where rho_c and r_c are the ensity and radius of the core, 
  !! and n (=7) is the core density power law index.
  !! The core radius is given by:
  !!   r_c = R*[1-x*(3-n)*M/(4*pi*rho_0*R^3)]^[1/(3-n)]
  !! while the core density by:
  !!   rho_c = 3*(1-x)*M/(4*pi*r_c^3)
  !! where x (=3/7) is the mass fraction of the core, rho_0 is the density
  !! at the edge of the core, and M and R are the total mass and radius of
  !! the SN ejecta.
  !! The velocity is taken as linearly proportional to radius:
  !!   v(r) = v0*r/R
  !! where v0 is the velocity at the outer boundary and is given by:
  !!   v0 = Ek^(1/2)*{2*pi*rho_c*r_c^5/(5*R^2) + 
  !!       + 2*pi*rho_0*R^3*[1-(R/r_c)^(n-5)]/(5-n)}
  !! Finally, pressure can be obtained from the thermal energy and volume:
  !!   P = (1/CV)*Eth/(4/3*PI*R^3)
  !! where Eth = E-Ek is the thermal energy (this is different from
  !! the method of Jun & Norman).
  !> @params snr_params A snrparams_type object containing the SNR parameters.
  !! See the snr module main documentation for further details.
  !> @param uvars Full flow variables array to be modified (U or UP)
  subroutine detonateSNRIa (snr_params, uvars)

    use parameters
    use globals
    implicit none
    type(snr_params_type), intent(inout) :: snr_params
    real, intent(inout) :: uvars(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)

    integer :: nb, bID, i, j, k
    real :: xc, yc, zc, RSN, MSN, ESN, chi, bx, by, bz, metal
    real :: x, y, z, dist
    real(kind=8) :: n, xm, rc, rhoc, rho0
    real(kind=8) :: Ek, Eth, pres, dens, vmax
    real(kind=8) :: vel, vx, vy, vz
    real :: primit(neqtot), flowvars(neqtot)
    real :: zone(6)
    integer :: zlevel

    write(logu,*) ""
    write(logu,'(1x,a)') "> Imposing Jun & Norman supernova remnant ..."

    ! Unpack SNR parameters
    xc = snr_params%xc
    yc = snr_params%yc
    zc = snr_params%zc
    RSN = snr_params%radius
    MSN = snr_params%mass
    ESN = snr_params%energy
    chi = snr_params%chi
    bx = snr_params%bx
    by = snr_params%by
    bz = snr_params%bz
    metal = snr_params%metal

    ! Ejecta model parameters
    n = 7.0
    xm = 3.0/7.0

    ! Calculate derived parameters in cgs
    Ek = chi*ESN
    Eth = (1.0-chi)*ESN
    rc = RSN*(1.0-xm*(3-n)*MSN/(4*PI*4*ism_dens*RSN**3))**(1.0/(3.0-n))
    rhoc = 3.0*(1.0-xm)*MSN/(4.0*PI*rc**3)
    rho0 = rhoc*(RSN/rc)**(-n)
    vmax = Ek**0.5*(2.0*PI*rhoc*rc**5/(5.0*RSN**2) + &
           2.0*PI*rho0*RSN**3*(1.0-(RSN/rc)**(n-5.0))/(5.0-n))**(-0.5)
    pres = (1.0/CV)*Eth/(4.0/3.0*PI*RSN**3)

    write(logu,'(1x,a)') "> SNR parameters (cgs):"
    write(logu,'(1x,a,es12.5)') "Ek   = ", Ek
    write(logu,'(1x,a,es12.5)') "Eth  = ", Eth
    write(logu,'(1x,a,es12.5)') "Etot = ", Ek+Eth
    write(logu,'(1x,a,es12.5)') "Core density = ", rhoc
    write(logu,'(1x,a,es12.5)') "Core radius  = ", rc
    write(logu,'(1x,a,es12.5)') "Outer density  = ", rho0
    write(logu,'(1x,a,es12.5)') "4*rho_ism      = ", 4*ism_dens
    write(logu,'(1x,a,es12.5)') "Outer velocity = ", vmax
    write(logu,'(1x,a,es12.5)') "Pres = ", pres
    write(logu,'(1x,a,es12.5,1x,es12.5,1x,es12.5)') "Location: ", xc, yc, zc

    ! A zone around the SNR will be refined to provide adequate resolution
    zone(1) = xc - RSN
    zone(2) = xc + RSN
    zone(3) = yc - RSN
    zone(4) = yc + RSN
    zone(5) = zc - RSN
    zone(6) = zc + RSN
    zlevel = maxlev
    call refineZone (zone, zlevel)

    ! Impose flow conditions
    do nb=1,nbMaxProc
      bID = localBlocks(nb)
      if (bID.ne.-1) then

        do i=nxmin,nxmax
          do j=nymin,nymax
            do k=nzmin,nzmax

              call cellPos (bID, i, j, k, x, y, z)
              x = x*l_sc;  y = y*l_sc;  z = z*l_sc
              dist = sqrt((x-xc)**2+(y-yc)**2+(z-zc)**2)   ! phys units

              ! DEBUG
!              write(logu,*) "CELL coords:", x/PC, y/PC, z/PC
!              write(logu,*) "dist to center of SN=", dist/PC
!              write(logu,*) "RSN=", RSN/PC

              if (dist.lt.RSN) then

                ! DEBUG
!                write(logu,*) "Cell", x, y, z, "in SNR region"

                ! Density
                if (dist.lt.rc) then
                  dens = rhoc
                else
                  dens = rhoc*(dist/rc)**(-n)
                end if

                ! Velocity components
                vel = vmax*(dist/RSN)
                vx = vel*(x-xc)/dist
                vy = vel*(y-yc)/dist
                vz = vel*(z-zc)/dist

                ! Scale primitives
                primit(1) = dens/d_sc
                primit(2) = vx/v_sc
                primit(3) = vy/v_sc
                primit(4) = vz/v_sc
                primit(5) = pres/p_sc

                ! Magnetic field (if applicable)
#ifdef PASBP
                primit(6) = bx
                primit(7) = by
                primit(8) = bz
#endif
                ! Passive scalar for metalicity
                if (cooling_type.eq.COOL_TABLE_METAL) then
                  primit(metalpas) = metal*primit(1)
                end if

                ! Convert primitives and set flow vars for this cell
                call prim2flow (primit, flowvars)
                uvars(nb,:,i,j,k) = flowvars(:)

              end if

            end do
          end do
        end do

      end if
    end do

    ! Disarm SNR
    snr_params%armed = .false.

  end subroutine detonateSNRIa

end module snr
