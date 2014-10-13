!===============================================================================
!> @file orbits.f90
!> @brief Computes orbits
!> @author Juan C. Toledo
!> @date 25/Apr/2013

! Copyright (c) 2014 Juan C. Toledo
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

!> @brief A module for computing orbits of bodies
!> @details At the moment contains code to simulate a binary system, 
!! giving the position of the two bodies as a function of time.
!! This implies solving the Kepler problem.

module orbits

  implicit none

  ! Some constants
  real, parameter :: YR   = 31556736.0
  real, parameter :: AU   = 1.4959787e13
  real, parameter :: MSUN = 1.9891e33
  real, parameter :: MEARTH = 5.97e27
  real, parameter :: KPS  = 1.0e5
  real, parameter :: GRAV = 6.67259e-8
  real, parameter :: PI   = acos(-1.0)

  ! ============================================================================
  ! USER-MODIFIABLE section
  ! ============================================================================
  ! This section defines a binary system, in the center of momentum frame
  ! (i.e., the barycentric velocity is assumed to be zero)

  ! Orbital parameters of the binary system
  ! All quantities in cgs
  real, parameter :: M1  = 4.4e6 * MSUN   ! Mass of primary
  real, parameter :: M2  = 3.0 * MEARTH   ! Mass of secondary
  real, parameter :: P   = 198 * YR       ! Period
  real, parameter :: ecc = 0.99           ! Eccentricity
  real, parameter :: Rx  = 0.0 * AU       ! Barycenter x coord
  real, parameter :: Ry  = 0.0 * AU       ! Barycenter y coord

contains

  ! ============================================
  !> @brief Computes the positions and velocities of the bodies 
  !!! at the given time
  !> @details Solves the reduced Kepler problema and returns the
  !! positions and velocities of the two bodies. Note that the orbital
  !! phase is defined as zero at periapsis and 0.5 at apoapsis. The
  !! phase is easily calculated from the time as:
  !!   phase = time/orbital_period + initial_phase
  !> @params phase The orbital phase at which the solution is to be computed
  !> @params x1,y1,x2,y2 Positions of the two bodies
  !> @params vx1,vy1,vx2,vy2 Velocities of the two bodies
  subroutine computeBinary (phase, x1, y1, x2, y2, vx1, vy1, vx2, vy2)

    implicit none

    real, intent(in) :: phase
    real, intent(out) :: x1, y1, x2, y2, vx1, vy1, vx2, vy2

    real, parameter :: PI = acos(-1.0)

    real :: mu, a, L, M, E
    real :: rad, theta, thetadot, raddot
    real :: radx, rady, vx, vy

    ! 1) Solve the reduced problem

    ! Compute reduced mass of the system, semi-major axis of reduced orbit
    ! and total angular momentum of the reduced mass
    mu = (M1*M2)/(M1+M2)
    a = (P*P*GRAV*(M1+M2)/(4*PI*PI))**(1.0/3.0)
    L = mu*(1-ecc)*a*sqrt((1+ecc)/(1-ecc)*GRAV*(M1+M2)/a)

    ! Mean anomaly M
    M = 2*PI*phase

    ! Obtain eccentric anomaly E by numerically solving Kepler's equation
    call solveKepler (M, ecc, E)

    ! Compute radial position and angular position (true anomaly) of reduced
    ! mass
    theta = 2.0*atan(sqrt((1+ecc)/(1-ecc))*tan(E/2))
    rad = a*(1.0-ecc*ecc)/(1+ecc*cos(theta))

    ! Compute radial and angular velocities of reduced mass
    thetadot = L/(mu*rad**2)
    raddot = a*ecc*sin(theta)*(1-ecc**2)/(1+ecc*cos(theta))**2 * thetadot

    ! 2) Convert back to two-body problem

    ! Compute position of bodies
    radx = rad*cos(theta)
    rady = rad*sin(theta)
    x1 = Rx - m2/(m1+m2)*radx
    y1 = Ry - m2/(m1+m2)*rady
    x2 = Rx + m1/(m1+m2)*radx
    y2 = Ry + m1/(m1+m2)*rady

    ! Compute velocities of bodies
    ! (in the center of momentum frame of reference)
    vx = raddot*cos(theta)-rad*thetadot*sin(theta)
    vy = raddot*sin(theta)+rad*thetadot*cos(theta)
    vx1 = -m2/(m1+m2)*vx
    vy1 = -m2/(m1+m2)*vy
    vx2 = m1/(m1+m2)*vx
    vy2 = m1/(m1+m2)*vy

!    print*, phase, x1/AU, y1/AU, sqrt(x1**2+y1**2)/AU, &
!            vx1/KPS, vy1/KPS, sqrt(vx1**2+vy1**2)/KPS, &
!            x2/AU, y2/AU, sqrt(x2**2+y2**2)/AU, &
!            vx2/KPS, vy2/KPS, sqrt(vx2**2+vy2**2)/KPS

  end subroutine computeBinary

  ! ============================================
  !> @brief Solves the Kepler equation, returning the eccentric anomaly
  !> @details Solves the Kepler equation, viz.,
  !!   M = E - e*sin(E)
  !! for E (eccentric anomaly), given M (mean anomaly). This uses a
  !! Newton-Raphson solver with the specified relative tolerance
  !> @params M Mean anomaly
  !> @params E Eccentric anomaly

  subroutine solveKepler (M, ecc, E)

    implicit none
  
    real, intent(in) :: M
    real, intent(in) :: ecc
    real, intent(out) :: E

   real, parameter :: tol = 1e-10
   real :: x, xn, rel_err

    if (M.eq.0.0) then
      E = 0.0
      return
    end if

    if (ecc>0.8) then
      x = sign(PI,M)
    else
      x = M
    end if

    rel_err = tol*2.0
    do while(rel_err.gt.tol)
      xn = x - (x-ecc*sin(x)-M)/(1.0-ecc*cos(x))
      if (x.ne.0.0) rel_err = abs((xn-x)/x)
      x = xn
    end do

    E = x
    return

  end subroutine solveKepler

end module orbits
