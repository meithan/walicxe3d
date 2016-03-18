!===============================================================================
!> @file godunov.f90
!> @brief Generic wrapper for Godunov schemes (HLL1, HLL, HLLC)
!> @author Juan C. Toledo
!> @date 7/Mar/2012

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

!> @brief Wrapper for Godunov-type numerical solvers
!> @details Computes intercell numerical fluxes using any Riemann solver and
!! advances the flow variables one timestep. The solver can be set to be
!! first or second order (in space and time). This is a high-level wrapper.
subroutine Godunov (order)

  use parameters
  use globals
  use tictoc
  implicit none

  integer, intent(in) :: order

  integer :: bIndx, bID, mark, bcount
  real :: dtp

  ! -----------------------------------
  ! 1st-order timestep

  ! Halve timestep for first part of second-order HLL
  if (order.eq.2) then
    dtp = dt/2.0
  else 
    dtp = dt
  end if

  ! Exchange 1-deep boundaries to fill ghost cells (all blocks)
  call tic(mark)
  write(logu,*) ""
  write(logu,'(1x,a)') "> Exchanging 1-deep boundary layers ..."
  call boundary (1, U)
  write(logu,*) "Boundaries exchanged in", nicetoc(mark)

  ! Update primitives in ghost cells
  call calcPrimsAll (U, PRIM, CELLS_GHOST)

  ! First-order integration of all local blocks
  bcount = 0
  call tic(mark)
  write(logu,*) ""
  if (order.eq.1) then
    write(logu,'(1x,a)') "> Integrating blocks ..."
  else if (order.eq.2) then
    write(logu,'(1x,a)') "> Integrating blocks (1st order half step) ..."
  end if

  do bIndx=1,nbMaxProc
    bID = localBlocks(bIndx)
    if (bID.ne.-1) then
    
      ! Compute numerical fluxes  for this block
      select case (solver_type)

        case (SOLVER_HLL1)
          call HLLfluxes (bIndx, 1)

        case (SOLVER_HLL)
          call HLLfluxes (bIndx, 1)

        case (SOLVER_HLLC)
          call HLLCfluxes (bIndx, 1)

      end select

      ! Apply conservative formula
      call upwindStep (bIndx, dtp)
      
      bcount = bcount + 1      

    end if
  end do

  write(logu,'(1x,a,i0,a,a)') "Integrated ", bcount, " blocks in ", nicetoc(mark)

  ! -----------------------------------
  ! 2nd-order full timestep (skipped in 1st-order schemes)

  if (order.eq.2) then

    ! Exchange 2-deep boundaries to fill ghost cells (all blocks)
    call tic(mark)
    write(logu,*) ""
    write(logu,'(1x,a)') "> Exchanging 2-deep boundary layers ..."
    call boundary (2, UP)
    write(logu,*) "Boundaries exchanged in", nicetoc(mark)

    ! Update primitives in ghost cells
    call calcPrimsAll (UP, PRIM, CELLS_ALL)

    ! Second-order integration of all local blocks
    bcount = 0
    call tic(mark)
    write(logu,*) ""
    write(logu,'(1x,a)') "> Integrating blocks (2nd order full step) ..."
    do bIndx=1,nbMaxProc

      bID = localBlocks(bIndx)
      if (bID.ne.-1) then

        ! Compute numerical fluxes for this block
        select case (solver_type)

          case (SOLVER_HLL)
            call HLLfluxes (bIndx, 2)

          case (SOLVER_HLLC)
            call HLLCfluxes (bIndx, 2)

        end select

        ! Apply conservative formula
        call upwindStep (bIndx, dt)

        ! Apply numerical viscosity
        call viscosity (bIndx, U, UP)

        bcount = bcount + 1      

      end if
    end do

    write(logu,'(1x,a,i0,a,a)') "Integrated ", bcount, " blocks in ", &
    nicetoc(mark)
  
  end if

end subroutine Godunov

!===============================================================================

!> @brief Applies conservative formula to calculate stepped U variables
!> @details This is the standard upwing Godunov step. Assumes the numerical
!! intercell fluxes FC have been calculated for this block
!> @param locIndx Block's local index
!> @param dtp Timestep to use
subroutine upwindStep (locIndx, dtp)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: locIndx
  real, intent(in) :: dtp
  
  integer :: i, j, k, bID, lev
  integer :: ieq   ! DEBUG
  real :: dtdx, dtdy, dtdz

  bID = localBlocks(locIndx)
  call meshlevel (bID, lev)

  dtdx = dtp/dx(lev)
  dtdy = dtp/dy(lev)
  dtdz = dtp/dz(lev)

  ! Apply conservative formula (10.2 of Toro) to obtain UPs
  do i=1,ncells_x
    do j=1,ncells_y
      do k=1,ncells_z
! DEBUG
! Check FC, GC and HC values for NaNs
do ieq=1,neqtot
if (FC(ieq,i-1,j,k).ne.FC(ieq,i-1,j,k)) then
  write(logu,'(a,i2,i2,i2)') "Nan in FC at ", i-1, j, k
end if
if (FC(ieq,i,j,k).ne.FC(ieq,i,j,k)) then
  write(logu,'(a,i2,i2,i2)') "Nan in FC at ", i, j, k
end if
if (GC(ieq,i,j-1,k).ne.GC(ieq,i,j-1,k)) then
  write(logu,'(a,i2,i2,i2)') "Nan in GC at ", i, j-1, k
end if
if (GC(ieq,i,j,k).ne.GC(ieq,i,j,k)) then
  write(logu,'(a,i2,i2,i2)') "Nan in GC at ", i, j, k
end if
if (HC(ieq,i,j,k-1).ne.HC(ieq,i,j,k-1)) then
  write(logu,'(a,i2,i2,i2)') "Nan in HC at ", i, j, k-1
end if
if (HC(ieq,i,j,k).ne.HC(ieq,i,j,k)) then
  write(logu,'(a,i2,i2,i2)') "Nan in HC at ", i, j, k
end if
end do
! DEBUG
        UP(locIndx,:,i,j,k)   &
          = U(locIndx,:,i,j,k)   &
          + dtdx*(FC(:,i-1,j,k)-FC(:,i,j,k))   &
          + dtdy*(GC(:,i,j-1,k)-GC(:,i,j,k))   &
          + dtdz*(HC(:,i,j,k-1)-HC(:,i,j,k))
      end do
    end do
  end do

end subroutine upwindStep

!===============================================================================

!> @brief Obtains wavespeed estimates given primitives at interface
!> @details Computes wavespeed estimates for the two outer shock waves and the
!! contact discontinuity (if needed) in the Riemann fan. This routine employs
!! the x-component of speed.
!> @param dl Density at "left" of interface
!> @param dr Density at "right" of interface
!> @param ul Velocity at "left" of interface
!> @param ur Velocity at "right" of interface
!> @param pl Pressure at "left" of interface
!> @param pr Pressure at "right" of interface
!> @param sl Wavespeed estimate for "left" moving wave
!> @param sr Wavespeed estimate for "right" moving wave
!> @param sst Wavespeed estimate for contact discontinuity
subroutine wavespeed (primL, primR, sl, sr, sst)

  use parameters
  implicit none

  real, intent(in) :: primL(neqtot)
  real, intent(in) :: primR(neqtot)
  real, intent(out) :: sl
  real, intent(out) :: sr
  real, intent(out) :: sst

  real :: dl1, ul1, pl1, cl
  real :: dr1, ur1, pr1, cr

  ! Unpack primitives
  dl1 = primL(1)
  dr1 = primR(1)
  ul1 = primL(2)
  ur1 = primR(2)
  pl1 = primL(5)
  pr1 = primR(5)
  
  ! Sound speeds
  call sound (primL,cl)
  call sound (primR,cr)

  ! Compute SL and SR
  ! Davis direct bounded, 10.38 of Toro
  sl = min( ul1-cl, ur1-cr )
  sr = max( ul1+cl, ur1+cr )

  ! Compute S* (HLLC only)
  if (solver_type.eq.SOLVER_HLLC) then
    sst = (pr1-pl1 + dl1*ul1*(sl-ul1) - dr1*ur1*(sr-ur1)) / &
          (dl1*(sl-ul1)-dr1*(sr-ur1))
  else
    sst = 0.0
  end if

  return

end subroutine wavespeed

!===============================================================================

!> @brief Applies the flux limiter to average left/right states
!> @details This routine uses the X-component of speed.
!> @param pll Pressure two cells "left" of interface
!> @param pl Pressure one cell "left" of interface
!> @param prr Pressure two cells "right" of interface
!> @param pr Pressure once cell "right" of interface
!> @param lim Limiter to use. See parameters.f90 for supported options.
!> @param neqs Number of equations
subroutine limiter (pll,pl,pr,prr,lim,neqs)

  implicit none

  integer, intent(in) :: neqs
  integer, intent(in) :: lim
  real, intent(in) :: pll(neqs)
  real, intent(in) :: prr(neqs)
  real, intent(inout) :: pl(neqs)
  real, intent(inout) :: pr(neqs)


  real :: dl, dm, dr, al, ar
  integer :: ieq

  do ieq=1,neqs
    dl = pl(ieq) - pll(ieq)
    dm = pr(ieq) - pl(ieq)
    dr = prr(ieq) - pr(ieq)
    al = average(dl, dm, lim)
    ar = average(dm, dr, lim)
    pl(ieq) = pl(ieq) + 0.5*al
    pr(ieq) = pr(ieq) - 0.5*ar
  end do

contains

  real function average (a,b,opt)

    use constants
    implicit none
    
    real, intent(in) :: a, b
    integer, intent(in) :: opt

    real :: s, c, d, eps

    select case (opt)
    
    case (LIMITER_NONE)
      average = 0.5*(a+b)

    case (LIMITER_VANLEER)
      if (a*b.le.0.0) then
        average = 0.0
      else
        average = a*b*(a+b)/(a*a+b*b)
      end if

    case (LIMITER_MINMOD)
      s = sign(1.0,a)
      average = s*max(0.0, min(abs(a), s*b))

    case (LIMITER_ALBADA)   ! NOT WORKING
      eps = 1.0e-7
      average = (a*(b*b+eps)+b*(a*a+eps))/(a*a+b*b*eps)

    case (LIMITER_UMIST)
      s = sign(1.0,a)
      c = 0.25*a + 0.75*b
      d = 0.75*a + 0.25*b
      average = min(2.0*abs(a), 2.0*s*b, s*c, s*d)
      average = s*max(0.0, average)

    case (LIMITER_WOODWARD)
      s = sign(1.0,a)
      c = 0.5*(a+b)
      average = min(2.0*abs(a), 2*s*b, s*c)
      average = s*max(0.0, average)

    case (LIMITER_SUPERBEE)
      s = sign(1.0,b)
      c = min(2.0*abs(b), s*a)
      d = min(abs(b),2.0*s*a)
      average = s*max(0.0,c,d)

    case default
      average = 0.0
      write(*,'(a)') "WARNING: no averaging in limiter!"
      write(*,'(a,i2)') "Passed limiter value: ", opt 

    end select

  end function average

end subroutine limiter

!===============================================================================
