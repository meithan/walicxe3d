!===============================================================================
!> @file hllc.f90
!> @brief HLLC Riemann solver
!> @author Juan C. Toledo
!> @date 22/Feb/2013

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

!> @brief Contact discontinuity-capturing modification of the
!! Harten, Lax, van Leer Approximate Riemann Solver (HLLC)
!> @details This is a modification of the standard HLL Riemann solver, designed
!! to incorporate information about the contact discontinuity. It computes the
!! intercell numerical fluxes for every cell interface in a block using the
!! HLLC solver. This routine assumes the global vector of primitives is
!! up-to-date (including enough ghost cells) and exits with computed values in
!! the intercell numerical fluxes (FC, GC, HC). Interfaces are defined to the
!! right, that is:
!! FC(:,i,j,k): intercell flux between cell (i,j,k) and cell (i+1,j,k)
!! GC(:,i,j,k): intercell flux between cell (i,j,k) and cell (i,j+1,k)
!! HC(:,i,j,k): intercell flux between cell (i,j,k) and cell (i,j,k+1)
!! The order of the solver must be specified. First-order employs piecewise
!! constant primitives in each cell, while second-order uses linear
!! reconstruction to interpolate the primitives at cell interfaces.
!> @param locIndx Local index of the block
!> @param order Order of solver (interpolation of primitives)
subroutine HLLCfluxes (locIndx, order)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: locIndx
  integer, intent(in) :: order
  
  integer :: i, j, k
  real :: pl(neqtot), pr(neqtot), pll(neqtot), prr(neqtot), ff(neqtot)
  logical :: valid

  ! First-order method
  select case (order)

  ! -------------------------------------------

  case (1)   ! First-order
  
    do i=0,ncells_x
      do j=0,ncells_y
        do k=0,ncells_z

          call validCell (i,j,k,valid)
          if (valid) then

! DEBUG
!write(*,*) i,j,k
!"if (rank.eq.master) then
!  write(*,'(a,i3,i3,i3)') "Order 1, cell:", i,j,k
!end if
! DEBUG
            ! X dimension
! DEBUG
!if (rank.eq.master) write(*,*) "X direction"
! DEBUG
            pl(:) = PRIM(locIndx,:,i,j,k)
            pr(:) = PRIM(locIndx,:,i+1,j,k)
            call primfhllc (pL, pR, ff)
            FC(:,i,j,k) = ff(:)

            ! Y dimension
! DEBUG
!if (rank.eq.master) write(*,*) "Y direction"
! DEBUG
            pL(:) = PRIM(locIndx,:,i,j,k)
            pR(:) = PRIM(locIndx,:,i,j+1,k)
            call swapxy (pL)
            call swapxy (pR)
            call primfhllc (pL, pR, ff)
            call swapxy (ff)
            GC(:,i,j,k) = ff(:)

            ! Z dimension
! DEBUG
!if (rank.eq.master) write(*,*) "Z direction"
! DEBUG
            pL(:) = PRIM(locIndx,:,i,j,k)
            pR(:) = PRIM(locIndx,:,i,j,k+1)
            call swapxz (pL)
            call swapxz (pR)
            call primfhllc (pL, pR, ff)
            call swapxz (ff)
            HC(:,i,j,k) = ff(:)
         
          end if
          
        end do
      end do
    end do

  ! -------------------------------------------

  case (2)  ! Second-order - requires limiter

    do i=0,ncells_x
      do j=0,ncells_y
        do k=0,ncells_z

          call validCell (i,j,k,valid)
          if (valid) then

! DEBUG
!if (rank.eq.master) then
!  write(*,'(a,i3,i3,i3)') "Order 2, cell:", i,j,k
!end if
!write(*,*) "PRIMS="
!write(*,*) PRIM(locIndx,:,i,j,k)
! DEBUG
            ! X dimension
! DEBUG
!if (rank.eq.master) write(*,*) "X direction"
! DEBUG
            pll(:) = PRIM(locIndx,:,i-1,j,k)
            pl(:)  = PRIM(locIndx,:,i,  j,k)
            pr(:)  = PRIM(locIndx,:,i+1,j,k)
            prr(:) = PRIM(locIndx,:,i+2,j,k)
            call limiter (pll,pl,pr,prr,limiter_type,neqtot)
            call primfhllc (pl, pr, ff)
            FC(:,i,j,k) = ff(:)

            ! Y dimension
! DEBUG
!if (rank.eq.master) write(*,*) "Y direction"
! DEBUG
            pll(:) = PRIM(locIndx,:,i,j-1,k)
            pl(:)  = PRIM(locIndx,:,i,j  ,k)
            pr(:)  = PRIM(locIndx,:,i,j+1,k)
            prr(:) = PRIM(locIndx,:,i,j+2,k)
            call swapxy (pll)
            call swapxy (pl)
            call swapxy (pr)
            call swapxy (prr)
            call limiter (pll,pl,pr,prr,limiter_type,neqtot)
            call primfhllc (pl, pr, ff)
            call swapxy (ff)
            GC(:,i,j,k) = ff(:)

            ! Z dimension
! DEBUG
!if (rank.eq.master) write(*,*) "Z direction"
! DEBUG
            pll(:) = PRIM(locIndx,:,i,j,k-1)
            pl(:)  = PRIM(locIndx,:,i,j,k  )
            pr(:)  = PRIM(locIndx,:,i,j,k+1)
            prr(:) = PRIM(locIndx,:,i,j,k+2)            
            call swapxz (pll)
            call swapxz (pl)
            call swapxz (pr)
            call swapxz (prr)
            call limiter (pll,pl,pr,prr,limiter_type,neqtot)            
            call primfhllc (pl, pr, ff)
            call swapxz (ff)
            HC(:,i,j,k) = ff(:)
         
          end if
          
        end do
      end do
    end do
    
  end select
  
end subroutine HLLCfluxes

!===============================================================================

!> @brief HLLC intercell fluxes along X using primitives at an interface
!> @param pL Vector of primitives at LEFT of interface
!> @param pR Vector of primitives at RIGHT of interface
!> @param ff Returned vector of HLLC intercell fluxes
subroutine primfhllc (pL, pR, ff)

  use parameters
!  use globals, only: rank
  implicit none
  
  real, intent(in) :: pL(neqtot)
  real, intent(in) :: pR(neqtot)
  real, intent(out) :: ff(neqtot)
  
  real :: sl, sr, sst, rhost, EL, ER
  real :: ustl(neqtot), ustr(neqtot)
  real :: uL(neqtot), uR(neqtot)
  real :: fL(neqtot), fR(neqtot)

  ! Calculate wavespeeds
  call wavespeed (pL, pR, sl, sr, sst)

! DEBUG
!if (rank.eq.master)  write(*,*) "sl,sr=", sl, sr
!if (rank.eq.master)  write(*,*) "sst=", sst
! DEBUG

  ! Obtain HLLC intercell fluxes as given by 10.34 of Toro
  if (sl.gt.0) then

!    write(*,*) "sl>0"
    call prim2fluxes (pL, DIM_X, fL)
    ff(:) = fL(:)
    return

  else if (sr.lt.0) then

!    write(*,*) "sr<0"
    call prim2fluxes (pR, DIM_X, fR)  
    ff(:) = fR(:)
    return

  else if (sst.ge.0) then

!    write(*,*) "sst>0"
    call prim2fluxes (pL, DIM_X, fL)
    call prim2flow (pL, uL)

    ! U-star state given by 10.33 of Toro
    rhost = pL(1)*(sl-pL(2))/(sl-sst)
    EL = 0.5*pL(1)*(pL(2)**2+pL(3)**2+pL(4)**2)+cv*pL(5)
    ustl(1) = rhost
    ustl(2) = rhost*sst
    ustl(3) = rhost*pL(3)
    ustl(4) = rhost*pL(4)
    ustl(5) = rhost*( EL/pL(1)+(sst-pL(2))*(sst+pL(5)/(pL(1)*(sl-pL(2)))) )
    if (npassive.ge.1) then
      ustl(6:neqtot) = rhost*pL(6:neqtot)/pL(1)
    end if

    ff(:) = fL + sl*(ustl(:)-uL(:))
    return
      
  else if (sst.le.0) then

!    write(*,*) "sst<0"
    call prim2fluxes (pR, DIM_X, fR)
    call prim2flow (pR, uR)

    ! U-star state given by 10.33 of Toro
    rhost = pR(1)*(sr-pR(2))/(sr-sst)
    ER = 0.5*pR(1)*(pR(2)**2+pR(3)**2+pR(4)**2)+cv*pR(5)
    ustr(1) = rhost
    ustr(2) = rhost*sst
    ustr(3) = rhost*pR(3)
    ustr(4) = rhost*pR(4)
    ustr(5) = rhost*( ER/pR(1)+(sst-pR(2))*(sst+pR(5)/(pR(1)*(sr-pR(2)))) )
    if (npassive.ge.1) then
      ustr(6:neqtot) = rhost*pR(6:neqtot)/pR(1)
    end if
    
    ff(:) = fR + sr*(ustr(:)-uR(:))
    return
  
  else

!    write(*,'(a)') "Error in primfhllc routine!"
!    write(*,'(a)') "***ABORTING***"
!    call clean_abort(0)
  
  end if

end subroutine primfhllc

!===============================================================================
