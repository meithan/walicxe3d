!===============================================================================
!> @file hll.f90
!> @brief HLL Riemann solver
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

!> @brief Harten, Lax, van Leer (HLL) Approximate Riemann Solver
!> @details Computes the intercell numerical fluxes for every cell interface
!! in a block using the HLL solver. This routine assumes the global vector of
!! primitives is up-to-date (including enough ghost cells) and exits with
!! computed values in the intercell numerical fluxes (FC, GC, HC). Interfaces
!! are defined to the right, that is:
!! FC(:,i,j,k): intercell flux between cell (i,j,k) and cell (i+1,j,k)
!! GC(:,i,j,k): intercell flux between cell (i,j,k) and cell (i,j+1,k)
!! HC(:,i,j,k): intercell flux between cell (i,j,k) and cell (i,j,k+1)
!! The order of the solver must be specified. First-order employs piece-wise
!! constant primitives in each cell, while second-ordr uses linear
!! reconstruction to interpolate the primitives at cell interfaces.
!> @param locIndx Local index of the block
!> @param order Order of solver (interpolation for primitives)
subroutine HLLfluxes (locIndx, order)

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

            ! X dimension
            pl(:) = PRIM(locIndx,:,i,j,k)
            pr(:) = PRIM(locIndx,:,i+1,j,k)
            call primfhll (pL, pR, ff)
            FC(:,i,j,k) = ff(:)

            ! Y dimension
            pL(:) = PRIM(locIndx,:,i,j,k)
            pR(:) = PRIM(locIndx,:,i,j+1,k)
            call swapxy (pL)
            call swapxy (pR)
            call primfhll (pL, pR, ff)
            call swapxy (ff)
            GC(:,i,j,k) = ff(:)

            ! Z dimension
            pL(:) = PRIM(locIndx,:,i,j,k)
            pR(:) = PRIM(locIndx,:,i,j,k+1)
            call swapxz (pL)
            call swapxz (pR)
            call primfhll (pL, pR, ff)
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

            ! X dimension
            pll(:) = PRIM(locIndx,:,i-1,j,k)
            pl(:)  = PRIM(locIndx,:,i,  j,k)
            pr(:)  = PRIM(locIndx,:,i+1,j,k)
            prr(:) = PRIM(locIndx,:,i+2,j,k)
            call limiter (pll,pl,pr,prr,limiter_type,neqtot)
            call primfhll (pl, pr, ff)
            FC(:,i,j,k) = ff(:)

            ! Y dimension
            pll(:) = PRIM(locIndx,:,i,j-1,k)
            pl(:)  = PRIM(locIndx,:,i,j  ,k)
            pr(:)  = PRIM(locIndx,:,i,j+1,k)
            prr(:) = PRIM(locIndx,:,i,j+2,k)
            call swapxy (pll)
            call swapxy (pl)
            call swapxy (pr)
            call swapxy (prr)
            call limiter (pll,pl,pr,prr,limiter_type,neqtot)
            call primfhll (pl, pr, ff)
            call swapxy (ff)
            GC(:,i,j,k) = ff(:)

            ! Z dimension
            pll(:) = PRIM(locIndx,:,i,j,k-1)
            pl(:)  = PRIM(locIndx,:,i,j,k  )
            pr(:)  = PRIM(locIndx,:,i,j,k+1)
            prr(:) = PRIM(locIndx,:,i,j,k+2)            
            call swapxz (pll)
            call swapxz (pl)
            call swapxz (pr)
            call swapxz (prr)
            call limiter (pll,pl,pr,prr,limiter_type,neqtot)            
            call primfhll (pl, pr, ff)
            call swapxz (ff)
            HC(:,i,j,k) = ff(:)
         
          end if
          
        end do
      end do
    end do
    
  end select
  
end subroutine HLLfluxes

!===============================================================================

!> @brief HLL intercell fluxes along X using primitives at an interface
!> @param pL Vector of primitives at LEFT of interface
!> @param pR Vector of primitives at RIGHT of interface
!> @param ff Returned vector of HLL intercell fluxes
subroutine primfhll (pL, pR, ff)

  use parameters
  implicit none
  
  real, intent(in) :: pL(neqtot)
  real, intent(in) :: pR(neqtot)
  real, intent(out) :: ff(neqtot)
  
  real :: sl, sr, sst
  real :: uL(neqtot), uR(neqtot)
  real :: fL(neqtot), fR(neqtot)

  ! Calculate wavespeeds
  call wavespeed (pL, pR, sl, sr, sst)

!  write(*,*) sl, sr

  ! Obtain intercell fluxes as given by 10.21 of Toro
  if (sl.ge.0) then
    call prim2fluxes (pL, DIM_X, fL)
    ff(:) = fL(:)
    return
  else if (sr.le.0) then
    call prim2fluxes (pR, DIM_X, fR)  
    ff(:) = fR(:)
    return
  else
    call prim2fluxes (pL, DIM_X, fL)
    call prim2fluxes (pR, DIM_X, fR)
    call prim2flow (pL, uL)
    call prim2flow (pR, uR)
    ff(:) = (sr*fl(:)-sl*fr(:)+sl*sr*(uR(:)-uL(:)))/(sr-sl)
    return
  end if

end subroutine primfhll

!===============================================================================
