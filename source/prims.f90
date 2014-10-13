!===============================================================================
!> @file prims.f90
!> @brief Calculation of primitives and conversion from/to flow variables
!> @author Juan C. Toledo
!> @date 30/Nov/2011

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

!> @brief Update primitives in *all* local blocks
!> @details This high-level wrapper routine updates primitives on all
!! cells of all local blocks, using the U flow vars array as source.
subroutine updatePrims ()

  use parameters
  use globals
  use tictoc
  implicit none
 
  integer :: mark
  
  write(logu,*) ""
  write(logu,'(1x,a)') "> Updating primitive variables ..."
  call tic(mark)
  call calcPrimsAll (U, PRIM, CELLS_ALL)
  write(logu,'(1x,a,a)') "> Updated primitives in", nicetoc(mark)

end subroutine updatePrims

!===============================================================================

!> @brief Wrapper routine to update primitives for all local blocks
!> @details The full vectors of flow variables and primitives must be given.
!! The 'cells' argument indicates the range of cells to update. See the
!! documentation of the calcPrimsBlock routine for further details.
subroutine calcPrimsAll (uvars, pvars, cells)

  use parameters
  use globals
  implicit none

  real, intent(in) :: uvars(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)
  real, intent(out) :: pvars(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)  
  integer, intent(in) :: cells

  integer :: nb, bID, badcells
  logical :: verbose

  verbose = .true.

  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then

      badcells = 0
      call calcPrimsBlock (uvars, pvars, nb, cells, badcells)

      if ((verbose).and.(badcells.ne.0)) then
        write(logu,'(1x,a,i0,a,i0)') "Warning: ", badcells, &
        " pressure corrections in block ", bID
      end if

    end if    
  end do


end subroutine calcPrimsAll

!===============================================================================

!> @brief Calculates and update primitives in one block
!> @details The full vectors of flow variables and primitives must be given.
!! The 'cells' argument indicates the range of cells to update. Currently
!! recognized options are:
!!  CELLS_ALL: update physical and ghost cells
!!  CELLS_PHYS: update only physical cells
!!  CELLS_GHOST: update only ghost cells 
!! The 'badcells' argument returns the number of cells in which 
!! pressure corrections were needed.
subroutine calcPrimsBlock (uvars, pvars, locIdx, cells, badcells)

  use parameters
  use globals
  implicit none

  real, intent(in) :: uvars(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)
  real, intent(out) :: pvars(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)  
  integer, intent(in) :: locIdx
  integer, intent(in) :: cells
  integer, intent(out) :: badcells

  integer :: i, j, k, istat

  badcells = 0


  ! Physical cells
  if ((cells.eq.CELLS_ALL).or.(cells.eq.CELLS_PHYS)) then

    do i=1,ncells_x
      do j=1,ncells_y
        do k=1,ncells_z
          call flow2prim (uvars(locIdx,:,i,j,k), pvars(locIdx,:,i,j,k), istat)
          if (istat.ne.0) badcells = badcells + 1
        end do
      end do
    end do     

  end if

  ! Ghost cells (will skip calculation if density is zero)
  if ((cells.eq.CELLS_ALL).or.(cells.eq.CELLS_GHOST)) then

    ! Left ghost cells
    do i=1-nghost,0
      do j=1,ncells_y
        do k=1,ncells_z
          if (uvars(locIdx,1,i,j,k).ne.0.0) then
            call flow2prim (uvars(locIdx,:,i,j,k), pvars(locIdx,:,i,j,k), istat)
          end if
        end do
      end do
    end do     

    ! Right ghost cells
    do i=ncells_x+1,ncells_x+2
      do j=1,ncells_y
        do k=1,ncells_z
          if (uvars(locIdx,1,i,j,k).ne.0.0) then
            call flow2prim (uvars(locIdx,:,i,j,k), pvars(locIdx,:,i,j,k), istat)
          end if
        end do
      end do
    end do     

    ! Front ghost cells
    do i=1,ncells_x
      do j=1-nghost,0
        do k=1,ncells_z
          if (uvars(locIdx,1,i,j,k).ne.0.0) then
            call flow2prim (uvars(locIdx,:,i,j,k), pvars(locIdx,:,i,j,k), istat)
          end if
        end do
      end do
    end do     

    ! Back ghost cells
    do i=1,ncells_x
      do j=ncells_y+1,ncells_y+2
        do k=1,ncells_z
          if (uvars(locIdx,1,i,j,k).ne.0.0) then
            call flow2prim (uvars(locIdx,:,i,j,k), pvars(locIdx,:,i,j,k), istat)
          end if
        end do
      end do
    end do     

    ! Bottom ghost cells
    do i=1,ncells_x
      do j=1,ncells_y
        do k=1-nghost,0
          if (uvars(locIdx,1,i,j,k).ne.0.0) then
            call flow2prim (uvars(locIdx,:,i,j,k), pvars(locIdx,:,i,j,k), istat)
          end if
        end do
      end do
    end do     

    ! Bottom ghost cells
    do i=1,ncells_x
      do j=1,ncells_y
        do k=ncells_z+1,ncells_z+2
          if (uvars(locIdx,1,i,j,k).ne.0.0) then
            call flow2prim (uvars(locIdx,:,i,j,k), pvars(locIdx,:,i,j,k), istat)
          end if
        end do
      end do
    end do     

  end if

  if ((cells.ne.CELLS_ALL).and.(cells.ne.CELLS_PHYS).and.(cells.ne.CELLS_GHOST)) then
    write(logu,*) "Invalid cell range passed to updatePrimBlock!"
    write(logu,*) "***Aborting!***"
    call clean_abort (ERROR_INVALID_CELL_RANGE)
  end if

end subroutine calcPrimsBlock

!===============================================================================

!> @brief Estimates temperature from primitives, in CGS
!> @details Will use the mean atomic mass for ionized gas if the estimated
!! temperature (using the neutral mean atomic mass) exceeds ion_thres
!> @param pvars Vector of primitives
subroutine calcTemp (pvars, temp)

  use parameters
  use globals
  implicit none

  real, intent(in) :: pvars(neqtot)
  real, intent(out) :: temp

  temp = pvars(5)/pvars(1)*(mu0*AMU*p_sc/d_sc/KB)

  if (temp.gt.ion_thres) then
    temp = pvars(5)/pvars(1)*(mui*AMU*p_sc/d_sc/KB)
  end if

  return

end subroutine calcTemp

!===============================================================================

!> @brief Calculates the vector of conserved (flow) vars given the vector of
!! primitive variables. Passive scalars are simply copied over.
!> @param pvars(neqtot) An input vector containing (rho,u,v,w,P,s1,s2,...)
!> @param uvars(neqtot) An output vector containing (rho,rho*u,rho*v,rho*w,E,s1,s2,...)
subroutine prim2flow (pvars, uvars)

  use parameters
  implicit none

  real, intent(in) :: pvars(neqtot)
  real, intent(out) :: uvars(neqtot)

  real :: v2

  uvars(1) = pvars(1)
  uvars(2) = pvars(1)*pvars(2)
  uvars(3) = pvars(1)*pvars(3)
  uvars(4) = pvars(1)*pvars(4)

  v2 = pvars(2)**2 + pvars(3)**2 + pvars(4)**2
  uvars(5) = 0.5*pvars(1)*v2 + CV*pvars(5)

#ifdef PASBP
  uvars(6) = pvars(6)
  uvars(7) = pvars(7)
  uvars(8) = pvars(8)
#endif

  if (npassive.ge.1) then
    uvars(firstpas:neqtot) = pvars(firstpas:neqtot)
  end if

  return

end subroutine prim2flow


!===============================================================================

!> @brief Calculates the vector of primitive variables given the vector of
!! conserved (flow) variables. Any passive scalars are simply copied over.
!> @param uvars(neqtot) An input vector containing (rho,rho*u,rho*v,rho*w,E,...)
!> @param pvars(neqtot) An output vector containing (rho,u,v,w,P,...)
!> @param istat An output integer indicating the success status of the
!! operation (0=all good, 1=pressure floored, 2=density floored).
subroutine flow2prim (uvars, pvars, istat)

  use parameters
  use globals   ! DEBUG
  implicit none

  real, intent(in) :: uvars(neqtot)
  real, intent(out) :: pvars(neqtot)
  integer, intent(out) :: istat

  real :: rhov2

  istat = 0

  if (uvars(1).eq.0.0) then
    write(logu,*) "Received zero density!!"
    call clean_abort(ERROR_DIVISION_BY_ZERO)
  end if

  pvars(1) = uvars(1)
  pvars(2) = uvars(2)/uvars(1)
  pvars(3) = uvars(3)/uvars(1)
  pvars(4) = uvars(4)/uvars(1)

  rhov2 = (uvars(2)**2 + uvars(3)**2 + uvars(4)**2)/uvars(1)
  pvars(5) = (uvars(5)-0.5*rhov2)/CV

  ! Floor on pressure
  if (pvars(5).lt.1.0e-30) then
    pvars(5) = 1.0e-30
    istat = 1
  end if

  ! Floor on density
  if (pvars(1).lt.1.0e-40) then
    pvars(1) = 1.0e-40
    istat = 2
  end if

#ifdef PASBP
  pvars(6) = uvars(6)
  pvars(7) = uvars(7)
  pvars(8) = uvars(8)
#endif
  
  if (npassive.ge.1) then
    pvars(firstpas:neqtot) = uvars(firstpas:neqtot)
  end if

  return
  
end subroutine flow2prim

!===============================================================================

!> @brief Calculates the vectors of physical fluxes along a specific dimension
!! given the vector of primitive variables. 
!> @param pvars(neqtot) An input vector containing primitives
!> @param dimf Dimension along which fluxes are to be calculated
!> @param F(neqtot) An output vector containing fluxes along dimf
subroutine prim2fluxes (pvars, dimens, flux)

  use parameters
  implicit none

  real, intent(in) :: pvars(neqtot)
  integer, intent(in) :: dimens
  real, intent(out) :: flux(neqtot)

  real :: v2, etot
  real :: pvars1(neqtot)

  ! Make a copy of the primitives so the original ones are not split
  pvars1(:) = pvars(:)

  ! Swap velocity components when dimf=2 (Y) or dimf=3(Z)
  if (dimens.eq.DIM_Y) then
    call swapxy(pvars1)
  else if (dimens.eq.DIM_Z) then
    call swapxz(pvars1)
  end if

  v2 = pvars1(2)**2 + pvars1(3)**2 + pvars1(4)**2
  etot = 0.5*pvars1(1)*v2 + CV*pvars1(5)

  ! Calculate Fluxes (formulae for X)
  flux(1) = pvars1(1)*pvars1(2)
  flux(2) = pvars1(5) + pvars1(1)*(pvars1(2)**2)
  flux(3) = pvars1(1)*pvars1(2)*pvars1(3)
  flux(4) = pvars1(1)*pvars1(2)*pvars1(4)
  flux(5) = pvars1(2)*(etot+pvars1(5))
#ifdef PASBP
  flux(6) = 0
  flux(7) = pvars1(2)*pvars1(7) - pvars1(3)*pvars1(6)
  flux(8) = pvars1(2)*pvars1(8) - pvars1(4)*pvars1(6)
#endif

  if (npassive.ge.1) then
    flux(firstpas:neqtot) = pvars1(2)*pvars1(firstpas:neqtot)
  end if

  ! Swap flux components to get correct fluxes for Y or Z
  if (dimens.eq.DIM_Y) then
    call swapxy(flux)
  else if (dimens.eq.DIM_Z) then
    call swapxz(flux)
  end if

  return

end subroutine prim2fluxes

!===============================================================================

!> @brief Swaps the x and y components of a vector of hydro variables
!> @param vec The vector to be swapped
subroutine swapxy (vec)

  use parameters
  implicit none

  real, intent(inout) :: vec(neqtot)

  real :: temp

  temp = vec(2)
  vec(2) = vec(3)
  vec(3) = temp

#ifdef PASBP
  temp = vec(6)
  vec(6) = vec(7)
  vec(7) = temp
#endif  

  return

end subroutine swapxy

!===============================================================================

!> @brief Swaps the x and z components of a vector of hydro variables
!> @param vec The vector to be swapped
subroutine swapxz (vec)

  use parameters
  implicit none

  real, intent(inout) :: vec(neqtot)

  real :: temp

  temp = vec(2)
  vec(2) = vec(4)
  vec(4) = temp

#ifdef PASBP
  temp = vec(6)
  vec(6) = vec(8)
  vec(8) = temp
#endif  

  return

end subroutine swapxz
