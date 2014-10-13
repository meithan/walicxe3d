!===============================================================================
!> @file basegrid.f90
!> @brief Constructs the adaptive mesh
!> @author Juan C. Toledo
!> @date 6/Jun/2011

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

!> @brief Initializes root blocks and refines the mesh to prepare it for ICs
!> @details This is a main-level routine. Calculates the mesh geometry
!! (how many root blocks) and number of levels to meet the required resolution,
!! and refines the IC zone.
subroutine basegrid ()

  use parameters
  use globals
  use tictoc
  implicit none

  ! Local variables
  integer :: maxlevx, maxlevy, maxlevz, smalldim
  integer :: nb, bID, ilev, x, y, z, mark
  real :: smallsize

  call tic(mark)

  write(logu,*) ""
  write(logu,*) "============================================"
  write(logu,'(1x,a)') " Initializing the base grid ..."
  write(logu,*) "============================================"
  write(logu,*) ""

  ! Every rank calculates the base grid geometry, based on the selected method
  if (mesh_method.eq.MESH_AUTO) then

    maxcells_x = p_maxcells_x
    maxcells_y = p_maxcells_y
    maxcells_z = p_maxcells_z

    ! Number of root blocks
    smalldim = min(maxcells_x,maxcells_y,maxcells_z)
    nbrootx = maxcells_x / smalldim
    nbrooty = maxcells_y / smalldim
    nbrootz = maxcells_z / smalldim

    ! Calculate number of levels
    maxlevx = int(log(float(maxcells_x)/nbrootx/ncells_x)/log(2.)+1)
    maxlevy = int(log(float(maxcells_y)/nbrooty/ncells_y)/log(2.)+1)
    maxlevz = int(log(float(maxcells_z)/nbrootz/ncells_z)/log(2.)+1)
    maxlev = max(maxlevx,maxlevy,maxlevz)


  else if (mesh_method.eq.MESH_MANUAL) then

    maxlev = p_maxlev
    nbrootx = p_nbrootx
    nbrooty = p_nbrooty
    nbrootz = p_nbrootz

    ! Number of max-resolution cells 
    maxcells_x = ncells_x * 2**(maxlev-1)
    maxcells_y = ncells_y * 2**(maxlev-1)
    maxcells_z = ncells_z * 2**(maxlev-1)

  end if

  ! The maximum allowable number of levels is 10 at the moment, due to the
  ! fact that bIDs are handled as 4-byte (signed) integers. With 10 levels,
  ! one can have up to 14 root blocks.
  if (maxlev.gt.10) then
    write(logu,*) ""
    write(logu,'(1x,a)') "Desired finest-level resolution would require more than 10 levels!"
    write(logu,'(1x,a)') "***ABORTING***"
    call clean_abort (ERROR_TOO_MANY_LEVS)
  end if

  ! Check if physical sizes conform to desired max-resolution cell numbers
  smallsize = min(xphystot, yphystot, zphystot)
  if ((maxcells_x/smalldim.ne.xphystot/smallsize).or.&
      (maxcells_y/smalldim.ne.yphystot/smallsize).or.&
      (maxcells_z/smalldim.ne.zphystot/smallsize)) then
    write(logu,*)
    write(logu,'(a)') "ERROR: The physical size of the simulation box has"//&
      " a different aspect ratio than the requested number of cells!"
    write(logu,'(a,f4.1,a,f4.1,a,f4.1)') "Physical size aspect ratio: ",&
      xphystot/smallsize, " : ", yphystot/smallsize, " : ", zphystot/smallsize
    write(logu,'(a,i3,a,i3,a,i3)') "Cell number aspect ratio:   ",&
      maxcells_x/smalldim, " : ", maxcells_y/smalldim, " : ", maxcells_z/smalldim
    write(logu,'(a)') "> Modify parameters.f90"
    call clean_abort (ERROR_BASEGRID_BAD_ASPECT)
  else
    write(logu,'(1x,a,i2,a,i2,a,i2)') "Grid geometry (root blocks): ", &
      nbrootx, " x ", nbrooty, " x ", nbrootz
    write(logu,'(1x,a,i3)') "Number of root blocks: ", nbrootx*nbrooty*nbrootz
    write(logu,'(1x,a,i2)') "Number of refinement levels: ", maxlev
  end if

  ! Allocate and initialize grid spacing at each level (code units)
  allocate( dx(maxlev) )
  allocate( dy(maxlev) )
  allocate( dz(maxlev) )
  do ilev=1,maxlev
     dx(ilev) = xphystot / (ncells_x * nbrootx * 2.0**(ilev-1)) / l_sc 
     dy(ilev) = yphystot / (ncells_y * nbrooty * 2.0**(ilev-1)) / l_sc
     dz(ilev) = zphystot / (ncells_z * nbrootz * 2.0**(ilev-1)) / l_sc
  end do

  ! Calculate total number of blocks per level
  allocate( nblockslev(maxlev) )
  do ilev=1,maxlev
    nblockslev(ilev) = nbrootx*nbrooty*nbrootz*8**(ilev-1)
  end do

  ! Calculate number of blocks along each dimension at each level
  allocate( nbx(maxlev) )
  allocate( nby(maxlev) )
  allocate( nbz(maxlev) )
  do ilev=1,maxlev
    nbx(ilev) = nbrootx*2**(ilev-1)
    nby(ilev) = nbrooty*2**(ilev-1)
    nbz(ilev) = nbrootz*2**(ilev-1)
  end do

  ! Calculate the bIDs of the first and last blocks at each level
  allocate( minID(maxlev) )
  allocate( maxID(maxlev) )
  minID(1) = 1
  maxID(1) = nblockslev(1)
  do ilev=2,maxlev
    minID(ilev) = maxID(ilev-1) + 1
    maxID(ilev) = maxID(ilev-1) + nblockslev(ilev)
  end do

  nbLocal = 0
  nbActive = nbrootx*nbrooty*nbrootz

  ! ==========================
  ! The following is only performed in cold starts

  if (.not.dowarm) then
  
    ! Activate root blocks and initialize block registry
    write(logu,'(1x,a)') "> Creating root blocks ..."

    ! Calculate bID of root blocks and register them in master's block list
    if (rank.eq.master) then
      nb = 1
      do z = 1,nbrootz
        do y = 1,nbrooty
          do x = 1,nbrootx

            ! bID is sensitive to numbering order
            bID = 1+(x-1)+(y-1)*nbrootx+(z-1)*nbrootx*nbrooty

            ! Register block
            localBlocks(nb) = bID

            nb = nb + 1
            nbLocal = nbLocal + 1

          end do
        end do
      end do
    end if

    ! Synchronize block lists
    call syncBlockLists ()

  else

    write(logu,*) ""
    write(logu,'(1x,a)') "> Skipping root block creation (warm start)"

  end if

  ! ==========================

  write(logu,*) ""
  write(logu,'(1x,a,a)') "> Created base grid in ", nicetoc(mark)
  write(logu,*) ""

end subroutine basegrid
