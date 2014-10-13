!===============================================================================
!> @file cut.f90
!> @brief Real-time data cut extractor
!> @author Juan C. Toledo
!> @date 8/May/2013

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

!> @brief Extracts 2D cuts (in VTK format) from 3D data files
subroutine extractCut (cut_axis, cut_location, nout)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: cut_axis
  real,    intent(in) :: cut_location
  integer, intent(in) :: nout

  character(*), parameter :: outputtpl = "CutD.YYYY"

  ! ==========================

  integer :: ilev, bID, blocksused, nb
  integer :: i, j, k, i1, j1, k1, ip, jp, i_off, j_off, i2, j2
  integer :: plane, nxmap, nymap, nx, ny, cell_count
  character(256) :: filename
  character(1) :: slash
  character(3) :: rankstr
  character(4) :: noutstr
  integer :: l

  real, allocatable :: outmap(:,:,:)
  real, allocatable :: outmap2(:,:,:)

  ! ==========================

  write(logu,'(a)') "> Writing 2D CUT to disk ..."

  ! Allocate output map
  if (cut_axis.eq.AXIS_X) then
    nxmap = nbrooty*2**(maxlev-1)*ncells_y
    nymap = nbrootz*2**(maxlev-1)*ncells_z
  else if (cut_axis.eq.AXIS_Y) then
    nxmap = nbrootx*2**(maxlev-1)*ncells_x
    nymap = nbrootz*2**(maxlev-1)*ncells_z
  else if (cut_axis.eq.AXIS_Z) then
    nxmap = nbrootx*2**(maxlev-1)*ncells_x
    nymap = nbrooty*2**(maxlev-1)*ncells_y
  end if
  allocate( outmap(neqtot,nxmap,nymap) )
  outmap(:,:,:) = 0.0

  if (rank.eq.master) then
    allocate( outmap2(neqtot,nxmap,nymap) )
    outmap2(:,:,:) = 0.0
  end if

  write(logu,'(a)') "Allocated map" ! DEBUG

  ! Loop over all local blocks. Determine if block intersects cut plane
  blocksused = 0
  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then

!      write(logu,'(a,i0)') "Checking local block ", bID  ! DEBUG

      call meshlevel (bID, ilev)
      call getCellPlane (bID, cut_axis, cut_location, plane)

      ! If block intersects, extract cell plane
      if (plane.ne.-1) then

!        write(logu,'(a)') "Block intersects!"  ! DEBUG
!        write(*,*) "Cutplane=", plane    

        blocksused = blocksused + 1
        if (cut_axis.eq.AXIS_X) then
          nx = ncells_y
          ny = ncells_z
        else if (cut_axis.eq.AXIS_Y) then
          nx = ncells_x
          ny = ncells_z
        else if (cut_axis.eq.AXIS_Z) then
          nx = ncells_x
          ny = ncells_y
        end if

        ! Go over cells that intersect cut plane
        do ip=1,nx
          do jp=1,ny

            ! Obtain block-relative coordinates
            if (cut_axis.eq.AXIS_X) then
              i = plane
              j = ip
              k = jp
            else if (cut_axis.eq.AXIS_Y) then
              i = ip
              j = plane
              k = jp
            else if (cut_axis.eq.AXIS_Z) then
              i = ip
              j = jp
              k = plane
            end if

            ! Calculate finest-mesh absolute coords
            call absCoords (bID,i,j,k,i1,j1,k1)

!            write(*,'(a,1x,i0,1x,i0,1x,i0,1x,a)') &
!              "Cell", i, j, k, "absolute coords:"
!            write(*,'(i0,1x,i0,1x,i0)') i1,j1,k1

            ! Reduce absolute coords to 2D
            if (cut_axis.eq.AXIS_X) then
              i1 = j1
              j1 = k1
            else if (cut_axis.eq.AXIS_Y) then
              i1 = i1
              j1 = k1
            else if (cut_axis.eq.AXIS_Z) then
              i1 = i1
              j1 = j1
            end if

            ! Copy data value into output map. Duplicate value 
            ! into multiple cells if block not at highest resolution
            do i_off=0,2**(maxlev-ilev)-1
              do j_off=0,2**(maxlev-ilev)-1

                ! Calculate outmap map coords
                i2 = i1 + i_off
                j2 = j1 + j_off

!                write(logu,*) "Setting cell", i2, j2
                outmap(:,i2,j2) = U(nb,:,i,j,k)
                cell_count = cell_count + 1

              end do
            end do

          end do
        end do

      end if
    end if

  end do

!  write(logu,*) "Density (before reduction):", minval(outmap(1,:,:)), maxval(outmap(1,:,:))
  write(logu,'(a)') "Done extracting data. Reducing across processes ..."

  ! Reduce all extracted data from all processors
  call mpi_reduce (outmap, outmap2, neqtot*nxmap*nymap, &
                   mpi_real_kind, MPI_SUM, master, mpi_comm_world, ierr)

!  if (rank.eq.master) write(logu,*) "Density (after reduction):", minval(outmap2(1,:,:)), maxval(outmap(1,:,:))

  write(logu,'(a)') "Output map successfully reduced."

  call mpi_barrier (mpi_comm_world, ierr)

  ! Master writes out resulting 2D cut
  if (rank.eq.master) then

    ! Generate filename
    l = len_trim(datadir)
    if (datadir(l:l).ne.'/') then
      slash = '/'
    else
      slash = ''
    end if
    write(rankstr,'(I3.3)') rank
    write(noutstr,'(I4.4)') nout
    filename = outputtpl
    call replace (filename, 'XXX', rankstr)
    call replace (filename, 'YYYY', noutstr)

    if (cut_axis.eq.AXIS_X) call replace(filename, 'D', 'X')
    if (cut_axis.eq.AXIS_Y) call replace(filename, 'D', 'Y')
    if (cut_axis.eq.AXIS_Z) call replace(filename, 'D', 'Z')

    write(filename,'(a)') trim(datadir) // trim(slash) // &
                          trim(filename) // ".vtk"

    ! Write output map to file
    write(*,*) "Writing cut map to file ", trim(filename)
    call write2DMapVTK (outmap2, nxmap, nymap, filename)

  end if

end subroutine extractCut

! ======================================================

!> @brief Returns the local cell plane of this block that intersects
!! the given plane, or -1 if it doesn't cut it
subroutine getCellPlane (bID, cut_axis, cut_loc, plane)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(in) :: cut_axis
  real,    intent(in) :: cut_loc
  integer, intent(out) :: plane

  real :: xl, xh, yl, yh, zl, zh, bl, bh
  integer :: ilev

  call bounds (bID, xl, xh, yl, yh, zl, zh)
  call meshlevel (bID, ilev)

  if (cut_axis.eq.AXIS_X) then
    bl = xl
    bh = xh
  else if (cut_axis.eq.AXIS_Y) then
    bl = yl
    bh = yh
  else if (cut_axis.eq.AXIS_Z) then
    bl = zl
    bh = zh
  else
    write(*,'(1x,a)') "Invalid cut axis!!"
  end if

  ! If block intersects cut plane, determine intersection cell plane.
  ! Otherwise, return -1.
  if ((cut_loc.ge.bl).and.(cut_loc.lt.bh)) then
    plane = int((cut_loc-bl)/(dx(ilev))) + 1
  else
    plane = -1
  end if

  return

end subroutine getCellPlane

! ======================================================

!> @brief Returns the bounding box of a block in physical coordinates
subroutine bounds(bID, xl, xh, yl, yh, zl, zh)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  real, intent(out) :: xl, xh, yl, yh, zl, zh

  integer :: ilev, bx, by, bz

  call meshlevel (bID, ilev)
  call bcoords(bID, bx, by, bz)

  xl = (bx-1)*xphystot/(nbrootx*2**(ilev-1))
  xh = bx*xphystot/(nbrootx*2**(ilev-1))
  yl = (by-1)*yphystot/(nbrooty*2**(ilev-1))
  yh = by*yphystot/(nbrooty*2**(ilev-1))
  zl = (bz-1)*zphystot/(nbrootz*2**(ilev-1))
  zh = bz*zphystot/(nbrootz*2**(ilev-1))

  return

end subroutine bounds

! ======================================================


