!===============================================================================
!> @file admesh.f90
!> @brief Adaptive mesh subroutines
!> @author Juan C. Toledo
!> @date 2/Dic/2011

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

!> @brief High-level routine that triggers mesh refinements
!> @details This routine is called in the main program, and is in charge of
!! deciding and triggering mesh refinement and coarsening, as needed based
!! on the physical data contained by the grid.
subroutine admesh ()

  use parameters
  use globals
  use tictoc
  implicit none

  integer :: nb, bID, flag
  integer :: refcount, crscount
  integer :: localFlags (nbMaxProc)
  integer :: mark
  logical :: test(nbMaxGlobal)  ! DEBUG

  logical :: verbose = .false.

  call tic(mark)
  write(logu,*) ""
  write(logu,'(1x,a)') "============================================"
  write(logu,'(1x,a)') " Updating Adaptive Mesh ..."
  write(logu,'(1x,a)') "============================================"
  write(logu,*) ""

  ! Each process inspects its own active blocks and flags those that meet
  ! the physical criteria for refinement or coarsening
  write(logu,'(1x,a)') "> Flagging blocks due to physical gradients ..."
  refcount = 0
  crscount = 0
  localFlags(:) = FLAG_NONE
  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then
      call markByPhysical (nb, bID, flag)
      localFlags(nb) = flag
      if (flag.eq.FLAG_REFINE) then
        refcount = refcount + 1
      else if (flag.eq.FLAG_COARSE) then
        crscount = crscount + 1
      end if
    end if
  end do

  write(logu,'(1x,i8,a)') refcount, " local blocks marked for refinement"
  write(logu,'(1x,i8,a)') crscount, " local blocks marked for coarsening"  

  ! Syncronize refinement flags
  call mpi_allgather (localFlags, nbMaxProc, mpi_integer, &
    refineFlags, nbMaxProc, mpi_integer, mpi_comm_world, ierr)

  if (verbose) then
    test = refineFlags.eq.FLAG_REFINE
    write(logu,'(1x,i8,a)') count(test), " global blocks marked for refinement"
    test = refineFlags.eq.FLAG_COARSE
    write(logu,'(1x,i8,a)') count(test), " global blocks marked for coarsening"
  end if

  ! Sweep all levels except the finest looking to:
  ! a) flag aditional blocks for refinement, and/or
  ! b) inhibit block coarsening
  ! when proximity to a block marked for refinement requires such change
  write(logu,'(1x,a)') "> Checking proximity criterion ..."
  call checkProximity (maxlev-1)

  ! Refine local blocks marked for refinement
  write(logu,'(1x,a)') "> Refining local blocks ..."  
  do nb=nbmin,nbmax
    bID = globalBlocks(nb)
    flag = refineFlags(nb)
    if (flag.eq.FLAG_REFINE) then
      write(logu,'(1x,a,i8)') "Refining block", bID
      call refineBlock(bID)
    end if
  end do

  ! Update global block list
  call syncBlockLists ()

  ! Do block coarsening
  write(logu,'(1x,a)') "> Coarsening block families ..."    
  call coarseBlocks ()

  ! Update global block list
  call syncBlockLists ()

  write(logu,*) ""
  write(logu,'(1x,a,a)') "> Mesh refined/coarsened in ", nicetoc(mark)

end subroutine admesh

!===============================================================================

!> @brief Marks a block for refinement or coarsening if it meets physical criteria
!> @details Note: this requires updated primitives
!> @param locInd Local index of the block to be checked
!> @param bID bID of the block
!> @param flag Refinement flag
subroutine markByPhysical (locIndx, bID, flag)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: locIndx
  integer, intent(in) :: bID
  integer, intent(out) :: flag
  integer :: i, j, k, ilev
  real :: grad, maxgrad, gradx, grady, gradz
  logical :: verbose

  ! DEBUG
  verbose = .false.

  call meshlevel (bID, ilev)
  grad = 0.0
  maxgrad = 0.0

  do i=1,ncells_x
    do j=1,ncells_y
      do k=1,ncells_z

        ! Pressure gradient check

        gradx = abs(PRIM(locIndx,5,i+1,j,k)-PRIM(locIndx,5,i-1,j,k)) / &
                PRIM(locIndx,5,i,j,k) / (2*dx(ilev))
        grady = abs(PRIM(locIndx,5,i,j+1,k)-PRIM(locIndx,5,i,j-1,k)) / &
                PRIM(locIndx,5,i,j,k) / (2*dx(ilev))
        gradz = abs(PRIM(locIndx,5,i,j,k+1)-PRIM(locIndx,5,i,j,k-1)) / &
                PRIM(locIndx,5,i,j,k) / (2*dx(ilev))
        grad = max( gradx, grady, gradz )
        maxgrad = max( maxgrad, grad )
        
        ! If the gradient is larger than the refinement threshold, mark the block
        ! unless it is already at the maximum level; exit subroutine regardless
        if (grad.ge.refineThres) then
          if (ilev.eq.maxlev) then
            if (verbose) then
              write(logu,'(a,i8,a)') "Block ", bID, " can't be refined any further - you might wanna increase the number of levels"
            end if
            flag = FLAG_NONE
            return
          else if (ilev.lt.maxlev) then       
            flag = FLAG_REFINE
            return
          end if
        end if

        ! Density gradient check
        
        gradx = abs(PRIM(locIndx,1,i+1,j,k)-PRIM(locIndx,1,i-1,j,k)) / &
                PRIM(locIndx,1,i,j,k) / (2*dx(ilev))
        grady = abs(PRIM(locIndx,1,i,j+1,k)-PRIM(locIndx,1,i,j-1,k)) / &
                PRIM(locIndx,1,i,j,k) / (2*dx(ilev))
        gradz = abs(PRIM(locIndx,1,i,j,k+1)-PRIM(locIndx,1,i,j,k-1)) / &
                PRIM(locIndx,1,i,j,k) / (2*dx(ilev))
        grad = max( gradx, grady, gradz )
        maxgrad = max( maxgrad, grad )

        ! If the gradient is larger than the refinement threshold, mark the block
        ! unless it is already at the maximum level; exit subroutine regardless
        if (grad.ge.refineThres) then
          if (ilev.eq.maxlev) then
            if (verbose) then
              write(logu,'(a,i8,a)') "Block ", bID, " can't be refined any further - you might wanna increase the number of levels"
            end if
            flag = FLAG_NONE
            return
          else if (ilev.lt.maxlev) then       
            flag = FLAG_REFINE
            return            
          end if
        end if

      end do
    end do
  end do
  
  ! If the maximum gradient is smaller than the coarsening threshold, mark
  ! it for coarsening, unless it is a root block
  if ((maxgrad.le.coarseThres).and.(ilev.gt.1)) then
    flag = FLAG_COARSE
  else
    flag = FLAG_NONE
  end if

  return

end subroutine markByPhysical

!===============================================================================

! NO LONGER NEEDED
!
!> @brief Marks a block's *neighbors* for refinement by proximity
!> @details The neighbors are checked.
!subroutine markByProximity (bID)

!  use parameters
!  use globals
!  implicit none

!  integer, intent(in) :: bID

!  integer :: dir, nID, fID, fowner, nbf
!  integer :: neighType, neighList(4)

!!  write(logu,*) "PROXIMITY CHECK for bID", bID

!  do dir=1,6
!    call neighbors(bID, dir, neighType, neighList)
!    ! DEBUG
!!    select case (dir)
!!    case (1)
!!      write(logu,*) "Direction:   LEFT"
!!    case (2)
!!      write(logu,*) "Direction:   RIGHT"
!!    case (3)
!!      write(logu,*) "Direction:   FRONT"            
!!    case (4)
!!      write(logu,*) "Direction:   BACK"
!!    case (5)
!!      write(logu,*) "Direction:   BOTTOM"
!!    case (6)
!!      write(logu,*) "Direction:   TOP"            
!!    end select
!!    write(logu,*) "Neighbor Type:", neighType
!!    write(logu,*) "Neighbor ID(s):", neighList
!    ! DEBUG
!  end do

!end subroutine markByProximity

!===============================================================================

!> @brief Refines a block into 8 child blocks
!> @details This requires the following steps:
!! 1) Obtaining the bIDs of all the child blocks
!! 2) Registering the children in the local block registry
!! 3) Copying the father's data octants into the children's data spaces
!! 4) Unregistering the father from the local block registry
!> @param bID The father block's (absolute) ID
subroutine refineBlock (fatherID)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: fatherID

  integer :: ilev, nb, fatherIndex, childIndex, childID
  integer :: c, i, j, k, ip, jp, kp, ieq, sx, sy, sz
  integer :: childList(8)

  ! TODO: THIS CHECK MIGHT BE ELIMINATED FOR EFFICIENCY
  call find (fatherID, localBlocks, nbMaxProc, fatherIndex)
  if (fatherIndex.eq.-1) then
    write(logu,*) ""
    write(logu,'(a,i5,a)') "Block ", fatherID, " is not a local block; can't refine it!"
    write(logu,'(1x,a)') "***ABORTING***"
    call clean_abort (ERROR_LOCAL_BID_NOT_FOUND)
  end if

  ! Don't refine past the highest mesh level
  ! TODO: THIS CHECK MIGHT BE REDUNDANT
  call meshlevel(fatherID, ilev)
  if (ilev.eq.maxlev) then
    write(logu,*) ""
    write(logu,'(a,i5,a)') "Trying to refine block ", fatherID, " past max mesh level! Aborting!"
    write(logu,'(a)') "***ABORTING***"
    call clean_abort (ERROR_ALREADY_MAX_LEV)
  end if

  ! Check if there are enough local block slots
  if (nbLocal+8>nbMaxProc) then
    write(logu,*) ""
    write(logu,'(a,i4)') "Not enough memory space for child blocks of bID ", fatherID
    write(logu,'(a)') "You might want to increase nbMaxProc."
    write(logu,'(a)') "***ABORTING***"
    call clean_abort (ERROR_INSUFICIENT_NBMAXPROC)
  end if

  ! If previous checks passed, go ahead

  ! Get list of children bIDs
  call children (fatherID, childList)

  ! Create children: assign children's bIDs to empty slots in localBlocks
  ! and copy corresponding father's data in children's data space
  do c=1,8
  
    childID = childList(c)

    ! Assign child's bID to empty local block slot
    call put (childID, localBlocks, nbMaxProc, childIndex)
    if (childIndex.eq.-1) then
      write(logu,*) ""
      write(logu,'(a)') "Child block couldn't be assigned to an empty slot!"
      write(logu,'(a)') "***ABORTING***"
      call clean_abort (ERROR_REGISTER_CHILD)
!    else
!      write(logu,'(a,i8,a,i5)') " Child block with bID ", childID, &
!                                " assigned to local slot ", childIndex
    end if
    
    ! Determine child's sibling coordinates (xs, ys, zs)
    call siblingCoords (childID, sx, sy, sz)

    ! Copy father's data into child's data space (flow vars and prims)
    ! For each (i,j,k) cell of the child, (ip,jp,kp) is the corresponding
    ! cell on the father's data
    do i=1,ncells_x
      do j=1,ncells_y
        do k=1,ncells_z
          ip = (ncells_x/2)*sx + (i+1)/2   ! INT division
          jp = (ncells_y/2)*sy + (j+1)/2   ! INT division
          kp = (ncells_z/2)*sz + (k+1)/2   ! INT division
          do ieq=1,neqtot
            U(childIndex,ieq,i,j,k) = U(fatherIndex,ieq,ip,jp,kp)
            PRIM(childIndex,ieq,i,j,k) = PRIM(fatherIndex,ieq,ip,jp,kp)
          end do
        end do
      end do
    end do

  end do

  ! Free father's slot
!  do nb=1,nbMaxProc
!    if (localBlocks(nb).eq.fatherID) then
!      localBlocks(nb) = -1
!!      write(logu,'(a,i8)') " Freeing slot of father block, bID ", fatherID
!      exit
!    end if
!  end do
  call pop (fatherID, localBlocks, nbMaxProc, nb)

  ! Update local block count
  nbLocal = nbLocal + 7

end subroutine refineBlock

!===============================================================================

!> @brief Synchronizes the global block registry
!> @details Simply merges the lists of local blocks
subroutine syncBlockLists()

  use parameters
  use globals
  implicit none

  integer :: nb

#ifndef MPIP

  ! TODO

#else

  write(logu,*) ""
  write(logu,'(1x,a)') "Synchronizing block lists ..."

  ! Reset global block registry
  globalBlocks(:) = -1

  ! Synchronize globalBlocks
  call mpi_allgather (localBlocks, nbMaxProc, mpi_integer, &
      globalBlocks, nbMaxProc, mpi_integer, mpi_comm_world, ierr)

  ! Update nbActive count
  nbActive = 0
  do nb=1,nbMaxGlobal
    if (globalBlocks(nb).ne.-1) then
      nbActive = nbActive + 1
    end if
  end do

  write(logu,'(1x,a,i0,a)') "There are ", nbActive, " active blocks globally."
  
  ! Everybody stops here
  call mpi_barrier (mpi_comm_world, ierr)

#endif  

end subroutine syncBlockLists

!===============================================================================

!> @brief Returns the neighbors of a block at the same level of refinement
!> @details Returns the bID, regardless of whether it is active or not
!> @param bID The (absolute) block ID of which one needs a neighbor
!> @param dir The direction of inquiry, one of the named constants
!! LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK
subroutine neighborLevel(bID, dir, neighID)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(in) :: dir
  integer, intent(out) :: neighID

  integer :: offsets
  integer :: x, y, z, ilev, nx, ny, nz
  integer :: xp, yp, zp

  call bcoords(bID, x, y, z)
  call meshlevel(bID, ilev)
  nx = nbx(ilev)
  ny = nby(ilev)
  nz = nbz(ilev)

  neighID = 0
  xp = x
  yp = y
  zp = z

  select case (dir)

  case (LEFT)
    if (x.eq.1) then
      neighID = -1
      return
    else
      xp = x-1
    end if

  case (RIGHT)
    if (x.eq.nx) then
      neighID = -1
      return
    else
      xp = x+1
    end if

  case (FRONT)
    if (y.eq.1) then
      neighID = -1
      return
    else
      yp = y-1
    end if

  case (BACK)
    if (y.eq.ny) then
      neighID = -1
      return
    else
      yp = y+1
    end if

  case (BOTTOM)
    if (z.eq.1) then
      neighID = -1
      return
    else
      zp = z-1
    end if

  case (TOP)
    if (z.eq.nz) then
      neighID = -1
      return
    else
      zp = z+1
    end if

  case DEFAULT
    print*, "INVALID neighbor direction:", dir
    return
    
  end select

  ! Calculate actual neighbor bID
  neighID = 1 + (xp-1) + (yp-1)*nx + (zp-1)*nx*ny + offsets(ilev)

end subroutine neighborLevel

!===============================================================================

!> @brief Returns the bID of actual active neighbors of a given block
!> @details Returns an integer neighType which describes the kind of
!! neighbors:
!! @n NEIGH_SAME = neighbors is as the same level of refinement
!! @n NEIGH_COARSER = neighbors is coarser
!! @n NEIGH_FINER = higher-level (finer) neighbors (four)
!! @n NEIGH_BOUNDARY = simulation box boundary
!> @param bID The (absolute) block ID of which one needs a neighbor
!> @param dir The direction of inquiry, one of the named constants
!! LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK
!> @param neighType The type of neighbor(s)
!! @param neighList List of neighbor(s) bID(s)
subroutine neighbors (bID, dir, neighType, neighList)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(in) :: dir
  integer, intent(out) :: neighType
  integer, intent(out) :: neighList(4)

  integer :: nID, fID, nb, childList(8)

  neighType = -1
  neighList(:) = -1

  ! Invalid bID
  if (bID.eq.-1) then
    neighType = -1
    return
  end if

  ! Get the direct neighbor's bID and its father's bID
  call neighborLevel(bID, dir, nID)
  call father(nID, fID)

  ! Neighbor is grid boundary
  if (nID.eq.-1) then
    neighType = NEIGH_BOUNDARY
    return
  end if

  ! Check which is active, if any
  do nb=1,nbMaxGlobal
    if (globalBlocks(nb).eq.nID) then
      neighType = NEIGH_SAME
      neighList(1) = nID
      return
    else if ((fID.ne.-1).and.(globalBlocks(nb).eq.fID)) then
      neighType = NEIGH_COARSER
      neighList(1) = fID
      return
    end if
  end do

  ! If neither is found, neighbors are four higher-level blocks, the children
  ! of the same-level neighbor
  neighType = NEIGH_FINER
  call children(nID, childList)

  ! Determine which four of the eight children are neighbors
  select case(dir)
 
    case (LEFT)
      neighList(1) = childList(2)
      neighList(2) = childList(4)
      neighList(3) = childList(6)
      neighList(4) = childList(8)

    case (RIGHT)
      neighList(1) = childList(1)
      neighList(2) = childList(3)
      neighList(3) = childList(5)
      neighList(4) = childList(7)
      
    case (FRONT)
      neighList(1) = childList(3)
      neighList(2) = childList(4)
      neighList(3) = childList(7)
      neighList(4) = childList(8)
      
    case (BACK)
      neighList(1) = childList(1)
      neighList(2) = childList(2)
      neighList(3) = childList(5)
      neighList(4) = childList(6)
      
    case (BOTTOM)
      neighList(1) = childList(5)
      neighList(2) = childList(6)
      neighList(3) = childList(7)
      neighList(4) = childList(8)
      
    case (TOP)
      neighList(1) = childList(1)
      neighList(2) = childList(2)
      neighList(3) = childList(3)
      neighList(4) = childList(4)
      
  end select
  return

end subroutine neighbors

!===============================================================================

!> @brief Returns the refinement level of a block
!> @param bID The (absolute) block ID of which one wants the mesh level
subroutine meshlevel(bID, level)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: level

  if (bID.eq.-1) then
    level = -1
    return
  end if

  level = -1
  do level=1,maxlev
    if ((bID.ge.minID(level)).and.(bID.le.maxID(level))) then
      return
    end if
  end do

end subroutine meshlevel

!===============================================================================

!> @brief Returns the number of blocks in all previous levels
!> @param level The mesh level for which one wants the block offset
integer function offsets(level)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: level
  
  integer :: ilev

  offsets = 0
  do ilev=2,level
    offsets = offsets + nbrootx*nbrooty*nbrootz*8**(ilev-2)
  end do

  return
end function offsets

!===============================================================================

!> @brief Returns the bID of a block's firstborn
!> @param bID The (absolute) block ID of the father
!> @param fbID The (absolute) block ID of the block's firstborn child
subroutine firstborn(bID, firstbornID)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: firstbornID

  integer :: offsets
  integer :: x, y, z, ilev, nx, ny, nz
  integer :: xp, yp, zp

  call bcoords(bID, x, y, z)
  call meshlevel(bID, ilev)
! DEBUG
!write(logu,*) "Calling firstborn..."
!write(logu,*) ilev
!if (ilev.gt.maxlev) then
!  write(logu,*) "Block", bID, "has firstborn beyond the max mesh level!!!"
!  write(logu,*) "Y U ASKING FOR IT?!"
!end if
! DEBUG    
  nx = nbx(ilev+1)
  ny = nby(ilev+1)
  nz = nbz(ilev+1)  

  xp = 2*x-1
  yp = 2*y-1
  zp = 2*z-1

  firstbornID = 1 + (xp-1) + (yp-1)*nx + (zp-1)*nx*ny + offsets(ilev+1)

end subroutine firstborn

!===============================================================================

!> @brief Returns a list of bIDs with all 8 siblings
!> @param bID The block ID of one of the siblings
!> @param sibling_list The list of the 8 (absolute) block IDs of the block's
!! siblings, including itself.
!> @details The siblings are returned in "natural order": first x, then y,
!! last z, producing the following sibling ordering:
!!    ---------
!!    | 7 | 8 |   
!!    ---------   z = top
!!    | 5 | 6 |   
!!    ---------
!!    ---------
!!    | 3 | 4 |   
!!    ---------   z = bottom
!!    | 1 | 2 |   
!! ^  ---------
!! |
!! y  x -->
subroutine siblings(bID, sibling_list)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: sibling_list(8)

  integer :: fatherID, firstbornID, ilev, nx, ny, nz, i

  call father(bID, fatherID)
  
  if (fatherID.eq.-1) then
    do i=1,8
      sibling_list(i) = -1
    end do
    return
  end if
  call firstborn(fatherID, firstbornID)
  call meshlevel(bID, ilev)

  nx = nbx(ilev)
  ny = nby(ilev)
  nz = nbz(ilev)

  sibling_list(1) = firstbornID
  sibling_list(2) = firstbornID + 1
  sibling_list(3) = firstbornID + nx
  sibling_list(4) = firstbornID + nx + 1

  sibling_list(5) = firstbornID + nx*ny
  sibling_list(6) = firstbornID + nx*ny + 1
  sibling_list(7) = firstbornID + nx*ny + nx
  sibling_list(8) = firstbornID + nx*ny + nx + 1
  
end subroutine siblings

!===============================================================================

! >> OBSOLETE
!> @brief Calculates the sibling ID of a block
!> @param bID The (absolute) block ID
!> @param sID The sibling ID of the block, following the "natural" order:
!!    ---------
!!    | 7 | 8 |   
!!    ---------   z = top
!!    | 5 | 6 |   
!!    ---------
!!    ---------
!!    | 3 | 4 |   
!!    ---------   z = bottom
!!    | 1 | 2 |   
!! ^  ---------
!! |   
!! y  x -->
subroutine siblingID(bID, sID)

  use parameters
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: sID

  integer :: sibling_list(8), i, fatherID

  call father(bID, fatherID)
!  print*, "Father of", bID, "is", fatherID  ! DEBUG
  if (fatherID.eq.-1) then
    ! return -1 if block is root block (no father or siblings)
    sID = -1
    return
  end if

  call siblings(bID, sibling_list)

  sID = -1
  do i=1,8
    if (bID.eq.sibling_list(i)) then
      sID = i
      exit
    end if
  end do

end subroutine siblingID
!===============================================================================

!> @brief Calculates the sibling position (coords) of a block
!> @param bID The (absolute) block ID
!> @param sx The sibling's x-coordinate (0 or 1)
!> @param sy The sibling's y-coordinate (0 or 1)
!> @param sz The sibling's z-coordinate (0 or 1)
subroutine siblingCoords (bID, sx, sy, sz)

  use parameters
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: sx, sy, sz

  integer :: sibling_list(8), i, fatherID, sID

  call father(bID, fatherID)
  if (fatherID.eq.-1) then
    ! return -1 if block is root block (no father or siblings)
    sx = -1;  sy = -1;  sz = -1;
    return
  end if

  call siblings(bID, sibling_list)

  sID = -1
  do i=1,8
    if (bID.eq.sibling_list(i)) then
      sID = i
      exit
    end if
  end do
  
  sx = 0
  sy = 0
  sz = 0
  if (mod(sID,2).eq.0) then
    sx = 1
  end if
  if ((sID.eq.3).or.(sID.eq.4).or.(sID.eq.7).or.(sID.eq.8)) then
    sy = 1
  end if
  if (sID.gt.4) then
    sz = 1
  end if

  return
  
end subroutine siblingCoords

!===============================================================================

!> @brief Returns a list of bIDs with all 8 children of a block
!> @param bID The (absolute) block ID of one of the siblings
!> @param sibling_list The list of the 8 (absolute) block IDs of the block's
!! siblings, including itself. These are given in "natural order": first x,
!! then y, then z
subroutine children(bID, children_list)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: children_list(8)

  integer :: firstbornID, ilev, nx, ny, nz

  call firstborn(bID, firstbornID)
  call meshlevel(bID, ilev)
  nx = nbx(ilev+1)
  ny = nby(ilev+1)
  nz = nbz(ilev+1)

  children_list(1) = firstbornID
  children_list(2) = firstbornID + 1
  children_list(3) = firstbornID + nx
  children_list(4) = firstbornID + nx + 1

  children_list(5) = firstbornID + nx*ny
  children_list(6) = firstbornID + nx*ny + 1
  children_list(7) = firstbornID + nx*ny + nx
  children_list(8) = firstbornID + nx*ny + nx + 1
    
end subroutine children

!===============================================================================

!> @brief Returns a block's father's bID
!> @param bID The block's (absolute) bID
!> @param fatherbID The block's father's (absolute) bID
subroutine father(bID, fatherID)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: fatherID

  integer :: offsets
  integer :: x, y, z, ilev, nx, ny, nz
  integer :: xp, yp, zp

  if (bID.eq.-1) then
    fatherID = -1
    return
  end if

  call meshlevel(bID, ilev)

  if (ilev.eq.1) then
  
    fatherID = -1
    return
    
  else
  
    call bcoords(bID, x, y, z)

    xp = floor((x+1)/2.0)
    yp = floor((y+1)/2.0)
    zp = floor((z+1)/2.0)

    nx = nbx(ilev-1)
    ny = nby(ilev-1)
    nz = nbz(ilev-1)

    fatherID = 1 + (xp-1) + (yp-1)*nx + (zp-1)*nx*ny + offsets(ilev-1)

!    print*, "Block", bID, ", x=", x, ", y=", y, ", z=", z  ! DEBUG
!    print*, bID, "'s father coords:", xp, yp, zp  ! DEBUG
!    print*, ''  ! DEBUG
    
  end if

end subroutine father


!===============================================================================

!> @brief Returns the (x,y,z) local integer coordinates of a block at
!! the block's mesh level
!> @param bID The (absolute) block ID of the block
!> @param x The level-specific integer x-coordinate
!> @param y The level-specific integer y-coordinate
!> @param z The level-specific integer z-coordinate
subroutine bcoords(bID, x, y, z)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: x, y, z

  integer :: offsets
  integer :: ilev, nx, ny, nz, localID

  call meshlevel(bID, ilev)
  nx = nbx(ilev)
  ny = nby(ilev)
  nz = nbz(ilev)

  localID = bID - offsets(ilev)
  x = mod(localID,nx)
  if (x.eq.0) x=nx
  y = mod(ceiling(localID*1.0/(nx)),ny)
  if (y.eq.0) y=ny
  z = mod(ceiling(localID*1.0/(nx*ny)),nz)
  if (z.eq.0) z=nz

end subroutine bcoords

!===============================================================================

!> @brief Returns the physical coordinates (in code units) of a block's
!! reference corner
!> @param bID The (absolute) block ID of the block
!> @param x The level-specific integer x-coordinate
!> @param y The level-specific integer y-coordinate
!> @param z The level-specific integer z-coordinate
subroutine getRefCorner(bID, xx, yy, zz)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  real, intent(out) :: xx, yy, zz

  integer :: x, y, z, ilev
  
  call meshlevel(bID, ilev)
  call bcoords(bID, x, y, z)

  xx = (x-1)*ncells_x*dx(ilev)
  yy = (y-1)*ncells_y*dy(ilev)
  zz = (z-1)*ncells_z*dz(ilev)

end subroutine getRefCorner

!===============================================================================

!> @brief Returns the physical position (in code units) of a single cell
!> @param bID The block ID of the host block
!> @param i The cell's integer x-position within the block (1 to ncells_x)
!> @param j The cell's integer y-position within the block (1 to ncells_y)
!> @param k The cell's integer z-position within the block (1 to ncells_z)
!> @param xx The x-position of the cell's physical center, in code units
!> @param yy The y-position of the cell's physical center, in code units
!> @param zz The z-position of the cell's physical center, in code units
subroutine cellPos (bID, i, j, k, x, y, z)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(in) :: i, j, k
  real, intent(out) :: x, y, z

  integer :: ilev, xb, yb, zb

  if (bID.eq.-1) then
    x = -1.0
    y = -1.0
    z = -1.0
    return
  end if

  call meshlevel(bID, ilev)
  call bcoords(bID, xb, yb, zb)

  x = ((xb-1)*ncells_x + (i-1) + 0.5) * dx(ilev)
  y = ((yb-1)*ncells_y + (j-1) + 0.5) * dy(ilev)
  z = ((zb-1)*ncells_z + (k-1) + 0.5) * dz(ilev)

  return

end subroutine cellPos

!===============================================================================

!> @brief Returns the rank of the process who owns this bID
!> @details Returns -1 if the bID is not in globalBlocks
!> @param bID The (absolute) block ID of the block
!> @param process The block owner's rank
subroutine getOwner (bID, owner)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: owner

  integer :: gIndx

  if (bID.eq.-1) then
    owner = -1
    return
  end if

  call find (bID, globalBlocks, nbMaxGlobal, gIndx)
  
  if (gIndx.ne.-1) then
    owner = (gIndx-1)/nbMaxProc
  else
    owner = -1
  end if
  return

end subroutine getOwner

!===============================================================================

!> @brief Ensure the proximity criterion is on all levels less than or equal
!! to finestLev
!> @details This routine will flag blocks for refinement or un-flag blocks
!! for coarsening if they would violate the level proximity condition. Only
!! level finestLev and any coarser levels are checked. Assumesthe refinement
!! flags were previously calculated and synchronized.  
!> @param startLev Initial (finer) level to check
subroutine checkProximity (finestLev)

  use globals
  use parameters
  implicit none

  integer, intent(in) :: finestLev

  integer :: level, nb, nb1, nb2 
  integer :: bID, nID, dir, i, bflag, nflag
  integer :: neighType, neighList(4), owner
  integer :: sib_list(8)  ! DEBUG
  logical :: inhibited

  ! DEBUG
!  write(logu,*) ""
!  write(logu,*) "Starting Proximity checks ..."
  ! DEBUG
  
  do level=finestLev,1,-1
    ! For every level from startLev up to root ...

    ! DEBUG
!    write(logu,*) "Checking blocks in level", level
!    write(logu,*) "GLOBAL refinement list:"
!    do nb=1,nbMaxProc
!      if (gToRefine(nb).ne.-1) then
!        write(logu,*) gToRefine(nb)
!      end if
!    end do
    ! DEBUG

    do nb=1,nbMaxGlobal
      bID = globalBlocks(nb)
      if ((bID.ne.-1).and.(bID.ge.minID(level)).and.(bID.le.maxID(level))) then

        ! This block's current refinement flag
        bflag = refineFlags(nb)

        ! Step 1: apply proximity refinement
        !
        ! For every block in this level currently marked for refinement, check if:
        ! a) an unmarked neighbor must be marked for refinement due to proximity
        ! b) a neighbor marked for coarsening must be inhibited and/or marked
        ! for refinement due to proximity
        ! Note: refinement is never inhibited
        !
        if (bflag.eq.FLAG_REFINE) then
          do dir=1,6
            ! Check all neighbors
          
            call neighbors (bID, dir, neighType, neighList)
            call getOwner (bID, owner)

            ! Inhibit coarsening of same-level neighbors
            if (neighType.eq.NEIGH_SAME) then

              nID = neighList(1)
              call find (nID, globalBlocks, nbMaxGlobal, nb1)
              nflag = refineFlags(nb1)
              
              if (nflag.eq.FLAG_COARSE) then
                refineFlags(nb1) = FLAG_NONE
              end if
              
            ! Force refinement of coarser neighbors
            else if (neighType.eq.NEIGH_COARSER) then

              nID = neighList(1)
              call find (nID, globalBlocks, nbMaxGlobal, nb2)
              nflag = refineFlags(nb2)

              if (nflag.ne.FLAG_REFINE) then
                refineFlags(nb2) = FLAG_REFINE
              end if
                
            end if
          
          end do
        end if

        ! Step 2: inhibit coarsening due to proximity
        !
        ! Now check all blocks marked for coarsening. Inhibit coarsening if:
        ! a) any finer neighbor is not marked for coarsening too
        !
        if (bflag.eq.FLAG_COARSE) then
        
          inhibited = .false.
          ! Check all neighbors
          do dir=1,6

!            call neighbors(bID, dir, neighType, neighList)         
!            if (neighType.eq.NEIGH_FINER) then
!              do i=1,4

!                nID = neighList(i)
!                call find (nID, globalBlocks, nbMaxGlobal, nb1)

!                if (nb1.eq.-1) then
!                  refineFlags(nb) = FLAG_NONE
!                  inhibited = .true.
!                  exit
!                end if

!              end do
!            end if

!            ! Skip remaining checks if already inhibited
!            if (inhibited) exit

            call neighbors(bID, dir, neighType, neighList)         
            if (neighType.eq.NEIGH_FINER) then

              call siblings(neighList(1), sib_list)
              do i=1,8

                nID = sib_list(i)
                call find (nID, globalBlocks, nbMaxGlobal, nb1)
                if (nb1.ne.-1) then
                  nflag = refineFlags(nb1)
                  if (nflag.ne.FLAG_COARSE) then
                    refineFlags(nb) = FLAG_NONE
                    inhibited = .true.
                    exit
                  end if
                end if

              end do
            end if

            ! Skip remaining checks if already inhibited
            if (inhibited) exit

          end do
        end if

      end if
      
    end do
  end do

!  write(logu,*) "Proximity checks done for all levels"   ! DEBUG
 

end subroutine checkProximity

!===============================================================================

!> @brief Refines a zone of the grid to the requested level
!> @details The rectangular zone to be refined is defined by the
!! provided zone array, which contains the bounding box of the zone:
!!  zone(1): left edge
!!  zone(2): right edge
!!  zone(3): front edge
!!  zone(4): back edge
!!  zone(5): bottom edge
!!  zone(6): top edge
!! These must be given in cgs units (cm). The second parameter, z_level,
!! specifies the level at which the zone is to be refined.
subroutine refineZone (zone, z_level)

  use parameters
  use globals
  implicit none

  real, intent(in) :: zone(6)  
  integer, intent(in) :: z_level

  integer :: ilev, nb, nb1, nb3, bID
  integer :: flag
  integer :: localFlags(nbMaxProc)
  real :: z_left, z_right, z_front, z_back, z_bottom, z_top
  real :: b_left, b_right, b_back, b_front, b_bottom, b_top

  ! Verbosity flag (for debugging, or if you're feeling lonely)
  logical :: verbose = .false.

  ! Sanity check
  if (z_level.gt.maxlev) then
    write(logu,*) "Impossible to refine IC zone to level", z_level, "!"
    write(logu,*) "The maximum refinement level is", maxlev
    write(logu,*) "***ABORTING***"
    call clean_abort (ERROR_REFINING_PAST_MAX)
  end if

  ! Unpack zone bounding box
  z_left = zone(1)
  z_right = zone(2)
  z_front = zone(3)
  z_back = zone(4)
  z_bottom = zone(5)
  z_top = zone(6)

  write(logu,*) ""
  write(logu,'(1x,a,i0)') "> Refining zone to level ", z_level
  write(logu,'(1x,a,es12.5,1x,es12.5)') "x-range: ", z_left, z_right
  write(logu,'(1x,a,es12.5,1x,es12.5)') "y-range: ", z_front, z_back
  write(logu,'(1x,a,es12.5,1x,es12.5)') "z-range: ", z_bottom, z_top
 
  ! For all refinement levels up to the requestes level
  do ilev=1,z_level-1

    if (verbose) then
      write(logu,*) ""
      write(logu,*) "Checking blocks in level", ilev
      write(logu,*) "Local blocks:"
      do nb=1,nbMaxProc
        if (localBlocks(nb).ne.-1) then
          write(logu,*) localBlocks(nb)
        end if
      end do
    end if

    ! Clear local refinement flags
    localFlags(:) = -1
 
    ! Flag for refinement all *local* blocks that contain the zone
    do nb=1,nbMaxProc
      bID = localBlocks(nb)
      if ((bID.ne.-1).and.(bID.ge.minID(ilev)).and.(bID.le.maxID(ilev))) then

        ! Get reference corner position and convert to cgs
        call getRefCorner(bID, b_left, b_front, b_bottom)
        b_right = (b_left + dx(ilev)*ncells_x) * l_sc
        b_back = (b_front + dy(ilev)*nCells_y) * l_sc
        b_top = (b_bottom + dz(ilev)*nCells_z) * l_sc
        b_left = b_left * l_sc
        b_front = b_front * l_sc
        b_bottom = b_bottom * l_sc

        if (verbose) then
          write(logu,*) ""
          write(logu,*) "Local block", bID, "has bounding volume"
          write(logu,*) "x:", b_left, b_right
          write(logu,*) "y:", b_front, b_back
          write(logu,*) "z:", b_bottom, b_top
        end if

        ! If block within bounding box, flag it for refinement         
        if ((z_left.lt.b_right).and.(z_right.gt.b_left).and.&
          (z_front.lt.b_back).and.(z_back.gt.b_front).and.&
          (z_bottom.lt.b_top).and.(z_top.gt.b_bottom)) then
          localFlags(nb) = FLAG_REFINE
          if (verbose) write(logu,*) "Block is contained - will be refined"
        else
          localFlags(nb) = FLAG_NONE
          if (verbose) write(logu,*) "Block is not contained"
        end if

      end if
    end do

    if (verbose) then
      write(logu,*) "localFlags after marking"
      write(logu,*) "         bID         flag"
      do nb3=1,nbMaxProc
        if (localFlags(nb3).ne.-1) then
          write(logu,*) localBlocks(nb3), localFlags(nb3)
        end if
      end do
    end if

    ! Syncronize refinement flags
    call mpi_allgather (localFlags, nbMaxProc, mpi_integer, &
    refineFlags, nbMaxProc, mpi_integer, mpi_comm_world, ierr)

    if (verbose) then
      write(logu,*) "refineFlags after synchronization"
      write(logu,*) "         bID         flag"
      do nb3=1,nbMaxGlobal
        if (refineFlags(nb3).ne.-1) then
          write(logu,*) globalBlocks(nb3), refineFlags(nb3)
        end if
      end do
    end if

    ! Check proximity criterion for all blocks in this level and coarser levels
    write(logu,'(1x,a)') "Doing proximity checks ..."
    call checkProximity(ilev)

    if (verbose) then
      write(logu,*) "         bID        owner        flag"
      do nb3=1,nbMaxGlobal
        if (refineFlags(nb3).ne.-1) then
          call getOwner(globalBlocks(nb3), nb1)
          write(logu,*) globalBlocks(nb3), nb1, refineFlags(nb3)
        end if
      end do
    end if
 
    ! Refine final list of local blocks
    write(logu,'(1x,a)') "Refining local blocks ..."
    do nb=nbmin,nbmax
      bID = globalBlocks(nb)
      flag = refineFlags(nb)
      if (flag.eq.FLAG_REFINE) then
      !if ((bID.ne.-1).and.(flag.eq.FLAG_REFINE)) then
        call refineBlock(bID)
      end if
    end do

    ! Update block registry
    call syncBlockLists ()

    ! Do load balance
    call doBalance ()
 
    ! Done for this level - continue with next
  end do

end subroutine refineZone

!===============================================================================

!> @brief Coarsens families of blocks into a single father block
!> @details Determines which families (sets of siblings) are all set to be
!! coarsened, then does the coarsening.
subroutine coarseBlocks ()

  use globals
  use parameters
  implicit none

  integer :: nb, nb1, bID, eID, sID, fID, fatherIndex, locIndx
  integer :: i, j, k, ip, jp, kp, i1, i2, j1, j2, k1, k2, sx, sy, sz
  integer :: next, flag, sib, sib_owner, elder_owner, ieq, badcells
  integer :: nData, mpistatus(MPI_STATUS_SIZE)
  integer :: family (nbActive,2), sib_list(8)
  real :: buf(neqtot, ncells_x/2, ncells_y/2, ncells_z/2)
  logical :: counted, verbose

  verbose = .false.

  ! 1. Build the list of families
  
  ! Sweep refineFlags and add every block marked for coarseningto its family.
  ! The family array has the following meaning:
  !  family(:,1) = bID of the elder block
  !  family(:,2) = number of member of this family marked for coarsening
  family(:,1) = -1
  family(:,2) = 0
  next = 1

  do nb=1,nbMaxGlobal
    flag = refineFlags(nb)
    bID = globalBlocks(nb)
    if ((bID.ne.-1).and.(flag.eq.FLAG_COARSE)) then

      if (verbose) write(logu,*) "Block", bID, "marked for coarsening"
      ! Obtain elder sibling ID
      call siblings (bID, sib_list)
      eID = sib_list(1)

      if (verbose) write(logu,*) "Family of", bID, ":", sib_list
      ! Check if this family is already registered. If so, increase its count.
      counted = .false.
      do nb1=1,next-1
        if (family(nb1,1).eq.eID) then
          if (verbose) write(logu,*) bID, "'s family exists; counting block"
          family(nb1,2) = family(nb1,2) + 1
          counted = .true.
          exit
        end if
      end do

      ! If it's not, add it
      if (.not.counted) then
        family(next,1) = eID
        family(next,2) = 1
        next = next + 1
        if (verbose) write(logu,*) bID, "'s family doesnt exist; adding new family"
      end if

    end if
  end do

  if (verbose) write (logu,*) "Done counting block families"

  ! 2. Apply coarsening for any family were all 8 siblings were marked
  
  do nb=1,nbActive
    
    ! Proceed if the family has 8 marked members
    eID = family(nb,1)
    if ((eID.ne.-1).and.(family(nb,2).eq.8)) then

      call getOwner (eID, elder_owner)
      call siblings (eID, sib_list)
      call father (eID, fID)
      write(logu,'(1x,a,i0)') "Coarsening children of block ", fID
      if (verbose) then
        write(logu,*) "Elder owner:", elder_owner
        write(logu,*) "family:", sib_list
        write(logu,*) "father's bID:", fID
      end if

      ! The owner of the elder sibling allocates space for the new father block
      if (rank.eq.elder_owner) then

        if (verbose) write(logu,*) "Registering new father block", fID
        call put (fID, localBlocks, nbMaxProc, fatherIndex)

        if (fatherIndex.eq.-1) then
          write(logu,*) ""
          write(logu,'(a)') "Father block couldn't be assigned an empty slot!"
          write(logu,'(a)') "***ABORTING***"
          call clean_abort (ERROR_REGISTER_FATHER)
        end if

        if ((fatherIndex.lt.1).or.(fatherIndex.gt.nbMaxProc)) then
          write(logu,*) "Invalid fatherIndex returned!"
          write(logu,*) "fatherIndex:", fatherIndex
          write(logu,*) "***ABORTING***"
          call clean_abort (ERROR_INVALID_LOC_INDEX)
        end if

      end if

      ! Now copy the data of each sibling into the father's data space
      ! The siblings might be local or on another process' memory
      do sib=1,8

        sID = sib_list(sib)
        call getOwner (sID, sib_owner)
        if (verbose) write(logu,*) "Sibling", sID, "is owned by", sib_owner

        ! SENDER-side operations

        ! If rank owns this sibling but does not own the elder, average
        ! the sibling's data and send it through MPI to the elder's owner
        if ((rank.ne.elder_owner).and.(rank.eq.sib_owner)) then

          ! Average sibling's data into buffer
          call find (sID, localBlocks, nbMaxProc, locIndx)
          do i=1,ncells_x/2
            do j=1,ncells_y/2
              do k=1,ncells_z/2
                ip = i*2-1
                jp = j*2-1
                kp = k*2-1
                do ieq=1,neqtot
                  buf(ieq,i,j,k) = sum( U(locIndx,ieq,ip:ip+1,jp:jp+1,kp:kp+1) ) / 8.0
                end do
              end do
            end do
          end do

          ! Send average data through MPI
          nData = neqtot*(ncells_x/2)*(ncells_y/2)*(ncells_z/2)
          call MPI_SEND ( buf, nData, mpi_real_kind, elder_owner, sID, &
            mpi_comm_world, ierr)
          if (verbose) write(logu,*) "Sent", nData, "values to process", elder_owner
          ! Free this sibling's local slot
          call pop (sID, localBlocks, nbMaxProc, nb1)

        end if

        ! RECEIVER-side operations

        ! The owner of the elder sibling either receives this sibling's data
        ! through MPI ...
        if ((rank.eq.elder_owner).and.(rank.ne.sib_owner)) then

          ! Received MPI data into appropriate memory space
          call siblingCoords (sID, sx, sy, sz)
          i1 = 1 + ncells_x/2*sx
          i2 = ncells_x/2*(1+sx)
          j1 = 1 + ncells_y/2*sy
          j2 = ncells_y/2*(1+sy)
          k1 = 1 + ncells_z/2*sz
          k2 = ncells_z/2*(1+sz)
          nData = neqtot*(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
          call MPI_RECV ( U(fatherIndex, 1:neqtot, i1:i2, j1:j2, k1:k2), nData, &
            mpi_real_kind, sib_owner, sID, mpi_comm_world, mpistatus, ierr)
          if (verbose) write(logu,*) "Received", nData, "values into local space", i1,i2,j1,j2,k1,k2

        ! ... or averages and copies it directly from local memory
        else if ((rank.eq.elder_owner).and.(rank.eq.sib_owner)) then

          call siblingCoords (sID, sx, sy, sz)
          call find (sID, localBlocks, nbMaxProc, locIndx)

          do i=1,ncells_x/2
            do j=1,ncells_y/2
              do k=1,ncells_z/2
                i1 = i*2-1
                j1 = j*2-1
                k1 = k*2-1
                ip = i + sx*ncells_x/2
                jp = j + sy*ncells_y/2
                kp = k + sz*ncells_z/2

                do ieq=1,neqtot

                if ((fatherIndex.lt.1).or.(fatherIndex.gt.nbMaxProc)) then
                  write(logu,*) "invalid fatherIndex !!"
                  write(logu,*) "fatherIndex=", fatherIndex
                  call clean_abort (ERROR_INVALID_LOC_INDEX)
                end if
                if (locIndx.eq.0) then
                  write(logu,*) "locIndx = 0 !!"
                  call clean_abort (ERROR_INVALID_LOC_INDEX)
                end if

                  U(fatherIndex,ieq,ip,jp,kp) = sum( U(locIndx,ieq,i1:i1+1,j1:j1+1,k1:k1+1) ) / 8.0

                end do

              end do
            end do
          end do

          if (verbose) write(logu,*) "Averaged and copied local data"

          ! Free this sibling's local slot
          call pop (sID, localBlocks, nbMaxProc, nb1)

        end if

        if (verbose) then
          if ((rank.ne.elder_owner).and.(rank.ne.sib_owner)) then
            write(logu,*) "No action required."
          end if
        end if

      end do

      ! Calculate primitives on new father block
      if (rank.eq.elder_owner) then
        call calcPrimsBlock (U, PRIM, fatherIndex, CELLS_PHYS, badcells)
      end if

    end if
  end do

end subroutine coarseBlocks

!===============================================================================

!> @brief Returns whether the cell coordinates are invalid
!> @details A cell is considered invalid if:
!! a) At least one of its indices exceeds the ranges of physical+ghost cells
!! b) At least two of its indices are in ghost cell range
subroutine validCell (i, j, k, isValid)

  use parameters
  implicit none

  integer, intent(in) :: i
  integer, intent(in) :: j
  integer, intent(in) :: k
  logical, intent(out) :: isValid

  integer :: counter
  
  if ((i.lt.nxmin).or.(i.gt.nxmax).or.   &
      (j.lt.nymin).or.(j.gt.nymax).or.   &
      (k.lt.nzmin).or.(k.gt.nzmax)) then
    isValid = .false.
    return
  end if

  counter = 0
  if ((i.lt.1).or.(i.gt.ncells_x)) counter = counter + 1
  if ((j.lt.1).or.(j.gt.ncells_y)) counter = counter + 1
  if ((k.lt.1).or.(k.gt.ncells_z)) counter = counter + 1  
  if (counter.le.1) then
    isValid = .true.
  else
    isValid = .false.
  end if

end subroutine

!===============================================================================

!> @brief Returns the global cell coords in max-resolution grid
!> @details Given the cell position i,j,k of block bID, returns the
!! ip,jp,kp global coordinates of that cell as if the whole simulation
!! box were at the maximum refinement level. Note that this will use
!! the reference corner (as opposed to the cell center) of the cell 
!! when the block is not at the highest refinement level.
!> @param i,j,k The block-level integer coordinates of the cell
!> @param ip,jp,kp The returned absolute integer coordinates of the cell
!! at the maximum refinement level
subroutine absCoords (bID, i, j, k, ip, jp, kp)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(in) :: i, j, k
  integer, intent(out) :: ip, jp, kp

  integer :: ilev, bx, by, bz

  ! Obtain block coords
  call meshlevel (bID, ilev)
  call bcoords(bID, bx, by, bz)

  ip = (bx-1)*ncells_x*2**(maxlev-ilev) + (i-1)*2**(maxlev-ilev) + 1
  jp = (by-1)*ncells_y*2**(maxlev-ilev) + (j-1)*2**(maxlev-ilev) + 1
  kp = (bz-1)*ncells_z*2**(maxlev-ilev) + (k-1)*2**(maxlev-ilev) + 1

end subroutine absCoords

!===============================================================================

