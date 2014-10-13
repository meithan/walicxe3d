!===============================================================================
!> @file loadbalance.f90
!> @brief Parallel execution load balancing
!> @author Juan C. Toledo
!> @date 13/Jan/2012

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

!> @brief Performs load balancing among all processes (wrapper)
!> @details Redistributes the blocks among all processes, with the purpose of
!! balancing the load. This is the high-level wrapper subroutine, intended
!! to be called by the main program.
subroutine loadBalance ()

  use parameters
  use globals
  use tictoc
  implicit none

  integer :: mark

  call tic(mark)

  write(logu,*) ""
  write(logu,'(1x,a)') "============================================"
  write(logu,'(1x,a)') " Performing Load Balance ..."
  write(logu,'(1x,a)') "============================================"
  write(logu,*) ""

  ! Main load balacing routine
  call doBalance ()

  ! Done
  write(logu,*) ""
  write(logu,'(1x,a,a)') "> Load balance completed in ", nicetoc(mark)
  write(logu,*) ""

end subroutine loadBalance

!===============================================================================

!> @brief Performs load balancing among all processors
!> @details It first builds a loadOrder of the blocks following a Hilbert walk.
!! Then, it builds a loadScheme, which contains a list of all the active
!! blocks and the current and new owner for each block. Finally, it applies
!! the computed loadScheme.
subroutine doBalance ()

  use parameters
  use globals
  use tictoc
  implicit none

  integer :: loadOrder(nbMaxGlobal)
  integer :: numBlocks(0:nProcs-1)
  integer :: loadScheme(nbActive, 3), tots(nProcs)
  integer :: nb, nbloc, p, c, bID, nData
  integer :: old_owner, new_owner, blocks, badcells
  integer :: mpistatus(MPI_STATUS_SIZE)
  integer :: mark
  integer :: sent, received
  logical :: verbose  

  verbose = .false.

  ! The goal is to build and apply the new load balancing scheme, which is
  ! containted in the loadScheme array. Its structure is:
  !   loadScheme(:,1) = block's bID
  !   loadScheme(:,2) = block's current owner
  !   loadScheme(:,3) = block's new owner
  ! After it's been built, it is used walked one entry at a time to transfer
  ! blocks between processes

  call tic(mark)

  ! First, build an ordered list of blocks according to its position in
  ! in a Hilbert walk at the finest level
  loadOrder(:) = -1
  call HilbertOrder (loadOrder)

  ! Determine how many blocks each process gets. Process 0 is assigned last.
  numBlocks(:) = floor(real(nbActive)/nProcs)
  do nb=1,mod(nbActive,nProcs)
    p = mod(nb,nProcs)
    numBlocks(p) = numBlocks(p) + 1
  end do

  ! Initialize the loadScheme list using loadOrder and the current owner
  c = 1
  do nb=1,nbMaxGlobal
    if (loadOrder(nb).ne.-1) then
      bID = loadOrder(nb)
      call getOwner (bID, old_owner)
      loadScheme(c,1) = bID
      loadScheme(c,2) = old_owner
      loadScheme(c,3) = -1
      c = c + 1
    end if
  end do

  ! Next, build the new loadScheme
  nb = 1
  do p=0,nProcs-1
    blocks = numBlocks(p)
    do c=1,blocks
      loadScheme(nb,3) = p
      nb = nb + 1
    end do
  end do

  ! DEBUG
  if (verbose) then
    ! Report the loadScheme
    write(logu,'(1x,a)') "The new load scheme is:"
    write(logu,'(6x,a)') "bID   old owner  ->   new owner"
    tots(:) = 0
    do nb=1,nbActive
      bID = loadScheme(nb,1)
      old_owner = loadScheme(nb,2)
      new_owner = loadScheme(nb,3)
      write(logu,'(1x,i8,i8,6x,a,i8)') bID, old_owner, "->", new_owner
      tots(new_owner+1) = tots(new_owner+1) + 1
    end do
    write(logu,'(1x,a)') "Summary:"
    do p=0,nProcs-1
      write(logu,'(2x,a,i8,a,i4,a)') "Process ", p, " gets ", tots(p+1), " blocks"
    end do
  end if  

  ! Everybody stop here
  call mpi_barrier (mpi_comm_world, ierr)

  ! DEBUG
  if (verbose) then
    write(logu,*) "There are", nbActive, "blocks globally"
    write(logu,*) "I have", nbLocal, "local blocks"
  end if
  ! DEBUG

  ! Finally, distribute the blocks via MPI messages
  !
  ! For every block, compare current and new owner.
  ! Then, one of four things will happen:
  ! 1) if the old owner is also the new owner, ignore this code
  !   (this block doesn't change process)
  ! 2) if they are different, but neither is my rank, ignore this code
  !   (I'm not involved in this exchange)
  ! 3) if they are different, and I'm the old owner, MPI_SEND the block to
  !    the new owner
  ! 4) if they are different, and I'm the new owner, MPI_RECV the block from
  !    the old owner
  call tic(mark)
  sent = 0
  received = 0
  do nb=1,nbActive

    bID = loadScheme(nb,1)
    old_owner = loadScheme(nb,2)
    new_owner = loadScheme(nb,3)

    if (old_owner.ne.new_owner) then
      if (rank.eq.old_owner) then

        ! Transmit this block (flow variables on physical cells)

        call find (bID, localBlocks, nbMaxProc, nbloc)
        if (nbloc.eq.-1) then
          write(logu,*) "Couldn't find bID", bID, "in local list!!"
          write(logu,*) "***ABORTING!***"
          call clean_abort (ERROR_LOADBAL_NO_LOCAL_BID)
        end if

        if (verbose) then        
          write(logu,'(1x,a,i8,a,i4,a,i4)') "Transmitting block ", bID, &
            " with local index ", nbloc, " to process ", new_owner
        end if
        nData = neqtot*ncells_x*ncells_y*ncells_z
        call MPI_SEND( U(nbloc, 1:neqtot, 1:ncells_x, 1:ncells_y, 1:ncells_z), &
             nData, mpi_real_kind, new_owner, bID, mpi_comm_world, ierr)
        localBlocks(nbloc) = -1

        sent = sent + 1
        
      else if (rank.eq.new_owner) then

        ! Receive this block (flow variables on physical cells)

        call find (-1, localBlocks, nbMaxProc, nbloc)        
        if (nbloc.eq.-1) then
          write(logu,*) "Couldn't find a free local slot!"
          write(logu,*) "***ABORTING!***"
          call clean_abort (ERROR_LOADBAL_INSUFFICIENT_RAM)
        end if

        if (verbose) then        
          write(logu,'(1x,a,i8,a,i4,a,i4)') "Receiving block ", bID, &
            " into free local index ", nbloc, " from process ", old_owner
        end if
        nData = neqtot*ncells_x*ncells_y*ncells_z
        call MPI_RECV( U(nbloc, 1:neqtot, 1:ncells_x, 1:ncells_y, 1:ncells_z), &
             nData, mpi_real_kind, old_owner, bID, mpi_comm_world, mpistatus, ierr)
        localBlocks(nbloc) = bID

        ! Update primitives of received block
        call calcPrimsBlock (U, PRIM, nbloc, CELLS_PHYS, badcells)

        received = received + 1
   
      end if
    end if

  end do
  
  ! Re-calculate number of local blocks
  nbLocal = 0
  do nb=1,nbMaxProc
    if (localBlocks(nb).ne.-1) then
      nbLocal = nbLocal + 1
    end if
  end do

  ! Report load balance results
  write(logu,'(1x,a,i0,a,i0,a)') "Sent ", sent, " blocks; received ", received, " blocks"
  write(logu,'(1x,a,i0,a)') "Now have ", nbLocal, " local blocks"

  ! Synchronize global block registry
  call syncBlockLists ()

  ! DEBUG
  if (verbose) then
    write(logu,*) "Local blocks after syncBlockLists:"
    do nb=1,nbMaxProc
      if (localBlocks(nb).ne.-1) then
        write(logu,*) nb, localBlocks(nb)
      end if
    end do
  end if
  
end subroutine doBalance

!===============================================================================

!> @brief Returns a load-balancing ordered list of blocks
!> @details The list is obtained by calculating the index of each block
!! along a Hilbert curve, and then sorting the block list using the indices
!> @param loadOrder The load-balanced ordered list of bIDs
subroutine HilbertOrder (loadOrder)

  use parameters
  use globals
  implicit none

  integer, intent(out) :: loadOrder(nbMaxGlobal)

  integer :: nb, next, bID, istat
  integer :: hk
  integer, allocatable :: blocks(:)
  integer, allocatable :: keys(:)

  ! Re-calculate nbActive, declare blocks(:) and keys(:), and fill them
  ! with blocks and Hilbert Keys

  nbActive = 0
  do nb=1,nbMaxGlobal
    if (globalBlocks(nb).ne.-1) then
      nbActive = nbActive + 1
    end if
  end do

  allocate( blocks(nbActive), stat=istat)
  allocate( keys(nbActive), stat=istat)
  blocks(:) = -1
  keys(:) = -1

  ! Build lists of blocks with corresponding Hilbert indices as keys
  next = 1
  do nb=1,nbMaxGlobal
    if (globalBlocks(nb).ne.-1) then
      ! Obtain this block's Hilbert key at the finest mesh level
      bID = globalBlocks(nb)
      call getHilbertKey(bID, hk)
      blocks(next) = bID
      keys(next) = hk     ! Implicit conversion to integer here
      next = next + 1    
    end if
  end do

  ! Sort 'blocks' list using 'keys' through QuickSort
  call QuickSort (nbActive, blocks, keys, 1, nbActive)

  ! Fill loadOrder list with sorted blocks list
  loadOrder(:) = -1
  do nb=1,nbActive
    loadOrder(nb) = blocks(nb)
  end do

  ! DEBUG
!  if (rank.eq.master) then
!    write(logu,*) "Hilbert indexing:"
!    do nb=1,nbActive
!      bID = loadOrder(nb)
!      call bcoords(bID, x, y, z)
!      write(logu,*) bID, x, y, z, keys(nb)
!    end do
!  end if
  ! DEBUG
  
  ! Don't forget to deallocate these arrays :)
  deallocate (blocks)
  deallocate (keys)

end subroutine HilbertOrder

!===============================================================================

!> @brief Returns a load-balancing ordered list of blocks (test routine)
!> @details The list is obtained by order in the globalBlocks list
!> @param loadOrder The load-balanced ordered list of bIDs
subroutine NaiveOrder (loadOrder)

  use parameters
  use globals
  implicit none

  integer, intent(out) :: loadOrder(nbMaxGlobal)

  integer :: nb, next

  loadOrder(:) = -1
  next = 1
  do nb=1,nbMaxGlobal
    if (globalBlocks(nb).ne.-1) then
      loadOrder(next) = globalBlocks(nb)
      next = next + 1
    end if
  end do

end subroutine NaiveOrder

!===============================================================================

!> @brief Quicksorts a list using an auxiliary list of keys (in-place)
!> @param l The size of the lists
!> @param list The list of items to be sorted
!> @param keys The list of keys used for the sort
!> @param first The index of the first element to be sorted
!> @param last The index of the last element to be sorted
recursive subroutine QuickSort (l, list, keys, first, last)

  implicit none
  integer, intent(in) :: l
  integer, intent(inout) :: list(l)
  integer, intent(inout) :: keys(l)
  integer, intent(in) :: first
  integer, intent(in) :: last

  integer :: pivotIndex, newPivotIndex

  ! Only do something if the lists have at least 2 elements
  if (first.lt.last) then

    ! Choose a pivot - middle item
    pivotIndex = first + (last-first+1)/2  ! Implicit conversion to integer
    
    ! Partition lists in-place
    call Partition (l, list, keys, first, last, pivotIndex, newPivotIndex)

    ! Recursively sort each sublist
    call QuickSort (l, list, keys, first, newPivotIndex - 1)
    call QuickSort (l, list, keys, newPivotIndex + 1, last)  

  end if

end subroutine QuickSort

!===============================================================================

!> @brief In-place auxiliary partition routine for QuickSort
!> @param l The size of the lists
!> @param list The list of items to be sorted
!> @param keys The list of keys used for the sort
subroutine Partition (l, list, keys, first, last, pivotIndex, storeIndex)

  implicit none
  integer, intent(in) :: l
  integer, intent(inout) :: list(l)
  integer, intent(inout) :: keys(l)
  integer, intent(in) :: first
  integer, intent(in) :: last
  integer, intent(in) :: pivotIndex
  integer, intent(out) :: storeIndex

  integer :: pivotKey
  integer :: temp, i

  ! Get the pivot's key
  pivotKey = keys(pivotIndex)

  ! Swap pivot value/key with item in last position
  temp = list(last)
  list(last) = list(pivotIndex)
  list(pivotIndex) = temp
  temp = keys(last)
  keys(last) = keys(pivotIndex)
  keys(pivotIndex) = temp

  ! Swap items to begginning if key smaller than pivot
  storeIndex = first
  do i=first,last-1
    if (keys(i).le.pivotKey) then
      temp = list(storeIndex)
      list(storeIndex) = list(i)
      list(i) = temp
      temp = keys(storeIndex)
      keys(storeIndex) = keys(i)
      keys(i) = temp
      storeIndex = storeIndex + 1
    end if
  end do

  ! Move pivot to its final position
  temp = list(last)
  list(last) = list(storeIndex)
  list(storeIndex) = temp
  temp = keys(last)
  keys(last) = keys(storeIndex)
  keys(storeIndex) = temp

end subroutine Partition

!===============================================================================

!> @brief Calculates the Hilbert Key of a bID at the finest mesh level
!> @details If the block's level is not the finest, 
!> @param bID The global ID of the block
!> @param hk The block's Hilbert key
subroutine getHilbertKey (bID, hk)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: hk

  integer :: maxnb, order, x, y, z, ilev, xp, yp, zp

  ! Determine from the base grid geometry the "bounding cube" for the
  ! Hilbert curve, and hence the order of the curve
  maxnb = max(nbx(maxlev), nby(maxlev), nbz(maxlev))
  order = ceiling(log(1.0*maxnb)/log(2.0))

  ! Obtain the block's local coords (own mesh level) and convert
  ! them to coords at the finest mesh level
  call bcoords(bID, x, y, z)
  call meshlevel(bID, ilev)
  xp = (x-1)*2**(maxlev-ilev) + 1
  yp = (y-1)*2**(maxlev-ilev) + 1
  zp = (z-1)*2**(maxlev-ilev) + 1

  ! Get the Hilbert Key from the subroutine in hilbert.f90
  hk = -1
  call HKey3(xp-1, yp-1, zp-1, order, hk)

end subroutine getHilbertKey
