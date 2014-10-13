!===============================================================================
!> @file boundary.f90
!> @brief Fills ghost cells by calculating or passing block boundaries
!> @author Juan C. Toledo
!> @date 20/Feb/2012

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
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.

!===============================================================================

!> @brief High-level wrapper routine for boundary conditions
!> @details This is a wrapper routine that calls the two boundary-passing
!! routines: normalBoundary, which applies standard boundary conditions,
!! and userBoundary, which applies user-defined custom conditions (this
!! routine is provided by user.f90, a file which must be modified
!! by the user on a per-problem basis).
!> @param depth The number of layers of ghost cells that are to be passed
!> @param uvars Vector of flow variables to be modified
subroutine boundary (depth, uvars)

  use parameters
  use userconds
  implicit none

  integer, intent(in) :: depth
  real, intent(inout) :: uvars(nbMaxProc, neqtot, &
                         nxmin:nxmax, nymin:nymax, nzmin:nzmax)
  
  ! Apply standard boundary conditions on ghost cells
  call normalBoundary (depth, uvars)
  call mpi_barrier (mpi_comm_world, ierr)  ! This barrier might not be needed

  ! Apply user-defined boundary conditions (provided by user.f90)
  call userBoundary (uvars)
  call mpi_barrier (mpi_comm_world, ierr)  ! This barrier might not be needed

end subroutine boundary

!===============================================================================

!> @brief Performs calculation and passing of block boundaries
!> @details This routine exchanges the boundaries between blocks prior to
!! numerical integration. Simulation box boundaries are calculated locally,
!! while boundaries betwen blocks are passed via message.
!> @param depth The number of layers of ghost cells that are to be passed
!> @param uvars Vector of flow variables to be modified
subroutine normalBoundary (depth, uvars)

  use globals
  use parameters
  use tictoc
  implicit none

  integer, intent(in) :: depth
  real, intent(inout) :: uvars(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)

  integer :: a, b, direction, nb, destID, destOwner, destInd, srcID, srcOwner
  integer :: srcInd, nData, src_face, depth1, ieq, bc_type
  integer :: i, j, k, ip, jp, kp, i1, i2, i3, i4, j1, j2, j3, j4, k1, k2, k3, k4
  integer :: sx, sy, sz
  integer :: ntype, neighs(4), mark
  integer :: mpistatus(MPI_STATUS_SIZE)
  real :: x_buf(neqtot, depth, ncells_y/2, ncells_z/2)
  real :: y_buf(neqtot, ncells_x/2, depth, ncells_z/2)
  real :: z_buf(neqtot, ncells_x/2, ncells_y/2, depth)
  logical :: verbose

  ! Debugging flag
  verbose = .false.

  call tic(mark)

  ! For all six directions, where direction indicates which ghost cell layer
  ! of the destination block if being set
  do a=1,6

    select case (a)
    case (1)
      direction = LEFT
      if (verbose) write(logu,'(1x,a)') "Step 1: setting LEFT ghost cells"
    case (2)
      direction = RIGHT
      if (verbose) write(logu,'(1x,a)') "Step 2: setting RIGHT ghost cells"      
    case (3)
      direction = FRONT
      if (verbose) write(logu,'(1x,a)') "Step 3: setting FRONT ghost cells"      
    case (4)
      direction = BACK
      if (verbose) write(logu,'(1x,a)') "Step 4: setting BACK ghost cells"      
    case (5)
      direction = BOTTOM
      if (verbose) write(logu,'(1x,a)') "Step 5: setting BOTTOM ghost cells"
    case (6)
      direction = TOP
      if (verbose) write(logu,'(1x,a)') "Step 6: setting TOP ghost cells"      
    end select

    ! For every active block ...
    do nb=1,nbMaxGlobal
      destID = globalBlocks(nb)
      if (destID.ne.-1) then

        ! The following will set destID's ghost cells on the face along the
        ! direction in turn. The block is called the "destination", while
        ! the neighbor(s) along this direction is(are) called the "source".

        ! Obtain the this block's owner and neighbor(s) along the given direction
        call getOwner (destID, destOwner)
        call neighbors (destID, direction, ntype, neighs)

        if (verbose) then
          write(logu,*) ""
          write(logu,'(a,i8)') "Destination block: ", destID       
          write(logu,'(1x,a,i4)') "Current owner: ", destOwner
          write(logu,'(1x,a,i2)') "Neighbor Type: ", ntype
          write(logu,*) "Neighbor(s): ", neighs
        end if

        ! ====================== !
        ! Sender-side operations !
        !  (MPI operations only) !
        ! ====================== !
        
        if (destOwner.ne.rank) then
          ! Only checked when destination block is *not* local
          !
          ! If it is, this rank will never send boundary cells through MPI
          ! (nonlocal neighbor data will get sent by someone else, and local
          ! neighbor data is copied from local memory)
          
          ! =================================

          ! Boundary to same-level block       
          
          if (ntype.eq.NEIGH_SAME) then
            srcID = neighs(1)
            call getOwner (srcID, srcOwner)
            if (srcOwner.eq.rank) then

              ! Send full boundary through MPI

              call opposite  (direction, src_face)
              call layerLimits (src_face, depth, .false., i1, i2, j1, j2, k1, k2)
              call find (srcID, localBlocks, nbMaxProc, srcInd)

              nData = neqtot*(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
              call MPI_SEND ( uvars(srcInd, 1:neqtot, i1:i2, j1:j2, k1:k2), &
               nData, mpi_real_kind, destOwner, srcID, mpi_comm_world, ierr)

              if (verbose) then
                write(logu,'(1x,i8,a,i3,i3,i3,i3,i3,i3,a,i6,a,i3)') srcID, &
                " is local: MPI sending boundary layer ", i1,i2,j1,j2,k1,k2, &
                " (", nData, " values) to rank ", destOwner
              end if

            else

              if (verbose) then
                write(logu,*) "I don't own neighbor", srcID, "-", srcOwner, " does"
              end if

            end if
          end if
          
          ! =================================

          ! Boundary to finer block (source is coarser than destination)
          
          if (ntype.eq.NEIGH_COARSER) then
            srcID = neighs(1)
            call getOwner (srcID, srcOwner)
            if (srcOwner.eq.rank) then
            
              ! Send one quarter of boundary layer through MPI

              call opposite (direction, src_face)
              depth1 = (depth+1)/2   ! INT division
              call siblingCoords (destID, sx, sy, sz)
              call quadrantLimits (src_face, depth1, .false., sx, sy, sz, i1, i2, j1, j2, k1, k2)
              call find (srcID, localBlocks, nbMaxProc, srcInd)

              nData = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)*neqtot
              call MPI_SEND ( uvars(srcInd, 1:neqtot, i1:i2, j1:j2, k1:k2), &
               nData, mpi_real_kind, destOwner, srcID, mpi_comm_world, ierr)

              if (verbose) then
                write(logu,'(1x,i8,a,i3,i3,i3,i3,i3,i3,a,i6,a,i3)') srcID, &
                " is local: MPI sending boundary quadrant ", i1,i2,j1,j2,k1,k2, &
                " (", nData, " values) to rank ", destOwner
              end if

            else
              if (verbose) then
                write(logu,*) "I don't own neighbor", srcID, "-", srcOwner, " does"
              end if

            end if
          end if

          ! =================================

          ! Boundary to coarser block (source is finer than destination)
          ! Will send up to 4 block boundaries
          
          if (ntype.eq.NEIGH_FINER) then
          
            do b=1,4 
              srcID = neighs(b)
              call getOwner (srcID, srcOwner)
              if (srcOwner.eq.rank) then
              
                ! Average groups of eight cells into buffer before sending

                call opposite (direction, src_face)
                depth1 = depth*2
                call layerLimits (src_face, depth1, .false., i1, i2, j1, j2, k1, k2)
                call find (srcID, localBlocks, nbMaxProc, srcInd)

                select case (direction)
                case (LEFT, RIGHT)
                  x_buf = 0.0
                case (FRONT, BACK)
                  y_buf = 0.0
                case (BOTTOM, TOP)
                  z_buf = 0.0
                end select

                do ieq=1,neqtot
                  do i=i1,i2
                    do j=j1,j2
                      do k=k1,k2
                        ip = (i+1)/2
                        jp = (j+1)/2
                        kp = (k+1)/2
                        select case (direction)
                        case (LEFT, RIGHT)
                          ip = (i-i1+2)/2
                          x_buf(ieq,ip,jp,kp) = x_buf(ieq,ip,jp,kp) + uvars(srcInd,ieq,i,j,k)/8.0
                        case (FRONT, BACK)
                          jp = (j-j1+2)/2
                          y_buf(ieq,ip,jp,kp) = y_buf(ieq,ip,jp,kp) + uvars(srcInd,ieq,i,j,k)/8.0
                        case (BOTTOM, TOP)
                          kp = (k-k1+2)/2
                          z_buf(ieq,ip,jp,kp) = z_buf(ieq,ip,jp,kp) + uvars(srcInd,ieq,i,j,k)/8.0
                        end select
                      end do
                    end do
                  end do
                end do

              ! Now send averaged cells in buffer through MPI
              
              nData = ((i2-i1+1)/2)*((j2-j1+1)/2)*((k2-k1+1)/2)*neqtot
              select case (direction)
              case (LEFT, RIGHT)
                call MPI_SEND ( x_buf, nData, mpi_real_kind, destOwner, &
                  srcID, mpi_comm_world, ierr)
                if (verbose) write(logu,*) "Averaged values range: ", minval(x_buf), maxval(x_buf)
              case (FRONT, BACK)
                call MPI_SEND ( y_buf, nData, mpi_real_kind, destOwner, &
                  srcID, mpi_comm_world, ierr)
                if (verbose) write(logu,*) "Averaged values range: ", minval(x_buf), maxval(x_buf)
              case (BOTTOM, TOP)
                call MPI_SEND ( z_buf, nData, mpi_real_kind, destOwner, &
                  srcID, mpi_comm_world, ierr)
                if (verbose) write(logu,*) "Averaged values range: ", minval(x_buf), maxval(x_buf)
              end select

              if (verbose) then
                write(logu,'(1x,i8,a,i3,i3,i3,i3,i3,i3,a,i6,a,i3)') srcID, &
                " is local: MPI sending averaged boundary layer ", i1,i2,j1,j2,k1,k2, &
                " (", nData, " values) to rank ", destOwner
              end if

            else
              if (verbose) then
                write(logu,*) "I don't own neighbor", srcID, "-", srcOwner, " does"
              end if

              end if
            end do
            
          end if

          ! No SENDER-side operations involved in ntype=NEIGH_BOUNDARY

        end if

        ! ======================== !
        ! Receiver-side operations !
        ! ======================== !

        ! Performed when destination block is LOCAL
        if (destOwner.eq.rank) then

          if (verbose) write(logu,'(1x,a)') "I own the destination block"

          ! =================================

          ! Boundary from same-level block
          
          if (ntype.eq.NEIGH_SAME) then

            srcID = neighs(1)
            call getOwner (srcID, srcOwner)
            call layerLimits (direction, depth, .true., i1, i2, j1, j2, k1, k2)
            call find (destID, localBlocks, nbMaxProc, destInd)

            if (srcOwner.eq.rank) then

              ! LOCAL boundary; ghost cells are copied one-to-one

              call opposite (direction, src_face)
              call layerLimits (src_face, depth, .false., i3, i4, j3, j4, k3, k4)
              call find (srcID, localBlocks, nbMaxProc, srcInd)
              uvars(destInd, 1:neqtot, i1:i2, j1:j2, k1:k2) = uvars(srcInd, 1:neqtot, i3:i4, j3:j4, k3:k4)

              if (verbose) then
                write(logu,'(1x,i8,a,i3,i3,i3,i3,i3,i3,a,i3,i3,i3,i3,i3,i3)') srcID, &
                " is local too: copying SRC boundary layer ", i3,i4,j3,j4,k3,k4, &
                " into DEST ghost layer ", i1,i2,j1,j2,k1,k2
              end if

            else
            
              ! NONLOCAL boundary; ghost cells are received one-to-one through MPI
              
              nData = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)*neqtot
              call MPI_RECV ( uvars(destInd, 1:neqtot, i1:i2, j1:j2, k1:k2), nData, &
                mpi_real_kind, srcOwner, srcID, mpi_comm_world, mpistatus, ierr)

              if (verbose) then
                write(logu,'(1x,i8,a,i8,a,i6,a,i3,i3,i3,i3,i3,i3)') srcID, &
                " is owned by rank ", srcOwner, " : receiving MPI data (", nData, &
                " values) into boundary layer ", i1,i2,j1,j2,k1,k2
              end if

            end if
            
          end if
          
          ! =================================
          
          ! Boundary from coarser block (source is coarser than destination)
          
          if (ntype.eq.NEIGH_COARSER) then

            srcID = neighs(1)
            call getOwner (srcID, srcOwner)
            call layerLimits (direction, depth, .true., i1, i2, j1, j2, k1, k2)
            call find (destID, localBlocks, nbMaxProc, destInd)

            if (srcOwner.eq.rank) then

              ! LOCAL boundary; ghost cells are calculated with duplicated data
              
              call siblingCoords (destID, sx, sy, sz)
              call find (srcID, localBlocks, nbMaxProc, srcInd)
              do ieq=1,neqtot
                do i=i1,i2
                  do j=j1,j2
                    do k=k1,k2
                      ip = (i+1)/2 + sx*ncells_x/2
                      jp = (j+1)/2 + sy*ncells_y/2
                      kp = (k+1)/2 + sz*ncells_z/2
                      select case (direction)
                      case (LEFT)
                        ip = ncells_x + (i+1)/2
                      case (RIGHT)
                        ip = (i-ncells_x+1)/2
                      case (FRONT)
                        jp = ncells_y + (j+1)/2
                      case (BACK)
                        jp = (j-ncells_y+1)/2
                      case (BOTTOM)
                        kp = ncells_z + (k+1)/2
                      case (TOP)
                        kp = (k-ncells_z+1)/2
                      end select
                      uvars(destInd, ieq, i, j, k) = uvars(srcInd, ieq, ip, jp, kp)
                    end do
                  end do
                end do
              end do

              if (verbose) then
                write(logu,'(1x,i8,a,i3,i3,i3,i3,i3,i3)') srcID, &
                " is local too: multiplying SRC boundary layer into DEST ghost layer ", &
                i1,i2,j1,j2,k1,k2
              end if

            else

              ! NONLOCAL boundary; ghost cells received into buffer and duplicated

              ! Receive a quadrant of boundary cells into appropriate buffer
              depth1 = (depth+1)/2
              select case (direction)
              case (LEFT, RIGHT)
                nData = neqtot*depth1*ncells_y/2*ncells_z/2
                call MPI_RECV ( x_buf(1:neqtot, 1:depth1, 1:ncells_y/2, 1:ncells_z/2), &
                  nData, mpi_real_kind, srcOwner, srcID, mpi_comm_world, mpistatus, ierr)
              case (FRONT, BACK)
                nData = neqtot*ncells_x/2*depth1*ncells_z/2
                call MPI_RECV ( y_buf(1:neqtot, 1:ncells_x/2, 1:depth1, 1:ncells_z/2), &
                  nData, mpi_real_kind, srcOwner, srcID, mpi_comm_world, mpistatus, ierr)
              case (BOTTOM, TOP)
                nData = neqtot*ncells_x/2*ncells_y/2*depth1
                call MPI_RECV ( z_buf(1:neqtot, 1:ncells_x/2, 1:ncells_y/2, 1:depth1), &
                  nData, mpi_real_kind, srcOwner, srcID, mpi_comm_world, mpistatus, ierr)
              end select

              ! Set ghost cells with duplicated data
              call siblingCoords (destID, sx, sy, sz)
              do ieq=1,neqtot
                do i=i1,i2
                  do j=j1,j2
                    do k=k1,k2
                      ip = (i+1)/2
                      jp = (j+1)/2
                      kp = (k+1)/2
                      select case (direction)
                      case (LEFT)
                        ip = (2-i)/2
                      case (RIGHT)
                        ip = (i-ncells_x+1)/2
                      case (FRONT)
                        jp = (2-j)/2
                      case (BACK)
                        jp = (j-ncells_y+1)/2
                      case (BOTTOM)
                        kp = (2-k)/2
                      case (TOP)
                        kp = (k-ncells_z+1)/2
                      end select
                      select case (direction)
                      case (LEFT, RIGHT)
                        uvars(destInd, ieq, i, j, k) = x_buf(ieq, ip, jp, kp)
                      case (FRONT, BACK)
                        uvars(destInd, ieq, i, j, k) = y_buf(ieq, ip, jp, kp)
                      case (BOTTOM, TOP)
                        uvars(destInd, ieq, i, j, k) = z_buf(ieq, ip, jp, kp)
                      end select
                    end do
                  end do
                end do
              end do

              if (verbose) then
                write(logu,'(1x,i8,a,i8,a,i6,a,i3,i3,i3,i3,i3,i3)') srcID, &
                " is owned by rank ", srcOwner, " : receiving MPI data (", nData, &
                " values) into buffer and duplicating into boundary layer ", &
                i1,i2,j1,j2,k1,k2
              end if

            end if

          end if

          ! =================================

          ! Boundary from four finer blocks (source is finer than destination)
          
          if (ntype.eq.NEIGH_FINER) then

            ! Set all ghost cells to zero
            call find (destID, localBlocks, nbMaxProc, destInd)
            call layerLimits (direction, depth, .true., i1, i2, j1, j2, k1, k2)
            uvars(destInd,:,i1:i2,j1:j2,k1:k2) = 0.0

            ! For each of the four neighbor blocks ...
            do b=1,4 

              srcID = neighs(b)
              call getOwner (srcID, srcOwner)
              
              if (srcOwner.eq.rank) then

                ! LOCAL boundary: SRC boundary cells are averaged into DEST ghost cells
                
                call siblingCoords (srcID, sx, sy, sz)
                call opposite (direction, src_face)
                depth1 = depth*2
                call layerLimits (src_face, depth1, .false., i3, i4, j3, j4, k3, k4)
                call find (srcID, localBlocks, nbMaxProc, srcInd)

                do ieq=1,neqtot
                  do i=i3,i4
                    do j=j3,j4
                      do k=k3,k4
                        ip = (i+1)/2 + sx*ncells_x/2
                        jp = (j+1)/2 + sy*ncells_y/2
                        kp = (k+1)/2 + sz*ncells_z/2
                        select case (direction)
                        case (LEFT)
                          ip = (i-ncells_x)/2
                        case (RIGHT)
                          ip = (i+1)/2 + ncells_x
                        case (FRONT)
                          jp = (j-ncells_y)/2
                        case (BACK)
                          jp = (j+1)/2 + ncells_y
                        case (BOTTOM)
                          kp = (k-ncells_z)/2
                        case (TOP)
                          kp = (k+1)/2 + ncells_z
                        end select
                        uvars(destInd,ieq,ip,jp,kp) = uvars(destInd,ieq,ip,jp,kp) + uvars(srcInd,ieq,i,j,k)/8.0
                      end do
                    end do
                  end do
                end do                

                if (verbose) then
                  write(logu,'(1x,i8,a,i3,i3,i3,i3,i3,i3,a)') srcID, &
                  " is local too: averaging SRC boundary layer ", i3,i4,j3,j4,k3,k4, &
                  " into DEST ghost layer"
                end if
                
              else

                ! NONLOCAL boundary; pre-averaged ghost cells received directly into DEST quadrant

                call siblingCoords (srcID, sx, sy, sz)
                call quadrantLimits (direction, depth, .true., sx, sy, sz, i1, i2, j1, j2, k1, k2)
                nData = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)*neqtot
                call MPI_RECV ( uvars(destInd, 1:neqtot, i1:i2, j1:j2, k1:k2), nData, &
                  mpi_real_kind, srcOwner, srcID, mpi_comm_world, mpistatus, ierr)

                if (verbose) then
                  write(logu,'(1x,i8,a,i8,a,i6,a,i3,i3,i3,i3,i3,i3)') srcID, &
                  " is owned by rank ", srcOwner, " : receiving MPI data (", nData, &
                  " values) into boundary layer ", i1,i2,j1,j2,k1,k2
                end if
  
              end if
          
            end do
            
          end if

          ! =================================

          ! Computational domain boundary
          ! Calculate ghost cells locally (except PERIODIC BCs)

          ! Pick the corresponding BC for this direction
          select case (direction)
          case (LEFT)
            bc_type = bc_left
          case (RIGHT)
            bc_type = bc_right
          case (FRONT)
            bc_type = bc_front
          case (BACK)
            bc_type = bc_back
          case (BOTTOM)
            bc_type = bc_bottom
          case (TOP)
            bc_type = bc_top
          case default
            write(logu,*) "INVALID BC TYPE=", bc_type
            exit
          end select
          
          if ((ntype.eq.NEIGH_BOUNDARY).and.(bc_type.ne.BC_PERIODIC)) then

            if (verbose) write(logu,'(1x,a)') "Computional domain boundary: calculating locally"

            call layerLimits (direction, depth, .true., i1, i2, j1, j2, k1, k2)
            call find (destID, localBlocks, nbMaxProc, destInd)

            ! For every ghost cell ...
            do i=i1,i2
              do j=j1,j2
                do k=k1,k2

                  ! Reflective (one component of velocity flips sign)
                  if (bc_type.eq.BC_REFLECTIVE) then
                  
                    ip = i
                    jp = j
                    kp = k
                    select case (direction)
                    case (LEFT)
                      ip = 1 - i
                    case (RIGHT)
                      ip = 2*ncells_x - i + 1
                    case (FRONT)
                      jp = 1 - j
                    case (BACK)
                      jp = 2*ncells_y - j + 1
                    case (BOTTOM)
                      kp = 1 - k
                    case (TOP)
                      kp = 2*ncells_z - k + 1
                    end select

                    do ieq=1,neqtot
                      if (((direction.eq.LEFT).or.(direction.eq.RIGHT)).and.(ieq.eq.2)) then
                        uvars(destInd,2,i,j,k) = -1.0 * uvars(destInd,2,ip,jp,kp)
                      else if (((direction.eq.FRONT).or.(direction.eq.BACK)).and.(ieq.eq.3)) then
                        uvars(destInd,3,i,j,k) = -1.0 * uvars(destInd,3,ip,jp,kp)
                      else if (((direction.eq.BOTTOM).or.(direction.eq.TOP)).and.(ieq.eq.4)) then
                        uvars(destInd,4,i,j,k) = -1.0 * uvars(destInd,4,ip,jp,kp)
                      else
                        uvars(destInd,ieq,i,j,k) = uvars(destInd,ieq,ip,jp,kp)
                      end if
                    end do

                  ! Free flow (zero gradient - duplicate the last cell)
                  else if (bc_type.eq.BC_FREE_FLOW) then

                    ip = i
                    jp = j
                    kp = k
                    select case (direction)
                    case (LEFT)
                      ip = 1
                    case (RIGHT)
                      ip = ncells_x
                    case (FRONT)
                      jp = 1
                    case (BACK)
                      jp = ncells_y
                    case (BOTTOM)
                      kp = 1
                    case (TOP)
                      kp = ncells_z
                    end select

                    uvars(destInd,:,i,j,k) = uvars(destInd,:,ip,jp,kp)

                  ! Periodic BCs
                  ! These are handled by the neighbors function, which
                  ! returns the correct neighbor in this case
                  
                  end if


                end do
              end do
            end do
            
          end if

          ! =================================

        end if

      if (verbose) write(logu,'(a,i8)') "Done with block", destID

      end if     
    end do
    
  end do

end subroutine normalBoundary

!===============================================================================

!> @brief Computes the limits (min/max i-j-k values) of boundary or ghost cell layers
!> @param face The face the layer is on
!> @param depth The number of boundary cell layers
!> @param ghost Are these limits meant for ghost cells?
!> @param imin The minimum i-value
!> @param imax The maximum i-value
!> @param jmin The minimum j-value
!> @param jmax The maximum j-value
!> @param kmin The minimum k-value
!> @param kmax The maximum k-value
subroutine layerLimits (face, depth, ghost, imin, imax, jmin, jmax, kmin, kmax)

  use parameters
  implicit none

  integer, intent(in) :: face
  integer, intent(in) :: depth
  logical, intent(in) :: ghost
  integer, intent(out) :: imin, imax, jmin, jmax, kmin, kmax

  ! Set full range initially ...
  imin = 1
  imax = ncells_x
  jmin = 1
  jmax = ncells_y
  kmin = 1
  kmax = ncells_z

  ! ... then crop along depth direction
  select case (face)
  case (LEFT)
    imax = depth
  case (RIGHT)
    imin = ncells_x - depth + 1
  case (FRONT)
    jmax = depth
  case (BACK)
    jmin = ncells_y - depth + 1
  case (BOTTOM)
    kmax = depth
  case (TOP)
    kmin = ncells_z - depth + 1
  case default
    ! Should abort execution
  end select

  ! If these are ghost cells, shift them accordingly
  if (ghost) then
    select case (face)
    case (LEFT)
      imin = imin - depth
      imax = imax - depth
    case (RIGHT)
      imin = imin + depth
      imax = imax + depth
    case (FRONT)
      jmin = jmin - depth
      jmax = jmax - depth
    case (BACK)
      jmin = jmin + depth
      jmax = jmax + depth      
    case (BOTTOM)
      kmin = kmin - depth
      kmax = kmax - depth
    case (TOP)
      kmin = kmin + depth
      kmax = kmax + depth
    case default
      ! Should abort execution
    end select
  end if
  
  return

end subroutine layerLimits

!====================================================================

!> @brief Computes the limits of a subquadrant of a layer of boundary or ghost cells
!> @param face The face the layer is on
!> @param depth The number of boundary cell layers
!> @param ghost Are these limits meant for ghost cells?
!> @param sx The shift of the quadrant along x (0 or 1)
!> @param sy The shift of the quadrant along y (0 or 1)
!> @param sz The shift of the quadrant along z (0 or 1)
!> @param imin The minimum i-value
!> @param imax The maximum i-value
!> @param jmin The minimum j-value
!> @param jmax The maximum j-value
!> @param kmin The minimum k-value
!> @param kmax The maximum k-value
subroutine quadrantLimits (face, depth, ghost, sx, sy, sz, imin, imax, jmin, jmax, kmin, kmax)

  use parameters
  implicit none

  integer, intent(in) :: face
  integer, intent(in) :: depth
  logical, intent(in) :: ghost
  integer, intent(in) :: sx, sy, sz
  integer, intent(out) :: imin, imax, jmin, jmax, kmin, kmax

  ! Set to correct full octant ...
  imin = 1 + (ncells_x/2)*sx
  imax = (ncells_x/2)*(sx+1)
  jmin = 1 + (ncells_y/2)*sy
  jmax = (ncells_y/2)*(sy+1)
  kmin = 1 + (ncells_z/2)*sz
  kmax = (ncells_z/2)*(sz+1)

  ! ... then crop in depth direction
  select case (face)
  case (LEFT)
    imin = 1
    imax = depth
  case (RIGHT)
    imin = ncells_x - depth + 1
    imax = ncells_x
  case (FRONT)
    jmin = 1
    jmax = depth
  case (BACK)
    jmin = ncells_y - depth + 1
    jmax = ncells_y
  case (BOTTOM)
    kmin = 1
    kmax = depth
  case (TOP)
    kmin = ncells_z - depth + 1
    kmax = ncells_z
  case default
    ! Should abort execution
  end select
  
  ! If these are ghost cells, shift them accordingly
  if (ghost) then
    select case (face)
    case (LEFT)
      imin = imin - depth
      imax = imax - depth
    case (RIGHT)
      imin = imin + depth
      imax = imax + depth
    case (FRONT)
      jmin = jmin - depth
      jmax = jmax - depth
    case (BACK)
      jmin = jmin + depth
      jmax = jmax + depth      
    case (BOTTOM)
      kmin = kmin - depth
      kmax = kmax - depth
    case (TOP)
      kmin = kmin + depth
      kmax = kmax + depth
    case default
      ! Should abort execution
    end select
  end if

end subroutine quadrantLimits

!===============================================================================

!> @brief Returns the direction opposite of 'direction'
!> @param direction Given direction
!> @param opposite Opposite direction
subroutine opposite (direction, oppDir)

  use parameters
  implicit none

  integer, intent(in) :: direction
  integer, intent(out) :: oppDir

  select case (direction)
  case (LEFT)
    oppDir = RIGHT
  case (RIGHT)
    oppDir = LEFT
  case (FRONT)
    oppDir = BACK
  case (BACK)
    oppDir = FRONT
  case (BOTTOM)
    oppDir = TOP
  case (TOP)
    oppDir = BOTTOM
  case default
    ! SHOULD ABORT EXECUTION
  end select

  return

end subroutine opposite

!===============================================================================
