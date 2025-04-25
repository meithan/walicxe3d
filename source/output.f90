!===============================================================================
!> @file output.f90
!> @brief Subroutines to read and write output files
!> @author Juan C. Toledo
!> @date 7/Jun/2011

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

!> @brief Output simulation data in any format required
!> param noutput The output number
subroutine output (noutput)

  use parameters
  use globals
  use tictoc
  implicit none

  integer, intent(in) :: noutput

  integer(8) :: mark, startmark
  integer :: p

  write(logu,*) ""
  write(logu,*) "============================================"
  write(logu,'(1x,a)') " Writing data ouput ..."
  write(logu,*) "============================================"
  write(logu,'(1x,a)') stamp()

  if (output_mode.eq.OUT_SIMULT) then
    write(logu,'(1x,a)') "Output mode: simultaneous"
  else if (output_mode.eq.OUT_TURNS) then
    write(logu,'(1x,a)') "Output mode: turn-based"
  else
    write(logu,'(1x,a,i0)') "Unrecognized data output mode: ", output_mode
    write(logu,'(1x,a)') "***ABORTING***"
    call clean_abort(ERROR_BAD_OUTPUT_MODE)
  end if

  if (.not.(output_bin.or.output_vtk)) then
    write(logu,'(1x,a)') "No output format was selected in parameters.f90!"
    write(logu,'(1x,a)') "Skipping data output."
    write(logu,*) ""
    return
  end if

  call tic(startmark)

  ! Write parameters file on first output
  if (noutput.eq.0) then
    if (rank.eq.0) then
      call writeParams ()
    end if
  end if

  ! Write data and state files
  do p=0,nprocs-1
    if (rank.eq.p) then

      call tic(mark)
      write(logu,'(1x,a,i0,a)') "> Writing output ", noutput, " to disk ..."
      write(logu,*) ""

      ! Binary data format - writes flow vars (Us)
      if (output_bin) then
        call writeBin (noutput)
      end if

      ! VTK format
      if (output_vtk) then
        call updatePrims ()     ! Make sure primitives are updated
        call write3DVTKBlocks (noutput)
        if (rank.eq.master) then
          call write3DVTKGrid (noutput)
        end if
      end if

      ! Write state file - master
      if (rank.eq.master) then
        call writeState (noutput)
      end if

      write(logu,'(1x,a,i0,a,a)') "> Completed output ", noutput, &
      " in ", nicetoc(mark)
      write(logu,'(1x,a)') stamp()
      write(logu,*) ""

    else

      if (output_mode.eq.OUT_TURNS) then
        write(logu,'(1x,a,i0,a)') "Process ", p, " writing output ..."
        write(logu,'(1x,a)') stamp()
        write(logu,*) ""
        call mpi_barrier (mpi_comm_world, ierr)
      end if

    end if
  end do

  call mpi_barrier (mpi_comm_world, ierr)
  write(logu,'(1x,a,i0,a,a)') "> All ranks completed output ", noutput, &
  " in ", nicetoc(startmark)
  write(logu,'(1x,a)') stamp()
  write(logu,*) ""

end subroutine output

!===============================================================================

!> @brief Dump data blocks in native binary format
subroutine writeBin (noutput)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: noutput

  character(100) :: blocksfile
  character(3) :: rankstr
  character(4) :: noutstr
  character(1) :: slash
  integer :: nb, l, unitout, istat, nblocks

  ! Generate filenames based on templates
  l = len_trim(datadir)
  if (datadir(l:l).ne.'/') then
    slash = '/'
  else
    slash = ''
  end if

  write(rankstr,'(I3.3)') rank
  write(noutstr,'(I4.4)') noutput
  blocksfile = blockstpl
  call replace (blocksfile, 'XXX', rankstr)
  call replace (blocksfile, 'YYYY', noutstr)
  write(blocksfile,'(a)') trim(datadir) // trim(slash) // trim(blocksfile) // ".bin"

  ! Open data file
  unitout = 10 + rank
  open (unit=unitout, file=blocksfile, status='replace', access='stream', iostat=istat)
  if (istat.ne.0) then
    write(logu,'(a,a,a)') "Could not open the file '", trim(blocksfile), "' !"
    write(logu,'(a,a,a)') "Does the datadir '", trim(datadir), "' exist?"
    close(unitout)
    call clean_abort (ERROR_OUTPUT_FILE_OPEN)
  end if

  ! Write data
  write(logu,'(1x,a,a,a)') "Writing data blocks to file ", trim(blocksfile), " ..."
  nblocks = 0
  write(unitout) nbLocal
  do nb=1,nbMaxProc
    if (localBlocks(nb).ne.-1) then
      write(unitout) localBlocks(nb)
      write(unitout) U(nb,:,1:ncells_x,1:ncells_y,1:ncells_z)
      nblocks = nblocks + 1
    end if
  end do

  close(unitout)
  write(logu,'(1x,a,i0,a)') "Dumping ", nblocks, " local blocks"
  write(logu,'(1x,a,a)') "Succesfully wrote data file."
  write(logu,*) ""

end subroutine writeBin

!===============================================================================

!> @brief Write state file in binary format (only one per output)
subroutine writeState (noutput)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: noutput

  character(100) :: statefile
  character(4) :: noutstr
  character(1) :: slash
  integer :: l, unitout, istat

  ! Generate filenames based on templates
  l = len_trim(datadir)
  if (datadir(l:l).ne.'/') then
    slash = '/'
  else
    slash = ''
  end if
  write(noutstr,'(I4.4)') noutput
  statefile = statetpl
  call replace (statefile, 'YYYY', noutstr)
  write(statefile,'(a)') trim(datadir) // trim(slash) // trim(statefile) // ".dat"

  write(logu,'(1x,a,a,a)') "Writing state file ", trim(statefile), " ..."

  ! Open data file
  unitout = 10 + rank
  open (unit=unitout, file=statefile, status='replace', iostat=istat)
  if (istat.ne.0) then
    write(logu,'(a,a,a)') "Could not open the file '", trim(statefile), "' !"
    write(logu,'(a,a,a)') "Does the datadir '", trim(datadir), "' exist?"
    close(unitout)
    call clean_abort (ERROR_OUTPUT_FILE_OPEN)
  end if

  ! Write simulation state variables and datadir location
  write(unitout,'(es22.15,i8,i5)') time*t_sc, it, noutput
  write(unitout,'(a)') datadir

  close(unitout)
  write(logu,'(1x,a,i3,a)') "Succesfully wrote state file."
  write(logu,*) ""

end subroutine writeState

!===============================================================================

!> @brief Write the simulations parameters file
subroutine writeParams ()

  use parameters
  use globals
  implicit none

  character(128) :: filepath
  integer :: unitout, istat

  ! Generate file path
  write(filepath,'(a)') trim(datadir) // "/" // "Params.dat"

  write(logu,'(1x,a,a,a)') "Writing state file ", trim(filepath), " ..."

  ! Open data file
  unitout = 100
  open (unit=unitout, file=filepath, status='replace', iostat=istat)
  if (istat.ne.0) then
    write(logu,'(a,a,a)') "Could not open the file '", trim(filepath), "' !"
    write(logu,'(a,a,a)') "Does the datadir '", trim(datadir), "' exist?"
    close(unitout)
    call clean_abort (ERROR_OUTPUT_FILE_OPEN)
  end if

  ! Write simulation parameters
  write(unitout,'(i0,1x,i0,1x,i0,1x,i0)') nprocs, neqtot, npassive, firstpas
  write(unitout,'(es22.15,1x,es22.15,1x,es22.15)') xphystot, yphystot, zphystot
  write(unitout,'(i0,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0)') nbrootx, nbrooty, nbrootz, maxlev, ncells_x, ncells_y, ncells_z
  write(unitout,'(es22.15,1x,es22.15,1x,es22.15,1x,es22.15)') gamma, mu0, mui, ion_thres
  write(unitout,'(es22.15,1x,es22.15,1x,es22.15,1x,i0)') l_sc, d_sc, v_sc, units_type

  close(unitout)
  write(logu,'(1x,a,i3,a)') "Succesfully wrote parameters file."
  write(logu,*) ""
  
end subroutine writeParams

!===============================================================================

!> @brief Visit-compatible binary VTK data files
subroutine write3DVTKBlocks (noutput)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: noutput

  character(100) :: blocksfile
  character(100) :: visitfile
  character(50) :: cbuffer
  character(3) :: rankstr
  character(4) :: noutstr
  character(1) :: slash
  character(1) :: lf = char(10)
  real :: xb, yb, zb
  real(kind=4) :: xx, yy, zz
  integer :: a, b, c, i, j, k, l, p, unitout, bID, pID, istat, nblocks
  integer :: nb, ncells, npoints, lvl, pas, pas0
  integer :: counter
!  real :: nextrep

  ! Check trailing slash in data dir path
  l = len_trim(datadir)
  if (datadir(l:l).ne.'/') then
    slash = '/'
  else
    slash = ''
  end if
  
  ! Write .visit metadata file - master only
  if (rank.eq.master) then
    write(noutstr,'(i4.4)') noutput
    write(visitfile,'(a)') trim(datadir) // trim(slash) // "master.visit"
    if ((.not.dowarm).and.(noutput.eq.0)) then
      open (unit=7, file=visitfile, status='replace', iostat=istat)
    else
      open (unit=7, file=visitfile, status='unknown', position='append', &
      iostat=istat)
    end if
    if (istat.ne.0) then
      write(logu,'(a,a,a)') "Could not open the file '", trim(visitfile), "' !"
      write(logu,'(a,a,a)') "Does the datadir '", trim(datadir), "' exist?"
      close(unit=7)
      call clean_abort (ERROR_OUTPUT_FILE_OPEN)
    else
      write(logu,*)
      write(logu,'(1x,a,a)') "Writing VisIt metadata file ", trim(visitfile)
      if (noutput.eq.0) then
        write (7, '(a,i5)') "!NBLOCKS ", nProcs
      end if
      do pID=0,nProcs-1
        blocksfile = blockstpl
        write(rankstr,'(I3.3)') pid
        call replace (blocksfile, 'XXX', rankstr)
        call replace (blocksfile, 'YYYY', noutstr)
        write(blocksfile,'(A)') trim(blocksfile) // ".vtk"
        write(7, '(A)') trim(blocksfile)
      end do
      close(unit=7)
    end if
  end if     

  ! Generate data filename based on templates
  write(rankstr,'(i3.3)') rank
  write(noutstr,'(i4.4)') noutput
  blocksfile = blockstpl
  call replace (blocksfile, 'XXX', rankstr)
  call replace (blocksfile, 'YYYY', noutstr)
  write(blocksfile,'(a)') trim(datadir) // trim(slash) // trim(blocksfile) // ".vtk"

  ! Write data file - all processes 
  unitout = 10 + rank
  open (unit=unitout, file=blocksfile, status='unknown', access='stream', &
        convert="BIG_ENDIAN", iostat=istat)
  if (istat.ne.0) then
    write(logu,'(a,a,a)') "Could not open the file '", trim(blocksfile), "' !"
    write(logu,'(a,a,a)') "Does the datadir '", trim(datadir), "' exist?"
    close(unitout)
    call clean_abort (ERROR_OUTPUT_FILE_OPEN)
  end if

  write(logu,*) ""
  write(logu,'(1x,a,i0,a,a,a)') "Rank ", rank, " writing BLOCKS data file ", trim(blocksfile), " ..."
     
  ! VTK Header
  write(unitout) "# vtk DataFile Version 4.2", lf
  write(unitout) "Output from Walicxe3D", lf
  write(unitout) "BINARY", lf
  write(unitout) "DATASET UNSTRUCTURED_GRID", lf

  ! Calculate total number of cells and vertices
  npoints = nbLocal*8*(ncells_x)*(ncells_y)*(ncells_z)
  ncells = nbLocal*(ncells_x)*(ncells_y)*(ncells_z)

  ! VTK cell vertices list
  ! DEBUG
!  nextrep = 0.01
  write(logu,'(1x,a)') "Writing VTK cell vertices list ..."
  ! DEBUG
  write(cbuffer,'(a,1x,i12,1x,a)') "POINTS", npoints, "float"
  write(unitout) trim(cbuffer), lf
  counter = 0
  nblocks = 0
  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then
      call meshlevel(bID, lvl)
      call getRefCorner(bID, xb, yb, zb)
      do k=1,ncells_z
        do j=1,ncells_y
          do i=1,ncells_x
            ! For each cell, dump all 8 vertices
!              write(logu,*) ''
!              write(logu,'(a,i3,a,i3,i3,i3)') "Block", bID, ", cell", i, j, k
!              write(logu,'(a,f8.3,f8.3,f8.3)') "Ref corner: ", xb, yb, zb
!              write(logu,'(a,f8.3,f8.3,f8.3)') "Cell coords: ", xb+i*dx(lvl), yb+j*dy(lvl), zb+k*dz(lvl)
            do a=0,1
              do b=0,1
                do c=0,1
                  if (units_type.eq.PHYS_UNITS) then
                    xx = real( ( xb + (i+c-1)*dx(lvl) ) * l_sc, 4 )
                    yy = real( ( yb + (j+b-1)*dy(lvl) ) * l_sc, 4 )
                    zz = real( ( zb + (k+a-1)*dz(lvl) ) * l_sc, 4 )
                  elseif (units_type.eq.CODE_UNITS) then
                    xx = real( ( xb + (i+c-1)*dx(lvl) ), 4 )
                    yy = real( ( yb + (j+b-1)*dy(lvl) ), 4 )
                    zz = real( ( zb + (k+a-1)*dz(lvl) ), 4 )
                  end if
                  write(unitout) xx, yy, zz
                  counter = counter + 1
                  ! DEBUG
!                  if (counter*1.0/npoints.ge.nextrep) then
!                    write(logu,'(a)',advance='no') ". "
!                    nextrep = nextrep + 0.01
!                  end if
                  ! DEBUG
                end do
              end do
            end do
          end do
        end do 
      end do
      nblocks = nblocks + 1
    end if
  end do
  write(unitout) lf
  if (counter.ne.npoints) then
    write(logu,'(a,i8,a,i8)') "POINTS expected:", npoints, ", ACTUAL:", counter
    write(logu,'(a)') "***ABORTING!***"
  end if
  if (nblocks.ne.nbLocal) then
    write(logu,'(a,i8,a,i8,a)') "Number of dumped blocks (", nblocks,") does not match nbLocal (", nbLocal, ")"
    write(logu,'(a)') "***ABORTING!***"
  end if

  ! VTK cell descriptions (i.e., which vertices define each cell)
  ! DEBUG
!  nextrep = 0.01
  write(logu,'(1x,a)') "Writing VTK cell descriptions ..."
  ! DEBUG
  write(cbuffer,'(a,1x,i12,1x,i12)') "CELLS", ncells, ncells*9
  write(unitout) trim(cbuffer), lf
  p = 0
  counter = 0
  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then
      do k=1,ncells_z
        do j=1,ncells_y
          do i=1,ncells_x
            write(unitout) 8, p, p+1, p+2, p+3, p+4, p+5, p+6, p+7
            p = p + 8
            counter = counter + 9
            ! DEBUG
!            if (counter*1.0/(ncells*9).ge.nextrep) then
!              write(logu,'(a)',advance='no') ". "
!              nextrep = nextrep + 0.01
!            end if
            ! DEBUG
          end do
        end do 
      end do
    end if
  end do
  write(unitout) lf
  if (counter.ne.ncells*9) then
    write(logu,'(a,i8,a,i8)') "CELLS*9 expected:", ncells*9, ", ACTUAL:", counter
    write(logu,'(a)') "***ABORTING!***"
  end if

  ! VTK cell types - all are voxels, type=11
  ! DEBUG
!  nextrep = 0.01
  write(logu,'(a)') " Writing VTK cell types ..."
  ! DEBUG
  write(cbuffer,'(a,1x,i12)') "CELL_TYPES", ncells
  write(unitout) trim(cbuffer), lf
  counter = 0
  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then
      do k=1,ncells_z
        do j=1,ncells_y
          do i=1,ncells_x
            write(unitout) 11
            counter = counter + 1
            ! DEBUG
!            if (counter*1.0/ncells.ge.nextrep) then
!              write(logu,'(a)',advance='no') ". "
!              nextrep = nextrep + 0.01
!            end if
            ! DEBUG
          end do
        end do
      end do
    end if
  end do    
  write(unitout) lf
  if (counter.ne.ncells) then
    write(logu,'(a,i8,a,i8)') "CELL TYPES expected:", ncells*9, ", ACTUAL:", counter
    write(logu,'(a)') "***ABORTING!***"
  end if
  
  ! VTK cell data dump begins here
  write(cbuffer,'(A,I12)') "CELL_DATA", ncells
  write(unitout) trim(cbuffer), lf

  ! Density
!  nextrep = 0.01
  counter = 0
  npoints = nbLocal*ncells_x*ncells_y*ncells_z
  write(logu,'(1x,a)') "Writing Density values ..."
  write(cbuffer,'(a)')  "SCALARS rho float 1"
  write(unitout) trim(cbuffer), lf
  write(cbuffer,'(a)')  "LOOKUP_TABLE default"
  write(unitout) trim(cbuffer), lf
  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then
      do k=1,ncells_z
        do j=1,ncells_y
          do i=1,ncells_x
            if (units_type.eq.PHYS_UNITS) then
              xx = real( PRIM(nb,1,i,j,k) * d_sc, 4)
            elseif (units_type.eq.CODE_UNITS) then
              xx = real( PRIM(nb,1,i,j,k), 4)
            end if
            write(unitout) xx
            counter = counter + 1
            ! DEBUG
!            if (counter*1.0/npoints.ge.nextrep) then
!              write(logu,'(a)',advance='no') ". "
!              nextrep = nextrep + 0.01
!            end if
            ! DEBUG
          end do
        end do
      end do
    end if
  end do
  write(unitout) lf

  ! Velocity
  ! DEBUG
!  nextrep = 0.01
  counter = 0
  npoints = nbLocal*ncells_x*ncells_y*ncells_z
  write(logu,'(1x,a)') "Writing Velocity values ..."
  ! DEBUG
  write(cbuffer,'(a)') "VECTORS vel float"
  write(unitout) trim(cbuffer), lf
  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then
      do k=1,nCells_z
        do j=1,ncells_y
          do i=1,ncells_x
            if (units_type.eq.PHYS_UNITS) then
              xx = real( PRIM(nb,2,i,j,k)*v_sc, 4)
              yy = real( PRIM(nb,3,i,j,k)*v_sc, 4)
              zz = real( PRIM(nb,4,i,j,k)*v_sc, 4)
            elseif (units_type.eq.CODE_UNITS) then
              xx = real( PRIM(nb,2,i,j,k), 4)
              yy = real( PRIM(nb,3,i,j,k), 4)
              zz = real( PRIM(nb,4,i,j,k), 4)
            end if
            write(unitout) xx, yy, zz
            ! DEBUG
!            counter = counter + 1
!            if (counter*1.0/npoints.ge.nextrep) then
!              write(logu,'(a)',advance='no') ". "
!              nextrep = nextrep + 0.01
!            end if
            ! DEBUG
          end do
        end do
      end do
    end if
  end do
  write(unitout) lf    

  ! Pressure
  ! DEBUG
!  nextrep = 0.01
  counter = 0
  npoints = nbLocal*ncells_x*ncells_y*ncells_z
  write(logu,'(1x,a)') "Writing Pressure values ..."
  ! DEBUG
  write(cbuffer,'(A)')  "SCALARS P float 1"
  write(unitout) trim(cbuffer), lf
  write(cbuffer,'(a)')  "LOOKUP_TABLE default"
  write(unitout) trim(cbuffer), lf
  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then
      do k=1,nCells_z
        do j=1,ncells_y
          do i=1,ncells_x
            bID = localBlocks(nb)
            if (units_type.eq.PHYS_UNITS) then
              xx = real( PRIM(nb,5,i,j,k)*p_sc, 4)
            elseif (units_type.eq.CODE_UNITS) then
              xx = real( PRIM(nb,5,i,j,k), 4)
            end if
            write(unitout) xx
            ! DEBUG
!            counter = counter + 1
!            if (counter*1.0/npoints.ge.nextrep) then
!              write(logu,'(a)',advance='no') ". "
!              nextrep = nextrep + 0.01
!            end if
            ! DEBUG
          end do
        end do
      end do
    end if
  end do
  write(unitout) lf 

#ifdef PASBP
  ! Velocity
  ! DEBUG
!  nextrep = 0.01
  counter = 0
  npoints = nbLocal*ncells_x*ncells_y*ncells_z
  write(logu,'(1x,a)') "Writing Magnetic Field values ..."
  ! DEBUG
  write(cbuffer,'(a)') "VECTORS B float"
  write(unitout) trim(cbuffer), lf
  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then
      do k=1,nCells_z
        do j=1,ncells_y
          do i=1,ncells_x
            if (units_type.eq.PHYS_UNITS) then
              xx = PRIM(nb,6,i,j,k)*v_sc
              yy = PRIM(nb,7,i,j,k)*v_sc
              zz = PRIM(nb,8,i,j,k)*v_sc
            elseif (units_type.eq.CODE_UNITS) then
              xx = PRIM(nb,6,i,j,k)
              yy = PRIM(nb,7,i,j,k)
              zz = PRIM(nb,8,i,j,k)
            end if
            write(unitout) xx, yy, zz
            ! DEBUG
!            counter = counter + 1
!            if (counter*1.0/npoints.ge.nextrep) then
!              write(logu,'(a)',advance='no') ". "
!              nextrep = nextrep + 0.01
!            end if
            ! DEBUG
          end do
        end do
      end do
    end if
  end do
  write(unitout) lf    
#endif

!    ! Temperature
!    write(cbuffer,'(A)')  "SCALARS T float 1"
!    write(unitout) trim(cbuffer), lf
!    write(cbuffer,'(A)')  "LOOKUP_TABLE default"
!    write(unitout) trim(cbuffer), lf
!    do nb=1,nbLocal
!      do i=1,ncells_x
!        do j=1,ncells_y
!          do k=1,ncells_z
!            bID = localBlocks(nb)
!            rho = PRIM(bID,1,i,j,k)*rho_sc
!            pres = PRIM(bID,5,i,j,k)*P_sc
!            call calcTemp (rho, pres, 0.0, T)
!            x = T
!            write(unitout) x
!          end do
!        end do
!      end do
!    end do
!    write(unitout) lf   

  ! Metallicity
  if (cooling_type.eq.COOL_TABLE_METAL) then
    write(logu,'(1x,a)') "Writing metallicity values ..."
    write(cbuffer,'(a)') "SCALARS Z float 1"
    write(unitout) trim(cbuffer), lf
    write(cbuffer,'(a)') "LOOKUP_TABLE default"
    write(unitout) trim(cbuffer), lf
    do nb=1,nbMaxProc
      bID = localBlocks(nb)
      if (bID.ne.-1) then
        do k=1,ncells_y
          do j=1,ncells_y
            do i=1,ncells_x
              bID = localBlocks(nb)
              xx = real( PRIM(nb,metalpas,i,j,k)/PRIM(nb,1,i,j,k), 4)
              write(unitout) xx
            end do
          end do
        end do
      end if
    end do
    write(unitout) lf
  end if

  ! Passive scalars -- except metallicity
  if (cooling_type.eq.COOL_TABLE_METAL) then
    pas0 = 2
  else
    pas0 = 1
  end if
  do pas=pas0,npassive
    write(logu,'(1x,a,i0,a)') "Writing Passive Scalar ", pas, " values ..."
    write(cbuffer,'(a,i0,a)') "SCALARS passive", pas, " float 1"
    write(unitout) trim(cbuffer), lf
    write(cbuffer,'(a)') "LOOKUP_TABLE default"
    write(unitout) trim(cbuffer), lf
    do nb=1,nbMaxProc
      bID = localBlocks(nb)
      if (bID.ne.-1) then
        do k=1,ncells_y
          do j=1,ncells_y
            do i=1,ncells_x
              bID = localBlocks(nb)
              xx = real( PRIM(nb,firstpas+pas-1,i,j,k)/PRIM(nb,1,i,j,k), 4)
              write(unitout) xx
            end do
          end do
        end do
      end if
    end do
    write(unitout) lf
  end do

  close(unitout)
  if (nbLocal.eq.0) then
    write(logu,'(1x,a)') "No local blocks to write. Wrote VTK file with no data."
  else
    write(logu,'(1x,a,i3,a)') "Succesfully wrote ", nblocks, " local blocks"
    write(logu,'(1x,a,a)') "Wrote VTK file."
  end if

end subroutine write3DVTKBlocks

!===============================================================================

!> @brief Visit-compatible VTK grid file (shows blocks, not cells)
subroutine write3DVTKGrid (noutput)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: noutput

  character(100) :: gridfile
  character(4) :: noutstr
  character(1) :: slash
  integer :: l
 
  character(50) :: cbuffer
  character(1) :: lf = char(10)
  real :: xb, yb, zb
  real(kind=4) :: xx, yy, zz
  integer :: i, j, k, p, unitout, bID, istat, owner
  integer :: nb, ncells, npoints, lvl, counter
  real(kind=8) :: hk
  real(kind=4) :: hilb

  ! DEBUG
  integer :: maxnb, order
  ! DEBUG

  ! Check trailing slash in data dir path
  l = len_trim(datadir)
  if (datadir(l:l).ne.'/') then
    slash = '/'
  else
    slash = ''
  end if
  
  ! Generate data filename based on templates
  write(noutstr,'(I4.4)') noutput
  gridfile = gridtpl
  call replace (gridfile, 'YYYY', noutstr)
  write(gridfile,'(a)') trim(datadir) // trim(slash) // trim(gridfile) // ".vtk"

  ! Write grid file
  unitout = 10
  open (unit=unitout, file=gridfile, status='unknown', access='stream', &
        convert="BIG_ENDIAN", iostat=istat)
  if (istat.ne.0) then
    write(logu,'(a,a,a)') "Could not open the file '", trim(gridfile), "' !"
    write(logu,'(a,a,a)') "Does the datadir '", trim(datadir), "' exist?"
    close(unitout)
    call clean_abort (ERROR_OUTPUT_FILE_OPEN)
  end if
  
  write(logu,*) ""
  write(logu,'(1x,a,i0,a,a,a)') "Rank ", rank, " writing GRID meta file ", trim(gridfile), " ..."
     
  ! VTK Header
  write(unitout) "# vtk DataFile Version 4.2", lf
  write(unitout) "Output from Walicxe3D", lf
  write(unitout) "BINARY", lf
  write(unitout) "DATASET UNSTRUCTURED_GRID", lf

  ! The number of vtk "cells" is the number of blocks (all processors)
  ! and there are 8 vtk "points" per block
  npoints = nbActive*8
  ncells = nbActive

  ! VTK cell vertices list
  write(cbuffer,'(a,1x,i12,1x,a)') "POINTS", npoints, "float"
  write(unitout) trim(cbuffer), lf
  counter = 0
  do nb=1,nbMaxGlobal
    bID = globalBlocks(nb)
    if (bID.ne.-1) then
!        print*, "Writing block", bID, "to grid"  !DEBUG
      counter = counter + 1
      call meshlevel(bID, lvl)
      call getRefCorner(bID, xb, yb, zb)
!        print*, "Writing block", bID
!        print*, "Ref corner:", xb, yb, zb
      ! For each block, write its eight corners
      do k=0,1
        do j=0,1
          do i=0,1
            if (units_type.eq.PHYS_UNITS) then
              xx = real( ( xb + i*ncells_x*dx(lvl) ) * l_sc, 4)
              yy = real( ( yb + j*ncells_y*dy(lvl) ) * l_sc, 4)
              zz = real( ( zb + k*ncells_z*dz(lvl) ) * l_sc, 4)
            else if (units_type.eq.CODE_UNITS) then
              xx = real( ( xb + i*ncells_x*dx(lvl) ), 4)
              yy = real( ( yb + j*ncells_y*dy(lvl) ), 4)
              zz = real( ( zb + k*ncells_z*dz(lvl) ), 4)
            end if
            write(unitout) xx, yy, zz
          end do
        end do
      end do
    end if
  end do
  write(unitout) lf
  write(logu,'(a,i6)') " Number of global active blocks: ", counter

  ! VTK cell descriptions (i.e., which vertices define each cell)
  write(cbuffer,'(a,1x,i12,1x,i12)') "CELLS", ncells, ncells*(8+1)
  write(unitout) trim(cbuffer), lf
!    print*, "CELLS ", ncells, ncells*(8+1)    !DEBUG
  p = 0
  do nb=1,nbMaxGlobal
    bID = globalBlocks(nb)
    if (bID.ne.-1) then
      write(unitout) 8, p, p+1, p+2, p+3, p+4, p+5, p+6, p+7
      p = p + 8
    end if
  end do
  write(unitout) lf

  ! VTK cell types - all are voxels, type=11
  write(cbuffer,'(a,1x,i12)') "CELL_TYPES", ncells
  write(unitout) trim(cbuffer), lf
  do nb=1,nbMaxGlobal
    bID = globalBlocks(nb)
    if (bID.ne.-1) then
      write(unitout) 11
    end if
  end do
  write(unitout) lf

  ! VTK cell "data" begins here: ranks and bIDs
  write(cbuffer,'(A,I12)') "CELL_DATA", ncells
  write(unitout) trim(cbuffer), lf

  ! Owner's rank
  write(cbuffer,'(A)')  "SCALARS rank integer 1"
  write(unitout) trim(cbuffer), lf
  write(cbuffer,'(a)')  "LOOKUP_TABLE default"
  write(unitout) trim(cbuffer), lf
  do nb=1,nbMaxGlobal
    bID = globalBlocks(nb)
    if (bID.ne.-1) then
      call getOwner (bID, owner)
      write(unitout) owner
    end if
  end do
  write(unitout) lf

  ! Block's bID
  write(cbuffer,'(A)')  "SCALARS bID integer 1"
  write(unitout) trim(cbuffer), lf
  write(cbuffer,'(a)')  "LOOKUP_TABLE default"
  write(unitout) trim(cbuffer), lf
  do nb=1,nbMaxGlobal
    bID = globalBlocks(nb)
    if (bID.ne.-1) then
      write(unitout) bID
    end if
  end do
  write(unitout) lf

  ! Hilbert Key
  maxnb = max(nbx(maxlev), nby(maxlev), nbz(maxlev))
  order = ceiling(log(1.0*maxnb)/log(2.0))
  write(cbuffer,'(A)')  "SCALARS hkey integer 1"
  write(unitout) trim(cbuffer), lf
  write(cbuffer,'(a)')  "LOOKUP_TABLE default"
  write(unitout) trim(cbuffer), lf
  do nb=1,nbMaxGlobal
    bID = globalBlocks(nb)
    if (bID.ne.-1) then
      call getHilbertKey(bID, hk)
      hilb = real(hk, 4)
      write(unitout) hilb
    end if
  end do
  write(unitout) lf

  close(unitout)
  write(logu,'(1x,a,a)') "Wrote GRID data file."

end subroutine write3DVTKGrid

!===============================================================================

!> @brief Writes a 2D daya array into a Visit-compatible 2D binary file
!> @param outmap The 2D data array to be written
!> @param nx,by The dimensions of the array
!> @param outfname The filename for the output data file
subroutine write2DMapVTK (outmap, nx, ny, outfname)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: nx, ny
  real, intent(in) :: outmap(neqtot,nx,ny)
  character(*), intent(in) :: outfname

  integer :: i, j
  integer :: istat, npoints, unitout
  real :: temp
  real(kind=4) :: xx, yy, zz
  character(100) :: cbuffer
  character(1) :: lf = char(10)

  ! Open data file
  unitout = 10
  open (unit=unitout, file=outfname, status='unknown', access='stream', &
        convert="BIG_ENDIAN", iostat=istat)
  if (istat.ne.0) then
    write (*,'(a,a,a)') "Could not open the file '", trim(outfname), "' !"
    close (unitout)
    stop
  end if

  ! Total number of points (cells)
  npoints = nx*ny

  ! VTK Header
  write(unitout) "# vtk DataFile Version 4.2", lf
  write(unitout) "Output from Walicxe3D", lf
  write(unitout) "BINARY", lf
  write(unitout) "DATASET STRUCTURED_POINTS", lf
  write(cbuffer,'(a,i6,1x,i6,1x,i6)') "DIMENSIONS ", nx, ny, 1
  write(unitout) trim(cbuffer),lf
  write(cbuffer,'(a,f10.5,1x,f10.5,1x,f10.5)') "ORIGIN ", 0.0, 0.0, 0.0
  write(unitout) trim(cbuffer),lf
  write(cbuffer,'(a,3(es12.5,1x))') "SPACING ", dx(maxlev), dx(maxlev), 0.0
  write(unitout) trim(cbuffer),lf
  write(cbuffer,'(a,i10)') 'POINT_DATA ',npoints
  write(unitout) trim(cbuffer),lf

  ! DENSITY
  write(logu,'(1x,a)') " Writing Density values ..."
  write(cbuffer,'(a)')  "SCALARS rho float 1"
  write(unitout) trim(cbuffer), lf
  write(cbuffer,'(a)')  "LOOKUP_TABLE default"
  write(unitout) trim(cbuffer), lf
  do j=1,ny
    do i=1,nx
      xx = real( outmap(1,i,j)*d_sc, 4)
      write(unitout) xx
    end do
  end do
  write(unitout) lf

  ! PRESSURE
  write(logu,'(1x,a)') " Writing Pressure values ..."
  write(cbuffer,'(a)')  "SCALARS P float 1"
  write(unitout) trim(cbuffer), lf
  write(cbuffer,'(a)')  "LOOKUP_TABLE default"
  write(unitout) trim(cbuffer), lf
  do j=1,ny
    do i=1,nx
      xx = real( outmap(5,i,j)*p_sc, 4)
      write(unitout) xx
    end do
  end do
  write(unitout) lf

  ! TEMPERATURE
  write(logu,'(1x,a)') " Writing Temperature values ..."
  write(cbuffer,'(a)')  "SCALARS T float 1"
  write(unitout) trim(cbuffer), lf
  write(cbuffer,'(a)')  "LOOKUP_TABLE default"
  write(unitout) trim(cbuffer), lf
  do j=1,ny
    do i=1,nx
      call calcTemp (outmap(:,i,j), temp)
      xx = real( temp, 4)
      write(unitout) xx
    end do
  end do
  write(unitout) lf

  ! VELOCITY
  write(logu,'(1x,a)') " Writing Velocity components ..."
  write(cbuffer,'(a)') 'VECTORS vel float'
  write(unitout) trim(cbuffer),lf
  do j=1,ny
    do i=1,nx
      xx = real( outmap(2,i,j)*v_sc, 4)
      yy = real( outmap(3,i,j)*v_sc, 4)
      zz = real( outmap(4,i,j)*v_sc, 4)
      write(unitout) xx, yy, zz
    end do
  end do

#ifdef PASBP
  ! MAGNETIC FIELD
  write(logu,'(1x,a)') " Writing Magnetic Field components ..."
  write(cbuffer,'(a)') 'VECTORS B float'
  write(unitout) trim(cbuffer),lf
  do j=1,ny
    do i=1,nx
      xx = real( outmap(6,i,j), 4)
      yy = real( outmap(7,i,j), 4)
      zz = real( outmap(8,i,j), 4)
      write(unitout) xx, yy, zz
    end do
  end do
#endif

  ! Missing passive scalars ...

  close(unitout)
  write(logu,'(1x,a,a)') "Successfully wrote VTK file."

end subroutine write2DMapVTK

!===============================================================================
