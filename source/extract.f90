!===============================================================================
!> @file extract.f90
!> @brief Data extractor for Walicxe3D ouput files
!> @author Juan C. Toledo
!> @date 18/Feb/2013

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
program extract

  use constants
  implicit none
  
  ! ============================================================================
  ! PROGRAM CONFIGURATION
  ! ============================================================================

  ! Output range to process
  integer, parameter :: noutmin = 0
  integer, parameter :: noutmax = 0

  ! Axis and location of cut
  ! cut_axis must be one of AXIS_X, AXIS_Y, AXIS_Z
  ! cut_location must be given in physical units (cgs)
  integer, parameter :: cut_axis = AXIS_Y
  real, parameter :: cut_location = 0.5*PC

  ! Path to data dir and name of parameters file (must be in the datadir)
  character(*), parameter :: datadir =      ! Path to data directory
  character(*), parameter :: params_file = "Params.dat"    ! Simulation parameters file

  ! Output format
  logical, parameter :: output_vtk = .false.   ! VTK output
  logical, parameter :: output_bin = .true.   ! Direct binary output

  ! ============================================================================
  !                    NO NEED TO MODIFY BELOW THIS POINT
  ! ============================================================================

  integer :: nprocs, neqtot, npassive, firstpas
  real :: xphystot, yphystot, zphystot
  integer :: nbrootx, nbrooty, nbrootz, maxlev, ncells_x, ncells_y,  ncells_z
  real :: gamma, mu0, mui, ion_thres, CV
  real :: l_sc, d_sc, v_sc, p_sc, e_sc, t_sc
  
  integer :: ilev, bID, blocksused, istat, nb, p
  integer :: i, j, k, i1, j1, k1, ip, jp, i_off, j_off, i2, j2
  integer :: mesh(7), unitin, nblocks, plane, nout
  integer :: nxmap, nymap, nx, ny, cell_count
  character(256) :: filename, filepath

  real, allocatable :: block(:,:,:,:)
  real, allocatable :: outmap(:,:,:)
  real, allocatable :: dx(:), pvars(:), uvars(:)
  
  ! Filename templates
  character(*), parameter :: tpl_blocks = "BlocksXXX.YYYY"  ! Data files template
  character(*), parameter :: tpl_state = "State.YYYY.dat"   ! State files template
  character(*), parameter :: tpl_outmap = "CutD.YYYY"  ! Output file template

  ! Uncomment this section, and comment parameters noutmin and noutmax above,
  ! to instead read these two params as command-line arguments
!  integer :: noutmin, noutmax
!  character(len=32) :: arg
!  call get_command_argument(1, arg)
!  read(arg, '(I10)') noutmin
!  call get_command_argument(2, arg)
!  read(arg, '(I10)') noutmax

  ! ====================================================

  ! Read simulation parameters from file
  write(filepath, "(a)") datadir // "/" // params_file
  open (unit=10, file=filepath, status='old', access='stream', form='formatted', iostat=istat)
  read(10,*) nprocs, neqtot, npassive, firstpas
  read(10,*) xphystot, yphystot, zphystot
  read(10,*) nbrootx, nbrooty, nbrootz, maxlev, ncells_x, ncells_y, ncells_z
  read(10,*) gamma, mu0, mui, ion_thres
  read(10,*) l_sc, d_sc, v_sc
  close(10)
  print*, nprocs, neqtot, npassive, firstpas
  print*, xphystot, yphystot, zphystot
  print*, nbrootx, nbrooty, nbrootz, maxlev, ncells_x, ncells_y, ncells_z
  print*, gamma, mu0, mui, ion_thres
  print*, l_sc, d_sc, v_sc

  ! Derived parameters
  CV = 1.0/(gamma-1.0)
  p_sc = d_sc*v_sc**2
  e_sc = p_sc
  t_sc = l_sc/v_sc
      
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

  ! Allocate data array for one block
  allocate( block(neqtot,ncells_x,ncells_y,ncells_z) )

  ! Allocate other variables
  allocate(dx(maxlev))
  allocate(pvars(neqtot))
  allocate(uvars(neqtot))

  ! Grid spacings - assumed EQUAL for all dimensions
  do ilev=1,maxlev
    dx(ilev) = xphystot/(ncells_x*nbrootx*2**(ilev-1))
  end do

  ! Pack mesh parameters
  mesh(1) = nbrootx
  mesh(2) = nbrooty
  mesh(3) = nbrootz
  mesh(4) = maxlev
  mesh(5) = ncells_x
  mesh(6) = ncells_y
  mesh(7) = ncells_z

  ! ====================================================

  do nout=noutmin,noutmax

  ! Reset arrays
  block(:,:,:,:) = 0.0
  outmap(:,:,:) = 0.0

  cell_count = 0

  ! Process blocks data files from all processes
  do p=0,nprocs-1

    ! Generate filename
    call genfname (p, nout, datadir, tpl_blocks, ".bin", filename)

    ! Open data file
    write(*,'(1x,a,a,a)') "Openinig data file '", trim(filename), "' ..."
    unitin = 10 + p
    open (unit=unitin,file=filename,status='old',access='stream',iostat=istat)
    if (istat.ne.0) then
      write(*,'(a,a,a)') "Could not open the file '", trim(filename), "' !"
      write(*,'(a,a,a)') "Does the datadir '", trim(datadir), "' exist?"
      close(unitin)
      stop
    end if

    ! Read file header
    blocksused = 0
    read(unitin) nblocks
    write(*,'(1x,a,i0,a)') "Processing ", nblocks, " blocks ..."

    ! Loop over all blocks, process one block at a time
    do nb=1,nblocks

      ! Read bID and determine if block intersecs cut plane
      read (unitin) bID
      read (unitin) block    
      call meshlevel (bID, mesh, ilev)        
      call getCellPlane (bID, mesh, plane)
      
      ! If block intersects, extract cell plane
      if (plane.ne.-1) then
  
        blocksused = blocksused + 1
!        write(*,*) "Cutplane=", plane    
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
            call absCoords (bID,i,j,k,mesh,i1,j1,k1)

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
!            write(*,'(a)') "Output map cells:"
            do i_off=0,2**(maxlev-ilev)-1
              do j_off=0,2**(maxlev-ilev)-1

                ! Calculate outmap map coords
                i2 = i1 + i_off
                j2 = j1 + j_off

!                write(*,'(i0,1x,i0)') i2, j2
                ! Calculate and de-scale primitives
                uvars(:) = block(:,i,j,k)
                call flow2prim (neqtot, uvars, pvars)
                pvars(1) = pvars(1)*d_sc
                pvars(2) = pvars(2)*v_sc
                pvars(3) = pvars(3)*v_sc
                pvars(4) = pvars(4)*v_sc
                pvars(5) = pvars(5)*p_sc
                outmap(:,i2,j2) = pvars(:)
                cell_count = cell_count + 1

              end do
            end do

          end do
        end do

      end if     
    end do

    write(*,'(1x,a,i0,a)') "Extracted data from ", blocksused, " blocks."

  end do

  write(*,*) ""
  write(*,*) "Done extracting 2D cut."
  write(*,*) "Total cells copied:", cell_count

  if (cell_count < nxmap*nymap) then
    write(*,*) "Warning: less cells copied than size of map!"
  end if

  write(*,*) ""
  write(*,'(1x,a)') "Range of values:"
  write(*,'(1x,a,es10.3,1x,es10.3)') "Density: ", &
    minval(outmap(1,:,:))*d_sc, maxval(outmap(1,:,:))*d_sc
  write(*,'(1x,a,es10.3,1x,es10.3)') "Pressure: ", &
    minval(outmap(5,:,:))*p_sc, maxval(outmap(5,:,:))*p_sc

  ! Write output map to disk
  if (output_vtk) then
    write(*,*) ""
    call genfname (0, nout, datadir, tpl_outmap, ".vtk", filename)
    write(*,*) "Writing output map to VTK file ", trim(filename)
    call write2DVTK (outmap, neqtot, nxmap, nymap, filename)
  end if

  if (output_bin) then
    write(*,*) ""
    call genfname (0, nout, datadir, tpl_outmap, ".bin", filename)
    write(*,*) "Writing output map to BIN file ", trim(filename)
    call write2Dbin (outmap, neqtot, nxmap, nymap, filename)
  end if

  write(*,*) ""

  end do

!===============================================================================

contains

!===============================================================================

subroutine getCellPlane (bID, mesh, plane)

  implicit none
  
  integer, intent(in) :: bID
  integer, intent(in) :: mesh(7)  
  integer, intent(out) :: plane

  real :: xl, xh, yl, yh, zl, zh, bl, bh
  integer :: ilev

  call bounds (bID, mesh, xl, xh, yl, yh, zl, zh)
  call meshlevel (bID, mesh, ilev)

!  write(*,*) bID, ilev
!  write(*,*) xl/PC, xh/PC
!  write(*,*) yl/PC, yh/PC
!  write(*,*) zl/PC, zh/PC

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
  if ((cut_location.ge.bl).and.(cut_location.lt.bh)) then
    plane = int((cut_location-bl)/(dx(ilev))) + 1
  else
    plane = -1
  end if

  return

end subroutine getCellPlane

!===============================================================================

subroutine bounds(bID, mesh, xl, xh, yl, yh, zl, zh)

  implicit none
  integer, intent(in) :: bID
  integer, intent(in) :: mesh(7)
  real, intent(out) :: xl, xh, yl, yh, zl, zh

  integer :: ilev, bx, by, bz

  call meshlevel (bID, mesh, ilev)
  call bcoords(bID, mesh, bx, by, bz)

  xl = (bx-1)*xphystot/(nbrootx*2**(ilev-1))
  xh = bx*xphystot/(nbrootx*2**(ilev-1))
  yl = (by-1)*yphystot/(nbrooty*2**(ilev-1))
  yh = by*yphystot/(nbrooty*2**(ilev-1))
  zl = (bz-1)*zphystot/(nbrootz*2**(ilev-1))
  zh = bz*zphystot/(nbrootz*2**(ilev-1))

  return

end subroutine bounds

!===============================================================================

subroutine bcoords(bID, mesh, ix, iy, iz)

  implicit none

  integer, intent(in) :: bID
  integer, intent(in) :: mesh(7)
  integer, intent(out) :: ix, iy, iz

  integer :: ilev, nx, ny, nz, localID, boffset

  call meshlevel (bID, mesh, ilev)
  call levelOffset (ilev, mesh, boffset)

  nx = mesh(1)*2**(ilev-1)
  ny = mesh(2)*2**(ilev-1)
  nz = mesh(3)*2**(ilev-1)

  localID = bID - boffset
  ix = mod(localID,nx)
  if (ix.eq.0) ix=nx
  iy = mod(ceiling(localID*1.0/(nx)),ny)
  if (iy.eq.0) iy=ny
  iz = mod(ceiling(localID*1.0/(nx*ny)),nz)
  if (iz.eq.0) iz=nz

end subroutine bcoords

!===============================================================================

subroutine meshlevel(bID, mesh, level)

  implicit none

  integer, intent(in) :: bID
  integer, intent(in) :: mesh(7)
  integer, intent(out) :: level

  integer :: minID, maxID

  if (bID.eq.-1) then
    level = -1
    return
  end if

  maxID = 0
  do level=1,mesh(4)
    minID = maxID + 1
    maxID = maxID + mesh(1)*mesh(2)*mesh(3)*8**(level-1)
    if ((bID.ge.minID).and.(bID.le.maxID)) then
      return
    end if
  end do

end subroutine meshlevel

!===============================================================================

subroutine levelOffset(level, mesh, boffset)

  implicit none

  integer, intent(in) :: level
  integer, intent(in) :: mesh(7)
  integer, intent(out) :: boffset

  integer :: ilev

  boffset = 0
  do ilev=1,level-1
    boffset = boffset + mesh(1)*mesh(2)*mesh(3)*8**(ilev-1)
  end do

end subroutine levelOffset

!===============================================================================

subroutine genfname (rank, nout, dir, template, ext, filename)

  implicit none

  integer, intent(in) :: rank
  integer, intent(in) :: nout
  character(*), intent(in) :: dir
  character(*), intent(in) :: template
  character(*), intent(in) :: ext
  character(256), intent(out) :: filename

  character(1) :: slash
  character(3) :: rankstr
  character(4) :: noutstr
  integer :: l

  l = len_trim(dir)
  if (datadir(l:l).ne.'/') then
    slash = '/'
  else
    slash = ''
  end if
  write(rankstr,'(I3.3)') rank
  write(noutstr,'(I4.4)') nout
  filename = template
  call replace (filename, 'XXX', rankstr)
  call replace (filename, 'YYYY', noutstr)

  if (cut_axis.eq.AXIS_X) call replace(filename, 'D', 'X')
  if (cut_axis.eq.AXIS_Y) call replace(filename, 'D', 'Y')
  if (cut_axis.eq.AXIS_Z) call replace(filename, 'D', 'Z')

  write(filename,'(a)') trim(dir) // trim(slash) // trim(filename) // trim(ext)

end subroutine genfname

!===============================================================================

! Calculate absolute cell coords in max-resolution grid
subroutine absCoords(bID,i,j,k,mesh,i1,j1,k1)

  implicit none

  integer, intent(in) :: bID
  integer, intent(in) :: i, j, k
  integer, intent(in) :: mesh(7)
  integer, intent(out) :: i1, j1, k1

  integer :: ilev, bx, by, bz
  integer :: maxlev, ncells_x, ncells_y, ncells_z

  maxlev = mesh(4)
  ncells_x = mesh(5)
  ncells_y = mesh(6)
  ncells_z = mesh(7)

  ! Obtain block coords
  call meshlevel (bID, mesh, ilev)
  call bcoords(bID, mesh, bx, by, bz)

  i1 = (bx-1)*ncells_x*2**(maxlev-ilev) + (i-1)*2**(maxlev-ilev) + 1
  j1 = (by-1)*ncells_y*2**(maxlev-ilev) + (j-1)*2**(maxlev-ilev) + 1
  k1 = (bz-1)*ncells_z*2**(maxlev-ilev) + (k-1)*2**(maxlev-ilev) + 1
  
  return

end subroutine absCoords

!===============================================================================

!> @brief Visit-compatible binary VTK data files (2D)
subroutine write2DVTK (outmap, neqtot, nx, ny, outfname)

  implicit none

  integer, intent(in) :: neqtot, nx, ny
  real, intent(in) :: outmap(neqtot,nx,ny)
  character(*), intent(in) :: outfname

  integer :: istat, npoints, unitout, ipas
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
  write(cbuffer,'(a,3(es12.5,1x))') "SPACING ", dx(maxlev), dx(maxlev), &
    dx(maxlev)
  write(unitout) trim(cbuffer),lf
  write(cbuffer,'(a,i10)') 'POINT_DATA ',npoints
  write(unitout) trim(cbuffer),lf

  ! DENSITY
  write(*,'(1x,a)') " Writing Density values ..."
  write(cbuffer,'(a)')  "SCALARS rho float 1"
  write(unitout) trim(cbuffer), lf
  write(cbuffer,'(a)')  "LOOKUP_TABLE default"
  write(unitout) trim(cbuffer), lf
  do j=1,ny
    do i=1,nx
      xx = outmap(1,i,j)
      write(unitout) xx
    end do
  end do
  write(unitout) lf  

  ! PRESSURE
  write(*,'(1x,a)') " Writing Pressure values ..."
  write(cbuffer,'(a)')  "SCALARS P float 1"
  write(unitout) trim(cbuffer), lf
  write(cbuffer,'(a)')  "LOOKUP_TABLE default"
  write(unitout) trim(cbuffer), lf
  do j=1,ny
    do i=1,nx
      xx = outmap(5,i,j)
      write(unitout) xx
    end do
  end do
  write(unitout) lf  

  ! TEMPERATURE
  write(*,'(1x,a)') " Writing Temperature values ..."
  write(cbuffer,'(a)')  "SCALARS T float 1"
  write(unitout) trim(cbuffer), lf
  write(cbuffer,'(a)')  "LOOKUP_TABLE default"
  write(unitout) trim(cbuffer), lf
  do j=1,ny
    do i=1,nx
      call calcTemp (neqtot, outmap(:,i,j), temp)
      xx = temp
      write(unitout) xx
    end do
  end do
  write(unitout) lf

  ! VELOCITY
  write(*,'(1x,a)') " Writing Velocity components ..."
  write(cbuffer,'(a)') 'VECTORS vel float'
  write(unitout) trim(cbuffer),lf
  do j=1,ny
    do i=1,nx
      xx = outmap(2,i,j)
      yy = outmap(3,i,j)
      zz = outmap(4,i,j)
      if (cut_axis.eq.AXIS_X) then
        write(unitout) yy, zz, xx 
      else if (cut_axis.eq.AXIS_Y) then
        write(unitout) xx, zz, yy
      else if (cut_axis.eq.AXIS_Z) then
        write(unitout) xx, yy, zz
      end if
    end do
  end do

  ! PASSIVE SCALARS
  if (npassive.ge.1) then
    do ipas=1,npassive
      write(*,'(1x,a,i0,a)') " Writing Passive Scalar ",ipas," ..."
      write(cbuffer,'(a,i0,a)')  "SCALARS pas",ipas," float 1"
      write(unitout) trim(cbuffer), lf
      write(cbuffer,'(a)')  "LOOKUP_TABLE default"
      write(unitout) trim(cbuffer), lf
      do j=1,ny
        do i=1,nx
          xx = outmap(firstpas+ipas-1,i,j)
          write(unitout) xx
        end do
      end do
      write(unitout) lf
    end do
  end if

#ifdef PASBP
  ! MAGNETIC FIELD
  write(*,'(1x,a)') " Writing Magnetic Field components ..."
  write(cbuffer,'(a)') 'VECTORS B float'
  write(unitout) trim(cbuffer),lf
  do j=1,ny
    do i=1,nx
      xx = outmap(6,i,j)
      yy = outmap(7,i,j)
      zz = outmap(8,i,j)
      if (cut_axis.eq.AXIS_X) then
        write(unitout) yy, zz, xx 
      else if (cut_axis.eq.AXIS_Y) then
        write(unitout) xx, zz, yy
      else if (cut_axis.eq.AXIS_Z) then
        write(unitout) xx, yy, zz
      end if
    end do
  end do
#endif

  close(unitout)
  write(*,'(1x,a,a)') "Successfully wrote VTK file."

end subroutine write2DVTK

!===============================================================================

subroutine write2Dbin (outmap, neqtot, nx, ny, outfname)

  integer, intent(in) :: neqtot, nx, ny
  real, intent(in) :: outmap(neqtot,nx,ny)
  character(*), intent(in) :: outfname

  integer :: unitout

  write(*,'(1x,a,a)') "Writing BIN output to file ", trim(outfname)
  write(*,'(2x,i0,a)') size(outmap), " values in array"

  unitout = 10
  open (unit=unitout,file=outfname,status='replace',access='stream',&
       iostat=istat)

  write(unitout) outmap

  close(unitout)
  write(*,'(1x,a,a)') "Successfully wrote BIN file."

end subroutine write2Dbin

!===============================================================================

! Calculates temperature from primitives, in CGS
subroutine calcTemp (neqtot, pvars, temp)

  implicit none

  integer, intent(in) :: neqtot
  real, intent(in) :: pvars(neqtot)
  real, intent(out) :: temp

!  temp = pvars(5)/pvars(1)*(mu0*AMU*p_sc/d_sc/KB)

!  if (temp.gt.ion_thres) then
!    temp = pvars(5)/pvars(1)*(mui*AMU*p_sc/d_sc/KB)
!  end if
  temp = pvars(5)/pvars(1)*(mu0*AMU/KB)

  if (temp.gt.ion_thres) then
    temp = pvars(5)/pvars(1)*(mui*AMU/KB)
  end if

  return

end subroutine calcTemp

!===============================================================================

subroutine flow2prim (neqtot, uvars, pvars)

  implicit none

  integer, intent(in) :: neqtot
  real, intent(in) :: uvars(neqtot)
  real, intent(out) :: pvars(neqtot)

  real :: rhov2

  pvars(1) = uvars(1)
  pvars(2) = uvars(2)/uvars(1)
  pvars(3) = uvars(3)/uvars(1)
  pvars(4) = uvars(4)/uvars(1)

  rhov2 = (uvars(2)**2 + uvars(3)**2 + uvars(4)**2)/uvars(1)
  pvars(5) = (uvars(5)-0.5*rhov2)/CV

  ! Floor on pressure
  if (pvars(5).lt.1.0e-30) then
!    write(logu,*) "PRESSURE FLOOR APPLIED!"
!    write(logu,*) "u(5)=", uvars(5)
!    write(logu,*) "p(5)=", pvars(5)
    pvars(5) = 1.0e-30
  end if

  ! Floor on density
  if (pvars(1).lt.1.0e-40) then
    pvars(1) = 1.0e-40
  end if

  if (neqtot.gt.5) then
    pvars(6:neqtot) = uvars(6:neqtot)
  end if

  return

end subroutine flow2prim

!===============================================================================

end program extract

