!===============================================================================
!> @file coldens.f90
!> @brief Line of sight integration facility
!> @author Juan C. Toledo
!> @date 2/Aug/2012

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

!> @brief Integrates 3D data along a line-of-sight (e.g., column density)
!> @details Version 2.0
program coldens

  implicit none

  ! Named constants -- don't modify

  integer, parameter :: AXIS_NONE = 0
  integer, parameter :: AXIS_X = 1
  integer, parameter :: AXIS_Y = 2
  integer, parameter :: AXIS_Z = 3

  integer, parameter :: INT_COLDENS = 1
  integer, parameter :: INT_XRAY = 2
  integer, parameter :: INT_SYNC = 3

  integer, parameter :: LOAD_BLOCK = 1
  integer, parameter :: LOAD_ALL = 2

  integer, parameter :: ACCEL_ISO = 0
  integer, parameter :: ACCEL_PAR = 1
  integer, parameter :: ACCEL_PER = 2

  ! Constants (7 significant digits)
  real, parameter :: PC  = 3.085680E+18
  real, parameter :: AMU = 1.660539E-24
  real, parameter :: KB  = 1.380658E-16
  real, parameter :: PI  = 3.141593
  
  ! ============================================================================
  ! PROGRAM CONFIGURATION
  ! ============================================================================

  ! Run Parameters
  integer, parameter :: outmin = 0     ! First output number to process
  integer, parameter :: outmax = 0     ! Last output number to process

  ! Integration Type
  ! Currently supported options:
  !  INT_COLDENS: column density integration
  !  INT_XRAY: thermal xray emission
  !  INT_SYNC: synchrotron emission
  integer, parameter :: int_type = INT_SYNC

  ! Data load method
  ! Currently supported options:
  !  LOAD_BLOCK: load and process one block at a time.
  !  LOAD_ALL: copy data from all blocks into one max-resolution buffer array
  ! Synchrotron emission requires LOAD_ALL; it will ignore this option.
  integer, parameter :: load_type = LOAD_BLOCK

  ! Rotation parameters
  ! The following parameters allow the user to pre-apply up to four
  ! 3D rotations to the the scene before doing the integration along the
  ! line of sight.
  ! The rot_center parameters specify the center of rotation, given in
  ! cell integer coordinates (same for all rotations)
  ! The rot_axis parameters must be one of the constants
  !   AXIS_X, AXIS_Y, AXIS_Z or AXIS_NONE (if no rotation is wanted)
  ! Positive rotation angles around a given axis are *counter-clockwise*
  ! when the arrow end of that axis points at the observer. 
  ! The scene is always projected along the Z-axis, with the tip of that
  ! axis pointed towards the observer, the +x axis pointing right, and
  ! the +y axis pointing up. All rotations are performed around these
  ! original axes (they remain *fixed* as the data rotates).
  !
  !        y
  !        ^
  !        |   Plane of the sky
  !        |                                
  !        #----- > x
  !       /
  !      /
  !     v
  !    z = line of sight
  !
  ! In other words, you might think of rotations around the AXIS_X, AXIS_Y
  ! and AXIS_Z axes as follows:
  ! Rotations around AXIS_Z are rotations around the line_of_sight
  ! Rotations around AXIS_X are rotations around a horizontal axis in the
  !   plane of the sky which points to the right
  ! Rotations around AXIS_Y are rotations around a vertical axis in the
  !   plane of the sky which points up
  ! Common shortcuts:
  ! If you want to project along the x-axis, set the following rotation:
  ! rot1_axis = AXIS_Y, rot1_angle = +90.0
  ! If you want to project along the y-axis, set the following rotation:
  ! rot1_axis = AXIS_X, rot1_angle = -90.0
  real, parameter :: rot_center_x = 256
  real, parameter :: rot_center_y = 256
  real, parameter :: rot_center_z = 256
  integer, parameter :: rot1_axis  = AXIS_X
  real,    parameter :: rot1_angle = -90.0
  integer, parameter :: rot2_axis  = AXIS_NONE
  real,    parameter :: rot2_angle = 0.0
  integer, parameter :: rot3_axis  = AXIS_NONE
  real,    parameter :: rot3_angle = 0.0
  integer, parameter :: rot4_axis  = AXIS_NONE
  real,    parameter :: rot4_angle = 0.0

  ! Output map size
  integer, parameter :: mapcells_x = 512
  integer, parameter :: mapcells_y = 512

  ! Output format
  logical, parameter :: output_bin = .true.    ! Direct binary output
  logical, parameter :: output_rg  = .true.    ! Graphic output

  ! Filenames
  character(*), parameter :: datadir    = "./BIN/"            ! Path to data dir
  character(*), parameter :: blockstpl  = "BlocksXXX.YYYY"    ! Data file template
  character(*), parameter :: blockstpl2 = "BlocksXXX.YYYY+100" ! Data file template (alt)
  character(*), parameter :: gridtpl    = "Grid.YYYY"         ! Grid file template
  character(*), parameter :: statetpl   = "State.YYYY"        ! State file template
!  character(*), parameter :: outtpl     = "Xray0.15_8kev2.6e22X.YYYY"  ! Output file template
  character(*), parameter :: outtpl     = "SyncZZZ+AA_new.YYYY"  ! Output file template
!  character(*), parameter :: outtpl     = "TestZZZ"  ! Output file template

  ! Physical box sizes (cgs)
  real, parameter :: xsize = 20*PC
  real, parameter :: ysize = 20*PC
  real, parameter :: zsize = 20*PC

  ! Mesh parameters
  integer, parameter :: nbrootx = 1
  integer, parameter :: nbrooty = 1
  integer, parameter :: nbrootz = 1
  integer, parameter :: maxlev = 6
  integer, parameter :: ncells_x = 16
  integer, parameter :: ncells_y = 16
  integer, parameter :: ncells_z = 16

  ! Simulation parameters
  integer, parameter :: nprocs = 24
  integer, parameter :: neqtot = 9

  ! Unit scalings
  real, parameter :: l_sc = 1.0*PC          !< length scale (cm)
  real, parameter :: d_sc = 1.0*1.3*AMU     !< density scale (g cm^-3)
  real, parameter :: v_sc = 3.2649e05       !< velocity scale (cm s^-1)
  real, parameter :: p_sc = d_sc*v_sc**2
  real, parameter :: e_sc = p_sc
  real, parameter :: t_sc = l_sc/v_sc

  ! ====================================================
  ! Xray calculation parameters
  ! ====================================================

  ! Xray emissivity coefficients
  character(*), parameter :: xray_coefs = '/home/claudio/coefs/coef0.15_8kev2.6e22.dat'

  ! Gas parameters: mean molecular masses for ionized and neutral gas,
  ! ionization threshold and gamma/heat capacity
  real, parameter :: mui = 1.3/2.1
  real, parameter :: mu0 = 1.3
  real, parameter :: ion_thres = 5000
  real, parameter :: gamma = 5.0/3.0
  real, parameter :: cv = 1.0/(gamma-1.0)

  ! ====================================================
  ! Synchrotron calculation parameters
  ! ====================================================

  ! Spectral index (S propto nu^(-alpha))
  real, parameter :: alpha = 0.6

  ! Entropy change threshold
  real, parameter :: deltaS = 0.0

  integer, parameter :: accel = ACCEL_PAR
  
  ! ============================================================================
  !                    NO NEED TO MODIFY BELOW THIS POINT
  ! ============================================================================

  integer :: nout, p, unitin, istat
  integer :: bID, blocksread, nb, nblocks
  integer :: bx, by, bz, i, j, k, i1, j1, k1, ip, jp, kp, ilev
  integer :: i2, j2, k2, i_off, j_off, k_off
  integer :: totcells_x, totcells_y, totcells_z
  integer :: mesh(7), numcoefs, counter
  real :: a, b, c, sigma, Kfact, prog, sumvalue
  real :: dx(maxlev), pvars(neqtot), uvars(neqtot)
  real :: line_of_sight(3)
  character(256) :: filename  

  real, allocatable :: block(:,:,:,:)
  real, allocatable :: prim(:,:,:,:)
  real, allocatable :: prim_old(:,:,:,:)
  real, allocatable :: outmap(:,:)
  real, allocatable :: xraytable(:,:)  

  ! ====================================================
  
  write(*,*) "Allocating data arrays ..."

  ! Allocate output map
  allocate ( outmap(mapcells_x, mapcells_y) )
  write(*,'(1x,a,f6.2,a)') "outmap: ", sizeof(outmap)/1024./1024., " MB"

  ! Allocate buffer data arrays
  if ((load_type.eq.LOAD_ALL).or.(int_type.eq.INT_SYNC)) then

    ! Allocate data array for whole simulation at max resolution
    totcells_x = nbrootx*ncells_x*2**(maxlev-1)
    totcells_y = nbrooty*ncells_y*2**(maxlev-1)
    totcells_z = nbrootz*ncells_z*2**(maxlev-1)

    allocate( prim(neqtot,totcells_x,totcells_y,totcells_z) )
    prim(:,:,:,:) = 0.0
    write(*,'(1x,a,f6.1,a)') "prim: ", sizeof(prim)/1024./1024., " MB"
    
    ! Second data array (for previous output) for synchrotron emission
    if (int_type.eq.INT_SYNC) then
      allocate( prim_old(neqtot,totcells_x,totcells_y,totcells_z) )
      prim_old(:,:,:,:) = 0.0
      write(*,'(1x,a,f6.1,a)') "prim_old: ", sizeof(prim_old)/1024./1024., " MB"
    end if

  else if (load_type.eq.LOAD_BLOCK) then

    ! Allocate data array for just one block
    allocate( block(neqtot,ncells_x,ncells_y,ncells_z) )
    block(:,:,:,:) = 0.0
    write(*,'(1x,a,f6.1,a)') "prim_old: ", sizeof(block)/1024./1024., " MB"

  else

    write(*,*) "Invalid LOAD_TYPE option:", load_type
    stop

  end if

  ! For xray emission, load coefficients from file
  if (int_type.eq.INT_XRAY) then

    write(*,'(2x,a,a,a)') "Loading xray coefficients from file ", &
      trim(xray_coefs), " ..."
    open (unit=99, file=xray_coefs, status='old', iostat=istat)
    if (istat.ne.0) then
      write(*,'(a,a,a)') "Could not open the file ", trim(xray_coefs), " !"
      write(*,*) "***ABORTING***"
      close(99)
      stop
    end if

    read(99,*) numcoefs
    allocate( xraytable(2,numcoefs) )

    do i=1,numcoefs
      read(99,*) a, b, c
      xraytable(1,i) = b
      xraytable(2,i) = c
      write(*,*) b, c
    end do

    close (unit=99)

  end if

  ! Grid spacings - assumed EQUAL for all dimensions
  do ilev=1,maxlev
    dx(ilev) = xsize/(ncells_x*nbrootx*2**(ilev-1))
  end do

  ! Compute line-of-sight vector (after rotations)
  call computeLOS (line_of_sight)

  ! Pack mesh parameters
  mesh(1) = nbrootx
  mesh(2) = nbrooty
  mesh(3) = nbrootz
  mesh(4) = maxlev
  mesh(5) = ncells_x
  mesh(6) = ncells_y
  mesh(7) = ncells_z

  ! ====================================================

  ! Iterate over all outputs
  do nout=outmin,outmax

    write(*,*) "=============================="
    write(*,'(1x,a,i0,a)') "Processing output ", nout, " ..."
    write(*,*) "=============================="
    write(*,*)

    if (int_type.eq.INT_COLDENS) then
      write(*,'(1x,a,1x)') "Integrating density and projecting ..."
    else if (int_type.eq.INT_XRAY) then
      write(*,'(1x,a,1x)') "Calculating x-ray emission and projecting ..."
    end if

    ! Delegate data processing to subroutine
    ! Synchrotron emission cannot be integrated on a block-by-block
    ! basis because it requires *two* different outputs and mesh may
    ! change in between
    if (int_type.eq.INT_SYNC) then
      call processAll ()
    else
      if (load_type.eq.LOAD_ALL) then
        call processAll ()
      else if (load_type.eq.LOAD_BLOCK) then
        call processByBlock ()
      end if
    end if

    ! Write output map to disk
    if (output_bin) then
      call genfname (0, nout, datadir, outtpl, ".bin", filename)
      call writebin (filename, mapcells_x, mapcells_y, outmap)
    end if
    if (output_rg) then
      call genfname (0, nout, datadir, outtpl, ".rg", filename)
      call writerg(filename, mapcells_x, mapcells_y, outmap)
    end if

    write(*,*) "=============================="

  end do

!===============================================================================

contains

!===============================================================================

subroutine processAll ()

  implicit none

  integer :: nextrep, idx
  character(3) :: progstr

  ! Reset arrays
  prim(:,:,:,:) = 0.0
  outmap(:,:) = 0.0
  if (int_type.eq.INT_SYNC) then
    prim_old(:,:,:,:) = 0.0
  end if
  
  ! Load data
  call readbin (nout,.false.,prim,neqtot,totcells_x,totcells_y,totcells_z)
  write(*,*) "Density:", minval(prim(1,:,:,:)), maxval(prim(1,:,:,:))
  if (int_type.eq.INT_SYNC) then
    call readbin (nout,.true.,prim_old,neqtot,totcells_x,totcells_y,totcells_z)
    write(*,*) "Density:", minval(prim_old(1,:,:,:)), maxval(prim_old(1,:,:,:))
  end if

  write(*,'(1x,a,1x)') "Calculating synchrotron emission and projecting ..."

  counter = 0  ! Progress counter
  nextrep = 2  ! Next report (percentage)
  write(*,'(3x)',advance='no')

  ! For each cell, compute desired quantity and project it into 2D map
  do i=1,totcells_x
    do j=1,totcells_y
      do k=1,totcells_z


        ! Compute physical quantity to integrate
        if (int_type.eq.INT_COLDENS) then

          ! COLUMN DENSITY CALCULATION
          sumvalue = prim(1,i,j,k)*dx(maxlev)

        else if (int_type.eq.INT_XRAY) then

          ! X-RAY CALCULATION
          call calcxray (prim(:,i,j,k), sumvalue)
          sumvalue = sumvalue*dx(maxlev)

        else if (int_type.eq.INT_SYNC) then

          ! SYNCHROTRON CALCULATION
          if ((i.eq.1).or.(i.eq.totcells_x).or.(j.eq.1).or. &
          (j.eq.totcells_y).or.(k.eq.1).or.(k.eq.totcells_z)) then
            ! Skip calculation for face cells
            sumvalue = 0.0
          else
            call calc_sigma (prim(:,i,j,k), prim_old(:,i,j,k), alpha, &
              deltaS, line_of_sight, sigma)
            call calc_K (prim, neqtot, totcells_x, totcells_y, totcells_z, &
              i, j, k, dx(maxlev), accel, Kfact)
            sumvalue = Kfact*sigma*dx(maxlev)
          end if

          if (sumvalue.ne.sumvalue) then
            write(*,'(a,i0,1x,i0,1x,i0,1x)') "NaN in sumvalue at", i, j, k
          end if

        end if

        ! Project this cell's value onto 2D output map
        call projectCell (sumvalue, i, j, k)
       
        ! Progress report
        counter = counter + 1
        prog = counter*1.0/(totcells_x*totcells_y*totcells_z)
        if (prog.ge.nextrep/100.0) then
          ! Fancy progress bar
          write(*,'(a)',advance="no") achar(13)
          write(*,'(1x,a)',advance="no") "["
          do idx=1,nextrep/2
            write(*,'(a)',advance="no") "="
          end do
          do idx=1,50-nextrep/2
            write(*,'(1x)',advance="no")
          end do
          write(progstr,'(i3)') nextrep
          write(*,'(a,a,a)',advance="no") "] ", trim(adjustl(progstr)), "%"
          write(*,*)
          !write(*,'(3a)',advance='no') char(8), char(8), char(8)
          !write(*,'(i2,a)',advance='no') int(nextrep*100), "%"
          nextrep = nextrep + 2
        end if

      end do
    end do
  end do

  write(*,'(a,a,a,a)',advance="no") achar(8), achar(8), achar(8), achar(8)
  write(*,'(1x,a)') "Done"

end subroutine processAll

!===============================================================================

subroutine processByBlock ()

  implicit none

  ! Reset arrays
  block(:,:,:,:) = 0.0
  outmap(:,:) = 0.0

  ! Read one data file at a time
  do p=0,nprocs-1

    ! Generate filename based on templates
    call genfname (p, nout, datadir, blockstpl, ".bin", filename)

    ! Open data file
    write(*,'(1x,a,a,a)') "Opening data file '", trim(filename), "' ..."
    unitin = 10 + p
    open (unit=unitin, file=filename, status='old', access='stream', iostat=istat)
    if (istat.ne.0) then
      write(*,'(a,a,a)') "Could not open the file '", trim(filename), "' !"
      write(*,'(a,a,a)') "Does the datadir '", trim(datadir), "' exist?"
      close(unitin)
      stop
    end if
    
    ! Read file header
    blocksread = 0
    read(unitin) nblocks
    write(*,'(1x,a,i0,a)') "File contains ", nblocks, " blocks."

    ! Loop over all blocks, process one block at a time
    do nb=1,nblocks

      ! Process this block's data

      read(unitin) bID
      read(unitin) block(:,:,:,:)
      blocksread = blocksread + 1
!      write(*,'(i0,a)',advance='no') bID, " ... "
      !write(*,*) minval(block(1,:,:,:)), maxval(block(1,:,:,:))
      
      call meshlevel (bID, mesh, ilev)
      call bcoords(bID, mesh, bx, by, bz)

!        write(*,'(1x,a,i4)') "level=", ilev
!        write(*,'(1x,a,i4,i4,i4)') "bcoords=", bx, by, bz

      ! For every cell ...
      do i=1,ncells_x
        do j=1,ncells_y
          do k=1,ncells_z

            ! Compute and de-scale primitives on this cell
            uvars = block(:,i,j,k)
            call flow2prim (uvars, pvars)
            pvars(1) = pvars(1)*d_sc
            pvars(2) = pvars(2)*v_sc
            pvars(3) = pvars(3)*v_sc
            pvars(4) = pvars(4)*v_sc
            pvars(5) = pvars(5)*p_sc

            ! 1) Compute physical quantity to be integrated

            if (int_type.eq.INT_COLDENS) then

              ! COLUMN DENSITY CALCULATION
              sumvalue = pvars(1)*dx(maxlev)

            else if (int_type.eq.INT_XRAY) then

              ! X-RAY CALCULATION
              call calcxray (pvars, sumvalue)
              sumvalue = sumvalue*dx(maxlev)

            else if (int_type.eq.INT_SYNC) then

              ! SYNCHROTRON CALCULATION CANNOT BE PERFORMED BLOCK-BY-BLOCK
              write(*,'(1x,a)') "Can't perform synchrotron calculation on a block-by-block basis"
              stop

            end if

            ! 2) Project cells into 2D output map

            call absCoords (bID,i,j,k,mesh,i1,j1,k1)
            do i_off=0,2**(maxlev-ilev)-1
              do j_off=0,2**(maxlev-ilev)-1
                do k_off=0,2**(maxlev-ilev)-1
                  i2 = i1 + i_off
                  j2 = j1 + j_off
                  k2 = k1 + k_off
                  call projectCell (sumvalue, i2, j2, k2)
                end do
              end do
            end do

          end do
        end do
      end do

!      write(*,'(4a)',advance='no') char(8), char(8), char(8), char(8)

    end do

    write(*,'(1x,a,i0,a)') "Successfully processed ", blocksread, " blocks."
    write(*,*) ""

    close (unitin)

  end do

end subroutine processByBlock

!===============================================================================

! Reads flow data from a Walicxe3D .bin data ouput and returns an array of primitives
subroutine readbin (nout,use_alt,prim,neq,nx,ny,nz)

  implicit none

  integer, intent(in) :: nout
  logical, intent(in) :: use_alt
  integer, intent(in) :: neq
  integer, intent(in) :: nx
  integer, intent(in) :: ny
  integer, intent(in) :: nz
  real, intent(out) :: prim(neq,nx,ny,nz)

  integer :: xoffset, yoffset, zoffset

  real :: block(neqtot,ncells_x,ncells_y,ncells_z)
  real :: uvars(neq), pvars(neq)

  ! Read data files
  do p=0,nprocs-1

    ! Generate filename
    if (use_alt) then
      call genfname (p, nout, datadir, blockstpl2, ".bin", filename)
    else
      call genfname (p, nout, datadir, blockstpl, ".bin", filename)
    end if

    ! Open data file
    write(*,'(1x,a,a,a)') "Openinig data file '", trim(filename), "' ..."
    unitin = 10 + p
    open (unit=unitin, file=filename, status='old', access='stream', iostat=istat)
    if (istat.ne.0) then
      write(*,'(a,a,a)') "Could not open the file '", trim(filename), "' !"
      write(*,'(a,a,a)') "Does the datadir '", trim(datadir), "' exist?"
      close(unitin)
      stop
    end if

    ! Read a block
    blocksread = 0
    read(unitin) nblocks
    write(*,'(1x,a,i0,a)') "File contains ", nblocks, " blocks. Reading blocks ..."
    do nb=1,nblocks

      read(unitin) bID
      write(*,'(1x,i0)',advance='no') bID
      read(unitin) block(:,:,:,:)
!      write(*,*) minval(block(1,:,:,:)), maxval(block(1,:,:,:))
      blocksread = blocksread + 1

      call meshlevel (bID, mesh, ilev)
      call bcoords (bID, mesh, bx, by, bz)

      ! Read and copy all cells of this cell into buffer array
      do i=1,ncells_x
        do j=1,ncells_y
          do k=1,ncells_z

            ! Calculate absolute cell coords in max-resolution grid
            i1 = (bx-1)*ncells_x*2**(maxlev-ilev) + (i-1)*2**(maxlev-ilev) + 1
            j1 = (by-1)*ncells_y*2**(maxlev-ilev) + (j-1)*2**(maxlev-ilev) + 1
            k1 = (bz-1)*ncells_z*2**(maxlev-ilev) + (k-1)*2**(maxlev-ilev) + 1

            ! Compute and de-scale primitives
            ! (Magnetic field isn't scaled at the moment)
            uvars = block(:,i,j,k)
            call flow2prim (uvars,pvars)
            pvars(1) = pvars(1)*d_sc
            pvars(2) = pvars(2)*v_sc
            pvars(3) = pvars(3)*v_sc
            pvars(4) = pvars(4)*v_sc
            pvars(5) = pvars(5)*p_sc

            ! Add the cell's data to the corresponding max-resolution cell(s)
            do xoffset=0,2**(maxlev-ilev)-1
              do yoffset=0,2**(maxlev-ilev)-1
                do zoffset=0,2**(maxlev-ilev)-1
                  ip = i1 + xoffset
                  jp = j1 + yoffset
                  kp = k1 + zoffset
                  prim(:,ip,jp,kp) = prim(:,ip,jp,kp) + pvars(:)
                end do
              end do
            end do
            
          end do
        end do
      end do

    end do
    write(*,'(/1x,a,i0,a)') "Successfully read ", blocksread, " blocks."
    write(*,*) ""

    close (unitin)

  end do

end subroutine readbin

!===============================================================================

subroutine projectCell (cellvalue, i, j, k)

  implicit none

  real, intent(in) :: cellvalue
  integer, intent(in) :: i, j, k

  real :: x, y, z, xp, yp, zp, xm, ym, x2, y2, a, b, weight
  integer :: im, jm, ii, jj, i2, j2, off_x, off_y

  ! Absolute mesh coordinates of cell's center in max-resolution grid
  x = i + 0.5
  y = j + 0.5
  z = k + 0.5

  ! Apply all rotations
  call transform (x, y, z, xp, yp, zp)

  ! Project transformed point into output map along the Z-axis
  ! The transformed point center is contained in the output map
  ! cell with center at (xm,ym)
  im = int(xp)
  jm = int(yp)
  xm = im + 0.5
  ym = jm + 0.5

  ! DEBUG
!  write(*,'(i6,i6,i6,a,i6,i6)') i,j,k,"->",im,jm

  ! Cell offsets for neighbor cells on output map
  if (xp.le.xm) then
    off_x = -1
  else
    off_x = +1
  end if
  if (yp.le.ym) then
    off_y = -1
  else
    off_y = +1
  end if

  ! Distribute value among 4 map cells using weights
  do ii=0,1
    do jj=0,1
    
      i2 = im + ii*off_x
      j2 = jm + jj*off_y

      if ((i2.ge.1).and.(i2.le.mapcells_x).and.\
      (j2.ge.1).and.(j2.le.mapcells_y)) then
        x2 = i2 + 0.5
        y2 = j2 + 0.5
        a = 1.0 - abs(xp-x2)
        b = 1.0 - abs(yp-y2)
        weight = a*b
        outmap(i2,j2) = outmap(i2,j2) + weight*cellvalue
      end if
      
    end do
  end do

end subroutine

!===============================================================================

subroutine computeLOS (los)

  implicit none

  real, intent(out) :: los(3)

  real :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4 

  ! Rotate z unit vector
  call rotatePoint (0.0, 0.0, 1.0, rot1_axis, rot1_angle, x1, y1, z1)
  call rotatePoint (x1, y1, z1, rot2_axis, rot2_angle, x2, y2, z2)
  call rotatePoint (x2, y2, z2, rot3_axis, rot3_angle, x3, y3, z3)
  call rotatePoint (x3, y3, z3, rot4_axis, rot4_angle, x4, y4, z4)

  los(1) = x4
  los(2) = y4
  los(3) = z4

end subroutine computeLOS

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

end subroutine absCoords

!===============================================================================

subroutine flow2prim (uvars, pvars)

  implicit none

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

subroutine calcxray (pvars,xrayj)

  implicit none

  real, intent(in) :: pvars(neqtot)
  real, intent(out) :: xrayj

  integer :: i
  real :: mintemp, maxtemp, temp, T1, T2, C1, C2, CX

  ! Calculate temperature (in cgs)
  temp = pvars(5)/pvars(1)*(mu0*AMU/KB)
  if (temp.gt.ion_thres) then
    temp = pvars(5)/pvars(1)*(mui*AMU/KB)
  end if

  ! Interpolate emission coefficient
  mintemp = xraytable(1,1)
  maxtemp = xraytable(1,numcoefs)
  if (temp.lt.mintemp) then
    CX = 0.0
  else if (temp.gt.maxtemp) then
    CX = xraytable(2,numcoefs)*(temp/maxtemp)**0.5
!    write(*,*) "T>1e8, CX=", CX
  else
    do i=2,numcoefs
      if (xraytable(1,i).gt.temp) then
        T1 = xraytable(1,i-1)
        T2 = xraytable(1,i)
        C1 = xraytable(2,i-1)
        C2 = xraytable(2,i)
        CX = C1 + (C2-C1)/(T2-T1)*(temp-T1)
        exit
      end if
    end do
!    write(*,*) "T=", temp, "CX=", CX
  end if

  ! Calculate X-ray emission, j = n^2 * Lambda(T)
  if (temp<ion_thres) then
    xrayj=CX*(pvars(1)/mu0/AMU)**2
  else
    xrayj=CX*(pvars(1)/mui/AMU)**2
  end if

  return

end subroutine calcxray

!===============================================================================

! Computes "raw" emissivity (in arbitrary units)
! Raw Emisivity (sigma) is calculated as:
!  sigma ~ rho^(1-2*alpha) * pres^(2*alpha) * (B_sky)^(alpha+1)
! where B_sky = B*sin(psi) is the component of the magnetic field
! in the plane of the sky (psi=angle between B and line of sight)
! Inputs:
!   prim: real vector of primitives for this output
!   prim_old: real vector of primitives for previous output
!   alpha: spectral index (S~nu**(-alpha))
!   deltaS: fractional entropy change threshold for nonzero emission
!   los: projection axis, one of AXIS_X, AXIS_Y or AXIS_Z
! Outputs:
!   sigma: calculated raw emissivity
subroutine calc_sigma (prim,prim_old,alpha,deltaS,los,sigma)

  implicit none

  real, intent(in) :: prim(neqtot)
  real, intent(in) :: prim_old(neqtot)
  real, intent(in) :: alpha
  real, intent(in) :: deltaS
  real, intent(in) :: los(3)
  real, intent(out) :: sigma

  real :: dens, dens_old, pres, pres_old, S_new, S_old, delta, mu, mgam
  real :: al1, al2, al3, B, Bx, By, Bz, nx, ny, nz, cospsi, sinpsi

  ! Gas parameters
  mgam = -5.0/3.0
  mu = 1.3

  dens = prim(1)
  pres = prim(5)
  dens_old = prim_old(1)
  pres_old = prim_old(5)
#ifdef PASBP
  Bx = prim(7)
  By = prim(8)
  Bz = prim(9)
#endif

  ! Obtain entropy change
  dens = dens/(mu*AMU)
  dens_old = dens_old/(mu*AMU)
  pres = pres/1.0e-8
  pres_old = pres_old/1.0e-8
  S_new = pres*dens**mgam
  S_old = pres_old*dens_old**mgam
  delta = abs((S_new-S_old)/S_old)
 
  ! Compute sigma only when entropy change is large enough
  if (delta.gt.deltaS) then

    ! Line-of-sight unit vector
    nx = los(1)
    ny = los(2)
    nz = los(3)
    if (abs(nx*nx+ny*ny+nz*nz-1.0).gt.1e-5) then
      write(*,*) "Passed line of sight is not a unit vector!"
      write(*,*) "LOS =", los
      write(*,*) "***ABORTING***"
      stop
    end if

    al1 = 1-2.0*alpha
    al2 = 2.0*alpha
    al3 = alpha+1.0

    B = sqrt(Bx*Bx+By*By+Bz*Bz)
    cospsi = (Bx*nx+By*ny+Bz*nz)/B
    sinpsi = sqrt(1.0-cospsi*cospsi)

    sigma = (dens**al1)*(pres**al2)*(B*sinpsi)**al3

  else
  
    sigma = 0.0
    
  end if

end subroutine calc_sigma

!===============================================================================

! Angle-dependent amplification factor for synchrotron emission
subroutine calc_K (prim, neq, nx, ny, nz, i, j, k, dx, accel, Kfact)

  implicit none

  integer, intent(in) :: neq
  integer, intent(in) :: nx
  integer, intent(in) :: ny
  integer, intent(in) :: nz
  real, intent(in) :: prim(neq,nx,ny,nz)
  integer, intent(in) :: i
  integer, intent(in) :: j
  integer, intent(in) :: k
  real, intent(in) :: dx
  integer, intent(in) :: accel
  real, intent(out) :: Kfact
  
  real :: gradpx, gradpy, gradpz, gradp2, bx, by, bz, b2, bdotgradp
  real :: c2phibn, c2phibs

  ! Calculate obliquity dependence factor of shock acceleration:
  !   C2PHIBN = cos^2(phibn)
  ! where phibn is the angle between the shock normal and
  ! the magnetic field. The shock normal is approximated using
  ! the direction of the local pressure gradient

  gradpx = (prim(5,i+1,j,k)-prim(5,i-1,j,k))/(2*dx)
  gradpy = (prim(5,i,j+1,k)-prim(5,i,j-1,k))/(2*dx)
  gradpz = (prim(5,i,j,k+1)-prim(5,i,j,k-1))/(2*dx)
  bx = prim(6,i,j,k)
  by = prim(7,i,j,k)
  bz = prim(8,i,j,k)
  b2 = bx**2 + by**2 + bz**2
  gradp2 = gradpx**2 + gradpy**2 + gradpz**2
  bdotgradp = bx*gradpx + by*gradpy + bz*gradpz
  
  if ((b2.lt.1E-60).or.(gradp2.lt.1E-60)) then
    c2phibn = 0.0
  else
    c2phibn = bdotgradp**2/(b2*gradp2)
  end if

!  print*, c2phibn

  ! Derive corresponding angle-dependent factor,
  !   K = constant     ... (isotropic)
  !   K = cos^2(phibs) ... (quasi-parallel)
  !   K = sin^2(phibs) ... (quasi-perpendicular)
  ! where
  !   cos^2(phibs) = cos^2(phibn)/(16-15*cos^2(phibn)),          
  ! according to particle acceleration method.

  c2phibs = c2phibn/(16-15*c2phibn)
  if (accel.eq.ACCEL_ISO) then
    Kfact = 1.0
  else if (accel.eq.ACCEL_PAR) then
    Kfact = c2phibs
  else if (accel.eq.ACCEL_PER) then
    Kfact = 1.0 - c2phibs
  else
    Kfact = 0.0
  end if

end subroutine calc_K

!===============================================================================

subroutine writebin (fname,nx,ny,outmap)

  character(*), intent(in) :: fname
  integer, intent(in) :: nx, ny
  real, intent(in) :: outmap(nx,ny)

  write(*,*) ""
  write(*,'(1x,a,a)') "Writing output file ", trim(fname)
  write(*,'(1x,a,i4,i4)') "Map dimensions:", nx, ny
  write(*,'(1x,a,i0)') "Map size: ", sizeof(outmap)
  write(*,'(1x,a,es10.3,es10.3)') "Range of values:", minval(outmap), maxval(outmap)
  write(*,*) ""
  open (unit=10, file=fname, status='replace', form='unformatted', iostat=istat)
  write (10) outmap(:,:)
  close (10)

end subroutine writebin

!===============================================================================

subroutine writerg (fname,nx,ny,outmap)

  character(*), intent(in) :: fname
  integer, intent(in) :: nx, ny
  real, intent(in) :: outmap(nx,ny)

  integer :: i, j

  write(*,*) ""
  write(*,'(1x,a,a)') "Writing output file ", trim(fname)
  write(*,'(1x,a,i4,i4)') "Map dimensions:", nx, ny
  write(*,'(1x,a,es10.3,es10.3)') "Range of values:", minval(outmap), maxval(outmap)
  write(*,*) ""
  open(unit=10, file=fname, status='replace', iostat=istat)
  write(10,*) nx,1,0,1
  write(10,'(a)') ' '
  write(10,*) ny,1,0,1
  write(10,'(a)') ' '
  write(10,'(10Z8.8)') ((real(outmap(i,j),4),i=1,nx),j=1,ny)
  close(10)

end subroutine writerg

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
  character(4) :: noutstr
  character(3) :: rankstr
  character(2) :: rotstr
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

  if (accel.eq.ACCEL_ISO) call replace(filename, 'ZZZ', 'ISO')
  if (accel.eq.ACCEL_PAR) call replace(filename, 'ZZZ', 'PAR')
  if (accel.eq.ACCEL_PER) call replace(filename, 'ZZZ', 'PER')
!  if (line_of_sight.eq.AXIS_X) call replace(filename, 'D', 'X')
!  if (line_of_sight.eq.AXIS_Y) call replace(filename, 'D', 'Y')
!  if (line_of_sight.eq.AXIS_Z) call replace(filename, 'D', 'Z')
!  write(rotstr,'(I2.2)') int(rot_beta)
  write(rotstr,'(I2.2)') int(rot2_angle)
  call replace(filename, 'AA', rotstr)

  write(filename,'(a)') trim(dir) // trim(slash) // trim(filename) // trim(ext)

end subroutine genfname

!===============================================================================

! Performs the rotations asked by the user on a point
subroutine transform (x, y, z, xp, yp, zp)

  implicit none
  
  real, intent(in) :: x, y, z
  real, intent(out) :: xp, yp, zp

  real :: x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4

  ! Translate so that rotation center is origin coord system
  call translatePoint (x, y, z,  &
    rot_center_x, rot_center_y, rot_center_z, x0, y0, z0)

  ! Apply rotations
  call rotatePoint (x0, y0, z0, rot1_axis, rot1_angle, x1, y1, z1)
  call rotatePoint (x1, y1, z1, rot2_axis, rot2_angle, x2, y2, z2)
  call rotatePoint (x2, y2, z2, rot3_axis, rot3_angle, x3, y3, z3)
  call rotatePoint (x3, y3, z3, rot4_axis, rot4_angle, x4, y4, z4)

  ! Translate back to original coordinates
  call translatePoint (x4, y4, z4, &
    -rot_center_x, -rot_center_y, -rot_center_z, xp, yp, zp)

end subroutine transform

!===============================================================================

! Rotate a point about a cartesian axis
subroutine rotatePoint (x, y, z, axis, theta, xp, yp, zp)

  implicit none
  real, intent(in) :: x, y, z
  integer, intent(in) :: axis
  real, intent(in) :: theta
  real, intent(out) :: xp, yp, zp

  real :: theta_rad

  ! Don't do anything for null rotations
  if ((axis.eq.AXIS_NONE).or.(theta.eq.0.0)) then
    xp = x
    yp = y
    zp = z
    return
  end if

  ! Apply rotation matrix
  theta_rad = theta*PI/180.0
  if (axis.eq.AXIS_X) then
    xp = x
    yp = cos(theta_rad)*y - sin(theta_rad)*z
    zp = sin(theta_rad)*y + cos(theta_rad)*z
    return
  else if (axis.eq.AXIS_Y) then
    xp = cos(theta_rad)*x + sin(theta_rad)*z
    yp = y
    zp = -sin(theta_rad)*x + cos(theta_rad)*z
    return
  else if (axis.eq.AXIS_Z) then
    xp = cos(theta_rad)*x - sin(theta_rad)*y
    yp = sin(theta_rad)*x + cos(theta_rad)*y
    zp = z
    return
  else
    write(*,*) "Invalid rotation axis!"
    write(*,*) "axis=", axis
    write(*,*) "***ABORTING***"
    stop
  end if

end subroutine rotatePoint

!===============================================================================

! Translate a point, applying shifts along +x, +y and +z
subroutine translatePoint (x, y, z, lx, ly, lz, xp, yp, zp)

  implicit none
  real, intent(in) :: x, y, z
  real, intent(in) :: lx, ly, lz
  real, intent(out) :: xp, yp, zp

  xp = x - lx
  yp = y - ly
  zp = z - lz

end subroutine translatePoint

!===============================================================================

end program coldens
