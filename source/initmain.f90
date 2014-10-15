!===============================================================================
!> @file initmain.f90
!> @brief Basic allocations and initializations
!> @author Juan C. Toledo
!> @date 3/Jun/2011

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

!> @brief Code-wide initializations
!> @details Initializes MPI, allocates big arrays, loads cooling data,
!! reports compilation parameters
subroutine initmain ()

  use parameters
  use globals
  use tictoc
  implicit none

  integer :: nps, istat, inext
  integer :: mark, l
  logical :: existing, success
  real :: totalsize
  character(3) :: rankstr
  character(1) :: slash

  ! =================================

  ! MPI Initialization
#ifdef MPIP
  call mpi_init (ierr)
  call mpi_comm_size (mpi_comm_world, nps, ierr)
  call mpi_comm_rank (mpi_comm_world, rank, ierr)
  if (nps.ne.nProcs) then
    if (rank.eq.master) then
      write(*,'(a,i0,a)') "Actual number of processors ( ", nps, " ) does not match value specified "
      write(*,'(a,i0,a)') "in parameters.f90 ( ", nProcs, " ) !"
      write(*,'(a)') "***ABORTING***"
      write(*,*) ""
    end if
    call clean_abort (ERROR_WRONG_NPROCS)
  end if
#else
  rank = master
#endif

  ! =================================

  ! Start logging
  if (logged) then

    ! Generate logfile name
    l = len_trim(logdir)
    if (logdir(l:l).ne.'/') then
      slash = '/'
    else
      slash = ''
    end if

    write(rankstr,'(i3.3)') rank
    logfile = trim(logdir) // trim(slash) // "rank" // rankstr // ".log"

    ! Check if file already exists
    success = .false.
    inext = 1
    do while(.not.success)
      inquire (file=logfile, exist=existing) 
      if (existing) then
        ! file exists - append a number, try again
        write(logfile,'(a,a,a,a,i0,a)') trim(logdir) // trim(slash), &
          "rank", rankstr, "-", inext, ".log"
        inext = inext + 1
      else
        ! Open logfile
        logu = 10 + nProcs + rank
        open (unit=logu, file=logfile, status='unknown', position="append", &
          iostat=istat)
        success = .true.
      end if   
    end do

    if (istat.ne.0) then
      if (rank.eq.master) then
        write(*,'(a)') "Could not open the log file!"
        write(*,'(a,a)') "Tried to open: ", logfile
        write(*,'(a,a,a)') "Does the logdir '", trim(logdir), "' exist ?"
        write(*,'(1x,a)') "***ABORTING***"
      end if
      close(logu)
      call clean_abort (ERROR_NO_LOGFILE)
    end if

  else

    logu = 6

  end if
 
  ! =================================
  call tic(mark)
  write(logu,*) ""
  write(logu,'(a)') "================================================================================"
  if (logged) then
    write(logu,*) ""
    write(logu,'(1x,a,a)') "> Started logfile ", trim(logfile)
  end if
  write(logu,*) ""
  write(logu,'(1x,a)') "> Initializing ... "

  if (dowarm) then
    write(logu,'(1x,a)') "A WARM START has been scheduled"
  end if

  ! Get hostname
  call HostNm(host)

  if(logged.or.rank.eq.master) then
    write(logu,*)
    write(logu,'(1x,a)') "***************************************************"  
    write(logu,'(1x,a)') "*   __    __      _ _               _____ ____    *"  
    write(logu,'(1x,a)') "*  / / /\ \ \__ _| (_) _____  _____|___ /|  _ \   *"  
    write(logu,'(1x,a)') "*  \ \/  \/ / _` | | |/ __\ \/ / _ \ |_ \| | | |  *"  
    write(logu,'(1x,a)') "*   \  /\  / (_| | | | (__ |  |  __/___) | |_| |  *"  
    write(logu,'(1x,a)') "*    \/  \/ \__,_|_|_|\___/_/\_\___|____/|____/   *"  
    write(logu,'(1x,a)') "*                                                 *"  
    write(logu,'(1x,a)') "*         Version 1.2 ($Revision:: 76  $)         *"  
    write(logu,'(1x,a)') "*                                                 *"
#ifdef MPIP
    write(logu,'(1x,a,i3,a)') "*        Running with MPI on ", nProcs , " processors       *"
#else
    write(logu,'(1x,a)') "*               Running serially                  *"
#endif
    write(logu,'(1x,a,a,a)') "*        Hostname: ", host, "                *"
    write(logu,'(1x,a,a,a)') "*        Start: ", STAMP(), "            *"
    write(logu,'(1x,a)') "*                                                 *"
    write(logu,'(1x,a)') "***************************************************"
  end if

  call tic (start_mark)

  write(logu,*) ""
  write(logu,'(1x,a,i0,a)') "Processor ", rank, " ready."

  ! =================================

  ! Report parameters. For warm starts, critical parameters are also verified.

  write(logu,*) ""
  write(logu,'(1x,a)') "============================================"
  write(logu,'(1x,a)') " Doing basic initializations ..."
  write(logu,'(1x,a)') "============================================"

  write(logu,*) ""
  write(logu,'(1x,a,es12.5,a,f7.3,a)') "Simulation box x-size:  ", xphystot, " cm / ", xphystot/PC, " pc"
  write(logu,'(1x,a,es12.5,a,f7.3,a)') "Simulation box y-size:  ", yphystot, " cm / ", yphystot/PC, " pc"
  write(logu,'(1x,a,es12.5,a,f7.3,a)') "Simulation box z-size:  ", zphystot, " cm / ", zphystot/PC, " pc"

  write(logu,*) ""
  if (mesh_method.eq.MESH_AUTO) then
    write(logu,'(1x,a)') "> Mesh initialization method: Automatic"
    write(logu,'(1x,a,i0)') "Max-level cells along x:  ", p_maxcells_x
    write(logu,'(1x,a,i0)') "Max-level cells along y:  ", p_maxcells_y
    write(logu,'(1x,a,i0)') "Max-level cells along z:  ", p_maxcells_z
  else if (mesh_method.eq.MESH_MANUAL) then
    write(logu,'(1x,a)') "> Mesh initialization method: Manual"
    write(logu,'(1x,a,i0)') " Number of root blocks along x:  ", p_nbrootx
    write(logu,'(1x,a,i0)') " Number of root blocks along y:  ", p_nbrooty
    write(logu,'(1x,a,i0)') " Number of root blocks along z:  ", p_nbrootz
    write(logu,'(1x,a,i0)') " Number of refinement levels:  ", p_maxlev
  end if

  write(logu,*) ""
  write(logu,'(1x,a,f7.1,a)') "Allowed RAM per process: ", RAM_per_proc, " MB"

  write(logu,*) ""
  write(logu,'(1x,a)') "> Mesh parameters"
  write(logu,'(1x,a,i0)') "Max blocks per processor:  ", nbMaxProc
  write(logu,'(1x,a,i0)') "Cells per block along x:   ", ncells_x
  write(logu,'(1x,a,i0)') "Cells per block along y:   ", ncells_y
  write(logu,'(1x,a,i0)') "Cells per block along z:   ", ncells_z
  write(logu,'(1x,a,i0)') "Ghost cell layer depth:    ", nghost
  write(logu,'(1x,a,f5.2)') "Refinement Threshold:  ", refineThres
  write(logu,'(1x,a,f5.2)') "Coarsening Threshold:  ", coarseThres

  write(logu,*) ""
  write(logu,'(1x,a)') "> Boundary Conditions"
  write(logu,'(1x,a,a)') "Left:   " , trim(bcname(bc_left))
  write(logu,'(1x,a,a)') "Right:  ", trim(bcname(bc_right))
  write(logu,'(1x,a,a)') "Front:  ", trim(bcname(bc_front))
  write(logu,'(1x,a,a)') "Back:   ", trim(bcname(bc_back))
  write(logu,'(1x,a,a)') "Bottom: ", trim(bcname(bc_bottom))
  write(logu,'(1x,a,a)') "Top:    ", trim(bcname(bc_top))

  write(logu,*) ""
  write(logu,'(1x,a)') "> Data Output"
  write(logu,'(1x,a,a)') "Data directory:     ", datadir
  write(logu,'(1x,a,a)') "Logging directory:  ", logdir
  write(logu,'(1x,a,a)') "Datafile template:  ", blockstpl
  write(logu,'(1x,a,a)') "Gridfile template:  ", gridtpl
  write(logu,'(1x,a,a)') "Statefile template: ", statetpl
  if (output_bin) then
    write(logu,'(1x,a)') "Output in internal format is ON"
  else
    write(logu,'(1x,a)') "Output in internal format is OFF"
  end if
  if (output_vtk) then
    write(logu,'(1x,a)') "Output in VisIt format is ON"
  else
    write(logu,'(1x,a)') "Output in VisIt format is OFF"
  end if
  if (units_type.eq.CODE_UNITS) write(logu,'(1x,a)') "Data dumped in Code Units"
  if (units_type.eq.PHYS_UNITS) write(logu,'(1x,a)') "Data dumped in Physical Units"

#ifdef PASBP
  write(logu,*) ""
  write(logu,'(1x,a)') "> Passive magnetic field ENABLED"
#endif

  write(logu,*) ""
  write(logu,'(1x,a)') "> Hydro Solver"
  write(logu,'(1x,a,a)') "Type: ", trim(solvername(solver_type))
  if ((solver_type.eq.SOLVER_HLL).or.(solver_type.eq.SOLVER_HLLC)) then
    ! Check that two ghost cells are used for second-order solvers
    if (nghost.ne.2) then
      write(logu,*) "This solver requires TWO ghost cells!"
      write(logu,*) "Modify parameters.f90"
      write(logu,*) "***ABORTING***"
    end if    
    write(logu,'(1x,a,a)') "Limiter: ", trim(limitername(limiter_type))
  end if
  write(logu,'(1x,a,i0)') "Hydro equations:  ", neqhydro
  write(logu,'(1x,a,i0)') "MHD equations:    ", neqmhd
  write(logu,'(1x,a,i0)') "Passive scalars:  ", npassive
  write(logu,'(1x,a,i0)') "Total equations:  ", neqtot
  write(logu,'(1x,a,f6.3)') "CFL parameter:  ", CFL
  write(logu,'(1x,a,f6.3)') "Artificial viscosity:  ", visc_eta

  write(logu,*) ""
  write(logu,'(1x,a)') "> Gas Parameters"
  write(logu,'(1x,a,f6.3)') "gamma = ", gamma
  write(logu,'(1x,a,f6.3)') "mu0   = ", mu0
  write(logu,'(1x,a,f6.3)') "mui   = ", mui
  write(logu,'(1x,a,f8.1)') "ion_thres = ", ion_thres

  ! Report unit scalings
  write(logu,*) ""
  write(logu,'(1x,a)') "> Unit scalings (1 code unit = ?)"
  write(logu,'(1x,a,es12.5,a)') "Length:   ", l_sc, " cm"
  write(logu,'(1x,a,es12.5,a)') "Density:  ", d_sc, " g cm^-3"
  write(logu,'(1x,a,es12.5,a)') "Velocity: ", v_sc, " cm s^-1"
  write(logu,'(1x,a,es12.5,a)') "Pressure: ", p_sc, " erg cm^-3"
  write(logu,'(1x,a,es12.5,a)') "Time:     ", t_sc, " s"

  ! Radiative cooling
  write(logu,*) ""
  if (cooling_type.eq.COOL_NONE) then
    write(logu,'(1x,a)') "> Radiative cooling is OFF"  
  else
    write(logu,'(1x,a)') "> Radiative cooling is ON"
    write(logu,'(2x,a,a)') "Cooling Table: ", trim(cooling_file)
    if (cooling_type.eq.COOL_TABLE) then
      call loadcooldata ()
    else if (cooling_type.eq.COOL_TABLE_METAL) then
      if (npassive.lt.1) then
        write(logu,'(1x,a)') "At least one passive scalar is needed for &
          &metallicity dependent cooling!"
        write(logu,'(1x,a)') "Set npassive to at least 1 in parameters.f90"
        write(logu,'(1x,a)') "***ABORTING***"
        call clean_abort(ERROR_NOT_ENOUGH_PASSIVES)
      end if
      call loadcooldata_metal ()
    end if
  end if

  ! Allocate memory and initialize big data arrays
  write(logu,*) ""
  write(logu,'(1x,a)') "> Allocating memory for big arrays ..."
  write(logu,*) ""

  ! Sanity check: abort if RAM will be insufficient to allocate big arrays
!  totMem = neqtot*(nxmax-nxmin+1)*(nymax-nymin+1)*(nymax-nymin+1)*(nbmaxProc*3+6)
!  if (mpi_real_kind.eq.MPI_DOUBLE_PRECISION) then
!    totMem = totMem*8
!  else if (mpi_real_kind.eq.MPI_REAL) then
!    totMem = totMem*4
!  end if
!  totMem = totMem/1024.0/1024.0
!  if (totMem.gt.maxMemProc) then
!    write(logu,'(a)') "Required memory would surpass maximum allowed memory!"
!    write(logu,'(a,f6.1,a)') "Required: ", totMem, " MB"
!    write(logu,'(a,f6.1,a)') "Allowed:  ", maxMemProc, " MB"
!    write(logu,'(a,f6.1,a)') "Either decrease nbMaxProc or increase maxMemProc in parameters.f90"
!    write(logu,*) "***ABORTING***"    
!    call clean_abort (ERROR_INSUFFICIENT_RAM)
!  else
!    write(logu,'(1x,a,f6.1,a)') "Estimated required memory for big arrays: ", totMem, " MB"
!  end if

  ! Big data arrays

  allocate( U(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax), stat=istat)
  U(:,:,:,:,:) = 0.0
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for U array!"
    write(logu,*) "***ABORTING***"    
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if

  allocate( UP(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax), stat=istat)
  UP(:,:,:,:,:) = 0.0
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for UP array!"
    write(logu,*) "***ABORTING***"    
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if

  allocate( PRIM(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax), stat=istat)
  PRIM(:,:,:,:,:) = 0.0
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for P array!"
    write(logu,*) "***ABORTING***"     
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if

  ! Allocate memory and initialize auxiliary data arrays
  
  allocate( F(neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax), stat=istat)
  F(:,:,:,:) = 0.0
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for F array!"
    write(logu,*) "***ABORTING***"     
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if
  
  allocate( G(neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax), stat=istat)
  G(:,:,:,:) = 0.0
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for G array!"
    write(logu,*) "***ABORTING***"     
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if
  
  allocate( H(neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax), stat=istat)
  H(:,:,:,:) = 0.0
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for H array!"
    write(logu,*) "***ABORTING***"     
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if
  
  allocate( FC(neqtot, nxmin:nxmax, nxmin:nxmax, nxmin:nxmax), stat=istat)
  FC(:,:,:,:) = 0.0
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for FC array!"
    write(logu,*) "***ABORTING***"     
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if
  
  allocate( GC(neqtot, nxmin:nxmax, nxmin:nxmax, nxmin:nxmax), stat=istat)
  GC(:,:,:,:) = 0.0
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for GC array!"
    write(logu,*) "***ABORTING***"     
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if
  
  allocate( HC(neqtot, nxmin:nxmax, nxmin:nxmax, nxmin:nxmax), stat=istat)
  HC(:,:,:,:) = 0.0
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for HC array!"
    write(logu,*) "***ABORTING***"     
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if

  ! Allocate and initialize global and local block registries

  allocate ( globalBlocks(nbMaxGlobal), stat=istat )
  globalBlocks(:) = -1
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for global block registry!"
    write(logu,*) "***ABORTING***"     
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if

  allocate ( localBlocks(nbMaxProc), stat=istat )
  localBlocks(:) = -1
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for local block registry!"
    write(logu,*) "***ABORTING***"     
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if

  ! Allocate and initialize global and local refinement flags lists

  allocate ( refineFlags(nbMaxGlobal), stat=istat )
  refineFlags(:) = -1
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for refinement flags!"
    write(logu,*) "***ABORTING***"     
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if

  ! Calculate global block index range
  nbmin = rank*nbMaxProc + 1
  nbmax = rank*nbMaxProc + nbMaxProc

  ! Report allocation results
  write(logu,'(1x,a)') "Successfully allocated big arrays."
  write(logu,'(1x,a,f6.1,a)') "U:      ", sizeof(U)/(1024.*1024.), " MB" 
  write(logu,'(1x,a,f6.1,a)') "UP:     ", sizeof(UP)/(1024.*1024.), " MB"
  write(logu,'(1x,a,f6.1,a)') "PRIM:   ", sizeof(PRIM)/(1024.*1024.), " MB"
  write(logu,'(1x,a,f6.1,a)') "Fluxes: ", (sizeof(F)+sizeof(G)+sizeof(H)+&
    sizeof(FC)+sizeof(GC)+sizeof(HC))/(1024.*1024.), " MB"
  totalsize = (sizeof(U)+sizeof(UP)+sizeof(PRIM)+sizeof(F)+sizeof(G)+&
    sizeof(H)+sizeof(FC)+sizeof(GC)+sizeof(HC))/(1024.*1024.)
  write(logu,'(1x,a,f7.1,a,f7.1,a)') "Total: ", totalsize, " MB / ", &
    RAM_per_proc, "MB"

  ! =================================
  
  ! Initialize simulation state - (warm start does this later)
  time = 0.0
  it = 0
  nextout = 0

  ! =================================

  write(logu,'(1x,a,a)') ""
  write(logu,'(1x,a,a)') "> Performed initializations and allocated big arrays in ", nicetoc(mark)
  write(logu,*) ""  

  ! Barrier
  call mpi_barrier (mpi_comm_world, ierr)

end subroutine initmain
