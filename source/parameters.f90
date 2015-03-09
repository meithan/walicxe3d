!===============================================================================
!> @file parameters.f90
!> @brief Code parameters module
!> @author Juan C. Toledo
!> @date 2/Jun/2011

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
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.

!===============================================================================

!> @brief Definition of code parameters
!> @details This file declares most global runtime parameters (constants) used
!! by the code, most of them user-tweakable. It is included by most subroutines.
module parameters

  use constants
#ifdef MPIP
  use mpi
#endif

  implicit none

  ! ============================================
  ! Execution parameters
  ! ============================================

  real, parameter :: tfin = 1000 * YR    !< Final integration time (s)
  real, parameter :: dtout = 100 * YR     !< Time between data dumps (s)

  !> Perform warm start?
  logical, parameter :: dowarm = .false.
  !> State file to use for warm start
  character(*), parameter :: warm_file = ""

  !> Number of MPI processes to launch
  integer, parameter :: nProcs = 4

  !> Available memory (RAM) *per process*, in MB
  ! This will determine the number of blocks allocated by the code
  ! Note that 1 MB = 1024 kB = 1024*1024 or 2^20 bytes
  real, parameter :: RAM_per_proc = 1024
 
  ! ============================================
  ! Adaptive Mesh parameters
  ! ============================================

  ! Simulation domain physical size (all in cgs)
  real, parameter :: xphystot = 10.0 * PC      !< Physical domain size along x
  real, parameter :: yphystot = 10.0 * PC      !< Physical domain size along y
  real, parameter :: zphystot = 10.0 * PC      !< Physical domain size along z

  ! Mesh Geometry
  
  !> Method used to specify mesh geometry
  ! There are currently two ways to specify the mesh geometry:
  ! 1) MESH_AUTO: the user specifies max-level resolution.
  !    This will automatically calculate the number of root blocks
  !    and number of refinement levels 
  ! 2) MESH_MANUAL: the user specifies the number and geometry of the
  !    root blocks and the number of refinement levels.
  !    The max-level resolution will then be given by:
  !    dx(max) = phys_size / (ncells_x*2^(maxlev-1))
  ! In both cases, the requested grid geometry must match the aspect 
  ! ratio of the physical box dimensions.
  integer, parameter :: mesh_method = MESH_AUTO

  ! Only one of the following group of parameters will be used:

  ! 1) MESH_AUTO
  ! Specify the maximum number of cells desired at the highest refinement
  ! level
  ! > MUST BE POWERS OF TWO! <
  integer, parameter :: p_maxcells_x = 128
  integer, parameter :: p_maxcells_y = 128
  integer, parameter :: p_maxcells_z = 128

  ! -- OR --

  ! 2) MESH_MANUAL
  ! Specify number of root blocks along each dimension and number of
  ! refinement levels
  integer, parameter :: p_nbrootx = 1
  integer, parameter :: p_nbrooty = 1
  integer, parameter :: p_nbrootz = 1
  integer, parameter :: p_maxlev = 5

  ! Other Mesh Parameters

  !> Number of cells per block (per dimension)
  ! This specifies the number of cells per block along each dimension
  ! (i.e., along each side of the cube that is a block)
  ! This is a tradeoff between 1) having a few large blocks, which reduces
  ! the amount of inter-process communication but also makes the mesh
  ! refinement inefficient (forcing the code to refine large volumes),
  ! and 2) having a lot of small blocks, which makes the AMR very efficient
  ! but increases the amount of communication.
  ! MUST be a power of two. A good value is 16.
  integer, parameter :: ncells_block = 16

  ! Refinement / Coarsening thresholds
  ! These are the maximum relative gradients required to mark a block for
  ! either refinement or coarsening. Note: refineThres must be larger than
  ! twice coarseThres. Good values are 0.3 and 0.03.
  real, parameter :: refineThres = 0.3           !< Refinement threshold
  real, parameter :: coarseThres = 0.03          !< Coarsening threshold
  
  ! ============================================
  ! Boundary conditions
  ! ============================================
  
  ! Currently recognized options are:
  !  BC_REFLECTIVE: Reflective (velocity is reversed)
  !  BC_FREE_FLOW: Free flow (zero-gradients at boundary)
  !  BC_PERIODIC: periodic (boundaries wrap around)
  integer, parameter :: bc_left   = BC_FREE_FLOW     !< BC on LEFT face
  integer, parameter :: bc_right  = BC_FREE_FLOW     !< BC on RIGHT face
  integer, parameter :: bc_front  = BC_FREE_FLOW     !< BC on FRONT face
  integer, parameter :: bc_back   = BC_FREE_FLOW     !< BC on BACK face
  integer, parameter :: bc_bottom = BC_FREE_FLOW     !< BC on BOTTOM face
  integer, parameter :: bc_top    = BC_FREE_FLOW     !< BC on TOP face


  ! ============================================
  ! Data output and logging
  ! ============================================

  ! Data output formats
  !> Output in native Walicxe3D binary format (required for warm starts)?
  logical, parameter :: output_bin = .true.
  !> Output in VisIt-compatible VTK format?
  logical, parameter :: output_vtk = .true.

  !> Output mode: simultaneous or turn-based output?
  !! Currently recognized options:
  !!  OUT_SIMULT: all ranks write output simultaneously
  !!  OUT_TURNS: ranks take turns to write output
  integer, parameter :: output_mode = OUT_SIMULT

  !> Write outputs using physical or code units?
  !! Currently recognized options:
  !!  CODE_UNITS: output in code (scaled) units
  !!  PHYS_UNITS: output in physical (CGS) units
  integer, parameter :: units_type = CODE_UNITS

  ! Data directory and file templates
  ! In these file templates, the sequence XXX will be substituted by
  ! the process number, while the sequence YYYY will be substituted by
  ! the output number. A file extension will be appended automatically
  ! depending on the selected format and should not be given here.
  !> Path to data directory
  character(*), parameter :: datadir = "./DATA/"
  !> Filename template for Blocks data files
  character(*), parameter :: blockstpl = "BlocksXXX.YYYY"
  !> Filename template for Grid data files
  character(*), parameter :: gridtpl   = "Grid.YYYY"
  !> Filename template for State files
  character(*), parameter :: statetpl  = "State.YYYY"
  
  !> Send everything output to stdout to a logfile?
  logical, parameter :: logged = .true.
  !> Directory for logfiles (may be data directory)
  character(*), parameter :: logdir = datadir

  ! ============================================
  ! Solver parameters
  ! ============================================

  !> Numerical Integrator
  !! Currently recognized options:
  !!  SOLVER_LAX: First order Lax-Friedrichs scheme (for testing)
  !!  SOLVER_HLL1: Simplified (first-order) HLL Riemann solver
  !!  SOLVER_HLL: Full (second-order) HLL Riemann solver
  !!  SOLVER_HLLC: HLLC Riemann solver (second-order)
  integer, parameter :: solver_type = SOLVER_HLLC

  !> Slope Limiter (to be used with the HLL/HLLC solvers)
  !! Currently recognized options:
  !!  LIMITER_NONE: use arithmetic average, i.e., no limiter
  !!  LIMITER_MINMOD: Minmod limiter - most diffusive
  !!  LIMITER_VANLEER: van Leer limiter
  !!  LIMITER_ALBADA: van Albada limiter
  !!  LIMITER_UMIST: UMIST limiter 
  !!  LIMITER_WOODWARD: Woodward limiter
  !!  LIMITER_SUPERBEE: Superbee limiter - least diffusive
  integer, parameter :: limiter_type = LIMITER_VANLEER
  
  !> Number of ghost cells (equal to order of solver)
  integer, parameter :: nghost = 2

  !> Number of extra passive scalars
  ! At least one is needed if metallicity-dependent cooling is to be used
  integer, parameter :: npassive = 0

  !> Courant-Friedrichs-Lewis parameter (0 < CFL < 1.0)
  real, parameter :: CFL = 0.4

  !> Artificial viscosity
  real, parameter :: visc_eta = 5.0E-3

  ! ============================================
  ! Radiative Cooling
  ! ============================================

  !> Radiative cooling Type
  ! Currently recognized options:
  !  COOL_NONE: no radiative cooling
  !  COOL_TABLE: tabulated cooling function (temperature only)
  !  COOL_TABLE_METAL: tabulated cooling function (temperature and metallicity)
  integer, parameter :: cooling_type = COOL_NONE

  !> Filename with table of cooling coefficients
  ! Several cooling tables are provided in the tables/ subdirectory.
  character(*), parameter :: cooling_file = "./tables/CHIANTIMazotta.dat"

  !> Maximum *fractional* thermal energy loss allowed in a single cell
  !! per timestep.
  ! This helps prevent negative pressure errors.
  real, parameter :: cooling_limit = 0.5

  ! ============================================
  ! ISM (base flow) Properties
  ! ============================================

  ! General gas parameters

  !> Heat capacity ratio (5/3 for monoatomic ideal gas)
  real, parameter :: gamma = 5.0/3.0
  !> Mean atomic mass of *neutral* gas (in AMUs)
  real, parameter :: mu0 = 1.3
  !> Mean atomic mass of *ionized* gas (in AMUs)
  real, parameter :: mui = 0.61
  !> Gas is considered ionized above this temp (in K)
  real, parameter :: ion_thres = 1.0e4

  ! ISM (base flow) parameters -- all in cgs
  !> ISM mean atomic mass per particle (amu)
  real, parameter :: ism_mu0  = mu0
  !> ISM mass density (g/cm^3)
  real, parameter :: ism_dens = 1.0 * ism_mu0 * AMU
  !> ISM temperature (K)
  real, parameter :: ism_temp = 1.0e3
  !> ISM x-velocity (cm/s)
  real, parameter :: ism_vx = 0.0
  !> ISM y-velocity (cm/s)
  real, parameter :: ism_vy = 0.0
  !> ISM z-velocity (cm/s)
  real, parameter :: ism_vz = 0.0
  !> ISM metallicity (ignored if cooling type is not COOL_TABLE_METAL)
  real, parameter :: ism_metal = 1.0

  ! Magnetic field (only needed when passive magnetic field is enabled)
  ! All values given in gauss
#ifdef PASBP
  real, parameter :: ism_bx = 0.0     !< ISM B_x field (G)
  real, parameter :: ism_by = 0.0     !< ISM B_y field (G)
  real, parameter :: ism_bz = 0.0     !< ISM B_z field (G)
#endif

  ! ============================================
  ! Unit Scalings
  ! ============================================

  ! Unit scaling factors for length, density and velocity.
  ! These are conversion factors between physical (CGS) and code units,
  ! which state the physical unit equivalent of 1 code unit, so that:
  !   physical units = code_units * scaling_factor
  ! All other unit conversions are derived from these three
  real, parameter :: l_sc = 1.0*PC          !< length scale (cm)
  real, parameter :: d_sc = 1.0*mu0*AMU     !< density scale (g cm^-3)
  real, parameter :: v_sc = 1.0e5           !< velocity scale (cm s^-1)
  real, parameter :: pas_sc = d_sc          !< passive scalars scale

  ! ============================================
  ! Additional USER Parameters
  ! ============================================

  ! Add HERE any additional parameters required

  ! ========================================================================== !
  !                                                                            !
  !       Derived paramaters follow - NO NEED TO MODIFY BELOW THIS LINE        !
  !                                                                            !
  ! ========================================================================== !

  !> Number of hydro equations
  integer, parameter :: neqhydro = 5

  !> Number of mhd equations
#ifdef PASBP
  integer, parameter :: neqmhd = 3
#else
  integer, parameter :: neqmhd = 0
#endif

  !> Passive scalar index used for metallicity
  ! Only applicable for COOL_TABLE_METAL cooling
  integer, parameter :: metalpas = neqhydro + neqmhd + min(npassive,1)

  !> Total number of equations to integrate
  !! (hydro variables + mhd variables + passive scalars)
  integer, parameter :: neqtot = neqhydro + neqmhd + npassive

  !> First passive scalar index  
  integer, parameter :: firstpas = neqhydro + neqmhd + 1 

  !> Number of bytes per real
#ifdef DOUBLEP
  integer, parameter :: bytes_per_real = 8
#else
  integer, parameter :: bytes_per_real = 4
#endif

  !> Memory size of one block, in bytes
  integer, parameter :: block_ram_size = &
    (ncells_block+2*nghost)**3 * neqtot * bytes_per_real

  !> Maximum number of blocks per process
  ! Computed automatically from the amount of RAM per process and the 
  ! memory size of each block. The fudge factor in the formula below
  ! leaves a bit of headroom for the code.
  integer, parameter :: nbMaxProc = &
    int((RAM_per_proc*1024*1024 - 6*block_ram_size)*0.95/(block_ram_size*3))
 
  !> Maximum number of blocks across all process
  integer, parameter :: nbMaxGlobal = nbMaxProc*nProcs

  !> Number of cells per block along x
  integer, parameter :: ncells_x = ncells_block 
  !> Number of cells per block along y
  integer, parameter :: ncells_y = ncells_block
  !> Number of cells per block along z
  integer, parameter :: ncells_z = ncells_block

  ! Data array bounds - includes ghost cells
  integer, parameter :: nxmin = 1-nghost          !< Data array bound, x low
  integer, parameter :: nxmax = ncells_x+nghost   !< Data array bound, x high
  integer, parameter :: nymin = 1-nghost          !< Data array bound, y low
  integer, parameter :: nymax = ncells_y+nghost   !< Data array bound, y high
  integer, parameter :: nzmin = 1-nghost          !< Data array bound, z low
  integer, parameter :: nzmax = ncells_z+nghost   !< Data array bound, z high

  !> Floating point precision for reals in MPI messages
#ifdef MPIP
#ifdef DOUBLEP
  integer, parameter :: mpi_real_kind = MPI_DOUBLE_PRECISION
#else
  integer, parameter :: mpi_real_kind = MPI_REAL
#endif
#endif

  !> Heat capacity at constant volume
  real, parameter :: CV = 1.0/(gamma-1.0)

  !> Rank of master process
  integer, parameter :: master = 0
 
  ! Derived unit scalings
  real, parameter :: p_sc = d_sc*v_sc**2
  real, parameter :: e_sc = p_sc
  real, parameter :: t_sc = l_sc/v_sc

end module parameters
