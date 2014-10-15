!===============================================================================
!> @file initflow.f90
!> @brief Initial flow conditions module
!> @author Juan C. Toledo
!> @date 2/Feb/2012

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

!> @brief Sets the initial flow conditions
!> @details This wrapper subroutine is called to set the initial flow
!! conditions on the mesh. It will first initialize the grid with a 
!! uniform IC and then call userIC(), where user-defined custom 
!! initial conditions can be specified.
subroutine initflow ()

  use parameters
  use globals
  use userconds
  implicit none

  write(logu,*) ""
  write(logu,*) "============================================"
  write(logu,'(1x,a)') " Setting Initial Conditions  ..."
  write(logu,*) "============================================"
  write(logu,*) ""

  write(logu,'(1x,a)') "> Initializing ambient medium ..."
  write(logu,*) ""

  ! Uniform ISM - base IC
  call uniformIC ()

  ! User-defined IC (defined in user.f90)
  call userInitialCondition (U)

  write(logu,*) ""
  write(logu,'(1x,a)') "> Done setting ICs"
  
end subroutine initflow

! ============================================================================

!> @brief Sets a uniform IC on all blocks
!> @details The initial medium is uniformly filled with gas having
!! the properties defined in the parameters module.
subroutine uniformIC ()

  use parameters
  use globals
  implicit none

  integer :: i, j, k, nb, bID
  real :: primit(neqtot), uvars(neqtot)
  real :: ism_pres

  ! ISM is assumed neutral
  ism_pres = ism_dens/(ism_mu0*AMU)*KB*ism_temp

  write(logu,'(1x,a)') "Ambient medium parameters:"    
  write(logu,'(1x,a,es12.5,a,es12.5,a)') "Dens= ", ism_dens, " g cm^-3   (", ism_dens/d_sc, " code units)"
  write(logu,'(1x,a,es12.5,a,es12.5,a)') "Pres= ", ism_pres, " dyn cm^-2 (", ism_pres/p_sc, ")"
  write(logu,'(1x,a,es12.5,a)') "Temp= ", ism_temp, " K         (no code units)"
  write(logu,'(1x,a,es12.5,a,es12.5,a)') "Velx= ", ism_vx, " cm s^-1   (", ism_vx/p_sc, ")"
  write(logu,'(1x,a,es12.5,a,es12.5,a)') "Vely= ", ism_vy, " cm s^-1   (", ism_vy/p_sc, ")"
  write(logu,'(1x,a,es12.5,a,es12.5,a)') "Velz= ", ism_vz, " cm s^-1   (", ism_vz/p_sc, ")"
#ifdef PASBP
  write(logu,'(1x,a,es12.5,a)') "Bx= ", ism_bx, " G"
  write(logu,'(1x,a,es12.5,a)') "By= ", ism_by, " G"
  write(logu,'(1x,a,es12.5,a)') "Bz= ", ism_bz, " G"
#endif  
  
  ! For all local blocks ...
  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then

      ! Set all physical cells
      do i=1,ncells_x
        do j=1,ncells_y
          do k=1,ncells_z

            ! Calculate scaled primitives for ISM
            primit(1) = ism_dens/d_sc
            primit(2) = ism_vx/v_sc
            primit(3) = ism_vy/v_sc
            primit(4) = ism_vz/v_sc
            primit(5) = ism_pres/p_sc
#ifdef PASBP
            primit(6) = ism_bx
            primit(7) = ism_by
            primit(8) = ism_bz
#endif

            ! Initialize passive scalars, if any
            if (npassive.ge.1) then
              primit(firstpas:neqtot) = 0.0
            end if

            ! Passive scalar for metalicity, for COOL_TABLE_METAL
            if (cooling_type.eq.COOL_TABLE_METAL) then
              primit(metalpas) = ism_metal*primit(1)
            end if

            ! Convert primitives and set flow vars for this cell
            call prim2flow (primit, uvars)
            U(nb,:,i,j,k) = uvars(:)
            
          end do
        end do
      end do

      ! DEBUG
!        write(logu,*) minval(U(nb,1,:,:,:)), maxval(U(nb,1,:,:,:))
      ! DEBUG

    end if
  end do

end subroutine uniformIC 
