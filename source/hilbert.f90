!===============================================================================
!> @file hilbert.f90
!> @brief Computes Hilbert keys
!> @author Juan C. Toledo
!> @date 30/Jan/2012

!! CHECK correct copyright notice for this subroutine

!! @n Note: this subroutine was adapted from a subroutine of the RAMSES code
!! RAMSES is developed by Romain Teyssier (the copyright holder is the
!! Commissariat a l'energie atomique et aux energies alternatives, CEA) and
!! is published under the CeCILL Free Software license, which can be 
!! in http://www.cecill.info/

!===============================================================================

!> @brief Returns the key along a 3D Hilbert Curve corresponding to a point
!> @details The list is obtained by calculating the index of each block
!! along a Hilbert curve, and then sorting the block list using the keys.
!> @param x The x-coordinate of the point (starts at 0)
!> @param y The y-coordinate of the point (starts at 0)
!> @param z The z-coordinate of the point (starts at 0)
!> @param order The order of the Hilbert curve
!> @param hk The key along the Hilbert curve (starts at 0)
subroutine HKey3 (x, y, z, order, hk)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: x
  integer, intent(in) :: y
  integer, intent(in) :: z
  integer, intent(in) :: order
  integer, intent(out) :: hk

  logical, dimension(0:3*order-1) :: i_bit_mask
  logical, dimension(0:1*order-1) :: x_bit_mask, y_bit_mask, z_bit_mask
  integer, dimension(0:7,0:1,0:11) :: state_diagram
  integer :: i, cstate, nstate, b0, b1, b2, sdigit, hdigit

  if (order > bit_size(order)) then
     write(logu,*) 'Maximum Hilbert curve order=', bit_size(order)
     call clean_abort (ERROR_HILBERT_ORDER)
  end if

  state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
                            &   0, 1, 3, 2, 7, 6, 4, 5,&
                            &   2, 6, 0, 7, 8, 8, 0, 7,&
                            &   0, 7, 1, 6, 3, 4, 2, 5,&
                            &   0, 9,10, 9, 1, 1,11,11,&
                            &   0, 3, 7, 4, 1, 2, 6, 5,&
                            &   6, 0, 6,11, 9, 0, 9, 8,&
                            &   2, 3, 1, 0, 5, 4, 6, 7,&
                            &  11,11, 0, 7, 5, 9, 0, 7,&
                            &   4, 3, 5, 2, 7, 0, 6, 1,&
                            &   4, 4, 8, 8, 0, 6,10, 6,&
                            &   6, 5, 1, 2, 7, 4, 0, 3,&
                            &   5, 7, 5, 3, 1, 1,11,11,&
                            &   4, 7, 3, 0, 5, 6, 2, 1,&
                            &   6, 1, 6,10, 9, 4, 9,10,&
                            &   6, 7, 5, 4, 1, 0, 2, 3,&
                            &  10, 3, 1, 1,10, 3, 5, 9,&
                            &   2, 5, 3, 4, 1, 6, 0, 7,&
                            &   4, 4, 8, 8, 2, 7, 2, 3,&
                            &   2, 1, 5, 6, 3, 0, 4, 7,&
                            &   7, 2,11, 2, 7, 5, 8, 5,&
                            &   4, 5, 7, 6, 3, 2, 0, 1,&
                            &  10, 3, 2, 6,10, 3, 4, 4,&
                            &   6, 1, 7, 0, 5, 2, 4, 3 /), &
                            & (/8 ,2, 12 /) )

  ! convert to binary
  do i=0,order-1
    x_bit_mask(i) = btest(x,i)
    y_bit_mask(i) = btest(y,i)
    z_bit_mask(i) = btest(z,i)
  enddo

  ! interleave bits
  do i=0,order-1
    i_bit_mask(3*i+2) = x_bit_mask(i)
    i_bit_mask(3*i+1) = y_bit_mask(i)
    i_bit_mask(3*i  ) = z_bit_mask(i)
  end do

  ! build Hilbert ordering using state diagram
  cstate=0
  do i=order-1,0,-1
    b2 = 0
    if (i_bit_mask(3*i+2)) b2 = 1
    b1 = 0
    if (i_bit_mask(3*i+1)) b1 = 1
    b0 = 0
    if (i_bit_mask(3*i)) b0 = 1
    sdigit = b2*4 + b1*2 + b0
    nstate = state_diagram (sdigit, 0, cstate)
    hdigit = state_diagram (sdigit, 1, cstate)
    i_bit_mask(3*i+2) = btest(hdigit, 2)
    i_bit_mask(3*i+1) = btest(hdigit, 1)
    i_bit_mask(3*i  ) = btest(hdigit, 0)
    cstate = nstate
  enddo

  ! save Hilbert key as an integer
  ! THIS MIGHT NOT BE ENOUGH
  hk = 0
  do i=0,3*order-1
    b0 = 0
    if (i_bit_mask(i)) b0 = 1
    hk = hk + b0*2**i
  end do

end subroutine HKey3
