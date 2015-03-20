! a module for storing the geometric information so we don't have to pass it
!
! This module provides the coordinate value for the left edge of a base-state
! zone (r_edge_loc) and the zone center (r_cc_loc).  As always, it is assumed that 
! the base state arrays begin with index 0, not 1.

module geometry

  use bl_types
  use ml_layout_module

  implicit none

  private
  public :: initialize_dx

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initialize_dx(dx,mba,num_levs)

    use probin_module, only: prob_lo, prob_hi

    real(dp_t)       , pointer     :: dx(:,:)
    type(ml_boxarray), intent(in ) :: mba
    integer          , intent(in ) :: num_levs
    
    integer :: n,d,dm

    dm = mba%dim

    write(*,*) 'allocating space for dx', prob_lo(1), prob_hi(1)
    allocate(dx(num_levs,dm))
    
    do d=1,dm
       dx(1,d) = (prob_hi(d)-prob_lo(d)) / real(extent(mba%pd(1),d),kind=dp_t)
    end do
    do n=2,num_levs
       dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
    end do

  end subroutine initialize_dx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module geometry
