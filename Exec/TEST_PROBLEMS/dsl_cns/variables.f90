!
! A module to provide integer indices into the various storage arrays
! for accessing the different variables by name.
!
module variables

  use bl_types

  implicit none

  integer, save :: rho_comp, rhoh_comp
  integer, save :: press_comp, div_comp
  integer, save :: foextrap_comp, hoextrap_comp
  integer, save :: icomp_vel, icomp_rho, icomp_rhoh
  integer, save :: icomp_magvel, icomp_mom, icomp_vort
  integer, save :: icomp_machno,icomp_cs,icomp_press
  integer, save :: icomp_proc

  ! divergence
  integer, save :: icomp_div

  ! the total number of plot components
  integer, save :: n_plot_comps = 0
  integer, save :: n_comp = 0

  integer, save :: nscal
  real(kind=dp_t), save :: rel_eps

contains

  function get_next_plot_index(num) result (next)

    ! return the next starting index for a plotfile quantity,
    ! and increment the counter of plotfile quantities by num
    integer :: num, next

    next = n_plot_comps + 1
    n_plot_comps = n_plot_comps + num

    return
  end function get_next_plot_index

  subroutine init_variables()

    use probin_module, only: dm_in

    rho_comp    = 1

    if (dm_in .eq. 3) then
       n_comp      = 5
       rhoh_comp   = 5
       press_comp  = 5
    else
       n_comp     = 4
       rhoh_comp  = 4
       press_comp = 4
    end if

    foextrap_comp = press_comp + 1
    hoextrap_comp = foextrap_comp + 1

    ! divergence
    div_comp      = hoextrap_comp + 1

  end subroutine init_variables

  subroutine init_plot_variables()

    use probin_module, only: dm_in, plot_processors

    icomp_rho      = get_next_plot_index(1)
    icomp_vel      = get_next_plot_index(dm_in)
    icomp_rhoh     = get_next_plot_index(1)
    icomp_magvel   = get_next_plot_index(1)
    icomp_mom      = get_next_plot_index(1)
    icomp_vort     = get_next_plot_index(1)
    icomp_press    = get_next_plot_index(1)

    if (plot_processors) then
       icomp_proc = get_next_plot_index(1)
    endif

    icomp_div         = get_next_plot_index(1)
    
  end subroutine init_plot_variables

end module variables
