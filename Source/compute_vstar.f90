! advance_premac is the driver routine that orchestrates the creation
! of the edge-centered, 1/2-time advective velocities, without enforcing
! the divergence constraint.

module vstar_module

  use bl_types, only: dp_t
  use multifab_module, only: multifab, multifab_build, multifab_build_edge, &
                             multifab_plus_plus_c, setval, &
                             get_layout, nghost, destroy, max_val, min_val
  use ml_layout_module, only: ml_layout
  use define_bc_module, only: bc_level

  implicit none

  private
  public :: compute_vstar

contains

  subroutine compute_vstar(uold,sold,ustar,dx,dt,the_bc_level,mla)

    use bl_prof_module, only: bl_prof_timer, build, destroy
    use velpred_module, only: velpred
    use addw0_module, only: addw0
    use bl_constants_module, only: ZERO, ONE
    use variables, only: rho_comp
    use fill_3d_module, only: put_1d_array_on_cart
    use probin_module, only: ppm_trace_forces

    use compute_rhs_module                , only: momentum_force
    use interpolate_face_velocities_module, only: interpolate_face_velocities

    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: ustar(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla

    ! local
    type(multifab) ::  cforce(mla%nlevel,mla%dim), dforce(mla%nlevel,mla%dim)
    type(multifab) :: utrans(mla%nlevel,mla%dim)
    type(multifab) ::  ufull(mla%nlevel)
    integer        :: n,comp,dm,nlevs
    logical        :: is_final_update

    type(bl_prof_timer), save :: bpt

    call build(bpt, "compute_vstar")

    dm = mla%dim
    nlevs = mla%nlevel

    !*************************************************************
    !     Interpolate the cell-centered velocities to the edges
    !*************************************************************
    call interpolate_face_velocities(uold, ustar, dx, dt, the_bc_level, mla)

    !*************************************************************
    !     Create the forcees for the Momentum equation
    !*************************************************************
    do n=1,nlevs
       do comp = 1,dm
          call multifab_build_edge(cforce(n,comp),get_layout(uold(n)),dm,1,comp)
          call multifab_build_edge(dforce(n,comp),get_layout(uold(n)),dm,1,comp)
       end do
    end do

    call momentum_force(cforce, dforce, uold, ustar, sold, rho_comp, dx, the_bc_level, mla)

    !*************************************************************
    !     Predict the (non-divergence free) edge velocities
    !*************************************************************
    call velpred(uold,ustar,cforce,dforce,dx,dt,the_bc_level,mla)

    do n = 1,nlevs
       do comp = 1,dm
          call destroy(cforce(n,comp))
          call destroy(dforce(n,comp))
       end do
    end do

    call destroy(bpt)

  end subroutine compute_vstar
end module vstar_module
