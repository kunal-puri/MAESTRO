! advance_premac is the driver routine that orchestrates the creation
! of the edge-centered, 1/2-time advective velocities, without enforcing
! the divergence constraint.

module vstar_module

  use bl_types, only: dp_t
  use multifab_module, only: multifab, multifab_build, multifab_build_edge, &
                             multifab_plus_plus_c, &
                             get_layout, nghost, destroy, max_val, min_val
  use ml_layout_module, only: ml_layout
  use define_bc_module, only: bc_level

  implicit none

  private
  public :: compute_vstar

contains

  subroutine compute_vstar(uold,sold,uface,gpi,normal,w0,w0mac, &
                           w0_force,w0_force_cart_vec,rho0,grav_cell,dx,dt, &
                           the_bc_level,mla)

    use bl_prof_module, only: bl_prof_timer, build, destroy
    use velpred_module, only: velpred
    use addw0_module, only: addw0
    use bl_constants_module, only: ONE
    use variables, only: rho_comp
    use fill_3d_module, only: put_1d_array_on_cart
    use probin_module, only: ppm_trace_forces

    use compute_rhs_module                , only: momentum_rhs
    use interpolate_face_velocities_module, only: interpolate_face_velocities

    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: uface(:,:)
    type(multifab) , intent(in   ) :: gpi(:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: w0_force(:,0:)
    type(multifab) , intent(in   ) :: w0_force_cart_vec(:)
    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    real(kind=dp_t), intent(in   ) :: grav_cell(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla

    ! local
    type(multifab) ::  force(mla%nlevel)
    type(multifab) :: utrans(mla%nlevel,mla%dim)
    type(multifab) ::  ufull(mla%nlevel)
    integer        :: n,comp,dm,nlevs
    logical        :: is_final_update

    type(bl_prof_timer), save :: bpt

    call build(bpt, "compute_vstar")

    dm = mla%dim
    nlevs = mla%nlevel

    do n=1,nlevs
       ! tracing needs more ghost cells
       call multifab_build(force(n),get_layout(uold(n)),dm,uold(n)%ng)
       call multifab_build(ufull(n),get_layout(uold(n)),dm,nghost(uold(n)))
    end do

    ! create fullu = uold + w0
    call put_1d_array_on_cart(w0,ufull,1,.true.,.true.,dx,the_bc_level,mla)
    do n=1,nlevs
       call multifab_plus_plus_c(ufull(n),1,uold(n),1,dm,nghost(uold(n)))
    end do    

    !*************************************************************
    !     Create utrans.
    !*************************************************************

    do n=1,nlevs
       do comp=1,dm
          call multifab_build_edge(utrans(n,comp), mla%la(n),1,1,comp)
       end do
    end do

    call interpolate_face_velocities(uold, utrans, dx, dt, the_bc_level, mla)

    !*************************************************************
    !     Create RHS for the Momentum equation
    !*************************************************************
    call momentum_rhs(force, uold, utrans, sold, rho_comp, dx, the_bc_level, mla)

    !*************************************************************
    !     Create the edge states to be used for the MAC velocity 
    !*************************************************************
    call velpred(uold,uface,utrans,force,dx,dt,the_bc_level,mla)
    !call velpred(uold,ufull,uface,utrans,force,w0,w0mac,dx,dt,the_bc_level,mla)

    do n = 1,nlevs
       call destroy(force(n))
       call destroy(ufull(n))
       do comp=1,dm
          call destroy(utrans(n,comp))
       end do
    end do

    call destroy(bpt)

  end subroutine compute_vstar
  
end module vstar_module
