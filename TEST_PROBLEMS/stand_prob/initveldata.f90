module init_vel_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_physbc_module
  use define_bc_module
  use multifab_module
  use fill_3d_module
  use eos_module
  use variables
  use network
  use geometry
  use ml_layout_module
  use ml_restriction_module
  use multifab_fill_ghost_module

  implicit none

  private
  public :: initveldata

contains

  subroutine initveldata(u,s0_init,p0_init,dx,bc,mla)

    type(multifab) , intent(inout) :: u(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: uop(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim),ng
    integer :: i,n,dm,nlevs
    
    dm = mla%dim
    nlevs = mla%nlevel

    ng = u(1)%ng

    do n=1,nlevs
       do i = 1, u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          uop => dataptr(u(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call initveldata_2d(uop(:,:,1,:), lo, hi, ng, dx(n,:), &
                                 s0_init(n,:,:), p0_init(n,:))
          case (3) 
             if (spherical .eq. 1) then
                call bl_error("ERROR: spherical not implemented in initdata")
             else
                call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx(n,:), &
                                   s0_init(n,:,:), p0_init(n,:))
             end if
          end select
       end do
    enddo

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(u(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(u(nlevs),1,1,dm,bc(nlevs))
    else
    
       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(u(n-1),u(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(u(n),u(n-1),ng,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),1,1,dm,fill_crse_input=.false.)
       enddo
       
    end if

  end subroutine initveldata

  subroutine initveldata_2d(u,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only : prob_lo, prob_hi, amp, yhigh, ylow, &
         perturb_type
    use mt19937_module

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind=dp_t)  , intent(in   ) :: s0_init(0:,:)
    real (kind=dp_t)  , intent(in   ) :: p0_init(0:)

    ! Local variables
    integer :: i, j, n
    real (kind=dp_t) :: x, y, Lx
    real (kind=dp_t) :: pert
    real (kind=dp_t) :: rand, ymid, ywidth

    if (perturb_type .eq. 1) then

       u = ZERO

    else if (perturb_type .eq. 2) then

       ymid = HALF*(yhigh + ylow)
       ywidth = yhigh-ymid
       Lx = (prob_hi(1) - prob_lo(1))

       ! initial the velocity
       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j)+HALF) * dx(2)
          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i)+HALF) * dx(1)

             u(i,j,1) = ZERO
             u(i,j,2) = amp*exp(-(y-ymid)**2/ywidth**2)*sin(2.0d0*M_PI*x/Lx)

          enddo
       enddo

    else

       u = ZERO

    end if

  end subroutine initveldata_2d

  subroutine initveldata_3d(u,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only : prob_lo, prob_hi, &
         perturb_type, prob_hi_z, n_cellz, n_celly
    use mt19937_module

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! Local variables
    integer :: i, j, k, iseed
    real (kind=dp_t) :: z, z_0, L_z
    real (kind=dp_t) :: pert
    real (kind=dp_t) :: rand
    real, allocatable :: alpha(:), phi(:)

    if ( perturb_type .eq. 1) then

      u = ZERO

    end if
    
  end subroutine initveldata_3d

end module init_vel_module
