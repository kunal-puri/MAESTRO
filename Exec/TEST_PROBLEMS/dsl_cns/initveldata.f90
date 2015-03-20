module init_vel_module

  use multifab_module
  use ml_layout_module
  use bl_constants_module
  use ml_restrict_fill_module

  implicit none

  private
  public :: initveldata

contains

  subroutine initveldata(u,w,dx,the_bc_level,mla,prob_lo)

    type(multifab) , intent(inout) :: u(:),w(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), intent(in   ) :: prob_lo(mla%dim)

    real(kind=dp_t), pointer:: uop(:,:,:,:), wop(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim),ng
    integer :: i,n,dm,nlevs
    
    dm = mla%dim
    nlevs = mla%nlevel
    ng = nghost(u(1))

    do n=1,nlevs

       do i = 1, nfabs(u(n))
          uop => dataptr(u(n),i)
          wop => dataptr(w(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call initveldata_2d(uop(:,:,1,:), wop(:,:,1,:),lo, hi, ng, dx(n,:), &
                                 prob_lo)
          end select
       end do

    enddo

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,u,mla%mba%rr,the_bc_level, &
                              icomp=1, &
                              bcomp=1, &
                              nc=multifab_ncomp(u(1)), &
                              ng=multifab_nghost(u(1)))

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,w,mla%mba%rr,the_bc_level, &
                              icomp=1, &
                              bcomp=1, &
                              nc=multifab_ncomp(w(1)), &
                              ng=multifab_nghost(w(1)))

  end subroutine initveldata

  subroutine initveldata_2d(u,w,lo,hi,ng,dx,prob_lo)
    use probin_module, only: deltaw, deltap, p0, rho0, gamma

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,:),w(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: prob_lo(2)

    ! Local variables
    integer           :: i, j
    real(kind=dp_t)   :: x, y, gamma1i

    gamma1i = 1./(gamma - 1.0d0)
    write(*, *) 'DELTAP, DELTAW', deltap, deltaw, gamma, rho0, p0

    ! initialize the velocity
    do j=lo(2), hi(2)
       y = prob_lo(2) + (dble(j) + 0.5d0) * dx(2)
       do i=lo(1), hi(1)
          x = prob_lo(1) + (dble(i) + 0.5d0) * dx(1)

          ! density
          u(i,j,1) = rho0
          w(i,j,1) = rho0
          
          ! Velocity & Momentum
          if ( y .le. 0.5d0 ) then
             w(i, j, 2) = tanh(deltaw * (y - 0.25d0))
          else
             w(i, j, 2) = tanh(deltaw * (0.75d0 - y))
          end if

          w(i, j, 3) = deltap * sin(2*M_PI*(x + 0.25d0))

          ! Momentum
          u(i, j, 2) = rho0*w(i, j, 2)
          u(i, j, 3) = rho0*w(i, j, 3)

          ! Pressure and Energy
          w(i,j,4) = p0
          u(i,j,4) = rho0*(0.5*(u(i,j,2)**2 + u(i,j,3)**2) + &
                           w(i,j,4)*gamma1i/rho0)

       end do
    end do
    
  end subroutine initveldata_2d

  ! subroutine initveldata_3d(u,lo,hi,ng,dx,s0_init,p0_init)

  !   integer           , intent(in   ) :: lo(:), hi(:), ng
  !   real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
  !   real (kind = dp_t), intent(in   ) :: dx(:)
  !   real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
  !   real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

  !   ! Local variables

  !   ! initial the velocity
  !   u = ZERO
    
  ! end subroutine initveldata_3d

end module init_vel_module
