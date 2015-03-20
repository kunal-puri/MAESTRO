module init_scalar_module

  use multifab_module
  use ml_layout_module
  use bl_constants_module
  use ml_restrict_fill_module
  use variables

  implicit none

  private
  public :: initscalardata, initscalardata_on_level

contains

  subroutine initscalardata(s,dx,the_bc_level,mla,prob_lo)

    use geometry, only: spherical

    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla
    real(kind=dp_t), intent(in   ) :: prob_lo(mla%dim)

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim),ng
    integer :: i,n,dm,nlevs
    
    dm = mla%dim
    nlevs = mla%nlevel

    ng = nghost(s(1))

    do n=1,nlevs
       do i = 1, nfabs(s(n))
          sop => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))

          select case (dm)
          case (2)
             call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx(n,:))
          ! case (3)
          !    call initscalardata_3d(sop(:,:,:,:), lo, hi, ng, dx(n,:), s0_init(n,:,:), &
          !         p0_init(n,:))
          end select
       end do
    enddo

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,s,mla%mba%rr,the_bc_level, &
                              icomp=rho_comp, &
                              bcomp=dm+rho_comp, &
                              nc=nscal, &
                              ng=s(1)%ng)

  end subroutine initscalardata

  subroutine initscalardata_on_level(n,s,dx,the_bc_level)

    use geometry, only: spherical

    integer        , intent(in   ) :: n
    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: the_bc_level

    ! local
    integer                  :: ng,i,dm
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    real(kind=dp_t), pointer :: sop(:,:,:,:)

    dm = get_dim(s)
    ng = nghost(s)

    do i = 1, nfabs(s)
       sop => dataptr(s,i)
       lo =  lwb(get_box(s,i))
       hi =  upb(get_box(s,i))
       select case (dm)
       case (2)
          call initscalardata_2d(sop(:,:,1,:),lo,hi,ng,dx)
       ! case (3)
       !    call initscalardata_3d(sop(:,:,:,:),lo,hi,ng,dx,s0_init,p0_init)
       end select
    end do

  end subroutine initscalardata_on_level

  subroutine initscalardata_2d(s,lo,hi,ng,dx)

    use probin_module, only: prob_lo, perturb_model
    use init_perturb_module

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer         :: i,j
    real(kind=dp_t) :: x,y
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

    ! initial the domain with the base state
    s = ZERO

    ! initialize the scalars
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          s(i,j,rho_comp)  = 1.0
          s(i,j,rhoh_comp) = 1.0
       enddo
    enddo
    
  end subroutine initscalardata_2d

  ! subroutine initscalardata_3d(s,lo,hi,ng,dx,s0_init,p0_init)

  !   use probin_module, only: prob_lo, perturb_model
  !   use init_perturb_module
    
  !   integer           , intent(in   ) :: lo(:),hi(:),ng
  !   real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
  !   real (kind = dp_t), intent(in   ) :: dx(:)
  !   real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
  !   real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

  !   !     Local variables
  !   integer         :: i,j,k
  !   real(kind=dp_t) :: x,y,z
  !   real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
  !   real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

  !   ! initial the domain with the base state
  !   s = ZERO
    
  !   ! initialize the scalars
  !   do k = lo(3), hi(3)
  !      do j = lo(2), hi(2)
  !         do i = lo(1), hi(1)
  !            s(i,j,k,rho_comp)  = s0_init(k,rho_comp)
  !            s(i,j,k,rhoh_comp) = s0_init(k,rhoh_comp)
  !            s(i,j,k,temp_comp) = s0_init(k,temp_comp)
  !            s(i,j,k,spec_comp:spec_comp+nspec-1) = &
  !                 s0_init(k,spec_comp:spec_comp+nspec-1)
  !            s(i,j,k,trac_comp:trac_comp+ntrac-1) = &
  !                 s0_init(k,trac_comp:trac_comp+ntrac-1)
  !         enddo
  !      enddo
  !   enddo

  !   if (perturb_model) then

  !      ! add an optional perturbation
  !      do k = lo(3), hi(3)
  !         z = prob_lo(3) + (dble(k)+HALF) * dx(3)

  !         do j = lo(2), hi(2)
  !            y = prob_lo(2) + (dble(j)+HALF) * dx(2)

  !            do i = lo(1), hi(1)
  !               x = prob_lo(1) + (dble(i)+HALF) * dx(1)

  !               call perturb_3d(x, y, z, p0_init(k), s0_init(k,:), &
  !                               dens_pert, rhoh_pert, rhoX_pert, temp_pert, trac_pert)

  !               s(i,j,k,rho_comp) = dens_pert
  !               s(i,j,k,rhoh_comp) = rhoh_pert
  !               s(i,j,k,temp_comp) = temp_pert
  !               s(i,j,k,spec_comp:spec_comp+nspec-1) = rhoX_pert(:)
  !               s(i,j,k,trac_comp:trac_comp+ntrac-1) = trac_pert(:)
  !            enddo
  !         enddo
  !      enddo
  !   endif

  ! end subroutine initscalardata_3d

end module init_scalar_module
