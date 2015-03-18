module compute_rhs_module

  ! Compute the force that appears in the velocity (or momentum)
  ! equations.  This is used both when predicting the interface
  ! states and in the final, conservative update.

  ! for the final conservative update of the velocity, we need to
  ! time-center the Coriolis term ( -2 omega x U ), which means we
  ! should use umac.  This is selected by setting is_final_update = T

  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private
  public :: momentum_flux

contains

  subroutine momentum_flux(cflux, dflux, uold, ustar, s, index_rho, dx, the_bc_level, mla)

    ! index_rho refers to the index into s where the density lives.
    ! Usually s will be the full state array, and index_rho would
    ! be rho_comp, but sometimes we pass in only a single-variable
    ! multifab array, so index_rho may be different.

    use bl_prof_module
    use geometry, only: spherical, nr_fine, dr
    use bl_constants_module
    use ml_restrict_fill_module
    use probin_module, only: evolve_base_state
    use fill_3d_module, only : put_1d_array_on_cart
    use variables, only : foextrap_comp

    type(multifab) , intent(inout) :: cflux(:,:), dflux(:,:)
    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(in   ) :: ustar(:,:)
    type(multifab) , intent(in   ) :: s(:)
    integer                        :: index_rho  
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla

    ! Local variables
    real(kind=dp_t), pointer ::  uop(:,:,:,:)
    real(kind=dp_t), pointer ::  usp(:,:,:,:)
    real(kind=dp_t), pointer ::  vsp(:,:,:,:)
    real(kind=dp_t), pointer ::  wsp(:,:,:,:)
    real(kind=dp_t), pointer ::   np(:,:,:,:)

    real(kind=dp_t), pointer ::   rp(:,:,:,:)

    integer                  :: i,r,lo(mla%dim),hi(mla%dim),dm,nlevs,comp
    integer                  :: ng_s,ng_f,n,ng_uo,ng_us

    real(kind=dp_t), pointer ::   fp_x(:,:,:,:)
    real(kind=dp_t), pointer ::   fp_y(:,:,:,:)
   
    type(bl_prof_timer), save :: bpt
    call build(bpt, "compute_flux")

    dm = mla%dim
    nlevs = mla%nlevel

    ! number of ghosts for the scalar
    ng_s  = nghost(s(1))

    ! number of ghosts for the convective flux
    ng_f  = nghost(cflux(1,1))

    do n = 1, nlevs
       do comp = 1,dm
          call setval(cflux(n,comp),ZERO,all=.true.)
          call setval(dflux(n,comp),ZERO,all=.true.)
       end do
    end do

    do n=1,nlevs
       do i=1,nfabs(s(n))

          ! pointer to the scalars
          rp  => dataptr(s(n),i)

          ! The scalars ans uold are cell-centered. lo & hi represents
          ! the domain indices over which we compute the fluxes
          lo = lwb(get_box(s(n),i))
          hi = upb(get_box(s(n),i))

          ! pointer to the edge-centered velocities
          usp => dataptr(ustar(n,1),i)
          ng_us = nghost(ustar(1,1))

          select case (dm)

          case (2)
             vsp => dataptr(ustar(n,2),i)
             
             fp_x => dataptr(cflux(n,1), i)
             fp_y => dataptr(cflux(n,2), i)

             call convective_fluxes_2d(fp_x(:,:,1,:), fp_y(:,:,1,:), ng_f, &
                                       usp(:,:,1,1), vsp(:,:,1,1), ng_us, &
                                       lo, hi, dx, n)
             
          ! case (3)
          !    uop => dataptr(uold(n),i)
          !    vmp => dataptr(umac(n,2),i)
          !    wmp => dataptr(umac(n,3),i)

          !    ng_uo = nghost(uold(1))

          !    call mk_vel_force_3d_cart(fp(:,:,:,:),ng_f,is_final_update, &
          !                              uop(:,:,:,:),ng_uo, &
          !                              ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),ng_um, &
          !                              w0(n,:), &
          !                              gpp(:,:,:,:),ng_gp,rp(:,:,:,index_rho),ng_s, &
          !                              rho0(n,:),grav(n,:),w0_force(n,:),lo,hi,n, &
          !                              do_add_utilde_force)
          end select
       end do
    enddo

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,cflux(:,1),mla%mba%rr,the_bc_level, &
                              icomp=1, &
                              bcomp=1, &
                              nc=dm,   &
                              ng=cflux(1,1)%ng)

    call ml_restrict_and_fill(nlevs,cflux(:,2),mla%mba%rr,the_bc_level, &
                              icomp=1, &
                              bcomp=1, &
                              nc=dm,   &
                              ng=cflux(1,1)%ng)

    call ml_restrict_and_fill(nlevs,dflux(:,1),mla%mba%rr,the_bc_level, &
                              icomp=1, &
                              bcomp=1, &
                              nc=dm,   &
                              ng=dflux(1,1)%ng)

    call ml_restrict_and_fill(nlevs,dflux(:,2),mla%mba%rr,the_bc_level, &
                              icomp=1, &
                              bcomp=1, &
                              nc=dm,   &
                              ng=dflux(1,1)%ng)

    call destroy(bpt)

  end subroutine momentum_flux

  subroutine convective_fluxes_2d(cflux_x, cflux_y, ng_f, ustar, vstar, ng_us, &
                                  lo, hi, dx, n)
    use bl_constants_module

    integer        , intent(in   ) ::  lo(:),hi(:),ng_f,ng_us, n
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(inout) :: cflux_x(lo(1)-ng_f:, lo(2)-ng_f:, :)
    real(kind=dp_t), intent(inout) :: cflux_y(lo(1)-ng_f:, lo(2)-ng_f:, :)
    real(kind=dp_t), intent(in   ) ::    ustar(lo(1)-ng_us:, lo(2)-ng_us:)
    real(kind=dp_t), intent(in   ) ::    vstar(lo(1)-ng_us:, lo(2)-ng_us:)

    !locals
    integer         :: i,j
    real(kind=dp_t) :: dxi, dyi

    ! initialize the flux to 0
    cflux_x = ZERO
    cflux_y = ZERO

    dxi = 1.0d0/dx(n, 1)
    dyi = 1.0d0/dx(n, 2)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute the Convective fluxes for the U-Momentum equation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do j = lo(2),hi(2)
       do i = lo(1), hi(1)+1
          
          cflux_x(i-1, j, 1) = 0.25*( (ustar(i-1,j)+ustar(i,j)) * &
                                      (ustar(i-1,j)+ustar(i,j)) )

          cflux_x(i, j, 1)   = 0.25*( (ustar(i,j)+ustar(i+1,j)) * &
                                      (ustar(i,j)+ustar(i+1,j)) )

          cflux_y(i-1, j, 1) = 0.25*( (vstar(i,j-1)+vstar(i-1,j-1)) * &
                                      (ustar(i,j)+ustar(i,j-1)) )
          
          cflux_y(i, j, 1)   = 0.25*( (vstar(i,j)+vstar(i-1,j)) * &
                                      (ustar(i,j)+ustar(i,j+1)) )

       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute the Convective fluxes for the V-Momentum equation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
          
          cflux_x(i, j-1, 2) = 0.25*( (ustar(i-1,j)+ustar(i-1,j-1)) * &
                                      (vstar(i-1,j-1)+vstar(i,j-1)) )

          cflux_x(i, j, 2)   = 0.25*( (ustar(i,j)+ustar(i,j-1)) * &
                                      (vstar(i,j-1)+vstar(i+1,j-1)) )

          cflux_y(i, j-1, 2) = 0.25*( (vstar(i,j-1)+vstar(i,j)) * &
                                      (vstar(i,j-1)+vstar(i,j)) )

          cflux_y(i, j, 2  ) = 0.25*( (vstar(i,j)+vstar(i,j+1)) * &
                                      (vstar(i,j)+vstar(i,j+1)) )   
       end do
    end do
    

  end subroutine convective_fluxes_2d

  subroutine mk_vel_force_1d(vel_force,ng_f,gpi,ng_gp, &
                             rho,ng_s, &
                             umac,ng_um, &
                             rho0,grav,w0,w0_force,lo,hi,n, &
                             do_add_utilde_force)

    use geometry, only: nr, dr
    use bl_constants_module
    use probin_module, only: base_cutoff_density, buoyancy_cutoff_factor

    integer        , intent(in   ) ::  lo(:),hi(:),ng_f,ng_gp,ng_s,ng_um
    real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f :)
    real(kind=dp_t), intent(in   ) ::     gpi(lo(1)-ng_gp:)
    real(kind=dp_t), intent(in   ) ::     rho(lo(1)-ng_s :)
    real(kind=dp_t), intent(in   ) ::    umac(lo(1)-ng_um:)
    real(kind=dp_t), intent(in   ) :: rho0(0:)
    real(kind=dp_t), intent(in   ) :: grav(0:)
    real(kind=dp_t), intent(in   ) :: w0(0:), w0_force(0:)
    logical        , intent(in   ) :: do_add_utilde_force
    integer        , intent(in   ) :: n
    integer         :: i
    real(kind=dp_t) :: rhopert

    vel_force = ZERO

    do i = lo(1),hi(1)

       rhopert = rho(i) - rho0(i)
       
       ! cutoff the buoyancy term if we are outside of the star
       if (rho(i) .lt. buoyancy_cutoff_factor*base_cutoff_density) then
          rhopert = 0.d0
       end if

       ! note: if use_alt_energy_fix = T, then gphi is already weighted
       ! by beta_0
       vel_force(i) =  rhopert / rho(i) * grav(i) - gpi(i) / rho(i) - w0_force(i)

    end do

  end subroutine mk_vel_force_1d

  ! subroutine mk_vel_force_2d(vel_force,ng_f,gpi,ng_gp, &
  !                            rho,ng_s, &
  !                            vmac, ng_um, &
  !                            rho0,grav,w0,w0_force,lo,hi,n, &
  !                            do_add_utilde_force)

  !   use geometry, only: nr, dr
  !   use bl_constants_module
  !   use probin_module, only: base_cutoff_density, buoyancy_cutoff_factor

  !   integer        , intent(in   ) ::  lo(:),hi(:),ng_f,ng_gp,ng_s,ng_um, n
  !   real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f :,lo(2)-ng_f :,:)
  !   real(kind=dp_t), intent(in   ) ::     gpi(lo(1)-ng_gp:,lo(2)-ng_gp:,:)
  !   real(kind=dp_t), intent(in   ) ::     rho(lo(1)-ng_s :,lo(2)-ng_s :)
  !   real(kind=dp_t), intent(in   ) ::    vmac(lo(1)-ng_um:,lo(2)-ng_um:)
  !   real(kind=dp_t), intent(in   ) :: rho0(0:)
  !   real(kind=dp_t), intent(in   ) :: grav(0:)
  !   real(kind=dp_t), intent(in   ) :: w0(0:),w0_force(0:)
  !   logical        , intent(in   ) :: do_add_utilde_force

  !   integer         :: i,j
  !   real(kind=dp_t) :: rhopert

  !   vel_force = ZERO

  !   do j = lo(2),hi(2)
  !      do i = lo(1),hi(1)

  !         rhopert = rho(i,j) - rho0(j)
          
  !         ! cutoff the buoyancy term if we are outside of the star
  !         if (rho(i,j) .lt. buoyancy_cutoff_factor*base_cutoff_density) then
  !            rhopert = 0.d0
  !         end if

  !         ! note: if use_alt_energy_fix = T, then gphi is already weighted
  !         ! by beta_0
  !         vel_force(i,j,1) = - gpi(i,j,1) / rho(i,j)
  !         vel_force(i,j,2) =  rhopert / rho(i,j) * grav(j) &
  !              - gpi(i,j,2) / rho(i,j) - w0_force(j)
  !      end do
  !   end do

  !   if (do_add_utilde_force) then

  !      do j=lo(2),hi(2)
  !         do i=lo(1),hi(1)

  !            if (j .le. -1) then
  !               ! do not modify force since dw0/dr=0                                                                       
  !            else if (j .ge. nr(n)) then
  !               ! do not modify force since dw0/dr=0                                                                       
  !            else
  !               vel_force(i,j,2) = vel_force(i,j,2) &
  !                    - (vmac(i,j+1)+vmac(i,j))*(w0(j+1)-w0(j)) / (TWO*dr(n))
  !            end if
          
  !         end do
  !      end do

  !   endif


  ! end subroutine mk_vel_force_2d

  ! subroutine mk_vel_force_3d_cart(vel_force,ng_f,is_final_update, &
  !                                 uold,ng_uo, &
  !                                 umac,vmac,wmac,ng_um, &
  !                                 w0, &
  !                                 gpi,ng_gp,rho,ng_s, &
  !                                 rho0,grav,w0_force,lo,hi,n, &
  !                                 do_add_utilde_force)

  !   use geometry,  only: sin_theta, cos_theta, omega, nr, dr
  !   use bl_constants_module
  !   use probin_module, only: base_cutoff_density, buoyancy_cutoff_factor, &
  !                            rotation_radius

  !   integer        , intent(in   ) ::  lo(:),hi(:),ng_f,ng_gp,ng_s, ng_uo, ng_um, n
  !   real(kind=dp_t), intent(inout) :: vel_force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
  !   logical        , intent(in   ) :: is_final_update
  !   real(kind=dp_t), intent(in   ) ::      uold(lo(1)-ng_uo:,lo(2)-ng_uo:,lo(3)-ng_uo:,:)
  !   real(kind=dp_t), intent(in   ) ::      umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
  !   real(kind=dp_t), intent(in   ) ::      vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
  !   real(kind=dp_t), intent(in   ) ::      wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
  !   real(kind=dp_t), intent(in   ) ::   w0(0:)
  !   real(kind=dp_t), intent(in   ) ::     gpi(lo(1)-ng_gp:,lo(2)-ng_gp:,lo(3)-ng_gp:,:)
  !   real(kind=dp_t), intent(in   ) ::       rho(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
  !   real(kind=dp_t), intent(in   ) :: rho0(0:)
  !   real(kind=dp_t), intent(in   ) :: grav(0:)
  !   real(kind=dp_t), intent(in   ) :: w0_force(0:)
  !   logical        , intent(in   ) :: do_add_utilde_force

  !   integer         :: i,j,k
  !   real(kind=dp_t) :: rhopert

  !   real(kind=dp_t) :: coriolis_term(3), centrifugal_term(3)

  !   vel_force = ZERO

  !   ! CURRENTLY for rotation in plane-parallel, we make the (bad) assumption 
  !   ! that all points within the patch have the same centrifugal forcing terms.
  !   !
  !   ! We assume the centrifugal term applies at a constant radius, 
  !   ! rotation_radius, for the patch.  In otherwords, the patch lives on the
  !   ! surface of a sphere of radius rotation_radius.
  !   !
  !   ! Furthermore, we assume the patch lives at longitude = 0.
  !   !
  !   ! Then the orientation of the patch is such that e_z is in the 
  !   ! outward radial direction of the star, e_x is in the co_latitude (polar) 
  !   ! angle direction and e_y is in the global y-direction.
  !   !
  !   ! centrifugal_term = omega x (omega x r) = (omega dot r) * omega
  !   !                                          - omega^2 * r
  !   ! where omega = (-|omega| sin_theta) e_x + (|omega| cos_theta) e_z
  !   !           r = rotation_radius e_z
  !   !
  !   ! See docs/rotation for derivation and figures.
  !   ! 

  !   centrifugal_term(1) = - omega**2 * rotation_radius * sin_theta * sin_theta
  !   centrifugal_term(2) = ZERO
  !   centrifugal_term(3) = omega**2 * rotation_radius * cos_theta * sin_theta &
  !                         - omega**2 * rotation_radius

  !   !$OMP PARALLEL DO PRIVATE(i,j,k,rhopert,coriolis_term)
  !   do k = lo(3),hi(3)
  !      do j = lo(2),hi(2)
  !         do i = lo(1),hi(1)

  !            rhopert = rho(i,j,k) - rho0(k)
             
  !            ! cutoff the buoyancy term if we are outside of the star
  !            if (rho(i,j,k) .lt. buoyancy_cutoff_factor*base_cutoff_density) then
  !               rhopert = 0.d0
  !            end if

  !            ! the coriolis term is:
  !            !    TWO * omega x U
  !            ! where omega is given above and U = (u, v, w) is the velocity

  !            if (is_final_update) then

  !               ! use umac so we are time-centered
  !               coriolis_term(1) = -TWO * omega * &
  !                    HALF*(vmac(i,j,k) + vmac(i,j+1,k)) * cos_theta

  !               coriolis_term(2) =  TWO * omega * &
  !                    (HALF*(wmac(i,j,k)   + w0(k) + &
  !                           wmac(i,j,k+1) + w0(k+1)) * sin_theta + &
  !                     HALF*(umac(i,j,k) + umac(i+1,j,k)) * cos_theta)

  !               coriolis_term(3) = -TWO * omega * &
  !                    HALF*(vmac(i,j,k) + vmac(i,j+1,k)) * sin_theta

  !            else
  !               coriolis_term(1) = -TWO * omega * uold(i,j,k,2) * cos_theta

  !               coriolis_term(2) =  TWO * omega * ((uold(i,j,k,3) + HALF*(w0(k) + w0(k+1))) * sin_theta + &
  !                                                  uold(i,j,k,1) * cos_theta)

  !               coriolis_term(3) = -TWO * omega * uold(i,j,k,2) * sin_theta
  !            endif

  !            ! note: if use_alt_energy_fix = T, then gphi is already
  !            ! weighted by beta_0
  !            vel_force(i,j,k,1) = -coriolis_term(1) - centrifugal_term(1) - &
  !                 gpi(i,j,k,1) / rho(i,j,k) 

  !            vel_force(i,j,k,2) = -coriolis_term(2) - centrifugal_term(2) - &
  !                 gpi(i,j,k,2) / rho(i,j,k) 

  !            vel_force(i,j,k,3) = -coriolis_term(3) - centrifugal_term(3) + &
  !                 ( rhopert * grav(k) - gpi(i,j,k,3) ) / rho(i,j,k) &
  !                 - w0_force(k)

  !         end do
  !      end do
  !   end do
  !   !$OMP END PARALLEL DO


  !   if (do_add_utilde_force) then
  !      !$OMP PARALLEL DO PRIVATE(i,j,k)
  !      do k=lo(3),hi(3)
  !         do j=lo(2),hi(2)
  !            do i=lo(1),hi(1)

  !               if (k .le. -1) then
  !                  ! do not modify force since dw0/dr=0
  !               else if (k .ge. nr(n)) then
  !                  ! do not modify force since dw0/dr=0
  !               else
  !                  vel_force(i,j,k,3) = vel_force(i,j,k,3) &
  !                       - (wmac(i,j,k+1)+wmac(i,j,k))*(w0(k+1)-w0(k)) / (TWO*dr(n))
  !               end if

  !            end do
  !         end do
  !      end do
  !      !$OMP END PARALLEL DO

  !   endif

  ! end subroutine mk_vel_force_3d_cart

end module compute_rhs_module
