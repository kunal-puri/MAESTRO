module plot_variables_module

  use bl_types
  use multifab_module
  use define_bc_module
  use probin_module, only: use_tfromp
  use bl_constants_module, only: HALF

  implicit none

  private

  public :: make_divergence, make_vorticity, make_magvel
  public :: make_processor_number

contains

  !---------------------------------------------------------------------------
  ! make_vorticity
  !---------------------------------------------------------------------------
  subroutine make_vorticity(vort,comp,u,dx,bc)

    use bl_prof_module

    integer        , intent(in   ) :: comp
    type(multifab) , intent(inout) :: vort
    type(multifab) , intent(in   ) :: u
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: bc

    real(kind=dp_t), pointer:: up(:,:,:,:)
    real(kind=dp_t), pointer:: vp(:,:,:,:)
    integer :: lo(get_dim(vort)),hi(get_dim(vort))
    integer :: i,ng_u,ng_v,dm

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_vort")

    dm = get_dim(vort)

    ng_u = nghost(u)
    ng_v = nghost(vort)

    do i = 1, nfabs(u)
       up  => dataptr(u, i)
       vp  => dataptr(vort, i)
       lo  =  lwb(get_box(u, i))
       hi  =  upb(get_box(u, i))
       select case (dm)
       case (1)
          call setval(vort,0.d0,comp,1,all=.true.)
       case (2)
          call makevort_2d(vp(:,:,1,comp),ng_v,up(:,:,1,:),ng_u,lo,hi,dx, &
                           bc%phys_bc_level_array(i,:,:))
       case (3)
          call makevort_3d(vp(:,:,:,comp),ng_v,up(:,:,:,:),ng_u,lo,hi,dx, &
                           bc%phys_bc_level_array(i,:,:))
       end select
    end do

    call destroy(bpt)

  end subroutine make_vorticity

  subroutine makevort_2d(vort,ng_v,u,ng_u,lo,hi,dx,bc)

    use bc_module
    use bl_constants_module

    integer           , intent(in   ) :: lo(:), hi(:), ng_v, ng_u
    real (kind = dp_t), intent(  out) :: vort(lo(1)-ng_v:,lo(2)-ng_v:)  
    real (kind = dp_t), intent(in   ) ::    u(lo(1)-ng_u:,lo(2)-ng_u:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    integer           , intent(in   ) :: bc(:,:)

    !     Local variables
    integer :: i, j
    real (kind = dp_t) :: vx,uy

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          vx = (u(i+1,j,3) - u(i-1,j,3)) / (2.d0*dx(1)) 
          uy = (u(i,j+1,2) - u(i,j-1,2)) / (2.d0*dx(2))
          vort(i,j) = vx - uy
       enddo
    enddo

    if (bc(1,1) .eq. INLET .or. bc(1,1) .eq. SLIP_WALL .or. bc(1,1) .eq. NO_SLIP_WALL) then
       i = lo(1)
       do j = lo(2), hi(2)
          vx = (u(i+1,j,3) + 3.d0*u(i,j,3) - 4.d0*u(i-1,j,3)) / dx(1)
          uy = (u(i,j+1,2) - u(i,j-1,2)) / (2.d0*dx(2))
          vort(i,j) = vx - uy
       end do
    end if

    if (bc(1,2) .eq. INLET .or. bc(1,2) .eq. SLIP_WALL .or. bc(1,2) .eq. NO_SLIP_WALL) then
       i = hi(1)
       do j = lo(2), hi(2)
          vx = -(u(i-1,j,3) + 3.d0*u(i,j,3) - 4.d0*u(i+1,j,3)) / dx(1)
          uy = (u(i,j+1,2) - u(i,j-1,2)) / (2.d0*dx(2))
          vort(i,j) = vx - uy
       end do
    end if

    if (bc(2,1) .eq. INLET .or. bc(2,1) .eq. SLIP_WALL .or. bc(2,1) .eq. NO_SLIP_WALL) then
       j = lo(2)
       do i = lo(1), hi(1)
          vx = (u(i+1,j,3) - u(i-1,j,3)) / (2.d0*dx(1)) 
          uy = (u(i,j+1,2) + 3.d0*u(i,j,2) - 4.d0*u(i,j-1,2)) / dx(2)
          vort(i,j) = vx - uy
       end do
    end if

    if (bc(2,2) .eq. INLET .or. bc(2,2) .eq. SLIP_WALL .or. bc(2,2) .eq. NO_SLIP_WALL) then
       j = hi(2)
       do i = lo(1), hi(1)
          vx =  (u(i+1,j,3) - u(i-1,j,3)) / (2.d0*dx(1)) 
          uy = -(u(i,j-1,2) + 3.d0*u(i,j,2) - 4.d0*u(i,j+1,2)) / dx(2)
          vort(i,j) = vx - uy
       end do
    end if

  end subroutine makevort_2d

  subroutine makevort_3d(vort,ng_v,u,ng_u,lo,hi,dx,bc)

    use bc_module
    use bl_constants_module

    integer           , intent(in   ) :: lo(:), hi(:), ng_v, ng_u
    real (kind = dp_t), intent(  out) :: vort(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:)
    real (kind = dp_t), intent(in   ) ::    u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    integer           , intent(in   ) :: bc(:,:)

    !     Local variables
    integer :: i, j, k
    logical :: fix_lo_x,fix_hi_x,fix_lo_y,fix_hi_y,fix_lo_z,fix_hi_z
    real (kind = dp_t) :: wy,vz,uz,wx,vx,uy

    !$OMP PARALLEL DO PRIVATE(i,j,k,uy,uz,vx,vz,wx,wy)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             uy = uycen(i,j,k)
             uz = uzcen(i,j,k)
             vx = vxcen(i,j,k)
             vz = vzcen(i,j,k)
             wx = wxcen(i,j,k)
             wy = wycen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    fix_lo_x = ( bc(1,1) .eq. INLET .or. bc(1,1) .eq. NO_SLIP_WALL )
    fix_hi_x = ( bc(1,2) .eq. INLET .or. bc(1,2) .eq. NO_SLIP_WALL )

    fix_lo_y = ( bc(2,1) .eq. INLET .or. bc(2,1) .eq. NO_SLIP_WALL )
    fix_hi_y = ( bc(2,2) .eq. INLET .or. bc(2,2) .eq. NO_SLIP_WALL )

    fix_lo_z = ( bc(3,1) .eq. INLET .or. bc(3,1) .eq. NO_SLIP_WALL )
    fix_hi_z = ( bc(3,2) .eq. INLET .or. bc(3,2) .eq. NO_SLIP_WALL )

    !
    !     First do all the faces
    !
    if (fix_lo_x) then
       i = lo(1)
       !$OMP PARALLEL DO PRIVATE(j,k,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(i)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             vx = vxlo(i,j,k)
             wx = wxlo(i,j,k)
             uy = uycen(i,j,k)
             wy = wycen(i,j,k)
             uz = uzcen(i,j,k)
             vz = vzcen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_hi_x) then
       i = hi(1)
       !$OMP PARALLEL DO PRIVATE(j,k,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(i)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             vx = vxhi(i,j,k)
             wx = wxhi(i,j,k)
             uy = uycen(i,j,k)
             wy = wycen(i,j,k)
             uz = uzcen(i,j,k)
             vz = vzcen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_lo_y) then
       j = lo(2)
       !$OMP PARALLEL DO PRIVATE(i,k,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(j)
       do k = lo(3),hi(3)
          do i = lo(1),hi(1)
             vx = vxcen(i,j,k)
             wx = wxcen(i,j,k)
             uy = uylo(i,j,k)
             wy = wylo(i,j,k)
             uz = uzcen(i,j,k)
             vz = vzcen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_hi_y) then
       j = hi(2)
       !$OMP PARALLEL DO PRIVATE(i,k,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(j)
       do k = lo(3),hi(3)
          do i = lo(1),hi(1)
             vx = vxcen(i,j,k)
             wx = wxcen(i,j,k)
             uy = uyhi(i,j,k)
             wy = wyhi(i,j,k)
             uz = uzcen(i,j,k)
             vz = vzcen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_lo_z) then
       k = lo(3)
       !$OMP PARALLEL DO PRIVATE(i,j,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(k)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             vx = vxcen(i,j,k)
             wx = wxcen(i,j,k)
             uy = uycen(i,j,k)
             wy = wycen(i,j,k)
             uz = uzlo(i,j,k)
             vz = vzlo(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_hi_z) then
       k = hi(3)
       !$OMP PARALLEL DO PRIVATE(i,j,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(k)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             vx = vxcen(i,j,k)
             wx = wxcen(i,j,k)
             uy = uycen(i,j,k)
             wy = wycen(i,j,k)
             uz = uzhi(i,j,k)
             vz = vzhi(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if
    !
    !     Next do all the edges
    !
    if (fix_lo_x .and. fix_lo_y) then
       i = lo(1)
       j = lo(2)
       do k = lo(3),hi(3)
          vx = vxlo(i,j,k)
          wx = wxlo(i,j,k)
          uy = uylo(i,j,k)
          wy = wylo(i,j,k)
          uz = uzcen(i,j,k)
          vz = vzcen(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_x .and. fix_lo_y) then
       i = hi(1)
       j = lo(2)
       do k = lo(3),hi(3)
          vx = vxhi(i,j,k)
          wx = wxhi(i,j,k)
          uy = uylo(i,j,k)
          wy = wylo(i,j,k)
          uz = uzcen(i,j,k)
          vz = vzcen(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_lo_x .and. fix_hi_y) then
       i = lo(1)
       j = hi(2)
       do k = lo(3),hi(3)
          vx = vxlo(i,j,k)
          wx = wxlo(i,j,k)
          uy = uyhi(i,j,k)
          wy = wyhi(i,j,k)
          uz = uzcen(i,j,k)
          vz = vzcen(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_x .and. fix_hi_y) then
       i = hi(1)
       j = hi(2)
       do k = lo(3),hi(3)
          vx = vxhi(i,j,k)
          wx = wxhi(i,j,k)
          uy = uyhi(i,j,k)
          wy = wyhi(i,j,k)
          uz = uzcen(i,j,k)
          vz = vzcen(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_lo_x .and. fix_lo_z) then
       i = lo(1)
       k = lo(3)
       do j = lo(2),hi(2)
          vx = vxlo(i,j,k)
          wx = wxlo(i,j,k)
          uy = uycen(i,j,k)
          wy = wycen(i,j,k)
          uz = uzlo(i,j,k)
          vz = vzlo(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_x .and. fix_lo_z) then
       i = hi(1)
       k = lo(3)
       do j = lo(2),hi(2)
          vx = vxhi(i,j,k)
          wx = wxhi(i,j,k)
          uy = uycen(i,j,k)
          wy = wycen(i,j,k)
          uz = uzlo(i,j,k)
          vz = vzlo(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_lo_x .and. fix_hi_z) then
       i = lo(1)
       k = hi(3)
       do j = lo(2),hi(2)
          vx = vxlo(i,j,k)
          wx = wxlo(i,j,k)
          uy = uycen(i,j,k)
          wy = wycen(i,j,k)
          uz = uzhi(i,j,k)
          vz = vzhi(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_x .and. fix_hi_z) then
       i = hi(1)
       k = hi(3)
       do j = lo(2),hi(2)
          vx = vxhi(i,j,k)
          wx = wxhi(i,j,k)
          uy = uycen(i,j,k)
          wy = wycen(i,j,k)
          uz = uzhi(i,j,k)
          vz = vzhi(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_lo_y .and. fix_lo_z) then
       j = lo(2)
       k = lo(3)
       do i = lo(1),hi(1)
          vx = vxcen(i,j,k)
          wx = wxcen(i,j,k)
          uy = uylo(i,j,k)
          wy = wylo(i,j,k)
          uz = uzlo(i,j,k)
          vz = vzlo(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_y .and. fix_lo_z) then
       j = hi(2)
       k = lo(3)
       do i = lo(1),hi(1)
          vx = vxcen(i,j,k)
          wx = wxcen(i,j,k)
          uy = uyhi(i,j,k)
          wy = wyhi(i,j,k)
          uz = uzlo(i,j,k)
          vz = vzlo(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_lo_y .and. fix_hi_z) then
       j = lo(2)
       k = hi(3)
       do i = lo(1),hi(1)
          vx = vxcen(i,j,k)
          wx = wxcen(i,j,k)
          uy = uylo(i,j,k)
          wy = wylo(i,j,k)
          uz = uzhi(i,j,k)
          vz = vzhi(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_y .and. fix_hi_z) then
       j = hi(2)
       k = hi(3)
       do i = lo(1),hi(1)
          vx = vxcen(i,j,k)
          wx = wxcen(i,j,k)
          uy = uyhi(i,j,k)
          wy = wyhi(i,j,k)
          uz = uzhi(i,j,k)
          vz = vzhi(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if
    !
    !     Finally do all the corners
    !
    if (fix_lo_x .and. fix_lo_y .and. fix_lo_z) then
       i = lo(1)
       j = lo(2)
       k = lo(3)
       vx = vxlo(i,j,k)
       wx = wxlo(i,j,k)
       uy = uylo(i,j,k)
       wy = wylo(i,j,k)
       uz = uzlo(i,j,k)
       vz = vzlo(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_hi_x .and. fix_lo_y .and. fix_lo_z) then
       i = hi(1)
       j = lo(2)
       k = lo(3)
       vx = vxhi(i,j,k)
       wx = wxhi(i,j,k)
       uy = uylo(i,j,k)
       wy = wylo(i,j,k)
       uz = uzlo(i,j,k)
       vz = vzlo(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_lo_x .and. fix_hi_y .and. fix_lo_z) then
       i = lo(1)
       j = hi(2)
       k = lo(3)
       vx = vxlo(i,j,k)
       wx = wxlo(i,j,k)
       uy = uyhi(i,j,k)
       wy = wyhi(i,j,k)
       uz = uzlo(i,j,k)
       vz = vzlo(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_hi_x .and. fix_hi_y .and. fix_lo_z) then
       i = hi(1)
       j = hi(2)
       k = lo(3)
       vx = vxhi(i,j,k)
       wx = wxhi(i,j,k)
       uy = uyhi(i,j,k)
       wy = wyhi(i,j,k)
       uz = uzlo(i,j,k)
       vz = vzlo(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_lo_x .and. fix_lo_y .and. fix_hi_z) then
       i = lo(1)
       j = lo(2)
       k = hi(3)
       vx = vxlo(i,j,k)
       wx = wxlo(i,j,k)
       uy = uylo(i,j,k)
       wy = wylo(i,j,k)
       uz = uzhi(i,j,k)
       vz = vzhi(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_hi_x .and. fix_lo_y .and. fix_hi_z) then
       i = hi(1)
       j = lo(2)
       k = hi(3)
       vx = vxhi(i,j,k)
       wx = wxhi(i,j,k)
       uy = uylo(i,j,k)
       wy = wylo(i,j,k)
       uz = uzhi(i,j,k)
       vz = vzhi(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_lo_x .and. fix_hi_y .and. fix_hi_z) then
       i = lo(1)
       j = hi(2)
       k = hi(3)
       vx = vxlo(i,j,k)
       wx = wxlo(i,j,k)
       uy = uyhi(i,j,k)
       wy = wyhi(i,j,k)
       uz = uzhi(i,j,k)
       vz = vzhi(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_hi_x .and. fix_hi_y .and. fix_hi_z) then
       i = hi(1)
       j = hi(2)
       k = hi(3)
       vx = vxhi(i,j,k)
       wx = wxhi(i,j,k)
       uy = uyhi(i,j,k)
       wy = wyhi(i,j,k)
       uz = uzhi(i,j,k)
       vz = vzhi(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

  contains

    function uycen(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i,j+1,k,1)-u(i,j-1,k,1))/dx(2)
    end function uycen

    function uylo(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i,j+1,k,1)+THREE*u(i,j,k,1)-FOUR*u(i,j-1,k,1))/(THREE*dx(2))
    end function uylo

    function uyhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = -(u(i,j-1,k,1)+THREE*u(i,j,k,1)-FOUR*u(i,j+1,k,1))/(THREE*dx(2))
    end function uyhi

    function uzcen(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i,j,k+1,1)-u(i,j,k-1,1))/dx(3)
    end function uzcen

    function uzlo(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i,j,k+1,1)+THREE*u(i,j,k,1)-FOUR*u(i,j,k-1,1))/(THREE*dx(3))
    end function uzlo

    function uzhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r =-(u(i,j,k-1,1)+THREE*u(i,j,k,1)-FOUR*u(i,j,k+1,1))/(THREE*dx(3))
    end function uzhi

    function vxcen(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i+1,j,k,2)-u(i-1,j,k,2))/dx(1)
    end function vxcen

    function vxlo(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i+1,j,k,2)+THREE*u(i,j,k,2)-FOUR*u(i-1,j,k,2))/(THREE*dx(1))
    end function vxlo

    function vxhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r =-(u(i-1,j,k,2)+THREE*u(i,j,k,2)-FOUR*u(i+1,j,k,2))/(THREE*dx(1))
    end function vxhi

    function vzcen(i,j,k) result(r) 
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i,j,k+1,2)-u(i,j,k-1,2))/dx(3)
    end function vzcen

    function vzlo(i,j,k) result(r) 
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i,j,k+1,2)+THREE*u(i,j,k,2)-FOUR*u(i,j,k-1,2))/(THREE*dx(3))
    end function vzlo

    function vzhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r =-(u(i,j,k-1,2)+THREE*u(i,j,k,2)-FOUR*u(i,j,k+1,2))/(THREE*dx(3))
    end function vzhi

    function wxcen(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i+1,j,k,3)-u(i-1,j,k,3))/dx(1)
    end function wxcen

    function wxlo(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i+1,j,k,3)+THREE*u(i,j,k,3)-FOUR*u(i-1,j,k,3))/(THREE*dx(1))
    end function wxlo

    function wxhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r =-(u(i-1,j,k,3)+THREE*u(i,j,k,3)-FOUR*u(i+1,j,k,3))/(THREE*dx(1))
    end function wxhi

    function wycen(i,j,k) result(r) 
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i,j+1,k,3)-u(i,j-1,k,3))/dx(2)
    end function wycen

    function wylo(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i,j+1,k,3)+THREE*u(i,j,k,3)-FOUR*u(i,j-1,k,3))/(THREE*dx(2))
    end function wylo

    function wyhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r =-(u(i,j-1,k,3)+THREE*u(i,j,k,3)-FOUR*u(i,j+1,k,3))/(THREE*dx(2))
    end function wyhi

    function vorfun(uy,uz,vx,vz,wx,wy) result(r)
      real(dp_t) :: uy,uz,vx,vz,wx,wy
      real(dp_t) :: r
      r = sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
    end function vorfun

  end subroutine makevort_3d


  !---------------------------------------------------------------------------
  ! make_magvel
  !---------------------------------------------------------------------------
  subroutine make_magvel(plotdata,comp_magvel,comp_mom,u)

    use bc_module
    use bl_constants_module
    use variables, only : rho_comp

    integer        , intent(in   ) :: comp_magvel, comp_mom
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: u

    real(kind=dp_t), pointer:: pp(:,:,:,:)
    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: up(:,:,:,:)
    real(kind=dp_t), pointer:: wxp(:,:,:,:)
    real(kind=dp_t), pointer:: wyp(:,:,:,:)
    real(kind=dp_t), pointer:: wzp(:,:,:,:)

    integer :: lo(get_dim(plotdata)),hi(get_dim(plotdata)),ng_p,ng_s,ng_u,ng_w
    integer :: i,dm

    dm = get_dim(plotdata)

    ng_u = nghost(u)
    ng_p = nghost(plotdata)

    do i = 1, nfabs(u)
       pp => dataptr(plotdata, i)
       up => dataptr(u, i)
       lo =  lwb(get_box(u, i))
       hi =  upb(get_box(u, i))
       select case (dm)
       case (2)
          call makemagvel_2d(pp(:,:,1,comp_magvel),pp(:,:,1,comp_mom),ng_p, &
                             up(:,:,1,:),ng_u,lo,hi)
       ! case (3)
       !    call makemagvel_3d_cart(pp(:,:,:,comp_magvel),pp(:,:,:,comp_mom),ng_p, &
       !                               sp(:,:,:,rho_comp),ng_s,up(:,:,:,:),ng_u,lo,hi)
       end select
    end do

  end subroutine make_magvel

  ! subroutine makemagvel_1d(magvel,mom,ng_p,rho,ng_s,u,ng_u,lo,hi)

  !   integer           , intent(in   ) :: lo(:), hi(:), ng_p, ng_u, ng_s
  !   real (kind = dp_t), intent(  out) :: magvel(lo(1)-ng_p:)
  !   real (kind = dp_t), intent(  out) ::    mom(lo(1)-ng_p:)
  !   real (kind = dp_t), intent(in   ) ::    rho(lo(1)-ng_s:)
  !   real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u:)

  !   !     Local variables
  !   integer :: i

  !   ! Recall w0 is edge-centered
  !   do i = lo(1), hi(1)
  !      magvel(i) = abs(u(i))
  !      mom(i) = rho(i)*magvel(i)
  !   enddo

  ! end subroutine makemagvel_1d

  subroutine makemagvel_2d(magvel,mom,ng_p,u,ng_u,lo,hi)

    integer           , intent(in   ) :: lo(:), hi(:), ng_p, ng_u
    real (kind = dp_t), intent(  out) :: magvel(lo(1)-ng_p:,lo(2)-ng_p:)
    real (kind = dp_t), intent(  out) ::    mom(lo(1)-ng_p:,lo(2)-ng_p:)
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,:)  

    !     Local variables
    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          magvel(i,j) = sqrt( u(i,j,2)**2 + u(i,j,3)**2 )
          mom(i,j) = magvel(i,j)
       enddo
    enddo

  end subroutine makemagvel_2d

  ! subroutine makemagvel_3d_cart(magvel,mom,ng_p,rho,ng_s,u,ng_u,lo,hi)

  !   integer           , intent(in   ) :: lo(:), hi(:), ng_p, ng_u, ng_s
  !   real (kind = dp_t), intent(  out) :: magvel(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
  !   real (kind = dp_t), intent(  out) ::    mom(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
  !   real (kind = dp_t), intent(in   ) ::    rho(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:) 
  !   real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:) 
  !   real (kind = dp_t), intent(in   ) :: w0(0:)

  !   !     Local variables
  !   integer :: i, j, k
  !   real (kind = dp_t) :: w0_cent

  !   ! Recall w0 is edge-centered

  !   !$OMP PARALLEL DO PRIVATE(i,j,k,w0_cent)
  !   do k = lo(3), hi(3)
  !      do j = lo(2), hi(2)
  !         do i = lo(1), hi(1)
  !            magvel(i,j,k) = sqrt(u(i,j,k,1)**2 + u(i,j,k,2)**2 + (u(i,j,k,3))**2)
  !            mom(i,j,k) = rho(i,j,k)*magvel(i,j,k)
  !         enddo
  !      enddo
  !   enddo
  !   !$OMP END PARALLEL DO

  ! end subroutine makemagvel_3d_cart

  !---------------------------------------------------------------------------
  ! make_processor_number
  !---------------------------------------------------------------------------
  subroutine make_processor_number(plotdata,comp_proc)

    type(multifab), intent(inout) :: plotdata
    integer,        intent(in   ) :: comp_proc
    
    real(kind=dp_t), pointer :: pp(:,:,:,:)
    integer :: i

    do i = 1, nfabs(plotdata)
       pp => dataptr(plotdata, i)
       pp(:,:,:,comp_proc) = parallel_myproc()
    enddo
  end subroutine make_processor_number

  !---------------------------------------------------------------------------
  ! make_divergence
  !---------------------------------------------------------------------------
  subroutine make_divergence(div,comp,u,dx,bc)

    use bl_prof_module

    integer        , intent(in   ) :: comp
    type(multifab) , intent(inout) :: div
    type(multifab) , intent(in   ) :: u
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: bc

    real(kind=dp_t), pointer:: up(:,:,:,:)
    real(kind=dp_t), pointer:: dp(:,:,:,:)
    integer :: lo(get_dim(div)),hi(get_dim(div))
    integer :: i,ng_u,ng_d,dm

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_divergence")

    dm = get_dim(div)

    ng_u = nghost(u)
    ng_d = nghost(div)

    do i = 1, nfabs(u)
       up  => dataptr(u, i)
       dp  => dataptr(div, i)
       lo  =  lwb(get_box(u, i))
       hi  =  upb(get_box(u, i))
       select case (dm)
       case (2)
          call divergence_2d(dp(:,:,1,comp),ng_d,up(:,:,1,:),ng_u,lo,hi,dx, &
                             bc%phys_bc_level_array(i,:,:))
       ! case (3)
       !    call makevort_3d(vp(:,:,:,comp),ng_v,up(:,:,:,:),ng_u,lo,hi,dx, &
       !                     bc%phys_bc_level_array(i,:,:))
       end select
    end do

    call destroy(bpt)

  end subroutine make_divergence

  subroutine divergence_2d(div,ng_d,u,ng_u,lo,hi,dx,bc)

    use bc_module
    use bl_constants_module

    integer           , intent(in   ) :: lo(:), hi(:), ng_d, ng_u
    real (kind = dp_t), intent(  out) :: div(lo(1)-ng_d:,lo(2)-ng_d:)  
    real (kind = dp_t), intent(in   ) :: u(lo(1)-ng_u:,lo(2)-ng_u:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    integer           , intent(in   ) :: bc(:,:)

    !     Local variables
    integer :: i, j
    real (kind = dp_t) :: ux, vy

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          ux = (u(i+1,j,2) - u(i-1,j,2)) / (2.0d0*dx(1))
          vy = (u(i,j+1,3) - u(i,j-1,3)) / (2.0d0*dx(2))
          div(i,j) = ux + vy
       enddo
    enddo

    if (bc(1,1) .eq. INLET .or. bc(1,1) .eq. SLIP_WALL .or. bc(1,1) .eq. NO_SLIP_WALL) then
       i = lo(1)
       do j = lo(2), hi(2)
          ux = (u(i+1,j,2) + 3.d0*u(i,j,2) - 4.d0*u(i-1,j,2)) / dx(1)
          vy = (u(i,j+1,3) - u(i,j-1,3)) / dx(2)
          div(i,j) = ux + vy
       end do
    end if

    if (bc(1,2) .eq. INLET .or. bc(1,2) .eq. SLIP_WALL .or. bc(1,2) .eq. NO_SLIP_WALL) then
       i = hi(1)
       do j = lo(2), hi(2)
          ux = -(u(i-1,j,2) + 3.d0*u(i,j,2) - 4.d0*u(i+1,j,2)) / dx(1)
          vy = (u(i,j+1,3) - u(i,j-1,3)) / (2.d0*dx(2))
          div(i,j) = ux + vy
       end do
    end if

    if (bc(2,1) .eq. INLET .or. bc(2,1) .eq. SLIP_WALL .or. bc(2,1) .eq. NO_SLIP_WALL) then
       j = lo(2)
       do i = lo(1), hi(1)
          ux = (u(i+1,j,2) - u(i-1,j,2)) / (2.d0*dx(1)) 
          vy = (u(i,j+1,3) + 3.d0*u(i,j,3) - 4.d0*u(i,j-1,3)) / dx(2)
          div(i,j) = ux + vy
       end do
    end if

    if (bc(2,2) .eq. INLET .or. bc(2,2) .eq. SLIP_WALL .or. bc(2,2) .eq. NO_SLIP_WALL) then
       j = hi(2)
       do i = lo(1), hi(1)
          ux =  (u(i+1,j,2) - u(i-1,j,2)) / (2.d0*dx(1)) 
          vy = -(u(i,j-1,3) + 3.d0*u(i,j,3) - 4.d0*u(i,j+1,3)) / dx(2)
          div(i,j) = ux + vy
       end do
    end if

  end subroutine divergence_2d


end module plot_variables_module

