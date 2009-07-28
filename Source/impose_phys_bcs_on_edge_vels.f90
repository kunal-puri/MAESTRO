module impose_phys_bcs_on_edges_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use define_bc_module

  implicit none

  private

  public :: impose_phys_bcs_on_edges

contains

  subroutine impose_phys_bcs_on_edges(u,uedge,the_bc_level)

    use bl_prof_module
    use geometry, only: dm, nlevs

    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(inout) :: uedge(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local variables
    real(kind=dp_t), pointer :: utp(:,:,:,:)
    real(kind=dp_t), pointer :: vtp(:,:,:,:)
    real(kind=dp_t), pointer :: wtp(:,:,:,:)
    integer                  :: lo(dm),hi(dm)
    integer                  :: i,n,ng_ut

    type(bl_prof_timer), save :: bpt

    call build(bpt, "impose_phys_bcs_on_edges")

    ng_ut = uedge(1,1)%ng

    do n=1,nlevs

       do i=1,u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          utp => dataptr(uedge(n,1),i)
          vtp => dataptr(uedge(n,2),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call impose_phys_bcs_2d( &
                              utp(:,:,1,1), vtp(:,:,1,1), ng_ut, lo,hi, &
                              the_bc_level(n)%phys_bc_level_array(i,:,:))
          case (3)
             wtp => dataptr(uedge(n,3), i)
             call impose_phys_bcs_3d( &
                              utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), ng_ut, &
                              lo, hi, &
                              the_bc_level(n)%phys_bc_level_array(i,:,:))
          end select
       end do

    end do

    call destroy(bpt)

  end subroutine impose_phys_bcs_on_edges

  subroutine impose_phys_bcs_2d(uedge,vedge,ng_ut,lo,hi,phys_bc)

    use bc_module

    integer,         intent(in   ) :: lo(:),hi(:),ng_ut
    real(kind=dp_t), intent(inout) :: uedge(lo(1)-ng_ut:,lo(2)-ng_ut:)
    real(kind=dp_t), intent(inout) :: vedge(lo(1)-ng_ut:,lo(2)-ng_ut:)
    integer        , intent(in   ) :: phys_bc(:,:)
    
    ! Local variables
    integer :: is,ie,js,je
    integer :: i,j

    is = lo(1)
    js = lo(2)
    ie = hi(1)
    je = hi(2)
    
    ! impose lo j side bc's on uedge
    if (phys_bc(2,1) .eq. INLET .or. phys_bc(2,1) .eq. SLIP_WALL .or. &
        phys_bc(2,1) .eq. NO_SLIP_WALL) then 
       uedge(is:ie+1,js-1) = ZERO
    else if (phys_bc(2,1) .eq. OUTLET) then
       uedge(is:ie+1,js-1) = uedge(is:ie+1,js)
    end if

    ! impose hi j side bc's on uedge
    if (phys_bc(2,2) .eq. INLET .or. phys_bc(2,2) .eq. SLIP_WALL .or. &
        phys_bc(2,2) .eq. NO_SLIP_WALL) then 
       uedge(is:ie+1,je+1) = ZERO
    else if (phys_bc(2,2) .eq. OUTLET) then
       uedge(is:ie+1,je+1) = uedge(is:ie+1,je)
    end if

    ! impose lo i side bc's on vedge
    if (phys_bc(1,1) .eq. INLET .or. phys_bc(1,1) .eq. SLIP_WALL .or. &
        phys_bc(1,1) .eq. NO_SLIP_WALL) then 
       vedge(is-1,js:je+1) = ZERO
    else if (phys_bc(1,1) .eq. OUTLET) then
       vedge(is-1,js:je+1) = vedge(is,js:je+1)
    end if

    ! impose hi i side bc's on vedge
    if (phys_bc(1,2) .eq. INLET .or. phys_bc(1,2) .eq. SLIP_WALL .or. &
        phys_bc(1,2) .eq. NO_SLIP_WALL) then 
       vedge(ie+1,js:je+1) = ZERO
    else if (phys_bc(2,2) .eq. OUTLET) then
       vedge(ie+1,js:je+1) = vedge(ie,js:je+1)
    end if

  end subroutine impose_phys_bcs_2d

  subroutine impose_phys_bcs_3d(uedge,vedge,wedge,ng_ut,lo,hi,phys_bc)

    use bc_module

    integer,         intent(in   ) :: lo(:),hi(:),ng_ut
    real(kind=dp_t), intent(inout) :: uedge(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(inout) :: vedge(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(inout) :: wedge(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    integer        , intent(in   ) :: phys_bc(:,:)
    
    ! Local variables
    integer :: is,ie,js,je,ks,ke
    integer :: i,j,k

    is = lo(1)
    js = lo(2)
    ks = lo(3)
    ie = hi(1)
    je = hi(2)
    ke = hi(3)
    
    ! impose lo k side bc's on uedge and vedge
    if (phys_bc(3,1) .eq. INLET .or. phys_bc(3,1) .eq. SLIP_WALL .or. &
        phys_bc(3,1) .eq. NO_SLIP_WALL) then 
       uedge(is  :ie+1,js-1:je+1,ks-1) = ZERO
       vedge(is-1:ie+1,js  :je+1,ks-1) = ZERO
    else if (phys_bc(3,1) .eq. OUTLET) then
       uedge(is  :ie+1,js-1:je+1,ks-1) = uedge(is  :ie+1,js-1:je+1,ks)
       vedge(is-1:ie+1,js  :je+1,ke-1) = vedge(is-1:ie+1,js  :je+1,ks)
    end if

    ! impose hi k side bc's on uedge and vedge
    if (phys_bc(3,2) .eq. INLET .or. phys_bc(3,2) .eq. SLIP_WALL .or. &
        phys_bc(3,2) .eq. NO_SLIP_WALL) then 
       uedge(is  :ie+1,js-1:je+1,ke+1) = ZERO
       vedge(is-1:ie+1,js  :je+1,ke+1) = ZERO
    else if (phys_bc(3,2) .eq. OUTLET) then
       uedge(is  :ie+1,js-1:je+1,ke+1) = uedge(is  :ie+1,js-1:je+1,ke)
       vedge(is-1:ie+1,js  :je+1,ke+1) = vedge(is-1:ie+1,js  :je+1,ke)
    end if
    
    ! impose lo j side bc's on uedge and wedge
    if (phys_bc(2,1) .eq. INLET .or. phys_bc(2,1) .eq. SLIP_WALL .or. &
        phys_bc(2,1) .eq. NO_SLIP_WALL) then 
       uedge(is  :ie+1,js-1,ks-1:ke+1) = ZERO
       wedge(is-1:ie+1,js-1,ks  :ke+1) = ZERO
    else if (phys_bc(2,1) .eq. OUTLET) then
       uedge(is  :ie+1,js-1,ks-1:ke+1) = uedge(is  :ie+1,js,ks-1:ke+1)
       wedge(is-1:ie+1,js-1,ks  :ke+1) = wedge(is-1:ie+1,js,ks  :ke+1)
    end if

    ! impose hi j side bc's on uedge and wedge
    if (phys_bc(2,2) .eq. INLET .or. phys_bc(2,2) .eq. SLIP_WALL .or. &
        phys_bc(2,2) .eq. NO_SLIP_WALL) then 
       uedge(is  :ie+1,je+1,ks-1:ke+1) = ZERO
       wedge(is-1:ie+1,je+1,ks  :ke+1) = ZERO
    else if (phys_bc(2,2) .eq. OUTLET) then
       uedge(is  :ie+1,je+1,ks-1:ke+1) = uedge(is  :ie+1,je,ks-1:ke+1)
       wedge(is-1:ie+1,je+1,ks  :ke+1) = wedge(is-1:ie+1,je,ks  :ke+1)
    end if

    ! impose lo i side bc's on vedge and wedge
    if (phys_bc(1,1) .eq. INLET .or. phys_bc(1,1) .eq. SLIP_WALL .or. &
        phys_bc(1,1) .eq. NO_SLIP_WALL) then 
       vedge(is-1,js  :je+1,ks-1:ke+1) = ZERO
       wedge(is-1,js-1:je+1,ks  :ke+1) = ZERO
    else if (phys_bc(1,1) .eq. OUTLET) then
       vedge(is-1,js  :je+1,ks-1:ke+1) = vedge(is,js  :je+1,ks-1:ke+1)
       wedge(is-1,js-1:je+1,ks  :ke+1) = wedge(is,js-1:je+1,ks  :ke+1)
    end if

    ! impose hi i side bc's on vedge and wedge
    if (phys_bc(1,2) .eq. INLET .or. phys_bc(1,2) .eq. SLIP_WALL .or. &
        phys_bc(1,2) .eq. NO_SLIP_WALL) then 
       vedge(ie+1,js  :je+1,ks-1:ke+1) = ZERO
       wedge(ie+1,js-1:je+1,ks  :ke+1) = ZERO
    else if (phys_bc(2,2) .eq. OUTLET) then
       vedge(ie+1,js  :je+1,ks-1:ke+1) = vedge(ie,js  :je+1,ks-1:ke+1)
       wedge(ie+1,js-1:je+1,ks  :ke+1) = wedge(ie,js-1:je+1,ks  :ke+1)
    end if

  end subroutine impose_phys_bcs_3d
  
end module impose_phys_bcs_on_edges_module