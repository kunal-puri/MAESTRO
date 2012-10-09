module tag_boxes_module

  use BoxLib
  use f2kcli
  use list_box_module
  use boxarray_module
  use ml_boxarray_module
  use layout_module
  use multifab_module
  use box_util_module
  use bl_IO_module
  use cluster_module
  use ml_layout_module

  implicit none 

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tag_boxes(mf,tagboxes,dx,lev,aux_tag_mf)

    use variables, only: rho_comp, spec_comp
    use geometry, only: nr_fine, nr

    type( multifab)          , intent(in   ) :: mf
    type(lmultifab)          , intent(inout) :: tagboxes
    real(dp_t)               , intent(in   ) :: dx
    integer                  , intent(in   ) :: lev
    type( multifab), optional, intent(in   ) :: aux_tag_mf

    real(kind = dp_t), pointer :: sp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, j, lo(get_dim(mf)), ng_s, dm
    logical           ::      radialtag(0:nr_fine-1)
    logical           :: radialtag_proc(0:nr_fine-1)
    integer, parameter :: npad = 4

    radialtag = .false.
    radialtag_proc = .false.

    dm = get_dim(mf)

    ng_s = mf%ng

    do i = 1, nfabs(mf)
       sp => dataptr(mf, i)
       lo =  lwb(get_box(tagboxes, i))
       select case (dm)
       case (2)
          call radialtag_2d(radialtag_proc,sp(:,:,1,spec_comp),sp(:,:,1,rho_comp),lo,ng_s,lev)
       case  (3)
          call radialtag_3d(radialtag_proc,sp(:,:,:,spec_comp),sp(:,:,:,rho_comp),lo,ng_s,lev)
       end select
    end do

    ! gather radialtag
    call parallel_reduce(radialtag, radialtag_proc, MPI_LOR)

    ! apply some padding                                                                                                  
    do j = 1, npad

       ! pad the start of a tagged region                                                                                 
       do i = 1, nr(lev)-1
          if (radialtag(i) .and. .not. radialtag(i-1)) then
             ! found start of a tagged region                                                                             
             radialtag(i-1) = .true.
          endif
       enddo

       ! pad the end of a tagged region                                                                                   
       do i = nr(lev)-1, 1, -1
          if (radialtag(i) .and. .not. radialtag(i+1)) then
             ! found end of a tagged region                                                                               
             radialtag(i+1) = .true.
          endif
          enddo

    enddo


    do i = 1, nfabs(mf)
       tp => dataptr(tagboxes, i)
       lo =  lwb(get_box(tagboxes, i))
       select case (dm)
       case (2)
          call tag_boxes_2d(tp(:,:,1,1),radialtag,lo)
       case  (3)
          call tag_boxes_3d(tp(:,:,:,1),radialtag,lo)
       end select
    end do

  end subroutine tag_boxes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine radialtag_2d(radialtag,he,rho,lo,ng,lev)

    use probin_module, only: tag_minval, tag_maxval

    integer          , intent(in   ) :: lo(:),ng
    logical          , intent(inout) :: radialtag(0:)
    real(kind = dp_t), intent(in   ) ::  he(lo(1)-ng:,lo(2)-ng:)
    real(kind = dp_t), intent(in   ) :: rho(lo(1)-ng:,lo(2)-ng:)
    integer, optional, intent(in   ) :: lev

    ! local
    integer :: i,j,nx,ny
    real(kind=dp_t) :: Xhe

    nx = size(he,dim=1) - 2*ng
    ny = size(he,dim=2) - 2*ng

    do j = lo(2),lo(2)+ny-1
       do i = lo(1),lo(1)+nx-1
          Xhe = he(i,j)/rho(i,j)
          if (Xhe .gt. tag_minval .and. Xhe .lt. tag_maxval) then
             radialtag(j) = .true.
           end if
        end do
     enddo

  end subroutine radialtag_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine radialtag_3d(radialtag,he,rho,lo,ng,lev)

    use probin_module, only: tag_minval, tag_maxval

    integer          , intent(in   ) :: lo(:),ng
    logical          , intent(inout) :: radialtag(0:)
    real(kind = dp_t), intent(in   ) ::  he(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(kind = dp_t), intent(in   ) :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer, optional, intent(in   ) :: lev

    ! local
    integer :: i,j,k,nx,ny,nz
    real(kind=dp_t) :: Xhe

    nx = size(he,dim=1) - 2*ng
    ny = size(he,dim=2) - 2*ng
    nz = size(he,dim=3) - 2*ng

    do k = lo(3),lo(3)+nz-1
       do j = lo(2),lo(2)+ny-1
          do i = lo(1),lo(1)+nx-1
             Xhe = he(i,j,k)/rho(i,j,k)
             if (Xhe .gt. tag_minval .and. Xhe .lt. tag_maxval) then
                radialtag(k) = .true.
             end if
          end do
       enddo
    end do

  end subroutine radialtag_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tag_boxes_2d(tagbox,radialtag,lo)

    integer          , intent(in   ) :: lo(:)
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):)
    logical          , intent(in   ) :: radialtag(0:)

    integer :: j,ny

    ny = size(tagbox,dim=2)

    tagbox = .false.

    ! tag all boxes with radialtag = .true
    do j = lo(2),lo(2)+ny-1
       tagbox(:,j) = radialtag(j)
    enddo

  end subroutine tag_boxes_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tag_boxes_3d(tagbox,radialtag,lo)

    integer          , intent(in   ) :: lo(:)
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):,lo(3):)
    logical          , intent(in   ) :: radialtag(0:)

    integer :: k,nz

    nz = size(tagbox,dim=3)

    tagbox = .false.

    ! tag all boxes with radialtag = .true.
    do k = lo(3),lo(3)+nz-1
       tagbox(:,:,k) = radialtag(k)
    end do

  end subroutine tag_boxes_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tag_boxes_module