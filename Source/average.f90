! Given a multifab of data (phi), average down to a base state quantity, phibar.
! If we are in plane-parallel, the averaging is at constant height.  
! If we are spherical, then the averaging is done at constant radius.  

module average_module
  
  use bl_types        , only: dp_t
  use box_module      , only: lwb, upb, box
  use multifab_module , only: multifab, nfabs, dataptr, &
                              get_layout, get_box, nghost, nboxes, get_dim
  use layout_module   , only: get_pd
  use ml_layout_module, only: ml_layout
  use parallel        , only: parallel_reduce, MPI_SUM
  use bl_error_module , only: bl_error

  implicit none

  private

  public :: average, average_irreg, average_one_level

contains

  subroutine average(mla,phi,phibar,dx,incomp)

    use geometry, only: nr_fine, nr_irreg, r_start_coord, r_end_coord, spherical, &
                        numdisjointchunks, dr
    use bl_prof_module, only: bl_prof_timer, build, destroy
    use bl_constants_module, only: ZERO, HALF
    use restrict_base_module, only: restrict_base, fill_ghost_base
    use probin_module, only: max_levs, drdxfac

    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: incomp
    type(multifab) , intent(in   ) :: phi(:)
    real(kind=dp_t), intent(inout) :: phibar(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    ! local
    real(kind=dp_t), pointer     :: pp(:,:,:,:)
    logical, pointer             :: mp(:,:,:,:)

    type(box)                    :: domain
    integer                      :: domlo(mla%dim),domhi(mla%dim),lo(mla%dim),hi(mla%dim)
    integer                      :: max_rcoord(mla%nlevel),rcoord(mla%nlevel)
    integer                      :: stencil_coord(mla%nlevel)
    integer                      :: dm,nlevs,i,j,r,n,ng,min_all,min_lev
    real(kind=dp_t)              :: radius

    integer, allocatable ::  ncell_proc(:,:)
    integer, allocatable ::       ncell(:,:)

    integer, allocatable :: which_lev(:)

    real(kind=dp_t), allocatable :: phisum_proc(:,:)
    real(kind=dp_t), allocatable ::      phisum(:,:)

    real(kind=dp_t), allocatable :: radii(:,:)

    type(bl_prof_timer), save :: bpt

    logical :: limit

    call build(bpt, "average")

    dm = mla%dim
    nlevs = mla%nlevel

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! NOTE: The indices for ncell_proc, etc., are switched in this subroutine
    !       to have the number of levels second because due to memory ordering, 
    !       it will run much faster for larger problems
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (spherical .eq. 1) then

       allocate(ncell_proc ( 0:nr_irreg,nlevs))
       allocate(ncell      (-1:nr_irreg,nlevs))

       allocate(which_lev(0:nr_fine-1))

       allocate(phisum_proc(0:nr_irreg,nlevs))
       allocate(phisum     (-1:nr_irreg,nlevs))
       allocate(radii      (-1:nr_irreg+1,nlevs))
       !
       ! radii contains every possible distance that a cell-center at the finest
       ! level can map into
       !
       !$OMP PARALLEL PRIVATE(r,n)
       do n=1,nlevs
          !$OMP DO
          do r=0,nr_irreg
             radii(r,n) = sqrt(0.75d0+2.d0*r)*dx(n,1)
          end do
          !$OMP END DO NOWAIT
       end do
       !$OMP END PARALLEL

       radii(nr_irreg+1,:) = 1.d99

    else

       allocate(ncell_proc (0:nr_fine-1,nlevs))
       allocate(ncell      (0:nr_fine-1,nlevs))
       allocate(phisum_proc(0:nr_fine-1,nlevs))
       allocate(phisum     (0:nr_fine-1,nlevs))

    end if

    ng = nghost(phi(1))

    phibar       = ZERO
    ncell        = 0
    ncell_proc   = 0
    phisum       = ZERO       
    phisum_proc  = ZERO

    if (spherical .eq. 0) then
       
       ! The plane-parallel case is straightforward.
       ! Simply average all the cells values at a particular height.

       do n=1,nlevs

          domain = get_pd(get_layout(phi(n)))
          domlo  = lwb(domain)
          domhi  = upb(domain)

          if (dm .eq. 1) then
             ncell(:,n) = 1
          else if (dm .eq. 2) then
             ncell(:,n) = domhi(1)-domlo(1)+1
          else if (dm .eq. 3) then
             ncell(:,n) = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
          end if

          do i=1, nfabs(phi(n))
             pp => dataptr(phi(n), i)
             lo =  lwb(get_box(phi(n), i))
             hi =  upb(get_box(phi(n), i))
             select case (dm)
             case (1)
                call sum_phi_1d(pp(:,1,1,:),phisum_proc(:,n),lo,hi,ng,incomp)
             case (2)
                call sum_phi_2d(pp(:,:,1,:),phisum_proc(:,n),lo,hi,ng,incomp)
             case (3)
                call sum_phi_3d(pp(:,:,:,:),phisum_proc(:,n),lo,hi,ng,incomp)
             end select
          end do

          call parallel_reduce(phisum(:,n), phisum_proc(:,n), MPI_SUM)

          ! compute phibar by normalizing phisum
          do i=1,numdisjointchunks(n)
             do r=r_start_coord(n,i),r_end_coord(n,i)
                phibar(n,r) = phisum(r,n) / dble(ncell(r,n))
             end do
          end do

       end do

       call restrict_base(phibar,.true.)
       call fill_ghost_base(phibar,.true.)

    else if(spherical .eq. 1) then

       ! For spherical, we construct a 1D array at each level, phisum(:,:), that has space
       ! allocated for every possible radius that a cell-center at each
       ! level can map into.  The radial locations have been precomputed and stored
       ! in radii(:,:).
       do n=nlevs,1,-1

          do i=1, nfabs(phi(n))
             pp => dataptr(phi(n), i)
             lo =  lwb(get_box(phi(n), i))
             hi =  upb(get_box(phi(n), i))

             if (n .eq. nlevs) then
                call sum_phi_3d_sphr(radii(0:,n),nr_irreg,pp(:,:,:,:),phisum_proc(:,n), &
                                     lo,hi,ng,dx(n,:),ncell_proc(:,n),incomp)
             else
                ! we include the mask so we don't double count; i.e., we only consider
                ! cells that are not covered by finer cells when constructing the sum
                mp => dataptr(mla%mask(n), i)
                call sum_phi_3d_sphr(radii(0:,n),nr_irreg,pp(:,:,:,:),phisum_proc(:,n), &
                                     lo,hi,ng,dx(n,:),ncell_proc(:,n),incomp, &
                                     mp(:,:,:,1))
             end if
          end do

          call parallel_reduce(ncell (0:,n), ncell_proc (:,n), MPI_SUM)
          call parallel_reduce(phisum(0:,n), phisum_proc(:,n), MPI_SUM)
       end do

       ! normalize phisum so it actually stores the average at a radius
       do n=1,nlevs
          do r=0,nr_irreg
             if (ncell(r,n) .ne. 0.d0) then
                phisum(r,n) = phisum(r,n) / dble(ncell(r,n))
             end if
          end do
       end do

       ! compute center point for the finest level
       phisum(-1,nlevs) = (11.d0/8.d0)*phisum(0,nlevs) - (3.d0/8.d0)*phisum(1,nlevs)
       radii (-1,nlevs) = 0.d0
       ncell (-1,nlevs) = 1

       ! choose which level to interpolate from
       do n=1,nlevs
          rcoord(n) = 0
       end do

       !$OMP PARALLEL DO PRIVATE(r,radius,n,j,min_all,min_lev) FIRSTPRIVATE(rcoord)
       do r=0,nr_fine-1

         radius = (dble(r)+HALF)*dr(1)

         ! for each level, find the closest coordinate
         do n=1,nlevs
            do j=rcoord(n),nr_irreg
               if (abs(radius-radii(j,n)) .lt. abs(radius-radii(j+1,n))) then
                  rcoord(n) = j
                  exit
               end if
            end do
         end do

         ! make sure closest coordinate is in bounds
         do n=1,nlevs-1
            rcoord(n) = max(rcoord(n),1)
         end do
         do n=1,nlevs
            rcoord(n) = min(rcoord(n),nr_irreg-1)
         end do
         
         ! choose the level with the largest min over the ncell interpolation points
         which_lev(r)=1
         min_all = min(ncell(rcoord(1)-1,1), &
                       ncell(rcoord(1)  ,1), &
                       ncell(rcoord(1)+1,1))

         do n=2,nlevs
            min_lev = min(ncell(rcoord(n)-1,n), &
                          ncell(rcoord(n)  ,n), &
                          ncell(rcoord(n)+1,n))
            
            if (min_lev .gt. min_all) then
               min_all = min_lev
               which_lev(r) = n
            end if
         end do

         ! if the min hit count at all levels is zero, we expand the search
         ! to find the closest instance of where the hitcount becomes nonzero
         j = 1
         do while (min_all .eq. 0)
            j = j+1
            do n=1,nlevs
               min_lev = max(ncell(max(1,rcoord(n)-j),n), &
                             ncell(min(rcoord(n)+j,nr_irreg-1),n))
               if (min_lev .ne. 0) then
                  which_lev(r) = n
                  min_all = min_lev
                  exit
               end if
            end do
         end do

      end do
      !$OMP END PARALLEL DO

      ! squish the list at each level down to exclude points with no contribution
      do n=1,nlevs
         j=0
         do r=0,nr_irreg
            do while(ncell(j,n) .eq. 0)
               j = j+1
               if (j .gt. nr_irreg) then
                  exit
               end if
            end do
            if (j .gt. nr_irreg) then
               phisum(r:nr_irreg,n)   = 1.d99
               radii (r:nr_irreg+1,n) = 1.d99
               max_rcoord(n) = r-1
               exit
            end if
            phisum(r,n) = phisum(j,n)
            radii (r,n) = radii (j,n)
            ncell (r,n) = ncell (j,n)
            j = j+1
            if (j .gt. nr_irreg) then
               max_rcoord(n) = r
               exit
            end if
         end do
      end do

      ! compute phibar
      stencil_coord = 0

      !$OMP PARALLEL DO PRIVATE(r,radius,j,limit) FIRSTPRIVATE(stencil_coord)
      do r=0,nr_fine-1

         radius = (dble(r)+HALF)*dr(1)

         ! find the closest coordinate
         do j=stencil_coord(which_lev(r)),max_rcoord(which_lev(r))
            if (abs(radius-radii(j  ,which_lev(r))) .lt. &
                abs(radius-radii(j+1,which_lev(r)))) then
               stencil_coord(which_lev(r)) = j
               exit
            end if
         end do

         ! make sure the interpolation points will be in bounds
         if (which_lev(r) .ne. nlevs) then
            stencil_coord(which_lev(r)) = max(stencil_coord(which_lev(r)),1)
         end if
         stencil_coord(which_lev(r)) = min(stencil_coord(which_lev(r)), &
                                           max_rcoord(which_lev(r))-1)

         if (r > nr_fine - 1 - drdxfac*2.d0**(max_levs-1)) then
            limit = .false. 
         else 
            limit = .true.
         end if

         call quad_interp(radius, &
                          radii(stencil_coord(which_lev(r))-1,which_lev(r)), &
                          radii(stencil_coord(which_lev(r))  ,which_lev(r)), &
                          radii(stencil_coord(which_lev(r))+1,which_lev(r)), &
                          phibar(1,r), &
                          phisum(stencil_coord(which_lev(r))-1,which_lev(r)), &
                          phisum(stencil_coord(which_lev(r))  ,which_lev(r)), &
                          phisum(stencil_coord(which_lev(r))+1,which_lev(r)), limit)

      end do
      !$OMP END PARALLEL DO

   end if

   call destroy(bpt)

 contains

   subroutine cubic_interp(x,x0,x1,x2,x3,y,y0,y1,y2,y3)

     real(kind=dp_t), intent(in   ) :: x,x0,x1,x2,x3,y0,y1,y2,y3
     real(kind=dp_t), intent(  out) :: y

     y = y0 + (y1-y0)/(x1-x0)*(x-x0) &
          + ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)*(x-x0)*(x-x1) &
          + ( ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0) &
          -((y3-y2)/(x3-x2)-(y2-y1)/(x2-x1))/(x3-x1) ) / (x3-x0) &
          *(x-x0)*(x-x1)*(x-x2)

     if (y .gt. max(y0,y1,y2,y3)) y = max(y0,y1,y2,y3)
     if (y .lt. min(y0,y1,y2,y3)) y = min(y0,y1,y2,y3)

   end subroutine cubic_interp

   subroutine quad_interp(x,x0,x1,x2,y,y0,y1,y2,limit)

     real(kind=dp_t), intent(in   ) :: x,x0,x1,x2,y0,y1,y2
     real(kind=dp_t), intent(  out) :: y
     logical,         intent(in   ) :: limit

     y = y0 + (y1-y0)/(x1-x0)*(x-x0) &
          + ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)*(x-x0)*(x-x1)


     if (limit) then
        if (y .gt. max(y0,y1,y2)) y = max(y0,y1,y2)
        if (y .lt. min(y0,y1,y2)) y = min(y0,y1,y2)
     end if

   end subroutine quad_interp

   subroutine lin_interp(x,x0,x1,y,y0,y1)

     real(kind=dp_t), intent(in   ) :: x,x0,x1,y0,y1
     real(kind=dp_t), intent(  out) :: y

     y = y0 + (y1-y0)/(x1-x0)*(x-x0)

   end subroutine lin_interp

 end subroutine average

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine average_irreg(mla,phi,phibar_irreg,dx,incomp)

    use geometry, only: nr_irreg, spherical
    use bl_prof_module, only: bl_prof_timer, build, destroy
    use bl_constants_module, only: ZERO
    use restrict_base_module, only: restrict_base, fill_ghost_base

    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: incomp
    type(multifab) , intent(in   ) :: phi(:)
    real(kind=dp_t), intent(inout) :: phibar_irreg(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    ! local
    real(kind=dp_t), pointer     :: pp(:,:,:,:)
    logical, pointer             :: mp(:,:,:,:)

    integer                      :: lo(mla%dim),hi(mla%dim),dm,nlevs
    integer                      :: i,r,n,ng

    integer, allocatable ::  ncell_proc(:,:)
    integer, allocatable ::       ncell(:,:)

    real(kind=dp_t), allocatable :: phisum_proc(:,:)
    real(kind=dp_t), allocatable ::      phisum(:,:)
    real(kind=dp_t), allocatable :: radii(:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "average_irreg")

    dm = mla%dim
    nlevs = mla%nlevel

    if (spherical .eq. 0) then
       call bl_error("average_irreg only written for spherical")
    end if
    
    allocate(ncell_proc ( 0:nr_irreg  ,nlevs))
    allocate(ncell      (-1:nr_irreg  ,nlevs))
    allocate(phisum_proc( 0:nr_irreg  ,nlevs))
    allocate(phisum     (-1:nr_irreg  ,nlevs))
    allocate(radii      (-1:nr_irreg+1,nlevs))
    !
    ! radii contains every possible distance that a cell-center at the finest
    ! level can map into
    !
    !$OMP PARALLEL PRIVATE(r,n)
    do n=1,nlevs
       !$OMP DO
       do r=0,nr_irreg
          radii(r,n) = sqrt(0.75d0+2.d0*r)*dx(n,1)
       end do
       !$OMP END DO NOWAIT
    end do
    !$OMP END PARALLEL

    radii(nr_irreg+1,:) = 1.d99

    ng = nghost(phi(1))

    phibar_irreg = ZERO
    ncell        = 0
    ncell_proc   = 0
    phisum       = ZERO       
    phisum_proc  = ZERO
    !
    ! For spherical, we construct a 1D array, phisum, that has space
    ! allocated for every possible radius that a cell-center at the finest
    ! level can map into.  The radius locations have been precomputed and stored
    ! in radii(:).
    
    ! For cells at the non-finest level, map a weighted contribution into the nearest
    ! bin in phisum.
    !
    do n=nlevs,1,-1

       do i=1, nfabs(phi(n))
          pp => dataptr(phi(n), i)
          lo =  lwb(get_box(phi(n), i))
          hi =  upb(get_box(phi(n), i))

          if (n .eq. nlevs) then
             call sum_phi_3d_sphr(radii(0:,n),nr_irreg,pp(:,:,:,:),phisum_proc(:,n), &
                                  lo,hi,ng,dx(n,:),ncell_proc(:,n),incomp)
          else
             ! we include the mask so we don't double count; i.e., we only consider
             ! cells that we can "see" when constructing the sum
             mp => dataptr(mla%mask(n), i)
             call sum_phi_3d_sphr(radii(0:,n),nr_irreg,pp(:,:,:,:),phisum_proc(:,n), &
                                  lo,hi,ng,dx(n,:),ncell_proc(:,n),incomp, &
                                  mp(:,:,:,1))
          end if
       end do

       call parallel_reduce(ncell (0:,n), ncell_proc (:,n), MPI_SUM)
       call parallel_reduce(phisum(0:,n), phisum_proc(:,n), MPI_SUM)
    end do

    ! compute phibar_irreg
    do n=1,nlevs
       do r=0,nr_irreg
          if (ncell(r,n) .ne. 0.d0) then
             phibar_irreg(n,r) = phisum(r,n) / dble(ncell(r,n))
          end if
       end do
    end do

   call destroy(bpt)

 end subroutine average_irreg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sum_phi_1d(phi,phisum,lo,hi,ng,incomp)

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,:)
    real (kind=dp_t), intent(inout) :: phisum(0:)

    integer :: i

    do i=lo(1),hi(1)
       phisum(i) = phi(i,incomp)
    end do

  end subroutine sum_phi_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sum_phi_2d(phi,phisum,lo,hi,ng,incomp)

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(inout) :: phisum(0:)

    integer :: i,j

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          phisum(j) = phisum(j) + phi(i,j,incomp)
       end do
    end do

  end subroutine sum_phi_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sum_phi_3d(phi,phisum,lo,hi,ng,incomp)

    integer         , intent(in   ) :: lo(:), hi(:), ng, incomp
    real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout) :: phisum(0:)

    integer :: i,j,k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             phisum(k) = phisum(k) + phi(i,j,k,incomp)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine sum_phi_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sum_phi_3d_sphr(radii,nr_irreg,phi,phisum,lo,hi,ng,dx,ncell,incomp,mask)

    use geometry, only: center
    use probin_module, only: prob_lo
    use bl_constants_module, only: HALF

    integer         , intent(in   )           :: lo(:), hi(:), ng, incomp, nr_irreg
    real (kind=dp_t), intent(in   )           :: radii(0:)
    real (kind=dp_t), intent(in   )           :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout)           :: phisum(0:)
    real (kind=dp_t), intent(in   )           :: dx(:)
    integer         , intent(inout)           :: ncell(0:)
    logical         , intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    ! local
    real (kind=dp_t) :: x, y, z, radius
    integer          :: i, j, k, index
    logical          :: cell_valid

    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k) + HALF)*dx(3) - center(3)
       
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j) + HALF)*dx(2) - center(2)
          
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i) + HALF)*dx(1) - center(1)
             
             ! make sure the cell isn't covered by finer cells
             cell_valid = .true.
             if ( present(mask) ) then
                if ( (.not. mask(i,j,k)) ) cell_valid = .false.
             end if
             
             if (cell_valid) then

                ! compute distance to center
                radius = sqrt(x**2 + y**2 + z**2)
                
                ! figure out which radii index this point maps into
                index = ((radius / dx(1))**2 - 0.75d0) / 2.d0
                
                ! due to roundoff error, need to ensure that we are in the proper radial bin
                if (index .lt. nr_irreg) then
                   if (abs(radius-radii(index)) .gt. abs(radius-radii(index+1))) then
                      index = index+1
                   end if
                end if
                
                phisum(index) = phisum(index) + phi(i,j,k,incomp)
                ncell(index)  = ncell(index) + 1

             end if
             
          end do
       end do
    end do

  end subroutine sum_phi_3d_sphr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine average_one_level(n,phi,phibar,incomp)

    use geometry, only: nr_fine, nr, spherical
    use bl_prof_module, only: bl_prof_timer, build, destroy
    use bl_constants_module, only: ZERO

    integer        , intent(in   ) :: n,incomp
    type(multifab) , intent(in   ) :: phi(:)
    real(kind=dp_t), intent(inout) :: phibar(:,0:)

    ! local
    real(kind=dp_t), pointer :: pp(:,:,:,:)
    type(box)                :: domain
    integer                  :: domlo(get_dim(phi(1))),domhi(get_dim(phi(1)))
    integer                  :: lo(get_dim(phi(1))),hi(get_dim(phi(1)))
    integer                  :: i,r,ng,ncell,dm

    real(kind=dp_t) ::   phisum_proc(0:nr_fine-1)
    real(kind=dp_t) ::        phisum(0:nr_fine-1)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "average_one_level")

    dm = get_dim(phi(1))
    ng = nghost(phi(1))

    phibar = ZERO

    ncell        = 0
    phisum       = ZERO       
    phisum_proc  = ZERO

    if (spherical .eq. 0) then
       
       domain = get_pd(get_layout(phi(n)))
       domlo  = lwb(domain)
       domhi  = upb(domain)

       if (dm .eq. 1) then
          ncell = 1
       else if (dm .eq. 2) then
          ncell = domhi(1)-domlo(1)+1
       else if (dm .eq. 3) then
          ncell = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
       end if

       do i=1, nfabs(phi(n))
          pp => dataptr(phi(n), i)
          lo =  lwb(get_box(phi(n), i))
          hi =  upb(get_box(phi(n), i))
          select case (dm)
          case (1)
             call sum_phi_1d(pp(:,1,1,:),phisum_proc,lo,hi,ng,incomp)
          case (2)
             call sum_phi_2d(pp(:,:,1,:),phisum_proc,lo,hi,ng,incomp)
          case (3)
             call sum_phi_3d(pp(:,:,:,:),phisum_proc,lo,hi,ng,incomp)
          end select
       end do

       call parallel_reduce(phisum, phisum_proc, MPI_SUM)

       do r=0,nr(n)-1
          phibar(n,r) = phisum(r) / dble(ncell)
       end do

    else if(spherical .eq. 1) then

       call bl_error("average_one_level not written for multilevel spherical")

    end if

    call destroy(bpt)

  end subroutine average_one_level

end module average_module
