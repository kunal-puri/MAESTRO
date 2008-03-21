module estdt_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use fill_3d_module
  use geometry
  use variables

  implicit none

contains

   subroutine estdt (istep, u, s, force, divU, normal, w0, p0, gam1, dx, cflfac, dtold, dt)

      integer        , intent(in ) :: istep
      type(multifab) , intent(in ) :: u
      type(multifab) , intent(in ) :: s
      type(multifab) , intent(in ) :: force
      type(multifab) , intent(in ) :: divU
      type(multifab) , intent(in ) :: normal
      real(kind=dp_t), intent(in ) :: w0(:), p0(:), gam1(:)
      real(kind=dp_t), intent(in ) :: dx(:)
      real(kind=dp_t), intent(in ) :: cflfac, dtold
      real(kind=dp_t), intent(out) :: dt

      real(kind=dp_t), pointer:: uop(:,:,:,:)
      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: fp(:,:,:,:)
      real(kind=dp_t), pointer:: np(:,:,:,:)
      real(kind=dp_t), pointer:: dUp(:,:,:,:)     
      integer :: lo(u%dim),hi(u%dim),ng,dm
      real(kind=dp_t) :: dt_grid,dt_proc
      real(kind=dp_t) :: dtchange
      integer         :: i

      real(kind=dp_t), parameter :: rho_min = 1.d-20

      ng = u%ng
      dm = u%dim

      dtchange = 1.1d0
      dt_grid  = 1.d20
      dt_proc  = 1.d20

      do i = 1, u%nboxes
         if ( multifab_remote(u, i) ) cycle
         uop => dataptr(u, i)
         sop => dataptr(s, i)
          fp => dataptr(force, i)
         dUp => dataptr(divU, i)

         lo =  lwb(get_box(u, i))
         hi =  upb(get_box(u, i))

         select case (dm)
            case (2)
              call estdt_2d(uop(:,:,1,:), sop(:,:,1,:), fp(:,:,1,:), dUp(:,:,1,1), &
                            w0, p0, gam1, lo, hi, ng, dx, rho_min, dt_grid)
            case (3)
              np => dataptr(normal, i)
              call estdt_3d(uop(:,:,:,:), sop(:,:,:,:), fp(:,:,:,:), dUp(:,:,:,1), np(:,:,:,:), &
                            w0, p0, gam1, lo, hi, ng, dx, rho_min, dt_grid)
         end select

         dt_proc = min(dt_proc,dt_grid)
      end do

      ! This sets dt to be the min of dt_proc over all processors.
      call parallel_reduce(dt,dt_proc,MPI_MIN)

      dt = dt * cflfac

      print *, 'advective dt = ', dt

      if (dtold .gt. 0.0D0 ) dt = min(dt,dtchange*dtold)

!     dt = 0.020_dp_t*(dx(1)/1562500.0_dp_t)
      dt = 0.005_dp_t*(dx(1)/ 281250.0_dp_t)
!     dt = 0.01_dp_t*(dx(1)/ 281250.0_dp_t)

      print *,'Computing dt at istep ',istep,' to be ',dt

   end subroutine estdt


   subroutine estdt_2d (u, s, force, divU, w0, p0, gam1, lo, hi, ng, dx, rho_min, dt)

     integer, intent(in) :: lo(:), hi(:), ng
     real (kind = dp_t), intent(in ) ::     u(lo(1)-ng:,lo(2)-ng:,:)  
     real (kind = dp_t), intent(in ) ::     s(lo(1)-ng:,lo(2)-ng:,:)  
     real (kind = dp_t), intent(in ) :: force(lo(1)- 1:,lo(2)- 1:,:)  
     real (kind = dp_t), intent(in ) ::  divU(lo(1):,lo(2):)
     real (kind = dp_t), intent( in) ::   w0(0:), p0(lo(2):), gam1(lo(2):)
     real (kind = dp_t), intent(in ) :: dx(:)
     real (kind = dp_t), intent(in ) :: rho_min
     real (kind = dp_t), intent(out) :: dt

!    Local variables
     real (kind = dp_t)  :: spdx, spdy, spdr
     real (kind = dp_t)  :: pforcex, pforcey
     real (kind = dp_t)  :: eps
     real (kind = dp_t)  :: dt_divU, denom, gradp0
     integer             :: i,j


     eps = 1.0e-8

     ! advective constraints
     spdx  = 0.0D0 
     spdy  = 0.0D0 

     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           spdx = max(spdx ,abs(u(i,j,1)))
           spdy = max(spdy ,abs(u(i,j,2)+w0(j)))
        enddo
     enddo

     spdx = spdx / dx(1)
     spdy = spdy / dx(2)

     spdr = ZERO 
     do j = lo(2),hi(2)
        spdr = max(spdr ,abs(w0(j)))
     enddo
     spdr = spdr / dx(2)

     if (spdx < eps .and. spdy < eps .and. spdr < eps) then
        dt = min(dx(1),dx(2))
     else
        dt = 1.0D0  / max(spdx,spdy,spdr)
     endif


     ! force constraints
     pforcex = 0.0D0 
     pforcey = 0.0D0 
     
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           pforcex = max(pforcex,abs(force(i,j,1)))
           pforcey = max(pforcey,abs(force(i,j,2)))
        enddo
     enddo

     if (pforcex > eps) then
        dt = min(dt,sqrt(2.0D0 *dx(1)/pforcex))
     endif

     if (pforcey > eps) then
        dt = min(dt,sqrt(2.0D0 *dx(2)/pforcey))
     endif


     ! divU constraint
     dt_divU = 1.d30

     do j = lo(2), hi(2)
        
        if (j .eq. lo(2)) then
           gradp0 = (p0(j+1) - p0(j))/dx(2)
        else if (j .eq. hi(2)) then
           gradp0 = (p0(j) - p0(j-1))/dx(2)
        else
           gradp0 = HALF*(p0(j+1) - p0(j-1))/dx(2)
        endif

        do i = lo(1), hi(1)

           denom = divU(i,j) - u(i,j,2)*gradp0/(gam1(j)*p0(j))

           if (denom > ZERO) then

              dt_divU = min(dt_divU, &
                            HALF*(ONE - rho_min/s(i,j,rho_comp))/denom)
           endif

        enddo
     enddo
     print *, '%%% timesteps (source, advective): ', dt_divU, dt

     dt = min(dt,dt_divU)

   end subroutine estdt_2d

   subroutine estdt_3d (u, s, force, divU, normal, w0, p0, gam1, lo, hi, ng, dx, rho_min, dt)

     integer, intent(in) :: lo(:), hi(:), ng
     real (kind = dp_t), intent(in ) ::      u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
     real (kind = dp_t), intent(in ) ::      s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
     real (kind = dp_t), intent(in ) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
     real (kind = dp_t), intent(in ) ::   divU(lo(1):,lo(2):,lo(3):)
     real (kind = dp_t), intent(in ) :: normal(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
     real (kind = dp_t), intent( in) ::   w0(0:), p0(lo(3):), gam1(lo(3):)
     real (kind = dp_t), intent(in ) :: dx(:)
     real (kind = dp_t), intent(in ) :: rho_min
     real (kind = dp_t), intent(out) :: dt

!    Local variables
     real (kind = dp_t), allocatable :: w0_cart(:,:,:,:)
     real (kind = dp_t)  :: spdx, spdy, spdz, spdr
     real (kind = dp_t)  :: pforcex, pforcey, pforcez
     real (kind = dp_t)  :: eps
     real (kind = dp_t)  :: dt_divU
     integer             :: i,j,k

     eps = 1.0e-8

     spdx    = ZERO
     spdy    = ZERO 
     spdz    = ZERO 
     
     if (spherical .eq. 0) then

        ! Limit dt based on velocity terms
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 spdx = max(spdx ,abs(u(i,j,k,1)))
                 spdy = max(spdy ,abs(u(i,j,k,2)))
                 spdz = max(spdz ,abs(u(i,j,k,3)+w0(k)))
              enddo
           enddo
        enddo
        
        spdr = ZERO 
        do k = lo(3),hi(3)
           spdr = max(spdr ,abs(w0(k)))
        enddo
        spdr = spdr / dx(3)
        
     else

        allocate(w0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3))
        call put_w0_on_3d_cells(w0(0:),w0_cart,normal,lo,hi,dx,0)

        ! Limit dt based on velocity terms
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 spdx = max(spdx ,abs(u(i,j,k,1)+w0_cart(i,j,k,1)))
                 spdy = max(spdy ,abs(u(i,j,k,2)+w0_cart(i,j,k,2)))
                 spdz = max(spdz ,abs(u(i,j,k,3)+w0_cart(i,j,k,3)))
              enddo
           enddo
        enddo

        deallocate(w0_cart)

        spdr = ZERO 
        do k = 0,size(w0,dim=1)-1
           spdr = max(spdr ,abs(w0(k)))
        enddo
        spdr = spdr / dr
        
     end if

     spdx = spdx / dx(1)
     spdy = spdy / dx(2)
     spdz = spdz / dx(3)
      
     if (spdx < eps .and. spdy < eps .and. spdz < eps .and. spdr < eps) then
        dt = min(dx(1),dx(2),dx(3))
     else
        dt = 1.0D0  / max(spdx,spdy,spdz,spdr)
     endif
     

     ! Limit dt based on forcing terms
     pforcex = ZERO 
     pforcey = ZERO 
     pforcez = ZERO 

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              pforcex = max(pforcex,abs(force(i,j,k,1)))
              pforcey = max(pforcey,abs(force(i,j,k,2)))
              pforcez = max(pforcez,abs(force(i,j,k,3)))
           enddo
        enddo
     enddo
     
     if (pforcex > eps) then
        dt = min(dt,sqrt(2.0D0*dx(1)/pforcex))
     endif
     
     if (pforcey > eps) then
        dt = min(dt,sqrt(2.0D0*dx(2)/pforcey))
     endif
     
     if (pforcez > eps) then
        dt = min(dt,sqrt(2.0D0*dx(3)/pforcez))
     endif


     ! divU constraint
     dt_divU = 1.d30

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              if (divU(i,j,k) > ZERO) then
                 dt_divU = min(dt_divU, &
                               (ONE - rho_min/s(i,j,k,rho_comp))/divU(i,j,k))
              endif

           enddo
        enddo
     enddo

     dt = min(dt,dt_divU)

   end subroutine estdt_3d

end module estdt_module
