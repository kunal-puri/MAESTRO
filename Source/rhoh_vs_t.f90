module rhoh_vs_t_module

  use bl_types
  use bl_constants_module
  use variables
  use geometry
  use network
  use probin_module
  use eos_module

  implicit none

  contains

   subroutine makeRhoHfromT_2d (sx,sy,s0_old,s0_edge_old,s0_new,s0_edge_new,lo,hi)

    implicit none
    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sx(lo(1):,lo(2):,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1):,lo(2):,:)
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:), s0_edge_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(0:,:), s0_edge_new(0:,:)

    !     Local variables
    integer :: i, j, n, nr
    real(kind=dp_t) qreact

    nr = size(s0_old,dim=1)

    do_diag = .false.

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)+1
         sx(i,j,rho_comp) = 0.d0
         do n = 1,nspec
           sx(i,j,rho_comp) = sx(i,j,rho_comp) + sx(i,j,spec_comp+n-1)
         end do
       end do
    end do

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)+1

         temp_row(1) = sx(i,j,temp_comp)
          den_row(1) = sx(i,j,rho_comp) + HALF * (s0_old(j,rho_comp) + s0_new(j,rho_comp))
          xn_row(1,:) = (sx(i,j,spec_comp:spec_comp+nspec-1)  + &
                       HALF * ( s0_old(j,spec_comp:spec_comp+nspec-1) + &
                                s0_new(j,spec_comp:spec_comp+nspec-1) ) ) /den_row(1) 

         call eos(eos_input_rt, den_row, temp_row, &
                  npts, nspec, &
                  xn_row, &
                  p_row, h_row, e_row, &
                  cv_row, cp_row, xne_row, eta_row, pele_row, &
                  dpdt_row, dpdr_row, dedt_row, dedr_row, &
                  dpdX_row, dhdX_row, &
                  gam1_row, cs_row, s_row, &
                  dsdt_row, dsdr_row, &
                  do_diag)

         sx(i,j,rhoh_comp) = den_row(1)*h_row(1)

         qreact = 0.0d0
         if(use_big_h) then
            do n=1,nspec
               qreact = qreact + ebin(n)*xn_row(1,n)
            enddo
            sx(i,j,rhoh_comp) = sx(i,j,rhoh_comp) + den_row(1) * qreact
         endif

         sx(i,j,rhoh_comp) = sx(i,j,rhoh_comp) - HALF * (s0_old(j,rhoh_comp) + s0_new(j,rhoh_comp))
          
       enddo
    enddo

    do j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
         sy(i,j,rho_comp) = 0.d0
         do n = 1,nspec
           sy(i,j,rho_comp) = sy(i,j,rho_comp) + sy(i,j,spec_comp+n-1)
         end do
       end do
    end do

    do j = lo(2), hi(2)+1
       do i = lo(1), hi(1)

         temp_row(1) = sy(i,j,temp_comp)
          den_row(1) = sy(i,j,rho_comp) + HALF * (s0_edge_old(j,rho_comp) + s0_edge_new(j,rho_comp))
          xn_row(1,:) = (sy(i,j,spec_comp:spec_comp+nspec-1)  + &
                      HALF * ( s0_edge_old(j,spec_comp:spec_comp+nspec-1) + &
                               s0_edge_new(j,spec_comp:spec_comp+nspec-1) ) ) /den_row(1) 

         call eos(eos_input_rt, den_row, temp_row, &
                  npts, nspec, &
                  xn_row, &
                  p_row, h_row, e_row, &
                  cv_row, cp_row, xne_row, eta_row, pele_row, &
                  dpdt_row, dpdr_row, dedt_row, dedr_row, &
                  dpdX_row, dhdX_row, &
                  gam1_row, cs_row, s_row, &
                  dsdt_row, dsdr_row, &
                  do_diag)

         sy(i,j,rhoh_comp) = den_row(1)*h_row(1) 

         qreact = 0.0d0
         if(use_big_h) then
            do n=1,nspec
               qreact = qreact + ebin(n)*xn_row(1,n)
            enddo
            sy(i,j,rhoh_comp) = sy(i,j,rhoh_comp) + den_row(1) * qreact
         endif

         sy(i,j,rhoh_comp) = sy(i,j,rhoh_comp) - HALF * (s0_edge_old(j,rhoh_comp) + s0_edge_new(j,rhoh_comp))
          
       enddo
    enddo

   end subroutine makeRhoHfromT_2d

   subroutine makeRhoHfromT_3d (sx,sy,sz,s0_old,s0_edge_old,s0_new,s0_edge_new,lo,hi)

    implicit none
    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: sx(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(inout) :: sy(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(inout) :: sz(lo(1):,lo(2):,lo(3):,:)
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:), s0_edge_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(0:,:), s0_edge_new(0:,:)

    !     Local variables
    integer :: i, j, k, n, nr
    real(kind=dp_t) qreact

    nr = size(s0_old,dim=1)

    do_diag = .false.

    if (spherical .eq. 1) then
       print *,'MAKERHOHFROMT_3D NOT YET SET UP FOR SPHERICAL '
       stop
    end if

    do k = lo(3), hi(3)
     do j = lo(2), hi(2)
       do i = lo(1), hi(1)+1
         sx(i,j,k,rho_comp) = 0.d0
         do n = 1,nspec
           sx(i,j,k,rho_comp) = sx(i,j,k,rho_comp) + sx(i,j,k,spec_comp+n-1)
         end do
       end do
      end do
    end do

    do k = lo(3), hi(3)
     do j = lo(2), hi(2)
       do i = lo(1), hi(1)+1

         temp_row(1) = sx(i,j,k,temp_comp)
          den_row(1) = sx(i,j,k,rho_comp) + HALF * (s0_old(k,rho_comp) + s0_new(k,rho_comp))
          xn_row(1,:) = (sx(i,j,k,spec_comp:spec_comp+nspec-1)  + &
                       HALF * ( s0_old(k,spec_comp:spec_comp+nspec-1) + &
                                s0_new(k,spec_comp:spec_comp+nspec-1) ) ) /den_row(1) 

         call eos(eos_input_rt, den_row, temp_row, &
                  npts, nspec, &
                  xn_row, &
                  p_row, h_row, e_row, &
                  cv_row, cp_row, xne_row, eta_row, pele_row, &
                  dpdt_row, dpdr_row, dedt_row, dedr_row, &
                  dpdX_row, dhdX_row, &
                  gam1_row, cs_row, s_row, &
                  dsdt_row, dsdr_row, &
                  do_diag)

         sx(i,j,k,rhoh_comp) = den_row(1)*h_row(1)

         qreact = 0.0d0
         if(use_big_h) then
            do n=1,nspec
               qreact = qreact + ebin(n)*xn_row(1,n)
            enddo
            sx(i,j,k,rhoh_comp) = sx(i,j,k,rhoh_comp) + sx(i,j,k,rho_comp) * qreact
         endif

         sx(i,j,k,rhoh_comp) = sx(i,j,k,rhoh_comp) - HALF * (s0_old(k,rhoh_comp) + s0_new(k,rhoh_comp))
          
       enddo
      enddo
    enddo

    do k = lo(3), hi(3)
     do j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
         sy(i,j,k,rho_comp) = 0.d0
         do n = 1,nspec
           sy(i,j,k,rho_comp) = sy(i,j,k,rho_comp) + sy(i,j,k,spec_comp+n-1)
         end do
       end do
      end do
    end do

    do k = lo(3), hi(3)
     do j = lo(2), hi(2)+1
       do i = lo(1), hi(1)

         temp_row(1) = sy(i,j,k,temp_comp)
          den_row(1) = sy(i,j,k,rho_comp) + HALF * (s0_old(k,rho_comp) + s0_new(k,rho_comp))
          xn_row(1,:) = (sy(i,j,k,spec_comp:spec_comp+nspec-1)  + &
                       HALF * ( s0_old(k,spec_comp:spec_comp+nspec-1) + &
                                s0_new(k,spec_comp:spec_comp+nspec-1) ) ) /den_row(1) 

         call eos(eos_input_rt, den_row, temp_row, &
                  npts, nspec, &
                  xn_row, &
                  p_row, h_row, e_row, &
                  cv_row, cp_row, xne_row, eta_row, pele_row, &
                  dpdt_row, dpdr_row, dedt_row, dedr_row, &
                  dpdX_row, dhdX_row, &
                  gam1_row, cs_row, s_row, &
                  dsdt_row, dsdr_row, &
                  do_diag)

         sy(i,j,k,rhoh_comp) = den_row(1)*h_row(1)

         qreact = 0.0d0
         if(use_big_h) then
            do n=1,nspec
               qreact = qreact + ebin(n)*xn_row(1,n)
            enddo
            sy(i,j,k,rhoh_comp) = sy(i,j,k,rhoh_comp) + sy(i,j,k,rho_comp) * qreact
         endif

         sy(i,j,k,rhoh_comp) = sy(i,j,k,rhoh_comp) - HALF * (s0_old(k,rhoh_comp) + s0_new(k,rhoh_comp))
          
       enddo
      enddo
    enddo

    do k = lo(3), hi(3)+1
     do j = lo(2), hi(2)
       do i = lo(1), hi(1)
         sz(i,j,k,rho_comp) = 0.d0
         do n = 1,nspec
           sz(i,j,k,rho_comp) = sz(i,j,k,rho_comp) + sz(i,j,k,spec_comp+n-1)
         end do
       end do
      end do
    end do

    do k = lo(3), hi(3)+1
     do j = lo(2), hi(2)
       do i = lo(1), hi(1)

         temp_row(1) = sz(i,j,k,temp_comp)
          den_row(1) = sz(i,j,k,rho_comp) + HALF * (s0_edge_old(k,rho_comp) + s0_edge_new(k,rho_comp))
          xn_row(1,:) = (sz(i,j,k,spec_comp:spec_comp+nspec-1)  + &
                        HALF * ( s0_edge_old(k,spec_comp:spec_comp+nspec-1) + &
                                 s0_edge_new(k,spec_comp:spec_comp+nspec-1) ) ) /den_row(1)

         call eos(eos_input_rt, den_row, temp_row, &
                  npts, nspec, &
                  xn_row, &
                  p_row, h_row, e_row, &
                  cv_row, cp_row, xne_row, eta_row, pele_row, &
                  dpdt_row, dpdr_row, dedt_row, dedr_row, &
                  dpdX_row, dhdX_row, &
                  gam1_row, cs_row, s_row, &
                  dsdt_row, dsdr_row, &
                  do_diag)

         sz(i,j,k,rhoh_comp) = den_row(1)*h_row(1)

         qreact = 0.0d0
         if(use_big_h) then
            do n=1,nspec
               qreact = qreact + ebin(n)*xn_row(1,n)
            enddo
            sz(i,j,k,rhoh_comp) = sz(i,j,k,rhoh_comp) + sz(i,j,k,rho_comp) * qreact
         endif

         sz(i,j,k,rhoh_comp) = sz(i,j,k,rhoh_comp) - HALF * (s0_edge_old(k,rhoh_comp) + s0_edge_new(k,rhoh_comp))
          
       enddo
      enddo
    enddo

   end subroutine makeRhoHfromT_3d

  subroutine makeTfromRhoH_2d (state,lo,hi,ng,t0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) ::  state(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  t0(0:)

    !     Local variables
    integer :: i, j, n
    real(kind=dp_t) qreact

    do_diag = .false.

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          ! (rho, H) --> T, p

          den_row(1)  = state(i,j,rho_comp)
          temp_row(1) = t0(j)
          xn_row(1,:) = state(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

          qreact = 0.0d0
          if(use_big_h) then
             do n=1,nspec
                qreact = qreact + ebin(n)*xn_row(1,n)
             enddo
             h_row(1) = state(i,j,rhoh_comp) / state(i,j,rho_comp) - qreact
          else
             h_row(1) = state(i,j,rhoh_comp) / state(i,j,rho_comp)
          endif

          call eos(eos_input_rh, den_row, temp_row, &
                   npts, nspec, &
                   xn_row, &
                   p_row, h_row, e_row, &
                   cv_row, cp_row, xne_row, eta_row, pele_row, &
                   dpdt_row, dpdr_row, dedt_row, dedr_row, &
                   dpdX_row, dhdX_row, &
                   gam1_row, cs_row, s_row, &
                   dsdt_row, dsdr_row, &
                   do_diag)

          state(i,j,temp_comp) = temp_row(1)

       enddo
    enddo

  end subroutine makeTfromRhoH_2d

  subroutine makeTfromRhoH_3d (state,lo,hi,ng,t0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) ::  state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  t0(0:)

    !     Local variables
    integer :: i, j, k, n
    real(kind=dp_t) qreact

    do_diag = .false.

    do k = lo(3), hi(3)
     do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          ! (rho, H) --> T, p

          den_row(1)  = state(i,j,k,rho_comp)
          temp_row(1) = t0(k)
          xn_row(1,:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

          qreact = 0.0d0
          if(use_big_h) then
             do n=1,nspec
                qreact = qreact + ebin(n)*xn_row(1,n)
             enddo
             h_row(1) = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp) - qreact
          else
             h_row(1) = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)
          endif

          call eos(eos_input_rh, den_row, temp_row, &
                   npts, nspec, &
                   xn_row, &
                   p_row, h_row, e_row, &
                   cv_row, cp_row, xne_row, eta_row, pele_row, &
                   dpdt_row, dpdr_row, dedt_row, dedr_row, &
                   dpdX_row, dhdX_row, &
                   gam1_row, cs_row, s_row, &
                   dsdt_row, dsdr_row, &
                   do_diag)

          state(i,j,k,temp_comp) = temp_row(1)

       enddo
      enddo
    enddo

  end subroutine makeTfromRhoH_3d

end module rhoh_vs_t_module
