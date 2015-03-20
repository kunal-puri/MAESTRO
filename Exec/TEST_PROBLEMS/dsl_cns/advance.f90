! advance_timestep advances the solution through one timestep,
! proceeding through the 12 steps described in the multilevel paper
! (Nonaka et al. 2010).

module advance_timestep_module

  use bl_types           , only: dp_t
  use bl_constants_module, only: ZERO, HALF, TWO
  use multifab_module 
  use ml_layout_module   , only: ml_layout
  use define_bc_module   , only: bc_tower
  use parallel           , only: parallel_IOProcessor, parallel_IOProcessorNode, &
                                 parallel_wtime, parallel_reduce, parallel_barrier, &
                                 MPI_MAX
  implicit none

  private
  public :: advance_timestep

contains
    
  subroutine advance_timestep(mla,uold,sold,unew,snew,dx,dt,&
                              dtold, the_bc_tower)

    use bl_prof_module  , only : bl_prof_timer, build, destroy
    use variables       , only : nscal, temp_comp, rho_comp, rhoh_comp, pi_comp
    use probin_module   , only : barrier_timers, verbose, ppm_trace_forces
    use time_module     , only : time
    
    logical,         intent(in   ) :: init_mode
    type(ml_layout), intent(inout) :: mla
    type(multifab),  intent(in   ) :: uold(:)
    type(multifab),  intent(inout) :: sold(:)
    type(multifab),  intent(inout) :: unew(:)
    type(multifab),  intent(inout) :: snew(:)

    real(dp_t)    ,  intent(in   ) :: dx(:,:),dt,dtold
    type(bc_tower),  intent(in   ) :: the_bc_tower

    ! local
    type(multifab) ::             rhohalf(mla%nlevel)
    type(multifab) ::       w0_force_cart(mla%nlevel)
    type(multifab) ::              macrhs(mla%nlevel)
    type(multifab) ::              macphi(mla%nlevel)
    type(multifab) ::           hgrhs_old(mla%nlevel)
    type(multifab) ::          Source_nph(mla%nlevel)
    type(multifab) ::            thermal1(mla%nlevel)
    type(multifab) ::              s2star(mla%nlevel)
    type(multifab) ::                  s1(mla%nlevel)
    type(multifab) ::                  s2(mla%nlevel)
    type(multifab) ::   delta_gamma1_term(mla%nlevel)
    type(multifab) ::        delta_gamma1(mla%nlevel)
    type(multifab) ::       rho_omegadot1(mla%nlevel)
    type(multifab) ::           rho_Hnuc1(mla%nlevel)
    type(multifab) ::        div_coeff_3d(mla%nlevel)
    type(multifab) :: div_coeff_cart_edge(mla%nlevel,mla%dim)
    type(multifab) ::              gamma1(mla%nlevel)
    type(multifab) ::          etarhoflux(mla%nlevel)
    type(multifab) ::                peos(mla%nlevel)
    type(multifab) ::        peosbar_cart(mla%nlevel)
    type(multifab) ::        delta_p_term(mla%nlevel)
    type(multifab) ::             Tcoeff1(mla%nlevel)
    type(multifab) ::             hcoeff1(mla%nlevel)
    type(multifab) ::            Xkcoeff1(mla%nlevel)
    type(multifab) ::             pcoeff1(mla%nlevel)
    type(multifab) ::             Tcoeff2(mla%nlevel)
    type(multifab) ::             hcoeff2(mla%nlevel)
    type(multifab) ::            Xkcoeff2(mla%nlevel)
    type(multifab) ::             pcoeff2(mla%nlevel)
    type(multifab) ::          scal_force(mla%nlevel)
    type(multifab) ::               pi_cc(mla%nlevel)
    type(multifab) ::           delta_chi(mla%nlevel)

    type(multifab) ::               w0mac(mla%nlevel,mla%dim)
    type(multifab) ::                umac(mla%nlevel,mla%dim)
    type(multifab) ::               sedge(mla%nlevel,mla%dim)
    type(multifab) ::               sflux(mla%nlevel,mla%dim)

    type(multifab) ::               ustar(mla%nlevel,mla%dim)

    real(kind=dp_t), allocatable ::        grav_cell_nph(:,:)
    real(kind=dp_t), allocatable ::        grav_cell_new(:,:)
    real(kind=dp_t), allocatable ::             rho0_nph(:,:)
    real(kind=dp_t), allocatable ::               p0_nph(:,:)
    real(kind=dp_t), allocatable ::     p0_minus_peosbar(:,:)
    real(kind=dp_t), allocatable ::              peosbar(:,:)
    real(kind=dp_t), allocatable ::             w0_force(:,:)
    real(kind=dp_t), allocatable ::                 Sbar(:,:)
    real(kind=dp_t), allocatable ::        div_coeff_nph(:,:)
    real(kind=dp_t), allocatable ::        gamma1bar_old(:,:)
    real(kind=dp_t), allocatable ::      gamma1bar_temp1(:,:)
    real(kind=dp_t), allocatable ::      gamma1bar_temp2(:,:)
    real(kind=dp_t), allocatable :: delta_gamma1_termbar(:,:)
    real(kind=dp_t), allocatable ::               w0_old(:,:)
    real(kind=dp_t), allocatable ::       div_coeff_edge(:,:)
    real(kind=dp_t), allocatable ::  rho0_predicted_edge(:,:)
    real(kind=dp_t), allocatable ::         delta_chi_w0(:,:)

    integer    :: i,n,comp,proj_type,nlevs,dm
    real(dp_t) :: halfdt

    ! need long int to store numbers greater than 2^31
    integer(kind=ll_t) :: numcell

    real(kind=dp_t) :: advect_time , advect_time_start , advect_time_max
    real(kind=dp_t) :: macproj_time, macproj_time_start, macproj_time_max
    real(kind=dp_t) :: ndproj_time , ndproj_time_start , ndproj_time_max
    real(kind=dp_t) :: thermal_time, thermal_time_start, thermal_time_max
    real(kind=dp_t) :: react_time  , react_time_start  , react_time_max
    real(kind=dp_t) :: misc_time   , misc_time_start   , misc_time_max

    integer :: nreduce
    real(kind=dp_t), allocatable :: times_local(:), times_global(:)
    
    type(bl_prof_timer), save :: bpt

    call build(bpt, "advance_timestep")

    nlevs = mla%nlevel
    dm = mla%dim

    advect_time  = 0.d0
    macproj_time = 0.d0
    ndproj_time  = 0.d0
    thermal_time = 0.d0
    react_time   = 0.d0
    misc_time    = 0.d0

    misc_time_start = parallel_wtime()

    if (verbose .ge. 1) then

       if (parallel_IOProcessor()) then
          do n = 1, nlevs
             print *, 'level: ', n
             !print *, '   number of boxes = ', nboxes(pi(n)%la)
             print *, '   maximum zones   = ', (extent(mla%mba%pd(n),i),i=1,dm)
          end do
       end if

       do n=1,nlevs
          numcell = multifab_volume(pi(n),.false.)
          if (parallel_IOProcessor()) then
             print*,'Number of valid cells at level        ',n,numcell
          end if
          numcell = multifab_volume(pi(n),.true.)
          if (parallel_IOProcessor()) then
             print*,'Number of valid + ghost cells at level',n,numcell
          end if
       end do
       if (parallel_IOProcessor()) then
          print*,''
       end if
    end if

    halfdt = half*dt

    if (barrier_timers) call parallel_barrier()
    misc_time = misc_time + parallel_wtime() - misc_time_start

    do n = 1, nlevs
       call setval(sold(n), ONE, all=.true.)
    end do
    
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !! STEP 3 -- construct the advective velocity
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!     if (parallel_IOProcessor() .and. verbose .ge. 1) then
!        write(6,*) '<<< STEP  3 : create MAC velocities>>> '
!     end if

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     !! Construct the face-centered velocities
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     do n=1,nlevs
!        do comp=1,dm
!           call multifab_build_edge(ustar(n,comp), mla%la(n),1,1,comp)
!           call setval(ustar(n,comp), ZERO, all=.true.)
!        end do
!     end do
    
!     ! Compute the non-divergence free velocities
!     call compute_vstar(uold,sold,ustar,dx,dt,the_bc_tower%bc_tower_array,mla)

!     do n=1,nlevs
!        call multifab_build(delta_gamma1_term(n), mla%la(n), 1, 0)
!        call multifab_build(macrhs(n),            mla%la(n), 1, 0)
!        call multifab_build(delta_chi(n),         mla%la(n), 1, 0)
!        call setval(delta_gamma1_term(n), ZERO, all=.true.)
!        call setval(delta_chi(n),         ZERO, all=.true.)
!     end do

!     if (barrier_timers) call parallel_barrier()
!     advect_time = advect_time + parallel_wtime() - advect_time_start

!     macproj_time_start = parallel_wtime()

!     ! Additional RHS for the Projection equation
!     call make_macrhs(macrhs,rho0_old,Source_nph,delta_gamma1_term,Sbar,div_coeff_old,dx, &
!                      gamma1bar_old,p0_old,delta_p_term,dt,delta_chi,.true.)

!     do n=1,nlevs
!        call destroy(delta_gamma1_term(n))
!        call destroy(Source_nph(n))
!        call destroy(delta_p_term(n))
!     end do

!     ! Phi is the variable used in the velocity projection step
!     do n=1,nlevs
!        call multifab_build(macphi(n), mla%la(n), 1, 1)
!        call setval(macphi(n), ZERO, all=.true.)
!     end do

!     ! Velocity Projection
!     call cell_to_edge(div_coeff_old,div_coeff_edge)
!     call macproject(mla,ustar,macphi,sold,dx,dt,the_bc_tower)

!     do n=1,nlevs
!        call destroy(macrhs(n))
!     end do

!     if (barrier_timers) call parallel_barrier()
!     macproj_time = macproj_time + (parallel_wtime() - macproj_time_start)

  end subroutine advance_timestep

end module advance_timestep_module
