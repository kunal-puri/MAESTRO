! advance_timestep advances the solution through one timestep,
! proceeding through the 12 steps described in the multilevel paper
! (Nonaka et al. 2010).

module advance_timestep_module

  use bl_types           , only: dp_t
  use bl_constants_module, only: ZERO, HALF, TWO
  use multifab_module !   , only: multifab, multifab_build, multifab_build_edge, &
                      !           destroy, setval, nghost, &
                      !           extent, multifab_volume, nboxes, &
                      !           multifab_copy, multifab_copy_c, &
                      !           multifab_sub_sub, multifab_div_div_s, multifab_plus_plus
  use ml_layout_module   , only: ml_layout
  use define_bc_module   , only: bc_tower
  use parallel           , only: parallel_IOProcessor, parallel_IOProcessorNode, &
                                 parallel_wtime, parallel_reduce, parallel_barrier, &
                                 MPI_MAX
  use particle_module

  implicit none

  private

  public :: advance_timestep

contains
    
  subroutine advance_timestep(init_mode,mla,uold,sold,unew,snew, &
                              gpi,pi,normal,rho0_old,rhoh0_old, &
                              rho0_new,rhoh0_new,p0_old,p0_new,tempbar,gamma1bar,w0, &
                              rho_omegadot2,rho_Hnuc2,rho_Hext,thermal2,&
                              div_coeff_old,div_coeff_new, &
                              grav_cell_old,dx,dt,dtold,the_bc_tower, &
                              dSdt,Source_old,Source_new,etarho_ec,etarho_cc, &
                              psi,sponge,hgrhs,tempbar_init,particles)

    use bl_prof_module              , only : bl_prof_timer, build, destroy
    use      pre_advance_module     , only : advance_premac
    use velocity_advance_module     , only : velocity_advance
    use  density_advance_module     , only : density_advance
    use enthalpy_advance_module     , only : enthalpy_advance
    use make_div_coeff_module       , only : make_div_coeff
    use make_w0_module              , only : make_w0
    use advect_base_module          , only : advect_base_dens, advect_base_enthalpy
    use react_state_module          , only : react_state
    use make_S_module               , only : make_S
    use average_module              , only : average
    use phihalf_module              , only : make_S_at_halftime, make_at_halftime
    use extraphalf_module           , only : extrap_to_halftime
    use thermal_conduct_module      , only : thermal_conduct
    use make_explicit_thermal_module, only : make_explicit_thermal, make_thermal_coeffs 
    use make_grav_module            , only : make_grav_cell
    use make_eta_module             , only : make_etarho_planar, make_etarho_spherical
    use make_psi_module             , only : make_psi_planar, make_psi_spherical
    use fill_3d_module              , only : put_1d_array_on_cart, make_w0mac, make_s0mac
    use cell_to_edge_module         , only : cell_to_edge
    use make_gamma_module           , only : make_gamma
    use rhoh_vs_t_module            , only : makePfromRhoH, makeTfromRhoP, makeTfromRhoH
    use diag_module                 , only : diag
    use sanity_module               , only : sanity_check
    use enforce_HSE_module          , only : enforce_HSE

    use macrhs_module               , only : make_macrhs
    use macproject_module           , only : macproject

    use hgrhs_module                , only : make_hgrhs, correct_hgrhs
    use hgproject_module            , only : hgproject
    use proj_parameters             , only : pressure_iters_comp, regular_timestep_comp

    use variables                   , only : nscal, temp_comp, rho_comp, rhoh_comp, pi_comp, &
                                             foextrap_comp
    use geometry                    , only : nlevs_radial, spherical, nr_fine, compute_cutoff_coords
    use network                     , only : nspec
    use probin_module               , only : barrier_timers, evolve_base_state, fix_base_state, &
                                             use_etarho, dpdt_factor, verbose, &
                                             use_tfromp, use_thermal_diffusion, &
                                             use_delta_gamma1_term, nodal, mach_max_abort, &
                                             prob_lo, prob_hi, use_particles, ppm_trace_forces
    use time_module                 , only : time
    use addw0_module                , only : addw0
    use make_pi_cc_module           , only : make_pi_cc
    
    logical,         intent(in   ) :: init_mode
    type(ml_layout), intent(inout) :: mla
    type(multifab),  intent(in   ) ::   uold(:)
    type(multifab),  intent(inout) ::   sold(:)
    type(multifab),  intent(inout) ::   unew(:)
    type(multifab),  intent(inout) ::   snew(:)
    type(multifab),  intent(inout) ::  gpi(:)
    type(multifab),  intent(inout) ::   pi(:)
    type(multifab),  intent(in   ) :: normal(:)
    real(dp_t)    ,  intent(inout) ::  rho0_old(:,0:)
    real(dp_t)    ,  intent(inout) :: rhoh0_old(:,0:)
    real(dp_t)    ,  intent(inout) ::  rho0_new(:,0:)
    real(dp_t)    ,  intent(inout) :: rhoh0_new(:,0:)
    real(dp_t)    ,  intent(inout) ::    p0_old(:,0:)
    real(dp_t)    ,  intent(inout) ::    p0_new(:,0:)
    real(dp_t)    ,  intent(inout) ::   tempbar(:,0:)
    real(dp_t)    ,  intent(inout) ::   tempbar_init(:,0:)
    real(dp_t)    ,  intent(inout) :: gamma1bar(:,0:)
    real(dp_t)    ,  intent(inout) ::        w0(:,0:)
    type(multifab),  intent(inout) :: rho_omegadot2(:)
    type(multifab),  intent(inout) :: rho_Hnuc2(:)
    type(multifab),  intent(inout) :: rho_Hext(:)
    type(multifab),  intent(inout) ::  thermal2(:)
    real(dp_t)    ,  intent(inout) :: div_coeff_old(:,0:)
    real(dp_t)    ,  intent(inout) :: div_coeff_new(:,0:)
    real(dp_t)    ,  intent(inout) :: grav_cell_old(:,0:)
    real(dp_t)    ,  intent(in   ) :: dx(:,:),dt,dtold
    type(bc_tower),  intent(in   ) :: the_bc_tower
    type(multifab),  intent(inout) ::       dSdt(:)
    type(multifab),  intent(inout) :: Source_old(:)
    type(multifab),  intent(inout) :: Source_new(:)
    real(dp_t)    ,  intent(inout) ::  etarho_ec(:,0:)
    real(dp_t)    ,  intent(inout) ::  etarho_cc(:,0:)
    real(dp_t)    ,  intent(inout) ::        psi(:,0:)
    type(multifab),  intent(in   ) :: sponge(:)
    type(multifab),  intent(inout) ::  hgrhs(:)
    type(particle_container), intent(inout) :: particles

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

    allocate(       grav_cell_nph(nlevs_radial,0:nr_fine-1))
    allocate(       grav_cell_new(nlevs_radial,0:nr_fine-1))
    allocate(            rho0_nph(nlevs_radial,0:nr_fine-1))
    allocate(              p0_nph(nlevs_radial,0:nr_fine-1))
    allocate(    p0_minus_peosbar(nlevs_radial,0:nr_fine-1))
    allocate(             peosbar(nlevs_radial,0:nr_fine-1))
    allocate(            w0_force(nlevs_radial,0:nr_fine-1))
    allocate(                Sbar(nlevs_radial,0:nr_fine-1))
    allocate(       div_coeff_nph(nlevs_radial,0:nr_fine-1))
    allocate(       gamma1bar_old(nlevs_radial,0:nr_fine-1))
    allocate(     gamma1bar_temp1(nlevs_radial,0:nr_fine-1))
    allocate(     gamma1bar_temp2(nlevs_radial,0:nr_fine-1))
    allocate(delta_gamma1_termbar(nlevs_radial,0:nr_fine-1))
    allocate(              w0_old(nlevs_radial,0:nr_fine))
    allocate(      div_coeff_edge(nlevs_radial,0:nr_fine))
    allocate( rho0_predicted_edge(nlevs_radial,0:nr_fine))
    allocate(        delta_chi_w0(nlevs_radial,0:nr_fine-1))

    if (verbose .ge. 1) then

       if (parallel_IOProcessor()) then
          do n = 1, nlevs
             print *, 'level: ', n
             print *, '   number of boxes = ', nboxes(pi(n)%la)
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

    ! Initialize these to previous values
    w0_old        = w0
    gamma1bar_old = gamma1bar

    halfdt = half*dt

    if (barrier_timers) call parallel_barrier()
    misc_time = misc_time + parallel_wtime() - misc_time_start

    do n = 1, nlevs
       call setval(sold(n), ONE, all=.true.)
    end do
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 1 -- react the full state and then base state through dt/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    react_time_start = parallel_wtime()

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  1 : react state     '
    end if

    do n=1,nlevs
       call multifab_build(s1(n),            mla%la(n), nscal, nghost(sold(n)))
       call multifab_build(rho_omegadot1(n), mla%la(n), nspec, 0)
       call multifab_build(rho_Hnuc1(n),     mla%la(n), 1,     0)
       call setval( s1(n), ZERO, all=.true. )
    end do

    do n=1,nlevs
       call destroy(rho_omegadot1(n))
       call destroy(rho_Hnuc1(n))
       call destroy(rho_Hext(n))
    end do

    if (barrier_timers) call parallel_barrier()
    react_time = react_time + parallel_wtime() - react_time_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 2 -- define average expansion at time n+1/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    advect_time_start = parallel_wtime()
    
    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< CALLING advance_timestep with dt =',dt 
       write(6,*) '<<< STEP  2 : make w0 >>> ', dpdt_factor, evolve_base_state
    end if
    
    do n=1,nlevs
       call multifab_build(Source_nph(n), mla%la(n), 1, 0)
       call setval(Source_nph(n), ZERO, all=.true.)
    end do

    do n=1,nlevs
       call multifab_build(delta_p_term(n), mla%la(n), 1, 0)
       call setval(delta_p_term(n),ZERO,all=.true.)
    end do

    ! this should have no effect if dpdt_factor .le. 0
    p0_minus_peosbar = ZERO

    ! these should have no effect if evolve_base_state = F
    w0_force = ZERO
    Sbar = ZERO
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 3 -- construct the advective velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  3 : create MAC velocities>>> '
    end if

    do n=1,nlevs
       do comp=1,dm
          call multifab_build_edge(umac(n,comp), mla%la(n),1,1,comp)
       end do
    end do
    
    call advance_premac(uold,sold,umac,gpi,normal,w0,w0mac,w0_force,w0_force_cart, &
                        rho0_old,grav_cell_old,dx,dt,the_bc_tower%bc_tower_array,mla)

    if (dm .eq. 3) then
       do n=1,nlevs
          call destroy(w0_force_cart(n))
       end do
    end if

    do n=1,nlevs
       call multifab_build(delta_gamma1_term(n), mla%la(n), 1, 0)
       call multifab_build(macrhs(n),            mla%la(n), 1, 0)
       call multifab_build(delta_chi(n),         mla%la(n), 1, 0)

       call setval(delta_gamma1_term(n), ZERO, all=.true.)
       call setval(delta_chi(n),         ZERO, all=.true.)
    end do

    if (barrier_timers) call parallel_barrier()
    advect_time = advect_time + parallel_wtime() - advect_time_start

    macproj_time_start = parallel_wtime()

    call make_macrhs(macrhs,rho0_old,Source_nph,delta_gamma1_term,Sbar,div_coeff_old,dx, &
                     gamma1bar_old,p0_old,delta_p_term,dt,delta_chi,.true.)

    do n=1,nlevs
       call destroy(delta_gamma1_term(n))
       call destroy(Source_nph(n))
       call destroy(delta_p_term(n))
    end do

    do n=1,nlevs
       call multifab_build(macphi(n), mla%la(n), 1, 1)
       call setval(macphi(n), ZERO, all=.true.)
    end do

    ! MAC projection !
    call cell_to_edge(div_coeff_old,div_coeff_edge)
    call macproject(mla,umac,macphi,sold,dx,the_bc_tower, &
                    macrhs,div_coeff_1d=div_coeff_old,div_coeff_1d_edge=div_coeff_edge)

    do n=1,nlevs
       call destroy(macrhs(n))
    end do

    if (barrier_timers) call parallel_barrier()
    macproj_time = macproj_time + (parallel_wtime() - macproj_time_start)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 4 -- advect the base state and full state through dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    advect_time_start = parallel_wtime()

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  4 : advect base        '
    end if
    
    rho0_new = rho0_old

    do n=1,nlevs
       call multifab_build(thermal1(n), mla%la(n), 1, 0)
       call setval(thermal1(n),ZERO,all=.true.)
    end do
    
    do n=1,nlevs
       call multifab_build(s2(n), mla%la(n), nscal, nghost(sold(n)))
       ! copy temperature into s2 for seeding eos calls only
       ! temperature will be overwritten later after enthalpy advance
       call multifab_copy_c(s2(n), temp_comp, s1(n), temp_comp, 1, nghost(sold(n)))
    end do

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '            :  density_advance >>> '
       write(6,*) '            :   tracer_advance >>> '
    end if

    do n=1,nlevs
       do comp = 1,dm
          call multifab_build_edge(sedge(n,comp),mla%la(n),nscal,0,comp)
          call multifab_build_edge(sflux(n,comp),mla%la(n),nscal,0,comp)
       end do
       if (ppm_trace_forces == 0) then
          call multifab_build(scal_force(n), mla%la(n), nscal, 1)
       else
          ! we need more ghostcells if we are tracing the forces
          call multifab_build(scal_force(n), mla%la(n), nscal, nghost(sold(n)))
       endif
       call multifab_build_edge(etarhoflux(n), mla%la(n), 1, 0, dm)
       call setval(etarhoflux(n),ZERO,all=.true.)
    end do

    p0_new = p0_old
    rhoh0_new = rhoh0_old

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '            : enthalpy_advance >>> '
    end if

    do n = 1, nlevs
       do comp = 1,dm
          call destroy(sedge(n,comp))
          call destroy(sflux(n,comp))
          call destroy(umac(n,comp))
       end do
       call destroy(scal_force(n))
    end do

    if (barrier_timers) call parallel_barrier()
    advect_time = advect_time + parallel_wtime() - advect_time_start
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 4a (Option I) -- Add thermal conduction (only enthalpy terms)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    thermal_time_start = parallel_wtime()

    if (barrier_timers) call parallel_barrier()
    thermal_time = thermal_time + (parallel_wtime() - thermal_time_start)

    misc_time_start = parallel_wtime()

    ! pass temperature through for seeding the temperature update eos call
    ! pi goes along for the ride
    do n=1,nlevs
       call multifab_copy_c(s2(n),temp_comp,s1(n),temp_comp,1,nghost(sold(n)))
       call multifab_copy_c(s2(n),pi_comp,  s1(n),pi_comp,  1,nghost(sold(n)))
    end do
    ! now update temperature
    if (use_tfromp) then
       call makeTfromRhoP(s2,p0_new,mla,the_bc_tower%bc_tower_array,dx)
    else
       call makeTfromRhoH(s2,p0_new,mla,the_bc_tower%bc_tower_array,dx)
    end if

    if (barrier_timers) call parallel_barrier()
    misc_time = misc_time + parallel_wtime() - misc_time_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 5 -- react the full state and then base state through dt/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    react_time_start = parallel_wtime()

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  5 : react state     '
    end if

    do n=1,nlevs
       call multifab_build(rho_Hext(n), mla%la(n), 1, 0)
    end do

    do n=1,nlevs
       call destroy(s2(n))
    end do

    if (barrier_timers) call parallel_barrier()
    react_time = react_time + parallel_wtime() - react_time_start
    
    misc_time_start = parallel_wtime()

    ! Just copy div_coeff_new from div_coeff_old if not evolving the base state
    div_coeff_new = div_coeff_old
    div_coeff_nph = HALF*(div_coeff_old + div_coeff_new)

    if (barrier_timers) call parallel_barrier()
    misc_time = misc_time + parallel_wtime() - misc_time_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 6 -- define a new average expansion rate at n+1/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    advect_time_start = parallel_wtime()
    
    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  6 : make new S and new w0 >>> '
    end if

    ! reset cutoff coordinates to old time value
    call compute_cutoff_coords(rho0_old)

    do n=1,nlevs
       call setval(thermal2(n),ZERO,all=.true.)
    end do

    do n=1,nlevs
       call multifab_build(delta_gamma1_term(n), mla%la(n), 1, 0)
       call multifab_build(delta_gamma1(n), mla%la(n), 1, 0)
       call setval(delta_gamma1_term(n), ZERO, all=.true.)
       call setval(delta_gamma1(n), ZERO, all=.true.)
    end do

    do n=1,nlevs
       call destroy(rho_Hext(n))
       call destroy(delta_gamma1(n))
    end do

    do n=1,nlevs
       call multifab_build(Source_nph(n), mla%la(n), 1, 0)
    end do

    do n=1,nlevs
       call multifab_build(delta_p_term(n), mla%la(n), 1, 0)
       call setval(delta_p_term(n),ZERO,all=.true.)
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 7 -- redo the construction of the advective velocity using the 
!! current w0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  7 : create MAC velocities >>> '
    end if

    do n=1,nlevs
       do comp=1,dm
          call multifab_build_edge(umac(n,comp),mla%la(n),1,1,comp)
       end do
    end do

    if (barrier_timers) call parallel_barrier()
    advect_time = advect_time + parallel_wtime() - advect_time_start

    macproj_time_start = parallel_wtime()

    do n=1,nlevs
       call multifab_build(macrhs(n), mla%la(n), 1, 0)
    end do

    do n=1,nlevs
       call destroy(delta_gamma1_term(n))
       call destroy(Source_nph(n))
       call destroy(delta_p_term(n))
       call destroy(delta_chi(n))
    end do

    do n=1,nlevs
       call multifab_build(rhohalf(n), mla%la(n), 1, 1)
    end do

    do n=1,nlevs
       call destroy(rhohalf(n))
       call destroy(macrhs(n))
       call destroy(macphi(n))
    end do

    if (barrier_timers) call parallel_barrier()
    macproj_time = macproj_time + (parallel_wtime() - macproj_time_start)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 8 -- advect the base state and full state through dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    advect_time_start = parallel_wtime()

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  8 : advect base   '
    end if

    rho0_new = rho0_old

    do n=1,nlevs
       call multifab_build(s2(n), mla%la(n), nscal, nghost(sold(n)))
       ! copy temperature into s2 for seeding eos calls only
       ! temperature will be overwritten later after enthalpy advance
       call multifab_copy_c(s2(n), temp_comp, s1(n), temp_comp, 1, nghost(sold(n)))
       call setval(etarhoflux(n),ZERO,all=.true.)
    end do

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '            :  density_advance >>>'
       write(6,*) '            :   tracer_advance >>>'
    end if

    ! Build the sedge array.
    do n=1,nlevs
       do comp = 1,dm
          call multifab_build_edge(sedge(n,comp),mla%la(n),nscal,0,comp)
          call multifab_build_edge(sflux(n,comp),mla%la(n),nscal,0,comp)
       end do
       if (ppm_trace_forces == 0) then
          call multifab_build(scal_force(n), mla%la(n), nscal, 1)
       else
          call multifab_build(scal_force(n), mla%la(n), nscal, nghost(sold(n)))
       endif
    end do

    do n=1,nlevs
       call destroy(etarhoflux(n))
    end do

    grav_cell_new = grav_cell_old
    rho0_nph = rho0_old
    grav_cell_nph = grav_cell_old

    p0_new = p0_old
    rhoh0_new = rhoh0_old

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '            : enthalpy_advance >>>'
    end if

    do n=1,nlevs
       do comp = 1,dm
          call destroy(sedge(n,comp))
          call destroy(sflux(n,comp))
       end do
       call destroy(scal_force(n))
       call destroy(thermal1(n))
    end do

    if (barrier_timers) call parallel_barrier()
    advect_time = advect_time + parallel_wtime() - advect_time_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 8a (Option I) -- Add thermal conduction (only enthalpy terms)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    thermal_time_start = parallel_wtime()

    if (barrier_timers) call parallel_barrier()
    thermal_time = thermal_time + (parallel_wtime() - thermal_time_start)

    misc_time_start = parallel_wtime()

    ! pass temperature through for seeding the temperature update eos call
    ! pi just goes along for the ride too
    do n=1,nlevs
       call multifab_copy_c(s2(n),temp_comp,s1(n),temp_comp,1,nghost(sold(n)))
       call multifab_copy_c(s2(n),pi_comp,  s1(n),pi_comp,  1,nghost(sold(n)))
    end do

    ! now update temperature
    if (use_tfromp) then
       call makeTfromRhoP(s2,p0_new,mla,the_bc_tower%bc_tower_array,dx)
    else
       call makeTfromRhoH(s2,p0_new,mla,the_bc_tower%bc_tower_array,dx)
    end if

    do n=1,nlevs
       call destroy(s1(n))
    end do

    if (barrier_timers) call parallel_barrier()
    misc_time = misc_time + parallel_wtime() - misc_time_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 9 -- react the full state and then base state through dt/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    react_time_start = parallel_wtime()

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  9 : react state '
    end if

    do n=1,nlevs
       call multifab_build(rho_Hext(n), mla%la(n), 1, 0)
    end do

    do n=1,nlevs
       call destroy(s2(n))
    end do

    if (barrier_timers) call parallel_barrier()
    react_time = react_time + parallel_wtime() - react_time_start

    misc_time_start = parallel_wtime()

    ! compute gamma1bar
    div_coeff_nph = HALF*(div_coeff_old+div_coeff_new)

    if (barrier_timers) call parallel_barrier()
    misc_time = misc_time + parallel_wtime() - misc_time_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 10 -- compute S^{n+1} for the final projection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ndproj_time_start = parallel_wtime()
    
    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP 10 : make new S >>>'
    end if
          
    do n=1,nlevs
       call setval(thermal2(n),ZERO,all=.true.)
    end do
    
    do n=1,nlevs
       call multifab_build(delta_gamma1_term(n), mla%la(n), 1, 0)
       call multifab_build(delta_gamma1(n), mla%la(n), 1, 0)
    end do

    ! p0 is only used for the delta_gamma1_term

    do n=1,nlevs
       call destroy(delta_gamma1(n))
    end do

    ! define dSdt = (Source_new - Source_old) / dt
    do n=1,nlevs
       call multifab_copy(dSdt(n),Source_new(n))
       call multifab_sub_sub(dSdt(n),Source_old(n))
       call multifab_div_div_s(dSdt(n),dt)
       call setval(dSdt(n),ZERO,all=.true.)
    end do

    if (barrier_timers) call parallel_barrier()
    ndproj_time = ndproj_time + parallel_wtime() - ndproj_time_start
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 11 -- update the velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    advect_time_start = parallel_wtime()

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP 11 : update and project new velocity >>>'
    end if
    
    ! Define rho at half time using the new rho from Step 8
    do n=1,nlevs
       call multifab_build(rhohalf(n), mla%la(n), 1, 1)
    end do
    !call make_at_halftime(rhohalf,sold,snew,rho_comp,1,the_bc_tower%bc_tower_array,mla)
    call velocity_advance(mla,uold,unew,sold,rhohalf,umac,gpi,normal,w0,w0mac,w0_force, &
                          w0_force_cart,rho0_old,rho0_nph,grav_cell_old,grav_cell_nph, &
                          dx,dt,the_bc_tower%bc_tower_array,sponge)

    do n=1,nlevs
       do comp=1,dm
          call destroy(umac(n,comp))
       end do
    end do

    if (barrier_timers) call parallel_barrier()
    advect_time = advect_time + parallel_wtime() - advect_time_start
       
    ndproj_time_start = parallel_wtime()

    ! Project the new velocity field.
    proj_type = regular_timestep_comp

    do n=1,nlevs
       call destroy(delta_gamma1_term(n))
    end do

    do n=1,nlevs
       call multifab_build(div_coeff_3d(n), mla%la(n), 1, 1)
    end do
       
    do n = 1, nlevs
       call multifab_build(pi_cc(n), mla%la(n), 1, 0)
       call setval(pi_cc(n), ZERO, all=.true.)
    enddo

    do n = 1, nlevs
       call multifab_copy_c(snew(n),pi_comp,pi_cc(n),1,1)
       call destroy(pi_cc(n))
    enddo

    do n=1,nlevs
       call destroy(div_coeff_3d(n))
       call destroy(rhohalf(n))
    end do
    
    if (barrier_timers) call parallel_barrier()
    ndproj_time = ndproj_time + (parallel_wtime() - ndproj_time_start)

    misc_time_start = parallel_wtime()

    if (.not. init_mode) then
       write(*, *) 'NOT INIT MODE'
       grav_cell_old = grav_cell_new

       if (.not. fix_base_state) then
          ! compute tempbar by "averaging"
          call average(mla,snew,tempbar,dx,temp_comp)
       end if

       ! output any runtime diagnostics
       ! pass in the new time value, time+dt
       call diag(time+dt,dt,dx,snew,rho_Hnuc2,rho_Hext,thermal2,rho_omegadot2,&
                 rho0_new,rhoh0_new,p0_new,tempbar, &
                 gamma1bar,div_coeff_new, &
                 unew,w0,normal, &
                 mla,the_bc_tower)


       ! perform sanity checks, if desired
       if (mach_max_abort > ZERO) then
          call sanity_check(time+dt,dx,snew, &
                 rho0_new,rhoh0_new,p0_new,tempbar, &
                 gamma1bar,div_coeff_new, &
                 unew,w0,normal, &
                 mla,the_bc_tower)
       endif

    end if

    if (barrier_timers) call parallel_barrier()
    misc_time = misc_time + parallel_wtime() - misc_time_start

    call destroy(bpt)

    call parallel_reduce(advect_time_max, advect_time, MPI_MAX, &
                         proc=parallel_IOProcessorNode())

    call parallel_reduce(macproj_time_max,   macproj_time, MPI_MAX, &
                         proc=parallel_IOProcessorNode())

    call parallel_reduce(ndproj_time_max,   ndproj_time, MPI_MAX, &
                         proc=parallel_IOProcessorNode())

    call parallel_reduce(react_time_max,  react_time, MPI_MAX, &
                         proc=parallel_IOProcessorNode())

    call parallel_reduce(misc_time_max,  misc_time, MPI_MAX, &
                         proc=parallel_IOProcessorNode())

    if(use_thermal_diffusion) then
      call parallel_reduce(thermal_time_max,  thermal_time, MPI_MAX, &
                           proc=parallel_IOProcessorNode())
    end if

    if (parallel_IOProcessor()) then
       print *, 'Timing summary:'
       print *, '   Advection       : ', advect_time_max , ' seconds'
       print *, '   MAC   Projection: ', macproj_time_max, ' seconds'
       print *, '   Nodal Projection: ', ndproj_time_max , ' seconds'
       if (use_thermal_diffusion) &
          print *, '   Thermal         : ', thermal_time_max, ' seconds'
       print *, '   Reactions       : ', react_time_max  , ' seconds'
       print *, '   Misc            : ', misc_time_max   , ' seconds'
       print *, ' '
    endif

    do n = 1, nlevs
       call setval(snew(n), ONE, all=.true.)
       call setval(sold(n), ONE, all=.true.)
    end do

  end subroutine advance_timestep

end module advance_timestep_module
