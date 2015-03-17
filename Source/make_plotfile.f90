module make_plotfile_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_boxarray_module

  implicit none

  private
  public :: get_plot_names, make_plotfile

contains

  subroutine get_plot_names(plot_names)

    use plot_variables_module
    use variables
    use network, only: nspec, short_spec_names
    use probin_module, only: plot_spec, plot_trac, plot_base, &
                             use_thermal_diffusion, plot_omegadot, plot_Hnuc, &
                             plot_Hext, plot_eta, plot_ad_excess, &
                             use_tfromp, plot_h_with_use_tfromp, plot_gpi, plot_cs, &
                             plot_sponge_fdamp, dm_in, use_particles, &
                             plot_processors, plot_pidivu
    use geometry, only: spherical

    character(len=20), intent(inout) :: plot_names(:)

    ! Local variables
    integer :: comp

    plot_names(icomp_vel  ) = "x_vel"
    if (dm_in > 1) then
       plot_names(icomp_vel+1) = "y_vel"
    end if
    if (dm_in > 2) then
       plot_names(icomp_vel+2) = "z_vel"
    end if
    plot_names(icomp_rho)  = "density"

    plot_names(icomp_magvel)      = "magvel"
    plot_names(icomp_mom)         = "momentum"
    plot_names(icomp_vort)        = "vort"
    plot_names(icomp_src)         = "S"

    if (plot_processors) then
       plot_names(icomp_proc) = "processor_number"
    endif

  end subroutine get_plot_names

  subroutine make_plotfile(dirname,mla,u,s,pi,gpi,rho_omegadot, &
                           rho_Hnuc,rho_Hext, &
                           thermal,Source,sponge,mba,plot_names,dx, &
                           the_bc_tower,w0,rho0,rhoh0,p0, &
                           tempbar,gamma1bar,etarho_cc, &
                           normal,dt,particles,write_pf_time)

    use bl_prof_module
    use fabio_module
    use variables
    use plot_variables_module
    use fill_3d_module
    use probin_module, only: nOutFiles, lUsingNFiles, plot_spec, plot_trac, & 
                             plot_base, plot_omegadot, plot_Hnuc, plot_Hext, &
                             plot_eta, plot_ad_excess, &
                             single_prec_plotfiles, &
                             do_smallscale, use_thermal_diffusion, &
                             evolve_base_state, prob_lo, prob_hi, &
                             use_tfromp, plot_h_with_use_tfromp, plot_gpi, &
                             plot_cs, sponge_kappa, plot_sponge_fdamp, use_particles, &
                             plot_processors, plot_pidivu, use_alt_energy_fix
    use geometry, only: spherical, nr_fine, nlevs_radial, numdisjointchunks, &
         r_start_coord, r_end_coord
    use average_module
    use ml_restrict_fill_module
    use bl_constants_module
    use network, only: nspec
    use time_module, only: time
    use particle_module, only: particle_container, make_particle_count
    use make_grav_module
    use make_div_coeff_module
    use make_pi_cc_module

    character(len=*) , intent(in   ) :: dirname
    type(ml_layout)  , intent(in   ) :: mla
    type(multifab)   , intent(in   ) :: u(:)
    type(multifab)   , intent(in   ) :: s(:)
    type(multifab)   , intent(in   ) :: pi(:)
    type(multifab)   , intent(in   ) :: gpi(:)
    type(multifab)   , intent(in   ) :: rho_omegadot(:)
    type(multifab)   , intent(in   ) :: rho_Hnuc(:)
    type(multifab)   , intent(in   ) :: rho_Hext(:)
    type(multifab)   , intent(in   ) :: thermal(:)
    type(multifab)   , intent(in   ) :: Source(:)
    type(multifab)   , intent(in   ) :: sponge(:)
    type(ml_boxarray), intent(in   ) :: mba
    character(len=20), intent(in   ) :: plot_names(:)
    real(dp_t)       , intent(in   ) :: dt,dx(:,:)
    type(bc_tower)   , intent(in   ) :: the_bc_tower
    real(dp_t)       , intent(in   ) :: w0(:,0:)
    real(dp_t)       , intent(in   ) :: rho0(:,0:)
    real(dp_t)       , intent(in   ) :: rhoh0(:,0:)
    real(dp_t)       , intent(in   ) :: p0(:,0:)
    real(dp_t)       , intent(in   ) :: tempbar(:,0:)
    real(dp_t)       , intent(in   ) :: gamma1bar(:,0:)
    real(dp_t)       , intent(in   ) :: etarho_cc(:,0:)
    type(multifab)   , intent(in   ) :: normal(:)
    type(particle_container), intent(inout) :: particles
    real(dp_t)       , intent(  out) :: write_pf_time

    type(multifab) :: plotdata(mla%nlevel)
    type(multifab) ::  tempfab(mla%nlevel)
    type(multifab) ::    w0mac(mla%nlevel,mla%dim)
    type(multifab) :: w0r_cart(mla%nlevel)
    type(multifab) ::    pi_cc(mla%nlevel)

    real(dp_t)  :: div_coeff(nlevs_radial,0:nr_fine-1)
    real(dp_t)  :: grav_cell(nlevs_radial,0:nr_fine-1)

    real(dp_t) :: entropybar(nlevs_radial,0:nr_fine-1)
    real(dp_t) ::         h0(nlevs_radial,0:nr_fine-1)

    real(dp_t) :: tempval

    integer :: n,r,j,n_1d,prec,comp,dm,nlevs

    type(bl_prof_timer), save :: bpt

    real(dp_t) :: writetime1, writetime2

    call build(bpt, "make_plotfile")

    dm = mla%dim
    nlevs = mla%nlevel

    if (single_prec_plotfiles) then
       prec = FABIO_SINGLE
    else
       prec = FABIO_DOUBLE
    endif

    do n = 1,nlevs

       call multifab_build(plotdata(n), mla%la(n), n_plot_comps, 0)
       call multifab_build(tempfab(n),  mla%la(n), dm,           1)
              
       ! VELOCITY 
       call multifab_copy_c(plotdata(n),icomp_vel,u(n),1,dm)

       ! DENSITY
       call multifab_copy_c(plotdata(n),icomp_rho,s(n),rho_comp,1)

    end do

    do n = 1,nlevs

       n_1d = n

       ! MAGVEL = |U + w0|
       call make_magvel(plotdata(n),icomp_magvel,icomp_mom,s(n),u(n),w0(n_1d,:),w0mac(n,:))

       ! VORTICITY
       call make_vorticity(plotdata(n),icomp_vort,u(n),dx(n,:), &
                           the_bc_tower%bc_tower_array(n))

       ! DIVU
       call multifab_copy_c(plotdata(n),icomp_src,Source(n),1,1)

    end do

    ! processor number
    if (plot_processors) then
       do n = 1, nlevs
          call make_processor_number(plotdata(n),icomp_proc)
       enddo
    endif
    
    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,tempfab,mla%mba%rr,the_bc_tower%bc_tower_array, &
                              icomp=1, &
                              bcomp=foextrap_comp, &
                              nc=1, &
                              ng=tempfab(1)%ng)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (parallel_IOProcessor()) then
       write(6,*) 'Writing state to plotfile ',trim(dirname)
    end if

    writetime1 = parallel_wtime()

    call fabio_ml_multifab_write_d(plotdata, mba%rr(:,1), dirname, plot_names, &
                                   mba%pd(1), prob_lo, prob_hi, time, dx(1,:), &
                                   nOutFiles = nOutFiles, &
                                   lUsingNFiles = lUsingNFiles, prec = prec)

    writetime2 = parallel_wtime() - writetime1
    call parallel_reduce(writetime1, writetime2, MPI_MAX, proc=parallel_IOProcessorNode())
    if (parallel_IOProcessor()) then
       print*,'Time to write plotfile: ',writetime1,' seconds'
       print*,''
    end if

    write_pf_time = writetime1

    do n = 1,nlevs
       call destroy(plotdata(n))
       call destroy(tempfab(n))
    end do

    call destroy(bpt)

  end subroutine make_plotfile

end module make_plotfile_module
