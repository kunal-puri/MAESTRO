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
    use probin_module, only: dm_in, plot_processors

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
    plot_names(icomp_rhoh) = "energy"

    plot_names(icomp_magvel)      = "magvel"
    plot_names(icomp_mom)         = "momentum"
    plot_names(icomp_vort)        = "vort"
    plot_names(icomp_press)       = "P"

    if (plot_processors) then
       plot_names(icomp_proc) = "processor_number"
    endif

    plot_names(icomp_div)      = "divergence"

  end subroutine get_plot_names

  subroutine make_plotfile(dirname,mla,u,w,mba,plot_names,dx,dt,the_bc_tower,&
                           write_pf_time)
    use bl_prof_module
    use fabio_module
    use variables
    use plot_variables_module
    use probin_module

    use ml_restrict_fill_module
    use bl_constants_module
    use time_module, only: time

    character(len=*) , intent(in   ) :: dirname
    type(ml_layout)  , intent(in   ) :: mla
    type(multifab)   , intent(in   ) :: u(:),w(:)
    type(ml_boxarray), intent(in   ) :: mba
    character(len=20), intent(in   ) :: plot_names(:)
    real(dp_t)       , intent(in   ) :: dt,dx(:,:)
    type(bc_tower)   , intent(in   ) :: the_bc_tower
    real(dp_t)       , intent(  out) :: write_pf_time

    type(multifab) :: plotdata(mla%nlevel)
    !type(multifab) ::  tempfab(mla%nlevel)

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
       !call multifab_build(tempfab(n),  mla%la(n), dm,           1)

       ! DENSITY
       write(*, *) 'DENSITY', icomp_rho, rho_comp
       call multifab_copy_c(plotdata(n),icomp_rho,u(n),rho_comp,1)
              
       ! VELOCITY 
       write(*, *) 'VELOCITY', icomp_vel, 2
       call multifab_copy_c(plotdata(n),icomp_vel,u(n),2,dm)

       ! ENERGY
       write(*, *) 'Energy', icomp_rhoh, rhoh_comp
       call multifab_copy_c(plotdata(n), icomp_rhoh,u(n),rhoh_comp,1)

    end do

    do n = 1,nlevs

       ! MAGVEL = |U + w0| AND MOMENTUM
       call make_magvel(plotdata(n),icomp_magvel,icomp_mom,w(n))

       ! VORTICITY
       call make_vorticity(plotdata(n),icomp_vort,w(n),dx(n,:), &
                          the_bc_tower%bc_tower_array(n))

       ! PRESSURE
       call multifab_copy_c(plotdata(n), icomp_press, w(n), press_comp,1)

       ! PROCESSOR NUMBER
       if (plot_processors) then
          call make_processor_number(plotdata(n),icomp_proc)
       end if

       ! DIVERGENCE
       call make_divergence(plotdata(n), icomp_div, w(n), dx(n,:), &
                            the_bc_tower%bc_tower_array(n))
    end do

    ! restrict data and fill all ghost cells
    ! call ml_restrict_and_fill(nlevs,tempfab,mla%mba%rr,the_bc_tower%bc_tower_array, &
    !                           icomp=1, &
    !                           bcomp=foextrap_comp, &
    !                           nc=1, &
    !                           ng=tempfab(1)%ng)
    
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
       !call destroy(tempfab(n))
    end do

    call destroy(bpt)

  end subroutine make_plotfile

end module make_plotfile_module
