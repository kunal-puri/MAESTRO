module initial_proj_module

  implicit none

  private
  public :: initial_proj

contains

  subroutine initial_proj(nlevs,uold,sold,pres,gpres,vel_force,normal,rho_Hext,Source_old, &
                          hgrhs,rho_omegadot1,div_coeff_3d,rhohalf, &
                          div_coeff_old,s0_old,p0_old,gam1,grav_cell,dx,the_bc_tower,mla)

    use variables, only: temp_comp, rho_comp, press_comp
    use define_bc_module
    use bl_constants_module
    use probin_module
    use geometry, only: spherical, nr
    use proj_parameters, only: initial_projection_comp
    use mk_vel_force_module
    use make_explicit_thermal_module
    use make_S_module
    use average_module
    use hgrhs_module
    use fill_3d_module
    use hgproject_module
    use multifab_module
    use ml_layout_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: uold(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: pres(:)
    type(multifab) , intent(inout) :: gpres(:)
    type(multifab) , intent(inout) :: vel_force(:)
    type(multifab) , intent(in   ) :: normal(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    type(multifab) , intent(inout) :: Source_old(:)
    type(multifab) , intent(inout) :: hgrhs(:)
    type(multifab) , intent(inout) :: rho_omegadot1(:)
    type(multifab) , intent(inout) :: div_coeff_3d(:)
    type(multifab) , intent(inout) :: rhohalf(:)
    real(kind=dp_t), intent(in   ) :: div_coeff_old(:,0:)
    real(kind=dp_t), intent(in   ) :: s0_old(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: gam1(:,0:)
    real(kind=dp_t), intent(in   ) :: grav_cell(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(ml_layout), intent(inout) :: mla

    ! local
    integer                     :: n
    real(dp_t)                  :: dt_temp
    real(dp_t), allocatable     :: Sbar(:,:,:)  
    type(multifab), allocatable :: gamma1_term(:)
    type(multifab), allocatable :: thermal(:)

    if ( parallel_IOProcessor() ) then
       print *, 'DOING THE INITIAL VELOCITY PROJECTION'
       print *, ' '
    end if
    
    call mk_vel_force(nlevs,vel_force,gpres,sold,normal,s0_old,grav_cell,dx, &
                      the_bc_tower%bc_tower_array,mla)
    
    do n = 1, nlevs
       ! we don't have a legit timestep yet, so set rho_omegadot1 = 0 to be safe
       call setval(rho_omegadot1(n), ZERO, all=.true.)
    end do

    allocate(thermal(nlevs))

    do n=1,nlevs
       call multifab_build(thermal(n), mla%la(n), 1, 1)
       call setval(thermal(n), 0.0_dp_t, all=.true.)
    end do
    
    if(use_thermal_diffusion) then
       call make_explicit_thermal(mla,dx,thermal,sold,p0_old,mg_verbose,cg_verbose, &
                                  the_bc_tower,temp_diffusion_formulation)
    end if
    
    allocate(gamma1_term(nlevs))

    do n=1,nlevs
       call multifab_build(gamma1_term(n), mla%la(n), 1, 0)
       call setval(gamma1_term(n), 0.0_dp_t, all=.true.)
    end do

    call make_S(nlevs,Source_old,gamma1_term,sold,rho_omegadot1,rho_Hext,thermal, &
                s0_old(:,:,temp_comp),gam1,dx)

    do n=1,nlevs
       call destroy(thermal(n))
    end do

    deallocate(thermal)
    
    allocate(Sbar(nlevs,nr(nlevs),1))

    call average(mla,Source_old,Sbar,dx,1,1)
    
    ! Note that we use rhohalf, filled with 1 at this point, as a temporary
    ! in order to do a constant-density initial projection.
    do n=1,nlevs
       call setval(rhohalf(n),ONE,rho_comp,1,all=.true.)
       call setval(hgrhs(n),ZERO,all=.true.)
    end do
    
    call make_hgrhs(nlevs,hgrhs,Source_old,gamma1_term,Sbar(:,:,1),div_coeff_old,dx)

    do n=1,nlevs
       call destroy(gamma1_term(n))
    end do

    deallocate(gamma1_term,Sbar)
    
    ! dt doesn't matter for the initial projection since we're throwing
    ! away the p and gpres anyway
    dt_temp = ONE

    if (spherical .eq. 1) then
       call fill_3d_data_wrapper(nlevs,div_coeff_3d,div_coeff_old,dx)
       call hgproject(initial_projection_comp,mla,uold,uold,rhohalf,pres,gpres,dx, &
                      dt_temp,the_bc_tower,verbose,mg_verbose,cg_verbose,press_comp, &
                      hgrhs,div_coeff_3d=div_coeff_3d,eps_in=1.d-10)
       
    else
       call hgproject(initial_projection_comp,mla,uold,uold,rhohalf,pres,gpres,dx, &
                      dt_temp,the_bc_tower,verbose,mg_verbose,cg_verbose,press_comp, &
                      hgrhs,div_coeff_1d=div_coeff_old)
    end if
    
    do n = 1,nlevs
       call setval( pres(n)  ,0.0_dp_t, all=.true.)
       call setval(gpres(n)  ,0.0_dp_t, all=.true.)
    end do
    
  end subroutine initial_proj

end module initial_proj_module
