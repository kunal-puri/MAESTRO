subroutine varden()
  use variables
  use geometry
  !use estdt_module
  !use firstdt_module
  !use advance_timestep_module
  use box_util_module
  use fabio_module
  !use checkpoint_module
  use make_plotfile_module
  use probin_module
  use runtime_init_module
  use bl_constants_module
  use define_bc_module
  use initialize_module
  use make_new_grids_module
  !use regrid_module
  use omp_module
  use time_module, only: time
  use cputime_module, only: start_cputime_clock
  use simple_log_module, only: simple_log_init, simple_log_finalize
  
  implicit none

  integer    :: init_step,istep
  integer    :: istep_divu_iter,istep_init_iter
  integer    :: i,n,r,numcell,dm,nlevs
  integer    :: last_plt_written,last_chk_written
  real(dp_t) :: smin,smax,smaxold,nuclear_dt_scalefac
  real(dp_t) :: umin,umax,vmin,vmax,wmin,wmax
  real(dp_t) :: dt,dtold

  type(ml_layout) :: mla
  type(bc_tower)  :: the_bc_tower

  real(dp_t)  , pointer     :: dx(:,:)

  type(multifab), allocatable :: unew(:)
  type(multifab), allocatable :: wnew(:)
  
  ! these are pointers because they need to be allocated and built within 
  ! another function
  type(multifab), pointer :: uold(:)
  type(multifab), pointer :: wold(:)

  type(multifab), pointer :: chkdata(:)

  character(len=5)               :: plot_index, check_index
  character(len=6)               :: plot_index6, check_index6
  character(len=7)               :: plot_index7, check_index7
  character(len=256)             :: plot_file_name, check_file_name
  character(len=20), allocatable :: plot_names(:)
  
  integer :: npartdata
  integer, allocatable :: index_partdata(:)
  character(len=16), allocatable :: names_partdata(:)

  real(dp_t), parameter :: SMALL = 1.d-13
  real(dp_t)            :: runtime1, runtime2

  logical :: init_mode

  logical :: dump_plotfile, dump_checkpoint, abort_maestro
  real(dp_t) :: write_pf_time

  integer :: numparticles

  logical :: have_overview

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! keep track of cputime
  call start_cputime_clock()

  last_plt_written = -1
  last_chk_written = -1

  call runtime_init()

  call init_variables()
  call init_plot_variables()

  allocate(plot_names(n_plot_comps))
  call get_plot_names(plot_names)

  ! a simple logging facility for writting both to the screen and a file
  call simple_log_init()

  !if (restart >= 0) then

  !    call initialize_from_restart(mla,restart,dt,pmask,dx,uold,sold,gpi,pi, &
  !                                 dSdt,Source_old,Source_new, &
  !                                 rho_omegadot2,rho_Hnuc2,rho_Hext,thermal2,the_bc_tower, &
  !                                 div_coeff_old,div_coeff_new,gamma1bar,gamma1bar_hold, &
  !                                 s0_init,rho0_old,rhoh0_old,rho0_new,rhoh0_new,p0_init, &
  !                                 p0_old,p0_new,w0,etarho_ec,etarho_cc,psi, &
  !                                 tempbar,tempbar_init,grav_cell)

  !    if (restart <= 99999) then
  !       write(unit=check_index,fmt='(i5.5)') restart
  !       check_file_name = trim(check_base_name) // check_index
  !    else if (restart <= 999999) then
  !       write(unit=check_index6,fmt='(i6.6)') restart
  !       check_file_name = trim(check_base_name) // check_index6
  !    else
  !       write(unit=check_index7,fmt='(i7.7)') restart
  !       check_file_name = trim(check_base_name) // check_index7
  !    endif

  !    if (use_particles) then
  !       call particle_container_restart(particles,check_file_name,&
  !                                       mla,dx,prob_lo)
  !    endif
     
  ! else if (test_set /= '') then
  !    call initialize_with_fixed_grids(mla,dt,pmask,dx,uold,sold,gpi,pi,dSdt, &
  !                                     Source_old,Source_new, &
  !                                     rho_omegadot2,rho_Hnuc2,rho_Hext,thermal2, &
  !                                     the_bc_tower, &
  !                                     div_coeff_old,div_coeff_new,gamma1bar, &
  !                                     gamma1bar_hold,s0_init,rho0_old,rhoh0_old, &
  !                                     rho0_new,rhoh0_new,p0_init,p0_old,p0_new,w0, &
  !                                     etarho_ec,etarho_cc,psi,tempbar,tempbar_init,grav_cell)

  !else
  write(*, *) 'Calling INITIALIZE'
  call initialize_with_adaptive_grids(mla,dt,pmask,dx,uold,wold,the_bc_tower)
  !end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! error checking
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dm = dm_in
  nlevs = mla%nlevel

  ! check to make sure dimensionality is consistent in the inputs file
  if (dm .ne. get_dim(mla%mba)) then 
     call bl_error('dm_in not properly set in inputs file')
  end if

  ! check to make sure our grid is square -- the solvers assume this
  if (dm == 2) then
     if (abs(dx(1,1) - dx(1,2)) > SMALL) then
        call bl_error('zones must be square')
     end if
  else if (dm == 3) then
     if (abs(dx(1,1) - dx(1,2)) > SMALL .OR. abs(dx(1,1) - dx(1,3)) > SMALL) then
        call bl_error('zones must be square')
     end if
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! print processor and grid info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (parallel_IOProcessor()) then
     print *, 'number of MPI processes = ', parallel_nprocs()
     print *, 'number of threads       = ', omp_get_max_threads()
     print *, ' '
     print *, 'number of dimensions    = ', dm
     do n = 1, nlevs
        print *, 'level: ', n
        print *, '   number of boxes = ', nboxes(uold(n)%la)
        print *, '   maximum zones   = ', (extent(mla%mba%pd(n),i),i=1,dm)
     end do
     print *, ' '
  end if

  if (verbose .ge. 1) then
     do n=1,nlevs
        numcell = multifab_volume(uold(n),.false.)
        if (parallel_IOProcessor()) then
           print*,"Number of valid cells at level        ",n,numcell
        end if
        numcell = multifab_volume(uold(n),.true.)
        if (parallel_IOProcessor()) then
           print*,"Number of valid + ghost cells at level",n,numcell
        end if
     end do
     if (parallel_IOProcessor()) then
        print*,""
     end if
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize all remaining arrays and multifabs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(*, *) 'Allocating and BUILDING UNEW/WNEW'
  allocate(unew(nlevs),wnew(nlevs))

  do n = 1,nlevs
     call multifab_build( unew(n), mla%la(n), n_comp, nghost(uold(n)) )
     call multifab_build( wnew(n), mla%la(n), n_comp, nghost(wold(n)) )

     call setval(      unew(n), ZERO, all=.true.)
     call setval(      wnew(n), ZERO, all=.true.)
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Print out some diagnostics and warnings about cutoff / sponge parameter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

1010 format(1x, "restarted from step ", i7, 2x, "MPI tasks = ", i6, 2x, "OMP threads = ", i3)

  ! write the maestro-overview.out file
  if (restart < 0) then
     !call write_job_info("", mla%mba, the_bc_tower, 0.0d0)
  else
     if (parallel_IOProcessor()) then
        inquire(file="maestro-overview.out", exist=have_overview)
        if (have_overview) then
           open(unit=99, file="maestro-overview.out", form = "formatted", &
                action="write", status="old", position="append")        
           write(99,1010) restart, parallel_nprocs(), omp_get_max_threads()
           close(unit=99)
        endif
     endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (restart < 0) then 
     istep = 0
  else 
     istep = restart
  end if

  if (stop_time >= 0.d0) then
     if (time+dt > stop_time) dt = min(dt, stop_time - time)
  end if

  write(*, *) 'RESTART MODE?'
  dtold = dt
  new_sim_mode: if (restart < 0) then
     !if ( chk_int > 0 ) then

        !-----------------------------------------------------------------------
        ! write a checkpoint file
        !-----------------------------------------------------------------------

     !    allocate(chkdata(nlevs))
     !    do n = 1,nlevs
     !       call multifab_build(chkdata(n), mla%la(n), 2*dm+nscal, 0)
     !       call multifab_copy_c(chkdata(n),1                ,uold(n), 1,dm)
     !       call multifab_copy_c(chkdata(n),rho_comp+dm      ,sold(n), 1,nscal)
     !    end do

     !    if (istep <= 99999) then
     !       write(unit=check_index,fmt='(i5.5)') istep
     !       check_file_name = trim(check_base_name) // check_index
     !    else if (istep <= 999999) then
     !       write(unit=check_index6,fmt='(i6.6)') istep
     !       check_file_name = trim(check_base_name) // check_index6
     !    else
     !       write(unit=check_index7,fmt='(i7.7)') istep
     !       check_file_name = trim(check_base_name) // check_index7
     !    endif

     !    call checkpoint_write(check_file_name, chkdata, &
     !                          pi, dSdt, Source_old, Source_new, &
     !                          rho_omegadot2, rho_Hnuc2, rho_Hext, thermal2, &
     !                          mla%mba%rr, dt)

     !    call write_base_state(istep, check_file_name, &
     !                          rho0_old, rhoh0_old, p0_old, gamma1bar, &
     !                          w0, etarho_ec, etarho_cc, &
     !                          div_coeff_old, psi, tempbar, tempbar_init, prob_lo(dm))

     !    call write_aux_data(istep, check_file_name)

     !    last_chk_written = istep

     !    do n = 1,nlevs
     !       call destroy(chkdata(n))
     !    end do
     !    deallocate(chkdata)

     !    if (use_particles) then
     !       call particle_container_checkpoint(particles,check_file_name,mla)
     !    endif
     ! end if

     if ( plot_int > 0 .or. plot_deltat > ZERO) then

        !-----------------------------------------------------------------------
        ! write a plotfile
        !-----------------------------------------------------------------------

        if (istep <= 99999) then
           write(unit=plot_index,fmt='(i5.5)') istep
           plot_file_name = trim(plot_base_name) // plot_index
        else if (istep <= 999999) then
           write(unit=plot_index6,fmt='(i6.6)') istep
           plot_file_name = trim(plot_base_name) // plot_index6
        else
           write(unit=plot_index7,fmt='(i7.7)') istep
           plot_file_name = trim(plot_base_name) // plot_index7
        endif

        write(*, *) 'WRITING PLOTFILE'
        call make_plotfile(plot_file_name,mla,uold,wold,mla%mba,plot_names,dx,dt, &
                           the_bc_tower,write_pf_time)

        !call write_job_info(plot_file_name, mla%mba, the_bc_tower, write_pf_time)
        last_plt_written = istep
     end if

  end if new_sim_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main evolution loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (restart < 0) then
     init_step = 1
  else
     init_step = restart+1
  end if

  if ( parallel_IOProcessor()) then
     print*,""
     print*,"BEGIN MAIN EVOLUTION LOOP WITH dt =",dt
     print*,""
  end if

  if ( (max_step >= init_step) .and. (time < stop_time .or. stop_time < 0.d0) ) then
     do istep = init_step,max_step
     
        if ( verbose .ge. 1 ) then
           if ( parallel_IOProcessor() ) then
              print *, 'MEMORY STATS AT START OF TIMESTEP ', istep
              print*, ' '
           end if
           call print(multifab_mem_stats(),    "    multifab")
           call print(fab_mem_stats(),         "         fab")
           call print(boxarray_mem_stats(),    "    boxarray")
           call print(layout_mem_stats(),      "      layout")
           call print(boxassoc_mem_stats(),    "    boxassoc")
           call print(fgassoc_mem_stats(),     "     fgassoc")
           call print(syncassoc_mem_stats(),   "   syncassoc")
           call print(copyassoc_mem_stats(),   "   copyassoc")
           call print(fluxassoc_mem_stats(),   "   fluxassoc")
           if ( parallel_IOProcessor() ) print*, ''

           do n = 1,nlevs
              umax = multifab_max_c(uold(n),1)
              umin = multifab_min_c(uold(n),1)
              if (dm.eq.1) then
                 if ( parallel_IOProcessor()) then
                    write(6,1001) n,umax
                    write(6,1011) n,umin
                 end if
              else if (dm .eq. 2) then
                 vmax = multifab_max_c(uold(n),2)
                 vmin = multifab_min_c(uold(n),2)
                 if ( parallel_IOProcessor()) then
                    write(6,1002) n,umax,vmax
                    write(6,1012) n,umin,vmin
                 end if
              
              else if (dm .eq. 3) then
                 vmax = multifab_max_c(uold(n),2)
                 vmin = multifab_min_c(uold(n),2)
                 wmax = multifab_max_c(uold(n),3)
                 wmin = multifab_min_c(uold(n),3)
                 if ( parallel_IOProcessor()) then
                    write(6,1003) n,umax,vmax,wmax
                    write(6,1013) n,umin,vmin,wmin
                 end if
              end if
           end do
        end if

1001 format('... level / max : old vels   ',i2,2x,e17.10)
1011 format('... level / min : old vels   ',i2,2x,e17.10)
1002 format('... level / max : old vels   ',i2,2x,e17.10,2x,e17.10)
1012 format('... level / min : old vels   ',i2,2x,e17.10,2x,e17.10)
1003 format('... level / max : old vels   ',i2,2x,e17.10,2x,e17.10,2x,e17.10)
1013 format('... level / min : old vels   ',i2,2x,e17.10,2x,e17.10,2x,e17.10)

        ! !---------------------------------------------------------------------
        ! ! BEGIN REGRIDDING
        ! !---------------------------------------------------------------------
        ! regridding: if (max_levs > 1 .and. regrid_int > 0 .and. (mod(istep-1,regrid_int) .eq. 0)) then

        !    ! Regrid psi, etarho_cc, etarho_ec, and w0.
        !    ! We do not regrid these in spherical since the base state array only
        !    ! contains 1 level of refinement.
        !    ! We do not regrid these if evolve_base_state=F since they are
        !    ! identically zero, set this way in initialize.
        !    regridding_base_state: if (spherical .eq. 0) then
        !       ! evolve_base_state == F and spherical == 0
              
        !       ! Here we want to fill in the rho0 array so there is
        !       ! valid data in any new grid locations that are created
        !       ! during the regrid.
              
        !       ! copy the coarsest level of the real arrays into the
        !       ! temp arrays
        !       rho0_temp(1,:) = rho0_old(1,:)
              
        !       ! piecewise linear interpolation to fill the cc temp arrays
        !       do n=2,max_levs
        !          do r=0,nr(n)-1
        !             if (r .eq. 0 .or. r .eq. nr(n)-1) then
        !                rho0_temp(n,r) = rho0_temp(n-1,r/2)
        !             else
        !                if (mod(r,2) .eq. 0) then
        !                   rho0_temp(n,r) = 0.75d0*rho0_temp(n-1,r/2) &
        !                        + 0.25d0*rho0_temp(n-1,r/2-1)
        !                else
        !                   rho0_temp(n,r) = 0.75d0*rho0_temp(n-1,r/2) &
        !                        + 0.25d0*rho0_temp(n-1,r/2+1)
        !                end if
        !             end if
        !          end do
        !       end do
              
        !       ! copy valid data into temp -- this way we haven't
        !       ! overwritten any of the original information.
        !       do n=2,nlevs_radial
        !          do i=1,numdisjointchunks(n)
        !             do r=r_start_coord(n,i),r_end_coord(n,i)
        !                rho0_temp(n,r) = rho0_old(n,r)
        !             end do
        !          end do
        !       end do
              
        !       ! don't copy rho0_temp back into the rho0 array just
        !       ! yet.  Wait until we regrid
              
        !       ! regardless of evolve_base_state, if new grids were
        !       ! created, we need to initialize tempbar_init there, in
        !       ! case drive_initial_convection = T

        !       ! copy the coarsest level of the real arrays into the 
        !       ! temp arrays
        !       tempbar_init_temp(1,:)  = tempbar_init(1,:)

        !       ! piecewise linear interpolation to fill the cc temp arrays
        !       do n=2,max_levs
        !          do r=0,nr(n)-1
        !             if (r .eq. 0 .or. r .eq. nr(n)-1) then
        !                tempbar_init_temp(n,r) = tempbar_init_temp(n-1,r/2)
        !             else
        !                if (mod(r,2) .eq. 0) then
        !                   tempbar_init_temp(n,r) = 0.75d0*tempbar_init_temp(n-1,r/2) &
        !                        + 0.25d0*tempbar_init_temp(n-1,r/2-1)
        !                else
        !                   tempbar_init_temp(n,r) = 0.75d0*tempbar_init_temp(n-1,r/2) &
        !                        + 0.25d0*tempbar_init_temp(n-1,r/2+1)
        !                end if
        !             end if
        !          end do
        !       end do

        !       ! copy valid data into temp
        !       do n=2,nlevs_radial
        !          do i=1,numdisjointchunks(n)
        !             do r=r_start_coord(n,i),r_end_coord(n,i)
        !                tempbar_init_temp(n,r) = tempbar_init(n,r)
        !             end do
        !          end do
        !       end do

        !       ! copy temp array back into the real thing
        !       tempbar_init = tempbar_init_temp

        !    end if regridding_base_state ! end regridding of base state
           
        !    ! we can pass as an auxillary tagging quantity either the
        !    ! energy generate rate (per volume) or tpert
        !    if (use_tpert_in_tagging) then
        !       do n = 1, nlevs
        !          call multifab_copy_c(tag_mf(n), 1, sold(n), temp_comp, 1)
        !       enddo
              
        !       call put_in_pert_form(mla,tag_mf,tempbar,dx,1, &
        !                             foextrap_comp,.true., &
        !                             the_bc_tower%bc_tower_array)

        !    else
        !       ! figure out if we are tagging off of heating or rxns
        !       if (do_heating) then
        !          do n = 1, nlevs
        !             call multifab_copy_c(tag_mf(n), 1, rho_Hext(n), 1, 1)
        !          enddo
        !       else
        !          do n = 1, nlevs
        !             call multifab_copy_c(tag_mf(n), 1, rho_Hnuc2(n), 1, 1)
        !          enddo
        !       endif
        !    endif

        !    do n=1,nlevs
        !       call multifab_destroy(unew(n))
        !       call multifab_destroy(snew(n))
        !       call multifab_destroy(sponge(n))
        !       call multifab_destroy(hgrhs(n))
        !       call multifab_destroy(Source_new(n))
        !       call multifab_destroy(rho_omegadot2(n))
        !       call multifab_destroy(rho_Hnuc2(n))
        !       call multifab_destroy(rho_Hext(n))
        !       call multifab_destroy(thermal2(n))
        !       if (dm .eq. 3) then
        !          call multifab_destroy(normal(n))
        !       end if
        !    end do

        !    ! create new grids and fill in data on those grids
        !    call regrid(istep,mla,uold,sold,gpi,pi,dSdt,Source_old,dx,the_bc_tower, &
        !                rho0_old,rhoh0_old,.false.,tag_mf)

        !    ! nlevs is local so we need to reset it
        !    nlevs = mla%nlevel

        !    if (nlevs .ne. max_levs) then
        !       call bl_error('varden.f90: nlevs .ne. max_levs not supported yet')
        !    end if

        !    call init_multilevel(sold)

        !    do n = 1,nlevs
        !       call multifab_build(      unew(n),    mla%la(n),    dm, nghost(uold(n)))
        !       call multifab_build(      snew(n),    mla%la(n), nscal, nghost(sold(n)))
        !       call multifab_build(    sponge(n),    mla%la(n),     1, 0)
        !       call multifab_build(     hgrhs(n),    mla%la(n),     1, 0, nodal)
        !       call multifab_build(Source_new(n),    mla%la(n),     1, 1)
        !       call multifab_build(rho_omegadot2(n), mla%la(n), nspec, 0)
        !       call multifab_build(    rho_Hnuc2(n), mla%la(n),     1, 0)
        !       call multifab_build(    rho_Hext(n), mla%la(n),     1, 0)
        !       call multifab_build(     thermal2(n), mla%la(n),     1, 1)
        !       if (dm .eq. 3) then
        !          call multifab_build(normal(n), mla%la(n),    dm, 1)
        !       end if
              
        !       call setval(      unew(n), ZERO, all=.true.)
        !       call setval(      snew(n), ZERO, all=.true.)
        !       call setval(    sponge(n), ONE,  all=.true.)
        !       call setval(     hgrhs(n), ZERO, all=.true.)
        !       call setval(Source_new(n), ZERO, all=.true.)
        !    end do

        !    ! Create normal now that we have defined center and dx
        !    call make_normal(normal,dx)

        !    if (spherical .eq. 0) then
        !       ! copy the old base state density with piecewise linear
        !       ! interpolated data in the new positions -- this is 
        !       ! only necessary for evolve_base_state = F and
        !       ! spherical = 0.
        !       rho0_old = rho0_temp

        !       ! zero out any data where there is no corresponding full
        !       ! state array
        !       do n=2,nlevs_radial
        !          do i=1,numdisjointchunks(n)
        !             if (i .eq. numdisjointchunks(n)) then
        !                do r=r_end_coord(n,i)+1,nr(n)-1
        !                   rho0_old(n,r) = 0.d0
        !                end do
        !             else
        !                do r=r_end_coord(n,i)+1,r_start_coord(n,i+1)-1
        !                   rho0_old(n,r) = 0.d0
        !                end do
        !             end if
        !          end do
        !       end do
        !    endif

        !    ! recompute p0 based on the new rho0 
        !    call compute_cutoff_coords(rho0_old)
           
        !    call make_grav_cell(grav_cell,rho0_old)

        !    ! enforce HSE
        !    call enforce_HSE(rho0_old,p0_old,grav_cell)

        !    if (use_tfromp) then
        !       ! compute full state T = T(rho,p0,X)
        !       call makeTfromRhoP(sold,p0_old,mla,the_bc_tower%bc_tower_array,dx)
        !    else
        !       ! compute full state T = T(rho,h,X)
        !       call makeTfromRhoH(sold,p0_old,mla,the_bc_tower%bc_tower_array,dx)
        !    end if

        !    ! force tempbar to be the average of temp
        !    call average(mla,sold,tempbar,dx,temp_comp)

        !    ! gamma1bar needs to be recomputed
        !    if (allocated(gamma1)) deallocate(gamma1)
        !    allocate(gamma1(nlevs))
           
        !    do n=1,nlevs
        !       call multifab_build(gamma1(n), mla%la(n), 1, 0)
        !    end do
           
        !    call make_gamma(mla,gamma1,sold,p0_old,dx)
        !    call average(mla,gamma1,gamma1bar,dx,1)
           
        !    do n=1,nlevs
        !       call destroy(gamma1(n))
        !    end do

        !    ! div_coeff_old needs to be recomputed
        !    call make_div_coeff(div_coeff_old,rho0_old,p0_old,gamma1bar,grav_cell)

        !    ! redistribute the particles to their new processor locations
        !    if (use_particles) call redistribute(particles,mla,dx,prob_lo)

        ! end if regridding ! end regridding
        ! !---------------------------------------------------------------------
        ! ! END REGRIDDING
        ! !---------------------------------------------------------------------

        !---------------------------------------------------------------------
        ! get the new timestep
        !---------------------------------------------------------------------
        dtold = dt

        if (istep > 1) then

           dt = 1.d20
           
           !call estdt(mla,the_bc_tower,uold,sold,gpi,Source_old,dSdt, &
           !           w0,rho0_old,p0_old,gamma1bar,grav_cell,dx,cflfac,dt)
           dt = fixed_dt

           if (parallel_IOProcessor() .and. verbose .ge. 1) then
              print*,''
              print*,"Call to estdt at beginning of step",istep
              print*,"gives dt =",dt
           end if

           if(dt .gt. max_dt) then
              if (parallel_IOProcessor() .and. verbose .ge. 1) then
                 print*,'max_dt limits the new dt =',max_dt
              end if
              dt = max_dt
           end if

           if(fixed_dt .ne. -1.0d0) then
              dt = fixed_dt
              if (parallel_IOProcessor()) then
                 print*, "Setting fixed dt =",dt
              end if
           end if

           if (parallel_IOProcessor() .and. verbose .ge. 1) then
              print*,''
           end if

        end if

        if (stop_time >= 0.d0) then
           if (time+dt > stop_time) then
              dt = stop_time - time 
              if (parallel_IOProcessor()) then
                 print*, "Stop time limits dt =",dt
              end if
           end if
        end if

        !---------------------------------------------------------------------
        ! Advance a single timestep at all levels.
        !---------------------------------------------------------------------
        runtime1 = parallel_wtime()

        ! call advance_timestep(init_mode,mla,uold,sold,unew,snew,gpi,pi,normal,rho0_old, &
        !                       rhoh0_old,rho0_new,rhoh0_new,p0_old,p0_new,tempbar,gamma1bar, &
        !                       w0,rho_omegadot2,rho_Hnuc2,rho_Hext,thermal2, &
        !                       div_coeff_old,div_coeff_new, &
        !                       grav_cell,dx,dt,dtold,the_bc_tower,dSdt,Source_old, &
        !                       Source_new,etarho_ec,etarho_cc,psi,sponge,hgrhs,tempbar_init, &
        !                       particles)

        runtime2 = parallel_wtime() - runtime1
        call parallel_reduce(runtime1, runtime2, MPI_MAX, proc = parallel_IOProcessorNode())
        if (parallel_IOProcessor()) print*,'Time to advance timestep: ',runtime1,' seconds'

        call print_and_reset_fab_byte_spread()

        time = time + dt

        if ( verbose .ge. 1 ) then

           if ( parallel_IOProcessor() ) then
              print *, 'MEMORY STATS AT END OF TIMESTEP ', istep
              print*, ' '
           end if
           call print(multifab_mem_stats(),    "    multifab")
           call print(fab_mem_stats(),         "         fab")
           call print(boxarray_mem_stats(),    "    boxarray")
           call print(layout_mem_stats(),      "      layout")
           call print(boxassoc_mem_stats(),    "    boxassoc")
           call print(fgassoc_mem_stats(),     "     fgassoc")
           call print(syncassoc_mem_stats(),   "   syncassoc")
           call print(copyassoc_mem_stats(),   "   copyassoc")
           call print(fluxassoc_mem_stats(),   "   fluxassoc")
           if ( parallel_IOProcessor() ) print*, ''

           do n = 1,nlevs
              if (parallel_IOProcessor()) write(6,1100) n

              smin = multifab_min_c(unew(n),1) 
              smax = multifab_max_c(unew(n),1)
              if (parallel_IOProcessor()) write(6,1101) smin,smax

              if (dm .ge. 2) then
                 smin = multifab_min_c(unew(n),2) 
                 smax = multifab_max_c(unew(n),2)
                 if (parallel_IOProcessor()) write(6,1102) smin,smax
              end if

              if (dm .eq. 3) then
                 smin = multifab_min_c(unew(n),3) 
                 smax = multifab_max_c(unew(n),3)
                 if (parallel_IOProcessor()) write(6,1103) smin,smax
              end if
              
              if (parallel_IOProcessor()) write(6,1107)
           end do
        end if

1100    format('At level ',i3)
1101    format('... min/max : x-velocity       ',e17.10,2x,e17.10)
1102    format('... min/max : y-velocity       ',e17.10,2x,e17.10)
1103    format('... min/max : z-velocity       ',e17.10,2x,e17.10)
1107    format(' ')

        if (parallel_IOProcessor()) then
           write(6,1000) istep,time,dt
           write(6,*)
        end if

        !---------------------------------------------------------------------
        ! save the old velocity, scalars and source terms for the next step
        !---------------------------------------------------------------------

        do n = 1,nlevs
           call multifab_copy_c(uold(n),      1,unew(n),      1,n_comp,   nghost(uold(n)))
           call multifab_copy_c(wold(n),      1,wnew(n),      1,n_comp,   nghost(uold(n)))
        end do

        !---------------------------------------------------------------------
        ! output
        !---------------------------------------------------------------------

        ! if the file .dump_checkpoint exists in our output directory, then
        ! automatically dump a plotfile
        inquire(file=".dump_checkpoint", exist=dump_checkpoint)

        ! if (chk_int > 0 .or. dump_checkpoint) then
        !    if (mod(istep,chk_int) .eq. 0 .or. dump_checkpoint) then

        !       ! write out any buffered diagnostic information
        !       call flush_diag()

        !       allocate(chkdata(nlevs))
        !       do n = 1,nlevs
        !          call multifab_build(chkdata(n), mla%la(n), 2*dm+nscal, 0)
        !          call multifab_copy_c(chkdata(n),1,unew(n),1,dm)
        !          call multifab_copy_c(chkdata(n),rho_comp+dm,snew(n),1,nscal)
        !       end do

        !       if (istep <= 99999) then
        !          write(unit=check_index,fmt='(i5.5)') istep
        !          check_file_name = trim(check_base_name) // check_index
        !       else if (istep <= 999999)  then
        !          write(unit=check_index6,fmt='(i6.6)') istep
        !          check_file_name = trim(check_base_name) // check_index6
        !       else
        !          write(unit=check_index7,fmt='(i7.7)') istep
        !          check_file_name = trim(check_base_name) // check_index7
        !       endif

        !       call checkpoint_write(check_file_name, chkdata, &
        !                             pi, dSdt, Source_old, Source_new, &
        !                             rho_omegadot2, rho_Hnuc2, rho_Hext, &
        !                             thermal2, mla%mba%rr, &
        !                             dt)

        !       call write_base_state(istep, check_file_name, &
        !                             rho0_new, rhoh0_new, p0_new, gamma1bar(:,:), &
        !                             w0, etarho_ec, etarho_cc, &
        !                             div_coeff_old, psi, tempbar, tempbar_init, prob_lo(dm))

        !       call write_aux_data(istep, check_file_name)

        !       last_chk_written = istep

        !       do n = 1,nlevs
        !          call destroy(chkdata(n))
        !       end do
        !       deallocate(chkdata)

        !       if (use_particles) then
        !          call particle_container_checkpoint(particles,check_file_name,mla)
        !       endif

        !    end if
        ! end if

        ! if the file .dump_plotfile exists in our output directory, then
        ! automatically dump a plotfile
        inquire(file=".dump_plotfile", exist=dump_plotfile)

        if (plot_int > 0 .or. plot_deltat > ZERO .or. dump_plotfile) then
           if ( (plot_int > 0 .and. mod(istep,plot_int) .eq. 0) .or. &
                (plot_deltat > ZERO .and. &
                mod(time - dt,plot_deltat) > mod(time,plot_deltat)) .or. &
                dump_plotfile) then

              if (istep <= 99999) then
                 write(unit=plot_index,fmt='(i5.5)') istep
                 plot_file_name = trim(plot_base_name) // plot_index
              else if (istep <= 999999) then
                 write(unit=plot_index6,fmt='(i6.6)') istep
                 plot_file_name = trim(plot_base_name) // plot_index6
              else
                 write(unit=plot_index7,fmt='(i7.7)') istep
                 plot_file_name = trim(plot_base_name) // plot_index7
              endif

              call make_plotfile(plot_file_name,mla,uold,wold,mla%mba,plot_names,dx,dt, &
                                 the_bc_tower,write_pf_time)

              !call write_job_info(plot_file_name, mla%mba, the_bc_tower, write_pf_time)
              last_plt_written = istep
           end if
        end if

        ! if the file .abort_maestro exists in our output directory, then
        ! automatically end the run.  This has the effect of also dumping
        ! a final checkpoint file.
        inquire(file=".abort_maestro", exist=abort_maestro)
        if (abort_maestro) exit

        ! have we reached the stop time?
        if (stop_time >= 0.d0) then
           if (time >= stop_time) goto 999
        end if

     end do  ! end main evolution loop

999  continue
     if (istep > max_step) istep = max_step


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write the final checkpoint and plotfile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

1000 format('STEP = ',i7,1x,' TIME = ',es16.10,1x,' DT = ',es16.10)

     ! if ( chk_int > 0 .and. last_chk_written .ne. istep ) then
     !    !       This writes a checkpoint file.

     !    ! write out any buffered diagnostic information
     !    call flush_diag()

     !    allocate(chkdata(nlevs))
     !    do n = 1,nlevs
     !       call multifab_build(chkdata(n), mla%la(n), 2*dm+nscal, 0)
     !       call multifab_copy_c(chkdata(n),1,unew(n),1,dm)
     !       call multifab_copy_c(chkdata(n),rho_comp+dm,snew(n),1,nscal)
     !       call multifab_copy_c(chkdata(n),rho_comp+dm+nscal,gpi(n),1,dm)
     !    end do

     !    if (istep <= 99999) then
     !       write(unit=check_index,fmt='(i5.5)') istep
     !       check_file_name = trim(check_base_name) // check_index
     !    else if (istep <= 999999) then
     !       write(unit=check_index6,fmt='(i6.6)') istep
     !       check_file_name = trim(check_base_name) // check_index6
     !    else
     !       write(unit=check_index7,fmt='(i7.7)') istep
     !       check_file_name = trim(check_base_name) // check_index7
     !    endif

     !    call checkpoint_write(check_file_name, chkdata, &
     !                          pi, dSdt, Source_old, Source_new, &
     !                          rho_omegadot2, rho_Hnuc2, rho_Hext, thermal2, &
     !                          mla%mba%rr, dt)

     !    call write_base_state(istep, check_file_name, &
     !                          rho0_new, rhoh0_new, p0_new, gamma1bar, &
     !                          w0, etarho_ec, etarho_cc, &
     !                          div_coeff_old, psi, tempbar, tempbar_init, prob_lo(dm))

     !    call write_aux_data(istep, check_file_name)

     !    do n = 1,nlevs
     !       call destroy(chkdata(n))
     !    end do
     !    deallocate(chkdata)

     !    if (use_particles) then
     !       call particle_container_checkpoint(particles,check_file_name,mla)
     !    endif

     ! end if

     if ( ( plot_int > 0 .or. &
            (plot_deltat > ZERO .and. &
             mod(time - dt,plot_deltat) > mod(time,plot_deltat)) ) .and. &
          last_plt_written .ne. istep ) then

        if (istep <= 99999) then
           write(unit=plot_index,fmt='(i5.5)') istep
           plot_file_name = trim(plot_base_name) // plot_index
        else if (istep <= 999999) then
           write(unit=plot_index6,fmt='(i6.6)') istep
           plot_file_name = trim(plot_base_name) // plot_index6
        else
           write(unit=plot_index7,fmt='(i7.7)') istep
           plot_file_name = trim(plot_base_name) // plot_index7
        endif

        call make_plotfile(plot_file_name,mla,uold,wold,mla%mba,plot_names,dx,dt, &
                           the_bc_tower,write_pf_time)

        !call write_job_info(plot_file_name, mla%mba, the_bc_tower, write_pf_time)
     end if
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! clean-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do n=1,nlevs
     call destroy(uold(n))
     call destroy(unew(n))
     call destroy(wold(n))
     call destroy(wnew(n))
  end do

  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)

  call simple_log_finalize()

  call runtime_close()

  deallocate(uold,wold)
  deallocate(dx)

end subroutine varden
