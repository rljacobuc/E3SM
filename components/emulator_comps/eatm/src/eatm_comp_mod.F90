module eatm_comp_mod

  ! !USES:

  use mct_mod
  use esmf
  use perf_mod
  use shr_const_mod
  use shr_sys_mod
  use shr_kind_mod   , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_file_mod   , only: shr_file_getunit, shr_file_freeunit
  use shr_cal_mod    , only: shr_cal_date2julian, shr_cal_ymdtod2string
  use shr_mpi_mod    , only: shr_mpi_bcast
  use seq_timemgr_mod, only: seq_timemgr_EClockGetData, seq_timemgr_RestartAlarmIsOn
  use eatmMod
  use eatmSpmdMod
  use eatmIO

  ! !PUBLIC TYPES:

  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: eatm_comp_init
  public :: eatm_comp_run
  public :: eatm_comp_final

  save

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !===============================================================================

  subroutine eatm_comp_init(Eclock, x2a, a2x, &
       seq_flds_x2a_fields, seq_flds_a2x_fields, &
       gsmap, ggrid, read_restart)

    ! !DESCRIPTION: initialize emulator atm model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)       , intent(in)    :: EClock
    type(mct_aVect)        , intent(inout) :: x2a, a2x
    character(len=*)       , intent(in)    :: seq_flds_x2a_fields ! fields from mediator
    character(len=*)       , intent(in)    :: seq_flds_a2x_fields ! fields to mediator
    type(mct_gsMap)        , pointer       :: gsMap               ! model global sep map (output)
    type(mct_gGrid)        , pointer       :: ggrid               ! model ggrid (output)
    logical                , intent(in)    :: read_restart        ! start from restart

    !--- local variables ---
    integer(IN)   :: n,k         ! generic counters
    integer(IN)   :: kmask       ! field reference
    integer(IN)   :: klat        ! field reference
    integer(IN)   :: kfld        ! fld index
    integer(IN)   :: cnt         ! counter
    integer(IN)   :: idt         ! integer timestep

    logical       :: exists      ! filename existance
    integer(IN)   :: nu          ! unit number
    integer(IN)   :: CurrentYMD  ! model date
    integer(IN)   :: CurrentTOD  ! model sec into model date
    integer(IN)   :: stepno      ! step number
    character(CL) :: calendar    ! calendar type
    character(CL) :: flds_strm

    ! TODO (AN): add this to the namelist
    character(CL)     :: ic_filename = "/global/cfs/cdirs/e3sm/anolan/ACE2-E3SMv3/initial_conditions/1971010100.nc"
    type(file_desc_t) :: ncid    ! netcdf file id

    !--- formats ---
    character(*), parameter :: F00   = "('(eatm_comp_init) ',8a)"
    character(*), parameter :: F0L   = "('(eatm_comp_init) ',a, l2)"
    character(*), parameter :: F01   = "('(eatm_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(eatm_comp_init) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(eatm_comp_init) ',a,i8,a)"
    character(*), parameter :: F04   = "('(eatm_comp_init) ',2a,2i8,'s')"
    character(*), parameter :: F05   = "('(eatm_comp_init) ',a,2f10.4)"
    character(*), parameter :: F90   = "('(eatm_comp_init) ',73('='))"
    character(*), parameter :: F91   = "('(eatm_comp_init) ',73('-'))"
    character(*), parameter :: subName = "(eatm_comp_init) "
    !-------------------------------------------------------------------------------

    call t_startf('EATM_INIT')

       call t_startf('eatm_initmctavs')
       if (masterproc) write(logunit_atm,F00) 'allocate AVs'
       call shr_sys_flush(logunit_atm)

       allocate(landfrac(lsize_x,lsize_y))
       allocate(ocnfrac(lsize_x,lsize_y))
       allocate(icefrac(lsize_x,lsize_y))
       allocate(phis(lsize_x,lsize_y))
       allocate(solin(lsize_x,lsize_y))

       allocate(ps(lsize_x,lsize_y))
       allocate(ts(lsize_x,lsize_y))
       allocate(t_0(lsize_x,lsize_y))
       allocate(specific_total_water_0(lsize_x,lsize_y))
       allocate(u_0(lsize_x,lsize_y))
       allocate(v_0(lsize_x,lsize_y))

       allocate(lhflx(lsize_x,lsize_y))
       allocate(shflx(lsize_x,lsize_y))
       allocate(surface_precipitation_rate(lsize_x,lsize_y))
       allocate(surface_upward_longwave_flux(lsize_x,lsize_y))
       allocate(flut(lsize_x,lsize_y))
       allocate(flds(lsize_x,lsize_y))
       allocate(fsds(lsize_x,lsize_y))
       allocate(surface_upward_shortwave_flux(lsize_x,lsize_y))
       allocate(top_of_atmos_upward_shortwave_flux(lsize_x,lsize_y))
       allocate(tendency_of_total_water_path_due_to_advection(lsize_x,lsize_y))

       ! determined from precip based surface temp
       call t_stopf('eatm_initmctavs')

       !----------------------------------------------------------------------------
       ! Read restart
       !----------------------------------------------------------------------------

       if (read_restart) then
          if (masterproc) then
             write(logunit_atm,F00) ' restart filename from rpointer'
             call shr_sys_flush(logunit_atm)
             inquire(file=trim(rpfile)//trim(inst_suffix),exist=exists)
             if (.not.exists) then
                write(logunit_atm,F00) ' ERROR: rpointer file does not exist'
                call shr_sys_abort(trim(subname)//' ERROR: rpointer file missing')
             endif
             nu = shr_file_getUnit()
             open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
             read(nu,'(a)') restart_file
             close(nu)
             call shr_file_freeUnit(nu)
          endif
          call shr_mpi_bcast(restart_file,mpicom_atm,'restart_file')

          call shr_sys_flush(logunit_atm)
       endif

       if (read_restart) then
          call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
          call seq_timemgr_EClockGetData( EClock, stepno=stepno, dtime=idt )
          call seq_timemgr_EClockGetData( EClock, calendar=calendar )
       endif

    !----------------------------------------------------------------------------
    ! Set initial atm state
    !----------------------------------------------------------------------------
    !JW IC file?
    ! populate all fields to zero OR from ACE IC if it's present
    call t_startf ('eatm_grid')

    call ncd_pio_openfile(ncid, trim(ic_filename), 0)

    ! AN: will do the actual readin in here. Need to pay attention to C vs. Fortran
    !     ordering for interoperability with FTorch
    call ncd_pio_closefile(ncid)


    call t_stopf ('eatm_grid')
    call t_stopf('EATM_INIT')

    return

  end subroutine eatm_comp_init

  !===============================================================================
  subroutine eatm_comp_run(EClock, x2a, a2x, gsmap, ggrid)

    ! !DESCRIPTION: run method for eatm model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)       , intent(in)    :: EClock
    type(mct_aVect)        , intent(inout) :: x2a
    type(mct_aVect)        , intent(inout) :: a2x
    type(mct_gsMap)        , pointer       :: gsMap
    type(mct_gGrid)        , pointer       :: ggrid

    !--- local ---
    integer(IN)   :: CurrentYMD        ! model date
    integer(IN)   :: CurrentTOD        ! model sec into model date
    integer(IN)   :: yy,mm,dd          ! year month day
    integer(IN)   :: n                 ! indices
    integer(IN)   :: idt               ! integer timestep
    real(R8)      :: dt                ! timestep
    logical       :: write_restart     ! restart now
    integer(IN)   :: nu                ! unit number
    integer(IN)   :: stepno            ! step number
    real(R8)      :: rday              ! elapsed day
    real(R8)      :: cosFactor         ! cosine factor
    real(R8)      :: factor            ! generic/temporary correction factor
    real(R8)      :: avg_alb           ! average albedo
    real(R8)      :: tMin              ! minimum temperature
    character(CL) :: calendar          ! calendar type

    character(len=18) :: date_str
    !--- temporaries
    real(R8)      :: uprime,vprime,swndr,swndf,swvdr,swvdf,ratio_rvrf
    real(R8)      :: tbot,pbot,rtmp,vp,ea,e,qsat,frac
    character(*), parameter :: F00   = "('(eatm_comp_run) ',8a)"
    character(*), parameter :: F04   = "('(eatm_comp_run) ',2a,2i8,'s')"
    character(*), parameter :: subName = "(eatm_comp_run) "
    !-------------------------------------------------------------------------------

    call t_startf('EATM_RUN')

    call t_startf('eatm_run1')

    call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
    call seq_timemgr_EClockGetData( EClock, curr_yr=yy, curr_mon=mm, curr_day=dd)
    call seq_timemgr_EClockGetData( EClock, stepno=stepno, dtime=idt)
    call seq_timemgr_EClockGetData( EClock, calendar=calendar)
    dt = idt * 1.0_r8
    write_restart = seq_timemgr_RestartAlarmIsOn(EClock)

    call t_stopf('eatm_run1')

    !--------------------
    ! ADVANCE ATM
    !--------------------

    call t_barrierf('eatm_BARRIER',mpicom_atm)
    call t_startf('eatm')

    !JW real work goes here


    call t_stopf('eatm_datamode')

    !--------------------
    ! Write restart
    !--------------------

    if (write_restart) then
       call t_startf('eatm_restart')
       call shr_cal_ymdtod2string(date_str, yy,mm,dd,currentTOD)

       write(restart_file,"(6a)") &
            trim(case_name), '.eatm',trim(inst_suffix),'.r.', trim(date_str), '.nc'
       if (masterproc) then
          nu = shr_file_getUnit()
          open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
          write(nu,'(a)') restart_file
          close(nu)
          call shr_file_freeUnit(nu)
       endif

       call shr_sys_flush(logunit_atm)
       call t_stopf('eatm_restart')
    endif

    call t_stopf('eatm')

    !----------------------------------------------------------------------------
    ! Log output for model date
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call t_startf('eatm_run2')
    if (masterproc) then
       write(logunit_atm,F04) trim(myModelName),': model date ', CurrentYMD,CurrentTOD
       call shr_sys_flush(logunit_atm)
    end if

    call t_stopf('eatm_run2')

    call t_stopf('EATM_RUN')

    return

  end subroutine eatm_comp_run

  !===============================================================================
  subroutine eatm_comp_final()

    ! !DESCRIPTION: finalize method for eatm model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    !--- formats ---
    character(*), parameter :: F00   = "('(eatm_comp_final) ',8a)"
    character(*), parameter :: F91   = "('(eatm_comp_final) ',73('-'))"
    character(*), parameter :: subName = "(eatm_comp_final) "
    !-------------------------------------------------------------------------------

    call t_startf('EATM_FINAL')

    ! deallocate arrays from init step
    deallocate(landfrac)
    deallocate(ocnfrac)
    deallocate(icefrac)
    deallocate(phis)
    deallocate(solin)

    deallocate(ps)
    deallocate(ts)
    deallocate(t_0)
    deallocate(specific_total_water_0)
    deallocate(u_0)
    deallocate(v_0)

    deallocate(lhflx)
    deallocate(shflx)
    deallocate(surface_precipitation_rate)
    deallocate(surface_upward_longwave_flux)
    deallocate(flut)
    deallocate(flds)
    deallocate(fsds)
    deallocate(surface_upward_shortwave_flux)
    deallocate(top_of_atmos_upward_shortwave_flux)
    deallocate(tendency_of_total_water_path_due_to_advection)

    if (masterproc) then
       write(logunit_atm,F91)
       write(logunit_atm,F00) trim(myModelName),': end of main integration loop'
       write(logunit_atm,F91)
    end if

    call t_stopf('EATM_FINAL')

  end subroutine eatm_comp_final
  !===============================================================================

end module eatm_comp_mod
