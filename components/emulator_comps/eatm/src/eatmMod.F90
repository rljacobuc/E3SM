module eatmMod

  ! !USES:

  use shr_kind_mod   , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL

  ! !PUBLIC TYPES:

  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public module data
  !--------------------------------------------------------------------------
  integer, public           :: gsize, lsize, lsize_x, lsize_y
  character(CL), public     :: restart_file
  character(CL), public     :: case_name      ! case name
  character(len=16), public :: inst_name
  character(len=16), public :: inst_suffix = ""    ! char string associated with instance (ie. "_0001" or "")

     !JW TODO: load up all arrays into a big 3D container? lots of 2D arrays?
     !         for now, using 2D arrays with EAM naming convention
     !         and EAM sign conventions

  ! variables that only go **into** ACE.
  real(kind=R8), dimension(:,:), allocatable, public :: landfrac               ! land area fraction
  real(kind=R8), dimension(:,:), allocatable, public :: ocnfrac                ! ocean area fraction
  real(kind=R8), dimension(:,:), allocatable, public :: icefrac                ! sea-ice area fraction
  real(kind=R8), dimension(:,:), allocatable, public :: phis                   ! surface geopotential
  real(kind=R8), dimension(:,:), allocatable, public :: solin                  ! solar insulation

  ! variables go **into** and **out** of ACE.  Excluding levels 1--8 fields for u, v, t, and specific_total_water
  real(kind=R8), dimension(:,:), allocatable, public :: ps                     ! surface pressure
  real(kind=R8), dimension(:,:), allocatable, public :: ts                     ! surface temperature (radiative)
  real(kind=R8), dimension(:,:), allocatable, public :: t_0                    ! temperature level-0
  real(kind=R8), dimension(:,:), allocatable, public :: specific_total_water_0 ! specific total water level-0
  real(kind=R8), dimension(:,:), allocatable, public :: u_0                    ! zonal wind level-0
  real(kind=R8), dimension(:,:), allocatable, public :: v_0                    ! meridional wind level-0

  ! variables that only come **out** of ACE.
  real(kind=R8), dimension(:,:), allocatable, public :: lhflx                  ! surface latent heat flux
  real(kind=R8), dimension(:,:), allocatable, public :: shflx                  ! surface sensible heat flux
  real(kind=R8), dimension(:,:), allocatable, public :: surface_precipitation_rate
  real(kind=R8), dimension(:,:), allocatable, public :: surface_upward_longwave_flux
  real(kind=R8), dimension(:,:), allocatable, public :: flut                   ! upwelling longwave flux at top of model
  real(kind=R8), dimension(:,:), allocatable, public :: flds                   ! downwelling longwave flux at top of model
  real(kind=R8), dimension(:,:), allocatable, public :: fsds                   ! downwelling solar flux at surace
  real(kind=R8), dimension(:,:), allocatable, public :: surface_upward_shortwave_flux
  real(kind=R8), dimension(:,:), allocatable, public :: top_of_atmos_upward_shortwave_flux
  real(kind=R8), dimension(:,:), allocatable, public :: tendency_of_total_water_path_due_to_advection

  character(CS), public :: myModelName = 'atm'   ! user defined model name

  character(len=*), parameter, public :: rpfile = 'rpointer.atm'

end module eatmMod
