
PROGRAM FLake1D

!------------------------------------------------------------------------------
!
! Description:
!
!  A Main Program to run FLake.
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================



!==============================================================================
! Modules used:
!==============================================================================

USE data_parameters , ONLY : &
    ireals,                  & ! KIND-type parameter for real variables
    iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_derivedtypes         ! Definitions of several derived TYPEs

USE flake_parameters           ! Thermodynamic parameters and dimensionless constants of FLake
USE flake_configure, ONLY: sediments_on => lflk_botsed_use, &  ! Switches that configure FLake
                           nmlfile

USE flake_albedo_ref           ! Reference values the albedo for the lake water, lake ice and snow

USE flake, ONLY :    & 
  h_Snow_min_flk   , & ! Minimum snow thickness [m]
  h_Ice_min_flk    , & ! Minimum ice thickness [m]
  Q_snow_flk       , & ! Heat flux through the air-snow interface [W m^{-2}]
  Q_ice_flk        , & ! Heat flux through the snow-ice or air-ice interface [W m^{-2}]
  Q_w_flk          , & ! Heat flux through the ice-water or air-water interface [W m^{-2}]
  Q_bot_flk        , & ! Heat flux through the water-bottom sediment interface [W m^{-2}]
  I_snow_flk       , & ! Radiation flux through the air-snow interface [W m^{-2}]
  I_ice_flk        , & ! Radiation flux through the snow-ice or air-ice interface [W m^{-2}]
  I_w_flk          , & ! Radiation flux through the ice-water or air-water interface [W m^{-2}
  I_h_flk          , & ! Radiation flux through the mixed-layer-thermocline interface [W m^{-2}]
!_gk_db: the integral radfluxes used only in output 
  I_intm_0_h_flk   , & ! Integral-mean radiation flux over the mixed layer [W m^{-1}]
  I_intm_h_D_flk   , & ! Integral-mean radiation flux over the thermocline [W m^{-1}]
!_gk_db end
  I_bot_flk        , & ! Radiation flux through the water-bottom sediment interface [W m^{-2}]
  Q_star_flk       , & ! A generalized heat flux scale [W m^{-2}]
  u_star_w_flk     , & ! Friction velocity in the surface layer of lake water [m s^{-1}]
  w_star_sfc_flk       ! Convective velocity scale, using a generalized heat flux scale [m s^{-1}]

USE SfcFlx, ONLY :         &  
  u_star_min_sf          , & ! Minimum value of friction velocity [m s^{-1}]
  tpsf_c_a_p             , & ! Specific heat of air at constant pressure [J kg^{-1} K^{-1}]
  tpsf_L_evap            , & ! Specific heat of evaporation [J kg^{-1}]
  tpsf_R_dryair          , & ! Gas constant for dry air [J kg^{-1} K^{-1}]
  tpsf_Rd_o_Rv           , & ! Ratio of gas constants (Rd/Rv)
  tpsf_alpha_q           , & ! Diemsnionless ratio
  c_Karman               , & ! The von Karman constant
  Sc_neutral             , & ! Turbulent Schmidt number at neutral static stability
  c_MO_q_stab            , & ! Constant of the MO theory (humidity, stable stratification)
  SfcFlx_rhoair          , & ! Function that returns air density
  SfcFlx_lwradatm        , & ! Function that returns long-wave radiation flux from the atmosperic 
  SfcFlx_satwvpres       , & ! Function that returns saturation water vapour pressure
  SfcFlx_wvpreswetbulb   , & ! Function that returns water vapour pressure using wet-buld temperature
  SfcFlx_spechum         , & ! Function that returns specific humidity
  SfcFlx_lwradwsfc       , & ! Function, returns the surface long-wave radiation flux
  SfcFlx_momsenlat       , & ! Subroutine, computes fluxes of momentum and of sensible and latent heat
  P_a_ref                , & ! Reference pressure [N m^{-2} = kg m^{-1} s^{-2}]
  u_star_a_sf            , & ! Friction velocity [m s^{-1}]
  Q_mom_a_sf             , & ! Momentum flux [N m^{-2}]
  Q_sens_a_sf            , & ! Sensible heat flux [W m^{-2}]
  Q_lat_a_sf             , & ! Laten heat flux [W m^{-2}]
  Q_watvap_a_sf              ! Flux of water vapout [kg m^{-2} s^{-1}]

!==============================================================================
! End of Modules used
!==============================================================================


!==============================================================================
! Declarations
!==============================================================================
IMPLICIT NONE

! Parameters of type INTEGER
INTEGER (KIND = iintegers), PARAMETER :: &
!! GK: array dimensions are static in FLake. It is a bug. The program will fail to work for long periods. 10^6 is used as a plug
  dim_out_arrays = 1000000_iintegers      , &
  n_lev_wind_m   = 1_iintegers           , & ! Number of levels where wind speed is measured
  n_lev_Taqa_m   = 1_iintegers           , & ! Number of levels where T_a and q_a (or T_webulb) are measured
  n_lev_Tw_m     = 1_iintegers               ! Number of levels where T_w is measured


! Parameters of type REAL
REAL (KIND = ireals), PARAMETER ::    &
  pi_value = 3.1415927              , & ! The value of \Pi
  hour_in_sec     = 3600._ireals    , & ! One hour in seconds [s] 
  day_length      = 86400._ireals   , & ! The number of seconds in 24 hours [s]
  small_number    = 1.E-10_ireals   , & ! A small number
  large_number    = 9999.0_ireals   , & ! A large number, used to mark wrong values (c/o FB)
  omega_earth     = 7.29E-05_ireals     ! The angular velocity of the earth's rotation [s^{-1}]



! Variables of type CHARACTER
CHARACTER (LEN = 10, KIND = 1)  ::    &  
  current_date                      , & ! Current date
  current_time						    ! Current time

CHARACTER (LEN = 300, KIND = 1)  ::    &  
  meteofile							, &	!Input file name 
  outputfile							!Output file name
!  Variables of type LOGICAL
LOGICAL ::          &
  l_save              ! Switch, TRUE = save results of calculation
 

! Variables of type INTEGER
INTEGER (KIND = iintegers) :: &
  ntsc                      , & ! The number of current time step
  time_step_number          , & ! The total number time steps
  n_out_save                , & ! The number of saved outputs
  save_interval_s           , & ! Saving interval [s]
  save_interval_n           , & ! Saving interval in time steps
  i, j                      , & ! Loop indices
  nband_optic				    ! Number of wave-length bands
  
! Variables of type REAL
REAL (KIND = ireals) ::    &
  current_time_s         , & ! Current simulation time [s]
  comput_length          , & ! The length of simulation [s]
! Optical characteristics of water, ice and snow 
  albedo_water           , & ! Water surface albedo with respect to the short-wave radiation
  albedo_ice             , & ! Ice surface albedo with respect to the short-wave radiation
  albedo_snow            , & ! Snow surface albedo with respect to the short-wave radiation
! Light extinction parameters
  frac_optic(10)		 , & ! Fractions of total radiation flux 
  extincoef_optic(10)    , & ! Extinction coefficients 
! Lake-specific parameters
  depth_w_lk             , & ! The lake depth [m]
  fetch_lk               , & ! Typical wind fetch [m] 
  depth_bs_lk            , & ! Depth of the thermally active layer of the bottom sediments [m]
  T_bs_lk                , & ! Temperature at the outer edge of 
                             ! the thermally active layer of bottom sediments [K]
  latitude_lk            , & ! Geographical latitude [degrees]
  par_Coriolis_lk        , & ! The Coriolis parameter [s^{-1}]
  del_time_lk            , & ! The model time step [s]
! Measurement heights
  height_u              , &  ! Height where wind is measured [m]
  height_tq             , &  ! Height where temperature and humidity are measured [m]
! Variables at the previous (in) time step, and updated variable (out) 
  T_snow_in, T_snow_out  , & ! Snow temperature [K]
  T_ice_in,  T_ice_out   , & ! Ice temperature [K]
  T_mnw_in,  T_mnw_out   , & ! Mean temperature of the water column [K]
  T_wML_in,  T_wML_out   , & ! Mixed-layer temperature [K]
  T_bot_in,  T_bot_out   , & ! Bottom temperature [K]
  T_B1_in,   T_B1_out    , & ! Temperature at the bottom of the upper layer of the sediment [K]
  C_T_in,    C_T_out     , & ! Shape factor (thermocline)
  h_snow_in, h_snow_out  , & ! Snow thickness [m]
  h_ice_in,  h_ice_out   , & ! Ice thickness [m]
  h_ML_in,   h_ML_out    , & ! Mixed-layer thickness [m] 
  H_B1_in,   H_B1_out    , & ! Thickness of the upper layer of the sediment [m]
  T_sfc_in,  T_sfc_out       ! Surface temperature [K]

! Arrays of type REAL
REAL (KIND = ireals), DIMENSION(0:dim_out_arrays) ::   &
  simul_time             , & ! Simulation time [h] 
  T_snow                 , & ! Snow temperature [K]
  T_ice                  , & ! Ice temperature [K]
  T_mnw                  , & ! Mean temperature of the water column [K]
  T_wML                  , & ! Mixed-layer temperature [K]
  T_bot                  , & ! Bottom temperature [K]
  T_B1                   , & ! Temperature at the bottom of the upper layer of the sediment [K]
  C_T                    , & ! Shape factor (thermocline)
  h_snow                 , & ! Snow thickness [m]
  h_ice                  , & ! Ice thickness [m]
  h_ML                   , & ! Mixed-layer thickness [m] 
  H_B1                   , & ! Thickness of the upper layer of the sediment [m]
  T_sfc                  , & ! Surface temperature [K]
  dMsnowdt               , & ! The rate of snow accumulation [kg m^{-2} s^{-1}]
  Q_mom                  , & ! Momentum flux [N m^{-2}]
  u_star_a               , & ! Friction velocity in the surface air layer [m s^{-1}]
  u_star_w               , & ! Friction velocity in the surface layer of lake water [m s^{-1}]
  w_star_gen             , & ! Convective velocity scale, using a generalized heat flux scale [m s^{-1}]
  Q_sen                  , & ! Sensible heat flux [W m^{-2}]
  Q_lat                  , & ! Latent heat flux [W m^{-2}]
  Q_wv                   , & ! Flux of water vapour [kg m^{-2} s^{-1}]
  Q_atm_lw               , & ! Long-wave radiation flux from the atmosphere [W m^{-2}]
  I_ice                  , & ! Radiation flux through the snow-ice or air-ice interface [W m^{-2}]
  I_w                    , & ! Radiation flux through the ice-water or air-water interface [W m^{-2}
  Q_ice                  , & ! Heat flux through the snow-ice or air-ice interface [W m^{-2}]
  Q_w                    , & ! Heat flux through the ice-water or air-water interface [W m^{-2}]
  Q_bot                      ! Heat flux through the water-bottom sediment interface [W m^{-2}]

! Arrays of type REAL
REAL (KIND = ireals) ::          &
  z_wind_m (1:n_lev_wind_m)    , & ! Levels where wind speed is measured [m]
  z_Taqa_m (1:n_lev_Taqa_m)    , & ! Levels where T_a and q_a (or T_webulb) are measured [m]
  z_Tw_m   (1:n_lev_Tw_m)          ! Levels where T_w is measured [m]

!Variables of derived types
TYPE (opticpar_medium) :: &
  opticpar_water        , & ! Optical characteristics of water
  opticpar_ice          , & ! Optical characteristics of ice
  opticpar_snow             ! Optical characteristics of snow

! ALLOCATABLE ARRAYS 
REAL (KIND = ireals), ALLOCATABLE  ::       &
  U_dir(:)                                , & ! Wind direction from sonic [dgr]
  u_star_m(:)                             , & ! Friction velocity from sonic [m s^{-1}]
  Q_sensible_m(:)                         , & ! Sensible heat flux from sonic (corrected) [W m^{-2}] 
  Q_sens_sonic_m(:)                       , & ! Sensible heat flux from sonic (not corrected) [W m^{-2}] 
  Q_latent_m(:)                           , & ! Latent heat flux [W m^(-2)]
  I_solar_in_m(:)                         , & ! Downward short-wave radiation flux [W m^{-2}]
  I_solar_out_m(:)                        , & ! Upward short-wave radiation flux [W m^{-2}]
  T_radsfc_m(:)                           , & ! Radiation surface temperature [dgr C, K]
  T_wmean_m (:)                           , & ! Mean water temperature over the upper 1 m [K]
  rho_a_m(:)                              , & ! Air density [kg m^{-3}]
  cl_low_m(:)                             , & ! Low-level cloud cover [0,1]
  Q_atm_lw_m(:)                           , & ! Long-wave radiation flux from the atmosphere [W m^{-2}]
  U_a_m  (:,:)							  , & ! Wind speed [m s^(-1)]
  T_a_m  (:,:)							  , & ! Air temperature [dgr C, K]
  T_wb_m (:,:)							  , & ! Wet buld temperature [dgr C, K]
  e_a    (:,:)							  , & ! Water vapour pressure [N m^{-2} = kg m^{-1} s^{-2}]
  q_a    (:,:)							  , & ! Specific humidity
  T_w_m  (:,:)							  , & ! Water temperature buld temperature [dgr C, K]
  P_a_m  (:)                              , & ! Air pressure [N m^{-2} = kg m^{-1} s^{-2}]
  cl_tot_m (:)                                ! Total cloud cover [0,1]

! Input namelist descriptions
NAMELIST /METEO/ z_wind_m, z_Taqa_m, z_Tw_m, meteofile, outputfile
NAMELIST /LAKE_PARAMS/ depth_w_lk, fetch_lk, depth_bs_lk, T_bs_lk, latitude_lk, sediments_on, c_relax_C, C_T_max
NAMELIST /SIMULATION_PARAMS/ del_time_lk,time_step_number, save_interval_n, &
	T_wML_in, T_bot_in, h_ML_in
NAMELIST /TRANSPARENCY/ nband_optic, frac_optic, extincoef_optic
!==============================================================================
!  END OF DECLARATIONS
!==============================================================================
!
!==============================================================================
!  START CALCULATIONS
!==============================================================================
!
!   
  T_wML_in = tpl_T_r-tpl_T_f ! initial temperature set to the maximum dens. val. (C)
                             ! (can be overwritten by reading the value from .nml file below)
  T_bot_in = tpl_T_r-tpl_T_f ! initial temperature set to the maximum dens. val. (C)
  h_ML_in =  3._ireals       ! initial mixed layer depth
!  Read configuration file (name from command line)
CALL GET_COMMAND_ARGUMENT(1,nmlfile)

OPEN(UNIT=11, FILE=nmlfile)
READ(UNIT=11, NML=SIMULATION_PARAMS)
READ(UNIT=11, NML=METEO)
READ(UNIT=11, NML=LAKE_PARAMS)
READ(UNIT=11, NML=TRANSPARENCY)
CLOSE(UNIT=11)

! Allocate arrays
ALLOCATE (U_dir(0:time_step_number)		   , &
  u_star_m(0:time_step_number)             , & ! Friction velocity from sonic [m s^{-1}]
  Q_sensible_m(0:time_step_number)         , & ! Sensible heat flux from sonic (corrected) [W m^{-2}] 
  Q_sens_sonic_m(0:time_step_number)       , & ! Sensible heat flux from sonic (not corrected) [W m^{-2}] 
  Q_latent_m(0:time_step_number)           , & ! Latent heat flux [W m^(-2)]
  I_solar_in_m(0:time_step_number)         , & ! Downward short-wave radiation flux [W m^{-2}]
  I_solar_out_m(0:time_step_number)		   , & ! Upward short-wave radiation flux [W m^{-2}]
  T_radsfc_m(0:time_step_number)           , & ! Radiation surface temperature [dgr C, K]
  T_wmean_m (0:time_step_number)		   , & ! Mean water temperature over the upper 1 m [K]
  rho_a_m(0:time_step_number)			   , & ! Air density [kg m^{-3}]
  cl_low_m(0:time_step_number)			   , & ! Low-level cloud cover [0,1]
  Q_atm_lw_m(0:time_step_number)		   , & ! Long-wave radiation flux from the atmosphere [W m^{-2}]
  U_a_m  (1:n_lev_wind_m, 0:time_step_number)    , & ! Wind speed [m s^(-1)]
  T_a_m  (1:n_lev_Taqa_m, 0:time_step_number)    , & ! Air temperature [dgr C, K]
  T_wb_m (1:n_lev_Taqa_m, 0:time_step_number)    , & ! Wet buld temperature [dgr C, K]
  e_a    (1:n_lev_Taqa_m, 0:time_step_number)    , & ! Water vapour pressure [N m^{-2} = kg m^{-1} s^{-2}]
  q_a    (1:n_lev_Taqa_m, 0:time_step_number)    , & ! Specific humidity
  T_w_m  (1:n_lev_Tw_m,   0:time_step_number)    , &  ! Water temperature buld temperature [dgr C, K]
  P_a_m  (0:time_step_number+12_iintegers)       , & ! Air pressure [N m^{-2} = kg m^{-1} s^{-2}]
  cl_tot_m(0:time_step_number+12_iintegers)        & ! Total cloud cover [0,1]
)


! Read input data
OPEN(UNIT=11, FILE=meteofile, STATUS='old')
DO i=0_iintegers, time_step_number-1_iintegers
  READ(11,*) I_solar_in_m(i),             & 
             T_a_m(1,i),                  &
             e_a(1,i),                    & 
             U_a_m(1,i),                  &
             cl_tot_m(i)
END DO 
CLOSE(UNIT=11)

! Compute miscellaneous quantities
comput_length = time_step_number*del_time_lk          ! The length of simulation [s]

par_Coriolis_lk =        &							  ! The Coriolis parameter [s^{-1}]
	2._ireals*omega_earth*SIN(latitude_lk/180._ireals*pi_value)

T_bs_lk = T_bs_lk + tpl_T_f							  ! Temperature at the outer edge of
													  ! the thermally active layer of the bottom sediments [convert to K]

opticpar_water = opticpar_medium(nband_optic,       & !water transparency
								   frac_optic,      &
								   extincoef_optic) 

!Initialize some arrays
DO i=0_iintegers, time_step_number-1_iintegers
  T_w_m(1,0_iintegers) = T_wML_in + tpl_T_f           !(Initial )Water temperature  
  P_a_m(i) = P_a_ref                                  ! Use reference pressure
  DO j=1_iintegers, n_lev_Taqa_m                     
    T_a_m(j,i)  = T_a_m(j,i)  + tpl_T_f               ! convert air temperature (°C to K)
	e_a(j,i)    = e_a(j,i)*100                        ! convert air humidity (mb to N/m^2)
  END DO 
  I_solar_in_m(i) = MAX(0._ireals, I_solar_in_m(i))   ! Security, remove negative values
END DO 

!  Set initial values (no snow-ice cover, mixing down to the bottom)
n_out_save            = 0_iintegers
T_snow (n_out_save)   = tpl_T_f
T_ice  (n_out_save)   = tpl_T_f
h_snow (n_out_save)   = 0._ireals
h_ice  (n_out_save)   = 0._ireals
T_wML  (n_out_save)   = T_w_m(1,n_out_save)   
h_ML   (n_out_save)   = h_ML_in
T_bot  (n_out_save)   = T_bot_in + tpl_T_f
C_T    (n_out_save)   = C_T_min 
T_mnw  (n_out_save)   = T_wML(n_out_save) - C_T(n_out_save)   & 
                      * (1._ireals-h_ML(n_out_save)/depth_w_lk)*(T_wML(n_out_save)-T_bot(n_out_save))
T_B1   (n_out_save)   = T_bs_lk 
H_B1   (n_out_save)   = depth_bs_lk
T_sfc(n_out_save)     = T_wML(n_out_save)

! Initialise FLake "*_in" variables
T_snow_in = T_snow  (n_out_save)
T_ice_in  = T_ice   (n_out_save)
T_mnw_in  = T_mnw   (n_out_save)
T_wML_in  = T_wML   (n_out_save)
T_bot_in  = T_bot   (n_out_save)
T_B1_in   = T_B1    (n_out_save)
C_T_in    = C_T     (n_out_save)
h_snow_in = h_snow  (n_out_save)
h_ice_in  = h_ice   (n_out_save)
h_ML_in   = h_ML    (n_out_save)
H_B1_in   = H_B1    (n_out_save)
T_sfc_in  = T_sfc   (n_out_save)


!------------------------------------------------------------------------------
!  Call "flake_driver" and other required procedures in a DO loop
!------------------------------------------------------------------------------

DO ntsc = 0_iintegers, time_step_number-1_iintegers
! Set current time 
  current_time_s = ntsc*del_time_lk    

! Compute water vapour pressure and air specific humidity 
  DO j=1_iintegers,n_lev_Taqa_m
    q_a(j,ntsc) = SfcFlx_spechum (e_a(j,ntsc), P_a_m(ntsc))
  END DO
! Compute air density, using measured values
  rho_a_m(ntsc) = SfcFlx_rhoair(T_a_m(1,ntsc), q_a(1,ntsc), P_a_m(ntsc))
! The rate of snow accumulation [kg m^{-2} s^{-1}]
  dMsnowdt(ntsc) = 0._ireals
! Long-wave radiation
!!GK 
!! here is the trick: 
!! if input data are in the range 0-1, it is cloudiness, and LW rdiation should be calculated
!! if input >> 1 (the atmospheric LW _IS_ >> 1, as long as we have our earth atmosphere), 
!! the radiation goes directy from the input to FLake 
!! NB: the trick is _BUGGY_, can produce weird behaviour with random input!
if(cl_tot_m(ntsc).gt. 1.0) then
	Q_atm_lw_m(ntsc) = cl_tot_m(ntsc)
else 
	Q_atm_lw_m(ntsc) = &
	    -SfcFlx_lwradatm(T_a_m(1,ntsc),&
						e_a(1,ntsc),&
						cl_tot_m(ntsc),&
						cl_tot_m(ntsc))
endif
! Call FLake
CALL flake_interface ( dMsnowdt(ntsc), I_solar_in_m(ntsc), Q_atm_lw_m(ntsc), z_wind_m(1), z_Taqa_m(1),  &
                       U_a_m(1,ntsc), T_a_m(1,ntsc), e_a(1,ntsc), q_a(1,ntsc), P_a_m(ntsc),             &

                       depth_w_lk, fetch_lk, depth_bs_lk, T_bs_lk, par_Coriolis_lk, del_time_lk,   &
                       T_snow_in,  T_ice_in,  T_mnw_in,  T_wML_in,  T_bot_in,  T_B1_in,            &
                       C_T_in,  h_snow_in,  h_ice_in,  h_ML_in,  H_B1_in, T_sfc_in,                &

                       albedo_water,   albedo_ice,   albedo_snow,                              &
                       opticpar_water, opticpar_ice, opticpar_snow,                            &

                       T_snow_out, T_ice_out, T_mnw_out, T_wML_out, T_bot_out, T_B1_out,   &
                       C_T_out, h_snow_out, h_ice_out, h_ML_out, H_B1_out, T_sfc_out )
! The FLake output is saved every "save_interval_n" time step
  l_save = MOD(ntsc, save_interval_n).EQ.0_iintegers
  IF(l_save) THEN                                         ! Save output
    n_out_save = ntsc/save_interval_n                     ! The number of the output to save 
    simul_time (n_out_save) = current_time_s/day_length   ! Time in days  
! FLake prognstic variables, some FLake diagnostic variables
    T_snow     (n_out_save) = T_snow_in 
    T_ice      (n_out_save) = T_ice_in 
    T_mnw      (n_out_save) = T_mnw_in 
    T_wML      (n_out_save) = T_wML_in 
    T_bot      (n_out_save) = T_bot_in 
    T_B1       (n_out_save) = T_B1_in 
    C_T        (n_out_save) = C_T_in 
    h_snow     (n_out_save) = h_snow_in 
    h_ice      (n_out_save) = h_ice_in 
    h_ML       (n_out_save) = h_ML_in 
    H_B1       (n_out_save) = H_B1_in 
    T_sfc      (n_out_save) = T_sfc_in 
! Fluxes and related quantities (heat fluxes are positive downwards)
    u_star_a   (n_out_save) = u_star_a_sf
    u_star_w   (n_out_save) = u_star_w_flk
    w_star_gen (n_out_save) = w_star_sfc_flk
    Q_sen      (n_out_save) = -Q_sens_a_sf
    Q_lat      (n_out_save) = -Q_lat_a_sf  
    I_ice      (n_out_save) = I_ice_flk
    I_w        (n_out_save) = I_w_flk
    Q_ice      (n_out_save) = Q_ice_flk
    Q_w        (n_out_save) = Q_w_flk
    Q_bot      (n_out_save) = Q_bot_flk
  END IF

! Assign "*_out" --> "*_in"
  T_snow_in = T_snow_out    
  T_ice_in  = T_ice_out    
  T_mnw_in  = T_mnw_out    
  T_wML_in  = T_wML_out    
  T_bot_in  = T_bot_out    
  T_B1_in   = T_B1_out    
  C_T_in    = C_T_out    
  h_snow_in = h_snow_out    
  h_ice_in  = h_ice_out    
  h_ML_in   = h_ML_out    
  H_B1_in   = H_B1_out    
  T_sfc_in  = T_sfc_out

! Increment time step counter and set current time
!  ntsc = ntsc + 1_iintegers 

END DO 

!------------------------------------------------------------------------------
!  Save output
!------------------------------------------------------------------------------


!_ Results from simulations
OPEN(UNIT=11, FILE=outputfile, STATUS='replace')
  WRITE(11,610)
  DO i=0_iintegers, n_out_save
    WRITE(11,520) i, simul_time(i),                                             &
                  T_wML(i)-tpl_T_f, T_mnw(i)-tpl_T_f, T_bot(i)-tpl_T_f,         &
                  u_star_a(i), u_star_w(i), w_star_gen(i),                      &
                  Q_w(i), Q_sen(i), Q_lat(i),                                   &
                  I_w(i), Q_atm_lw(i), -SfcFlx_lwradwsfc(T_sfc(i)),             &
                  h_ML(i), C_T(i), H_B1(i), T_B1(i)-tpl_T_f, Q_bot(i),			&
				  h_ice(i),h_snow(i),T_ice(i)-tpl_T_f,T_snow(i)-tpl_T_f
  END DO
CLOSE(UNIT=11)


WRITE(*, '(" End of calculation, number of time steps: ", F10.0)' ) current_time_s/day_length
WRITE(*, '(" Results are written to:  ",A30)') outputfile
!==============================================================================
!  End of calculations
!==============================================================================
DEALLOCATE(U_dir   , & ! Wind direction
  u_star_m		   , & ! Friction velocity from sonic [m s^{-1}]
  Q_sensible_m     , & ! Sensible heat flux from sonic (corrected) [W m^{-2}] 
  Q_sens_sonic_m   , & ! Sensible heat flux from sonic (not corrected) [W m^{-2}] 
  Q_latent_m       , & ! Latent heat flux [W m^(-2)]
  I_solar_in_m     , & ! Downward short-wave radiation flux [W m^{-2}]
  I_solar_out_m	   , & ! Upward short-wave radiation flux [W m^{-2}]
  T_radsfc_m       , & ! Radiation surface temperature [dgr C, K]
  T_wmean_m		   , & ! Mean water temperature over the upper 1 m [K]
  rho_a_m		   , & ! Air density [kg m^{-3}]
  cl_low_m		   , & ! Low-level cloud cover [0,1]
  Q_atm_lw_m	   , & ! Long-wave radiation flux from the atmosphere [W m^{-2}]
    U_a_m          , & ! Wind speed [m s^(-1)]
  T_a_m            , & ! Air temperature [dgr C, K]
  T_wb_m		   , & ! Wet buld temperature [dgr C, K]
  e_a			   , & ! Water vapour pressure [N m^{-2} = kg m^{-1} s^{-2}]
  q_a		       , & ! Specific humidity
  T_w_m			   , &  ! Water temperature buld temperature [dgr C, K]
  P_a_m            , & ! Air pressure [N m^{-2} = kg m^{-1} s^{-2}]
  cl_tot_m           & ! Total cloud cover [0,1]
)

!------------------------------------------------------------------------------
!  FORMATs
!------------------------------------------------------------------------------

520  FORMAT(1X, I8, 25(1X, G13.6))

610  FORMAT('Results from FLAKE simulations.', / &
            ' No   ' ,     '  time        ',                                 &
            ' Ts           ', ' Tm           ', ' Tb           ',            &
            ' ufr_a        ', ' ufr_w        ', ' Wconv        ', ' Qw           ',&
            ' Q_se         ', ' Q_la         ',                              &
            ' I_w          ', ' Q_lwa        ', ' Q_lww        ',            &
            ' h_ML         ', ' C_T          ',                              &
            ' H_B1         ', ' T_B1         ', ' Qbot         ',            &
            ' H_ice        ', ' H_snow       ', ' T_ice        ',' T_snow   ')


STOP
END PROGRAM FLake1D