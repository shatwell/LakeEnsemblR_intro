! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

MODULE SfcFlx

!------------------------------------------------------------------------------
!
! Description:
!
!  The main program unit of 
!  the atmospheric surface-layer parameterization scheme "SfcFlx".
!  "SfcFlx" is used to compute fluxes 
!  of momentum and of sensible and latent heat over lakes.
!  The surface-layer scheme developed by Mironov (1991) was used as the starting point.
!  It was modified and further developed to incorporate recent results as to 
!  the roughness lenghts for scalar quantities,
!  heat and mass transfer in free convection,
!  and the effect of limited fetch on the momentum transfer.
!  Apart from the momentum flux and sensible and latent heat fluxes,
!  the long-wave radiation flux from the water surface and
!  the long-wave radiation flux from the atmosphere can also be computed.
!  The atmospheric long-wave radiation flux is computed with simple empirical formulae,
!  where the atmospheric emissivity is taken to be dependent on 
!  the water vapour pressure and cloud fraction.
!
!  References:
!  Andreas, E., L., 2002: 
!     Parameterizing scalar trasnfer over ssnow and ice: a review.
!     J. Hydrometeorology, 3, 417-432. 
!  Fung, I. Y., D. E. Harrison, and A. A. Lacis, 1984: 
!     On the variability of the net long-wave radiation at the ocean surface.
!     Rev. Geophys. Space Physics, 22(2), 177-193.
!  Mironov, D. V., 1991:
!     Air-water interaction parameters over lakes.
!     Modelling Air-Lake Interaction. Physical Background,
!     S. S. Zilitinkevich, Ed., Springer-Verlag, Berlin, etc., 50-62.
!  Mironov, D., F. Beyrich, E. Heise, and M. Raschendorfer, 2002:
!     The water surface roughness lengths for temperature: an observational study.
!     COSMO Newsletter No. 2, February 2002, 146-148.
!     (available from the WWW site of the Consortium for Small Scale Modelling,
!     www.cosmo-model.org)
!  Mironov, D., F. Beyrich, E. Heise, and M. Raschendorfer, 2003:
!     The water surface roughness lengths for scalars: an observational study.
!     In preparation.
!  Zapadka, T., and S. B. Wozniak, 2000:
!     Preliminary comparison between various models of the long-wave radiation
!     budget of the sea and experimental data from the Baltic Sea.
!     Oceanologia, 42 (3), 359-369.
!  Zapadka, T., S. B. Wozniak, and B. Wozniak, 2000:
!     A simple formula for the net long-wave radiation flux 
!     in the southern Baltic Sea.
!     Oceanologia, 43 (3), 265-277.
!  Zilitinkevich, S. S., A. A. Grachev, and C. W. Fairall, 2000:
!     Scaling reasoning and field data on the sea-surface roughness lengths for scalars.
!     J. Atmos. Sci., 58, 320-325.
!
!
!  Contact Dmitrii Mironov 
!  German Weather Service, Referat FE14,
!  Frankfurter Str. 135, D-63067 Offenbach am Main, Germany.
!  E-mail: dmitrii.mironov@dwd.de 
!
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! N.NN       YYYY/MM/DD Dmitrii Mironov
!  Initial release
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE data_parameters  , ONLY :   &
  ireals                      , & ! KIND-type parameter for real variables
  iintegers                       ! KIND-type parameter for "normal" integer variables

USE flake_parameters , ONLY :   &
  tpl_grav                    , & ! Acceleration due to gravity [m s^{-2}]
  tpl_T_f                     , & ! Fresh water freezing point [K]
  tpl_rho_w_r                 , & ! Maximum density of fresh water [kg m^{-3}]
  tpl_c_w                     , & ! Specific heat of water [J kg^{-1} K^{-1}]
  tpl_L_f                         ! Latent heat of fusion [J kg^{-1}]

USE flake            , ONLY :   &
  h_Ice_min_flk                   ! Minimum ice thickness [m]

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Dimensionless constants in the Monin-Obukhov surface-layer 
!  similarity relations and in the expressions for the roughness lengths.
REAL (KIND = ireals), PARAMETER ::   &
  c_Karman      = 0.40             , & ! The von Karman constant 
  Pr_neutral    = 1.0              , & ! Turbulent Prandtl number at neutral static stability
  Sc_neutral    = 1.0              , & ! Turbulent Schmidt number at neutral static stability
  c_MO_u_stab   = 5.0              , & ! Constant of the MO theory (wind, stable stratification)
  c_MO_t_stab   = 5.0              , & ! Constant of the MO theory (temperature, stable stratification)
  c_MO_q_stab   = 5.0              , & ! Constant of the MO theory (humidity, stable stratification)
  c_MO_u_conv   = 15.0             , & ! Constant of the MO theory (wind, convection)
  c_MO_t_conv   = 15.0             , & ! Constant of the MO theory (temperature, convection)
  c_MO_q_conv   = 15.0             , & ! Constant of the MO theory (humidity, convection)
  c_MO_u_exp    = 0.25             , & ! Constant of the MO theory (wind, exponent)
  c_MO_t_exp    = 0.5              , & ! Constant of the MO theory (temperature, exponent)
  c_MO_q_exp    = 0.5              , & ! Constant of the MO theory (humidity, exponent)
  z0u_ice_rough = 1.0E-03          , & ! Aerodynamic roughness of the ice surface [m] (rough flow)
  c_z0u_smooth  = 0.1              , & ! Constant in the expression for z0u (smooth flow) 
  c_z0u_rough   = 1.23E-02         , & ! The Charnock constant in the expression for z0u (rough flow)
  c_z0u_rough_L = 1.00E-01         , & ! An increased Charnock constant (used as the upper limit)
  c_z0u_ftch_f  = 0.70             , & ! Factor in the expression for fetch-dependent Charnock parameter
  c_z0u_ftch_ex = 0.333333333      , & ! Exponent in the expression for fetch-dependent Charnock parameter
  c_z0t_rough_1 = 4.0              , & ! Constant in the expression for z0t (factor) 
  c_z0t_rough_2 = 3.2              , & ! Constant in the expression for z0t (factor)
  c_z0t_rough_3 = 0.5              , & ! Constant in the expression for z0t (exponent) 
  c_z0q_rough_1 = 4.0              , & ! Constant in the expression for z0q (factor)
  c_z0q_rough_2 = 4.2              , & ! Constant in the expression for z0q (factor)
  c_z0q_rough_3 = 0.5              , & ! Constant in the expression for z0q (exponent)
  c_z0t_ice_b0s = 1.250            , & ! Constant in the expression for z0t over ice
  c_z0t_ice_b0t = 0.149            , & ! Constant in the expression for z0t over ice
  c_z0t_ice_b1t = -0.550           , & ! Constant in the expression for z0t over ice
  c_z0t_ice_b0r = 0.317            , & ! Constant in the expression for z0t over ice
  c_z0t_ice_b1r = -0.565           , & ! Constant in the expression for z0t over ice
  c_z0t_ice_b2r = -0.183           , & ! Constant in the expression for z0t over ice
  c_z0q_ice_b0s = 1.610            , & ! Constant in the expression for z0q over ice
  c_z0q_ice_b0t = 0.351            , & ! Constant in the expression for z0q over ice
  c_z0q_ice_b1t = -0.628           , & ! Constant in the expression for z0q over ice
  c_z0q_ice_b0r = 0.396            , & ! Constant in the expression for z0q over ice
  c_z0q_ice_b1r = -0.512           , & ! Constant in the expression for z0q over ice
  c_z0q_ice_b2r = -0.180           , & ! Constant in the expression for z0q over ice
  Re_z0s_ice_t  = 2.5              , & ! Threshold value of the surface Reynolds number 
                                       ! used to compute z0t and z0q over ice (Andreas 2002)
  Re_z0u_thresh = 0.1                  ! Threshold value of the roughness Reynolds number 
                                       ! [value from Zilitinkevich, Grachev, and Fairall (200),
                                       ! currently not used] 

!  Dimensionless constants 
REAL (KIND = ireals), PARAMETER ::   &
  c_free_conv   = 0.14                 ! Constant in the expressions for fluxes in free convection

!  Dimensionless constants 
REAL (KIND = ireals), PARAMETER ::   &
  c_lwrad_emis  = 0.99                 ! Surface emissivity with respect to the long-wave radiation

!  Thermodynamic parameters
REAL (KIND = ireals), PARAMETER ::     &
  tpsf_C_StefBoltz    = 5.67E-08     , & ! The Stefan-Boltzmann constant [W m^{-2} K^{-4}]
  tpsf_R_dryair       = 2.8705E+02   , & ! Gas constant for dry air [J kg^{-1} K^{-1}]
  tpsf_R_watvap       = 4.6151E+02   , & ! Gas constant for water vapour [J kg^{-1} K^{-1}]
  tpsf_c_a_p          = 1.005E+03    , & ! Specific heat of air at constant pressure [J kg^{-1} K^{-1}]
  tpsf_L_evap         = 2.501E+06    , & ! Specific heat of evaporation [J kg^{-1}]
  tpsf_nu_u_a         = 1.50E-05     , & ! Kinematic molecular viscosity of air [m^{2} s^{-1}]
  tpsf_kappa_t_a      = 2.20E-05     , & ! Molecular temperature conductivity of air [m^{2} s^{-1}]
  tpsf_kappa_q_a      = 2.40E-05         ! Molecular diffusivity of air for water vapour [m^{2} s^{-1}]

!  Derived thermodynamic parameters
REAL (KIND = ireals), PARAMETER ::                        &
  tpsf_Rd_o_Rv  = tpsf_R_dryair/tpsf_R_watvap           , & ! Ratio of gas constants (Rd/Rv)
  tpsf_alpha_q  = (1._ireals-tpsf_Rd_o_Rv)/tpsf_Rd_o_Rv     ! Diemsnionless ratio 

!  Thermodynamic parameters
REAL (KIND = ireals), PARAMETER ::     &
  P_a_ref             = 1.0E+05          ! Reference pressure [N m^{-2} = kg m^{-1} s^{-2}]


!  The variables declared below
!  are accessible to all program units of the MODULE "SfcFlx"
!  and to the driving routines that use "SfcFlx".
!  These are basically the quantities computed by SfcFlx.
!  Apart from these quantities, there a few local scalars 
!  used by SfcFlx routines mainly for security reasons.
!  All variables declared below have a suffix "sf".

!  SfcFlx variables of type REAL

!  Roughness lengths
REAL (KIND = ireals) ::    &
  z0u_sf                 , & ! Roughness length with respect to wind velocity [m]
  z0t_sf                 , & ! Roughness length with respect to potential temperature [m]
  z0q_sf                     ! Roughness length with respect to specific humidity [m]

!  Fluxes in the surface air layer
REAL (KIND = ireals) ::    &
  u_star_a_sf            , & ! Friction velocity [m s^{-1}]
  Q_mom_a_sf             , & ! Momentum flux [N m^{-2}]
  Q_sens_a_sf            , & ! Sensible heat flux [W m^{-2}]
  Q_lat_a_sf             , & ! Laten heat flux [W m^{-2}]
  Q_watvap_a_sf              ! Flux of water vapout [kg m^{-2} s^{-1}]

!  Security constants
REAL (KIND = ireals), PARAMETER ::   &
  u_wind_min_sf  = 1.0E-02         , & ! Minimum wind speed [m s^{-1}]
  u_star_min_sf  = 1.0E-04         , & ! Minimum value of friction velocity [m s^{-1}]
  c_accur_sf     = 1.0E-07         , & ! A small number (accuracy)
  c_small_sf     = 1.0E-04             ! A small number (used to compute fluxes)

!  Useful constants
REAL (KIND = ireals), PARAMETER ::     &
  num_1o3_sf = 1._ireals/3._ireals       ! 1/3

!==============================================================================
! Procedures 
!==============================================================================

CONTAINS

!==============================================================================
!  The codes of the SfcFlx procedures are stored in separate "*.incf" files
!  and are included below.
!------------------------------------------------------------------------------

include 'src_SfcFlx_lwradatm.incf'

include 'src_SfcFlx_lwradwsfc.incf'

include 'src_SfcFlx_momsenlat.incf'

include 'src_SfcFlx_rhoair.incf'

include 'src_SfcFlx_roughness.incf'

include 'src_SfcFlx_satwvpres.incf'

include 'src_SfcFlx_spechum.incf'

include 'src_SfcFlx_wvpreswetbulb.incf'

!------------------------------------------------------------------------------
!  Procedures used for debugging purposes only
!------------------------------------------------------------------------------

!_nu SUBROUTINE prn_dbg_sf
!_nu  PRINT*, ' '
!_nu END SUBROUTINE prn_dbg_sf

!==============================================================================

END MODULE SfcFlx

