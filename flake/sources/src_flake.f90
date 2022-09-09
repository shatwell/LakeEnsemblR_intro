! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

MODULE flake

!------------------------------------------------------------------------------
!
! Description:
!
!  The main program unit of the FLake model,  
!  containing most of the FLake procedures.
!  Most FLake variables and local parameters are declared.
!
!  "FLake" (Fresh-water Lake) is a lake model capable of predicting the surface temperature 
!  in lakes of various depth on the time scales from a few hours to a year.
!  The model is based on a two-layer parametric representation of
!  the temperature profile, where the structure of the stratified layer between the
!  upper mixed layer and the basin bottom, the lake thermocline,
!  is described using the concept of self-similarity of the evolving temperature profile.
!  The concept was put forward by Kitaigorodskii and Miropolsky (1970) 
!  to describe the vertical temperature structure of the oceanic seasonal thermocline.
!  It has since been successfully used in geophysical applications.
!  The concept of self-similarity of the evolving temperature profile
!  is also used to describe the vertical structure of the thermally active upper layer 
!  of bottom sediments and of the ice and snow cover.
!
!  The lake model incorporates the heat budget equations
!  for the four layers in question, viz., snow, ice, water and bottom sediments,
!  developed with due regard for the vertically distributed character
!  of the short-wave radiation heating.
!  The entrainment equation that incorporates the Zilitinkevich (1975) spin-up term
!  is used to compute the depth of a convectively-mixed layer. 
!  A relaxation-type equation is used
!  to compute the wind-mixed layer depth in stable and neutral stratification,
!  where a multi-limit formulation for the equilibrium mixed-layer depth
!  proposed by Zilitinkevich and Mironov (1996)
!  accounts for the effects of the earth's rotation, of the surface buoyancy flux
!  and of the static stability in the thermocline.
!  The equations for the mixed-layer depth are developed with due regard for  
!  the volumetric character of the radiation heating.
!  Simple thermodynamic arguments are invoked to develop
!  the evolution equations for the ice thickness and for the snow thickness.
!  The heat flux through the water-bottom sediment interface is computed,
!  using a parameterization proposed by Golosov et al. (1998).
!  The heat flux trough the air-water interface 
!  (or through the air-ice or air-snow interface)
!  is provided by the driving atmospheric model.
!
!  Empirical constants and parameters of the lake model
!  are estimated, using independent empirical and numerical data.
!  They should not be re-evaluated when the model is applied to a particular lake.
!  The only lake-specific parameters are the optical characteristics of lake water,
!  the temperature at the bottom of the thermally active layer
!  of bottom sediments and the depth of this layer.
!
!  A detailed description of the lake model is given in
!  Mironov, D. V., et al., 2003:
!  Parameterization of Lakes in Numerical Weather Prediction.
!  Part 1: Description of a Lake Model.
!  Details of the ice-snow model are given in
!  Mironov, D., and B. Ritter, 2003:
!  A Thermodynamic Model of Ice and Snow over Water Bodies.
!  Manuscripts are available from the authors.
!  Contact Dmitrii Mironov 
!  German Weather Service, Referat FE14,
!  Frankfurter Str. 135, D-63067 Offenbach am Main, Germany.
!  E-mail: dmitrii.mironov@dwd.de 
!
!  Apart from serving as a lake module (parameterisation scheme) 
!  within a three-dimensional NWP or climate model,
!  "FLake" can serve as an independent single-column lake model.
!  The code is easily re-configured by setting a number of logical switches accordingly.
!  This is explained in comment lines where appropriate.
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

USE data_parameters , ONLY : &
  ireals                   , & ! KIND-type parameter for real variables
  iintegers                    ! KIND-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
!
!  The variables declared below
!  are accessible to all program units of the MODULE "flake"
!  and to the driving routines that use "flake".
!  These are basically the quantities computed by FLake.
!  Apart from these quantities, there are a few local scalars 
!  used by FLake routines mainly for security reasons.
!  These are declared as PRIVATE.
!  All variables declared below have a suffix "flk".

!  FLake variables of type REAL

!  Temperatures at the previous time step ("p") and the updated temperatures ("n") 
REAL (KIND = ireals) ::           &
  T_mnw_p_flk, T_mnw_n_flk      , & ! Mean temperature of the water column [K] 
  T_snow_p_flk, T_snow_n_flk    , & ! Temperature at the air-snow interface [K] 
  T_ice_p_flk, T_ice_n_flk      , & ! Temperature at the snow-ice or air-ice interface [K] 
  T_wML_p_flk, T_wML_n_flk      , & ! Mixed-layer temperature [K] 
  T_bot_p_flk, T_bot_n_flk      , & ! Temperature at the water-bottom sediment interface [K] 
  T_B1_p_flk, T_B1_n_flk            ! Temperature at the bottom of the upper layer of the sediments [K] 

!  Thickness of various layers at the previous time step ("p") and the updated values ("n") 
REAL (KIND = ireals) ::           &
  h_snow_p_flk, h_snow_n_flk    , & ! Snow thickness [m]
  h_ice_p_flk, h_ice_n_flk      , & ! Ice thickness [m]
  h_ML_p_flk, h_ML_n_flk        , & ! Thickness of the mixed-layer [m] 
  H_B1_p_flk, H_B1_n_flk            ! Thickness of the upper layer of bottom sediments [m] 

!  The shape factor(s) at the previous time step ("p") and the updated value(s) ("n") 
REAL (KIND = ireals) ::           &
  C_T_p_flk, C_T_n_flk          , & ! Shape factor (thermocline)
  C_TT_flk                      , & ! Dimensionless parameter (thermocline)
  C_Q_flk                       , & ! Shape factor with respect to the heat flux (thermocline)
  C_I_flk                       , & ! Shape factor (ice)
  C_S_flk                           ! Shape factor (snow)

!  Derivatives of the shape functions
REAL (KIND = ireals) ::           &
  Phi_T_pr0_flk                 , & ! d\Phi_T(0)/d\zeta   (thermocline)
  Phi_I_pr0_flk                 , & ! d\Phi_I(0)/d\zeta_I (ice)
  Phi_I_pr1_flk                 , & ! d\Phi_I(1)/d\zeta_I (ice)
  Phi_S_pr0_flk                     ! d\Phi_S(0)/d\zeta_S (snow)

!  Heat and radiation fluxes
REAL (KIND = ireals) ::           &
  Q_snow_flk                    , & ! Heat flux through the air-snow interface [W m^{-2}]
  Q_ice_flk                     , & ! Heat flux through the snow-ice or air-ice interface [W m^{-2}]
  Q_w_flk                       , & ! Heat flux through the ice-water or air-water interface [W m^{-2}]
  Q_bot_flk                     , & ! Heat flux through the water-bottom sediment interface [W m^{-2}]
  I_atm_flk                     , & ! Radiation flux at the lower boundary of the atmosphere [W m^{-2}],
                                    ! i.e. the incident radiation flux with no regard for the surface albedo.
  I_snow_flk                    , & ! Radiation flux through the air-snow interface [W m^{-2}]
  I_ice_flk                     , & ! Radiation flux through the snow-ice or air-ice interface [W m^{-2}]
  I_w_flk                       , & ! Radiation flux through the ice-water or air-water interface [W m^{-2}]
  I_h_flk                       , & ! Radiation flux through the mixed-layer-thermocline interface [W m^{-2}]
  I_bot_flk                     , & ! Radiation flux through the water-bottom sediment interface [W m^{-2}]
  I_intm_0_h_flk                , & ! Integral-mean radiation flux over the mixed layer [W m^{-1}]
  I_intm_h_D_flk                , & ! Integral-mean radiation flux over the thermocline [W m^{-1}]
  Q_star_flk                        ! A generalized heat flux scale [W m^{-2}]

!  Velocity scales
REAL (KIND = ireals) ::           &
  u_star_w_flk                  , & ! Friction velocity in the surface layer of lake water [m s^{-1}]
  w_star_sfc_flk                    ! Convective velocity scale, 
                                    ! using a generalized heat flux scale [m s^{-1}]

!  The rate of snow accumulation
REAL (KIND = ireals) ::           &
  dMsnowdt_flk                      ! The rate of snow accumulation [kg m^{-2} s^{-1}]

!  Security constants
REAL (KIND = ireals), PARAMETER :: &
  h_Snow_min_flk = 1.0E-5        , & ! Minimum snow thickness [m] 
  h_Ice_min_flk  = 1.0E-9        , & ! Minimum ice thickness [m] 
  h_ML_min_flk   = 1.0E-2        , & ! Minimum mixed-layer depth [m] 
  h_ML_max_flk   = 1.0E+3        , & ! Maximum mixed-layer depth [m] 
  H_B1_min_flk   = 1.0E-3        , & ! Minimum thickness of the upper layer of bottom sediments [m] 
  u_star_min_flk = 1.0E-6            ! Minimum value of the surface friction velocity [m s^{-1}]

!  Security constant(s)
REAL (KIND = ireals), PARAMETER, PRIVATE :: &
  c_small_flk    = 1.0E-10                    ! A small number 

!==============================================================================
! Procedures 
!==============================================================================

CONTAINS

!==============================================================================
!  The codes of the FLake procedures are stored in separate "*.incf" files
!  and are included below.
!------------------------------------------------------------------------------
! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

SUBROUTINE flake_driver ( depth_w, depth_bs, T_bs, par_Coriolis,       &
                          extincoef_water_typ,                         &
                          del_time, T_sfc_p, T_sfc_n )         

!------------------------------------------------------------------------------
!
! Description:
!
!  The main driving routine of the FLake model,
!  where computations are performed.
!  Advances the surface temperature
!  and other FLake variables one time step.
!  At the moment, the Euler explicit scheme is used.
!
!  Lines embraced with "!_tmp" contain temporary parts of the code.
!  These should be removed prior to using FLake in applications.
!  Lines embraced/marked with "!_dev" may be replaced
!  as improved parameterizations are developed and tested.
!  Lines embraced/marked with "!_dm" are DM's comments
!  that may be helpful to a user.
!  Lines embraced/marked with "!_dbg" are used 
!  for debugging purposes only.
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
! !VERSION!  !2005/11/29!     Kirillin

! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE flake_parameters           ! Thermodynamic parameters and dimensionless constants of FLake


USE flake_configure            ! Switches that configure FLake


!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations



!  Input (procedure arguments)

REAL (KIND = ireals), INTENT(IN) ::   &
  depth_w                           , & ! The lake depth [m]
  depth_bs                          , & ! Depth of the thermally active layer of the bottom sediments [m]
  T_bs                              , & ! Temperature at the outer edge of 
                                        ! the thermally active layer of the bottom sediments [K]
  par_Coriolis                      , & ! The Coriolis parameter [s^{-1}]
  extincoef_water_typ               , & ! "Typical" extinction coefficient of the lake water [m^{-1}],
                                        ! used to compute the equilibrium CBL depth
  del_time                          , & ! The model time step [s]
  T_sfc_p                               ! Surface temperature at the previous time step [K]  
                                        ! (equal to either T_ice, T_snow or to T_wML)

!  Output (procedure arguments)

REAL (KIND = ireals), INTENT(OUT) ::  &
  T_sfc_n                               ! Updated surface temperature [K] 
                                        ! (equal to the updated value of either T_ice, T_snow or T_wML)


!  Local variables of type LOGICAL
LOGICAL ::          &
  l_ice_create    , & ! Switch, TRUE = ice does not exist but should be created
  l_snow_exists   , & ! Switch, TRUE = there is snow above the ice
  l_ice_meltabove     ! Switch, TRUE = snow/ice melting from above takes place

!  Local variables of type INTEGER
INTEGER (KIND = iintegers) :: &
  i                             ! Loop index

!  Local variables of type REAL
REAL (KIND = ireals) ::    &
  d_T_mnw_dt             , & ! Time derivative of T_mnw [K s^{-1}] 
  d_T_ice_dt             , & ! Time derivative of T_ice [K s^{-1}] 
  d_T_bot_dt             , & ! Time derivative of T_bot [K s^{-1}] 
  d_T_B1_dt              , & ! Time derivative of T_B1 [K s^{-1}] 
  d_h_snow_dt            , & ! Time derivative of h_snow [m s^{-1}]
  d_h_ice_dt             , & ! Time derivative of h_ice [m s^{-1}]
  d_h_ML_dt              , & ! Time derivative of h_ML [m s^{-1}]
  d_H_B1_dt              , & ! Time derivative of H_B1 [m s^{-1}]
  d_C_T_dt                   ! Time derivative of C_T [s^{-1}]

!  Local variables of type REAL
REAL (KIND = ireals) ::    &
  N_T_mean               , & ! The mean buoyancy frequency in the thermocline [s^{-1}] 
  ZM_h_scale             , & ! The ZM96 equilibrium SBL depth scale [m] 
  conv_equil_h_scale         ! The equilibrium CBL depth scale [m]

!  Local variables of type REAL
REAL (KIND = ireals) :: &
  h_ice_threshold     , & ! If h_ice<h_ice_threshold, use quasi-equilibrium ice model 
  flk_str_1           , & ! Help storage variable
  flk_str_2           , & ! Help storage variable
  R_H_icesnow         , & ! Dimensionless ratio, used to store intermediate results
  R_rho_c_icesnow     , & ! Dimensionless ratio, used to store intermediate results
  R_TI_icesnow        , & ! Dimensionless ratio, used to store intermediate results
  R_Tstar_icesnow         ! Dimensionless ratio, used to store intermediate results


!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Compute fluxes, using variables from the previous time step.
!------------------------------------------------------------------------------

!_dm
! At this point, the heat and radiation fluxes, namely,
! Q_snow_flk, Q_ice_flk, Q_w_flk, 
! I_atm_flk, I_snow_flk, I_ice_flk, I_w_flk, I_h_flk, I_bot_flk,     
! the integral-mean radiation flux over the mixed layer, I_intm_0_h_flk, 
! and the integral-mean radiation flux over the thermocline, I_intm_h_D_flk, 
! should be known.
! They are computed within "flake_interface" (or within the driving model)
! and are available to "flake_driver"
! through the above variables declared in the MODULE "flake".
! In case a lake is ice-covered, Q_w_flk is re-computed below.
!_dm

! Heat flux through the ice-water interface
IF(h_ice_p_flk.GE.h_Ice_min_flk) THEN    ! Ice exists 
  IF(h_ML_p_flk.LE.h_ML_min_flk) THEN    ! Mixed-layer depth is zero, compute flux 
    Q_w_flk = -tpl_kappa_w*(T_bot_p_flk-T_wML_p_flk)/depth_w  ! Flux with linear T(z) 
    Phi_T_pr0_flk = Phi_T_pr0_1*C_T_p_flk-Phi_T_pr0_2         ! d\Phi(0)/d\zeta (thermocline)
    Q_w_flk = Q_w_flk*MAX(Phi_T_pr0_flk, 1._ireals)           ! Account for an increased d\Phi(0)/d\zeta 
  ELSE                    
    Q_w_flk = 0._ireals                  ! Mixed-layer depth is greater than zero, set flux to zero
  END IF   
END IF   

! A generalized heat flux scale 
Q_star_flk = Q_w_flk + I_w_flk + I_h_flk - 2._ireals*I_intm_0_h_flk

! Heat flux through the water-bottom sediment interface
IF(lflk_botsed_use) THEN
  Q_bot_flk = -tpl_kappa_w*(T_B1_p_flk-T_bot_p_flk)/MAX(H_B1_p_flk, H_B1_min_flk)*Phi_B1_pr0
ELSE  
  Q_bot_flk = 0._ireals   ! The bottom-sediment scheme is not used
END IF


!------------------------------------------------------------------------------
!  Check if ice exists or should be created.
!  If so, compute the thickness and the temperature of ice and snow.
!------------------------------------------------------------------------------

!_dm
! Notice that a quasi-equilibrium ice-snow model is used 
! to avoid numerical instability when the ice is thin.
! This is always the case when new ice is created.
!_dm

!_dev
! The dependence of snow density and of snow heat conductivity 
! on the snow thickness is accounted for parametrically.
! That is, the time derivatives of \rho_S and \kappa_S are neglected.
! The exception is the equation for the snow thickness 
! in case of snow accumulation and no melting, 
! where d\rho_S/dt is incorporated.
! Furthermore, some (presumably small) correction terms incorporating 
! the snow density and the snow heat conductivity are dropped out.
! Those terms may be included as better formulations 
! for \rho_S and \kappa_S are available.
!_dev

! Default values
l_ice_create    = .FALSE.  
l_ice_meltabove = .FALSE.  

Ice_exist: IF(h_ice_p_flk.LT.h_Ice_min_flk) THEN   ! Ice does not exist 

  l_ice_create = T_wML_p_flk.LE.(tpl_T_f+c_small_flk).AND.Q_w_flk.LT.0._ireals
  IF(l_ice_create) THEN                            ! Ice does not exist but should be created
    d_h_ice_dt = -Q_w_flk/tpl_rho_I/tpl_L_f                                  
    h_ice_n_flk = h_ice_p_flk + d_h_ice_dt*del_time                          ! Advance h_ice 
    T_ice_n_flk = tpl_T_f + h_ice_n_flk*Q_w_flk/tpl_kappa_I/Phi_I_pr0_lin    ! Ice temperature
    d_h_snow_dt = dMsnowdt_flk/tpl_rho_S_min 
    h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time                       ! Advance h_snow
    Phi_I_pr1_flk = Phi_I_pr1_lin                                    & 
                  + Phi_I_ast_MR*MIN(1._ireals, h_ice_n_flk/H_Ice_max)       ! d\Phi_I(1)/d\zeta_I (ice)
    R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/flake_snowheatconduct(h_snow_n_flk) &
                * h_snow_n_flk/MAX(h_ice_n_flk, h_Ice_min_flk)
    T_snow_n_flk = T_ice_n_flk + R_H_icesnow*(T_ice_n_flk-tpl_T_f)           ! Snow temperature
  END IF

ELSE Ice_exist                                     ! Ice exists

  l_snow_exists = h_snow_p_flk.GE.h_Snow_min_flk   ! Check if there is snow above the ice

  Melting: IF(T_snow_p_flk.GE.(tpl_T_f-c_small_flk)) THEN  ! T_sfc = T_f, check for melting from above
                                                           ! T_snow = T_ice if snow is absent 
    IF(l_snow_exists) THEN   ! There is snow above the ice
      flk_str_1 = Q_snow_flk + I_snow_flk - I_ice_flk        ! Atmospheric forcing
      IF(flk_str_1.GE.0._ireals) THEN  ! Melting of snow and ice from above
        l_ice_meltabove = .TRUE.
        d_h_snow_dt = (-flk_str_1/tpl_L_f+dMsnowdt_flk)/flake_snowdensity(h_snow_p_flk)
        d_h_ice_dt  = -(I_ice_flk - I_w_flk - Q_w_flk)/tpl_L_f/tpl_rho_I 
      END IF 
    ELSE                     ! No snow above the ice
      flk_str_1 = Q_ice_flk + I_ice_flk - I_w_flk - Q_w_flk  ! Atmospheric forcing + heating from the water
      IF(flk_str_1.GE.0._ireals) THEN  ! Melting of ice from above, snow accumulation may occur
        l_ice_meltabove = .TRUE.
        d_h_ice_dt  = -flk_str_1/tpl_L_f/tpl_rho_I 
        d_h_snow_dt = dMsnowdt_flk/tpl_rho_S_min
      END IF 
    END IF 
    IF(l_ice_meltabove) THEN  ! Melting from above takes place
      h_ice_n_flk  = h_ice_p_flk  + d_h_ice_dt *del_time  ! Advance h_ice
      h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time  ! Advance h_snow
      T_ice_n_flk  = tpl_T_f                              ! Set T_ice to the freezing point
      T_snow_n_flk = tpl_T_f                              ! Set T_snow to the freezing point
    END IF

  END IF Melting

  No_Melting: IF(.NOT.l_ice_meltabove) THEN                 ! No melting from above

    d_h_snow_dt = flake_snowdensity(h_snow_p_flk)  
    IF(d_h_snow_dt.LT.tpl_rho_S_max) THEN    ! Account for d\rho_S/dt
     flk_str_1 = h_snow_p_flk*tpl_Gamma_rho_S/tpl_rho_w_r
     flk_str_1 = flk_str_1/(1._ireals-flk_str_1)
    ELSE                                     ! Snow density is equal to its maximum value, d\rho_S/dt=0
     flk_str_1 = 0._ireals
    END IF
    d_h_snow_dt = dMsnowdt_flk/d_h_snow_dt/(1._ireals+flk_str_1)       ! Snow accumulation
    h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time                         ! Advance h_snow
    
    Phi_I_pr0_flk = h_ice_p_flk/H_Ice_max                              ! h_ice relative to its maximum value
    C_I_flk = C_I_lin - C_I_MR*(1._ireals+Phi_I_ast_MR)*Phi_I_pr0_flk  ! Shape factor (ice)
    Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*Phi_I_pr0_flk         ! d\Phi_I(1)/d\zeta_I (ice)
    Phi_I_pr0_flk = Phi_I_pr0_lin - Phi_I_pr0_flk                      ! d\Phi_I(0)/d\zeta_I (ice)

    h_ice_threshold = MAX(1._ireals, 2._ireals*C_I_flk*tpl_c_I*(tpl_T_f-T_ice_p_flk)/tpl_L_f)
    h_ice_threshold = Phi_I_pr0_flk/C_I_flk*tpl_kappa_I/tpl_rho_I/tpl_c_I*h_ice_threshold
    h_ice_threshold = SQRT(h_ice_threshold*del_time)                   ! Threshold value of h_ice
    h_ice_threshold = MIN(0.9_ireals*H_Ice_max, MAX(h_ice_threshold, h_Ice_min_flk))
                                                                       ! h_ice(threshold) < 0.9*H_Ice_max

    IF(h_ice_p_flk.LT.h_ice_threshold) THEN  ! Use a quasi-equilibrium ice model

      IF(l_snow_exists) THEN   ! Use fluxes at the air-snow interface
        flk_str_1 = Q_snow_flk + I_snow_flk - I_w_flk
      ELSE                     ! Use fluxes at the air-ice interface
        flk_str_1 = Q_ice_flk + I_ice_flk - I_w_flk
      END IF
      d_h_ice_dt = -(flk_str_1-Q_w_flk)/tpl_L_f/tpl_rho_I
      h_ice_n_flk = h_ice_p_flk + d_h_ice_dt *del_time                         ! Advance h_ice
      T_ice_n_flk = tpl_T_f + h_ice_n_flk*flk_str_1/tpl_kappa_I/Phi_I_pr0_flk  ! Ice temperature

    ELSE                                     ! Use a complete ice model

      d_h_ice_dt = tpl_kappa_I*(tpl_T_f-T_ice_p_flk)/h_ice_p_flk*Phi_I_pr0_flk
      d_h_ice_dt = (Q_w_flk+d_h_ice_dt)/tpl_L_f/tpl_rho_I
      h_ice_n_flk = h_ice_p_flk  + d_h_ice_dt*del_time                         ! Advance h_ice

      R_TI_icesnow = tpl_c_I*(tpl_T_f-T_ice_p_flk)/tpl_L_f         ! Dimensionless parameter
      R_Tstar_icesnow = 1._ireals - C_I_flk                        ! Dimensionless parameter
      IF(l_snow_exists) THEN  ! There is snow above the ice
        R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/flake_snowheatconduct(h_snow_p_flk) &
                    * h_snow_p_flk/h_ice_p_flk
        R_rho_c_icesnow = flake_snowdensity(h_snow_p_flk)*tpl_c_S/tpl_rho_I/tpl_c_I 
!_dev 
!_dm 
! These terms should be included as an improved understanding of the snow scheme is gained, 
! of the effect of snow density in particular. 
!_dm 
!_nu        R_Tstar_icesnow = R_Tstar_icesnow                                                           &
!_nu                        + (1._ireals+C_S_lin*h_snow_p_flk/h_ice_p_flk)*R_H_icesnow*R_rho_c_icesnow
!_dev

        R_Tstar_icesnow = R_Tstar_icesnow*R_TI_icesnow             ! Dimensionless parameter

!_dev
!_nu        R_Tstar_icesnow = R_Tstar_icesnow                                                         &
!_nu                        + (1._ireals-R_rho_c_icesnow)*tpl_c_I*T_ice_p_flk/tpl_L_f
!_dev
        flk_str_2 = Q_snow_flk+I_snow_flk-I_w_flk                  ! Atmospheric fluxes
        flk_str_1  = C_I_flk*h_ice_p_flk + (1._ireals+C_S_lin*R_H_icesnow)*R_rho_c_icesnow*h_snow_p_flk
        d_T_ice_dt = -(1._ireals-2._ireals*C_S_lin)*R_H_icesnow*(tpl_T_f-T_ice_p_flk)             & 
                   * tpl_c_S*dMsnowdt_flk                          ! Effect of snow accumulation
      ELSE                    ! No snow above the ice
        R_Tstar_icesnow = R_Tstar_icesnow*R_TI_icesnow             ! Dimensionless parameter
        flk_str_2 = Q_ice_flk+I_ice_flk-I_w_flk                    ! Atmospheric fluxes
        flk_str_1  = C_I_flk*h_ice_p_flk
        d_T_ice_dt = 0._ireals
      END IF 
      d_T_ice_dt = d_T_ice_dt + tpl_kappa_I*(tpl_T_f-T_ice_p_flk)/h_ice_p_flk*Phi_I_pr0_flk       &
                 * (1._ireals-R_Tstar_icesnow)                     ! Add flux due to heat conduction
      d_T_ice_dt = d_T_ice_dt - R_Tstar_icesnow*Q_w_flk            ! Add flux from water to ice
      d_T_ice_dt = d_T_ice_dt + flk_str_2                          ! Add atmospheric fluxes
      d_T_ice_dt = d_T_ice_dt/tpl_rho_I/tpl_c_I                    ! Total forcing
      d_T_ice_dt = d_T_ice_dt/flk_str_1                            ! dT_ice/dt 
      T_ice_n_flk = T_ice_p_flk + d_T_ice_dt*del_time                          ! Advance T_ice
    END IF

    Phi_I_pr1_flk = MIN(1._ireals, h_ice_n_flk/H_Ice_max)          ! h_ice relative to its maximum value
    Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*Phi_I_pr1_flk     ! d\Phi_I(1)/d\zeta_I (ice)
    R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/flake_snowheatconduct(h_snow_n_flk) &
                 *h_snow_n_flk/MAX(h_ice_n_flk, h_Ice_min_flk)
    T_snow_n_flk = T_ice_n_flk + R_H_icesnow*(T_ice_n_flk-tpl_T_f)             ! Snow temperature

  END IF No_Melting

END IF Ice_exist   

! Security, limit h_ice by its maximum value
h_ice_n_flk = MIN(h_ice_n_flk, H_Ice_max)      

! Security, limit the ice and snow temperatures by the freezing point 
T_snow_n_flk = MIN(T_snow_n_flk, tpl_T_f)  
T_ice_n_flk =  MIN(T_ice_n_flk,  tpl_T_f)    

!_tmp
! Security, avoid too low values (these constraints are used for debugging purposes)
  T_snow_n_flk = MAX(T_snow_n_flk, 73.15_ireals)  
  T_ice_n_flk =  MAX(T_ice_n_flk,  73.15_ireals)    
!_tmp

! Remove too thin ice and/or snow
IF(h_ice_n_flk.LT.h_Ice_min_flk)  THEN        ! Check ice
  h_ice_n_flk = 0._ireals       ! Ice is too thin, remove it, and
  T_ice_n_flk = tpl_T_f         ! set T_ice to the freezing point.
  h_snow_n_flk = 0._ireals      ! Remove snow when there is no ice, and
  T_snow_n_flk = tpl_T_f        ! set T_snow to the freezing point.
  l_ice_create = .FALSE.        ! "Exotic" case, ice has been created but proved to be too thin
ELSE IF(h_snow_n_flk.LT.h_Snow_min_flk) THEN  ! Ice exists, check snow
  h_snow_n_flk = 0._ireals      ! Snow is too thin, remove it, 
  T_snow_n_flk = T_ice_n_flk    ! and set the snow temperature equal to the ice temperature.
END IF


!------------------------------------------------------------------------------
!  Compute the mean temperature of the water column.
!------------------------------------------------------------------------------

IF(l_ice_create) Q_w_flk = 0._ireals     ! Ice has just been created, set Q_w to zero
d_T_mnw_dt = (Q_w_flk - Q_bot_flk + I_w_flk - I_bot_flk)/tpl_rho_w_r/tpl_c_w/depth_w
T_mnw_n_flk = T_mnw_p_flk + d_T_mnw_dt*del_time   ! Advance T_mnw
T_mnw_n_flk = MAX(T_mnw_n_flk, tpl_T_f)           ! Limit T_mnw by the freezing point 


!------------------------------------------------------------------------------
!  Compute the mixed-layer depth, the mixed-layer temperature, 
!  the bottom temperature and the shape factor
!  with respect to the temperature profile in the thermocline. 
!  Different formulations are used, depending on the regime of mixing. 
!------------------------------------------------------------------------------

HTC_Water: IF(h_ice_n_flk.GE.h_Ice_min_flk) THEN    ! Ice exists

  T_mnw_n_flk = MIN(T_mnw_n_flk, tpl_T_r) ! Limit the mean temperature under the ice by T_r 
  T_wML_n_flk = tpl_T_f                   ! The mixed-layer temperature is equal to the freezing point 

  IF(l_ice_create) THEN                  ! Ice has just been created 
    IF(h_ML_p_flk.GE.depth_w-h_ML_min_flk) THEN    ! h_ML=D when ice is created 
      h_ML_n_flk = 0._ireals                 ! Set h_ML to zero 
      C_T_n_flk = C_T_min                    ! Set C_T to its minimum value 
    ELSE                                          ! h_ML<D when ice is created 
      h_ML_n_flk = h_ML_p_flk                ! h_ML remains unchanged 
      C_T_n_flk = C_T_p_flk                  ! C_T (thermocline) remains unchanged 
    END IF 
    T_bot_n_flk = T_wML_n_flk - (T_wML_n_flk-T_mnw_n_flk)/C_T_n_flk/(1._ireals-h_ML_n_flk/depth_w)
                                             ! Update the bottom temperature 

  ELSE IF(T_bot_p_flk.LT.tpl_T_r) THEN   ! Ice exists and T_bot < T_r, molecular heat transfer 
    h_ML_n_flk = h_ML_p_flk                  ! h_ML remains unchanged 
    C_T_n_flk = C_T_p_flk                    ! C_T (thermocline) remains unchanged 
    T_bot_n_flk = T_wML_n_flk - (T_wML_n_flk-T_mnw_n_flk)/C_T_n_flk/(1._ireals-h_ML_n_flk/depth_w)
                                             ! Update the bottom temperature 

  ELSE                                   ! Ice exists and T_bot = T_r, convection due to bottom heating 
    T_bot_n_flk = tpl_T_r                      ! T_bot is equal to the temperature of maximum density 
    IF(h_ML_p_flk.GE.c_small_flk) THEN   ! h_ML > 0 
      C_T_n_flk = C_T_p_flk                     ! C_T (thermocline) remains unchanged 
      h_ML_n_flk = depth_w*(1._ireals-(T_wML_n_flk-T_mnw_n_flk)/(T_wML_n_flk-T_bot_n_flk)/C_T_n_flk)
      h_ML_n_flk = MAX(h_ML_n_flk, 0._ireals)   ! Update the mixed-layer depth  
    ELSE                                 ! h_ML = 0 
      h_ML_n_flk = h_ML_p_flk                   ! h_ML remains unchanged 
      C_T_n_flk = (T_wML_n_flk-T_mnw_n_flk)/(T_wML_n_flk-T_bot_n_flk) 
      C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min)) ! Update the shape factor (thermocline)  
    END IF 
  END IF 

  T_bot_n_flk = MIN(T_bot_n_flk, tpl_T_r)    ! Security, limit the bottom temperature by T_r 

ELSE HTC_Water                                      ! Open water

! Generalised buoyancy flux scale and convective velocity scale
  flk_str_1 = flake_buoypar(T_wML_p_flk)*Q_star_flk/tpl_rho_w_r/tpl_c_w                    
  IF(flk_str_1.LT.0._ireals) THEN       
    w_star_sfc_flk = (-flk_str_1*h_ML_p_flk)**(1._ireals/3._ireals)  ! Convection     
  ELSE 
    w_star_sfc_flk = 0._ireals                                       ! Neutral or stable stratification
  END IF 

!_dm
! The equilibrium depth of the CBL due to surface cooling with the volumetric heating
! is not computed as a solution to the transcendental equation.
! Instead, an algebraic formula is used
! that interpolates between the two asymptotic limits.
!_dm
  conv_equil_h_scale = -Q_w_flk/MAX(I_w_flk, c_small_flk)
  IF(conv_equil_h_scale.GT.0._ireals .AND. conv_equil_h_scale.LT.1._ireals  &
    .AND. T_wML_p_flk.GT.tpl_T_r) THEN   ! The equilibrium CBL depth scale is only used above T_r
    conv_equil_h_scale = SQRT(6._ireals*conv_equil_h_scale)                 &
                       + 2._ireals*conv_equil_h_scale/(1._ireals-conv_equil_h_scale)
    conv_equil_h_scale = MIN(depth_w, conv_equil_h_scale/extincoef_water_typ)
  ELSE
    conv_equil_h_scale = 0._ireals       ! Set the equilibrium CBL depth to zero
  END IF

! Mean buoyancy frequency in the thermocline
  N_T_mean = flake_buoypar(0.5_ireals*(T_wML_p_flk+T_bot_p_flk))*(T_wML_p_flk-T_bot_p_flk)
  IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN
    N_T_mean = SQRT(N_T_mean/(depth_w-h_ML_p_flk))  ! Compute N                   
  ELSE 
    N_T_mean = 0._ireals                            ! h_ML=D, set N to zero
  END IF 

! The rate of change of C_T
  d_C_T_dt = MAX(w_star_sfc_flk, u_star_w_flk, u_star_min_flk)**2_iintegers
  d_C_T_dt = N_T_mean*(depth_w-h_ML_p_flk)**2_iintegers       &
           / c_relax_C/d_C_T_dt                               ! Relaxation time scale for C_T
  d_C_T_dt = (C_T_max-C_T_min)/MAX(d_C_T_dt, c_small_flk)     ! Rate-of-change of C_T 

! Compute the shape factor and the mixed-layer depth, 
! using different formulations for convection and wind mixing

  C_TT_flk = C_TT_1*C_T_p_flk-C_TT_2         ! C_TT, using C_T at the previous time step
  C_Q_flk = 2._ireals*C_TT_flk/C_T_p_flk     ! C_Q using C_T at the previous time step

  Mixing_regime: IF(flk_str_1.LT.0._ireals) THEN  ! Convective mixing 

    C_T_n_flk = C_T_p_flk + d_C_T_dt*del_time                        ! Update C_T, assuming dh_ML/dt>0
    C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min))                ! Limit C_T 
    d_C_T_dt = (C_T_n_flk-C_T_p_flk)/del_time                        ! Re-compute dC_T/dt

    IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN       ! Compute dh_ML/dt
      IF(h_ML_p_flk.LE.h_ML_min_flk) THEN    ! Use a reduced entrainment equation (spin-up)
        d_h_ML_dt = c_cbl_1/c_cbl_2*MAX(w_star_sfc_flk, c_small_flk)

!_dbg
! WRITE*, ' FLake: reduced entrainment eq. D_time*d_h_ML_dt  = ', d_h_ML_dt*del_time
! WRITE*, '         w_*       = ', w_star_sfc_flk
! WRITE*, '         \beta*Q_* = ', flk_str_1
!_dbg

      ELSE                                   ! Use a complete entrainment equation 
        R_H_icesnow     = depth_w/h_ML_p_flk
        R_rho_c_icesnow = R_H_icesnow-1._ireals
        R_TI_icesnow    = C_T_p_flk/C_TT_flk
        R_Tstar_icesnow = (R_TI_icesnow/2._ireals-1._ireals)*R_rho_c_icesnow + 1._ireals
        d_h_ML_dt = -Q_star_flk*(R_Tstar_icesnow*(1._ireals+c_cbl_1)-1._ireals) - Q_bot_flk
        d_h_ML_dt = d_h_ML_dt/tpl_rho_w_r/tpl_c_w                        ! Q_* and Q_b flux terms
        flk_str_2 = (depth_w-h_ML_p_flk)*(T_wML_p_flk-T_bot_p_flk)*C_TT_2/C_TT_flk*d_C_T_dt 
        d_h_ML_dt = d_h_ML_dt + flk_str_2                                 ! Add dC_T/dt term
        flk_str_2 = I_bot_flk + (R_TI_icesnow-1._ireals)*I_h_flk - R_TI_icesnow*I_intm_h_D_flk
        flk_str_2 = flk_str_2 + (R_TI_icesnow-2._ireals)*R_rho_c_icesnow*(I_h_flk-I_intm_0_h_flk)
        flk_str_2 = flk_str_2/tpl_rho_w_r/tpl_c_w
        d_h_ML_dt = d_h_ML_dt + flk_str_2                                 ! Add radiation terms
        flk_str_2 = -c_cbl_2*R_Tstar_icesnow*Q_star_flk/tpl_rho_w_r/tpl_c_w/MAX(w_star_sfc_flk, c_small_flk)
        flk_str_2 = flk_str_2 + C_T_p_flk*(T_wML_p_flk-T_bot_p_flk)
        d_h_ML_dt = d_h_ML_dt/flk_str_2                                   ! dh_ML/dt = r.h.s.
      END IF 
!_dm
! Notice that dh_ML/dt may appear to be negative  
! (e.g. due to buoyancy loss to bottom sediments and/or
! the effect of volumetric radiation heating),
! although a negative generalized buoyancy flux scale indicates 
! that the equilibrium CBL depth has not yet been reached
! and convective deepening of the mixed layer should take place.
! Physically, this situation reflects an approximate character of the lake model.
! Using the self-similar temperature profile in the thermocline, 
! there is always communication between the mixed layer, the thermocline 
! and the lake bottom. As a result, the rate of change of the CBL depth
! is always dependent on the bottom heat flux and the radiation heating of the thermocline.
! In reality, convective mixed-layer deepening may be completely decoupled
! from the processes underneath. In order to account for this fact,
! the rate of CBL deepening is set to a small value
! if dh_ML/dt proves to be negative.
! This is "double insurance" however, 
! as a negative dh_ML/dt is encountered very rarely.
!_dm

!_dbg
! IF(d_h_ML_dt.LT.0._ireals) THEN 
!   WRITE*, 'FLake: negative d_h_ML_dt during convection, = ', d_h_ML_dt
!   WRITE*, '                d_h_ML_dt*del_time = ', MAX(d_h_ML_dt, c_small_flk)*del_time
!   WRITE*, '         u_*       = ', u_star_w_flk   
!   WRITE*, '         w_*       = ', w_star_sfc_flk
!   WRITE*, '         h_CBL_eqi = ', conv_equil_h_scale
!   WRITE*, '         ZM scale  = ', ZM_h_scale
!   WRITE*, '        h_ML_p_flk = ', h_ML_p_flk
! END IF
!   WRITE*, 'FLake: Convection, = ', d_h_ML_dt
!   WRITE*, '         Q_*       = ', Q_star_flk
!   WRITE*, '         \beta*Q_* = ', flk_str_1
!_dbg

      d_h_ML_dt = MAX(d_h_ML_dt, c_small_flk)    
      h_ML_n_flk = h_ML_p_flk + d_h_ML_dt*del_time                       ! Update h_ML 
      h_ML_n_flk = MAX(h_ML_min_flk, MIN(h_ML_n_flk, depth_w))           ! Security, limit h_ML
    ELSE                                              ! Mixing down to the lake bottom
      h_ML_n_flk = depth_w
    END IF

  ELSE Mixing_regime                              ! Wind mixing

    d_h_ML_dt = MAX(u_star_w_flk, u_star_min_flk)                        ! The surface friction velocity
    ZM_h_scale = (ABS(par_Coriolis)/c_sbl_ZM_n + N_T_mean/c_sbl_ZM_i)*d_h_ML_dt**2_iintegers
    ZM_h_scale = ZM_h_scale + flk_str_1/c_sbl_ZM_s
    ZM_h_scale = MAX(ZM_h_scale, c_small_flk)
    ZM_h_scale = d_h_ML_dt**3_iintegers/ZM_h_scale 
    ZM_h_scale = MAX(h_ML_min_flk, MIN(ZM_h_scale, h_ML_max_flk))        ! The ZM96 SBL depth scale 
    ZM_h_scale = MAX(ZM_h_scale, conv_equil_h_scale)                     ! Equilibrium mixed-layer depth 

!_dm 
! In order to avoid numerical discretization problems,
! an analytical solution to the evolution equation 
! for the wind-mixed layer depth is used.
! That is, an exponential relaxation formula is applied
! over the time interval equal to the model time step.
!_dm 

    d_h_ML_dt = c_relax_h*d_h_ML_dt/ZM_h_scale*del_time
    h_ML_n_flk = ZM_h_scale - (ZM_h_scale-h_ML_p_flk)*EXP(-d_h_ML_dt)    ! Update h_ML 
    h_ML_n_flk = MAX(h_ML_min_flk, MIN(h_ML_n_flk, depth_w))             ! Limit h_ML 
    d_h_ML_dt = (h_ML_n_flk-h_ML_p_flk)/del_time                         ! Re-compute dh_ML/dt

    IF(h_ML_n_flk.LE.h_ML_p_flk)           &
      d_C_T_dt = -d_C_T_dt                 ! Mixed-layer retreat or stationary state, dC_T/dt<0
    C_T_n_flk = C_T_p_flk + d_C_T_dt*del_time                            ! Update C_T
    C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min))                    ! Limit C_T 
    d_C_T_dt = (C_T_n_flk-C_T_p_flk)/del_time                            ! Re-compute dC_T/dt

!_dbg
! WRITE*, 'FLake: wind mixing: d_h_ML_dt*del_time = ', d_h_ML_dt*del_time
! WRITE*, '              h_CBL_eqi = ', conv_equil_h_scale
! WRITE*, '              ZM scale  = ', ZM_h_scale
! WRITE*, '              w_*       = ', w_star_sfc_flk
! WRITE*, '              u_*       = ', u_star_w_flk
! WRITE*, '             h_ML_p_flk = ', h_ML_p_flk
!_dbg

  END IF Mixing_regime

! Compute the time-rate-of-change of the the bottom temperature, 
! depending on the sign of dh_ML/dt 
! Update the bottom temperature and the mixed-layer temperature

  IF(h_ML_n_flk.LE.depth_w-h_ML_min_flk) THEN       ! Mixing did not reach the bottom 

    IF(h_ML_n_flk.GT.h_ML_p_flk) THEN   ! Mixed-layer deepening 
      R_H_icesnow     = h_ML_p_flk/depth_w
      R_rho_c_icesnow = 1._ireals-R_H_icesnow 
      R_TI_icesnow    = 0.5_ireals*C_T_p_flk*R_rho_c_icesnow+C_TT_flk*(2._ireals*R_H_icesnow-1._ireals)
      R_Tstar_icesnow = (0.5_ireals+C_TT_flk-C_Q_flk)/R_TI_icesnow
      R_TI_icesnow    = (1._ireals-C_T_p_flk*R_rho_c_icesnow)/R_TI_icesnow
     
      d_T_bot_dt = (Q_w_flk-Q_bot_flk+I_w_flk-I_bot_flk)/tpl_rho_w_r/tpl_c_w
      d_T_bot_dt = d_T_bot_dt - C_T_p_flk*(T_wML_p_flk-T_bot_p_flk)*d_h_ML_dt
      d_T_bot_dt = d_T_bot_dt*R_Tstar_icesnow/depth_w                   ! Q+I fluxes and dh_ML/dt term

      flk_str_2 = I_intm_h_D_flk - (1._ireals-C_Q_flk)*I_h_flk - C_Q_flk*I_bot_flk
      flk_str_2 = flk_str_2*R_TI_icesnow/(depth_w-h_ML_p_flk)/tpl_rho_w_r/tpl_c_w
      d_T_bot_dt = d_T_bot_dt + flk_str_2                               ! Add radiation-flux term

      flk_str_2 = (1._ireals-C_TT_2*R_TI_icesnow)/C_T_p_flk
      flk_str_2 = flk_str_2*(T_wML_p_flk-T_bot_p_flk)*d_C_T_dt
      d_T_bot_dt = d_T_bot_dt + flk_str_2                               ! Add dC_T/dt term
      
    ELSE                                ! Mixed-layer retreat or stationary state
      d_T_bot_dt = 0._ireals                                            ! dT_bot/dt=0
    END IF

    T_bot_n_flk = T_bot_p_flk + d_T_bot_dt*del_time                      ! Update T_bot  
    T_bot_n_flk = MAX(T_bot_n_flk, tpl_T_f)           ! Security, limit T_bot by the freezing point
    flk_str_2 = (T_bot_n_flk-tpl_T_r)*flake_buoypar(T_mnw_n_flk)
    IF(flk_str_2.LT.0._ireals) T_bot_n_flk = tpl_T_r  ! Security, avoid T_r crossover 
    T_wML_n_flk = C_T_n_flk*(1._ireals-h_ML_n_flk/depth_w)
    T_wML_n_flk = (T_mnw_n_flk-T_bot_n_flk*T_wML_n_flk)/(1._ireals-T_wML_n_flk)
    T_wML_n_flk = MAX(T_wML_n_flk, tpl_T_f)           ! Security, limit T_wML by the freezing point

  ELSE                                              ! Mixing down to the lake bottom 

    h_ML_n_flk = depth_w
    T_wML_n_flk = T_mnw_n_flk
    T_bot_n_flk = T_mnw_n_flk
    C_T_n_flk = C_T_min

  END IF

END IF HTC_Water


!------------------------------------------------------------------------------
!  Compute the depth of the thermally active layer of bottom sediments
!  and the temperature at this depth.
!------------------------------------------------------------------------------

Use_sediment: IF(lflk_botsed_use) THEN   ! The bottom-sediment scheme is used
  
  IF(H_B1_p_flk.GE.depth_bs-H_B1_min_flk) THEN   ! No T(z) maximum (no thermal wave) 
    H_B1_p_flk = 0._ireals                       ! Set H_B1_p to zero
    T_B1_p_flk = T_bot_p_flk                     ! Set T_B1_p to the bottom temperature
  END IF 

  flk_str_1 = 2._ireals*Phi_B1_pr0/(1._ireals-C_B1)*tpl_kappa_w/tpl_rho_w_r/tpl_c_w*del_time
  h_ice_threshold = SQRT(flk_str_1)                              ! Threshold value of H_B1
  h_ice_threshold = MIN(0.9_ireals*depth_bs, h_ice_threshold)    ! Limit H_B1
  flk_str_2 = C_B2/(1._ireals-C_B2)*(T_bs-T_B1_p_flk)/(depth_bs-H_B1_p_flk)

  IF(H_B1_p_flk.LT.h_ice_threshold) THEN  ! Use a truncated equation for H_B1(t)
    H_B1_n_flk = SQRT(H_B1_p_flk**2_iintegers+flk_str_1)  ! Advance H_B1
    d_H_B1_dt = (H_B1_n_flk-H_B1_p_flk)/del_time          ! Re-compute dH_B1/dt 
  ELSE                                    ! Use a full equation for H_B1(t)
    flk_str_1 = (Q_bot_flk+I_bot_flk)/H_B1_p_flk/tpl_rho_w_r/tpl_c_w
    flk_str_1 = flk_str_1 - (1._ireals-C_B1)*(T_bot_n_flk-T_bot_p_flk)/del_time
    d_H_B1_dt = (1._ireals-C_B1)*(T_bot_p_flk-T_B1_p_flk)/H_B1_p_flk + C_B1*flk_str_2
    d_H_B1_dt = flk_str_1/d_H_B1_dt
    H_B1_n_flk = H_B1_p_flk + d_H_B1_dt*del_time          ! Advance H_B1
  END IF 
  d_T_B1_dt = flk_str_2*d_H_B1_dt
  T_B1_n_flk = T_B1_p_flk + d_T_B1_dt*del_time            ! Advance T_B1

!_dbg
! WRITE*, 'BS module: '
! WRITE*, '  Q_bot   = ', Q_bot_flk
! WRITE*, '  d_H_B1_dt = ', d_H_B1_dt
! WRITE*, '  d_T_B1_dt = ', d_T_B1_dt
! WRITE*, '  H_B1    = ', H_B1_n_flk
! WRITE*, '    T_bot = ', T_bot_n_flk
! WRITE*, '  T_B1    = ', T_B1_n_flk
! WRITE*, '    T_bs  = ',  T_bs
!_dbg

!_dev 
! Use a very simplistic procedure, where only the upper layer profile is used, 
! H_B1 is always set to depth_bs, and T_B1 is always set to T_bs.
! Then, the time derivatives are zero, and the sign of the bottom heat flux depends on 
! whether T_bot is smaller or greater than T_bs.
! This is, of course, an oversimplified scheme.
!_nu  d_H_B1_dt = 0._ireals
!_nu  d_T_B1_dt = 0._ireals
!_nu  H_B1_n_flk = H_B1_p_flk + d_H_B1_dt*del_time   ! Advance H_B1
!_nu  T_B1_n_flk = T_B1_p_flk + d_T_B1_dt*del_time   ! Advance T_B1
!_dev 

  l_snow_exists = H_B1_n_flk.GE.depth_bs-H_B1_min_flk                    & ! H_B1 reached depth_bs, or
             .OR. H_B1_n_flk.LT.H_B1_min_flk                             & ! H_B1 decreased to zero, or
             .OR.(T_bot_n_flk-T_B1_n_flk)*(T_bs-T_B1_n_flk).LE.0._ireals   ! there is no T(z) maximum
  IF(l_snow_exists) THEN      
    H_B1_n_flk = depth_bs                     ! Set H_B1 to the depth of the thermally active layer
    T_B1_n_flk = T_bs                         ! Set T_B1 to the climatological temperature 
  END IF

ELSE Use_sediment                        ! The bottom-sediment scheme is not used

  H_B1_n_flk = rflk_depth_bs_ref              ! H_B1 is set to a reference value 
  T_B1_n_flk = tpl_T_r                        ! T_B1 is set to the temperature of maximum density

END IF Use_sediment


!------------------------------------------------------------------------------
!  Impose additional constraints.
!------------------------------------------------------------------------------

! In case of unstable stratification, force mixing down to the bottom
flk_str_2 = (T_wML_n_flk-T_bot_n_flk)*flake_buoypar(T_mnw_n_flk)
IF(flk_str_2.LT.0._ireals) THEN 

!_dbg
! WRITE*, 'FLake: inverse (unstable) stratification !!! '
! WRITE*, '       Mixing down to the bottom is forced.'
! WRITE*, '  T_wML_p, T_wML_n ', T_wML_p_flk-tpl_T_f, T_wML_n_flk-tpl_T_f
! WRITE*, '  T_mnw_p, T_mnw_n ', T_mnw_p_flk-tpl_T_f, T_mnw_n_flk-tpl_T_f
! WRITE*, '  T_bot_p, T_bot_n ', T_bot_p_flk-tpl_T_f, T_bot_n_flk-tpl_T_f
! WRITE*, '  h_ML_p,  h_ML_n  ', h_ML_p_flk,          h_ML_n_flk
! WRITE*, '  C_T_p,   C_T_n   ', C_T_p_flk,           C_T_n_flk
!_dbg

  h_ML_n_flk = depth_w
  T_wML_n_flk = T_mnw_n_flk
  T_bot_n_flk = T_mnw_n_flk
  C_T_n_flk = C_T_min

END IF


!------------------------------------------------------------------------------
!  Update the surface temperature.
!------------------------------------------------------------------------------

IF(h_snow_n_flk.GE.h_Snow_min_flk) THEN   
  T_sfc_n = T_snow_n_flk                   ! Snow exists, use the snow temperature
ELSE IF(h_ice_n_flk.GE.h_Ice_min_flk) THEN
  T_sfc_n = T_ice_n_flk                    ! Ice exists but there is no snow, use the ice temperature
ELSE 
  T_sfc_n = T_wML_n_flk                    ! No ice-snow cover, use the mixed-layer temperature
END IF

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END SUBROUTINE flake_driver


include 'src_flake_buoypar.incf'

include 'src_flake_snowdensity.incf'

include 'src_flake_snowheatconduct.incf'

!------------------------------------------------------------------------------
!  Procedures used for debugging purposes only
!------------------------------------------------------------------------------

!_nu SUBROUTINE prn_dbg
!_nu  PRINT*, ' '
!_nu END SUBROUTINE prn_dbg

!==============================================================================

END MODULE flake

