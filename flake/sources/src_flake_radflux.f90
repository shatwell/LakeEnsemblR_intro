! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

SUBROUTINE flake_radflux ( depth_w, albedo_water, albedo_ice, albedo_snow, & 
                           opticpar_water, opticpar_ice, opticpar_snow )       

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the radiation fluxes 
!  at the snow-ice, ice-water, air-water, 
!  mixed layer-thermocline and water column-bottom sediment interfaces,
!  the integral-mean radiation flux over the mixed layer,
!  and the integral-mean radiation flux over the thermocline.
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
    ireals,                  & ! KIND-type parameter for real variables
    iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_derivedtypes         ! Definitions of several derived TYPEs

USE flake           , ONLY : & 
  h_Snow_min_flk           , & ! Minimum snow thickness [m]
  h_Ice_min_flk            , & ! Minimum ice thickness [m]
  h_ML_min_flk             , & ! Minimum mixed-layer depth [m]
  h_snow_p_flk             , & ! Snow thickness [m]
  h_ice_p_flk              , & ! Ice thickness [m]
  h_ML_p_flk               , & ! Thickness of the mixed-layer [m]
                               ! All layer depths are from the previous time step ("p").
  I_atm_flk                , & ! Radiation flux at the lower boundary of the atmosphere [W m^{-2}],
                               ! i.e. the incident radiation flux with no regard for the surface albedo.
  I_snow_flk               , & ! Radiation flux through the air-snow interface [W m^{-2}]
  I_ice_flk                , & ! Radiation flux through the snow-ice or air-ice interface [W m^{-2}]
  I_w_flk                  , & ! Radiation flux through the ice-water or air-water interface [W m^{-2}
  I_h_flk                  , & ! Radiation flux through the mixed-layer-thermocline interface [W m^{-2}]
  I_bot_flk                , & ! Radiation flux through the water-bottom sediment interface [W m^{-2}]
  I_intm_0_h_flk           , & ! Integral-mean radiation flux over the mixed layer [W m^{-1}]
  I_intm_h_D_flk               ! Integral-mean radiation flux over the thermocline [W m^{-1}]

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Input (procedure arguments)

REAL (KIND = ireals), INTENT(IN) ::   &
  depth_w                           , & ! The lake depth [m]
  albedo_water                      , & ! Albedo of the water surface 
  albedo_ice                        , & ! Albedo of the water ice 
  albedo_snow                           ! Albedo of the water snow 

TYPE (opticpar_medium), INTENT(IN) :: & 
  opticpar_water                    , & ! Optical characteristics of water
  opticpar_ice                      , & ! Optical characteristics of ice
  opticpar_snow                         ! Optical characteristics of snow 

INTEGER (KIND = iintegers) :: & ! Help variable(s)
  i                             ! DO loop index

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

  IF(h_ice_p_flk.GE.h_Ice_min_flk) THEN            ! Ice exists
    IF(h_snow_p_flk.GE.h_Snow_min_flk) THEN        ! There is snow above the ice
      I_snow_flk = I_atm_flk*(1._ireals-albedo_snow) 
      I_bot_flk = 0._ireals
      DO i=1_iintegers, opticpar_snow%nband_optic
        I_bot_flk = I_bot_flk +                    & 
        opticpar_snow%frac_optic(i)*EXP(-opticpar_snow%extincoef_optic(i)*h_snow_p_flk) 
      END DO 
      I_ice_flk  = I_snow_flk*I_bot_flk
    ELSE                                           ! No snow above the ice 
      I_snow_flk = I_atm_flk  
      I_ice_flk  = I_atm_flk*(1._ireals-albedo_ice)
    END IF 
    I_bot_flk = 0._ireals
    DO i=1_iintegers, opticpar_ice%nband_optic
      I_bot_flk = I_bot_flk +                      & 
      opticpar_ice%frac_optic(i)*EXP(-opticpar_ice%extincoef_optic(i)*h_ice_p_flk) 
    END DO 
    I_w_flk      = I_ice_flk*I_bot_flk
  ELSE                                             ! No ice-snow cover
    I_snow_flk   = I_atm_flk  
    I_ice_flk    = I_atm_flk
    I_w_flk      = I_atm_flk*(1._ireals-albedo_water)
  END IF 

  IF(h_ML_p_flk.GE.h_ML_min_flk) THEN           ! Radiation flux at the bottom of the mixed layer
    I_bot_flk = 0._ireals
    DO i=1_iintegers, opticpar_water%nband_optic
      I_bot_flk = I_bot_flk +            & 
      opticpar_water%frac_optic(i)*EXP(-opticpar_water%extincoef_optic(i)*h_ML_p_flk) 
    END DO 
    I_h_flk = I_w_flk*I_bot_flk
  ELSE                                          ! Mixed-layer depth is less then a minimum value
    I_h_flk = I_w_flk
  END IF

  I_bot_flk = 0._ireals                         ! Radiation flux at the lake bottom
  DO i=1_iintegers, opticpar_water%nband_optic
    I_bot_flk = I_bot_flk +              & 
    opticpar_water%frac_optic(i)*EXP(-opticpar_water%extincoef_optic(i)*depth_w) 
  END DO 
  I_bot_flk = I_w_flk*I_bot_flk

  IF(h_ML_p_flk.GE.h_ML_min_flk) THEN           ! Integral-mean radiation flux over the mixed layer
    I_intm_0_h_flk = 0._ireals
    DO i=1_iintegers, opticpar_water%nband_optic
      I_intm_0_h_flk = I_intm_0_h_flk +                                &
      opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*  &
      (1._ireals - EXP(-opticpar_water%extincoef_optic(i)*h_ML_p_flk))
    END DO 
    I_intm_0_h_flk = I_w_flk*I_intm_0_h_flk/h_ML_p_flk
  ELSE
    I_intm_0_h_flk = I_h_flk
  END IF

  IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN   ! Integral-mean radiation flux over the thermocline
    I_intm_h_D_flk = 0._ireals 
    DO i=1_iintegers, opticpar_water%nband_optic
      I_intm_h_D_flk = I_intm_h_D_flk +                                &
      opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*  &
      ( EXP(-opticpar_water%extincoef_optic(i)*h_ML_p_flk)             &
      - EXP(-opticpar_water%extincoef_optic(i)*depth_w) )
    END DO 
    I_intm_h_D_flk = I_w_flk*I_intm_h_D_flk/(depth_w-h_ML_p_flk)
  ELSE
    I_intm_h_D_flk = I_h_flk
  END IF

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END SUBROUTINE flake_radflux

