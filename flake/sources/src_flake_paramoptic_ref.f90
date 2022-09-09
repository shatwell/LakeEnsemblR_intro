! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

MODULE flake_paramoptic_ref

!------------------------------------------------------------------------------
!
! Description:
!
!  This module contains "reference" values of the optical characteristics
!  of the lake water, lake ice and snow. These reference values may be used 
!  if no information about the optical characteristics of the lake in question 
!  is available. An exponential decay law for the radiation flux is assumed,
!  The extinction coefficient for the water is set to a large value,
!  leading to the absorption of 95% of the incoming radiation 
!  within the uppermost 1 m of the lake water. 
!  The extinction coefficients for ice and snow are taken from 
!  Launiainen and Cheng (1998). The estimates for the ice correspond 
!  to the uppermost 0.1 m of the ice layer and to the clear sky conditions 
!  (see Table 2 in op. cit.).
!  Very large values of the extinction coefficients for ice and snow ("opaque")
!  can be used to prevent penetration of the solar radiation 
!  through the snow-ice cover.
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

USE data_parameters, ONLY :      &
  ireals                       , & ! KIND-type parameter for real variables 
  iintegers                        ! KIND-type parameter for "normal" integer variables

USE flake_derivedtypes, ONLY :   &
  nband_optic_max              , & ! Maximum value of the wave-length bands
  opticpar_medium                  ! Derived TYPE 

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
INTEGER (KIND = iintegers), PRIVATE :: & ! Help variable(s)
  i                                      ! DO loop index

!  Set optical characteristics for water, ice and snow.
!  The simplest one-band approximation is assumed as a reference.
TYPE (opticpar_medium), PARAMETER ::                       & 
  opticpar_water_ref = opticpar_medium(1,                  & ! Water
    (/1., (0.,i=2_iintegers,nband_optic_max)/),            &
    (/3., (1.E+10,i=2_iintegers,nband_optic_max)/))      , &
  opticpar_water_trans = opticpar_medium(2,                  & ! Transparent Water
    (/0.10, 0.90, (0.,i=3_iintegers,nband_optic_max)/),      &
    (/2.0, 0.20, (1.E+10,i=3_iintegers,nband_optic_max)/)) , &
!_nu  opticpar_water_trans = opticpar_medium(1,                & ! Transparent Water
!_nu    (/1., (0.,i=2_iintegers,nband_optic_max)/),            &
!_nu    (/0.30, (1.E+10,i=2_iintegers,nband_optic_max)/))    , &
  opticpar_whiteice_ref = opticpar_medium(1,               & ! White ice
    (/1., (0.,i=2_iintegers,nband_optic_max)/),            &   
    (/17.1, (1.E+10,i=2_iintegers,nband_optic_max)/))    , &
  opticpar_blueice_ref = opticpar_medium(1,                & ! Blue ice
    (/1., (0.,i=2_iintegers,nband_optic_max)/),            &
    (/8.4, (1.E+10,i=2_iintegers,nband_optic_max)/))     , &
  opticpar_drysnow_ref = opticpar_medium(1,                & ! Dry snow 
    (/1., (0.,i=2_iintegers,nband_optic_max)/),            &
    (/25.0, (1.E+10,i=2_iintegers,nband_optic_max)/))    , &
  opticpar_meltingsnow_ref = opticpar_medium(1,            & ! Melting snow 
    (/1., (0.,i=2_iintegers,nband_optic_max)/),            &
    (/15.0, (1.E+10,i=2_iintegers,nband_optic_max)/))    , &
  opticpar_ice_opaque = opticpar_medium(1,                 & ! Opaque ice
    (/1., (0.,i=2_iintegers,nband_optic_max)/),            &
    (/1.0E+07, (1.E+10,i=2_iintegers,nband_optic_max)/)) , &
  opticpar_snow_opaque = opticpar_medium(1,                & ! Opaque snow
    (/1., (0.,i=2_iintegers,nband_optic_max)/),            &
    (/1.0E+07, (1.E+10,i=2_iintegers,nband_optic_max)/)) 

!==============================================================================

END MODULE flake_paramoptic_ref
