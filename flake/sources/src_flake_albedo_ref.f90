! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

MODULE flake_albedo_ref

!------------------------------------------------------------------------------
!
! Description:
!
!  This module contains "reference" values of the albedo for 
!  the lake water, lake ice and snow. 
!  As in "flake_paramoptic_ref", two ice categories, viz. white ice and blue ice,
!  and two snow categories, viz. dry snow and melting snow, are used.  
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

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
!  Set the albedo for water, ice and snow.
REAL (KIND = ireals), PARAMETER :: &
  albedo_water_ref       = 0.07  , & ! Water
  albedo_whiteice_ref    = 0.60  , & ! White ice
  albedo_blueice_ref     = 0.10  , & ! Blue ice
  albedo_drysnow_ref     = 0.60  , & ! Dry snow 
  albedo_meltingsnow_ref = 0.10      ! Melting snow 

!==============================================================================
!  Empirical parameters.
REAL (KIND = ireals), PARAMETER :: &
  c_albice_MR = 95.6                 ! Constant in the interpolation formula for
                                     ! the ice albedo (Mironov and Ritter 2004)


END MODULE flake_albedo_ref
