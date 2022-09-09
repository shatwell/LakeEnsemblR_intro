!subroutine Qaddinflow(Ts, Qinflow)

real (KIND = ireals) function Qinflow(Ts)
! external variables
USE data_parameters , ONLY : &
    ireals,                  & ! KIND-type parameter for real variables
    iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_configure, ONLY: &
    nmlfile
USE flake_parameters, ONLY: &
    tpl_c_w,&
    tpl_rho_w_r,&
    tpl_T_f
!------------------------------
implicit none
!------------------------------
! arguments
REAL (KIND = ireals) :: Ts
! REAL (KIND = ireals), INTENT(OUT) :: Qinflow
!------------------------------
! Internal variables:
REAL (KIND = ireals):: &
    DTinflow, &   ! Tinflow, K or (Tinflow-Ts), K 
                 ! depending on the value of 'DTfixed'
    QS           ! QS - inflow rate ( inflow[m^3/s] / lake surface[m^2] ) m/s
LOGICAL DTfixed  ! 

NAMELIST /inflow/ DTinflow, DTfixed, QS 
!------------------------------

! read variables from namelist
OPEN(UNIT=11, FILE=nmlfile)
!
READ(UNIT=11, NML=inflow)
CLOSE(UNIT=11)

if (.NOT.DTfixed) DTinflow = (DTinflow+tpl_T_f - Ts)
Qinflow = tpl_c_w*tpl_rho_w_r*DTinflow*QS
return
end