! fortran_free_source

module pfc1Dfuns

implicit none

contains



!============================================================================
! function to compute the conveyance coef at the western face
! the convention used here is that 'in' refers to cell 'i'
! and 'out' refers to cell 'i-1', which is the western cell

FUNCTION F_CC( xin,  dxin,  hin, zin, &
              xout, dxout, hout, zout )  Result( CC )
use shared, only: K_input, n_mann_input, b_pfc_input, h_pfc_min          
use pfc1dfuns2     
implicit none
REAL :: xin,  xout   ! coordinate of the ith cell and the WESTERN cell center
REAL :: dxin, dxout  ! cell sizes
REAL :: hin,  hout   ! thicknesses at cell center
REAL :: zin,  zout   ! elevations  at cell center 
REAL :: CC   !, F_hp_face, F_hs_face  <----these now in a module
! dummy vars
REAL :: hpw !thickness in the PAVEMENT at the western face
REAL :: hsw !thickness on the SURFACE at the western face
REAL :: Sfw !magnitude of hydraulic gradient at the western face
logical :: error
!-------------------------------------------


! Intermediate quantities
hpw =  F_hp_face(dxin, hin, zin, dxout, hout, zout, b_pfc_input)
hsw =  F_hs_face(dxin, hin, zin, dxout, hout, zout, b_pfc_input)
Sfw =  sqrt( ( ( hout + zout - hin - zin ) * 2. / &
                (dxout + dxin) ) ** 2 )

! Set hpw to small but positive and with enough range
!   left to allow further calcs.zero if negative
if( hpw .LT. TINY( hpw ) ) then
    hpw   =  h_pfc_min   
end if

!Conveyance coefficient itself
if( hsw .GT. 0.0 ) then
   CC  = 1. / abs( xin - xout ) *                       &
            (   K_input * hpw +                               &
              1./ n_mann_input * hsw ** (5./3.) /  sqrt(Sfw) )
else
! only PFC flow
   CC  = 1. / abs( xin - xout ) *                       &
            (   K_input * hpw  ) 
end if
   

! ERROR CHECKING FOR CONVEYANCE COEFS
if( (CC .GT. HUGE(CC) ) .OR. (CC .LT. -HUGE(CC) ) ) then
    error = .true.
else
    error = .false.
endif

!Output the parts of the calculation if the error is true
if( error .eqv. .true. ) then
    write(*,*) 'Problem with 1D conveyance coefficient!'
    print *, '       K = ', K_input
    print *, '      hp = ', hpw
    print *, '  n_mann = ', n_mann_input
    print *, '      hs = ', hsw
    print *, '      Sf = ', Sfw
    print *, '     xin = ', xin
    print *, '    xout = ', xout 
    print *, '      CC = ', CC
    write(*,*) 'Stopping Program'
    STOP
endif 

END FUNCTION 
!======================================================================    

!======================================================================
!Function to switch the porosity on/off if the
!  water is in/out of the pavement
FUNCTION F_por(h)
USE shared, only: b_pfc_input, por_input
IMPLICIT NONE
REAL h, F_por
if     ( h >= b_pfc_input ) then
                        F_por = 1.
ELSEIF ( h < b_pfc_input )  then
                        F_por = 1./por_input
end if
END function F_por
!======================================================================




!======================================================================
!       \\\\\\\\\\                                       //////////
                      END MODULE pfc1Dfuns
!       //////////                                       \\\\\\\\\\\
!======================================================================
