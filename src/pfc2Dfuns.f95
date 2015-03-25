! fortran_free_source

! This module holds external procedures (subroutine and functions)
!  for the pfc2D model (PERFCODE).
! Using module creates an explicit interface for the procedures


module pfc2Dfuns

implicit none

contains

!   1. F_LinearIndex
!   2. F_por
!   3. F_RHS_n
!   4. F_RHS_n1


!=============================================================
Function F_LinearIndex( i, j, jmax )
!   Converts grid index to one-dimensional storage location
implicit none
integer, intent( in ) :: i, j, jmax
integer               :: F_LinearIndex
F_LinearIndex = ( i - 1) * jmax + j
end Function F_LinearIndex
!==============================================================


!===================================================================================
!Function to switch the porosity on/off if the water is in/out of the pavement
FUNCTION F_por(h, b_pfc, por)
IMPLICIT NONE
REAL h, F_por, b_pfc, por
if     ( h >= b_pfc ) then
            F_por = 1.
ELSEIF ( h < b_pfc  )  then
            F_por = 1./por
end if
END Function F_por
!====================================================================================

Function F_RHS_n( i, j, Cw, Ce, Cs, Cn, rr, pf, dt ) Result( Fn )
!   Computes the RHS of the linear system for time level n
use shared, only: h_old, Z, imax, jmax
implicit none
! Arguments
integer, intent( in ) :: i, j
real   , intent( in ) :: Cw, Ce, Cs, Cn, rr, pf, dt
! Internal variables
!   added a bunch of dummy variables with if statements to have this function
!   also work at the boundaries.
real :: Fn 
real :: hw, he, hn, hs, Zw, Ze, Zs, Zn


! Thicknesses
if( i == 1 )   then; hw = 0.0; else; hw = h_old(i-1,j); endif
if( j == 1 )   then; hs = 0.0; else; hs = h_old(i,j-1); endif
if( j == jmax) then; hn = 0.0; else; hn = h_old(i,j+1); endif
if( i == imax) then; he = 0.0; else; he = h_old(i+1,j); endif
! Elevations
if( i == 1 )   then; Zw = 0.0; else; Zw = Z(i-1,j); endif
if( j == 1 )   then; Zs = 0.0; else; Zs = Z(i,j-1); endif
if( j == jmax) then; Zn = 0.0; else; Zn = Z(i,j+1); endif
if( i == imax) then; Ze = 0.0; else; Ze = Z(i+1,j); endif

!Compute the RHS from time level n
    Fn    = h_old(i,j)   +                        &
            pf * dt / 2. * (  Cw * hw  +  Cs * hs &
                            + Cn * hn  +  Ce * he &
                            + Cw * Zw  +  Cs * Zs &
                            + Cn * Zn  +  Ce * Ze &
                            -  (Cw + Cs + Cn + Ce) * h_old(i,j)       &  
                            -  (Cw + Cs + Cn + Ce) *     Z(i,j)       &
                            + rr                                        )
end function F_RHS_n
!====================================================================================


Function F_RHS_n1( i, j, Cw1, Ce1, Cs1, Cn1, rr, pf, dt ) Result (F1)
!   Computes the part of the RHS due to time level n+1
use shared, only: Z, imax, jmax
implicit none
!   Arguments
integer, intent( in ) :: i,j
real   , intent( in ) :: Cw1, Ce1, Cs1, Cn1, rr, pf, dt
!   Internal Variables
real :: F1
real :: Zw, Ze, Zs, Zn

! Elevations
if( i == 1 )   then; Zw = 0.0; else; Zw = Z(i-1,j); endif
if( j == 1 )   then; Zs = 0.0; else; Zs = Z(i,j-1); endif
if( j == jmax) then; Zn = 0.0; else; Zn = Z(i,j+1); endif
if( i == imax) then; Ze = 0.0; else; Ze = Z(i+1,j); endif

F1        = pf * dt / 2. * (  Cw1 * Zw  +  Cs1 * Zs &
                            + Cn1 * Zn  +  Ce1 * Ze &
                            - (Cw1 + Cs1 + Cn1 + Ce1)* Z(i,j)   &  
                            + rr                                   )
end function F_RHS_n1
!=================================================================================



end module pfc2Dfuns
