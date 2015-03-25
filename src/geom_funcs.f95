! fortran_free_source
!
! (c) Copyright 2010, 2015 Bradley J. Eck
! This module is part of PERFCODE 
!

!============================================================================
!   \\\\\\\\\\                                        //////////
                        MODULE geom_funcs
!   //////////                                        \\\\\\\\\\
                        implicit none

                        contains
!============================================================================
!   1.  F_L_xi
!   2.  unmap_x
!   3.  unmap_y

!===========================================================================
Function F_L_xi(xi, eta, seg) Result(L_xi) !xcc1, ycc1, dx, dy, R1, dR, W, theta1, dtheta)
!   Computes the METRIC COEFFICIENT for the length mapping.
!Function F_length_xi(xi, eta, xcc1, ycc1, dx, dy, R1, dR, W, theta1, dtheta)
! GEOMETRY MAPPING FUNCTIONS from Geometry.xlsb
use shared, only: CLSEG
implicit none
!   Arguments
real xi, eta  
type(CLSEG) :: seg
!   Result
real L_xi
!   Internal Variables
real angle, dx_dxi, dy_dxi
real xcc1, ycc1, dx, dy, R1, dR, W, theta1, dtheta
!------------------------------------------------------------------------
! Assign parts of the derived type to local variables
!  to keep the formulas cleaner
xcc1   = seg%xcc1
ycc1   = seg%ycc1
dx     = seg%dx
dy     = seg%dy
R1     = seg%R1
dR     = seg%dR
W      = seg%W
theta1 = seg%theta1
dtheta = seg%dtheta

! compute intermediate variables
Angle = theta1 + xi * dtheta

dx_dxi = dx + dR * Cos(Angle) - dtheta * Sin(Angle) * &
         (R1 + W * (eta - 0.5) + xi * dR)

dy_dxi = dx + dR * Sin(Angle) + dtheta * Cos(Angle) * &
         (R1 + W * (eta - 0.5) + xi * dR)
         
! Calculate metric coefficient
L_xi = sqrt( (dx_dxi ** 2 + dy_dxi ** 2) )

End Function F_L_xi
!=============================================================================

Function unmap_x(xi, eta, seg) Result( X )
!   
use shared, only: CLSEG
implicit none
!   Arguments
real xi, eta  
type(CLSEG) :: seg
!   Result
real X
!   Internal Variables
real xcc1, dx, R1, dR, W, theta1, dtheta
!------------------------------------------------------------------------
! Assign parts of the derived type to local variables
!  to keep the formulas cleaner
xcc1   = seg%xcc1
dx     = seg%dx
R1     = seg%R1
dR     = seg%dR
W      = seg%W
theta1 = seg%theta1
dtheta = seg%dtheta

! Compute the X coordinate
X = (xcc1 + xi * dx) + &
          (R1 + xi * dR + (eta - 0.5) * W) * Cos(theta1 + xi * dtheta)

end function unmap_x
!===============================================================================


Function unmap_y(xi, eta, seg) result( Y )

use shared, only: CLSEG
implicit none
real xi, eta  
type(CLSEG) :: seg
!   Result
real Y
!   Internal Variables
real ycc1, dy, R1, dR, W, theta1, dtheta
!------------------------------------------------------------------------
! Assign parts of the derived type to local variables
!  to keep the formulas cleaner
ycc1   = seg%ycc1
dy     = seg%dy
R1     = seg%R1
dR     = seg%dR
W      = seg%W
theta1 = seg%theta1
dtheta = seg%dtheta


Y = (ycc1 + xi * dy) + &
          (R1 + xi * dR + (eta - 0.5) * W) * Sin(theta1 + xi * dtheta)

end function unmap_y
!==============================================================================









!============================================================================
!   \\\\\\\\\                                       ///////////
                       END MODULE geom_funcs
!   /////////                                       \\\\\\\\\\\\
!============================================================================
