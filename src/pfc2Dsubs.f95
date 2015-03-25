
! fortran_free_source
!
! (c) Copyright 2010, 2015 Bradley J. Eck
! This module is part of PERFCODE 
!

! This module holds external procedures (subroutine and functions)
!  for the pfc2D model (PERFCODE).
! Using module creates an explicit interface for the procedures


!=========================================================================
!     \\\\\\\\\\\  B E G I N                          //////////
                                  MODULE pfc2Dsubs
!     //////////                                      \\\\\\\\\\
!=========================================================================

IMPLICIT NONE

CONTAINS


!===================================================================
SUBROUTINE set_ABCDEF( i, j, Cw1, Ce1, Cs1, Cn1, pf, dt, rr )
!   Fills the arrays of the linear system for Cell v
USE shared,    ONLY: A, B, C, D, E, Fn, F1, F, jmax
USE pfc2Dfuns, ONLY: F_LinearIndex, F_RHS_n1
implicit none
!   Arguments
integer, intent( in ) :: i, j
real   , intent( in ) :: Cw1, Ce1, Cs1, Cn1, pf, dt, rr
!   Internal Variables??
integer :: v
!real, external :: F_RHS_n1
!----------------------------------------------------------
! Linear Index
v   = F_LinearIndex( i, j, jmax)

! Bands of penta-diagonal matrix
A(v) = - dt / 2. * pf * Cw1
B(v) = - dt / 2. * pf * Cs1
C(v) =   dt / 2. * pf * ( Cw1 + Cs1 + Cn1 + Ce1 ) + 1.
D(v) = - dt / 2. * pf * Cn1
E(v) = - dt / 2. * pf * Ce1
! Right-hand-side
! Portion from time level n+1
F1(v) = F_RHS_n1( i, j, Cw1, Ce1, Cs1, Cn1, rr, pf, dt )
!The complete right hand side has contributions from 
! time level n and time level n+1
F(v) = Fn(v) + F1(v)
!---------------------------------------------------------------------
end subroutine set_ABCDEF
!========================================================================
!   \\\\\\\\\\   E N D    S U B R O U T I N E        //////////
!   //////////     S E T _ A B C D E F               \\\\\\\\\\\
!=========================================================================

!========================================================================
!   \\\\\\\\\\   B E G I N    S U B R O U T I N E    //////////
!   //////////     S E T _ X Y H                     \\\\\\\\\\\
!=========================================================================
subroutine SET_xyh( i , j, xx, yy , hh)
!   Assigns values to X, Y, and Z for pointing to the bi-linear
!    interpolatoin subroutine
use shared, only: jmax, h_old, Z, CV_Info
use pfc2dfuns, only: F_LinearIndex
!VARIABLE DECLARATIONS
!   Arguments
integer,intent( in ) :: i, j        ! Grid indices
real, intent( out ) :: xx, yy, hh   ! physical coordinates
!Internal variables
integer :: v
!-----------------------------------------------------------------------
v = F_LinearIndex( i, j, jmax )
xx = CV_Info( v ) % X
yy = CV_Info( v ) % Y
hh = h_old( i, j )

!----------------------------------------------------------------------
end subroutine SET_xyh
!========================================================================
!   \\\\\\\\\\   E N D        S U B R O U T I N E    //////////
!   //////////     S E T _ X Y H                     \\\\\\\\\\\
!=========================================================================


!=========================================================================
!     \\\\\\\\\\\                                     //////////
                         END MODULE pfc2Dsubs
!     //////////                                      \\\\\\\\\\
!=========================================================================













