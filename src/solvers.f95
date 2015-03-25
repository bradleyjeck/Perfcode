! fortran_free_source
!
! (c) Copyright 2010, 2015 Bradley J. Eck
! This module is part of PERFCODE 
!
!
!  This module contains subroutines for a few linear solvers
!============================================================================
!   \\\\\\\\\\                                        //////////
                        MODULE solvers 
!   //////////                                        \\\\\\\\\\
                        implicit none

                        contains
!============================================================================
!   Subroutines related to solving linear systems:
!   1. DIAGDOM_PENTA checks for diagonal dominance
!        given the bands of a penta-diagonal matrix
!   2. GAUSS_SEIDEL_PENTA uses the Gauss-Seidel method
!       for iterative solution of a penta-diagonal system 
!       of linear equations. 
!   3. THOMAS uses the tri-diagonal matrix algorithm to solve
!        a tri-diagonal linear system 



!===================================================================
!   \\\\\\\\\\       B E G I N    S U B R O U T I N E ///////////
!   //////////          D I A G D O M _ P E N T A     \\\\\\\\\\\
!===================================================================
!
!   Purpose:    Checks to see if a penta-diagonal matrix is
!               diagonally dominant.  Knowing this helps select
!               a solver.  The routine operates only on the bands
!               of the coefficent matrix.
subroutine diagdom_penta( A, B, C, D, E, n, LB, UB, diagdom)
!  A,B -- lower bands of the penta-diagonal matrix
!   C  -- main diagonal
!  D,E -- upper banks of the penta-diagonal matrix
!   n  -- number of unknowns (size of system)
!  LB  -- lower bandwidth
!  UB  -- upper bandwidth
! tolit-- iteration tolerence
! diagdom-- a logical that stores the result.
!-------------------------------------------------------------------
! VARIABLE DECLARATIONS

! Arguments
integer, intent(in)  :: n, LB, UB
real, intent( in )   :: A(n), B(n), C(n), D(n), E(n)
logical, intent(out) :: diagdom
! Internal Variables
integer              :: k
real                 :: T1, T2, T4, T5
real                 :: tot


!--------------------------------------------------------------

! set the logical to true, the following loop changes it if 
! a row is not diagonally dominant.
diagdom = .true.

do k = 1, n
    ! compute the magnitude of each term in the row of the matrix
    if( k-LB .LT. 1 ) then; T1 = 0. ; else; T1 = abs( A(k) ); endif
    if( k    .EQ. 1 ) then; T2 = 0. ; else; T2 = abs( B(k) ); endif
    if( k    .EQ. n ) then; T4 = 0. ; else; T4 = abs( D(k) ); endif
    if( k+UB .GT. n ) then; T5 = 0. ; else; T5 = abs( E(k) ); endif
    ! Test for diagonal dominance
    tot = T1 + T2 + T4 + T5
    if( tot .GT. abs( C(k) ) ) then
        write(*,*) 'Row ', k, 'of the matrix is not diagonally dominant'
        diagdom = .false.
    endif
enddo
!----------------------------------------------------------------
end subroutine diagdom_penta
!===================================================================
!   \\\\\\\\\\       E N D        S U B R O U T I N E ///////////
!   //////////          D I A G D O M _ P E N T A     \\\\\\\\\\\
!===================================================================

!===================================================================
!   \\\\\\\\\\       B E G I N    S U B R O U T I N E ///////////
!   //////////   G A U S S _ S E I D E L _ P E N T A  \\\\\\\\\\\
!===================================================================

! CAUTION -- This routine DOES NOT check convergence criteria
!            so it possible to converge to the wrong answer.


subroutine gauss_seidel_penta( A, B, C, D, E, F, n , LB, UB, &
                               tolit, maxit, Xold, Xnew, dev, numits )
!  A,B -- lower bands of the penta-diagonal matrix
!   C  -- main diagonal
!  D,E -- upper banks of the penta-diagonal matrix
!   F  -- right hand side (force vector) of linear system
!   n  -- number of unknowns (size of system)
!  LB  -- lower bandwidth
!  UB  -- upper bandwidth
! tolit-- iteration tolerence
! maxit-- maximum number of iterations allowed 
! Xold -- initial guess
! Xnew -- converged solution
! dev  -- device for outputting information from the solver
!numits-- number of iterations required to converge 
!----------------------------------------------------------------------
!DECLARATIONS
use utilities, only: F_L2_NORM   
!arguments
integer, intent(IN) :: n, LB, UB
real, intent(in )   :: A(n), B(n), C(n), D(n), E(n), F(n)
real, intent(in )   :: tolit
integer, intent(in) :: maxit
real, intent(in)    :: Xold(n) 
real, intent(out)   :: Xnew(n)
integer, intent(in) :: dev
integer, intent(out):: numits  ! number of iterations required
! internal variables
integer :: k    ! array index
integer :: m    ! iteration index
real :: T1, T2, T4, T5  ! Terms in the equation 
real :: relchng(n)  !relative change between iterations
real :: Xtmp(n)     !temporary array to store the progressive solutions

!------------------------------------------------------------------------

!store the starting guess in the temporary array
Xtmp = Xold

! Perform the iterative solution
do m = 1, maxit
!    write(*,*) ' iteration Number = ', m
!    WRITe(*,*) 'Row, T1, T2, T4, T5, Xnew'
    do k = 1, n
        ! compute terms in the expression using if statements to 
        ! sort out which terms apply based on the indices
        if( k-LB .LT. 1 ) then; T1 = 0. ; else; T1 = A(k)*Xnew(k-LB); endif
        if( k    .EQ. 1 ) then; T2 = 0. ; else; T2 = B(k)*Xnew(k-1 ); endif
        if( k    .EQ. n ) then; T4 = 0. ; else; T4 = D(k)*Xtmp(k+1 ); endif
        if( k+UB .GT. n ) then; T5 = 0. ; else; T5 = E(k)*Xtmp(k+UB); endif
        ! Compute 
        Xnew( k ) = 1./C(k) * ( F(k) - T1 - T2  - T4 - T5 )
!        write(*,10) k, T1, T2, T4, T5, Xnew( k )
    end do
    !compute relative change for this iteration
    do k = 1, n
        relchng(k) = ( Xnew(k) - Xtmp(k) ) / Xtmp(k)
    end do
    ! check for convergence
!    write(dev,*) 'GAUSS_SIEDEL_PENTA: Iteration', m, ', Max rel change:', maxval(abs(relchng)) 
    if( maxval( abs( relchng ) ) .LT. tolit .AND. &
          F_L2_Norm( relchng, n) .LT. tolit          ) then
!            write(dev,*) 'GAUSS_SIEDEL_PENTA: Iterations required to converge: ', m
            numits = m
            exit ! exit iteration loop
!    elseif( maxval( abs(relchng) ) .GT. tolit ) then
     else
            Xtmp = Xnew
    endif
!end iteration loop
end do
! Print message if maximum number of iterations was exceeded
if ( m .gt. maxit ) then
    write(*,*) 'GAUSS_SEIDEL_PENTA: maximum number of iterations &
                 &exceeded; program will terminate.'
    STOP
endif

10 format( I7, 10f12.7 )

!---------------------------------------------------------------------------
end subroutine gauss_seidel_penta
!============================================================================
!   \\\\\\\\\\          E N D    S U B R O U T I N E     //////////
!   //////////     G A U S S _ S E I D E L _ P E N T A   \\\\\\\\\\
!============================================================================
      



!======================================================================
!       \\\\\\\\\\   B E G I N   S U B R O U T I N E   //////////
!       //////////          T H O M A S                \\\\\\\\\\
!======================================================================
SUBROUTINE THOMAS(A,B,C,D,X,N)
integer :: N
REAL A(N), B(N), C(N), D(N), X(N), Q(n+1), G(n+1)
REAL :: WI
integer :: i,j 
!
!       Purpose:        Solve a system of linear equations that appear
!                       as a tri-diagonal matrix.
!
!       Source:         This algorithm was handed out in class.
!       Written By:     Brad Eck
!       Revision 0:     Original coding on 19 Feb 09
!
!       A -- Main diagonal
!       B -- Superdiagonal
!       C -- Subdiagonal
!       D -- RHS vector
!       X -- Solution vector
!       N -- number of unknowns
!-----------------------------------------------------------------------

WI=A(1)
G(1)=D(1)/WI
DO I=2,n
        Q(I-1) = B(I-1)/WI
        WI = A(I) - C(I) * Q(I-1)
        G(I) = ( D(I) -C(I) * G(I-1))/WI
END DO

X(N) = G(N)

DO I=2,n
        J = N - I + 1
        X(J) = G(J) - Q(J) * X(J+1)
END DO

END SUBROUTINE THOMAS

!=========================================================================
!       \\\\\\\\\\    E N D   S U B R O U T I N E    //////////
!       //////////         T H O M A S               \\\\\\\\\\
!========================================================================

!============================================================================
!   \\\\\\\\\\\\\\\\\\                    //////////////////////
                            END MODULE solvers
!   //////////////////                    \\\\\\\\\\\\\\\\\\\\\\
!============================================================================


