! fortran_free_source
!
! (c) Copyright 2010, 2015 Bradley J. Eck
! This module is part of PERFCODE 
!

!======================================================================
!======================================================================
!   \\\\\\\\\\   B E G I N   M O D U L E           ///////////
!   //////////                  U T I L I T I E S  \\\\\\\\\\\
!======================================================================
module utilities
implicit none
contains

! This module holds subroutines and functions for various jobs:
!   1.  Subroutine  GET_BANDS
!   2.  Subroutine  PUT_BANDS
!   3.  Subroutine  UNLINEARIZE
!   4.  Subroutine  BILINEAR_INTERP
!   5.  Function    F_LINTERP
!   6.  Function    F_L2_NORM
!   7.  Function    F_PYTHAGSUM
!   8.  Function    F_EXTRAPOLATE
!=======================================================================


!=======================================================================
!   \\\\\\\\\\      B E G I N    F U N C T I O N   ///////////
!   //////////       G E T _ B A N D S             \\\\\\\\\\\
!=======================================================================
!
!   PURPOSE:    Extracts the five bands from a penta-diagonal matrix.
!
SUBROUTINE GET_BANDS( COEF, N, LB, UB, A, B, C, D, E)
!
! COEF -- Penta-diagonal coefficient matrix.
!   N  -- number of unknowns (size of system)
!  LB  -- lower bandwidth
!  UB  -- upper bandwidth
!  A,B -- lower bands of the penta-diagonal matrix
!   C  -- main diagonal
!  D,E -- upper bands of the penta-diagonal matrix
!----------------------------------------------------------------------
!VARIABLE DECLARATIONS
!   Arguments
integer, intent( in ) :: N, LB, UB
real,    intent( in ) :: COEF( N, N )
real,    intent( out) :: A(N), B(N), C(N), D(N), E(N)
!   Internal variables
integer ::  i   !looping variable
!---------------------------------------------------------------------

! Lowermost subdiagonal 
do i = LB+1, n
    A(i) = coef(i,i-LB)
end do

! Subdiagonal
do i = 2, n
    B(i) = coef(i,i-1)
end do

! Main Diagonal
do i = 1, n
     C(i) = coef(i,i)
end do

! Super diagonal
do i = 1, n-1
     D(i) = coef(i,i+1)
end do


! Uppermost diagonal
do i = 1, n - UB
     E(i) = coef(i,i+ub)
end do
!---------------------------------------------------------------------
end subroutine GET_BANDS
!======================================================================
!   \\\\\\\\\\     E N D       S U B R O U T I N E ///////////
!   //////////       G E T _ B A N D S             \\\\\\\\\\\
!======================================================================

!======================================================================
!   \\\\\\\\\\    B E G I N    S U B R O U T I N E ///////////
!   //////////       P U T  _ B A N D S            \\\\\\\\\\\
!======================================================================
!
!   PURPOSE:    Puts the five bands into a square matrix.
!
SUBROUTINE PUT_BANDS(A, B, C, D, E, N, LB, UB, COEF)
!
!  A,B -- lower bands of the penta-diagonal matrix
!   C  -- main diagonal
!  D,E -- upper bands of the penta-diagonal matrix
!   N  -- number of unknowns (size of system)
!  LB  -- lower bandwidth
!  UB  -- upper bandwidth
! COEF -- Penta-diagonal coefficient matrix.
!----------------------------------------------------------------------
!VARIABLE DECLARATIONS
!   Arguments
integer, intent( in ) :: N, LB, UB
real,    intent( in  ) :: A(N), B(N), C(N), D(N), E(N)
real,    intent( out ) :: COEF( N, N )
!   Internal variables
integer ::  i   !looping variable
!---------------------------------------------------------------------
! Fill coefficient matrix with zeros

coef(:,:) = 0.0

! Lowermost subdiagonal 
do i = LB+1, n
    coef(i,i-LB) = A(i)
end do

! Subdiagonal
do i = 2, n
    coef(i,i-1) = B(i)
end do

! Main Diagonal
do i = 1, n
     coef(i,i) = C(i)
end do

! Super diagonal
do i = 1, n-1
     coef(i,i+1) = D(i)
end do


! Uppermost diagonal
do i = 1, n - UB
     coef(i,i+ub) = E(i)
end do
!---------------------------------------------------------------------
end subroutine PUT_BANDS
!======================================================================
!   \\\\\\\\\\     E N D       S U B R O U T I N E ///////////
!   //////////       P U T _ B A N D S             \\\\\\\\\\\
!======================================================================


!======================================================================
!   \\\\\\\\\\     B E G I N      S U B R O U T I N E ///////////
!   //////////       U N L I N E A R I Z E            \\\\\\\\\\\
!======================================================================
subroutine unlinearize( vector, imax, jmax, vmax, matrix )

!   Puts unknowns in linear (vector) form into matrix form
!   Assumes column-wise ordering from southwest corner of domain
use pfc2Dfuns, only: F_LinearIndex
implicit none
!Arguments
integer,                     intent( in ):: imax, jmax, vmax
real, dimension(vmax),       intent(in)  :: vector
real, dimension(imax, jmax), intent(out) :: matrix
!Internal Variables
integer                                  :: i, j, v

DO j = 1, jmax
    DO i = 1, imax
        v = F_LinearIndex( i, j, jmax )
        matrix(i,j) = vector( v )
    end do
end do

end subroutine unlinearize
!=======================================================================
!   \\\\\\\\\\      E N D         S U B R O U T I N E ///////////
!   //////////       U N L I N E A R I Z E            \\\\\\\\\\\
!=======================================================================

!=======================================================================
!   \\\\\\\\\\   B E G I N    S U B R O U T I N E     //////////
!   //////////       B I L I N E A R _ I N T E R P    \\\\\\\\\\
!=======================================================================
subroutine BILINEAR_INTERP       (  X , Y , Z ,  &
                                    x1, y1, z1,  &
                                    x2, y2, z2,  &
                                    x3, y3, z3,  &
                                    x4, y4, z4,  &
                                    dev, error       )

! Finds the value of Z at the point X,Y using Finite Element
! style interpolation with a Bi-linear element.  The physical
! coordinates (x,y,z) are mapped into ksi, eta space that
! ranges from -1 to 1.  The values of ksi and eta for the point
! X, Y are found by solving the non-linear system using the
! Newton-Raphson method.
!
!

implicit none
!VARIABLE DECLARATIONS
!   Arguments
real, intent( in  ) :: X, Y         ! Coordinates of point where Z is desired
real, intent( out ) ::       Z      ! Unknown function value
real, intent( in  ) :: x1, y1, z1   ! Coordinates of point 1
real, intent( in  ) :: x2, y2, z2   !   "    "       point 2
real, intent( in  ) :: x3, y3, z3   !   "    "       point 3
real, intent( in  ) :: x4, y4, z4   !   "    "       point 4
integer, optional   :: dev          ! output device for writing errors
logical, optional   :: error 

!   Internal variables
real :: ksi, eta                    ! mapped coordinates of XY
real :: X_guess, Y_guess            ! Values of X and Y computed from ksi and eta
real :: delta_ksi, delta_eta        ! incremental change in values over iteration
real :: J_11, J_12, J_21, J_22      ! elements of the jacobian matrix
real :: PSI_1, PSI_2, PSI_3, PSI_4  ! Shape functions for BiLinear element
real, parameter :: tolit = 1.e-5   ! iteration tolerance
integer, parameter :: qmax = 10    ! maximum number of iterations
integer            :: q             ! looping variable
integer            :: device        ! output device        

!--------------------------------------------------------------------------------

! Default values for output device
if( present( dev ) .EQV. .FALSE. ) then
        device = 6
else
        device = dev
end if
    

! STEP 1:  Find the value of ksi and eta that correspond to the point X,Y
!  initial guess for ksi and eta is in the middle of the element ( 0,0 )
ksi = 0.0
eta = 0.0

Map: do q = 1, qmax
        
        !   Values of the shape functions at the point (X, Y)
        PSI_1 = 0.25 * ( 1. - ksi ) * ( 1. - eta )
        PSI_2 = 0.25 * ( 1. + ksi ) * ( 1. - eta )
        PSI_3 = 0.25 * ( 1. + ksi ) * ( 1. + eta )
        PSI_4 = 0.25 * ( 1. - ksi ) * ( 1. + eta )
        
        !  figure out value of X and Y using ksi and eta
        X_guess = x1*PSI_1 + x2*PSI_2 + x3*PSI_3 + x4*PSI_4
        Y_guess = y1*PSI_1 + y2*PSI_2 + y3*PSI_3 + y4*PSI_4


        !compute values of jacobian
        !J_11 = d X_guess / d ksi
        J_11 =   x1 / 4. * ( eta - 1.  ) &
               + x2 / 4. * ( 1.  - eta ) &
               + x3 / 4. * ( eta + 1.  ) &
               - x4 / 4. * ( eta + 1.  )

        !J_12 = d X_guess / d eta
        J_12 =   x1 / 4. * ( ksi - 1.  ) &
               - x2 / 4. * ( ksi + 1.  ) &
               + x3 / 4. * ( ksi + 1.  ) &
               + x4 / 4. * (  1. - ksi )

        !J_21 = d Y_guess / d ksi
        J_21 =   y1 / 4. * ( eta - 1.  ) &
               + y2 / 4. * (  1. - eta ) &
               + y3 / 4. * ( eta + 1.  ) &
               - y4 / 4. * ( eta + 1.  )

        !J_22 = d Y_guess / d eta
        J_22 =   y1 / 4. * ( ksi - 1.  ) &
               - y2 / 4. * ( ksi + 1.  ) &
               + y3 / 4. * ( ksi + 1.  ) &
               + y4 / 4. * ( 1.  - ksi )

        !Manual solution of 2 x 2 system:  J * delta_ksi/eta = X/Y_guess - X/Y
        delta_ksi = (   (X_guess - X)*J_22       &
                      - (Y_guess - Y)*J_12   ) / &
                    (   J_11 * J_22              &
                      - J_12 * J_21          )

        delta_eta = (   (Y_guess - Y)*J_11       &
                      - (X_guess - X)*J_21   ) / &
                    (   J_11 * J_22              &
                      - J_12 * J_21          )


!write(device,*) 'BILINER_INTERP q=', q, 'delta_ksi =', delta_ksi, ' delta_eta', delta_eta

        ! update vales of ksi and eta
        ! rembeber delta = ksi_q - ksi_q+1
        ksi = ksi - delta_ksi
        eta = eta - delta_eta

        !Convergence Test
        if( abs( delta_ksi ) .LT. tolit .AND.   &
            abs( delta_eta ) .LT. tolit          ) then

                exit Map
        endif

end do Map



!report mapping result
!write( device, * ) 'BILINEAR_INTERP: Mapping result: ksi =', ksi, ' eta = ', eta

! assume no error and change if there is one
if( present( error) .eqv. .TRUE. ) then
    error = .FALSE.
end if

! Give Error if iteration fails to converge
if( q .GT. qmax ) then
    write( device, * ) 'BILINEAR_INTERP: Mapping iteration failed. ksi =', ksi, ' eta = ', eta

    if( present( error) .eqv. .TRUE. ) then   !assign an error if the variable was provided.
        error = .TRUE.
    end if
end if

                        
! Confirm that mapped point lies inside the range of the datapoints
if( abs( ksi ) .GT. 1. + tolit  .OR. &
    abs( eta ) .GT. 1. + tolit        )  then
    write( device, * ) 'BILINEAR_INTERP: Desired point lies outside &
                        & known points: ksi =', ksi, ' eta = ', eta 
    if( present( error) .eqv. .TRUE. ) then   !assign an error if the variable was provided.
        error = .TRUE.
    end if
end if

! STEP 2: Having found the values of ksi and eta that correspond
!         to the point (X, Y) compute the value of Z at that location.

!   Values of the shape functions at the point (X, Y)
PSI_1 = 0.25 * ( 1. - ksi ) * ( 1. - eta )
PSI_2 = 0.25 * ( 1. + ksi ) * ( 1. - eta )
PSI_3 = 0.25 * ( 1. + ksi ) * ( 1. + eta )
PSI_4 = 0.25 * ( 1. - ksi ) * ( 1. + eta )

! Value of Z at the point (X, Y)
Z = z1*PSI_1 + z2*PSI_2 + z3*PSI_3 + z4*PSI_4

!write( device, *) 'BILINEAR_INTERP: PSI_1=', PSI_1, &
!                                  ' PSI_2=', PSI_2, &
!                                   'PSI_3=', PSI_3, &
!                                   'PSI_4=', PSI_4, &
!                                      ' Z=', Z

!----------------------------------------------------------------------
end subroutine BILINEAR_INTERP
!=======================================================================
!   \\\\\\\\\\    E N D    S U B R O U T I N E         //////////
!   //////////        B I L I N E A R _ I N T E R P    \\\\\\\\\\
!======================================================================


!======================================================================
!   \\\\\\\\\\      B E G I N    F U N C T I O N   ///////////
!   //////////       F _ L I N T E R P             \\\\\\\\\\\
!======================================================================
! linear interpolation function
function F_linterp( X, known_X, known_Y, n) Result( Y )
!VARIABLE DECLARATIONS
!   Arguments
integer, intent( in ) :: n
real   , intent( in ) :: X
real, dimension(n), intent(in) :: known_X, known_Y
!   Internal Variables
real :: Y
integer :: i1, i2, j, im
!------------------------------------------------------------------
! bi-section method to find the right place in the table
! initialize indices
i1 = 1
i2 = n


if( known_X(n) .GT. known_X(1)  ) then
!ASCENDING ORDER
    do j = 1, 1000
        if ( i2 - i1 .gt. 1 ) then
            im = (i1+i2)/2      !midpoint
            if ( X .eq. known_X(im) ) then
                i1 = im
                i2 = im + 1
            elseif( X .gt. known_X(im) ) then
                i1 = im
            elseif( X .lt. known_X(im) ) then
                i2 = im
            endif    
        else
            exit
        end if
    end do

elseif( known_X(n) .LT. known_X(1) )then
!DESCENDING ORDER
    do j = 1, 1000
        if ( i2 - i1 .gt. 1 ) then
            im = (i1+i2)/2      !midpoint
            if ( X .eq. known_X(im) ) then
                i1 = im
                i2 = im + 1
            elseif( X .gt. known_X(im) ) then
                i2 = im
            elseif( X .lt. known_X(im) ) then
                i1 = im
            endif    
        else
            exit
        end if
    end do

end if


!        WRITE(*,*) 'j=', j, 'im=', im, 'i1=', i1, 'i2=', i2

if( j .eq. 1000 ) then
    write( *,* ) 'F_LINTERP: Arrays too large for this routine, &
                   & increase number of searching steps and recompile.'
endif


! bounds found; compute interpolated value
Y  = (X - Known_X(i1)) / &
     (Known_X(i2) - Known_X(i1)) * &
     (Known_Y(i2) - Known_Y(i1)) + Known_Y(i1)


!----------------------------------------------------------------------
end function  F_LINTERP          
!======================================================================
!   \\\\\\\\\\      E N D    F U N C T I O N       ///////////
!   //////////       F _ L I N T E R P             \\\\\\\\\\\
!======================================================================


!======================================================================
!   \\\\\\\\\\      B E G I N    F U N C T I O N   ///////////
!   //////////       F _ L 2 _ N O R M             \\\\\\\\\\\
!======================================================================
function F_L2_NORM( vector, n ) result( L2 )
!   Computes the L2 norm of a real-valued vector with n elements
!
! Variable Declarations
implicit none
! Arguments
integer,                 intent( in ) :: n
real   , dimension( n ), intent( in ) :: vector
! Internal Variables
real :: L2
real, dimension( n ) :: squares
integer :: i
!---------------------------------------------------------------------
do i = 1, n
    squares(i) = vector(i) ** 2
end do

L2 = sqrt( sum( squares(:) ) )
!---------------------------------------------------------------------
end function F_L2_NORM
!======================================================================
!   \\\\\\\\\\      E N D    F U N C T I O N       ///////////
!   //////////       F _ L 2 _ N O R M             \\\\\\\\\\\
!=====================================================================


!=============================================================
!   \\\\\\\\\\  B E G I N    F U N C T I O N       ///////////
!   //////////       F _ P y t h a g S u m         \\\\\\\\\\\
!=============================================================
Function F_PythagSum( x, y)
! Computes the pythagorean sum of twovariables
implicit none
REAL :: x, y, F_PythagSum
F_PythagSum = sqrt( x**2 + y**2 )
end function F_PythagSum
!==============================================================
!   \\\\\\\\\\   E N D       F U N C T I O N       ///////////
!   //////////       F _ P y t h a g S u m         \\\\\\\\\\\
!=============================================================


!=============================================================
!   \\\\\\\\\\  B E G I N    F U N C T I O N       ///////////
!   //////////       F _ E X T R A P O L A T E     \\\\\\\\\\\
!=============================================================
Function F_Extrapolate( X, x1, y1, x2, y2) RESULT(Y)
!   Finds the value of Y corresponding to the location X
!   on the line passing through ( x1, y1) and (x2, y2)
!
!   Called from:    convcoef@frictionslope


implicit none
REAL, intent(in) :: X, x1, y1, x2, y2
REAL :: Y
REAL :: slope, intercept
slope = (y2 - y1) / (x2 - x1)
intercept = y1 - slope * x1
Y = slope * X + intercept
end function F_Extrapolate
!==============================================================
!   \\\\\\\\\\   E N D       F U N C T I O N       ///////////
!   //////////       F _ E X T R A P O L A T E     \\\\\\\\\\\
!=============================================================


!==============================================================
end module utilities
!==============================================================
!   \\\\\\\\\\   E N D       M O D U L E           ///////////
!   //////////       U T I L I T I E S             \\\\\\\\\\\
!=============================================================
!==============================================================
