! fortran_free_source
! 
!   This module is part of PERFCODE, written by Bradley J Eck.
!
!   Purpose:   Holds subroutines used in the 2D model
!
!   History:   April 2010 -- Original coding
!              July  2011 -- added PFC_PARAMETERS and WATER_BUDGET
!

MODULE pfc2Dsubs
IMPLICIT NONE

CONTAINS

! This module holds
!   1. Subroutine CELLWISE_PFC_PARAMS
!   2. Subroutine CELLWISE_BCs
!   2. Subroutine SET_ABCDEF
!   3. Subroutine SET_XYH
!   4. Subroutine DRAIN_SLOPE 
!   5. Subroutine CALC_Qx
!

!===================================================================
!   \\\\\\\\\\   B E G I N    S U B R O U T I N E    //////////
!   //////////  C E L L W I S E _ P F C _ P A R A M S  \\\\\\\\\\
!===================================================================
subroutine CELLWISE_PFC_PARAMS( params_vary )
!   Sets up the arrays for each PFC parameter based on whether
!       or not the parameters vary by cell
use shared, only: K_input, por_input, b_pfc_input, n_mann_input, &
                  hyd_cond, por, b_pfc, n_mann,                  & 
                  imax, jmax
use outputs, only: write_flipped_matrix
use  inputs, only: read_flipped_matrix 
!-------------------------------------------------------------------
implicit none
! VARIABLES
! Arguments
logical, intent( in ) :: params_vary
! Internal variables
logical :: have_arrays
character(len=20) :: K_infile, por_infile, b_infile, n_infile
!-----------------------------------------------------------------
! Assign input values to arrays
hyd_cond(:,:) =      K_input
     por(:,:) =    por_input
   b_pfc(:,:) =  b_pfc_input
  n_mann(:,:) = n_mann_input

! go further if the params vary
if( params_vary .eqv. .true. ) then

  ! default filenames for the parameter arrays
    K_infile = 'hyd_cond.csv'
  por_infile = 'por.csv'
    b_infile = 'b_pfc.csv'
    n_infile = 'n_mann.csv'

  write(*,*) 'CELLWISE_PFC_PARAMS: Hydraulic conductivity, porosity, thickness and Mannings n'
  write(*,*) '                are set to vary by grid cell.'
  write(*,*) '                Do you have paramater arrays to read in? (T/F)'
  read (*,*) have_arrays

  if( have_arrays .eqv. .false. ) then
      write(*,*) 'CELLWISE_PFC_PARAMS: Writing parameter arrays . . .'
      call write_flipped_matrix( hyd_cond, imax, jmax,   K_infile    )
      call write_flipped_matrix(      por, imax, jmax, por_infile    )
      call write_flipped_matrix(    b_pfc, imax, jmax,   b_infile    )
      call write_flipped_matrix(   n_mann, imax, jmax,   n_infile    )
      write(*,*) 'CELLWISE_PFC_PARAMS:  . . . please modify values of the appropriate cells. '
  end if

  !Hydraulic conductivity
  write(*,*) 'CELLWISE_PFC_PARAMS: Enter filename for hydraulic conductivity array' 
  write(*,*) '                     or press / for ', K_infile
  read(*,*) K_infile
  call read_flipped_matrix( hyd_cond, imax, jmax, K_infile )

  ! Porosity
  write(*,*) 'CELLWISE_PFC_PARAMS: Enter filename for porosity array'
  write(*,*) '                     or press / for ', por_infile
  read(*,*) por_infile
  call read_flipped_matrix( por, imax, jmax, por_infile )

  ! b_pfc 
  write(*,*) 'CELLWISE_PFC_PARAMS: Enter filename for thickness array'
  write(*,*) '                     array or press / for ', b_infile
  read(*,*) b_infile
  call read_flipped_matrix( b_pfc, imax, jmax, b_infile )
  
  ! Manning's n 
  write(*,*) 'CELLWISE_PFC_PARAMS: Enter filename for mannings n array'
  write(*,*) '                     or press / for ', n_infile
  read(*,*) n_infile
  call read_flipped_matrix( n_mann, imax, jmax, n_infile )

  write(*,*) 'CELLWISE_PFC_PARAMS: Parameter arrays loaded!'
  write(*,*) '                     outputting arrays for verification...'
  call write_flipped_matrix( hyd_cond, imax, jmax, 'VERIFY ' //   K_infile    )
  call write_flipped_matrix(      por, imax, jmax, 'VERIFY ' // por_infile    )
  call write_flipped_matrix(    b_pfc, imax, jmax, 'VERIFY ' //   b_infile    )
  call write_flipped_matrix(   n_mann, imax, jmax, 'VERIFY ' //   n_infile    )
   
endif 

end subroutine
!===================================================================
!   \\\\\\\\\\    E N  D   S U B R O U T I N E         //////////
!   //////////  C E L L W I S E _ P F C _ P A R A M S  \\\\\\\\\\
!===================================================================



!===================================================================
!   \\\\\\\\\\  B E G I N      S U B R O U T I N E   //////////
!   //////////      C E L L W I S E _ B C S          \\\\\\\\\\
!===================================================================
subroutine CELLWISE_BCs( BCs_vary )
!   Sets logical flags for boundary and corner cells.  Assigns input 
!       boundary condtions to arrays, outputs a summary
!       for the user to modify, and reads this back in for use
!       in model calculations.
use shared, only: is_boundary, is_corner, imax, jmax,    &
                  north_bc, south_bc, east_bc, west_bc,  & 
                  north_bc_input, south_bc_input,        &
                   east_bc_input, west_bc_input
use inputs, only: READ_BCs
use outputs, only: WRITE_BCs_AND_GRID
implicit none
! Arguments
logical :: BCs_vary
! Internal variables
logical :: have_file
character(len=20) :: filename = 'BCs_and_grid.csv'
integer :: i, j

!-----------------------
! Boundary & corner logical flags
is_boundary(:,:) = .false.

do i = 1, imax
    is_boundary(i, 1)    = .true.
    is_boundary(i, jmax) = .true.
end do

do j = 1, jmax
    is_boundary(1   , j) = .true.
    is_boundary(imax, j) = .true.
end do


is_corner(:, :)       = .false.
is_corner(1   , 1)    = .true.
is_corner(imax, 1)    = .true.
is_corner(1   , jmax) = .true.
is_corner(imax, jmax) = .true.




! Assign input values to arrays
north_bc(:) = north_bc_input
south_bc(:) = south_bc_input
 east_bc(:) =  east_bc_input
 west_bc(:) =  west_bc_input

! Go further if the BCs vary

if( BCs_vary .eqv. .true. ) then

  write(*,*) 'CELLWISE_BCs: Boundary conditions are set to vary by boundary cell.'
  write(*,*) '              Do you have a boundary file to read in? (T/F)'
  read (*,*) have_file


  if( have_file .eqv. .false. ) then
    write(*,*) 'CELLWISE_BCs: Writing boundary file . . .'
    call write_BCs_and_grid( filename )  
    write(*,*) 'CELLWISE_BCs: . . . please modify BCs as appropriate.' 
  end if

 
  write(*,*) 'CELLWISE_BCs: Enter filename for boundary file' 
  write(*,*) '              or press / for ', filename
  read(*,*) filename
  call READ_BCs( filename )

  write(*,*) '' 
  write(*,*) 'CELLWISE_BCs: Boundary conditions read!'
  write(*,*) '              ... outputting for verification.'
  write(*,*) '' 

  call write_BCs_and_grid( 'VERIFY' // filename )

elseif( BCs_vary .eqv. .false. ) then

  call write_BCs_and_grid( 'VERIFY' // filename )

end if

end subroutine CELLWISE_BCs
!===================================================================
!   \\\\\\\\\\    E N  D   S U B R O U T I N E       //////////
!   //////////      C E L L W I S E _ B C S          \\\\\\\\\\
!===================================================================





!===================================================================
!   \\\\\\\\\\   B E G I N   S U B R O U T I N E       //////////
!   //////////          S E T _ A B C D E F            \\\\\\\\\\
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


!========================================================================
!   \\\\\\\\\\   B E G I N    S U B R O U T I N E    //////////
!   //////////     D R A I N _ S L O P E             \\\\\\\\\\\
!=========================================================================
subroutine DRAIN_SLOPE(i, j, side, S_drain)
!   Computes the drainage slope vector for a grid cell. Since this is
!       primarily used at the boundary of the domain the side must
!       be specified.  The output is a 2D vector with components in the
!       longidudinal (ksi) and transverse (eta) directions on the grid.
! 
!   HISTORY:    12 Jul 2011 -- moved code from Subroutine MOC_KIN BC
!
use shared, only: Z, jmax, north, south, east, west, cv_info, lng, wid
use pfc2Dfuns, only: F_LinearIndex
use utilities, only: F_PythagSum
implicit none
! Arguments
integer, intent(in) :: i, j, side
real, dimension(2), intent(out) :: S_drain
! Internal variables
integer :: v, vi1, vj1, vjm1, vim1
real, dimension( 2 ) :: ksi_ii1  !vector in the ksi direction from point i to point i+1
real, dimension( 2 ) :: eta_jj1  !vector in the eta direction from point j to point j+1
real, dimension( 2 ) :: S_ksi    ! slope vector in the ksi direction
real, dimension( 2 ) :: S_eta    ! slope vector in the eta direction
!-----------------------------------------------------------------------
! setup to figure out the slope components in 
! the i (ksi)  and j (eta) directions 
v    = F_LinearIndex( i  , j  , jmax ) 
vi1  = F_LinearIndex( i+1, j  , jmax )
vj1  = F_LinearIndex( i  , j+1, jmax )
vjm1 = F_LinearIndex( i  , j-1, jmax )
vim1 = F_LinearIndex( i-1, j  , jmax )

!---------------------------------------------------
! Compute unit vectors in the  longitudinal (ksi)
!    and tranverse (eta) directions.  If statements
!    are careful around the boundaries
!---------------------------------------------------


! LONGITUDINAL DIRECTION (ksi)
if( (side == north) .or. (side == south) .or. (side == west) ) then
    ksi_ii1 = (/  CV_info(vi1)%X - CV_Info(v)%X ,  &
                  CV_info(vi1)%Y - CV_Info(v)%Y      /)

elseif( side == east )  then
    ksi_ii1 = (/  CV_info(v)%X - CV_Info(vim1)%X ,  &
                  CV_info(v)%Y - CV_Info(vim1)%Y      /)

endif


! TRANSVERSE DIRECTION (eta)
if( side == south ) then

        eta_jj1 =   (/  CV_info(vj1)%X - CV_Info(v)%X ,  &
                        CV_info(vj1)%Y - CV_Info(v)%Y      /)

elseif( side == north ) then

        eta_jj1 = - (/  CV_info(vjm1)%X - CV_Info(v)%X,  &
                        CV_info(vjm1)%Y - CV_Info(v)%Y     /)

elseif( side == east .or. side == west ) then

    if( j  /= jmax ) then
        ! j+1 is OK
        eta_jj1 = (/  CV_info(vj1)%X - CV_Info(v)%X ,  &
                      CV_info(vj1)%Y - CV_Info(v)%Y      /)
    elseif( j == jmax ) then
        ! special treatment for jmax
        eta_jj1 = (/  CV_info(v)%X - CV_Info(vjm1)%X ,  &
                      CV_info(v)%Y - CV_Info(vjm1)%Y      /)
    endif

endif


!write(device,*) 'Direction Vectors: ksi_ii1 = ', ksi_ii1,  & 
!                                  ' eta_jj1 = ', eta_jj1

!Make the direction vectors of unit length
ksi_ii1 = ksi_ii1 / F_PythagSum( ksi_ii1(1), ksi_ii1(2) )
eta_jj1 = eta_jj1 / F_PythagSum( eta_jj1(1), eta_jj1(2) )

!write(device,*)'Direction UNIT Vectors: ksi_ii1 = ', ksi_ii1, &
!                                      ' eta_jj1 = ', eta_jj1


!---------------------------------------------------------
! Compute a slope vector for each direction by
!   estimating the magnitude and using the unit vectors
!   obtained above for the directions
!---------------------------------------------------------

! LONGITUDINAL DIRECTION (ksi)
if(  side == north .or. side == south .or. side == west ) then

    S_ksi = ksi_ii1 * ( Z( i+1, j   ) - Z( i, j ) ) / lng( i, j )

elseif( side == east ) then

    S_ksi = ksi_ii1 * ( Z( i, j   ) - Z( i-1, j ) ) / lng( i, j )

end if


! TRANSVERSE DIRECTION (eta)
if( side == south) then

    S_eta =   eta_jj1 * ( Z( i  , j+1 ) - Z( i, j ) ) / wid( i, j )

elseif( side == north ) then

    S_eta = - eta_jj1 * ( Z( i  , j-1 ) - Z( i, j ) ) / wid( i, j )

elseif( side == east .or. side == west ) then

    if( j /= jmax ) then
        
        S_eta =   eta_jj1 * ( Z( i  , j+1 ) - Z( i, j ) ) / wid( i, j )

    elseif( j == jmax ) then

        S_eta =  eta_jj1 * (  Z( i, j) - Z( i, j-1) ) / wid( i, j )

    endif

end if

! compute drainage slope vector
S_drain = S_ksi + S_eta

end subroutine DRAIN_SLOPE
!========================================================================
!   \\\\\\\\\\   E N D        S U B R O U T I N E    //////////
!   //////////     D R A I N _ S L O P E             \\\\\\\\\\
!=========================================================================






!============================================================================
!      \\\\\\\\\\  B E G I N   S U B R O U T I N E    //////////
!      //////////           C A L C  _ Q x            \\\\\\\\\\
!============================================================================
Subroutine CALC_Qx( i, j, face, Qx)
!   Computes the average flow rate out of a cell face during the current
!       time step.  This routine is used by Subroutine WATER_BUDGET
!       to compute the flow rate out of boundary cells based on the
!       solution at a time step.
!
!   HISTORY:    20 July 2011 -- Original coding
!
use shared, only: north, south, east, west, old, itr, &
                  h_old, h_new, z, area
use convcoef
implicit none
! Arguments
integer, intent(in)  :: i,j, face
real   , intent(out) :: Qx
! Internal variables
integer :: iout, jout
real :: cc, cc1
real :: head_diff_old, head_diff_new
real :: Qx_old, Qx_new
!-------------------------------------------

! Get indicees for the face
if( face == west ) then

    iout = i-1
    jout = j

elseif( face == north ) then

    iout = i
    jout = j+1

elseif( face == east ) then

    iout = i+1
    jout = j

elseif( face == south ) then

    iout = i
    jout = j-1

endif

! Conveyance coefs
call conveyance( face, old, i, j, cc )
call conveyance( face, itr, i, j, cc1 )

! Head difference across the face
head_diff_old = (h_old(iout,jout) - h_old(i,j)) + (z(iout,jout)  - z(i,j))
head_diff_new = (h_new(iout,jout) - h_new(i,j)) + (z(iout,jout)  - z(i,j))

! Flow rates
Qx_old = cc  * head_diff_old * area(i,j)
Qx_new = cc1 * head_diff_new * area(i,j)
Qx     = (Qx_old + Qx_new ) / 2.

!----------------------------------------------------------------------------
end subroutine
!============================================================================
!      \\\\\\\\\\     E N D    S U B R O U T I N E        //////////
!      //////////          C A L C _ Q x                  \\\\\\\\\\
!============================================================================



 

END MODULE pfc2Dsubs













