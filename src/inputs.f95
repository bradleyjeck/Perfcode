! fortran_free_source
!
! (c) Copyright 2010, 2015 Bradley J Eck 
! This module is part of PERFCODE

!   Purpose:    This module contains subroutines to read input files


module inputs

implicit none

contains

!   1. GET_PARAMETERS
!   2. GET_RAINFALL

!===================================================================
!   \\\\\\\\\\       B E G I N    S U B R O U T I N E ///////////
!   //////////          G E T _ P A R A M E T E R S   \\\\\\\\\\\
!==================================================================
!
!   Purpose:    This subroutine reads problem parameters
!               from a user selected input file.
subroutine GET_PARAMETERS( K, por, b_pfc, n_mann, g, dt_pfc, dt_sheet, max_time, &
                           dx, dy, qmax, maxit, h0, eps_matrix, eps_itr, eps_ss, &
                           relax, relax_tran,                                    & 
                           north_bc, south_bc, east_bc, west_bc,                 &
                           animate, dt_ani                                         )
!
!   K   --  Darcy Hydraulic Conductivity
!   por --  Effective porosity of the PFC
!  b_pfc--  Thickness of the pfc
! n_mann--  Manning's n
!   g   --  gravitational acceleration
!-------------------------------------------------------------------
! VARIABLE DECLARATIONS
!Arguments
REAL, intent( out ) :: K, por, b_pfc, n_mann, g, dt_pfc, dt_sheet, max_time, dx, dy  
integer, intent(out) :: qmax, maxit
real, intent( out ) :: h0, eps_matrix, eps_itr, eps_ss, relax, relax_tran
character(7), intent(out) :: north_bc, south_bc, east_bc, west_bc
logical, intent( out)  :: animate
real, intent( out ) :: dt_ani
! Internal variables
CHARACTER(20) infile    ! the file name to read parameters from
CHARACTER(5) :: dummy_line
!---------------------------------------------------------------
! Executable

! default value for input file
infile = 'parameters.dat'

! Prompt the user for the input file
WRITE(*,*) 'Enter filename or press / for ', infile
READ(*,*) infile

!read the file
OPEN( UNIT=8,   FILE = infile, ACTION = 'read', STATUS = 'old' )  

READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * ) dummy_line
! PFC Properties
READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * )   K
READ( unit=8, fmt = * )   por
READ( unit=8, fmt = * )   b_pfc
READ( unit=8, fmt = * )   n_mann

!Physical Constants
READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * )   g

!Timesteps
READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * )   dt_pfc
READ( unit=8, fmt = * )   dt_sheet
READ( unit=8, fmt = * ) max_time

! preliminary grid spacing
READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * ) dx
READ( unit=8, fmt = * ) dy

! tolerences
READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * ) qmax
READ( unit=8, fmt = * ) maxit
READ( unit=8, fmt = * ) eps_matrix
READ( unit=8, fmt = * ) eps_itr
READ( unit=8, fmt = * ) eps_ss
READ( unit=8, fmt = * ) relax
READ( unit=8, fmt = * ) relax_tran


! inital depth
READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * ) h0

! Boundary conditions
READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * ) north_bc
READ( unit=8, fmt = * ) south_bc
READ( unit=8, fmt = * ) east_bc
READ( unit=8, fmt = * ) west_bc

!Animation options
READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * ) animate
READ( unit=8, fmt = * ) dt_ani

close( 8 )

!------------------------------------------------------------------
end subroutine GET_PARAMETERS
!===================================================================
!   \\\\\\\\\\       E N D        S U B R O U T I N E ///////////
!   //////////          G E T _ P A R A M E T E R S   \\\\\\\\\\\
!==================================================================



!===================================================================
!   \\\\\\\\\\       B E G I N    S U B R O U T I N E ///////////
!   //////////          G E T _  R A I N F A L L      \\\\\\\\\\\
!==================================================================
!
!   Purpose:    This subroutine reads the rainfall record
!               from a user selected input file.

SUBROUTINE GET_RAINFALL(max_rec, rain_time, rain_rate, nrr)

INTEGER, intent( in) :: max_rec ! maximum allowable number of rainfall records
INTEGER, intent( out ) :: nrr   !Number of rainfall records
REAL, DIMENSION( max_rec), intent( inout) :: rain_time, rain_rate

! Internal variables
CHARACTER(30) infile    ! the file name to read parameters from
integer :: i, j     ! looping variables


!-------------------------
! read in the rainfall data 

! default value for input file
infile = 'rainfall.dat'

! Prompt the user for the input file
WRITE(*,*) 'Enter filename or press / for ', infile
READ(*,*) infile

OPEN( UNIT=8,   FILE = infile, ACTION = 'read', STATUS = 'old' )  

! Rainfall Rate
read( unit=8, fmt = * ) nrr

if ( nrr .gt. size ( rain_time ) ) then
    print *, 'GET_RAINFALL: Too many rainfall records--increase array size and recompile'
else
    do i = 1, nrr
        READ( unit = 8 , fmt = * ) j, rain_time(j), rain_rate(j)
    end do
end if

close( 8 )

!------------------------------------------------------------------
end subroutine GET_RAINFALL
!===================================================================
!   \\\\\\\\\\       E N D        S U B R O U T I N E ///////////
!   //////////          G E T _ R A I N F A L L       \\\\\\\\\\\
!==================================================================


end module inputs
