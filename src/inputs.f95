! fortran_free_source

!   Purpose:    This module contains subroutines to read input files


module inputs

implicit none

contains

!   1. Subroutine GET_PARAMETERS
!   2. Subroutine GET_RAINFALL
!   3. Subroutine GET_CROSS_SECTION
!   4. Subroutine GET_LONG_PROFILE
!   5. Subroutine READ_FLIPPED_MATRIX
!   6. Subroutine READ_BCs

!===================================================================
!   \\\\\\\\\\       B E G I N    S U B R O U T I N E ///////////
!   //////////          G E T _ P A R A M E T E R S   \\\\\\\\\\\
!==================================================================
!
subroutine GET_PARAMETERS(  )
!   Read parameter input file
USE shared, ONLY: params_vary, K_input, por_input, b_pfc_input,  &
                  n_mann_input, g, dt_pfc, dt_sheet, max_time,   & 
                  dx, dy, qmax, maxit, h0, eps_matrix, eps_itr,  &
                  relax, relax_tran, BCs_vary,                   &
                  north_bc_input, south_bc_input, east_bc_input, &
                  west_bc_input, animate, dt_ani, model_mode                               
implicit none
!-------------------------------------------------------------------
! VARIABLE DECLARATIONS
! Internal variables
CHARACTER(20):: infile      ! file name to read parameters from
CHARACTER(5) :: dummy_line  
!---------------------------------------------------------------


! default value for input file
infile = 'parameters.dat'

! Prompt the user for the input file
WRITE(*,*) 'GET_PARAMETERS: Enter filename or press / for ', infile
READ(*,*) infile

!read the file
OPEN( UNIT=8,   FILE = infile, ACTION = 'read', STATUS = 'old' )  

READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * ) dummy_line
! PFC Properties
READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * ) params_vary 
READ( unit=8, fmt = * ) K_input
READ( unit=8, fmt = * ) por_input
READ( unit=8, fmt = * ) b_pfc_input
READ( unit=8, fmt = * ) n_mann_input

!Physical Constants
READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * ) g

!Timesteps
READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * ) dt_pfc
READ( unit=8, fmt = * ) dt_sheet
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
READ( unit=8, fmt = * ) relax
READ( unit=8, fmt = * ) relax_tran


! inital depth
READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * ) h0

! Boundary conditions
READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * ) BCs_vary
READ( unit=8, fmt = * ) north_bc_input
READ( unit=8, fmt = * ) south_bc_input
READ( unit=8, fmt = * ) east_bc_input
READ( unit=8, fmt = * ) west_bc_input

!Animation options
READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * ) animate
READ( unit=8, fmt = * ) dt_ani

!Model mode
READ( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * ) model_mode

close( 8 )


if(  (params_vary .eqv. .true.) .and. (model_mode .eq. '1D') ) then
    write(*,*) 'GET_PARAMETERS:  Variable parameters not allowed for 1D model;'
    write(*,*) '                 please revise parameters file.'
    write(*,*) '                 STOPPING PROGRAM.'
    STOP
endif

if( (BCs_vary .eqv. .true. ) .and. (model_mode .eq. '1D') ) then
    write(*,*) 'GET_PARAMETERS:  Variable boundary conditions not allowed for 1D model;'
    write(*,*) '                 please revise parameters file.'
    write(*,*) '                 STOPPING PROGRAM.'
    STOP
endif





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



!===================================================================
!   \\\\\\\\\\   B E G I N        S U B R O U T I N E  ///////////
!   //////////    G E T _ C R O S S _ S E C T I O N    \\\\\\\\\\\
!==================================================================

subroutine GET_CROSS_SECTION( )
!   Reads the cross section input file and echos the results to 
!       the screen.
!
!   History:  July 6, 2011 -- moved this from SET_ELEVATIONS into
!                             its own subroutine so PARALLEL_GRID
!                             could use it as well.
use shared, only: nr_cs, slope_cs, wid_cs
implicit none
! Arguments
! Internal Variables
integer :: i, j
character(len=20) :: infile, dummy_line
!-----------------------------------------------------------------

! default value for input file
infile = 'CrossSection.dat'

! Prompt the user for the input file
write(*,*) '' 
WRITE(*,*) 'GET_CROSS_SECTION: Enter filename or / for ', infile
READ(*,*) infile

! Read the file
OPEN( UNIT=8,   FILE = infile, ACTION = 'read', STATUS = 'old' ) 

!Cross Seection geometry
read( unit=8, fmt = * ) dummy_line
READ( unit=8, fmt = * ) nr_cs
read( unit=8, fmt = * ) dummy_line
if( nr_cs .gt. size( slope_cs) ) then
    write(*,*) 'GET_CROSS_SECTION: Too many records in', infile, &
                'increase array size and recompile'

else
   do i = 1, nr_cs
       READ( unit=8, fmt = * ) j, slope_cs(j), wid_cs(j)
   end do
end if

! Close the input file
close(8)

! Echo to screen
PRINT *, 'CROSS SECTION INPUTS '
WRITE(*,*)  ' Segment    Slope       Width  '
WRITE(*,*)  '=================================='

DO j = 1, nr_cs
       WRITE(*,10) j, slope_cs(j), wid_cs(j)
END DO

write(*,*) ' '

10     FORMAT('   ', ( i3, '    '),  ( F10.3, '  ') , F10.6, i12  ) 

end subroutine GET_CROSS_SECTION
!===================================================================
!   \\\\\\\\\\    E N D       S U B R O U T I N E      ///////////
!   //////////    G E T _ C R O S S _ S E C T I O N    \\\\\\\\\\\
!==================================================================


!===================================================================
!   \\\\\\\\\\       B E G I N    S U B R O U T I N E   ///////////
!   //////////   G E T _ L O N G _ P R O F I L E        \\\\\\\\\\\
!==================================================================
subroutine GET_LONG_PROFILE( )
!   Reads the longitudinal profile input file and echos the result
!       to the screen
!
!   History: July 6, 2011  -- Moved from SET_ELEVATIONS to its own
!                             subroutine so that PARALLEL_GRID
!                             could use it.
use shared, only: nr_lp, dist_lp, Z_lp

! Arguments
! Internal variables
character(len=20) :: infile, dummy_line
integer :: i, j
!---------------------------------------------------------
! default value for input file
infile = 'LongProfile.dat'

! Prompt the user for the input file
WRITE(*,*) ''
WRITE(*,*) 'Enter filename or press / for ', infile
READ(*,*) infile

! Read the file
OPEN( UNIT=8,   FILE = infile, ACTION = 'read', STATUS = 'old' )  

read( unit=8, fmt = * ) dummy_line
read( unit=8, fmt = * ) nr_lp    !number of rows to define cross section
read( unit=8, fmt = * ) dummy_line

if ( nr_lp .gt. size ( dist_lp ) ) then
    print *, 'SET_ELEVATIONS: Too many records in', infile, &
               'increase array size and recompile'
else
    do i = 1, nr_lp
        READ( unit = 8 , fmt = * ) j, dist_lp(j), Z_lp(j)
    end do
end if


close( 8 )

! Echo to screen
PRINT *, 'LONGITUDINAL PROFILE '
WRITE(*,*)  ' Point    Distance     Elevation  '
WRITE(*,*)  '=================================='
!            '----|----|----|----|----|----|----|----|----|----|----|----|
!                 5   10   15   20   25   30   35   40   45   50   55   60
DO i = 1, nr_lp
       WRITE(*,10) i, dist_lp(i), Z_lp(i)
END DO

write(*,*) ''

10     FORMAT('   ', ( i3, '    '),  ( F10.3, '  ') , F10.6, i12   )

end subroutine GET_LONG_PROFILE
!===================================================================
!   \\\\\\\\\\   E N D      S U B R O U T I N E         ///////////
!   //////////   G E T _ L O N G _ P R O F I L E        \\\\\\\\\\\
!==================================================================




!===================================================================
!   \\\\\\\\\\       B E G I N    S U B R O U T I N E   ///////////
!   //////////   R E A D _ F L I P P E D _ M A T R I X  \\\\\\\\\\\
!==================================================================
subroutine read_flipped_matrix( array, imax, jmax, infile )
!   Reads the content of a file writtin by WRITE_FLIPPED_MATRIX
!       into an array variable.  Primary is use is for
!       pfc parameters that vary by cell.
!
!   History:  28 June 2011  Original coding by BJE
!------------------------------------------------------------------
implicit none
! Arguments
integer, intent( in ) :: imax, jmax
real, dimension(imax,jmax), intent(inout) :: array
character(len=*), intent(in) :: infile
! Internal variables
integer :: i, j
integer :: jlist( jmax )
character(len=10) dummy
!--------------------------------------

OPEN( UNIT=9,   FILE = infile,  STATUS = 'OLD', ACTION = 'READ' ) 

read(9,*) dummy

do j = jmax, 1, -1
    read(9, *) jlist(j), array(:,j)
end do

close(9)

end subroutine
!===================================================================
!   \\\\\\\\\\       E N D      S U B R O U T I N E     ///////////
!   //////////   R E A D _ F L I P P E D _ M A T R I X  \\\\\\\\\\\
!==================================================================



!===================================================================
!   \\\\\\\\\\  B E G I N       S U B R O U T I N E     ///////////
!   //////////   R E A D _ B C s                        \\\\\\\\\\\
!==================================================================
subroutine READ_BCs( infile )
!   Reads boundary conditions that vary cell by cell from 
!       an output file.
use shared, only: north_bc, south_bc, east_bc, west_bc, &
                  imax, jmax
implicit none
! Arguments
character( len=*) :: infile
! Internal vars
integer :: i, j
character( len=5) dummy
integer ::    vj( imax )  !CV numbers for constant j
!-------------------------------------------------------------

OPEN( UNIT=9,   FILE = infile,  STATUS = 'OLD', ACTION = 'READ' ) 

read(9,*) dummy  
read(9,*) dummy
read(9,*) dummy, dummy, north_bc(:)

do j = jmax, 1, -1
    read(9,*) i, west_bc(j), vj(:), east_bc(j)
end do

read(9, *) dummy, dummy, south_bc(:)

CLOSE(9)

end subroutine READ_BCs
!===================================================================
!   \\\\\\\\\\       E N D      S U B R O U T I N E     ///////////
!   //////////   R E A D _ B C s                        \\\\\\\\\\\
!==================================================================

end module inputs
