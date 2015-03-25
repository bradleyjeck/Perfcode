! fortran_free_source
!
! (c) Copyright 2010, 2015 Bradley J. Eck
! This module is part of PERFCODE 
!

!   Purpose:    This module contains subroutines to output information

!============================================================================
!   \\\\\\\\\\                                        //////////
                        MODULE outputs 
!   //////////                                        \\\\\\\\\\
                        implicit none

                        contains
!============================================================================
!   1. ECHO_INPUTS
!   2. WRITE_FLIPPED_MATRIX
!   3. WRITE_MATRIX
!   4. WRITE_VECTOR
!   5. WRITE_SYSTEM

!===================================================================
!   \\\\\\\\\\       B E G I N    S U B R O U T I N E ///////////
!   //////////          E C H O _ I N P U T S         \\\\\\\\\\\
!==================================================================
!
!   Purpose:    This subroutine echos the input data to the 
!               specified device in comma seperated values format.
subroutine ECHO_INPUTS( dev )
use shared, only: K, por, b_pfc, n_mann, g
!---------------------------------------------------------------------
! VARIABLE DECLARATIONS
!   Arguments
integer, intent( In ) :: dev    ! The device number that the output 
!integer, intent( in ) :: nrr    ! number of rainfall records
!REAL, intent( in ) :: K, por, b_pfc, n_mann, g
!REAL, intent( in ) :: rain_time(:), rain_rate(:)
!real, intent( in ) :: dt
!   Internal variables
!integer :: i

!----------------------------------------------------------------------
write(dev, *  ) 'SUMMARY OF INPUT DATA,'
write(dev, 200) 'Hydraulic Conductivity (m/s),', K
write(dev, 200) 'Effective Porosity,', por
write(dev, 200) 'PFC Thickness (m),', b_pfc
write(dev, 200) "Manning's n,", n_mann
write(dev, 200) 'Gravitational Acceleration (m/s/s),', g

!--------------------------------------------------------------------------------
! Format Statements
200     FORMAT (' ', A, ( F10.6, ',') )

!------------------------------------------------------------------
end subroutine ECHO_INPUTS
!===================================================================
!   \\\\\\\\\\      E N D         S U B R O U T I N E ///////////
!   //////////          E C H O _ I N P U T S         \\\\\\\\\\\
!==================================================================

!===================================================================
!   \\\\\\\\\\    B E G I N    S U B R O U T I N E     ///////////
!   ////////// W R I T E _ F L I P P E D _ M A T R I X \\\\\\\\\\\
!==================================================================
subroutine write_flipped_matrix( array, imax, jmax, outfile )
!   Writes matrix in 'flipped' form so it corresponds to
!       the physical geometry.  This means the (1,1) entry
!       appears at the bottom left corner of the ouput file.
!-----------------------------------------------------------
! VARIABLE DECLARATIONS
!   Arguments
integer imax, jmax
real array( imax, jmax)
character(len=*):: outfile  !assumed length specifier *
!   Internal variables
integer :: i, j
integer :: ilist( imax ), jlist( jmax )
!-------------------------------------------------------
! Create lists of indices
do i = 1, imax
    ilist(i) = i
end do

do j = 1,  jmax
    jlist(j) = j
end do

print *, 'WRITE_FLIPPED_MATRIX: writing the file ', outfile

! Write the length array in upside down form so it corresponds to
!   to the physical geometry.
OPEN( UNIT=9,   FILE = outfile,  STATUS = 'REPLACE' )  

! First line
write(9, 1) ' j \ i ', ilist(:) 

! and the rest
do j = jmax, 1, -1
        write(9, 2) jlist( j ), array(:,j) 
end do 

close( 9 )

!---------------------------------------------------------------
! Format statements
1       format( A, ',', 10000 ( I, ',') )
2       format( I, ',', 10000 ( F12.7, ',') )

!-------------------------------------------------------------------
end subroutine write_flipped_matrix
!===================================================================
!   \\\\\\\\\\    E N D        S U B R O U T I N E    ///////////
!   ////////// W R I T E _ F L I P P E D _ M A T R I X \\\\\\\\\\\
!==================================================================


!===================================================================
!   \\\\\\\\\\    B E G I N    S U B R O U T I N E     ///////////
!   //////////        W R I T E _ M A T R I X          \\\\\\\\\\\
!==================================================================
subroutine write_matrix( array, imax, jmax, outfile )
!   Writes matrix in usual form so the (1,1) entry
!       appears at the top left corner of the ouput file.
!-----------------------------------------------------------
! VARIABLE DECLARATIONS
!   Arguments
integer imax, jmax
real array( imax, jmax)
character(len=*):: outfile  !assumed length specifier *
!   Internal variables
integer :: i, j
integer :: ilist( imax ), jlist( jmax )
!-------------------------------------------------------
! Create lists of indices
do i = 1, imax
    ilist(i) = i
end do

do j = 1,  jmax
    jlist(j) = j
end do

print *, 'WRITE_MATRIX: writing the file ', outfile

! 
OPEN( UNIT=9,   FILE = outfile,  STATUS = 'REPLACE' )  


do i = 1, imax 
         write(9, 4) array(i,:)
end do 

close( 9 )

!---------------------------------------------------------------
! Format statements
1       format( A, ',', 10000 ( I, ',') )
2       format( I, ',', 10000 ( F12.7, ',') )
3       format(         10000 ( F12.7, ',') )
4       format(         10000 (   E  , ',') )

!-------------------------------------------------------------------
end subroutine write_matrix
!===================================================================
!   \\\\\\\\\\    E N D        S U B R O U T I N E    ///////////
!   //////////       W R I T E _ M A T R I X          \\\\\\\\\\\
!==================================================================

!===================================================================
!   \\\\\\\\\\    B E G I N    S U B R O U T I N E     ///////////
!   //////////        W R I T E _ V E C T O R          \\\\\\\\\\\
!==================================================================
subroutine write_vector( array, imax,  outfile )
!   Writes matrix in usual form so the (1,1) entry
!       appears at the top left corner of the ouput file.
!-----------------------------------------------------------
! VARIABLE DECLARATIONS
!   Arguments
integer imax 
real array( imax )
character(len=*):: outfile  !assumed length specifier *
!   Internal variables
integer :: i 
!-------------------------------------------------------
! Create lists of indices

print *, 'WRITE_VECTOR: writing the file ', outfile


! 
OPEN( UNIT=9,   FILE = outfile,  STATUS = 'REPLACE' )  

! and the rest
do i = 1, imax 
        write(9, 3 ) array(i) 
end do 

close( 9 )

!---------------------------------------------------------------
! Format statements
3       format( ( E, ',') )

!-------------------------------------------------------------------
end subroutine write_vector
!===================================================================
!   \\\\\\\\\\    E N D        S U B R O U T I N E    ///////////
!   //////////       W R I T E _ V E C T O R          \\\\\\\\\\\
!==================================================================


!======================================================================
subroutine WRITE_SYSTEM( A, B, C, D, E, F, n, outfile )
integer, intent( in ) :: n
real, dimension( n ), intent( in ) :: A, B, C, D, E, F
character( len=*) :: outfile

integer :: i


open( unit = 11, file = outfile, status = 'REPLACE' )
write( 11, *) 'v, A, B, C, D, E, F,'
do i = 1, n
    write( 11, 3 ) i, A(i), B(i), C(i), D(i), E(i), F(i)
end do
close( 11 )


!----------------------------------------------------------------------
! Format statements
3   format( (I, ','),  6( E, ',' ) )
!----------------------------------------------------------------------
end subroutine WRITE_SYSTEM
!======================================================================


!============================================================================
!   \\\\\\\\\                                       ///////////
                       END MODULE outputs
!   /////////                                       \\\\\\\\\\\\
!============================================================================
