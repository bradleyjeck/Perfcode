! fortran_free_source
!
!   This module is part of PERFCODE, written by Brad Eck
!
!   Purpose:    This module contains subroutines to output information

!============================================================================
!   \\\\\\\\\\                                        //////////
                        MODULE outputs 
!   //////////                                        \\\\\\\\\\
                        implicit none

                        contains
!============================================================================
!   1. ECHO_INPUTS (USE: shared)
!   2. WRITE_FLIPPED_MATRIX 
!   3. WRITE_MATRIX
!   4. WRITE_VECTOR
!   5. WRITE_SYSTEM
!   6. HEADER
!   7. WRITE_BCs_and_GRID
!   8. WRITE_CV_INFO
!   9. OUTPUTS_2D
!  10. WRITE_WATER_BUDGET

!===================================================================
!   \\\\\\\\\\       B E G I N    S U B R O U T I N E ///////////
!   //////////          E C H O _ I N P U T S         \\\\\\\\\\\
!==================================================================
!
!   Purpose:    This subroutine echos the input data to the 
!               specified device in comma seperated values format.
subroutine ECHO_INPUTS( dev )
use shared, only: K_input, por_input, b_pfc_input, n_mann_input, g
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
write(dev, 200) 'Hydraulic Conductivity (m/s),', K_input
write(dev, 200) 'Effective Porosity,', por_input
write(dev, 200) 'PFC Thickness (m),', b_pfc_input
write(dev, 200) "Manning's n,", n_mann_input
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

write(*,*) 'WRITE_FLIPPED_MATRIX: writing the file ', outfile

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
1       format( A, ',', 10000 ( I10, ',') )
2       format( I10, ',', 10000 ( F12.7, ',') )

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
1       format( A, ',', 10000 ( I10, ',') )
2       format( I10, ',', 10000 ( F14.7, ',') )
3       format(         10000 ( F14.7, ',') )
4       format(         10000 (   E14.7  , ',') )

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
!   Writes vector in usual form so the (1,1) entry
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

write(*,*) 'WRITE_VECTOR: writing the file ', outfile
 
OPEN( UNIT=9,   FILE = outfile,  STATUS = 'REPLACE' )  

! and the rest
do i = 1, imax 
        write(9, 3 ) array(i) 
end do 

close( 9 )

!---------------------------------------------------------------
! Format statements
3       format( ( E14.7, ',') )

!-------------------------------------------------------------------
end subroutine write_vector
!===================================================================
!   \\\\\\\\\\    E N D        S U B R O U T I N E    ///////////
!   //////////       W R I T E _ V E C T O R          \\\\\\\\\\\
!==================================================================


!===================================================================
!   \\\\\\\\\\    B E G I N      S U B R O U T I N E    ///////////
!   //////////       W R I T E _ S Y S T E M            \\\\\\\\\\\
!==================================================================
subroutine WRITE_SYSTEM( A, B, C, D, E, F, n, outfile )
!   Writes an output file of the equation system
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
3   format( (I10, ','),  6( E14.7, ',' ) )
!----------------------------------------------------------------------
end subroutine WRITE_SYSTEM
!===================================================================
!   \\\\\\\\\\     E N D      S U B R O U T I N E    ///////////
!   //////////       W R I T E _ S Y S T E M         \\\\\\\\\\\
!==================================================================


!===================================================================
!   \\\\\\\\\\    B E G I N      S U B R O U T I N E    ///////////
!   //////////            H E A D E R                   \\\\\\\\\\\
!==================================================================
subroutine HEADER( dev )
!   Writes a pretty header to the console
integer, intent (in ) :: dev    ! the device to write this on
!--------------------------------------------------------------------

write( dev, * ) ''
write( dev, * ) ''
write( dev, * ) ' PPP   EEEE  RRR   FFFF   CCC   OOO   DDD   EEEE '
write( dev, * ) ' P  P  E     R  R  F     C     O   O  D  D  E    '
write( dev, * ) ' P  P  E     R  R  F     C     O   O  D  D  E    '
write( dev, * ) ' PRP   EEEE  RRR   FFFF  C     O   O  D  D  EEEE '
write( dev, * ) ' P     E     R  R  F     C     O   O  D  D  E    '
write( dev, * ) ' P     E     R  R  F     C     O   O  D  D  E    '
write( dev, * ) ' P     EEEE  R  R  F      CCC   OOO   DDD   EEEE '
write( dev, * ) ''
write( dev, * ) '    PERmeable Friction COurse Drainage codE      '
write( dev, * ) ''
write( dev, * ) '        Written By: Brad Eck                     '
write( dev, * ) '                    bradleyjeck@gmail.com        '
write( dev, * ) ''
write( dev, * ) ''

!-------------------------------------------------------------------
end subroutine HEADER
!===================================================================
!   \\\\\\\\\\    E N D    S U B R O U T I N E    ///////////
!   //////////         H E A D E R                \\\\\\\\\\\
!==================================================================


!===================================================================
!   \\\\\\\\\\    B E G I N    S U B R O U T I N E   ///////////
!   //////////  W R I T E _ B C s _ A N D _ G R I D  \\\\\\\\\\\
!==================================================================
subroutine WRITE_BCs_AND_GRID(  outfile    )
!   Outputs a file that shows the boundary conditions along with
!       the grid numering scheme.  A companion routine reads the
!       file back into Perfcode to handle boundary conditions that
!       differ cell by cell. The grid is written in flipped form
!       as in WRITE_FLIIPED_MATRIX so that it appears consistent
!       with the physical geometry.
!
!    HISTORY:  30 June 2011 -- Original coding
!
use shared, only: north_bc, south_bc, east_bc, west_bc, &
                  imax, jmax
use pfc2Dfuns, only: F_LinearIndex
implicit none
! Arguments
character(len=*), intent(in) :: outfile
! Internal varieables
integer ::    vj( imax )  !CV numbers for constant j
integer :: ilist( imax )
integer :: i, j
character( len = imax*9 + 20 ) :: fmt1  ! character strings for formats 
character( len = imax*7 + 20 ) :: fmt2
character( len = imax*9 + 30 ) :: fmt3
!-------------------------------------------------------------------

! fortran i/o is messy.  Variables such as imax are not allowed
! in format statements.  However, formats are just character strings
! so we can write a long character string to do this instead.

! fmt1 <== 1 format(  2( A, ','), imax( I10, ',') )
fmt1 = "( A, ',', A, ',' "
do i = 1, imax   
    fmt1 = trim(fmt1) // "I10, ',' "
end do
fmt1 = trim(fmt1) // ")"

! fmt2 <== 2 format(  imax+2( A, ',') )
fmt2 = "("
do i = 1, imax+2
    fmt2 = trim( fmt2) // "A, ',' "
end do
fmt2 = trim(fmt2) // ")"

! fmt3 <== 3 format(  ( I10, ',' ), ( A, ','), imax( I10, ','), ( A, ',') )
fmt3 = "( I10, ',', A, ',' "
do i = 1, imax
    fmt3 = trim(fmt3) // "I10, ',' "
end do
fmt3 = trim(fmt3) // "A, ',' )"

! fill in ilist 
do i = 1, imax
    ilist(i) = i
end do

! With the format strings built, now write the actual file
open( unit = 22, file = outfile, status = 'REPLACE' )

write(22,*) 'Perfcode boundary conditions and grid numbering.'
write(22, fmt1) ' j \ i ', '  ', ilist(:)
write(22, fmt2) '   ', '  ', north_bc(:)

do j = jmax, 1, -1
    do i = 1, imax
        vj(i) = F_LinearIndex(i, j, jmax)
    end do

    write(22, fmt3) j, west_bc(j), vj(:), east_bc(j)

end do

write(22, fmt2) '  ', '  ', south_bc(:)

close( 22 )

! Screen message
write(*,*) 'WRITE_BCs_AND_GRID: Boundary conditions and grid numering'
write(*,*) '                    have been written to ', outfile
write(*,*)

end subroutine WRITE_BCs_AND_GRID
!===================================================================
!   \\\\\\\\\\    E N D    S U B R O U T I N E       ///////////
!   //////////  W R I T E _ B C s _ A N D _ G R I D  \\\\\\\\\\\
!==================================================================


!===================================================================
!   \\\\\\\\\\  B E G I N    S U B R O U T I N E  ///////////
!   //////////     W R I T E _ C V _ I N F O      \\\\\\\\\\\
!==================================================================
subroutine write_cv_info( )
use shared, only: cv_info, vmax, Z
implicit none
integer :: v

 open( unit=40, file = 'CV_info.csv', status = 'REPLACE')
  write(40,*) 'v,i,j,segment,xi,eta,X,Y,Z,'
  do v = 1, vmax
      WRITE(40,44) v, CV_info(v), Z( CV_info(v)%i, cv_info(v)%j )
  end do

  44  format( 4(I12, ','), 5(E14.7, ','))

end subroutine
!===================================================================
!   \\\\\\\\\\  E N D        S U B R O U T I N E  ///////////
!   //////////     W R I T E _ C V _ I N F O      \\\\\\\\\\\
!==================================================================


!===================================================================
!   \\\\\\\\\\  B E G I N    S U B R O U T I N E  ///////////
!   //////////     O U T P U T S _ 2 D            \\\\\\\\\\\
!==================================================================
subroutine OUTPUTS_2D(  )
use shared
use pfc2Dfuns, only: F_LinearIndex
!   Writes output files that are only appropriate to the 2D model
!-----------------------------------------------------------------

!-------------------------------------------------------------------------------
! Output grid numbering scheme to a file
! store grid numbering scheme and write it to a file
do j = 1, jmax
    do i = 1, imax
        grid( j, i) = F_LinearIndex( i, j, jmax )
    end do
end do

open( unit = 30, file = 'grid.csv', status = 'REPLACE' )
do j = jmax, 1, -1
    WRITE(30, 400 ) grid( j, : )
end do
close(30)

!-----------------------------------------------------------------------------------
! Output depth grid for last timestep

! an internal write statement to store the value of the REAL variable
! "time_simulated" in the CHARACTER variable "out_time"
write( out_time, 111 ) time_simulated

call write_flipped_matrix( h_old, imax, jmax, 'h_old'//out_time//' sec.csv' )

!------------------------------------------------------------------------------
! Output iteration history for the last time-step 

call write_matrix( h_temp_hist, vmax, qmax, 'h_temp_hist'//out_time//' sec.csv')


!-------------------------------------------------------------------------------
! file to show 1D solution along i = imax / 2; j = 1:jmax



OPEN(UNIT = 10, FILE = 'CrossSectionSoln.csv', STATUS='REPLACE')
WRITE(10,*) 'Output From PERFCODE.f95'
WRITE(10,*) 'Timestamp,', FILE_DATE,' ', FILE_TIME,','
do i = 1, 17
    write( 10, * ) input_variables(i), ',', input_values(i), ','
end do
write(10,*) 'north_bc,', north_bc       
write(10,*) 'south_bc,', south_bc 
write(10,*) 'east_bc,', east_bc
write(10,*) 'west_bc,', west_bc

WRITE(10,200) 'Average Rainfall Intensity (m/s),', &
               sum( rain(1:nlast )) / time_simulated
WRITE(10,200) 'Average Rainfall Intensity (cm/hr),',&
               sum( rain(1:nlast )) / time_simulated * 3600. * 100. 
WRITE(10,200) 'Final Time (sec),', time_simulated
WRITE(10,201) 'Number of cells longitudinally,', imax
WRITE(10,201) 'Number of cells transversly,'   , jmax
WRITE(10,201) 'Total Number of Grid Cells,', vmax
WRITE(10,200) 'CPU Time (seconds),', cputime
WRITE(10,200) 'Run Time (seconds),', &
               real(run_end_time - run_start_time)/real(count_rate)
WRITE(10, *) '************************** &
             & MODEL OUTPUT IN [ SI ] UNITS &
             &***********************************,'
i = imax / 2       
write(10,*) ' i = ', i,','
write(10, *) 'j,eta,Z,PFC_Surf,h,Head,Surf_Thk.mm,'
do j = 1, jmax
    v = F_LinearIndex( i, j, jmax)
    write(10, 2) j, CV_Info(v)%eta, Z(i,j), Z(i,j) + b_pfc(i,j),      &
                                h_old(i,j), Z(i,j) + h_old(i,j), &
                              ( h_old(i,j) - b_pfc(i,j) ) * 1000.
end do




!------------------------------------------------------------------------------------
! 3d plotting output for maximum depth 
!   ( contour plots of the resuls are made from  this file )
open( unit = 10, file = 'max_depth.csv', status = 'replace' )

write( 10, * ) 'v, X, Y, Z, h,'

do j = 1, jmax
    do i = 1, imax
        v = F_LinearIndex( i, j, jmax)
        write(10, 2) v, CV_Info (v ) % X,                    &
                        CV_Info( v ) % Y, Z(i,j), h_max(i,j)
    end do
end do

close( 10 )



!------------------------------------------------------------------------------------
! 3d plotting output for maximum discharge 
!   ( contour plots of the resuls are made from  this file )
open( unit = 10, file = 'max_Q.csv', status = 'replace' )

write( 10, * ) 'v, X, Y, Z, h,'

do j = 1, jmax
    do i = 1, imax
        v = F_LinearIndex( i, j, jmax)
        write(10, 2) v, CV_Info (v ) % X,                    &
                        CV_Info( v ) % Y, Z(i,j), h_Q_max(i,j)
    end do
end do

close( 10 )
!------------------------------------------------------------------------------------
! 3d plotting output for maximum mid-domain outlet depth 
!   ( contour plots of the resuls are made from  this file )
open( unit = 10, file = 'max_imidj1depth.csv', status = 'replace' )

write( 10, * ) 'v, X, Y, Z, h,'

do j = 1, jmax
    do i = 1, imax
        v = F_LinearIndex( i, j, jmax)
        write(10, 2) v, CV_Info (v ) % X,                    &
                        CV_Info( v ) % Y, Z(i,j), h_imid_j1_max(i,j)
    end do
end do

close( 10 )
!-----------------------------------------------------------------------------
! 3d plotting output for maximum mid-domain outlet depth 
!   ( contour plots of the resuls are made from  this file )
open( unit = 10, file = 'max_imiddepth.csv', status = 'replace' )

write( 10, * ) 'v, X, Y, Z, h,'

do j = 1, jmax
    do i = 1, imax
        v = F_LinearIndex( i, j, jmax)
        write(10, 2) v, CV_Info (v ) % X,                    &
                        CV_Info( v ) % Y, Z(i,j), h_imid_max(i,j)
    end do
end do

close( 10 )

!-----------------------------------------------------------------
!Format statements

2       FORMAT( I12, ',', 10000 ( E12.7, ',') )
10      FORMAT('   ', ( i3, '    '),  ( F10.3, '  ') , F10.6 ) 
111     FORMAT( f9.2 )
200     FORMAT ( A, ( E12.7, ',') ) 
201     FORMAT ( A, ( I12, ',') ) 
300     FORMAT ( 2 (     I12, ','),      F12.7, ',' , &    ! n, numit, maxdif
                         I12, ',' ,          E12.7, ',' , &    ! loc, L2_History
                 2 ( F12.8, ','), ( F12.3, ',' ), 3 ( F12.8 ,',')     )  ! rain, maxthk, time, Qout, h_imid_j1hist, h_imid_max_hist
400     FORMAT( 10000 ( I12, ',' ) )
401     FORMAT( (I12, ',') , 2( F12.7, ',' ) )
660     FORMAT(  2( I12, ','), 2( F12.7, ',') )

end subroutine OUTPUTS_2D
!===================================================================
!   \\\\\\\\\\   E N D    S U B R O U T I N E     ///////////
!   //////////     O U T P U T S _ 2 D            \\\\\\\\\\\
!==================================================================

!===================================================================
!   \\\\\\\\\\   B E G I N     S U B R O U T I N E     ///////////
!   //////////   W R I T E _ W A T E R _ B U D G E T   \\\\\\\\\\\
!==================================================================
Subroutine write_water_budget( nlast )
!   Writes an output file with water budget information.
!
!  History:  21 July 2011 -- Original coding
!
use shared, only: pfc_vol, surface_vol, rain_vol, north_vol, & 
                  south_vol, east_vol, west_vol, net_flow,    &
                   storage_change, vol_error, vol_error_frac, &
                  time
implicit none
! Arguments
integer, intent(in) :: nlast
! Internal variables
integer :: n
!------------------------------------------------------------------


open( unit = 9, file = 'WaterBudget.csv', status = 'replace')

write(9,*) 'PERFCODE GLOBAL WATER BUDGET,'
write(9,*) 'Volume output is in liters,'
write(9,*) ''
! Header line
write(9,*) 'Timestep, Time, pfc.vol.L, surface.vol.L, rain.vol.L, &
           & north.vol.L, south.vol.L, east.vol.L, west.vol.L, &
           & net.flow.L, storage.change.L, vol.error.L, &
           & vol.error.frac,'


do n = 1, nlast

    write(9, 92) n, time(n), pfc_vol(n) * 1000., &
                         surface_vol(n) * 1000., &
                            rain_vol(n) * 1000., &
                           north_vol(n) * 1000., &
                           south_vol(n) * 1000., &
                            east_vol(n) * 1000., &
                            west_vol(n) * 1000., &
                            net_flow(n) * 1000., &
                      storage_change(n) * 1000., &
                           vol_error(n) * 1000., &
                      vol_error_frac(n) * 1000. 
 
end do

close( 9 )

write(*,*) ''
write(*,*) 'WRITE_WATER_BUDGET:  WaterBudget.csv written to disk'
write(*,*) ''

! Format statements
92 format(  (I12  , ',' ), &   !time step
            (F14.7, ',' ), &   !time
          11(E14.7, ',' )    ) !all other values


!-------------------------------------------------------------------
end subroutine
!===================================================================
!   \\\\\\\\\\   E N D    S U B R O U T I N E          ///////////
!   //////////   W R I T E _ W A T E R _ B U D G E T   \\\\\\\\\\\
!==================================================================





END MODULE outputs
