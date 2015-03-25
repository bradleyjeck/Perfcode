! fortran_free_source
!
!  This module is part of PERFCODE, written by Bradley J. Eck.
!
!  Contents:  
!       1. Generate_Grid  Computes the length, width, and area 
!                          of each grid cell.
!       2. SET_ELEVATIONS  Gives an elevation to each grid cell center
!
MODULE GridGen 
implicit none
contains
!============================================================================
!   \\\\\\\\\        B E G I N      S U B R O U T I N E       ///////////
!   /////////            G E N E R A T E  _ G R I D           \\\\\\\\\\\\
!============================================================================
subroutine Generate_Grid( prelim_dx, prelim_dy )
!   This subroutine takes the centerline geometry file that is generated
!   mannualy and creates a curvilinear grid.
!
!       INPUTS:  Preliminary grid spacing
!       OUTPUTS: Size of computational domain (imax & jmax )
!                Length, width and area of each grid cell ( module SHARED)
!                Coordinates of each CV center 
!
USE shared,     ONLY: seg, area, lng, lng_south, lng_north, wid,   & 
                      imax, jmax, vmax, CV_Info, astat, gridcell, model_mode
USE outputs,    ONLY: WRITE_MATRIX, WRITE_FLIPPED_MATRIX, WRITE_CV_INFO
USE geom_funcs, ONLY: F_L_xi, UNMAP_X, UNMAP_Y
USE pfc2Dfuns,  ONLY: F_LinearIndex
!------------------------------------------------------------------------------
implicit none
! VARIABLE DECLARATIONS
!   Agruments
real, intent(in) :: prelim_dx, prelim_dy
!   Internal Variables
INTEGER, parameter :: max_rec=144 ! maximum allowable number of rainfall records
integer :: N_seg   
CHARACTER(len=20) infile    ! the file name to read parameters from
integer :: i, j, v     ! looping variables
character :: TRASH
real :: xi, eta, X, Y
real :: eta_s, eta_n
integer :: N_xi, N_eta, Seg_num
! for writing border info to file
integer :: factor = 5         ! how many times more border points than CVs?
integer :: ijf ! dummy variable for i or j times the factor
real, allocatable, dimension(:) :: NX, NY, SX, SY, EX, EY, WX, WY
! for writing grid to a file
integer :: res = 4
real, allocatable, dimension(:,:) :: X_gl_long, Y_gl_long, X_gl_tran, Y_gl_tran
!----------------------------------------------------------------------------------------------
! READ IN THE GEOMETRY DATA

! default value for input file
infile = 'CL_Segments.dat'

! Prompt the user for the input file
WRITE(*,*) 'GENERATE_GRID: Enter filename or press / for ', infile
READ(*,*) infile

OPEN( UNIT=8,   FILE = infile, ACTION = 'read', STATUS = 'old' )  

! Number of segments
read( unit=8, fmt = * ) N_seg
read( unit=8, fmt = * ) trash

! allocate enough space to read in the CL segments
allocate( seg( N_seg ) )

do i = 1, N_seg
    READ( unit = 8 , fmt = * ) j, seg(j)%xcc1, seg(j)%ycc1, seg(j)%dx, &
                                  seg(j)%dy,   seg(j)%R1,   seg(j)%dR, &
                                  seg(j)%W,    seg(j)%theta1, seg(j)%dtheta
end do

close( 8 )

! estimate the length of each segment by evaluating the metric coefficients
do i = 1, N_seg
    seg(i)%arclen = F_L_xi( 0.5, 0.5, seg(i) ) * 1.0   ! L_xi * Delta_xi 
end do


print *, 'Segment    ArcLength' 
do i = 1, N_seg
    print *, i, seg(i)%arclen
end do

!total length
print *, '  Total Length: ', sum( seg(:)%arclen )

! average length
print *, 'Average Length: ', sum( seg(:)%arclen ) / real( N_seg ) 


! Number of elements per segment
!nint = nearest integer
N_xi = nint(  sum( seg(:)%arclen)  / real( N_seg ) / prelim_dx )   
N_eta= nint(  sum( seg(:)%W     )  / real( N_seg ) / prelim_dy )  
print *, 'N_xi = ', N_xi, '    N_eta = ', N_eta
 
! size of computational domain
imax = N_xi * N_seg
jmax = N_eta
vmax = imax * jmax


!NOTE: the rest of this subroutine only applies if
!      the model mode is 2D
if2d: if( model_mode .eq. '2D' ) then 

  !--------------------------------------------------------
  ! ALLOCATE ARRAYS

  allocate(    lng( imax, jmax),  STAT = astat( 1) )
  allocate(    wid( imax, jmax),  STAT = astat( 2) )
  allocate(   area( imax, jmax),  STAT = astat( 3) ) 


  allocate(   lng_south( imax, jmax ),    STAT = astat( 6) )
  allocate(   lng_north( imax, jmax ),    STAT = astat( 7) )

  allocate ( CV_Info( vmax ), STAT = astat(17) ) 
  !---------------------------------------------------------------

  ! Now compute length, width, and area of each cell
  ! should confirm that all widths are the same

  do j = 1, jmax
  !    print *, 'j = ', j
  !    print *, 'i   Segment'
      do i = 1, imax
          ! Determine which segment we're in
          ! the intrinic function CEILING is like ROUNDUP in excel
          Seg_Num = ceiling( real(i) / real( imax ) * real( N_seg ) )
          ! Compute values of xi for the cell that we're in
          if ( i .LE. N_xi ) then
                  xi =   i * 1. / N_xi - 1. / N_xi / 2.
          else
                  xi = ( i - ( Seg_Num - 1 ) * N_xi ) * 1. / N_xi - 1. / N_xi / 2.
          end if
          ! value of eta
          eta = 1. / N_eta * j  -  1. / N_eta / 2.
          ! Physical Coordinates of CV
          X = unmap_x( xi = xi, eta = eta, seg = seg( Seg_Num ) )
          Y = unmap_y( xi = xi, eta = eta, seg = seg( Seg_Num ) ) 
          ! store the summary information for this cell
          v = F_LinearIndex( i, j, jmax ) 
          CV_Info(v) = gridcell( i, j, Seg_Num, xi, eta, X, Y ) 
          ! now compute the quantities of interest for each cell.  
  !        print *, i, Seg_Num
          lng (i,j) = F_L_xi( xi, eta, seg( Seg_Num ) ) * 1. / N_xi
          wid (i,j) = seg( Seg_Num)%W / N_eta
          area(i,j) = lng(i,j) * wid(i,j)
          ! compute lengths for north and south faces of cell
          !   south
          eta_s = eta - 1./N_eta * 1./2.
          lng_south(i,j) = F_L_xi( xi, eta_s, seg( Seg_Num) ) * 1./N_xi
          !   north
          eta_n = eta + 1./N_eta * 1./2.
          lng_north(i,j) = F_L_xi( xi, eta_n, seg( Seg_Num) ) * 1./N_xi
      end do
  end do

 


  !----------------------------------------------------------------------------
  !   >>>>>>>>>>     W R I T E   B O U N D A R Y   C O O R D S    <<<<<<<<<<
  !----------------------------------------------------------------------------
  ! was going to make this a subroutine, but it seemed easier to add it here

  ! Allocate arrays
  ijf = imax * factor
  allocate( NX( ijf) )
  allocate( NY( ijf) )
  allocate( SX( ijf) )
  allocate( SY( ijf) )

  ijf = jmax * factor + 1
  allocate( EX( ijf ) )
  allocate( EY( ijf ) )
  allocate( WX( ijf ) )
  allocate( WY( ijf ) )

  ! NORTH and SOUTH borders
  !re-calc N_xi so as not to change the following formula
  N_xi = N_xi * factor
  do i = 1, imax * factor  
      Seg_Num = ceiling( real(i) / real( imax * factor ) * real( N_seg ) )
      ! Compute values of xi 
      if ( i .LE. N_xi ) then
              xi =   i * 1. / N_xi - 1. / N_xi / 2.
      else
              xi = ( i - ( Seg_Num - 1 ) * N_xi ) * 1. / N_xi - 1. / N_xi / 2.
      end if
      ! NORTH -- Physical Coordinates on border 
      eta = 1.0
      NX(i) = unmap_x( xi = xi, eta = eta, seg = seg( Seg_Num ) )
      NY(i) = unmap_y( xi = xi, eta = eta, seg = seg( Seg_Num ) ) 
      ! SOUTH -- Physical Coordinates on border
      eta = 0.0
      SX(i) = unmap_x( xi = xi, eta = eta, seg = seg( Seg_Num ) )
      SY(i) = unmap_y( xi = xi, eta = eta, seg = seg( Seg_Num ) )    
  end do
      
  ! EAST and WEST borders
  do j = 1, jmax * factor + 1
      eta = ( j - 1. ) / ( jmax * factor ) 
      ! WEST
      xi = 0.0
      WX(j) = unmap_x( xi = xi, eta = eta, seg = seg( 1 ) )
      WY(j) = unmap_y( xi = xi, eta = eta, seg = seg( 1 ) )
      ! EAST
      xi = 1.0
      EX(j) = unmap_x( xi = xi, eta = eta, seg = seg( N_seg ) )
      EY(j) = unmap_Y( xi = xi, eta = eta, seg = seg( N_seg ) )    
  end do

  ! Write the borders to files
  open( unit = 40, file = 'NS_borders.csv', status = 'REPLACE')
  WRITE(40, *) 'NX, NY, SX, SY,'
  do i = 1, imax * factor
      write(40, 4) NX(i), NY(i), SX(i), SY(i)
  end do
  close(40)
       
  open( unit = 41, file = 'EW_borders.csv', status = 'REPLACE')
  write( 41, * ) 'EX, EY, WX, WY,'
  do j = 1, jmax * factor + 1
      write(41, 4) EX(j), EY(j), WX(j), WY(j)
  end do
  close(41)

  4 format( 4 (E12.7, ',') ) 


  !----------------------------------------------------------------------------
  !   >>>>>>>>>   WRITE GRID FOR PLOTTING                     <<<<<<<<<
  !----------------------------------------------------------------------------
  !
  ! envisioning the use of R's MATPLOT command, write a matrix of X coords
  ! and a matrix of Y coords
  ! we plot a bunch of points

  res = 4   !parameter to control the resolution of of the plotting
             ! this is analogous to the variable 'factor'  how many
             ! points do you want per CV?

  allocate( X_gl_long( imax * res, jmax + 1 ) )
  allocate( Y_gl_long( imax * res, jmax + 1 ) )

  ! redefine N_xi to reflect the amplified number of points for the grid
  N_xi = imax / N_seg * res 

  ! LONGITUDINAL GRID LINES
  do j = 1, jmax + 1
      !value of eta is constant for each j
      eta = ( j - 1. ) / jmax 
      do i = 1,  imax * res 
          ! figure out which segment we're in
          Seg_Num = ceiling( real(i) / real( imax * res ) * real( N_seg ) )
          ! Compute values of xi 
          if ( i .LE. N_xi ) then
              xi =  ( i - 1. ) / ( N_xi )
          else
              xi = ( i - ( ( Seg_Num - 1. ) * N_xi ) ) * 1. / N_xi
          end if
          X_gl_long(i, j)  = unmap_x( xi = xi, eta = eta, seg = seg( Seg_Num ) )
          Y_gl_long(i, j)  = unmap_y( xi = xi, eta = eta, seg = seg( Seg_Num ) )
      end do
  end do


  call WRITE_MATRIX( X_gl_long, imax * res  , jmax + 1, 'X_gl_long.csv' )
  call WRITE_MATRIX( Y_gl_long, imax * res  , jmax + 1, 'Y_gl_long.csv' )


  !TRANSVERSE GRID LINES (these are just straight and so only require two points)
  allocate( X_gl_tran( 2, imax + 1 ) )
  allocate( Y_gl_tran( 2, imax + 1 ) )
     
  ! put N_xi back to what its proper value
  N_xi = imax / N_seg

  do j = 1, 2
      eta = j - 1.
      do i = 1, imax + 1
          ! figure out which segment we're in
          Seg_Num = ceiling( real(i) / real( imax + 1 ) * real( N_seg ) )
          ! Compute values of xi 
          if ( i .LE. N_xi ) then
              xi =  ( i - 1. ) / ( N_xi )
          else
              xi = ( i - 1.  - ( Seg_Num - 1. ) * N_xi ) / N_xi
          end if
          X_gl_tran(j, i)  = unmap_x( xi = xi, eta = eta, seg = seg( Seg_Num ) )
          Y_gl_tran(j, i)  = unmap_y( xi = xi, eta = eta, seg = seg( Seg_Num ) ) 
      end do
  end do


  call WRITE_MATRIX( X_gl_tran, 2, imax + 1 , 'X_gl_tran.csv' )
  call WRITE_MATRIX( Y_gl_tran, 2, imax + 1, 'Y_gl_tran.csv' )



!----------------------------------------------------------------------------
endif if2d 
!---------------------------------------------------------------------------
end subroutine Generate_Grid
!============================================================================
!   \\\\\\\\\   E N D         S U B R O U T I N E   ///////////
!   /////////       G E N E R A T E _ G R I D       \\\\\\\\\\\\
!===========================================================================


!============================================================================
!   \\\\\\\\\   B E G I N     S U B R O U T I N E   ///////////
!   /////////       S E T _ E L E V A T I O N S     \\\\\\\\\\\\
!===========================================================================
subroutine Set_Elevations( ) 
!   Gives an elevation to each node of the grid
!
use SHARED, only: Z, imax, jmax, lng_south, CV_Info, seg,  &
                  nr_cs, slope_cs, wid_cs, eta_cs, Z_cs,   &
                  nr_lp, dist_lp, Z_lp, model_mode
use inputs, only: GET_LONG_PROFILE, GET_CROSS_SECTION 
use utilities, only: F_Linterp
use outputs,   only: write_flipped_matrix
use pfc2dfuns,  only: F_LinearIndex
implicit none
INTEGER :: i, j, v      ! looping variables
real :: tot_wid  ! total width of the cross section
real :: CL_wid
! Generating elevations
real :: dist_along_lp, Z_eta_0, Z_add

!------------------------------------------------------------------------
! C R O S S   S E C T I O N
!------------------------------------------------------------------------

call GET_CROSS_SECTION( )

! Compute eta and elevation from widths and slopes
!    Given:  widths and slopes  slope_cs   wid_cs, nr_cs
!    Find : Z vs eta
tot_wid = sum( wid_cs(1:nr_cs) )

CL_wid = seg(1) % W   

! Check tot_wid for consistency with the width given in CL segments
if(    abs(  tot_wid  -  CL_wid )   .GE. 1.e-3 ) then
    write( *,*) ' Cross Section Width =', tot_wid
    write( *,*) ' Centerline Width = ', CL_wid
    write(*,*) ' SET_ELEVATIONS: Total width specified in Cross Section file'
    write(*,*) '                 is inconsistent with the centerine geometry...Stopping Program'
    STOP

end if

eta_cs(1) = 0.
  z_cs(1) = 0.   !<---dummy value here, elevations are made relative to eta=0

! Compute etas and elevations 
do i = 2, nr_cs + 1
    eta_cs(i) = eta_cs(i-1) +   wid_cs(i-1) / tot_wid
      Z_cs(i) =   z_cs(i-1) - slope_cs(i-1) * wid_cs(i-1)
end do

! print the results to confirm
print *, ' CROSS SECTION POINTS '
WRITE(*,*) ' Point        Eta       Elevation '
WRITE(*,*)  '=================================='
do i = 1, nr_cs + 1
    write(*,10) i, eta_cs(i), Z_cs(i)
end do


!------------------------------------------------------------------------
!     L O N G I T D I N A L   P R O F I L E 
!------------------------------------------------------------------------
call GET_LONG_PROFILE


!--------------------------------------------------------------------------------
! I N T E R P O L A T E   E L E V A T I O N S  
!--------------------------------------------------------------------------------

!NOTE: the rest of this subroutine only applies if
!      the model mode is 2D
if2d: if( model_mode .eq. '2D' ) then 

    allocate(         Z( imax, jmax ) );           Z = 0.0  !,  STAT = astat( 4)

    do i = 1, imax
        ! Compute distance along longitudinal profile at eta = 0
        dist_along_lp = sum( lng_south(1:i,1 ) ) - lng_south(i,1)/2.
        ! Compute elevation at eta = 0 for this column of grid cells
        Z_eta_0 = F_Linterp(      x = dist_along_lp  ,   &
                            Known_X = dist_lp        ,   &
                            Known_Y =    Z_lp        ,   &
                                  n =   nr_lp               )
        do j = 1, jmax
            v = F_LinearIndex( i, j , jmax )
            Z_add  =        F_Linterp(       X = CV_Info(v)%eta , &
                                       known_X = eta_cs         , &
                                       known_Y = Z_cs           , & 
                                             n = (nr_cs + 1)        )
            Z(i,j) = Z_eta_0  + Z_add
        end do
    end do

endif if2d

!--------------------------------------------------------------------------------
! F O R M A T   S T A T E M E N T S 
!--------------------------------------------------------------------------------

10     FORMAT('   ', ( i3, '    '),  ( F10.3, '  ') , F10.6, i12   )


!---------------------------------------------------------------------------
end subroutine Set_Elevations 
!============================================================================
!   \\\\\\\\\   E N D      S U B R O U T I N E       ///////////
!   /////////           S E T _ E L E V A T I O N S  \\\\\\\\\\\\
!===========================================================================


END MODULE GridGen









