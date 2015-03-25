! fortran_free_source
!
! (c) Copyright 2010, 2015 Bradley J. Eck
! This module is part of PERFCODE 
!

module pfc1Dsubs

implicit none

contains





!===================================================================================
!       \\\\\\\\\\      B E G I N    S U B R O U T I N E      //////////
!       //////////     S E T U P    1 D   S E C T I O N       \\\\\\\\\\
!===================================================================================
SUBROUTINE setup_1D_section( )
!   does the setup work for looking at this as a 1D section
USE shared,    ONLY: Z_lp, nr_lp, dist_lp, slope_cs_1D,        &
                     wid_cs_1d, eta_cs_1d, nr_cs, long_slope,  &
                     slope_cs, wid_cs
USE utilities, ONLY: F_PythagSum
!------------------------------------------------------------------
integer :: i
!------------------------------------------------------------------
! Compute longitudinal slope 
!  ( assumed to be constant thorought the domain) 
long_slope = (    Z_lp( nr_lp ) -    Z_lp( 1 ) ) / &
             ( dist_lp( nr_lp ) - dist_lp( 1 ) )


! Using the longitudinal slope and cross slope,
! compute the slopes and segment widths for the 1D profile

allocate( slope_cs_1d( nr_cs) )
allocate(   wid_cs_1d( nr_cs) )
allocate(   eta_cs_1d( nr_cs+1) )

write(*,*) ''
write(*,*) '  1D CROSS SECTION   '
write(*,*) ' Segment   Slope   Width '
write(*,*) '============================'
do i = 1, nr_cs
    !For the slope, compute the magnitude of the resultant slope
    ! using pythagorean sum and then use the intrinsic SIGN function
    ! to give the resultant the same sign as the cross slope.
    slope_cs_1D(i) = SIGN (  F_PythagSum( long_slope, slope_cs(i) ) , &
                                                        slope_cs(i)       )
    ! Compute 1D width using similar triangles
    wid_cs_1D(i) = slope_cs_1D(i) / slope_cs(i) * wid_cs(i)  
    ! Print the results as we go
    write(*,10) i, slope_cs_1D(i), wid_cs_1d(i)
end do
write(*,*) ''



eta_cs_1d(1) = 0.

! Compute etas and elevations 
do i = 2, nr_cs + 1
    eta_cs_1d(i) = eta_cs_1d(i-1) +   wid_cs_1d(i-1) / sum( wid_cs_1D(1:nr_cs) )
end do

! print the results to confirm
print *, ' 1D CROSS SECTION POINTS '
WRITE(*,*) ' Point        Eta        '
WRITE(*,*)  '=================================='
do i = 1, nr_cs + 1
    write(*,10) i, eta_cs_1d(i)
end do


!--------------------------------------------------------------------------------
! Format Statements
10      FORMAT('   ', ( i3, '    '),  ( F10.3, '  ') , F10.6 )

!--------------------------------------------------------------------------------
end subroutine setup_1D_section
!===================================================================================
!       \\\\\\\\\\      E N D        S U B R O U T I N E   //////////
!       //////////     S E T U P  1 D  S E C T I O N       \\\\\\\\\\
!===================================================================================

!===================================================================================
!       \\\\\\\\\\      B E G I N    S U B R O U T I N E   //////////
!       //////////     G R I D   1 D   S E C T I O N       \\\\\\\\\\
!===================================================================================
SUBROUTINE grid_1d_section( slope_in, width_in, seg, dx )
!---------------------------------------------------------------------------
use shared, only: TNE, XCV, ZCV, EDX, etaCV
use pfc1Dfuns
use utilities, only: F_Linterp
use outputs, only: WRITE_MATRIX

!Define variables
    
    IMPLICIT NONE

    !CONSTANTS
    INTEGER, intent(in) :: seg
    real, intent(in) :: dx

    !ARRAYS
    REAL, dimension( seg ), intent(in) :: slope_in,  width_in 


   
!calculation variables
    real, dimension(seg) :: slope, width
    INTEGER, dimension( seg ) :: ne, ir
    
    INTEGER :: gb

    INTEGER :: i, n, s, start, finish

    REAL, ALLOCATABLE :: xface(:), zface(:)
    real, allocatable :: seg_X(:), seg_Z(:)

    REAL :: DX1




!---------------------------------------------------
! Need to reverse slope and width arrays based on
!   the design of this subroutine.
! Indices for Reverse arrays (uses an implied DO loop )
ir = (/ ( i, i = seg, 1, -1  ) /)


    do i = 1, seg
        slope(i) = slope_in( ir(i) )
        width(i) = width_in( ir(i) )
           ne(i) = NINT( width(i) / dx )
    end do

!-------------------------------------------------------------
!Compute derivative quanties   & allocate remaining arrays  

    gb = seg - 1    
    TNE = sum(ne) + gb

    allocate( XFACE(TNE+1), &
              ZFACE(TNE+1), &
                XCV(TNE)  , &
                ZCV(TNE)  , &
                EDX(TNE)  , &
              etaCV(TNE)         )

!---------------------------------------------------------------

! compute the points for the boundaries and the CV centers

    XFACE( 1 ) = 0. !could use a different starting point
    EDX(:) = dx     ! all elemnts are the same size

    do i = 1, TNE
        XFACE( i+1 ) = XFACE( i ) + EDX( i )
        XCV  ( i   ) = (  XFACE( i ) + XFACE( i+1 )  ) / 2. 
        etaCV( i   ) = 1. - XCV( i ) / sum( width(:) )
    end do
!------------------------------------------------------------------------
! interpolate elevations of the points
    
    allocate ( seg_X(seg+1), seg_Z(seg+1) )

    seg_x(1) = 0.
    seg_Z(1) = 10.

    do i = 1, seg
        seg_x(i+1) = seg_X(i) + width(i)
        seg_Z(i+1) = seg_Z(i) + width(i) * slope(i) 
    end do

    ! first cross section by interpolation
    ZFACE( 1 ) = seg_z( 1 ) 
    do i = 1, TNE
        ZFACE( i+1 ) =  F_linterp( XFACE( i+1 ) , seg_X, seg_Z, seg+1 )
        ZCV  ( i   ) =  F_linterp( XCV  ( i   ) , seg_X, seg_Z, seg+1 )
    end do
    
!------------------------------------------------------------------------
! output the resulting arrays
    
! VECTOR FORM
    OPEN( UNIT = 20, FILE = 'grid_1D_section.csv', STATUS = 'REPLACE' )
    WRITE(20,*) 'XFACE, ZFACE, CV, XCV, ZCV, EDX, etaCV'
    do i = 1, TNE
        WRITE(20,100) XFACE( i ), ZFACE( i ), i, XCV( i ), ZCV( i ), EDX(i), etaCV(i)
    end do
    WRITE(20,99) XFACE( TNE+1 ), ZFACE( TNE+1)
    close(20)

!--------------------------------------------------
!Format statements
90      FORMAT( i, F12.6 )
99      FORMAT( 2( F12.6, ',' ) )
100     FORMAT( 2( F12.6, ',' ), 1(I, ','), 5( F12.6, ',')  )


!------------------------------------------------------------------------------
    end subroutine GRID_1D_SECTION
!===================================================================================
!       \\\\\\\\\\          E N D   S U B R O U T I N E     //////////
!       //////////       G R I D  1 D  S E C T I O N        \\\\\\\\\\
!===================================================================================



!====================================================================================
!	\\\\\\\\\\\\        B E G I N      P R O G R A M     //////////
!	////////////            P F C 1 D I M P              \\\\\\\\\\
!====================================================================================
!       Purpose:        This program computes a 1D solution for unsteady
!                       drainage through a PFC.  The water THICKNESS in each
!                       cell is used as the primary variable.  
!       History:        Source code revised from previous program that used
!                       an explicit method, and revised again to use as a 
!                       subroutien within the 2D model
!       IC:             The depth on input to the subroutine
!       BCs:            Various
!       Linearization:  Picard Iteration (lag the coefficients)
!       Linear Solver:  Thomas Alogorithm used for this 1D case
!       
!       Externals:      1.
!                       2.


SUBROUTINE pfc1Dimp( h_old, dt, rain, tolit, qmax, &
                     h_new, imax, eta_0_BC, eta_1_BC )
!bring in related modules
use shared,    only: K, g, b_pfc, n_mann, por, XCV, ZCV, EDX, &
                     etaCV, relax, relax_tran, h_pfc_min, eta_0_hp2_max
use utilities, only: F_Linterp, F_L2_NORM
use solvers,   only: Thomas
use pfc1Dfuns
use outputs,   only: WRITE_MATRIX, write_vector

!-----------------------------------------------------------------------------------
!VARIABLE DECLARATIONS

implicit none
! NOTE:  Because of ONLY statement for modules, only the names variables
!        are used in this subroutine.  This allows re-using variables names
!        thus the variable 'h_old' in this subroutine is not the same as 
!        the 'h_old' matrix in PERFCODE.

! Arguments
integer,              intent(in) :: imax    ! imax == TNE from 1D gridding
real, dimension(imax), intent(in) :: h_old
real,                   intent(in) :: dt
real,                   intent(in) :: rain
real,                   intent(in) :: tolit
integer,                intent(in) :: qmax
real, dimension(imax), intent(out):: h_new
character( len = 7 ) :: eta_0_BC, eta_1_BC

!SCALERS
INTEGER :: i, n, q
REAL    :: PF, PF1, Cw, Ce, Ce1, Cw1
REAL    :: hp, hs
REAL    :: cputime
real, dimension(:), pointer ::  X, Z, DX 


!guess values for water thickness in cm
REAL    :: h_itr(imax)
!linear system (tri-diagonal)
REAL    :: main(imax), super(imax), sub(imax), RHS(imax)
!solution 
REAL    :: h_temp (imax) !,  h_new  (imax)
!convergence test
REAL    :: relchg(imax)
!post processing
REAL    :: head  (imax), hygrad (imax)              ! head and hydraulic gradient
REAL    :: q_pav (imax), q_surf (imax), q_tot(imax) !fluxes      
!summary info
integer, parameter :: nmax = 2
INTEGER :: numit  (nmax), loc(nmax)  ! number of iterations at each timestep & locaiton of max change
REAL    :: maxdiff(nmax)             ! max change and its location
!FUNCTIONS
!       REAL    :: F_por, F_CC
!CHARACTERS
CHARACTER(8) DATE
CHARACTER(10) TIME



logical :: transition
real, dimension( imax, qmax) :: h_temp_hist, h_itr_hist !for monitoring the solution through the iteration
real :: relaxation_factor   !  Relaxation Factor
real, dimension( imax ) :: residual

real :: b  ! an extra value for b_pfc


real :: hs1, hs2, ds  ! Sheet flow MOC
real :: hp1, hp2, dx_moc ! PFC flow MOC
real :: cross_slope
integer, dimension(1) :: i_Zmax   !index of point with highest elevation
real :: eps_itr_tol

integer :: nip    ! number of interpolation points
!----------------------------------------
! Set pointers

X => XCV
Z => ZCV
DX => EDX

b = b_pfc
n = 1
i_Zmax = maxloc( Z ) 
!-----------------------------------------------------------------------------------
!SPATIAL GRID  (from GRID_1D_SECTION


! INITIALIZE ARRAYS
    !iteration array
    h_itr = h_old
    
    ! Linear system
    main  = 0.
    super = 0.
    sub   = 0.
    rhs   = 0.

    ! Summary arrays
    numit   = 0.
    maxdiff = 0.
    loc     = 0.

!-----------------------------------------------------------------------------------
! Solution using Crank-Nicolson with tri-diagonal matrix algorithm
! main, super & sub are diagonals of the coefficient matrix
! RHS is the right hand side of the linear system

open( unit = 50, file = 'PF_smry.csv',  status = 'REPLACE')
write(50,54) 'n/i', (/ (i,  i = 1, imax) /)

open( unit = 110, file = '1DRunDetails.txt', status = 'REPLACE')

54 format( ( A, ','), 10000( I, ',') )



!iteration loop
do q = 1, qmax 

            
! UPSTREAM BOUNDARY
i = 1        
! First cell
 PF  = F_por( h_old(i) )
 PF1 = F_por( h_itr(i) )

BC1: if( eta_1_BC .EQ. 'NO_FLOW') then
  ! NO FLOW BOUNDARY  Cw ---> 0
    Cw  = 0.  
    Ce  = F_CC( X(i  ), DX(i  ), h_old(i  ), Z(i),     &
                X(i+1), DX(i+1), h_old(i+1), Z(i+1) )
    !these coefficints are updated as the iteration progresses
    Cw1 = 0. 
    Ce1 = F_CC( X(i  ), DX(i  ), h_itr(i  ), Z(i),     &
                X(i+1), DX(i+1), h_itr(i+1), Z(i+1) )
    !Diagonals of coefficient matrix
    main (i) =  1  +  dt / 2 * PF1 * Cw1 / DX(i) &
                   +  dt / 2 * PF1 * Ce1 / DX(i) 
    super(i) =      - dt / 2 * PF1 * Ce1 / DX(i)
    sub  (i) =      - dt / 2 * PF1 * Cw1 / DX(i)    ! will be zero
    !Right hand side of linearized system
               ! part that came from n level            
    RHS  (i) =    h_old(i) + dt / 2.  * PF *                   &
             !  (  Cw / DX(i) * (  h_old(i-1) + Z(i-1)         &  !enforce BC
             !                   - h_old(i  ) - Z(i  ) )       &  !enforce BC
               ( + Ce / DX(i) * (  h_old(i+1) + z(i+1)         &
                                - h_old(i  ) - Z(i  ) )        &
                + rain  )                                      &
                ! part from n+1 level
                + dt / 2. * PF1 *                              &
            !   (   Cw1 / DX(i) * ( Z(i-1) - Z(i) )            &   !enforce BC
                ( + Ce1 / DX(i) * ( Z(i+1) - Z(i) )            &
                 + rain   )                                    

elseif( eta_1_BC .EQ. 'MOC_KIN') then
        ! METHOD OF CHARCTERISTICS KINEMATIC BOUNDAARY
        ! cross slope is defined positive for use in SQRT
        cross_slope = ( Z(i+1) - Z(i) ) / ( X(i+1) - X(i) )   
        ! If it comes out negative, a different BC is needed
        if( cross_slope .LT. 0.0) then
            write(110,*) 'Different BC needed at eta = 1'
        endif 
        if( h_old(i) .LE. b_pfc) then
                ! PFC FLOW MOC BC
                dx_moc =  K * (cross_slope) * dt / por   
                !Interpolate up the drainage slope...make sure we have enough points
                nip = NINT(dx_moc / DX(1) ) + 2
                write(110,*) 'eta_1_bc X=', (xcv(1) + dx_moc), &
                                    'nip=', nip,               &
                                     'KX=',   XCV(1:nip),      &
                                     'KY=', h_old(1:nip)
                hp1 = F_Linterp(       X = (xcv(1) + dx_moc), &                   
                                 Known_X = XCV( 1:nip)      , & 
                                 Known_Y = h_old(1:nip)     , &
                                       n = nip                  )
                if( hp1 .LT. 0.0) then
                        write(100,*) 'PFC1DIMP: eta_1_bc hp1=', hp1
                        stop
                endif
                hp2 = hp1 + rain * dt / por
                if( rain .LT. TINY ( rain )) then
                    ! Rainfall rate is effectively zero
                    hp2 = hp2
                else
                    ! Rainfall is non-zero, set a maximum value for hp2
                    ! the total drainage distance is the sum from the highest point in 
                    ! the 1D domain to the end, this is why MAXLOC is used.
                    hp2 = min( hp2, sum( DX( i : i_Zmax(1) ))*rain/K/cross_slope)
                endif
                ! Fill in linear system
                main(i) = 1.
                RHS (i) = hp2
! write(100,*) 'PFC1DIMP: eta_BC = MOC_KIN i=',i, 'dx_moc=', dx_moc, 'hp1=', hp1, 'hp2=', hp2
         else
                !SHEET FLOW MOC BC
                hs2 = h_old(i) - b_pfc
                ! Handle zero rainfall
                if( rain .LT. TINY( rain ) ) then
                    ! no increase in flow rate along drainge path
                    ! ds is arbitray, so use the PFC value
                    ds =  K * (cross_slope) * dt / por
                else
                    ds = sqrt( cross_slope ) / n_mann / rain *    &
                         ( ( hs2 + rain * dt )**(5./3.) - hs2**(5./3.) )
                endif
                ! Interpolate up the slope to find hs1
                nip = NINT(ds / DX(1) ) + 2
                write(110,*) 'eta_1_bc X=',(xcv(1) + ds), 'nip=', nip, 'KX=',   XCV(1:nip), 'KY=', h_old(1:nip)
                hs1 = F_Linterp(       X = (xcv(1) + ds),   &                   
                                 Known_X = XCV( 1:nip)  , &    !(/ 0.      ,     DX(i)  /) ,&
                                 Known_Y =  h_old(1:nip), &                            !(/ h_old(i), h_old(i+1) /), &
                                       n = nip ) - b_pfc
                ! Handle return to sheet flow
                if( hs1 .GT. 0.0 ) then
                    ! we have sheet flow
                    main(i) = 1.
                    RHS (i) =  b_pfc + ( hs1**(5./3.) + ( hs2 + rain*dt )**(5./3.) - hs2**(5./3.) )**0.6
                else
                    ! upstream point has sheet flow
                    ! use the PFC characteristic
                    main(i) = 1.
                    RHS (i) = hs1 + b_pfc + rain * dt / por 
                end if
           end if 
endif BC1

!Interior of domain
do i = 2, imax - 1
    PF  = F_por( h_old(i) )
    PF1 = F_por( h_itr(i) )
!        FUNCTION F_CC(    xin,  dxin,  hin, zin, &
!                         xout, dxout, hout, zout )
    !these coefficients are stationary (time level n)
    Cw  = F_CC( X(i  ), DX(i  ), h_old(i  ), Z(i),     &
                X(i-1), DX(i-1), h_old(i-1), Z(i-1) )
    Ce  = F_CC( X(i  ), DX(i  ), h_old(i  ), Z(i),     &
                X(i+1), DX(i+1), h_old(i+1), Z(i+1) )
    !these coefficints are updated as the iteration progresses
    Cw1 = F_CC( X(i  ), DX(i  ), h_itr(i  ), Z(i),     &
                X(i-1), DX(i-1), h_itr(i-1), Z(i-1) )
    Ce1 = F_CC( X(i  ), DX(i  ), h_itr(i  ), Z(i),     &
                X(i+1), DX(i+1), h_itr(i+1), Z(i+1) )


    !Diagonals of coefficient matrix
    main (i) =  1  +  dt / 2 * PF1 * Cw1 / DX(i) &
                   +  dt / 2 * PF1 * Ce1 / DX(i) 
    super(i) =      - dt / 2 * PF1 * Ce1 / DX(i)
    sub  (i) =      - dt / 2 * PF1 * Cw1 / DX(i)
    !Right hand side of linearized system
                ! part that came from n level
    RHS  (i) =    h_old(i) + dt / 2.  * PF *                   &
               (  Cw / DX(i) * (  h_old(i-1) + z(i-1)          &
                                - h_old(i  ) - Z(i  ) )        &
                + Ce / DX(i) * (  h_old(i+1) + z(i+1)          &
                                - h_old(i  ) - Z(i  ) )        &
                + rain  )                                      &
                ! part from n+1 level
                + dt / 2. * PF1 *                              &
               (   Cw1 / DX(i) * ( z(i-1) - z(i) )             &
                 + Ce1 / DX(i) * ( z(i+1) - Z(i) )             &
                 + rain   )
 end do

! DOWNSTREAM BOUNDARY  
i = imax 
! use BC from input argument
BC0:if( eta_0_BC .EQ. 'NO_FLOW') then
    !NO FLOW BOUNDARY ---> Ce == 0
    !these coefficients are stationary (time level n)
    Cw  = F_CC( X(i  ), DX(i  ), h_old(i  ), Z(i),     &
                X(i-1), DX(i-1), h_old(i-1), Z(i-1) )
    Ce  = 0.0 
    !these coefficints are updated as the iteration progresses
    Cw1 = F_CC( X(i  ), DX(i  ), h_itr(i  ), Z(i),     &
                X(i-1), DX(i-1), h_itr(i-1), Z(i-1) )
    Ce1 = 0.0
    !Diagonals of coefficient matrix
    main (i) =  1  +  dt / 2 * PF1 * Cw1 / DX(i) &
                   +  dt / 2 * PF1 * Ce1 / DX(i) 
    super(i) =      - dt / 2 * PF1 * Ce1 / DX(i)
    sub  (i) =      - dt / 2 * PF1 * Cw1 / DX(i)
    !Right hand side of linearized system
                ! part that came from n level
    RHS  (i) =    h_old(i) + dt / 2.  * PF *                   &
               (  Cw / DX(i) * (  h_old(i-1) + z(i-1)          &
                                - h_old(i  ) - Z(i  ) )        &
             !   + Ce / DX(i) * (  h_old(i+1) + z(i+1)          &    !enforce BC
            !                    - h_old(i  ) - Z(i  ) )        &    !enforce BC
                + rain  )                                      &
                ! part from n+1 level
                + dt / 2. * PF1 *                               &
               (   Cw1 / DX(i) * ( z(i-1) - z(i) )             &
            !     + Ce1 / DX(i) * ( z(i+1) - Z(i) )             &    !enforce BC
                 + rain   )
  elseif( eta_0_BC .EQ. 'MOC_KIN') then
        ! METHOD OF CHARCTERISTICS KINEMATIC BOUNDAARY
        cross_slope = ( Z(i-1) - Z(i) ) / ( X(i) - X(i-1) )    ! cross slope is positive downwars for use in SQRT
        if( cross_slope .LT. 0.0) then
            write(110,*) 'Different BC needed at eta = 0'
        endif
        if( h_old(i) .LE. b_pfc) then
                ! PFC FLOW MOC BC
                dx_moc =  K * (cross_slope) * dt / por       !NEED MAGNITUDE OF SLOPE AT CENTER OF EACH CELL               
                nip = NINT(dx_moc / DX(imax) ) + 2
                write(110,*) 'eta_0_BC X=', (XCV(imax) - dx_moc)  ,  &
                                    'nip=', nip,                     &
                                     'KX=', XCV( imax-nip+1:imax) ,  &
                                     'KY=', h_old( imax-nip+1 : imax)

                hp1 = F_Linterp(       X = (XCV(imax) - dx_moc)   ,  &                   
                                 Known_X = (XCV( imax-nip+1:imax)),  & 
                                 Known_Y = h_old( imax-nip+1 : imax) , & 
                                       n = nip ) 
                if( hp1 .LT. 0.0) then
                        write(100,*) 'PFC1DIMP: eta_0_bc hp1=', hp1
                        stop
                endif
                hp2 = hp1 + rain * dt / por
                if( rain .LT. TINY ( rain )) then
                    ! Rainfall rate is effectively zero
                    hp2 = hp2
                else
                    ! Rainfall is non-zero, set a maximum value for hp2
                    ! the total drainage distance is the sum from the highest point in 
                    ! the 1D domain to the end, this is why MAXLOC is used.
                    eta_0_hp2_max = sum( DX( i_Zmax(1) : i ))*rain/K/cross_slope
                    hp2 = min( hp2, eta_0_hp2_max)
                endif
                ! Fill in linear system
                main(i) = 1.
                RHS (i) = hp2
         else
                !SHEET FLOW MOC BC
                hs2 = h_old(i) - b_pfc
                ! Handle zero rainfall
                if( rain .LT. TINY( rain ) ) then
                    ! no increase in flow rate along drainge path
                    ! ds is arbitray, so use the PFC value
                    ds =  K * (cross_slope) * dt / por
                else
                    ds = sqrt( cross_slope ) / n_mann / rain *    &
                         ( ( hs2 + rain * dt )**(5./3.) - hs2**(5./3.) )
                endif
                ! Interpolate up the slope to find hs1
                nip = NINT(ds / DX(imax) ) + 2
                write(110,*) 'eta_0_BC X=', (XCV(imax) - ds)         , &
                                    'nip=', nip                      , &
                                     'KX=', XCV( imax-nip+1:imax)    , &
                                     'KY=', h_old( imax-nip+1 : imax)

                hs1 = F_Linterp(       X = (XCV(imax) - ds)          , &                   
                                 Known_X = (XCV( imax-nip+1:imax))   , &
                                 Known_Y = h_old( imax-nip+1 : imax) , &
                                       n = nip ) - b_pfc
                ! Handle return to sheet flow
                if( hs1 .GT. 0.0 ) then
                    ! we have sheet floe
                    main(i) = 1.
                    RHS (i) =  b_pfc + ( hs1**(5./3.) + ( hs2 + rain*dt )**(5./3.) - hs2**(5./3.) )**0.6
                else
                    ! upstream point has sheet flow
                    ! use the PFC characteristic
                    main(i) = 1.
                    RHS (i) = hs1 + b_pfc + rain * dt / por 
                end if
           end if 
  end if BC0  



! TRANSITION CHECK
!   test to see if there is a transition to or from sheet flow
!   happening during this timestep.   Use under-relaxtion to
!   control oscillations during a transition timestep.

transition = .false.
do i = 1, imax
    pf = F_por( h_old(i) )
    pf1= F_por( h_itr(i) )
    if( pf .GT. pf1 .OR. pf .LT. pf1) then
        transition = .true.
        write(110,*) 'PERFCODE: transition for cell i=', i, 'pf=', pf, 'pf1=', pf1 
    endif
end do

if( transition .eqv. .true. ) then
    relaxation_factor = relax_tran
    eps_itr_tol = tolit * 10.
else
    relaxation_factor = relax
    eps_itr_tol = tolit
endif


        !Solve linear system
        CALL THOMAS(main, super, sub, RHS, h_temp, imax)



! Compute residual and relative change for this iteration.  This took
!  some careful though to handle both filling and draining cases.
!  relative change is used when the solution is far from zero
!  and absolute change (residual) is used near zero.

do i = 1, imax
 
    if( h_temp(i) .GT. TINY( h_temp(i) ) ) then

            ! Compute residual for this iteration
            residual(i) =  h_temp(i) - h_itr(i)

            ! Handle a result that is effectively zero by
            ! using an absolute tolerance instead of 
            ! a relative one
            if( h_temp(i) .LE. h_pfc_min  .and. &
                residual (i) .LE. eps_itr_tol     ) then
        
                    relchg(i) = 0.0
            
            else
                    relchg (i) =  residual (i) / h_itr(i)
            endif
            
    elseif( h_temp(i) .LE. TINY( h_temp(i) ) ) then

            ! the model is saying the cell is empty, 
            ! so force the solution to be zero      
            h_temp(i) =  0.0
            ! compute the residual
            residual(i) =  h_temp(i) - h_itr(i)
            ! For the zero case, use an absolute rather than
            ! relative tolerance by setting the value of relchng
            ! below the tolerance instead of computing it.
            if( abs( residual(i) ) .LE. eps_itr_tol ) then
                    
                    relchg(i) = 0.0
            endif
    endif

end do
         

if( maxval( h_temp) .LT. TINY( h_temp(1) )  ) then
       write(*,*) 'PFC1DIMP: Zeroed out.  Writing system and stopping program'
       open( unit = 10, file = '1Dsystem.csv', status = 'REPLACE' )
       write(10,*) 'i,sub,main,super,rhs,h_temp,'
       do i = 1, imax
             write(10, 10) i, sub(i), main(i), super(i), RHS(i), h_temp(i)
       end do
       close( 10 )

       call write_vector( h_old, imax, 'h_old_1d.csv')
       STOP
10  FORMAT (  (I, ','), 5(E, ',') ) 

end if


!perform usual iteration check
IF ( maxval ( ABS( relchg ) ) .le. eps_itr_tol .AND.  & 
      F_L2_NORM( relchg, imax )   .le. eps_itr_tol       ) then
        !  WRITE(*,*) 'Time step n = ', n,' converged in q = ', q, ' iterations.'
          EXIT
      end if    




! Smith page 32

h_itr = h_itr + relaxation_factor * residual

!Store the result of this iteration 
h_temp_hist( : , q ) = h_temp
h_itr_hist ( : , q ) = h_itr



 WRITE(110,*) 'ITERATION q=', q                     , &
             'Max Change of', maxval( abs( relchg) ), &
                     'at i=', maxloc( abs( relchg)) , &
                'h_temp(i)=', h_temp( maxloc( abs( relchg)) ) , &
                  'L2_Norm=', F_L2_NORM( relchg, imax) 



!end iteration loop
end do

!update the old and new solutions
!At the end of the iteration, we have found values for the
!next time step.

h_new  = h_temp


!Store summary info for this timestep
numit  ( n ) = q   
loc    ( n ) = maxloc ( abs( relchg ), dim=1 )     
maxdiff( n ) = relchg ( loc ( n ) ) 

!Give Error if Iteration fails to converge
if (q .gt. qmax) then
  WRITE(*,*)   'PFC1DIMP: Iteration failed to converge. '
  write(100,*)  'PFC1DIMP: Iteration failed to converge. '
  CALL WRITE_MATRIX( h_temp_hist, imax, qmax, 'h_temp_hist_1D.csv')
  CALL WRITE_MATRIX( h_itr_hist , imax, qmax, 'h_itr_hist_1D.csv')

!          EXIT
end if


close(50)   ! pf summary file
close(110)  ! Run details file
!-----------------------------------------------------------------------------------
! POST PROCESSING

!Compute head
head(:) = h_new(:) + Z(:)

!Compute hydraulic gradient and flux ( both positive downwards)

!Upstream boundary node (using a 1-sided approximation)
i = 1
hygrad(i) =  ( head(i) - head(i+1) )  /  ( X(i+1) - X(i) )
!In the pavement
hp = min( h_new(i), b )
q_pav (i) = K * hp * hygrad(i)
!on the surface
hs = max( 0., h_new(i) - b )
q_surf(i) = 1. / n_mann * hs ** (2./3.) * sqrt( abs( hygrad(i) ) ) * hs
q_tot(i) = q_pav(i) + q_surf(i)

do i = 2, imax - 1
    !hydrualic gradient
    hygrad(i) =  ( head(i-1) - head(i+1) )  /  ( X(i+1) - X(i-1) )
    !thickness in the pavement
    hp = min( h_new(i), b )
    q_pav (i) = K * hp * hygrad(i)
    !thickness on the surface
    hs = max( 0., h_new(i) - b )
    q_surf(i) = 1. / n_mann * hs**(2./3.) * sqrt( abs( hygrad(i) ) ) * hs
!             WRITE(*,*) 'i = ', i, 'hs = ', hs, 'q_surf =', q_surf(i)
    q_tot(i) = q_pav(i) + q_surf(i)
end do

! DOWNSTREAM BOUNDARY ( 1 sided approximation)
i = imax
!hydrualic gradient
hygrad(i) =  ( head(i-1) - head(i) )  /  ( X(i) - X(i-1) )
!thickness in the pavement
hp = min( h_new(i), b )
q_pav (i) = K * hp * hygrad(i)
!thickness on the surface
hs = max( 0., h_new(i) - b )
q_surf(i) = 1. / n_mann * hs ** (2./3.) * sqrt( abs( hygrad(i) ) ) * hs
q_tot(i) = q_pav(i) + q_surf(i)

!-----------------------------------------------------------------------------------
!Write results to a file

        call DATE_AND_TIME(DATE,TIME)
        call CPU_TIME(cputime)

        OPEN(UNIT = 10, FILE = 'pfc1Dimp.csv', STATUS='REPLACE')
        WRITE(10,*) 'Output From pfc1Dimp.f95'
        WRITE(10,*) 'Timestamp,', DATE,' ', TIME,','
        WRITE(10,*) 'Upstream Boundary == Fixed Value'
        WRITE(10,*) 'Downstream Boundary == Sf = So'
        WRITE(10,200) 'Hydraulic Conductivity (cm/s),', k
        WRITE(10,200) 'Rainfall Intensity (cm/hr),', rain * 3600.
        WRITE(10,200) 'PFC Thickness (cm),', b 
        WRITE(10,200) 'Final Time (sec),', ( n-1 ) * dt 
        WRITE(10,200) 'Time step (seconds),', dt 
        WRITE(10,200) 'Grid spacing (cm),', dx(5) 
        WRITE(10,*) 'Number of elements,', imax - 2
        WRITE(10,200) 'CPU Time (seconds),', cputime
        WRITE(10,*) '************************** &
                     &MODEL OUTPUT IN [ SI  ] UNITS &
                     &***********************************'
        WRITE(10,*) 'X, eta,    Z, PFC Surface, Thickness, Head,', &
                    'Hydraulic Gradient, Pavement Flux,'         , &
                    ' Surface Flux, Total Flux,'
        do i = 1, imax
                WRITE(10,100) X(i), etaCV(i), Z(i), Z(i) + b,  &
                              h_new(i), Head(i), hygrad(i),    &
                              q_pav(i), q_surf(i), q_tot(i)    
        END do
        CLOSE(10)

!-----------------------------------------------------------------------------------
!Write calculation summary to file

    OPEN( UNIT = 20, FILE = '1Ddetails.csv', STATUS='REPLACE')
    WRITE(20,*) 'Timestamp,', DATE, ' ', TIME, ','
    WRITE(20,*) '-----,'
    WRITE(20,*) 'Timestep, Iterations, MaxRelChng, MaxLocn'
    DO n = 1, nmax
        WRITE(20,300) n, numit(n), maxdiff(n), loc(n)
    end do
    close(20)    
        
!-----------------------------------------------------------------------------------
!Format statements

100     FORMAT (100 ( F14.7, ',' ) )        !  Formatting for the actual output
200     FORMAT ( A, F10.4 )
300     FORMAT ( 2 (I, ','), E, ',', I, ',' )

!-----------------------------------------------------------------------------------
        END subroutine pfc1Dimp
!===================================================================================
!       \\\\\\\\\\      E N D      S U B R O U T I N E  //////////
!       //////////          P F C 1 D I M P             \\\\\\\\\\
!===================================================================================







END MODULE pfc1Dsubs
