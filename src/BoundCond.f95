! fortran_free_source
!
!  This module is part of PERFCODE, written by Bradley J Eck
!
MODULE BoundCond 
implicit none
contains
!
!   1. Subroutine MOC_KIN_BC
!   2. Subroutine CHECK_BCs
!   3. Function   F_RHS_n_CURB_IN
!   4. Function   F_RHS_n1_CURB_IN
!   5. Function   F_Qout
!
!===================================================================
!   \\\\\\\\\\       B E G I N    S U B R O U T I N E ///////////
!   //////////          M O C _ K I N _ B C           \\\\\\\\\\\
!==================================================================
subroutine MOC_KIN_BC( i, j, rain, dt, side,  h_bound, dev)
!   Uses the method of characteristics (MOC) to solve the 1D 
!       equation under the assumption of kinematic flow and applies
!       this solution as a boundary condition.
!
!  History:   March 2010  -- Original coding
!             5 Jul 2011  -- changed variable 'side' to integer
!             10 Jul 2011 -- added code for western boundary
!             12 Jul 2011 -- moved drainage slope calc to its own
!                            subroutine in pfc2Dsubs so that it 
!                            could be used by other subroutines
!                            like the CURB_IN boundary condition
use shared   , only: hyd_cond, por, b_pfc, n_mann, CV_Info, wid, &
                     imax, jmax, lng, wid, h_old, Z, eta_0_hp2_max, &
                     north, south, east, west
use pfc2dsubs, only: set_xyh, DRAIN_SLOPE
use pfc2Dfuns, only: F_LinearIndex
use utilities, only: BILINEAR_INTERP, F_PythagSum 
implicit none
!---------------------------------------------------------------------------
! Arguments
integer, intent( in ) :: i, j
real, intent ( in )   :: dt      ! timestep
integer, intent(in) :: side ! which side of the domain are we working on
real, intent( in ) :: rain       ! rainfall rate for this timestep 
real, intent( out )   :: h_bound
integer, optional :: dev  !device for outputing errors 
! Internal variables
real, dimension( 2 ) :: S_drain
real, dimension( 2 ) :: S_drain_unit  ! slope vector for drainage slope
real :: drainage_slope  ! magnitude of drainge slope

integer :: v
integer :: vi1  !value of v for the cell i+1
integer :: vj1  !value of v for the cell j+1
integer :: vjm1 !value of v for the cell j-1
integer :: vim1 !value of v for the cell i-1


! Bilinear Interpolation
real :: XX, YY       ! Coordinates of point where depth is interpolated
real :: x1, y1, h1   ! Coordinates of point 1  Interpolation points
real :: x2, y2, h2   !   "    "       point 2
real :: x3, y3, h3   !   "    "       point 3
real :: x4, y4, h4   !   "    "       point 4

! Method of Characteristics
!   PFC
real :: dx_moc, hp1, hp2
!   Sheet flow
real :: ds, hs1, hs2

integer :: device
logical :: bilin_err

!---------------------------------------------------------------------

! Default values for output device
if( present( dev ) .EQV. .FALSE. ) then
        device = 6
else
        device = dev
end if

!--------------------------------------------------
! Compute vector for the drainage slope ( S_drain )
! and its magnitude, and unit vector for direction
!-------------------------------------------------

! compute drainage slope vector
call DRAIN_SLOPE( i, j, side, S_drain )

! ! and the magnitude
drainage_slope = F_PythagSum( S_drain(1), S_drain(2) )
! and a drainage slope unit vector
S_drain_unit = S_drain / drainage_slope

!write( device, *) 'Slope Vectors:  S_ksi = ', S_ksi, ' S_eta', S_eta


!-----------------------------------------------------------------------
!    I N T E R P O L A T I O N     P O I N T S
!-----------------------------------------------------------------------


! now we can figure out which points to use for the
! bilinear interploation routine. Points must be specified 
! counter-clockwise around the perimeter:
!
!       4-----3
!       |     |
!       |     |
!       1-----2

if    ( (side == south)  .AND. (S_drain_unit(1) .LE. 0.) )  then
        ! This is the southern boundary and 
        ! The domain slopes from left to right
        !Point 1
        call set_xyh( i-1, j  , x1, y1, h1 )
        !Point 2
        call set_xyh( i  , j  , x2, y2, h2 )
        !Point 3
        call set_xyh( i  , j+1, x3, y3, h3 )
        !Point 4
        call set_xyh( i-1, j+1, x4, y4, h4)

elseif( (side == south)  .AND.  (S_drain_unit(1) .GE. 0.0) ) then
        ! This is the southern boundary and 
        ! The domain slopes from right to left
        ! Point 1
        call set_xyh( i  , j  , x1, y1, h1 )
        ! Point 2
        call set_xyh( i+1, j  , x2, y2, h2 )
        !Point 3
        call set_xyh( i+1, j+1, x3, y3, h3 )
        !Point 4
        call set_xyh( i  , j+1, x4, y4, h4 )

elseif( (side == north) .AND. (S_drain_unit(1) .LE. 0.0) ) then
        ! This is the northern boundary and 
        ! The domain slopes from left to right
        call set_xyh( i-1, j-1, x1, y1, h1 )
        call set_xyh( i  , j-1, x2, y2, h2 )
        call set_xyh( i  , j  , x3, y3, h3 )
        call set_xyh( i-1, j  , x4, y4, h4 )

elseif( (side == north) .AND. (S_drain_unit(1) .GE. 0.0) ) then
        ! This is the northen boundary and 
        ! The domain slopes from right to left
        call set_xyh( i  , j-1, x1, y1, h1 )
        call set_xyh( i+1, j-1, x2, y2, h2 )
        call set_xyh( i+1, j  , x3, y3, h3 )
        call set_xyh( i  , j  , x4, h4, h4 )


elseif( (side == east)  .AND. (S_drain_unit(2) .GE. 0.0) ) then
        ! This is the eastern boundary and
        ! and uphill is the positive Y direction
        call set_xyh( i-1, j  , x1, y1, h1 )
        call set_xyh( i  , j  , x2, y2, h2 )
        call set_xyh( i  , j+1, x3, y3, h3 )
        call set_xyh( i-1, j+1, x4, y4, h4 )


elseif( (side == east) .AND. (S_drain_unit(2) .LT. 0.0) ) then
        ! This is the eastern boundary and
        ! and uphill is the negative Y direction
        call set_xyh( i-1, j-1, x1, y1, h1 )
        call set_xyh( i  , j-1, x2, y2, h2 )
        call set_xyh( i  , j  , x3, y3, h3 )
        call set_xyh( i-1, j  , x4, y4, h4 )


elseif( ( side == west) .and. (S_drain_unit(2) .gt. 0.0) ) then
        ! This is the western boundary and
        ! and uphill is the positive y direction
        call set_xyh( i  , j  , x1, y1, h1 )
        call set_xyh( i+1, j  , x2, y2, h2 )
        call set_xyh( i+1, j+1, x3, y3, h3 )
        call set_xyh( i  , j+1, x4, y4, h4 )

elseif( side == west .and. S_drain_unit(2)  .lt. 0.0 ) then
        ! This is the western boundary and
        ! and uphill is the negative Y  direction
        call set_xyh( i  , j-1  , x1, y1, h1 )
        call set_xyh( i+1, j-1  , x2, y2, h2 )
        call set_xyh( i+1, j    , x3, y3, h3 )
        call set_xyh( i  , j    , x4, y4, h4 )

endif




!----------------------------------------------------------------------
!   M E T H O D     O F     C H A R A C T E R I S T I C S 
!----------------------------------------------------------------------


! Reset v to confirm we're in the right cell
v = F_LinearIndex( i, j, jmax )


MOC:if( h_old(i,j) .LE. b_pfc(i,j) ) then
!------------------
!PFC FLOW
!-----------------
        ! Sheet flow has not started yet
        ! use MOC to estimate the solution at the next time step
        ! figure out how far up the drainage slope to go
        dx_moc = hyd_cond(i,j) * (drainage_slope) * dt / por(i,j)
!        write( device, * ) 'MOC_KIN_BC: i = ', i, ' j = ', j, 'pfc char len = ', dx_moc
        ! and the coordinates of this location
        XX = CV_Info( v ) % X  +  dx_moc * S_drain_unit( 1 )
        YY = CV_Info( v ) % Y  +  dx_moc * S_drain_unit( 2 )
        ! use bilinear interpolation to find the
        ! thickness (hp1) at this location
        call BILINEAR_INTERP(  XX, YY, hp1,  &
                               x1, y1, h1 ,  &
                               x2, y2, h2 ,  &
                               x3, y3, h3 ,  &
                               x4, y4, h4 ,  &
                               device, bilin_err     )  

        ! value at next time step
        hp2 = hp1 + rain * dt / por(i,j)
        ! set maximum value for hp2 (1D flow)      
        if( rain .LT. TINY ( rain )) then
                ! Rainfall rate is effectively zero
                hp2 = hp2   ! Eqv to hp2 = hp1 
        else
                ! Rainfall is non-zero, set a maximum value for hp2
                hp2 = min( hp2, sum(wid(i,:))*rain/hyd_cond(i,j)/drainage_slope)
                ! Use hp1 (basically zero rainfall) if there
                !  is a decrease in depth
                if( hp2 .LT. hp1 ) then
                        hp2 = hp1
                end if
        endif

       
!
!     ! Error checking for eastern boundary
!     if( i == imax ) then
!         if( j == jmax -5 .or. j == jmax/2 .or. j == 5 ) then
!             
!             write( device, *)   'MOC_KIN: i =', i ,          &
!                                          'j =', j ,          &
!                                    'S_drain =', S_drain,     & 
!                               ' drainage_slope =', drainage_slope, &
!                              ' S_drain_unit =', S_drain_unit
!             write(device,*) 'Bilinear Interpolation'
!             write(device,*) '       X,         Y,              h,' 
!             write(device,32) 0, XX, YY, hp1
!             write(device,32) 1, x1, y1,  h1
!             write(device,32) 2, x2, y2,  h2     
!             write(device,32) 3, x3, y3,  h3
!             write(device,32) 4, x4, y4,  h4  
!         end if
!     endif
!

        ! error checking for interpolation
        if( bilin_err .eqv. .true. ) then
            write(device,*) 'MOC_KIN_BC: Bilinear interpolation error &
                            & for grid cell i = ', i, ' j = ', j
            write(device,*)          'S_drain=', S_drain,     & 
                                ' drainage_slope=', drainage_slope, &
                               ' S_drain_unit=', S_drain_unit   
            write(device,*)              'hp2=', hp2   , &
                                      'dx_moc=', dx_moc, &
                                         'hp1=', hp1   , &
                                        'rain=', rain  , &
                                          'dt=', dt    , &
                                         'por=', por          

            write(device,* ) 'Interploation points/result:'
            write(device,* ) '       X,         Y,              h,' 
            write(device,32) 0, XX, YY, hp1
            write(device,32) 1, x1, y1,  h1
            write(device,32) 2, x2, y2,  h2     
            write(device,32) 3, x3, y3,  h3
            write(device,32) 4, x4, y4,  h4  

        end if





! ERROR CHECKING CODE for MOC BC
!        if( i == imax/2  .OR. j == jmax/2  ) then
!           write(device,*) 'PFC Flow MOC BC: i=', i     , &
!                                            'j=', j     , &
!                                          'hp2=', hp2   , &
!                                       'dx_moc=', dx_moc, &
!                                          'hp1=', hp1   , &
!                                         'rain=', rain  , &
!                                           'dt=', dt    , &
!                                          'por=', por
!        endif


        h_bound = hp2
else
!-----------------
!SHEET FLOW
!-----------------
        hs2 =  h_old(i,j) - b_pfc(i,j)
        ! Handle Zero Rainfall        
        if( rain .LT. TINY( rain ) ) then
            ! there is no increase in flow rate along the drainage path
            !  and ds becomes arbitrary so use the characteristic length for PFC flow
            ! b/c you might need it later
            ds = hyd_cond(i,j) * ( drainage_slope ) * dt / por(i,j)
        else
            ds = sqrt( drainage_slope ) / n_mann(i,j) / rain *      &
                 ( ( hs2 + rain*dt )**(5./3.) - hs2**(5./3.) )
        end if
       ! interpolate up the drainage slope to find hs1
        XX = CV_Info( v ) % X  +  ds * S_drain_unit( 1 )
        YY = CV_Info( v ) % Y  +  ds * S_drain_unit( 2 )
        ! use bilinear interpolation to find the thickness (hs1) at this location
        call  BILINEAR_INTERP(  XX, YY, hs1,  &
                                x1, y1, h1 ,  &
                                x2, y2, h2 ,  &
                                x3, y3, h3 ,  &
                                x4, y4, h4 ,  &
                                device, bilin_err ) 
        if( bilin_err .eqv. .true. ) then
            write(device,*) 'MOC_KIN_BC: Bilinear interpolation error &
                            & for grid cell i = ', i, ' j = ', j

              write(device,*) '       X,         Y,              h,' 
              write(device,32) 0, XX, YY, hs1
              write(device,32) 1, x1, y1,  h1
              write(device,32) 2, x2, y2,  h2     
              write(device,32) 3, x3, y3,  h3
              write(device,32) 4, x4, y4,  h4

        end if
        
        ! subtract off the pavement thickness
        hs1 = hs1 - b_pfc(i,j)
        !Handle return to PFC flow
        if( hs1 .GT. 0. ) then
            !we have sheet flow
            !Output some summary info
!            if( i == imax/2 .or. j == jmax / 2 ) then
!               write(device,*) 'Sheet Flow MOC BC: i=',   i, &
!                                                  'j=',   j, &
!                                                'hs2=', hs2, &
!                                                 'ds=',  ds, &
!                                                'hs1=',  hs1
!            endif
!            !checking for good values of inputs
!            if( hs1 .LT. 0.  .OR.  hs2 .LT. 0.  .OR. rain .LT. 0.) then
!                write(device,*) 'Sheet Flow MOC BC: i=',i, 'j=',j, &
!                                 'hs1=', hs1, 'hs2=',hs2, 'rain=',rain
!            end if
            ! return value for the boundary
            h_bound =  b_pfc(i,j) + (   hs1**(5./3.) +             & 
                                      ( hs2 + rain*dt )**(5./3.) - &
                                        hs2**(5./3.) )**0.6
        else
            ! the upstream point does not have sheet flow
            ! use PFC characterisitic
            h_bound = hs1 + b_pfc(i,j) + rain * dt / por(i,j)  
        end if 
                   
end if MOC



!------------------------------------------------------------------
! Format statements
31  format( 3( F12.7, '  ') )
32  format( I3, '  ', 3(F12.7, '  ') )

!-------------------------------------------------------------------
end subroutine MOC_KIN_BC
!===================================================================
!   \\\\\\\\\\    E N D       S U B R O U T I N E     ///////////
!   //////////          M O C _ K I N _ B C           \\\\\\\\\\\
!==================================================================


!===================================================================
!   \\\\\\\\\\  B E G I N     S U B R O U T I N E     ///////////
!   //////////          C H E C K _ B C s             \\\\\\\\\\\
!==================================================================
subroutine CHECK_BCs( model_mode )
!   Confirms that the BCs given as inputs are acceptable/have been
!       implemented in the model.  This subroutine should be modified
!       along with APPLY_BCs in module pfc2Dsubs.
use shared, only: north_bc_input, south_bc_input, east_bc_input, &
                  west_bc_input, north_bc, south_bc, east_bc,    &
                  west_bc, imax, jmax

implicit none
! Arguments
character(len=2), intent( in) :: model_mode
! Internal
integer :: i

!-----------------------------------------------------------------------


! 1D_FLOW not allowed for north or south boundary
if( model_mode .eq. '1D' ) then
    if( north_bc_input .eq. '1D_FLOW' .or. &
        south_bc_input .eq. '1D_FLOW'      ) then
        
        write(*,*) 'CHECK_BCs:  1D_FLOW not allowed for North or South boundary.'
        STOP
   
    end if


    if( north_bc_input .eq. 'CURB_IN' .or. &
        south_bc_input .eq. 'CURB_IN'      ) then
        
        write(*,*) 'CHECK_BCs:  CURB_IN not implemented for 1D model.'
        STOP
   
    end if


elseif( model_mode .eq. '2D' ) then

    do i = 1, imax
        if( north_bc(i) .eq. '1D_FLOW' .or. &
            south_bc(i) .eq. '1D_FLOW'      ) then

            write(*,*) 'CHECK_BCs: 1D_FLOW not allowed for North or South boundary.'
            STOP
       
        end if

    end do    
        
end if



write(*,*) '' 
write(*,*) 'CHECK_BCs:  Boundary conditions OK.'
write(*,*) '' 

end subroutine    
!===================================================================
!   \\\\\\\\\\    E N D       S U B R O U T I N E     ///////////
!   //////////          C H E C K _ B C s             \\\\\\\\\\\
!==================================================================


!===================================================================
!   \\\\\\\\\\    B E G I N       F U N C T I O N      ///////////
!   //////////   F _ R H S _ n _ C U R B _ I N         \\\\\\\\\\\
!==================================================================
Function F_RHS_n_CURB_IN(i, j, Cw, Ce, Cs, Cn, rr, pf, dt, side ) &
Result(Fn)
!   Computes the right hand side of the linear system at time level
!       n when the boundary condtion is CURB_IN.  See also the 
!       similar function in the module pfc2Dfuns.
!
!   HISTORY: July 10, 2011 -- Original coding
!
use shared, only: h_old, Z, g, lng_south, south, b_pfc, hyd_cond
use pfc2Dsubs, only: DRAIN_SLOPE
implicit none
! Arguments
integer, intent(in) :: i, j, side
real   , intent(in) :: Cw, Ce, Cs, Cn, rr, pf, dt
! Internal variables
real :: hw, he, hn, hs, Zw, Ze, Zs, Zn !values at center of adjacent cells
real :: h_curb ! sheet flow depth at curb inlet (hs)
real :: V_curb ! sheet flow velocity into curb inlet
real :: Qs_curb  ! sheet flow into curb inlet
real :: Qp_curb  ! pfc flow into curb inlet
real :: Fn
real, dimension(2) :: S_drain  ! drainage slope vector 
!----------------------------------------------------------------
! Depth & elevs for east and west sides 
! Depth
hw = h_old(i-1,j)
he = h_old(i+1,j)
! Elev
Zw = Z(i-1,j)
Ze = Z(i+1,j)


if( side .eq. south ) then

    ! Depths
    hn = h_old(i,j+1)
    ! Elevations
    Zn = Z(i,j+1)
    ! curb
    h_curb = h_old(i,j) - b_pfc(i,j) + ( Zn - Z(i,j) ) / 2.
    V_curb = sqrt( 2. * g * h_curb / 3. )
    Qs_curb = v_curb * (2. * h_curb / 3. ) * lng_south(i,j)

    ! Flow out of the pfc
    !   Drainage slope at the cell center
    call DRAIN_SLOPE(i, j, side, S_drain)
    !   Darcy's law to compute flow out the pfc at the face
    if( S_drain(2) .LT. 0.0 ) then
        write(*,*) 'F_RHS_n1_CURB_IN:  S_drain(2) is negative...Stopping program'
    end if
    Qp_curb = hyd_cond(i,j) * S_drain(2) * b_pfc(i,j) * lng_south(i,j)

    ! Fn
    Fn = h_old(i,j)   +                        &
         pf * dt / 2. * (  Cw * hw             &
                         + Cn * hn  +  Ce * he &
                         + Cw * Zw             &
                         + Cn * Zn  +  Ce * Ze &
                         -  (Cw      + Cn + Ce) * h_old(i,j)       &  
                         -  (Cw      + Cn + Ce) *     Z(i,j)       &
                         + rr                                      &  
                         - Qs_curb - Qp_curb )

end if

end function F_rhs_n_curb_in
!===================================================================
!   \\\\\\\\\\    E N D           F U N C T I O N      ///////////
!   //////////   F _ R H S _ n _ C U R B _ I N         \\\\\\\\\\\
!==================================================================



!===================================================================
!   \\\\\\\\\\    B E G I N       F U N C T I O N        ///////////
!   //////////   F _ R H S _ n 1 _ C U R B _ I N         \\\\\\\\\\\
!==================================================================
Function F_RHS_n1_CURB_IN(i, j, Cw1, Ce1, Cs1, Cn1, rr, pf, dt, side ) &
Result(F1)
!   Computes the right hand side of the linear system at time level
!       n+1 when the boundary condtion is CURB_IN.  See also the 
!       similar function in the module pfc2Dfuns.
!
!   HISTORY: July 10, 2011 -- Original coding
use shared, only: Z, g, lng_south, south, h_itr, b_pfc, hyd_cond
use pfc2Dsubs, only: DRAIN_SLOPE
implicit none
! Arguments
integer, intent(in) :: i, j, side
real   , intent(in) :: Cw1, Ce1, Cs1, Cn1, rr, pf, dt

! Internal variables
real :: Zw, Ze, Zs, Zn !values at center of adjacent cells
real :: h_curb ! depth at curb inlet (hs)
real :: V_curb ! velocity into curb inlet
real :: Qp_curb
real :: F1
real, dimension(2) :: S_drain
!-----------------------------------------------------------------------

Zw = Z(i-1,j)
Ze = Z(i+1,j)


if( side .eq. south ) then

    Zn = Z(i,j+1)

    ! curb
    h_curb = h_itr(i,j) - b_pfc(i,j) + ( Zn - Z(i,j) ) / 2.
    V_curb = sqrt( 2. * g * h_curb / 3. )

    !   Drainage slope at the cell center
    call DRAIN_SLOPE(i, j, side, S_drain)
    if( S_drain(2) .LT. 0.0 ) then
        write(*,*) 'F_RHS_n1_CURB_IN:  S_drain(2) is negative...Stopping program'
    end if
    !   Darcy's law to compute flow out the pfc at the face
    Qp_curb = hyd_cond(i,j) * S_drain(2) * b_pfc(i,j) * lng_south(i,j)
    

    F1 = pf * dt / 2. * (  Cw1 * Zw   &
                         + Cn1 * Zn  +  Ce1 * Ze &
                         - (Cw1 + Cn1 + Ce1)* Z(i,j)   &  
                         + rr &
                         - V_curb * lng_south(i,j) * 2./3. * (Zn - Z(i,j)) / 2. &
                         - Qp_curb )  

end if

end function F_rhs_n1_curb_in
!===================================================================
!   \\\\\\\\\\    E N D           F U N C T I O N      ///////////
!   //////////   F _ R H S _ n 1 _ C U R B _ I N         \\\\\\\\\\\
!==================================================================

!====================================================================
!   \\\\\\\\\\    B E G I N         F U N C T I O N      ///////////
!   //////////          F _ Q O U T _ B O U N D          \\\\\\\\\\\
!====================================================================
Function F_Qout_bound( i, j) result(Qout)
!   Computes the flow out of a boundary cell by difference
!      once the solution for a timestep has been found.
!      This is a post-processed quanitity at each timestep
!      and is used for the global water balance.
!
!       |----Q2---|      The boundary face is Qout
!       |         |      other faces numbered around
!      Q1    *    Q3     clockwise.
!       |         |
!       |---Qout--|
!   
use shared, only:  b_pfc, por, area, imax, jmax,       &
                   is_boundary, is_corner,             &
                   h_old, h_new, rain,                 &
                   dt, n, north, south, east, west
use pfc2Dsubs, only: CALC_Qx
implicit none
! Arguments
integer, intent(in) :: i, j
! Internal variables
real :: hp, hs
real :: Vold, Vnew, Qout, Qrain, Q1, Q2, Q3

!-------------------------------------------------------------
! error checking
if( is_boundary(i,j) .eqv. .false.) then
    write(*,*) 'F_Qout: Function called for non-boundary cell.'
    write(*,*) '        STOPPING PROGRAM'
    STOP
endif

! Compute cell volumes at old and new time levels
!   Vnew
hp = min( h_new(i, j),  b_pfc(i, j) )
hs = max( h_new(i, j) - b_pfc(i, j), 0.0 ) 
Vnew =  hp * area(i,j) * por(i, j)  +  hs * area(i,j) 
!  Vold
hp = min( h_old(i, j),  b_pfc(i, j) )
hs = max( h_old(i, j) - b_pfc(i, j), 0.0 ) 
Vold =  hp * area(i,j) * por(i, j)  +  hs * area(i,j) 

! Rainfall flow rate
Qrain = (rain(n) + rain(n-1))/2. * area(i,j)

! Compute flow rates across internal faces and
!  solve water balance to get Qout 
!-------------------------------------------
if( is_corner(i,j) .eqv. .false. ) then
!-----------------------------------------
  ! We compute the flow across three faces 
  ! and use these to find flow out the fourth face.
  if( j == 1 ) then

      !We are computing the flow out the SOUTHERN face
      ! Q1 = west
      ! Q2 = north
      ! Q3 = east
      call calc_Qx( i, j, west,  Q1 )
      call calc_Qx( i, j, north, Q2 )
      call calc_Qx( i, j, east,  Q3 )

  elseif( i == 1 ) then
      
      !We are computing the flow out the WESTERN face
      ! Q1 = north
      ! Q2 = east
      ! Q3 = south
      call calc_Qx( i, j, north, Q1 )
      call calc_Qx( i, j, east,  Q2 )
      call calc_Qx( i, j, south, Q3 )

  elseif( j == jmax ) then

      !We are computing the flow out the NORTHERN face
      ! Q1 = east
      ! Q2 = south
      ! Q3 = west
      call calc_Qx( i, j, east,  Q1 )
      call calc_Qx( i, j, south, Q2 )
      call calc_Qx( i, j, west,  Q3 )

  elseif( i == imax ) then

      !We are computing the flow out the EASTERN face
      ! Q1 = south
      ! Q2 = west
      ! Q3 = north
      call calc_Qx( i, j, south, Q1 )
      call calc_Qx( i, j, west,  Q2 )
      call calc_Qx( i, j, north, Q3 )

  endif
!---------------------------------------
elseif( is_corner(i,j) .eqv. .true. ) then
!---------------------------------------
  ! We compute flow across two interior faces and
  ! use these to find the outflow across the boundary
  Q3 = 0.0

  if( i == 1 .and. j == 1 ) then

      ! corner = SouthWest
      call calc_Qx( i, j, north, Q1 )
      call calc_Qx( i, j, east,  Q2 )  

  elseif( i == 1 .and. j == jmax ) then
      ! corner = NorthWest
      call calc_Qx( i, j, east , Q1 )
      call calc_Qx( i, j, south, Q2 )  

  elseif( i == imax .and. j == 1 ) then
      ! corner = SouthEast
      call calc_Qx( i, j, west , Q1 )
      call calc_Qx( i, j, north, Q2 )  

  elseif( i == imax .and. j == jmax ) then
      ! corner = NorthEast
      call calc_Qx( i, j, south, Q1 )
      call calc_Qx( i, j, west,  Q2 )  

  endif 
  
end if


! At last, compute Qout 
Qout = (Vnew - Vold)/dt - (Q1 + Q2 + Q3 + Qrain)
!-------------------------------------------------------------------
end function F_Qout_bound
!====================================================================
!   \\\\\\\\\\    B E G I N         F U N C T I O N      ///////////
!   //////////          F _ Q O U T _ B O U N D          \\\\\\\\\\\
!====================================================================

END MODULE BoundCond

