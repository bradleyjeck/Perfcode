! fortran_free_source
!
! (c) Copyright 2010, 2015 Bradley J. Eck
! This module is part of PERFCODE 
!

!============================================================================
!   \\\\\\\\\\                                        //////////
                        MODULE BoundCond 
!   //////////                                        \\\\\\\\\\
                        implicit none
                        contains
!============================================================================

!===================================================================
!   \\\\\\\\\\       B E G I N    S U B R O U T I N E ///////////
!   //////////          M O C _ K I N _ B C           \\\\\\\\\\\
!==================================================================

!   inputs:  everything
!  outputs:  the depth in the boundary cell


subroutine MOC_KIN_BC( i, j, rain, dt, side,  h_bound, dev)

use shared   , only: K, por, b_pfc, n_mann, CV_Info, wid, &
                     imax, jmax, lng, wid, h_old, Z, eta_0_hp2_max
use pfc2dsubs, only: set_xyh
use pfc2Dfuns, only: F_LinearIndex
use utilities, only: BILINEAR_INTERP, F_PythagSum 


integer, intent( in ) :: i, j
real, intent ( in )   :: dt      ! timestep
character(5), intent(in) :: side ! which side of the domain are we working on
real, intent( in ) :: rain       ! rainfall rate for this timestep 
real, intent( out )   :: h_bound
integer, optional :: dev  !device for outputing errors 

real, dimension( 2 ) :: ksi_ii1  !vector in the ksi direction from point i to point i+1
real, dimension( 2 ) :: eta_jj1  !vector in the eta direction from point j to point j+1
real, dimension( 2 ) :: S_ksi    ! slope vector in the ksi direction
real, dimension( 2 ) :: S_eta    ! slope vector in the eta direction
real, dimension( 2 ) :: S_drain
real, dimension( 2 ) :: S_drain_unit  ! slope vector for drainage slope
real :: drain_slope  ! magnitude of drainge slope

integer :: v
integer :: vi1  !value of v for the cell i+1
integer :: vj1  !value of v for the cell j+1
integer :: vjm1 !value of v for the cell j-1
integer :: vim1 !value of v for the cell i-1
integer :: v1, v2, v3, v4  ! global index for interpolation points
! Get a vector that points up the drainage slope from the point i,j


! Bilinear Interpolation
real :: XX, YY       ! Coordinates of point where depth is interpolated
real :: x1, y1, h1   ! Coordinates of point 1  Interpolation points
real :: x2, y2, h2   !   "    "       point 2
real :: x3, y3, h3   !   "    "       point 3
real :: x4, y4, h4   !   "    "       point 4

! Method of Characteristics
!   PFC
real :: dx_moc, hp1, hp2, hp2_max
!   Sheet flow
real :: ds, hs1, hs2


integer :: device
logical :: bilin_err

!-------------------------------------------------------------------------------------------

! Default values for output device
if( present( dev ) .EQV. .FALSE. ) then
        device = 6
else
        device = dev
end if


!----------------------------------------------------------------------
!   D R A I N A G E   S L O P E    C A L C U L A T I O N S 
!----------------------------------------------------------------------

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
if( side == 'north' .or. side == 'south' ) then
    ksi_ii1 = (/  CV_info(vi1)%X - CV_Info(v)%X ,  &
                  CV_info(vi1)%Y - CV_Info(v)%Y      /)

elseif( side == 'east ' )  then
    ksi_ii1 = (/  CV_info(v)%X - CV_Info(vim1)%X ,  &
                  CV_info(v)%Y - CV_Info(vim1)%Y      /)

endif


! TRANSVERSE DIRECTION (eta)

if( side == 'south' ) then

        eta_jj1 =   (/  CV_info(vj1)%X - CV_Info(v)%X ,  &
                        CV_info(vj1)%Y - CV_Info(v)%Y      /)

elseif( side == 'north' ) then

        eta_jj1 = - (/  CV_info(vjm1)%X - CV_Info(v)%X,  &
                        CV_info(vjm1)%Y - CV_Info(v)%Y     /)


elseif( side == 'east ' ) then

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
if(  side == 'north' .or. side == 'south' ) then

    S_ksi = ksi_ii1 * ( Z( i+1, j   ) - Z( i, j ) ) / lng( i, j )

elseif( side == 'east ' ) then

    S_ksi = ksi_ii1 * ( Z( i, j   ) - Z( i-1, j ) ) / lng( i, j )

end if


! TRANSVERSE DIRECTION (eta)
if( side == 'south') then

    S_eta =   eta_jj1 * ( Z( i  , j+1 ) - Z( i, j ) ) / wid( i, j )

elseif( side == 'north' ) then

    S_eta = - eta_jj1 * ( Z( i  , j-1 ) - Z( i, j ) ) / wid( i, j )

elseif( side == 'east ' ) then

    if( j /= jmax ) then
        
        S_eta =   eta_jj1 * ( Z( i  , j+1 ) - Z( i, j ) ) / wid( i, j )

    elseif( j == jmax ) then

        S_eta =  eta_jj1 * (  Z( i, j) - Z( i, j-1) ) / wid( i, j )

    endif

end if


!--------------------------------------------------
! Compute vector for the drainage slope ( S_drain )
! and its magnitude, and unit vector for direction
!-------------------------------------------------


! compute drainage slope vector
S_drain = S_ksi + S_eta
! and the magnitude
drain_slope = F_PythagSum( S_drain(1), S_drain(2) )
! and a drainage slope unit vector
S_drain_unit = S_drain / drain_slope


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

if    ( side == 'south'  .AND.  S_drain_unit(1) .LE. 0. )  then
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

elseif( side == 'south'  .AND.  S_drain_unit(1) .GE. 0.0 ) then
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

elseif( side == 'north' .AND. S_drain_unit(1) .LE. 0.0 ) then
        ! This is the northern boundary and 
        ! The domain slopes from left to right
        call set_xyh( i-1, j-1, x1, y1, h1 )
        call set_xyh( i  , j-1, x2, y2, h2 )
        call set_xyh( i  , j  , x3, y3, h3 )
        call set_xyh( i-1, j  , x4, y4, h4 )

elseif( side == 'north' .AND. S_drain_unit(1) .GE. 0.0 ) then
        ! This is the northen boundary and 
        ! The domain slopes from right to left
        call set_xyh( i  , j-1, x1, y1, h1 )
        call set_xyh( i+1, j-1, x2, y2, h2 )
        call set_xyh( i+1, j  , x3, y3, h3 )
        call set_xyh( i  , j  , x4, h4, h4 )


elseif( side == 'east ' .AND. S_drain_unit(2) .GE. 0.0 ) then
        ! This is the eastern boundary and
        ! and uphill is the positive Y direction
        call set_xyh( i-1, j  , x1, y1, h1 )
        call set_xyh( i  , j  , x2, y2, h2 )
        call set_xyh( i  , j+1, x3, y3, h3 )
        call set_xyh( i-1, j+1, x4, y4, h4 )


elseif( side == 'east ' .AND. S_drain_unit(2) .LT. 0.0 ) then
        ! This is the eastern boundary and
        ! and uphill is the negative Y direction
        call set_xyh( i-1, j-1, x1, y1, h1 )
        call set_xyh( i  , j-1, x2, y2, h2 )
        call set_xyh( i  , j  , x3, y3, h3 )
        call set_xyh( i-1, j  , x4, y4, h4 )

endif




!----------------------------------------------------------------------
!   M E T H O D     O F     C H A R A C T E R I S T I C S 
!----------------------------------------------------------------------


! Reset v to confirm we're in the right cell
v = F_LinearIndex( i, j, jmax )


MOC:if( h_old(i,j) .LE. b_pfc) then
!------------------
!PFC FLOW
!-----------------
        ! Sheet flow has not started yet
        ! use MOC to estimate the solution at the next time step
        ! figure out how far up the drainage slope to go
        dx_moc =  K * (drain_slope) * dt / por
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
        hp2 = hp1 + rain * dt / por
        ! set maximum value for hp2 (1D flow)      
        if( rain .LT. TINY ( rain )) then
                ! Rainfall rate is effectively zero
                hp2 = hp2   ! Eqv to hp2 = hp1 
        else
                ! Rainfall is non-zero, set a maximum value for hp2
                hp2 = min( hp2, sum(wid(i,:))*rain/K/drain_slope)
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
!                               ' drain_slope =', drain_slope, &
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
                                ' drain_slope=', drain_slope, &
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






        if( i == imax/2  .OR. j == jmax/2  ) then
           write(device,*) 'PFC Flow MOC BC: i=', i     , &
                                            'j=', j     , &
                                          'hp2=', hp2   , &
                                       'dx_moc=', dx_moc, &
                                          'hp1=', hp1   , &
                                         'rain=', rain  , &
                                           'dt=', dt    , &
                                          'por=', por
        endif
        h_bound = hp2
else
!-----------------
!SHEET FLOW
!-----------------
        hs2 =  h_old(i,j) - b_pfc
        ! Handle Zero Rainfall        
        if( rain .LT. TINY( rain ) ) then
            ! there is no increase in flow rate along the drainage path
            !  and ds becomes arbitrary so use the characteristic length for PFC flow
            ! b/c you might need it later
            ds = K * ( drain_slope ) * dt / por
        else
            ds = sqrt( drain_slope ) / n_mann / rain *      &
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
        hs1 = hs1 - b_pfc
        !Handle return to PFC flow
        if( hs1 .GT. 0. ) then
            !we have sheet flow
            !Output some summary info
            if( i == imax/2 .or. j == jmax / 2 ) then
               write(device,*) 'Sheet Flow MOC BC: i=',   i, &
                                                  'j=',   j, &
                                                'hs2=', hs2, &
                                                 'ds=',  ds, &
                                                'hs1=',  hs1
            endif
            !checking for good values of inputs
            if( hs1 .LT. 0.  .OR.  hs2 .LT. 0.  .OR. rain .LT. 0.) then
                write(device,*) 'Sheet Flow MOC BC: i=',i, 'j=',j, &
                                 'hs1=', hs1, 'hs2=',hs2, 'rain=',rain
            end if
            ! return value for the boundary
            h_bound =  b_pfc + (   hs1**(5./3.) +             & 
                                 ( hs2 + rain*dt )**(5./3.) - &
                                   hs2**(5./3.) )**0.6
        else
            ! the upstream point does not have sheet flow
            ! use PFC characterisitic
            h_bound = hs1 + b_pfc + rain * dt / por  
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





!============================================================================
!   \\\\\\\\\                                       ///////////
                       END MODULE BoundCond
!   /////////                                       \\\\\\\\\\\\
!============================================================================

