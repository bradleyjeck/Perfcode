! fortran_free_source


module ConvCoef

implicit none

contains


!   1. CONVEYANCE
!   2. FrictionSlope

!=======================================================================
!   \\\\\\\\\\  B E G I N     S U B R O U T I N E //////////
!   //////////      C O N V E Y A N C E           \\\\\\\\\\
!=======================================================================
SUBROUTINE CONVEYANCE( face, sol, i, j, CC )          
!   This subroutine computes the conveyance coefficient
!               for a given cell face...look out, its fancy!
!   HISTORY:    March 2010 -- Original coding
!               July 2011  -- Modified to use harmonic average K & n_mann
!                             in case values differ between cells.  
!                             Also added a check for outcell index to 
!                             allow calling at a boundary face
USE shared, ONLY: Sfw_old, Sfe_old, Sfs_old, Sfn_old,     &
                  Sfw_itr, Sfe_itr, Sfs_itr, Sfn_itr,     &
                  h_old, h_itr, wid, Z, hyd_cond, n_mann, & 
                  b_pfc, lng, lng_south, lng_north,       &
                  h_pfc_min, area,                        &
                  north, south, east, west, old, itr
USE utilities, ONLY: F_HarmonicAvg

implicit none
!VARIABLE DECLARATIONS
!   Arguments
integer, intent( in )  :: face     ! Which face?
integer, intent( in)   :: sol      ! Computed based on which solution?
integer, intent( in )  :: i, j     ! of which cell?
REAL, intent(out)      :: CC       ! the conveyance coefficient
!   Internal Variables
real :: K_equiv  ! equivalent hydraulic conductivity
real :: hp, hs   ! the thickness in the pavement and on the surface
REAL :: distin, distout  ! size of cells for scaling purposes (will be length or width depeding on which direction we're going.
REAL :: fluxdist    ! the distance (size) of the cell face that the flux applies t
REAL :: hin, hout    ! thickness at CV center
REAL :: zin, zout    ! elevation at CV center
REAL :: head_at_face, Zface !HEAD and ELEVATION at the face
REAL, POINTER, DIMENSION(:,:) :: h  ! pointer to the thickness array
REAL, POINTER, DIMENSION(:,:) :: Sfw, Sfe, Sfs, Sfn ! points to magnitude of friction slope at compass face.
REAL :: Sf !the friction slope for the particular face that we're working with
logical :: error !make sure the result is reasonable
!-------------------------------------------------------------------------
! Compute based on the old or iterative thickness?
if ( sol .EQ. old ) then
    ! thickness array
    h => h_old
    ! friction slope arrays
    Sfw => Sfw_old
    Sfe => Sfe_old
    Sfs => Sfs_old
    Sfn => Sfn_old
elseif( sol .eq. itr ) then
    h => h_itr
    Sfw => Sfw_itr
    Sfe => Sfe_itr
    Sfs => Sfs_itr
    Sfn => Sfn_itr
endif

! set internal/generic variables based on cell face
if ( face .EQ. west ) then
    distin = lng( i, j)
    distout= lng( i-1, j)
    hin = h(i,j)
    hout= h(i-1,j)
    zin = Z(i,j)
    zout= Z(i-1,j)
    fluxdist = wid(i,j)
    Sf = Sfw( i, j)
    K_equiv = F_HarmonicAvg( x=(/ hyd_cond(i,j), hyd_cond(i-1,j) /), imax=2 )
elseif( face .eq. east ) then
    distin = lng( i, j)
    distout= lng( i+1, j)
    hin = h(i,j)
    hout= h(i+1,j)
    zin = Z(i,j)
    zout= Z(i+1,j)
    fluxdist = wid(i,j)
    Sf = Sfe( i, j)
    K_equiv = F_HarmonicAvg( x=(/ hyd_cond(i,j), hyd_cond(i+1,j) /), imax=2 )
elseif( face .EQ. south ) then
    distin = wid( i, j)
    distout= wid( i, j-1)
    hin = h(i,j)
    hout= h(i,j-1)
    zin = Z(i,j)
    zout= Z(i,j-1)
    fluxdist = lng_south(i,j)
    Sf = Sfs(i,j)
    K_equiv = F_HarmonicAvg( x=(/ hyd_cond(i,j), hyd_cond(i,j-1) /), imax=2 )
elseif( face .EQ. north ) then
    distin = wid( i, j)
    distout= wid( i, j+1)
    hin = h(i,j)
    hout= h(i,j+1)
    zin = Z(i,j)
    zout= Z(i,j+1)
    fluxdist = lng_north(i,j)
    Sf = Sfn(i,j)
    K_equiv = F_HarmonicAvg( x=(/ hyd_cond(i,j), hyd_cond(i,j+1) /), imax=2 )
endif
!Compute the total head at the cell face
head_at_face = ( (hin+zin)*distout + (hout+zout)*distin )  &
                        / ( distin + distout)
!Elevation at the cell face
Zface        = ( zin*distout + zout*distin ) / ( distin + distout)
!compute the thicknesses
hp = MIN ( b_pfc(i,j), head_at_face - Zface             )
hs = MAX ( 0.        , head_at_face - Zface - b_pfc(i,j))


!Force hp to stay positive
if( hp .LT. 0.0 ) then
    hp = TINY(h_pfc_min)
end if

! Compute the Conveyance coefficient
! would really like to just one statement to calc the conv coef
!   but sqrt(Sf) sometimes gives problems, even when there is no
!   sheet flow, so this if block hopefully avoids the problem

if( hs .GT. 0.) then
    !Sheet flow occurs and compute CC as usual
    CC = ( K_equiv * hp + 1./n_mann(i,j)*hs**(5./3.)/sqrt(Sf) )  *  &
            ( 2.*fluxdist / ( distout + distin ) ) / Area(i,j)
else
    !Sheet flow does not occur and CC only depends on subsurface
    CC = ( K_equiv * hp   )  *  &     
            ( 2.*fluxdist / ( distout + distin ) ) / Area(i,j)
end if


! Set CC = 0 if Sf ==0
if( Sf == 0. ) then
    CC =  0.
end if


! ERROR CHECKING FOR CONVEYANCE COEFS
if( ( CC .GT. HUGE(CC) ) .OR. ( CC .LT. -HUGE(CC) ) ) then
    error = .true.
else
    error = .false.
endif

!Output the parts of the calculation if the error is true
if( error .eqv. .true. ) then
    write(*,*) 'Problem with conveyance coefficient!'
    print *, 'i = ', i, ' j = ', j, ' Face = ', face, ' Soln = ', sol
    print *, '       K = ', K_equiv
    print *, '      hp = ', hp
    print *, '  n_mann = ', n_mann(i,j)
    print *, '      hs = ', hs
    print *, '      Sf = ', Sf
    print *, 'fluxdist = ', fluxdist
    print *, ' distout = ', distout
    print *, '  distin = ', distin
    print *, '    Area = ', Area(i,j)
    print *, '      CC = ', CC
    write(*,*) 'Stopping Program'
    STOP
endif

END subroutine conveyance

!=======================================================================
!   \\\\\\\\\\  E N D         S U B R O U T I N E //////////
!   //////////      C O N V E Y A N C E           \\\\\\\\\\
!=======================================================================

!=======================================================================
!   \\\\\\\\\\  B E G I N     S U B R O U T I N E //////////
!   //////////     F R I C T I O N   S L O P E    \\\\\\\\\\
!=======================================================================
!   Purpose:    This subroutine computes the magnitude of the friction
!               slope at the cell faces.
!               The arguments specifcy whether to use the OLD or ITR
!               solution array in the calculations and the arrays 
!               for storing the results.
!
!
!       ---x---|---x---     Key:  * is CV Center
!       |      |      |           x normal component of friction slope
!       |  *   O   *  |             computed here by central difference
!       |      |      |           O the four normal components are
!       ---x---|---x---             tangent here and so are averaged
!
SUBROUTINE FrictionSlope( sol, Sfw, Sfe, Sfs, Sfn )
use SHARED,    only: h_old, h_itr, Z, lng, wid, imax, jmax, old, itr, dry                  
use outputs,   only: write_flipped_matrix
use utilities, only: F_PythagSum, F_Extrapolate
!-------------------------------------------------------------------------
!VARIABLE DECLARATIONS
implicit none
!   Arguments
integer, intent( in ) :: sol
REAL, DIMENSION(imax,jmax), intent(out), optional :: Sfw, Sfe, Sfs, Sfn
!   Internal Variables
REAL, DIMENSION(imax, jmax) :: HD      ! Total HEAD at cell centers
REAL, DIMENSION(:,:), pointer :: h     ! Pointer to array of thicknesses
REAL, DIMENSION(imax,jmax) :: Sf_norm_west, Sf_norm_east, Sf_norm_south, Sf_norm_north
REAL, DIMENSION(imax,jmax) :: Sf_tan_west, Sf_tan_east, Sf_tan_south, Sf_tan_north
REAL, ALLOCATABLE, DIMENSION(:,:), target :: h_dry  ! for computing pavement slopes
INTEGER :: i, j         !array idices

!----------------------------------------------------------------------------
! choose which thickness array to use for estimating the friction slope
if ( sol .eq. old ) then
    h => h_old
elseif( sol .eq. itr ) then
    h => h_itr
elseif( sol .eq. dry ) then
    allocate( h_dry( imax, jmax ) )
    h_dry = 0.0
    h => h_dry
endif

! compute the total head
HD = h + z !total head is thickness plus elevation

!initialize arrays to zero
!       Sf_norm_west = 0.0
!       Sf_norm_east = 0.0
!       Sf_norm_south= 0.0
!       Sf_norm_north= 0.0
!       Sf_tan_west = 0.0
!       Sf_tan_east = 0.0
!       Sf_tan_south= 0.0
!       Sf_tan_north= 0.0
!
!------------------------------------------------------------------------
!   C O M P O N E N T   N O R M A L   T O   E A C H    F A C E 
!------------------------------------------------------------------------

! 

! WEST
do j = 1, jmax
    ! Domain interior, by central differences 
    do i = 2, imax
        Sf_norm_west(i,j) = ( ( HD(i,j) - HD(i-1,j)  ) / 0.5 / ( lng(i-1,j) + lng(i,j) ) )
    end do
    ! Western boundary of domain by extrapolation
    i = 1
    Sf_norm_west(i,j) = F_Extrapolate( 0.,                                     &
                                     lng(i,j)           , Sf_norm_west(i+1,j), &
                                     lng(i,j)+lng(i+1,j), Sf_norm_west(i+2,j)  )
end do

!EAST
do j = 1, jmax
    ! Domain interior, by central differences
    do i = 1, imax - 1
        Sf_norm_east(i,j) = (  ( HD(i+1,j) - HD(i,j) ) / 0.5 / ( lng(i+1,j) + lng(i,j) ) )
    end do
    ! Eastern boundary of domain by extrapolation
    i = imax
    Sf_norm_east(i,j) = F_Extrapolate( 0.,                                     &
                                     lng(i,j)           , Sf_norm_east(i-1,j), &
                                     lng(i,j)+lng(i-1,j), Sf_norm_east(i-2,j)  )
end do


!SOUTH
do i = 1, imax
    ! Domain interior, by central differences
    do j = 2, jmax
        Sf_norm_south(i,j) = ( ( HD(i,j-1) - HD(i,j) ) / 0.5 / ( wid(i,j-1) + wid(i,j) ) )
    end do
    ! Southern boundary of domain by extrapolation
    j = 1
    Sf_norm_south(i,j) = F_Extrapolate( 0.,                                     &
                                     wid(i,j)           , Sf_norm_south(i,j+1), &
                                     wid(i,j)+wid(i,j+1), Sf_norm_south(i,j+2)  )

end do

!NORTH
do i = 1, imax
    ! Domain interior, by central differences
    do j = 1, jmax - 1
        Sf_norm_north(i,j) = (  ( HD(i,j) - HD(i,j+1) ) / 0.5 / ( wid(i,j+1) + wid(i,j) ) )
    end do
    ! Northern oundary of domain by extrapolation
    j = jmax
    Sf_norm_north(i,j) = F_Extrapolate( 0.,                                     &
                                     wid(i,j)           , Sf_norm_north(i,j-1), &
                                     wid(i,j)+wid(i,j-1), Sf_norm_north(i,j-2)  )
end do

!------------------------------------------------------------------------
!   C O M P O N E N T   T A N G E N T   T O   E A C H    F A C E 
!       A N D    M A G N I T U D E  A T   E A C H   F A C E
!------------------------------------------------------------------------
! component of friction slope that is TANGENT to each cell face
! computed by averaging the four nearest locations where the
! component is normal to a face.
!
!WEST
do j = 1, jmax  
    do i = 2, imax
        Sf_tan_west(i,j) =  ( (Sf_norm_north(i,j) + Sf_norm_south(i,j))*lng(i-1,j)   &
                             +(Sf_norm_north(i-1,j) + Sf_norm_south(i-1,j))*lng(i,j) ) & 
                            / ( 2. * ( lng(i,j) + lng(i-1, j) ) )
        Sfw(i,j) = F_PythagSum( Sf_norm_west (i,j), Sf_tan_west (i,j) )
    end do
end do 

!EAST
do j = 1, jmax 
    do i = 1, imax - 1
        Sf_tan_east(i,j) =  ( (Sf_norm_north(i,j) + Sf_norm_south(i,j))*lng(i+1,j)    &
                             +(Sf_norm_north(i+1,j) + Sf_norm_south(i+1,j))*lng(i,j) ) & 
                            / ( 2. * ( lng(i,j) + lng(i+1, j) ) )
        Sfe(i,j) = F_PythagSum( Sf_norm_east (i,j), Sf_tan_east (i,j) )
    end do
end do 

!SOUTH
do i = 1, imax
    do j = 2, jmax
        Sf_tan_south(i,j) =  ( ( Sf_norm_east(i,j) + Sf_norm_west(i,j) ) * wid(i,j-1)   &
                              +( Sf_norm_east(i,j-1)+Sf_norm_west(i,j-1))* wid(i,j  ) ) &
                            / ( 2. * ( wid(i,j) + wid(i,j-1) ) )
        Sfs(i,j) = F_PythagSum( Sf_norm_south(i,j), Sf_tan_south(i,j) )
    end do
end do 

!NORTH
do i = 1, imax
    do j = 1, jmax - 1
        Sf_tan_north(i,j) = ( ( Sf_norm_east(i,j) + Sf_norm_west(i,j) ) * wid(i,j+1)   &
                             +( Sf_norm_east(i,j+1)+Sf_norm_west(i,j+1))* wid(i,j  ) ) &
                            / ( 2. * ( wid(i,j) + wid(i,j+1) ) ) 
        Sfn(i,j) = F_PythagSum( Sf_norm_north(i,j), Sf_tan_north(i,j) )
    end do
end do


! deallocate space for h_dry
if( sol .eq. dry ) then
    deallocate( h_dry )
endif




!
!        i = 50; j = 51
!       write(100,*) 'i,j=', i, j
!       write(100,*) 'Sf_norm_east', Sf_norm_east(i,j)
!       write(100,*) 'Sf_tan_east', Sf_tan_east(i,j)
!       write(100,*) 'HD(i,j)', HD(i,j)
!       write(100,*) 'HD(i+1,j)', HD(i+1,j)
!       write(100,*) 'HD(i,j+1)', HD(i,j+1)
!       write(100,*) 'HD(i,j-1)', HD(i,j-1)
!       write(100,*) 'Sf_norm_north(i,j)', Sf_norm_north(i,j)
!       write(100,*) 'Sf_norm_north(i+1,j)', Sf_norm_north(i+1,j)
!       write(100,*) 'Sf_norm_south(i,j)', Sf_norm_south(i,j)
!       write(100,*) 'Sf_norm_south(i+1,j)', Sf_norm_south(i+1,j)
!
!

end subroutine FrictionSlope

!=======================================================================
!   \\\\\\\\\\     E N D      S U B R O U T I N E //////////
!   //////////     F R I C T I O N   S L O P E    \\\\\\\\\\
!=======================================================================





end module ConvCoef

