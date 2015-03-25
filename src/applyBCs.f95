! fortran_free_source
! This module is part of PERFCODE, written by Bradley J Eck
!
module applyBCs
implicit none
contains
!
!   1.  Subroutine  SET_BC
!
!========================================================================
!   \\\\\\\\\\  B E G I N     S U B R O U T I N E    //////////
!   //////////         S E T  _ B C                  \\\\\\\\\\\
!=========================================================================
subroutine SET_BC( sol )                                   
!  Applies boundary conditions to the grid cell (i, j).
!       The aim is to reduce the number of lines dedicated to applying
!       boundary conditions and make the logic clear(er). 
!       'Apply boundary condtions' means compute the
!       appropriate parts of the matrix system for that cell.
!       
!
use shared
use pfc2dfuns
use utilities, only: F_Linterp
use ConvCoef
use pfc2dsubs
use BoundCond
implicit none
! Arguments
integer, intent(in) :: sol  ! which solution are we talking about here?
! Internal variables
integer :: side, corner
character(len=7) :: condition 
character(len=7) :: EW_bc, NS_bc  ! BC for East/West or North/South side of corner cell
real :: h_curb, v_curb !depth and velocity into curb inlet
!------------------------------------------------------------------------
! NO CORNER 
!    First deal with boundary cells that are not corners
no_corner: if( is_corner(i,j) .eqv. .false. ) then

!What side of the domain and boundary condition?
if( j .eq. 1 ) then
  side      = south
  condition = south_bc(i)

elseif( j .eq. jmax ) then
  side      = north
  condition = north_bc(i)

elseif( i .eq. 1 ) then
  side      = west
  condition = west_bc(j)

elseif( i .eq. imax ) then
  side      = east
  condition = east_bc(j)

end if

! Some conditions contribute to both the stationary (old) 
!   and iterative (itr) while others do not
if( sol .eq. old ) then
  ! Implement BCs condition by condition
  if( condition .eq. 'NO_FLOW' ) then
        
      if( side .eq. north ) then
        Cn = 0.0
        CALL Conveyance( west , old, i, j, Cw )
        CALL Conveyance( east , old, i, j, Ce )
        CALL Conveyance( south, old, i, j, Cs )

      elseif( side .eq. south ) then 
        Cs = 0.0
        CALL Conveyance( west , old, i, j, Cw )
        CALL Conveyance( east , old, i, j, Ce )
        CALL Conveyance( north, old, i, j, Cn )

      elseif( side .eq. east ) then
        Ce = 0.0
        CALL Conveyance( west , old, i, j, Cw )
        CALL Conveyance( south, old, i, j, Cs )
        CALL Conveyance( north, old, i, j, Cn )

      elseif( side .eq. west ) then
        Cw = 0.0
        CALL Conveyance( east , old, i, j, Ce )
        CALL Conveyance( south, old, i, j, Cs )
        CALL Conveyance( north, old, i, j, Cn )

      endif

      ! Compute the old parts of the linear system
      pf    = F_por( h_old( i, j ), b_pfc(i, j), por(i, j) )    
      v     = F_LinearIndex( i, j, jmax)
      Fn(v) = F_RHS_n( i, j, Cw, Ce, Cs, Cn, rain(n-1), pf, dt )

  elseif( condition .eq. '1D_FLOW' ) then
      
      v = F_LinearIndex ( i, j, jmax)
      eta_1D = F_LINTERP(       X = CV_Info( v ) % eta ,  &
                          known_X = eta_cs             ,  &
                          known_Y = eta_cs_1D          ,  &
                                n = nr_cs + 1                  ) 
      h_bound= F_LINTERP(       X = eta_1D             ,  &
                          known_X = etaCV              ,  &
                          known_Y = h_new_1D           ,  &
                                n = TNE                        )
      C(v) = 1.0
      F(v) = h_bound

  elseif( condition .eq. 'MOC_KIN' ) then
      CALL MOC_KIN_BC( i, j, rain(n), dt, side, h_bound, 100 )
      v = F_LinearIndex( i, j, jmax )
      C(v) = 1.
      F(v) = h_bound
       
  elseif( ( side .eq. west) .and.      &
          ( condition .eq. 'eastKIN' ) ) then
      call MOC_KIN_BC( imax, j, rain(n), dt, east, h_bound, 100 )
      v = F_LinearIndex( i, j, jmax )
      C(v) = 1.
      F(v) = h_bound  

  elseif(  ( side .eq. east ) .and.      &
           ( condition .eq. 'westKIN' )  ) then
      call MOC_KIN_BC( 1, j, rain(n), dt, west, h_bound, 100 )
      v = F_LinearIndex( i, j, jmax )
      C(v) = 1.
      F(v) = h_bound    

  elseif( (side .eq. south) .and. (condition .eq. 'CURB_IN') ) then

    if( h_old(i,j) .lt. b_pfc(i,j) ) then
        ! water below the pfc use MOC KIN
      CALL MOC_KIN_BC( i, j, rain(n), dt, side, h_bound, 100 )
      v = F_LinearIndex( i, j, jmax )
      C(v) = 1.
      F(v) = h_bound

    else  ! water above the pfc, use curb inlet

      Cs = 0.0 ! since we pass this to the following function
      CALL Conveyance( west , old, i, j, Cw )
      CALL Conveyance( east , old, i, j, Ce )
      CALL Conveyance( north, old, i, j, Cn )
      pf    = F_por( h_old( i, j ), b_pfc(i, j), por(i, j) )    
      Fn(v) = F_RHS_n_CURB_IN( i, j, Cw, Ce, Cs, Cn, rain(n-1), pf, dt, south )       

    endif

  endif 

          
elseif( sol .eq. itr ) then

   ! Only NO_FLOW and CURB_IN boundaries have an itr part...use similar logic as above
  if( condition .eq. 'NO_FLOW' ) then
        
      if( side .eq. north ) then
        Cn1 = 0.0
        CALL Conveyance( west , itr, i, j, Cw1 )
        CALL Conveyance( east , itr, i, j, Ce1 )
        CALL Conveyance( south, itr, i, j, Cs1 )

      elseif( side .eq. south ) then 
        Cs1 = 0.0
        CALL Conveyance( west , itr, i, j, Cw1 )
        CALL Conveyance( east , itr, i, j, Ce1 )
        CALL Conveyance( north, itr, i, j, Cn1 )

      elseif( side .eq. east ) then
        Ce1 = 0.0
        CALL Conveyance( west , itr, i, j, Cw1 )
        CALL Conveyance( south, itr, i, j, Cs1 )
        CALL Conveyance( north, itr, i, j, Cn1 )

      elseif( side .eq. west ) then
        Cw1 = 0.0
        CALL Conveyance( east , itr, i, j, Ce1 )
        CALL Conveyance( south, itr, i, j, Cs1 )
        CALL Conveyance( north, itr, i, j, Cn1 )

      endif

      pf = F_por( h_itr(i, j) , b_pfc(i, j), por(i, j))
      CALL set_ABCDEF( i, j, Cw1, Ce1, Cs1, Cn1, pf, dt, rain(n) )

  elseif( ( condition .eq. 'CURB_IN' ) .and. &
          ( h_itr(i,j)    .gt. b_pfc(i,j))          ) then
    
    if(  side .eq. south ) then

        Cs1 = 0.0
        CALL Conveyance( west , itr, i, j, Cw1 )
        CALL Conveyance( east , itr, i, j, Ce1 )
        CALL Conveyance( north, itr, i, j, Cn1 )

        v   = F_LinearIndex( i, j, jmax)

        h_curb = h_itr(i,j) - b_pfc(i,j) + ( Z(i,j+1) - Z(i,j) ) / 2.
        V_curb = sqrt( 2. * g * h_curb / 3. )

        ! Bands of penta-diagonal matrix
        A(v) = - dt / 2. * pf * Cw1
        B(v) = 0.0 
        C(v) =   dt / 2. * pf * ( Cw1 + V_curb * lng_south(i,j)*2./3. + Cn1 + Ce1 ) + 1.
        D(v) = - dt / 2. * pf * Cn1
        E(v) = - dt / 2. * pf * Ce1
        ! Right-hand-side
        ! Portion from time level n+1
        F1(v) = F_RHS_n1_CURB_IN( i, j, Cw1, Ce1, Cs1, Cn1, rain(n), pf, dt, south )
        !The complete right hand side has contributions from 
        ! time level n and time level n+1
        F(v) = Fn(v) + F1(v)

    end if
    

  endif

end if

end if no_corner

!-------------------------------------------
! CORNER POINTS
! Boundary conditions are a real pain at the domain corners
! because there are many possible combinations and the number
! of possibilities grows each time a new condition is added.
! Further, we implement essentially two types of BCs:
!     1) Dirichlet where the answer at the cell center is
!        specified.  This is the 1D_FLOW or MOC_KIN condition.
!     2) Neumann(ish) where the BC is applied on the conveyance
!        coefficient.  In the case of a NO_FLOW boundary
!        the derivative is prescribed and it is a true Neumann 
!        condition.  In the case of CURB_IN we are specifying
!        the depth at the boundary and not the cell center 
!        and this is dealt with in the conveyance coef.
!
! In applying BCs at the corners we take the answer at the cell
! center if its available with 1D_FLOW gets precedence over 
! MOC_KIN. Otherwise we compute the conveyance
! coefficients for each face and solve the system as usual.
 
if( is_corner(i,j) .eqv. .true. ) then

! Which corner are we at?
! And what condition applies to the face east or west (EW) face
! and the north or south (NS)  
if( i == 1 .and. j == 1 ) then
    corner = SouthWest
    NS_bc = south_bc( i )
    EW_bc =  west_bc( j )        

elseif( i == 1 .and. j == jmax ) then
    corner = NorthWest
    NS_bc = north_bc( i )
    EW_bc =  west_bc( j )        

elseif( i == imax .and. j == 1 ) then
    corner = SouthEast
    NS_bc = south_bc( i )
    EW_bc =  east_bc( j )        

elseif( i == imax .and. j == jmax ) then
    corner = NorthEast
    NS_bc = north_bc( i )
    EW_bc =  east_bc( j )        

endif


!---------------------------------------------------------------
if( sol .eq. old ) then

  ! Take a soln at the cell center if its available
  if( EW_bc .eq. '1D_FLOW' ) then
      
      v = F_LinearIndex ( i, j, jmax)
      eta_1D = F_LINTERP(       X = CV_Info( v ) % eta ,  &
                          known_X = eta_cs             ,  &
                          known_Y = eta_cs_1D          ,  &
                                n = nr_cs + 1                  ) 
      h_bound= F_LINTERP(       X = eta_1D             ,  &
                          known_X = etaCV              ,  &
                          known_Y = h_new_1D           ,  &
                                n = TNE                        )
      C(v) = 1.0
      F(v) = h_bound

  elseif( NS_bc .eq. 'MOC_KIN' ) then !  .and.  &
!          EW_bc .eq. 'MOC_KIN'         ) then
    ! use the edge of pavement side
    if( corner .eq. southwest ) then 
        call MOC_KIN_BC( i, j, rain(n), dt, west, h_bound, 100 )

    elseif( corner .eq. southeast ) then
        call MOC_KIN_BC( i, j, rain(n), dt, east, h_bound, 100 )

    elseif( corner .eq. northwest ) then

        call MOC_KIN_BC( i, j, rain(n), dt, north, h_bound, 100 )

    elseif(  corner .eq. northeast ) then

        call MOC_KIN_BC( i, j, rain(n), dt, east, h_bound, 100 )

    end if
    
    v = F_LinearIndex( i, j, jmax )
    C(v) = 1.
    F(v) = h_bound            


  elseif( (NS_bc .eq. 'NO_FLOW') .and. &
          (EW_bc .eq. 'NO_FLOW')        ) then

      if( corner .eq. NorthEast ) then

          Cn = 0.0
          Ce = 0.0
          CALL Conveyance( west , old, i, j, Cw )
          CALL Conveyance( south, old, i, j, Cs )

      elseif( corner .eq. southeast ) then
        
          Cs = 0.0
          Ce = 0.0
          CALL Conveyance( west , old, i, j, Cw )
          CALL Conveyance( north, old, i, j, Cn )

      elseif( corner .eq. southwest ) then

          Cs = 0.0
          Cw = 0.0
          CALL Conveyance( east , old, i, j, Ce )
          CALL Conveyance( north, old, i, j, Cn )

      elseif( corner .eq. northwest ) then

          Cn = 0.0
          Cw = 0.0
          CALL Conveyance( east , old, i, j, Ce )
          CALL Conveyance( south, old, i, j, Cs )

      end if

      ! Compute the old parts of the linear system
      pf    = F_por( h_old( i, j ), b_pfc(i, j), por(i, j) )    
      v     = F_LinearIndex( i, j, jmax)
      Fn(v) = F_RHS_n( i, j, Cw, Ce, Cs, Cn, rain(n-1), pf, dt )
 

  elseif( (EW_bc .eq. 'eastKIN') .and. (i .eq. 1) ) then
    ! Figure the soln as the east side  
    call MOC_KIN_BC( imax, j, rain(n), dt, east, h_bound, 100 )
    v = F_LinearIndex( i, j, jmax )
    C(v) = 1.
    F(v) = h_bound


  elseif( EW_bc .eq. 'westKIN' .and. i .eq. imax ) then
    !Figure the soln as the west side
    if( corner .eq. northeast) then
        ! use one cell to the south
        call MOC_KIN_BC( 1, j-1, rain(n), dt, west, h_bound, 100 )
    else
        call MOC_KIN_BC( 1, j, rain(n), dt, west, h_bound, 100 )
    endif
    v = F_LinearIndex( i, j, jmax )
    C(v) = 1.
    F(v) = h_bound

!   elseif( NS_bc .eq. 'MOC_KIN' .and.  EW_bc .eq. 'NO_FLOW' )  then
!       ! run moc kin for the north or south side
!       if( corner .eq. northwest .or. corner .eq. northeast ) then
!           call MOC_KIN_BC(i, j, rain(n), dt, north, h_bound, 100)
!       elseif( corner .eq. southwest .or. corner .eq. southeast) then
!           call MOC_KIN_BC(i, j, rain(n), dt, south, h_bound, 100)
!       end if
!          
!       v = F_LinearIndex( i, j, jmax )
!       C(v) = 1.
!       F(v) = h_bound

  elseif( EW_bc .eq. 'MOC_KIN' ) then
      ! run moc kin for the east or west side
      if( corner .eq. southwest ) then
          call MOC_KIN_BC(i, j, rain(n), dt, west, h_bound, 100)
      elseif( corner .eq. northwest) then
          ! use one cell to the south for this corner
          call MOC_KIN_BC(i, j-1, rain(n), dt, west, h_bound, 100)
      elseif( corner .eq. northeast .or. corner .eq. southeast) then
          call MOC_KIN_BC(i, j, rain(n), dt, east, h_bound, 100)
      end if
       
      v = F_LinearIndex( i, j, jmax )
      C(v) = 1.
      F(v) = h_bound
  
  end if

!---------------------------------------------------------------
elseif( sol .eq. itr ) then

    if( (NS_bc .eq. 'NO_FLOW') .and. &
        (EW_bc .eq. 'NO_FLOW')       ) then

      if( corner .eq. NorthEast ) then

        Cn1 = 0.0
        Ce1 = 0.0
        CALL Conveyance( west , itr, i, j, Cw1 )
        CALL Conveyance( south, itr, i, j, Cs1 )

      elseif( corner .eq. southeast ) then
      
        Cs1 = 0.0
        Ce1 = 0.0
        CALL Conveyance( west , itr, i, j, Cw1 )
        CALL Conveyance( north, itr, i, j, Cn1 )

      elseif( corner .eq. southwest ) then

        Cs1 = 0.0
        Cw1 = 0.0
        CALL Conveyance( east , itr, i, j, Ce1 )
        CALL Conveyance( north, itr, i, j, Cn1 )

      elseif( corner .eq. northwest ) then

        Cn1 = 0.0
        Cw1 = 0.0
        CALL Conveyance( east , itr, i, j, Ce1 )
        CALL Conveyance( south, itr, i, j, Cs1 )

      end if

      pf = F_por( h_itr(i, j) , b_pfc(i, j), por(i, j))
      CALL set_ABCDEF( i, j, Cw1, Ce1, Cs1, Cn1, pf, dt, rain(n) )

   end if

end if

end if 

end subroutine SET_BC
!========================================================================
!   \\\\\\\\\\  E N D         S U B R O U T I N E    //////////
!   //////////         S E T _ B C                   \\\\\\\\\\\
!=========================================================================

end module
