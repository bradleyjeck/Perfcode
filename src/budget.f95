! fortran_free_source
!  
!  This module is part of Perfcode, written by Bradley J Eck
!
module budget
implicit none
contains
!
!	1. Subroutine WATER_BUDGET
!
!========================================================================
!   \\\\\\\\\\  B E G I N     S U B R O U T I N E    //////////
!   //////////     W A T E R _ B U D G E T           \\\\\\\\\\\
!=========================================================================
subroutine WATER_BUDGET(     old_volume, &
                             pfc_volume, & 
                         surface_volume, &
                            rain_volume, &
                       north_out_volume, &
                       south_out_volume, &
                        east_out_volume, &
                        west_out_volume, & 
                        net_flow_volume, &
                     storage_chg_volume, &
                          volume_error , &
                    volume_error_fraction       )
!   Computes a global water balance at each time step.
use shared, only: imax, jmax, area, por, b_pfc, h_new, rain, n, dt
use boundcond, only: F_Qout_bound
use pfc2Dsubs, only: CALC_Qx
implicit none
! Arguments
real, intent(in)  :: old_volume
real, intent(out) :: pfc_volume, surface_volume, rain_volume
real, intent(out) :: north_out_volume, south_out_volume
real, intent(out) ::  east_out_volume, west_out_volume
real, intent(out) :: net_flow_volume,  storage_chg_volume  
real, intent(out) :: volume_error, volume_error_fraction      
! Internal variables
integer :: i, j, v
real :: hp, hs
real, dimension(imax, jmax) :: vol  !volume array
!-----------------------------------------------------------------------
vol(:,:) = 0.0
! WATER BALANCE
! Volume of water in the pfc 
do j = 1, jmax
    do i = 1, imax
        hp = min( h_new(i, j), b_pfc(i, j) )
        vol(i, j) = hp * area(i,j) * por(i, j) 
    end do
end do

pfc_volume = sum( vol(:,:) )

! Volume of water on the surface
vol(:,:) = 0.0

do j = 1, jmax
    do i = 1, imax
        hs = max( h_new(i, j) - b_pfc(i, j), 0.0 ) 
        vol(i, j) = hs * area(i,j)  
    end do
end do

surface_volume = sum( vol(:,:) )

! Rainfall
rain_volume = rain(n) * dt * sum( area(:,:) )

! Outflow
vol(:,:) = 0.
! North & South (including all of corner cells)
do j = 1, jmax, jmax-1
    do i = 1, imax

        vol(i,j) = F_Qout_bound(i, j) * dt

    end do
end do

south_out_volume = sum( vol(:,1)    )
north_out_volume = sum( vol(:,jmax) )

! East & West    
do j = 2, jmax-1
    do i = 1, imax, imax-1

        vol(i,j) = F_Qout_bound(i, j) * dt
    
    end do
end do

west_out_volume = sum( vol( 1   , 2:jmax-1 ) )
east_out_volume = sum( vol( imax, 2:jmax-1 ) )

! Net flow
net_flow_volume = rain_volume + north_out_volume + south_out_volume & 
                              +  east_out_volume +  west_out_volume

! Storage Change
storage_chg_volume = pfc_volume + surface_volume - old_volume

! Error
volume_error = net_flow_volume - storage_chg_volume
volume_error_fraction = volume_error / storage_chg_volume

!-----------------------------------------------------------------------
end subroutine water_budget
!=======================================================================
!  \\\\\\\\\\   E N D        S U B R O U T I N E    //////////
!  //////////     W A T E R _ B U D G E T           \\\\\\\\\\\
!========================================================================     

end module 
