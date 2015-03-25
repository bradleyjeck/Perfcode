!fortran_free_source

! need a seperate module for these last two functions b/c
! the function F_CC calls both of them, and they cannot be 
! in the same module
module pfc1dfuns2
implicit none
contains

!==========================================================================
!function to determine thickness in the pavement at the cell face
FUNCTION F_hp_face(dxin, hin, zin, dxout, hout, zout, b)
implicit none
!INPUTS
REAL :: dxin, dxout  ! size of cells
REAL :: hin, hout    ! thickness at CV center
REAL :: zin, zout    ! elevation at CV center
REAL :: b            ! pavement thickness
REAL :: F_hp_face
!DUMMY
REAL :: head_at_face, Zface !HEAD and ELEVATION at the face
!
head_at_face = ( (hin+zin)*dxin + (hout+zout)*dxout )  &
                        / ( dxin + dxout)
Zface        = ( zin*dxin + zout*dxout ) / ( dxin + dxout)
F_hp_face = MIN ( b, head_at_face - Zface )
END function 
!===========================================================================
!function to determine the thickness on the surface at the cell face
FUNCTION F_hs_face(dxin, hin, zin, dxout, hout, zout, b)
implicit none
!INPUTS
REAL :: dxin, dxout  ! size of cells
REAL :: hin, hout    ! thickness at CV center
REAL :: zin, zout    ! elevation at CV center
REAL :: b            ! pavement thickness
REAL :: F_hs_face
!DUMMY
REAL :: head_at_face, Zface !HEAD and ELEVATION at the face
!
head_at_face = ( (hin+zin)*dxin + (hout+zout)*dxout )  &
                        / ( dxin + dxout)
Zface        = ( zin*dxin + zout*dxout ) / ( dxin + dxout)
F_hs_face = MAX ( 0., head_at_face - Zface - b )
END FUNCTION 
!=======================================================================

end module pfc1dfuns2

