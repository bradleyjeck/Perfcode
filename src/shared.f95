! fortran_free_source
!
!  This module is part of PERFCODE, written by Bradley J. Eck.
!
!   History:  5 April 2010  -- original coding
!                July 2011  -- revised
!
!   
!
!   Purpose:  This module declares variables to be used globally
!  
!     Notes:  - Variable organization tries to mirror program
!               the organization of the program
!             - See begining of main program for alphabetical
!               listing of variables with descriptions
!             - Use ONLY statement in subroutines to restrict 
!               access to variables in this module

MODULE SHARED  
implicit none
save  

!--------------------------------------------------
!PARAMETERS INPUT FILE
!--------------------------------------------------
!   PFC Properties
logical :: params_vary
REAL, TARGET :: K_input      ! Hydraulic Conductivity [m/s]
REAL, TARGET :: por_input
REAL, TARGET :: b_pfc_input  !PFC Thickness [ m ]
REAL, TARGET :: n_mann_input !Manning's n [ s / m ^(1/3) ]
!   Physical constants
REAL :: g      ! Gravitational Acceleration [m/s/s]
! Time Steps
REAL :: dt_pfc, dt_sheet, max_time 
! Grid Spacing
REAL :: dx, dy
!Tolerances
INTEGER :: qmax, maxit
REAL :: eps_matrix, eps_itr
REAL :: relax, relax_tran
!Initial Condition
real :: h0 ! initial depth in meters
!Boundary Conditions
logical :: BCs_vary
character( len=7 ) :: north_bc_input, south_bc_input
character( len=7 ) :: east_bc_input, west_bc_input
!Animation Options
logical :: animate ! at all and for this step
real :: dt_ani
!Mode of the model
character( len = 2 ) :: model_mode


!--------------------------
! OTHER PARAMETERS
!--------------------------
INTEGER, PARAMETER :: max_rec = 20000
REAL, PARAMETER :: h_pfc_min  = 1.e-10   ! use this instead of TINY

! compass directions as integers by azimuth or bearing
! for use in conveyance subroutine and elsewhere
!  using keywords as integers rather than characters speeds up
!   if/then statements
INTEGER, PARAMETER :: north = 360, south = 180, east = 90, west = 270
INTEGER, PARAMETER :: NorthEast =  45, SouthEast = 135
INTEGER, PARAMETER :: NorthWest = 315, SouthWest = 225
INTEGER, PARAMETER :: old = 0, itr = 1, dry = 2

!---------------------------
! RAINFALL
!---------------------------
INTEGER :: nrr  ! Number of rainfall records
REAL, DIMENSION( max_rec ) :: rain_time, rain_rate


!---------------------------------------------------------
!  GRID GENERATION  
!---------------------------------------------------------

!----------------------
!  Derived data types
type CLSEG  !describes a centerline segment
    real xcc1, ycc1, dx, dy, R1, dR, W, theta1, dtheta, arclen
end type CLSEG

type gridcell   ! Summary information for a grid cell
    integer :: i, j, segment
    real    :: xi, eta
    real    :: X , Y
end type gridcell

! allocatable variables of derived types
type(CLSEG), allocatable, dimension(:) :: seg
type(gridcell) , allocatable, dimension(:) :: CV_Info  !17
!-------------------------

!Array sizes
integer :: imax, jmax, vmax

! Grid numbering scheme
integer, allocatable, dimension(:,:) :: grid 

! Geometric Arrays
REAL, ALLOCATABLE, DIMENSION(:,:) :: lng, wid, area, Z  
REAL, ALLOCATABLE, DIMENSION(:,:) :: lng_south, lng_north


!-------------------------------------------------
!  ELEVATIONS
!------------------------------------------------

!CROSS SECTION  ( Transverse direction)
! input file
integer :: nr_cs
REAL ::  slope_cs(10), wid_cs(10)

! derived values
real, dimension( 11 ) :: eta_cs=0., Z_cs=0.

!LONGITUDINAL PROFILE
integer :: nr_lp
real, dimension(100) :: dist_lp, Z_lp
real :: long_slope  !longitudinal slope at each end of domain


!  1D GRID GENERATION
integer, TARGET :: TNE
REAL, ALLOCATABLE, DIMENSION(:), TARGET :: EDX, XCV, ZCV, etaCV
! 1D model/boudary conditions
real, allocatable, dimension(:) :: h_old_1d, h_new_1d, q_surf, q_pav
real, allocatable, dimension(:) :: slope_cs_1D, wid_cs_1d, eta_cs_1D

!-------------------------------
!  PFC PROPERTIES  & BOUNDARY CONDITIONS
!-----------------------------
real, allocatable, dimension(:, :) :: hyd_cond, por, b_pfc, n_mann
character(len=7), allocatable, dimension(:) :: north_bc, south_bc
character(len=7), allocatable, dimension(:) :: east_bc, west_bc
logical, allocatable, dimension(:,:) :: is_boundary, is_corner

!------------------------------------------------------------
! INTERMEDIATE VARIABLES
!-------------------------------------------------------------

! ARRAY INDICES AND LIMITING VALUES
integer :: i, j, v, q, n
integer :: ve
integer :: v_in !global index of 'inside' adjacent cell
integer :: nmax     ! maximum number of time steps
integer :: nlast    !the last timestep taken

! TIME STUFF
REAL :: dt
REAL, ALLOCATABLE, DIMENSION(:) :: rain  !rainfall depth for each time step
REAL, ALLOCATABLE, DIMENSION(:) :: time  
real :: time_simulated = 0.
character( len = 9 ) :: out_time   ! Characters to for internal writes to store
character( len = 9 ) :: sim_time   ! simulation time w/o floating point error 99999.99
character( len = 9 ) :: sim_time2

! FRICTION SLOPES, POROSITY FUNCTIONS, AND CONVEYANCE COEFFICIENTS
!  'old' means time level 'n'
!  'itr' or '1' means time level n+1
REAL, ALLOCATABLE, DIMENSION(:,:), TARGET :: Sfw_old, Sfe_old, Sfs_old, Sfn_old  
REAL, ALLOCATABLE, DIMENSION(:,:), TARGET :: Sfw_itr, Sfe_itr, Sfs_itr, Sfn_itr  
REAL :: pf, pf1                
REAL :: Cw , Ce , Cs , Cn     
REAL :: Cw1, Ce1, Cs1, Cn1   

! BOUNDARY CONDITION STUFF
real :: eta_1D
real :: hs1, hs2, ds     ! Sheet flow MOC
real :: hp1, hp2, dx_moc ! PFC flow MOC
real :: h_bound          ! depth at boundary (returnd by MOC_KIN or 1D_FLOW
real :: eta_0_hp2_max    ! max possible value for the MOC BC 

! CONVERGENCE TESTING 
logical :: transition
real :: relaxation_factor
REAL :: eps_itr_tol
integer :: pf_int, pf1_int  ! use integers to detect transition
REAL, ALLOCATABLE, DIMENSION(:) :: residual, relchng
real :: maxrelchng_ss 

! LINEAR SYSTEM
!  Bands
REAL, ALLOCATABLE, DIMENSION(:) :: A, B, C, D, E, Fn, F1, F  
!  Test for diagonal Dominance
logical diagdom 
!  Square matrix for outputting/use with library solvers

!---------------------------------------------------------
! THE SOLUTION (at various stages and in various formats)
!---------------------------------------------------------

! Vector Form
REAL, ALLOCATABLE, DIMENSION(:) :: h_itr_vec, h_tmp_vec  
REAL, ALLOCATABLE, DIMENSION(:) :: h_old_vec, h_new_vec

! Vector form, within a timestep (during an iteration)
real, allocatable, dimension(:,:) :: h_temp_hist 

! Vector form, at intervals for animation
!   rows --> grid cells
!   cols --> times
REAL, ALLOCATABLE, DIMENSION(:,:) :: h_vec_ani 

! Matrix Form
REAL, ALLOCATABLE, DIMENSION(:,:), TARGET :: h_old, h_itr, h_new  

! Matrix form, at special times
real, allocatable, dimension(:,:) :: h_max, h_Q_max
real, allocatable, dimension(:,:) :: h_imid_j1_max, h_imid_max

!------------------------------------------------------
! SUMMARY INFORMATION
!------------------------------------------------------

! Input variables and values
character( len=15), dimension(17) :: input_variables
real, dimension(17) :: input_values

! Information about each timestep
INTEGER, ALLOCATABLE, DIMENSION(:) :: numit, loc 
REAL, ALLOCATABLE, DIMENSION(:) :: maxdiff 
real, allocatable, dimension(:) :: maxthk 
integer, allocatable, dimension(:) :: matrix_numits
real, allocatable, dimension(:) :: Qout, L2_History
integer :: solver_numits, timestep_solver_numits

! time history of the depth at i=imax/2 j=1
real, allocatable, dimension( : ) :: h_imid_j1_hist, h_imid_max_hist 

! location of transition to sheet flow
real :: X_sheet
real, allocatable, dimension( : ) :: eta_sheet

! Water budget
real :: old_vol
real, allocatable, dimension(:) :: pfc_vol, surface_vol, rain_vol
real, allocatable, dimension(:) :: north_vol, south_vol, east_vol, west_vol
real, allocatable, dimension(:) :: net_flow, storage_change, vol_error, vol_error_frac

!---------------------------------------------------
!  MISCELLANEOUS  (gotta love this category)
!----------------------------------------------------

integer, dimension(200) :: astat=0     ! for keeping track of allocation statuses

CHARACTER(8)  FILE_DATE
CHARACTER(10) FILE_TIME

character (len=10) :: timechar, datechar, zonechar
integer, dimension(8) :: dtvals

! Routine timing
REAL    :: cputime
integer :: run_start_time, run_end_time, count_rate, count_max
integer :: report = 1  ! determine if we should write out the timestep. 

! For animation output
integer :: ani = 0  ! use this like 'report'
integer :: animax   ! maximum value of ani, compute from max_time / ani_step  
character( len = 10 ), allocatable, dimension(:) :: ani_lab   ! labels for animtaion output
real, allocatable, dimension(:) :: ani_time
character(len = 10) :: lab

END MODULE SHARED
