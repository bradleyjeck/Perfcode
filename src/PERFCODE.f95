! fortran_free_source
!
! (c) Copyright 2010, 2015 Bradley J Eck 
!
!
!
!   PPPP     EEEEEE   RRRR      FFFFFF      CCCC        OOOO     DDDDD     EEEEEE
!   P   P    E        R   R     F          C    C      O    O    D    D    E
!   P    P   E        R    R    F         C           O      O   D     D   E 
!   P    P   E        R    R    F        C           O        O  D      D  E     
!   P   P    E        R   R     F        C           O        O  D      D  E 
!   PPPP     EEEEEE   RRRR      FFFFFF   C           O        O  D      D  EEEEEE 
!   P        E        R   R     F        C           O        O  D      D  E 
!   P        E        R    R    F        C           O        O  D      D  E
!   P        E        R     R   F         C           O      O   D     D   E
!   P        E        R     R   F          C    C      O    O    D    D    E
!   P        EEEEEE   R     R   F           CCCC        OOOO     DDDDD     EEEEEE
!
!
!   P E R m e a b l e   F r i c t i o n   C O u r s e    D r a i n g e   c o d E
!
!       Written By:     Brad Eck
!
!             Date:     April 2010
!
!================================================================================
!       \\\\\\\\\\     P R O G R A M                            //////////
!       //////////                     D E S C R I P T I O N    \\\\\\\\\\
!================================================================================ 
!
!  
!       Purpose:        This program computes a 2D solution for unsteady
!                       drainage through a PFC.  The water THICKNESS in each
!                       cell is used as the primary variable.  
!       IC:             Specified in input file
!       BCs:            Specified in input file
!       Linearization:  Picard Iteration (lag the coefficients)
!       Linear Solver:  Gauss-Seidel iteration
!       
! Alphabetical list of variables used in the main program PERFCODE 
!       (variables used in subroutines are described there)
!
!   A         -- lowest band of penta diagonal matrix 
!   area      -- area of a grid cell
!   astat     -- array allocation statuses
!   B         -- subdiadonal band of penta diagonal matrix
!   b_pfc     -- thickness of the PFC layer
!   C         -- main diagonal of penta diagonal matrix
!   Ce        -- conveyance coefficient ( conv coef ) for the 
!                EASTtern cell face at time level n
!   Ce1       -- conv coef for EASTern cell face at time level n + 1
!   Cn        -- conv coef for the NORTHern cell face at time level n
!   Cn1       --      ''             ''           ''  at time level n + 1
!   Cs        -- conv coef for the SOUTHern cell face at time level n
!   Cs1       --      ''             ''           ''  at time level n + 1
!   CV_Info   -- information about each grid cell (aka Control Volume)
!   Cw        -- conv coef for the WESTern  cell face at time level n
!   Cw1       --      ''             ''           ''  at time level n + 1
!   D         -- superdiagonal band of penta diagonal matrix
!   dist_lp   -- distance along longitudinal profile
!   diagdom   -- logical flag for test of diagonal dominance
!   ds        -- distance up characteristic in sheet flow moc bc
!   dt        -- time step for the simulation
!   dt_pfc    -- time step for PFC flow
!   dt_sheet  -- time step for sheet flow
!   dx        -- prelim. grid size for longitudinal direction
!   dx_moc    -- distance up drainage path in pfc moc bc
!   dy        -- prelim. grid size for transverse direction
!   E         -- uppermost band of penta diagonal matrix
!   east_bc   -- condition for east boundary
!   eps_matrix-- tolerance (epsilon) for matrix solver
!   eps_itr   -- tolerance for an iteration
!   eps_itr_tol-- selected tolerance for the iteration (based on transition)
!   eps_ss    -- tolerance for steady state (not used)
!   eta_cs    -- values of eta along the cross slope
!   eta_0_hp2_max-- max possible value for pfc moc bc
!   eta_cs_1D -- values of eta for 1D model
!   eta1D     -- value of eta at each point in 1D domain
!   etaCV     -- value of eta at CV center for 1D grid
!   F_        -- the letter F with an underscore ( F_ ) denotes a 
!                    function call and NOT an array
!   F         -- right hand side of linear system in pentadiagonal matrix
!   F1        -- contribution to F from time level n+1
!   Fn        -- contribution to F from time level n
!   g         -- constant of gravitational acceleration  
!   grid      -- number of each grid cell
!   h0        -- initial depth (m)
!   h_bound   -- depth at boundary (returned by MOC_KIN or 1D_FLOW)
!   h_imid_j1_max-- solution when depth at middle of south boundary is max
!   h_imid_j1_max_hist
!   h_imid_max-- solution when depth in middle of domain is max
!   h_imid_max_hist
!   h_itr     -- matrix form of solution at level n+1
!   h_itr_vec -- vector form of solution at time level n+1 
!   h_max     -- solution at maximum depth
!   h_new_1d  -- solution at time level n+1 for 1D problem 
!   h_old     -- solution at time level n
!   h_old_1d  -- solution at time level n for 1D problem
!   h_old_vec --
!   h_pfc_min -- minimum value for pfc flow thickness
!   h_Q_max   -- solution at maximum flow
!   h_temp_hist -- history of solution during an iteration
!   h_tmp_vec --
!   hp1       -- depth at point 1 in pfc MOC bc
!   hp2       -- depth at point 2 in pfc MOC bc
!   hs1       -- sheet flow depth at point 1 in sheet flow moc bc
!   hs2       -- sheet flow depth at point 2 in sheet flow moc bc
!   i         -- array index ( longitudinally in the domain )
!   input_values-- array of values of the input variables
!   input_variables -- character array of input variables
!   imax      -- maximum value of the array index i
!   j         -- array index ( transverse in the domain )
!   jmax      -- maximum value of the array index j
!   K         -- the saturated hydraulic conductivity of the PFC  
!   L2_history -- value of the L2 norm for each timestep 
!   lng       -- curvilinear length of a grid cell at its center
!   lng_north -- curvilinear length of the northern face  
!   lng_south -- curvilinear length of the southern face
!   loc       -- the location of the largest relative change in a time step
!   long_slope -- overall longitudinal slope
!   max_rec   -- maximum number of records (for pre-allocating arrays
!                where values are read in from a file )
!   max_time  -- longest time to simulate
!   maxdiff   -- the change in head at location LOC for timestep n
!   maxit     -- maximum number of matrix iterations 
!   maxrelchng_ss-- maximum relative change for a timestep, for stdy state check
!   maxthk    -- maximum thickness fot the timestep
!   matrix_numits-- number of iterations to solve the matrix
!   n         -- index for time stepping
!   n_mann    -- Manning's roughness coefficient
!   north_bc  -- condition for north boundary
!   nlast     -- last timestep taken
!   nmax      -- maximum number of time steps in the simulation
!   numit     -- the number of iterations required for a timestep to converge
!   nr_cs     -- number of records in the cross slope file
!   nr_lp     -- number of records in the longitudinal profile file
!   nrr       -- number of rainfall records
!   out_time  -- 
!   pf        -- porosity factor ( includes effect of porosity 
!                                   when pavement is not saturated )
!   pf_int    -- porosity factor as an integer
!   pf1       -- porosity factor for time level n+1
!   pf1_int   --  "      "                      ""   as integer
!   por       -- the effective porosity of the PFC
!   q         -- iteration index
!   qmax      -- maximum number of iterations
!   rain      -- rainfall rate for each timestep of the simulation
!   Qout      -- flow rate out the southern boundary for a timestep
!   rain_rate -- rainfall rate for each time increment in the 
!                   rainfall input file
!   rain_time -- time column of rainfall input file
!   relax     -- relaxation factor for non-transition iterations
!   relaxation_factor -- underrelaxation factor for non-linear iteration
!   relax_tran -- relaxation factor for transition
!   relchng   -- the relative change between solutions for an iteration or timestep
!   residual  -- difference between old and itr solutions 
!   seg       -- properties of a centerline segment
!   Sfe_itr   -- friction slope at center of east  face at time level n+1  
!   Sfe_old   -- friction slope at center of east  face at time level n
!   Sfn_itr   -- friction slope at center of north face at time level n+1  
!   Sfn_old   -- friction slope at center of north face at time level n
!   Sfs_itr   -- friction slope at center of south face at time level n+1  
!   Sfs_old   -- friction slope at center of south face at time level n
!   Sfw_itr   -- friction slope at center of west  face at time level n+1  
!   Sfw_old   -- friction slope at center of west  face at time level n
!   slope_cs  -- slope column of cross section file
!   slope_cs_1d -- slope of 1D segment 
!   sim_tim   -- character variable for time simulated
!   solver_numits-- number of iterations for the solver
!   south_bc  -- condition for south boundary
!   time      -- time at each timestep
!   time_simulated-- the time simulated
!   timestep_solver_numits -- 
!   transition -- logical to see if we're in a transition timestep
!   tolit     -- tolerence for iterations, used for relative (fractional) changes
!   TNE       -- total number of elements for 1D grid
!   v         -- linear index for domain
!   ve        -- linear index for cell to the east
!   v_in      -- linear index of adjacent inside cell
!   vmax      -- number of unknowns in the domain   
!   west_bc   -- condition for west boundary
!   wid       -- curvilinear width of a grid cell at its center
!   wid_cs    -- width column of cross slope file
!   wid_cs_1d -- width of 1D segment
!   XCV       -- coordinate of CV center for 1D grid
!   Z         -- elevation at the cell center
!   Z_cs      -- elevation along the cross slope
!   Z_lp      -- elevation along longitudinal profile
!   ZCV       -- elevation of CV center for 1D grid

!===================================================================================
!       \\\\\\\\\\      B E G I N    P R O G R A M         //////////
!       //////////           P E R F C O D E               \\\\\\\\\\
!===================================================================================
program PERFCODE

!-------------------------------------------------------------------------
!   >>>>>>>>>>          M O D U L E S                 <<<<<<<<<<
!-------------------------------------------------------------------------
! Refer to the modules that are referred to by this code

USE SHARED      ! SHARED is used to store VARIABLES
USE INPUTS      ! INPUTS has subroutines
USE OUTPUTS     ! OUTPUTS has subroutines  
USE ConvCoef    ! computes conveyance coefficinnts
USE SOLVERS     ! linear solvers
USE Utilities
USE gridgen
use pfc1Dsubs
use pfc2Dsubs
use pfc2Dfuns
use BoundCond


!-------------------------------------------------------------------------
!   >>>>>>>>>>          V A R I A B L E S                  <<<<<<<<<<
!-------------------------------------------------------------------------

implicit none   

!  All variables are declared in module SHARED
      
!-------------------------------------------------------------------------------------
!   >>>>>>>>>>      P R O B L E M    S E T U P          <<<<<<<<<<
!-------------------------------------------------------------------------------------

!Create a file to store details of the run
open( unit = 100, file = 'PERFCODE_Run.txt', status = 'REPLACE' )

!----------------------------------
! Problem parameters file
!-----------------------------------
CALL GET_PARAMETERS( K, por, b_pfc, n_mann, g, dt_pfc, dt_sheet, max_time,  &
                     dx, dy, qmax, maxit, h0, eps_matrix, eps_itr, eps_ss,  &
                     relax, relax_tran,                                     & 
                     north_bc, south_bc, east_bc, west_bc,                  &
                     animate, dt_ani                                          )


!----------------------------------------------------
! Rainfall file & maximum number of timesteps
!----------------------------------------------------

call GET_RAINFALL( max_rec, rain_time, rain_rate, nrr )

nmax = ( maxval( rain_time(1:nrr) ) / min( dt_pfc, dt_sheet ) ) !* 100


!--------------------------------------------------------
! GRID GENERATION
!-------------------------------------------------------
!   This subroutine takes the centerline geometry file that is generated
!   mannualy and creates a curvilinear grid. 
!       INPUTS:  Preliminary grid spacing
!       OUTPUTS: Size of computational domain (imax & jmax )
!                Length, width and area of each grid cell ( module SHARED)
!                Coordinates of each CV center 
call GENERATE_GRID( prelim_dx = dx , prelim_dy = dy )


! Reads in cross section and longitudinal files and computes elevations
!  of CV Centers

CALL SET_ELEVATIONS()


! Creates a grid for a 1D section in case a 1D boundary condition is used
call setup_1d_section()
CALL grid_1d_section( slope_in = slope_cs_1D  ,  &
                      width_in =   wid_cs_1D  ,  &
                           seg =       nr_cs  ,  &
                            dx = (( dx+dy ) / 2.)   )  

!-----------------------------------------------------------
! inputs summary
!-----------------------------------------------------------
! make a list of input variables and values
input_variables = (/ 'K         ', &
                     'por       ', &
                     'b_pfc     ', &
                     'n_mann    ', &
                     'g         ', &
                     'dt_pfc    ', &
                     'dt_sheet  ', &
                     'max_time  ', &
                     'dx        ', &
                     'dy        ', &
                     'qmax      ', &
                     'maxit     ', &
                     'h0        ', &
                     'eps_matrix', &
                     'eps_itr   ', &
                     'eps_ss    ', &
                     'relax     ', &
                     'relax_tran' /) 

!also collect and store values of input variales
input_values =  (/ K, por, b_pfc, n_mann, g, dt_pfc, dt_sheet, &
                  max_time, dx, dy, real(qmax), real(maxit),   &
                  h0, eps_matrix, eps_itr, eps_ss, relax, relax_tran /)


! Echo inputs to the screen, unit 6 by default
CALL ECHO_INPUTS( dev = 6 ) 
!also echo to log file
CALL ECHO_INPUTS( dev = 100 )


!--------------------------------------
!  Animation setup
!---------------------------------------

if( animate .eqv. .TRUE. ) then
    
    animax = int( floor(max_time / dt_ani) ) 
    allocate( h_vec_ani ( vmax, animax  ) )
    allocate(   ani_lab (       animax  ) )
    allocate(   ani_time(       animax  ) )

endif

!-------------------------------------------------------------------------------
!   >>>>>>>>>>      A L L O C A T E   A R R A Y S        <<<<<<<<<<
!-------------------------------------------------------------------------------

! inialize as we go.

! VARIABLES IN MODULE SHARED

allocate(     h_old( imax, jmax ),  STAT = astat( 7) );       h_old = 0.0 
allocate(     h_itr( imax, jmax ),  STAT = astat( 8) );       h_itr = 0.0
allocate(   Sfw_old( imax, jmax ),  STAT = astat( 9) );     Sfw_old = 0.0
allocate(   Sfe_old( imax, jmax ),  STAT = astat(10) );     Sfe_old = 0.0
allocate(   Sfs_old( imax, jmax ),  STAT = astat(11) );     Sfs_old = 0.0
allocate(   Sfn_old( imax, jmax ),  STAT = astat(12) );     Sfn_old = 0.0
allocate(   Sfw_itr( imax, jmax ),  STAT = astat(13) );     Sfw_itr = 0.0
allocate(   Sfe_itr( imax, jmax ),  STAT = astat(14) );     Sfe_itr = 0.0
allocate(   Sfs_itr( imax, jmax ),  STAT = astat(15) );     Sfs_itr = 0.0
allocate(   Sfn_itr( imax, jmax ),  STAT = astat(16) );     Sfn_itr = 0.0
allocate(     h_max( imax, jmax )                    );       h_max = 0.0
allocate(   h_Q_max( imax, jmax )                    );     h_Q_max = 0.0
allocate( h_imid_j1_max ( imax, jmax )              );h_imid_j1_max = 0.0
allocate( h_imid_max_hist( nmax ) );                h_imid_max_hist = 0.0
allocate( h_imid_max( imax, jmax ) ) ;                   h_imid_max = 0.0

allocate( h_old_1d( TNE ) )
allocate( h_new_1d( TNE ) )


! Check allocation statuses
do i = 1, 20
    if( astat( i ) .NE. 0 ) then
        WRITE(100,*) 'PERFCODE: allocation problem!! &
                            & check shared variable:', i
    end if
end do

if( maxval(astat) .eq. 0 ) then
    WRITE(100,*) 'PERFCODE: allocation of shared variables sucessful'
endif


! VARIABLES IN THIS PROGRAM

allocate(   A( vmax ), stat = astat2( 1) );     A = 0.0
allocate(   B( vmax ), stat = astat2( 2) );     B = 0.0
allocate(   C( vmax ), stat = astat2( 3) );     C = 0.0
allocate(   D( vmax ), stat = astat2( 4) );     D = 0.0
allocate(   E( vmax ), stat = astat2( 5) );     E = 0.0
allocate(  Fn( vmax ), stat = astat2( 6) );     Fn= 0.0
allocate(  F1( vmax ), stat = astat2( 7) );     F1= 0.0
allocate(   F( vmax ), stat = astat2( 8) );     F = 0.0


allocate(   h_itr_vec( vmax ),  stat = astat2( 9) );    h_itr_vec = 0.0
allocate(   h_tmp_vec( vmax ),  stat = astat2(10) );    h_tmp_vec = 0.0
allocate(   h_old_vec( vmax ),  stat = astat2(11) );    h_old_vec = 0.0
allocate(   h_new_vec( vmax ),  stat = astat2(12) );    h_new_vec = 0.0
allocate(   relchng  ( vmax ),  stat = astat2(13) );    relchng   = 0.0

allocate(     numit( nmax ),    stat = astat2(14) );      numit   = 0
allocate(       loc( nmax ),    stat = astat2(15) );      loc     = 0
allocate(   maxdiff( nmax ),    stat = astat2(16) );      maxdiff = 0.0
allocate(      Qout( nmax )  ); Qout = 0.0
allocate(   matrix_numits(nmax) )
allocate(  L2_History( nmax ) ); L2_History = 0.0


! Set indices for rain so that n-1 always works.  This is b/c
!  in Crank-Nicolson half of the rainfall rate is from time level
!  n  and half is from time level n-1  
allocate(  rain( 0 : nmax-1 ),    stat = astat2( 0) )   
allocate(  time( nmax ),    stat = astat2(17) )


allocate(   grid( jmax, imax ),     stat = astat2(19) )

allocate( maxthk( nmax ), stat = astat2(20) );             maxthk = 0.0

allocate( residual( vmax ), stat = astat2(21) );         residual = 0.0



allocate( h_temp_hist( vmax, qmax) ); h_temp_hist = 0.0 

allocate( h_imid_j1_hist( nmax ), stat = astat(22) )
! Check allocation statuses
do i = 1, 29
    if( astat2( i ) .NE. 0 ) then
        WRITE(100,*) 'PERFCODE: allocation problem in main &
                                & program!!, check variable:', i
    end if
end do

if( maxval(astat2) .eq. 0 ) then
    WRITE(100,*) 'PERFCODE: allocation of main program variables sucessful'
endif


!-------------------------------------------------------------------------------
!   >>>>>>>>>>      P R O B L E M    S O L V I N G        <<<<<<<<<<
!-------------------------------------------------------------------------------

!  INITIAL CONDITIONS
!  set all all arrays to the initial depth value
h_old = h0     
h_itr = h0   ! added this after b/c the first iteration kept failing
h_old_vec = h0
h_itr_vec = h0

h_old_1D = h0   ! initial depth for 1D boundary condition
h_new_1D = h0


WRITE(*,*) 'PERFCODE: starting time stepping loop,&
            & max time = ', max_time, ' seconds'

CALL SYSTEM_CLOCK( RUN_START_TIME, count_rate, count_max)


 
!      !open a file to store each timestep
!          open( unit = 50, file = 'timesteps.csv', status = 'REPLACE' )
!          write(50,5) ' n / v,', (v, v=1, vmax) !implied DO loop
!      5   format( A, 10000( I, ','))
!


! Set rainfall rate for begining of simulation
n=0
rain(n) = F_Linterp( 0.0             , &
                     rain_time(1:nrr), &
                     rain_rate(1:nrr), &
                     nrr                  )

!----------------------------------------------
! BEGIN TIME STEPPING
!----------------------------------------------

time_stepping: do while (time_simulated .LT. max_time )

!increment n and store the largest n we've gotten so far
n = n + 1
nlast = n

! Select the time step 
if( maxval( h_old ) .GT. b_pfc * 0.95 ) then
    dt = dt_sheet
else
    dt = dt_pfc
endif

!Computed the time simulated  
!   Do the accumulation with an internal write/read to 
!    avoid accumulating the floating point errors

write( sim_time, 123 ) time_simulated
read( sim_time,  * ) time_simulated

123 format( F8.2 )

time_simulated = time_simulated + dt
time(n) = time_simulated
!Report which timestep we're in every 20 or so time steps
if( nint( real(n)/2. ) .gt. report ) then
    report = report + 1
    write(*,*) ' n = ', n, ' time = ', time_simulated,          &
                      'L_inf_norm = ', maxrelchng_ss,           &
                         'L2_norm = ', F_L2_Norm(relchng,vmax), &
                            'Qout = ', Qout(n-1)
endif


!Come up with the rainfall rate for this timestep
rain(n) = F_Linterp( time_simulated  , &
                     rain_time(1:nrr), &
                     rain_rate(1:nrr), &
                     nrr                  )


!PART OF NON-LINEAR SYSTEM FROM TIME LEVEL n
!   FRICTION SLOPE
!   Compute friction slope magnitudes based on the converged thicknesses
!   from the previous time step
CALL FrictionSlope( 'old', Sfw_old, Sfe_old, Sfs_old, Sfn_old )



       ! Compute solution for 1D model to use as a boundary condition

! only invoke the 1D solver if called for by the boundary conditions

if( west_bc .eq. '1D_FLOW' .or. &
    east_bc .eq. '1D_FLOW'        ) then

       h_old_1d = h_new_1d

       CALL PFC1DIMP(  h_old = h_old_1d, &
                          dt = dt      , &
                        rain = rain(n) , &   !  Should probably add rain(n-1)
                       tolit = eps_itr , &
                        qmax = qmax    , & 
                       h_new = h_new_1d, &
                        imax = TNE     , &
                    eta_0_BC = south_bc, &
                    eta_1_BC = north_bc   )


       ! Vet the solution to avoid a weird problem
       if( maxval( h_new_1d) .LT. TINY(h_new_1d(1) ) ) then
           write(100,*) 'PERFCODE: 1D Model zeroed out....stopping program'
           call write_vector( h_old_1d, TNE, 'h_old_1D.csv' )
           call write_vector( h_new_1d, TNE, 'h_new_1D.csv' )
           stop
       end if

endif

!-------------------------------------------
!   B O U N D A R Y   C O N D I T I O N S
!-------------------------------------------

! put east first so that west bc 'eastKIN' could copy it

!EASTERN BOUNDARY
i = imax
if( east_bc .eq. 'NO_FLOW' ) then
        do j = 2, jmax - 1
           pf = F_por( h_old( i, j ) )       
           CALL Conveyance( 'west ', 'old', i, j, Cw )
           Ce = 0.0                        !<---- NO FLOW BOUNDARY
           CALL Conveyance( 'south', 'old', i, j, Cs )
           CALL Conveyance( 'north', 'old', i, j, Cn )
           v     = F_LinearIndex( i, j, jmax)
           Fn(v) = F_RHS_n( i, j, Cw, Ce, Cs, Cn, rain(n-1), pf, dt )
        end do

elseif( east_bc .eq. '1D_FLOW') then

!open( unit = 66, file = 'eta_mapping.csv', status = 'REPLACE' )
!write( 66, * ) 'i,j,eta,eta_1D'

        do j = 1, jmax
            v = F_LinearIndex( i, j, jmax)
           eta_1D = F_LINTERP(        X = CV_Info( v ) % eta ,  &
                                known_X = eta_cs             ,  &
                                known_Y = eta_cs_1D          ,  &
                                      n = nr_cs + 1                  ) 
           h_bound= F_LINTERP(        X = eta_1D             ,  &
                                 known_X = etaCV              ,  &
                                 known_Y = h_new_1D           ,  &
                                       n = TNE                       )
            C(v) = 1.0
            F(v) = h_bound

!            write( 66, 660) i, j, CV_Info( v ) % eta, eta_1D
            
        end do 

!close(66)


elseif( east_bc .eq. 'MOC_KIN' ) then

        do j = 2, jmax - 1
            CALL MOC_KIN_BC( i, j, rain(n), dt, 'east ', h_bound, 100 )
            v = F_LinearIndex( i, j, jmax )
            C(v) = 1.
            F(v) = h_bound
!            write(100,*) 'PERFCODE: east bc i=',i, 'j=',j, 'h_bound=',h_bound
        end do   


end if



!WESTERN BOUNDARY 
i = 1
if( west_bc .eq. 'NO_FLOW' ) then
        do j = 2, jmax - 1
           ! Set porosity factor for this cell
           pf = F_por( h_old( i, j ) )       
           ! Set the conveyance coefficients 
           Cw = 0.0     !<---- NO FLOW BOUNDARY
           CALL Conveyance( 'east ', 'old', i, j, Ce )
           CALL Conveyance( 'south', 'old', i, j, Cs )
           CALL Conveyance( 'north', 'old', i, j, Cn )
           ! Compute the part if the right-hand-side that is from 
           ! time level n
           v     = F_LinearIndex( i, j, jmax)
           Fn(v) = F_RHS_n( i, j, Cw, Ce, Cs, Cn, rain(n-1), pf, dt )
        end do
elseif( west_bc .eq. '1D_FLOW' ) then
        do j = 1, jmax
            v = F_LinearIndex( i, j, jmax)
            eta_1D = F_LINTERP(        X = CV_Info( v ) % eta ,  &
                                 known_X = eta_cs             ,  &
                                 known_Y = eta_cs_1D          ,  &
                                       n = nr_cs + 1                ) 
            h_bound= F_LINTERP(        X = eta_1D             ,  &
                                 known_X = etaCV              ,  &
                                 known_Y = h_new_1D           ,  &
                                       n = TNE                      )
            C(v) = 1.0
            F(v) = h_bound
        end do 

elseif( west_bc .eq. 'MOC_KIN' ) then
        write(*,*) 'PERFCODE: Boundary condition ', west_bc, &
                             'not supported for western boundary'

elseif( west_bc .eq. 'eastKIN' ) then

        do j = 2, jmax-1
            v = F_LinearIndex( i, j, jmax )
            ! index of corresponding eastern cell
            ve = F_LinearIndex( imax, j, jmax )  
            ! Use solutions from east side on the west side
            C(v) = C(ve)
            F(v) = F(ve)
        end do

endif





!NORTHERN BOUNDARY
j = jmax
if( north_bc .eq. 'NO_FLOW' ) then
        do i = 2, imax - 1
           ! Set porosity factor for this cell
           pf = F_por( h_old( i, j ) )       
           ! Set the conveyance coefficients 
           CALL Conveyance( 'west ', 'old', i, j, Cw )
           CALL Conveyance( 'east ', 'old', i, j, Ce )
           CALL Conveyance( 'south', 'old', i, j, Cs )
           Cn = 0.0  !  <---- NO FLOW BOUNDARY 
           ! Compute the part if the right-hand-side that is from 
           ! time level n
           v     = F_LinearIndex( i, j, jmax)
           Fn(v) = F_RHS_n( i, j, Cw, Ce, Cs, Cn, rain(n-1), pf, dt )
        end do

elseif( north_bc .eq. 'MOC_KIN' ) then 
        do i = 2, imax - 1
            CALL MOC_KIN_BC( i, j, rain(n), dt, 'north', h_bound, 100)
            v = F_LinearIndex( i, j, jmax )
            C(v) = 1.
            F(v) = h_bound
        end do
!       ! Use the value of the next inside cell for cells
!       ! second from the end of the domain
!       i = 2
!       v = F_LinearIndex( i, j, jmax )
!       v_in = F_LinearIndex( i+1, j, jmax)
!       C(v) = 1.
!       F(v) = F( v_in)
!       i = imax - 1
!       v = F_LinearIndex( i, j, jmax )
!       v_in = F_LinearIndex( i-1, j, jmax)
!       C(v) = 1.
!       F(v) = F( v_in)
!
elseif( north_bc .eq. '1D_FLOW' ) then
        write(*,*) 'PERFCODE: Boundary condition ', north_bc, &
                             'not supported for northern boundary'

elseif( north_bc .eq. 'west_1D' .and. &
        west_bc .eq. '1D_FLOW'       ) then

      ! Put the answer for the northern most cell on the west end (i=1, j=jmax)
      !  in all of the northern cells

      do i = 2, imax - 1
             v = F_LinearIndex( i, j, jmax )
          v_in = F_LinearIndex( 1, jmax, jmax )
          C(v) = 1.
          F(v) = F(v_in)
      end do

end if


!SOUTHERN BOUNDARY
j = 1
if( south_bc .eq. 'NO_FLOW' ) then
        do i = 2, imax - 1
           pf = F_por( h_old( i, j ) )
           CALL Conveyance( 'west ', 'old', i, j, Cw )
           CALL Conveyance( 'east ', 'old', i, j, Ce )
           Cs = 0.0  !<------- NO FLOW BOUNDARY
           CALL Conveyance( 'north', 'old', i, j, Cn )
           v     = F_LinearIndex( i, j, jmax)
           Fn(v) = F_RHS_n( i, j, Cw, Ce, Cs, Cn, rain(n-1), pf, dt )
        end do

elseif( south_bc .eq. 'MOC_KIN' ) then      
        do i = 2, imax - 1
            CALL MOC_KIN_BC( i, j, rain(n), dt, 'south',  h_bound, 100)
            v = F_LinearIndex( i, j, jmax )
            C(v) = 1.
            F(v) = h_bound
        end do
!       ! Use the value of the next inside cell for cells
!       ! second from the west end of the domain
!       i = 2
!       v = F_LinearIndex( i, j, jmax )
!       v_in = F_LinearIndex( i+1, j, jmax)
!       C(v) = 1.
!       F(v) = F( v_in)
!       ! second from east end of domain
!       i = imax - 1
!       v = F_LinearIndex( i, j, jmax )
!       v_in = F_LinearIndex( i-1, j, jmax)
!       C(v) = 1.
!       F(v) = F( v_in )
       
elseif( south_bc .eq. '1D_FLOW' ) then
        write(*,*) 'PERFCODE: Boundary condition ', south_bc, &
                             'not supported for southern boundary'

elseif( south_bc .eq. 'west_1D'  .AND.  &
         west_bc .eq. '1D_FLOW'            ) then

        ! Put the answer for the southern most cell on the west end (v=1)
        !  in all of the southen cells

        do i = 2, imax - 1
            v = F_LinearIndex( i, j, jmax )
           v_in= F_LinearIndex( 1, 1, jmax )
            C(v) = 1.
            F(v) = F(v_in)
        end do
            
end if
               


!----------------------------
! C O R N E R  P O I N T S 
!----------------------------
! only the 1D_FLOW condition is already handled for the corner points

! NORTH EAST CORNER
i = imax; j = jmax
if( north_bc .eq. 'NO_FLOW' .AND. east_bc .eq. 'NO_FLOW' ) then
           ! Set porosity factor for this cell
           pf = F_por( h_old( i, j ) )       
           ! Set the conveyance coefficients 
           CALL Conveyance( 'west ', 'old', i, j, Cw )
           Ce = 0.0  !  <---- NO FLOW BOUNDARY
           CALL Conveyance( 'south', 'old', i, j, Cs )
           Cn = 0.0  !  <---- NO FLOW BOUNDARY 
           ! Compute the part of the right-hand-side that is from 
           ! time level n
           v     = F_LinearIndex( i, j, jmax)
           Fn(v) = F_RHS_n( i, j, Cw, Ce, Cs, Cn, rain(n-1), pf, dt )

elseif( north_bc .eq. 'MOC_KIN' .AND. east_bc .eq. 'NO_FLOW') then
           ! use the depth in the adjacent MOC_KIN cell
           v    = F_LinearIndex( i, j, jmax )
           v_in = F_LinearIndex( i-1, j, jmax )
           C(v) = 1.0
           F(v) = F( v_in )   

elseif( north_bc .eq. 'NO_FLOW' .AND.  &
         east_bc .eq. 'MOC_KIN'          )  then

        ! is a problem when there are no grade breaks
        ! just value of adjacent no flow cell ??
        v    = F_LinearIndex( i, j, jmax )
        A(v) = -1.
        C(v) = 1.
        F(v) = 0.

elseif( Z( imax, jmax) .GE. Z(imax, jmax-1)   .AND. &
              north_bc .NE. 'NO_FLOW'                 ) then

              write( 100, *) ' North east corner drains to the south &
                             &consider NO_FLOW boundary for the north &
                                                  &side of the domain. ' 

elseif( Z( imax, jmax) .LT. Z(imax, jmax-1) .AND.  &
               east_bc .eq. 'MOC_KIN'                 ) then
       
        ! drainage is to the north and MOC KIN will work
        call MOC_KIN_BC( i, j, rain(n), dt, 'east ', h_bound, 100 )
        v = F_LinearIndex( i, j, jmax )
        C(v) = 1.
        F(v) = h_bound

end if



!   NORTH WEST CORNER POINTS
i = 1; j = jmax
if( north_bc .eq. 'NO_FLOW' .AND. west_bc .eq. 'NO_FLOW' ) then
           ! Set porosity factor for this cell
           pf = F_por( h_old( i, j ) )       
           ! Set the conveyance coefficients 
           Cw = 0.0  !  <---- NO FLOW BOUNDARY 
           CALL Conveyance( 'east ', 'old', i, j, Ce )
           CALL Conveyance( 'south', 'old', i, j, Cs )
           Cn = 0.0  !  <---- NO FLOW BOUNDARY 
           ! Compute the part if the right-hand-side that is from 
           ! time level n
           v     = F_LinearIndex( i, j, jmax)
           Fn(v) = F_RHS_n( i, j, Cw, Ce, Cs, Cn, rain(n-1), pf, dt )

elseif( north_bc .eq. 'MOC_KIN' .AND. west_bc .eq. 'NO_FLOW' ) then
           !use the depth from the adjaent MOC_KIN cell
           v    = F_LinearIndex( i, j, jmax )
           v_in = F_LinearIndex( i+1, j, jmax )
           C(v) = 1.0
           F(v) = F( v_in )   

elseif( west_bc .eq. 'eastKIN' ) then

            v = F_LinearIndex( i, j, jmax )
            ! index of corresponding eastern cell
            ve = F_LinearIndex( imax, j, jmax )  
            ! Use solutions from east side on the west side
            C(v) = C(ve)
            F(v) = F(ve)

end if



! SOUTH EAST CORNER
i = imax; j = 1
if( south_bc .eq. 'NO_FLOW' .and. east_bc .eq. 'NO_FLOW' ) then
           pf = F_por( h_old( i, j ) )       
           ! Set the conveyance coefficients 
           CALL Conveyance( 'west ', 'old', i, j, Cw )
           Ce = 0.0  !  <---- NO FLOW BOUNDARY
           Cs = 0.0  !  <---- NO FLOW BOUNDARY
           CALL Conveyance( 'north', 'old', i, j, Cn )
           ! Compute the part of the right-hand-side that is from 
           ! time level n
           v     = F_LinearIndex( i, j, jmax)
           Fn(v) = F_RHS_n( i, j, Cw, Ce, Cs, Cn, rain(n-1), pf, dt )   

elseif( south_bc .eq. 'MOC_KIN' .AND. east_bc .eq. 'NO_FLOW' ) then
           ! use the depth in the adjacent MOC_KIN cell
           v    = F_LinearIndex( i, j, jmax )
           v_in = F_LinearIndex( i-1, j, jmax )
           C(v) = 1.0
           F(v) = F( v_in ) 

elseif( south_bc .eq. 'MOC_KIN' .AND. east_bc .eq. 'MOC_KIN' ) then

        call MOC_KIN_BC( i, j, rain(n), dt, 'east ', h_bound, 100 )
        v = F_LinearIndex( i, j, jmax )
        C(v) = 1.
        F(v) = h_bound
     
end if   
            
! SOUTHWEST CORNER
i = 1; j = 1
if( south_bc .eq. 'NO_FLOW' .AND. west_bc .eq. 'NO_FLOW' ) then
           pf = F_por( h_old( i, j ) )       
           ! Set the conveyance coefficients 
           Cw = 0.0  !  <---- NO FLOW BOUNDARY
           CALL Conveyance( 'east ', 'old', i, j, Ce )
           Cs = 0.0  !  <---- NO FLOW BOUNDARY
           CALL Conveyance( 'north', 'old', i, j, Cn )
           ! Compute the part of the right-hand-side that is from 
           ! time level n
           v     = F_LinearIndex( i, j, jmax)
           Fn(v) = F_RHS_n( i, j, Cw, Ce, Cs, Cn, rain(n-1), pf, dt )

elseif( south_bc .eq. 'MOC_KIN' .AND. west_bc .eq. 'NO_FLOW' ) then
           ! use the depth in the adjacent MOC_KIN cell
           v    = F_LinearIndex( i, j, jmax )
           v_in = F_LinearIndex( i+1, j, jmax )
           C(v) = 1.0
           F(v) = F( v_in )

elseif( west_bc .eq. 'eastKIN' ) then

            v = F_LinearIndex( i, j, jmax )
            ! index of corresponding eastern cell
            ve = F_LinearIndex( imax, j, jmax )  
            ! Use solutions from east side on the west side
            C(v) = C(ve)
            F(v) = F(ve)
    
 
end if


!-------------------------------------------
!  D O M A I N   I N T E R I O R  
!-------------------------------------------
!   Compute the part of the right hand side of the linear system
!   that is from time level n (the stationary part that does not
!   change as the iteration progresses)
do j = 2, jmax -1;  do i = 2, imax - 1

    ! Set porosity factor for this cell
    pf  = F_por( h_old( i, j ) )
    ! Set the conveyance coefficients
    CALL Conveyance( 'west ', 'old', i, j, Cw )
    CALL Conveyance( 'east ', 'old', i, j, Ce )
    CALL Conveyance( 'south', 'old', i, j, Cs )
    CALL Conveyance( 'north', 'old', i, j, Cn )
    ! Compute the part of the right-hand-side that is from 
    ! time level n
    v     = F_LinearIndex( i, j, jmax)
    Fn(v) = F_RHS_n( i, j, Cw, Ce, Cs, Cn, rain(n-1), pf, dt )

end do; end do


!----------------------------------------------
!ITERATIVE (LAGGED) PART OF NON-LINEAR SYSTEM 
!----------------------------------------------

!zero out matrix iteration counter
timestep_solver_numits = 0

iteration: do q = 1, qmax

!   FRICTION SLOPE
!   compute friction slope magnitudes based on the thickness
!   from the previous iteration
CALL FrictionSlope( 'itr', Sfw_itr, Sfe_itr, Sfs_itr, Sfn_itr )


!   BOUNDARY CELLS
!WESTERN BOUNDARY 
if( west_bc .eq. 'NO_FLOW' ) then
        i = 1
        do j = 2, jmax - 1
            pf = F_por( h_itr( i, j ) )       
            Cw1 = 0.0             !<---------------No flow boundary
            CALL Conveyance( 'east ', 'itr', i, j, Ce1 )
            CALL Conveyance( 'south', 'itr', i, j, Cs1 )
            CALL Conveyance( 'north', 'itr', i, j, Cn1 )
            CALL set_ABCDEF( i, j, Cw1, Ce1, Cs1, Cn1, pf, dt, rain(n) )
        end do
end if

!EASTERN BOUNDARY
if( east_bc .eq. 'NO_FLOW') then
        i = imax
        do j = 2, jmax - 1
            pf = F_por( h_itr( i, j ) )       
            CALL Conveyance( 'west ', 'itr', i, j, Cw1 )
            Ce1 = 0.0             !<---------------No flow boundary
            CALL Conveyance( 'south', 'itr', i, j, Cs1 )
            CALL Conveyance( 'north', 'itr', i, j, Cn1 )
            CALL set_ABCDEF( i, j, Cw1, Ce1, Cs1, Cn1, pf, dt, rain(n) )
        end do
end if

!SOUTHERN BOUNDARY
if( south_bc .eq. 'NO_FLOW') then
        j = 1
        do i = 2, imax - 1
           ! Set porosity factor for this cell
           pf = F_por( h_itr( i, j ) )       
           ! Set the conveyance coefficients 
           CALL Conveyance( 'west ', 'itr', i, j, Cw1 )
           CALL Conveyance( 'east ', 'itr', i, j, Ce1 )
           Cs1 = 0.0  !  <---- NO FLOW BOUNDARY
           CALL Conveyance( 'north', 'itr', i, j, Cn1 )
           ! Fill in the linear system
           CALL set_ABCDEF( i, j, Cw1, Ce1, Cs1, Cn1, pf, dt, rain(n) ) 
        end do
end if


!NORTHERN BOUNDARY
if( north_bc .eq. 'NO_FLOW') then
        j = jmax
        do i = 2, imax - 1
           ! Set porosity factor for this cell
           pf = F_por( h_itr( i, j ) )       
           ! Set the conveyance coefficients 
           CALL Conveyance( 'west ', 'itr', i, j, Cw1 )
           CALL Conveyance( 'east ', 'itr', i, j, Ce1 )
           CALL Conveyance( 'south', 'itr', i, j, Cs1 )
           Cn1 = 0.0  !  <---- NO FLOW BOUNDARY 
           ! Fill in the linear system
           CALL set_ABCDEF( i, j, Cw1, Ce1, Cs1, Cn1, pf, dt, rain(n) ) 
        end do
end if


!NORTH WEST CORNER
if( north_bc .eq. 'NO_FLOW' .AND. west_bc .eq. 'NO_FLOW' ) then
           i = 1; j = jmax
           ! Set porosity factor for this cell
           pf = F_por( h_itr( i, j ) )       
           ! Set the conveyance coefficients 
           Cw1 = 0.0  !  <---- NO FLOW BOUNDARY 
           CALL Conveyance( 'east ', 'itr', i, j, Ce1 )
           CALL Conveyance( 'south', 'itr', i, j, Cs1 )
           Cn1 = 0.0  !  <---- NO FLOW BOUNDARY 
           ! Fill in the linear system
           CALL set_ABCDEF( i, j, Cw1, Ce1, Cs1, Cn1, pf, dt, rain(n) )
end if 

!NORTH EAST CORNER
if( north_bc .eq. 'NO_FLOW' .AND. east_bc .eq. 'NO_FLOW' ) then
           i = imax; j = jmax
           ! Set porosity factor for this cell
           pf = F_por( h_itr( i, j ) )       
           ! Set the conveyance coefficients 
           CALL Conveyance( 'west ', 'itr', i, j, Cw1 )
           Ce1 = 0.0  !  <---- NO FLOW BOUNDARY
           CALL Conveyance( 'south', 'itr', i, j, Cs1 )
           Cn1 = 0.0  !  <---- NO FLOW BOUNDARY 
           ! Fill in the linear system
           CALL set_ABCDEF( i, j, Cw1, Ce1, Cs1, Cn1, pf, dt, rain(n) ) 
end if


! SOUTH WEST CORNER
if(  south_bc .eq. 'NO_FLOW'  .and. &
      west_bc .eq. 'NO_FLOW'          ) then

        i = 1; j = 1
        pf = F_por( h_itr( i, j ) )
        Cw1 = 0.0
        call conveyance( 'east ', 'itr', i, j, Ce1 )
        Cs1 = 0.0
        call Conveyance( 'north', 'itr', i, j, Cn1 )
        call set_ABCDEF( i, j, Cw1, Ce1, Cs1, Cn1, pf, dt, rain(n) )

end if


! SOUTH EAST CORNER
if(  south_bc .eq. 'NO_FLOW'  .and. &
      east_bc .eq. 'NO_FLOW'           ) then

        i = 1; j = 1
        pf = F_por( h_itr( i, j ) )
        call conveyance( 'east ', 'itr', i, j, Cw1 )
        Ce1 = 0.0
        Cs1 = 0.0
        call Conveyance( 'north', 'itr', i, j, Cn1 )
        call set_ABCDEF( i, j, Cw1, Ce1, Cs1, Cn1, pf, dt, rain(n) )

end if


! INTERIOR of DOMAIN
do j = 2, jmax - 1;     do i = 2, imax - 1  
   
    ! set porosity factor for this cell
    pf = F_por( h_itr( i, j ) )       
    ! These things Do change as the iteration progresses
    CALL Conveyance( 'west ', 'itr', i, j, Cw1 )
    CALL Conveyance( 'east ', 'itr', i, j, Ce1 )
    CALL Conveyance( 'south', 'itr', i, j, Cs1 )
    CALL Conveyance( 'north', 'itr', i, j, Cn1 )
    ! Fill in the linear system
    CALL set_ABCDEF( i, j, Cw1, Ce1, Cs1, Cn1, pf, dt, rain(n) )

end do;     end do


! TRANSITION CHECK
!   test to see if there is a transition to or from sheet flow
!   happening during this timestep.   Use under-relaxtion to
!   control oscillations during a transition timestep.

transition = .false.
do j = 1, jmax
    do i = 1, imax
        ! integers used to assure correct behavor when  equal      
        pf_int = nint( F_por( h_old(i,j) ) )  
        pf1_int= nint( F_por( h_itr(i,j) ) )
        if( pf_int .NE. pf1_int ) then
            transition = .true.
        endif
    end do
end do

if( transition .eqv. .true. ) then    
    relaxation_factor = relax_tran
    eps_itr_tol = eps_itr * 10.
else
    relaxation_factor = relax
    eps_itr_tol = eps_itr
endif

 
! CALL diagdom_penta( A, B, C, D, E, n, LB, UB, diagdom) 
! WRITE(100,*) 'Timestep ', n, 'Iteration ', q, &
!               'Is matrix diagonally dominant?', diagdom   


! Confirm that there is a value of C for all of the rows
! this is mostly a check to see that the corner points of
! the domain had values put in.
do v = 1, vmax
    if( abs( C(v) ) .LT.  TINY(  C(v) ) ) then
        write(100,*) ' No value of C: v = ', v, 'C(v)=', C(v)
        write(*,*) 'STOPPING PROGRAM'
        STOP
    end if
end do



!CALL SOLVER  
! gauss_seidel_penta(A,B,C,D,E,F,n,LB,UB,tolit,maxit,Xold,Xnew)  
CALL GAUSS_SEIDEL_penta( A, B, C, D, E, F, vmax, jmax, jmax, eps_matrix, maxit,&
                         h_itr_vec, h_tmp_vec, 100, solver_numits )


! Compute residual and relative change for this iteration.  This took
!  some careful thought to handle both filling and draining cases.
!  Relative change is used when the solution is far from zero
!  and absolute change (residual) is used near zero.


! Should put residual/ relchng computation block into a subroutine.

do v = 1, vmax
 
    if( h_tmp_vec(v) .GT. TINY( h_tmp_vec(v) ) ) then

            ! Compute residual for this iteration
            residual(v) =  h_tmp_vec(v) - h_itr_vec(v)

            ! Handle a result that is effectively zero by
            ! using an absolute tolerance instead of 
            ! a relative one
            if( h_tmp_vec(v) .LE. h_pfc_min  .and. &
                residual (v) .LE. eps_itr_tol     ) then
        
                    relchng(v) = 0.0
            
            else
                    relchng (v) =  residual (v) / h_itr_vec(v)
            endif
            
    elseif( h_tmp_vec(v) .LE. TINY( h_tmp_vec(v) ) ) then

            ! the model is saying the cell is empty, 
            ! so force the solution to be zero      
            h_tmp_vec(v) = 0.0
            ! compute the residual
            residual(v) =  h_tmp_vec(v) - h_itr_vec(v)
            ! For the zero case, use an absolute rather than
            ! relative tolerance by setting the value of relchng
            ! below the tolerance instead of computing it.
            if( abs( residual(v) ) .LE. eps_itr_tol ) then
                    
                    relchng(v) = 0.0
            endif
    endif

end do
            

! Store solution history during iteration in case 
!  the model fails to converge
h_temp_hist( :, q ) = h_tmp_vec

  
!Output the biggest change for this iteration
WRITE(100,*) 'PERFCODE: Iteration q =',          q                 , &
                'Solver Interations =',     solver_numits          , &
                       'L_inf_norm  =', maxval( abs( relchng ) )   , &
                        'At Cell v  =', maxloc( abs( relchng ) )   , &
                         ' L2 Norm  =', F_L2_NORM( relchng, vmax ) , &
                       'eps_itr_tol =', eps_itr_tol

! CONVERGENCE TEST
! Exit iteration loop if this timestep has converged
if( maxval( abs ( relchng ) )    .le. eps_itr_tol  .AND. &
    F_L2_NORM   ( relchng, vmax) .le. eps_itr_tol           ) then    
    WRITE(100,*) 'Time step n = ', n, time(n),'sec ' , &
                     'rain(n) = ', rain(n)           , &
             ' converged in q = ', q, ' iterations.' , &
             ' maxdepth=', maxval( h_tmp_vec ),        &
             ' max 1D =', maxval( h_new_1D ) , 'min 1D =', minval( h_new_1D)            

    WRITE(100,*)''
    !output results for each timestep for checking purposes
!    write(50,2) n, h_itr_vec(:)
    EXIT iteration
endif



!update iteration variables
h_itr_vec = h_itr_vec + relaxation_factor * residual



! un-linearize the thicknesses back to a matrix h_itr_vec ---> h_itr 
call unlinearize( h_itr_vec, imax, jmax, vmax, h_itr )


end do iteration


!Give Error if Iteration fails to converge and write some diagnostics
if (q .gt. qmax) then
    WRITE(*,*) ' Iteration failed to converge for time level n = ', n
    !output the coefficient matrix and main diagonal
    call write_system( A, B, C, D, E, F, vmax, 'ABCDEF.csv' )
    call write_flipped_matrix( h_old, imax, jmax, 'h_old.csv' )
    call write_matrix( h_temp_hist, vmax, qmax, 'h_temp_hist.csv')
    call WRITE_VECTOR( residual, vmax, 'residual_iteration.csv')
    call WRITE_VECTOR( relchng, vmax, 'relchng_iteration.csv')
!    call put_bands(a, b, c, d, e, vmax, lb, ub, amatrix)
!    call write_matrix( amatrix, vmax, vmax, 'amatrix.csv')
    EXIT time_stepping
end if


! Compute Change for this time step
!Time stepping residual (re-uses the arrays) 
residual = h_tmp_vec - h_old_vec

! compute relative change for this timestep
do v = 1, vmax
    if( abs(residual(v)) .LT. TINY(residual(v)) ) then
            ! The converged solution is zero
            relchng(v) = 0.0
    else
            !the solution is non-zero, compute as ususal
            relchng(v) = residual(v) / h_old_vec(v) 
    endif
end do

maxrelchng_ss = maxval ( ABS( relchng ) )

!call WRITE_VECTOR( relchng, vmax, 'relchng_time.csv')



!Update the old and new solutions
!At the end of the iteration, we have found values for the
!next time step. 
h_new_vec  = h_tmp_vec

!but when we go back to the top of the loop, the old is what we just found 
h_old_vec  = h_new_vec
!and now we need to unlinearize the h_old values
call unlinearize( h_old_vec, imax, jmax, vmax, h_old )






!--------------------------------------------------------------------
!  Summary Info for this timestep
!-------------------------------------------------------------------

numit         ( n ) = q   
loc           ( n ) = maxloc ( abs( relchng ), dim = 1 )     
maxdiff       ( n ) = relchng ( loc ( n ) ) 
maxthk        ( n ) = maxval( h_old_vec )
L2_History    ( n ) = F_L2_Norm( relchng, vmax )
h_imid_j1_hist( n ) = h_old( imax/2, 1 ) 
h_imid_max_hist(n)  = maxval( h_old( imax/2 , :) )

! Compute the flow into the southern boundary for this time step
!  (assume that we can neglect the drainage area of the last row)
j = 2
do i = 1, imax
    CALL Conveyance( 'south', 'itr', i, j, Cs1 )
    Qout(n) = Qout(n) + Cs1 * area(i,j) *       &
              (  ( h_itr(i, j-1) - h_itr(i,j) ) &
               + (     Z(i, j-1) -     Z(i,j) )   )
end do


!SELECTIVELY STORE MODEL RESULTS
! MAXIMUM DEPTH
! Check to see if this was the maximum time-step and store if so
if( maxval( h_old_vec) .GT. maxval( h_max ) ) then
        call unlinearize( h_old_vec, imax, jmax, vmax, h_max )
endif

! MAXIMUM DISCHARGE
if( Qout(n) .GT. maxval( Qout(1:n-1) ) ) then
        call unlinearize( h_old_vec, imax, jmax, vmax, h_Q_max )
endif

! MAXIMUM MID DOMAIN DISCHARGE DEPTH
if( h_imid_j1_hist( n ) .GT. maxval( h_imid_j1_hist(1:n-1) )) then
        call unlinearize( h_old_vec, imax, jmax, vmax, h_imid_j1_max )
endif

! MAXIMUM MID DOMAIN DISCHARGE DEPTH
if( h_imid_max_hist( n ) .GT. maxval( h_imid_max_hist(1:n-1) )) then
        call unlinearize( h_old_vec, imax, jmax, vmax, h_imid_max )
endif


! ANIMATION
!  Decide if the results from this timestep should be stored for
!  animation output.  Take the time, divide by the animation step,
!  round to the lowest integer and then convert to integer
if( animate .eqv. .true. ) then

    if( int( floor( time(n) / dt_ani ) )  .gt. ani ) then
        ! set the value of ani
        ani = ani + 1
print *, 'n = ', n, 'ani=', ani
        ! store the solution for this step
        h_vec_ani( :, ani ) = h_old_vec
        ! also store a label
        write( sim_time2, 123 ) time_simulated
        ani_lab( ani ) = 'h'//sim_time2//'s' 
        ani_time(ani ) = time_simulated
    endif

endif 



!STEADY-STATE CHECK  ( disabled in favor of setting
!                      the time for the simulation to run )
!IF (          maxrelchng_ss .le. eps_ss .AND. & 
!     F_L2_NORM( relchng, v) .le. eps_ss         ) then  
!  WRITE(*,*) 'Simulation reached steady state after', n, &
!               &'time steps or', time_simulated, 'seconds'
!  EXIT time_stepping
!end if  

end do time_stepping

!for outputting each timestep
! close(50)

!close log file
close(100)

!-------------------------------------------------------------------------------
!   >>>>>>>>>>      P O S T   P R O C E S S I N G       <<<<<<<<<<
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
!   >>>>>>>>>>     W R I T E    O U T P U T    F I L E S        <<<<<<<<<<
!-------------------------------------------------------------------------------
!Set date and time stamps

CALL SYSTEM_CLOCK( RUN_END_TIME, COUNT_RATE, COUNT_MAX )
call DATE_AND_TIME(FILE_DATE,FILE_TIME)
call CPU_TIME(cputime)

!-------------------------------------------------------------------------------
! file to show 1D solution along i = imax / 2; j = 1:jmax

OPEN(UNIT = 10, FILE = 'PERFCODE.csv', STATUS='REPLACE')
WRITE(10,*) 'Output From PERFCODE.f95'
WRITE(10,*) 'Timestamp,', FILE_DATE,' ', FILE_TIME,','
do i = 1, 18
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
             &1D MODEL OUTPUT IN [ SI ] UNITS &
             &***********************************,'
i = imax / 2       
write(10,*) ' i = ', i,','
write(10, *) 'j,eta,Z,PFC_Surf,h,Head,Surf_Thk.mm,'
do j = 1, jmax
    v = F_LinearIndex( i, j, jmax)
    write(10, 2) j, CV_Info(v)%eta, Z(i,j), Z(i,j) + b_pfc,      &
                                h_old(i,j), Z(i,j) + h_old(i,j), &
                              ( h_old(i,j) - b_pfc ) * 1000.
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



!--------------------------------------------------------------------------------
! Write parameters to a seperate file for convenicence
open( unit = 15, file = 'params.csv', status = 'REPLACE' )
write( 15, 155 ) input_variables(:), 'north_bc', 'south_bc', 'east_bc', 'west_bc'
write( 15, 156 ) input_values(:), north_bc, south_bc, east_bc, west_bc
close(15)

155  format (  22( A, ',') )
156  format (  18( E, ','), 4 ( A, ',') )
!-------------------------------------------------------------------------
!Write time history to a file
!  ( hydrographs and anything else time-dependant 
!    is plotted from this file )

OPEN( UNIT = 20, FILE = 'details.csv', STATUS='REPLACE')
WRITE(20,*) 'Timestamp,', FILE_DATE, ' ', FILE_TIME, ','
DO i = 1, 18 
    WRITE( 20, * ) input_variables(i), ',', input_values(i), ','
END DO
write(20,*) 'north_bc,', north_bc       
write(20,*) 'south_bc,', south_bc 
write(20,*) 'east_bc,', east_bc
write(20,*) 'west_bc,', west_bc
WRITE(20,*) 'imax,', imax, ','
WRITE(20,*) 'jmax,', jmax, ','
WRITe(20,*) 'vmax,', vmax, ','
WRITE(20,*) '-----,'
WRITE(20,*) 'Timestep,Iterations,MaxRelChng,MaxLocn,' , &
                                'L2_Norm,Rain.mmphr,' , &
                           'MaxThk.cm,Time,Qout.Lps,' , &
                                    'h_imid_j1_hist,' , &
                                    'h_imid_max_hist,' 
DO n = 1, nlast
WRITE(20,300) n, numit(n), maxdiff(n), loc(n)         , &
              L2_History(n), rain(n)*1000.*3600.      , &
              maxthk(n)*100., time(n), -Qout(n)*1000. , &
              h_imid_j1_hist(n), h_imid_max_hist(n)
end do
close(20)    


!-----------------------------------------------------------------------------------
! Output depth grid for last timestep

! an internal write statement to store the value of the REAL variable
! "time_simulated" in the CHARACTER variable "out_time"
write( out_time, 111 ) time_simulated

call write_flipped_matrix( h_old, imax, jmax, 'h_old'//out_time//' sec.csv' )

!------------------------------------------------------------------------------
! Output iteration history for the last time-step 

call write_matrix( h_temp_hist, vmax, qmax, 'h_temp_hist'//out_time//' sec.csv')

!-------------------------------------------------------------------------
! Animation output

if( animate .eqv. .TRUE. ) then

!Animation results
open( unit = 70, file = 'animate.csv', status = 'REPLACE' )
write( 70, 700) 'v,X,Y,Z,', ani_lab(:)
do j = 1, jmax
    do i = 1, imax
        v = F_LinearIndex( i, j, jmax)
        write(70, 2) v, CV_Info( v ) % X,                    &
                        CV_Info( v ) % Y, Z(i,j), h_vec_ani( v, :)
    end do
end do
close( 70 )

700 format( (A, 10000( A, ',') ) )

!Also sperately output the list of animation lables
open( unit = 71, file = 'ani_labs.csv', status = 'REPLACE' )
write( 71, *) 'ani,lab,time,'
do ani = 1, animax
    write( 71, 711 ) ani, ani_lab(ani), ani_time(ani)
end do
close( 71 )

end if


711 format( (I, ','), (A, ','), (F8.2, ',') )

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

!----------------------------------------------------------------------------------
!Format statements

2       FORMAT( I, ',', 10000 ( E, ',') )
10      FORMAT('   ', ( i3, '    '),  ( F10.3, '  ') , F10.6 ) 
111     FORMAT( f9.2 )
200     FORMAT ( A, ( E, ',') ) 
201     FORMAT ( A, ( I, ',') ) 
300     FORMAT ( 2 (     I, ','),      F12.7, ',' , &    ! n, numit, maxdif
                         I, ',' ,          E, ',' , &    ! loc, L2_History
                 2 ( F12.8, ','), ( F12.3, ',' ), 3 ( F12.8 ,',')     )  ! rain, maxthk, time, Qout, h_imid_j1hist, h_imid_max_hist
400     FORMAT( 10000 ( I, ',' ) )
401     FORMAT( (I, ',') , 2( F12.7, ',' ) )
660     FORMAT(  2( I, ','), 2( F12.7, ',') )
!------------------------------------------------------------------------------
end program PERFCODE
!===================================================================================
!       \\\\\\\\\\          E N D    P R O G R A M          //////////
!       //////////             P E R F C O D E              \\\\\\\\\\
!===================================================================================




