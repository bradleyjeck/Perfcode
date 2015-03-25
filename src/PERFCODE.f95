! fortran_free_source
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
!             Date:     April 2010  -- original coding
!                       July  2011  -- curb boundary condition, 
!                                      specify properties by cell
!                                      global water balance,
!                                      
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
! Alphabetical list of variables declared in the module SHARED.
!       (variables used in subroutines are described there)
!
!   A         -- lowest band of penta diagonal matrix 
!   area      -- area of a grid cell
!   astat     -- array allocation statuses
!   animate   -- logical flag for whether to store animation output
!   B         -- subdiadonal band of penta diagonal matrix
!   b_pfc     -- thickness of the PFC layer for each cell
!   b_pfc_input--PFC thickness read from input file
!   BCs_vary  -- logical noting whether BCs vary by cell
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
!   dry       -- integer code 
!   ds        -- distance up characteristic in sheet flow moc bc
!   dt        -- time step for the simulation
!   dt_ani    -- time step for storing animation results
!   dt_pfc    -- time step for PFC flow
!   dt_sheet  -- time step for sheet flow
!   dx        -- prelim. grid size for longitudinal direction
!   dx_moc    -- distance up drainage path in pfc moc bc
!   dy        -- prelim. grid size for transverse direction
!   E         -- uppermost band of penta diagonal matrix
!   east_bc   -- condition for each cell on east boundary
!   east_bc_input--single conditino for all cells on east boundary
!   east_vol  -- volume that leaves through east boundary during a time step
!   EDX       -- width of each cell in 1D model
!   eta_0_hp2_max-- max possible value for pfc moc bc
!   eta_cs    -- values of eta along the cross slope
!   eta_cs_1D -- values of eta for 1D model
!   eta_sheet -- location of the transition to sheet flow (1D model)
!   etaCV     -- value of eta at CV center for cells in 1D model
!   eps_itr   -- tolerance for an iteration
!   eps_itr_tol-- selected tolerance for the iteration (based on transition)
!   eps_matrix-- tolerance (epsilon) for matrix solver
!   eta_sheet -- the location of sheet flow for each time step
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
!   j_imid_j1_hist-- time history of the depth at the i=imax/2 j=1 cel 
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
!   hyd_cond  -- hydraulic conductivity for each grid cell (formarly K)
!   i         -- array index ( longitudinally in the domain )
!   input_values-- array of values of the input variables
!   input_variables -- character array of input variables
!   imax      -- maximum value of the array index i
!   is_boundary--logical matrix noting if a cell is on the domain boundary
!   is_corner -- logical matrix noting if a cell is on a domain corner
!   itr       -- integer code for time level n
!   j         -- array index ( transverse in the domain )
!   jmax      -- maximum value of the array index j
!   K_input   -- value of hydraulic conductivity read from parameters input file
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
!   model_mode-- whether the model runs in 1D or 2D mode
!   n         -- index for time stepping
!   n_mann    -- Manning's roughness coefficient for each grid cell
!   n_mann_input--Manning's n read from parameters input file
!   net_flow  -- net flow volume during a a time step
!   North     -- integer aziuth for north direction
!   NorthEast -- integer azimuth for direction
!   NorthWest -- integer azimuth for direction
!   north_bc  -- boundary condition for cells on north boundary of domain
!   north_bc_input--single condition for north BC, read from paramter input file
!   north_vol -- volume through the north boundary during a time step
!   nlast     -- last timestep taken
!   nmax      -- maximum number of time steps in the simulation
!   numit     -- the number of iterations required for a timestep to converge
!   nr_cs     -- number of records in the cross slope file
!   nr_lp     -- number of records in the longitudinal profile file
!   nrr       -- number of rainfall records
!   old       -- integer code for time level n-1
!   old_vol   -- system volume at previous time step
!   out_time  -- 
!   params_vary -- logical noting if pfc parameteters vary by grid cell
!   pf        -- porosity factor ( includes effect of porosity 
!                                   when pavement is not saturated )
!   pf_int    -- porosity factor as an integer
!   pf1       -- porosity factor for time level n+1
!   pf1_int   --  "      "                      ""   as integer
!   pfc_vol   -- volume of water in the pfc
!   por       -- the effective porosity of each cell in the grid  
!   por_input -- porosity read from paramters input file
!   q         -- iteration index
!   qmax      -- maximum number of iterations
!   rain      -- rainfall rate for each timestep of the simulation
!   q_pav     -- unit flux in the pavement for 1D model
!   q_surf    -- unit flux on the surface for 1D model
!   Qout      -- flow rate out the southern boundary for a timestep
!   rain_rate -- rainfall rate for each time increment in the 
!                   rainfall input file
!   rain_time -- time column of rainfall input file
!   rain_vol  -- volume of rainfall during a time step
!   relax     -- relaxation factor for non-transition iterations
!   relaxation_factor -- underrelaxation factor for non-linear iteration
!   relax_tran -- relaxation factor for transition
!   relchng   -- the relative change between solutions for an iteration or timestep
!   report    -- report results of this step to screen
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
!   South     -- integer azimuth for direction
!   SouthEast -- integer azimuth for direction
!   south_bc  -- condition for each cell south boundary
!   south_bc_input--single conditino for south boundary, read from param input
!   south_vol -- volume that leaves the south face during a timestep
!   storage_change-- volumetric change in storage during a time step
!   surface_vol--volume of water on the road surface
!   time      -- time at each timestep
!   time_simulated-- the time simulated
!   timestep_solver_numits -- 
!   transition -- logical to see if we're in a transition timestep
!   tolit     -- tolerence for iterations, used for relative (fractional) changes
!   TNE       -- total number of elements for 1D grid
!   v         -- linear index for domain
!   ve        -- linear index for cell to the east
!   vol_error -- difference between net_flow and storage_change volumes
!   vol_error_frac--volume error as a fraction of the storage_change
!   v_in      -- linear index of adjacent inside cell
!   vmax      -- number of unknowns in the domain   
!   west      -- integer azimuth for direction
!   west_bc   -- condition for each cell on west boundary
!   west_bc_input--single condition for west BC, read from param input file
!   west_vol  -- volume that leaves through west boundary during timestep
!   wid       -- curvilinear width of a grid cell at its center
!   wid_cs    -- width column of cross slope file
!   wid_cs_1d -- width of 1D segment
!   XCV       -- coordinate of CV center for 1D grid
!   X_sheet   -- location of transition to sheet flow (1D model)
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
USE applyBCs
use budget

!-------------------------------------------------------------------------
!   >>>>>>>>>>          V A R I A B L E S                  <<<<<<<<<<
!-------------------------------------------------------------------------

implicit none   

!  All variables are declared in module SHARED
      
!-------------------------------------------------------------------------------------
!   >>>>>>>>>>      P R O B L E M    S E T U P          <<<<<<<<<<
!-------------------------------------------------------------------------------------

! Disply program header 
call HEADER ( dev = 6 )

! Create a file to store details of the run
open( unit = 100, file = 'PERFCODE_Run.txt', status = 'REPLACE' )

! Problem parameters file
CALL GET_PARAMETERS(       )  

! Rainfall file & maximum number of timesteps
call GET_RAINFALL( max_rec, rain_time, rain_rate, nrr )
nmax = ( maxval( rain_time(1:nrr) ) / min( dt_pfc, dt_sheet ) ) 


! GRID GENERATION

! Read centerline geomtry file and generate horizontal grid
call GENERATE_GRID( prelim_dx = dx , prelim_dy = dy )

! Read cross section and longitudinal profile and computes elevations
call SET_ELEVATIONS( )


!Output grid information for 2D simulations
if( model_mode .eq. '2D' ) then

  call write_cv_info( )
     
  CALL WRITE_FLIPPED_MATRIX( Z, imax, jmax, 'Z.csv')

  call write_flipped_matrix(  lng, imax, jmax, 'length.csv' )

  call write_flipped_matrix(  wid, imax, jmax,  'width.csv' )

  call write_flipped_matrix( area, imax, jmax,   'area.csv' )

  call write_flipped_matrix( lng_south, imax, jmax, 'lng_south.csv')

end if

! Creates a grid for a 1D section in case a 1D boundary condition is used
call setup_1d_section( )
CALL grid_1d_section( slope_in = slope_cs_1D  ,  &
                      width_in =   wid_cs_1D  ,  &
                           seg =       nr_cs  ,  &
                            dx = (( dx+dy ) / 2.)   )  


! ALLOCATE ARRAYS
call ALLOCATE_ARRAYS( device = 6 )

!----------------------------------------------------------
! VARIABLE PARAMETERS & BOUNDARY CONDITIONS (2D only)
!-------------------------------------------------------

if( model_mode .eq. '2D' ) then
  ! Deal with pfc parameters that vary by cell
  call cellwise_pfc_params( params_vary )

  ! Deal with boundary conditions that vary by cell
  call cellwise_BCs( BCs_vary )

end if


call CHECK_BCs( model_mode )

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
                     'relax     ', &
                     'relax_tran' /) 

!also collect and store values of input variales
input_values =  (/ K_input, por_input, b_pfc_input, n_mann_input,      &
                   g, dt_pfc, dt_sheet,                                &
                   max_time, dx, dy, real(qmax), real(maxit),          &
                   h0, eps_matrix, eps_itr, relax, relax_tran /)


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

write(*,*) 'PERFCODE: Model mode is ', model_mode

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

!------------------------
! ONE-DIMENSIONAL MODEL
!-------------------------

if1D: if( model_mode .eq. '1D' ) then

write( *,*) 'entering 1D model', model_mode

stepping_1D: do while( time_simulated .LT. max_time )

!increment n and store the largest n we've gotten so far
n = n + 1
nlast = n

! Select the time step 
if( maxval( h_old_1D ) .GT. b_pfc_input * 0.95 ) then
    dt = dt_sheet
else
    dt = dt_pfc
endif

!Computed the time simulated  
!   Do the accumulation with an internal write/read to 
!    avoid accumulating the floating point errors
write( sim_time, 123 ) time_simulated
read( sim_time,  * ) time_simulated


time_simulated = time_simulated + dt
time(n) = time_simulated
!Report which timestep we're in every 20 or so time steps
if( nint( real(n)/20. ) .gt. report ) then
    report = report + 1
    write(*,*) ' n = ', n-1, ' time = ', time(n-1),          &
!                      'L_inf_norm = ', maxrelchng_ss,           &
!                         'L2_norm = ', F_L2_Norm(relchng,vmax), &
                            'Qout = ', Qout(n-1)
endif


!Come up with the rainfall rate for this timestep
rain(n) = F_Linterp( time_simulated  , &
                     rain_time(1:nrr), &
                     rain_rate(1:nrr), &
                     nrr                  )


h_old_1d = h_new_1d

CALL PFC1DIMP( h_old = h_old_1d      , &
                  dt = dt            , &
                rain = rain(n)       , & 
               tolit = eps_itr       , &
                qmax = qmax          , & 
               h_new = h_new_1d      , &
                imax = TNE           , &
            eta_0_BC = south_bc_input, &
            eta_1_BC = north_bc_input, & 
                Qout = Qout(n)       , &
             X_sheet = X_sheet       , &
               numit = numit(n)           )


! Check the solution for weird problems
if( maxval( h_new_1d) .LT. TINY(h_new_1d(1) ) ) then
   write(100,*) 'PERFCODE: 1D Model zeroed out....stopping program'
   call write_vector( h_old_1d, TNE, 'h_old_1D.csv' )
   call write_vector( h_new_1d, TNE, 'h_new_1D.csv' )
   stop
end if


call MIXED_RESIDUAL( old     = h_old_1D        , &
                     new     = h_new_1D        , &
                     vmax    = TNE             , & 
                     res     = residual(1:TNE) , &
                     mixres  = relchng(1:TNE)  , &
                     abs_tol = h_pfc_min       , &
                     rel_tol = eps_itr_tol         ) 

! Summary info for this 1D time step
loc           ( n ) = maxloc ( abs( relchng(1:TNE) ), dim = 1 )     
maxdiff       ( n ) = relchng ( loc ( n ) ) 
maxthk        ( n ) = maxval( h_new_1D )
L2_History    ( n ) = F_L2_Norm( relchng(1:TNE), TNE )
h_imid_j1_hist( n ) = h_new_1d( TNE ) 
h_imid_max_hist(n ) = maxval( h_new_1D ) 
eta_sheet     ( n ) = 1 - X_sheet / sum( wid_cs_1d(:) ) 

end do stepping_1D
end if if1D


!--------------------------
! TWO-DIMENSIONAL MODEL
!--------------------------

if2D: if( model_mode .eq. '2D' ) then
time_stepping: do while (time_simulated .LT. max_time )

!increment n and store the largest n we've gotten so far
n = n + 1
nlast = n

! Select the time step 
if( maxval( h_old  ) .GT.  b_pfc_input * 0.95 ) then
    dt = dt_sheet
else
    dt = dt_pfc
endif

!Computed the time simulated  
!   Do the accumulation with an internal write/read to 
!    avoid accumulating the floating point errors
write( sim_time, 123 ) time_simulated
read( sim_time,  * ) time_simulated

time_simulated = time_simulated + dt
time(n) = time_simulated

!Report values from the last timestep every 20 or so time steps
if( nint( real(n)/20. ) .gt. report ) then
    report = report + 1
    write(*,*) ' n = ', n-1, ' time = ', time(n-1),          &
                      'L_inf_norm = ', maxrelchng_ss,           &
                         'L2_norm = ', F_L2_Norm(relchng,vmax), &
                            'Qout = ', Qout(n-1)
endif


!Come up with the rainfall rate for this timestep
rain(n) = F_Linterp( time_simulated  , &
                     rain_time(1:nrr), &
                     rain_rate(1:nrr), &
                     nrr                  )

!Compute solution for 1D model only if if it was used as a BC 
if( west_bc_input .eq. '1D_FLOW' .or. &
    east_bc_input .eq. '1D_FLOW'        ) then

       h_old_1d = h_new_1d

       CALL PFC1DIMP(  h_old = h_old_1d      , &
                          dt = dt            , &
                        rain = rain(n)       , & 
                       tolit = eps_itr       , &
                        qmax = qmax          , & 
                       h_new = h_new_1d      , &
                        imax = TNE           , &
                    eta_0_BC = south_bc_input, &
                    eta_1_BC = north_bc_input    )


       ! Check the solution to avoid a weird problem
       if( maxval( h_new_1d) .LT. TINY(h_new_1d(1) ) ) then
           write(100,*) 'PERFCODE: 1D Model zeroed out....stopping program'
           call write_vector( h_old_1d, TNE, 'h_old_1D.csv' )
           call write_vector( h_new_1d, TNE, 'h_new_1D.csv' )
           stop
       end if

endif


!------------------
! STATIONARY PART 
!------------------
!PART OF NON-LINEAR SYSTEM FROM TIME LEVEL n
!   FRICTION SLOPE
!   Compute friction slope magnitudes based on the converged thicknesses
!   from the previous time step
CALL FrictionSlope( old, Sfw_old, Sfe_old, Sfs_old, Sfn_old )

!   BCs, CONVEYANCE COEFS and RHS   
!   Loop over the whole domain applying BCs and computing
!   conveyance coefs 
do j = 1, jmax
    do i = 1, imax

        if( is_boundary( i, j) .eqv. .true. ) then

            CALL SET_BC( old ) 

        elseif( is_boundary( i, j) .eqv. .false. ) then
            ! Set porosity factor for this cell
            pf  = F_por( h_old(i, j), b_pfc(i, j), por(i, j) )
           
            ! Set the conveyance coefficients
            CALL Conveyance( west , old, i, j, Cw )
            CALL Conveyance( east , old, i, j, Ce )
            CALL Conveyance( south, old, i, j, Cs )
            CALL Conveyance( north, old, i, j, Cn )
           
            ! Compute the part of the right-hand-side that is from 
            ! time level n
            v     = F_LinearIndex( i, j, jmax)
            Fn(v) = F_RHS_n( i, j, Cw, Ce, Cs, Cn, rain(n-1), pf, dt )

        end if
    end do
end do

!----------------------------------------------
!ITERATIVE (LAGGED) PART OF NON-LINEAR SYSTEM 
!----------------------------------------------

!zero out matrix iteration counter
timestep_solver_numits = 0

iteration: do q = 1, qmax

!   FRICTION SLOPE
!   compute friction slope magnitudes based on the thickness
!   from the previous iteration
CALL FrictionSlope( itr, Sfw_itr, Sfe_itr, Sfs_itr, Sfn_itr )

do j = 1, jmax
    do i = 1, imax   

        ! Figure out if this is a boundary cell
        if( is_boundary( i, j) .eqv. .true. ) then

            call SET_BC( itr ) 
        
        elseif( is_boundary( i, j) .eqv. .false. ) then

            ! These things Do change as the iteration progresses
            CALL Conveyance( west , itr, i, j, Cw1 )
            CALL Conveyance( east , itr, i, j, Ce1 )
            CALL Conveyance( south, itr, i, j, Cs1 )
            CALL Conveyance( north, itr, i, j, Cn1 )
        
            ! set porosity factor for this cell
            pf = F_por( h_itr( i, j ), b_pfc(i, j), por(i, j) )       

            ! Fill in the linear system
            CALL set_ABCDEF( i, j, Cw1, Ce1, Cs1, Cn1, pf, dt, rain(n) )

        end if
    end do
end do

! TRANSITION CHECK
!   test to see if there is a transition to or from sheet flow
!   happening during this timestep.   Use under-relaxtion to
!   control oscillations during a transition timestep.

transition = .false.
do j = 1, jmax
    do i = 1, imax
        ! integers used to assure correct behavor when  equal      
        pf_int = nint( F_por( h_old(i,j), b_pfc(i, j), por(i, j) ) )  
        pf1_int= nint( F_por( h_itr(i,j), b_pfc(i, j), por(i, j) ) )
        if( pf_int .NE. pf1_int ) then
            transition = .true.
        endif
    end do
end do

if( transition .eqv. .true. ) then    
    relaxation_factor = relax_tran
    eps_itr_tol = eps_itr  * 10.
else
    relaxation_factor = relax
    eps_itr_tol = eps_itr
endif

 
! Confirm that there is a value of C for all of the rows
! this is mostly a check to see that the corner points of
! the domain had values put in.
do v = 1, vmax
    if( abs( C(v) ) .LT.  TINY(  C(v) ) ) then
        write(*,*) 'PERFCODE: No value of C: v = ', v, 'C(v)=', C(v)
        write(*,*) '          STOPPING PROGRAM'
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
call MIXED_RESIDUAL( old     = h_itr_vec   , &
                     new     = h_tmp_vec   , &
                     vmax    = vmax        , & 
                     res     = residual    , &
                     mixres  = relchng     , &
                     abs_tol = h_pfc_min   , &
                     rel_tol = eps_itr_tol    ) 
          

! Store solution history during iteration in case 
!  the model fails to converge
h_temp_hist( :, q ) = h_tmp_vec


! Diagnostic output to the file PERFCODE_RUN.txt
!  
!        !Output the biggest change for this iteration
!        WRITE(100,*) 'PERFCODE: Iteration q =',          q                 , &
!                        'Solver Interations =',     solver_numits          , &
!                               'L_inf_norm  =', maxval( abs( relchng ) )   , &
!                                'At Cell v  =', maxloc( abs( relchng ) )   , &
!                                 ' L2 Norm  =', F_L2_NORM( relchng, vmax ) , &
!                               'eps_itr_tol =', eps_itr_tol
!


! CONVERGENCE TEST 
! Exit iteration loop if this timestep has converged
if( maxval( abs ( relchng ) )    .le. eps_itr_tol  .AND. &
    F_L2_NORM   ( relchng, vmax) .le. eps_itr_tol           ) then    

            EXIT iteration
endif

!update iteration variables
h_itr_vec = h_itr_vec + relaxation_factor * residual

! un-linearize the thicknesses back to a matrix h_itr_vec ---> h_itr 
call unlinearize( h_itr_vec, imax, jmax, vmax, h_itr )


! end of non-linear iteration loop
end do iteration

!Give Error if Iteration fails to converge and write some diagnostics
if (q .gt. qmax) then
    WRITE(*,*) ''
    WRITE(*,*) 'PERFCODE:       ITERATION FAILED TO CONVERGE   for time level n = ', n
    write(*,*) ''
    !output diagnostic information  
    call write_system( A, B, C, D, E, F, vmax, 'ABCDEF.csv' )
    call write_flipped_matrix( h_old, imax, jmax, 'h_old.csv' )
    call write_matrix( h_temp_hist, vmax, qmax, 'h_temp_hist.csv')
    call WRITE_VECTOR( residual, vmax, 'residual_iteration.csv')
    call WRITE_VECTOR( relchng, vmax, 'relchng_iteration.csv')
    EXIT time_stepping
end if


! The converged solution at this time step this the version of h_tmp_vec
!  the exited the iteration loop. 
h_new_vec = h_tmp_vec
call unlinearize( h_new_vec, imax, jmax, vmax, h_new )

! 
! ! Compute Change for this time step
! !Time stepping residual (re-uses the arrays) 
! residual = h_new_vec - h_old_vec
! 
! ! compute relative change for this timestep
! do v = 1, vmax
!     if( abs(residual(v)) .LT. TINY(residual(v)) ) then
!             ! The converged solution is zero
!             relchng(v) = 0.0
!     else
!             !the solution is non-zero, compute as ususal
!             relchng(v) = residual(v) / h_old_vec(v) 
!     endif
! end do
! 
! maxrelchng_ss = maxval ( ABS( relchng ) )


call MIXED_RESIDUAL( old     = h_old_vec   , &
                     new     = h_new_vec   , &
                     vmax    = vmax        , & 
                     res     = residual    , &
                     mixres  = relchng     , &
                     abs_tol = h_pfc_min   , &
                     rel_tol = eps_itr_tol    )

maxrelchng_ss = maxval( abs( relchng) ) 

!--------------------------------------------------------------------
!  Summary Info for this timestep
!-------------------------------------------------------------------

numit         ( n ) = q   
loc           ( n ) = maxloc ( abs( relchng ), dim = 1 )     
maxdiff       ( n ) = relchng ( loc ( n ) ) 
maxthk        ( n ) = maxval( h_new_vec )
L2_History    ( n ) = F_L2_Norm( relchng, vmax )
h_imid_j1_hist( n ) = h_new( imax/2, 1 ) 
h_imid_max_hist(n)  = maxval( h_new( imax/2 , :) )

! Compute the flow into the southern boundary for this time step
!  (assume that we can neglect the drainage area of the last row)
j = 2
do i = 1, imax
    CALL Conveyance( south, itr, i, j, Cs1 )
    Qout(n) = Qout(n) + Cs1 * area(i,j) *       &
              (  ( h_new(i, j-1) - h_new(i,j) ) &
               + (     Z(i, j-1) -     Z(i,j) )   )
end do

! GLOBAL WATER BUDGET
! set old volume
if( n == 1 ) then
    old_vol = sum(area(:,:) * h0 * por(:,:) )
else

    old_vol = pfc_vol(n-1) + surface_vol(n-1)

end if

CALL WATER_BUDGET(    old_volume =     old_vol      ,  &
                      pfc_volume =     pfc_vol   (n),  &
                  surface_volume = surface_vol   (n),  &
                     rain_volume =    rain_vol   (n),  &
                north_out_volume =   north_vol   (n),  &
                south_out_volume =   south_vol   (n),  &
                 east_out_volume =    east_vol   (n),  &
                 west_out_volume =    west_vol   (n),  &
                 net_flow_volume =     net_flow  (n),  &
              storage_chg_volume = storage_change(n),  &
                    volume_error =     vol_error (n),  &
           volume_error_fraction = vol_error_frac(n)      )



!SELECTIVELY STORE MODEL RESULTS
! MAXIMUM DEPTH
! Check to see if this was the maximum time-step and store if so
if( maxval( h_new_vec) .GT. maxval( h_max ) ) then
        call unlinearize( h_new_vec, imax, jmax, vmax, h_max )
endif

! MAXIMUM DISCHARGE
if( Qout(n) .GT. maxval( Qout(1:n-1) ) ) then
        call unlinearize( h_new_vec, imax, jmax, vmax, h_Q_max )
endif

! MAXIMUM MID DOMAIN DISCHARGE DEPTH
if( h_imid_j1_hist( n ) .GT. maxval( h_imid_j1_hist(1:n-1) )) then
        call unlinearize( h_new_vec, imax, jmax, vmax, h_imid_j1_max )
endif

! MAXIMUM MID DOMAIN DISCHARGE DEPTH
if( h_imid_max_hist( n ) .GT. maxval( h_imid_max_hist(1:n-1) )) then
        call unlinearize( h_new_vec, imax, jmax, vmax, h_imid_max )
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
        h_vec_ani( :, ani ) = h_new_vec
        ! also store a label
        write( sim_time2, 123 ) time_simulated
        ani_lab( ani ) = 'h'//sim_time2//'s' 
        ani_time(ani ) = time_simulated
    endif

endif 



! Very last piece of time stepping loop: 
!  update the 'old' soln
h_old_vec  = h_new_vec
!and now we need to unlinearize the h_old values
call unlinearize( h_old_vec, imax, jmax, vmax, h_old )

end do time_stepping
end if if2d

!for outputting each timestep
! close(50)

!close log file
close(100)

!-------------------------------------------------------------------------------
!   >>>>>>>>>>     W R I T E    O U T P U T    F I L E S        <<<<<<<<<<
!-------------------------------------------------------------------------------
!Set date and time stamps

CALL SYSTEM_CLOCK( RUN_END_TIME, COUNT_RATE, COUNT_MAX )
call DATE_AND_TIME(FILE_DATE,FILE_TIME)
call CPU_TIME(cputime)

!--------------------------------------------------------------------------------
! Write parameters to a seperate file for convenicence
open( unit = 15, file = 'params.csv', status = 'REPLACE' )
write( 15, 155 ) input_variables(:), 'north_bc', 'south_bc', 'east_bc', 'west_bc'
write( 15, 156 ) input_values(:), north_bc_input, south_bc_input, east_bc_input, west_bc_input
close(15)

155  format (  22( A, ',') )
156  format (  17( E12.7, ','), 4 ( A, ',') )
!-------------------------------------------------------------------------
!Write time history to a file
!  ( hydrographs and other time-dependant data
!    is plotted from this file )

OPEN( UNIT = 20, FILE = 'details.csv', STATUS='REPLACE')
WRITE(20,*) 'Timestamp,', FILE_DATE, ' ', FILE_TIME, ','
DO i = 1, 17 
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
                                    'h_imid_max_hist,', &
                                          'eta_sheet,' 
DO n = 1, nlast
WRITE(20,300) n, numit(n), maxdiff(n), loc(n)         , &
              L2_History(n), rain(n)*1000.*3600.      , &
              maxthk(n)*100., time(n), -Qout(n)*1000. , &
              h_imid_j1_hist(n), h_imid_max_hist(n)   , & 
              eta_sheet( n )
end do
close(20)    


!-------------------------------------------------------------------------
! Write other outputs appropriate for a 2D model

if( model_mode .eq. '2D') then

    CALL OUTPUTS_2D( )

    CALL WRITE_WATER_BUDGET( nlast )

end if




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


711 format( (I12, ','), (A, ','), (F8.2, ',') )

!----------------------------------------------------------------------------------
! Message for sucessful completion

call DATE_AND_TIME( timechar, datechar, zonechar, dtvals  )

write(*,*) " " 
write(*,*) "PERFCODE: Program completed sucessfully at  ", &
                      ! hour:min:sec.millisec (what a pain!)             
                      dtvals(5), ":", &
                      dtvals(6), ":", &
                      dtvals(7), ".", &
                      dtvals(8) 
write(*,*) "          Output files are in the current directory. "


!----------------------------------------------------------------------------------
!Format statements

2       FORMAT( I12, ',', 10000 ( E12.7, ',') )
10      FORMAT('   ', ( i3, '    '),  ( F10.3, '  ') , F10.6 ) 
111     FORMAT( f9.2 )
123     format( F9.2 )
200     FORMAT ( A, ( E12.7, ',') ) 
201     FORMAT ( A, ( I12, ',') ) 
300     FORMAT ( 2 (   I12, ','),   F12.7, ',' , &    ! n, numit, maxdif
                       I12, ',' ,   E12.7, ',' , &    ! loc, L2_History
                 2 ( F12.8, ','), ( F12.3, ',' ),&    ! rain, maxthk, time 
                 4 ( F12.8 ,',')                     )! Qout, h_imid_j1hist, h_imid_max_hist
400     FORMAT( 10000 ( I12, ',' ) )
401     FORMAT( (I12, ',') , 2( F12.7, ',' ) )
660     FORMAT(  2( I12, ','), 2( F12.7, ',') )
!------------------------------------------------------------------------------
end program PERFCODE
!===================================================================================
!       \\\\\\\\\\          E N D    P R O G R A M          //////////
!       //////////             P E R F C O D E              \\\\\\\\\\
!===================================================================================

