:: (c) Copyright 2015 Bradley J Eck
::
::windows batch file for compiling Perfcode with gfortran 
:: Step 1 compile the modules and in the proper order  
gfortran -c ..\src\shared.f95    
gfortran -c ..\src\pfc2Dfuns.f95  
gfortran -c ..\src\utilities.f95 
gfortran -c ..\src\inputs.f95    
gfortran -c ..\src\outputs.f95   
gfortran -c ..\src\geom_funcs.f95
gfortran -c ..\src\ConvCoef.f95  
gfortran -c ..\src\GridGen.f95   
gfortran -c ..\src\solvers.f95   
gfortran -c ..\src\pfc1dfuns2.f95
gfortran -c ..\src\pfc1Dfuns.f95 
gfortran -c ..\src\pfc1Dsubs.f95 
gfortran -c ..\src\pfc2Dsubs.f95 
gfortran -c ..\src\BoundCond.f95 
gfortran -c ..\src\applyBCs.f95 
gfortran -c ..\src\budget.f95 

:: Step 2 compile and link the main program 
gfortran -o ..\bin\Perfcode ..\src\PERFCODE.f95 shared.o pfc2Dfuns.o utilities.o inputs.o outputs.o geom_funcs.o ConvCoef.o GridGen.o solvers.o pfc1dfuns2.o pfc1dfuns.o pfc1dsubs.o pfc2dsubs.o boundcond.o applyBCs.o budget.o
