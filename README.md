# Perfcode : PERmeable Friction COurse Drainage codE

## Overview 
Perfcode is a numerical model for simulating drainage from roadways. The model
was developed to study the drainage behaviour of a porous asphalt known as
permeable friction course that is placed on top of regular impermeable
pavement.   Under small rainfall intensities, drainage is contained within the
porous layer; but under higher rainfall intensities, drainage occurs both
within and on top of the porous pavement.  Perfcode simulates this coupled
unsteady drainage process. 

The model predicts the variation of water depth within and on top of the porous
layer through time. Inputs are the rainfall pattern (hyetograph), geometric
information regarding the roadway layout, and hydraulic properties of the
pavement.  The porous layer is treated as an unconfined aquifer using Darcy’s
law and the Dupuit-Forchheimer assumptions. Surface flow is modeled using the
diffusion wave approximation to the shallow water equations. A mass balance
approach is used to couple surface and subsurface phases. Straight and curved
roadway geometries are accommodated via a curvilinear grid.

## References
Further details on the development, validation, and application of Perfcode are
available in the following references. Please cite these references when using
Perfcode.  

BJ Eck, ME Barrett, RJ Charbeneau (2012)
[Coupled Surface-Subsurface Model for Simulating Drainage
 from Permeable Friction Course Highways]
(http://ascelibrary.org/doi/abs/10.1061/%28ASCE%29HY.1943-7900.0000474)
Journal of Hydraulic Engineering 138 (1), 13-22

BJ Eck, RJ Charbeneau, ME Barrett (2010)
[Drainage hydraulics of porous pavement: Coupling surface and subsurface flow]
(http://www.crwr.utexas.edu/reports/2010/rpt10-2.shtml)
Center for Research in Water Resources, University of Texas at Austin

A video showing results of a Perfcode simulation is also available at http://bradeck.net/pfc-video/index.htm 

## Quick Start for Windows
The easiest way to become familiar with the inputs Perfcode needs and the
outputs it generates is to run one of the test examples included in the
repository. 

1. Download the zip file of the repository and Unzip to your preferred location
2. Open a command prompt in one of the \test- directories 
3. Invoke the run-test.bat batch file   
4. Refer to Appendix A of the report linked above for details on the input and output files.

## Compilation 
### Windows
A windows binary is included in the \bin directory. 

The application may also be compiled as follows:

1. Open a command prompt in the \build directory 
2. Confirm that gfortran is installed and can be invoked from the \build
   directory with command gfortran
3. Run the batch file compile.bat
4. Find the executable in the \bin directory   

### Linux 
ToDo



