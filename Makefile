# compile the executable and place into the bin directory
all: modules
	cd build; \
	gfortran -o ../bin/Perfcode ../src/PERFCODE.f95 shared.o pfc2Dfuns.o \
	Utilities.o inputs.o outputs.o geom_funcs.o ConvCoef.o GridGen.o \
	solvers.o pfc1Dfuns2.o pfc1Dfuns.o pfc1Dsubs.o pfc2Dsubs.o BoundCond.o \
	applyBCs.o budget.o

# compile all of the modules into the build directory 
modules: src/*.f95
	cd build; \
	gfortran -c  ../src/shared.f95     ;\
	gfortran -c  ../src/pfc2Dfuns.f95  ;\
	gfortran -c  ../src/Utilities.f95  ;\
	gfortran -c  ../src/inputs.f95     ;\
	gfortran -c  ../src/outputs.f95    ;\
	gfortran -c  ../src/geom_funcs.f95 ;\
	gfortran -c  ../src/ConvCoef.f95   ;\
	gfortran -c  ../src/GridGen.f95    ;\
	gfortran -c  ../src/solvers.f95    ;\
	gfortran -c  ../src/pfc1Dfuns2.f95 ;\
	gfortran -c  ../src/pfc1Dfuns.f95  ;\
	gfortran -c  ../src/pfc1Dsubs.f95  ;\
	gfortran -c  ../src/pfc2Dsubs.f95  ;\
	gfortran -c  ../src/BoundCond.f95  ;\
	gfortran -c  ../src/applyBCs.f95   ;\
	gfortran -c  ../src/budget.f95     


test: bin/Perfcode
	cd test-curved-road-2d; ../bin/Perfcode
	cd test-straight-road-2d; ../bin/Perfcode

.PHONY: clean cleanobj cleanmod
clean : cleanobj cleanmod
	rm -f bin/Perfcode

cleanobj :
	rm -f build/*.o

cleanmod : 
	rm -f build/*.mod
