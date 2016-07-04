all: modules
	gfortran -o bin/Perfcode src/PERFCODE.f95 build/shared.o build/pfc2Dfuns.o build/Utilities.o build/inputs.o build/outputs.o build/geom_funcs.o build/ConvCoef.o build/GridGen.o build/solvers.o build/pfc1Dfuns2.o build/pfc1Dfuns.o build/pfc1Dsubs.o build/pfc2Dsubs.o build/BoundCond.o build/applyBCs.o build/budget.o

modules: src/*.f95
	gfortran -c -o build/shared.o src/shared.f95
	gfortran -c -o build/pfc2Dfuns.o src/pfc2Dfuns.f95
	gfortran -c -o build/Utilities.o src/Utilities.f95
	gfortran -c -o build/inputs.o src/inputs.f95
	gfortran -c -o build/outputs.o src/outputs.f95
	gfortran -c -o build/geom_funcs.o src/geom_funcs.f95
	gfortran -c -o build/ConvCoef.o src/ConvCoef.f95
	gfortran -c -o build/GridGen.o src/GridGen.f95
	gfortran -c -o build/solvers.o src/solvers.f95
	gfortran -c -o build/pfc1Dfuns2.o src/pfc1Dfuns2.f95
	gfortran -c -o build/pfc1Dfuns.o src/pfc1Dfuns.f95
	gfortran -c -o build/pfc1Dsubs.o src/pfc1Dsubs.f95
	gfortran -c -o build/pfc2Dsubs.o src/pfc2Dsubs.f95
	gfortran -c -o build/BoundCond.o src/BoundCond.f95
	gfortran -c -o build/applyBCs.o src/applyBCs.f95
	gfortran -c -o build/budget.o src/budget.f95

clean:
	rm bin/Perfcode  
	rm build/*.o

test: bin/Perfcode
	cd test-curved-road-2d; ../bin/Perfcode
	cd test-straight-road-2d; ../bin/Perfcode
