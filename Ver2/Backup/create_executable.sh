#!/bin/sh
mpifort -c src_other/constants_mod.f90
mpifort -c src_main/tracerType_mod.f90
mpifort -c src_other/data_mod.f90
mpifort -c src_other/tracerInit_mod.f90
mpifort -c src_main/tracerSolver_mod.f90
mpifort -check bounds constants_mod.o tracerType_mod.o data_mod.o tracerInit_mod.o tracerSolver_mod.o src_other/tracer_driver.f90 -o tracer 
rm *.o
rm *.mod
#./tracer
