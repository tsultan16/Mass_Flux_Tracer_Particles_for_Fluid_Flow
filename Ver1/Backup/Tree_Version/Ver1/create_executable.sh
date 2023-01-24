#!/bin/sh
ifort -c src_other/constants_mod.f90
ifort -c src_main/tracerType_mod.f90
ifort -c src_other/data_mod.f90
ifort -c src_main/tracerSolver_mod.f90
ifort -check bounds constants_mod.o tracerType_mod.o data_mod.o tracerSolver_mod.o src_other/tracer_driver.f90 -o tracer 
rm *.o
rm *.mod
#./tracer
