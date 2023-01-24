#!/bin/sh
ifort -c constants_mod.f90
ifort -c tracerType_mod.f90
ifort -c data_mod.f90
ifort -c tracerSolver_mod.f90
ifort constants_mod.o tracerType_mod.o data_mod.o tracerSolver_mod.o tracer_driver.f90 -o tracer
rm *.o
rm *.mod
#./tracer
