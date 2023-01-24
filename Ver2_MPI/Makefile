OBJS = constants_mod.o tracerType_mod.o data_mod.o tracerInit_mod.o \
       domain_mod.o tracerSolver_mod.o tracer_driver.o 

tracer : $(OBJS)
	mpifort -check bounds -o tracer $(OBJS)
	
constants_mod.o : src_other/constants_mod.f90
	mpifort -c src_other/constants_mod.f90
	
tracerType_mod.o : src_main/tracerType_mod.f90
	mpifort -c src_main/tracerType_mod.f90
	
data_mod.o : src_other/data_mod.f90 tracerType_mod.o
	mpifort -c src_other/data_mod.f90
	
domain_mod.o : src_other/domain_mod.f90 constants_mod.o data_mod.o
	mpifort -c src_other/domain_mod.f90

tracerInit_mod.o : src_other/tracerInit_mod.f90 constants_mod.o tracerType_mod.o data_mod.o
	mpifort -c src_other/tracerInit_mod.f90
	
tracerSolver_mod.o : src_main/tracerSolver_mod.f90  data_mod.o tracerType_mod.o
	mpifort -c src_main/tracerSolver_mod.f90
	
tracer_driver.o : src_other/tracer_driver.f90  constants_mod.o tracerType_mod.o data_mod.o domain_mod.o tracerInit_mod.o tracerSolver_mod.o
	mpifort -c src_other/tracer_driver.f90
	
clean:
	rm  -f *.o *.mod 