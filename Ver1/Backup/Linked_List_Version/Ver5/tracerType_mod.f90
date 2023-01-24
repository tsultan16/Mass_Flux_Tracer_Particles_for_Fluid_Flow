MODULE tracertype_mod
	IMPLICIT NONE

	! tracer (derived)  types for 1,2 and 3D
	TYPE :: tracer_1d
        ! type data
        INTEGER :: id = 0							 ! label for distinguishing tracers
		INTEGER :: x = 0				 		 	 ! cell location indices
		REAL :: some_fluid_property = 0.0     		 ! fluid property recorded by tracer at it's current location
		TYPE (tracer_1d), POINTER :: next => null()  ! next pointer for linked list capability
	END TYPE tracer_1d

	TYPE :: tracer_2d
        ! type data
        INTEGER :: id = 0							 ! label for distinguishing tracers
		INTEGER :: x = 0				 		 	 ! cell location indices
        INTEGER :: y = 0
		REAL :: some_fluid_property = 0.0     		 ! fluid property recorded by tracer at it's current location
		TYPE (tracer_2d), POINTER :: next => null()  ! next pointer for linked list capability
	END TYPE tracer_2d

	TYPE :: tracer_3d
        ! type data
        INTEGER :: id = 0							 ! label for distinguishing tracers
		INTEGER :: x = 0				 		 	 ! cell location indices
        INTEGER :: y = 0
		INTEGER :: z = 0
		REAL :: some_fluid_property = 0.0     		 ! fluid property recorded by tracer at it's current location
		TYPE (tracer_3d), POINTER :: next => null()  ! next pointer for linked list capability
	END TYPE tracer_3d


	TYPE tracer_ptr_1d
		TYPE (tracer_1d), POINTER :: p => null()
    END TYPE tracer_ptr_1d

	TYPE tracer_ptr_2d
		TYPE (tracer_2d), POINTER :: p => null()
    END TYPE tracer_ptr_2d

    TYPE tracer_ptr_3d
		TYPE (tracer_3d), POINTER :: p => null()
    END TYPE tracer_ptr_3d


END MODULE tracertype_mod
