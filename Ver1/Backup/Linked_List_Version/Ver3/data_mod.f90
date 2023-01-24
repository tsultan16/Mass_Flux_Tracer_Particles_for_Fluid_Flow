MODULE data_mod
    USE tracertype_mod
	IMPLICIT NONE 
    
	INTEGER, ALLOCATABLE :: N_cell(:), N_cell2d(:,:)
	REAL, ALLOCATABLE :: flux(:) , rho(:), flux2d(:,:,:) , rho2d(:,:) 
	! define array of tracer pointers
    TYPE (tracer_ptr_1d), ALLOCATABLE :: tr1d(:)
    TYPE (tracer_ptr_2d), ALLOCATABLE :: tr2d(:)

    ! define array of head tracer pointers 
	! (each cell has it's own linked list of tracers
	! the "head tracer" in a given cell is the root element
	! of its linked list.)
    TYPE (tracer_ptr_1d), ALLOCATABLE :: cellHead1d(:)
	TYPE (tracer_ptr_2d), ALLOCATABLE :: cellHead2d(:,:)

END MODULE data_mod
