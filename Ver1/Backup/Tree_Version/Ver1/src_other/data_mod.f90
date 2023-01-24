MODULE data_mod
    USE tracertype_mod
	IMPLICIT NONE 
    
	INTEGER, ALLOCATABLE :: N_cell1d(:), N_cell2d(:,:), N_cell3d(:,:,:)
	REAL, ALLOCATABLE :: flux(:) , rho(:), flux2d(:,:,:) , rho2d(:,:), &
                         flux3d(:,:,:,:) , rho3d(:,:,:) 

	! define array of tracer pointers
    TYPE (tracer_ptr_1d), ALLOCATABLE :: tr1d(:)
    TYPE (tracer_ptr_2d), ALLOCATABLE :: tr2d(:)
    TYPE (tracer_ptr_3d), ALLOCATABLE :: tr3d(:)


    ! define array of head tracer pointers 
	! (each cell has it's own linked list of tracers
	! the "head tracer" in a given cell is the root element
	! of its linked list.)
    TYPE (tree_root_1d), ALLOCATABLE :: cellRoot1d(:)
	!TYPE (tracer_ptr_2d), ALLOCATABLE :: cellHead2d(:,:)
	!TYPE (tracer_ptr_3d), ALLOCATABLE :: cellHead3d(:,:,:)

END MODULE data_mod
