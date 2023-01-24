MODULE tracertype_mod
IMPLICIT NONE

! tracer data structures for 1,2 and 3D

TYPE :: tracer_1d

    ! type data
    INTEGER :: id = 0                            ! label for distinguishing tracers
    INTEGER :: x = 0                             ! cell location indices
    REAL :: some_fluid_property = 0.0            ! fluid property recorded by tracer at it's current location
    TYPE (tracer_1d), POINTER :: next => null()  ! next pointer for linked list capability

END TYPE tracer_1d

TYPE :: tracer_2d

    ! type data
    INTEGER :: id = 0                            ! label for distinguishing tracers
    INTEGER :: x = 0                             ! cell location indices
    INTEGER :: y = 0
    REAL :: some_fluid_property = 0.0            ! fluid property recorded by tracer at it's current location
    TYPE (tracer_2d), POINTER :: next => null()  ! next pointer for linked list capability

END TYPE tracer_2d

TYPE :: tracer_3d

    ! type data
    INTEGER :: id = 0                            ! label for distinguishing tracers
    INTEGER :: x = 0                             ! cell location indices
    INTEGER :: y = 0
    INTEGER :: z = 0
    REAL :: some_fluid_property = 0.0            ! fluid property recorded by tracer at it's current location
    TYPE (tracer_3d), POINTER :: next => null()  ! next pointer for linked list capability

END TYPE tracer_3d



! "tree" data structure for 1D

TYPE tree_node_1d

    !INTEGER :: occupancy
    TYPE(tracer_1d), POINTER :: leaf_L => null()       ! pointer to leftmost leaf/tracer in the list 
    TYPE(tracer_1d), POINTER :: leaf_R => null()       ! pointer to rightmost leaf/tracer in the list 
    TYPE(tree_node_1d), POINTER :: next => null()      ! pointer to the next node on the right

END TYPE tree_node_1d

TYPE tree_root_1d

   INTEGER, ALLOCATABLE :: node_occupancy(:)          ! array for storing tracer occupancy within each node
   TYPE(tree_node_1d), POINTER :: node_L => null()    ! pointer to leftmost node
   TYPE(tree_node_1d), POINTER :: node_R => null()    ! pointer to rightmost node      

END TYPE tree_root_1d



! some derived types for making pointer arrays

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