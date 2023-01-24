MODULE tracertype_mod
IMPLICIT NONE


! tracer data structures for 1,2 and 3D


TYPE :: tracer_1d ! need to get rid of the '_1d' label 

    ! type data
    INTEGER :: id = 0                            ! label for distinguishing tracers
    INTEGER :: x = 0                             ! cell location indices
    INTEGER :: y = 0   
    INTEGER :: z = 0   
    REAL :: some_fluid_property = 0.0            ! fluid property recorded by tracer at it's current location

END TYPE tracer_1d


! binary tree node derived data type

TYPE, EXTENDS(tracer_1d) :: node
    ! type data
    INTEGER :: bf = 0                          ! balance factor(= height of right sub-tree - height of left sub-tree)
    INTEGER :: rank = 1                        ! rank (i.e. 1+ #of nodes in left sub-tree = # of nodes with smaller key)
                                               ! can be used to determine relative position of a node within a sub-tree        
         
    TYPE(node), POINTER :: node_L => null()    ! pointer to left sub-tree/child
    TYPE(node), POINTER :: node_R => null()    ! pointer to right sub-tree/child

END TYPE node



! auxiliary path pointer
TYPE path_ptr
    INTEGER :: a = 0
    TYPE(node), POINTER :: node
    TYPE(path_ptr), POINTER :: next => null()
    TYPE(path_ptr), POINTER :: prev => null()
END TYPE path_ptr


! derived type for making pointer arrays

TYPE node_ptr
    TYPE(node), POINTER :: p
END TYPE node_ptr




END MODULE tracertype_mod