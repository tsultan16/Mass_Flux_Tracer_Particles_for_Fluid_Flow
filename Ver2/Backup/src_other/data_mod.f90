MODULE data_mod
    USE tracertype_mod
	IMPLICIT NONE 
    
	INTEGER, ALLOCATABLE :: N_cell_1d(:), N_cell_2d(:,:), N_cell_3d(:,:,:)
	REAL, ALLOCATABLE :: flux_1d(:)       , rho_1d(:),   &
    				     flux_2d(:,:,:)   , rho_2d(:,:), &
						 flux_3d(:,:,:,:) , rho_3d(:,:,:) 
 
                   

	! define array of tracer pointers
    TYPE (node_ptr), ALLOCATABLE :: tr(:)


    ! Define array of tree_head pointers. 
	! Each cell has it's own binary tree, 
    ! each node in the tree is represented by a tracer.
    TYPE (node_ptr), ALLOCATABLE :: cellHead_1d(:), cellHead_2d(:,:), cellHead_3d(:,:,:)  

    REAL :: randtime = 0.0, advecttime = 0.0

    REAL :: rate 
    INTEGER(8) :: cr,cm


END MODULE data_mod
