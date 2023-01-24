MODULE data_mod
    USE MPI
    USE tracertype_mod
	IMPLICIT NONE 
    
	! Grid variables
	INTEGER, ALLOCATABLE :: N_cell_1d(:), N_cell_2d(:,:), N_cell_3d(:,:,:)
	REAL, ALLOCATABLE :: flux_1d(:)       , rho_1d(:),   &
    				     flux_2d(:,:,:)   , rho_2d(:,:), &
						 flux_3d(:,:,:,:) , rho_3d(:,:,:) 
                    
	TYPE (node_ptr), ALLOCATABLE :: cellHead_1d(:), cellHead_2d(:,:), cellHead_3d(:,:,:), cellHead_bndry(:)  ! array of root-nodes
	TYPE (node_ptr), ALLOCATABLE :: tr(:) ! array of tracer pointers

    INTEGER :: xlow, xhi, ylow, yhi, zlow, zhi
    INTEGER :: my_particle_num

    ! MPI variables
	INTEGER :: comm2d, ierr, ndim, myrank, numprocs(1), dims(2), mycoord(2), req(4), &
           coltype, status(MPI_STATUS_SIZE), source, tag, destination, buffer_size
	LOGICAL :: isperiodic(2), reorder 
    REAL, ALLOCATABLE :: MPI_buffer_in(:), MPI_buffer_out(:)

    ! Timer variables
    REAL :: randtime = 0.0, advecttime = 0.0
    REAL :: rate 
    INTEGER(8) :: cr,cm


END MODULE data_mod
