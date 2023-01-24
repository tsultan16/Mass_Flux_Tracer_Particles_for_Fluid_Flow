MODULE domain_mod

USE MPI
USE constants_mod
USE data_mod

IMPLICIT NONE


CONTAINS


! domain decomposition setup
SUBROUTINE setup_domain()

    ! set cart communicator properties
	isperiodic(1) = .FALSE. ! periodic boundaries
	isperiodic(2) = .FALSE. ! periodic boundaries
	reorder = .TRUE. ! allow MPI rank reordering 
	dims(1) = nranks_x 
	dims(2) = nranks_y 
	
    ! create cart communicator handle
	CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, isperiodic, reorder, comm2d, ierr)


END SUBROUTINE setup_domain



SUBROUTINE compute_bound()

    xlow = 0 + mycoord(1)*(nx+2*nb)/dims(1)
    xhi = xlow +  (nx+2*nb)/dims(1) - 1
    ylow = 0 + mycoord(2)*(ny+2*nb)/dims(2)
    yhi = ylow + (ny+2*nb)/dims(2) - 1

END SUBROUTINE compute_bound








END MODULE domain_mod
