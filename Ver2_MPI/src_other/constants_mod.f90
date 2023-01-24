MODULE constants_mod
	IMPLICIT NONE	
    
	INTEGER, PARAMETER :: N = 100
    INTEGER, PARAMETER :: nx = 10
    INTEGER, PARAMETER :: ny = 10
    INTEGER, PARAMETER :: nz = 1
    INTEGER, PARAMETER :: nb = 1
    INTEGER, PARAMETER :: ndims = 2
    INTEGER, PARAMETER :: nt = 2
	
	INTEGER, PARAMETER :: nvars = 4, max_mpi_buffer_particles = 100
	INTEGER, PARAMETER :: nranks_x = 1
	INTEGER, PARAMETER :: nranks_y = 2
	INTEGER, PARAMETER :: nranks_z = 1 
	

END MODULE constants_mod
