PROGRAM tracer_test_run

USE MPI
USE constants_mod
USE tracertype_mod
USE data_mod
USE tracerInit_mod
USE tracersolver_mod
USE domain_mod

IMPLICIT NONE

REAL :: dt1, dx1, dy1, dz1
INTEGER :: i, ix, iy, iz
INTEGER ::c_i,c_f
REAL*8 :: t1,t2,t_tot
LOGICAL :: all_particles_buffered
TYPE(MASSFLUXTRACER) :: solver


!##########################################################################################

! Initialize MPI
CALL MPI_INIT(ierr)

! Get number of processes
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs(1), ierr)

! set up domain decomposition
CALL setup_domain()

! get rank
CALL MPI_COMM_RANK(comm2d, myrank, ierr)

! get rank co-ordinates
CALL MPI_CART_COORDS(comm2d, myrank, ndims, mycoord, ierr)

! compute domain boundary indices
CALL compute_bound()

!allocate memory for MPI buffers
buffer_size = 1+nvars*max_mpi_buffer_particles
ALLOCATE(MPI_buffer_in(1:buffer_size),MPI_buffer_out(1:buffer_size))


PRINT*,''
PRINT*,'My rank, coordinate, xlow, xhi, ylow, yhi =',myrank, mycoord, xlow, xhi, ylow, yhi
PRINT*,''
!##########################################################################################



! Initialize the system_clock
!CALL system_clock(count_rate=cr)
!CALL system_clock(count_max=cm)
!rate = REAL(cr)
!WRITE(*,*) "system_clock rate = ",rate

!Initialize grid variables
CALL data_init()

!set up solver object
solver = MASSFLUXTRACER(ndims,nb,N,dx1,dy1,dz1)
CALL solver%initialize_workpool()

! initialize tracers
IF(ndims .EQ. 1) CALL initialize_tracer_distribution_1d(N_cell_1d, cellHead_1d, cellHead_bndry_out)

IF(ndims .EQ. 2) CALL initialize_tracer_distribution_2d(N_cell_2d, cellHead_2d, cellHead_bndry_out)

IF(ndims .EQ. 3) CALL initialize_tracer_distribution_3d(N_cell_3d, cellHead_3d, cellHead_bndry_out)


!CALL output1d(N_cell1d)

CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

t1 = MPI_Wtime()

! Simulation loop    
DO i= 1, nt
    PRINT*,' '
    PRINT*,'TIME STEP, % complete = ',i, (i*100./(nt*1.))
    PRINT*,' '

    ! advect tracers
    CALL solver%solve(dt1)

    !################################
    CALL boundary_particle_exchange()
    !################################	


    IF(ndims .EQ. 1) THEN
	    !CALL output1d(N_cell_1d)
	END IF  
  
    IF(ndims .EQ. 2) THEN
	    !CALL output2d(N_cell_2d,i)
	END IF

    IF(i .EQ. 1) flux_2d = 0.d0
	
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
 
END DO
t2 = MPI_Wtime()
 
t_tot = t2-t1
 
 
IF(myrank .EQ. 0) THEN
	PRINT*,'Total simulation time (sec) =',t_tot
	PRINT*,'Advect Time (sec) = ',advecttime
	PRINT*,'Rand num generation Time (sec) = ',randtime
	PRINT*,'Rand time fraction =', randtime/t_tot
	PRINT*,'Advect time fraction =', advecttime/t_tot
END IF

CALL solver%destroy_workpool()

IF(ndims .EQ. 1) DEALLOCATE(cellHead_1d, rho_1d, flux_1d, N_cell_1d)
IF(ndims .EQ. 2) DEALLOCATE(cellHead_2d, rho_2d, flux_2d, N_cell_2d)
IF(ndims .EQ. 3) DEALLOCATE(cellHead_3d, rho_3d, flux_3d, N_cell_3d)
DEALLOCATE(cellHead_bndry_in, cellHead_bndry_out, N_cell_bndry_in, N_cell_bndry_out)
DEALLOCATE(MPI_buffer_in, MPI_buffer_out)

CLOSE(UNIT=1)


PRINT*,'RANK#',myrank,' DONE!'


!Terminate MPI
CALL MPI_FINALIZE(ierr)


!##########################################################################################


CONTAINS



SUBROUTINE boundary_particle_exchange()

	INTEGER:: bndry, bndry_in
	
    ! loop over boundaries, pack outgoing particles, unpack incoming particles	
	DO bndry = 1, 2*ndims
		all_particles_buffered = .FALSE.
		DO WHILE(.NOT. all_particles_buffered)
		    
			! clear buffers
    		 MPI_buffer_out = 0.0
			 MPI_buffer_in = 0.0
			 
			! start packing
			CALL solver%pack_MPI_buffer_out(MPI_buffer_out, N_cell_bndry_out, cellHead_bndry_out, bndry, max_mpi_buffer_particles, all_particles_buffered)
			
            ! exchange particle buffer with neighbor domain
			CALL MPI_communication(MPI_buffer_out, MPI_buffer_in, buffer_size, bndry)
			
			SELECT CASE(bndry)
				CASE(1)
			        bndry_in = 2
				CASE(2)
					bndry_in = 1
				CASE(3)
					bndry_in = 4
				CASE(4)
					bndry_in = 3
			END SELECT		
			
			! start unpacking
			CALL solver%unpack_MPI_buffer_in(MPI_buffer_in, cellHead_bndry_in, N_cell_bndry_in, bndry_in)			
			
		END DO	
	END DO
    !PRINT*,'Particle boundary packing and unpacking complete'
 


END SUBROUTINE boundary_particle_exchange



SUBROUTINE MPI_communication(buffer_out, buffer_in, buffer_size, bndry)

    REAL, INTENT(INOUT) :: buffer_out(:), buffer_in(:)
    INTEGER, INTENT(IN) :: buffer_size, bndry
	
    INTEGER :: bottom, top, left, right, neighbor_rank
	INTEGER :: i, tag
	
    ! get neighbor ranks
    CALL MPI_CART_SHIFT(comm2d, 0, 1, left, right, ierr) 
    CALL MPI_CART_SHIFT(comm2d, 1, 1, bottom, top, ierr)
     
    SELECT CASE(bndry)

    CASE(1) ! x+
        neighbor_rank = right
        PRINT*,'Exchanging with right neighbor, Rank=',neighbor_rank 
    CASE(2) ! x-
        neighbor_rank = left
        PRINT*,'Exchanging with left neighbor, Rank=',neighbor_rank
    CASE(3) ! y+
        neighbor_rank = top
        PRINT*,'Exchanging with top neighbor, Rank=',neighbor_rank
    CASE(4) ! y-
        neighbor_rank = bottom
        PRINT*,'Exchanging with bottom neighbor, Rank=',neighbor_rank
    END SELECT
	
	PRINT*,'Particles in out buffer = ',buffer_out(1)
       
    ! send and receive particle buffers from neighbor
	tag = 0
    CALL MPI_SENDRECV(buffer_out,            &  ! send buffer
                      buffer_size,           &  ! send count
                      MPI_REAL,              &  ! data type
                      neighbor_rank,         &  ! dest
                      tag,                   &  ! send tag
                      buffer_in,             &  ! recv buffer
                      buffer_size,           &  ! recv count
                      MPI_REAL,              &  ! data type
                      neighbor_rank,         &  ! source
                      tag,                   &  ! recv tag
                      comm2d, MPI_STATUS_IGNORE, ierr)

    !PRINT*,'Communication complete'

END SUBROUTINE MPI_communication
	

SUBROUTINE data_init()


IF(ndims .EQ. 1) THEN
    OPEN(UNIT=1, FILE='output1d.txt')	
   
    ALLOCATE(cellHead_1d(xlow:xhi))
	ALLOCATE(cellHead_bndry_in(1:2),cellHead_bndry_out(1:2))
    ALLOCATE(rho_1d(xlow:xhi), flux_1d(xlow:xhi-1))
    ALLOCATE(N_cell_1d(xlow:xhi))
	ALLOCATE(N_cell_bndry_in(1:2),N_cell_bndry_out(1:2))
    rho_1d = 0.0
    flux_1d = 0.0
    N_cell_1d = 0
    N_cell_bndry_in = 0
	N_cell_bndry_out = 0
	
ELSE IF(ndims .EQ. 2) THEN
    OPEN(UNIT=1, FILE='output1d.txt')	
   
    ALLOCATE(cellHead_2d(xlow:xhi,ylow:yhi))
	ALLOCATE(cellHead_bndry_in(1:4),cellHead_bndry_out(1:4))
    ALLOCATE(rho_2d(xlow:xhi,ylow:yhi), flux_2d(xlow:xhi-1,ylow:yhi-1,2))
    ALLOCATE(N_cell_2d(xlow:xhi,ylow:yhi))
  	ALLOCATE(N_cell_bndry_in(1:4),N_cell_bndry_out(1:4))

    rho_2d = 0.0
    flux_2d = 0.0
    N_cell_2d = 0
    N_cell_bndry_in = 0
	N_cell_bndry_out = 0
	
ELSE IF(ndims .EQ. 3) THEN
    OPEN(UNIT=1, FILE='output1d.txt')	
   
    ALLOCATE(cellHead_3d(xlow:xhi,ylow:yhi,zlow:zhi))
	ALLOCATE(cellHead_bndry_in(1:6),cellHead_bndry_out(1:6))
    ALLOCATE(rho_3d(xlow:xhi,ylow:yhi,zlow:zhi), flux_3d(xlow:xhi-1,ylow:yhi-1,zlow:zhi-1,3))
    ALLOCATE(N_cell_3d(xlow:xhi,ylow:yhi,zlow:zhi))
  	ALLOCATE(N_cell_bndry_in(1:6),N_cell_bndry_out(1:6))

    rho_3d = 0.0
    flux_3d = 0.0
    N_cell_3d = 0
    N_cell_bndry_in = 0
	N_cell_bndry_out = 0

ELSE
    PRINT*,'NEED TO HAVE 1 <= NDIMS <= 3'
	STOP

END IF  
  
! allocate root nodes for boundary cells
DO i= 1,ndims*2
	ALLOCATE(cellHead_bndry_in(i)%p)
END DO
	
	
dt1 = 0.1/5.0
dx1 = 1./real(nx)
dy1 = 1./real(ny)
dz1 = 1./real(nz)

IF(ndims .EQ. 1) THEN
	rho_1d(:) = 1.
	
	DO ix = xlow, xhi-1 
	   	flux_1d(ix) = 0.2501*dx1/dt1
	END DO
		
END IF

IF(ndims .EQ. 2) THEN
	rho_2d(:,:) = 1.
	
	DO ix = xlow, xhi-1
		DO iy = ylow, yhi-1 
		    !IF(ix < nb+nx/2) THEN
			!	flux_2d(ix,iy,1) = -0.2501*dx1*dy1/dt1
			!ELSE
			!	flux_2d(ix,iy,1) = 0.2501*dx1*dy1/dt1
			!END IF
			flux_2d(ix,iy,1) = 0.d0!0.2501*dx1*dy1/dt1
			flux_2d(ix,iy,2) = 0.2501*dx1*dy1/dt1
		END DO
	END DO
END IF


END SUBROUTINE data_init

SUBROUTINE output1d(N_cell)

    INTEGER, INTENT(INOUT) :: N_cell(:)
    INTEGER :: i
    REAL :: x

    DO i=1,nx
        x=i*dx1
        WRITE(1,*) x,N_Cell(i)
    END DO

END SUBROUTINE output1d


SUBROUTINE output2d(N_cell,ts)

	INTEGER, INTENT(IN) :: N_cell(1-nb:nx+nb,1-nb:ny+nb), ts
    INTEGER :: i,j
	REAL :: x,y
    CHARACTER(len=40)::filename
    CHARACTER(len=6)::uniti


	IF(ts < 10) THEN
		WRITE(uniti,'(I1.1)') ts
	ELSE IF(ts>=10 .and. ts<100) THEN
		WRITE(uniti,'(I2.2)') ts
	ELSE IF(ts>=100 .and. ts<1000) THEN
		WRITE(uniti,'(I3.3)') ts
	ELSE IF(ts>=1000 .and. ts<10000) THEN
		WRITE(uniti,'(I4.3)') ts
	ELSE IF(ts>=10000 .and. ts<100000) THEN
		WRITE(uniti,'(I5.3)') ts
	END IF
  
filename=trim('Output/t=')//trim(uniti)//trim('.txt')
!print*,'filename=',filename

OPEN(unit=12,file=filename)

DO j = 1-nb, ny+nb
    y = (j-1) * dy1
    DO i = 1-nb , nx+nb
       x = (i-1) * dx1
       WRITE(12,*) x,y,REAL(N_cell(i,j))
  END DO
END DO

CLOSE(12)

END SUBROUTINE output2d


		
END PROGRAM tracer_test_run

