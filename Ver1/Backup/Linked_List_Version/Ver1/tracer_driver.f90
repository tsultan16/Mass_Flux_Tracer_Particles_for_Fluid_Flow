PROGRAM tracer_test_run
    USE constants_mod
    USE data_mod
	USE tracertype_mod
	USE tracersolver_mod
	IMPLICIT NONE

	! define array of tracer pointers
    TYPE (tracer_ptr_1d), ALLOCATABLE :: tr(:)

    ! define array of head tracer pointers 
	! (each cell has it's own linked list of tracers
	! the "head tracer" in a given cell is the root element
	! of its linked list.)
    TYPE (tracer_ptr_1d), ALLOCATABLE :: cellHead1d(:)

	! define mass flux tracer object
    TYPE(MASSFLUXTRACER) :: solver


    REAL :: dt1, dx1

	! allocate memory for tracer array
    ALLOCATE(tr(1:N))
	! allocate memory for cell head array
	ALLOCATE(cellHead1d(0:nx+1))

    ALLOCATE(rho(1:nx), flux(0:nx))
    ALLOCATE(N_cell(0:nx+1))
  

    N_cell = 0

	! set up solver object
    solver = MASSFLUXTRACER(ndims,nx,ny,nz,N)
    CALL solver%initialize_workpool()


    ! initialize tracers
    CALL initialize_tracer_distribution(tr, cellHead1d, N_cell)
   
    ! advance by one advection step using sample rho and flux arrays
    dt1 = 0.1/2.
    dx1 = 1./real(nx) 
    rho = 1.
    flux = (/ -1., -1., -1., 1., 1., 1. /)
    
    
    CALL solver%solve(tr, cellHead1d, dt1, dx1)
	

    CALL solver%destroy_workpool()

    PRINT*,'DONE!'


	CONTAINS

		! This subroutine assigns a cell to each tracer 
        ! (according to a given distribution)
        ! For now, tracers will be randomly distributed
        ! across the cells.
		SUBROUTINE initialize_tracer_distribution(tr, cellHead1d, N_cell)

			TYPE(tracer_ptr_1d), INTENT(INOUT) :: tr(:),cellHead1d(:)
            INTEGER, INTENT(INOUT) :: N_cell(:)

			TYPE(tracer_1d), POINTER :: current
     		INTEGER :: i, cellnum
            REAL :: p 
        

            ! loop over tracers, assign them to cells
            DO i = 1, N

                ! allocate memory for tracer pointer
                ALLOCATE(tr(i)%p) 

            	! draw a random integer between 1 and nx
                CALL RANDOM_NUMBER(p)
                cellnum= 1+FLOOR(p*nx)

				! assign id and cell number to tracer 
                tr(i)%p%id = i
                tr(i)%p%x = cellnum
                N_cell(cellnum+1) = N_cell(cellnum+1)+1

                ! if cell is empty, designate this tracer as head
                IF(.NOT. ASSOCIATED(cellHead1d(cellnum+1)%p)) THEN
                   cellHead1d(cellnum+1)%p => tr(i)%p

				! otherwise add current tracer to the cell's linked list                
				ELSE
	              	current => cellHead1d(cellnum+1)%p

                     
                    DO WHILE(ASSOCIATED(current%next))
                        ! traverse through list until last tracer is reached
						current => current%next			
                    END DO
                    current%next => tr(i)%p
				END IF				     

			END DO


        	PRINT*, 'Tracer initialization completed.'

		END SUBROUTINE initialize_tracer_distribution


END PROGRAM tracer_test_run

