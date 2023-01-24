PROGRAM tracer_test_run
    USE constants_mod
	USE tracertype_mod
    USE data_mod
	USE tracersolver_mod
	IMPLICIT NONE


	! define mass flux tracer object
    TYPE(MASSFLUXTRACER) :: solver
    REAL :: dt1, dx1, dy1
    INTEGER :: i

    OPEN(UNIT=1, FILE='output1d.txt')	


	IF(ndims .EQ. 1) THEN
		ALLOCATE(tr1d(1:N))
		ALLOCATE(cellHead1d(0:nx+1))
		ALLOCATE(rho(1:nx), flux(0:nx))
		ALLOCATE(N_cell(0:nx+1))
  
    ELSE IF(ndims .EQ. 2) THEN
		ALLOCATE(tr2d(1:N))
		ALLOCATE(cellHead2d(0:nx+1,0:ny+1))
		ALLOCATE(rho2d(1:nx,1:ny), flux2d(0:nx,0:ny,2))
		ALLOCATE(N_cell2d(0:nx+1,0:ny+1))
    END IF


    N_cell = 0

	! set up solver object
    solver = MASSFLUXTRACER(ndims,nx,ny,nz,N)
    CALL solver%initialize_workpool()


    ! initialize tracers
	IF(ndims .EQ. 1) THEN
		CALL initialize_tracer_distribution_1d(tr1d, cellHead1d, N_cell)
	ELSE IF(ndims .EQ. 2) THEN
		CALL initialize_tracer_distribution_2d(tr2d, cellHead2d, N_cell2d)
	END IF

    ! advance by one advection step using sample rho and flux arrays
    dt1 = 0.1/5.
    dx1 = 1./real(nx)
    dy1 = 1./real(ny) 

	IF(ndims .EQ. 1) THEN
		rho = 1.
		flux(0:nx/2) = -1./4.
        flux(1+nx/2:nx) = 1./4.
    ELSE IF(ndims .EQ. 2) THEN
		rho2d = 1.0
		flux2d(:,:,1) = 5.0
        flux2d(:,:,2) = 5.0
	END IF


    ! Simulation loop    
    DO i= 1, nt
        PRINT*,'TIME STEP, % complete = ',i, (i*100./(nt*1.))
        CALL solver%solve(dt1, dx1, dy1)
        CALL fileOut1d()	

    END DO
    



    CALL solver%destroy_workpool()



	IF(ndims .EQ. 1) THEN
		DEALLOCATE(tr1d,cellHead1d,rho,flux,N_cell)
  
    ELSE IF(ndims .EQ. 2) THEN
		DEALLOCATE(tr2d,cellHead2d,rho2d,flux2d,N_cell2d)
    END IF

    CLOSE(UNIT=1)

    PRINT*,'DONE!'




	CONTAINS

		! This subroutine assigns a cell to each tracer 
        ! (according to a given distribution)
        ! For now, tracers will be randomly distributed
        ! across the cells.
		SUBROUTINE initialize_tracer_distribution_1d(tr, cellHead1d, N_cell)

			TYPE(tracer_ptr_1d), INTENT(INOUT) :: tr(:),cellHead1d(:)
            INTEGER, INTENT(INOUT) :: N_cell(:)

			TYPE(tracer_1d), POINTER :: current
     		INTEGER :: i, cellnum
            REAL :: p 
        

            ! loop over tracers, assign them to cells
            DO i = 1, N

                ! allocate memory for tracer pointer
                ALLOCATE(tr(i)%p) 

            	! draw a random integer between nx/3 and nx*2/3
                CALL RANDOM_NUMBER(p)
                cellnum= nx/3+FLOOR(p*nx/3)

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

		END SUBROUTINE initialize_tracer_distribution_1d


		SUBROUTINE initialize_tracer_distribution_2d(tr, cellHead, N_cell)

			TYPE(tracer_ptr_2d), INTENT(INOUT) :: tr(:),cellHead(:,:)
            INTEGER, INTENT(INOUT) :: N_cell(:,:)

			TYPE(tracer_2d), POINTER :: current
     		INTEGER :: i, j, cellnumx, cellnumy
            REAL :: p 
        

            ! loop over tracers, assign them to cells
            DO i = 1, N

                ! allocate memory for tracer pointer
                ALLOCATE(tr(i)%p) 

            	! draw a random integer between 1 and nx
                CALL RANDOM_NUMBER(p)
                cellnumx= 1+FLOOR(p*nx)
                CALL RANDOM_NUMBER(p)
				cellnumy= 1+FLOOR(p*ny)

				! assign id and cell number to tracer 
                tr(i)%p%id = i
                tr(i)%p%x = cellnumx
				tr(i)%p%y = cellnumy
                N_cell(cellnumx+1,cellnumy+1) = N_cell(cellnumx+1,cellnumy+1)+1

                ! if cell is empty, designate this tracer as head
                IF(.NOT. ASSOCIATED(cellHead(cellnumx+1,cellnumy+1)%p)) THEN
                   cellHead(cellnumx+1,cellnumy+1)%p => tr(i)%p

				! otherwise add current tracer to the cell's linked list                
				ELSE
	              	current => cellHead(cellnumx+1,cellnumy+1)%p

                     
                    DO WHILE(ASSOCIATED(current%next))
                        ! traverse through list until last tracer is reached
						current => current%next			
                    END DO
                    current%next => tr(i)%p
				END IF				     

			END DO


        	PRINT*, 'Tracer initialization completed.'

		END SUBROUTINE initialize_tracer_distribution_2d


        SUBROUTINE fileout1d()

            REAL::x
            INTEGER::i

                DO i=1,nx   
                    WRITE(1,*) i*dx, N_cell(i) 
                END DO

        END SUBROUTINE fileout1d


END PROGRAM tracer_test_run

