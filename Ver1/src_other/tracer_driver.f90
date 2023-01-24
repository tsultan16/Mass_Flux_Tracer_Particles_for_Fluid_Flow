PROGRAM tracer_test_run
    USE constants_mod
	USE tracertype_mod
    USE data_mod
	USE tracersolver_mod
    USE MPI

	IMPLICIT NONE

	! define mass flux tracer object
    TYPE(MASSFLUXTRACER) :: solver
    REAL :: dt1, dx1, dy1, dz1
    INTEGER :: i

    REAL :: mfrac(4)
    REAL :: time_start, time_end, time_init, time_sim(nt), time_del, t1, t2, t_tot


    IF(ndims .EQ. 1) THEN
        OPEN(UNIT=1, FILE='output1d.txt')	
    ELSE IF(ndims .EQ. 2) THEN
        OPEN(UNIT=1, FILE='output2d.txt')
    END IF

	IF(ndims .EQ. 1) THEN
		ALLOCATE(tr1d(1:N))
		ALLOCATE(cellRoot1d(0:nx+1))
		ALLOCATE(rho(1:nx), flux(0:nx))
		ALLOCATE(N_cell1d(0:nx+1))
  
        rho = 0.0
        flux = 0.0
        N_cell1d = 0

    ELSE IF(ndims .EQ. 2) THEN
		!ALLOCATE(tr2d(1:N))
		!ALLOCATE(cellHead2d(0:nx+1,0:ny+1))
		!ALLOCATE(rho2d(1:nx,1:ny), flux2d(0:nx,0:ny,2))
		!ALLOCATE(N_cell2d(0:nx+1,0:ny+1))

        !rho2d = 0.0
        !flux2d = 0.0
        !N_cell2d = 0

    ELSE IF(ndims .EQ. 3) THEN
		!ALLOCATE(tr3d(1:N))
		!ALLOCATE(cellHead3d(0:nx+1,0:ny+1,0:nz+1))
		!ALLOCATE(rho3d(1:nx,1:ny,1:nz), flux3d(0:nx,0:ny,0:nz,3))
		!ALLOCATE(N_cell3d(0:nx+1,0:ny+1,0:nz+1))
   
        !rho3d = 0.0
        !flux3d = 0.0
        !N_cell3d = 0
    END IF

    dt1 = 0.1/5.0
    dx1 = 1./real(nx)
    dy1 = 1./real(ny) 
    dz1 = 1./real(ny) 

	IF(ndims .EQ. 1) THEN
		rho = 1.
		!flux(0:nx/2) = -0.251*dx1/dt1
        !flux(1+nx/2:nx) = 0.251*dx1/dt1
        flux(:) = 0.51*dx1/dt1
        
    ELSE IF(ndims .EQ. 2) THEN
		!rho2d = 1.0
		!flux2d(0:nx/2,:,1) = -0.25*dx1/dt1
		!flux2d(1+nx/2:nx,:,1) = 0.25*dx1/dt1
        !flux2d(:,0:ny/2,2) = -0.5*dx1/dt1
        !flux2d(:,1+ny/2:ny,2) = 0.5*dx1/dt1
        !flux2d(:,:,1)=2.*dx1/dt1
        !flux2d(:,:,2)=0.0

    ELSE IF(ndims .EQ. 3) THEN
		!rho3d = 1.0
		!flux3d(:,:,:,1) = 0.0
        !flux3d(:,:,:,2) = 0.0
        !flux3d(:,:,:,3) = 100.0
	END IF
  

	! set up solver object
    solver = MASSFLUXTRACER(ndims,nx,ny,nz,N,dx1,dy1,dz1)
    CALL solver%initialize_workpool()

    ! initialize tracers
	IF(ndims .EQ. 1) THEN
		CALL initialize_tracer_distribution_1d(tr1d, cellRoot1d, N_cell1d)
        CALL output1d(N_cell1d)
	ELSE IF(ndims .EQ. 2) THEN
		!CALL initialize_tracer_distribution_2d(tr2d, cellHead2d, N_cell2d)
    ELSE IF(ndims .EQ. 3) THEN
		!CALL initialize_tracer_distribution_3d(tr3d, cellHead3d, N_cell3d)
	END IF

    ! Simulation loop  
    t1 = MPI_Wtime()  
    DO i= 1, nt
        PRINT*,' '
        PRINT*,'TIME STEP, % complete = ',i, (i*100./(nt*1.))
        PRINT*,' '

        CALL CPU_TIME(time_init)
        CALL solver%solve(dt1)
        CALL CPU_TIME(time_end)
        time_sim(i) = time_end-time_init
        !CALL output1d(N_cell1d)
        !IF(ndims .EQ. 1)THEN
        !    CALL fileOut1d()
        !ELSE IF(ndims .EQ. 2)THEN
        !    CALL fileOut2d(i)	
        !ELSE IF(ndims .EQ. 3)THEN
        !    CALL fileOut3d_yz(i)	
        !END IF
 

    END DO
    t2 = MPI_Wtime()

    t_tot = t2-t1
    PRINT*,'Total simulation time (sec) =',t_tot
    PRINT*,'Rand num generation Time (sec) = ',randtime
    PRINT*,'Rand time fraction =', randtime/t_tot

    CALL solver%destroy_workpool()



	IF(ndims .EQ. 1) THEN
		DEALLOCATE(tr1d,cellRoot1d,rho,flux,N_cell1d)
    ELSE IF(ndims .EQ. 2) THEN
		!DEALLOCATE(tr2d,cellHead2d,rho2d,flux2d,N_cell2d)
    ELSE IF(ndims .EQ. 3) THEN
		!DEALLOCATE(tr3d,cellHead3d,rho3d,flux3d,N_cell3d)
    END IF

    CLOSE(UNIT=1)

    PRINT*,'DONE!'

	CONTAINS

		! This subroutine assigns a cell to each tracer 
        ! (according to a given distribution). The number of tracers
        ! inside each tracer needs to be specified. 
        ! Each cell contains a tree data structure that has two nodes. Each
        ! node is a linked list. Tracers are divided equaly among these two nodes.
        
		SUBROUTINE initialize_tracer_distribution_1d(tr, cellRoot, N_cell)

			TYPE(tracer_ptr_1d), INTENT(INOUT) :: tr(:)
            TYPE(tree_root_1d), INTENT(INOUT) :: cellRoot(:)
            INTEGER, INTENT(INOUT) :: N_cell(:)

			TYPE(tree_node_1d), POINTER :: current_node
            TYPE(tracer_1d), POINTER :: current_tracer
     		INTEGER :: i, ilow, l, cellnum, num_node1, num_node2
            INTEGER :: tr_counter
            REAL :: p 
        
            ilow=1

            ! Assign a sample tracer distribution for testing
            N_cell(:)=0
            N_cell(ilow+nx/2)=N/2
            N_cell(ilow+1+nx/2)=N-N/2


            ! initialize counter
            tr_counter=1

            ! loop over cells
            DO i = ilow+0, ilow+nx+1



                ! Allocate memory for two nodes
                ALLOCATE(cellRoot(i)%node_L)
                ALLOCATE(cellRoot(i)%node_R)
                ALLOCATE(cellRoot(i)%node_L%next)
                cellRoot(i)%node_R => cellRoot(i)%node_L%next ! for two nodes only, next node = rightmost node
                ALLOCATE(cellRoot(i)%node_occupancy(2))


                PRINT*,'CELL#,  N_cell ',i-ilow, N_cell(i)
  
                IF (N_cell(i) .GT. 0) THEN

                ! Compute number of tracers to be placed inside each node
                num_node1 = N_cell(i)/2
                num_node2 = N_cell(i) - num_node1

                cellRoot(i)%node_occupancy(1) = num_node1
                cellRoot(i)%node_occupancy(2) = num_node2


                !PRINT*,'# of tracers in node 1=',cellRoot(i)%node_occupancy(1)
                !PRINT*,'# of tracers in node 2=',cellRoot(i)%node_occupancy(2)

                !**********************
                ! link tracers to node1
                !**********************

                ! set pointer to first node in the cell
                current_node =>  cellRoot(i)%node_L

                ! allocate memory for first tracer in this node
                ALLOCATE(tr(tr_counter)%p) 
                current_tracer => tr(tr_counter)%p
            
				! assign id and cell number to this tracer 
                tr(tr_counter)%p%id = tr_counter
                tr(tr_counter)%p%x = i-ilow
            
                !PRINT*,'Placing Tracer#',tr(tr_counter)%p%id,' in node1.'

                ! link this tracer to the node
                current_node%leaf_L => tr(tr_counter)%p ! the first tracer in the node is the lefttmost leaf
                tr_counter = tr_counter+1

                ! link with remaining tracers in the node
                DO l = 1, num_node1-1
               
                    ! allocate memory for next tracer in the first node
                    ALLOCATE(tr(tr_counter)%p) 
                    
                    ! assign id and cell number to thid tracer 
                    tr(tr_counter)%p%id = tr_counter
                    tr(tr_counter)%p%x = i-ilow

                    !PRINT*,'Placing Tracer#',tr(tr_counter)%p%id,' in node1.'

                    ! link this tracer to the node
                    current_tracer%next => tr(tr_counter)%p
                    current_tracer =>current_tracer%next 
                    tr_counter = tr_counter+1
                END DO    
                current_node%leaf_R => current_tracer ! the last tracer in the node is the rightmost leaf
                  

                !**********************
                ! link tracers to node2
                !**********************

                ! set pointer to next node in the cell
                current_node =>  cellRoot(i)%node_L%next  ! OR current_node =>  cellRoot(i)%node_R

                ! allocate memory for first tracer pointer in this node
                ALLOCATE(tr(tr_counter)%p) 
                current_tracer => tr(tr_counter)%p
            
				! assign id and cell number to this tracer 
                tr(tr_counter)%p%id = tr_counter
                tr(tr_counter)%p%x = i-ilow

                !PRINT*,'Placing Tracer#',tr(tr_counter)%p%id,' in node2.'

                ! link this tracer to the node
                current_node%leaf_L => tr(tr_counter)%p ! the first tracer in the node is the lefttmost leaf
                tr_counter = tr_counter+1

                ! link with remaining tracers in the node
                DO l = 1, num_node2-1
                
                    ! allocate memory for next tracer pointer in the first node
                    ALLOCATE(tr(tr_counter)%p) 

                    ! assign id and cell number to thid tracer 
                    tr(tr_counter)%p%id = tr_counter
                    tr(tr_counter)%p%x = i-ilow

                    !PRINT*,'Placing Tracer#',tr(tr_counter)%p%id,' in node2.'


                    ! link this tracer to the node
                    current_tracer%next => tr(tr_counter)%p
                    current_tracer =>current_tracer%next 
                    tr_counter = tr_counter+1    
                END DO
                current_node%leaf_R => current_tracer ! the last tracer in the node is the rightmost leaf
               
                END IF
 
			END DO


        	PRINT*, 'Tracer initialization completed.'

		END SUBROUTINE initialize_tracer_distribution_1d


SUBROUTINE output1d(N_cell)

    INTEGER, INTENT(INOUT) :: N_cell(:)
    INTEGER :: i
    REAL :: x

    DO i=1,nx
        x=i*dx1
        WRITE(1,*) x,N_Cell(i)
    END DO

END SUBROUTINE output1d




		
END PROGRAM tracer_test_run

