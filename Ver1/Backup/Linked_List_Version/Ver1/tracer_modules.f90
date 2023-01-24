MODULE constants_mod
	IMPLICIT NONE	
    
	INTEGER, PARAMETER :: N = 50
    INTEGER, PARAMETER :: nx = 5
    INTEGER, PARAMETER :: ny = 1
    INTEGER, PARAMETER :: nz = 1
    INTEGER, PARAMETER :: nb = 3
    INTEGER, PARAMETER :: ndims = 1

END MODULE constants_mod

MODULE data_mod
	IMPLICIT NONE 
    
	INTEGER, ALLOCATABLE :: N_cell(:)
	REAL, ALLOCATABLE :: flux(:) , rho(:)
	
END MODULE data_mod


MODULE tracertype_mod
	IMPLICIT NONE

	! tracer (derived)  types for 1,2 and 3D
	TYPE :: tracer_1d
        ! type data
        INTEGER :: id = 0							 ! label for distinguishing tracers
		INTEGER :: x = 0				 		 	 ! cell location indices
		REAL :: some_fluid_property = 0.0     		 ! fluid property recorded by tracer at it's current location
		TYPE (tracer_1d), POINTER :: next => null()  ! next pointer for linked list capability
	END TYPE tracer_1d

	TYPE, EXTENDS(tracer_1d) :: tracer_2d
		! type data
		INTEGER :: y = 0    						 ! cell location indices
	END TYPE tracer_2d

    TYPE, EXTENDS(tracer_2d) :: tracer_3d
		! type data
		INTEGER :: z = 0    						 ! cell location indices
	END TYPE tracer_3d

	TYPE tracer_ptr_1d
		TYPE (tracer_1d), POINTER :: p => null()
    END TYPE tracer_ptr_1d

END MODULE tracertype_mod


MODULE tracersolver_mod
	USE tracertype_mod
    USE data_mod
	IMPLICIT NONE
	
	! parameters
	INTEGER, PARAMETER, PRIVATE :: debug = 1
	REAL :: dt, dx

    ! work pool variables (1d arrays, for ndims>1 these will hold flattened data)
	REAL, ALLOCATABLE, PRIVATE :: work_flux(:) , work_rho(:)
    INTEGER, ALLOCATABLE, PRIVATE :: work_N_cell(:),work_N_cell_temp(:) 

	! tracer advection solver object
	TYPE MASSFLUXTRACER
        ! type data
        INTEGER, PRIVATE :: ndims, nx, ny, nz, N    

        ! type bound procedures		
		CONTAINS
            PROCEDURE, PASS(SELF), PUBLIC :: initialize_workpool
            PROCEDURE, PASS(SELF), PUBLIC :: destroy_workpool
            PROCEDURE, PASS(SELF), PUBLIC :: solve
			PROCEDURE, PASS(SELF), PUBLIC :: advect_1d
			PROCEDURE, PASS(SELF), PRIVATE :: cell_transfer_1d


    END TYPE MASSFLUXTRACER

	INTERFACE MASSFLUXTRACER
		MODULE PROCEDURE constructor
	END INTERFACE 

    
	CONTAINS

		FUNCTION constructor(ndims,nx,ny,nz,N) RESULT(SELF)

			TYPE(MASSFLUXTRACER) :: SELF
            INTEGER, INTENT(IN) :: ndims, nx, ny, nz, N

			SELF%ndims = ndims
			SELF%nx = nx
			SELF%ny = ny
			SELF%nz = nz
            SELF%N = N 

			RETURN

		END FUNCTION constructor

        SUBROUTINE initialize_workpool(SELF)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF

			INTEGER :: ndims,nx,ny,nz

			ndims = SELF%ndims
            nx = SELF%nx
            ny = SELF%ny
            nz = SELF%nz
			

			ALLOCATE(work_rho(1:MAX(1,nx)*MAX(1,ny)*MAX(1,nz)))
 
			IF(ndims .EQ. 1) THEN
				ALLOCATE(work_flux(1:nx+1))
				ALLOCATE(work_N_cell(1:nx+2),work_N_cell_temp(1:nx+2))
			END IF


            IF(ndims .EQ. 2) THEN
				ALLOCATE(work_flux(1:(nx+1)*(ny+1)))
				ALLOCATE(work_N_cell(1:(nx+2)*(ny+2)),work_N_cell_temp(1:(nx+2)*(ny+2)) )
			END IF

            IF(ndims .EQ. 3) THEN
				ALLOCATE(work_flux(1:(nx+1)*(ny+1)*(nz+1)))
				ALLOCATE(work_N_cell(1:(nx+2)*(ny+2)*(nz+2)),work_N_cell_temp(1:(nx+2)*(ny+2)*(nz+2)))
			END IF

			work_rho = 0.0
			work_flux = 0.0
			work_N_cell = 0
			work_N_cell_temp = 0	  

		END SUBROUTINE initialize_workpool


		SUBROUTINE destroy_workpool(SELF)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF

			DEALLOCATE(work_rho, work_flux, &
                       work_N_cell, work_N_cell_temp)

		END SUBROUTINE destroy_workpool


		! Top-level interface for mass-flux tracer advection
		SUBROUTINE solve(SELF, tr, cellHead1d, dt_in, dx_in)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
			TYPE(tracer_ptr_1d), INTENT(INOUT) :: tr(:),cellHead1d(:)
            REAL, INTENT(IN) :: dt_in, dx_in
			
            dt = dt_in
            dx = dx_in


            IF(SELF%ndims .EQ. 1) Then
				PRINT*,'CHECKPOINT1'
				CALL SELF%advect_1d(tr,cellHead1d,N_cell,flux,rho)
			END IF

        END SUBROUTINE solve

        ! Possible Improvements: Instead of sending density and flux arrays,
        ! have the mhd solver compute the fluid outgoing mass fraction through 
        ! cell faces and send that.
		SUBROUTINE advect_1d(SELF, tr, cellHead1d, N_cell, flux, rho)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF 
			TYPE(tracer_ptr_1d), INTENT(INOUT) :: tr(:),cellHead1d(:)
            INTEGER, INTENT(INOUT) :: N_cell(:)
            REAL, INTENT(IN) :: flux(:), rho(:)
			
			INTEGER :: nx           
			TYPE(tracer_1d), POINTER :: current
     		INTEGER :: i, j, ntr, move_count
            REAL :: p, m_frac, m0

			nx = SELF%nx

			PRINT*,'CHECKPOINT2'


            ! copy values into work pool arrays
            work_N_cell(:) = N_cell(:)
            work_N_cell_temp(:) = N_cell(:)
			work_flux(:) = flux(:)
            work_rho(:) = rho(:) 
            

			! Loop over cells, compute outgoing fluid mass fractions (m_out) through 
            ! each face and transfer outgoing tracers into neighboring cells.
            !
            ! Use Monte Carlo sampling to determine whether or not to move a tracer:
			! For each cell face: for each tracer (originally in the cell) draw a 
			! random number p from (0,1). If p < m_out, then move 
			! that tracer across the cell face into neighbor cell, otherwise don't.
	
            DO i = 1, nx

				current => cellHead1d(i+1)%p
               
                ntr= 0
                
                IF(debug .EQ. 1) THEN
                    PRINT*, ' '
					PRINT*,'Cell #, N_cell:',i, work_N_cell(i+1)
					PRINT*,'F_L, F_R = ', work_flux(i), work_flux(i+1)
                    PRINT*, ' '
                END IF
				   
				! non-empty cell 
                IF(ASSOCIATED(current)) THEN 

			       	! calculate fluid mass in cell
                	m0=work_rho(i)*dx

                   IF(debug .EQ. 1) THEN
						PRINT*,'Tracers initially occupying this cell:'
						PRINT*,'Tracer # ',current%id
                        DO WHILE(ASSOCIATED(current%next))
							current => current%next
                            PRINT*,'Tracer # ',current%id
						END DO
                        current => cellHead1d(i+1)%p
	        	    END IF
 
                    move_count = 0

                	! Outgoing through right face
                	IF(SIGN(1.0,work_flux(i+1)) .GT. 0.0) THEN
                		
						m_frac = work_flux(i+1)*dt/m0

                        IF(debug .EQ. 1) PRINT*,'OUTGOING RIGHT: m0, m_frac=',m0,m_frac 

                        m0=m0*(1.0-m_frac)

						! Loop over each tracer originally in the cell
						DO j = 1, work_N_cell(i+1)

               				CALL RANDOM_NUMBER(p)
                             
                            IF(debug .EQ. 1) THEN
							    PRINT*,' Tracer #, p, m_frac = ',current%id,p,m_frac
	        	            END IF

               				IF(p .LT. m_frac)THEN

								IF(debug .EQ. 1) THEN
								    PRINT*,' Moving Tracer # ',current%id,' into cell ',i+1
                           			IF(debug .EQ. 1) current => current%next 
		        	            END IF

								CALL SELF%cell_transfer_1d(tr,cellHead1d,work_N_cell_temp(i+1),&
 												   	 work_N_cell_temp(i+2),j-move_count,i,i+1)

                                work_N_cell_temp(i+1) = work_N_cell_temp(i+1)-1 
                                work_N_cell_temp(i+2) = work_N_cell_temp(i+2)+1
                                move_count = move_count+1
                                ntr = ntr+1
                                IF(debug .EQ. 1) PRINT*,'ntr=',ntr

 							END IF

                            IF((p .GT. m_frac) .AND. (debug .EQ. 1)) current => current%next
                         END DO
                	END IF			

                    ! reset counter
                	move_count = 0
	
					! Outgoing through left face
                	IF(SIGN(1.0,work_flux(i)) .LT. 0.0) THEN

                		m_frac = abs(work_flux(i))*dt/m0

                        IF(debug .EQ. 1) THEN
							PRINT*,'OUTGOING LEFT: m0, m_frac=',m0,m_frac
                        END IF

                        IF(debug .EQ. 1) current => cellHead1d(i+1)%p

                     	! Loop over each tracer in the cell  
						DO j = 1, work_N_cell(i+1)-ntr

							CALL RANDOM_NUMBER(p)

							IF(debug .EQ. 1) THEN
							    PRINT*,' Tracer #, p, m_frac = ',current%id,p,m_frac
	        	            END IF

               				IF(p .LT. m_frac)THEN

								IF(debug .EQ. 1) THEN
								    PRINT*,' Moving Tracer # ',current%id,' into cell ',i-1 
                                    current => current%next
		        	            END IF

								CALL SELF%cell_transfer_1d(tr,cellHead1d, work_N_cell_temp(i+1), &
                                                     work_N_cell_temp(i),j-move_count,i,i-1)

								work_N_cell_temp(i+1) = work_N_cell_temp(i+1)-1 
                                work_N_cell_temp(i) = work_N_cell_temp(i)+1
                                move_count = move_count+1

                            END IF
                            IF((p .GT. m_frac) .AND. (debug .EQ. 1)) current => current%next
						END DO

               	    END IF
				END IF

			END DO		
			
			! update N_cell array
            N_cell(:) = work_N_cell_temp(:)


		END SUBROUTINE advect_1d


		SUBROUTINE cell_transfer_1d(SELF, tr, cellHead1d, N_src, N_des, tr_num ,x_src, x_des)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
			TYPE(tracer_ptr_1d), INTENT(INOUT) :: tr(:),cellHead1d(:)
			INTEGER, INTENT(IN) :: tr_num, N_src, N_des, x_src, x_des ! tr_num: tracer number in linked list 
														! x_src: source cell
														! x_des: destination cell								
            INTEGER :: i
            TYPE(tracer_1d), POINTER :: current_src, current_des, temp, &
										current_src_above => null(), current_src_below => null() 


            IF(debug .EQ. 1) THEN
                temp => cellHead1d(x_src+1)%p 
				PRINT*,'Tracers currently in source cell:'
				IF(ASSOCIATED(temp)) THEN
					 PRINT*,'Tracer # ',temp%id
    	            DO WHILE(ASSOCIATED(temp%next))
						temp => temp%next
        	      		PRINT*,'Tracer # ',temp%id
					END DO
				END IF

                temp => cellHead1d(x_des+1)%p 
				PRINT*,'Tracers currently in destination cell:'
				IF(ASSOCIATED(temp)) THEN
                	PRINT*,'Tracer # ',temp%id
               	 	DO WHILE(ASSOCIATED(temp%next))
						temp => temp%next
                    	PRINT*,'Tracer # ',temp%id
					END DO
				END IF
	        END IF


            
            !#####################################################################
            ! CORE STEP: Link the last tracer in the destination cell to the
			! outgoing tracer from the source cell. Then nullify the next pointer 
			! of this outgoing tracer.
            !#####################################################################
             
			! set a temp pointer to head tracer in source cell
			current_src => cellHead1d(x_src+1)%p

            ! traverse to the outgoing tracer location within link
			DO i = 1, tr_num-1
                current_src_above => current_src
				current_src => current_src%next 
				IF(debug .EQ. 1) PRINT*,'Traversing list. Currently &
                                         pointing at Tracer #',current_src%id                    

       	    END DO
               
            ! set a temp pointer to tracer next in list after the outgoing tracer
            IF(ASSOCIATED(current_src%next)) current_src_below => current_src%next

			! set a temp pointer to head tracer in destination cell
			current_des => cellHead1d(x_des+1)%p

			! if destination cell empty, then incoming tracer is assigned as head	
            IF(.NOT. ASSOCIATED(current_des))THEN
				cellHead1d(x_des+1)%p => current_src
				NULLIFY(current_src%next)

				IF(debug .EQ. 1) PRINT*,'Destination cell empty.'

            ! otherwise, traverse to last tracer of destination cell
            ELSE
				
            	DO WHILE(ASSOCIATED(current_des%next))
            	  	current_des => current_des%next
                END DO
                
               	! link last tracer of destination cell to incoming
               	! tracer from source cell
               	ALLOCATE(current_des%next)
               	current_des%next => current_src
	            NULLIFY(current_src%next)

			END IF


            !#####################################################################
            ! AUXILIARY STEPS: Reattach broken link left behind in the source cell
			! by the outgoing tracer.
			!
			! 	Caveats:  1) If outoing tracer is the head of source cell, then 
			!				 need to promote the second tracer in the source 
            !				 cell to new head, assuming there is a second tracer.
			!				 (one-sided link re-attach)
			!
			!			  2) If outgoing tracer is tail of source cell, then 
			! 				 nulify next pointer of second to last tracer in source
			!                cell.
			! 				(once-sided link re-attach)
			!
			!		      3) If outgoing tracer is neither head nor tail, then
			!				 link the tracer that is one above it in the source cell
			!				 to the tracer that is one below it.
			!				(two-sided link re-attach)
            !#####################################################################

			IF(tr_num .EQ. 1) THEN	

				IF(debug .EQ. 1) PRINT*,'Case 1. Moving Head tracer.'

				! if more than 1 tracer in cell, promote next tracer to head
                IF(ASSOCIATED(current_src_below))THEN
                	cellHead1d(x_src+1)%p => current_src_below
                    IF(debug .EQ. 1) PRINT*,'Tracer # ', current_src_below%id ,' promoted to head.'
				! otherwise destroy head tracer pointer
                ELSE
                	NULLIFY(cellHead1d(x_src+1)%p)
				END IF

			END IF

			IF(tr_num .GT. 1 .AND. tr_num .EQ. N_src) THEN

				IF(debug .EQ. 1) PRINT*,'Case 2. Moving Tail tracer.'

				NULLIFY(current_src_above%next)
			END IF
			
			IF(tr_num .GT. 1 .AND. tr_num .LT. N_src) THEN

				IF(debug .EQ. 1) PRINT*,'Case 3. Neither Head nor Tail.'
	
				current_src_above%next => current_src_below
			END IF

			!#####################################################################

            IF(debug .EQ. 1) THEN
                temp => cellHead1d(x_src+1)%p 
				PRINT*,'Tracers currently in source cell:'
				IF(ASSOCIATED(temp)) THEN 
					PRINT*,'Tracer # ',temp%id
                	DO WHILE(ASSOCIATED(temp%next))
						temp => temp%next
              			PRINT*,'Tracer # ',temp%id
					END DO
				END IF
                temp => cellHead1d(x_des+1)%p 
				PRINT*,'Tracers currently in destination cell:'
				IF(ASSOCIATED(temp)) THEN
				PRINT*,'Tracer # ',temp%id
	                DO WHILE(ASSOCIATED(temp%next))
						temp => temp%next
	              		PRINT*,'Tracer # ',temp%id
					END DO
				END IF
	        END IF
  

		END SUBROUTINE cell_transfer_1d

END MODULE tracersolver_mod

