MODULE tracersolver_mod
	USE tracertype_mod
    USE data_mod
	IMPLICIT NONE
	
	! parameters
	INTEGER, PARAMETER, PRIVATE :: debug = 1
	REAL :: dt, dx, dy, dz

    ! work pool variables 
	REAL, ALLOCATABLE, PRIVATE :: work_flux(:) , work_rho(:)
    REAL, ALLOCATABLE, PRIVATE :: work_flux2d(:,:,:) , work_rho2d(:,:)
    REAL, ALLOCATABLE, PRIVATE :: work_flux3d(:,:,:,:) , work_rho3d(:,:,:) 
    INTEGER, ALLOCATABLE, PRIVATE :: work_N_cell(:),work_N_cell_temp(:)
    INTEGER, ALLOCATABLE, PRIVATE :: work_N_cell2d(:,:),work_N_cell_temp2d(:,:) 
    INTEGER, ALLOCATABLE, PRIVATE :: work_N_cell3d(:,:,:),work_N_cell_temp3d(:,:,:)

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

			PROCEDURE, PASS(SELF), PUBLIC :: advect_2d
			PROCEDURE, PASS(SELF), PRIVATE :: cell_transfer_2d

			PROCEDURE, PASS(SELF), PUBLIC :: advect_3d
			PROCEDURE, PASS(SELF), PRIVATE :: cell_transfer_3d


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
			
			IF(ndims .EQ. 1) THEN
                ALLOCATE(work_rho(1:nx))
				ALLOCATE(work_flux(1:nx+1))
				ALLOCATE(work_N_cell(1:nx+2),work_N_cell_temp(1:nx+2))

				work_rho = 0.0
				work_flux = 0.0
				work_N_cell = 0
				work_N_cell_temp = 0	  

			END IF

            IF(ndims .EQ. 2) THEN
                ALLOCATE(work_rho2d(1:nx,1:ny))
				ALLOCATE(work_flux2d(1:nx+1,1:ny+1,2))
				ALLOCATE(work_N_cell2d(1:nx+2,1:ny+2),work_N_cell_temp2d(1:nx+2,1:ny+2))

				work_rho2d = 0.0
				work_flux2d = 0.0
				work_N_cell2d = 0
				work_N_cell_temp2d = 0	  

			END IF

            IF(ndims .EQ. 3) THEN
				ALLOCATE(work_rho3d(1:nx,1:ny,1:nz))
				ALLOCATE(work_flux3d(1:nx+1,1:ny+1,1:nz+1,3))
				ALLOCATE(work_N_cell3d(1:nx+2,1:ny+2,1:nz+2),work_N_cell_temp3d(1:nx+2,1:ny+2,1:nz+2))

				work_rho3d = 0.0
				work_flux3d = 0.0
				work_N_cell3d = 0
				work_N_cell_temp3d = 0	  

			END IF

			

		END SUBROUTINE initialize_workpool


		SUBROUTINE destroy_workpool(SELF)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF

            IF(SELF%ndims .EQ. 1) THEN
				DEALLOCATE(work_rho, work_flux, &
						   work_N_cell, work_N_cell_temp)
			ELSE IF(SELF%ndims .EQ. 2) THEN
				DEALLOCATE(work_rho2d, work_flux2d, &
						   work_N_cell2d, work_N_cell_temp2d)
            ELSE IF(SELF%ndims .EQ. 3) THEN
				DEALLOCATE(work_rho3d, work_flux3d, &
						   work_N_cell3d, work_N_cell_temp3d)
			END IF


		END SUBROUTINE destroy_workpool


		! Top-level routine for mass-flux tracer advection

		SUBROUTINE solve(SELF, dt_in, dx_in, dy_in, dz_in)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
            REAL, INTENT(IN) :: dt_in, dx_in, dy_in, dz_in
			
            dt = dt_in
            dx = dx_in
            dy = dy_in
            dz = dz_in

            SELECT CASE(SELF%ndims) 
            
            CASE(1) 
				CALL SELF%advect_1d(tr1d,cellHead1d,N_cell,flux,rho)
			CASE(2)
                CALL SELF%advect_2d(tr2d,cellHead2d,N_cell2d,flux2d,rho2d)
            CASE(3)
                CALL SELF%advect_3d(tr3d,cellHead3d,N_cell3d,flux3d,rho3d)

			END SELECT


        END SUBROUTINE solve



        ! Routine for 3d tracer advection
		!
		! Loop over cells. Cycle through the cell faces (in the order x+, 
        ! x-, y+, y-, z+, z-). For each face: 
		!
		! 1) Compute outgoing fluid mass fractions (m_frac).
		!
		! e.g. through right face, m_frac = mass_flux*dy*dz*dt/m0
        !      where m0 is the mass of fluid currently in the cell (which is
		!	   updated after each cell face cycle: m0_new = m0_old*(1-m_frac) ).
		!
        ! 2) Use Monte Carlo sampling to determine whether or not to move a tracer
		! 	 through each cell face. Cycle through every tracer currently occupying  
		!    that cell. For each tracer, draw a random number p from [0,1).  
		!    If p < m_frac, then move that tracer, otherwise don't.	

		SUBROUTINE advect_3d(SELF, tr, cellHead, N_cell, flux, rho)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF 
			TYPE(tracer_ptr_3d), INTENT(INOUT) :: tr(:),cellHead(:,:,:)
            INTEGER, INTENT(INOUT) :: N_cell(:,:,:)
            REAL, INTENT(IN) :: flux(:,:,:,:), rho(:,:,:)
			           
			TYPE(tracer_3d), POINTER :: current
            INTEGER :: nx, ny, nz
     		INTEGER :: i, j, k, l, ntr, move_count, face_count
            REAL :: m0, m_frac, p

            ! clear temp pointer
            current => null()

			nx = SELF%nx
            ny = SELF%ny
            nz = SELF%nz     

            ! copy values into work pool arrays
            work_N_cell3d(:,:,:) = N_cell(:,:,:)
            work_N_cell_temp3d(:,:,:) = N_cell(:,:,:)
			work_flux3d(:,:,:,:) = flux(:,:,:,:)
            work_rho3d(:,:,:) = rho(:,:,:) 
            
			! loop across cells in our 2d grid
            DO i = 1, nx
				DO j = 1, ny
                    DO k = 1, nz
                    
                    ! reset counters
					ntr= 0
					move_count = 0
					face_count = 0
               
					! calculate initial fluid mass in cell
					m0 = work_rho3d(i,j,k)*dx*dy*dz

					IF(debug .EQ. 1) THEN
                        PRINT*, ' '
                        PRINT*,'Cell #:',i,j,k,' , N_cell: ',work_N_cell3d(i+1,j+1,k+1)
                        PRINT*, ' '
                    END IF
				   
				    ! cycle through all six cell faces
					DO WHILE(face_count .LT. 6 .AND. work_N_cell_temp3d(i+1,j+1,k+1) .GT. 0 &
							 .AND. work_N_cell3d(i+1,j+1,k+1) .GT. 0)

						current => cellHead(i+1,j+1,k+1)%p

						IF(debug .EQ. 1) THEN
							PRINT*,'Tracers initially occupying this cell:'
							PRINT*,'Tracer # ',current%id
							DO WHILE(ASSOCIATED(current%next))
								current => current%next
								PRINT*,'Tracer # ',current%id
							END DO
							current => cellHead(i+1,j+1,k+1)%p
						END IF
 

						SELECT CASE(face_count)

						CASE(0)	
							! Outgoing through x+ face
							IF(SIGN(1.0,work_flux3d(i+1,j+1,k+1,1)) .GT. 0.0) THEN
                		
								m_frac = work_flux3d(i+1,j+1,k+1,1)*dt*dy*dz/m0

								IF(debug .EQ. 1) PRINT*,'OUTGOING x+: m0, m_frac=',m0,m_frac 

								m0=m0*(1.0-m_frac)

								! Loop over each tracer originally in the cell
								DO l = 1, work_N_cell3d(i+1,j+1,k+1)

									CALL RANDOM_NUMBER(p)
                  
									IF(debug .EQ. 1) THEN
										PRINT*,' Tracer #, p, m_frac = ',current%id,p,m_frac
									END IF

									IF(p .LT. m_frac)THEN

										IF(debug .EQ. 1) THEN
											PRINT*,' Moving Tracer # ',current%id,' into cell ',i+1,j,k
										END IF

										current%x = i+1
										current => current%next
										CALL SELF%cell_transfer_3d(tr,cellHead,work_N_cell_temp3d(i+1,j+1,k+1),&
 												   	 work_N_cell_temp3d(i+2,j+1,k+1),l-move_count,i,i+1,j,j,k,k)

										work_N_cell_temp3d(i+1,j+1,k+1) = work_N_cell_temp3d(i+1,j+1,k+1)-1 
										work_N_cell_temp3d(i+2,j+1,k+1) = work_N_cell_temp3d(i+2,j+1,k+1)+1
										move_count = move_count+1
									ELSE
										current => current%next
									END IF

								END DO

                                ntr= ntr + move_count
								! reset counter
								move_count = 0
							END IF			
	
						CASE(1)
							! Outgoing through x- face
							IF(SIGN(1.0,work_flux3d(i,j+1,k+1,1)) .LT. 0.0) THEN

								m_frac = abs(work_flux3d(i,j+1,k+1,1))*dt*dy*dz/m0

								IF(debug .EQ. 1) THEN
									PRINT*,'OUTGOING x-: m0, m_frac=',m0,m_frac
								END IF

								current => cellHead(i+1,j+1,k+1)%p

								! Loop over each tracer in the cell  
								DO l = 1, work_N_cell3d(i+1,j+1,k+1)-ntr

									CALL RANDOM_NUMBER(p)

									IF(debug .EQ. 1) THEN
										PRINT*,' Tracer #, p, m_frac = ',current%id,p,m_frac
									END IF

									IF(p .LT. m_frac)THEN

										IF(debug .EQ. 1) THEN
											PRINT*,' Moving Tracer # ',current%id,' into cell ',i-1,j,k 
										END IF

										current%x = i-1
										current => current%next

										CALL SELF%cell_transfer_3d(tr,cellHead, work_N_cell_temp3d(i+1,j+1,k+1), &
                                                     work_N_cell_temp3d(i,j+1,k+1),l-move_count,i,i-1,j,j,k,k)

										work_N_cell_temp3d(i+1,j+1,k+1) = work_N_cell_temp3d(i+1,j+1,k+1)-1 
										work_N_cell_temp3d(i,j+1,k+1) = work_N_cell_temp3d(i,j+1,k+1)+1
										move_count = move_count+1                                        
									ELSE
										current => current%next
									END IF

								END DO
  
                                ntr = ntr+move_count
								! reset counter
								move_count = 0
							END IF

						CASE(2)
							! Outgoing through y+ face
							IF(SIGN(1.0,work_flux3d(i+1,j+1,k+1,2)) .GT. 0.0) THEN
                		
								m_frac = work_flux3d(i+1,j+1,k+1,2)*dt*dx*dz/m0

								IF(debug .EQ. 1) PRINT*,'OUTGOING y+: m0, m_frac=',m0,m_frac 

								current => cellHead(i+1,j+1,k+1)%p
								m0=m0*(1.0-m_frac)

								! Loop over each tracer originally in the cell
								DO l = 1, work_N_cell3d(i+1,j+1,k+1)-ntr

									CALL RANDOM_NUMBER(p)
                             
									IF(debug .EQ. 1) THEN
										PRINT*,' Tracer #, p, m_frac = ',current%id,p,m_frac
									END IF

									IF(p .LT. m_frac)THEN

										IF(debug .EQ. 1) THEN
											PRINT*,' Moving Tracer # ',current%id,' into cell ',i,j+1,k
										END IF

										current%y = j+1
										current => current%next

										CALL SELF%cell_transfer_3d(tr,cellHead,work_N_cell_temp3d(i+1,j+1,k+1),&
 												   	 work_N_cell_temp3d(i+1,j+2,k+1),l-move_count,i,i,j,j+1,k,k)

										work_N_cell_temp3d(i+1,j+1,k+1) = work_N_cell_temp3d(i+1,j+1,k+1)-1 
										work_N_cell_temp3d(i+1,j+2,k+1) = work_N_cell_temp3d(i+1,j+2,k+1)+1
										move_count = move_count+1
									ELSE
										current => current%next
									END IF

								END DO

                                ntr = ntr+move_count 
								! reset counter
								move_count = 0
							END IF		

						CASE(3)
							! Outgoing through y- face
							IF(SIGN(1.0,work_flux3d(i+1,j,k+1,2)) .LT. 0.0) THEN
                		
								m_frac = abs(work_flux3d(i+1,j,k+1,2))*dt*dx*dz/m0

								IF(debug .EQ. 1) PRINT*,'OUTGOING y-: m0, m_frac=',m0,m_frac 

								current => cellHead(i+1,j+1,k+1)%p

								m0=m0*(1.0-m_frac)

								! Loop over each tracer originally in the cell
								DO l = 1, work_N_cell3d(i+1,j+1,k+1)-ntr

									CALL RANDOM_NUMBER(p)
                             
									IF(debug .EQ. 1) THEN
										PRINT*,' Tracer #, p, m_frac = ',current%id,p,m_frac
									END IF

									IF(p .LT. m_frac)THEN

										IF(debug .EQ. 1) THEN
											PRINT*,' Moving Tracer # ',current%id,' into cell ',i,j-1,k
										END IF

										current%y = j-1
										current => current%next

										CALL SELF%cell_transfer_3d(tr,cellHead,work_N_cell_temp3d(i+1,j+1,k+1),&
 												   	 work_N_cell_temp3d(i+1,j,k+1),l-move_count,i,i,j,j-1,k,k)

										work_N_cell_temp3d(i+1,j+1,k+1) = work_N_cell_temp3d(i+1,j+1,k+1)-1 
										work_N_cell_temp3d(i+1,j,k+1) = work_N_cell_temp3d(i+1,j,k+1)+1
										move_count = move_count+1
									ELSE
										current => current%next
									END IF

								END DO

                                ntr = ntr+move_count 
								! reset counter
								move_count = 0
							END IF			

                        CASE(4)
							! Outgoing through z+ face
							IF(SIGN(1.0,work_flux3d(i+1,j+1,k+1,3)) .GT. 0.0) THEN
                		
								m_frac = work_flux3d(i+1,j+1,k+1,3)*dt*dx*dy/m0

								IF(debug .EQ. 1) PRINT*,'OUTGOING z+: m0, m_frac=',m0,m_frac 

								current => cellHead(i+1,j+1,k+1)%p
								m0=m0*(1.0-m_frac)

								! Loop over each tracer originally in the cell
								DO l = 1, work_N_cell3d(i+1,j+1,k+1)-ntr

									CALL RANDOM_NUMBER(p)
                             
									IF(debug .EQ. 1) THEN
										PRINT*,' Tracer #, p, m_frac = ',current%id,p,m_frac
									END IF

									IF(p .LT. m_frac)THEN

										IF(debug .EQ. 1) THEN
											PRINT*,' Moving Tracer # ',current%id,' into cell ',i,j,k+1
										END IF

										current%z = k+1
										current => current%next

										CALL SELF%cell_transfer_3d(tr,cellHead,work_N_cell_temp3d(i+1,j+1,k+1),&
 												   	 work_N_cell_temp3d(i+1,j+1,k+2),l-move_count,i,i,j,j,k,k+1)

										work_N_cell_temp3d(i+1,j+1,k+1) = work_N_cell_temp3d(i+1,j+1,k+1)-1 
										work_N_cell_temp3d(i+1,j+1,k+2) = work_N_cell_temp3d(i+1,j+1,k+2)+1
										move_count = move_count+1
									ELSE
										current => current%next
									END IF

								END DO

                                ntr = ntr+move_count 
								! reset counter
								move_count = 0
							END IF		

						CASE(5)
							! Outgoing through z- face
							IF(SIGN(1.0,work_flux3d(i+1,j+1,k,2)) .LT. 0.0) THEN
                		
								m_frac = abs(work_flux3d(i+1,j+1,k,2))*dt*dx*dy/m0

								IF(debug .EQ. 1) PRINT*,'OUTGOING z-: m0, m_frac=',m0,m_frac 

								current => cellHead(i+1,j+1,k+1)%p

								m0=m0*(1.0-m_frac)

								! Loop over each tracer originally in the cell
								DO l = 1, work_N_cell3d(i+1,j+1,k+1)-ntr

									CALL RANDOM_NUMBER(p)
                             
									IF(debug .EQ. 1) THEN
										PRINT*,' Tracer #, p, m_frac = ',current%id,p,m_frac
									END IF

									IF(p .LT. m_frac)THEN

										IF(debug .EQ. 1) THEN
											PRINT*,' Moving Tracer # ',current%id,' into cell ',i,j,k-1
										END IF

										current%z = k-1
										current => current%next

										CALL SELF%cell_transfer_3d(tr,cellHead,work_N_cell_temp3d(i+1,j+1,k+1),&
 												   	 work_N_cell_temp3d(i+1,j+1,k),l-move_count,i,i,j,j,k,k-1)

										work_N_cell_temp3d(i+1,j+1,k+1) = work_N_cell_temp3d(i+1,j+1,k+1)-1 
										work_N_cell_temp3d(i+1,j+1,k) = work_N_cell_temp3d(i+1,j+1,k)+1
										move_count = move_count+1
									ELSE
										current => current%next
									END IF

								END DO
							END IF	
						END SELECT

                        face_count = face_count +1

                    END DO

                    
                    END DO
				END DO
			END DO		
			
			! update N_cell array
            N_cell(:,:,:) = work_N_cell_temp3d(:,:,:)


		END SUBROUTINE advect_3d




        ! routine for performing a single tracer exchange between two cells on 3d grid

		SUBROUTINE cell_transfer_3d(SELF, tr, cellHead, N_src, N_des, tr_num ,x_src, x_des, &
                                    y_src, y_des, z_src, z_des)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
			TYPE(tracer_ptr_3d), INTENT(INOUT) :: tr(:),cellHead(:,:,:)
			INTEGER, INTENT(IN) :: tr_num, N_src, N_des, x_src, x_des, y_src, y_des, z_src, z_des 
                                                        ! tr_num: tracer position in linked list 
                                                        ! N_src : num of tracers in source cell
                                                        ! N_des : num of tracers in destination cell
														! x_src, y_src, z_src : source cell indices
														! x_des, y_des, z_des : destination cell indices	
							
            INTEGER :: i
            TYPE(tracer_3d), POINTER :: current_src, current_des, temp, &
										current_src_above => null(), current_src_below => null() 


            ! clear temp pointers
            current_src => null()
            current_src_above => null()
            current_src_below => null()
            current_des => null()
            temp => null()



            IF(debug .EQ. 1) THEN
                temp => cellHead(x_src+1, y_src+1, z_src+1)%p 
				PRINT*,'Tracers currently in source cell:'
				IF(ASSOCIATED(temp)) THEN
					 PRINT*,'Tracer # ',temp%id
    	            DO WHILE(ASSOCIATED(temp%next))
						temp => temp%next
        	      		PRINT*,'Tracer # ',temp%id
					END DO
				END IF

                temp => cellHead(x_des+1,y_des+1, z_des+1)%p 
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
			current_src => cellHead(x_src+1,y_src+1,z_src+1)%p

            ! traverse to the outgoing tracer location within link
			DO i = 1, tr_num-1
                current_src_above => current_src
				current_src => current_src%next 
				!IF(debug .EQ. 1) PRINT*,'Traversing list. Currently &
                !                         pointing at Tracer #',current_src%id                    

       	    END DO
               
            ! set a temp pointer to tracer next in list after the outgoing tracer
            IF(ASSOCIATED(current_src%next)) current_src_below => current_src%next

			! set a temp pointer to head tracer in destination cell
			current_des => cellHead(x_des+1,y_des+1,z_des+1)%p

			! if destination cell empty, then incoming tracer is assigned as head	
            IF(.NOT. ASSOCIATED(current_des))THEN
				cellHead(x_des+1,y_des+1,z_des+1)%p => current_src
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
			! 				 (once-sided link re-attach)
			!
			!		      3) If outgoing tracer is neither head nor tail, then
			!				 link the tracer that is one above it in the source cell
			!				 to the tracer that is one below it.
			!				 (two-sided link re-attach)
            !#####################################################################

			IF(tr_num .EQ. 1) THEN	

				IF(debug .EQ. 1) PRINT*,'Case 1. Moving Head tracer.'

				! if more than 1 tracer in cell, promote next tracer to head
                IF(ASSOCIATED(current_src_below))THEN
                	cellHead(x_src+1,y_src+1,z_src+1)%p => current_src_below
                    IF(debug .EQ. 1) PRINT*,'Tracer # ', current_src_below%id ,' promoted to head.'
				! otherwise destroy head tracer pointer
                ELSE
                	NULLIFY(cellHead(x_src+1,y_src+1,z_src+1)%p)
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
                temp => cellHead(x_src+1,y_src+1,z_src+1)%p 
				PRINT*,'Tracers currently in source cell:'
				IF(ASSOCIATED(temp)) THEN 
					PRINT*,'Tracer # ',temp%id
                	DO WHILE(ASSOCIATED(temp%next))
						temp => temp%next
              			PRINT*,'Tracer # ',temp%id
					END DO
				END IF
                temp => cellHead(x_des+1,y_des+1,z_des+1)%p 
				PRINT*,'Tracers currently in destination cell:'
				IF(ASSOCIATED(temp)) THEN
				PRINT*,'Tracer # ',temp%id
	                DO WHILE(ASSOCIATED(temp%next))
						temp => temp%next
	              		PRINT*,'Tracer # ',temp%id
					END DO
				END IF
	        END IF
  

		END SUBROUTINE cell_transfer_3d



        ! Routine for 2d tracer advection
		!
		! Loop over cells. Cycle through the cell faces (in the order x+, 
        ! x-, y+, y-). For each face: 
		!
		! 1) Compute outgoing fluid mass fractions (m_frac).
		!
		! e.g. through right face, m_frac = mass_flux*dy*dt/m0
        !      where m0 is the mass of fluid currently in the cell (which is
		!	   updated after each cell face cycle: m0_new = m0_old*(1-m_frac) ).
		!
        ! 2) Use Monte Carlo sampling to determine whether or not to move a tracer
		! 	 through each cell face. Cycle through every tracer currently occupying  
		!    that cell. For each tracer, draw a random number p from [0,1).  
		!    If p < m_frac, then move that tracer, otherwise don't.	

		SUBROUTINE advect_2d(SELF, tr, cellHead, N_cell, flux, rho)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF 
			TYPE(tracer_ptr_2d), INTENT(INOUT) :: tr(:),cellHead(:,:)
            INTEGER, INTENT(INOUT) :: N_cell(:,:)
            REAL, INTENT(IN) :: flux(:,:,:), rho(:,:)
			           
			TYPE(tracer_2d), POINTER :: current
            INTEGER :: nx, ny
     		INTEGER :: i, j, l, ntr, move_count, face_count
            REAL :: m0, m_frac, p

            ! clear temp pointer
            current => null()

			nx = SELF%nx
            ny = SELF%ny

            ! copy values into work pool arrays
            work_N_cell2d(:,:) = N_cell(:,:)
            work_N_cell_temp2d(:,:) = N_cell(:,:)
			work_flux2d(:,:,:) = flux(:,:,:)
            work_rho2d(:,:) = rho(:,:) 
            
			! loop across cells in our 2d grid
            DO i = 1, nx
				DO j = 1, ny

                    ! reset counters
					ntr= 0
					move_count = 0
					face_count = 0
               
					! calculate initial fluid mass in cell
					m0 = work_rho2d(i,j)*dx*dy

					IF(debug .EQ. 1) THEN
                        PRINT*, ' '
                        PRINT*,'Cell #:',i,j,' , N_cell: ',work_N_cell2d(i+1,j+1)
                        !PRINT*,'Flux(left, right, bottom, top) :',work_flux2d(i,j+1,1), &
                        !		work_flux2d(i+1,j+1,1),work_flux2d(i+1,j,2),work_flux2d(i+1,j+1,2)
                        PRINT*, ' '
                    END IF
				   
				    ! cycle through all four cell faces
					DO WHILE(face_count .LT. 4 .AND. work_N_cell_temp2d(i+1,j+1) .GT. 0 &
							 .AND. work_N_cell2d(i+1,j+1) .GT. 0)

						current => cellHead(i+1,j+1)%p

						IF(debug .EQ. 1) THEN
							PRINT*,'Tracers initially occupying this cell:'
							PRINT*,'Tracer # ',current%id
							DO WHILE(ASSOCIATED(current%next))
								current => current%next
								PRINT*,'Tracer # ',current%id
							END DO
							current => cellHead(i+1,j+1)%p
						END IF
 

						SELECT CASE(face_count)

						CASE(0)	
							! Outgoing through x+ face
							IF(SIGN(1.0,work_flux2d(i+1,j+1,1)) .GT. 0.0) THEN
                		
								m_frac = work_flux2d(i+1,j+1,1)*dt*dy/m0

								IF(debug .EQ. 1) PRINT*,'OUTGOING x+: m0, m_frac=',m0,m_frac 

								m0=m0*(1.0-m_frac)

								! Loop over each tracer originally in the cell
								DO l = 1, work_N_cell2d(i+1,j+1)

									CALL RANDOM_NUMBER(p)
                  
									IF(debug .EQ. 1) THEN
										PRINT*,' Tracer #, p, m_frac = ',current%id,p,m_frac
									END IF

									IF(p .LT. m_frac)THEN

										IF(debug .EQ. 1) THEN
											PRINT*,' Moving Tracer # ',current%id,' into cell ',i+1,j
										END IF

										current%x = i+1
										current => current%next
										CALL SELF%cell_transfer_2d(tr,cellHead,work_N_cell_temp2d(i+1,j+1),&
 												   	 work_N_cell_temp2d(i+2,j+1),l-move_count,i,i+1,j,j)

										work_N_cell_temp2d(i+1,j+1) = work_N_cell_temp2d(i+1,j+1)-1 
										work_N_cell_temp2d(i+2,j+1) = work_N_cell_temp2d(i+2,j+1)+1
										move_count = move_count+1
									ELSE
										current => current%next
									END IF

								END DO

                                ntr= ntr + move_count
								! reset counter
								move_count = 0
							END IF			
	
						CASE(1)
							! Outgoing through x- face
							IF(SIGN(1.0,work_flux2d(i,j+1,1)) .LT. 0.0) THEN

								m_frac = abs(work_flux2d(i,j+1,1))*dt*dy/m0

								IF(debug .EQ. 1) THEN
									PRINT*,'OUTGOING x-: m0, m_frac=',m0,m_frac
								END IF

								current => cellHead(i+1,j+1)%p

								! Loop over each tracer in the cell  
								DO l = 1, work_N_cell2d(i+1,j+1)-ntr

									CALL RANDOM_NUMBER(p)

									IF(debug .EQ. 1) THEN
										PRINT*,' Tracer #, p, m_frac = ',current%id,p,m_frac
									END IF

									IF(p .LT. m_frac)THEN

										IF(debug .EQ. 1) THEN
											PRINT*,' Moving Tracer # ',current%id,' into cell ',i-1, j 
										END IF

										current%x = i-1
										current => current%next

										CALL SELF%cell_transfer_2d(tr,cellHead, work_N_cell_temp2d(i+1,j+1), &
                                                     work_N_cell_temp2d(i,j+1),l-move_count,i,i-1,j,j)

										work_N_cell_temp2d(i+1,j+1) = work_N_cell_temp2d(i+1,j+1)-1 
										work_N_cell_temp2d(i,j+1) = work_N_cell_temp2d(i,j+1)+1
										move_count = move_count+1                                        
									ELSE
										current => current%next
									END IF

								END DO
  
                                ntr = ntr+move_count
								! reset counter
								move_count = 0
							END IF

						CASE(2)
							! Outgoing through y+ face
							IF(SIGN(1.0,work_flux2d(i+1,j+1,2)) .GT. 0.0) THEN
                		
								m_frac = work_flux2d(i+1,j+1,2)*dt*dx/m0

								IF(debug .EQ. 1) PRINT*,'OUTGOING y+: m0, m_frac=',m0,m_frac 

								current => cellHead(i+1,j+1)%p
								m0=m0*(1.0-m_frac)

								! Loop over each tracer originally in the cell
								DO l = 1, work_N_cell2d(i+1,j+1)-ntr

									CALL RANDOM_NUMBER(p)
                             
									IF(debug .EQ. 1) THEN
										PRINT*,' Tracer #, p, m_frac = ',current%id,p,m_frac
									END IF

									IF(p .LT. m_frac)THEN

										IF(debug .EQ. 1) THEN
											PRINT*,' Moving Tracer # ',current%id,' into cell ',i,j+1
										END IF

										current%y = j+1
										current => current%next

										CALL SELF%cell_transfer_2d(tr,cellHead,work_N_cell_temp2d(i+1,j+1),&
 												   	 work_N_cell_temp2d(i+1,j+2),l-move_count,i,i,j,j+1)

										work_N_cell_temp2d(i+1,j+1) = work_N_cell_temp2d(i+1,j+1)-1 
										work_N_cell_temp2d(i+1,j+2) = work_N_cell_temp2d(i+1,j+2)+1
										move_count = move_count+1
									ELSE
										current => current%next
									END IF

								END DO

                                ntr = ntr+move_count 
								! reset counter
								move_count = 0
							END IF		

						CASE(3)
							! Outgoing through y- face
							IF(SIGN(1.0,work_flux2d(i+1,j,2)) .LT. 0.0) THEN
                		
								m_frac = abs(work_flux2d(i+1,j,2))*dt*dx/m0

								IF(debug .EQ. 1) PRINT*,'OUTGOING y-: m0, m_frac=',m0,m_frac 

								current => cellHead(i+1,j+1)%p

								m0=m0*(1.0-m_frac)

								! Loop over each tracer originally in the cell
								DO l = 1, work_N_cell2d(i+1,j+1)-ntr

									CALL RANDOM_NUMBER(p)
                             
									IF(debug .EQ. 1) THEN
										PRINT*,' Tracer #, p, m_frac = ',current%id,p,m_frac
									END IF

									IF(p .LT. m_frac)THEN

										IF(debug .EQ. 1) THEN
											PRINT*,' Moving Tracer # ',current%id,' into cell ',i,j-1
										END IF

										current%y = j-1
										current => current%next

										CALL SELF%cell_transfer_2d(tr,cellHead,work_N_cell_temp2d(i+1,j+1),&
 												   	 work_N_cell_temp2d(i+1,j),l-move_count,i,i,j,j-1)

										work_N_cell_temp2d(i+1,j+1) = work_N_cell_temp2d(i+1,j+1)-1 
										work_N_cell_temp2d(i+1,j) = work_N_cell_temp2d(i+1,j)+1
										move_count = move_count+1
									ELSE
										current => current%next
									END IF

								END DO
							END IF			

						END SELECT

                        face_count = face_count +1

                    END DO

				END DO
			END DO		
			
			! update N_cell array
            N_cell(:,:) = work_N_cell_temp2d(:,:)


		END SUBROUTINE advect_2d



        ! routine for performing a single tracer exchange between two cells on 2d grid

		SUBROUTINE cell_transfer_2d(SELF, tr, cellHead, N_src, N_des, tr_num ,x_src, x_des, y_src, y_des)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
			TYPE(tracer_ptr_2d), INTENT(INOUT) :: tr(:),cellHead(:,:)
			INTEGER, INTENT(IN) :: tr_num, N_src, N_des, x_src, x_des, y_src, y_des 
                                                        ! tr_num: tracer position in linked list 
                                                        ! N_src : num of tracers in source cell
                                                        ! N_des : num of tracers in destination cell
														! x_src, y_src : source cell indices
														! x_des, y_des : destination cell indices	
							
            INTEGER :: i
            TYPE(tracer_2d), POINTER :: current_src, current_des, temp, &
										current_src_above => null(), current_src_below => null() 


            ! clear temp pointers
            current_src => null()
            current_src_above => null()
            current_src_below => null()
            current_des => null()
            temp => null()



            IF(debug .EQ. 1) THEN
                temp => cellHead(x_src+1, y_src+1)%p 
				PRINT*,'Tracers currently in source cell:'
				IF(ASSOCIATED(temp)) THEN
					 PRINT*,'Tracer # ',temp%id
    	            DO WHILE(ASSOCIATED(temp%next))
						temp => temp%next
        	      		PRINT*,'Tracer # ',temp%id
					END DO
				END IF

                temp => cellHead(x_des+1,y_des+1)%p 
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
			current_src => cellHead(x_src+1, y_src+1)%p

            ! traverse to the outgoing tracer location within link
			DO i = 1, tr_num-1
                current_src_above => current_src
				current_src => current_src%next 
				!IF(debug .EQ. 1) PRINT*,'Traversing list. Currently &
                !                         pointing at Tracer #',current_src%id                    

       	    END DO
               
            ! set a temp pointer to tracer next in list after the outgoing tracer
            IF(ASSOCIATED(current_src%next)) current_src_below => current_src%next

			! set a temp pointer to head tracer in destination cell
			current_des => cellHead(x_des+1,y_des+1)%p

			! if destination cell empty, then incoming tracer is assigned as head	
            IF(.NOT. ASSOCIATED(current_des))THEN
				cellHead(x_des+1,y_des+1)%p => current_src
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
			! 				 (once-sided link re-attach)
			!
			!		      3) If outgoing tracer is neither head nor tail, then
			!				 link the tracer that is one above it in the source cell
			!				 to the tracer that is one below it.
			!				 (two-sided link re-attach)
            !#####################################################################

			IF(tr_num .EQ. 1) THEN	

				IF(debug .EQ. 1) PRINT*,'Case 1. Moving Head tracer.'

				! if more than 1 tracer in cell, promote next tracer to head
                IF(ASSOCIATED(current_src_below))THEN
                	cellHead(x_src+1,y_src+1)%p => current_src_below
                    IF(debug .EQ. 1) PRINT*,'Tracer # ', current_src_below%id ,' promoted to head.'
				! otherwise destroy head tracer pointer
                ELSE
                	NULLIFY(cellHead(x_src+1,y_src+1)%p)
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
                temp => cellHead(x_src+1,y_src+1)%p 
				PRINT*,'Tracers currently in source cell:'
				IF(ASSOCIATED(temp)) THEN 
					PRINT*,'Tracer # ',temp%id
                	DO WHILE(ASSOCIATED(temp%next))
						temp => temp%next
              			PRINT*,'Tracer # ',temp%id
					END DO
				END IF
                temp => cellHead(x_des+1,y_des+1)%p 
				PRINT*,'Tracers currently in destination cell:'
				IF(ASSOCIATED(temp)) THEN
				PRINT*,'Tracer # ',temp%id
	                DO WHILE(ASSOCIATED(temp%next))
						temp => temp%next
	              		PRINT*,'Tracer # ',temp%id
					END DO
				END IF
	        END IF
  

		END SUBROUTINE cell_transfer_2d



        ! Routine for 1d tracer advection
		!
		! Loop over cells. Cycle through the cell faces (in the order x+, 
        ! x-). For each face: 
		!
		! 1) Compute outgoing fluid mass fractions (m_frac).
		!
		! e.g. through right face, m_frac = mass_flux*dt/m0
        !      where m0 is the mass of fluid currently in the cell (which is
		!	   updated after each cell face cycle: m0_new = m0_old*(1-m_frac) ).
		!
        ! 2) Use Monte Carlo sampling to determine whether or not to move a tracer
		! 	 through each cell face. Cycle through every tracer currently occupying  
		!    that cell. For each tracer, draw a random number p from [0,1).  
		!    If p < m_frac, then move that tracer, otherwise don't.

		SUBROUTINE advect_1d(SELF, tr, cellHead, N_cell, flux, rho)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF 
			TYPE(tracer_ptr_1d), INTENT(INOUT) :: tr(:),cellHead(:)
            INTEGER, INTENT(INOUT) :: N_cell(:)
            REAL, INTENT(IN) :: flux(:), rho(:)
			
			INTEGER :: nx           
			TYPE(tracer_1d), POINTER :: current
     		INTEGER :: i, l, ntr, move_count, face_count
            REAL :: p, m_frac, m0


            ! clear temp pointers
            current => null()

			nx = SELF%nx


            ! copy values into work pool arrays
            work_N_cell(:) = N_cell(:)
            work_N_cell_temp(:) = N_cell(:)
			work_flux(:) = flux(:)
            work_rho(:) = rho(:) 
            
			! loop across cells in our 1d grid

            DO i = 1, nx

               ! reset counters
				ntr= 0
				move_count = 0
				face_count = 0

                ! calculate fluid mass in cell
                m0=work_rho(i)*dx
              
                IF(debug .EQ. 1) THEN
                    PRINT*, ' '
					PRINT*,'Cell #, N_cell:',i, work_N_cell(i+1)
					PRINT*,'F_L, F_R = ', work_flux(i), work_flux(i+1)
                    PRINT*, ' '
                END IF
				  
                ! cycle through both cell faces
				DO WHILE(face_count .LT. 2 .AND. work_N_cell_temp(i+1) .GT. 0 &
					 .AND. work_N_cell(i+1) .GT. 0)

                    current => cellHead(i+1)%p
			       	
                    IF(debug .EQ. 1) THEN
                        PRINT*,'face_index=',face_count+1
                        PRINT*,'Tracers initially occupying this cell:'
						PRINT*,'Tracer # ',current%id
                        DO WHILE(ASSOCIATED(current%next))
							current => current%next
                            PRINT*,'Tracer # ',current%id
						END DO
                        current => cellHead(i+1)%p
	        	    END IF
 
                    SELECT CASE(face_count)

                    CASE(0)
                        ! Outgoing through x+ face
                        IF(SIGN(1.0,work_flux(i+1)) .GT. 0.0) THEN
                		
                            m_frac = work_flux(i+1)*dt/m0

                            IF(debug .EQ. 1) PRINT*,'OUTGOING x+: m0, m_frac=',m0,m_frac 

                            m0=m0*(1.0-m_frac)

                            ! Loop over each tracer originally in the cell
                            DO l = 1, work_N_cell(i+1)

                                CALL RANDOM_NUMBER(p)
                             
                                IF(debug .EQ. 1) THEN
                                    PRINT*,' Tracer #, p, m_frac = ',current%id,p,m_frac
                                END IF

                                IF(p .LT. m_frac)THEN

                                    IF(debug .EQ. 1) THEN
                                        PRINT*,' Moving Tracer # ',current%id,' into cell ',i+1
                                    END IF

                                    current%x=i+1
                                    current => current%next

                                    CALL SELF%cell_transfer_1d(tr,cellHead,work_N_cell_temp(i+1),&
 												   	 work_N_cell_temp(i+2),l-move_count,i,i+1)

                                    work_N_cell_temp(i+1) = work_N_cell_temp(i+1)-1 
                                    work_N_cell_temp(i+2) = work_N_cell_temp(i+2)+1
                                    move_count = move_count+1
                                ELSE
                                    current => current%next
                                END IF

                            END DO

                            ntr = ntr+move_count
                            ! reset counter
                            move_count = 0

                        END IF			

                    CASE(1)
                        ! Outgoing through x- face
                        IF(SIGN(1.0,work_flux(i)) .LT. 0.0) THEN

                            m_frac = abs(work_flux(i))*dt/m0

                            IF(debug .EQ. 1) THEN
                                PRINT*,'OUTGOING x-: m0, m_frac=',m0,m_frac
                            END IF

                            current => cellHead(i+1)%p

                            ! Loop over each tracer in the cell  
                            DO l = 1, work_N_cell(i+1)-ntr
                            
                                CALL RANDOM_NUMBER(p)

                                IF(debug .EQ. 1) THEN
                                    PRINT*,' Tracer #, p, m_frac = ',current%id,p,m_frac
                                END IF

                                IF(p .LT. m_frac)THEN

                                    IF(debug .EQ. 1) THEN
                                        PRINT*,' Moving Tracer # ',current%id,' into cell ',i-1 
                                    END IF

                                    current%x = i-1
                                    current => current%next

                                    CALL SELF%cell_transfer_1d(tr,cellHead, work_N_cell_temp(i+1), &
                                                     work_N_cell_temp(i),l-move_count,i,i-1)

                                    work_N_cell_temp(i+1) = work_N_cell_temp(i+1)-1 
                                    work_N_cell_temp(i) = work_N_cell_temp(i)+1
                                    move_count = move_count+1
                                ELSE
                                    current => current%next
                                END IF
                            END DO
                        END IF

                    END SELECT

                    face_count = face_count +1

                END DO

			END DO		
			
			! update N_cell array
            N_cell(:) = work_N_cell_temp(:)


		END SUBROUTINE advect_1d



        ! routine for performing a single tracer exchange between two cells on 1d grid

		SUBROUTINE cell_transfer_1d(SELF, tr, cellHead, N_src, N_des, tr_num ,x_src, x_des)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
			TYPE(tracer_ptr_1d), INTENT(INOUT) :: tr(:),cellHead(:)
			INTEGER, INTENT(IN) :: tr_num, N_src, N_des, x_src, x_des ! tr_num: tracer number in linked list 
														! x_src: source cell
														! x_des: destination cell								
            INTEGER :: i
            TYPE(tracer_1d), POINTER :: current_src, current_des, temp, &
										current_src_above, current_src_below


            ! clear temp pointers
            current_src => null()
            current_src_above => null()
            current_src_below => null()
            current_des => null()
            temp => null()

            IF(debug .EQ. 1) THEN
                temp => cellHead(x_src+1)%p 
				PRINT*,'Tracers currently in source cell:'
				IF(ASSOCIATED(temp)) THEN
					 PRINT*,'Tracer # ',temp%id
    	            DO WHILE(ASSOCIATED(temp%next))
						temp => temp%next
        	      		PRINT*,'Tracer # ',temp%id
					END DO
				END IF

                temp => cellHead(x_des+1)%p 
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
			current_src => cellHead(x_src+1)%p

    
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
			current_des => cellHead(x_des+1)%p

			! if destination cell empty, then incoming tracer is assigned as head	
            IF(.NOT. ASSOCIATED(current_des))THEN
				cellHead(x_des+1)%p => current_src
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
                	cellHead(x_src+1)%p => current_src_below
                    IF(debug .EQ. 1) PRINT*,'Tracer # ', current_src_below%id ,' promoted to head.'
				! otherwise destroy head tracer pointer
                ELSE
                	NULLIFY(cellHead(x_src+1)%p)
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
                temp => cellHead(x_src+1)%p 
				PRINT*,'Tracers currently in source cell:'
				IF(ASSOCIATED(temp)) THEN 
					PRINT*,'Tracer # ',temp%id
                	DO WHILE(ASSOCIATED(temp%next))
						temp => temp%next
              			PRINT*,'Tracer # ',temp%id
					END DO
				END IF
                temp => cellHead(x_des+1)%p 
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

