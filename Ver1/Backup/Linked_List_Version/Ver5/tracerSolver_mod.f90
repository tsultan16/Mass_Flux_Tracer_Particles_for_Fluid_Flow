MODULE tracersolver_mod
	USE tracertype_mod
    USE data_mod
	IMPLICIT NONE

	! tracer advection solver object

    PUBLIC :: MASSFLUXTRACER

    INTEGER, PRIVATE :: ndims, nx, ny, nz, N    
    REAL, PRIVATE :: dt, dx, dy, dz

    REAL,PRIVATE :: ti,tf,t_ex!cell exchange time

    ! work pool variables 
    REAL, ALLOCATABLE, PRIVATE :: work_flux(:), work_rho(:), &
                                  work_flux2d(:,:,:), work_rho2d(:,:), &
                                  work_flux3d(:,:,:,:), work_rho3d(:,:,:)

    INTEGER, ALLOCATABLE, PRIVATE :: work_N_cell1d(:), work_N_cell_temp(:), &
                                     work_N_cell2d(:,:), work_N_cell_temp2d(:,:), & 
                                     work_N_cell3d(:,:,:), work_N_cell_temp3d(:,:,:)
        
    TYPE (tracer_ptr_1d), ALLOCATABLE, PRIVATE :: work_cellHead1d(:)
    TYPE (tracer_ptr_2d), ALLOCATABLE, PRIVATE :: work_cellHead2d(:,:)
    TYPE (tracer_ptr_3d), ALLOCATABLE, PRIVATE :: work_cellHead3d(:,:,:)


	TYPE :: MASSFLUXTRACER
        ! type data

        ! type bound procedures		
		CONTAINS

        PROCEDURE, PASS(SELF), PUBLIC :: initialize_workpool
        PROCEDURE, PASS(SELF), PUBLIC :: destroy_workpool
        PROCEDURE, PASS(SELF), PUBLIC :: solve

		PROCEDURE, PASS(SELF), PRIVATE :: advect_1d
		PROCEDURE, PASS(SELF), PRIVATE :: cell_transfer_1d

		PROCEDURE, PASS(SELF), PRIVATE :: advect_2d
		PROCEDURE, PASS(SELF), PRIVATE :: cell_transfer_2d

		PROCEDURE, PASS(SELF), PRIVATE :: advect_3d
		PROCEDURE, PASS(SELF), PRIVATE :: cell_transfer_3d

    END TYPE MASSFLUXTRACER

	INTERFACE MASSFLUXTRACER
		MODULE PROCEDURE constructor
	END INTERFACE 

    
	CONTAINS

		FUNCTION constructor(ndims_in,nx_in,ny_in,nz_in,N_in,dx_in,dy_in,dz_in) RESULT(SELF)

			TYPE(MASSFLUXTRACER) :: SELF
            INTEGER, INTENT(IN) :: ndims_in, nx_in, ny_in, nz_in, N_in
            REAL, INTENT(IN) :: dx_in, dy_in, dz_in

			ndims = ndims_in
			nx = nx_in
			ny = ny_in
			nz = nz_in
            N = N_in
            dx = dx_in
			dy = dy_in
			dz = dz_in
            
			RETURN

		END FUNCTION constructor

        SUBROUTINE initialize_workpool(SELF)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
			
			IF(ndims .EQ. 1) THEN
                ALLOCATE(work_rho(1:nx))
				ALLOCATE(work_flux(1:nx+1))
				ALLOCATE(work_N_cell1d(1:nx+2),work_N_cell_temp(1:nx+2))
                ALLOCATE(work_cellHead1d(1:nx+2))

				work_rho(:) = 0.0
				work_flux(:) = 0.0
				work_N_cell1d(:) = 0
				work_N_cell_temp(:) = 0	  

			END IF

            IF(ndims .EQ. 2) THEN
                ALLOCATE(work_rho2d(1:nx,1:ny))
				ALLOCATE(work_flux2d(1:nx+1,1:ny+1,2))
				ALLOCATE(work_N_cell2d(1:nx+2,1:ny+2),work_N_cell_temp2d(1:nx+2,1:ny+2))
                ALLOCATE(work_cellHead2d(1:nx+2,1:ny+2))

				work_rho2d(:,:) = 0.0
				work_flux2d(:,:,:) = 0.0
				work_N_cell2d(:,:) = 0
				work_N_cell_temp2d(:,:) = 0

			END IF

            IF(ndims .EQ. 3) THEN
				ALLOCATE(work_rho3d(1:nx,1:ny,1:nz))
				ALLOCATE(work_flux3d(1:nx+1,1:ny+1,1:nz+1,3))
				ALLOCATE(work_N_cell3d(1:nx+2,1:ny+2,1:nz+2),work_N_cell_temp3d(1:nx+2,1:ny+2,1:nz+2))
                ALLOCATE(work_cellHead3d(1:nx+2,1:ny+2,1:nz+2))

				work_rho3d(:,:,:) = 0.0
				work_flux3d(:,:,:,:) = 0.0
				work_N_cell3d(:,:,:) = 0
				work_N_cell_temp3d(:,:,:) = 0	  

			END IF
	

		END SUBROUTINE initialize_workpool


		SUBROUTINE destroy_workpool(SELF)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF

            IF(ndims .EQ. 1) THEN

				DEALLOCATE(work_rho, work_flux, work_N_cell1d, &
                           work_N_cell_temp, work_cellHead1d)

			ELSE IF(ndims .EQ. 2) THEN

				DEALLOCATE(work_rho2d, work_flux2d,work_N_cell2d, &
                           work_N_cell_temp2d, work_cellHead2d)

            ELSE IF(ndims .EQ. 3) THEN

				DEALLOCATE(work_rho3d, work_flux3d, work_N_cell3d, &
                           work_N_cell_temp3d, work_cellHead3d)

			END IF


		END SUBROUTINE destroy_workpool


		! Top-level routine for mass-flux tracer advection

		SUBROUTINE solve(SELF, dt_in)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
            REAL, INTENT(IN) :: dt_in
			
            dt = dt_in

            SELECT CASE(ndims) 
            
            CASE(1) 
				CALL SELF%advect_1d(cellHead1d,N_cell1d,flux,rho)
			CASE(2)
                CALL SELF%advect_2d(cellHead2d,N_cell2d,flux2d,rho2d)
            CASE(3)
                CALL SELF%advect_3d(cellHead3d,N_cell3d,flux3d,rho3d)

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

		SUBROUTINE advect_3d(SELF, cellHead, N_cell, flux, rho)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF 
			TYPE(tracer_ptr_3d), INTENT(INOUT) :: cellHead(:,:,:)
            INTEGER, INTENT(INOUT) :: N_cell(:,:,:)
            REAL, INTENT(IN) :: flux(:,:,:,:), rho(:,:,:)
			           
			TYPE(tracer_3d), POINTER :: current
     		INTEGER :: i, j, k, l, ntr, move_count, face_count
            REAL :: m0, m_frac, p
         

            ! clear temp pointer
            current => null()

            ! copy values into work pool arrays
            work_cellHead3d(:,:,:) = cellHead(:,:,:)
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
                    
				    ! cycle through all six cell faces
					DO WHILE(face_count .LT. 6 .AND. work_N_cell_temp3d(i+1,j+1,k+1) .GT. 0 &
							 .AND. work_N_cell3d(i+1,j+1,k+1) .GT. 0)

						SELECT CASE(face_count)

						CASE(0)	
							! Outgoing through x+ face
							IF(SIGN(1.0,work_flux3d(i+1,j+1,k+1,1)) .GT. 0.0) THEN

                                current => work_cellHead3d(i+1,j+1,k+1)%p                		
								m_frac = work_flux3d(i+1,j+1,k+1,1)*dt*dy*dz/m0
								m0=m0*(1.0-m_frac)

								! Loop over each tracer originally in the cell
								DO l = 1, work_N_cell3d(i+1,j+1,k+1)

									CALL RANDOM_NUMBER(p)
                  
									IF(p .LT. m_frac)THEN

                                        current%x = i+1
										current => current%next
										CALL SELF%cell_transfer_3d(work_cellHead3d,work_N_cell_temp3d(i+1,j+1,k+1),&
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

								current => work_cellHead3d(i+1,j+1,k+1)%p
								m_frac = abs(work_flux3d(i,j+1,k+1,1))*dt*dy*dz/m0
								m0=m0*(1.0-m_frac)

								! Loop over each tracer in the cell  
								DO l = 1, work_N_cell3d(i+1,j+1,k+1)-ntr

									CALL RANDOM_NUMBER(p)
							
									IF(p .LT. m_frac)THEN

										current%x = i-1
										current => current%next

										CALL SELF%cell_transfer_3d(work_cellHead3d, work_N_cell_temp3d(i+1,j+1,k+1), &
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

								current => work_cellHead3d(i+1,j+1,k+1)%p                		
								m_frac = work_flux3d(i+1,j+1,k+1,2)*dt*dx*dz/m0
								m0=m0*(1.0-m_frac)

								! Loop over each tracer originally in the cell
								DO l = 1, work_N_cell3d(i+1,j+1,k+1)-ntr

									CALL RANDOM_NUMBER(p)
                             
									IF(p .LT. m_frac)THEN

										current%y = j+1
										current => current%next

										CALL SELF%cell_transfer_3d(work_cellHead3d,work_N_cell_temp3d(i+1,j+1,k+1),&
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
                		
								current => work_cellHead3d(i+1,j+1,k+1)%p
								m_frac = abs(work_flux3d(i+1,j,k+1,2))*dt*dx*dz/m0
								m0=m0*(1.0-m_frac)

								! Loop over each tracer originally in the cell
								DO l = 1, work_N_cell3d(i+1,j+1,k+1)-ntr

									CALL RANDOM_NUMBER(p)
                             
									IF(p .LT. m_frac)THEN

										current%y = j-1
										current => current%next

										CALL SELF%cell_transfer_3d(work_cellHead3d,work_N_cell_temp3d(i+1,j+1,k+1),&
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
                		
                                current => work_cellHead3d(i+1,j+1,k+1)%p
								m_frac = work_flux3d(i+1,j+1,k+1,3)*dt*dx*dy/m0
								m0=m0*(1.0-m_frac)

								! Loop over each tracer originally in the cell
								DO l = 1, work_N_cell3d(i+1,j+1,k+1)-ntr

									CALL RANDOM_NUMBER(p)
                         
									IF(p .LT. m_frac)THEN

                                        current%z = k+1
										current => current%next

										CALL SELF%cell_transfer_3d(work_cellHead3d,work_N_cell_temp3d(i+1,j+1,k+1),&
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
                		
								current => work_cellHead3d(i+1,j+1,k+1)%p
								m_frac = abs(work_flux3d(i+1,j+1,k,2))*dt*dx*dy/m0
								m0=m0*(1.0-m_frac)

								! Loop over each tracer originally in the cell
								DO l = 1, work_N_cell3d(i+1,j+1,k+1)-ntr

									CALL RANDOM_NUMBER(p)

									IF(p .LT. m_frac)THEN

										current%z = k-1
										current => current%next

										CALL SELF%cell_transfer_3d(work_cellHead3d,work_N_cell_temp3d(i+1,j+1,k+1),&
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
			
			! copy work array data into storage arrays
            N_cell(:,:,:) = work_N_cell_temp3d(:,:,:)
            cellHead(:,:,:) = work_cellHead3d(:,:,:)

		END SUBROUTINE advect_3d



        ! routine for performing a single tracer exchange between two cells on 3d grid

		SUBROUTINE cell_transfer_3d(SELF, cellHead, N_src, N_des, tr_num ,x_src, x_des, &
                                    y_src, y_des, z_src, z_des)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
			TYPE(tracer_ptr_3d), INTENT(INOUT) :: cellHead(:,:,:)
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
       	    END DO
               
            ! set a temp pointer to tracer next in list after the outgoing tracer
            IF(ASSOCIATED(current_src%next)) current_src_below => current_src%next

			! set a temp pointer to head tracer in destination cell
			current_des => cellHead(x_des+1,y_des+1,z_des+1)%p

			! if destination cell empty, then incoming tracer is assigned as head	
            IF(.NOT. ASSOCIATED(current_des))THEN
				cellHead(x_des+1,y_des+1,z_des+1)%p => current_src
				NULLIFY(current_src%next)

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
			!				 (one-sided re-attach)
			!
			!			  2) If outgoing tracer is tail of source cell, then 
			! 				 nulify next pointer of second to last tracer in source
			!                cell.
			! 				 (once-sided re-attach)
			!
			!		      3) If outgoing tracer is neither head nor tail, then
			!				 link the tracer that is one above it in the source cell
			!				 to the tracer that is one below it.
			!				 (two-sided re-attach)
            !#####################################################################

			IF(tr_num .EQ. 1) THEN	
				! if more than 1 tracer in cell, promote next tracer to head
                IF(ASSOCIATED(current_src_below))THEN
                	cellHead(x_src+1,y_src+1,z_src+1)%p => current_src_below
				! otherwise destroy head tracer pointer
                ELSE
                	NULLIFY(cellHead(x_src+1,y_src+1,z_src+1)%p)
				END IF

			END IF

			IF(tr_num .GT. 1 .AND. tr_num .EQ. N_src) THEN
				NULLIFY(current_src_above%next)
			END IF
			
			IF(tr_num .GT. 1 .AND. tr_num .LT. N_src) THEN
				current_src_above%next => current_src_below
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

		SUBROUTINE advect_2d(SELF, cellHead, N_cell, flux, rho)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF 
			TYPE(tracer_ptr_2d), INTENT(INOUT) :: cellHead(:,:)
            INTEGER, INTENT(INOUT) :: N_cell(:,:)
            REAL, INTENT(IN) :: flux(:,:,:), rho(:,:)
			           
			TYPE(tracer_2d), POINTER :: current
     		INTEGER :: i, j, l, ntr, move_count, face_count
            REAL :: m0, m_frac, p

            ! clear temp pointer
            current => null()

            ! copy values into work pool arrays
            work_N_cell2d(:,:) = N_cell(:,:)
            work_N_cell_temp2d(:,:) = N_cell(:,:)
			work_flux2d(:,:,:) = flux(:,:,:)
            work_rho2d(:,:) = rho(:,:) 
            work_cellHead2d(:,:) = cellHead(:,:)

			! loop across cells in our 2d grid
            DO i = 1, nx
				DO j = 1, ny
		  
                PRINT*,'CELL#',i,j


                    ! reset counters
                    ntr= 0
                    move_count = 0
                    face_count = 0

                    ! calculate initial fluid mass in cell
                    m0 = work_rho2d(i,j)*dx*dy

				    ! cycle through all four cell faces
					DO WHILE(face_count .LT. 4 .AND. work_N_cell_temp2d(i+1,j+1) .GT. 0 &
							 .AND. work_N_cell2d(i+1,j+1) .GT. 0)

     					SELECT CASE(face_count)

						CASE(0)	
							! Outgoing through x+ face
							IF(SIGN(1.0,work_flux2d(i+1,j+1,1)) .GT. 0.0) THEN


                                CALL CPU_TIME(ti)                                


                                current => work_cellHead2d(i+1,j+1)%p
                                m_frac = work_flux2d(i+1,j+1,1)*dt*dy/m0
                                m0=m0*(1.0-m_frac)

								! Loop over each tracer originally in the cell
								DO l = 1, work_N_cell2d(i+1,j+1)

									CALL RANDOM_NUMBER(p)
                                    
                                    IF(p .LT. m_frac)THEN
                        
                                        current%x = i+1
										current => current%next

										CALL SELF%cell_transfer_2d(work_cellHead2d,work_N_cell_temp2d(i+1,j+1),&
 												   	 work_N_cell_temp2d(i+2,j+1),l-move_count,i,i+1,j,j)

										work_N_cell_temp2d(i+1,j+1) = work_N_cell_temp2d(i+1,j+1)-1 
										work_N_cell_temp2d(i+2,j+1) = work_N_cell_temp2d(i+2,j+1)+1
										move_count = move_count+1

									ELSE
										current => current%next
									END IF

								END DO

                                CALL CPU_TIME(tf)
                                t_ex=tf-ti
                                PRINT*,'Moving ', move_count,' tracer between cells', i,j,i+1,j 
                                PRINT*,'Time taken=',t_ex


                                ntr= ntr + move_count
								! reset counter
								move_count = 0
							END IF			
	
						CASE(1)
							! Outgoing through x- face
							IF(SIGN(1.0,work_flux2d(i,j+1,1)) .LT. 0.0) THEN

                                current => work_cellHead2d(i+1,j+1)%p
								m_frac = abs(work_flux2d(i,j+1,1))*dt*dy/m0
                                m0=m0*(1.0-m_frac)
								
								! Loop over each tracer in the cell  
								DO l = 1, work_N_cell2d(i+1,j+1)-ntr

									CALL RANDOM_NUMBER(p)

                                    IF(p .LT. m_frac)THEN

                                        current%x = i-1
										current => current%next

										CALL SELF%cell_transfer_2d(work_cellHead2d, work_N_cell_temp2d(i+1,j+1), &
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
                		
								current => work_cellHead2d(i+1,j+1)%p
								m_frac = work_flux2d(i+1,j+1,2)*dt*dx/m0
								m0=m0*(1.0-m_frac)

								! Loop over each tracer originally in the cell
								DO l = 1, work_N_cell2d(i+1,j+1)-ntr

									CALL RANDOM_NUMBER(p)
                             
									IF(p .LT. m_frac)THEN

										current%y = j+1
										current => current%next

										CALL SELF%cell_transfer_2d(work_cellHead2d,work_N_cell_temp2d(i+1,j+1),&
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
                		
								current => work_cellHead2d(i+1,j+1)%p
								m_frac = abs(work_flux2d(i+1,j,2))*dt*dx/m0
								m0=m0*(1.0-m_frac)

								! Loop over each tracer originally in the cell
								DO l = 1, work_N_cell2d(i+1,j+1)-ntr

									CALL RANDOM_NUMBER(p)
                             
									IF(p .LT. m_frac)THEN

                                        current%y = j-1
										current => current%next

										CALL SELF%cell_transfer_2d(work_cellHead2d,work_N_cell_temp2d(i+1,j+1),&
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
			
			!copy work array data into storage arrays
            N_cell(:,:) = work_N_cell_temp2d(:,:)
            cellHead(:,:) = work_cellHead2d(:,:)

		END SUBROUTINE advect_2d



        ! routine for performing a single tracer exchange between two cells on 2d grid

		SUBROUTINE cell_transfer_2d(SELF, cellHead, N_src, N_des, tr_num ,x_src, x_des, y_src, y_des)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
			TYPE(tracer_ptr_2d), INTENT(INOUT) :: cellHead(:,:)
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
       	    END DO
               
            ! set a temp pointer to tracer next in list after the outgoing tracer
            IF(ASSOCIATED(current_src%next)) current_src_below => current_src%next

			! set a temp pointer to head tracer in destination cell
			current_des => cellHead(x_des+1,y_des+1)%p

			! if destination cell empty, then incoming tracer is assigned as head	
            IF(.NOT. ASSOCIATED(current_des))THEN
				cellHead(x_des+1,y_des+1)%p => current_src
				NULLIFY(current_src%next)
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
				! if more than 1 tracer in cell, promote next tracer to head
                IF(ASSOCIATED(current_src_below))THEN
                	cellHead(x_src+1,y_src+1)%p => current_src_below
				! otherwise destroy head tracer pointer
                ELSE
                	NULLIFY(cellHead(x_src+1,y_src+1)%p)
				END IF

			END IF

			IF(tr_num .GT. 1 .AND. tr_num .EQ. N_src) THEN
				NULLIFY(current_src_above%next)
			END IF
			
			IF(tr_num .GT. 1 .AND. tr_num .LT. N_src) THEN
				current_src_above%next => current_src_below
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

		SUBROUTINE advect_1d(SELF, cellHead, N_cell, flux, rho)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF 
			TYPE(tracer_ptr_1d), INTENT(INOUT) :: cellHead(:)
            INTEGER, INTENT(INOUT) :: N_cell(:)
            REAL, INTENT(IN) :: flux(:), rho(:)
			
			TYPE(tracer_1d), POINTER :: current
     		INTEGER :: i, l, ntr, move_count, face_count
            REAL :: p, m_frac, m0

            ! clear temp pointers
            current => null()

	        ! copy values into work pool arrays
            work_N_cell1d(:) = N_cell(:)
            work_N_cell_temp(:) = N_cell(:)
			work_flux(:) = flux(:)
            work_rho(:) = rho(:) 
            work_cellHead1d(:) = cellHead(:)
 
			! loop across cells in our 1d grid

            DO i = 1, nx

                ! reset counters
                ntr= 0
                move_count = 0
                face_count = 0

                ! calculate fluid mass in cell
                m0=work_rho(i)*dx
 				  
                ! cycle through both cell faces
				DO WHILE(face_count .LT. 2 .AND. work_N_cell_temp(i+1) .GT. 0 &
					 .AND. work_N_cell1d(i+1) .GT. 0)
             
                    SELECT CASE(face_count)

                    CASE(0)
                        ! Outgoing through x+ face
                        IF(SIGN(1.0,work_flux(i+1)) .GT. 0.0) THEN
                		
                            current => work_cellHead1d(i+1)%p
			                m_frac = work_flux(i+1)*dt/m0
                            m0=m0*(1.0-m_frac)

                            ! Loop over each tracer originally in the cell
                            DO l = 1, work_N_cell1d(i+1)

                                CALL RANDOM_NUMBER(p)
                            
                                IF(p .LT. m_frac)THEN

                                    current%x=i+1
                                    current => current%next

                                    CALL SELF%cell_transfer_1d(work_cellHead1d,work_N_cell_temp(i+1),&
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

                            current => work_cellHead1d(i+1)%p
                            m_frac = abs(work_flux(i))*dt/m0
                      
                            ! Loop over each tracer in the cell  
                            DO l = 1, work_N_cell1d(i+1)-ntr
                            
                                CALL RANDOM_NUMBER(p)

                                IF(p .LT. m_frac)THEN

                                    current%x = i-1
                                    current => current%next

                                    CALL SELF%cell_transfer_1d(work_cellHead1d, work_N_cell_temp(i+1), &
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
			
			! copy work array data into storage arrays
            N_cell(:) = work_N_cell_temp(:)
            cellHead(:) = work_cellHead1d(:)

		END SUBROUTINE advect_1d



        ! routine for performing a single tracer exchange between two cells on 1d grid

		SUBROUTINE cell_transfer_1d(SELF, cellHead, N_src, N_des, tr_num ,x_src, x_des)

			CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
			TYPE(tracer_ptr_1d), INTENT(INOUT) :: cellHead(:)
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
       	    END DO
               
            ! set a temp pointer to tracer next in list after the outgoing tracer
            IF(ASSOCIATED(current_src%next)) current_src_below => current_src%next

			! set a temp pointer to head tracer in destination cell
			current_des => cellHead(x_des+1)%p

			! if destination cell empty, then incoming tracer is assigned as head	
            IF(.NOT. ASSOCIATED(current_des))THEN
				cellHead(x_des+1)%p => current_src
				NULLIFY(current_src%next)
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
                ! if more than 1 tracer in cell, promote next tracer to head
                IF(ASSOCIATED(current_src_below))THEN
                	cellHead(x_src+1)%p => current_src_below
				! otherwise destroy head tracer pointer
                ELSE
                	NULLIFY(cellHead(x_src+1)%p)
				END IF

			END IF

			IF(tr_num .GT. 1 .AND. tr_num .EQ. N_src) THEN
				NULLIFY(current_src_above%next)
			END IF
			
			IF(tr_num .GT. 1 .AND. tr_num .LT. N_src) THEN
				current_src_above%next => current_src_below
			END IF


		END SUBROUTINE cell_transfer_1d



END MODULE tracersolver_mod

