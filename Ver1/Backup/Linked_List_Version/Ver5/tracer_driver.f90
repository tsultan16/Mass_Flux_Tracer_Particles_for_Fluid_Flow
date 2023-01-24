PROGRAM tracer_test_run
    USE constants_mod
	USE tracertype_mod
    USE data_mod
	USE tracersolver_mod
	IMPLICIT NONE

	! define mass flux tracer object
    TYPE(MASSFLUXTRACER) :: solver
    REAL :: dt1, dx1, dy1, dz1
    INTEGER :: i

    REAL :: mfrac(4)
    REAL :: time_start, time_end, time_init, time_sim(nt), time_del


    IF(ndims .EQ. 1) THEN
        OPEN(UNIT=1, FILE='output1d.txt')	
    ELSE IF(ndims .EQ. 2) THEN
        OPEN(UNIT=1, FILE='output2d.txt')
    END IF

	IF(ndims .EQ. 1) THEN
		ALLOCATE(tr1d(1:N))
		ALLOCATE(cellHead1d(0:nx+1))
		ALLOCATE(rho(1:nx), flux(0:nx))
		ALLOCATE(N_cell1d(0:nx+1))
  
        rho = 0.0
        flux = 0.0
        N_cell1d = 0

    ELSE IF(ndims .EQ. 2) THEN
		ALLOCATE(tr2d(1:N))
		ALLOCATE(cellHead2d(0:nx+1,0:ny+1))
		ALLOCATE(rho2d(1:nx,1:ny), flux2d(0:nx,0:ny,2))
		ALLOCATE(N_cell2d(0:nx+1,0:ny+1))

        rho2d = 0.0
        flux2d = 0.0
        N_cell2d = 0

    ELSE IF(ndims .EQ. 3) THEN
		ALLOCATE(tr3d(1:N))
		ALLOCATE(cellHead3d(0:nx+1,0:ny+1,0:nz+1))
		ALLOCATE(rho3d(1:nx,1:ny,1:nz), flux3d(0:nx,0:ny,0:nz,3))
		ALLOCATE(N_cell3d(0:nx+1,0:ny+1,0:nz+1))
   
        rho3d = 0.0
        flux3d = 0.0
        N_cell3d = 0
    END IF

    dt1 = 0.1/5.0
    dx1 = 1./real(nx)
    dy1 = 1./real(ny) 
    dz1 = 1./real(ny) 

	IF(ndims .EQ. 1) THEN
		rho = 1.
		!flux(0:nx/2) = -1./4.
        !flux(1+nx/2:nx) = 1./4.
        flux(:)=0.251*dx1/dt1

    ELSE IF(ndims .EQ. 2) THEN
		rho2d = 1.0
		!flux2d(0:nx/2,:,1) = -0.25*dx1/dt1
		!flux2d(1+nx/2:nx,:,1) = 0.25*dx1/dt1
        !flux2d(:,0:ny/2,2) = -0.5*dx1/dt1
        !flux2d(:,1+ny/2:ny,2) = 0.5*dx1/dt1
        flux2d(:,:,1)=2.*dx1/dt1
        flux2d(:,:,2)=0.0

    ELSE IF(ndims .EQ. 3) THEN
		rho3d = 1.0
		flux3d(:,:,:,1) = 0.0
        flux3d(:,:,:,2) = 0.0
        flux3d(:,:,:,3) = 100.0
	END IF
  

	! set up solver object
    solver = MASSFLUXTRACER(ndims,nx,ny,nz,N,dx1,dy1,dz1)
    CALL solver%initialize_workpool()

    ! initialize tracers
	IF(ndims .EQ. 1) THEN
		CALL initialize_tracer_distribution_1d(tr1d, cellHead1d, N_cell1d)
	ELSE IF(ndims .EQ. 2) THEN
		CALL initialize_tracer_distribution_2d(tr2d, cellHead2d, N_cell2d)
    ELSE IF(ndims .EQ. 3) THEN
		CALL initialize_tracer_distribution_3d(tr3d, cellHead3d, N_cell3d)
	END IF

    
    ! Simulation loop    
    DO i= 1, nt
        PRINT*,'TIME STEP, % complete = ',i, (i*100./(nt*1.))
        
        CALL CPU_TIME(time_init)
        CALL solver%solve(dt1)
        CALL CPU_TIME(time_end)
        time_sim(i) = time_end-time_init

        !IF(ndims .EQ. 1)THEN
        !    CALL fileOut1d()
        !ELSE IF(ndims .EQ. 2)THEN
        !    CALL fileOut2d(i)	
        !ELSE IF(ndims .EQ. 3)THEN
        !    CALL fileOut3d_yz(i)	
        !END IF
 

    END DO

    PRINT*,'N=' ,N, 't(sec)=',time_sim 


    CALL solver%destroy_workpool()



	IF(ndims .EQ. 1) THEN
		DEALLOCATE(tr1d,cellHead1d,rho,flux,N_cell1d)
    ELSE IF(ndims .EQ. 2) THEN
		DEALLOCATE(tr2d,cellHead2d,rho2d,flux2d,N_cell2d)
    ELSE IF(ndims .EQ. 3) THEN
		DEALLOCATE(tr3d,cellHead3d,rho3d,flux3d,N_cell3d)
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
     		INTEGER :: i, l, cellnum, counter, ilow
            REAL :: p 
        

            N_cell(:) = 0
            N_cell(nx/2) = N
            counter = 1

            ilow = 1
            ! loop over cells
            DO  i = ilow+1,ilow+nx

                IF(N_cell(i) .GE. 1)  THEN
                    ! allocate memory for head tracer 
                    ALLOCATE(tr(counter)%p) 
                    cellHead1d(i)%p => tr(counter)%p
                    current => tr(counter)%p
                    ! assign id and cell number to tracer 
                    tr(counter)%p%id = counter
                    tr(counter)%p%x = i
                    counter = counter +1
                END IF

                ! loop over tracers
                DO l = 1, N_cell(i)-1
                    ALLOCATE(tr(counter)%p)
                    ! assign id and cell number to tracer 
                    tr(counter)%p%id = counter
                    tr(counter)%p%x = i

                    current%next => tr(counter)%p
                    current => tr(counter)%p
                    counter = counter+1                 
                END DO

            END DO
        

            GO TO 100

            ! loop over tracers, assign them to cells
            DO i = 1, N

                ! allocate memory for tracer pointer
                ALLOCATE(tr(i)%p) 

            	! draw a random integer between nx/3 and nx*2/3
                CALL RANDOM_NUMBER(p)
                cellnum= nx/2!nx/3+FLOOR(p*nx/3)

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

            100 CONTINUE


        	PRINT*, 'Tracer initialization completed.'

		END SUBROUTINE initialize_tracer_distribution_1d


		SUBROUTINE initialize_tracer_distribution_2d(tr, cellHead, N_cell)

			TYPE(tracer_ptr_2d), INTENT(INOUT) :: tr(:),cellHead(:,:)
            INTEGER, INTENT(INOUT) :: N_cell(:,:)

			TYPE(tracer_2d), POINTER :: current
     		INTEGER :: i, j, cellnumx, cellnumy
            REAL :: p 
        
            current => null()

            ! loop over tracers, assign them to cells
            DO i = 1, N

                ! allocate memory for tracer pointer
                ALLOCATE(tr(i)%p) 

                CALL RANDOM_NUMBER(p)
                cellnumx= 1!45+10*p!3*nx/4+FLOOR(p*nx/10)
                CALL RANDOM_NUMBER(p)
				cellnumy= 1!45+10*p!3*ny/4+FLOOR(p*ny/10)

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



		SUBROUTINE initialize_tracer_distribution_3d(tr, cellHead, N_cell)

			TYPE(tracer_ptr_3d), INTENT(INOUT) :: tr(:),cellHead(:,:,:)
            INTEGER, INTENT(INOUT) :: N_cell(:,:,:)

			TYPE(tracer_3d), POINTER :: current
     		INTEGER :: i, cellx, celly, cellz
            REAL :: p 
        
            current => null()

            ! loop over tracers, assign them to cells
            DO i = 1, N

                ! allocate memory for tracer pointer
                ALLOCATE(tr(i)%p) 

            	! draw a random integer between 1 and nx
                CALL RANDOM_NUMBER(p)
                cellx= 45+10*p!2*nx/5+FLOOR(p*nx/4)
                CALL RANDOM_NUMBER(p)
				celly= 45+10*p!2*ny/5+FLOOR(p*ny/4)
                CALL RANDOM_NUMBER(p)
				cellz= 45+10*p!2*nz/5+FLOOR(p*nz/4)

				! assign id and cell number to tracer 
                tr(i)%p%id = i
                tr(i)%p%x = cellx
				tr(i)%p%y = celly
                tr(i)%p%z = cellz

                N_cell(cellx+1,celly+1,cellz+1) = N_cell(cellx+1,celly+1,cellz+1)+1

                ! if cell is empty, designate this tracer as head
                IF(.NOT. ASSOCIATED(cellHead(cellx+1,celly+1,cellz+1)%p)) THEN
                   cellHead(cellx+1,celly+1,cellz+1)%p => tr(i)%p

				! otherwise add current tracer to the cell's linked list                
				ELSE
	              	current => cellHead(cellx+1,celly+1,cellz+1)%p

                     
                    DO WHILE(ASSOCIATED(current%next))
                        ! traverse through list until last tracer is reached
						current => current%next			
                    END DO
                    current%next => tr(i)%p
				END IF				     

			END DO


        	PRINT*, 'Tracer initialization completed.'

		END SUBROUTINE initialize_tracer_distribution_3d



        SUBROUTINE fileout1d()

            REAL::x
            INTEGER::i

                DO i=1,nx   
                    WRITE(1,*) i*dx1, N_cell1d(i) 
                END DO

        END SUBROUTINE fileout1d


        SUBROUTINE fileout2d(ts)
        
            INTEGER :: ts,i,j
            REAL :: x,y
            CHARACTER (len=40) :: filename
            CHARACTER (len=6) :: uniti


            IF(ts<10) THEN
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
  
            filename = trim('Output/t=')//trim(uniti)//trim('.txt')
    
            OPEN (UNIT=12,FILE=filename)

            DO j=1,ny
                y=j*dy1
                DO i=1,nx
                    x=i*dx1
                    WRITE(12,*) x,y,N_cell2d(i,j)
                END DO
            END DO

            CLOSE(12)

        END SUBROUTINE fileout2d

        SUBROUTINE fileout3d_xy(ts)
        
            INTEGER :: ts,i,j
            REAL :: x,y
            CHARACTER (len=40) :: filename
            CHARACTER (len=6) :: uniti


            IF(ts<10) THEN
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
  
            filename = trim('Output/t_xy=')//trim(uniti)//trim('.txt')
    
            OPEN (UNIT=13,FILE=filename)

            DO j=1,ny
                y=j*dy1
                DO i=1,nx
                    x=i*dx1
                    WRITE(13,*) x,y,N_cell3d(i,j,nz/2)
                END DO
            END DO

            CLOSE(12)

        END SUBROUTINE fileout3d_xy

        SUBROUTINE fileout3d_yz(ts)
        
            INTEGER :: ts,i,j
            REAL :: x,y
            CHARACTER (len=40) :: filename
            CHARACTER (len=6) :: uniti


            IF(ts<10) THEN
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
  
            filename = trim('Output/t_yz=')//trim(uniti)//trim('.txt')
    
            OPEN (UNIT=13,FILE=filename)

            DO j=1,ny
                y=j*dy1
                DO i=1,nx
                    x=i*dx1
                    WRITE(13,*) x,y,N_cell3d(nx/2,i,j)
                END DO
            END DO

            CLOSE(12)

        END SUBROUTINE fileout3d_yz

END PROGRAM tracer_test_run

