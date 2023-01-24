!  MASS-FLUX Tracer Particle Advection Solver (based
!  on Genel et al., 2013, MNRAS, 435, 1426)  
! 
!  Data structure used for handling the tracers is a 
!  self-balancing "AVL" Tree. 
!  (Using AVL Tree Algorithm from 'The Art of Computer 
!   Programming, Vol. 3, Knuth')
!


MODULE tracersolver_mod

USE tracertype_mod
USE data_mod
USE MPI

IMPLICIT NONE


! tracer advection solver object

PUBLIC :: MASSFLUXTRACER


INTEGER, PRIVATE :: ndims, nb, N    
REAL, PRIVATE :: dt, dx, dy, dz
INTEGER(KIND=4), PRIVATE :: s1, s2, s3

! work pool variables 
REAL, ALLOCATABLE, PRIVATE :: work_flux(:),          work_rho(:),      &
							  work_flux_2d(:,:,:),   work_rho_2d(:,:), &
						      work_flux_3d(:,:,:,:), work_rho_3d(:,:,:)

INTEGER, ALLOCATABLE, PRIVATE :: work_N_cell_1d(:),  work_N_cell_2d(:,:),  work_N_cell_3d(:,:,:)


TYPE :: MASSFLUXTRACER
    ! type data

    ! type bound procedures
    CONTAINS

    PROCEDURE, PASS(SELF), PUBLIC :: initialize_workpool
    PROCEDURE, PASS(SELF), PUBLIC :: destroy_workpool
    PROCEDURE, PASS(SELF), PUBLIC :: solve

    PROCEDURE, PASS(SELF), PRIVATE :: advect_1d
	PROCEDURE, PASS(SELF), PRIVATE :: advect_bndry_1d
    PROCEDURE, PASS(SELF), PRIVATE :: cell_face_cycle_1d

    PROCEDURE, PASS(SELF), PRIVATE :: advect_2d
	PROCEDURE, PASS(SELF), PRIVATE :: advect_bndry_2d
    PROCEDURE, PASS(SELF), PRIVATE :: cell_face_cycle_2d

    PROCEDURE, PASS(SELF), PRIVATE :: advect_3d
	PROCEDURE, PASS(SELF), PRIVATE :: advect_bndry_3d
    PROCEDURE, PASS(SELF), PRIVATE :: cell_face_cycle_3d

    PROCEDURE, PASS(SELF), PRIVATE :: insert_node
	PROCEDURE, PASS(SELF), PRIVATE :: insert_node_bndry
    PROCEDURE, PASS(SELF), PRIVATE :: check_balance_insert
    PROCEDURE, PASS(SELF), PRIVATE :: rebalance_insert

    PROCEDURE, PASS(SELF), PRIVATE :: delete_node
    PROCEDURE, PASS(SELF), PRIVATE :: check_balance_delete
    PROCEDURE, PASS(SELF), PRIVATE :: rebalance_delete

    PROCEDURE, PASS(SELF), PRIVATE :: create_tracer
    PROCEDURE, PASS(SELF), PRIVATE :: destroy_tracer

	PROCEDURE, PASS(SELF), PUBLIC :: pack_MPI_buffer_out
	PROCEDURE, PASS(SELF), PUBLIC :: unpack_MPI_buffer_in	

    PROCEDURE, PASS(SELF), PRIVATE :: rand_taus


END TYPE MASSFLUXTRACER

INTERFACE MASSFLUXTRACER
    MODULE PROCEDURE constructor
END INTERFACE 

    
CONTAINS



FUNCTION constructor(ndims_in, nb_in, N_in,dx_in,dy_in,dz_in) RESULT(SELF)

    TYPE(MASSFLUXTRACER) :: SELF
    INTEGER, INTENT(IN) :: ndims_in, nb_in, N_in
    REAL, INTENT(IN) :: dx_in, dy_in, dz_in

    ndims = ndims_in
    nb = nb_in
    N = N_in
    dx = dx_in
    dy = dy_in
    dz = dz_in
            
    ! set taus88 PRNG non-zero seed values (can take these from constuctor inputs, 
    ! but hardcoded for now.)
    s1 = 1234
    s2 = 5678
    s3 = 9123

    RETURN 
END FUNCTION constructor



SUBROUTINE initialize_workpool(SELF)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF

    SELECT CASE(ndims) 
            
    CASE(1) 
        ALLOCATE(work_rho(xlow:xhi))
        ALLOCATE(work_flux(xlow:xhi-1))
        ALLOCATE(work_N_cell_1d(xlow:xhi))
        work_rho(:) = 0.0
        work_flux(:) = 0.0
        work_N_cell_1d(:) = 0
        
    CASE(2)
        ALLOCATE(work_rho_2d(xlow:xhi,ylow:yhi))
        ALLOCATE(work_flux_2d(xlow:xhi-1,ylow:yhi-1,2))
        ALLOCATE(work_N_cell_2d(xlow:xhi,ylow:yhi))
        work_rho_2d(:,:) = 0.0
        work_flux_2d(:,:,:) = 0.0
        work_N_cell_2d(:,:) = 0
 
    CASE(3)
        ALLOCATE(work_rho_3d(xlow:xhi,ylow:yhi,zlow:zhi))
        ALLOCATE(work_flux_3d(xlow:xhi-1,ylow:yhi-1,zlow:zhi-1,3))
        ALLOCATE(work_N_cell_3d(xlow:xhi,ylow:yhi,zlow:zhi))
        work_rho_3d(:,:,:) = 0.0
        work_flux_3d(:,:,:,:) = 0.0
        work_N_cell_3d(:,:,:) = 0		
		
    END SELECT

END SUBROUTINE initialize_workpool



SUBROUTINE destroy_workpool(SELF)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF

    SELECT CASE(ndims) 
            
    CASE(1) 
        DEALLOCATE(work_rho, work_flux, work_N_cell_1d)
    CASE(2)
        DEALLOCATE(work_rho_2d, work_flux_2d, work_N_cell_2d)
    CASE(3)
        DEALLOCATE(work_rho_3d, work_flux_3d, work_N_cell_3d)	
    END SELECT


END SUBROUTINE destroy_workpool



! Top-level routine for mass-flux tracer advection

SUBROUTINE solve(SELF, dt_in)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
    REAL, INTENT(IN) :: dt_in

    dt = dt_in

    SELECT CASE(ndims) 
            
    CASE(1) 
        CALL SELF%advect_1d(cellHead_1d,N_cell_1d,flux_1d,rho_1d)
    CASE(2)
        CALL SELF%advect_2d(cellHead_2d,N_cell_2d,flux_2d,rho_2d)
    CASE(3)
        CALL SELF%advect_3d(cellHead_3d,N_cell_3d,flux_3d,rho_3d)
    END SELECT


END SUBROUTINE solve



! Routine for 1d tracer advection
!
! Loop over cells. 
! Cycle through the cell faces (in the order x+,x-). 
! For each face: 
!     1) Compute outgoing fluid mass fractions (m_frac).
!     2) Compute number of outgoing tracers, Nout (= m_frac * total num. of tracers in the cell).
!     3) Pick Nout tracers randomly and transfer them. 

SUBROUTINE advect_1d(SELF, cellHead, N_cell, flux, rho)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF 
    TYPE(node_ptr), INTENT(INOUT) :: cellHead(:)
    INTEGER, INTENT(INOUT) :: N_cell(:)
    REAL, INTENT(IN) :: flux(:), rho(:)

    INTEGER :: i, l, Ntr, Nout, face_count, j, ix
    REAL :: mfrac, m0
    REAL*8 :: p, t1, t2

    ! move particles from boundary cells (i.e. from neighboring domain) into interior cells
    CALL SELF%advect_bndry_1d(cellHead, cellHead_bndry_in, N_cell, N_cell_bndry_in)

    ! copy values into work pool arrays
    work_N_cell_1d(:) = N_cell(:)
    work_flux(:) = flux(:)
    work_rho(:) = rho(:) 
	
    ! loop across cells in our 1d grid and advect tracers
    DO i = xlow+nb, xhi-nb

        ! calculate fluid mass in cell
        m0 = work_rho(i)*dx

        ! reset counter
        face_count = 0

        Ntr = N_cell(i)

        ! cycle through both cell faces
        DO WHILE(face_count .LT. 2 .AND. Ntr .GT. 0)
             
            SELECT CASE(face_count)

            CASE(0)

                ! Outgoing through x+ face
                IF(SIGN(1.0,work_flux(i)) .GT. 0.0) THEN

                    ! outgoing fluid mass fraction
                    mfrac = MIN(0.9999999,work_flux(i)*dt/m0)                   
                    m0 = m0*(1.0-mfrac)
				ELSE
					mfrac = 0.d0
                END IF

				ix = 1

            CASE(1)

                ! Outgoing through x- face
                IF(SIGN(1.0,work_flux(i-1)) .LT. 0.0) THEN

                    ! outgoing fluid mass fraction
                    mfrac = MIN(0.9999999,ABS(work_flux(i-1))*dt/m0)
				ELSE		
					mfrac = 0.d0                   
                END IF

				ix = -1

            END SELECT

			! compute number of outgoing tracers
			t1 = MPI_Wtime()
			p = SELF%rand_taus()
			t2 = MPI_Wtime()
			randtime = randtime + t2 - t1

			! decimal part of outgoing tracer number is randomly set to either 0 or 1
			IF(p .LT. (Ntr*mfrac-FLOOR(Ntr*mfrac)))THEN
				Nout = mfrac*Ntr+1
			ELSE
				Nout = mfrac*Ntr
			END IF

            IF(Nout .GT. 0) THEN
				CALL SELF%cell_face_cycle_1d(Ntr,Nout,cellHead(i)%p,cellHead(i+ix)%p,mfrac,i,i+ix)
			END IF
			
            face_count = face_count +1

        END DO

    END DO
   
    ! copy work array data into storage arrays
    N_cell(:) = work_N_cell_1d(:)

	! update boundary cell count array
	N_cell_bndry_out(1) = N_cell(xhi)
	N_cell_bndry_out(2) = N_cell(xlow)
	


END SUBROUTINE advect_1d



! Subroutine for moving particles from boundary cells into interior cells.
! Boundary cells contain only particles that are exchanged with neighboring domains.

SUBROUTINE advect_bndry_1d(SELF, cellHead, cellHead_bndry_in, N_cell, N_cell_bndry_in)	

	CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF 
    TYPE(node_ptr), INTENT(INOUT) :: cellHead(xlow:xhi), cellHead_bndry_in(:)
    INTEGER, INTENT(INOUT) :: N_cell(xlow:xhi)
    INTEGER, INTENT(INOUT) :: N_cell_bndry_in(:)
	
    INTEGER :: bndry, n_bndry, i_p
	INTEGER :: x_des, y_des, new_rank
	
    TYPE(node), POINTER :: transfer_node, temp

	! loop acrosss all boundaries
	DO bndry = 1,2
	
	    ! find out how many tracers in boundary  
	    n_bndry = N_cell_bndry_in(bndry) 
		
	    DO i_p = 1, n_bndry
		
		    ! pull tracer out from the top of boundary list
		    transfer_node => cellHead_bndry_in(bndry)%p%node_R
			cellHead_bndry_in(bndry)%p%node_R => transfer_node%node_R
			transfer_node%node_R => null()
			
			! find out where to put the tracer
		    SELECT CASE(bndry)		   
				CASE(1)
					x_des = transfer_node%x + 2 
				CASE(2)	
					x_des = transfer_node%x - 2
			END SELECT
			
			! Update tracer state variables
			transfer_node%x = x_des
			
            !CALL print_tree(cellHead(x_des)%p)

			! insert into destination tree
			new_rank = N_cell(x_des) + 1
			CALL SELF%insert_node(cellHead(x_des)%p,transfer_node,new_rank)
         		
			! update cell counts	
			N_cell(x_des) = N_cell(x_des) + 1
		    N_cell_bndry_in(bndry) = N_cell_bndry_in(bndry) - 1  
		
            !CALL print_tree(cellHead(x_des)%p)
		
		END DO
	    	
	END DO


END SUBROUTINE advect_bndry_1d



! routine for carrying out a single cell-face cycle 

SUBROUTINE cell_face_cycle_1d(SELF, Ntr, Nout, head_src, head_des, mfrac, x_src, x_des)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
    TYPE(node), POINTER ,INTENT(INOUT) :: head_src, head_des
    REAL, INTENT(IN) :: mfrac  
    INTEGER, INTENT(INOUT) :: Ntr 
	INTEGER, INTENT(IN) :: Nout
    INTEGER, INTENT(IN) :: x_src, x_des    ! x_src: source cell
                                           ! x_des: destination cell
    INTEGER :: i, Nmove, new_rank, bndry
    TYPE(node), POINTER :: transfer_node        
    REAL*8 :: p, t1, t2

   
    ! Pick Nout tracers randomly from the cell and move them one by one
    DO i = 1, Nout

        ! generate a random integer between [1,Ntr]
        t1 = MPI_Wtime()
        p = SELF%rand_taus()
        t2 = MPI_Wtime()
        randtime = randtime + t2 - t1
		
        Nmove = 1+p*Ntr
		
        ! clear temp pointer
        transfer_node => null()
        
        t1 = MPI_Wtime()

        ! delete node from source tree
        CALL SELF%delete_node(head_src,transfer_node,Nmove)
      
        ! update tracer position (and other state variables...)
        transfer_node%x = x_des    

        ! insert node into destination tree
		IF(x_des .EQ. xlow .OR. x_des .EQ. xhi) THEN
            CALL SELF%insert_node_bndry(head_des,transfer_node)
		ELSE	
		    new_rank = work_N_cell_1d(x_des)+1
			CALL SELF%insert_node(head_des,transfer_node,new_rank)
        END IF


        t2 = MPI_Wtime()
        advecttime = advecttime + t2 - t1

        Ntr = Ntr - 1
        work_N_cell_1d(x_src) = work_N_cell_1d(x_src) - 1 
        work_N_cell_1d(x_des) = work_N_cell_1d(x_des) + 1

    END DO  

       
END SUBROUTINE cell_face_cycle_1d



! Routine for 2d tracer advection
!
! Loop over cells. 
! Cycle through the cell faces (in the order x+,x-,y+,y-). 
! For each face: 
!     1) Compute outgoing fluid mass fractions (m_frac).
!     2) Compute number of outgoing tracers, Nout (= m_frac * total num. of tracers in the cell).
!     3) Pick Nout tracers randomly and transfer them. 


SUBROUTINE advect_2d(SELF, cellHead, N_cell, flux, rho)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF 
    TYPE(node_ptr), INTENT(INOUT) :: cellHead(xlow:xhi,ylow:yhi)
    INTEGER, INTENT(INOUT) :: N_cell(xlow:xhi,ylow:yhi)
    REAL, INTENT(IN) :: flux(xlow:xhi-1,ylow:yhi-1,2), rho(xlow:xhi,ylow:yhi)

    INTEGER :: i, j, l, Ntr, Nout, face_count, ix, iy
    REAL :: mfrac, m0
    REAL*8 :: p, t1,t2


    ! move particles from boundary cells (i.e. from neighboring domain) into interior cells
    CALL SELF%advect_bndry_2d(cellHead, cellHead_bndry_in, N_cell, N_cell_bndry_in)

    ! copy values into work pool arrays
    work_N_cell_2d(:,:) = N_cell(:,:)
    work_flux_2d(:,:,:) = flux(:,:,:)
    work_rho_2d(:,:) = rho(:,:) 


    ! loop across cells in our 2d grid and advect tracers
    DO i = xlow+nb, xhi-nb
		DO j = ylow+nb, yhi-nb

	
			! calculate fluid mass in cell
			m0 = work_rho_2d(i,j)*dx*dy

			! reset counter
			face_count = 0

			Ntr = N_cell(i,j)

			! cycle through both cell faces
			DO WHILE(face_count .LT. 4 .AND. Ntr .GT. 0)
             
				SELECT CASE(face_count)

				CASE(0)

					! Outgoing through x+ face
					IF(SIGN(1.0,work_flux_2d(i,j,1)) .GT. 0.0) THEN
										
						! outgoing fluid mass fraction
						mfrac = MIN(0.9999999,work_flux_2d(i,j,1)*dt/m0)                   
						m0 = m0*(1.0-mfrac)
					ELSE	
						mfrac = 0.d0
					END IF

					ix = 1
					iy = 0
					
				CASE(1)

					! Outgoing through x- face
					IF(SIGN(1.0,work_flux_2d(i-1,j,1)) .LT. 0.0) THEN

						! outgoing fluid mass fraction
						mfrac = MIN(0.9999999,ABS(work_flux_2d(i-1,j,1))*dt/m0)
						m0 = m0*(1.0-mfrac)
					ELSE	
						mfrac = 0.d0                                     
					END IF
				
					ix = -1
					iy = 0  

				CASE(2)

					! Outgoing through y+ face
					IF(SIGN(1.0,work_flux_2d(i,j,2)) .GT. 0.0) THEN

						! outgoing fluid mass fraction
						mfrac = MIN(0.9999999,work_flux_2d(i,j,2)*dt/m0) 
						m0 = m0*(1.0-mfrac)
					ELSE	
						mfrac = 0.d0
					END IF
				
					ix = 0
					iy = 1
					
            CASE(3)

                ! Outgoing through y- face
                IF(SIGN(1.0,work_flux_2d(i,j-1,2)) .LT. 0.0) THEN

                    ! outgoing fluid mass fraction
                    mfrac = MIN(0.9999999,ABS(work_flux_2d(i,j-1,2))*dt/m0)

				ELSE
					mfrac = 0.d0
                END IF
                
				ix = 0
				iy = -1

            END SELECT

			! compute number of outgoing tracers
			t1 = MPI_Wtime()
			p = SELF%rand_taus()
			t2 = MPI_Wtime()
			randtime = randtime + t2 - t1

			! decimal part of outgoing tracer number is randomly set to either 0 or 1
			IF(p .LT. (Ntr*mfrac-REAL(FLOOR(Ntr*mfrac))))THEN
				Nout = mfrac*Ntr+1
			ELSE
				Nout = mfrac*Ntr
			END IF
	
            ! move tracers
            IF(Nout .GT. 0) THEN
			    PRINT*,'CELL#,NTR=',i,j,Ntr
			    PRINT*,'Moving ',Nout,' tracers from',i,j,' to ',i+ix,j+iy 			
				CALL SELF%cell_face_cycle_2d(Ntr,Nout,cellHead(i,j)%p,cellHead(i+ix,j+iy)%p,mfrac,i,i+ix,j,j+iy) 
			END IF
			
            face_count = face_count +1

        END DO

		END DO
    END DO

    ! copy work array data into storage arrays
    N_cell(:,:) = work_N_cell_2d(:,:)
   
   	! update boundary cell count array
	N_cell_bndry_out = 0
	
	DO j = ylow, yhi 
		N_cell_bndry_out(1) = N_cell_bndry_out(1) + N_cell(xhi ,j)
		N_cell_bndry_out(2) = N_cell_bndry_out(2) + N_cell(xlow,j)
	END DO
	
	DO i = xlow, xhi 
		N_cell_bndry_out(3) = N_cell_bndry_out(3) + N_cell(i,yhi)
		N_cell_bndry_out(4) = N_cell_bndry_out(4) + N_cell(i,ylow)
	END DO
	
	N_cell(xlow,:) = 0
	N_cell(xhi, :) = 0
	N_cell(:,ylow) = 0
	N_cell(:, yhi) = 0
    

	PRINT*,''
	DO j = yhi , ylow, -1
		DO i = xlow, xhi 
			WRITE(*,FMT='(i4)', ADVANCE = 'NO'), N_cell(i,j)
		END DO
		PRINT*,''
	END DO
	PRINT*,''
    	
	
END SUBROUTINE advect_2d



! Subroutine for moving particles from boundary cells into interior cells.
! Boundary cells contain only particles that are exchanged with neighboring domains.

SUBROUTINE advect_bndry_2d(SELF, cellHead, cellHead_bndry_in, N_cell, N_cell_bndry_in)	

	CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF 
    TYPE(node_ptr), INTENT(INOUT) :: cellHead(xlow:xhi,ylow:yhi), cellHead_bndry_in(:)
    INTEGER, INTENT(INOUT) :: N_cell(xlow:xhi,ylow:yhi)
    INTEGER, INTENT(INOUT) :: N_cell_bndry_in(:)
	
    INTEGER :: bndry, n_bndry, i_p
	INTEGER :: x_des, y_des, new_rank
	
    TYPE(node), POINTER :: transfer_node, temp

	! loop acrosss all boundaries
	DO bndry = 1,4	
	    
	    ! find out how many tracers in boundary  
	    n_bndry = N_cell_bndry_in(bndry) 
		
		PRINT*,'Boundary advection: Myrank, bndry, n_bndry = ',myrank,bndry,n_bndry
		

	    DO i_p = 1, n_bndry
		
		    ! pull tracer out from the top of boundary list
		    transfer_node => cellHead_bndry_in(bndry)%p%node_R
			cellHead_bndry_in(bndry)%p%node_R => transfer_node%node_R
			transfer_node%node_R => null()
			
			! find out where to put the tracer
			SELECT CASE(bndry)		   
				CASE(1)
					x_des = transfer_node%x + 2
					y_des = transfer_node%y 			
				CASE(2)	
					x_des = transfer_node%x - 2
					y_des = transfer_node%y 
				CASE(3)
					x_des = transfer_node%x 
					y_des = transfer_node%y + 2			
				CASE(4)	
					x_des = transfer_node%x 
					y_des = transfer_node%y - 2	
			END SELECT			

            !CALL print_tree(cellHead(x_des,y_des)%p)

            ! update tracer state variables
            transfer_node%x = x_des 
			transfer_node%y = y_des
			
			! insert into destination tree
			new_rank = N_cell(x_des,y_des) + 1
			CALL SELF%insert_node(cellHead(x_des,y_des)%p,transfer_node,new_rank)
         		
			! update cell counts	
			N_cell(x_des,y_des) = N_cell(x_des,y_des) + 1
		    N_cell_bndry_in(bndry) = N_cell_bndry_in(bndry) - 1  
		
            !CALL print_tree(cellHead(x_des,y_des)%p)
		
		END DO
	    	
	END DO


END SUBROUTINE advect_bndry_2d



! routine for moving tracers across a cell-face in 2d

SUBROUTINE cell_face_cycle_2d(SELF, Ntr, Nout ,head_src, head_des, mfrac, x_src, x_des, y_src, y_des)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
    TYPE(node), POINTER ,INTENT(INOUT) :: head_src, head_des
    REAL, INTENT(IN) :: mfrac  
    INTEGER, INTENT(INOUT) :: Ntr
    INTEGER, INTENT(IN) :: Nout	
    INTEGER, INTENT(IN) :: x_src, x_des, &    ! x_src, y_src: source cell
                           y_src, y_des       ! x_des, y_des: destination cell
    INTEGER :: i, Nmove, new_rank, ix, iy
    TYPE(node), POINTER :: transfer_node        
    REAL*8 :: p, t1,t2

	
    ! Pick Nout tracers randomly from the cell and move them one by one
    DO i = 1, Nout

        ! generate a random integer between [1,Ntr]
        t1 = MPI_Wtime()
        p = SELF%rand_taus()
        t2 = MPI_Wtime()
        randtime = randtime + t2 - t1


        Nmove = 1+p*Ntr

        ! clear temp pointer
        transfer_node => null()
        
        t1 = MPI_Wtime()
				
        ! delete node from source tree
        CALL SELF%delete_node(head_src,transfer_node,Nmove)

        ! update tracer position (and other state variables...)
        transfer_node%x = x_des 
        transfer_node%y = y_des

        ! insert node into destination tree
		IF(x_des .EQ. xlow .OR. x_des .EQ. xhi .OR. y_des .EQ. ylow .OR. y_des .EQ. yhi) THEN
            CALL SELF%insert_node_bndry(head_des,transfer_node)
		ELSE	
			new_rank = work_N_cell_2d(x_des,y_des)+1
			CALL SELF%insert_node(head_des,transfer_node,new_rank)
        END IF


        t2 = MPI_Wtime()
        advecttime = advecttime + t2 - t1

        Ntr = Ntr - 1
        work_N_cell_2d(x_src,y_src) = work_N_cell_2d(x_src,y_src) - 1 
        work_N_cell_2d(x_des,y_des) = work_N_cell_2d(x_des,y_des) + 1

    END DO  

       
END SUBROUTINE cell_face_cycle_2d



! Routine for 3d tracer advection
!
! Loop over cells. 
! Cycle through the cell faces (in the order x+,x-,y+,y-,z+,z-). 
! For each face: 
!     1) Compute outgoing fluid mass fractions (m_frac).
!     2) Compute number of outgoing tracers, Nout (= m_frac * total num. of tracers in the cell).
!     3) Pick Nout tracers randomly and transfer them. 


SUBROUTINE advect_3d(SELF, cellHead, N_cell, flux, rho)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF 
    TYPE(node_ptr), INTENT(INOUT) :: cellHead(xlow:xhi,ylow:yhi,zlow:zhi)
    INTEGER, INTENT(INOUT) :: N_cell(xlow:xhi,ylow:yhi,zlow:zhi)
    REAL, INTENT(IN) :: flux(xlow:xhi-1,ylow:yhi-1,zlow:zhi-1,3), rho(xlow:xhi,ylow:yhi,zlow:zhi)

    INTEGER :: i, j, k, l, Ntr, Nout, face_count, ix, iy, iz
    REAL :: mfrac, m0
    REAL*8 :: p, t1,t2


    ! move particles from boundary cells (i.e. from neighboring domain) into interior cells
    CALL SELF%advect_bndry_3d(cellHead, cellHead_bndry_in, N_cell, N_cell_bndry_in)

	! copy values into work pool arrays
    work_N_cell_3d(:,:,:) = N_cell(:,:,:)
    work_flux_3d(:,:,:,:) = flux(:,:,:,:)
    work_rho_3d(:,:,:) = rho(:,:,:) 

    ! loop across cells in our 2d grid and advect tracers
    DO i = xlow+nb, xhi-nb
		DO j = ylow+nb, yhi-nb
			DO k = zlow+nb, zhi-nb
			
		! calculate fluid mass in cell
				m0 = work_rho_3d(i,j,k)*dx*dy*dz

				! reset counter
				face_count = 0

				Ntr = N_cell(i,j,k)

				! cycle through both cell faces
				DO WHILE(face_count .LT. 6 .AND. Ntr .GT. 0)
             
					SELECT CASE(face_count)

					CASE(0)

						! Outgoing through x+ face
						IF(SIGN(1.0,work_flux_3d(i,j,k,1)) .GT. 0.0) THEN
										
							! outgoing fluid mass fraction
							mfrac = MIN(0.9999999,work_flux_3d(i,j,k,1)*dt/m0)                   
							m0 = m0*(1.0-mfrac)
						ELSE	
							mfrac = 0.d0
						END IF

						ix = 1
						iy = 0
						iz = 0
						
					CASE(1)

						! Outgoing through x- face
						IF(SIGN(1.0,work_flux_3d(i-1,j,k,1)) .LT. 0.0) THEN

							! outgoing fluid mass fraction
							mfrac = MIN(0.9999999,ABS(work_flux_3d(i-1,j,k,1))*dt/m0)
							m0 = m0*(1.0-mfrac)
						ELSE	
							mfrac = 0.d0                                     
						END IF
				
						ix = -1
						iy = 0
						iz = 0		

					CASE(2)

						! Outgoing through y+ face
						IF(SIGN(1.0,work_flux_3d(i,j,k,2)) .GT. 0.0) THEN

							! outgoing fluid mass fraction
							mfrac = MIN(0.9999999,work_flux_3d(i,j,k,2)*dt/m0) 
							m0 = m0*(1.0-mfrac)
						ELSE	
							mfrac = 0.d0
						END IF
				
						ix = 0
						iy = 1
						iz = 0
						
					CASE(3)

						! Outgoing through y- face
						IF(SIGN(1.0,work_flux_3d(i,j-1,k,2)) .LT. 0.0) THEN

							! outgoing fluid mass fraction
							mfrac = MIN(0.9999999,ABS(work_flux_3d(i,j-1,k,2))*dt/m0)
							m0 = m0*(1.0-mfrac)

						ELSE
							mfrac = 0.d0
						END IF
                
						ix = 0
						iy = -1
						iz = 0
	
					CASE(4)

						! Outgoing through z+ face
						IF(SIGN(1.0,work_flux_3d(i,j,k,3)) .GT. 0.0) THEN

							! outgoing fluid mass fraction
							mfrac = MIN(0.9999999,work_flux_3d(i,j,k,3)*dt/m0) 
							m0 = m0*(1.0-mfrac)
						ELSE	
							mfrac = 0.d0
						END IF
				
						ix = 0
						iy = 0
						iz = 1
						
					CASE(5)

						! Outgoing through z- face
						IF(SIGN(1.0,work_flux_3d(i,j-1,k,3)) .LT. 0.0) THEN

							! outgoing fluid mass fraction
							mfrac = MIN(0.9999999,ABS(work_flux_3d(i,j,k-1,3))*dt/m0)

						ELSE
							mfrac = 0.d0
						END IF
                
						ix = 0
						iy = 0
						iz = -1
	
					END SELECT

					! compute number of outgoing tracers
					t1 = MPI_Wtime()
					p = SELF%rand_taus()
					t2 = MPI_Wtime()
					randtime = randtime + t2 - t1

					! decimal part of outgoing tracer number is randomly set to either 0 or 1
					IF(p .LT. (Ntr*mfrac-REAL(FLOOR(Ntr*mfrac))))THEN
						Nout = mfrac*Ntr+1
					ELSE
						Nout = mfrac*Ntr
					END IF
	
					! move tracers
					IF(Nout .GT. 0) THEN
						CALL SELF%cell_face_cycle_3d(Ntr,Nout,cellHead(i,j,k)%p,cellHead(i+ix,j+iy,k+iz)%p,mfrac,i,i+ix,j,j+iy,k,k+iz) 
					END IF
			
					face_count = face_count +1

				END DO
				
			END DO
		END DO
    END DO

    ! copy work array data into storage arrays
    N_cell(:,:,:) = work_N_cell_3d(:,:,:)
   
   	! update boundary cell count array
	N_cell_bndry_out = 0
	
	DO k = zlow, zhi
		DO j = ylow, yhi 
			N_cell_bndry_out(1) = N_cell_bndry_out(1) + N_cell(xhi ,j,k)
			N_cell_bndry_out(2) = N_cell_bndry_out(2) + N_cell(xlow,j,k)
		END DO		
	END DO
	
	DO k = zlow, zhi	
		DO i = xlow, xhi 
			N_cell_bndry_out(3) = N_cell_bndry_out(3) + N_cell(i,yhi,k)
			N_cell_bndry_out(4) = N_cell_bndry_out(4) + N_cell(i,ylow,k)
		END DO	
	END DO
	
	DO j = ylow, yhi	
		DO i = xlow, xhi 
			N_cell_bndry_out(5) = N_cell_bndry_out(5) + N_cell(i,j,zhi)
			N_cell_bndry_out(6) = N_cell_bndry_out(6) + N_cell(i,j,zlow)
		END DO	
	END DO
	
	
	N_cell(xlow,:,:) = 0
	N_cell(xhi, :,:) = 0
	N_cell(:,ylow,:) = 0
	N_cell(:, yhi,:) = 0
	N_cell(:,:,zlow) = 0
	N_cell(:,:, zhi) = 0    

	PRINT*,''
	DO j = yhi , ylow, -1
		DO i = xlow, xhi 
			WRITE(*,FMT='(i4)', ADVANCE = 'NO'), N_cell(i,j,zlow)
		END DO
		PRINT*,''
	END DO
	PRINT*,''
    	
	
END SUBROUTINE advect_3d



! Subroutine for moving particles from boundary cells into interior cells.
! Boundary cells contain only particles that are exchanged with neighboring domains.

SUBROUTINE advect_bndry_3d(SELF, cellHead, cellHead_bndry_in, N_cell, N_cell_bndry_in)	

	CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF 
    TYPE(node_ptr), INTENT(INOUT) :: cellHead(xlow:xhi,ylow:yhi,zlow:zhi), cellHead_bndry_in(:)
    INTEGER, INTENT(INOUT) :: N_cell(xlow:xhi,ylow:yhi,zlow:zhi)
    INTEGER, INTENT(INOUT) :: N_cell_bndry_in(:)
	
    INTEGER :: bndry, n_bndry, i_p
	INTEGER :: x_des, y_des, z_des, new_rank
	
    TYPE(node), POINTER :: transfer_node, temp

	! loop acrosss all boundaries
	DO bndry = 1,6	
	    
	    ! find out how many tracers in boundary  
	    n_bndry = N_cell_bndry_in(bndry) 
		
		PRINT*,'Boundary advection: Myrank, bndry, n_bndry = ',myrank,bndry,n_bndry
		

	    DO i_p = 1, n_bndry
		
		    ! pull tracer out from the top of boundary list
		    transfer_node => cellHead_bndry_in(bndry)%p%node_R
			cellHead_bndry_in(bndry)%p%node_R => transfer_node%node_R
			transfer_node%node_R => null()
			
			! find out where to put the tracer
			SELECT CASE(bndry)		   
				CASE(1)
					x_des = transfer_node%x + 2
					y_des = transfer_node%y 			
				CASE(2)	
					x_des = transfer_node%x - 2
					y_des = transfer_node%y 
				CASE(3)
					x_des = transfer_node%x 
					y_des = transfer_node%y + 2			
				CASE(4)	
					x_des = transfer_node%x 
					y_des = transfer_node%y - 2	
				CASE(5)
					x_des = transfer_node%x 
					y_des = transfer_node%y 
					z_des = transfer_node%z + 2					
				CASE(6)	
					x_des = transfer_node%x 
					y_des = transfer_node%y
					z_des = transfer_node%z - 2					
			END SELECT			

            !CALL print_tree(cellHead(x_des,y_des,z_des)%p)

            ! update tracer state variables
            transfer_node%x = x_des 
			transfer_node%y = y_des
			transfer_node%y = z_des
			
			! insert into destination tree
			new_rank = N_cell(x_des,y_des,z_des) + 1
			CALL SELF%insert_node(cellHead(x_des,y_des,z_des)%p,transfer_node,new_rank)
         		
			! update cell counts	
			N_cell(x_des,y_des,z_des) = N_cell(x_des,y_des,z_des) + 1
		    N_cell_bndry_in(bndry) = N_cell_bndry_in(bndry) - 1  
		
            !CALL print_tree(cellHead(x_des,y_des,z_des)%p)
		
		END DO
	    	
	END DO


END SUBROUTINE advect_bndry_3d



! routine for moving tracers across a cell-face in 3d

SUBROUTINE cell_face_cycle_3d(SELF, Ntr, Nout ,head_src, head_des, mfrac, x_src, x_des, y_src, y_des, z_src, z_des)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
    TYPE(node), POINTER ,INTENT(INOUT) :: head_src, head_des
    REAL, INTENT(IN) :: mfrac  
    INTEGER, INTENT(INOUT) :: Ntr
    INTEGER, INTENT(IN) :: Nout	
    INTEGER, INTENT(IN) :: x_src, y_src, z_src, &    ! x_src, y_src, z_src: source cell
                           x_des, y_des, z_des       ! x_des, y_des, z_des: destination cell
						   
    INTEGER :: i, Nmove, new_rank, ix, iy, iz
    TYPE(node), POINTER :: transfer_node        
    REAL*8 :: p, t1,t2

	
    ! Pick Nout tracers randomly from the cell and move them one by one
    DO i = 1, Nout

        ! generate a random integer between [1,Ntr]
        t1 = MPI_Wtime()
        p = SELF%rand_taus()
        t2 = MPI_Wtime()
        randtime = randtime + t2 - t1


        Nmove = 1+p*Ntr

        ! clear temp pointer
        transfer_node => null()
        
        t1 = MPI_Wtime()
				
        ! delete node from source tree
        CALL SELF%delete_node(head_src,transfer_node,Nmove)

        ! update tracer position (and other state variables...)
        transfer_node%x = x_des 
        transfer_node%y = y_des
		transfer_node%z = z_des

        ! insert node into destination tree
		IF(x_des .EQ. xlow .OR. x_des .EQ. xhi .OR. y_des .EQ. ylow .OR. y_des .EQ. yhi &
		   .OR. z_des .EQ. zlow .OR. z_des .EQ. zhi) THEN
            CALL SELF%insert_node_bndry(head_des,transfer_node)
		ELSE	
			new_rank = work_N_cell_3d(x_des,y_des,z_des)+1
			CALL SELF%insert_node(head_des,transfer_node,new_rank)
        END IF


        t2 = MPI_Wtime()
        advecttime = advecttime + t2 - t1

        Ntr = Ntr - 1
		work_N_cell_3d(x_src,y_src,z_src) = work_N_cell_3d(x_src,y_src,z_src) - 1 
        work_N_cell_3d(x_des,y_des,z_des) = work_N_cell_3d(x_des,y_des,z_des) + 1

    END DO  

       
END SUBROUTINE cell_face_cycle_3d



!****************************************************************************************
!                                                                                       *
!                                   AVL TREE ROUTINES                                   *
!                                                                                       *
!****************************************************************************************



! Routine for a standard Binary Search Tree node insertion
! Note: Balance factors are adjuated along the way to the insertion location
! Only need to update balance factor for sub-trees that lie on the traversal path). 

SUBROUTINE insert_node(SELF,head,Q,k)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
    TYPE(node), POINTER, INTENT(INOUT) :: head, Q
    INTEGER, INTENT(IN) ::  k ! (node rank) k = N+1 (where N = # number of nodes in tree)
    

    TYPE(node), POINTER :: S,P,R,T
    INTEGER :: M, U    

    ! Set temp pointers
    T => head
    S => head
    P => head
    R => null()

    IF(ASSOCIATED(T%node_R)) THEN
        S => T%node_R
        P => T%node_R
        R => P%node_R
    END IF

    U = k   
    M = k

    IF(ASSOCIATED(P)) M = M - P%rank
	

    
    ! Traverse to insert location (insert is going to be the new right-most leaf)
    DO WHILE(ASSOCIATED(R))
    
        ! advance temp pointers
        IF(R%bf .NE. 0) THEN
            T => P
            S => R
            U = M
        END IF
      
        P => R
        M = M - P%rank
        R => P%node_R    

    END DO

    ! link newly inserted node
    P%node_R => Q
    ! set rank and balance factor
    Q%rank = 1
    Q%bf = 0    
    
    !update balance factors and check balance (only if the tree was initially non-empty)
	IF(k .GT. 1) THEN
		CALL SELF%check_balance_insert(T,S,P,Q,R,U)
	END IF

END SUBROUTINE insert_node



! Separate subroutine for insertions into boundary cells. Instead of using a tree in a boundary cell, we use 
! a linked list to keep things simple (since all particles in boundary cells will need to be packed up for MPI 
! communications, there is no performance benefit by using a tree over a linked list).

SUBROUTINE insert_node_bndry(SELF,head,Q)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
    TYPE(node), POINTER, INTENT(INOUT) :: head, Q    

    TYPE(node), POINTER :: temp


    IF(ASSOCIATED(head%node_R)) THEN
		temp => head%node_R

        ! link new node to the top of the list
		head%node_R => Q
		Q%node_R => temp
		
	ELSE
        head%node_R => Q
    END IF
		

END SUBROUTINE insert_node_bndry



! this subroutine adjusts the balance factors between nodes S and Q after an insert,
! then checks for imbalanced nodes on insertion path

! S points to imbalanced node
! R points to sub-tree of S in which the insertion occured
! T points to the parent of S
! P, Q both point to the newly inserted node 
SUBROUTINE check_balance_insert(SELF,T,S,P,Q,R,U)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
    TYPE(node), POINTER, INTENT(INOUT) :: T,S,P,Q,R
    INTEGER, INTENT(IN) :: U

    INTEGER :: M


    M = U
    R => S%node_R
    P => S%node_R 

    DO WHILE(.NOT. ASSOCIATED(P,Q)) ! loop exits when P=Q  
        P%bf = 1       
        M = M - P%rank
        P => P%node_R
    END DO

    ! S is imbalanced only if it already has a balance factor of 1 
    IF(S%bf .EQ. 0) THEN
        S%bf = 1   
    ELSE IF (S%bf .EQ. -1) THEN  
		S%bf = 0
    ELSE IF (S%bf .EQ. 1) THEN       
       CALL SELF%rebalance_insert(T,S,P,Q,R,U)
    END IF

END SUBROUTINE check_balance_insert



! this subroutine re-balances the tree after an insert

SUBROUTINE rebalance_insert(SELF,T,S,P,Q,R,U)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
    TYPE(node), POINTER, INTENT(INOUT) :: T,S,P,Q,R
    INTEGER, INTENT(IN) :: U


    ! Perform single rotation
    
    P => R        
    S%node_R => R%node_L
    R%node_L => S 
       
    ! update balance factors and ranks
    S%bf = 0
    R%bf = 0
    R%rank = R%rank + S%rank
       
    ! Assign new sub-tree root
    T%node_R => P

END SUBROUTINE rebalance_insert



! subroutine for deleting a node from a tree

SUBROUTINE delete_node(SELF,head,del_node,k)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
    TYPE(node), POINTER, INTENT(IN) :: head
    TYPE(node), POINTER, INTENT(INOUT) :: del_node
    INTEGER, INTENT(IN) :: k  ! rank of deletion node

    TYPE(node), POINTER:: T,S,P,Q,R  
    TYPE(path_ptr), POINTER:: p0,current, temp
    INTEGER :: M, a

    INTEGER:: temp_id

    ! initialize temp pointers
    T => null()
    S => null()
    P => head
    Q => head%node_R    
    R => null()
    a = 1

    M = k

    ALLOCATE(p0)
    p0%node => P
    p0%a = 1
    IF(ASSOCIATED(P%node_R)) THEN
        ALLOCATE(p0%next)
        current => p0%next
        current%prev => p0
    END IF

    P => Q

    ! search for delete node
    DO WHILE(M .NE. P%rank)

        ! Go right if rank is greater
        IF(M .GT. P%rank) THEN
            
            Q => P%node_R
            a = 1
            M = M - P%rank
              
        ! Go left if rank is smaller
        ELSE IF(M .LT. P%rank) THEN

            Q => P%node_L
            a = -1

            ! update rank of predecessor node after every left turn
            P%rank = P%rank-1

        END IF

        current%node => P     
        current%a = a
        ALLOCATE(current%next)
        current%next%prev => current
        current => current%next

        P => Q

    END DO
	
    !PRINT*,'Found delete node.'

    T => Q  ! Q and T both point to the node that will get deleted

    current%node => Q
   
    ! Case 1: Q is a leaf node
    IF(.NOT. ASSOCIATED(T%node_L) .AND. .NOT. ASSOCIATED(T%node_R)) THEN
               
        IF(a .EQ. -1)THEN
            current%prev%node%node_L => null()
        ELSE IF(a .EQ. 1)THEN
            current%prev%node%node_R => null()
        END IF

    ! Case 2: Q has empty right node
    !         Choose left node as sucecssor
    ELSE IF(.NOT. ASSOCIATED(T%node_R)) THEN
 
        R => T%node_L
      
        ! exchange contents of successor node with the original deletion node
        temp_id = Q%id
        Q%id = R%id
        Q%bf = R%bf
        Q%rank = R%rank 
        Q%node_R => R%node_R
        Q%node_L => R%node_L
        R%id = temp_id
        R%node_L => null()
        R%node_R => null()

        ! point T at the successor node
        T => R
        
        current%a = -a         
       
    ! Case 3: Q has empty left node
    ELSE IF(.NOT. ASSOCIATED(T%node_L)) THEN

        R => T%node_R
      
        ! exchange contents of successor node with the original deletion node
        temp_id = Q%id
        Q%id = R%id
        Q%bf = R%bf
        Q%rank = R%rank 
        Q%node_R => R%node_R
        Q%node_L => R%node_L
        R%id = temp_id
        R%node_L => null()
        R%node_R => null()

        ! point T at the successor node
        T => R

        current%a = -a
       
    ! Case 4: Q has non-empty left and right nodes. 
    ! Find a successor node.
    ELSE

        ! find a successor node from the right sub-tree of Q
        ! (i.e. a node on the left branch of right sub-tree of Q, that has only one or no child)
        R => T%node_R

        current%a = 1
        ALLOCATE(current%next)
        current%next%prev => current
        current => current%next
        current%node => R

        ! check for null left sub-tree of R 
        ! If null, then choose right node of R as the successor
        IF(.NOT. ASSOCIATED(R%node_L)) THEN

            ! copy contents of successor node into the original deletion node
            ! (leave rank unchanged)
            temp_id = Q%id
            Q%id = R%id
            Q%node_R => R%node_R
            R%id = temp_id
            R%node_R => null()

            ! point T at the successor node
            T => R

            current%a = -1

        ! Otherwise, find successor from left sub-tree of R
        ELSE
            ! traverse to nearest left node with a non-empty subtree
            S => R%node_L

            current%a = -1
            ALLOCATE(current%next)
            current%next%prev => current
            current => current%next
            current%node => S  

            ! update rank of R
            R%rank = R%rank -1            

            DO WHILE(ASSOCIATED(S%node_L))
                R => S
                S => S%node_L

                current%a = -1
                ALLOCATE(current%next)
                current%next%prev => current
                current => current%next
                current%node => S                  

                ! update rank of R
                R%rank = R%rank -1            

            END DO       
 
            R%node_L => S%node_R

            ! exchange contents of successor node with the original deletion node
            temp_id = Q%id
            Q%id = S%id             
            S%id = temp_id
            S%node_R => null()

            ! point T at the successor node
            T => S

            current%a = -1

        END IF

    END IF
    del_node => T

    ! check balance    
    CALL SELF%check_balance_delete(current,p0)

    ! Destroy temp pointers
    current => p0
    DO WHILE(ASSOCIATED(current%next))
        current => current%next
        DEALLOCATE(current%prev)
    END DO
    DEALLOCATE(current)    


END SUBROUTINE delete_node



! subroutine for checking balance after a node deletion

SUBROUTINE check_balance_delete(SELF,P_k, p0)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
    TYPE(path_ptr), POINTER, INTENT(INOUT):: P_k, p0

    TYPE(node), POINTER :: T,S,P,R
    LOGICAL :: CASE3_FLAG 

    ! clear temp pointers
    T => null()
    P => null()
    R => null()
    S => null()
 
    ! reset case 3 flag
    CASE3_FLAG = .FALSE.


    ! Adjust balance factors and ranks and rebalance if an unbalance node is found.
    ! Bottom up, starting from P_l-1 ... P_1
   
    P_k => P_k%prev
    DO WHILE(ASSOCIATED(P_k%prev) .AND. .NOT. CASE3_FLAG)

        IF(P_k%node%bf .EQ. P_k%a) THEN
            P_k%node%bf = 0

        ELSE IF(P_k%node%bf .EQ. 0 ) THEN
            P_k%node%bf = -P_k%a
            EXIT

        ELSE IF(P_k%node%bf .EQ. -P_k%a) THEN
         
            ! set temp pointers
            ! S points to imbalanced node
            ! R points to opposite sub-tree of S to where th deletion occured
            ! T points to the parent of S
            
            S => P_k%node
            IF(P_k%a .EQ. 1) THEN
                R => S%node_L
            ELSE IF(P_k%a .EQ. -1) THEN
                R => S%node_R
            END IF
            T => P_k%prev%node

            CALL SELF%rebalance_delete(T,S,R,P,P_k%a, CASE3_FLAG) 

            ! clear temp pointers
            T => null()
            P => null()
            R => null()
            S => null()

        END IF

        P_k => P_k%prev

    END DO


END SUBROUTINE check_balance_delete



! subroutine for rebalancing tree after a node deletion

SUBROUTINE rebalance_delete(SELF,T,S,R,P,a,CASE3_FLAG)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
    TYPE(node), POINTER, INTENT(INOUT) :: T,S,R,P
    INTEGER, INTENT(IN) :: a
    LOGICAL, INTENT(INOUT) :: CASE3_FLAG

    ! Case 1: Perform single rotation
    IF(R%bf .EQ. -a) THEN
        P => R
        IF(a .EQ. -1) THEN
            S%node_R => R%node_L
            R%node_L => S 
        ELSE IF(a .EQ. 1)THEN
            S%node_L => R%node_R
            R%node_R => S
        END IF

        ! update balance factors
        S%bf = 0
        R%bf = 0

        ! update ranks
        IF(a .EQ. -1) THEN
            R%rank = R%rank + S%rank
        ELSE
            S%rank = S%rank - R%rank
        END IF


    ! Case 2: Perform double rotation
    ELSE IF(R%bf .EQ. a) THEN

        IF(a .EQ. -1) THEN
            P => R%node_L
            R%node_L => P%node_R 
            P%node_R => R
            S%node_R => P%node_L
            P%node_L => S
        ELSE IF(a .EQ. 1)THEN
            P => R%node_R
            R%node_R => P%node_L 
            P%node_L => R
            S%node_L => P%node_R
            P%node_R => S
        END IF

        IF(P%bf .EQ. -a) THEN
            S%bf = a
            R%bf = 0
        ELSE IF(P%bf .EQ. 0) THEN
            S%bf = 0
            R%bf = 0
        ELSE IF(P%bf .EQ. a) THEN
            S%bf = 0
            R%bf = -a
        END IF

        P%bf = 0

        ! update ranks
        IF(a .EQ. -1) THEN
            R%rank = R%rank - P%rank
            P%rank = P%rank + S%rank
        ELSE
            P%rank = P%rank + R%rank
            S%rank = S%rank - P%rank
        END IF


    ! Case 3: Sub-tree balanced: Perform single rotation
    !        (same as Case 1, except balance factor of S remains unchanged)
    !         Rebalancing sequence needs to be terminated for Case 3.
    ELSE IF(R%bf .EQ. 0) THEN
       
        P => R

        IF(a .EQ. -1) THEN
            S%node_R => R%node_L
            R%node_L => S 
        ELSE IF(a .EQ. 1)THEN
            S%node_L => R%node_R
            R%node_R => S
        END IF

        ! update balance factors
        R%bf = a

        ! update ranks
        IF(a .EQ. -1) THEN
            R%rank = R%rank + S%rank
        ELSE
            S%rank = S%rank - R%rank
        END IF

        CASE3_FLAG = .TRUE.

    END IF

    ! Assign new sub-tree root
    IF(ASSOCIATED(S,T%node_R)) THEN
        T%node_R => P
    ELSE IF(ASSOCIATED(S,T%node_L)) THEN
        T%node_L => P
    END IF

END SUBROUTINE rebalance_delete



! Returns (double precision) random number between 0 and 1
! Adapted from Taus88 C-code in L’ECUYER(1996)

FUNCTION rand_taus(SELF) RESULT(r)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF

    REAL*8 :: r
    INTEGER(KIND = 4) :: b   ! Signed 32-bit integer


    b = ISHFT(IEOR(ISHFT(s1,13),s1),-19)    
    s1 = IEOR(ISHFT(IAND(s1,4294967294),12),b)
    b = ISHFT(IEOR(ISHFT(s2,2),s2),-25)
    s2 = IEOR(ISHFT(IAND(s2,4294967288),4),b)
    b = ISHFT(IEOR(ISHFT(s3,3),s3),-11)
    s3 = IEOR(ISHFT(IAND(s3,4294967280),17),b)

    r = 0.5D0+IEOR(s1,IEOR(s2,s3))*2.3283064365D-10  
    !CALL RANDOM_NUMBER(r)
    
    IF(r .GT. 0.9999999999D0 .OR. r .LT. 0.0000000001D0) THEN
		PRINT*,'Random number generator failed..r=',r
		STOP
	END IF

END FUNCTION rand_taus



!****************************************************************************************
!                                                                                       *
!                              MPI BUFFER PREP ROUTINES                                 *
!                                                                                       *
!****************************************************************************************



SUBROUTINE create_tracer(SELF, tracer)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
    TYPE(node), POINTER, INTENT(INOUT) :: tracer

    ALLOCATE(tracer)

END SUBROUTINE create_tracer



SUBROUTINE destroy_tracer(SELF, tracer)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
    TYPE(node), POINTER, INTENT(INOUT) :: tracer

    DEALLOCATE(tracer)

END SUBROUTINE destroy_tracer


! Subroutine for packing particles residing in boundary cells of our rectangular grid into an MPI buffer array. 
! For 1D, there are two boundaries: x+ and x-
! For 2D, there are four boundaries: x+, x-, y+, y-
! For 3D, there are six boundaries: x+, x-, y+, y-, z+, z-
!
! The boundary is specified by the 'bndry' index. The maximum number of particles that the buffer can hold is also specified.
! If there are more particles than will fit in the buffer, this subroutine needs to be called as many times as necessary 
! until all particles have been copied, as indicated by the 'all_particles_buffered' flag.

SUBROUTINE pack_MPI_buffer_out(SELF, buffer, N_cell_bndry_out, cellHead_bndry_out, bndry, max_particles, all_particles_buffered)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
	REAL, INTENT(INOUT) :: buffer(:)
    INTEGER, INTENT(IN) :: bndry, max_particles    ! bndry = 1 = x+, 2 = x-, 3 = y+, 4 = y-, 5 = z+, 6 = z-
    INTEGER, INTENT(INOUT) :: N_cell_bndry_out(:)
    TYPE (node_ptr), INTENT(INOUT) :: cellHead_bndry_out(:)
	LOGICAL, INTENT(INOUT) :: all_particles_buffered
	
    INTEGER :: i_p, counter, n_bndry, i_b, i_v, &
               i_bndry, n_remaining
    REAL(8) :: particle_data(4)
    LOGICAL :: buffer_full
    TYPE(node), POINTER :: top, temp


    ! find out number of particles currently in boundary cell
	n_bndry = N_cell_bndry_out(bndry)
 
    !####################################################################
    ! copy particle data into buffer (particles removed from top of list)
	!####################################################################

    buffer(1) =  MIN(max_particles, n_bndry)  ! first entry contains number of particles in buffer (not to exceed the buffer limit)


    i_b = 2
	
    DO i_p = 1, MIN(max_particles, n_bndry)

        ! pull out particle at the top   
        top => cellHead_bndry_out(bndry)%p%node_R

        particle_data(1) = top%id   
        particle_data(2) = top%x
        particle_data(3) = top%y
        particle_data(4) = top%z

		DO i_v = 1, 4
            buffer(i_b) = particle_data(i_v)
            i_b = i_b + 1    
        END DO

        ! link root with the next particle 
        cellHead_bndry_out(bndry)%p%node_R => top%node_R


		! destroy tracer from memory after its data has been copied
        CALL SELF%destroy_tracer(top)
		
	END DO	
	
    ! update particle number in boundary cells
    N_cell_bndry_out(bndry) = N_cell_bndry_out(bndry) - MIN(max_particles, n_bndry)

	! set particle buffer completion flag
    IF(n_bndry .LT. max_particles) THEN
        all_particles_buffered = .TRUE.
	END IF
	
	
	GO TO 101
    IF(myrank .EQ. 0) THEN
	
	PRINT*,''
	PRINT*,'myrank=',myrank
	PRINT*,'bndry=',bndry
	PRINT*,'MPI_OUT_PACKUP Status:'
	PRINT*,'Number of particles copied into buffer = ',MIN(max_particles, n_bndry)
	PRINT*,'Number of particles left to copy = ',N_cell_bndry_out(bndry)
	PRINT*,'all_particles_buffered = ',all_particles_buffered
	PRINT*,'STATUS(cellHead_bndry(bndry)%p%node_R)=',ASSOCIATED(cellHead_bndry_out(bndry)%p%node_R)
    PRINT*,''
    PRINT*,'Buffer contents = '
	PRINT*,buffer(1)
	DO i_p = 1, max_particles
		PRINT*,buffer(4*(i_p-1)+2),buffer(4*(i_p-1)+3),buffer(4*(i_p-1)+4),buffer(4*(i_p-1)+5)
	END DO
    PRINT*,''
	END IF
	101 CONTINUE
	
	
END SUBROUTINE pack_MPI_buffer_out



! Subroutine for unpacking particles from an MPI buffer array and placing them in boundary cells.

SUBROUTINE unpack_MPI_buffer_in(SELF, buffer, cellHead_bndry_in, N_cell_bndry_in, bndry)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
	REAL, INTENT(IN) :: buffer(:)
    INTEGER, INTENT(IN) :: bndry                  ! bndry = 1 = x+, 2 = x-, 3 = y+, 4 = y-, 5 = z+, 6 = z-
    TYPE (node_ptr), INTENT(INOUT) :: cellHead_bndry_in(:)
    INTEGER, INTENT(INOUT) :: N_cell_bndry_in(:)
	
    INTEGER :: i_p, counter, n_bndry, i_b, i_v, &
               i_bndry, n_remaining
    REAL(8) :: particle_data(4)
    TYPE(node), POINTER :: new_particle, temp


    !##########################################################
    ! Copy particle data out of MPI buffer into new particles. 
	! Place these new particles in cell boundary list
	!##########################################################

    ! find out how many particles in buffer
    n_bndry = buffer(1) 

    new_particle => null()
	
    DO i_p = 1, n_bndry

		! create a new particle
		CALL SELF%create_tracer(new_particle)

        ! copy particle data
        new_particle%id = buffer(4*(i_p-1)+2)   
        new_particle%x  = buffer(4*(i_p-1)+3)
        new_particle%y  = buffer(4*(i_p-1)+4)
        new_particle%z  = buffer(4*(i_p-1)+5)

        ! link the new particle to the cell boundary list (inserted at the top)
        new_particle%node_R => cellHead_bndry_in(bndry)%p%node_R
        cellHead_bndry_in(bndry)%p%node_R => new_particle		
		        		
	END DO	
	
    ! update particle number in boundary cells
    N_cell_bndry_in(bndry) = N_cell_bndry_in(bndry) + n_bndry

	GO TO 102
    IF(myrank .EQ. 0) THEN

	PRINT*,''
	PRINT*,'myrank=',myrank
	PRINT*,'bndry=',bndry
	PRINT*,'MPI_IN_UNPACK Status:'
	PRINT*,'Number of particles copied out of buffer = ',n_bndry
	PRINT*,'STATUS(cellHead_bndry_in(bndry)%p%node_R)=',ASSOCIATED(cellHead_bndry_in(bndry)%p%node_R)
    PRINT*,''
    PRINT*,'Buffer contents = '
	PRINT*,buffer(1)
	DO i_p = 1, n_bndry
		PRINT*,buffer(4*(i_p-1)+2),buffer(4*(i_p-1)+3),buffer(4*(i_p-1)+4),buffer(4*(i_p-1)+5)
	END DO
    PRINT*,''
	END IF
	102 CONTINUE

END SUBROUTINE unpack_MPI_buffer_in


!************************************************************************************8
! Subroutine for displaying tree contents on the terminal
! (Use for debugging purposes)

SUBROUTINE print_tree(head)

TYPE(node), POINTER, INTENT(IN) :: head

TYPE(node_ptr) :: temp0, temp1(2), temp2(4), temp3(8), temp4(16), temp5(32), temp6(64) 

INTEGER :: L0, L1(2), L2(4), L3(8), L4(16), L5(32), L6(64)
INTEGER :: B0, B1(2), B2(4), B3(8), B4(16), B5(32), B6(64) 
INTEGER :: R0, R1(2), R2(4), R3(8), R4(16), R5(32), R6(64) 
INTEGER :: i,j,k
LOGICAL :: node_stat_0, node_stat_1(2), node_stat_2(4), node_stat_3(8), node_stat_4(16), node_stat_5(32)


PRINT*,' '
PRINT*,' '

L0 = 0
L1 = 0
L2 = 0
L3 = 0 
L4 = 0
L5 = 0
L6 = 0 
B0 = 0
B1 = 0
B2 = 0
B3 = 0 
B4 = 0
B5 = 0
B6 = 0 
R0 = 0
R1 = 0
R2 = 0
R3 = 0 
R4 = 0
R5 = 0 
R6 = 0
node_stat_0 = .FALSE.
node_stat_1 = .FALSE.
node_stat_2 = .FALSE.
node_stat_3 = .FALSE.
node_stat_4 = .FALSE.
node_stat_5 = .FALSE.



temp0%p => head%node_R
node_stat_0 = ASSOCIATED(temp0%p)
IF(ASSOCIATED(temp0%p)) THEN
    L0 = temp0%p%id
    B0 = temp0%p%bf
    R0 = temp0%p%rank
    temp1(1)%p => temp0%p%node_L
    temp1(2)%p => temp0%p%node_R
ELSE
    temp1(1)%p => null()
    temp1(2)%p => null()
END IF


! Get Tree Level 1 values
DO i = 1,2
    IF(ASSOCIATED(temp1(i)%p))THEN 
        L1(i) = temp1(i)%p%id
        B1(i) = temp1(i)%p%bf
        R1(i) = temp1(i)%p%rank
		node_stat_1(i) = ASSOCIATED(temp1(i)%p)
    END IF
END DO


! Get Tree Level 2 values
DO i = 1,4
    j = MOD(i,2)+i/2
    temp2(i)%p => null()
    IF(MOD(i,2) .NE. 0)THEN
        IF(ASSOCIATED(temp1(j)%p))THEN
            temp2(i)%p => temp1(j)%p%node_L
        END IF
    ELSE
        IF(ASSOCIATED(temp1(j)%p))THEN
            temp2(i)%p => temp1(j)%p%node_R
        END IF
    END IF

    IF(ASSOCIATED(temp2(i)%p))THEN
        L2(i) = temp2(i)%p%id
        B2(i) = temp2(i)%p%bf
        R2(i) = temp2(i)%p%rank
		node_stat_2(i) = ASSOCIATED(temp2(i)%p)		
    END IF
    
END DO


! Get Tree Level 3 values
DO i = 1,8
    j = MOD(i,2)+i/2
    temp3(i)%p => null()
    IF(MOD(i,2) .NE. 0)THEN
        IF(ASSOCIATED(temp2(j)%p)) THEN
            temp3(i)%p => temp2(j)%p%node_L
        END IF
    ELSE
        IF(ASSOCIATED(temp2(j)%p)) THEN
            temp3(i)%p => temp2(j)%p%node_R
        END IF    
    END IF

    IF(ASSOCIATED(temp3(i)%p))THEN
        L3(i) = temp3(i)%p%id
        B3(i) = temp3(i)%p%bf  
        R3(i) = temp3(i)%p%rank                   
        node_stat_3(i) = ASSOCIATED(temp3(i)%p)		
    END IF
    
END DO


! Get Tree Level 4 values
DO i = 1, 16
    j = MOD(i,2)+i/2
    temp4(i)%p => null()
    IF(MOD(i,2) .NE. 0)THEN
        IF(ASSOCIATED(temp3(j)%p)) THEN
            temp4(i)%p => temp3(j)%p%node_L
        END IF
    ELSE
        IF(ASSOCIATED(temp3(j)%p)) THEN
            temp4(i)%p => temp3(j)%p%node_R
        END IF    
    END IF

    IF(ASSOCIATED(temp4(i)%p))THEN
        L4(i) = temp4(i)%p%id
        B4(i) = temp4(i)%p%bf                 
        R4(i) = temp4(i)%p%rank                 
		node_stat_4(i) = ASSOCIATED(temp4(i)%p)
    END IF
    
END DO

! Get Tree Level 5 values
DO i = 1, 32
    j = MOD(i,2)+i/2
    temp5(i)%p => null()
    IF(MOD(i,2) .NE. 0)THEN
        IF(ASSOCIATED(temp4(j)%p)) THEN
            temp5(i)%p => temp4(j)%p%node_L
        END IF
    ELSE
        IF(ASSOCIATED(temp4(j)%p)) THEN
            temp5(i)%p => temp4(j)%p%node_R
        END IF    
    END IF

    IF(ASSOCIATED(temp5(i)%p))THEN
        L5(i) = temp5(i)%p%id
        B5(i) = temp5(i)%p%bf
        R5(i) = temp5(i)%p%rank
		node_stat_5(i) = ASSOCIATED(temp5(i)%p)		
    END IF
    
END DO

! Get Tree Level 6 values
DO i = 1, 64
    j = MOD(i,2)+i/2
    temp6(i)%p => null()

    IF(MOD(i,2) .NE. 0)THEN
        IF(ASSOCIATED(temp5(j)%p)) THEN
            temp6(i)%p => temp5(j)%p%node_L
        END IF
    ELSE
        IF(ASSOCIATED(temp5(j)%p)) THEN
            temp6(i)%p => temp5(j)%p%node_R
        END IF    
    END IF

    IF(ASSOCIATED(temp6(i)%p))THEN
        L6(i) = temp6(i)%p%id
        B6(i) = temp6(i)%p%bf
        R6(i) = temp6(i)%p%rank       
    END IF
    
END DO



PRINT*,'IDs:'
WRITE(*,FMT ='(a60)',ADVANCE='no') ,''
WRITE(*,FMT =('(i5)')) L0
PRINT*,' '
PRINT*,' '

!WRITE(*,FMT ='(a1)',ADVANCE='no') ,''
DO i= 1,2
    WRITE(*,FMT =('(i44)'),ADVANCE='no') L1(i)
END DO
PRINT*,' '
PRINT*,' ' 

WRITE(*,FMT ='(a7)',ADVANCE='no') ,''
DO i=1,4
    WRITE(*,FMT =('(i24)'),ADVANCE='no') L2(i)
END DO
PRINT*,' '
PRINT*,' '

WRITE(*,FMT ='(a12)',ADVANCE='no') ,''
DO i = 1,8  
    WRITE(*,FMT =('(i12)'),ADVANCE='no') L3(i)
END DO
PRINT*,' '
PRINT*,' '

WRITE(*,FMT ='(a14)',ADVANCE='no') ,''
DO i = 1,16  
    WRITE(*,FMT =('(i6)'),ADVANCE='no') L4(i)
END DO
PRINT*,' '
PRINT*,' '


WRITE(*,FMT ='(a15)',ADVANCE='no') ,''
DO i = 1,32  
    WRITE(*,FMT =('(i3)'),ADVANCE='no') L5(i)
END DO
PRINT*,' '
PRINT*,' '
PRINT*,' '


!DO i = 1,64  
!    WRITE(*,FMT =('(i3)'),ADVANCE='no') L6(i)
!END DO
!PRINT*,' '
!PRINT*,' '
!PRINT*,' '


GO TO 100
!************************************************
PRINT*,'Balance Factors:'
WRITE(*,FMT ='(a95)',ADVANCE='no') ,''
WRITE(*,FMT =('(i3)')) B0
PRINT*,' '
PRINT*,' '

WRITE(*,FMT ='(a2)',ADVANCE='no') ,''
DO i= 1,2
    WRITE(*,FMT =('(i64)'),ADVANCE='no') B1(i)
END DO
PRINT*,' ' 
PRINT*,' '

WRITE(*,FMT ='(a17)',ADVANCE='no') ,''
DO i=1,4
    WRITE(*,FMT =('(i32)'),ADVANCE='no') B2(i)
END DO
PRINT*,' '
PRINT*,' '

WRITE(*,FMT ='(a25)',ADVANCE='no') ,''
DO i = 1,8  
    WRITE(*,FMT =('(i16)'),ADVANCE='no') B3(i)
END DO
PRINT*,' '
PRINT*,' '

WRITE(*,FMT ='(a29)',ADVANCE='no') ,''
DO i = 1,16  
    WRITE(*,FMT =('(i8)'),ADVANCE='no') B4(i)
END DO
PRINT*,' '
PRINT*,' '

WRITE(*,FMT ='(a31)',ADVANCE='no') ,''
DO i = 1,32  
    WRITE(*,FMT =('(i4)'),ADVANCE='no') B5(i)
END DO
PRINT*,' '
PRINT*,' '

WRITE(*,FMT ='(a32)',ADVANCE='no') ,''
DO i = 1,64  
    WRITE(*,FMT =('(i2)'),ADVANCE='no') B6(i)
END DO
PRINT*,' '
PRINT*,' '
PRINT*,' '

!************************************************
PRINT*,'Ranks:'
WRITE(*,FMT ='(a95)',ADVANCE='no') ,''
WRITE(*,FMT =('(i3)')) R0
PRINT*,' '
PRINT*,' '

WRITE(*,FMT ='(a2)',ADVANCE='no') ,''
DO i= 1,2
    WRITE(*,FMT =('(i64)'),ADVANCE='no') R1(i)
END DO
PRINT*,' ' 
PRINT*,' '

WRITE(*,FMT ='(a17)',ADVANCE='no') ,''
DO i=1,4
    WRITE(*,FMT =('(i32)'),ADVANCE='no') R2(i)
END DO
PRINT*,' '
PRINT*,' '

WRITE(*,FMT ='(a25)',ADVANCE='no') ,''
DO i = 1,8  
    WRITE(*,FMT =('(i16)'),ADVANCE='no') R3(i)
END DO
PRINT*,' '
PRINT*,' '

WRITE(*,FMT ='(a29)',ADVANCE='no') ,''
DO i = 1,16  
    WRITE(*,FMT =('(i8)'),ADVANCE='no') R4(i)
END DO
PRINT*,' '
PRINT*,' '

WRITE(*,FMT ='(a31)',ADVANCE='no') ,''
DO i = 1,32  
    WRITE(*,FMT =('(i4)'),ADVANCE='no') R5(i)
END DO
PRINT*,' '
PRINT*,' '

WRITE(*,FMT ='(a32)',ADVANCE='no') ,''
DO i = 1,64  
    WRITE(*,FMT =('(i2)'),ADVANCE='no') R6(i)
END DO
PRINT*,' '
PRINT*,' '
PRINT*,' '

!************************************************
PRINT*,'Node Pointer Status:'
WRITE(*,FMT ='(a60)',ADVANCE='no') ,''
WRITE(*,FMT =('(l5)')) node_stat_0
PRINT*,' '
PRINT*,' '

!WRITE(*,FMT ='(a1)',ADVANCE='no') ,''
DO i= 1,2
    WRITE(*,FMT =('(l44)'),ADVANCE='no') node_stat_1(i)
END DO
PRINT*,' '
PRINT*,' ' 

WRITE(*,FMT ='(a7)',ADVANCE='no') ,''
DO i=1,4
    WRITE(*,FMT =('(l24)'),ADVANCE='no') node_stat_2(i)
END DO
PRINT*,' '
PRINT*,' '

WRITE(*,FMT ='(a12)',ADVANCE='no') ,''
DO i = 1,8  
    WRITE(*,FMT =('(l12)'),ADVANCE='no') node_stat_3(i)
END DO
PRINT*,' '
PRINT*,' '

WRITE(*,FMT ='(a14)',ADVANCE='no') ,''
DO i = 1,16  
    WRITE(*,FMT =('(l6)'),ADVANCE='no') node_stat_4(i)
END DO
PRINT*,' '
PRINT*,' '


WRITE(*,FMT ='(a15)',ADVANCE='no') ,''
DO i = 1,32  
    WRITE(*,FMT =('(l3)'),ADVANCE='no') node_stat_5(i)
END DO
PRINT*,' '
PRINT*,' '
PRINT*,' '

100 CONTINUE


END SUBROUTINE print_tree



END MODULE tracersolver_mod