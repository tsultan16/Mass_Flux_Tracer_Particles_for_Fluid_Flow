MODULE tracersolver_mod

USE tracertype_mod
USE data_mod
USE MPI

IMPLICIT NONE

! tracer advection solver object

PUBLIC :: MASSFLUXTRACER

INTEGER, PRIVATE :: ndims, nx, ny, nz, N    
REAL, PRIVATE :: dt, dx, dy, dz

! work pool variables 
REAL, ALLOCATABLE, PRIVATE :: work_flux(:), work_rho(:)
INTEGER, ALLOCATABLE, PRIVATE :: work_N_cell1d(:), work_N_cell_temp(:)
TYPE (tree_root_1d), ALLOCATABLE, PRIVATE :: work_cellRoot1d(:)

TYPE :: MASSFLUXTRACER
    ! type data

    ! type bound procedures
    CONTAINS

    PROCEDURE, PASS(SELF), PUBLIC :: initialize_workpool
    PROCEDURE, PASS(SELF), PUBLIC :: destroy_workpool
    PROCEDURE, PASS(SELF), PUBLIC :: solve

    PROCEDURE, PASS(SELF), PRIVATE :: advect_1d
    PROCEDURE, PASS(SELF), PRIVATE :: cell_face_cycle_1d
    PROCEDURE, PASS(SELF), PRIVATE :: cell_transfer_1d


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

    SELECT CASE(ndims) 
            
    CASE(1) 

        ALLOCATE(work_rho(1:nx))
        ALLOCATE(work_flux(1:nx+1))
        ALLOCATE(work_N_cell1d(1:nx+2),work_N_cell_temp(1:nx+2))
        ALLOCATE(work_cellRoot1d(1:nx+2))
        work_rho(:) = 0.0
        work_flux(:) = 0.0
        work_N_cell1d(:) = 0
        work_N_cell_temp(:) = 0	  

    CASE(2)
    
    CASE(3)

    END SELECT

END SUBROUTINE initialize_workpool


SUBROUTINE destroy_workpool(SELF)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF

    SELECT CASE(ndims) 
            
    CASE(1) 

        DEALLOCATE(work_rho, work_flux, work_N_cell1d, &
                   work_N_cell_temp,work_cellRoot1d)

    CASE(2)
    
    CASE(3)

    END SELECT


END SUBROUTINE destroy_workpool


! Top-level routine for mass-flux tracer advection

SUBROUTINE solve(SELF, dt_in)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
    REAL, INTENT(IN) :: dt_in

    dt = dt_in

    SELECT CASE(ndims) 
            
    CASE(1) 

        CALL SELF%advect_1d(cellRoot1d,N_cell1d,flux,rho)

    CASE(2)
    
    CASE(3)

    END SELECT


END SUBROUTINE solve


! Routine for 1d tracer advection
!
! Loop over cells. Cycle through the cell faces (in the order x+, 
! x-). For each face: 
!
! 1) Compute outgoing fluid mass fractions (m_frac).
! 2) Compute number of outgoing tracers (Nout = m_frac * total num. of tracers in the cell).
! 3) Pick Nout tracers randomly and transfer them. 

SUBROUTINE advect_1d(SELF, cellRoot, N_cell, flux, rho)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF 
    TYPE(tree_root_1d), INTENT(INOUT) :: cellRoot(:)
    INTEGER, INTENT(INOUT) :: N_cell(:)
    REAL, INTENT(IN) :: flux(:), rho(:)

    TYPE(tracer_1d), POINTER :: current
    INTEGER :: i, ilow, l, ntr, move_count, face_count
    REAL :: p, mfrac, m0


    ilow = 1

    ! clear temp pointers
    current => null()

    ! copy values into work pool arrays
    work_N_cell1d(:) = N_cell(:)
    work_N_cell_temp(:) = N_cell(:)
    work_flux(:) = flux(:)
    work_rho(:) = rho(:) 
    work_cellRoot1d(:)=cellRoot(:)

    ! loop across cells in our 1d grid

    DO i = ilow+1, ilow+nx

        !PRINT*,' '
        !PRINT*,'CELL#, N_cell ',i-ilow,N_cell(i)
        !PRINT*,' '

        ! calculate fluid mass in cell
        m0=work_rho(i)*dx

        ! reset counter
        face_count = 0
  
        ! cycle through both cell faces
        DO WHILE(face_count .LT. 2 .AND. N_cell(i) .GT. 0)
             
            !PRINT*,'face_count=',face_count

            SELECT CASE(face_count)

            CASE(0)

                !PRINT*,'flux_x+ = ',work_flux(i)

                ! Outgoing through x+ face
                IF(SIGN(1.0,work_flux(i)) .GT. 0.0) THEN

                   ! PRINT*,'Outgoing through x+ face'

                    ! outgoing fluid mass fraction
                    mfrac = MIN(0.9999999,work_flux(i)*dt/m0)                   
                    m0 = m0*(1.0-mfrac)
                   !PRINT*,'mfrac=',mfrac               
                    CALL SELF%cell_face_cycle_1d(cellRoot,mfrac,i,i+1)

                END IF

            CASE(1)

                !PRINT*,'flux_x- = ',work_flux(i)

                ! Outgoing through x- face
                IF(SIGN(1.0,work_flux(i-1)) .LT. 0.0) THEN

                    !PRINT*,'Outgoing through x- face'

                    ! outgoing fluid mass fraction
                    mfrac = MIN(0.9999999,abs(work_flux(i-1))*dt/m0)
                    !PRINT*,'mfrac=',mfrac                
                    CALL SELF%cell_face_cycle_1d(cellRoot,mfrac,i,i-1)
                                       
                END IF

            END SELECT

            face_count = face_count +1

        END DO

    END DO

    ! copy work array data into storage arrays
    N_cell(:) = work_N_cell_temp(:)
   


END SUBROUTINE advect_1d



! routine for carrying out a single cell-face cycle 

SUBROUTINE cell_face_cycle_1d(SELF,cellRoot,mfrac, x_src, x_des)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
    TYPE(tree_root_1d), INTENT(INOUT) :: cellRoot(:)
    REAL, INTENT(IN) :: mfrac  
    INTEGER, INTENT(IN) :: x_src, x_des     ! x_src: source cell
                                            ! x_des: destination cell
    INTEGER :: i, Ntr, Nout, Nmove, N_br1, N_br2
    REAL :: p
    TYPE(tracer_1d), POINTER :: temp         
    
    REAL :: t1,t2

    ! compute number of tracers intially inside the cell
    N_br1 = work_cellRoot1d(x_src)%node_occupancy(1)
    N_br2 = work_cellRoot1d(x_src)%node_occupancy(2)
    Ntr = N_br1 + N_br2

    !PRINT*,'Number of tracers occupying source cell: ',Ntr
    !PRINT*,'Branch 1, Branch 2: : ',N_br1,N_br2

    

    ! compute number of outgoing tracers
    t1 = MPI_Wtime()
    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER(p)
    t2 = MPI_Wtime()
    randtime= randtime + t2 - t1

    IF(p .LT. (Ntr*mfrac-FLOOR(Ntr*mfrac)))THEN
        Nout = mfrac*Ntr+1
    ELSE
        Nout = mfrac*Ntr
    END IF

    work_N_cell_temp(x_src) = work_N_cell_temp(x_src) - Nout 
    work_N_cell_temp(x_des) = work_N_cell_temp(x_des) + Nout

    PRINT*,'Cell# ', x_src-1,Nout

    GO TO 100
        !*********************************************************************************
        PRINT*,' '
        PRINT*,'BEFORE Transfer: '
        PRINT*,' '
        PRINT*,'Tracers currently occupying source tree:'
        PRINT*,'Branch 1. No. of tracers=',cellRoot(x_src)%node_occupancy(1)
        IF(ASSOCIATED(cellRoot(x_src)%node_L%leaf_L)) THEN  
            temp => cellRoot(x_src)%node_L%leaf_L
            PRINT*,'Tracer# ',temp%id
            DO WHILE(ASSOCIATED(temp%next))
                temp => temp%next
                PRINT*,'Tracer# ',temp%id
            END DO
            NULLIFY(temp)
        END IF

        PRINT*,'Branch 2. No. of tracers=',cellRoot(x_src)%node_occupancy(2)
        IF(ASSOCIATED(cellRoot(x_src)%node_R%leaf_L)) THEN
            temp => cellRoot(x_src)%node_R%leaf_L
            PRINT*,'Tracer# ',temp%id
            DO WHILE(ASSOCIATED(temp%next))
                temp => temp%next
                PRINT*,'Tracer# ',temp%id
            END DO
            NULLIFY(temp)
        END IF
        

        PRINT*,'Tracers currently occupying destination tree:'
        PRINT*,'Branch 1. No. of tracers=',cellRoot(x_des)%node_occupancy(1)
        IF(ASSOCIATED(cellRoot(x_des)%node_L%leaf_L)) THEN
            temp => cellRoot(x_des)%node_L%leaf_L
            PRINT*,'Tracer# ',temp%id
            DO WHILE(ASSOCIATED(temp%next))
                temp => temp%next
                PRINT*,'Tracer# ',temp%id
            END DO
            NULLIFY(temp)
        END IF
        
        PRINT*,'Branch 2. No. of tracers=',cellRoot(x_des)%node_occupancy(2)
        IF(ASSOCIATED(cellRoot(x_des)%node_R%leaf_L)) THEN
            temp => cellRoot(x_des)%node_R%leaf_L
            PRINT*,'Tracer# ',temp%id
            DO WHILE(ASSOCIATED(temp%next))
                temp => temp%next
                PRINT*,'Tracer# ',temp%id
            END DO
            NULLIFY(temp)  
        END IF
        
        !*********************************************************************************
    100 CONTINUE
    
    ! Pick Nout tracers randomly from the cell and move them one by one
    DO i = 1, Nout

        ! generate a random intger in [1,Ntr]
        t1 = MPI_Wtime()
        CALL RANDOM_NUMBER(p)
        t2 = MPI_Wtime()
        randtime= randtime + t2 - t1
        Nmove = 1+p*Ntr

        !PRINT*,'Tracer ',Nmove, ' in the list is going to be moved'

        ! Find out which branch in the tree is occupied by this tracer, then move it
        IF(Nmove .LE. N_br1)THEN
            !PRINT*,'Occupied: Branch 1, Brach Position:', Nmove

            CALL SELF%cell_transfer_1d(cellRoot,1,Nmove,x_src, x_des)
            N_br1 = N_br1 - 1 
            work_cellRoot1d(x_src)%node_occupancy(1) = work_cellRoot1d(x_src)%node_occupancy(1)-1 ! update for other cell-face passes

        ELSE
            !PRINT*,'Occupied: Branch 2, Brach Position:', Nmove-N_br1

            CALL SELF%cell_transfer_1d(cellRoot,2,Nmove-N_br1,x_src, x_des)
            work_cellRoot1d(x_src)%node_occupancy(2) = work_cellRoot1d(x_src)%node_occupancy(2)-1 ! update for other cell-face passes

        END IF

        Ntr=Ntr-1

    END DO

    GO TO 101
      !*********************************************************************************
        PRINT*,' '
        PRINT*,'AFTER Transfer: '
        PRINT*,' '
        PRINT*,'Tracers currently occupying source tree:'
        PRINT*,'Branch 1. No. of tracers=',cellRoot(x_src)%node_occupancy(1)
        IF(ASSOCIATED(cellRoot(x_src)%node_L%leaf_L)) THEN
            temp => cellRoot(x_src)%node_L%leaf_L
            PRINT*,'Tracer# ',temp%id
            DO WHILE(ASSOCIATED(temp%next))
                temp => temp%next
                PRINT*,'Tracer# ',temp%id
            END DO
        END IF
        NULLIFY(temp)
        

        PRINT*,'Branch 2. No. of tracers=',cellRoot(x_src)%node_occupancy(2)
        IF(ASSOCIATED(cellRoot(x_src)%node_R%leaf_L)) THEN
            temp => cellRoot(x_src)%node_R%leaf_L
            PRINT*,'Tracer# ',temp%id
            DO WHILE(ASSOCIATED(temp%next))
                temp => temp%next
                PRINT*,'Tracer# ',temp%id
            END DO
        END IF
        NULLIFY(temp)       
       

        PRINT*,'Tracers currently occupying destination tree:'
        PRINT*,'Branch 1. No. of tracers=',cellRoot(x_des)%node_occupancy(1)
        IF(ASSOCIATED(cellRoot(x_des)%node_L%leaf_L)) THEN
            temp => cellRoot(x_des)%node_L%leaf_L
            PRINT*,'Tracer# ',temp%id
            DO WHILE(ASSOCIATED(temp%next))
                temp => temp%next
                PRINT*,'Tracer# ',temp%id
            END DO
            NULLIFY(temp)
        END IF
        
        PRINT*,'Branch 2. No. of tracers=',cellRoot(x_des)%node_occupancy(2)
        IF(ASSOCIATED(cellRoot(x_des)%node_R%leaf_L)) THEN  
            temp => cellRoot(x_des)%node_R%leaf_L
            PRINT*,'Tracer# ',temp%id
            DO WHILE(ASSOCIATED(temp%next))
                temp => temp%next
                PRINT*,'Tracer# ',temp%id
            END DO
        NULLIFY(temp) 
        END IF
               
        !*********************************************************************************
        101 CONTINUE
    
END SUBROUTINE cell_face_cycle_1d




! routine for transferring a single tracer between a source cell and a destination cell

SUBROUTINE cell_transfer_1d(SELF, cellRoot, branch, branch_pos, x_src, x_des)

    CLASS(MASSFLUXTRACER), INTENT(IN) :: SELF
    TYPE(tree_root_1d), INTENT(INOUT) :: cellRoot(:)
	INTEGER, INTENT(IN) :: branch, branch_pos, x_src, x_des   ! branch_pos is the outgoing tracer position within the branch   

    TYPE(tree_node_1d), POINTER :: current_node
    TYPE(tracer_1d), POINTER :: current_src, current_src_above, current_src_below
    INTEGER :: i
            
    ! clear temp pointers
    current_src => null()
    current_src_above => null()
    current_src_below => null()

    
    !*************************************************   
    ! Link outgoing tracer to tree in destination cell
    !*************************************************

    ! traverse to outgoing tracer location within src branch

    IF(branch .EQ. 1)THEN
        current_node => cellRoot(x_src)%node_L
        current_src => current_node%leaf_L
    ELSE
        current_node => cellRoot(x_src)%node_R 
        current_src => current_node%leaf_L
    END IF

    DO i= 1, branch_pos-1
        current_src_above => current_src
        current_src => current_src%next
    END DO
    IF(ASSOCIATED(current_src%next)) current_src_below => current_src%next
  

    !PRINT*,'Outgoing: Tracer# ',current_src%id

    ! ->go to location within destination tree where incoming tracer will be linked
    ! ->pick the branch with the fewest tracers (to avoid creating or worsening an imbalance)
    ! ->link incoming tracer to this branch

    !PRINT*,'Destination cell. Tracer number in Branch 1, 2 = ',cellRoot(x_des)%node_occupancy(1),&
    !                                                           cellRoot(x_des)%node_occupancy(2)

    IF(cellRoot(x_des)%node_occupancy(1) .LE. cellRoot(x_des)%node_occupancy(2))THEN
        !PRINT*,'Placing in Destination Branch 1.'
        ! check if destination branch is empty
        IF(cellRoot(x_des)%node_occupancy(1) .GT. 0) THEN
            ! if not empty, then link incoming tracer to tail of destination branch, and assign as new tail
            cellRoot(x_des)%node_L%leaf_R%next => current_src
            cellRoot(x_des)%node_L%leaf_R => current_src
        ELSE 
            !PRINT*,'Branch Empty'
            ! if empty, then incoming tracer is assigned both head and tail of destination branch
            cellRoot(x_des)%node_L%leaf_L => current_src
            cellRoot(x_des)%node_L%leaf_R => current_src
        END IF
        cellRoot(x_des)%node_occupancy(1) = cellRoot(x_des)%node_occupancy(1) + 1
    ELSE
        !PRINT*,'Placing in Destination Branch 2.'
        ! check if destination branch is empty
        IF(cellRoot(x_des)%node_occupancy(2) .GT. 0) THEN
            ! if not empty, then link incoming tracer to tail of destination branch, and assign as new tail
            cellRoot(x_des)%node_R%leaf_R%next => current_src
            cellRoot(x_des)%node_R%leaf_R => current_src
        ELSE 
            !PRINT*,'Branch Empty'
            ! if empty, then incoming tracer is assigned both head and tail of destination branch
            cellRoot(x_des)%node_R%leaf_L => current_src
            cellRoot(x_des)%node_R%leaf_R => current_src
        END IF
        cellRoot(x_des)%node_occupancy(2) = cellRoot(x_des)%node_occupancy(2) + 1
    END IF
    NULLIFY(current_src%next)

   
    ! update state variables of incoming tracer

    current_src%x = x_des

    cellRoot(x_src)%node_occupancy(branch) = cellRoot(x_src)%node_occupancy(branch) - 1 


    !*****************************************************   
    ! re-attach broken link left behind in the source tree
    !*****************************************************

    ! Case 1: Outgoing tracer was the only occupant of the branch
    !         -> nullify pointer to branch head/tail
    IF(.NOT. ASSOCIATED(current_src_above) .AND. .NOT. ASSOCIATED(current_src_below)) THEN

        !PRINT*,'Case 1. Outgoing tracer is sole ocupant of src branch'

        NULLIFY(current_node%leaf_L)
        NULLIFY(current_node%leaf_R)

    ! Case 2: Outgoing tracer was not the only occupant of the branch
    ELSE    

        !PRINT*,'Case 1. Src branch has >1 tracers.'

        ! Sub-case 1: Outgoing tracer was head of the branch
        !             -> Re-assign branch head
        IF(.NOT. ASSOCIATED(current_src_above))THEN
            !PRINT*,'Sub-case 1. Outgoing tracer was head of src branch'
            current_node%leaf_L => current_src_below
        END IF

        ! Sub-case 2: Outgoing tracer was neither head nor tail
        !             -> Re-attach break in link
        IF(ASSOCIATED(current_src_above) .AND. ASSOCIATED(current_src_below)) THEN
            !PRINT*,'Sub-case 2. Outgoing tracer was neither head nor tail of src branch'
            current_src_above%next => current_src_below
        END IF

        ! Sub-case 3: Outgoing tracer was branch tail
        !             -> Re-assign branch tail
        IF(.NOT. ASSOCIATED(current_src_below)) THEN
            !PRINT*,'Sub-case 3. Outgoing tracer was tail of src branch'
            NULLIFY(current_src_above%next)
            current_node%leaf_R => current_src_above
        END IF

    END IF

END SUBROUTINE cell_transfer_1d


END MODULE tracersolver_mod