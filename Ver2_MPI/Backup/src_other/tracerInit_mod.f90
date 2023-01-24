MODULE tracerInit_mod

USE MPI
USE constants_mod
USE tracerType_mod
USE data_mod

IMPLICIT NONE

CONTAINS

! initial particle distribution function
FUNCTION fn(i,j) RESULT(fx)
    INTEGER, INTENT(IN) :: i,j
    INTEGER :: fx
    
    IF(i .GE. nx/2 - 1 .AND. i .LE. 2+nx/2 .AND. j .EQ. 1+ny/2) THEN
        fx = N/4   
    ELSE
        fx = 0
    END IF
   
   

END FUNCTION fn




SUBROUTINE initialize_tracer_distribution_1d(N_cell_1d, cellHead_1d)

	INTEGER, INTENT(INOUT) :: N_cell_1d(xlow:xhi)
	TYPE (node_ptr), INTENT(INOUT) :: cellHead_1d(xlow:xhi)
    INTEGER:: i, ilow, tr_counter, my_start_id, my_end_id



 
    my_particle_num = 0
		
    ! intiliaze tracer distribution	
	DO i = xlow, xhi
	    N_cell_1d(i) = fn(i,1)	
	    my_particle_num = my_particle_num + N_cell_1d(i)
	END DO

	! compute start ID for my particles
    IF(myrank .EQ. 0) THEN

       my_start_id = 1
       my_end_id = my_start_id + my_particle_num - 1
       CALL MPI_SEND(my_end_id, 1, MPI_INTEGER, myrank+1, myrank, MPI_COMM_WORLD, ierr)

    ELSE IF(myrank .EQ. numprocs(1)-1) THEN

       CALL MPI_RECV(my_start_id, 1, MPI_INTEGER, myrank-1, myrank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
       my_start_id = my_start_id + 1
       my_end_id = my_start_id + my_particle_num - 1
       
    ELSE

       CALL MPI_RECV(my_start_id, 1, MPI_INTEGER, myrank-1, myrank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
       my_start_id = my_start_id + 1
       my_end_id = my_start_id + my_particle_num - 1
       CALL MPI_SEND(my_end_id, 1, MPI_INTEGER, myrank+1, myrank, MPI_COMM_WORLD, ierr) 
   
    END IF
    
    ! initialize counter
    tr_counter = my_start_id
	
    DO i = xlow, xhi

        ALLOCATE(cellHead_1d(i)%p)
		cellHead_1d(i)%p%id = -1
		cellHead_1d(i)%p%rank = 1
		cellHead_1d(i)%p%bf = 0

		
        PRINT*,'Cell #, N = ',i, N_cell_1d(i)
        IF(N_cell_1d(i) .GT. 0) THEN            
            CALL tree_init(cellHead_1d(i)%p,N_cell_1d(i),tr_counter,i-ilow,-1,-1)    
        END IF
	END DO

	PRINT*,''
    DO i = xlow, xhi 
	   WRITE(*,FMT='(i4)', ADVANCE = 'NO'), N_cell_1d(i)
    END DO
	
	PRINT*,''
  

    PRINT*, 'Tracer initialization completed.'

END SUBROUTINE initialize_tracer_distribution_1d



SUBROUTINE initialize_tracer_distribution_2d(N_cell_2d, cellHead_2d)

 	INTEGER, INTENT(INOUT) :: N_cell_2d(xlow:xhi,ylow:yhi)
	TYPE (node_ptr), INTENT(INOUT) :: cellHead_2d(xlow:xhi,ylow:yhi)

    INTEGER:: i, ilow, j, jlow, tr_counter, my_start_id, my_end_id


    my_particle_num = 0

    !Initialize tracer distribution
	DO i = xlow, xhi
		DO j = ylow, yhi
			N_cell_2d(i,j) = fn(i,j)
			my_particle_num = my_particle_num + N_cell_2d(i,j)
		END DO
    END DO

	! compute start ID for my particles
    IF(myrank .EQ. 0) THEN

       my_start_id = 1
       my_end_id = my_start_id + my_particle_num - 1
       CALL MPI_SEND(my_end_id, 1, MPI_INTEGER, myrank+1, myrank, MPI_COMM_WORLD, ierr)

    ELSE IF(myrank .EQ. numprocs(1)-1) THEN

       CALL MPI_RECV(my_start_id, 1, MPI_INTEGER, myrank-1, myrank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
       my_start_id = my_start_id + 1
       my_end_id = my_start_id + my_particle_num - 1
       
    ELSE

       CALL MPI_RECV(my_start_id, 1, MPI_INTEGER, myrank-1, myrank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
       my_start_id = my_start_id + 1
       my_end_id = my_start_id + my_particle_num - 1
       CALL MPI_SEND(my_end_id, 1, MPI_INTEGER, myrank+1, myrank, MPI_COMM_WORLD, ierr) 
   
    END IF

    PRINT*,''
	PRINT*,'Tracer: start_id, end_id = ',my_start_id,my_end_id
	PRINT*,''
	

   
    ! initialize counter
    tr_counter = my_start_id
	

    ! loop over all cells
    DO i = xlow, xhi
		DO j = ylow, yhi

			ALLOCATE(cellHead_2d(i,j)%p)
			cellHead_2d(i,j)%p%id = -1
			cellHead_2d(i,j)%p%rank = 1
			cellHead_2d(i,j)%p%bf = 0
			cellHead_2d(i,j)%p%x = i
			cellHead_2d(i,j)%p%y = j
		
			!PRINT*,'Cell #, N = ',i,j, N_cell_2d(i,j)
			IF(N_cell_2d(i,j) .GT. 0) THEN            
				CALL tree_init(cellHead_2d(i,j)%p,N_cell_2d(i,j),tr_counter,i,j,-1)    
			END IF
		END DO
    END DO

    PRINT*, 'Tracer initialization completed.'


	DO j = yhi , ylow, -1
		DO i = xlow, xhi 
			WRITE(*,FMT='(i4)', ADVANCE = 'NO'), N_cell_2d(i,j)
		END DO
		PRINT*,''
	END DO
	PRINT*,''
    
END SUBROUTINE initialize_tracer_distribution_2d



! this subroutine inserts N (unique) nodes into an empty tree 

SUBROUTINE tree_init(head, N, tr_counter, tr_x, tr_y, tr_z)
    
    TYPE(node), POINTER, INTENT(INOUT) :: head
    INTEGER, INTENT(IN) :: N, tr_x, tr_y, tr_z 
    INTEGER, INTENT(INOUT) :: tr_counter

    TYPE(node), POINTER :: T ! temp pointers
    INTEGER :: counter   

    ! clear temp pointers
    T => null()
    
    ! reset counters
    counter =  0

    ! Allocate memory for root node
    CALL create_tracer(head%node_R) 

    ! sey root balance factor
    head%node_R%bf = 0
    ! set root node rank
    head%node_R%rank = 1
    
    ! set tracer id and position
    head%node_R%id = tr_counter
    head%node_R%x = tr_x
    head%node_R%y = tr_y
    head%node_R%z = tr_z

    tr(tr_counter)%p => head%node_R

    tr_counter= tr_counter + 1        

    ! insert N-1 nodes into the tree
    DO WHILE(counter .LT. N-1)

        T => head
 
        !CALL print_tree(head)

        !PRINT*,'Inserting new node with Rank =',counter+2


        CALL insert_node(T,counter+2,tr_counter, tr_x, tr_y, tr_z)

        
        !PRINT*,'Finished inserting new node.'

               
        counter = counter + 1
        tr_counter = tr_counter + 1
    END DO

    PRINT*,'Tree has been populated with tracers.'
    !CALL print_tree(head)


END SUBROUTINE tree_init


SUBROUTINE create_tracer(tracer)

    TYPE(node), POINTER, INTENT(INOUT) :: tracer

    ALLOCATE(tracer)

END SUBROUTINE create_tracer


SUBROUTINE destroy_tracer(tracer)

    TYPE(node), POINTER, INTENT(INOUT) :: tracer

    DEALLOCATE(tracer)

END SUBROUTINE destroy_tracer





! Routine for a standard Binary Search Tree node insertion
! Note: Balance factors are adjuated along the way to the insertion location
! Only need to update balance factor for sub-trees that lie on the traversal path). 

SUBROUTINE insert_node(head,k,tr_id, tr_x, tr_y, tr_z)

    TYPE(node), POINTER, INTENT(INOUT) :: head
    INTEGER, INTENT(IN) ::  k, tr_id, tr_x, tr_y, tr_z ! (node rank) k = N+1 (where N = # number of nodes in tree)
    

    TYPE(node), POINTER :: S,P,Q,R,T
    INTEGER :: M, U


    T => head

    ! Reset temp pointers
    R => null()
    S => T%node_R
    P => T%node_R
    Q => null()

    U = k   
    M = k


    M = M - P%rank
    R => P%node_R
    
    ! Traverse to insert location (insert is going to the new right-most leaf)

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

    ! allocate memory for the new node
    ALLOCATE(Q)
    ! link newly created node
    P%node_R => Q
    ! set rank and balance factor
    Q%rank = 1
    Q%bf = 0
    ! set tracer id and position
    Q%id = tr_id
    Q%x = tr_x  
    Q%y = tr_y  
    Q%z = tr_z    


    M = U

    tr(tr_id)%p => Q  
  
    !PRINT*,'Node inserted.'
  
  
    !update balance factors
    CALL adjust_bf(T,S,P,Q,R,U)


END SUBROUTINE insert_node



! this subroutine adjusts the balance factors between nodes S and Q after an insert

SUBROUTINE adjust_bf(T,S,P,Q,R,U)

    TYPE(node), POINTER, INTENT(INOUT) :: T,S,P,Q,R
    INTEGER, INTENT(IN) :: U

    INTEGER :: M

    !PRINT*,'Updating balance factors...'

    M = U
    R => S%node_R
    P => S%node_R 

    !PRINT*,'Rank(S),Rank(P)=',S%rank,P%rank

    DO WHILE(.NOT. ASSOCIATED(P,Q)) ! loop exits when P=Q  
        P%bf = 1       
        M = M - P%rank
        P => P%node_R
    END DO

    !PRINT*,'M,Rank(P)=',M,P%rank

    CALL check_balance_insert(T,S,P,Q,R,U)

END SUBROUTINE adjust_bf



! this subroutine checks for imbalanced nodes

! S points to imbalanced node
! R points to sub-tree of S in which the insertion occured
! T points to the parent of S
! P, Q both point to the newly inserted node 
SUBROUTINE check_balance_insert(T,S,P,Q,R,U)

    TYPE(node), POINTER, INTENT(INOUT) :: T,S,P,Q,R
    INTEGER, INTENT(IN) :: U

    !PRINT*,'Checking for imbalance...'

    ! no imbalances 
    IF(S%bf .EQ. 0) THEN
        S%bf = 1   
    ELSE IF (S%bf .EQ. 1) THEN       
       ! PRINT*,'Found imbalance! Rebalancing sub-tree...'
       CALL rebalance_insert(T,S,P,Q,R,U)
    END IF

END SUBROUTINE check_balance_insert


! this subroutine re-balances the tree after an insert

SUBROUTINE rebalance_insert(T,S,P,Q,R,U)

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


! Subroutine for displaying tree contents on the terminal
! (Use for debugging purposes)

SUBROUTINE print_tree(head)

TYPE(node), POINTER, INTENT(IN) :: head

TYPE(node_ptr) :: temp0, temp1(2), temp2(4), temp3(8), temp4(16), temp5(32), temp6(64) 

INTEGER :: L0, L1(2), L2(4), L3(8), L4(16), L5(32), L6(64)
INTEGER :: B0, B1(2), B2(4), B3(8), B4(16), B5(32), B6(64) 
INTEGER :: R0, R1(2), R2(4), R3(8), R4(16), R5(32), R6(64) 
INTEGER :: i,j,k

LOGICAL, PARAMETER :: bf_option = .FALSE., rank_option = .FALSE.

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

temp0%p => head%node_R
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


DO i = 1,64  
    WRITE(*,FMT =('(i3)'),ADVANCE='no') L6(i)
END DO
PRINT*,' '
PRINT*,' '
PRINT*,' '


!************************************************
IF(rank_option) THEN

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


END IF

!************************************************
IF(bf_option) THEN

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

END IF

END SUBROUTINE print_tree





END MODULE tracerInit_mod
