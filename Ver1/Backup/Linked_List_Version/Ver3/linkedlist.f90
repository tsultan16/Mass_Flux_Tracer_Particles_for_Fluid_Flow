MODULE link_mod

	TYPE link
		REAL :: n
		TYPE(link) , pointer :: next => null()
	END TYPE link

END MODULE link_mod


PROGRAM linked_list_test
	USE link_mod
	IMPLICIT NONE

	! define pointers for linked list root and current elements
	TYPE(link),pointer :: root, current 


	INTEGER :: i = 0, m
	INTEGER :: io_stat_number = 0

	! define an allocatable array of real numbers
	REAL, ALLOCATABLE :: x(:)

	! define character array for file name
	CHARACTER (LEN = 80) :: fname = 'num.txt'

	OPEN(UNIT = 1, FILE = fname, status = 'old')

	! Allocate memory for linked list root element
	ALLOCATE(root)
	
    ! read in first number from file and put it in root
	READ(UNIT=1, FMT=*, iostat=io_stat_number) root%n

	! check end of file (iostat = 0 means not end of file)
    IF(io_stat_number .EQ. 0) THEN
		i = i + 1
        ALLOCATE(root%next)		
    END IF

	current => root
    ! read in the remaining numbers in the file (1 number per line)
    DO WHILE(ASSOCIATED(current%next))
    	current => current%next
		READ(UNIT=1, FMT=*, iostat=io_stat_number) current%n
		! check end of file 
    	IF(io_stat_number .EQ. 0) THEN
			i = i + 1
       		ALLOCATE(current%next)		
        END IF
    END DO

    m = i
    
    ! allocate memory for real number array
    ALLOCATE(x(1:m))

    i=1
   
    ! traverse through the linked list and copy the numbers into the array
    current => root
    
    DO WHILE(ASSOCIATED(current%next))
		x(i)=current%n
        current => current%next
        i = i + 1
    END DO
    
    PRINT*, m, ' numbers read from file.'
    PRINT*, m, ' elements of array x are:'
    
    DO i = 1 , m
    	PRINT*,x(i)
    END DO 

END PROGRAM linked_list_test
