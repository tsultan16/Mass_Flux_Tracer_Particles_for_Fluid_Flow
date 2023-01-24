MODULE flat_mod

	TYPE flat_base
		! type data
        
        ! type bound procedures		
		CONTAINS
            PROCEDURE, PASS(SELF), PUBLIC :: flatten_2d
            

	END TYPE flat_base

	CONTAINS
	
	SUBROUTINE flatten_2d(self, flux, rho, N_cell, work_flux, work_rho, work_N_cell, nx, ny)
    
    CLASS(flt_base) :: self

    REAL, INTENT(IN) :: flux(:,:,:), rho(:,:), N_cell(:,:)
    INTEGER INTENT(IN) :: nx, ny 
    REAL, INTENT(OUT) :: work_flux(:,:), work_rho(:), work_N_cell(:)

    INTEGER(WID) :: i,j, idx, dj
    
    ! Clean target arrays

    work_flux(:,:)  = 0.0
	work_rho(:)  = 0.0
    work_N_cell(:)  = 0
        
    dj = nx+2

    ! Always flatten running grids 
    DO j = 1, ny + 2 
        DO i = 1, nx + 2
            idx = (j-1)*dj + i
                 
            q(idx,1:8) = grid2d(i,j,:)
            b(idx,:)     = bface2d(i,j,:)
                
        END DO
    END DO
    
END SUBROUTINE flatten_2d




END MODULE flat_mod
