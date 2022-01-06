SUBROUTINE Analyze(A,c, x_opt, STAT, b, x0)

REAL, INTENT(IN) :: A(:,:), c(:), x_opt(:)
LOGICAL, INTENT(INOUT) :: STAT
REAL, OPTIONAL :: b(:), x0(:)

REAL :: x_sci(size(x_opt)), A_dual(size(A,2),size(A,1)), b_dual(size(c)), c_dual(size(A,1)), x_dual(size(A,1)), eps, error


IF((.NOT. PRESENT(b)) .OR. (PRESENT(b) .AND. PRESENT(x0))) THEN

	error = DOT(c,x_opt)
	eps = 1.E-3

	IF( ABS(error) < eps ) THEN
		STAT = .TRUE.
	ELSE
		STAT = .FALSE.
	END IF

ELSE

	CALL Dual(A, b, c, A_dual, b_dual, c_dual)

	x_dual = LPSolver(A_dual,c_dual, b_dual)
	
	x_sci = ScilabLPSolver(A, b, c)	

	error = (ABS(DOT(c,x_opt)+DOT(c_dual,x_dual))*100)/ABS(DOT(c,x_opt)-DOT(c_dual,x_dual))

	eps = 1.E-01

	IF( error < eps ) THEN
		STAT = .TRUE.
	ELSE
		STAT = .FALSE.	
	END IF

END IF

END SUBROUTINE Analyze
