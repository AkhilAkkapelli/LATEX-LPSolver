FUNCTION LPSolver(A,c, b, x0) RESULT(x_opt)


REAL, INTENT(IN) :: A(:,:), c(:)
REAL, OPTIONAL :: b(:), x0(:)

REAL :: x_opt(size(A,2))

REAL :: A_can(size(A,1)+size(A,2)+1,2*(size(A,1)+size(A,2)+1)), c_can(2*(size(A,1)+size(A,2)+1)), &
		 x_can(2*(size(A,1)+size(A,2)+1)), x(2*(size(A,1)+size(A,2))+1), &
			 A_art(size(A,1),size(A,2)+1), c_art(size(c)+1), x_art(size(A,2)+1), eps


IF(.NOT. PRESENT(b)) THEN

	x_opt = Karmarkar(A, c)

ELSE IF( PRESENT(x0) ) THEN

	eps = 1.E-04

	!!!Feasibility Check
	IF(ANY(ABS(GEMV(A,x0)-b) >= eps))  STOP "Initial Point doesn't lie in Feasible Region"

	CALL Transform(A,b,c, A_art,c_art, x0)

	x_art = Karmarkar(A_art, c_art)

	x_opt = InvTransform(x_art, x0)

ELSE

	CALL StdToCan(A,b,c, A_can,c_can)

	x_can = Karmarkar(A_can,c_can)

	x = InvTransform(x_can)

	x_opt(:) = x(:size(A,2))

END IF

END FUNCTION LPSolver
