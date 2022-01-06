SUBROUTINE StdToCan(A,b,c ,A_can,c_can)

REAL, INTENT(IN) :: A(:,:), b(:), c(:)
REAL, INTENT(OUT) :: A_can(size(A,1)+size(A,2)+1,2*(size(A,1)+size(A,2)+1)), c_can(2*(size(A,1)+size(A,2)+1))

REAL :: A_pd(size(A,1)+size(A,2)+1,size(A,1)+size(A,2)), &
		b_pd(size(b)+size(c)+1),c_pd(size(c)+size(b))
REAL :: A_slk(size(A_pd,1),2*size(A_pd,2)), &
		b_slk(size(b_pd)),c_slk(2*size(c_pd))
REAL :: A_art(size(A_slk,1),size(A_slk,2)+1),b_art(size(b_slk)),c_art(size(c_slk)+1), x0(size(c_slk)+1)

INTEGER :: m,n


m = size(A,1); n = size(A,2)

!!!STEP 1 Combine Primal and Dual problems
CALL Combine(A,b,c, A_pd,b_pd,c_pd, m,n)


!!!STEP 2 Introduce Slack Variables
CALL Slack(A_pd,b_pd,c_pd, A_slk,b_slk,c_slk, m,n)


!!!STEP 3 Introdue an Artificial variable to create an interior starting point.
CALL Artificial(A_slk,b_slk,c_slk, A_art,b_art,c_art, m,n, x0)


!!!STEP 4 A Projective transformation of the positive orthant into a simplex.
CALL Transform(A_art,b_art,c_art, A_can,c_can, x0)

END SUBROUTINE StdToCan
