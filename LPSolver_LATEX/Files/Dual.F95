SUBROUTINE Dual(A, b, c, A_dual, b_dual, c_dual)

REAL, INTENT(IN) :: A(:,:), b(:), c(:)

REAL, INTENT(OUT) :: A_dual(size(A,2),size(A,1)), b_dual(size(c)), c_dual(size(b))


c_dual = -b

A_dual = -TRANSPOSE(A)

b_dual = -c

END SUBROUTINE Dual
