SUBROUTINE Artificial(A_slk,b_slk,c_slk, A_art,b_art,c_art, m,n, a0)

REAL, INTENT(IN) :: A_slk(:,:),b_slk(:),c_slk(:)
INTEGER, INTENT(IN) :: m,n

REAL, INTENT(OUT) :: A_art(m+n+1,2*m+2*n+1),b_art(m+n+1),c_art(2*m+2*n+1), a0(2*m+2*n+1)

REAL :: x0(n), y0(m), u0(m), v0(n), lambda_0


x0 = 1.; y0 = 1.; u0 = 1.; v0 = 1.; lambda_0 = 1.

A_art(:,:2*(m+n)) = A_slk
A_art(:m,2*m+2*n+1) = -A_slk(m+n+1,n+1:m+n) - GEMV(A_slk(:m,:n),x0) + y0
A_art(m+1:m+n,2*m+2*n+1) =  A_slk(m+n+1,:n) - GEMV(A_slk(m+1:m+n,n+1:m+n),u0) - v0
A_art(m+n+1,2*m+2*n+1) = DOT((-A_slk(m+n+1,:n) .VAUG. -A_slk(m+n+1,n+1:m+n)),(x0 .VAUG. u0))

b_art = 0.
b_art(:size(b_slk)) = b_slk

c_art = 0.
c_art(2*m+2*n+1) = 1.

a0 = x0 .VAUG. y0 .VAUG. u0 .VAUG. v0 .VAUG. lambda_0

END SUBROUTINE Artificial
