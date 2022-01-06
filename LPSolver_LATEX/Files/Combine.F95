SUBROUTINE Combine(A_std,b_std,c_std, A_pd,b_pd,c_pd, m,n)

REAL, INTENT(IN) :: A_std(:,:),b_std(:),c_std(:)
INTEGER, INTENT(IN) :: m,n

REAL, INTENT(OUT) :: A_pd(m+n+1,m+n), b_pd(m+n+1), c_pd(m+n)


A_pd = 0.
A_pd(:m,:n) = A_std
A_pd(m+1:m+n,n+1:m+n) = TRANSPOSE(A_std)
A_pd(m+n+1,:) = c_std .VAUG. -b_std

b_pd(:m) = b_std
b_pd(m+1:m+n) = c_std
b_pd(m+n+1) = 0

c_pd = 0

END SUBROUTINE Combine
