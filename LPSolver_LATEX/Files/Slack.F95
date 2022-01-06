SUBROUTINE Slack(A_pd,b_pd,c_pd, A_slk,b_slk,c_slk, m,n)

REAL, INTENT(IN) :: A_pd(:,:),b_pd(:),c_pd(:)
INTEGER, INTENT(IN) :: m,n

REAL, INTENT(OUT) :: A_slk(m+n+1,2*m+2*n), b_slk(m+n+1), c_slk(2*m+2*n)


A_slk = 0.
A_slk(:,:m+n) = A_pd
A_slk(:m,m+n+1:2*m+n) = -Ones(m)
A_slk(m+1:m+n,2*m+n+1:) = Ones(n)

b_slk = b_pd

c_slk=0.
c_slk(:size(c_pd)) = c_pd

END SUBROUTINE Slack
