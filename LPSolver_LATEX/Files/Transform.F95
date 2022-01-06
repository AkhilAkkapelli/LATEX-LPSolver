SUBROUTINE Transform(A,b,c, A_can,c_can, a0)

REAL, INTENT(IN) :: A(:,:),b(:),c(:), a0(:)

REAL, INTENT(OUT) :: A_can(size(A,1),size(A,2)+1), c_can(size(c)+1)


A_can(:size(A,1),:size(A,2)) = ColMult(a0,A)
A_can(:,size(A,2)+1) = -b

c_can(:size(c)) = a0*c
c_can(size(c)+1) = 0.

END SUBROUTINE Transform
