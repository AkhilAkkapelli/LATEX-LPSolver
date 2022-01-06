SUBROUTINE Equality

UREAL :: Aeq(row, column), beq(row), ceq(column), x0(column), xeq(column), eps, &
			Acan(size(Aeq,1),size(Aeq,2)+1), ccan(size(ceq)+1), xcan(size(Aeq,2)+1), &
			Acan1(size(Aeq,1),size(Aeq,2)+2), ccan1(size(ceq)+2), xcan1(size(Aeq,2)+2), &
				A(row,column+1), b(row), c(column+1), x01(column+1), xopt(column), xopt1(column+1)
INTEGER :: i, j


OPEN(unit = 10, file= './1.IO/Aeq.txt')
DO i=1,row
READ(10, *) (Aeq(i,j) ,j=1,column)
END DO
CLOSE(10)  

OPEN(unit = 20, file= './1.IO/beq.txt')
DO i=1,row
     READ(20, *) beq(i)
END DO
CLOSE(20)

OPEN(unit = 30, file= './1.IO/ceq.txt')
DO j=1,column
     READ(30, *) ceq(j)
END DO
CLOSE(30) 

x0 = one

x01 = x0 .VAUG. one

b = beq + 0*(beq - GEMV(Aeq,x0))

A(:,:column) = Aeq
A(:,column+1) = 1*(beq - GEMV(Aeq,x0))

c(:column) = 0
c(column+1) = one

CALL ProjectiveTransform(A,b,c, Acan1,ccan1, x01)

!ccan = ccan - MAX(MAXVAL(ccan),0.)

xcan1 = ProjectiveScale(Acan1, ccan1, optimumstop, potentialratio)

xopt1 = InvProjTransform(xcan1, x01)

x0 = xopt1(:column) 
print*, 'x0', xopt1
IF(x0(column+1) < 1.Q0-10) THEN

CALL ProjectiveTransform(Aeq,beq,ceq, Acan,ccan, x0)

xcan = ProjectiveScale(Acan, ccan, potentialstop, potentialratio)

xopt = InvProjTransform(xcan, x0)
print*, 'xopt', xopt
END IF

END SUBROUTINE Equality
