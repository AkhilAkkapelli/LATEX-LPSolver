SUBROUTINE Canonical

UREAL :: Acan(row,column), ccan(column), xopt(column)
INTEGER :: i, j


OPEN(unit = 10, file= './1.IO/Acan.txt')
  DO i=1,row
    READ(10, *) (Acan(i,j) ,j=1,column)
  END DO
CLOSE(10) 

OPEN(unit = 30, file= './1.IO/ccan.txt')
  DO j=1,column
    READ(30, *) ccan(j)
  END DO
CLOSE(30)


!CALL ProbGen(Acan)
!CALL ProbGen(ccan)

xopt = ProjectiveScale(Acan, ccan, potentialstop, potentialratio)

END SUBROUTINE Canonical
