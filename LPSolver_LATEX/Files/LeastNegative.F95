SUBROUTINE LeastNegative

UREAL :: Aln(row, column), bln(row), cln(column), y0(row), x0(column), xln(column), eps, e(column), &
		      Aaug(row,column+1), baug(row), caug(column+1), x0aug(column+1), xaug(column+1), &
		    	Aaugcan(size(Aln,1),size(Aln,2)+2), caugcan(size(cln)+2), xaugcan(size(Aln,2)+2), &
		  		Acan(size(Aln,1),size(Aln,2)+1), ccan(size(cln)+1), xcan(size(Aln,2)+1)	
			
INTEGER :: i, j


CALL ProbGen(Aln)
CALL ProbGen(bln)
CALL ProbGen(cln)

e=one

y0 = PDSOL(GEMM(Aln,TRANS(Aln)), bln)

x0 = GEMV(TRANS(Aln), y0)

IF(ANY(x0 < 0)) THEN

eps = -2*minval(x0)

Aaug = Aln .HAUG. -GEMV(Aln,e)

baug = bln 

x0aug = (x0 + eps) .VAUG. eps

caug = 0
caug(column+1) = 1

CALL ProjectiveTransform(Aaug,baug,caug, Aaugcan, caugcan, x0aug)

xaugcan = ProjectiveScale(Aaugcan, caugcan, potentialstop, potentialratio)

xaug = InvProjTransform(xaugcan, x0aug)

x0 = xaug(:column)

END IF

IF(xaug(column+1) < 1.0Q-30) THEN

CALL ProjectiveTransform(Aln,bln,cln, Acan,ccan, x0)

ccan = ccan - MAX(MAXVAL(ccan),0.)

xcan = ProjectiveScale(Acan, ccan, potentialstop, potentialratio)

xln = InvProjTransform(xcan, x0)

ELSE

eps = xaug(column +1)

END IF

END SUBROUTINE LeastNegative
