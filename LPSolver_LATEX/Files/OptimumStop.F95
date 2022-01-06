LOGICAL PURE FUNCTION optimumstop(n, x,xp, iter, c) RESULT(stp)

INTEGER, INTENT(IN) :: iter, n
UREAL, INTENT(IN) :: x(n), xp(n), c(n)

INTEGER, PARAMETER :: iterlimit = 10000
UREAL :: obj, objp

stp = .FALSE. 
obj = DOT(c,x)
objp = DOT(c,xp)

IF(iter >= iterlimit .OR. objp - obj < real(1.Q-100,realKind) ) stp = .TRUE.

END FUNCTION optimumstop
