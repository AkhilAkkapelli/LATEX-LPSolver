alpha = (1.Q0-beta)/(n*c_unit(n-1))
beta = real(1.Q-1,realKind)

DO idx=1,n
IF(c_unit(idx)<0 .OR. idx == n-1) CYCLE
alpha = min((1.Q0-beta)/(n*c_unit(idx)),alpha)
END DO

x = x0 - alpha*c_unit
