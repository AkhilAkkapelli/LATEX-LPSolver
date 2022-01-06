#include "Preprocessor.F90"

MODULE AUGOperator
IMPLICIT NONE
!!! AUGMENT OPERATION on ANY SCALAR/VECTOR/MATRIX AND ANY SCALAR/VECTOR/MATRIX


PRIVATE :: sHAugv, vHaugs, vHAugv, MHAugs, sHAugM, MHAugv, vHAugM, MHAugM, sVAugs, sVAugv, vVAugs, vVAugv, vVAugM, MVAugv, MVAugM

PUBLIC  :: OPERATOR(.HAUG.), OPERATOR(.VAUG.)
! .HAUG. Operates Horizontal Augmentation and .VAUG. Operates Vertical Augmentation


INTERFACE OPERATOR(.HAUG.)
  MODULE PROCEDURE sHAugv
  MODULE PROCEDURE vHAugs
  MODULE PROCEDURE vHAugv
  MODULE PROCEDURE MHAugs
  MODULE PROCEDURE sHAugM
  MODULE PROCEDURE MHAugv
  MODULE PROCEDURE vHAugM
  MODULE PROCEDURE MHAugM
END INTERFACE

INTERFACE OPERATOR(.VAUG.)
  MODULE PROCEDURE sVAugs
  MODULE PROCEDURE sVAugv
  MODULE PROCEDURE vVAugs
  MODULE PROCEDURE vVAugv
  MODULE PROCEDURE vVAugM
	MODULE PROCEDURE MVAugv
	MODULE PROCEDURE MVAugM
END INTERFACE


CONTAINS


PURE FUNCTION sHAugv(s,v) RESULT(sv)
! Horizontal Augment SCALAR and VECTOR : sv = [s v]

UREAL, INTENT(IN) :: s, v(:)
UREAL             :: sv(size(v),2)
	
sv = 0
sv(1,1) = s
sv(:,2) = v(:)

END FUNCTION sHAugv

PURE FUNCTION vHAugs(v,s) RESULT(vs)
! Horizontal Augment VECTOR and SCALAR : vs = [v s]

UREAL, INTENT(IN) :: v(:), s
UREAL             :: vs(size(v),2)
	

vs = 0
vs(:,1) = v(:)
vs(1,2) = s

END FUNCTION vHAugs

PURE FUNCTION vHAugv(u,v) RESULT(uv)
! Horizontal Augment VECTOR and VECTOR : uv = [u v]

UREAL, INTENT(IN) :: u(:), v(:)
UREAL             :: uv(max(size(u),size(v)),2)

	
uv = 0
uv(:size(u),1) = u(:)
uv(:size(v),2) = v(:)

END FUNCTION vHAugv

PURE FUNCTION MHAugs(M,s) RESULT(Ms)
! Horizontal Augment MATRIX and SCALAR : Ms = [M s]

UREAL, INTENT(IN) :: M(:,:), s
UREAL             :: Ms(size(M,1),size(M,2)+1)


Ms = 0
Ms(:,:size(M,2)) = M
Ms(1,size(M,2)+1) = s

END FUNCTION MHAugs

PURE FUNCTION sHAugM(s,M) RESULT(sM)
! Horizontal Augment SCALAR and MATRIX : sM = [s M]

UREAL, INTENT(IN) :: s, M(:,:)
UREAL             :: sM(size(M,1),size(M,2)+1)


sM = 0
sM(1,1) = s
sM(:,2:) = M	

END FUNCTION sHAugM

PURE FUNCTION MHAugv(M,v) RESULT(Mv)
! Horizontal Augment MATRIX and VECTOR : Mv = [M v]

UREAL, INTENT(IN) :: M(:,:), v(:)
UREAL             :: Mv(max(size(M,1),size(v)),size(M,2)+1)


Mv = 0
Mv(:size(M,1),:size(M,2)) = M
Mv(:size(v),size(M,2)+1) = v

END FUNCTION MHAugv

PURE FUNCTION vHAugM(v,M) RESULT(vM)
! Horizontal Augment VECTOR and MATRIX : vM = [v M]

UREAL, INTENT(IN) :: v(:), M(:,:)
UREAL             :: vM(max(size(v),size(M,1)),size(M,2)+1)


vM = 0
vM(:size(M),1) = v
vM(:size(M,1),2:) = M	

END FUNCTION vHAugM

PURE FUNCTION MHAugM(M,N) RESULT(MN)
! Horizontal Augment MATRIX and MATRIX : MN = [M N]

UREAL, INTENT(IN) :: M(:,:), N(:,:)
UREAL             :: MN(max(size(M,1),size(N,1)),size(M,2)+size(N,2))


MN = 0
MN(:size(M,1),:size(M,2)) = M
MN(:size(N,1),size(M,2)+1:size(M,2)+size(N,2)) = N

END FUNCTION MHAugM


PURE FUNCTION sVAugs(s,t) RESULT(st)
! Vertical Augment SCALAR and SCALAR : st = [s t]^T

UREAL, INTENT(IN) :: s,t
UREAL             :: st(2)


st(1) = s
st(2) = t

END FUNCTION sVAugs

PURE FUNCTION sVAugv(s,v) RESULT(sv)
! Vertical Augment SCALAR and VECTOR : sv = [s v]^T

UREAL, INTENT(IN) :: s, v(:)
UREAL             :: sv(1+size(v))


sv(1) = s
sv(2:) = v

END FUNCTION sVAugv

PURE FUNCTION vVAugs(v,s) RESULT(vs)
! Vertical Augment VECTOR and SCALAR : vs = [v s]^T

UREAL, INTENT(IN) :: v(:), s
UREAL             :: vs(size(v)+1)


vs(:size(v)) = v
vs(size(v)+1) = s

END FUNCTION vVAugs


PURE FUNCTION vVAugv(u,v) RESULT(uv)
! Vertical Augment VECTOR and VECTOR : uv = [u v]^T

UREAL, INTENT(IN) :: u(:), v(:)
UREAL             :: uv(size(u)+size(v))


uv(:size(u)) = u
uv(size(u)+1:) = v

END FUNCTION vVAugv

PURE FUNCTION vVAugM(v,M) RESULT(vM)
! Vertical Augment VECTOR and MATRIX : vM = [v^T M]^T

UREAL, INTENT(IN) :: v(:), M(:,:)
UREAL :: vM(size(M,1)+1,max(size(v),size(M,2)))


vM = 0
vM(1,:size(v)) = v
vM(2: ,:size(M,2)) = M

END FUNCTION vVAugM

PURE FUNCTION MVAugv(M,v) RESULT(Mv)
! Vertical Augment VECTOR and MATRIX : Mv = [M v^T]^T

UREAL, INTENT(IN) :: M(:,:), v(:)
UREAL :: Mv(size(M,1)+1,max(size(M,2),size(v)))


Mv = 0
Mv(:,:size(M,2)) = M
Mv(size(M,1)+1,:size(v)) = v

END FUNCTION MVAugv

PURE FUNCTION MVAugM(M,N) RESULT(MN)
! Vertical Augment VECTOR and MATRIX : MN = [M N]^T

UREAL, INTENT(IN) :: M(:,:), N(:,:)
UREAL :: MN(size(M,1)+size(N,1),max(size(M,2),size(N,2)))


MN = 0
MN(:size(M,1),:size(M,2)) = M
MN(size(M,1)+1:size(M,1)+size(N,1),:size(N,2)) = N

END FUNCTION MVAugM

END MODULE AUGOperator
