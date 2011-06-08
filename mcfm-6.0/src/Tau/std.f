      SUBROUTINE STD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX SPL(10,10),SMN(10,10),C23(10)
      DOUBLE PRECISION ROOT(10)
      COMMON/CSTD/SPL,SMN
      COMMON/MOM/PLAB(4,10)
      DO 5 K=1,10
      SPL(K,K)=(0.D0,0.D0)
      SMN(K,K)=SPL(K,K)
      ROOT(K)=DSQRT(PLAB(4,K)-PLAB(1,K))
   5  C23(K)=DCMPLX(PLAB(2,K),PLAB(3,K))
      DO 10 I=2,10
      DO 10 J=1,I-1
      SPL(I,J)=C23(I)*ROOT(J)/ROOT(I) - C23(J)*ROOT(I)/ROOT(J)
      SPL(J,I)=-SPL(I,J)
      SMN(I,J)=-DCONJG(SPL(I,J))
      SMN(J,I)=-SMN(I,J)
  10  CONTINUE
      RETURN
      END
C
