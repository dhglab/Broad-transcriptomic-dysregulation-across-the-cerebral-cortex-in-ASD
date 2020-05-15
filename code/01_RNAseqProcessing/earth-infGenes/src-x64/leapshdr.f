C Copied from Thomas Lumley's leaps 2.9 package for earth 3.2-6 to avoid
C use of external routines in earth causing complaints from CRAN check.
C Original Fortran code by Alan Miller.
C
      SUBROUTINE MAKEQR(NP,NN,WEIGHTS,TXMAT,YVEC,D,RBAR,THETAB,
     $     SSERR,IER)
C     Calls INCLUD to construct Banachiewicz factorisation
C
C
      INTEGER NP, NN, IER
      DOUBLE PRECISION WEIGHTS(NN), TXMAT(*), YVEC(NN), D(NP), RBAR(*),
     +     THETAB(NP), SSERR
C     local variables
      INTEGER I, NRBAR
      IER=0
      NRBAR=(NP*(NP-1))/2
      DO 10 I=0, NN-1
         CALL INCLUD(NP,NRBAR,WEIGHTS(I+1),TXMAT(I*NP+1),YVEC(I+1),D,
     $        RBAR,THETAB, SSERR,IER)
         IF (IER .NE. 0) RETURN
 10   CONTINUE
      RETURN
      END
