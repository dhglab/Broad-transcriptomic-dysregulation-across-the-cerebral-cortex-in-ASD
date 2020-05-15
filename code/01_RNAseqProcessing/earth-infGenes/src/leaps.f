C Copied from Thomas Lumley's leaps 2.9 package for earth 3.2-6 to avoid
C use of external routines in earth causing complaints from CRAN check.
C Code which writes to the console has been commented out.
C Original Fortran code from Alan Miller.
C
C      PROGRAM START
C
C     This is the starting program for the SUBSETS package of programs.
C     It forms the upper-triangular Banachiewicz factorization of the
C     input data.
C     Free-format input is assumed, i.e. with data fields separated by
C     spaces, CR's, tabs or commas.   N.B. Some Fortran compilers will
C     not accept tabs and/or commas as delimiters.
C     Warning: Some Fortran compilers will not allow free format input
C     of character data.   This program inputs the names of variables
C     in free format.
C
C     Latest revision - 16 August 1992
C
C     Stephen Milborrow 5 Sep 2016: tweaked code to eliminate warnings:
C       Warning: Obsolescent feature: DO termination statement which is
C       not END DO or CONTINUE with label
C
c___c     IMPLICIT NONE
c___      integer npmax, dimu
c___      parameter (npmax=50, dimu=npmax*(npmax+1)/2)
c___      DOUBLE PRECISION U(dimu), EL(0:npmax), RHS(0:npmax), X(0:npmax),
c___     +                 WT, ONE, Y, RESSQ
c___      CHARACTER ANS, FNAME*20, VNAME(0:npmax)*8, YNAME*8, TEXT*79
c___      INTEGER LIN, YPOS, IPOS, I, K, ICONST, NCOLS, NOBS, NRBAR, IER,
c___     +        LINE1, LOUT
c___      LOGICAL OK, LSEL
c___      DATA WT/1.D0/, ONE/1.D0/, LSEL/.FALSE./
c___
c___C
c___C     Set unit numbers for I/O in the data statement below.
c___C
c___      DATA LIN/5/, LOUT/6/
c___C
c___C     Ask for details of the data file.
c___C
c___   10 WRITE(LOUT, 900)
c___  900 FORMAT(' Name of data file = ? ')
c___      READ(LIN, *) FNAME
c___C
c___C     Add extension .dat if none has been entered,
c___C     detected by the lack of a '.'
c___C
c___      IF (INDEX(FNAME, '.') .EQ. 0) THEN
c___    IPOS = INDEX(FNAME, ' ')
c___    FNAME = FNAME(1:IPOS-1) // '.dat'
c___      END IF
c___C
c___C     Check that file exists.
c___C
c___      INQUIRE(FILE=FNAME, EXIST=OK)
c___      IF (.NOT. OK) THEN
c___    WRITE(*, 910) FNAME
c___  910   FORMAT(' *** File not found - ', a, ' **')
c___    GO TO 10
c___      END IF
c___C
c___C     Display first part of file.
c___C
c___      OPEN(10, FILE=FNAME, STATUS='OLD')
c___      WRITE(*, *)'Start of your data file follows'
c___      DO 20 I = 1, 12
c___    READ(10, '(A)') TEXT
c___    WRITE(*, '(1X, A)') TEXT
c___   20 CONTINUE
c___      REWIND 10
c___C
c___      WRITE(LOUT, 920)
c___  920 FORMAT(' How many X-variables ? ')
c___      READ(LIN, *) K
c___      WRITE(LOUT, 930)
c___  930 FORMAT('+Do you want a constant in the model ? ')
c___      READ(LIN, *) ANS
c___      ICONST = 0
c___      IF(ANS.EQ.'Y' .OR. ANS .EQ. 'y') ICONST = 1
c___      NCOLS = K + ICONST
c___      NRBAR = NCOLS * (NCOLS - 1) / 2
c___C
c___C     Get position of dependant variable.
c___C
c___      WRITE(*, *)'Is dependant variable at end ? (Y/N): '
c___      READ(*, *) ANS
c___      IF (ANS .EQ. 'Y' .OR. ANS .EQ. 'y') THEN
c___    YPOS = K+1
c___      ELSE
c___    WRITE(*, *)'Enter no. of position of dependant variable: '
c___    READ(*, *) YPOS
c___    IF (YPOS .LT. 1) YPOS = 1
c___    IF (YPOS .GT. K) YPOS = K + 1
c___      END IF
c___C
c___C     Enter variable names, read them from file, or set defaults.
c___C
c___      VNAME(0) = 'Constant'
c___      WRITE(*, *)'Are variable names in data file ? (Y/N): '
c___      READ(*, *) ANS
c___      IF (ANS .EQ. 'Y' .OR. ANS .EQ. 'y') THEN
c___    WRITE(*, *)'Which line do names start on ? '
c___    READ(*, *) LINE1
c___    IF (LINE1 .GT. 1) THEN
c___      DO 30 I = 1, LINE1-1
c___   30     READ(10, *)
c___    END IF
c___    IF (YPOS .GT. K) THEN
c___      READ(10, *) (VNAME(I),I=1,K), YNAME
c___    ELSE IF (YPOS .EQ. 1) THEN
c___      READ(10, *) YNAME, (VNAME(I),I=1,K)
c___    ELSE
c___      READ(10, *) (VNAME(I),I=1,YPOS-1), YNAME,
c___     +                        (VNAME(I),I=YPOS,K)
c___    END IF
c___    REWIND 10
c___      ELSE
c___    WRITE(*, *)'Do you want to name variables ? (Y/N): '
c___    READ(*, '(a)') ANS
c___    IF (ANS .EQ. 'Y' .OR. ANS .EQ. 'y') THEN
c___      WRITE(*, *)'Variable names may contain up to 8 characters'
c___      WRITE(*, *)'Name for dependant (Y) variable = ? '
c___      READ(*, '(a)') YNAME
c___      DO 40 I = 1, K
c___        WRITE(*, *)'Name for variable', I, ' = ? '
c___        READ(*, '(a)') VNAME(I)
c___   40     CONTINUE
c___    ELSE
c___      DO 50 I = 1, K
c___        WRITE(VNAME(I), 940) I
c___  940       FORMAT('XVAR(', I2, ')')
c___   50     CONTINUE
c___      YNAME = 'Dept.var'
c___    END IF
c___      END IF
c___C
c___      WRITE(*, *)'Which line does the data start on ? '
c___      READ(*, *) LINE1
c___      IF (LINE1 .GT. 1) THEN
c___    DO 60 I = 1, LINE1-1
c___   60   READ(10, *)
c___      END IF
c___C
c___C     Read in data and form the upper-triangular factorization.
c___C
c___      IF (ICONST .EQ. 1) THEN
c___    CALL CLEAR(NCOLS, NRBAR, EL, U, RHS, RESSQ, IER)
c___      ELSE
c___    CALL CLEAR(NCOLS, NRBAR, EL(1), U, RHS(1), RESSQ, IER)
c___      END IF
c___      NOBS = 1
c___      X(0) = ONE
c___C
c___C     Case is skipped if spurious characters are found (e.g. for
c___C     missing values).
c___C
c___   70 CONTINUE
c___      IF (YPOS .GT. K) THEN
c___    READ(10, *, ERR=70, END=80) (X(I),I=1,K), Y
c___      ELSE IF (YPOS .EQ. 1) THEN
c___    READ(10, *, ERR=70, END=80) Y, (X(I),I=1,K)
c___      ELSE
c___    READ(10, *, ERR=70, END=80) (X(I),I=1,YPOS-1), Y,
c___     +                              (X(I),I=YPOS,K)
c___      END IF
c___      IF (ICONST .EQ. 1) THEN
c___    CALL INCLUD(NCOLS, NRBAR, WT, X, Y, EL, U, RHS, RESSQ, IER)
c___      ELSE
c___    CALL INCLUD(NCOLS, NRBAR, WT, X(1), Y, EL(1), U, RHS(1), RESSQ,
c___     +              IER)
c___      END IF
c___      NOBS = NOBS + 1
c___      GO TO 70
c___C
c___C     Change extension to .red for output file.
c___C
c___   80 IPOS = INDEX(FNAME, '.')
c___      FNAME(IPOS+1:IPOS+3) = 'red'
c___      NOBS = NOBS - 1
c___C
c___C     Write U, EL, RHS & the residual sum of squares (RESSQ) to disk.
c___C
c___      OPEN(9, FILE=FNAME, STATUS='NEW', ACCESS='SEQUENTIAL',
c___     +          FORM='UNFORMATTED')
c___      WRITE(9) K, ICONST, NCOLS, NOBS, NRBAR, LSEL
c___      IF (ICONST .EQ. 0) THEN
c___    WRITE(9) YNAME, (VNAME(I),I=1,K)
c___    WRITE(9) (U(I),I=1,NRBAR), (EL(I),I=1,K), (RHS(I),I=1,K), RESSQ
c___      ELSE
c___    WRITE(9) YNAME, (VNAME(I),I=0,K)
c___    WRITE(9) (U(I),I=1,NRBAR), (EL(I),I=0,K), (RHS(I),I=0,K), RESSQ
c___      END IF
c___      ENDFILE 9
c___C
c___      END


      SUBROUTINE CLEAR(NP, NRBAR, D, RBAR, THETAB, SSERR, IER)
C
C     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
C
C     Sets arrays to zero prior to calling INCLUD
C
      INTEGER NP, NRBAR, IER
      DOUBLE PRECISION D(NP), RBAR(*), THETAB(NP), SSERR
C
C     Local variables
C
      INTEGER I
      DOUBLE PRECISION ZERO
C
      DATA ZERO/0.D0/
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (IER .NE. 0) RETURN
C
      DO 10 I = 1, NP
        D(I) = ZERO
        THETAB(I) = ZERO
   10 CONTINUE
      DO 20 I = 1, NRBAR
      RBAR(I) = ZERO
   20 END DO
      SSERR = ZERO
      RETURN
      END


c___      SUBROUTINE INCLUD(NP, NRBAR, WEIGHT, XROW, YELEM, D,
c___     +      RBAR, THETAB, SSERR, IER)
c___C
c___C     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
c___C     Modified from algorithm AS 75.1
c___C
c___C     Calling this routine updates d, rbar, thetab and sserr by the
c___C     inclusion of xrow, yelem with the specified weight.   The number
c___C     of columns (variables) may exceed the number of rows (cases).
c___C
c___C**** WARNING: The elements of XROW are overwritten  ****
c___C
c___      INTEGER NP, NRBAR, IER
c___      DOUBLE PRECISION WEIGHT, XROW(NP), YELEM, D(NP), RBAR(*),
c___     +    THETAB(NP), SSERR
c___C
c___C     Local variables
c___C
c___      INTEGER I, K, NEXTR
c___      DOUBLE PRECISION ZERO, W, Y, XI, DI, WXI, DPI, CBAR, SBAR, XK
c___C
c___      DATA ZERO/0.D0/
c___C
c___C     Some checks.
c___C
c___      IER = 0
c___      IF (NP .LT. 1) IER = 1
c___      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
c___      IF (IER .NE. 0) RETURN
c___C
c___      W = WEIGHT
c___      Y = YELEM
c___      NEXTR = 1
c___      DO 30 I = 1, NP
c___C
c___C     Skip unnecessary transformations.   Test on exact zeros must be
c___C     used or stability can be destroyed.
c___C
c___    IF (W .EQ. ZERO) RETURN
c___    XI = XROW(I)
c___    IF (XI .EQ. ZERO) THEN
c___      NEXTR = NEXTR + NP - I
c___      GO TO 30
c___    END IF
c___    DI = D(I)
c___    WXI = W * XI
c___    DPI = DI + WXI*XI
c___    CBAR = DI / DPI
c___    SBAR = WXI / DPI
c___    W = CBAR * W
c___    D(I) = DPI
c___    IF (I .EQ. NP) GO TO 20
c___    DO 10 K = I+1, NP
c___      XK = XROW(K)
c___      XROW(K) = XK - XI * RBAR(NEXTR)
c___      RBAR(NEXTR) = CBAR * RBAR(NEXTR) + SBAR * XK
c___      NEXTR = NEXTR + 1
c___   10   CONTINUE
c___   20   XK = Y
c___    Y = XK - XI * THETAB(I)
c___    THETAB(I) = CBAR * THETAB(I) + SBAR * XK
c___   30 CONTINUE
c___C
c___C     Y * SQRT(W) is now equal to Brown & Durbin's recursive residual.
c___C
c___      SSERR = SSERR + W * Y * Y
c___C
c___      RETURN
c___      END
      SUBROUTINE ADD1(NP, NRBAR, D, RBAR, THETAB, FIRST, LAST, TOL, SS,
     +         SXX, SXY, SMAX, JMAX, IER)
C
C     Calculate the reduction in residual sum of squares when one
C     variable, selected from those in positions FIRST .. LAST, is
C     added, given that the variables in positions 1 .. FIRST-1 (if
C     any) are already included.
C
      INTEGER NP, NRBAR, FIRST, LAST, JMAX, IER
      DOUBLE PRECISION D(NP), RBAR(NRBAR), THETAB(NP), TOL(NP), SS(NP),
     +        SXX(NP), SXY(NP), SMAX
C
C     Local variables
C
      INTEGER J, INC, POS, ROW, COL
      DOUBLE PRECISION ZERO, DIAG, DY, SSQX
      DATA ZERO/0.D0/
C
C     Check call arguments
C
      JMAX = 0
      SMAX = ZERO
      IER = 0
      IF (FIRST .GT. NP) IER = 1
      IF (LAST .LT. FIRST) IER = IER + 2
      IF (FIRST .LT. 1) IER = IER + 4
      IF (LAST .GT. NP) IER = IER + 8
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 16
      IF (IER .NE. 0) RETURN
C
C     Accumulate sums of squares & products from row FIRST
C
      DO 10 J = FIRST, LAST
        SXX(J) = ZERO
        SXY(J) = ZERO
   10 CONTINUE
      INC = NP - LAST
      POS = (FIRST-1) * (NP+NP-FIRST)/2 + 1
      DO 30 ROW = FIRST, LAST
        DIAG = D(ROW)
        DY = DIAG * THETAB(ROW)
        SXX(ROW) = SXX(ROW) + DIAG
        SXY(ROW) = SXY(ROW) + DY
        DO 20 COL = ROW+1, LAST
          SXX(COL) = SXX(COL) + DIAG * RBAR(POS)**2
          SXY(COL) = SXY(COL) + DY * RBAR(POS)
          POS = POS + 1
   20   CONTINUE
        POS = POS + INC
   30 CONTINUE
C
C     Incremental sum of squares for a variable = SXY * SXY / SXX.
C     Calculate whenever sqrt(SXX) > TOL for that variable.
C
      DO 40 J = FIRST, LAST
        SSQX = SXX(J)
        IF (SQRT(SSQX) .GT. TOL(J)) THEN
          SS(J) = SXY(J)**2 / SXX(J)
          IF (SS(J) .GT. SMAX) THEN
            SMAX = SS(J)
            JMAX = J
          END IF
        ELSE
          SS(J) = ZERO
        END IF
   40 CONTINUE
C
      RETURN
      END

      SUBROUTINE BAKWRD(NP, NRBAR, D, RBAR, THETAB, FIRST, LAST,
     *   VORDER, TOL, RSS, BOUND, NVMAX, RESS, IR, NBEST, LOPT, IL,
     *   WK, IWK, IER)
C
C     Backward elimination from variables in positions FIRST .. LAST.
C     If FIRST > 1, variables in positions prior to this are forced in.
C     If LAST < NP, variables in positions after this are forced out.
C     On exit, the array VORDER contains the numbers of the variables
C     in the order in which they were deleted.
C
      INTEGER NP, NRBAR, FIRST, LAST, VORDER(NP), NVMAX, IR, NBEST,
     *   IL, LOPT(IL, *), IWK, IER
      DOUBLE PRECISION D(NP), RBAR(NRBAR), THETAB(NP), TOL(NP), RSS(NP),
     *   BOUND(NVMAX), RESS(IR, *), WK(IWK)
C
C     Local variables
C
      INTEGER NEED, POS, J1, JMIN, I
      DOUBLE PRECISION SMIN
C
C     Check call arguments
C
      IER = 0
      IF (FIRST .GE. NP) IER = 1
      IF (LAST .LE. 1) IER = IER + 2
      IF (FIRST .LT. 1) IER = IER + 4
      IF (LAST .GT. NP) IER = IER + 8
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 16
      IF (IWK .LT. 2*LAST) IER = IER + 32
      IF (NBEST .GT. 0) THEN
        NEED = NVMAX*(NVMAX+1)/2
        IF (IR .LT. NVMAX) IER = IER + 64
        IF (IL .LT. NEED) IER = IER + 128
      END IF
      IF (IER .NE. 0) RETURN
C
C     For POS = LAST, ..., FIRST+1 call DROP1 to find best variable to
C     find which variable to drop next.
C
      J1 = LAST + 1
      DO 20 POS = LAST, FIRST+1, -1
        CALL DROP1(NP, NRBAR, D, RBAR, THETAB, FIRST, POS, TOL, WK,
     *            WK(J1), SMIN, JMIN, IER)
        IF (JMIN .GT. 0 .AND. JMIN .LT. POS) THEN
          CALL VMOVE(NP, NRBAR, VORDER, D, RBAR, THETAB, RSS, JMIN, POS,
     *            TOL, IER)
          IF (NBEST .GT. 0) THEN
            DO 10 I = JMIN, POS-1
            CALL REPORT(I, RSS(I), BOUND, NVMAX, RESS, IR, NBEST, LOPT,
     *            IL, VORDER)
   10       END DO
          END IF
        END IF
   20 CONTINUE
C
      RETURN
      END

      SUBROUTINE DROP1(NP, NRBAR, D, RBAR, THETAB, FIRST, LAST, TOL,
     *    SS, WK, SMIN, JMIN, IER)
C
C     Calculate the increase in the residual sum of squares when
C     variable J is dropped from the model, for J = FIRST, ..., LAST.
C
      INTEGER NP, NRBAR, FIRST, LAST, JMIN, IER
      DOUBLE PRECISION D(NP), RBAR(NRBAR), THETAB(NP), TOL(NP),
     *    SS(LAST), WK(LAST), SMIN
C
C     Local variables
C
      INTEGER J, POS1, INC, POS, ROW, COL, I
      DOUBLE PRECISION LARGE, ZERO, D1, RHS, D2, X
      DATA LARGE/1.D+35/, ZERO/0.D0/
C
C     Check call arguments
C
      JMIN = 0
      SMIN = LARGE
      IER = 0
      IF (FIRST .GT. NP) IER = 1
      IF (LAST .LT. FIRST) IER = IER + 2
      IF (FIRST .LT. 1) IER = IER + 4
      IF (LAST .GT. NP) IER = IER + 8
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 16
      IF (IER .NE. 0) RETURN
C
C     POS1 = position of first element of row FIRST in RBAR.
C
      POS1 = (FIRST - 1) * (NP + NP - FIRST)/2 + 1
      INC = NP - LAST
C
C     Start of outer cycle for the variable to be dropped.
C
      DO 60 J = FIRST, LAST
        D1 = D(J)
        IF (SQRT(D1) .LT. TOL(J)) THEN
          SS(J) = ZERO
          SMIN = ZERO
          JMIN = J
          GO TO 50
        END IF
        RHS = THETAB(J)
        IF (J .EQ. LAST) GO TO 40
C
C     Copy row J of RBAR into WK.
C
        POS = POS1
        DO 10 I = J+1, LAST
          WK(I) = RBAR(POS)
          POS = POS + 1
   10   CONTINUE
        POS = POS + INC
C
C     Lower the variable past each row.
C
        DO 30 ROW = J+1, LAST
          X = WK(ROW)
          D2 = D(ROW)
          IF (ABS(X) * SQRT(D1) .LT. TOL(ROW) .OR. D2 .EQ. ZERO) THEN
            POS = POS + NP - ROW
            GO TO 30
          END IF
          D1 = D1 * D2 / (D2 + D1 * X**2)
          DO 20 COL = ROW+1, LAST
            WK(COL) = WK(COL) - X * RBAR(POS)
            POS = POS + 1
   20     CONTINUE
          RHS = RHS - X * THETAB(ROW)
          POS = POS + INC
   30   CONTINUE
   40   SS(J) = RHS * D1 * RHS
        IF (SS(J) .LT. SMIN) THEN
          JMIN = J
          SMIN = SS(J)
        END IF
C
C     Update position of first element in row of RBAR.
C
   50   IF (J .LT. LAST) POS1 = POS1 + NP - J
C
   60 CONTINUE
C
      RETURN
      END

      SUBROUTINE EXADD1(IVAR, RSS, BOUND, NVMAX, RESS, IR, NBEST,
     1 LOPT, IL, VORDER, SMAX, JMAX, SS, WK, LAST)
C
C     Update the NBEST subsets of IVAR variables found from a call
C     to subroutine ADD1.
C
      INTEGER IVAR, NVMAX, IR, NBEST, IL, LOPT(IL, NBEST), LAST,
     *   VORDER(LAST), JMAX
      DOUBLE PRECISION RSS(LAST), BOUND(NVMAX), RESS(IR, NBEST), SMAX,
     *   SS(LAST), WK(LAST)
C
C     Local variables
C
      DOUBLE PRECISION ZERO, SSBASE, SM, TEMP
      INTEGER I, J, LTEMP, JM
      DATA ZERO/0.D0/
C
      IF (JMAX .EQ. 0) RETURN
      IF (IVAR .LE. 0) RETURN
      IF (IVAR .GT. NVMAX) RETURN
      LTEMP = VORDER(IVAR)
      JM = JMAX
      SM = SMAX
      IF (IVAR .GT. 1) SSBASE= RSS(IVAR-1)
      IF (IVAR .EQ. 1) SSBASE= RSS(IVAR) + SS(1)
      DO 10 J = IVAR, LAST
      WK(J) = SS(J)
   10 END DO
C
      DO 30 I = 1, NBEST
        TEMP = SSBASE - SM
        IF (TEMP .GE. BOUND(IVAR)) GO TO 40
        VORDER(IVAR) = VORDER(JM)
        IF (JM .EQ. IVAR) VORDER(IVAR) = LTEMP
        CALL REPORT(IVAR, TEMP, BOUND, NVMAX, RESS, IR, NBEST, LOPT, IL,
     *              VORDER)
        IF (I .GE. NBEST) GO TO 40
        WK(JM) = ZERO
        SM = ZERO
        JM = 0
        DO 20 J = IVAR, LAST
          IF (WK(J) .LE. SM) GO TO 20
          JM = J
          SM = WK(J)
   20   CONTINUE
        IF (JM .EQ. 0) GO TO 40
   30 CONTINUE
C
C     Restore VORDER(IVAR)
C
   40 VORDER(IVAR) = LTEMP
C
      RETURN
      END

      SUBROUTINE FORWRD(NP, NRBAR, D, RBAR, THETAB, FIRST, LAST,
     *   VORDER, TOL, RSS, BOUND, NVMAX, RESS, IR, NBEST, LOPT, IL,
     *   WK, IWK, IER)
C
C     Forward selection from variables in positions FIRST .. LAST.
C     If FIRST > 1, variables in positions prior to this are forced in.
C     If LAST < NP, variables in positions after this are forced out.
C     On exit, the array VORDER contains the numbers of the variables
C     in the order in which they were added.
C
      INTEGER NP, NRBAR, FIRST, LAST, VORDER(NP), NVMAX, IR, NBEST,
     *   IL, LOPT(IL, *), IWK, IER
      DOUBLE PRECISION D(NP), RBAR(NRBAR), THETAB(NP), TOL(NP), RSS(NP),
     *   BOUND(NVMAX), RESS(IR, *), WK(IWK)
C
C     Local variables
C
      INTEGER NEED, POS, J1, J2, JMAX
      DOUBLE PRECISION SMAX
C
C     Check call arguments
C
      IER = 0
      IF (FIRST .GE. NP) IER = 1
      IF (LAST .LE. 1) IER = IER + 2
      IF (FIRST .LT. 1) IER = IER + 4
      IF (LAST .GT. NP) IER = IER + 8
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 16
      IF (IWK .LT. 3*LAST) IER = IER + 32
      IF (NBEST .GT. 0) THEN
        NEED = NVMAX*(NVMAX+1)/2
        IF (IR .LT. NVMAX) IER = IER + 64
        IF (IL .LT. NEED) IER = IER + 128
      END IF
      IF (IER .NE. 0) RETURN
C
C     For POS = FIRST .. LAST-1, call ADD1 to find best variable to put
C     into position POS.
C
      J1 = LAST + 1
      J2 = LAST + J1
      DO 10 POS = FIRST, LAST-1
          CALL ADD1(NP, NRBAR, D, RBAR, THETAB, POS, LAST, TOL, WK,
     *            WK(J1), WK(J2), SMAX, JMAX, IER)
          IF (NBEST .GT. 0) CALL EXADD1(POS, RSS, BOUND, NVMAX, RESS,
     *        IR, NBEST, LOPT, IL, VORDER, SMAX, JMAX, WK, WK(J1), LAST)
C
C     Move the best variable to position POS.
C
          IF (JMAX .GT. POS) CALL VMOVE(NP, NRBAR, VORDER, D, RBAR,
     *           THETAB, RSS, JMAX, POS, TOL, IER)
   10 CONTINUE
C
      RETURN
      END

      SUBROUTINE INITR(NP, NVMAX, NBEST, BOUND, RESS, IR, LOPT, IL,
     *          VORDER, RSS, IER)
C
C     Initialize the recording of best subsets
C
      INTEGER NP, NVMAX, NBEST, IR, IL, LOPT(IL,NBEST), VORDER(NP), IER
      DOUBLE PRECISION BOUND(NP), RESS(IR,NBEST), RSS(NP)
C
C     Local variables
C
      INTEGER BEST, POS, NVAR, I
      DOUBLE PRECISION LARGE
      DATA LARGE/1.D+35/
C
C     Check call arguments
C
      IER = 0
      IF (NBEST .LE. 0) IER = 1
      IF (NVMAX .LE. 0) IER = IER + 2
      IF (NVMAX .GT. NP) IER = IER + 4
      IF (IR .LT. NVMAX) IER = IER + 8
      IF (IL .LT. NVMAX*(NVMAX+1)/2) IER = IER + 16
      IF (IER .NE. 0) RETURN
C
C     Initialize arrays BOUND, RESS & LOPT
C
      DO 30 BEST = 1, NBEST
        POS = 1
        DO 20 NVAR = 1, NVMAX
          IF (BEST .EQ. 1) THEN
            RESS(NVAR,BEST) = RSS(NVAR)
          ELSE
            RESS(NVAR,BEST) = LARGE
          END IF
          IF (BEST .EQ. NBEST) BOUND(NVAR) = RESS(NVAR,NBEST)
          DO 10 I = 1, NVAR
            IF (BEST .EQ. 1) THEN
              LOPT(POS,BEST) = VORDER(I)
            ELSE
              LOPT(POS,BEST) = 0
            END IF
            POS = POS + 1
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C
      RETURN
      END

      SUBROUTINE REPORT(POS, SSQ, BOUND, NVMAX, RESS, IR, NBEST, LOPT,
     *   IL, VORDER)
C
C     Update record of the best NBEST subsets of POS variables, if
C     necessary, using SSQ.
C
      INTEGER POS, NVMAX, IR, IL, NBEST, LOPT(IL,NBEST), VORDER(POS)
      DOUBLE PRECISION SSQ, BOUND(NVMAX), RESS(IR,NBEST)
C
C     Local variables
C
      INTEGER RANK, L0, JJ, J, LISTJ, L, I, K
      DOUBLE PRECISION UNDER1, OVER1
      DATA UNDER1/0.9999D0/, OVER1/1.0001D0/
C
C     If residual sum of squares (SSQ) for the new subset > the
C     appropriate bound, return.
C
      IF (POS .GT. NVMAX) RETURN
      IF (SSQ .GE. BOUND(POS)) RETURN
C
C     Find rank of the new subset
C
      DO 30 RANK = 1,NBEST
        IF (SSQ .LE. RESS(POS,RANK)) GO TO 40
   30 CONTINUE
   40 L0 = (POS*(POS-1))/2
C
C     Check that the subset is not a duplicate of one which has already
C     been recorded.
C
      JJ = RANK
      IF (SSQ .GT. UNDER1*RESS(POS,RANK)) GO TO 50
      IF (RANK .EQ. 1) GO TO 90
      IF (SSQ .GT. OVER1*RESS(POS,RANK-1)) GO TO 90
      JJ = RANK-1
   50 DO 70 J = 1, POS
        LISTJ = VORDER(J)
        L = L0
        DO 60 I = 1, POS
          L = L + 1
          IF (LISTJ .EQ. LOPT(L,JJ)) GO TO 70
   60   CONTINUE
        GO TO 80
   70 CONTINUE
      RETURN
   80 JJ = JJ - 1
      IF (JJ .GT. 0 .AND. JJ .EQ. RANK-1) GO TO 50
C
C     Record new subset, and move down the other records.
C
   90 IF (RANK .EQ. NBEST) GO TO 110
      J = NBEST - RANK
      DO 101 I = 1, J
        JJ = NBEST - I
        RESS(POS,JJ+1) = RESS(POS,JJ)
        L = L0
        DO 100 K = 1, POS
          L = L + 1
          LOPT(L,JJ+1) = LOPT(L,JJ)
  100   END DO
  101 END DO
  110 RESS(POS,RANK) = SSQ
      L = L0
      DO 120 K = 1, POS
        L = L + 1
        LOPT(L,RANK) = VORDER(K)
  120 CONTINUE
      BOUND(POS) = RESS(POS,NBEST)
      END

      SUBROUTINE SEQREP(NP, NRBAR, D, RBAR, THETAB, FIRST, LAST,
     *   VORDER, TOL, RSS, BOUND, NVMAX, RESS, IR, NBEST, LOPT, IL,
     *   WK, IWK, IER)
C
C     Sequential replacement algorithm applied to the variables in
C     positions FIRST, ..., LAST.
C     If FIRST > 1, variables in positions prior to this are forced in.
C     If LAST < NP, variables in positions after this are forced out.
C
      INTEGER NP, NRBAR, FIRST, LAST, VORDER(NP), NVMAX, IR, NBEST,
     *   IL, LOPT(IL, *), IWK, IER
      DOUBLE PRECISION D(NP), RBAR(NRBAR), THETAB(NP), TOL(NP), RSS(NP),
     *   BOUND(NVMAX), RESS(IR, *), WK(IWK)
C
C     Local variables
C
      INTEGER NEED, J1, J2, NV, SIZE, START, BEST, FROM, I, JMAX, COUNT
      DOUBLE PRECISION ZERO, SSRED, SMAX
      DATA ZERO/0.D0/
C
C     Check call arguments
C
      IER = 0
      IF (FIRST .GE. NP) IER = 1
      IF (LAST .LE. 1) IER = IER + 2
      IF (FIRST .LT. 1) IER = IER + 4
      IF (LAST .GT. NP) IER = IER + 8
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 16
      IF (IWK .LT. 3*LAST) IER = IER + 32
      IF (NBEST .GT. 0) THEN
        NEED = NVMAX*(NVMAX+1)/2
        IF (IR .LT. NVMAX) IER = IER + 64
        IF (IL .LT. NEED) IER = IER + 128
      END IF
      IF (IER .NE. 0 .OR. NBEST .LE. 0) RETURN
C
      J1 = 1 + LAST
      J2 = J1 + LAST
      NV = MIN(NVMAX, LAST-1)
C
C     Outer loop; SIZE = current size of subset being considered.
C
      DO 30 SIZE = FIRST, NV
        COUNT = 0
        START = FIRST
   10   SSRED = ZERO
        BEST = 0
        FROM = 0
C
C     Find the best variable from those in positions SIZE+1, ..., LAST
C     to replace the one in position SIZE.   Then rotate variables in
C     positions START, ..., SIZE.
C
        DO 20 I = START, SIZE
          CALL ADD1(NP, NRBAR, D, RBAR, THETAB, SIZE, LAST, TOL, WK,
     *                    WK(J1), WK(J2), SMAX, JMAX, IER)
          IF (JMAX .GT. SIZE) THEN
            CALL EXADD1(SIZE, RSS, BOUND, NVMAX, RESS, IR, NBEST,
     *                  LOPT, IL, VORDER, SMAX, JMAX, WK, WK(J1), LAST)
            IF (SMAX .GT. SSRED) THEN
              SSRED = SMAX
              BEST = JMAX
              IF (I .LT. SIZE) THEN
                FROM = SIZE + START - I - 1
              ELSE
                FROM = SIZE
              END IF
            END IF
          END IF
          IF (I .LT. SIZE) CALL VMOVE(NP, NRBAR, VORDER, D, RBAR,
     *                  THETAB, RSS, SIZE, START, TOL, IER)
   20   CONTINUE
C
C     If any replacement reduces the RSS, make the best one.
C     Move variable from position FROM to SIZE.
C     Move variable from position BEST to FIRST.
C
        IF (BEST .GT. SIZE) THEN
          IF (FROM .LT. SIZE) CALL VMOVE(NP, NRBAR, VORDER, D, RBAR,
     *                  THETAB, RSS, FROM, SIZE, TOL, IER)
          CALL VMOVE(NP, NRBAR, VORDER, D, RBAR, THETAB, RSS, BEST,
     *                  FIRST, TOL, IER)
          COUNT = 0
          START = FIRST + 1
        ELSE
          COUNT = COUNT + 1
        END IF
C
C     Repeat until COUNT = SIZE - START + 1
C
        IF (COUNT .LE. SIZE - START) GO TO 10
   30 CONTINUE
C
      RETURN
      END

      SUBROUTINE XHAUST(NP, NRBAR, D, RBAR, THETAB, FIRST, LAST,
     *   VORDER, TOL, RSS, BOUND, NVMAX, RESS, IR, NBEST, LOPT, IL,
     *   WK, DIMWK, IWK, DIMIWK, IER)
C
C     Exhaustive search algorithm, using leaps and bounds, applied to
C     the variables in positions FIRST, ..., LAST.
C     If FIRST > 1, variables in positions prior to this are forced in.
C     If LAST < NP, variables in positions after this are forced out.
C
      INTEGER NP, NRBAR, FIRST, LAST, VORDER(NP), NVMAX, IR, NBEST,
     *   IL, LOPT(IL, *), DIMWK, DIMIWK, IWK(DIMIWK), IER
      DOUBLE PRECISION D(NP), RBAR(NRBAR), THETAB(NP), TOL(NP), RSS(NP),
     *   BOUND(NVMAX), RESS(IR, *), WK(DIMWK)
C
C     Local variables
C
      INTEGER NEED, J1, J2, ROW, I, JMAX, IPT, NEWPOS
      DOUBLE PRECISION SMAX, TEMP
C
C     Check call arguments
C
      IER = 0
      IF (FIRST .GE. NP) IER = 1
      IF (LAST .LE. 1) IER = IER + 2
      IF (FIRST .LT. 1) IER = IER + 4
      IF (LAST .GT. NP) IER = IER + 8
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 16
      IF (DIMWK .LT. 3*LAST .OR. DIMIWK .LT. NVMAX) IER = IER + 32
      IF (NBEST .GT. 0) THEN
        NEED = NVMAX*(NVMAX+1)/2
        IF (IR .LT. NVMAX) IER = IER + 64
        IF (IL .LT. NEED) IER = IER + 128
      END IF
      IF (IER .NE. 0 .OR. NBEST .LE. 0) RETURN
C
      J1 = 1 + LAST
      J2 = J1 + LAST
C
C     Record subsets contained in the initial ordering, including check
C     for variables which are linearly related to earlier variables.
C     This should be redundant if the user has first called SING and
C     INITR.
C
      DO 10 ROW = FIRST, NVMAX
        IF (D(ROW) .LE. TOL(ROW)) THEN
          IER = -999
          RETURN
        END IF
        CALL REPORT(ROW, RSS(ROW), BOUND, NVMAX, RESS, IR, NBEST, LOPT,
     *        IL, VORDER)
   10 CONTINUE
C
C     IWK(I) contains the upper limit for the I-th simulated DO-loop for
C     I = FIRST, ..., NVMAX-1.
C     IPT points to the current DO loop.
C
      DO 20 I = FIRST, NVMAX
      IWK(I) = LAST
   20 END DO
C
C     Innermost loop.
C     Find best possible variable for position NVMAX from those in
C     positions NVMAX, .., IWK(NVMAX).
C
   30 CALL ADD1(NP, NRBAR, D, RBAR, THETAB, NVMAX, IWK(NVMAX), TOL, WK,
     *          WK(J1), WK(J2), SMAX, JMAX, IER)
      CALL EXADD1(NVMAX, RSS, BOUND, NVMAX, RESS, IR, NBEST, LOPT, IL,
     *            VORDER, SMAX, JMAX, WK, WK(J1), IWK(NVMAX))
C
C     Move to next lower numbered loop which has not been exhausted.
C
      IPT = NVMAX - 1
   40 IF (IPT .GE. IWK(IPT)) THEN
        IPT = IPT - 1
        IF (IPT .GE. FIRST) GO TO 40
        RETURN
      END IF
C
C     Lower variable from position IPT to position IWK(IPT).
C     Record any good new subsets found by the move.
C
      NEWPOS = IWK(IPT)
      CALL VMOVE(NP, NRBAR, VORDER, D, RBAR, THETAB, RSS, IPT, NEWPOS,
     *           TOL, IER)
      DO 50 I = IPT, MIN(NVMAX, NEWPOS-1)
      CALL REPORT(I, RSS(I), BOUND, NVMAX, RESS, IR, NBEST, LOPT, IL,
     *             VORDER)
   50 END DO
C
C     Reset all ends of loops for I >= IPT.
C
      DO 60 I = IPT, NVMAX
      IWK(I) = NEWPOS - 1
   60 END DO
C
C     If residual sum of squares for all variables above position NEWPOS
C     is greater than BOUND(I), no better subsets of size I can be found
C     inside the current loop.
C
      TEMP = RSS(NEWPOS-1)
      DO 70 I = IPT, NVMAX
        IF (TEMP .GT. BOUND(I)) GO TO 80
   70 CONTINUE
      IF (IWK(NVMAX) .GT. NVMAX) GO TO 30
      IPT = NVMAX - 1
      GO TO 40
   80 IPT = I - 1
      IF (IPT .LT. FIRST) RETURN
      GO TO 40
C
      END
C
C___      SUBROUTINE EFROYM(NP, NRBAR, D, RBAR, THETAB, FIRST, LAST,
C___     *   FIN, FOUT, SIZE, NOBS, VORDER, TOL, RSS, BOUND, NVMAX, RESS,
C___     *   IR, NBEST, LOPT, IL, WK, IWK, IER)
C___C
C___C     Efroymson's stepwise regression from variables in positions FIRST,
C___C     ..., LAST.  If FIRST > 1, variables in positions prior to this are
C___C     forced in.  If LAST < NP, variables in positions after this are
C___C     forced out.
C___C
C___c     IMPLICIT NONE
C___      INTEGER NP, NRBAR, FIRST, LAST, SIZE, NOBS, VORDER(NP), NVMAX, IR,
C___     *   NBEST, IL, LOPT(IL, *), IWK, IER
C___      DOUBLE PRECISION D(NP), RBAR(NRBAR), THETAB(NP), FIN, FOUT,
C___     *   TOL(NP), RSS(NP), BOUND(NVMAX), RESS(IR, *), WK(IWK)
C___C
C___C     Local variables
C___C
C___      INTEGER NEED, J1, J2, JMAX, JMIN, I
C___      DOUBLE PRECISION ONE, EPS, ZERO, SMAX, BASE, VAR, F, SMIN
C___      DATA ONE/1.D0/, EPS/1.D-16/, ZERO/0.D0/
C___C
C___C     Check call arguments
C___C
C___      IER = 0
C___      IF (FIRST .GE. NP) IER = 1
C___      IF (LAST .LE. 1) IER = IER + 2
C___      IF (FIRST .LT. 1) IER = IER + 4
C___      IF (LAST .GT. NP) IER = IER + 8
C___      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 16
C___      IF (IWK .LT. 3*LAST) IER = IER + 32
C___      IF (NBEST .GT. 0) THEN
C___               NEED = NVMAX*(NVMAX+1)/2
C___          IF (IR .LT. NVMAX) IER = IER + 64
C___          IF (IL .LT. NEED) IER = IER + 128
C___      END IF
C___      IF (FIN .LT. FOUT .OR. FIN .LE. ZERO) IER = IER + 256
C___      IF (NOBS .LE. NP) IER = IER + 512
C___      IF (IER .NE. 0) RETURN
C___C
C___C     EPS approximates the smallest quantity such that the calculated
C___C     value of (1 + EPS) is > 1.   It is used to test for a perfect fit
C___C     (RSS = 0).
C___C
C___   10 IF (ONE + EPS .LE. ONE) THEN
C___          EPS = EPS + EPS
C___          GO TO 10
C___      END IF
C___C
C___C     SIZE = number of variables in the current subset
C___C
C___      SIZE = FIRST - 1
C___      J1 = LAST + 1
C___      J2 = LAST + J1
C___C
C___C     Find the best variable to add next
C___C
C___   20 CALL ADD1(NP, NRBAR, D, RBAR, THETAB, SIZE+1, LAST, TOL, WK,
C___     *   WK(J1), WK(J2), SMAX, JMAX, IER)
C___      IF (NBEST .GT. 0) CALL EXADD1(SIZE+1, RSS, BOUND, NVMAX, RESS,
C___     *    IR, NBEST, LOPT, IL, VORDER, SMAX, JMAX, WK, WK(J1), LAST)
C___      write(*, *) 'Best variable to add: ', VORDER(JMAX)
C___C
C___C     Calculate 'F-to-enter' value
C___C
C___      IF (SIZE .GT. 0) THEN
C___          BASE = RSS(SIZE)
C___      ELSE
C___          BASE = RSS(1) + WK(1)
C___      END IF
C___      VAR = (BASE - SMAX) / (NOBS - SIZE - 1)
C___      IF (VAR .LT. EPS*BASE) THEN
C___          IER = -1
C___          F = ZERO
C___      ELSE
C___          F = SMAX / VAR
C___      END IF
C___      write(*, 900) F
C___  900 format(' F-to-enter = ', f10.2)
C___C
C___C     Exit if F < FIN or IER < 0 (perfect fit)
C___C
C___      IF (F .LT. FIN .OR. IER .LT. 0) RETURN
C___C
C___C     Add the variable to the subset (in position FIRST).
C___C
C___      SIZE = SIZE + 1
C___      IF (JMAX .GT. FIRST) CALL VMOVE(NP, NRBAR, VORDER, D, RBAR,
C___     *    THETAB, RSS, JMAX, FIRST, TOL, IER)
C___C
C___C     See whether a variable entered earlier can be deleted now.
C___C
C___   30 IF (SIZE .LE. FIRST) GO TO 20
C___      CALL DROP1(NP, NRBAR, D, RBAR, THETAB, FIRST+1, SIZE, TOL, WK,
C___     *    WK(J1), SMIN, JMIN, IER)
C___      VAR = RSS(SIZE) / (NOBS - SIZE)
C___      F = SMIN / VAR
C___      write(*, 910) VORDER(JMIN), F
C___  910 format(' F-to-drop variable: ', i4, ' = ', f10.2)
C___      IF (F .LT. FOUT) THEN
C___        CALL VMOVE(NP, NRBAR, VORDER, D, RBAR, THETAB, RSS, JMIN, SIZE,
C___     *        TOL, IER)
C___        IF (NBEST .GT. 0) THEN
C___          DO 40 I = JMIN, SIZE-1
C___   40     CALL REPORT(I, RSS(I), BOUND, NVMAX, RESS, IR, NBEST, LOPT,
C___     *        IL, VORDER)
C___        END IF
C___        SIZE = SIZE - 1
C___        GO TO 30
C___      END IF
C___C
C___      GO TO 20
C___      END

      SUBROUTINE REGCF(NP, NRBAR, D, RBAR, THETAB, TOL, BETA,
     +     NREQ, IER)
C
C     ALGORITHM AS274  APPL. STATIST. (1992) VOL 41, NO. x
C
C     Modified version of AS75.4 to calculate regression coefficients
C     for the first NREQ variables, given an orthogonal reduction from
C     AS75.1.
C
      INTEGER NP, NRBAR, NREQ, IER
      DOUBLE PRECISION D(NP), RBAR(*), THETAB(NP), TOL(NP),
     +     BETA(NP)
C
C     Local variables
C
      INTEGER I, J, NEXTR
      DOUBLE PRECISION ZERO
C
      DATA ZERO/0.D0/
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (NREQ .LT. 1 .OR. NREQ .GT. NP) IER = IER + 4
      IF (IER .NE. 0) RETURN
C
      DO 20 I = NREQ, 1, -1
      IF (SQRT(D(I)) .LT. TOL(I)) THEN
      BETA(I) = ZERO
      D(I) = ZERO
      GO TO 20
      END IF
      BETA(I) = THETAB(I)
      NEXTR = (I-1) * (NP+NP-I)/2 + 1
      DO 10 J = I+1, NREQ
      BETA(I) = BETA(I) - RBAR(NEXTR) * BETA(J)
      NEXTR = NEXTR + 1
   10 CONTINUE
   20 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE SING(NP, NRBAR, D, RBAR, THETAB, SSERR, TOL,
     +   LINDEP, WORK, IER)
C
C     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
C
C     Checks for singularities, reports, and adjusts orthogonal
C     reductions produced by AS75.1.
C
      INTEGER NP, NRBAR, IER
      DOUBLE PRECISION D(NP), RBAR(NRBAR), THETAB(NP), SSERR,
     +      TOL(NP), WORK(NP)
      LOGICAL LINDEP(NP)
C
C     Local variables
C
      DOUBLE PRECISION ZERO, TEMP
      INTEGER COL, POS, ROW, NC2, POS2
C
      DATA ZERO/0.D0/
C
C     Check input parameters
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (IER .NE. 0) RETURN
C
      DO 10 COL = 1, NP
      WORK(COL) = SQRT(D(COL))
   10 END DO
C
      DO 40 COL = 1, NP
C
C     Set elements within RBAR to zero if they are less than TOL(COL) in
C     absolute value after being scaled by the square root of their row
C     multiplier.
C
      TEMP = TOL(COL)
      POS = COL - 1
      DO 30 ROW = 1, COL-1
      IF (ABS(RBAR(POS)) * WORK(ROW) .LT. TEMP) RBAR(POS) = ZERO
      POS = POS + NP - ROW - 1
   30 CONTINUE
C
C     If diagonal element is near zero, set it to zero, set appropriate
C     element of LINDEP, and use INCLUD to augment the projections in
C     the lower rows of the orthogonalization.
C
      LINDEP(COL) = .FALSE.
      IF (WORK(COL) .LE. TEMP) THEN
      LINDEP(COL) = .TRUE.
      IER = IER - 1
      IF (COL .LT. NP) THEN
        NC2 = NP - COL
        POS2 = POS + NP - COL + 1
        CALL INCLUD(NC2, NC2*(NC2-1)/2, D(COL), RBAR(POS+1),
     +            THETAB(COL), D(COL+1), RBAR(POS2), THETAB(COL+1),
     +            SSERR, IER)
      ELSE
        SSERR = SSERR + D(COL) * THETAB(COL)**2
      END IF
      D(COL) = ZERO
      WORK(COL) = ZERO
      THETAB(COL) = ZERO
      END IF
   40 CONTINUE
      RETURN
      END
C
      SUBROUTINE SSLEAPS(NP, D, THETAB, SSERR, RSS, IER)
C
C     ALGORITHM AS274  APPL. STATIST. (1992) VOL. 41, NO. 2
C
C     Calculates partial residual sums of squares from an orthogonal
C     reduction from AS75.1.
C
      INTEGER NP, IER
      DOUBLE PRECISION D(NP), THETAB(NP), SSERR, RSS(NP)
C
C     Local variables
C
      INTEGER I
      DOUBLE PRECISION SUM
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (IER .NE. 0) RETURN
C
      SUM = SSERR
      RSS(NP) = SSERR
      DO 10 I = NP, 2, -1
      SUM = SUM + D(I) * THETAB(I)**2
      RSS(I-1) = SUM
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE TOLSET(NP, NRBAR, D, RBAR, TOL, WORK, IER)
C
C     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
C
C     Sets up array TOL for testing for zeros in an orthogonal
C     reduction formed using AS75.1.
C
      INTEGER NP, NRBAR, IER
      DOUBLE PRECISION D(NP), RBAR(*), TOL(NP), WORK(NP)
C
C     Local variables.
C
      INTEGER COL, ROW, POS
      DOUBLE PRECISION EPS, SUM, ZERO
C
C     EPS is a machine-dependent constant.   For compilers which use
C     the IEEE format for floating-point numbers, recommended values
C     are 1.E-06 for single precision and 1.D-12 for double precision.
C
c     changed EPS from 10^-12 to 5x10^-10 to try to fix a bug
      DATA EPS/5.D-10/, ZERO/0.D0/
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (IER .NE. 0) RETURN
C
C     Set TOL(I) = sum of absolute values in column I of RBAR after
C     scaling each element by the square root of its row multiplier.
C
      DO 10 ROW = 1, NP
      WORK(ROW) = SQRT(D(ROW))
   10 END DO
      DO 30 COL = 1, NP
      POS = COL - 1
      IF (COL .LE. NP) THEN
      SUM = WORK(COL)
      ELSE
      SUM = ZERO
      END IF
      DO 20 ROW = 1, MIN(COL-1, NP)
      SUM = SUM + ABS(RBAR(POS)) * WORK(ROW)
      POS = POS + NP - ROW - 1
  20  CONTINUE
      TOL(COL) = EPS * SUM
  30  CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PCORR(NP, NRBAR, D, RBAR, THETAB, SSERR, IN,
     +      WORK, CORMAT, DIMC, YCORR, IER)
C
C     ALGORITHM AS274  APPL. STATIST. (1992) VOL. 41, NO. 2
C
C     Calculate partial correlations after the first IN variables
C     have been forced into the regression.
C
C     Auxiliary routine called: COR
C
      INTEGER NP, NRBAR, IN, DIMC, IER
      DOUBLE PRECISION D(NP), RBAR(*), THETAB(NP), SSERR,
     +     WORK(NP), CORMAT(*), YCORR(NP)
C
C     Local variables.
C
      INTEGER START, IN1, I
      DOUBLE PRECISION ZERO
C
      DATA ZERO/0.D0/
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (IN .LT. 0 .OR. IN .GT. NP-1) IER = IER + 4
      IF (DIMC .LT. (NP-IN)*(NP-IN-1)/2) IER = IER + 8
      IF (IER .NE. 0) RETURN
C
      START = IN * (NP+NP-IN-1)/2 + 1
      IN1 = IN + 1
      CALL COR(NP-IN, D(IN1), RBAR(START), THETAB(IN1),
     +  SSERR, WORK, CORMAT, YCORR)
C
C     Check for zeros.
C
      DO 10 I = 1, NP-IN
      IF (WORK(I) .LE. ZERO) IER = -I
   10 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE COR(NP, D, RBAR, THETAB, SSERR, WORK, CORMAT,
     +    YCORR)
C
C     ALGORITHM AS274  APPL. STATIST. (1992) VOL. 41, NO. 2
C
C     Calculate correlations from an orthogonal reduction.   This
C     routine will usually be called from PCORR, which will have
C     removed the appropriate number of rows at the start.
C
      INTEGER NP
      DOUBLE PRECISION D(NP), RBAR(*), THETAB(NP), SSERR,
     +     WORK(NP), CORMAT(*), YCORR(NP)
C
C     Local variables.
C
      INTEGER ROW, POS, COL1, POS1, COL2, POS2, DIFF
      DOUBLE PRECISION SUMY, SUM, ZERO
C
      DATA ZERO/0.D0/
C
C     Process by columns, including the projections of the dependent
C     variable (THETAB).
C
      SUMY = SSERR
      DO 10 ROW = 1, NP
      SUMY = SUMY + D(ROW) * THETAB(ROW)**2
   10 END DO
      SUMY = SQRT(SUMY)
      POS = NP*(NP-1)/2
      DO 70 COL1 = NP, 1, -1
C
C     Calculate the length of column COL1.
C
      SUM = D(COL1)
      POS1 = COL1 - 1
      DO 20 ROW = 1, MIN(COL1-1, NP)
      SUM = SUM + D(ROW) * RBAR(POS1)**2
      POS1 = POS1 + NP - ROW - 1
  20  CONTINUE
      WORK(COL1) = SQRT(SUM)
C
C     If SUM = 0, set all correlations with this variable to zero.
C
      IF (SUM .EQ. ZERO) THEN
      YCORR(COL1) = ZERO
      DO 30 COL2 = NP, COL1+1, -1
      CORMAT(POS) = ZERO
      POS = POS - 1
   30 CONTINUE
      GO TO 70
      END IF
C
C     Form cross-products, then divide by product of column lengths.
C
      SUM = D(COL1) * THETAB(COL1)
      POS1 = COL1 - 1
      DO 40 ROW = 1, MIN(COL1-1, NP)
      SUM = SUM + D(ROW) * RBAR(POS1) * THETAB(ROW)
      POS1 = POS1 + NP - ROW - 1
   40 CONTINUE
      YCORR(COL1) = SUM / (SUMY * WORK(COL1))
C
      DO 60 COL2 = NP, COL1+1, -1
      IF (WORK(COL2) .GT. ZERO) THEN
        POS1 = COL1 - 1
        POS2 = COL2 - 1
        DIFF = COL2 - COL1
        SUM = ZERO
        DO 50 ROW = 1, MIN(COL1-1, NP)
          SUM = SUM + D(ROW) * RBAR(POS1) * RBAR(POS2)
          POS1 = POS1 + NP - ROW - 1
          POS2 = POS1 + DIFF
   50   CONTINUE
        SUM = SUM + D(COL1) * RBAR(POS2)
        CORMAT(POS) = SUM / (WORK(COL1) * WORK(COL2))
        ELSE
        CORMAT(POS) = ZERO
        END IF
        POS = POS - 1
   60   CONTINUE
   70   CONTINUE
C
      RETURN
      END
C
      SUBROUTINE VMOVE(NP, NRBAR, VORDER, D, RBAR, THETAB,
     +    RSS, FROM, TO, TOL, IER)
C
C     ALGORITHM AS274 APPL. STATIST. (1992) VOL.41, NO. 2
C
C     Move variable from position FROM to position TO in an
C     orthogonal reduction produced by AS75.1.
C
      INTEGER NP, NRBAR, VORDER(NP), FROM, TO, IER
      DOUBLE PRECISION D(NP), RBAR(*), THETAB(NP), RSS(NP),
     +    TOL(NP)
C
C     Local variables
C
      DOUBLE PRECISION ZERO, D1, D2, X, ONE, D1NEW, D2NEW, CBAR, SBAR, Y
      INTEGER M, FIRST, LAST, INC, M1, M2, MP1, COL, POS, ROW
C
      DATA ZERO/0.D0/, ONE/1.D0/
C
C     Check input parameters
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (FROM .LT. 1 .OR. FROM .GT. NP) IER = IER + 4
      IF (TO .LT. 1 .OR. TO .GT. NP) IER = IER + 8
      IF (IER .NE. 0) RETURN
C
      IF (FROM .EQ. TO) RETURN
C
      IF (FROM .LT. TO) THEN
      FIRST = FROM
      LAST = TO - 1
       INC = 1
      ELSE
      FIRST = FROM - 1
      LAST = TO
      INC = -1
      END IF
      DO 70 M = FIRST, LAST, INC
C
C     Find addresses of first elements of RBAR in rows M and (M+1).
C
      M1 = (M-1)*(NP+NP-M)/2 + 1
      M2 = M1 + NP - M
      MP1 = M + 1
      IF (M .LE. NP) THEN
      D1 = D(M)
      IF (MP1 .LE. NP) THEN
        D2 = D(MP1)
          ELSE
        D2 = ZERO
      END IF
      ELSE
      D1 = ZERO
      D2 = ZERO
      END IF
C
C     Special cases.
C
      IF (D1 .EQ. ZERO .AND. D2 .EQ. ZERO) GO TO 40
       X = RBAR(M1)
      IF (ABS(X) * SQRT(D1) .LT. TOL(MP1)) THEN
       X = ZERO
      END IF
      IF (D1 .EQ. ZERO .OR. X .EQ. ZERO) THEN
      D(M) = D2
      D(MP1) = D1
      RBAR(M1) = ZERO
      DO 10 COL = M+2, NP
        M1 = M1 + 1
        X = RBAR(M1)
        RBAR(M1) = RBAR(M2)
        RBAR(M2) = X
        M2 = M2 + 1
   10     CONTINUE
       X = THETAB(M)
       THETAB(M) = THETAB(MP1)
       THETAB(MP1) = X
       GO TO 40
      ELSE IF (D2 .EQ. ZERO) THEN
       D(M) = D1 * X**2
       RBAR(M1) = ONE / X
       DO 20 COL = M+2, NP
        M1 = M1 + 1
        RBAR(M1) = RBAR(M1) / X
   20   CONTINUE
       THETAB(M) = THETAB(M) / X
       GO TO 40
       END IF
C
C     Planar rotation in regular case.
C
      D1NEW = D2 + D1*X**2
      CBAR = D2 / D1NEW
      SBAR = X * D1 / D1NEW
      D2NEW = D1 * CBAR
      D(M) = D1NEW
      D(MP1) = D2NEW
      RBAR(M1) = SBAR
      DO 30 COL = M+2, NP
        M1 = M1 + 1
        Y = RBAR(M1)
        RBAR(M1) = CBAR*RBAR(M2) + SBAR*Y
        RBAR(M2) = Y - X*RBAR(M2)
        M2 = M2 + 1
   30  CONTINUE
       Y = THETAB(M)
       THETAB(M) = CBAR*THETAB(MP1) + SBAR*Y
       THETAB(MP1) = Y - X*THETAB(MP1)
C
C     Swap columns M and (M+1) down to row (M-1).
C
   40 IF (M .EQ. 1) GO TO 60
       POS = M
      DO 50 ROW = 1, M-1
      X = RBAR(POS)
      RBAR(POS) = RBAR(POS-1)
      RBAR(POS-1) = X
      POS = POS + NP - ROW - 1
   50   CONTINUE
C
C     Adjust variable order (VORDER), the tolerances (TOL) and
C     the vector of residual sums of squares (RSS).
C
   60   M1 = VORDER(M)
      VORDER(M) = VORDER(MP1)
      VORDER(MP1) = M1
      X = TOL(M)
       TOL(M) = TOL(MP1)
       TOL(MP1) = X
      RSS(M) = RSS(MP1) + D(MP1) * THETAB(MP1)**2
   70 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE REORDR(NP, NRBAR, VORDER, D, RBAR, THETAB,
     +    RSS, TOL, LIST, N, POS1, IER)
C
C     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
C
C     Re-order the variables in an orthogonal reduction produced by
C     AS75.1 so that the N variables in LIST start at position POS1,
C     though will not necessarily be in the same order as in LIST.
C     Any variables in VORDER before position POS1 are not moved.
C
C     Auxiliary routine called: VMOVE
C
      INTEGER NP, NRBAR, VORDER(NP), N, LIST(N), POS1, IER
      DOUBLE PRECISION D(NP), RBAR(NRBAR), THETAB(NP), RSS(NP),
     +    TOL(NP)
C
C     Local variables.
C
      INTEGER NEXT, I, L, J
C
C     Check N.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (N .LT. 1 .OR. N .GT. NP+1-POS1) IER = IER + 4
      IF (IER .NE. 0) RETURN
C
C     Work through VORDER finding variables which are in LIST.
C
      NEXT = POS1
      I = POS1
   10 L = VORDER(I)
      DO 20 J = 1, N
      IF (L .EQ. LIST(J)) GO TO 40
   20 CONTINUE
   30 I = I + 1
      IF (I .LE. NP) GO TO 10
C
C     If this point is reached, one or more variables in LIST has not
C     been found.
C
      IER = NEXT - N - 1
      RETURN
C
C     Variable L is in LIST; move it up to position NEXT if it is not
C     already there.
C
   40 IF (I .GT. NEXT) CALL VMOVE(NP, NRBAR, VORDER, D, RBAR,
     +   THETAB, RSS, I, NEXT, TOL, IER)
      NEXT = NEXT + 1
      IF (NEXT .LT. N+POS1) GO TO 30
C
      RETURN
      END


      SUBROUTINE INCLUD(NP, NRBAR, WEIGHT, XROW, YELEM, D,
     +      RBAR, THETAB, SSERR, IER)
C
C     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
C     Modified from algorithm AS 75.1
C
C     Calling this routine updates d, rbar, thetab and sserr by the
C     inclusion of xrow, yelem with the specified weight.   The number
C     of columns (variables) may exceed the number of rows (cases).
C
C**** WARNING: The elements of XROW are overwritten  ****
C
      INTEGER NP, NRBAR, IER
      DOUBLE PRECISION WEIGHT, XROW(NP), YELEM, D(NP), RBAR(*),
     +    THETAB(NP), SSERR
C
C     Local variables
C
      INTEGER I, K, NEXTR
      DOUBLE PRECISION ZERO, W, Y, XI, DI, WXI, DPI, CBAR, SBAR, XK
C
      DATA ZERO/0.D0/
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (IER .NE. 0) RETURN
C
      W = WEIGHT
      Y = YELEM
      NEXTR = 1
      DO 30 I = 1, NP
C
C     Skip unnecessary transformations.   Test on exact zeros must be
C     used or stability can be destroyed.
C
      IF (W .EQ. ZERO) RETURN
       XI = XROW(I)
      IF (XI .EQ. ZERO) THEN
       NEXTR = NEXTR + NP - I
      GO TO 30
      END IF
      DI = D(I)
      WXI = W * XI
      DPI = DI + WXI*XI
      CBAR = DI / DPI
      SBAR = WXI / DPI
      W = CBAR * W
      D(I) = DPI
      IF (I .EQ. NP) GO TO 20
      DO 10 K = I+1, NP
        XK = XROW(K)
        XROW(K) = XK - XI * RBAR(NEXTR)
        RBAR(NEXTR) = CBAR * RBAR(NEXTR) + SBAR * XK
        NEXTR = NEXTR + 1
   10   CONTINUE
   20   XK = Y
       Y = XK - XI * THETAB(I)
       THETAB(I) = CBAR * THETAB(I) + SBAR * XK
   30  CONTINUE
C
C     Y * SQRT(W) is now equal to Brown & Durbin's recursive residual.
C
      SSERR = SSERR + W * Y * Y
C
      RETURN
      END
c___C      PROGRAM SUBSET
c___C
c___C     Interactive program to perform regressions on subsets of
c___C     variables.   Max. no. of variables, excl. constant = 50.
c___C
c___C     Subroutines called:-
c___C     ADD1, BAKWRD, EFROYM, XHAUST, FORWRD, EXADD1, INITR, REGCF,
c___C     REORDR, LSORT, PCORR, REPORT, SEQREP, SHELL, SS, TOLSET.
c___C
c___C     Latest revision - 10 November 1993
c___C
c___C     IMPLICIT NONE
c___      INTEGER MAXCOL, MAXSUB, MAXBST, MAXL, MAXR
c___      PARAMETER (MAXCOL=50, MAXSUB=25, MAXBST=20, MAXL=1000,
c___     +           MAXR=MAXSUB*MAXBST)
c___      INTEGER UDIM, IIW, IW
c___      PARAMETER (UDIM=MAXCOL*(MAXCOL+1)/2, IIW=3*MAXCOL, IW=UDIM+IIW)
c___      CHARACTER FNAME*30, ANS, OPTION(22), VNAME(0:MAXCOL)*8, YNAME*8
c___      LOGICAL LSEL, OK
c___      INTEGER LIN, LOUT, LPR, LOPT(MAXL), IWK(IIW), K, IRTN, LINE,
c___     +        ICONST, NCOLS, I1, IPOS, I, NOBS, VORDER(0:MAXCOL), IL,
c___     +        NRBAR, NVMAX, NVMX, NBEST, IOPT, IER, NDF, J, NV, NB,
c___     +        IPRINT, IR, L, IPROC, FIRST, LAST, SIZE, M, ILNB
c___      DOUBLE PRECISION U(UDIM), EL(0:MAXCOL), RHS(0:MAXCOL), RESSQ,
c___     +        SSQ(0:MAXCOL), TOL(0:MAXCOL), BOUND(MAXSUB), RESS(MAXR),
c___     +        WK(IW), TEMP, FIN, FOUT
c___      REAL VAR
c___      DATA OPTION/'C', 'c', 'F', 'f', 'B', 'b', 'R', 'r', 'E', 'e',
c___     +    'P', 'p', 'I', 'i', 'O', 'o', 'L', 'l', 'X', 'x', 'Q', 'q'/
c___C
c___C     Set unit numbers for I/O in LIN & LOUT below.
c___C
c___      DATA LIN/5/, LOUT/6/
c___C
c___C     Ask for name of the data set.
c___C
c___   10 WRITE(LOUT,9000)
c___ 9000 FORMAT(' Enter name of .RED file for data (e.g. B:myfile): ')
c___      READ(LIN,8000) FNAME
c___ 8000 FORMAT(A)
c___C
c___C     Add the .RED extension if necessary.
c___C
c___      IF (INDEX(FNAME, '.RED') .EQ. 0) THEN
c___    IPOS = INDEX(FNAME, ' ')
c___    IF (IPOS .EQ. 0 .OR. IPOS .GT. 11) THEN
c___      WRITE(LOUT, 9010) FNAME
c___ 9010     FORMAT(' ** Illegal filename entered - ', A, ' **')
c___      GO TO 10
c___    END IF
c___    FNAME(IPOS: IPOS+3) = '.RED'
c___      END IF
c___C
c___C     Check that file exists.
c___C
c___      INQUIRE(FILE=FNAME, EXIST=OK)
c___      IF (.NOT. OK) THEN
c___    WRITE(LOUT, 9020) FNAME
c___ 9020   FORMAT(' ** File not found - ', A, ' **')
c___    GO TO 10
c___      END IF
c___      OPEN(9, FILE=FNAME, STATUS='OLD', ACCESS='SEQUENTIAL',
c___     +        FORM='UNFORMATTED')
c___C
c___C     Read contents of file.
c___C
c___      READ(9) K, ICONST, NCOLS, NOBS, NRBAR, LSEL
c___      IF (ICONST .EQ. 0) THEN
c___    READ(9) YNAME, (VNAME(I),I=1,K)
c___    READ(9) (U(I),I=1,NRBAR), (EL(I),I=1,K), (RHS(I),I=1,K), RESSQ
c___      ELSE
c___    READ(9) YNAME, (VNAME(I),I=0,K)
c___    READ(9) (U(I),I=1,NRBAR), (EL(I),I=0,K), (RHS(I),I=0,K), RESSQ
c___      END IF
c___      I1 = 1 + ICONST
c___      WRITE(LOUT, 9030) K, NOBS, YNAME
c___ 9030 FORMAT(' No. of predictors = ', I3, 5X, 'No. of cases = ', I5/
c___     +       ' Dependant variable is ', A)
c___      WRITE(LOUT, 9930) (I, VNAME(I),I=1,K)
c___C
c___C     Initially, all variables except the constant (if there is one)
c___C     are considered candidates for either inclusion or exclusion
c___C
c___      FIRST = I1
c___      LAST = NCOLS
c___C
c___C     Set up array VORDER.
c___C
c___      DO 30 I = 0, K
c___    VORDER(I) = I
c___   30 CONTINUE
c___C
c___C     Ask for values of NVMAX & NBEST.
c___C
c___   50 WRITE(LOUT, 9040)
c___ 9040 FORMAT(' Enter max. size of subsets (excl. constant): ')
c___      READ(LIN, 8010) NVMAX
c___ 8010 FORMAT(I3)
c___      NVMX = NVMAX + ICONST
c___      IF(NVMX .LE. MAXBST) GO TO 70
c___   60 WRITE(LOUT, 9050)
c___ 9050 FORMAT(' *** Too many, sorry, try again')
c___      GO TO 50
c___   70 IL = NVMX*(NVMX + 1)/2
c___      L = MIN(MAXL/IL, MAXR/NVMX, MAXBST)
c___      WRITE(LOUT, 9060) L, NVMAX
c___ 9060 FORMAT('+How many subsets of each size to be recorded ?'/
c___     +       ' Max. = ', I4, ' with NVMAX =', I3, ' : ')
c___      READ(LIN, 8010) NBEST
c___      IF(NBEST .GT. L) GO TO 60
c___C
c___C     Call TOLSET, SS & INITR to initialize arrays.
c___C
c___      IF (ICONST .EQ. 1) THEN
c___    CALL TOLSET(NCOLS, NRBAR, EL, U, TOL, WK, IER)
c___    CALL SS(NCOLS, EL, RHS, RESSQ, SSQ, IER)
c___    CALL INITR(NCOLS, NVMX, NBEST, BOUND, RESS, NVMX, LOPT, IL,
c___     +             VORDER, SSQ, IER)
c___      ELSE
c___    CALL TOLSET(NCOLS, NRBAR, EL(1), U, TOL(1), WK, IER)
c___    CALL SS(NCOLS, EL(1), RHS(1), RESSQ, SSQ(1), IER)
c___    CALL INITR(NCOLS, NVMX, NBEST, BOUND, RESS, NVMX, LOPT, IL,
c___     +             VORDER(1), SSQ(1), IER)
c___      END IF
c___      WRITE(LOUT, 9065) NCOLS, RESSQ
c___ 9065 FORMAT(' Initially NCOLS = ', I4,'  RESSQ = ', G13.5)
c___      IF (NOBS .GT. NCOLS) THEN
c___    NDF = NOBS - NCOLS
c___    VAR = RESSQ / NDF
c___    WRITE(*, 9068) VAR, NDF
c___ 9068   FORMAT(' Resid. variance estimate = ', g11.4, ' with ', i4,
c___     +         ' deg. of freedom'/)
c___      END IF
c___      IPROC = 0
c___C
c___C     Display menu & ask for choice.
c___C
c___  100 WRITE(LOUT, 9070)
c___ 9070 FORMAT(' Options:-'/
c___     1 ' C Corrlns. & partial corrlns.    F Forward selection'/
c___     2 ' B Backward elimination           R Sequential replacement'/
c___     3 ' E Efroymson stepwise             P Print summary of subsets'/
c___     4 ' I Specify IN variables           O Specify OUT variables'/
c___     5 ' L Least-squares regn.coeffs.     X Exhaustive search'/
c___     6 ' Q Quit           ENTER YOUR OPTION : ')
c___      READ(LIN, *) ANS
c___C
c___C     Compare ANS with currently available options.
c___C
c___      DO 110 IOPT = 1,22
c___    IF(ANS .EQ. OPTION(IOPT)) GO TO 120
c___  110 CONTINUE
c___      WRITE(LOUT, 9080) ANS
c___ 9080 FORMAT(' Option ', A, ' not available')
c___      GO TO 100
c___  120 L = (IOPT + 1)/2
c___C
c___C             C    F    B    R    E    P    I    O    L    X    Q
c___      GO TO (200, 300, 400, 500, 550, 700, 800, 900, 250, 600, 850), L
c___C
c___C-----------------------------------------------------------------------
c___C
c___C     Option 1. Correlations.
c___C
c___  200 WRITE(LOUT, 9200)
c___ 9200 FORMAT('+Do you want partial correlations ? (Y or N) ')
c___      NV = 0
c___      READ(LIN, *) ANS
c___      IF(ANS .EQ. 'N' .OR. ANS .EQ. 'n') GO TO 210
c___      IF(ANS .NE. 'Y' .AND. ANS .NE. 'y') GO TO 200
c___      ASSIGN 210 TO IRTN
c___      WRITE(LOUT, 9210)
c___ 9210 FORMAT(' Partial corrlns. on how many variables (excl.const.) ? ')
c___      READ(LIN, 8010) NV
c___      IF(NV .GT. 0) GO TO 1000
c___  210 WRITE(LOUT, 9220)
c___ 9220 FORMAT('+Correlations amongst all variables (A) or with Y only',
c___     +       1X,'(Y) ? ')
c___      IOPT = 0
c___      READ(LIN, *) ANS
c___      IF(ANS .EQ. 'A' .OR. ANS .EQ. 'a') IOPT = 1
c___      NB = NV + ICONST
c___      CALL PCORR(NCOLS, NRBAR, EL, U, RHS, RESSQ, NB, WK(UDIM+NCOLS+1),
c___     +           WK, IW, WK(UDIM+1), IER)
c___C
c___C     Display the (partial) correlations.
c___C     Correlations amongst the X-variables start at WK(1); correlations
c___C     with Y start at WK(UDIM+1).
c___C
c___      CALL PRINTC(NCOLS, NB, WK, UDIM, WK(UDIM+1), VORDER, VNAME(1),
c___     +            YNAME, IOPT, LOUT, IER)
c___      GO TO 100
c___C
c___C-----------------------------------------------------------------------
c___C
c___C     Option 9. Least - squares regression coefficients.
c___C
c___  250 WRITE(LOUT, 9850)
c___      READ(LIN, 8010) NV
c___      ASSIGN 260 TO IRTN
c___      GO TO 1000
c___  260 IF (ICONST .EQ. 1) THEN
c___    CALL REGCF(NCOLS, NRBAR, EL, U, RHS, TOL, WK, NV, IER)
c___      ELSE
c___    CALL REGCF(NCOLS, NRBAR, EL(1), U, RHS(1), TOL(1), WK, NV, IER)
c___      END IF
c___      IER = -IER
c___      IF(IER .NE. 0) WRITE(LOUT, 9250) IER
c___ 9250 FORMAT(' Variables linearly dependant, rank deficiency =',I4)
c___      WRITE(LOUT, 9260)(VORDER(I-ICONST),WK(I),I=1,NV)
c___ 9260 FORMAT(' Least-squares regn.coeffs.',
c___     +  7(/1X, I5, G13.5, 2X, I5, G13.5, 2X, I5, G13.5, 2X, I5, G13.5))
c___      WRITE(LOUT, 9270) SSQ(NV-ICONST)
c___ 9270 FORMAT(' Resid. sum of sq. =',G13.5/)
c___      GO TO 100
c___C
c___C-----------------------------------------------------------------------
c___C
c___C     Option 2. Forward selection.
c___C
c___  300 IF (ICONST .EQ. 1) THEN
c___    CALL FORWRD(NCOLS, NRBAR, EL, U, RHS, FIRST, LAST, VORDER, TOL,
c___     +      SSQ, BOUND, NVMX, RESS, NVMX, NBEST, LOPT, IL, WK, IW, IER)
c___      ELSE
c___    CALL FORWRD(NCOLS, NRBAR, EL(1), U, RHS(1), FIRST, LAST,
c___     +      VORDER(1), TOL(1), SSQ(1), BOUND, NVMX, RESS, NVMX, NBEST,
c___     +      LOPT, IL, WK, IW, IER)
c___      END IF
c___      NV = NVMX
c___      IF (IPROC .EQ. 2*(IPROC/2)) IPROC = IPROC + 1
c___      GO TO 1100
c___C
c___C-----------------------------------------------------------------------
c___C
c___C     Option 3. Backward elimination.
c___C
c___  400 IF (ICONST .EQ. 1) THEN
c___    CALL BAKWRD(NCOLS, NRBAR, EL, U, RHS, FIRST, LAST, VORDER, TOL,
c___     +     SSQ, BOUND, NVMX, RESS, NVMX, NBEST, LOPT, IL, WK, IW, IER)
c___      ELSE
c___    CALL BAKWRD(NCOLS, NRBAR, EL(1), U, RHS(1), FIRST, LAST,
c___     +     VORDER(1), TOL(1), SSQ(1), BOUND, NVMX, RESS, NVMX, NBEST,
c___     +     LOPT, IL, WK, IW, IER)
c___      END IF
c___      NV = LAST
c___      I = IPROC/2
c___      IF (I .EQ. 2*(I/2)) IPROC = IPROC + 2
c___      GO TO 1100
c___C
c___C-----------------------------------------------------------------------
c___C
c___C     Option 4. Sequential replacement.
c___C
c___  500 IF (ICONST .EQ. 1) THEN
c___    CALL SEQREP(NCOLS, NRBAR, EL, U, RHS, FIRST, LAST, VORDER, TOL,
c___     +     SSQ, BOUND, NVMX, RESS, NVMX, NBEST, LOPT, IL, WK, IW, IER)
c___      ELSE
c___    CALL SEQREP(NCOLS, NRBAR, EL(1), U, RHS(1), FIRST, LAST,
c___     +     VORDER(1), TOL(1), SSQ(1), BOUND, NVMX, RESS, NVMX, NBEST,
c___     +     LOPT, IL, WK, IW, IER)
c___      END IF
c___      I = IPROC/8
c___      IF (I .EQ. 2*(I/2)) IPROC = IPROC + 8
c___      NV = NVMX
c___      GO TO 1100
c___C
c___C-----------------------------------------------------------------------
c___C
c___C     Option 5. Efroymson (stepwise)
c___C
c___  550 WRITE(LOUT, 9550)
c___ 9550 FORMAT(' Enter F-to-enter value : ')
c___      READ(LIN, 8550) FIN
c___ 8550 FORMAT(F10.0)
c___      WRITE(LOUT, 9560)
c___ 9560 FORMAT(' Enter F-to-remove value : ')
c___      READ(LIN, 8550) FOUT
c___      IF (ICONST .EQ. 1) THEN
c___    CALL EFROYM(NCOLS, NRBAR, EL, U, RHS, FIRST, LAST, FIN, FOUT,
c___     +    SIZE, NOBS, VORDER, TOL, SSQ, BOUND, NVMX, RESS, NVMX, NBEST,
c___     +    LOPT, IL, WK, IW, IER)
c___      ELSE
c___    CALL EFROYM(NCOLS, NRBAR, EL(1), U, RHS(1), FIRST, LAST, FIN,
c___     +   FOUT, SIZE, NOBS, VORDER(1), TOL(1), SSQ(1), BOUND, NVMX, RESS,
c___     +   NVMX, NBEST, LOPT, IL, WK, IW, IER)
c___      END IF
c___      IF (IER .NE. 0) THEN
c___    WRITE(LOUT, 9570) IER
c___ 9570   FORMAT(' Error code',I4,' returned by EFROYM')
c___    GO TO 100
c___      ELSE
c___    NV = SIZE
c___    I = IPROC/4
c___    IPROC = IPROC + 4
c___    GO TO 1100
c___      END IF
c___C
c___C-----------------------------------------------------------------------
c___C
c___C     Option 10. Exhaustive search.
c___C
c___  600  IF (ICONST .EQ. 1) THEN
c___     CALL XHAUST(NCOLS, NRBAR, EL, U, RHS, FIRST, LAST, VORDER, TOL,
c___     +      SSQ, BOUND, NVMX, RESS, NVMX, NBEST, LOPT, IL, WK, IW, IWK,
c___     +      IIW, IER)
c___      ELSE
c___     CALL XHAUST(NCOLS, NRBAR, EL(1), U, RHS(1), FIRST, LAST,
c___     +      VORDER(1), TOL(1), +   SSQ(1), BOUND, NVMX, RESS, NVMX,
c___     +      NBEST, LOPT, IL, WK, IW, IWK, IIW, IER)
c___      END IF
c___      IF (IPROC .LT. 16) IPROC = IPROC + 16
c___      GO TO 100
c___C
c___C-----------------------------------------------------------------------
c___C
c___C     Option 6. Print summary of best subsets found so far.
c___C
c___  700 CALL LSORT(LOPT, IL, NBEST, NVMX)
c___      L = FIRST*(FIRST-1)/2 + 1
c___      LINE = 1
c___      M = FIRST - ICONST
c___      DO 730 NV = FIRST, NVMX
c___    WRITE(LOUT,9700) M
c___ 9700   FORMAT(20X,'Best subsets found of',I3,' variables')
c___    LINE = LINE + 1
c___    DO 720 NB = 1,NBEST
c___      J = (NB-1)*NVMX + NV
c___      TEMP = RESS(J)
c___      IF(TEMP .GT. 1.E+35) GO TO 720
c___      IPOS = L
c___      DO 710 I = 1,NV
c___        J = (NB-1)*IL + IPOS
c___        IWK(I) = LOPT(J)
c___        IPOS = IPOS + 1
c___  710     CONTINUE
c___        WRITE(LOUT,9710) TEMP,(IWK(I),I=FIRST,NV)
c___      LINE = LINE + 1 + (NV-1)/10
c___ 9710     FORMAT(' RSS =',G14.6,3X,'Variables:',10I4,4(/10X,10I4))
c___  720   CONTINUE
c___    IF (LINE .GE. 25 - NB) THEN
c___      PAUSE
c___      LINE = 1
c___    END IF
c___    L = L + NV
c___    M = M + 1
c___  730 CONTINUE
c___      GO TO 100
c___C
c___C----------------------------------------------------------------------
c___C
c___C     Option 7. Force variables into models.
c___C
c___  800 WRITE(LOUT, 9800)
c___ 9800 FORMAT('+How many variables, excl. constant ? ')
c___      READ(LIN, 8010) NV
c___      ASSIGN 810 TO IRTN
c___      GO TO 1000
c___  810 GO TO 1100
c___C
c___C----------------------------------------------------------------------
c___C
c___C     Option 11. Exit.
c___C
c___  850 IF (IPROC .EQ. 0) STOP
c___      WRITE(LOUT, 9860)
c___ 9860 FORMAT(' Do you want to save the best subsets found ? (Y/N) ')
c___      READ(LIN, *) ANS
c___      IF (ANS .EQ. 'Y' .OR. ANS .EQ. 'y') THEN
c___    REWIND(9)
c___    CALL LSORT(LOPT, IL, NBEST, NVMX)
c___    READ(9) K, ICONST, NCOLS, NOBS, NRBAR, LSEL
c___    IF (ICONST .EQ. 0) THEN
c___      READ(9) YNAME, (VNAME(I),I=1,K)
c___      READ(9) (U(I),I=1,NRBAR), (EL(I),I=1,K), (RHS(I),I=1,K), RESSQ
c___    ELSE
c___      READ(9) YNAME, (VNAME(I),I=0,K)
c___      READ(9) (U(I),I=1,NRBAR), (EL(I),I=0,K), (RHS(I),I=0,K), RESSQ
c___    END IF
c___    LSEL = .TRUE.
c___    REWIND(9)
c___    ILNB = IL*NBEST
c___    IR = NVMX*NBEST
c___    WRITE(9) K, ICONST, NCOLS, NOBS, NRBAR, LSEL
c___    IF (ICONST .EQ. 0) THEN
c___      WRITE(9) YNAME, (VNAME(I),I=1,K)
c___      WRITE(9) (U(I),I=1,NRBAR), (EL(I),I=1,K), (RHS(I),I=1,K),
c___     +              RESSQ
c___    ELSE
c___      WRITE(9) YNAME, (VNAME(I),I=0,K)
c___      WRITE(9) (U(I),I=1,NRBAR), (EL(I),I=0,K), (RHS(I),I=0,K),
c___     +              RESSQ
c___    END IF
c___    WRITE(9) NVMAX, NBEST, IL, ILNB, IR, IPROC
c___    WRITE(9) (LOPT(L),L=1,ILNB)
c___    WRITE(9) (RESS(I),I=1,IR)
c___      END IF
c___      STOP
c___C
c___C----------------------------------------------------------------------
c___C
c___C     Simulated subroutine to force variables into starting positions.
c___C     NV = no. of variables to be forced in.
c___C
c___ 1000 WRITE(LOUT, 9930) (I, VNAME(I),I = 1,K)
c___ 9930 FORMAT('+Variables & their numbers:', 10(/1X, 5(I3, 1X, A8, 3X)))
c___      IF(NV .LE. 0) GO TO 100
c___      WRITE(LOUT, 9920)
c___ 9920 FORMAT(' List variable nos. : ')
c___      READ(LIN, *) (IWK(I),I = 1,NV)
c___C
c___C     Find variables in VORDER which are in the input list and move up
c___C     to the next available position.
c___C
c___      IF (ICONST .EQ. 1) THEN
c___    CALL REORDR(NCOLS, NRBAR, VORDER, EL, U, RHS, SSQ, TOL, IWK, NV,
c___     +            2, IER)
c___      ELSE
c___    CALL REORDR(NCOLS, NRBAR, VORDER(1), EL(1), U, RHS(1), SSQ(1),
c___     +              TOL(1), IWK, NV, 1, IER)
c___      END IF
c___      NV = NV + ICONST
c___      FIRST = NV + 1
c___      GO TO IRTN,(210, 260, 810)
c___C
c___C----------------------------------------------------------------------
c___C
c___C     Option 8. Force variables out of models.
c___C
c___  900 WRITE(LOUT, 9850)
c___ 9850 FORMAT('+How many variables ? ')
c___      READ(LIN, 8010) NV
c___      WRITE(LOUT, 9920)
c___      DO 910 I = 1, NV
c___  910 READ(LIN, *) IWK(I)
c___      LAST = NCOLS
c___      J = LAST
c___  920 L = VORDER(J)
c___      DO 930 M = 1, NV
c___    IF(L .EQ. IWK(M)) GO TO 940
c___  930 CONTINUE
c___      GO TO 960
c___  940 IF(J .EQ. LAST) GO TO 950
c___      CALL VMOVE(NCOLS, NRBAR, VORDER, EL, U, RHS, SSQ, J, LAST, TOL,
c___     +    IER)
c___  950 LAST = LAST - 1
c___      IF(J .LT. FIRST) FIRST = FIRST - 1
c___  960 J = J - 1
c___      IF(J .GT. 0) GO TO 920
c___      GO TO 100
c___C
c___C----------------------------------------------------------------------
c___C
c___C     Print current order of the first NV variables and their RSS's.
c___C
c___ 1100 WRITE(LOUT, 9900)
c___ 9900 FORMAT(' Order  Variable   Resid.sumsq.')
c___      DO 1110 I = 1-ICONST, NV-ICONST
c___    J = VORDER(I)
c___    WRITE(LOUT, 9910) I, VNAME(J), SSQ(I)
c___ 9910   FORMAT(I5, 3X, A8, 1X, G14.6)
c___ 1110 CONTINUE
c___      GO TO 100
c___      END
c___C
c___C
c___C
c___C
      SUBROUTINE LSORT(LOPT, IL, NBEST, NVMX)
C
C      Sort the variable numbers in LOPT into increasing order.
C
C      Latest revision - 12 February 1986
C
      DIMENSION LOPT(IL, NBEST)
      INTEGER COL, SIZE, TEMP, START
C
      IF (NVMX .LT. 2) RETURN
      DO 20 COL = 1, NBEST
        TEMP = LOPT(2, COL)
        IF (TEMP .GT. LOPT(3,COL)) THEN
          LOPT(2,COL) = LOPT(3,COL)
          LOPT(3,COL) = TEMP
        END IF
        IF (IL .LE. 3) GO TO 20
        START = 4
        DO 10 SIZE = 3, NVMX
          CALL SHELL(LOPT(START,COL), SIZE)
          START = START + SIZE
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
C
C
C
C
      SUBROUTINE SHELL(L, N)
C
C      Perform a SHELL-sort on integer array L, sorting into
C      increasing order.
C
C      Latest revision - 12 February 1986
C
      DIMENSION L(N)
      INTEGER START, END, TEMP
C
      INCR = N
   10 INCR = INCR/3
      IF (INCR .EQ. 2*(INCR/2)) INCR = INCR + 1
      DO 50 START = 1, INCR
        END = N
C
C      TEMP contains the element being compared; IT holds its current
C      location.   It is compared with the elements in locations
C      IT+INCR, IT+2.INCR, ... until a larger element is found.   All
C      smaller elements move INCR locations towards the start.   After
C      each time through the sequence, the END is decreased by INCR
C      until END <= INCR.
C
   20   I1 = START
        TEMP = L(I1)
        IT = I1
C
C      I2 = location of element NEW to be compared with TEMP.
C      Test I2 <= END.
C
   30   I2 = I1 + INCR
        IF (I2 .GT. END) THEN
          IF (I1 .GT. IT) L(I1) = TEMP
          END = END - INCR
          GO TO 40
        END IF
        NEW = L(I2)
C
C      If TEMP > NEW, move NEW to lower-numbered position.
C
        IF (TEMP .GT. NEW) THEN
          L(I1) = NEW
          I1 = I2
          GO TO 30
        END IF
C
C      TEMP <= NEW so do not swap.
C      Use NEW as the next TEMP.
C
        IF (I1 .GT. IT) L(I1) = TEMP
        I1 = I2
        TEMP = NEW
        IT = I1
        GO TO 30
C
C      Repeat until END <= INCR.
C
   40   IF (END .GT. INCR) GO TO 20
   50 CONTINUE
C
C      Repeat until INCR = 1.
C
      IF (INCR .GT. 1) GO TO 10
      RETURN
      END
C
C
C
c___      SUBROUTINE PRINTC(NP, IN, CORMAT, DIMC, YCORR, VORDER, VNAME,
c___     +              YNAME, IOPT, LOUT, IER)
c___C
c___C     Print (partial) correlations calculated using PCORR.
c___C     If IOPT = 0, print correlations with the Y-variable only.
c___C
c___C     IMPLICIT NONE
c___      INTEGER NP, IN, DIMC, VORDER(NP), IOPT, LOUT, IER
c___      DOUBLE PRECISION CORMAT(DIMC), YCORR(NP)
c___      CHARACTER VNAME(NP)*8, YNAME*8
c___C
c___C     Local variables.
c___C
c___      INTEGER NROWS, J1, J2, J, I1, I2, I, ROW, UPOS, TPOS, LAST
c___      CHARACTER TEXT*74, EMPTY*65, CHAR1*9
c___
c___      DATA EMPTY/' '/, CHAR1/' 1.0'/
c___C
c___C     Check validity of arguments
c___C
c___      IER = 0
c___      IF (IN .GE. NP) IER = 1
c___      IF (NP .LE. 1) IER = IER + 2
c___      NROWS = NP - IN
c___      IF (DIMC .LE. NROWS*(NROWS-1)/2) IER = IER + 4
c___      IF (IER .NE. 0) RETURN
c___C
c___C     If IOPT.NE.0 output heading
c___C
c___      IF(IOPT .EQ. 0) GO TO 30
c___      WRITE(LOUT, 900)
c___  900 FORMAT(/5X, 'Correlation matrix')
c___      J1 = IN + 1
c___   10 J2 = MIN(J1+6, NP)
c___      I1 = J1 - IN
c___      I2 = J2 - IN
c___      WRITE(LOUT, 910)(VNAME(VORDER(J)), J=J1,J2)
c___  910 FORMAT(11X, 7(A8, 1X))
c___C
c___C     Print correlations for rows 1 to I2, columns I1 to I2.
c___C
c___      DO 20 ROW = 1, I2
c___      TEXT = ' ' // VNAME(VORDER(ROW+IN)) // EMPTY
c___      IF (I1 .GT. ROW) THEN
c___    UPOS = (ROW-1) * (NROWS+NROWS-ROW) /2 + (I1-ROW)
c___    LAST = UPOS + I2 - I1
c___    WRITE(TEXT(12:74), '(7(F8.5,1X))')(CORMAT(I),I=UPOS,LAST)
c___      ELSE
c___    UPOS = (ROW-1) * (NROWS+NROWS-ROW) /2 + 1
c___    TPOS = 12 + 9*(ROW-I1)
c___    TEXT(TPOS:TPOS+8) = CHAR1
c___    LAST = UPOS + I2 - ROW - 1
c___    IF (ROW .LT. I2) WRITE(TEXT(TPOS+9:74), '(6(F8.5, 1X))')
c___     +                 (CORMAT(I),I=UPOS, LAST)
c___      END IF
c___      WRITE(LOUT, '(A)') TEXT
c___   20 CONTINUE
c___C
c___C     Move onto the next block of columns.
c___C
c___      J1 = J2 + 1
c___      IF (J1 .LE. NP) GO TO 10
c___C
c___C     Correlations with the Y-variable.
c___C
c___   30 WRITE(LOUT, 920) YNAME
c___  920 FORMAT(/5X, 'Correlations with the dependent variable: ', A)
c___      J1 = IN + 1
c___   40 J2 = MIN(J1+7, NP)
c___      I1 = J1 - IN
c___      I2 = J2 - IN
c___      WRITE(LOUT, 930)(VNAME(VORDER(J)), J=J1,J2)
c___  930 FORMAT(/1X, 8(A8, 1X))
c___      WRITE(LOUT, 940)(YCORR(I),I=I1,I2)
c___  940 FORMAT(1X, 8(F8.5, 1X))
c___      J1 = J2 + 1
c___      IF (J1 .LE. NP) GO TO 40
c___C
c___C     Put extra blank line into output
c___C
c___      WRITE(LOUT, *)
c___C
c___      RETURN
c___      END

c     include 'subs.for'
