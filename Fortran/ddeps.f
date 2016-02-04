* COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils

C Original date 29 Jan 2001
C 29 January 2001. Modified from MC49 to be threadsafe.

C 12th July 2004 Version 1.0.0. Version numbering added.
C 28 February 2008. Version 1.0.1. Comments flowed to column 72.
C 21 September 2009. Version 1.0.2. Minor change to documentation.

      SUBROUTINE MC59AD(ICNTL,NC,NR,NE,IRN,LJCN,JCN,LA,A,LIP,IP,
     &                  LIW,IW,INFO)
C
C To sort the sparsity pattern of a matrix to an ordering by columns.
C There is an option for ordering the entries within each column by
C increasing row indices and an option for checking the user-supplied
C matrix entries for indices which are out-of-range or duplicated.
C
C ICNTL:  INTEGER array of length 10. Intent(IN). Used to specify
C         control parameters for the subroutine.
C ICNTL(1): indicates whether the user-supplied matrix entries are to
C           be checked for duplicates, and out-of-range indices.
C           Note  simple checks are always performed.
C           ICNTL(1) = 0, data checking performed.
C           Otherwise, no data checking.
C ICNTL(2): indicates the ordering requested.
C           ICNTL(2) = 0, input is by rows and columns in arbitrary
C           order and the output is sorted by columns.
C           ICNTL(2) = 1, the output is also row ordered
C           within each column.
C           ICNTL(2) = 2, the input is already ordered by
C           columns and is to be row ordered within each column.
C           Values outside the range 0 to 2 are flagged as an error.
C ICNTL(3): indicates whether matrix entries are also being ordered.
C           ICNTL(3) = 0, matrix entries are ordered.
C           Otherwise, only the sparsity pattern is ordered
C           and the array A is not accessed by the routine.
C ICNTL(4): the unit number of the device to
C           which error messages are sent. Error messages
C           can be suppressed by setting ICNTL(4) < 0.
C ICNTL(5): the unit number of the device to
C           which warning messages are sent. Warning
C           messages can be suppressed by setting ICNTL(5) < 0.
C ICNTL(6)  indicates whether matrix symmetric. If unsymmetric, ICNTL(6)
C           must be set to 0.
C           If ICNTL(6) = -1 or 1, symmetric and only the lower
C           triangular part of the reordered matrix is returned.
C           If ICNTL(6) = -2 or 2, Hermitian and only the lower
C           triangular part of the reordered matrix is returned.
C           If error checks are performed (ICNTL(1) = 0)
C           and ICNTL(6)> 1 or 2, the values of duplicate
C           entries are added together; if ICNTL(6) < -1 or -2, the
C           value of the first occurrence of the entry is used.
C ICNTL(7) to ICNTL(10) are not currently accessed by the routine.
C
C NC:      INTEGER variable. Intent(IN). Must be set by the user
C          to the number of columns in the matrix.
C NR:      INTEGER variable. Intent(IN). Must be set by the user
C          to the number of rows in the matrix.
C NE:      INTEGER variable. Intent(IN). Must be set by the user
C          to the number of entries in the matrix.
C IRN: INTEGER array of length NE. Intent (INOUT). Must be set by the
C            user to hold the row indices of the entries in the matrix.
C          If ICNTL(2).NE.2, the entries may be in any order.
C          If ICNTL(2).EQ.2, the entries in column J must be in
C            positions IP(J) to IP(J+1)-1 of IRN. On exit, the row
C            indices are reordered so that the entries of a single
C            column are contiguous with column J preceding column J+1, J
C            = 1, 2, ..., NC-1, with no space between columns.
C          If ICNTL(2).EQ.0, the order within each column is arbitrary;
C            if ICNTL(2) = 1 or 2, the order within each column is by
C            increasing row indices.
C LJCN:    INTEGER variable. Intent(IN). Defines length array
C JCN:     INTEGER array of length LJCN. Intent (INOUT).
C          If ICNTL(2) = 0 or 1, JCN(K) must be set by the user
C          to the column index of the entry
C          whose row index is held in IRN(K), K = 1, 2, ..., NE.
C          On exit, the contents of this array  will have been altered.
C          If ICNTL(2) = 2, the array is not accessed.
C LA:      INTEGER variable. Intent(IN). Defines length of array
C          A.
C A:       is a REAL (DOUBLE PRECISION in the D version, INTEGER in
C          the I version, COMPLEX in the C version,
C          or COMPLEX"*"16 in the Z version) array of length LA.
C          Intent(INOUT).
C          If ICNTL(3).EQ.0, A(K) must be set by the user to
C          hold the value of the entry with row index IRN(K),
C          K = 1, 2, ..., NE. On exit, the array will have been
C          permuted in the same way as the array IRN.
C          If ICNTL(3).NE.0, the array is not accessed.
C LIP:     INTEGER variable. Intent(IN). Defines length of array
C          IP.
C IP:      INTEGER array of length LIP. Intent(INOUT). IP
C          need only be set by the user if ICNTL(2) = 2.
C          In this case, IP(J) holds the position in
C          the array IRN of the first entry in column J, J = 1, 2,
C          ..., NC, and IP(NC+1) is one greater than the number of
C          entries in the matrix.
C          In all cases, the array IP will have this meaning on exit
C          from the subroutine and is altered when ICNTL(2) = 2 only
C          when ICNTL(1) =  0 and there are out-of-range
C          indices or duplicates.
C LIW:     INTEGER variable. Intent(IN). Defines length of array
C          IW.
C IW:      INTEGER array of length LIW. Intent(OUT). Used by the
C          routine as workspace.
C INFO:    INTEGER array of length 10.  Intent(OUT). On exit,
C          a negative value of INFO(1) is used to signal a fatal
C          error in the input data, a positive value of INFO(1)
C          indicates that a warning has been issued, and a
C          zero value is used to indicate a successful call.
C          In cases of error, further information is held in INFO(2).
C          For warnings, further information is
C          provided in INFO(3) to INFO(6).  INFO(7) to INFO(10) are not
C          currently used and are set to zero.
C          Possible nonzero values of INFO(1):
C         -1 -  The restriction ICNTL(2) = 0, 1, or 2 violated.
C               Value of ICNTL(2) is given by INFO(2).
C         -2 -  NC.LE.0. Value of NC is given by INFO(2).
C         -3 -  Error in NR. Value of NR is given by INFO(2).
C         -4 -  NE.LE.0. Value of NE is given by INFO(2).
C         -5 -  LJCN too small. Min. value of LJCN is given by INFO(2).
C         -6 -  LA too small. Min. value of LA is given by INFO(2).
C         -7 -  LIW too small. Value of LIW is given by INFO(2).
C         -8 -  LIP too small. Value of LIP is given by INFO(2).
C         -9 -  The entries of IP not monotonic increasing.
C        -10 -  For each I, IRN(I) or JCN(I) out-of-range.
C        -11 -  ICNTL(6) is out of range.
C         +1 -  One or more duplicated entries. One copy of
C               each such entry is kept and, if ICNTL(3) = 0 and
C               ICNTL(6).GE.0, the values of these entries are
C               added together. If  ICNTL(3) = 0 and ICNTL(6).LT.0,
C               the value of the first occurrence of the entry is used.
C               Initially INFO(3) is set to zero. If an entry appears
C               k times, INFO(3) is incremented by k-1 and INFO(6)
C               is set to the revised number of entries in the
C               matrix.
C         +2 - One or more of the entries in IRN out-of-range. These
C               entries are removed by the routine.`INFO(4) is set to
C               the number of entries which were out-of-range and
C               INFO(6) is set to the revised number of entries in the
C               matrix.
C         +4 - One or more of the entries in JCN out-of-range. These
C               entries are removed by the routine. INFO(5) is set to
C               the number of entries which were out-of-range and
C               INFO(6) is set to the revised number of entries in the
C               matrix. Positive values of INFO(1) are summed so that
C               the user can identify all warnings.
C
C     .. Scalar Arguments ..
      INTEGER LA,LIP,LIW,LJCN,NC,NE,NR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER ICNTL(10),IP(LIP),INFO(10),IRN(NE),IW(LIW),JCN(LJCN)
C     ..
C     .. Local Scalars ..
      INTEGER I,ICNTL1,ICNTL2,ICNTL3,ICNTL6,LAA
      INTEGER IDUP,IOUT,IUP,JOUT,LP,MP,KNE,PART
      LOGICAL LCHECK
C     ..
C     .. External Subroutines ..
      EXTERNAL MC59BD,MC59CD,MC59DD,MC59ED,MC59FD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
C     .. Executable Statements ..

C Initialise
      DO 10 I = 1,10
         INFO(I) = 0
   10 CONTINUE

      ICNTL1 = ICNTL(1)
      ICNTL2 = ICNTL(2)
      ICNTL3 = ICNTL(3)
      ICNTL6 = ICNTL(6)
      LCHECK = (ICNTL1.EQ.0)
C Streams for errors/warnings
      LP = ICNTL(4)
      MP = ICNTL(5)

C  Check the input data
      IF (ICNTL2.GT.2 .OR. ICNTL2.LT.0) THEN
         INFO(1) = -1
         INFO(2) = ICNTL2
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9010) ICNTL2
         END IF
         GO TO 70
      END IF

      IF (ICNTL6.GT.2 .OR. ICNTL6.LT.-2) THEN
         INFO(1) = -11
         INFO(2) = ICNTL6
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9150) ICNTL6
         END IF
         GO TO 70
      END IF
C For real matrices, symmetric = Hermitian so only
C have to distinguish between unsymmetric (ICNTL6 = 0) and
C symmetric (ICNTL6.ne.0)

      IF (NC.LT.1) THEN
        INFO(1) = -2
        INFO(2) = NC
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9020) NC
        END IF
        GO TO 70
      END IF

      IF (NR.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9030) NR
        END IF
        GO TO 70
      END IF

      IF (ICNTL6.NE.0 .AND. NR.NE.NC) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9035) NC,NR
        END IF
        GO TO 70
      END IF

      IF (NE.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NE
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9040) NE
        END IF
        GO TO 70
      END IF

      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.1) THEN
        IF (LJCN.LT.NE) THEN
          INFO(1) = -5
          INFO(2) = NE
        END IF
      ELSE
        IF (LJCN.LT.1) THEN
          INFO(1) = -5
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-5) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9050) LJCN,INFO(2)
         END IF
         GO TO 70
      END IF

      IF (ICNTL3.EQ.0) THEN
        IF (LA.LT.NE) THEN
          INFO(1) = -6
          INFO(2) = NE
        END IF
      ELSE
        IF (LA.LT.1) THEN
          INFO(1) = -6
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-6) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9060) LA,INFO(2)
         END IF
         GO TO 70
      END IF

      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.2) THEN
        IF (LIP.LT.NC+1) THEN
          INFO(1) = -7
          INFO(2) = NC+1
        END IF
      ELSE IF (LIP.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -7
        INFO(2) = MAX(NR,NC)+1
      END IF
      IF (INFO(1).EQ.-7) THEN
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9065) LIP,INFO(2)
        END IF
        GO TO 70
      END IF

C Check workspace sufficient
      IF (LIW.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -8
        INFO(2) = MAX(NR,NC)+1
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9070) LIW,INFO(2)
        END IF
        GO TO 70
      END IF

      LAA = NE
      IF (ICNTL3.NE.0) LAA = 1
C Initialise counts of number of out-of-range entries and duplicates
      IOUT = 0
      JOUT = 0
      IDUP = 0
      IUP = 0

C PART is used by MC59BD to indicate if upper or lower or
C all of matrix is required.
C PART =  0 : unsymmetric case, whole matrix wanted
C PART =  1 : symmetric case, lower triangular part of matrix wanted
C PART = -1 : symmetric case, upper triangular part of matrix wanted
      PART = 0
      IF (ICNTL6.NE.0) PART = 1

      IF (ICNTL2.EQ.0) THEN

C Order directly by columns
C On exit from MC59BD, KNE holds number of entries in matrix
C after removal of out-of-range entries. If no data checking, KNE = NE.
        CALL MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LAA,A,IP,IW,
     +              IOUT,JOUT,KNE)
C Return if ALL entries out-of-range.
        IF (KNE.EQ.0) GO TO 50

C Check for duplicates
        IF (LCHECK) CALL MC59ED(NC,NR,NE,IRN,LIP,IP,LAA,A,IW,IDUP,
     &                          KNE,ICNTL6)

      ELSE IF (ICNTL2.EQ.1) THEN

C First order by rows.
C Interchanged roles of IRN and JCN, so set PART = -1
C if matrix is symmetric case
        IF (ICNTL6.NE.0) PART = -1
        CALL MC59BD(LCHECK,PART,NR,NC,NE,JCN,IRN,LAA,A,IW,IP,
     +              JOUT,IOUT,KNE)
C Return if ALL entries out-of-range.
        IF (KNE.EQ.0) GO TO 50

C At this point, JCN and IW hold column indices and row pointers
C Optionally, check for duplicates.
        IF (LCHECK) CALL MC59ED(NR,NC,NE,JCN,NR+1,IW,LAA,A,IP,
     &                          IDUP,KNE,ICNTL6)

C Now order by columns and by rows within each column
        CALL MC59CD(NC,NR,KNE,IRN,JCN,LAA,A,IP,IW)

      ELSE IF (ICNTL2.EQ.2) THEN
C Input is using IP, IRN.
C Optionally check for duplicates and remove out-of-range entries
        IF (LCHECK) THEN
          CALL MC59FD(NC,NR,NE,IRN,NC+1,IP,LAA,A,LIW,IW,IDUP,
     +                IOUT,IUP,KNE,ICNTL6,INFO)
C Return if IP not monotonic.
          IF (INFO(1).EQ.-9) GO TO 40
C Return if ALL entries out-of-range.
          IF (KNE.EQ.0) GO TO 50
        ELSE
           KNE = NE
        END IF

C  Order by rows within each column
        CALL MC59DD(NC,KNE,IRN,IP,LAA,A)

      END IF

      INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(6) = KNE
      INFO(7) = IUP

C Set warning flag if out-of-range /duplicates found
      IF (IDUP.GT.0) INFO(1) = INFO(1) + 1
      IF (IOUT.GT.0) INFO(1) = INFO(1) + 2
      IF (JOUT.GT.0) INFO(1) = INFO(1) + 4
      IF (INFO(1).GT.0 .AND. MP.GT.0) THEN
        WRITE (MP,FMT=9080) INFO(1)
        IF (IOUT.GT.0) WRITE (MP,FMT=9090) IOUT
        IF (JOUT.GT.0) WRITE (MP,FMT=9110) JOUT
        IF (IDUP.GT.0) WRITE (MP,FMT=9100) IDUP
        IF (IUP.GT.0)  WRITE (MP,FMT=9130) IUP
      END IF
      GO TO 70

   40 INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(7) = IUP
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9140)
      END IF
      GO TO 70

   50 INFO(1) = -10
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(2) = IOUT + JOUT
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9120)
      END IF
   70 RETURN

 9000 FORMAT (/,' *** Error return from MC59AD *** INFO(1) = ',I3)
 9010 FORMAT (1X,'ICNTL(2) = ',I2,' is out of range')
 9020 FORMAT (1X,'NC = ',I6,' is out of range')
 9030 FORMAT (1X,'NR = ',I6,' is out of range')
 9035 FORMAT (1X,'Symmetric case. NC = ',I6,' but NR = ',I6)
 9040 FORMAT (1X,'NE = ',I10,' is out of range')
 9050 FORMAT (1X,'Increase LJCN from ',I10,' to at least ',I10)
 9060 FORMAT (1X,'Increase LA from ',I10,' to at least ',I10)
 9065 FORMAT (1X,'Increase LIP from ',I8,' to at least ',I10)
 9070 FORMAT (1X,'Increase LIW from ',I8,' to at least ',I10)
 9080 FORMAT (/,' *** Warning message from MC59AD *** INFO(1) = ',I3)
 9090 FORMAT (1X,I8,' entries in IRN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9100 FORMAT (1X,I8,' duplicate entries were supplied by the user')
 9110 FORMAT (1X,I8,' entries in JCN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9120 FORMAT (1X,'All entries out of range')
 9130 FORMAT (1X,I8,' of these entries were in the upper triangular ',
     +       /,'       part of matrix')
 9140 FORMAT (1X,'Entries in IP are not monotonic increasing')
 9150 FORMAT (1X,'ICNTL(6) = ',I2,' is out of range')
      END
C***********************************************************************
      SUBROUTINE MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LA,A,IP,IW,IOUT,
     +                  JOUT,KNE)
C
C   To sort a sparse matrix from arbitrary order to
C   column order, unordered within each column. Optionally
C   checks for out-of-range entries in IRN,JCN.
C
C LCHECK - logical variable. Intent(IN). If true, check
C          for out-of-range indices.
C PART -   integer variable. Intent(IN)
C PART =  0 : unsymmetric case, whole matrix wanted
C PART =  1 : symmetric case, lower triangular part of matrix wanted
C             (ie IRN(K) .ge. JCN(K) on exit)
C PART = -1 : symmetric case, upper triangular part of matrix wanted
C             (ie IRN(K) .le. JCN(K) on exit)
C   NC - integer variable. Intent(IN)
C      - on entry must be set to the number of columns in the matrix
C   NR - integer variable. Intent(IN)
C      - on entry must be set to the number of rows in the matrix
C   NE - integer variable. Intent(IN)
C      - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(INOUT)
C      - on entry set to contain the row indices of the nonzeros
C        in arbitrary order.
C      - on exit, the entries in IRN are reordered so that the row
C        indices for column 1 precede those for column 2 and so on,
C        but the order within columns is arbitrary.
C  JCN - integer array of length NE. Intent(INOUT)
C      - on entry set to contain the column indices of the nonzeros
C      - JCN(K) must be the column index of
C        the entry in IRN(K)
C      - on exit, JCN(K) is the column index for the entry with
C        row index IRN(K) (K=1,...,NE).
C  LA  - integer variable which defines the length of the array A.
C        Intent(IN)
C   A  - real (double precision/complex/complex*16) array of length LA
C        Intent(INOUT)
C      - if LA > 1, the array must be of length NE, and A(K)
C        must be set to the value of the entry in (IRN(K), JCN(K));
C        on exit A is reordered in the same way as IRN
C      - if LA = 1, the array is not accessed
C  IP  - integer array of length NC+1. Intent(INOUT)
C      - not set on entry
C      - on exit, IP(J) contains the position in IRN (and A) of the
C        first entry in column J (J=1,...,NC)
C      - IP(NC+1) is set to NE+1
C  IW  - integer array of length NC+1.  Intent(INOUT)
C      - the array is used as workspace
C      - on exit IW(I) = IP(I) (so IW(I) points to the beginning
C        of column I).
C IOUT - integer variable. Intent(OUT). On exit, holds number
C        of entries in IRN found to be out-of-range
C JOUT - integer variable. Intent(OUT). On exit, holds number
C        of entries in JCN found to be out-of-range
C  KNE - integer variable. Intent(OUT). On exit, holds number
C        of entries in matrix after removal of out-of-range entries.
C        If no data checking, KNE = NE.

C     .. Scalar Arguments ..
      INTEGER LA,NC,NE,NR,IOUT,JOUT,KNE,PART
      LOGICAL LCHECK
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NC+1),JCN(NE)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,JCE,JCEP,K,L,LOC
C     ..
C     .. Executable Statements ..

C Initialise IW
      DO 10 J = 1,NC + 1
        IW(J) = 0
   10 CONTINUE

      KNE = 0
      IOUT = 0
      JOUT = 0
C Count the number of entries in each column and store in IW.
C We also allow checks for out-of-range indices
      IF (LCHECK) THEN
C Check data.
C Treat case of pattern only separately.
        IF (LA.GT.1) THEN
          IF (PART.EQ.0) THEN
C Unsymmetric
            DO 20 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                A(KNE) = A(K)
                IW(J) = IW(J) + 1
              END IF
   20       CONTINUE
          ELSE IF (PART.EQ.1) THEN
C Symmetric, lower triangle
            DO 21 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
C Lower triangle ... swap if necessary
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   21       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
C Symmetric, upper triangle
            DO 22 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
C Upper triangle ... swap if necessary
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   22       CONTINUE
          END IF
        ELSE
C Pattern only
          IF (PART.EQ.0) THEN
            DO 25 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                IW(J) = IW(J) + 1
              END IF
   25       CONTINUE
          ELSE IF (PART.EQ.1) THEN
            DO 26 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
C Lower triangle ... swap if necessary
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   26       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
            DO 27 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
C Upper triangle ... swap if necessary
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   27       CONTINUE
          END IF
        END IF
C Return if ALL entries out-of-range.
        IF (KNE.EQ.0) GO TO 130

      ELSE

C No checks
        KNE = NE
        IF (PART.EQ.0) THEN
          DO 30 K = 1,NE
            J = JCN(K)
            IW(J) = IW(J) + 1
   30     CONTINUE
        ELSE IF (PART.EQ.1) THEN
          DO 35 K = 1,NE
            I = IRN(K)
            J = JCN(K)
C Lower triangle ... swap if necessary
            IF (I.LT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   35     CONTINUE
        ELSE IF (PART.EQ.-1) THEN
          DO 36 K = 1,NE
            I = IRN(K)
            J = JCN(K)
C Upper triangle ... swap if necessary
            IF (I.GT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   36     CONTINUE
        END IF
      END IF

C KNE is now the number of nonzero entries in matrix.

C Put into IP and IW the positions where each column
C would begin in a compressed collection with the columns
C in natural order.

      IP(1) = 1
      DO 37 J = 2,NC + 1
        IP(J) = IW(J-1) + IP(J-1)
        IW(J-1) = IP(J-1)
   37 CONTINUE

C Reorder the elements into column order.
C Fill in each column from the front, and as a new entry is placed
C in column K increase the pointer IW(K) by one.

      IF (LA.EQ.1) THEN
C Pattern only
        DO 70 L = 1,NC
          DO 60 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            DO 40 J = 1,NE
              IF (JCE.EQ.L) GO TO 50
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
   40       CONTINUE
   50       JCN(K) = JCE
            IRN(K) = ICE
   60     CONTINUE
   70   CONTINUE
      ELSE

        DO 120 L = 1,NC
          DO 110 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            ACE = A(K)
            DO 90 J = 1,NE
              IF (JCE.EQ.L) GO TO 100
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
              ACEP = A(LOC)
              A(LOC) = ACE
              ACE = ACEP
   90       CONTINUE
  100       JCN(K) = JCE
            IRN(K) = ICE
            A(K) = ACE
  110     CONTINUE
  120   CONTINUE
      END IF

  130 CONTINUE

      RETURN
      END
C
C**********************************************************
      SUBROUTINE MC59CD(NC,NR,NE,IRN,JCN,LA,A,IP,IW)
C
C   To sort a sparse matrix stored by rows,
C   unordered within each row, to ordering by columns, with
C   ordering by rows within each column.
C
C   NC - integer variable. Intent(IN)
C      - on entry must be set to the number of columns in the matrix
C   NR - integer variable. Intent(IN)
C      - on entry must be set to the number of rows in the matrix
C  NE - integer variable. Intent(IN)
C      - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(OUT).
C      - not set on entry.
C      - on exit,  IRN holds row indices with the row
C        indices for column 1 preceding those for column 2 and so on,
C        with ordering by rows within each column.
C  JCN - integer array of length NE. Intent(INOUT)
C      - on entry set to contain the column indices of the nonzeros
C        with indices for column 1 preceding those for column 2
C        and so on, with the order within columns is arbitrary.
C      - on exit, contents destroyed.
C  LA  - integer variable which defines the length of the array A.
C        Intent(IN)
C   A  - real (double precision/complex/complex*16) array of length LA
C        Intent(INOUT)
C      - if LA > 1, the array must be of length NE, and A(K)
C        must be set to the value of the entry in JCN(K);
C        on exit A, A(K) holds the value of the entry in IRN(K).
C      - if LA = 1, the array is not accessed
C  IP  - integer array of length NC+1. Intent(INOUT)
C      - not set on entry
C      - on exit, IP(J) contains the position in IRN (and A) of the
C        first entry in column J (J=1,...,NC)
C      - IP(NC+1) is set to NE+1
C  IW  - integer array of length NR+1.  Intent(IN)
C      - on entry, must be set on entry so that IW(J) points to the
C        position in JCN of the first entry in row J, J=1,...,NR, and
C        IW(NR+1) must be set to NE+1
C
C     .. Scalar Arguments ..
      INTEGER LA,NC,NE,NR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NR+1),JCN(NE)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,J1,J2,K,L,LOC,LOCP
C     ..
C     .. Executable Statements ..

C  Count the number of entries in each column

      DO 10 J = 1,NC
        IP(J) = 0
   10 CONTINUE

      IF (LA.GT.1) THEN

        DO 20 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
          IRN(K) = JCN(K)
   20   CONTINUE
        IP(NC+1) = NE + 1

C  Set IP so that IP(I) points to the first entry in column I+1

        IP(1) = IP(1) + 1
        DO 30 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
   30   CONTINUE

        DO 50 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 40 J = J1,J2
            K = IRN(J)
            L = IP(K) - 1
            JCN(J) = L
            IRN(J) = I
            IP(K) = L
   40     CONTINUE
   50   CONTINUE
        IP(NC+1) = NE + 1
        DO 70 J = 1,NE
          LOC = JCN(J)
          IF (LOC.EQ.0) GO TO 70
          ICE = IRN(J)
          ACE = A(J)
          JCN(J) = 0
          DO 60 K = 1,NE
            LOCP = JCN(LOC)
            ICEP = IRN(LOC)
            ACEP = A(LOC)
            JCN(LOC) = 0
            IRN(LOC) = ICE
            A(LOC) = ACE
            IF (LOCP.EQ.0) GO TO 70
            ICE = ICEP
            ACE = ACEP
            LOC = LOCP
   60     CONTINUE
   70   CONTINUE
      ELSE

C Pattern only

C  Count the number of entries in each column

        DO 90 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
   90   CONTINUE
        IP(NC+1) = NE + 1

C  Set IP so that IP(I) points to the first entry in column I+1

        IP(1) = IP(1) + 1
        DO 100 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
  100   CONTINUE

        DO 120 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 110 J = J1,J2
            K = JCN(J)
            L = IP(K) - 1
            IRN(L) = I
            IP(K) = L
  110     CONTINUE
  120   CONTINUE

      END IF

      RETURN
      END

C**********************************************************

      SUBROUTINE MC59DD(NC,NE,IRN,IP,LA,A)
C
C To sort from arbitrary order within each column to order
C by increasing row index. Note: this is taken from MC20B/BD.
C
C   NC - integer variable. Intent(IN)
C      - on entry must be set to the number of columns in the matrix
C   NE - integer variable. Intent(IN)
C      - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(INOUT)
C      - on entry set to contain the row indices of the nonzeros
C        ordered so that the row
C        indices for column 1 precede those for column 2 and so on,
C        but the order within columns is arbitrary.
C        On exit, the order within each column is by increasing
C        row indices.
C   LA - integer variable which defines the length of the array A.
C        Intent(IN)
C    A - real (double precision/complex/complex*16) array of length LA
C        Intent(INOUT)
C      - if LA > 1, the array must be of length NE, and A(K)
C        must be set to the value of the entry in IRN(K);
C        on exit A is reordered in the same way as IRN
C      - if LA = 1, the array is not accessed
C  IP  - integer array of length NC. Intent(IN)
C      - on entry, IP(J) contains the position in IRN (and A) of the
C        first entry in column J (J=1,...,NC)
C     . .
C     .. Scalar Arguments ..
      INTEGER LA,NC,NE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(NC)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ACE
      INTEGER ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Executable Statements ..

C Jump if pattern only.
      IF (LA.GT.1) THEN
        KMAX = NE
        DO 50 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 40
          KOR = KMAX
          DO 30 KDUMMY = KLO,KMAX
C Items KOR, KOR+1, .... ,KMAX are in order
            ACE = A(KOR-1)
            ICE = IRN(KOR-1)
            DO 10 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 20
              IRN(K-1) = IK
              A(K-1) = A(K)
   10       CONTINUE
            K = KMAX + 1
   20       IRN(K-1) = ICE
            A(K-1) = ACE
            KOR = KOR - 1
   30     CONTINUE
C Next column
   40     KMAX = KLO - 2
   50   CONTINUE
      ELSE

C Pattern only.
        KMAX = NE
        DO 150 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 140
          KOR = KMAX
          DO 130 KDUMMY = KLO,KMAX
C Items KOR, KOR+1, .... ,KMAX are in order
            ICE = IRN(KOR-1)
            DO 110 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 120
              IRN(K-1) = IK
  110       CONTINUE
            K = KMAX + 1
  120       IRN(K-1) = ICE
            KOR = KOR - 1
  130     CONTINUE
C Next column
  140     KMAX = KLO - 2
  150   CONTINUE
      END IF
      END
C***********************************************************************

      SUBROUTINE MC59ED(NC,NR,NE,IRN,LIP,IP,LA,A,IW,IDUP,KNE,ICNTL6)

C Checks IRN for duplicate entries.
C On exit, IDUP holds number of duplicates found and KNE is number
C of entries in matrix after removal of duplicates
C     . .
C     .. Scalar Arguments ..
      INTEGER ICNTL6,IDUP,KNE,LIP,LA,NC,NR,NE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(LIP),IW(NR)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,K,KSTART,KSTOP,NZJ

      IDUP = 0
      KNE = 0
C Initialise IW
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE

      KSTART = IP(1)
      IF (LA.GT.1) THEN
C Matrix entries considered
        NZJ = 0
        DO 30 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
C We have a duplicate in column J
              IDUP = IDUP + 1
C If requested, sum duplicates
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE

      ELSE

C Pattern only
        DO 50 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 40 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
C  We have a duplicate in column J
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF

      RETURN
      END
C***********************************************************************

      SUBROUTINE MC59FD(NC,NR,NE,IRN,LIP,IP,LA,A,LIW,IW,IDUP,IOUT,
     +                  IUP,KNE,ICNTL6,INFO)

C Checks IRN for duplicate and out-of-range entries.
C For symmetric matrix, also checks NO entries lie in upper triangle.
C Also checks IP is monotonic.
C On exit:
C IDUP holds number of duplicates found
C IOUT holds number of out-of-range entries
C For symmetric matrix, IUP holds number of entries in upper
C triangular part.
C KNE holds number of entries in matrix after removal of
C out-of-range and duplicate entries.
C Note: this is similar to MC59ED except it also checks IP is
C monotonic and removes out-of-range entries in IRN.
C     . .
C     .. Scalar Arguments ..
      INTEGER ICNTL6,IDUP,IOUT,IUP,KNE,LA,LIP,LIW,NC,NR,NE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(LIP),IW(LIW),INFO(2)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,K,KSTART,KSTOP,NZJ,LOWER

      IDUP = 0
      IOUT = 0
      IUP = 0
      KNE = 0
C Initialise IW
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE

      KSTART = IP(1)
      LOWER = 1
      IF (LA.GT.1) THEN
        NZJ = 0
        DO 30 J = 1,NC
C In symmetric case, entries out-of-range if they lie
C in upper triangular part.
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
C Check for out-of-range
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
C In symmetric case, check if entry is out-of-range because
C it lies in upper triangular part.
              IF (ICNTL6.NE.0 .AND. I.LT.J) IUP = IUP + 1
            ELSE IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
C  We have a duplicate in column J
              IDUP = IDUP + 1
C If requested, sum duplicates
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE

      ELSE

C Pattern only
        DO 50 J = 1,NC
C In symmetric case, entries out-of-range if lie
C in upper triangular part.
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO  40 K = KSTART,KSTOP - 1
            I = IRN(K)
C Check for out-of-range
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
              IF (ICNTL6.NE.0 .AND. I.GT.1) IUP = IUP + 1
            ELSE IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
C  We have a duplicate in column J
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF

      RETURN
      END

* COPYRIGHT (c) 1977 AEA Technology
* Original date 8 Oct 1992
C######8/10/92 Toolpack tool decs employed.
C######8/10/92 D version created by name change only.
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC21AD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUMNZ
C     ..
C     .. Array Arguments ..
      INTEGER ICN(LICN),IP(N),IPERM(N),IW(N,4),LENR(N)
C     ..
C     .. External Subroutines ..
      EXTERNAL MC21BD
C     ..
C     .. Executable Statements ..
      CALL MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),
     +            IW(1,3),IW(1,4))
      RETURN
C
      END
      SUBROUTINE MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
C   PR(I) IS THE PREVIOUS ROW TO I IN THE DEPTH FIRST SEARCH.
C IT IS USED AS A WORK ARRAY IN THE SORTING ALGORITHM.
C   ELEMENTS (IPERM(I),I) I=1, ... N  ARE NON-ZERO AT THE END OF THE
C ALGORITHM UNLESS N ASSIGNMENTS HAVE NOT BEEN MADE.  IN WHICH CASE
C (IPERM(I),I) WILL BE ZERO FOR N-NUMNZ ENTRIES.
C   CV(I) IS THE MOST RECENT ROW EXTENSION AT WHICH COLUMN I
C WAS VISITED.
C   ARP(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED WHEN LOOKING FOR A CHEAP ASSIGNMENT.
C   OUT(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED DURING ONE PASS THROUGH THE MAIN LOOP.
C
C   INITIALIZATION OF ARRAYS.
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUMNZ
C     ..
C     .. Array Arguments ..
      INTEGER ARP(N),CV(N),ICN(LICN),IP(N),IPERM(N),LENR(N),OUT(N),PR(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,II,IN1,IN2,IOUTK,J,J1,JORD,K,KK
C     ..
C     .. Executable Statements ..
      DO 10 I = 1,N
        ARP(I) = LENR(I) - 1
        CV(I) = 0
        IPERM(I) = 0
   10 CONTINUE
      NUMNZ = 0
C
C
C   MAIN LOOP.
C   EACH PASS ROUND THIS LOOP EITHER RESULTS IN A NEW ASSIGNMENT
C OR GIVES A ROW WITH NO ASSIGNMENT.
      DO 100 JORD = 1,N
        J = JORD
        PR(J) = -1
        DO 70 K = 1,JORD
C LOOK FOR A CHEAP ASSIGNMENT
          IN1 = ARP(J)
          IF (IN1.LT.0) GO TO 30
          IN2 = IP(J) + LENR(J) - 1
          IN1 = IN2 - IN1
          DO 20 II = IN1,IN2
            I = ICN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
C   NO CHEAP ASSIGNMENT IN ROW.
          ARP(J) = -1
C   BEGIN LOOKING FOR ASSIGNMENT CHAIN STARTING WITH ROW J.
   30     CONTINUE
          OUT(J) = LENR(J) - 1
C INNER LOOP.  EXTENDS CHAIN BY ONE OR BACKTRACKS.
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENR(J) - 1
            IN1 = IN2 - IN1
C FORWARD SCAN.
            DO 40 II = IN1,IN2
              I = ICN(II)
              IF (CV(I).EQ.JORD) GO TO 40
C   COLUMN I HAS NOT YET BEEN ACCESSED DURING THIS PASS.
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
C
   40       CONTINUE
C
C   BACKTRACKING STEP.
   50       CONTINUE
            J = PR(J)
            IF (J.EQ.-1) GO TO 100
   60     CONTINUE
C
   70   CONTINUE
C
C   NEW ASSIGNMENT IS MADE.
   80   CONTINUE
        IPERM(I) = J
        ARP(J) = IN2 - II - 1
        NUMNZ = NUMNZ + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 100
          II = IP(J) + LENR(J) - OUT(J) - 2
          I = ICN(II)
          IPERM(I) = J
   90   CONTINUE
C
  100 CONTINUE
C
C   IF MATRIX IS STRUCTURALLY SINGULAR, WE NOW COMPLETE THE
C PERMUTATION IPERM.
      IF (NUMNZ.EQ.N) RETURN
      DO 110 I = 1,N
        ARP(I) = 0
  110 CONTINUE
      K = 0
      DO 130 I = 1,N
        IF (IPERM(I).NE.0) GO TO 120
        K = K + 1
        OUT(K) = I
        GO TO 130
C
  120   CONTINUE
        J = IPERM(I)
        ARP(J) = I
  130 CONTINUE
      K = 0
      DO 140 I = 1,N
        IF (ARP(I).NE.0) GO TO 140
        K = K + 1
        IOUTK = OUT(K)
        IPERM(IOUTK) = I
  140 CONTINUE
      RETURN
C
      END
* COPYRIGHT (c) 1976 AEA Technology
* Original date 21 Jan 1993
C 8 August 2000: CONTINUEs given to DOs.
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC22AD(N,ICN,A,NZ,LENROW,IP,IQ,IW,IW1)
C     .. Scalar Arguments ..
      INTEGER N,NZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NZ)
      INTEGER ICN(NZ),IP(N),IQ(N),IW(N,2),IW1(NZ),LENROW(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AVAL
      INTEGER I,ICHAIN,IOLD,IPOS,J,J2,JJ,JNUM,JVAL,LENGTH,NEWPOS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS
C     ..
C     .. Executable Statements ..
      IF (NZ.LE.0) GO TO 1000
      IF (N.LE.0) GO TO 1000
C SET START OF ROW I IN IW(I,1) AND LENROW(I) IN IW(I,2)
      IW(1,1) = 1
      IW(1,2) = LENROW(1)
      DO 10 I = 2,N
        IW(I,1) = IW(I-1,1) + LENROW(I-1)
        IW(I,2) = LENROW(I)
   10 CONTINUE
C PERMUTE LENROW ACCORDING TO IP.  SET OFF-SETS FOR NEW POSITION
C     OF ROW IOLD IN IW(IOLD,1) AND PUT OLD ROW INDICES IN IW1 IN
C     POSITIONS CORRESPONDING TO THE NEW POSITION OF THIS ROW IN A/ICN.
      JJ = 1
      DO 20 I = 1,N
        IOLD = IP(I)
        IOLD = IABS(IOLD)
        LENGTH = IW(IOLD,2)
        LENROW(I) = LENGTH
        IF (LENGTH.EQ.0) GO TO 20
        IW(IOLD,1) = IW(IOLD,1) - JJ
        J2 = JJ + LENGTH - 1
        DO 15 J = JJ,J2
          IW1(J) = IOLD
   15   CONTINUE
        JJ = J2 + 1
   20 CONTINUE
C SET INVERSE PERMUTATION TO IQ IN IW(.,2).
      DO 30 I = 1,N
        IOLD = IQ(I)
        IOLD = IABS(IOLD)
        IW(IOLD,2) = I
   30 CONTINUE
C PERMUTE A AND ICN IN PLACE, CHANGING TO NEW COLUMN NUMBERS.
C
C ***   MAIN LOOP   ***
C EACH PASS THROUGH THIS LOOP PLACES A CLOSED CHAIN OF COLUMN INDICES
C     IN THEIR NEW (AND FINAL) POSITIONS ... THIS IS RECORDED BY
C     SETTING THE IW1 ENTRY TO ZERO SO THAT ANY WHICH ARE SUBSEQUENTLY
C     ENCOUNTERED DURING THIS MAJOR SCAN CAN BE BYPASSED.
      DO 200 I = 1,NZ
        IOLD = IW1(I)
        IF (IOLD.EQ.0) GO TO 200
        IPOS = I
        JVAL = ICN(I)
C IF ROW IOLD IS IN SAME POSITIONS AFTER PERMUTATION GO TO 150.
        IF (IW(IOLD,1).EQ.0) GO TO 150
        AVAL = A(I)
C **  CHAIN LOOP  **
C EACH PASS THROUGH THIS LOOP PLACES ONE (PERMUTED) COLUMN INDEX
C     IN ITS FINAL POSITION  .. VIZ. IPOS.
        DO 100 ICHAIN = 1,NZ
C NEWPOS IS THE ORIGINAL POSITION IN A/ICN OF THE ELEMENT TO BE PLACED
C IN POSITION IPOS.  IT IS ALSO THE POSITION OF THE NEXT ELEMENT IN
C     THE CHAIN.
          NEWPOS = IPOS + IW(IOLD,1)
C IS CHAIN COMPLETE ?
          IF (NEWPOS.EQ.I) GO TO 130
          A(IPOS) = A(NEWPOS)
          JNUM = ICN(NEWPOS)
          ICN(IPOS) = IW(JNUM,2)
          IPOS = NEWPOS
          IOLD = IW1(IPOS)
          IW1(IPOS) = 0
C **  END OF CHAIN LOOP  **
  100   CONTINUE
  130   A(IPOS) = AVAL
  150   ICN(IPOS) = IW(JVAL,2)
C ***   END OF MAIN LOOP   ***
  200 CONTINUE
C
 1000 RETURN

      END
* COPYRIGHT (c) 1987 AEA Technology
* Original date 10 Feb 1993
C       Toolpack tool decs employed.
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC34AD(N,IRN,JCOLST,YESA,A,IW)
C THIS SUBROUTINE ACCEPTS AS INPUT THE STANDARD DATA STRUCTURE FOR
C     A SYMMETRIC MATRIX STORED AS A LOWER TRIANGLE AND PRODUCES
C     AS OUTPUT THE SYMMETRIC MATRIX HELD IN THE SAME DATA
C     STRUCTURE AS A GENERAL MATRIX.
C N IS AN INTEGER VARIABLE THAT MUST BE SET BY THE USER TO THE
C     ORDER OF THE MATRIX. NOT ALTERED BY THE ROUTINE
C     RESTRICTION (IBM VERSION ONLY): N LE 32767.
C IRN IS AN INTEGER (INTEGER*2 IN IBM VERSION) ARRAY THAT
C     MUST BE SET BY THE USER TO HOLD THE ROW INDICES OF THE LOWER
C     TRIANGULAR PART OF THE SYMMETRIC MATRIX.  THE ENTRIES OF A
C     SINGLE COLUMN ARE CONTIGUOUS. THE ENTRIES OF COLUMN J
C     PRECEDE THOSE OF COLUMN J+1 (J_=_1, ..., N-1), AND THERE IS
C     NO WASTED SPACE BETWEEN COLUMNS. ROW INDICES WITHIN A COLUMN
C     MAY BE IN ANY ORDER.  ON EXIT IT WILL HAVE THE SAME MEANING
C     BUT WILL BE CHANGED TO HOLD THE ROW INDICES OF ENTRIES IN
C     THE EXPANDED STRUCTURE.  DIAGONAL ENTRIES NEED NOT BE
C     PRESENT. THE NEW ROW INDICES ADDED IN THE UPPER TRIANGULAR
C     PART WILL BE IN ORDER FOR EACH COLUMN AND WILL PRECEDE THE
C     ROW INDICES FOR THE LOWER TRIANGULAR PART WHICH WILL REMAIN
C     IN THE INPUT ORDER.
C JCOLST IS AN INTEGER ARRAY OF LENGTH N+1 THAT MUST BE SET BY
C     THE USER SO THAT JCOLST(J) IS THE POSITION IN ARRAYS IRN AND
C     A OF THE FIRST ENTRY IN COLUMN J (J_=_1, ..., N).
C     JCOLST(N+1) MUST BE SET TO ONE MORE THAN THE TOTAL NUMBER OF
C     ENTRIES.  ON EXIT, JCOLST(J) WILL HAVE THE SAME MEANING BUT
C     WILL BE CHANGED TO POINT TO THE POSITION OF THE FIRST ENTRY
C     OF COLUMN J IN THE EXPANDED STRUCTURE. THE NEW VALUE OF
C     JCOLST(N+1) WILL BE ONE GREATER THAN THE NUMBER OF ENTRIES
C     IN THE EXPANDED STRUCTURE.
C YESA IS A LOGICAL VARIABLE THAT MUST BE SET TO .TRUE. IF THE
C     USER DESIRES TO GENERATE THE EXPANDED FORM FOR THE VALUES ALSO.
C     IF YESA IS .FALSE., THE ARRAY A WILL NOT BE REFERENCED.  IT IS
C     NOT ALTERED BY THE ROUTINE.
C A IS A REAL (DOUBLE PRECISION IN THE D VERSION) ARRAY THAT
C     CAN BE SET BY THE USER SO THAT A(K) HOLDS THE VALUE OF THE
C     ENTRY IN POSITION K OF IRN, {K = 1, _..._ JCOLST(N+1)-1}.
C     ON EXIT, IF YESA IS .TRUE., THE ARRAY WILL HOLD THE VALUES
C     OF THE ENTRIES IN THE EXPANDED STRUCTURE CORRESPONDING TO
C     THE OUTPUT VALUES OF IRN.   IF YESA IS .FALSE., THE ARRAY IS
C     NOT ACCESSED BY THE SUBROUTINE.
C IW IS AN INTEGER (INTEGER*2 IN IBM VERSION) ARRAY OF LENGTH
C     N THAT WILL BE USED AS WORKSPACE.
C
C CKP1 IS A LOCAL VARIABLE USED AS A RUNNING POINTER.
C OLDTAU IS NUMBER OF ENTRIES IN SYMMETRIC STORAGE.
C     .. Scalar Arguments ..
      INTEGER N
      LOGICAL YESA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*)
      INTEGER IRN(*),IW(*),JCOLST(*)
C     ..
C     .. Local Scalars ..
      INTEGER CKP1,I,I1,I2,II,IPKP1,IPOS,J,JSTART,LENK,NDIAG,NEWTAU,
     +        OLDTAU
C     ..
C     .. Executable Statements ..
C
      OLDTAU = JCOLST(N+1) - 1
C INITIALIZE WORK ARRAY
      DO 5 I = 1,N
        IW(I) = 0
    5 CONTINUE
C
C IW(J) IS SET EQUAL TO THE TOTAL NUMBER OF ENTRIES IN COLUMN J
C     OF THE EXPANDED SYMMETRIC MATRIX.
C NDIAG COUNTS NUMBER OF DIAGONAL ENTRIES PRESENT
      NDIAG = 0
      DO 20 J = 1,N
        I1 = JCOLST(J)
        I2 = JCOLST(J+1) - 1
        IW(J) = IW(J) + I2 - I1 + 1
        DO 10 II = I1,I2
          I = IRN(II)
          IF (I.NE.J) THEN
            IW(I) = IW(I) + 1

          ELSE
            NDIAG = NDIAG + 1
          END IF

   10   CONTINUE
   20 CONTINUE
C
C NEWTAU IS NUMBER OF ENTRIES IN EXPANDED STORAGE.
      NEWTAU = 2*OLDTAU - NDIAG
C IPKP1 POINTS TO POSITION AFTER END OF COLUMN BEING CURRENTLY
C     PROCESSED
      IPKP1 = OLDTAU + 1
C CKP1 POINTS TO POSITION AFTER END OF SAME COLUMN IN EXPANDED
C     STRUCTURE
      CKP1 = NEWTAU + 1
C GO THROUGH THE ARRAY IN THE REVERSE ORDER PLACING LOWER TRIANGULAR
C     ELEMENTS IN THE APPROPRIATE SLOTS.
      DO 40 J = N,1,-1
        I1 = JCOLST(J)
        I2 = IPKP1
C LENK IS NUMBER OF ENTRIES IN COLUMN J OF ORIGINAL STRUCTURE
        LENK = I2 - I1
C JSTART IS RUNNING POINTER TO POSITION IN NEW STRUCTURE
        JSTART = CKP1
C SET IKP1 FOR NEXT COLUMN
        IPKP1 = I1
        I2 = I2 - 1
C RUN THROUGH COLUMNS IN REVERSE ORDER
C LOWER TRIANGULAR PART OF COLUMN MOVED TO END OF SAME COLUMN IN
C     EXPANDED FORM
        DO 30 II = I2,I1,-1
          JSTART = JSTART - 1
          IF (YESA) A(JSTART) = A(II)
          IRN(JSTART) = IRN(II)
   30   CONTINUE
C JCOLST IS SET TO POSITION OF FIRST ENTRY IN LOWER TRIANGULAR PART OF
C     COLUMN J IN EXPANDED FORM
        JCOLST(J) = JSTART
C SET CKP1 FOR NEXT COLUMN
        CKP1 = CKP1 - IW(J)
C RESET IW(J) TO NUMBER OF ENTRIES IN LOWER TRIANGLE OF COLUMN.
        IW(J) = LENK
   40 CONTINUE
C
C AGAIN SWEEP THROUGH THE COLUMNS IN THE REVERSE ORDER, THIS
C     TIME WHEN ONE IS HANDLING COLUMN J THE UPPER TRIANGULAR
C     ELEMENTS A(J,I) ARE PUT IN POSITION.
      DO 80 J = N,1,-1
        I1 = JCOLST(J)
        I2 = JCOLST(J) + IW(J) - 1
C RUN DOWN COLUMN IN ORDER
C NOTE THAT I IS ALWAYS GREATER THAN OR EQUAL TO J
        DO 60 II = I1,I2
          I = IRN(II)
          IF (I.EQ.J) GO TO 60
          JCOLST(I) = JCOLST(I) - 1
          IPOS = JCOLST(I)
          IF (YESA) A(IPOS) = A(II)
          IRN(IPOS) = J
   60   CONTINUE
   80 CONTINUE
      JCOLST(N+1) = NEWTAU + 1
      RETURN

      END
C COPYRIGHT (c) 1995 Timothy A. Davis, Patrick Amestoy and
C           Council for the Central Laboratory of the Research Councils
C Original date 30 November 1995
C  April 2001: call to MC49 changed to MC59 to make routine threadsafe
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.
C 23 May 2007 Version 1.1.0. Absolute value of hash taken to cover the
C            case of integer overflow.
C            Comments with character in column 2 corrected.
C 2 August 2007 Version 2.0.0 Dense row handling added, error & warning
C            messages added, iovflo added, interface changed. MC47I/ID
C            added.
C 31 October 2007 Version 2.1.0 Corrected tree formation when handling
C            full variables
C

      SUBROUTINE MC47ID(ICNTL)
      INTEGER ICNTL(10)
C ICNTL is an INTEGER array of length 10 that contains control
C     parameters and must be set by the user. Default values are set
C     by MA57ID.
C
C     ICNTL(1) is the stream number for error messages. Printing is
C     suppressed if ICNTL(1)<0. The default is 6.
C
C     ICNTL(2) is the stream number for warning messages. Printing is
C     suppressed if ICNTL(2)<0. The default is 6.
C
C     ICNTL(3) is the stream number for printing matrix data on entry
C     to and exit from MC47A/AD. Printing is suppressed if ICNTL(3)<0.
C     The default is -1.
C
C     ICNTL(4) controls the choice of AMD algorithm
C   =-1 Classical MC47B (AMD) algorithm (no dense row detection)
C   = 0 Only exactly dense rows in the reduced matrix are selected.
C   = 1 Corresponds to automatic setting of the minimum density
C     requirement.
C     The default value is 1.
C
C     ICNTL(5) defines the largest positive
C     integer that your computer can represent (-iovflo should also
C     be representable). HUGE(1) in Fortran 95. The default value is
C     2139062143
C
C     ICNTL(6) to ICNTL(10) are set to zero by MA57ID but are not
C     currently used by MC47.
C
C Local variables
      INTEGER I

      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = -1
      ICNTL(4) = 1
      ICNTL(5) = 2139062143

      DO 100 I=6,10
        ICNTL(I) = 0
 100  CONTINUE
      RETURN
      END


      SUBROUTINE MC47AD(N, NE, PE, IW, IWLEN,
     *      ICNTL,INFO, RINFO)
      INTEGER N, NE, PE(N+1), IWLEN, IW(IWLEN), INFO(10)
      INTEGER ICNTL(10)
      DOUBLE PRECISION RINFO(10)
C N is an INTEGER variable that must be set by the user to the
C     order of the matrix A.  It is not altered by the subroutine.
C     Restriction: N >= 1.
C NE is an INTEGER variable that must be set by the user to the
C     number of entries in the matrix A.
C     It is not altered by the subroutine.
C PE is an INTEGER  array of size N+1 that must be set by the user.
C     If the user is supplying the entries by columns, then PE(i) must
C     hold the index in IW of the start of column i, i=1, ..., N and
C     PE(N+1) must be equal to NE+1.  If the user is supplying row and
C     column indices for the matrix, PE(1) must be negative.
C     On exit, PE will hold information on the matrix factors.
C IW is an INTEGER  array of length IWLEN that must be set by the user
C     to hold the pattern of the matrix A.  If PE(1) is positive,
C     IW(PE(J)), ..., IW(PE(J+1)-1) must hold the row indices of entries
C     in column J, J = 1, ..., N.
C     The entries within a column need not be in order.
C     If PE(1) is negative, then (IW(k), IW(NE+k)), k = 1, ..., NE, must
C     hold the row and column index of an entry.
C     Duplicates, out-of-range entries, diagonal entries, and entries
C     in upper triangle are ignored.
C     IW is used as workspace by the subroutine.  On exit, the
C     permutation generated by MC47 is held in IW(IWLEN-I+1), I=1, ... N
C     and is such that the kth column of the permuted matrix is column
C     IW(IWLEN-N+k) of the original matrix.
C     The inverse permutation is held in positions IWLEN-2N+1 to
C     IWLEN-N of IW, preceded by further information on the structure
C     of the factors
C IWLEN is an INTEGER variable. It must be set by the user to the length
C     of array IW and is not altered by the subroutine.
C     We recommend IWLEN > 2NE + 9N.   Restriction: IWLEN >= 2NE + 8N.
C ICNTL is an INTEGER array of length 10 that contains control
C     parameters and must be set by the user. Default values for the
C     components may be set by a call to MC47ID.
C INFO is an INTEGER array of length 8 that need not be set by the user.
C     On return from MC47A/AD, a value of zero for INFO(1) indicates
C     that the subroutine has performed successfully.  Negative values
C     for INFO(1) signify a fatal error.  Possible values are:
C     -1    N < 1
C     -2    IWLEN < 2NE + 8N.
C     -3    Error in PE when input is by columns.
C     -4    Matrix is null (usually because all entries in upper
C           triangle.
C     There is one warning indicated by a positive value for INFO(1)
C     +1    Out-of-range index, duplicate, diagonal or entry in upper
C           triangle in input. Action taken is to ignore these entries.
C     The other entries of INFO give information to the user.
C     INFO(2) gives the number of compresses performed on the array IW.
C             A large value for this indicates that the ordering could
C             be found more quickly if IWLEN were increased.
C     INFO(3) gives the minimum necessary value for IWLEN for a
C             successful run of MC47A/AD on the same matrix as has just
C             been analysed, without the need for any compresses of IW.
C     INFO(4) gives the number of entries with row or column indices
C             that are out of range.  Any such entry is ignored.
C     INFO(5) gives the number of duplicate entries.  Any such entry is
C             ignored.
C     INFO(6) gives the number of entries in upper triangle.  Any such
C             entry is ignored.
C     INFO(7) gives the number of diagonal entries.  Any such entry is
C             ignored.
C     INFO(8) gives the number of restarts performed.
C RINFO is an INTEGER array of length 10 that need not be set by the
C     user. The other entries of RINFO give information to the user.
C     RINFO(1) gives forecast number of reals to hold the factorization
C     RINFO(2) gives the forecast number of flops required by the
C     factorization if no pivoting is performed.

C Local variables
      INTEGER DEGREE
      DOUBLE PRECISION DUMMY(1)
      INTEGER ELEN,HEAD,I,II,I1,I2,J,LAST,LEN,LENIW,LP,MP,
     *        NEXT,NV,PFREE,W,WP
      INTEGER ICT59(10),INFO59(10),IOUT,JOUT,IDUP,JNFO(10)
C DEGREE is used to subdivide array IW
C DUMMY  is dummy real for call to MC34A/AD
C ELEN   is used to subdivide array IW
C HEAD   is used to subdivide array IW
C I      is DO loop variable
C IFLAG  is eror return from MC59A/AD
C II     is running index for entries of IW in column
C I1     is first location for entries of IW in column
C I2     is flast location for entries of IW in column
C J      is column index and DO loop variable
C LAST   is used to subdivide array IW
C LEN    is used to subdivide array IW
C LENIW  is space left in IW after allocation of work vectors
C LP     is a local copy of ICNTL(1)
C MP     is a local copy of ICNTL(3)
C NEXT   is used to subdivide array IW
C NV     is used to subdivide array IW
C PFREE  marks active length of IW for call to MC47B/BD
C W      is used to subdivide array IW
C WP     is a local copy of ICNTL(2)
C IOUT   number of indices out-of-range
C JOUT   number of column indices out-of-range (as detected by MC59A/AD)
C IDUP   number of duplicates


C Subroutines called
      EXTERNAL MC59AD,MC34AD,MC47BD

C Initialize info
      DO 5 J = 1,10
         INFO(J) = 0
 5    CONTINUE

C Set LP,WP,MP
      LP = ICNTL(1)
      WP = ICNTL(2)
      MP = ICNTL(3)

C Check value of N
      IF (N.LT.1) THEN
        INFO(1) = -1
        IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10)')
     +       '**** Error return from MC47AD **** INFO(1) =',INFO(1),
     +       'N has value ',N
        GO TO 1000
      ENDIF

C Error check to see if enough space in IW to get started
      IF (PE(1).LT.1) THEN
        IF (2*NE+N.GT.IWLEN) THEN
          INFO(1) = -2
         IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10)')
     +        '**** Error return from MC47AD **** INFO(1) =',INFO(1),
     +        'IWLEN has value ',IWLEN
          GO TO 1000
        ENDIF
      ELSE
        IF (NE+N.GT.IWLEN) THEN
          INFO(1) = -2
         IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10)')
     +        '**** Error return from MC47AD **** INFO(1) =',INFO(1),
     +        'IWLEN has value ',IWLEN
          GO TO 1000
        ENDIF
      ENDIF

C Diagnostic print
      IF (MP.GE.0) THEN
        WRITE(MP,'(/A)') 'Entry to MC47A/AD'
        WRITE(MP,'(A,I10,A,I10,A)') 'Matrix of order',N,' with',NE,
     *                            ' entries'
        IF (PE(1).LT.0)  THEN
          WRITE(MP,'(A)') 'Matrix input in coordinate form'
          WRITE(MP,'(A/(4(I8,I8)))') 'Row and column indices',
     *          (IW(I),IW(NE+I),I=1,NE)
        ELSE
          WRITE(MP,'(A)') 'Matrix input by columns'
          DO 10 J=1,N
            WRITE(MP,'(A,I4/(10I8))') 'Column',J,
     *                                (IW(I),I=PE(J),PE(J+1)-1)
   10     CONTINUE
        ENDIF
      ENDIF

C Divide workspace
      LAST   = IWLEN  - N + 1
      ELEN   = LAST   - N
      NV     = ELEN   - N
      W      = NV     - N
      DEGREE = W      - N
      HEAD   = DEGREE - N
      NEXT   = HEAD   - N
      LEN    = NEXT   - N
      LENIW = LEN-1

C Set counters for number of upper triangular entries and diagonals
C     present in input matrix.
C These will be removed for later processing.
      INFO(6) = 0
      INFO(7) = 0

      IF (PE(1).LT.0) THEN
C First remove diagonals (if present) and all entries in upper triangle
C Note that this may give out-of-range entries from MC59.
        DO 20 I=1,NE
          IF (IW(I).LE.IW(NE+I)) THEN
            IF (IW(I).EQ.IW(NE+I) .AND. IW(I).NE.0) THEN
              INFO(7) = INFO(7) + 1
            ELSE
              IF (IW(I).GT.0) INFO(6) = INFO(6) + 1
            ENDIF
            IW(I)=0
          ENDIF
   20   CONTINUE

C Call sort routine
        ICT59(1) = 0
        ICT59(2) = 1
        ICT59(3) = 1
        ICT59(4) = LP
        ICT59(5) = -1
        ICT59(6) = 0
        CALL MC59AD(ICT59,N,N,NE,IW,NE,IW(NE+1),1,DUMMY,
     *              N+1,PE,N+1,IW(2*NE+1),INFO59)
C        IFLAG = INFO59(1)
        IDUP  = INFO59(3)
        IOUT  = INFO59(4)
        JOUT  = INFO59(5)
      ELSE

C Matrix already sorted by columns.
C Remove duplicates, out-of-range indices, and entries in upper
C       triangle.  First initialize counts.
        IDUP = 0
        IOUT = 0
        JOUT = 0

C Set array used to find duplicates
        DO 30 I = 1,N
          IW(NE+I) = 0
   30   CONTINUE

        DO 50 J=1,N
          I1 = PE(J)
          PE(J) = I1-(IOUT+IDUP)
          I2 = PE(J+1)-1
          IF (I2.LT.I1-1) THEN
            INFO(1) = -3
            GO TO 1000
          ENDIF
          DO 40 II = I1,I2
            I = IW(II)
            IF (I.LE.J .OR. I.GT.N) THEN
              IF (I.EQ.J) INFO(7) = INFO(7) + 1
              IF (I.GT.0 .AND. I.LT.J) INFO(6) = INFO(6) + 1
              IOUT = IOUT + 1
            ELSE
              IF (IW(NE+I).EQ.J) THEN
C Duplicate found
                IDUP = IDUP + 1
              ELSE
                IW(NE+I)=J
                IW(II-(IOUT+IDUP)) = I
              ENDIF
            ENDIF
   40     CONTINUE
   50   CONTINUE
        PE(N+1) = NE - (IOUT+IDUP) + 1
      ENDIF

C Set flags for duplicates or out-of-range entries
C Check if there were duplicates
      IF (IDUP.GT.0) THEN
        INFO(1) = 1
        INFO(4) = IDUP
        IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
     +       '**** Warning from MC47AD **** INFO(1) =',INFO(1),
     +       'Number of duplicates found: ',INFO(4)
      ELSE
        INFO(4) = 0
      ENDIF
C Check for out of range entries
      IF (IOUT+ JOUT - INFO(7) .GT.0 ) THEN
        INFO(1) = 1
        INFO(5) = IOUT + JOUT - INFO(7)
        IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
     +       '**** Warning from MC47AD **** INFO(1) =',INFO(1),
     +       'Number of out of range entries found and ignored: ',
     +       INFO(5)
      ELSE
        INFO(5) = 0
      ENDIF

C Check for entries in upper triangle
      IF (INFO(6).GT.0) THEN
         INFO(1) = 1
        IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
     +       '**** Warning from MC47AD **** INFO(1) =',INFO(1),
     +       'Number of entries in upper triangle found and ignored: ',
     +        INFO(6)
      ENDIF

C Check for entries in diagonals
      IF (INFO(7).GT.0) THEN
         INFO(1) = 1
        IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
     +       '**** Warning from MC47AD **** INFO(1) =',INFO(1),
     +       'Number of entries in diagonals found and ignored: ',
     +        INFO(7)
      ENDIF

C Check for null matrix
C Usually happens if wrong triangle is input
      IF (NE-(IOUT+IDUP).EQ.0) THEN
        INFO(1) = -4
        IF (LP.GE.0) WRITE(LP,'(/A,I3/A)')
     +       '**** Error return from MC47AD **** INFO(1) =',INFO(1),
     +       'Matrix is null'
        GO TO 1000
      ENDIF

C Generate matrix in expanded form
C First check there is sufficient space in IW
      IF (LENIW.LT.2*(PE(N+1)-1)) THEN
        INFO(1) = -2
         IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10/A,I10)')
     +        '**** Error return from MC47AD **** INFO(1) =',INFO(1),
     +        'IWLEN has value ',IWLEN,
     +        'Should be at least', 2*(PE(N+1)-1)+8*N
        GO TO 1000
      ENDIF

      CALL MC34AD(N,IW,PE,.FALSE.,DUMMY,IW(W))
      PFREE = PE(N+1)

C Set length array for MC47B/BD
      DO 60 I=1,N
        IW(LEN+I-1) = PE(I+1) - PE(I)
   60 CONTINUE

C Call to approximate minimum degree subroutine.
      CALL MC47BD(N,LENIW,PE,PFREE,IW(LEN),IW,IW(NV),
     *            IW(ELEN),IW(LAST),IW(DEGREE),
     *            IW(HEAD),IW(NEXT),IW(W), ICNTL,JNFO, RINFO)

      INFO(2) = JNFO(1)
      INFO(3) = PFREE+8*N
      INFO(8) = JNFO(2)

C Print diagnostics
      IF (MP.GE.0) THEN
        WRITE(MP,'(/A)') 'Exit from MC47A/AD'
        WRITE(MP,'(A/(7I10))') 'INFO(1-10):',(INFO(I),I=1,10)
        WRITE(MP,'(A/(8I10))') 'Parent array',(PE(I),I=1,N)
        WRITE(MP,'(A/(8I10))') 'Permutation',(IW(ELEN+I-1),I=1,N)
        WRITE(MP,'(A/(8I10))') 'Inverse permutation',
     *                         (IW(LAST+I-1),I=1,N)
        WRITE(MP,'(A/(8I10))') 'Degree array',(IW(NV+I-1),I=1,N)
      ENDIF

 1000 RETURN
      END

C ====================================================================
C ====================================================================
C ====================================================================

      SUBROUTINE MC47BD (N, IWLEN, PE, PFREE, LEN, IW, NV,
     $                   ELEN, LAST, DEGREE,
     $                   HEAD, DENXT, W, ICNTL, JNFO, RJNFO)

      INTEGER N, IWLEN, PE(N), PFREE, LEN(N), IW(IWLEN), NV(N),
     $        ELEN(N), LAST(N),  DEGREE(N),
     $         HEAD(N), DENXT(N), W(N), ICNTL(10), JNFO(10)

      DOUBLE PRECISION RJNFO(10)

C -------------------------------------------------------------------
C AMDD is a modified version of
C     MC47B:  Approximate Minimum (UMFPACK/MA38-style, external) Degree
C          ordering algorithm, with aggresive absorption
C designed to automatically detect and exploit dense
C rows in the reduced matrix at any step of the minimum degree.
C
C We use the term Le to denote the set of all supervariables in element
C E.
C **** Reword below***
C A row is declared as full if none of its entries can be guaranteed
C     to be zero.
C A row is quasi dense if at most N-THRESM-1 of its entries can be
C     guaranteed to be zero at it has not been recognized as being full
C     in the calculation so far.
C A row is dense if it is either full or quasi dense.
C A row is sparse if it is not dense.
C -------------------------------------------------------------------
C
C
C N must be set to the matrix order. It is not altered.
C     Restriction:  N .ge. 1
C
C IWLEN must be set to the length of IW. It is not altered. On input,
C     the matrix is stored in IW (1..PFREE-1).
C     *** We do not recommend running this algorithm with ***
C     ***      IWLEN .LT. PFREE + N.                      ***
C     *** Better performance will be obtained if          ***
C     ***      IWLEN .GE. PFREE + N                       ***
C     *** or better yet                                   ***
C     ***      IWLEN .GT. 1.2 * PFREE                     ***
C     Restriction: IWLEN .GE. PFREE-1
C
C PE(i) must be set to the the index in IW of the start of row I, or be
C     zero if row I has no off-diagonal entries. During execution,
C     it is used for both supervariables and elements:
C       * Principal supervariable I:  index into IW of the
C               list of supervariable I.  A supervariable
C               represents one or more rows of the matrix
C               with identical pattern.
C       * Non-principal supervariable I:  if I has been absorbed
C               into another supervariable J, then PE(I) = -J.
C               That is, J has the same pattern as I.
C               Note that J might later be absorbed into another
C               supervariable J2, in which case PE(I) is still -J,
C               and PE(J) = -J2.
C       * Unabsorbed element E:  the index into IW of the list
C               of element E.  Element E is created when
C               the supervariable of the same name is selected as
C               the pivot.
C       * Absorbed element E:  if element E is absorbed into element
C               E2, then PE(E) = -E2.  This occurs when one of its
C               variables is eliminated and when the pattern of
C               E (that is, Le) is found to be a subset of the pattern
C               of E2 (that is, Le2).  If element E is "null" (it has
C               no entries outside its pivot block), then PE(E) = 0.
C
C     On output, PE holds the assembly tree/forest, which implicitly
C     represents a pivot order with identical fill-in as the actual
C     order (via a depth-first search of the tree). If NV(I) .GT. 0,
C     then I represents a node in the assembly tree, and the parent of
C     I is -PE(I), or zero if I is a root. If NV(I)=0, then (I,-PE(I))
C     represents an edge in a subtree, the root of which is a node in
C     the assembly tree.
C
C PFREE must be set to the position in IW of the first free variable.
C     During execution, additional data is placed in IW, and PFREE is
C     modified so that components  of IW from PFREE are free.
C     On output, PFREE is set equal to the size of IW that would have
C     caused no compressions to occur.  If NCMPA is zero, then
C     PFREE (on output) is less than or equal to IWLEN, and the space
C     IW(PFREE+1 ... IWLEN) was not used. Otherwise, PFREE (on output)
C     is greater than IWLEN, and all the memory in IW was used.
C
C LEN(I) must be set to hold the number of entries in row I of the
C     matrix, excluding the diagonal.  The contents of LEN(1..N) are
C     undefined on output.
C
C IW(1..PFREE-1) must be set to  hold the patterns of the rows of
C     the matrix.  The matrix must be symmetric, and both upper and
C     lower triangular parts must be present.  The diagonal must not be
C     present.  Row I is held as follows:
C               IW(PE(I)...PE(I) + LEN(I) - 1) must hold the list of
C               column indices for entries in row I (simple
C               supervariables), excluding the diagonal.  All
C               supervariables start with one row/column each
C               (supervariable I is just row I). If LEN(I) is zero on
C               input, then PE(I) is ignored on input. Note that the
C               rows need not be in any particular order, and there may
C               be empty space between the rows.
C     During execution, the supervariable I experiences fill-in. This
C     is represented by constructing a list of the elements that cause
C     fill-in in supervariable I:
C               IE(PE(i)...PE(I) + ELEN(I) - 1) is the list of elements
C               that contain I. This list is kept short by removing
C               absorbed elements. IW(PE(I)+ELEN(I)...PE(I)+LEN(I)-1)
C               is the list of supervariables in I. This list is kept
C               short by removing nonprincipal variables, and any entry
C               J that is also contained in at least one of the
C               elements in the list for I.
C     When supervariable I is selected as pivot, we create an element E
C     of the same name (E=I):
C               IE(PE(E)..PE(E)+LEN(E)-1) is the list of supervariables
C                in element E.
C     An element represents the fill-in that occurs when supervariable
C     I is selected as pivot.
C     CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
C     The contents of IW are undefined on output.
C
C NV(I) need not be set. During execution, ABS(NV(I)) is equal to the
C     number of rows represented by the principal supervariable I. If I
C     is a nonprincipal variable, then NV(I) = 0. Initially, NV(I) = 1
C     for all I.  NV(I) .LT. 0 signifies that I is a principal variable
C     in the pattern Lme of the current pivot element ME. On output,
C     NV(E) holds the true degree of element E at the time it was
C     created (including the diagonal part).
C
C ELEN(I) need not be set. See the description of IW above. At the
C     start of execution, ELEN(I) is set to zero. For a supervariable,
C     ELEN(I) is the number of elements in the list for supervariable
C     I. For an element, ELEN(E) is the negation of the position in the
C     pivot sequence of the supervariable that generated it. ELEN(I)=0
C     if I is nonprincipal.
C     On output ELEN(1..N) holds the inverse permutation (the same
C     as the 'INVP' argument in Sparspak). That is, if K = ELEN(I),
C     then row I is the Kth pivot row.  Row I of A appears as the
C     (ELEN(I))-th row in the permuted matrix, PAP^T.
C
C LAST(I) need not be set on input. In a degree list, LAST(I) is the
C     supervariable preceding I, or zero if I is the head of the list.
C     In a hash bucket, LAST(I) is the hash key for I. LAST(HEAD(HASH))
C     is also used as the head of a hash bucket if HEAD(HASH) contains
C     a degree list (see HEAD, below).
C     On output, LAST(1..N) holds the permutation (the same as the
C     'PERM' argument in Sparspak). That is, if I = LAST(K), then row I
C     is the Kth pivot row.  Row LAST(K) of A is the K-th row in the
C     permuted matrix, PAP^T.
C
C
C DEGREE need not be set on input. If I is a supervariable and sparse,
C     then DEGREE(I) holds the current approximation of the external
C     degree of row I (an upper bound). The external degree is the
C     number of entries in row I, minus ABS(NV(I)) (the diagonal
C     part). The bound is equal to the external degree if ELEN(I) is
C     less than or equal to two. We also use the term "external degree"
C     for elements E to refer to |Le \ Lme|. If I is full in the reduced
C     matrix, then DEGREE(I)=N+1. If I is dense in the reduced matrix,
C     then DEGREE(I)=N+1+last_approximate_external_deg of I.
C     All dense rows are stored in the list pointed by HEAD(N).
C     Quasi dense rows are stored first, and are followed by full rows
C     in the reduced matrix. LASTD holds the last row in
C     this list of dense rows or is zero if the list is empty.
C
C HEAD(DEG) need not be set on input. HEAD is used for degree lists.
C     HEAD(DEG) is the first supervariable in a degree list (all
C     supervariables I in a degree list DEG have the same approximate
C     degree, namely, DEG = DEGREE(I)). If the list DEG is empty then
C     HEAD(DEG) = 0.
C     During supervariable detection HEAD(HASH) also serves as a
C     pointer to a hash bucket.
C     If HEAD(HASH) .GT. 0, there is a degree list of degree HASH. The
C     hash bucket head pointer is LAST(HEAD(HASH)).
C     If HEAD(HASH) = 0, then the degree list and hash bucket are
C     both empty.
C     If HEAD(HASH) .LT. 0, then the degree list is empty, and
C     -HEAD(HASH) is the head of the hash bucket.
C     After supervariable detection is complete, all hash buckets are
C     empty, and the (LAST(HEAD(HASH)) = 0) condition is restored for
C     the non-empty degree lists.
C
C DENXT(I) need not be set on input. For supervariable I, DENXT(I) is
C     the supervariable following I in a link list, or zero if I is
C     the last in the list. Used for two kinds of lists: degree lists
C     and hash buckets (a supervariable can be in only one kind of
C     list at a time). For element E, DENXT(E) is the number of
C     variables with dense or full rows in the element E.
C
C W(I) need not be set on input. The flag array W determines the status
C     of elements and variables, and the external degree of elements.
C     For elements:
C          if W(E) = 0, then the element E is absorbed.
C          if W(E) .GE. WFLG, then W(E)-WFLG is the size of the set
C               |Le \ Lme|, in terms of nonzeros (the sum of ABS(NV(I))
C               for each principal variable I that is both in the
C               pattern of element E and NOT in the pattern of the
C               current pivot element, ME).
C          if WFLG .GT. WE(E) .GT. 0, then E is not absorbed and has
C               not yet been seen in the scan of the element lists in
C               the computation of |Le\Lme| in loop 150 below.
C               ***SD: change comment to remove reference to label***
C     For variables:
C          during supervariable detection, if W(J) .NE. WFLG then J is
C          not in the pattern of variable I.
C     The W array is initialized by setting W(I) = 1 for all I, and by
C     setting WFLG = 2. It is reinitialized if WFLG becomes too large
C     (to ensure that WFLG+N does not cause integer overflow).
C
C ICNTL is an INTEGER array of length 10 that contains control
C     parameters and must be set by the user. Default values for the
C     components may be set by a call to MC47ID.
C
C RJNFO is an REAL (DOUBLE PRECISION in  D version) array of length 7
C     that need not be set by the user. This array supplies information
C     on the execution of MC47BD.
C     RJNFO(1) gives forecast number of reals to hold the factorization
C     RJNFO(2) gives the forecast number of flops required by the
C     factorization if no pivoting is performed.
C
C Local variables:
C ---------------
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, HASH, HMOD, I,
     $     IDUMMY, ILAST, INEXT, IOVFLO,J, JDUMMY, JLAST, JNEXT, K,
     $     KNT1, KNT2, KNT3, LASTD,  LENJ, LN, MAXMEM, ME,
     $     MEM, MINDEG, NBD, NCMPA, NDME, NEL, NELME, NEWMEM,
     $     NFULL, NLEFT, NRLADU, NVI, NVJ, NVPIV, P, P1, P2, P3, PDST,
     $     PEE, PEE1, PEND, PJ, PME, PME1, PME2, PN, PSRC, RSTRT,
     $     SLENME, THRESH, THRESM, WE, WFLG, WNVI,X
     $
      DOUBLE PRECISION RELDEN, SM, STD, OPS
      LOGICAL IDENSE
C
C DEG:        the degree of a variable or element
C DEGME:      size (no. of variables), |Lme|, of the current element,
C             ME (= DEGREE(ME))
C DEXT:       external degree, |Le \ Lme|, of some element E
C DMAX:       largest |Le| seen so far
C E:          an element
C ELENME:     the length, ELEN(ME), of element list of pivotal var.
C ELN:        the length, ELEN(...), of an element list
C HASH:       the computed value of the hash function
C HMOD:       the hash function is computed modulo HMOD = MAX(1,N-1)
C I:          a supervariable
C IDUMMY:     loop counter
C ILAST:      the entry in a link list preceding I
C INEXT:      the entry in a link list following I
C IOVFLO:     local copy of ICNTL(5)
C J:          a supervariable
C JDUMMY:     loop counter
C JLAST:      the entry in a link list preceding J
C JNEXT:      the entry in a link list, or path, following J
C K:          the pivot order of an element or variable
C KNT1:       loop counter used during element construction
C KNT2:       loop counter used during element construction
C KNT3:       loop counter used during element construction
C LASTD:      index of the last row in the list of dense rows
C LENJ:       LEN(J)
C LN:         length of a supervariable list
C MAXMEM:     amount of memory needed for no compressions
C ME:         current supervariable being eliminated, and the
C                     current element created by eliminating that
C                     supervariable
C MEM:        memory in use assuming no compressions have occurred
C MINDEG:     current approximate minimum degree
C NBD:        total number of dense rows selected
C NCMPA:      counter for the number of times IW was compressed
C NDME  :     number of dense rows adjacent to me
C NEL:        number of pivots selected so far
C NELME:      number of pivots selected when reaching the root
C NEWMEM:     amount of new memory needed for current pivot element
C NFULL:      total number of full rows detected.
C NLEFT:      N-NEL, the number of nonpivotal rows/columns remaining
C NRLADU:     counter for the forecast number of reals in matrix factor
C NVI:        the number of variables in a supervariable I (= NV(I))
C NVJ:        the number of variables in a supervariable J (= NV(J))
C NVPIV:      number of pivots in current element
C P:          pointer into lots of things
C P1:         pe (i) for some variable i (start of element list)
C P2:         pe (i) + elen (i) -  1 for some var. i (end of el. list)
C P3:         index of first supervariable in clean list
C PJ:         pointer into an element or variable
C PDST:       destination pointer, for compression
C PEE:        pointer into element E
C PEE1:       pointer into element E
C PEND:       end of memory to compress
C PME:        pointer into the current element (PME1...PME2)
C PME1:       the current element, ME, is stored in IW(PME1...PME2)
C PME2:       the end of the current element
C PN:         pointer into a "clean" variable, also used to compress
C PSRC:       source pointer, for compression
C RSTRT:      counter for the number of restarts carried out
C SLENME:     number of variables in variable list of pivotal variable
C THRESH:     local copy of ICNTL(4)
C THRESM :    local integer holding the threshold used to detect quasi
C             dense rows. When quasi dense rows are reintegrated in the
C             graph to be processed then THRESM is modified.
C WE:         W(E)
C WFLG:       used for flagging the W array.  See description of W.
C WNVI:       WFLG-NV(I)
C X:          either a supervariable or an element
C
C OPS:        counter for forecast number of flops
C RELDEN :    holds average density to set THRESM automatically
C SM:         counter used for forming standard deviation
C STD:        standard deviation
C
C IDENSE is true if supervariable I is dense
C
C -------------------------------------------------------------------
C  FUNCTIONS CALLED:
C -------------------------------------------------------------------
      INTRINSIC MAX, MIN, MOD
C ====================================================================
C  INITIALIZATIONS
C ====================================================================

      DO 2 I = 1,10
         RJNFO(I) = 0.0
         JNFO(I) = 0
 2    CONTINUE
      DMAX = 0
      HMOD = MAX (1, N-1)
      IOVFLO = ICNTL(5)
      LASTD = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      MINDEG = 1
      NBD   = 0
      NCMPA = 0
      NEL = 0
      NFULL  = 0
      NRLADU = 0
      RSTRT = 0
      OPS = 0.00
      THRESH = ICNTL(4)
      WFLG = 2

C     ------------------------------------------------------
C     Experiments with automatic setting of parameter THRESH.
C     ------------------------------------------------------
      IF (THRESH.GT.0) THEN
         THRESM  = 0
         RELDEN = 0.0
         SM = 0
C        ----------------------------------------------------------
C        initialize arrays and eliminate rows with no off-diag. nz.
C        ----------------------------------------------------------
         DO 5 I=1,N
            THRESM = MAX(THRESM, LEN(I))
            IF (LEN(I).GT.0) THEN
               RELDEN = RELDEN + LEN(I)
               SM = SM + (LEN(I) * LEN(I))
            END IF
            LAST (I) = 0
            HEAD (I) = 0
            NV (I) = 1
            DEGREE (I) = LEN (I)
            IF (DEGREE(I) .EQ. 0) THEN
               NEL = NEL + 1
               ELEN (I) = -NEL
               PE (I) = 0
               W (I) = 0
               NRLADU = NRLADU + 1
               OPS = OPS + 1
            ELSE
               W (I) = 1
               ELEN (I) = 0
            ENDIF
 5       CONTINUE
         IF (N .EQ. NEL) GOTO 265

         RELDEN = RELDEN/(N-NEL)
C        RELDEN holds average row length
         SM = SM/(N-NEL-NFULL) - RELDEN*RELDEN
         STD = SQRT(ABS(SM))
C        STD holds standard deviation of the row lengths
         IF (STD .LE. RELDEN) THEN
            THRESM = -1
         ELSE
            THRESM = INT(9*RELDEN + 0.5*STD*((STD/(RELDEN+0.01))**1.5)+
     *           2*RELDEN*RELDEN/(STD+0.01) +1)
         END IF
C     ------------------------------------------------------
C     end automatic setting of THRESM
C     ------------------------------------------------------

      ELSE
         THRESM = THRESH
         DO 10 I = 1, N
            LAST (I) = 0
            HEAD (I) = 0
            NV (I) = 1
            DEGREE (I) = LEN (I)
            IF (DEGREE(I) .EQ. 0) THEN
               NEL = NEL + 1
               ELEN (I) = -NEL
               PE (I) = 0
               W (I) = 0
               NRLADU = NRLADU + 1
               OPS = OPS + 1
            ELSE
               W (I) = 1
               ELEN (I) = 0
            ENDIF
 10      CONTINUE
      ENDIF
      IF (THRESM.GE.0) THEN
         IF (THRESM.GE.N) THEN
C           full rows only
            THRESM = -1
         ELSE IF (THRESM.EQ.0) THEN
            THRESM = N
         ENDIF
      ENDIF

C     ----------------------------------------------------------------
C     initialize degree lists
C     ----------------------------------------------------------------
      DO 20 I = 1, N
         DEG = DEGREE (I)
         IF (DEG .GT. 0) THEN
C           ----------------------------------------------------------
C           place i in the degree list corresponding to its degree
C           or in the dense row list if i is dense
C           ----------------------------------------------------------
C           test for row density
            IF ( (THRESM.GE.0) .AND.
     &           (DEG+1.GE.THRESM.OR.DEG+1.GE.N-NEL )) THEN
C              I is dense and will be inserted in the degree
C              list of N
               NBD = NBD+1
               IF (DEG+1.NE.N-NEL) THEN
C                 I is quasi dense
                  DEGREE(I) = DEGREE(I)+N+1
C                 insert I at the beginning of degree list of n
                  DEG = N
                  INEXT = HEAD (DEG)
                  IF (INEXT .NE. 0) LAST (INEXT) = I
                  DENXT (I) = INEXT
                  HEAD (DEG) = I
                  LAST(I)  = 0
                  IF (LASTD.EQ.0) THEN
                     LASTD=I
                  END IF
               ELSE
C                 I is full
                  NFULL = NFULL+1
                  DEGREE(I) = N+1
C                 insert I at the end of degree list of n
                  DEG = N
                  IF (LASTD.EQ.0) THEN
C                    degree list is empty
                     LASTD     = I
                     HEAD(DEG) = I
                     DENXT(I)   = 0
                     LAST(I)   = 0
                  ELSE
C                     IF (NFULL.EQ.1) THEN
C                       First full row encountered
                        DENXT(LASTD) = I
                        LAST(I)     = LASTD
                        LASTD       = I
                        DENXT(I)     = 0
C                     ELSE
C                       Absorb I into LASTD (first full row found)
C                        PE(I) = - LASTD
C                        NV(LASTD) = NV(LASTD) + NV(I)
C                        NV(I) = 0
C                        ELEN(I) = 0
C                     END IF
                  ENDIF
               ENDIF
            ELSE
C              place i in the degree list corresponding to its degree
               INEXT = HEAD (DEG)
               IF (INEXT .NE. 0) LAST (INEXT) = I
               DENXT (I) = INEXT
               HEAD (DEG) = I
            ENDIF
         ENDIF
 20   CONTINUE

C     We suppress dense row selection if none of them was found in A
C     in the 1st pass
      IF (NBD.EQ.0 .AND. THRESH.GT.0) THEN
         THRESM = -1
      END IF
C
C ====================================================================
C  WHILE (selecting pivots) DO
C ====================================================================

 30   IF (NEL .LT. N) THEN

C ==================================================================
C  GET PIVOT OF MINIMUM APPROXIMATE DEGREE
C ==================================================================
C       -------------------------------------------------------------
C       find next supervariable for elimination
C       -------------------------------------------------------------
         DO 40 DEG = MINDEG, N
            ME = HEAD (DEG)
            IF (ME .GT. 0) GO TO 50
 40      CONTINUE
 50      MINDEG = DEG
         IF (DEG.LT.N)  THEN
C       -------------------------------------------------------------
C       remove chosen variable from linked list
C       -------------------------------------------------------------
            INEXT = DENXT (ME)
            IF (INEXT .NE. 0) LAST (INEXT) = 0
            HEAD (DEG) = INEXT
         ELSE
            IF (DEGREE(ME).EQ.N+1) GO TO 263
C DEGREE(ME).GT.N+1 so ME is quasi dense
C RESTARTING STRATEGY
C         FOR EACH  quasi dense row d
C           1/ insert d in the degree list according to the
C            value degree(d)-(N+1) (updating MINDEG)
C           2/ Build the adjacency list of d in the quotient graph
C              update DENXT(e_me)= DENXT(e_me)-NV(ME)
C           4/ get back to min degree process
C
C           THRESM > 0 because quasi dense rows were selected
C           While loop: ME is the current dense row
C           make sure that WFLG is not too large
            RSTRT = RSTRT + 1
            RELDEN = 0.0
            SM = 0
            IF (WFLG .GT. IOVFLO-NBD-1) THEN
               DO  51 X = 1, N
                  IF (W (X) .NE. 0) W (X) = 1
 51            CONTINUE
               WFLG = 2
            END IF
            WFLG = WFLG + 1
            DO 57 IDUMMY = 1,N

C           ---------------------------------------------------------
C           remove chosen variable from link list
C           ---------------------------------------------------------
               INEXT = DENXT (ME)
               IF (INEXT .NE. 0) THEN
                  LAST (INEXT) = 0
               ELSE
                  LASTD = 0
               ENDIF
C           ----------------------------------------------------------
c           build adjacency list of ME in quotient graph
C           and calculate its external degree in ndense(me)
C           ----------------------------------------------------------
               DENXT(ME) = 0
C              Flag ME as having been considered in this calculation
               W(ME)      = WFLG
               P1 = PE(ME)
               P2 = P1 + LEN(ME) -1
C           LN-1 holds the pointer in IW to last elt/var in adj list
C              of ME.  LEN(ME) will then be set to LN-P1
C           ELN-1 hold the pointer in IW to  last elt in in adj list
C              of ME.  ELEN(ME) will then be set to ELN-P1
C           element adjacent to ME
               LN       = P1
               ELN      = P1
               DO 55 P=P1,P2
                  E= IW(P)
                  IF (W(E).EQ.WFLG) GO TO 55
                  W(E) = WFLG
C             -------------------------------------------
C             Ensure that E is an unabsorbed element or a quasi dense
C                 row and flag it
C             -------------------------------------------
                  DO 52 JDUMMY = 1,N
                     IF ( PE(E) .GE. 0 ) GOTO 53
                     E = -PE(E)
                     IF (W(E) .EQ.WFLG) GOTO 55
                     W(E) = WFLG
 52               CONTINUE
 53               IF (ELEN(E).LT.0) THEN
C                    E is a new element in adj(ME)
                     DENXT(E) = DENXT(E) - NV(ME)
C                    Move first entry in ME's list of adjacent variables
C                    to the end
                     IW(LN) = IW(ELN)
C                    Place E at end of ME's list of adjacent elements
                     IW(ELN) = E
                     LN  = LN+1
                     ELN = ELN + 1
C                    update ndense of ME with all unflagged dense
C                    rows in E
                     PEE1 = PE(E)
                     DO 54 PEE = PEE1, PEE1+LEN(E)-1
                        X = IW(PEE)
                        IF ((ELEN(X).GE.0).AND.(W(X).NE.WFLG)) THEN
C                          X is a dense row
                           DENXT(ME) = DENXT(ME) + NV(X)
                           W(X) = WFLG
                        ENDIF
 54                  CONTINUE
                  ELSE
C                    E is a dense row
                     DENXT(ME) = DENXT(ME) + NV(E)
C                    Place E at end of ME's list of adjacent variables
                     IW(LN)=E
                     LN = LN+1
                  ENDIF
 55            CONTINUE

C           ----------------------------------------------
C           DEGREE(ME)-(N+1) holds last external degree computed
C           when ME was detected as dense
C           DENXT(ME) is the exact external degree of ME
C           ----------------------------------------------
               WFLG     = WFLG + 1
               LEN(ME)  = LN-P1
               ELEN(ME) = ELN- P1
               NDME = DENXT(ME)+NV(ME)
C            If we want to select ME as full (NDME.EQ.NBD)
C            or quasi dense (NDME.GE.THRESM) then
C            denxt(of elements adjacent to ME) should be updated
               IF (DENXT(ME).EQ.0) DENXT(ME) =1
C           ---------------------------------------------------------
C           place ME in the degree list of DENXT(ME), update DEGREE
C           ---------------------------------------------------------
C               IF (DEGREE(ME)+NV(ME) .LT. NBD  ) THEN
               IF (NDME .LT. NBD) THEN
C                 ME is not full
                  RELDEN = RELDEN + NV(ME)*NDME
                  SM = SM + NV(ME)*NDME*NDME
                  DEGREE(ME) = DENXT(ME)
                  DEG = DEGREE(ME)
                  MINDEG = MIN(DEG,MINDEG)
                  JNEXT = HEAD(DEG)
                  IF (JNEXT.NE. 0) LAST (JNEXT) = ME
                  DENXT(ME) = JNEXT
                  HEAD(DEG) = ME
               ELSE
C                 ME is full
                  DEGREE(ME) = N+1
                  DEG = DENXT(ME)
                  MINDEG = MIN(DEG,MINDEG)
                  DEG = N

C                 Update DENXT of all elements in the list of elements
C                 adjacent to ME
                  P1 = PE(ME)
                  P2 = P1 + ELEN(ME) - 1
                  DO 56 PJ=P1,P2
                     E= IW(PJ)
                     DENXT (E) = DENXT(E) + NV(ME)
 56               CONTINUE
C                  insert ME in the list of dense rows
                  DEG = N
C                 ME at the end of the list
                  NFULL = NFULL +NV(ME)
                  IF (LASTD.EQ.0) THEN
C                        degree list is empty
                     LASTD     = ME
                     HEAD(N) = ME
                     DENXT(ME)   = 0
                     LAST(ME)   = 0
                     IF (INEXT.EQ.0) INEXT = LASTD
                  ELSE
C                     IF (NFULL.EQ.NV(ME)) THEN
C                           First full row encountered
                        DENXT(LASTD) = ME
                        LAST(ME)     = LASTD
                        LASTD        = ME
                        DENXT(ME)     = 0
                        IF (INEXT.EQ.0) INEXT = LASTD
C                     ELSE
C                   Absorb ME into LASTD (first full row found)
C                        PE(ME) = - LASTD
C                        NV(LASTD) = NV(LASTD) + NV(ME)
C                        NV(ME) = 0
C                        ELEN(ME) = 0
C                     END IF
                  ENDIF
               END IF

C           ------------------------------
C           process next quasi dense row
C           ------------------------------
               ME    = INEXT
               IF (ME.EQ.0) GO TO 58
               IF (DEGREE(ME).LE.(N+1) ) GOTO 58
 57         CONTINUE
 58         HEAD (N) = ME
C           ---------------------------------------
C           update dense row selection strategy
C           -------------------------------------
            IF (NBD.EQ.NFULL) THEN
               RELDEN = 0
               SM = 0
            ELSE
               RELDEN = (RELDEN + NFULL*NBD)/(NBD)
               SM = (SM + NFULL*NBD*NBD)/(NBD) - RELDEN*RELDEN
            END IF
            STD = SQRT(ABS(SM))
            THRESM = INT(9*RELDEN+0.5*STD*((STD/(RELDEN + 0.01))**1.5)
     *           + 2*RELDEN*RELDEN/(STD+0.01) +1)
            THRESM = MIN(THRESM,NBD)
            IF (THRESM.GE.NBD) THEN
               THRESM = N
            END IF
            NBD = NFULL
C           get back to min degree elimination loop
            GOTO 30
C         -------------------------------------------------------------
C         -------------------------------------------------------------
         ENDIF
C       -------------------------------------------------------------
C       me represents the elimination of pivots nel+1 to nel+nv(me).
C       place me itself as the first in this set.  It will be moved
C       to the nel+nv(me) position when the permutation vectors are
C       computed.
C       -------------------------------------------------------------
         ELENME = ELEN (ME)
         ELEN (ME) = - (NEL + 1)
         NVPIV = NV (ME)
         NEL = NEL + NVPIV
         DENXT(ME) = 0

C ====================================================================
C  CONSTRUCT NEW ELEMENT
C ====================================================================
C
C       -------------------------------------------------------------
C       At this point, me is the pivotal supervariable.  It will be
C       converted into the current element.  Scan list of the
C       pivotal supervariable, me, setting tree pointers and
C       constructing new list of supervariables for the new element,
C       me.  p is a pointer to the current position in the old list.
C       -------------------------------------------------------------
C
C       flag the variable "me" as being in the front by negating nv(me)
         NV (ME) = -NVPIV
         DEGME = 0
         IF (ELENME .EQ. 0) THEN
C         ----------------------------------------------------------
C         There are no elements involved.
C         Construct the new element in place.
C         ----------------------------------------------------------
            PME1 = PE (ME)
            PME2 = PME1 - 1
            DO 60 P = PME1, PME1 + LEN (ME) - 1
               I = IW (P)
               NVI = NV (I)
               IF (NVI .GT. 0) THEN
C             ----------------------------------------------------
C             i is a principal variable not yet placed in the
C             generated element. Store i in new list
C             ----------------------------------------------------
                  DEGME = DEGME + NVI
C                 flag i as being in Lme by negating nv (i)
                  NV (I) = -NVI
                  PME2 = PME2 + 1
                  IW (PME2) = I

C             ----------------------------------------------------
C             remove variable i from degree list.
C             ----------------------------------------------------
C             only done for sparse rows
                  IF (DEGREE(I).LE.N) THEN
                     ILAST = LAST (I)
                     INEXT = DENXT (I)
                     IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                     IF (ILAST .NE. 0) THEN
                        DENXT (ILAST) = INEXT
                     ELSE
C                       i is at the head of the degree list
                        HEAD (DEGREE (I)) = INEXT
                     ENDIF
                  ELSE
C                    Dense rows remain dense so do not remove from list
                     DENXT(ME) = DENXT(ME) + NVI
                  ENDIF
               ENDIF
 60         CONTINUE
C           this element takes no new memory in iw:
            NEWMEM = 0
         ELSE
C         ----------------------------------------------------------
C         construct the new element in empty space, iw (pfree ...)
C         ----------------------------------------------------------
            P  = PE (ME)
            PME1 = PFREE
            SLENME = LEN (ME) - ELENME
            DO 120 KNT1 = 1, ELENME
C              search the elements in me.
               E = IW (P)
               P = P + 1
               PJ = PE (E)
               LN = LEN (E)
C           -------------------------------------------------------
C           search for different supervariables and add them to the
C           new list, compressing when necessary.
C           -------------------------------------------------------
               DO 110 KNT2 = 1, LN
                  I = IW (PJ)
                  PJ = PJ + 1
                  NVI = NV (I)
                  IF (NVI .GT. 0) THEN
C               -------------------------------------------------
C               compress iw, if necessary
C               -------------------------------------------------
                     IF (PFREE .GT. IWLEN) THEN
C                    prepare for compressing iw by adjusting
C                    pointers and lengths so that the lists being
C                    searched in the inner and outer loops contain
C                    only the remaining entries.
C                    ***** SD: Seperate compression subroutine tried
C                      but found to be inefficient in comparison ****
                        PE (ME) = P
                        LEN (ME) = LEN (ME) - KNT1
C                       Check if anything left in supervariable ME
                        IF (LEN (ME) .EQ. 0) PE (ME) = 0
                        PE (E) = PJ
                        LEN (E) = LN - KNT2
C                       Check if anything left in element E
                        IF (LEN (E) .EQ. 0) PE (E) = 0
                        NCMPA = NCMPA + 1
C                       store first item in pe
C                       set first entry to -item
                        DO 70 J = 1, N
                           PN = PE (J)
                           IF (PN .GT. 0) THEN
                              PE (J) = IW (PN)
                              IW (PN) = -J
                           ENDIF
 70                     CONTINUE

C                       psrc/pdst point to source/destination
                        PDST = 1
                        PSRC = 1
                        PEND = PME1 - 1

C                       while loop:
                        DO 91 IDUMMY = 1, IWLEN
                           IF (PSRC .GT. PEND) GO TO 95
C                          search for next negative entry
                           J = -IW (PSRC)
                           PSRC = PSRC + 1
                           IF (J .GT. 0) THEN
                              IW (PDST) = PE (J)
                              PE (J) = PDST
                              PDST = PDST + 1
C                     copy from source to destination
                              LENJ = LEN (J)
                              DO 90 KNT3 = 0, LENJ - 2
                                 IW (PDST + KNT3) = IW (PSRC + KNT3)
 90                           CONTINUE
                              PDST = PDST + LENJ - 1
                              PSRC = PSRC + LENJ - 1
                           ENDIF
 91                     END DO

C                       move the new partially-constructed element
 95                     P1 = PDST
                        DO 100 PSRC = PME1, PFREE - 1
                           IW (PDST) = IW (PSRC)
                           PDST = PDST + 1
 100                    CONTINUE
                        PME1 = P1
                        PFREE = PDST
                        PJ = PE (E)
                        P = PE (ME)
                     ENDIF

C               -------------------------------------------------
C               i is a principal variable not yet placed in Lme
C               store i in new list
C               -------------------------------------------------
                     DEGME = DEGME + NVI
C                    flag i as being in Lme by negating nv (i)
                     NV (I) = -NVI
                     IW (PFREE) = I
                     PFREE = PFREE + 1

C               -------------------------------------------------
C               remove variable i from degree link list
C               -------------------------------------------------
C                    only done for sparse rows
                     IF (DEGREE(I).LE.N) THEN
                        ILAST = LAST (I)
                        INEXT = DENXT (I)
                        IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                        IF (ILAST .NE. 0) THEN
                           DENXT (ILAST) = INEXT
                        ELSE
C                         i is at the head of the degree list
                           HEAD (DEGREE (I)) = INEXT
                        ENDIF
                     ELSE
C                    Dense rows remain dense so do not remove from list
                        DENXT(ME) = DENXT(ME) + NVI
                     ENDIF
                  ENDIF
 110           CONTINUE

C                set tree pointer and flag to indicate element e is
C                absorbed into new element me (the parent of e is me)
                  PE (E) = -ME
                  W (E) = 0
 120        CONTINUE

C           search the supervariables in me.
            KNT1 = ELENME + 1
            E = ME
            PJ = P
            LN = SLENME

C           -------------------------------------------------------
C           search for different supervariables and add them to the
C           new list, compressing when necessary.
C           -------------------------------------------------------
            DO 126 KNT2 = 1, LN
               I = IW (PJ)
               PJ = PJ + 1
               NVI = NV (I)
               IF (NVI .GT. 0) THEN
C               -------------------------------------------------
C               compress iw, if necessary
C               -------------------------------------------------
                  IF (PFREE .GT. IWLEN) THEN
C                 prepare for compressing iw by adjusting
C                 pointers and lengths so that the lists being
C                 searched in the inner and outer loops contain
C                 only the remaining entries.
                     PE (ME) = P
                     LEN (ME) = LEN (ME) - KNT1
C                    Check if anything left in supervariable ME
                     IF (LEN (ME) .EQ. 0) PE (ME) = 0
                     PE (E) = PJ
                     LEN (E) = LN - KNT2
C                    Check if anything left in element E
                     IF (LEN (E) .EQ. 0) PE (E) = 0
                     NCMPA = NCMPA + 1
C                    store first item in pe
C                    set first entry to -item
                     DO 121 J = 1, N
                        PN = PE (J)
                        IF (PN .GT. 0) THEN
                           PE (J) = IW (PN)
                           IW (PN) = -J
                        ENDIF
 121                 CONTINUE

C                    psrc/pdst point to source/destination
                     PDST = 1
                     PSRC = 1
                     PEND = PME1 - 1

C                 while loop:
C 122              CONTINUE
                     DO 123 IDUMMY = 1,IWLEN
                        IF (PSRC .GT. PEND) GO TO 124
C                       search for next negative entry
                        J = -IW (PSRC)
                        PSRC = PSRC + 1
                        IF (J .GT. 0) THEN
                           IW (PDST) = PE (J)
                           PE (J) = PDST
                           PDST = PDST + 1
C                          copy from source to destination
                           LENJ = LEN (J)
                           DO 122 KNT3 = 0, LENJ - 2
                              IW (PDST + KNT3) = IW (PSRC + KNT3)
 122                       CONTINUE
                           PDST = PDST + LENJ - 1
                           PSRC = PSRC + LENJ - 1
                        ENDIF
 123                 END DO

C                 move the new partially-constructed element
 124                 P1 = PDST
                     DO 125 PSRC = PME1, PFREE - 1
                        IW (PDST) = IW (PSRC)
                        PDST = PDST + 1
 125                 CONTINUE
                     PME1 = P1
                     PFREE = PDST
                     PJ = PE (E)
                     P = PE (ME)
                  END IF

C               -------------------------------------------------
C               i is a principal variable not yet placed in Lme
C               store i in new list
C               -------------------------------------------------
                  DEGME = DEGME + NVI
C                 flag i as being in Lme by negating nv (i)
                  NV (I) = -NVI
                  IW (PFREE) = I
                  PFREE = PFREE + 1

C               -------------------------------------------------
C               remove variable i from degree link list
C               -------------------------------------------------
C                 only done for sparse rows
                  IF (DEGREE(I).LE.N) THEN
                     ILAST = LAST (I)
                     INEXT = DENXT (I)
                     IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                     IF (ILAST .NE. 0) THEN
                        DENXT (ILAST) = INEXT
                     ELSE
C                 i is at the head of the degree list
                        HEAD (DEGREE (I)) = INEXT
                     ENDIF
                  ELSE
C                    Dense rows remain dense so do not remove from list
                     DENXT(ME) = DENXT(ME) + NVI
                  ENDIF
               ENDIF
 126           CONTINUE

            PME2 = PFREE - 1
C           this element takes newmem new memory in iw (possibly zero)
            NEWMEM = PFREE - PME1
            MEM = MEM + NEWMEM
            MAXMEM = MAX (MAXMEM, MEM)
         ENDIF

C       -------------------------------------------------------------
C       me has now been converted into an element in iw (pme1..pme2)
C       -------------------------------------------------------------
C       degme holds the external degree of new element
         DEGREE (ME) = DEGME
         PE (ME) = PME1
         LEN (ME) = PME2 - PME1 + 1

C       -------------------------------------------------------------
C       make sure that wflg is not too large.  With the current
C       value of wflg, wflg+n must not cause integer overflow
C       -------------------------------------------------------------
         IF (WFLG .GT. IOVFLO-N) THEN
            DO 130 X = 1, N
               IF (W (X) .NE. 0) W (X) = 1
 130        CONTINUE
            WFLG = 2
         ENDIF

C ====================================================================
C   COMPUTE (w(e) - wflg) = |Le(G')\Lme(G')| FOR ALL ELEMENTS
C   where G' is the subgraph of G containing just the sparse rows)
C ====================================================================
C       -------------------------------------------------------------
C       Scan 1:  compute the external degrees of elements touched
C       with respect to the current element.  That is:
C            (w (e) - wflg) = |Le \ Lme|
C       for each element e involving a supervariable in Lme.
C       The notation Le refers to the pattern (list of
C       supervariables) of a previous element e, where e is not yet
C       absorbed, stored in iw (pe (e) + 1 ... pe (e) + iw (pe (e))).
C       The notation Lme refers to the pattern of the current element
C       (stored in iw (pme1..pme2)).
C       aggressive absorption is possible only if DENXT(ME) = NBD
C       which is true when only full rows have been selected.
C       -------------------------------------------------------------
         IF (NBD.GT.0) THEN
C           Dense rows have been found
            DO 150 PME = PME1, PME2
               I = IW (PME)
c              skip dense rows
               IF (DEGREE(I).GT.N) GOTO 150
               ELN = ELEN (I)
               IF (ELN .GT. 0) THEN
C              note that nv (i) has been negated to denote i in Lme:
                  NVI = -NV (I)
                  WNVI = WFLG - NVI
                  DO 140 P = PE (I), PE (I) + ELN - 1
                     E = IW (P)
                     WE = W (E)
                     IF (WE .GE. WFLG) THEN
C                    unabsorbed element e has been seen in this loop
                        WE = WE - NVI
                     ELSE IF (WE .NE. 0) THEN
C                    e is an unabsorbed element - this is
C                    the first we have seen e in all of Scan 1
                        WE = DEGREE (E) + WNVI - DENXT(E)
                     ENDIF
                     W (E) = WE
 140              CONTINUE
               ENDIF
 150        CONTINUE
         ELSE
C           No dense rows have been found
            DO 152 PME = PME1, PME2
               I = IW (PME)
               ELN = ELEN (I)
               IF (ELN .GT. 0) THEN
C           note that nv (i) has been negated to denote i in Lme:
                  NVI = -NV (I)
                  WNVI = WFLG - NVI
                  DO 151 P = PE (I), PE (I) + ELN - 1
                     E = IW (P)
                     WE = W (E)
                     IF (WE .GE. WFLG) THEN
C                    unabsorbed element e has been seen in this loop
                        WE = WE - NVI
                     ELSE IF (WE .NE. 0) THEN
C                    e is an unabsorbed element - this is
C                    the first we have seen e in all of Scan 1
                        WE = DEGREE (E) + WNVI
                     ENDIF
                     W (E) = WE
 151              CONTINUE
               ENDIF
 152        CONTINUE
         END IF

C ====================================================================
C  DEGREE UPDATE AND ELEMENT ABSORPTION
C ====================================================================

C       -------------------------------------------------------------
C       Scan 2:  for each sparse i in Lme, sum up the external degrees
C       of each Le for the elements e appearing within i, plus the
C       supervariables in i.  Place i in hash list.
C       -------------------------------------------------------------
         IF (NBD.GT.0) THEN
C           Dense rows have been found
            DO 180 PME = PME1, PME2
               I = IW (PME)
C              skip dense rows
               IF (DEGREE(I).GT.N) GOTO 180
C              remove absorbed elements from the list for i
               P1 = PE (I)
               P2 = P1 + ELEN (I) - 1
               PN = P1
               HASH = 0
               DEG = 0

C         ----------------------------------------------------------
C         scan the element list associated with supervariable i
C         ----------------------------------------------------------
               DO 160 P = P1, P2
                  E = IW (P)
C                 dext = | Le | - | (Le \cap Lme)\D | - DENXT(e)
                  DEXT = W (E) - WFLG
                  IF (DEXT .GT. 0) THEN
                     DEG = DEG + DEXT
                     IW (PN) = E
                     PN = PN + 1
                     HASH = HASH+E
                  ELSE IF ((DEXT .EQ. 0) .AND.
     &                    (DENXT(ME).EQ.NBD)) THEN
C             aggressive absorption: e is not adjacent to me, but
C             |Le(G') \ Lme(G')| is 0 and all dense rows
C             are in me, so absorb it into me
                     PE (E) = -ME
                     W (E)  = 0
                  ELSE IF (DEXT.EQ.0) THEN
                     IW(PN) = E
                     PN     = PN+1
                     HASH = HASH + E
                  ENDIF
 160           CONTINUE

C           count the number of elements in i (including me):
               ELEN (I) = PN - P1 + 1

C         ----------------------------------------------------------
C         scan the supervariables in the list associated with i
C         ----------------------------------------------------------
               P3 = PN
               DO 170 P = P2 + 1, P1 + LEN (I) - 1
                  J = IW (P)
                  NVJ = NV (J)
                  IF (NVJ .GT. 0) THEN
C                 j is unabsorbed, and not in Lme.
C                 add to degree and add to new list
C                 add degree only of sparse rows.
                     IF (DEGREE(J).LE.N) DEG=DEG+NVJ
                     IW (PN) = J
                     PN = PN + 1
                     HASH = HASH + J
                  ENDIF
 170           CONTINUE

C         ----------------------------------------------------------
C         update the degree and check for mass elimination
C         ----------------------------------------------------------
               IF ((DEG .EQ. 0).AND.(DENXT(ME).EQ.NBD)) THEN
C              mass elimination only possible when all dense rows
C              are in ME
C           -------------------------------------------------------
C           mass elimination - supervariable i can be eliminated
C           -------------------------------------------------------
                  PE (I) = -ME
                  NVI = -NV (I)
                  DEGME = DEGME - NVI
                  NVPIV = NVPIV + NVI
                  NEL = NEL + NVI
                  NV (I) = 0
                  ELEN (I) = 0
               ELSE
C           -------------------------------------------------------
C           update the upper-bound degree of i
C           A bound for the new external degree is the old bound plus
C           the size of the generated element
C           -------------------------------------------------------
C           the following degree does not yet include the size
C           of the current element, which is added later:
                  DEGREE(I) = MIN (DEG+NBD-DENXT(ME), DEGREE(I))

C           -------------------------------------------------------
C           add me to the list for i
C           -------------------------------------------------------
C              move first supervariable to end of list
                  IW (PN) = IW (P3)
C              move first element to end of element part of list
                  IW (P3) = IW (P1)
C              add new element to front of list.
                  IW (P1) = ME
C              store the new length of the list in len (i)
                  LEN (I) = PN - P1 + 1

C           -------------------------------------------------------
C           place in hash bucket.  Save hash key of i in last (i).
C           -------------------------------------------------------
                  HASH = ABS(MOD (HASH, HMOD)) + 1
                  J = HEAD (HASH)
                  IF (J .LE. 0) THEN
C                the degree list is empty, hash head is -j
                     DENXT (I) = -J
                     HEAD (HASH) = -I
                  ELSE
C                degree list is not empty - has j as its head
C                last is hash head
                     DENXT (I) = LAST (J)
                     LAST (J) = I
                  ENDIF
                  LAST (I) = HASH
               ENDIF
 180        CONTINUE
         ELSE
C           No dense rows have been found
            DO 183 PME = PME1, PME2
               I = IW (PME)
C              remove absorbed elements from the list for i
               P1 = PE (I)
               P2 = P1 + ELEN (I) - 1
               PN = P1
               HASH = 0
               DEG = 0

C              -------------------------------------------------------
C              scan the element list associated with supervariable i
C              -------------------------------------------------------
               DO 181 P = P1, P2
                  E = IW (P)
C                 dext = | Le | - | (Le \cap Lme)\D | - DENXT(e)
                  DEXT = W (E) - WFLG
                  IF (DEXT .GT. 0) THEN
                     DEG = DEG + DEXT
                     IW (PN) = E
                     PN = PN + 1
                     HASH = HASH + E
                  ELSE IF (DEXT .EQ. 0) THEN
C                  aggressive absorption: e is not adjacent to me, but
C                  |Le(G') \ Lme(G')| is 0, so absorb it into me
                     PE (E) = -ME
                     W (E)  = 0
                  ENDIF
 181           CONTINUE

C           count the number of elements in i (including me):
               ELEN (I) = PN - P1 + 1

C         ----------------------------------------------------------
C         scan the supervariables in the list associated with i
C         ----------------------------------------------------------
               P3 = PN
               DO 182 P = P2 + 1, P1 + LEN (I) - 1
                  J = IW (P)
                  NVJ = NV (J)
                  IF (NVJ .GT. 0) THEN
C                 j is unabsorbed, and not in Lme.
C                 add to degree and add to new list
                     DEG=DEG+NVJ
                     IW (PN) = J
                     PN = PN + 1
                     HASH = HASH + J
                  ENDIF
 182           CONTINUE

C         ----------------------------------------------------------
C         update the degree and check for mass elimination
C         ----------------------------------------------------------
               IF (DEG .EQ. 0) THEN
C           -------------------------------------------------------
C           mass elimination - supervariable i can be eliminated
C           -------------------------------------------------------
                  PE (I) = -ME
                  NVI = -NV (I)
                  DEGME = DEGME - NVI
                  NVPIV = NVPIV + NVI
                  NEL = NEL + NVI
                  NV (I) = 0
                  ELEN (I) = 0
               ELSE
C           -------------------------------------------------------
C           update the upper-bound degree of i
C           A bound for the new external degree is the old bound plus
C           the size of the generated element
C           -------------------------------------------------------
C
C           the following degree does not yet include the size
C           of the current element, which is added later:
                  DEGREE(I) = MIN (DEG,  DEGREE(I))

C           -------------------------------------------------------
C           add me to the list for i
C           -------------------------------------------------------
C              move first supervariable to end of list
                  IW (PN) = IW (P3)
C              move first element to end of element part of list
                  IW (P3) = IW (P1)
C              add new element to front of list.
                  IW (P1) = ME
C              store the new length of the list in len (i)
                  LEN (I) = PN - P1 + 1

C           -------------------------------------------------------
C           place in hash bucket.  Save hash key of i in last (i).
C           -------------------------------------------------------
                  HASH = ABS(MOD (HASH, HMOD)) + 1
                  J = HEAD (HASH)
                  IF (J .LE. 0) THEN
C                the degree list is empty, hash head is -j
                     DENXT (I) = -J
                     HEAD (HASH) = -I
                  ELSE
C                degree list is not empty - has j as its head
C                last is hash head
                     DENXT (I) = LAST (J)
                     LAST (J) = I
                  ENDIF
                  LAST (I) = HASH
               ENDIF
 183        CONTINUE
         END IF
         DEGREE (ME) = DEGME

C       -------------------------------------------------------------
C       Clear the counter array, w (...), by incrementing wflg.
C       -------------------------------------------------------------
         DMAX = MAX (DMAX, DEGME)
         WFLG = WFLG + DMAX

C        make sure that wflg+n does not cause integer overflow
         IF (WFLG .GE. IOVFLO - N) THEN
            DO 190 X = 1, N
               IF (W (X) .NE. 0) W (X) = 1
 190        CONTINUE
            WFLG = 2
         ENDIF
C        at this point, w (1..n) .lt. wflg holds

C ====================================================================
C  SUPERVARIABLE DETECTION
C ====================================================================
         DO 250 PME = PME1, PME2
            I = IW (PME)
            IF ( (NV(I).GE.0) .OR. (DEGREE(I).GT.N) ) GO TO 250
C           only done for sparse rows
C           replace i by head of its hash bucket, and set the hash
C           bucket header to zero

C           -------------------------------------------------------
C           examine all hash buckets with 2 or more variables.  We
C           do this by examing all unique hash keys for super-
C           variables in the pattern Lme of the current element, me
C           -------------------------------------------------------
            HASH = LAST (I)
C           let i = head of hash bucket, and empty the hash bucket
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
C             degree list is empty
               I = -J
               HEAD (HASH) = 0
            ELSE
C             degree list is not empty, restore last () of head
               I = LAST (J)
               LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250

C           while loop:
            DO 247 JDUMMY = 1,N
               IF (DENXT (I) .EQ. 0) GO TO 250
C             ----------------------------------------------------
C             this bucket has one or more variables following i.
C             scan all of them to see if i can absorb any entries
C             that follow i in hash bucket.  Scatter i into w.
C             ----------------------------------------------------
               LN = LEN (I)
               ELN = ELEN (I)
C              do not flag the first element in the list (me)
               DO 210 P = PE (I) + 1, PE (I) + LN - 1
                  W (IW (P)) = WFLG
 210           CONTINUE

C             ----------------------------------------------------
C             scan every other entry j following i in bucket
C             ----------------------------------------------------
               JLAST = I
               J = DENXT (I)

C             while loop:
               DO 245 IDUMMY=1,N
                  IF (J .EQ. 0) GO TO 246

C               -------------------------------------------------
C               check if j and i have identical nonzero pattern
C               -------------------------------------------------
C               jump if i and j do not have same size data structure
                  IF (LEN (J) .NE. LN) GO TO 240
C               jump if i and j do not have same number adj elts
                  IF (ELEN (J) .NE. ELN) GO TO 240
C               do not flag the first element in the list (me)

                  DO 230 P = PE (J) + 1, PE (J) + LN - 1
C                 jump if an entry (iw(p)) is in j but not in i
                     IF (W (IW (P)) .NE. WFLG) GO TO 240
 230              CONTINUE

C                 -------------------------------------------------
C                 found it!  j can be absorbed into i
C                 -------------------------------------------------
                  PE (J) = -I
C                 both nv (i) and nv (j) are negated since they
C                 are in Lme, and the absolute values of each
C                 are the number of variables in i and j:
                  NV (I) = NV (I) + NV (J)
                  NV (J) = 0
                  ELEN (J) = 0
C                 delete j from hash bucket
                  J = DENXT (J)
                  DENXT (JLAST) = J
                  GO TO 245

C                 -------------------------------------------------
 240              CONTINUE
C               j cannot be absorbed into i
C               -------------------------------------------------
                  JLAST = J
                  J = DENXT (J)
 245           CONTINUE

C             ----------------------------------------------------
C             no more variables can be absorbed into i
C             go to next i in bucket and clear flag array
C             ----------------------------------------------------
 246           WFLG = WFLG + 1
               I = DENXT (I)
               IF (I .EQ. 0) GO TO 250
 247         CONTINUE
 250      CONTINUE

C ====================================================================
C  RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMENT
C  Squeeze out absorbed variables
C ====================================================================
          P = PME1
          NLEFT = N - NEL
          DO 260 PME = PME1, PME2
             I = IW (PME)
             NVI = -NV (I)
             IF (NVI .LE. 0) GO TO 260
C            i is a principal variable in Lme
C            restore nv (i) to signify that i is principal
             NV (I) = NVI
             IF (DEGREE(I).GT.N) GO TO 258
C           -------------------------------------------------------
C           compute the external degree (add size of current elem)
C           -------------------------------------------------------
             DEG = MIN (DEGREE (I)+ DEGME - NVI, NLEFT - NVI)
             DEGREE (I) = DEG
             IDENSE = .FALSE.
C           -------------------
C           Dense row detection
C           -------------------
             IF (THRESM.GE.0) THEN
C           DEGME is exact external degree of pivot ME |Le\Ve|,
C           DEG is is approx external degree of I
                IF ((DEG+NVI .GE. THRESM).OR.
     &               (DEG+NVI .GE. NLEFT)) THEN
                   IF (THRESM.EQ.N) THEN
C                    We must be sure that I is full in reduced matrix
                      IF ((ELEN(I).LE.2) .AND.((DEG+NVI).EQ.NLEFT)
     &                     .AND. NBD.EQ.NFULL ) THEN
C                        DEG approximation is exact and I is dense
                         DEGREE(I) = N+1
                         IDENSE = .TRUE.
                      ENDIF
                   ELSE
C                     relaxed dense row detection
                      IDENSE = .TRUE.
                      IF ((ELEN(I).LE.2).AND. ((DEG+NVI).EQ.NLEFT)
     &                     .AND. NBD.EQ.NFULL ) THEN
                         DEGREE(I) = N+1
                      ELSE
                         DEGREE(I) = N+1+DEGREE(I)
                      ENDIF
                   ENDIF
                ENDIF
                IF (IDENSE) THEN
C                  update DENXT of all elements in the list of element
C                  adjacent to I (including ME).
                   P1 = PE(I)
                   P2 = P1 + ELEN(I) - 1
                   DO 255 PJ=P1,P2
                      E= IW(PJ)
                      DENXT (E) = DENXT(E) + NVI
 255               CONTINUE
C                  insert I in the list of dense rows
                   NBD = NBD+NVI
                   DEG = N
                   IF (DEGREE(I).EQ.N+1) THEN
c                     insert I at the end of the list
                      NFULL = NFULL +NVI
                      IF (LASTD.EQ.0) THEN
C                        degree list is empty
                         LASTD     = I
                         HEAD(DEG) = I
                         DENXT(I)   = 0
                         LAST(I)   = 0
                      ELSE
C                         IF (NFULL.EQ.NVI) THEN
C                           First full row encountered
                            DENXT(LASTD) = I
                            LAST(I)     = LASTD
                            LASTD       = I
                            DENXT(I)     = 0
C                         ELSE
C                           Absorb I into LASTD (first full row found)
C                            PE(I) = - LASTD
C                            NV(LASTD) = NV(LASTD) + NV(I)
C                            NV(I) = 0
C                            ELEN(I) = 0
C                         END IF
                      ENDIF
                   ELSE
C                     insert I at the beginning of the list
                      INEXT = HEAD(DEG)
                      IF (INEXT .NE. 0) LAST (INEXT) = I
                      DENXT (I) = INEXT
                      HEAD (DEG) = I
                      LAST(I)    = 0
                      IF (LASTD.EQ.0) LASTD=I
                   ENDIF
C               end of IDENSE=true
                ENDIF
C            end of THRESM>0
             ENDIF

             IF (.NOT.IDENSE) THEN
C           -------------------------------------------------------
C           place the supervariable at the head of the degree list
C           -------------------------------------------------------
                INEXT = HEAD (DEG)
                IF (INEXT .NE. 0) LAST (INEXT) = I
                DENXT (I) = INEXT
                LAST (I) = 0
                HEAD (DEG) = I
             ENDIF
C           -------------------------------------------------------
C           save the new degree, and find the minimum degree
C           -------------------------------------------------------
             MINDEG = MIN (MINDEG, DEG)
 258         CONTINUE
C           -------------------------------------------------------
C           place the supervariable in the element pattern
C           -------------------------------------------------------
             IW (P) = I
             P = P + 1
 260      CONTINUE

C =====================================================================
C  FINALIZE THE NEW ELEMENT
C =====================================================================
          OPS = OPS + DEGME*NVPIV + DEGME * NVPIV*NVPIV +
     *         DEGME*DEGME*NVPIV + NVPIV*NVPIV*NVPIV/3 +
     *         NVPIV*NVPIV/2 + NVPIV/6 + NVPIV
          NRLADU = NRLADU + (NVPIV*(NVPIV+1))/2 + (DEGME*NVPIV)
          NV (ME) = NVPIV + DEGME
C         nv (me) is now the degree of pivot (including diagonal part)
C         save the length of the list for the new element me
          LEN (ME) = P - PME1
          IF (LEN (ME) .EQ. 0) THEN
C            there is nothing left of the current pivot element
             PE (ME) = 0
             W (ME) = 0
          ENDIF
          IF (NEWMEM .NE. 0) THEN
C            element was not constructed in place: deallocate part
C            of it (final size is less than or equal to newmem,
C            since newly nonprincipal variables have been removed).
             PFREE = P
             MEM = MEM - NEWMEM + LEN (ME)
          ENDIF

C =====================================================================
C       END WHILE (selecting pivots)
          GO TO 30
       ENDIF
C =====================================================================
       GO TO 265

C      We have only full rows that we amalgamate at the root
C      node and ME = LASTD
C      Perform mass elimination of all full rows
 263   NELME    = -(NEL+1)
       DO 264 X=1,N
          IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
C            X is an unabsorbed element
             PE(X) = -ME
          ELSEIF (DEGREE(X).EQ.N+1) THEN
C            X is a dense row, absorb it in ME (mass elimination)
             NEL   = NEL + NV(X)
             PE(X) = -ME
             ELEN(X) = 0
             NV(X) = 0
          ENDIF
 264   CONTINUE
C      ME is the root node
       ELEN(ME) = NELME
       NV(ME)   = NBD
       NRLADU = NRLADU + (NBD*(NBD+1))/2
       OPS = OPS + NBD*NBD*NBD/3 + NBD*NBD/2 + NBD/6 + NBD
       PE(ME)   = 0

 265   CONTINUE
C ===================================================================
C  COMPUTE THE PERMUTATION VECTORS
C ===================================================================

C     ----------------------------------------------------------------
C     The time taken by the following code is O(n).  At this
C     point, elen (e) = -k has been done for all elements e,
C     and elen (i) = 0 has been done for all nonprincipal
C     variables i.  At this point, there are no principal
C     supervariables left, and all elements are absorbed.
C     ----------------------------------------------------------------
C
C     ----------------------------------------------------------------
C     compute the ordering of unordered nonprincipal variables
C     ----------------------------------------------------------------

       DO 290 I = 1, N
          IF (ELEN (I) .EQ. 0) THEN
C         ----------------------------------------------------------
C         i is an un-ordered row.  Traverse the tree from i until
C         reaching an element, e.  The element, e, was the
C         principal supervariable of i and all nodes in the path
C         from i to when e was selected as pivot.
C         ----------------------------------------------------------
             J = -PE (I)
C         while (j is a variable) do:
             DO 270 JDUMMY = 1,N
                IF (ELEN (J) .LT. 0) GO TO 275
                J = -PE (J)
 270         CONTINUE
 275         E = J
C           ----------------------------------------------------------
C           get the current pivot ordering of e
C           ----------------------------------------------------------
             K = -ELEN (E)

C           ----------------------------------------------------------
C           traverse the path again from i to e, and compress the
C           path (all nodes point to e).  Path compression allows
C           this code to compute in O(n) time.  Order the unordered
C           nodes in the path, and place the element e at the end.
C           ----------------------------------------------------------
             J = I
C            while (j is a variable) do:
             DO 280 IDUMMY = 1,N
                IF (ELEN (J) .LT. 0) GO TO 285
                JNEXT = -PE (J)
                PE (J) = -E
                IF (ELEN (J) .EQ. 0) THEN
C                  j is an unordered row
                   ELEN (J) = K
                   K = K + 1
                ENDIF
                J = JNEXT
 280         CONTINUE
C            leave elen (e) negative, so we know it is an element
 285         ELEN (E) = -K
          ENDIF
 290   CONTINUE

C     ----------------------------------------------------------------
C     reset the inverse permutation (elen (1..n)) to be positive,
C     and compute the permutation (last (1..n)).
C     ----------------------------------------------------------------
       DO 300 I = 1, N
          K = ABS (ELEN (I))
          LAST (K) = I
          ELEN (I) = K
 300   CONTINUE

C ====================================================================
C  RETURN THE MEMORY USAGE IN IW AND SET INFORMATION ARRAYS
C ====================================================================
C     If maxmem is less than or equal to iwlen, then no compressions
C     occurred, and iw (maxmem+1 ... iwlen) was unused.  Otherwise
C     compressions did occur, and iwlen would have had to have been
C     greater than or equal to maxmem for no compressions to occur.
C     Return the value of maxmem in the pfree argument.

       RJNFO(1) = OPS
       RJNFO(2) = NRLADU
       JNFO(1) = NCMPA
       JNFO(2) = RSTRT
       PFREE = MAXMEM

       RETURN
       END
C ====================================================================
C ====================================================================
C ====================================================================
! COPYRIGHT (c) 1999 Council for the Central Laboratory
!                    of the Research Councils
! Original date July 1999
! AUTHORS Iain Duff (i.duff@rl.ac.uk) and
!         Jacko Koster (jacko.koster@uninett.no)
!
! Version 1.6.0
! See ChangeLog for version history.
!
      SUBROUTINE MC64ID(ICNTL)
      IMPLICIT NONE

C  Purpose
C  =======
C
C  The components of the array ICNTL control the action of MC64A/AD.
C  Default values for these are set in this subroutine.
C
C  Parameters
C  ==========
C
      INTEGER ICNTL(10)
C
C  Local variables
      INTEGER I
C
C    ICNTL(1) has default value 6.
C     It is the output stream for error messages. If it
C     is negative, these messages will be suppressed.
C
C    ICNTL(2) has default value 6.
C     It is the output stream for warning messages.
C     If it is negative, these messages are suppressed.
C
C    ICNTL(3) has default value -1.
C     It is the output stream for monitoring printing.
C     If it is negative, these messages are suppressed.
C
C    ICNTL(4) has default value 0.
C     If left at the defaut value, the incoming data is checked for
C     out-of-range indices and duplicates.  Setting ICNTL(4) to any
C     other will avoid the checks but is likely to cause problems
C     later if out-of-range indices or duplicates are present.
C     The user should only set ICNTL(4) non-zero, if the data is
C     known to avoid these problems.
C
C    ICNTL(5) to ICNTL(10) are not used by MC64A/AD but are set to
C     zero in this routine.

C Initialization of the ICNTL array.
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = -1
      DO 10 I = 4,10
        ICNTL(I) = 0
   10 CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64AD(JOB,N,NE,IP,IRN,A,NUM,CPERM,LIW,IW,LDW,DW,
     &           ICNTL,INFO)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
C  Purpose
C  =======
C
C This subroutine attempts to find a column permutation for an NxN
C sparse matrix A = {a_ij} that makes the permuted matrix have N
C entries on its diagonal.
C If the matrix is structurally nonsingular, the subroutine optionally
C returns a column permutation that maximizes the smallest element
C on the diagonal, maximizes the sum of the diagonal entries, or
C maximizes the product of the diagonal entries of the permuted matrix.
C For the latter option, the subroutine also finds scaling factors
C that may be used to scale the matrix so that the nonzero diagonal
C entries of the permuted matrix are one in absolute value and all the
C off-diagonal entries are less than or equal to one in absolute value.
C The natural logarithms of the scaling factors u(i), i=1..N, for the
C rows and v(j), j=1..N, for the columns are returned so that the
C scaled matrix B = {b_ij} has entries b_ij = a_ij * EXP(u_i + v_j).
C
C  Parameters
C  ==========
C
      INTEGER JOB,N,NE,NUM,LIW,LDW
      INTEGER IP(N+1),IRN(NE),CPERM(N),IW(LIW),ICNTL(10),INFO(10)
      DOUBLE PRECISION A(NE),DW(LDW)
C
C JOB is an INTEGER variable which must be set by the user to
C control the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   1 Compute a column permutation of the matrix so that the
C     permuted matrix has as many entries on its diagonal as possible.
C     The values on the diagonal are of arbitrary size. HSL subroutine
C     MC21A/AD is used for this. See [1].
C   2 Compute a column permutation of the matrix so that the smallest
C     value on the diagonal of the permuted matrix is maximized.
C     See [3].
C   3 Compute a column permutation of the matrix so that the smallest
C     value on the diagonal of the permuted matrix is maximized.
C     The algorithm differs from the one used for JOB = 2 and may
C     have quite a different performance. See [2].
C   4 Compute a column permutation of the matrix so that the sum
C     of the diagonal entries of the permuted matrix is maximized.
C     See [3].
C   5 Compute a column permutation of the matrix so that the product
C     of the diagonal entries of the permuted matrix is maximized
C     and vectors to scale the matrix so that the nonzero diagonal
C     entries of the permuted matrix are one in absolute value and
C     all the off-diagonal entries are less than or equal to one in
C     absolute value. See [3].
C  Restriction: 1 <= JOB <= 5.
C
C N is an INTEGER variable which must be set by the user to the
C   order of the matrix A. It is not altered by the subroutine.
C   Restriction: N >= 1.
C
C NE is an INTEGER variable which must be set by the user to the
C   number of entries in the matrix. It is not altered by the
C   subroutine.
C   Restriction: NE >= 1.
C
C IP is an INTEGER array of length N+1.
C   IP(J), J=1..N, must be set by the user to the position in array IRN
C   of the first row index of an entry in column J. IP(N+1) must be set
C   to NE+1. It is not altered by the subroutine.
C
C IRN is an INTEGER array of length NE.
C   IRN(K), K=1..NE, must be set by the user to hold the row indices of
C   the entries of the matrix. Those belonging to column J must be
C   stored contiguously in the positions IP(J)..IP(J+1)-1. The ordering
C   of the row indices within each column is unimportant. Repeated
C   entries are not allowed. The array IRN is not altered by the
C   subroutine.
C
C A is a DOUBLE PRECISION array of length NE.
C   The user must set A(K), K=1..NE, to the numerical value of the
C   entry that corresponds to IRN(K).
C   It is not used by the subroutine when JOB = 1.
C   It is not altered by the subroutine.
C
C NUM is an INTEGER variable that need not be set by the user.
C   On successful exit, NUM will be the number of entries on the
C   diagonal of the permuted matrix.
C   If NUM < N, the matrix is structurally singular.
C
C CPERM is an INTEGER array of length N that need not be set by the
C   user. On successful exit, CPERM contains the column permutation.
C   Column ABS(CPERM(J)) of the original matrix is column J in the
C   permuted matrix, J=1..N. For the N-NUM  entries of CPERM that are
C   not matched the permutation array is set negative so that a full
C   permutation of the matrix is obtained even in the structurally
C   singular case.
C
C LIW is an INTEGER variable that must be set by the user to
C   the dimension of array IW. It is not altered by the subroutine.
C   Restriction:
C     JOB = 1 :  LIW >= 5N
C     JOB = 2 :  LIW >= 4N
C     JOB = 3 :  LIW >= 10N + NE
C     JOB = 4 :  LIW >= 5N
C     JOB = 5 :  LIW >= 5N
C
C IW is an INTEGER array of length LIW that is used for workspace.
C
C LDW is an INTEGER variable that must be set by the user to the
C   dimension of array DW. It is not altered by the subroutine.
C   Restriction:
C     JOB = 1 :  LDW is not used
C     JOB = 2 :  LDW >= N
C     JOB = 3 :  LDW >= NE
C     JOB = 4 :  LDW >= 2N + NE
C     JOB = 5 :  LDW >= 3N + NE
C
C DW is a DOUBLE PRECISION array of length LDW
C   that is used for workspace. If JOB = 5, on return,
C   DW(i) contains u_i, i=1..N, and DW(N+j) contains v_j, j=1..N.
C
C ICNTL is an INTEGER array of length 10. Its components control the
C   output of MC64A/AD and must be set by the user before calling
C   MC64A/AD. They are not altered by the subroutine.
C
C   ICNTL(1) must be set to specify the output stream for
C   error messages. If ICNTL(1) < 0, messages are suppressed.
C   The default value set by MC46I/ID is 6.
C
C   ICNTL(2) must be set by the user to specify the output stream for
C   warning messages. If ICNTL(2) < 0, messages are suppressed.
C   The default value set by MC46I/ID is 6.
C
C   ICNTL(3) must be set by the user to specify the output stream for
C   diagnostic messages. If ICNTL(3) < 0, messages are suppressed.
C   The default value set by MC46I/ID is -1.
C
C   ICNTL(4) must be set by the user to a value other than 0 to avoid
C   checking of the input data.
C   The default value set by MC46I/ID is 0.
C
C INFO is an INTEGER array of length 10 which need not be set by the
C   user. INFO(1) is set non-negative to indicate success. A negative
C   value is returned if an error occurred, a positive value if a
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : successful entry (for structurally singular matrix).
C   +2 : the returned scaling factors are large and may cause
C        overflow when used to scale the matrix.
C        (For JOB = 5 entry only.)
C   -1 : JOB < 1 or JOB > 5.  Value of JOB held in INFO(2).
C   -2 : N < 1.  Value of N held in INFO(2).
C   -3 : NE < 1. Value of NE held in INFO(2).
C   -4 : the defined length LIW violates the restriction on LIW.
C        Value of LIW required given by INFO(2).
C   -5 : the defined length LDW violates the restriction on LDW.
C        Value of LDW required given by INFO(2).
C   -6 : entries are found whose row indices are out of range. INFO(2)
C        contains the index of a column in which such an entry is found.
C   -7 : repeated entries are found. INFO(2) contains the index of a
C        column in which such entries are found.
C  INFO(3) to INFO(10) are not currently used and are set to zero by
C        the routine.
C
C References:
C  [1]  I. S. Duff, (1981),
C       "Algorithm 575. Permutations for a zero-free diagonal",
C       ACM Trans. Math. Software 7(3), 387-390.
C  [2]  I. S. Duff and J. Koster, (1998),
C       "The design and use of algorithms for permuting large
C       entries to the diagonal of sparse matrices",
C       SIAM J. Matrix Anal. Appl., vol. 20, no. 4, pp. 889-901.
C  [3]  I. S. Duff and J. Koster, (1999),
C       "On algorithms for permuting large entries to the diagonal
C       of sparse matrices",
C       Technical Report RAL-TR-1999-030, RAL, Oxfordshire, England.

C Local variables and parameters
      INTEGER I,J,K
      DOUBLE PRECISION FACT,ZERO,RINF
      PARAMETER (ZERO=0.0D+00)
C External routines and functions
      EXTERNAL MC21AD,MC64BD,MC64RD,MC64SD,MC64WD
C Intrinsic functions
      INTRINSIC ABS,LOG

C Set RINF to largest positive real number (infinity)
      RINF = HUGE(RINF)

C Check value of JOB
      IF (JOB.LT.1 .OR. JOB.GT.5) THEN
        INFO(1) = -1
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
C Check value of N
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
C Check value of NE
      IF (NE.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NE
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NE',NE
        GO TO 99
      ENDIF
C Check LIW
      IF (JOB.EQ.1) K = 5*N
      IF (JOB.EQ.2) K = 4*N
      IF (JOB.EQ.3) K = 10*N + NE
      IF (JOB.EQ.4) K = 5*N
      IF (JOB.EQ.5) K = 5*N
      IF (LIW.LT.K) THEN
        INFO(1) = -4
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
C Check LDW
C If JOB = 1, do not check
      IF (JOB.GT.1) THEN
        IF (JOB.EQ.2) K = N
        IF (JOB.EQ.3) K = NE
        IF (JOB.EQ.4) K = 2*N + NE
        IF (JOB.EQ.5) K = 3*N + NE
        IF (LDW.LT.K) THEN
          INFO(1) = -5
          INFO(2) = K
          IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
          GO TO 99
        ENDIF
      ENDIF
      IF (ICNTL(4).EQ.0) THEN
C Check row indices. Use IW(1:N) as workspace
        DO 3 I = 1,N
          IW(I) = 0
    3   CONTINUE
        DO 6 J = 1,N
          DO 4 K = IP(J),IP(J+1)-1
            I = IRN(K)
C Check for row indices that are out of range
            IF (I.LT.1 .OR. I.GT.N) THEN
              INFO(1) = -6
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
C Check for repeated row indices within a column
            IF (IW(I).EQ.J) THEN
              INFO(1) = -7
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
    4     CONTINUE
    6   CONTINUE
      ENDIF

C Print diagnostics on input
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,N,NE
        WRITE(ICNTL(3),9021) (IP(J),J=1,N+1)
        WRITE(ICNTL(3),9022) (IRN(J),J=1,NE)
        IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(J),J=1,NE)
      ENDIF

C Set components of INFO to zero
      DO 8 I=1,10
        INFO(I) = 0
    8 CONTINUE

C Compute maximum matching with MC21A/AD
      IF (JOB.EQ.1) THEN
C Put length of column J in IW(J)
        DO 10 J = 1,N
          IW(J) = IP(J+1) - IP(J)
   10   CONTINUE
C IW(N+1:5N) is workspace
        CALL MC21AD(N,IRN,NE,IP,IW(1),CPERM,NUM,IW(N+1))
        GO TO 90
      ENDIF

C Compute bottleneck matching
      IF (JOB.EQ.2) THEN
C IW(1:5N), DW(1:N) are workspaces
        CALL MC64BD(N,NE,IP,IRN,A,CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),DW)
        GO TO 90
      ENDIF

C Compute bottleneck matching
      IF (JOB.EQ.3) THEN
C Copy IRN(K) into IW(K), ABS(A(K)) into DW(K), K=1..NE
        DO 20 K = 1,NE
          IW(K) = IRN(K)
          DW(K) = ABS(A(K))
   20   CONTINUE
C Sort entries in each column by decreasing value.
        CALL MC64RD(N,NE,IP,IW,DW)
C IW(NE+1:NE+10N) is workspace
        CALL MC64SD(N,NE,IP,IW(1),DW,CPERM,NUM,IW(NE+1),
     &     IW(NE+N+1),IW(NE+2*N+1),IW(NE+3*N+1),IW(NE+4*N+1),
     &     IW(NE+5*N+1),IW(NE+6*N+1))
        GO TO 90
      ENDIF

      IF (JOB.EQ.4) THEN
        DO 50 J = 1,N
          FACT = ZERO
          DO 30 K = IP(J),IP(J+1)-1
            IF (ABS(A(K)).GT.FACT) FACT = ABS(A(K))
   30     CONTINUE
          DO 40 K = IP(J),IP(J+1)-1
            DW(2*N+K) = FACT - ABS(A(K))
   40     CONTINUE
   50   CONTINUE
C B = DW(2N+1:2N+NE); IW(1:5N) and DW(1:2N) are workspaces
        CALL MC64WD(N,NE,IP,IRN,DW(2*N+1),CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
     &     DW(1),DW(N+1))
        GO TO 90
      ENDIF

      IF (JOB.EQ.5) THEN
        DO 75 J = 1,N
          FACT = ZERO
          DO 60 K = IP(J),IP(J+1)-1
            DW(3*N+K) = ABS(A(K))
            IF (DW(3*N+K).GT.FACT) FACT = DW(3*N+K)
   60     CONTINUE
          DW(2*N+J) = FACT
          IF (FACT.NE.ZERO) THEN
            FACT = LOG(FACT)
          ELSE
            FACT = RINF/N
          ENDIF
          DO 70 K = IP(J),IP(J+1)-1
            IF (DW(3*N+K).NE.ZERO) THEN
              DW(3*N+K) = FACT - LOG(DW(3*N+K))
            ELSE
              DW(3*N+K) = RINF/N
            ENDIF
   70     CONTINUE
   75   CONTINUE
C B = DW(3N+1:3N+NE); IW(1:5N) and DW(1:2N) are workspaces
        CALL MC64WD(N,NE,IP,IRN,DW(3*N+1),CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
     &     DW(1),DW(N+1))
        IF (NUM.EQ.N) THEN
          DO 80 J = 1,N
            IF (DW(2*N+J).NE.ZERO) THEN
              DW(N+J) = DW(N+J) - LOG(DW(2*N+J))
            ELSE
              DW(N+J) = ZERO
            ENDIF
   80     CONTINUE
        ENDIF
C Check size of scaling factors
        FACT = 0.5*LOG(RINF)
        DO 86 J = 1,N
          IF (DW(J).LT.FACT .AND. DW(N+J).LT.FACT) GO TO 86
          INFO(1) = 2
          GO TO 90
   86   CONTINUE
C       GO TO 90
      ENDIF

   90 IF (INFO(1).EQ.0 .AND. NUM.LT.N) THEN
C Matrix is structurally singular, return with warning
        INFO(1) = 1
        IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9011) INFO(1)
      ENDIF
      IF (INFO(1).EQ.2) THEN
C Scaling factors are large, return with warning
        IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9012) INFO(1)
      ENDIF

C Print diagnostics on output
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,2)
        WRITE(ICNTL(3),9031) NUM
        WRITE(ICNTL(3),9032) (CPERM(J),J=1,N)
        IF (JOB.EQ.5) THEN
          WRITE(ICNTL(3),9033) (DW(J),J=1,N)
          WRITE(ICNTL(3),9034) (DW(N+J),J=1,N)
        ENDIF
      ENDIF

C Return from subroutine.
   99 RETURN

 9001 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2,
     &        ' because ',(A),' = ',I10)
 9004 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains two or more entries with row index ',I8)
 9011 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
     &        '        The matrix is structurally singular.')
 9012 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
     &        '        Some scaling factors may be too large.')
 9020 FORMAT (' ****** Input parameters for MC64A/AD:'/
     &        ' JOB = ',I8/' N   = ',I8/' NE  = ',I8)
 9021 FORMAT (' IP(1:N+1)  = ',8I8/(14X,8I8))
 9022 FORMAT (' IRN(1:NE)  = ',8I8/(14X,8I8))
 9023 FORMAT (' A(1:NE)    = ',4(1PD14.4)/(14X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC64A/AD:'/
     &        ' INFO(1:2)  = ',2I8)
 9031 FORMAT (' NUM        = ',I8)
 9032 FORMAT (' CPERM(1:N) = ',8I8/(14X,8I8))
 9033 FORMAT (' DW(1:N)    = ',5(F11.3)/(14X,5(F11.3)))
 9034 FORMAT (' DW(N+1:2N) = ',5(F11.3)/(14X,5(F11.3)))
      END

C**********************************************************************
      SUBROUTINE MC64BD(N,NE,IP,IRN,A,IPERM,NUM,JPERM,PR,Q,L,D)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(N),JPERM(N),PR(N),Q(N),L(N)
      DOUBLE PRECISION A(NE),D(N)

C N, NE, IP, IRN are described in MC64A/AD.
C A is a DOUBLE PRECISION array of length
C   NE. A(K), K=1..NE, must be set to the value of the entry
C   that corresponds to IRN(K). It is not altered.
C IPERM is an INTEGER array of length N. On exit, it contains the
C    matching: IPERM(I) = 0 or row I is matched to column IPERM(I).
C NUM is INTEGER variable. On exit, it contains the cardinality of the
C    matching stored in IPERM.
C IW is an INTEGER work array of length 4N.
C DW is a DOUBLE PRECISION array of length N.

C Local variables
      INTEGER I,II,J,JJ,JORD,Q0,QLEN,IDUM,JDUM,ISP,JSP,
     &        K,KK,KK1,KK2,I0,UP,LOW,LPOS
      DOUBLE PRECISION CSP,DI,DNEW,DQ0,AI,A0,BV
C Local parameters
      DOUBLE PRECISION RINF,ZERO,MINONE
      PARAMETER (ZERO=0.0D+0,MINONE=-1.0D+0)
C Intrinsic functions
      INTRINSIC ABS,MIN
C External subroutines and/or functions
      EXTERNAL MC64DD,MC64ED,MC64FD


C Set RINF to largest positive real number
      RINF = HUGE(RINF)

C Initialization
      NUM = 0
      BV = RINF
      DO 10 K = 1,N
        IPERM(K) = 0
        JPERM(K) = 0
        PR(K) = IP(K)
        D(K) = ZERO
   10 CONTINUE
C Scan columns of matrix;
      DO 20 J = 1,N
        A0 = MINONE
        DO 30 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.GT.D(I)) D(I) = AI
          IF (JPERM(J).NE.0) GO TO 30
          IF (AI.GE.BV) THEN
            A0 = BV
            IF (IPERM(I).NE.0) GO TO 30
            JPERM(J) = I
            IPERM(I) = J
            NUM = NUM + 1
          ELSE
            IF (AI.LE.A0) GO TO 30
            A0 = AI
            I0 = I
          ENDIF
   30   CONTINUE
        IF (A0.NE.MINONE .AND. A0.LT.BV) THEN
          BV = A0
          IF (IPERM(I0).NE.0) GO TO 20
          IPERM(I0) = J
          JPERM(J) = I0
          NUM = NUM + 1
        ENDIF
   20 CONTINUE
C Update BV with smallest of all the largest maximum absolute values
C of the rows.
      DO 25 I = 1,N
        BV = MIN(BV,D(I))
   25 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
C Rescan unassigned columns; improve initial assignment
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        DO 50 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.LT.BV) GO TO 50
          IF (IPERM(I).EQ.0) GO TO 90
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 50
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).NE.0) GO TO 70
            IF (ABS(A(KK)).GE.BV) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   50   CONTINUE
        GO TO 95
   80   JPERM(JJ) = II
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = I
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000

C Prepare for main loop
      DO 99 I = 1,N
        D(I) = MINONE
        L(I) = 0
   99 CONTINUE

C Main loop ... each pass round this loop is similar to Dijkstra's
C algorithm for solving the single source shortest path problem

      DO 100 JORD = 1,N

        IF (JPERM(JORD).NE.0) GO TO 100
        QLEN = 0
        LOW = N + 1
        UP = N + 1
C CSP is cost of shortest path to any unassigned row
C ISP is matrix position of unassigned row element in shortest path
C JSP is column index of unassigned row element in shortest path
        CSP = MINONE
C Build shortest path tree starting from unassigned column JORD
        J = JORD
        PR(J) = -1

C Scan column J
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = ABS(A(K))
          IF (CSP.GE.DNEW) GO TO 115
          IF (IPERM(I).EQ.0) THEN
C Row I is unassigned; update shortest path info
            CSP = DNEW
            ISP = I
            JSP = J
            IF (CSP.GE.BV) GO TO 160
          ELSE
            D(I) = DNEW
            IF (DNEW.GE.BV) THEN
C Add row I to Q2
              LOW = LOW - 1
              Q(LOW) = I
            ELSE
C Add row I to Q, and push it
              QLEN = QLEN + 1
              L(I) = QLEN
              CALL MC64DD(I,N,Q,D,L,1)
            ENDIF
            JJ = IPERM(I)
            PR(JJ) = J
          ENDIF
  115   CONTINUE

        DO 150 JDUM = 1,NUM
C If Q2 is empty, extract new rows from Q
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (CSP.GE.D(I)) GO TO 160
            BV = D(I)
            DO 152 IDUM = 1,N
              CALL MC64ED(QLEN,N,Q,D,L,1)
              L(I) = 0
              LOW = LOW - 1
              Q(LOW) = I
              IF (QLEN.EQ.0) GO TO 153
              I = Q(1)
              IF (D(I).NE.BV) GO TO 153
  152       CONTINUE
C End of dummy loop; this point is never reached
          ENDIF
C Move row Q0
  153     UP = UP - 1
          Q0 = Q(UP)
          DQ0 = D(Q0)
          L(Q0) = UP
C Scan column that matches with row Q0
          J = IPERM(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
C Update D(I)
            IF (L(I).GE.UP) GO TO 155
            DNEW = MIN(DQ0,ABS(A(K)))
            IF (CSP.GE.DNEW) GO TO 155
            IF (IPERM(I).EQ.0) THEN
C Row I is unassigned; update shortest path info
              CSP = DNEW
              ISP = I
              JSP = J
              IF (CSP.GE.BV) GO TO 160
            ELSE
              DI = D(I)
              IF (DI.GE.BV .OR. DI.GE.DNEW) GO TO 155
              D(I) = DNEW
              IF (DNEW.GE.BV) THEN
C Delete row I from Q (if necessary); add row I to Q2
                IF (DI.NE.MINONE) THEN
                  LPOS = L(I)
                  CALL MC64FD(LPOS,QLEN,N,Q,D,L,1)
                ENDIF
                L(I) = 0
                LOW = LOW - 1
                Q(LOW) = I
              ELSE
C Add row I to Q (if necessary); push row I up Q
                IF (DI.EQ.MINONE) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64DD(I,N,Q,D,L,1)
              ENDIF
C Update tree
              JJ = IPERM(I)
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE

C If CSP = MINONE, no augmenting path is found
  160   IF (CSP.EQ.MINONE) GO TO 190
C Update bottleneck value
        BV = MIN(BV,CSP)
C Find augmenting path by tracing backward in PR; update IPERM,JPERM
        NUM = NUM + 1
        I = ISP
        J = JSP
        DO 170 JDUM = 1,NUM+1
          I0 = JPERM(J)
          JPERM(J) = I
          IPERM(I) = J
          J = PR(J)
          IF (J.EQ.-1) GO TO 190
          I = I0
  170   CONTINUE
C End of dummy loop; this point is never reached
  190   DO 191 KK = UP,N
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  191   CONTINUE
        DO 192 KK = LOW,UP-1
          I = Q(KK)
          D(I) = MINONE
  192   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  193   CONTINUE

  100 CONTINUE
C End of main loop

C BV is bottleneck value of final matching
      IF (NUM.EQ.N) GO TO 1000

C Matrix is structurally singular, complete IPERM.
C JPERM, PR are work arrays
      DO 300 J = 1,N
        JPERM(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          PR(K) = I
        ELSE
          J = IPERM(I)
          JPERM(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 I = 1,N
        IF (JPERM(I).NE.0) GO TO 320
        K = K + 1
        JDUM = PR(K)
        IPERM(JDUM) = - I
  320 CONTINUE

 1000 RETURN
      END

C**********************************************************************
      SUBROUTINE MC64DD(I,N,Q,D,L,IWAY)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER I,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)

C Variables N,Q,D,L are described in MC64B/BD
C IF IWAY is equal to 1, then
C node I is pushed from its current position upwards
C IF IWAY is not equal to 1, then
C node I is pushed from its current position downwards

C Local variables and parameters
      INTEGER IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DI

      POS = L(I)
      IF (POS.LE.1) GO TO 20
      DI = D(I)
C POS is index of current position of I in the tree
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
C End of dummy loop; this point is never reached
      ELSE
        DO 15 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   15   CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C End of dummy if; this point is never reached
   20 Q(POS) = I
      L(I) = POS

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64ED(QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)

C Variables QLEN,N,Q,D,L are described in MC64B/BD (IWAY = 1) or
C     MC64W/WD (IWAY = 2)
C The root node is deleted from the binary heap.

C Local variables and parameters
      INTEGER I,IDUM,K,POS,POSK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI

C Move last element to begin of Q
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = 1
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 20
C Exchange old last element with larger priority child
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   10   CONTINUE
C End of dummy loop; this point is never reached
      ELSE
        DO 15 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 20
C Exchange old last element with smaller child
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   15   CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C End of dummy if; this point is never reached
   20 Q(POS) = I
      L(I) = POS

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64FD(POS0,QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER POS0,QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)

C Variables QLEN,N,Q,D,L are described in MC64B/BD (IWAY = 1) or
C     MC64WD (IWAY = 2).
C Move last element in the heap

      INTEGER I,IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI

C Quick return, if possible
      IF (QLEN.EQ.POS0) THEN
        QLEN = QLEN - 1
        RETURN
      ENDIF

C Move last element from queue Q to position POS0
C POS is current position of node I in the tree
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = POS0
      IF (IWAY.EQ.1) THEN
        IF (POS.LE.1) GO TO 20
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
C End of dummy loop; this point is never reached
   20   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 30 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   30   CONTINUE
C End of dummy loop; this point is never reached
      ELSE
        IF (POS.LE.1) GO TO 34
        DO 32 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 34
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 34
   32   CONTINUE
C End of dummy loop; this point is never reached
   34   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 36 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   36   CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C End of dummy if; this point is never reached
   40 Q(POS) = I
      L(I) = POS

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64RD(N,NE,IP,IRN,A)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER N,NE
      INTEGER IP(N+1),IRN(NE)
      DOUBLE PRECISION A(NE)

C This subroutine sorts the entries in each column of the
C sparse matrix (defined by N,NE,IP,IRN,A) by decreasing
C numerical value.

C Local constants
      INTEGER THRESH,TDLEN
      PARAMETER (THRESH=15,TDLEN=50)
C Local variables
      INTEGER J,IPJ,K,LEN,R,S,HI,FIRST,MID,LAST,TD
      DOUBLE PRECISION HA,KEY
C Local arrays
      INTEGER TODO(TDLEN)

      DO 100 J = 1,N
        LEN = IP(J+1) - IP(J)
        IF (LEN.LE.1) GO TO 100
        IPJ = IP(J)

C Sort array roughly with partial quicksort
        IF (LEN.LT.THRESH) GO TO 400
        TODO(1) = IPJ
        TODO(2) = IPJ + LEN
        TD = 2
  500   CONTINUE
        FIRST = TODO(TD-1)
        LAST = TODO(TD)
C KEY is the smallest of two values present in interval [FIRST,LAST)
        KEY = A((FIRST+LAST)/2)
        DO 475 K = FIRST,LAST-1
          HA = A(K)
          IF (HA.EQ.KEY) GO TO 475
          IF (HA.GT.KEY) GO TO 470
          KEY = HA
          GO TO 470
  475   CONTINUE
C Only one value found in interval, so it is already sorted
        TD = TD - 2
        GO TO 425

C Reorder interval [FIRST,LAST) such that entries before MID are gt KEY
  470   MID = FIRST
        DO 450 K = FIRST,LAST-1
          IF (A(K).LE.KEY) GO TO 450
          HA = A(MID)
          A(MID) = A(K)
          A(K) = HA
          HI = IRN(MID)
          IRN(MID) = IRN(K)
          IRN(K) = HI
          MID = MID + 1
  450   CONTINUE
C Both subintervals [FIRST,MID), [MID,LAST) are nonempty
C Stack the longest of the two subintervals first
        IF (MID-FIRST.GE.LAST-MID) THEN
          TODO(TD+2) = LAST
          TODO(TD+1) = MID
          TODO(TD) = MID
C          TODO(TD-1) = FIRST
        ELSE
          TODO(TD+2) = MID
          TODO(TD+1) = FIRST
          TODO(TD) = LAST
          TODO(TD-1) = MID
        ENDIF
        TD = TD + 2

  425   CONTINUE
        IF (TD.EQ.0) GO TO 400
C There is still work to be done
        IF (TODO(TD)-TODO(TD-1).GE.THRESH) GO TO 500
C Next interval is already short enough for straightforward insertion
        TD = TD - 2
        GO TO 425

C Complete sorting with straightforward insertion
  400   DO 200 R = IPJ+1,IPJ+LEN-1
          IF (A(R-1) .LT. A(R)) THEN
            HA = A(R)
            HI = IRN(R)
            A(R) = A(R-1)
            IRN(R) = IRN(R-1)
            DO 300 S = R-1,IPJ+1,-1
              IF (A(S-1) .LT. HA) THEN
                A(S) = A(S-1)
                IRN(S) = IRN(S-1)
              ELSE
                A(S) = HA
                IRN(S) = HI
                GO TO 200
              END IF
  300       CONTINUE
            A(IPJ) = HA
            IRN(IPJ) = HI
          END IF
  200   CONTINUE

  100 CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64SD(N,NE,IP,IRN,A,IPERM,NUMX,
     &           W,LEN,LENL,LENH,FC,IW,IW4)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER N,NE,NUMX
      INTEGER IP(N+1),IRN(NE),IPERM(N),
     &        W(N),LEN(N),LENL(N),LENH(N),FC(N),IW(N),IW4(4*N)
      DOUBLE PRECISION A(NE)

C N, NE, IP, IRN, are described in MC64A/AD.
C A is a DOUBLE PRECISION array of length NE.
C   A(K), K=1..NE, must be set to the value of the entry that
C   corresponds to IRN(k). The entries in each column must be
C   non-negative and ordered by decreasing value.
C IPERM is an INTEGER array of length N. On exit, it contains the
C   bottleneck matching: IPERM(I) - 0 or row I is matched to column
C   IPERM(I).
C NUMX is an INTEGER variable. On exit, it contains the cardinality
C   of the matching stored in IPERM.
C IW is an INTEGER work array of length 10N.

C FC is an integer array of length N that contains the list of
C   unmatched columns.
C LEN(J), LENL(J), LENH(J) are integer arrays of length N that point
C   to entries in matrix column J.
C   In the matrix defined by the column parts IP(J)+LENL(J) we know
C   a matching does not exist; in the matrix defined by the column
C   parts IP(J)+LENH(J) we know one exists.
C   LEN(J) lies between LENL(J) and LENH(J) and determines the matrix
C   that is tested for a maximum matching.
C W is an integer array of length N and contains the indices of the
C   columns for which LENL ne LENH.
C WLEN is number of indices stored in array W.
C IW is integer work array of length N.
C IW4 is integer work array of length 4N used by MC64U/UD.

      INTEGER NUM,NVAL,WLEN,II,I,J,K,L,CNT,MOD,IDUM1,IDUM2,IDUM3
      DOUBLE PRECISION BVAL,BMIN,BMAX,RINF
      EXTERNAL MC64QD,MC64UD

C BMIN and BMAX are such that a maximum matching exists for the input
C   matrix in which all entries smaller than BMIN are dropped.
C   For BMAX, a maximum matching does not exist.
C BVAL is a value between BMIN and BMAX.
C CNT is the number of calls made to MC64U/UD so far.
C NUM is the cardinality of last matching found.

C Set RINF to largest positive real number
      RINF = HUGE(RINF)

C Compute a first maximum matching from scratch on whole matrix.
      DO 20 J = 1,N
        FC(J) = J
        IW(J) = 0
        LEN(J) = IP(J+1) - IP(J)
   20 CONTINUE
C The first call to MC64U/UD
      CNT = 1
      MOD = 1
      NUMX = 0
      CALL MC64UD(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUMX,N,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))

C IW contains a maximum matching of length NUMX.
      NUM = NUMX

      IF (NUM.NE.N) THEN
C Matrix is structurally singular
        BMAX = RINF
      ELSE
C Matrix is structurally nonsingular, NUM=NUMX=N;
C Set BMAX just above the smallest of all the maximum absolute
C values of the columns
        BMAX = RINF
        DO 30 J = 1,N
          BVAL = 0.0
          DO 25 K = IP(J),IP(J+1)-1
            IF (A(K).GT.BVAL) BVAL = A(K)
   25     CONTINUE
          IF (BVAL.LT.BMAX) BMAX = BVAL
   30   CONTINUE
        BMAX = 1.001 * BMAX
      ENDIF

C Initialize BVAL,BMIN
      BVAL = 0.0
      BMIN = 0.0
C Initialize LENL,LEN,LENH,W,WLEN according to BMAX.
C Set LEN(J), LENH(J) just after last entry in column J.
C Set LENL(J) just after last entry in column J with value ge BMAX.
      WLEN = 0
      DO 48 J = 1,N
        L = IP(J+1) - IP(J)
        LENH(J) = L
        LEN(J) = L
        DO 45 K = IP(J),IP(J+1)-1
          IF (A(K).LT.BMAX) GO TO 46
   45   CONTINUE
C Column J is empty or all entries are ge BMAX
        K = IP(J+1)
   46   LENL(J) = K - IP(J)
C Add J to W if LENL(J) ne LENH(J)
        IF (LENL(J).EQ.L) GO TO 48
        WLEN = WLEN + 1
        W(WLEN) = J
   48 CONTINUE

C Main loop
      DO 90 IDUM1 = 1,NE
        IF (NUM.EQ.NUMX) THEN
C We have a maximum matching in IW; store IW in IPERM
          DO 50 I = 1,N
            IPERM(I) = IW(I)
   50     CONTINUE
C Keep going round this loop until matching IW is no longer maximum.
          DO 80 IDUM2 = 1,NE
            BMIN = BVAL
            IF (BMAX .EQ. BMIN) GO TO 99
C Find splitting value BVAL
            CALL MC64QD(IP,LENL,LEN,W,WLEN,A,NVAL,BVAL)
            IF (NVAL.LE.1) GO TO 99
C Set LEN such that all matrix entries with value lt BVAL are
C discarded. Store old LEN in LENH. Do this for all columns W(K).
C Each step, either K is incremented or WLEN is decremented.
            K = 1
            DO 70 IDUM3 = 1,N
              IF (K.GT.WLEN) GO TO 71
              J = W(K)
              DO 55 II = IP(J)+LEN(J)-1,IP(J)+LENL(J),-1
                IF (A(II).GE.BVAL) GO TO 60
                I = IRN(II)
                IF (IW(I).NE.J) GO TO 55
C Remove entry from matching
                IW(I) = 0
                NUM = NUM - 1
                FC(N-NUM) = J
   55         CONTINUE
   60         LENH(J) = LEN(J)
C IP(J)+LEN(J)-1 is last entry in column ge BVAL
              LEN(J) = II - IP(J) + 1
C If LENH(J) = LENL(J), remove J from W
              IF (LENL(J).EQ.LENH(J)) THEN
                W(K) = W(WLEN)
                WLEN = WLEN - 1
              ELSE
                K = K + 1
              ENDIF
   70       CONTINUE
   71       IF (NUM.LT.NUMX) GO TO 81
   80     CONTINUE
C End of dummy loop; this point is never reached
C Set mode for next call to MC64U/UD
   81     MOD = 1
        ELSE
C We do not have a maximum matching in IW.
          BMAX = BVAL
C BMIN is the bottleneck value of a maximum matching;
C for BMAX the matching is not maximum, so BMAX>BMIN
C          IF (BMAX .EQ. BMIN) GO TO 99
C Find splitting value BVAL
          CALL MC64QD(IP,LEN,LENH,W,WLEN,A,NVAL,BVAL)
          IF (NVAL.EQ.0. OR. BVAL.EQ.BMIN) GO TO 99
C Set LEN such that all matrix entries with value ge BVAL are
C inside matrix. Store old LEN in LENL. Do this for all columns W(K).
C Each step, either K is incremented or WLEN is decremented.
          K = 1
          DO 87 IDUM3 = 1,N
            IF (K.GT.WLEN) GO TO 88
            J = W(K)
            DO 85 II = IP(J)+LEN(J),IP(J)+LENH(J)-1
              IF (A(II).LT.BVAL) GO TO 86
   85       CONTINUE
   86       LENL(J) = LEN(J)
            LEN(J) = II - IP(J)
            IF (LENL(J).EQ.LENH(J)) THEN
              W(K) = W(WLEN)
              WLEN = WLEN - 1
            ELSE
              K = K + 1
            ENDIF
   87     CONTINUE
C End of dummy loop; this point is never reached
C Set mode for next call to MC64U/UD
   88     MOD = 0
        ENDIF
        CNT = CNT + 1
        CALL MC64UD(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUM,NUMX,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))

C IW contains maximum matching of length NUM
   90 CONTINUE
C End of dummy loop; this point is never reached

C BMIN is bottleneck value of final matching
   99 IF (NUMX.EQ.N) GO TO 1000
C The matrix is structurally singular, complete IPERM
C W, IW are work arrays
      DO 300 J = 1,N
        W(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          IW(K) = I
        ELSE
          J = IPERM(I)
          W(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 J = 1,N
        IF (W(J).NE.0) GO TO 320
        K = K + 1
        IDUM1 = IW(K)
        IPERM(IDUM1) =  - J
  320 CONTINUE

 1000 RETURN
      END

C**********************************************************************
      SUBROUTINE MC64QD(IP,LENL,LENH,W,WLEN,A,NVAL,VAL)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER WLEN,NVAL
      INTEGER IP(*),LENL(*),LENH(*),W(*)
      DOUBLE PRECISION A(*),VAL

C This routine searches for at most XX different numerical values
C in the columns W(1:WLEN). XX>=2.
C Each column J is scanned between IP(J)+LENL(J) and IP(J)+LENH(J)-1
C until XX values are found or all columns have been considered.
C On output, NVAL is the number of different values that is found
C and SPLIT(1:NVAL) contains the values in decreasing order.
C If NVAL > 0, the routine returns VAL = SPLIT((NVAL+1)/2).
C
      INTEGER XX,J,K,II,S,POS
      PARAMETER (XX=10)
      DOUBLE PRECISION SPLIT(XX),HA

C Scan columns in W(1:WLEN). For each encountered value, if value not
C already present in SPLIT(1:NVAL), insert value such that SPLIT
C remains sorted by decreasing value.
C The sorting is done by straightforward insertion; therefore the use
C of this routine should be avoided for large XX (XX < 20).
      NVAL = 0
      DO 10 K = 1,WLEN
        J = W(K)
        DO 15 II = IP(J)+LENL(J),IP(J)+LENH(J)-1
          HA = A(II)
          IF (NVAL.EQ.0) THEN
            SPLIT(1) = HA
            NVAL = 1
          ELSE
C Check presence of HA in SPLIT
            DO 20 S = NVAL,1,-1
              IF (SPLIT(S).EQ.HA) GO TO 15
              IF (SPLIT(S).GT.HA) THEN
                POS = S + 1
                GO TO 21
              ENDIF
  20        CONTINUE
            POS = 1
C The insertion
  21        DO 22 S = NVAL,POS,-1
              SPLIT(S+1) = SPLIT(S)
  22        CONTINUE
            SPLIT(POS) = HA
            NVAL = NVAL + 1
          ENDIF
C Exit loop if XX values are found
          IF (NVAL.EQ.XX) GO TO 11
  15    CONTINUE
  10  CONTINUE
C Determine VAL
  11  IF (NVAL.GT.0) VAL = SPLIT((NVAL+1)/2)

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64UD(ID,MOD,N,IRN,LIRN,IP,LENC,FC,IPERM,NUM,NUMX,
     &           PR,ARP,CV,OUT)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER ID,MOD,N,LIRN,NUM,NUMX
      INTEGER ARP(N),CV(N),IRN(LIRN),IP(N),
     &        FC(N),IPERM(N),LENC(N),OUT(N),PR(N)

C PR(J) is the previous column to J in the depth first search.
C   Array PR is used as workspace in the sorting algorithm.
C Elements (I,IPERM(I)) I=1,..,N are entries at the end of the
C   algorithm unless N assignments have not been made in which case
C   N-NUM pairs (I,IPERM(I)) will not be entries in the matrix.
C CV(I) is the most recent loop number (ID+JORD) at which row I
C   was visited.
C ARP(J) is the number of entries in column J which have been scanned
C   when looking for a cheap assignment.
C OUT(J) is one less than the number of entries in column J which have
C   not been scanned during one pass through the main loop.
C NUMX is maximum possible size of matching.

      INTEGER I,II,IN1,IN2,J,J1,JORD,K,KK,LAST,NFC,
     &        NUM0,NUM1,NUM2,ID0,ID1

      IF (ID.EQ.1) THEN
C The first call to MC64U/UD.
C Initialize CV and ARP; parameters MOD, NUMX are not accessed
        DO 5 I = 1,N
          CV(I) = 0
          ARP(I) = 0
    5   CONTINUE
        NUM1 = N
        NUM2 = N
      ELSE
C Not the first call to MC64U/UD.
C Re-initialize ARP if entries were deleted since last call to MC64U/UD
        IF (MOD.EQ.1) THEN
          DO 8 I = 1,N
            ARP(I) = 0
    8     CONTINUE
        ENDIF
        NUM1 = NUMX
        NUM2 = N - NUMX
      ENDIF
      NUM0 = NUM

C NUM0 is size of input matching
C NUM1 is maximum possible size of matching
C NUM2 is maximum allowed number of unassigned rows/columns
C NUM is size of current matching

C Quick return if possible
C      IF (NUM.EQ.N) GO TO 199
C NFC is number of rows/columns that could not be assigned
      NFC = 0
C Integers ID0+1 to ID0+N are unique numbers for call ID to MC64U/UD,
C so 1st call uses 1..N, 2nd call uses N+1..2N, etc
      ID0 = (ID-1)*N

C Main loop. Each pass round this loop either results in a new
C assignment or gives a column with no assignment

      DO 100 JORD = NUM0+1,N

C Each pass uses unique number ID1
        ID1 = ID0 + JORD
C J is unmatched column
        J = FC(JORD-NUM0)
        PR(J) = -1
        DO 70 K = 1,JORD
C Look for a cheap assignment
          IF (ARP(J).GE.LENC(J)) GO TO 30
          IN1 = IP(J) + ARP(J)
          IN2 = IP(J) + LENC(J) - 1
          DO 20 II = IN1,IN2
            I = IRN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
C No cheap assignment in row
          ARP(J) = LENC(J)
C Begin looking for assignment chain starting with row J
   30     OUT(J) = LENC(J) - 1
C Inner loop.  Extends chain by one or backtracks
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENC(J) - 1
            IN1 = IN2 - IN1
C Forward scan
            DO 40 II = IN1,IN2
              I = IRN(II)
              IF (CV(I).EQ.ID1) GO TO 40
C Column J has not yet been accessed during this pass
              J1 = J
              J = IPERM(I)
              CV(I) = ID1
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
C Backtracking step.
   50       J1 = PR(J)
            IF (J1.EQ.-1) THEN
C No augmenting path exists for column J.
              NFC = NFC + 1
              FC(NFC) = J
              IF (NFC.GT.NUM2) THEN
C A matching of maximum size NUM1 is not possible
                LAST = JORD
                GO TO 101
              ENDIF
              GO TO 100
            ENDIF
            J = J1
   60     CONTINUE
C End of dummy loop; this point is never reached
   70   CONTINUE
C End of dummy loop; this point is never reached

C New assignment is made.
   80   IPERM(I) = J
        ARP(J) = II - IP(J) + 1
        NUM = NUM + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 95
          II = IP(J) + LENC(J) - OUT(J) - 2
          I = IRN(II)
          IPERM(I) = J
   90   CONTINUE
C End of dummy loop; this point is never reached

   95   IF (NUM.EQ.NUM1) THEN
C A matching of maximum size NUM1 is found
          LAST = JORD
          GO TO 101
        ENDIF
C
  100 CONTINUE

C All unassigned columns have been considered
      LAST = N

C Now, a transversal is computed or is not possible.
C Complete FC before returning.
  101 DO 110 JORD = LAST+1,N
        NFC = NFC + 1
        FC(NFC) = FC(JORD-NUM0)
  110 CONTINUE

C  199 RETURN
      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64WD(N,NE,IP,IRN,A,IPERM,NUM,
     &           JPERM,OUT,PR,Q,L,U,D)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...                                   ***
C     Iain Duff (I.Duff@rl.ac.uk) or                                ***
C     Jacko Koster (jacko.koster@uninett.no)                        ***
C
      INTEGER N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(N),
     &        JPERM(N),OUT(N),PR(N),Q(N),L(N)
      DOUBLE PRECISION A(NE),U(N),D(N)

C N, NE, IP, IRN are described in MC64A/AD.
C A is a DOUBLE PRECISION array of length NE.
C   A(K), K=1..NE, must be set to the value of the entry that
C   corresponds to IRN(K). It is not altered.
C   All values A(K) must be non-negative.
C IPERM is an INTEGER array of length N. On exit, it contains the
C   weighted matching: IPERM(I) = 0 or row I is matched to column
C   IPERM(I).
C NUM is an INTEGER variable. On exit, it contains the cardinality of
C   the matching stored in IPERM.
C IW is an INTEGER work array of length 5N.
C DW is a DOUBLE PRECISION array of length 2N.
C   On exit, U = D(1:N) contains the dual row variable and
C   V = D(N+1:2N) contains the dual column variable. If the matrix
C   is structurally nonsingular (NUM = N), the following holds:
C      U(I)+V(J) <= A(I,J)  if IPERM(I) |= J
C      U(I)+V(J)  = A(I,J)  if IPERM(I)  = J
C      U(I) = 0  if IPERM(I) = 0
C      V(J) = 0  if there is no I for which IPERM(I) = J

C Local variables
      INTEGER I,I0,II,J,JJ,JORD,Q0,QLEN,JDUM,ISP,JSP,
     &        K,K0,K1,K2,KK,KK1,KK2,UP,LOW,LPOS
      DOUBLE PRECISION CSP,DI,DMIN,DNEW,DQ0,VJ
C Local parameters
      DOUBLE PRECISION RINF,ZERO
      PARAMETER (ZERO=0.0D+0)
C External subroutines and/or functions
      EXTERNAL MC64DD,MC64ED,MC64FD


C Set RINF to largest positive real number
      RINF = HUGE(RINF)

C Initialization
      NUM = 0
      DO 10 K = 1,N
        U(K) = RINF
        D(K) = ZERO
        IPERM(K) = 0
        JPERM(K) = 0
        PR(K) = IP(K)
        L(K) = 0
   10 CONTINUE
C Initialize U(I)
      DO 30 J = 1,N
        DO 20 K = IP(J),IP(J+1)-1
          I = IRN(K)
          IF (A(K).GT.U(I)) GO TO 20
          U(I) = A(K)
          IPERM(I) = J
          L(I) = K
   20   CONTINUE
   30 CONTINUE
      DO 40 I = 1,N
        J = IPERM(I)
        IF (J.EQ.0) GO TO 40
C Row I is not empty
        IPERM(I) = 0
        IF (JPERM(J).NE.0) GO TO 40
C Don't choose cheap assignment from dense columns
        IF (IP(J+1)-IP(J) .GT. N/10 .AND. N.GT.50) GO TO 40
C Assignment of column J to row I
        NUM = NUM + 1
        IPERM(I) = J
        JPERM(J) = L(I)
   40 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
C Scan unassigned columns; improve assignment
      DO 95 J = 1,N
C JPERM(J) ne 0 iff column J is already assigned
        IF (JPERM(J).NE.0) GO TO 95
        K1 = IP(J)
        K2 = IP(J+1) - 1
C Continue only if column J is not empty
        IF (K1.GT.K2) GO TO 95
C       VJ = RINF
C Changes made to allow for NaNs
        I0 = IRN(K1)
        VJ = A(K1) - U(I0)
        K0 = K1
        DO 50 K = K1+1,K2
          I = IRN(K)
          DI = A(K) - U(I)
          IF (DI.GT.VJ) GO TO 50
          IF (DI.LT.VJ .OR. DI.EQ.RINF) GO TO 55
          IF (IPERM(I).NE.0 .OR. IPERM(I0).EQ.0) GO TO 50
   55     VJ = DI
          I0 = I
          K0 = K
   50   CONTINUE
        D(J) = VJ
        K = K0
        I = I0
        IF (IPERM(I).EQ.0) GO TO 90
        DO 60 K = K0,K2
          I = IRN(K)
          IF (A(K)-U(I).GT.VJ) GO TO 60
          JJ = IPERM(I)
C Scan remaining part of assigned column JJ
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 60
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).GT.0) GO TO 70
            IF (A(KK)-U(II).LE.D(JJ)) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   60   CONTINUE
        GO TO 95
   80   JPERM(JJ) = KK
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = K
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000

C Prepare for main loop
      DO 99 I = 1,N
        D(I) = RINF
        L(I) = 0
   99 CONTINUE

C Main loop ... each pass round this loop is similar to Dijkstra's
C algorithm for solving the single source shortest path problem

      DO 100 JORD = 1,N

        IF (JPERM(JORD).NE.0) GO TO 100
C JORD is next unmatched column
C DMIN is the length of shortest path in the tree
        DMIN = RINF
        QLEN = 0
        LOW = N + 1
        UP = N + 1
C CSP is the cost of the shortest augmenting path to unassigned row
C IRN(ISP). The corresponding column index is JSP.
        CSP = RINF
C Build shortest path tree starting from unassigned column (root) JORD
        J = JORD
        PR(J) = -1

C Scan column J
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = A(K) - U(I)
          IF (DNEW.GE.CSP) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = K
            JSP = J
          ELSE
            IF (DNEW.LT.DMIN) DMIN = DNEW
            D(I) = DNEW
            QLEN = QLEN + 1
            Q(QLEN) = K
          ENDIF
  115   CONTINUE
C Initialize heap Q and Q2 with rows held in Q(1:QLEN)
        Q0 = QLEN
        QLEN = 0
        DO 120 KK = 1,Q0
          K = Q(KK)
          I = IRN(K)
          IF (CSP.LE.D(I)) THEN
            D(I) = RINF
            GO TO 120
          ENDIF
          IF (D(I).LE.DMIN) THEN
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
          ELSE
            QLEN = QLEN + 1
            L(I) = QLEN
            CALL MC64DD(I,N,Q,D,L,2)
          ENDIF
C Update tree
          JJ = IPERM(I)
          OUT(JJ) = K
          PR(JJ) = J
  120   CONTINUE

        DO 150 JDUM = 1,NUM

C If Q2 is empty, extract rows from Q
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (D(I).GE.CSP) GO TO 160
            DMIN = D(I)
  152       CALL MC64ED(QLEN,N,Q,D,L,2)
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
            IF (QLEN.EQ.0) GO TO 153
            I = Q(1)
            IF (D(I).GT.DMIN) GO TO 153
            GO TO 152
          ENDIF
C Q0 is row whose distance D(Q0) to the root is smallest
  153     Q0 = Q(UP-1)
          DQ0 = D(Q0)
C Exit loop if path to Q0 is longer than the shortest augmenting path
          IF (DQ0.GE.CSP) GO TO 160
          UP = UP - 1

C Scan column that matches with row Q0
          J = IPERM(Q0)
          VJ = DQ0 - A(JPERM(J)) + U(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
            IF (L(I).GE.UP) GO TO 155
C DNEW is new cost
            DNEW = VJ + A(K)-U(I)
C Do not update D(I) if DNEW ge cost of shortest path
            IF (DNEW.GE.CSP) GO TO 155
            IF (IPERM(I).EQ.0) THEN
C Row I is unmatched; update shortest path info
              CSP = DNEW
              ISP = K
              JSP = J
            ELSE
C Row I is matched; do not update D(I) if DNEW is larger
              DI = D(I)
              IF (DI.LE.DNEW) GO TO 155
              IF (L(I).GE.LOW) GO TO 155
              D(I) = DNEW
              IF (DNEW.LE.DMIN) THEN
                LPOS = L(I)
                IF (LPOS.NE.0)
     *            CALL MC64FD(LPOS,QLEN,N,Q,D,L,2)
                LOW = LOW - 1
                Q(LOW) = I
                L(I) = LOW
              ELSE
                IF (L(I).EQ.0) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64DD(I,N,Q,D,L,2)
              ENDIF
C Update tree
              JJ = IPERM(I)
              OUT(JJ) = K
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE

C If CSP = RINF, no augmenting path is found
  160   IF (CSP.EQ.RINF) GO TO 190
C Find augmenting path by tracing backward in PR; update IPERM,JPERM
        NUM = NUM + 1
        I = IRN(ISP)
        IPERM(I) = JSP
        JPERM(JSP) = ISP
        J = JSP
        DO 170 JDUM = 1,NUM
          JJ = PR(J)
          IF (JJ.EQ.-1) GO TO 180
          K = OUT(J)
          I = IRN(K)
          IPERM(I) = JJ
          JPERM(JJ) = K
          J = JJ
  170   CONTINUE
C End of dummy loop; this point is never reached

C Update U for rows in Q(UP:N)
  180   DO 185 KK = UP,N
          I = Q(KK)
          U(I) = U(I) + D(I) - CSP
  185   CONTINUE
  190   DO 191 KK = LOW,N
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  191   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  193   CONTINUE

  100 CONTINUE
C End of main loop


C Set dual column variable in D(1:N)
 1000 DO 200 J = 1,N
        K = JPERM(J)
        IF (K.NE.0) THEN
          D(J) = A(K) - U(IRN(K))
        ELSE
          D(J) = ZERO
        ENDIF
        IF (IPERM(J).EQ.0) U(J) = ZERO
  200 CONTINUE

      IF (NUM.EQ.N) GO TO 1100

C The matrix is structurally singular, complete IPERM.
C JPERM, OUT are work arrays
      DO 300 J = 1,N
        JPERM(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          OUT(K) = I
        ELSE
          J = IPERM(I)
          JPERM(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 J = 1,N
        IF (JPERM(J).NE.0) GO TO 320
        K = K + 1
        JDUM = OUT(K)
        IPERM(JDUM) = - J
  320 CONTINUE
 1100 RETURN
      END


* COPYRIGHT (c) 1988 AEA Technology and
* Council for the Central Laboratory of the Research Councils
C Original date 14 June 2001
C  June 2001: threadsafe version of MC41
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC71AD(N,KASE,X,EST,W,IW,KEEP)
C
C      MC71A/AD ESTIMATES THE 1-NORM OF A SQUARE MATRIX A.
C      REVERSE COMMUNICATION IS USED FOR EVALUATING
C      MATRIX-VECTOR PRODUCTS.
C
C
C         N       INTEGER
C                 THE ORDER OF THE MATRIX.  N .GE. 1.
C
C         KASE    INTEGER
C                 SET INITIALLY TO ZERO . IF N .LE. 0 SET TO -1
C                 ON INTERMEDIATE RETURN
C                 = 1 OR 2.
C                ON FINAL RETURN
C                 =  0  ,IF SUCCESS
C                 = -1  ,IF N .LE.0
C
C         X       DOUBLE PRECISION ARRAY OF DIMENSION (N)
C                 IF 1-NORM IS REQUIRED
C                 MUST BE OVERWRITTEN BY
C
C                      A*X,             IF KASE=1,
C                      TRANSPOSE(A)*X,  IF KASE=2,
C
C                 AND MC71 MUST BE RE-CALLED, WITH ALL THE OTHER
C                 PARAMETERS UNCHANGED.
C                 IF INFINITY-NORM IS REQUIRED
C                 MUST BE OVERWRITTEN BY
C
C                      TRANSPOSE(A)*X,  IF KASE=1,
C                      A*X,             IF KASE=2,
C
C                 AND MC71 MUST BE RE-CALLED, WITH ALL THE OTHER
C                 PARAMETERS UNCHANGED.
C
C         EST     DOUBLE PRECISION
C                 CONTAINS AN ESTIMATE (A LOWER BOUND) FOR NORM(A).
C
C         W       DOUBLE PRECISION ARRAY OF DIMENSION (N)
C                 = A*V,   WHERE  EST = NORM(W)/NORM(V)
C                          (V  IS NOT RETURNED).
C         IW      INTEGER(N) USED AS WORKSPACE.
C
C         KEEP    INTEGER ARRAY LENGTH 5 USED TO PRESERVE PRIVATE
C                 DATA, JUMP, ITER, J AND JLAST BETWEEN CALLS,
C                 KEEP(5) IS SPARE.
C
C      REFERENCE
C      N.J. HIGHAM (1987) FORTRAN CODES FOR ESTIMATING
C      THE ONE-NORM OF A
C      REAL OR COMPLEX MATRIX, WITH APPLICATIONS
C      TO CONDITION  ESTIMATION, NUMERICAL ANALYSIS REPORT NO. 135,
C      UNIVERSITY OF MANCHESTER, MANCHESTER M13 9PL, ENGLAND.
C
C      SUBROUTINES AND FUNCTIONS
C
C
C      INTERNAL VARIABLES
C
C     .. Parameters ..
      INTEGER ITMAX
      PARAMETER (ITMAX=5)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EST
      INTEGER KASE,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION W(*),X(*)
      INTEGER IW(*),KEEP(5)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALTSGN,TEMP
      INTEGER I,ITER,J,JLAST,JUMP
C     ..
C     .. External Functions ..
      INTEGER IDAMAX
      EXTERNAL IDAMAX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN,NINT,DBLE
C     ..
C     .. Executable Statements ..
C
      IF (N.LE.0) THEN
        KASE = -1
        RETURN

      END IF

      IF (KASE.EQ.0) THEN
        DO 10 I = 1,N
          X(I) = ONE/DBLE(N)
   10   CONTINUE
        KASE = 1
        JUMP = 1
        KEEP(1) = JUMP
        KEEP(2) = 0
        KEEP(3) = 0
        KEEP(4) = 0
        RETURN

      END IF
C
      JUMP  = KEEP(1)
      ITER  = KEEP(2)
      J     = KEEP(3)
      JLAST = KEEP(4)
C
      GO TO (100,200,300,400,500) JUMP
C
C      ................ ENTRY   (JUMP = 1)
C
  100 CONTINUE
      IF (N.EQ.1) THEN
        W(1) = X(1)
        EST = ABS(W(1))
C         ... QUIT
        GO TO 510

      END IF
C
      DO 110 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  110 CONTINUE
      KASE = 2
      JUMP = 2
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 2)
C
  200 CONTINUE
      J = IDAMAX(N,X,1)
      ITER = 2
C
C      MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
C
  220 CONTINUE
      DO 230 I = 1,N
        X(I) = ZERO
  230 CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 3)
C
  300 CONTINUE
C
C      COPY X INTO W
C
      DO 310 I = 1,N
        W(I) = X(I)
  310 CONTINUE
      DO 320 I = 1,N
        IF (NINT(SIGN(ONE,X(I))).NE.IW(I)) GO TO 330
  320 CONTINUE
C
C      REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 410
C
  330 CONTINUE
      DO 340 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  340 CONTINUE
      KASE = 2
      JUMP = 4
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 4)
C
  400 CONTINUE
      JLAST = J
      J = IDAMAX(N,X,1)
      IF ((ABS(X(JLAST)).NE.ABS(X(J))) .AND. (ITER.LT.ITMAX)) THEN
        ITER = ITER + 1
        GO TO 220

      END IF
C
C      ITERATION COMPLETE.  FINAL STAGE.
C
  410 CONTINUE
      EST = ZERO
      DO 420 I = 1,N
        EST = EST + ABS(W(I))
  420 CONTINUE
C
      ALTSGN = ONE
      DO 430 I = 1,N
        X(I) = ALTSGN* (ONE+DBLE(I-1)/DBLE(N-1))
        ALTSGN = -ALTSGN
  430 CONTINUE
      KASE = 1
      JUMP = 5
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 5)
C
  500 CONTINUE
      TEMP = ZERO
      DO 520 I = 1,N
        TEMP = TEMP + ABS(X(I))
  520 CONTINUE
      TEMP = 2.0*TEMP/DBLE(3*N)
      IF (TEMP.GT.EST) THEN
C
C      COPY X INTO W
C
        DO 530 I = 1,N
          W(I) = X(I)
  530   CONTINUE
        EST = TEMP
      END IF
C
  510 KASE = 0
C
 1010 CONTINUE
      KEEP(1) = JUMP
      KEEP(2) = ITER
      KEEP(3) = J
      KEEP(4) = JLAST
      RETURN
C
      END
