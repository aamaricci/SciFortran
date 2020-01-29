c$$$     To find the minimum of a function F(x) of n variables x. It is assumed that the function is differentiable, although it is not necessary to supply a formula for the derivatives. The method used is a quasi-Newton method in which derivatives are estimated by differences and is described in R. Fletcher, FORTRAN subroutines for minimization by quasi-Newton methods, AERE Report R7125 (1972). The subroutine complements VA04 but some comparisons (R. Fletcher, loc. cit) indicate that VA04 is less efficient than VA10 and more affected by round off error. VA04 also uses twice as much storage as VA10. It is therefore suggested that VA10 be used in the first instance on any problem. VA10 should not be used when explicit expressions are available for derivatives (use VA09) nor when the function is a sum of squares (use VA05, VA02 or one of the NS routines as appropriate).

c$$$
c$$$ CALL VA10AD(FUNCT,N,X,F,G,H,W,DFN,XM,HH,EPS, MODE,MAXFN,IPRINT,IEXIT)
c$$$ FUNCT is the name of the subroutine provided by the user which is considered in section 3. It must be declared in an EXTERNAL statement.
c$$$ N is an INTEGER to be set to the number of variables (N ≥ 2).
c$$$ X is a REAL (DOUBLE PRECISION in the D version) array in which the solution is stored. An initial approximation
c$$$must be set in X on entry to VA10 and the best estimate obtained will be returned on exit.
c$$$ F is a REAL (DOUBLE PRECISION in the D version) number in which the best value of F(x) corresponding to X
c$$$above will be returned.
c$$$ G is a REAL (DOUBLE PRECISION in the D version) array of N elements which is used to store an estimate of the
c$$$gradient vector ∂F(x)/∂x. This array need not be set on entry.
c$$$ H is a REAL (DOUBLE PRECISION in the D version) array of dimension N(N+1)/2 elements in which an estimate of the Hessian matrix ∂ 2 F/(∂x i ∂x j ) is stored. The matrix is represented in the product form LDLT where L is a lower triangular matrix with unit diagonals and D is a diagonal matrix. The lower triangle of L is stored by columns in H excepting that the unit diagonal elements are replaced by the corresponding elements of D. The setting of H on entry is controlled by the parameter MODE (q.v.).
c$$$ A is a REAL (DOUBLE PRECISION in the D version) array of 3N elements used as working space.
c$$$ DFN isaREAL(DOUBLEPRECISIONintheDversion)numberwhichmustbesetsoastogiveVA10anestimateof the likely reduction to be obtained in F(x). DFN is used only on the first iteration so an order of magnitude estimate will suffice. The information can be provided in different ways depending upon the sign of DFN which should be set in one of the following ways:
c$$$
c$$$ if DFN > 0 the setting of DFN itself will be taken as the likely reduction to be obtained in F(x).
c$$$ if DFN=0 it will be assumed that an estimate of the minimum value of F(x) has been set in argument F, and the
c$$$likely reduction in F(x) will be computed according to the initial function value.
c$$$ if DFN < 0 a multiple |DFN| of the modulus of the initial function value will be taken as an estimate of the likely
c$$$reduction.
c$$$ XM is a REAL (DOUBLE PRECISION in the D version) array of N elements to be set on entry so that XM(I) > 0 contains an indication of the magnitude of X(I). This quantity need not be set precisely as it is merely used in scaling the problem.
c$$$ HH is a REAL (DOUBLE PRECISION in the D version) number to be set so that HH XM(I) contains a step length to be usedincalculatingG(I)bydifferences.SetHHequalto2−t/2 wheretisthenumberofsignificantbinarydigits in the calculation of F. If F contains only small errors the setting HH=1.E-3 is appropriate for VA10A and HH=1.0D-6 for VA10AD.
c$$$ EPS is a REAL (DOUBLE PRECISION in the D version) number to be set on entry so that the accuracy required in X(I) is EPS XM(I) for all I, (EPS > 0).
c$$$ MODE is an INTEGER which controls the setting of the initial estimate of the Hessian matrix in the parameter H. The following settings of MODE are permitted:
c$$$ if MODE=1 an estimate corresponding to a unit matrix is set in H by VA10A.
c$$$ if MODE=2 VA10A/AD assumes that the Hessian matrix itself has been set in H by columns of its lower triangle,
c$$$ T
c$$$and the conversion to LDL form is carried out by VA10. The Hessian matrix must be positive definite.
c$$$i f MODE=3 VA10A/AD assumes that the Hessian matrix has been set in H in product form. This is convenient when using the H matrix from one problem as an initial estimate for another, in which case the contents of H are passed on unchanged.
c$$$ MAXFN is an INTEGER set to the maximum number of calls of FUNCT permitted. Up to 2N more calls may be taken if the limit is exceeded whilst evaluating a gradient vector by differences.
c$$$ IPRINT An is an INTEGER controlling printing. Printing occurs every |IPRINT| iterations and also on exit, in the form Iteration No. No of calls of FUNCT,IEXIT (on exit only) Function only X(1),X(2),...,X(N) 8 to a line (5 in VA10AD) G(1),G(2),...,G(N) 8 to a line (5 in VA10AD)
c$$$The values of X and G can be suppressed on intermediate iterations by setting IPRINT < 0. All intermediate printing can be suppressed by setting IPRINT=MAXFN+1. All printing can be suppressed by setting IPRINT=0.
c$$$ IEXIT is an INTEGER giving the reason for exit from VA10A/AD. This will be set by VA10A/AD as follows if IEXIT=0 (MODE=2 only) the estimate of the Hessian matrix is not positive definite.
c$$$ if IEXIT=1 a normal exit has been made in which |DX(I)| < EPS(I) for all I=1,2,...,N, where DX(I) is the change in X on an iteration.
c$$$ T
c$$$ if IEXIT=2 G DX > 0. This is an error exit, either due to rounding errors because EPS is set too small for the computer word length, or because the truncation error in the finite difference formula for G is dominant.
c$$$ if IEXIT=3 FUNCT has been called MAXFN times.

      
      
*######DATE 9 February 1994 COPYRIGHT AEA Technology
C       Toolpack tool decs employed.
C       Arg dimensions made *.
C
C HSL VERSION CONTROL RECORD
      subroutine minimize_sascha_(funct, n, x, f, g, h, w, dfn, xm,
     $  hh, eps, mode, maxfn, iprint, iexit,itn)
c     SUBROUTINE VA10AD(FUNCT,N,X,F,G,H,W,DFN,XM,HH,EPS,MODE,MAXFN,
c    +                  IPRINT,IEXIT)
      DOUBLE PRECISION DFN,EPS,F,HH
      INTEGER IEXIT,IPRINT,MAXFN,MODE,N
      DOUBLE PRECISION G(*),H(*),W(*),X(*),XM(*)
      EXTERNAL FUNCT
      DOUBLE PRECISION AEPS,ALPHA,DF,DGS,EPSMCH,F1,F2,FF,GS0,GYS,SIG,
     +                 TOT,Z,ZZ
      INTEGER I,IDIFF,IFN,IG,IGG,IJ,INT,IR,IS,ITN,J,LINK,NN
      DOUBLE PRECISION FD05AD
      EXTERNAL FD05AD
      EXTERNAL MC11AD,MC11BD,MC11ED
      INTRINSIC DABS,MOD
      EPSMCH = FD05AD(1)*10.0D0
      IF (IPRINT.NE.0) WRITE (6,FMT=1000)
 1000 FORMAT ('1ENTRY TO VA10AD',/)
      NN = N* (N+1)/2
      IG = N
      IGG = N + N
      IS = IGG
      IDIFF = 1
      IEXIT = 0
      IR = N
      IF (MODE.EQ.3) GO TO 15
      IF (MODE.EQ.2) GO TO 10
      IJ = NN + 1
      DO 5 I = 1,N
        DO 6 J = 1,I
          IJ = IJ - 1
    6   H(IJ) = 0.D0
    5 H(IJ) = 1.D0
      GO TO 15
   10 CONTINUE
      CALL MC11BD(H,N,IR)
      IF (IR.LT.N) RETURN
   15 CONTINUE
      Z = F
      ITN = 0
      CALL FUNCT(N,X,F)
      IFN = 1
      DF = DFN
      IF (DFN.EQ.0D0) DF = F - Z
      IF (DFN.LT.0D0) DF = DABS(DF*F)
      IF (DF.LE.0D0) DF = 1.D0
   17 CONTINUE
      LINK = 1
      IF (IDIFF-1) 100,100,110
   18 CONTINUE
      IF (IFN.GE.MAXFN) GO TO 90
   20 CONTINUE
      IF (IPRINT.EQ.0) GO TO 21
      IF (MOD(ITN,IPRINT).NE.0) GO TO 21
      WRITE (6,FMT=1001) ITN,IFN
 1001 FORMAT (24I5)
      WRITE (6,FMT=1002) F
 1002 FORMAT (5D24.16)
      IF (IPRINT.LT.0) GO TO 21
      WRITE (6,FMT=1002) (X(I),I=1,N)
      WRITE (6,FMT=1002) (W(IG+I),I=1,N)
   21 CONTINUE
      ITN = ITN + 1
      DO 22 I = 1,N
   22 W(I) = -W(IG+I)
      CALL MC11ED(H,N,W,G,IR)
      Z = 0.D0
      GS0 = 0.D0
      DO 29 I = 1,N
        W(IS+I) = W(I)
        IF (Z*XM(I).GE.DABS(W(I))) GO TO 29
        Z = DABS(W(I))/XM(I)
   29 GS0 = GS0 + W(IG+I)*W(I)
      AEPS = EPS/Z
      IEXIT = 2
      IF (GS0.GE.0D0) GO TO 92
      ALPHA = -2D0*DF/GS0
      IF (ALPHA.GT.1D0) ALPHA = 1D0
      FF = F
      TOT = 0.D0
      INT = 0
      IEXIT = 1
   30 CONTINUE
      IF (IFN.GE.MAXFN) GO TO 90
      DO 31 I = 1,N
   31 W(I) = X(I) + ALPHA*W(IS+I)
      CALL FUNCT(N,W,F1)
      IFN = IFN + 1
      IF (F1.GE.F) GO TO 40
      F2 = F
      TOT = TOT + ALPHA
   32 CONTINUE
      DO 33 I = 1,N
   33 X(I) = W(I)
      F = F1
      IF (INT-1) 35,49,50
   35 CONTINUE
      IF (IFN.GE.MAXFN) GO TO 90
      DO 34 I = 1,N
   34 W(I) = X(I) + ALPHA*W(IS+I)
      CALL FUNCT(N,W,F1)
      IFN = IFN + 1
      IF (F1.GE.F) GO TO 50
      IF (F1+F2.GE.F+F .AND. 7D0*F1+5D0*F2.GT.12D0*F) INT = 2
      TOT = TOT + ALPHA
      ALPHA = 2D0*ALPHA
      GO TO 32
   40 CONTINUE
      IF (ALPHA.LT.AEPS) GO TO 92
      IF (IFN.GE.MAXFN) GO TO 90
      ALPHA = .5D0*ALPHA
      DO 41 I = 1,N
   41 W(I) = X(I) + ALPHA*W(IS+I)
      CALL FUNCT(N,W,F2)
      IFN = IFN + 1
      IF (F2.GE.F) GO TO 45
      TOT = TOT + ALPHA
      F = F2
      DO 42 I = 1,N
   42 X(I) = W(I)
      GO TO 49
   45 CONTINUE
      Z = .1D0
      IF (F1+F.GT.F2+F2) Z = 1D0 + .5D0* (F-F1)/ (F+F1-F2-F2)
      IF (Z.LT..1D0) Z = .1D0
      ALPHA = Z*ALPHA
      INT = 1
      GO TO 30
   49 CONTINUE
      IF (TOT.LT.AEPS) GO TO 92
   50 CONTINUE
      ALPHA = TOT
      DO 56 I = 1,N
   56 W(I) = W(IG+I)
      LINK = 2
      IF (IDIFF-1) 100,100,110
   54 CONTINUE
      IF (IFN.GE.MAXFN) GO TO 90
      GYS = 0.D0
      DO 55 I = 1,N
        GYS = GYS + W(IG+I)*W(IS+I)
   55 W(IGG+I) = W(I)
      DF = FF - F
      DGS = GYS - GS0
      IF (DGS.LE.0D0) GO TO 20
      IF (DGS+ALPHA*GS0.GT.0D0) GO TO 70
      SIG = 1D0/GS0
      IR = -IR
      CALL MC11AD(H,N,W,SIG,G,IR,1,0D0)
      DO 60 I = 1,N
   60 G(I) = W(IG+I) - W(IGG+I)
      SIG = 1D0/ (ALPHA*DGS)
      IR = -IR
      CALL MC11AD(H,N,G,SIG,W,IR,0,0D0)
      GO TO 20
   70 CONTINUE
      ZZ = ALPHA/ (DGS-ALPHA*GS0)
      SIG = -ZZ
      CALL MC11AD(H,N,W,SIG,G,IR,1,EPSMCH)
      Z = DGS*ZZ - 1.D0
      DO 71 I = 1,N
   71 G(I) = W(IG+I) + Z*W(IGG+I)
      SIG = 1.D0/ (ZZ*DGS**2)
      CALL MC11AD(H,N,G,SIG,W,IR,0,0D0)
      GO TO 20
   90 CONTINUE
      IEXIT = 3
      GO TO 94
   92 CONTINUE
      IF (IDIFF.EQ.2) GO TO 94
      IDIFF = 2
      GO TO 17
   94 CONTINUE
      DO 95 I = 1,N
   95 G(I) = W(IG+I)
      IF (IPRINT.EQ.0) RETURN
      WRITE (6,FMT=1001) ITN,IFN,IEXIT
      WRITE (6,FMT=1002) F
      WRITE (6,FMT=1002) (X(I),I=1,N)
      WRITE (6,FMT=1002) (G(I),I=1,N)
      RETURN
  100 CONTINUE
      DO 101 I = 1,N
        Z = HH*XM(I)
        ZZ = X(I)
        X(I) = ZZ + Z
        CALL FUNCT(N,X,F1)
        W(IG+I) = (F1-F)/Z
  101 X(I) = ZZ
      IFN = IFN + N
  102 GO TO (18,54) LINK
  110 CONTINUE
      DO 111 I = 1,N
        Z = HH*XM(I)
        ZZ = X(I)
        X(I) = ZZ + Z
        CALL FUNCT(N,X,F1)
        X(I) = ZZ - Z
        CALL FUNCT(N,X,F2)
        W(IG+I) = (F1-F2)/ (2D0*Z)
  111 X(I) = ZZ
      IFN = IFN + N + N
      GO TO 102
      END
*######DATE 1 Feb 1993 COPYRIGHT AEA Technology
C       Toolpack tool decs employed.
C       Arg dimensions set to *.
C
      SUBROUTINE MC11AD(A,N,Z,SIG,W,IR,MK,EPS)
      DOUBLE PRECISION EPS,SIG
      INTEGER IR,MK,N
      DOUBLE PRECISION A(*),W(*),Z(*)
      DOUBLE PRECISION AL,B,GM,R,TI,TIM,V,Y
      INTEGER I,IJ,IP,J,MM,NP
      IF (N.GT.1) GO TO 1
      A(1) = A(1) + SIG*Z(1)**2
      IR = 1
      IF (A(1).GT.0.0D0) RETURN
      A(1) = 0.D0
      IR = 0
      RETURN
    1 CONTINUE
      NP = N + 1
      IF (SIG.GT.0.0D0) GO TO 40
      IF (SIG.EQ.0.0D0 .OR. IR.EQ.0) RETURN
      TI = 1.0D0/SIG
      IJ = 1
      IF (MK.EQ.0) GO TO 10
      DO 7 I = 1,N
        IF (A(IJ).NE.0.0D0) TI = TI + W(I)**2/A(IJ)
    7 IJ = IJ + NP - I
      GO TO 20
   10 CONTINUE
      DO 11 I = 1,N
   11 W(I) = Z(I)
      DO 15 I = 1,N
        IP = I + 1
        V = W(I)
        IF (A(IJ).GT.0.0D0) GO TO 12
        W(I) = 0.D0
        IJ = IJ + NP - I
        GO TO 15
   12   CONTINUE
        TI = TI + V**2/A(IJ)
        IF (I.EQ.N) GO TO 14
        DO 13 J = IP,N
          IJ = IJ + 1
   13   W(J) = W(J) - V*A(IJ)
   14   IJ = IJ + 1
   15 CONTINUE
   20 CONTINUE
      IF (IR.LE.0) GO TO 21
      IF (TI.GT.0.0D0) GO TO 22
      IF (MK-1) 40,40,23
   21 TI = 0.D0
      IR = -IR - 1
      GO TO 23
   22 TI = EPS/SIG
      IF (EPS.EQ.0.0D0) IR = IR - 1
   23 CONTINUE
      MM = 1
      TIM = TI
      DO 30 I = 1,N
        J = NP - I
        IJ = IJ - I
        IF (A(IJ).NE.0.0D0) TIM = TI - W(J)**2/A(IJ)
        W(J) = TI
   30 TI = TIM
      GO TO 41
   40 CONTINUE
      MM = 0
      TIM = 1.0D0/SIG
   41 CONTINUE
      IJ = 1
      DO 66 I = 1,N
        IP = I + 1
        V = Z(I)
        IF (A(IJ).GT.0.0D0) GO TO 53
        IF (IR.GT.0 .OR. SIG.LT.0.0D0 .OR. V.EQ.0.0D0) GO TO 52
        IR = 1 - IR
        A(IJ) = V**2/TIM
        IF (I.EQ.N) RETURN
        DO 51 J = IP,N
          IJ = IJ + 1
   51   A(IJ) = Z(J)/V
        RETURN
   52   CONTINUE
        TI = TIM
        IJ = IJ + NP - I
        GO TO 66
   53   CONTINUE
        AL = V/A(IJ)
        IF (MM) 54,54,55
   54   TI = TIM + V*AL
        GO TO 56
   55   TI = W(I)
   56   CONTINUE
        R = TI/TIM
        A(IJ) = A(IJ)*R
        IF (R.EQ.0.0D0) GO TO 70
        IF (I.EQ.N) GO TO 70
        B = AL/TI
        IF (R.GT.4.0D0) GO TO 62
        DO 61 J = IP,N
          IJ = IJ + 1
          Z(J) = Z(J) - V*A(IJ)
   61   A(IJ) = A(IJ) + B*Z(J)
        GO TO 64
   62   GM = TIM/TI
        DO 63 J = IP,N
          IJ = IJ + 1
          Y = A(IJ)
          A(IJ) = B*Z(J) + Y*GM
   63   Z(J) = Z(J) - V*Y
   64   CONTINUE
        TIM = TI
        IJ = IJ + 1
   66 CONTINUE
   70 CONTINUE
      IF (IR.LT.0) IR = -IR
      RETURN
      END
      SUBROUTINE MC11BD(A,N,IR)
      INTEGER IR,N
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION AA,V
      INTEGER I,II,IJ,IK,IP,JK,NI,NP
      IR = N
      IF (N.GT.1) GO TO 100
      IF (A(1).GT.0.0D0) RETURN
      A(1) = 0.D0
      IR = 0
      RETURN
  100 CONTINUE
      NP = N + 1
      II = 1
      DO 104 I = 2,N
        AA = A(II)
        NI = II + NP - I
        IF (AA.GT.0.0D0) GO TO 101
        A(II) = 0.D0
        IR = IR - 1
        II = NI + 1
        GO TO 104
  101   CONTINUE
        IP = II + 1
        II = NI + 1
        JK = II
        DO 103 IJ = IP,NI
          V = A(IJ)/AA
          DO 102 IK = IJ,NI
            A(JK) = A(JK) - A(IK)*V
  102     JK = JK + 1
  103   A(IJ) = V
  104 CONTINUE
      IF (A(II).GT.0.0D0) RETURN
      A(II) = 0.D0
      IR = IR - 1
      RETURN
      END
      SUBROUTINE MC11CD(A,N)
      INTEGER N
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION AA,V
      INTEGER II,IJ,IK,IP,JK,NI,NIP,NP
      IF (N.EQ.1) RETURN
      NP = N + 1
      II = N*NP/2
      DO 202 NIP = 2,N
        JK = II
        NI = II - 1
        II = II - NIP
        AA = A(II)
        IP = II + 1
        IF (AA.GT.0.0D0) GO TO 203
        DO 204 IJ = IP,NI
  204   A(IJ) = 0.D0
        GO TO 202
  203   CONTINUE
        DO 201 IJ = IP,NI
          V = A(IJ)*AA
          DO 200 IK = IJ,NI
            A(JK) = A(JK) + A(IK)*V
  200     JK = JK + 1
  201   A(IJ) = V
  202 CONTINUE
      RETURN
      END
      SUBROUTINE MC11DD(A,N,Z,W)
      INTEGER N
      DOUBLE PRECISION A(*),W(*),Z(*)
      DOUBLE PRECISION Y
      INTEGER I,II,IJ,IP,J,K,N1,NP
      IF (N.GT.1) GO TO 300
      Z(1) = Z(1)*A(1)
      W(1) = Z(1)
      RETURN
  300 CONTINUE
      NP = N + 1
      II = 1
      N1 = N - 1
      DO 303 I = 1,N1
        Y = Z(I)
        IF (A(II).EQ.0.0D0) GO TO 302
        IJ = II
        IP = I + 1
        DO 301 J = IP,N
          IJ = IJ + 1
  301   Y = Y + Z(J)*A(IJ)
  302   Z(I) = Y*A(II)
        W(I) = Z(I)
  303 II = II + NP - I
      Z(N) = Z(N)*A(II)
      W(N) = Z(N)
      DO 311 K = 1,N1
        I = N - K
        II = II - NP + I
        IF (Z(I).EQ.0.0D0) GO TO 311
        IP = I + 1
        IJ = II
        Y = Z(I)
        DO 310 J = IP,N
          IJ = IJ + 1
  310   Z(J) = Z(J) + A(IJ)*Z(I)
  311 CONTINUE
      RETURN
      END
      SUBROUTINE MC11ED(A,N,Z,W,IR)
      INTEGER IR,N
      DOUBLE PRECISION A(*),W(*),Z(*)
      DOUBLE PRECISION V
      INTEGER I,I1,II,IJ,IP,J,NIP,NP
      IF (IR.LT.N) RETURN
      W(1) = Z(1)
      IF (N.GT.1) GO TO 400
      Z(1) = Z(1)/A(1)
      RETURN
  400 CONTINUE
      DO 402 I = 2,N
        IJ = I
        I1 = I - 1
        V = Z(I)
        DO 401 J = 1,I1
          V = V - A(IJ)*Z(J)
  401   IJ = IJ + N - J
        W(I) = V
  402 Z(I) = V
      Z(N) = Z(N)/A(IJ)
      NP = N + 1
      DO 411 NIP = 2,N
        I = NP - NIP
        II = IJ - NIP
        V = Z(I)/A(II)
        IP = I + 1
        IJ = II
        DO 410 J = IP,N
          II = II + 1
  410   V = V - A(II)*Z(J)
  411 Z(I) = V
      RETURN
      END
      SUBROUTINE MC11FD(A,N,IR)
      INTEGER IR,N
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION AA,V
      INTEGER I,I1,II,IJ,IK,IP,J,JK,K,N1,NI,NIP,NP
      IF (IR.LT.N) RETURN
      A(1) = 1.0D0/A(1)
      IF (N.EQ.1) RETURN
      NP = N + 1
      N1 = N - 1
      II = 2
      DO 511 I = 2,N
        A(II) = -A(II)
        IJ = II + 1
        IF (I.EQ.N) GO TO 502
        DO 501 J = I,N1
          IK = II
          JK = IJ
          V = A(IJ)
          DO 500 K = I,J
            JK = JK + NP - K
            V = V + A(IK)*A(JK)
  500     IK = IK + 1
          A(IJ) = -V
  501   IJ = IJ + 1
  502   CONTINUE
        A(IJ) = 1.0D0/A(IJ)
        II = IJ + 1
        AA = A(IJ)
        IJ = I
        IP = I + 1
        NI = N - I
        DO 511 J = 2,I
          V = A(IJ)*AA
          IK = IJ
          K = IJ - IP + J
          I1 = IJ - 1
          NIP = NI + IJ
          DO 510 JK = K,I1
            A(JK) = A(JK) + V*A(IK)
  510     IK = IK + NIP - JK
          A(IJ) = V
  511 IJ = IJ + NP - J
      RETURN
      END
C#####################################################
      DOUBLE PRECISION FUNCTION FD05AD( INUM )
      INTEGER INUM
      DOUBLE PRECISION DC( 5 )
C
C  REAL CONSTANTS (DOUBLE PRECISION ARITHMETIC).
C
C  OBTAINED FROM H.S.L. SUBROUTINE ZE02AM.
C  NICK GOULD AND SID MARLOW, HARWELL, JULY 1988.
C
C  DC(1) THE 'SMALLEST' POSITIVE NUMBER: 1 + DC(1) > 1.
C  DC(2) THE 'SMALLEST' POSITIVE NUMBER: 1 - DC(2) < 1.
C  DC(3) THE SMALLEST NONZERO +VE REAL NUMBER.
C  DC(4) THE SMALLEST FULL PRECISION +VE REAL NUMBER.
C  DC(5) THE LARGEST FINITE +VE REAL NUMBER.
C
      SAVE DC
      DATA DC( 1 ) /    2.220446049250314D-016 /
      DATA DC( 2 ) /    1.110223024625158D-016 /
      DATA DC( 3 ) /    2.225073858507202D-308 /
      DATA DC( 4 ) /    2.225073858507202D-308 /
      DATA DC( 5 ) /    1.797693134862315D+308 /
      IF ( INUM .LE. 0 .OR. INUM .GE. 6 ) THEN
         WRITE(6, 2000) INUM
         STOP
      ELSE
         FD05AD = DC( INUM )
      ENDIF
      RETURN
 2000 FORMAT( ' INUM =', I3, ' OUT OF RANGE IN FD05AD.',
     *        ' EXECUTION TERMINATED.' )
      END
