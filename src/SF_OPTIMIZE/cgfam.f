
C
C     --------------------------------------------------------------------
C     Conjugate Gradient methods for solving unconstrained nonlinear
C     optimization problems, as described in the paper:
C
C     Gilbert, J.C. and Nocedal, J. (1992). "Global Convergence Properties 
C     of Conjugate Gradient Methods", SIAM Journal on Optimization, Vol. 2,
C     pp. 21-42. 
C
C     A web-based Server which solves unconstrained nonlinear optimization
C     problems using this Conjugate Gradient code can be found at:
C
C       http://www-neos.mcs.anl.gov/neos/solvers/UCO:CGPLUS/
C
C     --------------------------------------------------------------------
C
      SUBROUTINE CGFAM(N,X,F,G,D,GOLD,IPRINT,EPS,W,
     *                  IFLAG,IREST,METHOD,FINISH,ITER,NFUN)
C
C Subroutine parameters
C
      DOUBLE PRECISION X(N),G(N),D(N),GOLD(N),W(N),F,EPS
      INTEGER N,IPRINT(2),IFLAG,IREST,METHOD,IM,NDES
C
C     N      =  NUMBER OF VARIABLES
C     X      =  ITERATE
C     F      =  FUNCTION VALUE
C     G      =  GRADIENT VALUE
C     GOLD   =  PREVIOUS GRADIENT VALUE
C     IPRINT =  FREQUENCY AND TYPE OF PRINTING
C               IPRINT(1) < 0 : NO OUTPUT IS GENERATED
C               IPRINT(1) = 0 : OUTPUT ONLY AT FIRST AND LAST ITERATION
C               IPRINT(1) > 0 : OUTPUT EVERY IPRINT(1) ITERATIONS
C               IPRINT(2)     : SPECIFIES THE TYPE OF OUTPUT GENERATED;
C                               THE LARGER THE VALUE (BETWEEN 0 AND 3),
C                               THE MORE INFORMATION
C               IPRINT(2) = 0 : NO ADDITIONAL INFORMATION PRINTED
C 		IPRINT(2) = 1 : INITIAL X AND GRADIENT VECTORS PRINTED
C		IPRINT(2) = 2 : X VECTOR PRINTED EVERY ITERATION
C		IPRINT(2) = 3 : X VECTOR AND GRADIENT VECTOR PRINTED 
C				EVERY ITERATION 
C     EPS    =  CONVERGENCE CONSTANT
C     W      =  WORKING ARRAY
C     IFLAG  =  CONTROLS TERMINATION OF CODE, AND RETURN TO MAIN
C               PROGRAM TO EVALUATE FUNCTION AND GRADIENT
C               IFLAG = -3 : IMPROPER INPUT PARAMETERS
C               IFLAG = -2 : DESCENT WAS NOT OBTAINED
C               IFLAG = -1 : LINE SEARCH FAILURE
C               IFLAG =  0 : INITIAL ENTRY OR 
C                            SUCCESSFUL TERMINATION WITHOUT ERROR   
C               IFLAG =  1 : INDICATES A RE-ENTRY WITH NEW FUNCTION VALUES
C               IFLAG =  2 : INDICATES A RE-ENTRY WITH A NEW ITERATE
C     IREST  =  0 (NO RESTARTS); 1 (RESTART EVERY N STEPS)
C     METHOD =  1 : FLETCHER-REEVES 
C               2 : POLAK-RIBIERE
C               3 : POSITIVE POLAK-RIBIERE ( BETA=MAX{BETA,0} )
C
C Local variables
C
      DOUBLE PRECISION GTOL,ONE,ZERO,GNORM,DDOT,STP1,FTOL,XTOL,STPMIN,
     .       STPMAX,STP,BETA,BETAFR,BETAPR,DG0,GG,GG0,DGOLD,
     .       DGOUT,DG,DG1
      INTEGER MP,LP,ITER,NFUN,MAXFEV,INFO,I,NFEV,NRST,IDES
      LOGICAL NEW,FINISH
C
C     THE FOLLOWING PARAMETERS ARE PLACED IN COMMON BLOCKS SO THEY
C     CAN BE EASILY ACCESSED ANYWHERE IN THE CODE
C
C     MP = UNIT NUMBER WHICH DETERMINES WHERE TO WRITE REGULAR OUTPUT
C     LP = UNIT NUMBER WHICH DETERMINES WHERE TO WRITE ERROR OUPUT
C      COMMON /CGDD/MP,LP
C
C     ITER: KEEPS TRACK OF THE NUMBER OF ITERATIONS
C     NFUN: KEEPS TRACK OF THE NUMBER OF FUNCTION/GRADIENT EVALUATIONS
C      COMMON /RUNINF/ITER,NFUN
      SAVE
      DATA ONE,ZERO/1.0D+0,0.0D+0/
C
C IFLAG = 1 INDICATES A RE-ENTRY WITH NEW FUNCTION VALUES
      IF(IFLAG.EQ.1) GO TO 72
C
C IFLAG = 2 INDICATES A RE-ENTRY WITH A NEW ITERATE
      IF(IFLAG.EQ.2) GO TO 80
C
C     INITIALIZE
C     ----------
C
C
C     IM =   NUMBER OF TIMES BETAPR WAS NEGATIVE FOR METHOD 2 OR
C            NUMBER OF TIMES BETAPR WAS 0 FOR METHOD 3
C
C     NDES = NUMBER OF LINE SEARCH ITERATIONS AFTER WOLFE CONDITIONS
C            WERE SATISFIED
C
      ITER= 0
      IF(N.LE.0) GO TO 96
      NFUN= 1
      NEW=.TRUE.
      NRST= 0
      IM=0
      NDES=0
C
      DO 5 I=1,N
 5    D(I)= -G(I)
      GNORM= DSQRT(DDOT(N,G,1,G,1))
      STP1= ONE/GNORM
C
C     PARAMETERS FOR LINE SEARCH ROUTINE
C     ----------------------------------
C
C     FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. TERMINATION
C       OCCURS WHEN THE SUFFICIENT DECREASE CONDITION AND THE
C       DIRECTIONAL DERIVATIVE CONDITION ARE SATISFIED.
C
C     XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
C       WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
C       IS AT MOST XTOL.
C
C     STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH
C       SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP.
C
C     MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
C       OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
C       MAXFEV BY THE END OF AN ITERATION.

      FTOL= 1.0D-4
      GTOL= 1.0D-1
      IF(GTOL.LE.1.D-04) THEN
        IF(6.GT.0) WRITE(6,145)
        GTOL=1.D-02
      ENDIF
      XTOL= 1.0D-17
      STPMIN= 1.0D-20
      STPMAX= 1.0D+20
      MAXFEV= 40
C
      IF(IPRINT(1).GE.0) CALL CGBD(IPRINT,ITER,NFUN,
     *   GNORM,N,X,F,G,STP,FINISH,NDES,IM,BETAFR,BETAPR,BETA)
C
C     MAIN ITERATION LOOP
C    ---------------------
C
 8    ITER= ITER+1
C     WHEN NRST>N AND IREST=1 THEN RESTART
      NRST= NRST+1
      INFO=0
C
C
C     CALL THE LINE SEARCH ROUTINE OF MOR'E AND THUENTE
C     (modified for our CG method)
C     -------------------------------------------------
C
C       JJ Mor'e and D Thuente, "Linesearch Algorithms with Guaranteed
C       Sufficient Decrease". ACM Transactions on Mathematical
C       Software 20 (1994), pp 286-307.
C
      NFEV=0
      DO 70 I=1,N
  70  GOLD(I)= G(I)
      DG= DDOT(N,D,1,G,1)
      DGOLD=DG
      STP=ONE
C
C Shanno-Phua's Formula For Trial Step
C
      IF(.NOT.NEW) STP= DG0/DG
      IF (ITER.EQ.1) STP=STP1
      IDES=0
      new=.false.
  72  CONTINUE
C
C     write(6,*) 'step= ', stp
C
C Call to the line search subroutine
C
      CALL CVSMOD(N,X,F,G,D,STP,FTOL,GTOL,
     *            XTOL,STPMIN,STPMAX,MAXFEV,INFO,NFEV,W,DG,DGOUT)

C       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
C         INFO = 0  IMPROPER INPUT PARAMETERS.
C         INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.
C         INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
C                   DIRECTIONAL DERIVATIVE CONDITION HOLD.
C         INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
C                   IS AT MOST XTOL.
C         INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
C         INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
C         INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
C         INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
C                   THERE MAY NOT BE A STEP WHICH SATISFIES THE
C                   SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
C                   TOLERANCES MAY BE TOO SMALL.

      IF (INFO .EQ. -1) THEN
C       RETURN TO FETCH FUNCTION AND GRADIENT
        IFLAG=1
        RETURN
      ENDIF
      IF (INFO .NE. 1) GO TO 90
C
C     TEST IF DESCENT DIRECTION IS OBTAINED FOR METHODS 2 AND 3
C     ---------------------------------------------------------
C
      GG= DDOT(N,G,1,G,1)
      GG0= DDOT(N,G,1,GOLD,1)
      BETAPR= (GG-GG0)/GNORM**2
      IF (IREST.EQ.1.AND.NRST.GT.N) THEN
        NRST=0
        NEW=.TRUE.
        GO TO 75
      ENDIF 
C
      IF (METHOD.EQ.1) THEN
        GO TO 75
      ELSE
        DG1=-GG + BETAPR*DGOUT
        IF (DG1.lt. 0.0d0 ) GO TO 75
        IF (IPRINT(1).GE.0) write(6,*) 'no descent'
        IDES= IDES + 1
        IF(IDES.GT.5) GO TO 95
        GO TO 72
      ENDIF
C
C     DETERMINE CORRECT BETA VALUE FOR METHOD CHOSEN
C     ----------------------------------------------
C
C     IM =   NUMBER OF TIMES BETAPR WAS NEGATIVE FOR METHOD 2 OR
C            NUMBER OF TIMES BETAPR WAS 0 FOR METHOD 3
C
C     NDES = NUMBER OF LINE SEARCH ITERATIONS AFTER WOLFE CONDITIONS
C            WERE SATISFIED
C
  75  NFUN= NFUN + NFEV
      NDES= NDES + IDES
      BETAFR= GG/GNORM**2
      IF (NRST.EQ.0) THEN
        BETA= ZERO
      ELSE
        IF (METHOD.EQ.1) BETA=BETAFR
        IF (METHOD.EQ.2) BETA=BETAPR
        IF ((METHOD.EQ.2.OR.METHOD.EQ.3).AND.BETAPR.LT.0) IM=IM+1
        IF (METHOD.EQ.3) BETA=MAX(ZERO,BETAPR)
      ENDIF
C
C     COMPUTE THE NEW DIRECTION
C     --------------------------
C
      DO 78 I=1,N
  78  D(I) = -G(I) +BETA*D(I)
      DG0= DGOLD*STP
C
C     RETURN TO DRIVER FOR TERMINATION TEST
C     -------------------------------------
C
      GNORM=DSQRT(DDOT(N,G,1,G,1))
      IFLAG=2
      RETURN

  80  CONTINUE
C
C Call subroutine for printing output
C
      IF(IPRINT(1).GE.0) CALL CGBD(IPRINT,ITER,NFUN,
     *     GNORM,N,X,F,G,STP,FINISH,NDES,IM,BETAFR,BETAPR,BETA)
      IF (FINISH) THEN
         IFLAG = 0
         RETURN
      END IF
      GO TO 8
C
C     ----------------------------------------
C     END OF MAIN ITERATION LOOP. ERROR EXITS.
C     ----------------------------------------
C
  90  IFLAG=-1
      IF(6.GT.0) WRITE(6,100) INFO
      RETURN
  95  IFLAG=-2
      IF(6.GT.0) WRITE(6,135) I
      RETURN
  96  IFLAG= -3
      IF(6.GT.0) WRITE(6,140)
C
C     FORMATS
C     -------
C
 100  FORMAT(/' IFLAG= -1 ',/' LINE SEARCH FAILED. SEE'
     .          ' DOCUMENTATION OF ROUTINE CVSMOD',/' ERROR RETURN'
     .          ' OF LINE SEARCH: INFO= ',I2,/
     .          ' POSSIBLE CAUSE: FUNCTION OR GRADIENT ARE INCORRECT')
 135  FORMAT(/' IFLAG= -2',/' DESCENT WAS NOT OBTAINED')
 140  FORMAT(/' IFLAG= -3',/' IMPROPER INPUT PARAMETERS (N',
     .       ' IS NOT POSITIVE)')
 145  FORMAT(/'  GTOL IS LESS THAN OR EQUAL TO 1.D-04',
     .       / ' IT HAS BEEN RESET TO 1.D-02')
      RETURN
      END
C
C     LAST LINE OF ROUTINE CGFAM
C     ***************************
C
C
C**************************************************************************
      SUBROUTINE CGBD(IPRINT,ITER,NFUN,
     *           GNORM,N,X,F,G,STP,FINISH,NDES,IM,BETAFR,BETAPR,BETA)
C
C     ---------------------------------------------------------------------
C     THIS ROUTINE PRINTS MONITORING INFORMATION. THE FREQUENCY AND AMOUNT
C     OF OUTPUT ARE CONTROLLED BY IPRINT.
C     ---------------------------------------------------------------------
C
      DOUBLE PRECISION X(N),G(N),F,GNORM,STP,BETAFR,BETAPR,BETA
      INTEGER IPRINT(2),ITER,NFUN,LP,MP,N,NDES,IM,I
      LOGICAL FINISH
C      COMMON /CGDD/MP,LP
C
      IF (ITER.EQ.0)THEN
           PRINT*
           WRITE(6,10)
           WRITE(6,20) N
           WRITE(6,30) F,GNORM
                 IF (IPRINT(2).GE.1)THEN
                     WRITE(6,40)
                     WRITE(6,50) (X(I),I=1,N)
                     WRITE(6,60)
                     WRITE(6,50) (G(I),I=1,N)
                 ENDIF
           WRITE(6,10)
           WRITE(6,70)
      ELSE
          IF ((IPRINT(1).EQ.0).AND.(ITER.NE.1.AND..NOT.FINISH))RETURN
          IF (IPRINT(1).NE.0)THEN
               IF(MOD(ITER-1,IPRINT(1)).EQ.0.OR.FINISH)THEN
                     IF(IPRINT(2).GT.1.AND.ITER.GT.1) WRITE(6,70)
                     WRITE(6,80)ITER,NFUN,F,GNORM,STP,BETA
               ELSE
                     RETURN
               ENDIF
          ELSE
               IF( IPRINT(2).GT.1.AND.FINISH) WRITE(6,70)
               WRITE(6,80)ITER,NFUN,F,GNORM,STP,BETA
          ENDIF
          IF (IPRINT(2).EQ.2.OR.IPRINT(2).EQ.3)THEN
                  WRITE(6,40)
                  WRITE(6,50)(X(I),I=1,N)
              IF (IPRINT(2).EQ.3)THEN
                  WRITE(6,60)
                  WRITE(6,50)(G(I),I=1,N)
              ENDIF
          ENDIF
          IF (FINISH) WRITE(6,100)
      ENDIF
C
 10   FORMAT('*************************************************')
 20   FORMAT(' N=',I5,//,'INITIAL VALUES:')
 30   FORMAT(' F= ',1PD10.3,'   GNORM= ',1PD10.3)
 40   FORMAT(/,' VECTOR X= ')
 50   FORMAT(6(2X,1PD10.3/))
 60   FORMAT(' GRADIENT VECTOR G= ')
 70   FORMAT(/'   I  NFN',4X,'FUNC',7X,'GNORM',6X,
     *   'STEPLEN',4x,'BETA',/,
     *   ' ----------------------------------------------------')
 80   FORMAT(I4,1X,I3,2X,2(1PD10.3,2X),1PD8.1,2x,1PD8.1)
100   FORMAT(/' SUCCESSFUL CONVERGENCE (NO ERRORS).'
     *          ,/,' IFLAG = 0')
C
      RETURN
      END
C     
C     LAST LINE OF CGBD
C*************************************************************************


