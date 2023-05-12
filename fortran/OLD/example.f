C *Examples:
C     The following is a simple example problem, with the coding needed
C     for its solution by SLSODE. The problem is from chemical kinetics,
C     and consists of the following three rate equations:
C
C        dy1/dt = -.04*y1 + 1.E4*y2*y3
C        dy2/dt = .04*y1 - 1.E4*y2*y3 - 3.E7*y2**2
C        dy3/dt = 3.E7*y2**2
C
C     on the interval from t = 0.0 to t = 4.E10, with initial conditions
C     y1 = 1.0, y2 = y3 = 0. The problem is stiff.
C
C     The following coding solves this problem with SLSODE, using 
C     MF = 21 and printing results at t = .4, 4., ..., 4.E10.  It uses 
C     ITOL = 2 and ATOL much smaller for y2 than for y1 or y3 because y2 
C     has much smaller values.  At the end of the run, statistical 
C     quantities of interest are printed.
C
       EXTERNAL  FEX, JEX
       INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(23), LIW, LRW,
     *         MF, NEQ
       REAL  ATOL(3), RTOL, RWORK(58), T, TOUT, Y(3)
       NEQ = 3
       Y(1) = 1
       Y(2) = 0
       Y(3) = 0
       T = 0
       TOUT = .4
       ITOL = 2
       RTOL = 1.E-4
       ATOL(1) = 1.E-6
       ATOL(2) = 1.E-10
       ATOL(3) = 1.E-6
       ITASK = 1
       ISTATE = 1
       IOPT = 0
       LRW = 58
       LIW = 23
       MF = 21
       DO 40 IOUT = 1,12
       CALL SLSODE (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     *               ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
       WRITE(6,20)  T, Y(1), Y(2), Y(3)
   20  FORMAT(' At t =',E12.4,'   y =',3E14.6)
       IF (ISTATE .LT. 0)  GO TO 80
   40  TOUT = TOUT*10.
       WRITE(6,60)  IWORK(11), IWORK(12), IWORK(13)
   60  FORMAT(/' No. steps =',i4,',  No. f-s =',i4,',  No. J-s =',i4)
       STOP
   80  WRITE(6,90)  ISTATE
   90  FORMAT(///' Error halt.. ISTATE =',I3)
       STOP
       END

       SUBROUTINE  FEX (NEQ, T, Y, YDOT)
       INTEGER  NEQ
       REAL  T, Y(3), YDOT(3)
       YDOT(1) = -.04*Y(1) + 1.E4*Y(2)*Y(3)
       YDOT(3) = 3.E7*Y(2)*Y(2)
       YDOT(2) = -YDOT(1) - YDOT(3)
       RETURN
       END

       SUBROUTINE  JEX (NEQ, T, Y, ML, MU, PD, NRPD)
       INTEGER  NEQ, ML, MU, NRPD
       REAL  T, Y(3), PD(NRPD,3)
       PD(1,1) = -.04
       PD(1,2) = 1.E4*Y(3)
       PD(1,3) = 1.E4*Y(2)
       PD(2,1) = .04
       PD(2,3) = -PD(1,3)
       PD(3,2) = 6.E7*Y(2)
       PD(2,2) = -PD(1,2) - PD(3,2)
       RETURN
       END

