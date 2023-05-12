      INTEGER  K, N, L, NQ, NYH, IFLAG
       REAL  T, TN, H, HU, UROUND, DKY(3), YH(3,2)
      COMMON /SLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   IOWND(6), IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU

C     Ininitalize arguments
       UROUND = RUMACH()
       K = 2
       N = 3
       L = 2
       NQ = 4
       IFLAG = 0
       NYH = N
       DKY(1) = 1.0E0
       DKY(2) = 0.0E0
       DKY(3) = 0.0E0
       T = 0.4
       TN = T
       H = 1.0E0
       HU = 0.0E0
       YH(1,1) = 1.0E0
       YH(1,2) = 0.0E0
       YH(2,1) = 0.0E0
       YH(2,2) = 1.0E0
       YH(3,1) = 1.0E0
       YH(3,2) = 1.0E0

C     Run SINTDY
       CALL SINTDY (T, K, YH, NYH, DKY, IFLAG)

C     Print DKY
       WRITE(6,20)  DKY(1), DKY(2), DKY(3)
   20  FORMAT('   DKY =',3E14.6)
       STOP
       END
