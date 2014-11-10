**********************************************************************
      SUBROUTINE PEGASUS(X,DX,ABSE,RELE,NMAX,LM,LO,LN,X0)
**********************************************************************
C--  PEGASUS-Verfahren bei lm = 1 ,
C--  ANDERSEN-BJOERNSEN-KING-Verfahren bei lm = 2 .
C--  X   : startvalue for x
C--  DX  : initial stepsize in x
C--  ABSE: absolut error limit
C--  RELE: relativ error limit
C--  NMAX: max. count of steps
C--  LM  : (1/2) chooses algorithm (see above)
C--  LO  : (0/1) print nothing/print result of each step
C--  LN  : on output: 0=convergence, 1:no convergence
C--  X0  : on output: calculated root
C--
C--  PEGASUS calls an external function named
C--  REAL*8 FUNCTION FN(RA) to determine the function value
C---------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-G,O-Z)
C
      IF( ABSE.NE.0.D0 .AND. RELE.NE.0.D0 ) THEN
         WRITE(*,*) 'ABSERR OR RELERR SHOULD BE ZERO.'
         STOP
      ENDIF
C
C--   assuming the startvalue X is already fairly well:
      DX=MAX(DX,X/100)
      DX=MIN(DX,X)
C
      X1=X
      X2=X+DX
      F1=FN(X1)
      F2=FN(X2)
      N=2  
      IF( LO.EQ.1 ) THEN
        WRITE(*,'(A)') ' PEGASUS.F: STEP   1'
        WRITE(*,10010) X1 ,F1
        WRITE(*,'(A)') ' PEGASUS.F: STEP   2'
        WRITE(*,10010) X2 ,F2
      ENDIF

C
5100  FF=F1*F2
C
      NO=1
      IF( FF.GT.0.D0 ) THEN
         ND=1
         DO 5200 I=3,100
C
            ND=-ND
            X2=X1 + INT(I/2 + 0.001) * DX * ND
            IF( X2.LE.0.D0) THEN
               IF(NO.EQ.1) THEN
                  X2=0.d0
      	          NO=-1
               ELSEIF(NO.EQ.-1) THEN
                  GOTO 5200
               ENDIF
            ENDIF

            F2=FN(X2)
            N=N+1
            IF( LO.EQ.1 ) THEN
               WRITE(*,'(A,I4)') ' PEGASUS.F: STEP',N
               WRITE(*,10010) X2 ,F2
            ENDIF
C
            FF=F1*F2
            IF( FF.LT.0 ) THEN
C
               IF( I.GT.3 ) THEN 
      	          X1=X1 + ( INT(I/2 + 0.001) - 1 ) * DX * ND
      		  F1=FN(X1)
                  N=N+1
                 IF( LO.EQ.1 ) THEN
                    WRITE(*,'(A,I4)') ' PEGASUS.F: STEP',N
                    WRITE(*,10010) X2 ,F2
                 ENDIF
               ENDIF
C
               IF( X2.LT.X1 ) THEN
                  CALL PTAUSCH(X2,F2,X1,F1)
               ENDIF
C
               GOTO 5100
C
            ENDIF
C
            IF( N.GT.NMAX ) THEN
               LN=1
               GOTO 20000
            ENDIF
C
5200     CONTINUE
C
      ELSE
C
         IF( DABS(X1-X2).LT.(Dabs(X1)*RELE+ABSE) ) THEN
            X0=DABS(X1)
            LN=0 
            GOTO 20000
         ENDIF
C
         IF( N.GT.NMAX ) THEN
            LN=1
            GOTO 20000
         ENDIF
C
         FD=F2-F1
         X3=X2-F2*(X2-X1)/FD
c
         F3=FN(X3)
         N=N+1
C
         IF( LO.EQ.1 ) THEN
            WRITE(*,'(A,I4)') ' PEGASUS.F: STEP',N
            WRITE(*,10010) X3,F3
         ENDIF
         FF=F2*F3
         IF( FF.LT.0.D0 ) THEN
c
      	    CALL PTAUSCH(X1,F1,X2,F2)            
      	    CALL PTAUSCH(X2,F2,X3,F3)            
c
         ELSE
C
            IF( LM.EQ.1 ) THEN
C----------- PEGASUSVERFAHREN 
               QF=F2/(F2+F3)
            ELSEIF( LM.EQ.2 ) THEN
c----------- ABK-VERFAHREN
      	       QF=1.D0-F2/F3
               IF( QF.LT.0.D0 ) QF=0.5D0
            ENDIF
C
            F1=F1*QF
C
            CALL PTAUSCH(X2,F2,X3,F3)
C
	 ENDIF
      ENDIF
C
      GOTO 5100
C
10010 FORMAT(1X,' X = ',G20.13,'  F(X) = ',G20.13)
C
20000 continue
      RETURN
      END
C
C--------------------------------------------------------------
C-- END OF PEGASUS 
C--------------------------------------------------------------
c
c
****************************************************************
      SUBROUTINE PTAUSCH(X1,Y1,X2,Y2)
****************************************************************
C--  Vertauscht Werte 1 und 2 
C---------------------------------------------------------------
C
      IMPLICIT REAL*8(A-G,O-Z)
C
      XH=X1
      YH=Y1
      X1=X2
      Y1=Y2
      X2=XH
      Y2=YH
C
      RETURN
      END
C
C----------------------------------------------------------------------
C-- END OF PTAUSCH 
c---------------------------------------------------------------------






