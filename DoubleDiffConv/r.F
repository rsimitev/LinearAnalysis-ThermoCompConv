C-----------------------------------------------------------------------
      FUNCTION R(TRII,NEI,II,JI,KI)
C-----------------------------------------------------------------------
C     BERECHNET DIE RADIALINTEGRALE FUER DAS DYNAMOPROBLEM
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      CHARACTER*3 TRI,TRII
      PARAMETER(DPI=3.141592653589793D0)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
C
      TRI=TRII
      NE=NEI
      I=II
      J=JI
      K=KI
      CALL TAUSCH(TRI,I,J,K)
      IF(NE.EQ.0) THEN
         IF(TRI.EQ.'SS ')THEN 
            IF(I*J.EQ.0) THEN
               RINT=0D0
      	    ELSE             
               RINT=1D0/2D0 * (C0(I-J)-C0(I+J))
            ENDIF
         ELSEIF(TRI.EQ.'SC ') THEN
            IF(I.EQ.0) THEN
      	       RINT=0D0
      	    ELSE
      	       RINT=1D0/2D0 * (S0(I-J)+S0(I+J))
      	    ENDIF
         ELSEIF(TRI.EQ.'CC ') THEN
            RINT=1D0/2D0 * (C0(I+J)+C0(I-J))
         ELSEIF(TRI.EQ.'SSS') THEN
            IF(I*J*K.EQ.0) THEN
      	       RINT=0D0
      	    ELSE
      	       RINT=1D0/4D0 * (S0(I-J+K)+S0(I+J-K)-S0(I-J-K)-S0(I+J+K))
      	    ENDIF
         ELSEIF(TRI.EQ.'SSC') THEN
            IF(I*J.NE.0) THEN
      	       RINT=1D0/4D0 * (C0(I-J-K)-C0(I+J+K)+C0(I-J+K)-C0(I+J-K))
      	    ELSE
               RINT=0D0
      	    ENDIF
         ELSEIF(TRI.EQ.'SCC') THEN
            IF(I.EQ.0) THEN
               RINT=0D0
      	    ELSE
      	       RINT=1D0/4D0 * (S0(I+J-K)+S0(I+J+K)+S0(I-J+K)+S0(I-J-K))
      	    ENDIF
         ELSEIF(TRI.EQ.'CCC') THEN
            RINT= 1D0/4D0 * (C0(I+J+K)+C0(I+J-K)+C0(I-J+K)+C0(I-J-K))
         ENDIF           
      ELSEIF(NE.EQ.1) THEN
         IF(TRI.EQ.'SS ')THEN 
            IF(I*J.EQ.0) THEN
               RINT=0D0
      	    ELSE             
               RINT=1D0/2D0 * (C1(I-J)-C1(I+J))
            ENDIF
         ELSEIF(TRI.EQ.'SC ') THEN
            IF(I.EQ.0) THEN
      	       RINT=0D0
      	    ELSE
      	       RINT=1D0/2D0 * (S1(I-J)+S1(I+J))
      	    ENDIF
         ELSEIF(TRI.EQ.'CC ') THEN
            RINT=1D0/2D0 * (C1(I+J)+C1(I-J))
         ELSEIF(TRI.EQ.'SSS') THEN
            IF(I*J*K.EQ.0) THEN
      	       RINT=0D0
      	    ELSE
      	       RINT=1D0/4D0 * (S1(I-J+K)+S1(I+J-K)-S1(I-J-K)-S1(I+J+K))
      	    ENDIF
         ELSEIF(TRI.EQ.'SSC') THEN
            IF(I*J.NE.0) THEN
      	       RINT=1D0/4D0 * (C1(I-J-K)-C1(I+J+K)+C1(I-J+K)-C1(I+J-K))
      	    ELSE
               RINT=0D0
      	    ENDIF
         ELSEIF(TRI.EQ.'SCC') THEN
            IF(I.EQ.0) THEN
               RINT=0D0
      	    ELSE
      	       RINT=1D0/4D0 * (S1(I+J-K)+S1(I+J+K)+S1(I-J+K)+S1(I-J-K))
      	    ENDIF
         ELSEIF(TRI.EQ.'CCC') THEN
            RINT= 1D0/4D0 * (C1(I+J+K)+C1(I+J-K)+C1(I-J+K)+C1(I-J-K))
         ENDIF
      ELSEIF(NE.EQ.2) THEN
         IF(TRI.EQ.'SS ')THEN 
            IF(I*J.EQ.0) THEN
               RINT=0D0
      	    ELSE             
               RINT=1D0/2D0 * (C2(I-J)-C2(I+J))
            ENDIF
         ELSEIF(TRI.EQ.'SC ') THEN
            IF(I.EQ.0) THEN
      	       RINT=0D0
      	    ELSE
      	       RINT=1D0/2D0 * (S2(I-J)+S2(I+J))
      	    ENDIF
         ELSEIF(TRI.EQ.'CC ') THEN
            RINT=1D0/2D0 * (C2(I+J)+C2(I-J))
         ELSEIF(TRI.EQ.'SSS') THEN
            IF(I*J*K.EQ.0) THEN
      	       RINT=0D0
      	    ELSE
      	       RINT=1D0/4D0 * (S2(I-J+K)+S2(I+J-K)-S2(I-J-K)-S2(I+J+K))
      	    ENDIF
         ELSEIF(TRI.EQ.'SSC') THEN
            IF(I*J.NE.0) THEN
      	       RINT=1D0/4D0 * (C2(I-J-K)-C2(I+J+K)+C2(I-J+K)-C2(I+J-K))
      	    ELSE
               RINT=0D0
      	    ENDIF
         ELSEIF(TRI.EQ.'SCC') THEN
            IF(I.EQ.0) THEN
               RINT=0D0
      	    ELSE
      	       RINT=1D0/4D0 * (S2(I+J-K)+S2(I+J+K)+S2(I-J+K)+S2(I-J-K))
      	    ENDIF
         ELSEIF(TRI.EQ.'CCC') THEN
            RINT= 1D0/4D0 * (C2(I+J+K)+C2(I+J-K)+C2(I-J+K)+C2(I-J-K))
         ENDIF           
      ELSEIF(NE.EQ.3) THEN
         IF(TRI.EQ.'SS')THEN 
            IF(I*J.EQ.0) THEN
               RINT=0D0
      	    ELSE             
               RINT=1D0/2D0 * (C3(I-J)-C3(I+J))
            ENDIF
         ELSEIF(TRI.EQ.'SC') THEN
            IF(I.EQ.0) THEN
      	       RINT=0D0
      	    ELSE
      	       RINT=1D0/2D0 * (S3(I-J)+S3(I+J))
      	    ENDIF
         ELSEIF(TRI.EQ.'CC') THEN
            RINT=1D0/2D0 * (C3(I+J)+C3(I-J))
         ELSEIF(TRI.EQ.'SSS') THEN
            IF(I*J*K.EQ.0) THEN
      	       RINT=0D0
      	    ELSE
      	       RINT=1D0/4D0 * (S3(I-J+K)+S3(I+J-K)-S3(I-J-K)-S3(I+J+K))
      	    ENDIF
         ELSEIF(TRI.EQ.'SSC') THEN
            IF(I*J.NE.0) THEN
      	       RINT=1D0/4D0 * (C3(I-J-K)-C3(I+J+K)+C3(I-J+K)-C3(I+J-K))
      	    ELSE
               RINT=0D0
      	    ENDIF
         ELSEIF(TRI.EQ.'SCC') THEN
            IF(I.EQ.0) THEN
               RINT=0D0
      	    ELSE
      	       RINT=1D0/4D0 * (S3(I+J-K)+S3(I+J+K)+S3(I-J+K)+S3(I-J-K))
      	    ENDIF
         ELSEIF(TRI.EQ.'CCC') THEN
            RINT= 1D0/4D0 * (C3(I+J+K)+C3(I+J-K)+C3(I-J+K)+C3(I-J-K))
         ENDIF           
      ELSEIF(NE.EQ.4) THEN
         IF(TRI.EQ.'SS')THEN 
            IF(I*J.EQ.0) THEN
               RINT=0D0
      	    ELSE             
               RINT=1D0/2D0 * (C4(I-J)-C4(I+J))
            ENDIF
         ELSEIF(TRI.EQ.'SC') THEN
            IF(I.EQ.0) THEN
      	       RINT=0D0
      	    ELSE
      	       RINT=1D0/2D0 * (S4(I-J)+S4(I+J))
      	    ENDIF
         ELSEIF(TRI.EQ.'CC') THEN
            RINT=1D0/2 * (C4(I+J)+C4(I-J))
         ELSEIF(TRI.EQ.'SSS') THEN
            IF(I*J*K.EQ.0) THEN
      	       RINT=0D0
      	    ELSE
      	       RINT=1D0/4D0 * (S4(I-J+K)+S4(I+J-K)-S4(I-J-K)-S4(I+J+K))
      	    ENDIF
         ELSEIF(TRI.EQ.'SSC') THEN
            IF(I*J.NE.0) THEN
      	       RINT=1D0/4D0 * (C4(I-J-K)-C4(I+J+K)+C4(I-J+K)-C4(I+J-K))
      	    ELSE
               RINT=0D0
      	    ENDIF
         ELSEIF(TRI.EQ.'SCC') THEN
            IF(I.EQ.0) THEN
               RINT=0D0
      	    ELSE
      	       RINT=1D0/4D0 * (S4(I+J-K)+S4(I+J+K)+S4(I-J+K)+S4(I-J-K))
      	    ENDIF
         ELSEIF(TRI.EQ.'CCC') THEN
            RINT= 1D0/4D0 * (C4(I+J+K)+C4(I+J-K)+C4(I-J+K)+C4(I-J-K))
         ENDIF           
      ELSEIF(NE.EQ.-1) THEN
         IF(TRI.EQ.'SS ') THEN
            IF(I*J.EQ.0) THEN
               RINT=0D0
      	    ELSE             
               RINT=1D0/2D0 * (CM1(I-J)-CM1(I+J))
            ENDIF
         ELSEIF(TRI.EQ.'SC ') THEN
            IF(I.EQ.0) THEN
      	       RINT=0D0
      	    ELSE
      	       RINT=1D0/2D0 * (SM1(I-J)+SM1(I+J))
      	    ENDIF
         ELSEIF(TRI.EQ.'CC ') THEN
            RINT=1D0/2D0 * (CM1(I+J)+CM1(I-J))
         ELSEIF(TRI.EQ.'SSS') THEN
            IF(I*J*K.EQ.0) THEN
      	       RINT=0D0
      	    ELSE
      	     RINT=1D0/4D0*(SM1(I-J+K)+SM1(I+J-K)-SM1(I-J-K)-SM1(I+J+K))
      	    ENDIF
         ELSEIF(TRI.EQ.'SSC') THEN
            IF(I*J.NE.0) THEN
      	     RINT=1D0/4D0*(CM1(I-J-K)-CM1(I+J+K)+CM1(I-J+K)-CM1(I+J-K))
      	    ELSE
               RINT=0D0
      	    ENDIF
         ELSEIF(TRI.EQ.'SCC') THEN
            IF(I.EQ.0) THEN
               RINT=0D0
      	    ELSE
      	     RINT=1D0/4D0*(SM1(I+J-K)+SM1(I+J+K)+SM1(I-J+K)+SM1(I-J-K))
      	    ENDIF
         ELSEIF(TRI.EQ.'CCC') THEN
           RINT= 1D0/4D0*(CM1(I+J+K)+CM1(I+J-K)+CM1(I-J+K)+CM1(I-J-K))
         ENDIF           
      ELSEIF(NE.EQ.-2) THEN
         IF(TRI.EQ.'SS ')THEN 
            IF(I*J.EQ.0) THEN
               RINT=0D0
      	    ELSE             
               RINT=1D0/2D0 * (CM2(I-J)-CM2(I+J))
            ENDIF
         ELSEIF(TRI.EQ.'SC ') THEN
            IF(I.EQ.0) THEN
      	       RINT=0D0
      	    ELSE
      	       RINT=1D0/2D0 * (SM2(I-J)+SM2(I+J))
      	    ENDIF
         ELSEIF(TRI.EQ.'CC ') THEN
            RINT=1D0/2D0 * (CM2(I+J)+CM2(I-J))
         ELSEIF(TRI.EQ.'SSS') THEN
            IF(I*J*K.EQ.0) THEN
      	       RINT=0D0
      	    ELSE
      	     RINT=1D0/4D0*(SM2(I-J+K)+SM2(I+J-K)-SM2(I-J-K)-SM2(I+J+K))
      	    ENDIF
         ELSEIF(TRI.EQ.'SSC') THEN
            IF(I*J.NE.0) THEN
      	     RINT=1D0/4D0*(CM2(I-J-K)-CM2(I+J+K)+CM2(I-J+K)-CM2(I+J-K))
      	    ELSE
               RINT=0D0
      	    ENDIF
         ELSEIF(TRI.EQ.'SCC') THEN
            IF(I.EQ.0) THEN
               RINT=0D0
      	    ELSE
      	     RINT=1D0/4D0*(SM2(I+J-K)+SM2(I+J+K)+SM2(I-J+K)+SM2(I-J-K))
      	    ENDIF
         ELSEIF(TRI.EQ.'CCC') THEN
          RINT= 1D0/4D0*(CM2(I+J+K)+CM2(I+J-K)+CM2(I-J+K)+CM2(I-J-K))
         ENDIF           
      ENDIF
      R=RINT
C      WRITE(12,'(I4,A3,3I4,D14.4)') NE,TRI,I,J,K,R
      END
C-----------------------------------------------------------------------
C     END OF RI
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION NGU(N)
      NGU=-1
      IF( MOD(N,2).EQ.0 ) NGU=1      
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      SUBROUTINE TAUSCH(TRI,I,J,K)
      IMPLICIT REAL*8(A-H,O-Y)
      CHARACTER*3 TRI
      IF(TRI.EQ.'CS ') THEN
         N=I
         I=J
         J=N
         TRI='SC '
      ELSEIF(TRI.EQ.'SCS') THEN
         N=J
         J=K
         K=N
         TRI='SSC'
      ELSEIF(TRI.EQ.'CSS') THEN
         N=I
         I=K
         K=N
         TRI='SSC'
      ELSEIF(TRI.EQ.'CCS') THEN
         N=I
         I=K
         K=N
         TRI='SCC'
      ELSEIF(TRI.EQ.'CSC') THEN
         N=I
         I=J
         J=N
         TRI='SCC'
      ENDIF            
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION S0(N)
      IMPLICIT REAL*8(A-H,O-Y)
      PARAMETER(DPI=3.141592653589793D0)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      IF(NGU(N).EQ.1) THEN
         S0=0D0
      ELSE
         S0=2D0/N/DPI
      ENDIF
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION C0(N)
      IMPLICIT REAL*8(A-H,O-Y)
      PARAMETER(DPI=3.141592653589793D0)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      IF(N.EQ.0) THEN
         C0=1D0
      ELSE
         C0=0D0
      ENDIF
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION S1(N)
      IMPLICIT REAL*8(A-H,O-Y)
      PARAMETER(DPI=3.141592653589793D0)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      IF(NGU(N).EQ.1) THEN
         IF(N.EQ.0) THEN
      	    S1=0D0
      	 ELSE
            S1=-1D0/N/DPI
         ENDIF
      ELSE
         S1=(1+2*RI)/N/DPI
      ENDIF
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION C1(N)
      IMPLICIT REAL*8(A-H,O-Y)
      PARAMETER(DPI=3.141592653589793D0)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      IF(NGU(N).EQ.1) THEN
         IF(N.EQ.0) THEN
            C1=RI+1D0/2
         ELSE
            C1=0D0
         ENDIF
      ELSE
         C1=-2D0/N**2/DPI**2
      ENDIF
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION S2(N)
      IMPLICIT REAL*8(A-H,O-Y)
      PARAMETER(DPI=3.141592653589793D0)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      IF(NGU(N).EQ.1) THEN
      	 IF(N.EQ.0) THEN
            S2=0D0
         ELSE
            S2=-(2*RI+1)/DPI/N
         ENDIF
      ELSE
         S2=(1+2*RI*(RI+1))/DPI/N - 4D0/N**3/DPI**3
      ENDIF
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION C2(N)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y)
      PARAMETER(DPI=3.141592653589793D0)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      IF(NGU(N).EQ.1) THEN
         IF(N.EQ.0) THEN
            C2=RI*(RI+1)+1D0/3
         ELSE
            C2=2D0/N**2/DPI**2 
         ENDIF
      ELSE
         C2=-2D0*(2*RI+1)/N**2/DPI**2
      ENDIF
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION S3(N)
      IMPLICIT REAL*8(A-H,O-Y)
      PARAMETER(DPI=3.141592653589793D0)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      IF(NGU(N).EQ.1) THEN
      	 IF(N.EQ.0) THEN
            S3=0D0
         ELSE
            S3=-(1+3*RI+3*RI**2)/DPI/N + 6/DPI**3/N**3
         ENDIF
      ELSE
         S3=(1+3*RI+3*RI**2+2*RI**3)/DPI/N - (6+12*RI)/DPI**3/N**3
      ENDIF
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION C3(N)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y)
      PARAMETER(DPI=3.141592653589793D0)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      IF(NGU(N).EQ.1) THEN
         IF(N.EQ.0) THEN
            C3=1D0/4+RI+3D0/2*RI**2+RI**3
         ELSE
            C3=(3+6*RI)/DPI**2/N**2
         ENDIF
      ELSE
         C3=-(3+6*RI+6*RI**2)/DPI**2/N**2+12/DPI**4/N**4
      ENDIF
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION S4(N)
      IMPLICIT REAL*8(A-H,O-Y)
      PARAMETER(DPI=3.141592653589793D0)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      IF(NGU(N).EQ.1) THEN
      	 IF(N.EQ.0) THEN
            S4=0D0
         ELSE
            S4=-(1+4*RI+6*RI**2+4*RI**3)/DPI/N+(12+24*RI)/DPI**3/N**3
         ENDIF
      ELSE
         S4=(1+4*RI+6*RI**2+4*RI**3+2*RI**4)/DPI/N - 
     -	     (12+24*RI+24*RI**2)/DPI**3/N**3 + 48D0/DPI**5/N**5
      ENDIF
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION C4(N)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y)
      PARAMETER(DPI=3.141592653589793D0)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      IF(NGU(N).EQ.1) THEN
         IF(N.EQ.0) THEN
            C4=1D0/5+RI+2*RI**2+2*RI**3+RI**4
         ELSE
            C4=(4+12*RI+12*RI**2)/DPI**2/N**2-24/DPI**4/N**4
         ENDIF
      ELSE
         C4=-(4+12*RI+12*RI**2+8*RI**3)/DPI**2/N**2 + 
     +       (24+48*RI)/DPI**4/N**4 
      ENDIF
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION SM1(N)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y)
      PARAMETER(DPI=3.141592653589793D0)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      IF(N.EQ.0) THEN
         SM1=0D0
      ELSE
         SM1= DCOS(N*RI*DPI)*DIS(RI,RI+1,DPI*N) - 
     -        DSIN(N*RI*DPI)*DIC(RI,RI+1,N*DPI)
      ENDIF
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION CM1(N)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y)
      PARAMETER(DPI=3.141592653589793D0)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      IF(N.EQ.0) THEN
         CM1=DLOG((RI+1)/RI)
      ELSE
         CM1= DCOS(N*RI*DPI)*DIC(RI,RI+1,N*DPI) + 
     +        DSIN(N*RI*DPI)*DIS(RI,RI+1,N*DPI)
      ENDIF
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION SM2(N)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y)
      PARAMETER(DPI=3.141592653589793D0)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      IF(N.EQ.0) THEN
         SM2=0D0
      ELSE
      SM2=DCOS(N*RI*DPI)*( DSIN(N*DPI*RI)/RI-DSIN(N*DPI*(RI+1))/(RI+1) +
     + 			    N*DPI*DIC(RI,RI+1,N*DPI) ) -
     -    DSIN(N*RI*DPI)*( DCOS(N*DPI*RI)/RI-DCOS(N*DPI*(RI+1))/(RI+1) -
     -			    N*DPI*DIS(RI,RI+1,N*DPI) )
      ENDIF
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION CM2(N)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y)
      PARAMETER(DPI=3.141592653589793D0)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      IF(N.EQ.0) THEN
         CM2=1D0/RI-1D0/(1+RI)
      ELSE
      CM2=DSIN(N*RI*DPI)*( DSIN(N*DPI*RI)/RI-DSIN(N*DPI*(RI+1))/(RI+1) +
     + 			    N*DPI*DIC(RI,RI+1,N*DPI) ) +
     +    DCOS(N*RI*DPI)*( DCOS(N*DPI*RI)/RI-DCOS(N*DPI*(RI+1))/(RI+1) -
     -			    N*DPI*DIS(RI,RI+1,N*DPI) )
      ENDIF
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION DIS(XMIN,XMAX,A)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y)
C
      X=DABS(A*XMAX)
      IF(X.LT.1) THEN
         DISMAX=SIS(X)
      ELSE
         DISMAX=SIA(X)
      ENDIF
      IF(A*XMAX.LT.0) DISMAX=-DISMAX
      X=DABS(A*XMIN)
      IF(X.LT.1) THEN
         DISMIN=SIS(X)
      ELSE
         DISMIN=SIA(X)
      ENDIF
      IF(A*XMIN.LT.0) DISMIN=-DISMIN
      DIS= DISMAX-DISMIN
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION SIA(X)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y)
      DIMENSION AF(4),BF(4),AG(4),BG(4)
      PARAMETER(DPI=3.141592653589793D0)
      DATA AF/38.027264D0,265.187033D0,335.677320D0,38.102495D0/
     &     BF/40.021433D0,322.624911D0,570.236280D0,157.105423D0/
     &     AG/42.242855D0,302.757865D0,352.018498D0,21.821899D0/
     &     BG/48.196927D0,482.485984D0,1114.978885D0,449.690326D0/
      F= (X**8+AF(1)*X**6+AF(2)*X**4+AF(3)*X**2+AF(4)) /
     /   (X**8+BF(1)*X**6+BF(2)*X**4+BF(3)*X**2+BF(4)) / X
      G= (X**8+AG(1)*X**6+AG(2)*X**4+AG(3)*X**2+AG(4)) /
     /   (X**8+BG(1)*X**6+BG(2)*X**4+BG(3)*X**2+BG(4)) / X**2
      SIA=DPI/2D0 - F*DCOS(X) - G*DSIN(X)
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION SIS(X)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y)
      GENAU=1D-7
C      
      SISO=X
      DO 200 I=1,200000
         JZ=0
         JN=0
         SISZ=1D0
         SISN=2D0*I+1D0
10       DO 20 J=JZ+1,2*I+1
            SISZ=SISZ*X
            IF(SISZ.GT.1D20) THEN
               GOTO 30
            ENDIF
20       CONTINUE
30	 JZ=J
         DO 40 J=JN+1,2*I+1
            SISN=SISN*J
            IF(SISN.GT.1D20) THEN
               GOTO 50
            ENDIF
40       CONTINUE
50	 JN=J
	 SISZ=SISZ/SISN
         SISN=1
         IF((JZ.LT.2*I+1).OR.(JN.LT.2*I+1)) GOTO 10
         SISN=2D0*I+1D0
         JN=0
C         
         SIS=SISO + (-1)**I * SISZ
         IF(I.GT.1) THEN
            IF(DABS(1D0-SIS/SISO).LE.GENAU) GOTO 300
         ENDIF
         SISO=SIS
200   CONTINUE
300   END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION DIC(XMIN,XMAX,A)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y)
C
      X=DABS(A*XMAX)
      IF(X.LT.1) THEN
         DICMAX=CIS(X)
      ELSE
         DICMAX=CIA(X)
      ENDIF
      X=DABS(A*XMIN)
      IF(X.LT.1) THEN
         DICMIN=CIS(X)
      ELSE
         DICMIN=CIA(X)
      ENDIF
      DIC= DICMAX-DICMIN
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION CIA(X)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y)
      DIMENSION AF(4),BF(4),AG(4),BG(4)
      DATA AF/38.027264D0,265.187033D0,335.677320D0,38.102495D0/
     &     BF/40.021433D0,322.624911D0,570.236280D0,157.105423D0/
     &     AG/42.242855D0,302.757865D0,352.018498D0,21.821899D0/
     &     BG/48.196927D0,482.485984D0,1114.978885D0,449.690326D0/
      F= (X**8+AF(1)*X**6+AF(2)*X**4+AF(3)*X**2+AF(4)) /
     /   (X**8+BF(1)*X**6+BF(2)*X**4+BF(3)*X**2+BF(4)) / X
      G= (X**8+AG(1)*X**6+AG(2)*X**4+AG(3)*X**2+AG(4)) /
     /   (X**8+BG(1)*X**6+BG(2)*X**4+BG(3)*X**2+BG(4)) / X**2
      CIA=F*DSIN(X) - G*DCOS(X)
      END
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
      FUNCTION CIS(X)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y)
      PARAMETER(C=0.5772156649D0)
      GENAU=1D-7
C      
      CISO=DLOG(X)+C
      DO 200 I=1,200000
         JZ=0
         JN=0
         CISZ=1D0
         CISN=2D0*I
10       DO 20 J=JZ+1,2*I
            CISZ=CISZ*X
            IF(CISZ.GT.1D20) THEN
               GOTO 30
            ENDIF
20       CONTINUE
30	 JZ=J
         DO 40 J=JN+1,2*I
            CISN=CISN*J
            IF(CISN.GT.1D20) THEN
               GOTO 50
            ENDIF
40       CONTINUE
50	 JN=J
	 CISZ=CISZ/CISN
         CISN=1D0
         IF((JZ.LT.2*I).OR.(JN.LT.2*I)) GOTO 10
C         
         CIS=CISO + (-1)**I * CISZ
         IF(I.GT.1) THEN
            IF(DABS(1D0-CIS/CISO).LE.GENAU) GOTO 300
         ENDIF
         CISO=CIS
200   CONTINUE
300   END
C-----------------------------------------------------------------------
