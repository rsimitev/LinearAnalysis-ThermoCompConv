***************************************************************
      SUBROUTINE READX(STARTFILE,NUDSR,TIMER,XR,LIC,LCALCI)
*****************************************************************
C
C----------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CF,CFR,CFS
      CHARACTER*2 CRR,CRRR
      CHARACTER*40 STARTFILE
      PARAMETER (NM=5000,NMT=50)
C
      DIMENSION XR(*),LR(NM),MR(NM),NR(NM),KR(NM),CFR(NM),CRRR(NM)
C
      COMMON/LOG/LCALC,LWRITE,L3,LVAR,LDY,LT,L7,LREAD,L9,L10
      COMMON/NUM/RELE,EPS,ALPH,STOER,NITMAX,NJA
      COMMON/TINT/DT,TIME0,NST,NTW,NDTW,NDTPW,NKWT,IKWT(NMT)
C
      COMMON/DIMI/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/DIMSI/NDVS,NDWS,NDTS,NDHS,NDGS,NDS
      COMMON/PAR/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/PARI/RAI,TAI,PRI,PMI,ETAI,CI,OMI,FTWI,FTGI,MFI
      COMMON/NPAR/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD
      COMMON/NPARI/M0I,NTVI,NTHI,LTVI,LTHI,KTVI,KTHI,LEVI,LRBI,LDI
C
      COMMON/QNUVI/NUCI,NUOMI
      COMMON/QNUSI/NMQC,LS(NM),MS(NM),NS(NM),CFS(NM)
      COMMON/QNUI/NMC,L(NM),M(NM),N(NM),K(NM),CF(NM),CRR(NM)
      COMMON/CALLS/NCNONLIN,NCLIN,NCSETBL,NCINPUT
C
      NCINPUT=NCINPUT+1
C
      IF( NMC.NE.NM .OR. NMQC.NE.NM ) THEN
	 WRITE(*,*) 'WRONG DIMENSION NM IN READX: ',NM,NMC
	 STOP
      ENDIF
C
C-- COUNTER , DIMENSIONS AND LOGICALS SET TO ZERO:
C
      OPEN(12,FILE=STARTFILE,STATUS='old')
      REWIND(12)
C
C-- LST IS FORMAT PARAMETER ( L=1 FOR HIRSCHING FORMAT ) ,
C                           ( LST=0 FOR OLD FORMAT ) ,
C                           ( LST=2 FOR NEW FORMAT ) .
C   LCALCI TELLS HOW THE INPUTFILE HAS BEEN CALCULATED:
      READ(12,*) LST,LCALCI
C
C-- IF LTR.EQ.1 THE INPUTFILE IS TIMEINTEGRATION:
      IF( LCALCI.EQ.5 .OR. LCALCI.EQ.6 ) THEN
         LTR=1
      ELSE
         LTR=0
      ENDIF
C
      IF( LST.EQ.0 .OR. LST.EQ.2 ) THEN
C
C-- READH READS THE HEADER OF THE INPUTFILE AND DETERMINS WETHER
C   THE DATASET (NUDS,TIME) IS THE RIGHT ONE (LDR=1):
         CALL READH(12,LTR,NUDSR,TIMER,NUDS,TIME,LDR)
         IF( MFI.NE.0 ) THEN
            WRITE(*,*) 'SORRY , FLOQUET PAR. OF INPUT HAS TO BE 0.'
            STOP
         ENDIF
C
C-- LOOKING FOR THE RIGHT DATASET:
         DO 1000 I=1,1000
C
C----- READD READS FULL SET OF COEFFITIENTS:
            CALL READD(12,LDR,NK,XR,CFR,CRRR,LR,MR,NR,KR,
     & 	         EVPM,EVPF,EVTM,EVTF,DNU,EMPM,EMPF,EMTM,EMTF)
C
            IF( LDR.EQ.1 ) GOTO 1200
C
C---- READP READS THE DATASETNUMBER, TIME IF GIVEN (LT=1) AND 
C     NEW PARAMETERS FOR LT.NE.1 , LDR IS DETERMINED:
            CALL READP(12,LST,LTR,NUDSR,TIMER,NUDS,TIME,LDR,LEND)
            IF( LEND.EQ.1 ) GOTO 1100
C
1000     CONTINUE
1100     CONTINUE
         WRITE(*,'('' SORRY, SET '',I4,D14.4,'' NOT FOUND'')')
     &		  NUDSR,TIMER
         STOP
1200     CONTINUE

      ELSEIF( LST.EQ.1 ) THEN
C
C-- READW READS PARAMETERS AND SET OF KOEFFITIENTS FOR HIRSCHING INPUT:
         CALL READW(12,NK,XR,CFR,CRRR,LR,MR,NR,KR)
	 LCALCI=LCALC
C
      ELSE
         WRITE(*,*) 'WRONG FORMAT OF INPUTFILE.'
         STOP
      ENDIF
C
      CALL GETSYM(LCALCI,NK,CFR,LR,MR,M0I,LEVI,LRBI,LDI)
      CALL INPUTCHECK(LCALCI,LIC,TIME0)
C
      IF( MF.NE.0 ) THEN
C
C----- LOOKING FOR M<0:
         DO 1300 I=1,NK
1300     IF( MR(I).LT.0 ) GOTO 1500
C
C----- PRODUCING OF M<0:
         NDX=NK
         DO 1400 I=1,NDX
            IF( MR(I).NE.0 ) THEN
               NK=NK+1
               IF( CRRR(I).EQ.'RR' .OR. CRRR(I).EQ.'RI' ) THEN
                  XR(NK)=XR(I)
               ELSEIF( CRRR(I).EQ.'IR' .OR. CRRR(I).EQ.'II' ) THEN
                  XR(NK)=-XR(I)
               ENDIF
               CFR(NK)=CFR(I)
               CRRR(NK)=CRRR(I)
               LR(NK)=LR(I)
               MR(NK)=-MR(I)
               NR(NK)=NR(I)
               KR(NK)=KR(I)
            ENDIF
1400     CONTINUE
      ENDIF
1500  CONTINUE   
C
C-- SORT INPUT:
      LSX=1
      CALL SORTK(NK,LSX,XR,CFR,CRRR,LR,MR,NR,KR,NUCI,NUOMI)
C
C-- DETERMINE DIMENSION OF INPUT AND POSITION OF VARIABLE C AND OM:
      ND=NK
      IF( ND.GT.NM ) THEN
	 WRITE(*,*) 'TOO SMALL DIMENSION NM IN READX.',NM,ND
	 STOP
      ENDIF
      DO 1600 I=1,ND
         CF(I)=CFR(I)
	 CRR(I)=CRRR(I)
	 L(I)=LR(I)
	 M(I)=MR(I)
	 N(I)=NR(I)
	 K(I)=KR(I)
1600  CONTINUE
      CALL RDIM(ND,CF,CRR,L,M,N,K,CFS,LS,MS,NS,
     &           NDV,NDW,NDT,NDH,NDG,NDVS,NDWS,NDTS,NDHS,NDGS,NDS)
C
C-- FOR KINEMATIC DYNAMO MAGNETIC FIELD OF INPUT IS NOT NEEDED:	
      IF( LCALC.EQ.2 ) THEN
      	 FTW=FTWI
      	 ND=ND-NDH-NDG
         NDH=0
         NDG=0
         NDS=NDS-NDHS-NDGS
         NDHS=0
         NDGS=0
      ENDIF
C
      CALL KOEFF(1)
C
      IF( DABS(CI).GT.0.D0 ) C=CI
      IF( DABS(OMI).GT.0.D0) OM=OMI
      IF( KTV.EQ.0 .AND. KTH.EQ.0 ) OM=0.D0
C
      CLOSE(12)
C
      RETURN
      END
C
C-------------------------------------------------------------------------
C-- END OF SUBROUTINE READX
C-------------------------------------------------------------------------
C
C
*****************************************************************
      SUBROUTINE READH(NOUT,LTR,NUDSR,TIMER,NUDS,TIME,LDR)
*****************************************************************
C
C----------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*2 C2
      PARAMETER(NMT=50)
C
      COMMON/LOG/LCALC,LWRITE,L3,LVAR,LDY,LT,L7,L8,L9,L10
      COMMON/PARI/RAI,TAI,PRI,PMI,ETAI,CI,OMI,FTWI,FTGI,MFI
      COMMON/NPARI/M0I,NTVI,NTHI,LTVI,LTHI,KTVI,KTHI,LEVI,LRBI,LDI
      COMMON/TINT/DT,TIME0,NST,NTW,NDTW,NDTPW,NKWT,IKWT(NMT)
      COMMON/QNUVI/NUCI,NUOMI
C
      READ(NOUT,*) M0I,NTVI,NTHI,KTVI,KTHI,LTVI,LTHI,FTWI,FTGI
      READ(NOUT,*) NUDS,TAI,RAI,PRI,PMI,ETAI
      READ(NOUT,*) CI,OMI,NUCI,NUOMI,MFI
      TAI=DSQRT(TAI)
C
      IF( LTR.EQ.1 ) THEN
         READ(NOUT,9100) C2,NUDS,TIME
      ENDIF
      IF( ( LT.EQ.1 .AND. TIME.EQ.TIMER ) .OR.
     &	  ( LT.EQ.0 .AND. NUDS.EQ.NUDSR )  ) THEN
         LDR=1
         IF( LTR.EQ.1 ) THEN
            TIME0=TIME
         ELSE
            TIME0=0.D0
         ENDIF
      ELSE
         LDR=0
      ENDIF
C
9100  FORMAT(1X,A2,I6,D16.6)
C
      RETURN
      END
C
C-------------------------------------------------------------------------
C-- END OF SUBROUTINE READH
C-------------------------------------------------------------------------
C
C
*****************************************************************
      SUBROUTINE READD(NOUT,LDR,NKR,XR,CFR,CRRR,LR,MR,NR,KR,
     &	               EVPM,EVPF,EVTM,DNU,EVTF,EMPM,EMPF,EMTM,EMTF)
*****************************************************************
C
C----------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CFI,CFR
      CHARACTER*2 C2,CRRR
C
      DIMENSION XR(*),LR(*),MR(*),NR(*),KR(*),CFR(*),CRRR(*)
C
      COMMON/LOG/LCALC,LWRITE,L3,LVAR,LDY,LT,L7,LREAD,L9,L10
      COMMON/NPARI/M0I,NTVI,NTHI,LTVI,LTHI,KTVI,KTHI,LEVI,LRBI,LDI
C
      NKR=0
      DO 1000 I=1,5000
	 READ(NOUT,9100,END=1100,ERR=1100) C2
	 IF( C2.EQ.'V ' .OR. C2.EQ.'W ' .OR. C2.EQ.'T ' .OR.
     &       C2.EQ.'G ' .OR. C2.EQ.'H ' ) THEN 

            IF( LDR.EQ.1 ) THEN
               BACKSPACE(NOUT)
               READ(NOUT,9300) CFI,LI,MI,NI,KI,DRR,DIR,DRI,DII
               NKR=NKR+1
               CFR(NKR)=CFI
               CRRR(NKR)='RR'
      	       LR(NKR)=LI
      	       MR(NKR)=MI
      	       NR(NKR)=NI
      	       KR(NKR)=KI
      	       XR(NKR)=DRR
               IF( MI.NE.0 ) THEN
                  NKR=NKR+1
                  CFR(NKR)=CFI
                  CRRR(NKR)='IR'
      	          LR(NKR)=LI
      		  MR(NKR)=MI
      	          NR(NKR)=NI
      	          KR(NKR)=KI
      		  XR(NKR)=DIR
               ENDIF
               IF( KI.NE.0 ) THEN
                  NKR=NKR+1
                  CFR(NKR)=CFI
                  CRRR(NKR)='RI'
      	          LR(NKR)=LI
      		  MR(NKR)=MI
      	          NR(NKR)=NI
      	          KR(NKR)=KI
      		  XR(NKR)=DRI
               ENDIF
               IF( MI.NE.0 .AND. KI.NE.0 ) THEN
                  NKR=NKR+1
                  CFR(NKR)=CFI
                  CRRR(NKR)='II'
      	          LR(NKR)=LI
      		  MR(NKR)=MI
      	          NR(NKR)=NI
      	          KR(NKR)=KI
      		  XR(NKR)=DII
               ENDIF
            ENDIF                 
	 ELSEIF( C2.EQ.'EV' ) THEN
            IF( LDR.EQ.1 ) THEN
               BACKSPACE(NOUT)
               READ(NOUT,9400) C2,EVPM,EVPF,EVTM,EVTF,DNU
     	    ENDIF
	 ELSEIF( C2.EQ.'EM' ) THEN
            IF( LDR.EQ.1 ) THEN
               BACKSPACE(NOUT) 
               READ(NOUT,9500) C2,EMPM,EMPF,EMTM,EMTF
     	    ENDIF
	 ELSE
	    BACKSPACE(NOUT) 
	    GOTO 1100
	 ENDIF
1000  CONTINUE
1100  CONTINUE
C
9100  FORMAT(1X,A2)
9200  FORMAT(1X,A2,I6,D16.6)
9300  FORMAT(1X,A1,4I3,4D16.8)
9400  FORMAT(1X,A2,5D14.6)
9500  FORMAT(1X,A2,4D14.6)
C 
      RETURN
      END
C
C-------------------------------------------------------------------------
C-- END OF SUBROUTINE READD
C-------------------------------------------------------------------------
C
C
*****************************************************************
      SUBROUTINE READP(NOUT,LST,LTR,NUDSR,TIMER,NUDS,TIME,LDR,LEND)
*****************************************************************
C
C----------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*2 C2
      PARAMETER(NMT=50)
C
      COMMON/LOG/LCALC,LWRITE,L3,LVAR,LDY,LT,L7,L8,L9,L10
      COMMON/PARI/RAI,TAI,PRI,PMI,ETAI,CI,OMI,FTWI,FTGI,MFI
      COMMON/NPARI/M0I,NTVI,NTHI,LTVI,LTHI,KTVI,KTHI,LEVI,LRBI,LDI
      COMMON/TINT/DT,TIME0,NST,NTW,NDTW,NDTPW,NKWT,IKWT(NMT)
      COMMON/QNUVI/NUCI,NUOMI
C
      IF( LTR.EQ.1 ) THEN
         READ(NOUT,9100,END=2000) C2,NUDS,TIME
      ELSEIF( LTR.EQ.0 ) THEN
         IF( LST.EQ.2 ) READ(NOUT,*,END=2000) 
     &       M0I,NTVI,NTHI,KTVI,KTHI,LTVI,LTHI,FTWI,FTGI
         READ(NOUT,*,END=2000) NUDS,TAI,RAI,PRI,PMI,ETAI
         READ(NOUT,*) CI,OMI,NUCI,NUOMI
         TAI=DSQRT(TAI)
      ENDIF
      IF( ( LT.EQ.1 .AND. TIME.EQ.TIMER ) .OR.
     &	  ( LT.EQ.0 .AND. NUDS.EQ.NUDSR )  ) THEN
         LDR=1
         IF( LTR.EQ.1 ) THEN
            TIME0=TIME
         ELSE
            TIME0=0.D0
      	 ENDIF
      ELSE
         LDR=0
      ENDIF
C
      LEND=0
      RETURN
C
2000  CONTINUE
      LEND=1
C
9100  FORMAT(1X,A2,I6,D16.6)
C
      RETURN
      END
C
C-------------------------------------------------------------------------
C-- END OF SUBROUTINE READp
C-------------------------------------------------------------------------
C
C
*****************************************************************
      SUBROUTINE READW(NOUT,NXD,XR,CFR,CRRR,LR,MR,NR,KR)
*****************************************************************
C
C----------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CFR,CFI
      CHARACTER*2 CRRR
C
      DIMENSION XR(*),LR(*),MR(*),NR(*),KR(*),CFR(*),CRRR(*)
C
      COMMON/LOG/LCALC,LWRITE,LTR,LVAR,LDY,LT,L7,LREAD,L9,L10
      COMMON/PARI/RAI,TAI,PRI,PMI,ETAI,CI,OMI,FTWI,FTGI,MFI
      COMMON/NPARI/M0I,NTVI,NTHI,LTVI,LTHI,KTVI,KTHI,LEVI,LRBI,LDI
      COMMON/QNUVI/NUCI,NUOMI
C
      WRITE(*,*) 'CALL OF READW TO READ HIRSCHING INPUT.'
C
      NXD=0
      L2B=0
      NCR=0
      NCI=0
      NUOMI=0
      NVMAX=1
      NHMAX=1
C
      ETAI=0.4D0
      READ(NOUT,*) M0I,I1,I2,NC,SR,SI
      READ(NOUT,*) TAI,RAI,PRI,PMI
      DO 3000 I=1,3000
         READ(NOUT,*,END=3100,ERR=3100) NF,LI,NUI,NI,DRR,DIR
         IF( LI.EQ.0 .AND. NUI.EQ.0 .AND. NI.EQ.0 ) GOTO 3100
      	 IF( NF.EQ.1 ) THEN
       	    CFI='V'
      	 ELSEIF( NF.EQ.3 ) THEN
       	    CFI='W'
      	 ELSEIF( NF.EQ.5 ) THEN
       	    CFI='T'
      	 ELSEIF( NF.EQ.7 ) THEN
       	    CFI='H'
      	 ELSEIF( NF.EQ.9 ) THEN
       	    CFI='G'
         ENDIF
      	 DRI=0.D0
         DII=0.D0
         MI=IABS(NUI)*M0I
      	 KI=0
         IF( I.EQ.1 ) L1=LI
         IF( LI.NE.L1 .AND. L2B.EQ.0 ) THEN
      	    L2=LI
            L2B=1
      	 ENDIF
         IF(  NI.GT.NVMAX .AND. ( NF.EQ.1 .OR. NF.EQ.3 ) ) THEN
      	    NVMAX=NI
         ELSEIF(  NI.GT.NHMAX .AND. ( NF.EQ.7 .OR. NF.EQ.9 ) ) THEN
      	    NHMAX=NI
         ENDIF
      	 NXD=NXD+1
         CFR(NXD)=CFI
         CRRR(NXD)='RR'
         XR(NXD)=DRR
         LR(NXD)=LI
      	 MR(NXD)=MI
      	 NR(NXD)=NI
      	 KR(NXD)=KI
         IF( NF.EQ.1 ) THEN
      	    NCR=NCR+1
      	 ENDIF
         IF( MI.NE.0 ) THEN
      	    NXD=NXD+1
            CFR(NXD)=CFI
            CRRR(NXD)='IR'
            XR(NXD)=DIR
      	    LR(NXD)=LI
      	    MR(NXD)=MI
      	    NR(NXD)=NI
      	    KR(NXD)=KI
            IF( NF.EQ.1 ) THEN
      	       NCI=NCI+1
      	    ENDIF
         ENDIF
         IF( KI.NE.0 ) THEN
      	    NXD=NXD+1
            CFR(NXD)=CFI
      	    CRRR(NXD)='RI'
            XR(NXD)=DRI
      	    LR(NXD)=LI
      	    MR(NXD)=MI
      	    NR(NXD)=NI
      	    KR(NXD)=KI
         ENDIF
         IF( MI.NE.0 .AND. KI.NE.0 ) THEN
      	    NXD=NXD+1
            CFR(NXD)=CFI
      	    CRRR(NXD)='II'
            XR(NXD)=DII
      	    LR(NXD)=LI
      	    MR(NXD)=MI
      	    NR(NXD)=NI
      	    KR(NXD)=KI
         ENDIF
3000  CONTINUE
3100  CONTINUE
C
      TAI=DSQRT(TAI)
C
      NUCI=2*NC-NCR-NCI
      CI=XR(NUCI)
C
      WRITE(*,*) 'DRIFT IN HIRSCHING FILE: ',CI
C
      IF( LCALC.EQ.2 ) XR(NUCI)=0.D0
      OMI=0.D0
      NTVI=NVMAX-1
      NTHI=NHMAX-1
      KTVI=0
      KTHI=0
      LTVI=0
      LTHI=0
C
      RETURN
      END
C
C-------------------------------------------------------------------------
C-- END OF SUBROUTINE READW
C-------------------------------------------------------------------------
C
C
*************************************************************************
      SUBROUTINE READKH(NOUT,NK,NUDS,TIME,CFR,CRRR,LR,MR,NR)
*************************************************************************
C
C------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NKM=5)
      CHARACTER*1 CFR,CFI
      CHARACTER*2 C2,CRRR,CRRI
C
      DIMENSION CFR(*),CRRR(*),LR(*),MR(*),NR(*)
      DIMENSION CFI(NKM),CRRI(NKM),LI(NKM),MI(NKM),NI(NKM)
C
      COMMON/PAR/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPAR/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD
C
      READ(NOUT,*) M0,NTV,NTH,KTV,KTH,LTV,LTH,FTW,FTG
      READ(NOUT,*) NUST,TA,RA,PR,PM,ETA
      READ(NOUT,*) C,OM,NUC,NUOM,MF
      READ(NOUT,*) NK
C
      NKZ=(NK+4)/5
C
      NKR=0
      DO 1000 I=1,NKZ
         READ(NOUT,9100) CK,(CFI(J),CRRI(J),LI(J),MI(J),NI(J),J=1,NKM)
         DO 100 J=1,NKM
            NKR=NKR+1
            CFR(NKR)=CFI(J)
      	    CRRR(NKR)=CRRI(J)
            LR(NKR)=LI(J)
            MR(NKR)=MI(J)
      	    NR(NKR)=NI(J)
100	 CONTINUE
1000  CONTINUE
C
      READ(NOUT,'(1X,A2,I8,D16.6)') C2,NUDS,TIME
C
9100  FORMAT(1X,A2,1X,A1,1X,A2,3I3,2X,A1,1X,A2,3I3,2X,A1,1X,A2,3I3,2X,
     &	        A1,1X,A2,3I3,2X,A1,1X,A2,3I3)
C
      RETURN
      END
C
C----------------------------------------------------------------------
C  END OF SUBROUTINE READKH
C----------------------------------------------------------------------
C
C
*************************************************************************
      SUBROUTINE READKX(NOUT,NK,XR,LEM,
     &			EVPM,EVPF,EVTM,EVTF,DNU,EMPM,EMPF,EMTM,EMTF)
*************************************************************************
C
C------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*2 C2
      PARAMETER (NKM=5)
C
      DIMENSION XR(*),XI(NKM)
C
      NKZ=(NK+4)/5
      NKO=0
C
      DO 1000 I=1,NKZ
      	 READ(NOUT,'(1X,A2,D14.6,4D15.6)') C2,(XI(J),J=1,NKM)
         DO 100 J=1,NKM
      	    NKO=NKO+1
            XR(NKO)=XI(J)
100	 CONTINUE
1000  CONTINUE
C
      READ(NOUT,9100) C2,EVPM,EVPF,EVTM,EVTF,DNU
      IF( LEM.EQ.1 ) THEN
         READ(NOUT,9200) C2,EMPM,EMPF,EMTM,EMTF
      ENDIF
C
9100  FORMAT(1X,A2,5D14.6)
9200  FORMAT(1X,A2,4D14.6)
C
C
      RETURN
      END
C
C----------------------------------------------------------------------
C  END OF SUBROUTINE READK            
C----------------------------------------------------------------------
C
C
*************************************************************************
      SUBROUTINE READKP(NOUT,NUDS,TIME,LEND)
*************************************************************************
C
C------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*2 C2
C
      READ(NOUT,'(1X,A2,I8,D16.6)',END=1000) C2,NUDS,TIME
C
      LEND=0
      RETURN
C
1000  CONTINUE
      LEND=1
      RETURN
C
      END
C
C----------------------------------------------------------------------
C  END OF SUBROUTINE READKP
C----------------------------------------------------------------------
C
C
*******************************************************************
      SUBROUTINE ROTX(X,CRR,M,K,NX,ND)
*******************************************************************
C-- ROTATES X SO THAT X(NX) = 0 .
C------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*2 CRR
      DIMENSION X(*),CRR(*),M(*),K(*)
C
      IF( DABS(X(NX)).LT.1.D-6 ) RETURN
C
      PI=3.141592654D0
C
      IF( CRR(NX).EQ.'IR' ) THEN
C----- ROTATION IN PHI
         LR=1         
	 MC=M(NX)
         IF( MC.EQ.0 ) THEN
            WRITE(*,*) 'WRONG MC IN ROTX.',MC,NX
            STOP
         ENDIF
      ELSEIF( CRR(NX).EQ.'RI' ) THEN
C----- ROTATION IN TIME
         LR=2
	 KOM=K(NX)
         IF( KOM.EQ.0 ) THEN
            WRITE(*,*) 'WRONG KOM IN ROTX.',KOM,NX
            STOP
         ENDIF
      ELSE
         WRITE(*,*) 'NO ROTATION POSSIBLE IN ROTX.'
         STOP
      ENDIF
C
      IF( LR.EQ.1 ) THEN
         XR=X(NX-1)
         XI=X(NX)
         IF( DABS(XR).LT.1.D-6 ) THEN
            PHI=PI/2/MC
         ELSE
	    PHI=DATAN(XI/XR)/MC
         ENDIF
      ELSEIF( LR.EQ.2 ) THEN
         XR=X(NX-2)
         XI=X(NX)
         IF( DABS(XR).LT.1.D-6 ) THEN
            PHI=PI/2/KOM
         ELSE
	    PHI=DATAN(XI/XR)/KOM
         ENDIF
      ENDIF
C
      NXC=0
      DO 100 I=1,ND
         NXC=NXC+1
      	 NXRR=NXC
      	 NXROT=NXC
         XRR=X(NXC)
         XIR=0.D0
         XRI=0.D0
         XII=0.D0
	 MX=M(NXRR)
	 KX=K(NXRR)
         IF( MX.NE.0 ) THEN
            NXC=NXC+1
            XIR=X(NXC)
         ENDIF
         IF( KX.NE.0 ) THEN
            NXC=NXC+1
            XRI=X(NXC)
         ENDIF
         IF( MX.NE.0 .AND. KX.NE.0 ) THEN
            NXC=NXC+1
            XII=X(NXC)
         ENDIF
	 IF( LR.EQ.1 .AND. MX.NE.0 ) THEN
	    X(NXROT)=XRR*DCOS(MX*PHI)+XIR*DSIN(MX*PHI)
      	    NXROT=NXROT+1
            X(NXROT)=XIR*DCOS(MX*PHI)-XRR*DSIN(MX*PHI)
      	    IF( KX.NE.0 ) THEN
      	       NXROT=NXROT+1
	       X(NXROT)=XRI*DCOS(MX*PHI)+XII*DSIN(MX*PHI)
      	       NXROT=NXROT+1
               X(NXROT)=XII*DCOS(MX*PHI)-XRI*DSIN(MX*PHI)
      	    ENDIF                                     
         ELSEIF( LR.EQ.2 .AND. KX.NE.0 ) THEN
	    X(NXROT)=XRR*DCOS(KX*PHI)+XRI*DSIN(KX*PHI)
	    IF( MX.NE.0 ) THEN
      	       NXROT=NXROT+1
	       X(NXROT)=XIR*DCOS(KX*PHI)+XII*DSIN(KX*PHI)
	    ENDIF
      	    NXROT=NXROT+1
            X(NXROT)=XRI*DCOS(KX*PHI)-XRR*DSIN(KX*PHI)
      	    IF( MX.NE.0 ) THEN
      	       NXROT=NXROT+1
               X(NXROT)=XII*DCOS(KX*PHI)-XIR*DSIN(KX*PHI)
      	    ENDIF        
         ENDIF
         IF( NXC.EQ.ND ) GOTO 200
100   CONTINUE
200   CONTINUE
C
      RETURN
      END
C
C-------------------------------------------------------------------
C
C
*****************************************************************
      SUBROUTINE RDIM(ND,CF,CRR,L,M,N,K,CFS,LS,MS,NS,
     &      NDV,NDW,NDT,NDH,NDG,NDVS,NDWS,NDTS,NDHS,NDGS,NDS)
*****************************************************************
C----------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CF,CFS
      CHARACTER*2 CRR
C
      DIMENSION CF(*),CRR(*),L(*),M(*),N(*),K(*)
      DIMENSION CFS(*),LS(*),MS(*),NS(*)
C
      NDV=0
      NDW=0
      NDT=0
      NDH=0
      NDG=0
      NDVS=0
      NDWS=0
      NDTS=0
      NDHS=0
      NDGS=0
      NDS=0
C
      DO 1000 I=1,ND
C
         IF( CRR(I).EQ.'RR' .AND. K(I).EQ.0 ) THEN
            NDS=NDS+1
            CFS(NDS)=CF(I)
            LS(NDS)=L(I)
            MS(NDS)=M(I)
            NS(NDS)=N(I)
         ENDIF
         IF( CF(I).EQ.'V' ) THEN
            NDV=NDV+1
            IF( CRR(I).EQ.'RR' .AND. K(I).EQ.0 ) NDVS=NDVS+1
         ELSEIF( CF(I).EQ.'W' ) THEN
            NDW=NDW+1
            IF( CRR(I).EQ.'RR' .AND. K(I).EQ.0 ) NDWS=NDWS+1
         ELSEIF( CF(I).EQ.'T' ) THEN
            NDT=NDT+1
            IF( CRR(I).EQ.'RR' .AND. K(I).EQ.0 ) NDTS=NDTS+1
         ELSEIF( CF(I).EQ.'H' ) THEN
            NDH=NDH+1
            IF( CRR(I).EQ.'RR' .AND. K(I).EQ.0 ) NDHS=NDHS+1
         ELSEIF( CF(I).EQ.'G' ) THEN
            NDG=NDG+1
            IF( CRR(I).EQ.'RR' .AND. K(I).EQ.0 ) NDGS=NDGS+1
         ENDIF
C
1000  CONTINUE
C
      RETURN
      END
C
C-------------------------------------------------------------------    
C
C
*****************************************************************
      SUBROUTINE KOEFF(LOU)
*****************************************************************
C----------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CF,CFS,CFI,CFSI
      CHARACTER*2 CRR,CRRI
      PARAMETER(NM=5000,NMT=50)
C
      COMMON/LOG/LCALC,LWRITE,L3,LVAR,LDY,LT,L7,LREAD,L9,L10
      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/DIMS/NDVS,NDWS,NDTS,NDHS,NDGS,NDS
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),CF(NM),CRR(NM)
      COMMON/QNUS/NMQC,LS(NM),MS(NM),NS(NM),CFS(NM)
      COMMON/DIMI/NDVI,NDWI,NDTI,NDHI,NDGI,NDI
      COMMON/DIMSI/NDVSI,NDWSI,NDTSI,NDHSI,NDGSI,NDSI
      COMMON/QNUI/NMCI,LI(NM),MI(NM),NI(NM),KI(NM),CFI(NM),CRRI(NM)
      COMMON/QNUSI/NMSCI,LSI(NM),MSI(NM),NSI(NM),CFSI(NM)
      COMMON/TINT/DT,TIME0,NST,NTW,NDTW,NDTPW,NKWT,IKWT(NMT)
C
      IF( NMC.NE.NM ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN KOEFF.'
         STOP
      ENDIF
C
      IF( LOU.EQ.0 ) THEN
C
         IF( ND.GT.NM ) THEN
            WRITE(*,*) 'DIMENSION NM OF CALCULATION TOO SMALL.',ND
            STOP
         ENDIF
         IF( NDS.GT.NM ) THEN
            WRITE(*,*) 'DIMENSION NMS OF CALCULATION TOO SMALL.',NDS
            STOP
         ENDIF
         IF( LCALC.EQ.2 ) THEN
            WRITE(*,*) 'DIMENSION OF CALCULATION ND = ',NDH+NDG
         ELSE
            WRITE(*,*) 'DIMENSION OF CALCULATION ND = ',ND
         ENDIF
C
         OPEN(21,FILE='koeff.dat',STATUS='unknown')
C
         WRITE(21,*) 'DIMENSIONS OF CALCULATION: '
         WRITE(21,*) 'ND = ',ND
         WRITE(21,*) 'NDV= ',NDV
         WRITE(21,*) 'NDW= ',NDW
         WRITE(21,*) 'NDT= ',NDT
         WRITE(21,*) 'NDH= ',NDH
         WRITE(21,*) 'NDG= ',NDG
         WRITE(21,*) 'NDS = ',NDS
         WRITE(21,*) 'NDVS= ',NDVS
         WRITE(21,*) 'NDWS= ',NDWS
         WRITE(21,*) 'NDTS= ',NDTS
         WRITE(21,*) 'NDHS= ',NDHS
         WRITE(21,*) 'NDGS= ',NDGS
C
C-- AUSSCHREIBEN DER BENUTZTEN KOEFFITIENTEN
         WRITE(21,*) 'BENUTZTE KOEFFITIENTEN FUER RECHNUNG:'
         DO 100 I=1,ND
100      WRITE(21,'(I5,A2,A3,4I4)') I,CF(I),CRR(I),L(I),M(I),N(I),K(I)
         WRITE(21,*)
         WRITE(21,*) 'REDUZIERTE KOEFFITIENTEN:'
         DO 200 I=1,NDS
200      WRITE(21,'(I5,A2,3I4)') I,CFS(I),LS(I),MS(I),NS(I)
C
C-- SETZEN DER NUMMERN DER AUSZUSCHREIBENDEN KOEFFITIENTEN:
         IF( LCALC.EQ.5 .OR. LCALC.EQ.6 ) THEN
            CALL SETWT(ND)
            WRITE(21,*)
            WRITE(21,*) 'NUMMERN DER AUSGESCHRIEBENEN KOEFFITIENTEN:'
            DO 300 I=1,NKWT
 300        WRITE(21,*) I,IKWT(I)      
         ENDIF
         CLOSE(21)
C
      ELSEIF( LOU.EQ.1 ) THEN
C
         IF( NDI.GT.NM ) THEN
            WRITE(*,*) 'DIMENSION NM OF INPUT TOO SMALL.',NDI
            STOP
         ELSE
            WRITE(*,*) 'DIMENSIONS OF INPUT ND = ',NDI
         ENDIF
         IF( NDSI.GT.NM ) THEN
            WRITE(*,*) 'DIMENSION NMS OF INPUT TOO SMALL.',NDSI
            STOP
         ENDIF
C
         OPEN(21,FILE='koeffi.dat',STATUS='unknown')
C
         WRITE(21,*) 'DIMENSIONS OF INPUT: '
         WRITE(21,*) 'ND = ',NDI
         WRITE(21,*) 'NDV= ',NDVI
         WRITE(21,*) 'NDW= ',NDWI
         WRITE(21,*) 'NDT= ',NDTI
         WRITE(21,*) 'NDH= ',NDHI
         WRITE(21,*) 'NDG= ',NDGI
         WRITE(21,*) 'NDS = ',NDSI
         WRITE(21,*) 'NDVS= ',NDVSI
         WRITE(21,*) 'NDWS= ',NDWSI
         WRITE(21,*) 'NDTS= ',NDTSI
         WRITE(21,*) 'NDHS= ',NDHSI
         WRITE(21,*) 'NDGS= ',NDGSI
C
C-- AUSSCHREIBEN DER BENUTZTEN KOEFFITIENTEN
         WRITE(21,*) 'BENUTZTE KOEFFITIENTEN DES INPUTS:'
         DO 400 I=1,NDI
400      WRITE(21,'(I5,A2,A3,4I4)') 
     &			I,CFI(I),CRRI(I),LI(I),MI(I),NI(I),KI(I)
         WRITE(21,*)
         WRITE(21,*) 'REDUZIERTE KOEFFITIENTEN DES INPUTS:'
         DO 500 I=1,NDSI
500      WRITE(21,'(I5,A2,3I4)') I,CFSI(I),LSI(I),MSI(I),NSI(I)
C
      ENDIF
C
      RETURN
      END
C
C---------------------------------------------------------------------
C
C
*****************************************************************
      SUBROUTINE SETWT(ND)
*****************************************************************
C----------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CF,CFWT
      CHARACTER*2 CRR,CRRWT
      PARAMETER (NM=5000,NMT=50)
C
      COMMON/QNU/NMQC,L(NM),M(NM),N(NM),K(NM),CF(NM),CRR(NM)
      COMMON/TINT/DT,TIME0,NST,NTW,NDTW,NDTPW,NKWT,IKWT(NMT)
      COMMON/COEW/NMKTC,NKT,LWT(NMT),MWT(NMT),NWT(NMT),KWT(NMT),
     &		  CFWT(NMT),CRRWT(NMT)
C
      IF( NMT.NE.NMKTC ) THEN 
	 WRITE(*,*) 'WRONG DIMENSION NMT IN SETWT (SK).'
	 STOP
      ENDIF
C
      DO 200 J=1,NKT
      DO 100 I=1,ND
         IF( CF(I).EQ.CFWT(J) .AND. CRR(I).EQ.CRRWT(J) .AND. 
     &	     L(I).EQ.LWT(J) .AND. M(I).EQ.MWT(J) .AND. 
     &			N(I).EQ.NWT(J) .AND. K(I).EQ.KWT(J) )  THEN
      	    NKWT=NKWT+1
      	    IKWT(NKWT)=I
            GOTO 200
         ENDIF
100   CONTINUE
200   CONTINUE
C
      RETURN
      END
C
C--------------------------------------------------------------------
C-- END OF SUBROUTINE SETWT
C--------------------------------------------------------------------
C
C
*****************************************************************
      SUBROUTINE DIMEN
*****************************************************************
C----------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CF,CFQ,CFS
      CHARACTER*2 CRR,CRRQ
      PARAMETER (NM=5000)
C
      DIMENSION X(NM)
C
      COMMON/LOG/LCALC,LWRITE,LTR,LVAR,LDY,L6,L7,L8,L9,L10
      COMMON/NUM/RELE,EPS,ALPH,STOER,NITMAX,NJA
      COMMON/PAR/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPAR/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LDI
      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/DIMS/NDVS,NDWS,NDTS,NDHS,NDGS,NDS
      COMMON/QNUV/NUC,NUOM
C
      COMMON/QNU/NMQC,LQ(NM),MQ(NM),NQ(NM),KQ(NM),CFQ(NM),CRRQ(NM)
      COMMON/QNUS/NMSC,LS(NM),MS(NM),NS(NM),CFS(NM)
C
      IF( NM.NE.NMQC ) THEN 
	 WRITE(*,*) 'WRONG DIMENSION NM IN DIMEN (INPUT).'
	 STOP
      ENDIF
      IF( NM.NE.NMSC ) THEN 
	 WRITE(*,*) 'WRONG DIMENSION NM IN DIMEN (INPUT).'
	 STOP
      ENDIF
C
      ND=0
C
      IF( LCALC.NE.2 ) THEN
         CF='V'
         DO 100 M=MMIN(CF,LRB,M0,MF,NTV,LTV,LTR,LCALC),
     &			MMAX(M0,MF,NTV,LTV,LTR,LCALC),M0
         DO 100 L=LMIN(CF,M,LEV,LCALC),LMAX(M,M0,NTV,LTV,LTR),LDI 
         DO 100 N=1,NMAX(CF,FTW,FTG,L,M,M0,NTV,LTR)
         DO 100 K=0,KTV
            CRR='RR'
      	    CALL SETQU(ND,CF,CRR,L,M,N,K)
            IF( M.NE.0 ) THEN 
               CRR='IR'
      	       CALL SETQU(ND,CF,CRR,L,M,N,K)
            ENDIF
            IF( K.NE.0 ) THEN 
               CRR='RI'
      	       CALL SETQU(ND,CF,CRR,L,M,N,K)
            ENDIF
            IF( M.NE.0 .AND. K.NE.0 ) THEN 
               CRR='II'
      	       CALL SETQU(ND,CF,CRR,L,M,N,K)
            ENDIF
100      CONTINUE
C
         IF( LCALC.EQ.1 .AND. TA.LT.EPS ) GOTO 250
C-------  ONSET IN NON ROTATING CASE HAS NO TOROIDAL COMPONENT.
C
         CF='W'
         DO 200 M=MMIN(CF,LRB,M0,MF,NTV,LTV,LTR,LCALC),
     &			MMAX(M0,MF,NTV,LTV,LTR,LCALC),M0
         DO 200 L=LMIN(CF,M,LEV,LCALC),LMAX(M,M0,NTV,LTV,LTR),LDI
         DO 200 N=1,NMAX(CF,FTW,FTG,L,M,M0,NTV,LTR)
         DO 200 K=0,KTV
         IF( LRIBOUN(CF,L,M,N).EQ.1 ) GOTO 200
            CRR='RR'
      	    CALL SETQU(ND,CF,CRR,L,M,N,K)
            IF( M.NE.0 ) THEN 
               CRR='IR'
      	       CALL SETQU(ND,CF,CRR,L,M,N,K)
            ENDIF
            IF( K.NE.0 ) THEN 
               CRR='RI'
      	       CALL SETQU(ND,CF,CRR,L,M,N,K)
            ENDIF
            IF( M.NE.0 .AND. K.NE.0 ) THEN 
               CRR='II'
      	       CALL SETQU(ND,CF,CRR,L,M,N,K)
            ENDIF
200      CONTINUE
C
250      CONTINUE
C
         CF='T'
         DO 300 M=MMIN(CF,LRB,M0,MF,NTV,LTV,LTR,LCALC),
     &			MMAX(M0,MF,NTV,LTV,LTR,LCALC),M0
         DO 300 L=LMIN(CF,M,LEV,LCALC),LMAX(M,M0,NTV,LTV,LTR),LDI
         DO 300 N=1,NMAX(CF,FTW,FTG,L,M,M0,NTV,LTR)
         DO 300 K=0,KTV
            CRR='RR'
      	    CALL SETQU(ND,CF,CRR,L,M,N,K)
            IF( M.NE.0 ) THEN 
               CRR='IR'
      	       CALL SETQU(ND,CF,CRR,L,M,N,K)
            ENDIF
            IF( K.NE.0 ) THEN 
               CRR='RI'
      	       CALL SETQU(ND,CF,CRR,L,M,N,K)
            ENDIF
            IF( M.NE.0 .AND. K.NE.0 ) THEN 
               CRR='II'
      	       CALL SETQU(ND,CF,CRR,L,M,N,K)
            ENDIF
300      CONTINUE
      ENDIF
C
      IF( LCALC.EQ.2 .OR. LCALC.EQ.4 .OR. LCALC.EQ.6 .OR. 
     &					  LCALC.EQ.8 ) THEN
         CF='H'
         DO 400 M=MMIN(CF,LRB,M0,MF,NTH,LTV,LTR,LCALC),
     &			MMAX(M0,MF,NTH,LTH,LTR,LCALC),M0
         DO 400 L=LMIN(CF,M,LEV,LCALC),LMAX(M,M0,NTH,LTH,LTR),LDI
         DO 400 N=1,NMAX(CF,FTW,FTG,L,M,M0,NTH,LTR)
         DO 400 K=0,KTH
            CRR='RR'
      	    CALL SETQU(ND,CF,CRR,L,M,N,K)
            IF( M.NE.0 ) THEN 
               CRR='IR'
      	       CALL SETQU(ND,CF,CRR,L,M,N,K)
            ENDIF
            IF( K.NE.0 ) THEN 
               CRR='RI'
      	       CALL SETQU(ND,CF,CRR,L,M,N,K)
            ENDIF
            IF( M.NE.0 .AND. K.NE.0 ) THEN 
               CRR='II'
      	       CALL SETQU(ND,CF,CRR,L,M,N,K)
            ENDIF
400      CONTINUE
C
         CF='G'
         DO 500 M=MMIN(CF,LRB,M0,MF,NTH,LTV,LTR,LCALC),
     &			MMAX(M0,MF,NTH,LTH,LTR,LCALC),M0
         DO 500 L=LMIN(CF,M,LEV,LCALC),LMAX(M,M0,NTH,LTH,LTR),LDI
         DO 500 N=1,NMAX(CF,FTW,FTG,L,M,M0,NTH,LTR)
         DO 500 K=0,KTH
            CRR='RR'
      	    CALL SETQU(ND,CF,CRR,L,M,N,K)
            IF( M.NE.0 ) THEN 
               CRR='IR'
      	       CALL SETQU(ND,CF,CRR,L,M,N,K)
            ENDIF
            IF( K.NE.0 ) THEN 
               CRR='RI'
      	       CALL SETQU(ND,CF,CRR,L,M,N,K)
            ENDIF
            IF( M.NE.0 .AND. K.NE.0 ) THEN 
               CRR='II'
      	       CALL SETQU(ND,CF,CRR,L,M,N,K)
            ENDIF
500      CONTINUE
C
      ENDIF
C
      LSX=0
      CALL SORTK(ND,LSX,X,CFQ,CRRQ,LQ,MQ,NQ,KQ,NUC,NUOM)
      CALL RDIM(ND,CFQ,CRRQ,LQ,MQ,NQ,KQ,CFS,LS,MS,NS,
     &           NDV,NDW,NDT,NDH,NDG,NDVS,NDWS,NDTS,NDHS,NDGS,NDS)
C
      RETURN
      END
C
C---------------------------------------------------------------------
C-- END OF SUBROUTINE DIMEN
C---------------------------------------------------------------------
C
C
*****************************************************************
      SUBROUTINE SETQU(ND,CFI,CRRI,LI,MI,NI,KI)
*****************************************************************
C----------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CF,CFI
      CHARACTER*2 CRR,CRRI
      PARAMETER (NM=5000)
C
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),CF(NM),CRR(NM)
C
      IF( NM.NE.NMC ) THEN 
	 WRITE(*,*) 'WRONG DIMENSION NM IN SETQU (SK).'
	 STOP
      ENDIF
C
      IF( LI.LT.IABS(MI) ) RETURN
C
      ND=ND+1
C
      IF( ND.LE.NM ) THEN
         CF(ND)=CFI
         CRR(ND)=CRRI
         L(ND)=LI
         M(ND)=MI
         N(ND)=NI
         K(ND)=KI
      ENDIF
C
      RETURN
      END
C
C--------------------------------------------------------------------
C-- END OF SUBROUTINE SETQU
C--------------------------------------------------------------------
C
C
C
***********************************************************************
      FUNCTION LRIBOUN(CT,L,M,N)
***********************************************************************
C  DETERMINES WETHER (L,M,N) IS THE POSTION OF W(L=1,M=0,N=1),
C  DESCRIBING RIGID BODY ROTATION.
C----------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CT
C
      COMMON/LOG/LCALC,LWRITE,LTR,LVAR,LDY,L6,L7,L8,L9,L10
C
C      IF( CT.EQ.'W' .AND. L.EQ.1 .AND. M.EQ.0 .AND. N.EQ.1 .AND.
C     &	  LCALC.NE.1 .AND. LCALC.NE.2 ) THEN
      IF( CT.EQ.'W' .AND. L.EQ.1 .AND. N.EQ.1 ) THEN
      	 LRIBOUN=1
      ELSE
      	 LRIBOUN=0
      ENDIF
C
      RETURN
      END
C
C--------------------------------------------------------------------
C-- END OF FUNCTION LRIBOUN
C--------------------------------------------------------------------
C
C
C--------------------------------------------------------------------------
      FUNCTION NMAX(CF,FTW,FTG,L,M,M0,NT,LTR)
C--------------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CF
C
      IF( LTR.EQ.0 ) THEN
         IF( CF.EQ.'W' ) THEN      
            NMAX=IDINT(FTW*NT)
         ELSEIF( CF.EQ.'G' ) THEN      
            NMAX=IDINT(FTG*NT)
         ELSE
            NMAX=NT
         ENDIF
      ELSEIF( LTR.EQ.1 ) THEN
c----- DIESE BEDINGUNG GIBT KANN FUER (3-L+IABS(M))/2 ANDERE ERGEBNISSE 
C      LIEFERN.
        NMAX=NT-IABS(M)/M0-L/2+IABS(M)/2+1
      ELSEIF( LTR.EQ.2 ) THEN
         IF( CF.EQ.'W' ) THEN
            NMAX=IDINT( FTW*( NT-IABS(M)/M0+IABS(M)/2-L/2+1 ) )
         ELSEIF( CF.EQ.'G' ) THEN
            NMAX=IDINT( FTG*( NT-IABS(M)/M0+IABS(M)/2-L/2+1 ) )
         ELSE
            NMAX=NT-IABS(M)/M0+IABS(M)/2-L/2+1
      	 ENDIF
      ELSEIF( LTR.EQ.3 ) THEN
         IF( CF.EQ.'W' ) THEN
            NMAX=IDINT( FTW*( NT-L/M0+1 ) )
         ELSEIF( CF.EQ.'G' ) THEN
            NMAX=IDINT( FTG*( NT-L/M0+1 ) )
         ELSE
            NMAX=NT-L/M0+1
      	 ENDIF
      ENDIF
C
      RETURN
      END
C--------------------------------------------------------------------------
C
C
C--------------------------------------------------------------------------
      FUNCTION LMIN(CF,M,LEV,LCALC)
C--------------------------------------------------------------------------
      CHARACTER*1 CF
C
      IF( LEV.EQ.0 ) THEN
C----- KEINE SYMMETRIE AUSGEZEICHNET:
         IF( CF.EQ.'T' ) THEN
C-------- NUR FUER T AUCH L=0 ERLAUBT:
      	    LMIN=IABS(M)
         ELSE
            LMIN=MAX0(IABS(M),1)
         ENDIF
      ELSE
C----- DIFFERENZIERTE SYMMETRIE:
	 IF( CF.EQ.'V' ) THEN
C-------- POLOIDALES GESCHWINDIGKEITSFELD:
            IF( LCALC.EQ.1 .AND. LEV.EQ.1 ) THEN
C----------- ASYMMETRIE NUR FUER ONSET DER KONVEKTION MOEGLICH:
               LMIN=IABS(M)+1
            ELSE
C----------- SYMMETRISCHES ODER GEMISCHTES FELD:
	       IF( M.NE.0 ) THEN
	          LMIN=IABS(M)
	       ELSE
	          LMIN=2
	       ENDIF
      	    ENDIF
	 ELSEIF( CF.EQ.'W' ) THEN
C-------- TOROIDALES GESCHWINDIGKEITSFELD:
            IF( LEV.EQ.0 .OR. ( LCALC.EQ.1 .AND. LEV.EQ.1 ) ) THEN
C----------- SYMMETRISCHES ( FUER ONSET ) ODER GEMISCHTES FELD:
	       IF( M.NE.0 ) THEN
	          LMIN=IABS(M)
	       ELSEIF( M.EQ.0 ) THEN
C-------------- M=0 NICHT ERLAUBT:
	          LMIN=2
	       ENDIF
            ELSE
C----------- ASYMMETRISCHES FELD:
	       LMIN=IABS(M)+1
      	    ENDIF
	 ELSEIF( CF.EQ.'T' ) THEN     
C-------- TEMPERATUR FELD:
            IF( LCALC.EQ.1 .AND. LEV.EQ.1 ) THEN
C----------- ASYMMETRIE NUR FUER ONSET DER KONVEKTION MOEGLICH:
	       LMIN=IABS(M)+1
      	    ELSE
C----------- SYMMETRISCHES ODER GEMISCHTES FELD , L=0 MOEGLICH:
	       LMIN=IABS(M)    
      	    ENDIF
	 ELSEIF( CF.EQ.'H' ) THEN
C-------- POLOIDALES MAGNETFELD:
            IF( LEV.EQ.1 ) THEN
C----------- ASYMMETRISCHES FELD:
               LMIN=IABS(M)+1
            ELSEIF( LEV.EQ.2 ) THEN
C----------- SYMMETRISCHES FELD:
	       IF( M.NE.0 ) THEN
	          LMIN=IABS(M)
	       ELSEIF( M.EQ.0 ) THEN
C-------------- M=0 NICHT ERLAUBT:
	          LMIN=2
	       ENDIF
      	    ELSE
C----------- KEINE SYMMETRIE AUSGEZEICHNET , M=0 VERBOTEN:
      	       LMIN=MAX0(IABS(M),1)
	    ENDIF                
	 ELSEIF( CF.EQ.'G' ) THEN
C-------- TOROIDALES MAGNETFELD:
            IF( LEV.EQ.2 ) THEN
C----------- ASYMMETRISCHES FELD:
               LMIN=IABS(M)+1
            ELSEIF( LEV.EQ.1 ) THEN
C----------- SYMMETRISCHES FELD , M=0 NICHT ERLAUBT:
	       IF( M.NE.0 ) THEN
	          LMIN=IABS(M)
	       ELSE
	          LMIN=2
	       ENDIF
      	    ELSE
C----------- KEINE SYMMETRIE AUSGEZEICHNET , M=0 VERBOTEN:
      	       LMIN=MAX0(IABS(M),1)
	    ENDIF                
	 ENDIF
      ENDIF
C
      RETURN
      END
C
C--------------------------------------------------------------------------
C
C
C--------------------------------------------------------------------------
      FUNCTION LMAX(M,M0,NT,LT,LTR)
C--------------------------------------------------------------------------
C
      IF( LTR.EQ.0 ) THEN
         LMAX=2*(LT/2)+1
      ELSEIF( LTR.EQ.1 ) THEN
         LMAX=IABS(M)+2*(NT-IABS(M)/M0)+1
      ELSEIF( LTR.EQ.2 ) THEN
         LMAX=IABS(M)+2*(LT/M0-IABS(M)/M0)+1
      ELSEIF( LTR.EQ.3 ) THEN
         LMAX=LT+1
      ENDIF
C
      RETURN
      END
C
C--------------------------------------------------------------------------
C
C
C--------------------------------------------------------------------------
      FUNCTION MMIN(CF,LRB,M0,MF,NT,LT,LTR,LCALC)
C--------------------------------------------------------------------------
      CHARACTER*1 CF
C
      IF( LCALC.EQ.1 ) THEN
C----- ONSET OF CONVECTION
         MMIN=M0
      ELSE
         IF( MF.EQ.0 ) THEN
	    IF( ( CF.EQ.'H' .OR. CF.EQ.'G' ) .AND.
     &		   LRB.EQ.1 .AND. MOD(M0,2).EQ.0 ) THEN
C----------- SUBHARMONISCHES MAGNETFELD:
	       MMIN=M0/2
	    ELSE
C----------- HARMONISCHES MAGNETFELD UND ANDERE FELDER:
               MMIN=0
	    ENDIF
         ELSE
C--------- FLOQUET STOERUNG , SUBTRAKTION VON M0 ERMOEGLICHT VERGLEICH
C          ZWISCHEN (M0=2,MF=0) UND (M0=2,MF=2)
	    IF( ( CF.EQ.'H' .OR. CF.EQ.'G' ) .AND.
     &		  LRB.EQ.1 .AND. MOD(M0,2).EQ.0 ) THEN 
C----------- SUBHARMONISCHES MAGNETFELD:
               IF( LTR.EQ.0 .OR. LTR.EQ.2 .OR. LTR.EQ.3 ) THEN
      		  MMIN=-LT-M0+MF-M0/2
      	       ELSEIF( LTR.EQ.1 ) THEN
                  MMIN=-NT*M0-M0+MF-M0/2
      	       ENDIF
	    ELSE
C----------- HARMONISCHES MAGNETFELD:
      	       IF( LTR.EQ.0 .OR. LTR.EQ.2 .OR. LTR.EQ.3 ) THEN	
      		  MMIN=-LT-M0+MF
      	       ELSEIF( LTR.EQ.1 ) THEN
                  MMIN=-NT*M0-M0+MF
      	       ENDIF
	    ENDIF
         ENDIF
      ENDIF
C
      RETURN
      END
C
C--------------------------------------------------------------------------
C
C
C--------------------------------------------------------------------------
      FUNCTION MMAX(M0,MF,NT,LT,LTR,LCALC)
C--------------------------------------------------------------------------
C
      IF( LCALC.EQ.1 ) THEN
C----- ONSET OF CONVECTION
         MMAX=M0
      ELSE
         IF( LTR.EQ.0 .OR. LTR.EQ.2 .OR. LTR.EQ.3 ) THEN
            MMAX=LT+MF
      	 ELSEIF( LTR.EQ.1 ) THEN
      	    MMAX=NT*M0+MF
         ENDIF
      ENDIF
C
      RETURN
      END
C
C--------------------------------------------------------------------------
C
C
*******************************************************************
      SUBROUTINE SORTK(NK,LSX,X,CF,CRR,L,M,N,K,NUC,NUOM)
*******************************************************************
C-- SORTS THE KOEFFITIENTS IN THE REQUIRED ORDER.
C   IF MF.NE.0 THERE HAVE TO BE KOEFFITIENTS WITH M<0 , IF THERE
C   ARE NOT , THE KOEFFITIENTS FOR M>0 ARE STORED AS KOEFFITIENTS
C   M<0 WITH CC X AS WELL.
C   THE X VALUES ARE ONLY TOUCHED FOR LSX.EQ.1 .
C------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CF,CFH
      CHARACTER*2 CRR,CRRH
      DIMENSION X(*),CF(*),CRR(*),L(*),M(*),N(*),K(*)
C
      DO 2000 I=2,NK
         DO 1000 J=1,I-1
            IF( ( CF(I).EQ.'V' .AND. CF(J).NE.'V' ) 		.OR.
     &		( CF(I).EQ.'W' .AND. 
     &		  ( CF(J).NE.'V' .AND. CF(J).NE.'W' ) )		.OR.
     &	        ( CF(I).EQ.'T' .AND. 
     &		  ( CF(J).EQ.'H' .OR. CF(J).EQ.'G' ) )		.OR.
     &		( CF(I).EQ.'H' .AND. CF(J).EQ.'G' ) 		.OR.
     &      	( CF(J).EQ.CF(I) .AND. L(J).GT.L(I) ) 		.OR.
     &		( CF(J).EQ.CF(I) .AND. L(J).EQ.L(I) .AND. 
     &				       M(J).GT.M(I) )	 	.OR.
     &		( CF(J).EQ.CF(I) .AND. L(J).EQ.L(I) .AND. 
     &		   M(J).EQ.M(I) .AND. N(J).GT.N(I) ) 		.OR.
     &          ( CF(J).EQ.CF(I) .AND. L(J).EQ.L(I) .AND. 
     &		  M(J).EQ.M(I) .AND. N(J).EQ.N(I)  .AND. 
     &				     K(J).GT.K(I) ) 		.OR.
     &          ( CF(J).EQ.CF(I) .AND. L(J).EQ.L(I) .AND. 
     &		  M(J).EQ.M(I) .AND. N(J).EQ.N(I)  .AND. 
     &	             		     K(J).EQ.K(I)  .AND. 
     &		  ( CRR(I).EQ.'RR'  .OR. CRR(J).EQ.'II'  .OR.
     &		    ( CRR(J).EQ.'RI' .AND. CRR(I).EQ.'IR' ) ) ) ) THEN
               IF( LSX.EQ.1 ) XH=X(I)
      	       CFH=CF(I)	
      	       CRRH=CRR(I)
      	       LH=L(I)
      	       MH=M(I)
      	       NH=N(I)
      	       KH=K(I)
      	       IF( LSX.EQ.1 ) X(I)=X(J)
      	       CF(I)=CF(J)	
      	       CRR(I)=CRR(J)
      	       L(I)=L(J)
      	       M(I)=M(J)
      	       N(I)=N(J)
      	       K(I)=K(J)
      	       IF( LSX.EQ.1 ) X(J)=XH
      	       CF(J)=CFH
      	       CRR(J)=CRRH
      	       L(J)=LH
      	       M(J)=MH
      	       N(J)=NH
      	       K(J)=KH
	       IF( LSX.EQ.1 ) THEN
                  IF( I.EQ.NUC ) THEN
                     NUC=J
                  ELSEIF( J.EQ.NUC ) THEN
      	             NUC=I
                  ENDIF
                  IF( I.EQ.NUOM ) THEN
                     NUOM=J
                  ELSEIF( J.EQ.NUOM ) THEN
      	             NUOM=I
                  ENDIF
	       ENDIF
      	    ENDIF
1000     CONTINUE
2000  CONTINUE
C
      RETURN
      END
C
C----------------------------------------------------------------------
C
C      		
C----------------------------------------------------------------------
      SUBROUTINE GETSYM(LCALC,ND,CF,L,M,M0,LEV,LRB,LD)
C----------------------------------------------------------------------
      CHARACTER*1 CF
      DIMENSION CF(*),L(*),M(*)
C
C-- BESTIMMUNG DER SYMMETRIE DES EINGABEVEKTORS:
C   EINE SEPARATION IN E-SYMMETRISCHES ODER E-ASYMMETRISCHES MAGNETFELD
C   IST NUR MOEGLICH, WENN DAS GESCHWINDIGKEITSFELD REIN E-SYMMETRISCH IST.
C   SONST MUESSEN SOWOHL IM MAGNETFELD ALS AUCH IN V ALLE SYMMETRISCHEN WIE
C   AUCH ANTISYMMETRISCHEN KOMPONENTEN MITGENOMMEN WERDEN.
C   FUER GERADES M0 KANN DAS MAGNETFELD ROTATIONSSYMMETRISCH (LRB=2 ) ODER 
C   ROTATIONSANTISYMMETRISCH (LRB=1 ) SEIN. IST KEIN MAGNETFELD VORHANDEN,
C   SO WIRD LRB=0 GESETZ. FERNER GIBT DANN LEV DIE SYMMETRIE DES V-FELDES
C   AN MIT LEV=1 :E-ASYMMETRISCH UND LEV=2:E-SYMMETRISCH SOWIE LEV=0:
C   KEINE SYMMETRIE AUSGEZEICHNET. MIT MAGNETFELD WIRD IM FALLE LEV>0 
C   IN E-SYMMETRISCHES V-FELD VORRAUSGESETZT, LEV GIBT DANN AUSKUNFT UEBER 
C   DIE SYMMETRIE DES MAGNETFELDES: LEV=1:E-ASYMMETRISCH USW. .
      LVSYM=0
      LVASYM=0
      LHSYM=0
      LHASYM=0
      MMIN=M0
      DO 100 I=1,ND
	 IF( CF(I).EQ.'V' ) THEN
            IF( MOD(L(I)+M(I),2).EQ.0 ) THEN
               LVSYM=1
            ELSEIF( MOD(L(I)+M(I),2).EQ.1 ) THEN
               LVASYM=1
            ENDIF                             
	 ELSEIF( CF(I).EQ.'H' ) THEN
            IF( MOD(L(I)+M(I),2).EQ.0 ) THEN
               LHSYM=1
            ELSEIF( MOD(L(I)+M(I),2).EQ.1 ) THEN
               LHASYM=1
            ENDIF                             
	    IF( CF(I).EQ.'H' .AND. 
     &		 M(I).LT.MMIN .AND. M(I).GT.0 ) MMIN=M(I)
	 ENDIF
100   CONTINUE
C
      IF( LVSYM.EQ.1 .AND. LVASYM.EQ.0 ) THEN
	 IF( LHSYM.EQ.1 .AND. LHASYM.EQ.0 ) THEN
	    LEV=2
	    LD=2
	 ELSEIF( LHSYM.EQ.0 .AND. LHASYM.EQ.1 ) THEN
	    LEV=1
	    LD=2
	 ELSEIF( LHSYM.EQ.1 .AND. LHASYM.EQ.1 ) THEN
	    LEV=0
	    LD=1
	 ELSEIF( LHSYM.EQ.0 .AND. LHASYM.EQ.0 ) THEN
	    LEV=2
	    LD=2
	 ENDIF
      ELSEIF( LHSYM.EQ.0 .AND. LHASYM.EQ.1 ) THEN
	 IF( LCALC.EQ.1 ) THEN
	    LEV=1
	    LD=2
	 ELSE
	    LEV=0
	    LD=1
	 ENDIF
      ELSE
	 LEV=0
	 LD=1
      ENDIF 
      IF( MMIN.EQ.M0/2 .AND. MOD(M0,2).EQ.0 ) THEN
	 LRB=1
      ELSE
	 LRB=2
      ENDIF
      IF( LCALC.EQ.3 .OR. LCALC.EQ.5 .OR. LCALC.EQ.7 ) LRB=0
C
      RETURN
      END
C
C-----------------------------------------------------------------------
C
C
*****************************************************************
      SUBROUTINE INPUTCHECK(LCALCI,LIC,TIME0)
*****************************************************************
C
C----------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/LOG/LCALC,LWRITE,L3,LVAR,LDY,LT,L7,LREAD,L9,L10
      COMMON/PAR/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/PARI/RAI,TAI,PRI,PMI,ETAI,CI,OMI,FTWI,FTGI,MFI
      COMMON/NPAR/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD
      COMMON/NPARI/M0I,NTVI,NTHI,LTVI,LTHI,KTVI,KTHI,LEVI,LRBI,LDI
      COMMON/GAL/LL,LLC,LLM,LNL,LNLC,LNLM,LNLS,LDJ,LDK,LZ
      COMMON/CALLS/NCNONLIN,NCLIN,NCSETBL,NCINPUT
C
C-- LDIFF INDICATES WETHER A DIFFERENCE IS FOUND ORE NOT:
      LDIFFPAR=0
      LDIFFTRUNC=0
C
C-- FOR STAIBILITY ANALYSIS THE PARAMETERS OF THE INPUTFILE ARE
C   NOT INTERISTING:
      IF( ( LCALC.EQ.7 .OR. LCALC.EQ.8 ) .AND. NCINPUT.EQ.1 ) GOTO 10
C
C-- ARE THE PARAMETERS CONCERNING THE MAGNETIC FIELD
C   HAVE TO BE CHECKED?
      IF( ( LCALCI.EQ.2 .OR. LCALCI.EQ.4 .OR. LCALCI.EQ.6 ) .AND.
     &    ( LCALC.EQ.2 .OR. LCALC.EQ.4 .OR. LCALC.EQ.6 ) ) THEN
         LMAG=1
C
C--DO WE CHANGE FROM NONMAGNETIC TO MAGNETIC CASE OR VICE VERSA?
      ELSEIF( ( ( LCALCI.EQ.2 .OR. LCALCI.EQ.4 .OR. LCALCI.EQ.6 ) .AND. 
     &          ( LCALC.EQ.1 .OR. LCALC.EQ.3 .OR. LCALC.EQ.5 ) ) .OR.
     &        ( ( LCALC.EQ.2 .OR. LCALC.EQ.4 .OR. LCALC.EQ.6 ) .AND. 
     &          ( LCALCI.EQ.1 .OR. LCALCI.EQ.3 .OR. LCALCI.EQ.5 ) ) 
     &		) THEN
         LMAG=-1
      ELSE
         LMAG=0
      ENDIF
C
C-- THE TRUNCATION PARAMETERS ARE ONLY CHECKED IF LREAD.EQ.0 ,
C   OTHERWISE THEY WILL BE OVERWRITTEN BY THE TRUNCATION OF THE
C   INPUTFILE ANYWAY.
C
C-- FOR LIC=1 THE PROGRAM IS STOPED IF A DIFFERENCE IS FOUND.
C
      IF( M0.NE.M0I ) THEN
         WRITE(*,*) 'NEW GROUND MODE (OLD,NEW):.',M0,M0I
	 LDIFFTRUNC=1
      ENDIF
C
      IF( TA.NE.TAI ) THEN
         WRITE(*,*) 'NEW TAYLOR NUMBER (OLD,NEW):',TAI,TA
	 LDIFFPAR=1
      ENDIF               
      IF( RA.NE.RAI ) THEN
         WRITE(*,*) 'NEW RAYLEIGH NUMBER (OLD,NEW):',RAI,RA
	 LDIFFPAR=1
      ENDIF               
      IF( PR.NE.PRI ) THEN
         WRITE(*,*) 'NEW PRANDTL NUMBER (OLD,NEW):',PRI,PR
	 LDIFFPAR=1
      ENDIF               
      IF( PM.NE.PMI .AND. LMAG.EQ.1 ) THEN
         WRITE(*,*) 'NEW MAGNETIC PRANDTL NUMBER (OLD,NEW):',PMI,PM
	 LDIFFPAR=1
      ENDIF               
      IF( ETA.NE.ETAI ) THEN
         WRITE(*,*) 'NEW ASPECT RATIO (OLD,NEW):',ETAI,ETA
	 LDIFFPAR=1
      ENDIF                 
C
      IF( NTV.NE.NTVI .AND. LREAD.EQ.0 ) THEN
         WRITE(*,*) 'NEW RADIAL TRUNCATION OF V (OLD,NEW):',
     &							NTVI,NTV
	 LDIFFTRUNC=1
      ENDIF               
      IF( LTV.NE.LTVI .AND. LREAD.EQ.0 ) THEN
         WRITE(*,*) 'NEW SPHERICAL TRUNCATION OF V (OLD,NEW):',
     &							LTVI,LTV
	 LDIFFTRUNC=1
      ENDIF                 
      IF( NTH.NE.NTHI .AND. LMAG.EQ.1 ) THEN
         WRITE(*,*) 'NEW RADIAL TRUNCATION OF H (OLD,NEW):',
     &							NTHI,NTH
	 LDIFFTRUNC=1
      ENDIF                 
      IF( LTH.NE.LTHI .AND. LMAG.EQ.1 ) THEN
         WRITE(*,*) 'NEW SPHERICAL TRUNCATION OF H (OLD,NEW):',
     &							LTHI,LTH
	 LDIFFTRUNC=1
      ENDIF                 
      IF( KTV.NE.KTVI .AND. LREAD.EQ.0 ) THEN
         WRITE(*,*) 'NEW TIME TRUNCATION OF V (OLD,NEW):',
     &							KTVI,KTV
	 LDIFFTRUNC=1
      ENDIF                 
      IF( KTH.NE.KTHI .AND. LMAG.EQ.1 ) THEN
         WRITE(*,*) 'NEW TIME TRUNCATION OF H (OLD,NEW):',
     &							KTHI,KTH
	 LDIFFTRUNC=1
      ENDIF                 
      IF( LMAG.EQ.-1 .AND. LEV.EQ.LEVI .AND. LEV*LEVI.NE.0 ) THEN
C
C-- GOING FROM THE NON MAGNETIC TO THE MAGNETIC CASE, ORE VICE
C   VERSA, LEV CHANGES IT MEANING, SO A CHANGE IS NECESARRY IF
C   LEV.NE.0 .
         WRITE(*,*) 'WRONG EQUATORIAL SYMMETRY (OLD,NEW):',
     &                                                  LEV,LEVI
         STOP    
      ELSEIF( LMAG.NE.-1 .AND. LEV.NE.LEVI ) THEN
         WRITE(*,*) 'WRONG EQUATORIAL SYMMETRY (OLD,NEW):',
     &                                                  LEV,LEVI
         STOP    
      ENDIF
      IF( LRB.NE.LRBI .AND. LMAG.EQ.1 ) THEN
         WRITE(*,*) 'NEW SYMMETRY OF H (OLD,NEW):',LRB,LRBI
	 LDIFFTRUNC=1
      ENDIF                 
      IF( FTW.NE.FTWI .AND. LREAD.EQ.0 ) THEN
         WRITE(*,*) 'NEW TRUNCATION FACTOR FOR W (OLD,NEW):',
     &							FTWI,FTW
	 LDIFFTRUNC=1
      ENDIF                 
      IF( FTG.NE.FTGI .AND. LMAG.EQ.1 ) THEN
         WRITE(*,*) 'NEW TRUNCATION FACTOR FOR G (OLD,NEW):',
     &							FTGI,FTG
	 LDIFFTRUNC=1
      ENDIF                 
C
C-- FOR STABILITY ANALYSIS THE PARAMETERS OF THE INPUTFILE ARE TAKEN:
10    CONTINUE
      IF( LCALC.EQ.7 .OR. LCALC.EQ.8 ) THEN
	 ETA=ETAI
         M0=M0I
         NTV=NTVI
         LTV=LTVI
         KTV=KTVI
         NTH=NTHI
         LTH=LTHI
         KTH=KTHI
         LEV=LEVI
         LRB=LRBI
	 FTW=FTWI
	 FTG=FTGI
         TA=TAI
         RA=RAI	
         PR=PRI
         PM=PMI
      ENDIF
C
C-- IF THE DIMENSION CHANGES FOR A NEW INPUTFILE, ALL TERMS 
C   HAVE TO BE RECALCULATED, BECAUSE THEIR ORDER IN THE VECTOR x
C   CHANGES, FOR A CHANGE IN THE PARAMETERS ONLY LINEAR TERMS 
C   MUST BE CALCULATED NEW:
      IF( NCINPUT.GT.1 ) THEN
	 IF( LDIFFTRUNC.EQ.1 ) THEN
            LL=1
            LNL=1
            LNLS=1
	 ENDIF
	 IF( LDIFFPAR.EQ.1 ) THEN
            LL=1
	 ENDIF
      ENDIF
      IF( LDIFFPAR.EQ.1 .OR. LDIFFTRUNC.EQ.1 ) THEN
	 IF( LIC.EQ.1 ) STOP
	 WRITE(*,*) 'NEW DATASET'
	 TIME0=0.D0                            
      ENDIF
C
      RETURN
      END
C
C---------------------------------------------------------------------------- 
C
