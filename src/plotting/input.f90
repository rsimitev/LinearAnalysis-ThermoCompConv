   subroutine READH(unit_out,LTR,NUDSR,TIMER,NUDS,TIME,LDR)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*2 C2
      PARAMETER(NMT=50)

      COMMON/LOG/LCALC,LWRITE,L3,LVAR,LDY,LT,L7,L8,L9,L10
      COMMON/PARI/RAI,TAI,PRI,PMI,ETAI,CI,OMI,FTWI,FTGI,MFI
      COMMON/NPARI/M0I,NTVI,NTHI,LTVI,LTHI,KTVI,KTHI,LEVI,LRBI,LDI
      COMMON/TINT/DT,TIME0,NST,NTW,NDTW,NDTPW,NKWT,IKWT(NMT)
      COMMON/QNUVI/NUCI,NUOMI

      READ(unit_out,*) M0I,NTVI,NTHI,KTVI,KTHI,LTVI,LTHI,FTWI,FTGI
      READ(unit_out,*) NUDS,TAI,RAI,PRI,PMI,ETAI
      READ(unit_out,*) CI,OMI,NUCI,NUOMI,MFI
      TAI=DSQRT(TAI)

      IF( LTR.EQ.1 ) THEN
         READ(unit_out,'(1X,A2,I6,D16.6)') C2,NUDS,TIME
      ENDIF
      IF( ( LT.EQ.1 .AND. TIME.EQ.TIMER ) .OR.
     &    ( LT.EQ.0 .AND. NUDS.EQ.NUDSR )  ) THEN
         LDR=1
         IF( LTR.EQ.1 ) THEN
            TIME0=TIME
         ELSE
            TIME0=0.D0
         ENDIF
      ELSE
         LDR=0
      ENDIF
   END subroutine

   !----------------------------------------------------------------
   subroutine READD(unit_out,LDR,NKR,XR,CFR,CRRR,LR,MR,NR,KR,
     &                       EVPM,EVPF,EVTM,DNU,EVTF,EMPM,EMPF,EMTM,EMTF)


      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CFI,CFR
      CHARACTER*2 C2,CRRR

      DIMENSION XR(*),LR(*),MR(*),NR(*),KR(*),CFR(*),CRRR(*)

      COMMON/LOG/LCALC,LWRITE,L3,LVAR,LDY,LT,L7,LREAD,L9,L10
      COMMON/NPARI/M0I,NTVI,NTHI,LTVI,LTHI,KTVI,KTHI,LEVI,LRBI,LDI

      NKR=0
      DO I=1, 5000
         READ(unit_out,'(1X,A2)',END=1100,ERR=1100) C2
         select case (C2)
            case('V ','W ','T ','G ','H ')
               IF( LDR.EQ.1 ) THEN
                  BACKSPACE(unit_out)
                  READ(unit_out,'(1X,A1,4I3,4D16.8)') CFI,LI,MI,NI,KI,DRR,DIR,DRI,DII
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
            case('EV')
               IF( LDR.EQ.1 ) THEN
                  BACKSPACE(unit_out)
                  READ(unit_out,'(1X,A2,5D14.6)') C2,EVPM,EVPF,EVTM,EVTF,DNU
               ENDIF
            case('EM')
               IF( LDR.EQ.1 ) THEN
                  BACKSPACE(unit_out)
                  READ(unit_out,'(1X,A2,4D14.6)') C2,EMPM,EMPF,EMTM,EMTF
               ENDIF
            case default
               BACKSPACE(unit_out)
               exit 
          end select
       enddo
1100  CONTINUE
   END subroutine

   !-------------------------------------------------------------------------
   subroutine READP(unit_out,LST,LTR,NUDSR,TIMER,NUDS,TIME,LDR,LEND)

      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*2 C2
      PARAMETER(NMT=50)

      COMMON/LOG/LCALC,LWRITE,L3,LVAR,LDY,LT,L7,L8,L9,L10
      COMMON/PARI/RAI,TAI,PRI,PMI,ETAI,CI,OMI,FTWI,FTGI,MFI
      COMMON/NPARI/M0I,NTVI,NTHI,LTVI,LTHI,KTVI,KTHI,LEVI,LRBI,LDI
      COMMON/TINT/DT,TIME0,NST,NTW,NDTW,NDTPW,NKWT,IKWT(NMT)
      COMMON/QNUVI/NUCI,NUOMI

      IF( LTR.EQ.1 ) THEN
         READ(unit_out,'(1X,A2,I6,D16.6)',END=2000) C2,NUDS,TIME
      ELSEIF( LTR.EQ.0 ) THEN
         IF( LST.EQ.2 ) READ(unit_out,*,END=2000) M0I,NTVI,NTHI,KTVI,KTHI,LTVI,LTHI,FTWI,FTGI
         READ(unit_out,*,END=2000) NUDS,TAI,RAI,PRI,PMI,ETAI
         READ(unit_out,*) CI,OMI,NUCI,NUOMI
         TAI=DSQRT(TAI)
      ENDIF
      IF( ( LT.EQ.1 .AND. TIME.EQ.TIMER ) .OR.
     &    ( LT.EQ.0 .AND. NUDS.EQ.NUDSR )  ) THEN
         LDR=1
         IF( LTR.EQ.1 ) THEN
            TIME0=TIME
         ELSE
            TIME0=0.D0
         ENDIF
      ELSE
         LDR=0
      ENDIF

      LEND=0
      RETURN

2000  CONTINUE
      LEND=1

   END


   !----------------------------------------------------------------
   subroutine RDIM(ND,CF,CRR,L,M,N,K,CFS,LS,MS,NS,
     &      NDV,NDW,NDT,NDH,NDG,NDVS,NDWS,NDTS,NDHS,NDGS,NDS)

      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CF,CFS
      CHARACTER*2 CRR

      DIMENSION CF(*),CRR(*),L(*),M(*),N(*),K(*)
      DIMENSION CFS(*),LS(*),MS(*),NS(*)

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

      DO I=1, ND
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
      enddo
   END subroutine

   !-------------------------------------------------------------------
   subroutine KOEFF(LOU)

      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CF,CFS,CFI,CFSI
      CHARACTER*2 CRR,CRRI
      PARAMETER(NM=5000,NMT=50)

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

      IF( NMC.NE.NM ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN KOEFF.'
         STOP
      ENDIF

      IF( LOU.EQ.0 ) THEN

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

         OPEN(21,FILE='koeff.dat',STATUS='unknown')

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

!-- AUSSCHREIBEN DER BENUTZTEN KOEFFITIENTEN
         WRITE(21,*) 'BENUTZTE KOEFFITIENTEN FUER RECHNUNG:'
         DO 100 I=1,ND
100      WRITE(21,'(I5,A2,A3,4I4)') I,CF(I),CRR(I),L(I),M(I),N(I),K(I)
         WRITE(21,*)
         WRITE(21,*) 'REDUZIERTE KOEFFITIENTEN:'
         DO 200 I=1,NDS
200      WRITE(21,'(I5,A2,3I4)') I,CFS(I),LS(I),MS(I),NS(I)

!-- SETZEN DER NUMMERN DER AUSZUSCHREIBENDEN KOEFFITIENTEN:
         IF( LCALC.EQ.5 .OR. LCALC.EQ.6 ) THEN
            CALL SETWT(ND)
            WRITE(21,*)
            WRITE(21,*) 'NUMMERN DER AUSGESCHRIEBENEN KOEFFITIENTEN:'
            DO 300 I=1,NKWT
 300        WRITE(21,*) I,IKWT(I)
         ENDIF
         CLOSE(21)

      ELSEIF( LOU.EQ.1 ) THEN

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

         OPEN(21,FILE='koeffi.dat',STATUS='unknown')

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

!-- AUSSCHREIBEN DER BENUTZTEN KOEFFITIENTEN
         WRITE(21,*) 'BENUTZTE KOEFFITIENTEN DES INPUTS:'
         DO 400 I=1,NDI
400      WRITE(21,'(I5,A2,A3,4I4)')
     &                        I,CFI(I),CRRI(I),LI(I),MI(I),NI(I),KI(I)
         WRITE(21,*)
         WRITE(21,*) 'REDUZIERTE KOEFFITIENTEN DES INPUTS:'
         DO 500 I=1,NDSI
500      WRITE(21,'(I5,A2,3I4)') I,CFSI(I),LSI(I),MSI(I),NSI(I)

      ENDIF

      RETURN
      END

   !---------------------------------------------------------------------
   subroutine SETWT(ND)

      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CF,CFWT
      CHARACTER*2 CRR,CRRWT
      PARAMETER (NM=5000,NMT=50)

      COMMON/QNU/NMQC,L(NM),M(NM),N(NM),K(NM),CF(NM),CRR(NM)
      COMMON/TINT/DT,TIME0,NST,NTW,NDTW,NDTPW,NKWT,IKWT(NMT)
      COMMON/COEW/NMKTC,NKT,LWT(NMT),MWT(NMT),NWT(NMT),KWT(NMT),
     &                  CFWT(NMT),CRRWT(NMT)

      IF( NMT.NE.NMKTC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NMT IN SETWT (SK).'
         STOP
      ENDIF

      DO 200 J=1,NKT
      DO 100 I=1,ND
         IF( CF(I).EQ.CFWT(J) .AND. CRR(I).EQ.CRRWT(J) .AND.
     &             L(I).EQ.LWT(J) .AND. M(I).EQ.MWT(J) .AND.
     &                        N(I).EQ.NWT(J) .AND. K(I).EQ.KWT(J) )  THEN
                  NKWT=NKWT+1
                  IKWT(NKWT)=I
            GOTO 200
         ENDIF
100   CONTINUE
200   CONTINUE

      RETURN
      END

   !--------------------------------------------------------------------
   subroutine DIMEN

      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CF,CFQ,CFS
      CHARACTER*2 CRR,CRRQ
      PARAMETER (NM=5000)

      DIMENSION X(NM)

      COMMON/LOG/LCALC,LWRITE,LTR,LVAR,LDY,L6,L7,L8,L9,L10
      COMMON/NUM/RELE,EPS,ALPH,STOER,NITMAX,NJA
      COMMON/PAR/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPAR/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LDI
      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/DIMS/NDVS,NDWS,NDTS,NDHS,NDGS,NDS
      COMMON/QNUV/NUC,NUOM

      COMMON/QNU/NMQC,LQ(NM),MQ(NM),NQ(NM),KQ(NM),CFQ(NM),CRRQ(NM)
      COMMON/QNUS/NMSC,LS(NM),MS(NM),NS(NM),CFS(NM)

      IF( NM.NE.NMQC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN DIMEN (INPUT).'
         STOP
      ENDIF
      IF( NM.NE.NMSC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN DIMEN (INPUT).'
         STOP
      ENDIF

      ND=0

      IF( LCALC.NE.2 ) THEN
         CF='V'
         DO 100 M=MMIN(CF,LRB,M0,MF,NTV,LTV,LTR,LCALC),
     &                        MMAX(M0,MF,NTV,LTV,LTR,LCALC),M0
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

         IF( LCALC.EQ.1 .AND. TA.LT.EPS ) GOTO 250
   !-------  ONSET IN NON ROTATING CASE HAS NO TOROIDAL COMPONENT.

         CF='W'
         DO 200 M=MMIN(CF,LRB,M0,MF,NTV,LTV,LTR,LCALC),
     &                        MMAX(M0,MF,NTV,LTV,LTR,LCALC),M0
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

250      CONTINUE

         CF='T'
         DO 300 M=MMIN(CF,LRB,M0,MF,NTV,LTV,LTR,LCALC),
     &                        MMAX(M0,MF,NTV,LTV,LTR,LCALC),M0
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

      IF( LCALC.EQ.2 .OR. LCALC.EQ.4 .OR. LCALC.EQ.6 .OR.
     &                                          LCALC.EQ.8 ) THEN
         CF='H'
         DO 400 M=MMIN(CF,LRB,M0,MF,NTH,LTV,LTR,LCALC),
     &                        MMAX(M0,MF,NTH,LTH,LTR,LCALC),M0
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

         CF='G'
         DO 500 M=MMIN(CF,LRB,M0,MF,NTH,LTV,LTR,LCALC),
     &                        MMAX(M0,MF,NTH,LTH,LTR,LCALC),M0
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

      ENDIF

      LSX=0
      CALL SORTK(ND,LSX,X,CFQ,CRRQ,LQ,MQ,NQ,KQ,NUC,NUOM)
      CALL RDIM(ND,CFQ,CRRQ,LQ,MQ,NQ,KQ,CFS,LS,MS,NS,
     &           NDV,NDW,NDT,NDH,NDG,NDVS,NDWS,NDTS,NDHS,NDGS,NDS)

      RETURN
      END

   !----------------------------------------------------------------
   subroutine SETQU(ND,CFI,CRRI,LI,MI,NI,KI)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CF,CFI
      CHARACTER*2 CRR,CRRI
      PARAMETER (NM=5000)

      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),CF(NM),CRR(NM)

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN SETQU (SK).'
         STOP
      ENDIF

      IF( LI.LT.IABS(MI) ) RETURN

      ND=ND+1

      IF( ND.LE.NM ) THEN
         CF(ND)=CFI
         CRR(ND)=CRRI
         L(ND)=LI
         M(ND)=MI
         N(ND)=NI
         K(ND)=KI
      ENDIF

      RETURN
   END

   !--------------------------------------------------------------------
!  DETERMINES WETHER (L,M,N) IS THE POSTION OF W(L=1,M=0,N=1),
!  DESCRIBING RIGID BODY ROTATION.
   FUNCTION LRIBOUN(CT,L,M,N)

      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CT

      COMMON/LOG/LCALC,LWRITE,LTR,LVAR,LDY,L6,L7,L8,L9,L10

!      IF( CT.EQ.'W' .AND. L.EQ.1 .AND. M.EQ.0 .AND. N.EQ.1 .AND.
!     &          LCALC.NE.1 .AND. LCALC.NE.2 ) THEN
      IF( CT.EQ.'W' .AND. L.EQ.1 .AND. N.EQ.1 ) THEN
               LRIBOUN=1
      ELSE
               LRIBOUN=0
      ENDIF
   END function

   !--------------------------------------------------------------------------
   FUNCTION NMAX(CF,FTW,FTG,L,M,M0,NT,LTR)
   !--------------------------------------------------------------------------

      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CF

      IF( LTR.EQ.0 ) THEN
         IF( CF.EQ.'W' ) THEN
            NMAX=IDINT(FTW*NT)
         ELSEIF( CF.EQ.'G' ) THEN
            NMAX=IDINT(FTG*NT)
         ELSE
            NMAX=NT
         ENDIF
      ELSEIF( LTR.EQ.1 ) THEN
   !----- DIESE BEDINGUNG GIBT KANN FUER (3-L+IABS(M))/2 ANDERE ERGEBNISSE
!      LIEFERN.
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

      RETURN
      END
   !--------------------------------------------------------------------------


   !--------------------------------------------------------------------------
      FUNCTION LMIN(CF,M,LEV,LCALC)
   !--------------------------------------------------------------------------
      CHARACTER*1 CF

      IF( LEV.EQ.0 ) THEN
   !----- KEINE SYMMETRIE AUSGEZEICHNET:
         IF( CF.EQ.'T' ) THEN
   !-------- NUR FUER T AUCH L=0 ERLAUBT:
                  LMIN=IABS(M)
         ELSE
            LMIN=MAX0(IABS(M),1)
         ENDIF
      ELSE
   !----- DIFFERENZIERTE SYMMETRIE:
         IF( CF.EQ.'V' ) THEN
   !-------- POLOIDALES GESCHWINDIGKEITSFELD:
            IF( LCALC.EQ.1 .AND. LEV.EQ.1 ) THEN
   !----------- ASYMMETRIE NUR FUER ONSET DER KONVEKTION MOEGLICH:
               LMIN=IABS(M)+1
            ELSE
   !----------- SYMMETRISCHES ODER GEMISCHTES FELD:
               IF( M.NE.0 ) THEN
                  LMIN=IABS(M)
               ELSE
                  LMIN=2
               ENDIF
                  ENDIF
         ELSEIF( CF.EQ.'W' ) THEN
   !-------- TOROIDALES GESCHWINDIGKEITSFELD:
            IF( LEV.EQ.0 .OR. ( LCALC.EQ.1 .AND. LEV.EQ.1 ) ) THEN
   !----------- SYMMETRISCHES ( FUER ONSET ) ODER GEMISCHTES FELD:
               IF( M.NE.0 ) THEN
                  LMIN=IABS(M)
               ELSEIF( M.EQ.0 ) THEN
   !-------------- M=0 NICHT ERLAUBT:
                  LMIN=2
               ENDIF
            ELSE
   !----------- ASYMMETRISCHES FELD:
               LMIN=IABS(M)+1
                  ENDIF
         ELSEIF( CF.EQ.'T' ) THEN
   !-------- TEMPERATUR FELD:
            IF( LCALC.EQ.1 .AND. LEV.EQ.1 ) THEN
   !----------- ASYMMETRIE NUR FUER ONSET DER KONVEKTION MOEGLICH:
               LMIN=IABS(M)+1
                  ELSE
   !----------- SYMMETRISCHES ODER GEMISCHTES FELD , L=0 MOEGLICH:
               LMIN=IABS(M)
                  ENDIF
         ELSEIF( CF.EQ.'H' ) THEN
   !-------- POLOIDALES MAGNETFELD:
            IF( LEV.EQ.1 ) THEN
   !----------- ASYMMETRISCHES FELD:
               LMIN=IABS(M)+1
            ELSEIF( LEV.EQ.2 ) THEN
   !----------- SYMMETRISCHES FELD:
               IF( M.NE.0 ) THEN
                  LMIN=IABS(M)
               ELSEIF( M.EQ.0 ) THEN
   !-------------- M=0 NICHT ERLAUBT:
                  LMIN=2
               ENDIF
                  ELSE
   !----------- KEINE SYMMETRIE AUSGEZEICHNET , M=0 VERBOTEN:
                     LMIN=MAX0(IABS(M),1)
            ENDIF
         ELSEIF( CF.EQ.'G' ) THEN
   !-------- TOROIDALES MAGNETFELD:
            IF( LEV.EQ.2 ) THEN
   !----------- ASYMMETRISCHES FELD:
               LMIN=IABS(M)+1
            ELSEIF( LEV.EQ.1 ) THEN
   !----------- SYMMETRISCHES FELD , M=0 NICHT ERLAUBT:
               IF( M.NE.0 ) THEN
                  LMIN=IABS(M)
               ELSE
                  LMIN=2
               ENDIF
                  ELSE
   !----------- KEINE SYMMETRIE AUSGEZEICHNET , M=0 VERBOTEN:
                     LMIN=MAX0(IABS(M),1)
            ENDIF
         ENDIF
      ENDIF

      RETURN
   END function

   !--------------------------------------------------------------------------
   FUNCTION LMAX(M,M0,NT,LT,LTR)
      IF( LTR.EQ.0 ) THEN
         LMAX=2*(LT/2)+1
      ELSEIF( LTR.EQ.1 ) THEN
         LMAX=IABS(M)+2*(NT-IABS(M)/M0)+1
      ELSEIF( LTR.EQ.2 ) THEN
         LMAX=IABS(M)+2*(LT/M0-IABS(M)/M0)+1
      ELSEIF( LTR.EQ.3 ) THEN
         LMAX=LT+1
      ENDIF

      RETURN
   END function

   !--------------------------------------------------------------------------
   FUNCTION MMIN(CF,LRB,M0,MF,NT,LT,LTR,LCALC)
      CHARACTER*1 CF

      IF( LCALC.EQ.1 ) THEN
   !----- ONSET OF CONVECTION
         MMIN=M0
      ELSE
         IF( MF.EQ.0 ) THEN
            IF( ( CF.EQ.'H' .OR. CF.EQ.'G' ) .AND.
     &                   LRB.EQ.1 .AND. MOD(M0,2).EQ.0 ) THEN
   !----------- SUBHARMONISCHES MAGNETFELD:
               MMIN=M0/2
            ELSE
   !----------- HARMONISCHES MAGNETFELD UND ANDERE FELDER:
               MMIN=0
            ENDIF
         ELSE
   !--------- FLOQUET STOERUNG , SUBTRAKTION VON M0 ERMOEGLICHT VERGLEICH
!          ZWISCHEN (M0=2,MF=0) UND (M0=2,MF=2)
            IF( ( CF.EQ.'H' .OR. CF.EQ.'G' ) .AND.
     &                  LRB.EQ.1 .AND. MOD(M0,2).EQ.0 ) THEN
   !----------- SUBHARMONISCHES MAGNETFELD:
               IF( LTR.EQ.0 .OR. LTR.EQ.2 .OR. LTR.EQ.3 ) THEN
                        MMIN=-LT-M0+MF-M0/2
                     ELSEIF( LTR.EQ.1 ) THEN
                  MMIN=-NT*M0-M0+MF-M0/2
                     ENDIF
            ELSE
   !----------- HARMONISCHES MAGNETFELD:
                     IF( LTR.EQ.0 .OR. LTR.EQ.2 .OR. LTR.EQ.3 ) THEN
                        MMIN=-LT-M0+MF
                     ELSEIF( LTR.EQ.1 ) THEN
                  MMIN=-NT*M0-M0+MF
                     ENDIF
            ENDIF
         ENDIF
      ENDIF
   END function

   !--------------------------------------------------------------------------
   FUNCTION MMAX(M0,MF,NT,LT,LTR,LCALC)
      IF( LCALC.EQ.1 ) THEN
   !----- ONSET OF CONVECTION
         MMAX=M0
      ELSE
         IF( LTR.EQ.0 .OR. LTR.EQ.2 .OR. LTR.EQ.3 ) THEN
            MMAX=LT+MF
               ELSEIF( LTR.EQ.1 ) THEN
                  MMAX=NT*M0+MF
         ENDIF
      ENDIF
   END function

   !--------------------------------------------------------------------------
!-- SORTS THE KOEFFITIENTS IN THE REQUIRED ORDER.
!   IF MF.NE.0 THERE HAVE TO BE KOEFFITIENTS WITH M<0 , IF THERE
!   ARE NOT , THE KOEFFITIENTS FOR M>0 ARE STORED AS KOEFFITIENTS
!   M<0 WITH CC X AS WELL.
!   THE X VALUES ARE ONLY TOUCHED FOR LSX.EQ.1 .
   subroutine SORTK(NK,LSX,X,CF,CRR,L,M,N,K,NUC,NUOM)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CF,CFH
      CHARACTER*2 CRR,CRRH
      DIMENSION X(*),CF(*),CRR(*),L(*),M(*),N(*),K(*)

      DO 2000 I=2,NK
         DO 1000 J=1,I-1
            IF( ( CF(I).EQ.'V' .AND. CF(J).NE.'V' )                 .OR.
     &                ( CF(I).EQ.'W' .AND.
     &                  ( CF(J).NE.'V' .AND. CF(J).NE.'W' ) )                .OR.
     &                ( CF(I).EQ.'T' .AND.
     &                  ( CF(J).EQ.'H' .OR. CF(J).EQ.'G' ) )                .OR.
     &                ( CF(I).EQ.'H' .AND. CF(J).EQ.'G' )                 .OR.
     &              ( CF(J).EQ.CF(I) .AND. L(J).GT.L(I) )                 .OR.
     &                ( CF(J).EQ.CF(I) .AND. L(J).EQ.L(I) .AND.
     &                                       M(J).GT.M(I) )                 .OR.
     &                ( CF(J).EQ.CF(I) .AND. L(J).EQ.L(I) .AND.
     &                   M(J).EQ.M(I) .AND. N(J).GT.N(I) )                 .OR.
     &          ( CF(J).EQ.CF(I) .AND. L(J).EQ.L(I) .AND.
     &                  M(J).EQ.M(I) .AND. N(J).EQ.N(I)  .AND.
     &                                     K(J).GT.K(I) )                 .OR.
     &          ( CF(J).EQ.CF(I) .AND. L(J).EQ.L(I) .AND.
     &                  M(J).EQ.M(I) .AND. N(J).EQ.N(I)  .AND.
     &                                          K(J).EQ.K(I)  .AND.
     &                  ( CRR(I).EQ.'RR'  .OR. CRR(J).EQ.'II'  .OR.
     &                    ( CRR(J).EQ.'RI' .AND. CRR(I).EQ.'IR' ) ) ) ) THEN
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

      RETURN
   END subroutine

   !----------------------------------------------------------------------


   !----------------------------------------------------------------------
   subroutine GETSYM(LCALC,ND,CF,L,M,M0,LEV,LRB,LD)
   !----------------------------------------------------------------------
      CHARACTER*1 CF
      DIMENSION CF(*),L(*),M(*)

!-- BESTIMMUNG DER SYMMETRIE DES EINGABEVEKTORS:
!   EINE SEPARATION IN E-SYMMETRISCHES ODER E-ASYMMETRISCHES MAGNETFELD
!   IST NUR MOEGLICH, WENN DAS GESCHWINDIGKEITSFELD REIN E-SYMMETRISCH IST.
!   SONST MUESSEN SOWOHL IM MAGNETFELD ALS AUCH IN V ALLE SYMMETRISCHEN WIE
!   AUCH ANTISYMMETRISCHEN KOMPONENTEN MITGENOMMEN WERDEN.
!   FUER GERADES M0 KANN DAS MAGNETFELD ROTATIONSSYMMETRISCH (LRB=2 ) ODER
!   ROTATIONSANTISYMMETRISCH (LRB=1 ) SEIN. IST KEIN MAGNETFELD VORHANDEN,
!   SO WIRD LRB=0 GESETZ. FERNER GIBT DANN LEV DIE SYMMETRIE DES V-FELDES
!   AN MIT LEV=1 :E-ASYMMETRISCH UND LEV=2:E-SYMMETRISCH SOWIE LEV=0:
!   KEINE SYMMETRIE AUSGEZEICHNET. MIT MAGNETFELD WIRD IM FALLE LEV>0
!   IN E-SYMMETRISCHES V-FELD VORRAUSGESETZT, LEV GIBT DANN AUSKUNFT UEBER
!   DIE SYMMETRIE DES MAGNETFELDES: LEV=1:E-ASYMMETRISCH USW. .
      LVSYM=0
      LVASYM=0
      LHSYM=0
      LHASYM=0
      MMIN=M0
      DO I=1,ND
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
            IF( CF(I).EQ.'H' .AND. M(I).LT.MMIN .AND. M(I).GT.0 ) MMIN=M(I)
         ENDIF
      enddo

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

      RETURN
      END

   !-----------------------------------------------------------------------
   subroutine INPUTCHECK(LCALCI,LIC,TIME0)
      IMPLICIT REAL*8(A-H,O-Z)

      COMMON/LOG/LCALC,LWRITE,L3,LVAR,LDY,LT,L7,LREAD,L9,L10
      COMMON/PAR/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/PARI/RAI,TAI,PRI,PMI,ETAI,CI,OMI,FTWI,FTGI,MFI
      COMMON/NPAR/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD
      COMMON/NPARI/M0I,NTVI,NTHI,LTVI,LTHI,KTVI,KTHI,LEVI,LRBI,LDI
      COMMON/GAL/LL,LLC,LLM,LNL,LNLC,LNLM,LNLS,LDJ,LDK,LZ
      COMMON/CALLS/NCNONLIN,NCLIN,NCSETBL,NCINPUT

!-- LDIFF INDICATES WETHER A DIFFERENCE IS FOUND ORE NOT:
      LDIFFPAR=0
      LDIFFTRUNC=0

!-- FOR STAIBILITY ANALYSIS THE PARAMETERS OF THE INPUTFILE ARE
!   NOT INTERISTING:
      IF( ( LCALC.EQ.7 .OR. LCALC.EQ.8 ) .AND. NCINPUT.EQ.1 ) GOTO 10

!-- ARE THE PARAMETERS CONCERNING THE MAGNETIC FIELD
!   HAVE TO BE CHECKED?
      IF( ( LCALCI.EQ.2 .OR. LCALCI.EQ.4 .OR. LCALCI.EQ.6 ) .AND.
     &    ( LCALC.EQ.2  .OR. LCALC.EQ.4  .OR. LCALC.EQ.6  ) ) THEN
         LMAG=1

!--DO WE CHANGE FROM NONMAGNETIC TO MAGNETIC CASE OR VICE VERSA?
      ELSEIF( ( ( LCALCI.EQ.2 .OR. LCALCI.EQ.4 .OR. LCALCI.EQ.6 ) .AND.
     &          ( LCALC.EQ.1  .OR. LCALC.EQ.3  .OR. LCALC.EQ.5  ) ) .OR.
     &        ( ( LCALC.EQ.2  .OR. LCALC.EQ.4  .OR. LCALC.EQ.6  ) .AND.
     &          ( LCALCI.EQ.1 .OR. LCALCI.EQ.3 .OR. LCALCI.EQ.5 ) )
     &                ) THEN
         LMAG=-1
      ELSE
         LMAG=0
      ENDIF

!-- THE TRUNCATION PARAMETERS ARE ONLY CHECKED IF LREAD.EQ.0 ,
!   OTHERWISE THEY WILL BE OVERWRITTEN BY THE TRUNCATION OF THE
!   INPUTFILE ANYWAY.

!-- FOR LIC=1 THE PROGRAM IS STOPED IF A DIFFERENCE IS FOUND.

      IF( M0.NE.M0I ) THEN
         WRITE(*,*) 'NEW GROUND MODE (OLD,NEW):.',M0,M0I
         LDIFFTRUNC=1
      ENDIF

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

      IF( NTV.NE.NTVI .AND. LREAD.EQ.0 ) THEN
         WRITE(*,*) 'NEW RADIAL TRUNCATION OF V (OLD,NEW):', NTVI,NTV
         LDIFFTRUNC=1
      ENDIF
      IF( LTV.NE.LTVI .AND. LREAD.EQ.0 ) THEN
         WRITE(*,*) 'NEW SPHERICAL TRUNCATION OF V (OLD,NEW):', LTVI,LTV
         LDIFFTRUNC=1
      ENDIF
      IF( NTH.NE.NTHI .AND. LMAG.EQ.1 ) THEN
         WRITE(*,*) 'NEW RADIAL TRUNCATION OF H (OLD,NEW):', NTHI,NTH
         LDIFFTRUNC=1
      ENDIF
      IF( LTH.NE.LTHI .AND. LMAG.EQ.1 ) THEN
         WRITE(*,*) 'NEW SPHERICAL TRUNCATION OF H (OLD,NEW):', LTHI,LTH
         LDIFFTRUNC=1
      ENDIF
      IF( KTV.NE.KTVI .AND. LREAD.EQ.0 ) THEN
         WRITE(*,*) 'NEW TIME TRUNCATION OF V (OLD,NEW):', KTVI,KTV
         LDIFFTRUNC=1
      ENDIF
      IF( KTH.NE.KTHI .AND. LMAG.EQ.1 ) THEN
         WRITE(*,*) 'NEW TIME TRUNCATION OF H (OLD,NEW):', KTHI,KTH
         LDIFFTRUNC=1
      ENDIF
      IF( LMAG.EQ.-1 .AND. LEV.EQ.LEVI .AND. LEV*LEVI.NE.0 ) THEN

!-- GOING FROM THE NON MAGNETIC TO THE MAGNETIC CASE, ORE VICE
!   VERSA, LEV CHANGES IT MEANING, SO A CHANGE IS NECESARRY IF
!   LEV.NE.0 .
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
     &                                                        FTWI,FTW
         LDIFFTRUNC=1
      ENDIF
      IF( FTG.NE.FTGI .AND. LMAG.EQ.1 ) THEN
         WRITE(*,*) 'NEW TRUNCATION FACTOR FOR G (OLD,NEW):',
     &                                                        FTGI,FTG
         LDIFFTRUNC=1
      ENDIF

!-- FOR STABILITY ANALYSIS THE PARAMETERS OF THE INPUTFILE ARE TAKEN:
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

!-- IF THE DIMENSION CHANGES FOR A NEW INPUTFILE, ALL TERMS
!   HAVE TO BE RECALCULATED, BECAUSE THEIR ORDER IN THE VECTOR x
!   CHANGES, FOR A CHANGE IN THE PARAMETERS ONLY LINEAR TERMS
!   MUST BE CALCULATED NEW:
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

      RETURN
   END

! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
