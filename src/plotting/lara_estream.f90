!************************************************************************
!-- PROGRAM TO read data from Galerkinprogram of J.W. and convert it for IDL.
!--
!-- Input:   stdin  (short version of LARA.EXE for DISSPLA)
!--
!-- Output:
!--          3 files: idl.z, idl.x, idl.y
!--
!------------------------------------------------------------------------
PROGRAM LARA
   implicit none
   double precision, parameter:: PI = 3.14159265358979D0
   integer, PARAMETER:: NM = 5500, NAM = 400, nPlotsMAX = 9, nSubPlotsMAX = 4, NPAM = nPlotsMAX*nSubPlotsMAX
   integer, PARAMETER:: NLMA = 100
   integer, PARAMETER:: NMX = 65, NMY = 128
   integer:: NMC
   integer:: NCPLOT, LR, NQ, NR, n, l
   integer:: drawPlotNum, drawFrame, drawHeader, drawTime
   logical:: countourParIsNumber
   integer:: dataSetNumber
   CHARACTER*40 INPUTFILE,OUTPUTFILE
   CHARACTER*30 CTEXT1
   CHARACTER*10 CTEXT2
   CHARACTER*1 whatToPlot,CFS
   CHARACTER*2 CRR
   integer:: timeSeriesControl
   INTEGER:: i, j
   double precision:: THETA(NMC)
   double precision:: dt
   double precision:: headerSpaceY

   double precision:: DX(NM), constantCoordinateValue(nPlotsMAX,nSubPlotsMAX),normRadiusMax(nPlotsMAX,nSubPlotsMAX)
   double precision:: TIME(nPlotsMAX), contourPar(nPlotsMAX,nSubPlotsMAX)
   integer:: nSubPlots(nPlotsMAX)
   character(len=1):: constantCoordinate(nPlotsMAX,nSubPlotsMAX), thisPlotconstantCoordinate(NPAM)
   character(len=2):: whatToPlot(nPlotsMAX,nSubPlotsMAX), thisPlotWhatToPlot(NPAM)
   character(len=2):: domain(nPlotsMAX), quadrant(nPlotsMAX,nSubPlotsMAX), thisPlotQuadrant(NPAM), XCP(NPAM)
   character(len=3):: subPlotLabel(nPlotsMAX,nSubPlotsMAX), ABCN(NPAM)
   double precision:: zdo
   DIMENSION ZDP(NPAM),TIMEP(NPAM)
   DIMENSION XOR(NPAM),YOR(NPAM),XAR(NPAM),YAR(NPAM)
   DIMENSION XROCM(NPAM),XRICM(NPAM),XRMCM(NPAM),XRM(NPAM)

   COMMON/LOG/LCALC,LWRITE,LTR,LVAR,LDY,LT,L7,L8,L9,L10
   COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
   COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),whatToPlot(NM),CRR(NM)
   COMMON/QNUS/NMSC,LS(NM),MS(NM),NS(NM),CFS(NM)
   COMMON/PAR/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
   COMMON/PARI/RAI,TAI,PRI,PMI,ETAI,CI,OMI,FTWI,FTGI,MFI
   COMMON/NPAR/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD
   COMMON/NPARI/M0I,NTVI,NTHI,LTVI,LTHI,KTVI,KTHI,LEVI,LRBI,LDI
   COMMON/LNMAX/NLMAC,NL,LC(NLMA),NMAXC(NLMA),NMABC
   COMMON/AB/A(NAM),B(NAM),NAMC
   COMMON/NUM/RELE,EPS,ALPH,STOER,NITMAX,NJA
   COMMON/POLE/XP,YP

   NCPLOT = 0
   ZDO = 0.E0

   !-- INPUT:
   LR = 0
   NQ = 0
   NR = 0

   READ(*,*)
   READ(*,*) INPUTFILE,OUTPUTFILE,dataSetNumber,driftRate
   READ(*,*)
   READ(*,*) timeSeriesControl, drawHeader, drawPlotNum, drawTime, &
             plotSize, countourParIsNumber, drawFrame

   !-- drawHeader CONTROLLS THE HEAD:
   !   drawHeader = 0 : NHEAD,
   !   drawHeader = 1 : NHEAD WRITTEN,
   !   drawPlotNum = 0 : NO PLOTNUMBERS
   !   drawPlotNum = 1 : PLOTNUMBERS WITH DESCRIPTION,
   !   drawPlotNum = 2 : PLOTNUMBERS WITHOUT DESCRIPTION,
   !   drawPlotNum = 3 : PLOTNUMBERS GIVEN BY subPlotLabel,
   !   drawTime CONTROLLS WETHER THE TIME IS WRITTEN (0/1).
   !   plotSize DETERMINS THE SIZE OF THE PLOT:
   !   plotSize = 0 : SMALL
   !   plotSize = 1 : MEDIUM
   !   plotSize = 2 : BIG
   !   countourParIsNumber CONTROLLS WETHER contourPar IS THE NUMBER OF CONTOURLINES (countourParIsNumber = .true.) OR
   !   THE DIFFERENCE BETWEEN THE CONTOUR LEVELS (countourParIsNumber = .false.).
   !   FOR drawFrame = 1 A FRAME IS DRAWN AROUND EACH SUBPLOT.
   !   timeSeriesControl CONTROLLS TIMESERIES
   !   timeSeriesControl = 0  : NORMAL
   !   timeSeriesControl = -1 : TIMESERIES OF 6 PLOTS WITH TIME GIVEN INDIVIDUALLY,
   !   timeSeriesControl = 1  : TIMESERIES OF 6 PLOTS WITH TIME GIVEN BY OM,
   !   timeSeriesControl = -2 : TIMESERIES OF 8 PLOTS WITH TIME GIVEN INDIVIDUALLY,
   !   timeSeriesControl = 2  : TIMESERIES OF 8 PLOTS WITH TIME GIVEN BY OM.
   IF( drawPlotNum.LT.0 .OR. drawPlotNum.GT.3 ) THEN
      WRITE(*,*) 'WRONG INPUT OF drawPlotNum: ',drawPlotNum
      STOP
   ENDIF
   IF( plotSize.LT.0 .OR. plotSize.GT.2 ) THEN
      WRITE(*,*) 'WRONG INPUT OF plotSize: ',plotSize
      STOP
   ENDIF
   IF( drawHeader.NE.0 .AND. drawHeader.NE.1 ) THEN
      WRITE(*,*) 'WRONG INPUT OF drawHeader: ',drawHeader
      STOP
   ENDIF
   IF( drawFrame.NE.0 .AND. drawFrame.NE.1 ) THEN
      WRITE(*,*) 'WRONG INPUT OF drawFrame: ',drawFrame
      STOP
   ENDIF
   IF( timeSeriesControl.LT.-2 .OR. timeSeriesControl.GT.2 ) THEN
      WRITE(*,*) 'WRONG INPUT OF timeSeriesControl: ',timeSeriesControl
      STOP
   ENDIF

   OPEN(14,FILE = OUTPUTFILE,STATUS = 'unknown')
   write(14,*) 'Inputfile,dataSetNumber ',INPUTFILE,dataSetNumber

   RELE = 1.D-9
   EPS = 1.D-13

   !-- nPlots IS NUMBER OF PLOTS, XP AND YP ARE LATITUDE AND LONGITUDE OF
   !   THE POLE FOR PROJECTION OF A SPHERE ON A CIRCLE ( quadrant = 'PL','PR','PS' ) .
   READ(*,*)
   READ(*,*) nPlots,XP,YP

   if(nPlots.ne.1.) then
      write(*,*) 'wrong number of plots.'
      stop
   endif

   IF( timeSeriesControl.GT.0 ) THEN
      plotSize = 0
      nPlots = 1
   ENDIF

   IF( ( plotSize.EQ.0 .AND. nPlots.GT.nPlotsMAX ) .OR. &
         ( plotSize.EQ.1 .AND. nPlots.GT.4 )   .OR. &
         ( plotSize.EQ.2 .AND. nPlots.GT.2 ) ) THEN
      WRITE(*,*) 'TOO BIG NUMBER OF PLOTS nPlots: ',nPlots
      STOP
   ENDIF

   DO I = 1,nPlots
      !----- nSubPlots IS NUMBER OF SUBPLOTS, domain DESTINGUISHES BETWEEN
      !      QUADRANT (domain = 'QU'), HALFSPHERE (domain = 'HS') AND FULL SPHERE (domain = 'SP').
      !      TIME IS THE TIME OF THE PLOTTED FIELD FOR TIME DEPENDENCE.
      READ(*,*)
      READ(*,*) domain(I),TIME(I),nSubPlots(I)
      if(nSubPlots(I).ne.1) then
         write(*,*) 'wrong number of plots.'
         stop
      endif

      IF( nSubPlots(I).GT.nSubPlotsMAX ) THEN
         WRITE(*,*) 'TOO BIG NUMBER OF SUBPLOTS nSubPlots:',nSubPlots(I)
         STOP
      ENDIF
      IF( domain(I).EQ.'HS' ) THEN
         NR = NR+1
      ELSE
         NQ = NQ+1
      ENDIF

      !----- quadrant DESTINGUSHES BETWEEN QUADRANTS (quadrant = 'Q1','Q2','Q3','Q4') ,
      !      HALF SPHERES ( quadrant = 'HL','HR','HU','HO') ,SPHERE ( quadrant = 'SP' )
      !      PROJECTION ON A SPHERE ( quadrant = ' PS','PL','PR' ) .
      !      constantCoordinate DETERMINS WETHER R = constantCoordinateValue (constantCoordinate = 'R') , PHI = constantCoordinateValue (constantCoordinate = 'P') OR
      !      THETA = constantCoordinateValue (constantCoordinate = 'T') IS KEPT CONSTANT ,
      !      whatToPlot DETERMINS THE FIELD TO BE PLOTTED:
      !      'VS' : STREAMFUNCTIONS OF VELOCITY FIELD IN BUSSE NOTATION,
      !      'BS' : STREAMFUNCTIONS OF MAGNETIC FIELD IN BUSSE NOTATION,
      !      'JS' : STREAMFUNCTIONS OF ELECTRIC CURRENT IN BUSSE NOTATION,
      !      'VR' : RADIAL VELOCITY FIELD,
      !      'BR' : RADIAL MAGNETIC FIELD,
      !      'TE' : TEMPERATURE FIELD Theta,
      !      'ZF' : ZONAL FLOW ( Mean part of phi comp. of velocity),
      !      'MF' : MERIDIONAL FLOW ( MEAN STREAM FUNCTION IN PHI = CONST. PLANE ),
      !      'MT' : MEAN TOROIDAL MAGNETIC FIELD FOR PHI = CONST,
      !      'BT' : TOROIDAL MAGNETIC FIELD FOR PHI = CONST,
      !      'MP' : STREAMLINES OF MEAN POLOIDAL MAGNETIC FIELD FOR PHI = CONST,
      !      'MJ' : STREAMLINES OF MEAN ELECTRIC CURRENT FOR PHI = CONST,
      !      'MC' : CONTOUR LINES OF MEAN PHI COMPONENT OF ELECTRIC CURRENT FOR PHI = CONST,
      !      'TT' : Temperature field Theta+Ts,
      !      'UP' : Phi component of velocity,
      !      'NU' : local Nusselt number for r = ri.
      !
      !      normRadiusMax IS A MULTIPLIER FOR THE LARGEST RADIUS TO BE PLOTTED: RM = normRadiusMax*RO.
      !      contourPar IS THE STEP FOR THE CONTOURS FOR countourParIsNumber = .false. OR
      !      THE NUMBER OF CONTPUR LINES FOR Z>0 OR Z<0.
      ! Next two lines  are repeated for the number of subplots
      !| SUBPL | PLANE(RPT) | CONST | FIELD |(MAX RAD)/RO|contourPar/STEP|PlotNR|
      !   'SP'      'T'        90      'VS'     1.E0          9    '000'
      DO J = 1,nSubPlots(I)
         READ(*,*)
         READ(*,*) quadrant(I,J), constantCoordinate(I,J), &
                   constantCoordinateValue(I,J), whatToPlot(I,J), &
                   normRadiusMax(I,J),contourPar(I,J),subPlotLabel(I,J)
      enddo
   enddo
   !-- END OF PARAMETER INPUT.

   !-- INPUT CHECK:
   DO I = 1,nPlots
      DO J = 1,nSubPlots(I)
         IF( domain(I).NE.'QU' .AND. domain(I).NE.'HS' .AND.  domain(I).NE.'SP' ) THEN
            WRITE(*,*) 'WRONG INPUT OF domain.'
            WRITE(*,*) 'CHOOSE BETWEEN QUADRANT : domain = QU ,'
            WRITE(*,*) '            HALF SPHERE : domain = HS ,'
            WRITE(*,*) '                 SPHERE : domain = SP .'
            WRITE(*,'('' PLOT '',I3,'' , SUBPLOT '',I3)') I,J
            STOP
         ENDIF
         IF( quadrant(I,J).NE.'Q1' .AND. quadrant(I,J).NE.'Q2' .AND. &
            quadrant(I,J).NE.'Q3' .AND. quadrant(I,J).NE.'Q4' .AND. &
            quadrant(I,J).NE.'HU' .AND. quadrant(I,J).NE.'HO' .AND. &
            quadrant(I,J).NE.'HL' .AND. quadrant(I,J).NE.'HR' .AND. &
            quadrant(I,J).NE.'PL' .AND. quadrant(I,J).NE.'PR' .AND. &
            quadrant(I,J).NE.'SP' .AND. quadrant(I,J).NE.'PS' ) THEN
            WRITE(*,*) 'WRONG INPUT OF quadrant.'
            WRITE(*,*) '  CHOOSE BETWEEN QUADRANTS : quadrant = Q1,Q2,Q3,Q4 ,'
            WRITE(*,*) '              HALF SPHERES : quadrant = HL,HR,HO,HU ,'
            WRITE(*,*) '               FULL SPHERE : quadrant = SP ,'
            WRITE(*,*) '   PROJECTION OF LEFT HALF : quadrant = PL ,'
            WRITE(*,*) '  PROJECTION OF RIGHT HALF : quadrant = PL ,'
            WRITE(*,*) ' PROJECTION OF FULL SPHERE : quadrant = PL .'
            WRITE(*,'('' PLOT '',I3,'' , SUBPLOT '',I3)') I,J
            STOP
         ENDIF
         IF( constantCoordinate(I,J).NE.'P' .AND. constantCoordinate(I,J).NE.'T' .AND. constantCoordinate(I,J).NE.'R' ) THEN
            WRITE(*,*) 'WRONG INPUT OF CONSTANT COORDINATE constantCoordinate.'
            WRITE(*,*) '          CHOOSE BETWEEN PHI : constantCoordinate = P ,'
            WRITE(*,*) '                       THETA : constantCoordinate = T ,'
            WRITE(*,*) '                           R : constantCoordinate = R ,'
            WRITE(*,*) ' RADIAL FIELD FOR CONSTANT R : constantCoordinate = R .'
            WRITE(*,'('' PLOT '',I3,'' , SUBPLOT '',I3)') I,J
            STOP
         ENDIF
         IF( constantCoordinate(I,J).EQ.'R' ) normRadiusMax(I,J) = 1.E0
!
         IF( constantCoordinate(I,J).EQ.'P' .AND. constantCoordinateValue(I,J).GT.360.E0 ) THEN
            WRITE(*,*) 'PHI SHOULD BE GIVEN IN DEGREES < =  360.'
            WRITE(*,'('' PLOT '',I3,'' , SUBPLOT '',I3)') I,J
            STOP
         ELSEIF( constantCoordinate(I,J).EQ.'T' .AND. constantCoordinateValue(I,J).GT.180.E0 ) THEN
            WRITE(*,*) 'THETA SHOULD BE GIVEN IN DEGREES < =  180.'
            WRITE(*,'('' PLOT '',I3,'' , SUBPLOT '',I3)') I,J
            STOP
         ELSEIF(  constantCoordinate(I,J).EQ.'R' .AND. ( constantCoordinateValue(I,J).LT.0.E0 .OR. constantCoordinateValue(I,J).GT.1.E0 )  )THEN
            WRITE(*,*) 'RREL SHOULD BE > = 0 , < =  1 .'
            WRITE(*,*) 'RREL IS DEFINED AS: R = RI+RREL*(RO-RI).'
            WRITE(*,'('' PLOT '',I3,'' , SUBPLOT '',I3)') I,J
            STOP
         ENDIF
         IF(  constantCoordinate(I,J).EQ.'P' ) THEN
            IF( constantCoordinateValue(I,J).LT.0.E0 ) constantCoordinateValue(I,J) = constantCoordinateValue(I,J)+360.E0
            IF(( quadrant(I,J).EQ.'Q1' .OR. &
                 quadrant(I,J).EQ.'Q4' .OR. &
                 quadrant(I,J).EQ.'HR' ) .AND. &
                 constantCoordinateValue(I,J).GT.180.E0 .AND. &
                 constantCoordinateValue(I,J).LT.360.E0 ) THEN
                 constantCoordinateValue(I,J) = constantCoordinateValue(I,J) - 180.E0
            ELSEIF(( quadrant(I,J).EQ.'Q2' .OR. &
                     quadrant(I,J).EQ.'Q3' .OR. &
                     quadrant(I,J).EQ.'HL' ) .AND. &
                     constantCoordinateValue(I,J).LT.180.E0 .AND. &
                     constantCoordinateValue(I,J).GT.0.E0 ) THEN
               constantCoordinateValue(I,J) = constantCoordinateValue(I,J) + 180.E0
            ENDIF
         ENDIF
         IF( constantCoordinate(I,J).NE.'P' .AND. whatToPlot(I,J).EQ.'MT' ) THEN
            WRITE(*,*) 'FOR MT PHI HAS TO BE KEPT CONSTANT.'
            constantCoordinate(I,J) = 'P'
         ENDIF
         IF( constantCoordinate(I,J).NE.'P' .AND. ( whatToPlot(I,J).EQ.'MP' .OR. whatToPlot(I,J).EQ.'BT' ) ) THEN
            WRITE(*,*) 'FOR MP AND BT PHI HAS TO BE KEPT CONSTANT.'
            constantCoordinate(I,J) = 'P'
         ENDIF
         IF( constantCoordinate(I,J).NE.'P' .AND. whatToPlot(I,J).EQ.'MJ' ) THEN
            WRITE(*,*) 'FOR  MJ PHI HAS TO BE KEPT CONSTANT.'
            constantCoordinate(I,J) = 'P'
         ENDIF
         IF( constantCoordinate(I,J).NE.'P' .AND. whatToPlot(I,J).EQ.'MC' ) THEN
            WRITE(*,*) 'FOR  MC PHI HAS TO BE KEPT CONSTANT.'
            constantCoordinate(I,J) = 'P'
         ENDIF
         IF( constantCoordinate(I,J).NE.'P' .AND. whatToPlot(I,J).EQ.'ZF' ) THEN
            WRITE(*,*) 'FOR ZONAL FLOW PHI HAS TO BE KEPT CONSTANT.'
            constantCoordinate(I,J) = 'P'
         ENDIF
         IF( constantCoordinate(I,J).NE.'P' .AND. whatToPlot(I,J).EQ.'MF' ) THEN
            WRITE(*,*) 'FOR MERIDIONAL FLOW PHI HAS TO BE KEPT CONSTANT.'
            constantCoordinate(I,J) = 'P'
         ENDIF
         IF( constantCoordinate(I,J).NE.'R' .AND. whatToPlot(I,J).EQ.'NU' ) THEN
            WRITE(*,*) 'FOR NUSSELT NUMBER R HAS TO BE KEPT CONSTANT.'
            constantCoordinate(I,J) = 'R'
         ENDIF
         if( constantCoordinate(I,J).EQ.'R' .and. quadrant(I,J)(:1).NE.'P' ) then
            write(*,*) 'For R = const the subplot must be a projection.'
            stop
         endif
         IF( domain(I).EQ.'QU' ) THEN
            IF( nSubPlots(I).GT.1 ) THEN
               WRITE(*,*) 'ONLY ONE SUBPLOT ALLOWED FOR domain = QU.'
               WRITE(*,*) 'nSubPlots ASSUMED TO BE 1.'
               nSubPlots(I) = 1
            ENDIF
            IF( constantCoordinate(I,J).NE.'P' .AND. constantCoordinate(I,J).NE.'T' ) THEN
               WRITE(*,*) 'FOR domain = QU ONLY constantCoordinate = P OR constantCoordinate = T VALID.'
               WRITE(*,'('' PLOT '',I3,'' , SUBPLOT '',I3)') I,J
               STOP
            ENDIF
            IF( quadrant(I,J).NE.'Q1' .AND. quadrant(I,J).NE.'Q2' .AND. &
               quadrant(I,J).NE.'Q3' .AND. quadrant(I,J).NE.'Q4' ) THEN
               WRITE(*,*) 'FOR domain = QU ONLY quadrant = Q1,Q2,Q3,Q4 VALID.'
               WRITE(*,'('' PLOT '',I3,'' , SUBPLOT '',I3)') I,J
               STOP
            ENDIF
         ELSEIF( domain(I).EQ.'HS' ) THEN
            IF( nSubPlots(I).GT.2 ) THEN
               WRITE(*,*) 'ONLY TWO SUBPLOT MAXIMUM ALLOWED FOR domain = HS.'
               WRITE(*,*) 'nSubPlots ASSUMED TO BE 2.'
               nSubPlots(I) = 2
            ENDIF
            IF( constantCoordinate(I,J).NE.'P' .AND. constantCoordinate(I,J).NE.'T' ) THEN
               WRITE(*,*) 'FOR domain = HS ONLY constantCoordinate = P OR constantCoordinate = T VALID.'
               WRITE(*,'('' PLOT '',I3,'' , SUBPLOT '',I3)') I,J
               STOP
            ENDIF
            IF( quadrant(I,J).NE.'Q1' .AND. quadrant(I,J).NE.'Q2' .AND. &
                quadrant(I,J).NE.'Q3' .AND. quadrant(I,J).NE.'Q4' .AND. &
                quadrant(I,J).NE.'HO' .AND. quadrant(I,J).NE.'HU' ) THEN
               WRITE(*,*) 'FOR domain = HS ONLY quadrant = Q1,Q2,Q3,Q4,HO,HU VALID.'
               WRITE(*,'('' PLOT '',I3,'' , SUBPLOT '',I3)') I,J
               STOP
            ENDIF
            IF( quadrant(I,J).EQ.'HO' .AND. constantCoordinate(I,J).EQ.'P' ) THEN
               quadrant(I,J) = 'Q1'
               IF( constantCoordinateValue(I,J).GT.180.E0 ) THEN
                  constantCoordinateValue(I,J) = constantCoordinateValue(I,J)-180.E0
               ENDIF
               nSubPlots(I) = nSubPlots(I)+1
               quadrant(I,nSubPlots(I)) = 'Q2'
               constantCoordinate(I,nSubPlots(I)) = constantCoordinate(I,J)
               IF( constantCoordinateValue(I,J).LT.180.E0 ) THEN
                  constantCoordinateValue(I,nSubPlots(I)) = constantCoordinateValue(I,nSubPlots(I))+180.E0
               ENDIF
               whatToPlot(I,nSubPlots(I)) = whatToPlot(I,J)
               normRadiusMax(I,nSubPlots(I)) = normRadiusMax(I,J)
               contourPar(I,nSubPlots(I)) = contourPar(I,J)
            ELSEIF( quadrant(I,J).EQ.'HU' .AND. constantCoordinate(I,J).EQ.'P' ) THEN
               quadrant(I,J) = 'Q4'
                     IF( constantCoordinateValue(I,J).GT.180.E0 ) THEN
                  constantCoordinateValue(I,J) = constantCoordinateValue(I,J)-180.E0
               ENDIF
               nSubPlots(I) = nSubPlots(I)+1
               quadrant(I,nSubPlots(I)) = 'Q3'
               constantCoordinate(I,nSubPlots(I)) = constantCoordinate(I,J)
                     IF( constantCoordinateValue(I,J).LT.180.E0 ) THEN
                  constantCoordinateValue(I,nSubPlots(I)) = constantCoordinateValue(I,nSubPlots(I))+180.E0
               ENDIF
               whatToPlot(I,nSubPlots(I)) = whatToPlot(I,J)
               normRadiusMax(I,nSubPlots(I)) = normRadiusMax(I,J)
               contourPar(I,nSubPlots(I)) = contourPar(I,J)
            ENDIF
         ELSEIF( domain(I).EQ.'SP' ) THEN
            IF( nSubPlots(I).GT.4 ) THEN
            WRITE(*,*) 'ONLY FOUR SUBPLOT MAXIMUM ALLOWED FOR domain = SP.'
               WRITE(*,*) 'nSubPlots ASSUMED TO BE 4.'
               nSubPlots(I) = 4
            ENDIF
            IF( quadrant(I,J).EQ.'HO' .AND. constantCoordinate(I,J).EQ.'P' ) THEN
               quadrant(I,J) = 'Q1'
               IF( constantCoordinate(I,J).EQ.'P' .AND. constantCoordinateValue(I,J).GT.180.E0 ) constantCoordinateValue(I,J) = constantCoordinateValue(I,J)-180.E0
               nSubPlots(I) = nSubPlots(I)+1
               quadrant(I,nSubPlots(I)) = 'Q2'
               constantCoordinate(I,nSubPlots(I)) = constantCoordinate(I,J)
               IF( constantCoordinate(I,J).EQ.'P' .AND. constantCoordinateValue(I,J).LT.180.E0 ) constantCoordinateValue(I,nSubPlots(I)) = constantCoordinateValue(I,J)+180.E0
               whatToPlot(I,nSubPlots(I)) = whatToPlot(I,J)
               normRadiusMax(I,nSubPlots(I)) = normRadiusMax(I,J)
               contourPar(I,nSubPlots(I)) = contourPar(I,J)
            ELSEIF( quadrant(I,J).EQ.'HU' .AND. constantCoordinate(I,J).EQ.'P' ) THEN
               quadrant(I,J) = 'Q4'
               IF( constantCoordinate(I,J).EQ.'P' .AND. constantCoordinateValue(I,J).GT.180.E0 ) constantCoordinateValue(I,J) = constantCoordinateValue(I,J)-180.E0
               nSubPlots(I) = nSubPlots(I)+1
               quadrant(I,nSubPlots(I)) = 'Q3'
               constantCoordinate(I,nSubPlots(I)) = constantCoordinate(I,J)
               IF( constantCoordinate(I,J).EQ.'P' .AND. constantCoordinateValue(I,J).LT.180.E0 ) constantCoordinateValue(I,nSubPlots(I)) = constantCoordinateValue(I,J)+180.E0
               whatToPlot(I,nSubPlots(I)) = whatToPlot(I,J)
               normRadiusMax(I,nSubPlots(I)) = normRadiusMax(I,J)
               contourPar(I,nSubPlots(I)) = contourPar(I,J)
            ENDIF
         ENDIF
      enddo
   enddo
   !-- SETTING OF CONSTANTS:
   LTR = 1
   NMC = NM
   NMSC = NM
   NAMC = NAM
   NLMAC = NLMA

   IF( timeSeriesControl.GT.0 ) TIME(1) = 0.D0

   !-- READLA READS THE SET OF COEFFITIENTS TO BE PLOTTED ,
   !   IT IS NECESSARY TO CALL IS HERE TO GET PARAMETERS.
   TIMEO = TIME(1)
   write(*,'(A,A,I3,D9.2)') 'reading data from ',INPUTFILE,dataSetNumber,TIME(1),'...'
   CALL READLA(INPUTFILE,dataSetNumber,TIME(1),DX)
   write(*,*) '...done'
   RA = RAI
   TA = TAI
   PR = PRI
   PM = PMI
   ETA = ETAI
   C = CI
   OM = OMI
   FTW = FTWI
   FTG = FTGI
   MF = 0
   M0 = M0I
   NTV = NTVI
   NTH = NTHI
   LTV = LTVI
   LTH = LTHI
   KTV = KTVI
   KTH = KTHI
   LD = LDI
   LEV = LEVI
   LRB = LRBI

   IF( timeSeriesControl.GT.0 ) THEN
      NSUBP = nSubPlots(1)
      IF( LCALC.NE.3 .AND. LCALC.NE.4 ) THEN
         WRITE(*,*) 'TIMESERIES WITH timeSeriesControl.GT.0 ONLY FOR TIME EXPANSION.'
         STOP
      ENDIF
      IF( OM.LT.0.D-4 ) THEN
         WRITE(*,*) 'OM TOO SMALL FOR timeSeriesControl.GT.0.'
         STOP
      ENDIF
      IF( timeSeriesControl.EQ.1 ) THEN
         nPlots = 6
      ELSEIF( timeSeriesControl.EQ.2 ) THEN
         nPlots = 8
      ENDIF
      TPERIOD = 2*PI/OM
      DT = TPERIOD/nPlots
      drawPlotNum = 3
      TIME(1) = 0.D0
      subPlotLabel(1,1) = '(a)'
      DO J = 2,nSubPlots(1)
         subPlotLabel(1,J) = '   '
      enddo
      DO I = 2,nPlots
         domain(I) = domain(1)
         nSubPlots(I) = nSubPlots(1)
         IF( domain(I).EQ.'HS' ) THEN
            NR = NR+1
         ELSE
            NQ = NQ+1
         ENDIF
         IF( nPlots.EQ.6 ) THEN
            select case(I)
               case(2)
                  TIME(I) = 1*DT
               case(3)
                  TIME(I) = 5*DT
               case(4)
                  TIME(I) = 2*DT
               case(5)
                  TIME(I) = 4*DT
               case(6)
                  TIME(I) = 3*DT
            end select
         ELSEIF( nPlots.EQ.8 ) THEN
            select case(I)
               case(2)
                  TIME(I) = 1*DT
               case(3)
                  TIME(I) = 2*DT
               case(4)
                  TIME(I) = 7*DT
               case(5)
                  TIME(I) = 3*DT
               case(6)
                  TIME(I) = 6*DT
               case(7)
                  TIME(I) = 5*DT
               case(8)
                  TIME(I) = 4*DT
            end select
         ENDIF
         DO J = 1,nSubPlots(1)
            quadrant(I,J) = quadrant(1,J)
            constantCoordinate(I,J) = constantCoordinate(1,J)
            constantCoordinateValue(I,J) = constantCoordinateValue(1,J)
            whatToPlot(I,J) = whatToPlot(1,J)
            normRadiusMax(I,J) = normRadiusMax(1,J)
            contourPar(I,J) = -contourPar(1,J)
            IF( J.EQ.1 ) THEN
               IF( nPlots.EQ.6 ) THEN
                  select case(I)
                     case(2)
                        subPlotLabel(I,J) = '(b)'
                     case(3)
                        subPlotLabel(I,J) = '(f)'
                     case(4)
                        subPlotLabel(I,J) = '(c)'
                     case(5)
                        subPlotLabel(I,J) = '(e)'
                     case(6)
                        subPlotLabel(I,J) = '(d)'
                     case default
                        subPlotLabel(I,J) = '   '
                  end select
               ELSEIF( nPlots.EQ.8 ) THEN
                  select case(I)
                     case(2)
                        subPlotLabel(I,J) = '(b)'
                     case(3)
                        subPlotLabel(I,J) = '(c)'
                     case(4)
                        subPlotLabel(I,J) = '(h)'
                     case(5)
                        subPlotLabel(I,J) = '(d)'
                     case(6)
                        subPlotLabel(I,J) = '(g)'
                     case(7)
                        subPlotLabel(I,J) = '(f)'
                     case(8)
                        subPlotLabel(I,J) = '(e)'
                     case default
                        subPlotLabel(I,J) = '   '
                  end select
               ENDIF
            endif
         enddo
      enddo
   ENDIF

   !-- CALCULATION OF INNER AND OUTER RADIUS:
   RI = ETA/(1.D0-ETA)
   RO = 1.D0+RI
   XRI = DBLE(RI)
   XRO = DBLE(RO)

   !-- ABG CALCULATES THE ALPHAS AND BETAS IN THE RADIAL FUNCTION
   !   OF THE POLOIDAL MAGNETIC FIELD:
   IF( LCALC.EQ.2 .OR. LCALC.EQ.4 .OR. LCALC.EQ.6 ) CALL ABG(ND,whatToPlot,L,N)

   XLRAND = 3.0D0
   XRRAND = 3.0D0
   NROWR = NR
   NROWQ = NQ/2
   IF( MOD(NQ,2).NE.0 ) NROWQ = NROWQ+1
   IF( timeSeriesControl.NE.0 ) THEN
      NROWR = 0
      NROWQ = 3
   ELSE
      NROWR = NR
      NROWQ = NQ/2
      IF( MOD(NQ,2).NE.0 ) NROWQ = NROWQ+1
   ENDIF
   NROW = NROWR+NROWQ
   IF( plotSize.EQ.0 .AND. NROW.GT.3 ) THEN
      WRITE(*,*) 'TOO MANY ROWS, ONLY 3 ALLOWED FOR plotSize = 0.',NROW
      STOP
   ELSEIF( plotSize.EQ.1 .AND. NROW.GT.2 ) THEN
      WRITE(*,*) 'TOO MANY ROWS, ONLY 2 ALLOWED FOR plotSize = 1.',NROW
      STOP
   ELSEIF( plotSize.EQ.2 .AND. NROW.GT.1 ) THEN
      WRITE(*,*) 'TOO MANY ROWS, ONLY 1 ALLOWED FOR plotSize = 2.',NROW
      STOP
   ENDIF
   IF( NQ.LE.1 ) THEN
      NCOL = 1
   ELSE
      NCOL = 2
   ENDIF
   IF( plotSize.EQ.2 .AND. NCOL.EQ.2 ) THEN
      WRITE(*,*) 'TOO MANY COLUMNS, ONLY 1 ALLOWED FOR plotSize = 2.',NCOL
      STOP
   ENDIF
   IF( IABS(timeSeriesControl).EQ.2 ) NCOL = 3
   XTEXT = XLRAND

   !-- GROESSE DER PLOTS:
   IF( plotSize.EQ.0 ) THEN
      YHR = 5.5D0
      YHQ = 5.5D0
      XLR = 2.0D0
      XLQ = 1.5D0
      XINTER = 1.D0
      YINTER = 1.D0
   ELSEIF( plotSize.EQ.1 ) THEN
      YHR = 6.75D0
      YHQ = 6.75D0
      XLR = 0.75D0
      XLQ = 0.0D0
      XINTER = 1.5D0
      YINTER = 1.5D0
   ELSEIF( plotSize.EQ.2 ) THEN
      YHR = 6.75D0
      YHQ = 13.0
      XLR = 1.0D0
      XINTER = 0.0D0
      YINTER = 0.0D0
   ENDIF
   XBQ = YHQ
   XBR = 2*YHR

   IF( drawHeader.EQ.0 ) THEN
      headerSpaceY = 0.0D0
   ELSEIF( drawHeader.EQ.1 ) THEN
      headerSpaceY = 3.D0
   ENDIF
   IF( drawPlotNum.EQ.1 ) THEN
      plotNumSpaceY = 5.0D0
   ELSE
      plotNumSpaceY = 0.0D0
   ENDIF
   YHPG = NROWR*YHR+NROWQ*YHQ+(NROW-1)*YINTER
   IF( NQ.GT.0 ) THEN
      XBPG = 2*XLQ+NCOL*XBQ+(NCOL-1)*XINTER
   ELSE
      XBPG = 2*XLR+XBR
   ENDIF
   XAREA = XLRAND+XBPG+XRRAND
   YAREA = plotNumSpaceY+YHPG+headerSpaceY

   !-- NZEI ZAEHLT ZEILEN , NSPA SPALTEN UND NP ZAHL DER PLOTS.
   !   NQT ZAEHLT DIE ZAHL DER QUADRATE.
   NP = 0
   NQT = 0

   !-- DIE DATEN FUER DIE EINZELNEN PLOTS WERDEN FESTGELEGT UND LINEAR
   !   ABGESPEICHERT: URSPRUNG IN CM = (XORIG,YORIG) ,
   !   PLOTGEBIET IN CM = (XAR,YAR) , RADIEN IN CM = (XRICM,XROCM,XRMCM).
   DO I = 1,nPlots
      IF( I.EQ.1 ) THEN
         NSPA = 1
         NZEI = 1
         IF( domain(I).EQ.'HS' ) THEN
            YHPLOT = YHR
            XBPLOT = XBR
            XLPLOT = XLR
         ELSE
            YHPLOT = YHQ
            XBPLOT = XBQ
            XLPLOT = XLQ
         ENDIF
      ELSEIF( domain(I).EQ.'HS' .OR. domain(I-1).EQ.'HS' ) THEN
         NSPA = 1
         NZEI = NZEI+1
         YHPLOT = YHR
         XBPLOT = XBR
         XLPLOT = XLR
      ELSE
         IF( NSPA.EQ.NCOL ) THEN
            NSPA = 1
            NZEI = NZEI+1
         ELSE
            NSPA = NSPA+1
            IF( IABS(timeSeriesControl).EQ.2 .AND. I.EQ.5 ) NSPA = NSPA+1
         ENDIF
         YHPLOT = YHQ
         XBPLOT = XBQ
         XLPLOT = XLQ
      ENDIF
      XORIG = XLRAND+XLPLOT+(NSPA-1)*(XBPLOT+XINTER)
      YORIG = plotNumSpaceY+YHPG-NZEI*YHPLOT-(NZEI-1)*YINTER

      IF( domain(I).EQ.'QU' .OR. domain(I).EQ.'SP' .AND. NCOL.GT.1 ) THEN
         NQT = NQT+1
         IF( NSPA.EQ.1 .AND. NQT.EQ.NQ ) XLQ = XLQ+(XBQ+XINTER)/2
      ENDIF
      DO J = 1,nSubPlots(I)
         NP = NP+1
         thisPlotQuadrant(NP)           = quadrant(I,J)
         thisPlotWhatToPlot(NP)         = whatToPlot(I,J)
         thisPlotconstantCoordinate(NP) = constantCoordinate(I,J)
         ABCN(NP) = subPlotLabel(I,J)
         IF( constantCoordinate(I,J).EQ.'R' ) THEN
            XCP(NP) = XRI+constantCoordinateValue(I,J)
         ELSE
            XCP(NP) = constantCoordinateValue(I,J)
         ENDIF
         ZDP(NP)   = contourPar(I,J)
         TIMEP(NP) = TIME(I)
         IF( domain(I).EQ.'HS' ) THEN
            IF( quadrant(I,J).EQ.'HO' .OR. quadrant(I,J).EQ.'HU' ) THEN
               XOR(NP) = XORIG
               XAR(NP) = XBR
               YOR(NP) = YORIG
               YAR(NP) = YHR
               XRMCM(NP) = XBR/2
            ELSEIF( quadrant(I,J).EQ.'Q1' ) THEN
               XOR(NP) = XORIG+XBR/2
               XAR(NP) = XBR/2
               YOR(NP) = YORIG
               YAR(NP) = YHR
               XRMCM(NP) = XBR/2
            ELSEIF( quadrant(I,J).EQ.'Q2' ) THEN
               XOR(NP) = XORIG
               XAR(NP) = XBR/2
               YOR(NP) = YORIG
               YAR(NP) = YHR
               XRMCM(NP) = XBR/2
            ELSEIF( quadrant(I,J).EQ.'Q3' ) THEN
               XOR(NP) = XORIG
               XAR(NP) = XBR/2
               YOR(NP) = YORIG
               YAR(NP) = YHR
               XRMCM(NP) = XBR/2
            ELSEIF( quadrant(I,J).EQ.'Q4' ) THEN
               XOR(NP) = XORIG+XBR/2
               XAR(NP) = XBR/2
               YOR(NP) = YORIG
               YAR(NP) = YHR
               XRMCM(NP) = XBR/2
            ENDIF
         ELSEIF( domain(I).EQ.'QU' ) THEN
            XOR(NP) = XORIG
            XAR(NP) = XBQ
            YOR(NP) = YORIG
            YAR(NP) = YHQ
            XRMCM(NP) = XBQ
         ELSEIF( domain(I).EQ.'SP' ) THEN
            IF( quadrant(I,J).EQ.'Q1' ) THEN
               XOR(NP) = XORIG+XBQ/2
               XAR(NP) = XBQ/2
               YOR(NP) = YORIG+YHQ/2
               YAR(NP) = YHQ/2
               XRMCM(NP) = XBQ/2
            ELSEIF( quadrant(I,J).EQ.'Q2' ) THEN
               XOR(NP) = XORIG
               XAR(NP) = XBQ/2
               YOR(NP) = YORIG+YHQ/2
               YAR(NP) = YHQ/2
               XRMCM(NP) = XBQ/2
            ELSEIF( quadrant(I,J).EQ.'Q3' ) THEN
               XOR(NP) = XORIG
               XAR(NP) = XBQ/2
               YOR(NP) = YORIG
               YAR(NP) = YHQ/2
               XRMCM(NP) = XBQ/2
            ELSEIF( quadrant(I,J).EQ.'Q4' ) THEN
               XOR(NP) = XORIG+XBQ/2
               XAR(NP) = XBQ/2
               YORIG = YORIG
               YOR(NP) = YORIG
               YAR(NP) = YHQ/2
               XRMCM(NP) = XBQ/2
            ELSEIF( quadrant(I,J).EQ.'HU' ) THEN
               XOR(NP) = XORIG
               XAR(NP) = XBQ
               YOR(NP) = YORIG
               YAR(NP) = YHQ/2
               XRMCM(NP) = XBQ/2
            ELSEIF( quadrant(I,J).EQ.'HO' ) THEN
               XOR(NP) = XORIG
               XAR(NP) = XBQ
               YOR(NP) = YORIG+YHQ/2
               YAR(NP) = YHQ/2
               XRMCM(NP) = XBQ/2
            ELSEIF( quadrant(I,J).EQ.'HL' ) THEN
               XOR(NP) = XORIG
               XAR(NP) = XBQ/2
               YOR(NP) = YORIG
               YAR(NP) = YHQ
               XRMCM(NP) = XBQ/2
            ELSEIF( quadrant(I,J).EQ.'HR' ) THEN
               XOR(NP) = XORIG+XBQ/2
               XAR(NP) = XBQ/2
               YOR(NP) = YORIG
               YAR(NP) = YHQ
               XRMCM(NP) = XBQ/2
            ELSEIF( quadrant(I,J).EQ.'SP' .OR. quadrant(I,J).EQ.'PS' .OR. quadrant(I,J).EQ.'PL' ) THEN
               XOR(NP) = XORIG
               XAR(NP) = XBQ
               YOR(NP) = YORIG
               YAR(NP) = YHQ
               XRMCM(NP) = XBQ/2
            ELSEIF( quadrant(I,J).EQ.'PR' ) THEN
               XOR(NP) = XORIG
               XAR(NP) = XBQ
               YOR(NP) = YORIG
               YAR(NP) = YHQ
               XRMCM(NP) = XBQ/2
            ENDIF
         ENDIF
         XRM(NP) = normRadiusMax(I,J)*XRO
         XROCM(NP) = XRMCM(NP)*XRO/XRM(NP)
         XRICM(NP) = XRMCM(NP)*XRI/XRM(NP)
      enddo
   enddo

   YTEXT = plotNumSpaceY-1.0E0

   !-- PLO FUEHRT DIE EINZELNEN SUBPLOTS AUS:

   DO I = 1,NP
      write(14,*) 'Plot Nr. ',I,':'

      IF( TIMEP(I).NE.TIMEO .AND. LT.EQ.1 ) CALL READLA(INPUTFILE,dataSetNumber,TIMEP(I),DX)

      CALL PLO(I,NSUBP,driftRate,DX,countourParIsNumber,plotSize,drawFrame,&
                     ZDP(I),TIMEP(I),thisPlotQuadrant(I),thisPlotWhatToPlot(I), &
                     thisPlotconstantCoordinate(I),XCP(I),        &
                     XOR(I),YOR(I),XAR(I),YAR(I),                       &
                     XRI,XRO,XRM(I),XRICM(I),XROCM(I),XRMCM(I),XP,YP)

      !-- BESCHRIFTUNG SUBPLOT: BUCHSTABE IM PLOT AN POSITION (XNUM,YNUM) ,
      !   BESCHRIFTUNG UNTER BUCHSTABE UNTEN AUF SEITE.
      IF( drawPlotNum.GT.0 ) THEN
         IF( thisPlotQuadrant(I).EQ.'Q1' ) THEN
            XNUM = XOR(I)+XAR(I)-0.9E0
            YNUM = YOR(I)+YAR(I)-0.5E0
            XTIME = XNUM-2.E0
            YTIME = YNUM+0.8E0
         ELSEIF( thisPlotQuadrant(I).EQ.'Q2' ) THEN
            XNUM = XOR(I)+0.25E0
            YNUM = YOR(I)+YAR(I)-0.5E0
            XTIME = XNUM
            YTIME = YNUM+0.8E0
         ELSEIF( thisPlotQuadrant(I).EQ.'Q3' ) THEN
            XNUM = XOR(I)+0.25E0
            YNUM = YOR(I)+0.2E0
            XTIME = XNUM
            YTIME = YNUM-1.0E0
         ELSEIF( thisPlotQuadrant(I).EQ.'Q4' ) THEN
            XNUM = XOR(I)+XAR(I)-0.9E0
            YNUM = YOR(I)+0.3E0
            XTIME = XNUM-2.E0
            YTIME = YNUM-1.0E0
         ELSEIF( thisPlotQuadrant(I).EQ.'HO' ) THEN
            XNUM = XOR(I)+0.25E0
            YNUM = YOR(I)+YAR(I)-0.5E0
            XTIME = XNUM
            YTIME = YNUM+0.8E0
         ELSEIF( thisPlotQuadrant(I).EQ.'HU' ) THEN
            XNUM = XOR(I)+0.25E0
            YNUM = YOR(I)+0.2E0
            XTIME = XNUM
            YTIME = YNUM-1.0E0
         ELSEIF( thisPlotQuadrant(I).EQ.'HL' .OR. thisPlotQuadrant(I).EQ.'PL' ) THEN
            XNUM = XOR(I)+0.25E0
            YNUM = YOR(I)+YAR(I)-0.5E0
            XTIME = XNUM
            YTIME = YNUM+0.8E0
         ELSEIF( thisPlotQuadrant(I).EQ.'HR' .OR. thisPlotQuadrant(I).EQ.'PR' ) THEN
            XNUM = XOR(I)+XAR(I)-0.9E0
            YNUM = YOR(I)+YAR(I)-0.5E0
            XTIME = XNUM-2.E0
            YTIME = YNUM+0.8E0
         ELSEIF( thisPlotQuadrant(I).EQ.'SP' .OR. thisPlotQuadrant(I).EQ.'PS' ) THEN
            XNUM = XOR(I)+0.25E0
            YNUM = YOR(I)+YAR(I)-0.5E0
            XTIME = XNUM
            YTIME = YNUM+0.8E0
         ENDIF
         IF( drawPlotNum.EQ.1 ) THEN
            IF( thisPlotWhatToPlot(I).EQ.'BR' ) THEN
               CTEXT1 = '  Contours of radial magnetic field for '
            ELSEIF( thisPlotWhatToPlot(I).EQ.'BR' ) THEN
               CTEXT1 = '  Contours of toriodal magnetic field for '
            ELSEIF( thisPlotWhatToPlot(I).EQ.'BS' ) THEN
               CTEXT1 = '  Magnetic field lines in the plane '
            ELSEIF( thisPlotWhatToPlot(I).EQ.'VR' ) THEN
               CTEXT1 = '  Contours of radial velocity field for '
            ELSEIF( thisPlotWhatToPlot(I).EQ.'VS' ) THEN
               CTEXT1 = '  Streamlines in the plane  '
            ELSEIF( thisPlotWhatToPlot(I).EQ.'VS' ) THEN
               CTEXT1 = '  Streamlines of electric current in the plane  '
            ELSEIF( thisPlotWhatToPlot(I).EQ.'TE' ) THEN
               CTEXT1 = '  Temperature field '
            ELSEIF( thisPlotWhatToPlot(I).EQ.'ZF' ) THEN
               CTEXT1 = '  Mean zonal flow '
            ELSEIF( thisPlotWhatToPlot(I).EQ.'MF' ) THEN
               CTEXT1 = '  Mean meridional flow '
            ELSEIF( thisPlotWhatToPlot(I).EQ.'MT' ) THEN
               CTEXT1 = '  Mean toroidal magnetic field '
            ELSEIF( thisPlotWhatToPlot(I).EQ.'MP' ) THEN
               CTEXT1 = '  Fieldlines of mean poloidal magnetic field '
            ELSEIF( thisPlotWhatToPlot(I).EQ.'MJ' ) THEN
               CTEXT1 = '  Fieldlines of mean electric current'
            ELSEIF( thisPlotWhatToPlot(I).EQ.'MC' ) THEN
               CTEXT1 = '  Contourlines of mean phi comp. of elec. curr.'
            ELSEIF( thisPlotWhatToPlot(I).EQ.'TT' ) THEN
               CTEXT1 = '  Temperature field including basic state for '
            ELSEIF( thisPlotWhatToPlot(I).EQ.'UP' ) THEN
               CTEXT1 = '  Contours of U{M7}v{M0} for '
            ELSEIF( thisPlotWhatToPlot(I).EQ.'NU' ) THEN
               CTEXT1 = '  Contours of local nusselt number for '
               XCP(I) = REAL(ETA/(1.D0-ETA))
            ENDIF
            IF( thisPlotconstantCoordinate(I).EQ.'P') THEN
               CTEXT2 = 'phi  = '
            ELSEIF( thisPlotconstantCoordinate(I).EQ.'T') THEN
               CTEXT2 = 'theta  = '
            ELSEIF( thisPlotconstantCoordinate(I).EQ.'R' ) THEN
               CTEXT2 = 'r  = '
            ENDIF
            write(14,*) CTEXT1,CTEXT2,XCP(I)
            write(*,*) CTEXT1,CTEXT2,XCP(I)
         ENDIF
         YTEXT = YTEXT-0.5E0
      ENDIF
      TIMEO = TIMEP(I)
   enddo
   CLOSE(14)

contains

   !------------------------------------------------------------------------
   !     calculates the field Z and makes one subplot.
   SUBROUTINE PLO(NPLOT,NSUBP,driftRate,DX,countourParIsNumber,plotSize,drawFrame,contourPar,TIME,quadrant,whatToPlot,constantCoordinate,constantCoordinateValue, &
             XOR,YOR,XAR,YAR,XRI,XRO,XRM,XRICM,XROCM,XRMCM,XP,YP)
      IMPLICIT REAL*8(A-H,O-W)
      IMPLICIT REAL*8(X,Y,Z)
      CHARACTER*1 constantCoordinate,CCC
      CHARACTER*2 quadrant,whatToPlot,CPC
      character*20 filez,filex,filey

      DIMENSION DX(*),THETA(NMY),XIDL(NMX,NMY),YIDL(NMX,NMY)
      DIMENSION Z(NMX,NMY),XML(2),YML(2),ZDS(4)

      COMMON/PARI/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/CNULL/ZNULL
      integer on_a_sphere

      !-- COUNTER
      NCPLOT = NCPLOT+1
      WRITE(14,'(2X,''TIME =  '',2D16.6)') TIME

      !-- UEBERGABE AN COMMONBLOCK FUER TRANS:
      CCC = constantCoordinate
      CPC = quadrant
      XCC = constantCoordinateValue

      !-- INITIALISIERUNG VON DISSPLA UND ZEICHNEN EINES RAHMENS (FRAME):
      DXY = XRO/100

      !-- FESTLEGEN DER X BZW Y ACHSE UND ZEICHNEN DES INNEREN UND
      !   AEUSSEREN KERNS MIT ARC:
      IF( quadrant.EQ.'Q1' ) THEN
         XML(1) = XRI
         YML(1) = 0.E0
         XML(2) = XRO
         YML(2) = 0.E0
         XML(1) = 0.E0
         YML(1) = XRI
         XML(2) = 0.E0
         YML(2) = XRO
         XMIN = XRI
         XMAX = XRM
         IF( constantCoordinate.EQ.'T' ) THEN
            YMIN = 90.E0
            YMAX = 180.E0
         ELSEIF( constantCoordinate.EQ.'P' ) THEN
            YMIN = 0.E0
            YMAX = 90.E0
         ENDIF
      ELSEIF( quadrant.EQ.'Q2' ) THEN
         XMIN = XRI
         XMAX = XRM
         IF( constantCoordinate.EQ.'T' ) THEN
            YMIN = 180.E0
            YMAX = 270.E0
         ELSEIF( constantCoordinate.EQ.'P' ) THEN
            YMIN = 0.E0
            YMAX = 90.E0
         ENDIF
         XML(1) = -XRO
         YML(1) = 0.E0
         XML(2) = -XRI
         YML(2) = 0.E0
         XML(1) = 0.E0
         YML(1) = XRI
         XML(2) = 0.E0
         YML(2) = XRO
      ELSEIF( quadrant.EQ.'Q3' ) THEN
         XMIN = XRI
         XMAX = XRM
         IF( constantCoordinate.EQ.'T' ) THEN
            YMIN = 270.E0
            YMAX = 360.E0
         ELSEIF( constantCoordinate.EQ.'P' ) THEN
            YMIN = 90.E0
            YMAX = 180.E0
         ENDIF
         XML(1) = -XRO
         YML(1) = 0.E0
         XML(2) = -XRI
         YML(2) = 0.E0
         XML(1) = 0.E0
         YML(1) = -XRI
         XML(2) = 0.E0
         YML(2) = -XRO
      ELSEIF( quadrant.EQ.'Q4' ) THEN
         XMIN = XRI
         XMAX = XRM
         IF( constantCoordinate.EQ.'T' ) THEN
            YMIN = 0.E0
            YMAX = 90.E0
         ELSEIF( constantCoordinate.EQ.'P' ) THEN
            YMIN = 90.E0
            YMAX = 180.E0
         ENDIF
         XML(1) = XRI
         YML(1) = 0.E0
         XML(2) = XRO
         YML(2) = 0.E0
         XML(1) = 0.E0
         YML(1) = -XRI
         XML(2) = 0.E0
         YML(2) = -XRO
      ELSEIF( quadrant.EQ.'HO' ) THEN
         XMIN = XRI
         XMAX = XRM
         IF( constantCoordinate.EQ.'T' ) THEN
            YMIN = 90.E0-XP
            YMAX = 270.E0-XP
         ELSEIF( constantCoordinate.EQ.'P' ) THEN
            YMIN = 0.E0
            YMAX = 90.E0
         ENDIF
         XML(1) = -XRO
         YML(1) = 0.E0
         XML(2) = -XRI
         YML(2) = 0.E0
         XML(1) = XRI
         YML(1) = 0.E0
         XML(2) = XRO
         YML(2) = 0.E0
      ELSEIF( quadrant.EQ.'HU' ) THEN
         XMIN = XRI
         XMAX = XRM
         IF( constantCoordinate.EQ.'T' ) THEN
            YMIN = 270.E0-XP
            YMAX = 450.E0-XP
         ELSEIF( constantCoordinate.EQ.'P' ) THEN
            YMIN = 90.E0
            YMAX = 180.E0
         ENDIF
         XML(1) = -XRO
         YML(1) = 0.E0
         XML(2) = -XRI
         YML(2) = 0.E0
         XML(1) = XRI
         YML(1) = 0.E0
         XML(2) = XRO
         YML(2) = 0.E0
      ELSEIF( quadrant.EQ.'HL' ) THEN
         XMIN = XRI
         XMAX = XRM
         IF( constantCoordinate.EQ.'T' ) THEN
            YMIN = 180.E0
            YMAX = 360.E0
         ELSEIF( constantCoordinate.EQ.'P' ) THEN
            YMIN = 0.E0
            YMAX = 180.E0
         ENDIF
         XML(1) = 0.D0
         YML(1) = -XRO
         XML(2) = 0.D0
         YML(2) = -XRI
         XML(1) = 0.E0
         YML(1) = XRI
         XML(2) = 0.E0
         YML(2) = XRO
      ELSEIF( quadrant.EQ.'HR' ) THEN
         XMIN = XRI
         XMAX = XRM
         IF( constantCoordinate.EQ.'T' ) THEN
            YMIN = 0.E0
            YMAX = 180.E0
         ELSEIF( constantCoordinate.EQ.'P' ) THEN
            YMIN = 0.E0
            YMAX = 180.E0
         ENDIF
         XML(1) = 0.E0
         YML(1) = -XRO
         XML(2) = 0.E0
         YML(2) = -XRI
         XML(1) = 0.E0
         YML(1) = XRI
         XML(2) = 0.E0
         YML(2) = XRO
      ELSEIF( quadrant.EQ.'SP' ) THEN
         XMIN = XRI
         XMAX = XRM
         IF( constantCoordinate.EQ.'T' ) THEN
            YMIN = 0.E0
            YMAX = 360.E0
         ELSEIF( constantCoordinate.EQ.'P' ) THEN
            YMIN = 0.E0
            YMAX = 180.E0
         ENDIF
      ELSEIF( quadrant.EQ.'PS' ) THEN
         XMIN = -180.E0-XP
         XMAX = 180.E0-XP
         YMIN = 180.E0-YP
         YMAX = 0.E0-YP
      ELSEIF( quadrant.EQ.'PL' ) THEN
         XMIN = -180.E0-XP
         XMAX = 180.E0-XP
         YMIN = 180.E0-YP
         YMAX = 0.E0-YP
         XML(1) = 0.5E0
         YML(1) = 0.E0
         XML(2) = 0.5E0
         YML(2) = 1.E0
      ELSEIF( quadrant.EQ.'PR' ) THEN
         XMIN = -180.E0-XP
         XMAX = 180.E0-XP
         YMIN = 180.E0-YP
         YMAX = 0.E0-YP
         XML(1) = 0.5E0
         YML(1) = 0.E0
         XML(2) = 0.5E0
         YML(2) = 1.E0
      ENDIF

      !-- IST DER MAXIMALE RADIUS XRM GROESSER ALS DER AEUSSERE RADIUS RO
      !   UND EXISTIERT ABER NUR FUER R< = RO EIN FELD , SO MUESSEN DAS
      !   PLOTGEBIET UND DER URSPRUNG ENTSPRECHEND ANGEPASST WERDEN.
      !   TEILWEISE WIRD ZUDEM DIE X-ACHSE AUF R< = RO EINGESCHRAENKT.
      IF( ( ( whatToPlot.NE.'BS' .AND. whatToPlot.NE.'MP' ) .OR. &
            ( whatToPlot.EQ.'BS' .AND. constantCoordinate.EQ.'R' ) ) .AND. &
            XROCM.NE.XRMCM ) THEN
         IF( quadrant(:1).EQ.'Q' ) THEN
            XAR = XROCM
            YAR = XROCM
            XMAX = XRO
            IF( quadrant.EQ.'Q2' ) THEN
               XOR = XOR+XRMCM-XROCM
            ELSEIF( quadrant.EQ.'Q3' ) THEN
               XOR = XOR+XRMCM-XROCM
               YOR = YOR+XRMCM-XROCM
            ELSEIF( quadrant.EQ.'Q4' ) THEN
               YOR = YOR+XRMCM-XROCM
            ENDIF
         ELSEIF( quadrant.EQ.'HL' ) THEN
            XAR = XROCM
            YAR = 2*XROCM
            XMAX = XRO
            XOR = XOR+XRMCM-XROCM
            YOR = YOR+XRMCM-XROCM
         ELSEIF( quadrant.EQ.'HR' ) THEN
            XAR = XROCM
            YAR = 2*XROCM
            XMAX = XRO
            YOR = YOR+XRMCM-XROCM
         ELSEIF( quadrant.EQ.'HO' ) THEN
            XAR = 2*XROCM
            YAR = XROCM
            XMAX = XRO
            XOR = XOR+XRMCM-XROCM
         ELSEIF( quadrant.EQ.'HU' ) THEN
            XAR = 2*XROCM
            YAR = XROCM
            XMAX = XRO
            XOR = XOR+XRMCM-XROCM
            YOR = YOR+XRMCM-XROCM
         ELSEIF( quadrant.EQ.'SP' .OR. quadrant(:1).EQ.'P' ) THEN
            XAR = 2*XROCM
            YAR = 2*XROCM
            IF( quadrant.EQ.'SP' ) XMAX = XRO
            XOR = XOR+XRMCM-XROCM
            YOR = YOR+XRMCM-XROCM
         ENDIF
      ENDIF
      XARC = XAR
      YARC = YAR

      write(*,*) 'computing the fields...'

      !-- BERECHNEN DER Z-WERTE FUER EIN RASTER MIT JE NXM PUNKTEN IN
      !   X-RICHTUNG UND NYM PUNKTEN IN Y-RICHTUNG:
      !   THETA WIRD EIN INTEGER NTHETA ZUGEORDNET UNTER DEM PLM(THETA)
      !   ABGESPEICHERT WIRD, NMTHETA IST DIE ANZAHL DER BENOETIGTEN THETA.
      IF( constantCoordinate.EQ.'T' ) THEN
         NMTHETA = 1
         THETA(NMTHETA) = DBLE(constantCoordinateValue)
      ELSEIF( constantCoordinate.EQ.'P' ) THEN
         PHI = DBLE(constantCoordinateValue)
      ELSEIF( constantCoordinate.EQ.'R' ) THEN
         R = DBLE(constantCoordinateValue)
      ENDIF
      XD = (XMAX-XMIN)/(NMX-1)
      YD = (YMAX-YMIN)/(NMY-1)
      IF( constantCoordinate.NE.'T' ) THEN
         NMTHETA = NMY
         DO I = 1, NMTHETA
            THETA(I) = DBLE(YMIN+(I-1)*YD)
         enddo
      ENDIF

      !-- BESTIMMUNG DER PLM(THETA) , ABSPEICHERUNG:
      CALL STOREPLM(THETA,NMTHETA)

      ZMIN = 1.E10
      ZMAX = -1.E10
      DO I = 1,NMX
         X = XMIN+(I-1)*XD
         DO J = 1,NMY
            Y = YMIN+(J-1)*YD
            IF( constantCoordinate.EQ.'T' ) THEN
               R = DBLE(X)
               PHI = DBLE(Y)
               NTHETA = 1
            ELSEIF( constantCoordinate.EQ.'P' ) THEN
               R = DBLE(X)
               NTHETA = J
            ELSEIF( constantCoordinate.EQ.'R' ) THEN
               PHI = DBLE(X)
               NTHETA = J
            ENDIF
            !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ! IDL
            !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if( constantCoordinate.eq.'T' ) then
               XIDL(I,J) = R*COS(pi*PHI/180.d0)
               YIDL(I,J) = R*SIN(pi*PHI/180.d0)
               on_a_sphere = 0
            elseif( constantCoordinate.eq.'P' ) then
               XIDL(I,J) = R*COS(pi*(THETA(J)-90.d0)/180.d0)
               YIDL(I,J) = R*SIN(pi*(THETA(J)-90.d0)/180.d0)
               on_a_sphere = 0
            elseif( constantCoordinate.eq.'R' ) then
               on_a_sphere = 1
            else
               on_a_sphere = 0
               write(*,*) 'wrong constant variable: ',cc
               stop
            endif
            !-------- R,PHI UND THETA SIND DIE KUGELKOORDINATEN:
            IF( whatToPlot.EQ.'VS' .OR. whatToPlot.EQ.'BS' .OR. whatToPlot.EQ.'JS' ) THEN
               IF(  constantCoordinate.EQ.'T' ) THEN
                  Z(I,J) = REAL(FT(DX,whatToPlot,R,PHI,NTHETA,TIME,driftRate))
               ELSEIF( constantCoordinate.EQ.'P' ) THEN
                  Z(I,J) = REAL(FP(DX,whatToPlot,R,PHI,NTHETA,TIME,driftRate))
               ELSEIF( constantCoordinate.EQ.'R' ) THEN
                  Z(I,J) = REAL(FR(DX,whatToPlot,R,PHI,NTHETA,TIME,driftRate))
               ENDIF
            ELSEIF( whatToPlot.EQ.'VR' .OR. whatToPlot.EQ.'BR' ) THEN
               Z(I,J) = REAL(RF(DX,whatToPlot,R,PHI,NTHETA,TIME,driftRate))
            ELSEIF( whatToPlot.EQ.'TE' ) THEN
               Z(I,J) = REAL(TEMP(DX,whatToPlot,R,PHI,NTHETA,TIME,driftRate))
            ELSEIF( whatToPlot.EQ.'ZF' ) THEN
               Z(I,J) = REAL(FZONAL(DX,R,NTHETA,TIME))
            ELSEIF( whatToPlot.EQ.'MF' ) THEN
               Z(I,J) = REAL(DMERI(DX,R,NTHETA,TIME))
            ELSEIF( whatToPlot.EQ.'MT' ) THEN
               Z(I,J) = REAL(DMTOR(DX,R,NTHETA,TIME))
            ELSEIF( whatToPlot.EQ.'MP' .OR. whatToPlot.EQ.'MJ' ) THEN
               Z(I,J) = REAL(DMPJ(DX,whatToPlot,R,NTHETA,TIME))
            ELSEIF( whatToPlot.EQ.'BT' ) THEN
               Z(I,J) = REAL(DBT(DX,R,PHI,NTHETA,TIME,driftRate))
            ELSEIF( whatToPlot.EQ.'MC' ) THEN
               Z(I,J) = REAL(DMC(DX,R,NTHETA,TIME))
            ELSEIF( whatToPlot.EQ.'TT' ) THEN
               Z(I,J) = REAL(TT(DX,whatToPlot,R,PHI,NTHETA,TIME,driftRate))
            ELSEIF( whatToPlot.EQ.'UP' ) THEN
               Z(I,J) = REAL(UP(DX,whatToPlot,R,PHI,NTHETA,TIME,driftRate))
            ELSEIF( whatToPlot.EQ.'NU' ) THEN
               Z(I,J) = REAL(FNU(DX,whatToPlot,R,PHI,NTHETA,TIME,driftRate))
            ELSE
               WRITE(*,*) 'WRONG INPUT OF whatToPlot.'
               STOP
            ENDIF
            IF( Z(I,J).GT.ZMAX ) ZMAX = Z(I,J)
            IF( Z(I,J).LT.ZMIN ) ZMIN = Z(I,J)
         enddo
      enddo

      range = MAX(ABS(ZMIN),ABS(ZMAX))
      ZNULL = 1.E-11*range
      ZNULLM = 1.E-11
      ZANULL = 1.E-13
      ZSCALE = 1.E0

      IF( contourPar.GT.0.E0 ) THEN
         IF( range.LT.ZANULL ) THEN
            WRITE(14,*) 'ZMAX AND ZMIN CLOSE TO ZERO: ',ZMAX,ZMIN
            WRITE(14,*) 'NO PLOT POSSIBLE.'
            GOTO 9000
         ELSEIF( ZNULL.LE.ZNULLM ) THEN
            ZSCALE = 1.E0/ZNULLM
            WRITE(14,*) 'SCALED BY ',ZSCALE
            ZMIN = ZSCALE*ZMIN
            ZMAX = ZSCALE*ZMAX
            ZNULL = ZSCALE*ZNULL
            DO IX = 1,NMX
               DO IY = 1,NMY
                  Z(IX,IY) = ZSCALE*Z(IX,IY)
               enddo
            enddo
         ENDIF
      ELSEIF( contourPar.LT.0.E0 ) THEN
         IF( NCPLOT.GT.NSUBP ) THEN
            NZD = MOD(NCPLOT,NSUBP)
            IF( NZD.EQ.0 ) NZD = NSUBP
            contourPar = ZSCALE*ZDS(NZD)
            countourParIsNumber = .true.
         ELSE
            contourPar = ZSCALE*ABS(contourPar)
         ENDIF
      ENDIF
      IF( countourParIsNumber ) THEN
         NCL = AINT(contourPar+0.1E0)
         IF( ZMIN.GT.-ZNULL .OR. ZMAX.LT.ZNULL ) THEN
            contourPar = ((ZMAX-ZMIN)-ZNULL)/(NCL-1)
         ELSE
            contourPar = (MAX(ABS(ZMIN),ABS(ZMAX))-ZNULL)/(NCL-1)
         ENDIF
         contourPar = contourPar-contourPar/100
      ELSE
         IF( contourPar.GT.ABS(ZMAX) .AND. contourPar.GT.ABS(ZMIN) ) THEN
            WRITE(*,*) 'TOO LARGE contourPar , ZMIN,ZMAX ARE ONLY: ',ZMIN,ZMAX
            STOP
         ENDIF
         IF( ZMIN.GT.-ZNULL .OR. ZMAX.LT.ZNULL ) THEN
            NCL = AINT( ABS(ZMAX-ZMIN)/contourPar+0.1E0 )+1
         ELSE
            NCL = AINT(MAX(ABS(ZMAX),ABS(ZMIN))/contourPar+0.1E0)+1
         ENDIF
      ENDIF
      IF( NCPLOT.LE.NSUBP ) ZDS(NCPLOT) = contourPar/ZSCALE
      IF( ZMIN.GT.-ZNULL .OR. ZMAX.LT.ZNULL ) THEN
         ZMINP = ZMIN
         ZMAXP = ZMAX
      ELSE
         ZMINP = contourPar*AINT( ZMIN/contourPar )-ZNULL
         ZMAXP = contourPar*AINT( ZMAX/contourPar )+ZNULL
      ENDIF
      IF( ABS(ZMINP).LT.ZNULL ) ZMINP = ZNULL
      IF( ABS(ZMAXP).LT.ZNULL ) ZMAXP = ZNULL
      WRITE(14,*) 'DIFFERENCE BETWEEN CONTOURLINES contourPar =  ',contourPar
      WRITE(14,*) 'NUMBER OF CONTOURLINES NCL =  ',NCL
      WRITE(14,*) 'ZMAX,ZMIN =  ',ZMAX,ZMIN
      WRITE(14,*) 'ZMAXP,ZMINP =  ',ZMAXP,ZMINP

      TENSN = 0.D0
9000  CONTINUE
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! IDL
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*) 'writing files idl.z, idl.x, idl.y ...'

      filez = 'idl.z'
      filex = 'idl.x'
      filey = 'idl.y'

      open(21,file = filez,STATUS =  'UNKNOWN')
      open(22,file = filex,STATUS =  'UNKNOWN')
      open(23,file = filey,STATUS =  'UNKNOWN')

      if (on_a_sphere.eq.1) then
         DO I = 1,NMX
            X = XMIN+(I-1)*XD
            write(22,*) DBLE(X) + 180.
         enddo
         DO J = 1,NMY
            write(23,*) THETA(J)-90.
         enddo
         do i = 1,nmx
            do j = 1,nmy
               write(21,*) z(i,j)
            enddo
         enddo
      ELSE
         do i = 1,nmx
            do j = 1,nmy
               write(21,*) z(i,j)
               write(22,*) xidl(i,j)
               write(23,*) yidl(i,j)
            enddo
         enddo
      ENDIF

      close(21)
      close(22)
      close(23)
9999  CONTINUE
   END subroutine plo

   !------------------------------------------------------------------------
   !   Stromfunktion fuer theta = konstant:
   !      F_theta = r dphi v             (Busse: r/sin(th) d/dphi v )
   !   Fuer den elektrischen Strom:
   !              F_theta = r dphi g
   !
   !     optimized for K = 0.
   FUNCTION FT(X,whatToPlot,R,PHI,NTHETA,TIME,driftRate)
      IMPLICIT REAL*8(A-H,O-Z)
      double precision, intent(in):: x(:)
      character(len=2), intent(in):: whatToPlot
      integer, intent(in):: NTHETA
      CHARACTER*1 whatToPlot
      CHARACTER*2 CRR

      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),whatToPlot(NM),CRR(NM)
      COMMON/PARI/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPARI/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD
      COMMON/AB/A(NAM),B(NAM),NAMC

      IF( NAM.NE.NAMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NAM IN FT.'
         STOP
      ENDIF
      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN FT.'
         STOP
      ENDIF

      PPHI = PHI*PI/180.D0
      RI = ETA/(1.D0-ETA)
      RO = RI+1.D0
      FT = 0.D0
      IF( whatToPlot.EQ.'VS' ) THEN
         NDOMIN = 1
         NDOMAX = NDV
      ELSEIF( whatToPlot.EQ.'BS' ) THEN
         NDOMIN = NDV+NDW+NDT+1
         NDOMAX = NDV+NDW+NDT+NDH
      ELSEIF( whatToPlot.EQ.'JS' ) THEN
         NDOMIN = NDV+NDW+NDT+NDH+1
         NDOMAX = NDV+NDW+NDT+NDH+NDG
      ELSE
         WRITE(*,*) 'WRONG whatToPlot IN FT, SHOULD BE VS OR BS OR JS BUT IS: ',whatToPlot
         STOP
      ENDIF
      DO  I = NDOMIN,NDOMAX
         IF( .NOT.( ( whatToPlot(I).EQ.'V' .AND. whatToPlot.EQ.'VS' ) .OR. &
                    ( whatToPlot(I).EQ.'H' .AND. whatToPlot.EQ.'BS' ) .OR. &
                    ( whatToPlot(I).EQ.'G' .AND. whatToPlot.EQ.'JS' )  )  ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN FT, SHOULD BE V OR H OR G BUT IS: ',whatToPlot(I)
            STOP
         ENDIF
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         ENDIF
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         ENDIF
         FTT = EPSM*EPSK*M(I)*PLMS(L(I),M(I),NTHETA)*R
         IF( whatToPlot(I).EQ.'V' .OR. whatToPlot(I).EQ.'G' ) THEN
            FTT = FTT*DSIN( N(I)*PI*(R-RI) )
         ELSEIF( whatToPlot(I).EQ.'H' ) THEN
            NR = NAB(L(I),N(I))
            IF( R.LE.RO ) THEN
               IF( A(NR).EQ.0.D0 .OR. B(NR).EQ.0.D0 ) THEN
                  WRITE(*,*) 'ALPHA AND BETA NOT CALCULATED.'
                  STOP
               ENDIF
               FTT = -FTT*DCOS( A(NR)*R-B(NR) )
            ELSE
               IF( A(NR).EQ.0.D0 .OR. B(NR).EQ.0.D0 ) THEN
                  WRITE(*,*) 'ALPHA AND BETA NOT CALCULATED.'
                  STOP
               ENDIF
               FTT = -FTT * (RO/R)**(L(I)+1) * DCOS( A(NR)*RO-B(NR) )
            ENDIF
         ENDIF

         if(K(I).EQ.0) then
            IF( CRR(I).EQ.'RR' ) THEN
               FTT = -FTT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               FTT = -FTT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) )
            ELSE
               FTT = 0.D0
            ENDIF
         else
            IF( CRR(I).EQ.'RR' ) THEN
               FTT = -FTT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               FTT = -FTT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               FTT = FTT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'II' ) THEN
               FTT = FTT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
            ENDIF
         endif
         FT = FT-FTT
      enddo
   END function ft

   !------------------------------------------------------------------------
   ! Stromfunktion fuer phi = konstant:
   !              F_phi = r sin(theta) dtheta v  (like Busse)
   ! Fuer den elektrischen Strom:
   !              F_phi = r sin(theta) dtheta g
   !
   !     optimized for K = 0.
   FUNCTION FP(X,whatToPlot,R,PHI,NTHETA,TIME,driftRate)
      IMPLICIT REAL*8(A-H,O-Z)
      double precision, intent(in):: x(:)
      integer, intent(in):: NTHETA
      CHARACTER*1 whatToPlot
      CHARACTER*2 CRR,whatToPlot

      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),whatToPlot(NM),CRR(NM)
      COMMON/PARI/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
!cc   COMMON/NPARI/M0,NE,NTV,NTH,LTV,LTH,KTV,KTH,LD
      COMMON/NPARI/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD
      COMMON/AB/A(NAM),B(NAM),NAMC

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN FP.'
         STOP
      ENDIF
      IF( NAM.NE.NAMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NAM IN FP.'
         STOP
      ENDIF

      PPHI = PHI*PI/180.D0
      RI = ETA/(1.D0-ETA)
      RO = RI+1.D0
      FP = 0.D0
      IF( whatToPlot.EQ.'VS' ) THEN
         NDOMIN = 1
         NDOMAX = NDV

      ELSEIF( whatToPlot.EQ.'BS' ) THEN
         NDOMIN = NDV+NDW+NDT+1
         NDOMAX = NDV+NDW+NDT+NDH
      ELSEIF( whatToPlot.EQ.'JS' ) THEN
         NDOMIN = NDV+NDW+NDT+NDH+1
         NDOMAX = NDV+NDW+NDT+NDH+NDG
      ELSE
         WRITE(*,*) 'WRONG whatToPlot IN FP, SHOULD BE VS OR BS OR JS BUT IS: ',whatToPlot
         STOP
      ENDIF
      DO I = NDOMIN,NDOMAX
         IF( .NOT.( ( whatToPlot(I).EQ.'V' .AND. whatToPlot.EQ.'VS' ) .OR.&
                   ( whatToPlot(I).EQ.'H' .AND. whatToPlot.EQ.'BS' ) .OR. &
                   ( whatToPlot(I).EQ.'G' .AND. whatToPlot.EQ.'JS' )  )  ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN FP, SHOULD BE V OR H OR G BUT IS: ',whatToPlot(I)
            STOP
         ENDIF
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         ENDIF
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         ENDIF
         FPT = EPSM*EPSK*R * ( &
                 DBLE(L(I))*DSQRT( DBLE( (L(I)-M(I)+1)*(L(I)+M(I)+1) ) /  &
                 DBLE( (2*L(I)+1)*(2*L(I)+3) ) ) * PLMS(L(I)+1,M(I),NTHETA) -     &
                 DBLE(L(I)+1)*DSQRT( DBLE( (L(I)-M(I))*(L(I)+M(I)) ) /    &
                 DBLE( (2*L(I)+1)*(2*L(I)-1) ) ) * PLMS(L(I)-1,M(I),NTHETA)  )

         IF( whatToPlot(I).EQ.'V' .OR. whatToPlot(I).EQ.'G' ) THEN
            FPT = FPT*DSIN( N(I)*PI*(R-RI) )
         ELSEIF( whatToPlot(I).EQ.'H' ) THEN
            NR = NAB(L(I),N(I))
            IF( R.LE.RO ) THEN
               IF( A(NR).EQ.0.D0 .OR. B(NR).EQ.0.D0 ) THEN
                  WRITE(*,*) 'ALPHA AND BETA NOT CALCULATED.'
                  STOP
               ENDIF
               FPT = FPT*DCOS( A(NR)*R-B(NR) )
            ELSE
               IF( A(NR).EQ.0.D0 .OR. B(NR).EQ.0.D0 ) THEN
                  WRITE(*,*) 'ALPHA AND BETA NOT CALCULATED.'
                  STOP
               ENDIF
               FPT = FPT * (RO/R)**(L(I)+1) * DCOS( A(NR)*RO-B(NR) )
            ENDIF
         ENDIF

         if(K(I).EQ.0) then
            IF( CRR(I).EQ.'RR' ) THEN
               FPT = FPT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               FPT = -FPT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) )
            ELSE
               FPT = 0.D0
            ENDIF
         else
            IF( CRR(I).EQ.'RR' ) THEN
               FPT = FPT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               FPT = -FPT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               FPT = -FPT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'II' ) THEN
               FPT = FPT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
            ENDIF
         endif
         FP = FP+FPT
      enddo
   END function FP

   !------------------------------------------------------------------------
   !   Stromfunktion fuer r = konstant:
   !                     F_r = w      (like Busse, Hirsching: rw )
   !   Stromfunktion fuer r = konstant des elektrischen Stroms:
   !                     F_r = - laplace h
   !
   !     optimized for K = 0.
   FUNCTION FR(X,whatToPlot,R,PHI,NTHETA,TIME,driftRate)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*1 whatToPlot
      CHARACTER*2 CRR,whatToPlot

      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),whatToPlot(NM),CRR(NM)
      COMMON/PARI/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
!cc   COMMON/NPARI/M0,NE,NTV,NTH,LTV,LTH,KTV,KTH,LD
      COMMON/NPARI/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD
      COMMON/AB/A(NAM),B(NAM),NAMC

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN FR.'
         STOP
      ENDIF

      PPHI = PHI*PI/180.D0
      RI = ETA/(1.D0-ETA)
      FR = 0.D0
      IF( whatToPlot.EQ.'VS' ) THEN
         NDOMIN = NDV+1
         NDOMAX = NDV+NDW
      ELSEIF( whatToPlot.EQ.'JS' ) THEN
         NDOMIN = NDV+NDW+NDT+1
         NDOMAX = NDV+NDW+NDT+NDH
      ELSEIF( whatToPlot.EQ.'BS' ) THEN
         NDOMIN = NDV+NDW+NDT+NDH+1
         NDOMAX = NDV+NDW+NDT+NDH+NDG
      ELSE
         WRITE(*,*) 'WRONG whatToPlot IN FR, SHOULD BE VS OR BS OR JS BUT IS: ',whatToPlot
         STOP
      ENDIF
      DO I = NDOMIN,NDOMAX
         IF(  .NOT.( ( whatToPlot.EQ.'VS' .AND. whatToPlot(I).EQ.'W' ) .OR.&
                     ( whatToPlot.EQ.'BS' .AND. whatToPlot(I).EQ.'G' ) .OR.&
                     ( whatToPlot.EQ.'JS' .AND. whatToPlot(I).EQ.'H' )  )  ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN FR, SHOULD BE W OR G OR H BUT IS: ',whatToPlot(I)
            STOP
         ENDIF
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         ENDIF
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         ENDIF
         FRT = EPSM*EPSK*PLMS(L(I),M(I),NTHETA)
         IF( whatToPlot(I).EQ.'W' ) THEN
            FRT = FRT*R*DCOS( (N(I)-1)*PI*(R-RI) )
         ELSEIF( whatToPlot(I).EQ.'G' ) THEN
            FRT = FRT*DSIN( N(I)*PI*(R-RI) )
         ELSEIF( whatToPlot(I).EQ.'H' ) THEN
            NR = NAB(L(I),N(I))
            IF( A(NR).EQ.0.D0 .OR. B(NR).EQ.0.D0 ) THEN
               WRITE(*,*) 'ALPHA AND BETA NOT CALCULATED.'
               STOP
            ENDIF
            FRT = FRT*( ( A(NR)*A(NR)+DBLE(L(I)*(L(I)+1))/(R*R) ) * DCOS( A(NR)*R-B(NR) ) + 2*A(NR)/R * DSIN( A(NR)*R-B(NR) )  )
         ENDIF
         if(K(I).EQ.0) then
            IF( CRR(I).EQ.'RR' ) THEN
               FRT = FRT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               FRT = -FRT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) )
            ELSE
               FRT = 0.D0
            ENDIF
         else
            IF( CRR(I).EQ.'RR' ) THEN
               FRT = FRT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               FRT = -FRT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               FRT = -FRT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'II' ) THEN
               FRT = FRT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
            ENDIF
         endif
         FR = FR+FRT
      enddo
   END function fr

   !------------------------------------------------------------------------
   ! Radiales Geschw.feld: U_r = L_2/r v
   !
   !     optimized for K = 0.
   FUNCTION RF(X,whatToPlot,R,PHI,NTHETA,TIME,driftRate)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*1 whatToPlot
      CHARACTER*2 CRR,whatToPlot

      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),whatToPlot(NM),CRR(NM)
      COMMON/PARI/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPARI/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD
      COMMON/AB/A(NAM),B(NAM),NAMC

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN RF.'
         STOP
      ENDIF
      IF( NAM.NE.NAMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NAM IN RF.'
         STOP
      ENDIF

      PPHI = PHI*PI/180.D0
      RI = ETA/(1.D0-ETA)
      RF = 0.D0
      IF( whatToPlot.EQ.'VR' ) THEN
         NDOMIN = 1
         NDOMAX = NDV
      ELSEIF( whatToPlot.EQ.'BR' ) THEN
         NDOMIN = NDV+NDW+NDT+1
         NDOMAX = NDV+NDW+NDT+NDH
      ELSE
         WRITE(*,*) 'WRONG whatToPlot IN RF, SHOULD BE V OR H BUT IS: ',whatToPlot
         STOP
      ENDIF
      DO I = NDOMIN,NDOMAX
         IF( .NOT.( ( whatToPlot(I).EQ.'V' .AND. whatToPlot.EQ.'VR' ) .OR. &
                    ( whatToPlot(I).EQ.'H' .AND. whatToPlot.EQ.'BR' )  ) ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN RF, SHOULD BE V OR H BUT IS: ', whatToPlot(I)
            STOP
         ENDIF
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         ENDIF
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         ENDIF
         RFT = EPSM*EPSK*L(I)*(L(I)+1) * PLMS(L(I),M(I),NTHETA) / R
         IF( whatToPlot(I).EQ.'V' ) THEN
            RFT = RFT*DSIN( N(I)*PI*(R-RI) )
         ELSEIF( whatToPlot(I).EQ.'H' ) THEN
            NR = NAB(L(I),N(I))
            IF( A(NR).EQ.0.D0 .OR. B(NR).EQ.0.D0 ) THEN
               WRITE(*,*) 'ALPHA AND BETA NOT CALCULATED.'
               STOP
            ENDIF
            RFT = RFT*DCOS( A(NR)*R-B(NR) )
         ENDIF

         if(K(I).EQ.0) then
            IF( CRR(I).EQ.'RR' ) THEN
               RFT = RFT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               RFT = -RFT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) )
            ELSE
               RFT = 0.D0
            ENDIF
         else
            IF( CRR(I).EQ.'RR' ) THEN
               RFT = RFT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               RFT = -RFT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               RFT = -RFT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'II' ) THEN
               RFT = RFT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
            ENDIF
         endif
         RF = RF+RFT
      enddo
   END function rf

   !------------------------------------------------------------------------
   !   Temperaturfeld Theta ( =  Abweichung vom Grundzust.)
   !   optimized for K = 0.
   FUNCTION TEMP(X,whatToPlot,R,PHI,NTHETA,TIME,driftRate)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 whatToPlot
      CHARACTER*2 CRR,whatToPlot
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)

      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),whatToPlot(NM),CRR(NM)
      COMMON/PARI/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPARI/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD
      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN TEMP.'
         STOP
      ENDIF

      TEMP = 0.D0
      PPHI = PHI*PI/180.D0
      RI = ETA/(1.D0-ETA)
      IF( whatToPlot.EQ.'TE' ) THEN
         NDOMIN = 1+NDV+NDW
         NDOMAX = NDV+NDW+NDT
      ELSE
         WRITE(*,*) 'WRONG whatToPlot IN TEMP, SHOULD BE TE BUT IS: ',whatToPlot
         STOP
      ENDIF

      DO I = NDOMIN, NDOMAX
         IF( whatToPlot(I).NE.'T' ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN TEMP, SHOULD BE T BUT IS: ', whatToPlot(I)
            STOP
         ENDIF
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         ENDIF
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         ENDIF
         TEM = EPSM*EPSK*PLMS(L(I),M(I),NTHETA)*DSIN( N(I)*PI*(R-RI) )

         if(K(I).EQ.0) then
            IF( CRR(I).EQ.'RR' ) THEN
               TEM = TEM * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               TEM = -TEM * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) )
            ELSE
               TEM = 0.D0
            ENDIF
         else
            IF( CRR(I).EQ.'RR' ) THEN
               TEM = TEM * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               TEM = -TEM * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               TEM = -TEM * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'II' ) THEN
               TEM = TEM * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
            ENDIF
         endif
         TEMP = TEMP+TEM
      enddo
   END function temp

   !------------------------------------------------------------------------
   !     temperature field Theta + Ts
   !     optimized for K = 0.
   FUNCTION TT(X,whatToPlot,R,PHI,NTHETA,TIME,driftRate)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*1 whatToPlot
      CHARACTER*2 CRR,whatToPlot

      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),whatToPlot(NM),CRR(NM)
      COMMON/PARI/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPARI/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN TT.'
         STOP
      ENDIF

      PPHI = PHI*PI/180.D0
      RI = ETA/(1.D0-ETA)
      T = 0.D0
      IF( whatToPlot.EQ.'TT' ) THEN
         NDOMIN = NDV+NDW+1
         NDOMAX = NDV+NDW+NDT
      ELSE
         WRITE(*,*) 'WRONG whatToPlot IN TT, SHOULD TT BUT IS: ',whatToPlot
         STOP
      ENDIF
      DO I = NDOMIN,NDOMAX
         IF( whatToPlot(I).NE.'T' ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN T, SHOULD BE T BUT IS: ', whatToPlot(I)
            STOP
         ENDIF
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         ENDIF
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         ENDIF
         TT = EPSM*EPSK*PLMS(L(I),M(I),NTHETA)
         TT = TT*DSIN( N(I)*PI*(R-RI) )

         IF(K(I).EQ.0) THEN
          IF( CRR(I).EQ.'RR' ) THEN
            TT = TT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) )
          ELSEIF( CRR(I).EQ.'IR' ) THEN
            TT = -TT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) )
          ELSEIF( CRR(I).EQ.'RI' ) THEN
            TT = 0.0D0
          ELSEIF( CRR(I).EQ.'II' ) THEN
            TT = 0.0D0
          ENDIF
         ELSE
          IF( CRR(I).EQ.'RR' ) THEN
            TT = TT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) )  * DCOS(K(I)*OM*TIME)
          ELSEIF( CRR(I).EQ.'IR' ) THEN
            TT = -TT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
          ELSEIF( CRR(I).EQ.'RI' ) THEN
            TT = -TT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
          ELSEIF( CRR(I).EQ.'II' ) THEN
            TT = TT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
          ENDIF
         ENDIF
         T = T+TT
      enddo

      !  add basic temperature field Ts:
      T = T - R * R / ( 2.D0 * PR )
      TT = T
   END

   !------------------------------------------------------------------------
   !   local Nusselt number NU(r = ri)
   !   optimized for K = 0.
   FUNCTION FNU(X,whatToPlot,R,PHI,NTHETA,TIME,driftRate)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*1 whatToPlot
      CHARACTER*2 CRR,whatToPlot

      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),whatToPlot(NM),CRR(NM)
      COMMON/PARI/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPARI/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD
      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN FNU.'
         STOP
      ENDIF

      FNU = 0.D0
      PPHI = PHI*PI/180.D0
      RI = ETA/(1.D0-ETA)
      IF( whatToPlot.EQ.'NU' ) THEN
         NDOMIN = 1+NDV+NDW
         NDOMAX = NDV+NDW+NDT
      ELSE
         WRITE(*,*) 'WRONG whatToPlot IN FNU, SHOULD BE NU BUT IS: ',whatToPlot
         STOP
      ENDIF

      DO I = NDOMIN,NDOMAX
         IF( whatToPlot(I).NE.'T' ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN TEMP, SHOULD BE T BUT IS: ', whatToPlot(I)
            STOP
         ENDIF
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         ENDIF
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         ENDIF
         FNUT = EPSM*EPSK*PLMS(L(I),M(I),NTHETA)*DBLE(N(I))*PI

         if(K(I).EQ.0) then
            IF( CRR(I).EQ.'RR' ) THEN
               FNUT = FNUT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               FNUT = -FNUT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) )
            ELSE
               FNUT = 0.D0
            ENDIF
         else
            IF( CRR(I).EQ.'RR' ) THEN
               FNUT = FNUT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               FNUT = -FNUT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               FNUT = -FNUT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'II' ) THEN
               FNUT = FNUT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
            ENDIF
         endif
         FNU = FNU+FNUT
      enddo
      FNU = 1.D0 - PR/RI*FNU
   END function fnu

   !------------------------------------------------------------------------
   !   Zonaler Fluss = gemittelte phi-Komponente der Geschwindigkeit:
   !          < u_phi > = - dtheta w   (m = 0)
   !
   !     optimized for K = 0.
   FUNCTION FZONAL(X,R,NTHETA,TIME)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*1 whatToPlot
      CHARACTER*2 CRR

      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),whatToPlot(NM),CRR(NM)
      COMMON/PARI/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPARI/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN FZONAL.'
         STOP
      ENDIF

      FZONAL = 0.D0
      RI = ETA/(1.D0-ETA)
      NDOMIN = 1+NDV
      NDOMAX = NDV+NDW

      DO I = NDOMIN,NDOMAX
         IF( whatToPlot(I).NE.'W' ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN FZONAL, SHOULD BE W BUT IS: ', whatToPlot(I)
            STOP
         ENDIF

         IF( M(I).NE.0 ) cycle

         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         ENDIF
         ZON = EPSK*DSQRT(DBLE(L(I)*(L(I)+1))) * PLMS(L(I),1,NTHETA) * R * DCOS( (N(I)-1)*PI*(R-RI) )

         if(K(I).EQ.0) then
            IF( CRR(I).EQ.'RR' ) THEN
               ZON = ZON * X(I)
            ELSE
               ZON = 0.D0
            ENDIF
         else
            IF( CRR(I).EQ.'RR' ) THEN
               ZON = ZON * X(I) * DCOS(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               ZON = -ZON * X(I) * DSIN(K(I)*OM*TIME)
            ELSE
               ZON = 0.D0
            ENDIF
         endif
         FZONAL = FZONAL+ZON
      enddo
   END function FZONAL

   !------------------------------------------------------------------------
   !     Uphi = 1/(r*sinphi) d^2/drdph rv - d/dth w
   !
   !     optimized for K = 0.
   FUNCTION UP(X,whatToPlot,R,PHI,NTHETA,TIME,driftRate)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*1 whatToPlot
      CHARACTER*2 CRR,whatToPlot
      DIMENSION THETA(NMY)

      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),whatToPlot(NM),CRR(NM)
      COMMON/PARI/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPARI/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN UP.'
         STOP
      ENDIF

      THETAR = PI*THETA(NTHETA)/180.D0
      SINTH = DSIN(THETAR)

      PPHI = PHI*PI/180.D0
      RI = ETA/(1.D0-ETA)
      UP = 0.D0

      IF( whatToPlot.NE.'UP' ) THEN
        WRITE(*,*) 'WRONG whatToPlot IN UP, SHOULD BE UP BUT IS: ',whatToPlot
        STOP
      ENDIF

      NDOMIN = NDV+1
      NDOMAX = NDV+NDW
      !------- toroidal part: --------------------------
      DO I = NDOMIN,NDOMAX
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         ENDIF
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         ENDIF

         IF(SINTH.EQ.0.D0) THEN
            UPT = 0.D0
         ELSE
            DL = DBLE(L(I))
            DM = DBLE(M(I))
            DLPM = DL+DM
            DLMM = DL-DM
            !-------               -d/dth w       ----------------
            UPT = EPSM*EPSK/SINTH * ( (DL+1.D0)*PLMS(L(I)-1,M(I),NTHETA) * &
               DSQRT(DLPM*DLMM/((2.D0*DL-1)*(2D0*DL+1D0))) - &
               DL*PLMS(L(I)+1,M(I),NTHETA) * &
               DSQRT((DLMM+1.D0)*(DLPM+1.D0)/((2D0*DL+3D0)*(2D0*DL+1D0))) )
         ENDIF

         IF( whatToPlot(I).EQ.'W' ) THEN
            UPT = UPT*R*DCOS( (N(I)-1)*PI*(R-RI) )
         ELSEIF( whatToPlot(I).EQ.'G' ) THEN
            UPT = UPT*DSIN( N(I)*PI*(R-RI) )
         ENDIF
         IF(K(I).EQ.0) THEN
            IF( CRR(I).EQ.'RR' ) THEN
               UPT = UPT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               UPT = -UPT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) )
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               UPT = 0.0D0
            ELSEIF( CRR(I).EQ.'II' ) THEN
               UPT = 0.0D0
            ENDIF
         ELSE
            IF( CRR(I).EQ.'RR' ) THEN
               UPT = UPT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               UPT = -UPT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               UPT = -UPT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'II' ) THEN
               UPT = UPT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
            ENDIF
         ENDIF

         UP = UP+UPT
      enddo

      !------- poloidal part: --------------------------
      NDOMIN = 1
      NDOMAX = NDV
      DO I = NDOMIN,NDOMAX
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
!           poloidal part is 0 for M = 0:
            cycle
         ELSE
            EPSM = 2.D0
         ENDIF
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         ENDIF

         !------- 1/(rsinth) d^2/drdphi (rv) -------------
         IF( SINTH.EQ.0.D0 .OR. M(I).EQ.0 ) THEN
          UPT = 0.D0
         ELSE
            UPT = EPSM*EPSK * M(I) * PLMS(L(I),M(I),NTHETA)  / (R*SINTH)

            UPT = UPT*(R*N(I)*PI*DCOS(N(I)*PI*(R-RI))+DSIN(N(I)*PI*(R-RI)))

            IF(K(I).EQ.0) THEN
               IF( CRR(I).EQ.'RR' ) THEN
                  UPT = -UPT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) )
               ELSEIF( CRR(I).EQ.'IR' ) THEN
                  UPT = -UPT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) )
               ELSE
                  UPT = 0.D0
               ENDIF
            ELSE
               IF( CRR(I).EQ.'RR' ) THEN
                  UPT = -UPT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
               ELSEIF( CRR(I).EQ.'IR' ) THEN
                  UPT = -UPT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
               ELSEIF( CRR(I).EQ.'RI' ) THEN
                  UPT = UPT * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
               ELSEIF( CRR(I).EQ.'II' ) THEN
                  UPT = UPT * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
               ENDIF
            ENDIF
         ENDIF

         UP = UP+UPT
      enddo
   END function UP

   !------------------------------------------------------------------------
   !   Meridionale Zirkulation = phi-gemittelt Stromlinien fuer phi = kostant:
   !        < F_phi > = < r sin(theta) dtheta v>
   !
   !     optimized for K = 0.
   FUNCTION DMERI(X,R,NTHETA,TIME)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*1 whatToPlot
      CHARACTER*2 CRR

      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),whatToPlot(NM),CRR(NM)
      COMMON/PARI/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPARI/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN DMERI.'
         STOP
      ENDIF
      RI = ETA/(1.D0-ETA)
      DMERI = 0.D0
      NDOMIN = 1
      NDOMAX = NDV

      DO I = NDOMIN, NDOMAX
         IF( whatToPlot(I).NE.'V' ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN FP, SHOULD BE V BUT IS: ', whatToPlot(I)
            STOP
         ENDIF

         IF( M(I).NE.0 ) cycle

         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         ENDIF

         DMER = EPSK*R * DBLE(L(I)*(L(I)+1)) * ( &
           1.D0/DSQRT( DBLE( (2*L(I)+1)*(2*L(I)+3) ) ) * PLMS(L(I)+1,0,NTHETA) - &
           1.D0/DSQRT( DBLE( (2*L(I)+1)*(2*L(I)-1) ) ) * PLMS(L(I)-1,0,NTHETA)  )
         DMER = DMER*DSIN( N(I)*PI*(R-RI) )

         if(K(I).EQ.0) then
            IF( CRR(I).EQ.'RR' ) THEN
               DMER = DMER * X(I)
            ELSE
               DMER = 0.D0
            ENDIF
         else
            IF( CRR(I).EQ.'RR' ) THEN
               DMER = DMER * X(I) * DCOS(K(I)*OM*TIME)
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               DMER = -DMER * X(I) * DSIN(K(I)*OM*TIME)
            ELSE
               DMER = 0.D0
            ENDIF
         endif
         DMERI = DMERI+DMER
      enddo
   END function DMERI

   !------------------------------------------------------------------------
   !   phi-gemittelte phi-Komponente der Toroidalfeldes:
   !           < B_phi > = - dtheta g (m = 0)
   !
   !     optimized for K = 0.
   FUNCTION DMTOR(X,R,NTHETA,TIME)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*1 whatToPlot
      CHARACTER*2 CRR

      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),whatToPlot(NM),CRR(NM)
      COMMON/PARI/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPARI/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN DMTOR.'
         STOP
      ENDIF

      RI = ETA/(1.D0-ETA)
      DMTOR = 0.D0
      NDOMIN = NDV+NDW+NDT+NDH+1
      NDOMAX = NDV+NDW+NDT+NDH+NDG
      DO I = NDOMIN,NDOMAX
         IF( whatToPlot(I).NE.'G' ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN DMTOR, SHOULD BE G BUT IS: ', whatToPlot(I)
            STOP
         ENDIF
         IF( M(I).NE.0 ) cycle
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         ENDIF
         DMT = EPSK*DSQRT(DBLE(L(I)*(L(I)+1))) * PLMS(L(I),1,NTHETA) * DSIN( N(I)*PI*(R-RI) )
        if(K(I).EQ.0) then
         IF( CRR(I).EQ.'RR' ) THEN
            DMT = DMT * X(I)
         ELSE
            DMT = 0.D0
         ENDIF
        else
         IF( CRR(I).EQ.'RR' ) THEN
            DMT = DMT * X(I) * DCOS(K(I)*OM*TIME)
         ELSEIF( CRR(I).EQ.'RI' ) THEN
            DMT = -DMT * X(I) * DSIN(K(I)*OM*TIME)
         ELSE
            DMT = 0.D0
         ENDIF
        endif
         DMTOR = DMTOR-DMT
      enddo
   END function DMTOR

   !------------------------------------------------------------------------
   !  phi-gemittelte Stomlinien des Poloidalfeldes fuer phi = konstant:
   !            < F_phi > = r sin(theta) dtheta h (m = 0)
   !  oder des elektrischen Stromes:
   !            < F_phi > = r sin(theta) dtheta g (m = 0)
   !
   !     optimized for K = 0.
   FUNCTION DMPJ(X,whatToPlot,R,NTHETA,TIME)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*1 whatToPlot
      CHARACTER*2 CRR,whatToPlot

      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),whatToPlot(NM),CRR(NM)
      COMMON/PARI/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPARI/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD
      COMMON/AB/A(NAM),B(NAM),NAMC

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN DMPJ.'
         STOP
      ENDIF
      IF( NAM.NE.NAMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NAM IN DMPJ.'
         STOP
      ENDIF

      RI = ETA/(1.D0-ETA)
      RO = RI+1.D0
      DMPJ = 0.D0
      IF( whatToPlot.EQ.'MP' ) THEN
         NDOMIN = NDV+NDW+NDT+1
         NDOMAX = NDV+NDW+NDT+NDH
      ELSEIF( whatToPlot.EQ.'MJ' ) THEN
         NDOMIN = NDV+NDW+NDT+NDH+1
         NDOMAX = NDV+NDW+NDT+NDH+NDG
      ELSE
         WRITE(*,*) 'WRONG whatToPlot IN DMPJ, SHOULD BE MP OR MJ BUT IS: ',whatToPlot
         STOP
      ENDIF
      DO I = NDOMIN,NDOMAX
         IF( .NOT.( ( whatToPlot(I).EQ.'H' .AND. whatToPlot.EQ.'MP' ) .OR. &
                   ( whatToPlot(I).EQ.'G' .AND. whatToPlot.EQ.'MJ' )  ) ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN DMPJ, SHOULD BE H OR G BUT IS: ', whatToPlot(I)
            STOP
         ENDIF
         IF( M(I).NE.0 ) cycle
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         ENDIF
         DMP = EPSK*R * DBLE(L(I)*(L(I)+1)) * ( &
              PLMS(L(I)+1,M(I),NTHETA) / DSQRT( DBLE( (2*L(I)+1)*(2*L(I)+3) ) )  - &
              PLMS(L(I)-1,M(I),NTHETA) / DSQRT( DBLE( (2*L(I)+1)*(2*L(I)-1) ) )   )
         IF( whatToPlot(I).EQ.'H' ) THEN
            NR = NAB(L(I),N(I))
            IF( R.LE.RO ) THEN
               IF( A(NR).EQ.0.D0 .OR. B(NR).EQ.0.D0 ) THEN
                  WRITE(*,*) 'ALPHA AND BETA NOT CALCULATED.'
                  STOP
               ENDIF
               DMP = DMP*DCOS( A(NR)*R-B(NR) )
            ELSE
               IF( A(NR).EQ.0.D0 .OR. B(NR).EQ.0.D0 ) THEN
                  WRITE(*,*) 'ALPHA AND BETA NOT CALCULATED.'
                  STOP
               ENDIF
               DMP = DMP * (RO/R)**(L(I)+1) * DCOS( A(NR)*RO-B(NR) )
            ENDIF
         ELSEIF( whatToPlot(I).EQ.'G' ) THEN
            DMP = DMP*DSIN( N(I)*PI*(R-RI) )
         ENDIF
        if(K(I).EQ.0) then
         IF( CRR(I).EQ.'RR' ) THEN
            DMP = DMP * X(I)
         ELSE
            DMP = 0.D0
         ENDIF
        else
         IF( CRR(I).EQ.'RR' ) THEN
            DMP = DMP * X(I) * DCOS(K(I)*OM*TIME)
         ELSEIF( CRR(I).EQ.'RI' ) THEN
            DMP = -DMP * X(I) * DSIN(K(I)*OM*TIME)
         ELSE
            DMP = 0.D0
         ENDIF
        endif
         DMPJ = DMPJ+DMP
      enddo
   END function dmpj

   !------------------------------------------------------------------------
   ! Ueber Phi gemittelte Phi-Komponente des elektrischen Stromes:
   !            dtheta laplace h  (m = 0).
   !
   !     optimized for K = 0.
   FUNCTION DMC(X,R,NTHETA,TIME)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*1 whatToPlot
      CHARACTER*2 CRR

      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),whatToPlot(NM),CRR(NM)
      COMMON/PARI/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPARI/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD
      COMMON/AB/A(NAM),B(NAM),NAMC

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN DMC.'
         STOP
      ENDIF
      IF( NAM.NE.NAMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NAM IN DMC.'
         STOP
      ENDIF
      RI = ETA/(1.D0-ETA)
      DMC = 0.D0
      NDOMIN = NDV + NDW + NDT + 1
      NDOMAX = NDV + NDW + NDT + NDH
      DO I = NDOMIN,NDOMAX
         IF( whatToPlot(I).NE.'H' ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN DMC, SHOULD BE H BUT IS: ', whatToPlot(I)
            STOP
         ENDIF
         IF( M(I).NE.0 ) cycle
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         ENDIF
         DM = EPSK*DSQRT(DBLE(L(I)*(L(I)+1))) * PLMS(L(I),1,NTHETA)
         NR = NAB(L(I),N(I))
         IF( A(NR).EQ.0.D0 .OR. B(NR).EQ.0.D0 ) THEN
            WRITE(*,*) 'ALPHA AND BETA NOT CALCULATED.'
            STOP
         ENDIF
         DM = -DM*( ( A(NR)*A(NR)+DBLE(L(I)*(L(I)+1))/(R*R) )*DCOS(A(NR)*R-B(NR) ) + &
                    2*A(NR)/R * DSIN( A(NR)*R-B(NR) )  )

        if(K(I).EQ.0) then
           IF( CRR(I).EQ.'RR' ) THEN
              DM = DM * X(I)
           ELSE
              DM = 0.D0
           ENDIF
        else
           IF( CRR(I).EQ.'RR' ) THEN
              DM = DM * X(I) * DCOS(K(I)*OM*TIME)
           ELSEIF( CRR(I).EQ.'RI' ) THEN
              DM = -DM * X(I) * DSIN(K(I)*OM*TIME)
           ELSE
              DM = 0.D0
           ENDIF
        endif
         DMC = DMC-DM
      enddo
   END function dmc

   !------------------------------------------------------------------------
   !  THIS PROGRAM FINDS THE A'S AND B'S OF THE POLODIAL MAGNETIC
   !  FIELD TO FULLFILL THE BOUNDARY CONDITIONS:
   !  A(I)*TAN(A(I)*RO-B(I))-(L+1)/RO = 0  AND
   !  A(I)*TAN(A(I)*RI-B(I))+L/RI = 0 WITH A PRCISSION OF 1D-13.
   !  THE A'S AND B'S ARE STORED LINEARLY IN THE ARRAYS, NAB(L,N)
   !  DETERMINS THE POSITION IN THE ARRAY.
   !  NEEDS FUNCTIONS AMIN,NAB .
   SUBROUTINE ABG(ND,whatToPlot,LA,NA)
      IMPLICIT REAL*8(A-H,O-Y)
      CHARACTER(len=1):: whatToPlot(:)
      integer:: ND
      double precision:: LA(:),NA(:)
      double precision:: ri, ro

      COMMON/LOG/LCALC,LWRITE,LTR,LVAR,LDY,L6,L7,L8,L9,L10
      COMMON/PAR/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPAR/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD
      COMMON/NUM/RELE,EPS,ALPH,STOER,NITMAX,NJA

      COMMON/AB/A(NAM),B(NAM),NAMC
      COMMON/ABMIN/RIAB,ROAB,RELEAB,LAB
      COMMON/LNMAX/NLMAC,NL,LC(100),NMAXC(100),NMABC

      IF( NAM.NE.NAMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NAM IN ABG.'
         STOP
      ENDIF
      IF( NLMA.NE.NLMAC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NLMA IN ABG.'
         STOP
      ENDIF

      CALL CALCNMAX(ND,whatToPlot,LA,NA)

      RI = ETA/(1.D0-ETA)
      RO = RI+1.D0
      RIAB = RI
      ROAB = RO
      DAX = 1.D-3

      IF( DAX.LT.RELE*1.D3 ) THEN
         RELEAB = RELE*1.D-4
      ELSE
         RELEAB = RELE
      ENDIF
      RELEAB = DMAX1(RELEAB,EPS)

      IA = 1
      DO 1000 NI = 1,NL
         L = IABS(LC(NI))
         NMAX = NMAXC(NI)

         IF( NMAX.LE.0 ) THEN
            cycle
         ELSEIF( NMAX.GT.NAM ) THEN
            WRITE(*,*) 'TOO SMALL NAM IN DABG.'
            STOP
         ENDIF
         N = 1
         LAB = L
         IF( RI.EQ.0 ) THEN
            DO I = 0,2000
               IF( I.EQ.0 ) THEN
                  AXMIN = DAX
               ELSE
                  AXMIN = (I-0.5D0)*DPI+DAX
               ENDIF
               AXMAX = (I+0.5D0)*DPI-DAX
               IF( IA.GT.NAM ) THEN
                  WRITE(*,*) 'TOO SMALL DIMENSION NAM IN ABG.'
                  STOP
               ENDIF
               AGUESS = AXMAX
               LBT = 0
90             A(IA) = AMINB(AGUESS)
               IF( LBT.EQ.0 .AND. ( A(IA).LT.AXMIN .OR. A(IA).GT.AXMAX ) ) THEN
                  LBT = 1
                  AGUESS = AXMIN
                  GOTO 90
               ELSEIF( LBT.EQ.1 .AND. ( A(IA).LT.AXMIN .OR. A(IA).GT.AXMAX ) ) THEN
                  WRITE(*,*) 'WRONG ALPHA!!!!'
                  WRITE(*,'(X,'' IA,L,ALPHAMIN,ALPHA,ALPHAMAX: '',2I4,3D16.6)') IA,L,AXMIN,A(IA),AXMAX
                  STOP
               ENDIF
               B(IA) = 0.0D0
               IA = IA+1
               N = N+1
               IF(N.GT.NMAX) GOTO 1000
            enddo
         ELSE
            CD = DSQRT(L*(L+1)/RI/RO)
            AXMIN = DAX
            AXMAX = 0.5D0*DPI-DAX
            IF(AXMAX.GT.CD) THEN
               IF( IA.GT.NAM ) THEN
                  WRITE(*,*) 'TOO SMALL DIMENSION NAM IN ABG.'
                  STOP
               ENDIF
               AGUESS = AXMAX
               LBT = 0
190            A(IA) = AMIN(AGUESS)
               IF( LBT.EQ.0 .AND. ( A(IA).LT.AXMIN .OR. A(IA).GT.AXMAX ) ) THEN
                  LBT = 1
                  AGUESS = AXMIN
                  GOTO 190
               ELSEIF( LBT.EQ.1 .AND.  ( A(IA).LT.AXMIN .OR. A(IA).GT.AXMAX ) ) THEN
                  WRITE(*,*) 'WRONG ALPHA!!!!'
                  WRITE(*,'(X,'' IA,L,ALPHAMIN,ALPHA,ALPHAMAX: '',2I4,3D16.6)') IA,L,AXMIN,A(IA),AXMAX
                  STOP
               ENDIF
               B(IA) = A(IA)*RI+DATAN(L/A(IA)/RI)
               IA = IA+1
               N = N+1
               IF(N.GT.NMAX) GOTO 1000
            ENDIF
            DO I = 1,2000
               DAX = 1D-3
               AXMIN = (I-0.5D0)*DPI+DAX
               AXMAX = (I+0.5D0)*DPI-DAX
               IF(AXMIN.LT.CD .AND. AXMAX.GT.CD ) THEN
                  IF( IA.GT.NAM ) THEN
                     WRITE(*,*) 'TOO SMALL DIMENSION NAM IN ABG.'
                     STOP
                  ENDIF
                  AGUESS = AXMIN
                  LBT = 0
290               A(IA) = AMIN(AGUESS)
                  IF( LBT.EQ.0 .AND. ( A(IA).LT.AXMIN .OR. A(IA).GT.AXMAX ) ) THEN
                     LBT = 1
                     AGUESS = AXMAX
                     GOTO 290
                  ELSEIF( LBT.EQ.1 .AND. ( A(IA).LT.AXMIN .OR. A(IA).GT.AXMAX ) ) THEN
                     WRITE(*,*) 'WRONG ALPHA!!!!'
                     WRITE(*,'(X,'' IA,L,ALPHAMIN,ALPHA,ALPHAMAX: '',2I4,3D16.6)') IA,L,AXMIN,A(IA),AXMAX
                     STOP
                  ENDIF
                  B(IA) = A(IA)*RI+DATAN(L/A(IA)/RI)
                  IA = IA+1
                  N = N+1
                  IF(N.GT.NMAX) GOTO 1000
               ENDIF
150            CONTINUE
               IF( IA.GT.NAM ) THEN
                  WRITE(*,*) 'TOO SMALL DIMENSION NAM IN ABG.'
                  STOP
               ENDIF
               AGUESS = AXMAX
               LBT = 0
390            A(IA) = AMIN(AGUESS)
               IF( LBT.EQ.0 .AND. ( A(IA).LT.AXMIN .OR. A(IA).GT.AXMAX ) ) THEN
                  LBT = 1
                  AGUESS = AXMIN
                  GOTO 390
               ELSEIF( LBT.EQ.1 .AND. ( A(IA).LT.AXMIN .OR. A(IA).GT.AXMAX ) ) THEN
                  WRITE(*,*) 'WRONG ALPHA!!!!'
                  WRITE(*,'(X,'' IA,L,ALPHAMIN,ALPHA,ALPHAMAX: '',2I4,3D16.6)') IA,L,AXMIN,A(IA),AXMAX
                  STOP
               ENDIF
               B(IA) = A(IA)*RI+DATAN(L/A(IA)/RI)
               IA = IA+1
               N = N+1
               IF(N.GT.NMAX) exit
            enddo
         ENDIF
1000  CONTINUE
      DO I = 1,IA-1
         IF( I.GT.1 .AND. ( A(I).GT.A(I-1)-RELE .AND. A(I).LT.A(I-1)+RELE ) ) THEN
            WRITE(*,*) 'TWO ALPHAS EQUAL: ',A(I-1),A(I)
            STOP
         ENDIF
         DO J = 1,100
            B(I) = B(I)-DPI
            IF(B(I).LT.0.D0) THEN
               B(I) = B(I)+DPI
               exit
            ENDIF
         enddo
      enddo

      IA = IA-1
      WRITE(*,*) IA,' ALPHA AND BETA CALCULATED.'
      NMABC = IA
      DO I = 1,IA
         WRITE(*,'(2X,I4,2D14.6)') I,A(I),B(I)
      enddo
   END SUBROUTINE abg

   !------------------------------------------------------------------------
   !  FINDS THE MINIMUM FOR THE FUNCTION IN LINE 5 WITH A NEWTON METHOD.
   FUNCTION AMIN(AX)
      IMPLICIT REAL*8(A-H,O-Y)
      COMMON/ABMIN/RI,RO,RELE,L
      ICOUNT = 0

5     FA = DTAN(AX)-(L*RO+(L+1)*RI)*AX/(RI*RO*AX**2-L*(L+1))
      FAA = 1D0/DCOS(AX)**2-( (L*RO+(L+1)*RI)*(RI*RO*AX**2-L*(L+1)) - &
            AX*(L*RO+(L+1)*RI)*2*RI*RO*AX )/(RI*RO*AX**2-L*(L+1))**2
      IF(FAA.EQ.0) THEN
         AX = AX+RELE
         GOTO 5
      ENDIF
      DA = FA/FAA
      AOX = AX
      AX = AX-DA
      IF(DABS(1-DABS(AOX/AX)).LT.RELE) THEN
         AMIN = AX
         RETURN
      ENDIF
      ICOUNT = ICOUNT+1
      IF(ICOUNT.GT.100) THEN
         WRITE(*,*) 'NO ZERO FOUND IN DABG/AMIN.'
         STOP
      ENDIF
      GOTO 5
   END function amin

   !------------------------------------------------------------------------
   !   FINDS THE MINIMUM FOR THE FUNCTION IN LINE 5 WITH A NEWTON METHOD.
   FUNCTION AMINB(AX)
      IMPLICIT REAL*8(A-H,O-Y)
      COMMON/ABMIN/RI,RO,RELE,L
      do
         do
            FA = DTAN(AX*RO)-(L+1)/AX/RO
            FAA = RO/DCOS(AX*RO)**2+(L+1)/AX**2/RO
            IF(FAA.EQ.0) THEN
               AX = AX+RELE
            else
               exit
            ENDIF
         enddo
         DA = FA/FAA
         AOX = AX
         AX = AX-DA
         IF(DABS(1-DABS(AOX/AX)).LT.RELE) THEN
            AMINB = AX
            RETURN
         ENDIF
      enddo
   END function aminb

   !------------------------------------------------------------------------
   !  DETERMINS THE POSITION OF AN A OR B IN THE ARRAY A(I),B(I)
   !  DEPENDING ON L AND N.
   integer FUNCTION NAB(L,N)
      IMPLICIT REAL*8(A-H,O-Y)
      integer, intent(in):: l,n
      integer:: lmin, lmax, NI, ll

      COMMON/LOG/LCALC,LWRITE,LTR,LVAR,LDY,L6,L7,L8,L9,L10
      COMMON/PAR/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPAR/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD
      COMMON/LNMAX/NLMAC,NL,LC(NLMA),NMAXC(NLMA),NMABC

      IF( NLMA.NE.NLMAC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NLMA IN NAB.'
         STOP
      ENDIF

      LMIN = LC(1)
      LMAX = LC(NL)
      IF( L.GT.LMAX .OR. L.LT.LMIN ) THEN
         WRITE(*,*) 'WRONG L IN NAB',L,LMIN,LMAX
         STOP
      ENDIF
      NAB = 0
      DO NI = 1,NL
         LL = LC(NI)
         NMAX = NMAXC(NI)
         IF( LL.LT.L ) THEN
            IF( NMAX.GT.0 ) NAB = NAB+NMAX
         ELSEIF( LL.EQ.L ) THEN
            IF( NMAX.LT.N ) THEN
               WRITE(*,*) 'WRONG N IN NAB',N,NMAX
               STOP
            ELSE
               NAB = NAB+N
               RETURN
            ENDIF
         ENDIF
      enddo

      IF( N.GT.NMABC ) THEN
         WRITE(*,*) 'N LARGER THE CALCULATED NUMBER OF A,B IN NAB: ',N,NMABC
         STOP
      ENDIF
   END function nab

   !------------------------------------------------------------------------
   !-- CALCULATES THE MAXIMUM N FOR EACH L.
   !   THIS IS USED FOR CALCULATING THE RADIAL FUNCTION OF H.
   SUBROUTINE CALCNMAX(NK,whatToPlot,L,N)
      IMPLICIT REAL*8(A-H,O-Z)
      integer:: NK
      CHARACTER(len=1):: whatToPlot(:)
      double precision:: L(*),N(*)
      integer:: i, lold

      COMMON/LNMAX/NLMAC,NL,LC(NLMA),NMAX(NLMA),NMABC

      NLMAC = NLMA

      !-- BESTIMMMUNG VON NMAX FUER JEDES L , NOTWENDIG IN ABG:
      LOLD = 10000
      NL = 0
      DO I = 1,NK
         IF( whatToPlot(I).EQ.'H' ) THEN
            IF( L(I).NE.LOLD ) THEN
               NL = NL+1
               IF( NL.GT.NLMA ) THEN
                  WRITE(*,*) 'TOO SMALL DIMENSION NLMA IN CALCNMAX.'
                  STOP
               ENDIF
               LC(NL) = L(I)
               NMAX(NL) = N(I)
               LOLD = L(I)
            ELSEIF( L(I).EQ.LOLD .AND. N(I).GT.NMAX(NL) ) THEN
               NMAX(NL) = N(I)
            ENDIF
         ENDIF
      enddo
   END SUBROUTINE CALCNMAX

   !------------------------------------------------------------------------
   ! Phi -Komponente des Toroidalfeldes: - dtheta g
   FUNCTION DBT(X,R,PHI,NTHETA,TIME,driftRate)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*1 whatToPlot
      CHARACTER*2 CRR

      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),whatToPlot(NM),CRR(NM)
      COMMON/PARI/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPARI/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN DMPJ.'
         STOP
      ENDIF

      PPHI = PHI*PI/180.D0
      RI = ETA/(1.D0-ETA)
      DBT = 0.D0
      NDOMIN = NDV+NDW+NDT+NDH+1
      NDOMAX = NDV+NDW+NDT+NDH+NDG
      DO I = NDOMIN,NDOMAX
         IF( whatToPlot(I).NE.'G' ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN DBT, SHOULD BE G BUT IS: ', whatToPlot(I)
            STOP
         ENDIF
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         ENDIF
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         ENDIF
         DB = -0.5D0*EPSM*EPSK* ( &
              DSQRT(DBLE((L(I)-M(I)+1)*(L(I)+M(I)))) * PLMS(L(I),M(I)-1,NTHETA) -&
              DSQRT(DBLE((L(I)+M(I)+1)*(L(I)-M(I)))) * PLMS(L(I),M(I)+1,NTHETA) )
         DB = DB*DSIN( N(I)*PI*(R-RI) )
         IF( CRR(I).EQ.'RR' ) THEN
            DB = DB * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
         ELSEIF( CRR(I).EQ.'IR' ) THEN
            DB = -DB * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DCOS(K(I)*OM*TIME)
         ELSEIF( CRR(I).EQ.'RI' ) THEN
            DB = -DB * X(I) * DCOS( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
         ELSEIF( CRR(I).EQ.'II' ) THEN
            DB = DB * X(I) * DSIN( M(I)*(PPHI-driftRate*TIME) ) * DSIN(K(I)*OM*TIME)
         ENDIF
         DBT = DBT+DB
      enddo
   END function dbt

   !----------------------------------------------------------------
   SUBROUTINE READLA(STARTFILE,NUDSR,TIMER,X)
      IMPLICIT REAL*8(A-H,O-Z)
      integer:: NUDSR
      CHARACTER*1 whatToPlot,CFS
      CHARACTER*2 CRR
      CHARACTER*40 STARTFILE
      double precision, intent(in):: x(:)

      COMMON/LOG/LCALC,LWRITE,L3,LVAR,LDY,LT,L7,L8,L9,L10
      COMMON/PARI/RA,TA,PR,PM,ETA,C,OM,FTW,FTG,MF
      COMMON/NPARI/M0,NTV,NTH,LTV,LTH,KTV,KTH,LEV,LRB,LD
      COMMON/NUM/RELE,EPS,ALPH,STOER,NITMAX,NJA
      COMMON/DIM/NDV,NDW,NDT,NDH,NDG,ND
      COMMON/DIMS/NDVS,NDWS,NDTS,NDHS,NDGS,NDS
      COMMON/QNU/NMC,L(NM),M(NM),N(NM),K(NM),whatToPlot(NM),CRR(NM)
      COMMON/QNUS/NMSC,LS(NM),MS(NM),NS(NM),CFS(NM)
      COMMON/QNUVI/NUC,NUOM
      COMMON/TINT/DT,TIME0,NST,NTW,NDTW,NDTPW,NKWT,IKWT(NMT)

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN READP.'
         STOP
      ENDIF

      !-- IF LRT.EQ.0 WE ARE LOOKING FOR THE RIGHT NUMBER OF DATASET NUDSR,
      !-- IF LRT.EQ.1 WE ARE LOOKING FOR THE RIGHT TIME.
      OPEN(12,FILE=STARTFILE,STATUS='old')

      !-- COUNTER AND LOGICALS SET TO ZERO:
      NKZR=0

      !-- LST IS FORMAT PARAMETER ( L=1 FOR HIRSCHING FORMAT ) ,
      !   LCALCI TELLS HOW THE INPUTFILE HAS BEEN CALCULATED:
      READ(12,*) LST,LCALC

      !-- FUER LCALC.EQ.5 ODER LCALC.EQ.6 LIEGT ZEITINTEGRATION VOR ,
      !   D.H. VERSCHIEDENEN ZEITEN MUESSEN IM INPUTFILE DURCH
      !   READLA GESUCHT WERDEN. (LT=1), LTR=1 ZEIGT AN , DASS
      !   ZEITINTEGRATION VORLIEGT UND DIE ZEIT GELESEN WERDEN MUSS:
      IF( LCALC.EQ.5 .OR. LCALC.EQ.6 ) THEN
         LT=1
         LTR=1
      ELSE
         LT=0
         LTR=0
      ENDIF

      !-- READH READS THE HEADER OF THE INPUTFILE AND DETERMINS WETHER
      !   THE DATASET (dataSetNumber,TIME) IS THE RIGHT ONE (LDR=1):
      CALL READH(12,LTR,NUDSR,TIMER,dataSetNumber,TIME,LDR)

      !-- LOOKING FOR THE RIGHT DATASET:
      DO I=1,1000
         !----- READD READS FULL SET OF COEFFITIENTS:
         CALL READD(12,LDR,ND,X,whatToPlot,CRR,L,M,N,K, EVPM,EVPF,EVTM,DNU,EVTF,EMPM,EMPF,EMTM,EMTF)
         IF( LDR.EQ.1 ) exit
      enddo

      if(ND.GT.NM) then
        write(*,*) 'To small dimension NM in READLA'
        stop
      endif
      TA=TA**2

      LSX=1
      CALL SORTK(ND,LSX,X,whatToPlot,CRR,L,M,N,K,NUC,NUOM)
      CALL RDIM(ND,whatToPlot,CRR,L,M,N,K,CFS,LS,MS,NS, &
                       NDV,NDW,NDT,NDH,NDG,NDVS,NDWS,NDTS,NDHS,NDGS,NDS)
      CLOSE(12)
   END SUBROUTINE READLA
END PROGRAM LARA
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
