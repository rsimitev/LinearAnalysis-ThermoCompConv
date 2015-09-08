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
   use parameters
   implicit none
   double precision, parameter:: PI = 3.14159265358979D0
   integer, PARAMETER:: NM = 5500, NAM = 400, nPlotsMAX = 9, nSubPlotsMAX = 4, NPAM = nPlotsMAX*nSubPlotsMAX
   integer, PARAMETER:: NLMA = 100
   integer, PARAMETER:: NMX = 65, NMY = 128
   integer:: NMC
   integer:: NCPLOT, LR, NQ, NR, n, l
   integer:: drawPlotNum, drawFrame, drawHeader, drawTime
   logical:: countourParIsNumber
   CHARACTER*40 INPUTFILE,OUTPUTFILE
   CHARACTER*1 CFS
   CHARACTER*2 CRR
   INTEGER:: i, j
   double precision:: THETA(NMC)
   double precision:: dt

   double precision:: DX(NM), constantCoordinateValue(nPlotsMAX,nSubPlotsMAX),normRadiusMax(nPlotsMAX,nSubPlotsMAX)
   double precision:: TIME(nPlotsMAX), contourPar(nPlotsMAX,nSubPlotsMAX)
   integer:: nSubPlots(nPlotsMAX)
   character(len=1):: constantCoordinate(nPlotsMAX,nSubPlotsMAX), thisPlotconstantCoordinate(NPAM)
   character(len=2):: whatToPlot(nPlotsMAX,nSubPlotsMAX), thisPlotWhatToPlot(NPAM)
   character(len=2):: domain(nPlotsMAX), quadrant(nPlotsMAX,nSubPlotsMAX), thisPlotQuadrant(NPAM), XCP(NPAM)
   double precision:: zdo
   DIMENSION ZDP(NPAM),TIMEP(NPAM)
   DIMENSION XRICM(NPAM),XRMCM(NPAM),XRM(NPAM)

   NCPLOT = 0
   ZDO = 0.E0

   RELE = 1.D-9
   EPS = 1.D-13

   !-- INPUT:
   LR = 0
   NQ = 0
   NR = 0

   READ(*,*)
   READ(*,*) INPUTFILE,OUTPUTFILE

   OPEN(14,FILE = OUTPUTFILE,STATUS = 'unknown')
   write(14,*) 'Inputfile: ',INPUTFILE

   !-- SETTING OF CONSTANTS:
   LTR = 1
   NMC = NM
   NMSC = NM
   NAMC = NAM
   NLMAC = NLMA

   !-- READLA READS THE SET OF COEFFITIENTS TO BE PLOTTED ,
   !   IT IS NECESSARY TO CALL IS HERE TO GET PARAMETERS.
   TIMEO = TIME(1)
   write(*,'(A,A,I3,D9.2)') 'reading data from ',INPUTFILE,TIME(1),'...'
   CALL READLA(INPUTFILE,TIME(1),DX)
   write(*,*) '...done'
   RA = RAI
   TA = TAI
   PR = PRI
   PM = PMI
   ETA = ETAI
   C = CI
   OM = OMI
   MF = 0
   M0 = M0I
   NTH = NTHI
   LTH = LTHI
   KTV = KTVI
   KTH = KTHI
   LD = LDI
   LEV = LEVI
   LRB = LRBI

   !-- CALCULATION OF INNER AND OUTER RADIUS:
   RI = ETA/(1.D0-ETA)
   RO = 1.D0+RI
   XRI = DBLE(RI)
   XRO = DBLE(RO)

   !-- ABG CALCULATES THE ALPHAS AND BETAS IN THE RADIAL function
   !   OF THE POLOIDAL MAGNETIC FIELD:
   IF( LCALC.EQ.2 .OR. LCALC.EQ.4 .OR. LCALC.EQ.6 ) CALL ABG(ND,whatToPlot,L,N)

   XLRAND = 3.0D0
   XRRAND = 3.0D0
   NROWR = NR
   NROWQ = NQ/2
   IF( MOD(NQ,2).NE.0 ) NROWQ = NROWQ+1
   NROWR = NR
   NROWQ = NQ/2
   IF( MOD(NQ,2).NE.0 ) NROWQ = NROWQ+1
   NROW = NROWR+NROWQ

   XBR = 2*YHR

   YHPG = NROWR*YHR+NROWQ*YHQ+(NROW-1)*YINTER
   IF( NQ.GT.0 ) THEN
      XBPG = 2*XLQ+NCOL*XBQ+(NCOL-1)*XINTER
   ELSE
      XBPG = 2*XLR+XBR
   endif
   XAREA = XLRAND+XBPG+XRRAND
   YAREA = plotNumSpaceY+YHPG

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
         endif
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
         endif
         YHPLOT = YHQ
         XBPLOT = XBQ
         XLPLOT = XLQ
      endif
      XORIG = XLRAND+XLPLOT+(NSPA-1)*(XBPLOT+XINTER)
      YORIG = plotNumSpaceY+YHPG-NZEI*YHPLOT-(NZEI-1)*YINTER

      IF( domain(I).EQ.'QU' .OR. domain(I).EQ.'SP' .AND. NCOL.GT.1 ) THEN
         NQT = NQT+1
         IF( NSPA.EQ.1 .AND. NQT.EQ.NQ ) XLQ = XLQ+(XBQ+XINTER)/2
      endif
      DO J = 1,nSubPlots(I)
         NP = NP+1
         thisPlotQuadrant(NP)           = quadrant(I,J)
         thisPlotWhatToPlot(NP)         = whatToPlot(I,J)
         thisPlotconstantCoordinate(NP) = constantCoordinate(I,J)
         IF( constantCoordinate(I,J).EQ.'R' ) THEN
            XCP(NP) = XRI+constantCoordinateValue(I,J)
         ELSE
            XCP(NP) = constantCoordinateValue(I,J)
         endif
         ZDP(NP)   = contourPar(I,J)
         TIMEP(NP) = TIME(I)
      enddo
   enddo

   YTEXT = plotNumSpaceY-1.0E0

   !-- PLO FUEHRT DIE EINZELNEN SUBPLOTS AUS:

   CALL READLA(INPUTFILE,DX)
   CALL PLO(DX)

   CLOSE(14)

contains

   !------------------------------------------------------------------------
   !     calculates the field Z and makes one subplot.
   SUBROUTINE PLO()
      character*20 filez,filex,filey

      DIMENSION DX(*),THETA(NMY),XIDL(NMX,NMY),YIDL(NMX,NMY)
      DIMENSION Z(NMX,NMY),ZDS(4)

      write(*,*) 'computing the fields...'

      IF( constantCoordinate.EQ.'T' ) THEN
         NMTHETA = 1
         THETA(NMTHETA) = DBLE(constantCoordinateValue)
      ELSEIF( constantCoordinate.EQ.'P' ) THEN
         PHI = DBLE(constantCoordinateValue)
      ELSEIF( constantCoordinate.EQ.'R' ) THEN
         R = DBLE(constantCoordinateValue)
      endif
      XD = (XMAX-XMIN)/(NMX-1)
      YD = (YMAX-YMIN)/(NMY-1)
      IF( constantCoordinate.NE.'T' ) THEN
         NMTHETA = NMY
         DO I = 1, NMTHETA
            THETA(I) = DBLE(YMIN+(I-1)*YD)
         enddo
      endif

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
            endif
            !-------------------------------------------------------------------
            ! IDL- generate x and y coordinates.
            !-------------------------------------------------------------------
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
            Z(I,J) = flow_r(DX,R,PHI,THETA)
            IF( Z(I,J).GT.ZMAX ) ZMAX = Z(I,J)
            IF( Z(I,J).LT.ZMIN ) ZMIN = Z(I,J)
         enddo
      enddo

      range = MAX(ABS(ZMIN),ABS(ZMAX))
      ZNULL = 1.E-11*range
      ZNULLM = 1.E-11
      ZANULL = 1.E-13
      ZSCALE = 1.E0

      !-------------------------------------------------------------------------
      ! IDL
      !-------------------------------------------------------------------------
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
      endif

      close(21)
      close(22)
      close(23)
   end subroutine plo


   !------------------------------------------------------------------------
   ! Radiales Geschw.feld: U_r = L_2/r v
   !
   !     optimized for K = 0.
   function flow_r(X,R,PHI,THETA)
      IMPLICIT none
      double precision, intent(in):: x(:), R(:), theta(:), phi(:)

      PPHI = PHI*PI/180.D0
      RI = ETA/(1.D0-ETA)
      RF = 0.D0

      DO I = NDOMIN,NDOMAX
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         endif
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         endif

         RFT = EPSM*EPSK*L(I)*(L(I)+1) * PLMS(L(I),M(I),NTHETA) / R

         IF( whatToPlot(I).EQ.'V' ) THEN
            RFT = RFT*DSIN( N(I)*PI*(R-RI) )
         ELSEIF( whatToPlot(I).EQ.'H' ) THEN
            NR = NAB(L(I),N(I))
            IF( A(NR).EQ.0.D0 .OR. B(NR).EQ.0.D0 ) THEN
               WRITE(*,*) 'ALPHA AND BETA NOT CALCULATED.'
               STOP
            endif
            RFT = RFT*DCOS( A(NR)*R-B(NR) )
         endif

         if(K(I).EQ.0) then
            IF( CRR(I).EQ.'RR' ) THEN
               RFT = RFT * X(I) * DCOS( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               RFT = -RFT * X(I) * DSIN( M(I)*PPHI )
            ELSE
               RFT = 0.D0
            endif
         else
            IF( CRR(I).EQ.'RR' ) THEN
               RFT = RFT * X(I) * DCOS( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               RFT = -RFT * X(I) * DSIN( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               RFT = 0
            ELSEIF( CRR(I).EQ.'II' ) THEN
               RFT = 0
            endif
         endif
         flow_r = flow_r + RFT
      enddo
   end function

   !------------------------------------------------------------------------
   !   Temperaturfeld Theta ( =  Abweichung vom Grundzust.)
   !   optimized for K = 0.
   function TEMP(X,whatToPlot,R,PHI,NTHETA)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*2 CRR,whatToPlot(:)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN TEMP.'
         STOP
      endif

      TEMP = 0.D0
      PPHI = PHI*PI/180.D0
      RI = ETA/(1.D0-ETA)
      IF( whatToPlot.EQ.'TE' ) THEN
         NDOMIN = 1+NDV+NDW
         NDOMAX = NDV+NDW+NDT
      ELSE
         WRITE(*,*) 'WRONG whatToPlot IN TEMP, SHOULD BE TE BUT IS: ',whatToPlot
         STOP
      endif

      DO I = NDOMIN, NDOMAX
         IF( whatToPlot(I).NE.'T' ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN TEMP, SHOULD BE T BUT IS: ', whatToPlot(I)
            STOP
         endif
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         endif
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         endif
         TEM = EPSM*EPSK*PLMS(L(I),M(I),NTHETA)*DSIN( N(I)*PI*(R-RI) )

         if(K(I).EQ.0) then
            IF( CRR(I).EQ.'RR' ) THEN
               TEM = TEM * X(I) * DCOS( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               TEM = -TEM * X(I) * DSIN( M(I)*PPHI )
            ELSE
               TEM = 0.D0
            endif
         else
            IF( CRR(I).EQ.'RR' ) THEN
               TEM = TEM * X(I) * DCOS( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               TEM = -TEM * X(I) * DSIN( M(I)*PPHI )
            endif
         endif
         TEMP = TEMP+TEM
      enddo
   end function temp

   !------------------------------------------------------------------------
   !   Stromfunktion fuer theta = konstant:
   !      F_theta = r dphi v             (Busse: r/sin(th) d/dphi v )
   !   Fuer den elektrischen Strom:
   !              F_theta = r dphi g
   !
   !     optimized for K = 0.
   function FT(X,whatToPlot,R,PHI,NTHETA)
      IMPLICIT REAL*8(A-H,O-Z)
      double precision, intent(in):: x(:)
      character(len=2), intent(in):: whatToPlot(:)
      integer, intent(in):: NTHETA
      CHARACTER*2 CRR

      IF( NAM.NE.NAMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NAM IN FT.'
         STOP
      endif
      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN FT.'
         STOP
      endif

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
      endif
      DO  I = NDOMIN,NDOMAX
         IF( .NOT.( ( whatToPlot(I).EQ.'V' .AND. whatToPlot.EQ.'VS' ) .OR. &
                    ( whatToPlot(I).EQ.'H' .AND. whatToPlot.EQ.'BS' ) .OR. &
                    ( whatToPlot(I).EQ.'G' .AND. whatToPlot.EQ.'JS' )  )  ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN FT, SHOULD BE V OR H OR G BUT IS: ',whatToPlot(I)
            STOP
         endif
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         endif
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         endif
         FTT = EPSM*EPSK*M(I)*PLMS(L(I),M(I),NTHETA)*R
         IF( whatToPlot(I).EQ.'V' .OR. whatToPlot(I).EQ.'G' ) THEN
            FTT = FTT*DSIN( N(I)*PI*(R-RI) )
         ELSEIF( whatToPlot(I).EQ.'H' ) THEN
            NR = NAB(L(I),N(I))
            IF( R.LE.RO ) THEN
               IF( A(NR).EQ.0.D0 .OR. B(NR).EQ.0.D0 ) THEN
                  WRITE(*,*) 'ALPHA AND BETA NOT CALCULATED.'
                  STOP
               endif
               FTT = -FTT*DCOS( A(NR)*R-B(NR) )
            ELSE
               IF( A(NR).EQ.0.D0 .OR. B(NR).EQ.0.D0 ) THEN
                  WRITE(*,*) 'ALPHA AND BETA NOT CALCULATED.'
                  STOP
               endif
               FTT = -FTT * (RO/R)**(L(I)+1) * DCOS( A(NR)*RO-B(NR) )
            endif
         endif

         if(K(I).EQ.0) then
            IF( CRR(I).EQ.'RR' ) THEN
               FTT = -FTT * X(I) * DSIN( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               FTT = -FTT * X(I) * DCOS( M(I)*PPHI )
            ELSE
               FTT = 0.D0
            endif
         else
            IF( CRR(I).EQ.'RR' ) THEN
               FTT = -FTT * X(I) * DSIN( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               FTT = -FTT * X(I) * DCOS( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'RI' ) THEN
            ELSEIF( CRR(I).EQ.'II' ) THEN
            endif
         endif
         FT = FT-FTT
      enddo
   end function ft

   !------------------------------------------------------------------------
   ! Stromfunktion fuer phi = konstant:
   !              F_phi = r sin(theta) dtheta v  (like Busse)
   ! Fuer den elektrischen Strom:
   !              F_phi = r sin(theta) dtheta g
   !
   !     optimized for K = 0.
   function FP(X,whatToPlot,R,PHI,NTHETA)
      IMPLICIT REAL*8(A-H,O-Z)
      double precision, intent(in):: x(:)
      integer, intent(in):: NTHETA
      CHARACTER*2:: CRR,whatToPlot(:)

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN FP.'
         STOP
      endif
      IF( NAM.NE.NAMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NAM IN FP.'
         STOP
      endif

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
      endif
      DO I = NDOMIN,NDOMAX
         IF( .NOT.( ( whatToPlot(I).EQ.'V' .AND. whatToPlot.EQ.'VS' ) .OR.&
                   ( whatToPlot(I).EQ.'H' .AND. whatToPlot.EQ.'BS' ) .OR. &
                   ( whatToPlot(I).EQ.'G' .AND. whatToPlot.EQ.'JS' )  )  ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN FP, SHOULD BE V OR H OR G BUT IS: ',whatToPlot(I)
            STOP
         endif
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         endif
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         endif
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
               endif
               FPT = FPT*DCOS( A(NR)*R-B(NR) )
            ELSE
               IF( A(NR).EQ.0.D0 .OR. B(NR).EQ.0.D0 ) THEN
                  WRITE(*,*) 'ALPHA AND BETA NOT CALCULATED.'
                  STOP
               endif
               FPT = FPT * (RO/R)**(L(I)+1) * DCOS( A(NR)*RO-B(NR) )
            endif
         endif

         if(K(I).EQ.0) then
            IF( CRR(I).EQ.'RR' ) THEN
               FPT = FPT * X(I) * DCOS( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               FPT = -FPT * X(I) * DSIN( M(I)*PPHI )
            ELSE
               FPT = 0.D0
            endif
         else
            IF( CRR(I).EQ.'RR' ) THEN
               FPT = FPT * X(I) * DCOS( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               FPT = -FPT * X(I) * DSIN( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               FPT = 0
            ELSEIF( CRR(I).EQ.'II' ) THEN
               FPT = 0
            endif
         endif
         FP = FP+FPT
      enddo
   end function FP

   !------------------------------------------------------------------------
   !   Stromfunktion fuer r = konstant:
   !                     F_r = w      (like Busse, Hirsching: rw )
   !   Stromfunktion fuer r = konstant of the electric currents:
   !                     F_r = - laplace h
   !
   !     optimized for K = 0.
   function FR(X,whatToPlot,R,PHI,NTHETA)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*2 CRR,whatToPlot(:)

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN FR.'
         STOP
      endif

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
      endif
      DO I = NDOMIN,NDOMAX
         IF(  .NOT.( ( whatToPlot.EQ.'VS' .AND. whatToPlot(I).EQ.'W' ) .OR.&
                     ( whatToPlot.EQ.'BS' .AND. whatToPlot(I).EQ.'G' ) .OR.&
                     ( whatToPlot.EQ.'JS' .AND. whatToPlot(I).EQ.'H' )  )  ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN FR, SHOULD BE W OR G OR H BUT IS: ',whatToPlot(I)
            STOP
         endif
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         endif
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         endif
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
            endif
            FRT = FRT*( ( A(NR)*A(NR)+DBLE(L(I)*(L(I)+1))/(R*R) ) * DCOS( A(NR)*R-B(NR) ) + 2*A(NR)/R * DSIN( A(NR)*R-B(NR) )  )
         endif
         if(K(I).EQ.0) then
            IF( CRR(I).EQ.'RR' ) THEN
               FRT = FRT * X(I) * DCOS( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               FRT = -FRT * X(I) * DSIN( M(I)*PPHI )
            ELSE
               FRT = 0.D0
            endif
         else
            IF( CRR(I).EQ.'RR' ) THEN
               FRT = FRT * X(I) * DCOS( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               FRT = -FRT * X(I) * DSIN( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               FRT = 0
            ELSEIF( CRR(I).EQ.'II' ) THEN
               FRT = 0
            endif
         endif
         FR = FR+FRT
      enddo
   end function fr


   !------------------------------------------------------------------------
   !     temperature field Theta + Ts
   !     optimized for K = 0.
   function TT(X,whatToPlot,R,PHI,NTHETA)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*2 CRR,whatToPlot(:)

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN TT.'
         STOP
      endif

      PPHI = PHI*PI/180.D0
      RI = ETA/(1.D0-ETA)
      T = 0.D0
      IF( whatToPlot.EQ.'TT' ) THEN
         NDOMIN = NDV+NDW+1
         NDOMAX = NDV+NDW+NDT
      ELSE
         WRITE(*,*) 'WRONG whatToPlot IN TT, SHOULD TT BUT IS: ',whatToPlot
         STOP
      endif
      DO I = NDOMIN,NDOMAX
         IF( whatToPlot(I).NE.'T' ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN T, SHOULD BE T BUT IS: ', whatToPlot(I)
            STOP
         endif
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         endif
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         endif
         TT = EPSM*EPSK*PLMS(L(I),M(I),NTHETA)
         TT = TT*DSIN( N(I)*PI*(R-RI) )

         IF(K(I).EQ.0) THEN
          IF( CRR(I).EQ.'RR' ) THEN
            TT = TT * X(I) * DCOS( M(I)*PPHI )
          ELSEIF( CRR(I).EQ.'IR' ) THEN
            TT = -TT * X(I) * DSIN( M(I)*PPHI )
          ELSEIF( CRR(I).EQ.'RI' ) THEN
            TT = 0.0D0
          ELSEIF( CRR(I).EQ.'II' ) THEN
            TT = 0.0D0
          endif
         ELSE
          IF( CRR(I).EQ.'RR' ) THEN
            TT = TT * X(I) * DCOS( M(I)*PPHI )
          ELSEIF( CRR(I).EQ.'IR' ) THEN
            TT = -TT * X(I) * DSIN( M(I)*PPHI )
          ELSEIF( CRR(I).EQ.'RI' ) THEN
            TT = 0
          ELSEIF( CRR(I).EQ.'II' ) THEN
            TT = 0
          endif
         endif
         T = T+TT
      enddo

      !  add basic temperature field Ts:
      T = T - R * R / ( 2.D0 * PR )
      TT = T
   end

   !------------------------------------------------------------------------
   !   local Nusselt number NU(r = ri)
   !   optimized for K = 0.
   function localNusselt(X,whatToPlot,R,PHI,NTHETA)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      double precision:: localNusselt
      CHARACTER*2 CRR,whatToPlot(:)

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN localNusselt.'
         STOP
      endif

      localNusselt = 0.D0
      PPHI = PHI*PI/180.D0
      RI = ETA/(1.D0-ETA)
      IF( whatToPlot.EQ.'NU' ) THEN
         NDOMIN = 1+NDV+NDW
         NDOMAX = NDV+NDW+NDT
      ELSE
         WRITE(*,*) 'WRONG whatToPlot IN localNusselt, SHOULD BE NU BUT IS: ',whatToPlot
         STOP
      endif

      DO I = NDOMIN,NDOMAX
         IF( whatToPlot(I).NE.'T' ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN TEMP, SHOULD BE T BUT IS: ', whatToPlot(I)
            STOP
         endif
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         endif
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         endif
         localNusseltT = EPSM*EPSK*PLMS(L(I),M(I),NTHETA)*DBLE(N(I))*PI

         if(K(I).EQ.0) then
            IF( CRR(I).EQ.'RR' ) THEN
               localNusseltT = localNusseltT * X(I) * DCOS( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               localNusseltT = -localNusseltT * X(I) * DSIN( M(I)*PPHI )
            ELSE
               localNusseltT = 0.D0
            endif
         else
            IF( CRR(I).EQ.'RR' ) THEN
               localNusseltT = localNusseltT * X(I) * DCOS( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               localNusseltT = -localNusseltT * X(I) * DSIN( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               localNusseltT = 0
            ELSEIF( CRR(I).EQ.'II' ) THEN
               localNusseltT = 0
            endif
         endif
         localNusselt = localNusselt+localNusseltT
      enddo
      localNusselt = 1.D0 - PR/RI*localNusselt
   end function localNusselt

   !------------------------------------------------------------------------
   !   Zonaler Fluss = gemittelte phi-Komponente der Geschwindigkeit:
   !          < u_phi > = - dtheta w   (m = 0)
   !
   !     optimized for K = 0.
   function flow_p_zonal(X,R,NTHETA)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*2:: whatToPlot(:)
      CHARACTER*2:: CRR

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN flow_p_zonal.'
         STOP
      endif

      flow_p_zonal = 0.D0
      RI = ETA/(1.D0-ETA)
      NDOMIN = 1+NDV
      NDOMAX = NDV+NDW

      DO I = NDOMIN,NDOMAX
         IF( whatToPlot(I).NE.'W' ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN flow_p_zonal, SHOULD BE W BUT IS: ', whatToPlot(I)
            STOP
         endif

         IF( M(I).NE.0 ) cycle

         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         endif
         ZON = EPSK*DSQRT(DBLE(L(I)*(L(I)+1))) * PLMS(L(I),1,NTHETA) * R * DCOS( (N(I)-1)*PI*(R-RI) )

         if(K(I).EQ.0) then
            IF( CRR(I).EQ.'RR' ) THEN
               ZON = ZON * X(I)
            ELSE
               ZON = 0.D0
            endif
         else
            IF( CRR(I).EQ.'RR' ) THEN
               ZON = ZON * X(I)
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               ZON = 0
            ELSE
               ZON = 0.D0
            endif
         endif
         flow_p_zonal = flow_p_zonal+ZON
      enddo
   end function

   !------------------------------------------------------------------------
   !     Uphi = 1/(r*sinphi) d^2/drdph rv - d/dth w
   !
   !     optimized for K = 0.
   function flow_p(X,whatToPlot,R,PHI,NTHETA)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*2:: CRR,whatToPlot(:)
      DIMENSION:: THETA(NMY)

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN UP.'
         STOP
      endif

      THETAR = PI*THETA(NTHETA)/180.D0
      SINTH = DSIN(THETAR)

      PPHI = PHI*PI/180.D0
      RI = ETA/(1.D0-ETA)
      UP = 0.D0

      IF( whatToPlot.NE.'UP' ) THEN
        WRITE(*,*) 'WRONG whatToPlot IN UP, SHOULD BE UP BUT IS: ',whatToPlot
        STOP
      endif

      NDOMIN = NDV+1
      NDOMAX = NDV+NDW
      !------- toroidal part: --------------------------
      DO I = NDOMIN,NDOMAX
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         endif
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         endif

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
         endif

         IF( whatToPlot(I).EQ.'W' ) THEN
            UPT = UPT*R*DCOS( (N(I)-1)*PI*(R-RI) )
         ELSEIF( whatToPlot(I).EQ.'G' ) THEN
            UPT = UPT*DSIN( N(I)*PI*(R-RI) )
         endif
         IF(K(I).EQ.0) THEN
            IF( CRR(I).EQ.'RR' ) THEN
               UPT = UPT * X(I) * DCOS( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               UPT = -UPT * X(I) * DSIN( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               UPT = 0.0D0
            ELSEIF( CRR(I).EQ.'II' ) THEN
               UPT = 0.0D0
            endif
         ELSE
            IF( CRR(I).EQ.'RR' ) THEN
               UPT = UPT * X(I) * DCOS( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'IR' ) THEN
               UPT = -UPT * X(I) * DSIN( M(I)*PPHI )
            ELSEIF( CRR(I).EQ.'RI' ) THEN
               UPT = 0
            ELSEIF( CRR(I).EQ.'II' ) THEN
               UPT = 0
            endif
         endif

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
         endif
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         endif

         !------- 1/(rsinth) d^2/drdphi (rv) -------------
         IF( SINTH.EQ.0.D0 .OR. M(I).EQ.0 ) THEN
          UPT = 0.D0
         ELSE
            UPT = EPSM*EPSK * M(I) * PLMS(L(I),M(I),NTHETA)  / (R*SINTH)

            UPT = UPT*(R*N(I)*PI*DCOS(N(I)*PI*(R-RI))+DSIN(N(I)*PI*(R-RI)))

            IF(K(I).EQ.0) THEN
               IF( CRR(I).EQ.'RR' ) THEN
                  UPT = -UPT * X(I) * DSIN( M(I)*PPHI )
               ELSEIF( CRR(I).EQ.'IR' ) THEN
                  UPT = -UPT * X(I) * DCOS( M(I)*PPHI )
               ELSE
                  UPT = 0.D0
               endif
            ELSE
               IF( CRR(I).EQ.'RR' ) THEN
                  UPT = -UPT * X(I) * DSIN( M(I)*PPHI )
               ELSEIF( CRR(I).EQ.'IR' ) THEN
                  UPT = -UPT * X(I) * DCOS( M(I)*PPHI )
               ELSEIF( CRR(I).EQ.'RI' ) THEN
                  UPT = 0
               ELSEIF( CRR(I).EQ.'II' ) THEN
                  UPT = 0
               endif
            endif
         endif

         UP = UP+UPT
      enddo
   end function


   !------------------------------------------------------------------------
   !  phi-gemittelte Stomlinien des Poloidalfeldes fuer phi = konstant:
   !            < F_phi > = r sin(theta) dtheta h (m = 0)
   !  oder des elektrischen Stromes:
   !            < F_phi > = r sin(theta) dtheta g (m = 0)
   !
   !     optimized for K = 0.
   function DMPJ(X,whatToPlot,R,NTHETA)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*2:: CRR,whatToPlot(:)

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN DMPJ.'
         STOP
      endif
      IF( NAM.NE.NAMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NAM IN DMPJ.'
         STOP
      endif

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
      endif
      DO I = NDOMIN,NDOMAX
         IF( .NOT.( ( whatToPlot(I).EQ.'H' .AND. whatToPlot.EQ.'MP' ) .OR. &
                   ( whatToPlot(I).EQ.'G' .AND. whatToPlot.EQ.'MJ' )  ) ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN DMPJ, SHOULD BE H OR G BUT IS: ', whatToPlot(I)
            STOP
         endif
         IF( M(I).NE.0 ) cycle
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         endif
         DMP = EPSK*R * DBLE(L(I)*(L(I)+1)) * ( &
              PLMS(L(I)+1,M(I),NTHETA) / DSQRT( DBLE( (2*L(I)+1)*(2*L(I)+3) ) )  - &
              PLMS(L(I)-1,M(I),NTHETA) / DSQRT( DBLE( (2*L(I)+1)*(2*L(I)-1) ) )   )
         IF( whatToPlot(I).EQ.'H' ) THEN
            NR = NAB(L(I),N(I))
            IF( R.LE.RO ) THEN
               IF( A(NR).EQ.0.D0 .OR. B(NR).EQ.0.D0 ) THEN
                  WRITE(*,*) 'ALPHA AND BETA NOT CALCULATED.'
                  STOP
               endif
               DMP = DMP*DCOS( A(NR)*R-B(NR) )
            ELSE
               IF( A(NR).EQ.0.D0 .OR. B(NR).EQ.0.D0 ) THEN
                  WRITE(*,*) 'ALPHA AND BETA NOT CALCULATED.'
                  STOP
               endif
               DMP = DMP * (RO/R)**(L(I)+1) * DCOS( A(NR)*RO-B(NR) )
            endif
         ELSEIF( whatToPlot(I).EQ.'G' ) THEN
            DMP = DMP*DSIN( N(I)*PI*(R-RI) )
         endif
        if(K(I).EQ.0) then
         IF( CRR(I).EQ.'RR' ) THEN
            DMP = DMP * X(I)
         ELSE
            DMP = 0.D0
         endif
        else
         IF( CRR(I).EQ.'RR' ) THEN
            DMP = DMP * X(I)
         ELSEIF( CRR(I).EQ.'RI' ) THEN
            DMP = 0
         ELSE
            DMP = 0.D0
         endif
        endif
         DMPJ = DMPJ+DMP
      enddo
   end function dmpj

   !------------------------------------------------------------------------
   ! Ueber Phi gemittelte Phi-Komponente des elektrischen Stromes:
   !            dtheta laplace h  (m = 0).
   !
   !     optimized for K = 0.
   function DMC(X,R,NTHETA)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in):: NTHETA
      double precision, intent(in):: x(:)
      CHARACTER*1:: whatToPlot(:)
      CHARACTER*2:: CRR

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN DMC.'
         STOP
      endif
      IF( NAM.NE.NAMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NAM IN DMC.'
         STOP
      endif
      RI = ETA/(1.D0-ETA)
      DMC = 0.D0
      NDOMIN = NDV + NDW + NDT + 1
      NDOMAX = NDV + NDW + NDT + NDH
      DO I = NDOMIN,NDOMAX
         IF( whatToPlot(I).NE.'H' ) THEN
            WRITE(*,*) 'WRONG whatToPlot IN DMC, SHOULD BE H BUT IS: ', whatToPlot(I)
            STOP
         endif
         IF( M(I).NE.0 ) cycle
         IF( K(I).EQ.0 ) THEN
            EPSK = 1.D0
         ELSE
            EPSK = 2.D0
         endif
         DM = EPSK*DSQRT(DBLE(L(I)*(L(I)+1))) * PLMS(L(I),1,NTHETA)
         NR = NAB(L(I),N(I))
         IF( A(NR).EQ.0.D0 .OR. B(NR).EQ.0.D0 ) THEN
            WRITE(*,*) 'ALPHA AND BETA NOT CALCULATED.'
            STOP
         endif
         DM = -DM*( ( A(NR)*A(NR)+DBLE(L(I)*(L(I)+1))/(R*R) )*DCOS(A(NR)*R-B(NR) ) + &
                    2*A(NR)/R * DSIN( A(NR)*R-B(NR) )  )

        if(K(I).EQ.0) then
           IF( CRR(I).EQ.'RR' ) THEN
              DM = DM * X(I)
           ELSE
              DM = 0.D0
           endif
        else
           IF( CRR(I).EQ.'RR' ) THEN
              DM = DM * X(I)
           ELSEIF( CRR(I).EQ.'RI' ) THEN
              DM = 0
           ELSE
              DM = 0.D0
           endif
        endif
         DMC = DMC-DM
      enddo
   end function dmc

   !------------------------------------------------------------------------
   !  THIS PROGRAM FINDS THE A'S AND B'S OF THE POLODIAL MAGNETIC
   !  FIELD TO FULLFILL THE BOUNDARY CONDITIONS:
   !  A(I)*TAN(A(I)*RO-B(I))-(L+1)/RO = 0  AND
   !  A(I)*TAN(A(I)*RI-B(I))+L/RI = 0 WITH A PRCISSION OF 1D-13.
   !  THE A'S AND B'S ARE STORED LINEARLY IN THE ARRAYS, NAB(L,N)
   !  DETERMINS THE POSITION IN THE ARRAY.
   !  NEEDS functionS AMIN,NAB .
   SUBROUTINE ABG(ND,whatToPlot,LA,NA)
      IMPLICIT REAL*8(A-H,O-Y)
      CHARACTER(len=1):: whatToPlot(:)
      integer:: ND
      double precision:: LA(:),NA(:)
      double precision:: ri, ro

      IF( NAM.NE.NAMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NAM IN ABG.'
         STOP
      endif
      IF( NLMA.NE.NLMAC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NLMA IN ABG.'
         STOP
      endif

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
      endif
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
         endif
         N = 1
         LAB = L
         IF( RI.EQ.0 ) THEN
            DO I = 0,2000
               IF( I.EQ.0 ) THEN
                  AXMIN = DAX
               ELSE
                  AXMIN = (I-0.5D0)*DPI+DAX
               endif
               AXMAX = (I+0.5D0)*DPI-DAX
               IF( IA.GT.NAM ) THEN
                  WRITE(*,*) 'TOO SMALL DIMENSION NAM IN ABG.'
                  STOP
               endif
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
               endif
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
               endif
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
               endif
               B(IA) = A(IA)*RI+DATAN(L/A(IA)/RI)
               IA = IA+1
               N = N+1
               IF(N.GT.NMAX) GOTO 1000
            endif
            DO I = 1,2000
               DAX = 1D-3
               AXMIN = (I-0.5D0)*DPI+DAX
               AXMAX = (I+0.5D0)*DPI-DAX
               IF(AXMIN.LT.CD .AND. AXMAX.GT.CD ) THEN
                  IF( IA.GT.NAM ) THEN
                     WRITE(*,*) 'TOO SMALL DIMENSION NAM IN ABG.'
                     STOP
                  endif
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
                  endif
                  B(IA) = A(IA)*RI+DATAN(L/A(IA)/RI)
                  IA = IA+1
                  N = N+1
                  IF(N.GT.NMAX) GOTO 1000
               endif
150            CONTINUE
               IF( IA.GT.NAM ) THEN
                  WRITE(*,*) 'TOO SMALL DIMENSION NAM IN ABG.'
                  STOP
               endif
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
               endif
               B(IA) = A(IA)*RI+DATAN(L/A(IA)/RI)
               IA = IA+1
               N = N+1
               IF(N.GT.NMAX) exit
            enddo
         endif
1000  CONTINUE
      DO I = 1,IA-1
         IF( I.GT.1 .AND. ( A(I).GT.A(I-1)-RELE .AND. A(I).LT.A(I-1)+RELE ) ) THEN
            WRITE(*,*) 'TWO ALPHAS EQUAL: ',A(I-1),A(I)
            STOP
         endif
         DO J = 1,100
            B(I) = B(I)-DPI
            IF(B(I).LT.0.D0) THEN
               B(I) = B(I)+DPI
               exit
            endif
         enddo
      enddo

      IA = IA-1
      WRITE(*,*) IA,' ALPHA AND BETA CALCULATED.'
      NMABC = IA
      DO I = 1,IA
         WRITE(*,'(2X,I4,2D14.6)') I,A(I),B(I)
      enddo
   end SUBROUTINE abg

   !------------------------------------------------------------------------
   !  FINDS THE MINIMUM FOR THE function IN LINE 5 WITH A NEWTON METHOD.
   function AMIN(AX)
      IMPLICIT REAL*8(A-H,O-Y)
      COMMON/ABMIN/RI,RO,RELE,L
      ICOUNT = 0

5     FA = DTAN(AX)-(L*RO+(L+1)*RI)*AX/(RI*RO*AX**2-L*(L+1))
      FAA = 1D0/DCOS(AX)**2-( (L*RO+(L+1)*RI)*(RI*RO*AX**2-L*(L+1)) - &
            AX*(L*RO+(L+1)*RI)*2*RI*RO*AX )/(RI*RO*AX**2-L*(L+1))**2
      IF(FAA.EQ.0) THEN
         AX = AX+RELE
         GOTO 5
      endif
      DA = FA/FAA
      AOX = AX
      AX = AX-DA
      IF(DABS(1-DABS(AOX/AX)).LT.RELE) THEN
         AMIN = AX
         RETURN
      endif
      ICOUNT = ICOUNT+1
      IF(ICOUNT.GT.100) THEN
         WRITE(*,*) 'NO ZERO FOUND IN DABG/AMIN.'
         STOP
      endif
      GOTO 5
   end function amin

   !------------------------------------------------------------------------
   !   FINDS THE MINIMUM FOR THE function IN LINE 5 WITH A NEWTON METHOD.
   function AMINB(AX)
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
            endif
         enddo
         DA = FA/FAA
         AOX = AX
         AX = AX-DA
         IF(DABS(1-DABS(AOX/AX)).LT.RELE) THEN
            AMINB = AX
            RETURN
         endif
      enddo
   end function aminb

   !------------------------------------------------------------------------
   !  DETERMINS THE POSITION OF AN A OR B IN THE ARRAY A(I),B(I)
   !  DEPendING ON L AND N.
   integer function NAB(L,N)
      IMPLICIT REAL*8(A-H,O-Y)
      integer, intent(in):: l,n
      integer:: lmin, lmax, NI, ll

      IF( NLMA.NE.NLMAC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NLMA IN NAB.'
         STOP
      endif

      LMIN = LC(1)
      LMAX = LC(NL)
      IF( L.GT.LMAX .OR. L.LT.LMIN ) THEN
         WRITE(*,*) 'WRONG L IN NAB',L,LMIN,LMAX
         STOP
      endif
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
            endif
         endif
      enddo

      IF( N.GT.NMABC ) THEN
         WRITE(*,*) 'N LARGER THE CALCULATED NUMBER OF A,B IN NAB: ',N,NMABC
         STOP
      endif
   end function nab

   !------------------------------------------------------------------------
   !-- CALCULATES THE MAXIMUM N FOR EACH L.
   !   THIS IS USED FOR CALCULATING THE RADIAL function OF H.
   SUBROUTINE CALCNMAX(NK,whatToPlot,L,N)
      IMPLICIT REAL*8(A-H,O-Z)
      integer:: NK
      CHARACTER(len=1):: whatToPlot(:)
      double precision:: L(*),N(*)
      integer:: i, lold

      NLMAC = NLMA

      !-- BESTIMMMUNG VON NMAX FUER JEDES L , NOTWendIG IN ABG:
      LOLD = 10000
      NL = 0
      DO I = 1,NK
         IF( whatToPlot(I).EQ.'H' ) THEN
            IF( L(I).NE.LOLD ) THEN
               NL = NL+1
               IF( NL.GT.NLMA ) THEN
                  WRITE(*,*) 'TOO SMALL DIMENSION NLMA IN CALCNMAX.'
                  STOP
               endif
               LC(NL) = L(I)
               NMAX(NL) = N(I)
               LOLD = L(I)
            ELSEIF( L(I).EQ.LOLD .AND. N(I).GT.NMAX(NL) ) THEN
               NMAX(NL) = N(I)
            endif
         endif
      enddo
   end SUBROUTINE CALCNMAX

   !----------------------------------------------------------------
   SUBROUTINE READLA(STARTFILE,NUDSR,X)
      IMPLICIT REAL*8(A-H,O-Z)
      integer:: NUDSR
      CHARACTER*1:: whatToPlot(:),CFS
      CHARACTER*2:: CRR
      CHARACTER*40 STARTFILE
      double precision, intent(in):: x(:)

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN READP.'
         STOP
      endif

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
      endif

      !-- READH READS THE HEADER OF THE INPUTFILE AND DETERMINS WETHER
      CALL READH(12,LTR,NUDSR,LDR)

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
   end SUBROUTINE READLA

end PROGRAM LARA
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
