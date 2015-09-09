!************************************************************************
!-- PROGRAM TO read data from Galerkinprogram of J.W. and convert it for IDL.
!--
!-- Input:   stdin  (short version of LARA.EXE for DISSPLA)
!--
!-- Output:
!--
!------------------------------------------------------------------------
PROGRAM LARA
   use parameters
   implicit none
   double precision, parameter:: PI = 3.14159265358979D0
   integer, PARAMETER:: NM = 5500, NAM = 400
   integer, PARAMETER:: NLMA = 100
   integer, PARAMETER:: NMX = 65, NMY = 128
   integer:: NMC
   integer:: NCPLOT, LR, NQ, NR, n, l
   logical:: countourParIsNumber
   CHARACTER*40:: INPUTFILE,OUTPUTFILE
   CHARACTER*1:: CFS
   CHARACTER*2:: CRR
   INTEGER:: i, j
   double precision:: THETA(NMC)
   double precision:: dt

   double precision:: DX(NM)

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
   write(*,'(A,A,I3,D9.2)') 'reading data from ',INPUTFILE,'...'
   CALL READLA(INPUTFILE,DX)
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

   CALL READLA(INPUTFILE,DX)
   CALL PLO(DX)

   CLOSE(14)

contains

   !------------------------------------------------------------------------
   !     calculates the field Z and makes one subplot.
   SUBROUTINE PLO(dx)
      integer, parameter:: nr=41, nt=360, np=720
      double precision, intent(in):: DX(:)
      double precision:: THETA(Nt), r, phi, dtheta

      write(*,*) 'computing the fields...'

      ! Avoid the poles
      dtheta  = 180.d0/(nt+1)
      DO I = 1, Nt
         THETA(I) = (I-0.5d0)*dtheta
      enddo

      !-- BESTIMMUNG DER PLM(THETA) , ABSPEICHERUNG:
      CALL STOREPLM(THETA, Nt)

      open(20,file ='glo-render.dat',STATUS =  'UNKNOWN')
      do j=1, nt
         do k=0, np
            phi = k*(360.0/np)
            do i=0, nr
               r = ri + dble(i)/dble(nr)
               Z = flow_r(DX, r, PHI, j)
               
               Write(20,*) '#,nmr,x,nmp,x,nmt'
               Write(20,*) '#',nmr,'x',nmp,'x',nmt
               Write(20,*) '#========================'
               Write(20,*) r*sin((theta(j)-90.)*pi/180.0)*cos(phi*pi/180.0), &
                           r*sin((theta(j)-90.)*pi/180.0)*sin(phi*pi/180.0), &
                           r*cos((theta(j)-90.)*pi/180.0), &
                           z
            enddo
         enddo
      enddo
      close(20)
   end subroutine plo


   !------------------------------------------------------------------------
   !> Radiales flow: U_r = L_2/r v
   double precision function flow_r(X,R,PHI,iTHETA)
      IMPLICIT none
      double precision, intent(in):: x(:), R, iTheta, phi

      PPHI = PHI*PI/180.D0
      RI = ETA/(1.D0-ETA)
      RF = 0.D0

      DO I = NDOMIN,NDOMAX
         IF( M(I).EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         endif

         RFT = EPSM*L(I)*(L(I)+1)*PLMS(L(I),M(I),iTheta)*DSIN( N(I)*PI*(R-RI) ) / R

         RFT = RFT

         IF( CRR(I).EQ.'RR' ) THEN
            RFT = RFT * X(I) * DCOS( M(I)*PPHI )
         ELSEIF( CRR(I).EQ.'IR' ) THEN
            RFT = -RFT * X(I) * DSIN( M(I)*PPHI )
         ELSE
            RFT = 0.D0
         endif

         flow_r = flow_r + RFT
      enddo
   end function

   !------------------------------------------------------------------------
   !   Temperaturfeld Theta ( =  Abweichung vom Grundzust.)
   !   optimized for K = 0.
   double precision function TEMP(X,whatToPlot,R,PHI,NTHETA)
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
   double precision function FT(X,whatToPlot,R,PHI,NTHETA)
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
   double precision function FP(X,whatToPlot,R,PHI,NTHETA)
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
   double precision function FR(X,whatToPlot,R,PHI,NTHETA)
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
   double precision function TT(X,whatToPlot,R,PHI,NTHETA)
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
   double precision function localNusselt(X,whatToPlot,R,PHI,NTHETA)
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
   double precision function flow_p_zonal(X,R,NTHETA)
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
   double precision function flow_p(X,whatToPlot,R,PHI,NTHETA)
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
   double precision function DMPJ(X,whatToPlot,R,NTHETA)
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
   double precision function DMC(X,R,NTHETA)
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
   SUBROUTINE READLA(STARTFILE,X)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CF,CFS
      CHARACTER*2 CRR
      CHARACTER*40 STARTFILE
      character(len=20):: title
      double precision, intent(out):: x(:)

      IF( NM.NE.NMC ) THEN
         WRITE(*,*) 'WRONG DIMENSION NM IN READP.'
         STOP
      ENDIF

      OPEN(12,FILE=STARTFILE,STATUS='old')
      READ(12,*) title

      LT=0

      !-- READH READS THE HEADER OF THE INPUTFILE AND DETERMINS WETHER
      !   THE DATASET (dataSetNumber,TIME) IS THE RIGHT ONE (LDR=1):
      CALL READH(12)

      !-- LOOKING FOR THE RIGHT DATASET:
      DO I=1,1000
         !----- READD READS FULL SET OF COEFFITIENTS:
         CALL READD(12,ND,X,CF,CRR,L,M,N,K)
         IF( LDR.EQ.1 ) exit
      enddo

      if(ND.GT.NM) then
        write(*,*) 'To small dimension NM in READLA'
        stop
      endif
      TA=TA**2

      LSX=1
      CALL SORTK(ND,LSX,X,CF,CRR,L,M,N,K,NUC,NUOM)
      CALL RDIM(ND,CF,CRR,L,M,N,K,CFS,LS,MS,NS, &
                       NDV,NDW,NDT,NDH,NDG,NDVS,NDWS,NDTS,NDHS,NDGS,NDS)
      CLOSE(12)
   END SUBROUTINE READLA

   !--------------------------------------------------------------------------
   subroutine READH(unit_in)
      IMPLICIT none
      READ(unit_in,*) M0, Truncation
      READ(unit_in,*) TA, Rt, Rc, Pt, Pc, eta
   END subroutine readh

   !--------------------------------------------------------------------------
   subroutine READD(unit_in,NKR,XR,CFR,CRRR,LR,MR,NR,KR)
      IMPLICIT none
      CHARACTER(len=1):: CFR(:)
      CHARACTER(len=2):: CRRR(:)
      integer:: err

      DIMENSION XR(*),LR(*),MR(*),NR(*),KR(*)

      NKR=0

      do
         READ(unit_in,'(1X,A1,3I3,2D16.8)',status=err) CFI,LI,MI,NI,DRR,DIR
         if (err.ne.0) exit
         NKR=NKR+1
         CFR(NKR)=CFI
         CRRR(NKR)='RR'
         LR(NKR)=LI
         MR(NKR)=MI
         NR(NKR)=NI
         XR(NKR)=DRR
         IF( MI.NE.0 ) THEN
            NKR=NKR+1
            CFR(NKR)=CFI
            CRRR(NKR)='IR'
            LR(NKR)=LI
            MR(NKR)=MI
            NR(NKR)=NI
            XR(NKR)=DIR
         ENDIF
       enddo
   END subroutine

! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab



end PROGRAM LARA
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
