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
   integer, parameter:: nr=41, nt=360, np=720
   integer, PARAMETER:: NM = 5500, NAM = 400
   integer, PARAMETER:: NLMA = 100
   integer, PARAMETER:: NMX = 65, NMY = 128
   integer:: NMC
   integer:: LR, NQ, NR, n, l
   CHARACTER*40:: INPUTFILE,OUTPUTFILE
   CHARACTER*1:: CFS
   CHARACTER*2:: CRR
   INTEGER:: i, j
   double precision:: THETA(NMC)
   double precision:: dt

   double precision:: DX(NM)
   type eigenElement
      integer:: n
      integer:: l
      integer:: m
      character:: fieldCode
      double precision:: realPart
      double precision:: imagPart
   end type
   

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
      
      double precision, intent(in):: DX(:)
      double precision:: THETA(Nt), r, phi, dtheta

      write(*,*) 'computing the fields...'

      ! Avoid the poles
      dtheta  = 180.d0/(nt+1)
      DO I = 1, Nt
         THETA(I) = (I-0.5d0)*dtheta
      enddo

      !-- BESTIMMUNG DER PLM(THETA) , ABSPEICHERUNG:
      CALL STOREPLM(THETA, Nt_s)

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
   subroutine READD(unit_in,NKR,FieldValue,FieldCode,RealOrImag,LR,MR,NR,KR)
      IMPLICIT none
      CHARACTER(len=1):: FieldCode(:)
      CHARACTER(len=2):: RealOrImag(:)
      integer:: err

      DIMENSION FieldValue(*),LR(*),MR(*),NR(*),KR(*)

      NKR=0

      do
         READ(unit_in,'(1X,A1,3I3,2D16.8)',status=err) QuantI,LI,MI,NI,RealPart,ImagPart
         if (err.ne.0) exit
         NKR=NKR+1
         FieldCode(NKR)  = QuantI
         LR(NKR)         = LI
         MR(NKR)         = MI
         NR(NKR)         = NI
         FieldValue(NKR)=RealPart
         IF( MI.NE.0 ) THEN
            NKR=NKR+1
            FieldCode(NKR)=QuantI
            LR(NKR)=LI
            MR(NKR)=MI
            NR(NKR)=NI
            FieldValue(NKR)=ImagPart
         ENDIF
       enddo
   END subroutine

end PROGRAM LARA
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
