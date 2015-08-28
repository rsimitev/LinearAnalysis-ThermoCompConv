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

      !-- IF LRT.EQ.0 WE ARE LOOKING FOR THE RIGHT NUMBER OF DATASET NUDSR,
      !-- IF LRT.EQ.1 WE ARE LOOKING FOR THE RIGHT TIME.
      OPEN(12,FILE=STARTFILE,STATUS='old')

      !-- LST IS FORMAT PARAMETER ( L=1 FOR HIRSCHING FORMAT ) ,
      !   LCALCI TELLS HOW THE INPUTFILE HAS BEEN CALCULATED:
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
      IMPLICIT REAL*8(A-H,O-Z)
      
      WRITE(unitOut,'(3I5,2D16.8,'' M0, TRUNC'')')  M0, Truncation
      READ(unit_in,*) M0, Truncation
      READ(unit_in,*) TA, Rt, Rc, Pt, Pc, eta

      TAI=DSQRT(TAI)
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
