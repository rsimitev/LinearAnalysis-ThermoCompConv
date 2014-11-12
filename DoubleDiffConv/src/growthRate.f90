module GrowthRateMod
#include "losub-inc.h"
   use parameters
   implicit none
   double precision, private:: ZACOPY(NMAX,NMAX),ZBCOPY(NMAX,NMAX)
   double precision, parameter, private:: DPI=3.141592653589793D0
contains

   subroutine GrowthRate_init()
      implicit none
      ZACOPY = dcmplx(0.0d0,0.0d0)
      ZBCOPY = dcmplx(0.0d0,0.0d0)
   end subroutine

   !***********************************************************************
   !> Computes the maximum imaginary part of the frequency,
   !! that is, the maximum growth rate for all eigen modes.
   double precision FUNCTION MaxGrowthRate(Ra)
      IMPLICIT none
      double precision, intent(in):: Ra
      double precision:: RtOld
      double complex:: ZEW(NMAX)

      RtOld = Rt
      Rt    = Ra
      Zew = dcmplx(0.0d0,0.0d0)
      call computeGrowthRateModes(.false., zew)
      ! search for lowest imaginary part:
      MaxGrowthRate = minval(DIMAG(ZEW),1)
      Rt=RtOld
   end function

   !> Computes the complex frequency for which the growth rate is maximum.
   double complex FUNCTION MaxGrowthRateCmplx(Ra)
      IMPLICIT none
      double precision, intent(in):: Ra
      double precision:: RtOld
      double complex:: ZEW(NMAX)
      integer:: imin

      RtOld = Rt
      Rt    = Ra
      Zew = dcmplx(0.0d0,0.0d0)
      call computeGrowthRateModes(.false., zew)
      ! search for lowest imaginary part:
      IMIN = minloc(DIMAG(ZEW),1)
      MaxGrowthRateCmplx = ZEW(IMIN)
      Rt=RtOld
   end function

   subroutine computeGrowthRateModes(sort, zew, zeval)
      implicit none
      !> .True. Sort the eigenvalues and eigenvectors if computed.
      logical, intent(in):: sort
      double complex, intent(out):: ZEW(NMAX)
      double complex, intent(out), optional::ZEVAL(NMAX,NMAX)
      double complex:: ZA(NMAX,NMAX),ZB(NMAX,NMAX)
      double complex:: ZEWA(NMAX),ZEWB(NMAX)
      double complex:: ZEVEC(NMAX), ZSAVE, ZWORK(3*NMAX)
      double precision:: RWORK(8*NMAX)
      integer:: i, j, k, info

      ! - MAT SETS THE complex(8) MATRICES ZA AND ZB SETTING OF MATRIX:
      CALL MAT(ZA,ZB,NMAX)

!       SUBROUTINE zGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,
!     $                  VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
      IF(present(zeval)) THEN ! Compute eigen values and vectors.
        call zggev('N', 'V',NMAX,ZA,NMAX,ZB,NMAX, ZEWA, ZEWB, ZEVAL, NMAX, ZEVAL, NMAX, ZWORK, 3*NMAX, rwork, info)
      ELSE ! Only compute eigenvalues
        call zggev('N', 'N',NMAX,ZA,NMAX,ZB,NMAX, ZEWA, ZEWB, ZEVAL, NMAX, ZEVAL, NMAX, ZWORK, 3*NMAX, rwork, info)
      endIF

      ZEW(:)=ZEWA(:)/ZEWB(:)

      ! sort eigenvalues:
      if(sort) then
         DO I=1,ND
            DO J=I,ND
               IF( DIMAG(ZEW(J)).LT.DIMAG(ZEW(I)) ) THEN
                  ZSAVE = ZEW(J)
                  ZEW(J) = ZEW(I)
                  ZEW(I) = ZSAVE
                  IF(present(zeval)) THEN
                     DO K=1,ND
                        ZEVEC(K)    = ZEVAL(K,J)
                        ZEVAL(K,J) = ZEVAL(K,I)
                        ZEVAL(K,I) = ZEVEC(K)
                     enddo
                  endIF
               endIF
            enddo
         enddo
      endif
   end subroutine computeGrowthRateModes


   !************************************************************************
   ! - SETS THE complex(8) MATRICES ZA AND ZB.
   SUBROUTINE MAT(ZA,ZB,NDIM)
      implicit none
      double complex:: ZA(ndim,ndim),ZB(ndim,ndim)
      integer:: ndim, lmax
      integer:: ni, i, nimax, li, lpi, lti
      integer:: nj, j, njmax, lj, lpj, ltj


      ZA(:,:)=DCMPLX(0D0,0D0)
      ZB(:,:)=DCMPLX(0D0,0D0)

      I=0
      LMAX=2*NT+M0-1

      DO 4000 LI=LMIN,LMAX,LD
         LPI=LI

   !        Determine L for toroidal (w) field:
         IF( NE.EQ.2 ) THEN
            LTI=LI+1
         ELSEIF( NE.EQ.1 ) THEN
            LTI=LI-1
         ELSEIF( NE.EQ.0 ) THEN
            LTI=LI
         endIF

         NIMAX=INT( DBLE(2*NT+1-LI+M0)/2 )
         DO 3000 NI=1,NIMAX
            J=0
            DO 2000 LJ=LMIN,LMAX,LD
               LPJ=LJ
               IF( NE.EQ.2 ) THEN
                  LTJ=LJ+1
               ELSEIF( NE.EQ.1 ) THEN
                  LTJ=LJ-1
               ELSEIF( NE.EQ.0 ) THEN
                  LTJ=LJ
               endIF
               NJMAX=INT( DBLE(2*NT+1-LJ+M0)/2 )

   !  ******************** I: Equation (Line) ******************
   !  ******************** J: Variable (Column) ****************
   !  ******************** I+1: v (poloidal)  ******************
   !  ******************** I+2: theta         ******************
   !  ******************** I+3: w (toroidal)  ******************
   ! new****************** I+4: gamma (concentration) **********
               DO 1000 NJ=1,NJMAX
               IF(J+3.GT.NDIM .OR. I+3.GT.NDIM) THEN
                  write(*,*) 'MAT(): NDIM too small.'
                  stop
               endIF
               IF( LI.EQ.LJ ) THEN
                     ZB(I+1,J+1)=DCMPLX(0.D0,-DIII2(NI,NJ,LPI,1))
                     ZA(I+1,J+1)=DCMPLX(DIII1(NI,NJ,LPI),DIII3(NI,NJ,LPI,1))
                     ZA(I+1,J+2)=DCMPLX(DIII5(NI,NJ,LPI),0.D0)
                     ! -- concentration driving
                     ZA(I+1,J+4)=DCMPLX(DIII5conc(NI,NJ,LPI),0.D0)

                     ZB(I+2,J+2)=DCMPLX(0.D0,-DI1(NI,NJ,1))
                     ZA(I+2,J+1)=DCMPLX(DI3(NI,NJ,LPI),0.D0)
                     ZA(I+2,J+2)=DCMPLX(DI2(NI,NJ,LPI),0.D0)
                     ZB(I+3,J+3)=DCMPLX(0.D0,-DII2(NI,NJ,LTI,1))
                     ZA(I+3,J+3)=DCMPLX(DII1(NI,NJ,LTI),DII3(NI,NJ,1))
                     ! -- concentration equation
                     ZB(I+4,J+4)=DCMPLX(0.D0,-DI1(NI,NJ,1))
                     ZA(I+4,J+1)=DCMPLX(DI3(NI,NJ,LPI),0.D0)
                     ZA(I+4,J+4)=DCMPLX(1.D0/Le * DI2(NI,NJ,LPI),0.D0)
                  endIF
                  IF( LPI.EQ.LTJ+1 ) THEN
                      ZA(I+1,J+3)=DCMPLX(DIII4A(NI,NJ,LPI,1),0.D0)
                  ELSEIF( LPI.EQ.LTJ-1 ) THEN
                      ZA(I+1,J+3)=DCMPLX(DIII4B(NI,NJ,LPI,1),0.D0)
                  endIF
                  IF( LTI.EQ.LPJ+1 ) THEN
                      ZA(I+3,J+1)=DCMPLX(DII4A(NI,NJ,LTI,1),0.D0)
                  ELSEIF( LTI.EQ.LPJ-1 ) THEN
                      ZA(I+3,J+1)=DCMPLX(DII4B(NI,NJ,LTI,1),0.D0)
                  endIF
                  J=J+4
1000           CONTINUE
2000        CONTINUE
             I=I+4
3000     CONTINUE
4000  CONTINUE
   end subroutine

   !************************************************************************
   ! - GALERKIN TERMS:
   !************************************************************************
   double precision function DI1(N1,N2,NU1)
   ! ---- HEAT EQUATION, TIME DERIVATIVE
      implicit none
      integer:: N1, N2, NU1
      DI1=Pt*NU1*R('SS ',2,N1,N2,0)
   end function

   double precision function DI2(N1,N2,L1)
   ! ---- HEAT EQUATION , DISSIPATION
      implicit none
      integer:: N1, N2, l1
      DI2=N2**2*DPI**2*R('SS ',2,N1,N2,0) - 2*N2*DPI*R('SC ',1,N1,N2,0) + DL(L1)*R('SS ',0,N1,N2,0)
   end function

   double precision function DI3(N1,N2,L1)
   ! ---- HEAT EQUATION , SOURCE
      implicit none
      integer:: N1, N2, l1
      DI3 =-DL(L1)*R('SS ',2,N1,N2,0)
   end function

   double precision function DII1(N1,N2,L1)
   ! ---- TOROIDAL EQUATION , DISSIPATION
      implicit none
      integer:: N1, N2, l1
      DII1=DL(L1)*( (N2-1)**2*DPI**2*R('CC ',4,N1-1,N2-1,0) + 4*(N2-1)*DPI*R('CS ',3,N1-1,N2-1,0) + (DL(L1)-2)*R('CC ',2,N1-1,N2-1,0) )
   end function

   double precision function DII2(N1,N2,L1,NU1)
   ! ---- TOROIDAL EQUATION , TIME DERIVATIVE
      implicit none
      integer:: N1, N2, l1, NU1
      DII2=NU1*DL(L1)*R('CC ',4,N1-1,N2-1,0)
   end function

   double precision function DII3(N1,N2,NU1)
   ! ---- TOROIDAL EQUATION , CORRIOLIS
      implicit none
      integer:: N1, N2, NU1
      DII3=-TAU*NU1*M0*R('CC ',4,N1-1,N2-1,0)
   end function

   double precision function DII4A(N1,N2,L1,NU1)
   ! ---- TOROIADL EQUATION , Q-TERM 1 (L1=L3+1)
      implicit none
      integer:: N1, N2, l1, NU1
      DII4A= TAU * DSQRT( DBLE(L1-NU1*M0)*(L1+NU1*M0)/(2*L1-1)/(2*L1+1) ) * ( (L1**2-1)*(L1-1)*R('CS ',2,N1-1,N2,0)  - (L1+1)*(L1-1)*N2*DPI*R('CC ',3,N1-1,N2,0)    )
   end function

   double precision function DII4B(N1,N2,L1,NU1)
   ! ---- TOROIADL EQUATION , Q-TERM 1 (L1=L3-1)
      implicit none
      integer:: N1, N2, l1, NU1
      DII4B= TAU * DSQRT( DBLE(L1-NU1*M0+1)*(L1+NU1*M0+1)/(2*L1+1)/(2*L1+3) ) * ( (1-(L1+1)**2)*(L1+2)*R('CS ',2,N1-1,N2,0)  - L1*(L1+2)*N2*DPI*R('CC ',3,N1-1,N2,0)  )
   end function

   double precision function DIII1(N1,N2,L1)
   ! ---- POLOIDAL EQUOATION , DISSIPATION
      implicit none
      integer:: N1, N2, l1
      DIII1=DL(L1)* ( N2**4*DPI**4*R('SS ',2,N1,N2,0) - 4*N2**3*DPI**3*R('SC ',1,N1,N2,0) + 2*DL(L1)*N2**2*DPI**2*R('SS ',0,N1,N2,0) + (DL(L1)**2-2*DL(L1))*R('SS ',-2,N1,N2,0) )
   end function

   double precision function DIII2(N1,N2,L1,NU1)
   ! ---- POLOIDAL EQUATION , TIME DERIVATIVE
      implicit none
      integer:: N1, N2, l1, NU1
      DIII2= -NU1*DL(L1)*(-N2**2*DPI**2*R('SS ',2,N1,N2,0) + 2*N2*DPI*R('SC ',1,N1,N2,0)-DL(L1)*R('SS ',0,N1,N2,0) )
   end function

   double precision function DIII3(N1,N2,L1,NU1)
   ! ---- POLOIDAL EQUATION , CORRIOLIS
      implicit none
      integer:: N1, N2, l1, NU1
      DIII3= TAU*NU1*M0*( -N2**2*DPI**2*R('SS ',2,N1,N2,0) + 2*N2*DPI*R('SC ',1,N1,N2,0) - DL(L1)*R('SS ',0,N1,N2,0) )
   end function

   double precision function DIII4A(N1,N2,L1,NU1)
   ! ---- POLOIDAL EUQUATION , Q-TERM 1 (L1=L3+1)
      implicit none
      integer:: N1, N2, l1, NU1
      DIII4A= TAU * DSQRT(DBLE(L1-M0*NU1)*(L1+M0*NU1)/(2*L1-1)/(2*L1+1)) * ( (L1*(L1-1)-2)*(L1-1)*R('SC ',2,N1,N2-1,0) + (L1+1)*(L1-1)*(N2-1)*DPI*R('SS ',3,N1,N2-1,0) )
   end function

   double precision function DIII4B(N1,N2,L1,NU1)
   ! ---- POLOIDAL EQUATION , Q-TERM 2 (L1=L3-1)
      implicit none
      integer:: N1, N2, l1, NU1
      DIII4B= TAU * DSQRT( DBLE(L1-M0*NU1+1)*(L1+M0*NU1+1)/(2*L1+1)/(2*L1+3) ) * ( (L1+2)*(2-(L1+1)*(L1+2))*R('SC ',2,N1,N2-1,0) + L1*(L1+2)*(N2-1)*DPI*R('SS ',3,N1,N2-1,0) )
   end function

   double precision function DIII5(N1,N2,L1)
   ! ---- POLOIDAL EQUATION ,
      implicit none
      integer:: N1, N2, l1
      DIII5=-Rt*DL(L1)*R('SS ',2,N1,N2,0)
   end function

   double precision function DIII5conc(N1,N2,L1)
   ! ---- POLOIDAL EQUATION ,
      implicit none
      integer:: N1, N2, l1
      DIII5conc=-Rc*DL(L1)*R('SS ',2,N1,N2,0)
   end function

   !**************************************************************************
   !-- SUBROUTINES:
   !**************************************************************************
   double precision function DL(L)
      implicit none
      integer:: l
      DL = DBLE(L*(L+1))
   end function

   !***************************************************************************
   SUBROUTINE DIMENSION(LMIN,LD,NT,M0,ND)
   !***************************************************************************
      implicit none
      integer:: NT,M0,ND,LMIN,LD
      integer:: L
   ! - DETERMINATION OF DIMENSION:
   ! - for each value of L the number of possible N-values is added
   !         print*, "Triangular truncation (2.12)"
   !         print*, LMIN, "...", 2*NT+M0-1,LD
      ND=0
      DO L = LMIN, 2*NT+M0-1, LD
   !         print*, L, 1, "...", INT( DBLE(2*NT+1-L+M0)/2 )
   ! cccccccc18    ND=ND+3*DINT( DBLE(2*NT+1-L+M0)/2 )
         ND = ND + 4*INT( DBLE(2*NT+1-L+M0)/2 )
      endDO

      IF(ND.GT.NMAX) THEN
         WRITE(*,*) 'DIMENSION OF MATRIX TOO SMALL:',ND,'>',NMAX
         STOP DIM_TO_SMALL
      endIF
   end subroutine

!-----------------------------------------------------------------------
   double precision function R (TRII, NEI, II, JI, KI)
!-----------------------------------------------------------------------
!     BERECHNET DIE RADIALINTEGRALE FUER DAS DYNAMOPROBLEM
!-----------------------------------------------------------------------
      IMPLICIT double precision (A - H, O - Y)
      IMPLICIT complex (8) (Z)
      CHARACTER(3) TRI, TRII
      integer:: i,j,k,NEI, ii, ji, ki
!
      TRI = TRII
      NE = NEI
      I = II
      J = JI
      K = KI
      CALL TAUSCH (TRI, I, J, K)
      IF (NE.EQ.0) THEN
         IF (TRI.EQ.'SS ') THEN
            IF (I * J.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0 * (C0 (I - J) - C0 (I + J) )
            ENDIF
         ELSEIF (TRI.EQ.'SC ') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0 * (S0 (I - J) + S0 (I + J) )
            ENDIF
         ELSEIF (TRI.EQ.'CC ') THEN
            RINT = 1D0 / 2D0 * (C0 (I + J) + C0 (I - J) )
         ELSEIF (TRI.EQ.'SSS') THEN
            IF (I * J * K.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0 * (S0 (I - J + K) + S0 (I + J - K)      &
               - S0 (I - J - K) - S0 (I + J + K) )
            ENDIF
         ELSEIF (TRI.EQ.'SSC') THEN
            IF (I * J.NE.0) THEN
               RINT = 1D0 / 4D0 * (C0 (I - J - K) - C0 (I + J + K)      &
               + C0 (I - J + K) - C0 (I + J - K) )
            ELSE
               RINT = 0D0
            ENDIF
         ELSEIF (TRI.EQ.'SCC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0 * (S0 (I + J - K) + S0 (I + J + K)      &
               + S0 (I - J + K) + S0 (I - J - K) )
            ENDIF
         ELSEIF (TRI.EQ.'CCC') THEN
            RINT = 1D0 / 4D0 * (C0 (I + J + K) + C0 (I + J - K) + C0 (I &
            - J + K) + C0 (I - J - K) )
         ENDIF
      ELSEIF (NE.EQ.1) THEN
         IF (TRI.EQ.'SS ') THEN
            IF (I * J.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0 * (C1 (I - J) - C1 (I + J) )
            ENDIF
         ELSEIF (TRI.EQ.'SC ') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0 * (S1 (I - J) + S1 (I + J) )
            ENDIF
         ELSEIF (TRI.EQ.'CC ') THEN
            RINT = 1D0 / 2D0 * (C1 (I + J) + C1 (I - J) )
         ELSEIF (TRI.EQ.'SSS') THEN
            IF (I * J * K.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0 * (S1 (I - J + K) + S1 (I + J - K)      &
               - S1 (I - J - K) - S1 (I + J + K) )
            ENDIF
         ELSEIF (TRI.EQ.'SSC') THEN
            IF (I * J.NE.0) THEN
               RINT = 1D0 / 4D0 * (C1 (I - J - K) - C1 (I + J + K)      &
               + C1 (I - J + K) - C1 (I + J - K) )
            ELSE
               RINT = 0D0
            ENDIF
         ELSEIF (TRI.EQ.'SCC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0 * (S1 (I + J - K) + S1 (I + J + K)      &
               + S1 (I - J + K) + S1 (I - J - K) )
            ENDIF
         ELSEIF (TRI.EQ.'CCC') THEN
            RINT = 1D0 / 4D0 * (C1 (I + J + K) + C1 (I + J - K) + C1 (I &
            - J + K) + C1 (I - J - K) )
         ENDIF
      ELSEIF (NE.EQ.2) THEN
         IF (TRI.EQ.'SS ') THEN
            IF (I * J.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0 * (C2 (I - J) - C2 (I + J) )
            ENDIF
         ELSEIF (TRI.EQ.'SC ') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0 * (S2 (I - J) + S2 (I + J) )
            ENDIF
         ELSEIF (TRI.EQ.'CC ') THEN
            RINT = 1D0 / 2D0 * (C2 (I + J) + C2 (I - J) )
         ELSEIF (TRI.EQ.'SSS') THEN
            IF (I * J * K.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0 * (S2 (I - J + K) + S2 (I + J - K)      &
               - S2 (I - J - K) - S2 (I + J + K) )
            ENDIF
         ELSEIF (TRI.EQ.'SSC') THEN
            IF (I * J.NE.0) THEN
               RINT = 1D0 / 4D0 * (C2 (I - J - K) - C2 (I + J + K)      &
               + C2 (I - J + K) - C2 (I + J - K) )
            ELSE
               RINT = 0D0
            ENDIF
         ELSEIF (TRI.EQ.'SCC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0 * (S2 (I + J - K) + S2 (I + J + K)      &
               + S2 (I - J + K) + S2 (I - J - K) )
            ENDIF
         ELSEIF (TRI.EQ.'CCC') THEN
            RINT = 1D0 / 4D0 * (C2 (I + J + K) + C2 (I + J - K) + C2 (I &
            - J + K) + C2 (I - J - K) )
         ENDIF
      ELSEIF (NE.EQ.3) THEN
         IF (TRI.EQ.'SS') THEN
            IF (I * J.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0 * (C3 (I - J) - C3 (I + J) )
            ENDIF
         ELSEIF (TRI.EQ.'SC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0 * (S3 (I - J) + S3 (I + J) )
            ENDIF
         ELSEIF (TRI.EQ.'CC') THEN
            RINT = 1D0 / 2D0 * (C3 (I + J) + C3 (I - J) )
         ELSEIF (TRI.EQ.'SSS') THEN
            IF (I * J * K.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0 * (S3 (I - J + K) + S3 (I + J - K)      &
               - S3 (I - J - K) - S3 (I + J + K) )
            ENDIF
         ELSEIF (TRI.EQ.'SSC') THEN
            IF (I * J.NE.0) THEN
               RINT = 1D0 / 4D0 * (C3 (I - J - K) - C3 (I + J + K)      &
               + C3 (I - J + K) - C3 (I + J - K) )
            ELSE
               RINT = 0D0
            ENDIF
         ELSEIF (TRI.EQ.'SCC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0 * (S3 (I + J - K) + S3 (I + J + K)      &
               + S3 (I - J + K) + S3 (I - J - K) )
            ENDIF
         ELSEIF (TRI.EQ.'CCC') THEN
            RINT = 1D0 / 4D0 * (C3 (I + J + K) + C3 (I + J - K) + C3 (I &
            - J + K) + C3 (I - J - K) )
         ENDIF
      ELSEIF (NE.EQ.4) THEN
         IF (TRI.EQ.'SS') THEN
            IF (I * J.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0 * (C4 (I - J) - C4 (I + J) )
            ENDIF
         ELSEIF (TRI.EQ.'SC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0 * (S4 (I - J) + S4 (I + J) )
            ENDIF
         ELSEIF (TRI.EQ.'CC') THEN
            RINT = 1D0 / 2 * (C4 (I + J) + C4 (I - J) )
         ELSEIF (TRI.EQ.'SSS') THEN
            IF (I * J * K.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0 * (S4 (I - J + K) + S4 (I + J - K)      &
               - S4 (I - J - K) - S4 (I + J + K) )
            ENDIF
         ELSEIF (TRI.EQ.'SSC') THEN
            IF (I * J.NE.0) THEN
               RINT = 1D0 / 4D0 * (C4 (I - J - K) - C4 (I + J + K)      &
               + C4 (I - J + K) - C4 (I + J - K) )
            ELSE
               RINT = 0D0
            ENDIF
         ELSEIF (TRI.EQ.'SCC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0 * (S4 (I + J - K) + S4 (I + J + K)      &
               + S4 (I - J + K) + S4 (I - J - K) )
            ENDIF
         ELSEIF (TRI.EQ.'CCC') THEN
            RINT = 1D0 / 4D0 * (C4 (I + J + K) + C4 (I + J - K) + C4 (I &
            - J + K) + C4 (I - J - K) )
         ENDIF
      ELSEIF (NE.EQ. - 1) THEN
         IF (TRI.EQ.'SS ') THEN
            IF (I * J.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0 * (CM1 (I - J) - CM1 (I + J) )
            ENDIF
         ELSEIF (TRI.EQ.'SC ') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0 * (SM1 (I - J) + SM1 (I + J) )
            ENDIF
         ELSEIF (TRI.EQ.'CC ') THEN
            RINT = 1D0 / 2D0 * (CM1 (I + J) + CM1 (I - J) )
         ELSEIF (TRI.EQ.'SSS') THEN
            IF (I * J * K.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0 * (SM1 (I - J + K) + SM1 (I + J - K)    &
               - SM1 (I - J - K) - SM1 (I + J + K) )
            ENDIF
         ELSEIF (TRI.EQ.'SSC') THEN
            IF (I * J.NE.0) THEN
               RINT = 1D0 / 4D0 * (CM1 (I - J - K) - CM1 (I + J + K)    &
               + CM1 (I - J + K) - CM1 (I + J - K) )
            ELSE
               RINT = 0D0
            ENDIF
         ELSEIF (TRI.EQ.'SCC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0 * (SM1 (I + J - K) + SM1 (I + J + K)    &
               + SM1 (I - J + K) + SM1 (I - J - K) )
            ENDIF
         ELSEIF (TRI.EQ.'CCC') THEN
            RINT = 1D0 / 4D0 * (CM1 (I + J + K) + CM1 (I + J - K)       &
            + CM1 (I - J + K) + CM1 (I - J - K) )
         ENDIF
      ELSEIF (NE.EQ. - 2) THEN
         IF (TRI.EQ.'SS ') THEN
            IF (I * J.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0 * (CM2 (I - J) - CM2 (I + J) )
            ENDIF
         ELSEIF (TRI.EQ.'SC ') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0 * (SM2 (I - J) + SM2 (I + J) )
            ENDIF
         ELSEIF (TRI.EQ.'CC ') THEN
            RINT = 1D0 / 2D0 * (CM2 (I + J) + CM2 (I - J) )
         ELSEIF (TRI.EQ.'SSS') THEN
            IF (I * J * K.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0 * (SM2 (I - J + K) + SM2 (I + J - K)    &
               - SM2 (I - J - K) - SM2 (I + J + K) )
            ENDIF
         ELSEIF (TRI.EQ.'SSC') THEN
            IF (I * J.NE.0) THEN
               RINT = 1D0 / 4D0 * (CM2 (I - J - K) - CM2 (I + J + K)    &
               + CM2 (I - J + K) - CM2 (I + J - K) )
            ELSE
               RINT = 0D0
            ENDIF
         ELSEIF (TRI.EQ.'SCC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0 * (SM2 (I + J - K) + SM2 (I + J + K)    &
               + SM2 (I - J + K) + SM2 (I - J - K) )
            ENDIF
         ELSEIF (TRI.EQ.'CCC') THEN
            RINT = 1D0 / 4D0 * (CM2 (I + J + K) + CM2 (I + J - K)       &
            + CM2 (I - J + K) + CM2 (I - J - K) )
         ENDIF
      ENDIF
      R = RINT
!      WRITE(12,'(I4,A3,3I4,D14.4)') NE,TRI,I,J,K,R
   END FUNCTION R
!-----------------------------------------------------------------------
!     END OF RI
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   integer function NGU (N)
        implicit none
        integer:: n
        NGU = - 1
        IF (MOD (N, 2) .EQ.0) NGU = 1
   END FUNCTION NGU
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   SUBROUTINE TAUSCH (TRI, I, J, K)
      implicit none
      CHARACTER(3) TRI
      integer:: i, j, k, N
      IF (TRI.EQ.'CS ') THEN
         N = I
         I = J
         J = N
         TRI = 'SC '
      ELSEIF (TRI.EQ.'SCS') THEN
         N = J
         J = K
         K = N
         TRI = 'SSC'
      ELSEIF (TRI.EQ.'CSS') THEN
         N = I
         I = K
         K = N
         TRI = 'SSC'
      ELSEIF (TRI.EQ.'CCS') THEN
         N = I
         I = K
         K = N
         TRI = 'SCC'
      ELSEIF (TRI.EQ.'CSC') THEN
         N = I
         I = J
         J = N
         TRI = 'SCC'
      ENDIF
      END SUBROUTINE TAUSCH
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function S0 (N)
      IMPLICIT double precision (A - H, O - Y)
      integer:: N
      IF (NGU (N) .EQ.1) THEN
         S0 = 0D0
      ELSE
         S0 = 2D0 / N / DPI
      ENDIF
      END FUNCTION S0
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function C0 (N)
      IMPLICIT double precision (A - H, O - Y)
      integer:: N
      IF (N.EQ.0) THEN
         C0 = 1D0
      ELSE
         C0 = 0D0
      ENDIF
      END FUNCTION C0
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function S1 (N)
      IMPLICIT double precision (A - H, O - Y)
      integer:: N
      IF (NGU (N) .EQ.1) THEN
         IF (N.EQ.0) THEN
            S1 = 0D0
         ELSE
            S1 = - 1D0 / N / DPI
         ENDIF
      ELSE
         S1 = (1 + 2 * RI) / N / DPI
      ENDIF
      END FUNCTION S1
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function C1 (N)
      IMPLICIT double precision (A - H, O - Y)
      integer:: N
      IF (NGU (N) .EQ.1) THEN
         IF (N.EQ.0) THEN
            C1 = RI + 1D0 / 2
         ELSE
            C1 = 0D0
         ENDIF
      ELSE
         C1 = - 2D0 / N**2 / DPI**2
      ENDIF
      END FUNCTION C1
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function S2 (N)
      IMPLICIT double precision (A - H, O - Y)
      integer:: N
      IF (NGU (N) .EQ.1) THEN
         IF (N.EQ.0) THEN
            S2 = 0D0
         ELSE
            S2 = - (2 * RI + 1) / DPI / N
         ENDIF
      ELSE
         S2 = (1 + 2 * RI * (RI + 1) ) / DPI / N - 4D0 / N**3 / DPI**3
      ENDIF
      END FUNCTION S2
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function C2 (N)
!-----------------------------------------------------------------------
      IMPLICIT double precision (A - H, O - Y)
      integer:: N
      IF (NGU (N) .EQ.1) THEN
         IF (N.EQ.0) THEN
            C2 = RI * (RI + 1) + 1D0 / 3
         ELSE
            C2 = 2D0 / N**2 / DPI**2
         ENDIF
      ELSE
         C2 = - 2D0 * (2 * RI + 1) / N**2 / DPI**2
      ENDIF
      END FUNCTION C2
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function S3 (N)
      IMPLICIT double precision (A - H, O - Y)
      integer:: N
      IF (NGU (N) .EQ.1) THEN
         IF (N.EQ.0) THEN
            S3 = 0D0
         ELSE
            S3 = - (1 + 3 * RI + 3 * RI**2) / DPI / N + 6 / DPI**3 / N**3
         ENDIF
      ELSE
         S3 = (1 + 3 * RI + 3 * RI**2 + 2 * RI**3) / DPI / N - (6 + 12 * RI) / DPI**3 / N**3
      ENDIF
      END FUNCTION S3
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function C3 (N)
!-----------------------------------------------------------------------
      IMPLICIT double precision (A - H, O - Y)
      integer:: N
      IF (NGU (N) .EQ.1) THEN
         IF (N.EQ.0) THEN
            C3 = 1D0 / 4 + RI + 3D0 / 2 * RI**2 + RI**3
         ELSE
            C3 = (3 + 6 * RI) / DPI**2 / N**2
         ENDIF
      ELSE
         C3 = - (3 + 6 * RI + 6 * RI**2) / DPI**2 / N**2 + 12 / DPI**4 / N**4
      ENDIF
      END FUNCTION C3
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function S4 (N)
      IMPLICIT double precision (A - H, O - Y)
      integer:: N
      IF (NGU (N) .EQ.1) THEN
         IF (N.EQ.0) THEN
            S4 = 0D0
         ELSE
            S4 = - (1 + 4 * RI + 6 * RI**2 + 4 * RI**3) / DPI / N + (12 + 24 * RI) / DPI**3 / N**3
         ENDIF
      ELSE
         S4 = (1 + 4 * RI + 6 * RI**2 + 4 * RI**3 + 2 * RI**4) / DPI /  &
         N - (12 + 24 * RI + 24 * RI**2) / DPI**3 / N**3 + 48D0 / DPI**5 / N**5
      ENDIF
      END FUNCTION S4
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function C4 (N)
!-----------------------------------------------------------------------
      IMPLICIT double precision (A - H, O - Y)
      integer:: N
      IF (NGU (N) .EQ.1) THEN
         IF (N.EQ.0) THEN
            C4 = 1D0 / 5 + RI + 2 * RI**2 + 2 * RI**3 + RI**4
         ELSE
            C4 = (4 + 12 * RI + 12 * RI**2) / DPI**2 / N**2 - 24 / DPI**&
            4 / N**4
         ENDIF
      ELSE
         C4 = - (4 + 12 * RI + 12 * RI**2 + 8 * RI**3) / DPI**2 / N**2 +&
         (24 + 48 * RI) / DPI**4 / N**4
      ENDIF
      END FUNCTION C4
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function SM1 (N)
!-----------------------------------------------------------------------
      IMPLICIT double precision (A - H, O - Y)
      integer:: N
      IF (N.EQ.0) THEN
         SM1 = 0D0
      ELSE
         SM1 = DCOS (N * RI * DPI) * DIS (RI, RI + 1, DPI * N) - DSIN ( &
         N * RI * DPI) * DIC (RI, RI + 1, N * DPI)
      ENDIF
      END FUNCTION SM1
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function CM1 (N)
!-----------------------------------------------------------------------
      IMPLICIT double precision (A - H, O - Y)
      integer:: N
      IF (N.EQ.0) THEN
         CM1 = DLOG ( (RI + 1) / RI)
      ELSE
         CM1 = DCOS (N * RI * DPI) * DIC (RI, RI + 1, N * DPI) + DSIN ( &
         N * RI * DPI) * DIS (RI, RI + 1, N * DPI)
      ENDIF
      END FUNCTION CM1
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function SM2 (N)
!-----------------------------------------------------------------------
      IMPLICIT double precision (A - H, O - Y)
      integer:: N
      IF (N.EQ.0) THEN
         SM2 = 0D0
      ELSE
         SM2 = DCOS (N * RI * DPI) * (DSIN (N * DPI * RI) / RI - DSIN ( &
         N * DPI * (RI + 1) ) / (RI + 1) + N * DPI * DIC (RI, RI + 1, N &
         * DPI) ) - DSIN (N * RI * DPI) * (DCOS (N * DPI * RI) / RI -   &
         DCOS (N * DPI * (RI + 1) ) / (RI + 1) - N * DPI * DIS (RI, RI +&
         1, N * DPI) )
      ENDIF
      END FUNCTION SM2
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function CM2 (N)
!-----------------------------------------------------------------------
      IMPLICIT double precision (A - H, O - Y)
      integer:: N
      IF (N.EQ.0) THEN
         CM2 = 1D0 / RI - 1D0 / (1 + RI)
      ELSE
         CM2 = DSIN (N * RI * DPI) * (DSIN (N * DPI * RI) / RI - DSIN ( &
         N * DPI * (RI + 1) ) / (RI + 1) + N * DPI * DIC (RI, RI + 1, N &
         * DPI) ) + DCOS (N * RI * DPI) * (DCOS (N * DPI * RI) / RI -   &
         DCOS (N * DPI * (RI + 1) ) / (RI + 1) - N * DPI * DIS (RI, RI +&
         1, N * DPI) )
      ENDIF
      END FUNCTION CM2
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function DIS (XMIN, XMAX, A)
!-----------------------------------------------------------------------
      IMPLICIT double precision (A - H, O - Y)
!
      X = DABS (A * XMAX)
      IF (X.LT.1) THEN
         DISMAX = SIS (X)
      ELSE
         DISMAX = SIA (X)
      ENDIF
      IF (A * XMAX.LT.0) DISMAX = - DISMAX
      X = DABS (A * XMIN)
      IF (X.LT.1) THEN
         DISMIN = SIS (X)
      ELSE
         DISMIN = SIA (X)
      ENDIF
      IF (A * XMIN.LT.0) DISMIN = - DISMIN
      DIS = DISMAX - DISMIN
      END FUNCTION DIS
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function SIA (X)
!-----------------------------------------------------------------------
      IMPLICIT double precision (A - H, O - Y)
      DIMENSION AF (4), BF (4), AG (4), BG (4)
      PARAMETER (DPI = 3.141592653589793D0)
      DATA AF / 38.027264D0, 265.187033D0, 335.677320D0, 38.102495D0 /  &
      BF / 40.021433D0, 322.624911D0, 570.236280D0, 157.105423D0 / AG / &
      42.242855D0, 302.757865D0, 352.018498D0, 21.821899D0 / BG /       &
      48.196927D0, 482.485984D0, 1114.978885D0, 449.690326D0 /
      F = (X**8 + AF (1) * X**6 + AF (2) * X**4 + AF (3) * X**2 + AF (4)&
      ) / (X**8 + BF (1) * X**6 + BF (2) * X**4 + BF (3) * X**2 + BF (4)&
      ) / X
      G = (X**8 + AG (1) * X**6 + AG (2) * X**4 + AG (3) * X**2 + AG (4)&
      ) / (X**8 + BG (1) * X**6 + BG (2) * X**4 + BG (3) * X**2 + BG (4)&
      ) / X**2
      SIA = DPI / 2D0 - F * DCOS (X) - G * DSIN (X)
      END FUNCTION SIA
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function SIS (X)
!-----------------------------------------------------------------------
      IMPLICIT double precision (A - H, O - Y)
      integer:: Jz, Jn, i, j
      GENAU = 1D-7
!
      SISO = X
      DO 200 I = 1, 200000
         JZ = 0
         JN = 0
         SISZ = 1D0
         SISN = 2D0 * I + 1D0
   10    DO 20 J = JZ + 1, 2 * I + 1
            SISZ = SISZ * X
            IF (SISZ.GT.1D20) THEN
               GOTO 30
            ENDIF
   20    END DO
   30    JZ = J
         DO 40 J = JN + 1, 2 * I + 1
            SISN = SISN * J
            IF (SISN.GT.1D20) THEN
               GOTO 50
            ENDIF
   40    END DO
   50    JN = J
         SISZ = SISZ / SISN
         SISN = 1
         IF ( (JZ.LT.2 * I + 1) .OR. (JN.LT.2 * I + 1) ) GOTO 10
         SISN = 2D0 * I + 1D0
         JN = 0
!
         SIS = SISO + ( - 1) **I * SISZ
         IF (I.GT.1) THEN
            IF (DABS (1D0 - SIS / SISO) .LE.GENAU) GOTO 300
         ENDIF
         SISO = SIS
  200 END DO
  300 END FUNCTION SIS
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function DIC (XMIN, XMAX, A)
!-----------------------------------------------------------------------
      IMPLICIT double precision (A - H, O - Y)
!
      X = DABS (A * XMAX)
      IF (X.LT.1) THEN
         DICMAX = CIS (X)
      ELSE
         DICMAX = CIA (X)
      ENDIF
      X = DABS (A * XMIN)
      IF (X.LT.1) THEN
         DICMIN = CIS (X)
      ELSE
         DICMIN = CIA (X)
      ENDIF
      DIC = DICMAX - DICMIN
      END FUNCTION DIC
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function CIA (X)
!-----------------------------------------------------------------------
      IMPLICIT double precision (A - H, O - Y)
      DIMENSION AF (4), BF (4), AG (4), BG (4)
      DATA AF / 38.027264D0, 265.187033D0, 335.677320D0, 38.102495D0 /  &
      BF / 40.021433D0, 322.624911D0, 570.236280D0, 157.105423D0 / AG / &
      42.242855D0, 302.757865D0, 352.018498D0, 21.821899D0 / BG /       &
      48.196927D0, 482.485984D0, 1114.978885D0, 449.690326D0 /
      F = (X**8 + AF (1) * X**6 + AF (2) * X**4 + AF (3) * X**2 + AF (4)&
      ) / (X**8 + BF (1) * X**6 + BF (2) * X**4 + BF (3) * X**2 + BF (4)&
      ) / X
      G = (X**8 + AG (1) * X**6 + AG (2) * X**4 + AG (3) * X**2 + AG (4)&
      ) / (X**8 + BG (1) * X**6 + BG (2) * X**4 + BG (3) * X**2 + BG (4)&
      ) / X**2
      CIA = F * DSIN (X) - G * DCOS (X)
      END FUNCTION CIA
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   double precision function CIS (X)
!-----------------------------------------------------------------------
      IMPLICIT double precision (A - H, O - Y)
      PARAMETER (C = 0.5772156649D0)
      integer:: Jz, Jn, i, j
      GENAU = 1D-7
!
      CISO = DLOG (X) + C
      DO 200 I = 1, 200000
         JZ = 0
         JN = 0
         CISZ = 1D0
         CISN = 2D0 * I
   10    DO 20 J = JZ + 1, 2 * I
            CISZ = CISZ * X
            IF (CISZ.GT.1D20) THEN
               GOTO 30
            ENDIF
   20    END DO
   30    JZ = J
         DO 40 J = JN + 1, 2 * I
            CISN = CISN * J
            IF (CISN.GT.1D20) THEN
               GOTO 50
            ENDIF
   40    END DO
   50    JN = J
         CISZ = CISZ / CISN
         CISN = 1D0
         IF ( (JZ.LT.2 * I) .OR. (JN.LT.2 * I) ) GOTO 10
!
         CIS = CISO + ( - 1) **I * CISZ
         IF (I.GT.1) THEN
            IF (DABS (1D0 - CIS / CISO) .LE.GENAU) GOTO 300
         ENDIF
         CISO = CIS
  200 END DO
  300 END FUNCTION CIS
end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
