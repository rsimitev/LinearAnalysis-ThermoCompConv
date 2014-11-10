************************************************************************
*       LOSUB.F for LO.C 
*
* Program to calc. lin. onset of conv. in rotating spherical shell.
*
* Rev 1.0 05/94                                               J.W.
* Rev 1.1 05/94 added loop for increment auf TAU              M.A.
* Rev 2.0 07/94 added loop for calcul. of minimal wavenumber
*              M (LCALC=3)                                    M.A.
* Rev 2.1 07/94 added part for calcul. of the eigenvector
*              at onset (LCALC=4)                             M.A.
* Rev 2.2 12/04/94 changed sign of drift c for consistency
*                 with LC.F                                   M.A.
* Rev 3.0 29/03/07  - double diffusive convection             R.S.
************************************************************************
* Parameters:
* RA=RAYLEIGH NUMBER, TA= TAYLOR NUMBER, PR= PRANTEL NUMBER,
* ETA=RATIO OF RADII, NT= TRUNCATION, M0=WAVE NUMBER,
* pL = Lewis number,  RAC=Rayleigh number due to concentration
* NE= SYMMETRIE PARAMETER, NE=0 : UNDEFINED SYMMETRIE,
* NE=2 : EQUATORIAL SYMMETRIE, NE=1 : EQUATORIAL ANTISYMMETRIE.
* DRIFT C IS DEFINED LIKE (PHI+C*T).
* Rev. 2.2: Drift is now def. as (phi-c*t).
*
* LCALC=1 : Eigenvalues are determined for const. parameters
* LCALC=2 : Onset determined for constant wavenumber M
*           (by searching root of grothrate in R, using pegasus.f).
* LCALC=3 : Onset determined by variing Rayleighnumber R and
*           wavenumber M. 
* LCALC=4 : Eigenvector determined for one set of parameters
*           at onset
*
* LO.F calculates R (crit. Rayleighn.) and Omega (and M) in the 
* range TAU=TTA to TAU=TTF.
* It only calculates the mode with minimal value of R.
*
* Main prg:    'lo.c'
* Subroutines: from modules
*              'losub.f'
*              'r.f'
*              'pegasus.f'
*              'imsl.f'   (some routines of the IMSL-Library)
*
* Compile and link:
*           make lo   (uses "makefile")
*
* Start:    lo inputfilename outputfilename
*
************************************************************************

******* Dimension of Matrix: *******************************************
*       NT=10:  NMAX=165  
*       NT=12:  NMAX=234  
*       NT=14:  NMAX=316  
*       NT=16:  NMAX=408  
*       NT=18:  NMAX=513  
*       NT=20:  NMAX=630
*       NT=22:  NMAX=759
#ifndef NMAX
#define NMAX 144
#endif

#ifdef __hpux
#define fdate fdate_
#endif

********* limit for CPU time in sec, for C routine clock(): *************
#ifdef __hpux
*******  medium queue:   4 h = 14000 sec
*******  long   queue:  16 h = 57600 sec
#define CPUTIMELIM         52000.0D0
#else
#define CPUTIMELIM         40000.0D0
#endif
#define CLOCKS_PER_SECOND  1000000.0D0
#define NCPUMAXVALUE       2147483647

******* DQS 3.0: Signal, when process exceeds CPU time limit: ********** 
******* (doesn't work actually, because signal is send only to 
*******  shell, not to the working process.
*******  Only the Bourneshell can catch the signal, but reacts
*******  not before one shell command is finished (especially the
*******  working process).
#define SIGRESTART  80

*******  defines from <sys/signal.h> *********************
#  define SIGABRT       6       /* Process abort signal */
#  define SIGUSR1       16      /* user defined signal 1 */
#  define SIGUSR2       17      /* user defined signal 2 */

******* Process return values (2...255) ********************************
#define NO_INFILE      100
#define ERR_IN_INFILE  101
#define DIM_TO_SMALL   102
#define ERR_WRT_OUTFILE 103
#define NO_RA_FOUND    120
#define START_NEXT_RUN 190
#define FINISHED       199

#define DOUBLEMAXVALUE  1.0D306
***********************************************************************
      SUBROUTINE losub(inputfile,outputfile)
      IMPLICIT REAL*8(A-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      COMMON/PAR2/ETA
      COMMON/DIM/NT,M0,ND,LMIN,LD,NE
      COMMON/NUM/TTA,TTF,TTSTEP,DRA,ABSE,RELE,NSMAX,LCALC
      CHARACTER*40 infile,outfile
      COMMON/FILES/infile,outfile
      CHARACTER*40 inputfile,outputfile
      REAL*8    RACI(3),OMI(3)
      INTEGER   M0I(3),LLI(3),LMINI(3)
      INTEGER nloop
      COMPLEX*16 ZEVEC(NMAX)
      COMMON/MODE/ZEVEC

      infile = inputfile
      outfile = outputfile

C-----Default values:
      NE = 2
      LCALC = 2
      RA = 4.D3
      Rconc = 4.D3
      RAOLD = 4.D3
      RconcOLD = 4.D3
      TAU = 100.0D0
      TAU1 = 100.0D0
      TAU2 = 100.0D0
      PR = 1.D-1
      pL = 1.D00
      ETA = 0.4D0
      NT = 3
      M0 = 6
      DRA = RA/10.0D0
      DRconc = Rconc/10.0D0
      ABSE = 0.D0
      RELE = 1.D-6
      NSMAX = 100
  
C-----INPUT:
      OPEN(15,FILE=infile,STATUS='OLD',ERR=10)
      GOTO 11
10    WRITE(*,*) 'LOSUB.F: Error while reading inputfile!'
      STOP   NO_INFILE
11    CONTINUE
      READ(15,'(A)',END=15) 
      READ(15,*,END=15) NE,LCALC
      READ(15,'(A)',END=15) 
      READ(15,*,END=15) RA,TTA,PR,ETA,pL,Rconc
      READ(15,'(A)',END=15) 
      READ(15,*,END=15) NT,M0
      READ(15,'(A)',END=15) 
      READ(15,*,END=15) DRA,ABSE,RELE,NSMAX
      READ(15,'(A)',END=15) 
      READ(15,*,END=15) TTSTEP,TTF
      CLOSE(15)
      GOTO 16
15    WRITE(*,*) 'Error in inputfile ',inputfile
      STOP ERR_IN_INFILE
16    CONTINUE

C-----LOSUB.F doesn't work for M=0 !!!!!!
      IF(M0.LT.1) THEN
        write(*,*) 'The code does not work for M0<1.', M0, ' --> 1'
        M0 = 1
      ENDIF

      TAU  = TTA
      TAU1 = TTA
      TAU2 = TTA
      TA = TAU*TAU
      RAOLD = RA
      RconcOLD = Rconc

C-----OUTPUT:
      OPEN(16,FILE=outfile,STATUS='UNKNOWN')
      IF(LCALC.GT.0 .AND. LCALC.LT.4 .or. LCALC.eq.5.or.LCALC.eq.6) THEN
       WRITE(16,*)  '### Output of Program lo.f Ver.2.1:        ###'
       WRITE(16,*)  '### Lin. Onset of Conv. via Galerkinmethod ###'
       WRITE(16,*)  '# P     ', PR, '#'
       WRITE(16,*)  '# Lewis ', pL, '#'
       WRITE(16,*)  '# TAU   ', TAU, '#'
       WRITE(16,*)  '# R     ', RA, '#'
       WRITE(16,*)  '# RC     ', Rconc, '#'
       WRITE(16,*)  '# ETA   ', ETA, '#'
       WRITE(16,*)  '# m     ', M0, '   #'
       WRITE(16,*)  '# cvar     TAU     #'
       WRITE(16,*)  '# NE    ', NE, '   #'
       WRITE(16,*)  '# TTA   ', TTA, '#'
       WRITE(16,*)  '# TTF   ', TTF, '#'
       WRITE(16,*)  '# TTSTEP ',TTSTEP, '#'
       WRITE(16,*)  '# NT    ', NT, '   #'
       WRITE(16,*)  '# see definition of LCALC for output. LCALC:',
     &            LCALC,'   #'
       WRITE(16,*)  '#                                      #'
      ENDIF
      
      CLOSE(16)
C
      RI=ETA/(1-ETA)
      RO=1.D0+RI
C
      IF(NE.EQ.0) THEN
C-- UNDEFINED SYMMETRIE:
         LMIN=M0
         LD=1
      ELSEIF(NE.EQ.1) THEN
C-- EQUATORIAL ANTISYMMETRIE (L+M ODD):
         LMIN=M0+1
         LD=2
      ELSEIF(NE.EQ.2) THEN
C-- EQUATORIAL SYMMETRIE (L+M EVEN);
         LMIN=M0
         LD=2
      ENDIF
C
      CALL DIMENSION
      write(*,*) 'DIMENSION OF MATRIX:',ND
C

      IF( LCALC.EQ.-1) THEN
c---LCALC=-1 : most basic case: find the most unstable growth rate at all parameters fixed
         if (.true.) then
c---  plot FN(x)
         Ra0=RA
         niter=(TTF-Ra0)/TTSTEP
         do nloop=0,niter
            RA=Ra0+nloop*TTSTEP
            GROR=FN(RA)
            WRITE(*,*) RA,GROR
            call open_file_at_end(16,outfile)
            WRITE(16,'(3D16.8)') RA,GROR
            close(16)
         enddo
         endif

         GROR=FN(RA)
         WRITE(*,*) 'R=',RA,' TAU=',TAU,' P=',PR,' M0=',M0,' eta=',ETA
         WRITE(*,* ) 'most unstable growth rate', GROR
         WRITE(*,* ) 'if growth rate < 0 then above onset'
         stop ' growth rate at fixed other parametres'

      else if( LCALC.EQ.0 ) THEN
c---LCALC=0 : Critical Ra, for constant other parameters
         call open_file_at_end(16,outfile)
         CALL PEGASUS(RA,DRA,ABSE,RELE,NSMAX,1,1,LL,RAC)
         GROR=FN(RAC)
         OMEGA=C*M0
         WRITE(*,*) 'TAU=',TAU,' P=',PR,' M0=',M0,' eta=',ETA
         WRITE(*,*) 'Lewis=',pL,' Rconc=',Rconc
         WRITE(*,* ) 'R_crit=',RAC, '  (growth rate =',GROR
         stop '             Ra_crit at fixed other parametres'

         else if( LCALC.EQ.1 ) THEN
C-----eigenvalues determined for this value of RA: ------------------
         GROR=FN(RA)
         RETURN

      ELSE IF (LCALC.EQ.4) THEN
C-----calculate the critical eigenvector. Print for plotting.
         open(16,FILE=outfile,STATUS='UNKNOWN')
C--      searching for zero grothrate:
         CALL PEGASUS(RA,DRA,ABSE,RELE,NSMAX,1,1,LL,RAC)
         GROR=FN(RAC)
         OMEGA=C*M0

         IF(LL.EQ.0) THEN
C--       print eigenvector
C---------Fileformat for outputfile:
C---------LST=0 formatted Wicht
C---------LST=1 formatted Hirsching
C---------LST=3 unformatted
          LST=0
          LCALC=1
          WRITE(16,'(2I2,'' LINEAR ONSET '')')  LST,LCALC
          NTH=0
          KTV=0
          KTH=0
          LTV=0
          LTH=0
          GRR=0.D0
          GRI=0.D0
          WRITE(16,'(I2,7I3,2D16.8,'' M0,TRUNC,LD,GROTH,DRIFT'')')
     &      M0,NT,NTH,KTV,KTH,LTV,LTH,LD,GRR,GRI
          NUDS=1
          PM=0.D0
          WRITE(16,
     &   '(I5,2D14.6,D9.2,D13.6,D9.2,'' I,TA,RA,PR,PM,E'')')
     &                                      NUDS,TA,RAC,PR,PM,ETA
          C0 = OMEGA/M0
          OMM=0.D0
          NUC=0
          NUOM=0
          MF=0
          WRITE(16,9100) C0,OMM,NUC,NUOM,MF
9100      FORMAT(2D17.10,3I4,'    C,OM, WHERE?,FLOQUET')

          LMAX=2*NT+M0-1
          I=0
          DO 3000 LI=LMIN,LMAX,LD
C           L for poloidal (v) field:
            LPI = LI
            NIMAX=DINT( DBLE(2*NT+1-LI+M0)/2 )
            DO 3000 NI=1,NIMAX
              IF( LST.EQ.0 ) THEN
               WRITE(16,9200) 'V',LPI,M0,NI,0,
     &         DREAL(ZEVEC(I+1)),DIMAG(ZEVEC(I+1)),0.D0,0.D0
              ELSEIF(LST.EQ.3) THEN
               WRITE(16,'(A,4I3,A,2F11.7,A)') ' ''V ''',LPI,M0,NI,0,' ',
     &         DREAL(ZEVEC(I+1)),DIMAG(ZEVEC(I+1)),' .0D+00 .0D+00 '
              ENDIF
             I=I+4
3000      CONTINUE
C         
          I=0
          DO 3200 LI=LMIN,LMAX,LD
C           L for toroidal (w) field:
            IF( NE.EQ.2 ) THEN
             LTI=LI+1
            ELSEIF( NE.EQ.1 ) THEN
             LTI=LI-1
            ELSEIF( NE.EQ.0 ) THEN
             LTI=LI
            ENDIF
            NIMAX=DINT( DBLE(2*NT+1-LI+M0)/2 )
            DO 3200 NI=1,NIMAX
              IF( LST.EQ.0 ) THEN
               WRITE(16,9200) 'W',LTI,M0,NI,0,
     &         DREAL(ZEVEC(I+3)),DIMAG(ZEVEC(I+3)),0.D0,0.D0
              ELSEIF(LST.EQ.3) THEN
               WRITE(16,'(A,4I3,A,2F11.7,A)') ' ''W ''',LTI,M0,NI,0,' ',
     &         DREAL(ZEVEC(I+3)),DIMAG(ZEVEC(I+3)),' .0D+00 .0D+00 '
              ENDIF
              I=I+4
3200      CONTINUE
C         
          I=0
          DO 3400 LI=LMIN,LMAX,LD
            NIMAX=DINT( DBLE(2*NT+1-LI+M0)/2 )
            DO 3400 NI=1,NIMAX
              IF( LST.EQ.0 ) THEN
               WRITE(16,9200) 'T',LI,M0,NI,0,
     &         DREAL(ZEVEC(I+2)),DIMAG(ZEVEC(I+2)),0.D0,0.D0
              ELSEIF(LST.EQ.3) THEN
               WRITE(16,'(A,4I3,A,2F11.7,A)') ' ''T ''',LI,M0,NI,0,' ',
     &         DREAL(ZEVEC(I+2)),DIMAG(ZEVEC(I+2)),' .0D+00 .0D+00 '
              ENDIF
             I=I+4
3400      CONTINUE
C         
          I=0
          DO 3600 LI=LMIN,LMAX,LD
            NIMAX=DINT( DBLE(2*NT+1-LI+M0)/2 )
            DO 3600 NI=1,NIMAX
              IF( LST.EQ.0 ) THEN
               WRITE(16,9200) 'G',LI,M0,NI,0,
     &         DREAL(ZEVEC(I+4)),DIMAG(ZEVEC(I+4)),0.D0,0.D0
              ELSEIF(LST.EQ.3) THEN
               WRITE(16,'(A,4I3,A,2F11.7,A)') ' ''G ''',LI,M0,NI,0,' ',
     &         DREAL(ZEVEC(I+4)),DIMAG(ZEVEC(I+4)),' .0D+00 .0D+00 '
              ENDIF
             I=I+4
3600      CONTINUE
C         
9200            FORMAT(1X,A1,4I3,4D16.8)
         ELSE
          WRITE(16,*) 'NO CRITICAL RAYLEIGH NUMBER FOUND.'
          STOP NO_RA_FOUND
       ENDIF
      endif

      if(( LCALC.EQ.2 ).or.( LCALC.EQ.3 )) then
C     LCALC-s which requiree increment of TAU:
      nloop=0
20    continue
      TA = TAU*TAU

      IF( LCALC.EQ.2 ) THEN
C-----searching for zero grothrate by varying RA: ---------------------- 
         call open_file_at_end(16,outfile)
         CALL PEGASUS(RA,DRA,ABSE,RELE,NSMAX,1,1,LL,RAC)
         GROR=FN(RAC)
         OMEGA=C*M0
         IF(LL.EQ.0) THEN
            WRITE(16,'(1P,3G13.6,I4)') TAU,RAC,OMEGA,M0
            write( *,'(1P,3G13.6,I4)') TAU,RAC,OMEGA,M0
            NTRYCOUNT = 0
         ELSE IF(NTRYCOUNT.GE.3) THEN 
            WRITE(16,*) 'NO CRITICAL RAYLEIGH NUMBER FOUND.'
            STOP NO_RA_FOUND
         ELSE
            print*, LL
            write(*,*) 'NO CRIT. RAYLEIGH NUMBER FOUND. Trying again.'
            NTRYCOUNT = NTRYCOUNT + 1
         ENDIF

      ELSE IF( LCALC.EQ.3 ) THEN
C-----searching for zero grothrate by variing RA and M0: ---------------
         call open_file_at_end(16,outfile)
         IF(M0.GT.1) THEN
           M0I(1) = M0-1
           LMINI(1) = LMIN-1
         ELSE
           M0I(1) = M0 
           LMINI(1) = LMIN
         ENDIF
         M0I(2) = M0
         M0I(3) = M0+1 
         LMINI(2) = LMIN
         LMINI(3) = LMIN+1
         DO 100 II=1,3
           M0  = M0I(II)
           LMIN = LMINI(II)
           call dimension
           DRA = RA/10.0D0
           CALL PEGASUS(RA,DRA,ABSE,RELE,NSMAX,1,0,LL,RAC)
           GROR=FN(RAC)
           OMI(II) = C*M0
           RACI(II)= RAC
           LLI(II) = LL
           write( *,'(1P,3G13.6,I4,A3,2G13.6)') 
     & TAU,RAC,OMI(II),M0,' | ',
     & RAC*(1/(1-eta))**4*(2/TAU),(OMI(II)*2.0/TAU)
100      CONTINUE

         INDEX = 0
         IF(LLI(1).EQ.0 .AND. LLI(2).EQ.0) THEN
           IF( RACI(1).LT.RACI(2)) THEN
               INDEX = 1
           ELSE
               INDEX = 2
           ENDIF
         ELSEIF(LLI(1).EQ.0) THEN
               INDEX = 1
         ELSEIF(LLI(2).EQ.0) THEN 
               INDEX = 2
         ENDIF
         IF(LLI(3).EQ.0 .AND. INDEX.GT.0) THEN
           IF( RACI(3).LT.RACI(INDEX) ) THEN
             INDEX = 3
           ENDIF
         ELSEIF(LLI(3).EQ.0) THEN
             INDEX = 3
         ENDIF
          
         IF( INDEX.GT.0 ) THEN
            RAC=RACI(INDEX)
            OMEGA=OMI(INDEX)
            M0=M0I(INDEX)
            LMIN=LMINI(INDEX)
           write( 16,'(1P,3G13.6,I4)') 
     & TAU,RAC,OMEGA,M0
           write( *,'(1P,3G13.6,I4,A3,2G13.6)') 
     & TAU,RAC,OMEGA,M0,' | ',
     & RAC*(1/(1-eta))**4*(2/TAU),(OMEGA*2.0/TAU)
            write(*,*)
            NTRYCOUNT = 0
         ELSE IF(NTRYCOUNT.GE.3) THEN 
            WRITE(16,*) 'NO CRITICAL RAYLEIGH NUMBER FOUND.'
            STOP NO_RA_FOUND
         ELSE
            write(*,*) 'NO CRIT. RAYLEIGH NUMBER FOUND. Trying again.'
            NTRYCOUNT = NTRYCOUNT + 1
         ENDIF
      ENDIF 
C-----------------------------------------------------------------------
      CLOSE(16)

C--   increment TAU:
      TAU0 = TAU1
      TAU1 = TAU
      IF( DABS(TTSTEP) .LT. DABS(TAU*0.1D0) ) THEN
         TAU   = TAU1 + TTSTEP
      ELSE
         TAU = TAU1 + TTSTEP/DABS(TTSTEP)*TAU1*0.1D0
      ENDIF 
C--   interpolate new startingvalue for RA:
      IF(TAU1.NE.TAU0) THEN
        RA = RAC + (RAC-RAOLD)/(TAU1-TAU0)*(TAU-TAU1)
      ELSE
        RA = 2.0D0*RAC - RAOLD 
      ENDIF
      RAOLD = RAC

C--   endvalue of TAU reached?
      IF ( ((TAU.GT.TTF).AND.(TTF.GT.TTA)).OR.
     &      ((TAU.LT.TTF).AND.(TTF.LT.TTA)) )  THEN
c      WRITE(*,*) 'LOSUB.F: finished at ',fdate()
      STOP FINISHED
      ENDIF

C-----count loops:
      nloop=nloop+1
      GOTO 20
C--   End of tau Loop
      end if

      if (LCALC.EQ.5) THEN
C-----vary m and calculate critical R at fixed P, tau, eta.
            M0A=M0
         do M0=M0A,TTF, INT(TTSTEP)

            IF(NE.EQ.0) THEN
C--   UNDEFINED SYMMETRIE:
               LMIN=M0
               LD=1
            ELSEIF(NE.EQ.1) THEN
C--   EQUATORIAL ANTISYMMETRIE (L+M ODD):
               LMIN=M0+1
               LD=2
            ELSEIF(NE.EQ.2) THEN
C--   EQUATORIAL SYMMETRIE (L+M EVEN);
               LMIN=M0
               LD=2
            ENDIF
            call dimension
            DRA = RA/10.0D0
            CALL PEGASUS(RA,DRA,ABSE,RELE,NSMAX,1,0,LL,RAC)
            WRITE(*,* ) M0,RAC
            call open_file_at_end(16,outfile)
            WRITE(16,*) M0,RAC
            close(16)
         enddo

      else if (LCALC.EQ.6) THEN
C-----vary Le and calculate critical R at fixed P, tau, eta, M
         pL0=pL
         niter=(TTF-pL0)/TTSTEP
         do nloop=0,niter
            pL=pL0+nloop*TTSTEP
            DRA = RA/10.0D0
            CALL PEGASUS(RA,DRA,ABSE,RELE,NSMAX,1,0,LL,RAC)
            GROR=FN(RAC)
            WRITE(*,*) pL,RAC,GROR
            call open_file_at_end(16,outfile)
            WRITE(16,'(3D16.8)') pL,RAC,GROR
            close(16)
         enddo
      end if

      END

************************************************************************
      FUNCTION FN(RA)
************************************************************************
C-- FUNCTION CALLED BY PEGASUS, FN=GROTHRATE.
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      DIMENSION ZA(NMAX,NMAX),ZB(NMAX,NMAX)
      DIMENSION ZACOPY(NMAX,NMAX),ZBCOPY(NMAX,NMAX)
      DIMENSION ZEW(NMAX),ZEWA(NMAX),ZEWB(NMAX)
      DIMENSION ZEVALL(NMAX,NMAX),ZEVEC(NMAX)
      COMMON/PAR/TAU,RAL,PR,RI,C,pL,Rconc
      COMMON/DIM/NT,M0,ND,LMIN,LD,NE
      COMMON/MODE/ZEVEC
      COMMON/NUM/TTA,TTF,TTSTEP,DRA,ABSE,RELE,NSMAX,LCALC

      RAL=RA

C-- MAT SETS THE COMPLEX MATRICES ZA AND ZB SETTING OF MATRIX:
      CALL MAT(ZA,ZB,NMAX)

      DO 50 I=1,NMAX
       ZEWA(I)=DCMPLX(0D0,0D0)
       ZEWB(I)=DCMPLX(0D0,0D0)
       DO 50 J=1,NMAX
       ZACOPY(I,J)=DCMPLX(0D0,0D0)
       ZBCOPY(I,J)=DCMPLX(0D0,0D0)
50    CONTINUE
C
C-- DGVLCG IS AN IMSL ROUTINE THAT CALCULATES THE EIGENVALUES ZEWA/ZEWB:
C-- Problem: generalized complex eigensystem A*x = lam*B*x
C--          Input:  A: ZA, B: ZB
C--          Output: complex ZEWA(ND), ZEWB(ND)
C      CALL DGVLCG(ND,ZA,NMAX,ZB,NMAX,ZEWA,ZEWB)
      IF(LCALC.EQ.4) THEN
        CALL DG2CCG (ND,ZA,NMAX,ZB,NMAX,ZEWA,ZEWB,ZEVALL,
     &                           NMAX,ZACOPY,ZBCOPY)
      ELSE
       CALL DG2LCG(ND,ZA,NMAX,ZB,NMAX,ZEWA,ZEWB,ZACOPY,ZBCOPY)
      ENDIF
C
      DO 100 I=1,ND
100   ZEW(I)=ZEWA(I)/ZEWB(I)
      IMIN=1
C---- search for lowest imaginary part:
      EWMIN=DIMAG(ZEW(IMIN))
      DO 200 I=2,ND
         IF( DIMAG(ZEW(I)).LT.EWMIN ) THEN 
            IMIN=I
            EWMIN=DIMAG(ZEW(IMIN))
         ENDIF
200   CONTINUE
    
c      IF( EWMIN.LT.0.D0 ) WRITE(*,*) 'CONVECTION REGIME !',RAL
      C=-DREAL(ZEW(IMIN))/M0
C
      IF(LCALC.EQ.1) THEN
C---- sort eigenvalues:
       write(*,*) '     Frequ.(exp(+iwt))   -Grothrate  '
       DO 320 I=1,ND
        DO 310 J=I,ND
         IF( DIMAG(ZEW(J)).LT.DIMAG(ZEW(I)) ) THEN
            ZSAVE = ZEW(J)
            ZEW(J) = ZEW(I)
            ZEW(I) = ZSAVE
            DO 305 K=1,ND
             ZEVEC(K)    = ZEVALL(K,J)
             ZEVALL(K,J) = ZEVALL(K,I)
             ZEVALL(K,I) = ZEVEC(K)
305         continue
         ENDIF 
310     continue
320    WRITE(*,'(I4,2D16.6)') I,ZEW(I)     
      ENDIF

C------save crit. eigenvector:
      IF(LCALC.EQ.4) THEN
        DO 400 J=1,ND
         ZEVEC(J) = ZEVALL(J,IMIN)
400     CONTINUE
      ENDIF

      FN=EWMIN
c      if ( fn.lt.0.0) then
c         print*,  ' Convection onset found :) :  Ra= ...  ',RAL
c      else
c         print*,  ' Stability :( :  Ra= ...  ',RAL
c      endif         

C
      END
C
C------------------------------------------------------------------------
C
C
*************************************************************************
      SUBROUTINE MAT(ZA,ZB,NDIM)
*************************************************************************
C-- SETS THE COMPLEX MATRICES ZA AND ZB.
C------------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
C
      DIMENSION ZA(NDIM,*),ZB(NDIM,*)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      COMMON/DIM/NT,M0,ND,LMIN,LD,NE
C
      DO 200 J=1,NDIM
       DO 200 I=1,NDIM
         ZA(I,J)=DCMPLX(0D0,0D0)
         ZB(I,J)=DCMPLX(0D0,0D0)
200   CONTINUE
C
      I=0
      LMAX=2*NT+M0-1
C
      DO 4000 LI=LMIN,LMAX,LD
         LPI=LI
C
C        Determine L for toroidal (w) field:
         IF( NE.EQ.2 ) THEN 
            LTI=LI+1
         ELSEIF( NE.EQ.1 ) THEN 
            LTI=LI-1
         ELSEIF( NE.EQ.0 ) THEN
            LTI=LI
         ENDIF
C
         NIMAX=DINT( DBLE(2*NT+1-LI+M0)/2 )
C
         DO 3000 NI=1,NIMAX
C
            J=0
            DO 2000 LJ=LMIN,LMAX,LD
               LPJ=LJ
C
               IF( NE.EQ.2 ) THEN
                  LTJ=LJ+1
               ELSEIF( NE.EQ.1 ) THEN 
                  LTJ=LJ-1
               ELSEIF( NE.EQ.0 ) THEN
                  LTJ=LJ
               ENDIF
C
               NJMAX=DINT( DBLE(2*NT+1-LJ+M0)/2 )
C
C  ******************** I: Equation (Line) ******************
C  ******************** J: Variable (Column) ****************
C  ******************** I+1: v (poloidal)  ******************
C  ******************** I+2: theta         ******************
C  ******************** I+3: w (toroidal)  ******************
C-new****************** I+4: gamma (concentration) ********** 
               DO 1000 NJ=1,NJMAX
C
               IF(J+3.GT.NDIM .OR. I+3.GT.NDIM) THEN
                  write(*,*) 'MAT(): NDIM too small.'
                  stop
               ENDIF
C
                  IF( LI.EQ.LJ ) THEN
      	             ZB(I+1,J+1)=DCMPLX(0.D0,-DIII2(NI,NJ,LPI,1))
                     ZA(I+1,J+1)=
     &			DCMPLX(DIII1(NI,NJ,LPI),DIII3(NI,NJ,LPI,1))
                     ZA(I+1,J+2)=DCMPLX(DIII5(NI,NJ,LPI),0.D0)
c--- concentration driving
                     ZA(I+1,J+4)=DCMPLX(DIII5conc(NI,NJ,LPI),0.D0)
c
      	             ZB(I+2,J+2)=DCMPLX(0.D0,-DI1(NI,NJ,1))
                     ZA(I+2,J+1)=DCMPLX(DI3(NI,NJ,LPI),0.D0)
                     ZA(I+2,J+2)=DCMPLX(DI2(NI,NJ,LPI),0.D0)
      	             ZB(I+3,J+3)=DCMPLX(0.D0,-DII2(NI,NJ,LTI,1))
                     ZA(I+3,J+3)=
     &			DCMPLX(DII1(NI,NJ,LTI),DII3(NI,NJ,1))
c--- concentration equation
      	             ZB(I+4,J+4)=DCMPLX(0.D0,-DI1(NI,NJ,1))
                     ZA(I+4,J+1)=DCMPLX(DI3(NI,NJ,LPI),0.D0)
                     ZA(I+4,J+4)=DCMPLX(1.D0/pL * DI2(NI,NJ,LPI),0.D0)
                  ENDIF
                  IF( LPI.EQ.LTJ+1 ) THEN
      	             ZA(I+1,J+3)=DCMPLX(DIII4A(NI,NJ,LPI,1),0.D0)      
                  ELSEIF( LPI.EQ.LTJ-1 ) THEN
      	             ZA(I+1,J+3)=DCMPLX(DIII4B(NI,NJ,LPI,1),0.D0)      
      	          ENDIF
                  IF( LTI.EQ.LPJ+1 ) THEN
      	             ZA(I+3,J+1)=DCMPLX(DII4A(NI,NJ,LTI,1),0.D0)      
                  ELSEIF( LTI.EQ.LPJ-1 ) THEN
      	             ZA(I+3,J+1)=DCMPLX(DII4B(NI,NJ,LTI,1),0.D0)      
                  ENDIF
C
c                  print*, LI,NI,LJ,NJ,I,J

                  J=J+4
1000           CONTINUE
C
2000        CONTINUE
C
      	    I=I+4
3000     CONTINUE
C
4000  CONTINUE
C
c      print*, I, J
c      stop
      END
C
C-----------------------------------------------------------------------
C
C
*************************************************************************
C-- GALERKIN TERMS:
*************************************************************************
C
      FUNCTION DI1(N1,N2,NU1)
C----- HEAT EQUATION, TIME DERIVATIVE
         IMPLICIT REAL*8(A-H,O-Z)
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
         DI1=PR*NU1*R('SS ',2,N1,N2,0)
      END
C
      FUNCTION DI2(N1,N2,L1)
C----- HEAT EQUATION , DISSIPATION
         IMPLICIT REAL*8(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
         DI2=N2**2*DPI**2*R('SS ',2,N1,N2,0) - 
     -           2*N2*DPI*R('SC ',1,N1,N2,0) +
     +	           DL(L1)*R('SS ',0,N1,N2,0)  
      END
C
      FUNCTION DI3(N1,N2,L1)
C----- HEAT EQUATION , SOURCE 
         IMPLICIT REAL*8(A-H,O-Z)
         DI3 =-DL(L1)*R('SS ',2,N1,N2,0)
      END
C
      FUNCTION DII1(N1,N2,L1)
C----- TOROIDAL EQUATION , DISSIPATION
         IMPLICIT REAL*8(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
      	 DII1=DL(L1)*( (N2-1)**2*DPI**2*R('CC ',4,N1-1,N2-1,0) +
     +                     4*(N2-1)*DPI*R('CS ',3,N1-1,N2-1,0) + 
     +                       (DL(L1)-2)*R('CC ',2,N1-1,N2-1,0) )
      END
C
      FUNCTION DII2(N1,N2,L1,NU1)
C----- TOROIDAL EQUATION , TIME DERIVATIVE
         IMPLICIT REAL*8(A-H,O-Z)
      	 DII2=NU1*DL(L1)*R('CC ',4,N1-1,N2-1,0)
      END
C
      FUNCTION DII3(N1,N2,NU1)
C----- TOROIDAL EQUATION , CORRIOLIS
         IMPLICIT REAL*8(A-H,O-Z)
         COMMON/DIM/NT,M0,ND,LMIN,LD,NE
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
         DII3=-TAU*NU1*M0*R('CC ',4,N1-1,N2-1,0) 
      END
C
      FUNCTION DII4A(N1,N2,L1,NU1)
C----- TOROIADL EQUATION , Q-TERM 1 (L1=L3+1)
         IMPLICIT REAL*8(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
         COMMON/DIM/NT,M0,ND,LMIN,LD,NE
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      	 DII4A= TAU * 
     *        DSQRT( DBLE(L1-NU1*M0)*(L1+NU1*M0)/(2*L1-1)/(2*L1+1) ) * 
     *	             (    (L1**2-1)*(L1-1)*R('CS ',2,N1-1,N2,0)  - 
     -		      (L1+1)*(L1-1)*N2*DPI*R('CC ',3,N1-1,N2,0)    ) 
      END
C
      FUNCTION DII4B(N1,N2,L1,NU1)
C----- TOROIADL EQUATION , Q-TERM 1 (L1=L3-1)
         IMPLICIT REAL*8(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
         COMMON/DIM/NT,M0,ND,LMIN,LD,NE
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      	 DII4B= TAU * 
     *	  DSQRT( DBLE(L1-NU1*M0+1)*(L1+NU1*M0+1)/(2*L1+1)/(2*L1+3) ) * 
     *	             ( (1-(L1+1)**2)*(L1+2)*R('CS ',2,N1-1,N2,0)  - 
     -		           L1*(L1+2)*N2*DPI*R('CC ',3,N1-1,N2,0)  ) 
      END
C
      FUNCTION DIII1(N1,N2,L1)
C----- POLOIDAL EQUOATION , DISSIPATION
         IMPLICIT REAL*8(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
         DIII1=DL(L1)* ( N2**4*DPI**4*R('SS ',2,N1,N2,0) - 
     -	               4*N2**3*DPI**3*R('SC ',1,N1,N2,0) +
     +		2*DL(L1)*N2**2*DPI**2*R('SS ',0,N1,N2,0) +
     +           (DL(L1)**2-2*DL(L1))*R('SS ',-2,N1,N2,0) )
      END
C
      FUNCTION DIII2(N1,N2,L1,NU1)
C----- POLOIDAL EQUATION , TIME DERIVATIVE
         IMPLICIT REAL*8(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
         DIII2=        -NU1*DL(L1)*(-N2**2*DPI**2*R('SS ',2,N1,N2,0) +
     +         2*N2*DPI*R('SC ',1,N1,N2,0)-DL(L1)*R('SS ',0,N1,N2,0)  ) 
      END
C
      FUNCTION DIII3(N1,N2,L1,NU1)
C----- POLOIDAL EQUATION , CORRIOLIS
         IMPLICIT REAL*8(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
         COMMON/DIM/NT,M0,ND,LMIN,LD,NE
      	 DIII3= TAU*NU1*M0*( -N2**2*DPI**2*R('SS ',2,N1,N2,0) +
     +	                          2*N2*DPI*R('SC ',1,N1,N2,0) -
     -			            DL(L1)*R('SS ',0,N1,N2,0)  ) 
      END
C
      FUNCTION DIII4A(N1,N2,L1,NU1)
C----- POLOIDAL EUQUATION , Q-TERM 1 (L1=L3+1)
         IMPLICIT REAL*8(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
         COMMON/DIM/NT,M0,ND,LMIN,LD,NE
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      	 DIII4A= TAU * 
     *		DSQRT(DBLE(L1-M0*NU1)*(L1+M0*NU1)/(2*L1-1)/(2*L1+1)) *
     *          (     (L1*(L1-1)-2)*(L1-1)*R('SC ',2,N1,N2-1,0) +
     +	          (L1+1)*(L1-1)*(N2-1)*DPI*R('SS ',3,N1,N2-1,0)    )
      END
C
      FUNCTION DIII4B(N1,N2,L1,NU1)
C----- POLOIDAL EQUATION , Q-TERM 2 (L1=L3-1)
         IMPLICIT REAL*8(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
         COMMON/DIM/NT,M0,ND,LMIN,LD,NE
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
         DIII4B= TAU * 
     *	     DSQRT( DBLE(L1-M0*NU1+1)*(L1+M0*NU1+1)/(2*L1+1)/(2*L1+3) )*
     *	           (    (L1+2)*(2-(L1+1)*(L1+2))*R('SC ',2,N1,N2-1,0) + 
     +  	     L1*(L1+2)*(N2-1)*DPI*R('SS ',3,N1,N2-1,0)   )
      END
C
      FUNCTION DIII5(N1,N2,L1)
C----- POLOIDAL EQUATION , 
         IMPLICIT REAL*8(A-H,O-Z)
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
         DIII5=-RA*DL(L1)*R('SS ',2,N1,N2,0)
      END

      FUNCTION DIII5conc(N1,N2,L1)
C----- POLOIDAL EQUATION , 
         IMPLICIT REAL*8(A-H,O-Z)
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
         DIII5conc=-Rconc*DL(L1)*R('SS ',2,N1,N2,0)
      END
C
C--------------------------------------------------------------------------
C
C
***************************************************************************
*-- SUBROUTINES:
***************************************************************************
C
      FUNCTION DL(L)
         IMPLICIT REAL*8(A-H,O-Z)
         DL = DBLE(L*(L+1))
      END
C
C-- DUMMY FUNCTIONS:
C
      FUNCTION F2(X)
         F2=X
      END
C
      FUNCTION F3(X)
         F3=X
      END
C
      FUNCTION F4(X)
         F4=X
      END

****************************************************************************
      SUBROUTINE DIMENSION
****************************************************************************
      IMPLICIT REAL*8(A-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      COMMON/DIM/NT,M0,ND,LMIN,LD,NE
C-- DETERMINATION OF DIMENSION:
C-- for each value of L the number of possible N-values is added
      ND=0 
c         print*, "Triangular truncation (2.12)"
c         print*, LMIN, "...", 2*NT+M0-1,LD
      DO 18 L=LMIN,2*NT+M0-1,LD
c         print*, L, 1, "...", INT( DBLE(2*NT+1-L+M0)/2 )
cccccccccc18    ND=ND+3*DINT( DBLE(2*NT+1-L+M0)/2 )
18    ND=ND+4*DINT( DBLE(2*NT+1-L+M0)/2 )
C
      IF(ND.GT.NMAX) THEN
         WRITE(*,*) 'DIMENSION OF MATRIX TOO SMALL:',ND,'>',NMAX
         STOP DIM_TO_SMALL
      ENDIF

      END
*************************************************************************
      SUBROUTINE open_file_at_end(NHANDLE,filename)
*************************************************************************
*     opens file <filename> and puts the filepointer at EOF
*************************************************************************
      INTEGER NHANDLE
      INTEGER*4 II
      CHARACTER*40 filename

      OPEN(NHANDLE,FILE=filename,STATUS='OLD',ERR=990)
      DO 100 II=1,1000000000
         READ(NHANDLE,'(A)',END=999,ERR=990)
100   CONTINUE
      GOTO 999
990   WRITE(*,*) 'Error reading ',filename
      STOP ERR_WRT_OUTFILE
999   CONTINUE
      END


*************************************************************************
      SUBROUTINE abort
*     this function should replace the internal FORTRAN routine
*     'abort' which is e.g. called from the IMSL routines
*************************************************************************
      IMPLICIT REAL*8(A-H,O-Y)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      COMMON/NUM/TTA,TTF,TTSTEP,DRA,ABSE,RELE,NSMAX,LCALC

      write(*,*) 'LOSUB.F: abort(fortran) called.'
C-----increment TAU:
      IF( DABS(TTSTEP) .LT. DABS(TAU*0.1D0) ) THEN
         TAU   = TAU + TTSTEP
      ELSE
         TAU = TAU + TAU*0.1D0
      ENDIF 
      
C--   endvalue of TAU reached?
      IF ( ((TAU.GT.TTF).AND.(TTF.GT.TTA)).OR.
     &      ((TAU.LT.TTF).AND.(TTF.LT.TTA)) )  THEN
c      WRITE(*,*) 'LOSUB.F: finished at ',fdate()
      STOP FINISHED
      ENDIF

      call ende(SIGRESTART)
      RETURN
      END

*************************************************************************
      SUBROUTINE ende(II)
*     ende() is called explicit in the program or implicit when the
*     process receives a signal.
*     II indicates the behaviour of ende() or the received signal number.
*************************************************************************
      IMPLICIT REAL*8(A-H,O-Y)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      COMMON/PAR2/ETA
      COMMON/DIM/NT,M0,ND,LMIN,LD,NE
      COMMON/NUM/TTA,TTF,TTSTEP,DRA,ABSE,RELE,NSMAX,LCALC
      CHARACTER*40 infile,outfile
      COMMON/FILES/infile,outfile
      INTEGER II,NHANDLE
      LOGICAL LOP
c      CHARACTER*24 fdate
c      EXTERNAL     fdate

      write(*,*) 'ende(',II,'):  TAU = ',TAU

      IF(II.EQ.SIGUSR1) THEN
c       write(*,'(A,A)') ' ende(): received signal SIGUSR1 at ',fdate()
      ELSEIF(II.EQ.SIGABRT) THEN
c       write(*,'(A,A)') ' ende(): received signal SIGABRT at ',fdate()
      ENDIF

      INQUIRE(FILE=outfile,OPENED=LOP,NUMBER=NHANDLE)
      IF(LOP) THEN
        CLOSE(16)
      ENDIF

      IF(II.EQ.SIGRESTART .OR. II.EQ.SIGUSR1) THEN
         WRITE(*,*) 'LOSUB.F: Terminating this run and starting next.'
         OPEN(15,FILE=infile,STATUS='UNKNOWN')
         WRITE(15,*) ' NE (0/1/2) | LCALC (1/2/3/4) |'
         WRITE(15,'(A,2I12)') ' ',NE,LCALC
         WRITE(15,*) '|  RAYLEIGH  |  TAU     |  PRANTEL  |  ETA  |'
         WRITE(15,'(1P,G13.6,A,G13.6,A,G13.6,A,G13.6)') 
     &         RA,' ',TAU,' ',PR,' ',ETA
         WRITE(15,*) '|   NTRUNC (>=1) | MODE |'
         WRITE(15,'(A,2I12)') ' ',NT,M0
         WRITE(15,*) '|   DRA   | ABSERR  |  RELERR  | NMAX |'
         WRITE(15,'(1PG12.6,A,1PG11.5,A,1PG11.5,A,I4)') 
     &           DRA,' ',ABSE,' ',RELE,' ',NSMAX
         WRITE(15,*) '|   TAU_STEP | TAU_END'
         WRITE(15,'(1P,2G11.4)') TTSTEP,TTF
         CLOSE(15)
C--------return value START_NEXT_RUN to the shell script: 
c         call exit(START_NEXT_RUN) 
      ELSE
c         call exit(1)
      ENDIF

      END
