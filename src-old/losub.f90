!***********************************************************************
!       LOSUB.F for LO.C
!
! Program to calc. lin. onset of conv. in rotating spherical shell.
!
! Rev 1.0 05/94                                               J.W.
! Rev 1.1 05/94 added loop for increment auf TAU              M.A.
! Rev 2.0 07/94 added loop for calcul. of minimal wavenumber
!              M (LCALC=3)                                    M.A.
! Rev 2.1 07/94 added part for calcul. of the eigenvector
!              at onset (LCALC=4)                             M.A.
! Rev 2.2 12/04/94 changed sign of drift c for consistency
!                 with LC.F                                   M.A.
! Rev 3.0 29/03/07  - double diffusive convection             R.S.
!***********************************************************************
! Parameters:
! RA=RAYLEIGH NUMBER, TA= TAYLOR NUMBER, PR= PRANTEL NUMBER,
! ETA=RATIO OF RADII, NT= TRUNCATION, M0=WAVE NUMBER,
! pL = Lewis number,  RAC=Rayleigh number due to concentration
! NE= SYMMETRIE PARAMETER, NE=0 : UNDEFINED SYMMETRIE,
! NE=2 : EQUATORIAL SYMMETRIE, NE=1 : EQUATORIAL ANTISYMMETRIE.
! DRIFT C IS DEFINED LIKE (PHI+C*T).
! Rev. 2.2: Drift is now def. as (phi-c*t).
!
! LCALC=1 : Eigenvalues are determined for const. parameters
! LCALC=2 : Onset determined for constant wavenumber M
!           (by searching root of grothrate in R, using pegasus.f).
! LCALC=3 : Onset determined by variing Rayleigh number R and
!           wavenumber M.
! LCALC=4 : Eigenvector determined for one set of parameters
!           at onset
!
! LO.F calculates R (crit. Rayleighn.) and Omega (and M) in the
! range TAU=TTA to TAU=TTF.
! It only calculates the mode with minimal value of R.
!
! Main prg:    'lo.c'
! Subroutines: from modules
!              'losub.f'
!              'r.f'
!              'pegasus.f'
!              'imsl.f'   (some routines of the IMSL-Library)
!
! Compile and link:
!           make lo   (uses "makefile")
!
! Start:    lo inputfilename outputfilename
!
!***********************************************************************

!****** Dimension of Matrix: *******************************************
!       NT=10:  NMAX=165
!       NT=12:  NMAX=234
!       NT=14:  NMAX=316
!       NT=16:  NMAX=408
!       NT=18:  NMAX=513
!       NT=20:  NMAX=630
!       NT=22:  NMAX=759
#ifndef NMAX
#define NMAX 144
#endif

#ifdef __hpux
#define fdate fdate_
#endif

!******** limit for CPU time in sec, for C routine clock(): ************
#ifdef __hpux
!******  medium queue:   4 h = 14000 sec
!******  long   queue:  16 h = 57600 sec
#define CPUTIMELIM         52000.0D0
#else
#define CPUTIMELIM         40000.0D0
#endif
#define CLOCKS_PER_SECOND  1000000.0D0
#define NCPUMAXVALUE       2147483647

!****** DQS 3.0: Signal, when process exceeds CPU time limit: **********
!****** (doesn't work actually, because signal is send only to
!******  shell, not to the working process.
!******  Only the Bourneshell can catch the signal, but reacts
!******  not before one shell command is finished (especially the
!******  working process).
#define SIGRESTART  80

!******  defines from <sys/signal.h> *********************
#  define SIGABRT       6       /* Process abort signal */
#  define SIGUSR1       16      /* user defined signal 1 */
#  define SIGUSR2       17      /* user defined signal 2 */

!****** Process return values (2...255) ********************************
#define NO_INFILE      100
#define ERR_IN_INFILE  101
#define DIM_TO_SMALL   102
#define ERR_WRT_OUTFILE 103
#define NO_RA_FOUND    120
#define START_NEXT_RUN 190
#define FINISHED       199

#define DOUBLEMAXVALUE  1.0D306
!**********************************************************************
      SUBROUTINE setDefaults ()
      IMPLICIT double precision (A - H, O - Y)
      IMPLICIT complex (8) (Z)
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      COMMON / PAR2 / ETA
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      COMMON / NUM / TTA, TTF, TTSTEP, DRA, ABSE, RELE, NSMAX, LCALC
      NE = 2
      RA = 4.D3
      Rconc = 4.D3
      TAU = 100.0D0
      PR = 1.D-1
      pL = 1.D00
      ETA = 0.4D0
      NT = 3
      M0 = 6
      DRA = RA / 10.0D0
      ABSE = 0.D0
      RELE = 1.D-6
      NSMAX = 100
      LCALC = 2
      endsubroutine

!**********************************************************************
      SUBROUTINE readInputFile (inputfile)
      IMPLICIT double precision (A - H, O - Y)
      IMPLICIT complex (8) (Z)
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      COMMON / PAR2 / ETA
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      COMMON / NUM / TTA, TTF, TTSTEP, DRA, ABSE, RELE, NSMAX, LCALC
      CHARACTER (len = * ) inputfile
      OPEN (15, FILE = inputfile, STATUS = 'OLD', ERR = 10)
      GOTO 11
   10 WRITE ( * , * ) 'LOSUB.F: Error while reading inputfile!'
      STOP NO_INFILE
   11 CONTINUE
      READ (15, '(A)', END = 15)
      READ (15, *, END = 15) NE, LCALC
      READ (15, '(A)', END = 15)
      READ (15, *, END = 15) RA, TTA, PR, ETA, pL, Rconc
      READ (15, '(A)', END = 15)
      READ (15, *, END = 15) NT, M0
      READ (15, '(A)', END = 15)
      READ (15, *, END = 15) DRA, ABSE, RELE, NSMAX
      READ (15, '(A)', END = 15)
      READ (15, *, END = 15) TTSTEP, TTF
      CLOSE (15)
      GOTO 16
   15 WRITE ( * , * ) 'Error in inputfile ', inputfile
      STOP ERR_IN_INFILE
   16 CONTINUE
      endsubroutine

!**********************************************************************
      SUBROUTINE writeOutputHeader (outputfile)
      IMPLICIT double precision (A - H, O - Y)
      IMPLICIT complex (8) (Z)
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      COMMON / PAR2 / ETA
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      COMMON / NUM / TTA, TTF, TTSTEP, DRA, ABSE, RELE, NSMAX, LCALC
      CHARACTER (len = * ) outputfile
      OPEN (16, FILE = outputfile, STATUS = 'UNKNOWN')
      IF (LCALC.GT.0.AND.LCALC.LT.4.or.LCALC.eq.5.or.LCALC.eq.6) THEN
      WRITE (16,  * ) '### Output of Program lo.f Ver.2.1:        #'
      WRITE (16,  * ) '### Lin. Onset of Conv. via Galerkinmethod #'
      WRITE (16, '(A11,E12.5,A2)') '# P     ', PR, '#'
         WRITE (16, '(A11,E12.5,A2)') '# Lewis ', pL, '#'
      WRITE (16, '(A11,E12.5,A2)') '# TAU   ', TAU, '#'
      WRITE (16, '(A11,E12.5,A2)') '# R     ', RA, '#'
      WRITE (16, '(A11,E12.5,A2)') '# RC    ', Rconc, '#'
      WRITE (16, '(A11,E12.5,A2)') '# ETA   ', ETA, '#'
      WRITE (16, '(A11,G12.5,A2)') '# m     ', M0, '#'
      WRITE (16, '(A11,A12,A2)') '# cvar  ', 'TAU', '#'
      WRITE (16, '(A11,I12,A2)') '# NE    ', NE, '#'
      WRITE (16, '(A11,E12.5,A2)') '# TTA   ', TTA, '#'
      WRITE (16, '(A11,E12.5,A2)') '# TTF   ', TTF, '#'
         WRITE (16, '(A11,E12.5,A2)') '# TTSTEP', TTSTEP, '#'
      WRITE (16, '(A11,I12,A2)') '# NT    ', NT, '#'
      WRITE (16,  * ) '# see definition of LCALC for output. LCALC:', &
               LCALC, '   #'
      WRITE (16,  * ) '#                                      #'
      END IF

      CLOSE (16)
      endsubroutine

!**********************************************************************
      SUBROUTINE fixedParGrowthRate (outputfile)
      IMPLICIT double precision (A - H, O - Y)
      IMPLICIT complex (8) (Z)
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      COMMON / PAR2 / ETA
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      COMMON / NUM / TTA, TTF, TTSTEP, DRA, ABSE, RELE, NSMAX, LCALC
      CHARACTER (len = * ) outputfile
      INTEGER nloop
      IF (.true.) then
         Ra0 = RA
         niter = int ( (TTF - Ra0) / TTSTEP)
         CALL open_file_at_end (16, outputfile)
         DO nloop = 0, niter
         RA = Ra0 + nloop * TTSTEP
         GROR = FN (RA)
         WRITE (*,*) RA, GROR
         WRITE (16, '(3D16.8)') RA, GROR
         enddo
         CLOSE (16)
      END IF

      GROR = FN (RA)
      WRITE ( * ,  * ) 'R=', RA, ' TAU=', TAU, ' P=', PR, ' M0=', M0, ' &
     &eta=', ETA
      WRITE ( * , * ) 'Most unstable growth rate', GROR
      WRITE ( * , * ) 'If growth rate < 0 then above onset'
      endsubroutine

!**********************************************************************
      SUBROUTINE fixedParCriticalRa (outputfile)
      IMPLICIT double precision (A - H, O - Y)
      IMPLICIT complex (8) (Z)
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      COMMON / PAR2 / ETA
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      COMMON / NUM / TTA, TTF, TTSTEP, DRA, ABSE, RELE, NSMAX, LCALC
      CHARACTER (len = * ) outputfile
      CALL open_file_at_end (16, outputfile)
      CALL PEGASUS (RA, DRA, ABSE, RELE, NSMAX, 1, 1, LL, RAC)
      GROR = FN (RAC)
      OMEGA = C * M0
      WRITE ( * , * ) 'TAU=', TAU, ' P=', PR, ' M0=', M0, ' eta=', ETA
      WRITE ( * , * ) 'Lewis=', pL, ' Rconc=', Rconc
      WRITE ( * ,  * ) 'R_crit=', RAC, '  (growth rate =', GROR, ')'
      CLOSE (16)
      endsubroutine

!**********************************************************************
      SUBROUTINE fixedParCriticalRaAndM (outputfile, RAC, NTRYCOUNT)
      IMPLICIT double precision (A - H, O - Y)
      IMPLICIT complex (8) (Z)
      CHARACTER (len = * ) outputfile
      DOUBLE PRECISION RAC
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      COMMON / PAR2 / ETA
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      COMMON / NUM / TTA, TTF, TTSTEP, DRA, ABSE, RELE, NSMAX, LCALC
      INTEGER NTRYCOUNT

      CALL open_file_at_end (16, outputfile)
      CALL PEGASUS (RA, DRA, ABSE, RELE, NSMAX, 1, 1, LL, RAC)
      GROR = FN (RAC)
      OMEGA = C * M0
      IF (LL.EQ.0) THEN
         WRITE (16, '(1P,3E17.6,I4)') TAU, RAC, OMEGA, M0
         WRITE ( * , '(1P,3E17.6,I4)') TAU, RAC, OMEGA, M0
         NTRYCOUNT = 0
      ELSEIF (NTRYCOUNT.GE.3) THEN
         WRITE (16, * ) 'NO CRITICAL RAYLEIGH NUMBER FOUND.'
         STOP NO_RA_FOUND
      ELSE
         PRINT *, LL
      WRITE ( * ,  * ) 'NO CRIT. RAYLEIGH NUMBER FOUND. Trying again.'
         NTRYCOUNT = NTRYCOUNT + 1
      END IF
      CLOSE (16)
      endsubroutine

!**********************************************************************
      SUBROUTINE fixedParCriticalRaAndM0 (outputfile, RAC, NTRYCOUNT)
      IMPLICIT double precision (A - H, O - Y)
      IMPLICIT complex (8) (Z)
      CHARACTER (len = * ) outputfile
      DOUBLE PRECISION RAC
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      COMMON / PAR2 / ETA
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      COMMON / NUM / TTA, TTF, TTSTEP, DRA, ABSE, RELE, NSMAX, LCALC
      DOUBLE PRECISION RACI (3), OMI (3)
      INTEGER M0I (3), LLI (3), LMINI (3)
      INTEGER NTRYCOUNT
      COMPLEX (8) ZEVEC (NMAX)
      COMMON / MODE / ZEVEC

!     Start a search around the previous m0
      IF (M0.GT.1) THEN
         M0I (1) = M0 - 1
         LMINI (1) = LMIN - 1
      ELSE
         M0I (1) = M0
         LMINI (1) = LMIN
      END IF
      M0I (2) = M0
      M0I (3) = M0 + 1
      LMINI (2) = LMIN
      LMINI (3) = LMIN + 1
      RAC = RA
      DO II = 1, 3
         M0 = M0I (II)
         LMIN = LMINI (II)
         DRA = RA / 10.0D0
         CALL dimension (LMIN, LD, NT, M0, ND)
!            Print*, Ra, dRa, abse, rele, nsmax, ll, rac
!           LL is either 0 for success or 1 for failure
         CALL PEGASUS (RA, DRA, ABSE, RELE, NSMAX, 1, 0, LL, RAC)
         GROR = FN (RAC)
         OMI (II) = C * M0
         RACI (II) = RAC
         LLI (II) = LL
         WRITE ( * , '(1X,1P,4E17.6,I4,A3,2E17.6)') TAU, RAC, OMI (II) ,   &
         GROR, M0, ' | ', RAC * (1 / (1 - eta) ) **4 * (2 / TAU) , (OMI (  &
         II) * 2.0 / TAU)
      enddo

!       INDEX = 0
!       IF(LLI(1).EQ.0 .AND. LLI(2).EQ.0) THEN
!         IF( RACI(1).LT.RACI(2)) THEN
!             INDEX = 1
!         ELSE
!             INDEX = 2
!         END IF
!       ELSEIF(LLI(1).EQ.0) THEN
!             INDEX = 1
!       ELSEIF(LLI(2).EQ.0) THEN
!             INDEX = 2
!       END IF
!       IF(LLI(3).EQ.0 .AND. INDEX.GT.0) THEN
!         IF( RACI(3).LT.RACI(INDEX) ) THEN
!           INDEX = 3
!         END IF
!       ELSEIF(LLI(3).EQ.0) THEN
!           INDEX = 3
!       END IF
      IF (any (lli.eq.0) ) then
         where (lli.ne.0) RACI = 1.0d100
         index = minloc (RACI, 1)
      ELSE
         index = 0
      END IF

      CALL open_file_at_end (16, outputfile)
      IF (INDEX.GT.0) THEN
         RAC = RACI (INDEX)
         OMEGA = OMI (INDEX)
         M0 = M0I (INDEX)
         LMIN = LMINI (INDEX)
         WRITE (16, '(1P,3E17.6,I4)') TAU, RAC, OMEGA, M0
         WRITE ( * , '(">",1P,3E17.6,I4,A3,2E17.6)') TAU, RAC, OMEGA,   &
         M0, ' | ', RAC * (1 / (1 - eta) ) **4 * (2 / TAU) , (OMEGA *   &
         2.0 / TAU)
         WRITE (*,*)
         NTRYCOUNT = 0
      ELSEIF (NTRYCOUNT.GE.3) THEN
         WRITE (16, * ) 'NO CRITICAL RAYLEIGH NUMBER FOUND.'
         STOP NO_RA_FOUND
      ELSE
      WRITE ( * ,  * ) 'NO CRIT. RAYLEIGH NUMBER FOUND. Trying again.'
         NTRYCOUNT = NTRYCOUNT + 1
      END IF

      CLOSE (16)
      end subroutine

!**********************************************************************
      SUBROUTINE fixedParCriticalEigenVector (outputfile)
      IMPLICIT double precision (A - H, O - Y)
      IMPLICIT complex (8) (Z)
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      COMMON / PAR2 / ETA
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      COMMON / NUM / TTA, TTF, TTSTEP, DRA, ABSE, RELE, NSMAX, LCALC
      CHARACTER (len = * ) outputfile
      COMPLEX (8) ZEVEC (NMAX)
      COMMON / MODE / ZEVEC
      OPEN (16, FILE = outputfile, STATUS = 'UNKNOWN')
!--      searching for zero grothrate:
      CALL PEGASUS (RA, DRA, ABSE, RELE, NSMAX, 1, 1, LL, RAC)
      GROR = FN (RAC)
      OMEGA = C * M0

      IF (LL.EQ.0) THEN
!--       print eigenvector
!---------Fileformat for outputfile:
!---------LST=0 formatted Wicht
!---------LST=1 formatted Hirsching
!---------LST=3 unformatted
         LST = 0
         LCALC = 1
         WRITE (16, '(2I2,'' LINEAR ONSET '')') LST, LCALC
         NTH = 0
         KTV = 0
         KTH = 0
         LTV = 0
         LTH = 0
         GRR = 0.D0
         GRI = 0.D0
         WRITE (16, '(I2,7I3,2D16.8,'' M0,TRUNC,LD,GROTH,DRIFT'')') M0, &
         NT, NTH, KTV, KTH, LTV, LTH, LD, GRR, GRI
         NUDS = 1
         PM = 0.D0
         WRITE (16, '(I5,2D14.6,D9.2,D13.6,D9.2,'' I,TA,RA,PR,PM,E'')') &
         NUDS, TA, RAC, PR, PM, ETA
         C0 = OMEGA / M0
         OMM = 0.D0
         NUC = 0
         NUOM = 0
         MF = 0
         WRITE (16, 9100) C0, OMM, NUC, NUOM, MF

         LMAX = 2 * NT + M0 - 1
         I = 0
         DO LI = LMIN, LMAX, LD
!           L for poloidal (v) field:
         LPI = LI
         NIMAX = DINT (DBLE (2 * NT + 1 - LI + M0) / 2)
         DO NI = 1, NIMAX
         IF (LST.EQ.0) THEN
            WRITE (16, 9200) 'V', LPI, M0, NI, 0, DBLE (ZEVEC (I + 1) ) &
            , DIMAG (ZEVEC (I + 1) ) , 0.D0, 0.D0
         ELSEIF (LST.EQ.3) THEN
      WRITE (16, '(A,4I3,A,2F11.7,A)') ' V ', LPI, M0, NI, 0, ' ', DBLE &
     &(ZEVEC (I + 1) ) , DIMAG (ZEVEC (I + 1) ) , ' .0D+00 .0D+00 '
         END IF
         I = I + 4
         enddo
         enddo
!
         I = 0
         DO LI = LMIN, LMAX, LD
!           L for toroidal (w) field:
         IF (NE.EQ.2) THEN
            LTI = LI + 1
         ELSEIF (NE.EQ.1) THEN
            LTI = LI - 1
         ELSEIF (NE.EQ.0) THEN
            LTI = LI
         END IF
         NIMAX = DINT (DBLE (2 * NT + 1 - LI + M0) / 2)
         DO NI = 1, NIMAX
         IF (LST.EQ.0) THEN
            WRITE (16, 9200) 'W', LTI, M0, NI, 0, DBLE (ZEVEC (I + 3) ) &
            , DIMAG (ZEVEC (I + 3) ) , 0.D0, 0.D0
         ELSEIF (LST.EQ.3) THEN
      WRITE (16, '(A,4I3,A,2F11.7,A)') ' W ', LTI, M0, NI, 0, ' ', DBLE &
     &(ZEVEC (I + 3) ) , DIMAG (ZEVEC (I + 3) ) , ' .0D+00 .0D+00 '
         END IF
         I = I + 4
         enddo
         enddo
!
         I = 0
         DO LI = LMIN, LMAX, LD
         NIMAX = DINT (DBLE (2 * NT + 1 - LI + M0) / 2)
         DO NI = 1, NIMAX
         IF (LST.EQ.0) THEN
            WRITE (16, 9200) 'T', LI, M0, NI, 0, DBLE (ZEVEC (I + 2) ) ,&
            DIMAG (ZEVEC (I + 2) ) , 0.D0, 0.D0
         ELSEIF (LST.EQ.3) THEN
      WRITE (16, '(A,4I3,A,2F11.7,A)') ' ''T ''', LI, M0, NI, 0, ' ', DB&
     &LE (ZEVEC (I + 2) ) , DIMAG (ZEVEC (I + 2) ) , ' .0D+00 .0D+00 '
         END IF
         I = I + 4
         enddo
         enddo
!
         I = 0
         DO LI = LMIN, LMAX, LD
         NIMAX = DINT (DBLE (2 * NT + 1 - LI + M0) / 2)
         DO NI = 1, NIMAX
         IF (LST.EQ.0) THEN
            WRITE (16, 9200) 'G', LI, M0, NI, 0, DBLE (ZEVEC (I + 4) ) ,&
            DIMAG (ZEVEC (I + 4) ) , 0.D0, 0.D0
         ELSEIF (LST.EQ.3) THEN
      WRITE (16, '(A,4I3,A,2F11.7,A)') ' ''G ''', LI, M0, NI, 0, ' ', DB&
     &LE (ZEVEC (I + 4) ) , DIMAG (ZEVEC (I + 4) ) , ' .0D+00 .0D+00 '
         END IF
         I = I + 4
         enddo
         enddo
!
 9100 FORMAT   (2D17.10,3I4,'    C,OM, WHERE?,FLOQUET')
 9200 FORMAT   (1X,A1,4I3,4D16.8)
      ELSE
         WRITE (16, * ) 'NO CRITICAL RAYLEIGH NUMBER FOUND.'
         STOP NO_RA_FOUND
      END IF
      endsubroutine

!**********************************************************************
      SUBROUTINE losub (inputfile, outputfile)
      IMPLICIT double precision (A - H, O - Y)
      IMPLICIT complex (8) (Z)
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      COMMON / PAR2 / ETA
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      COMMON / NUM / TTA, TTF, TTSTEP, DRA, ABSE, RELE, NSMAX, LCALC
      CHARACTER (len = 40) infile, outfile
      COMMON / FILES / infile, outfile
      CHARACTER (len = * ) inputfile, outputfile
      INTEGER nloop
      COMPLEX (8) ZEVEC (NMAX)
      COMMON / MODE / ZEVEC

      infile = inputfile
      outfile = outputfile

!-----Default values:
      CALL setDefaults ()
      RAOLD = 4.D3
      RconcOLD = 4.D3
      TAU1 = 100.0D0
      TAU2 = 100.0D0

!-----INPUT:
      CALL readInputFile (inputfile)

!-----LOSUB.F doesn't work for M=0 !!!!!!
      IF (M0.LT.1) THEN
      WRITE ( * ,  * ) 'The code does not work for M0<1.', M0, ' --> 1'
         M0 = 1
      END IF

      TAU = TTA
      TAU1 = TTA
      TAU2 = TTA
      TA = TAU * TAU
      RAOLD = RA
      RconcOLD = Rconc

!-----OUTPUT:
      CALL writeOutputHeader (outputfile)
!
      RI = ETA / (1 - ETA)
      RO = 1.D0 + RI
!
      IF (NE.EQ.0) THEN
!-- UNDEFINED SYMMETRIE:
         LMIN = M0
         LD = 1
      ELSEIF (NE.EQ.1) THEN
!-- EQUATORIAL ANTISYMMETRIE (L+M ODD):
         LMIN = M0 + 1
         LD = 2
      ELSEIF (NE.EQ.2) THEN
!-- EQUATORIAL SYMMETRIE (L+M EVEN);
         LMIN = M0
         LD = 2
      END IF
!
      CALL DIMENSION (LMIN, LD, NT, M0, ND)
      WRITE ( * , * ) 'DIMENSION OF MATRIX:', ND
!
!---LCALC=-1 : most basic case: find the most unstable growth rate at al
      IF (LCALC.EQ. - 1) THEN
         CALL fixedParGrowthRate (outputfile)
         STOP ' Growth rate at fixed other parametres'
!---LCALC=0 : Critical Ra, for constant other parameters
      ELSEIF (LCALC.EQ.0) THEN
         CALL fixedParCriticalRa (outputfile)
         STOP ' Ra_crit at fixed other parametres'
!-----eigenvalues determined for this value of RA: ------------------
      ELSEIF (LCALC.EQ.1) THEN
         GROR = FN (RA)
         RETURN
!     LCALC-s which requiree increment of TAU:
      ELSEIF ( (LCALC.EQ.2) .or. (LCALC.EQ.3) ) then
         NTRYCOUNT = 0
         nloop = 0
         DO
         TA = TAU * TAU
!-----------searching for zero grothrate by varying RA: ----------------
         IF (LCALC.EQ.2) THEN
            CALL fixedParCriticalRaAndM (outputfile, RAC, NTRYCOUNT)
!-----------searching for zero grothrate by variing RA and M0: ---------
         ELSEIF (LCALC.EQ.3) THEN
            CALL fixedParCriticalRaAndM0 (outputfile, RAC, NTRYCOUNT)
         END IF
!--         increment TAU:
         TAU0 = TAU1
         TAU1 = TAU
         IF (DABS (TTSTEP) .LT.DABS (TAU * 0.1D0) ) THEN
            TAU = TAU1 + TTSTEP
         ELSE
            TAU = TAU1 + TTSTEP / DABS (TTSTEP) * TAU1 * 0.1D0
         END IF
!--         interpolate new startingvalue for RA:
         IF (TAU1.NE.TAU0) THEN
            RA = RAC + (RAC - RAOLD) / (TAU1 - TAU0) * (TAU - TAU1)
         ELSE
            RA = 2.0D0 * RAC - RAOLD
         END IF
         RAOLD = RAC
!--      endvalue of TAU reached?
         IF ( ( (TAU.GT.TTF) .AND. (TTF.GT.TTA) ) .OR. ( (TAU.LT.TTF)   &
         .AND. (TTF.LT.TTA) ) ) THEN
!              WRITE(*,*) 'LOSUB.F: finished at ',fdate()
            STOP FINISHED
         END IF

!--------count loops:
         nloop = nloop + 1
!--      End of tau Loop
         enddo
!-----calculate the critical eigenvector. Print for plotting.
      ELSEIF (LCALC.EQ.4) THEN
         CALL fixedParCriticalEigenVector (outputfile)
!-----vary m and calculate critical R at fixed P, tau, eta.
      ELSEIF (LCALC.EQ.5) THEN
         M0A = M0
         DO M0 = M0A, TTF, INT (TTSTEP)
!--      UNDEFINED SYMMETRIE:
            IF (NE.EQ.0) THEN
               LMIN = M0
               LD = 1
!--      EQUATORIAL ANTISYMMETRIE (L+M ODD):
            ELSEIF (NE.EQ.1) THEN
               LMIN = M0 + 1
               LD = 2
!--      EQUATORIAL SYMMETRIE (L+M EVEN);
            ELSEIF (NE.EQ.2) THEN
               LMIN = M0
               LD = 2
            END IF
            CALL dimension (LMIN, LD, NT, M0, ND)
            DRA = RA / 10.0D0
            CALL PEGASUS (RA, DRA, ABSE, RELE, NSMAX, 1, 0, LL, RAC)
            WRITE (*,*) M0, RAC
            CALL open_file_at_end (16, outputfile)
            WRITE (16,*) M0, RAC
            CLOSE (16)
         enddo
      ELSEIF (LCALC.EQ.6) THEN
!-----vary Le and calculate critical R at fixed P, tau, eta, M
         pL0 = pL
         niter = (TTF - pL0) / TTSTEP
         DO nloop = 0, niter
            pL = pL0 + nloop * TTSTEP
            DRA = RA / 10.0D0
            CALL PEGASUS (RA, DRA, ABSE, RELE, NSMAX, 1, 0, LL, RAC)
            GROR = FN (RAC)
            WRITE(*,*) pL, RAC, GROR
            CALL open_file_at_end (16, outputfile)
            WRITE(16,'(3D16.8)') pL, RAC, GROR
            CLOSE(16)
         enddo
      END IF
      END SUBROUTINE losub

!***********************************************************************
!-- FUNCTION CALLED BY PEGASUS, FN=GROTHRATE.
      FUNCTION FN (RA)
      IMPLICIT double precision (A - H, O - Y)
      IMPLICIT complex (8) (Z)
      DIMENSION ZA (NMAX, NMAX), ZB (NMAX, NMAX)
      DIMENSION ZACOPY (NMAX, NMAX), ZBCOPY (NMAX, NMAX)
      DIMENSION ZEW (NMAX), ZEWA (NMAX), ZEWB (NMAX)
      DIMENSION ZEVALL (NMAX, NMAX), ZEVEC (NMAX)
      COMMON / PAR / TAU, RAL, PR, RI, C, pL, Rconc
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      COMMON / MODE / ZEVEC
      COMMON / NUM / TTA, TTF, TTSTEP, DRA, ABSE, RELE, NSMAX, LCALC

      RAL = RA

!-- MAT SETS THE complex(8) MATRICES ZA AND ZB SETTING OF MATRIX:
      CALL MAT (ZA, ZB, NMAX)

      ZEWA   = DCMPLX (0D0, 0D0)
      ZEWB   = DCMPLX (0D0, 0D0)
      ZACOPY = DCMPLX (0D0, 0D0)
      ZBCOPY = DCMPLX (0D0, 0D0)
!
!-- DGVLCG IS AN IMSL ROUTINE THAT CALCULATES THE EIGENVALUES ZEWA/ZEWB:
!-- Problem: generalized complex(8) eigensystem A*x = lam*B*x
!--          Input:  A: ZA, B: ZB
!--          Output: complex(8) ZEWA(ND), ZEWB(ND)
!      CALL DGVLCG(ND,ZA,NMAX,ZB,NMAX,ZEWA,ZEWB)
      IF (LCALC.EQ.4) THEN
         CALL DG2CCG(ND, ZA, NMAX, ZB, NMAX, ZEWA, ZEWB, ZEVALL, NMAX, &
         ZACOPY, ZBCOPY)
      ELSE
         CALL DG2LCG(ND, ZA, NMAX, ZB, NMAX, ZEWA, ZEWB, ZACOPY, ZBCOPY)
      END IF
!
      DO I = 1, ND
         ZEW (I) = ZEWA (I) / ZEWB (I)
      enddo
      IMIN = 1
!---- search for lowest imaginary part:
      EWMIN = DIMAG (ZEW (IMIN) )
      DO I = 2, ND
         IF (DIMAG (ZEW (I) ) .LT.EWMIN) THEN
            IMIN = I
            EWMIN = DIMAG (ZEW (IMIN) )
         END IF
      END DO

!      IF( EWMIN.LT.0.D0 ) WRITE(*,*) 'CONVECTION REGIME !',RAL
      C = - DBLE (ZEW (IMIN) ) / M0
!
      IF (LCALC.EQ.1) THEN
!---- sort eigenvalues:
         WRITE ( * ,  * ) '     Frequ.(exp(+iwt))   -Grothrate  '
         DO I = 1, ND
            DO J = I, ND
               IF (DIMAG (ZEW (J) ) .LT.DIMAG (ZEW (I) ) ) THEN
                  ZSAVE = ZEW (J)
                  ZEW (J) = ZEW (I)
                  ZEW (I) = ZSAVE
                  DO K = 1, ND
                     ZEVEC (K) = ZEVALL (K, J)
                     ZEVALL (K, J) = ZEVALL (K, I)
                     ZEVALL (K, I) = ZEVEC (K)
                  END DO
               END IF
            END DO
         END DO
         WRITE (*,'(I4,2D16.6)') I, ZEW (I)
      END IF

!------save crit. eigenvector:
      IF (LCALC.EQ.4) THEN
         DO J = 1, ND
            ZEVEC (J) = ZEVALL (J, IMIN)
         END DO
      END IF

      FN = EWMIN
!      if ( fn.lt.0.0) then
!         print*,  ' Convection onset found :) :  Ra= ...  ',RAL
!      else
!         print*,  ' Stability :( :  Ra= ...  ',RAL
!      endif
      END FUNCTION FN
!
!***********************************************************************
!-- SETS THE complex(8) MATRICES ZA AND ZB.
      SUBROUTINE MAT (ZA, ZB, NDIM)
      IMPLICIT double precision (A - H, O - Y)
      IMPLICIT complex (8) (Z)
!
      DIMENSION ZA (NDIM, * ), ZB (NDIM, * )
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
!
      DO J = 1, NDIM
      DO I = 1, NDIM
      ZA (I, J) = DCMPLX (0D0, 0D0)
      ZB (I, J) = DCMPLX (0D0, 0D0)
      enddo
      enddo
!
      I = 0
      LMAX = 2 * NT + M0 - 1
!
      DO LI = LMIN, LMAX, LD
      LPI = LI
!
!        Determine L for toroidal (w) field:
      IF (NE.EQ.2) THEN
         LTI = LI + 1
      ELSEIF (NE.EQ.1) THEN
         LTI = LI - 1
      ELSEIF (NE.EQ.0) THEN
         LTI = LI
      END IF
!
      NIMAX = DINT (DBLE (2 * NT + 1 - LI + M0) / 2)
!
      DO NI = 1, NIMAX
!
      J = 0
      DO LJ = LMIN, LMAX, LD
      LPJ = LJ
!
      IF (NE.EQ.2) THEN
         LTJ = LJ + 1
      ELSEIF (NE.EQ.1) THEN
         LTJ = LJ - 1
!               ELSEIF( NE.EQ.0 ) THEN
      ELSE
         LTJ = LJ
      END IF
!
      NJMAX = DINT (DBLE (2 * NT + 1 - LJ + M0) / 2)
!
!  ******************** I: Equation (Line) ******************
!  ******************** J: Variable (Column) ****************
!  ******************** I+1: v (poloidal)  ******************
!  ******************** I+2: theta         ******************
!  ******************** I+3: w (toroidal)  ******************
!-new****************** I+4: gamma (concentration) **********
      DO NJ = 1, NJMAX
!
      IF (J + 3.GT.NDIM.OR.I + 3.GT.NDIM) THEN
         WRITE ( * , * ) 'MAT(): NDIM too small.'
         STOP
      END IF
!
      IF (LI.EQ.LJ) THEN
         ZB (I + 1, J + 1) = DCMPLX (0.D0, - DIII2 (NI, NJ, LPI, 1) )
         ZA (I + 1, J + 1) = DCMPLX (DIII1 (NI, NJ, LPI), DIII3 (NI, NJ,&
         LPI, 1) )
         ZA (I + 1, J + 2) = DCMPLX (DIII5 (NI, NJ, LPI), 0.D0)
!--- concentration driving
         ZA (I + 1, J + 4) = DCMPLX (DIII5conc (NI, NJ, LPI), 0.D0)
!
         ZB (I + 2, J + 2) = DCMPLX (0.D0, - DI1 (NI, NJ, 1) )
         ZA (I + 2, J + 1) = DCMPLX (DI3 (NI, NJ, LPI), 0.D0)
         ZA (I + 2, J + 2) = DCMPLX (DI2 (NI, NJ, LPI), 0.D0)
         ZB (I + 3, J + 3) = DCMPLX (0.D0, - DII2 (NI, NJ, LTI, 1) )
         ZA (I + 3, J + 3) = DCMPLX (DII1 (NI, NJ, LTI), DII3 (NI, NJ,  &
         1) )
!--- concentration equation
         ZB (I + 4, J + 4) = DCMPLX (0.D0, - DI1 (NI, NJ, 1) )
         ZA (I + 4, J + 1) = DCMPLX (DI3 (NI, NJ, LPI), 0.D0)
         ZA (I + 4, J + 4) = DCMPLX (1.D0 / pL * DI2 (NI, NJ, LPI),     &
         0.D0)
      END IF
      IF (LPI.EQ.LTJ + 1) THEN
         ZA (I + 1, J + 3) = DCMPLX (DIII4A (NI, NJ, LPI, 1), 0.D0)
      ELSEIF (LPI.EQ.LTJ - 1) THEN
         ZA (I + 1, J + 3) = DCMPLX (DIII4B (NI, NJ, LPI, 1), 0.D0)
      END IF
      IF (LTI.EQ.LPJ + 1) THEN
         ZA (I + 3, J + 1) = DCMPLX (DII4A (NI, NJ, LTI, 1), 0.D0)
      ELSEIF (LTI.EQ.LPJ - 1) THEN
         ZA (I + 3, J + 1) = DCMPLX (DII4B (NI, NJ, LTI, 1), 0.D0)
      END IF
      J = J + 4
      END DO
!
      END DO
!
      I = I + 4
      END DO
!
      END DO
      END SUBROUTINE MAT
!
!***********************************************************************
!-- GALERKIN TERMS:
!***********************************************************************
      FUNCTION DI1 (N1, N2, NU1)
!----- HEAT EQUATION, TIME DERIVATIVE
      IMPLICIT double precision (A - H, O - Z)
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      DI1 = PR * NU1 * R ('SS ', 2, N1, N2, 0)
      END FUNCTION DI1
!
      FUNCTION DI2 (N1, N2, L1)
!----- HEAT EQUATION , DISSIPATION
      IMPLICIT double precision (A - H, O - Z)
      PARAMETER (DPI = 3.141592653589793D0)
      DI2 = N2**2 * DPI**2 * R ('SS ', 2, N1, N2, 0) - 2 * N2 * DPI * R &
      ('SC ', 1, N1, N2, 0) + DL (L1) * R ('SS ', 0, N1, N2, 0)
      END FUNCTION DI2
!
!----- HEAT EQUATION , SOURCE
      FUNCTION DI3 (N1, N2, L1)
      IMPLICIT double precision (A - H, O - Z)
      DI3 = - DL (L1) * R ('SS ', 2, N1, N2, 0)
      END FUNCTION DI3
!
!----- TOROIDAL EQUATION , DISSIPATION
      FUNCTION DII1 (N1, N2, L1)
      IMPLICIT double precision (A - H, O - Z)
      PARAMETER (DPI = 3.141592653589793D0)
      DII1 = DL (L1) * ( (N2 - 1) **2 * DPI**2 * R ('CC ', 4, N1 - 1,   &
      N2 - 1, 0) + 4 * (N2 - 1) * DPI * R ('CS ', 3, N1 - 1, N2 - 1, 0) &
      + (DL (L1) - 2) * R ('CC ', 2, N1 - 1, N2 - 1, 0) )
      END FUNCTION DII1
!
!----- TOROIDAL EQUATION , TIME DERIVATIVE
      FUNCTION DII2 (N1, N2, L1, NU1)
      IMPLICIT double precision (A - H, O - Z)
      DII2 = NU1 * DL (L1) * R ('CC ', 4, N1 - 1, N2 - 1, 0)
      END FUNCTION DII2
!
!----- TOROIDAL EQUATION , CORRIOLIS
      FUNCTION DII3 (N1, N2, NU1)
      IMPLICIT double precision (A - H, O - Z)
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      DII3 = - TAU * NU1 * M0 * R ('CC ', 4, N1 - 1, N2 - 1, 0)
      END FUNCTION DII3
!
!----- TOROIADL EQUATION , Q-TERM 1 (L1=L3+1)
      FUNCTION DII4A (N1, N2, L1, NU1)
      IMPLICIT double precision (A - H, O - Z)
      PARAMETER (DPI = 3.141592653589793D0)
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      DII4A = TAU * DSQRT (DBLE (L1 - NU1 * M0) * (L1 + NU1 * M0)       &
      / (2 * L1 - 1) / (2 * L1 + 1) ) * ( (L1**2 - 1) * (L1 - 1)        &
      * R ('CS ', 2, N1 - 1, N2, 0) - (L1 + 1) * (L1 - 1) * N2 * DPI *  &
      R ('CC ', 3, N1 - 1, N2, 0) )
      END FUNCTION DII4A
!
!----- TOROIADL EQUATION , Q-TERM 1 (L1=L3-1)
      FUNCTION DII4B (N1, N2, L1, NU1)
      IMPLICIT double precision (A - H, O - Z)
      PARAMETER (DPI = 3.141592653589793D0)
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      DII4B = TAU * DSQRT (DBLE (L1 - NU1 * M0 + 1) * (L1 + NU1 * M0 +  &
      1) / (2 * L1 + 1) / (2 * L1 + 3) ) * ( (1 - (L1 + 1) **2) *       &
      (L1 + 2) * R ('CS ', 2, N1 - 1, N2, 0) - L1 * (L1 + 2) * N2 * DPI &
      * R ('CC ', 3, N1 - 1, N2, 0) )
      END FUNCTION DII4B
!
!----- POLOIDAL EQUOATION , DISSIPATION
      FUNCTION DIII1 (N1, N2, L1)
      IMPLICIT double precision (A - H, O - Z)
      PARAMETER (DPI = 3.141592653589793D0)
      DIII1 = DL (L1) * (N2**4 * DPI**4 * R ('SS ', 2, N1, N2, 0)       &
      - 4 * N2**3 * DPI**3 * R ('SC ', 1, N1, N2, 0) + 2 * DL (L1)      &
      * N2**2 * DPI**2 * R ('SS ', 0, N1, N2, 0) + (DL (L1) **2 - 2 *   &
      DL (L1) ) * R ('SS ', - 2, N1, N2, 0) )
      END FUNCTION DIII1
!
!----- POLOIDAL EQUATION , TIME DERIVATIVE
      FUNCTION DIII2 (N1, N2, L1, NU1)
      IMPLICIT double precision (A - H, O - Z)
      PARAMETER (DPI = 3.141592653589793D0)
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      DIII2 = - NU1 * DL (L1) * ( - N2**2 * DPI**2 * R ('SS ', 2, N1,   &
      N2, 0) + 2 * N2 * DPI * R ('SC ', 1, N1, N2, 0) - DL (L1) * R (   &
      'SS ', 0, N1, N2, 0) )
      END FUNCTION DIII2
!
!----- POLOIDAL EQUATION , CORRIOLIS
      FUNCTION DIII3 (N1, N2, L1, NU1)
      IMPLICIT double precision (A - H, O - Z)
      PARAMETER (DPI = 3.141592653589793D0)
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      DIII3 = TAU * NU1 * M0 * ( - N2**2 * DPI**2 * R ('SS ', 2, N1, N2,&
      0) + 2 * N2 * DPI * R ('SC ', 1, N1, N2, 0) - DL (L1) * R ('SS ', &
      0, N1, N2, 0) )
      END FUNCTION DIII3
!
      FUNCTION DIII4A (N1, N2, L1, NU1)
!----- POLOIDAL EUQUATION , Q-TERM 1 (L1=L3+1)
      IMPLICIT double precision (A - H, O - Z)
      PARAMETER (DPI = 3.141592653589793D0)
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      DIII4A = TAU * DSQRT (DBLE (L1 - M0 * NU1) * (L1 + M0 * NU1)      &
      / (2 * L1 - 1) / (2 * L1 + 1) ) * ( (L1 * (L1 - 1) - 2) * (L1 - 1)&
      * R ('SC ', 2, N1, N2 - 1, 0) + (L1 + 1) * (L1 - 1) * (N2 - 1)    &
      * DPI * R ('SS ', 3, N1, N2 - 1, 0) )
      END FUNCTION DIII4A
!
      FUNCTION DIII4B (N1, N2, L1, NU1)
!----- POLOIDAL EQUATION , Q-TERM 2 (L1=L3-1)
      IMPLICIT double precision (A - H, O - Z)
      PARAMETER (DPI = 3.141592653589793D0)
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      DIII4B = TAU * DSQRT (DBLE (L1 - M0 * NU1 + 1) * (L1 + M0 * NU1 + &
      1) / (2 * L1 + 1) / (2 * L1 + 3) ) * ( (L1 + 2) * (2 - (L1 + 1)   &
      * (L1 + 2) ) * R ('SC ', 2, N1, N2 - 1, 0) + L1 * (L1 + 2)        &
      * (N2 - 1) * DPI * R ('SS ', 3, N1, N2 - 1, 0) )
      END FUNCTION DIII4B
!
      FUNCTION DIII5 (N1, N2, L1)
!----- POLOIDAL EQUATION ,
      IMPLICIT double precision (A - H, O - Z)
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      DIII5 = - RA * DL (L1) * R ('SS ', 2, N1, N2, 0)
      END FUNCTION DIII5

      FUNCTION DIII5conc (N1, N2, L1)
!----- POLOIDAL EQUATION ,
      IMPLICIT double precision (A - H, O - Z)
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      DIII5conc = - Rconc * DL (L1) * R ('SS ', 2, N1, N2, 0)
      END FUNCTION DIII5conc
!
!***********************************************************************
!-- SUBROUTINE S:
!***********************************************************************
      FUNCTION DL (L)
      IMPLICIT none
      DOUBLE PRECISION DL
      INTEGER L
      DL = DBLE (L * (L + 1) )
      END FUNCTION DL
!
!-- DUMMY FUNCTION S:
!
      FUNCTION F2 (X)
      F2 = X
      END FUNCTION F2
!
      FUNCTION F3 (X)
      F3 = X
      END FUNCTION F3
!
      FUNCTION F4 (X)
      F4 = X
      END FUNCTION F4

!***********************************************************************
      SUBROUTINE DIMENSION (LMIN, LD, NT, M0, ND)
!***********************************************************************
      INTEGER NT, M0, ND, LMIN, LD
      INTEGER L
!-- DETERMINATION OF DIMENSION:
!-- for each value of L the number of possible N-values is added
!         print*, "Triangular truncation (2.12)"
!         print*, LMIN, "...", 2*NT+M0-1,LD
      ND = 0
      DO L = LMIN, 2 * NT + M0 - 1, LD
!         print*, L, 1, "...", INT( DBLE(2*NT+1-L+M0)/2 )
!ccccccccc18    ND=ND+3*DINT( DBLE(2*NT+1-L+M0)/2 )
      ND = ND+4 * DINT (DBLE (2 * NT + 1 - L + M0) / 2)
      END DO
!
      IF (ND.GT.NMAX) THEN
         WRITE ( * , * ) 'DIMENSION OF MATRIX TOO SMALL:', ND, '>',     &
         NMAX
         STOP DIM_TO_SMALL
      END IF

      END SUBROUTINE DIMENSION

!***********************************************************************
      SUBROUTINE open_file_at_end (NHANDLE, filename)
!***********************************************************************
!     opens file <filename> and puts the filepointer at EOF
!***********************************************************************
      INTEGER NHANDLE
      CHARACTER (len = * ) filename

      OPEN (NHANDLE, FILE = filename, STATUS = 'OLD', POSITION =        &
      'APPEND ', ERR = 990)
      GOTO 999
  990 WRITE ( * , * ) 'Error reading ', filename
      STOP ERR_WRT_OUTFILE
  999 CONTINUE
      END SUBROUTINE open_file_at_end

!***********************************************************************
      SUBROUTINE abort1
!     this function should replace the internal FORTRAN routine
!     'abort' which is e.g. CALLed from the IMSL routines
!***********************************************************************
      IMPLICIT double precision (A - H, O - Y)
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      COMMON / NUM / TTA, TTF, TTSTEP, DRA, ABSE, RELE, NSMAX, LCALC

      WRITE ( * , * ) 'LOSUB.F: abort(fortran) CALLed.'
!-----increment TAU:
      IF (DABS (TTSTEP) .LT.DABS (TAU * 0.1D0) ) THEN
         TAU = TAU + TTSTEP
      ELSE
         TAU = TAU + TAU * 0.1D0
      END IF

!--   endvalue of TAU reached?
      IF ( ( (TAU.GT.TTF) .AND. (TTF.GT.TTA) ) .OR. ( (TAU.LT.TTF)      &
      .AND. (TTF.LT.TTA) ) ) THEN
!      WRITE(*,*) 'LOSUB.F: finished at ',fdate()
         STOP FINISHED
      END IF

      CALL ende (SIGRESTART)
      RETURN
      END SUBROUTINE abort1

!***********************************************************************
      SUBROUTINE ende (II)
!     ende() is CALLed explicit in the program or implicit when the
!     process receives a signal.
!     II indicates the behaviour of ende() or the received signal number
!***********************************************************************
      IMPLICIT double precision (A - H, O - Y)
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc
      COMMON / PAR2 / ETA
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE
      COMMON / NUM / TTA, TTF, TTSTEP, DRA, ABSE, RELE, NSMAX, LCALC
      CHARACTER (40) infile, outfile
      COMMON / FILES / infile, outfile
      INTEGER II, NHANDLE
      LOGICAL LOP
!      CHARACTER*24 fdate
!      EXTERNAL     fdate

      WRITE ( * ,  * ) 'ende(', II, '):  TAU = ', TAU

!      IF(II.EQ.SIGUSR1) THEN
!        write(*,'(A,A)') ' ende(): received signal SIGUSR1 at ',fdate()
!      ELSEIF(II.EQ.SIGABRT) THEN
!        write(*,'(A,A)') ' ende(): received signal SIGABRT at ',fdate()
!      END IF

      INQUIRE (FILE = trim (outfile), OPENED = LOP, NUMBER = NHANDLE)
      IF (LOP) THEN
         CLOSE (NHANDLE)
      END IF

      IF (II.EQ.SIGRESTART.OR.II.EQ.SIGUSR1) THEN
      WRITE ( * ,  * ) 'LOSUB.F: Terminating this run and starting next.&
     &'
         OPEN (15, FILE = infile, STATUS = 'UNKNOWN')
         WRITE (15, * ) ' NE (0/1/2) | LCALC (1/2/3/4) |'
         WRITE (15, '(A,2I12)') ' ', NE, LCALC
      WRITE (15,  * ) '|  RAYLEIGH  |  TAU     |  PRANTEL  |  ETA  |    &
     &      Lewis |   Rconc   |'
         WRITE (15, '(1P,E17.6,5(A,E17.6))') RA, ' ', TAU, ' ', PR, ' ',&
         ETA, ' ', pL, ' ', Rconc
      WRITE (15,  * ) '|   NTRUNC (>=1) | MODE |'
         WRITE (15, '(A,2I12)') ' ', NT, M0
      WRITE (15,  * ) '|   DRA   | ABSERR  |  RELERR  | NMAX |'
         WRITE (15, '(1PG12.6,A,1PG11.5,A,1PG11.5,A,I4)') DRA, ' ',     &
         ABSE, ' ', RELE, ' ', NSMAX
      WRITE (15,  * ) '|   TAU_STEP | TAU_END '
         WRITE (15, '(1P,2G11.4)') TTSTEP, TTF
         CLOSE (15)
         STOP 10
!--------return value START_NEXT_RUN to the shell script:
!         CALL exit(START_NEXT_RUN)
      ELSE
!         CALL exit(1)
      END IF
      END SUBROUTINE ende
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
