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
#include "losub-inc.h"
***********************************************************************
      subroutine setDefaults()
      IMPLICIT double precision(A-H,O-Y)
      IMPLICIT complex(8)(Z)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      COMMON/PAR2/ETA
      COMMON/DIM/NT,M0,ND,LMIN,LD,NE
      COMMON/NUM/TTA,TTF,TTSTEP,DRA,ABSE,RELE,NSMAX,LCALC
      NE = 2
      RA = 4.D3
      Rconc = 4.D3
      TAU = 100.0D0
      PR = 1.D-1
      pL = 1.D00
      ETA = 0.4D0
      NT = 3
      M0 = 6
      DRA = RA/10.0D0
      ABSE = 0.D0
      RELE = 1.D-6
      NSMAX = 100
      LCALC = 2
      end subroutine
  
***********************************************************************
      subroutine readInputFile(inputfile)
      IMPLICIT double precision(A-H,O-Y)
      IMPLICIT complex(8)(Z)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      COMMON/PAR2/ETA
      COMMON/DIM/NT,M0,ND,LMIN,LD,NE
      COMMON/NUM/TTA,TTF,TTSTEP,DRA,ABSE,RELE,NSMAX,LCALC
      CHARACTER(len=*) inputfile
      OPEN(15,FILE=inputfile,STATUS='OLD',ERR=10)
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
      end subroutine

***********************************************************************
      subroutine writeOutputHeader(outputfile)
      IMPLICIT double precision(A-H,O-Y)
      IMPLICIT complex(8)(Z)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      COMMON/PAR2/ETA
      COMMON/DIM/NT,M0,ND,LMIN,LD,NE
      COMMON/NUM/TTA,TTF,TTSTEP,DRA,ABSE,RELE,NSMAX,LCALC
      CHARACTER(len=*) outputfile
      OPEN(16,FILE=outputfile,STATUS='UNKNOWN')
      IF(LCALC.GT.0 .AND. LCALC.LT.4 .or. LCALC.eq.5.or.LCALC.eq.6) THEN
       WRITE(16,*)  '### Output of Program lo.f Ver.2.1:        ###'
       WRITE(16,*)  '### Lin. Onset of Conv. via Galerkinmethod ###'
       WRITE(16,'(A11,E12.5,A2)') '# P     ', PR,     '#'
       WRITE(16,'(A11,E12.5,A2)') '# Lewis ', pL,     '#'
       WRITE(16,'(A11,E12.5,A2)') '# TAU   ', TAU,    '#'
       WRITE(16,'(A11,E12.5,A2)') '# R     ', RA,     '#'
       WRITE(16,'(A11,E12.5,A2)') '# RC    ', Rconc,  '#'
       WRITE(16,'(A11,E12.5,A2)') '# ETA   ', ETA,    '#'
       WRITE(16,'(A11,G12.5,A2)') '# m     ', M0,     '#'
       WRITE(16,'(A11,A12,A2)')   '# cvar  ','TAU',   '#'
       WRITE(16,'(A11,I12,A2)')   '# NE    ', NE,     '#'
       WRITE(16,'(A11,E12.5,A2)') '# TTA   ', TTA,    '#'
       WRITE(16,'(A11,E12.5,A2)') '# TTF   ', TTF,    '#'
       WRITE(16,'(A11,E12.5,A2)') '# TTSTEP',TTSTEP,  '#'
       WRITE(16,'(A11,I12,A2)')   '# NT    ', NT,     '#'
       WRITE(16,*)  '# see definition of LCALC for output. LCALC:',
     &            LCALC,'   #'
       WRITE(16,*)  '#                                      #'
      ENDIF
      
      CLOSE(16)
      end subroutine

***********************************************************************
      subroutine fixedParGrowthRate(outputfile)
      IMPLICIT double precision(A-H,O-Y)
      IMPLICIT complex(8)(Z)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      COMMON/PAR2/ETA
      COMMON/DIM/NT,M0,ND,LMIN,LD,NE
      COMMON/NUM/TTA,TTF,TTSTEP,DRA,ABSE,RELE,NSMAX,LCALC
      CHARACTER(len=*) outputfile
      INTEGER aux, factor

      Ra0=RA

      if (Ra0.gt.TTF) then
         factor = -1.0d0
      else
         factor = 1.0d0
      endif
      CALL open_file_at_end(16,outputfile)
      do while(abs(Ra-ttf) > abs(0.11*ttf))
         aux = int(log10(abs(Ra)))
         if(aux.ge.2) then
            TTSTEP = factor*10.0d0**(aux-1)
         else
            TTSTEP = factor*10.0d0
         endif
         RA=Ra+TTSTEP
         GROR=FN(RA)
         WRITE(*,*) RA,GROR
         Write(16, *) Ra, GroR
      enddo
      close(16)

      WRITE(*,*) 'R=',RA,' TAU=',TAU,' P=',PR,' M0=',M0,' eta=',ETA
      WRITE(*,*) 'Most unstable growth rate', GROR
      WRITE(*,*) 'If growth rate < 0 then above onset'
      end subroutine

***********************************************************************
      subroutine fixedParCriticalRa(outputfile)
      IMPLICIT double precision(A-H,O-Y)
      IMPLICIT complex(8)(Z)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      COMMON/PAR2/ETA
      COMMON/DIM/NT,M0,ND,LMIN,LD,NE
      COMMON/NUM/TTA,TTF,TTSTEP,DRA,ABSE,RELE,NSMAX,LCALC
      CHARACTER(len=*) outputfile
      CALL open_file_at_end(16,outputfile)
         CALL pegasus(RA,DRA,ABSE,RELE,NSMAX,1,QUIET,LL,RAC)
         GROR=FN(RAC)
         OMEGA=C*M0
         WRITE(*,*) 'TAU=',TAU,' P=',PR,' M0=',M0,' eta=',ETA
         WRITE(*,*) 'Lewis=',pL,' Rconc=',Rconc
      WRITE(*,*) 'R_crit=',RAC, '  (growth rate =',GROR,')'
      close(16)
      end subroutine

***********************************************************************
      subroutine fixedParCriticalRaAndM(outputfile, RAC, NTRYCOUNT)
      IMPLICIT double precision(A-H,O-Y)
      IMPLICIT complex(8)(Z)
      CHARACTER(len=*) outputfile
      double precision RAC
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      COMMON/PAR2/ETA
      COMMON/DIM/NT,M0,ND,LMIN,LD,NE
      COMMON/NUM/TTA,TTF,TTSTEP,DRA,ABSE,RELE,NSMAX,LCALC
      integer NTRYCOUNT

      CALL open_file_at_end(16,outputfile)
      CALL pegasus(RA,DRA,ABSE,RELE,NSMAX,1,QUIET,LL,RAC)
      GROR=FN(RAC)
      OMEGA=C*M0
      IF(LL.EQ.0) THEN
         WRITE(16,'(1P,3E17.6,I4)') TAU,RAC,OMEGA,M0
         WRITE( *,'(1P,3E17.6,I4)') TAU,RAC,OMEGA,M0
         NTRYCOUNT = 0
      ELSE IF(NTRYCOUNT.GE.3) THEN
         WRITE(16,*) 'NO CRITICAL RAYLEIGH NUMBER FOUND.'
         STOP NO_RA_FOUND
      ELSE
         write(*,*) 'NO CRIT. RAYLEIGH NUMBER FOUND. Trying again.'
         NTRYCOUNT = NTRYCOUNT + 1
      ENDIF
      close(16)
      end subroutine

***********************************************************************
      subroutine fixedParCriticalRaAndM0(outputfile, RAC, NTRYCOUNT)
      IMPLICIT double precision(A-H,O-Y)
      IMPLICIT complex(8)(Z)
      CHARACTER(len=*) outputfile
      double precision RAC
      parameter(nm0=3)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      COMMON/DIM/NT,M0,ND,LMIN,LD,NE
      COMMON/NUM/TTA,TTF,TTSTEP,DRA,ABSE,RELE,NSMAX,LCALC
      double precision RACI(nm0),OMI(nm0),GRORi(nm0)
      INTEGER M0I(nm0),LLI(nm0),LMINI(nm0)
      INTEGER NTRYCOUNT
      complex(8) ZEVEC(NMAX)
      COMMON/MODE/ZEVEC
      integer first_m0, first_lmin, i, ii

      do i=0, nm0/2
         first_m0   = m0 - i
         first_lmin = lmin-i
         if (first_m0==1) exit
      enddo
c     Start a search around the previous m0
      do i=1, nm0
        M0I(i)   = first_m0+i-1
        LMINI(i) = first_lmin+i-1
      enddo
      DO II=1, nm0
         M0   = M0I(II)
         LMIN = LMINI(II)
C        Search for the new Rac in steps of 1/10 of the previous Rac.         
         DRA = RA/10.0D0
         CALL dimension(LMIN,LD,NT,M0,ND)
c        LL is either 0 for success or 1 for failure
         CALL pegasus(RA,DRA,ABSE,RELE,NSMAX,1,VERBOSE,LL,RAC)
C        Compute the grouth rate for each case so we can output it for
C        debug purposes.
         GRORi(ii)    = FN(RAC)
c        C is set inside FN 
         OMI(II) = C*M0
         RACI(II)= RAC
         LLI(II) = LL
         write( *,'(1X,1P,4E17.6,I4,A3,2E17.6)')
     &      TAU, RAC, OMI(II), GRORi(ii), M0,' | ',
     &          RAC*(1/(1-eta))**4*(2/TAU),(OMI(II)*2.0/TAU)
      enddo

C      Let's not allow growth rates that are too high
      where(grori.gt.1.0d0) lli = 1
      if(any(lli.eq.0)) then
         where(lli.ne.0) RACI = 1.0d100
         index = minloc(RACI,1)
      else
         index = 0
      endif

      CALL open_file_at_end(16,outputfile)
      IF( INDEX.GT.0 ) THEN
         RAC   = RACI(INDEX)
         OMEGA = OMI(INDEX)
         M0    = M0I(INDEX)
         LMIN  = LMINI(INDEX)
         write( 16,'(1P,3E17.6,I4)') TAU, RAC, OMEGA, M0
         write( *,'(">",1P,3E17.6,I4,A3,2E17.6)')
     &     TAU,RAC,                        OMEGA,    M0,' | ',
     &         RAC*(1/(1-eta))**4*(2/TAU),(OMEGA*2.0/TAU)
         write(*,*)
         NTRYCOUNT = 0
      ELSE IF(NTRYCOUNT.GE.3) THEN
         WRITE(16,*) 'NO CRITICAL RAYLEIGH NUMBER FOUND.'
         STOP NO_RA_FOUND
      ELSE
         write(*,*) 'NO CRIT. RAYLEIGH NUMBER FOUND. Trying again.'
         NTRYCOUNT = NTRYCOUNT + 1
      ENDIF

      close(16)
      end subroutine

***********************************************************************
      subroutine fixedParCriticalEigenVector(outputfile)
      IMPLICIT double precision(A-H,O-Y)
      IMPLICIT complex(8)(Z)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      COMMON/PAR2/ETA
      COMMON/DIM/NT,M0,ND,LMIN,LD,NE
      COMMON/NUM/TTA,TTF,TTSTEP,DRA,ABSE,RELE,NSMAX,LCALC
      CHARACTER(len=*) outputfile
      complex(8) ZEVEC(NMAX)
      COMMON/MODE/ZEVEC
         open(16,FILE=outputfile,STATUS='UNKNOWN')
C--      searching for zero grothrate:
         CALL pegasus(RA,DRA,ABSE,RELE,NSMAX,1,VERBOSE,LL,RAC)
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
      end subroutine

***********************************************************************
      SUBROUTINE losub(inputfile,outputfile)
      IMPLICIT double precision(A-H,O-Y)
      IMPLICIT complex(8)(Z)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      COMMON/PAR2/ETA
      COMMON/DIM/NT,M0,ND,LMIN,LD,NE
      COMMON/NUM/TTA,TTF,TTSTEP,DRA,ABSE,RELE,NSMAX,LCALC
      CHARACTER(len=40) infile,outfile
      COMMON/FILES/infile,outfile
      CHARACTER(len=*) inputfile,outputfile
      INTEGER nloop
      complex(8) ZEVEC(NMAX)
      COMMON/MODE/ZEVEC

      infile = inputfile
      outfile = ''
      outfile = outputfile

C-----Default values:
      CALL setDefaults()
      RAOLD = 4.D3
      RconcOLD = 4.D3
      TAU1 = 100.0D0
      TAU2 = 100.0D0

C-----INPUT:
      CALL readInputFile(inputfile)

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
      CALL writeOutputHeader(outputfile)
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
      CALL DIMENSION(LMIN,LD,NT,M0,ND)
      write(*,*) 'DIMENSION OF MATRIX:',ND
C
c---LCALC=-1 : most basic case: find the most unstable growth rate at all parameters fixed
      IF( LCALC.EQ.-1) THEN
         CALL fixedParGrowthRate(outputfile)
         stop ' Growth rate at fixed other parametres'
c---LCALC=0 : Critical Ra, for constant other parameters
      ELSE IF( LCALC.EQ.0 ) THEN
         CALL fixedParCriticalRa(outputfile)
         stop ' Ra_crit at fixed other parametres'
C-----eigenvalues determined for this value of RA: ------------------
      ELSE IF( LCALC.EQ.1 ) THEN
         GROR=FN(RA)
         RETURN
C     LCALC-s which requiree increment of TAU:
      ELSE IF(( LCALC.EQ.2 ).or.( LCALC.EQ.3 )) then
         NTRYCOUNT = 0
         nloop = 0
         do
            TA = TAU*TAU
C-----------searching for zero grothrate by varying RA: ----------------------
            IF( LCALC.EQ.2 ) THEN
               CALL fixedParCriticalRaAndM(outputfile, RAC, NTRYCOUNT)
C-----------searching for zero grothrate by variing RA and M0: ---------------
            ELSE IF( LCALC.EQ.3 ) THEN
               CALL fixedParCriticalRaAndM0(outputfile, RAC, NTRYCOUNT)
            ENDIF 
C--         increment TAU:
            TAU0 = TAU1
            TAU1 = TAU
            IF( DABS(TTSTEP) .LT. DABS(TAU*0.1D0) ) THEN
               TAU   = TAU1 + TTSTEP
            ELSE
               TAU = TAU1 + TTSTEP/DABS(TTSTEP)*TAU1*0.1D0
            ENDIF 
C--         interpolate new startingvalue for RA:
            IF(TAU1.NE.TAU0) THEN
C              RA = RAC + (RAC-RAOLD)/(TAU1-TAU0)*(TAU-TAU1)
C             Assume Rac grows with Ta**(4/3)               
              RA = RAC + (4.0/3.0)*Tau1**(1.0/3.0)*(TAU-TAU1)
            ELSE
              RA = 2.0D0*RAC - RAOLD 
            ENDIF
            RAOLD = RAC
          
C--         endvalue of TAU reached?
            IF ( ((TAU.GT.TTF).AND.(TTF.GT.TTA)).OR.
     &           ((TAU.LT.TTF).AND.(TTF.LT.TTA)) )  THEN
c            WRITE(*,*) 'LOSUB.F: finished at ',fdate()
               STOP FINISHED
            ENDIF
C--------------count loops:
            nloop=nloop+1
C--         End of tau Loop
         enddo
C-----calculate the critical eigenvector. Print for plotting.
      ELSE IF (LCALC.EQ.4) THEN
         CALL fixedParCriticalEigenVector(outputfile)
C-----vary m and calculate critical R at fixed P, tau, eta.
      ELSE IF (LCALC.EQ.5) THEN
         M0A=M0
         do M0=M0A, int(TTF), INT(TTSTEP)
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
            CALL dimension(LMIN,LD,NT,M0,ND)
            DRA = RA/10.0D0
            CALL pegasus(RA,DRA,ABSE,RELE,NSMAX,1,QUIET,LL,RAC)
            WRITE(*,* ) M0,RAC
            CALL open_file_at_end(16,outputfile)
            WRITE(16,*) M0,RAC
            close(16)
         enddo
      ELSE IF (LCALC.EQ.6) THEN
C-----vary Le and calculate critical R at fixed P, tau, eta, M
         pL0=pL
         niter=(TTF-pL0)/TTSTEP
         do nloop=0,niter
            pL=pL0+nloop*TTSTEP
            DRA = RA/10.0D0
            CALL pegasus(RA,DRA,ABSE,RELE,NSMAX,1,QUIET,LL,RAC)
            GROR=FN(RAC)
            WRITE(*,*) pL,RAC,GROR
            CALL open_file_at_end(16,outputfile)
            WRITE(16,'(3D16.8)') pL,RAC,GROR
            close(16)
         enddo
      end if

      END

************************************************************************
C-- FUNCTION CALLED BY PEGASUS, FN=GROTHRATE.
      FUNCTION FN(RA)
      IMPLICIT double precision(A-H,O-Y)
      IMPLICIT complex(8)(Z)
      DIMENSION ZA(NMAX,NMAX),ZB(NMAX,NMAX)
      DIMENSION ZACOPY(NMAX,NMAX),ZBCOPY(NMAX,NMAX)
      DIMENSION ZEW(NMAX),ZEWA(NMAX),ZEWB(NMAX)
      DIMENSION ZEVALL(NMAX,NMAX),ZEVEC(NMAX)
      COMMON/PAR/TAU,RAL,PR,RI,C,pL,Rconc
      COMMON/DIM/NT,M0,ND,LMIN,LD,NE
      COMMON/MODE/ZEVEC
      COMMON/NUM/TTA,TTF,TTSTEP,DRA,ABSE,RELE,NSMAX,LCALC

      RAL=RA

C-- MAT SETS THE complex(8) MATRICES ZA AND ZB SETTING OF MATRIX:
      CALL MAT(ZA,ZB,NMAX)

      ZEWA   = DCMPLX(0D0,0D0)
      ZEWB   = DCMPLX(0D0,0D0)
      ZACOPY = DCMPLX(0D0,0D0)
      ZBCOPY = DCMPLX(0D0,0D0)
C
C-- DGVLCG IS AN IMSL ROUTINE THAT CALCULATES THE EIGENVALUES ZEWA/ZEWB:
C-- Problem: generalized complex(8) eigensystem A*x = lam*B*x
C--          Input:  A: ZA, B: ZB
C--          Output: complex(8) ZEWA(ND), ZEWB(ND)
C      CALL DGVLCG(ND,ZA,NMAX,ZB,NMAX,ZEWA,ZEWB)
      IF(LCALC.EQ.4) THEN
        CALL DG2CCG (ND,ZA,NMAX,ZB,NMAX,ZEWA,ZEWB,ZEVALL,
     &                           NMAX,ZACOPY,ZBCOPY)
      ELSE
        CALL DG2LCG(ND,ZA,NMAX,ZB,NMAX,ZEWA,ZEWB,ZACOPY,ZBCOPY)
      ENDIF
C
      ZEW(:)=ZEWA(:)/ZEWB(:)
C---- search for lowest imaginary part:
      IMIN  = minloc(DIMAG(ZEW),1)
      EWMIN = DIMAG(ZEW(IMIN))
    
c      IF( EWMIN.LT.0.D0 ) WRITE(*,*) 'CONVECTION REGIME !',RAL
      C=-DREAL(ZEW(IMIN))/M0
      FN = EWMIN
C
      IF(LCALC.EQ.1) THEN
C------- sort eigenvalues:
         write(*,*) '     Frequ.(exp(+iwt))   -Grothrate  '
         DO I=1,ND
            DO J=I,ND
               IF( DIMAG(ZEW(J)).LT.DIMAG(ZEW(I)) ) THEN
                  ZSAVE = ZEW(J)
                  ZEW(J) = ZEW(I)
                  ZEW(I) = ZSAVE
                  DO K=1,ND
                     ZEVEC(K)    = ZEVALL(K,J)
                     ZEVALL(K,J) = ZEVALL(K,I)
                     ZEVALL(K,I) = ZEVEC(K)
                  enddo
               ENDIF 
            enddo
            WRITE(*,'(I4,2D16.6)') I,ZEW(I)     
         enddo
      ELSEIF(LCALC.EQ.4) THEN
C---------save crit. eigenvector:
         DO J=1,ND
            ZEVEC(J) = ZEVALL(J,IMIN)
         enddo
      ENDIF
      END
C
*************************************************************************
C-- SETS THE complex(8) MATRICES ZA AND ZB.
      SUBROUTINE MAT(ZA,ZB,NDIM)
      IMPLICIT double precision(A-H,O-Y)
      IMPLICIT complex(8)(Z)
      DIMENSION ZA(ndim,ndim),ZB(ndim,ndim)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      COMMON/DIM/NT,M0,ND,LMIN,LD,NE
C
      ZA(:,:)=DCMPLX(0D0,0D0)
      ZB(:,:)=DCMPLX(0D0,0D0)
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

         NIMAX=DINT( DBLE(2*NT+1-LI+M0)/2 )
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
               ENDIF
               NJMAX=DINT( DBLE(2*NT+1-LJ+M0)/2 )
C
C  ******************** I: Equation (Line) ******************
C  ******************** J: Variable (Column) ****************
C  ******************** I+1: v (poloidal)  ******************
C  ******************** I+2: theta         ******************
C  ******************** I+3: w (toroidal)  ******************
C-new****************** I+4: gamma (concentration) ********** 
               DO 1000 NJ=1,NJMAX
               IF(J+3.GT.NDIM .OR. I+3.GT.NDIM) THEN
                  write(*,*) 'MAT(): NDIM too small.'
                  stop
               ENDIF
                  IF( LI.EQ.LJ ) THEN
                      ZB(I+1,J+1)=DCMPLX(0.D0,-DIII2(NI,NJ,LPI,1))
                ZA(I+1,J+1)= DCMPLX(DIII1(NI,NJ,LPI),DIII3(NI,NJ,LPI,1))
                     ZA(I+1,J+2)=DCMPLX(DIII5(NI,NJ,LPI),0.D0)
c--- concentration driving
                     ZA(I+1,J+4)=DCMPLX(DIII5conc(NI,NJ,LPI),0.D0)
c
                      ZB(I+2,J+2)=DCMPLX(0.D0,-DI1(NI,NJ,1))
                     ZA(I+2,J+1)=DCMPLX(DI3(NI,NJ,LPI),0.D0)
                     ZA(I+2,J+2)=DCMPLX(DI2(NI,NJ,LPI),0.D0)
                      ZB(I+3,J+3)=DCMPLX(0.D0,-DII2(NI,NJ,LTI,1))
                     ZA(I+3,J+3)=DCMPLX(DII1(NI,NJ,LTI),DII3(NI,NJ,1))
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
                  J=J+4
1000           CONTINUE
2000        CONTINUE
             I=I+4
3000     CONTINUE
4000  CONTINUE
      END
C
*************************************************************************
C-- GALERKIN TERMS:
*************************************************************************
      FUNCTION DI1(N1,N2,NU1)
C----- HEAT EQUATION, TIME DERIVATIVE
         IMPLICIT double precision(A-H,O-Z)
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
         DI1=PR*NU1*R('SS ',2,N1,N2,0)
      END
C
      FUNCTION DI2(N1,N2,L1)
C----- HEAT EQUATION , DISSIPATION
         IMPLICIT double precision(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
         DI2=N2**2*DPI**2*R('SS ',2,N1,N2,0) - 
     -           2*N2*DPI*R('SC ',1,N1,N2,0) +
     +              DL(L1)*R('SS ',0,N1,N2,0)  
      END
C
      FUNCTION DI3(N1,N2,L1)
C----- HEAT EQUATION , SOURCE 
         IMPLICIT double precision(A-H,O-Z)
         DI3 =-DL(L1)*R('SS ',2,N1,N2,0)
      END
C
      FUNCTION DII1(N1,N2,L1)
C----- TOROIDAL EQUATION , DISSIPATION
         IMPLICIT double precision(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
          DII1=DL(L1)*( (N2-1)**2*DPI**2*R('CC ',4,N1-1,N2-1,0) +
     +                     4*(N2-1)*DPI*R('CS ',3,N1-1,N2-1,0) + 
     +                       (DL(L1)-2)*R('CC ',2,N1-1,N2-1,0) )
      END
C
      FUNCTION DII2(N1,N2,L1,NU1)
C----- TOROIDAL EQUATION , TIME DERIVATIVE
         IMPLICIT double precision(A-H,O-Z)
          DII2=NU1*DL(L1)*R('CC ',4,N1-1,N2-1,0)
      END
C
      FUNCTION DII3(N1,N2,NU1)
C----- TOROIDAL EQUATION , CORRIOLIS
         IMPLICIT double precision(A-H,O-Z)
         COMMON/DIM/NT,M0,ND,LMIN,LD,NE
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
         DII3=-TAU*NU1*M0*R('CC ',4,N1-1,N2-1,0) 
      END
C
      FUNCTION DII4A(N1,N2,L1,NU1)
C----- TOROIADL EQUATION , Q-TERM 1 (L1=L3+1)
         IMPLICIT double precision(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
         COMMON/DIM/NT,M0,ND,LMIN,LD,NE
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
          DII4A= TAU * 
     *        DSQRT( DBLE(L1-NU1*M0)*(L1+NU1*M0)/(2*L1-1)/(2*L1+1) ) * 
     *                (    (L1**2-1)*(L1-1)*R('CS ',2,N1-1,N2,0)  - 
     -            (L1+1)*(L1-1)*N2*DPI*R('CC ',3,N1-1,N2,0)    ) 
      END
C
      FUNCTION DII4B(N1,N2,L1,NU1)
C----- TOROIADL EQUATION , Q-TERM 1 (L1=L3-1)
         IMPLICIT double precision(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
         COMMON/DIM/NT,M0,ND,LMIN,LD,NE
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
          DII4B= TAU * 
     *     DSQRT( DBLE(L1-NU1*M0+1)*(L1+NU1*M0+1)/(2*L1+1)/(2*L1+3) ) * 
     *                ( (1-(L1+1)**2)*(L1+2)*R('CS ',2,N1-1,N2,0)  - 
     -                 L1*(L1+2)*N2*DPI*R('CC ',3,N1-1,N2,0)  ) 
      END
C
      FUNCTION DIII1(N1,N2,L1)
C----- POLOIDAL EQUOATION , DISSIPATION
         IMPLICIT double precision(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
         DIII1=DL(L1)* ( N2**4*DPI**4*R('SS ',2,N1,N2,0) - 
     -                  4*N2**3*DPI**3*R('SC ',1,N1,N2,0) +
     +      2*DL(L1)*N2**2*DPI**2*R('SS ',0,N1,N2,0) +
     +           (DL(L1)**2-2*DL(L1))*R('SS ',-2,N1,N2,0) )
      END
C
      FUNCTION DIII2(N1,N2,L1,NU1)
C----- POLOIDAL EQUATION , TIME DERIVATIVE
         IMPLICIT double precision(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
         DIII2=        -NU1*DL(L1)*(-N2**2*DPI**2*R('SS ',2,N1,N2,0) +
     +         2*N2*DPI*R('SC ',1,N1,N2,0)-DL(L1)*R('SS ',0,N1,N2,0)  ) 
      END
C
      FUNCTION DIII3(N1,N2,L1,NU1)
C----- POLOIDAL EQUATION , CORRIOLIS
         IMPLICIT double precision(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
         COMMON/DIM/NT,M0,ND,LMIN,LD,NE
          DIII3= TAU*NU1*M0*( -N2**2*DPI**2*R('SS ',2,N1,N2,0) +
     +                             2*N2*DPI*R('SC ',1,N1,N2,0) -
     -                     DL(L1)*R('SS ',0,N1,N2,0)  ) 
      END
C
      FUNCTION DIII4A(N1,N2,L1,NU1)
C----- POLOIDAL EUQUATION , Q-TERM 1 (L1=L3+1)
         IMPLICIT double precision(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
         COMMON/DIM/NT,M0,ND,LMIN,LD,NE
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
          DIII4A= TAU * 
     *      DSQRT(DBLE(L1-M0*NU1)*(L1+M0*NU1)/(2*L1-1)/(2*L1+1)) *
     *          (     (L1*(L1-1)-2)*(L1-1)*R('SC ',2,N1,N2-1,0) +
     +             (L1+1)*(L1-1)*(N2-1)*DPI*R('SS ',3,N1,N2-1,0)    )
      END
C
      FUNCTION DIII4B(N1,N2,L1,NU1)
C----- POLOIDAL EQUATION , Q-TERM 2 (L1=L3-1)
         IMPLICIT double precision(A-H,O-Z)
         PARAMETER(DPI=3.141592653589793D0)
         COMMON/DIM/NT,M0,ND,LMIN,LD,NE
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
         DIII4B= TAU * 
     *       DSQRT( DBLE(L1-M0*NU1+1)*(L1+M0*NU1+1)/(2*L1+1)/(2*L1+3) )*
     *             (    (L1+2)*(2-(L1+1)*(L1+2))*R('SC ',2,N1,N2-1,0) + 
     +         L1*(L1+2)*(N2-1)*DPI*R('SS ',3,N1,N2-1,0)   )
      END
C
      FUNCTION DIII5(N1,N2,L1)
C----- POLOIDAL EQUATION , 
         IMPLICIT double precision(A-H,O-Z)
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
         DIII5=-RA*DL(L1)*R('SS ',2,N1,N2,0)
      END

      FUNCTION DIII5conc(N1,N2,L1)
C----- POLOIDAL EQUATION , 
         IMPLICIT double precision(A-H,O-Z)
         COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
         DIII5conc=-Rconc*DL(L1)*R('SS ',2,N1,N2,0)
      END
C
***************************************************************************
*-- SUBROUTINES:
***************************************************************************
      FUNCTION DL(L)
         IMPLICIT double precision(A-H,O-Z)
         DL = DBLE(L*(L+1))
      END

****************************************************************************
      SUBROUTINE DIMENSION(LMIN,LD,NT,M0,ND)
****************************************************************************
      integer NT,M0,ND,LMIN,LD
      integer L
C-- DETERMINATION OF DIMENSION:
C-- for each value of L the number of possible N-values is added
c         print*, "Triangular truncation (2.12)"
c         print*, LMIN, "...", 2*NT+M0-1,LD
      ND=0 
      DO L = LMIN, 2*NT+M0-1, LD
c         print*, L, 1, "...", INT( DBLE(2*NT+1-L+M0)/2 )
cccccccccc18    ND=ND+3*DINT( DBLE(2*NT+1-L+M0)/2 )
         ND = ND + 4*DINT( DBLE(2*NT+1-L+M0)/2 )
      ENDDO
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
      CHARACTER(len=*) filename

      OPEN(NHANDLE,FILE=filename,STATUS='OLD',POSITION='APPEND',ERR=990)
      GOTO 999
990   WRITE(*,*) 'Error reading ',filename
      STOP ERR_WRT_OUTFILE
999   CONTINUE
      END

*************************************************************************
      SUBROUTINE abort1
*     this function should replace the internal FORTRAN routine
*     'abort' which is e.g. CALLed from the IMSL routines
*************************************************************************
      IMPLICIT double precision(A-H,O-Y)
      COMMON/PAR/TAU,RA,PR,RI,C,pL,Rconc
      COMMON/NUM/TTA,TTF,TTSTEP,DRA,ABSE,RELE,NSMAX,LCALC

      write(*,*) 'LOSUB.F: abort(fortran) CALLed.'
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

      CALL ende(SIGRESTART)
      RETURN
      END

*************************************************************************
      SUBROUTINE ende(II)
*     ende() is CALLed explicit in the program or implicit when the
*     process receives a signal.
*     II indicates the behaviour of ende() or the received signal number.
*************************************************************************
      IMPLICIT double precision(A-H,O-Y)
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

c      IF(II.EQ.SIGUSR1) THEN
c       write(*,'(A,A)') ' ende(): received signal SIGUSR1 at ',fdate()
c      ELSEIF(II.EQ.SIGABRT) THEN
c       write(*,'(A,A)') ' ende(): received signal SIGABRT at ',fdate()
c      ENDIF

      INQUIRE(FILE=trim(outfile),OPENED=LOP,NUMBER=NHANDLE)
      IF(LOP) THEN
        CLOSE(NHANDLE)
      ENDIF

      IF(II.EQ.SIGRESTART .OR. II.EQ.SIGUSR1) THEN
         WRITE(*,*) 'LOSUB.F: Terminating this run and starting next.'
         OPEN(15,FILE=infile,STATUS='UNKNOWN')
         WRITE(15,*) ' NE (0/1/2) | LCALC (1/2/3/4) |'
         WRITE(15,'(A,2I12)') ' ',NE,LCALC
         WRITE(15,*) '|  RAYLEIGH  |  TAU     |  PRANTEL  |  ETA  |
     &    Lewis |   Rconc   |'
         WRITE(15,'(1P,E17.6,5(A,E17.6))') 
     &         RA,' ',TAU,' ',PR,' ',ETA,' ',pL,' ',Rconc
         WRITE(15,*) '|   NTRUNC (>=1) | MODE |'
         WRITE(15,'(A,2I12)') ' ',NT,M0
         WRITE(15,*) '|   DRA   | ABSERR  |  RELERR  | NMAX |'
         WRITE(15,'(1PG12.6,A,1PG11.5,A,1PG11.5,A,I4)') 
     &           DRA,' ',ABSE,' ',RELE,' ',NSMAX
         WRITE(15,*) '|   TAU_STEP | TAU_END'
         WRITE(15,'(1P,2G11.4)') TTSTEP,TTF
         CLOSE(15)
         stop 10
C--------return value START_NEXT_RUN to the shell script: 
c         CALL exit(START_NEXT_RUN) 
      ELSE
c         CALL exit(1)
      ENDIF

      END subroutine
