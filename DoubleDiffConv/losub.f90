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
!#ifndef NMAX                                                           
!#define NMAX 144                                                       
!#endif                                                                 
                                                                        
!#ifdef __hpux                                                          
!#define fdate fdate_                                                   
!#endif                                                                 
                                                                        
!******** limit for CPU time in sec, for C routine clock(): ************
!#ifdef __hpux                                                          
!******  medium queue:   4 h = 14000 sec                                
!******  long   queue:  16 h = 57600 sec                                
!#define CPUTIMELIM         52000.0D0                                   
!#else                                                                  
!#define CPUTIMELIM         40000.0D0                                   
!#endif                                                                 
!#define CLOCKS_PER_SECOND  1000000.0D0                                 
!#define NCPUMAXVALUE       2147483647                                  
                                                                        
!****** DQS 3.0: Signal, when process exceeds CPU time limit: **********
!****** (doesn't work actually, because signal is send only to          
!******  shell, not to the working process.                             
!******  Only the Bourneshell can catch the signal, but reacts          
!******  not before one shell command is finished (especially the       
!******  working process).                                              
!#define SIGRESTART  80                                                 
                                                                        
!******  defines from <sys/signal.h> *********************              
!#  define SIGABRT       6       /* Process abort signal */             
!#  define SIGUSR1       16      /* user defined signal 1 */            
!#  define SIGUSR2       17      /* user defined signal 2 */            
                                                                        
!****** Process return values (2...255) ********************************
!#define NO_INFILE      100                                             
!#define ERR_IN_INFILE  101                                             
!#define DIM_TO_SMALL   102                                             
!#define ERR_WRT_OUTFILE 103                                            
!#define NO_RA_FOUND    120                                             
!#define START_NEXT_RUN 190                                             
!#define FINISHED       199                                             
                                                                        
!#define DOUBLEMAXVALUE  1.0D306                                        
!********************************************************************** 
                                                                        
!********************************************************************** 
                                                                        
!********************************************************************** 
      SUBROUTINE losub (inputfile, outputfile) 
      use parameters
      use io
      implicit none
      CHARACTER (len = * ) inputfile, outputfile 
      INTEGER:: nloop 
      COMPLEX(8):: ZEVEC (NMAX) 
      COMMON / MODE / ZEVEC 
                                                                        
!-----Default values:                                                   
      CALL setDefaults () 
                                                                        
!-----INPUT:                                                            
      CALL readInputFile (inputfile) 
                                                                        
!-----LOSUB.F doesn't work for M=0 !!!!!!                               
      IF (M0.LT.1) THEN 
      WRITE ( * ,  * ) 'The code does not work for M0<1.', M0, ' --> 1' 
         M0 = 1 
      ENDIF 
                                                                        
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
      ENDIF 
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
            CALL fixedParCriticalRaAndM (outputfile, NTRYCOUNT) 
!-----------searching for zero grothrate by variing RA and M0: ---------
         ELSEIF (LCALC.EQ.3) THEN 
            CALL fixedParCriticalRaAndM0 (outputfile, NTRYCOUNT) 
         ENDIF 
!--         increment TAU:                                              
         TAU0 = TAU1 
         TAU1 = TAU 
         IF (DABS (TTSTEP) .LT.DABS (TAU * 0.1D0) ) THEN 
            TAU = TAU1 + TTSTEP 
         ELSE 
            TAU = TAU1 + TTSTEP / DABS (TTSTEP) * TAU1 * 0.1D0 
         ENDIF 
!--         interpolate new startingvalue for RA:                       
         IF (TAU1.NE.TAU0) THEN 
            RA = RAC + (RAC - RAOLD) / (TAU1 - TAU0) * (TAU - TAU1) 
         ELSE 
            RA = 2.0D0 * RAC - RAOLD 
         ENDIF 
         RAOLD = RAC 
                                                                        
!--      endvalue of TAU reached?                                       
         IF ( ( (TAU.GT.TTF) .AND. (TTF.GT.TTA) ) .OR. ( (TAU.LT.TTF)   &
         .AND. (TTF.LT.TTA) ) ) THEN                                    
!              WRITE(*,*) 'LOSUB.F: finished at ',fdate()               
            STOP FINISHED 
         ENDIF 
                                                                        
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
         IF (NE.EQ.0) THEN 
!--   UNDEFINED SYMMETRIE:                                              
            LMIN = M0 
            LD = 1 
         ELSEIF (NE.EQ.1) THEN 
!--   EQUATORIAL ANTISYMMETRIE (L+M ODD):                               
            LMIN = M0 + 1 
            LD = 2 
         ELSEIF (NE.EQ.2) THEN 
!--   EQUATORIAL SYMMETRIE (L+M EVEN);                                  
            LMIN = M0 
            LD = 2 
         ENDIF 
         CALL dimension (LMIN, LD, NT, M0, ND) 
         DRA = RA / 10.0D0 
         CALL PEGASUS (RA, DRA, ABSE, RELE, NSMAX, 1, 0, LL, RAC) 
         WRITE ( *, * ) M0, RAC 
         CALL open_file_at_end (16, outfile) 
         WRITE (16, * ) M0, RAC 
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
         WRITE ( *, * ) pL, RAC, GROR 
         CALL open_file_at_end (16, outfile) 
         WRITE (16, '(3D16.8)') pL, RAC, GROR 
         CLOSE (16) 
         enddo 
      ENDIF 
                                                                        
                                                                        
                                                                        
                                                                        
      END SUBROUTINE losub                          
                                                                        
!***********************************************************************
!-- FUNCTION CALLED BY PEGASUS, FN=GROTHRATE.                           
      FUNCTION FN (RA) 
      IMPLICIT doubleprecision (A - H, O - Y) 
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
                                                                        
      DO 50 I = 1, NMAX 
         ZEWA (I) = DCMPLX (0D0, 0D0) 
         ZEWB (I) = DCMPLX (0D0, 0D0) 
         DO 50 J = 1, NMAX 
            ZACOPY (I, J) = DCMPLX (0D0, 0D0) 
            ZBCOPY (I, J) = DCMPLX (0D0, 0D0) 
   50 CONTINUE 
!                                                                       
!-- DGVLCG IS AN IMSL ROUTINE THAT CALCULATES THE EIGENVALUES ZEWA/ZEWB:
!-- Problem: generalized complex(8) eigensystem A*x = lam*B*x           
!--          Input:  A: ZA, B: ZB                                       
!--          Output: complex(8) ZEWA(ND), ZEWB(ND)                      
!      CALL DGVLCG(ND,ZA,NMAX,ZB,NMAX,ZEWA,ZEWB)                        
      IF (LCALC.EQ.4) THEN 
         CALL DG2CCG (ND, ZA, NMAX, ZB, NMAX, ZEWA, ZEWB, ZEVALL, NMAX, &
         ZACOPY, ZBCOPY)                                                
      ELSE 
         CALL DG2LCG (ND, ZA, NMAX, ZB, NMAX, ZEWA, ZEWB, ZACOPY,       &
         ZBCOPY)                                                        
      ENDIF 
!                                                                       
      DO 100 I = 1, ND 
  100 ZEW (I) = ZEWA (I) / ZEWB (I) 
      IMIN = 1 
!---- search for lowest imaginary part:                                 
      EWMIN = DIMAG (ZEW (IMIN) ) 
      DO 200 I = 2, ND 
         IF (DIMAG (ZEW (I) ) .LT.EWMIN) THEN 
            IMIN = I 
            EWMIN = DIMAG (ZEW (IMIN) ) 
         ENDIF 
  200 END DO 
                                                                        
!      IF( EWMIN.LT.0.D0 ) WRITE(*,*) 'CONVECTION REGIME !',RAL         
      C = - DREAL (ZEW (IMIN) ) / M0 
!                                                                       
      IF (LCALC.EQ.1) THEN 
!---- sort eigenvalues:                                                 
      WRITE ( * ,  * ) '     Frequ.(exp(+iwt))   -Grothrate  ' 
         DO 320 I = 1, ND 
            DO 310 J = I, ND 
               IF (DIMAG (ZEW (J) ) .LT.DIMAG (ZEW (I) ) ) THEN 
                  ZSAVE = ZEW (J) 
                  ZEW (J) = ZEW (I) 
                  ZEW (I) = ZSAVE 
                  DO K = 1, ND 
                  ZEVEC (K) = ZEVALL (K, J) 
                  ZEVALL (K, J) = ZEVALL (K, I) 
                  ZEVALL (K, I) = ZEVEC (K) 
                  enddo 
               ENDIF 
  310       END DO 
  320    WRITE ( * , '(I4,2D16.6)') I, ZEW (I) 
      ENDIF 
                                                                        
!------save crit. eigenvector:                                          
      IF (LCALC.EQ.4) THEN 
         DO J = 1, ND 
         ZEVEC (J) = ZEVALL (J, IMIN) 
         enddo 
      ENDIF 
                                                                        
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
      IMPLICIT doubleprecision (A - H, O - Y) 
      IMPLICIT complex (8) (Z) 
!                                                                       
      DIMENSION ZA (NDIM, * ), ZB (NDIM, * ) 
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc 
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE 
!                                                                       
      DO 200 J = 1, NDIM 
         DO 200 I = 1, NDIM 
            ZA (I, J) = DCMPLX (0D0, 0D0) 
            ZB (I, J) = DCMPLX (0D0, 0D0) 
  200 CONTINUE 
!                                                                       
      I = 0 
      LMAX = 2 * NT + M0 - 1 
!                                                                       
      DO 4000 LI = LMIN, LMAX, LD 
         LPI = LI 
!                                                                       
!        Determine L for toroidal (w) field:                            
         IF (NE.EQ.2) THEN 
            LTI = LI + 1 
         ELSEIF (NE.EQ.1) THEN 
            LTI = LI - 1 
         ELSEIF (NE.EQ.0) THEN 
            LTI = LI 
         ENDIF 
!                                                                       
         NIMAX = DINT (DBLE (2 * NT + 1 - LI + M0) / 2) 
!                                                                       
         DO 3000 NI = 1, NIMAX 
!                                                                       
            J = 0 
            DO 2000 LJ = LMIN, LMAX, LD 
               LPJ = LJ 
!                                                                       
               IF (NE.EQ.2) THEN 
                  LTJ = LJ + 1 
               ELSEIF (NE.EQ.1) THEN 
                  LTJ = LJ - 1 
!               ELSEIF( NE.EQ.0 ) THEN                                  
               ELSE 
                  LTJ = LJ 
               ENDIF 
!                                                                       
               NJMAX = DINT (DBLE (2 * NT + 1 - LJ + M0) / 2) 
!                                                                       
!  ******************** I: Equation (Line) ******************           
!  ******************** J: Variable (Column) ****************           
!  ******************** I+1: v (poloidal)  ******************           
!  ******************** I+2: theta         ******************           
!  ******************** I+3: w (toroidal)  ******************           
!-new****************** I+4: gamma (concentration) **********           
               DO 1000 NJ = 1, NJMAX 
!                                                                       
                  IF (J + 3.GT.NDIM.OR.I + 3.GT.NDIM) THEN 
                     WRITE ( * , * ) 'MAT(): NDIM too small.' 
                     STOP 
                  ENDIF 
!                                                                       
                  IF (LI.EQ.LJ) THEN 
                     ZB (I + 1, J + 1) = DCMPLX (0.D0, - DIII2 (NI, NJ, &
                     LPI, 1) )                                          
                     ZA (I + 1, J + 1) = DCMPLX (DIII1 (NI, NJ, LPI),   &
                     DIII3 (NI, NJ, LPI, 1) )                           
                     ZA (I + 1, J + 2) = DCMPLX (DIII5 (NI, NJ, LPI),   &
                     0.D0)                                              
!--- concentration driving                                              
                     ZA (I + 1, J + 4) = DCMPLX (DIII5conc (NI, NJ, LPI)&
                     , 0.D0)                                            
!                                                                       
                     ZB (I + 2, J + 2) = DCMPLX (0.D0, - DI1 (NI, NJ, 1)&
                     )                                                  
                     ZA (I + 2, J + 1) = DCMPLX (DI3 (NI, NJ, LPI),     &
                     0.D0)                                              
                     ZA (I + 2, J + 2) = DCMPLX (DI2 (NI, NJ, LPI),     &
                     0.D0)                                              
                     ZB (I + 3, J + 3) = DCMPLX (0.D0, - DII2 (NI, NJ,  &
                     LTI, 1) )                                          
                     ZA (I + 3, J + 3) = DCMPLX (DII1 (NI, NJ, LTI),    &
                     DII3 (NI, NJ, 1) )                                 
!--- concentration equation                                             
                     ZB (I + 4, J + 4) = DCMPLX (0.D0, - DI1 (NI, NJ, 1)&
                     )                                                  
                     ZA (I + 4, J + 1) = DCMPLX (DI3 (NI, NJ, LPI),     &
                     0.D0)                                              
                     ZA (I + 4, J + 4) = DCMPLX (1.D0 / pL * DI2 (NI,   &
                     NJ, LPI), 0.D0)                                    
                  ENDIF 
                  IF (LPI.EQ.LTJ + 1) THEN 
                     ZA (I + 1, J + 3) = DCMPLX (DIII4A (NI, NJ, LPI, 1)&
                     , 0.D0)                                            
                  ELSEIF (LPI.EQ.LTJ - 1) THEN 
                     ZA (I + 1, J + 3) = DCMPLX (DIII4B (NI, NJ, LPI, 1)&
                     , 0.D0)                                            
                  ENDIF 
                  IF (LTI.EQ.LPJ + 1) THEN 
                     ZA (I + 3, J + 1) = DCMPLX (DII4A (NI, NJ, LTI, 1),&
                     0.D0)                                              
                  ELSEIF (LTI.EQ.LPJ - 1) THEN 
                     ZA (I + 3, J + 1) = DCMPLX (DII4B (NI, NJ, LTI, 1),&
                     0.D0)                                              
                  ENDIF 
                  J = J + 4 
 1000          END DO 
!                                                                       
 2000       END DO 
!                                                                       
            I = I + 4 
 3000    END DO 
!                                                                       
 4000 END DO 
      END SUBROUTINE MAT                            
!                                                                       
!***********************************************************************
!-- GALERKIN TERMS:                                                     
!***********************************************************************
      FUNCTION DI1 (N1, N2, NU1) 
!----- HEAT EQUATION, TIME DERIVATIVE                                   
      IMPLICIT doubleprecision (A - H, O - Z) 
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc 
      DI1 = PR * NU1 * R ('SS ', 2, N1, N2, 0) 
      END FUNCTION DI1                              
!                                                                       
      FUNCTION DI2 (N1, N2, L1) 
!----- HEAT EQUATION , DISSIPATION                                      
      IMPLICIT doubleprecision (A - H, O - Z) 
      PARAMETER (DPI = 3.141592653589793D0) 
      DI2 = N2**2 * DPI**2 * R ('SS ', 2, N1, N2, 0) - 2 * N2 * DPI * R &
      ('SC ', 1, N1, N2, 0) + DL (L1) * R ('SS ', 0, N1, N2, 0)         
      END FUNCTION DI2                              
!                                                                       
      FUNCTION DI3 (N1, N2, L1) 
!----- HEAT EQUATION , SOURCE                                           
      IMPLICIT doubleprecision (A - H, O - Z) 
      DI3 = - DL (L1) * R ('SS ', 2, N1, N2, 0) 
      END FUNCTION DI3                              
!                                                                       
      FUNCTION DII1 (N1, N2, L1) 
!----- TOROIDAL EQUATION , DISSIPATION                                  
      IMPLICIT doubleprecision (A - H, O - Z) 
      PARAMETER (DPI = 3.141592653589793D0) 
      DII1 = DL (L1) * ( (N2 - 1) **2 * DPI**2 * R ('CC ', 4, N1 - 1,   &
      N2 - 1, 0) + 4 * (N2 - 1) * DPI * R ('CS ', 3, N1 - 1, N2 - 1, 0) &
      + (DL (L1) - 2) * R ('CC ', 2, N1 - 1, N2 - 1, 0) )               
      END FUNCTION DII1                             
!                                                                       
      FUNCTION DII2 (N1, N2, L1, NU1) 
!----- TOROIDAL EQUATION , TIME DERIVATIVE                              
      IMPLICIT doubleprecision (A - H, O - Z) 
      DII2 = NU1 * DL (L1) * R ('CC ', 4, N1 - 1, N2 - 1, 0) 
      END FUNCTION DII2                             
!                                                                       
      FUNCTION DII3 (N1, N2, NU1) 
!----- TOROIDAL EQUATION , CORRIOLIS                                    
      IMPLICIT doubleprecision (A - H, O - Z) 
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE 
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc 
      DII3 = - TAU * NU1 * M0 * R ('CC ', 4, N1 - 1, N2 - 1, 0) 
      END FUNCTION DII3                             
!                                                                       
      FUNCTION DII4A (N1, N2, L1, NU1) 
!----- TOROIADL EQUATION , Q-TERM 1 (L1=L3+1)                           
      IMPLICIT doubleprecision (A - H, O - Z) 
      PARAMETER (DPI = 3.141592653589793D0) 
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE 
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc 
      DII4A = TAU * DSQRT (DBLE (L1 - NU1 * M0) * (L1 + NU1 * M0)       &
      / (2 * L1 - 1) / (2 * L1 + 1) ) * ( (L1**2 - 1) * (L1 - 1)        &
      * R ('CS ', 2, N1 - 1, N2, 0) - (L1 + 1) * (L1 - 1) * N2 * DPI *  &
      R ('CC ', 3, N1 - 1, N2, 0) )                                     
      END FUNCTION DII4A                            
!                                                                       
      FUNCTION DII4B (N1, N2, L1, NU1) 
!----- TOROIADL EQUATION , Q-TERM 1 (L1=L3-1)                           
      IMPLICIT doubleprecision (A - H, O - Z) 
      PARAMETER (DPI = 3.141592653589793D0) 
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE 
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc 
      DII4B = TAU * DSQRT (DBLE (L1 - NU1 * M0 + 1) * (L1 + NU1 * M0 +  &
      1) / (2 * L1 + 1) / (2 * L1 + 3) ) * ( (1 - (L1 + 1) **2) *       &
      (L1 + 2) * R ('CS ', 2, N1 - 1, N2, 0) - L1 * (L1 + 2) * N2 * DPI &
      * R ('CC ', 3, N1 - 1, N2, 0) )                                   
      END FUNCTION DII4B                            
!                                                                       
      FUNCTION DIII1 (N1, N2, L1) 
!----- POLOIDAL EQUOATION , DISSIPATION                                 
      IMPLICIT doubleprecision (A - H, O - Z) 
      PARAMETER (DPI = 3.141592653589793D0) 
      DIII1 = DL (L1) * (N2**4 * DPI**4 * R ('SS ', 2, N1, N2, 0)       &
      - 4 * N2**3 * DPI**3 * R ('SC ', 1, N1, N2, 0) + 2 * DL (L1)      &
      * N2**2 * DPI**2 * R ('SS ', 0, N1, N2, 0) + (DL (L1) **2 - 2 *   &
      DL (L1) ) * R ('SS ', - 2, N1, N2, 0) )                           
      END FUNCTION DIII1                            
!                                                                       
      FUNCTION DIII2 (N1, N2, L1, NU1) 
!----- POLOIDAL EQUATION , TIME DERIVATIVE                              
      IMPLICIT doubleprecision (A - H, O - Z) 
      PARAMETER (DPI = 3.141592653589793D0) 
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc 
      DIII2 = - NU1 * DL (L1) * ( - N2**2 * DPI**2 * R ('SS ', 2, N1,   &
      N2, 0) + 2 * N2 * DPI * R ('SC ', 1, N1, N2, 0) - DL (L1) * R (   &
      'SS ', 0, N1, N2, 0) )                                            
      END FUNCTION DIII2                            
!                                                                       
      FUNCTION DIII3 (N1, N2, L1, NU1) 
!----- POLOIDAL EQUATION , CORRIOLIS                                    
      IMPLICIT doubleprecision (A - H, O - Z) 
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
      IMPLICIT doubleprecision (A - H, O - Z) 
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
      IMPLICIT doubleprecision (A - H, O - Z) 
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
      IMPLICIT doubleprecision (A - H, O - Z) 
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc 
      DIII5 = - RA * DL (L1) * R ('SS ', 2, N1, N2, 0) 
      END FUNCTION DIII5                            
                                                                        
      FUNCTION DIII5conc (N1, N2, L1) 
!----- POLOIDAL EQUATION ,                                              
      IMPLICIT doubleprecision (A - H, O - Z) 
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc 
      DIII5conc = - Rconc * DL (L1) * R ('SS ', 2, N1, N2, 0) 
      END FUNCTION DIII5conc                        
!                                                                       
!***********************************************************************
!-- SUBROUTINES:                                                        
!***********************************************************************
      FUNCTION DL (L) 
      IMPLICIT doubleprecision (A - H, O - Z) 
      DL = DBLE (L * (L + 1) ) 
      END FUNCTION DL                               
!                                                                       
!-- DUMMY FUNCTIONS:                                                    
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
      ENDDO 
!                                                                       
      IF (ND.GT.NMAX) THEN 
         WRITE ( * , * ) 'DIMENSION OF MATRIX TOO SMALL:', ND, '>',     &
         NMAX                                                           
         STOP DIM_TO_SMALL 
      ENDIF 
                                                                        
      END SUBROUTINE DIMENSION                      
                                                                        
!***********************************************************************
      SUBROUTINE open_file_at_end (NHANDLE, filename) 
!***********************************************************************
!     opens file <filename> and puts the filepointer at EOF             
!***********************************************************************
      INTEGER NHANDLE 
      CHARACTER(40) filename 
                                                                        
      OPEN (NHANDLE, FILE = trim (filename) , STATUS = 'OLD', POSITION =&
      'APPEND', ERR = 990)                                              
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
      IMPLICIT doubleprecision (A - H, O - Y) 
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc 
      COMMON / NUM / TTA, TTF, TTSTEP, DRA, ABSE, RELE, NSMAX, LCALC 
                                                                        
      WRITE ( * , * ) 'LOSUB.F: abort(fortran) CALLed.' 
!-----increment TAU:                                                    
      IF (DABS (TTSTEP) .LT.DABS (TAU * 0.1D0) ) THEN 
         TAU = TAU + TTSTEP 
      ELSE 
         TAU = TAU + TAU * 0.1D0 
      ENDIF 
                                                                        
!--   endvalue of TAU reached?                                          
      IF ( ( (TAU.GT.TTF) .AND. (TTF.GT.TTA) ) .OR. ( (TAU.LT.TTF)      &
      .AND. (TTF.LT.TTA) ) ) THEN                                       
!      WRITE(*,*) 'LOSUB.F: finished at ',fdate()                       
         STOP FINISHED 
      ENDIF 
                                                                        
      CALL ende (SIGRESTART) 
      RETURN 
      END SUBROUTINE abort1                         
                                                                        
!***********************************************************************
      SUBROUTINE ende (II) 
!     ende() is CALLed explicit in the program or implicit when the     
!     process receives a signal.                                        
!     II indicates the behaviour of ende() or the received signal number
!***********************************************************************
      IMPLICIT doubleprecision (A - H, O - Y) 
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc 
      COMMON / PAR2 / ETA 
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE 
      COMMON / NUM / TTA, TTF, TTSTEP, DRA, ABSE, RELE, NSMAX, LCALC 
      CHARACTER(40) infile, outfile 
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
!      ENDIF                                                            
                                                                        
      INQUIRE (FILE = trim (outfile), OPENED = LOP, NUMBER = NHANDLE) 
      IF (LOP) THEN 
         CLOSE (NHANDLE) 
      ENDIF 
                                                                        
      IF (II.EQ.SIGRESTART.OR.II.EQ.SIGUSR1) THEN 
      WRITE ( * ,  * ) 'LOSUB.F: Terminating this run and starting next.&
     &'                                                                 
         OPEN (15, FILE = infile, STATUS = 'UNKNOWN') 
         WRITE (15, * ) ' NE (0/1/2) | LCALC (1/2/3/4) |' 
         WRITE (15, '(A,2I12)') ' ', NE, LCALC 
      WRITE (15,  * ) '|  RAYLEIGH  |  TAU     |  PRANTEL  |  ETA  |    &
     &     Lewis |   Rconc   |'                                         
         WRITE (15, '(1P,E17.6,5(A,E17.6))') RA, ' ', TAU, ' ', PR, ' ',&
         ETA, ' ', pL, ' ', Rconc                                       
      WRITE (15,  * ) '|   NTRUNC (>=1) | MODE |' 
         WRITE (15, '(A,2I12)') ' ', NT, M0 
      WRITE (15,  * ) '|   DRA   | ABSERR  |  RELERR  | NMAX |' 
         WRITE (15, '(1PG12.6,A,1PG11.5,A,1PG11.5,A,I4)') DRA, ' ',     &
         ABSE, ' ', RELE, ' ', NSMAX                                    
      WRITE (15,  * ) '|   TAU_STEP | TAU_END' 
         WRITE (15, '(1P,2G11.4)') TTSTEP, TTF 
         CLOSE (15) 
         STOP 10 
!--------return value START_NEXT_RUN to the shell script:               
!         CALL exit(START_NEXT_RUN)                                     
      ELSE 
!         CALL exit(1)                                                  
      ENDIF 
                                                                        
      END SUBROUTINE ende                           
