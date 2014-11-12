module linearProblem
  use parameters
  implicit none
contains
   SUBROUTINE fixedParGrowthRate (outputfile) 
      IMPLICIT none 
      CHARACTER (len = * ):: outputfile 
      INTEGER:: nloop, niter 

      Ra0 = RA 
      niter = int ( (TTF - Ra0) / TTSTEP) 
      CALL open_file_at_end (16, outputfile) 
      DO nloop = 0, niter 
         RA = Ra0 + nloop * TTSTEP 
         GROR = FN (RA) 
         WRITE ( *, * ) RA, GROR 
         WRITE (16, '(3D16.8)') RA, GROR 
      enddo 
      CLOSE (16) 
                                                                        
      GROR = FN (RA) 
      WRITE (*,*) 'R=', RA
      WRITE (*,*) 'TAU=', TAU, ' P=', PR, ' M0=', M0, 'eta=', ETA                                                        
      WRITE (*,*) 'Most unstable growth rate', GROR 
      WRITE (*,*) 'If growth rate < 0 then above onset' 
   end subroutine 
                                                                        
!********************************************************************** 
   SUBROUTINE fixedParCriticalRa (outputfile) 
      IMPLICIT none
      CHARACTER(len=*)::outputfile 
      CALL PEGASUS (RA, DRA, ABSE, RELE, NSMAX, 1, 1, LL, RAC) 
      GROR = FN(RAC) 
      OMEGA = C * M0 
      CALL open_file_at_end (16, outputfile) 
      WRITE (*,*) 'TAU=', TAU, ' P=', PR, ' M0=', M0, ' eta=', ETA 
      WRITE (*,*) 'Lewis=', pL, ' Rconc=', Rconc 
      WRITE (*,*) 'R_crit=', RAC, '  (growth rate =', GROR, ')' 
      CLOSE (16) 
   endsubroutine 
                                                                        
!********************************************************************** 
   SUBROUTINE fixedParCriticalRaAndM (outputfile, NTRYCOUNT) 
      IMPLICIT doubleprecision (A - H, O - Y) 
      IMPLICIT complex (8) (Z) 
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc 
      COMMON / PAR2 / ETA 
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE 
      COMMON / NUM / TTA, TTF, TTSTEP, DRA, ABSE, RELE, NSMAX, LCALC 
      INTEGER NTRYCOUNT 
      CHARACTER (len = * ) outputfile 
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
      ENDIF 
      CLOSE (16) 
      endsubroutine 
                                                                        
!********************************************************************** 
      SUBROUTINE fixedParCriticalRaAndM0 (outputfile, NTRYCOUNT) 
      IMPLICIT doubleprecision (A - H, O - Y) 
      IMPLICIT complex (8) (Z) 
      COMMON / PAR / TAU, RA, PR, RI, C, pL, Rconc 
      COMMON / PAR2 / ETA 
      COMMON / DIM / NT, M0, ND, LMIN, LD, NE 
      COMMON / NUM / TTA, TTF, TTSTEP, DRA, ABSE, RELE, NSMAX, LCALC 
      CHARACTER (len = * ) outputfile 
      DOUBLEPRECISION RACI (3), OMI (3) 
      INTEGER M0I (3), LLI (3), LMINI (3) 
      INTEGER NTRYCOUNT 
      COMPLEX (8) ZEVEC (NMAX) 
      COMMON / MODE / ZEVEC 
                                                                        
      IF (M0.GT.1) THEN 
         M0I (1) = M0 - 1 
         LMINI (1) = LMIN - 1 
      ELSE 
         M0I (1) = M0 
         LMINI (1) = LMIN 
      ENDIF 
      M0I (2) = M0 
      M0I (3) = M0 + 1 
      LMINI (2) = LMIN 
      LMINI (3) = LMIN + 1 
      RAC = RA 
      DO II = 1, 3 
      M0 = M0I (II) 
      LMIN = LMINI (II) 
      DRA = RA / 10.0D0 
!         Print*, Ra, dRa, abse, rele, nsmax, ll, rac                   
!        LL is either 0 for success or 1 for failure                    
      CALL PEGASUS (RA, DRA, ABSE, RELE, NSMAX, 1, 0, LL, RAC) 
      GROR = FN (RAC) 
      OMI (II) = C * M0 
      RACI (II) = RAC 
      LLI (II) = LL 
      WRITE ( * , '(X,1P,4E17.6,I4,A3,2E17.6)') TAU, RAC, OMI (II) ,    &
      GROR, M0, ' | ', RAC * (1 / (1 - eta) ) **4 * (2 / TAU) , (OMI (  &
      II) * 2.0 / TAU)                                                  
      enddo 
!      CALL dimension(LMIN,LD,NT,M0,ND)                                 
                                                                        
      INDEX = 0 
      IF (LLI (1) .EQ.0.AND.LLI (2) .EQ.0) THEN 
         IF (RACI (1) .LT.RACI (2) ) THEN 
            INDEX = 1 
         ELSE 
            INDEX = 2 
         ENDIF 
      ELSEIF (LLI (1) .EQ.0) THEN 
         INDEX = 1 
      ELSEIF (LLI (2) .EQ.0) THEN 
         INDEX = 2 
      ENDIF 
      IF (LLI (3) .EQ.0.AND.INDEX.GT.0) THEN 
         IF (RACI (3) .LT.RACI (INDEX) ) THEN 
            INDEX = 3 
         ENDIF 
      ELSEIF (LLI (3) .EQ.0) THEN 
         INDEX = 3 
      ENDIF 
                                                                        
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
         WRITE ( *, * ) 
         NTRYCOUNT = 0 
      ELSEIF (NTRYCOUNT.GE.3) THEN 
         WRITE (16, * ) 'NO CRITICAL RAYLEIGH NUMBER FOUND.' 
         STOP NO_RA_FOUND 
      ELSE 
      WRITE ( * ,  * ) 'NO CRIT. RAYLEIGH NUMBER FOUND. Trying again.' 
         NTRYCOUNT = NTRYCOUNT + 1 
      ENDIF 
                                                                        
      CLOSE (16) 
   endsubroutine 
                                                                        
!********************************************************************** 
   SUBROUTINE fixedParCriticalEigenVector (outputfile) 
      IMPLICIT doubleprecision (A - H, O - Y) 
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
         Write(16,*) 'C               OM, WHERE?, FLOQUET'
         WRITE (16, '(2D17.10,3I4)') C0, OMM, NUC, NUOM, MF 
 9100    FORMAT 
                                                                        
         LMAX = 2 * NT + M0 - 1 
         I = 0 
         DO LI = LMIN, LMAX, LD 
!           L for poloidal (v) field:                                   
            LPI = LI 
            NIMAX = DINT (DBLE (2 * NT + 1 - LI + M0) / 2) 
            DO NI = 1, NIMAX 
               IF (LST.EQ.0) THEN 
                  WRITE (16, '(1X,A1,4I3,4D16.8)') 'V', LPI, M0, NI, 0, DREAL (ZEVEC (I + 1) ) , DIMAG (ZEVEC (I + 1) ) , 0.D0, 0.D0          
               ELSEIF (LST.EQ.3) THEN 
                  WRITE (16, '(A,4I3,A,2F11.7,A)') ' ''V ''', LPI, M0, NI, 0, ' ', DREAL (ZEVEC (I + 1) ) , DIMAG (ZEVEC (I + 1) ) , ' .0D+00 .0D+00 '                              
               ENDIF 
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
            ENDIF 
            NIMAX = DINT (DBLE (2 * NT + 1 - LI + M0) / 2) 
            DO NI = 1, NIMAX 
               IF (LST.EQ.0) THEN 
                  WRITE (16, '(1X,A1,4I3,4D16.8)') 'W', LTI, M0, NI, 0, DREAL (ZEVEC (I + 3) ) , DIMAG (ZEVEC (I + 3) ) , 0.D0, 0.D0          
               ELSEIF (LST.EQ.3) THEN 
                  WRITE (16, '(A,4I3,A,2F11.7,A)') ' ''W ''', LTI, M0,  NI, 0, ' ', DREAL (ZEVEC (I + 3) ) , DIMAG (ZEVEC (I + 3) ) , ' .0D+00 .0D+00 '                              
               ENDIF 
               I = I + 4 
            enddo
         enddo
!                                                                       
         I = 0 
         DO LI = LMIN, LMAX, LD 
            NIMAX = DINT (DBLE (2 * NT + 1 - LI + M0) / 2) 
            DO NI = 1, NIMAX 
               IF (LST.EQ.0) THEN 
                  WRITE (16, '(1X,A1,4I3,4D16.8)') 'T', LI, M0, NI, 0, DREAL (ZEVEC (I + 2) ) , DIMAG (ZEVEC (I + 2) ) , 0.D0, 0.D0            
               ELSEIF (LST.EQ.3) THEN 
                  WRITE (16, '(A,4I3,A,2F11.7,A)') ' ''T ''', LI, M0,   NI, 0, ' ', DREAL (ZEVEC (I + 2) ) , DIMAG (ZEVEC (I + 2) ) , ' .0D+00 .0D+00 '                              
               ENDIF 
               I = I + 4 
            enddo
         enddo
!                                                                       
         I = 0 
         DO LI = LMIN, LMAX, LD 
            NIMAX = DINT (DBLE (2 * NT + 1 - LI + M0) / 2) 
            DO NI = 1, NIMAX 
               IF (LST.EQ.0) THEN 
                  WRITE (16, '(1X,A1,4I3,4D16.8)') 'G', LI, M0, NI, 0, DREAL (ZEVEC (I + 4) ) , DIMAG (ZEVEC (I + 4) ) , 0.D0, 0.D0            
               ELSEIF (LST.EQ.3) THEN 
                  WRITE (16, '(A,4I3,A,2F11.7,A)') ' ''G ''', LI, M0,   NI, 0, ' ', DREAL (ZEVEC (I + 4) ) , DIMAG (ZEVEC (I + 4) ) , ' .0D+00 .0D+00 '                              
               ENDIF 
               I = I + 4 
            enddo
         enddo
      ELSE 
         WRITE (16, * ) 'NO CRITICAL RAYLEIGH NUMBER FOUND.' 
         STOP NO_RA_FOUND 
      ENDIF 
      endsubroutine 

end module linearProblem
