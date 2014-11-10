!***********************************************************************
!  IMSL.F                                                               
!  collection of IMSL routines for LOSUB.F                              
!  attention: routine DG4CCG is modified    MA 23.02.95                 
!********************************************************************** 
!                                                                       
!-----------------------------------------------------------------------
!  IMSL Name:  G2CCG/DG2CCG (Single/Double precision version)           
!                                                                       
!  Computer:   sgruxs/DOUBLE                                            
!                                                                       
!  Revised:    September 25, 1985                                       
!                                                                       
!  Purpose:    Compute all of the eigenvalues and eigenvectors of a     
!              generalized complex(8) eigensystem A*z = w*B*z.          
!                                                                       
!  Usage:      CALL G2CCG (N, A, LDA, B, LDB, ALPHA, BETA, EVEC, LDEVEC,
!                          ACOPY, BCOPY)                                
!                                                                       
!  Arguments:  (See GVCCG)                                              
!                                                                       
!  Chapter:    MATH/LIBRARY Eigensystem Analysis                        
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE DG2CCG (N, A, LDA, B, LDB, ALPHA, BETA, EVEC, LDEVEC,  &
      ACOPY, BCOPY)                                                     
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N, LDA, LDB, LDEVEC 
      COMPLEX (8) A (LDA, * ), B (LDB, * ), ALPHA ( * ), BETA ( * ),    &
      EVEC (LDEVEC, * ), ACOPY (N, * ), BCOPY (N, * )                   
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, J 
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL ZSCAL, E1MES, E1POP, E1PSH, E1STI, DCCGCG, DG3CCG,       &
      DG4CCG, DG5CCG                                                    
!                                  SPECIFICATIONS FOR FUNCTIONS         
      EXTERNAL IZAMAX, N1RTY 
      INTEGER IZAMAX, N1RTY 
!      INTRINSIC  DIMAG                                                 
!                                                                       
      CALL E1PSH ('DG2CCG ') 
!                                  Check N                              
      IF (N.LT.1) THEN 
         CALL E1STI (1, N) 
      CALL E1MES (5, 1, 'The argument N = %(I1).  The '//'order of the m&
     &atrix must be at least 1.')                                       
         GOTO 9000 
      ENDIF 
!                                  Check LDA                            
      IF (LDA.LT.N) THEN 
         CALL E1STI (1, LDA) 
         CALL E1STI (2, N) 
      CALL E1MES (5, 2, 'The argument LDA = %(I1).  The '//'leading dime&
     &nsion of the matrix A must be at '//'least equal to the order, N =&
     & %(I2).')                                                         
      ENDIF 
!                                  Check LDB                            
      IF (LDB.LT.N) THEN 
         CALL E1STI (1, LDB) 
         CALL E1STI (2, N) 
      CALL E1MES (5, 3, 'The argument LDB = %(I1).  The '//'leading dime&
     &nsion of the matrix B must be at '//'least equal to the order N = &
     &%(I2).')                                                          
      ENDIF 
!                                  Check LDEVEC                         
      IF (LDEVEC.LT.N) THEN 
         CALL E1STI (1, LDEVEC) 
         CALL E1STI (2, N) 
      CALL E1MES (5, 4, 'The argument LDEVEC = %(I1).  The '//'order of &
     &the eigenvector matrix must be at '//'least equal to the order N =&
     & %(I2).')                                                         
      ENDIF 
      IF (N1RTY (0) .GT.0) GOTO 9000 
!                                  Copy A and B                         
      CALL DCCGCG (N, A, LDA, ACOPY, N) 
      CALL DCCGCG (N, B, LDB, BCOPY, N) 
!                                  Reduce                               
      CALL DG3CCG (N, ACOPY, BCOPY, 1, EVEC, LDEVEC) 
!                                  Find ALPHA, BETA and EVEC            
      CALL DG4CCG (N, ACOPY, BCOPY, EVEC, LDEVEC, 1, ALPHA, BETA) 
!                                  Order ALPHA and BETA                 
      CALL DG5CCG (N, N, ALPHA, BETA, EVEC, LDEVEC, .TRUE., ACOPY (1, 1)&
      )                                                                 
!                                  Normalize eigenvectors               
      DO 10 J = 1, N 
         I = IZAMAX (N, EVEC (1, J), 1) 
         IF (DBLE (EVEC (I, J) ) .NE.0.0D0.OR.DIMAG (EVEC (I, J) )      &
         .NE.0.0D0) THEN                                                
            CALL ZSCAL (N, 1.0D0 / EVEC (I, J), EVEC (1, J), 1) 
         ENDIF 
   10 END DO 
!                                                                       
 9000 CALL E1POP ('DG2CCG ') 
      RETURN 
      END SUBROUTINE DG2CCG                         
!-----------------------------------------------------------------------
!  IMSL Name:  G2LCG/DG2LCG (Single/Double precision version)           
!                                                                       
!  Computer:   sgruxs/DOUBLE                                            
!                                                                       
!  Revised:    September 26, 1985                                       
!                                                                       
!  Purpose:    Compute all of the eigenvalues of a generalized complex(8
!              eigensystem A*z = w*B*z.                                 
!                                                                       
!  Usage:      CALL G2LCG (N, A, LDA, B, LDB, ALPHA, BETA, ACOPY, BCOPY)
!                                                                       
!  Arguments:  (See GVLCG)                                              
!                                                                       
!  Chapter:    MATH/LIBRARY Eigensystem Analysis                        
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE DG2LCG (N, A, LDA, B, LDB, ALPHA, BETA, ACOPY, BCOPY) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N, LDA, LDB 
      COMPLEX (8) A (LDA, * ), B (LDB, * ), ALPHA ( * ), BETA ( * ),    &
      ACOPY (N, * ), BCOPY (N, * )                                      
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      COMPLEX (8) EVEC (1, 1) 
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL E1MES, E1POP, E1PSH, E1STI, DCCGCG, DG3CCG, DG4CCG,      &
      DG5CCG                                                            
!                                  SPECIFICATIONS FOR FUNCTIONS         
      EXTERNAL N1RTY 
      INTEGER N1RTY 
!      INTRINSIC  DIMAG                                                 
                                                                        
!                                                                       
      CALL E1PSH ('DG2LCG ') 
!                                  Check N                              
      IF (N.LT.1) THEN 
         CALL E1STI (1, N) 
      CALL E1MES (5, 1, 'The argument N = %(I1).  The '//'order of the m&
     &atrix must be at least 1.')                                       
         GOTO 9000 
      ENDIF 
!                                  Check LDA                            
      IF (LDA.LT.N) THEN 
         CALL E1STI (1, LDA) 
         CALL E1STI (2, N) 
      CALL E1MES (5, 2, 'The argument LDA = %(I1).  The '//'leading dime&
     &nsion of the matrix A must be at '//'least equal to the order, N =&
     & %(I2).')                                                         
      ENDIF 
!                                  Check LDB                            
      IF (LDB.LT.N) THEN 
         CALL E1STI (1, LDB) 
         CALL E1STI (2, N) 
      CALL E1MES (5, 3, 'The argument LDB = %(I1).  The '//'leading dime&
     &nsion of the matrix B must be at '//'least equal to the order, N =&
     & %(I2).')                                                         
      ENDIF 
      IF (N1RTY (0) .GT.0) GOTO 9000 
!                                  Copy A and B                         
      CALL DCCGCG (N, A, LDA, ACOPY, N) 
      CALL DCCGCG (N, B, LDB, BCOPY, N) 
!                                  Reduce                               
      CALL DG3CCG (N, ACOPY, BCOPY, 0, EVEC, 1) 
!                                  Find ALPHA and BETA                  
      CALL DG4CCG (N, ACOPY, BCOPY, EVEC, 1, 0, ALPHA, BETA) 
!                                  Order ALPHA and BETA                 
      CALL DG5CCG (N, N, ALPHA, BETA, EVEC, 1, .FALSE., ACOPY (1, 1) ) 
!                                                                       
 9000 CALL E1POP ('DG2LCG ') 
      RETURN 
      END SUBROUTINE DG2LCG                         
!-----------------------------------------------------------------------
!  IMSL Name:  CCGCG/DCCGCG (Single/Double precision version)           
!                                                                       
!  Computer:   sgruxs/DOUBLE                                            
!                                                                       
!  Revised:    June 5, 1985                                             
!                                                                       
!  Purpose:    Copy a complex(8) general matrix.                        
!                                                                       
!  Usage:      CALL CCGCG (N, A, LDA, B, LDB)                           
!                                                                       
!  Arguments:                                                           
!     N      - Order of the matrices A and B.  (Input)                  
!     A      - complex(8) matrix of order N.  (Input)                   
!     LDA    - Leading dimension of A exactly as specified in the       
!              dimension statement of the calling program.  (Input)     
!     B      - complex(8) matrix of order N containing a copy of A.     
!              (Output)                                                 
!     LDB    - Leading dimension of B exactly as specified in the       
!              dimension statement of the calling program.  (Input)     
!                                                                       
!  GAMS:       D1b8                                                     
!                                                                       
!  Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations              
!                                                                       
!  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE DCCGCG (N, A, LDA, B, LDB) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N, LDA, LDB 
      COMPLEX (8) A (LDA, * ), B (LDB, * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER J 
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL ZCOPY, E1MES, E1POP, E1PSH, E1STI 
!                                                                       
      CALL E1PSH ('DCCGCG ') 
!                                  Check N                              
      IF (N.LT.1) THEN 
         CALL E1STI (1, N) 
      CALL E1MES (5, 1, 'The argument N = %(I1).  It must be at '//'leas&
     &t 1.')                                                            
         GOTO 9000 
      ENDIF 
!                                  Check LDA                            
      IF (LDA.LT.N) THEN 
         CALL E1STI (1, LDA) 
         CALL E1STI (2, N) 
      CALL E1MES (5, 2, 'The argument LDA = %(I1).  It must be at '//'le&
     &ast as large as N = %(I2).')                                      
         GOTO 9000 
      ENDIF 
!                                  Check LDB                            
      IF (LDB.LT.N) THEN 
         CALL E1STI (1, LDB) 
         CALL E1STI (2, N) 
      CALL E1MES (5, 3, 'The argument LDB = %(I1).  It must be at '//'le&
     &ast as large as N = %(I2).')                                      
         GOTO 9000 
      ENDIF 
!                                  Copy                                 
      IF (LDA.EQ.N.AND.LDB.EQ.N) THEN 
         CALL ZCOPY (N * N, A, 1, B, 1) 
      ELSEIF (LDA.GE.LDB) THEN 
         DO 10 J = 1, N 
            CALL ZCOPY (N, A (1, J), 1, B (1, J), 1) 
   10    END DO 
      ELSE 
         DO 20 J = N, 1, - 1 
            CALL ZCOPY (N, A (1, J), - 1, B (1, J), - 1) 
   20    END DO 
      ENDIF 
!                                                                       
 9000 CONTINUE 
      CALL E1POP ('DCCGCG ') 
      RETURN 
      END SUBROUTINE DCCGCG                         
!-----------------------------------------------------------------------
!  IMSL Name:  G3CCG/DG3CCG (Single/Double precision version)           
!                                                                       
!  Computer:   sgruxs/DOUBLE                                            
!                                                                       
!  Revised:    May 21, 1985                                             
!                                                                       
!  Purpose:    Reduce the complex(8) matrix A to upper Hessenberg form a
!              the complex(8) matrix B to upper triangular form.        
!                                                                       
!  Usage:      CALL G3CCG (N, A, B, IJOB, EVEC, LDEVEC)                 
!                                                                       
!  Arguments:                                                           
!     N      - Order of the system.  (Input)                            
!     A      - complex(8) matrix of order and leading dimension N.      
!              (Input/Output)                                           
!     B      - complex(8) matrix of order and leading dimension N.      
!              (Input/Output)                                           
!     IJOB   - If IJOB = 1 then the transformation is accumulated in    
!              EVEC.  Otherwise EVEC is not used.  (Input)              
!     EVEC   - complex(8) matrix of order N.  If IJOB = 1 then it contai
!              the reduction transformation.  (Output)                  
!     LDEVEC - Leading dimension of EVEC exactly as specified in the    
!              dimension statement of the calling program.  (Input)     
!                                                                       
!  Chapter:    MATH/LIBRARY Eigensystem Analysis                        
!                                                                       
!  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE DG3CCG (N, A, B, IJOB, EVEC, LDEVEC) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N, IJOB, LDEVEC 
      COMPLEX (8) A (N, * ), B (N, * ), EVEC (LDEVEC, * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, II, J 
      COMPLEX (8) Y, ZDUM 
!                                  SPECIFICATIONS FOR INTRINSICS        
!     INTRINSIC  DABS,DIMAG,DBLE                                        
!      INTRINSIC  DABS, DIMAG, DBLE                                     
!      DOUBLE PRECISION DABS, DIMAG, DBLE                               
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL ZAXPY, ZSET, ZSWAP 
!                                  SPECIFICATIONS FOR FUNCTIONS         
      EXTERNAL IZAMAX 
      INTEGER IZAMAX 
      DOUBLEPRECISION CABS1 
!                                                                       
      CABS1 (ZDUM) = DABS (DBLE (ZDUM) ) + DABS (DIMAG (ZDUM) ) 
!                                                                       
      IF (N.EQ.1) GOTO 9000 
!                                  Reduce B to triangular form using    
!                                  elementary transformations           
      DO 20 I = 1, N - 1 
         II = I + IZAMAX (N - I, B (I + 1, I), 1) 
         IF (DBLE (B (II, I) ) .EQ.0.0D0.AND.DIMAG (B (II, I) )         &
         .EQ.0.0D0) GOTO 20                                             
         IF (CABS1 (B (II, I) ) .GT.CABS1 (B (I, I) ) ) THEN 
!                                  Must interchange                     
            CALL ZSWAP (N, A (I, 1), N, A (II, 1), N) 
            CALL ZSWAP (N - I + 1, B (I, I), N, B (II, I), N) 
         ENDIF 
         DO 10 J = I + 1, N 
            Y = B (J, I) / B (I, I) 
            IF (DBLE (Y) .NE.0.0D0.OR.DIMAG (Y) .NE.0.0D0) THEN 
               CALL ZAXPY (N, - Y, A (I, 1), N, A (J, 1), N) 
               CALL ZAXPY (N - I, - Y, B (I, I + 1), N, B (J, I + 1),   &
               N)                                                       
            ENDIF 
   10    END DO 
         B (I + 1, I) = 0.0D0 
   20 END DO 
!                                  Initialize EVEC                      
      IF (IJOB.EQ.1) THEN 
         DO 30 J = 1, N 
            CALL ZSET (N, (0.0D0, 0.0D0), EVEC (1, J), 1) 
   30    END DO 
         CALL ZSET (N, (1.0D0, 0.0D0), EVEC, LDEVEC + 1) 
      ENDIF 
!                                  Reduce A to upper Hessenberg form    
      DO 50 J = 1, N - 2 
         DO 40 I = N, J + 2, - 1 
            IF (CABS1 (A (I, J) ) .GT.CABS1 (A (I - 1, J) ) ) THEN 
!                                  Must interchange rows                
               CALL ZSWAP (N - J + 1, A (I, J), N, A (I - 1, J),        &
               N)                                                       
               CALL ZSWAP (N - I + 2, B (I, I - 1), N, B (I - 1, I - 1),&
               N)                                                       
            ENDIF 
            IF (DBLE (A (I, J) ) .NE.0.0D0.OR.DIMAG (A (I, J) )         &
            .NE.0.0D0) THEN                                             
               Y = A (I, J) / A (I - 1, J) 
               CALL ZAXPY (N - J, - Y, A (I - 1, J + 1), N, A (I, J + 1)&
               , N)                                                     
               CALL ZAXPY (N - I + 2, - Y, B (I - 1, I - 1), N, B (I, I &
               - 1), N)                                                 
            ENDIF 
!                                  Transformation from the right        
            IF (CABS1 (B (I, I - 1) ) .GT.CABS1 (B (I, I) ) ) THEN 
!                                  Must interchange columns             
               CALL ZSWAP (I, B (1, I), 1, B (1, I - 1), 1) 
               CALL ZSWAP (N, A (1, I), 1, A (1, I - 1), 1) 
               IF (IJOB.EQ.1) THEN 
                  CALL ZSWAP (N - I + J + 1, EVEC (I - J, I), 1, EVEC ( &
                  I - J, I - 1), 1)                                     
               ENDIF 
            ENDIF 
            IF (DBLE (B (I, I - 1) ) .NE.0.0D0.OR.DIMAG (B (I, I - 1) ) &
            .NE.0.0D0) THEN                                             
               Y = B (I, I - 1) / B (I, I) 
               CALL ZAXPY (I - 1, - Y, B (1, I), 1, B (1, I - 1),       &
               1)                                                       
               B (I, I - 1) = 0.0D0 
               CALL ZAXPY (N, - Y, A (1, I), 1, A (1, I - 1), 1) 
               IF (IJOB.EQ.1) THEN 
                  CALL ZAXPY (N - I + J + 1, - Y, EVEC (I - J, I),      &
                  1, EVEC (I - J, I - 1), 1)                            
               ENDIF 
            ENDIF 
   40    END DO 
         A (J + 2, J) = 0.0D0 
   50 END DO 
 9000 CONTINUE 
      RETURN 
      END SUBROUTINE DG3CCG                         
!-----------------------------------------------------------------------
!  IMSL Name:  G4CCG/DG4CCG (Single/Double precision version)           
!                                                                       
!********************************************************************** 
!  modified version :                                                   
!  24.02.95 MA       calling ende(NSIGRESTART) when iteration failed    
!********************************************************************** 
!  Computer:   sgruxs/DOUBLE                                            
!                                                                       
!  Revised:    May 21, 1985                                             
!                                                                       
!  Purpose:    Compute eigenvalues and eigenvectors of a generalized    
!              complex(8) eigensystem A*z = w*B*z given its reduction to
!              an upper Hessenberg and upper triangular matrix.         
!                                                                       
!  Usage:      CALL G4CCG (N, A, B, EVEC, LDEVEC, IJOB, ALPHA, BETA)    
!                                                                       
!  Arguments:                                                           
!     N      - Order of the system.  (Input)                            
!     A      - complex(8) upper Hessenberg matrix.  (Input)             
!     B      - complex(8) upper triangular matrix.  (Input)             
!     EVEC   - complex(8) matrix of order N.  If IJOB = 1, on input, the
!              transformation matrix.  On output, containing the        
!              eigenvectors.   (Input/Output)                           
!     LDEVEC - Leading dimension of EVEC exactly as specified in the    
!              dimension statement of the calling program.  (Input)     
!     IJOB   - If IJOB = 0 then only the eigenvalues are computed.  If  
!              IJOB = 1 then the eigenvectors are also computed.        
!              (Input)                                                  
!     ALPHA  - complex(8) vector of length N containing the numerators o
!              the eigenvalues.  (Output)                               
!     BETA   - complex(8) vector of length N containing the demoninators
!              of the eigenvalues.  (Output)                            
!                                                                       
!  Chapter:    MATH/LIBRARY Eigensystem Analysis                        
!                                                                       
!  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE DG4CCG (N, A, B, EVEC, LDEVEC, IJOB, ALPHA, BETA) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      PARAMETER (NSIGRESTART = 80) 
!                                                                       
      INTEGER N, LDEVEC, IJOB 
      COMPLEX (8) A (N, * ), B (N, * ), EVEC (LDEVEC, * ), ALPHA ( * ), &
      BETA ( * )                                                        
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, INFER, ITS, J, L, LOR1, M, NN, NNORN 
      DOUBLEPRECISION ANI, ANORM, BNI, BNORM, D0, D1, D2, E0, E1, EPS,  &
      EPSA, EPSB, R, SS                                                 
      COMPLEX (8) ALFM, ANM1M1, ANNM1, BETM, D, DEN, NUM, S, SL, W, Y,  &
      ZC, ZDUM, TEMP1                                                   
!                                  SPECIFICATIONS FOR INTRINSICS        
!      INTRINSIC  DABS,DIMAG,DMAX1,DCMPLX,CDSQRT,DBLE                   
!      INTRINSIC  DABS, DIMAG, DMAX1, DCMPLX, CDSQRT, DBLE              
!      DOUBLE PRECISION DABS, DIMAG, DMAX1, DBLE                        
!      complex(8) DCMPLX, CDSQRT                                        
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL ZAXPY, ZSWAP, E1MES, E1POP, E1PSH, E1STI 
!                                  SPECIFICATIONS FOR FUNCTIONS         
      EXTERNAL DMACH, ZDOTU, IZAMAX, DZASUM 
      EXTERNAL ende 
      INTEGER IZAMAX 
      DOUBLEPRECISION DMACH, DZASUM 
      COMPLEX (8) ZDOTU 
      DOUBLEPRECISION CABS1 
!                                                                       
      CABS1 (ZDUM) = DABS (DBLE (ZDUM) ) + DABS (DIMAG (ZDUM) ) 
!                                  First executable statement           
      CALL E1PSH ('DG4CCG ') 
!                                                                       
      EPS = DMACH (4) 
      INFER = 0 
      NN = N 
!                                  Compute the machine precision times  
!                                  the norm of A and B                  
      ANORM = 0.0D0 
      BNORM = 0.0D0 
      DO 10 I = 1, N 
         IF (I.EQ.1) THEN 
            ANI = DZASUM (N - I + 1, A (I, I), N) 
         ELSE 
            ANI = ANI + CABS1 (A (I, I - 1) ) + DZASUM (N - I + 1, A (I,&
            I), N)                                                      
         ENDIF 
         BNI = DZASUM (N - I + 1, B (I, I), N) 
         ANORM = DMAX1 (ANI, ANORM) 
         BNORM = DMAX1 (BNI, BNORM) 
   10 END DO 
      IF (ANORM.EQ.0.0D0) ANORM = 1.0D0 
      IF (BNORM.EQ.0.0D0) BNORM = 1.0D0 
      EPSB = EPS * BNORM 
      EPSA = EPS * ANORM 
      IF (N.EQ.1) THEN 
         ALPHA (1) = A (1, 1) 
         BETA (1) = B (1, 1) 
         GOTO 100 
      ENDIF 
   20 ITS = 0 
!                                  Check for negligible                 
!                                  subdiagonal elements                 
   30 D2 = CABS1 (A (NN, NN) ) 
      DO 40 L = NN, 2, - 1 
         SS = D2 
         D2 = CABS1 (A (L - 1, L - 1) ) 
         SS = SS + D2 
         R = SS + CABS1 (A (L, L - 1) ) 
         IF (R.EQ.SS) GOTO 50 
   40 END DO 
      L = 1 
   50 CONTINUE 
!                                                                       
      IF (L.EQ.NN) THEN 
   55    ALPHA (NN) = A (NN, NN) 
         BETA (NN) = B (NN, NN) 
         IF (NN.EQ.1) GOTO 100 
         NN = NN - 1 
         IF (NN.GT.1) THEN 
            GOTO 20 
         ELSE 
            GOTO 55 
         ENDIF 
      ENDIF 
!                                                                       
      IF (ITS.GE.30) THEN 
         IF (INFER.EQ.0) THEN 
            INFER = NN 
            CALL E1STI (1, NN) 
            IF ( (IJOB.EQ.0) .or. (IJOB.EQ.1) ) THEN 
               WRITE ( * , * ) 'DG4CCG: iteration failed, aborting!' 
               WRITE ( * , * ) 'DG4CCG: IJOB = ', IJOB 
               CALL ende (NSIGRESTART) 
            ENDIF 
         ENDIF 
         IF (CABS1 (A (NN, NN - 1) ) .GT.0.8D0 * CABS1 (ANNM1) ) GOTO   &
         9000                                                           
      ENDIF 
!                                                                       
      IF (ITS.EQ.10.OR.ITS.EQ.20) THEN 
!                                  Ad-hoc shift                         
         NUM = DCMPLX (CABS1 (ANNM1), CABS1 (A (NN - 1, NN - 2) ) ) 
         DEN = DCMPLX (1.0D0, 0.0D0) 
      ELSE 
!                                  Compute shift as eigenvalue          
!                                  of lower 2 by 2                      
         ANNM1 = A (NN, NN - 1) 
         ANM1M1 = A (NN - 1, NN - 1) 
         S = A (NN, NN) * B (NN - 1, NN - 1) - ANNM1 * B (NN - 1, NN) 
         W = ANNM1 * B (NN, NN) * (A (NN - 1, NN) * B (NN - 1, NN - 1)  &
         - B (NN - 1, NN) * ANM1M1)                                     
         Y = (ANM1M1 * B (NN, NN) - S) / 2.0D0 
         ZC = CDSQRT (Y * Y + W) 
         IF (CABS1 (ZC) .GT.0.0D0) THEN 
            TEMP1 = Y / ZC 
            D0 = DBLE (TEMP1) 
            IF (D0.LT.0.0D0) ZC = - ZC 
         ENDIF 
         DEN = (Y + ZC) * B (NN - 1, NN - 1) * B (NN, NN) 
         IF (CABS1 (DEN) .EQ.0.0D0) DEN = DCMPLX (EPSA, 0.0D0) 
         NUM = (Y + ZC) * S - W 
      ENDIF 
!                                  Check for 2 consecutive small        
!                                  subdiagonal elements                 
      IF (NN.GT.L + 1) THEN 
         D2 = CABS1 (A (NN - 1, NN - 1) ) 
         E1 = CABS1 (ANNM1) 
         D1 = CABS1 (A (NN, NN) ) 
         DO 60 M = NN - 1, L + 1, - 1 
            E0 = E1 
            E1 = CABS1 (A (M, M - 1) ) 
            D0 = D1 
            D1 = D2 
            D2 = CABS1 (A (M - 1, M - 1) ) 
            D0 = (D0 + D1 + D2) * CABS1 (A (M, M) * DEN - B (M, M)      &
            * NUM)                                                      
            E0 = E0 * E1 * CABS1 (DEN) + D0 
            IF (E0.EQ.D0) GOTO 70 
   60    END DO 
      ENDIF 
      M = L 
   70 CONTINUE 
      ITS = ITS + 1 
      W = A (M, M) * DEN - B (M, M) * NUM 
      ZC = A (M + 1, M) * DEN 
      D1 = CABS1 (ZC) 
      D2 = CABS1 (W) 
!                                  Find L and M and set A=LAM and B=LBM 
      LOR1 = L 
      NNORN = NN 
      IF (IJOB.EQ.1) THEN 
         LOR1 = 1 
         NNORN = N 
      ENDIF 
      DO 90 I = M, NN - 1 
         J = I + 1 
!                                  Find row transformations to restore  
!                                  a to upper Hessenberg form.          
!                                  Apply transformations to A and B     
         IF (I.NE.M) THEN 
            W = A (I, I - 1) 
            ZC = A (J, I - 1) 
            D1 = CABS1 (ZC) 
            D2 = CABS1 (W) 
            IF (D1.EQ.0.0D0) GOTO 30 
         ENDIF 
         IF (D2.LE.D1) THEN 
!                                  Must interchange rows                
            CALL ZSWAP (NNORN - I + 1, A (I, I), N, A (J, I), N) 
            CALL ZSWAP (NNORN - I + 1, B (I, I), N, B (J, I), N) 
            IF (I.GT.M) A (I, I - 1) = A (J, I - 1) 
            IF (D2.EQ.0.0D0) GOTO 80 
!                                  The scaling of W and EVEC is         
!                                  designed to avoid a division by      
!                                  zero when the denominator is small   
            Y = DCMPLX (DBLE (W) / D1, DIMAG (W) / D1) / DCMPLX (DBLE ( &
            ZC) / D1, DIMAG (ZC) / D1)                                  
         ELSE 
            Y = DCMPLX (DBLE (ZC) / D2, DIMAG (ZC) / D2) / DCMPLX (DBLE &
            (W) / D2, DIMAG (W) / D2)                                   
         ENDIF 
         CALL ZAXPY (NNORN - I + 1, - Y, A (I, I), N, A (J, I), N) 
         CALL ZAXPY (NNORN - I + 1, - Y, B (I, I), N, B (J, I), N) 
   80    IF (I.GT.M) A (J, I - 1) = 0.0D0 
!                                  Perform transformations from right   
!                                  to restore B to triangular form      
!                                  apply transformations to A           
         ZC = B (J, I) 
         W = B (J, J) 
         D2 = CABS1 (W) 
         D1 = CABS1 (ZC) 
         IF (D1.EQ.0.0D0) GOTO 30 
         IF (D2.LE.D1) THEN 
!                                  Must interchange columns             
            CALL ZSWAP (J - LOR1 + 1, A (LOR1, J), 1, A (LOR1, I),      &
            1)                                                          
            CALL ZSWAP (J - LOR1 + 1, B (LOR1, J), 1, B (LOR1, I),      &
            1)                                                          
            IF (I.NE.NN - 1) THEN 
               CALL ZSWAP (1, A (J + 1, J), 1, A (J + 1, I), 1) 
            ENDIF 
            IF (IJOB.EQ.1) THEN 
               CALL ZSWAP (N, EVEC (1, J), 1, EVEC (1, I), 1) 
            ENDIF 
            B (J, I) = 0.0D0 
            IF (D2.EQ.0.0D0) GOTO 90 
            ZC = DCMPLX (DBLE (W) / D1, DIMAG (W) / D1) / DCMPLX (DBLE (&
            ZC) / D1, DIMAG (ZC) / D1)                                  
         ELSE 
            ZC = DCMPLX (DBLE (ZC) / D2, DIMAG (ZC) / D2) / DCMPLX (    &
            DBLE (W) / D2, DIMAG (W) / D2)                              
         ENDIF 
         CALL ZAXPY (J - LOR1 + 1, - ZC, A (LOR1, J), 1, A (LOR1, I),   &
         1)                                                             
         CALL ZAXPY (J - LOR1 + 1, - ZC, B (LOR1, J), 1, B (LOR1, I),   &
         1)                                                             
         B (J, I) = 0.0D0 
         IF (I.LT.NN - 1) A (J + 1, I) = A (J + 1, I) - ZC * A (J + 1,  &
         J)                                                             
         IF (IJOB.EQ.1) THEN 
            CALL ZAXPY (N, - ZC, EVEC (1, J), 1, EVEC (1, I), 1) 
         ENDIF 
   90 END DO 
      GOTO 30 
!                                  Find eigenvectors using B for        
!                                  intermediate storage                 
  100 IF (IJOB.NE.1) GOTO 9000 
      DO 130 M = N, 1, - 1 
         ALFM = A (M, M) 
         BETM = B (M, M) 
         B (M, M) = DCMPLX (1.0D0, 0.0D0) 
         DO 120 L = M - 1, 1, - 1 
            SL = 0.0D0 
            DO 110 J = L + 1, M 
               SL = SL + (BETM * A (L, J) - ALFM * B (L, J) ) * B (J, M) 
  110       END DO 
            Y = BETM * A (L, L) - ALFM * B (L, L) 
            IF (CABS1 (Y) .EQ.0.0D0) Y = DCMPLX ( (EPSA + EPSB) / 2.0D0,&
            0.0D0)                                                      
            B (L, M) = - SL / Y 
  120    END DO 
  130 END DO 
!                                  Transform to original coordinate     
!                                  system                               
      DO 150 M = N, 1, - 1 
         DO 140 I = 1, N 
            EVEC (I, M) = ZDOTU (M, EVEC (I, 1), LDEVEC, B (1, M),      &
            1)                                                          
  140    END DO 
  150 END DO 
!                                  Normalize so that largest            
!                                  component = 1.0                      
      DO 170 M = N, 1, - 1 
         I = IZAMAX (N, EVEC (1, M), 1) 
         D = EVEC (I, M) 
         IF (DBLE (D) .NE.0.0D0.OR.DIMAG (D) .NE.0.0D0) THEN 
            DO 160 I = 1, N 
               EVEC (I, M) = EVEC (I, M) / D 
  160       END DO 
         ENDIF 
  170 END DO 
!                                                                       
 9000 CALL E1POP ('DG4CCG ') 
      RETURN 
      END SUBROUTINE DG4CCG                         
!-----------------------------------------------------------------------
!  IMSL Name:  G5CCG/DG5CCG (Single/Double precision version)           
!                                                                       
!  Computer:   sgruxs/DOUBLE                                            
!                                                                       
!  Revised:    February 9, 1986                                         
!                                                                       
!  Purpose:    Lexicographical sort on complex(8) eigenvalues and       
!              corresponding eigenvectors of a complex(8) general       
!              eigensystem.                                             
!                                                                       
!  Usage:      CALL G5CCG (N, NEVEC, ALPHA, BETA, EVEC, LDEVEC, VECTOR, 
!                          EVAL)                                        
!                                                                       
!  Arguments:                                                           
!     N      - Order of the matrix.  (Input)                            
!     NEVEC  - Number of eigenvalues (eigenvectors) to be sorted.       
!              (Input)                                                  
!     ALPHA  - complex(8) vector of length N.  (Input/Output)           
!     BETA   - complex(8) vector of length N.  (Input/Output)           
!              The J-th eigenvalue is ALPHA(J)/BETA(J), assuming        
!              BETA(J) is zero then the eigenvalue is to be regarded    
!              as infinite.                                             
!     EVEC   - complex(8) matrix of order N.  (Input/Output)            
!              The J-th eigenvector, corresponding to ALPH(J)/BETA(J)   
!              is stored in the J-th column.                            
!     LDEVEC - Leading dimension of EVEC exactly as specified in the    
!              dimension statement of the calling program.  (Input)     
!     VECTOR - Logical parameter specifiying if corresponding           
!              eigenvectors to eigenvalues exist.  (Input)              
!     EVAL   - complex(8) work array of length N.  (Output)             
!              Used to store the computed eigenvalues to be sorted.     
!                                                                       
!  Chapter:    MATH/LIBRARY Eigensystem Analysis                        
!                                                                       
!  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE DG5CCG (N, NEVEC, ALPHA, BETA, EVEC, LDEVEC, VECTOR,   &
      EVAL)                                                             
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N, NEVEC, LDEVEC 
      COMPLEX (8) ALPHA ( * ), BETA ( * ), EVEC (LDEVEC, * ), EVAL ( * ) 
      LOGICAL VECTOR 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, J 
      DOUBLEPRECISION AJ, AJP1, BIG, RJ, RJP1 
      COMPLEX (8) TEMP 
!                                  SPECIFICATIONS FOR INTRINSICS        
!     INTRINSIC  DIMAG,DMIN1,CDABS,DCMPLX,DBLE                          
!      INTRINSIC  DIMAG, DMIN1, CDABS, DCMPLX, DBLE                     
!      DOUBLE PRECISION DIMAG, DMIN1, CDABS, DBLE                       
!      complex(8) DCMPLX                                                
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL ZSWAP 
!                                  SPECIFICATIONS FOR FUNCTIONS         
      EXTERNAL DMACH 
      DOUBLEPRECISION DMACH 
!                                  Calculate eigenvalues                
      BIG = DMACH (2) 
      DO 10 I = 1, NEVEC 
         IF (CDABS (ALPHA (I) ) .GT.DMIN1 (CDABS (BETA (I) ), 1.0D0)    &
         * BIG) THEN                                                    
            EVAL (I) = DCMPLX (BIG, 0.0D0) 
         ELSE 
            EVAL (I) = ALPHA (I) / BETA (I) 
         ENDIF 
   10 END DO 
!                                  Sort eigenvalues and swap according  
      DO 30 I = 1, NEVEC - 1 
         DO 20 J = 1, NEVEC - I 
            RJ = DBLE (EVAL (J) ) 
            RJP1 = DBLE (EVAL (J + 1) ) 
            AJ = DIMAG (EVAL (J) ) 
            AJP1 = DIMAG (EVAL (J + 1) ) 
            IF (RJ.GT.RJP1.OR. (RJ.EQ.RJP1.AND.AJ.GT.AJP1) ) THEN 
               TEMP = EVAL (J) 
               EVAL (J) = EVAL (J + 1) 
               EVAL (J + 1) = TEMP 
               TEMP = ALPHA (J) 
               ALPHA (J) = ALPHA (J + 1) 
               ALPHA (J + 1) = TEMP 
               TEMP = BETA (J) 
               BETA (J) = BETA (J + 1) 
               BETA (J + 1) = TEMP 
               IF (VECTOR) THEN 
                  CALL ZSWAP (NEVEC, EVEC (1, J), 1, EVEC (1, J + 1),   &
                  1)                                                    
               ENDIF 
            ENDIF 
   20    END DO 
   30 END DO 
!                                                                       
      RETURN 
      END SUBROUTINE DG5CCG                         
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  IMSL Name:  E1PSH                                                    
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    March 2, 1984                                            
!                                                                       
!  Purpose:    To push a subroutine name onto the error control stack.  
!                                                                       
!  Usage:      CALL E1PSH(NAME)                                         
!                                                                       
!  Arguments:                                                           
!     NAME   - A character string of length six specifing the name of   
!              the subroutine.  (Input)                                 
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE E1PSH (NAME) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      CHARACTER NAME * ( * ) 
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      INTEGER IFINIT 
      SAVE IFINIT 
!                                  SPECIFICATIONS FOR SPECIAL CASES     
!                              SPECIFICATIONS FOR COMMON /ERCOM1/       
      INTEGER CALLVL, MAXLEV, MSGLEN, ERTYPE (51), ERCODE (51), PRINTB (&
      7), STOPTB (7), PLEN, IFERR6, IFERR7, IALLOC (51), HDRFMT (7),    &
      TRACON (7)                                                        
      COMMON / ERCOM1 / CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, PRINTB, &
      STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, TRACON              
      SAVE / ERCOM1 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM2/       
      CHARACTER MSGSAV (255), PLIST (300), RNAME (51) * 6 
      COMMON / ERCOM2 / MSGSAV, PLIST, RNAME 
      SAVE / ERCOM2 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM3/       
      DOUBLEPRECISION ERCKSM 
      COMMON / ERCOM3 / ERCKSM 
      SAVE / ERCOM3 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM4/       
      LOGICAL ISUSER (51) 
      COMMON / ERCOM4 / ISUSER 
      SAVE / ERCOM4 / 
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL E1INIT, E1MES, E1STI 
!                                  SPECIFICATIONS FOR FUNCTIONS         
      EXTERNAL I1KST 
      INTEGER I1KST 
!                                                                       
      DATA IFINIT / 0 / 
!                                  INITIALIZE ERROR TABLE IF NECESSARY  
      IF (IFINIT.EQ.0) THEN 
         CALL E1INIT 
         IFINIT = 1 
      ENDIF 
      IF (CALLVL.GE.MAXLEV) THEN 
         CALL E1STI (1, MAXLEV) 
      CALL E1MES (5, 1, 'Error condition in E1PSH.  Push would '//'cause&
     & stack level to exceed %(I1). ')                                  
         STOP 
      ELSE 
!                                  STORE ALLOCATION LEVEL               
         IALLOC (CALLVL) = I1KST (1) 
!                                  INCREMENT THE STACK POINTER BY ONE   
         CALLVL = CALLVL + 1 
!                                  PUT SUBROUTINE NAME INTO STACK       
         RNAME (CALLVL) = NAME 
!                                  SET ERROR TYPE AND ERROR CODE        
         ERTYPE (CALLVL) = 0 
         ERCODE (CALLVL) = 0 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE E1PSH                          
!-----------------------------------------------------------------------
!  IMSL Name:  E1STI                                                    
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    March 6, 1984                                            
!                                                                       
!  Purpose:    To store an integer for subsequent use within an error   
!              message.                                                 
!                                                                       
!  Usage:      CALL E1STI(II, IVALUE)                                   
!                                                                       
!  Arguments:                                                           
!     II     - Integer specifying the substitution index.  II must be   
!              between 1 and 9.  (Input)                                
!     IVALUE - The integer to be stored.  (Input)                       
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE E1STI (II, IVALUE) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER II, IVALUE 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER IBEG, IER, ILEN 
      CHARACTER ARRAY (14) 
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      INTEGER IFINIT 
      CHARACTER BLANK (1) 
      SAVE BLANK, IFINIT 
!                                  SPECIFICATIONS FOR SPECIAL CASES     
!                              SPECIFICATIONS FOR COMMON /ERCOM1/       
      INTEGER CALLVL, MAXLEV, MSGLEN, ERTYPE (51), ERCODE (51), PRINTB (&
      7), STOPTB (7), PLEN, IFERR6, IFERR7, IALLOC (51), HDRFMT (7),    &
      TRACON (7)                                                        
      COMMON / ERCOM1 / CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, PRINTB, &
      STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, TRACON              
      SAVE / ERCOM1 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM2/       
      CHARACTER MSGSAV (255), PLIST (300), RNAME (51) * 6 
      COMMON / ERCOM2 / MSGSAV, PLIST, RNAME 
      SAVE / ERCOM2 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM3/       
      DOUBLEPRECISION ERCKSM 
      COMMON / ERCOM3 / ERCKSM 
      SAVE / ERCOM3 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM4/       
      LOGICAL ISUSER (51) 
      COMMON / ERCOM4 / ISUSER 
      SAVE / ERCOM4 / 
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL C1TIC, E1INIT, E1INPL 
!                                  SPECIFICATIONS FOR FUNCTIONS         
      EXTERNAL I1ERIF 
      INTEGER I1ERIF 
!                                                                       
      DATA BLANK / ' ' / , IFINIT / 0 / 
!                                  INITIALIZE IF NECESSARY              
      IF (IFINIT.EQ.0) THEN 
         CALL E1INIT 
         IFINIT = 1 
      ENDIF 
      CALL C1TIC (IVALUE, ARRAY, 14, IER) 
      IBEG = I1ERIF (ARRAY, 14, BLANK, 1) 
      IF (II.GE.1.AND.II.LE.9.AND.IER.EQ.0) THEN 
         ILEN = 15 - IBEG 
         CALL E1INPL ('I', II, ILEN, ARRAY (IBEG) ) 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE E1STI                          
!-----------------------------------------------------------------------
!  IMSL Name:  E1MES                                                    
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    March 2, 1984                                            
!                                                                       
!  Purpose:    Set an error state for the current level in the stack.   
!              The message is printed immediately if the error type is  
!              5, 6, or 7 and the print attribute for that type is YES. 
!                                                                       
!  Usage:      CALL E1MES(IERTYP,IERCOD,MSGPKD)                         
!                                                                       
!  Arguments:                                                           
!     IERTYP - Integer specifying the error type.  (Input)              
!                IERTYP=1,  informational/note                          
!                IERTYP=2,  informational/alert                         
!                IERTYP=3,  informational/warning                       
!                IERTYP=4,  informational/fatal                         
!                IERTYP=5,  terminal                                    
!                IERTYP=6,  PROTRAN/warning                             
!                IERTYP=7,  PROTRAN/fatal                               
!     IERCOD - Integer specifying the error code.  (Input)              
!     MSGPKD - A character string containing the message.               
!              (Input)  Within the message, any of following may appear 
!                %(A1),%(A2),...,%(A9) for character arrays             
!                %(C1),%(C2),...,%(C9) for complex(8) numbers           
!                %(D1),%(D2),...,%(D9) for double precision numbers     
!                %(I1),%(I2),...,%(I9) for integer numbers              
!                %(K1),%(K2),...,%(K9) for keywords                     
!                %(L1),%(L2),...,%(L9) for literals (strings)           
!                %(R1),%(R2),...,%(R9) for real numbers                 
!                %(Z1),%(Z2),...,%(Z9) for complex(8) numbers           
!              This provides a way to insert character arrays, strings, 
!              numbers, and keywords into the message.  See remarks     
!              below.                                                   
!                                                                       
!  Remarks:                                                             
!     The number of characters in the message after the insertion of    
!     the corresponding strings, etc. should not exceed 255.  If the    
!     limit is exceeded, only the first 255 characters will be used.    
!     The appropriate strings, etc. need to have been previously stored 
!     in common via calls to E1STA, E1STD, etc.  Line breaks may be     
!     specified by inserting the two characters '%/' into the message   
!     at the desired locations.                                         
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE E1MES (IERTYP, IERCOD, MSGPKD) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER IERTYP, IERCOD 
      CHARACTER MSGPKD * ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER ERTYP2, I, IER, IPLEN, ISUB, LAST, LEN2, LOC, M, MS, NLOC,&
      NUM, PBEG                                                         
      CHARACTER MSGTMP (255) 
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      INTEGER IFINIT, NFORMS 
      CHARACTER BLNK, DBB (3), FIND (4), FORMS (9), INREF (25), LPAR,   &
      NCHECK (3), PERCNT, RPAR                                          
      SAVE BLNK, DBB, FIND, FORMS, IFINIT, INREF, LPAR, NCHECK, NFORMS, &
      PERCNT, RPAR                                                      
!                                  SPECIFICATIONS FOR SPECIAL CASES     
!                              SPECIFICATIONS FOR COMMON /ERCOM1/       
      INTEGER CALLVL, MAXLEV, MSGLEN, ERTYPE (51), ERCODE (51), PRINTB (&
      7), STOPTB (7), PLEN, IFERR6, IFERR7, IALLOC (51), HDRFMT (7),    &
      TRACON (7)                                                        
      COMMON / ERCOM1 / CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, PRINTB, &
      STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, TRACON              
      SAVE / ERCOM1 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM2/       
      CHARACTER MSGSAV (255), PLIST (300), RNAME (51) * 6 
      COMMON / ERCOM2 / MSGSAV, PLIST, RNAME 
      SAVE / ERCOM2 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM3/       
      DOUBLEPRECISION ERCKSM 
      COMMON / ERCOM3 / ERCKSM 
      SAVE / ERCOM3 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM4/       
      LOGICAL ISUSER (51) 
      COMMON / ERCOM4 / ISUSER 
      SAVE / ERCOM4 / 
!                                  SPECIFICATIONS FOR INTRINSICS        
!     INTRINSIC  LEN,MIN0                                               
!      INTRINSIC  LEN, MIN0                                             
!      INTEGER    LEN, MIN0                                             
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL C1TCI, E1INIT, E1PRT, E1UCS, M1VE, M1VECH 
!                                  SPECIFICATIONS FOR FUNCTIONS         
      EXTERNAL I1DX 
      INTEGER I1DX 
!                                                                       
      DATA FORMS / 'A', 'C', 'D', 'I', 'K', 'L', 'R', 'S', 'Z' / ,      &
      NFORMS / 9 /                                                      
      DATA PERCNT / '%' / , LPAR / '(' / , RPAR / ')' / , BLNK / ' ' / 
      DATA INREF / ' ', 'i', 'n', ' ', 'r', 'e', 'f', 'e', 'r', 'e',    &
      'n', 'c', 'e', ' ', 't', 'o', ' ', 'k', 'e', 'y', 'w', 'o', 'r',  &
      'd', ' ' /                                                        
      DATA NCHECK / 'N', '1', '*' / , DBB / '.', ' ', ' ' / 
      DATA FIND / '*', ' ', ' ', '*' / 
      DATA IFINIT / 0 / 
!                                  INITIALIZE ERROR TABLE IF NECESSARY  
      IF (IFINIT.EQ.0) THEN 
         CALL E1INIT 
         IFINIT = 1 
      ENDIF 
!                                  CHECK AND SET ERROR TYPE IF NECESSARY
      IF (IERTYP.NE. - 1) THEN 
         ERTYPE (CALLVL) = IERTYP 
      ELSEIF (IERTYP.LT. - 1.OR.IERTYP.GT.7) THEN 
         MSGLEN = 51 
      CALL M1VECH ('.  Error from E1MES.  Illegal error type'//' specifi&
     &ed. ', MSGLEN, MSGSAV, MSGLEN)                                    
         CALL E1PRT 
         STOP 
      ENDIF 
!                                                                       
      ERTYP2 = ERTYPE (CALLVL) 
!                                  SET ERROR CODE IF NECESSARY          
      IF (IERCOD.GT. - 1) ERCODE (CALLVL) = IERCOD 
      LEN2 = LEN (MSGPKD) 
!                                                                       
      IF (IERTYP.EQ.0.OR.IERCOD.EQ.0) THEN 
!                                  REMOVE THE ERROR STATE               
         MSGLEN = 0 
      ELSEIF (LEN2.EQ.0.OR. (LEN2.EQ.1.AND.MSGPKD (1:1) .EQ.BLNK) )     &
      THEN                                                              
         IF (ERTYP2.EQ.6) IFERR6 = 1 
         IF (ERTYP2.EQ.7) IFERR7 = 1 
!                                  UPDATE CHECKSUM PARAMETER ERCKSM     
         CALL E1UCS 
!                                  PRINT MESSAGE IF NECESSARY           
         IF (ERTYP2.GE.5.AND.PRINTB (ERTYP2) .EQ.1) CALL E1PRT 
      ELSE 
!                                  FILL UP MSGSAV WITH EXPANDED MESSAGE 
         LEN2 = MIN0 (LEN2, 255) 
         DO 10 I = 1, LEN2 
            MSGTMP (I) = MSGPKD (I:I) 
   10    END DO 
         MS = 0 
         M = 0 
!                                  CHECK PLIST FOR KEYWORD NAME         
         NLOC = I1DX (PLIST, PLEN, NCHECK, 3) 
         IF (NLOC.GT.0.AND.HDRFMT (ERTYP2) .EQ.3) THEN 
!                                  M1VE INREF INTO MSGSAV               
            CALL M1VE (INREF, 1, 25, 25, MSGSAV, 1, 25, 25, IER) 
!                                  GET LENGTH OF KEYWORD NAME           
            CALL C1TCI (PLIST (NLOC + 3), 3, IPLEN, IER) 
            PBEG = NLOC + 3 + IER 
!                                  M1VE KEYWORD NAME INTO MSGSAV        
            CALL M1VE (PLIST, PBEG, PBEG + IPLEN - 1, PLEN, MSGSAV, 26, &
            IPLEN + 25, 255, IER)                                       
!                                  UPDATE POINTER                       
            MS = IPLEN + 25 
         ENDIF 
!                                  INSERT DOT, BLANK, BLANK             
         CALL M1VE (DBB, 1, 3, 3, MSGSAV, MS + 1, MS + 3, 255, IER) 
         MS = MS + 3 
!                                  LOOK AT NEXT CHARACTER               
   20    M = M + 1 
         ISUB = 0 
         IF (M.GT.LEN2 - 4) THEN 
            LAST = LEN2 - M + 1 
            DO 30 I = 1, LAST 
   30       MSGSAV (MS + I) = MSGTMP (M + I - 1) 
            MSGLEN = MS + LAST 
            GOTO 40 
         ELSEIF (MSGTMP (M) .EQ.PERCNT.AND.MSGTMP (M + 1)               &
         .EQ.LPAR.AND.MSGTMP (M + 4) .EQ.RPAR) THEN                     
            CALL C1TCI (MSGTMP (M + 3), 1, NUM, IER) 
            IF (IER.EQ.0.AND.NUM.NE.0.AND.I1DX (FORMS, NFORMS, MSGTMP ( &
            M + 2), 1) .NE.0) THEN                                      
!                                  LOCATE THE ITEM IN THE PARAMETER LIST
               CALL M1VE (MSGTMP (M + 2), 1, 2, 2, FIND, 2, 3, 4, IER) 
               LOC = I1DX (PLIST, PLEN, FIND, 4) 
               IF (LOC.GT.0) THEN 
!                                  SET IPLEN = LENGTH OF STRING         
                  CALL C1TCI (PLIST (LOC + 4), 4, IPLEN, IER) 
                  PBEG = LOC + 4 + IER 
!                                  ADJUST IPLEN IF IT IS TOO BIG        
                  IPLEN = MIN0 (IPLEN, 255 - MS) 
!                                  M1VE STRING FROM PLIST INTO MSGSAV   
                  CALL M1VE (PLIST, PBEG, PBEG + IPLEN - 1, PLEN,       &
                  MSGSAV, MS + 1, MS + IPLEN, 255, IER)                 
                  IF (IER.GE.0.AND.IER.LT.IPLEN) THEN 
!                                  UPDATE POINTERS                      
                     M = M + 4 
                     MS = MS + IPLEN - IER 
!                                  BAIL OUT IF NO MORE ROOM             
                     IF (MS.GE.255) THEN 
                        MSGLEN = 255 
                        GOTO 40 
                     ENDIF 
!                                  SET FLAG TO SHOW SUBSTITION WAS MADE 
                     ISUB = 1 
                  ENDIF 
               ENDIF 
            ENDIF 
         ENDIF 
         IF (ISUB.EQ.0) THEN 
            MS = MS + 1 
            MSGSAV (MS) = MSGTMP (M) 
         ENDIF 
         GOTO 20 
   40    ERTYP2 = ERTYPE (CALLVL) 
         IF (ERTYP2.EQ.6) IFERR6 = 1 
         IF (ERTYP2.EQ.7) IFERR7 = 1 
!                                  UPDATE CHECKSUM PARAMETER ERCKSM     
         CALL E1UCS 
!                                  PRINT MESSAGE IF NECESSARY           
         IF (ERTYP2.GE.5.AND.PRINTB (ERTYP2) .EQ.1) CALL E1PRT 
      ENDIF 
!                                  CLEAR PARAMETER LIST                 
      PLEN = 1 
!                                                                       
      RETURN 
      END SUBROUTINE E1MES                          
!-----------------------------------------------------------------------
!  IMSL Name:  ZSCAL (Double precision version)                         
!                                                                       
!  Computer:   sgruxs/DOUBLE                                            
!                                                                       
!  Revised:    August 9, 1986                                           
!                                                                       
!  Purpose:    Multiply a vector by a scalar, y = ay, both double       
!              complex(8).                                              
!                                                                       
!  Usage:      CALL ZSCAL (N, ZA, ZX, INCX)                             
!                                                                       
!  Arguments:                                                           
!     N      - Length of vectors X.  (Input)                            
!     ZA     - complex(8) scalar.  (Input)                              
!     ZX     - complex(8) vector of length MAX(N*IABS(INCX),1).         
!                 (Input/Output)                                        
!              ZSCAL replaces X(I) with ZA*X(I) for I = 1,...,N.        
!              X(I) refers to a specific element of ZX.                 
!     INCX   - Displacement between elements of ZX.  (Input)            
!              X(I) is defined to be ZX(1+(I-1)*INCX). INCX must be     
!              greater than 0.                                          
!                                                                       
!  GAMS:       D1a6                                                     
!                                                                       
!  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations              
!              STAT/LIBRARY Mathematical Support                        
!                                                                       
!  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE ZSCAL (N, ZA, ZX, INCX) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N, INCX 
      COMPLEX (8) ZA, ZX ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, IX 
!                                                                       
      IF (N.GT.0) THEN 
         IF (INCX.NE.1) THEN 
!                                  CODE FOR INCREMENTS NOT EQUAL TO 1   
            IX = 1 
            IF (INCX.LT.0) IX = ( - N + 1) * INCX + 1 
            DO 10 I = 1, N 
               ZX (IX) = ZA * ZX (IX) 
               IX = IX + INCX 
   10       END DO 
         ELSE 
!                                  CODE FOR INCREMENTS EQUAL TO 1       
            DO 20 I = 1, N 
               ZX (I) = ZA * ZX (I) 
   20       END DO 
         ENDIF 
      ENDIF 
      RETURN 
      END SUBROUTINE ZSCAL                          
!-----------------------------------------------------------------------
!  IMSL Name:  E1POP                                                    
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    March 13, 1984                                           
!                                                                       
!  Purpose:    To pop a subroutine name from the error control stack.   
!                                                                       
!  Usage:      CALL E1POP(NAME)                                         
!                                                                       
!  Arguments:                                                           
!     NAME   - A character string of length six specifying the name     
!              of the subroutine.  (Input)                              
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE E1POP (NAME) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      CHARACTER NAME * ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER IERTYP, IR 
!                                  SPECIFICATIONS FOR SPECIAL CASES     
!                              SPECIFICATIONS FOR COMMON /ERCOM1/       
      INTEGER CALLVL, MAXLEV, MSGLEN, ERTYPE (51), ERCODE (51), PRINTB (&
      7), STOPTB (7), PLEN, IFERR6, IFERR7, IALLOC (51), HDRFMT (7),    &
      TRACON (7)                                                        
      COMMON / ERCOM1 / CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, PRINTB, &
      STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, TRACON              
      SAVE / ERCOM1 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM2/       
      CHARACTER MSGSAV (255), PLIST (300), RNAME (51) * 6 
      COMMON / ERCOM2 / MSGSAV, PLIST, RNAME 
      SAVE / ERCOM2 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM3/       
      DOUBLEPRECISION ERCKSM 
      COMMON / ERCOM3 / ERCKSM 
      SAVE / ERCOM3 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM4/       
      LOGICAL ISUSER (51) 
      COMMON / ERCOM4 / ISUSER 
      SAVE / ERCOM4 / 
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL E1MES, E1PRT, E1PSH, E1STI, E1STL, I1KRL 
!                                  SPECIFICATIONS FOR FUNCTIONS         
      EXTERNAL I1KST 
      INTEGER I1KST 
!                                                                       
      IF (CALLVL.LE.1) THEN 
         CALL E1PSH ('E1POP ') 
         CALL E1STL (1, NAME) 
      CALL E1MES (5, 1, 'Error condition in E1POP.  Cannot pop '//'from &
     &%(L1) because stack is empty.')                                   
         STOP 
      ELSEIF (NAME.NE.RNAME (CALLVL) ) THEN 
         CALL E1STL (1, NAME) 
         CALL E1STL (2, RNAME (CALLVL) ) 
      CALL E1MES (5, 2, 'Error condition in E1POP.  %(L1) does '//'not m&
     &atch the name %(L2) in the stack.')                               
         STOP 
      ELSE 
         IERTYP = ERTYPE (CALLVL) 
         IF (IERTYP.NE.0) THEN 
!                                  M1VE ERROR TYPE AND ERROR CODE TO    
!                                    PREVIOUS LEVEL FOR ERROR TYPES 2-7 
            IF (IERTYP.GE.2.AND.IERTYP.LE.7) THEN 
               ERTYPE (CALLVL - 1) = ERTYPE (CALLVL) 
               ERCODE (CALLVL - 1) = ERCODE (CALLVL) 
            ENDIF 
!                                  CHECK PRINT TABLE TO DETERMINE       
!                                    WHETHER TO PRINT STORED MESSAGE    
            IF (IERTYP.LE.4) THEN 
               IF (ISUSER (CALLVL - 1) .AND.PRINTB (IERTYP) .EQ.1) CALL &
               E1PRT                                                    
            ELSE 
               IF (PRINTB (IERTYP) .EQ.1) CALL E1PRT 
            ENDIF 
!                                  CHECK STOP TABLE AND ERROR TYPE TO   
!                                    DETERMINE WHETHER TO STOP          
            IF (IERTYP.LE.4) THEN 
               IF (ISUSER (CALLVL - 1) .AND.STOPTB (IERTYP) .EQ.1) THEN 
                  STOP 
               ENDIF 
            ELSEIF (IERTYP.EQ.5) THEN 
               IF (STOPTB (IERTYP) .EQ.1) THEN 
                  STOP 
               ENDIF 
            ELSEIF (HDRFMT (IERTYP) .EQ.1) THEN 
               IF (ISUSER (CALLVL - 1) ) THEN 
                  IF (N1RGB (0) .NE.0) THEN 
                     STOP 
                  ENDIF 
               ENDIF 
            ENDIF 
         ENDIF 
!                                  SET ERROR TYPE AND CODE              
         IF (CALLVL.LT.MAXLEV) THEN 
            ERTYPE (CALLVL + 1) = - 1 
            ERCODE (CALLVL + 1) = - 1 
         ENDIF 
!                                  SET IR = AMOUNT OF WORKSPACE         
!                                  ALLOCATED AT THIS LEVEL              
         IR = I1KST (1) - IALLOC (CALLVL - 1) 
         IF (IR.GT.0) THEN 
!                                  RELEASE WORKSPACE                    
            CALL I1KRL (IR) 
            IALLOC (CALLVL) = 0 
         ELSEIF (IR.LT.0) THEN 
            CALL E1STI (1, CALLVL) 
            CALL E1STI (2, IALLOC (CALLVL - 1) ) 
            CALL E1STI (3, I1KST (1) ) 
      CALL E1MES (5, 3, 'Error condition in E1POP. '//' The number of wo&
     &rkspace allocations at '//'level %(I1) is %(I2).  However, the tot&
     &al '//'number of workspace allocations is %(I3).')                
            STOP 
         ENDIF 
!                                  DECREASE THE STACK POINTER BY ONE    
         CALLVL = CALLVL - 1 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE E1POP                          
!-----------------------------------------------------------------------
!  IMSL Name:  E1INIT                                                   
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    March 13, 1984                                           
!                                                                       
!  Purpose:    Initialization.                                          
!                                                                       
!  Usage:      CALL E1INIT                                              
!                                                                       
!  Arguments:  None                                                     
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE E1INIT 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER L 
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      INTEGER ISINIT 
      SAVE ISINIT 
!                                  SPECIFICATIONS FOR SPECIAL CASES     
!                              SPECIFICATIONS FOR COMMON /ERCOM1/       
      INTEGER CALLVL, MAXLEV, MSGLEN, ERTYPE (51), ERCODE (51), PRINTB (&
      7), STOPTB (7), PLEN, IFERR6, IFERR7, IALLOC (51), HDRFMT (7),    &
      TRACON (7)                                                        
      COMMON / ERCOM1 / CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, PRINTB, &
      STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, TRACON              
      SAVE / ERCOM1 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM2/       
      CHARACTER MSGSAV (255), PLIST (300), RNAME (51) * 6 
      COMMON / ERCOM2 / MSGSAV, PLIST, RNAME 
      SAVE / ERCOM2 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM3/       
      DOUBLEPRECISION ERCKSM 
      COMMON / ERCOM3 / ERCKSM 
      SAVE / ERCOM3 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM4/       
      LOGICAL ISUSER (51) 
      COMMON / ERCOM4 / ISUSER 
      SAVE / ERCOM4 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM8/       
      INTEGER PROLVL, XXLINE (10), XXPLEN (10), ICALOC (10), INALOC (10) 
      COMMON / ERCOM8 / PROLVL, XXLINE, XXPLEN, ICALOC, INALOC 
      SAVE / ERCOM8 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM9/       
      CHARACTER XXPROC (10) * 31 
      COMMON / ERCOM9 / XXPROC 
      SAVE / ERCOM9 / 
!                                                                       
      DATA ISINIT / 0 / 
!                                                                       
      IF (ISINIT.EQ.0) THEN 
!                                  INITIALIZE                           
         CALLVL = 1 
         ERCODE (1) = 0 
         ERTYPE (1) = 0 
         IALLOC (1) = 0 
         ISUSER (1) = .TRUE. 
         IFERR6 = 0 
         IFERR7 = 0 
         PLEN = 1 
         MAXLEV = 50 
         DO 10 L = 2, 51 
            ERTYPE (L) = - 1 
            ERCODE (L) = - 1 
            IALLOC (L) = 0 
            ISUSER (L) = .FALSE. 
   10    END DO 
         DO 20 L = 1, 7 
            HDRFMT (L) = 1 
            TRACON (L) = 1 
   20    END DO 
         PROLVL = 1 
         DO 30 L = 1, 10 
   30    ICALOC (L) = 0 
         XXLINE (1) = 0 
         XXPLEN (1) = 1 
         XXPROC (1) = '?' 
         RNAME (1) = 'USER' 
         PRINTB (1) = 0 
         PRINTB (2) = 0 
         DO 40 L = 3, 7 
   40    PRINTB (L) = 1 
         STOPTB (1) = 0 
         STOPTB (2) = 0 
         STOPTB (3) = 0 
         STOPTB (4) = 1 
         STOPTB (5) = 1 
         STOPTB (6) = 0 
         STOPTB (7) = 1 
         ERCKSM = 0.0D0 
!                                  SET FLAG TO INDICATE THAT            
!                                    INITIALIZATION HAS OCCURRED        
         ISINIT = 1 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE E1INIT                         
!-----------------------------------------------------------------------
!  IMSL Name:  I1KST                                                    
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    August 9, 1983                                           
!                                                                       
!  Purpose:    Return control information about the workspace stack.    
!                                                                       
!  Usage:      I1KST(NFACT)                                             
!                                                                       
!  Arguments:                                                           
!     NFACT  - Integer value between 1 and 6 inclusive returns the      
!                 following information: (Input)                        
!                   NFACT = 1 - LOUT: number of current allocations     
!                               excluding permanent storage. At the     
!                               end of a run, there should be no        
!                               active allocations.                     
!                   NFACT = 2 - LNOW: current active length             
!                   NFACT = 3 - LTOTAL: total storage used thus far     
!                   NFACT = 4 - LMAX: maximum storage allowed           
!                   NFACT = 5 - LALC: total number of allocations made  
!                               by I1KGT thus far                       
!                   NFACT = 6 - LNEED: number of numeric storage units  
!                               by which the stack size must be         
!                               increased for all past allocations      
!                               to succeed                              
!     I1KST  - Integer function. (Output) Returns a workspace stack     
!              statistic according to value of NFACT.                   
!                                                                       
!  Copyright:  1983 by IMSL, Inc.  All Rights Reserved                  
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      INTEGER FUNCTION I1KST (NFACT) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER NFACT 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER ISTATS (7) 
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      LOGICAL FIRST 
      SAVE FIRST 
!                                  SPECIFICATIONS FOR SPECIAL CASES     
!                                  SPECIFICATIONS FOR COMMON /WORKSP/   
      REAL RWKSP (5000) 
      REAL RDWKSP (5000) 
      DOUBLEPRECISION DWKSP (2500) 
      COMPLEX (8) CWKSP (2500) 
      COMPLEX (8) CZWKSP (2500) 
      COMPLEX (8) ZWKSP (1250) 
      INTEGER IWKSP (5000) 
      LOGICAL LWKSP (5000) 
      EQUIVALENCE (DWKSP (1), RWKSP (1) ) 
      EQUIVALENCE (CWKSP (1), RWKSP (1) ), (ZWKSP (1), RWKSP (1) ) 
      EQUIVALENCE (IWKSP (1), RWKSP (1) ), (LWKSP (1), RWKSP (1) ) 
      EQUIVALENCE (RDWKSP (1), RWKSP (1) ), (CZWKSP (1), RWKSP (1) ) 
      COMMON / WORKSP / DWKSP 
!                                  SPECIFICATIONS FOR EQUIVALENCE       
      EQUIVALENCE (ISTATS (1), IWKSP (1) ) 
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL E1MES, IWKIN 
!                                                                       
      DATA FIRST / .TRUE. / 
!                                                                       
      IF (FIRST) THEN 
!                                  INITIALIZE WORKSPACE IF NEEDED       
         FIRST = .FALSE. 
         CALL IWKIN (0) 
      ENDIF 
!                                                                       
      IF (NFACT.LE.0.OR.NFACT.GE.7) THEN 
      CALL E1MES (5, 9, 'Error from subroutine I1KST:  Argument'//' for &
     &I1KST must be between 1 and 6 inclusive.')                        
      ELSEIF (NFACT.EQ.1) THEN 
!                                  LOUT                                 
         I1KST = ISTATS (1) 
      ELSEIF (NFACT.EQ.2) THEN 
!                                  LNOW + PERMANENT                     
         I1KST = ISTATS (2) + (ISTATS (5) - ISTATS (4) + 1) 
      ELSEIF (NFACT.EQ.3) THEN 
!                                  LUSED + PERMANENT                    
         I1KST = ISTATS (3) + (ISTATS (5) - ISTATS (4) + 1) 
      ELSEIF (NFACT.EQ.4) THEN 
!                                  LMAX                                 
         I1KST = ISTATS (5) 
      ELSEIF (NFACT.EQ.5) THEN 
!                                  LALC                                 
         I1KST = ISTATS (6) 
      ELSEIF (NFACT.EQ.6) THEN 
!                                  LNEED                                
         I1KST = ISTATS (7) 
      ENDIF 
!                                                                       
      RETURN 
      END FUNCTION I1KST                            
!-----------------------------------------------------------------------
!  IMSL Name:  C1TIC                                                    
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    March 9, 1984                                            
!                                                                       
!  Purpose:    Convert an integer to its corresponding character form.  
!              (Right justified)                                        
!                                                                       
!  Usage:      CALL C1TIC(NUM, CHRSTR, SLEN, IER)                       
!                                                                       
!  Arguments:                                                           
!     NUM    - Integer number.  (Input)                                 
!     CHRSTR - Character array that receives the result.  (Output)      
!     SLEN   - Length of the character array.  (Input)                  
!     IER    - Completion code.  (Output) Where                         
!                 IER < 0  indicates that SLEN <= 0,                    
!                 IER = 0  indicates normal completion,                 
!                 IER > 0  indicates that the character array is too    
!                       small to hold the complete number.  IER         
!                       indicates how many significant digits are       
!                       being truncated.                                
!                                                                       
!  Remarks:                                                             
!  1. The character array is filled in a right justified manner.        
!  2. Leading zeros are replaced by blanks.                             
!  3. Sign is inserted only for negative number.                        
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE C1TIC (NUM, CHRSTR, SLEN, IER) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER NUM, SLEN, IER 
      CHARACTER CHRSTR ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, J, K, L 
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      CHARACTER BLANK (1), DIGIT (10), MINUS (1) 
      SAVE BLANK, DIGIT, MINUS 
!                                  SPECIFICATIONS FOR INTRINSICS        
!     INTRINSIC  IABS                                                   
!      INTRINSIC  IABS                                                  
!      INTEGER    IABS                                                  
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL M1VE 
!                                                                       
      DATA DIGIT / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      DATA BLANK / ' ' / , MINUS / '-' / 
!                                  CHECK SLEN                           
      IF (SLEN.LE.0) THEN 
         IER = - 1 
         RETURN 
      ENDIF 
!                                  THE NUMBER IS ZERO                   
      IF (NUM.EQ.0) THEN 
         CALL M1VE (BLANK, 1, 1, 1, CHRSTR, 1, SLEN - 1, SLEN, I) 
         CHRSTR (SLEN) = DIGIT (1) 
         IER = 0 
         RETURN 
      ENDIF 
!                                  CONVERT NUMBER DIGIT BY DIGIT TO     
!                                  CHARACTER FORM                       
      J = SLEN 
      K = IABS (NUM) 
   10 IF (K.GT.0.AND.J.GE.1) THEN 
         L = K 
         K = K / 10 
         L = L - K * 10 
         CHRSTR (J) = DIGIT (L + 1) 
         J = J - 1 
         GOTO 10 
      ENDIF 
!                                                                       
   20 IF (K.EQ.0) THEN 
         IF (NUM.LT.0) THEN 
            CALL M1VE (MINUS, 1, 1, 1, CHRSTR, J, J, SLEN, I) 
            IF (I.NE.0) THEN 
               IER = 1 
               RETURN 
            ENDIF 
            J = J - 1 
         ENDIF 
         IER = 0 
         CALL M1VE (BLANK, 1, 1, 1, CHRSTR, 1, J, SLEN, I) 
         RETURN 
      ENDIF 
!                                  DETERMINE THE NUMBER OF SIGNIFICANT  
!                                  DIGITS BEING TRUNCATED               
      I = 0 
   30 IF (K.GT.0) THEN 
         K = K / 10 
         I = I + 1 
         GOTO 30 
      ENDIF 
!                                                                       
      IF (NUM.LT.0) I = I + 1 
      IER = I 
!                                                                       
      RETURN 
      END SUBROUTINE C1TIC                          
!-----------------------------------------------------------------------
!  IMSL Name:  I1ERIF                                                   
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    March 13, 1984                                           
!                                                                       
!  Purpose:    Return the position of the first element of a given      
!              character array which is not an element of another       
!              character array.                                         
!                                                                       
!  Usage:      I1ERIF(STR1, LEN1, STR2, LEN2)                           
!                                                                       
!  Arguments:                                                           
!     STR1   - Character array to be searched.  (Input)                 
!     LEN1   - Length of STR1.  (Input)                                 
!     STR2   - Character array to be searched for.  (Input)             
!     LEN2   - Length of STR2.  (Input)                                 
!     I1ERIF - Integer function.  (Output)                              
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      INTEGER FUNCTION I1ERIF (STR1, LEN1, STR2, LEN2) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER LEN1, LEN2 
      CHARACTER STR1 ( * ), STR2 ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I 
!                                  SPECIFICATIONS FOR FUNCTIONS         
      EXTERNAL I1X 
      INTEGER I1X 
!                              FIRST EXECUTABLE STATEMENT               
      IF (LEN1.LE.0.OR.LEN2.LE.0) THEN 
         I1ERIF = 1 
      ELSE 
         DO 10 I = 1, LEN1 
            IF (I1X (STR2, LEN2, STR1 (I), 1) .EQ.0) THEN 
               I1ERIF = I 
               RETURN 
            ENDIF 
   10    END DO 
         I1ERIF = 0 
      ENDIF 
!                                                                       
      RETURN 
      END FUNCTION I1ERIF                           
!-----------------------------------------------------------------------
!  IMSL Name:  E1INPL                                                   
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    March 2, 1984                                            
!                                                                       
!  Purpose:    To store a character string in the parameter list PLIST  
!              for use by the error message handler.                    
!                                                                       
!  Usage:      CALL E1INPL(FORM,NUM,SLEN,STRUP)                         
!                                                                       
!  Arguments:                                                           
!     FORM   - A character string of length one to be inserted into     
!              PLIST which specifies the form of the string.  (Input)   
!              For example, 'L' for string, 'A' for character array,    
!              'I' for integer, 'K' for keyword (PROTRAN only).  An     
!              asterisk is inserted into PLIST preceding FORM.          
!     NUM    - Integer to be inserted as a character into PLIST         
!              immediately following FORM.  (Input)  NUM must be between
!              1 and 9.                                                 
!     SLEN   - The number of characters in STRUP.  (Input)  LEN must be 
!              less than or equal to 255.  The character representation 
!              of SLEN is inserted into PLIST after NUM and an asterisk.
!     STRUP  - A character string of length LEN which is to be inserted 
!              into PLIST.  (Input)  Trailing blanks are ignored.       
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE E1INPL (FORM, NUM, SLEN, STRUP) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER NUM, SLEN 
      CHARACTER FORM, STRUP ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER IER, L, LEN2, LENCK, LOC, NLEN, NNUM 
      CHARACTER STRNCH (3) 
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      CHARACTER BLANK, PRCNT (1), TEMP (4) 
      SAVE BLANK, PRCNT, TEMP 
!                                  SPECIFICATIONS FOR SPECIAL CASES     
!                              SPECIFICATIONS FOR COMMON /ERCOM1/       
      INTEGER CALLVL, MAXLEV, MSGLEN, ERTYPE (51), ERCODE (51), PRINTB (&
      7), STOPTB (7), PLEN, IFERR6, IFERR7, IALLOC (51), HDRFMT (7),    &
      TRACON (7)                                                        
      COMMON / ERCOM1 / CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, PRINTB, &
      STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, TRACON              
      SAVE / ERCOM1 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM2/       
      CHARACTER MSGSAV (255), PLIST (300), RNAME (51) * 6 
      COMMON / ERCOM2 / MSGSAV, PLIST, RNAME 
      SAVE / ERCOM2 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM3/       
      DOUBLEPRECISION ERCKSM 
      COMMON / ERCOM3 / ERCKSM 
      SAVE / ERCOM3 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM4/       
      LOGICAL ISUSER (51) 
      COMMON / ERCOM4 / ISUSER 
      SAVE / ERCOM4 / 
!                                  SPECIFICATIONS FOR INTRINSICS        
!     INTRINSIC  IABS                                                   
!      INTRINSIC  IABS                                                  
!      INTEGER    IABS                                                  
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL C1TIC, M1VE 
!                                                                       
      DATA TEMP / '*', ' ', ' ', '*' / , PRCNT / '%' / , BLANK / ' ' / 
!                                                                       
      NNUM = IABS (NUM) 
      LENCK = PLEN + SLEN + 8 
      IF (NNUM.GE.1.AND.NNUM.LE.9.AND.LENCK.LE.300) THEN 
         TEMP (2) = FORM 
         CALL C1TIC (NNUM, TEMP (3), 1, IER) 
         LOC = PLEN + 1 
         IF (LOC.EQ.2) LOC = 1 
         CALL M1VE (TEMP, 1, 4, 4, PLIST (LOC), 1, 4, 262, IER) 
         LOC = LOC + 4 
         IF (NUM.LT.0) THEN 
            LEN2 = SLEN 
         ELSE 
            DO 10 L = 1, SLEN 
               LEN2 = SLEN - L + 1 
               IF (STRUP (LEN2) .NE.BLANK) GOTO 20 
   10       END DO 
            LEN2 = 1 
   20       CONTINUE 
         ENDIF 
         NLEN = 1 
         IF (LEN2.GE.10) NLEN = 2 
         IF (LEN2.GE.100) NLEN = 3 
         CALL C1TIC (LEN2, STRNCH, NLEN, IER) 
         CALL M1VE (STRNCH, 1, NLEN, 3, PLIST (LOC), 1, NLEN, 262, IER) 
         LOC = LOC + NLEN 
         CALL M1VE (PRCNT, 1, 1, 1, PLIST (LOC), 1, 1, 262, IER) 
         LOC = LOC + 1 
         CALL M1VE (STRUP, 1, LEN2, LEN2, PLIST (LOC), 1, LEN2, 262,    &
         IER)                                                           
         PLEN = LOC + LEN2 - 1 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE E1INPL                         
!-----------------------------------------------------------------------
!  IMSL Name:  M1VECH                                                   
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    December 31, 1984                                        
!                                                                       
!  Purpose:    Character substring assignment.                          
!                                                                       
!  Usage:      CALL M1VECH (STR1, LEN1, STR2, LEN2)                     
!                                                                       
!  Arguments:                                                           
!     STR1   - Source substring.  (Input)                               
!              The source substring is STR1(1:LEN1).                    
!     LEN1   - Length of STR1.  (Input)                                 
!     STR2   - Destination substring.  (Output)                         
!              The destination substring is STR2(1:LEN2).               
!     LEN2   - Length of STR2.  (Input)                                 
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE M1VECH (STR1, LEN1, STR2, LEN2) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER LEN1, LEN2 
      CHARACTER STR1 * ( * ), STR2 * ( * ) 
!                                                                       
      STR2 (1:LEN2) = STR1 (1:LEN1) 
!                                                                       
      RETURN 
      END SUBROUTINE M1VECH                         
!-----------------------------------------------------------------------
!  IMSL Name:  E1PRT                                                    
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    March 14, 1984                                           
!                                                                       
!  Purpose:    To print an error message.                               
!                                                                       
!  Usage:      CALL E1PRT                                               
!                                                                       
!  Arguments:  None                                                     
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE E1PRT 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER ALL, I, IBEG, IBLOC, IBLOC2, IEND, IER, IHDR, J, LERTYP,  &
      LOC, LOCM1, LOCX, MAXLOC, MAXTMP, MLOC, MOD, NCBEG, NLOC, NOUT    
      CHARACTER MSGTMP (70), STRING (10) 
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      CHARACTER ATLINE (9), BLANK (1), DBB (3), FROM (6), MSGTYP (8, 7),&
      PERSLA (2), QMARK, UNKNOW (8)                                     
!                                  SPECIFICATIONS FOR SPECIAL CASES     
!                              SPECIFICATIONS FOR COMMON /ERCOM1/       
      INTEGER CALLVL, MAXLEV, MSGLEN, ERTYPE (51), ERCODE (51), PRINTB (&
      7), STOPTB (7), PLEN, IFERR6, IFERR7, IALLOC (51), HDRFMT (7),    &
      TRACON (7)                                                        
      COMMON / ERCOM1 / CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, PRINTB, &
      STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, TRACON              
      SAVE / ERCOM1 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM2/       
      CHARACTER MSGSAV (255), PLIST (300), RNAME (51) * 6 
      COMMON / ERCOM2 / MSGSAV, PLIST, RNAME 
      SAVE / ERCOM2 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM3/       
      DOUBLEPRECISION ERCKSM 
      COMMON / ERCOM3 / ERCKSM 
      SAVE / ERCOM3 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM4/       
      LOGICAL ISUSER (51) 
      COMMON / ERCOM4 / ISUSER 
      SAVE / ERCOM4 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM8/       
      INTEGER PROLVL, XXLINE (10), XXPLEN (10), ICALOC (10), INALOC (10) 
      COMMON / ERCOM8 / PROLVL, XXLINE, XXPLEN, ICALOC, INALOC 
      SAVE / ERCOM8 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM9/       
      CHARACTER XXPROC (10) * 31 
      COMMON / ERCOM9 / XXPROC 
      SAVE / ERCOM9 / 
      SAVE ATLINE, BLANK, DBB, FROM, MSGTYP, PERSLA, QMARK, UNKNOW 
!                                  SPECIFICATIONS FOR INTRINSICS        
!     INTRINSIC  MIN0                                                   
!      INTRINSIC  MIN0                                                  
!      INTEGER    MIN0                                                  
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL C1TIC, M1VE, UMACH 
!                                  SPECIFICATIONS FOR FUNCTIONS         
      EXTERNAL I1DX, I1ERIF 
      INTEGER I1DX, I1ERIF 
!                                                                       
      DATA MSGTYP / 'N', 'O', 'T', 'E', ' ', ' ', ' ', ' ', 'A', 'L',   &
      'E', 'R', 'T', ' ', ' ', ' ', 'W', 'A', 'R', 'N', 'I', 'N', 'G',  &
      ' ', 'F', 'A', 'T', 'A', 'L', ' ', ' ', ' ', 'T', 'E', 'R', 'M',  &
      'I', 'N', 'A', 'L', 'W', 'A', 'R', 'N', 'I', 'N', 'G', ' ', 'F',  &
      'A', 'T', 'A', 'L', ' ', ' ', ' ' /                               
      DATA UNKNOW / 'U', 'N', 'K', 'N', 'O', 'W', 'N', ' ' / 
      DATA ATLINE / ' ', 'a', 't', ' ', 'l', 'i', 'n', 'e', ' ' / 
      DATA BLANK / ' ' / , FROM / ' ', 'f', 'r', 'o', 'm', ' ' / 
      DATA DBB / '.', ' ', ' ' / , PERSLA / '%', '/' / 
      DATA QMARK / '?' / 
!                                                                       
      IF (MSGLEN.LE.0) RETURN 
      CALL UMACH (2, NOUT) 
      MAXTMP = 70 
      MOD = 0 
      LERTYP = ERTYPE (CALLVL) 
      IHDR = HDRFMT (LERTYP) 
      IF (IHDR.EQ.3) THEN 
         IF (XXPROC (PROLVL) (1:1) .EQ.QMARK.AND.XXLINE (PROLVL) .EQ.0) &
         THEN                                                           
            IHDR = 1 
         ENDIF 
      ENDIF 
      IEND = 0 
      IF (IHDR.EQ.1.AND.ERTYPE (CALLVL) .LE.4) THEN 
         MSGTMP (1) = BLANK (1) 
         IEND = 1 
!                                  CONVERT ERROR CODE INTO CHAR STRING  
         CALL C1TIC (ERCODE (CALLVL), STRING, 10, IER) 
!                                  LOCATE START OF NON-BLANK CHARACTERS 
         IBEG = I1ERIF (STRING, 10, BLANK, 1) 
!                                  M1VE IT TO MSGTMP                    
         CALL M1VE (STRING, IBEG, 10, 10, MSGTMP, IEND+1, IEND+11 -     &
         IBEG, MAXTMP, IER)                                             
         IEND = IEND+11 - IBEG 
      ENDIF 
      IF (IHDR.NE.2) THEN 
         CALL M1VE (FROM, 1, 6, 6, MSGTMP, IEND+1, IEND+6, MAXTMP, IER) 
         IEND = IEND+6 
      ENDIF 
      IF (IHDR.EQ.3) THEN 
!                                  THIS IS A PROTRAN RUN TIME ERROR MSG.
!                                  RETRIEVE THE PROCEDURE NAME          
         CALL M1VE (XXPROC (PROLVL), 1, XXPLEN (PROLVL), 31, MSGTMP,    &
         IEND+1, IEND+XXPLEN (PROLVL), MAXTMP, IER)                     
         MLOC = IEND+XXPLEN (PROLVL) + 1 
         MSGTMP (MLOC) = BLANK (1) 
         IEND = IEND+I1DX (MSGTMP (IEND+1), XXPLEN (PROLVL) + 1, BLANK, &
         1) - 1                                                         
         IF (XXLINE (PROLVL) .GT.0) THEN 
!                                  INSERT ATLINE                        
            CALL M1VE (ATLINE, 1, 9, 9, MSGTMP, IEND+1, IEND+9, MAXTMP, &
            IER)                                                        
            IEND = IEND+9 
!                                  CONVERT PROTRAN GLOBAL LINE NUMBER   
            CALL C1TIC (XXLINE (PROLVL), STRING, 10, IER) 
!                                  LOCATE START OF NON-BLANK CHARACTERS 
            IBEG = I1ERIF (STRING, 10, BLANK, 1) 
!                                  M1VE GLOBAL LINE NUMBER TO MSGTMP    
            CALL M1VE (STRING, IBEG, 10, 10, MSGTMP, IEND+1, IEND+11 -  &
            IBEG, MAXTMP, IER)                                          
            IEND = IEND+11 - IBEG 
         ENDIF 
      ELSE 
!                                  THIS IS EITHER A LIBRARY ERROR MSG   
!                                  OR A PROTRAN PREPROCESSOR ERROR MSG  
         IF (IHDR.EQ.1) THEN 
!                                  THIS IS A LIBRARY ERROR MESSAGE.     
!                                  RETRIEVE ROUTINE NAME                
            CALL M1VE (RNAME (CALLVL), 1, 6, 6, MSGTMP, IEND+1, IEND+6, &
            MAXTMP, IER)                                                
            MSGTMP (IEND+7) = BLANK (1) 
            IEND = IEND+I1DX (MSGTMP (IEND+1), 7, BLANK, 1) - 1 
         ENDIF 
!                                  ADD DOT, BLANK, BLANK IF NEEDED      
         IF (I1DX (MSGSAV, 3, DBB, 3) .NE.1) THEN 
            CALL M1VE (DBB, 1, 3, 3, MSGTMP, IEND+1, IEND+3, MAXTMP,    &
            IER)                                                        
            IEND = IEND+3 
            MOD = 3 
         ENDIF 
      ENDIF 
!                                  MSGTMP AND MSGSAV NOW CONTAIN THE    
!                                   ERROR MESSAGE IN FINAL FORM.        
      NCBEG = 59 - IEND-MOD 
      ALL = 0 
      IBLOC = I1DX (MSGSAV, MSGLEN, PERSLA, 2) 
      IF (IBLOC.NE.0.AND.IBLOC.LT.NCBEG) THEN 
         LOCM1 = IBLOC - 1 
         LOC = IBLOC + 1 
      ELSEIF (MSGLEN.LE.NCBEG) THEN 
         LOCM1 = MSGLEN 
         ALL = 1 
      ELSE 
         LOC = NCBEG 
!                                  CHECK FOR APPROPRIATE PLACE TO SPLIT 
   10    CONTINUE 
         IF (MSGSAV (LOC) .NE.BLANK (1) ) THEN 
            LOC = LOC - 1 
            IF (LOC.GT.1) GOTO 10 
            LOC = NCBEG + 1 
         ENDIF 
         LOCM1 = LOC - 1 
      ENDIF 
!                                  NO BLANKS FOUND IN FIRST NCBEG CHARS 
      IF (LERTYP.GE.1.AND.LERTYP.LE.7) THEN 
         WRITE (NOUT, 99995) (MSGTYP (I, LERTYP), I = 1, 8), (MSGTMP (I)&
         , I = 1, IEND), (MSGSAV (I), I = 1, LOCM1)                     
      ELSE 
         WRITE (NOUT, 99995) (UNKNOW (I), I = 1, 8), (MSGTMP (I),       &
         I = 1, IEND), (MSGSAV (I), I = 1, LOCM1)                       
      ENDIF 
      IF (ALL.EQ.0) THEN 
!                                  PREPARE TO WRITE CONTINUATION OF     
!                                    MESSAGE                            
!                                                                       
!                                  FIND WHERE TO BREAK MESSAGE          
!                                    LOC = NUMBER OF CHARACTERS OF      
!                                          MESSAGE WRITTEN SO FAR       
   20    LOCX = LOC + 64 
         NLOC = LOC + 1 
         IBLOC2 = IBLOC 
         MAXLOC = MIN0 (MSGLEN - LOC, 64) 
         IBLOC = I1DX (MSGSAV (NLOC), MAXLOC, PERSLA, 2) 
         IF (MSGSAV (NLOC) .EQ.BLANK (1) .AND.IBLOC2.EQ.0) NLOC = NLOC +&
         1                                                              
         IF (IBLOC.GT.0) THEN 
!                                  PAGE BREAK FOUND AT IBLOC            
            LOCX = NLOC + IBLOC - 2 
            WRITE (NOUT, 99996) (MSGSAV (I), I = NLOC, LOCX) 
            LOC = NLOC + IBLOC 
            GOTO 20 
!                                  DON'T BOTHER LOOKING FOR BLANK TO    
!                                    BREAK AT IF LOCX .GE. MSGLEN       
         ELSEIF (LOCX.LT.MSGLEN) THEN 
!                                  CHECK FOR BLANK TO BREAK THE LINE    
   30       CONTINUE 
            IF (MSGSAV (LOCX) .EQ.BLANK (1) ) THEN 
!                                  BLANK FOUND AT LOCX                  
               WRITE (NOUT, 99996) (MSGSAV (I), I = NLOC, LOCX) 
               LOC = LOCX 
               GOTO 20 
            ENDIF 
            LOCX = LOCX - 1 
            IF (LOCX.GT.NLOC) GOTO 30 
            LOCX = LOC + 64 
!                                  NO BLANKS FOUND IN NEXT 64 CHARS     
            WRITE (NOUT, 99996) (MSGSAV (I), I = NLOC, LOCX) 
            LOC = LOCX 
            GOTO 20 
         ELSE 
!                                  ALL THE REST WILL FIT ON 1 LINE      
            LOCX = MSGLEN 
            WRITE (NOUT, 99996) (MSGSAV (I), I = NLOC, LOCX) 
         ENDIF 
      ENDIF 
!                                  SET LENGTH OF MSGSAV AND PLEN        
!                                    TO SHOW THAT MESSAGE HAS           
!                                    ALREADY BEEN PRINTED               
 9000 MSGLEN = 0 
      PLEN = 1 
      IF (TRACON (LERTYP) .EQ.1.AND.CALLVL.GT.2) THEN 
!                                  INITIATE TRACEBACK                   
         WRITE (NOUT, 99997) 
         DO 9005 J = CALLVL, 1, - 1 
            IF (J.GT.1) THEN 
               IF (ISUSER (J - 1) ) THEN 
                  WRITE (NOUT, 99998) RNAME (J), ERTYPE (J), ERCODE (J) 
               ELSE 
                  WRITE (NOUT, 99999) RNAME (J), ERTYPE (J), ERCODE (J) 
               ENDIF 
            ELSE 
               WRITE (NOUT, 99998) RNAME (J), ERTYPE (J), ERCODE (J) 
            ENDIF 
 9005    END DO 
      ENDIF 
!                                                                       
      RETURN 
99995 FORMAT (/, ' *** ', 8A1, ' ERROR', 59A1) 
99996 FORMAT (' *** ', 9X, 64A1) 
99997 FORMAT (14X, 'Here is a traceback of subprogram calls',           &
     &       ' in reverse order:', /, 14X, '      Routine    Error ',   &
     &       'type    Error code', /, 14X, '      -------    ',         &
     &       '----------    ----------')                                
99998 FORMAT (20X, A6, 5X, I6, 8X, I6) 
99999 FORMAT (20X, A6, 5X, I6, 8X, I6, 4X, '(Called internally)') 
      END SUBROUTINE E1PRT                          
!-----------------------------------------------------------------------
!  IMSL Name:  E1UCS                                                    
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    March 8, 1984                                            
!                                                                       
!  Purpose:    To update the checksum number for error messages.        
!                                                                       
!  Usage:      CALL E1UCS                                               
!                                                                       
!  Arguments:  None                                                     
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE E1UCS 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, IBEG, IBEG2, IEND, ILOC, IPOS, JLOC, NCODE, NLEN 
      DOUBLEPRECISION DNUM 
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      DOUBLEPRECISION DMAX 
      CHARACTER BLANK (1), COMMA (1), EQUAL (1), LPAR (1) 
      SAVE BLANK, COMMA, DMAX, EQUAL, LPAR 
!                                  SPECIFICATIONS FOR SPECIAL CASES     
!                              SPECIFICATIONS FOR COMMON /ERCOM1/       
      INTEGER CALLVL, MAXLEV, MSGLEN, ERTYPE (51), ERCODE (51), PRINTB (&
      7), STOPTB (7), PLEN, IFERR6, IFERR7, IALLOC (51), HDRFMT (7),    &
      TRACON (7)                                                        
      COMMON / ERCOM1 / CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, PRINTB, &
      STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, TRACON              
      SAVE / ERCOM1 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM2/       
      CHARACTER MSGSAV (255), PLIST (300), RNAME (51) * 6 
      COMMON / ERCOM2 / MSGSAV, PLIST, RNAME 
      SAVE / ERCOM2 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM3/       
      DOUBLEPRECISION ERCKSM 
      COMMON / ERCOM3 / ERCKSM 
      SAVE / ERCOM3 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM4/       
      LOGICAL ISUSER (51) 
      COMMON / ERCOM4 / ISUSER 
      SAVE / ERCOM4 / 
!                                  SPECIFICATIONS FOR INTRINSICS        
!     INTRINSIC  DMOD                                                   
!      INTRINSIC  DMOD                                                  
!      DOUBLE PRECISION DMOD                                            
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL S1ANUM 
!                                  SPECIFICATIONS FOR FUNCTIONS         
      EXTERNAL ICASE, I1X 
      INTEGER ICASE, I1X 
!                                                                       
      DATA BLANK (1) / ' ' / , COMMA (1) / ',' / , LPAR (1) / '(' / 
      DATA EQUAL (1) / '=' / , DMAX / 1.0D+9 / 
!                                                                       
      IF (MSGLEN.GT.1) THEN 
         IPOS = 0 
         IBEG2 = 1 
   10    IBEG = IBEG2 
         IEND = MSGLEN 
!                                  LOOK FOR BLANK, COMMA, LEFT PAREN.,  
!                                  OR EQUAL SIGN                        
         ILOC = I1X (MSGSAV (IBEG), IEND-IBEG + 1, BLANK, 1) 
         JLOC = I1X (MSGSAV (IBEG), IEND-IBEG + 1, COMMA, 1) 
         IF (ILOC.EQ.0.OR. (JLOC.GT.0.AND.JLOC.LT.ILOC) ) ILOC = JLOC 
         JLOC = I1X (MSGSAV (IBEG), IEND-IBEG + 1, LPAR, 1) 
         IF (ILOC.EQ.0.OR. (JLOC.GT.0.AND.JLOC.LT.ILOC) ) ILOC = JLOC 
         JLOC = I1X (MSGSAV (IBEG), IEND-IBEG + 1, EQUAL, 1) 
         IF (ILOC.EQ.0.OR. (JLOC.GT.0.AND.JLOC.LT.ILOC) ) ILOC = JLOC 
         IF (ILOC.GE.1) THEN 
            CALL S1ANUM (MSGSAV (IBEG + ILOC), IEND-IBEG - ILOC + 1,    &
            NCODE, NLEN)                                                
            IF (NCODE.EQ.2.OR.NCODE.EQ.3) THEN 
!                                  FLOATING POINT NUMBER FOUND.         
!                                  SET POINTERS TO SKIP OVER IT         
               IBEG2 = IBEG + ILOC + NLEN 
               IF (IBEG2.LE.MSGLEN) THEN 
                  CALL S1ANUM (MSGSAV (IBEG2), IEND-IBEG2 + 1, NCODE,   &
                  NLEN)                                                 
      IF ( (MSGSAV (IBEG2) .EQ.'+'.OR.MSGSAV (IBEG2) .EQ.'-') .AND. (NCO&
     &DE.EQ.1.OR.NCODE.EQ.2) ) THEN                                     
!                                  INTEGER IMMEDIATELY FOLLOWS A REAL AS
!                                  WITH SOME CDC NOS. LIKE 1.2345678+123
!                                  SET POINTERS TO SKIP OVER IT         
                     IF (NCODE.EQ.2.AND.MSGSAV (IBEG2 + NLEN - 1)       &
                     .EQ.'.') THEN                                      
!                                  DO NOT SKIP AN END-OF-SENTENCE PERIOD
                        IBEG2 = IBEG2 + NLEN - 1 
                     ELSE 
                        IBEG2 = IBEG2 + NLEN 
                     ENDIF 
                  ENDIF 
               ENDIF 
            ELSE 
               IBEG2 = IBEG + ILOC 
            ENDIF 
            IEND = IBEG + ILOC - 1 
         ENDIF 
!                                  UPDATE CKSUM USING PART OF MESSAGE   
         DO 20 I = IBEG, IEND 
            IPOS = IPOS + 1 
            DNUM = ICASE (MSGSAV (I) ) 
            ERCKSM = DMOD (ERCKSM + DNUM * IPOS, DMAX) 
   20    END DO 
!                                  GO BACK FOR MORE IF NEEDED           
         IF (IEND.LT.MSGLEN.AND.IBEG2.LT.MSGLEN) GOTO 10 
!                                  UPDATE CKSUM USING ERROR TYPE        
         DNUM = ERTYPE (CALLVL) 
         ERCKSM = DMOD (ERCKSM + DNUM * (IPOS + 1), DMAX) 
!                                  UPDATE CKSUM USING ERROR CODE        
         DNUM = ERCODE (CALLVL) 
         ERCKSM = DMOD (ERCKSM + DNUM * (IPOS + 2), DMAX) 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE E1UCS                          
!-----------------------------------------------------------------------
!  IMSL Name:  I1DX (Single precision version)                          
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    September 9, 1985                                        
!                                                                       
!  Purpose:    Determine the array subscript indicating the starting    
!              element at which a key character sequence begins.        
!              (Case-insensitive version)                               
!                                                                       
!  Usage:      I1DX(CHRSTR, I1LEN, KEY, KLEN)                           
!                                                                       
!  Arguments:                                                           
!     CHRSTR - Character array to be searched.  (Input)                 
!     I1LEN  - Length of CHRSTR.  (Input)                               
!     KEY    - Character array that contains the key sequence.  (Input) 
!     KLEN   - Length of KEY.  (Input)                                  
!     I1DX   - Integer function.  (Output)                              
!                                                                       
!  Remarks:                                                             
!  1. Returns zero when there is no match.                              
!                                                                       
!  2. Returns zero if KLEN is longer than ISLEN.                        
!                                                                       
!  3. Returns zero when any of the character arrays has a negative or   
!     zero length.                                                      
!                                                                       
!  GAMS:       N5c                                                      
!                                                                       
!  Chapter:    MATH/LIBRARY Utilities                                   
!                                                                       
!  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      INTEGER FUNCTION I1DX (CHRSTR, I1LEN, KEY, KLEN) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER I1LEN, KLEN 
      CHARACTER CHRSTR ( * ), KEY ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, II, J 
!                                  SPECIFICATIONS FOR FUNCTIONS         
      EXTERNAL ICASE, I1CSTR 
      INTEGER ICASE, I1CSTR 
!                                                                       
      I1DX = 0 
      IF (KLEN.LE.0.OR.I1LEN.LE.0) GOTO 9000 
      IF (KLEN.GT.I1LEN) GOTO 9000 
!                                                                       
      I = 1 
      II = I1LEN - KLEN + 1 
   10 IF (I.LE.II) THEN 
         IF (ICASE (CHRSTR (I) ) .EQ.ICASE (KEY (1) ) ) THEN 
            IF (KLEN.NE.1) THEN 
               J = KLEN - 1 
               IF (I1CSTR (CHRSTR (I + 1), J, KEY (2), J) .EQ.0) THEN 
                  I1DX = I 
                  GOTO 9000 
               ENDIF 
            ELSE 
               I1DX = I 
               GOTO 9000 
            ENDIF 
         ENDIF 
         I = I + 1 
         GOTO 10 
      ENDIF 
!                                                                       
 9000 RETURN 
      END FUNCTION I1DX                             
!-----------------------------------------------------------------------
!  IMSL Name:  M1VE                                                     
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    March 5, 1984                                            
!                                                                       
!  Purpose:    Move a subset of one character array to another.         
!                                                                       
!  Usage:      CALL M1VE(INSTR, INBEG, INEND, INLEN, OUTSTR, OUTBEG,    
!                         OUTEND, OUTLEN, IER)                          
!                                                                       
!  Arguments:                                                           
!     INSTR  - Source character array.  (Input)                         
!     INBEG  - First element of INSTR to be moved.  (Input)             
!     INEND  - Last element of INSTR to be moved.  (Input)              
!              The source subset is INSTR(INBEG),...,INSTR(INEND).      
!     INLEN  - Length of INSTR.  (Input)                                
!     OUTSTR - Destination character array.  (Output)                   
!     IUTBEG - First element of OUTSTR destination.  (Input)            
!     IUTEND - Last element of OUTSTR  destination.  (Input)            
!              The destination subset is OUTSRT(IUTBEG),...,            
!              OUTSTR(IUTEND).                                          
!     IUTLEN - Length of OUTSTR.  (Input)                               
!     IER    - Completion code.  (Output)                               
!              IER = -2  indicates that the input parameters, INBEG,    
!                        INEND, INLEN, IUTBEG, IUTEND are not           
!                        consistent.  One of the conditions             
!                        INBEG.GT.0, INEND.GE.INBEG, INLEN.GE.INEND,    
!                        IUTBEG.GT.0, or IUTEND.GE.IUTBEG is not        
!                        satisfied.                                     
!              IER = -1  indicates that the length of OUTSTR is         
!                        insufficient to hold the subset of INSTR.      
!                        That is, IUTLEN is less than IUTEND.           
!              IER =  0  indicates normal completion                    
!              IER >  0  indicates that the specified subset of OUTSTR, 
!                        OUTSTR(IUTBEG),...,OUTSTR(IUTEND) is not long  
!                        enough to hold the subset INSTR(INBEG),...,    
!                        INSTR(INEND) of INSTR.  IER is set to the      
!                        number of characters that were not moved.      
!                                                                       
!  Remarks:                                                             
!  1. If the subset of OUTSTR is longer than the subset of INSTR,       
!     trailing blanks are moved to OUTSTR.                              
!  2. If the subset of INSTR is longer than the subset of OUTSTR,       
!     the shorter subset is moved to OUTSTR and IER is set to the number
!     of characters that were not moved to OUTSTR.                      
!  3. If the length of OUTSTR is insufficient to hold the subset,       
!     IER is set to -2 and nothing is moved.                            
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE M1VE (INSTR, INBEG, INEND, INLEN, OUTSTR, IUTBEG,      &
      IUTEND, IUTLEN, IER)                                              
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER INBEG, INEND, INLEN, IUTBEG, IUTEND, IUTLEN, IER 
      CHARACTER INSTR ( * ), OUTSTR ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER IUTLAS, KI, KO 
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      CHARACTER BLANK 
      SAVE BLANK 
!                                  SPECIFICATIONS FOR INTRINSICS        
!     INTRINSIC  MIN0                                                   
!      INTRINSIC  MIN0                                                  
!      INTEGER    MIN0                                                  
!                                                                       
      DATA BLANK / ' ' / 
!                                  CHECK INBEG, INEND, INLEN, IUTBEG,   
!                                  AND IUTEND                           
!                                                                       
      IF (INBEG.LE.0.OR.INEND.LT.INBEG.OR.INLEN.LT.INEND.OR.IUTBEG.LE.0.&
     &OR.IUTEND.LT.IUTBEG) THEN                                         
         IER = - 2 
         RETURN 
      ELSEIF (IUTLEN.LT.IUTEND) THEN 
         IER = - 1 
         RETURN 
      ENDIF 
!                                  DETERMINE LAST CHARACTER TO M1VE     
      IUTLAS = IUTBEG + MIN0 (INEND-INBEG, IUTEND-IUTBEG) 
!                                  M1VE CHARACTERS                      
      KI = INBEG 
      DO 10 KO = IUTBEG, IUTLAS 
         OUTSTR (KO) = INSTR (KI) 
         KI = KI + 1 
   10 END DO 
!                                   SET IER TO NUMBER OF CHARACTERS THAT
!                                   WHERE NOT MOVED                     
      IER = KI - INEND-1 
!                                   APPEND BLANKS IF NECESSARY          
      DO 20 KO = IUTLAS + 1, IUTEND 
         OUTSTR (KO) = BLANK 
   20 END DO 
!                                                                       
      RETURN 
      END SUBROUTINE M1VE                           
!-----------------------------------------------------------------------
!  IMSL Name:  C1TCI                                                    
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    August 13, 1984                                          
!                                                                       
!  Purpose:    Convert character string into corresponding integer      
!              form.                                                    
!                                                                       
!  Usage:      CALL C1TCI (CHRSTR, SLEN, NUM, IER)                      
!                                                                       
!  Arguments:                                                           
!   CHRSTR  - Character array that contains the number description.     
!             (Input)                                                   
!   SLEN    - Length of the character array.  (Input)                   
!   NUM     - The answer.  (Output)                                     
!   IER     - Completion code.  (Output)  Where                         
!                IER =-2  indicates that the number is too large to     
!                         be converted;                                 
!                IER =-1  indicates that SLEN <= 0;                     
!                IER = 0  indicates normal completion;                  
!                IER > 0  indicates that the input string contains a    
!                         nonnumeric character.  IER is the index of    
!                         the first nonnumeric character in CHRSTR.     
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE C1TCI (CHRSTR, SLEN, NUM, IER) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER SLEN, NUM, IER 
      CHARACTER CHRSTR ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER COUNT, I, IMACH5, J, N, S, SIGN 
      CHARACTER ZERO 
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      CHARACTER BLANK, DIGIT * 10, MINUS, PLUS 
      SAVE BLANK, DIGIT, MINUS, PLUS 
!                                  SPECIFICATIONS FOR EQUIVALENCE       
      EQUIVALENCE (DIGIT, ZERO) 
!                                  SPECIFICATIONS FOR INTRINSICS        
!     INTRINSIC  INDEX                                                  
!      INTRINSIC  INDEX                                                 
!     INTEGER    INDEX                                                  
!                                  SPECIFICATIONS FOR FUNCTIONS         
      EXTERNAL IMACH 
      INTEGER IMACH 
!                                                                       
      DATA DIGIT / '0123456789' / 
      DATA BLANK / ' ' / , MINUS / '-' / , PLUS / '+' / 
!                                                                       
!                                  CHECK SLEN                           
      NUM = 0 
      IF (SLEN.LE.0) THEN 
         IER = - 1 
         GOTO 50 
      ENDIF 
!                                  HANDLE LEADING BLANKS                
      SIGN = 1 
      I = 1 
   10 IF (I.LE.SLEN) THEN 
         IF (CHRSTR (I) .EQ.BLANK) THEN 
            I = I + 1 
            GOTO 10 
         ENDIF 
      ELSE 
         IER = 1 
         GOTO 50 
      ENDIF 
!                                  CHECK FOR SIGN, IF ANY               
      S = I 
      IF (CHRSTR (I) .EQ.MINUS) THEN 
         SIGN = - 1 
         I = I + 1 
      ELSEIF (CHRSTR (I) .EQ.PLUS) THEN 
         I = I + 1 
      ENDIF 
   20 IF (I.LE.SLEN) THEN 
         IF (CHRSTR (I) .EQ.BLANK) THEN 
            I = I + 1 
            GOTO 20 
         ENDIF 
      ELSE 
         IER = S 
         GOTO 50 
      ENDIF 
!                                  SKIP LEADING ZERO                    
      J = I 
   30 IF (I.LE.SLEN) THEN 
         IF (CHRSTR (I) .EQ.ZERO) THEN 
            I = I + 1 
            GOTO 30 
         ENDIF 
      ELSE 
         IER = 0 
         GOTO 50 
      ENDIF 
!                                  CHECK FIRST NONBLANK CHARACTER       
      COUNT = 0 
!                                  CHECK NUMERIC CHARACTERS             
      IMACH5 = IMACH (5) 
   40 N = INDEX (DIGIT, CHRSTR (I) ) 
      IF (N.NE.0) THEN 
         COUNT = COUNT + 1 
         IF (NUM.GT. ( (IMACH5 - N) + 1) / 10) THEN 
            IER = - 2 
            GOTO 50 
         ELSE 
            NUM = NUM * 10 - 1 + N 
            I = I + 1 
            IF (I.LE.SLEN) GOTO 40 
         ENDIF 
      ENDIF 
!                                                                       
      IF (COUNT.EQ.0) THEN 
         IF (I.GT.J) THEN 
            IER = I 
         ELSE 
            IER = S 
         ENDIF 
      ELSEIF (I.GT.SLEN) THEN 
         NUM = SIGN * NUM 
         IER = 0 
      ELSE 
         NUM = SIGN * NUM 
         IER = I 
      ENDIF 
!                                                                       
   50 CONTINUE 
      RETURN 
      END SUBROUTINE C1TCI                          
!-----------------------------------------------------------------------
!  IMSL Name:  E1STL                                                    
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    November 8, 1985                                         
!                                                                       
!  Purpose:    To store a string for subsequent use within an error     
!              message.                                                 
!                                                                       
!  Usage:      CALL E1STL(IL,STRING)                                    
!                                                                       
!  Arguments:                                                           
!     IL     - Integer specifying the substitution index.  IL must be   
!              between 1 and 9.  (Input)                                
!     STRING - A character string.  (Input)                             
!                                                                       
!  Copyright:  1985 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE E1STL (IL, STRING) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER IL 
      CHARACTER STRING * ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, LEN2 
      CHARACTER STRGUP (255) 
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      INTEGER IFINIT 
      SAVE IFINIT 
!                                  SPECIFICATIONS FOR SPECIAL CASES     
!                              SPECIFICATIONS FOR COMMON /ERCOM1/       
      INTEGER CALLVL, MAXLEV, MSGLEN, ERTYPE (51), ERCODE (51), PRINTB (&
      7), STOPTB (7), PLEN, IFERR6, IFERR7, IALLOC (51), HDRFMT (7),    &
      TRACON (7)                                                        
      COMMON / ERCOM1 / CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, PRINTB, &
      STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, TRACON              
      SAVE / ERCOM1 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM2/       
      CHARACTER MSGSAV (255), PLIST (300), RNAME (51) * 6 
      COMMON / ERCOM2 / MSGSAV, PLIST, RNAME 
      SAVE / ERCOM2 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM3/       
      DOUBLEPRECISION ERCKSM 
      COMMON / ERCOM3 / ERCKSM 
      SAVE / ERCOM3 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM4/       
      LOGICAL ISUSER (51) 
      COMMON / ERCOM4 / ISUSER 
      SAVE / ERCOM4 / 
!                                  SPECIFICATIONS FOR INTRINSICS        
!     INTRINSIC  IABS,LEN,MIN0                                          
!      INTRINSIC  IABS, LEN, MIN0                                       
!      INTEGER    IABS, LEN, MIN0                                       
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL E1INIT, E1INPL 
!                                                                       
      DATA IFINIT / 0 / 
!                                  INITIALIZE IF NECESSARY              
      IF (IFINIT.EQ.0) THEN 
         CALL E1INIT 
         IFINIT = 1 
      ENDIF 
      LEN2 = LEN (STRING) 
      LEN2 = MIN0 (LEN2, 255) 
      DO 10 I = 1, LEN2 
         STRGUP (I) = STRING (I:I) 
   10 END DO 
      IF (IABS (IL) .GE.1.AND.IABS (IL) .LE.9) THEN 
         CALL E1INPL ('L', IL, LEN2, STRGUP) 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE E1STL                          
!-----------------------------------------------------------------------
!  IMSL Name:  N1RGB                                                    
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    March 2, 1984                                            
!                                                                       
!  Purpose:    Return a positive number as a flag to indicated that a   
!              stop should occur due to one or more global errors.      
!                                                                       
!  Usage:      N1RGB(IDUMMY)                                            
!                                                                       
!  Arguments:                                                           
!     IDUMMY - Integer scalar dummy argument.                           
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      INTEGER FUNCTION N1RGB (IDUMMY) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER IDUMMY 
!                                  SPECIFICATIONS FOR SPECIAL CASES     
!                              SPECIFICATIONS FOR COMMON /ERCOM1/       
      INTEGER CALLVL, MAXLEV, MSGLEN, ERTYPE (51), ERCODE (51), PRINTB (&
      7), STOPTB (7), PLEN, IFERR6, IFERR7, IALLOC (51), HDRFMT (7),    &
      TRACON (7)                                                        
      COMMON / ERCOM1 / CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, PRINTB, &
      STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, TRACON              
      SAVE / ERCOM1 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM2/       
      CHARACTER MSGSAV (255), PLIST (300), RNAME (51) * 6 
      COMMON / ERCOM2 / MSGSAV, PLIST, RNAME 
      SAVE / ERCOM2 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM3/       
      DOUBLEPRECISION ERCKSM 
      COMMON / ERCOM3 / ERCKSM 
      SAVE / ERCOM3 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM4/       
      LOGICAL ISUSER (51) 
      COMMON / ERCOM4 / ISUSER 
      SAVE / ERCOM4 / 
!                                  INITIALIZE FUNCTION                  
      N1RGB = 0 
!                                  CHECK FOR GLOBAL ERROR TYPE 6        
      IF (IFERR6.GT.0) THEN 
         N1RGB = STOPTB (6) 
         IFERR6 = 0 
      ENDIF 
!                                  CHECK FOR GLOBAL ERROR TYPE 7        
      IF (IFERR7.GT.0) THEN 
         N1RGB = STOPTB (7) 
         IFERR7 = 0 
      ENDIF 
!                                                                       
      RETURN 
      END FUNCTION N1RGB                            
!-----------------------------------------------------------------------
!  IMSL Name:  I1KRL                                                    
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    August 9, 1983                                           
!                                                                       
!  Purpose:    Deallocate the last N allocations made in the workspace. 
!              stack by I1KGT                                           
!                                                                       
!  Usage:      CALL I1KRL(N)                                            
!                                                                       
!  Arguments:                                                           
!     N      - Number of allocations to be released top down (Input)    
!                                                                       
!  Copyright:  1983 by IMSL, Inc.  All Rights Reserved                  
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE I1KRL (N) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, IN, LALC, LBND, LBOOK, LMAX, LNEED, LNOW, LOUT, LUSED, &
      NDX, NEXT                                                         
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      LOGICAL FIRST 
      SAVE FIRST 
!                                  SPECIFICATIONS FOR SPECIAL CASES     
!                                  SPECIFICATIONS FOR COMMON /WORKSP/   
      REAL RWKSP (5000) 
      REAL RDWKSP (5000) 
      DOUBLEPRECISION DWKSP (2500) 
      COMPLEX (8) CWKSP (2500) 
      COMPLEX (8) CZWKSP (2500) 
      COMPLEX (8) ZWKSP (1250) 
      INTEGER IWKSP (5000) 
      LOGICAL LWKSP (5000) 
      EQUIVALENCE (DWKSP (1), RWKSP (1) ) 
      EQUIVALENCE (CWKSP (1), RWKSP (1) ), (ZWKSP (1), RWKSP (1) ) 
      EQUIVALENCE (IWKSP (1), RWKSP (1) ), (LWKSP (1), RWKSP (1) ) 
      EQUIVALENCE (RDWKSP (1), RWKSP (1) ), (CZWKSP (1), RWKSP (1) ) 
      COMMON / WORKSP / DWKSP 
!                                  SPECIFICATIONS FOR EQUIVALENCE       
      EQUIVALENCE (LOUT, IWKSP (1) ) 
      EQUIVALENCE (LNOW, IWKSP (2) ) 
      EQUIVALENCE (LUSED, IWKSP (3) ) 
      EQUIVALENCE (LBND, IWKSP (4) ) 
      EQUIVALENCE (LMAX, IWKSP (5) ) 
      EQUIVALENCE (LALC, IWKSP (6) ) 
      EQUIVALENCE (LNEED, IWKSP (7) ) 
      EQUIVALENCE (LBOOK, IWKSP (8) ) 
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL E1MES, E1STI, IWKIN 
!                                                                       
      DATA FIRST / .TRUE. / 
!                                                                       
      IF (FIRST) THEN 
!                                  INITIALIZE WORKSPACE IF NEEDED       
         FIRST = .FALSE. 
         CALL IWKIN (0) 
      ENDIF 
!                                  CALLING I1KRL(0) WILL CONFIRM        
!                                  INTEGRITY OF SYSTEM AND RETURN       
      IF (N.LT.0) THEN 
      CALL E1MES (5, 10, 'Error from subroutine I1KRL:  Attempt'//' to r&
     &elease a negative number of workspace'//' allocations. ')         
         GOTO 9000 
      ENDIF 
!                                  BOOKKEEPING OVERWRITTEN              
      IF (LNOW.LT.LBOOK.OR.LNOW.GT.LUSED.OR.LUSED.GT.LMAX.OR.LNOW.GE.LBN&
     &D.OR.LOUT.GT.LALC) THEN                                           
      CALL E1MES (5, 11, 'Error from subroutine I1KRL:  One or '//'more &
     &of the first eight bookkeeping locations '//'in IWKSP have been ov&
     &erwritten.  ')                                                    
         GOTO 9000 
      ENDIF 
!                                  CHECK ALL THE POINTERS IN THE        
!                                  PERMANENT STORAGE AREA.  THEY MUST   
!                                  BE MONOTONE INCREASING AND LESS THAN 
!                                  OR EQUAL TO LMAX, AND THE INDEX OF   
!                                  THE LAST POINTER MUST BE LMAX+1.     
      NDX = LBND 
      IF (NDX.NE.LMAX + 1) THEN 
         DO 10 I = 1, LALC 
            NEXT = IWKSP (NDX) 
            IF (NEXT.EQ.LMAX + 1) GOTO 20 
!                                                                       
            IF (NEXT.LE.NDX.OR.NEXT.GT.LMAX) THEN 
      CALL E1MES (5, 12, 'Error from subroutine I1KRL:  '//'A pointer in&
     & permanent storage has been '//' overwritten. ')                  
               GOTO 9000 
            ENDIF 
            NDX = NEXT 
   10    END DO 
      CALL E1MES (5, 13, 'Error from subroutine I1KRL:  A '//'pointer in&
     & permanent storage has been '//'overwritten. ')                   
         GOTO 9000 
      ENDIF 
   20 IF (N.GT.0) THEN 
         DO 30 IN = 1, N 
            IF (LNOW.LE.LBOOK) THEN 
      CALL E1MES (5, 14, 'Error from subroutine I1KRL:  '//'Attempt to r&
     &elease a nonexistant '//'workspace  allocation. ')                
               GOTO 9000 
            ELSEIF (IWKSP (LNOW) .LT.LBOOK.OR.IWKSP (LNOW) .GE.LNOW - 1)&
            THEN                                                        
!                                  CHECK TO MAKE SURE THE BACK POINTERS 
!                                  ARE MONOTONE.                        
               CALL E1STI (1, LNOW) 
      CALL E1MES (5, 15, 'Error from subroutine I1KRL:  '//'The pointer &
     &at IWKSP(%(I1)) has been '//'overwritten.  ')                     
               GOTO 9000 
            ELSE 
               LOUT = LOUT - 1 
               LNOW = IWKSP (LNOW) 
            ENDIF 
   30    END DO 
      ENDIF 
!                                                                       
 9000 RETURN 
      END SUBROUTINE I1KRL                          
!-----------------------------------------------------------------------
!  IMSL Name:  N1RTY                                                    
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    March 6, 1984                                            
!                                                                       
!  Purpose:    Retrieve an error type.                                  
!                                                                       
!  Usage:      N1RTY(IOPT)                                              
!                                                                       
!  Arguments:                                                           
!     IOPT   - Integer specifying the level.  (Input)                   
!              If IOPT=0 the error type for the current level is        
!              returned.  If IOPT=1 the error type for the most         
!              recently called routine (last pop) is returned.          
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      INTEGER FUNCTION N1RTY (IOPT) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER IOPT 
!                                  SPECIFICATIONS FOR SPECIAL CASES     
!                              SPECIFICATIONS FOR COMMON /ERCOM1/       
      INTEGER CALLVL, MAXLEV, MSGLEN, ERTYPE (51), ERCODE (51), PRINTB (&
      7), STOPTB (7), PLEN, IFERR6, IFERR7, IALLOC (51), HDRFMT (7),    &
      TRACON (7)                                                        
      COMMON / ERCOM1 / CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, PRINTB, &
      STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, TRACON              
      SAVE / ERCOM1 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM2/       
      CHARACTER MSGSAV (255), PLIST (300), RNAME (51) * 6 
      COMMON / ERCOM2 / MSGSAV, PLIST, RNAME 
      SAVE / ERCOM2 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM3/       
      DOUBLEPRECISION ERCKSM 
      COMMON / ERCOM3 / ERCKSM 
      SAVE / ERCOM3 / 
!                              SPECIFICATIONS FOR COMMON /ERCOM4/       
      LOGICAL ISUSER (51) 
      COMMON / ERCOM4 / ISUSER 
      SAVE / ERCOM4 / 
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL E1PRT, M1VECH 
!                                                                       
      IF (IOPT.NE.0.AND.IOPT.NE.1) THEN 
         ERTYPE (CALLVL) = 5 
         ERCODE (CALLVL) = 1 
         MSGLEN = 47 
      CALL M1VECH ('.  The argument passed to N1RTY must be 0 or '//'1. &
     &', MSGLEN, MSGSAV, MSGLEN)                                        
         CALL E1PRT 
         STOP 
      ELSE 
         N1RTY = ERTYPE (CALLVL + IOPT) 
      ENDIF 
!                                                                       
      RETURN 
      END FUNCTION N1RTY                            
!-----------------------------------------------------------------------
!  IMSL Name:  IZAMAX (Single precision version)                        
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    August 9, 1986                                           
!                                                                       
!  Purpose:    Find the smallest index of the component of a            
!              double-complex(8) vector having maximum magnitude.       
!                                                                       
!  Usage:      IZAMAX(N, ZX, INCX)                                      
!                                                                       
!  Arguments:                                                           
!     N      - Length of vector X.  (Input)                             
!     ZX     - complex(8) vector of length N*INCX.  (Input)             
!     INCX   - Displacement between elements of ZX.  (Input)            
!              X(I) is defined to be ZX(1+(I-1)*INCX). INCX must be     
!              greater than zero.                                       
!     IZAMAX - The smallest index I such that DCABS(X(I)) is the maximum
!              of DCABS(X(J)) for J=1 to N.  (Output)                   
!              X(I) refers to a specific element of ZX.                 
!                                                                       
!  GAMS:       D1a2                                                     
!                                                                       
!  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations              
!              STAT/LIBRARY Mathematical Support                        
!                                                                       
!  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      INTEGER FUNCTION IZAMAX (N, ZX, INCX) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N, INCX 
      COMPLEX (8) ZX ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, IX 
      DOUBLEPRECISION SMAX 
      COMPLEX (8) ZDUM 
!                                  SPECIFICATIONS FOR INTRINSICS        
!     INTRINSIC  DABS,DBLE,DIMAG                                        
!      INTRINSIC  DABS, DBLE, DIMAG                                     
!      DOUBLE PRECISION DABS, DBLE, DIMAG                               
      DOUBLEPRECISION ZABS1 
!                                                                       
      ZABS1 (ZDUM) = DABS (DBLE (ZDUM) ) + DABS (DIMAG (ZDUM) ) 
!                                                                       
      IZAMAX = 0 
      IF (N.GE.1) THEN 
         IZAMAX = 1 
         IF (N.NE.1) THEN 
            IF (INCX.NE.1) THEN 
!                                  CODE FOR INCREMENTS NOT EQUAL TO 1   
               IX = 1 
               IF (INCX.LT.0) IX = ( - N + 1) * INCX + 1 
               SMAX = ZABS1 (ZX (IX) ) 
               IX = IX + INCX 
               DO 10 I = 2, N 
                  IF (ZABS1 (ZX (IX) ) .GT.SMAX) THEN 
                     IZAMAX = I 
                     SMAX = ZABS1 (ZX (IX) ) 
                  ENDIF 
                  IX = IX + INCX 
   10          END DO 
            ELSE 
!                                  CODE FOR INCREMENTS EQUAL TO 1       
               SMAX = ZABS1 (ZX (1) ) 
               DO 20 I = 2, N 
                  IF (ZABS1 (ZX (I) ) .GT.SMAX) THEN 
                     IZAMAX = I 
                     SMAX = ZABS1 (ZX (I) ) 
                  ENDIF 
   20          END DO 
            ENDIF 
         ENDIF 
      ENDIF 
      RETURN 
      END FUNCTION IZAMAX                           
!-----------------------------------------------------------------------
!  IMSL Name:  ZCOPY (Double precision version)                         
!                                                                       
!  Computer:   sgruxs/DOUBLE                                            
!                                                                       
!  Revised:    August 9, 1986                                           
!                                                                       
!  Purpose:    Copy a vector X to a vector Y, both complex(8).          
!                                                                       
!  Usage:      CALL ZCOPY (N, ZX, INCX, ZY, INCY)                       
!                                                                       
!  Arguments:                                                           
!     N      - Length of vectors X and Y.  (Input)                      
!     ZX     - complex(8) vector of length MAX(N*IABS(INCX),1).         
!              (Input)                                                  
!     INCX   - Displacement between elements of ZX.  (Input)            
!              X(I) is defined to be                                    
!                 ZX(1+(I-1)*INCX) if INCX.GE.0  or                     
!                 ZX(1+(I-N)*INCX) if INCX.LT.0.                        
!     ZY     - complex(8) vector of length MAX(N*IABS(INCY),1).         
!              (Output)                                                 
!              ZCOPY copies X(I) to Y(I) for I = 1,...N. X(I) and Y(I)  
!              refer to specific elements of ZX and ZY.                 
!     INCY   - Displacement between elements of ZY.  (Input)            
!              Y(I) is defined to be                                    
!                 ZY(1+(I-1)*INCY) if INCY.GE.0  or                     
!                 ZY(1+(I-N)*INCY) if INCY.LT.0.                        
!                                                                       
!  GAMS:       D1a                                                      
!                                                                       
!  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations              
!              STAT/LIBRARY Mathematical Support                        
!                                                                       
!  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE ZCOPY (N, ZX, INCX, ZY, INCY) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N, INCX, INCY 
      COMPLEX (8) ZX ( * ), ZY ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, IX, IY 
!                                                                       
      IF (N.GT.0) THEN 
         IF (INCX.NE.1.OR.INCY.NE.1) THEN 
!                                  CODE FOR UNEQUAL INCREMENTS OR EQUAL 
!                                    INCREMENTS NOT EQUAL TO 1          
            IX = 1 
            IY = 1 
            IF (INCX.LT.0) IX = ( - N + 1) * INCX + 1 
            IF (INCY.LT.0) IY = ( - N + 1) * INCY + 1 
            DO 10 I = 1, N 
               ZY (IY) = ZX (IX) 
               IX = IX + INCX 
               IY = IY + INCY 
   10       END DO 
         ELSE 
!                                  CODE FOR BOTH INCREMENTS EQUAL TO 1  
            DO 20 I = 1, N 
               ZY (I) = ZX (I) 
   20       END DO 
         ENDIF 
      ENDIF 
      RETURN 
      END SUBROUTINE ZCOPY                          
!-----------------------------------------------------------------------
!  IMSL Name:  ZSWAP (Double precision version)                         
!                                                                       
!  Computer:   sgruxs/DOUBLE                                            
!                                                                       
!  Revised:    August 9, 1986                                           
!                                                                       
!  Purpose:    Interchange vectors X and Y, both complex(8).            
!                                                                       
!  Usage:      CALL ZSWAP (N, ZX, INCX, ZY, INCY)                       
!                                                                       
!  Arguments:                                                           
!     N      - Length of vectors X and Y.  (Input)                      
!     ZX     - complex(8) vector of length MAX(N*IABS(INCX),1).         
!              (Input/Output)                                           
!     INCX   - Displacement between elements of ZX.  (Input)            
!              X(I) is defined to be                                    
!                ZX(1+(I-1)*INCX) if INCX.GE.0  or                      
!                ZX(1+(I-N)*INCX) if INCX.LT.0.                         
!     ZY     - complex(8) vector of length MAX(N*IABS(INCY),1).         
!              (Input/Output)                                           
!     INCY   - Displacement between elements of ZY.  (Input)            
!              Y(I) is defined to be                                    
!                ZY(1+(I-1)*INCY) if INCY.GE.0  or                      
!                ZY(1+(I-N)*INCY) if INCY.LT.0.                         
!                                                                       
!  Keyword:    Level 1 BLAS; ZSWAP; Swap; Exchange                      
!                                                                       
!  GAMS:       D1a5                                                     
!                                                                       
!  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations              
!              STAT/LIBRARY Mathematical Support                        
!                                                                       
!  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE ZSWAP (N, ZX, INCX, ZY, INCY) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N, INCX, INCY 
      COMPLEX (8) ZX ( * ), ZY ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, IX, IY 
      COMPLEX (8) ZTEMP 
!                                                                       
      IF (N.GT.0) THEN 
         IF (INCX.NE.1.OR.INCY.NE.1) THEN 
!                                  CODE FOR UNEQUAL INCREMENTS OR EQUAL 
!                                    INCREMENTS NOT EQUAL TO 1          
            IX = 1 
            IY = 1 
            IF (INCX.LT.0) IX = ( - N + 1) * INCX + 1 
            IF (INCY.LT.0) IY = ( - N + 1) * INCY + 1 
            DO 10 I = 1, N 
               ZTEMP = ZX (IX) 
               ZX (IX) = ZY (IY) 
               ZY (IY) = ZTEMP 
               IX = IX + INCX 
               IY = IY + INCY 
   10       END DO 
         ELSE 
!                                  CODE FOR BOTH INCREMENTS EQUAL TO 1  
            DO 20 I = 1, N 
               ZTEMP = ZX (I) 
               ZX (I) = ZY (I) 
               ZY (I) = ZTEMP 
   20       END DO 
         ENDIF 
      ENDIF 
      RETURN 
      END SUBROUTINE ZSWAP                          
!-----------------------------------------------------------------------
!  IMSL Name:  ZAXPY (Double precision version)                         
!                                                                       
!  Computer:   sgruxs/DOUBLE                                            
!                                                                       
!  Revised:    August 9, 1986                                           
!                                                                       
!  Purpose:    Compute the scalar times a vector plus a vector,         
!              y = ax + y, all complex(8).                              
!                                                                       
!  Usage:      CALL ZAXPY (N, ZA, ZX, INCX, ZY, INCY)                   
!                                                                       
!  Arguments:                                                           
!     N      - Length of vectors X and Y.  (Input)                      
!     ZA     - complex(8) scalar.  (Input)                              
!     ZX     - complex(8) vector of length MAX(N*IABS(INCX),1).  (Input)
!     INCX   - Displacement between elements of ZX.  (Input)            
!              X(I) is defined to be                                    
!                 ZX(1+(I-1)*INCX) if INCX.GE.0  or                     
!                 ZX(1+(I-N)*INCX) if INCX.LT.0.                        
!     ZY     - complex(8) vector of length MAX(N*IABS(INCY),1).         
!                 (Input/Output)                                        
!              ZAXPY replaces Y(I) with ZA*X(I) + Y(I) for I = 1,...N.  
!              X(I) and Y(I) refer to specific elements of ZX and ZY.   
!     INCY   - Displacement between elements of ZY.  (Input)            
!              Y(I) is defined to be                                    
!                 ZY(1+(I-1)*INCY) if INCY.GE.0  or                     
!                 ZY(1+(I-N)*INCY) if INCY.LT.0.                        
!                                                                       
!  GAMS:       D1a7                                                     
!                                                                       
!  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations              
!              STAT/LIBRARY Mathematical Support                        
!                                                                       
!  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE ZAXPY (N, ZA, ZX, INCX, ZY, INCY) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N, INCX, INCY 
      COMPLEX (8) ZA, ZX ( * ), ZY ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, IX, IY 
      DOUBLEPRECISION ZZABS 
!                                  SPECIFICATIONS FOR INTRINSICS        
!     INTRINSIC  DABS,DBLE,DIMAG                                        
!      INTRINSIC  DABS, DBLE, DIMAG                                     
!      DOUBLE PRECISION DABS, DBLE, DIMAG                               
!                                                                       
      IF (N.GT.0) THEN 
         ZZABS = DABS (DBLE (ZA) ) + DABS (DIMAG (ZA) ) 
         IF (ZZABS.NE.0.0D0) THEN 
            IF (INCX.NE.1.OR.INCY.NE.1) THEN 
!                                  CODE FOR UNEQUAL INCREMENTS OR EQUAL 
!                                    INCREMENTS NOT EQUAL TO 1          
               IX = 1 
               IY = 1 
               IF (INCX.LT.0) IX = ( - N + 1) * INCX + 1 
               IF (INCY.LT.0) IY = ( - N + 1) * INCY + 1 
               DO 10 I = 1, N 
                  ZY (IY) = ZY (IY) + ZA * ZX (IX) 
                  IX = IX + INCX 
                  IY = IY + INCY 
   10          END DO 
            ELSE 
!                                  CODE FOR BOTH INCREMENTS EQUAL TO 1  
               DO 20 I = 1, N 
                  ZY (I) = ZY (I) + ZA * ZX (I) 
   20          END DO 
            ENDIF 
         ENDIF 
      ENDIF 
      RETURN 
      END SUBROUTINE ZAXPY                          
!-----------------------------------------------------------------------
!  IMSL Name:  ZSET (Double precision version)                          
!                                                                       
!  Computer:   sgruxs/DOUBLE                                            
!                                                                       
!  Revised:    August 9, 1986                                           
!                                                                       
!  Purpose:    Set the components of a vector to a scalar, all double   
!              complex(8).                                              
!                                                                       
!  Usage:      CALL ZSET (N, ZA, ZX, INCX)                              
!                                                                       
!  Arguments:                                                           
!     N      - Length of vector X.  (Input)                             
!     ZA     - complex(8) scalar.  (Input)                              
!     ZX     - complex(8) vector of length N*INCX.  (Input/Output)      
!              ZSET replaces X(I) with ZA for I=1,...,N. X(I) refers to 
!              a specific element of ZX. See INCX argument description. 
!     INCX   - Displacement between elements of ZX.  (Input)            
!              X(I) is defined to be ZX(1+(I-1)*INCX). INCX must be     
!              greater than zero.                                       
!                                                                       
!  GAMS:       D1a1                                                     
!                                                                       
!  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations              
!              STAT/LIBRARY Mathematical Support                        
!                                                                       
!  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE ZSET (N, ZA, ZX, INCX) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N, INCX 
      COMPLEX (8) ZA, ZX ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, NINCX 
!                                                                       
      IF (N.GT.0) THEN 
         IF (INCX.NE.1) THEN 
!                                  CODE FOR INCREMENT NOT EQUAL TO 1    
            NINCX = N * INCX 
            DO 10 I = 1, NINCX, INCX 
               ZX (I) = ZA 
   10       END DO 
         ELSE 
!                                  CODE FOR INCREMENT EQUAL TO 1        
            DO 20 I = 1, N 
               ZX (I) = ZA 
   20       END DO 
         ENDIF 
      ENDIF 
      RETURN 
      END SUBROUTINE ZSET                           
!-----------------------------------------------------------------------
!  IMSL Name:  DMACH (Double precision version)                         
!                                                                       
!  Computer:   sgruxs/DOUBLE                                            
!                                                                       
!  Revised:    March 15, 1984                                           
!                                                                       
!  Purpose:    Retrieve double precision machine constants.             
!                                                                       
!  Usage:      DMACH(N)                                                 
!                                                                       
!  Arguments:                                                           
!     N      - Index of desired constant.  (Input)                      
!     DMACH  - Machine constant.  (Output)                              
!              DMACH(1) = B**(EMIN-1), the smallest positive magnitude. 
!              DMACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude. 
!              DMACH(3) = B**(-T), the smallest relative spacing.       
!              DMACH(4) = B**(1-T), the largest relative spacing.       
!              DMACH(5) = LOG10(B), the log, base 10, of the radix.     
!              DMACH(6) = not-a-number.                                 
!              DMACH(7) = positive machine infinity.                    
!              DMACH(8) = negative machine infinity.                    
!                                                                       
!  GAMS:       R1                                                       
!                                                                       
!  Chapters:   MATH/LIBRARY Reference Material                          
!              STAT/LIBRARY Reference Material                          
!              SFUN/LIBRARY Reference Material                          
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      DOUBLEPRECISION FUNCTION DMACH (N) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N 
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      DOUBLEPRECISION RMACH (8) 
      SAVE RMACH 
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL E1MES, E1POP, E1PSH, E1STI 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER IRMACH (16) 
!                                                                       
      EQUIVALENCE (RMACH, IRMACH) 
!                                  DEFINE CONSTANTS                     
      DATA RMACH (1) / 2.22559D-308 / 
      DATA RMACH (2) / 1.79728D+308 / 
      DATA RMACH (3) / 1.11048D-16 / 
      DATA RMACH (4) / 2.22096D-16 / 
      DATA RMACH (5) / .3010299956639811952137388947245D0 / 
      DATA IRMACH (11) / 2146959360 / 
      DATA IRMACH (12) / 0 / 
      DATA IRMACH (13) / 2146435072 / 
      DATA IRMACH (14) / 0 / 
      DATA IRMACH (15) / - 1048576 / 
      DATA IRMACH (16) / 0 / 
!                                                                       
      IF (N.LT.1.OR.N.GT.8) THEN 
         CALL E1PSH ('DMACH ') 
         DMACH = RMACH (6) 
         CALL E1STI (1, N) 
      CALL E1MES (5, 5, 'The argument must be between 1 '//'and 8 inclus&
     &ive. N = %(I1)')                                                  
         CALL E1POP ('DMACH ') 
      ELSE 
         DMACH = RMACH (N) 
      ENDIF 
!                                                                       
      RETURN 
      END FUNCTION DMACH                            
!-----------------------------------------------------------------------
!  IMSL Name:  DZASUM (Double precision version)                        
!                                                                       
!  Computer:   sgruxs/DOUBLE                                            
!                                                                       
!  Revised:    August 9, 1986                                           
!                                                                       
!  Purpose:    Sum the absolute values of the real part together with   
!              the absolute values of the imaginary part of the         
!              components of a double-complex(8) vector.                
!                                                                       
!  Usage:      DZASUM(N, ZX, INCX)                                      
!                                                                       
!  Arguments:                                                           
!     N      - Length of vectors X.  (Input)                            
!     ZX     - complex(8) vector of length N*INCX.  (Input)             
!     INCX   - Displacement between elements of ZX.  (Input)            
!              X(I) is defined to be ZX(1+(I-1)*INCX).  INCX must be    
!              greater than 0.                                          
!     DZASUM - Sum from I=1 to N of DABS(REAL(X(I)))+DABS(AIMAG(X(I)))).
!              (Output)                                                 
!              X(I) refers to a specific element of ZX.                 
!                                                                       
!  Keyword:    Level 1 BLAS                                             
!                                                                       
!  GAMS:       D1a3a                                                    
!                                                                       
!  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations              
!              STAT/LIBRARY Mathematical Support                        
!                                                                       
!  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      DOUBLEPRECISION FUNCTION DZASUM (N, ZX, INCX) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N, INCX 
      COMPLEX (8) ZX ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, IX 
      DOUBLEPRECISION STEMP 
!                                  SPECIFICATIONS FOR INTRINSICS        
!     INTRINSIC  DABS,DBLE,DIMAG                                        
!      INTRINSIC  DABS, DBLE, DIMAG                                     
!      DOUBLE PRECISION DABS, DBLE, DIMAG                               
!                                                                       
      DZASUM = 0.0D0 
      STEMP = 0.0D0 
      IF (N.GT.0) THEN 
         IF (INCX.NE.1) THEN 
!                                  CODE FOR INCREMENTS NOT EQUAL TO 1   
            IX = 1 
            IF (INCX.LT.0) IX = ( - N + 1) * INCX + 1 
            DO 10 I = 1, N 
               STEMP = STEMP + DABS (DIMAG (ZX (IX) ) ) + DABS (DBLE (  &
               ZX (IX) ) )                                              
               IX = IX + INCX 
   10       END DO 
            DZASUM = STEMP 
         ELSE 
!                                  CODE FOR INCREMENTS EQUAL TO 1       
            DO 20 I = 1, N 
               STEMP = STEMP + DABS (DIMAG (ZX (I) ) ) + DABS (DBLE (ZX &
               (I) ) )                                                  
   20       END DO 
            DZASUM = STEMP 
         ENDIF 
      ENDIF 
      RETURN 
      END FUNCTION DZASUM                           
!-----------------------------------------------------------------------
!  IMSL Name:  ZDOTU (Double precision version)                         
!                                                                       
!  Computer:   SGRUXS/DOUBLE                                            
!                                                                       
!  Revised:    August 9, 1986                                           
!                                                                       
!  Purpose:    Compute the double-complex(8) dot product x*y.           
!                                                                       
!  Usage:      ZDOTU(N, ZX, INCX, ZY, INCY)                             
!                                                                       
!  Arguments:                                                           
!     N      - Length of vectors X and Y.  (Input)                      
!     ZX     - complex(8) vector of length MAX(N*IABS(INCX),1).         
!              (Input)                                                  
!     INCX   - Displacement between elements of ZX.  (Input)            
!              X(I) is defined to be                                    
!                 ZX(1+(I-1)*INCX) if INCX.GE.0  or                     
!                 ZX(1+(I-N)*INCX) if INCX.LT.0.                        
!     ZY     - complex(8) vector of length MAX(N*IABS(INCY),1).         
!              (Input)                                                  
!     INCY   - Displacement between elements of ZY.  (Input)            
!              Y(I) is defined to be                                    
!                 ZY(1+(I-1)*INCY) if INCY.GE.0  or                     
!                 ZY(1+(I-N)*INCY) if INCY.LT.0.                        
!     ZDOTU  - complex(8) sum from I=1 to N of X(I)*Y(I).  (Output)     
!              X(I) and Y(I) refer to specific elements of ZX and ZY    
!              respectively.                                            
!                                                                       
!  Keyword:    Level 1 BLAS; Inner product; Scalar product              
!                                                                       
!  GAMS:       D1a4                                                     
!                                                                       
!  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations              
!              STAT/LIBRARY Mathematical Support                        
!                                                                       
!  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      COMPLEX (8) FUNCTIONZDOTU (N, ZX, INCX, ZY, INCY) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N, INCX, INCY 
      COMPLEX (8) ZX ( * ), ZY ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, IX, IY 
      COMPLEX (8) ZTEMP 
!                                                                       
      ZDOTU = (0.0D0, 0.0D0) 
      IF (N.GT.0) THEN 
         ZTEMP = (0.0D0, 0.0D0) 
         IF (INCX.NE.1.OR.INCY.NE.1) THEN 
!                                  CODE FOR UNEQUAL INCREMENTS OR EQUAL 
!                                    INCREMENTS NOT EQUAL TO 1          
            IX = 1 
            IY = 1 
            IF (INCX.LT.0) IX = ( - N + 1) * INCX + 1 
            IF (INCY.LT.0) IY = ( - N + 1) * INCY + 1 
            DO 10 I = 1, N 
               ZTEMP = ZTEMP + ZX (IX) * ZY (IY) 
               IX = IX + INCX 
               IY = IY + INCY 
   10       END DO 
            ZDOTU = ZDOTU 
         ELSE 
!                                  CODE FOR BOTH INCREMENTS EQUAL TO 1  
            DO 20 I = 1, N 
               ZTEMP = ZTEMP + ZX (I) * ZY (I) 
   20       END DO 
         ENDIF 
         ZDOTU = ZTEMP 
      ENDIF 
      RETURN 
      END                                           
!-----------------------------------------------------------------------
!  IMSL Name:  IWKIN (Single precision version)                         
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    January 17, 1984                                         
!                                                                       
!  Purpose:    Initialize bookkeeping locations describing the          
!              workspace stack.                                         
!                                                                       
!  Usage:      CALL IWKIN (NSU)                                         
!                                                                       
!  Argument:                                                            
!     NSU    - Number of numeric storage units to which the workspace   
!              stack is to be initialized                               
!                                                                       
!  GAMS:       N4                                                       
!                                                                       
!  Chapters:   MATH/LIBRARY Reference Material                          
!              STAT/LIBRARY Reference Material                          
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE IWKIN (NSU) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER NSU 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER ISIZE (6), LALC, LBND, LBOOK, LMAX, LNEED, LNOW, LOUT,    &
      LUSED, MELMTS, MTYPE                                              
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      LOGICAL FIRST 
      SAVE FIRST 
!                                  SPECIFICATIONS FOR SPECIAL CASES     
!                                  SPECIFICATIONS FOR COMMON /WORKSP/   
      REAL RWKSP (5000) 
      REAL RDWKSP (5000) 
      DOUBLEPRECISION DWKSP (2500) 
      COMPLEX (8) CWKSP (2500) 
      COMPLEX (8) CZWKSP (2500) 
      COMPLEX (8) ZWKSP (1250) 
      INTEGER IWKSP (5000) 
      LOGICAL LWKSP (5000) 
      EQUIVALENCE (DWKSP (1), RWKSP (1) ) 
      EQUIVALENCE (CWKSP (1), RWKSP (1) ), (ZWKSP (1), RWKSP (1) ) 
      EQUIVALENCE (IWKSP (1), RWKSP (1) ), (LWKSP (1), RWKSP (1) ) 
      EQUIVALENCE (RDWKSP (1), RWKSP (1) ), (CZWKSP (1), RWKSP (1) ) 
      COMMON / WORKSP / DWKSP 
!                                  SPECIFICATIONS FOR EQUIVALENCE       
      EQUIVALENCE (LOUT, IWKSP (1) ) 
      EQUIVALENCE (LNOW, IWKSP (2) ) 
      EQUIVALENCE (LUSED, IWKSP (3) ) 
      EQUIVALENCE (LBND, IWKSP (4) ) 
      EQUIVALENCE (LMAX, IWKSP (5) ) 
      EQUIVALENCE (LALC, IWKSP (6) ) 
      EQUIVALENCE (LNEED, IWKSP (7) ) 
      EQUIVALENCE (LBOOK, IWKSP (8) ) 
      EQUIVALENCE (ISIZE (1), IWKSP (11) ) 
!                                  SPECIFICATIONS FOR INTRINSICS        
!     INTRINSIC  MAX0                                                   
!      INTRINSIC  MAX0                                                  
!      INTEGER    MAX0                                                  
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL E1MES, E1STI 
!                                                                       
      DATA FIRST / .TRUE. / 
!                                                                       
      IF (.NOT.FIRST) THEN 
         IF (NSU.NE.0) THEN 
            CALL E1STI (1, LMAX) 
      CALL E1MES (5, 100, 'Error from subroutine IWKIN:  '//'Workspace s&
     &tack has previously been '//'initialized to %(I1). Correct by maki&
     &ng the '//'call to IWKIN the first executable '//'statement in the&
     & main program.  ')                                                
!                                                                       
            STOP 
!                                                                       
         ELSE 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
      IF (NSU.EQ.0) THEN 
!                                  IF NSU=0 USE DEFAULT SIZE 5000       
         MELMTS = 5000 
      ELSE 
         MELMTS = NSU 
      ENDIF 
!                                  NUMBER OF ITEMS .LT. 0               
      IF (MELMTS.LE.0) THEN 
         CALL E1STI (1, MELMTS) 
      CALL E1MES (5, 1, 'Error from subroutine IWKIN:  Number '//'of num&
     &eric storage units is not positive. NSU '//'= %(I1) ')            
      ELSE 
!                                                                       
         FIRST = .FALSE. 
!                                  HERE TO INITIALIZE                   
!                                                                       
!                                  SET DATA SIZES APPROPRIATE FOR A     
!                                  STANDARD CONFORMING FORTRAN SYSTEM   
!                                  USING THE FORTRAN                    
!                                  *NUMERIC STORAGE UNIT* AS THE        
!                                  MEASURE OF SIZE.                     
!                                                                       
!                                  TYPE IS REAL                         
         MTYPE = 3 
!                                  LOGICAL                              
         ISIZE (1) = 1 
!                                  INTEGER                              
         ISIZE (2) = 1 
!                                  REAL                                 
         ISIZE (3) = 1 
!                                  DOUBLE PRECISION                     
         ISIZE (4) = 2 
!                                  complex(8)                           
         ISIZE (5) = 2 
!                                  complex(8)                           
         ISIZE (6) = 4 
!                                  NUMBER OF WORDS USED FOR BOOKKEEPING 
         LBOOK = 16 
!                                  CURRENT ACTIVE LENGTH OF THE STACK   
         LNOW = LBOOK 
!                                  MAXIMUM VALUE OF LNOW ACHIEVED THUS  
!                                  FAR                                  
         LUSED = LBOOK 
!                                  MAXIMUM LENGTH OF THE STORAGE ARRAY  
         LMAX = MAX0 (MELMTS, ( (LBOOK + 2) * ISIZE (2) + ISIZE (3)     &
         - 1) / ISIZE (3) )                                             
!                                  LOWER BOUND OF THE PERMANENT STORAGE 
!                                  WHICH IS ONE WORD MORE THAN THE      
!                                  MAXIMUM ALLOWED LENGTH OF THE STACK  
         LBND = LMAX + 1 
!                                  NUMBER OF CURRENT ALLOCATIONS        
         LOUT = 0 
!                                  TOTAL NUMBER OF ALLOCATIONS MADE     
         LALC = 0 
!                                  NUMBER OF WORDS BY WHICH THE ARRAY   
!                                  SIZE MUST BE INCREASED FOR ALL PAST  
!                                  ALLOCATIONS TO SUCCEED               
         LNEED = 0 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE IWKIN                          
!-----------------------------------------------------------------------
!  IMSL Name:  I1X (Single precision version)                           
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    August 30, 1985                                          
!                                                                       
!  Purpose:    Determine the array subscript indicating the starting    
!              element at which a key character sequence begins.        
!              (Case-sensitive version)                                 
!                                                                       
!  Usage:      I1X(CHRSTR, I1LEN, KEY, KLEN)                            
!                                                                       
!  Arguments:                                                           
!     CHRSTR - Character array to be searched.  (Input)                 
!     I1LEN  - Length of CHRSTR.  (Input)                               
!     KEY    - Character array that contains the key sequence.  (Input) 
!     KLEN   - Length of KEY.  (Input)                                  
!     I1X    - Integer function.  (Output)                              
!                                                                       
!  Remarks:                                                             
!  1. Returns zero when there is no match.                              
!                                                                       
!  2. Returns zero if KLEN is longer than ISLEN.                        
!                                                                       
!  3. Returns zero when any of the character arrays has a negative or   
!     zero length.                                                      
!                                                                       
!  GAMS:       N5c                                                      
!                                                                       
!  Chapter:    MATH/LIBRARY Utilities                                   
!                                                                       
!  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      INTEGER FUNCTION I1X (CHRSTR, I1LEN, KEY, KLEN) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER I1LEN, KLEN 
      CHARACTER CHRSTR ( * ), KEY ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, II, J 
!                                                                       
      I1X = 0 
      IF (KLEN.LE.0.OR.I1LEN.LE.0) GOTO 9000 
      IF (KLEN.GT.I1LEN) GOTO 9000 
!                                                                       
      I = 1 
      II = I1LEN - KLEN + 1 
   10 IF (I.LE.II) THEN 
         IF (CHRSTR (I) .EQ.KEY (1) ) THEN 
            DO 20 J = 2, KLEN 
               IF (CHRSTR (I + J - 1) .NE.KEY (J) ) GOTO 30 
   20       END DO 
            I1X = I 
            GOTO 9000 
   30       CONTINUE 
         ENDIF 
         I = I + 1 
         GOTO 10 
      ENDIF 
!                                                                       
 9000 RETURN 
      END FUNCTION I1X                              
!-----------------------------------------------------------------------
!  IMSL Name:  UMACH (Single precision version)                         
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    March 21, 1984                                           
!                                                                       
!  Purpose:    Set or retrieve input or output device unit numbers.     
!                                                                       
!  Usage:      CALL UMACH (N, NUNIT)                                    
!                                                                       
!  Arguments:                                                           
!     N      - Index of desired unit.  (Input)                          
!              The values of N are defined as follows:                  
!              N = 1, corresponds to the standard input unit.           
!              N = 2, corresponds to the standard output unit.          
!     NUNIT  - I/O unit.  (Input or Output)                             
!              If the value of N is negative, the unit corresponding    
!              to the index is reset to the value given in NUNIT.       
!              Otherwise, the value corresponding to the index is       
!              returned in NUNIT.                                       
!                                                                       
!  GAMS:       R1                                                       
!                                                                       
!  Chapters:   MATH/LIBRARY Reference Material                          
!              STAT/LIBRARY Reference Material                          
!              SFUN/LIBRARY Reference Material                          
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE UMACH (N, NUNIT) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N, NUNIT 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER NN, NOUT 
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      INTEGER UNIT (2) 
      SAVE UNIT 
!                                  SPECIFICATIONS FOR INTRINSICS        
!     INTRINSIC  IABS                                                   
!      INTRINSIC  IABS                                                  
!      INTEGER    IABS                                                  
!                                                                       
      DATA UNIT (1) / 5 / 
      DATA UNIT (2) / 6 / 
!                                                                       
      NN = IABS (N) 
      IF (NN.NE.1.AND.NN.NE.2) THEN 
!                                  ERROR.  INVALID RANGE FOR N.         
         NOUT = UNIT (2) 
         WRITE (NOUT, 99999) NN 
99999 FORMAT    (/, ' *** TERMINAL ERROR 5 from UMACH.  The absolute',  &
     &          /, ' ***          value of the index variable must be'  &
     &          , /, ' ***          1 or 2.  IABS(N) = ', I6,           &
     &          '.', /)                                                 
         STOP 
!                                  CHECK FOR RESET OR RETRIEVAL         
      ELSEIF (N.LT.0) THEN 
!                                  RESET                                
         UNIT (NN) = NUNIT 
      ELSE 
!                                  RETRIEVE                             
         NUNIT = UNIT (N) 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE UMACH                          
!-----------------------------------------------------------------------
!  IMSL Name:  S1ANUM                                                   
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    March 28, 1984                                           
!                                                                       
!  Purpose:    Scan a token and identify it as follows: integer, real   
!              number (single/double), FORTRAN relational operator,     
!              FORTRAN logical operator, or FORTRAN logical constant.   
!                                                                       
!  Usage:      CALL S1ANUM(INSTR, SLEN, CODE, OLEN)                     
!                                                                       
!  Arguments:                                                           
!     INSTR  - Character string to be scanned.  (Input)                 
!     SLEN   - Length of INSTR.  (Input)                                
!     CODE   - Token code.  (Output)  Where                             
!                 CODE =  0  indicates an unknown token,                
!                 CODE =  1  indicates an integer number,               
!                 CODE =  2  indicates a (single precision) real number,
!                 CODE =  3  indicates a (double precision) real number,
!                 CODE =  4  indicates a logical constant (.TRUE. or    
!                               .FALSE.),                               
!                 CODE =  5  indicates the relational operator .EQ.,    
!                 CODE =  6  indicates the relational operator .NE.,    
!                 CODE =  7  indicates the relational operator .LT.,    
!                 CODE =  8  indicates the relational operator .LE.,    
!                 CODE =  9  indicates the relational operator .GT.,    
!                 CODE = 10  indicates the relational operator .GE.,    
!                 CODE = 11  indicates the logical operator .AND.,      
!                 CODE = 12  indicates the logical operator .OR.,       
!                 CODE = 13  indicates the logical operator .EQV.,      
!                 CODE = 14  indicates the logical operator .NEQV.,     
!                 CODE = 15  indicates the logical operator .NOT..      
!     OLEN   - Length of the token as counted from the first character  
!              in INSTR.  (Output)  OLEN returns a zero for an unknown  
!              token (CODE = 0).                                        
!                                                                       
!  Remarks:                                                             
!  1. Blanks are considered significant.                                
!  2. Lower and upper case letters are not significant.                 
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All rights reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      SUBROUTINE S1ANUM (INSTR, SLEN, CODE, OLEN) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER SLEN, CODE, OLEN 
      CHARACTER INSTR ( * ) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER I, IBEG, IIBEG, J 
      LOGICAL FLAG 
      CHARACTER CHRSTR (6) 
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      INTEGER TABPTR (16), TDCNST, TICNST, TOKEN (13), TRCNST, TZERR 
      CHARACTER DIGIT (10), LETTER (52), MINUS, PERIOD, PLUS, TABLE (38) 
      SAVE DIGIT, LETTER, MINUS, PERIOD, PLUS, TABLE, TABPTR, TDCNST,   &
      TICNST, TOKEN, TRCNST, TZERR                                      
!                                  SPECIFICATIONS FOR FUNCTIONS         
      EXTERNAL I1X, I1CSTR 
      INTEGER I1X, I1CSTR 
!                                                                       
      DATA TOKEN / 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 4, 4 / 
      DATA TABLE / 'D', 'E', 'E', 'Q', 'N', 'E', 'L', 'T', 'L', 'E',    &
      'G', 'T', 'G', 'E', 'A', 'N', 'D', 'O', 'R', 'E', 'Q', 'V', 'N',  &
      'E', 'Q', 'V', 'N', 'O', 'T', 'T', 'R', 'U', 'E', 'F', 'A', 'L',  &
      'S', 'E' /                                                        
      DATA TABPTR / 1, 2, 3, 5, 7, 9, 11, 13, 15, 18, 20, 23, 27, 30,   &
      34, 39 /                                                          
      DATA DIGIT / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      DATA LETTER / 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',   &
      'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W',  &
      'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',  &
      'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w',  &
      'x', 'y', 'z' /                                                   
      DATA PERIOD / '.' / , PLUS / '+' / , MINUS / '-' / 
      DATA TZERR / 0 /, TICNST / 1 / 
      DATA TRCNST / 2 /, TDCNST / 3 / 
!                                                                       
      IF (SLEN.LE.0) THEN 
         CODE = 0 
         OLEN = 0 
         RETURN 
      ENDIF 
!                                  STATE 0 - ASSUME ERROR TOKEN         
      IBEG = 1 
      CODE = TZERR 
!                                  CHECK SIGN                           
      IF (INSTR (IBEG) .EQ.MINUS.OR.INSTR (IBEG) .EQ.PLUS) THEN 
         FLAG = .TRUE. 
         IIBEG = IBEG 
         IBEG = IBEG + 1 
      ELSE 
         FLAG = .FALSE. 
      ENDIF 
!                                  STATE 1 - ASSUME INTEGER CONSTANT    
      IF (I1X (DIGIT, 10, INSTR (IBEG), 1) .NE.0) THEN 
         CODE = TICNST 
         IIBEG = IBEG 
         IBEG = IBEG + 1 
!                                                                       
   10    IF (IBEG.LE.SLEN) THEN 
!                                                                       
            IF (I1X (DIGIT, 10, INSTR (IBEG), 1) .NE.0) THEN 
               IIBEG = IBEG 
               IBEG = IBEG + 1 
               GOTO 10 
!                                                                       
            ENDIF 
!                                                                       
         ELSE 
            GOTO 80 
!                                                                       
         ENDIF 
!                                                                       
         IF (INSTR (IBEG) .NE.PERIOD) GOTO 80 
      ENDIF 
!                                  STATE 2 - ASSUME REAL CONSTANT       
      IF (CODE.EQ.TICNST) THEN 
         CODE = TRCNST 
         IIBEG = IBEG 
         IBEG = IBEG + 1 
         IF (IBEG.GT.SLEN) GOTO 80 
      ELSEIF (INSTR (IBEG) .EQ.PERIOD.AND.SLEN.GE.2) THEN 
         IF (I1X (DIGIT, 10, INSTR (IBEG + 1), 1) .NE.0) THEN 
            CODE = TRCNST 
            IIBEG = IBEG + 1 
            IBEG = IBEG + 2 
            IF (IBEG.GT.SLEN) GOTO 80 
         ENDIF 
      ENDIF 
!                                                                       
      IF (I1X (DIGIT, 10, INSTR (IBEG), 1) .NE.0) THEN 
         CODE = TRCNST 
         IIBEG = IBEG 
         IBEG = IBEG + 1 
!                                                                       
   20    IF (IBEG.LE.SLEN) THEN 
!                                                                       
            IF (I1X (DIGIT, 10, INSTR (IBEG), 1) .NE.0) THEN 
               IIBEG = IBEG 
               IBEG = IBEG + 1 
               GOTO 20 
!                                                                       
            ENDIF 
!                                                                       
         ELSE 
            GOTO 80 
!                                                                       
         ENDIF 
!                                                                       
      ENDIF 
!                                                                       
      IF (CODE.EQ.TZERR) THEN 
         IF (INSTR (IBEG) .NE.PERIOD) GOTO 80 
         IBEG = IBEG + 1 
         IF (IBEG.GT.SLEN) GOTO 80 
      ENDIF 
!                                                                       
      IF (I1X (LETTER, 52, INSTR (IBEG), 1) .EQ.0) GOTO 80 
      CHRSTR (1) = INSTR (IBEG) 
!                                                                       
      DO 30 I = 2, 6 
         IBEG = IBEG + 1 
         IF (IBEG.GT.SLEN) GOTO 80 
         IF (I1X (LETTER, 52, INSTR (IBEG), 1) .EQ.0) GOTO 40 
         CHRSTR (I) = INSTR (IBEG) 
   30 END DO 
!                                                                       
      GOTO 80 
!                                                                       
   40 CONTINUE 
!                                                                       
      DO 50 J = 1, 15 
         IF (I1CSTR (CHRSTR, I - 1, TABLE (TABPTR (J) ), TABPTR (J + 1) &
         - TABPTR (J) ) .EQ.0) GOTO 60                                  
   50 END DO 
!                                                                       
      GOTO 80 
!                                  STATE 4 - LOGICAL OPERATOR           
   60 IF (J.GT.2) THEN 
!                                                                       
         IF (CODE.EQ.TRCNST) THEN 
!                                                                       
            IF (INSTR (IBEG) .EQ.PERIOD) THEN 
               CODE = TICNST 
               IIBEG = IIBEG - 1 
            ENDIF 
!                                                                       
            GOTO 80 
!                                                                       
         ELSEIF (INSTR (IBEG) .NE.PERIOD) THEN 
            GOTO 80 
!                                                                       
         ELSEIF (FLAG) THEN 
            GOTO 80 
!                                                                       
         ELSE 
            CODE = TOKEN (J - 2) 
            IIBEG = IBEG 
            GOTO 80 
!                                                                       
         ENDIF 
!                                                                       
      ENDIF 
!                                  STATE 5 - DOUBLE PRECISION CONSTANT  
      IF (CODE.NE.TRCNST) GOTO 80 
      IF (INSTR (IBEG) .EQ.MINUS.OR.INSTR (IBEG) .EQ.PLUS) IBEG = IBEG +&
      1                                                                 
      IF (IBEG.GT.SLEN) GOTO 80 
!                                                                       
      IF (I1X (DIGIT, 10, INSTR (IBEG), 1) .EQ.0) THEN 
         GOTO 80 
!                                                                       
      ELSE 
         IIBEG = IBEG 
         IBEG = IBEG + 1 
!                                                                       
   70    IF (IBEG.LE.SLEN) THEN 
!                                                                       
            IF (I1X (DIGIT, 10, INSTR (IBEG), 1) .NE.0) THEN 
               IIBEG = IBEG 
               IBEG = IBEG + 1 
               GOTO 70 
!                                                                       
            ENDIF 
!                                                                       
         ENDIF 
!                                                                       
      ENDIF 
!                                                                       
      IF (J.EQ.1) CODE = TDCNST 
!                                                                       
   80 CONTINUE 
!                                                                       
      IF (CODE.EQ.TZERR) THEN 
         OLEN = 0 
!                                                                       
      ELSE 
         OLEN = IIBEG 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE S1ANUM                         
!-----------------------------------------------------------------------
!  IMSL Name:  ICASE (Single precision version)                         
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    September 9, 1985                                        
!                                                                       
!  Purpose:    Convert from character to the integer ASCII value without
!              regard to case.                                          
!                                                                       
!  Usage:      ICASE(CH)                                                
!                                                                       
!  Arguments:                                                           
!     CH     - Character to be converted.  (Input)                      
!     ICASE  - Integer ASCII value for CH without regard to the case    
!              of CH.  (Output)                                         
!              ICASE returns the same value as IMSL routine IACHAR for  
!              all but lowercase letters.  For these, it returns the    
!              IACHAR value for the corresponding uppercase letter.     
!                                                                       
!  GAMS:       N3                                                       
!                                                                       
!  Chapter:    MATH/LIBRARY Utilities                                   
!              STAT/LIBRARY Utilities                                   
!                                                                       
!  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      INTEGER FUNCTION ICASE (CH) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      CHARACTER CH 
!                                  SPECIFICATIONS FOR FUNCTIONS         
!      EXTERNAL   IACHAR                                                
!      INTEGER    IACHAR                                                
!                                                                       
      ICASE = IACHAR (CH) 
      IF (ICASE.GE.97.AND.ICASE.LE.122) ICASE = ICASE-32 
!                                                                       
      RETURN 
      END FUNCTION ICASE                            
!-----------------------------------------------------------------------
!  IMSL Name:  I1CSTR (Single precision version)                        
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    September 10, 1985                                       
!                                                                       
!  Purpose:    Case insensitive comparison of two character arrays.     
!                                                                       
!  Usage:      I1CSTR(STR1, LEN1, STR2, LEN2)                           
!                                                                       
!  Arguments:                                                           
!     STR1   - First character array.  (Input)                          
!     LEN1   - Length of STR1.  (Input)                                 
!     STR2   - Second character array.  (Input)                         
!     LEN2   - Length of STR2.  (Input)                                 
!     I1CSTR - Integer function.  (Output) Where                        
!              I1CSTR = -1  if STR1 .LT. STR2,                          
!              I1CSTR =  0  if STR1 .EQ. STR2,                          
!              I1CSTR =  1  if STR1 .GT. STR2.                          
!                                                                       
!  Remarks:                                                             
!  1. If the two arrays, STR1 and STR2,  are of unequal length, the     
!     shorter array is considered as if it were extended with blanks    
!     to the length of the longer array.                                
!                                                                       
!  2. If one or both lengths are zero or negative the I1CSTR output is  
!     based on comparison of the lengths.                               
!                                                                       
!  GAMS:       N5c                                                      
!                                                                       
!  Chapter:    MATH/LIBRARY Utilities                                   
!                                                                       
!  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      INTEGER FUNCTION I1CSTR (STR1, LEN1, STR2, LEN2) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER LEN1, LEN2 
      CHARACTER STR1 (LEN1), STR2 (LEN2) 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER IC1, IC2, ICB, IS, L, LENM 
!                                  SPECIFICATIONS FOR INTRINSICS        
!     INTRINSIC  ISIGN,MIN0                                             
!      INTRINSIC  ISIGN, MIN0                                           
!      INTEGER    ISIGN, MIN0                                           
!                                  SPECIFICATIONS FOR FUNCTIONS         
      EXTERNAL ICASE 
      INTEGER ICASE 
!                                                                       
      IF (LEN1.GT.0.AND.LEN2.GT.0) THEN 
!                                  COMPARE FIRST LENM CHARACTERS        
         LENM = MIN0 (LEN1, LEN2) 
         DO 10 L = 1, LENM 
            IC1 = ICASE (STR1 (L) ) 
            IC2 = ICASE (STR2 (L) ) 
            IF (IC1.NE.IC2) THEN 
               I1CSTR = ISIGN (1, IC1 - IC2) 
               RETURN 
            ENDIF 
   10    END DO 
      ENDIF 
!                                  COMPARISON BASED ON LENGTH OR        
!                                  TRAILING BLANKS                      
      IS = LEN1 - LEN2 
      IF (IS.EQ.0) THEN 
         I1CSTR = 0 
      ELSE 
         IF (LEN1.LE.0.OR.LEN2.LE.0) THEN 
!                                  COMPARISON BASED ON LENGTH           
            I1CSTR = ISIGN (1, IS) 
         ELSE 
!                                  COMPARISON BASED ON TRAILING BLANKS  
!                                  TO EXTEND SHORTER ARRAY              
            LENM = LENM + 1 
            ICB = ICASE (' ') 
            IF (IS.GT.0) THEN 
!                                  EXTEND STR2 WITH BLANKS              
               DO 20 L = LENM, LEN1 
                  IC1 = ICASE (STR1 (L) ) 
                  IF (IC1.NE.ICB) THEN 
                     I1CSTR = ISIGN (1, IC1 - ICB) 
                     RETURN 
                  ENDIF 
   20          END DO 
            ELSE 
!                                  EXTEND STR1 WITH BLANKS              
               DO 30 L = LENM, LEN2 
                  IC2 = ICASE (STR2 (L) ) 
                  IF (ICB.NE.IC2) THEN 
                     I1CSTR = ISIGN (1, ICB - IC2) 
                     RETURN 
                  ENDIF 
   30          END DO 
            ENDIF 
!                                                                       
            I1CSTR = 0 
         ENDIF 
      ENDIF 
!                                                                       
      RETURN 
      END FUNCTION I1CSTR                           
!-----------------------------------------------------------------------
!  IMSL Name:  IMACH (Single precision version)                         
!                                                                       
!  Computer:   sgruxs/SINGLE                                            
!                                                                       
!  Revised:    March 26, 1984                                           
!                                                                       
!  Purpose:    Retrieve integer machine constants.                      
!                                                                       
!  Usage:      IMACH(N)                                                 
!                                                                       
!  Arguments:                                                           
!     N      - Index of desired constant.  (Input)                      
!     IMACH  - Machine constant.  (Output)                              
!                                                                       
!  Remark:                                                              
!     Following is a description of the assorted integer machine        
!     constants.                                                        
!                                                                       
!     Words                                                             
!                                                                       
!        IMACH( 1) = Number of bits per integer storage unit.           
!        IMACH( 2) = Number of characters per integer storage unit.     
!                                                                       
!     Integers                                                          
!                                                                       
!        Assume integers are represented in the S-DIGIT, BASE-A form    
!        SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )                 
!        where 0 .LE. X(I) .LT. A for I=0,...,S-1.  Then                
!                                                                       
!        IMACH( 3) = A, the base.                                       
!        IMACH( 4) = S, number of BASE-A digits.                        
!        IMACH( 5) = A**S - 1, largest magnitude.                       
!                                                                       
!     Floating-point numbers                                            
!                                                                       
!        Assume floating-point numbers are represented in the T-DIGIT,  
!        BASE-B form SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )       
!        where 0 .LE. X(I) .LT. B for I=1,...,T,                        
!        0 .LT. X(1), and EMIN .LE. E .LE. EMAX.  Then                  
!                                                                       
!        IMACH( 6) = B, the base.                                       
!                                                                       
!        Single precision                                               
!                                                                       
!           IMACH( 7) = T, number of BASE-B digits.                     
!           IMACH( 8) = EMIN, smallest exponent E.                      
!           IMACH( 9) = EMAX, largest exponent E.                       
!                                                                       
!        Double precision                                               
!                                                                       
!           IMACH(10) = T, number of BASE-B digits.                     
!           IMACH(11) = EMIN, smallest exponent E.                      
!           IMACH(12) = EMAX, largest exponent E.                       
!                                                                       
!  GAMS:       R1                                                       
!                                                                       
!  Chapters:   MATH/LIBRARY Reference Material                          
!              STAT/LIBRARY Reference Material                          
!              SFUN/LIBRARY Reference Material                          
!                                                                       
!  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.                 
!                                                                       
!  Warranty:   IMSL warrants only that IMSL testing has been applied    
!              to this code.  No other warranty, expressed or implied,  
!              is applicable.                                           
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      INTEGER FUNCTION IMACH (N) 
!                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER N 
!                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER NOUT 
!                                  SPECIFICATIONS FOR SAVE VARIABLES    
      INTEGER IMACHV (12) 
      SAVE IMACHV 
!                                  SPECIFICATIONS FOR SUBROUTINES       
      EXTERNAL UMACH 
!                                  DEFINE CONSTANTS                     
      DATA IMACHV (1) / 32 / 
      DATA IMACHV (2) / 4 / 
      DATA IMACHV (3) / 2 / 
      DATA IMACHV (4) / 31 / 
      DATA IMACHV (5) / 2147483647 / 
      DATA IMACHV (6) / 2 / 
      DATA IMACHV (7) / 24 / 
      DATA IMACHV (8) / - 125 / 
      DATA IMACHV (9) / 128 / 
      DATA IMACHV (10) / 53 / 
      DATA IMACHV (11) / - 1021 / 
      DATA IMACHV (12) / 1024 / 
!                                                                       
      IF (N.LT.1.OR.N.GT.12) THEN 
!                                  ERROR.  INVALID RANGE FOR N.         
         CALL UMACH (2, NOUT) 
         WRITE (NOUT, 99999) N 
99999 FORMAT    (/, ' *** TERMINAL ERROR 5 from IMACH.  The argument',  &
     &          /, ' ***          must be between 1 and 12 inclusive.'  &
     &          , /, ' ***          N = ', I6, '.', /)                  
         IMACH = 0 
         STOP 
!                                                                       
      ELSE 
         IMACH = IMACHV (N) 
      ENDIF 
!                                                                       
      RETURN 
      END FUNCTION IMACH                            
