      PROGRAM lolo 
      IMPLICIT none 
      CHARACTER(60) infile, outfile 
                                                                        
!---------------------------------------------------------              
!     arg #1 - filename or usage ?                                      
      CALL getarg (1, infile) 
      IF (infile.eq.' ') then 
         PRINT * , 'Usage : ' 
         PRINT * , 'lo <in file> <out file>' 
         STOP 
      ENDIF 
      IF (infile.eq.'-h') then 
         PRINT * , 'Usage : ' 
         PRINT * , 'lo <in file> <out file>' 
         STOP 
      ENDIF 
                                                                        
      CALL getarg (2, outfile) 
      PRINT * , trim (infile) , ' - ', trim (outfile) 
!     call losub('in','out')                                            
      CALL losub (trim (infile), trim (outfile) ) 
      END PROGRAM lolo                              
