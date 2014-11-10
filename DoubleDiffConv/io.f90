module io
   use parameters
   implicit none
contains

      subroutine readInputFile(inputfile)
      implicit none
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
      IMPLICIT none
      CHARACTER(len=*):: outputfile
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
*************************************************************************
*     opens file <filename> and puts the filepointer at EOF
      subroutine open_file_at_end(NHANDLE,filename)
      INTEGER:: NHANDLE
      CHARACTER(len=*):: filename

      OPEN(NHANDLE,FILE=trim(filename),STATUS='OLD', POSITION='APPEND',ERR=990)
      GOTO 999
990   WRITE(*,*) 'Error reading ',filename
      STOP ERR_WRT_OUTFILE
999   CONTINUE
      END subroutine
end module io
