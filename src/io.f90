module io
#include "errorcodes.h"
   use parameters
   use parser
   implicit none
contains

   subroutine readConfigFile(inputfile)
      IMPLICIT none
      CHARACTER(len=*) inputfile
      OPEN(15,FILE=inputfile,STATUS='OLD',ERR=10)
      GOTO 11
10    WRITE(*,*) 'LOSUB.F: Error while reading inputfile!'
      STOP   NO_INFILE
11    CONTINUE
      READ(15,'(A)',END=15) 
      READ(15,*,END=15) Symmetry,LCALC
      READ(15,'(A)',END=15) 
      READ(15,*,END=15) Rt,Tau,Pt,ETA,Le,Rc
      READ(15,'(A)',END=15) 
      READ(15,*,END=15) Truncation,M0
      READ(15,'(A)',END=15) 
      READ(15,*,END=15) DRt,ABSE,RELE,NSMAX
      READ(15,'(A)',END=15) 
      READ(15,*,END=15) StepSize, UpperLimit
      CLOSE(15)
      GOTO 16
15    WRITE(*,*) 'Error in inputfile ',inputfile
      STOP ERR_IN_INFILE
16    CONTINUE
   end subroutine

   !**********************************************************************
   subroutine readConfigFileNew(inputfile)
      implicit none
      CHARACTER(len=*) inputfile
      CHARACTER(len=60) varname
      CHARACTER(len=256) line
      integer:: err
      OPEN(15,FILE=inputfile,STATUS='OLD', iostat=err)
      if(err.ne.0) then
         WRITE(*,*) 'Error opening input file!'
         STOP   NO_INFILE
      endif
      do
         call parse(15, varname, line, err)
         if (err.ne.0) exit
         select case(varname)
            case('Calculation')
               call read_val(line, LCALC)
            case('VariablePar')
               call read_val(line, VariablePar)
            case('Symmetry')
               call read_val(line, Symmetry)
            case('Truncation')
               call read_val(line, Truncation)
            case('Rt')
               call read_val(line, Rt)
            case('Rc')
               call read_val(line, Rc)
            case('Pt')
               call read_val(line, Pt)
            case('Le')
               call read_val(line, Le)
            case('m0')
               call read_val(line, m0)
            case('eta')
               call read_val(line, eta)
            case('tau')
               call read_val(line, tau)
            case('StepSize')
               call read_val(line, StepSize)
            case('UpperLimit')
               call read_val(line, UpperLimit)
            case('AbsParameterError')
               call read_val(line, ABSE)
            case('RelativeGREror')
               call read_val(line, RELE)
            case('MaxIterations')
               call read_val(line, NSMAX)
            case default
               cycle
         end select
      enddo
      close(15)
   end subroutine
   !**********************************************************************
   subroutine writeOutputHeader(unitOut)
      use parameters
      IMPLICIT none
      integer, intent(in):: unitOut
      IF(LCALC.GT.0 .AND. LCALC.LT.4 .or. LCALC.eq.5.or.LCALC.eq.6) THEN
         WRITE(unitOut,*)  '### Output of Program lo.f Ver.2.1:        ###'
         WRITE(unitOut,*)  '### Lin. Onset of Conv. via Galerkinmethod ###'
         WRITE(unitOut,'(A11,E12.5,A2)') '# P     ', Pt,     '#'
         WRITE(unitOut,'(A11,E12.5,A2)') '# Lewis ', Le,     '#'
         WRITE(unitOut,'(A11,E12.5,A2)') '# TAU   ', TAU,    '#'
         WRITE(unitOut,'(A11,E12.5,A2)') '# R     ', Rt,     '#'
         WRITE(unitOut,'(A11,E12.5,A2)') '# RC    ', Rc,     '#'
         WRITE(unitOut,'(A11,E12.5,A2)') '# ETA   ', ETA,    '#'
         WRITE(unitOut,'(A11,G12.5,A2)') '# m     ', M0,     '#'
         WRITE(unitOut,'(A11,I12,A2)')   '# Symmetry     ', Symmetry,   '#'
         WRITE(unitOut,'(A11,E12.5,A2)') '# LowerLimit   ', LowerLimit, '#'
         WRITE(unitOut,'(A11,E12.5,A2)') '# UpperLimit   ', UpperLimit, '#'
         WRITE(unitOut,'(A11,E12.5,A2)') '# StepSize     ', StepSize,   '#'
         WRITE(unitOut,'(A11,I12,A2)')   '# Truncation   ', Truncation, '#'
         WRITE(unitOut,*)  '# see definition of LCALC for output. LCALC:', LCALC,'   #'
         WRITE(unitOut,*)  '#                                      #'
      ENDIF
   end subroutine

! *************************************************************************
!     opens file <filename> and puts the filepointer at EOF
      subroutine open_file_at_end(NHANDLE,filename)
         implicit none
         INTEGER, intent(in):: NHANDLE
         CHARACTER(len=*):: filename
         
         OPEN(NHANDLE,FILE=trim(filename),STATUS='OLD', POSITION='APPEND',ERR=990)
         GOTO 999
990      WRITE(*,*) 'Error reading ',filename
         STOP ERR_WRT_OUTFILE
999      CONTINUE
      END subroutine

! *************************************************************************
      subroutine writeConfigFile(outputfile)
         implicit none
         CHARACTER(len=*), intent(in):: outputfile

         OPEN(99,FILE=outputfile,STATUS='UNKNOWN')
         WRITE(99,*) ' Symmetry (0/1/2) | LCALC (1/2/3/4) |'
         WRITE(99,'(A,2I12)') ' ',Symmetry, LCALC
         WRITE(99,*) '|  RAYLEIGH  |  TAU     |  PRANTEL  |  ETA  | Lewis |   Rconc   |'
         WRITE(99,'(1P,E17.6,5(A,E17.6))') Rt,' ',TAU,' ',Pt,' ',ETA,' ',Le,' ',Rc
         WRITE(99,*) '|   NTRUNC (>=1) | MODE |'
         WRITE(99,'(A,2I12)') ' ',Truncation, M0
         WRITE(99,*) '|   DRA   | ABSERR  |  RELERR  | NSMAX |'
         WRITE(99,'(1PG13.6,A,1PG12.5,A,1PG12.5,A,I4)') DRt,' ',ABSE,' ',RELE,' ',NSMAX
         WRITE(99,*) '|  StepSize  | UpperLimit'
         WRITE(99,'(1P,2G11.4)') StepSize, UpperLimit
         CLOSE(99)
      end subroutine writeConfigFile
end module io
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
