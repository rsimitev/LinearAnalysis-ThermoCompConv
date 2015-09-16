module io
#include "errorcodes.h"
#include "version.h"
   use parameters
   use parser
   implicit none
contains

   !**********************************************************************
   subroutine readConfigFileNew(inputfile)
      implicit none
      CHARACTER(len=*) inputfile
      CHARACTER(len=60) varname
      CHARACTER(len=256) line
      integer:: err
      logical:: aaSet, RaSet, RtSet, RcSet
      aaSet=.false.
      RaSet=.false.
      RtSet=.false.
      RcSet=.false.
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
               RtSet=.true.
            case('Rc')
               call read_val(line, Rc)
               RcSet=.true.
            case('aa')
               call read_val(line, alpha)
               aaSet=.true.
            case('Ra')
               call read_val(line, Ra)
               RaSet=.true.
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
      ! Validation
      if (aaSet.and.RaSet)then
         if((.not.RtSet).or.(.not.RcSet)) then
            Rt = Ra*cos(alpha)
            Rc = Ra*sin(alpha)
         endif
      endif
   end subroutine
   !**********************************************************************
   subroutine writeOutputHeader(unitOut)
      use parameters
      IMPLICIT none
      integer, intent(in):: unitOut
      WRITE(unitOut,*)  '### Output of Program glo  Ver.', VERSION,':   ###'
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
      WRITE(unitOut,'(A11,A3,A2)')    '# Variable par ',VariablePar , '         #'
      WRITE(unitOut,*)  '# see definition of LCALC for output. LCALC:', LCALC,'   #'
      WRITE(unitOut,*)  '#                                      #'
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

end module io
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
