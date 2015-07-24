!***********************************************************************
!
!***********************************************************************
#include "errorcodes.h"
#include "version.h"
program CriticalRaEff
   use parameters
   use growthRateMod
   use io
   implicit none
   character*60 infile,outfile
   integer, parameter:: unitOut=16

!---------------------------------------------------------
!  arg #1 - filename or usage ?
   call getarg(1,infile)
   if (infile.eq.' ') then
      print*, 'Usage : '
      print*, 'CriticalRaEff <in file> <out file>'
      stop
   endif
   if (infile.eq.'-h') then
      print*, 'Usage : '
      print*, 'CriticalRaEff <in file> <out file>'
      stop
   endif

   call getarg(2,outfile)
   print*,  trim(infile),' - ',trim(outfile)

   call init(trim(infile),trim(outfile))


   close(unitOut)
contains

   !**********************************************************************
   !> Initialises things.
   SUBROUTINE init(inputfile,outputfile)
      implicit none
      CHARACTER(len=*) inputfile,outputfile

      ! ----Default values:
      call setDefaults()
      ! ----INPUT:
      call readConfigFileNew(inputfile)

      ! ---- doesn't work for M=0 !!!!!!
      IF(M0.LT.1) THEN
        write(*,*) 'The code does not work for M0<1. ', M0, ' --> 1'
        M0 = 1
      ENDIF

      call GrowthRateInit(Rt, Rc, Pt, Le, tau, eta, m0, Symmetry, Truncation)
      call setVariableParam(VariablePar)

      ! ----OUTPUT:
      OPEN(unitOut,FILE=outputfile,STATUS='UNKNOWN')
      call writeOutputHeader(unitOut)
   END subroutine

   !**********************************************************************
   !> Computes the lowest critical thermal Rayleigh number of all m's
   !! for all other parameters fixed.
   subroutine fixedParCriticalRaAndM0()
      implicit none
      double precision, parameter:: DPI=3.141592653589793D0
      double precision:: CriticalRa
      double precision:: Omega, alpha
      double precision:: RaMin, RaMax, gr1,gr2
      double complex:: frequency
      integer:: i, ii
      integer:: info, idx

      DO m0=int(LowerLimit), int(UpperLimit)
         info=0
         call GrowthRateUpdatePar(m=m0)
         alpha = 0.0d0
         ! Increase the interval, until we find a critical Ra
         RaMin = 0
         RtMax = 10*CriticalRa
         do
            gr1 = MaxGrowthRate(RaMin)
            gr2 = MaxGrowthRate(RaMax)
            if (gr1*gr2.gt.0.0d0) then
               RaMin = RaMax
               RaMax = 10*RaMin
            endif
         enddo
         ! Now that we found an interval find the critical value for Ra.
         call minimizer(MaxGrowthRate, RaMin, RaMax, RELE ,ABSE, NSMAX, CriticalRa, info)
         ! Cache this value for future use.
         CriticalRaAlpha0 = CriticalRa
         do
            call minimizer(MaxGrowthRate, RtMin, RtMax, RELE ,ABSE, NSMAX, CriticalRt, info)
            if (info.NE.2) exit
         enddo
         if(info.NE.0) then
            Write(*,*) 'Failed to find roots: error:', info
            stop 5
         endif
         ! Compute the grouth rate for each case so we can output it for
         ! debug purposes.
         frequency = MaxGrowthRateCmplx(CriticalRt)
         GRORi(ii) = dimag(frequency)
         OMI(II) = -dble(frequency)
         RACI(II)= CriticalRt
         LLI(II) = info
         write( *,'(1X,1P,4E17.6,I4,A3,2E17.6)') TAU, CriticalRt, OMI(II), GRORi(ii), M0,' | ', CriticalRt*(1/(1-eta))**4*(2/TAU),(OMI(II)*2.0/TAU)
      enddo
end program
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
