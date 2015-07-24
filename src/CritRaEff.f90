!***********************************************************************
!
!***********************************************************************
#include "errorcodes.h"
#include "version.h"
program CriticalRaEff
   use parameters
   use GrowthRateRaEffMod
   use CritRaEff_io
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
   Print*, 'Out of init()'

   call computeCriticalCurve()
   
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

      call GrowthRateInit(Ra, alpha, Pt, Le, tau, eta, m0, Symmetry, Truncation)
      call setVariableParam(VariablePar)

      ! ----OUTPUT:
      OPEN(unitOut,FILE=outputfile,STATUS='UNKNOWN')
      call writeOutputHeader(unitOut)
   END subroutine

   !**********************************************************************
   !> Computes the lowest critical thermal Rayleigh number of all m's
   !! for all other parameters fixed.
   subroutine computeCriticalCurve()
      implicit none
      double precision, parameter:: DPI=3.141592653589793D0
      double precision:: CriticalRa, CriticalRaAlpha0
      double precision:: alpha, dalpha
      double precision:: RaMin, RaMax, gr1,gr2
      double complex:: frequency
      integer:: i
      integer:: info

      print*, 'In computeCriticalCurve()'
      info=0
      dalpha = 2.d0*dpi/3600.d0
      ! Start at alpha=0
      alpha = 0.0d0
      ! At this point a critical Ra is certain to exist so,
      ! increase the interval, until we find it.
      RaMin = 0
      RaMax = 10*Ra
      do
         gr1 = MaxGrowthRate(RaMin)
         gr2 = MaxGrowthRate(RaMax)
         Print*, 'RaMin=',RaMin,' -> gr1 = ', gr1, 'RaMax=', Ramax, ' -> gr2 = ', gr2
         if (gr1*gr2.gt.0.0d0) then
            RaMin = RaMax
            RaMax = 2*RaMax
	 else
            exit
         endif
      enddo
      ! Now that we found an interval find the critical value for Ra.
      call minimizer(MaxGrowthRate, RaMin, RaMax, RELE ,ABSE, NSMAX, CriticalRa, info)
      ! Cache this value for future use.
      CriticalRaAlpha0 = CriticalRa
      print*, 'First minimisation found CriticalRaAlpha0 =',  CriticalRa
      do i=1, 1800
         alpha = alpha + dalpha
         call GrowthRateUpdatePar(Ra=CriticalRa, alpha=alpha)
         RaMin = 0.5d0*CriticalRa
         RaMax = 1.5d0*CriticalRa
         call minimizer(MaxGrowthRate, RaMin, RaMax, RELE ,ABSE, NSMAX, CriticalRa, info)
         if (info.NE.0) exit
         frequency = MaxGrowthRateCmplx(CriticalRa)
         WRITE(*,*)        alpha, CriticalRa, dimag(frequency), dble(frequency)
         Write(unitOut, *) alpha, CriticalRa, dimag(frequency), dble(frequency)
      enddo
      alpha = 0.0d0
      CriticalRa = CriticalRaAlpha0
      do i=1, 1800
         alpha = alpha - dalpha
         call GrowthRateUpdatePar(Ra=CriticalRa, alpha=alpha)
         RaMin = 0.5d0*CriticalRa
         RaMax = 1.5d0*CriticalRa
         call minimizer(MaxGrowthRate, RaMin, RaMax, RELE ,ABSE, NSMAX, CriticalRa, info)
         if (info.NE.0) exit
         frequency = MaxGrowthRateCmplx(CriticalRa)
         WRITE(*,*)        alpha, CriticalRa, dimag(frequency), dble(frequency)
         Write(unitOut, *) alpha, CriticalRa, dimag(frequency), dble(frequency)
      enddo
   end subroutine
end program
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
