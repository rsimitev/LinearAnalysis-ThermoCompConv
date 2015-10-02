!***********************************************************************
!
!***********************************************************************
#include "errorcodes.h"
#include "version.h"
program CriticalRaEff
   use parameters
   use GrowthRateMod
   use CritRaEff_io
   use glo_constants
   implicit none
   type GlobalCrit
      double precision:: alpha
      double precision:: Ra
      double precision:: w
      integer:: m
   end type
   integer, parameter:: NN = 900  ! This is one point every 0.4 degrees.
   character*60:: infile,outfile
   character(len=2):: par
   integer, parameter:: unitOut=16
   double precision:: alphas(NN)
   type(GlobalCrit):: crit(NN)
   logical:: recompute=.false.
   logical:: restart=.false.

   !---------------------------------------------------------
   !  arg #1 - filename or usage ?
   call getarg(1,infile)
   if (trim(infile).eq.'' .or. infile.eq.'-h') then
      call usage()
      stop
   endif

   call getarg(2,outfile)
   if (trim(outfile).eq.'' .or. outfile.eq.'-h') then
      call usage()
      stop
   endif
   print*,  trim(infile),' - ',trim(outfile)

   call getarg(3,par)
   if (par.eq.'-h') then
      call usage()
      stop
   elseif(par.eq.'-r') then
      recompute = .true.
   elseif(par.eq.'-c') then
      restart = .true.
   endif
   
   call init(trim(infile))
   
   call createAlphas(alphas)
   select case(LCALC)
      case(0)
         ! Single m
         call computeCriticalCurveSingleM(alphas,crit)
      case(1)
         ! Global
         call computeCriticalCurve(alphas,crit)
   end select

contains

   !**********************************************************************
   !> Prints a usage message.
   subroutine usage()
      implicit none
      print*, 'Usage : '
      print*, 'CriticalRaEff -h'
      print*, '   Prints this message.'
      print*, ''
      print*, 'CriticalRaEff <in file> <out file> [-r|-c]'
      print*, '   Parameters must be used in this order.'
      print*, '   -r ignores any previous calculation and recomputes everything.'
      print*, '   -c tries to continue any previously started computation.'
      print*, '   default is to compute nothing if  any previous result is found.'
   end subroutine

   !**********************************************************************
   !> Initialises things.
   subroutine init(inputfile)
      implicit none
      CHARACTER(len=*), intent(in)::  inputfile

      ! ----Default values:
      call setDefaults()
      ! ----INPUT:
      call readConfigFileNew(inputfile)

      ! ---- doesn't work for M=0 !!!!!!
      IF(M0.LT.1) THEN
        write(*,*) 'The code does not work for M0<1. ', M0, ' --> 1'
        M0 = 1
      ENDIF

      call GrowthRateInitAlpha(Ra, alpha, Pt, Le, tau, eta, m0, Symmetry, Truncation)
      call setVariableParam('Ra ')

   end subroutine

   !**********************************************************************
   subroutine createAlphas(alphas)
      implicit none
      double precision, intent(out):: alphas(:)
      double precision:: dalpha
      integer:: n, i
      n=size(alphas,1)
      dalpha=2.0d0*dpi/n
      alphas(1)=-dpi

      do i=2, n
         alphas(i) = alphas(i-1) + dalpha
      enddo
   end subroutine

   !**********************************************************************
   !> Computes the lowest critical effective Rayleigh number as a function of alpha
   !! for all other parameters fixed.
   subroutine computeCriticalCurveM(alpha, crit)
      implicit none
      double precision, intent(in):: alpha(:)
      double precision, intent(out):: crit(:,:)
      double precision:: CriticalRa, CriticalRaAlpha0
      double precision:: RaMin, RaMax, gr1,gr2
      integer:: i, N, HalfN
      integer:: info, counter

      N     = size(alpha,1)
      HalfN = N/2
      Write(*,*) N, HalfN
      info  = 0
      crit(:,1) = huge(1.0d0)
      crit(:,2) = 0.0d0

      ! At alpha=0 a critical Ra is certain to exist so,
      ! change the interval, until we find it.
      RaMin = 0
      RaMax = 10*Ra
      call GrowthRateUpdateParAlpha(alpha=0.0d0)
      do
         gr1 = MaxGrowthRate(RaMin)
         gr2 = MaxGrowthRate(RaMax)
         if (gr1*gr2.gt.0.0d0) then
            RaMin = RaMax
            RaMax = 2*RaMax
         else
            exit
         endif
      enddo
      ! Now that we found an interval find the critical value for Ra.
      counter=0
      CriticalRa = (RaMin+RaMax)*0.5
      call computeSinglePoint(alpha, HalfN, CriticalRa, crit, counter)
      ! Cache this value for future use.
      CriticalRaAlpha0 = CriticalRa
      ! Compute the positive half of the alphas
      counter=0
      do i=HalfN+1, N
         call computeSinglePoint(alpha, i, CriticalRa, crit, counter)
         if (mod(i,5)==0 .and. counter.ne.3) call writeCriticalCurveSingleM(alpha, crit, m0)
      enddo

      ! and the negative half
      counter=0
      CriticalRa = CriticalRaAlpha0
      do i=HalfN-1, 1, -1
         call computeSinglePoint(alpha, i, CriticalRa, crit, counter)
         if (mod(i,5)==0 .and. counter.ne.3) call writeCriticalCurveSingleM(alpha, crit, m0)
      enddo
   end subroutine

   !**********************************************************************
   !>
   subroutine computeSinglePoint(alpha, i, CriticalRa, crit, counter)
      implicit none
      double precision, intent(in):: alpha(:)
      integer, intent(in):: i
      integer, intent(inout):: counter
      double precision, intent(inout):: CriticalRa, crit(:,:)
      double precision:: RaMax, Ramin
      integer:: info
      info=0
      
      if (counter.ge.3) then
         crit(i,1) = huge(1.0d0)
         crit(i,2) = 0.0d0
         return
      else
         ! if we are to restart the calclulation, use any previous result
         ! that is not a huge number and return.
         if (restart.and.crit(i,1).lt.huge(1.0d0)) return
         ! Otherwhise, do the actual leg work.
         call GrowthRateUpdateParAlpha(Ra=CriticalRa, alpha=alpha(i))
         RaMin = 0.05d0*CriticalRa
         RaMax = 100.d0*CriticalRa
         call minimizer(MaxGrowthRate, RaMin, RaMax, RELE ,ABSE, NSMAX, CriticalRa, info)
         if (info.NE.0) then
            if (counter.ne.3) then
               counter = counter + 1
               CriticalRa = (RaMin+RaMax)/2.0
               crit(i,1) = huge(1.0d0)
               crit(i,2) = 0.0d0
               return
            endif
         else
            counter = 0
         endif
         call logSinglePoint(alpha(i), CriticalRa)
         crit(i,1) = CriticalRa
         crit(i,2) = dble(MaxGrowthRateCmplx(CriticalRa))
      endif
   end subroutine
   
   !**********************************************************************
   !>
   subroutine logSinglePoint(alpha, CriticalRa)
      implicit none
      double precision, intent(in):: alpha, CriticalRa
      character(8)  :: date
      character(10) :: time

      call date_and_time(DATE=date,TIME=time)
      Write(*,*) '[',date,'-',time,']', ' alpha = ', alpha, CriticalRa
   end subroutine
   
   !**********************************************************************
   !>
   subroutine computeCriticalCurve(alpha, crit)
      implicit none
      double precision, intent(in):: alpha(:)
      type(GlobalCrit), intent(out):: crit(:)
      double precision, allocatable:: crit_new(:,:)
      integer:: m, N, i

      N = size(alpha,1)
      allocate(crit_new(N,2))
      crit_new(:,1) = huge(1.0d0)
      crit_new(:,2) = 0.0d0
      do i=1, N
         crit(i)%alpha = alpha(i)
         crit(i)%m  = Huge(1)
         crit(i)%w  = Huge(1.0d0)
         crit(i)%Ra = Huge(1.0d0)
      enddo

      do m=1, m0
         call GrowthRateUpdateParAlpha(m=m)
         Write(*,*) 'Computing critical value for m =', m
         if(wasPreviouslyComputed(m)) then
            call readCriticalCurveSingleM(alpha,crit_new,m)
            if (recompute.or.restart) then
               call computeCriticalCurveM(alpha, crit_new)
               call writeCriticalCurveSingleM(alpha, crit_new, m)
            endif
         else
            call computeCriticalCurveM(alpha, crit_new)
            call writeCriticalCurveSingleM(alpha, crit_new, m)
         endif
         do i=1,N
            if(crit_new(i,1).lt.crit(i)%Ra) then
               crit(i)%Ra = crit_new(i,1)
               crit(i)%w  = crit_new(i,2)
               crit(i)%m  = m
            endif
         enddo
         call writeCriticalCurve(crit)
      enddo
      deallocate(crit_new)
   end subroutine

   !**********************************************************************
   !>
   subroutine computeCriticalCurveSingleM(alpha, crit)
      implicit none
      double precision, intent(in):: alpha(:)
      type(GlobalCrit), intent(out):: crit(:)
      double precision, allocatable:: crit_new(:,:)
      integer:: N, i
      integer:: info

      N     = size(alpha,1)
      allocate(crit_new(N,2))
      info  = 0
      crit_new(:,1) = huge(1.0d0)
      crit_new(:,2) = 0.0d0
      call GrowthRateUpdateParAlpha(m=m0)
      Write(*,*) 'Computing critical value for m =', m0
      if(wasPreviouslyComputed(m0)) then
         call readCriticalCurveSingleM(alpha,crit_new,m0)
         if (recompute) then
            call computeCriticalCurveM(alpha, crit_new)
            call writeCriticalCurveSingleM(alpha, crit_new, m0)
         endif
      else
         call computeCriticalCurveM(alpha, crit_new)
         call writeCriticalCurveSingleM(alpha, crit_new, m0)
      endif
      do i=1,N
         crit(i)%alpha = alpha(i)
         crit(i)%Ra = crit_new(i,1)
         crit(i)%w  = crit_new(i,2)
         crit(i)%m  = m0
      enddo
      deallocate(crit_new)
   end subroutine

   !**********************************************************************
   !>
   subroutine writeCriticalCurve(crit)
      implicit none
      type(GlobalCrit), intent(in):: crit(:)
      integer:: N, i
      ! ----OUTPUT:
      OPEN(unitOut,FILE=outfile,STATUS='UNKNOWN')
      call writeOutputHeader(unitOut)
      N = size(crit,1)
      do i=1, N
         if (crit(i)%Ra.ge.huge(1.0d0)) then
            Write(unitOut,*) crit(i)%alpha, ' NaN ', ' NaN ' , ' NaN'
         else
            Write(unitOut,*) crit(i)%alpha, crit(i)%Ra, crit(i)%m, crit(i)%w  
         endif
      enddo
      close(unitOut)
   end subroutine

   !**********************************************************************
   !>
   subroutine writeCriticalCurveSingleM(alpha, crit, m)
      implicit none
      double precision, intent(in):: alpha(:)
      double precision, intent(in):: crit(:,:)
      integer, intent(in):: m
      character(len=3):: num
      integer:: N, i
      integer, parameter:: unitm=999
      N = size(crit,1)
      Write(num,'(I3.3)') m
      open(unit=unitm,file=trim(outfile)//'.'//trim(num), status='UNKNOWN')
      do i=1, N
         Write(unitm,*) alpha(i), crit(i,1), crit(i,2)
      enddo
      close(unitm)
   end subroutine
   
   !**********************************************************************
   !>
   subroutine readCriticalCurveSingleM(alpha, crit, m)
      implicit none
      double precision, intent(out):: crit(:,:)
      double precision, intent(in):: alpha(:)
      double precision, allocatable:: aa2(:), crit2(:,:), trans(:,:)
      integer, intent(in):: m
      double precision:: a1, a2, GR, OM
      character(len=3):: num
      integer:: N, N2, i
      integer, parameter:: unitm=999

      N = size(crit,1)
      Write(num,'(I3.3)') m
      open(unit=unitm,file=trim(outfile)//'.'//trim(num), status='OLD')
      ! read the first entry and check that the dalpha is the same
      read(unitm,*) a1, GR, OM
      read(unitm,*) a2, GR, OM
      if (a2.ne.alpha(2)) then
         N2 = int(2*dpi/(a2-a1))
      else
         N2 = N
      endif
      rewind(unitm)
      ! Allocate buffer arrays with the appropriate size.
      allocate(aa2(N2))
      allocate(crit2(N2,2))
      ! Store the first entry.
      ! Read the rest
      do i=1, N2
         read(unitm,*) aa2(i), crit2(i,1), crit2(i,2)
      enddo
      close(unitm)
      ! Interpolate or just copy if sizes are the same.
      if (N2.eq.N) then
         crit(:,1) = crit2(:,1)
         crit(:,2) = crit2(:,2)
      else
         allocate(trans(N,N2))
         call createTrans(alpha,aa2,trans)
         crit(:,1) = matmul(trans(:,:),crit2(:,1))
         crit(:,2) = matmul(trans(:,:),crit2(:,2))
         deallocate(trans)
      endif
      deallocate(aa2,crit2)
   end subroutine
   
   !**********************************************************************
   !> Creates the transformation matrix that converts from 
   !! one resolution to the other. At the moment a simple linear 
   !! interpolation between nearest neighbours is used.
   subroutine createTrans(aaOut,aaIn,trans)
      implicit none
      double precision, intent(in):: aaOut(:), aaIn(:)
      double precision, intent(out):: trans(:,:)
      integer:: N, M, i, j

      N = size(aaOut,1)
      M = size(aaIn,1)
      trans = 0.0d0
      ! Take care of the initial point
      if (aaOut(1).lt.aaIn(1)) then
         trans(1,1) = 1.0d0
      else
         trans(1,1)   =  (aaOut(1) - aaIn(1+1))/(aaIn(1)-aaIn(1+1))
         trans(1,1+1) = -(aaOut(1) - aaIn(1)  )/(aaIn(1)-aaIn(1+1))
      endif
      ! Use a linear interpolation where necessary.
      do i=2, N-1
         do j=1, M-1
            if( (aaOut(i).ge.aaIn(j)) .and. (aaOut(i).lt.aaIn(j+1)) ) then
               trans(i,j)   =  (aaOut(i) - aaIn(j+1))/(aaIn(j)-aaIn(j+1))
               trans(i,j+1) = -(aaOut(i) - aaIn(j)  )/(aaIn(j)-aaIn(j+1))
            endif
         enddo
      enddo
      ! take cae of the last point
      if (aaOut(N).gt.aaIn(M)) then
         trans(1,1) = 1.0d0
      else
         trans(N,M-1) =  (aaOut(N) - aaIn(M)   )/(aaIn(M-1)-aaIn(M))
         trans(N,M  ) = -(aaOut(N) - aaIn(M-1) )/(aaIn(M-1)-aaIn(M))
      endif
   end subroutine
   
   !**********************************************************************
   !> 
   logical function wasPreviouslyComputed(m)
      implicit none
      integer, intent(in):: m
      character(len=3):: num
      Write(num,'(I3.3)') m
      inquire(file=trim(outfile)//'.'//trim(num), EXIST=wasPreviouslyComputed)
   end function
end program
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
