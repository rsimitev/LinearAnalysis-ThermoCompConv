!***********************************************************************
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
! Rt=RAYLEIGH NUMBER, TA= TAYLOR NUMBER, Pt= PRANTEL NUMBER,
! ETA=RATIO OF RADII, Truncation= TRUNCATION, M0=WAVE NUMBER,
! Le = Lewis number,  CriticalRt=Rayleigh number due to concentration
! Symmetry= SYMMETRIE PARAMETER, Symmetry=0 : UNDEFINED SYMMETRIE,
! Symmetry=2 : EQUATORIAL SYMMETRIE, Symmetry=1 : EQUATORIAL ANTISYMMETRIE.
! DRIFT C IS DEFINED LIKE (PHI+C*T).
! Rev. 2.2: Drift is now def. as (phi-c*t).
!
! LCALC=1 : Eigenvalues are determined for const. parameters
! LCALC=2 : Onset determined for constant wavenumber M
!           (by searching root of grothrate in R, using pegasus.f).
! LCALC=3 : Onset determined by variing Rayleighnumber R and
!           wavenumber M.
! LCALC=4 : Eigenvector determined for one set of parameters
!           at onset
!
! LO.F calculates R (crit. Rayleighn.) and Omega (and M) in the
! range TAU=LowerLimit to TAU=UpperLimit.
! It only calculates the mode with minimal value of R.
!
! Start:    glo inputfilename outputfilename
!
#include "losub-inc.h"
program linearOnset
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
      print*, 'glo <in file> <out file>'
      stop
   endif
   if (infile.eq.'-h') then
      print*, 'Usage : '
      print*, 'glo <in file> <out file>'
      stop
   endif

   call getarg(2,outfile)
   print*,  trim(infile),' - ',trim(outfile)

   call init(trim(infile),trim(outfile))

   select case(LCALC)
      case(-1)! most basic case: find the most unstable growth rate at all parameters fixed
         CALL fixedParGrowthRate()
      case(0)! Critical Rt, for constant other parameters
         CALL fixedParCriticalRa()
      case(1) ! eigenvalues determined for all parameters fixed.
         call fixedParEigenValues()
      case(2)
         LowerLimit = tau
         call varyTauCriticalRt()
      case(3)
         LowerLimit = tau
         call varyTauCriticalState()
      case(4) ! calculate the critical eigenvector. Print for plotting.
         CALL fixedParCriticalEigenVector()
      case(5) ! vary m and calculate critical R at fixed P, tau, eta.
         LowerLimit = m0
         call varyMCriticalRt()
      case(6) ! vary Le and calculate critical R at fixed P, tau, eta, M
         LowerLimit = Le
         call varyLeCriticalRt()
      case default
         Write(*,*) 'Unknown computation type:', LCALC
   end select
   close(unitOut)
contains

   !**********************************************************************
   SUBROUTINE init(inputfile,outputfile)
      implicit none
      CHARACTER(len=*) inputfile,outputfile

      ! ----Default values:
      CALL setDefaults()
      ! ----INPUT:
      CALL readConfigFile(inputfile)

      ! ---- doesn't work for M=0 !!!!!!
      IF(M0.LT.1) THEN
        write(*,*) 'The code does not work for M0<1. ', M0, ' --> 1'
        M0 = 1
      ENDIF

      ! ----OUTPUT:
      OPEN(unitOut,FILE=outputfile,STATUS='UNKNOWN')
      CALL writeOutputHeader(unitOut)

      RI = ETA/(1.0d0-ETA)
      RO = 1.0D0 + RI

      IF(Symmetry.EQ.0) THEN
         ! - UNDEFINED SYMMETRIE:
         LMIN=M0
         LD=1
      ELSEIF(Symmetry.EQ.1) THEN
         ! - EQUATORIAL ANTISYMMETRIE (L+M ODD):
         LMIN=M0+1
         LD=2
      ELSEIF(Symmetry.EQ.2) THEN
         ! - EQUATORIAL SYMMETRIE (L+M EVEN);
         LMIN=M0
         LD=2
      ENDIF

      CALL DIMENSION(LMIN,LD,Truncation,M0,NEigenmodes)
      write(*,*) 'DIMENSION OF MATRIX:',NEigenmodes
   END subroutine

   !**********************************************************************
   subroutine fixedParEigenValues()
      implicit none
      double complex:: ZEW(NMAX)
      integer:: i
      call computeGrowthRateModes(.true., zew)
      write(*,*) 'n     Frequ.(exp(+iwt))   -Grothrate  '
      write(unitOut,*) 'n     Frequ.(exp(+iwt))   -Grothrate  '
      do i=1, NEigenmodes
         WRITE(*,'(I4,2D16.6)') I,ZEW(I)
         WRITE(unitOut,'(I4,2D16.6)') I,ZEW(I)
      enddo
   end subroutine

!**********************************************************************
   subroutine fixedParGrowthRate()
      implicit none
      INTEGER:: aux
      double precision:: Rt0, GroR, factor

      Rt0=Rt

      if (Rt.gt.UpperLimit) then
         factor = -1.0d0
      else
         factor = 1.0d0
      endif
      do while(abs(Rt-UpperLimit) > abs(0.11*UpperLimit))
         aux = int(log10(abs(Rt)))
         if(aux.ge.2) then
            StepSize = factor*10.0d0**(aux-1)
         else
            StepSize = factor*10.0d0
         endif
         Rt   = Rt + StepSize
         GROR = MaxGrowthRate(Rt)
         WRITE(*,*) Rt, GROR
         Write(unitOut, *) Rt, GroR
      enddo

      WRITE(*,*) 'R=',Rt,' TAU=',TAU,' P=',Pt,' M0=',M0,' eta=',ETA
      WRITE(*,*) 'Most unstable growth rate', GROR
      WRITE(*,*) 'If growth rate < 0 then above onset'
   end subroutine

!**********************************************************************
   subroutine fixedParCriticalRa()
      implicit none
      double Complex:: frequency
      double precision:: CriticalRt, Omega, GroR
      integer:: info
      CALL minimizer(MaxGrowthRate,Rt/10, Rt*10,RELE,ABSE,NSMAX,CriticalRt, info)
      frequency = MaxGrowthRateCmplx(CriticalRt)
      GROR  = dimag(frequency)
      OMEGA = -dble(frequency)
      Write(unitOut,*) CriticalRt, Omega, GroR
      WRITE(*,*) 'TAU=',TAU,' Pt=',Pt,' M0=',M0,' eta=',ETA
      WRITE(*,*) 'Lewis=',Le,' Rc=',Rc
      WRITE(*,*) 'R_crit=',CriticalRt, '  (growth rate =',GROR,')'
   end subroutine


!**********************************************************************
      subroutine fixedParCriticalRaAndM0(CriticalRt)
      implicit none
      double precision, intent(out):: CriticalRt
      integer, parameter::nm0=3
      double precision:: Omega
      double precision:: RACI(nm0),OMI(nm0),GRORi(nm0)
      double precision:: RtMin, RtMax
      double Complex:: frequency
      INTEGER:: M0I(nm0),LLI(nm0),LMINI(nm0)
      INTEGER, save:: NTRYCOUNT=0
      integer:: first_m0, first_lmin, i, ii
      integer:: info, idx

      info=0
      CriticalRt = Rt
      do i=0, nm0/2
         first_m0   = m0 - i
         first_lmin = lmin-i
         if (first_m0==1) exit
      enddo
!     Start a search around the previous m0
      do i=1, nm0
        M0I(i)   = first_m0+i-1
        LMINI(i) = first_lmin+i-1
      enddo
      DO II=1, nm0
         M0   = M0I(II)
         LMIN = LMINI(II)
         CALL dimension(LMIN,LD,Truncation,M0,NEigenmodes)
         ! Increase the interval, in case we did not find anything.
         do i=1, 3
            RtMin=Rt/(2.0d0**dble(i))
            RtMax=Rt*(2.0d0**dble(i))
            ! Write(*,*)  i, Rt, RtMin, RtMax, RELE ,ABSE, NSMAX
            CALL minimizer(MaxGrowthRate, RtMin, RtMax, RELE ,ABSE, NSMAX, CriticalRt, info)
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

      ! Let's not allow growth rates that are too high
      if(any(lli.eq.0)) then
         where(lli.NE.0) RACI = 1.0d100
         idx = minloc(RACI,1)
      else
         idx = 0
      endif

      IF( idx.GT.0 ) THEN
         CriticalRt   = RACI(idx)
         OMEGA = OMI(idx)
         M0    = M0I(idx)
         LMIN  = LMINI(idx)
         write( unitOut,'(1P,3E17.6,I4)') TAU, CriticalRt, OMEGA, M0
         write( *,'(">",1P,3E17.6,I4,A3,2E17.6)') TAU,CriticalRt,OMEGA, M0,' | ', CriticalRt*(1/(1-eta))**4*(2/TAU),(OMEGA*2.0/TAU)
         write(*,*)
         NTRYCOUNT = 0
      ELSE IF(NTRYCOUNT.GE.3) THEN
         WRITE(unitOut,*) 'NO CRITICAL RAYLEIGH NUMBER FOUND.'
         STOP NO_RA_FOUND
      ELSE
         write(*,*) 'NO CRIT. RAYLEIGH NUMBER FOUND. Trying again.'
         NTRYCOUNT = NTRYCOUNT + 1
      ENDIF
   end subroutine

!**********************************************************************
   subroutine fixedParCriticalEigenVector()
      implicit none
      double precision:: GroR, Omega, Ta
      complex(8):: ZEVEC(NMAX), zew(NMAX), ZEVAL(NMAX,NMAX)
      integer:: info, i, ni, li, lti, lpi
      integer:: NTH, KTV, KTH, LTV, LTH, lst, NUC
      double precision:: GRR, GRI, Pm, C0,CriticalRt, OMM
      integer:: nuds, nuom, MF, LMAX, NIMAX

      TA = TAU*TAU
      ! find the zero grothrate:
      CriticalRt = Rt
      CALL minimizer(MaxGrowthRate, Rt/10, Rt*10, RELE ,ABSE, NSMAX, CriticalRt, info)
      Rt = CriticalRt
      call computeGrowthRateModes(.TRUE.,zew,zeval)
      GROR  = DIMAG(zew(1))
      Omega = -dble(zew(1))
      zevec(:) = zeval(:,1)

      IF(info.EQ.0) THEN
         ! print eigenvector
         ! Fileformat for outputfile:
         ! LST=0 formatted Wicht
         ! LST=1 formatted Hirsching
         ! LST=3 unformatted
          LST=0
          WRITE(unitOut,'(2I2,'' LINEAR ONSET '')')  LST,LCALC
          NTH=0
          KTV=0
          KTH=0
          LTV=0
          LTH=0
          GRR=0.D0
          GRI=0.D0
          WRITE(unitOut,'(I2,7I3,2D16.8,'' M0,TRUNC,LD,GROTH,DRIFT'')')  M0,Truncation,NTH,KTV,KTH,LTV,LTH,LD,GRR,GRI
          NUDS=1
          PM=0.D0
          WRITE(unitOut, '(I5,2D14.6,D9.2,D13.6,D9.2,'' I,TA,Rt,Pt,PM,E'')') NUDS,TA,CriticalRt,Pt,PM,ETA
          C0 = OMEGA/M0
          OMM=0.D0
          NUC=0
          NUOM=0
          MF=0
          WRITE(unitOut,9100) C0,OMM,NUC,NUOM,MF
9100      FORMAT(2D17.10,3I4,'    C,OM, WHERE?,FLOQUET')

          LMAX=2*Truncation+M0-1
          I=0
          DO 3000 LI=LMIN,LMAX,LD
!           L for poloidal (v) field:
            LPI = LI
            NIMAX=INT( DBLE(2*Truncation+1-LI+M0)/2 )
            DO 3000 NI=1,NIMAX
              IF( LST.EQ.0 ) THEN
               WRITE(unitOut,9200) 'V',LPI,M0,NI,0, DREAL(ZEVEC(I+1)),DIMAG(ZEVEC(I+1)),0.D0,0.D0
              ELSEIF(LST.EQ.3) THEN
               WRITE(unitOut,'(A,4I3,A,2F11.7,A)') ' ''V ''',LPI,M0,NI,0,' ', DREAL(ZEVEC(I+1)),DIMAG(ZEVEC(I+1)),' .0D+00 .0D+00 '
              ENDIF
             I=I+4
3000      CONTINUE
!
          I=0
          DO 3200 LI=LMIN,LMAX,LD
!           L for toroidal (w) field:
            IF( Symmetry.EQ.2 ) THEN
             LTI=LI+1
            ELSEIF( Symmetry.EQ.1 ) THEN
             LTI=LI-1
            ELSEIF( Symmetry.EQ.0 ) THEN
             LTI=LI
            ENDIF
            NIMAX=INT( DBLE(2*Truncation+1-LI+M0)/2 )
            DO 3200 NI=1,NIMAX
              IF( LST.EQ.0 ) THEN
               WRITE(unitOut,9200) 'W',LTI,M0,NI,0, DREAL(ZEVEC(I+3)),DIMAG(ZEVEC(I+3)),0.D0,0.D0
              ELSEIF(LST.EQ.3) THEN
               WRITE(unitOut,'(A,4I3,A,2F11.7,A)') ' ''W ''',LTI,M0,NI,0,' ', DREAL(ZEVEC(I+3)),DIMAG(ZEVEC(I+3)),' .0D+00 .0D+00 '
              ENDIF
              I=I+4
3200      CONTINUE
!
          I=0
          DO 3400 LI=LMIN,LMAX,LD
            NIMAX=INT( DBLE(2*Truncation+1-LI+M0)/2 )
            DO 3400 NI=1,NIMAX
              IF( LST.EQ.0 ) THEN
               WRITE(unitOut,9200) 'T',LI,M0,NI,0, DBLE(ZEVEC(I+2)),DIMAG(ZEVEC(I+2)),0.D0,0.D0
              ELSEIF(LST.EQ.3) THEN
               WRITE(unitOut,'(A,4I3,A,2F11.7,A)') ' ''T ''',LI,M0,NI,0,' ', DBLE(ZEVEC(I+2)),DIMAG(ZEVEC(I+2)),' .0D+00 .0D+00 '
              ENDIF
             I=I+4
3400      CONTINUE
!
          I=0
          DO 3600 LI=LMIN,LMAX,LD
            NIMAX=INT( DBLE(2*Truncation+1-LI+M0)/2 )
            DO 3600 NI=1,NIMAX
              IF( LST.EQ.0 ) THEN
               WRITE(unitOut,9200) 'G',LI,M0,NI,0, DREAL(ZEVEC(I+4)),DIMAG(ZEVEC(I+4)),0.D0,0.D0
              ELSEIF(LST.EQ.3) THEN
               WRITE(unitOut,'(A,4I3,A,2F11.7,A)') ' ''G ''',LI,M0,NI,0,' ', DREAL(ZEVEC(I+4)),DIMAG(ZEVEC(I+4)),' .0D+00 .0D+00 '
              ENDIF
             I=I+4
3600      CONTINUE
!
9200            FORMAT(1X,A1,4I3,4D16.8)
         ELSE
          WRITE(unitOut,*) 'NO CRITICAL RAYLEIGH NUMBER FOUND.'
          STOP NO_RA_FOUND
       ENDIF
      end subroutine

!**********************************************************************
   subroutine varyTauCriticalRt()
      implicit none
      integer:: NTRYCOUNT, info
      double complex:: frequency
      double precision:: Ta, CriticalRt, GroR, RtOld
      double precision:: Tau0, Tau1, omega, RtMin, RtMax

      NTRYCOUNT = 0
      LowerLimit = Tau
      TAU0 = LowerLimit
      TAU1 = LowerLimit
      TA = TAU*TAU
      do
         TA = TAU**2
         ! searching for zero grothrate by varying Rt: ----------------------
         RtMin = Rt/5.0
         RtMax = Rt*5
         RtOld = Rt
         CALL minimizer(MaxGrowthRate, RtMin, RtMax, RELE, ABSE, 50, CriticalRt, info)
         frequency = MaxGrowthRateCmplx(CriticalRt)
         GROR  = dimag(frequency)
         OMEGA = -dble(frequency)
         IF(info.EQ.0) THEN
            WRITE(unitOut,'(1P,3E17.6,I4)') TAU,CriticalRt,OMEGA,M0
            WRITE( *,'(1P,3E17.6,I4)') TAU,CriticalRt,OMEGA,M0
            NTRYCOUNT = 0
         ELSE
            write(*,*) '# NO CRIT. RAYLEIGH NUMBER FOUND. Trying again.'
            NTRYCOUNT = NTRYCOUNT + 1
         ENDIF
!--      increment TAU:
         TAU0 = TAU1
         TAU1 = TAU
         IF(NTRYCOUNT==0) THEN
            IF( DABS(StepSize) .LT. TAU*0.1D0 ) THEN
               TAU = TAU1 + StepSize
            ELSE
               TAU = TAU1 + 10.0d0**(int(log10(abs(tau1)))-1)*StepSize/dabs(StepSize)
            ENDIF
!--         interpolate new startingvalue for Rt:
!           Assume CriticalRt grows with tau**(4/3)
            Rt = CriticalRt + (4.0/3.0)*Tau1**(1.0/3.0)*(TAU-TAU1)
            RtOld = CriticalRt
         ELSE
            Rt = CriticalRt
         ENDIF

         if (NTRYCOUNT==3) exit
         ! end value of TAU reached?
         if(UpperLimit.GT.LowerLimit) then
            if (TAU.GT.UpperLimit) exit
         else
            if (TAU.LT.UpperLimit) exit
         endif
      enddo !--         End of tau Loop
   end subroutine

   subroutine varyTauCriticalState()
      implicit none
      double precision:: CriticalRt, RtOld
      double precision:: Tau0, Tau1
      Rtold = Rt
      LowerLimit = Tau
      TAU0 = LowerLimit
      TAU1 = LowerLimit
      do
!--------searching for zero grothrate by varying Rt and M0: ---------------
         CALL fixedParCriticalRaAndM0(CriticalRt)
!--      increment TAU:
         TAU0 = TAU1
         TAU1 = TAU
         IF( DABS(StepSize) .LT. DABS(TAU*0.1D0) ) THEN
            TAU = TAU1 + StepSize
         ELSE
            TAU = TAU1 + 10.0d0**(int(log10(abs(tau1)))-1)*StepSize/dabs(StepSize)
         ENDIF
!--      interpolate new startingvalue for Rt:
         IF(dabs(TAU1-TAU0).le.1.0d-10) THEN
!          Assume Rac grows with Ta**(4/3)
           Rt = CriticalRt + (4.0/3.0)*Tau1**(1.0/3.0)*(TAU-TAU1)
         ELSE
           Rt = 2.0D0*CriticalRt - RtOld
         ENDIF
         RtOld = CriticalRt

         ! end value of TAU reached?
         if(UpperLimit.GT.LowerLimit) then
            if (TAU.GT.UpperLimit) exit
         else
            if (TAU.LT.UpperLimit) exit
         endif
      enddo
   end subroutine

   subroutine varyMCriticalRt()
      implicit none
      double precision:: CriticalRt
      integer:: info
      CriticalRt=Rt
      do M0=int(LowerLimit), int(UpperLimit), INT(StepSize)
         IF(Symmetry.EQ.0) THEN
! -   UNDEFINED SYMMETRIE:
            LMIN=M0
            LD=1
         ELSEIF(Symmetry.EQ.1) THEN
! -   EQUATORIAL ANTISYMMETRIE (L+M ODD):
            LMIN=M0+1
            LD=2
         ELSEIF(Symmetry.EQ.2) THEN
! -   EQUATORIAL SYMMETRIE (L+M EVEN);
            LMIN=M0
            LD=2
         ENDIF
         CALL dimension(LMIN,LD,Truncation,M0,NEigenmodes)
         CALL minimizer(MaxGrowthRate,Rt/10, Rt*10,RELE,ABSE,NSMAX,CriticalRt, info)
         WRITE(*,*) M0,CriticalRt
         WRITE(unitOut,*) M0, CriticalRt
      enddo
   end subroutine

   !> Varies the Lewis number and computes the critical
   !! thermal Rayleigh number for fixed other parameters.
   subroutine varyLeCriticalRt()
      implicit none
      double precision:: LeOld, GroR
      double precision:: CriticalRt
      integer:: niter, i, info
      LeOld = Le
      niter = int((UpperLimit-LeOld)/StepSize)
      do i=0, niter
         Le  = LeOld + i*StepSize
         CALL minimizer(MaxGrowthRate,Rt/10, Rt*10,RELE,ABSE,NSMAX,CriticalRt, info)
         GROR = MaxGrowthRate(CriticalRt)
         WRITE(*,*) Le, CriticalRt, GROR
         WRITE(unitOut,'(3D16.8)') Le, CriticalRt, GROR
      enddo
   end subroutine


end program
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
