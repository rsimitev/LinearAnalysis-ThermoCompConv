!***********************************************************************
! Program to calc. lin. onset of conv. in rotating spherical shell.
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
program LinearModePlotter
   use parameters
   use growthRateMod
   use io
   implicit none
   character*60 infile,outfilebase
   integer, parameter:: unitFlowPolOut=16
   integer, parameter:: unitFlowTorOut=17
   integer, parameter:: unitFlowTempOut=18
   integer, parameter:: unitFlowCompOut=19

!---------------------------------------------------------
!  arg #1 - filename or usage ?
   call getarg(1,infile)
   if (infile.eq.' '.or.infile.eq.'-h') then
      print*, 'Usage : '
      print*, 'modePlotter <in file> <out files base name>'
      stop
   endif

   call getarg(2,outfilebase)
   print*,  trim(infile),' - ',trim(outfilebase)

   call init(trim(infile),trim(outfilebase))

   CALL fixedParCriticalEigenVector()

   close(unitOut)
contains

   !**********************************************************************
   !> Initialises things.
   SUBROUTINE init(inputfile,outputfilebase)
      implicit none
      CHARACTER(len=*) inputfile,outputfilebase

      ! ----Default values:
      CALL setDefaults()
      ! ----INPUT:
      CALL readConfigFile(inputfile)

      ! ----OUTPUT:
      OPEN(unitFlowPolOut,FILE=outputfile//'.pol',STATUS='UNKNOWN')
      OPEN(unitFlowTorOut,FILE=outputfile//'.tor',STATUS='UNKNOWN')
      OPEN(unitTempOut,FILE=outputfile//'.temp',STATUS='UNKNOWN')
      OPEN(unitCompOut,FILE=outputfile//'.thetac',STATUS='UNKNOWN')
      OPEN(unitParOut,FILE=outputfile//'.par',STATUS='UNKNOWN')

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

   END subroutine

   !**********************************************************************
   !> Computes and Writes out the eigen vector (equatorial render)
   !! Corresponding to the critical value of Rt with all other parameters 
   !! fixed.
   subroutine fixedParCriticalEigenVector()
      implicit none
      double precision:: GroR, Omega
      complex(8), allocatable:: ZEVEC(:), zew(:), ZEVAL(:,:)
      integer:: info, i, ni, li, lti, lpi
      integer:: NTH, KTV, KTH, LTV, LTH, lst, NUC
      double precision:: GRR, GRI, Pm, C0, OMM
      integer:: nuds, nuom, MF, LMAX, NIMAX

      ! find the zero grothrate:
      call fixedParCriticalParAndM0_v2()
      Nmodes = getEigenProblemSize()
      allocate( ZEVEC(Nmodes), zew(Nmodes), ZEVAL(Nmodes,Nmodes))
      call computeGrowthRateModes(.TRUE.,zew,zeval)
      GROR  = -DIMAG(zew(1))
      Omega = dble(zew(1))
      zevec(:) = zeval(:,1)

      ! print eigenvector

      LMAX=2*Truncation+M0-1

      I=0
      DO  LI=LMIN,LMAX,LD
         ! L for toroidal (w) field:
         IF( Symmetry.EQ.2 ) THEN
             LTI=LI+1
         ELSEIF( Symmetry.EQ.1 ) THEN
             LTI=LI-1
         ELSEIF( Symmetry.EQ.0 ) THEN
             LTI=LI
         ENDIF
         LPI = LI
         NIMAX=INT( DBLE(2*Truncation+1-LI+M0)/2 )
         DO  NI=1,NIMAX
            WRITE(unitFlowPolOut, '(3I5,2D18.8)') LPI,M0,NI, DBLE(ZEVEC(I+1)),DIMAG(ZEVEC(I+1))
            WRITE(unitFlowTorOut, '(3I5,2D18.8)') LTI,M0,NI, DBLE(ZEVEC(I+3)),DIMAG(ZEVEC(I+3))
            WRITE(unitTempOut,    '(3I5,2D18.8)') LI, M0,NI, DBLE(ZEVEC(I+2)),DIMAG(ZEVEC(I+2))
            WRITE(unitCompOut,    '(3I5,2D18.8)') LI, M0,NI, DBLE(ZEVEC(I+4)),DIMAG(ZEVEC(I+4))
            I=I+4
         enddo
      enddo
   end subroutine

   !****************************************************************
   !     subroutine write_parp(...)
   !
   !>     writes the parameter file 'file'.par
   !****************************************************************
   subroutine write_parp(unit_out)
      implicit none
      integer, intent(in):: unit_out
      double precision:: TaOut
      TaOut=Ta*Ta
      write(unit_out,'(A)') '| magic | lformat | drs_calc_type | tempBC | compBC | flowBC | magBC |'
      write(unit_out,'(7I8)')    0,     0,         9,           0,       0,       0,       0

      write(unit_out,'(A)') '|  eta  | Therm. Prandtl | Compo. Prandtl | Taylor | Therm. Rayleigh | Compo. Rayleigh |  Magn. Prandtl  |'
      write(unit_out,'(7D14.5)') eta,     Pt,                 Pc,         TaOut,       Ra_t,             Ra_c,           1
      write(unit_out,'(A)') '| Nr | Nt | Np | Nr_s | Nt_s |Np_s | lsymm | m0 |'
      write(unit_out,'(8I8)')  Nr,   Nt, Np,  Nr_s,  NIMAX, Np_s, Symmetry,   m0
      write(unit_out,'(A)') '|delta_t|nsteps|transient|sampling_rate|step|time|drift|'
      write(unit_out,'(D11.3,I7,I7,I6,I8,D14.5,D11.3)') &
      h, stepmax, transient, sampling_rate, steps, time, drift
      write(unit_out,*) '''',comment,''''

   end subroutine write_parp
end program
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
