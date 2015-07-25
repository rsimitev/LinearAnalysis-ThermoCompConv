module GrowthRateRaEffMod
#include "errorcodes.h"
   implicit none
   private
   double precision, parameter:: DPI=3.141592653589793D0
   !> Internally used parameters
   double precision:: Le_i, Pt_i, tau_i, Ra_i, alpha_i
   double precision:: ri, ro, eta_i
   character(len=3):: variable='Rt'
   integer:: Symmetry_i, mm_i, NEigenmodes, lmin, ld, Truncation_i
   public:: MaxGrowthRate, MaxGrowthRateCmplx, computeGrowthRateModes
   public:: setEigenProblemSize, getEigenProblemSize, GrowthRateInit, setVariableParam
   public:: testMAT, setLminAndLD, GrowthRateUpdatePar
   public:: setParameterValue,saveParameterValue
contains

   !***********************************************************************
   !> Initializes the module with fixed parameters.
   subroutine GrowthRateInit(Ra, alpha, Pt, Le, tau, eta, m, Symmetry, truncation)
      implicit none
      double precision, intent(in)::Ra, alpha, Pt, Le, tau, eta
      integer, intent(in):: m, Symmetry, truncation
      Ra_i  = Ra
      alpha_i  = alpha
      Pt_i  = Pt
      Le_i  = Le
      tau_i = tau
      mm_i  = m
      eta_i = eta
      Symmetry_i = Symmetry
      truncation_i = truncation
      ri = eta_i/(1.0d0-eta_i)
      ro = 1.0D0 + ri
      call setLminAndLD(Symmetry, m, LMIN, LD)
      call setEigenProblemSize(LMIN,LD,truncation,M)
   end subroutine

   !***********************************************************************
   !> Updates the specified parameters.
   !> @par Ra is the effective Rayleigh number
   !> @par alpha is the Rayleigh angle
   subroutine GrowthRateUpdatePar(Ra, alpha, Pt, Le, tau, eta, m, Symmetry, truncation)
      implicit none
      double precision, optional, intent(in)::Pt, Le, tau, eta, Ra, alpha
      integer,optional, intent(in):: m, Symmetry, truncation
      if(present(Pt))  Pt_i  = Pt
      if(present(Le))  Le_i  = Le
      if(present(tau)) tau_i = tau
      if(present(Ra))  Ra_i  = Ra
      if(present(alpha))  alpha_i  = alpha
      if(present(eta)) then
         eta_i = eta
         ri = eta_i/(1.0d0-eta_i)
         ro = 1.0D0 + ri
      endif
      if(present(m))          mm_i         = m
      if(present(Symmetry))   Symmetry_i   = Symmetry
      if(present(truncation)) truncation_i = truncation
      if(present(m) .or. present(Symmetry) .or.  present(truncation)) then
         call setLminAndLD(Symmetry_i, mm_i, LMIN, LD)
         call setEigenProblemSize(LMIN,LD,truncation_i,mm_i)
      endif
   end subroutine

   !***********************************************************************
   subroutine writeAllParameters()
      implicit none
      Write(*,*) '---BEGIN'
      Write(*,*) 'WAP1: ', Ra_i, alpha_i, Pt_i
      Write(*,*) 'WAP2: ', Le_i, tau_i, eta_i
      Write(*,*) 'WAP3: ', ri, ro
      Write(*,*) 'WAP4: ', mm_i, Symmetry_i, truncation_i
      Write(*,*) '---END'
   end subroutine

   !***********************************************************************
   subroutine setLminAndLD(Sym, mm, LMIN, LD)
      implicit none
      integer, intent(in):: Sym, mm
      integer, intent(out):: LMIN, LD
      IF(Sym.EQ.0) THEN
         ! - UNDEFINED SYMMETRIE:
         LMIN = mm
         LD   = 1
      ELSEIF(Sym.EQ.1) THEN
         ! - EQUATORIAL ANTISYMMETRIE (L+M ODD):
         LMIN = mm+1
         LD   = 2
      ELSEIF(Sym.EQ.2) THEN
         ! - EQUATORIAL SYMMETRIE (L+M EVEN);
         LMIN = mm
         LD   = 2
      ENDIF
   end subroutine

   !***************************************************************************
   ! - DETERMINATION OF DIMENSION of the problem:
   ! - for each value of L the number of possible N-values is added
   SUBROUTINE setEigenProblemSize(LMIN,LD,truncation,M)
      implicit none
      integer, intent(in):: LMIN,LD,truncation,M
      integer:: L
      ! - DETERMINATION OF DIMENSION:
      ! - for each value of L the number of possible N-values is added
      !         print*, "Triangular truncation (2.12)"
      !         print*, LMIN, "...", 2*truncation_i+M-1,LD
      NEigenmodes=0
      DO L = LMIN, 2*truncation+M-1, LD
         NEigenmodes = NEigenmodes + 4*INT( DBLE(2*truncation+1-L+M)/2 )
      endDO
   end subroutine

   !***************************************************************************
   integer function getEigenProblemSize()
      implicit none
      getEigenProblemSize = NEigenModes
   end function


   !***********************************************************************
   !> Sets which parameter is going to be varied.
   !! Possible values for par are Ra, aa, Rt, Rc, Pt, Le, tau, eta and m.\n
   !! The default value is 'Rt'.
   subroutine setVariableParam(par)
      implicit none
      character(len=3), intent(in):: par
      select case(trim(par))
         case ('Ra','aa','Pt','Le','tau','eta','m')
            variable = par
         case default
            variable = 'Ra'
      end select
   end subroutine

   !***********************************************************************
   !> Sets the internal parameter selected by setVariableParam(par) to the value
   !! \a val
   subroutine setParameterValue(val)
      implicit none
      double precision, intent(in):: val
      select case(trim(variable))
         case ('Ra')
            Ra_i = val
         case ('aa')
            alpha_i = val
         case ('Pt')
            Pt_i = val
         case ('Le')
            Le_i = val
         case ('tau')
            tau_i = val
         case ('eta')
            eta_i = val
         case ('m')
            mm_i = nint(val)
         case default
            Ra_i = val
      end select
   end subroutine

   !***********************************************************************
   !> Saves the value of the internal parameter selected by setVariableParam(par) in
   !! \a val
   subroutine saveParameterValue(val)
      implicit none
      double precision, intent(out):: val
      select case(trim(variable))
         case ('Ra')
            val = Ra_i
         case ('aa')
            val = alpha_i
         case ('Pt')
            val = Pt_i
         case ('Le')
            val = Le_i
         case ('tau')
            val = tau_i
         case ('eta')
            val = eta_i
         case ('m')
            val = dble(mm_i)
         case default
            val = Ra_i
      end select
   end subroutine

   !***********************************************************************
   !> Computes the maximum imaginary part of the frequency,
   !! that is, the maximum growth rate for all eigen modes.
   double precision FUNCTION MaxGrowthRate(val)
      implicit none
      double precision, intent(in):: val
      double precision:: valOld
      double complex:: ZEW(NEigenmodes)
      integer:: imin

      call saveParameterValue(valOld)
      call setParameterValue(val)
      call computeGrowthRateModes(.false., zew)
      ! search for lowest imaginary part:
      IMIN = minloc(DIMAG(ZEW),1)
      MaxGrowthRate = DIMAG(ZEW(imin))
      call setParameterValue(valOld)
   end function

   !***********************************************************************
   !> Computes the complex frequency for which the growth rate is maximum.
   double complex FUNCTION MaxGrowthRateCmplx(val)
      implicit none
      double precision, intent(in):: val
      double precision:: valOld
      double complex:: ZEW(NEigenmodes)
      integer:: imin

      call saveParameterValue(valOld)
      call setParameterValue(val)
      Zew = dcmplx(0.0d0,0.0d0)
      call computeGrowthRateModes(.false., zew)
      ! search for lowest imaginary part:
      IMIN = minloc(DIMAG(ZEW),1)
      MaxGrowthRateCmplx = ZEW(IMIN)
      call setParameterValue(valOld)
   end function

   !***********************************************************************
   !> Computes the eigenvalues and optionally the eigen modes of the
   !! algebraic problem. The real part of the eigen value is the oscillation
   !! frequency and the imaginary part is the symmetric of the growth rate.
   subroutine computeGrowthRateModes(sort, zew, zeval)
      implicit none
      !> .True. Sort the eigenvalues and eigenvectors if computed.
      logical, intent(in):: sort
      !> Stores the eigenvalues.
      double complex, intent(out):: ZEW(NEigenmodes)
      !> If present, stores the eigenvectors ordered as the eigenvalues.
      double complex, intent(out), optional::ZEVAL(NEigenmodes,NEigenmodes)
      double complex:: ZA(NEigenmodes,NEigenmodes),ZB(NEigenmodes,NEigenmodes)
      double complex:: ZEWA(NEigenmodes),ZEWB(NEigenmodes),ZEVALL(NEigenmodes,NEigenmodes)
      double complex:: ZEVEC(NEigenmodes), ZSAVE, ZWORK(3*NEigenmodes)
      double precision:: RWORK(8*NEigenmodes)
      integer:: i, j, k, info
      
      ! - MAT SETS THE complex(8) MATRICES ZA AND ZB SETTING OF MATRIX:
      CALL MAT(tau_i, Ra_i, alpha_i, Pt_i, Le_i, mm_i, Symmetry_i, ZA,ZB, NEigenmodes)

!       SUBROUTINE zGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,
!     $                  VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
      IF(.not.present(zeval)) THEN ! Compute eigen values and vectors.
        call zggev('N', 'N',NEigenmodes,ZA,NEigenmodes,ZB,NEigenmodes, ZEWA, ZEWB, ZEVALL, NEigenmodes, ZEVALL, NEigenmodes, ZWORK, 3*NEigenmodes, rwork, info)
      ELSE ! Only compute eigenvalues
        call zggev('N', 'V',NEigenmodes,ZA,NEigenmodes,ZB,NEigenmodes, ZEWA, ZEWB, ZEVALL, NEigenmodes, ZEVALL, NEigenmodes, ZWORK, 3*NEigenmodes, rwork, info)
        zeval = zevall
      endIF

      ZEW(:)=ZEWA(:)/ZEWB(:)

      ! sort eigenvalues:
      if(sort) then
         DO I=1,NEigenmodes
            DO J=I,NEigenmodes
               IF( DIMAG(ZEW(J)).LT.DIMAG(ZEW(I)) ) THEN
                  ZSAVE = ZEW(J)
                  ZEW(J) = ZEW(I)
                  ZEW(I) = ZSAVE
                  IF(present(zeval)) THEN
                     DO K=1,NEigenmodes
                        ZEVEC(K)   = ZEVAL(K,J)
                        ZEVAL(K,J) = ZEVAL(K,I)
                        ZEVAL(K,I) = ZEVEC(K)
                     enddo
                  endIF
               endIF
            enddo
         enddo
      endif
   end subroutine computeGrowthRateModes

   !************************************************************************
   ! - SETS THE complex(8) MATRICES ZA AND ZB.
   SUBROUTINE MAT(tau_i, Ra_i, alpha_i, Pt_i, Le_i, mm_i, Symmetry_i, ZA,ZB,NDIM)
      implicit none
      integer, intent(in):: ndim, mm_i, Symmetry_i
      double complex, intent(out):: ZA(ndim,ndim),ZB(ndim,ndim)
      double precision, intent(in):: tau_i, Ra_i, alpha_i, Pt_i, Le_i
      double precision:: Rt_i, Rc_i
      integer:: lmax
      integer:: ni, i, nimax, li, lpi, lti
      integer:: nj, j, njmax, lj, lpj, ltj

      Rt_i = Ra_i*cos(alpha_i)
      Rc_i = Ra_i/Le_i*sin(alpha_i)

      ZA(:,:)=DCMPLX(0D0,0D0)
      ZB(:,:)=DCMPLX(0D0,0D0)

      I=0
      LMAX=2*truncation_i+mm_i-1
      !Write(*,*) 'MAT():', mm_i, LMIN, LMAX, LD
      !Write(*,*) 'MAT():', Symmetry_i, truncation_i, NDIM
      DO LI=LMIN,LMAX,LD
         LPI=LI
         ! Determine L for toroidal (w) field:
         IF( Symmetry_i.EQ.2 ) THEN
            LTI=LI+1
         ELSEIF( Symmetry_i.EQ.1 ) THEN
            LTI=LI-1
         ELSEIF( Symmetry_i.EQ.0 ) THEN
            LTI=LI
         endIF
         NIMAX=INT( DBLE(2*truncation_i+1-LI+mm_i)/2 )

         DO NI=1,NIMAX
            J=0
            DO LJ=LMIN,LMAX,LD
               LPJ=LJ
               IF( Symmetry_i.EQ.2 ) THEN
                  LTJ=LJ+1
               ELSEIF( Symmetry_i.EQ.1 ) THEN
                  LTJ=LJ-1
               ELSEIF( Symmetry_i.EQ.0 ) THEN
                  LTJ=LJ
               endIF
               NJMAX=INT( DBLE(2*truncation_i+1-LJ+mm_i)/2 )

               !  ******************** I: Equation (Line) ******************
               !  ******************** J: Variable (Column) ****************
               !  ******************** I+1: v (poloidal)  ******************
               !  ******************** I+2: theta         ******************
               !  ******************** I+3: w (toroidal)  ******************
               ! new****************** I+4: gamma (concentration) **********
               DO NJ=1,NJMAX
                  IF(J+3.GT.NDIM .OR. I+3.GT.NDIM) THEN
                     write(*,*) 'MAT(): NDIM too small.'
                     Write(*,*) 'i =',i,'j =', j, 'NDIM =', NDIM
                     stop
                  endIF
                  IF( LI.EQ.LJ ) THEN
                     ! these are the terms on the LHS
                     ZB(I+1,J+1) = DCMPLX(0.D0, -pol_LHS(          NI, NJ, LPI, 1))
                     ZB(I+2,J+2) = DCMPLX(0.D0, -scalar_LHS( Pt_i, NI, NJ, 1)) ! temperature
                     ZB(I+3,J+3) = DCMPLX(0.D0, -tor_LHS   (       NI, NJ, LTI, 1))
                     ZB(I+4,J+4) = DCMPLX(0.D0, -scalar_LHS( Pt_i, NI, NJ, 1)) ! concentration

                     ! NS
                     ZA(I+1,J+1) = DCMPLX(Diffusion_pol(NI,NJ,LPI),         Coriolis_pol2pol(tau_i, mm_i, NI, NJ, LPI, 1))
                     ZA(I+3,J+3) = DCMPLX(Diffusion_tor(NI,NJ,LTI),         Coriolis_tor2tor(tau_i, mm_i, NI, NJ, 1))
                     ZA(I+1,J+2) = DCMPLX(Buoyancy_temp(Rt_i, NI, NJ, LPI), 0.D0)
                     ZA(I+1,J+4) = DCMPLX(Buoyancy_conc(Rc_i, NI, NJ, LPI), 0.D0)
                     ! HE
                     ZA(I+2,J+1) = DCMPLX(Advection_scalar(NI, NJ, LPI),    0.D0)
                     ZA(I+2,J+2) = DCMPLX(Diffusion_scalar(NI, NJ, LPI),    0.D0)
                     ! CE
                     ZA(I+4,J+1) = DCMPLX(Advection_scalar(NI,NJ,LPI),      0.D0)
                     ZA(I+4,J+4) = DCMPLX(Diffusion_scalar(NI,NJ,LPI)/Le_i, 0.D0)
                  endIF
                  ! The remaining cross-pol/tor terms of the Coriolis force.
                  IF( LPI.EQ.LTJ+1 ) THEN
                     ZA(I+1,J+3) = DCMPLX(DIII4A(tau_i, mm_i, NI, NJ, LPI, 1),  0.D0)
                  ELSEIF( LPI.EQ.LTJ-1 ) THEN
                     ZA(I+1,J+3) = DCMPLX(DIII4B(tau_i, mm_i, NI, NJ, LPI, 1),  0.D0)
                  endIF
                  IF( LTI.EQ.LPJ+1 ) THEN
                     ZA(I+3,J+1) = DCMPLX(DII4A(tau_i, mm_i, NI, NJ, LTI, 1),   0.D0)
                  ELSEIF( LTI.EQ.LPJ-1 ) THEN
                     ZA(I+3,J+1) = DCMPLX(DII4B(tau_i, mm_i, NI, NJ, LTI, 1),   0.D0)
                  endIF
                  J=J+4
               enddo
            enddo
             I=I+4
         enddo
      enddo
   end subroutine

   subroutine testMAT(ZA_refFile, ZB_refFile, error)
      implicit none
      integer, parameter:: Nmodes=220
      character(len=*), intent(in):: ZA_refFile, ZB_refFile
      integer, intent(out):: error
      double complex:: ZA(Nmodes,Nmodes),ZB(Nmodes,Nmodes)
      double complex:: ZA_ref,ZB_ref
      double precision:: Rt_i, Rc_i
      integer:: i,j,k1, k2
      integer:: lmax
      integer:: ni, nimax, li, lpi, lti
      integer:: nj, njmax, lj, lpj, ltj
      tau_i=1.0d5
      Rt_i=3.0d6
      Rc_i=1.0d2
      Pt_i=1.0d0
      Le_i=0.30d0
      eta_i=0.35d0
      mm_i=18
      Symmetry_i=2
      truncation_i=10
      Write(*,*) ZA_refFile,' ', ZB_refFile
      Write(*,*) tau_i, Rt_i, Rc_i, Pt_i, Le_i, mm_i, Symmetry_i, Nmodes
      RI = ETA_i/(1.0d0-ETA_i)
      RO = 1.0D0 + RI
      call setLminAndLD(Symmetry_i, mm_i, LMIN, LD)
      call setEigenProblemSize(LMIN,LD,truncation_i,Mm_i)
      Write(*,*) 'ri =', ri,'ro =', ro, Nmodes, NEigenmodes
      CALL MAT(tau_i, Rt_i, Rc_i, Pt_i, Le_i, mm_i, Symmetry_i, ZA,ZB, NEigenmodes)
      Write(*,*) 'Done with MAT.'

      I=0
      error=0
      LMAX = 2*truncation_i+mm_i-1
      open(unit=444, file=trim(ZA_refFile), status='OLD')
      open(unit=445, file=trim(ZB_refFile), status='OLD')
      DO LI=LMIN,LMAX,LD
         LPI=LI
         ! Determine L for toroidal (w) field:
         IF( Symmetry_i.EQ.2 ) THEN
            LTI=LI+1
         ELSEIF( Symmetry_i.EQ.1 ) THEN
            LTI=LI-1
         ELSEIF( Symmetry_i.EQ.0 ) THEN
            LTI=LI
         endIF
         NIMAX=INT( DBLE(2*truncation_i+1-LI+mm_i)/2 )

         DO NI=1,NIMAX
            J=0
            DO LJ=LMIN,LMAX,LD
               LPJ=LJ
               IF( Symmetry_i.EQ.2 ) THEN
                  LTJ=LJ+1
               ELSEIF( Symmetry_i.EQ.1 ) THEN
                  LTJ=LJ-1
               ELSEIF( Symmetry_i.EQ.0 ) THEN
                  LTJ=LJ
               endIF
               NJMAX=INT( DBLE(2*truncation_i+1-LJ+mm_i)/2 )

               !  ******************** I: Equation (Line) ******************
               !  ******************** J: Variable (Column) ****************
               !  ******************** I+1: v (poloidal)  ******************
               !  ******************** I+2: theta         ******************
               !  ******************** I+3: w (toroidal)  ******************
               ! new****************** I+4: gamma (concentration) **********
               DO NJ=1,NJMAX
                  do k1=1, 4
                     do k2=1, 4
                        Read(444,*) ZA_ref
                        Write(*,*) 'ZA(',i+k1,',',j+k2,')=',ZA(i+k1,j+k2),'.vs.',ZA_ref
                        if(abs(ZA(i+k1,j+k2)-ZA_ref)/abs(ZA_ref).gt.1.0e-10) then
                           error=1
                           return
                        endif
                        Read(445,*) ZB_ref
                        Write(*,*) 'ZB(',i+k1,',',j+k2,')=',ZB(i+k1,j+k2),'.vs.',ZB_ref
                        if(abs(ZB(i+k1,j+k2)-ZB_ref)/abs(ZB_ref).gt.1.0e-10) then
                           error=2
                           return
                        endif
                     enddo
                  enddo
                  J=J+4
               enddo
            enddo
             I=I+4
         enddo
      enddo
      close(444)
      close(445)
   end subroutine

   !************************************************************************
   ! - GALERKIN TERMS:
   !************************************************************************
   double precision function scalar_LHS(Pt, N1,N2,NU1)
   ! ---- HEAT EQUATION, TIME DERIVATIVE
      implicit none
      double precision, intent(in):: Pt
      integer, intent(in):: N1, N2, NU1
      scalar_LHS=Pt*NU1*R('SS ',2,N1,N2,0)
   end function

   double precision function Diffusion_scalar(N1,N2,L1)
   ! ---- HEAT EQUATION , DISSIPATION
      implicit none
      integer:: N1, N2, l1
      Diffusion_scalar=N2**2*DPI**2*R('SS ',2,N1,N2,0) - 2*N2*DPI*R('SC ',1,N1,N2,0) + llp1(L1)*R('SS ',0,N1,N2,0)
   end function

   double precision function Advection_scalar(N1,N2,L1)
   ! ---- HEAT EQUATION , SOURCE
      implicit none
      integer:: N1, N2, l1
      Advection_scalar =-llp1(L1)*R('SS ',2,N1,N2,0)
   end function

   double precision function Diffusion_tor(N1,N2,L1)
   ! ---- TOROIDAL EQUATION , DISSIPATION
      implicit none
      integer:: N1, N2, l1
      Diffusion_tor=llp1(L1)*( (N2-1)**2*DPI**2*R('CC ',4,N1-1,N2-1,0)+4*(N2-1)*DPI*R('CS ',3,N1-1,N2-1,0)+(llp1(L1)-2)*R('CC ',2,N1-1,N2-1,0) )
   end function

   double precision function tor_LHS(N1,N2,L1,NU1)
   ! ---- TOROIDAL EQUATION , TIME DERIVATIVE
      implicit none
      integer:: N1, N2, l1, NU1
      tor_LHS=NU1*llp1(L1)*R('CC ',4,N1-1,N2-1,0)
   end function

   double precision function Coriolis_tor2tor(tau, mm, N1,N2,NU1)
   ! ---- TOROIDAL EQUATION , CORRIOLIS
      implicit none
      double precision, intent(in):: tau
      integer, intent(in):: N1, N2, NU1, mm
      Coriolis_tor2tor=-tau*NU1*mm*R('CC ',4,N1-1,N2-1,0)
   end function

   double precision function DII4A(tau,mm, N1,N2,L1,NU1)
   ! ---- TOROIADL EQUATION , Q-TERM 1 (L1=L3+1)
      implicit none
      double precision, intent(in):: tau
      integer, intent(in):: N1, N2, l1, NU1, mm
      DII4A= tau*DSQRT( DBLE(L1-NU1*mm)*(L1+NU1*mm)/(2*L1-1)/(2*L1+1) )*( (L1**2-1)*(L1-1)*R('CS ',2,N1-1,N2,0)  - (L1+1)*(L1-1)*N2*DPI*R('CC ',3,N1-1,N2,0)    )
   end function

   double precision function DII4B(tau,mm,  N1,N2,L1,NU1)
   ! ---- TOROIADL EQUATION , Q-TERM 1 (L1=L3-1)
      implicit none
      double precision, intent(in):: tau
      integer, intent(in):: N1, N2, l1, NU1, mm
      DII4B= tau*DSQRT( DBLE(L1-NU1*mm+1)*(L1+NU1*mm+1)/(2*L1+1)/(2*L1+3) )*( (1-(L1+1)**2)*(L1+2)*R('CS ',2,N1-1,N2,0)  - L1*(L1+2)*N2*DPI*R('CC ',3,N1-1,N2,0)  )
   end function

   double precision function Diffusion_pol(N1,N2,L1)
   ! ---- POLOIDAL EQUOATION , DISSIPATION
      implicit none
      integer:: N1, N2, l1
      Diffusion_pol=llp1(L1)* ( N2**4*DPI**4*R('SS ',2,N1,N2,0) - 4*N2**3*DPI**3*R('SC ',1,N1,N2,0)+2*llp1(L1)*N2**2*DPI**2*R('SS ',0,N1,N2,0)+(llp1(L1)**2-2*llp1(L1))*R('SS ',-2,N1,N2,0) )
   end function

   ! ---- POLOIDAL EQUATION , TIME DERIVATIVE
   double precision function pol_LHS(N1,N2,L1,NU1)
      implicit none
      integer:: N1, N2, l1, NU1
      pol_LHS= -NU1*llp1(L1)*( -N2**2*DPI**2*R('SS ',2,N1,N2,0) + 2*N2*DPI*R('SC ',1,N1,N2,0) - llp1(L1)*R('SS ',0,N1,N2,0) )
   end function

   double precision function Coriolis_pol2pol(tau, mm, N1,N2,L1,NU1)
   ! ---- POLOIDAL EQUATION , CORRIOLIS
      implicit none
      double precision, intent(in):: tau
      integer, intent(in):: N1, N2, l1, NU1, mm
      Coriolis_pol2pol= tau*NU1*mm*( -N2**2*DPI**2*R('SS ',2,N1,N2,0)+2*N2*DPI*R('SC ',1,N1,N2,0) - llp1(L1)*R('SS ',0,N1,N2,0) )
   end function

   double precision function DIII4A(tau,mm,  N1, N2, L1, NU1)
   ! ---- POLOIDAL EQUATION , Q-TERM 1 (L1=L3+1)
      implicit none
      double precision, intent(in):: tau
      integer, intent(in):: N1, N2, l1, NU1, mm
      DIII4A= tau*DSQRT( DBLE(L1-mm*NU1)*(L1+mm*NU1)/(2*L1-1)/(2*L1+1) ) * ( (L1*(L1-1)-2)*(L1-1)*R('SC ',2,N1,N2-1,0) + (L1+1)*(L1-1)*(N2-1)*DPI*R('SS ',3,N1,N2-1,0) )
   end function

   double precision function DIII4B(tau,mm, N1, N2, L1, NU1)
   ! ---- POLOIDAL EQUATION , Q-TERM 2 (L1=L3-1)
      implicit none
      double precision, intent(in):: tau
      integer, intent(in):: N1, N2, l1, NU1, mm
      DIII4B= tau*DSQRT( DBLE(L1-mm*NU1+1)*(L1+mm*NU1+1)/(2*L1+1)/(2*L1+3) ) * ( (L1+2)*(2-(L1+1)*(L1+2))*R('SC ',2,N1,N2-1,0) + L1*(L1+2)*(N2-1)*DPI*R('SS ',3,N1,N2-1,0) )
   end function

   double precision function Buoyancy_temp(Ra,N1,N2,L1)
   ! ---- POLOIDAL EQUATION ,
      implicit none
      integer:: N1, N2, l1
      double precision:: Ra
      Buoyancy_temp=-Ra*llp1(L1)*R('SS ',2,N1,N2,0)
   end function

   double precision function Buoyancy_conc(Ra,N1,N2,L1)
   ! ---- POLOIDAL EQUATION ,
      implicit none
      integer:: N1, N2, l1
      double precision:: Ra
      Buoyancy_conc=-Ra*llp1(L1)*R('SS ',2,N1,N2,0)
   end function

   !**************************************************************************
   !-- SUBROUTINES:
   !**************************************************************************
   pure double precision function llp1(L)
      implicit none
      integer, intent(in):: l
      llp1 = DBLE(L*(L+1))
   end function

!-----------------------------------------------------------------------
   pure double precision function R (TRII, Symmetry_i, II, JI, KI)
!-----------------------------------------------------------------------
!     BERECHNET DIE RADIALINTEGRALE FUER DAS DYNAMOPROBLEM
!-----------------------------------------------------------------------
      implicit none
      CHARACTER(3), intent(in):: TRII
      integer, intent(in):: Symmetry_i, ii, ji, ki

      CHARACTER(3):: TRI
      integer:: i,j,k
      double precision:: Rint
!
      TRI = TRII
      I = II
      J = JI
      K = KI
      CALL tausch (TRI, I, J, K)
      ! TODO: Convert to select case.
      IF (Symmetry_i.EQ.0) THEN
         IF (TRI.EQ.'SS ') THEN
            IF (I*J.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0*(delta(I - J) - delta(I+J) )
            ENDIF
         ELSEIF (TRI.EQ.'SC ') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0*(two_over_nPi_odd (I - J)+two_over_nPi_odd (I+J) )
            ENDIF
         ELSEIF (TRI.EQ.'CC ') THEN
            RINT = 1D0 / 2D0*(delta (I+J)+delta (I - J) )
         ELSEIF (TRI.EQ.'SSS') THEN
            IF (I*J*K.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0*( two_over_nPi_odd(I-J+K) + two_over_nPi_odd(I+J-K) - two_over_nPi_odd(I-J-K) - two_over_nPi_odd(I+J+K) )
            ENDIF
         ELSEIF (TRI.EQ.'SSC') THEN
            IF (I*J.NE.0) THEN
               RINT = 1D0 / 4D0*(delta (I - J - K) - delta (I+J+K) +delta (I - J+K) - delta (I+J - K) )
            ELSE
               RINT = 0D0
            ENDIF
         ELSEIF (TRI.EQ.'SCC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0*(two_over_nPi_odd (I+J - K)+two_over_nPi_odd (I+J+K) +two_over_nPi_odd (I - J+K)+two_over_nPi_odd (I - J - K) )
            ENDIF
         ELSEIF (TRI.EQ.'CCC') THEN
            RINT = 1D0 / 4D0*( delta (I+J+K) + delta(I+J-K) + delta(I-J+K) + delta (I-J-K) )
         ENDIF
      ELSEIF (Symmetry_i.EQ.1) THEN
         IF (TRI.EQ.'SS ') THEN
            IF (I*J.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0*(C1 (I - J) - C1 (I+J) )
            ENDIF
         ELSEIF (TRI.EQ.'SC ') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0*(S1 (I - J)+S1 (I+J) )
            ENDIF
         ELSEIF (TRI.EQ.'CC ') THEN
            RINT = 1D0 / 2D0*(C1 (I+J)+C1 (I - J) )
         ELSEIF (TRI.EQ.'SSS') THEN
            IF (I*J*K.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0*(S1 (I - J+K)+S1 (I+J - K) - S1 (I - J - K) - S1 (I+J+K) )
            ENDIF
         ELSEIF (TRI.EQ.'SSC') THEN
            IF (I*J.NE.0) THEN
               RINT = 1D0 / 4D0*(C1 (I - J - K) - C1 (I+J+K) +C1 (I - J+K) - C1 (I+J - K) )
            ELSE
               RINT = 0D0
            ENDIF
         ELSEIF (TRI.EQ.'SCC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0*(S1 (I+J - K)+S1 (I+J+K) +S1 (I - J+K)+S1 (I - J - K) )
            ENDIF
         ELSEIF (TRI.EQ.'CCC') THEN
            RINT = 1D0 / 4D0*(C1 (I+J+K)+C1 (I+J - K)+C1 (I - J+K)+C1 (I - J - K) )
         ENDIF
      ELSEIF (Symmetry_i.EQ.2) THEN
         IF (TRI.EQ.'SS ') THEN
            IF (I*J.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0*(C2 (I - J) - C2 (I+J) )
            ENDIF
         ELSEIF (TRI.EQ.'SC ') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0*(S2 (I - J)+S2 (I+J) )
            ENDIF
         ELSEIF (TRI.EQ.'CC ') THEN
            RINT = 1D0 / 2D0*(C2 (I+J)+C2 (I - J) )
         ELSEIF (TRI.EQ.'SSS') THEN
            IF (I*J*K.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0*(S2 (I - J+K)+S2 (I+J - K) - S2 (I - J - K) - S2 (I+J+K) )
            ENDIF
         ELSEIF (TRI.EQ.'SSC') THEN
            IF (I*J.NE.0) THEN
               RINT = 1D0 / 4D0*(C2 (I - J - K) - C2 (I+J+K) +C2 (I - J+K) - C2 (I+J - K) )
            ELSE
               RINT = 0D0
            ENDIF
         ELSEIF (TRI.EQ.'SCC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0*(S2 (I+J - K)+S2 (I+J+K) +S2 (I - J+K)+S2 (I - J - K) )
            ENDIF
         ELSEIF (TRI.EQ.'CCC') THEN
            RINT = 1D0 / 4D0*(C2 (I+J+K)+C2 (I+J - K)+C2 (I - J+K)+C2 (I - J - K) )
         ENDIF
      ELSEIF (Symmetry_i.EQ.3) THEN
         IF (TRI.EQ.'SS') THEN
            IF (I*J.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0*(C3 (I - J) - C3 (I+J) )
            ENDIF
         ELSEIF (TRI.EQ.'SC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0*(S3 (I - J)+S3 (I+J) )
            ENDIF
         ELSEIF (TRI.EQ.'CC') THEN
            RINT = 1D0 / 2D0*(C3 (I+J)+C3 (I - J) )
         ELSEIF (TRI.EQ.'SSS') THEN
            IF (I*J*K.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0*(S3 (I - J+K)+S3 (I+J - K) - S3 (I - J - K) - S3 (I+J+K) )
            ENDIF
         ELSEIF (TRI.EQ.'SSC') THEN
            IF (I*J.NE.0) THEN
               RINT = 1D0 / 4D0*(C3 (I - J - K) - C3 (I+J+K) +C3 (I - J+K) - C3 (I+J - K) )
            ELSE
               RINT = 0D0
            ENDIF
         ELSEIF (TRI.EQ.'SCC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0*(S3 (I+J - K)+S3 (I+J+K) +S3 (I - J+K)+S3 (I - J - K) )
            ENDIF
         ELSEIF (TRI.EQ.'CCC') THEN
            RINT = 1D0 / 4D0*(C3 (I+J+K)+C3 (I+J - K)+C3 (I - J+K)+C3 (I - J - K) )
         ENDIF
      ELSEIF (Symmetry_i.EQ.4) THEN
         IF (TRI.EQ.'SS') THEN
            IF (I*J.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0*(C4 (I - J) - C4 (I+J) )
            ENDIF
         ELSEIF (TRI.EQ.'SC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0*(S4 (I - J)+S4 (I+J) )
            ENDIF
         ELSEIF (TRI.EQ.'CC') THEN
            RINT = 1D0 / 2*(C4 (I+J)+C4 (I - J) )
         ELSEIF (TRI.EQ.'SSS') THEN
            IF (I*J*K.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0*(S4 (I - J+K)+S4 (I+J - K)      &
               - S4 (I - J - K) - S4 (I+J+K) )
            ENDIF
         ELSEIF (TRI.EQ.'SSC') THEN
            IF (I*J.NE.0) THEN
               RINT = 1D0 / 4D0*(C4 (I - J - K) - C4 (I+J+K)      &
              +C4 (I - J+K) - C4 (I+J - K) )
            ELSE
               RINT = 0D0
            ENDIF
         ELSEIF (TRI.EQ.'SCC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0*(S4 (I+J - K)+S4 (I+J+K)      &
              +S4 (I - J+K)+S4 (I - J - K) )
            ENDIF
         ELSEIF (TRI.EQ.'CCC') THEN
            RINT = 1D0 / 4D0*(C4 (I+J+K)+C4 (I+J - K)+C4 (I &
            - J+K)+C4 (I - J - K) )
         ENDIF
      ELSEIF (Symmetry_i.EQ. - 1) THEN
         IF (TRI.EQ.'SS ') THEN
            IF (I*J.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0*(CM1 (I - J) - CM1 (I+J) )
            ENDIF
         ELSEIF (TRI.EQ.'SC ') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0*(SM1 (I - J)+SM1 (I+J) )
            ENDIF
         ELSEIF (TRI.EQ.'CC ') THEN
            RINT = 1D0 / 2D0*(CM1 (I+J)+CM1 (I - J) )
         ELSEIF (TRI.EQ.'SSS') THEN
            IF (I*J*K.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0*(SM1 (I - J+K)+SM1 (I+J - K)    &
               - SM1 (I - J - K) - SM1 (I+J+K) )
            ENDIF
         ELSEIF (TRI.EQ.'SSC') THEN
            IF (I*J.NE.0) THEN
               RINT = 1D0 / 4D0*(CM1 (I - J - K) - CM1 (I+J+K)    &
              +CM1 (I - J+K) - CM1 (I+J - K) )
            ELSE
               RINT = 0D0
            ENDIF
         ELSEIF (TRI.EQ.'SCC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0*(SM1 (I+J - K)+SM1 (I+J+K)    &
              +SM1 (I - J+K)+SM1 (I - J - K) )
            ENDIF
         ELSEIF (TRI.EQ.'CCC') THEN
            RINT = 1D0 / 4D0*(CM1 (I+J+K)+CM1 (I+J - K)       &
           +CM1 (I - J+K)+CM1 (I - J - K) )
         ENDIF
      ELSEIF (Symmetry_i.EQ. - 2) THEN
         IF (TRI.EQ.'SS ') THEN
            IF (I*J.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0*(CM2 (I - J) - CM2 (I+J) )
            ENDIF
         ELSEIF (TRI.EQ.'SC ') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 2D0*(SM2 (I - J)+SM2 (I+J) )
            ENDIF
         ELSEIF (TRI.EQ.'CC ') THEN
            RINT = 1D0 / 2D0*(CM2 (I+J)+CM2 (I - J) )
         ELSEIF (TRI.EQ.'SSS') THEN
            IF (I*J*K.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0*(SM2 (I - J+K)+SM2 (I+J - K)    &
               - SM2 (I - J - K) - SM2 (I+J+K) )
            ENDIF
         ELSEIF (TRI.EQ.'SSC') THEN
            IF (I*J.NE.0) THEN
               RINT = 1D0 / 4D0*(CM2 (I - J - K) - CM2 (I+J+K)    &
              +CM2 (I - J+K) - CM2 (I+J - K) )
            ELSE
               RINT = 0D0
            ENDIF
         ELSEIF (TRI.EQ.'SCC') THEN
            IF (I.EQ.0) THEN
               RINT = 0D0
            ELSE
               RINT = 1D0 / 4D0*(SM2 (I+J - K)+SM2 (I+J+K)    &
              +SM2 (I - J+K)+SM2 (I - J - K) )
            ENDIF
         ELSEIF (TRI.EQ.'CCC') THEN
            RINT = 1D0 / 4D0*(CM2 (I+J+K)+CM2 (I+J - K)       &
           +CM2 (I - J+K)+CM2 (I - J - K) )
         ENDIF
      ENDIF
      R = RINT
   END FUNCTION R
   !-----------------------------------------------------------------------
   !     END OF ri
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !> Returns 1 if N is even or -1 if N is odd
   pure integer function even_odd (N)
        implicit none
        integer, intent(in):: n
        even_odd = - 1
        IF (MOD (N, 2) .EQ.0) even_odd = 1
   END FUNCTION even_odd
   !-----------------------------------------------------------------------
   !> Swaps the three indices according to a 3 letter code
   pure SUBROUTINE tausch (TRI, I, J, K)
      implicit none
      CHARACTER(len=3), intent(inout):: TRI
      integer, intent(inout):: i, j, k
      integer:: N
      IF (TRI.EQ.'CS ') THEN
         N = I
         I = J
         J = N
         TRI = 'SC '
      ELSEIF (TRI.EQ.'SCS') THEN
         N = J
         J = K
         K = N
         TRI = 'SSC'
      ELSEIF (TRI.EQ.'CSS') THEN
         N = I
         I = K
         K = N
         TRI = 'SSC'
      ELSEIF (TRI.EQ.'CCS') THEN
         N = I
         I = K
         K = N
         TRI = 'SCC'
      ELSEIF (TRI.EQ.'CSC') THEN
         N = I
         I = J
         J = N
         TRI = 'SCC'
      ENDIF
      END SUBROUTINE TAUSCH
!-----------------------------------------------------------------------
!
!> Returns 2/(n*pi) in case n is odd.
!-----------------------------------------------------------------------
   pure double precision function two_over_nPi_odd (N)
      implicit none
      integer, intent(in):: N
      IF (even_odd(N) .EQ.1) THEN
         two_over_nPi_odd = 0D0
      ELSE
         two_over_nPi_odd = 2D0 / (N * DPI)
      ENDIF
   END FUNCTION two_over_nPi_odd
!-----------------------------------------------------------------------
!
!> Returns 1 id N is equal to 0. Returns 0 otherwise.
!-----------------------------------------------------------------------
   pure double precision function delta (N)
      implicit none
      integer, intent(in):: N
      delta = 0D0
      IF (N.EQ.0) delta = 1D0
   END FUNCTION delta
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function S1 (N)
      implicit none
      integer, intent(in):: N
      IF (even_odd (N) .EQ.1) THEN
         IF (N.EQ.0) THEN
            S1 = 0D0
         ELSE
            S1 = - 1D0 / (N * DPI)
         ENDIF
      ELSE
         S1 = (ri+ro) / (N * DPI)
      ENDIF
   END FUNCTION S1
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function C1 (N)
      implicit none
      integer, intent(in):: N
      IF (even_odd (N) .EQ.1) THEN
         IF (N.EQ.0) THEN
            C1 = RI+1D0 / 2
         ELSE
            C1 = 0D0
         ENDIF
      ELSE
         C1 = - 2D0 / N**2 / DPI**2
      ENDIF
   END FUNCTION C1
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function S2 (N)
      implicit none
      integer, intent(in):: N
      IF (even_odd (N) .EQ.1) THEN
         IF (N.EQ.0) THEN
            S2 = 0D0
         ELSE
            S2 = - (ri+ro) / DPI / N
         ENDIF
      ELSE
         S2 = (1+2*ri*ro ) / DPI / N - 4D0 / N**3 / DPI**3
      ENDIF
   END FUNCTION S2
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function C2 (N)
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in):: N
      IF (even_odd (N) .EQ.1) THEN
         IF (N.EQ.0) THEN
            C2 = ri*ro + 1D0 / 3
         ELSE
            C2 = 2D0 / N**2 / DPI**2
         ENDIF
      ELSE
         C2 = - 2D0*(ri+ro) / N**2 / DPI**2
      ENDIF
      END FUNCTION C2
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function S3 (N)
      implicit none
      integer, intent(in):: N
      IF (even_odd (N) .EQ.1) THEN
         IF (N.EQ.0) THEN
            S3 = 0D0
         ELSE
            S3 = - (1+3*ri+3*ri**2) / DPI / N+6 / DPI**3 / N**3
         ENDIF
      ELSE
         S3 = (1+3*ri+3*ri**2+2*ri**3) / DPI / N - (6+12*ri) / DPI**3 / N**3
      ENDIF
      END FUNCTION S3
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function C3 (N)
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in):: N
      IF (even_odd (N) .EQ.1) THEN
         IF (N.EQ.0) THEN
            C3 = 1D0 / 4+ri+3D0 / 2*ri**2+ri**3
         ELSE
            C3 = (3+6*ri) / DPI**2 / N**2
         ENDIF
      ELSE
         C3 = - (3+6*ri+6*ri**2) / DPI**2 / N**2+12 / DPI**4 / N**4
      ENDIF
      END FUNCTION C3
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function S4 (N)
      implicit none
      integer, intent(in):: N
      IF (even_odd (N) .EQ.1) THEN
         IF (N.EQ.0) THEN
            S4 = 0D0
         ELSE
            S4 = - (1+4*ri+6*ri**2+4*ri**3) / DPI / N+(12+24*ri) / DPI**3 / N**3
         ENDIF
      ELSE
         S4 = (1+4*ri+6*ri**2+4*ri**3+2*ri**4) / DPI /  &
         N - (12+24*ri+24*ri**2) / DPI**3 / N**3+48D0 / DPI**5 / N**5
      ENDIF
      END FUNCTION S4
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function C4(N)
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in):: N
      IF (even_odd(N).EQ.1) THEN
         IF (N.EQ.0) THEN
            C4 = 1D0 / 5 + ri + 2*ri**2 + 2*ri**3 + ri**4
         ELSE
            C4 = (4+12*RI+12*RI**2) / DPI**2 / N**2 - 24 / DPI**4 / N**4
         ENDIF
      ELSE
         C4 = - (4+12*ri+12*ri**2+8*ri**3) / DPI**2 / N**2 + (24+48*ri) / DPI**4 / N**4
      ENDIF
      END FUNCTION C4
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function SM1 (N)
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in):: N
      IF (N.EQ.0) THEN
         SM1 = 0D0
      ELSE
         SM1 = DCOS (N*ri*DPI)*DIS (ri, ro, DPI*N) - DSIN ( &
         N*ri*DPI)*DIC (ri, ro, N*DPI)
      ENDIF
      END FUNCTION SM1
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function CM1 (N)
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in):: N
      IF (N.EQ.0) THEN
         CM1 = DLOG ( ro / ri)
      ELSE
         CM1 = DCOS (N*ri*DPI)*DIC (ri, ro, N*DPI)+DSIN ( &
         N*ri*DPI)*DIS (ri, ro, N*DPI)
      ENDIF
      END FUNCTION CM1
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function SM2 (N)
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in):: N
      double precision:: NPi, NPiri, NPiro
      IF (N.EQ.0) THEN
         SM2 = 0D0
      ELSE
         NPi   = N*DPI
         NPiri = N*ri*DPI
         NPiro = N*ro*DPI
         SM2 = DCOS(NPiri)*( DSIN(NPiri)/ri - DSIN(NPiro)/ro + NPI*DIC(ri,ro,NPI) ) - &
               DSIN(NPiri)*( DCOS(NPiri)/ri - DCOS(NPiro)/ro - NPI*DIS(ri,ro,NPI) )
      ENDIF
      END FUNCTION SM2
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function CM2 (N)
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in):: N
      double precision:: NPi, NPiri, NPiro
      IF (N.EQ.0) THEN
         CM2 = 1D0/ri - 1D0/ro
      ELSE
         NPi   = N*DPI
         NPiri = N*ri*DPI
         NPiro = N*ro*DPI
         CM2 = DSIN(NPiri)*( DSIN(NPIRI)/ri - DSIN(NPiro)/ro + NPi*DIC (ri, ro, NPi) ) + &
               DCOS(NPiri)*( DCOS(NPIRI)/ri - DCOS(NPIro)/ro - NPi*DIS (ri, ro, NPi) )
      ENDIF
      END FUNCTION CM2
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function DIS (XMIN, XMAX, A)
!-----------------------------------------------------------------------
      implicit none
      double precision, intent(in):: XMIN, XMAX, A
      double precision:: X, dismin, dismax
!
      X = DABS (A*XMAX)
      IF (X.LT.1) THEN
         DISMAX = SIS (X)
      ELSE
         DISMAX = SIA (X)
      ENDIF
      IF (A*XMAX.LT.0) DISMAX = - DISMAX
      X = DABS (A*XMIN)
      IF (X.LT.1) THEN
         DISMIN = SIS (X)
      ELSE
         DISMIN = SIA (X)
      ENDIF
      IF (A*XMIN.LT.0) DISMIN = - DISMIN
      DIS = DISMAX - DISMIN
   END FUNCTION DIS
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function SIA (X)
!-----------------------------------------------------------------------
      implicit none
      double precision, intent(in):: x
      double precision, parameter, DIMENSION(4):: AF=(/ 38.027264D0, 265.187033D0, 335.677320D0, 38.102495D0 /)
      double precision, parameter, DIMENSION(4):: BF=(/ 40.021433D0, 322.624911D0, 570.236280D0, 157.105423D0 /)
      double precision, parameter, DIMENSION(4):: AG=(/ 42.242855D0, 302.757865D0, 352.018498D0, 21.821899D0 /)
      double precision, parameter, DIMENSION(4):: BG=(/ 48.196927D0, 482.485984D0, 1114.978885D0, 449.690326D0 /)
      double precision:: F,G
      F = (X**8 + AF(1)*X**6 + AF(2)*X**4 + AF(3)*X**2 + AF(4)&
      ) / (X**8 + BF(1)*X**6 + BF(2)*X**4 + BF(3)*X**2 + BF(4)&
      ) / X
      G = (X**8 + AG(1)*X**6 + AG(2)*X**4 + AG(3)*X**2 + AG(4)&
      ) / (X**8 + BG(1)*X**6 + BG(2)*X**4 + BG(3)*X**2 + BG(4)&
      ) / X**2
      SIA = DPI / 2D0 - F*DCOS (X) - G*DSIN (X)
   END FUNCTION SIA
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function SIS (X)
!-----------------------------------------------------------------------
      implicit none
      double precision, intent(in):: x
      double precision, parameter:: GENAU = 1.0D-7
      integer:: Jz, Jn, i, j
      double precision:: SISO,SISN, SISZ

      SISO = X
      DO I = 1, 200000
         JZ = 0
         JN = 0
         SISZ = 1D0
         SISN = 2D0*I+1D0
         do
            DO J = JZ+1, 2*I+1
               SISZ = SISZ*X
               IF (SISZ.GT.1D20) exit
            END DO
            JZ = J
            DO J = JN+1, 2*I+1
               SISN = SISN*J
               IF (SISN.GT.1D20) exit
            END DO
            JN = J
            SISZ = SISZ / SISN
            SISN = 1
            IF ( (JZ.ge.2*I+1) .and. (JN.ge.2*I+1) ) exit
         enddo
         SISN = 2D0*I+1D0
         JN = 0
!
         SIS = SISO + ((-1)**I)*SISZ
         IF (I.GT.1) THEN
            IF (DABS (1D0 - SIS / SISO) .LE.GENAU) exit
         ENDIF
         SISO = SIS
      END DO
   END FUNCTION SIS
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function DIC (XMIN, XMAX, A)
!-----------------------------------------------------------------------
      implicit none
      double precision, intent(in):: xmin, xmax, a
      double precision:: DICMAX, DICMIN, x
!
      X = DABS (A*XMAX)
      IF (X.LT.1.0d0) THEN
         DICMAX = CIS (X)
      ELSE
         DICMAX = CIA (X)
      ENDIF
      X = DABS (A*XMIN)
      IF (X.LT.1.0d0) THEN
         DICMIN = CIS (X)
      ELSE
         DICMIN = CIA (X)
      ENDIF
      DIC = DICMAX - DICMIN
   END FUNCTION DIC
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function CIA (X)
!-----------------------------------------------------------------------
      implicit none
      double precision, intent(in):: x
      double precision:: F, G
      double precision, parameter, dimension(4):: AF=(/ 38.027264D0, 265.187033D0, 335.677320D0, 38.102495D0 /)
      double precision, parameter, dimension(4):: BF=(/ 40.021433D0, 322.624911D0, 570.236280D0, 157.105423D0 /)
      double precision, parameter, dimension(4):: AG=(/ 42.242855D0, 302.757865D0, 352.018498D0, 21.821899D0 /)
      double precision, parameter, dimension(4):: BG=(/ 48.196927D0, 482.485984D0, 1114.978885D0, 449.690326D0 /)
      F = ( X**8 + AF(1)*X**6 + AF(2)*X**4 + AF(3)*X**2 + AF(4) ) / &
          ( X**8 + BF(1)*X**6 + BF(2)*X**4 + BF(3)*X**2 + BF(4) ) / X
      G = ( X**8 + AG(1)*X**6 + AG(2)*X**4 + AG(3)*X**2 + AG(4) ) / &
          ( X**8 + BG(1)*X**6 + BG(2)*X**4 + BG(3)*X**2 + BG(4) ) / X**2
      CIA = F*DSIN (X) - G*DCOS (X)
      END FUNCTION CIA
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
   pure double precision function CIS (X)
      implicit none
      double precision, PARAMETER:: C = 0.5772156649D0, GENAU = 1D-7
      integer:: Jz, Jn, i, j
      double precision, intent(in):: x
      double precision:: CISN, CISO, CISZ

      CISO = DLOG (X)+C
      DO I = 1, 200000
         JZ = 0
         JN = 0
         CISZ = 1D0
         CISN = 2D0*I
         do
            DO J = JZ+1, 2*I
               CISZ = CISZ*X
               IF (CISZ.GT.1D20) exit
            END DO
            JZ = J
            DO J = JN+1, 2*I
               CISN = CISN*J
               IF (CISN.GT.1D20) exit
            END DO
            JN = J
            CISZ = CISZ / CISN
            CISN = 1D0
            IF ( (JZ.ge.2*I) .and. (JN.ge.2*I) ) exit
         enddo
!
         CIS = CISO+(-1)**I*CISZ
         IF (I.GT.1) THEN
            IF (DABS (1D0 - CIS / CISO) .LE.GENAU) exit
         ENDIF
         CISO = CIS
      END DO
   END FUNCTION CIS
end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
