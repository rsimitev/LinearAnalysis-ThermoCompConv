!************************************************************************
!-- PROGRAM TO read data from Galerkinprogram of J.W. and convert it for IDL.
!--
!-- Input:   stdin  (short version of LARA.EXE for DISSPLA)
!--
!-- Output:
!--
!------------------------------------------------------------------------
PROGRAM simplePlot
   use parameters
   use growthRateMod
   use plms_mod
   use io
   implicit none
   double precision, parameter:: PI = 3.14159265358979D0
   integer:: nr=31, nt=180, np=361
   CHARACTER(len=60):: INFILE
   character(len=2):: domain

   type:: eigenElement
      integer:: n
      integer:: l
      integer:: m
      character:: fieldCode
      double precision:: realPart
      double precision:: imagPart
   end type
   type(eigenElement), target, allocatable:: eigenVector(:)
!---------------------------------------------------------
!  arg #1 - filename or usage ?
   call getarg(1,infile)
   if (trim(infile).eq.'' .or. trim(infile).eq.'-h') then
      print*, 'Usage : '
      print*, 'glo_plot <in file> <domain>'
      print*, '  <domain> can be one of 2D or 3D.'
      stop
   endif
   
   call getarg(2,domain)
   call init(trim(infile))
   call computeModes()

   Write(*,*) 'Plotting...'
   select case (domain)
      case('2D')
         CALL Plot_2D()
      case('3D')
         CALL Plot_3D()
      case default
         stop
   end select
   Write(*,*) 'Done!'

contains

   !**********************************************************************
   !> Initialises things.
   SUBROUTINE init(inputfile)
      implicit none
      CHARACTER(len=*), intent(in):: inputfile

      ! ----Default values:
      call setDefaults()
      ! ----INPUT:
      call readConfigFileNew(inputfile)

      ! ---- doesn't work for M=0 !!!!!!
      IF(M0.LT.1) THEN
        write(*,*) 'The code does not work for M0<1. ', M0, ' --> 1'
        M0 = 1
      ENDIF
      if(4*m0.gt.Np) then
         Np = 4*m0+1
         Nt = 2*m0
         Write(*,*) "Warning: This file is going to be huge!"
      endif
      Nr = max(Nr,floor(Nt*acos(eta)/pi)+1)

      call GrowthRateInit(Rt, Rc, Pt, Le, tau, eta, m0, Symmetry, Truncation)
      call setVariableParam(VariablePar)
   END subroutine

   !------------------------------------------------------------------------
   !     calculates the field Z and makes one subplot.
   SUBROUTINE Plot_2D()
      implicit none
      double precision:: THETA, r, phi
      double precision:: ri
      integer:: i, k
      !-- CALCULATION OF INNER AND OUTER RADIUS:
      RI = ETA/(1.D0-ETA)

      write(*,*) 'computing the fields...'
      ! Plot the equator
      theta  = 90.0d0
      
      ! Find the highest l or m that we need to compute
      ! Reuse k for that so we don't have to create a new variable
      k=0
      DO i=1, size(eigenVector,1)
         if (eigenVector(i)%m.gt.k) k = eigenVector(i)%m
         if (eigenVector(i)%l.gt.k) k = eigenVector(i)%l
      enddo

      !-- BESTIMMUNG DER PLM(THETA) , ABSPEICHERUNG:
      
      CALL STOREPLM( (/THETA/), 1, k )

      open(20,file ='glo-render_2D.general',STATUS =  'UNKNOWN')
      Write(20,*) 'grid = ',np+1,'x',nr+1
      Write(20,*) 'format = ascii'
      Write(20,*) 'interleaving = field'
      Write(20,*) 'majority = column'
      Write(20,*) 'field = locations, field0'
      Write(20,*) 'structure = 2-vector, scalar'
      Write(20,*) 'type = float, float'
      Write(20,*) ''
      Write(20,*) 'end'
      do i=0, nr
         r = ri + dble(i)/dble(nr)
         do k=0, np
            phi = k*(360.0/np)
            Write(20,*) r*cos(phi*pi/180.0), &
                        r*sin(phi*pi/180.0), &
                        flow_r(eigenVector, r, PHI, 1)
         enddo
      enddo
      close(20)
   end subroutine Plot_2D

   !------------------------------------------------------------------------
   !     calculates the field Z and makes one subplot.
   SUBROUTINE Plot_3D()
      implicit none
      double precision:: THETA(Nt), r, phi, dtheta
      double precision:: ri
      integer:: i, j, k
      !-- CALCULATION OF INNER AND OUTER RADIUS:
      RI = ETA/(1.D0-ETA)

      write(*,*) 'computing the fields...'
      ! Avoid the poles
      dtheta  = 180.d0/(nt+1)
      DO I = 1, Nt
         THETA(I) = (I-0.5d0)*dtheta
      enddo
      
      ! Find the highest l or m that we need to compute
      ! Reuse k for that so we don't have to create a new variable
      k=0
      DO i=1, size(eigenVector)
         if (eigenVector(i)%m.gt.k) k = eigenVector(i)%m
         if (eigenVector(i)%l.gt.k) k = eigenVector(i)%l
      enddo

      !-- BESTIMMUNG DER PLM(THETA) , ABSPEICHERUNG:
      
      CALL STOREPLM( THETA, Nt, k )

      open(20,file ='glo-render_3D.general',STATUS =  'UNKNOWN')
      Write(20,*) 'grid = ',nt,'x',np+1,'x',nr+1
      Write(20,*) 'format = ascii'
      Write(20,*) 'interleaving = field'
      Write(20,*) 'majority = column'
      Write(20,*) 'field = locations, field0'
      Write(20,*) 'structure = 3-vector, scalar'
      Write(20,*) 'type = float, float'
      Write(20,*) ''
      Write(20,*) 'end'
      do i=0, nr
         r = ri + dble(i)/dble(nr)
         do k=0, np
            phi = k*(360.0/np)
            do j=1, nt
               Write(20,*) r*sin((theta(j))*pi/180.0)*cos(phi*pi/180.0), &
                           r*sin((theta(j))*pi/180.0)*sin(phi*pi/180.0), &
                           r*cos((theta(j))*pi/180.0), &
                           flow_r(eigenVector, r, PHI, j)
            enddo
         enddo
      enddo
      close(20)
   end subroutine Plot_3D

   !------------------------------------------------------------------------
   !> Radial flow: U_r = L_2/r v
   double precision function flow_r(eigenVector,R,PHI,iTHETA)
      IMPLICIT none
      type(eigenElement)::eigenVector(:)
      double precision, intent(in):: R, phi
      double precision:: ri, epsm, pphi
      double precision:: RFT, RFTR, RFTI
      integer, intent(in):: iTheta
      integer:: i
      integer:: l, m, n

      PPHI = PHI*PI/180.D0
      !-- CALCULATION OF INNER AND OUTER RADIUS:
      RI = ETA/(1.D0-ETA)

      flow_r = 0.0d0
      DO i=1, size(eigenVector)
         if (eigenVector(i)%fieldCode.ne.'V') cycle
         ! Prefactor for Legendre Associated Polynomials
         IF( eigenVector(i)%M.EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         endif

         l = eigenVector(i)%l
         m = eigenVector(i)%m
         n = eigenVector(i)%n
         RFT = EPSM*l*(l+1)*PLMS(L,M,iTheta)*DSIN( N*PI*(R-RI) ) / R

         RFTR =  RFT * eigenVector(I)%realPart * DCOS( M*PPHI )
         RFTI = -RFT * eigenVector(I)%ImagPart * DSIN( M*PPHI )

         flow_r = flow_r + RFTR+RFTI
      enddo
   end function

   !------------------------------------------------------------------------
   SUBROUTINE computeModes()
      IMPLICIT none
      integer:: nElements
      complex(8), allocatable:: ZEVEC(:), zew(:), ZEVAL(:,:)
      integer:: i, ni, li, lti, lpi, ld, ii
      integer:: LMAX, NIMAX, LMIN, eqdim


      call setLminAndLD(Symmetry, m0, LMIN, LD) 
      call setEigenProblemSize(LMIN,LD,truncation,M0)
      nElements = getEigenProblemSize()
      allocate( ZEVEC(NElements), zew(NElements), ZEVAL(NElements,NElements))
      allocate (eigenVector(nElements))
      ! Recoompute critical state modes and eigenvectors
      call computeGrowthRateModes(sort=.TRUE., zew=zew, zeval=zeval)
      ! Most unstable mode will be the first
      zevec(:) = zeval(:,1)

      eqdim=nElements/4

      LMAX=2*Truncation+M0-1
      ! poloidal flow:
      I=1
      do LI=LMIN,LMAX,LD
         LPI = LI
         NIMAX=INT( DBLE(2*Truncation+1-LI+M0)/2 )
         do NI=1,NIMAX
            eigenVector(i)%fieldCode = 'V'
            eigenVector(i)%l         = LPI
            eigenVector(i)%m         = M0
            eigenVector(i)%n         = NI
            eigenVector(i)%RealPart  = DREAL(ZEVEC(I+0*eqdim))
            eigenVector(i)%ImagPart  = DIMAG(ZEVEC(I+0*eqdim))
            I=I+1
         enddo
      enddo

      i=1
      ii=eqdim+1
      DO LI=LMIN,LMAX,LD
         IF( Symmetry.EQ.2 ) THEN
         LTI=LI+1
         ELSEIF( Symmetry.EQ.1 ) THEN
         LTI=LI-1
         ELSEIF( Symmetry.EQ.0 ) THEN
         LTI=LI
         ENDIF
         NIMAX=INT( DBLE(2*Truncation+1-LI+M0)/2 )
         DO NI=1,NIMAX
            eigenVector(ii)%fieldCode = 'W'
            eigenVector(ii)%l         = LTI
            eigenVector(ii)%m         = M0
            eigenVector(ii)%n         = NI
            eigenVector(ii)%RealPart  = DREAL(ZEVEC(I+1*eqdim))
            eigenVector(ii)%ImagPart  = DIMAG(ZEVEC(I+1*eqdim))
            I=I+1
            Ii=Ii+1
         enddo
      enddo

      ! temperature
      i=1
      ii=2*eqdim+1
      DO LI=LMIN,LMAX,LD
         NIMAX=INT( DBLE(2*Truncation+1-LI+M0)/2 )
         DO NI=1,NIMAX
            eigenVector(ii)%fieldCode = 'T'
            eigenVector(ii)%l         = LI
            eigenVector(ii)%m         = M0
            eigenVector(ii)%n         = NI
            eigenVector(ii)%RealPart  = DREAL(ZEVEC(I+2*eqdim))
            eigenVector(ii)%ImagPart  = DIMAG(ZEVEC(I+2*eqdim))
            I=I+1
            Ii=Ii+1
         enddo
      enddo

      ! composition
      i=1
      ii=3*eqdim+1
      DO LI=LMIN,LMAX,LD
         NIMAX=INT( DBLE(2*Truncation+1-LI+M0)/2 )
         DO NI=1,NIMAX
            eigenVector(ii)%fieldCode = 'G'
            eigenVector(ii)%l         = LI
            eigenVector(ii)%m         = M0
            eigenVector(ii)%n         = NI
            eigenVector(ii)%RealPart  = DREAL(ZEVEC(I+3*eqdim))
            eigenVector(ii)%ImagPart  = DIMAG(ZEVEC(I+3*eqdim))
            I=I+1
            Ii=iI+1
         enddo
      enddo
      
   END SUBROUTINE computeModes

end PROGRAM simplePlot
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
