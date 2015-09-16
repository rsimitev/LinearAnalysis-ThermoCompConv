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
   implicit none
   double precision, parameter:: PI = 3.14159265358979D0
   integer, parameter:: nr=33, nt=180, np=360
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
      print*, 'simplePlot <in file> <domain>'
      print*, '  <domain> can be one of 2D or 3D.'
      stop
   endif
   
   call getarg(2,domain)
   print*, 'Will plot ', trim(infile),' in ',domain,'.'

   write(*,*) 'Reading data from ',INFILE,'...'
   
   CALL READLA(INFILE)
   Write(*,*) 'Plotting...'
   select case (domain)
      case('2D')
         CALL Plot_2D()
      case('3D')
         CALL Plot_3D()
   end select
   Write(*,*) 'Done!'

contains

   !------------------------------------------------------------------------
   !     calculates the field Z and makes one subplot.
   SUBROUTINE Plot_2D()
      implicit none
      double precision:: THETA, r, phi
      double precision:: ri, z
      integer:: i, j, k
      !-- CALCULATION OF INNER AND OUTER RADIUS:
      RI = ETA/(1.D0-ETA)

      write(*,*) 'computing the fields...'
      ! Plot the equator
      theta  = 90.0d0
      
      ! Find the highest l or m that we need to compute
      ! Reuse k for that so we don't have to create a new variable
      k=0
      DO i=1, size(eigenVector)
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
   SUBROUTINE READLA(STARTFILE)
      IMPLICIT none
      CHARACTER(len=*)::STARTFILE
      character(len=128):: title
      integer, parameter:: unit_in=12
      integer:: LMIN, LD, nElements, i

      OPEN(unit_in,FILE=STARTFILE,STATUS='old')
      title='#'
      do while (title(1:1)=='#')
         READ(unit_in,*) title
         title = adjustl(title)
      enddo
      backspace(unit_in)
      READ(unit_in,*) M0, Truncation, symmetry
      READ(unit_in,*) tau, Rt, Rc, Pt, Le, eta

      call setLminAndLD(Symmetry, m0, LMIN, LD) 
      call setEigenProblemSize(LMIN,LD,truncation,M0)
      nElements = getEigenProblemSize()
      allocate (eigenVector(nElements))
      do i=1, nElements
         read(unit_in,*) eigenVector(i)%fieldCode, &
                         eigenVector(i)%l, &
                         eigenVector(i)%m, &
                         eigenVector(i)%n, &
                         eigenVector(i)%RealPart, &
                         eigenVector(i)%ImagPart
      enddo
      close(unit_in)
   END SUBROUTINE READLA

end PROGRAM simplePlot
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
