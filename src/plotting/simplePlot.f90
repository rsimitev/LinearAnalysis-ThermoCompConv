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
   CHARACTER*40:: INFILE,OUTFILE

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
   if (infile.eq.' ') then
      print*, 'Usage : '
      print*, 'lara <in file> <out file>'
      stop
   endif
   if (infile.eq.'-h') then
      print*, 'Usage : '
      print*, 'lara <in file> <out file>'
      stop
   endif

   call getarg(2,outfile)
   print*,  trim(infile),' - ',trim(outfile)

   write(*,*) 'reading data from ',INFILE,'...'
   CALL READLA(INFILE)
   write(*,*) '...done'
   
   CALL Plot()

contains

   !------------------------------------------------------------------------
   !     calculates the field Z and makes one subplot.
   SUBROUTINE Plot()
      implicit none
      double precision:: THETA(Nt), r, phi, dtheta
      double precision:: ri, z
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

      open(20,file ='glo-render.dat',STATUS =  'UNKNOWN')
      Write(20,*) '#,nt,x,np,x,nr'
      Write(20,*) '#',nt,'x',np+1,'x',nr+1
      Write(20,*) '#========================'
      do i=0, nr
         r = ri + dble(i)/dble(nr)
         do k=0, np
            phi = k*(360.0/np)
            do j=1, nt
               Z = flow_r(eigenVector, r, PHI, j)
               Write(20,*) r*sin((theta(j))*pi/180.0)*cos(phi*pi/180.0), &
                           r*sin((theta(j))*pi/180.0)*sin(phi*pi/180.0), &
                           r*cos((theta(j))*pi/180.0), &
                           z
            enddo
         enddo
      enddo
      close(20)
   end subroutine plot

   !------------------------------------------------------------------------
   !> Radial flow: U_r = L_2/r v
   double precision function flow_r(eigenVector,R,PHI,iTHETA)
      IMPLICIT none
      type(eigenElement)::eigenVector(:)
      double precision, intent(in):: R, phi
      double precision:: ri, rft, epsm, pphi
      integer, intent(in):: iTheta
      integer:: i
      integer:: l, m, n

      PPHI = PHI*PI/180.D0
      !-- CALCULATION OF INNER AND OUTER RADIUS:
      RI = ETA/(1.D0-ETA)

      flow_r = 0.0d0
      DO i=1, size(eigenVector)
         if (eigenVector(i)%fieldCode.ne.'V') cycle
         IF( eigenVector(i)%M.EQ.0 ) THEN
            EPSM = 1.D0
         ELSE
            EPSM = 2.D0
         endif

         l = eigenVector(i)%l
         m = eigenVector(i)%m
         n = eigenVector(i)%n
         RFT = EPSM*l*(l+1)*PLMS(L,M,iTheta)*DSIN( N*PI*(R-RI) ) / R

         RFT = RFT * eigenVector(I)%realPart * DCOS( M*PPHI )
         RFT = RFT - RFT * eigenVector(I)%ImagPart * DSIN( M*PPHI )

         flow_r = flow_r + RFT
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
