**********************************************************************
*     PLMS.FOR
*     R1.0 JW
*     R1.1 MA optimization of STOREPLM() and PLMS(), comments.
**********************************************************************
      FUNCTION PLMS(L,M,NTHETA)
**********************************************************************
*     provides the value of PLM at angle theta = theta(NTHETA)
*     plms() reads it from COMMON/PLMC/PLMST(NMT,NMLM), where it was
*     stored by storeplm().
**********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      PARAMETER(NMT=512,NMLM=4000)
C
      COMMON/PLMC/PLMST(NMT,NMLM),NMTC,NMLMC,NMTS,NMLMS,NML
C
      IF( L.GT.NML ) THEN
         WRITE(*,*) 'SORRY, PLM FOR THIS L HAS NOT BEEN STORED.',L
         WRITE(*,*) 'PLEASE INCREASE NML IN SUBROUTINE STOREPLM.'
         STOP
      ENDIF
C
      IF( L.LT.M ) THEN
         PLMS=0.D0
         RETURN
      ENDIF

C
C--    calculate positi.on LM:
C--    first jump to the right M:
C      LM=0
C      DO 100 I=0,M-1
C100   LM=LM+NML-I+1
C--   or as closed formula (MA):
      LM = M*(NML+1) - (M*(M-1))/2
C--   add offset for the right L:
      LM = LM+L-M+1
      IF( LM.GT.NMLMS ) THEN
         WRITE(*,*) 'SORRY, PLM HAS NOT BEEN STORED FOR L,M = ',L,M
         STOP
      ENDIF
C
      PLMS=PLMST(NTHETA,LM)
C
      RETURN
      END function
C
c--------------------------------------------------------------------
C   stores the values PLM in PLMST(nmt,nmlm) for theta=theta(1..nmtheta).
C   storeplm() has to be called once before PLMS() is called.
c-- THETA HAS TO BE GIVEN IN DEGREES.
c   The plm's here are fully normalised!
C   There is no sign change.
      subroutine storeplm(theta,nmtheta)
      implicit real*8(a-h,o-z)
      dimension theta(*)
      PARAMETER(NMT=512,NMLM=4000)
c
      common/plmc/plmst(nmt,nmlm),nmtc,nmlc,nmts,nmls,nml
c
      pi   = 3.141592653589793d0
      nmtc = nmt
      nmlc = nmlm
c     max. value of M and L:
      nml=80
c
      if(nmtheta.GT.nmt) then
        write(*,*) 'STOREPLM: nmtheta.GT.nmt!'
        stop
      endif

c     i is the index of the theta point
      do i=1, nmtheta
c        theta in radians
         thetar = pi*theta(i)/180.d0
         sinth  = dsin(thetar)
         costh  = dcos(thetar)
c        sin(th)**m, m=0:
         sinthm = 1.d0
         lm = 0
c
c        indices run on l => m faster and then on m=0...nml
         do m=0, nml
            fac=1.d0
            do j=3,2*m+1,2
               fac=fac*dble(j)/dble(j-1)
            enddo
            plm=dsqrt(fac)
            if( sinth.ne.0.d0 ) then
               plm = plm*sinthm
c              sin(th)**m for next m:
               sinthm = sinthm*sinth
            elseif( m.ne.0 ) then
               plm=0.d0
            endif
c
c-- plm for l=m: P_l^l =
            lm = lm + 1
            if( lm.gt.nmlm ) then
               write(*,*) 'Too small dimension nmlm in storelpm.',nmlm
               stop
            endif
            plmst(i,lm) = plm
            plm1 = 0.d0
c
c-- plm for l>m: P_l^m =
            do l=m+1, nml
               plm2 = plm1
               plm1 = plm
               plm  = costh * dsqrt( dble( (2*l-1)*(2*l+1) ) / dble( (l-m)*(l+m) )  ) * plm1 -
     -                     dsqrt( dble( (2*l+1)*(l+m-1)*(l-m-1) ) / dble( (2*l-3)*(l-m)*(l+m) ) ) * plm2
               lm=lm+1
               if( lm.gt.nmlm ) then
                  write(*,*) 'Too small dimension nmlm in storelpm.', nmlm
                  stop
               endif
               plmst(i,lm)=plm
            enddo
         enddo
      enddo
c
      nmts = nmtheta
      nmls = lm
c
      return
      end subroutine
c
c------------------------------------------------------------------------
