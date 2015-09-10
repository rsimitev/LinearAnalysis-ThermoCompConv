module plms_mod

   private
   double precision, parameter:: pi   = 3.141592653589793d0
   double precision, allocatable:: plmst(:,:)
   integer:: Nlm_i, maxl_i
   public:: plms, storeplm

contains

   !**********************************************************************
   !*     provides the value of PLM at angle theta = theta(NTHETA)
   !*     plms() reads it from COMMON/PLMC/PLMST(NMT,NMLM), where it was
   !*     stored by storeplm().
   !**********************************************************************
   double precision FUNCTION PLMS(L,M,iTHETA)
      IMPLICIT none
      integer, intent(in):: l, m, itheta
      integer:: lm

      IF( L .GT. maxl_i ) THEN
         WRITE(*,*) 'SORRY, PLM FOR THIS L HAS NOT BEEN STORED.',L
         STOP
      ENDIF

      IF( L.LT.M ) THEN
         PLMS=0.D0
         RETURN
      ENDIF

      LM = M*(maxl_i+1) - (M*(M-1))/2
      !-- add offset for the right L:
      LM = LM+L-M+1
      IF( LM.GT.Nlm_i ) THEN
         WRITE(*,*) 'SORRY, PLM HAS NOT BEEN STORED FOR L,M = ',L,M
         STOP
      ENDIF

      PLMS=PLMST(iTHETA,LM)
   END function

   !--------------------------------------------------------------------
   !   stores the values PLM in PLMST(nmt,nmlm) for theta=theta(1..nmtheta).
   !   storeplm() has to be called once before PLMS() is called.
   !-- THETA HAS TO BE GIVEN IN DEGREES.
   !   The plm's here are fully normalised!
   !   There is no sign change.
   subroutine storeplm(theta, ntheta, maxml)
      implicit none
      integer, intent(in):: ntheta, maxml
      double precision, intent(in):: theta(ntheta)
      double precision:: thetar, sinth, sinthm, costh, fac
      double precision:: plm, plm1, plm2
      integer:: lm, i, m, l, j

      maxl_i = maxml
      Nlm_i  = (maxl_i+2)*(maxl_i+1)/2
      if(allocated(plmst)) deallocate(plmst)
      allocate(plmst(ntheta,Nlm_i))

      !i is the index of the theta point
      do i=1, ntheta
         ! theta in radians
         thetar = pi*theta(i)/180.d0
         sinth  = dsin(thetar)
         costh  = dcos(thetar)
         ! sin(th)**m, m=0:
         sinthm = 1.d0

         lm = 0
         ! indices run on l => m faster and then on m=0...nml
         do m=0, maxl_i
            fac=1.d0
            do j=3,2*m+1,2
               fac=fac*dble(j)/dble(j-1)
            enddo
            plm=dsqrt(fac)
            if( sinth.ne.0.d0 ) then
               plm = plm*sinthm
               ! sin(th)**m for next m:
               sinthm = sinthm*sinth
            elseif( m.ne.0 ) then
               plm=0.d0
            endif

            !-- plm for l=m: P_l^l =
            lm = lm + 1
            plmst(i,lm) = plm
            plm1 = 0.d0

            !-- plm for l>m: P_l^m =
            do l=m+1, maxl_i
               plm2 = plm1
               plm1 = plm
               plm  = costh * dsqrt( dble( (2*l-1)*(2*l+1) ) / dble( (l-m)*(l+m) )  ) * plm1 - dsqrt( dble( (2*l+1)*(l+m-1)*(l-m-1) ) / dble( (2*l-3)*(l-m)*(l+m) ) ) * plm2
               lm=lm+1
               plmst(i,lm)=plm
            enddo
         enddo
      enddo

      return
   end subroutine

end module