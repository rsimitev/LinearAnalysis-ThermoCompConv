!*********************************************************************
!>    fn     : double precision function of only one variable
!!    Xmin   : min of interval
!!    Xmax   : max of interval
!!    rdx    : relative error in x
!!    df     : absolute error limit in f
!!    stepMax   : max. count of steps
!!    X0     : on output: calculated root
!!    info   : on output: 0: convergence; 
!!                        1: no convergence;
!!                        2: no zero in the interval;
!!                        3: overflow;
!!                        4: NaN;
   SUBROUTINE minimizer (fn, Xmin, xmax, rdx, df, stepMax, X0, info)
      implicit none
      double precision, intent(in):: xmin, xmax, rdx, df
      integer, intent(in):: stepMax !< max. count of steps
!>    0: convergence;\n 
!!    1: no convergence;\n
!!    2: no zero in the interval;
      integer, intent(out):: info
      double precision, intent(out):: x0 !< Calculated root
      double precision:: x1, x2, x3
      double precision:: f1, f2, f3, ff, dxdf
      integer:: n

      interface
         function fn(x)
            double precision:: fn
            double precision:: x
         end function
      end interface

      info = 0
      if(xmin>xmax) then
         X1 = Xmax
         X2 = Xmin
      else
         X1 = Xmin
         X2 = Xmax
      endif
      if (X1.ge.1.0d100.or.X2.ge.1.0d100) then
         info = 3
         X0 = huge(1.0d0)
         return
      endif
      F1 = FN (X1)
      if (F1.ne.F1) then
         info = 4
         x0 = 0.0d0
         return
      endif
      F2 = FN (X2)
      if (F2.ne.F2) then
         info = 4
         x0 = 0.0d0
         return
      endif
      ! Write(*,*) 'In minimizer: X1,F1, X2,F2',X1, F1, X2, F2
      N = 0
!
      FF = F1 * F2
!
      if (ff.gt.0.0d0) then
         info = 2
         x0 = huge(1.0d0)
         return
      endif
      ! Try our luck!
      if(abs(f1).lt.df) then
         x0 = x1
         info = 0
         return
      elseif(abs(f2).lt.df) then
         x0 = x2
         info = 0
         return
      endif

      do !mainloop
         if(n==0.and.x1*x2<0.0d0) then ! Be slick and try 0
            x3 = 0.0d0
         ! Use bisection rule at the start and every other step.
         elseif (mod(n-1,2)==0.or.n<6) then
            x3 = (x2 + x1)/2.0d0
         else
            ! Project along the line joining the end points
            dxdf = (x2-x1)/(f2-f1)
            x3 = x1 - f1 * dxdf
         endif
         if (X3.ge.1.0d100) then
            info = 3
            X0 = huge(1.0d0)
            return
         endif
         F3 = FN (X3)
         if (F3.ne.F3) then
            info = 4
            x0 = 0.0d0
            return
         endif
         ! Write(*,*) 'In minimizer: X3,F3',X3, F3
         N = N + 1
         FF = F2 * F3
         ! Shift end points appropriately
         IF (FF.LT.0.D0) THEN
            CALL swap (X1, F1, X3, F3)
            if (n.gt.3) f1=2.0d0*f1
         ELSE
            CALL swap (X2, F2, X3, F3)
            if (n.gt.3) f2=2.0d0*f2
         ENDIF
         ! If any of the limits is close enough to zero, exit
         if(abs(f1).lt.df) then
            x0 = x1
            info = 0
            exit
         elseif(abs(f2).lt.df) then
            x0 = x2
            info = 0
            exit
         endif
         ! If the relative distance between points is lower than
         ! this limit, assume the zero is in the middle and exit.
         if (abs((x2-x1)/(x2+x1)).lt.rdx) then
            x0 = (x2+x1)/2.0d0
            info = 0
            exit
         endif

         ! If we have gone over stepMax iterations,
         ! assume the zero is in the middle and exit.
         IF (N.GT.stepMax) THEN
            x0 = (x2+x1)/2.0d0
            info = 1
            exit
         ENDIF
      enddo !mainloop
   END SUBROUTINE minimizer

   !***************************************************************
   SUBROUTINE swap (X1, Y1, X2, Y2)
      implicit none
      double precision, intent(inout):: x1, y1, x2, y2
      double precision:: XH, YH
!
      XH = X1
      YH = Y1
      X1 = X2
      Y1 = Y2
      X2 = XH
      Y2 = YH
   END SUBROUTINE swap
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
