program testMinimizer
   implicit none
   !> The arguments to this program are the test parameters.
   !! Valid tests are:
   !! -\a m find the correct minimum in an interval
   !! -\a n correctly find there is no minimum in an interval
   character:: test
   integer:: info
   double precision:: X, Y

   info = 0
   call getarg(1,test)
   select case (test)
      case('m') 
         ! Test if we can find a minimium in an interval where the function
         ! crosses 0.
         call minimizer(my_f, 0.0d0, 20.0d0, 1.0d-5, 1.0d-5, 1000, X, info)
         Print*, X, my_f(X), info
         if (info==0) then
            if (abs(x-6.213).lt.1.0d-5) then
               Y = my_f(x)
               if (abs(y).lt.1.0d-5) then
                  print*, 'Pass.'
               else
                  print*, 'Fail.'
               endif
            else
               print*, 'Fail.'
            endif
         else
            print*, 'Fail.'
         endif
      case('n') 
         ! Test fo rwhen the interval does not contain a zeo
         call minimizer(my_f, 20.0d0, 30.0d0, 1.0d-5, 1.0d-5, 1000, X, info)
         Print*, X, my_f(X), info
         if(info.ne.2) then
            print*, 'Fail.'
         else
            print*, 'Pass.'
         endif
      case default
         print*, 'Usage : '
         print*, 'testMinimizer <test>'
         stop -1
   end select
contains 
   !*****************************************************
   !> A quadratic function that has a minimum at x=7.213,
   !! y=0
   double precision function my_f(x)
      implicit none
      double precision, intent(in):: x
      my_f = (x-7.213)**3+1
   end function
end program 
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
