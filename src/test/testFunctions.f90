module testFunctions
      
contains 
   !*****************************************************
   !> A quadratic function that has a minimum at x=7.213,
   !! y=0
   double precision function my_f(x)
      implicit none
      double precision, intent(in):: x
      my_f = (x-7.213)**3+1
   end function
end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
