program testMinimizer
   use parameters
   use growthRateMod
   implicit none
   !> The arguments to this program are the test parameters.
   !! Valid tests are:
   !! -\a m find the correct minimum in an interval
   !! -\a n correctly find there is no minimum in an interval
   character(len=256):: file1, file2
   integer:: error

   call getarg(1,file1)
   call getarg(2,file2)
   Write(*,*) trim(file1), trim(file2)
   call testMAT(trim(file1), trim(file2), error)
   if(error.eq.0) then
      Write(*,*) 'Pass.'
   else
      Write(*,*) 'Fail.'
   endif
end program 
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
