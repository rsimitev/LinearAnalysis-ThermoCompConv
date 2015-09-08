! Reads stdin containing four columns:
! - the value of the first parameter being varied;
! - the value of the second parameter being varied;
! - the growth rate for this set of parameters;
! - the m that corresponds to this set of parameters.
! Values of the second column are expected to vary faster.
! All other parameters are considered fixed.
! 
! Outputs to stdout the values of the first, second and last column
! for the point imediately after the third column changes sign.
program giveMeM
   implicit none
   real:: tau, Rtc, gr, tau_old, Rtc_old, gr_old
   integer:: m, error
   error = 0
   read(*,*) tau_old, Rtc_old, gr_old, m
   do while (error==0)
      read(*,*, iostat=error) tau, Rtc, gr, m
      if(tau==tau_old) then
         if (gr*gr_old.le.0.0e0) then
            Write(*,*) tau, Rtc, m
         endif
      else
         tau_old = tau
      endif
      gr_old  = gr
   enddo

end program
