      program lolo
      implicit none
      character*60 infile,outfile

c---------------------------------------------------------
c     arg #1 - filename or usage ? 
      call getarg(1,infile)
      if (infile.eq.' ') then
         print*, 'Usage : '
         print*, 'lo <in file> <out file>'
         stop
      endif
      if (infile.eq.'-h') then
         print*, 'Usage : '
         print*, 'lo <in file> <out file>'
         stop
      endif

      call getarg(2,outfile)
      print*,  trim(infile),' - ',trim(outfile)
c     call losub('in','out')
      call losub(trim(infile),trim(outfile))
      end program
