c        program lolo
 
      implicit none
        character*60 infile,outfile
        integer lr,lr1
c	external iargc
c	external getarg

c---------------------------------------------------------
c       arg #1 - filename or usage ? 
        call getarg(1,infile)
        if (infile.eq.' ') then
           print*, 'Usage : '
	   print*, 'modify -h <in file> <out file> '
           stop
	endif
	if (infile.eq.'-h') then
           print*, 'Usage : '
	   print*, 'modify -h <in file> <out file>'
           stop
	endif

c       name of files with initial coefficients : argument #1
	lr=1
 	do while (infile(lr:).ne.' ') 
	   lr=lr+1
	enddo
c       name of files with modified coefficients : argument #2
	call getarg(2,outfile)
        lr1=1
 	do while (outfile(lr1:).ne.' ') 
	   lr1=lr1+1
	enddo

        print*,  infile(1:lr-1),' - ',outfile(1:lr1-1)

c        call losub('in','out')
        call losub(infile(1:lr-1),outfile(1:lr1-1))


        end











