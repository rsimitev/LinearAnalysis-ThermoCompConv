######## Makefile for lo ##################################
######## attention: the command must be preceeded by a "tab"
 
FCFLAGS:=-O3 -finit-local-zero -finit-real=zero -finit-integer=0 -fbounds-check -g -Wall -pedantic 
FCFLAGS:=-O3 -ffree-line-length-none -finit-local-zero -finit-real=zero -finit-integer=0 -cpp
#FCFLAGS:=-O2 -g -debug full -check all
LDFLAGS= 
LDFLAGS= 

FC=gfortran

default: lo

attention: 
	cat 0.README

OBJECTS:=parameters.o growthRate.o minimizer.o

lapack_OBJS:=$(addprefix lapack/lapack_routine/,disnan.o ieeeck.o zgeqr2.o zlacgv.o zlarf.o ztgevc.o dlabad.o ilaenv.o zgeqrf.o zlacpy.o zlarft.o zung2r.o dladiv.o ilazlc.o zggbak.o zladiv.o zlartg.o zungqr.o dlaisnan.o ilazlr.o zggbal.o zlange.o zlascl.o zunm2r.o dlamch.o iparmq.o zggev.o zlanhs.o zlaset.o zunmqr.o dlapy2.o lsame.o zgghrd.o zlarfb.o zlassq.o dlapy3.o xerbla.o zhgeqz.o zlarfg.o zrot.o)

lo.o: 

growthRate.o: parameters.o

%.o:%.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

####### Main program ###############################################
lo: main.f90 $(OBJECTS) $(lapack_OBJS)
	$(FC) $(FCFLAGS) -o $@ $(LDFLAGS) -Llib -lblas $^


clean: 
	rm *~ *.o *.mod $(lapack_OBJS)

realclean: clean
	lo

