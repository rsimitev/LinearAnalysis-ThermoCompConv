######## Makefile for lo ##################################
######## attention: the command must be preceeded by a "tab"

CCFLAGS=-c -O2 -o
LDFLAGS= 

FC=gfortran

attention: 
	cat 0.README

losub.o: losub.F
	$(FC) -DNMAX=760 $(CCFLAGS) losub.o losub.F

lo.o: main.f90
	$(FC) $(CCFLAGS) lo.o main.f90

r.o: r.F
	$(FC) $(CCFLAGS) r.o r.F

imsl.o: imsl.F
	$(FC) $(CCFLAGS) imsl.o imsl.F

pegasus.o: pegasus.F
	$(FC) $(CCFLAGS) pegasus.o pegasus.F

####### Main program ###############################################
lo:     attention losub.o lo.o r.o pegasus.o imsl.o 
	$(FC) $(LDFLAGS) -o lo lo.o losub.o r.o pegasus.o imsl.o 
	rm *.o


clean: 
	rm *~ *.o lo lo_in plot.in

