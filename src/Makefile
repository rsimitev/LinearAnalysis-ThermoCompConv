######## Makefile for lo ##################################
######## attention: the command must be preceeded by a "tab"

FCFLAGS=-O3 -finit-local-zero -finit-real=zero -finit-integer=0 -g -Wall -pedantic 
#FCFLAGS:=-O2 -g -debug full -check all
LDFLAGS= 

FC:=gfortran
#FC:=ifort

default: lo

attention: 
	cat 0.README

OBJECTS := lolo.o r.o imsl.o pegasus.o

%.o:%.F
	$(FC) -c $(FCFLAGS) -o $@ $^

####### Main program ###############################################
lo: attention $(OBJECTS)
	$(FC) -DNMAX=760 $(FCFLAGS) -o $@ $(OBJECTS) $(LDFLAGS) losub.F

lara: lara_estream.o readla.o plms.o input.o
	$(FC) $(FCFLAGS) -o lara lara_estream.o readla.o plms.o input.o

clean: 
	rm -vf *~ *.o 
	
distclean: clean
	rm -vf lo lara lo_in plot.in
	rm -vf idl.x idl.y idl.z idl.batch plot.ps in out

