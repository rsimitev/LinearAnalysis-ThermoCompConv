######## Makefile for lo ##################################
######## attention: the command must be preceeded by a "tab"
 
FCFLAGS:=-O3 -finit-local-zero -finit-real=zero -finit-integer=0 -fbounds-check -g -Wall -pedantic 
FCFLAGS:=-O3 -ffree-line-length-none -finit-local-zero -finit-real=zero -finit-integer=0 -cpp
#FCFLAGS:=-O2 -g -debug full -check all
LDFLAGS:= 

FC:=gfortran

default: $(ALL)

attention: 
	cat 0.README

include src/include.mk

%.o:%.f90
	$(FC) $(FCFLAGS) -c -o $@ $<


clean: 
	rm -vf *~ *.mod $(OBJS)

realclean: clean
	rm -vf lo

