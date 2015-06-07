FC=gfortran

FFLAGS=-O2 -g3 --fixed-line-length-none

default:lara

lara_estream.o: lara_estream.f90
	$(FC) $(FFLAGS) -c lara_estream.f90

lara: lara_estream.o readla.o plms.o input.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o lara lara_estream.o readla.o plms.o input.o

clean:
	@rm -fv *.o *~ idl.x idl.y idl.z idl.batch plot.ps in out

realclean: clean
	@rm -fv lara
