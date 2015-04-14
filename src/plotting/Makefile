FC=gfortran

FFLAGS=-O2 -g3 --fixed-line-length-none

lara_estream.o: lara_estream.f
	$(FC) $(FFLAGS) -c lara_estream.f

lara: lara_estream.o readla.o plms.o input.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o lara lara_estream.o readla.o plms.o input.o

clean:
	rm *.o *~ idl.x idl.y idl.z idl.batch plot.ps in out
