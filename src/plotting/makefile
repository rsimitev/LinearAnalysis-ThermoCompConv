FC=gfortran

FFLAGS=-O2  -g3
LDFLAGS=-O2  -g3

lara_estream.o: lara_estream.f90
	$(FC) $(FFLAGS) -c $<

lara: plms.o input.o lara_estream.o
	$(FC) $(LDFLAGS) -o $@ $^

clean:
	rm *.o *~ idl.x idl.y idl.z idl.batch plot.ps in out
