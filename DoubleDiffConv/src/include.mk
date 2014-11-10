VPATH+=src/

include src/lapack/include.mk

lo_BIN:=lo

lo_OBJS:=parameters.o io.o growthRate.o minimizer.o

growthRate.o: parameters.o
io.o: parameters.o

####### Main program ###############################################
$(lo_BIN): main.f90 $(lo_OBJS) $(lapack_OBJS)
	$(FC) $(FCFLAGS) -o $@ $(LDFLAGS) -Llib -lblas $^

ALL+=$(lo_BIN)
OBJS+=$(lo_OBJS)
