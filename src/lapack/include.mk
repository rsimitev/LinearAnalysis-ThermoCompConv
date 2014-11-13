VPATH+=src/lapack/

lapack_OBJS:=disnan.o ieeeck.o zgeqr2.o zlacgv.o zlarf.o ztgevc.o dlabad.o ilaenv.o zgeqrf.o zlacpy.o zlarft.o zung2r.o dladiv.o ilazlc.o zggbak.o zladiv.o zlartg.o zungqr.o dlaisnan.o ilazlr.o zggbal.o zlange.o zlascl.o zunm2r.o dlamch.o iparmq.o zggev.o zlanhs.o zlaset.o zunmqr.o dlapy2.o lsame.o zgghrd.o zlarfb.o zlassq.o dlapy3.o xerbla.o zhgeqz.o zlarfg.o zrot.o

OBJS+=$(lapack_OBJS)
