set(MYLAPACK_SRCS 
   disnan.f 
   dlabad.f 
   dladiv.f 
   dlaisnan.f 
   dlamch.f 
   dlapy2.f 
   dlapy3.f 
   ieeeck.f 
   ilaenv.f 
   ilazlc.f 
   ilazlr.f 
   iparmq.f 
   lsame.f 
   xerbla.f 
   zgeqr2.f 
   zgeqrf.f 
   zggbak.f 
   zggbal.f 
   zggev.f 
   zgghrd.f 
   zhgeqz.f 
   zlacgv.f 
   zlacpy.f 
   zladiv.f 
   zlange.f 
   zlanhs.f 
   zlarf.f 
   zlarfb.f 
   zlarfg.f 
   zlarft.f 
   zlartg.f 
   zlascl.f 
   zlaset.f 
   zlassq.f 
   zrot.f
   ztgevc.f 
   zung2r.f 
   zungqr.f 
   zunm2r.f 
   zunmqr.f)

add_library( mylapack ${MYLAPACK_SRCS} )
target_link_libraries(mylapack ${BLAS_LIBRARIES})

