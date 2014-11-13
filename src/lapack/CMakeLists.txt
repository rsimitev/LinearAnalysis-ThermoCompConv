set(MYLAPACK_SRCS 
   disnan.f 
   ieeeck.f 
   zgeqr2.f 
   zlacgv.f 
   zlarf.f 
   ztgevc.f 
   dlabad.f 
   ilaenv.f 
   zgeqrf.f 
   zlacpy.f 
   zlarft.f 
   zung2r.f 
   dladiv.f 
   ilazlc.f 
   zggbak.f 
   zladiv.f 
   zlartg.f 
   zungqr.f 
   dlaisnan.f 
   ilazlr.f 
   zggbal.f 
   zlange.f 
   zlascl.f 
   zunm2r.f 
   dlamch.f 
   iparmq.f 
   zggev.f 
   zlanhs.f 
   zlaset.f 
   zunmqr.f 
   dlapy2.f 
   lsame.f 
   zgghrd.f 
   zlarfb.f 
   zlassq.f 
   dlapy3.f 
   xerbla.f 
   zhgeqz.f 
   zlarfg.f 
   zrot.f )

add_library( mylapack ${MYLAPACK_SRCS} )
target_link_libraries(mylapack ${BLAS_LIBRARIES})

