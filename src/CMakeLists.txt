add_subdirectory(lapack)

set(GLO_SRCS
   growthRate.f90
   parser.f90
   io.f90
   glo.f90
   minimizer.f90
   parameters.f90
   )

add_executable(glo ${GLO_SRCS})
target_link_libraries(glo mylapack ${BLAS_LIBRARIES})
set_target_properties(glo PROPERTIES LINK_FLAGS "${BLAS_LINKER_FLAGS}")
install(TARGETS glo RUNTIME DESTINATION bin)

set(CriticalRaEff_SRCS
    CritRaEff.f90
    CritRaEff_io.f90
    growthRate_RaEff.f90
    parameters.f90
    parser.f90
    minimizer.f90)

add_executable(CritRaEff ${CriticalRaEff_SRCS})
target_link_libraries(CritRaEff mylapack ${BLAS_LIBRARIES})
set_target_properties(CritRaEff PROPERTIES LINK_FLAGS "${BLAS_LINKER_FLAGS}")

if(BUILD_UTILS)
   add_subdirectory(modeTracker)
endif(BUILD_UTILS)
if(BUILD_TESTS)
   add_subdirectory(test)
endif(BUILD_TESTS)

