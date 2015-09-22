add_subdirectory(lapack)

set (linearProblem_SRCS
     growthRate.f90
     parser.f90
     minimizer.f90
     parameters.f90
     glo_constants.f90
   )
    
add_library(linearProblem ${linearProblem_SRCS}) 
target_link_libraries(linearProblem mylapack ${BLAS_LIBRARIES})
set_target_properties(linearProblem PROPERTIES LINK_FLAGS "${BLAS_LINKER_FLAGS}")
   
set(GLO_SRCS io.f90 glo.f90 )
add_executable(glo ${GLO_SRCS})
target_link_libraries(glo linearProblem)
install(TARGETS glo RUNTIME DESTINATION bin)

set (CriticalRaEff_SRCS CritRaEff.f90 CritRaEff_io.f90)
add_executable(CritRaEff ${CriticalRaEff_SRCS})
target_link_libraries(CritRaEff linearProblem)
install(TARGETS CritRaEff RUNTIME DESTINATION bin)

if(BUILD_UTILS)
   add_subdirectory(modeTracker)
   add_subdirectory(plotting)
endif(BUILD_UTILS)
if(BUILD_TESTS)
   add_subdirectory(test)
endif(BUILD_TESTS)

