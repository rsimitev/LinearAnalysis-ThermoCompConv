set(GLO_SRCS
   growthRate.f90  
   io.f90  
   main.f90  
   minimizer.f90  
   parameters.f90
   )

add_executable(glo ${GLO_SRCS})
target_link_libraries(glo ${LAPACK_LIBRARIES})
set_target_properties(glo PROPERTIES LINK_FLAGS "${LAPACK_LINKER_FLAGS}")
install(TARGETS glo RUNTIME DESTINATION bin)
