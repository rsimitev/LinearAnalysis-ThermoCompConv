include_directories(${CMAKE_CURRENT_BINARY_DIR}/..)
set (simplePlot_SRCS
     plms.f90
     simplePlot.f90
   )
add_executable(simplePlot ${simplePlot_SRCS})
target_link_libraries(simplePlot linearProblem)
install(TARGETS simplePlot RUNTIME DESTINATION bin)