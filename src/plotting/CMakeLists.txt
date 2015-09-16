include_directories(${CMAKE_CURRENT_BINARY_DIR}/..)
set (glo_plot_SRCS
     plms.f90
     glo_plot.f90
   )
add_executable(glo_plot ${glo_plot_SRCS})
target_link_libraries(glo_plot linearProblem)
install(TARGETS glo_plot RUNTIME DESTINATION bin)

set (glo_compute_and_plot_SRCS
     plms.f90
     glo_compute_and_plot.f90
   )
add_executable(glo_compute_and_plot ${glo_compute_and_plot_SRCS})
target_link_libraries(glo_compute_and_plot linearProblem)
install(TARGETS glo_compute_and_plot RUNTIME DESTINATION bin)
