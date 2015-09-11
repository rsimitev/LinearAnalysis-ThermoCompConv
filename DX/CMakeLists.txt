include_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

configure_file(glo-ur-eq.net.in glo-ur-eq.net @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/glo-ur-eq.net
              DESTINATION ${DXPROGS})
