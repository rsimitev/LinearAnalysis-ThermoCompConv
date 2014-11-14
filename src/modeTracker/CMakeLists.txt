add_subdirectories(modeTracker)

add_executable(modeTracker main.f90)
install(TARGETS modeTracker RUNTIME DESTINATION bin)
