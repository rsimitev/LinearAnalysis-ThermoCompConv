include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_executable(testMinimizer ../minimizer.f90 testMinimizer.f90)

add_test(FindMinimum testMinimizer m)
set_tests_properties(FindMinimum PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")

add_test(CorrectlyIgnoreInterval testMinimizer n)
set_tests_properties(CorrectlyIgnoreInterval PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")
