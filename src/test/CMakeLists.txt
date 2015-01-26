include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(testMinimizer ../minimizer.f90 testFunctions.f90 testMinimizer.f90)

add_test(NAME FindMinimum COMMAND testMinimizer m)
set_tests_properties(FindMinimum PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")

add_test(NAME CorrectlyIgnoreInterval COMMAND testMinimizer n)
set_tests_properties(CorrectlyIgnoreInterval PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")
