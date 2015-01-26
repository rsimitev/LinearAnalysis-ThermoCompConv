include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(testMinimizer ../minimizer.f90 testFunctions.f90 testMinimizer.f90)

add_test(NAME FindMinimum COMMAND testMinimizer m)
set_tests_properties(FindMinimum PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")

add_test(NAME CorrectlyIgnoreInterval COMMAND testMinimizer n)
set_tests_properties(CorrectlyIgnoreInterval PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")


add_executable(testMatrixGeneration ../parameters.f90 ../growthRate.f90 testMatrixGeneration.f90)
target_link_libraries(testMatrixGeneration mylapack ${BLAS_LIBRARIES})
add_test(NAME CorrectMatrixGeneration COMMAND 
   testMatrixGeneration ${CMAKE_CURRENT_SOURCE_DIR}/data/ZA_ref.txt
   ${CMAKE_CURRENT_SOURCE_DIR}/data/ZB_ref.txt)
set_tests_properties(CorrectMatrixGeneration PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")
