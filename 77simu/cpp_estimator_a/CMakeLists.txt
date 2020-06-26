cmake_minimum_required(VERSION 3.0.0)
project(cpp_estimator_a VERSION 0.1.0)

include(CTest)
enable_testing()

add_definitions("-O2")
add_executable(
    cpp_estimator_a
    main.cpp
    around_csv.cpp
    value_io.cpp
    calculator.cpp
    progress.cpp
)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
