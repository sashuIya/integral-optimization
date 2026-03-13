cmake_minimum_required(VERSION 3.14)
project(IntegralOptimization VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -O3")
endif()

# Directories
include_directories(include)

# Executable (currently monolithic)
add_executable(integral_optimization main.cpp src/system_state.cc src/ode_solver.cc src/shooting_method.cc)

# Google Test
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(baseline_test tests/baseline_test.cc)
target_link_libraries(baseline_test gtest_main)

add_executable(system_state_test tests/system_state_test.cc src/system_state.cc)
target_link_libraries(system_state_test gtest_main)

add_executable(ode_solver_test tests/ode_solver_test.cc src/ode_solver.cc src/system_state.cc)
target_link_libraries(ode_solver_test gtest_main)

add_executable(shooting_method_test tests/shooting_method_test.cc src/shooting_method.cc src/ode_solver.cc src/system_state.cc)
target_link_libraries(shooting_method_test gtest_main)


include(GoogleTest)
gtest_discover_tests(baseline_test)
gtest_discover_tests(system_state_test)
gtest_discover_tests(ode_solver_test)
gtest_discover_tests(shooting_method_test)
