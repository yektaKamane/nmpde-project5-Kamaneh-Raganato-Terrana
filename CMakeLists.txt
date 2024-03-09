cmake_minimum_required(VERSION 3.12.0)
project(project_title LANGUAGES CXX C)

include(./common/cmake-common.cmake)

add_executable(solver src/main.cpp src/Fisher_Kolmogorov_solver.cpp)
deal_ii_setup_target(solver)
