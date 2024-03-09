cmake_minimum_required(VERSION 3.12.0)
project(project_title LANGUAGES CXX C)

include(common/cmake-common.cmake)

add_executable(test2D src/test2D.cpp src/FisherKolmogorov2D.cpp)
deal_ii_setup_target(test2D)
