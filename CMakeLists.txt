cmake_minimum_required(VERSION 3.12.0)
project(project_title LANGUAGES CXX C)

include(../cmake-common.cmake)

add_executable(lab-01 src/main.cpp src/FisherKolmogorov.cpp)
deal_ii_setup_target(lab-01)
