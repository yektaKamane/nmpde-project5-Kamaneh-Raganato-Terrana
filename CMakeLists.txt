cmake_minimum_required(VERSION 3.12.0)
project(project_title LANGUAGES CXX C)

include(common/cmake-common.cmake)

# add_executable(test2D src/test2D.cpp src/FisherKolmogorov2D.cpp)
# deal_ii_setup_target(test2D)

# add_executable(test3D src/test3D.cpp src/FisherKolmogorov3D.cpp)
# deal_ii_setup_target(test3D)

# add_executable(test3D_B src/test3D_B.cpp src/Fisher_Kolmogorov_solver.cpp)
# deal_ii_setup_target(test3D_B)

# add_executable(test2D_B src/test2D_B.cpp src/Fisher_Kolmogorov_solver.cpp)
# deal_ii_setup_target(test2D_B)

add_executable(test2D_B src/test2D_B.cpp src/Fisher_Kolmogorov_solver_convergence.cpp)
deal_ii_setup_target(test2D_B)