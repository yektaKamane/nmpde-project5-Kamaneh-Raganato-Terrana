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

# add_executable(test2D_B src/test2D_B.cpp src/Fisher_Kolmogorov_solver_convergence.cpp)
# deal_ii_setup_target(test2D_B)

add_executable(test1D_convergence src/test1D_convergence.cpp src/FisherKolmogorov1D_convergence.cpp)
deal_ii_setup_target(test1D_convergence)

add_executable(test2D_convergence src/test2D_convergence.cpp src/Fisher_Kolmogorov_solver_convergence.cpp)
deal_ii_setup_target(test2D_convergence)

add_executable(test3D_convergence src/test3D_convergence.cpp src/Fisher_Kolmogorov_solver_convergence.cpp)
deal_ii_setup_target(test3D_convergence)

# add_executable(test1D src/test1D.cpp src/FisherKolmogorov1D.cpp)
# deal_ii_setup_target(test1D)

# Working on axonal direction in 2D and 3D
# add_executable(main src/main.cpp src/Fisher_Kolmogorov_solver_convergence.cpp)
# deal_ii_setup_target(main)