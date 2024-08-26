#include <fstream>
#include <iostream>
#include <vector>

#include "Fisher_Kolmogorov_solver.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  const unsigned int               mpi_rank =
    Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

    // Check if the correct number of command-line arguments is provided
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " [parameter_file]" << std::endl;
    return 1;
  }



  // Read the optional parameter file name or use a default value
  std::string parameter_file;
  if (argc >= 3)
  {
    parameter_file = argv[1];
  }
  else
  {
    parameter_file = "../input/test1D_convergence.prm"; // Default parameter file
  }

    unsigned const int  convergence_test = 0;
    FisherKol<1> problem("", parameter_file, convergence_test);

    problem.setup();
    problem.solve();

  return 0;
}
