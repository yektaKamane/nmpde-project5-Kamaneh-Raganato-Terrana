#include <fstream>
#include <iostream>
#include <vector>

#include "Fisher_Kolmogorov_solver.hpp"

// Main function.
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  // Check if the correct number of command-line arguments is provided
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " <dimension> <mesh_file> [parameter_file]" << std::endl;
    return 1;
  }

  // Read the mandatory dimension from command line arguments and convert to unsigned int
  unsigned int dim;
  try
  {
    dim = std::stoul(argv[1]);
  }
  catch (const std::invalid_argument &)
  {
    std::cerr << "Invalid dimension: " << argv[1] << ". Must be a positive integer." << std::endl;
    return 1;
  }
  catch (const std::out_of_range &)
  {
    std::cerr << "Dimension out of range: " << argv[1] << std::endl;
    return 1;
  }

   // Read the mandatory mesh file name from command line arguments
  std::string mesh_file = argv[2];

  // Read the optional parameter file name or use a default value
  std::string parameter_file;
  if (argc >= 4)
  {
    parameter_file = argv[3];
  }
  else
  {
    parameter_file = "../input/test1.prm"; // Default parameter file
  }

  unsigned const int  convergence_test = 0;
  // Instantiate the problem with the specified dimension
  if (dim == 2)
  {
    FisherKol<2> problem(mesh_file, parameter_file, convergence_test);
    problem.setup();
    problem.solve();
  }
  else if (dim == 3)
  {
    FisherKol<3> problem(mesh_file, parameter_file, convergence_test);
    problem.setup();
    problem.solve();
  }
  else
  {
    std::cerr << "Unsupported dimension: " << dim << ". Only 2 and 3 are supported." << std::endl;
    return 1;
  }


  return 0;
}