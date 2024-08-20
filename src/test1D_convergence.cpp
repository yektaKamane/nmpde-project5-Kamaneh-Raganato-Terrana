#include <fstream>
#include <iostream>
#include <vector>

#include "Fisher_Kolmogorov_solver.hpp"
#include "convergence.hpp"

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

    FisherKol<1> problem("", parameter_file);

    problem.setup();
    problem.solve();
    // non so, perché vorremmo salvare tutto, oppure trovare il valore finale di c e vedere quando converge e quando no?
    // poi qualcosa tipo switch case per i tre valori di alpha e d
  // }

  return 0;
}
