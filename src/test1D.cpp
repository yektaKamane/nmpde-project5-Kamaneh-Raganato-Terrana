#include <deal.II/base/convergence_table.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "FisherKolmogorov1D.hpp"

// Main function.
int
main(int argc, char * argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  
  const unsigned int N = 200;
  const unsigned int r = 1;
  const double T       = 20.0;
  const double deltat  = 0.1;

  FisherKol problem(N, r, T, deltat);

  problem.setup();
  problem.solve();

  return 0;
}



