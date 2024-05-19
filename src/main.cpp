#include <deal.II/base/convergence_table.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "Fisher_Kolmogorov_solver.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  // const unsigned int               mpi_rank =
  //   Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  const unsigned int degree = 1;

  const double T      = 2.0;
  const double deltat = 0.1;
  const unsigned int dim = 2;

  FisherKol<dim> problem("../mesh/mesh-square-40.msh", degree, T, deltat, "../input/test1.prm");

  problem.setup();
  problem.solve();

  return 0;
}