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

  const unsigned int dim = 2;

  FisherKol<dim> problem("../mesh/mesh-square-5.msh", "../input/test1.prm");

  problem.setup();
  problem.solve();


  return 0;
}